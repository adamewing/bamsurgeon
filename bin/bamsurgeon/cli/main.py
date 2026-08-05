#!/usr/bin/env python

'''
One entry point for both engines.

Takes a single VCF (or a legacy varfile) describing SNVs, indels and SVs
together, applies them, and writes one truth VCF.

Structural variants are applied first, then small variants on top of the SV
output. That ordering is not an implementation convenience -- it is the only
one that works. The SV engine excludes every read over its interval and
replaces them with reads simulated from a reference template, which carries
no small variants, so anything spiked in first is reverted. Running the other
way round, the read-editing engine happily edits the simulated reads, and a
small variant inside a duplication or inversion survives.

The alternative considered was merging both donors and calling replace_reads
once. That is one pass instead of two, but it destroys every small variant
falling inside an SV's excluded read set, which is precisely the interaction
this ordering avoids.
'''

import os
import sys
import argparse
import logging

from argparse import Namespace
from uuid import uuid4

import pysam

import bamsurgeon.makevcf as makevcf

from bamsurgeon.aligners import SUPPORTED_ALIGNERS
from bamsurgeon.snvindel import cluster_small, resolve_haplosize, run_spikein
from bamsurgeon.svengine import run_sv_spikein
from bamsurgeon.varinput import read_variants, sample_read_length, warn_sv_overlap

FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def sv_args(args, out_bam):
    ''' the Namespace bamsurgeon.svengine expects '''
    return Namespace(
        bamFileName=args.bamFileName, refFasta=args.refFasta,
        outBamFile=out_bam, varFileName=args.varFileName,
        svfrac=args.vaf, maxlibsize=args.maxlibsize, pad=args.pad,
        readlen=args.readlen, mindepth=args.sv_mindepth,
        maxdepth=args.sv_maxdepth, maxdfrac=args.maxdfrac, minctglen=None,
        maxmuts=None, cnvfile=args.cnvfile, donorbam=args.donorbam,
        ismean=args.ismean, issd=args.issd, simerr=args.simerr,
        procs=args.procs, inslib=args.inslib, aligner=args.aligner,
        alignopts=args.alignopts, alignerthreads=args.alignerthreads,
        tagreads=args.tagreads, skipmerge=False,
        keepsecondary=args.keepsecondary, debug=args.debug,
        tmpdir=args.tmpdir + '.sv', seed=args.seed, allowN=args.allowN,
        maxopen=args.maxopen, vcf='',
    )


def small_args(args, in_bam, out_bam):
    ''' the Namespace bamsurgeon.snvindel expects '''
    return Namespace(
        bamFileName=in_bam, refFasta=args.refFasta, outBamFile=out_bam,
        varFileName=args.varFileName, snvfrac=args.snvfrac,
        mutfrac=args.vaf, numsnvs=0, cnvfile=args.cnvfile,
        coverdiff=args.coverdiff, indel_coverdiff=args.indel_coverdiff,
        haplosize=args.haplosize, procs=args.procs,
        mindepth=args.mindepth,
        maxdepth=args.maxdepth, minmutreads=args.minmutreads,
        avoidreads=args.avoidreads, nomut=args.nomut,
        ignoreref=args.ignoreref, ignoresnps=args.ignoresnps,
        insane=args.insane, force=args.force, single=args.single,
        maxopen=args.maxopen, requirepaired=args.requirepaired,
        tagreads=args.tagreads, skipmerge=False,
        ignorepileup=args.ignorepileup, aligner=args.aligner,
        alignerthreads=args.alignerthreads, alignopts=args.alignopts,
        tmpdir=args.tmpdir + '.small', seed=args.seed, vcf='',
    )


def sort_and_index(bam):
    '''
    The SV stage appends its simulated reads after the target reads, so its
    output is not coordinate ordered. The small-variant stage reads it with
    fetch and pileup, which need a sorted, indexed BAM.
    '''
    tmp = bam + '.sorting.bam'
    pysam.sort('-o', tmp, bam)
    os.replace(tmp, bam)
    pysam.index(bam)


def main(args):
    logger.info("starting %s called with args: %s" % (sys.argv[0], ' '.join(sys.argv)))

    requests = read_variants(args.varFileName, 'any', default_vaf=args.vaf)

    if not requests:
        sys.exit('no variants read from %s' % args.varFileName)

    svs = [r for r in requests if r.is_sv]
    small = [r for r in requests if not r.is_sv]

    logger.info('read %d SV and %d small variant request(s)' % (len(svs), len(small)))

    # reported for information: applying SVs first means these survive, but
    # a small variant inside a deletion still has nothing left to sit on
    if svs and small:
        warn_sv_overlap(requests)

    records = []
    current = args.bamFileName
    intermediate = None

    if svs:
        intermediate = '%s.sv.%s.bam' % (args.outBamFile.rsplit('.bam', 1)[0], uuid4().hex[:8])
        logger.info('applying %d structural variant(s) -> %s' % (len(svs), intermediate))
        records += run_sv_spikein(sv_args(args, intermediate), requests=svs)
        sort_and_index(intermediate)
        current = intermediate

    if small:
        haplosize, read_length = resolve_haplosize(
            Namespace(bamFileName=current, refFasta=args.refFasta,
                      haplosize=args.haplosize))
        clusters = cluster_small(small, haplosize, read_length)
        logger.info('applying %d small variant(s) in %d cluster(s)'
                    % (len(small), len(clusters)))
        records += run_spikein(small_args(args, current, args.outBamFile),
                               clusters, 'bamsurgeon')
    elif svs:
        os.rename(intermediate, args.outBamFile)
        intermediate = None

    if intermediate and not args.debug:
        for f in (intermediate, intermediate + '.bai'):
            if os.path.exists(f):
                os.remove(f)

    var_basename = '.'.join(os.path.basename(args.varFileName).split('.')[:-1])
    bam_basename = '.'.join(os.path.basename(args.outBamFile).split('.')[:-1])

    vcf_fn = args.vcf + bam_basename + '.bamsurgeon.' + var_basename + '.vcf'

    makevcf.write_vcf(records, args.refFasta, vcf_fn)

    logger.info('%d variant(s) spiked in; truth VCF written to %s'
                % (len(records), vcf_fn))


def run():
    parser = argparse.ArgumentParser(
        prog='bamsurgeon',
        description='spike SNVs, indels and SVs into a BAM from one VCF')

    parser.add_argument('-v', '--varfile', dest='varFileName', required=True,
                        help='VCF (or legacy varfile) describing the variants to add')
    parser.add_argument('-f', '--bamfile', dest='bamFileName', required=True,
                        help='sam/bam/cram file from which to obtain reads')
    parser.add_argument('-r', '--reference', dest='refFasta', required=True,
                        help='reference genome, indexed with bwa index and samtools faidx')
    parser.add_argument('-o', '--outbam', dest='outBamFile', required=True,
                        help='.bam file name for output')

    # -s is deliberately absent: it means --snvfrac in addsnv/addindel and
    # --svfrac in addsv, so a single letter cannot carry both here
    parser.add_argument('--vaf', default=0.5, type=float,
                        help='default variant allele fraction, overridden by INFO/VAF or FORMAT/AF (default 0.5)')
    parser.add_argument('--snvfrac', default=1, type=float,
                        help='maximum allowable linked SNP MAF, for avoiding haplotypes (default 1)')

    parser.add_argument('-p', '--procs', default=1, type=int,
                        help='split into multiple processes (default 1)')
    parser.add_argument('--seed', default=None, type=int,
                        help='seed random number generation')
    parser.add_argument('--picardjar', default=None, help=argparse.SUPPRESS)  # deprecated, no effect
    parser.add_argument('-c', '--cnvfile', default=None,
                        help='tabix-indexed list of absolute copy number values, used to adjust VAFs')

    parser.add_argument('--aligner', default='backtrack',
                        help='supported aligners: ' + ','.join(sorted(SUPPORTED_ALIGNERS)))
    parser.add_argument('--alignopts', default=None,
                        help='aligner-specific options as option1:value1,option2:value2,...')
    parser.add_argument('--alignerthreads', default=1, type=int,
                        help='threads used per realignment (default 1)')

    # depth cutoffs measure different things at different times, so the two
    # engines keep their own
    small = parser.add_argument_group('small variants (SNV, indel)')
    small.add_argument('--mindepth', default=10, type=int,
                       help='minimum read depth at the site (default 10)')
    small.add_argument('--maxdepth', default=2000, type=int,
                       help='maximum read depth at the site (default 2000)')
    small.add_argument('-d', '--coverdiff', default=0.9, type=float,
                       help='allowed input/output coverage ratio after realignment, SNVs (default 0.9)')
    small.add_argument('--indel-coverdiff', default=0.1, type=float,
                       help='same for indels, which lose more coverage on realignment (default 0.1)')
    small.add_argument('-z', '--haplosize', default=0,
                       help="group sites within this distance, or 'auto' for the sampled read length (default 0)")
    small.add_argument('--minmutreads', default=3, type=int,
                       help='minimum number of mutated reads to output per site')
    small.add_argument('--avoidreads', default=None,
                       help='file of read names to avoid')
    small.add_argument('--ignoreref', action='store_true', default=False,
                       help='make mutations even if the ALT matches the reference')
    small.add_argument('--ignoresnps', action='store_true', default=False,
                       help='skip the scan for confounding SNPs sharing reads')
    small.add_argument('--insane', action='store_true', default=False,
                       help='skip the coverage sanity check')
    small.add_argument('--ignorepileup', action='store_true', default=False,
                       help='do not check pileup depth in mutation regions')
    small.add_argument('--requirepaired', action='store_true', default=False,
                       help='skip mutations if unpaired reads are present')
    small.add_argument('--nomut', action='store_true', default=False, help='dry run')

    sv = parser.add_argument_group('structural variants')
    sv.add_argument('--sv-mindepth', default=10, type=int,
                    help='minimum read depth at a breakend (default 10)')
    sv.add_argument('--sv-maxdepth', default=2000, type=int,
                    help='maximum read depth at a breakend (default 2000)')
    sv.add_argument('-l', '--maxlibsize', default=600, type=int,
                    help='maximum fragment length of the sequencing library')
    sv.add_argument('--pad', default=None, type=int,
                    help='reference flank each side of a breakpoint (default 2x --maxlibsize)')
    sv.add_argument('--readlen', default=None, type=int,
                    help='simulated read length (default: modal read length from the BAM)')
    sv.add_argument('--maxdfrac', default=0.1, type=float,
                    help='maximum discordant read fraction (default 0.1)')
    sv.add_argument('--ismean', default=300, help='mean insert size (default 300)')
    sv.add_argument('--issd', default=70, help='insert size standard deviation (default 70)')
    sv.add_argument('--simerr', default=0.0, help='error rate for simulated reads')
    sv.add_argument('--inslib', default=None,
                    help='FASTA library of insertion sequences')
    sv.add_argument('--donorbam', default=None,
                    help='bam supplying duplication interior reads; without it the interior is simulated')
    sv.add_argument('--allowN', action='store_true', default=False,
                    help='allow N in the reference slice, replace with A and warn')
    sv.add_argument('--keepsecondary', action='store_true', default=False,
                    help='keep secondary reads in final BAM')

    parser.add_argument('--force', action='store_true', default=False,
                        help='force mutations regardless of nearby SNPs or low coverage')
    parser.add_argument('--single', action='store_true', default=False,
                        help='input BAM is single-ended')
    parser.add_argument('--tagreads', action='store_true', default=False,
                        help='add BS tag to altered reads')
    parser.add_argument('--maxopen', default=1000, type=int,
                        help='maximum number of open files during merge (default 1000)')
    parser.add_argument('--tmpdir', default='bamsurgeon.tmp',
                        help='temporary directory (default bamsurgeon.tmp)')
    parser.add_argument('--debug', action='store_true', default=False,
                        help='retain all intermediates')
    parser.add_argument('--vcf', default='',
                        help='directory prefix for the output truth VCF')

    args = parser.parse_args()

    if args.pad is None:
        args.pad = 2 * int(args.maxlibsize)

    if args.picardjar is not None:
        logger.warning('--picardjar is deprecated and has no effect: BAM to FASTQ '
                       'conversion no longer uses picard')

    main(args)


if __name__ == '__main__':
    run()

#!/usr/bin/env python

'''
The structural variant spike-in engine.

Builds the mutated haplotype from reference slices (see svtemplate), simulates
reads over it, and replaces the corresponding reads in the target BAM.

bin/addsv.py is argparse plus a call into here, matching addsnv and addindel.
'''

import re
import os
import sys
import random
import traceback
import subprocess
import argparse
import pysam

import bamsurgeon.replace_reads as rr
import bamsurgeon.svtemplate as svt
import bamsurgeon.makevcf as makevcf

from bamsurgeon.aligners import remap_fastq, SUPPORTED_ALIGNERS
from bamsurgeon.records import MutationRecord, MutationResult
from bamsurgeon.varinput import read_variants
from bamsurgeon.common import *
from uuid import uuid4
from shutil import move
from concurrent.futures import ProcessPoolExecutor

import logging
FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


INTERVAL_KINDS = ('DEL', 'DUP', 'INV', 'INS')


def _read_in_region(read, chrom, start, end):
    ''' primary, properly-mated, and fully contained in the interval '''
    read_end = read.reference_start + read.query_length
    pair_end = read.next_reference_start + read.query_length
    return not (read.is_duplicate or read.is_secondary or read.is_supplementary or
                read.is_unmapped or read.mate_is_unmapped or
                read.next_reference_name != chrom or
                pair_end > end or read.next_reference_start < start or
                read_end > end or read.reference_start < start)


def get_reads(bam_file, fasta_ref, chrom, start, end, svfrac, salt):
    '''
    Yield the svfrac subset of reads in the interval.

    Two passes: selection is by exact count over the whole candidate set, so
    the population has to be known before anything can be emitted. Selecting
    on qname rather than per-read keeps both mates on the same side of the
    decision, as the previous per-read hash also did.
    '''
    bam = pysam.AlignmentFile(bam_file, reference_filename=fasta_ref)

    qnames = set()
    for read in bam.fetch(chrom, start, end):
        if _read_in_region(read, chrom, start, end):
            qnames.add(read.query_name)

    selected = select_qnames(qnames, svfrac, salt)

    for read in bam.fetch(chrom, start, end):
        if read.query_name in selected and _read_in_region(read, chrom, start, end):
            yield read

    bam.close()


def runwgsim(newseq, readlen, pemean, pesd, tmpdir, nsimreads, mutid='null',
             err_rate=0.0, seed=None):
    ''' wrapper function for wgsim, could swap out to support other reads simulators (future work?) '''

    basefn = tmpdir + '/' + mutid + ".wgsimtmp." + str(uuid4())
    fasta = basefn + ".fasta"
    fq1 = basefn + ".1.fq"
    fq2 = basefn + ".2.fq"

    with open(fasta, 'w') as fout:
        fout.write(">" + mutid + "\n" + newseq + "\n")

    logger.info("%s template len: %d" % (mutid, len(newseq)))
    logger.info("%s num. sim. reads: %d" % (mutid, nsimreads))
    logger.info("%s read length: %d" % (mutid, readlen))
    logger.info("%s PE mean outer distance: %f" % (mutid, pemean))
    logger.info("%s PE outer distance SD: %f" % (mutid, pesd))
    logger.info("%s error rate: %f" % (mutid, err_rate))

    wgsim_args = ['wgsim', '-e', str(err_rate), '-d', str(pemean), '-s', str(pesd),
                  '-N', str(nsimreads), '-1', str(readlen), '-2', str(readlen),
                  '-r', '0', '-R', '0', fasta, fq1, fq2]

    seed = 1 if seed == 0 else seed # Fix for wgsim thinking 0 is no seed
    if seed is not None: wgsim_args += ['-S', str(seed)]

    logger.info(str(wgsim_args))
    subprocess.check_call(wgsim_args)

    os.remove(fasta)

    return (fq1, fq2)


def singleseqfa(file, mutid='null'):
    with open(file, 'r') as fasta:
        header = None
        seq = ''
        for line in fasta:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    logger.warning("%s multiple entries found in %s only using the first" % (mutid, file))
                header = line.lstrip('>')
            else:
                seq += line
    return seq


def load_inslib(infa):
    seqdict = {}

    with open(infa, 'r') as fa:
        seqid = ''
        seq   = ''
        for line in fa:
            if line.startswith('>'):
                if seq != '':
                    seqdict[seqid] = seq
                seqid = line.lstrip('>').strip()
                seq   = ''
            else:
                assert seqid != ''
                seq = seq + line.strip()

    if seqid not in seqdict and seq != '':
        seqdict[seqid] = seq

    return seqdict


def discordant_fraction(bamfile, fasta_ref, chrom, start, end):
    r = 0
    d = 0
    bam = pysam.AlignmentFile(bamfile, reference_filename=fasta_ref)
    for read in bam.fetch(chrom, start, end):
        r += 1
        if not read.is_proper_pair:
            d += 1

    if r > 0:
        return float(d)/float(r)
    else:
        return 0.0


def add_donor_reads(args, mutid, tmpbamfn, chrom, left_bnd, right_bnd, svfrac):
    ''' interior coverage for a duplication, taken from --donorbam '''
    tmpbam = pysam.AlignmentFile(tmpbamfn, reference_filename=args.refFasta)

    outbamfn = '%s/%s.%s.dup.merged.bam' % (args.tmpdir, mutid, str(uuid4()))
    outbam = pysam.AlignmentFile(outbamfn, 'wb', template=tmpbam)
    for read in tmpbam.fetch(until_eof=True):
        outbam.write(read)

    with pysam.AlignmentFile(args.donorbam, reference_filename=args.refFasta) as donorbam:
        cover_donor = donorbam.count(contig=chrom, start=left_bnd, end=right_bnd) / float(right_bnd-left_bnd)
    with pysam.AlignmentFile(args.bamFileName, reference_filename=args.refFasta) as origbam:
        cover_orig = origbam.count(contig=chrom, start=left_bnd, end=right_bnd) / float(right_bnd-left_bnd)

    donor_norm_factor = cover_orig * svfrac / cover_donor
    if donor_norm_factor > 1.0:
        logger.warning('%s: donor_norm_factor %f > 1.0. This means donor bam has less coverage than required.' % (mutid, donor_norm_factor))

    logger.info('%s: DUP donor coverage normalisation factor: %f' % (mutid, donor_norm_factor))
    logger.info('%s: fetch donor reads from %s-%d-%d' % (mutid, chrom, left_bnd, right_bnd))

    nreads = 0
    for read in get_reads(args.donorbam, args.refFasta, chrom, left_bnd, right_bnd, donor_norm_factor, args.salt):
        read.query_name = read.query_name + '_donor_' + mutid
        outbam.write(read)
        nreads += 1

    outbam.close()

    logger.info('%s: using %d donor reads from %s' % (mutid, nreads, args.donorbam))

    return outbamfn


def resolve_insertion(args, req, mutid):
    ''' turn the INS spec into literal sequence '''
    spec = req.insseqfile

    if spec is None:
        return req.insseq, None

    if spec == 'RND':
        assert args.inslib is not None, 'INS RND requires --inslib'
        name = random.choice(sorted(args.inslib.keys()))
        logger.info("%s chose sequence from insertion library: %s" % (mutid, name))
        return args.inslib[name], name

    if spec.startswith('INSLIB:'):
        assert args.inslib is not None, 'INSLIB: requires --inslib'
        name = spec.split(':', 1)[1]
        assert name in args.inslib, '%s not found in insertion library' % name
        logger.info("%s specify sequence from insertion library: %s" % (mutid, name))
        return args.inslib[name], name

    return singleseqfa(spec, mutid=mutid), os.path.basename(spec)


def check_depth(bamfile, mutid, chrom, positions, mindepth, maxdepth):
    for pos in positions:
        depth = bamfile.count(chrom, pos-1, pos)
        if depth < mindepth:
            logger.warning('%s skipping due to insufficient depth at %s:%d (%d)' % (mutid, chrom, pos, depth))
            return False
        if depth > maxdepth:
            logger.warning('%s skipping due to excessive depth at %s:%d (%d)' % (mutid, chrom, pos, depth))
            return False
    return True


def makemut(args, req, alignopts):
    if args.seed is not None:
        random.seed(args.seed + int(req.start))

    kind = req.kind
    chrom = req.chrom
    start = int(req.start)
    end = int(req.end)
    vaf = float(req.vaf)

    # the index keeps mutids unique when a varfile repeats an interval.
    # wgsim names simulated reads after the mutid, so a collision would make
    # two mutations' reads indistinguishable to replace_reads.
    mutid = '_'.join(map(str, (chrom, start, end, kind, req.index)))

    bamfile = pysam.AlignmentFile(args.bamFileName, reference_filename=args.refFasta)
    reffile = pysam.Fastafile(args.refFasta)

    logdir = 'addsv_logs_' + os.path.basename(args.outBamFile)
    logfile = open(os.path.join(logdir, os.path.basename(args.outBamFile) + '_' + mutid + '.log'), 'w')

    result = MutationResult()

    try:
        pad = int(args.pad)

        # copy number adjusts the requested VAF, as before
        if args.cnvfile:
            cnv = pysam.Tabixfile(args.cnvfile, 'r')
            if chrom in cnv.contigs:
                for cnregion in cnv.fetch(chrom, start, end):
                    cn = float(cnregion.strip().split()[3])
                    logger.info("%s copy number in sv region: %f" % (mutid, cn))
                    vaf = vaf / cn
                    assert vaf <= 1.0, 'copy number from %s must be at least 1: %s' % (args.cnvfile, cnregion.strip())
                    logger.info("%s adjusted VAF: %f" % (mutid, vaf))

        breakends = [start, end] if kind != 'BND' else [start]
        if kind == 'BND':
            if not check_depth(bamfile, mutid, req.mate_chrom, [req.mate_pos],
                               args.mindepth, args.maxdepth):
                return result

        if not check_depth(bamfile, mutid, chrom, breakends, args.mindepth, args.maxdepth):
            return result

        dfrac = discordant_fraction(args.bamFileName, args.refFasta, chrom,
                                    max(0, start - pad), end + pad)
        logger.info("%s discordant fraction: %f" % (mutid, dfrac))
        if dfrac > args.maxdfrac:
            logger.warning("%s discordant fraction %f > %f aborting mutation!" % (mutid, dfrac, args.maxdfrac))
            return result

        # ---- build the mutated haplotype from reference slices ----
        ins_id = None
        literal_insseq = None
        donor_interior = False

        if kind == 'BND':
            template = svt.build_breakend(
                reffile, chrom, start, req.mate_chrom, req.mate_pos, pad,
                mutid=mutid, allow_n=args.allowN,
                flip_left=req.flip_left, flip_right=req.flip_right)

        elif kind == 'DUP' and args.donorbam is not None:
            # interior copy number comes from real reads instead of simulation
            donor_interior = True
            template = svt.build_dup_junction(reffile, chrom, start, end, pad,
                                              mutid=mutid, allow_n=args.allowN)

        else:
            insseq = ''
            if kind == 'INS':
                insseq, ins_id = resolve_insertion(args, req, mutid)
                # a library entry round-trips by name; an inline sequence has
                # to carry itself
                if ins_id is None:
                    literal_insseq = insseq

            template = svt.build_interval(
                kind, reffile, chrom, start, end, pad, mutid=mutid,
                allow_n=args.allowN, insseq=insseq, tsdlen=int(req.tsdlen),
                ndups=int(req.ndups), ins_motif=req.ins_motif,
                maxlibsize=int(args.maxlibsize))

        logger.info("%s template length %d, replacing %d bp of reference (ratio %.4f)"
                    % (mutid, len(template.seq), template.exclude_length, template.reads_ratio()))

        # ---- collect the reads this replaces ----
        excl_names = set()
        for excl_chrom, excl_start, excl_end in template.exclude:
            for read in get_reads(args.bamFileName, args.refFasta, excl_chrom,
                                  max(0, excl_start), excl_end, vaf, args.salt):
                excl_names.add(read.query_name)

        if not excl_names:
            logger.warning("%s no reads to replace, skipping site." % mutid)
            return result

        nsimreads = int(round(len(excl_names) * template.reads_ratio()))
        if nsimreads < 1:
            logger.warning("%s simulated read count rounds to zero, skipping site." % mutid)
            return result

        readlen = args.readlen
        if readlen is None:
            readlen = svt.sample_read_length(bamfile, chrom, max(0, start - pad), end + pad)
        if not readlen:
            logger.warning("%s could not determine read length, skipping site." % mutid)
            return result

        exclfile = os.path.join(args.tmpdir, '.'.join((mutid, 'exclude', str(uuid4()), 'txt')))
        with open(exclfile, 'w') as exclude:
            for name in sorted(excl_names):
                exclude.write(name + "\n")

        # ---- simulate and realign ----
        outbam_mutsfile = os.path.join(args.tmpdir, '.'.join((mutid, str(uuid4()), "muts.bam")))

        fq1, fq2 = runwgsim(template.seq, int(readlen), float(args.ismean), float(args.issd),
                            args.tmpdir, nsimreads, mutid=mutid,
                            err_rate=float(args.simerr), seed=args.seed)

        remap_fastq(args.aligner, fq1, args.refFasta, outbam_mutsfile, alignopts,
                    fq2=fq2, mutid=mutid, threads=int(args.alignerthreads))

        if not check_min_read_count(outbam_mutsfile, args.refFasta, 0):
            logger.warning("%s outbam %s has no mapped reads!" % (mutid, outbam_mutsfile))
            return result

        if donor_interior:
            prev = outbam_mutsfile
            outbam_mutsfile = add_donor_reads(args, mutid, prev, chrom, start, end, vaf)
            os.remove(prev)
            if os.path.exists(prev + '.bai'):
                os.remove(prev + '.bai')

        # ---- report ----
        if kind == 'BND':
            record = MutationRecord(
                kind='BND', chrom=chrom, pos=start, vaf=vaf, mutid=mutid,
                mate_chrom=req.mate_chrom, mate_pos=req.mate_pos,
                flip_left=req.flip_left, flip_right=req.flip_right)
        else:
            # an insertion reports where it landed inside the requested
            # interval, and how much sequence it added; the others span
            # exactly what was asked for
            pos = template.event_pos if template.event_pos else start
            record = MutationRecord(
                kind=kind, chrom=chrom, pos=pos,
                end=pos if kind == 'INS' else end,
                svlen=template.event_len if kind == 'INS' else end - start,
                vaf=vaf, mutid=mutid, ins_id=ins_id, tsdlen=int(req.tsdlen),
                insseq=literal_insseq, ins_motif=req.ins_motif,
                ndups=int(req.ndups), donor_interior=donor_interior)

        logfile.write('%s\t%s\t%d\t%d\tvaf=%s\ttemplate_len=%d\texcluded=%d\tsimulated=%d\n'
                      % (kind, chrom, start, end, vaf, len(template.seq),
                         len(excl_names), nsimreads))

        logger.info("%s temporary bam: %s" % (mutid, outbam_mutsfile))

        result.bamfile = outbam_mutsfile
        result.excludefile = exclfile
        result.records = [record]

    except svt.TemplateError as e:
        logger.warning('%s %s' % (mutid, e))
    except Exception:
        logger.error('%s failed:' % mutid)
        traceback.print_exc(file=sys.stderr)
    finally:
        logfile.close()
        bamfile.close()

    return result


def run_sv_spikein(args, requests=None):
    '''
    Apply SV requests and return the MutationRecords.

    Returns records rather than writing a VCF so that a run combining both
    engines can emit one truth file.
    '''
    tmpbams = [] # temporary BAMs, each holds the realigned reads for one mutation
    exclfns = [] # 'exclude' files store reads to be removed from the original BAM due to deletions
    records = [] # MutationRecords, become the truth VCF

    if (args.bamFileName.endswith('.bam') and not os.path.exists(args.bamFileName + '.bai')) or \
        (args.bamFileName.endswith('.cram') and not os.path.exists(args.bamFileName + '.crai')):
        logger.error("input file must be indexed, not .bai or .crai file found for %s" % args.bamFileName)
        sys.exit(1)

    alignopts = {}
    if args.alignopts is not None:
        alignopts = dict([o.split(':') for o in args.alignopts.split(',')])

    # Must happen before any pool.submit: workers inherit this by pickle, and
    # deriving it per-process is what made read selection vary with --procs.
    args.salt = make_salt(args.seed)
    logger.info("read selection salt: %s" % args.salt)

    if args.minctglen is not None:
        logger.warning('--minctglen is deprecated and ignored; use --pad to set the flank size')

    # load insertion library if present
    try:
        if args.inslib is not None:
            logger.info("loading insertion library from %s" % args.inslib)
            args.inslib = load_inslib(args.inslib)
    except Exception:
        logger.error("failed to load insertion library %s" % args.inslib)
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)

    results = []
    pool = ProcessPoolExecutor(max_workers=int(args.procs))

    logdir = 'addsv_logs_' + os.path.basename(args.outBamFile)

    for d in (args.tmpdir, logdir):
        if not os.path.exists(d):
            os.mkdir(d)
            logger.info("created directory: %s" % d)

    assert os.path.exists(logdir), "could not create output directory!"
    assert os.path.exists(args.tmpdir), "could not create temporary directory!"

    if requests is None:
        requests = read_variants(args.varFileName, 'sv',
                                 default_vaf=float(args.svfrac),
                                 maxmuts=int(args.maxmuts) if args.maxmuts else 0)

    if not requests:
        logger.error("no SV requests to make")
        sys.exit(1)

    for nmuts, req in enumerate(requests):
        req.index = nmuts
        results.append(pool.submit(makemut, args, req, alignopts))

    ## process the results of mutation jobs
    for result in results:
        res = result.result()

        if res.ok() and os.path.exists(res.bamfile) and os.path.exists(res.excludefile):
            if check_min_read_count(res.bamfile, args.refFasta, 0):
                tmpbams.append(res.bamfile)
                exclfns.append(res.excludefile)
                records.extend(res.records)
            else:
                os.remove(res.bamfile)
                os.remove(res.excludefile)

    if len(tmpbams) == 0:
        logger.error("no succesful mutations")
        sys.exit(1)

    logger.info("tmpbams: %s" % tmpbams)
    logger.info("exclude: %s" % exclfns)

    excl_merged = 'addsv.exclude.final.' + str(uuid4()) + '.txt'
    mergedtmp = 'addsv.mergetmp.final.' + str(uuid4()) + '.bam'

    logger.info("merging exclude files into %s" % excl_merged)
    with open(excl_merged, 'w') as exclout:
        for exclfn in exclfns:
            with open(exclfn, 'r') as excl:
                for line in excl:
                    exclout.write(line)
            if not args.debug:
                os.remove(exclfn)

    if len(tmpbams) == 1:
        logger.info("only one bam: %s renaming to %s" % (tmpbams[0], mergedtmp))
        os.rename(tmpbams[0], mergedtmp)
    elif len(tmpbams) > 1:
        logger.info("merging bams into %s" % mergedtmp)
        mergebams(tmpbams, mergedtmp, maxopen=int(args.maxopen), debug=args.debug)

    if args.skipmerge:
        logger.info("final merge skipped, please merge manually: %s" % mergedtmp)
        logger.info("exclude file to use: %s" % excl_merged)
    else:
        if args.tagreads:
            from bamsurgeon.markreads import markreads
            tmp_tag_bam = 'tag.%s.bam' % str(uuid4())
            markreads(mergedtmp, args.refFasta, tmp_tag_bam)
            move(tmp_tag_bam, mergedtmp)
            logger.info("tagged reads.")

        logger.info("writing to %s" % args.outBamFile)
        rr.replace_reads(args.bamFileName, mergedtmp, args.outBamFile, args.refFasta, excludefile=excl_merged, allreads=True, keepsecondary=args.keepsecondary, seed=args.seed, quiet=True)
        if not args.debug:
            os.remove(excl_merged)
            os.remove(mergedtmp)

        logger.info("done.")

    if not args.debug:
        for tmpbam in tmpbams:
            if os.path.isfile(tmpbam):
                os.remove(tmpbam)
            if os.path.isfile(tmpbam + '.bai'):
                os.remove(tmpbam + '.bai')

    return records


def main(args):
    logger.info("starting %s called with args: %s" % (sys.argv[0], ' '.join(sys.argv)))

    records = run_sv_spikein(args)

    var_basename = '.'.join(os.path.basename(args.varFileName).split('.')[:-1])
    bam_basename = '.'.join(os.path.basename(args.outBamFile).split('.')[:-1])

    vcf_fn = args.vcf + bam_basename + '.addsv.' + var_basename + '.vcf'

    makevcf.write_vcf(records, args.refFasta, vcf_fn, salt=args.salt)

    logger.info('vcf output written to ' + vcf_fn)

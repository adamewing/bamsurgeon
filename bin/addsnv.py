#!/usr/bin/env python

import sys
import pysam
import argparse
import random
import subprocess
import os
import bamsurgeon.replace_reads as rr
import bamsurgeon.aligners as aligners
import bamsurgeon.mutation as mutation
import bamsurgeon.makevcf as makevcf

from operator import itemgetter
from bamsurgeon.common import *
from uuid import uuid4
from shutil import move
from re import sub
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict as dd

import logging
FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def mut(base, altbase):
    """ change base to something different
    """

    bases = ('A','T','C','G')
    base = base.upper()
    if base not in bases or (altbase is not None and altbase not in ['A','T','G','C']):
        raise ValueError("ERROR base passed to mut(): " + str(base) + " not one of (A,T,C,G)\n")

    if altbase is not None:
        return altbase

    else:
        alt = base
        while alt == base:
            alt = bases[int(random.uniform(0,4))]
        return alt


def makemut(args, hc, avoid, alignopts):
    mutid_list = []
    for site in hc:
        mutid_list.append(site['chrom'] + '_' + str(site['start']) + '_' + str(site['end']) + '_' + str(site['vaf']) + '_' + str(site['altbase']))


    if args.seed is not None: random.seed(int(args.seed) + int(hc[0]['start']))

    bamfile = pysam.AlignmentFile(args.bamFileName)
    bammate = pysam.AlignmentFile(args.bamFileName) # use for mates to avoid iterator problems
    reffile = pysam.Fastafile(args.refFasta)
    tmpbams = []

    #snvfrac = float(args.snvfrac)

    chrom = None
    vaf   = None

    mutpos_list = []
    altbase_list = []
    
    for site in hc:
        if chrom is None:
            chrom = site['chrom']
        else:
            assert chrom == site['chrom'], "haplotype clusters cannot span multiple chromosomes!"

        if vaf is None:
            vaf = site['vaf']
            
        elif vaf != site['vaf']:
            logger.warning("multiple VAFs for single haplotype, using first encountered VAF: %f" % vaf)

        mutpos = int(random.uniform(site['start'],site['end']+1)) # position of mutation in genome
        mutpos_list.append(mutpos)
        altbase_list.append(site['altbase'])

    mutbase_list = []
    refbase_list = []
    mutstr_list  = []

    for n, mutpos in enumerate(mutpos_list):
        refbase = reffile.fetch(chrom,mutpos-1,mutpos)
        altbase = altbase_list[n]
        refbase_list.append(refbase)

        if altbase == refbase.upper() and not args.ignoreref:
            logger.warning("%s specified ALT base matches reference, skipping mutation" % mutid_list[n])
            return None

        try:
            mutbase = mut(refbase, altbase)
            mutbase_list.append(mutbase)

        except ValueError as e:
            logger.warning(mutid_list[n] + " " + ' '.join(("skipped site:",chrom,str(hc[n]['start']),str(hc[n]['end']),"due to N base:",str(e),"\n")))
            return None

        mutstr_list.append(refbase + "-->" + str(mutbase))

    # optional CNV file
    cnv = None
    if (args.cnvfile):
        cnv = pysam.Tabixfile(args.cnvfile, 'r')

    hapstr = "_".join(('haplo',chrom,str(min(mutpos_list)),str(max(mutpos_list))))
    log = open('addsnv_logs_' + os.path.basename(args.outBamFile) + '/' + os.path.basename(args.outBamFile) + "." + hapstr + ".log",'w')

    tmpoutbamname = args.tmpdir + "/" + hapstr + ".tmpbam." + str(uuid4()) + ".bam"
    logger.info("%s creating tmp bam: %s" % (hapstr, tmpoutbamname))
    outbam_muts = pysam.AlignmentFile(tmpoutbamname, 'wb', template=bamfile)

    mutfail, hasSNP, maxfrac, outreads, mutreads, mutmates = mutation.mutate(args, log, bamfile, bammate, chrom, min(mutpos_list), max(mutpos_list)+1, mutpos_list, avoid=avoid, mutid_list=mutid_list, is_snv=True, mutbase_list=mutbase_list, reffile=reffile)

    if mutfail:
        outbam_muts.close()
        os.remove(tmpoutbamname)
        return None

    # pick reads to change
    readlist = []
    for extqname,read in outreads.items():
        if read.seq != mutreads[extqname]:
            readlist.append(extqname)

    logger.info("%s len(readlist): %s" % (hapstr, str(len(readlist))))
    readlist.sort()
    random.shuffle(readlist)

    if len(readlist) < int(args.mindepth):
        logger.warning("%s too few reads in region (%s) skipping..." % (hapstr, str(len(readlist))))
        outbam_muts.close()
        os.remove(tmpoutbamname)
        return None

    if vaf is None:
        vaf = float(args.mutfrac) # default minor allele freq if not otherwise specified

    if cnv: # cnv file is present
        if chrom in cnv.contigs:
            for cnregion in cnv.fetch(chrom,min(mutpos_list),max(mutpos_list)+1):
                cn = float(cnregion.strip().split()[3]) # expect chrom,start,end,CN
                logger.info(hapstr + "\t" + ' '.join(("copy number in snp region:",chrom,str(min(mutpos_list)),str(max(mutpos_list)),"=",str(cn))))
                if float(cn) > 0.0:
                    vaf = vaf/float(cn)
                else:
                    vaf = 0.0
                logger.info("%s adjusted VAF: %f" % (hapstr, vaf))
    else:
        logger.info("%s selected VAF: %f" % (hapstr, vaf))

    lastread = int(len(readlist)*vaf)

    # pick at least args.minmutreads if possible
    if lastread < int(args.minmutreads):
        if len(readlist) > int(args.minmutreads):
            lastread = int(args.minmutreads)
            logger.warning("%s forced %d reads." % (hapstr, lastread))
        else:
            logger.warning("%s dropped site with fewer reads than --minmutreads" % hapstr)
            os.remove(tmpoutbamname)
            return None
    elif lastread <= 0:
        logger.warning("%s dropped mutation with 0 reads" % hapstr)
        os.remove(tmpoutbamname)
        return None

    readtrack = dd(list)

    for readname in readlist:
        orig_name, readpos, pairend = readname.split(',')
        readtrack[orig_name].append('%s,%s' % (readpos, pairend))

    usedreads = 0
    newreadlist = []

    for orig_name in readtrack:
        for read_instance in readtrack[orig_name]:
            newreadlist.append(orig_name + ',' + read_instance)
            usedreads += 1

        if usedreads >= lastread:
            break

    readlist = newreadlist

    logger.info("%s picked: %d" % (hapstr, len(readlist)))

    wrote = 0
    nmut = 0
    mut_out = {}
    # change reads from .bam to mutated sequences
    for extqname,read in outreads.items():
        if read.seq != mutreads[extqname]:
            if not args.nomut and extqname in readlist:
                qual = read.qual # changing seq resets qual (see pysam API docs)
                read.seq = mutreads[extqname] # make mutation
                read.qual = qual
                nmut += 1
        if (not hasSNP) or args.force:
            wrote += 1
            mut_out[extqname] = read

    muts_written = {}

    for extqname in mut_out:
        if extqname not in muts_written:
            outbam_muts.write(mut_out[extqname])
            muts_written[extqname] = True

            if mutmates[extqname] is not None:
                # is mate also in mutated list?
                mate_read = mutmates[extqname]

                pairname = 'F' # read is first in pair
                if mate_read.is_read2:
                    pairname = 'S' # read is second in pair
                if not mate_read.is_paired:
                    pairname = 'U' # read is unpaired

                mateqname = ','.join((mate_read.qname,str(mate_read.pos),pairname))

                if mateqname in mut_out:
                    # yes: output mutated mate
                    outbam_muts.write(mut_out[mateqname])
                    muts_written[mateqname] = True

                else:
                    # no: output original mate
                    outbam_muts.write(mate_read)

    logger.info("%s wrote: %d, mutated: %d" % (hapstr,wrote,nmut))

    if not hasSNP or args.force:
        outbam_muts.close()

        aligners.remap_bam(args.aligner, tmpoutbamname, args.refFasta, alignopts, threads=int(args.alignerthreads), mutid=hapstr, paired=(not args.single), picardjar=args.picardjar, insane=args.insane)

        outbam_muts = pysam.AlignmentFile(tmpoutbamname)
        coverwindow = 1

        avgincover = get_avg_coverage(bamfile, chrom,min(mutpos_list)-coverwindow,max(mutpos_list)+coverwindow)
        avgoutcover = get_avg_coverage(outbam_muts, chrom,min(mutpos_list)-coverwindow,max(mutpos_list)+coverwindow)

        logger.info("%s avgincover: %f, avgoutcover: %f" % (hapstr, avgincover, avgoutcover))

        spikein_snvfrac = 0.0
        if wrote > 0:
            spikein_snvfrac = float(nmut)/float(wrote)

        # qc cutoff for final snv depth 
        if (avgoutcover > 0 and avgincover > 0 and avgoutcover/avgincover >= float(args.coverdiff)) or args.force:
            tmpbams.append(tmpoutbamname)
            for n,site in enumerate(hc):
                snvstr = chrom + ":" + str(site['start']) + "-" + str(site['end']) + " (VAF=" + str(vaf) + ")"
                log.write("\t".join(("snv",snvstr,str(mutpos_list[n]),mutstr_list[n],str(avgincover),str(avgoutcover),str(spikein_snvfrac),str(maxfrac)))+"\n")
        else:

            outbam_muts.close()
            os.remove(tmpoutbamname)
            if os.path.exists(tmpoutbamname + '.bai'):
                os.remove(tmpoutbamname + '.bai')
            logger.warning("%s dropped for outcover/incover < %s" % (hapstr, str(args.coverdiff)))
            return None

    outbam_muts.close()
    bamfile.close()
    bammate.close()
    log.close() 

    return tmpbams

def main(args):
    logger.info("starting %s called with args: %s" % (sys.argv[0], ' '.join(sys.argv)))
    bedfile = open(args.varFileName, 'r')
    reffile = pysam.Fastafile(args.refFasta)

    if (args.bamFileName.endswith('.bam') and not os.path.exists(args.bamFileName + '.bai')) or \
        (args.bamFileName.endswith('.cram') and not os.path.exists(args.bamFileName + '.crai')):
        logger.error("input file must be indexed, not .bai or .crai file found for %s" % args.bamFileName)
        sys.exit(1)

    alignopts = {}
    if args.alignopts is not None:
        alignopts = dict([o.split(':') for o in args.alignopts.split(',')])

    aligners.checkoptions(args.aligner, alignopts, args.picardjar)

    # load readlist to avoid, if specified
    avoid = None
    if args.avoidreads is not None:
        avoid = dictlist(args.avoidreads)

    # make a temporary file to hold mutated reads
    outbam_mutsfile = "addsnv." + str(uuid4()) + ".muts.bam"
    bamfile = pysam.AlignmentFile(args.bamFileName)
    outbam_muts = pysam.AlignmentFile(outbam_mutsfile, 'wb', template=bamfile)
    outbam_muts.close()
    bamfile.close()
    tmpbams = []

    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)
        logger.info("created tmp directory: %s" % args.tmpdir)

    if not os.path.exists('addsnv_logs_' + os.path.basename(args.outBamFile)):
        os.mkdir('addsnv_logs_' + os.path.basename(args.outBamFile))
        logger.info("created directory: addsnv_logs_%s" % os.path.basename(args.outBamFile))

    assert os.path.exists('addsnv_logs_' + os.path.basename(args.outBamFile)), "could not create output directory!"
    assert os.path.exists(args.tmpdir), "could not create temporary directory!"

    pool = ProcessPoolExecutor(max_workers=int(args.procs))
    results = []

    ntried = 0

    targets = []
    for bedline in bedfile:
        if ntried < int(args.numsnvs) or int(args.numsnvs) == 0:
            c = bedline.strip().split()
            target = {
            'chrom'   : c[0],
            'start'   : int(c[1]),
            'end'     : int(c[2]),
            'vaf'     : None,
            'altbase' : None
            }

            # VAF is 4th column, if present
            if len(c) > 3:
                target['vaf'] = float(c[3])

            # ALT is 5th column, if present
            if len(c) == 5:
                altbase = c[4].upper()
                assert altbase in ['A','T','C','G'], "ERROR: ALT " + altbase + " not A, T, C, or G!\n"
                target['altbase'] = altbase

            targets.append(target)
            ntried += 1

    targets = sorted(targets, key=itemgetter('chrom', 'start')) # sort list of dicts by chrom, start

    haploclusters = []

    hc = []
    lastchrom = None
    laststart = None

    hapsize = int(args.haplosize)
    for target in targets:
        if lastchrom is None:
            lastchrom = target['chrom']
            laststart = target['start']
            hc.append(target)

        elif target['chrom'] == lastchrom:
            if laststart is None:
                laststart = target['start']
                hc.append(target)
            elif target['start'] - laststart < hapsize:
                hc.append(target)
            else:
                haploclusters.append(hc)
                hc = []
                hc.append(target)

        elif target['chrom'] != lastchrom:
            haploclusters.append(hc)
            hc = []
            laststart = None
            hc.append(target)

    haploclusters.append(hc)

    for hc in haploclusters:
        # make mutation (submit job to thread pool)
        result = pool.submit(makemut, args, hc, avoid, alignopts)
        results.append(result)


    for result in results:
        tmpbamlist = result.result()
        if tmpbamlist is not None:
            for tmpbam in tmpbamlist:
                if os.path.exists(tmpbam):
                    tmpbams.append(tmpbam)

    if len(tmpbams) == 0:
        logger.error("no succesful mutations")
        sys.exit(1)        

    # merge tmp bams
    if len(tmpbams) == 1:
        move(tmpbams[0],outbam_mutsfile)
    elif len(tmpbams) > 1:
        mergebams(tmpbams,outbam_mutsfile,maxopen=int(args.maxopen))

    bedfile.close()

    # cleanup
    for bam in tmpbams:
        if os.path.exists(bam):
            os.remove(bam)
        if os.path.exists(bam + '.bai'):
            os.remove(bam + '.bai')

    if args.skipmerge:
        logger.info("skipping merge, plase merge reads from %s manually." % outbam_mutsfile)
    else:
        if args.tagreads:
            from bamsurgeon.markreads import markreads
            tmp_tag_bam = 'tag.%s.bam' % str(uuid4())
            markreads(outbam_mutsfile, tmp_tag_bam)
            move(tmp_tag_bam, outbam_mutsfile)
            logger.info("tagged reads.")

        logger.info("done making mutations, merging mutations into %s --> %s" % (args.bamFileName, args.outBamFile))
        rr.replace_reads(args.bamFileName, outbam_mutsfile, args.outBamFile, keepqual=True, seed=args.seed)

        #cleanup
        os.remove(outbam_mutsfile)

    var_basename = '.'.join(os.path.basename(args.varFileName).split('.')[:-1])
    bam_basename = '.'.join(os.path.basename(args.outBamFile).split('.')[:-1])

    vcf_fn = bam_basename + '.addsnv.' + var_basename + '.vcf'

    makevcf.write_vcf_snv('addsnv_logs_' + os.path.basename(args.outBamFile), args.refFasta, vcf_fn)

    logger.info('vcf output written to ' + vcf_fn)


def run():
    parser = argparse.ArgumentParser(description='adds SNVs to reads, outputs modified reads as .bam along with mates')
    parser.add_argument('-v', '--varfile', dest='varFileName', required=True, help='Target regions to try and add a SNV, as BED')
    parser.add_argument('-f', '--bamfile', dest='bamFileName', required=True, help='sam/bam file from which to obtain reads')
    parser.add_argument('-r', '--reference', dest='refFasta', required=True, help='reference genome, fasta indexed with bwa index _and_ samtools faidx')
    parser.add_argument('-o', '--outbam', dest='outBamFile', required=True, help='.bam file name for output')
    parser.add_argument('-s', '--snvfrac', dest='snvfrac', default=1, help='maximum allowable linked SNP MAF (for avoiding haplotypes) (default = 1)')
    parser.add_argument('-m', '--mutfrac', dest='mutfrac', default=0.5, help='allelic fraction at which to make SNVs (default = 0.5)')
    parser.add_argument('-n', '--numsnvs', dest='numsnvs', default=0, help="maximum number of mutations to try (default: entire input)")
    parser.add_argument('-c', '--cnvfile', dest='cnvfile', default=None, help="tabix-indexed list of genome-wide absolute copy number values (e.g. 2 alleles = no change)")
    parser.add_argument('-d', '--coverdiff', dest='coverdiff', default=0.9, help="allow difference in input and output coverage (default=0.9)")
    parser.add_argument('-z', '--haplosize', default=0, help='haplotype size (default = 0)')
    parser.add_argument('-p', '--procs', dest='procs', default=1, help="split into multiple processes (default=1)")
    parser.add_argument('--picardjar', default=None, help='path to picard.jar, required for most aligners')
    parser.add_argument('--mindepth', default=10, help='minimum read depth to make mutation (default = 10)')
    parser.add_argument('--maxdepth', default=2000, help='maximum read depth to make mutation (default = 2000)')
    parser.add_argument('--minmutreads', default=3, help='minimum number of mutated reads to output per site')
    parser.add_argument('--avoidreads', default=None, help='file of read names to avoid (mutations will be skipped if overlap)')
    parser.add_argument('--nomut', action='store_true', default=False, help="dry run")
    parser.add_argument('--ignoresnps', action='store_true', default=False, help="make mutations even if there are non-reference alleles sharing the relevant reads")
    parser.add_argument('--ignoreref', action='store_true', default=False, help="make mutations even if the mutation is back to the reference allele")
    parser.add_argument('--force', action='store_true', default=False, help="force mutation to happen regardless of nearby SNP or low coverage")
    parser.add_argument('--insane', action='store_true', default=False, help="ignore sanity check enforcing input read count = output read count in realignment")
    parser.add_argument('--single', action='store_true', default=False, help="input BAM is simgle-ended (default is paired-end)")
    parser.add_argument('--maxopen', dest='maxopen', default=1000, help="maximum number of open files during merge (default 1000)")
    parser.add_argument('--requirepaired', action='store_true', default=False, help='skip mutations if unpaired reads are present')
    parser.add_argument('--tagreads', action='store_true', default=False, help='add BS tag to altered reads')
    parser.add_argument('--skipmerge', action='store_true', default=False, help="final output is tmp file to be merged")
    parser.add_argument('--ignorepileup', action='store_true', default=False, help="do not check pileup depth in mutation regions")
    parser.add_argument('--aligner', default='backtrack', help='supported aligners: ' + ','.join(aligners.supported_aligners_bam))
    parser.add_argument('--alignerthreads', default=1, help='threads used per realignment (default = 1)')
    parser.add_argument('--alignopts', default=None, help='aligner-specific options as comma delimited list of option1:value1,option2:value2,...')
    parser.add_argument('--tmpdir', default='addsnv.tmp', help='temporary directory (default=addsnv.tmp)')
    parser.add_argument('--seed', default=None, help='seed random number generation')
    args = parser.parse_args()

    if 'BAMSURGEON_PICARD_JAR' in os.environ:
        args.picardjar = os.environ['BAMSURGEON_PICARD_JAR']

    main(args)

if __name__ == '__main__':
    run()

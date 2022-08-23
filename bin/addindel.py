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

from bamsurgeon.common import *
from uuid import uuid4
from re import sub
from shutil import move
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict as dd

import logging
FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_mutstr(chrom, start, end, ins, ref):
    return 'FIX get_mutstr'


def makemut(args, chrom, start, end, vaf, ins, avoid, alignopts):
    ''' is ins is a sequence, it will is inserted at start, otherwise delete from start to end'''

    if args.seed is not None: random.seed(int(args.seed) + int(start))

    mutid = chrom + '_' + str(start) + '_' + str(end) + '_' + str(vaf)
    if ins is None:
        mutid += ':DEL'
    else:
        mutid += ':INS:' + ins


    bamfile = pysam.AlignmentFile(args.bamFileName)
    bammate = pysam.AlignmentFile(args.bamFileName) # use for mates to avoid iterator problems
    reffile = pysam.Fastafile(args.refFasta)
    tmpbams = []

    is_insertion = ins is not None
    is_deletion  = ins is None

    snvfrac = float(args.snvfrac)

    mutstr = get_mutstr(chrom, start, end, ins, reffile)

    del_ln = 0
    if is_deletion:
        del_ln = end-start

    mutpos = start
    mutpos_list = [start]

    # optional CNV file
    cnv = None
    if (args.cnvfile):
        cnv = pysam.Tabixfile(args.cnvfile, 'r')

    log = open('addindel_logs_' + os.path.basename(args.outBamFile) + '/' + os.path.basename(args.outBamFile) + "." + "_".join((chrom,str(start),str(end))) + ".log",'w')

    tmpoutbamname = args.tmpdir + "/" + mutid + ".tmpbam." + str(uuid4()) + ".bam"
    logger.info("%s creating tmp bam: %s" % (mutid ,tmpoutbamname))
    outbam_muts = pysam.AlignmentFile(tmpoutbamname, 'wb', template=bamfile)

    mutfail, hasSNP, maxfrac, outreads, mutreads, mutmates = mutation.mutate(args, log, bamfile, bammate, chrom, mutpos, mutpos+del_ln+1, mutpos_list, avoid=avoid, mutid_list=[mutid], is_insertion=is_insertion, is_deletion=is_deletion, ins_seq=ins, reffile=reffile, indel_start=start, indel_end=end)

    if mutfail:
        outbam_muts.close()
        os.remove(tmpoutbamname)
        return None

    # pick reads to change
    readlist = []
    for extqname,read in outreads.items():
        if read.seq != mutreads[extqname]:
            readlist.append(extqname)

    logger.info("%s len(readlist): %d" % (mutid, len(readlist)))
    readlist.sort()
    random.shuffle(readlist)

    if len(readlist) < int(args.mindepth):
        logger.warning("%s skipped, too few reads in region: %d" % (mutid, len(readlist)))
        outbam_muts.close()
        os.remove(tmpoutbamname)
        return None

    if vaf is None:
        vaf = float(args.mutfrac) # default minor allele freq if not otherwise specified
        
    if cnv: # cnv file is present
        if chrom in cnv.contigs:
            for cnregion in cnv.fetch(chrom,start,end):
                cn = float(cnregion.strip().split()[3]) # expect chrom,start,end,CN
                logger.info(mutid + "\t" + ' '.join(("copy number in snp region:",chrom,str(start),str(end),"=",str(cn))))
                if float(cn) > 0.0:
                    vaf = vaf/float(cn)
                else:
                    vaf = 0.0
                logger.info("%s adjusted VAF: %f" % (mutid, vaf))
    else:
        logger.info("%s selected VAF: %f" % (mutid, vaf))

    lastread = int(len(readlist)*vaf)

    # pick at least args.minmutreads if possible
    if lastread < int(args.minmutreads):
        if len(readlist) > int(args.minmutreads):
            lastread = int(args.minmutreads)
            logger.warning("%s forced %d reads" % (mutid, lastread))
        else:
            logger.warning("%s dropped site with fewer reads than --minmutreads" % mutid)
            os.remove(tmpoutbamname)
            return None
    elif lastread <= 0:
        logger.warning("%s dropped mutation with 0 reads" % mutid)
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

    logger.info("%s picked: %d reads" % (mutid, len(readlist)))

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
        if not hasSNP or args.force:
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

    logger.info("%s wrote: %d, mutated: %d" % (mutid,wrote,nmut))

    if not hasSNP or args.force:
        outbam_muts.close()
        aligners.remap_bam(args.aligner, tmpoutbamname, args.refFasta, alignopts, threads=int(args.alignerthreads), mutid=mutid, paired=(not args.single), picardjar=args.picardjar, insane=args.insane)

        outbam_muts = pysam.AlignmentFile(tmpoutbamname)
        coverwindow = 1
        avgincover = get_avg_coverage(bamfile, chrom,mutpos-coverwindow,mutpos+del_ln+coverwindow)
        avgoutcover = get_avg_coverage(outbam_muts, chrom,mutpos-coverwindow,mutpos+del_ln+coverwindow)
        spikein_frac = 0.0
        if wrote > 0:
            spikein_frac = float(nmut)/float(wrote)

        # qc cutoff for final snv depth 
        if (avgoutcover > 0 and avgincover > 0 and avgoutcover/avgincover >= float(args.coverdiff)) or args.force:
            tmpbams.append(tmpoutbamname)
            indelstr = ''
            if is_insertion:
                indelstr = ':'.join(('INS', chrom, str(start), ins))
            else:
                indelstr = ':'.join(('DEL', chrom, str(start), str(end)))

            snvstr = chrom + ":" + str(start) + "-" + str(end) + " (VAF=" + str(vaf) + ")"
            log.write("\t".join(("indel",indelstr,str(mutpos),mutstr,str(avgincover),str(avgoutcover),str(spikein_frac),str(maxfrac)))+"\n")
        else:
            outbam_muts.close()
            os.remove(tmpoutbamname)
            if os.path.exists(tmpoutbamname + '.bai'):
                os.remove(tmpoutbamname + '.bai')
                
            logger.warning("%s dropped for outcover/incover < %s" % (mutid, str(args.coverdiff)))
            return None

    outbam_muts.close()
    bamfile.close()
    bammate.close()
    log.close() 

    return sorted(tmpbams)


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
    outbam_mutsfile = "addindel." + str(uuid4()) + ".muts.bam"
    bamfile = pysam.AlignmentFile(args.bamFileName)
    outbam_muts = pysam.AlignmentFile(outbam_mutsfile, 'wb', template=bamfile)
    outbam_muts.close()
    bamfile.close()
    tmpbams = []

    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)
        logger.info("created tmp directory: %s" % args.tmpdir)

    if not os.path.exists('addindel_logs_' + os.path.basename(args.outBamFile)):
        os.mkdir('addindel_logs_' + os.path.basename(args.outBamFile))
        logger.info("created directory: addindel_logs_%s" % os.path.basename(args.outBamFile))

    assert os.path.exists('addindel_logs_' + os.path.basename(args.outBamFile)), "could not create output directory!"
    assert os.path.exists(args.tmpdir), "could not create temporary directory!"

    pool = ProcessPoolExecutor(max_workers=int(args.procs))
    results = []

    ntried = 0
    for bedline in bedfile:
        if ntried < int(args.numsnvs) or int(args.numsnvs) == 0:
            c = bedline.strip().split()
            chrom = c[0]
            start = int(c[1])
            end   = int(c[2])
            vaf   = float(c[3])
            type  = c[4]
            ins   = None

            assert type in ('INS', 'DEL')
            if type == 'INS':
                ins = c[5]

            # make mutation (submit job to thread pool)
            result = pool.submit(makemut, args, chrom, start, end, vaf, ins, avoid, alignopts)
            results.append(result)
            ntried += 1

    for result in results:
        tmpbamlist = result.result()
        if tmpbamlist is not None:
            for tmpbam in tmpbamlist:
                if os.path.exists(tmpbam):
                    tmpbams.append(tmpbam)


    if len(tmpbams) == 0:
        logger.error("no succesful mutations")
        sys.exit(1)

    tmpbams.sort()

    # merge tmp bams
    if len(tmpbams) == 1:
        os.rename(tmpbams[0],outbam_mutsfile)
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

    vcf_fn = bam_basename + '.addindel.' + var_basename + '.vcf'

    makevcf.write_vcf_indel('addindel_logs_' + os.path.basename(args.outBamFile), args.refFasta, vcf_fn)

    logger.info('vcf output written to ' + vcf_fn)
    
def run():
    # run this script
    parser = argparse.ArgumentParser(description='adds INDELs to reads, outputs modified reads as .bam along with mates')
    parser.add_argument('-v', '--varfile', dest='varFileName', required=True, help='Target regions to try and add a SNV, as BED')
    parser.add_argument('-f', '--bamfile', dest='bamFileName', required=True, help='sam/bam file from which to obtain reads')
    parser.add_argument('-r', '--reference', dest='refFasta', required=True, help='reference genome, fasta indexed with bwa index _and_ samtools faidx')
    parser.add_argument('-o', '--outbam', dest='outBamFile', required=True, help='.bam file name for output')
    parser.add_argument('-s', '--snvfrac', dest='snvfrac', default=1, help='maximum allowable linked SNP MAF (for avoiding haplotypes) (default = 1)')
    parser.add_argument('-m', '--mutfrac', dest='mutfrac', default=0.5, help='allelic fraction at which to make SNVs (default = 0.5)')
    parser.add_argument('-n', '--numsnvs', dest='numsnvs', default=0, help="maximum number of mutations to try (default: entire input)")
    parser.add_argument('-c', '--cnvfile', dest='cnvfile', default=None, help="tabix-indexed list of genome-wide absolute copy number values (e.g. 2 alleles = no change)")
    parser.add_argument('-d', '--coverdiff', dest='coverdiff', default=0.1, help="allow difference in input and output coverage (default=0.1)")
    parser.add_argument('-p', '--procs', dest='procs', default=1, help="split into multiple processes (default=1)")
    parser.add_argument('--picardjar', default=None, help='path to picard.jar')
    parser.add_argument('--mindepth', default=10, help='minimum read depth to make mutation (default = 10)')
    parser.add_argument('--maxdepth', default=2000, help='maximum read depth to make mutation (default = 2000)')
    parser.add_argument('--minmutreads', default=3, help='minimum number of mutated reads to output per site')
    parser.add_argument('--avoidreads', default=None, help='file of read names to avoid (mutations will be skipped if overlap)')
    parser.add_argument('--nomut', action='store_true', default=False, help="dry run")
    parser.add_argument('--det', action='store_true', default=False, help="deterministic base changes: make transitions only")
    parser.add_argument('--force', action='store_true', default=False, help="force mutation to happen regardless of nearby SNP or low coverage")
    parser.add_argument('--insane', action='store_true', default=False, help="ignore sanity check enforcing input read count = output read count in realignment")
    parser.add_argument('--single', action='store_true', default=False, help="input BAM is simgle-ended (default is paired-end)")
    parser.add_argument('--maxopen', dest='maxopen', default=1000, help="maximum number of open files during merge (default 1000)")
    parser.add_argument('--requirepaired', action='store_true', default=False, help='skip mutations if unpaired reads are present')
    parser.add_argument('--aligner', default='backtrack', help='supported aligners: ' + ','.join(aligners.supported_aligners_bam))
    parser.add_argument('--alignopts', default=None, help='aligner-specific options as comma delimited list of option1:value1,option2:value2,...')
    parser.add_argument('--alignerthreads', default=1, help='threads used per realignment (default = 1)')
    parser.add_argument('--tagreads', action='store_true', default=False, help='add BS tag to altered reads')
    parser.add_argument('--skipmerge', action='store_true', default=False, help="final output is tmp file to be merged")
    parser.add_argument('--ignorepileup', action='store_true', default=False, help="do not check pileup depth in mutation regions")
    parser.add_argument('--tmpdir', default='addindel.tmp', help='temporary directory (default=addindel.tmp)')
    parser.add_argument('--seed', default=None, help='seed random number generation')
    args = parser.parse_args()

    if 'BAMSURGEON_PICARD_JAR' in os.environ:
        args.picardjar = os.environ['BAMSURGEON_PICARD_JAR']

    main(args)

if __name__ == '__main__':
    run()

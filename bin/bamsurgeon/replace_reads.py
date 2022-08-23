#!/usr/bin/env python

import sys
import pysam
import argparse
import random
from collections import defaultdict

from bamsurgeon.common import rc

import logging
FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def cleanup(read,orig,RG):
    '''
    fixes unmapped reads that are marked as 'reverse'
    fill in read group at random from existing RGs if 
    RG tags are present in .bam header 
    '''

    if read.is_unmapped and read.is_reverse:
        read.is_reverse = False
        qual = read.qual
        read.seq  = rc(read.seq)
        read.qual = qual[::-1]
    if read.mate_is_unmapped and read.mate_is_reverse:
        read.mate_is_reverse = False
        # mate seq/qual should be caught by the above logic

    if RG:
        hasRG = False
        if read.tags is not None:
            for tag in read.tags:
                if tag[0] == 'RG':
                    hasRG = True

        # use RG from original read if it exists
        if orig is not None:
            if not hasRG and orig.tags is not None:
                if orig.has_tag('RG'):
                    read.set_tag('RG', orig.get_tag('RG'), "Z")
                    hasRG = True

        if not hasRG:
            # give up and add random read group from list in header (e.g. for simulated reads)
            newRG = RG[random.randint(0,len(RG)-1)]
            read.set_tag('RG', newRG, "Z")
            
    return read

def get_RGs(bam):
    '''return list of RG IDs'''
    RG = []
    if 'RG' in bam.header:
        for headRG in bam.header['RG']:
            RG.append(headRG['ID'])
    return RG

def get_excluded_reads(file):
    '''read list of excluded reads into a set'''
    ex = set()
    with open(file) as f:
        for line in f:
            ex.add(line.strip())
    return ex

def compare_ref(targetbam, donorbam):
    ''' if targetbam and donorbam are aligned to different references 
        and the references are in a different order it's a problem
    '''
    for ref in targetbam.references:
        if ref not in donorbam.references or donorbam.gettid(ref) != targetbam.gettid(ref):
            logger.error("contig mismatch: %s\n" % ref)
            return False
    return True
    

def replace_reads(origbamfile, mutbamfile, outbamfile, nameprefix=None, excludefile=None, allreads=False, keepqual=False, progress=False, keepsecondary=False, keepsupplementary=False, seed=None, quiet=False):
    ''' outputbam must be writeable and use targetbam as template
        read names in excludefile will not appear in final output
    '''
    targetbam = pysam.AlignmentFile(origbamfile)
    donorbam  = pysam.AlignmentFile(mutbamfile)
    write_mode = 'wc' if outbamfile.endswith('.cram') else 'wb'
    outputbam  = pysam.AlignmentFile(outbamfile, write_mode, template=targetbam)

    if seed is not None: random.seed(int(seed))

    # check whether references are compatible
    if not compare_ref(targetbam, donorbam):
        sys.exit("Target and donor are aligned to incompatable reference genomes!")

    RG = get_RGs(targetbam) # read groups

    exclude = {}
    if excludefile:
        exclude = get_excluded_reads(excludefile)

    # load reads from donorbam into dict 
    logger.info("loading donor reads into dictionary...\n")

    rdict = defaultdict(lambda: [None, None, None])
    secondary = defaultdict(list) # track secondary alignments, if specified
    supplementary = defaultdict(list) # track supplementary alignments, if specified    
    excount = 0 # number of excluded reads
    nullcount = 0 # number of null reads
    nr = 0
    for read in donorbam.fetch(until_eof=True):
        if read.seq is None:
            # skip unmapped reads
            nullcount += 1
            continue
        if read.qname in exclude:
            excount += 1
            continue
        if nameprefix:
            qual = read.qual # temp
            read.qname = nameprefix + read.qname # must set name _before_ setting quality (see pysam docs)
            read.qual = qual
        extqname = read.qname
        if not read.is_secondary and not read.is_supplementary:
            rlist = rdict[extqname]
            # 0: first pair, 1: second pair, 2: unpaired
            rlist[0 if read.is_read1 else 1 if read.is_read2 else 2] = read
            nr += 1
        elif keepsecondary and read.is_secondary:
            secondary[extqname].append(read)
        elif keepsupplementary and read.is_supplementary:
            supplementary[extqname].append(read)
        else:
            excount += 1

    logger.info('secondary reads count:'+ str(sum([len(v) for k,v in secondary.items()])))
    logger.info('supplementary reads count:'+ str(sum([len(v) for k,v in supplementary.items()])))
    logger.info("loaded " + str(nr) + " reads, (" + str(excount) + " excluded, " + str(nullcount) + " null or secondary or supplementary--> ignored)\n")
    excount = 0
    recount = 0 # number of replaced reads
    used = set() # set of used reads
    prog = 0
    ignored_target = 0 # number of supplemental / secondary reads in original

    for read in targetbam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary:
            ignored_target += 1
            continue
        prog += 1
        if progress and prog % 10000000 == 0:
            logger.info("processed " + str(prog) + " reads.\n")
        if read.qname in exclude:
            excount += 1
            continue
        if nameprefix:
            qual = read.qual # temp
            read.qname = nameprefix + read.qname
            read.qual = qual
        extqname = read.qname
        newReads = []
        if extqname in rdict:
            newRead = rdict[extqname][0 if read.is_read1 else 1 if read.is_read2 else 2]
            if keepqual:
                newRead.qual = read.qual
            newReads = [newRead]
            used.add(extqname)
            recount += 1
        if keepsecondary and extqname in secondary:
                newReads.extend(secondary[extqname])
                used.add(extqname)
                recount += len(secondary[extqname])
        if keepsupplementary and extqname in supplementary:
                newReads.extend(supplementary[extqname])
                used.add(extqname)
                recount += len(supplementary[extqname])
        #non of the above, then write the original read back
        elif len(newReads) == 0:
            newReads = [read]
        for newRead in newReads:
            newRead = cleanup(newRead,read,RG)
            outputbam.write(newRead)
            
    if not quiet:
        logger.info("replaced " + str(recount) + " reads (" + str(excount) + " excluded )\n")
        logger.info("kept " + str(sum([len(v) for k,v in secondary.items()])) + " secondary reads.\n")
        logger.info("kept " + str(sum([len(v) for k,v in supplementary.items()])) + " supplementary reads.\n") 
        logger.info("ignored %d non-primary reads in target BAM.\n" % ignored_target) 

    nadded = 0
    # dump the unused reads from the donor if requested with --all
    if allreads:
        for extqname in rdict.keys():
            if extqname not in used and extqname not in exclude:
                for read in rdict[extqname]:
                    if read is None: continue
                    read = cleanup(read,None,RG)
                    outputbam.write(read)
                    nadded += 1
        logger.info("added " + str(nadded) + " reads due to --all\n")

    targetbam.close()
    donorbam.close()
    outputbam.close()


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='replaces aligned reads in bamfile1 with aligned reads from bamfile2')
    parser.add_argument('-b', '--bam', dest='targetbam', required=True, help='original .bam')
    parser.add_argument('-r', '--replacebam', dest='donorbam', required=True, help='.bam with reads to replace original bam')
    parser.add_argument('-o', '--outputbam', dest='outputbam', required=True, help="name for new .bam output")
    parser.add_argument('-n', '--namechange', dest='namechange', default=None, help="change all read names by prepending string (passed as -n [string])")
    parser.add_argument('-x', '--exclude', dest='exclfile', default=None, help="file containing a list of read names to ignore (exclude from output)")
    parser.add_argument('--all', action='store_true', default=False, help="append reads that don't match target .bam")
    parser.add_argument('--keepqual', action='store_true', default=False, help="keep original quality scores, replace read and mapping only for primary reads")
    parser.add_argument('--progress', action='store_true', default=False, help="output progress every 10M reads")
    parser.add_argument('--keepsecondary', action='store_true', default=False, help='keep secondary reads in final BAM')
    parser.add_argument('--keepsupplementary', action='store_true', default=False, help='keep supplementary reads in final BAM')    
    args = parser.parse_args()

    replace_reads(args.targetbam, args.donorbam, args.outputbam, args.namechange, args.exclfile, args.all, args.keepqual, args.progress, args.keepsecondary,args.keepsupplementary)

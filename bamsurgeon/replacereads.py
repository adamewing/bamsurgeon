#!/usr/bin/env python

import sys
import pysam
import argparse
import random
from collections import defaultdict
from string import maketrans

def rc(dna):
    ''' reverse complement '''
    complements = maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]

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
                for tag in orig.tags:
                    if tag[0] == 'RG':
                        read.tags = read.tags + [tag] 
                        hasRG = True

        if not hasRG:
            # give up and add random read group from list in header (e.g. for simulated reads)
            newRG = RG[random.randint(0,len(RG)-1)]
            read.tags = read.tags + [("RG",newRG)]
    return read

def getRGs(bam):
    '''return list of RG IDs'''
    RG = []
    if 'RG' in bam.header:
        for headRG in bam.header['RG']:
            RG.append(headRG['ID'])
    return RG

def getExcludedReads(file):
    '''read list of excluded reads into a dictionary'''
    ex = {}
    f = open(file,'r')
    for line in f:
        line = line.strip()
        ex[line] = True
    f.close()
    return ex

def compare_ref(targetbam, donorbam):
    ''' if targetbam and donorbam are aligned to different references 
        and the references are in a different order it's a problem
    '''
    for ref in targetbam.references:
        if ref not in donorbam.references or donorbam.gettid(ref) != targetbam.gettid(ref):
            sys.stderr.write("contig mismatch: %s\n" % ref)
            return False
    return True
    

def replaceReads(targetbam, donorbam, outputbam, nameprefix=None, excludefile=None, allreads=False, keepqual=False, progress=False, keepsecondary=False, keepsupplementary=False, seed=None):
    ''' targetbam, donorbam, and outputbam are pysam.Samfile objects
        outputbam must be writeable and use targetbam as template
        read names in excludefile will not appear in final output
    '''

    if seed is not None: random.seed(int(seed))

    # check whether references are compatible
    if not compare_ref(targetbam, donorbam):
        sys.exit("Target and donor are aligned to incompatable reference genomes!")

    RG = getRGs(targetbam) # read groups

    exclude = {}
    if excludefile:
        exclude = getExcludedReads(excludefile)

    # load reads from donorbam into dict 
    sys.stdout.write("loading donor reads into dictionary...\n")

    #rdict = defaultdict(list)
    rdict = {}
    secondary = defaultdict(list) # track secondary alignments, if specified
    supplementary = defaultdict(list) # track supplementary alignments, if specified    
    excount = 0 # number of excluded reads
    nullcount = 0 # number of null reads
    nr = 0
    for read in donorbam.fetch(until_eof=True):
        if read.seq is not None:
            if read.qname not in exclude:
                pairname = 'F' # read is first in pair
                if read.is_read2:
                    pairname = 'S' # read is second in pair
                if not read.is_paired:
                    pairname = 'U' # read is unpaired
                if nameprefix:
                    qual = read.qual # temp
                    read.qname = nameprefix + read.qname # must set name _before_ setting quality (see pysam docs)
                    read.qual = qual
                extqname = ','.join((read.qname,pairname))
                if not read.is_secondary and not read.is_supplementary:
                    rdict[extqname] = read
                    nr += 1
                elif read.is_secondary and keepsecondary:
                    secondary[extqname].append(read)
                elif read.is_supplementary and keepsupplementary:
                    supplementary[extqname].append(read)
            else: # no seq!
                excount += 1
        else:
            nullcount += 1

    print 'secondary reads count:'+ str(sum([len(v) for k,v in secondary.iteritems()]))
    print 'supplementary reads count:'+ str(sum([len(v) for k,v in supplementary.iteritems()]))
    sys.stdout.write("loaded " + str(nr) + " reads, (" + str(excount) + " excluded, " + str(nullcount) + " null or secondary or supplementary--> ignored)\n")
    excount = 0
    recount = 0 # number of replaced reads
    used = {}
    prog = 0
    ignored_target = 0 # number of supplemental / secondary reads in original

    for read in targetbam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary:
            ignored_target += 1
            continue

        prog += 1
        if progress and prog % 10000000 == 0:
            sys.stdout.write("processed " + str(prog) + " reads.\n")

        if read.qname not in exclude:
            pairname = 'F' # read is first in pair
            if read.is_read2:
                pairname = 'S' # read is second in pair
            if not read.is_paired:
                pairname = 'U' # read is unpaired
            if nameprefix:
                qual = read.qual # temp
                read.qname = nameprefix + read.qname
                read.qual = qual

            extqname = ','.join((read.qname,pairname))            
            #check if this read has been processed already. If so, skip to the next read
            if used.get(extqname): continue
            newReads = []
            if extqname in rdict:
                if keepqual:
                    try:
                        rdict[extqname].qual = read.qual
                    except ValueError as e:
                        sys.stdout.write("error replacing quality score for read: " + str(rdict[extqname].qname) + " : " + str(e) + "\n")
                        sys.stdout.write("donor:  " + str(rdict[extqname]) + "\n")
                        sys.stdout.write("target: " + str(read) + "\n")
                        sys.exit(1)
                newReads = [rdict[extqname]]
                used[extqname] = True
                recount += 1
            if extqname in secondary and keepsecondary:
                    newReads.extend(secondary[extqname])
                    used[extqname] = True
                    recount += len(secondary[extqname])
            if extqname in supplementary and keepsupplementary:
                    newReads.extend(supplementary[extqname])
                    used[extqname] = True
                    recount += len(supplementary[extqname])
            #non of the above, then write the original read back
            elif len(newReads) == 0:
                newReads = [read]
            assert(len(newReads) != 0)
            for newRead in newReads:
                newRead = cleanup(newRead,read,RG)
                outputbam.write(newRead)
        else:
            excount += 1
    sys.stdout.write("replaced " + str(recount) + " reads (" + str(excount) + " excluded )\n")
    sys.stdout.write("kept " + str(sum([len(v) for k,v in secondary.iteritems()])) + " secondary reads.\n")
    sys.stdout.write("kept " + str(sum([len(v) for k,v in supplementary.iteritems()])) + " supplementary reads.\n") 
    sys.stdout.write("ignored %d non-primary reads in target BAM.\n" % ignored_target) 

    nadded = 0
    # dump the unused reads from the donor if requested with --all
    if allreads:
        for extqname in rdict.keys():
            if extqname not in used and extqname not in exclude:
                rdict[extqname] = cleanup(rdict[extqname],None,RG)
                outputbam.write(rdict[extqname])
                nadded += 1
        sys.stdout.write("added " + str(nadded) + " reads due to --all\n")

def main(args):
    targetbam = pysam.Samfile(args.targetbam, 'rb')
    donorbam  = pysam.Samfile(args.donorbam, 'rb')
    outputbam = pysam.Samfile(args.outputbam, 'wb', template=targetbam)

    replaceReads(targetbam, donorbam, outputbam, args.namechange, args.exclfile, args.all, args.keepqual, args.progress, args.keepsecondary,args.keepsupplementary)

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
    main(args)

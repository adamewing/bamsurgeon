#!/bin/env python

import sys,pysam,argparse
from random import randint

def cleanup(read,RG):
    '''
    fixes unmapped reads that are marked as 'reverse'
    fill in read group at random from existing RGs if 
    RG tags are present in .bam header 
    '''
    if read.is_unmapped and read.is_reverse:
        read.is_reverse = False

    if RG:
        hasRG = False
        for tag in read.tags:
            if tag[0] == 'RG':
                hasRG = True

        if not hasRG:
            # add random read group from list in header
            newRG = RG[randint(0,len(RG)-1)]
            read.tags = read.tags + [("RG",newRG)]
    return read

def getRGs(bam):
    '''return list of RG IDs'''
    RG = []
    if 'RG' in bam.header:
        for headRG in bam.header['RG']:
            RG.append(headRG['ID'])
    return RG

def main(args):
    targetbam = pysam.Samfile(args.targetbam, 'rb')
    donorbam  = pysam.Samfile(args.donorbam, 'rb')
    outputbam = pysam.Samfile(args.outputbam, 'wb', template=targetbam)

    RG = getRGs(targetbam) # read groups

    namechange = None

    if args.namechange:
        namechange = args.namechange

    # load reads from donorbam into dict 
    sys.stderr.write("loading donor reads into dictionary...\n")
    nr = 0
    rdict = {}
    for read in donorbam.fetch(until_eof=True):
        pairname = 'F' # read is first in pair
        if read.is_read2:
            pairname = 'S' # read is second in pair
        if not read.is_paired:
            pairname = 'U' # read is unpaired
        if namechange:
            qual = read.qual # temp
            read.qname = args.namechange + read.qname # must set name _before_ setting quality (see pysam docs)
            read.qual = qual
        extqname = ','.join((read.qname,pairname))
        rdict[extqname] = read
        nr += 1

    sys.stderr.write("loaded " + str(nr) + " reads, reading " + args.targetbam + "...\n")


    used = {}
    for read in targetbam.fetch(until_eof=True):
        pairname = 'F' # read is first in pair
        if read.is_read2:
            pairname = 'S' # read is second in pair
        if not read.is_paired:
            pairname = 'U' # read is unpaired
        if namechange:
            qual = read.qual # temp
            read.qname = args.namechange + read.qname
            read.qual = qual

        extqname = ','.join((read.qname,pairname))
        if extqname in rdict: 
            if args.keepqual:
                rdict[extqname].qual = read.qual
            rdict[extqname] = cleanup(rdict[extqname],RG)
            outputbam.write(rdict[extqname])  # write read from donor .bam
            used[extqname] = True
        else:
            read = cleanup(read,RG)
            outputbam.write(read) # write read from target .bam

    nadded = 0
    # dump the unused reads from the donor if requested with --all
    if args.all:
        for extqname in rdict.keys():
            if extqname not in used:
                rdict[extqname] = cleanup(rdict[extqname],RG)
                outputbam.write(rdict[extqname])
                nadded += 1
        sys.stderr.write("added " + str(nadded) + " reads due to --all\n")

    targetbam.close()
    donorbam.close()
    outputbam.close()

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='replaces aligned reads in bamfile1 with aligned reads from bamfile2')
    parser.add_argument('-b', '--bam', dest='targetbam', required=True,
                        help='original .bam')
    parser.add_argument('-r', '--replacebam', dest='donorbam', required=True,
                        help='.bam with reads to replace original bam')
    parser.add_argument('-o', '--outputbam', dest='outputbam', required=True)
    parser.add_argument('-i', '--ignoresize', dest='ignoresize', default=0, help="don't replace reads with reads having insert size > ignoresize")
    parser.add_argument('-n', '--namechange', dest='namechange', default=None, help="change all read names by prepending string (passed as -n [string])")
    parser.add_argument('--all', action='store_true', default=False, help="append reads that don't match target .bam")
    parser.add_argument('--keepqual', action='store_true', default=False, help="keep original quality scores, replace read and mapping only")
    args = parser.parse_args()
    main(args)

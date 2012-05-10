#!/bin/env python

import sys,pysam, argparse

def main(args):
    targetbam = pysam.Samfile(args.targetbam, 'rb')
    donorbam  = pysam.Samfile(args.donorbam, 'rb')
    outputbam = pysam.Samfile(args.outputbam, 'wb', template=targetbam)

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

        extqname = ','.join((read.qname,str(read.pos),pairname))
        rdict[extqname] = read
        nr += 1

    sys.stderr.write("loaded " + str(nr) + " reads, reading " + args.targetbam + "...\n")
    for read in targetbam.fetch(until_eof=True):
        pairname = 'F' # read is first in pair
        if read.is_read2:
            pairname = 'S' # read is second in pair
        if not read.is_paired:
            pairname = 'U' # read is unpaired

        extqname = ','.join((read.qname,str(read.pos),pairname))
        if extqname in rdict: # experiment with try/except
            outputbam.write(rdict[extqname])
        else:
            outputbam.write(read)

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
    args = parser.parse_args()
    main(args)

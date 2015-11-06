#!/usr/bin/env python

import pysam
import random
import argparse
import sys
from re import sub
from collections import defaultdict as dd

def split(inbamfn, keep_secondary=False):
    assert inbamfn.endswith('.bam')
    outbam1fn = sub('.bam$', '.pick1.bam', inbamfn)
    outbam2fn = sub('.bam$', '.pick2.bam', inbamfn)

    inbam = pysam.Samfile(inbamfn, 'rb')
    outbam1 = pysam.Samfile(outbam1fn, 'wb', template=inbam)
    outbam2 = pysam.Samfile(outbam2fn, 'wb', template=inbam)

    lastread = None

    reads = dd(list) 

    for read in inbam.fetch(until_eof=True):
        if not read.is_secondary or args.secondary:
            if lastread is not None and read.qname == lastread.qname:
                reads[read.qname].append(read)

            else:
                if lastread is not None:
                    ob = None
                    if random.random() > 0.5:
                        ob = outbam1
                    else:
                        ob = outbam2

                    for outread in reads[lastread.qname]:
                        ob.write(outread)

                    del reads[lastread.qname]
                reads[read.qname].append(read)
            lastread = read
        else:
            print "skipped secondary alignment: ", read.qname

    ob = None
    if random.random() > 0.5:
        ob = outbam1
    else:
        ob = outbam2

    for read in reads[lastread.qname]:
        ob.write(read)

    outbam1.close()
    outbam2.close()
    inbam.close()

def main(args):
    split(args.bam, keep_secondary=args.secondary)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='split one bam into two')
    parser.add_argument('-b', '--bam', required=True, help='input BAM: must be sorted by readname (e.g. samtools sort -n)')
    parser.add_argument('--secondary', default=False, action='store_true', help='keep secondary alignments')
    args = parser.parse_args()
    main(args)

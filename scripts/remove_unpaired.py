#!/usr/bin/env python

import pysam
import sys
from re import sub

if len(sys.argv) == 2:
    assert sys.argv[1].endswith('.bam')
    outname = sub('bam$', 'fixed.bam', sys.argv[1])
    print("output bam: " + outname)

    inbam  = pysam.AlignmentFile(sys.argv[1])
    outbam = pysam.AlignmentFile(outname, 'wb', template=inbam)

    reads = {}

    for read in inbam.fetch(until_eof=True):
        if not read.is_secondary and not read.is_supplementary:
            if read.qname in reads:
                reads[read.qname].append(read)
            else:
                reads[read.qname] = []
                reads[read.qname].append(read)

            if len(reads[read.qname]) == 2:
                for outread in reads[read.qname]:
                    outbam.write(outread)
                del reads[read.qname]

    inbam.close()
    outbam.close()

else:
    sys.exit("usage: %s <BAM file>" % sys.argv[0])

#!/usr/bin/env python

import sys

if len(sys.argv) == 2:
    with open(sys.argv[1], 'r') as sites:
        lastend = 0
        lastchrom = None
        for line in sites:
            chrom, start, end, vaf = line.strip().split()
            start = int(start)
            end   = int(end)
            if chrom != lastchrom:
                lastend = 0

            assert end >= lastend

            if end-lastend > 1000:
                print line.strip()

            lastend = end
            lastchrom = chrom
else:
    print "usage:", sys.argv[0], "<mutation list (chrom, start, end, vaf) MUST BE SORTED>"

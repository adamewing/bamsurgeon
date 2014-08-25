#!/usr/bin/env python

import sys

if len(sys.argv) == 2:
    with open(sys.argv[1], 'r') as sites:
        lastend = 0
        lastchrom = None
        for line in sites:
            chrom, start, end = line.strip().split()[:3]
            start = int(start)
            end   = int(end)
            if chrom != lastchrom:
                lastend = 0

            if end-lastend > 2000:
                print line.strip()

            lastend = end
            lastchrom = chrom
else:
    print "usage:", sys.argv[0], "<mutation list (chrom, start, end) MUST BE SORTED>"

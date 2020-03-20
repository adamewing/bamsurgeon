#!/usr/bin/env python

import sys
import os

from collections import defaultdict as dd


class Site:
    def __init__(self, infn, line):
        self.line  = line.strip().split()
        self.chrom = self.line[0]
        self.start = int(self.line[1])
        self.end   = int(self.line[2])
        self.info  = self.line[3:]
        self.infn  = infn

    def overlap(self, other, sep):
        start = self.start - sep
        end   = self.end + sep

        if min(end, other.end) - max(start, other.start) >= 0:
            return True
        return False

    def __str__(self):
        return '\t'.join(self.line) 

    def __lt__(self, other):
        if self.chrom != other.chrom:
            return self.chrom < other.chrom
        else:
            return self.start < other.start


if len(sys.argv) > 1:
    sep = 2000
    sites = []
    for fn in sys.argv[1:]:
        assert os.path.exists(fn), "file not found: " + fn

        with open(fn, 'r') as insites:
            for line in insites:
                sites.append(Site(fn, line))
    sites.sort()

    sites_by_fn = dd(list)
    lastsite = sites[0]

    for site in sites[1:]:
        if not site.overlap(lastsite, sep):
            sites_by_fn[lastsite.infn].append(lastsite)
            lastsite = site

    if not site.overlap(lastsite, sep):
        sites_by_fn[lastsite.infn].append(lastsite)

    for fn in sites_by_fn:
        outfn = fn + '.sep.txt'
        with open(outfn, 'w') as out:
            for site in sites_by_fn[fn]:
                assert site.infn == fn
                out.write(str(site) + '\n')

else:
    sys.exit("usage: %s <mutation lists (chrom start end ... ) can input multiple files (file1.txt file2.txt ... )>" % sys.argv[0])

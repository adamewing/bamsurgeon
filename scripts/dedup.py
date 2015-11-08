#!/bin/env python

import sys

pad = 10000

class Interval:
    def __init__(self,chr,start,end,etc):
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
        self.info = etc
    def cmp(self,ivl):
        if self.chr == ivl.chr:
            if (self.start <= ivl.start and self.end >= ivl.start) or (self.start <= ivl.end and self.end >= ivl.end):
                if self.eq(ivl):
                    return False
                return True
            else:
                return False
        else:
            return False

    def eq(self,ivl):
        if self.chr == ivl.chr and self.start == ivl.start and self.end == ivl.end and self.info == ivl.info:
            return True
        return False

    def __str__(self):
        out = "\t".join((self.chr,str(self.start+pad),str(self.end-pad))) + "\t" + " ".join(self.info)
        return out

if len(sys.argv) == 2:
    f = open(sys.argv[1])
    ivl_list = []
    for line in f:
        c = line.strip().split()
        chrom = c[0]
        start = int(c[1]) - pad
        end = int(c[2]) + pad
        ivl = Interval(chrom,start,end,c[3:len(c)])
        ivl_list.append(ivl)

    dup = {}
    for i1 in ivl_list:
        for i2 in ivl_list:
            if i1.cmp(i2):
                dup[i1] = True

    for ivl in ivl_list:
        if ivl not in dup:
            print ivl

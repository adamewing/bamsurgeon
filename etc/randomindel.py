#!/bin/env python

from random import *
import sys

vaf = '0.5'

def randomindel(loc):
    chrom, start, end = loc.strip().split()
    i = randint(0,1)
    l = int(expovariate(10)*100)+1 # indel size, exponential dist to bias towards shorter indels
    if i == 0: # DEL
        end = str(int(start) + l)
        return '\t'.join((chrom,start,end,vaf,'DEL'))
    if i == 1: # INS
        iseq = ''
        for s in range(l):
            bp = ['A','T','G','C']
            iseq = iseq + bp[randint(0,3)]
        return '\t'.join((chrom,start,end,vaf,'INS',iseq))

if len(sys.argv) > 1:
    f = open(sys.argv[1],'r')
    for line in f:
        print randomindel(line)
else:
    print "usage:",sys.argv[0],"<BED-3>"

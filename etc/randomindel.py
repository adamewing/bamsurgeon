#!/bin/env python

from random import *
import sys

def randomindel():
    i = randint(0,1)
    l = int(expovariate(10)*100)+1 # indel size, exponential dist to bias towards shorter indels
    if i == 0: # DEL
        return "DEL"
    if i == 1: # INS
        iseq = ''
        for s in range(l):
            bp = ['A','T','G','C']
            iseq = iseq + bp[randint(0,3)]
        return "INS " + iseq 

if len(sys.argv) > 1:
    f = open(sys.argv[1],'r')
    for line in f:
        rndsv = randomindel()
        print line.strip() + "\t" + rndsv
else:
    print "usage:",sys.argv[0],"<BED-3>"

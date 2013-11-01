#!/bin/env python

from random import *
from os import getcwd
import sys

def randomsv():
    pre = os.path.dirname(os.path.realpath(__file__)) + '/ins_seqs/'
    i = randint(0,3)
    if i == 0: # DEL
        dfrac = uniform(0.5,1)
        return "DEL " + str(dfrac)
    if i == 1: # INS
        files = ['alu.fa','line.fa','line_trunc.fa','pseudo.fa']
        insfile = pre + files[randint(0,len(files)-1)]
        tsdlen = randint(0,30)
        return ' '.join(('INS',insfile,str(tsdlen)))
    if i == 2: # INV
        return "INV"
    if i == 3: # DUP
        ndups = randint(1,4)
        return "DUP " + str(ndups)

if len(sys.argv) > 1:
    f = open(sys.argv[1],'r')
    for line in f:
        rndsv = randomsv()
        print line.strip() + "\t" + rndsv

else:
    print "usage:",sys.argv[0],"<BED-3>"

#!/usr/bin/env python

import sys,os

def header():
    print '##fileformat=VCFv4.1'
    print '##INDIVIDUAL=TRUTH'
    print '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    print '##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation in primary">'
    print '##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Frequency">'
    print '##INFO=<ID=DPR,Number=1,Type=Float,Description="Avg Depth in Region (+/- 1bp)">'
    print '##SAMPLE=<ID=TRUTH,Individual="TRUTH",Description="Validation data"'
    print '\t'.join(('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','TRUTH'))

if len(sys.argv) == 2:
    header()
    logdir_files = os.listdir(sys.argv[1])

    for filename in logdir_files:
        if filename.endswith('.log'):
            with open(sys.argv[1] + '/' + filename, 'r') as infile:
                for line in infile:
                    if line.startswith('snv'):
                        #chrom, pos, mut = line.strip().split()
                        c = line.strip().split()
                        chrom = c[1].split(':')[0]
                        pos = c[3]
                        mut = c[4]
                        dpr = c[6]
                        vaf = c[7]
                        ref,alt = mut.split('-->')
                        print "\t".join((chrom,pos,'.',ref,alt,'100','PASS','SOMATIC;VAF=' + vaf + ';DPR=' + dpr,'GT','0/1'))
else:
    print "usage:", sys.argv[0], "<log directory>"

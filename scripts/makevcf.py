#!/usr/bin/env python

import sys,os
import textwrap

def print_header():
    print textwrap.dedent("""\
    ##fileformat=VCFv4.1
    ##phasing=none
    ##INDIVIDUAL=TRUTH
    ##SAMPLE=<ID=TRUTH,Individual="TRUTH",Description="bamsurgeon spike-in">
    ##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
    ##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
    ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
    ##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
    ##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation in primary">
    ##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Frequency">
    ##INFO=<ID=DPR,Number=1,Type=Float,Description="Avg Depth in Region (+/- 1bp)">
    ##INFO=<ID=MATEID,Number=1,Type=String,Description="Breakend mate">
    ##ALT=<ID=INV,Description="Inversion">
    ##ALT=<ID=DUP,Description="Duplication">
    ##ALT=<ID=DEL,Description="Deletion">
    ##ALT=<ID=INS,Description="Insertion">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSPIKEIN""")

if len(sys.argv) == 2:
    print_header()
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

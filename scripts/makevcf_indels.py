#!/usr/bin/env python

import sys
import pysam
import textwrap
import os

def print_header():
    print(textwrap.dedent("""\
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
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSPIKEIN"""))

if len(sys.argv) == 3:

    fa = pysam.Fastafile(sys.argv[2])
    print_header()

    logdir_files = os.listdir(sys.argv[1])
    for filename in logdir_files:
        if filename.endswith('.log'):
            with open(sys.argv[1] + '/' + filename, 'r') as logfile:
        
                for line in logfile:
                    ref = ''
                    alt = ''
                    start = ''
                    info = []
                    info.append('SOMATIC')
        
                    if line.startswith('indel'):
                        indelinfo = line.strip().split()[1].split(':')
                        info.append('VAF=%f' % float(line.strip().split()[7]))
                        if indelinfo[0] == 'DEL':
                            chrom, start, end = indelinfo[1:4]
                            ref = fa.fetch(chrom, int(start)-1, int(end)).upper()
                            alt = ref[0].upper()
        
                        if indelinfo[0] == 'INS':
                            chrom, start, seq = indelinfo[1:4]
                            ref = fa.fetch(chrom, int(start)-1, int(start)).upper()
                            alt = ref.upper() + seq.upper()
        
                        assert ref != '' and alt != '' and start != ''
        
                        print('\t'.join((chrom, start, '.', ref, alt, '100', 'PASS', ';'.join(info), 'GT', '0/1')))
        
else:
    sys.exit("usage: %s <indel log directory> <samtools faidx indexed reference>" % sys.argv[0])

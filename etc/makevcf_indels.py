#!/usr/bin/env python

import sys
import pysam
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

if len(sys.argv) == 3:

    fa = pysam.Fastafile(sys.argv[2])
    print_header()

    with open(sys.argv[1], 'r') as logfile:
        for line in logfile:
            ref = ''
            alt = ''
            start = ''
            info = []
            info.append('SOMATIC')

            if line.startswith('indel'):
                indelinfo = line.strip().split()[1].split(':')
                if indelinfo[0] == 'DEL':
                    chrom, start, end = indelinfo[1:4]
                    ref = fa.fetch(chrom, int(start)-1, int(end))
                    alt = ref[0]

                if indelinfo[0] == 'INS':
                    chrom, start, seq = indelinfo[1:4]
                    ref = fa.fetch(chrom, int(start)-1, int(start))
                    alt = ref + seq

                assert ref != '' and alt != '' and start != ''

                print '\t'.join((chrom, start, '.', ref, alt, '100', 'PASS', ';'.join(info), 'GT', '0/1'))

else:
    print "usage:",sys.argv[0],"<bamsurgeon indel .log file> <samtools faidx indexed reference>"

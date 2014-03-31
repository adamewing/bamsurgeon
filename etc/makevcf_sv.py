#!/usr/bin/env python

import sys
import os
import subprocess
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
    ##ALT=<ID=IGN,Description="Ignore SNVs in Interval">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSPIKEIN""")

def printvcf(chrom, bnd1, bnd2, precise, type, svlen, ref):
    base1 = ref.fetch(chrom, bnd1, bnd1+1) 
    base2 = ref.fetch(chrom, bnd2, bnd2+1)

    alt = '<' + type.upper() + '>'

    info = []
    if not precise:
        info.append('IMPRECISE')
        info.append('CIPOS=-100,100')
        info.append('CIEND=-100,100')

    info.append('SOMATIC')
    info.append('SVTYPE=' + type.upper())
    info.append('END=' + str(bnd2))
    info.append('SVLEN=' + str(svlen))
    
    infostr = ';'.join(info)

    print '\t'.join((chrom, str(bnd1), '.', base1, alt, '100', 'PASS', infostr, 'GT', './.'))

    '''
    Finish this later: we'd like to output breakends rather than events but it'll require special consideration of each SV type

    if type == 'del':
        id1 = '_'.join((type,str(n),'A'))
        id2 = '_'.join((type,str(n),'B'))
        alt1 = base2 + '[' + chrom + ':' + str(bnd2) + '['
        alt2 = base1 + '[' + chrom + ':' + str(bnd1) + '['
    '''

def precise_interval(mutline, ref):
    m = mutline.split()
    chrom, refstart, refend = m[1:4]
    refstart = int(refstart)
    refend   = int(refend)

    if m[0] == 'ins':
        contigstart = int(m[6])
        contigend   = int(m[6])+1
    else:
        contigstart = int(m[6])
        contigend   = int(m[7])

    precise = True 
    bnd1 = refstart + contigstart
    bnd2 = refstart + contigend

    assert bnd1 < bnd2

    printvcf(chrom, bnd1, bnd2, precise, m[0], bnd2-bnd1, ref)


def ignore_interval(mutline, ref):
    m = mutline.split()
    chrom, refstart, refend = m[1:4]
    refstart = int(refstart)
    refend   = int(refend)

    assert refstart < refend

    printvcf(chrom, refstart, refend, True, 'IGN', refend-refstart, ref)


if len(sys.argv) == 3:
    print_header()

    ref = pysam.Fastafile(sys.argv[2])

    with open(sys.argv[1], 'r') as log:
        for line in log:
            for mutype in ('ins', 'del', 'inv', 'dup'):
                if line.startswith(mutype):
                    precise_interval(line.strip(), ref)
                    ignore_interval(line.strip(), ref)

else:
    print "usage:",sys.argv[0],"<bamsurgeon SV .log file> <samtools indexed fasta reference>"

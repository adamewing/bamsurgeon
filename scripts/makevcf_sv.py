#!/usr/bin/env python

import sys
import os
import subprocess
import pysam
import textwrap
import argparse

from uuid import uuid4

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
    ##ALT=<ID=IGN,Description="Ignore SNVs in Interval">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSPIKEIN"""))

def printvcf(chrom, bnd1, bnd2, precise, type, svlen, ref, id, svfrac):
    base1 = ref.fetch(chrom, bnd1-1, bnd1) 
    base2 = ref.fetch(chrom, bnd2-1, bnd2)

    alt = '<' + type.upper() + '>'

    info = []
    if not precise:
        info.append('IMPRECISE')
        info.append('CIPOS=-100,100')
        info.append('CIEND=-100,100')

    info.append('SOMATIC')
    info.append('SVTYPE=' + type.upper())
    info.append('END=' + str(bnd2))
    info.append('SVLEN=' + str(abs(svlen)))
    info.append('VAF=' + svfrac)
    
    infostr = ';'.join(info)

    print('\t'.join((chrom, str(bnd1), id, base1, alt, '100', 'PASS', infostr, 'GT', './.')))


def precise_interval(mutline, ref):
    precise = True 
    m = mutline.split()
    chrom, refstart, refend = m[1:4]
    svfrac = m[-1]
    refstart = int(refstart)
    refend   = int(refend)
    id = '.'

    if m[0].startswith('big'):
        bnd1 = int(m[2])+int(m[5])
        bnd2 = int(m[7])+int(m[9])
        m[0] = m[0].replace('big', '')

    elif m[0] == 'ins':
        bnd1 = refstart+int(m[6])
        bnd2 = refstart+int(m[6])+1
        id = m[7]

    elif m[0] == 'trn':
        chr1 = chrom
        chr2 = m[6]

        bnd1 = int(m[2])+int(m[5])
        bnd2 = int(m[7])+int(m[9])

        id1 = str(uuid4()).split('-')[0]
        id2 = str(uuid4()).split('-')[0]

        base1 = ref.fetch(chr1, bnd1-1, bnd1) 
        base2 = ref.fetch(chr2, bnd2-1, bnd2)

        alt1 = '%s[%s:%d[' % (base1, chr2, bnd2) 
        alt2 = ']%s:%d]%s' % (chr1, bnd1, base2) 

        print('\t'.join((chr1, str(bnd1), id1, base1, alt1, '100', 'PASS', 'SOMATIC;SVTYPE=BND;PRECISE;MATEID='+id2+';VAF='+svfrac, 'GT', './.')))
        print('\t'.join((chr2, str(bnd2), id2, base2, alt2, '100', 'PASS', 'SOMATIC;SVTYPE=BND;PRECISE;MATEID='+id1+';VAF='+svfrac, 'GT', './.')))

    else:
        contigstart = int(m[6])
        contigend   = int(m[7])

        bnd1 = refstart + contigstart
        bnd2 = refstart + contigend

    if m[0] != 'trn': printvcf(chrom, bnd1, bnd2, precise, m[0], bnd2-bnd1, ref, id, svfrac)


def ignore_interval(mutline, ref):
    m = mutline.split()
    chrom, refstart, refend = m[1:4]
    refstart = int(refstart)
    refend   = int(refend)

    assert refstart < refend

    printvcf(chrom, refstart, refend, True, 'IGN', refend-refstart, ref, '.', '0.0')
    if m[0] == 'trn': printvcf(m[6], int(m[7]), int(m[8]), True, 'IGN', int(m[8])-int(m[7]), ref)


def main(args):
    print_header()

    ref = pysam.Fastafile(args.ref)

    logdir_files = os.listdir(args.logdir)

    for filename in logdir_files:
        if filename.endswith('.log'):
            with open(args.logdir + '/' + filename, 'r') as log:
                for line in log:
                    for mutype in ('ins', 'del', 'inv', 'dup', 'trn', 'bigdup', 'biginv', 'bigdel'):
                        if line.startswith(mutype):
                            precise_interval(line.strip(), ref)
                            if args.mask:
                                ignore_interval(line.strip(), ref)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="make VCF 'truth' file given log file (hint: concatenate them) from addsv.py")
    parser.add_argument('-r', '--ref', dest='ref', required=True, help="reference indexed with samtools faidx")
    parser.add_argument('-l', '--logdir', dest='logdir', required=True, help="log directory from addsv.py")
    parser.add_argument('--mask', action="store_true", default=False, help="output contig intervals, used to mask accidental SNVs if combining mutation types in one BAM")
    args = parser.parse_args()
    main(args)

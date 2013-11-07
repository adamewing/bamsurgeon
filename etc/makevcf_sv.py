#!/usr/bin/env python

import sys
import os
import random
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
    ##INFO=<ID=MATEID,Number=1,Type=String,Description="Breakend mate">
    ##ALT=<ID=INV,Description="Inversion">
    ##ALT=<ID=DUP,Description="Duplication">
    ##ALT=<ID=DEL,Description="Deletion">
    ##ALT=<ID=INS,Description="Insertion">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSPIKEIN""")

def printvcflines(chrom, bnd1, bnd2, precise, ref, n, type, svlen):
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

def process_block(lines, ref):
    before_str, before_seq, mut_str, after_str, after_seq = lines
    assert before_seq.startswith(after_seq[:100])
    mut = mut_str.split()

    m = mut_str.split()
    chrom, refstart, refend = m[1:4]
    refstart = int(refstart)
    refend   = int(refend)
    refseq   = ref.fetch(chrom, refstart, refend)

    rnd = str(random.random())
    tgtfa = 'tmp.' + rnd + '.tgt.fa'
    qryfa = 'tmp.' + rnd + '.qry.fa'

    tgt = open(tgtfa, 'w')
    qry = open(qryfa, 'w')

    tgt.write('>ref' + '\n' + refseq + '\n')
    qry.write('>qry' + '\n' + before_seq[:100] + '\n')

    tgt.close()
    qry.close()

    cmd = ['exonerate', '--bestn', '1', '--showalignment','0', '--ryo', 'SUMMARY\t%s\t%qab\t%qae\t%tab\t%tae\n', qryfa, tgtfa]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    best = []
    topscore = 0

    for pline in p.stdout.readlines():
        if pline.startswith('SUMMARY'):
            c = pline.strip().split()
            if int(c[1]) > topscore:
                topscore = int(c[1])
                best = c

    bnd1 = None
    bnd2 = None
    precise = False

    if best:
        asmstart = refstart + int(best[4])
        bnd1 = asmstart + int(mut[6])

        if mut[0] != 'ins':
            bnd2 = asmstart + int(mut[7])

        else:
            bnd2 = bnd1
            bnd1 -= int(mut[8])

        if topscore == 500:
            precise = True

    p.stdout.close()

    os.remove(tgtfa)
    os.remove(qryfa)

    if None in (bnd1, bnd2):
        precise = False
        bnd1 = refstart
        bnd2 = refend

    return chrom, bnd1, bnd2, precise, mut[0], len(before_seq)-len(after_seq)

if len(sys.argv) == 3:
    assert sys.argv[1].endswith('.log')
    assert os.path.exists(sys.argv[2] + '.fai')

    ref = pysam.Fastafile(sys.argv[2])

    ln = 0
    lines = []
    print_header()

    with open(sys.argv[1], 'r') as log:
        for line in log:
            ln += 1
            lines.append(line.strip())
            if ln % 5 == 0:
                chrom, bnd1, bnd2, precise, type, svlen = process_block(lines, ref)
                if None not in (bnd1, bnd2):
                    printvcflines(chrom, bnd1, bnd2, precise, ref, ln, type, svlen)
                lines = []
else:
    print "usage:",sys.argv[0],"<bamsurgeon SV .log file> <samtools faidx indexed ref fasta>"

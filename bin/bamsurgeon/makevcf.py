#!/usr/bin/env python

import os
import pysam
import datetime
import hashlib


def vcf_header(ref_fa):
    contigs = ''
    fasta_file = pysam.Fastafile(ref_fa)
    for contig in fasta_file.references:
        contigs += '##contig=<ID=%s,length=%d>' % (contig, fasta_file.get_reference_length(contig)) + '\n'

    header = '##fileformat=VCFv4.1' + '\n'
    header += '##fileDate=%s' % (str(datetime.date.today())) + '\n'
    header += '##phasing=none' + '\n'
    header += contigs
    header += '##INDIVIDUAL=TRUTH' + '\n'
    header += '##SAMPLE=<ID=TRUTH,Individual="TRUTH",Description="BAMSurgeon spike-in">' + '\n'
    header += '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">' + '\n'
    header += '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">' + '\n'
    header += '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">' + '\n'
    header += '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' + '\n'
    header += '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">' + '\n'
    header += '##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation in primary">' + '\n'
    header += '##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">' + '\n'
    header += '##INFO=<ID=DPR,Number=1,Type=Float,Description="Density of reads supporting the variant">' + '\n'
    header += '##INFO=<ID=MATEID,Number=1,Type=String,Description="Breakend mate">' + '\n'
    header += '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">' + '\n'
    header += '##ALT=<ID=INV,Description="Inversion">' + '\n'
    header += '##ALT=<ID=DEL,Description="Deletion">' + '\n'
    header += '##ALT=<ID=INS,Description="Insertion">' + '\n'
    header += '##ALT=<ID=DUP,Description="Duplication">' + '\n'
    header += '##ALT=<ID=IGN,Description="Ignore SNVs in Interval">' + '\n'
    header += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + '\n'
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSPIKEIN'
    return header


def write_vcf_snv(logdir, ref_fa, vcf_fn):
    header = vcf_header(ref_fa)
    logdir_files = os.listdir(logdir)

    with open(vcf_fn, 'w') as vcf_out:
        vcf_out.write(header + '\n')
        for filename in logdir_files:
            if filename.endswith('.log'):
                with open(logdir + '/' + filename, 'r') as infile:
                    for line in infile:
                        if line.startswith('snv'):
                            c = line.strip().split()
                            chrom = c[1].split(':')[0]
                            pos = c[3]
                            mut = c[4]
                            dpr = c[6]
                            vaf = c[7]
                            ref, alt = mut.split('-->')
                            vcf_out.write('\t'.join((chrom, pos, '.', ref.upper(), alt.upper(), '100', 'PASS',
                                          'SOMATIC;VAF=' + vaf + ';DPR=' + dpr, 'GT', '0/1')) + '\n')


def write_vcf_indel(logdir, ref_fa, vcf_fn):
    fa = pysam.Fastafile(ref_fa)
    header = vcf_header(ref_fa)

    logdir_files = os.listdir(logdir)
    with open(vcf_fn, 'w') as vcf_out:
        vcf_out.write(header + '\n')

        for filename in logdir_files:
            if filename.endswith('.log'):
                with open(logdir + '/' + filename, 'r') as logfile:

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

                            vcf_out.write('\t'.join((chrom, start, '.', ref, alt, '100',
                                          'PASS', ';'.join(info), 'GT', '0/1')) + '\n')


def sv_vcf_line(chrom, bnd1, bnd2, precise, type, svlen, ref, id, svfrac, fh):
    base1 = ref.fetch(chrom, bnd1-1, bnd1).upper()

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

    fh.write('\t'.join((chrom, str(bnd1), id, base1, alt, '100', 'PASS', infostr, 'GT', './.')) + '\n')


def sv_vcf_precise_interval(mutline, ref, fh):
    precise = True
    m = mutline.split()
    chrom, refstart, refend = m[1:4]
    svfrac = m[-1]
    refstart = int(refstart)
    refend = int(refend)
    id_ = '.'

    if m[0].startswith('big'):
        bnd1 = int(m[2])+int(m[5])
        bnd2 = int(m[7])+int(m[9])
        bnd1, bnd2 = sorted([bnd1, bnd2])
        m[0] = m[0].replace('big', '')
    elif m[0] == 'ins':
        bnd1 = refstart+int(m[6])
        length = int(m[9])
        bnd2 = bnd1 + length
        id_ = m[7] if m[7] != 'None' else '.'
    elif m[0] == 'trn':
        chr1 = chrom
        chr2 = m[6]

        bnd1 = int(m[2])+int(m[5])
        bnd2 = int(m[7])+int(m[9])

        id1 = hashlib.md5(mutline.encode()).hexdigest()
        id2 = hashlib.md5((mutline+mutline).encode()).hexdigest()

        base1 = ref.fetch(chr1, bnd1-1, bnd1).upper()
        base2 = ref.fetch(chr2, bnd2-1, bnd2).upper()

        flip_left = True if m[10] == 'True' else False
        flip_right = True if m[11] == 'True' else False
        bracket1 = '[' if not flip_right else ']'
        prefix1, suffix1 = (base1, '') if not flip_left else ('', base1)
        bracket2 = '[' if flip_left else ']'
        prefix2, suffix2 = (base2, '') if flip_right else ('', base2)

        alt1 = prefix1 + bracket1 + chr2 + ':' + str(bnd2) + bracket1 + suffix1
        alt2 = prefix2 + bracket2 + chr1 + ':' + str(bnd1) + bracket2 + suffix2

        fh.write('\t'.join((chr1, str(bnd1), id1, base1, alt1, '100', 'PASS',
                 'SOMATIC;SVTYPE=BND;PRECISE;MATEID='+id2+';VAF='+svfrac, 'GT', './.')) + '\n')
        fh.write('\t'.join((chr2, str(bnd2), id2, base2, alt2, '100', 'PASS',
                 'SOMATIC;SVTYPE=BND;PRECISE;MATEID='+id1+';VAF='+svfrac, 'GT', './.')) + '\n')
    else:
        contigstart = int(m[6])
        contigend = int(m[7])
        bnd1 = refstart + contigstart
        bnd2 = refstart + contigend

    if m[0] != 'trn':
        sv_vcf_line(chrom, bnd1, bnd2, precise, m[0], bnd2-bnd1, ref, id_, svfrac, fh)


def write_vcf_sv(logdir, ref_fa, vcf_fn):
    header = vcf_header(ref_fa)

    ref = pysam.Fastafile(ref_fa)

    logdir_files = os.listdir(logdir)

    with open(vcf_fn, 'w') as vcf_out:
        vcf_out.write(header + '\n')

        for filename in logdir_files:
            if filename.endswith('.log'):
                with open(logdir + '/' + filename, 'r') as log:
                    for line in log:
                        for mutype in ('ins', 'del', 'inv', 'dup', 'trn', 'bigdup', 'biginv', 'bigdel'):
                            if line.startswith(mutype):
                                sv_vcf_precise_interval(line.strip(), ref, vcf_out)

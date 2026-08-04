#!/usr/bin/env python

import os
import pysam
import datetime
import hashlib

from bamsurgeon.records import INTERVAL_KINDS


# small variants are written as heterozygous spike-ins, SVs as uncalled --
# matching what the log-scraping writers they replace emitted
HET = (0, 1)
UNCALLED = (None, None)


SV_HEADER_LINES = (
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
    '##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Somatic mutation in primary">',
    '##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">',
    '##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">',
    '##INFO=<ID=MATEID,Number=1,Type=String,Description="Breakend mate">',
    '##INFO=<ID=NDUPS,Number=1,Type=Integer,Description="Number of additional tandem copies">',
    '##INFO=<ID=TSDLEN,Number=1,Type=Integer,Description="Target site duplication length">',
    '##INFO=<ID=DPR,Number=1,Type=Float,Description="Read depth in the spiked region after realignment">',
    '##ALT=<ID=DEL,Description="Deletion">',
    '##ALT=<ID=DUP,Description="Duplication">',
    '##ALT=<ID=INV,Description="Inversion">',
    '##ALT=<ID=INS,Description="Insertion">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
)


def sv_variant_header(ref_fa, salt=None):
    ''' VariantHeader built from the reference .fai '''
    header = pysam.VariantHeader()

    fa = pysam.Fastafile(ref_fa)
    for contig in fa.references:
        header.add_line('##contig=<ID=%s,length=%d>'
                        % (contig, fa.get_reference_length(contig)))

    header.add_line('##fileDate=%s' % str(datetime.date.today()))
    header.add_line('##phasing=none')
    if salt is not None:
        header.add_line('##bamsurgeonSalt=%s' % salt)
    header.add_line('##INDIVIDUAL=TRUTH')
    header.add_line('##SAMPLE=<ID=TRUTH,Individual="TRUTH",Description="BAMSurgeon spike-in">')

    for line in SV_HEADER_LINES:
        header.add_line(line)

    header.add_sample('SPIKEIN')

    return header


def _breakend_alts(rec, base1, base2):
    '''
    The four bracket forms. Preserved from the previous writer, which built
    these by hand from a log line.
    '''
    bracket1 = ']' if rec.flip_right else '['
    prefix1, suffix1 = ('', base1) if rec.flip_left else (base1, '')
    bracket2 = '[' if rec.flip_left else ']'
    prefix2, suffix2 = (base2, '') if rec.flip_right else ('', base2)

    alt1 = '%s%s%s:%d%s%s' % (prefix1, bracket1, rec.mate_chrom, rec.mate_pos,
                              bracket1, suffix1)
    alt2 = '%s%s%s:%d%s%s' % (prefix2, bracket2, rec.chrom, rec.pos,
                              bracket2, suffix2)
    return alt1, alt2


def write_vcf(records, ref_fa, vcf_fn, salt=None):
    '''
    Write truth records of any kind, sorted by coordinate.

    The previous writer walked os.listdir(logdir) and re-parsed per-mutation
    log lines, so output order was arbitrary and the result could not be
    tabix-indexed -- meaning bamsurgeon could not produce input that its own
    evaluator.py would accept.
    '''
    header = sv_variant_header(ref_fa, salt=salt)
    ref = pysam.Fastafile(ref_fa)

    order = {c: i for i, c in enumerate(header.contigs)}

    def base_at(chrom, pos):
        return ref.fetch(chrom, pos - 1, pos).upper() or 'N'

    rows = []

    for rec in records:
        if rec.kind == 'BND':
            base1 = base_at(rec.chrom, rec.pos)
            base2 = base_at(rec.mate_chrom, rec.mate_pos)
            alt1, alt2 = _breakend_alts(rec, base1, base2)

            id1 = hashlib.md5(rec.key().encode()).hexdigest()
            id2 = hashlib.md5((rec.key() + '|mate').encode()).hexdigest()

            rows.append((rec.chrom, rec.pos, id1, base1, alt1,
                         {'SVTYPE': 'BND', 'PRECISE': True, 'SOMATIC': True,
                          'MATEID': id2, 'VAF': rec.vaf}, None, UNCALLED))
            rows.append((rec.mate_chrom, rec.mate_pos, id2, base2, alt2,
                         {'SVTYPE': 'BND', 'PRECISE': True, 'SOMATIC': True,
                          'MATEID': id1, 'VAF': rec.vaf}, None, UNCALLED))

        elif rec.is_small:
            info = {'SOMATIC': True, 'VAF': rec.vaf}
            if rec.dpr is not None:
                info['DPR'] = rec.dpr

            rows.append((rec.chrom, rec.pos, '.', rec.ref_allele.upper(),
                         rec.alt_allele.upper(), info, None, HET))

        elif rec.kind in INTERVAL_KINDS:
            info = {'SVTYPE': rec.kind, 'PRECISE': True, 'SOMATIC': True,
                    'SVLEN': int(rec.svlen), 'VAF': rec.vaf}
            if rec.kind == 'DUP' and rec.ndups:
                info['NDUPS'] = int(rec.ndups)
            if rec.kind == 'INS' and rec.tsdlen:
                info['TSDLEN'] = int(rec.tsdlen)

            rows.append((rec.chrom, rec.pos, rec.ins_id or '.',
                         base_at(rec.chrom, rec.pos), '<%s>' % rec.kind,
                         info, rec.end, UNCALLED))

    rows.sort(key=lambda r: (order.get(r[0], len(order)), r[1]))

    with pysam.VariantFile(vcf_fn, 'w', header=header) as out:
        for chrom, pos, rid, refbase, alt, info, end, gt in rows:
            vrec = out.new_record(
                contig=chrom,
                start=pos - 1,
                stop=end if end is not None else pos,
                alleles=(refbase, alt),
                id=None if rid == '.' else rid,
                qual=100,
                filter='PASS',
            )
            for k, v in info.items():
                vrec.info[k] = v
            vrec.samples['SPIKEIN']['GT'] = gt
            out.write(vrec)


def vcf_header(ref_fa, salt=None):
    contigs = ''
    fasta_file = pysam.Fastafile(ref_fa)
    for contig in fasta_file.references:
        contigs += '##contig=<ID=%s,length=%d>' % (contig, fasta_file.get_reference_length(contig)) + '\n'

    header = '##fileformat=VCFv4.1' + '\n'
    header += '##fileDate=%s' % (str(datetime.date.today())) + '\n'
    header += '##phasing=none' + '\n'
    # which reads were selected is a function of this salt; recording it makes
    # a run reproducible even when --seed was not given
    if salt is not None:
        header += '##bamsurgeonSalt=%s' % salt + '\n'
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


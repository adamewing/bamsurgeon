'''
Tests for scripts/evaluator.py.

These need only pysam and the reference FASTA -- no aligner, no picard -- so
unlike test_combined.py they run anywhere the package imports.

The expected numbers are hand-countable from the fixtures each test builds.
'''

import os
import subprocess
import sys

import pysam
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EVALUATOR = os.path.join(REPO, 'scripts', 'evaluator.py')
REF = os.path.join(REPO, 'test_data', 'Homo_sapiens_chr22_assembly19.fasta')

HEADER = '''##fileformat=VCFv4.1
##contig=<ID=22,length=51304566>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="CI around POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="CI around END">
##ALT=<ID=DEL,Description="Deletion">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
'''


def write_vcf(path, rows):
    with open(path, 'w') as fh:
        fh.write(HEADER)
        for row in rows:
            fh.write('\t'.join(str(f) for f in row) + '\n')
    return path


def indexed(path):
    ''' bgzip + tabix, which evaluator requires of the truth VCF '''
    return pysam.tabix_index(path, preset='vcf', force=True)


def run_evaluator(sub, truth, mutype, ref=None, extra=()):
    cmd = [sys.executable, EVALUATOR, '-v', sub, '-t', truth, '-m', mutype]
    if ref:
        cmd += ['-f', ref]
    cmd += list(extra)

    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr[-2000:]

    counts, scores = None, None
    for line in proc.stdout.splitlines():
        if line.startswith('tpcount'):
            counts = tuple(int(x) for x in line.split(':')[1].split(','))
        if line.startswith('precision'):
            scores = tuple(float(x) for x in line.split(':')[1].split(','))

    assert counts is not None and scores is not None, proc.stdout
    return counts, scores


def snv(pos, ref, alt):
    return ('22', pos, '.', ref, alt, 100, 'PASS', '.', 'GT', '0/1')


def test_snv_precision_recall(tmp_path):
    ''' three shared sites, two submission-only, two truth-only '''
    truth = write_vcf(str(tmp_path / 'truth.vcf'), [
        snv(33714965, 'A', 'G'), snv(33769483, 'T', 'A'), snv(33908770, 'C', 'T'),
        snv(34166720, 'T', 'C'), snv(34264774, 'G', 'A'),
    ])
    sub = write_vcf(str(tmp_path / 'sub.vcf'), [
        snv(33714965, 'A', 'G'), snv(33769483, 'T', 'A'), snv(33908770, 'C', 'T'),
        snv(33999999, 'A', 'T'), snv(34100000, 'C', 'G'),
    ])

    counts, scores = run_evaluator(sub, indexed(truth), 'SNV')

    assert counts == (3, 2, 5, 5)                  # tp, fp, sub, tru
    assert scores == pytest.approx((0.6, 0.6, 0.6))


def test_snv_tp_fp_fn_output(tmp_path):
    ''' the three output VCFs name the right sites '''
    truth = write_vcf(str(tmp_path / 'truth.vcf'), [
        snv(33714965, 'A', 'G'), snv(34166720, 'T', 'C'),
    ])
    sub = write_vcf(str(tmp_path / 'sub.vcf'), [
        snv(33714965, 'A', 'G'), snv(33999999, 'A', 'T'),
    ])

    tp = str(tmp_path / 'tp.vcf')
    fp = str(tmp_path / 'fp.vcf')
    fn = str(tmp_path / 'fn.vcf')

    run_evaluator(sub, indexed(truth), 'SNV',
                  extra=['--tp', tp, '--fp', fp, '--fn', fn])

    def positions(path):
        with pysam.VariantFile(path) as v:
            return sorted(r.pos for r in v)

    assert positions(tp) == [33714965]
    assert positions(fp) == [33999999]
    assert positions(fn) == [34166720]


def test_sv_overlap_matches(tmp_path):
    '''
    Two deletions that overlap count as a match.

    SV mode used to be unreachable: expand_sv_ends subscripted a scalar
    INFO/END, and the handler for the resulting TypeError referenced an
    undefined `logger`, so every SV evaluation died with a NameError.
    '''
    truth = write_vcf(str(tmp_path / 'truth.vcf'), [
        ('22', 33871043, '.', 'G', '<DEL>', 100, 'PASS',
         'SVTYPE=DEL;END=33884754', 'GT', './.'),
    ])
    sub = write_vcf(str(tmp_path / 'sub.vcf'), [
        ('22', 33871100, '.', 'G', '<DEL>', 100, 'PASS',
         'SVTYPE=DEL;END=33884700', 'GT', './.'),
    ])

    counts, scores = run_evaluator(sub, indexed(truth), 'SV')

    assert counts == (1, 0, 1, 1)
    assert scores == pytest.approx((1.0, 1.0, 1.0))


def test_sv_non_overlap_does_not_match(tmp_path):
    truth = write_vcf(str(tmp_path / 'truth.vcf'), [
        ('22', 33871043, '.', 'G', '<DEL>', 100, 'PASS',
         'SVTYPE=DEL;END=33884754', 'GT', './.'),
    ])
    sub = write_vcf(str(tmp_path / 'sub.vcf'), [
        ('22', 34100000, '.', 'G', '<DEL>', 100, 'PASS',
         'SVTYPE=DEL;END=34110000', 'GT', './.'),
    ])

    counts, _ = run_evaluator(sub, indexed(truth), 'SV')

    assert counts == (0, 1, 1, 1)


def _homopolymer(chrom='22', start=33700000, end=33720000, minlen=6):
    ''' first run of >= minlen identical bases, as a 1-based position '''
    import re
    fa = pysam.Fastafile(REF)
    seq = fa.fetch(chrom, start, end).upper()
    m = re.search(r'(A{%d,}|C{%d,}|G{%d,}|T{%d,})' % ((minlen,) * 4), seq)
    assert m, 'no homopolymer found in the test region'
    return start + m.start() + 1, m.group()[0]


def test_indel_haplotype_aware_rescue(tmp_path):
    '''
    Two insertions at different positions in a homopolymer produce the same
    haplotype, so they match only when --ref enables the allele-aware
    comparison. This is the behaviour have_identical_haplotypes() exists for.
    '''
    pos, base = _homopolymer()

    truth = write_vcf(str(tmp_path / 'truth.vcf'),
                      [('22', pos, '.', base, base + base, 100, 'PASS', '.', 'GT', '0/1')])
    sub = write_vcf(str(tmp_path / 'sub.vcf'),
                    [('22', pos + 3, '.', base, base + base, 100, 'PASS', '.', 'GT', '0/1')])

    truth_idx = indexed(truth)

    # positions differ, so a naive comparison finds nothing
    counts, _ = run_evaluator(sub, truth_idx, 'INDEL')
    assert counts == (0, 1, 1, 1)

    # with the reference, the identical haplotype is recognised
    counts, scores = run_evaluator(sub, truth_idx, 'INDEL', ref=REF)
    assert counts == (1, 0, 1, 1)
    assert scores == pytest.approx((1.0, 1.0, 1.0))


def test_a_truth_site_is_only_used_once(tmp_path):
    '''
    Two identical submission records cannot both match one truth record.

    Tracking which truth record was consumed used to follow the loop
    variable rather than the record that matched.
    '''
    truth = write_vcf(str(tmp_path / 'truth.vcf'), [snv(33714965, 'A', 'G')])
    sub = write_vcf(str(tmp_path / 'sub.vcf'), [
        snv(33714965, 'A', 'G'), snv(33714965, 'A', 'G'),
    ])

    counts, _ = run_evaluator(sub, indexed(truth), 'SNV')

    tp, fp = counts[0], counts[1]
    assert tp == 1, 'one truth site matched twice'
    assert fp == 1

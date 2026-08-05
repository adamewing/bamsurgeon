'''
Copy-number awareness.

Two things that both follow from duplicated sequence carrying more reads:
--maxdepth is an absolute cap and must rise with local copy number, and the
validator can predict a duplication's depth rather than only checking that it
moved the right way.

The unit tests here need nothing but pysam. The integration test needs the
toolchain from scripts/setup_env.sh.
'''

import os
import shutil
import subprocess
import sys

import pysam
import pytest

from conftest import BIN, REF, data, have, sort_and_index, truth_vcf_path

from bamsurgeon.copynumber import (depth_scale, dup_depth_scale,  # noqa: E402
                                   dup_intervals)
from bamsurgeon.records import MutationRecord  # noqa: E402
from bamsurgeon.makevcf import write_vcf  # noqa: E402
from bamsurgeon.validate import validate_vcf, Thresholds, OK, LOW_VAF  # noqa: E402


BAM = data('testregion_realign.bam')

# the duplication test_sv.txt uses, and a position well inside it
DUP_START, DUP_END = 34271043, 34284754
SNV_POS = 34278000


def test_dup_depth_scale():
    ''' (1 - vaf) + vaf * (ndups + 1) '''
    assert dup_depth_scale(1.0, 1) == pytest.approx(2.0)
    assert dup_depth_scale(1.0, 3) == pytest.approx(4.0)
    assert dup_depth_scale(0.9, 3) == pytest.approx(3.7)
    assert dup_depth_scale(0.5, 1) == pytest.approx(1.5)

    # nothing carries the duplication, so depth is unchanged
    assert dup_depth_scale(0.0, 3) == pytest.approx(1.0)


def test_depth_scale_from_intervals():
    intervals = [('22', 100, 200, 3.7)]

    assert depth_scale('22', 150, 160, intervals=intervals) == 3.7
    assert depth_scale('22', 199, 260, intervals=intervals) == 3.7

    # abutting, not overlapping
    assert depth_scale('22', 200, 260, intervals=intervals) == 1.0
    assert depth_scale('21', 150, 160, intervals=intervals) == 1.0


def test_overlapping_intervals_take_the_maximum():
    ''' not the product: they describe one locus, not stacked amplifications '''
    intervals = [('22', 100, 200, 2.0), ('22', 150, 300, 3.0)]

    assert depth_scale('22', 160, 170, intervals=intervals) == 3.0


def test_depth_scale_from_cnvfile():
    ''' the fixture declares absolute CN 3 across the region, so 1.5x diploid '''
    scale = depth_scale('22', 33900000, 33900100,
                        cnvfile=data('test_cnvlist.txt.gz'))

    assert scale == pytest.approx(1.5)


def test_dup_intervals_converts_to_zero_based_half_open():
    records = [
        MutationRecord('DUP', '22', 101, 200, vaf=1.0, ndups=1),
        MutationRecord('DEL', '22', 300, 400, vaf=1.0),
    ]

    assert dup_intervals(records) == [('22', 100, 200, 2.0)]


def write_dup_and_snv(path):
    ''' a duplication and an SNV inside it, as one VCF for the unified CLI '''
    fa = pysam.Fastafile(REF)
    ref = fa.fetch('22', SNV_POS - 1, SNV_POS).upper()

    write_vcf([
        MutationRecord('DUP', '22', DUP_START, DUP_END, svlen=DUP_END - DUP_START,
                       vaf=0.9, ndups=3),
        MutationRecord('SNV', '22', SNV_POS, SNV_POS, ref_allele=ref,
                       alt_allele='G' if ref != 'G' else 'C', vaf=0.5),
    ], REF, path)

    return path


requires_toolchain = pytest.mark.skipif(
    not (have('bwa', 'samtools', 'wgsim')),
    reason='needs bwa, samtools and wgsim; see scripts/setup_env.sh')


@requires_toolchain
def test_maxdepth_rises_with_local_copy_number(workdir):
    '''
    An SNV inside a duplication is not dropped for being too deep.

    A unified run applies SVs first, so by the time the SNV is placed the
    locus carries 3.7x the input depth -- entirely legitimately. Against the
    raw --maxdepth this aborts with "depth at site is greater than cutoff"
    and the whole run exits with no successful mutations.
    '''
    varfile = write_dup_and_snv(os.path.join(workdir, 'dupsnv.vcf'))
    outbam = os.path.join(workdir, 'dupsnv.bam')

    cmd = [sys.executable, '-O', '-m', 'bamsurgeon.cli.main',
           '-v', varfile, '-f', BAM, '-r', REF, '-o', outbam,
           '-p', '4', '--seed', '1234', '--aligner', 'mem',
           # input depth is ~55x, so this cap is only exceeded because of
           # the duplication
           '--maxdepth', '100',
           '--tmpdir', os.path.join(workdir, 'bs.tmp'),
           '--vcf', workdir + os.sep]

    proc = subprocess.run(cmd, cwd=workdir, env=dict(os.environ, PYTHONPATH=BIN),
                          capture_output=True, text=True)

    assert proc.returncode == 0, proc.stderr[-3000:]
    assert 'local copy number' in proc.stderr, proc.stderr[-2000:]

    sort_and_index(outbam)

    truth = truth_vcf_path(workdir, outbam, 'bamsurgeon', varfile)
    reports = validate_vcf(truth, outbam, REF, orig_bam=BAM,
                           thresholds=Thresholds())

    kinds = sorted(r.kind for r in reports)
    assert kinds == ['DUP', 'SNV'], kinds

    failures = [(r.kind, r.status, r.observation.note)
                for r in reports if r.status != OK]
    assert not failures, failures


@requires_toolchain
def test_a_duplication_must_reach_its_predicted_depth(workdir):
    '''
    Direction alone passes a 4-copy duplication that arrived at 1.02x.

    The spike-in below is honest and validates. Reasserting the same BAM
    against a truth VCF that claims twice the copies must not, or the check
    is only testing that depth went up at all.
    '''
    varfile = write_dup_and_snv(os.path.join(workdir, 'dupsnv.vcf'))
    outbam = os.path.join(workdir, 'dupsnv.bam')

    cmd = [sys.executable, '-O', '-m', 'bamsurgeon.cli.main',
           '-v', varfile, '-f', BAM, '-r', REF, '-o', outbam,
           '-p', '4', '--seed', '1234', '--aligner', 'mem',
           '--tmpdir', os.path.join(workdir, 'bs.tmp'),
           '--vcf', workdir + os.sep]

    proc = subprocess.run(cmd, cwd=workdir, env=dict(os.environ, PYTHONPATH=BIN),
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr[-3000:]

    sort_and_index(outbam)

    truth = truth_vcf_path(workdir, outbam, 'bamsurgeon', varfile)

    honest = [r for r in validate_vcf(truth, outbam, REF, orig_bam=BAM,
                                      thresholds=Thresholds())
              if r.kind == 'DUP']
    assert len(honest) == 1
    assert honest[0].status == OK, honest[0].observation

    overclaimed = os.path.join(workdir, 'overclaimed.vcf')
    with open(truth) as src, open(overclaimed, 'w') as dst:
        for line in src:
            dst.write(line.replace('NDUPS=3', 'NDUPS=9'))

    reports = [r for r in validate_vcf(overclaimed, outbam, REF, orig_bam=BAM,
                                       thresholds=Thresholds())
               if r.kind == 'DUP']

    assert reports[0].status == LOW_VAF, reports[0].observation
    assert 'predicted' in reports[0].observation.note

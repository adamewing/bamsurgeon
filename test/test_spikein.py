'''
Runs every case in cases.CASES and asserts the result.

This replaces test/*.sh. Those ended in a `samtools mpileup | bcftools call`
that printed calls and always exited 0, so none of them could fail; and they
all wrote into test_data/, so they overwrote each other and some depended on
having been run in a particular order.
'''

import os
import subprocess
import sys

import pysam
import pytest

from cases import CASES, Case
from conftest import (BIN, PICARD, REF, data, have, have_picard,
                      sort_and_index, truth_vcf_path)

from bamsurgeon.validate import validate_vcf, Thresholds, OK  # noqa: E402


def resolve(case):
    ''' expand test_data basenames in the argument list '''
    out = []
    for a in case.args:
        candidate = data(a)
        out.append(candidate if os.path.exists(candidate) else a)
    return out


def missing_requirements(case):
    if not have(*case.needs):
        return 'needs ' + ', '.join(case.needs)
    if case.tool != 'addsv' and not have_picard():
        return 'needs picard (BAMSURGEON_PICARD_JAR)'
    for path in (data(case.varfile), data(case.bam)):
        if not os.path.exists(path):
            return 'missing fixture %s' % os.path.basename(path)
    return None


def run_case(case, workdir):
    outbam = os.path.join(workdir, case.name + '.bam')

    cmd = [sys.executable, '-O', os.path.join(BIN, case.tool + '.py'),
           '-v', data(case.varfile), '-f', data(case.bam), '-r', REF,
           '-o', outbam, '--aligner', case.aligner, '--seed', '1234',
           '--tmpdir', os.path.join(workdir, 'tmp'),
           '--vcf', workdir + os.sep]
    cmd += resolve(case)

    # addsv realigns from FASTQ and never calls picard
    if case.tool != 'addsv' and not case.use_env_picard:
        cmd += ['--picardjar', PICARD]

    env = dict(os.environ, PYTHONPATH=BIN)
    if case.use_env_picard:
        env['BAMSURGEON_PICARD_JAR'] = PICARD

    proc = subprocess.run(cmd, cwd=workdir, env=env,
                          capture_output=True, text=True)
    return proc, outbam


def _param(case):
    marks = [pytest.mark.xfail(reason=case.xfail, strict=False)] if case.xfail else []
    return pytest.param(case, id=case.name, marks=marks)


@pytest.mark.parametrize('case', [_param(c) for c in CASES])
def test_spikein(case, workdir):
    skip = missing_requirements(case)
    if skip:
        pytest.skip(skip)

    proc, outbam = run_case(case, workdir)

    if case.expect == 'no_variants':
        # every site is excluded by the thresholds, so refusing is correct
        assert proc.returncode != 0, 'expected no successful mutations'
        assert 'no succesful mutations' in proc.stderr, proc.stderr[-2000:]
        return

    assert proc.returncode == 0, proc.stderr[-3000:]

    if case.expect == 'skipmerge':
        # the donor BAM is left unmerged for the caller to handle
        assert 'skipping merge' in proc.stderr, proc.stderr[-2000:]
        donors = [f for f in os.listdir(workdir) if f.endswith('.muts.bam')]
        assert donors, 'no unmerged donor BAM was left behind'
        return

    assert os.path.exists(outbam), proc.stderr[-3000:]
    sort_and_index(outbam)

    truth = truth_vcf_path(workdir, outbam, case.tool, data(case.varfile))
    assert os.path.exists(truth), 'no truth VCF at %s' % truth

    reports = validate_vcf(truth, outbam, REF, orig_bam=data(case.bam),
                           thresholds=Thresholds())
    assert reports, 'truth VCF has no records'

    passed = sum(1 for r in reports if r.status == OK)
    rate = passed / float(len(reports))

    failures = [(r.chrom, r.pos, r.kind, r.status, r.observation.note)
                for r in reports if r.status != OK]

    assert rate >= case.min_pass_rate, \
        '%s: %d/%d validated (need %.2f); failures: %s' % (
            case.name, passed, len(reports), case.min_pass_rate, failures)


@pytest.mark.skipif(not (have('bwa', 'samtools') and have_picard()),
                    reason='needs bwa, samtools and picard')
def test_avoidreads(workdir):
    '''
    --avoidreads drops mutations that would touch a listed read.

    The shell version of this used test_data/testregion_mut.bam as both input
    and output, and only worked if test_snv.sh had been run first to create
    it. Here it runs against the pristine input BAM and compares against the
    same spike-in without the avoid list, which is the actual claim.
    '''
    base = Case('avoid_off', 'addsnv', 'random_snvs.txt',
                args=('-n', '10', '-c', 'test_cnvlist.txt.gz'))
    avoid = Case('avoid_on', 'addsnv', 'random_snvs.txt',
                 args=('-n', '10', '-c', 'test_cnvlist.txt.gz',
                       '--avoidreads', 'test_avoid_snv1.txt'))

    counts = {}
    for case in (base, avoid):
        proc, outbam = run_case(case, workdir)
        assert proc.returncode == 0, proc.stderr[-3000:]
        truth = truth_vcf_path(workdir, outbam, case.tool, data(case.varfile))
        with open(truth) as fh:
            counts[case.name] = sum(1 for line in fh if not line.startswith('#'))

    assert counts['avoid_on'] < counts['avoid_off'], \
        'avoid list did not drop any site: %s' % counts


SV_COORDINATE_CASES = [c for c in CASES
                       if c.tool == 'addsv' and c.expect == 'variants'
                       and not c.xfail and c.aligner == 'mem']


@pytest.mark.parametrize('case', SV_COORDINATE_CASES, ids=lambda c: c.name)
def test_sv_coordinates_match_the_request(case, workdir):
    '''
    Output POS equals the requested start, and END the requested end.

    This is the property the reference-templated engine exists for. The
    assembly-based engine derived coordinates from an exonerate alignment of
    a velvet contig, so they drifted: a DEL requested at 33871043 came out at
    33870241.

    Insertions are excluded because a varfile gives an interval to place the
    insertion within, and the reported position is where it actually landed.
    '''
    skip = missing_requirements(case)
    if skip:
        pytest.skip(skip)

    from bamsurgeon.varinput import read_variants

    proc, outbam = run_case(case, workdir)
    assert proc.returncode == 0, proc.stderr[-3000:]

    requested = {}
    for req in read_variants(data(case.varfile), 'sv', default_vaf=1.0):
        if req.kind == 'INS':
            continue
        requested[(req.chrom, req.start)] = req

    truth = truth_vcf_path(workdir, outbam, case.tool, data(case.varfile))

    seen = 0
    with pysam.VariantFile(truth) as vcf:
        for rec in vcf:
            key = (rec.chrom, rec.pos)
            if key not in requested:
                continue
            req = requested[key]
            seen += 1

            if req.kind == 'BND':
                continue  # END is meaningless for a breakend

            assert rec.stop == req.end, (
                '%s %s:%d END is %d, requested %d'
                % (req.kind, rec.chrom, rec.pos, rec.stop, req.end))

    assert seen == len(requested), (
        'only %d of %d requested non-INS variants appear at their requested '
        'position in %s' % (seen, len(requested), truth))

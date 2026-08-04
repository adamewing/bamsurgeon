'''
End-to-end test of the unified entry point: SNVs, indels and SVs from one VCF.

Unlike the shell scripts in this directory, these assert. They run a spike-in
and then check the output BAM against the truth VCF the run emitted, using
the same validator a user would.

Requires the toolchain from scripts/setup_env.sh and either --picardjar or
$BAMSURGEON_PICARD_JAR.
'''

import os
import shutil
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, 'bin')

sys.path.insert(0, BIN)

from bamsurgeon.validate import validate_vcf, Thresholds, OK  # noqa: E402

REF = os.path.join(REPO, 'test_data', 'Homo_sapiens_chr22_assembly19.fasta')
BAM = os.path.join(REPO, 'test_data', 'testregion_realign.bam')
COMBINED_VCF = os.path.join(REPO, 'test_data', 'test_combined.vcf')

PICARD = os.environ.get('BAMSURGEON_PICARD_JAR')


requires_toolchain = pytest.mark.skipif(
    not (PICARD and os.path.exists(PICARD) and shutil.which('bwa')
         and shutil.which('samtools') and shutil.which('wgsim')),
    reason='needs bwa, samtools, wgsim and picard; see scripts/setup_env.sh')


def run_bamsurgeon(varfile, outbam, tmpdir, extra=()):
    ''' invoke the unified CLI as a subprocess, as a user would '''
    env = dict(os.environ, PYTHONPATH=BIN)

    cmd = [sys.executable, '-m', 'bamsurgeon.cli.main',
           '-v', varfile, '-f', BAM, '-r', REF, '-o', outbam,
           '-p', '4', '--seed', '1234', '--aligner', 'mem',
           '--picardjar', PICARD,
           '--tmpdir', os.path.join(tmpdir, 'bs.tmp'),
           '--vcf', tmpdir + os.sep]
    cmd.extend(extra)

    proc = subprocess.run(cmd, cwd=tmpdir, env=env,
                          capture_output=True, text=True)
    return proc


def sort_and_index(bam):
    subprocess.run(['samtools', 'sort', '-o', bam + '.s', bam], check=True)
    os.replace(bam + '.s', bam)
    subprocess.run(['samtools', 'index', bam], check=True)


def truth_vcf_for(tmpdir, outbam, varfile):
    stem = os.path.basename(outbam).rsplit('.bam', 1)[0]
    var = '.'.join(os.path.basename(varfile).split('.')[:-1])
    return os.path.join(tmpdir, '%s.bamsurgeon.%s.vcf' % (stem, var))


@requires_toolchain
def test_combined_snv_indel_sv_from_vcf(tmp_path):
    '''
    Two SNVs, two indels and two SVs in a single VCF, applied in one run.

    Every one must appear in the output BAM. This exercises the whole path:
    VCF parsing routes explicit alleles to the read-editing engine and
    symbolic ALTs to the SV engine, both run against the same BAM, and their
    records merge into one truth VCF.
    '''
    tmpdir = str(tmp_path)
    outbam = os.path.join(tmpdir, 'combined.bam')

    proc = run_bamsurgeon(COMBINED_VCF, outbam, tmpdir)
    assert proc.returncode == 0, proc.stderr[-3000:]
    assert os.path.exists(outbam), proc.stderr[-3000:]

    sort_and_index(outbam)

    truth = truth_vcf_for(tmpdir, outbam, COMBINED_VCF)
    assert os.path.exists(truth), 'no truth VCF at %s' % truth

    reports = validate_vcf(truth, outbam, REF, orig_bam=BAM,
                           thresholds=Thresholds())

    kinds = sorted(r.kind for r in reports)
    assert kinds == ['DEL', 'INDEL', 'INDEL', 'INV', 'SNV', 'SNV'], kinds

    failures = [(r.chrom, r.pos, r.kind, r.status)
                for r in reports if r.status != OK]
    assert not failures, 'sites did not validate: %s' % failures


@requires_toolchain
def test_small_variant_inside_an_sv_survives(tmp_path):
    '''
    An SNV inside an inversion still lands.

    This is why the unified run applies SVs first: the SV engine replaces
    every read over its interval with reads simulated from reference, so a
    small variant spiked in beforehand would be reverted. Applied afterwards,
    it edits the simulated reads instead.
    '''
    tmpdir = str(tmp_path)
    outbam = os.path.join(tmpdir, 'overlap.bam')

    # an SNV placed well inside the inversion used by the combined fixture
    varfile = os.path.join(tmpdir, 'overlap.vcf')
    _write_overlap_vcf(varfile)

    proc = run_bamsurgeon(varfile, outbam, tmpdir)
    assert proc.returncode == 0, proc.stderr[-3000:]

    # the overlap is reported rather than silently accepted
    assert 'falls inside' in proc.stderr, proc.stderr[-2000:]

    sort_and_index(outbam)

    truth = truth_vcf_for(tmpdir, outbam, varfile)
    reports = validate_vcf(truth, outbam, REF, orig_bam=BAM,
                           thresholds=Thresholds())

    snvs = [r for r in reports if r.kind == 'SNV']
    assert len(snvs) == 1, [r.kind for r in reports]
    assert snvs[0].status == OK, snvs[0].observation


def _write_overlap_vcf(path):
    import pysam
    from bamsurgeon.records import MutationRecord
    from bamsurgeon.makevcf import write_vcf

    fa = pysam.Fastafile(REF)
    pos = 34078000                      # inside 34071043-34084754
    ref = fa.fetch('22', pos - 1, pos).upper()

    write_vcf([
        MutationRecord('INV', '22', 34071043, 34084754, svlen=13711, vaf=0.5),
        MutationRecord('SNV', '22', pos, pos, ref_allele=ref,
                       alt_allele='G' if ref != 'G' else 'C', vaf=0.5),
    ], REF, path)

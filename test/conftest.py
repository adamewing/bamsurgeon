'''
Shared fixtures and helpers for the spike-in tests.

Each test gets its own tmp_path, which is the point of the conversion: the
shell scripts these replace all wrote into test_data/, so they overwrote each
other's output and some depended on having been run in a particular order.
'''

import os
import shutil
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, 'bin')
TEST_DATA = os.path.join(REPO, 'test_data')

sys.path.insert(0, BIN)

REF = os.path.join(TEST_DATA, 'Homo_sapiens_chr22_assembly19.fasta')

PICARD = os.environ.get('BAMSURGEON_PICARD_JAR')


def have(*executables):
    return all(shutil.which(e) for e in executables)


def have_picard():
    return bool(PICARD and os.path.exists(PICARD))


def data(name):
    return os.path.join(TEST_DATA, name)


def sort_and_index(bam):
    subprocess.run(['samtools', 'sort', '-o', bam + '.s', bam], check=True,
                   capture_output=True)
    os.replace(bam + '.s', bam)
    subprocess.run(['samtools', 'index', bam], check=True, capture_output=True)


def truth_vcf_path(outdir, outbam, tool, varfile):
    ''' the name each tool derives for its truth VCF '''
    stem = os.path.basename(outbam).rsplit('.bam', 1)[0]
    var = '.'.join(os.path.basename(varfile).split('.')[:-1])
    return os.path.join(outdir, '%s.%s.%s.vcf' % (stem, tool, var))


@pytest.fixture
def workdir(tmp_path):
    ''' an isolated output directory, one per test '''
    return str(tmp_path)

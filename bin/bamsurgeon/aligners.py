#!/usr/bin/env python


import os
import subprocess
import shutil
from uuid import uuid4

from bamsurgeon.common import bamtofastq, double_fastq_to_interleaved


import logging
FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

SUPPORTED_ALIGNERS = frozenset(['backtrack', 'mem', 'novoalign', 'gsnap', 'STAR', 'bowtie2', 'tmap', 'bwakit', 'minimap2'])
INTERLEAVED_FASTQ_ALIGNERS = frozenset(['minimap2'])


def _checkoptions(name, options):
    ''' checks if necessary options have been specified '''
    if name == 'novoalign':
        if 'novoref' not in options:
            raise ValueError("ERROR '--aligner novoalign' requires '--alignopts novoref:path_to_novoalign_reference' to be set\n")
    elif name == 'gsnap':
        if 'gsnaprefdir' not in options or 'gsnaprefname' not in options:
            raise ValueError("ERROR '--aligner gsnap' requires '--alignopts gsnaprefdir:GSNAP_reference_dir,gsnaprefname:GSNAP_ref_name\n")
    elif name == 'STAR':
        if 'STARrefdir' not in options:
            raise ValueError("ERROR '--aligner STAR' requires '--alignopts STARrefdir:/path/to/STAR_reference_dir\n")
    elif name == 'bowtie2':
        if 'bowtie2ref' not in options:
            raise ValueError("ERROR '--aligner bowtie2' requires '--alignopts bowtie2ref:bowtie2_ref_basepath\n")
    elif name == 'minimap2':
        if 'x' not in options:
            raise ValueError("ERROR '--aligner minimap2' requires '--alignopts x:map-pb or x:map-ont\n")
    elif name not in SUPPORTED_ALIGNERS:
        raise ValueError("ERROR unsupported aligner: " + name)


def remap_bam(name, bamfn, fastaref, options, mutid='null', threads=1, paired=True, picardjar=None, deltmp=True):
    ''' remap bam file with supported alignment method. "options" param is a dict of aligner-specific required options '''
    fastqs = bamtofastq(bamfn, picardjar, threads=threads, paired=paired, twofastq=(name not in INTERLEAVED_FASTQ_ALIGNERS))

    if len(fastqs) == 2:
        fq1, fq2 = fastqs
    else:
        fq1 = fastqs[0]
        fq2 = None

    remap_fastq(name, fq1, fastaref, bamfn, options, fq2=fq2, mutid=mutid, threads=threads, deltmp=deltmp)


def remap_fastq(name, fq1, fastaref, outbam, options, fq2=None, mutid='null', threads=1, deltmp=True):
    ''' remap bam file with supported alignment method. "options" param is a dict of aligner-specific required options '''

    _checkoptions(name, options)

    basefn = "tmp." + str(uuid4())
    sam_out = basefn + '.sam'
    sort_tmp = basefn + '.sort.bam'

    logger.info(mutid + " aligning " + fq1 + ',' + fq2 + " with " + name)
    if name == 'backtrack':
        _run_backtrack(fq1, fastaref, sam_out, fq2=fq2, threads=threads)
    elif name == 'mem':
        _run_bwa_mem(fq1, fastaref, sam_out, fq2=fq2, threads=threads)
    elif name == 'novoalign':
        _run_novoalign(fq1, sam_out, options['novoref'], fq2=fq2)
    elif name == 'gsnap':
        _run_gsnap(fq1, sam_out, options['gsnaprefdir'], options['gsnaprefname'], fq2=fq2, threads=threads)
    elif name == 'STAR':
        _run_STAR(fq1, sam_out, options['STARrefdir'], fq2=fq2)
    elif name == 'bowtie2':
        _run_bowtie2(fq1, options['bowtie2ref'], sam_out, fq2=fq2)
    elif name == 'tmap':
        _run_tmap(fq1, fastaref, sam_out, fq2=fq2, threads=threads)
    elif name == 'bwakit':
        _run_bwakit(fq1, fastaref, sam_out, fq2=fq2, threads=threads)
    elif name == 'minimap2':
        # minimap2 requires interleaved fastq
        if fq2 is not None:
            double_fastq_to_interleaved(fq1, fq2, fq1 + '.interleaved')
            os.remove(fq1)
            os.remove(fq2)
            fq1 = fq1 + '.interleaved'
        _run_minimap2(fq1, fastaref, sam_out, options['x'], threads=threads)
    else:
        raise ValueError("ERROR: " + name + " is not a supported aligner for remapping.")

    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-T', sort_tmp, '-o', outbam, sam_out]
    idx_cmd = ['samtools', 'index', outbam]
    subprocess.check_call(sort_cmd)
    subprocess.check_call(idx_cmd)

    if deltmp:
        os.remove(fq1)
        os.remove(sam_out)
        if fq2 is not None:
            os.remove(fq2)

# Aligner functions


def _run_backtrack(fq1, fastaref, sam_out, fq2=None, threads=1):
    sai1fn = sam_out + ".1.sai"
    sai2fn = sam_out + ".2.sai"

    sai1args = ['bwa', 'aln', '-q', '5', '-l', '32', '-k', '2', '-t', str(threads), '-o', '1', '-f', sai1fn, fastaref, fq1]
    subprocess.check_call(sai1args)
    if fq2 is not None:
        # paired-end
        sai2args = ['bwa', 'aln', '-q', '5', '-l', '32', '-k', '2', '-t', str(threads), '-o', '1', '-f', sai2fn, fastaref, fq2]
        subprocess.check_call(sai2args)
        samargs = ['bwa', 'sampe', '-P', '-f', sam_out, fastaref, sai1fn, sai2fn, fq1, fq2]
        subprocess.check_call(samargs)
        os.remove(sai2fn)
    else:
        # single-end
        samargs = ['bwa', 'samse', '-f', sam_out, fastaref, sai1fn, fq1]
        subprocess.check_call(samargs)
    os.remove(sai1fn)


def _run_bwa_mem(fq1, fastaref, sam_out, fq2=None, threads=1):
    if fq2 is None:
        sam_cmd = ['bwa', 'mem', '-t', str(threads), '-M', '-Y', fastaref, fq1]
    else:
        sam_cmd = ['bwa', 'mem', '-t', str(threads), '-M', '-Y', fastaref, fq1, fq2]
    with open(sam_out, 'w') as sam:
        subprocess.check_call(sam_cmd, stdout=sam)


def _run_novoalign(fq1, sam_out, novoref, fq2=None):
    sam_cmd = ['novoalign', '-F', 'STDFQ', '-r', 'Random', '-d', novoref, '-oSAM', '-f', fq1]
    if fq2 is not None:
        sam_cmd.append(fq2)
    with open(sam_out, 'w') as sam:
        subprocess.check_call(sam_cmd, stdout=sam)


def _run_gsnap(fq1, sam_out, gsnaprefdir, gsnaprefname, fq2=None, threads=1):
    sam_cmd = ['gsnap', '-D', gsnaprefdir, '-d', gsnaprefname, '-t', str(threads), '--quality-protocol=sanger',
               '-M', '2', '-n', '10', '-B', '2', '-i', '1', '--pairmax-dna=1000', '--terminal-threshold=1000',
               '--gmap-mode=none', '--clip-overlap', '-A', 'sam', '-a', 'paired', fq1]
    if fq2 is not None:
        sam_cmd.append(fq2)
    with open(sam_out, 'w') as sam:
        subprocess.check_call(sam_cmd, stdout=sam)


def _run_STAR(fq1, sam_out, STARrefdir, fq2=None):
    sam_cmd = ['STAR', '--genomeLoad', 'LoadAndKeep', '--genomeDir', STARrefdir, '--outFileNamePrefix', sam_out, '--readFilesIn', fq1]
    if fq2 is not None:
        sam_cmd.append(fq2)
    subprocess.check_call(sam_cmd)
    shutil.move(sam_out + 'Aligned.out.sam', sam_out)


def _run_bowtie2(fq1, bowtie2ref, sam_out, fq2=None):
    sam_cmd = ['bowtie2', '-x', bowtie2ref, '-S', sam_out, '-1', fq1]
    if fq2 is not None:
        sam_cmd += ['-2', fq2]
    subprocess.check_call(sam_cmd)


def _run_tmap(fq1, fastaref, sam_out, fq2=None, threads=1):
    if fq2 is not None:
        raise ValueError('tmap is not supported for paired-end reads.')
    sam_cmd = ['tmap', 'mapall', '-f', fastaref, '-r', fq1, '-n', str(threads), '-v', '-u', '-o', '0', 'stage1', 'map4']
    with open(sam_out, 'w') as sam:
        subprocess.check_call(sam_cmd, stdout=sam)


def _run_bwakit(fq1, fastaref, sam_out, fq2=None, threads=1):
    sam_cmd = ['run-bwamem', '-t', str(threads), '-o', sam_out, '-H', fastaref, fq1]
    if fq2 is not None:
        sam_cmd.append(fq2)
    # It is not a SAM file, but it should work
    shutil.move(sam_out + '.aln.bam', sam_out)


def _run_minimap2(fq1, fastaref, sam_out, x, threads=1):
    sam_cmd = ['minimap2', '-ax', x, '-t', str(threads), fastaref, fq1]
    with open(sam_out, 'w') as sam:
        subprocess.check_call(sam_cmd, stdout=sam)

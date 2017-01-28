#!/usr/bin/env python


import os
import sys
import subprocess

from common import *
from shutil import move
from re import sub
from uuid import uuid4


#
# Remapping functions bam --> bam
#


supported_aligners_bam   = ['backtrack', 'mem', 'novoalign', 'gsnap', 'STAR', 'bowtie2', 'tmap','bwakit']
supported_aligners_fastq = ['backtrack', 'mem', 'novoalign']

def checkoptions(name, options, picardjar, sv=False):
    ''' checks if necessary options have been specified '''

    if sv:
        if name not in supported_aligners_fastq:
            raise ValueError("ERROR\tunsupported aligner: " + name + "\n")
    else:
        if name not in supported_aligners_bam:
            raise ValueError("ERROR\tunsupported aligner: " + name + "\n")

    if name != 'backtrack' and not sv:
        if picardjar is None:
            raise ValueError("ERROR\t'--aligner " + name + "' requires '--picardjar' to be specified\n")

    if name == 'novoalign':
        if 'novoref' not in options:
            raise ValueError("ERROR\t'--aligner novoalign' requires '--alignopts novoref:path_to_novoalign_reference' to be set\n")

    if name == 'gsnap':
        if 'gsnaprefdir' not in options or 'gsnaprefname' not in options:
            raise ValueError("ERROR\t'--aligner gsnap' requires '--alignopts gsnaprefdir:GSNAP_reference_dir,gsnaprefname:GSNAP_ref_name\n")

    if name == 'STAR':
        if 'STARrefdir' not in options:
            raise ValueError("ERROR\t'--aligner STAR' requires '--alignopts STARrefdir:/path/to/STAR_reference_dir\n")

    if name == 'bowtie2':
        if 'bowtie2ref' not in options:
            raise ValueError("ERROR\t'--aligner bowtie2' requires '--alignopts bowtie2ref:bowtie2_ref_basepath\n")


def remap_bam(name, bamfn, fastaref, options, mutid='null', threads=1, paired=True, picardjar=None):
    ''' remap bam file with supported alignment method. "options" param is a dict of aligner-specific required options '''

    checkoptions(name, options, picardjar)

    assert os.path.exists(bamfn), 'cannot locate bam: %s' % bamfn

    if picardjar is not None:
        assert os.path.exists(picardjar), 'cannot locate picard.jar: %s' % picardjar

    if name == 'backtrack':
        remap_backtrack_bam(bamfn, threads, fastaref, mutid=mutid, paired=paired)

    if name == 'mem':
        remap_bwamem_bam(bamfn, threads, fastaref, picardjar, mutid=mutid, paired=paired)

    if name == 'bwakit':
        remap_bwakit_bam(bamfn, threads, fastaref, picardjar, mutid=mutid, paired=paired)

    if name == 'novoalign':
        remap_novoalign_bam(bamfn, threads, fastaref, picardjar, options['novoref'], mutid=mutid, paired=paired)

    if name == 'gsnap':
        remap_gsnap_bam(bamfn, threads, fastaref, picardjar, options['gsnaprefdir'], options['gsnaprefname'], mutid=mutid, paired=paired)

    if name == 'STAR':
        remap_STAR_bam(bamfn, threads, fastaref, picardjar, options['STARrefdir'], mutid=mutid, paired=paired)

    if name == 'bowtie2':
        remap_bowtie2_bam(bamfn, threads, fastaref, picardjar, options['bowtie2ref'], mutid=mutid, paired=paired)

    if name == 'tmap':
        remap_tmap_bam(bamfn, threads, fastaref, picardjar, mutid=mutid, paired=paired)


def remap_backtrack_bam(bamfn, threads, fastaref, mutid='null', paired=True):
    """ call bwa/samtools to remap .bam
    """
    if paired:
        sai1fn = bamfn + ".1.sai"
        sai2fn = bamfn + ".2.sai"
        samfn  = bamfn + ".sam"
        refidx = fastaref + ".fai"

        sai1args = ['bwa', 'aln', '-q', '5', '-l', '32', '-k', '3', '-t', str(threads), '-o', '1', '-f', sai1fn, '-b1', fastaref, bamfn]
        sai2args = ['bwa', 'aln', '-q', '5', '-l', '32', '-k', '3', '-t', str(threads), '-o', '1', '-f', sai2fn, '-b2', fastaref, bamfn]
        samargs  = ['bwa', 'sampe', '-P', '-f', samfn, fastaref, sai1fn, sai2fn, bamfn, bamfn]
        bamargs  = ['samtools', 'view', '-bt', refidx, '-o', bamfn, samfn] 

        print "INFO\t" + now() + "\t" + mutid + "\tmapping 1st end, cmd: " + " ".join(sai1args)
        subprocess.call(sai1args)
        print "INFO\t" + now() + "\t" + mutid + "\tmapping 2nd end, cmd: " + " ".join(sai2args)
        subprocess.call(sai2args)
        print "INFO\t" + now() + "\t" + mutid + "\tpairing ends, building .sam, cmd: " + " ".join(samargs)
        subprocess.call(samargs)
        print "INFO\t" + now() + "\t" + mutid + "\tsam --> bam, cmd: " + " ".join(bamargs)
        subprocess.call(bamargs)

        sortfn = bamfn + ".sort.bam"
        sortargs = ['samtools','sort','-m','10000000000', '-T', sortfn ,'-o',sortfn, bamfn]
        print "INFO\t" + now() + "\t" + mutid + "\tsorting, cmd: " + " ".join(sortargs)
        subprocess.call(sortargs)
        os.rename(sortfn, bamfn)

        indexargs = ['samtools','index',bamfn]
        print "INFO\t" + now() + "\t" + mutid + "\tindexing, cmd: " + " ".join(indexargs)
        subprocess.call(indexargs)

        # cleanup
        os.remove(sai1fn)
        os.remove(sai2fn)
        os.remove(samfn)

    else:
        saifn = bamfn + ".sai"
        samfn  = bamfn + ".sam"
        refidx = fastaref + ".fai"

        saiargs = ['bwa', 'aln', '-q', '5', '-l', '32', '-k', '3', '-t', str(threads), '-o', '1', '-f', saifn, '-b1', fastaref, bamfn]
        samargs  = ['bwa', 'samse', '-f', samfn, fastaref, saifn, bamfn]
        bamargs  = ['samtools', 'view', '-bt', refidx, '-o', bamfn, samfn] 

        print "INFO\t" + now() + "\t" + mutid + "\tmapping, cmd: " + " ".join(saiargs)
        subprocess.call(saiargs)
        print "INFO\t" + now() + "\t" + mutid + "\tpairing ends, building .sam, cmd: " + " ".join(samargs)
        subprocess.call(samargs)
        print "INFO\t" + now() + "\t" + mutid + "\tsam --> bam, cmd: " + " ".join(bamargs)
        subprocess.call(bamargs)

        sortfn = bamfn + ".sort.bam"
        sortargs = ['samtools','sort','-m','10000000000', '-T', sortfn, '-o',sortfn, bamfn]
        print "INFO\t" + now() + "\t" + mutid + "\tsorting, cmd: " + " ".join(sortargs)
        subprocess.call(sortargs)
        os.rename(sortfn,bamfn)

        indexargs = ['samtools','index',bamfn]
        print "INFO\t" + now() + "\t" + mutid + "\tindexing, cmd: " + " ".join(indexargs)
        subprocess.call(indexargs)

        # cleanup
        os.remove(saifn)
        os.remove(samfn)


def remap_bwamem_bam(bamfn, threads, fastaref, picardjar, mutid='null', paired=True, deltmp=True):
    """ call bwa mem and samtools to remap .bam
    """
    assert os.path.exists(picardjar)
    assert bamreadcount(bamfn) > 0
    if paired:
        assert bamreadcount(bamfn) > 1 

    sam_out  = bamfn + '.realign.sam'
    sort_out = bamfn + '.realign.sorted.bam'

    print "INFO\t" + now() + "\t" + mutid + "\tconverting " + bamfn + " to fastq\n"
    fastq = bamtofastq(bamfn, picardjar, threads=threads, paired=paired)[0]

    sam_cmd = []

    if paired:
        sam_cmd  = ['bwa', 'mem', '-t', str(threads), '-M', '-p', fastaref, fastq] # interleaved
    else:
        sam_cmd  = ['bwa', 'mem', '-t', str(threads), '-M', fastaref, fastq] # single-end

    assert len(sam_cmd) > 0

    bam_cmd  = ['samtools', 'view', '-bt', fastaref + '.fai', '-o', bamfn, sam_out]
    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-m', '10000000000', '-T', sort_out, '-o', sort_out, bamfn]
    idx_cmd  = ['samtools', 'index', bamfn]

    print "INFO\t" + now() + "\t" + mutid + "\taligning " + fastq + " with bwa mem\n"
    with open(sam_out, 'w') as sam:
        p = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE)
        for line in p.stdout:
            sam.write(line)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\twriting " + sam_out + " to BAM...\n")
    subprocess.call(bam_cmd)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tdeleting SAM: " + sam_out + "\n")
    if deltmp: os.remove(sam_out)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tsorting output: " + ' '.join(sort_cmd) + "\n")
    subprocess.call(sort_cmd)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremove original bam:" + bamfn + "\n")
    if deltmp: os.remove(bamfn)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\trename sorted bam: " + sort_out + " to original name: " + bamfn + "\n")
    move(sort_out, bamfn)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tindexing: " + ' '.join(idx_cmd) + "\n")
    subprocess.call(idx_cmd)

    # check if BAM readcount looks sane
    if bamreadcount(bamfn) < fastqreadcount(fastq): 
        raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tbam readcount < fastq readcount, alignment sanity check failed!\n")

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fastq + "\n")
    if deltmp: os.remove(fastq)


def remap_bwakit_bam(bamfn, threads, fastaref, picardjar, mutid='null', paired=True, deltmp=True):
    """ call bwa kit and samtools to remap .bam
    """
    assert os.path.exists(picardjar)
    assert bamreadcount(bamfn) > 0
    if paired:
        assert bamreadcount(bamfn) > 1 
        print 'paird'

    bam_out  = sub('.bam$','',bamfn) + '.aln.bam'
    sort_out = bamfn + '.realign.sorted.bam'
    out_prefix = sub('.bam$','',bamfn)

    print "INFO\t" + now() + "\t" + mutid + "\tconverting " + bamfn + " to fastq\n"
    kit_cmd1 = []
    fastq = bamtofastq(bamfn, picardjar, threads=threads, paired=paired, twofastq=True)
    if paired:
        kit_cmd1  = ['run-bwamem', '-t', '4', '-o', out_prefix, '-H', fastaref, fastq[0], fastq[1]]
    else:
        kit_cmd1  = ['run-bwamem', '-t', '4', '-o', out_prefix, '-H', fastaref, fastq[0]]

    print kit_cmd1
    assert len(kit_cmd1) > 0
    print "INFO\t" + now() + "\t" + mutid + "\taligning " + ', '.join(fastq)  + " with bwa kit\n"
    kit_cmd2 = ['sh']
    #robust way : open two processes and pipe them together
    process_runkit = subprocess.Popen(kit_cmd1, stdout=subprocess.PIPE, shell=False)
    process_sh = subprocess.Popen(kit_cmd2, stdin=process_runkit.stdout, stdout=subprocess.PIPE, shell=False)
    process_runkit.stdout.close()
    process_sh.communicate()[0]

    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-m', '10000000000', '-T', sort_out, '-o', sort_out, bamfn]
    idx_cmd  = ['samtools', 'index', bamfn]

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\trename aligned bam: " + bam_out + " to original name: " + bamfn + "\n")
    move(bam_out, bamfn)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tsorting output: " + ' '.join(sort_cmd) + "\n")
    subprocess.call(sort_cmd)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremove original bam:" + bamfn + "\n")
    if deltmp: os.remove(bamfn)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\trename sorted bam: " + sort_out + " to original name: " + bamfn + "\n")
    move(sort_out, bamfn)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tindexing: " + ' '.join(idx_cmd) + "\n")
    subprocess.call(idx_cmd)

    # check if BAM readcount looks sane
    if paired:
        if bamreadcount(bamfn) < fastqreadcount(fastq[0]) + fastqreadcount(fastq[1]): 
            raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tbam readcount < fastq readcount, alignment sanity check failed!\n")
    else:
        if bamreadcount(bamfn) < fastqreadcount(fastq[0]): 
            raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tbam readcount < fastq readcount, alignment sanity check failed!\n")

    if paired:
        for fq in fastq:
            sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fq + "\n")
            os.remove(fq)
    else:
        sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fastq[0] + "\n")
        os.remove(fastq[0])


def remap_novoalign_bam(bamfn, threads, fastaref, picardjar, novoref, mutid='null', paired=True):
    """ call novoalign and samtools to remap .bam
    """
    assert os.path.exists(picardjar)
    assert os.path.exists(novoref)
    assert bamreadcount(bamfn) > 0 

    sam_out  = bamfn + '.realign.sam'
    sort_out = bamfn + '.realign.sorted.bam'

    print "INFO\t" + now() + "\t" + mutid + "\tconverting " + bamfn + " to fastq\n"
    fastq = bamtofastq(bamfn, picardjar, threads=threads, paired=paired, twofastq=True)

    sam_cmd = []

    if paired:
        sam_cmd  = ['novoalign', '--mmapOff', '-F', 'STDFQ', '-f', fastq[0], fastq[1], '-r', 'Random', '-d', novoref, '-oSAM'] # interleaved
        #sam_cmd  = ['novoalign', '-F', 'STDFQ', '-f', fastq[0], fastq[1], '-r', 'Random', '-d', novoref, '-oSAM'] # uncomment for unlicensed
    else:
        sam_cmd  = ['novoalign', '--mmapOff', '-F', 'STDFQ', '-f', fastq[0], '-r', 'Random', '-d', novoref, '-oSAM'] # interleaved
        #sam_cmd  = ['novoalign', '-F', 'STDFQ', '-f', fastq, '-r', 'Random', '-d', novoref, '-oSAM'] # uncomment for unlicensed

    assert len(sam_cmd) > 0

    bam_cmd  = ['samtools', 'view', '-bt', fastaref + '.fai', '-o', bamfn, sam_out]
    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-m', '10000000000', '-T', sort_out, '-o', sort_out, bamfn]
    idx_cmd  = ['samtools', 'index', bamfn]

    print "INFO\t" + now() + "\t" + mutid + "\taligning " + str(fastq) + " with novoalign\n"
    with open(sam_out, 'w') as sam:
        p = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE)
        for line in p.stdout:
            sam.write(line)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\twriting " + sam_out + " to BAM...\n")
    subprocess.call(bam_cmd)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tdeleting SAM: " + sam_out + "\n")
    os.remove(sam_out)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tsorting output: " + ' '.join(sort_cmd) + "\n")
    subprocess.call(sort_cmd)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremove original bam:" + bamfn + "\n")
    os.remove(bamfn)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\trename sorted bam: " + sort_out + " to original name: " + bamfn + "\n")
    move(sort_out, bamfn)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tindexing: " + ' '.join(idx_cmd) + "\n")
    subprocess.call(idx_cmd)

    # check if BAM readcount looks sane
    if paired:
        if bamreadcount(bamfn) < fastqreadcount(fastq[0]) + fastqreadcount(fastq[1]): 
            raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tbam readcount < fastq readcount, alignment sanity check failed!\n")
    else:
        if bamreadcount(bamfn) < fastqreadcount(fastq[0]): 
            raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tbam readcount < fastq readcount, alignment sanity check failed!\n")

    if paired:
        for fq in fastq:
            sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fq + "\n")
            os.remove(fq)
    else:
        sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fastq[0] + "\n")
        os.remove(fastq[0])


def remap_gsnap_bam(bamfn, threads, fastaref, picardjar, gsnaprefdir, gsnaprefname, mutid='null', paired=True):
    """ call gsnap and samtools to remap .bam
    """
    assert os.path.exists(picardjar)
    assert os.path.exists(gsnaprefdir)
    assert bamreadcount(bamfn) > 0 

    sam_out  = bamfn + '.realign.sam'
    sort_out = bamfn + '.realign.sorted.bam'

    print "INFO\t" + now() + "\t" + mutid + "\tconverting " + bamfn + " to fastq\n"
    fastq = bamtofastq(bamfn, picardjar, threads=threads, paired=paired, twofastq=True)

    sam_cmd = []

    if paired:
        sam_cmd = ['gsnap', '-D', gsnaprefdir, '-d', gsnaprefname, '-t', str(threads), '--quality-protocol=sanger', 
                   '-M', '2', '-n', '10', '-B', '2', '-i', '1', '--pairmax-dna=1000', '--terminal-threshold=1000', 
                   '--gmap-mode=none', '--clip-overlap', '-A', 'sam', '-a', 'paired', fastq[0], fastq[1]]
    else:
        sam_cmd = ['gsnap', '-D', gsnaprefdir, '-d', gsnaprefname, '-t', str(threads), '--quality-protocol=sanger', 
                   '-M', '2', '-n', '10', '-B', '2', '-i', '1', '--terminal-threshold=1000', '--gmap-mode=none', 
                   '--clip-overlap', '-A', 'sam', fastq[0]]

    assert len(sam_cmd) > 0

    bam_cmd  = ['samtools', 'view', '-bt', fastaref + '.fai', '-o', bamfn, sam_out]
    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-m', '10000000000', '-T', sort_out, '-o', sort_out, bamfn]
    idx_cmd  = ['samtools', 'index', bamfn]

    print "INFO\t" + now() + "\t" + mutid + "\taligning " + str(fastq) + " with gsnap\n"
    with open(sam_out, 'w') as sam:
        p = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE)
        for line in p.stdout:
            sam.write(line)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\twriting " + sam_out + " to BAM...\n")
    subprocess.call(bam_cmd)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tdeleting SAM: " + sam_out + "\n")
    os.remove(sam_out)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tsorting output: " + ' '.join(sort_cmd) + "\n")
    subprocess.call(sort_cmd)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremove original bam:" + bamfn + "\n")
    os.remove(bamfn)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\trename sorted bam: " + sort_out + " to original name: " + bamfn + "\n")
    move(sort_out, bamfn)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tindexing: " + ' '.join(idx_cmd) + "\n")
    subprocess.call(idx_cmd)

    # check if BAM readcount looks sane
    if paired:
        if bamreadcount(bamfn) < fastqreadcount(fastq[0]) + fastqreadcount(fastq[1]): 
            raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tbam readcount < fastq readcount, alignment sanity check failed!\n")
    else:
        if bamreadcount(bamfn) < fastqreadcount(fastq[0]): 
            raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tbam readcount < fastq readcount, alignment sanity check failed!\n")

    if paired:
        for fq in fastq:
            sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fq + "\n")
            os.remove(fq)
    else:
        sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fastq[0] + "\n")
        os.remove(fastq[0])

def remap_STAR_bam(bamfn, threads, fastaref, picardjar, STARrefdir, mutid='null', paired=True):
    """ call gsnap and samtools to remap .bam
    """
    assert os.path.exists(picardjar)
    assert os.path.exists(STARrefdir)
    assert bamreadcount(bamfn) > 0 

    sam_out  = bamfn + 'Aligned.out.sam'
    sort_out = bamfn + '.realign.sorted.bam'

    print "INFO\t" + now() + "\t" + mutid + "\tconverting " + bamfn + " to fastq\n"
    fastq = bamtofastq(bamfn, picardjar, threads=threads, paired=paired, twofastq=True)

    sam_cmd = []

    sam_cmd = ['STAR', '--genomeLoad', 'LoadAndKeep', '--genomeDir', STARrefdir, '--outFileNamePrefix', bamfn, '--readFilesIn'] + fastq

    assert len(sam_cmd) > 0

    bam_cmd  = ['samtools', 'view', '-b', '-o', bamfn, sam_out]
    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-m', '10000000000', '-T', sort_out, '-o', sort_out, bamfn]
    idx_cmd  = ['samtools', 'index', bamfn]

    print "INFO\t" + now() + "\t" + mutid + "\taligning " + str(fastq) + " with STAR\n"

    p = subprocess.call(sam_cmd)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\twriting " + sam_out + " to BAM...\n")
    subprocess.call(bam_cmd)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tdeleting SAM: " + sam_out + "\n")
    os.remove(sam_out)
    os.remove(bamfn + "Log.final.out")
    os.remove(bamfn + "Log.out")
    os.remove(bamfn + "Log.progress.out")
    os.remove(bamfn + "SJ.out.tab")

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tsorting output: " + ' '.join(sort_cmd) + "\n")
    subprocess.call(sort_cmd)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremove original bam:" + bamfn + "\n")
    os.remove(bamfn)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\trename sorted bam: " + sort_out + " to original name: " + bamfn + "\n")
    move(sort_out, bamfn)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tindexing: " + ' '.join(idx_cmd) + "\n")
    subprocess.call(idx_cmd)

    # check if BAM readcount looks sane
    # STAR may not report unmapped reads in BAM

    if paired:
        if bamreadcount(bamfn) < .5*(fastqreadcount(fastq[0]) + fastqreadcount(fastq[1])): 
            raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tbam readcount < fastq readcount, alignment sanity check failed!\n")
    else:
        if bamreadcount(bamfn) < .5*fastqreadcount(fastq[0]): 
            raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tbam readcount < fastq readcount, alignment sanity check failed!\n")

    if paired:
        for fq in fastq:
            sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fq + "\n")
            os.remove(fq)
    else:
        sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fastq[0] + "\n")
        os.remove(fastq[0])

def remap_bowtie2_bam(bamfn, threads, fastaref, picardjar, bowtie2ref, mutid='null', paired=True):
    """ call bowtie2 and samtools to remap .bam
    """

    assert bamreadcount(bamfn) > 0 

    sam_out  = bamfn + '.realign.sam'
    sort_out = bamfn + '.realign.sorted.bam'

    print "INFO\t" + now() + "\t" + mutid + "\tconverting " + bamfn + " to fastq\n"
    fastq = bamtofastq(bamfn, picardjar, threads=threads, paired=paired, twofastq=True)

    sam_cmd = []

    if paired:
        sam_cmd = ['bowtie2', '-x', bowtie2ref, '-1', fastq[0], '-2', fastq[1], '-S', sam_out]
    else:
        sam_cmd = ['bowtie2', '-x', bowtie2ref, '-U', fastq[0], '-S', sam_out]

    assert len(sam_cmd) > 0

    bam_cmd  = ['samtools', 'view', '-bt', fastaref + '.fai', '-o', bamfn, sam_out]
    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-m', '10000000000', '-T', sort_out, '-o', sort_out, bamfn]
    idx_cmd  = ['samtools', 'index', bamfn]

    print "INFO\t" + now() + "\t" + mutid + "\taligning " + str(fastq) + " with bowtie2\n"
    with open(sam_out, 'w') as sam:
        p = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE)
        for line in p.stdout:
            sam.write(line)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\twriting " + sam_out + " to BAM...\n")
    subprocess.call(bam_cmd)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tdeleting SAM: " + sam_out + "\n")
    os.remove(sam_out)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tsorting output: " + ' '.join(sort_cmd) + "\n")
    subprocess.call(sort_cmd)


    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremove original bam:" + bamfn + "\n")
    os.remove(bamfn)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\trename sorted bam: " + sort_out + " to original name: " + bamfn + "\n")
    move(sort_out, bamfn)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tindexing: " + ' '.join(idx_cmd) + "\n")
    subprocess.call(idx_cmd)

    # check if BAM readcount looks sane
    if paired:
        if bamreadcount(bamfn) < fastqreadcount(fastq[0]) + fastqreadcount(fastq[1]): 
            raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tbam readcount < fastq readcount, alignment sanity check failed!\n")
    else:
        if bamreadcount(bamfn) < fastqreadcount(fastq[0]): 
            raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tbam readcount < fastq readcount, alignment sanity check failed!\n")

    if paired:
        for fq in fastq:
            sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fq + "\n")
            os.remove(fq)
    else:
        sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fastq[0] + "\n")
        os.remove(fastq[0])


def remap_tmap_bam(bamfn, threads, fastaref, picardjar, mutid='null', paired=False):
    """ call bowtie2 and samtools to remap .bam
    """

    assert bamreadcount(bamfn) > 0 

    sam_out  = bamfn + '.realign.sam'
    sort_out = bamfn + '.realign.sorted.bam'

    print "INFO\t" + now() + "\t" + mutid + "\tconverting " + bamfn + " to fastq\n"
    fastq = bamtofastq(bamfn, picardjar, threads=threads, paired=paired, twofastq=True)

    sam_cmd = []

    if paired:
        raise ValueError('tmap only supported in --single mode.')
    else:
        sam_cmd = ['tmap', 'mapall', '-f', fastaref, '-r', fastq[0], '-n', str(threads), '-v', '-u', '-o', '0', 'stage1', 'map4']

    assert len(sam_cmd) > 0

    bam_cmd  = ['samtools', 'view', '-bt', fastaref + '.fai', '-o', bamfn, sam_out]
    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-m', '10000000000', '-T', sort_out, '-o',  sort_out, bamfn]
    idx_cmd  = ['samtools', 'index', bamfn]

    print "INFO\t" + now() + "\t" + mutid + "\taligning " + str(fastq) + " with tmap\n"
    with open(sam_out, 'w') as sam:
        p = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE)
        for line in p.stdout:
            sam.write(line)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\twriting " + sam_out + " to BAM...\n")
    subprocess.call(bam_cmd)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tdeleting SAM: " + sam_out + "\n")
    os.remove(sam_out)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tsorting output: " + ' '.join(sort_cmd) + "\n")
    subprocess.call(sort_cmd)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremove original bam:" + bamfn + "\n")
    os.remove(bamfn)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\trename sorted bam: " + sort_out + " to original name: " + bamfn + "\n")
    move(sort_out, bamfn)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tindexing: " + ' '.join(idx_cmd) + "\n")
    subprocess.call(idx_cmd)

    # check if BAM readcount looks sane
    if paired:
        if bamreadcount(bamfn) < fastqreadcount(fastq[0]) + fastqreadcount(fastq[1]): 
            raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tbam readcount < fastq readcount, alignment sanity check failed!\n")
    else:
        if bamreadcount(bamfn) < fastqreadcount(fastq[0]): 
            raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tbam readcount < fastq readcount, alignment sanity check failed!\n")

    if paired:
        for fq in fastq:
            sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fq + "\n")
            os.remove(fq)
    else:
        sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fastq[0] + "\n")
        os.remove(fastq[0])


#
# Remapping functions paired fastq --> BAM
#


def remap_fastq(name, fq1, fq2, fastaref, outbam, options, mutid='null', threads=1, deltmp=True):
    ''' remap bam file with supported alignment method. "options" param is a dict of aligner-specific required options '''

    checkoptions(name, options, None, sv=True)

    if name == 'backtrack':
        return remap_backtrack_fastq(fq1, fq2, threads, fastaref, outbam, deltmp=deltmp, mutid=mutid)


    if name == 'mem':
        return remap_bwamem_fastq(fq1, fq2, threads, fastaref, outbam, deltmp=deltmp, mutid=mutid)


    if name == 'novoalign':
        return remap_novoalign_fastq(fq1, fq2, threads, fastaref, options['novoref'], outbam, deltmp=deltmp, mutid=mutid)


def remap_bwamem_fastq(fq1, fq2, threads, fastaref, outbam, deltmp=True, mutid='null'):
    """ call bwa mem and samtools to remap .bam
    """

    basefn   = "bwatmp." + str(uuid4())
    sam_out  = basefn + '.sam'
    sort_out = basefn + '.sorted.bam'

    sam_cmd  = ['bwa', 'mem', '-t', str(threads), '-M', '-Y', fastaref, fq1, fq2]

    bam_cmd  = ['samtools', 'view', '-bt', fastaref + '.fai', '-o', outbam, sam_out]
    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-m', '10000000000', '-T', sort_out, '-o',  sort_out, outbam]
    idx_cmd  = ['samtools', 'index', outbam]

    print "INFO\t" + now() + "\t" + mutid + "\taligning " + fq1 + ',' + fq2 + " with bwa mem"
    with open(sam_out, 'w') as sam:
        p = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE)
        for line in p.stdout:
            sam.write(line)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\twriting " + sam_out + " to BAM...\n")
    subprocess.call(bam_cmd)

    if deltmp:
        sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tdeleting SAM: " + sam_out + "\n")
        os.remove(sam_out)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tsorting output: " + ' '.join(sort_cmd) + "\n")
    subprocess.call(sort_cmd)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremove original bam:" + outbam + "\n")
    os.remove(outbam)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\trename sorted bam: " + sort_out + " to original name: " + outbam + "\n")
    move(sort_out, outbam)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tindexing: " + ' '.join(idx_cmd) + "\n")
    subprocess.call(idx_cmd)

    if deltmp:
        sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fq1 + "\n")
        os.remove(fq1)
        sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fq2 + "\n")
        os.remove(fq2)

    return bamreadcount(outbam)


def remap_novoalign_fastq(fq1, fq2, threads, fastaref, novoref, outbam, deltmp=True, mutid='null'):
    """ call novoalign and samtools to remap .bam
    """

    basefn   = "novotmp." + str(uuid4())
    sam_out  = basefn + '.sam'
    sort_out = basefn + '.sorted.bam'

    sam_cmd  = ['novoalign', '-F', 'STDFQ', '-f', fq1, fq2, '-r', 'Random', '-d', novoref, '-oSAM']
    bam_cmd  = ['samtools', 'view', '-bt', fastaref + '.fai', '-o', outbam, sam_out]
    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-m', '10000000000', '-T', sort_out, '-o', sort_out, outbam]
    idx_cmd  = ['samtools', 'index', outbam]

    print "INFO\t" + now() + "\t" + mutid + "\taligning " + fq1 + ',' + fq2 + " with novoalign"
    with open(sam_out, 'w') as sam:
        p = subprocess.Popen(sam_cmd, stdout=subprocess.PIPE)
        for line in p.stdout:
            sam.write(line)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\twriting " + sam_out + " to BAM...\n")
    subprocess.call(bam_cmd)

    if deltmp:
        sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tdeleting SAM: " + sam_out + "\n")
        os.remove(sam_out)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tsorting output: " + ' '.join(sort_cmd) + "\n")
    subprocess.call(sort_cmd)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremove original bam:" + outbam + "\n")
    os.remove(outbam)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\trename sorted bam: " + sort_out + " to original name: " + outbam + "\n")
    move(sort_out, outbam)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tindexing: " + ' '.join(idx_cmd) + "\n")
    subprocess.call(idx_cmd)

    if deltmp:
        sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fq1 + "\n")
        os.remove(fq1)
        sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fq2 + "\n")
        os.remove(fq2)

    return bamreadcount(outbam)


def remap_backtrack_fastq(fq1, fq2, threads, fastaref, outbam, deltmp=True, mutid='null'):
    """ call bwa/samtools to remap .bam and merge with existing .bam
    """
    basefn = "bwatmp." + str(uuid4())
    sai1fn = basefn + ".1.sai"
    sai2fn = basefn + ".2.sai"
    samfn  = basefn + ".sam"
    refidx = fastaref + ".fai"
    tmpbam = basefn + ".bam"
    tmpsrt = basefn + ".sort.bam"

    sai1args = ['bwa', 'aln', '-q', '5', '-l', '32', '-k', '2', '-t', str(threads), '-o', '1', '-f', sai1fn, fastaref, fq1]
    sai2args = ['bwa', 'aln', '-q', '5', '-l', '32', '-k', '2', '-t', str(threads), '-o', '1', '-f', sai2fn, fastaref, fq2]
    samargs  = ['bwa', 'sampe', '-P', '-f', samfn, fastaref, sai1fn, sai2fn, fq1, fq2]
    bamargs  = ['samtools', 'view', '-bt', refidx, '-o', tmpbam, samfn]
    sortargs = ['samtools', 'sort', '-T', tmpsrt, '-o', tmpsrt, tmpbam]

    print "INFO\t" + now() + "\t" + mutid + "\tmapping 1st end, cmd: " + " ".join(sai2args)
    p = subprocess.Popen(sai2args, stderr=subprocess.STDOUT)
    p.wait()

    print "INFO\t" + now() + "\t" + mutid + "\tmapping 2nd end, cmd: " + " ".join(sai1args)
    p = subprocess.Popen(sai1args, stderr=subprocess.STDOUT)
    p.wait()

    print "INFO\t" + now() + "\t" + mutid + "\tpairing ends, building .sam, cmd: " + " ".join(samargs)
    p = subprocess.Popen(samargs, stderr=subprocess.STDOUT)
    p.wait()

    print "INFO\t" + now() + "\t" + mutid + "\tsam --> bam, cmd: " + " ".join(bamargs)
    p = subprocess.Popen(bamargs, stderr=subprocess.STDOUT)
    p.wait()

    print "INFO\t" + now() + "\t" + mutid + "\tsorting, cmd: " + " ".join(sortargs)
    p = subprocess.Popen(sortargs, stderr=subprocess.STDOUT)
    p.wait()

    print "INFO\t" + now() + "\t" + mutid + "\trename " + tmpsrt + ".bam --> " + tmpbam
    os.remove(tmpbam)
    os.rename(tmpsrt, tmpbam)

    print "INFO\t" + now() + "\t" + mutid + "\trename " + tmpbam + " --> " + outbam
    os.rename(tmpbam, outbam)

    # cleanup
    if deltmp:
        os.remove(sai1fn)
        os.remove(sai2fn)
        os.remove(samfn)
        os.remove(fq1)
        os.remove(fq2)

    return bamreadcount(outbam)

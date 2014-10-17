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


def remap_backtrack_bam(bamfn, threads, bwaref, mutid='null', paired=True):
    """ call bwa/samtools to remap .bam
    """
    if paired:
        sai1fn = bamfn + ".1.sai"
        sai2fn = bamfn + ".2.sai"
        samfn  = bamfn + ".sam"
        refidx = bwaref + ".fai"

        sai1args = ['bwa', 'aln', '-q', '5', '-l', '32', '-k', '3', '-t', str(threads), '-o', '1', '-f', sai1fn, '-b1', bwaref, bamfn]
        sai2args = ['bwa', 'aln', '-q', '5', '-l', '32', '-k', '3', '-t', str(threads), '-o', '1', '-f', sai2fn, '-b2', bwaref, bamfn]
        samargs  = ['bwa', 'sampe', '-P', '-f', samfn, bwaref, sai1fn, sai2fn, bamfn, bamfn]
        bamargs  = ['samtools', 'view', '-bt', refidx, '-o', bamfn, samfn] 

        print "INFO\t" + now() + "\t" + mutid + "\tmapping 1st end, cmd: " + " ".join(sai1args)
        subprocess.call(sai1args)
        print "INFO\t" + now() + "\t" + mutid + "\tmapping 2nd end, cmd: " + " ".join(sai2args)
        subprocess.call(sai2args)
        print "INFO\t" + now() + "\t" + mutid + "\tpairing ends, building .sam, cmd: " + " ".join(samargs)
        subprocess.call(samargs)
        print "INFO\t" + now() + "\t" + mutid + "\tsam --> bam, cmd: " + " ".join(bamargs)
        subprocess.call(bamargs)

        sortbase = bamfn + ".sort"
        sortfn   = sortbase + ".bam"
        sortargs = ['samtools','sort','-m','10000000000',bamfn,sortbase]
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
        refidx = bwaref + ".fai"

        saiargs = ['bwa', 'aln', '-q', '5', '-l', '32', '-k', '3', '-t', str(threads), '-o', '1', '-f', saifn, '-b1', bwaref, bamfn]
        samargs  = ['bwa', 'samse', '-f', samfn, bwaref, saifn, bamfn]
        bamargs  = ['samtools', 'view', '-bt', refidx, '-o', bamfn, samfn] 

        print "INFO\t" + now() + "\t" + mutid + "\tmapping, cmd: " + " ".join(saiargs)
        subprocess.call(saiargs)
        print "INFO\t" + now() + "\t" + mutid + "\tpairing ends, building .sam, cmd: " + " ".join(samargs)
        subprocess.call(samargs)
        print "INFO\t" + now() + "\t" + mutid + "\tsam --> bam, cmd: " + " ".join(bamargs)
        subprocess.call(bamargs)

        sortbase = bamfn + ".sort"
        sortfn   = sortbase + ".bam"
        sortargs = ['samtools','sort','-m','10000000000',bamfn,sortbase]
        print "INFO\t" + now() + "\t" + mutid + "\tsorting, cmd: " + " ".join(sortargs)
        subprocess.call(sortargs)
        os.rename(sortfn,bamfn)

        indexargs = ['samtools','index',bamfn]
        print "INFO\t" + now() + "\t" + mutid + "\tindexing, cmd: " + " ".join(indexargs)
        subprocess.call(indexargs)

        # cleanup
        os.remove(saifn)
        os.remove(samfn)


def remap_bwamem_bam(bamfn, threads, bwaref, samtofastq, mutid='null', paired=True):
    """ call bwa mem and samtools to remap .bam
    """
    assert os.path.exists(samtofastq)
    assert bamreadcount(bamfn) > 0
    if paired:
        assert bamreadcount(bamfn) > 1 

    sam_out  = bamfn + '.realign.sam'
    sort_out = bamfn + '.realign.sorted'

    print "INFO\t" + now() + "\t" + mutid + "\tconverting " + bamfn + " to fastq\n"
    fastq = bamtofastq(bamfn, samtofastq, threads=threads, paired=paired)

    sam_cmd = []

    if paired:
        sam_cmd  = ['bwa', 'mem', '-t', str(threads), '-M', '-p', bwaref, fastq] # interleaved
    else:
        sam_cmd  = ['bwa', 'mem', '-t', str(threads), '-M', bwaref, fastq] # single-end

    assert len(sam_cmd) > 0

    bam_cmd  = ['samtools', 'view', '-bt', bwaref + '.fai', '-o', bamfn, sam_out]
    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-m', '10000000000', bamfn, sort_out]
    idx_cmd  = ['samtools', 'index', bamfn]

    print "INFO\t" + now() + "\t" + mutid + "\taligning " + fastq + " with bwa mem\n"
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

    sort_out += '.bam'

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremove original bam:" + bamfn + "\n")
    os.remove(bamfn)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\trename sorted bam: " + sort_out + " to original name: " + bamfn + "\n")
    move(sort_out, bamfn)

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tindexing: " + ' '.join(idx_cmd) + "\n")
    subprocess.call(idx_cmd)

    # check if BAM readcount looks sane
    # check if BAM readcount looks sane
    if bamreadcount(bamfn) < fastqreadcount(fastq): 
        raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tbam readcount < fastq readcount, alignment sanity check failed!\n")

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fastq + "\n")
    os.remove(fastq)


def remap_novoalign_bam(bamfn, threads, bwaref, samtofastq, novoref, mutid='null', paired=True):
    """ call novoalign and samtools to remap .bam
    """
    assert os.path.exists(samtofastq)
    assert os.path.exists(novoref)
    assert bamreadcount(bamfn) > 0 

    sam_out  = bamfn + '.realign.sam'
    sort_out = bamfn + '.realign.sorted'

    print "INFO\t" + now() + "\t" + mutid + "\tconverting " + bamfn + " to fastq\n"
    fastq = bamtofastq(bamfn, samtofastq, threads=threads, paired=paired, twofastq=True)

    sam_cmd = []

    if paired:
        sam_cmd  = ['novoalign', '--mmapOff', '-F', 'STDFQ', '-f', fastq[0], fastq[1], '-r', 'Random', '-d', novoref, '-oSAM'] # interleaved
    else:
        sam_cmd  = ['novoalign', '--mmapOff', '-F', 'STDFQ', '-f', fastq, '-r', 'Random', '-d', novoref, '-oSAM'] # interleaved

    assert len(sam_cmd) > 0

    bam_cmd  = ['samtools', 'view', '-bt', bwaref + '.fai', '-o', bamfn, sam_out]
    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-m', '10000000000', bamfn, sort_out]
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

    sort_out += '.bam'

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
        if bamreadcount(bamfn) < fastqreadcount(fastq): 
            raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tbam readcount < fastq readcount, alignment sanity check failed!\n")

    if paired:
        for fq in fastq:
            sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fq + "\n")
            os.remove(fq)
    else:
        sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fastq + "\n")
        os.remove(fastq)


def remap_gsnap_bam(bamfn, threads, bwaref, samtofastq, gsnaprefdir, gsnaprefname, mutid='null', paired=True):
    """ call gsnap and samtools to remap .bam
    """
    assert os.path.exists(samtofastq)
    assert os.path.exists(gsnaprefdir)
    assert bamreadcount(bamfn) > 0 

    sam_out  = bamfn + '.realign.sam'
    sort_out = bamfn + '.realign.sorted'

    print "INFO\t" + now() + "\t" + mutid + "\tconverting " + bamfn + " to fastq\n"
    fastq = bamtofastq(bamfn, samtofastq, threads=threads, paired=paired, twofastq=True)

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

    bam_cmd  = ['samtools', 'view', '-bt', bwaref + '.fai', '-o', bamfn, sam_out]
    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-m', '10000000000', bamfn, sort_out]
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

    sort_out += '.bam'

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
        if bamreadcount(bamfn) < fastqreadcount(fastq): 
            raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tbam readcount < fastq readcount, alignment sanity check failed!\n")

    if paired:
        for fq in fastq:
            sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fq + "\n")
            os.remove(fq)
    else:
        sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fastq + "\n")
        os.remove(fastq)


#
# Remapping functions paired fastq --> BAM
#


def remap_bwamem_fastq(fq1, fq2, threads, bwaref, outbam, deltmp=True, mutid='null'):
    """ call bwa mem and samtools to remap .bam
    """

    basefn   = "bwatmp." + str(uuid4())
    sam_out  = basefn + '.sam'
    sort_out = basefn + '.sorted'

    sam_cmd  = ['bwa', 'mem', '-t', str(threads), '-M', '-Y', bwaref, fq1, fq2]

    bam_cmd  = ['samtools', 'view', '-bt', bwaref + '.fai', '-o', outbam, sam_out]
    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-m', '10000000000', outbam, sort_out]
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

    sort_out += '.bam'

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


def remap_novoalign_fastq(fq1, fq2, threads, bwaref, novoref, outbam, deltmp=True, mutid='null'):
    """ call novoalign and samtools to remap .bam
    """

    basefn   = "novotmp." + str(uuid4())
    sam_out  = basefn + '.sam'
    sort_out = basefn + '.sorted'

    sam_cmd  = ['novoalign', '-F', 'STDFQ', '-f', fq1, fq2, '-r', 'Random', '-d', novoref, '-oSAM']
    bam_cmd  = ['samtools', 'view', '-bt', bwaref + '.fai', '-o', outbam, sam_out]
    sort_cmd = ['samtools', 'sort', '-@', str(threads), '-m', '10000000000', outbam, sort_out]
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

    sort_out += '.bam'

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


def remap_backtrack_fastq(fq1, fq2, threads, bwaref, outbam, deltmp=True, mutid='null'):
    """ call bwa/samtools to remap .bam and merge with existing .bam
    """
    basefn = "bwatmp." + str(uuid4())
    sai1fn = basefn + ".1.sai"
    sai2fn = basefn + ".2.sai"
    samfn  = basefn + ".sam"
    refidx = bwaref + ".fai"
    tmpbam = basefn + ".bam"
    tmpsrt = basefn + ".sort"

    sai1args = ['bwa', 'aln', '-q', '5', '-l', '32', '-k', '2', '-t', str(threads), '-o', '1', '-f', sai1fn, bwaref, fq1]
    sai2args = ['bwa', 'aln', '-q', '5', '-l', '32', '-k', '2', '-t', str(threads), '-o', '1', '-f', sai2fn, bwaref, fq2]
    samargs  = ['bwa', 'sampe', '-P', '-f', samfn, bwaref, sai1fn, sai2fn, fq1, fq2]
    bamargs  = ['samtools', 'view', '-bt', refidx, '-o', tmpbam, samfn]
    sortargs = ['samtools', 'sort', tmpbam, tmpsrt]

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
    os.rename(tmpsrt + ".bam", tmpbam)

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
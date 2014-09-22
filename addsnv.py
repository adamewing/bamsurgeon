#!/usr/bin/env python

import sys
import pysam
import argparse
import random
import subprocess
import os
import bs.replacereads as rr
import datetime
import traceback

from uuid import uuid4
from shutil import move
from re import sub
from multiprocessing import Pool, Value, Lock
from collections import Counter


sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 0)

def now():
    return str(datetime.datetime.now())

def majorbase(basepile):
    """returns tuple: (major base, count)
    """
    return Counter(basepile).most_common()[0]

def minorbase(basepile):
    """returns tuple: (minor base, count)
    """
    c = Counter(basepile)
    if len(list(c.elements())) > 1:
        return c.most_common(2)[-1]
    else:
        return c.most_common()[0]

def mut(base, altbase):
    """ change base to something different
    """

    bases = ('A','T','C','G')
    base = base.upper()
    if base not in bases or (altbase is not None and altbase not in ['A','T','G','C']):
        raise ValueError("ERROR\t" + now() + "\tbase passed to mut(): " + str(base) + " not one of (A,T,C,G)\n")

    if altbase is not None:
        return altbase

    else:
        alt = base
        while alt == base:
            alt = bases[int(random.uniform(0,4))]
        return alt

def countReadCoverage(bam,chrom,start,end):
    """ calculate coverage of aligned reads over region
    """

    coverage = []
    start = int(start)
    end = int(end)
    for i in range(end-start+1):
        coverage.append(0.0)

    i = 0
    if chrom in bam.references:
        for pcol in bam.pileup(chrom,start,end):
            n = 0
            if pcol.pos >= start and pcol.pos <= end:
                for read in pcol.pileups:
                    if read.alignment.mapq >= 0 and not read.alignment.is_duplicate:
                        n += 1
                coverage[i] = n
                i += 1

    return coverage


def countBaseAtPos(bamfile, chrom, pos, mutid='null'):
    """ return list of bases at position chrom,pos
    """
    locstr = chrom + ":" + str(pos) + "-" + str(pos)
    args = ['samtools','mpileup',bamfile,'-r',locstr]

    p = subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    p.wait()
    pout = p.stdout.readlines()

    pileup = None 

    for line in pout:
        try:
            c = line.strip().split()
            assert len(c) > 5
            pileup = c[4].upper()
        except AssertionError:
            sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\tmpileup failed, no coverage for base: " + chrom + ":" + str(pos))
            return []
    bases = []
    if pileup:
        for b in pileup:
            if b in ['A','T','C','G']:
                bases.append(b)

    return bases

def mergebams(bamlist, outbamfn, maxopen=1000):
    """ call samtools to merge two .bams
    """

    assert outbamfn.endswith('.bam')
    print "INFO\t" + now() + "\tlen(bamlist)", len(bamlist)

    if len(bamlist) == 1:
        print "INFO\t" + now() + "\tonly one BAM to merge, renaming",bamlist[0],"-->",outbamfn
        move(bamlist[0], outbamfn)
    else:
        nmerge = 1
        mergenum = 0
        merge_sublists = {}
        for tmpbam in bamlist:
            mergetmp = "tmp.merging." + str(mergenum) + "." + outbamfn
            if mergetmp in merge_sublists:
                merge_sublists[mergetmp].append(tmpbam)
            else:
                merge_sublists[mergetmp] = []
                merge_sublists[mergetmp].append(tmpbam)
            if nmerge % maxopen == 0:
                mergenum += 1
            nmerge += 1

        for submergefn, tmpbams in merge_sublists.iteritems():
            if len(tmpbams) == 1:
                move(tmpbams[0], submergefn)
                print "INFO\t" + now() + "\trenamed:",tmpbams[0],"-->",submergefn
            else:
                args = ['samtools','merge','-f',submergefn] + tmpbams 
                print "INFO\t" + now() + "\tmerging, cmd: ",args
                subprocess.call(args)

        if len(merge_sublists.keys()) == 1:
            print "INFO\t" + now() + "\tmerge finished, renaming:",merge_sublists.keys()[0],"-->",outbamfn
            move(merge_sublists.keys()[0], outbamfn)
        else:
            args = ['samtools','merge','-f',outbamfn] + merge_sublists.keys() 
            print "INFO\t" + now() + "\tfinal merge, cmd: ",args
            subprocess.call(args)

        for submergefn in merge_sublists.keys():
            if os.path.exists(submergefn):
                os.remove(submergefn)

    for bamfile in bamlist:
        if os.path.exists(bamfile):
            os.remove(bamfile)
            os.remove(bamfile + '.bai')

def remap_paired(bamfn, threads, bwaref, mutid='null'):
    """ call bwa/samtools to remap .bam
    """
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
    move(sortfn,bamfn)

    indexargs = ['samtools','index',bamfn]
    print "INFO\t" + now() + "\t" + mutid + "\tindexing, cmd: " + " ".join(indexargs)
    subprocess.call(indexargs)

    # cleanup
    os.remove(sai1fn)
    os.remove(sai2fn)
    os.remove(samfn)

############

def remap_bwamem(bamfn, threads, bwaref, samtofastq, mutid='null', paired=True):
    """ call bwa mem and samtools to remap .bam
    """
    assert os.path.exists(samtofastq)
    assert bamreadcount(bamfn) > 0 

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
    if bamreadcount(bamfn) < fastqreadcount(fastq): 
        raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tbam readcount < fastq readcount, alignment sanity check failed!\n")

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tremoving " + fastq + "\n")
    os.remove(fastq)


def remap_novoalign(bamfn, threads, bwaref, samtofastq, novoref, mutid='null', paired=True):
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
        sam_cmd  = ['novoalign', '-F', 'STDFQ', '-f', fastq[0], fastq[1], '-r', 'Random', '-d', novoref, '-oSAM'] # interleaved
    else:
        sam_cmd  = ['novoalign', '-F', 'STDFQ', '-f', fastq, '-r', 'Random', '-d', novoref, '-oSAM'] # interleaved

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


def remap_gsnap(bamfn, threads, bwaref, samtofastq, gsnaprefdir, gsnaprefname, mutid='null', paired=True):
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


def bamtofastq(bam, samtofastq, threads=1, paired=True, twofastq=False):
    ''' if twofastq is True output two fastq files instead of interleaved (default) for paired-end'''
    assert os.path.exists(samtofastq)
    assert bam.endswith('.bam')

    outfq = None
    outfq_pair = None

    cmd = ['java', '-XX:ParallelGCThreads=' + str(threads), '-Xmx4g', '-jar', samtofastq, 'VALIDATION_STRINGENCY=SILENT', 'INPUT=' + bam]
    if paired:
        if twofastq: # two-fastq paired end
            outfq_pair = [sub('bam$', '1.fastq', bam), sub('bam$', '2.fastq', bam)]
            cmd.append('F=' + outfq_pair[0])
            cmd.append('F2=' + outfq_pair[1])
        else: # interleaved paired-end
            outfq = sub('bam$', 'fastq', bam)
            cmd.append('FASTQ=' + outfq)
            cmd.append('INTERLEAVE=true')
    else:
        outfq = sub('bam$', 'fastq', bam)
        cmd.append('FASTQ=' + outfq)

    sys.stdout.write("INFO\t" + now() + "\tconverting BAM " + bam + " to FASTQ\n")
    subprocess.call(cmd)

    if outfq is not None:
        assert os.path.exists(outfq) # conversion failed
        return outfq

    if outfq_pair is not None:
        assert os.path.exists(outfq_pair[0]) and os.path.exists(outfq_pair[1])
        return outfq_pair

    return None


def bamreadcount(bamfile):
    bam = pysam.Samfile(bamfile, 'rb')
    if os.path.exists(bamfile + '.bai'):
        return bam.mapped + bam.unmapped
    else:
        return(list(bam.fetch(until_eof=True)))


def fastqreadcount(fastqfile):
    assert not fastqfile.endswith('gz') # not supported yet
    return sum(1 for line in open(fastqfile))/4


#########


def remap_single(bamfn, threads, bwaref, mutid='null'):
    """ call bwa/samtools to remap .bam
    """
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
    move(sortfn,bamfn)

    indexargs = ['samtools','index',bamfn]
    print "INFO\t" + now() + "\t" + mutid + "\tindexing, cmd: " + " ".join(indexargs)
    subprocess.call(indexargs)

    # cleanup
    os.remove(saifn)
    os.remove(samfn)


def replace(origbamfile, mutbamfile, outbamfile):
    ''' open .bam file and call replacereads
    '''
    origbam = pysam.Samfile(origbamfile, 'rb')
    mutbam  = pysam.Samfile(mutbamfile, 'rb')
    outbam  = pysam.Samfile(outbamfile, 'wb', template=origbam)

    rr.replaceReads(origbam, mutbam, outbam, keepqual=True)

    origbam.close()
    mutbam.close()
    outbam.close()


def dictlist(fn):
    d = {}
    with open(fn, 'r') as inlist:
        for name in inlist:
            d[name.strip()] = True
    return d


def makemut(args, chrom, start, end, vaf, altbase, avoid):
    mutid = chrom + '_' + str(start) + '_' + str(end) + '_' + str(vaf) + '_' + str(altbase)
    try:
        bamfile = pysam.Samfile(args.bamFileName, 'rb')
        bammate = pysam.Samfile(args.bamFileName, 'rb') # use for mates to avoid iterator problems
        reffile = pysam.Fastafile(args.refFasta)
        tmpbams = []

        snvfrac = float(args.snvfrac)

        mutpos = int(random.uniform(start,end+1)) # position of mutation in genome
        refbase = reffile.fetch(chrom,mutpos-1,mutpos)

        if altbase == refbase.upper():
            sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\tspecified ALT base matches reference, skipping mutation\n")
            return None

        try:
            mutbase = mut(refbase, altbase)
        except ValueError as e:
            sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\t" + ' '.join(("skipped site:",chrom,str(start),str(end),"due to N base:",str(e),"\n")))
            return None

        mutstr = refbase + "-->" + mutbase

        # optional CNV file
        cnv = None
        if (args.cnvfile):
            cnv = pysam.Tabixfile(args.cnvfile, 'r')

        log = open('addsnv_logs_' + os.path.basename(args.outBamFile) + '/' + os.path.basename(args.outBamFile) + "." + "_".join((chrom,str(start),str(end))) + ".log",'w')

        # keep a list of reads to modify - use hash to keep unique since each
        # read will be visited as many times as it has bases covering the region
        outreads = {}
        mutreads = {} # same keys as outreads
        mutmates = {} # same keys as outreads, keep track of mates
        numunmap = 0
        hasSNP = False
        tmpoutbamname = args.tmpdir + "/" + mutid + ".tmpbam." + str(uuid4()) + ".bam"
        print "INFO\t" + now() + "\t" + mutid + "\tcreating tmp bam: ",tmpoutbamname
        outbam_muts = pysam.Samfile(tmpoutbamname, 'wb', template=bamfile)
        maxfrac = 0.0

        for pcol in bamfile.pileup(reference=chrom,start=mutpos,end=mutpos+1):
            # this will include all positions covered by a read that covers the region of interest
            if pcol.pos: #> start and pcol.pos <= end:
                refbase = reffile.fetch(chrom,pcol.pos-1,pcol.pos)
                basepile = ''
                for pread in pcol.pileups:
                    if avoid is not None and pread.alignment.qname in avoid:
                        print "WARN\t" + now() + "\t" + mutid + "\tdropped mutation due to read in --avoidlist", pread.alignment.qname
                        os.remove(tmpoutbamname)
                        return None

                    if not pread.alignment.is_secondary: # only consider primary alignments
                        basepile += pread.alignment.seq[pread.qpos-1]
                        pairname = 'F' # read is first in pair
                        if pread.alignment.is_read2:
                            pairname = 'S' # read is second in pair
                        if not pread.alignment.is_paired:
                            pairname = 'U' # read is unpaired

                        extqname = ','.join((pread.alignment.qname,str(pread.alignment.pos),pairname))

                        if pcol.pos == mutpos:
                            if not pread.alignment.is_secondary and not pread.alignment.mate_is_unmapped:
                                outreads[extqname] = pread.alignment
                                mutbases = list(pread.alignment.seq)
                                mutbases[pread.qpos-1] = mutbase
                                mutread = ''.join(mutbases)
                                mutreads[extqname] = mutread
                                mate = None
                                if not args.single:
                                    try:
                                        mate = bammate.mate(pread.alignment)
                                    except:
                                        print "WARN\t" + now() + "\t" + mutid + "\twarning: no mate for", pread.alignment.qname
                                        if args.requirepaired:
                                            print "WARN\t" + now() + "\t" + mutid + "\tskipped mutation due to --requirepaired"
                                            os.remove(tmpoutbamname)
                                            return None
                                mutmates[extqname] = mate
                                log.write(" ".join(('read',extqname,mutread,"\n")))
                            else:
                                numunmap += 1

                            if len(mutreads) > 200: # FIXME maxdepth should be a parameter
                                sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\tdepth at site is greater than cutoff, aborting mutation.\n")
                                outbam_muts.close()
                                os.remove(tmpoutbamname)
                                return None

                # make sure region doesn't have any changes that are likely SNPs
                # (trying to avoid messing with haplotypes)
                    
                basepile = countBaseAtPos(args.bamFileName,chrom,pcol.pos,mutid=mutid)
                if basepile:
                    majb = majorbase(basepile)
                    minb = minorbase(basepile)

                    frac = float(minb[1])/(float(majb[1])+float(minb[1]))
                    if minb[0] == majb[0]:
                        frac = 0.0
                    if frac > maxfrac:
                        maxfrac = frac
                    if (not args.ignoresnps) and (frac > snvfrac or majb[0].upper() != refbase.upper()):
                        print "WARN\t" + now() + "\t" + mutid + "\tdropped for proximity to SNP, nearby SNP MAF:",frac,"maxfrac:",snvfrac
                        hasSNP = True
                else:
                    print "WARN\t" + now() + "\t" + mutid + "\tcould not pileup for region:",chrom,pcol.pos
                    hasSNP = True

        # pick reads to change
        readlist = []
        for extqname,read in outreads.iteritems():
            if read.seq != mutreads[extqname]:
                readlist.append(extqname)

        print "INFO\t" + now() + "\t" + mutid + "\tlen(readlist): " + str(len(readlist))
        random.shuffle(readlist)

        if len(readlist) < int(args.mindepth):
            print "WARN\t" + now() + "\t" + mutid + "\ttoo few reads in region (" + str(len(readlist)) + ") skipping..."
            outbam_muts.close()
            os.remove(tmpoutbamname)
            return None

        if vaf is None:
            vaf = float(args.mutfrac) # default minor allele freq if not otherwise specified
        if cnv: # cnv file is present
            if chrom in cnv.contigs:
                for cnregion in cnv.fetch(chrom,start,end):
                    cn = float(cnregion.strip().split()[3]) # expect chrom,start,end,CN
                    print "INFO\t" + now() + "\t" + mutid + "\t" + ' '.join(("copy number in snp region:",chrom,str(start),str(end),"=",str(cn))) + "\n"
                    if float(cn) > 0.0:
                        vaf = 1.0/float(cn)
                    else:
                        vaf = 0.0
                    print "adjusted VAF: " + str(vaf) + "\n"
        else:
            print "INFO\t" + now() + "\t" + mutid + "\tselected VAF: " + str(vaf) + "\n"

        ######
        lastread = int(len(readlist)*vaf)

        # pick at least args.minmutreads if possible
        if lastread < int(args.minmutreads):
            if len(readlist) > int(args.minmutreads):
                lastread = int(args.minmutreads)
                sys.stdout.write("WARN\t" + now() + "\t" + mutid + "\tforced " + str(lastread) + " reads.\n")
            else:
                print "WARN\t" + now() + "\t" + mutid + "\tdropped site with fewer reads than --minmutreads"
                os.remove(tmpoutbamname)
                return None

        readlist = readlist[0:lastread] 

        ######

        print "INFO\t" + now() + "\t" + mutid + "\tpicked:",str(len(readlist))

        wrote = 0
        nmut = 0
        # change reads from .bam to mutated sequences
        for extqname,read in outreads.iteritems():
            if read.seq != mutreads[extqname]:
                if not args.nomut and extqname in readlist:
                    qual = read.qual # changing seq resets qual (see pysam API docs)
                    read.seq = mutreads[extqname] # make mutation
                    read.qual = qual
                    nmut += 1
            if not hasSNP or args.force:
                wrote += 1
                outbam_muts.write(read)
                if args.single:
                    outbam_muts.write(mutmates[extqname])
                elif mutmates[extqname] is not None:
                    outbam_muts.write(mutmates[extqname])
        print "INFO\t" + now() + "\t" + mutid + "\twrote: ",wrote,"mutated:",nmut

        if not hasSNP or args.force:
            outbam_muts.close()
            if args.single:
                if args.bwamem:
                    remap_bwamem(tmpoutbamname, 1, args.refFasta, args.samtofastq, mutid=mutid, paired=False)
                elif args.novoalign:
                    remap_novoalign(tmpoutbamname, 1, args.refFasta, args.samtofastq, args.novoref, mutid=mutid, paired=False)
                elif args.gsnap:
                    remap_gsnap(tmpoutbamname, 1, args.refFasta, args.samtofastq, args.gsnaprefdir, args.gsnaprefname, mutid=mutid, paired=False)
                else:
                    remap_single(tmpoutbamname, 1, args.refFasta, mutid=mutid)
            else:
                if args.bwamem:
                    remap_bwamem(tmpoutbamname, 1, args.refFasta, args.samtofastq, mutid=mutid, paired=True)
                elif args.novoalign:
                    remap_novoalign(tmpoutbamname, 1, args.refFasta, args.samtofastq, args.novoref, mutid=mutid, paired=True)
                elif args.gsnap:
                    remap_gsnap(tmpoutbamname, 1, args.refFasta, args.samtofastq, args.gsnaprefdir, args.gsnaprefname, mutid=mutid, paired=True)
                else:
                    remap_paired(tmpoutbamname, 1, args.refFasta, mutid=mutid)

            outbam_muts = pysam.Samfile(tmpoutbamname,'rb')
            coverwindow = 1
            incover  = countReadCoverage(bamfile,chrom,mutpos-coverwindow,mutpos+coverwindow)
            outcover = countReadCoverage(outbam_muts,chrom,mutpos-coverwindow,mutpos+coverwindow)

            avgincover  = float(sum(incover))/float(len(incover)) 
            avgoutcover = float(sum(outcover))/float(len(outcover))

            print "INFO\t" + now() + "\t" + mutid + "\tavgincover: " + str(avgincover) + " avgoutcover: " + str(avgoutcover)

            spikein_snvfrac = 0.0
            if wrote > 0:
                spikein_snvfrac = float(nmut)/float(wrote)

            # qc cutoff for final snv depth 
            if (avgoutcover > 0 and avgincover > 0 and avgoutcover/avgincover >= float(args.coverdiff)) or args.force:
                tmpbams.append(tmpoutbamname)
                snvstr = chrom + ":" + str(start) + "-" + str(end) + " (VAF=" + str(vaf) + ")"
                log.write("\t".join(("snv",snvstr,str(mutpos),mutstr,str(avgoutcover),str(avgoutcover),str(spikein_snvfrac),str(maxfrac)))+"\n")
            else:
                outbam_muts.close()
                os.remove(tmpoutbamname)
                if os.path.exists(tmpoutbamname + '.bai'):
                    os.remove(tmpoutbamname + '.bai')

                return None

        outbam_muts.close()
        bamfile.close()
        bammate.close()
        log.close() 

        return tmpbams

    except Exception, e:
        sys.stderr.write("*"*60 + "\nERROR\t" + now() + "\tencountered error in mutation spikein: " + mutid + "\n")
        traceback.print_exc(file=sys.stdout)
        sys.stderr.write("*"*60 + "\n")
        if os.path.exists(tmpoutbamname):
            os.remove(tmpoutbamname)
        if os.path.exists(tmpoutbamname + '.bai'):
            os.remove(tmpoutbamname + '.bai')
        return None


def main(args):
    print "INFO\t" + now() + "\tstarting " + sys.argv[0] + " called with args: " + ' '.join(sys.argv) + "\n"
    bedfile = open(args.varFileName, 'r')
    reffile = pysam.Fastafile(args.refFasta)

    if (args.bwamem or args.novoalign or args.gsnap) and args.samtofastq is None:
        sys.stderr.write("ERROR\t" + now() + "\t --samtofastq must be specified with --bwamem or --novoalign option\n")
        sys.exit(1)

    if sum((args.bwamem, args.novoalign, args.gsnap)) > 1:
        sys.stderr.write("ERROR\t" + now() + "\t only one aligner allowed e.g. --bwamem and --novoalign cannot be specified together\n")
        sys.exit(1)

    if args.novoalign and args.novoref is None:
        sys.stderr.write("ERROR\t" + now() + "\t --novoref must be be specified with --novoalign\n")
        sys.exit(1)

    if args.gsnap and (args.gsnaprefname is None or args.gsnaprefdir is None):
        sys.stderr.write("ERROR\t" + now() + "\t --gsnaprefname and --gsnaprefdir must be be specified with --gsnap\n")
        sys.exit(1)

    # load readlist to avoid, if specified
    avoid = None
    if args.avoidreads is not None:
        avoid = dictlist(args.avoidreads)

    # make a temporary file to hold mutated reads
    outbam_mutsfile = "addsnv." + str(uuid4()) + ".muts.bam"
    bamfile = pysam.Samfile(args.bamFileName, 'rb')
    outbam_muts = pysam.Samfile(outbam_mutsfile, 'wb', template=bamfile)
    outbam_muts.close()
    bamfile.close()
    tmpbams = []

    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)
        print "INFO\t" + now() + "\tcreated tmp directory: " + args.tmpdir

    if not os.path.exists('addsnv_logs_' + os.path.basename(args.outBamFile)):
        os.mkdir('addsnv_logs_' + os.path.basename(args.outBamFile))
        print "INFO\t" + now() + "\tcreated directory: addsnv_logs_" + os.path.basename(args.outBamFile)

    assert os.path.exists('addsnv_logs_' + os.path.basename(args.outBamFile)), "could not create output directory!"
    assert os.path.exists(args.tmpdir), "could not create temporary directory!"

    pool = Pool(processes=int(args.procs))
    results = []

    ntried = 0
    for bedline in bedfile:
        if ntried < int(args.numsnvs) or int(args.numsnvs) == 0:
            c = bedline.strip().split()
            chrom   = c[0]
            start   = int(c[1])
            end     = int(c[2])
            vaf     = None
            altbase = None

            # VAF is 4th column, if present
            if len(c) > 3:
                vaf = float(c[3])

            # ALT is 5th column, if present
            if len(c) == 5:
                altbase = c[4].upper()
                assert altbase in ['A','T','C','G'], "ERROR:\t" + now() + "\tALT " + altbase + " not A, T, C, or G!\n"

            # make mutation (submit job to thread pool)
            result = pool.apply_async(makemut, [args, chrom, start, end, vaf, altbase, avoid])
            results.append(result)
            ntried += 1

    for result in results:
        tmpbamlist = result.get()
        if tmpbamlist is not None:
            for tmpbam in tmpbamlist:
                tmpbams.append(tmpbam)

    # merge tmp bams
    if len(tmpbams) == 1:
        move(tmpbams[0],outbam_mutsfile)
    elif len(tmpbams) > 1:
        mergebams(tmpbams,outbam_mutsfile,maxopen=int(args.maxopen))

    bedfile.close()

    # cleanup
    for bam in tmpbams:
        if os.path.exists(bam):
            os.remove(bam)
        if os.path.exists(bam + '.bai'):
            os.remove(bam + '.bai')

    if args.skipmerge:
        print "INFO\t" + now() + "\tskipping merge, plase merge reads from", outbam_mutsfile, "manually."
    else:
        print "INFO\t" + now() + "\tdone making mutations, merging mutations into", args.bamFileName, "-->", args.outBamFile
        replace(args.bamFileName, outbam_mutsfile, args.outBamFile)

        #cleanup
        os.remove(outbam_mutsfile)


def run():
    parser = argparse.ArgumentParser(description='adds SNVs to reads, outputs modified reads as .bam along with mates')
    parser.add_argument('-v', '--varfile', dest='varFileName', required=True, help='Target regions to try and add a SNV, as BED')
    parser.add_argument('-f', '--bamfile', dest='bamFileName', required=True, help='sam/bam file from which to obtain reads')
    parser.add_argument('-r', '--reference', dest='refFasta', required=True, help='reference genome, fasta indexed with bwa index -a stdsw _and_ samtools faidx')
    parser.add_argument('-o', '--outbam', dest='outBamFile', required=True, help='.bam file name for output')
    parser.add_argument('-s', '--snvfrac', dest='snvfrac', default=1, help='maximum allowable linked SNP MAF (for avoiding haplotypes) (default = 1)')
    parser.add_argument('-m', '--mutfrac', dest='mutfrac', default=0.5, help='allelic fraction at which to make SNVs (default = 0.5)')
    parser.add_argument('-n', '--numsnvs', dest='numsnvs', default=0, help="maximum number of mutations to try (default: entire input)")
    parser.add_argument('-c', '--cnvfile', dest='cnvfile', default=None, help="tabix-indexed list of genome-wide absolute copy number values (e.g. 2 alleles = no change)")
    parser.add_argument('-d', '--coverdiff', dest='coverdiff', default=0.9, help="allow difference in input and output coverage (default=0.9)")
    parser.add_argument('-p', '--procs', dest='procs', default=1, help="split into multiple processes (default=1)")
    parser.add_argument('--samtofastq', default=None, help='path to picard SamToFastq.jar')
    parser.add_argument('--mindepth', default=10, help='minimum read depth to make mutation')
    parser.add_argument('--minmutreads', default=3, help='minimum number of mutated reads to output per site')
    parser.add_argument('--avoidreads', default=None, help='file of read names to avoid (mutations will be skipped if overlap)')
    parser.add_argument('--nomut', action='store_true', default=False, help="dry run")
    parser.add_argument('--ignoresnps', action='store_true', default=False, help="make mutations even if there are non-reference alleles sharing the relevant reads")
    parser.add_argument('--force', action='store_true', default=False, help="force mutation to happen regardless of nearby SNP or low coverage")
    parser.add_argument('--single', action='store_true', default=False, help="input BAM is simgle-ended (default is paired-end)")
    parser.add_argument('--maxopen', dest='maxopen', default=1000, help="maximum number of open files during merge (default 1000)")
    parser.add_argument('--requirepaired', action='store_true', default=False, help='skip mutations if unpaired reads are present')
    parser.add_argument('--skipmerge', action='store_true', default=False, help="final output is tmp file to be merged")
    parser.add_argument('--bwamem', action='store_true', default=False, help='realignment with BWA MEM (instead of backtrack)')
    parser.add_argument('--novoalign', action='store_true', default=False, help='realignment with novoalign')
    parser.add_argument('--novoref', default=None, help='novoalign reference, must be specified with --novoalign')
    parser.add_argument('--gsnap', action='store_true', default=False, help='realignment with gsnap')
    parser.add_argument('--gsnaprefname', default=None, help='gsnap reference directory')
    parser.add_argument('--gsnaprefdir', default=None, help='gsnap reference name')
    parser.add_argument('--tmpdir', default='addsnv.tmp', help='temporary directory (default=addsnv.tmp)')
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
    run()

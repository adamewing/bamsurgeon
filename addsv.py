#!/usr/bin/env python

import re, os, sys, random
import subprocess
import traceback
import argparse
import pysam
import bs.replacereads as rr
import bs.asmregion as ar
import bs.mutableseq as ms
import datetime

from uuid import uuid4
from time import sleep
from shutil import move
from math import sqrt
from itertools import izip
from collections import Counter
from multiprocessing import Pool

sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 0)


def now():
    return str(datetime.datetime.now())


def remap_bwamem(fq1, fq2, threads, bwaref, outbam, deltmp=True, mutid='null'):
    """ call bwa mem and samtools to remap .bam
    """

    basefn   = "bwatmp." + str(uuid4())
    sam_out  = basefn + '.sam'
    sort_out = basefn + '.sorted'

    sam_cmd  = ['bwa', 'mem', '-t', str(threads), '-M', bwaref, fq1, fq2]

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

    return bamreads(outbam)


def remap_novoalign(fq1, fq2, threads, bwaref, novoref, outbam, deltmp=True, mutid='null'):
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

    return bamreads(outbam)


def remap(fq1, fq2, threads, bwaref, outbam, deltmp=True, mutid='null'):
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

    return bamreads(outbam)


def bamreads(bamfn):
    assert os.path.exists(bamfn) and bamfn.endswith('.bam')
    if not os.path.exists(bamfn + '.bai'):
        args = ['samtools', 'index', bamfn]
        subprocess.call(args)

    bam = pysam.Samfile(bamfn, 'rb')
    return bam.mapped


def runwgsim(contig, newseq, svfrac, exclude, pemean, pesd, mutid='null'):
    ''' wrapper function for wgsim
    '''
    namecount = Counter(contig.reads.reads)

    basefn = "wgsimtmp." + str(uuid4())
    fasta = basefn + ".fasta"
    fq1 = basefn + ".1.fq"
    fq2 = basefn + ".2.fq"

    fout = open(fasta,'w')
    fout.write(">target\n" + newseq + "\n")
    fout.close()

    totalreads = len(contig.reads.reads)
    paired = 0
    single = 0
    discard = 0
    pairednames = []
    # names with count 2 had both pairs in the contig
    for name,count in namecount.items():
        #print name,count
        if count == 1:
            single += 1
        elif count == 2:
            paired += 1 
            pairednames.append(name) 
        else:
            discard += 1

    print "INFO\t" + now() + "\t" + mutid + "\tpaired  reads :", paired
    print "INFO\t" + now() + "\t" + mutid + "\tsingle  reads :", single
    print "INFO\t" + now() + "\t" + mutid + "\tdiscard reads :", discard
    print "INFO\t" + now() + "\t" + mutid + "\ttotal   reads :", totalreads

    # adjustment factor for length of new contig vs. old contig
    lenfrac = float(len(newseq))/float(len(contig.seq))

    print "INFO\t" + now() + "\t" + mutid + "\told ctg len:", len(contig.seq)
    print "INFO\t" + now() + "\t" + mutid + "\tnew ctg len:", len(newseq)
    print "INFO\t" + now() + "\t" + mutid + "\tadj. factor:", lenfrac

    # number of paried reads to simulate
    nsimreads = int((paired + (single/2)) * svfrac * lenfrac)

    print "INFO\t" + now() + "\t" + mutid + "\tnum. sim. reads:", nsimreads 
    print "INFO\t" + now() + "\t" + mutid + "\tPE mean outer distance:", pemean
    print "INFO\t" + now() + "\t" + mutid + "\tPE outer distance SD:", pesd

    # length of quality score comes from original read, used here to set length of read
    maxqlen = 0
    for qual in (contig.rquals + contig.mquals):
        if len(qual) > maxqlen:
            maxqlen = len(qual)

    args = ['wgsim','-e','0','-d',str(pemean),'-s',str(pesd),'-N',str(nsimreads),'-1',str(maxqlen),'-2','101','-r','0','-R','0',fasta,fq1,fq2]
    print args
    subprocess.call(args)

    os.remove(fasta)

    fqReplaceList(fq1,pairednames,contig.rquals,svfrac,exclude, mutid=mutid)
    fqReplaceList(fq2,pairednames,contig.mquals,svfrac,exclude, mutid=mutid)

    return (fq1,fq2)


def fqReplaceList(fqfile,names,quals,svfrac,exclude,mutid='null'):
    '''
    Replace seq names in paired fastq files from a list until the list runs out
    (then stick with original names). fqfile = fastq file, names = list

    'exclude' is a filehandle, the exclude file contains read names that should
    not appear in the final output BAM

    '''
    fqin = open(fqfile,'r')

    ln = 0
    namenum = 0
    newnames = []
    seqs = []
    usednames = {}

    for fqline in fqin:
        if ln == 0:
            if len(names) > namenum:
                newnames.append(names[namenum])
            else:
                simname = fqline.strip().lstrip('@')
                simname = re.sub('/1$','',simname)  #wgsim
                simname = re.sub('/2$','',simname)  #wgsim
                newnames.append(simname) 
            namenum += 1
            ln += 1
        elif ln == 1:
            seqs.append(fqline.strip())
            ln += 1
        elif ln == 2:
            ln += 1
        elif ln == 3: # quals are passed as a list
            ln = 0
        else:
            raise ValueError("ERROR\t" + now() + "\t" + mutid + "\tfastq iteration problem\n")

    fqin.close()
    os.remove(fqfile)

    # make sure there's enough (bogus) quality scores
    while len(seqs) > len(quals):
        i = random.randint(0,len(quals)-1)
        quals.append(quals[i])

    # write .fq with new names
    fqout = open(fqfile,'w')
    for i in range(namenum):
        fqout.write("@" + newnames[i] + "\n")

        # make sure quality strings are the same length as the sequences
        while len(seqs[i]) > len(quals[i]):
            quals[i] = quals[i] + 'B'

        if len(seqs[i]) < len(quals[i]):
            quals[i] = quals[i][:len(seqs[i])]

        fqout.write(seqs[i] + "\n+\n" + quals[i] + "\n")
        if newnames[i] in usednames:
            print "INFO\t" + now() + "\t" + mutid + "\twarning, used read name: " + newnames[i] + " in multiple pairs"
        usednames[newnames[i]] = True
        
    # burn off excess
    if len(seqs) > 0:
        for name in names:
            if name not in usednames:
                if random.uniform(0,1) < svfrac:  # this controls deletion depth
                    exclude.write(name + "\n")

    fqout.close()


def singleseqfa(file,mutid='null'):
    with open(file, 'r') as fasta:
        header = None
        seq = ''
        for line in fasta:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\tmultiple entries found in " + file + " only using the first\n")
                header = line.lstrip('>')
            else:
                seq += line
    return seq


def pickseq(inslib, mutid='null'):
    pick = random.choice(inslib.keys())
    print "INFO\t" + now() + "\t" + mutid + "\tchose sequence from insertion library: " + pick
    return inslib[pick]


def load_inslib(infa):
    seqdict = {}

    with open(infa, 'r') as fa:
        seqid = ''
        seq   = ''
        for line in fa:
            if line.startswith('>'):
                if seq != '':
                    seqdict[seqid] = seq
                seqid = line.lstrip('>').strip()
                seq   = ''
            else:
                assert seqid != ''
                seq = seq + line.strip()

    if seqid not in seqdict and seq != '':
        seqdict[seqid] = seq

    return seqdict


def mergebams(bamlist, outbamfn, maxopen=100):
    ''' call samtools to merge a list of bams hierarchically '''

    assert outbamfn.endswith('.bam')
    print "INFO\t" + now() + "\tlen(bamlist)", len(bamlist)

    if len(bamlist) == 1:
        print "INFO\t" + now() + "\tonly one BAM to merge, renaming",bamlist[0],"-->",outbamfn
        os.rename(bamlist[0], outbamfn)
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
                os.rename(tmpbams[0], submergefn)
                print "INFO\t" + now() + "\trenamed:",tmpbams[0], " --> ", submergefn
            else:
                args = ['samtools','merge','-f',submergefn] + tmpbams 
                print "INFO\t" + now() + "\tmerging, cmd: ", args
                subprocess.call(args)

        if len(merge_sublists.keys()) == 1:
            print "INFO\t" + now() + "\tmerge finished, renaming: ", merge_sublists.keys()[0]," --> ", outbamfn
            os.rename(merge_sublists.keys()[0], outbamfn)
        else:
            args = ['samtools','merge','-f',outbamfn] + merge_sublists.keys()
            print "INFO\t" + now() + "\tfinal merge, cmd: ", args
            subprocess.call(args)

        for submergefn in merge_sublists.keys():
            if os.path.exists(submergefn):
                os.remove(submergefn)

    for bamfile in bamlist:
        if os.path.exists(bamfile):
            os.remove(bamfile)
            os.remove(bamfile + '.bai')


def align(qryseq, refseq):
    rnd = str(uuid4())
    tgtfa = 'tmp.' + rnd + '.tgt.fa'
    qryfa = 'tmp.' + rnd + '.qry.fa'

    tgt = open(tgtfa, 'w')
    qry = open(qryfa, 'w')

    tgt.write('>ref' + '\n' + refseq + '\n')
    qry.write('>qry' + '\n' + qryseq + '\n')

    tgt.close()
    qry.close()

    cmd = ['exonerate', '--bestn', '1', '-m', 'ungapped', '--showalignment','0', '--ryo', 'SUMMARY\t%s\t%qab\t%qae\t%tab\t%tae\n', qryfa, tgtfa]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    best = []
    topscore = 0

    for pline in p.stdout.readlines():
        if pline.startswith('SUMMARY'):
            c = pline.strip().split()
            if int(c[1]) > topscore:
                topscore = int(c[1])
                best = c

    os.remove(tgtfa)
    os.remove(qryfa)

    return best


def replace(origbamfile, mutbamfile, outbamfile, excludefile):
    ''' open .bam file and call replacereads
    '''
    origbam = pysam.Samfile(origbamfile, 'rb')
    mutbam  = pysam.Samfile(mutbamfile, 'rb')
    outbam  = pysam.Samfile(outbamfile, 'wb', template=origbam)

    rr.replaceReads(origbam, mutbam, outbam, excludefile=excludefile, allreads=True)

    origbam.close()
    mutbam.close()
    outbam.close()


def discordant_fraction(bamfile, chrom, start, end):
    r = 0
    d = 0
    bam = pysam.Samfile(bamfile, 'rb')
    for read in bam.fetch(chrom, start, end):
        r += 1
        if not read.is_proper_pair:
            d += 1

    if r > 0:
        return float(d)/float(r)
    else:
        return 0.0


def makemut(args, bedline):
    mutid = ':'.join(map(str, bedline.strip().split()))
    try:
        bamfile = pysam.Samfile(args.bamFileName, 'rb')
        reffile = pysam.Fastafile(args.refFasta)
        logfn   = '_'.join(map(os.path.basename, bedline.strip().split())) + ".log"
        logfile = open('addsv_logs_' + os.path.basename(args.outBamFile) + '/' + os.path.basename(args.outBamFile) + '_' + logfn, 'w')
        exclfile = 'exclude.' + str(uuid4()) + '.txt'
        exclude = open(exclfile, 'w')

        # optional CNV file
        cnv = None
        if (args.cnvfile):
            cnv = pysam.Tabixfile(args.cnvfile, 'r')

        # temporary file to hold mutated reads
        outbam_mutsfile = "tmp." + str(uuid4()) + ".muts.bam"

        c = bedline.strip().split()
        chrom    = c[0]
        start  = int(c[1])
        end    = int(c[2])
        araw   = c[3:len(c)] # INV, DEL, INS seqfile.fa TSDlength, DUP
 
        actions = map(lambda x: x.strip(),' '.join(araw).split(','))

        svfrac = float(args.svfrac) # default, can be overridden by cnv file

        if cnv: # CNV file is present
            if chrom in cnv.contigs:
                for cnregion in cnv.fetch(chrom,start,end):
                    cn = float(cnregion.strip().split()[3]) # expect chrom,start,end,CN
                    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\t" + ' '.join(("copy number in snp region:",chrom,str(start),str(end),"=",str(cn))) + "\n")
                    svfrac = 1.0/float(cn)
                    assert svfrac < 1.0
                    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tadjusted MAF: " + str(svfrac) + "\n")

        print "INFO\t" + now() + "\t" + mutid + "\tinterval:", c
        print "INFO\t" + now() + "\t" + mutid + "\tlength:", end-start
        # modify start and end if interval is too long
        maxctglen = int(args.maxctglen)
        assert maxctglen > 3*int(args.maxlibsize) # maxctglen is too short
        if end-start > maxctglen:
            adj   = (end-start) - maxctglen
            rndpt = random.randint(0,adj)
            start = start + rndpt
            end   = end - (adj-rndpt)
            print "INFO\t" + now() + "\t" + mutid + "\tnote: interval size too long, adjusted:",chrom,start,end

        dfrac = discordant_fraction(args.bamFileName, chrom, start, end)
        print "INFO\t" + now() + "\t" + mutid + "\tdiscordant fraction:", dfrac

        maxdfrac = 0.1 # FIXME make a parameter
        if dfrac > .1: 
            sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\tdiscordant fraction > " + str(maxdfrac) + " aborting mutation!\n")
            return None, None

        contigs = ar.asm(chrom, start, end, args.bamFileName, reffile, int(args.kmersize), args.noref, args.recycle, mutid=mutid)

        # find the largest contig        
        maxlen = 0
        maxcontig = None
        for contig in contigs:
            if contig.len > maxlen:
                maxlen = contig.len
                maxcontig = contig

        if maxcontig is None:
            sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\tmaxcontig has length 0, aborting mutation!\n")
            return None, None

        # trim contig to get best ungapped aligned region to ref.
        refseq = reffile.fetch(chrom,start,end)
        alignstats = align(maxcontig.seq, refseq)
        
        if len(alignstats) < 6:
            sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\talignstats:" + str(alignstats) + "\n")
            sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\tNo good alignment between mutated contig and original, aborting mutation!\n")
            return None, None
        
        qrystart, qryend = map(int, alignstats[2:4])
        tgtstart, tgtend = map(int, alignstats[4:6])

        refseq = refseq[tgtstart:tgtend]

        print "INFO\t" + now() + "\t" + mutid + "\tbest contig length:", maxlen
        print "INFO\t" + now() + "\t" + mutid + "\talignment result:", alignstats

        maxcontig.trimseq(qrystart, qryend)
        print "INFO\t" + now() + "\t" + mutid + "\ttrimmed contig length:", maxcontig.len

        refstart = start + tgtstart
        refend = start + tgtend

        if refstart > refend:
            refstart, refend = refend, refstart
    
        print "INFO\t" + now() + "\t" + mutid + "\tstart, end, tgtstart, tgtend, refstart, refend:", start, end, tgtstart, tgtend, refstart, refend

        # is there anough room to make mutations?
        if maxcontig.len > 3*int(args.maxlibsize):
            # make mutation in the largest contig
            mutseq = ms.MutableSeq(maxcontig.seq)

            # support for multiple mutations
            for actionstr in actions:
                a = actionstr.split()
                action = a[0]

                print "INFO\t" + now() + "\t" + mutid + "\taction: ", actionstr, action

                insseqfile = None
                insseq = ''
                tsdlen = 0  # target site duplication length
                ndups = 0   # number of tandem dups
                dsize = 0.0 # deletion size fraction
                dlen = 0
                if action == 'INS':
                    assert len(a) > 1 # insertion syntax: INS <file.fa> [optional TSDlen]
                    insseqfile = a[1]
                    if not (os.path.exists(insseqfile) or insseqfile == 'RND'): # not a file... is it a sequence? (support indel ins.)
                        assert re.search('^[ATGCatgc]*$',insseqfile) # make sure it's a sequence
                        insseq = insseqfile.upper()
                        insseqfile = None
                    if len(a) > 2:
                        tsdlen = int(a[2])

                if action == 'DUP':
                    if len(a) > 1:
                        ndups = int(a[1])
                    else:
                        ndups = 1

                if action == 'DEL':
                    if len(a) > 1:
                        dsize = float(a[1])
                        if dsize >= 1.0: # if DEL size is not a fraction, interpret as bp
                            # since DEL 1 is default, if DEL 1 is specified, interpret as 1 bp deletion
                            dlen = int(dsize)
                            dsize = 1.0
                    else:
                        dsize = 1.0

                logfile.write(">" + chrom + ":" + str(refstart) + "-" + str(refend) + " BEFORE\n" + str(mutseq) + "\n")

                if action == 'INS':
                    if insseqfile: # seq in file
                        if insseqfile == 'RND':
                            assert args.inslib is not None # insertion library needs to exist
                            mutseq.insertion(mutseq.length()/2,pickseq(args.inslib, mutid=mutid),tsdlen)
                        else:
                            mutseq.insertion(mutseq.length()/2,singleseqfa(insseqfile, mutid=mutid),tsdlen)
                    else: # seq is input
                        mutseq.insertion(mutseq.length()/2,insseq,tsdlen)
                    logfile.write("\t".join(('ins',chrom,str(refstart),str(refend),action,str(mutseq.length()),str(mutseq.length()/2),str(insseqfile),str(tsdlen))) + "\n")

                elif action == 'INV':
                    invstart = int(args.maxlibsize)
                    invend = mutseq.length() - invstart
                    mutseq.inversion(invstart,invend)
                    logfile.write("\t".join(('inv',chrom,str(refstart),str(refend),action,str(mutseq.length()),str(invstart),str(invend))) + "\n")

                elif action == 'DEL':
                    delstart = int(args.maxlibsize)
                    delend = mutseq.length() - delstart
                    if dlen == 0: # bp size not specified, delete fraction of contig
                        dlen = int((float(delend-delstart) * dsize)+0.5) 

                    dadj = delend-delstart-dlen
                    if dadj < 0:
                        dadj = 0
                        sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\twarning: deletion of length 0\n")
    
                    delstart += dadj/2
                    delend   -= dadj/2

                    mutseq.deletion(delstart,delend)
                    logfile.write("\t".join(('del',chrom,str(refstart),str(refend),action,str(mutseq.length()),str(delstart),str(delend),str(dlen))) + "\n")

                elif action == 'DUP':
                    dupstart = int(args.maxlibsize)
                    dupend = mutseq.length() - dupstart
                    mutseq.duplication(dupstart,dupend,ndups)
                    logfile.write("\t".join(('dup',chrom,str(refstart),str(refend),action,str(mutseq.length()),str(dupstart),str(dupend),str(ndups))) + "\n")

                else:
                    raise ValueError("ERROR\t" + now() + "\t" + mutid + "\t: mutation not one of: INS,INV,DEL,DUP\n")

                logfile.write(">" + chrom + ":" + str(refstart) + "-" + str(refend) +" AFTER\n" + str(mutseq) + "\n")

            pemean, pesd = float(args.ismean), float(args.issd) 
            print "INFO\t" + now() + "\t" + mutid + "\tset paired end mean distance: " + str(args.ismean)
            print "INFO\t" + now() + "\t" + mutid + "\tset paired end distance stddev: " + str(args.issd)

            # simulate reads
            (fq1, fq2) = runwgsim(maxcontig, mutseq.seq, svfrac, exclude, pemean, pesd, mutid=mutid)

            # remap reads
            if args.bwamem:
                outreads = remap_bwamem(fq1, fq2, 1, args.refFasta, outbam_mutsfile, mutid=mutid)
            elif args.novoalign:
                outreads = remap_novoalign(fq1, fq2, 1, args.refFasta, args.novoref, outbam_mutsfile, mutid=mutid)
            else:
                outreads = remap(fq1, fq2, 1, args.refFasta, outbam_mutsfile, mutid=mutid)

            if outreads == 0:
                sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\toutbam " + outbam_mutsfile + " has no mapped reads!\n")
                return None, None

        else:
            sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\tbest contig too short to make mutation!\n")
            return None, None

        print "INFO\t" + now() + "\t" + mutid + "\ttemporary bam: " + outbam_mutsfile

        exclude.close()
        bamfile.close()

        return outbam_mutsfile, exclfile

    except Exception, e:
        sys.stderr.write("*"*60 + "\nencountered error in mutation spikein: " + bedline + "\n")
        traceback.print_exc(file=sys.stderr)
        sys.stderr.write("*"*60 + "\n")
        return None, None


def main(args):
    print "INFO\t" + now() + "\tstarting " + sys.argv[0] + " called with args: " + ' '.join(sys.argv) + "\n"
    tmpbams = [] # temporary BAMs, each holds the realigned reads for one mutation
    exclfns = [] # 'exclude' files store reads to be removed from the original BAM due to deletions

    if args.bwamem and args.novoalign:
        sys.stderr.write("ERROR\t" + now() + "\t --bwamem and --novoalign cannot be specified together")
        sys.exit(1)

    if args.novoalign and args.novoref is None:
        sys.stderr.write("ERROR\t" + now() + "\t --novoref must be be specified with --novoalign")
        sys.exit(1)

    # load insertion library if present
    try:
        if args.inslib is not None:
            print "INFO\t" + now() + "\tloading insertion library from " + args.inslib
            args.inslib = load_inslib(args.inslib)
    except Exception, e:
        sys.stderr.write("ERROR\t" + now() + "\tfailed to load insertion library " + args.inslib + "\n")
        traceback.print_exc(file=sys.stderr)
        sys.stderr.write("\n")
        sys.exit(1)

    results = []
    pool = Pool(processes=int(args.procs))

    nmuts = 0

    if not os.path.exists('addsv_logs_' + os.path.basename(args.outBamFile)):
        os.mkdir('addsv_logs_' + os.path.basename(args.outBamFile))
        print "INFO\t" + now() + "\tcreated log directory: addsv_logs_" + os.path.basename(args.outBamFile)

    with open(args.varFileName, 'r') as varfile:
        for bedline in varfile:
            if re.search('^#',bedline):
                continue
            if args.maxmuts and nmuts >= int(args.maxmuts):
                break
            
            # submit each mutation as its own thread                
            result = pool.apply_async(makemut, [args, bedline]) 
            results.append(result)                              

            nmuts += 1
            if args.delay is not None:
                sleep(int(args.delay))

    ## process the results of multithreaded mutation jobs
    for result in results:
        tmpbam = None
        exclfn = None

        tmpbam, exclfn = result.get()

        if tmpbam is not None and exclfn is not None:
            tmpbams.append(tmpbam)
            exclfns.append(exclfn)

    print "INFO\t" + now() + "\ttmpbams:",tmpbams
    print "INFO\t" + now() + "\texclude:",exclfns

    excl_merged = 'addsv.exclude.final.' + str(uuid4()) + '.txt'
    mergedtmp = 'addsv.mergetmp.final.' + str(uuid4()) + '.bam'

    print "INFO\t" + now() + "\tmerging exclude files into", excl_merged, "..."
    exclout = open(excl_merged, 'w')
    for exclfn in exclfns:
        with open(exclfn, 'r') as excl:
            for line in excl:
                exclout.write(line)
    exclout.close()

    if len(tmpbams) == 1:
        print "INFO\t" + now() + "\tonly one bam:", tmpbams[0], "renaming to", mergedtmp
        os.rename(tmpbams[0], mergedtmp)

    elif len(tmpbams) > 1:
        print "INFO\t" + now() + "\tmerging bams into", mergedtmp, "..."
        mergebams(tmpbams, mergedtmp)

    if args.skipmerge:
        print "INFO\t" + now() + "\tfinal merge skipped, please merge manually:", mergedtmp
        print "INFO\t" + now() + "\texclude file to use:", excl_merged
        print "INFO\t" + now() + "\tcleaning up..."

        for exclfn in exclfns:
            if os.path.isfile(exclfn):
                os.remove(exclfn)

        for tmpbam in tmpbams:
            if os.path.isfile(tmpbam):
                os.remove(tmpbam)
            if os.path.isfile(tmpbam + '.bai'):
                os.remove(tmpbam + '.bai')

    else:
        print "INFO\t" + now() + "\tswapping reads into original and writing to", args.outBamFile
        replace(args.bamFileName, mergedtmp, args.outBamFile, excl_merged)

        os.remove(excl_merged)
        os.remove(mergedtmp)

        for exclfn in exclfns:
            if os.path.isfile(exclfn):
                os.remove(exclfn)

        for tmpbam in tmpbams:
            if os.path.isfile(tmpbam):
                os.remove(tmpbam)
            if os.path.isfile(tmpbam + '.bai'):
                os.remove(tmpbam + '.bai')

        print "INFO\t" + now() + "\tdone."

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='adds SNVs to reads, outputs modified reads as .bam along with mates')
    parser.add_argument('-v', '--varfile', dest='varFileName', required=True,
                        help='whitespace-delimited target regions to try and add a SNV: chrom,start,stop,action,seqfile (if insertion),TSDlength (if insertion)')
    parser.add_argument('-f', '--bamfile', dest='bamFileName', required=True,
                        help='sam/bam file from which to obtain reads')
    parser.add_argument('-r', '--reference', dest='refFasta', required=True,
                        help='reference genome, fasta indexed with bwa index -a stdsw _and_ samtools faidx')
    parser.add_argument('-o', '--outbam', dest='outBamFile', required=True,
                        help='.bam file name for output')
    parser.add_argument('-l', '--maxlibsize', dest='maxlibsize', default=600, help="maximum fragment length of seq. library")
    parser.add_argument('-k', '--kmer', dest='kmersize', default=31, 
                        help="kmer size for assembly (default = 31)")
    parser.add_argument('-s', '--svfrac', dest='svfrac', default=1.0, 
                        help="allele fraction of variant (default = 1.0)")
    parser.add_argument('--maxctglen', dest='maxctglen', default=32000, 
                        help="maximum contig length for assembly - can increase if velvet is compiled with LONGSEQUENCES")
    parser.add_argument('-n', dest='maxmuts', default=None,
                        help="maximum number of mutations to make")
    parser.add_argument('-c', '--cnvfile', dest='cnvfile', default=None, 
                        help="tabix-indexed list of genome-wide absolute copy number values (e.g. 2 alleles = no change)")
    parser.add_argument('--ismean', dest='ismean', default=300, help="mean insert size (default = estimate from region)")
    parser.add_argument('--issd', dest='issd', default=70, help="insert size standard deviation (default = estimate from region)")
    parser.add_argument('-p', '--procs', dest='procs', default=1, help="split into multiple processes (default=1)")
    parser.add_argument('--inslib', default=None, help='FASTA file containing library of possible insertions, use INS RND instead of INS filename to pick one')
    parser.add_argument('--delay', default=None, help='time delay between jobs (try to avoid thrashing disks)')
    parser.add_argument('--nomut', action='store_true', default=False, help="dry run")
    parser.add_argument('--noremap', action='store_true', default=False, help="dry run")
    parser.add_argument('--noref', action='store_true', default=False, 
                        help="do not perform reference based assembly")
    parser.add_argument('--recycle', action='store_true', default=False)
    parser.add_argument('--bwamem', action='store_true', default=False, help='realign with bwa mem (original shuld be aligned with mem as well!)')
    parser.add_argument('--novoalign', action='store_true', default=False, help='realignment with novoalign')
    parser.add_argument('--novoref', default=None, help='novoalign reference, must be specified with --novoalign')
    parser.add_argument('--skipmerge', action='store_true', default=False)
    args = parser.parse_args()
    main(args)


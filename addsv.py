#!/usr/bin/env python

import re, os, sys, random
import subprocess
import traceback
import argparse
import pysam
import tempfile
import bs.replacereads as rr
import bs.asmregion as ar
import bs.mutableseq as ms
from itertools import izip
from collections import Counter
from multiprocessing import Pool

def remap(fq1, fq2, threads, bwaref, outbam, deltmp=True):
    """ call bwa/samtools to remap .bam and merge with existing .bam
    """
    basefn = "bwatmp" + str(random.random())
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

    print "mapping 1st end, cmd: " + " ".join(sai2args)
    p = subprocess.Popen(sai2args, stderr=subprocess.STDOUT)
    p.wait()

    print "mapping 2nd end, cmd: " + " ".join(sai1args)
    p = subprocess.Popen(sai1args, stderr=subprocess.STDOUT)
    p.wait()

    print "pairing ends, building .sam, cmd: " + " ".join(samargs)
    p = subprocess.Popen(samargs, stderr=subprocess.STDOUT)
    p.wait()

    print "sam --> bam, cmd: " + " ".join(bamargs)
    p = subprocess.Popen(bamargs, stderr=subprocess.STDOUT)
    p.wait()

    print "sorting, cmd: " + " ".join(sortargs)
    p = subprocess.Popen(sortargs, stderr=subprocess.STDOUT)
    p.wait()

    print "rename " + tmpsrt + ".bam --> " + tmpbam
    os.remove(tmpbam)
    os.rename(tmpsrt + ".bam", tmpbam)

    print "rename " + tmpbam + " --> " + outbam
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

def set_isize(isizedistr):
    print(isizedistr)
    if isizedistr == ['']: #use default wgsim isize
        return None
    else:
        isizeopt=[]  #each file distr is consisted of components [mixrate,mean,sd] triples
        for distr in isizedistr:
            distr=distr.split(',')
            distropt=[]
            for mode in distr:
                #print mode, distr
                distropt.append([float(v) for v in mode.split('_')])
            isizeopt.append(distropt)
    #print isizeopt
    return isizeopt            

def runwgsim(contig,newseq,svfrac,exclude,readnames=None,isizedistr=None):
    ''' wrapper function for wgsim
    '''
    assert len(exclude) == len(readnames)
    nbf=len(exclude) #number of bamfiles

    try:
        print isizedistr
        isizeopt = set_isize(isizedistr)
        if isizeopt == None: #thst means only one default distr for each bam
            isizeopt = [ [ [1.0, 500, 50] ] for x in xrange(nbf) ]
    except Exception, e:
        print >>sys.stderr,"parsing isizedistr error!"

    assert len(isizeopt) == len(readnames)
    #print readnames

    namecount = Counter(contig.reads.reads)

    #need to identify how many each type are from which bam
    totalreads = len(contig.reads.reads)
    paired = [ 0 ] * nbf
    single = [ 0 ] * nbf
    discard = [ 0 ] * nbf
    pairednames = [ [] for x in xrange(nbf) ]
    # names with count 2 had both pairs in the contig
    for name,count in namecount.items():
        #print name,count
        origin = [ rni for rni in xrange(len(readnames)) if name in readnames[rni] ]
        try:
            assert len(origin) == 1 #otherwise where this read come from?
        except AssertionError:
            print "warning: some read in contig not in bam!"
            continue
        origin = origin[0]
        if count == 1:
            single[origin] += 1
        elif count == 2:
            paired[origin] += 1 
            pairednames[origin].append(name) 
        else:
            discard[origin] += 1

    basefn = "wgsimtmp" + str(random.random()) 
    fasta = basefn + ".fasta"
    fout = open(fasta,'w')
    fout.write(">target\n" + newseq + "\n")
    fout.close()
    fq1 = [ basefn + "." + str(bi) + ".1.fq" for bi in xrange(nbf) ]
    fq2 = [ basefn + "." + str(bi) + ".2.fq" for bi in xrange(nbf) ]
    tfq1 = []
    tfq2 = []
    for bi in xrange(nbf):
        nfq = len(isizeopt[bi])
        tfq1.append([ basefn + "." + str(bi) + "." + str(fqi) + ".fq1" for fqi in xrange(nfq) ])
        tfq2.append([ basefn + "." + str(bi) + "." + str(fqi) + ".fq2" for fqi in xrange(nfq) ])


    print "paired : " + str(sum(paired))
    print "single : " + str(sum(single))
    print "discard: " + str(sum(discard))
    print "total  : " + str(totalreads)

    # adjustment factor for length of new contig vs. old contig
    lenfrac = float(len(newseq))/float(len(contig.seq))

    print "old ctg len: " + str(len(contig.seq))
    print "new ctg len: " + str(len(newseq))
    print "adj. factor: " + str(lenfrac)

    # number of paried reads to simulate
    nsimreads = int((sum(paired) + (sum(single)/2)) * svfrac * lenfrac)

    print "num. sim. reads: " + str(nsimreads) 

    # length of quality score comes from original read, used here to set length of read
    maxqlen = 0
    for qual in (contig.rquals + contig.mquals):
        if len(qual) > maxqlen:
            maxqlen = len(qual)

    print "paired reads per bam:", paired
    print "single reads per bam:", single
    print "isize option per bam:", isizeopt

    for bi in xrange(nbf):
        bireads = nsimreads*(float(paired[bi] + single[bi]/2)/float(sum(paired)+sum(single)/2))
        print "sim reads in", bi, "-th bam:", bireads
        nfq = len(isizeopt[bi]) #number of modes in isize
        for fqi in xrange(nfq):
            fqreads = round(isizeopt[bi][fqi][0] * bireads)
            fqisize = isizeopt[bi][fqi][1]
            fqisizesd = isizeopt[bi][fqi][2]
            tfq1out = tfq1[bi][fqi]
            tfq2out = tfq2[bi][fqi]
            args = ['wgsim','-e','0','-N',str(fqreads),'-1',str(maxqlen),'-2','100','-r','0','-R','0','-d',str(fqisize),'-s',str(fqisizesd),fasta,tfq1out,tfq2out]
            print args
            subprocess.call(args)

    os.remove(fasta)

    for bi in xrange(nbf):
        merge_fq(fq1[bi],tfq1[bi]) #temporary fq files are merged to one
        merge_fq(fq2[bi],tfq2[bi])
        fqReplaceList(fq1[bi],pairednames[bi],contig.rquals,svfrac,exclude[bi])
        fqReplaceList(fq2[bi],pairednames[bi],contig.mquals,svfrac,exclude[bi])

    return (fq1,fq2)

def merge_fq(dest,source):
    if os.path.exists(dest):
        os.remove(dest)
    for idx in xrange(len(source)):
        args = ['cat',source[idx],'>>'+dest ]
        #subprocess.Popen(" ".join(args),shell=True)
        os.system(" ".join(args))

def fqReplaceList(fqfile,names,quals,svfrac,exclude):
    '''
    Replace seq names in paired fastq files from a list until the list runs out
    (then stick with original names). fqfile = fastq file, names = list

    'exclude' is a filehandle, the exclude file contains read names that should
    not appear in the final output BAM

    '''
    print fqfile
    print os.path.exists(fqfile)
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
            raise ValueError("fastq iteration problem")

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
            print "warning, used read name: " + newnames[i] + " in multiple pairs"
            #charlie: according to Adam's coding, this warning shouldn't happen at all
        usednames[newnames[i]] = True
        
    # burn off excess
    if len(seqs) > 0:
        for name in names:
            if name not in usednames:
                if random.uniform(0,1) < svfrac:  # this controls deletion depth
                    exclude.write(name + "\n")

    fqout.close()

def singleseqfa(file):
    with open(file, 'r') as fasta:
        header = None
        seq = ''
        for line in fasta:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    sys.stderr.write("multiple entries found in " + file + " only using the first\n")
                header = line.lstrip('>')
            else:
                seq += line
    return seq

def mergebams(bamlist, outbamfn, maxopen=100):
    ''' call samtools to merge a list of bams hierarchically '''

    assert outbamfn.endswith('.bam')
    print "DEBUG: len(bamlist)", len(bamlist)

    if len(bamlist) == 1:
        print "only one BAM to merge, renaming",bamlist[0],"-->",outbamfn
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
                print "renamed:",tmpbams[0],"-->",submergefn
            else:
                args = ['samtools','merge','-f',submergefn] + tmpbams 
                print "merging, cmd: ",args
                subprocess.call(args)

        if len(merge_sublists.keys()) == 1:
            print "merge finished, renaming:",merge_sublists.keys()[0],"-->",outbamfn
            os.rename(merge_sublists.keys()[0], outbamfn)
        else:
            args = ['samtools','merge','-f',outbamfn] + merge_sublists.keys() 
            print "final merge, cmd: ",args
            subprocess.call(args)

        for submergefn in merge_sublists.keys():
            if os.path.exists(submergefn):
                os.remove(submergefn)

    for bamfile in bamlist:
        if os.path.exists(bamfile):
            os.remove(bamfile)
            os.remove(bamfile + '.bai')

def align(qryseq, refseq):
    rnd = str(random.random())
    tgtfa = 'tmp.' + rnd + '.tgt.fa'
    qryfa = 'tmp.' + rnd + '.qry.fa'
    #print tgtfa
    #print qryfa

    tgt = open(tgtfa, 'w')
    qry = open(qryfa, 'w')

    tgt.write('>ref' + '\n' + refseq + '\n')
    qry.write('>qry' + '\n' + qryseq + '\n')

    tgt.close()
    qry.close()

    cmd = ['exonerate', '--bestn', '1', '-m', 'ungapped', '--showalignment','0', '--ryo', 'SUMMARY\t%s\t%qab\t%qae\t%tab\t%tae\n', '--query', qryfa, '--target', tgtfa]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print " ".join(cmd)

    best = []
    topscore = 0

    for pline in p.stdout.readlines():
        print pline
        if pline.startswith('SUMMARY'):
            c = pline.strip().split()
            if int(c[1]) > topscore:
                topscore = int(c[1])
                best = c

    #os.remove(tgtfa)
    #os.remove(qryfa)

    return best

def check_asmvariants(bamfile, contigseq, reffile, chrom, refstart, refend):
    ''' check contig for variants, add full depth variants if they're present in the original bam but not the local assembly '''

    refbases = list(reffile.fetch(chrom, refstart, refend))
    ctgbases = list(contigseq)

    offset = 0
    numchanges = 0
    changepos = []

    region = chrom + ':' + str(refstart) + '-' + str(refend)
    lastpos  = refstart-1
    mpargs = ['samtools', 'mpileup', '-q','10','-r', region, '-f', reffile.filename, bamfile]

    print "pileup:", mpargs

    p = subprocess.Popen(mpargs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in p.stdout.readlines():
        mpileupstr=''
        c = line.strip().split()
        if len(c) == 6:
            refpos     = int(c[1])
            refbase    = c[2].upper()
            depth      = int(c[3])
            mpileupstr = c[4]

            assert refpos > lastpos

            baselist = []
            for base in list(mpileupstr):
                base = base.upper()
                if base in ['A', 'T', 'G', 'C']:
                    baselist.append(base)

            basecount = Counter(baselist)

            for base, count in basecount.items():
                if count > depth/2:
                    ctgbases[offset] = base
                    numchanges += 1
                    changepos.append(refpos)

            offset += (refpos-lastpos)
            lastpos = refpos

    if numchanges > len(contigseq) * .10: # could set as a parameter
        sys.stderr.write("*** Whoops, I mangled a contig! Returning original... *** " + region + "\n")
        return contigseq 

    print "contig editing successful, made", numchanges, "changes."
    print "change positions: " + ','.join(map(str, changepos))
    return ''.join(ctgbases) 

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

def makemut(args, bedline, pct=1):
    #varfile = open(args.varFileName, 'r')

    try:
        #bamfile = pysam.Samfile(args.bamFileName, 'rb')
        #exclfile = 'exclude.' + str(random.random()) + '.txt'
        #exclude = open(exclfile, 'w')
        #logfile = open('addsv_logs_' + os.path.basename(args.outBamFile) + '/' + os.path.basename(args.outBamFile) + '_' + logfn, 'w')
        bamfilenames = args.bamFileName.split(':')
        isizedistr = args.isizeDistr.split(':')
        bamfile = [ pysam.Samfile(bfn, 'rb') for bfn in bamfilenames ]
        reffile = pysam.Fastafile(args.refFasta)
        logfn   = '_'.join(map(os.path.basename, bedline.strip().split())) + ".log"
        outbamfilenames = args.outBamFile.split(':')
        logfile = [ open('addsv_logs_' + os.path.basename(obfn) + '/' + os.path.basename(obfn) + '_' + logfn, 'w') for obfn in outbamfilenames ]
        exclfile = [ 'exclude.' + str(random.random()) + '.txt' for bfn in bamfilenames ]
        exclude = [ open(ef, 'w') for ef in exclfile ]
        assert len(bamfile) == len(outbamfilenames)
        nbf = len(bamfile)
        #charlie: bamfile is a list now

        # temporary file to hold mutated reads
        #outbam_mutsfile = "tmp." + str(random.random()) + ".muts.bam"
        outbam_mutsfile = [ "tmp." + str(random.random()) + ".muts.bam" for bfn in bamfilenames ]

        c = bedline.strip().split()
        chrom    = c[0]
        start  = int(c[1])
        end    = int(c[2])
        araw   = c[3:len(c)] # INV, DEL delen, INS seq/seqfile.fa, DUP times
 
        actions = map(lambda x: x.strip(),' '.join(araw).split(','))

        if pct is not None:
            svfrac = float(pct)
        else:
            svfrac = float(args.svfrac) # default, can be overridden by cnv file
        assert svfrac>0


        print "interval:", c
        print "length:", end-start
        # modify start and end if interval is too long
        maxctglen = int(args.maxctglen)
        assert maxctglen > 3*int(args.maxlibsize) # maxctglen is too short
        if end-start > maxctglen:
            adj   = (end-start) - maxctglen
            rndpt = random.randint(0,adj)
            start = start + rndpt
            end   = end - (adj-rndpt)
            print "note: interval size too long, adjusted:",chrom,start,end

        #dfrac = discordant_fraction(args.bamFileName, chrom, start, end)
        dfrac = [ discordant_fraction(bfn, chrom, start, end) for bfn in bamfilenames ]
        #charlie: dfrac is a list now
        print "discordant fraction:",dfrac

        maxdfrac = 0.1 # FIXME make a parameter
        #if dfrac > .1: 
        if any([df>maxdfrac for df in dfrac]):
            print "discordant fraction >", maxdfrac, "aborting mutation!"
            return None, None

        #need to keep track of all read names in the region
        readnames=[ [] for x in xrange(nbf) ]
        for bi in xrange(nbf):
            for read in bamfile[bi].fetch(chrom,start,end):
                readnames[bi].append(read.qname)

        #contigs = ar.asm(chrom, start, end, args.bamFileName, reffile, int(args.kmersize), args.noref, args.recycle)
        contigs = ar.asm(chrom, start, end, bamfilenames, reffile, int(args.kmersize), args.noref, args.recycle)
        #charlie: need assemble of multiple bam files and return one set of contigs
        #charlie: we still need to find a way to test each bam's contribution portion of reads

        # find the largest contig        
        maxlen = 0
        maxcontig = None
        for contig in contigs:
            if contig.len > maxlen:
                maxlen = contig.len
                maxcontig = contig

        if maxcontig is None:
            print "maxcontig has length 0, aborting mutation!"
            return None, None

        # trim contig to get best ungapped aligned region to ref.
        refseq = reffile.fetch(chrom,start,end)
        if refseq=="":
            print "empty reference seq, maybe seqname err"

        alignstats = align(maxcontig.seq, refseq)
        print alignstats

        if len(alignstats) < 6:
            print "alignstats:", alignstats
            print "No good alignment between mutated contig and original, aborting mutation!"
            return None, None
  
        qrystart, qryend = map(int, alignstats[2:4])
        tgtstart, tgtend = map(int, alignstats[4:6])

        refseq = refseq[tgtstart:tgtend]

        print "best contig length:", maxlen
        print "alignment result:", alignstats

        maxcontig.trimseq(qrystart, qryend)
        print "trimmed contig length:", maxcontig.len

        refstart = start + tgtstart
        refend = start + tgtend

        if refstart > refend:
            refstart, refend = refend, refstart
    
        print 'start, end, tgtstart, tgtend, refstart, refend:', start, end, tgtstart, tgtend, refstart, refend

        #fixedseq = check_asmvariants(args.bamFileName, maxcontig.seq, reffile, chrom, refstart, refend)
        fixedseq = maxcontig.seq  # FIXME

        # is there anough room to make mutations?
        if maxcontig.len > 3*int(args.maxlibsize):
            # make mutation in the largest contig
            mutseq = ms.MutableSeq(fixedseq)

            # support for multiple mutations
            for actionstr in actions:
                a = actionstr.split()
                action = a[0]

                print actionstr,action

                insseqfile = None
                insseq = ''
                tsdlen = 0  # target site duplication length
                ndups = 0   # number of tandem dups
                dsize = 0.0 # deletion size fraction
                dlen = 0
                if action == 'INS':
                    assert len(a) > 1 # insertion syntax: INS <file.fa> [optional TSDlen]
                    insseqfile = a[1]
                    if not os.path.exists(insseqfile): # not a file... is it a sequence? (support indel ins.)
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

                for lgf in logfile:
                    lgf.write(">" + chrom + ":" + str(refstart) + "-" + str(refend) + " BEFORE\n" + str(mutseq) + "\n")

                if action == 'INS':
                    if insseqfile: # seq in file
                        mutseq.insertion(mutseq.length()/2,singleseqfa(insseqfile),tsdlen)
                    else: # seq is input
                        mutseq.insertion(mutseq.length()/2,insseq,tsdlen)

                    for lgf in logfile:
                        lgf.write("\t".join(('ins',chrom,str(refstart),str(refend),action,str(mutseq.length()),str(mutseq.length()/2),str(insseqfile),str(tsdlen),str(svfrac)))+ "\n")

                elif action == 'INV':
                    invstart = int(args.maxlibsize)
                    invend = mutseq.length() - invstart
                    mutseq.inversion(invstart,invend)
                    for lgf in logfile:
                        lgf.write("\t".join(('inv',chrom,str(refstart),str(refend),action,str(mutseq.length()),str(invstart),str(invend),str(invend-invstart),str(svfrac))) + "\n")

                elif action == 'DEL':
                    delstart = int(args.maxlibsize)
                    delend = mutseq.length() - delstart
                    if dlen == 0: # bp size not specified, delete fraction of contig
                        dlen = int((float(delend-delstart) * dsize)+0.5) 

                    dadj = delend-delstart-dlen
                    if dadj < 0:
                        dadj = 0
                        print "warning: deletion of length 0"
    
                    delstart += dadj/2
                    delend   -= dadj/2

                    mutseq.deletion(delstart,delend)
                    for lgf in logfile:
                        lgf.write("\t".join(('del',chrom,str(refstart),str(refend),action,str(mutseq.length()),str(delstart),str(delend),str(dlen),str(svfrac))) + "\n")

                elif action == 'DUP':
                    dupstart = int(args.maxlibsize)
                    dupend = mutseq.length() - dupstart
                    mutseq.duplication(dupstart,dupend,ndups)
                    for lgf in logfile:
                        lgf.write("\t".join(('dup',chrom,str(refstart),str(refend),action,str(mutseq.length()),str(dupstart),str(dupend),str(ndups),str(svfrac))) + "\n")

                else:
                    raise ValueError(bedline.strip() + ": mutation not one of: INS,INV,DEL,DUP")

                for lgf in logfile:
                    lgf.write(">" + chrom + ":" + str(refstart) + "-" + str(refend) +" AFTER\n" + str(mutseq) + "\n")

            # simulate reads
            (fq1, fq2) = runwgsim(maxcontig, mutseq.seq, svfrac, exclude, readnames, isizedistr)
            # charlie: here we need to runwgsim for multiple bamfile
            # charlie: fq1, fq2 is list of len(bamfile)
            # charlie: exclude is list of len(bamfile)

            # remap reads
            for bi in xrange(nbf):
                outreads = remap(fq1[bi], fq2[bi], 4, args.refFasta, outbam_mutsfile[bi], deltmp=False)
            # charlie: this remap function also need update to reflect multi bams
                if outreads == 0:
                    print "outbam", outbam_mutsfile[bi], "has no mapped reads!"
                    return None, None

        else:
            print "best contig too short to make mutation: ",bedline.strip()
            return None, None

        sys.stderr.write("temporary bam: " + str(outbam_mutsfile) + "\n")

        for bi in xrange(nbf):
            exclude[bi].close()
            bamfile[bi].close()

        return outbam_mutsfile, exclfile

    except Exception, e:
        sys.stderr.write("*"*60 + "\nencountered error in mutation spikein: " + bedline + "\n")
        traceback.print_exc(file=sys.stdout)
        sys.stderr.write("*"*60 + "\n")
        return None, None

def main(args):

    outbamfilenames = args.outBamFile.split(':')
    bamfilenames = args.bamFileName.split(':')
    nbf = len(bamfilenames)

    tmpbams = [ [] for x in xrange(nbf) ]
    exclfns = [ [] for x in xrange(nbf) ] # 'exclude' files

    results = []
    pool = Pool(processes=int(args.procs))

    nmuts = 0

    for bi in xrange(nbf):
        if not os.path.exists('addsv_logs_' + os.path.basename(outbamfilenames[bi])):
            os.mkdir('addsv_logs_' + os.path.basename(outbamfilenames[bi]))
            print "created log directory: addsv_logs_" + os.path.basename(outbamfilenames[bi])
    
    pct=[]; npct=0
    if args.varPctFile is not None:
        with open(args.varPctFile, 'r') as pctfile:
            for pctline in pctfile:
                if re.search('^#',pctline):
                    continue
                if args.maxmuts and npct >= int(args.maxmuts):
                    break
                pct.append(pctline.rstrip())
                npct += 1

    with open(args.varFileName, 'r') as varfile:
        warnstring=None
        for bedline in varfile:
            if re.search('^#',bedline):
                continue
            if args.maxmuts and nmuts >= int(args.maxmuts):
                break
            #tmpbam, exclude = makemut(args, bedline)
            try:
                bedpct=pct[nmuts]
            except Exception,e:
                warnstring="WARNING: running out of percentage for some mutations, forcing None!"
                bedpct=None
            result = pool.apply_async(makemut, [args, bedline, bedpct]) #MT
            results.append(result)                              #MT
            #tmpbams.append(tmpbam)
            #exclfns.append(exclude)
            nmuts += 1
        if warnstring is not None:
            print warnstring
    
    # charlie:
    # optional CNV file, this basically covered by DUP
    # guess no longer needed
    # should mark bosolete later
    # or we need to generate a bedline of DUP and feed into makemut here 
    #cnv = None
    #if (args.cnvfile):
    #    cnv = pysam.Tabixfile(args.cnvfile, 'r')

    #if cnv: # CNV file is present
    #    if chrom in cnv.contigs:
    #        for cnregion in cnv.fetch(chrom,start,end):
    #            cn = float(cnregion.strip().split()[3]) # expect chrom,start,end,CN,pct
    #            try:
    #                pct = float(cnregion.strip().split()[4]) # expect chrom,start,end,CN,pct
    #                pct = cn * pct 
    #            except Exception,e: #not specified
    #                pct = float(args.svfrac) * cn
    #            sys.stderr.write(' '.join(("copy number in snp region:",chrom,str(start),str(end),"=",str(cn))),"pct=", str(pct)+ "\n")
    #            #charlie: what is svfrac?
    #            #charlie: is this a bug?
    #            #a CNV of 3 should have 3 times more reads in the same region
    #            #so the svfrace should be 3 not 1/3?
    #            sys.stderr.write("adjusted MAF: " + str(svfrac) + "\n")

    ## MT ##
    for result in results:
        tmpbam = None
        exclfn = None

        tmpbam, exclfn = result.get()

        print tmpbam
        print exclfn
        if tmpbam is not None and exclfn is not None:
            for bi in xrange(nbf):
                #print tmpbams[bi],tmpbam[bi]
                tmpbams[bi].append(tmpbam[bi])
                exclfns[bi].append(exclfn[bi])
    ########

    print "bams:",str(tmpbams)
    print "excl:",str(exclfns)

    def final_merge(tmpbams,exclfns,bamfile,outbamfile):

        excl_merged = 'exclude.final.' + str(random.random()) + '.txt'
        mergedtmp = 'mergetmp.final.' + str(random.random()) + '.bam'

        print "merging exclude files into", excl_merged, "..."
        exclout = open(excl_merged, 'w')
        for exclfn in exclfns:
            with open(exclfn, 'r') as excl:
                for line in excl:
                    exclout.write(line)
        exclout.close()

        if len(tmpbams) == 1:
            print "only one bam:", tmpbams[0], "renaming to", mergedtmp
            os.rename(tmpbams[0], mergedtmp)

        elif len(tmpbams) > 1:
            print "merging bams into", mergedtmp, "..."
            mergebams(tmpbams, mergedtmp)

        if args.skipmerge:
            print "final merge skipped, please merge manually:",mergedtmp
            print "exclude file to use:",excl_merged
            print "cleaning up..."

            for exclfn in exclfns:
                if os.path.isfile(exclfn):
                    os.remove(exclfn)

            for tmpbam in tmpbams:
                if os.path.isfile(tmpbam):
                    os.remove(tmpbam)
                if os.path.isfile(tmpbam + '.bai'):
                    os.remove(tmpbam + '.bai')

        else:
            print "swapping reads into original and writing to", outbamfile
            print bamfile, mergedtmp, outbamfile, excl_merged
            replace(bamfile, mergedtmp, outbamfile, excl_merged)

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
            print "done."

    for bi in xrange(nbf):
        if tmpbams[bi] != []:
            final_merge(tmpbams[bi],exclfns[bi],bamfilenames[bi],outbamfilenames[bi])
        else:
            print "no bam to merge for",bi,"-th lib bam, all sv spike-in failed on it?\n" 

def run():
    parser = argparse.ArgumentParser(description='adds SNVs to reads, outputs modified reads as .bam along with mates')
    parser.add_argument('-v', '--varfile', dest='varFileName', required=True,
                        help='whitespace-delimited target regions to try and add a SNV: chrom,start,stop,action,seqfile if insertion,TSDlength if insertion')
    parser.add_argument('-f', '--sambamfile', dest='bamFileName', required=True,
                        help='sam/bam files from which to obtain reads, f1.bam:f2.bam:...')
    parser.add_argument('-i', '--isizedistr', dest='isizeDistr', required=False, default="",
                        help='insert size distribution for each bam files, d11_p111_p112,d12_p121_p122:d21_p211_p222:...')
    parser.add_argument('-t', '--varpctfile', dest='varPctFile', required=False, default=None,
                        help='a file to specify variants percentage')
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
    parser.add_argument('-p', '--procs', dest='procs', default=1, help="split into multiple processes (default=1)")
    parser.add_argument('--nomut', action='store_true', default=False, help="dry run")
    parser.add_argument('--noremap', action='store_true', default=False, help="dry run")
    parser.add_argument('--noref', action='store_true', default=False, 
                        help="do not perform reference based assembly")
    parser.add_argument('--recycle', action='store_true', default=False)
    parser.add_argument('--skipmerge', action='store_true', default=False)
    args = parser.parse_args()
    main(args)

if __name__ == '__main__':
    run()

#!/usr/bin/env python

import re, os, sys, random
import subprocess
import traceback
import argparse
import pysam
import bs.replacereads as rr
import bs.asmregion as ar
import bs.mutableseq as ms

from math import sqrt
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

    sai1args = ['bwa', 'aln', bwaref, '-q', '5', '-l', '32', '-k', '2', '-t', str(threads), '-o', '1', '-f', sai1fn, fq1]
    sai2args = ['bwa', 'aln', bwaref, '-q', '5', '-l', '32', '-k', '2', '-t', str(threads), '-o', '1', '-f', sai2fn, fq2]
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


def runwgsim(contig, newseq, svfrac, exclude, pemean, pesd):
    ''' wrapper function for wgsim
    '''
    namecount = Counter(contig.reads.reads)

    basefn = "wgsimtmp" + str(random.random())
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

    print "paired :", paired
    print "single :", single
    print "discard:", discard
    print "total  :", totalreads

    # adjustment factor for length of new contig vs. old contig
    lenfrac = float(len(newseq))/float(len(contig.seq))

    print "old ctg len:", len(contig.seq)
    print "new ctg len:", len(newseq)
    print "adj. factor:", lenfrac

    # number of paried reads to simulate
    nsimreads = int((paired + (single/2)) * svfrac * lenfrac)

    print "num. sim. reads:", nsimreads 
    print "PE mean outer distance:", pemean
    print "PE outer distance SD:", pesd

    # length of quality score comes from original read, used here to set length of read
    maxqlen = 0
    for qual in (contig.rquals + contig.mquals):
        if len(qual) > maxqlen:
            maxqlen = len(qual)

    args = ['wgsim','-e','0','-d',str(pemean),'-s',str(pesd),'-N',str(nsimreads),'-1',str(maxqlen),'-2','100','-r','0','-R','0',fasta,fq1,fq2]
    print args
    subprocess.call(args)

    os.remove(fasta)

    fqReplaceList(fq1,pairednames,contig.rquals,svfrac,exclude)
    fqReplaceList(fq2,pairednames,contig.mquals,svfrac,exclude)

    return (fq1,fq2)


def fqReplaceList(fqfile,names,quals,svfrac,exclude):
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

#TODO: replace with median ... mean is almost always driven by a few read pairs with wacky distances marked proper...
#def mean(d):
#    return float(sum(d)) / float(len(d))


#def sd(d):
#    return sqrt(mean(map(lambda x: (x - mean(d))**2.0, d)))


def estimate_pedist(bamfile, chrom, start, end, window=10000, setmean=None, setsd=None):
    start -= window
    end   += window

    if start < 0:
        start = 0

    if setmean is not None and setsd is not None:
        return float(setmean), float(setsd)

    outerlen = [read.tlen + 2*len(read.seq) for read in bamfile.fetch(chrom, start-window, end+window) if read.is_proper_pair]
    if setmean is None:
        setmean = mean(outerlen)
    if setsd is None:
        setsd = sd(outerlen)
    return float(setmean), float(setsd)


def makemut(args, bedline):
    try:
        bamfile = pysam.Samfile(args.bamFileName, 'rb')
        reffile = pysam.Fastafile(args.refFasta)
        logfn   = '_'.join(map(os.path.basename, bedline.strip().split())) + ".log"
        logfile = open('addsv_logs_' + os.path.basename(args.outBamFile) + '/' + os.path.basename(args.outBamFile) + '_' + logfn, 'w')
        exclfile = 'exclude.' + str(random.random()) + '.txt'
        exclude = open(exclfile, 'w')

        # optional CNV file
        cnv = None
        if (args.cnvfile):
            cnv = pysam.Tabixfile(args.cnvfile, 'r')

        # temporary file to hold mutated reads
        outbam_mutsfile = "tmp." + str(random.random()) + ".muts.bam"

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
                    sys.stderr.write(' '.join(("copy number in snp region:",chrom,str(start),str(end),"=",str(cn))) + "\n")
                    svfrac = 1.0/float(cn)
                    sys.stderr.write("adjusted MAF: " + str(svfrac) + "\n")

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

        dfrac = discordant_fraction(args.bamFileName, chrom, start, end)
        print "discordant fraction:",dfrac

        maxdfrac = 0.1 # FIXME make a parameter
        if dfrac > .1: 
            print "discordant fraction >", maxdfrac, "aborting mutation!"
            return None, None

        contigs = ar.asm(chrom, start, end, args.bamFileName, reffile, int(args.kmersize), args.noref, args.recycle)

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
        alignstats = align(maxcontig.seq, refseq)
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

                logfile.write(">" + chrom + ":" + str(refstart) + "-" + str(refend) + " BEFORE\n" + str(mutseq) + "\n")

                if action == 'INS':
                    if insseqfile: # seq in file
                        mutseq.insertion(mutseq.length()/2,singleseqfa(insseqfile),tsdlen)
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
                        print "warning: deletion of length 0"
    
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
                    raise ValueError(bedline.strip() + ": mutation not one of: INS,INV,DEL,DUP")

                logfile.write(">" + chrom + ":" + str(refstart) + "-" + str(refend) +" AFTER\n" + str(mutseq) + "\n")

            # estimate paired-end distribution
            print "estimating paired-end insert size mean, stdev..."
            pemean, pesd = estimate_pedist(bamfile, chrom, start, end, window=10000, setmean=args.ismean, setsd=args.issd)
            # simulate reads
            (fq1, fq2) = runwgsim(maxcontig, mutseq.seq, svfrac, exclude, pemean, pesd)

            # remap reads
            outreads = remap(fq1, fq2, 4, args.refFasta, outbam_mutsfile)

            if outreads == 0:
                print "outbam", outbam_mutsfile, "has no mapped reads!"
                return None, None

        else:
            print "best contig too short to make mutation: ",bedline.strip()
            return None, None

        sys.stderr.write("temporary bam: " + outbam_mutsfile + "\n")

        exclude.close()
        bamfile.close()

        return outbam_mutsfile, exclfile

    except Exception, e:
        sys.stderr.write("*"*60 + "\nencountered error in mutation spikein: " + bedline + "\n")
        traceback.print_exc(file=sys.stdout)
        sys.stderr.write("*"*60 + "\n")
        return None, None

def main(args):
    tmpbams = []
    exclfns = [] # 'exclude' files

    results = []
    pool = Pool(processes=int(args.procs))

    nmuts = 0

    if not os.path.exists('addsv_logs_' + os.path.basename(args.outBamFile)):
        os.mkdir('addsv_logs_' + os.path.basename(args.outBamFile))
        print "created log directory: addsv_logs_" + os.path.basename(args.outBamFile)

    with open(args.varFileName, 'r') as varfile:
        for bedline in varfile:
            if re.search('^#',bedline):
                continue
            if args.maxmuts and nmuts >= int(args.maxmuts):
                break
            #tmpbam, exclude = makemut(args, bedline)
            result = pool.apply_async(makemut, [args, bedline]) #MT
            results.append(result)                              #MT
            #tmpbams.append(tmpbam)
            #exclfns.append(exclude)
            nmuts += 1

    ## MT ##
    for result in results:
        tmpbam = None
        exclfn = None

        tmpbam, exclfn = result.get()

        if tmpbam is not None and exclfn is not None:
            tmpbams.append(tmpbam)
            exclfns.append(exclfn)
    ########

    print "bams:",tmpbams
    print "excl:",exclfns

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
        print "swapping reads into original and writing to", args.outBamFile
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

        print "done."

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='adds SNVs to reads, outputs modified reads as .bam along with mates')
    parser.add_argument('-v', '--varfile', dest='varFileName', required=True,
                        help='whitespace-delimited target regions to try and add a SNV: chrom,start,stop,action,seqfile if insertion,TSDlength if insertion')
    parser.add_argument('-f', '--sambamfile', dest='bamFileName', required=True,
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
    parser.add_argument('--nomut', action='store_true', default=False, help="dry run")
    parser.add_argument('--noremap', action='store_true', default=False, help="dry run")
    parser.add_argument('--noref', action='store_true', default=False, 
                        help="do not perform reference based assembly")
    parser.add_argument('--recycle', action='store_true', default=False)
    parser.add_argument('--skipmerge', action='store_true', default=False)
    args = parser.parse_args()
    main(args)


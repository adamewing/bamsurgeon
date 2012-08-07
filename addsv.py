#!/bin/env python

import re, os, sys, random
import subprocess
import collections
import asmregion
import mutableseq
import argparse
import pysam

def remap(fq1, fq2, threads, bwaref, outbam):

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

    print "mapping 1st end, cmd: " + " ".join(sai1args)
    subprocess.call(sai1args)
    print "mapping 2nd end, cmd: " + " ".join(sai2args)
    subprocess.call(sai2args)
    print "pairing ends, building .sam, cmd: " + " ".join(samargs)
    subprocess.call(samargs)
    print "sam --> bam, cmd: " + " ".join(bamargs)
    subprocess.call(bamargs)
    print "sorting, cmd: " + " ".join(sortargs)
    subprocess.call(sortargs)
    print "rename " + tmpsrt + ".bam --> " + tmpbam
    os.remove(tmpbam)
    os.rename(tmpsrt + ".bam", tmpbam)

    if os.path.isfile(outbam):
        tmpmerge  = basefn + ".merge.bam"
        mergeargs = ['samtools','merge',tmpmerge,tmpbam,outbam]
        print outbam + " exists, merging: " + " ".join(mergeargs)
        subprocess.call(mergeargs)
        os.remove(outbam)
        os.remove(tmpbam)
        print "rename " + tmpmerge + " --> " + outbam
        os.rename(tmpmerge, outbam)
    else:
        print "rename " + tmpbam + " --> " + outbam
        os.rename(tmpbam, outbam)

    # cleanup
    os.remove(sai1fn)
    os.remove(sai2fn)
    os.remove(samfn)
    os.remove(fq1)
    os.remove(fq2)

def runwgsim(contig,newseq,svfrac):
    '''
    wrapper function for wgsim
    '''
    namecount = collections.Counter(contig.reads.reads)

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

    sys.stderr.write("paired : " + str(paired) + "\n" +
                     "single : " + str(single) + "\n" +
                     "discard: " + str(discard) + "\n" +
                     "total  : " + str(totalreads) + "\n")

    # adjustment factor for length of new contig vs. old contig
    lenfrac = float(len(newseq))/float(len(contig.seq))

    sys.stderr.write("old ctg len: " + str(len(contig.seq)) + "\n" +
                     "new ctg len: " + str(len(newseq)) + "\n" +
                     "adj. factor: " + str(lenfrac) + "\n")

    # number of paried reads to simulate
    nsimreads = int((paired + (single/2)) * svfrac * lenfrac)

    sys.stderr.write("num. sim. reads: " + str(nsimreads) + "\n")

    # length of quality score comes from original read, used here to set length of read
    maxqlen = 0
    for qual in (contig.rquals + contig.mquals):
        if len(qual) > maxqlen:
            maxqlen = len(qual)

    args = ['wgsim','-e','0','-N',str(nsimreads),'-1',str(maxqlen),'-2','100','-r','0','-R','0',fasta,fq1,fq2]
    print args
    subprocess.call(args)

    os.remove(fasta)

    fqReplaceList(fq1,pairednames,contig.rquals)
    fqReplaceList(fq2,pairednames,contig.mquals)

    return (fq1,fq2)

def fqReplaceList(fqfile,names,quals):
    '''
    Replace seq names in paired fastq files from a list until the list runs out
    (then stick with original names). fqfile = fastq file, names = list
    if there are more names in the list than needed assign the remainder a null
    sequence ("NNNNNNNNNNNN...N") so that they still replace a read in the original
    .bam when reads are replaced with output 
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
                newnames.append(fqline.strip().lstrip('@').rstrip('/1').rstrip('/2'))
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
        fqout.write(seqs[i] + "\n+\n" + quals[i] + "\n")
        usednames[newnames[i]] = True

    # burn off excess
    nullseq  = 'N'*len(seqs[0])
    nullqual = '#'*len(seqs[0])
    for name in names:
        if name not in usednames:
            fqout.write("@" + name + "\n")
            fqout.write(nullseq + "\n+\n" + nullqual + "\n")

    fqout.close()

def singleseqfa(file):
    print file
    f = open(file, 'r')
    seq = ""
    for line in f:
        if not re.search ('^>',line):
            seq += line.strip()
    return seq


def main(args):
    varfile = open(args.varFileName, 'r')
    bamfile = pysam.Samfile(args.bamFileName, 'rb')
    bammate = pysam.Samfile(args.bamFileName, 'rb') # use for mates to avoid iterator problems
    reffile = pysam.Fastafile(args.refFasta)
    logfile = open(args.outBamFile + ".log", 'w')

    svfrac = float(args.svfrac)

    if os.path.isfile(args.outBamFile):
        raise ValueError(args.outBamFile + " exists, delete or rename before running.\n")

    for bedline in varfile:
        if re.search('^#',bedline):
            continue
    
        c = bedline.strip().split()
        chr    = c[0]
        start  = int(c[1])
        end    = int(c[2])
        action = c[3] # INV, DEL, INS seqfile.fa TSDlength, DUP
        insseqfile = None
        tsdlen = 0
        ndups = 0
        if action == 'INS':
            insseqfile = c[4]
            if len(c) > 5:
                tsdlen = int(c[5])

        if action == 'DUP':
            if len(c) > 4:
                ndups = int(c[4])
            else:
                ndups = 1

        contigs = asmregion.asm(chr, start, end, args.bamFileName, reffile, int(args.kmersize), args.noref, args.recycle)

        # find the largest contig        
        maxlen = 0
        maxcontig = None
        for contig in contigs:
            if contig.len > maxlen:
                maxlen = contig.len
                maxcontig = contig

        # is there anough room to make mutations?
        if maxlen > 3*int(args.maxlibsize):
            # make mutation in the largest contig
            mutseq = mutableseq.MutableSeq(maxcontig.seq)

            print "BEFORE:",mutseq

            if action == 'INS':
                mutseq.insertion(mutseq.length()/2,singleseqfa(insseqfile),tsdlen)
                logfile.write("\t".join(('ins',str(mutseq.length()/2),inseqfile,str(tsdlen))) + "\n")
            elif action == 'INV':
                invstart = int(args.maxlibsize)
                invend = mutseq.length() - invstart
                mutseq.inversion(invstart,invend)
                logfile.write("\t".join(('inv',str(invstart),str(invend))) + "\n")
            elif action == 'DEL':
                delstart = int(args.maxlibsize)
                delend = mutseq.length() - delstart
                mutseq.deletion(delstart,delend)
                logfile.write("\t".join(('del',str(delstart),str(delend))) + "\n")
            elif action == 'DUP':
                dupstart = int(args.maxlibsize)
                dupend = mutseq.length() - dupstart
                mutseq.duplication(dupstart,dupend,ndups)
                logfile.write("\t".join(('dup',str(dupstart),str(dupend),str(ndups))) + "\n")
            else:
                raise ValueError(bedline.strip() + ": mutation not one of: INS,INV,DEL,DUP")

            print "AFTER:",mutseq

            # simulate reads
            (fq1,fq2) = runwgsim(maxcontig,mutseq.seq,svfrac)

            # remap reads
            remap(fq1,fq2,4,args.refFasta,args.outBamFile)
        else:
            print "best contig too short to make mutation: ",bedline.strip()

    varfile.close()
    bamfile.close()
    bammate.close()
    logfile.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='adds SNVs to reads, outputs modified reads as .bam along with mates')
    parser.add_argument('-v', '--varfile', dest='varFileName', required=True,
                        help='whitespace-delimited target regions to try and add a SNV: chr,start,stop,action,seqfile if insertion,TSDlength if insertion')
    parser.add_argument('-f', '--sambamfile', dest='bamFileName', required=True,
                        help='sam/bam file from which to obtain reads')
    parser.add_argument('-r', '--reference', dest='refFasta', required=True,
                        help='reference genome, fasta indexed with bwa index -a stdsw _and_ samtools faidx')
    parser.add_argument('-o', '--outbam', dest='outBamFile', required=True,
                        help='.bam file name for output')
    parser.add_argument('-l', '--maxlibsize', dest='maxlibsize', default=600)
    parser.add_argument('-k', '--kmer', dest='kmersize', default=31, help="kmer size for assembly (default = 31)")
    parser.add_argument('-s', '--svfrac', dest='svfrac', default=1.0, help="allele fraction of variant (default = 1.0)")
    parser.add_argument('--nomut', action='store_true', default=False)
    parser.add_argument('--noremap', action='store_true', default=False)
    parser.add_argument('--noref', action='store_true', default=False)
    parser.add_argument('--recycle', action='store_true', default=False)
    args = parser.parse_args()
    main(args)


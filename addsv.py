#!/usr/bin/env python

import re, os, sys, random
import subprocess
import argparse
import pysam
import bs.replacereads as rr
import bs.asmregion as ar
import bs.mutableseq as ms
from collections import Counter

def remap(fq1, fq2, threads, bwaref, outbam):
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

def runwgsim(contig,newseq,svfrac,exclude):
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

    print "paired : " + str(paired)
    print "single : " + str(single)
    print "discard: " + str(discard)
    print "total  : " + str(totalreads)

    # adjustment factor for length of new contig vs. old contig
    lenfrac = float(len(newseq))/float(len(contig.seq))

    print "old ctg len: " + str(len(contig.seq))
    print "new ctg len: " + str(len(newseq))
    print "adj. factor: " + str(lenfrac)

    # number of paried reads to simulate
    nsimreads = int((paired + (single/2)) * svfrac * lenfrac)

    print "num. sim. reads: " + str(nsimreads) 

    # length of quality score comes from original read, used here to set length of read
    maxqlen = 0
    for qual in (contig.rquals + contig.mquals):
        if len(qual) > maxqlen:
            maxqlen = len(qual)

    args = ['wgsim','-e','0','-N',str(nsimreads),'-1',str(maxqlen),'-2','100','-r','0','-R','0',fasta,fq1,fq2]
    print args
    subprocess.call(args)

    os.remove(fasta)

    fqReplaceList(fq1,pairednames,contig.rquals,svfrac,exclude)
    fqReplaceList(fq2,pairednames,contig.mquals,svfrac,exclude)

    return (fq1,fq2)

def fqReplaceList(fqfile,names,quals,svfrac,exclude):
    """
    Replace seq names in paired fastq files from a list until the list runs out
    (then stick with original names). fqfile = fastq file, names = list
    if there are more names in the list than needed assign the remainder a null
    sequence ("NNNNNNNNNNNN...N") so that they still replace a read in the original
    .bam when reads are replaced with output

    'exclude' is a filehandle

    """
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
        fqout.write(seqs[i] + "\n+\n" + quals[i] + "\n")
        if newnames[i] in usednames:
            print "warning, used read name: " + newnames[i] + " in multiple pairs"
        usednames[newnames[i]] = True
        
    # burn off excess
    nullseq  = 'N'*len(seqs[0])
    nullqual = '#'*len(seqs[0])
    for name in names:
        if name not in usednames:
            if random.uniform(0,1) < svfrac:
                fqout.write("@" + name + "\n")
                fqout.write(nullseq + "\n+\n" + nullqual + "\n")
                exclude.write(name + "\n")

    fqout.close()

def singleseqfa(file):
    print file
    f = open(file, 'r')
    seq = ""
    for line in f:
        if not re.search ('^>',line):
            seq += line.strip().upper()
    return seq

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

def main(args):
    """ needs refactoring
    """
    varfile = open(args.varFileName, 'r')
    bamfile = pysam.Samfile(args.bamFileName, 'rb')
    reffile = pysam.Fastafile(args.refFasta)
    logfile = open(args.outBamFile + ".log", 'w')
    exclude = open(args.exclfile, 'w')

    # temporary file to hold mutated reads
    outbam_mutsfile = "tmp." + str(random.random()) + ".muts.bam"

    svfrac = float(args.svfrac)

    nmuts = 0

    for bedline in varfile:
        if re.search('^#',bedline):
            continue
   
        if args.maxmuts and nmuts >= int(args.maxmuts):
            break
 
        c = bedline.strip().split()
        chr    = c[0]
        start  = int(c[1])
        end    = int(c[2])
        araw   = c[3:len(c)] # INV, DEL, INS seqfile.fa TSDlength, DUP
        actions = map(lambda x: x.strip(),' '.join(araw).split(','))

        print "interval:",c
        # modify start and end if interval is too long
        maxctglen = int(args.maxctglen)
        assert maxctglen > 3*int(args.maxlibsize) # maxctglen is too short
        if end-start > maxctglen:
            adj   = (end-start) - maxctglen
            rndpt = random.randint(0,adj)
            start = start + rndpt
            end   = end - (adj-rndpt)
            print "note: interval size too long, adjusted:",chr,start,end

        contigs = ar.asm(chr, start, end, args.bamFileName, reffile, int(args.kmersize), args.noref, args.recycle)

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
            mutseq = ms.MutableSeq(maxcontig.seq)

            # if we're this far along, we're making a mutation
            nmuts += 1 

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

                print "BEFORE:",mutseq

                if action == 'INS':
                    if insseqfile: # seq in file
                        mutseq.insertion(mutseq.length()/2,singleseqfa(insseqfile),tsdlen)
                    else: # seq is input
                        mutseq.insertion(mutseq.length()/2,insseq,tsdlen)
                    logfile.write("\t".join(('ins',chr,str(start),str(end),action,str(mutseq.length()),str(mutseq.length()/2),str(insseqfile),str(tsdlen))) + "\n")

                elif action == 'INV':
                    invstart = int(args.maxlibsize)
                    invend = mutseq.length() - invstart
                    mutseq.inversion(invstart,invend)
                    logfile.write("\t".join(('inv',chr,str(start),str(end),action,str(mutseq.length()),str(invstart),str(invend))) + "\n")

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
                    logfile.write("\t".join(('del',chr,str(start),str(end),action,str(mutseq.length()),str(delstart),str(delend),str(dlen))) + "\n")

                elif action == 'DUP':
                    dupstart = int(args.maxlibsize)
                    dupend = mutseq.length() - dupstart
                    mutseq.duplication(dupstart,dupend,ndups)
                    logfile.write("\t".join(('dup',chr,str(start),str(end),action,str(mutseq.length()),str(dupstart),str(dupend),str(ndups))) + "\n")

                else:
                    raise ValueError(bedline.strip() + ": mutation not one of: INS,INV,DEL,DUP")

                print "AFTER:",mutseq

            # simulate reads
            (fq1, fq2) = runwgsim(maxcontig, mutseq.seq, svfrac, exclude)

            # remap reads
            remap(fq1, fq2, 4, args.refFasta, outbam_mutsfile)

        else:
            print "best contig too short to make mutation: ",bedline.strip()

    print "addsv.py finished, made", nmuts, "mutations"

    exclude.close()
    varfile.close()
    bamfile.close()
    logfile.close()

    print "merging mutations into", args.bamFileName, "-->", args.outBamFile
    replace(args.bamFileName, outbam_mutsfile, args.outBamFile, args.exclfile)

    # cleanup
    os.remove(outbam_mutsfile)

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
    parser.add_argument('-l', '--maxlibsize', dest='maxlibsize', default=600, help="maximum fragment length of seq. library")
    parser.add_argument('-k', '--kmer', dest='kmersize', default=31, 
                        help="kmer size for assembly (default = 31)")
    parser.add_argument('-s', '--svfrac', dest='svfrac', default=1.0, 
                        help="allele fraction of variant (default = 1.0)")
    parser.add_argument('-x', '--excluded', dest='exclfile', default='excluded.txt',
                        help="output excluded (e.g. from a deletion) read names to file (default=excluded.txt)")
    parser.add_argument('--maxctglen', dest='maxctglen', default=32000, 
                        help="maximum contig length for assembly - can increase if velvet is compiled with LONGSEQUENCES")
    parser.add_argument('-n', dest='maxmuts', default=None,
                        help="maximum number of mutations to make")
    parser.add_argument('--nomut', action='store_true', default=False, help="dry run")
    parser.add_argument('--noremap', action='store_true', default=False, help="dry run")
    parser.add_argument('--noref', action='store_true', default=False, 
                        help="do not perform reference based assembly")
    parser.add_argument('--recycle', action='store_true', default=False)
    args = parser.parse_args()
    main(args)


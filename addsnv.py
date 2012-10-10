#!/usr/bin/env python

import sys, pysam, argparse, collections, random, subprocess, os

def majorbase(basepile):
    """returns tuple: (major base, count)
    """
    return collections.Counter(basepile).most_common()[0]

def minorbase(basepile):
    """returns tuple: (minor base, count)
    """
    c = collections.Counter(basepile)
    if len(list(c.elements())) > 1:
        return c.most_common(2)[-1]
    else:
        return c.most_common()[0]

def mut(base,det=False):
    """ change base to something different
        if 'det' (deterministic) is true, mutations will be made in a predictable pattern:
        A-->G, G-->A, T-->C, C-->T (transitions)
    """

    bases = ('A','T','C','G')
    base = base.upper()
    if base not in bases:
        raise ValueError("base passed to mut(): " + str(base) + " not one of (A,T,C,G)")

    if det:
        if base == 'A':
            return 'T'
        elif base == 'T':
            return 'A'
        elif base == 'G':
            return 'C'
        elif base == 'C':
            return 'G'

    else:
        mut = base
        while mut == base:
            mut = bases[int(random.uniform(0,4))]
        return mut

def countReadCoverage(bam,chrom,start,end,strand=None):
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
                    if strand == '+':
                        if not read.alignment.is_reverse and read.alignment.mapq >= 0 and not read.alignment.is_duplicate:
                            n += 1
                    elif strand == '-':
                        if read.alignment.is_reverse and read.alignment.mapq >= 0 and not read.alignment.is_duplicate:
                            n += 1
                    else:
                        if read.alignment.mapq >= 0 and not read.alignment.is_duplicate:
                            n += 1
                coverage[i] = n
                i += 1

    return coverage


def countBaseAtPos(bamfile,chrom,pos):
    """ return list of bases at position chrom,pos
    """
    locstr = chrom + ":" + str(pos) + "-" + str(pos)
    args = ['samtools','mpileup',bamfile,'-r',locstr]

    p = subprocess.Popen(args,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    p.wait()
    pout = p.stdout.readlines()

    pileup = None 

    for line in pout:
        c = line.strip().split()
        pileup = c[4].upper()

    bases = []
    if pileup:
        for b in pileup:
            if b in ['A','T','C','G']:
                bases.append(b)

    return bases

def mergebams(bamlist,outbamfn):
    """ call samtools to merge two .bams
    """
    args = ['samtools','merge','-f',outbamfn] + bamlist
    print "merging, cmd: ",args
    subprocess.call(args)

    for bamfile in bamlist:
        os.remove(bamfile)
        os.remove(bamfile + '.bai')

def remap(bamfn, threads, bwaref):
    """ call bwa/samtools to remap .bam
    """
    sai1fn = bamfn + ".1.sai"
    sai2fn = bamfn + ".2.sai"
    samfn  = bamfn + ".sam"
    refidx = bwaref + ".fai"

    sai1args = ['bwa', 'aln', bwaref, '-q', '5', '-l', '32', '-k', '3', '-t', str(threads), '-o', '1', '-f', sai1fn, '-b1', bamfn]
    sai2args = ['bwa', 'aln', bwaref, '-q', '5', '-l', '32', '-k', '3', '-t', str(threads), '-o', '1', '-f', sai2fn, '-b2', bamfn]
    samargs  = ['bwa', 'sampe', '-P', '-f', samfn, bwaref, sai1fn, sai2fn, bamfn, bamfn]
    bamargs  = ['samtools', 'view', '-bt', refidx, '-o', bamfn, samfn] 

    print "mapping 1st end, cmd: " + " ".join(sai1args)
    subprocess.call(sai1args)
    print "mapping 2nd end, cmd: " + " ".join(sai2args)
    subprocess.call(sai2args)
    print "pairing ends, building .sam, cmd: " + " ".join(samargs)
    subprocess.call(samargs)
    print "sam --> bam, cmd: " + " ".join(bamargs)
    subprocess.call(bamargs)

    sortbase = bamfn + ".sort"
    sortfn   = sortbase + ".bam"
    sortargs = ['samtools','sort','-m','10000000000',bamfn,sortbase]
    print "sorting, cmd: " + " ".join(sortargs)
    subprocess.call(sortargs)
    os.rename(sortfn,bamfn)

    indexargs = ['samtools','index',bamfn]
    print "indexing, cmd: " + " ".join(indexargs)
    subprocess.call(indexargs)

    # cleanup
    os.remove(sai1fn)
    os.remove(sai2fn)
    os.remove(samfn)

def main(args):
    """ needs refactoring
    """
    bedfile = open(args.varFileName, 'r')
    bamfile = pysam.Samfile(args.bamFileName, 'rb')
    bammate = pysam.Samfile(args.bamFileName, 'rb') # use for mates to avoid iterator problems
    reffile = pysam.Fastafile(args.refFasta)
    outbam  = pysam.Samfile(args.outBamFile, 'wb', template=bamfile)
    outbam.close()
    tmpbams = []

    log = open(args.outBamFile + ".log",'w')

    snvfrac = float(args.snvfrac)

    for bedline in bedfile:
        if len(tmpbams) < int(args.numsnvs) or int(args.numsnvs) == 0:
            c = bedline.strip().split()
            chrom   = c[0]
            start = int(c[1])
            end   = int(c[2])

            gmutpos = int(random.uniform(start,end+1)) # position of mutation in genome
            refbase = reffile.fetch(chrom,gmutpos-1,gmutpos)
            mutbase = mut(refbase,args.det)
            mutstr = refbase + "-->" + mutbase

            # keep a list of reads to modify - use hash to keep unique since each
            # read will be visited as many times as it has bases covering the region
            outreads = {}
            mutreads = {} # same keys as outreads
            mutmates = {} # same keys as outreads, keep track of mates
            numunmap = 0
            hasSNP = False
            tmpoutbamname = "tmpbam" + str(random.random()) + ".bam"
            print "creating tmp bam: ",tmpoutbamname #DEBUG
            outbam = pysam.Samfile(tmpoutbamname, 'wb', template=bamfile)
            maxfrac = 0.0

            for pcol in bamfile.pileup(reference=chrom,start=gmutpos,end=gmutpos+1):
                # this will include all positions covered by a read that covers the region of interest
                if pcol.pos: #> start and pcol.pos <= end:
                    refbase = reffile.fetch(chrom,pcol.pos-1,pcol.pos)
                    basepile = ''
                    for pread in pcol.pileups:
                        basepile += pread.alignment.seq[pread.qpos-1]
                        pairname = 'F' # read is first in pair
                        if pread.alignment.is_read2:
                            pairname = 'S' # read is second in pair
                        if not pread.alignment.is_paired:
                            pairname = 'U' # read is unpaired

                        extqname = ','.join((pread.alignment.qname,str(pread.alignment.pos),pairname))

                        if pcol.pos == gmutpos:
                            if not pread.alignment.mate_is_unmapped:
                                outreads[extqname] = pread.alignment
                                mutbases = list(pread.alignment.seq)
                                mutbases[pread.qpos-1] = mutbase
                                mutread = ''.join(mutbases)
                                mutreads[extqname] = mutread
                                mate = None
                                try:
                                    mate = bammate.mate(pread.alignment)
                                except:
                                    print "warning: no mate for",pread.alignment.qname
                                mutmates[extqname] = mate
                                log.write(" ".join(('read',extqname,mutread,"\n")))
                            else:
                                numunmap += 1

                    # make sure region doesn't have any changes that are likely SNPs
                    # (trying to avoid messing with haplotypes)
                
                    basepile = countBaseAtPos(args.bamFileName,chrom,pcol.pos)
                    if basepile:
                        majb = majorbase(basepile)
                        minb = minorbase(basepile)

                        frac = float(minb[1])/(float(majb[1])+float(minb[1]))
                        if minb[0] == majb[0]:
                            frac = 0.0
                        if frac > maxfrac:
                            maxfrac = frac
                        if frac > snvfrac:
                            hasSNP = True
                    else:
                        hasSNP = True

            # pick reads to change
            readlist = []
            for extqname,read in outreads.iteritems():
                if read.seq != mutreads[extqname]:
                    readlist.append(extqname)

            print "len(readlist):",str(len(readlist))
            random.shuffle(readlist)
            readlist = readlist[0:int(len(readlist)*float(args.mutfrac))] 
            print "picked:",str(len(readlist))

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
                    outbam.write(read)
                    if mutmates[extqname]:
                        outbam.write(mutmates[extqname])
            print "wrote: ",wrote,"mutated:",nmut

            if not hasSNP or args.force:
                outbam.close()
                remap(tmpoutbamname, 4, args.refFasta)

                outbam = pysam.Samfile(tmpoutbamname,'rb')
                coverwindow = 1
                incover  = countReadCoverage(bamfile,chrom,gmutpos-coverwindow,gmutpos+coverwindow)
                outcover = countReadCoverage(outbam,chrom,gmutpos-coverwindow,gmutpos+coverwindow)

                avgincover  = float(sum(incover))/float(len(incover)) 
                avgoutcover = float(sum(outcover))/float(len(outcover))
                snvfrac = 0.0
                if wrote > 0:
                    snvfrac = float(nmut)/float(wrote)

                # qc cutoff for final snv depth 
                if (avgoutcover > 0 and avgincover > 0 and avgoutcover/avgincover >= 0.9) or args.force:
                    tmpbams.append(tmpoutbamname)
                    log.write("\t".join(("snv",bedline.strip(),str(gmutpos),mutstr,str(avgoutcover),str(avgoutcover),str(snvfrac),str(maxfrac)))+"\n")

            outbam.close()

    # merge tmp bams
    if len(tmpbams) == 1:
        os.rename(tmpbams[0],args.outBamFile)
    elif len(tmpbams) > 1:
        mergebams(tmpbams,args.outBamFile)

    bedfile.close()
    bamfile.close()
    bammate.close()
    log.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='adds SNVs to reads, outputs modified reads as .bam along with mates')
    parser.add_argument('-v', '--varfile', dest='varFileName', required=True,
                        help='Target regions to try and add a SNV, as BED')
    parser.add_argument('-f', '--sambamfile', dest='bamFileName', required=True,
                        help='sam/bam file from which to obtain reads')
    parser.add_argument('-r', '--reference', dest='refFasta', required=True,
                        help='reference genome, fasta indexed with bwa index -a stdsw _and_ samtools faidx')
    parser.add_argument('-o', '--outbam', dest='outBamFile', required=True,
                        help='.bam file name for output')
    parser.add_argument('-s', '--snvfrac', dest='snvfrac', default=1, 
                        help='maximum allowable linked SNP MAF (for avoiding haplotypes) (default = 1)')
    parser.add_argument('-m', '--mutfrac', dest='mutfrac', default=0.5, 
                        help='allelic fraction at which to make SNVs (default = 0.5)')
    parser.add_argument('-n', '--numsnvs', dest='numsnvs', default=0, 
                        help="maximum number of mutations to make (default: entire input)")
    parser.add_argument('--nomut', action='store_true', default=False, help="dry run")
    parser.add_argument('--det', action='store_true', default=False, help="deterministic base changes: make transitions only")
    parser.add_argument('--force', action='store_true', default=False, help="force mutation to happen regardless of nearby SNP or low coverage")
    args = parser.parse_args()
    main(args)

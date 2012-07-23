#!/bin/env python

import sys, pysam, argparse, collections, random, subprocess, os

def majorbase(basepile):
    """returns tuple: (major base, count)"""
    return collections.Counter(basepile).most_common()[0]

def minorbase(basepile):
    """returns tuple: (minor base, count)"""
    c = collections.Counter(basepile)
    if len(list(c.elements())) > 1:
        return c.most_common(2)[-1]
    else:
        return c.most_common()[0]

def mut(base):
    """change base to something different"""
    bases = ('A','T','C','G')
    mut = base
    while mut == base:
        mut = bases[int(random.uniform(0,4))]
    return mut

def remap(bamfn, threads, bwaref):
    sai1fn = bamfn + ".1.sai"
    sai2fn = bamfn + ".2.sai"
    samfn  = bamfn + ".sam"
    refidx = bwaref + ".fai"

    sai1args = ['bwa', 'aln', bwaref, '-q', '5', '-l', '32', '-k', '2', '-t', str(threads), '-o', '1', '-f', sai1fn, '-b1', bamfn]
    sai2args = ['bwa', 'aln', bwaref, '-q', '5', '-l', '32', '-k', '2', '-t', str(threads), '-o', '1', '-f', sai2fn, '-b2', bamfn]
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

    # cleanup
    os.remove(sai1fn)
    os.remove(sai2fn)
    os.remove(samfn)

def main(args):
    bedfile = open(args.bedFileName, 'r')
    bamfile = pysam.Samfile(args.bamFileName, 'rb')
    bammate = pysam.Samfile(args.bamFileName, 'rb') # use for mates to avoid iterator problems
    reffile = pysam.Fastafile(args.refFasta)
    outbam  = pysam.Samfile(args.outBamFile, 'wb', template=bamfile)

    snpfrac = float(args.snpfrac)

    for bedline in bedfile:
        c = bedline.strip().split()
        chr   = c[0]
        start = int(c[1])
        end   = int(c[2])

        gmutpos = int(random.uniform(start,end+1)) # position of mutation in genome
        refbase = reffile.fetch(chr,gmutpos-1,gmutpos)
        mutbase = mut(refbase)
        mutstr = refbase + "-->" + mutbase

        print "\t".join((bedline.strip(),str(gmutpos),mutstr)) # debug

        # keep a list of reads to modify - use hash to keep unique since each
        # read will be visited as many times as it has bases covering the region
        outreads = {}
        mutreads = {} # same keys as outreads
        mutmates = {} # same keys as outreads, keep track of mates
        numunmap = 0
        hasSNP = False
        for pcol in bamfile.pileup(reference=chr,start=gmutpos,end=gmutpos+1):
            # this will include all positions covered by a read that covers the region of interest
            if pcol.pos: #> start and pcol.pos <= end:
                refbase = reffile.fetch(chr,pcol.pos-1,pcol.pos)
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
                            mate = bammate.mate(pread.alignment)
                            mutmates[extqname] = mate
                        else:
                            numunmap += 1

                # make sure region doesn't have any changes that are likely SNPs
                # (trying to avoid messing with haplotypes)
                majb = majorbase(basepile)
                minb = minorbase(basepile)
                frac = float(minb[1])/(float(majb[1])+float(minb[1]))
                if minb[0] == majb[0]:
                    frac = 0.0
                if frac > snpfrac:
                    hasSNP = True
                print " ".join((refbase,basepile,str(pcol.pos),str(majb),str(minb),str(hasSNP),str(frac))) #debug

        '''
        # pick reads to change
        readlist = []
        for extqname,read in outreads.iteritems():
            if read.seq != mutreads[extqname]:
                readlist.append(
        '''

        # change reads from .bam to mutated sequences
        for extqname,read in outreads.iteritems():
            if read.seq != mutreads[extqname]:
                #print read.seq
                #print mutreads[extqname]
                mutprob = random.uniform(0,1) # choose reads to mutate at random
                if not args.nomut and mutprob < float(args.mutfrac):
                    read.seq = mutreads[extqname] # make mutation
            if not hasSNP:
                outbam.write(read)
                outbam.write(mutmates[extqname])

    bedfile.close()
    bamfile.close()
    bammate.close()
    outbam.close()

    if not args.nomut and not args.noremap:
        remap(args.outBamFile, 4, args.refFasta)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='adds SNVs to reads, outputs modified reads as .bam along with mates')
    parser.add_argument('-b', '--bedfile', dest='bedFileName', required=True,
                        help='Target regions to try and add a SNV, as BED')
    parser.add_argument('-f', '--sambamfile', dest='bamFileName', required=True,
                        help='sam/bam file from which to obtain reads')
    parser.add_argument('-r', '--reference', dest='refFasta', required=True,
                        help='reference genome, fasta indexed with bwa index -a stdsw _and_ samtools faidx')
    parser.add_argument('-o', '--outbam', dest='outBamFile', required=True,
                        help='.bam file name for output')
    parser.add_argument('-s', '--snpfrac', dest='snpfrac', default=1)
    parser.add_argument('-m', '--mutfrac', dest='mutfrac', default=0.5)
    parser.add_argument('--nomut', action='store_true', default=False)
    parser.add_argument('--noremap', action='store_true', default=False)
    args = parser.parse_args()
    main(args)

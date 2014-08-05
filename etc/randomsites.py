#!/usr/bin/env python

import argparse
import random
import pysam
import re
import subprocess

def main(args):

    genome = pysam.Fastafile(args.fastaFile)

    maptabix = None
    if args.maptabix:
        maptabix = pysam.Tabixfile(args.maptabix)

    minlen = int(args.minlen)
    maxlen = int(args.maxlen)

    assert minlen <= maxlen

    mask = None
    if args.mask is not None:
        mask = pysam.Tabixfile(args.mask)

    chrlen = {}
    for line in open(args.fastaFile + '.fai', 'r'):
        c = line.strip().split()
        chr = c[0]
        size = int(c[1])
        chrlen[chr] = size

    # calculate offsets
    offset = 0
    chroffset = {}
    for chrom in sorted(chrlen.iterkeys()):
        offset += int(chrlen[chrom])
        chroffset[chrom] = offset

    genomelen = offset

    # random picks
    n = 0
    while n < int(args.numpicks):
        rndloc = random.randint(0,genomelen)
        lastoffset = 0
        rndchr = None
        for chrom in sorted(chrlen.iterkeys()):
            offset = int(chroffset[chrom])
            assert lastoffset < offset
            if rndloc >= lastoffset and rndloc < offset:
                rndchr = chrom
                rndloc -= lastoffset
            lastoffset = offset

        # handle single chromosome option
        if args.chrom and args.chrom != rndchr:
            continue

        fraglen = int(random.uniform(minlen,maxlen))
        fragstart = rndloc
        fragend   = rndloc + int(fraglen)

        # handle mask
        if mask is not None and rndchr in mask.contigs and len(list(mask.fetch(rndchr, fragstart, fragend))) > 0:
            continue

        # handle pmin/pmax args
        if args.minpos is not None and fragstart < int(args.minpos):
            continue

        if args.maxpos is not None and fragend > int(args.maxpos):
            continue

        # handle mappability option
        if maptabix:
            reject = False 

            if rndchr not in maptabix.contigs and 'chr' + rndchr in maptabix.contigs:
                mchrom = 'chr' + rndchr
            else:
                mchrom = rndchr 

            if mchrom in maptabix.contigs:
                for mapline in maptabix.fetch(mchrom, fragstart, fragend):
                    m = mapline.strip().split()
                    mscore = float(m[3])

                    if mscore < float(args.minmap):
                        reject = True
            if reject:
                continue

        # handle bamfile depth option
        if args.bamfile is not None and int(args.mindepth) > 0:
            reject = False
    
            # pileup 
            region = rndchr + ':' + str(fragstart) + '-' + str(fragend)
            mpargs = ['samtools', 'mpileup', '-r', region, args.bamfile]
            p = subprocess.Popen(mpargs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            lastpos = fragstart-1
            for line in p.stdout.readlines():
                c = line.strip().split()
                if len(c) == 6:
                    pos   = int(c[1])
                    depth = int(c[3])
                    if depth < int(args.mindepth):
                        reject = True
                    if pos - lastpos > 1:
                        reject = True
                    lastpos = pos

            if reject:
                continue

        # don't pick sites on contigs if --nocontigs is set
        if (rndchr.startswith('chr') and len(rndchr) > 5) or (not rndchr.startswith('chr') and len(rndchr) > 2):
            if args.nocontigs:
                continue

        if fragend < offset:
            if genome:
                seq = genome.fetch(rndchr,fragstart,fragend)
                assert seq

                if args.requireseq:
                    if re.search('[ATGCatgc]',seq):
                        print "\t".join((rndchr,str(fragstart),str(fragend),seq))
                        n += 1
                else:
                    print "\t".join((rndchr,str(fragstart),str(fragend),seq))
                    n += 1
            else:
                print "\t".join((rndchr,str(fragstart),str(fragend)))
                n += 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="pick random sites from a samtools-indexed genome")
    parser.add_argument('-f', '--fasta', dest='fastaFile', required=True, help="reference fasta, must be indexed with samtools faidx")
    parser.add_argument('-n', '--num', dest='numpicks', required=True, help="number of sites to pick")
    parser.add_argument('-c', '--chr', dest='chrom', default=None, help='make mutations on one chromosome only')
    parser.add_argument('-b', '--bamfile', dest='bamfile', default=None, help='bamfile to use with -d/--mindepth')
    parser.add_argument('-d', '--mindepth', dest='mindepth', default=0, help='mindepth to use with -b/--bamfile')
    parser.add_argument('--maptabix', dest='maptabix', default=None, help='mappability tabix, required for --minmap')
    parser.add_argument('--minmap', dest='minmap', default=0.8, help='only select regions above mappability threshold (default 0.8)')
    parser.add_argument('--lmin', dest='minlen', default=1, help='minimum fragment length (default=1)')
    parser.add_argument('--lmax', dest='maxlen', default=1, help='maximum fragment length (default=1)')
    parser.add_argument('--pmin', dest='minpos', default=None, help='minimum position')
    parser.add_argument('--pmax', dest='maxpos', default=None, help='maximum position')
    parser.add_argument('--mask', dest='mask', default=None, help='mask (tabix indexed BED-3)')
    parser.add_argument('--requireseq', action="store_true", help="do not select hits in unsequenced regions, requires fasta file")
    parser.add_argument('--nocontigs', action="store_true", help="exclude contigs")
    args = parser.parse_args()
    main(args)

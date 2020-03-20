#!/usr/bin/env python

import os
import random
import pysam
import argparse

'''
Tool for generating random sites for bamsurgeon input

'''

class Genome:
    def __init__(self, gfn, bedfile=None):
        ''' gfn = genome file name '''
        self.chrlen = {} # length of each chromosome
        self.chrmap = [] # used for picking chromosomes
        self.ref    = pysam.Fastafile(gfn)

        bp = 0
        bins = 100000

        with open(gfn + '.fai', 'r') as g:
            for line in g:
                chrom, length = line.strip().split()[:2]
                self.chrlen[chrom] = int(length)
                bp += int(length)

        for chrom, length in self.chrlen.items():
            self.chrmap += [chrom] * int(float(length) / float(bp) * bins)

        self.bed = []
        self.bedfile = None

        if bedfile is not None:
            self.bedfile = bedfile
            with open(self.bedfile, 'r') as bed:
                for line in bed:
                    chrom, start, end = line.strip().split()[:3]
                    start = int(start)
                    end   = int(end)
                    assert start < end, "start >= end in BED line: " + line.strip()
                    self.bed.append((chrom, start, end))


    def pick(self, mutlen, avoidN=False, usebed=False):
        ''' return a random chromosome and position '''
        goodmut = False

        if usebed:
            seg = random.choice(self.bed)
            rpos = int(random.uniform(seg[1], seg[2]))
            return seg[0], rpos, rpos + mutlen

        else:
            goodmut = False
            while not goodmut:
                rchrom = random.choice(self.chrmap)
                rpos   = int(random.uniform(1, self.chrlen[rchrom]))

                if avoidN and 'N' not in self.ref.fetch(rchrom, rpos-1, rpos+mutlen):
                    goodmut = True

                if not avoidN:
                    goodmut = True

            return rchrom, rpos, rpos + mutlen


def randomseq(len):
    ''' make a random DNA sequence '''
    return ''.join([random.choice(['A','T','G','C']) for _ in range(int(len))])


def randomsv():
    ''' random SV information '''
    i = random.randint(0,4)
    if i == 0: # DEL
        dfrac = random.uniform(0.5,1)
        return 'DEL ' + str(dfrac)
    if i == 1: # INS
        tsdlen = random.randint(0,30)
        return 'INS RND ' + str(tsdlen)
    if i == 2: # INV
        return 'INV'
    if i == 3: # DUP
        ndups = random.randint(1,4)
        return 'DUP ' + str(ndups)
    if i == 4: # TRN
        return 'TRN'


def betafunc(a,b):
    ''' return appropriate function from random class '''
    return lambda : random.betavariate(float(a), float(b))


def scalefunc(min, max):
    ''' return func for rescaling value from range 0,1 to min,max '''
    min, max = float(min), float(max)
    return lambda x : float(x) * (max-min) + min


def run_snv(g, args):
    ''' generate input for addsnv.py '''
    vaf = betafunc(args.vafbeta1, args.vafbeta2)
    vafscale = scalefunc(args.minvaf, args.maxvaf)

    usebed = args.bed is not None

    for _ in range(int(args.numpicks)):
        rchrom, rstart, rend = g.pick(0, avoidN=args.avoidN, usebed=usebed)
        info = [rchrom, rstart, rend, vafscale(vaf())]
        print('\t'.join(map(str, info)))


def run_indel(g, args):
    ''' generate input for addindel.py '''
    vaf = betafunc(args.vafbeta1, args.vafbeta2)
    len = betafunc(args.lenbeta1, args.lenbeta2)

    vafscale = scalefunc(args.minvaf, args.maxvaf)
    lenscale = scalefunc(args.minlen, args.maxlen)

    usebed = args.bed is not None

    for _ in range(int(args.numpicks)):
        mutlen = int(lenscale(len()))
        rchrom, rstart, rend = g.pick(mutlen, avoidN=args.avoidN, usebed=usebed)
        if random.uniform(0,1) < 0.5: # deletion
            info = [rchrom, rstart, rend, vafscale(vaf()), 'DEL']
        else: # insertion
            info = [rchrom, rstart, rend, vafscale(vaf()), 'INS', randomseq(mutlen)]

        print('\t'.join(map(str, info)))


def run_sv(g, args):
    ''' generate input for addsv.py '''
    vaf = betafunc(args.vafbeta1, args.vafbeta2)
    len = betafunc(args.lenbeta1, args.lenbeta2)

    vafscale = scalefunc(args.minvaf, args.maxvaf)
    lenscale = scalefunc(args.minlen, args.maxlen)

    usebed = args.bed is not None

    with open(args.cnvfile, 'w') as cnv:
        for _ in range(int(args.numpicks)):
            mutlen = int(lenscale(len()))

            rchrom, rstart, rend = g.pick(mutlen, avoidN=args.avoidN, usebed=usebed)

            info = [rchrom, rstart, rend, randomsv()]

            if info[-1] == 'TRN':
                tsd_partner = list(g.pick(mutlen, avoidN=args.avoidN, usebed=usebed))
                info = info[:3] + ['TRN'] + tsd_partner
                partner_cnv = tsd_partner + [1.0/(vafscale(vaf()))]

            cnvinfo = [rchrom, rstart, rend, 1.0/(vafscale(vaf()))]

            print('\t'.join(map(str, info)))

            cnv.write('\t'.join(map(str, cnvinfo)) + '\n')

def main(args):
    if args.seed is not None:
        random.seed(int(args.seed))

    assert args.minvaf < args.maxvaf, "--minvaf( = %s) is higher than --maxvaf ( = %s)" % (args.minvaf, args.maxvaf)
    assert os.path.exists(args.genome + '.fai'), "reference FASTA not indexed: " + args.genome

    g = Genome(args.genome, bedfile=args.bed)

    args.func(g, args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='make random sites for bamsurgeon')
    parser.add_argument('-g', '--genome', required=True, help='Genome FASTA, indexed with samtools faidx (expects .fai file exists)')
    parser.add_argument('-b', '--bed', default=None, help='use segments in specified BED file to choose mutations')
    parser.add_argument('-s', '--seed', default=None, help='use a seed to make reproducible random picks')
    parser.add_argument('-n', '--numpicks', default=1000, help='number of sites to generate (default = 1000)')
    parser.add_argument('--avoidN', default=False, action='store_true', help='avoid picking sites with N characters in ref')
    parser.add_argument('--minvaf', default=0.25, help='minimum variant allele fraction (default = 0.25)')
    parser.add_argument('--maxvaf', default=0.5, help='maximum variant allele fraction (default = 0.5)')
    parser.add_argument('--vafbeta1', default=2.0, help='left shape parameter for beta distribution of VAFs (default = 2.0)')
    parser.add_argument('--vafbeta2', default=2.0, help='right shape parameter for beta distribution of VAFs (default = 2.0)')
    subparsers = parser.add_subparsers(title="mode")

    parser_snv = subparsers.add_parser('snv')
    parser_snv.set_defaults(func=run_snv)

    parser_indel = subparsers.add_parser('indel')
    parser_indel.add_argument('--minlen', default=1,  help='minimum SV contig length (default = 1)')
    parser_indel.add_argument('--maxlen', default=90, help='maximum SV contig length (default = 90)')
    parser_indel.add_argument('--lenbeta1', default=0.5, help='left shape parameter for beta dist. of indel lengths (default = 0.5)')
    parser_indel.add_argument('--lenbeta2', default=4.0, help='right shape parameter for beta dist. of indel lengths (default = 4.0)')
    parser_indel.set_defaults(func=run_indel)

    parser_sv = subparsers.add_parser('sv')
    parser_sv.add_argument('--minlen', default=3000,  help='minimum SV contig length (default = 3000)')
    parser_sv.add_argument('--maxlen', default=30000, help='maximum SV contig length (default = 30000)')
    parser_sv.add_argument('--lenbeta1', default=1.0, help='left shape parameter for beta dist. of indel lengths (default = 1.0)')
    parser_sv.add_argument('--lenbeta2', default=1.0, help='right shape parameter for beta dist. of indel lengths (default = 1.0)')
    parser_sv.add_argument('--cnvfile', default='cnvs.txt', help='output file for CNV information (used for SV VAF)')
    parser_sv.set_defaults(func=run_sv)

    args = parser.parse_args()
    main(args)

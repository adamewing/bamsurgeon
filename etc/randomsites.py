#!/usr/bin/env python

import random
import argparse

'''
http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
'''

class Genome:
    def __init__(self, gfn):
        ''' gfn = genome file name '''
        self.chrlen = {} # length of each chromosome
        self.chrmap = [] # used for picking chromosomes

        bp = 0
        bins = 100000

        with open(gfn, 'r') as g:
            for line in g:
                chrom, len = line.strip().split()
                self.chrlen[chrom] = int(len)
                bp += int(len)

        for chrom, len in self.chrlen.iteritems():
            self.chrmap += [chrom] * int(float(len) / float(bp) * bins)

    def pick(self):
        ''' return a random chromosome and position '''
        rchrom = random.choice(self.chrmap)
        rpos   = int(random.uniform(1, self.chrlen[rchrom]))

        return rchrom, rpos


def randomseq(len):
    ''' make a random DNA sequence '''
    return ''.join([random.choice(['A','T','G','C']) for _ in range(int(len))])
    

def randomsv():
    ''' random SV information '''
    i = random.randint(0,3)
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

    for _ in range(int(args.numpicks)):
        rchrom, rstart = g.pick()
        info = [rchrom, rstart, rstart, vafscale(vaf())]
        print '\t'.join(map(str, info))


def run_indel(g, args):
    ''' generate input for addindel.py '''
    vaf = betafunc(args.vafbeta1, args.vafbeta2)
    len = betafunc(args.lenbeta1, args.lenbeta2)

    vafscale = scalefunc(args.minvaf, args.maxvaf)
    lenscale = scalefunc(args.minlen, args.maxlen)

    for _ in range(int(args.numpicks)):
        rchrom, rstart = g.pick()
        if random.uniform(0,1) < 0.5: # deletion
            info = [rchrom, rstart, rstart + int(lenscale(len())), vafscale(vaf()), 'DEL']
        else: # insertion
            info = [rchrom, rstart, rstart, vafscale(vaf()), 'INS', randomseq(int(lenscale(len())))]

        print '\t'.join(map(str, info))


def run_sv(g, args):
    ''' generate input for addsv.py '''
    vaf = betafunc(args.vafbeta1, args.vafbeta2)
    len = betafunc(args.lenbeta1, args.lenbeta2)

    vafscale = scalefunc(args.minvaf, args.maxvaf)
    lenscale = scalefunc(args.minlen, args.maxlen)

    with open(args.cnvfile, 'w') as cnv:
        for _ in range(int(args.numpicks)):
            rchrom, rstart = g.pick()
            rend = rstart + int(lenscale(len()))
            info = [rchrom, rstart, rend, randomsv()]
            cnvinfo = [rchrom, rstart, rend, 1.0/(vafscale(vaf()))]
            print '\t'.join(map(str, info))
            cnv.write('\t'.join(map(str, cnvinfo)) + '\n')

def main(args):
    if args.seed is not None:
        random.seed(int(args.seed))

    g = Genome(args.genome)

    args.func(g, args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='make random sites for bamsurgeon')
    parser.add_argument('-g', '--genome', required=True, help='genome description: chromosome name <TAB/SPACE> chromosome length')
    parser.add_argument('-s', '--seed', default=None, help='use a seed to make reproducible random picks')
    parser.add_argument('-n', '--numpicks', default=1000, help='number of sites to generate (default = 1000)')
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

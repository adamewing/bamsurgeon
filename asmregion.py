#!/bin/env python

"""
try to do ref-directed assembly for paired reads in a region of a .bam file
"""

import pysam,tempfile,argparse,subprocess,sys,shutil,os,re
import parseamos

def velvetContigs(dir):
    assert os.path.exists(dir)
    fh = open(dir + "/contigs.fa", 'r')
    contigs = []
    name = None
    seq  = None 
    for line in fh:
        if re.search("^>", line):
            if name and seq:
                contigs.append(Contig(name,seq,dir))
            name = line.lstrip('>').strip()
            seq = ''
        else:
            if name:
                seq += line.strip()
            else:
                raise ValueError("invalid fasta format: " + fastaFile)
    if name and seq:
        contigs.append(Contig(name,seq,dir))
    return contigs

class Contig:
    def __init__(self,name,seq,dir):
        self.name = name
        self.seq = seq
        self.len = len(seq)

        amosfile = dir + "/velvet_asm.afg"
        seqfile = dir + "/Sequences"
        inseq = parseamos.InputSeqs(seqfile)  # contains ALL input seqs...
        contigreadmap = parseamos.contigreadmap(amosfile,inseq)

        namefields = name.split('_')
        self.iid = namefields[1]

        self.reads = contigreadmap[self.iid] 
        self.rquals = [] # meaningless, used for filler later rather than have uniform quality
        self.mquals = [] # meaningless, used for filler later rather than have uniform quality

        assert self.name
        assert self.seq
        assert self.len
    def __str__(self):
        return ">" + self.name + "\n" + self.seq 

def median(list):
    list.sort()
    i = len(list)/2
    if len(list) % 2 == 0:
        return (list[i] + list[i+1])/2
    else:
        return list[i]

def n50(contigs):
    ln = map(lambda x: x.len, contigs)
    nlist = []
    for n in ln:
        for i in range(n):
            nlist.append(n)
    return median(nlist)

def runVelvet(reads,refseqname,refseq,kmer,isPaired=True,long=False, inputContigs=False, cov_cutoff=False, noref=False):
    """
    reads is either a dictionary of ReadPair objects, (if inputContigs=False) or a list of 
    Contig objects (if inputContigs=True), refseq is a single sequence, kmer is an odd int
    """
    readsFasta  = tempfile.NamedTemporaryFile(delete=False,dir='.')
    refseqFasta = tempfile.NamedTemporaryFile(delete=False,dir='.')

    if inputContigs:
        for contig in reads:
            readsFasta.write(str(contig) + "\n")
    else:
        for readpair in reads.values():
            readsFasta.write(readpair.fasta())

    if refseq:
        refseqFasta.write(">%s\n%s\n" % (refseqname,refseq))

    readsFasta.flush()
    refseqFasta.flush()

    readsFN  = readsFasta.name
    refseqFN = refseqFasta.name

    tmpdir = tempfile.mkdtemp(dir='.')

    print tmpdir

    if noref:
        if long:
            argsvelveth = ['velveth', tmpdir, str(kmer), '-long', readsFN]
        else:
            argsvelveth = ['velveth', tmpdir, str(kmer), '-shortPaired', readsFN]
    else:
        if long: 
            argsvelveth = ['velveth', tmpdir, str(kmer), '-reference', refseqFN, '-long', readsFN]
        else:
            argsvelveth = ['velveth', tmpdir, str(kmer), '-reference', refseqFN, '-shortPaired', readsFN]

    if cov_cutoff:
        argsvelvetg = ['velvetg', tmpdir, '-exp_cov', 'auto', '-cov_cutoff', 'auto', '-unused_reads', 'yes', '-read_trkg', 'yes', '-amos_file', 'yes']
    else:
        argsvelvetg = ['velvetg', tmpdir, '-unused_reads', 'yes', '-read_trkg', 'yes', '-amos_file', 'yes']
        
    subprocess.call(argsvelveth)
    subprocess.call(argsvelvetg)

    vcontigs = velvetContigs(tmpdir)

    # cleanup
    shutil.rmtree(tmpdir)
    os.unlink(readsFN)
    os.unlink(refseqFN)

    return vcontigs


class ReadPair:
    def __init__(self,read,mate):
        assert read.is_read1 != mate.is_read1

        self.read1 = None
        self.read2 = None

        if read.is_read1:
            self.read1 = read
            self.read2 = mate
        else:
            self.read1 = mate
            self.read2 = read

    def fasta(self):
        return ">" + self.read1.qname + "\n" + self.read1.seq + "\n>" + self.read2.qname + "\n" + self.read2.seq + "\n"

    def __str__(self):
        r1map = "mapped"
        r2map = "mapped"

        if self.read1.is_unmapped:
            r1map = "unmapped"
        if self.read2.is_unmapped:
            r2map = "unmapped"

        output = " ".join(("read1:", self.read1.qname, self.read1.seq, r1map, "read2:", self.read2.qname, self.read2.seq, r2map))
        return output

def asm(chr, start, end, bamfilename, reffile, kmersize, noref=False, recycle=False):
    bamfile  = pysam.Samfile(bamfilename,'rb')
    matefile = pysam.Samfile(bamfilename,'rb')

    readpairs = {}
    nreads = 0
    rquals = []
    mquals = []
    for read in bamfile.fetch(chr,start,end):
        if not read.mate_is_unmapped and read.is_paired:
            mate = matefile.mate(read)
            readpairs[read.qname] = ReadPair(read,mate)
            nreads += 1
            if read.is_read1:
                if read.is_reverse:
                    rquals.append(read.qual[::-1])
                    mquals.append(mate.qual)
                else:
                    rquals.append(read.qual)
                    mquals.append(mate.qual[::-1])
            else:
                if read.is_reverse:
                    rquals.append(mate.qual)
                    mquals.append(read.qual[::-1])
                else:
                    rquals.append(mate.qual[::-1])
                    mquals.append(read.qual)

    sys.stderr.write("found " + str(nreads) + " reads in region.\n")
    refseq = None
    if reffile:
        refseq = reffile.fetch(chr,start,end)

    region = chr + ":" + str(start) + "-" + str(end)

    contigs = runVelvet(readpairs, region, refseq, kmersize, cov_cutoff=True, noref=noref)
    newcontigs = None

    if recycle:
        if len(contigs) > 1:
            newcontigs = runVelvet(contigs, region, refseq, kmersize, long=True, inputContigs=True, noref=noref)

        if newcontigs and n50(newcontigs) > n50(contigs):
            contigs = newcontigs

    for contig in contigs:
        contig.rquals = rquals
        contig.mquals = mquals
    return contigs

def main(args):
    '''
    this is here for testing/debugging
    '''
    reffile  = None

    if not args.noref:
        if not args.refFasta:
            raise ValueError("no reference given and --noref not set")
        reffile  = pysam.Fastafile(args.refFasta)

    (chr,coords) = args.regionString.split(':')
    (start,end) = coords.split('-')
    start = re.sub(',','',start)
    end   = re.sub(',','',end)
    start = int(start)
    end   = int(end)

    contigs = asm(chr, start, end, args.bamFileName, reffile, int(args.kmersize), args.noref, args.recycle)

    maxlen = 0
    maxiid = None
    for contig in contigs:
        print contig
        if contig.len > maxlen:
            maxlen = contig.len
            maxiid = contig.iid

    print "read names in longest contig (" + str(maxiid) + "):"
    for contig in contigs:
        if contig.iid == maxiid:
            contig.reads.infodump()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='parse the output of pickreads.py')
    parser.add_argument('-f', '--fastaref', dest='refFasta', required=False,
                        help='indexed reference fo ref-directed assembly')
    parser.add_argument('-r', '--region', dest='regionString', required=True,
                        help='format: chrN:startbasenum-endbasenum')
    parser.add_argument('-b', '--bamfile', dest='bamFileName', required=True,
                        help='target .bam file, region specified in -r must exist')
    parser.add_argument('-k', '--kmersize', dest='kmersize', default=31,
                        help='kmer size for velvet, default=31')
    parser.add_argument('--noref', action="store_true")
    parser.add_argument('--recycle', action="store_true")
    args = parser.parse_args()
    main(args)

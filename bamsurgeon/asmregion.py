#!/usr/bin/env python

"""
try to do ref-directed assembly for paired reads in a region of a .bam file
"""

import pysam,argparse,subprocess,sys,shutil,os,re
import parseamos
import datetime

from uuid import uuid4


def now():
    return str(datetime.datetime.now())


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
        self.eid = namefields[1]

        self.reads = contigreadmap[self.eid]
        self.rquals = [] # meaningless, used for filler later rather than have uniform quality
        self.mquals = [] # meaningless, used for filler later rather than have uniform quality

        assert self.name
        assert self.seq
        assert self.len

    def subseq(self, start, end):
        start = int(start)
        end   = int(end)
        assert start < end
        return self.seq[start:end]

    def trimseq(self, start, end):
        self.seq = self.subseq(start, end)
        self.len = len(self.seq)

        # modify read list
        new_reads = parseamos.ContigReads(self.eid)
        for src, read in self.reads.reads.iteritems():
            if int(read.off) > int(start) and int(read.off) < int(end):
                new_reads.reads[read.src] = read
        assert new_reads is not None
        self.reads = new_reads

    def __len__(self): return self.len

    def __gt__(self, other):
        return self.len > other.len

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


def runVelvet(reads,refseqname,refseq,kmer,addsv_tmpdir,isPaired=True,long=False,cov_cutoff=False,mutid='null',debug=False):
    """
    reads is either a dictionary of ReadPair objects, kmer is an odd int
    """
    reads_fasta_tmpfn = addsv_tmpdir + '/' + '.'.join((mutid, 'tmpreads', 'fasta'))

    reads_fasta = open(reads_fasta_tmpfn, 'w')


    for readpair in reads.values():
        reads_fasta.write(readpair.fasta())

    reads_fasta.close()

    tmpdir = addsv_tmpdir + '/' + mutid + '.' + str(uuid4()).split('-')[0]

    print tmpdir

    if long:
        argsvelveth = ['velveth', tmpdir, str(kmer), '-long', reads_fasta_tmpfn]
    else:
        argsvelveth = ['velveth', tmpdir, str(kmer), '-shortPaired', reads_fasta_tmpfn]

    if cov_cutoff:
        argsvelvetg = ['velvetg', tmpdir, '-exp_cov', 'auto', '-cov_cutoff', 'auto', '-unused_reads', 'yes', '-read_trkg', 'yes', '-amos_file', 'yes']
    else:
        argsvelvetg = ['velvetg', tmpdir, '-unused_reads', 'yes', '-read_trkg', 'yes', '-amos_file', 'yes']
        
    subprocess.call(argsvelveth)
    subprocess.call(argsvelvetg)

    vcontigs = velvetContigs(tmpdir)

    # cleanup
    if not debug:
        shutil.rmtree(tmpdir)
        os.unlink(reads_fasta_tmpfn)

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


def asm(chrom, start, end, bamfilename, reffile, kmersize, tmpdir, mutid='null', debug=False):
    bamfile  = pysam.Samfile(bamfilename,'rb')
    matefile = pysam.Samfile(bamfilename,'rb')

    readpairs = {}
    nreads = 0
    ndisc  = 0 # track discordant reads
    rquals = []
    mquals = []


    for read in bamfile.fetch(chrom,start,end):
        if not read.mate_is_unmapped and read.is_paired:
            try:
                mate = matefile.mate(read)
                readpairs[read.qname] = ReadPair(read,mate)
                nreads += 1
                if not read.is_proper_pair:
                    ndisc  += 1
                if nreads % 1000 == 0:
                    print "INFO\t" + now() + "\t" + mutid + "\tfound mates for", nreads, "reads,", float(ndisc)/float(nreads), "discordant."
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
            except ValueError:
                sys.stderr.write("WARN\t" + now() + "\t" + mutid + "\tcannot find mate for read marked paired: " + read.qname + "\n")

    sys.stdout.write("INFO\t" + now() + "\t" + mutid + "\tfound " + str(nreads) + " reads in region.\n")

    if nreads == 0:
        return []

    refseq = None
    if reffile:
        refseq = reffile.fetch(chrom,start,end)

    region = "%s:%d-%d" % (chrom, start, end)

    contigs = runVelvet(readpairs, region, refseq, kmersize, tmpdir, cov_cutoff=True, mutid=mutid, debug=debug)
    newcontigs = None

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

    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)
        print "INFO\t" + now() + "\tcreated tmp directory: " + args.tmpdir

    contigs = asm(chr, start, end, args.bamFileName, reffile, int(args.kmersize), args.tmpdir, args.noref, args.recycle)

    maxlen = 0
    maxeid = None
    for contig in contigs:
        print contig
        if contig.len > maxlen:
            maxlen = contig.len
            maxeid = contig.eid

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
    parser.add_argument('-t', '--tmpdir', dest='tmpdir', default='addsv.tmp')
    parser.add_argument('--noref', action="store_true")
    parser.add_argument('--recycle', action="store_true")
    args = parser.parse_args()
    main(args)

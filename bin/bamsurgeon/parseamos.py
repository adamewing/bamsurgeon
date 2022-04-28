#!/usr/bin/env python

import sys, re

import logging
FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class Read:
    def __init__(self, src):
        self.src  = src  # read number
        self.off  = None # offset in contig
        self.name = None
        self.seq  = None
    def __str__(self):
        return '\t'.join(map(str, (self.src, self.off, self.name, '\n'+self.seq)))


class ContigReads:
    def __init__(self, eid):
        self.eid = eid
        self.reads = {} # read names

    def addread(self,src):
        #if src in self.reads:
        #    logger.error("AMOS PROBLEM: src already in readlist\n")
        self.reads[src] = Read(src)

    def getReadNames(self,seqs):
        for src, read in self.reads.items():
            self.reads[src].name = seqs.srcread[src]
            self.reads[src].seq  = seqs.seqmap[seqs.srcread[src]]

    def __str__(self):
        output = []
        for eid, read in self.reads.items():
            output.append(str(read))
        return '\n'.join(output)


class InputSeqs:
    def __init__(self,seqfile):
        self.srcread = {}  # read number --> read name
        self.seqmap  = {} # read name --> read seq

        readnames = []
        readseqs  = []

        f = open(seqfile,'r')
        read = None
        src  = None
        seq  = ""
        for line in f:
            if re.search('^>',line):
                if read:
                    assert seq != ""
                    readseqs.append(seq)
                    seq = ""
                (read,src,n) = re.sub('^>','',line).strip().split()
                assert read and src
                self.srcread[src] = read
                readnames.append(read)
            else:
                seq += line.strip()

        self.srcread[src] = read
        readseqs.append(seq)

        self.seqmap = dict(zip(readnames,readseqs))

    def __str__(self):
        output = ""
        for name, seq in self.seqmap.items():
            output += ">" + name + "\n" + seq + "\n"
        return output


def contigreadmap(amosfile,seqs):
    '''
    seqs is an InputSeqs object
    '''
    f = open(amosfile, 'r')
    inCTGblock = False
    inTLEblock = False
    CTGeid = None
    TLEsrc = None
    TLEoff = None
    contig = None
    contigs = {} 

    for line in f:
        line = line.strip()
        if inCTGblock:
            if not inTLEblock:
                if re.search('^eid:',line):
                    CTGeid = re.sub('eid:','',line)
                    CTGeid = re.sub('-0$', '', CTGeid)
                    contig = ContigReads(CTGeid)
                    #print "debug: CTGeid =",CTGeid

                if re.search('}', line):
                    contig.getReadNames(seqs)
                    contigs[CTGeid] = contig
                    inCTGblock = False
                    CTGeid = None
                    contig = None
                    #print "debug: left CTG block"

                if re.search('{TLE',line):
                    inTLEblock = True
                    #print "debug: in TLE block"

            else: # in TLE block
                if re.search('^src:',line):
                    TLEsrc = re.sub('src:','',line)
                    contig.addread(TLEsrc)
                    #print "debug TLEsrc =",TLEsrc

                if re.search('^off:',line):
                    TLEoff = re.sub('off:', '', line)
                    contig.reads[TLEsrc].off = int(TLEoff) # set read offset in contig

                if re.search('}', line):
                    inTLEblock = False
                    TLEsrc = None
                    TLEoff = None
                    #print "debug: left TLE block"

        else: # not in CTG block
            if re.search('{CTG',line):
                #print "debug: in CTG block"
                inCTGblock = True

    return contigs


if __name__ == '__main__':
    '''
    This is here for testing/debugging
    '''
    if len(sys.argv) == 2:
        amosfile = sys.argv[1].strip() + "/velvet_asm.afg"
        seqfile  = sys.argv[1].strip() + "/Sequences"
        inputseqs = InputSeqs(seqfile)
        contigmap = contigreadmap(amosfile,inputseqs)
        for eid,contig in contigmap.items():
            print(contig)

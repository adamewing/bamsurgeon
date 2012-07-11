#!/bin/env python

import sys,re

class ContigReads:
    def __init__(self, iid):
        self.iid = iid 
        self.srcs = [] # read numbers from velvet
        self.reads = [] # read names
    def addsrc(self,src):
        self.srcs.append(src)
    def getReadNames(self,seqs):
        for src in self.srcs:
            read = seqs.srcread[src]
            self.reads.append(read)
    def infodump(self):
        for i in range(len(self.srcs)):
            print self.srcs[i],self.reads[i]

class InputSeqs:
    def __init__(self,seqfile):
        self.srcread = {}  # read number --> read name
        self.readnames = []  # read names
        self.readseqs = []

        f = open(seqfile,'r')
        read = None
        src  = None
        seq  = ""
        for line in f:
            if re.search('^>',line):
                if read:
                    assert seq != ""
                    self.readseqs.append(seq)
                    seq = ""
                (read,src,n) = re.sub('^>','',line).strip().split()
                assert read and src
                self.srcread[src] = read
                self.readnames.append(read)
            else:
                seq += line.strip()
        self.srcread[src] = read
        self.readseqs.append(seq)

    def __str__(self):
        output = ""
        for i in range(len(self.readnames)):
            output += ">" + self.readnames[i] + "\n" + self.readseqs[i] + "\n"
        return output

def contigreadmap(amosfile,seqs):
    '''
    seqs is an InputSeqs object
    '''
    f = open(amosfile, 'r')
    inCTGblock = False
    inTLEblock = False
    CTGiid = None
    TLEsrc = None
    contig = None
    contigs = {} 

    for line in f:
        line = line.strip()
        if inCTGblock:
            if not inTLEblock:
                if re.search('^iid:',line):
                    CTGiid = re.sub('iid:','',line)
                    contig = ContigReads(CTGiid)
                    #print "debug: CTGiid =",CTGiid

                if re.search('}', line):
                    contig.getReadNames(seqs)
                    contigs[CTGiid] = contig
                    inCTGblock = False
                    CTGiid = None
                    contig = None
                    #print "debug: left CTG block"

                if re.search('{TLE',line):
                    inTLEblock = True
                    #print "debug: in TLE block"

            else: # in TLE block
                if re.search('^src:',line):
                    TLEsrc = re.sub('src:','',line)
                    contig.addsrc(TLEsrc)
                    #print "debug TLEsrc =",TLEsrc

                if re.search('}', line):
                    inTLEblock = False
                    TLEsrc = None
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
    amosfile = sys.argv[1].strip() + "/velvet_asm.afg"
    seqfile  = sys.argv[1].strip() + "/Sequences"
    inputseqs = InputSeqs(seqfile)
    contigmap = contigreadmap(amosfile,inputseqs)
    for iid,contig in contigmap.iteritems():
        contig.infodump()

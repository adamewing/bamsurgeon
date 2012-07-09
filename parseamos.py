#!/bin/env python

import sys,re

class ContigReads:
    def __init__(self, iid):
        self.iid = iid 
        self.srcdict = {} # TLEsrc --> read 
        self.reads   = {} # read --> TLEsrc
    def addsrc(self,src):
        self.srcdict[src] = 1
    def getReadNames(self,seqs):
        for src in self.srcdict.keys():
            read = seqs.srcread[src]
            self.srcdict[src] = read
            self.reads[read]  = src
    def infodump(self):
        for srcdict, read in self.srcdict.iteritems():
            print self.iid,srcdict,read


class InputSeqs:
    def __init__(self,seqfile):
        self.srcread = {}  # read number --> read name
        self.readseq = {}  # read name --> read seq

        f = open(seqfile,'r')
        read = None
        src  = None
        seq  = ""
        for line in f:
            if re.search('^>',line):
                if read:
                    assert seq != ""
                    self.readseq[read] = seq
                    seq = ""
                (read,src,n) = re.sub('^>','',line).strip().split()
                assert read and src
                self.srcread[src] = read
            else:
                seq += line.strip()

    def __str__(self):
        output = ""
        for read,seq in self.readseq.iteritems():
            output += ">" + read + "\n" + seq + "\n"
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
    amosfile = sys.argv[1].strip() + "/velvet_asm.afg"
    seqfile  = sys.argv[1].strip() + "/Sequences"
    inputseqs = InputSeqs(seqfile)
    contigmap = contigreadmap(amosfile,inputseqs)
    for iid,contig in contigmap.iteritems():
        contig.infodump()

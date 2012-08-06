#!/bin/env python

'''
Methods for making mutations in a sequence
'''

import string

def rc(seq):
    seq = seq[::-1]
    seq = seq.translate(string.maketrans("ATGC","TACG"))
    return seq

class MutableSeq:
    def __init__(self,seq):
        self.seq = string.upper(seq.strip())

    def __str__(self):
        return self.seq

    def length(self):
        return len(self.seq)

    def subseq(self, start, end):
        start = int(start)
        end   = int(end)
        assert start < end
        return self.seq[start:end]

    def deletion(self, start, end):
        """
        deletes between start and end, bases at positions corresponding to start and end are kept
        """
        start = int(start)
        end   = int(end)
        assert start < end
        self.seq = self.seq[:start] + self.seq[end:]

    def insertion(self, loc, seq, tsdlen=0):
        """
        inserts seq after position loc, adds taret site duplication (tsd) if tsdlen > 0
        """
        tsd = self.seq[loc:loc+tsdlen]
        self.seq = self.seq[:loc] + tsd + seq + self.seq[loc+1:]

    def inversion(self, start, end):
        """
        inverts sequence between start and end, bases at start and end positions are not affected
        """
        start = int(start)
        end   = int(end)
        assert start < end
        invseq = rc(self.seq[start:end])
        self.seq = self.seq[:start] + invseq + self.seq[end:]

    def duplication(self,start,end,fold=1):
        """
        duplicates sequence between start and end, bases at start and end positions are not included
        """
        start = int(start)
        end   = int(end)
        assert start < end
        dupseq = self.seq[start:end]
        for i in range(fold):
            dupseq = dupseq + self.seq[start:end]
        self.seq = self.seq[:start] + dupseq + self.seq[end:]
        
        

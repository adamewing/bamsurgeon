#!i/usr/bin/env python

'''
Methods for making mutations in a sequence
'''

import operator
from bamsurgeon.common import rc


def dist(seq1, seq2):
    ''' Hamming distance '''
    assert len(seq1) == len(seq2)
    if seq1 == seq2:
        return 0

    return sum([c1 != c2 for c1, c2 in zip(seq1, seq2)])


class MutableSeq:
    def __init__(self,seq):
        self.seq = str.upper(seq.strip())

    def __str__(self):
        return self.seq

    def find_site(self, site, left_trim=0, right_trim=0):
        ''' find closest occurance of site in self.seq '''
        left, right = site.split('^')
        site = left+right

        assert len(site) <= len(self.seq)

        d = [dist(self.seq[start:start+len(site)], site) for start in range(left_trim,len(self.seq)-len(site)-right_trim+1)]
        min_i, min_d = min(enumerate(d), key=operator.itemgetter(1))

        return min_i + len(left)

    def length(self):
        return len(self.seq)

    def subseq(self, start, end):
        start = int(start)
        end   = int(end)
        assert start < end
        return self.seq[start:end]

    def deletion(self, start, end):
        ''' deletes between start and end, bases at positions corresponding to start and end are kept '''
        start = int(start)
        end   = int(end)
        assert start < end
        self.seq = self.seq[:start] + self.seq[end:]

    def insertion(self, loc, seq, tsdlen=0):
        ''' inserts seq after position loc, adds taret site duplication (tsd) if tsdlen > 0 '''
        loc = int(loc)
        tsd = self.seq[loc:loc+tsdlen]
        self.seq = self.seq[:loc] + tsd + seq + self.seq[loc:]

    def inversion(self, start, end):
        ''' inverts sequence between start and end, bases at start and end positions are not affected '''
        start = int(start)
        end   = int(end)
        assert start < end
        invseq = rc(self.seq[start:end])
        self.seq = self.seq[:start] + invseq + self.seq[end:]

    def duplication(self,start,end,fold=1):
        ''' duplicates sequence between start and end, bases at start and end positions are not included '''
        start = int(start)
        end   = int(end)
        assert start < end
        dupseq = self.seq[start:end]
        for i in range(fold):
            dupseq = dupseq + self.seq[start:end]
        self.seq = self.seq[:start] + dupseq + self.seq[end:]

    def fusion(self, loc1, other, loc2, flip1=False, flip2=False):
        loc1 = int(loc1)
        loc2 = int(loc2)

        loc1_seq = self.seq[:loc1]
        loc2_seq = other.seq[loc2:]

        if flip1:
            loc1_seq = rc(self.seq[loc1:])

        if flip2:
            loc2_seq = rc(other.seq[:loc2])

        self.seq = loc1_seq + loc2_seq 

        
        

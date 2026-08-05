#!/usr/bin/env python


import subprocess
import pysam
import random
import hashlib
import string
import os

from collections import Counter
from shutil import move
from re import sub

import logging
FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def make_salt(seed=None):
    '''
    Build the salt that read selection is hashed against.

    This must be created once in the parent process and passed down to
    workers. It used to be a lazily-initialised module global, which meant
    every process in a ProcessPoolExecutor drew its own, so the set of
    replaced reads changed with --procs.
    '''
    if seed is not None:
        return hashlib.md5(('bamsurgeon:' + str(seed)).encode()).hexdigest()[:10]

    alphabet = string.ascii_lowercase + string.ascii_uppercase + string.digits
    return ''.join(random.choice(alphabet) for _ in range(10))


def read_hash_fraction(query_name, salt):
    ''' stable pseudo-random value in [0,1) for a read name '''
    read_hash = int(hashlib.md5((salt + query_name).encode()).hexdigest(), 16)
    return (read_hash % 1000000) / 1000000.0


def select_qnames(qnames, frac, salt):
    '''
    Pick exactly round(frac * N) of qnames, deterministically.

    The previous approach kept each read whose hash fell below frac, which
    makes the count binomially distributed around the target rather than
    equal to it -- a direct contributor to spiked depth not matching the
    requested VAF. Ranking by hash and taking a prefix gives the same
    pseudo-random membership with an exact count.
    '''
    if frac >= 1.0:
        return set(qnames)
    if frac <= 0.0:
        return set()

    # qname breaks ties so the ordering is total, not just stable
    ranked = sorted(qnames, key=lambda q: (read_hash_fraction(q, salt), q))
    return set(ranked[:int(round(frac * len(ranked)))])


def get_avg_coverage(alignment_file, chrom, start, end):
    split_coverage = alignment_file.count_coverage(chrom, start, end, quality_threshold=0)
    base_sum = [sum(x) for x in split_coverage]
    return sum(base_sum) / float(end - start)


def rc(dna):
    ''' reverse complement '''
    complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


def majorbase(basepile):
    """returns tuple: (major base, count)
    """
    return Counter(basepile).most_common()[0]


def minorbase(basepile):
    """returns tuple: (minor base, count)
    """
    c = Counter(basepile)
    if len(list(c.elements())) > 1:
        return c.most_common(2)[-1]
    else:
        return c.most_common()[0]


def mergebams(bamlist, outbamfn, maxopen=100, debug=False):
    ''' call samtools to merge a list of bams hierarchically '''

    assert outbamfn.endswith('.bam')
    logger.info("len(bamlist): %d" % len(bamlist))

    if len(bamlist) == 1:
        logger.info("only one BAM to merge, renaming " + bamlist[0] + " --> " + outbamfn)
        move(bamlist[0], outbamfn)
    else:
        nmerge = 1
        mergenum = 0
        merge_sublists = {}
        for tmpbam in bamlist:
            mergetmp = "tmp.merging." + str(mergenum) + "." + outbamfn
            if mergetmp in merge_sublists:
                merge_sublists[mergetmp].append(tmpbam)
            else:
                merge_sublists[mergetmp] = []
                merge_sublists[mergetmp].append(tmpbam)
            if nmerge % maxopen == 0:
                mergenum += 1
            nmerge += 1

        for submergefn, tmpbams in merge_sublists.items():
            if len(tmpbams) == 1:
                move(tmpbams[0], submergefn)
                logger.info("renamed: " + tmpbams[0] + " --> " + submergefn)
            else:
                args = ['samtools', 'merge', '-f', submergefn] + tmpbams
                logger.info("merging, cmd: " + ' '.join(args))
                subprocess.check_call(args)

        if len(merge_sublists.keys()) == 1:
            logger.info("merge finished, renaming: " + list(merge_sublists.keys())[0] + " --> " + outbamfn)
            move(list(merge_sublists.keys())[0], outbamfn)
        else:
            args = ['samtools', 'merge', '-f', outbamfn] + list(merge_sublists.keys())
            logger.info("final merge, cmd: " + ' '.join(args))
            subprocess.check_call(args)

        for submergefn in merge_sublists.keys():
            if os.path.exists(submergefn):
                os.remove(submergefn)

    if not debug:
        for bamfile in bamlist:
            if os.path.exists(bamfile):
                os.remove(bamfile)
            if os.path.exists(bamfile + '.bai'):
                os.remove(bamfile + '.bai')


def double_fastq_to_interleaved(fq1, fq2, outfq):
    def get_fastq_read_iterator(fq):
        f = open(fq, 'r')
        while True:
            name = f.readline()
            if not name:
                break
            other = [f.readline() for i in range(3)]
            yield [name] + other
        f.close()
    f_1_iter = get_fastq_read_iterator(fq1)
    f_2_iter = get_fastq_read_iterator(fq2)
    f_out = open(outfq, 'w')
    while True:
        try:
            f_out.write(''.join(next(f_1_iter)))
            f_out.write(''.join(next(f_2_iter)))
        except StopIteration:
            break
    f_out.close()


def _fastq_entry(read, suffix=''):
    '''
    One FASTQ record for an aligned read, in original sequencing orientation.

    A reverse-strand alignment stores the reverse complement of what the
    sequencer read, so the sequence has to be complemented back and the
    quality string reversed alongside it. Getting only one of the two right
    is silent: the realignment still succeeds, with the qualities attached to
    the wrong end of the read.
    '''
    seq = read.query_sequence
    if not seq:
        return None

    quals = read.query_qualities
    if quals is None:
        qual = 'I' * len(seq)
    else:
        qual = pysam.qualities_to_qualitystring(quals)

    if read.is_reverse:
        seq = rc(seq)
        qual = qual[::-1]

    return '@%s%s\n%s\n+\n%s\n' % (read.query_name, suffix, seq, qual)


def bamtofastq(bam, paired=True, twofastq=False, fasta_ref=None):
    '''
    Convert a BAM to FASTQ. Two files if twofastq, otherwise interleaved.

    This used to shell out to picard SamToFastq, which was the only thing in
    bamsurgeon that needed java. The BAMs handed to it are per-mutation temp
    files of at most a few thousand reads, so doing it in-process is both
    faster and one fewer dependency.

    Secondary and supplementary alignments are excluded, matching picard's
    INCLUDE_NON_PRIMARY_ALIGNMENTS=false: their sequences are duplicates or
    hard-clipped fragments, and either would corrupt the read on remapping.

    Record order is first-encounter order of the completing mate, which is
    deterministic for a given input BAM but is not picard's order, so FASTQ
    files and the BAMs realigned from them will not be byte-identical to
    those produced before this change.
    '''
    assert bam.endswith('.bam')

    logger.info("converting BAM " + bam + " to FASTQ")

    inbam = pysam.AlignmentFile(bam, 'rb', reference_filename=fasta_ref)

    if not paired:
        outfq = sub('bam$', 'fastq', bam)
        with open(outfq, 'w') as out:
            for read in inbam.fetch(until_eof=True):
                if read.is_secondary or read.is_supplementary:
                    continue
                entry = _fastq_entry(read)
                if entry is not None:
                    out.write(entry)
        inbam.close()
        return [outfq]

    if twofastq:
        outfqs = [sub('bam$', '1.fastq', bam), sub('bam$', '2.fastq', bam)]
        out1 = open(outfqs[0], 'w')
        out2 = open(outfqs[1], 'w')
    else:
        outfqs = [sub('bam$', 'fastq', bam)]
        out1 = out2 = open(outfqs[0], 'w')

    # mates can be arbitrarily far apart in a coordinate-sorted BAM, so hold
    # each until its partner turns up and write the pair together
    pending = {}

    for read in inbam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary:
            continue

        mate = pending.pop(read.query_name, None)
        if mate is None:
            pending[read.query_name] = read
            continue

        first, second = (mate, read) if mate.is_read1 else (read, mate)

        entries = (_fastq_entry(first, '/1'), _fastq_entry(second, '/2'))
        if None in entries:
            continue

        out1.write(entries[0])
        out2.write(entries[1])

    out1.close()
    if out2 is not out1:
        out2.close()
    inbam.close()

    if pending:
        # picard raised on this; there is nothing useful to align a lone mate
        # against, so report it and carry on rather than losing the mutation
        logger.warning("%s: %d reads had no mate in the file and were skipped"
                       % (bam, len(pending)))

    return outfqs


def check_min_read_count(bamfile, fasta_ref, min_num_reads=0):
    # Returns True if the BAM has at least min_num_reads reads
    bam = pysam.AlignmentFile(bamfile, reference_filename=fasta_ref)
    # If the BAM is indexed, we can just check the header
    try:
        if bam.is_bam and bam.check_index():
            return bam.mapped + bam.unmapped > min_num_reads
    except AttributeError:
        pass
    except ValueError:
        pass
    # If the BAM is not indexed, we have to count the reads manually
    for i, _ in enumerate(bam.fetch(until_eof=True)):
        if i + 1 > min_num_reads:
            return True
    return False


def dictlist(fn):
    d = {}
    with open(fn, 'r') as inlist:
        for name in inlist:
            d[name.strip()] = True
    return d

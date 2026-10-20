'''
Tests for bamsurgeon/replace_reads.py.

These need only pysam, no aligner, no subprocess -- so like test_evaluator.py
they run anywhere the package imports.
'''

import os

import pysam

from bamsurgeon.replace_reads import replace_reads

HEADER = {'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chr1', 'LN': 1000}]}


def make_read(name, is_read1, qual_char, pos):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = 'ACGT'
    a.flag = 1 | (64 if is_read1 else 128)
    a.reference_id = 0
    a.reference_start = pos
    a.mapping_quality = 60
    a.cigar = [(0, 4)]
    a.query_qualities = pysam.qualitystring_to_array(qual_char * 4)
    return a


def write_bam(path, reads):
    with pysam.AlignmentFile(path, 'wb', header=HEADER) as bam:
        for r in reads:
            bam.write(r)
    pysam.index(path)


def test_orphaned_mate_keeps_original_read_instead_of_crashing(tmp_path):
    '''
    replace_reads() used to crash with AttributeError: 'NoneType' object has
    no attribute 'qual' when the donor bam only contains one mate of a pair
    (the other mate fell outside the mutation window). It should instead
    keep the original read for the mate that has no donor copy.
    '''
    target_path = str(tmp_path / 'target.bam')
    donor_path = str(tmp_path / 'donor.bam')
    out_path = str(tmp_path / 'out.bam')

    # target has a full pair: "readA" mate1 (pos 10) + mate2 (pos 500, far away)
    write_bam(target_path, [
        make_read('readA', True, '#', 10),
        make_read('readA', False, '#', 500),
    ])
    # donor only covers mate1, realigned to pos 50 by the mutation step;
    # mate2 fell outside the mutation window and has no donor copy
    write_bam(donor_path, [
        make_read('readA', True, 'F', 50),
    ])

    replace_reads(target_path, donor_path, out_path, fasta_ref=None,
                  keepqual=True, quiet=True)

    with pysam.AlignmentFile(out_path, 'rb') as out:
        recs = [(r.is_read1, r.reference_start) for r in out.fetch(until_eof=True)]

    # mate1 was replaced by the donor's realigned copy (pos 50); mate2 had no
    # donor copy so it must fall back to the original target read (pos 500)
    # instead of crashing on a None donor read. Checked as a list, not a set,
    # so a duplicate write of either mate would be caught.
    assert sorted(recs) == [(False, 500), (True, 50)]

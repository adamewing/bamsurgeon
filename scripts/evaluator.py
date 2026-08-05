#!/usr/bin/env python
'''
Evaluate VCFs against BAMSurgeon "Truth" VCFs
Adam Ewing, adam.ewing@mater.uq.edu.au

Requires pysam. This used to require PyVCF, which has been unmaintained since
2018; the record classification helpers below (is_snv/is_indel/is_sv,
passfilter) reimplement the PyVCF record properties it relied on.
'''

import sys, os
import argparse
import pysam
from collections import OrderedDict
import logging

FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


# ---------------------------------------------------------------------------
# record classification
#
# pysam exposes the raw fields; PyVCF had is_snp/is_indel/is_sv properties.
# ---------------------------------------------------------------------------

def alts(rec):
    return [a for a in (rec.alts or ()) if a not in (None, '.')]


def info_get(rec, key, default=None):
    '''
    INFO lookup that tolerates undeclared keys.

    pysam raises ValueError('Invalid header') rather than returning None for
    a key the file's header does not define, and a submission VCF routinely
    lacks fields the truth VCF declares.
    '''
    try:
        value = rec.info.get(key, default)
    except (ValueError, KeyError):
        return default
    return default if value is None else value


def is_sv(rec):
    ''' SVTYPE, a symbolic ALT, or a breakend '''
    if info_get(rec, 'SVTYPE') is not None:
        return True
    return any(a.startswith('<') or '[' in a or ']' in a for a in alts(rec))


def is_snv(rec):
    a = alts(rec)
    return bool(a) and not is_sv(rec) and len(rec.ref) == 1 and all(len(x) == 1 for x in a)


def is_indel(rec):
    a = alts(rec)
    return bool(a) and not is_sv(rec) and any(len(x) != len(rec.ref) for x in a)


def of_type(rec, vtype):
    return {'SNV': is_snv, 'INDEL': is_indel, 'SV': is_sv}[vtype](rec)


def reckey(rec):
    '''
    Identity of a record, for tracking which truth sites have been consumed.

    This used to be str(record), whose meaning differed between the PyVCF
    repr and a full VCF line.
    '''
    return (rec.chrom, rec.pos, rec.ref, tuple(alts(rec)))


def match(subrec, trurec, vtype='SNV'):
    assert vtype in ('SNV', 'SV', 'INDEL')

    if vtype == 'SNV' and is_snv(subrec) and is_snv(trurec):
        if subrec.pos == trurec.pos and subrec.ref == trurec.ref and alts(subrec) == alts(trurec):
            return True

    if vtype == 'INDEL' and is_indel(subrec) and is_indel(trurec):
        if subrec.pos == trurec.pos and subrec.ref == trurec.ref and alts(subrec) == alts(trurec):
            return True

    if vtype == 'SV' and is_sv(subrec) and is_sv(trurec):
        trustart, truend = expand_sv_ends(trurec)
        substart, subend = expand_sv_ends(subrec)

        # check for overlap
        if min(truend, subend) - max(trustart, substart) > 0:
            return True

    return False


def expand_sv_ends(rec):
    '''
    Assign start and end positions to SV calls, widening by the confidence
    intervals when present.

    htslib folds INFO/END into rec.stop, so no INFO lookup is needed for the
    end. The previous version did `int(rec.INFO.get('END')[0])`, which
    subscripts a scalar and raised TypeError for every SV; the handler that
    was meant to report that referenced an undefined `logger`, so SV
    evaluation always died with a NameError instead.
    '''
    startpos, endpos = rec.start, rec.stop

    cipos = info_get(rec, 'CIPOS')
    if cipos and cipos[0] < 0:
        startpos += cipos[0]

    ciend = info_get(rec, 'CIEND')
    if ciend and len(ciend) > 1 and ciend[1] > 0:
        endpos += ciend[1]

    if startpos > endpos:
        endpos, startpos = startpos, endpos

    return startpos, endpos


def relevant(rec, vtype, ignorechroms):
    ''' Return true if a record matches the type of variant being investigated '''
    return of_type(rec, vtype) and (ignorechroms is None or rec.chrom not in ignorechroms)


def passfilter(rec, disabled=False):
    ''' Return true if a record is unfiltered or has PASS in the filter field '''
    if disabled:
        return True
    keys = list(rec.filter.keys())
    return not keys or keys == ['PASS'] or keys == ['.']


def svmask(rec, vcfh, truchroms):
    ''' mask snv calls in sv regions '''
    if is_snv(rec) and rec.chrom in truchroms:
        for overlap_rec in vcfh.fetch(rec.chrom, rec.pos-1, rec.pos):
            if is_sv(overlap_rec):
                return True
    return False


def var_dist(v1, v2):
    """compute absolute distance between two variants
    """
    assert v1.chrom == v2.chrom
    return abs(v1.pos-v2.pos)


def get_close_matches(var, vcf_fh, win, indels_only=True):
    """Find close matches for variant var in file vcf_fh within given window
    win and return as list of item,distance tuples, sorted ascendingly by
    distance.

    """

    matches = list(vcf_fh.fetch(var.chrom, max(0, var.pos-win), var.pos+1+win))
    if indels_only:
        matches = [m for m in matches if is_indel(m)]
    if len(matches) == 0:
        return []

    dist_map = [(m, var_dist(m, var)) for m in matches]
    return sorted(dist_map, key=lambda x: x[1])


def have_identical_haplotypes(v1, v2, ref):
    """Check if two variant produce the same haplotype / variant sequence.

    - v1 and v2: variant records to compare
    - ref: PySAM FastaFile
    """

    assert (is_indel(v1) or is_snv(v1)) and (is_indel(v2) or is_snv(v2))

    if v1.chrom != v2.chrom:
        return False

    v1_alts, v2_alts = alts(v1), alts(v2)

    if is_snv(v1) and is_snv(v2):
        assert v1.ref.upper() == v2.ref.upper()
        return v1_alts[0].upper() == v2_alts[0].upper()
    if is_snv(v1) or is_snv(v2):
        # one snp one indel: can't produce identical results
        return False

    assert is_indel(v1) and is_indel(v2)
    # only one allele per variant allowed
    assert len(v1_alts) == 1 and len(v2_alts) == 1, "Can't handle multi-allelic entries"

    # get the sequence context which fully overlaps both variants.
    # note: VCF is one-based, but start and end are zero-based half-open
    start = min([v1.pos, v2.pos])-1
    end = max([v1.pos + max([len(v1.ref), len(v1_alts[0])]),
               v2.pos + max([len(v2.ref), len(v2_alts[0])])
               ])
    chrom = v1.chrom  # made sure before they are identical before
    seq = list(ref.fetch(chrom, start, end).upper())

    if len(seq) != end-start:
        logger.warning("Couldn't fetch full sequence window. Skipping"
                       " allele-aware comparison, otherwise indices would"
                       " be off")
        raise NotImplementedError

    v1_offset = v1.pos-1-start
    v2_offset = v2.pos-1-start
    # lower() in replacement for debugging purposes only
    v1_seq = seq[:v1_offset] + list(v1_alts[0].lower()) + seq[v1_offset+len(v1.ref):]
    v2_seq = seq[:v2_offset] + list(v2_alts[0].lower()) + seq[v2_offset+len(v2.ref):]

    assert seq[v1_offset] == v1.ref[0].upper()
    assert seq[v2_offset] == v2.ref[0].upper()
    assert len(v1_seq) == len(seq) - len(v1.ref) + len(v1_alts[0])
    assert len(v2_seq) == len(seq) - len(v2.ref) + len(v2_alts[0])

    return ''.join(v1_seq).upper() == ''.join(v2_seq).upper()


def evaluate(submission, truth, vtype='SNV', reffa=None, ignorechroms=None, ignorepass=False,
             fp_vcf=None, fn_vcf=None, tp_vcf=None,
             debug=False):
    ''' return stats on sensitivity, specificity, balanced accuracy '''

    assert vtype in ('SNV', 'SV', 'INDEL')

    subvcfh = pysam.VariantFile(submission)

    # one handle to scan the truth file end to end, a second to fetch regions
    # from it while that scan's results are in hand
    truvcfh = pysam.VariantFile(truth)
    trufetch = pysam.VariantFile(truth)

    fpvcfh = fnvcfh = tpvcfh = None
    if fp_vcf:
        fpvcfh = pysam.VariantFile(fp_vcf, 'w', header=subvcfh.header)
    if tp_vcf:
        tpvcfh = pysam.VariantFile(tp_vcf, 'w', header=subvcfh.header)
    if fn_vcf:
        # false negatives are truth records, so they need the truth header;
        # writing them against the submission header would drop or mangle
        # INFO fields the submission does not declare
        fnvcfh = pysam.VariantFile(fn_vcf, 'w', header=truvcfh.header)

    reffa_fh = None
    if reffa:
        reffa_fh = pysam.Fastafile(reffa)
        if debug:
            print("DEBUG: Using haplotype aware indel comparison")

    tpcount = 0
    fpcount = 0
    subrecs = 0
    trurecs = 0

    truchroms = {}
    fns = OrderedDict()

    ''' count records in truth vcf, track contigs/chromosomes '''
    for trurec in truvcfh:
        if relevant(trurec, vtype, ignorechroms):
            trurecs += 1
            truchroms[trurec.chrom] = True
            fns[reckey(trurec)] = trurec

    used_truth = {}  # keep track of 'truth' sites used, they should only be usable once

    ''' parse submission vcf, compare to truth '''
    for subrec in subvcfh:
        if passfilter(subrec, disabled=ignorepass):
            if is_snv(subrec) and vtype == 'SNV':
                if not svmask(subrec, trufetch, truchroms):
                    subrecs += 1
            if is_sv(subrec) and vtype == 'SV':
                subrecs += 1
            if is_indel(subrec) and vtype == 'INDEL':
                subrecs += 1

        # the truth record this submission matched, if any. The previous
        # version tracked the loop variable instead, so it consumed and
        # deleted whichever truth record the fetch happened to yield last
        # rather than the one that actually matched.
        matched_tru = None

        startpos, endpos = subrec.start, subrec.stop

        if vtype == 'SV' and is_sv(subrec):
            startpos, endpos = expand_sv_ends(subrec)

        try:
            if relevant(subrec, vtype, ignorechroms) and passfilter(subrec, disabled=ignorepass) and subrec.chrom in truchroms:
                for trurec in trufetch.fetch(subrec.chrom, max(0, startpos), endpos):
                    if match(subrec, trurec, vtype=vtype) and reckey(trurec) not in used_truth:
                        matched_tru = trurec
                        break

                if matched_tru is None and is_indel(subrec) and reffa_fh:  # try haplotype aware comparison
                    window = 100
                    for (trurec, _) in get_close_matches(subrec, trufetch, window, indels_only=True):
                        if reckey(trurec) in used_truth:
                            continue
                        if have_identical_haplotypes(subrec, trurec, reffa_fh):
                            matched_tru = trurec
                            if debug:
                                print("DEBUG: Rescuing %s which has same haplotype as %s" % (subrec, trurec))
                            break

            if matched_tru is not None:
                used_truth[reckey(matched_tru)] = True

        except (ValueError, NotImplementedError) as e:
            logger.error("Warning: " + str(e))

        if matched_tru is not None:
            tpcount += 1
            if tpvcfh:
                tpvcfh.write(subrec)
            fns.pop(reckey(matched_tru), None)
        else:
            if relevant(subrec, vtype, ignorechroms) and passfilter(subrec, disabled=ignorepass) and not svmask(subrec, trufetch, truchroms):
                fpcount += 1  # FP counting method needs to change for real tumors
                if fpvcfh:
                    fpvcfh.write(subrec)

    if fnvcfh:
        for fn in fns.values():
            fnvcfh.write(fn)

    print(f"tpcount, fpcount, subrecs, trurecs: {tpcount},{fpcount},{subrecs},{trurecs}")

    recall = float(tpcount) / float(trurecs) if trurecs else 0.0
    if tpcount+fpcount > 0:
        precision = float(tpcount) / float(tpcount + fpcount)
    else:
        precision = 0.0
    f1score = 0.0 if tpcount == 0 else 2.0*(precision*recall)/(precision+recall)

    for fh in [fpvcfh, fnvcfh, tpvcfh]:
        if fh:
            fh.close()

    subvcfh.close()
    truvcfh.close()
    trufetch.close()

    return precision, recall, f1score


def main(args):

    chromlist = None
    if args.chromlist is not None:
        chromlist = args.chromlist.split(',')

    if not args.subvcf.endswith('.vcf') and not args.subvcf.endswith('.vcf.gz'):
        logger.error("submission VCF filename does not end in .vcf or .vcf.gz")
        sys.exit(1)

    if not os.path.exists(args.truvcf):
        logger.error("truth VCF does not exist.")
        sys.exit(1)
    if not os.path.exists(args.truvcf + '.tbi') and not os.path.exists(args.truvcf + '.csi'):
        logger.error("truth VCF does not appear to be indexed. bgzip + tabix index required.")
        sys.exit(1)

    if args.mutype not in ('SV', 'SNV', 'INDEL'):
        logger.error("-m/--mutype must be either SV, SNV, or INDEL")
        sys.exit(1)

    result = evaluate(args.subvcf, args.truvcf, vtype=args.mutype,
                      reffa=args.reffa, ignorechroms=chromlist, ignorepass=args.nonpass,
                      fp_vcf=args.fp_vcf, fn_vcf=args.fn_vcf, tp_vcf=args.tp_vcf,
                      debug=args.debug)

    print("precision, recall, F1 score: " + ','.join(map(str, result)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="check vcf output against a 'truth' vcf")
    parser.add_argument('-v',  '--vcf',    dest='subvcf', required=True, help="VCF being submitted for evaluation")
    parser.add_argument('-t',  '--truth',  dest='truvcf', required=True, help="'Truth' VCF containing true positives")
    parser.add_argument('-f',  '--ref', dest='reffa', help="Reference fasta file (enables haplotype-ware indel comparison)")
    parser.add_argument('-m', '--mutype', dest='mutype', required=True, help="Mutation type: must be either SNV, SV, or INDEL")
    parser.add_argument('--ignore', dest='chromlist', default=None, help="(optional) comma-seperated list of chromosomes to ignore")
    parser.add_argument('--nonpass', dest='nonpass', action="store_true", help="evaluate all records (not just PASS records) in VCF")
    parser.add_argument('--fp', dest='fp_vcf', help="print false positive positions to this vcf-file")
    parser.add_argument('--tp', dest='tp_vcf', help="print true positive positions to this file")
    parser.add_argument('--fn', dest='fn_vcf', help="print false negatives positions to this file")
    parser.add_argument('--debug', dest='debug', action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    main(args)

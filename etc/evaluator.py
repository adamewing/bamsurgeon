#!/usr/bin/env python

import sys, os
import vcf
import argparse
import string
import pysam
from collections import OrderedDict

'''
Evaluate VCFs against BAMSurgeon "Truth" VCFs
Adam Ewing, adam.ewing@mater.uq.edu.au
Requires PyVCF (https://github.com/jamescasbon/PyVCF)
and pysam.
'''

def match(subrec, trurec, vtype='SNV'):
    assert vtype in ('SNV', 'SV', 'INDEL')

    if vtype == 'SNV' and subrec.is_snp and trurec.is_snp:
        if subrec.POS == trurec.POS and subrec.REF == trurec.REF and subrec.ALT == trurec.ALT:
            return True

    if vtype == 'INDEL' and subrec.is_indel and trurec.is_indel:
        if subrec.POS == trurec.POS and subrec.REF == trurec.REF and subrec.ALT == trurec.ALT:
            return True

    if vtype == 'SV' and subrec.is_sv and trurec.is_sv:
        trustart, truend = expand_sv_ends(trurec)
        substart, subend = expand_sv_ends(subrec)

        # check for overlap
        if min(truend, subend) - max(trustart, substart) > 0:
            return True

    return False


def expand_sv_ends(rec):
    ''' assign start and end positions to SV calls using conf. intervals if present '''
    startpos, endpos = rec.start, rec.end
    assert rec.is_sv

    try:
        endpos = int(rec.INFO.get('END')[0])

        if rec.INFO.get('CIPOS'):
            ci = map(int, rec.INFO.get('CIPOS'))
            if ci[0] < 0:
                startpos += ci[0]

        if rec.INFO.get('CIEND'):
            ci = map(int, rec.INFO.get('CIEND')) 
            if ci[0] > 0:
                endpos += ci[0]

    except TypeError as e:
        sys.stderr.write("error expanding sv interval: " + str(e) + " for record: " + str(rec) + "\n")

    if startpos > endpos:
        endpos, startpos = startpos, endpos

    return startpos, endpos


def relevant(rec, vtype, ignorechroms):
    ''' Return true if a record matches the type of variant being investigated '''
    rel = (rec.is_snp and vtype == 'SNV') or (rec.is_sv and vtype == 'SV') or (rec.is_indel and vtype == 'INDEL')
    return rel and (ignorechroms is None or rec.CHROM not in ignorechroms)

def passfilter(rec, disabled=False):
    ''' Return true if a record is unfiltered or has 'PASS' in the filter field (pyvcf sets FILTER to None) '''
    if disabled:
        return True
    if rec.FILTER is None or rec.FILTER == '.' or not rec.FILTER:
        return True
    return False


def svmask(rec, vcfh, truchroms):
    ''' mask snv calls in sv regions '''
    if rec.is_snp and rec.CHROM in truchroms:
        for overlap_rec in vcfh.fetch(rec.CHROM, rec.POS-1, rec.POS):
            if overlap_rec.is_sv:
                    return True
    return False


def var_dist(v1, v2):
    """compute absolute distance between two variants
    """
    assert v1.CHROM == v2.CHROM
    return abs(v1.POS-v2.POS)


def get_close_matches(var, vcf_fh, win, indels_only=True):
    """Find close matches for variant (PyVCF record var) in file (PyVCF
    VCFReader vcf_fh) within given window win and return as list of
    item,distance tuples, sorted ascendingly by distance.
    
    """
    
    matches = list(vcf_fh.fetch(var.CHROM, var.POS-win, var.POS+1+win))
    if indels_only:
        matches = [m for m in matches if m.is_indel]
    if len(matches) == 0:
        return []
        
    dist_map = [(m, var_dist(m, var)) for m in matches]
    return sorted(dist_map, key=lambda x: x[1])


def have_identical_haplotypes(v1, v2, ref):
    """Check if two variant produce the same haplotype / variant sequence.

    - v1 and v2: PyVCF variants to compare
    - ref: PySAM FastaFile
    """

    assert (v1.is_indel or v1.is_snp) and (v2.is_indel or v2.is_snp)

    if v1.CHROM != v2.CHROM:
        return False

    if v1.is_snp and v2.is_snp:
        assert v1.REF.upper() == v2.REF.upper()
        return str(v1.ALT[0]).upper() == str(v2.ALT[0]).upper()
    if v1.is_snp or v2.is_snp:
        # one snp one indel: can't produce identical results
        return False

    assert v1.is_indel and v2.is_indel
    # only on allele per variant allowed
    assert len(v1.ALT) == 1 and len(v2.ALT) == 1, (
        "Can't handle multi-allelic entries")

    # get the sequence context whic fully overlaps both variants.
    # note: pyvcf is one-based, but start and end are zero-based half-open
    start = min([v1.POS, v2.POS])-1
    end = max([v1.POS + max([len(v1.REF), len(v1.ALT[0])]),
               v2.POS + max([len(v2.REF), len(v2.ALT[0])])
           ])
    chrom = v1.CHROM# made sure before they are identical before
    seq = list(ref.fetch(chrom, start, end).upper())

    if len(seq) != end-start:
        # FIXME how to handle
        LOG.critical("WARN: Couldn't fetch full sequence window. Skipping"
                     " allele-aware comparison, otherwise indices would"
                     " be off\n")
        raise NotImplementedError

    v1_offset = v1.POS-1-start
    v2_offset = v2.POS-1-start
    # lower() in replacement for debugging purposes only
    v1_seq = seq[:v1_offset] + list(str(v1.ALT[0]).lower()) + seq[v1_offset+len(v1.REF):]
    v2_seq = seq[:v2_offset] + list(str(v2.ALT[0]).lower()) + seq[v2_offset+len(v2.REF):]
    if False:
        print "reference sequence context\t%s" % (''.join(seq))
        print "v1 (offset %d) %s\t%s" % (v1_offset, v1, ''.join(v1_seq))
        print "v2 (offset %d) %s\t%s" % (v2_offset, v2, ''.join(v2_seq))
        print

    try:
        assert seq[v1_offset] == v1.REF[0].upper()
        assert seq[v2_offset] == v2.REF[0].upper()
        assert len(v1_seq) == len(seq) - len(v1.REF) + len(v1.ALT[0])
        assert len(v2_seq) == len(seq) - len(v2.REF) + len(v2.ALT[0])
    except AssertionError:
        #import pdb; pdb.set_trace()
        raise

    #if ''.join(v1_seq).upper() == ''.join(v2_seq).upper():
    #    print ''.join(v1_seq).upper()
    return ''.join(v1_seq).upper() == ''.join(v2_seq).upper()


def vars_identical(v1, v2):
    """Return if two variants are "identical" which here means, same
    chrom, same pos, same alt base
    
    FIXME check if pyvcf record cmp is the same
    """
    
    assert len(v1.ALT) == 1 and len(v2.ALT) == 1
    return all([v1.CHROM == v2.CHROM, v1.POS == v2.POS, v1.REF == v2.REF, v1.ALT[0] == v2.ALT[0]])


def write_var(var, fh):
    '''write pyvcf var to fh'''
    assert len(var.ALT)==1
    fh.write("%s\t%d\t.\t%s\t%s\t%s\t%s\t%s\n" % (
        var.CHROM, var.POS, var.REF,
        str(var.ALT[0]), var.QUAL if var.QUAL else ".",
        ';'.join(var.FILTER) if len(var.FILTER) else ".",
        string.translate(str(var.INFO), None, "{}'[] ").replace(":", "=").replace('=True', '')))
    
def evaluate(submission, truth, vtype='SNV', reffa=None, ignorechroms=None, ignorepass=False, 
             fpfh=False, tpfh=False, fnfh=False):
    ''' return stats on sensitivity, specificity, balanced accuracy '''

    assert vtype in ('SNV', 'SV', 'INDEL')
    subvcfh = vcf.Reader(filename=submission)
    truvcfh = vcf.Reader(filename=truth)

    if reffa:
        reffa_fh  = pysam.Fastafile(reffa)
    else:
        reffa_fh = None

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
            truchroms[trurec.CHROM] = True
            fns[str(trurec)] = trurec
    used_truth = {} # keep track of 'truth' sites used, they should only be usable once

    ''' parse submission vcf, compare to truth '''
    for subrec in subvcfh:
        if passfilter(subrec, disabled=ignorepass):
            if subrec.is_snp and vtype == 'SNV':
                if not svmask(subrec, truvcfh, truchroms):
                    subrecs += 1
            if subrec.is_sv and vtype == 'SV':
                subrecs += 1
            if subrec.is_indel and vtype == 'INDEL':
                subrecs += 1

        matched = False

        startpos, endpos = subrec.start, subrec.end

        if vtype == 'SV' and subrec.is_sv:
            startpos, endpos = expand_sv_ends(subrec)
        try:
            if relevant(subrec, vtype, ignorechroms) and passfilter(subrec, disabled=ignorepass) and subrec.CHROM in truchroms:
                if subrec.is_indel and reffa_fh:
                    for (trurec, d) in get_close_matches(subrec, truvcfh, 100, indels_only=True):
                        if vars_identical(trurec, subrec):
                            matched = True
                            break
                        elif have_identical_haplotypes(subrec, trurec, reffa_fh):
                            #sys.stderr.write("Rescueing %s which has same haplotype as %s\n" % (subrec, trurec))
                            matched = True
                            break
                    if matched:
                        used_truth[str(trurec)] = True
                else:
                    for trurec in truvcfh.fetch(subrec.CHROM, startpos, end=endpos):
                        if match(subrec, trurec, vtype=vtype) and str(trurec) not in used_truth:
                            matched = True
                            used_truth[str(trurec)] = True
                    
        except ValueError as e:
            sys.stderr.write("Warning: " + str(e) + "\n")

        if matched:
            if fns.has_key(str(trurec)):
                del fns[str(trurec)]
            tpcount += 1
            if tpfh:
                #import pdb; pdb.set_trace()
                assert not subrec.ID
                write_var(subrec, tpfh)
        else:
            if relevant(subrec, vtype, ignorechroms) and passfilter(subrec, disabled=ignorepass) and not svmask(subrec, truvcfh, truchroms): 
                fpcount += 1 # FP counting method needs to change for real tumors
                if fpfh:
                    #import pdb; pdb.set_trace()
                    assert not subrec.ID
                    #assert not subrec.FILTER
                    write_var(subrec, fpfh)

    if fnfh:
        for fn in fns.values():
            assert len(fn.ALT)==1
            write_var(fn, fnfh)
            
            
    print "tpcount, fpcount, subrecs, trurecs:"
    print tpcount, fpcount, subrecs, trurecs

    recall    = float(tpcount) / float(trurecs)
    if tpcount + fpcount > 0:
        precision = float(tpcount) / float(tpcount + fpcount)
    else:
        precision = 0
    fdr       = 1.0 - float(fpcount) / float(subrecs)
    f1score   = 0.0 if tpcount == 0 else 2.0*(precision*recall)/(precision+recall)

    return precision, recall, f1score

def main(args):

    chromlist = None
    if args.chromlist is not None:
        chromlist = args.chromlist.split(',')

    if not args.subvcf.endswith('.vcf') and not args.subvcf.endswith('.vcf.gz'):
        sys.stderr.write("submission VCF filename does not end in .vcf or .vcf.gz\n")
        sys.exit(1)

    if not os.path.exists(args.truvcf + '.tbi'):
        sys.stderr.write("truth VCF does not appear to be indexed. bgzip + tabix index required.\n")
        sys.exit(1)

    if args.mutype not in ('SV', 'SNV', 'INDEL'):
        sys.stderr.write("-m/--mutype must be either SV, SNV, or INDEL\n")
        sys.exit(1)

    if args.print_fp_to:
        if args.print_fp_to == "-":
            fpfh = sys.stdout
        else:
            fpfh = open(args.print_fp_to, 'w')
        print "Writing FPs to %s" % args.print_fp_to
    else:
        fpfh = False
    if args.print_tp_to:
        if args.print_tp_to == "-":
            tpfh = sys.stdout
        else:
            tpfh = open(args.print_tp_to, 'w')
        print "Writing TPs to %s" % args.print_tp_to
    else:
        tpfh = False
    if args.print_fn_to:
        if args.print_fn_to == "-":
            fnfh = sys.stdout
        else:
            fnfh = open(args.print_fn_to, 'w')
        print "Writing FNs to %s" % args.print_fn_to
    else:
        fnfh = False

    result = evaluate(args.subvcf, args.truvcf, vtype=args.mutype, reffa=args.reffa,
                      ignorechroms=chromlist, ignorepass=args.nonpass, 
                      fpfh=fpfh, tpfh=tpfh, fnfh=fnfh)

    print "precision, recall, F1 score: " + ','.join(map(str, result))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="check vcf output against a 'truth' vcf")
    parser.add_argument('-f',  '--ref', dest='reffa', help="Reference fasta file (enables haplotype-ware indel comparison and avoid problems due to positional amibiguity)")
    parser.add_argument('-v',  '--vcf', dest='subvcf', required=True, help="VCF being submitted for evaluation")
    parser.add_argument('-t',  '--truth', dest='truvcf', required=True, help="'Truth' VCF containing true positives")
    parser.add_argument('-m', '--mutype', dest='mutype', required=True, help="Mutation type: must be either SNV, SV, or INDEL")
    parser.add_argument('--ignore', dest='chromlist', default=None, help="(optional) comma-seperated list of chromosomes to ignore")
    parser.add_argument('--nonpass', dest='nonpass', action="store_true", help="evaluate all records (not just PASS records) in VCF")
    parser.add_argument('--print-fp-to', dest='print_fp_to', help="output false positive positions to this file")
    parser.add_argument('--print-tp-to', dest='print_tp_to', help="output true positive positions to this file")
    parser.add_argument('--print-fn-to', dest='print_fn_to', help="output false negatives positions to this file")
    args = parser.parse_args()
    main(args)

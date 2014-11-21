#!/usr/bin/env python

import sys, os
import vcf
import argparse
import string

'''
Evaluate VCFs against BAMSurgeon "Truth" VCFs
Adam Ewing, adam.ewing@mater.uq.edu.au
Requires PyVCF (https://github.com/jamescasbon/PyVCF)
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


def evaluate(submission, truth, vtype='SNV', ignorechroms=None, ignorepass=False, 
             fpfh=False, tpfh=False):
    ''' return stats on sensitivity, specificity, balanced accuracy '''

    assert vtype in ('SNV', 'SV', 'INDEL')
    subvcfh = vcf.Reader(filename=submission)
    truvcfh = vcf.Reader(filename=truth)

    tpcount = 0
    fpcount = 0
    subrecs = 0
    trurecs = 0

    truchroms = {}

    ''' count records in truth vcf, track contigs/chromosomes '''
    for trurec in truvcfh:
        if relevant(trurec, vtype, ignorechroms):
            trurecs += 1
            truchroms[trurec.CHROM] = True

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
                for trurec in truvcfh.fetch(subrec.CHROM, startpos, end=endpos):
                    if match(subrec, trurec, vtype=vtype) and str(trurec) not in used_truth:
                        matched = True
                        used_truth[str(trurec)] = True

        except ValueError as e:
            sys.stderr.write("Warning: " + str(e) + "\n")

        if matched:
            tpcount += 1
            if tpfh:
                #import pdb; pdb.set_trace()
                assert not subrec.ID
                assert not subrec.FILTER
                assert len(subrec.ALT)==1
                tpfh.write("%s\t%d\t.\t%s\t%s\t%s\t.\t%s\n" % (
                    subrec.CHROM, subrec.POS, subrec.REF,
                    str(subrec.ALT[0]), subrec.QUAL if subrec.QUAL else ".",
                    string.translate(str(subrec.INFO), None, "{}'[] ").replace(":", "=").replace('=True', '')))
        else:
            if relevant(subrec, vtype, ignorechroms) and passfilter(subrec, disabled=ignorepass) and not svmask(subrec, truvcfh, truchroms): 
                fpcount += 1 # FP counting method needs to change for real tumors
                if fpfh:
                    #import pdb; pdb.set_trace()
                    assert not subrec.ID
                    assert not subrec.FILTER
                    assert len(subrec.ALT)==1
                    fpfh.write("%s\t%d\t.\t%s\t%s\t%s\t.\t%s\n" % (
                        subrec.CHROM, subrec.POS, subrec.REF,
                        str(subrec.ALT[0]), subrec.QUAL if subrec.QUAL else ".",
                        string.translate(str(subrec.INFO), None, "{}'[] ").replace(":", "=").replace('=True', '')))

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
    result = evaluate(args.subvcf, args.truvcf, vtype=args.mutype, ignorechroms=chromlist, ignorepass=args.nonpass, 
                      fpfh=fpfh, tpfh=tpfh)

    print "precision, recall, F1 score: " + ','.join(map(str, result))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="check vcf output against a 'truth' vcf")
    parser.add_argument('-v',  '--vcf',    dest='subvcf', required=True, help="VCF being submitted for evaluation")
    parser.add_argument('-t',  '--truth',  dest='truvcf', required=True, help="'Truth' VCF containing true positives")
    parser.add_argument('-m', '--mutype', dest='mutype', required=True, help="Mutation type: must be either SNV, SV, or INDEL")
    parser.add_argument('--ignore', dest='chromlist', default=None, help="(optional) comma-seperated list of chromosomes to ignore")
    parser.add_argument('--nonpass', dest='nonpass', action="store_true", help="evaluate all records (not just PASS records) in VCF")
    parser.add_argument('--print-fp-to', dest='print_fp_to', help="output false positive positions to this file")
    parser.add_argument('--print-tp-to', dest='print_tp_to', help="output true positive positions to this file")
    args = parser.parse_args()
    main(args)

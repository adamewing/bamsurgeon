#!/usr/bin/env python

'''
One implementation of the read-editing spike-in, shared by addsnv and addindel.

Their makemut() bodies were the same ~120 lines apart from about nine
cosmetic differences, and both converged on mutation.mutate(), which already
took is_snv/is_insertion/is_deletion flags. That call was the seam; this
module is what sits on top of it.

The addsnv version is the one kept, because it was the more featured of the
two: copy-number adjustment of the requested VAF, --minmutreads forcing,
read regrouping so that mates are picked together, and the post-realignment
--coverdiff check. addindel had all of that except the haplotype grouping.

A haplotype cluster generalises: an indel is a cluster of one.
'''

import os
import random
import logging

import pysam

import bamsurgeon.mutation as mutation

from collections import defaultdict as dd
from dataclasses import dataclass
from uuid import uuid4

from bamsurgeon.aligners import remap_bam
from bamsurgeon.common import get_avg_coverage
from bamsurgeon.records import MutationRecord

FORMAT = '%(levelname)s %(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


BASES = ('A', 'T', 'C', 'G')


@dataclass
class VariantSite:
    ''' one requested edit; a cluster of these shares a haplotype '''
    chrom: str
    start: int
    end: int
    vaf: float = None
    kind: str = 'SNV'      # SNV, INS or DEL
    altbase: str = None    # SNV: explicit ALT, otherwise chosen at random
    insseq: str = None     # INS: sequence to insert

    @property
    def is_snv(self):
        return self.kind == 'SNV'


def mut(base, altbase):
    """ change base to something different
    """
    base = base.upper()
    if base not in BASES or (altbase is not None and altbase not in BASES):
        raise ValueError("ERROR base passed to mut(): " + str(base) + " not one of (A,T,C,G)\n")

    if altbase is not None:
        return altbase

    alt = base
    while alt == base:
        alt = BASES[int(random.uniform(0, 4))]
    return alt


def get_mutstr(site, reffile):
    '''
    Human-readable description of the edit, for the log.

    addindel's version of this returned the literal string 'FIX get_mutstr'.
    '''
    if site.is_snv:
        return 'SNV'
    if site.kind == 'INS':
        return 'INS:' + site.insseq
    return 'DEL:%d' % (site.end - site.start)


def resolve_snv_alt(sites, chrom, mutpos_list, reffile, mutid_list, ignoreref):
    '''
    Choose the ALT base for each site in the cluster.

    Returns (mutbase_list, refbase_list) or None if the cluster must be
    dropped.
    '''
    mutbase_list = []
    refbase_list = []

    for n, mutpos in enumerate(mutpos_list):
        refbase = reffile.fetch(chrom, mutpos-1, mutpos)
        altbase = sites[n].altbase
        refbase_list.append(refbase)

        if altbase == refbase.upper() and not ignoreref:
            logger.warning("%s specified ALT base matches reference, skipping mutation" % mutid_list[n])
            return None, None

        try:
            mutbase_list.append(mut(refbase, altbase))
        except ValueError as e:
            logger.warning("%s skipped site: %s %d %d due to N base: %s" %
                           (mutid_list[n], chrom, sites[n].start, sites[n].end, e))
            return None, None

    return mutbase_list, refbase_list


def resolve_indel(site):
    ''' the mutate() kwargs describing one indel '''
    is_insertion = site.kind == 'INS'
    return {
        'is_insertion': is_insertion,
        'is_deletion': not is_insertion,
        'ins_seq': site.insseq,
        'indel_start': site.start,
        'indel_end': site.end,
    }


def _small_record(site, chrom, pos, reffile, mutbase, vaf, dpr, mutid):
    ''' MutationRecord with explicit REF/ALT, as VCF wants for small variants '''
    if site.is_snv:
        ref = reffile.fetch(chrom, pos-1, pos).upper()
        return MutationRecord(kind='SNV', chrom=chrom, pos=pos, end=pos,
                              ref_allele=ref, alt_allele=mutbase,
                              vaf=vaf, dpr=dpr, mutid=mutid)

    if site.kind == 'INS':
        anchor = reffile.fetch(chrom, site.start-1, site.start).upper()
        return MutationRecord(kind='INDEL', chrom=chrom, pos=site.start,
                              end=site.start, svlen=len(site.insseq),
                              ref_allele=anchor,
                              alt_allele=anchor + site.insseq.upper(),
                              vaf=vaf, dpr=dpr, mutid=mutid)

    ref = reffile.fetch(chrom, site.start-1, site.end).upper()
    return MutationRecord(kind='INDEL', chrom=chrom, pos=site.start,
                          end=site.end, svlen=site.end - site.start,
                          ref_allele=ref, alt_allele=ref[0],
                          vaf=vaf, dpr=dpr, mutid=mutid)


def makemut(args, sites, avoid, alignopts):
    '''
    Apply one haplotype cluster of edits and return (tmpbams, records).

    Every site in a cluster shares one set of reads, so they must be applied
    together or they overwrite each other during the merge.
    '''
    if not sites:
        return None, []

    kinds = set(s.kind for s in sites)
    if len(kinds) > 1:
        # Composing makeins/makedel coordinate shifts with SNV query-position
        # offsets is a research problem, not an implementation detail. Refuse
        # rather than emit quietly-wrong reads.
        logger.error('cannot mix %s in one haplotype cluster at %s:%d; '
                     'run them as separate spike-ins'
                     % ('/'.join(sorted(kinds)), sites[0].chrom, sites[0].start))
        return None, []

    is_snv_cluster = sites[0].is_snv

    chrom = sites[0].chrom
    for site in sites:
        assert chrom == site.chrom, "haplotype clusters cannot span multiple chromosomes!"

    if args.seed is not None:
        random.seed(int(args.seed) + int(sites[0].start))

    mutid_list = ['%s_%s_%s_%s_%s' % (s.chrom, s.start, s.end, s.vaf, s.altbase)
                  for s in sites]

    bamfile = pysam.AlignmentFile(args.bamFileName, reference_filename=args.refFasta)
    bammate = pysam.AlignmentFile(args.bamFileName, reference_filename=args.refFasta) # use for mates to avoid iterator problems
    reffile = pysam.Fastafile(args.refFasta)

    tmpbams = []
    records = []

    vaf = None
    mutpos_list = []

    for site in sites:
        if vaf is None:
            vaf = site.vaf
        elif site.vaf is not None and vaf != site.vaf:
            logger.warning("multiple VAFs for single haplotype, using first encountered VAF: %f" % vaf)

        if is_snv_cluster:
            # position is drawn from within the requested interval
            mutpos_list.append(int(random.uniform(site.start, site.end+1)))
        else:
            mutpos_list.append(site.start)

    mutbase_list = None
    mut_kwargs = {}

    if is_snv_cluster:
        mutbase_list, _ = resolve_snv_alt(sites, chrom, mutpos_list, reffile,
                                          mutid_list, args.ignoreref)
        if mutbase_list is None:
            return None, []
        mut_kwargs = {'is_snv': True, 'mutbase_list': mutbase_list}
    else:
        mut_kwargs = resolve_indel(sites[0])

    hapstr = "_".join(('haplo', chrom, str(min(mutpos_list)), str(max(mutpos_list))))
    logdir = args.logdir
    log = open(os.path.join(logdir, os.path.basename(args.outBamFile) + "." + hapstr + ".log"), 'w')

    tmpoutbamname = os.path.join(args.tmpdir, hapstr + ".tmpbam." + str(uuid4()) + ".bam")
    logger.info("%s creating tmp bam: %s" % (hapstr, tmpoutbamname))
    outbam_muts = pysam.AlignmentFile(tmpoutbamname, 'wb', template=bamfile)

    # indels span start..end; SNVs are points inside the cluster
    if is_snv_cluster:
        mutstart, mutend = min(mutpos_list), max(mutpos_list) + 1
    else:
        site = sites[0]
        del_ln = (site.end - site.start) if site.kind == 'DEL' else 0
        mutstart, mutend = site.start, site.start + del_ln + 1

    res = mutation.mutate(args, log, bamfile, bammate, chrom, mutstart, mutend,
                          mutpos_list, avoid=avoid, mutid_list=mutid_list,
                          reffile=reffile, **mut_kwargs)

    if res.failed:
        outbam_muts.close()
        os.remove(tmpoutbamname)
        return None, []

    # pick reads to change
    readlist = [extqname for extqname, read in res.outreads.items()
                if read.seq != res.mutreads[extqname]]

    logger.info("%s len(readlist): %d" % (hapstr, len(readlist)))
    readlist.sort()
    random.shuffle(readlist)

    if len(readlist) < int(args.mindepth):
        logger.warning("%s too few reads in region (%d) skipping..." % (hapstr, len(readlist)))
        outbam_muts.close()
        os.remove(tmpoutbamname)
        return None, []

    if vaf is None:
        vaf = float(args.mutfrac) # default minor allele freq if not otherwise specified

    if args.cnvfile:
        cnv = pysam.Tabixfile(args.cnvfile, 'r')
        if chrom in cnv.contigs:
            for cnregion in cnv.fetch(chrom, min(mutpos_list), max(mutpos_list)+1):
                cn = float(cnregion.strip().split()[3]) # expect chrom,start,end,CN
                logger.info("%s copy number in snv region: %s %d %d = %f"
                            % (hapstr, chrom, min(mutpos_list), max(mutpos_list), cn))
                vaf = vaf / cn if cn > 0.0 else 0.0
                logger.info("%s adjusted VAF: %f" % (hapstr, vaf))
    else:
        logger.info("%s selected VAF: %f" % (hapstr, vaf))

    lastread = int(len(readlist)*vaf)

    # pick at least args.minmutreads if possible
    if lastread < int(args.minmutreads):
        if len(readlist) > int(args.minmutreads):
            lastread = int(args.minmutreads)
            logger.warning("%s forced %d reads." % (hapstr, lastread))
        else:
            logger.warning("%s dropped site with fewer reads than --minmutreads" % hapstr)
            outbam_muts.close()
            os.remove(tmpoutbamname)
            return None, []
    elif lastread <= 0:
        logger.warning("%s dropped mutation with 0 reads" % hapstr)
        outbam_muts.close()
        os.remove(tmpoutbamname)
        return None, []

    # regroup by original read name so both instances of a name travel together
    readtrack = dd(list)
    for readname in readlist:
        orig_name, readpos, pairend = readname.split(',')
        readtrack[orig_name].append('%s,%s' % (readpos, pairend))

    usedreads = 0
    newreadlist = []

    for orig_name in readtrack:
        for read_instance in readtrack[orig_name]:
            newreadlist.append(orig_name + ',' + read_instance)
            usedreads += 1

        if usedreads >= lastread:
            break

    readlist = newreadlist

    logger.info("%s picked: %d" % (hapstr, len(readlist)))

    wrote = 0
    nmut = 0
    mut_out = {}

    # change reads from .bam to mutated sequences
    for extqname, read in res.outreads.items():
        if read.seq != res.mutreads[extqname]:
            if not args.nomut and extqname in readlist:
                qual = read.qual # changing seq resets qual (see pysam API docs)
                read.seq = res.mutreads[extqname] # make mutation
                read.qual = qual
                nmut += 1
        if (not res.has_snp) or args.force:
            wrote += 1
            mut_out[extqname] = read

    muts_written = {}

    for extqname in mut_out:
        if extqname in muts_written:
            continue

        outbam_muts.write(mut_out[extqname])
        muts_written[extqname] = True

        if res.mutmates[extqname] is None:
            continue

        # is mate also in mutated list?
        mate_read = res.mutmates[extqname]

        pairname = 'F' # read is first in pair
        if mate_read.is_read2:
            pairname = 'S' # read is second in pair
        if not mate_read.is_paired:
            pairname = 'U' # read is unpaired

        mateqname = ','.join((mate_read.qname, str(mate_read.pos), pairname))

        if mateqname in mut_out:
            # yes: output mutated mate
            outbam_muts.write(mut_out[mateqname])
            muts_written[mateqname] = True
        else:
            # no: output original mate
            outbam_muts.write(mate_read)

    logger.info("%s wrote: %d, mutated: %d" % (hapstr, wrote, nmut))

    if res.has_snp and not args.force:
        # nothing was written; the temp BAM used to be left behind here
        outbam_muts.close()
        os.remove(tmpoutbamname)
        logger.warning("%s dropped for proximity to SNP" % hapstr)
        log.close()
        bamfile.close()
        bammate.close()
        return None, []

    outbam_muts.close()

    remap_bam(args.aligner, tmpoutbamname, args.refFasta, alignopts,
              threads=int(args.alignerthreads), mutid=hapstr,
              paired=(not args.single), picardjar=args.picardjar)

    outbam_muts = pysam.AlignmentFile(tmpoutbamname, reference_filename=args.refFasta)
    coverwindow = 1

    # SNVs measure across the cluster; an indel measures across its own span.
    # These are the windows the two tools used separately -- widening either
    # of them shifts the reported DPR and the --coverdiff decision.
    incover_start = min(mutpos_list) - coverwindow
    if is_snv_cluster:
        incover_end = max(mutpos_list) + coverwindow
    else:
        site = sites[0]
        del_ln = (site.end - site.start) if site.kind == 'DEL' else 0
        incover_end = site.start + del_ln + coverwindow

    avgincover = get_avg_coverage(bamfile, chrom, incover_start, incover_end)
    avgoutcover = get_avg_coverage(outbam_muts, chrom, incover_start, incover_end)

    logger.info("%s avgincover: %f, avgoutcover: %f" % (hapstr, avgincover, avgoutcover))

    spikein_frac = float(nmut)/float(wrote) if wrote > 0 else 0.0

    # qc cutoff for final depth
    if not ((avgoutcover > 0 and avgincover > 0 and
             avgoutcover/avgincover >= float(args.coverdiff)) or args.force):
        outbam_muts.close()
        os.remove(tmpoutbamname)
        if os.path.exists(tmpoutbamname + '.bai'):
            os.remove(tmpoutbamname + '.bai')
        logger.warning("%s dropped for outcover/incover < %s" % (hapstr, str(args.coverdiff)))
        log.close()
        bamfile.close()
        bammate.close()
        return None, []

    tmpbams.append(tmpoutbamname)

    for n, site in enumerate(sites):
        mutbase = mutbase_list[n] if mutbase_list else None
        records.append(_small_record(site, chrom, mutpos_list[n], reffile,
                                     mutbase, spikein_frac, avgoutcover,
                                     mutid_list[n]))
        log.write("\t".join((site.kind, chrom, str(mutpos_list[n]),
                             get_mutstr(site, reffile), str(avgincover),
                             str(avgoutcover), str(spikein_frac),
                             str(res.maxfrac))) + "\n")

    outbam_muts.close()
    bamfile.close()
    bammate.close()
    log.close()

    return tmpbams, records


def cluster_sites(targets, hapsize):
    '''
    Group sites that are close enough to share reads.

    Sites in the same cluster are mutated together; sites in different
    clusters are mutated independently and merged, so two clusters touching
    the same reads will overwrite each other. That is what --haplosize
    controls.
    '''
    clusters = []
    current = []
    lastchrom = None
    laststart = None

    for target in targets:
        if lastchrom is None:
            lastchrom = target.chrom
            laststart = target.start
            current.append(target)

        elif target.chrom == lastchrom:
            if laststart is None:
                laststart = target.start
                current.append(target)
            elif target.start - laststart < hapsize:
                current.append(target)
            else:
                clusters.append(current)
                current = [target]

        else:
            clusters.append(current)
            current = [target]
            lastchrom = target.chrom
            laststart = None

    clusters.append(current)

    return [c for c in clusters if c]


def run_spikein(args, clusters, toolname):
    '''
    Shared driver: farm clusters out to the pool, merge, replace, write VCF.

    addsnv and addindel had a copy of this each, differing only in the log
    directory name and the temp file prefix.
    '''
    import sys
    import bamsurgeon.replace_reads as rr
    import bamsurgeon.makevcf as makevcf

    from shutil import move
    from concurrent.futures import ProcessPoolExecutor
    from bamsurgeon.common import mergebams, dictlist

    logger.info("starting %s called with args: %s" % (sys.argv[0], ' '.join(sys.argv)))

    if (args.bamFileName.endswith('.bam') and not os.path.exists(args.bamFileName + '.bai')) or \
       (args.bamFileName.endswith('.cram') and not os.path.exists(args.bamFileName + '.crai')):
        logger.error("input file must be indexed, not .bai or .crai file found for %s" % args.bamFileName)
        sys.exit(1)

    alignopts = {}
    if args.alignopts is not None:
        alignopts = dict([o.split(':') for o in args.alignopts.split(',')])

    avoid = dictlist(args.avoidreads) if args.avoidreads is not None else None

    args.logdir = '%s_logs_%s' % (toolname, os.path.basename(args.outBamFile))

    for d in (args.tmpdir, args.logdir):
        if not os.path.exists(d):
            os.mkdir(d)
            logger.info("created directory: %s" % d)

    assert os.path.exists(args.logdir), "could not create output directory!"
    assert os.path.exists(args.tmpdir), "could not create temporary directory!"

    outbam_mutsfile = "%s.%s.muts.bam" % (toolname, str(uuid4()))
    with pysam.AlignmentFile(args.bamFileName, reference_filename=args.refFasta) as bamfile:
        pysam.AlignmentFile(outbam_mutsfile, 'wb', template=bamfile).close()

    pool = ProcessPoolExecutor(max_workers=int(args.procs))
    results = [pool.submit(makemut, args, cluster, avoid, alignopts)
               for cluster in clusters]

    tmpbams = []
    records = []

    for result in results:
        tmpbamlist, recs = result.result()
        if tmpbamlist:
            for tmpbam in tmpbamlist:
                if os.path.exists(tmpbam):
                    tmpbams.append(tmpbam)
            records.extend(recs)

    if len(tmpbams) == 0:
        logger.error("no succesful mutations")
        sys.exit(1)

    tmpbams.sort()

    if len(tmpbams) == 1:
        move(tmpbams[0], outbam_mutsfile)
    else:
        mergebams(tmpbams, outbam_mutsfile, maxopen=int(args.maxopen))

    for bam in tmpbams:
        if os.path.exists(bam):
            os.remove(bam)
        if os.path.exists(bam + '.bai'):
            os.remove(bam + '.bai')

    if args.skipmerge:
        logger.info("skipping merge, plase merge reads from %s manually." % outbam_mutsfile)
    else:
        if args.tagreads:
            from bamsurgeon.markreads import markreads
            tmp_tag_bam = 'tag.%s.bam' % str(uuid4())
            markreads(outbam_mutsfile, args.refFasta, tmp_tag_bam)
            move(tmp_tag_bam, outbam_mutsfile)
            logger.info("tagged reads.")

        logger.info("done making mutations, merging mutations into %s --> %s"
                    % (args.bamFileName, args.outBamFile))
        rr.replace_reads(args.bamFileName, outbam_mutsfile, args.outBamFile,
                         args.refFasta, keepqual=True, seed=args.seed)
        os.remove(outbam_mutsfile)

    var_basename = '.'.join(os.path.basename(args.varFileName).split('.')[:-1])
    bam_basename = '.'.join(os.path.basename(args.outBamFile).split('.')[:-1])

    vcf_fn = args.vcf + bam_basename + '.' + toolname + '.' + var_basename + '.vcf'

    makevcf.write_vcf(records, args.refFasta, vcf_fn)

    logger.info('vcf output written to ' + vcf_fn)


def read_snv_targets(path, maxmuts, haplosize):
    ''' BED-like: chrom start end [vaf] [altbase] '''
    from operator import attrgetter

    targets = []
    ntried = 0

    with open(path, 'r') as bedfile:
        for line in bedfile:
            if maxmuts and ntried >= maxmuts:
                break
            c = line.strip().split()
            if not c:
                continue

            site = VariantSite(chrom=c[0], start=int(c[1]), end=int(c[2]), kind='SNV')

            if len(c) > 3:
                site.vaf = float(c[3])
            if len(c) == 5:
                altbase = c[4].upper()
                assert altbase in BASES, "ERROR: ALT " + altbase + " not A, T, C, or G!\n"
                site.altbase = altbase

            targets.append(site)
            ntried += 1

    targets.sort(key=attrgetter('chrom', 'start'))

    return cluster_sites(targets, haplosize)


def read_indel_targets(path, maxmuts):
    '''
    BED-like: chrom start end vaf INS|DEL [seq]

    Indels are not grouped: each line becomes a cluster of one, preserving
    file order, which is what addindel did before.
    '''
    clusters = []
    ntried = 0

    with open(path, 'r') as bedfile:
        for line in bedfile:
            if maxmuts and ntried >= maxmuts:
                break
            c = line.strip().split()
            if not c:
                continue

            kind = c[4]
            assert kind in ('INS', 'DEL'), 'indel type must be INS or DEL: %s' % kind

            site = VariantSite(chrom=c[0], start=int(c[1]), end=int(c[2]),
                               vaf=float(c[3]), kind=kind,
                               insseq=c[5] if kind == 'INS' else None)

            clusters.append([site])
            ntried += 1

    return clusters

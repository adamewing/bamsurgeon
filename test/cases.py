'''
The spike-in test matrix, one entry per case that test/*.sh used to cover.

`expect` says what success looks like:

  'variants'     the run succeeds and every emitted site validates
  'no_variants'  the run correctly refuses, because its thresholds exclude
                 every site; exiting nonzero is the pass condition
  'skipmerge'    the run succeeds but leaves an unmerged donor BAM, so there
                 is no output BAM to validate

`needs` lists executables the case cannot run without. Cases naming an
aligner that is not installed are skipped rather than failed.
'''

from dataclasses import dataclass


@dataclass
class Case:
    name: str
    tool: str                  # addsnv | addindel | addsv
    varfile: str
    bam: str = 'testregion_realign.bam'
    aligner: str = 'mem'
    args: tuple = ()
    expect: str = 'variants'
    needs: tuple = ('bwa', 'samtools')
    # some cases legitimately lose sites; require this fraction to validate
    min_pass_rate: float = 1.0
    # a known limitation: the case should work and does not. Runs anyway, so
    # that pytest reports XPASS if it is ever fixed.
    xfail: str = None


ONT = dict(bam='testregion_ont.bam', aligner='minimap2',
           needs=('minimap2', 'samtools'))

# --coverdiff 0.9 (the addsnv default) drops every long-read site, because
# minimap2 remapping of ONT reads loses more coverage than that allows
ONT_SNV_ARGS = ('--alignopts', 'x:map-ont', '--ignorepileup', '--single',
                '--mindepth', '4', '--coverdiff', '0.1')
ONT_INDEL_ARGS = ('--alignopts', 'x:map-ont', '--ignorepileup', '--single',
                  '--mindepth', '4')


CASES = [
    # ---- SNVs -------------------------------------------------------------
    Case('snv', 'addsnv', 'random_snvs.txt',
         args=('-n', '5', '-c', 'test_cnvlist.txt.gz'),
         # two of the first five sites are 1bp apart and land in separate
         # haplotype clusters at the default -z 0, so one is overwritten
         min_pass_rate=0.8),

    Case('snv_haplosize_auto', 'addsnv', 'random_snvs.txt',
         args=('-n', '5', '-c', 'test_cnvlist.txt.gz', '-z', 'auto')),

    Case('snv_alt', 'addsnv', 'snv_alttest.txt',
         args=('-n', '10', '-c', 'test_cnvlist.txt.gz', '-p', '2')),

    Case('snv_minmutreads', 'addsnv', 'random_snvs.txt',
         args=('-n', '10', '-c', 'test_cnvlist.txt.gz', '-p', '2',
               '--minmutreads', '25'),
         # includes the same colliding pair, 1 of 10 sites
         min_pass_rate=0.9),

    Case('snv_maxdepth', 'addsnv', 'random_snvs.txt',
         args=('-n', '10', '--maxdepth', '5'),
         expect='no_variants'),

    Case('snv_skipmerge', 'addsnv', 'random_snvs.txt',
         args=('-n', '5', '-c', 'test_cnvlist.txt.gz', '--skipmerge'),
         expect='skipmerge'),

    Case('snv_ont', 'addsnv', 'test_ont_snvs.txt',
         args=('-n', '5',) + ONT_SNV_ARGS, **ONT),

    # ---- indels -----------------------------------------------------------
    Case('indel', 'addindel', 'test_indels.txt'),

    Case('indel_minmutreads', 'addindel', 'test_indels.txt',
         args=('-c', 'test_cnvlist.txt.gz', '-p', '2', '--minmutreads', '50')),

    Case('indel_maxdepth', 'addindel', 'test_indels.txt',
         args=('-c', 'test_cnvlist.txt.gz', '-p', '2', '--maxdepth', '20'),
         expect='no_variants'),

    Case('indel_ins_ont', 'addindel', 'test_ont_indel_ins.txt',
         args=('-n', '5',) + ONT_INDEL_ARGS, **ONT),

    Case('indel_del_ont', 'addindel', 'test_ont_indel_del.txt',
         args=('-n', '5',) + ONT_INDEL_ARGS, **ONT,
         xfail='deleting 100bp from 27-45kb ONT reads collapses coverage over '
               'the site from ~26x to ~1x on remapping, so every site is '
               'dropped by --coverdiff. The insertion equivalent works.'),

    # ---- structural variants ---------------------------------------------
    Case('sv', 'addsv', 'test_sv.txt',
         args=('-p', '4', '--keepsecondary', '--inslib', 'test_inslib.fa')),

    Case('sv_inslib', 'addsv', 'test_sv_inslib.txt',
         args=('-p', '4', '--keepsecondary', '--inslib', 'test_inslib.fa')),

    Case('sv_empty_cols', 'addsv', 'test_sv_empty_cols.txt',
         args=('-p', '4', '--keepsecondary')),

    Case('sv_simerr', 'addsv', 'test_sv.txt',
         args=('-p', '4', '--keepsecondary', '--inslib', 'test_inslib.fa',
               '--simerr', '0.05'),
         # simulated sequencing error costs some breakend support
         min_pass_rate=0.75),

    Case('big_sv', 'addsv', 'test_big_sv.txt',
         args=('-p', '4', '--keepsecondary',
               '--donorbam', 'testregion_realign.bam')),

    # ---- aligners not installed here -------------------------------------
    Case('snv_novoalign', 'addsnv', 'random_snvs.txt', aligner='novoalign',
         args=('-n', '5'), needs=('novoalign', 'samtools')),
    Case('snv_gsnap', 'addsnv', 'random_snvs.txt', aligner='gsnap',
         args=('-n', '5'), needs=('gsnap', 'samtools')),
    Case('snv_tmap', 'addsnv', 'random_snvs.txt', aligner='tmap',
         bam='testregion_tmap.bam', args=('-n', '5', '--single'),
         needs=('tmap', 'samtools')),
    Case('indel_novoalign', 'addindel', 'test_indels.txt', aligner='novoalign',
         needs=('novoalign', 'samtools')),
    Case('sv_novoalign', 'addsv', 'test_sv.txt', aligner='novoalign',
         args=('-p', '1', '--keepsecondary'), needs=('novoalign', 'samtools')),
]


# Cases from test/*.sh with no counterpart here, and why:
#
#   test_trn.sh            test_data/test_trn.txt is not in the repo
#   test_snv_haplo.sh      test_data/random_snvs_haplopairs.txt is not in the repo
#   test_snv_bowtie2.sh    test_data/testregion_bt2.bam is not in the repo
#   test_indel_bowtie2.sh  same
#   test_sv_exact.sh       --require_exact no longer exists; breakpoints are
#                          exact by construction, which test_spikein asserts
#   test_snv_avoid.sh      replaced by test_avoidreads below, which does not
#                          use one file as both input and output
#   test_snv_environ.sh    tested the $BAMSURGEON_PICARD_JAR fallback; picard
#                          is no longer used at all

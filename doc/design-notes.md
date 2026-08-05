# Design notes

Decisions taken during the SV-engine and read-editing rework, and the work
still outstanding. Kept in the repo because most of this is reasoning that is
expensive to reconstruct from the diff alone.

## Decisions taken

**Chained (`;`-delimited) compound mutations are not supported.** Few real
use cases, and the same result is reachable with serial spike-ins. Lines
containing `;` are reported and skipped rather than partially applied. No
fixture used them.

**`BIG*` is not a separate implementation.** DEL, DUP and INV are handled
identically at any size. `BIGDEL`/`BIGDUP`/`BIGINV` are accepted as deprecated
spellings. The 10 kbp promotion and 5 kbp demotion thresholds are gone. Note
that `BIGDUP`'s VAF used to be read only when the line had six fields, so
`BIGDUP 1.0` silently fell back to `--svfrac`; it is honoured now.

**Duplication interior is a choice, not a size class.** Simulated by default;
taken from `--donorbam` when one is given, which keeps the template small for
large duplications and gives the interior a real error profile.

**Mixed SNV and indel in one haplotype cluster is refused.** Composing
`makeins`/`makedel` coordinate shifts with SNV query-position offsets is a
research problem, not an implementation detail.

**`--ignoresnps` and `--insane` are `--force` split in two.** `--force`
continues to bypass both the SNP-proximity guard and the coverage check;
the other two flags bypass one each.

**The unified entry point applies SVs first, then small variants on the SV
output.** Not an implementation convenience: the SV engine replaces every read
over its interval with reads simulated from a reference template, which
carries no small variants, so anything spiked in beforehand is reverted. The
other order works because the read-editing engine happily edits simulated
reads, and it means a small variant inside a duplication or inversion
survives rather than needing the composition machinery described below. A
small variant inside a *deletion* still has nothing to sit on, which is
correct.

The alternative considered, and rejected, was merging both donors and calling
`replace_reads` once. One pass instead of two, but it destroys exactly the
small variants the ordering saves.

**`--coverdiff` stays per variant class.** addsnv defaulted to 0.9 and
addindel to 0.1, and that is tuning rather than drift: a deletion removes
bases and its reads may soft-clip, so the input/output coverage ratio
legitimately falls further. Forcing one value showed it immediately -- a 10bp
deletion that passes at 0.1 is dropped at 0.9. The unified CLI therefore
takes `--coverdiff` (0.9) and `--indel-coverdiff` (0.1). Tools that set only
one keep using it for both.

**`-s` is retired in the unified CLI only.** It means `--snvfrac` in
addsnv/addindel and `--svfrac` in addsv, so one letter cannot carry both
where the engines meet. Retiring it in the legacy shims would break every
existing invocation for no benefit, since within one tool it was never
ambiguous.

## Stage 4 — done

`varinput.py` gives both engines one `VariantRequest`, produced by a VCF
reader or by the legacy varfile readers (relocated verbatim). The unified
`bamsurgeon` entry point takes one VCF containing every class. Proximity and
overlap warnings, `-z auto`, `--seed type=int` everywhere and addsv's missing
`--maxopen` all landed.

Remaining from the original plan, deliberately not done: nothing. The
`--coverdiff` and `-s` items were resolved as described above rather than as
originally sketched.

### SVs overwrite small variants (background)

The ordering above avoids this in a single unified run. It still applies to
anyone running the tools separately in the wrong order, and to the general
composition problem.

`get_reads()` excludes every read in an SV region and `replace_reads(allreads=True)`
substitutes wgsim reads simulated from the template, which knows nothing about
variants spiked in earlier. Under the old assembler a prior SNV sometimes
survived into the contig depending on whether `-cov_cutoff auto` popped it;
under a reference template it is deterministically reverted.

The reference template makes the common case tractable — the template is a
known interval at known offsets, so known small variants can be applied to it
with the existing `MutableSeq` ops before simulating. What stays hard is the
modelling, not the plumbing: conditional VAF (an SNV at 0.3 inside an SV at
0.5 — what fraction of simulated SV reads carry it?) and phase (SNV on the SV
allele or the reference allele; VCF cannot say without phasing).

**This is what the unified CLI does automatically.** Running the tools by hand,
the same rule applies: SVs first, then SNVs/indels on the SV output BAM.

## Stage 5 — done

- ~~Port `scripts/evaluator.py` and `scripts/bamregions_from_vcf.py` off
  PyVCF onto `pysam.VariantFile`.~~ Done. `test/test_evaluator.py` covers it
  and needs only pysam, so it runs without an aligner or picard.

  Not mechanical in the end: SV mode had never worked. `expand_sv_ends()`
  did `int(rec.INFO.get('END')[0])`, which subscripts a scalar, and the
  `except TypeError` meant to report that referenced a `logger` the module
  never defined, so every SV evaluation died with a NameError. Separately,
  the record that had matched was tracked by following the fetch loop
  variable, so the wrong truth site was marked used and deleted from the
  false-negative set whenever a window returned more than one record.
- ~~Convert `test/*.sh` to a pytest table asserting via
  `validate_spikein.py`.~~ Done. `test/cases.py` is the matrix,
  `test/test_spikein.py` runs it, and the shell scripts are gone. Each case
  gets its own tmp directory, which the shell versions did not: they all wrote
  into `test_data/` and some only worked if another had run first.

  Four cases were dropped because their fixtures are not in the repo
  (`test_trn.txt`, `random_snvs_haplopairs.txt`, `testregion_bt2.bam`), and
  `test_sv_exact.sh` because `--require_exact` no longer exists. The reason
  for each is recorded in `test/cases.py`.

  `test_sv_coordinates_match_the_request` is the permanent guard on the Stage
  2 property: output POS equals the requested start and END the requested end,
  for every non-insertion SV.

### Known limitation: deletions in long reads

`indel_del_ont` is marked xfail. Deleting 100bp from 27-45kb ONT reads
collapses coverage over the site from ~26x to ~1x on remapping, so every site
is dropped by `--coverdiff` at any threshold. The insertion equivalent
(`indel_ins_ont`) works, so this is specific to `makedel` on long reads rather
than to long reads generally. It runs rather than being skipped, so pytest
reports XPASS if it is ever fixed.
- ~~**Remove picard, and with it java.**~~ Done. `bamtofastq()` is pysam-native:
  pairs are collected in a dict keyed by qname and written when the mate
  arrives, reverse-strand reads are complemented back to original orientation
  with their quality strings reversed alongside, and secondary and
  supplementary alignments are excluded — which is what
  `INCLUDE_NON_PRIMARY_ALIGNMENTS=false` did. Names carry the `/1` and `/2`
  suffixes picard appended, so read names in realigned output are unchanged.

  Checked against picard 2.27.3 on a slice of the test BAM: the two FASTQs
  are identical record-for-record once sorted. Record *order* differs, so
  realigned BAMs are not byte-identical to earlier runs — golden checksum
  comparison across this change does not transfer, and the pytest suite is
  what validates it.

  Unpaired reads (a mate excluded by the region slice, or by a filter) are
  dropped with a warning, which is what picard did silently.

  `--picardjar` and `$BAMSURGEON_PICARD_JAR` are still accepted so existing
  invocations do not break; the flag is hidden from `--help` and warns that it
  has no effect. `check_java()` and `default-jre` are gone from the dependency
  check, setup script and Dockerfile.

## Scale: what is untested, and where it will bite first

Everything here has only been exercised on a ~770 kb chr22 region at ~50x.
It needs stress-testing on whole genomes before it can be trusted at that
scale. Several choices deliberately traded headroom for simplicity, and these
are the ones to look at first.

**`get_reads()` is no longer streaming.** Exact-count selection has to know the
whole candidate population before it can emit anything, so it now holds a set
of qnames for the interval and makes two passes. Memory is O(reads in
interval) where it used to be O(1), and `select_qnames()` sorts that set. Fine
for a 14 kb interval; a chromosome-arm deletion at 30x is a different
proposition.

**Interval templates fetch the whole span as a Python string.**
`build_interval()` fetches `ref[start-pad : end+pad]` and lets `MutableSeq`
slice it. A 1 Mb deletion is a 1 MB string; a chromosome-arm deletion is
~100 MB, and the deletion case only ever needs the two flanks. Worth
special-casing DEL to fetch flanks only if large deletions become common.

**Duplication defaults do not scale.** The simulated interior is
`(ndups + 1)` copies of the interval, so a 1 Mb DUP with `ndups 3` is a 4 MB
template and proportionally many simulated reads. `--donorbam` exists exactly
for this and is the WGS path; the simulate-by-default choice suits test-scale
intervals.

**`replace_reads()` loads the entire donor BAM into a dict** keyed by qname,
holding read objects. This predates the rework and is untouched, but with
thousands of mutations merged into one donor it is the most likely place to
run out of memory. It is also the one step that sees every read in the target
BAM.

**`mergebams` in the SV path still inherits `maxopen=100`** rather than 1000,
because addsv never grew a `--maxopen`. With thousands of temp BAMs that is
many more hierarchical merge rounds than necessary. Listed under Stage 4.

**The validator does a fetch or pileup per record.** Acceptable for a
truth VCF of tens of thousands of sites, but it has never been run against a
whole-genome truth set and makes no attempt to batch by region.

## Fixture limitations worth knowing

`test_data/testregion_realign.bam` has no column with more than one distinct
base over the regions examined — no sequencing errors, no heterozygous sites.
The SNP-proximity guard therefore cannot fire on the shipped fixtures at any
`--snvfrac`, which is how several bugs in it survived: a per-column reset of
`hasSNP`, an off-by-one in `countBaseAtPos`, and a `TypeError` in its warning
path. Any test intended to exercise `--snvfrac`, `--ignoresnps` or the linked-SNP
logic needs a fixture with real allelic variation.

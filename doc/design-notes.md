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

## Outstanding: Stage 4 — VCF input and unified CLI

- `varinput.py`: one `VariantRequest` produced by both a VCF reader and the
  frozen legacy varfile readers. Legacy readers move verbatim; they are the
  compatibility contract.
- **Proximity warning.** Warn when two requested sites fall in different
  haplotype clusters but within a read length of each other, comparing target
  *intervals* rather than the drawn position. Warn harder when their VAFs
  differ, because clustering collapses a cluster to a single VAF.
- **`-z auto`** as an opt-in value alongside the integer. Deliberately not the
  default: auto-expanding clusters would trade a visible failure (a dropped
  mutation, detectable as BAM/VCF divergence) for an invisible one (a silently
  overridden VAF).
- **SV / small-variant overlap detection.** An SNV inside an SV's excluded read
  set is destroyed. Detection and warning only; composition is out of scope,
  see below.
- CLI collisions: `-s` means `--snvfrac` in two tools and `--svfrac` in the
  third; `--coverdiff` defaults 0.9 vs 0.1; `--seed` is `type=int` only in
  addsv; addsv has no `--maxopen` so it inherits `mergebams`'s default of 100
  rather than 1000.

### SVs overwrite small variants

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

**Ordering workaround, valid today:** SVs first, then SNVs/indels on the SV
output BAM. `addsnv` mutates whatever reads it finds, including simulated ones.
It is SNV-then-SV that reverts.

## Outstanding: Stage 5 — cleanup

- Port `scripts/evaluator.py` and `scripts/bamregions_from_vcf.py` off PyVCF
  (`import vcf`; unmaintained since 2018) onto `pysam.VariantFile`. Keep
  `have_identical_haplotypes()`.
- Convert `test/*.sh` to a pytest table asserting via `validate_spikein.py`.
  Four cases reference fixtures absent from the repo (`test_trn.sh`,
  `test_snv_haplo.sh`, `test_snv_bowtie2.sh`, `test_indel_bowtie2.sh`);
  `test_snv_avoid.sh` uses one file as both input and output and depends on
  `test_snv.sh` having run first.
- **Remove picard, and with it java.** `bamtofastq()` in `common.py` is
  picard's only use anywhere in the codebase; `check_java()` and the
  `default-jre` package exist solely to support it. Doing this drops a ~200 MB
  JRE from the Dockerfile and setup script.

  Do it pysam-native rather than shelling to `samtools fastq`: the
  per-mutation temp BAMs are small, so collecting pairs in a dict keyed by
  qname and emitting in first-encounter order is deterministic, avoids a fork,
  and works on CRAM. What has to be right is reverse-complementing
  reverse-strand reads back to original orientation, reversing their quality
  strings with them, and excluding secondary and supplementary alignments —
  which is what `INCLUDE_NON_PRIMARY_ALIGNMENTS=false` does today.

  Expect FASTQ record order to differ from picard's, so realignment output
  changes byte-wise. Not a correctness change, but it invalidates golden
  checksum comparison against earlier runs; validate with
  `validate_spikein.py` before and after instead.

## Fixture limitations worth knowing

`test_data/testregion_realign.bam` has no column with more than one distinct
base over the regions examined — no sequencing errors, no heterozygous sites.
The SNP-proximity guard therefore cannot fire on the shipped fixtures at any
`--snvfrac`, which is how several bugs in it survived: a per-column reset of
`hasSNP`, an off-by-one in `countBaseAtPos`, and a `TypeError` in its warning
path. Any test intended to exercise `--snvfrac`, `--ignoresnps` or the linked-SNP
logic needs a fixture with real allelic variation.

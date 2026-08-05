# BAMSurgeon

BAMSurgeon adds simulated mutations to existing BAM files. It is built for
benchmarking variant callers: you take a real alignment, spike in mutations you
chose, and you now have a dataset whose right answer you know exactly.

Everything is driven from one command and one VCF:

```
bamsurgeon -v variants.vcf -f input.bam -r reference.fasta -o mutated.bam
```

which writes the mutated BAM and a truth VCF describing what should have been modified.

---

## Contents

- [How it works](#how-it-works)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Input](#input)
- [Command reference](#command-reference)
- [Output](#output)
- [Validating a spike-in](#validating-a-spike-in)
- [Evaluating a variant caller](#evaluating-a-variant-caller)
- [The single-class tools](#the-single-class-tools)
- [Notes and limitations](#notes-and-limitations)
- [Testing](#testing)

---

## How it works

There are two engines, and which one handles a variant depends on how you
wrote it.

**The read-editing engine** handles SNVs and small indels. It pulls the reads
over the site, edits the bases in place, realigns them, and merges them back.
The reads stay the sample's own reads, with their own error profile.

**The SV engine** handles deletions, duplications, inversions, insertions and
breakends. It builds the mutated haplotype directly from reference sequence,
simulates reads over it with `wgsim`, and swaps them for the reads it removed.
Breakpoints are exact by construction: the output VCF's `POS` and `END` are
the ones you asked for, because the engine never has to work out where
anything landed.

### Ordering

A single run applies **SVs first, then small variants on the SV output**. This
is not an implementation convenience. The SV engine replaces every read over
its interval with reads simulated from a reference template, which carries no
small variants — so anything spiked in beforehand is reverted. The other order
works because the read-editing engine is happy to edit simulated reads, and it
means an SNV inside a duplication or an inversion survives.

An SNV inside a *deletion* still has nothing to sit on, which is correct.
Overlaps are reported when the input is read, so you will see them.

---

## Requirements

| | |
|---|---|
| Python | 3.6+ with `pysam` |
| Always | `samtools` |
| Realignment | `bwa`, or another supported aligner |
| Structural variants | `wgsim` (ships with samtools) |
| Optional | `bcftools`, `tabix` for the surrounding workflow |

No Java, and no picard — BAM to FASTQ conversion is done in pysam.

`scripts/setup_env.sh` installs all of it on a Debian/Ubuntu box and is safe
to re-run. `scripts/check_dependencies.py` reports what is missing.

The reference FASTA needs both indexes:

```
samtools faidx reference.fasta
bwa index reference.fasta
```

---

## Installation

```
git clone https://github.com/adamewing/bamsurgeon.git
cd bamsurgeon
pip install -e .
```

That puts `bamsurgeon` on your `$PATH`. Editable installs need pip 21.3 or
newer and setuptools 61 or newer; if pip reports *"File setup.py not found,
directory cannot be installed in editable mode"*, upgrade them:

```
pip install --upgrade pip setuptools
```

You can also run from the source tree without installing:

```
PYTHONPATH=bin python3 -m bamsurgeon.cli.main --help
```

---

## Quick start

Two SNVs, an indel and an inversion, in one file:

```
##fileformat=VCFv4.2
##contig=<ID=22,length=51304566>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele fraction">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
22	33714965	.	A	G	.	PASS	VAF=0.5
22	33769483	.	T	A	.	PASS	VAF=0.25
22	33944542	.	GA	G	.	PASS	VAF=0.5
22	34071043	.	C	<INV>	.	PASS	SVTYPE=INV;END=34084754;VAF=0.5
```

Run it:

```
bamsurgeon \
  -v variants.vcf \
  -f input.bam \
  -r reference.fasta \
  -o mutated.bam \
  -p 8 --seed 1234
```

Then sort, index, and check the result:

```
samtools sort -o mutated.sorted.bam mutated.bam
samtools index mutated.sorted.bam

python3 scripts/validate_spikein.py \
  --vcf mutated.bamsurgeon.variants.vcf \
  --bam mutated.sorted.bam \
  --orig-bam input.bam \
  --ref reference.fasta
```

---

## Input

The unified entry point takes **a VCF**. A file is treated as a VCF if its
name ends in `.vcf`, `.vcf.gz` or `.bcf`, or if it starts with
`##fileformat=VCF`. Legacy varfiles work only with the single-class tools —
they cannot express more than one variant class, so combined input needs a
VCF.

### The routing rule

How you write the ALT decides which engine gets it:

| ALT | Engine | Becomes |
|---|---|---|
| single base, e.g. `A` → `G` | read-editing | SNV |
| ALT longer than REF | read-editing | insertion of the extra bases |
| REF longer than ALT | read-editing | deletion of the extra bases |
| symbolic, e.g. `<DEL>` `<DUP>` `<INV>` `<INS>` | SV | that SV |
| bracket notation, e.g. `C[22:34314981[` | SV | breakend |
| any ALT with `SVTYPE` in INFO | SV | that SVTYPE |

MNPs and complex alleles (REF and ALT the same length but longer than 1) are
skipped with a warning. `<CNV>` and other unsupported `SVTYPE`s are skipped
with a warning.

### Allele fraction

Taken from the first of these that is present:

1. `INFO/VAF`
2. `FORMAT/AF`, if the VCF has exactly one sample
3. `--vaf` on the command line (default 0.5)

With `-c/--cnvfile`, the VAF is divided by the local absolute copy number.

### Structural variant details

Anything the VCF spec covers is read from where the spec puts it:

| Field | Meaning |
|---|---|
| `END` | end of the interval, for `DEL`, `DUP`, `INV` |
| `SVLEN` | used to derive `END` when `END` is absent |
| `NDUPS` | additional tandem copies for `DUP` (default 1, so two copies total) |
| `TSDLEN` | target site duplication length for `INS` |
| `MATEID` | de-duplicates the two records of a breakend pair |

Simulation directives with no natural VCF home use a `BS_` prefix:

| Key | Meaning |
|---|---|
| `BS_INSSEQ` | the literal sequence to insert |
| `BS_INSLIB` | name of an entry in the `--inslib` FASTA |
| `BS_MOTIF` | preferred cut site for an insertion, written `AA^TTTT` |

An insertion needs one of `BS_INSSEQ` or `BS_INSLIB`. For compatibility with
older output, the record `ID` is used as a library entry name if neither is
present.

Breakends use standard bracket notation, and all four orientations are
supported. Only one record of a mate pair is needed — the other end is derived
from the ALT — but a pair sharing a `MATEID` is read correctly and applied
once.

```
22	33607508	bnd_a	C	C[22:34314981[	.	PASS	SVTYPE=BND;MATEID=bnd_b
```

### Legacy varfiles

The BED-like formats still work with `addsnv.py`, `addindel.py` and
`addsv.py`:

```
addsnv    chrom  start  end  [vaf]  [altbase]
addindel  chrom  start  end  vaf    INS|DEL  [seq]
addsv     chrom  start  end  TYPE   [type-specific fields...]
```

For `addsv`, the type-specific fields are:

| Type | Fields |
|---|---|
| `DEL`, `INV` | `[vaf]` |
| `DUP` | `[ndups] [vaf]` |
| `INS` | `seq\|file\|RND\|INSLIB:name [tsdlen] [motif] [vaf]` |
| `TRN` | `mate_chrom mate_pos mate_end orientation [vaf]` |

Splitting is on single whitespace characters, not runs, because empty columns
carry meaning: `DUP\t\t0.9` is "default copy count, VAF 0.9".

`BIGDEL`, `BIGDUP` and `BIGINV` are accepted as deprecated spellings of `DEL`,
`DUP` and `INV`. There is no separate large-variant implementation any more,
and no size threshold — the same code handles an event at any size.

---

## Command reference

### Required

| Option | Meaning |
|---|---|
| `-v, --varfile` | VCF describing the variants to add |
| `-f, --bamfile` | input SAM/BAM/CRAM |
| `-r, --reference` | reference FASTA, `samtools faidx` and `bwa index`ed |
| `-o, --outbam` | output BAM |

### General

| Option | Default | Meaning |
|---|---|---|
| `--vaf` | 0.5 | default allele fraction, overridden per record |
| `--snvfrac` | 1 | maximum linked-SNP MAF before a site is dropped |
| `-p, --procs` | 1 | parallel spike-ins |
| `--seed` | none | seed the RNG; required for reproducible output |
| `-c, --cnvfile` | none | tabix-indexed absolute copy number, adjusts VAF |
| `--aligner` | backtrack | `STAR`, `backtrack`, `bowtie2`, `bwakit`, `gsnap`, `mem`, `minimap2`, `novoalign`, `tmap` |
| `--alignopts` | none | `option1:value1,option2:value2,...` |
| `--alignerthreads` | 1 | threads per realignment |
| `--force` | off | ignore nearby SNPs and low coverage |
| `--single` | off | input is single-ended |
| `--tagreads` | off | add a `BS` tag to every altered read |
| `--maxopen` | 1000 | open file limit during merge |
| `--tmpdir` | `bamsurgeon.tmp` | scratch directory |
| `--debug` | off | keep intermediates |
| `--vcf` | `''` | path *prefix* for the truth VCF; end it in `/` for a directory |

### Small variants

| Option | Default | Meaning |
|---|---|---|
| `--mindepth` | 10 | minimum depth at the site |
| `--maxdepth` | 2000 | maximum depth at the site |
| `-d, --coverdiff` | 0.9 | allowed output/input depth ratio after realignment, SNVs |
| `--indel-coverdiff` | 0.1 | the same for indels |
| `-z, --haplosize` | 0 | group sites within this distance onto one haplotype; `auto` uses the sampled read length |
| `--minmutreads` | 3 | minimum mutated reads per site |
| `--avoidreads` | none | file of read names to leave alone |
| `--ignoreref` | off | mutate even when ALT matches the reference |
| `--ignoresnps` | off | skip the confounding-SNP scan |
| `--insane` | off | skip the coverage sanity check |
| `--ignorepileup` | off | do not require a pileup across the region |
| `--requirepaired` | off | skip sites with unpaired reads |
| `--nomut` | off | dry run |

`--coverdiff` is per class on purpose. A deletion removes bases and its reads
may soft-clip, so its input/output ratio legitimately falls further than an
SNV's. Forcing one value for both drops small deletions.

### Structural variants

| Option | Default | Meaning |
|---|---|---|
| `--sv-mindepth` | 10 | minimum depth at a breakend |
| `--sv-maxdepth` | 2000 | maximum depth at a breakend |
| `-l, --maxlibsize` | 600 | maximum fragment length of the library |
| `--pad` | 2× maxlibsize | reference flank either side of a breakpoint |
| `--readlen` | modal read length in the BAM | simulated read length |
| `--maxdfrac` | 0.1 | maximum discordant read fraction before a site is refused |
| `--ismean` | 300 | mean insert size |
| `--issd` | 70 | insert size standard deviation |
| `--simerr` | 0 | error rate for simulated reads |
| `--inslib` | none | FASTA of insertion sequences |
| `--donorbam` | none | BAM supplying duplication interiors |
| `--allowN` | off | tolerate N in the reference slice, replacing with A |
| `--keepsecondary` | off | keep secondary alignments in the output |

`--ismean` and `--issd` should be measured from your input BAM.
`samtools stats` reports both, as *insert size average* and *insert size
standard deviation*.

`--donorbam` matters for large duplications. Without it the duplicated
interior is simulated, which means a template of `(ndups + 1)` copies of the
interval — fine at test scale, wasteful at megabase scale. With it, the
interior comes from real reads with a real error profile.

### A note on `-s`

`-s` is not available in the unified CLI. It means `--snvfrac` in
`addsnv`/`addindel` and `--svfrac` in `addsv`, and one letter cannot carry
both where the engines meet. It still works in the single-class tools, where
it was never ambiguous.

---

## Output

### The BAM

**The output BAM is not sorted, despite what its header says.** It inherits
the input's `@HD SO:coordinate` but reads are appended out of order, so
`samtools index` will refuse it. Sort it first:

```
samtools sort -o mutated.sorted.bam mutated.bam
samtools index mutated.sorted.bam
```

With `--tagreads`, every altered read carries a `BS` tag, which
`scripts/remove_non_BS.py` can filter on.

### The truth VCF

Named from the output BAM and the input variant file:

```
<--vcf prefix><outbam basename>.bamsurgeon.<varfile basename>.vcf
```

so `-o mutated.bam -v variants.vcf --vcf results/` gives
`results/mutated.bamsurgeon.variants.vcf`.

It is sorted and tabix-indexable, it describes what actually landed rather
than what was requested, and **it is valid input to another run** — the `BS_`
extensions exist so that a truth VCF round-trips.

Small variants are written with `GT 0/1`, SVs with `GT ./.`. Every record
carries `SOMATIC` and the achieved `VAF`; SVs also carry `PRECISE` and `DPR`,
the depth in the spiked region after realignment.

---

## Validating a spike-in

`scripts/validate_spikein.py` checks a truth VCF against the BAM that
produced it, and exits nonzero on failure. It replaces eyeballing a pileup.

```
python3 scripts/validate_spikein.py \
  --vcf truth.vcf --bam mutated.sorted.bam --ref reference.fasta \
  --orig-bam input.bam \
  [--tsv report.tsv]
```

`--orig-bam` is optional but worth giving: it baselines SV interval depth
against the input, which is what makes the depth checks possible.

Per-site output looks like:

```
chrom  pos       kind   req_vaf  depth  alt  obs_vaf  status  note
22     33714965  SNV    0.500    54     27   0.500    OK
22     34271043  DUP    0.900    103    48   0.466    OK      depth ratio 3.42 of 3.70 predicted
```

| Status | Meaning |
|---|---|
| `OK` | found, at the right place, at roughly the right fraction |
| `LOW_VAF` | present but under-represented |
| `ABSENT` | no supporting reads |
| `NO_COVERAGE` | too few reads to say anything |
| `SHIFTED` | breakend found, but not where it was requested |

Thresholds are all flags — `--vaf-ratio-min`, `--sv-exact-tolerance`,
`--pass-rate`, `--sv-pass-rate` and others — with generous defaults, so a pass
means something and a failure is worth looking at.

SVs are not graded on VAF. The supporting fraction reported for one is clipped
reads over local depth, which is not an allele fraction: only reads physically
spanning a junction can clip, so even a VAF of 1.0 gives a fraction well below
1. Deletions and duplications are instead checked on interval depth against
the original BAM, which is a sound signal — and for a duplication the expected
ratio, `(1 - vaf) + vaf * (ndups + 1)`, is checked for magnitude and not only
direction.

---

## Evaluating a variant caller

`scripts/evaluator.py` compares a caller's VCF against a truth VCF and reports
precision, recall and F1.

```
python3 scripts/evaluator.py \
  -v submission.vcf -t truth.vcf.gz -m SNV \
  [-f reference.fasta] [--tp tp.vcf --fp fp.vcf --fn fn.vcf]
```

`-m` is `SNV`, `INDEL` or `SV`. The truth VCF must be bgzipped and tabix
indexed.

Give `-f` for indels. Left-alignment ambiguity in repeat context means two
records at different positions can describe the same haplotype; with the
reference, the comparison recognises that and stops reporting false negatives
that are not.

---

## The single-class tools

`addsnv.py`, `addindel.py` and `addsv.py` still exist and still take legacy
varfiles. They share the engines the unified CLI uses, so behaviour matches.
Reach for them when you have an existing varfile-based pipeline; use
`bamsurgeon` for anything new.

If you do run them separately, **apply SVs first and then spike small variants
into the SV output**, for the reason described under [Ordering](#ordering).

`--picardjar` and `$BAMSURGEON_PICARD_JAR` are still accepted everywhere and
warn that they have no effect.

---

## Notes and limitations

**Reproducibility.** Pass `--seed`. With it, output is byte-identical across
different `-p` values; without it, read selection is seeded from the system
RNG and will differ run to run.

**Long reads.** Supported, but they want different settings — realignment
loses more coverage, so the default `--coverdiff` of 0.9 drops everything:

```
bamsurgeon -v variants.vcf -f ont.bam -r ref.fasta -o out.bam \
  --aligner minimap2 --alignopts x:map-ont \
  --single --ignorepileup --coverdiff 0.1 --mindepth 4
```

**Copy number.** `--maxdepth` is an absolute cap, and amplified sequence
carries more reads legitimately. It is raised in proportion to the local copy
number, taken from `--cnvfile` and from any duplication applied earlier in the
same run. `--mindepth` is not scaled.

**Deletions and coverage.** The coverage sanity check measures depth at a
mutation's breakends, not across it. Measuring inside a deletion measures the
deletion against itself: a read that carries it contributes no aligned bases
there, so a homozygous deletion would score ~0 and be dropped for having
worked.

**Chained mutations.** `;`-delimited compound mutations are not supported.
Lines containing `;` are reported and skipped. The same result is reachable
with serial spike-ins.

**Mixed SNV and indel in one haplotype cluster** is refused. Composing indel
coordinate shifts with SNV query-position offsets is a research problem.

**Scale.** Everything here has been exercised on a ~770 kb region at ~50x.
Whole-genome behaviour is untested, and several choices traded headroom for
simplicity — see `doc/design-notes.md` for where it is most likely to bite
first.

---

## Testing

```
python3 -m pytest
```

`test/cases.py` is the matrix of spike-in configurations and
`test/test_spikein.py` runs each one, then checks the output BAM against the
truth VCF that run emitted. Each case gets its own temporary directory.

Cases naming an aligner that is not installed are skipped rather than failed,
so a run reporting only skips is not evidence of anything — check that the
passed count is what you expect.

`test/test_evaluator.py` and the unit tests in `test/test_copynumber.py` need
only pysam and run anywhere. `test/test_combined.py` exercises the unified
entry point on a VCF containing all three classes.

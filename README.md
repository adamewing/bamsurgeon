# Bamsurgeon<!-- omit in toc -->
[![DOI](https://zenodo.org/badge/4290471.svg)](https://zenodo.org/badge/latestdoi/4290471)

*tools for adding mutations to .bam files, used for testing mutation callers*

Please see doc/Manual.pdf for instructions and contact adam.ewing@gmail.com with questions.
  
**Before using simulation-based methods, consider whether there are validated benchmarking data available for your use case**. For example the [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) or the [PrecisionFDA Truth Challenge](https://precision.fda.gov/challenges/truth/results) may provide a means to test on non-simulated data.

## Table of contents<!-- omit in toc -->
- [Dependencies](#dependencies)
- [Getting started](#getting-started)

## Dependencies

The following presumes $HOME/bin is in your $PATH

Samtools / BCFTools / wgsim:

```
git clone --recurse-submodules --remote-submodules https://github.com/samtools/htslib.git
make -C htslib

git clone https://github.com/samtools/samtools.git
make -C samtools
cp samtools/samtools $HOME/bin
cp samtools/misc/wgsim $HOME/bin

git clone https://github.com/samtools/bcftools.git
make -C bcftools
cp bcftools/bcftools $HOME/bin
```

BWA (or another aligner, currently supported aligners include bwa, novoalign, gsnap, bowtie2, and tmap)

```
git clone https://github.com/lh3/bwa.git
make -C bwa
cp bwa/bwa $HOME/bin
```

Picard tools

```
wget https://github.com/broadinstitute/picard/releases/download/1.131/picard-tools-1.131.zip
unzip picard-tools-1.131.zip
```

Exonerate

```
git clone https://github.com/adamewing/exonerate.git
cd exonerate
git checkout v2.4.0
autoreconf -i
./configure && make && make check && make install
```

Velvet

```
wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz
tar xvzf velvet_1.2.10.tgz
make -C velvet_1.2.10
cp velvet_1.2.10/velvetg $HOME/bin
cp velvet_1.2.10/velveth $HOME/bin
```

Pysam

```
pip install pysam
```

## Getting started

It is not necessary to install BAMSurgeon. It works as any other standalone Python script. Use case example:

```
python3 bin/addsv.py -p 1 -v test_data/test_sv.txt -f test_data/testregion_realign.bam -r $REF -o test_data/testregion_sv_mut.bam --aligner mem --keepsecondary --seed 1234 --inslib test_data/test_inslib.fa
```

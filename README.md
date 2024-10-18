# Bamsurgeon<!-- omit in toc -->
[![DOI](https://zenodo.org/badge/4290471.svg)](https://zenodo.org/badge/latestdoi/4290471)

*Tools (addsnv.py, addindel.py and addsv.py) for adding genomic variants to SAM/BAM/CRAM files, used for testing variant callers*

Please see doc/Manual.pdf for instructions and contact adam.ewing@gmail.com with questions.
  
**Before using simulation-based methods, consider whether there are validated benchmarking data available for your use case**. For example the [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) or the [PrecisionFDA Truth Challenge](https://precision.fda.gov/challenges/truth/results) may provide a means to test on non-simulated data.

## Installation

BAMSurgeon can be installed via bioconda: [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/bamsurgeon/README.html) (thanks @kubranarci!)

## Dependencies

The following presumes $HOME/bin is in your $PATH

SAMtools / wgsim:

```
git clone --recurse-submodules --remote-submodules https://github.com/samtools/htslib.git
make -C htslib

git clone https://github.com/samtools/samtools.git
make -C samtools
cp samtools/samtools $HOME/bin
cp samtools/misc/wgsim $HOME/bin
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
git clone https://github.com/dzerbino/velvet.git
cd velvet
make
cp velvetg $HOME/bin
cp velveth $HOME/bin
```

Pysam

```
pip install pysam
```

To check if all dependencies are correctly installed, you may execute the following script:
```
python3 -O scripts/check_dependencies.py
```

## Getting started

Alternately, it is not necessary to install BAMSurgeon. It works as any other standalone Python script. Usage example:

```
python3 -O bin/addsv.py -p 1 -v test_data/test_sv.txt -f test_data/testregion_realign.bam -r test_data/Homo_sapiens_chr22_assembly19.fasta -o test_data/testregion_sv_mut.bam --aligner mem --keepsecondary --seed 1234 --inslib test_data/test_inslib.fa
```

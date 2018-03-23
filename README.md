
[![Build Status](https://travis-ci.org/adamewing/bamsurgeon.svg?branch=master)](https://travis-ci.org/adamewing/bamsurgeon)

[![GSR Certified](https://popmodels.cancercontrol.cancer.gov/gsr/static//img/gsr-certified.jpg)](https://popmodels.cancercontrol.cancer.gov/gsr/certification/)

## Bamsurgeon:
*tools for adding mutations to .bam files, used for testing mutation callers*

Please see doc/Manual.pdf for instructions and contact adam.ewing@gmail.com with questions.

# Installation

```
python setup.py build
python setup.py install
```

# External Prerequisites:

The following presumes $HOME/bin is in your $PATH

Samtools / BCFTools / wgsim:

```
git clone https://github.com/samtools/htslib.git
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
pip install cython
pip install pysam
```

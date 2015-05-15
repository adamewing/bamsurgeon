FROM ubuntu:14.04
MAINTAINER Adam Ewing <adam.ewing@gmail.com>

WORKDIR ~/
RUN apt-get update && apt-get install -y software-properties-common

RUN add-apt-repository -y ppa:scipy/ppa

RUN apt-get install -y \
    python \
    python-numpy \
    python-scipy \
    python-pip \
    zlib1g-dev \
    git \
    wget \
    libncurses5-dev \
    unzip

RUN pip install cython
RUN pip install pysam

RUN mkdir $HOME/bin
RUN export PATH=$PATH:$HOME/bin

RUN wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz
RUN tar xvzf velvet_1.2.10.tgz
RUN make -C velvet_1.2.10
RUN cp velvet_1.2.10/velvetg $HOME/bin
RUN cp velvet_1.2.10/velveth $HOME/bin

RUN git clone https://github.com/lh3/bwa.git
RUN make -C bwa
RUN cp bwa/bwa $HOME/bin

RUN git clone https://github.com/samtools/htslib.git
RUN make -C htslib

RUN git clone https://github.com/samtools/samtools.git
RUN make -C samtools
RUN cp samtools/samtools $HOME/bin
RUN cp samtools/misc/wgsim $HOME/bin

RUN git clone https://github.com/samtools/bcftools.git
RUN make -C bcftools
RUN cp bcftools/bcftools $HOME/bin

RUN wget https://github.com/broadinstitute/picard/releases/download/1.131/picard-tools-1.131.zip
RUN unzip picard-tools-1.131.zip

RUN wget https://www.ebi.ac.uk/~guy/exonerate/exonerate-2.2.0-x86_64.tar.gz
RUN tar xvzf exonerate-2.2.0-x86_64.tar.gz
RUN cp exonerate-2.2.0-x86_64/bin/* $HOME/bin

RUN git clone https://github.com/adamewing/bamsurgeon.git


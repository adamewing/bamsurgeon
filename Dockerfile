FROM ubuntu:20.04
LABEL author="Adam Ewing <adam.ewing@gmail.com>"

ENV PATH=$PATH:$HOME/bin
ARG DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    python3 \
    python3-pip \
    git \
    wget \
    build-essential \
    libz-dev \
    libglib2.0-dev \
    libbz2-dev \
    liblzma-dev \
    default-jre \
    autoconf \
    samtools \
    bwa


RUN mkdir $HOME/bin

RUN wget https://github.com/dzerbino/velvet/archive/refs/tags/v1.2.10.tar.gz && tar -xvzf v1.2.10.tar.gz
RUN make -C velvet-1.2.10
RUN cp velvet-1.2.10/velvetg $HOME/bin && cp velvet-1.2.10/velveth $HOME/bin

RUN git clone https://github.com/adamewing/exonerate.git
RUN cd exonerate && autoreconf -fi  && ./configure && make && make install

RUN wget https://github.com/broadinstitute/picard/releases/download/2.27.3/picard.jar
RUN chmod +x picard.jar
RUN export BAMSURGEON_PICARD_JAR=$HOME/picard.jar

RUN pip install pysam

RUN git clone https://github.com/adamewing/bamsurgeon.git
RUN export PATH=$PATH:$HOME/bin && cd bamsurgeon

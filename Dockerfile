FROM ubuntu:19.04
MAINTAINER Adam Ewing <adam.ewing@gmail.com>

ENV PATH=$PATH:$HOME/bin

WORKDIR ~/

#install the bareminimum and remove the cache
RUN apt-get update && apt-get install -y --no-install-recommends \
    python \
    python-dev \
    python-numpy \
    python-scipy \
    python-pip \
    python-setuptools \
    python-wheel \
    zlib1g-dev \
    libbz2-dev \
    git \
    wget \
    libncurses5-dev \
    liblzma-dev \
    pkg-config \
    automake \
    autoconf \
    gcc \
    libglib2.0-dev \
    default-jre \
    samtools \
    bcftools \
    bwa \
    && rm -rf /var/lib/apt/lists/*


RUN mkdir $HOME/bin

RUN wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz && tar -xvzf velvet_1.2.10.tgz
RUN make -C velvet_1.2.10
RUN cp velvet_1.2.10/velvetg $HOME/bin && cp velvet_1.2.10/velveth $HOME/bin

RUN git clone https://github.com/adamewing/exonerate.git
RUN cd exonerate && autoreconf -fi  && ./configure && make && make install

#have to do it with the & because pysam needs cython
RUN pip install cython && pip install pysam==0.12.0

RUN git clone https://github.com/adamewing/bamsurgeon.git
RUN export PATH=$PATH:$HOME/bin && cd bamsurgeon && python setup.py install

CMD []

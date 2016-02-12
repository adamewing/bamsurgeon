#cd ..
#git clone https://github.com/samtools/htslib.git
#sudo make -C htslib
#git clone https://github.com/samtools/samtools.git
#sudo make -C samtools
#sudo cp samtools/samtools $HOME/bin
#sudo cp samtools/misc/wgsim $HOME/bin
#git clone https://github.com/samtools/bcftools.git
#sudo make -C bcftools
#sudo cp bcftools/bcftools $HOME/bin
#git clone https://github.com/lh3/bwa.git
#sudo make -C bwa
#sudo cp bwa/bwa $HOME/bin
mkdir bin
WORKDIR=`pwd`
wget https://github.com/broadinstitute/picard/releases/download/1.131/picard-tools-1.131.zip
unzip picard-tools-1.131.zip
wget http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-2.2.0-x86_64.tar.gz
tar -xvf exonerate-2.2.0-x86_64.tar.gz
sudo cp -rf exonerate-2.2.0-x86_64/ $WORKDIR/bin/
wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz
tar xvzf velvet_1.2.10.tgz
sudo make -C velvet_1.2.10
sudo cp velvet_1.2.10/velvetg $WORKDIR/bin/
sudo cp velvet_1.2.10/velveth $WORKDIR/bin/
sudo pip install --install-option="--prefix=$WORKDIR/bin" cython
sudo pip install --install-option="--prefix=$WORKDIR/bin" pysam
wget ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Assembly/GRCh37-HG19_Broad_variant/Homo_sapiens_assembly19.fasta
git clone https://github.com/lh3/bwa.git
sudo make -C bwa
sudo cp bwa/bwa $WORKDIR/bin
git clone https://github.com/samtools/htslib.git
sudo make -C htslib
git clone https://github.com/samtools/samtools.git
sudo make -C samtools
sudo cp samtools/samtools $WORKDIR/bin/
sudo cp samtools/misc/wgsim $WORKDIR/bin/
git clone https://github.com/samtools/bcftools.git
sudo make -C bcftools
sudo cp bcftools/bcftools $WORKDIR/bin/
export PATH=$WORKDIR/bin/:$PATH


#cd bamsurgeon


#bundle install --path .bundle --quiet --without=development

exit 0

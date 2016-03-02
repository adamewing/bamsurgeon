#wget ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Assembly/GRCh37-HG19_Broad_variant/Homo_sapiens_assembly19.fasta
#git clone https://github.com/lh3/bwa.git
#sudo make -C bwa
#sudo cp bwa/bwa $HOME/bin
bwa index Homo_sapiens_assembly19.fasta
samtools faidx Homo_sapiens_assembly19.fasta

exit 0

#!/bin/bash
 
echo -e "\n\n### Get NCBI RefSeq microbial FASTA sequence files ###"
mkdir -p ./blastdb &&
cd ./blastdb &&
wget ftp://ftp.ncbi.nih.gov/refseq/release/microbial/microbial.*.genomic.fna.gz &&
#wget ftp://ftp.ncbi.nih.gov/refseq/release/microbial/microbial.?.*genomic.fna.gz &&
gunzip -f -v *.gz &&
cat *.fna > microbial_all.fa &&

echo -e "\n\n### Make NCBI BLAST database ###"
makeblastdb -in microbial_all.fa -dbtype nucl -logfile blastdbmake.log &&

ls -alh . &&
cd ..

echo -e "\n### Done!"

# EOF

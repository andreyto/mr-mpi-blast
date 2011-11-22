#!/bin/bash

if [ $MRMPIBLAST_PREFIX == "" ]; then
    echo "The environment variable, MRMPIBLAST_PREFIX is not set"
    exit
fi

TUTORIAL_PREFIX=$MRMPIBLAST_PREFIX/tutorial

echo -e "\n\n### Get NCBI RefSeq microbial FASTA sequence files ###"
mkdir -p $TUTORIAL_PREFIX/blastdb &&
cd $TUTORIAL_PREFIX/blastdb &&
wget ftp://ftp.ncbi.nih.gov/refseq/release/microbial/microbial.*.genomic.fna.gz &&
#wget ftp://ftp.ncbi.nih.gov/refseq/release/microbial/microbial.?.*genomic.fna.gz &&
gunzip -f -v *.gz &&
cat *.fna > microbial_all.fa &&

echo -e "\n\n### Make NCBI BLAST database ###"
$MRMPIBLAST_PREFIX/bin/makeblastdb -in microbial_all.fa -dbtype nucl -logfile blastdbmake.log &&

ls -alh $TUTORIAL_PREFIX/blastdb &&
cd $TUTORIAL_PREFIX

echo -e "\n### Done!"

# EOF

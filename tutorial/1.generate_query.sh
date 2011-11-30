#!/bin/bash

if [ -z "$MRMPIBLAST_PREFIX" ]; then
    echo "The environment variable, MRMPIBLAST_PREFIX is not set"
    exit
fi

TUTORIAL_PREFIX=$MRMPIBLAST_PREFIX/tutorial

echo -e  "\n\n### Get NCBI RefSeq viral FASTA sequence files ###"
mkdir -p $TUTORIAL_PREFIX/query &&
cd $TUTORIAL_PREFIX/query &&
wget ftp://ftp.ncbi.nih.gov/refseq/release/viral/*.fna.gz &&
gunzip -f -v *.gz &&
cat *.fna > viral_all.fa &&

echo -e "\n\n### Generate query from the viral FASTA"
python $TUTORIAL_PREFIX/tools/splitter.py viral_all.fa 1000 100 viral_all_query.fa &&

echo -e "\n\n### Indexing query FASTA file for mr-mpi-blast"
python $TUTORIAL_PREFIX/tools/seqindexer.py -i viral_all_query.fa -o viral_all_query.fa.idx -d viral_all_query.fa.def -u 0 -s 0 -b 1 &&

ls -alh $TUTORIAL_PREFIX/query
cd $TUTORIAL_PREFIX

echo -e "\n### Done!"

# EOF


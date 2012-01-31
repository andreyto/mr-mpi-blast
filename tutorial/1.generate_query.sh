#!/bin/bash

echo -e  "\n\n### Get NCBI RefSeq viral FASTA sequence files ###"
mkdir -p ./query &&
cd ./query &&
wget ftp://ftp.ncbi.nih.gov/refseq/release/viral/*.fna.gz &&
gunzip -f -v *.gz &&
cat *.fna > viral_all.fa &&

echo -e "\n\n### Generate queries from the viral FASTA"
splitter.py viral_all.fa 1000 100 viral_all_query.fa &&

echo -e "\n\n### Indexing query FASTA file"
seqindexer.py -i viral_all_query.fa -o viral_all_query.fa.idx -d viral_all_query.fa.def -s 1 -b 1 &&

ls -alh . &&
cd ..

echo -e "\nDone!"

# EOF


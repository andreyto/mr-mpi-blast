#!/bin/bash

echo -e "\n### Indexing query FASTA file for mr-mpi-blast"

## Option descriptions
# -s 1: the generated sequential number will be starting from 1
# -b 1: save full deflines in .def file
mkdir -p ./query &&
cd ./query &&
seqindexer.py -i 100_simul_seqs.fa -o 100_simul_seqs.fa.idx -d 100_simul_seqs.fa.def -s 1 -b 1 &&
cd ..
echo -e "\nDone!"

# EOF


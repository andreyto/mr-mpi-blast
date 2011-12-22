#!/bin/bash

echo -e "\n### Indexing query FASTA file for mr-mpi-blast"

## Option descriptions
# -u 0: use generated sequential number as unique query id  
# -s 1: the generated sequential number will be starting from 1
# -b 1: save full deflines in .def file
seqindexer.py -i 100_simul_seqs.fa -o 100_simul_seqs.fa.idx -d 100_simul_seqs.fa.def -u 0 -s 1 -b 1


echo -e "\nDone!"

# EOF


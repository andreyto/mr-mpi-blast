#!/bin/bash

echo -e "\n### Indexing query FASTA file for mr-mpi-blast"

## Option descriptions
# -u 1: use gi in the input file as unique query id  
# -b 0: save full deflines in .def file
seqindexer.py -i 30_real_seqs.fa -o 30_real_seqs.fa.idx -d 30_real_seqs.fa.def -u 1 -b 0


echo -e "\nDone!"

# EOF


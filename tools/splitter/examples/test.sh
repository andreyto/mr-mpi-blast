#!/bin/bash

echo -e "\n\n### Generate queries from a FASTA file"
../splitter.py 30_real_seqs.fa 1000 100 out.fa 

echo -e "Done"

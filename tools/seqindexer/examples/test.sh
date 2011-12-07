#!/bin/bash 

python ../seqindexer.py -i test.query -o test.query.idx -d test.query.def -u 0 -s 0 -b 1 &&
python ../seqindexer.py -i 100.query -o 100.query.idx -d 100.query.def -u 0 -s 0 -b 1 &&
python ../seqindexer.py -i 10seqs.query -o 10seqs.query.idx -d 10seqs.query.def -u 1 -b 1 

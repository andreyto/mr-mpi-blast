#!/bin/bash 

mkdir -p generated
python ../seqindexer.py -i test.query -o ./generated/test.query.idx -d ./generated/test.query.def -u 0 -s 0 -b 1 &&
python ../seqindexer.py -i 100.query -o ./generated/100.query.idx -d ./generated/100.query.def -u 0 -s 0 -b 1 &&
python ../seqindexer.py -i 10seqs.query -o ./generated/10seqs.query.idx -d ./generated/10seqs.query.def -u 1 -b 1 

echo "Done"

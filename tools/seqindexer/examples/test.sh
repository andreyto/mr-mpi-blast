#!/bin/bash 

mkdir -p generated
python ../seqindexer.py -i test.query -o ./generated/test.query.idx -d ./generated/test.query.def -b 1 &&
python ../seqindexer.py -i 100.query -o ./generated/100.query.whole.idx -d ./generated/100.query.whole.def -s 1 -b 1 &&
python ../seqindexer.py -i 100.query -o ./generated/100.query.part.idx -d ./generated/100.query.part.def -s 100 -b 0 

echo "Done"

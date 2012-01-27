#!/bin/bash

echo -e "\n### Make NCBI BLAST database ###"
mkdir -p ./blastdb &&
cd ./blastdb &&
makeblastdb -in test.fa -out test.db -dbtype nucl -logfile blastdbmake.log &&
cd ..
echo -e "\nDone!"

# EOF

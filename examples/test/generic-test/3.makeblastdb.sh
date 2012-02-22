#!/bin/bash

echo -e "\n### Make NCBI BLAST database ###"
mkdir -p ./blastdb &&
cd ./blastdb &&
makeblastdb -in test.fa -out test.db -dbtype nucl -logfile blastdbmake_nucl.log &&
makeblastdb -in test_prot.fa -dbtype prot -out test_prot.db -logfile blastdbmake_prot.log &&

cd ..
echo -e "\nDone!"

# EOF

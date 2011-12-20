#!/bin/bash

echo -e "\n### Make NCBI BLAST database ###"
makeblastdb -in test.fa -out test.db -dbtype nucl -logfile blastdbmake.log && 
echo -e "\nDone!"

# EOF

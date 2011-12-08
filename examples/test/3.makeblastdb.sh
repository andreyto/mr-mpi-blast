#!/bin/bash

echo -e "\n### Make NCBI BLAST database ###"
makeblastdb -in test.fa -dbtype nucl -logfile blastdbmake.log && 
echo -e "\n### Done!"

# EOF

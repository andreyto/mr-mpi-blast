#!/bin/bash 

 

echo -e "\n\n### Convert hit files into CSV file ###"
cd ./blastdb &&
mkdir -p hits &&
#mv *log.txt ./logs &&
mv *.bin ./hits &&
#python ../tools/load_hd5.py ./hits/ hits 1 &&
#python ../tools/load_sql.py ./hits/ hits 1 &&
load_csv.py ./hits/ hits &&

echo -e "\n### Done!"

# EOF

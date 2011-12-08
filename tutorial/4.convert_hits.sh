#!/bin/bash 


echo -e "\n### Convert *.bin to CSV files: ###" &&
echo -e "    hits.csv : only qid" &&
echo -e "    hits_w_defline .csv : qid + original defline" &&
cd ./blastdb &&
mkdir -p hits &&
#mv *log.txt ./logs &&
mv *.bin ./hits &&
load_csv.py -b ./hits/ -o hits &&
load_csv.py -b ./hits/ -o hits_w_defline -d 1 -i ../query/viral_all_query.fa.def  &&

echo -e "\n### Done!"

# EOF

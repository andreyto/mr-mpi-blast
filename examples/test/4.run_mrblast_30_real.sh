#!/bin/bash

export BLASTDB=$(pwd)/blastdb
 
cp ./mrblast_30_real.ini mrblast.ini &&

echo -e "\n### Run mr-mpi-blast ###" &&
mpirun -np 4 mrblast -evalue 10 &&
mkdir -p ./30-real-hits &&
mv *.bin ./30-real-hits &&
echo -e "\n### Convert *.bin to CSV files: ###" &&
echo -e "    30_real_seq.csv : only qid" &&
echo -e "    30_real_seq_w_defline.csv : qid + defline" &&
load_csv.py -b ./30-real-hits -o 30_real_seqs &&
load_csv.py -b ./30-real-hits -o 30_real_seqs_w_defline -d 1 -i ./query/30_real_seqs.fa.def &&

echo -e "    30_real_seq.sqlite : only qid" &&
echo -e "    30_real_seq_w_defline.sqlite : qid + defline" &&
load_sql.py -b ./30-real-hits -o 30_real_seqs &&
load_sql.py -b ./30-real-hits -o 30_real_seqs_w_defline -d 1 -i ./query/30_real_seqs.fa.def &&

echo -e "    30_real_seq.hd5 : only qid" &&
load_hd5.py -b ./30-real-hits -o 30_real_seqs 


echo -e "\nDone!"

# EOF

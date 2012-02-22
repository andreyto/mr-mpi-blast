#!/bin/bash

export BLASTDB=$(pwd)/blastdb
#export MRBLASTDIR=/home/ssul/work/mrblast-build

cp ./mrblast_30_real_w_numiter.ini mrblast.ini &&

echo -e "\n### Run mr-mpi-blast ###" &&
mpirun -np 4 mrblast -evalue 10 -task megablast &&
 
mkdir -p ./30-real-hits-numiter &&
mv *.bin ./30-real-hits-numiter &&
mkdir -p ./30-real-hits-numiter-output &&
echo -e "\n### Convert *.bin to CSV files: ###" &&
echo -e "    30_real_seq_numiter.csv : only qid" &&
#echo -e "    30_real_seq_w_defline_numiter.csv : qid + defline" &&
load_csv_classifier.py -b ./30-real-hits-numiter -o ./30-real-hits-numiter-output/30_real_seqs_numiter &&
#load_csv_classifier.py -b ./30-real-hits-numiter -o ./30-real-hits-numiter-output/30_real_seqs_w_defline_numiter -d 1 -i ./query/30_real_seqs.fa.def &&

echo -e "    30_real_seq.sqlite : only qid" &&
#echo -e "    30_real_seq_w_defline.sqlite : qid + defline" &&
load_sql_classifier.py -b ./30-real-hits-numiter -o ./30-real-hits-numiter-output/30_real_seqs_numiter &&
#load_sql_classifier.py -b ./30-real-hits-numiter -o ./30-real-hits-numiter-output/30_real_seqs_w_defline_numiter -d 1 -i ./query/30_real_seqs.fa.def &&

echo -e "    30_real_seq.hd5 : only qid" &&
#echo -e "    30_real_seq_w_defline.hd5 : qid + defline" &&
load_hd5_classifier.py -b ./30-real-hits-numiter -o ./30-real-hits-numiter-output/30_real_seqs_numiter &&
#load_hd5_classifier.py -b ./30-real-hits-numiter -o ./30-real-hits-numiter-output/30_real_seqs_w_defline -d 1 -i ./query/30_real_seqs.fa.def 

echo -e "\nDone!"

# EOF

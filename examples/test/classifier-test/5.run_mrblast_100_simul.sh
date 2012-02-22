#!/bin/bash

export BLASTDB=$(pwd)/blastdb
#export MRBLASTDIR=/home/ssul/work/mrblast-build

cp ./mrblast_100_simul.ini mrblast.ini &&

echo -e "\n### Run mr-mpi-blast with 100 simulated queries ###" &&
mpirun -np 4 mrblast -evalue 10 -task megablast  &&
mkdir -p ./100-simul-hits &&
mv *.bin ./100-simul-hits &&
mkdir -p ./100-simul-hits-output &&
echo -e "\n### Convert *.bin to CSV files: ###" &&
echo -e "    100_simul_seq.csv : only qid" &&
echo -e "    100_simul_seq_w_defline.csv : qid + defline" &&
load_csv_classifier.py -b ./100-simul-hits -o ./100-simul-hits-output/100_simul_seqs &&
load_csv_classifier.py -b ./100-simul-hits -o ./100-simul-hits-output/100_simul_seqs_w_defline -d 1 -i ./query/100_simul_seqs.fa.def &&

echo -e "    100_simul_seq.sqlite : only qid" &&
echo -e "    100_simul_seq.sqlite : qid + defline" &&
load_sql_classifier.py -b ./100-simul-hits -o ./100-simul-hits-output/100_simul_seqs &&
load_sql_classifier.py -b ./100-simul-hits -o ./100-simul-hits-output/100_simul_seqs_w_defline -d 1 -i ./query/100_simul_seqs.fa.def &&

echo -e "    100_simul_seq.hd5 : only qid" &&
echo -e "    100_simul_seq.hd5 : qid + defline" &&
load_hd5_classifier.py -b ./100-simul-hits -o ./100-simul-hits-output/100_simul_seqs &&
load_hd5_classifier.py -b ./100-simul-hits -o ./100-simul-hits-output/100_simul_seqs_w_defline -d 1 -i ./query/100_simul_seqs.fa.def &&

echo -e "\nDone!"


# EOF

#!/bin/bash

export BLASTDB=$(pwd)/blastdb

cp ./mrblast_100_simul.ini mrblast.ini &&

echo -e "\n### Run mr-mpi-blast with 100 simulated queries ###" &&
mpirun -np 4 mrblast -evalue 10 &&
mkdir -p ./100-simul-hits &&
mv *.bin ./100-simul-hits &&
echo -e "\n### Convert *.bin to CSV files: ###" &&
echo -e "    100_simul_seq.csv : only qid" &&
echo -e "    100_simul_seq_w_defline.csv : qid + defline" &&
load_csv.py -b ./100-simul-hits -o 100_simul_seqs &&
load_csv.py -b ./100-simul-hits -o 100_simul_seqs_w_defline -d 1 -i ./query/100_simul_seqs.fa.def &&

load_sql.py -b ./100-simul-hits -o 100_simul_seqs &&
load_sql.py -b ./100-simul-hits -o 100_simul_seqs_w_defline -d 1 -i ./query/100_simul_seqs.fa.def &&

load_hd5.py -b ./100-simul-hits -o 100_simul_seqs &&

echo -e "\nDone!"


# EOF

#!/bin/bash
 
cp ./mrblast_30_real_w_numiter.ini mrblast.ini &&

echo -e "\n### Run mr-mpi-blast ###" &&
mpirun -np 4 mrblast -db - -evalue 10 &&
mkdir -p ./30-real-hits-numiter &&
mv *.bin ./30-real-hits-numiter &&
echo -e "\n### Convert *.bin to CSV files: ###" &&
echo -e "    30_real_seq_numiter.csv : only qid" &&
echo -e "    30_real_seq_w_defline_numiter.csv : qid + original defline" &&
load_csv.py -b ./30-real-hits-numiter -o 30_real_seqs_numiter &&
load_csv.py -b ./30-real-hits-numiter -o 30_real_seqs_w_defline_numiter -d 1 -i 30_real_seqs.fa.def &&

load_sql.py -b ./30-real-hits-numiter -o 30_real_seqs_numiter &&
load_sql.py -b ./30-real-hits-numiter -o 30_real_seqs_w_defline_numiter -d 1 -i 30_real_seqs.fa.def &&

load_hd5.py -b ./30-real-hits-numiter -o 30_real_seqs_numiter &&


echo -e "\n### Done!"

# EOF

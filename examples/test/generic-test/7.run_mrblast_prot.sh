#!/bin/bash

export BLASTDB=$(pwd)/blastdb
export MRBLASTDIR=/home/ssul/work/mrblast-build

echo -e "\n### Indexing protein query file ###" &&
seqindexer.py -i ./query/prot_seqs.fa -o ./query/prot_seqs.fa.idx -d ./query/prot_seqs.fa.def -s 1 -b 1 &&

cp ./mrblast_prot.ini mrblast.ini &&

echo -e "\n### Run mr-mpi-blast ###" &&

mpirun -np 4 $MRBLASTDIR/mrblast -evalue 10 -task blastp &&

mkdir -p ./prot-hits &&
mv *.bin ./prot-hits &&
mkdir -p ./prot-hits-output &&

echo -e "\n### Convert *.bin to CSV files: ###" &&
echo -e "    prot_seq.csv : only qid" &&
echo -e "    prot_seq_w_defline.csv : qid + defline" &&
load_csv.py -b ./prot-hits -o ./prot-hits-output/prot_seqs &&
load_csv.py -b ./prot-hits -o ./prot-hits-output/prot_seqs_w_defline -d 1 -i ./query/prot_seqs.fa.def &&

echo -e "    prot_seq.sqlite : only qid" &&
echo -e "    prot_seq_w_defline.sqlite : qid + defline" &&
load_sql.py -b ./prot-hits -o ./prot-hits-output/prot_seqs &&
load_sql.py -b ./prot-hits -o ./prot-hits-output/prot_seqs_w_defline -d 1 -i ./query/prot_seqs.fa.def &&

echo -e "    prot_seq.hd5 : only qid" &&
echo -e "    prot_seq_w_defline.hd5 : qid + defline" &&
load_hd5.py -b ./prot-hits -o ./prot-hits-output/prot_seqs 
load_hd5.py -b ./prot-hits -o ./prot-hits-output/prot_seqs_w_defline -d 1 -i ./query/prot_seqs.fa.def 


echo -e "\nDone!"

# EOF

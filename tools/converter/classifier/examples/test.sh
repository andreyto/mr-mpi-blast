#!/bin/bash

mkdir -p ./converted &&
echo -e "\n### Convert *.bin to CSV files: ###" &&
echo -e "    100_simul_seq.csv : only qid" &&
echo -e "    100_simul_seq_w_defline.csv : qid + defline" &&
load_csv_classifier.py -b ./hits -o ./converted/100_simul_seqs &&
load_csv_classifier.py -b ./hits -o ./converted/100_simul_seqs_w_defline -d 1 -i ./100_simul_seqs.fa.def &&

echo -e "    100_simul_seq.sqlite : only qid" &&
echo -e "    100_simul_seq.sqlite : qid + defline" &&
load_sql_classifier.py -b ./hits -o ./converted/100_simul_seqs &&
load_sql_classifier.py -b ./hits -o ./converted/100_simul_seqs_w_defline -d 1 -i ./100_simul_seqs.fa.def &&

echo -e "    100_simul_seq.hd5 : only qid" &&
echo -e "    100_simul_seq.hd5 : qid + defline" &&
load_hd5_classifier.py -b ./hits -o ./converted/100_simul_seqs &&
load_hd5_classifier.py -b ./hits -o ./converted/100_simul_seqs_w_defline -d 1 -i ./100_simul_seqs.fa.def &&

echo -e "\nDone!"


# EOF

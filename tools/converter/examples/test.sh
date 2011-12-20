#!/bin/bash 

mkdir -p converted
echo -e "\n### Convert *.bin to CSV files: ###" &&
echo -e "    30_real_seq.csv : only qid" &&
echo -e "    30_real_seq_w_defline.csv : qid + defline" &&
../load_csv.py -b ./hits -o ./converted/30_real_seqs &&
../load_csv.py -b ./hits -o ./converted/30_real_seqs_w_defline -d 1 -i 30_real_seqs.fa.def &&

echo -e "    30_real_seq.sqlite : only qid" &&
echo -e "    30_real_seq_w_defline.sqlite : qid + defline" &&
../load_sql.py -b ./hits -o ./converted/30_real_seqs &&
../load_sql.py -b ./hits -o ./converted/30_real_seqs_w_defline -d 1 -i 30_real_seqs.fa.def &&

echo -e "    30_real_seq.hd5 : only qid" &&
../load_hd5.py -b ./hits -o ./converted/30_real_seqs 

echo -e "Done"

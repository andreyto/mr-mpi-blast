#!/bin/bash

cp ./mrblast_30_real.ini mrblast.ini &&

echo -e "\n### Run mr-mpi-blast ###" &&
mpirun -np 4 mrblast -db - -evalue 10 &&
echo -e "\n### Convert *.bin to 30_real_seq.csv ###" &&
python ../tools/load-hd5-tables/load_csv.py . 30_real_seqs &&
echo -e "\n### Done!"


# EOF

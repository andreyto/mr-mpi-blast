#!/bin/bash

cp ./mrblast_100_simul.ini mrblast.ini &&

echo -e "\n### Run mr-mpi-blast with 100 simulated queries ###" &&
mpirun -np 4 mrblast -db - -evalue 10 &&
echo -e "\n### Convert *.bin to 100_simul_seq.csv ###" &&
python ../tools/load-hd5-tables/load_csv.py . 100_simul_seqs &&
echo -e "\n### Done!"


# EOF

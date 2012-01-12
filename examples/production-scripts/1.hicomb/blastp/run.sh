#!/bin/sh

mpirun -np 4 ../../mrblast -i 1seq_prot.query -d 1seq_prot.query.idx -n 1 -b 1 -s new_blastp_option.out -p 1


#!/bin/bash


if [ -z "$SGE_ACCOUNT" ]; then
    echo "Please set SGE_ACCOUNT environment variable."
    exit
fi


cd ./blastdb &&
cp ../dblist.txt . &&
cp ../mrblast.ini . &&
cp ../mrblast_job.sh . &&

echo -e "\n\n### qsub mr-mpi-blast job ###"
qsub mrblast_job.sh &&



# EOF


#!/bin/bash

if [ -z "$MRMPIBLAST_PREFIX" ]; then
    echo "The environment variable, MRMPIBLAST_PREFIX is not set."
    exit
fi

if [ -z "$SGE_ACCOUNT" ]; then
    echo "Please set SGE_ACCOUNT environment variable."
    exit
fi


cd $MRMPIBLAST_PREFIX/tutorial/blastdb &&
cp $MRMPIBLAST_PREFIX/tutorial/tools/dblist.txt . &&
cp $MRMPIBLAST_PREFIX/tutorial/tools/mrblast.ini . &&
cp $MRMPIBLAST_PREFIX/tutorial/tools/mrblast_job.sh . &&

echo -e "\n\n### qsub mr-mpi-blast job ###"
qsub mrblast_job.sh &&



# EOF


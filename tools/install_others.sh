#!/bin/bash 

if [ -z "$MRMPIBLAST_PREFIX" ]; then
    echo "The environment variable, MRMPIBLAST_PREFIX is not set."
    exit
fi

 
echo -e "\n\n### Copy other files to MRMPIBLAST_PREFIX"
cp -R ./src/app/mr-mpi-blast/tutorial $MRMPIBLAST_PREFIX &&

#cp -R ./src/app/mr-mpi-blast/tools/seqindexer/seqindexer.py $MRMPIBLAST_PREFIX/tutorial/tools &&
#cp -R ./src/app/mr-mpi-blast/tools/splitter/splitter.py $MRMPIBLAST_PREFIX/tutorial/tools &&
#cp -R ./src/app/mr-mpi-blast/tools/load-hd5-tables/*.py $MRMPIBLAST_PREFIX/tutorial/tools &&

cp -R ./src/app/mr-mpi-blast/tools/seqindexer/seqindexer.py $MRMPIBLAST_PREFIX/bin &&
cp -R ./src/app/mr-mpi-blast/tools/splitter/splitter.py $MRMPIBLAST_PREFIX/bin &&
cp -R ./src/app/mr-mpi-blast/tools/load-hd5-tables/*.py $MRMPIBLAST_PREFIX/bin &&

cp -R ./src/app/mr-mpi-blast/examples $MRMPIBLAST_PREFIX &&
cp -R ./src/app/mr-mpi-blast/docs $MRMPIBLAST_PREFIX &&


echo -e "\n### Done!"
# EOF

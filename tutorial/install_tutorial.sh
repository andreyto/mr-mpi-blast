#!/bin/bash 

if [ $MRMPIBLAST_PREFIX == "" ]; then
    echo "The environment variable, MRMPIBLAST_PREFIX is not set."
    exit
fi

if [ $WORKING_DIR == "" ]; then
    echo "The environment variable, WORKING_DIR is not set."
    exit
fi

echo -e "\n\n### Copy tutorial files to MRMPIBLAST_PREFIX/tutorial"
cp -R $WORKING_DIR/ncbi_cxx--7_0_0/src/app/mr-mpi-blast/tutorial $MRMPIBLAST_PREFIX/ &&
cp -R $WORKING_DIR/ncbi_cxx--7_0_0/src/app/mr-mpi-blast/tools/seqindexer/seqindexer.py $MRMPIBLAST_PREFIX/tutorial/tools &&
cp -R $WORKING_DIR/ncbi_cxx--7_0_0/src/app/mr-mpi-blast/tools/splitter/splitter.py $MRMPIBLAST_PREFIX/tutorial/tools &&
cp -R $WORKING_DIR/ncbi_cxx--7_0_0/src/app/mr-mpi-blast/tools/load-hd5-tables/*.py $MRMPIBLAST_PREFIX/tutorial/tools

echo -e "\n### Done!"
# EOF

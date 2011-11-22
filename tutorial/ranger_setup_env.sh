#!/bin/bash 

#source ~/.environ
#source ./.profile_user


if [ $MRMPIBLAST_PREFIX == "" ]; then
    echo "The environment variable, MRMPIBLAST_PREFIX is not set."
    exit
fi

if [ $WORKING_DIR == "" ]; then
    echo "The environment variable, WORKING_DIR is not set."
    exit
fi

echo -e "\n\n### Load modules ###"
module load python/2.7.1
module swap pgi gcc
module load cmake
module unload openmpi
module load mvapich
module load boost
module list


echo -e "\n\n### Copy make related files ###"
cp -R $WORKING_DIR/ncbi_cxx--7_0_0/src/app/mr-mpi-blast/tutorial/makefiles/build_mrmpi.sh $WORKING_DIR/ncbi_cxx--7_0_0/src/app/mr-mpi-blast/ &&

cp -R $WORKING_DIR/ncbi_cxx--7_0_0/src/app/mr-mpi-blast/tutorial/makefiles/Makefile.in $WORKING_DIR/ncbi_cxx--7_0_0/src/app/mr-mpi-blast/ &&

cp -R $WORKING_DIR/ncbi_cxx--7_0_0/src/app/mr-mpi-blast/tutorial/makefiles/Makefile.in.root $WORKING_DIR/ncbi_cxx--7_0_0/src/app/Makefile.in &&

cp -R $WORKING_DIR/ncbi_cxx--7_0_0/src/app/mr-mpi-blast/tutorial/makefiles/Makefile.mpicxx $WORKING_DIR/ncbi_cxx--7_0_0/src/app/mr-mpi-blast/mrmpi/MAKE/ &&

cp -R $WORKING_DIR/ncbi_cxx--7_0_0/src/app/mr-mpi-blast/tutorial/makefiles/Makefile.mrblast.app $WORKING_DIR/ncbi_cxx--7_0_0/src/app/mr-mpi-blast/

echo -e "\n### Done!"
# EOF

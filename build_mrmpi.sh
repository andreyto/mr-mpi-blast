#!/bin/bash
set -ex #abort if any command fails and echo all commands

################################################################################
MRMPI_HOME=$WORKING_DIR/ncbi_cxx--7_0_0/src/app/mr-mpi-blast/mrmpi
################################################################################

topdir=$(pwd)
cd $MRMPI_HOME &&
make clean-all &&
make mpicxx &&
#make mpicc &&
cd ${topdir}

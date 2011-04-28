#!/bin/bash
set -ex #abort if any command fails and echo all commands

################################################################################
MRMPI_HOME=/home/ssul/work/distros2/ncbi_cxx/ncbi_cxx--Jun_15_2010/src/app/mrblast/mrmpi
#MRMPI_HOME=/work/01471/ssul/work/distros3/ncbi_cxx/ncbi_cxx--Jun_15_2010/src/app/mr-mpi-blast/mrmpi
################################################################################

topdir=$(pwd)
cd $MRMPI_HOME &&
make clean-all &&
### Ranger
#make mpicxx &&
### Local
make mpicc &&
cd ${topdir}

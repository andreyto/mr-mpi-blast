#!/bin/bash
set -ex #abort if any command fails and echo all commands
topdir=$(pwd)

### Ranger
#cd ../../../../src/app/mr-mpi-blast/mrmpi &&
#make clean-all &&
#make mpicxx &&

### Local
cd ../../../../src/app/mrblast/mrmpi &&
make clean-all &&
make mpicc &&

cd ${topdir}

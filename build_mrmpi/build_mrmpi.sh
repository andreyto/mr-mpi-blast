#!/bin/bash
set -ex #abort if any command fails and echo all commands
topdir=$(pwd)
cd ../../../../src/app/mrblast/mrmpi &&
make clean-all &&
### Ranger
#make mpicxx &&
### Local
make mpicc &&
cd ${topdir}

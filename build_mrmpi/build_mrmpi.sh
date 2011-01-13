#!/bin/bash
set -ex #abort if any command fails and echo all commands
topdir=$(pwd)
cd ../../../../src/app/mrblast/mrmpi &&
rm -rf *.d *.o &&
#make -f Makefile.mvapich &&
make -f Makefile.mpicc &&
cd ${topdir}

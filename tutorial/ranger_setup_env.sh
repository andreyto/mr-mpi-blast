#!/bin/bash 
  
echo -e "\n\n### Load modules"
module load python
module swap pgi gcc
module load cmake
module unload openmpi
module load mvapich
module load boost
module load numpy
module list
 
echo -e "\nDone!"
# EOF

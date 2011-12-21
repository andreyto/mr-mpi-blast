#!/bin/bash 
  
echo -e "\n\n### Load modules"
module swap pgi gcc
module unload openmpi
module load mvapich
module load boost
module load python
module list
 
echo -e "\nDone!"
# EOF

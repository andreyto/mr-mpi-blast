#!/bin/bash 

#source ~/.environ
#source ./.profile_user

 
echo -e "\n\n### Load modules"
module load python
module swap pgi gcc
module load cmake
module unload openmpi
module load mvapich
module load boost
module list
 
echo -e "\n### Done!"
# EOF

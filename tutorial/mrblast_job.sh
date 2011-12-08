#!/bin/bash
#$ -A $SGE_ACCOUNT
#$ -V # Inherit the submission environment
#$ -cwd # Start job in submission directory
#$ -N mrblast # Job Name
#$ -j y # Combine stderr and stdout
#$ -o $JOB_NAME.o$JOB_ID # Name of the output file (eg. myMPI.oJobID)
#$ -pe 16way 64  # Requests 16 tasks/node, 32 cores total
#$ -q development # Queue name "normal"
#$ -l h_rt=01:00:00 # Run time (hh:mm:ss) 
##$ -M # Use email notification address
##$ -m be # Email at Begin and End of job

ibrun mrblast -db - -evalue 1e-4

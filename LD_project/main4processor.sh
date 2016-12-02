#!/bin/tcsh

# Instructions for new Lightningsmp Cluster users:
#  To use this script:
#   1) Save this script as a file named myscript on Lightningsmp.its
#   2) On lightningsmp, Issue                   
#       qsub myscript    to submit the job 
#        Use qstat -a to see job status, 
#         Use qdel jobname to delete one of your jobs
#         jobnames are of the form 1234.hpc4smp 

###########################################
# Output goes to file BATCH_OUTPUT.
# Error output goes to file BATCH_ERRORS.
# If you want the output to go to another file, change BATCH_OUTPUT 
# or BATCH_ERRORS in the following lines to the full path of that file. 
# HAO: In lightning, 16 cpus is one node. processor means cpus

#PBS  -o BATCH_OUTPUT 
#PBS  -e BATCH_ERRORS 

#PBS -lvmem=64GB,pmem=8Gb,mem=64Gb,nodes=1:ppn=16:ib,cput=80:00:00,walltime=10:00:00

# Change to directory from which qsub command was issued 
   cd $PBS_O_WORKDIR

#argv recombination/mutation/generation
mpirun -np 8 ./main4processor 0.002 0.000000001 2000 

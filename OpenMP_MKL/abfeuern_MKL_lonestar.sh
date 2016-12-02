
qsub -v MP=1 run_lightning_1node4mpi.sh
qsub -v MP=48 run_lightning_1node4mpi.sh
qsub -v MP=48 run_lonestar_2node4mpi.sh

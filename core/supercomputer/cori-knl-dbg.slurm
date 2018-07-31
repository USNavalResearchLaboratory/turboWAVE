#!/bin/bash
#SBATCH -C knl,quad,cache
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --job-name=wake
#SBATCH --output=tw.out
#SBATCH --error=tw.err
#SBATCH --account=m1689

#Enable core dumps (not needed for ATP)
#ulimit -c unlimited

#Enable ATP core dump analysis
export ATP_ENABLED=1

#OpenMP settings:
export OMP_NUM_THREADS=17
export OMP_PLACES=threads
export OMP_PROC_BIND=FALSE
#export KMP_AFFINITY=verbose

#Load valgrind for debugging purposes
#module load valgrind

#run the application with valgrind:
#srun --ntasks 16 --cpus-per-task 17 --cpu_bind=no valgrind ./tw3d
#run the application without valgrind (still produces ATP report):
srun --ntasks 16 --cpus-per-task 17 --cpu_bind=no ./tw3d

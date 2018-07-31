#!/bin/bash
#SBATCH -C knl,quad,cache
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --time=00:05:00
#SBATCH --job-name=dirac
#SBATCH --output=tw.out
#SBATCH --error=tw.err
#SBATCH --account=m1689

#Enable core dumps
#ulimit -c unlimited

#Enable ATP core dump analysis
#export ATP_ENABLED=1

#OpenMP settings:
export OMP_NUM_THREADS=8
export OMP_PLACES=cores
export OMP_PROC_BIND=true

#Load valgrind for debugging purposes
#module load valgrind

#Run the application.
#Ask for all logical cores even if not all are used.
srun --ntasks 17 --cpus-per-task 16 --cpu-bind=cores ./tw3d
#!/bin/bash -l
#SBATCH -C haswell
#SBATCH --partition=regular
#SBATCH --nodes=1
#SBATCH --time=00:05:00
#SBATCH --job-name=compression
#SBATCH --output=tw.out
#SBATCH --error=tw.err
#SBATCH --account=m1689

export OMP_NUM_THREADS=4
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

srun --ntasks 16 --cpus-per-task 4 ./tw3d

#!/bin/csh

#PBS -A NRLDC08284338
#PBS -l walltime=12:00:00
#PBS -l select=911:ncpus=36:mpiprocs=36
#PBS -q standard
#PBS -N turbowave
#PBS -o tw.out
#PBS -e tw.err

cd $WORKDIR
mpiexec_mpt -np 32768 ./tw3d
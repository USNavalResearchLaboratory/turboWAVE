#!/bin/csh

#PBS -A NRLDC08284338
#PBS -l walltime=12:00:00
#PBS -l select=32:ncpus=28:mpiprocs=1:ngpus=1
#PBS -q standard
#PBS -N turbowave
#PBS -o tw.out
#PBS -e tw.err

cd $WORKDIR
mpiexec_mpt -np 32 ./tw3d
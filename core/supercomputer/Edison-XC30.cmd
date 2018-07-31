#!/bin/csh
#PBS -q regular
#PBS -l mppwidth=32
#PBS -l walltime=00:30:00
#PBS -A m1689
#PBS -e tw.err
#PBS -o tw.out
#PBS -N turbowave
setenv MPICH_UNEX_BUFFER_SIZE 300M
setenv MPICH_PTL_UNEX_EVENTS 40960
cd $SCRATCH
# j option is for SMT : use 2 virtual cores per physical core
# physical cores appears in mppwidth line above
# virtual cores appears in -n option below
# uncomment the following to generate core files
# limit coredumpsize unlimited
aprun -j 2 -n 64 ./tw3d

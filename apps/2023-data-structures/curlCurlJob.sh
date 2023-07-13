#!/bin/bash

#PBS -N gmgConvTest
#PBS -l select=2:node_type=rome:mpiprocs=128
#PBS -l walltime=04:00:00

module load petsc

cd $PBS_O_WORKDIR
mpirun -n 256 ./curlCurlConvergence > out 2>&1

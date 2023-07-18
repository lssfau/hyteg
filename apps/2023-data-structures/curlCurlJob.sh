#!/bin/bash

#PBS -N gmgConvTest
#PBS -l select=5:node_type=rome:mpiprocs=128
#PBS -l walltime=08:00:00

module load petsc

cd $PBS_O_WORKDIR
mpirun -n 640 ./curlCurlConvergence > out 2>&1

#!/bin/bash
#PBS -N tokamak-4nodes
#PBS -l select=4:node_type=rome:mpiprocs=128
#PBS -l walltime=00:05:00
#PBS -m abe
#PBS -M dominik.thoennes@fau.de

export MPI_LAUNCH_TIMEOUT=1000
export MPI_IB_CONGESTED=1

module load gcc
module load mpt
module load cmake
module load petsc
module list


mpirun -np 512 --map-by core --bind-to core ./Tokamak Tokamak.prm -Parameters.maxLevel=7 -Parameters.refineCoarseMesh=1

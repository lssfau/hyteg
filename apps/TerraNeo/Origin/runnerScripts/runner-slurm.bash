#!/bin/bash
##
## This is a sample script for running
## programs on slurm systems
## This specific one was used on tethys-3g (geophysik)
##
## Absolute paths are better for cluster jobs

#SBATCH --output=./slurm-out/slurm-%j.out
#SBATCH --chdir=./
#SBATCH --job-name="tn-app"
#SBATCH --ntasks=240
#SBATCH --time=71:59:59

## Load all modules that were used to compile the executable, if required
# module -s load intel/2024
# module -s load mpi.intel/2024

mpirun.openmpi TerraNeo parameters.prm

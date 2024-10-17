#!/bin/bash
#SBATCH --output=/SCRATCH/pponkumar/sphbench/slurm-out/slurm-%j.out
#SBATCH --chdir=/SCRATCH/pponkumar/
#SBATCH --job-name="sph64re5p1"
#SBATCH --ntasks=128
#SBATCH --time=71:59:59

# module -s load intel/2024
# module -s load mpi.intel/2024

mpirun.openmpi /home/pponkumar/hyteg/hyteg-build/apps/2024-convbench/SphericalShellBenchRotation /home/pponkumar/hyteg/hyteg/apps/2024-convbench/SphericalShellBenchRotation.prm


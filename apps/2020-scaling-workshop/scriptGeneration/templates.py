param_template = """Parameters
{{
  discretization p2p1;

  minLevel 0;
  maxLevel {maxLevel};

  numEdgesPerSide {numEPS};

  scenario 1;
  domainInfoOnly true;

  vtk true;
  dbFile Scaling_Workshop_01_Cube.db;

  preSmooth 2;
  postSmooth 2;
  incSmooth 2;
  fmgInnerIterations 1;
  numCycles 1;
  absoluteResidualTolerance 1e-12;

  numGSVelocity 3;
  symmGSVelocity false;

  estimateOmega false;
  omega 0.2;
  omegaEstimationLevel 5;
  omegaEstimationIterations 100;

  coarseGridSolverType 1;
  maxIterations 10000;
  coarseGridAbsoluteResidualTolerance 1e-12;
}}"""

job_template = """#!/bin/bash
#PBS -N {job_name}
#PBS -l select={num_nodes}:node_type=rome:mpiprocs={num_cores}
#PBS -l walltime={walltime}
#PBS -m abe
#PBS -M dominik.thoennes@fau.de

export MPI_LAUNCH_TIMEOUT=1000
export MPI_IB_CONGESTED=1

module load gcc
module load mpt
module load cmake
module load petsc
module list

cd $PBS_O_WORKDIR
cd ..
pwd
ls -lha

mpirun -np {total_num_procs} omplace -c 0-{num_cores}:st=4 ./Scaling_Workshop_01_Cube hawk/{paramfile_name}"""
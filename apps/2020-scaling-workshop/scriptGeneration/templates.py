
def parameter_file_01_cube(scenario: int, max_level: int, num_edges_per_side: int, db_file: str, **kwargs):
    return f"""Parameters
{{
  discretization p2p1;

  minLevel 0;
  maxLevel {max_level};

  numEdgesPerSide {num_edges_per_side};

  scenario {scenario};
  domainInfoOnly false;

  vtk false;
  dbFile {db_file};

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
}}
"""


def job_file_hawk(job_name: str, num_nodes: int, num_cores: int, walltime: str, total_num_procs: int, paramfile_name: str, **kwargs):
    return f"""#!/bin/bash
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

mpirun -np {total_num_procs} omplace -c 0-{num_cores}:st=4 ./Scaling_Workshop_01_Cube hawk/{paramfile_name}
"""

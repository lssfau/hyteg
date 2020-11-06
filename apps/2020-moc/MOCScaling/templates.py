

def create_parameter_file(level: int, db_file_name: str, num_edges_per_side: int, **kwargs):
    return f"""Parameters
{{
    level {level};
    numTimeSteps 10;
    setTimeStepSizeManually true;
    manualTimeStepSize {1.5 / (2**level * num_edges_per_side)};
    threeDim true;

    enableGaussianCone false;
    enableLinearCone false;
    enableCylinder false;

    resetParticles false;
    resetParticlesInterval 1;
    adjustedAdvection false;
    vtk false;
    printInterval 1;
    vtkInterval 1;
    dbFile {db_file_name};

    numEdgesPerSide {num_edges_per_side};
    }}
"""



def job_file_hawk(job_name: str, binary_name: str, num_nodes: int, num_omp_threads_per_mpi_proc: int, walltime: str, total_num_procs: int, paramfile_name: str,
                  path: str, **kwargs):
    return f"""#!/bin/bash
#PBS -N {job_name}
#PBS -l select={num_nodes}:node_type=rome:mpiprocs={num_omp_threads_per_mpi_proc}
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

mpirun -np {total_num_procs} omplace -c 0-128:st={int(128 / num_omp_threads_per_mpi_proc)} {path}{binary_name} {path}hawk/{paramfile_name}
"""


def job_file_supermuc(job_name: str, binary_name: str, num_nodes: int, num_mpi_procs_per_node: int, walltime: str, paramfile_name: str, **kwargs):

    def partition(num_nodes):
        if num_nodes <= 16:
            return "micro"
        elif num_nodes <= 768:
            return "general"
        elif num_nodes <= 3072:
            return "large"

    constraint = ""
    if num_nodes <= 792:
        constraint = "#SBATCH --constraint=[i01|i02|i03|i04|i05|i06|i07|i08]"

    return f"""#!/bin/bash
# Job Name and Files (also --job-name)
#SBATCH -J {job_name}
#Output and error (also --output, --error):
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#Initial working directory (also --chdir):
#SBATCH -D ./
#Notification and type
#SBATCH --mail-type=END
#SBATCH --mail-user=nils.kohl@fau.de
# Wall clock limit:
#SBATCH --time={walltime}
#SBATCH --no-requeue
#Setup of execution environment
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pr86ma
 
#SBATCH --ear=off
#SBATCH --partition={partition(num_nodes)}
#Number of nodes and MPI tasks per node:
#SBATCH --nodes={num_nodes}
#SBATCH --ntasks-per-node={num_mpi_procs_per_node}
{constraint}

source load_modules_supermuc.sh

cd ..
pwd
ls -lha

module list

mpiexec -n $SLURM_NTASKS ./{binary_name} MOCScaling/{paramfile_name}
"""

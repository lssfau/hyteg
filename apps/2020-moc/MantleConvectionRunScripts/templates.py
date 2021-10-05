def create_parameter_file(max_level: int, ra: float, output_directory: str, base_name: str, max_num_time_steps: int,
                          uzawa_omega: float, cfl: float, uzawa_pre: int, uzawa_post: int, uzawa_inner: int,
                          stokes_rel_tol: float, stokes_abs_tol: float, ntan: int, nrad: int, vtk_interval: int,
                          vtk_vertex_dofs: bool, **kwargs):
    return f"""Parameters
{{
    level {max_level};
    minLevel 0;

    rMin 0.55;
    rMax 1.0;

    // nTan 12;
    // nRad 2;
    // threeDim false;

    nTan {ntan};
    nRad {nrad};
    threeDim true;

    // PETSC_MUMPS         = 0,
    // PETSC_MINRES_JACOBI = 1,
    // PETSC_MINRES_BOOMER = 2,
    // HYTEG_MINRES        = 3,
    // HYTEG_MINRES_GMG    = 4,
    // HYTEG_UZAWA_V       = 5
    // HYTEG_UZAWA_FMG     = 6

    stokesSolverType 5;

    stokesMaxNumIterations 20;
    stokesAbsoluteResidualUTolerance {stokes_abs_tol};
    stokesRelativeResidualUTolerance {stokes_rel_tol};
    uzawaOmega {uzawa_omega};
    uzawaPreSmooth {uzawa_pre};
    uzawaPostSmooth {uzawa_post};
    uzawaInnerIterations {uzawa_inner};

    // PETSC_MINRES = 0,
    // HYTEG_CG     = 1
    // HYTEG_GMG    = 2
    diffusionSolverType 1;

    diffusionMaxNumIterations 10000;
    diffusionAbsoluteResidualUTolerance 1e-12;

    maxNumTimeSteps {max_num_time_steps};
    simulationTime 10;
    cflMax {cfl};
    fixedTimeStep false;
    dtConstant 1e-4;
    rayleighNumber {ra};
    vtk true;
    vtkOutputVelocity false;
    vtkOutputInterval {vtk_interval};
    vtkOutputVertexDoFs {vtk_vertex_dofs};

    outputDirectory {output_directory};
    outputBaseName {base_name};

    verbose true;
    }}
"""


def job_file_hawk(job_name: str, binary_name: str, num_nodes: int, num_mpi_procs_per_node: int, walltime: str,
                  paramfile_name: str, **kwargs):
    petsc_detail_string = "-ksp_view -ksp_monitor -log_view -mat_mumps_icntl_4 2"

    return f"""#!/bin/bash
#PBS -N {job_name}
#PBS -l select={num_nodes}:node_type=rome:mpiprocs={num_mpi_procs_per_node}
#PBS -l walltime={walltime}
#PBS -m abe
#PBS -M nils.kohl@fau.de

export MPI_LAUNCH_TIMEOUT=1000
export MPI_IB_CONGESTED=1

source load_modules_hawk.sh

cd $PBS_O_WORKDIR
cd ..
pwd
ls -lha

mpirun -np {num_nodes*num_mpi_procs_per_node} omplace -c 0-128:st={int(128 / num_mpi_procs_per_node)} ./{binary_name} MantleConvectionRunScripts/{paramfile_name} {petsc_detail_string}
"""


def job_file_supermuc(job_name: str, binary_name: str, num_nodes: int, num_mpi_procs_per_node: int, walltime: str,
                      paramfile_name: str, **kwargs):
    petsc_detail_string = "-ksp_view -ksp_monitor -log_view -mat_mumps_icntl_4 2"

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

mpiexec -n $SLURM_NTASKS ./{binary_name} MantleConvectionRunScripts/{paramfile_name} {petsc_detail_string}
"""

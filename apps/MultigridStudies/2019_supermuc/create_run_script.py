import time
import datetime


def supermuc_scaling_prm_file_string(discretization="P2", mesh_spherical_shell=False, num_faces_per_side=1, ntan=2, nrad=2,
                                     num_cycles=1, fmg_r=1, omega=0.2, pre=3, post=3, num_gs_velocity=1, max_level=3, timing_file="timing.json", db_file="database.db"):

    base_config = """
Parameters
{{
    equation stokes;
    dim 3;
    numFacesPerSide {num_faces_per_side};
    discretization {discretization};

    // number of tets = 60 * (ntan-1) * (ntan-1) * (nrad-1)
    meshSphericalShell {mesh_spherical_shell};
    shellNTan {ntan};
    shellNRad {nrad};
    shellRMin 0.55;
    shellRMax 1.0;

    meshLayout CRISSCROSS;
    symmetricCuboidMesh true;
    numCycles {num_cycles};
    cycleType V;
    fmgInnerCycles {fmg_r}; // 0 == no fmg

    // CRISSCROSS: ~0.4
    // CRISS:    : P1: ? , P2: ~0.72
    sorRelax {omega};
    velocitySorRelax 1.0;

    symmGSVelocity false;
    numGSVelocity {num_gs_velocity};
    symmGSPressure false;
    numGSPressure 1;

    preSmoothingSteps {pre};
    postSmoothingSteps {post};
    smoothingIncrement 2;
    minLevel 2;
    maxLevel {max_level}; // P1 level, P2 level is automatically reduced by 1
    skipCyclesForAvgConvRate 0;
    L2residualTolerance 1e-16;
    projectPressureAfterRestriction true;
    calculateDiscretizationError false;
    coarseGridMaxIterations 2000;
    coarseGridResidualTolerance 1e-04;

    cyclesBeforeDC 0;
    postDCPreSmoothingSteps 3;
    postDCPostSmoothingSteps 3;
    postDCSmoothingIncrement 2;

    outputVTK false;
    outputTiming false;
    outputTimingJSON true;
    outputTimingJSONFile {timing_file};
    outputSQL true;
    outputSQLFile {db_file};
}}
""".format(discretization=discretization, ntan=ntan, nrad=nrad,
           num_cycles=num_cycles, fmg_r=fmg_r, omega=omega, pre=pre,
           post=post, max_level=max_level, timing_file=timing_file,
           db_file=db_file, mesh_spherical_shell=mesh_spherical_shell,
           num_faces_per_side=num_faces_per_side, num_gs_velocity=num_gs_velocity)
    return base_config


def supermuc_job_file_string(job_name="hyteg_job", wall_clock_limit="1:00:00", prm_file="parameter_file.prm", num_nodes=1):

    def partition(num_nodes):
        if num_nodes <= 16:
            return "micro"
        elif num_nodes <= 768:
            return "general"
        elif num_nodes <= 3072:
            return "large"

    base_config = """#!/bin/bash
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
#SBATCH --time={wall_clock_limit}
#SBATCH --no-requeue
#Setup of execution environment
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pr86ma
 
#SBATCH --ear=off
#SBATCH --partition={partition}
#Number of nodes and MPI tasks per node:
#SBATCH --nodes={num_nodes}
#SBATCH --ntasks-per-node=48

module load slurm_setup

source load_modules.sh

module list

cd ..
pwd
ls -lha

#Run the program:
mpiexec -n $SLURM_NTASKS ./MultigridStudies 2019_supermuc/{prm_file}

""".format(job_name=job_name, wall_clock_limit=wall_clock_limit, num_nodes=num_nodes, prm_file=prm_file, partition=partition(num_nodes))
    return base_config


def supermuc_scaling(cube_scaling=True):

    shell_base_config = {"discretization": "P2", "fmg_r": 1, "max_level": 6,
                         "num_cycles": 2, "pre": 6, "post": 6}

    ppn = 48

    node_dep_parameters_shell = {
        1: {"ntan": 3, "nrad":  7, "omega": 0.2},
        2: {"ntan": 3, "nrad": 13, "omega": 0.2}
    }

    cube_base_config = {
        "discretization": "P2",
        "fmg_r": 0,
        "max_level": 8,
        "num_cycles": 10,
        "pre": 4,
        "post": 4,
        "num_gs_velocity": 2,
        "omega": 0.55
    }

    node_dep_parameters_cube = {
        1: {"num_faces_per_side": 1},
        2: {"num_faces_per_side": 2},
        6: {"num_faces_per_side": 3},
        12: {"num_faces_per_side": 4},
        24: {"num_faces_per_side": 5},
        48: {"num_faces_per_side": 6},
        96: {"num_faces_per_side": 8},
        192: {"num_faces_per_side": 10},
        384: {"num_faces_per_side": 13},
        768: {"num_faces_per_side": 16},
        1536: {"num_faces_per_side": 20},
        3072: {"num_faces_per_side": 26},
    }

    if cube_scaling:
        base_config = cube_base_config
        node_dep_parameters = node_dep_parameters_cube

    for num_nodes, prms in node_dep_parameters.items():
        # some_id = str(uuid4())
        some_id = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
        job_name = "mg_studies_{}nodes_{}".format(num_nodes, some_id)
        prm_file_name = job_name + ".prm"
        job_file_name = job_name + ".job"

        timing_file = job_name + ".json"
        db_file = job_name + ".db"

        prm_string_prm_dict = {}
        prm_string_prm_dict.update(base_config)
        prm_string_prm_dict.update(prms)
        prm_string_prm_dict["timing_file"] = timing_file
        prm_string_prm_dict["db_file"] = db_file

        prm_string = supermuc_scaling_prm_file_string(**prm_string_prm_dict)
        job_string = supermuc_job_file_string(job_name=job_name, wall_clock_limit="1:00:00",
                                              num_nodes=num_nodes, prm_file=prm_file_name)

        with open(prm_file_name, "w") as f:
            f.write(prm_string)
        with open(job_file_name, "w") as f:
            f.write(job_string)

if __name__ == "__main__":
    supermuc_scaling()

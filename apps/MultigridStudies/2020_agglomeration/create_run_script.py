import time
import datetime
import os


def supermuc_scaling_prm_file_string(discretization="P2", mesh_spherical_shell=False, num_faces_per_side=1, ntan=2, nrad=2,
                                     num_cycles=1, fmg_r=1, omega=0.2, pre=3, post=3, num_gs_velocity=1, max_level=3,
                                     timing_file="timing.json", agglomeration_timing_file="timing_agglomeration.json",
                                     db_file="database.db", coarse_grid_tol=1e-10,
                                     coarse_grid_solver_type=1, coarse_grid_preconditioner_type=1,
                                     agglomeration=False, agglomeration_strategy='bulk',
                                     agglomeration_num_processes=4, agglomeration_interval=48,
                                     block_low_rank=False, block_low_rank_tolerance=1e-3,
                                     mesh_type='symmetricCube', t_domain_diameter=2, t_domain_height=100,
                                     t_domain_width=25, project_pressure_after_restriction=True):

    base_config = """
Parameters
{{
    equation stokes;
    
    // meshes
    // dim == 2:
    //  square:        square
    // dim == 3:
    //  symmetricCube:  symmetric cube with 24 tets
    //  cube:           cube with 6 tets
    //  sphericalShell: spherical shell
    //  tDomain:        T-shaped domain
    meshType {mesh_type};

    
    numEdgesPerSide {num_faces_per_side};
    discretization {discretization};

    // number of tets = 60 * (ntan-1) * (ntan-1) * (nrad-1)
    shellNTan {ntan};
    shellNRad {nrad};
    shellRMin 0.55;
    shellRMax 1.0;
    
    // parameters for T-domain
    //
    // top down view
    //
    //       out
    //       #
    // in    #     diam   = 1
    // #######     height = 6 // height minus junction
    //       #     width  = 2 // one sided distance from junction
    //       #
    //       out
    //
    //     out
    //     ##
    //     ##
    // in  ##     diam   = 2
    // ######     height = 4 // height minus junction
    // ######     width  = 3 // one sided distance from junction
    //     ##
    //     ##
    //     ##
    //     out
    //
    // height == width == 0 gives only junction with specified diameter (cube w/ diam^3 cubes)
    //
    tDomainDiameter {t_domain_diameter}; // number of cubes in "diameter", > 0 please
    tDomainHeight   {t_domain_height};
    tDomainWidth    {t_domain_width};

    meshLayout CRISSCROSS;
    
    numCycles {num_cycles};
    cycleType V;
    fmgInnerCycles {fmg_r}; // 0 == no fmg

    sorRelax {omega};
    sorRelaxEstimationIterations 20;
    sorRelaxEstimationLevel 3;
    velocitySorRelax 1.0;

    symmGSVelocity false;
    numGSVelocity {num_gs_velocity};
    symmGSPressure false;
    numGSPressure 1;

    preSmoothingSteps {pre};
    postSmoothingSteps {post};
    smoothingIncrement 2;
    minLevel 0;
    maxLevel {max_level}; // P1 level, P2 level is automatically reduced by 1
    skipCyclesForAvgConvRate 0;
    L2residualTolerance 1e-16;
    projectPressureAfterRestriction {project_pressure_after_restriction};
    calculateDiscretizationError false;
    coarseGridMaxIterations 100000;
    coarseGridResidualTolerance {coarse_grid_tol};
    
    // 0: MUMPS                          (PETSc)
    // 1: block preconditioned MINRES    (PETSc)
    // 2: MINRES                         (HyTeG)
    // 3: pressure preconditioned MINRES (HyTeG)
    // 4: SuperLU (dist)                 (PETSc)
    coarseGridSolverType {coarse_grid_solver_type};

    // for solver type 1:
    // 0: PCGAMG
    // 1: PCJACOBI
    // 2: Schur complement
    // 3: HYPRE
    coarseGridSolverVelocityPreconditionerType {coarse_grid_preconditioner_type};

    // BLR options for MUMPS
    blockLowRank {block_low_rank};
    blockLowRankTolerance {block_low_rank_tolerance};

    agglomeration {agglomeration};
    agglomerationStrategy {agglomeration_strategy};   // dedicated, bulk, or interval
    agglomerationNumProcesses {agglomeration_num_processes};
    agglomerationInterval {agglomeration_interval};

    cyclesBeforeDC 0;
    postDCPreSmoothingSteps 3;
    postDCPostSmoothingSteps 3;
    postDCSmoothingIncrement 2;

    outputBaseDirectory /hppfs/work/pr86ma/di36vuv2;
    outputVTK false;
    outputTiming false;
    outputTimingJSON true;
    outputTimingJSONFile {timing_file};
    outputAgglomerationTimingJSONFile {agglomeration_timing_file};
    outputSQL true;
    outputParallelSQL false;
    outputSQLFile {db_file};
}}
""".format(discretization=discretization, ntan=ntan, nrad=nrad,
           num_cycles=num_cycles, fmg_r=fmg_r, omega=omega, pre=pre,
           post=post, max_level=max_level, timing_file=timing_file, agglomeration_timing_file=agglomeration_timing_file,
           db_file=db_file, mesh_spherical_shell=mesh_spherical_shell,
           num_faces_per_side=num_faces_per_side, num_gs_velocity=num_gs_velocity,
           coarse_grid_tol=coarse_grid_tol, coarse_grid_solver_type=coarse_grid_solver_type,
           coarse_grid_preconditioner_type=coarse_grid_preconditioner_type,
           agglomeration=agglomeration, agglomeration_strategy=agglomeration_strategy,
           agglomeration_num_processes=agglomeration_num_processes, agglomeration_interval=agglomeration_interval,
           block_low_rank=block_low_rank, block_low_rank_tolerance=block_low_rank_tolerance,
           mesh_type=mesh_type, t_domain_diameter=t_domain_diameter, t_domain_height=t_domain_height,
           t_domain_width=t_domain_width, project_pressure_after_restriction=project_pressure_after_restriction)
    return base_config


def supermuc_job_file_string(job_name="hyteg_job", wall_clock_limit="1:00:00", prm_file="parameter_file.prm", num_nodes=1, ppn=48, petsc_detail=False, script_dir="."):

    petsc_detail_string = ""
    if petsc_detail:
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
#SBATCH --ntasks-per-node={ppn}
{constraint}

cd ..

module load slurm_setup

source load_modules.sh

module list

cd ..
pwd
ls -lha

#Run the program:
mpiexec -n $SLURM_NTASKS -print-rank-map ./MultigridStudies 2020_agglomeration/{script_dir}/{prm_file} {petsc_detail_string}

""".format(job_name=job_name, wall_clock_limit=wall_clock_limit, num_nodes=num_nodes, prm_file=prm_file, partition=partition(num_nodes),
           constraint=constraint, ppn=ppn, petsc_detail_string=petsc_detail_string, script_dir=script_dir)
    return base_config


def supermuc_scaling():

    cube_base_config_fmg = {
        "weak": {
            "P2": {
                "discretization": "P2",
                "fmg_r": 0,
                "max_level": 8,
                "num_cycles": 5,
                "pre": 3,
                "post": 3,
                "num_gs_velocity": 2,
                "omega": 0.65,
            }
        },

        "strong": {
            "P2": {
                "discretization": "P2",
                "fmg_r": 1,
                "max_level": 6,
                "num_cycles": 1,
                "pre": 3,
                "post": 3,
                "num_gs_velocity": 2,
                "omega": 0.65,
            }
        }
    }

    cube_base_config_fmg["weak_fast"] = cube_base_config_fmg["weak"]
    cube_base_config_fmg["weak_large"] = cube_base_config_fmg["weak"]
    cube_base_config_fmg["weak_large_3072_27fps"] = cube_base_config_fmg["weak"]
    cube_base_config_fmg["weak_large_3072_28fps"] = cube_base_config_fmg["weak"]
    cube_base_config_fmg["weak_large_3072_29fps"] = cube_base_config_fmg["weak"]
    cube_base_config_fmg["weak_fast_tdomain"] = cube_base_config_fmg["weak"]

    node_dep_parameters_cube = {
        "weak": {
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
        },
        "weak_large": {
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
        },
        "weak_large_3072_27fps": {
            3072: {"num_faces_per_side": 27},
        },
        "weak_large_3072_28fps": {
            3072: {"num_faces_per_side": 28},
        },
        "weak_large_3072_29fps": {
            3072: {"num_faces_per_side": 29},
        },
        "weak_fast": {
            2: {"num_faces_per_side": 1},
            6: {"num_faces_per_side": 2},
            12: {"num_faces_per_side": 3},
            24: {"num_faces_per_side": 4},
            48: {"num_faces_per_side": 5},
            96: {"num_faces_per_side": 6},
            192: {"num_faces_per_side": 8},
            384: {"num_faces_per_side": 10},
            768: {"num_faces_per_side": 13},
            1536: {"num_faces_per_side": 16},
            3072: {"num_faces_per_side": 20},
        },
        "strong": {
            1: {"num_faces_per_side": 6},
            2: {"num_faces_per_side": 6},
            6: {"num_faces_per_side": 6},
            12: {"num_faces_per_side": 6},
            24: {"num_faces_per_side": 6},
            48: {"num_faces_per_side": 6},
            96: {"num_faces_per_side": 6},
        },
        "weak_fast_tdomain": {
            192: {"mesh_type": 'tDomain', "t_domain_diameter": 2, "t_domain_height": 100, "t_domain_width": 25,
                  "project_pressure_after_restriction": False},
            384: {"mesh_type": 'tDomain', "t_domain_diameter": 2, "t_domain_height": 200, "t_domain_width": 50,
                  "project_pressure_after_restriction": False},
            768: {"mesh_type": 'tDomain', "t_domain_diameter": 2, "t_domain_height": 400, "t_domain_width": 100,
                  "project_pressure_after_restriction": False},
            1536: {"mesh_type": 'tDomain', "t_domain_diameter": 2, "t_domain_height": 800, "t_domain_width": 200,
                   "project_pressure_after_restriction": False},
            3072: {"mesh_type": 'tDomain', "t_domain_diameter": 2, "t_domain_height": 1600, "t_domain_width": 400,
                   "project_pressure_after_restriction": False},
        }
    }

    some_id = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
    os.mkdir(some_id)

    coarse_grid_solver_string = {
        (0, 0): "MUMPS",
        (1, 1): "MinRes-Jacobi",
        (1, 3): "MinRes-HYPRE",
        (4, 0): "SuperLU"
    }

    blr_settings = {
        0: [(False, 0), (True, 1e-3), (True, 1e-5)],
        1: [(False, 0)],
        4: [(False, 0)],
    }

    agglomeration_parameters = [
        {"agglomeration": False},
        {"agglomeration": True, "agglomeration_strategy": "bulk", "agglomeration_num_processes": 48},
        {"agglomeration": True, "agglomeration_strategy": "interval", "agglomeration_num_processes": 48},
        {"agglomeration": True, "agglomeration_strategy": "dedicated", "agglomeration_num_processes": 48},
        {"agglomeration": True, "agglomeration_strategy": "bulk", "agglomeration_num_processes": 96},
        {"agglomeration": True, "agglomeration_strategy": "interval", "agglomeration_num_processes": 96},
        {"agglomeration": True, "agglomeration_strategy": "dedicated", "agglomeration_num_processes": 96},
        {"agglomeration": True, "agglomeration_strategy": "bulk", "agglomeration_num_processes": 192},
        {"agglomeration": True, "agglomeration_strategy": "interval", "agglomeration_num_processes": 192},
        {"agglomeration": True, "agglomeration_strategy": "dedicated", "agglomeration_num_processes": 192},
        {"agglomeration": True, "agglomeration_strategy": "bulk", "agglomeration_num_processes": 384},
        {"agglomeration": True, "agglomeration_strategy": "interval", "agglomeration_num_processes": 384},
        {"agglomeration": True, "agglomeration_strategy": "dedicated", "agglomeration_num_processes": 384},
    ]

    for discretization in ["P2"]:
        for scaling_type in ["weak_fast_tdomain"]:
            for coarse_grid_tol in [1e-05, 1e-12]:
                for coarse_grid_solver_type, coarse_grid_preconditioner_type in coarse_grid_solver_string.keys():
                    for agglomeration_parameter_set in agglomeration_parameters:
                        for blr, blr_tol in blr_settings[coarse_grid_solver_type]:
                            if "weak_large" in scaling_type:
                                ppn = 24
                            else:
                                ppn = 48

                            base_config = cube_base_config_fmg[scaling_type][discretization]
                            base_config["coarse_grid_tol"] = coarse_grid_tol
                            base_config["coarse_grid_solver_type"] = coarse_grid_solver_type
                            base_config["coarse_grid_preconditioner_type"] = coarse_grid_preconditioner_type
                            base_config["block_low_rank"] = blr
                            base_config["block_low_rank_tolerance"] = blr_tol
                            base_config.update(agglomeration_parameter_set)
                            agglomeration_string = "agg_" + (agglomeration_parameter_set["agglomeration_strategy"] + \
                                                             str(agglomeration_parameter_set["agglomeration_num_processes"]) if agglomeration_parameter_set["agglomeration"] else "none")
                            node_dep_parameters = node_dep_parameters_cube[scaling_type]

                            blr_string = ""
                            if coarse_grid_solver_type == 0:
                                blr_string = "_BLR_{:.2e}".format(blr_tol)

                            for num_nodes, prms in node_dep_parameters.items():
                                # some_id = str(uuid4())

                                job_name = "mg_studies_{}_{}_cgstype_{}_cgtol_{:.2e}_{}nodes_{}_{}".format(
                                    scaling_type, discretization, coarse_grid_solver_string[(coarse_grid_solver_type, coarse_grid_preconditioner_type)] + blr_string,
                                    coarse_grid_tol, num_nodes, agglomeration_string, some_id)
                                prm_file_name = job_name + ".prm"
                                job_file_name = job_name + ".job"

                                timing_file = job_name + ".json"
                                agglomeration_timing_file = job_name + "_agglomeration.json"
                                db_file = job_name + ".db"

                                prm_string_prm_dict = {}
                                prm_string_prm_dict.update(base_config)
                                prm_string_prm_dict.update(prms)
                                prm_string_prm_dict["timing_file"] = timing_file
                                prm_string_prm_dict["agglomeration_timing_file"] = agglomeration_timing_file
                                prm_string_prm_dict["db_file"] = db_file

                                prm_string = supermuc_scaling_prm_file_string(**prm_string_prm_dict)
                                job_string = supermuc_job_file_string(job_name=job_name, wall_clock_limit="0:30:00",
                                                                      num_nodes=num_nodes, prm_file=prm_file_name, ppn=ppn, petsc_detail=True,
                                                                      script_dir=some_id)

                                with open(os.path.join(some_id, prm_file_name), "w") as f:
                                    f.write(prm_string)
                                with open(os.path.join(some_id, job_file_name), "w") as f:
                                    f.write(job_string)

if __name__ == "__main__":
    supermuc_scaling()

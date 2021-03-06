import time
import datetime

def parameter_string(vtkDirectory, vtkBaseName, ntan, num_layers):

    if num_layers == 3:
        layer_string = """  innerBoundary 0.55; // 3481 km from center (outer core boundary)
  layer1 0.775;
  outerBoundary 1.0; // 6371 km from center"""
    elif num_layers == 5:
        layer_string = """  innerBoundary 0.55; // 3481 km from center (outer core boundary)
  layer0 0.6625;
  layer1 0.775;
  layer2 0.8875;
  outerBoundary 1.0; // 6371 km from center"""

    return """Parameters
{{
  // domain
  ntan {ntan};
  minLevel 2;
  maxLevel 5;

  // simulation
  timeSteps 10000;
  diffusivity -1; // set negative to optimize for zero diffusivity
  dt 1e-2;
  rhsScaleFactor 100.0;

  // Stokes solver
  stokesResidual 1e-12; // after each timestep we iterate until we reach this residual
  stokesMaxNumVCycles 5;
  stokesSolveInterval 10;
  numDiffusionVCycles 1;

  /// output
  printTiming true;
  timingFile {vtkBaseName}.json;
  writeDomainVTK true;
  vtkOutput true;
  vtkFrequency 10;
  vtkBaseFile {vtkBaseName};
  vtkDirectory {vtkDirectory};
  vtkOutputLevel 4; // P1 level, max output level <= maxLevel but will output only the vertex unknowns (cannot output all DoFs)
  exitAfterWriteDomain false; // true to exit application after coarse domain output (to improve mesh offline)
}}

/// Layers for the spherical shell generator
/// the keys can be arbitrary but need to different
/// the values have to be sorted ascending
Layers
{{
{layer_string}
}}
""".format(vtkBaseName=vtkBaseName, vtkDirectory=vtkDirectory, ntan=ntan, layer_string=layer_string)


def supermuc_job_file_string(job_name="hyteg_job", wall_clock_limit="1:00:00", prm_file="parameter_file.prm", num_nodes=1, ppn=48):

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
 
#SBATCH --partition={partition}
#Number of nodes and MPI tasks per node:
#SBATCH --nodes={num_nodes}
#SBATCH --ntasks-per-node={ppn}
{constraint}

module load slurm_setup

source load_modules_supermucng.sh

module list

mkdir -p /hppfs/work/pr86ma/di36vuv2/stokes_sphere_transport

pwd
ls -lha

#Run the program:
mpiexec -n $SLURM_NTASKS ./StokesSphereTransport {prm_file}

""".format(job_name=job_name, wall_clock_limit=wall_clock_limit, num_nodes=num_nodes, prm_file=prm_file, partition=partition(num_nodes),
           constraint=constraint, ppn=ppn)
    return base_config


def generate_files():

    some_id = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')

    config_0 = {
        "name": "config_8_nodes_1920cells",
        "num_nodes": 8,
        "ntan": 5,
        "num_layers": 3
    }

    config_1 = {
        "name": "config_96nodes_15360cells",
        "num_nodes": 96,
        "ntan": 9,
        "num_layers": 5
    }

    for config in [config_0, config_1]:
        job_name = "stokes_{}_{}".format(config["name"], some_id)

        prm_string = parameter_string("/hppfs/work/pr86ma/di36vuv2/stokes_sphere_transport/", job_name, config["ntan"], config["num_layers"])
        prm_file_name = job_name + ".prm"
        job_file_name = job_name + ".job"

        job_string = supermuc_job_file_string(job_name=job_name, wall_clock_limit="4:00:00",
                                              num_nodes=config["num_nodes"], prm_file=prm_file_name, ppn=48)

        with open(prm_file_name, "w") as f:
            f.write(prm_string)
        with open(job_file_name, "w") as f:
            f.write(job_string)

if __name__ == "__main__":
    generate_files()

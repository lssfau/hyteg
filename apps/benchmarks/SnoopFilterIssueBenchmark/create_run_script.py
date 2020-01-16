import time
import datetime
import os

def supermuc_scaling_prm_file_string(level=7, num_cells_per_process=1, num_outer_iterations=1, num_inner_iterations=1, db_file="database.db"):

    base_config = """
Parameters
{{
    level {};
    numCellPerProcess {};
    numOuterSORIterations {};
    numInnerSORIterations {};
    dbFile {};
}}
""".format(level, num_cells_per_process, num_outer_iterations, num_inner_iterations, db_file)
    return base_config


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
 
#SBATCH --ear=off
#SBATCH --partition={partition}
#Number of nodes and MPI tasks per node:
#SBATCH --nodes={num_nodes}
#SBATCH --ntasks-per-node={ppn}
{constraint}

module load slurm_setup

source load_modules.sh

module list

cd ..
pwd
ls -lha

#Run the program:
mpiexec -n $SLURM_NTASKS ./SnoopFilterIssueBenchmark {prm_file}

""".format(job_name=job_name, wall_clock_limit=wall_clock_limit, num_nodes=num_nodes, prm_file=prm_file, partition=partition(num_nodes),
           constraint=constraint, ppn=ppn)
    return base_config


def supermuc_scaling():

    some_id = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
    os.mkdir(some_id)

    ppn = 48
    outer_iterations = 2
    inner_iterations = 10
    cells_pp = 2
    level = 7

    for num_nodes in [1, 2, 6, 12, 24, 48, 96, 192, 384, 768, 1536, 3072]:
        job_name = "snoop_n{}_ppn{}_level{}_cellspp{}_inner{}_outer{}_{}".format(num_nodes, ppn, level, cells_pp, inner_iterations, outer_iterations, some_id)
        prm_file_name = job_name + ".prm"
        job_file_name = job_name + ".job"
        db_file = job_name + ".db"

        prm_string = supermuc_scaling_prm_file_string(level=level, num_cells_per_process=cells_pp, num_outer_iterations=outer_iterations,
                                                      num_inner_iterations=inner_iterations, db_file=os.path.join(some_id, db_file))
        job_string = supermuc_job_file_string(job_name=job_name, wall_clock_limit="0:10:00",
                                              num_nodes=num_nodes, prm_file=os.path.join(some_id, prm_file_name), ppn=ppn)

        with open(os.path.join(some_id, prm_file_name), "w") as f:
            f.write(prm_string)
        with open(os.path.join(some_id, job_file_name), "w") as f:
            f.write(job_string)

if __name__ == "__main__":
    supermuc_scaling()

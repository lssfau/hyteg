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


def emmy_job_file_string(job_name="hyteg_job", wall_clock_limit="1:00:00", prm_string="parameter_file.prm", num_nodes=1, ppn=20):

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

    base_config = """#!/bin/bash -l
#
#PBS -l nodes={num_nodes}:ppn=40,walltime={wall_clock_limit}
#
# job name 
#PBS -N {job_name}
#
# first non-empty non-comment line ends PBS options

# jobs always start in $HOME - 
# change to work directory
cd /home/hpc/iwia/iwia007h/hyteg-build-release/apps/benchmarks/SnoopFilterIssueBenchmark

# load required modules (compiler, MPI, ...)
source load_modules_emmy.sh

# run, using only physical cores
mpirun -n {num_cores} ./SnoopFilterIssueBenchmark {prm_string}

""".format(job_name=job_name, wall_clock_limit=wall_clock_limit, num_nodes=num_nodes, prm_string=prm_string, partition=partition(num_nodes),
           num_cores=ppn*num_nodes)
    return base_config


def emmy_scaling():

    some_id = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
    os.mkdir(some_id)

    ppn = 20
    outer_iterations = 2
    inner_iterations_list = [1, 10, 100]
    cells_pp_list = [2, 10]
    level_list = [6, 7]

    for num_nodes in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]:
        for inner_iterations in inner_iterations_list:
            for cells_pp in cells_pp_list:
                for level in level_list:
                    job_name = "snoop_n{}_ppn{}_level{}_cellspp{}_inner{}_outer{}_{}".format(num_nodes, ppn, level, cells_pp, inner_iterations, outer_iterations, some_id)
                    prm_file_name = job_name + ".prm"
                    job_file_name = job_name + ".job"
                    db_file = job_name + ".csv"

                    # prm_string = supermuc_scaling_prm_file_string(level=level, num_cells_per_process=cells_pp, num_outer_iterations=outer_iterations,
                    #                                               num_inner_iterations=inner_iterations, db_file=os.path.join(some_id, db_file))
                    prm_string = "{} {} {} {} {}".format(level, cells_pp, outer_iterations, inner_iterations, db_file)
                    job_string = emmy_job_file_string(job_name=job_name, wall_clock_limit="0:15:00",
                                                      num_nodes=num_nodes, prm_string=prm_string, ppn=ppn)

                    # with open(os.path.join(some_id, prm_file_name), "w") as f:
                        # f.write(prm_string)
                    with open(os.path.join(some_id, job_file_name), "w") as f:
                        f.write(job_string)

if __name__ == "__main__":
    emmy_scaling()

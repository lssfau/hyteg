
import datetime
import time
import os
import subprocess

executables = {
    1: 'Benchmark_01_Cube', # discr error
    2: 'Benchmark_01_Cube'  # FMG table, and convergence plot, p2p1
}

parameterizations_base = {
    1: {'discretization': 'p2p1', 'minLevel': 0, 'numEdgesPerSide': 1, 'vtk': False, 'numCycles': 100, 'absoluteResidualTolerance': 1e-12,
        'estimateOmega': True, 'omegaEstimationIterations': 20,
        'preSmooth': 10, 'postSmooth': 10, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True,
        'coarseGridSolverType': 0, 'normCalculationLevelIncrement': 1},

    2: {'discretization': 'p2p1', 'minLevel': 0, 'maxLevel': 6, 'numEdgesPerSide': 1, 'vtk': False, 'numCycles': 1, 'absoluteResidualTolerance': 1e-12,
        'estimateOmega': True, 'omegaEstimationLevel': 6, 'omegaEstimationIterations': 20,
        'coarseGridSolverType': 0, 'normCalculationLevelIncrement': 1},
}

parameterizations = {
    1: [
        {'maxLevel': 0, 'omegaEstimationLevel': 0},
        {'maxLevel': 1, 'omegaEstimationLevel': 1},
        {'maxLevel': 2, 'omegaEstimationLevel': 2},
        {'maxLevel': 3, 'omegaEstimationLevel': 3},
        {'maxLevel': 4, 'omegaEstimationLevel': 4},
        {'maxLevel': 5, 'omegaEstimationLevel': 5},
        {'maxLevel': 6, 'omegaEstimationLevel': 6},
    ],
    2: [
        {'preSmooth': 0, 'postSmooth': 2, 'incSmooth': 1, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': False},
        {'preSmooth': 1, 'postSmooth': 1, 'incSmooth': 1, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': False},
        {'preSmooth': 2, 'postSmooth': 2, 'incSmooth': 1, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': False},
        {'preSmooth': 3, 'postSmooth': 3, 'incSmooth': 1, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': False},
        {'preSmooth': 0, 'postSmooth': 2, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': False},
        {'preSmooth': 1, 'postSmooth': 1, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': False},
        {'preSmooth': 2, 'postSmooth': 2, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': False},
        {'preSmooth': 3, 'postSmooth': 3, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': False},

        {'preSmooth': 0, 'postSmooth': 2, 'incSmooth': 1, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True},
        {'preSmooth': 1, 'postSmooth': 1, 'incSmooth': 1, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True},
        {'preSmooth': 2, 'postSmooth': 2, 'incSmooth': 1, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True},
        {'preSmooth': 3, 'postSmooth': 3, 'incSmooth': 1, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True},
        {'preSmooth': 0, 'postSmooth': 2, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True},
        {'preSmooth': 1, 'postSmooth': 1, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True},
        {'preSmooth': 2, 'postSmooth': 2, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True},
        {'preSmooth': 3, 'postSmooth': 3, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True},
    ],
}


def run_all_configs(benchmark_id):

    print("Running benchmark " + str(benchmark_id))

    date_id = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
    id = '_'.join([date_id, 'benchmark', str(benchmark_id)])

    # create new db dir
    db_dir = "db/db_" + id
    print("Creating directory for databases: " + db_dir)
    os.mkdir(db_dir)

    for run_id, parameterization in enumerate(parameterizations[benchmark_id]):
        db_file = os.path.join(db_dir, '_'.join([id, str(run_id)]) + '.db')
        cmd = 'mpirun -np 4 ./{exe} {exe}.prm '.format(exe=executables[benchmark_id])
        cmd += '-Parameters.dbFile={db_file} '.format(db_file=db_file)
        for prm_key, prm_val in parameterizations_base[benchmark_id].items():
            cmd += '-Parameters.{prm_key}={prm_val} '.format(prm_key=prm_key, prm_val=prm_val)
        for prm_key, prm_val in parameterization.items():
            cmd += '-Parameters.{prm_key}={prm_val} '.format(prm_key=prm_key, prm_val=prm_val)
        print(cmd)
        completed_process = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE, universal_newlines=True)
        out_file = os.path.join(db_dir, '_'.join([id, str(run_id)]) + '.out')
        with open(out_file, 'w') as f:
            f.write(completed_process.stdout)

run_all_configs(1)
run_all_configs(2)

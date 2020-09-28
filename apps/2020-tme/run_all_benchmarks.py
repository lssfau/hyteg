
import datetime
import time
import os
import subprocess

executables = {
    1: 'Benchmark_01_Cube', # discr error P2-P1
    2: 'Benchmark_01_Cube', # discr error P1-P1
    3: 'Benchmark_01_Cube', # FMG table, and convergence plot, P2-P1
    4: 'Benchmark_01_Cube', # v-cycles, P2-P1
    5: 'Benchmark_01_Cube', # FMG table, and convergence plot, P1-P1
    6: 'Benchmark_02_Y-Pipe',
}

# omega
# P2-P1, level 5, 100 iterations: 0.448872
# P1-P1, level 6, 100 iterations: 0.570751 -> bad choice, v-cycles diverge, but seems to work for FMG
# P1-P1, level 6, guess:          0.3      -> seems to work for FMG and v-cycles

parameterizations_base = {
    1: {'discretization': 'p2p1', 'minLevel': 0, 'numEdgesPerSide': 1, 'vtk': False, 'numCycles': 100, 'absoluteResidualTolerance': 1e-12,
        'estimateOmega': True, 'omegaEstimationIterations': 20,
        'preSmooth': 10, 'postSmooth': 10, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True,
        'coarseGridSolverType': 0, 'normCalculationLevelIncrement': 1},

    2: {'discretization': 'p1p1', 'minLevel': 0, 'numEdgesPerSide': 1, 'vtk': False, 'numCycles': 100, 'absoluteResidualTolerance': 1e-12,
        'estimateOmega': False, 'omega': 0.3,
        'preSmooth': 10, 'postSmooth': 10, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True,
        'coarseGridSolverType': 0, 'normCalculationLevelIncrement': 1},

    3: {'discretization': 'p2p1', 'minLevel': 0, 'maxLevel': 5, 'numEdgesPerSide': 1, 'vtk': False, 'numCycles': 1, 'absoluteResidualTolerance': 1e-12,
        'estimateOmega': False, 'omega': 0.448872, 'omegaEstimationLevel': 6, 'omegaEstimationIterations': 20,
        'coarseGridSolverType': 0, 'normCalculationLevelIncrement': 1},

    4: {'discretization': 'p2p1', 'minLevel': 0, 'maxLevel': 5, 'numEdgesPerSide': 1, 'vtk': False, 'numCycles': 100, 'absoluteResidualTolerance': 1e-12,
        'estimateOmega': False, 'omega': 0.448872, 'omegaEstimationLevel': 6, 'omegaEstimationIterations': 20,
        'coarseGridSolverType': 0, 'normCalculationLevelIncrement': 1, 'fmgInnerIterations': 0},

    5: {'discretization': 'p1p1', 'minLevel': 0, 'maxLevel': 6, 'numEdgesPerSide': 1, 'vtk': False, 'numCycles': 1, 'absoluteResidualTolerance': 1e-12,
        'estimateOmega': False, 'omega': 0.570751, 'omegaEstimationLevel': 6, 'omegaEstimationIterations': 20,
        'coarseGridSolverType': 0, 'normCalculationLevelIncrement': 1},

    6: {'discretization': 'p2p1', 'minLevel': 0, 'maxLevel': 4, 'numEdgesPerSide': 1, 'vtk': False, 'numCycles': 1,
        'absoluteResidualTolerance': 1e-12,
        'estimateOmega': False, 'omega': 0.2,
        'coarseGridSolverType': 0, 'normCalculationLevelIncrement': 1,
        'calculateDiscretizationError': True,
        'DiscretizationErrorSolver.preSmooth': 7,
        'DiscretizationErrorSolver.postSmooth': 7,
        'DiscretizationErrorSolver.incSmooth': 2,
        'DiscretizationErrorSolver.fmgInnerIterations': 2,
        'DiscretizationErrorSolver.numCycles': 50,
        'DiscretizationErrorSolver.absoluteResidualTolerance': 1e-12,
        'DiscretizationErrorSolver.numGSVelocity': 4,
        'DiscretizationErrorSolver.symmGSVelocity': False,
        'DiscretizationErrorSolver.estimateOmega': False,
        'DiscretizationErrorSolver.omega': 0.2,
        'DiscretizationErrorSolver.coarseGridSolverType': 0}
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
        {'maxLevel': 0},
        {'maxLevel': 1},
        {'maxLevel': 2},
        {'maxLevel': 3},
        {'maxLevel': 4},
        {'maxLevel': 5},
        {'maxLevel': 6},
        {'maxLevel': 7},
    ],
    3: [],
    4: [],
    5: [],

    6: [{'preSmooth': 1, 'postSmooth': 3, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 3, 'symmGSVelocity': False},
        {'preSmooth': 1, 'postSmooth': 2, 'incSmooth': 1, 'fmgInnerIterations': 1, 'numGSVelocity': 3, 'symmGSVelocity': False},
        {'preSmooth': 0, 'postSmooth': 2, 'incSmooth': 1, 'fmgInnerIterations': 1, 'numGSVelocity': 3, 'symmGSVelocity': False},]
}

for pre in range(4):
    for post in range(4):
        for inc in range(4):
            for kappa in range(1, 3):
                for gs, symm in [(1, True), (2, True), (1, False), (2, False), (3, False), (4, False)]:
                    parameterizations[3].append({'preSmooth': pre, 'postSmooth': post, 'incSmooth': inc, 'fmgInnerIterations': kappa, 'numGSVelocity': gs, 'symmGSVelocity': symm})
                    parameterizations[5].append({'preSmooth': pre, 'postSmooth': post, 'incSmooth': inc, 'fmgInnerIterations': kappa, 'numGSVelocity': gs, 'symmGSVelocity': symm})

for pre in range(4):
    for post in range(4):
        for inc in range(4):
            for gs, symm in [(1, True), (1, False), (2, False), (3, False)]:
                parameterizations[4].append({'preSmooth': pre, 'postSmooth': post, 'incSmooth': inc, 'numGSVelocity': gs, 'symmGSVelocity': symm})

def run_all_configs(benchmark_id):

    print("Running benchmark " + str(benchmark_id))

    date_id = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
    id = '_'.join([date_id, 'benchmark', str(benchmark_id)])

    # create new db dir
    db_dir = "db/db_" + id
    print("Creating directory for databases: " + db_dir)
    os.mkdir(db_dir)

    print("Running ", len(parameterizations[benchmark_id]), " benchmarks.")

    for run_id, parameterization in enumerate(parameterizations[benchmark_id]):
        db_file = os.path.join(db_dir, '_'.join([id, str(run_id)]) + '.db')
        cmd = 'mpirun -np 4 ./{exe} {exe}.prm '.format(exe=executables[benchmark_id])
        cmd += '-Parameters.dbFile={db_file} '.format(db_file=db_file)
        if benchmark_id == 6:
            cmd += '-Parameters.DiscretizationErrorSolver.dbFile={db_file} '.format(db_file=db_file + 'discr.db')
        for prm_key, prm_val in parameterizations_base[benchmark_id].items():
            cmd += '-Parameters.{prm_key}={prm_val} '.format(prm_key=prm_key, prm_val=prm_val)
        for prm_key, prm_val in parameterization.items():
            cmd += '-Parameters.{prm_key}={prm_val} '.format(prm_key=prm_key, prm_val=prm_val)
        print(cmd)
        completed_process = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE, universal_newlines=True)
        out_file = os.path.join(db_dir, '_'.join([id, str(run_id)]) + '.out')
        with open(out_file, 'w') as f:
            f.write(completed_process.stdout)

# run_all_configs(1)
# run_all_configs(2)
# run_all_configs(3)
# run_all_configs(4)
# run_all_configs(5)
run_all_configs(6)

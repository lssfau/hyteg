
import datetime
import time
import os
import subprocess
import numpy as np

executables = {

    # discretization errors

    'discr-cube-p2p1': 'Benchmark_01_Cube',
    'discr-cube-p1p1': 'Benchmark_01_Cube',

    'discr-pipe-p2p1': 'Benchmark_02_Y-Pipe',
    'discr-pipe-p1p1': 'Benchmark_02_Y-Pipe',

    # fmg studies

    'fmg-cube-p2p1': 'Benchmark_01_Cube',
    'fmg-cube-p1p1': 'Benchmark_01_Cube',

    'fmg-pipe-p2p1': 'Benchmark_02_Y-Pipe',
    'fmg-pipe-p1p1': 'Benchmark_02_Y-Pipe',
}

for prm in ['omega-cube-p1p1-1-1-2-1symm', 'omega-cube-p1p1-2-2-2-1symm', 'omega-cube-p1p1-3-3-2-1symm', 'omega-cube-p1p1-1-1-2-3fwd', 'omega-cube-p1p1-2-2-2-3fwd', 'omega-cube-p1p1-3-3-2-3fwd',
            'omega-cube-p2p1-1-1-2-1symm', 'omega-cube-p2p1-2-2-2-1symm', 'omega-cube-p2p1-3-3-2-1symm', 'omega-cube-p2p1-1-1-2-3fwd', 'omega-cube-p2p1-2-2-2-3fwd', 'omega-cube-p2p1-3-3-2-3fwd']:
    executables[prm] = 'Benchmark_01_Cube'

for prm in ['omega-pipe-p1p1-1-1-2-1symm', 'omega-pipe-p1p1-2-2-2-1symm', 'omega-pipe-p1p1-3-3-2-1symm', 'omega-pipe-p1p1-1-1-2-3fwd', 'omega-pipe-p1p1-2-2-2-3fwd', 'omega-pipe-p1p1-3-3-2-3fwd',
            'omega-pipe-p2p1-1-1-2-1symm', 'omega-pipe-p2p1-2-2-2-1symm', 'omega-pipe-p2p1-3-3-2-1symm', 'omega-pipe-p2p1-1-1-2-3fwd', 'omega-pipe-p2p1-2-2-2-3fwd', 'omega-pipe-p2p1-3-3-2-3fwd']:
    executables[prm] = 'Benchmark_02_Y-Pipe'

# omega

# from omega sweep:
# cube:
# - P1-P1-stab: 0.1
# - P2-P1:      0.2
# pipe:
# - P1-P1-stab: 0.3
# - P2-P1:      0.4

# old:
# P2-P1, level 5, 100 iterations: 0.448872
# P1-P1, level 6, 100 iterations: 0.570751 -> bad choice, v-cycles diverge, but seems to work for FMG
# P1-P1, level 6, guess:          0.3      -> seems to work for FMG and v-cycles

omega_base = {
    'minLevel': 0,
    'numEdgesPerSide': 1,
    'vtk': False,
    'numCycles': 1,
    'absoluteResidualTolerance': 1e-12,
    'estimateOmega': False,
    'coarseGridSolverType': 1,
    'normCalculationLevelIncrement': 1,
    'calculateDiscretizationError': False
}

discr_base =  {'minLevel': 0, 'numEdgesPerSide': 1, 'vtk': False, 'numCycles': 100, 'absoluteResidualTolerance': 1e-12,
               'estimateOmega': False, 'preSmooth': 10, 'postSmooth': 10, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 4, 'symmGSVelocity': False,
               'coarseGridSolverType': 1, 'normCalculationLevelIncrement': 1}

parameterizations_base = {

    # Omega estiamtion cube

    'omega-cube-p1p1-1-1-2-1symm' : {**omega_base, **{
        'discretization': 'p1p1', 'maxLevel': 5, 'preSmooth': 1, 'postSmooth': 1, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True,
    }},

    'omega-cube-p1p1-2-2-2-1symm' : {**omega_base, **{
        'discretization': 'p1p1', 'maxLevel': 5, 'preSmooth': 2, 'postSmooth': 2, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True,
    }},

    'omega-cube-p1p1-3-3-2-1symm' : {**omega_base, **{
        'discretization': 'p1p1', 'maxLevel': 5, 'preSmooth': 3, 'postSmooth': 3, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True,
    }},

    'omega-cube-p1p1-1-1-2-3fwd' : {**omega_base, **{
        'discretization': 'p1p1', 'maxLevel': 5, 'preSmooth': 1, 'postSmooth': 1, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 3, 'symmGSVelocity': False,
    }},

    'omega-cube-p1p1-2-2-2-3fwd' : {**omega_base, **{
        'discretization': 'p1p1', 'maxLevel': 5, 'preSmooth': 2, 'postSmooth': 2, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 3, 'symmGSVelocity': False,
    }},

    'omega-cube-p1p1-3-3-2-3fwd' : {**omega_base, **{
        'discretization': 'p1p1', 'maxLevel': 5, 'preSmooth': 3, 'postSmooth': 3, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 3, 'symmGSVelocity': False,
    }},

    'omega-cube-p2p1-1-1-2-1symm' : {**omega_base, **{
        'discretization': 'p2p1', 'maxLevel': 4, 'preSmooth': 1, 'postSmooth': 1, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True,
    }},

    'omega-cube-p2p1-2-2-2-1symm' : {**omega_base, **{
        'discretization': 'p2p1', 'maxLevel': 4, 'preSmooth': 2, 'postSmooth': 2, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True,
    }},

    'omega-cube-p2p1-3-3-2-1symm' : {**omega_base, **{
        'discretization': 'p2p1', 'maxLevel': 4, 'preSmooth': 3, 'postSmooth': 3, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True,
    }},

    'omega-cube-p2p1-1-1-2-3fwd' : {**omega_base, **{
        'discretization': 'p2p1', 'maxLevel': 4, 'preSmooth': 1, 'postSmooth': 1, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 3, 'symmGSVelocity': False,
    }},

    'omega-cube-p2p1-2-2-2-3fwd' : {**omega_base, **{
        'discretization': 'p2p1', 'maxLevel': 4, 'preSmooth': 2, 'postSmooth': 2, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 3, 'symmGSVelocity': False,
    }},

    'omega-cube-p2p1-3-3-2-3fwd' : {**omega_base, **{
        'discretization': 'p2p1', 'maxLevel': 4, 'preSmooth': 3, 'postSmooth': 3, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 3, 'symmGSVelocity': False,
    }},


    # Omega estimation y-pipe

    'omega-pipe-p1p1-1-1-2-1symm' : {**omega_base, **{
        'discretization': 'p1p1', 'maxLevel': 4, 'preSmooth': 1, 'postSmooth': 1, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True,
    }},

    'omega-pipe-p1p1-2-2-2-1symm' : {**omega_base, **{
        'discretization': 'p1p1', 'maxLevel': 4, 'preSmooth': 2, 'postSmooth': 2, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True,
    }},

    'omega-pipe-p1p1-3-3-2-1symm' : {**omega_base, **{
        'discretization': 'p1p1', 'maxLevel': 4, 'preSmooth': 3, 'postSmooth': 3, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True,
    }},

    'omega-pipe-p1p1-1-1-2-3fwd' : {**omega_base, **{
        'discretization': 'p1p1', 'maxLevel': 4, 'preSmooth': 1, 'postSmooth': 1, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 3, 'symmGSVelocity': False,
    }},

    'omega-pipe-p1p1-2-2-2-3fwd' : {**omega_base, **{
        'discretization': 'p1p1', 'maxLevel': 4, 'preSmooth': 2, 'postSmooth': 2, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 3, 'symmGSVelocity': False,
    }},

    'omega-pipe-p1p1-3-3-2-3fwd' : {**omega_base, **{
        'discretization': 'p1p1', 'maxLevel': 4, 'preSmooth': 3, 'postSmooth': 3, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 3, 'symmGSVelocity': False,
    }},

    'omega-pipe-p2p1-1-1-2-1symm' : {**omega_base, **{
        'discretization': 'p2p1', 'maxLevel': 3, 'preSmooth': 1, 'postSmooth': 1, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True,
    }},

    'omega-pipe-p2p1-2-2-2-1symm' : {**omega_base, **{
        'discretization': 'p2p1', 'maxLevel': 3, 'preSmooth': 2, 'postSmooth': 2, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True,
    }},

    'omega-pipe-p2p1-3-3-2-1symm' : {**omega_base, **{
        'discretization': 'p2p1', 'maxLevel': 3, 'preSmooth': 3, 'postSmooth': 3, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': True,
    }},

    'omega-pipe-p2p1-1-1-2-3fwd' : {**omega_base, **{
        'discretization': 'p2p1', 'maxLevel': 3, 'preSmooth': 1, 'postSmooth': 1, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 3, 'symmGSVelocity': False,
    }},

    'omega-pipe-p2p1-2-2-2-3fwd' : {**omega_base, **{
        'discretization': 'p2p1', 'maxLevel': 3, 'preSmooth': 2, 'postSmooth': 2, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 3, 'symmGSVelocity': False,
    }},

    'omega-pipe-p2p1-3-3-2-3fwd' : {**omega_base, **{
        'discretization': 'p2p1', 'maxLevel': 3, 'preSmooth': 3, 'postSmooth': 3, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 3, 'symmGSVelocity': False,
    }},

    ########################### discr error ###########################

    'discr-cube-p1p1': {**discr_base, **{'discretization': 'p1p1', 'omega': 0.1}},
    'discr-cube-p2p1': {**discr_base, **{'discretization': 'p2p1', 'omega': 0.2}},

    ########################

    'fmg-cube-p1p1': {'discretization': 'p1p1', 'minLevel': 0, 'maxLevel': 6, 'numEdgesPerSide': 1, 'vtk': False, 'numCycles': 1, 'absoluteResidualTolerance': 1e-12,
                      'estimateOmega': False, 'omega': 0.1, 'coarseGridSolverType': 1, 'normCalculationLevelIncrement': 1, 'fmgInnerIterations': 0},

    'fmg-cube-p2p1': {'discretization': 'p2p1', 'minLevel': 0, 'maxLevel': 5, 'numEdgesPerSide': 1, 'vtk': False, 'numCycles': 1, 'absoluteResidualTolerance': 1e-12,
                      'estimateOmega': False, 'omega': 0.2, 'coarseGridSolverType': 1, 'normCalculationLevelIncrement': 1, 'fmgInnerIterations': 0},

    'fmg-pipe-p1p1': {'discretization': 'p1p1', 'minLevel': 0, 'maxLevel': 4, 'numEdgesPerSide': 1, 'vtk': False, 'numCycles': 1, 'absoluteResidualTolerance': 1e-12,
                      'estimateOmega': False, 'omega': 0.3, 'coarseGridSolverType': 1, 'normCalculationLevelIncrement': 1, 'fmgInnerIterations': 0},

    'fmg-pipe-p2p1': {'discretization': 'p2p1', 'minLevel': 0, 'maxLevel': 3, 'numEdgesPerSide': 1, 'vtk': False, 'numCycles': 1, 'absoluteResidualTolerance': 1e-12,
                      'estimateOmega': False, 'omega': 0.4, 'coarseGridSolverType': 1, 'normCalculationLevelIncrement': 1, 'fmgInnerIterations': 0},
}

parameterizations = {
    'discr-cube-p1p1': [
        {'maxLevel': 0},
        {'maxLevel': 1},
        {'maxLevel': 2},
        {'maxLevel': 3},
        {'maxLevel': 4},
        {'maxLevel': 5},
        {'maxLevel': 6},
    ],
    'discr-cube-p2p1': [
        {'maxLevel': 0},
        {'maxLevel': 1},
        {'maxLevel': 2},
        {'maxLevel': 3},
        {'maxLevel': 4},
        {'maxLevel': 5},
    ],
    'discr-pipe-p1p1': [
        {'maxLevel': 0},
        {'maxLevel': 1},
        {'maxLevel': 2},
        {'maxLevel': 3},
        {'maxLevel': 4},
        {'maxLevel': 5},
    ],
    'discr-pipe-p2p1': [
        {'maxLevel': 0},
        {'maxLevel': 1},
        {'maxLevel': 2},
        {'maxLevel': 3},
        {'maxLevel': 4},
    ],

    'fmg-cube-p1p1' : [],
    'fmg-cube-p2p1' : [],

    'fmg-pipe-p1p1': [{'preSmooth': 1, 'postSmooth': 0, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 2, 'symmGSVelocity': False},
                      {'preSmooth': 1, 'postSmooth': 2, 'incSmooth': 2, 'fmgInnerIterations': 1, 'numGSVelocity': 2, 'symmGSVelocity': False},
                      {'preSmooth': 1, 'postSmooth': 1, 'incSmooth': 1, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': False, 'solveWithCoarseGridSolverOnEachFMGLevel': True}],

    'fmg-pipe-p2p1': [{'preSmooth': 1, 'postSmooth': 2, 'incSmooth': 3, 'fmgInnerIterations': 1, 'numGSVelocity': 2, 'symmGSVelocity': False},
                      {'preSmooth': 2, 'postSmooth': 3, 'incSmooth': 3, 'fmgInnerIterations': 1, 'numGSVelocity': 2, 'symmGSVelocity': False},
                      {'preSmooth': 1, 'postSmooth': 1, 'incSmooth': 1, 'fmgInnerIterations': 1, 'numGSVelocity': 1, 'symmGSVelocity': False, 'solveWithCoarseGridSolverOnEachFMGLevel': True}],
}

for prm in ['omega-cube-p1p1-1-1-2-1symm', 'omega-cube-p1p1-2-2-2-1symm', 'omega-cube-p1p1-3-3-2-1symm', 'omega-cube-p1p1-1-1-2-3fwd', 'omega-cube-p1p1-2-2-2-3fwd', 'omega-cube-p1p1-3-3-2-3fwd',
            'omega-cube-p2p1-1-1-2-1symm', 'omega-cube-p2p1-2-2-2-1symm', 'omega-cube-p2p1-3-3-2-1symm', 'omega-cube-p2p1-1-1-2-3fwd', 'omega-cube-p2p1-2-2-2-3fwd', 'omega-cube-p2p1-3-3-2-3fwd']:
    parameterizations[prm] = []

for prm in ['omega-pipe-p1p1-1-1-2-1symm', 'omega-pipe-p1p1-2-2-2-1symm', 'omega-pipe-p1p1-3-3-2-1symm', 'omega-pipe-p1p1-1-1-2-3fwd', 'omega-pipe-p1p1-2-2-2-3fwd', 'omega-pipe-p1p1-3-3-2-3fwd',
            'omega-pipe-p2p1-1-1-2-1symm', 'omega-pipe-p2p1-2-2-2-1symm', 'omega-pipe-p2p1-3-3-2-1symm', 'omega-pipe-p2p1-1-1-2-3fwd', 'omega-pipe-p2p1-2-2-2-3fwd', 'omega-pipe-p2p1-3-3-2-3fwd']:
    parameterizations[prm] = []

for omega in np.linspace(0.1, 1.0, 10):

    for prm in ['omega-cube-p1p1-1-1-2-1symm', 'omega-cube-p1p1-2-2-2-1symm', 'omega-cube-p1p1-3-3-2-1symm', 'omega-cube-p1p1-1-1-2-3fwd', 'omega-cube-p1p1-2-2-2-3fwd', 'omega-cube-p1p1-3-3-2-3fwd',
                'omega-cube-p2p1-1-1-2-1symm', 'omega-cube-p2p1-2-2-2-1symm', 'omega-cube-p2p1-3-3-2-1symm', 'omega-cube-p2p1-1-1-2-3fwd', 'omega-cube-p2p1-2-2-2-3fwd', 'omega-cube-p2p1-3-3-2-3fwd']:
        parameterizations[prm].append({'omega': omega})

    for prm in ['omega-pipe-p1p1-1-1-2-1symm', 'omega-pipe-p1p1-2-2-2-1symm', 'omega-pipe-p1p1-3-3-2-1symm', 'omega-pipe-p1p1-1-1-2-3fwd', 'omega-pipe-p1p1-2-2-2-3fwd', 'omega-pipe-p1p1-3-3-2-3fwd',
                'omega-pipe-p2p1-1-1-2-1symm', 'omega-pipe-p2p1-2-2-2-1symm', 'omega-pipe-p2p1-3-3-2-1symm', 'omega-pipe-p2p1-1-1-2-3fwd', 'omega-pipe-p2p1-2-2-2-3fwd', 'omega-pipe-p2p1-3-3-2-3fwd']:
        parameterizations[prm].append({'omega': omega})

for pre in range(4):
    for post in range(4):
        for inc in range(4):
            for kappa in range(1, 3):
                for gs, symm in [(1, True), (2, True), (1, False), (2, False), (3, False), (4, False)]:
                    parameterizations['fmg-cube-p1p1'].append({'preSmooth': pre, 'postSmooth': post, 'incSmooth': inc, 'fmgInnerIterations': kappa, 'numGSVelocity': gs, 'symmGSVelocity': symm})
                    parameterizations['fmg-cube-p2p1'].append({'preSmooth': pre, 'postSmooth': post, 'incSmooth': inc, 'fmgInnerIterations': kappa, 'numGSVelocity': gs, 'symmGSVelocity': symm})


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


# for prm in ['omega-cube-p1p1-1-1-2-1symm', 'omega-cube-p1p1-2-2-2-1symm', 'omega-cube-p1p1-3-3-2-1symm', 'omega-cube-p1p1-1-1-2-3fwd', 'omega-cube-p1p1-2-2-2-3fwd', 'omega-cube-p1p1-3-3-2-3fwd',
#             'omega-cube-p2p1-1-1-2-1symm', 'omega-cube-p2p1-2-2-2-1symm', 'omega-cube-p2p1-3-3-2-1symm', 'omega-cube-p2p1-1-1-2-3fwd', 'omega-cube-p2p1-2-2-2-3fwd', 'omega-cube-p2p1-3-3-2-3fwd']:
#     run_all_configs(prm)

# for prm in ['omega-pipe-p1p1-1-1-2-1symm', 'omega-pipe-p1p1-2-2-2-1symm', 'omega-pipe-p1p1-3-3-2-1symm', 'omega-pipe-p1p1-1-1-2-3fwd', 'omega-pipe-p1p1-2-2-2-3fwd', 'omega-pipe-p1p1-3-3-2-3fwd',
#             'omega-pipe-p2p1-1-1-2-1symm', 'omega-pipe-p2p1-2-2-2-1symm', 'omega-pipe-p2p1-3-3-2-1symm', 'omega-pipe-p2p1-1-1-2-3fwd', 'omega-pipe-p2p1-2-2-2-3fwd', 'omega-pipe-p2p1-3-3-2-3fwd']:
#     run_all_configs(prm)

# for prm in ['discr-cube-p1p1', 'discr-cube-p2p1']:
#     run_all_configs(prm)

# for prm in ['fmg-cube-p1p1', 'fmg-cube-p2p1']:
#     run_all_configs(prm)

for prm in ['fmg-pipe-p1p1', 'fmg-pipe-p2p1']:
    run_all_configs(prm)
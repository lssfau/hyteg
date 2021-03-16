
import datetime
import time
import os
import subprocess

executables = {
    1: 'Benchmark_01_CircularAdvection',
    2: 'Benchmark_02_SwirlingAdvection',
    4: 'Benchmark_04_BlendedAdvectionDiffusion',
    5: 'Benchmark_01_CircularAdvection',
}

parameterizations = {
    1: [
        # comparison with FCT
        {'threeDim': False, 'spaceDiscretization': 'P1', 'level': 7, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': False, 'spaceDiscretization': 'P1', 'level': 7, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 10, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': False, 'spaceDiscretization': 'P1', 'level': 7, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 100, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': False, 'spaceDiscretization': 'P1', 'level': 7, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 1000, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': False, 'spaceDiscretization': 'P1', 'level': 7, 'numTimeSteps': 6283, 'resetParticles': False, 'resetParticlesInterval': 10000, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': False, 'spaceDiscretization': 'P1', 'level': 7, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': True, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},

        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 10, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 100, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 1000, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 6283, 'resetParticles': False, 'resetParticlesInterval': 10000, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': True, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},

        # time step convergence (basically RK-4 benchmark)
        {'threeDim': False, 'spaceDiscretization': 'P1', 'level': 7, 'numTimeSteps': 3141, 'resetParticles': False, 'resetParticlesInterval': 10000, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': False, 'spaceDiscretization': 'P1', 'level': 7, 'numTimeSteps': 1570, 'resetParticles': False, 'resetParticlesInterval': 10000, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': False, 'spaceDiscretization': 'P1', 'level': 7, 'numTimeSteps': 785, 'resetParticles': False, 'resetParticlesInterval': 10000, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': False, 'spaceDiscretization': 'P1', 'level': 7, 'numTimeSteps': 392, 'resetParticles': False, 'resetParticlesInterval': 10000, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': False, 'spaceDiscretization': 'P1', 'level': 7, 'numTimeSteps': 196, 'resetParticles': False, 'resetParticlesInterval': 10000, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': False, 'spaceDiscretization': 'P1', 'level': 7, 'numTimeSteps': 98, 'resetParticles': False, 'resetParticlesInterval': 10000, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
    ],

    5: [
        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False, 'enableGaussianCone': True, 'enableLinearCone': False, 'enableCylinder': False},
        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 10, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False, 'enableGaussianCone': True, 'enableLinearCone': False, 'enableCylinder': False},
        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 62, 'resetParticles': False, 'resetParticlesInterval': 100, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False, 'enableGaussianCone': True, 'enableLinearCone': False, 'enableCylinder': False},

        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 628, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False, 'enableGaussianCone': True, 'enableLinearCone': False, 'enableCylinder': False},
        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 628, 'resetParticles': True, 'resetParticlesInterval': 10, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False, 'enableGaussianCone': True, 'enableLinearCone': False, 'enableCylinder': False},
        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 628, 'resetParticles': True, 'resetParticlesInterval': 100, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False, 'enableGaussianCone': True, 'enableLinearCone': False, 'enableCylinder': False},
        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 628, 'resetParticles': False, 'resetParticlesInterval': 1000, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False, 'enableGaussianCone': True, 'enableLinearCone': False, 'enableCylinder': False},

        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False, 'enableGaussianCone': True, 'enableLinearCone': False, 'enableCylinder': False},
        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 10, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False, 'enableGaussianCone': True, 'enableLinearCone': False, 'enableCylinder': False},
        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 100, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False, 'enableGaussianCone': True, 'enableLinearCone': False, 'enableCylinder': False},
        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 1000, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False, 'enableGaussianCone': True, 'enableLinearCone': False, 'enableCylinder': False},
        {'threeDim': False, 'spaceDiscretization': 'P2', 'level': 6, 'numTimeSteps': 6283, 'resetParticles': False, 'resetParticlesInterval': 10000, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False, 'enableGaussianCone': True, 'enableLinearCone': False, 'enableCylinder': False},
    ],

    2: [
        {'threeDim': True, 'spaceDiscretization': 'P1', 'level': 5, 'numTimeSteps': 15, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': True, 'spaceDiscretization': 'P1', 'level': 6, 'numTimeSteps': 15, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': True, 'spaceDiscretization': 'P1', 'level': 7, 'numTimeSteps': 15, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},

        {'threeDim': True, 'spaceDiscretization': 'P1', 'level': 5, 'numTimeSteps': 30, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': True, 'spaceDiscretization': 'P1', 'level': 6, 'numTimeSteps': 30, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': True, 'spaceDiscretization': 'P1', 'level': 7, 'numTimeSteps': 30, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},

        {'threeDim': True, 'spaceDiscretization': 'P1', 'level': 5, 'numTimeSteps': 60, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': True, 'spaceDiscretization': 'P1', 'level': 6, 'numTimeSteps': 60, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
        {'threeDim': True, 'spaceDiscretization': 'P1', 'level': 7, 'numTimeSteps': 60, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False, 'globalMaxLimiter': False},
    ],

    4: [
        {'threeDim': False, 'level': 4, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'diffusivity': 1e-3, 'printInterval': 1, 'vtk': False, 'strangSplitting': False},
        {'threeDim': False, 'level': 5, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'diffusivity': 1e-3, 'printInterval': 1, 'vtk': False, 'strangSplitting': False},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'diffusivity': 1e-3, 'printInterval': 1, 'vtk': False, 'strangSplitting': False},

        {'threeDim': False, 'level': 4, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'diffusivity': 1e-5, 'printInterval': 1, 'vtk': False, 'strangSplitting': False},
        {'threeDim': False, 'level': 5, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'diffusivity': 1e-5, 'printInterval': 1, 'vtk': False, 'strangSplitting': False},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'diffusivity': 1e-5, 'printInterval': 1, 'vtk': False, 'strangSplitting': False},

        {'threeDim': False, 'level': 4, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'diffusivity': 1e-7, 'printInterval': 1, 'vtk': False, 'strangSplitting': False},
        {'threeDim': False, 'level': 5, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'diffusivity': 1e-7, 'printInterval': 1, 'vtk': False, 'strangSplitting': False},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'diffusivity': 1e-7, 'printInterval': 1, 'vtk': False, 'strangSplitting': False},
    ]
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
        for prm_key, prm_val in parameterization.items():
            cmd += '-Parameters.{prm_key}={prm_val} '.format(prm_key=prm_key, prm_val=prm_val)
        print(cmd)
        completed_process = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE, universal_newlines=True)
        out_file = os.path.join(db_dir, '_'.join([id, str(run_id)]) + '.out')
        with open(out_file, 'w') as f:
            f.write(completed_process.stdout)

#run_all_configs(1)
#run_all_configs(2)
#run_all_configs(4)
run_all_configs(5)

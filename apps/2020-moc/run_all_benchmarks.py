
import datetime
import time
import os
import subprocess

executables = {
    1: 'Benchmark_01_CircularAdvection',
    2: 'Benchmark_02_SwirlingAdvection',
    4: 'Benchmark_04_BlendedAdvectionDiffusion',
}

parameterizations = {
    1: [
        # h and dt convergence measurements, adj. adv. tests
        {'threeDim': False, 'level': 5, 'numTimeSteps': 62, 'resetParticles': False, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 62, 'resetParticles': False, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},

        {'threeDim': False, 'level': 5, 'numTimeSteps': 628, 'resetParticles': False, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 628, 'resetParticles': False, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},

        {'threeDim': False, 'level': 5, 'numTimeSteps': 6283, 'resetParticles': False, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 6283, 'resetParticles': False, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},

        {'threeDim': False, 'level': 5, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},

        {'threeDim': False, 'level': 5, 'numTimeSteps': 628, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 628, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},

        {'threeDim': False, 'level': 5, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},

        {'threeDim': False, 'level': 5, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': True, 'printInterval': 1, 'vtk': False},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': True, 'printInterval': 1, 'vtk': False},

        {'threeDim': False, 'level': 5, 'numTimeSteps': 628, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': True, 'printInterval': 1, 'vtk': False},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 628, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': True, 'printInterval': 1, 'vtk': False},

        {'threeDim': False, 'level': 5, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': True, 'printInterval': 1, 'vtk': False},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 6283, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': True, 'printInterval': 1, 'vtk': False},

        # particle reset interval test, these maybe need to be repeated for VTK
        {'threeDim': False, 'level': 5, 'numTimeSteps': 628, 'resetParticles': True, 'resetParticlesInterval': 10, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': True},
        {'threeDim': False, 'level': 5, 'numTimeSteps': 628, 'resetParticles': True, 'resetParticlesInterval': 100, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': True},
    ],

    2: [
        # maybe do a time-step convergence study for this method
        {'threeDim': True, 'level': 4, 'numTimeSteps': 15, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},
        {'threeDim': True, 'level': 5, 'numTimeSteps': 15, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},

        {'threeDim': True, 'level': 4, 'numTimeSteps': 15, 'resetParticles': True, 'adjustedAdvection': True, 'printInterval': 1, 'vtk': False},
        {'threeDim': True, 'level': 5, 'numTimeSteps': 15, 'resetParticles': True, 'adjustedAdvection': True, 'printInterval': 1, 'vtk': False},

        {'threeDim': True, 'level': 4, 'numTimeSteps': 30, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},
        {'threeDim': True, 'level': 5, 'numTimeSteps': 30, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},

        {'threeDim': True, 'level': 4, 'numTimeSteps': 30, 'resetParticles': True, 'adjustedAdvection': True, 'printInterval': 1, 'vtk': False},
        {'threeDim': True, 'level': 5, 'numTimeSteps': 30, 'resetParticles': True, 'adjustedAdvection': True, 'printInterval': 1, 'vtk': False},

        {'threeDim': True, 'level': 4, 'numTimeSteps': 60, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},
        {'threeDim': True, 'level': 5, 'numTimeSteps': 60, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1, 'vtk': False},

        {'threeDim': True, 'level': 4, 'numTimeSteps': 60, 'resetParticles': True, 'adjustedAdvection': True, 'printInterval': 1, 'vtk': False},
        {'threeDim': True, 'level': 5, 'numTimeSteps': 60, 'resetParticles': True, 'adjustedAdvection': True, 'printInterval': 1, 'vtk': False},
    ],

    4: [
        {'threeDim': False, 'level': 5, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': True, 'diffusivity': 1e-3, 'printInterval': 1, 'vtk': False},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': True, 'diffusivity': 1e-3, 'printInterval': 1, 'vtk': False},
        {'threeDim': False, 'level': 5, 'numTimeSteps': 628, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': True, 'diffusivity': 1e-3, 'printInterval': 1, 'vtk': False},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 628, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': True, 'diffusivity': 1e-3, 'printInterval': 1, 'vtk': False},

        {'threeDim': False, 'level': 5, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': True, 'diffusivity': 1e-7, 'printInterval': 1, 'vtk': False},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 62, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': True, 'diffusivity': 1e-7, 'printInterval': 1, 'vtk': False},
        {'threeDim': False, 'level': 5, 'numTimeSteps': 628, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': True, 'diffusivity': 1e-7, 'printInterval': 1, 'vtk': False},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 628, 'resetParticles': True, 'resetParticlesInterval': 1, 'adjustedAdvection': True, 'diffusivity': 1e-7, 'printInterval': 1, 'vtk': False},
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
run_all_configs(4)
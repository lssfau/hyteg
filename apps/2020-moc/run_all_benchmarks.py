
import datetime
import time
import os
import subprocess

executables = {
    1: 'Benchmark_01_CircularAdvection',
}

parameterizations = {
    1: [
        {'threeDim': False, 'level': 5, 'numTimeSteps': 63, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 63, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1},
        {'threeDim': False, 'level': 7, 'numTimeSteps': 63, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1},

        {'threeDim': False, 'level': 5, 'numTimeSteps': 628, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 628, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1},
        {'threeDim': False, 'level': 7, 'numTimeSteps': 628, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1},

        {'threeDim': False, 'level': 5, 'numTimeSteps': 6283, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 6283, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1},
        {'threeDim': False, 'level': 7, 'numTimeSteps': 6283, 'resetParticles': False, 'adjustedAdvection': False, 'printInterval': 1},

        {'threeDim': False, 'level': 5, 'numTimeSteps': 63, 'resetParticles': True, 'adjustedAdvection': False, 'printInterval': 1},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 63, 'resetParticles': True, 'adjustedAdvection': False, 'printInterval': 1},
        {'threeDim': False, 'level': 7, 'numTimeSteps': 63, 'resetParticles': True, 'adjustedAdvection': False, 'printInterval': 1},

        {'threeDim': False, 'level': 5, 'numTimeSteps': 628, 'resetParticles': True, 'adjustedAdvection': False, 'printInterval': 1},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 628, 'resetParticles': True, 'adjustedAdvection': False, 'printInterval': 1},
        {'threeDim': False, 'level': 7, 'numTimeSteps': 628, 'resetParticles': True, 'adjustedAdvection': False, 'printInterval': 1},

        {'threeDim': False, 'level': 5, 'numTimeSteps': 6283, 'resetParticles': True, 'adjustedAdvection': False, 'printInterval': 1},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 6283, 'resetParticles': True, 'adjustedAdvection': False, 'printInterval': 1},
        {'threeDim': False, 'level': 7, 'numTimeSteps': 6283, 'resetParticles': True, 'adjustedAdvection': False, 'printInterval': 1},

        {'threeDim': False, 'level': 5, 'numTimeSteps': 63, 'resetParticles': True, 'adjustedAdvection': True, 'printInterval': 1},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 63, 'resetParticles': True, 'adjustedAdvection': True, 'printInterval': 1},
        {'threeDim': False, 'level': 7, 'numTimeSteps': 63, 'resetParticles': True, 'adjustedAdvection': True, 'printInterval': 1},

        {'threeDim': False, 'level': 5, 'numTimeSteps': 628, 'resetParticles': True, 'adjustedAdvection': True, 'printInterval': 1},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 628, 'resetParticles': True, 'adjustedAdvection': True, 'printInterval': 1},
        {'threeDim': False, 'level': 7, 'numTimeSteps': 628, 'resetParticles': True, 'adjustedAdvection': True, 'printInterval': 1},

        {'threeDim': False, 'level': 5, 'numTimeSteps': 6283, 'resetParticles': True, 'adjustedAdvection': True, 'printInterval': 1},
        {'threeDim': False, 'level': 6, 'numTimeSteps': 6283, 'resetParticles': True, 'adjustedAdvection': True, 'printInterval': 1},
        {'threeDim': False, 'level': 7, 'numTimeSteps': 6283, 'resetParticles': True, 'adjustedAdvection': True, 'printInterval': 1},
    ]
}

levels = {
    1: [3, 4, 5, 6, 7],
}

num_time_steps = {
    1: [62, 628, 6283, 62831, 628318],
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

run_all_configs(1)
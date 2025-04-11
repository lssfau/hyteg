import os
import sys
import subprocess

level = 4
nodes_procpernode = [(2, 120), (8, 120), (16, 120), (32, 120), (64, 120), (128, 120), (256, 120), (512, 120), (1024, 120), (2048, 120)]
ntan_nrad_level = [(3, 2, level), (5, 2, level), (5, 3, level), (5, 5, level), (9, 3, level), (9, 5, level), (9, 9, level), (17, 5, level), (17, 9, level), (17, 17, level)]

if len(sys.argv) < 2:
    raise ValueError("You need to pass an output directory")

outputDir = sys.argv[1]

outputDirLevel = f"{outputDir}/level_{level}"
prmDirLevel = f"{outputDir}/level_{level}/prmdir"
runnersDirLevel = f"{outputDir}/level_{level}/runnersdir"
pbsOutDirLevel = f"{outputDir}/level_{level}/pbsoutdir"
timingTreeDirLevel = f"{outputDir}/level_{level}/timingtreedir"

if f'level_{level}' not in os.listdir(outputDir):
    subprocess.call(["mkdir", "-p", outputDirLevel])

if 'prmdir' not in os.listdir(outputDirLevel):
    subprocess.call(["mkdir", "-p", prmDirLevel])

if 'runnersdir' not in os.listdir(outputDirLevel):
    subprocess.call(["mkdir", "-p", runnersDirLevel])

if 'pbsoutdir' not in os.listdir(outputDirLevel):
    subprocess.call(["mkdir", "-p", pbsOutDirLevel])

if 'timingtreedir' not in os.listdir(outputDirLevel):
    subprocess.call(["mkdir", "-p", timingTreeDirLevel])

for i in range(10):
    prmfile = f'{prmDirLevel}/SphericalShellBenchRotationMinimalMinimal_{i+1}{i+1}_level{level}.prm'
    runnerfile = f'{runnersDirLevel}/runner_{i+1}{i+1}_level{level}.pbs'
    timingtreefile = f'{timingTreeDirLevel}/timingTree_{i+1}{i+1}_level{level}.json'
    timingtreefile = timingtreefile.replace('/', '\/')
    logfile = f'{pbsOutDirLevel}/logfile_{i+1}{i+1}_level_{level}.out'
    logfile = logfile.replace('/', '\/') 
    pbsOutDirLevel = pbsOutDirLevel.replace('/', '\/')
    
    subprocess.call(["cp", "SphericalShellBenchRotationMinimalMinimal.prm", prmfile])
    subprocess.call(["sed", "-i",  f's/nTan .*;/nTan {ntan_nrad_level[i][0]};/', prmfile])
    subprocess.call(["sed", "-i",  f's/nRad .*;/nRad {ntan_nrad_level[i][1]};/', prmfile])
    subprocess.call(["sed", "-i",  f's/maxLevel [0-9];/maxLevel {ntan_nrad_level[i][2]};/', prmfile])
    subprocess.call(["sed", "-i",  's/minLevel [0-9];/minLevel 0;/', prmfile])
    subprocess.call(["sed", "-i",  f's/logFilename .*;/logFilename {logfile};/', prmfile])
    subprocess.call(["sed", "-i",  f's/timingTreeFilename .*;/timingTreeFilename {timingtreefile};/', prmfile])

    subprocess.call(["cp", "runner-scaling.pbs", runnerfile])
    # subprocess.call(["sed", "-i", f's/PBS -o .*/PBS -o {pbsOutDirLevel}/', runnerfile])
    # subprocess.call(["sed", "-i", f's/PBS -e .*/PBS -e {pbsOutDirLevel}/', runnerfile])
    subprocess.call(["sed", "-i", f's/PBS -N sph/PBS -N ttAvSwNWP1_{i+1}{i+1}L{level}/', runnerfile])
    subprocess.call(["sed", "-i", f's/select=.*:node_type=rome/select={nodes_procpernode[i][0]}:node_type=rome/', runnerfile])
    subprocess.call(["sed", "-i", f's/mpiprocs=.*:/mpiprocs={nodes_procpernode[i][1]}:/', runnerfile])
    subprocess.call(["sed", "-i", f's/mpirun -np \([^\/]*\)/mpirun -np {nodes_procpernode[i][0]*nodes_procpernode[i][1]} /', runnerfile])
    subprocess.call(["sed", "-i", 's/SphericalShellBenchRotationMinimalMinimal.prm/' + prmfile.replace('/', '\/') + '/', runnerfile])

    pbsOutDirLevel = pbsOutDirLevel.replace('\/', '/')
    # subprocess.call(f"qsub {runnerfile}", shell=True, cwd=pbsOutDirLevel)


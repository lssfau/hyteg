#!/bin/bash
#SBATCH --output=/SCRATCH/pponkumar/analyticalBenchmarks/minimal/slurm-out/slurm-%j.out
#SBATCH --chdir=/SCRATCH/pponkumar/
#SBATCH --job-name="3dsfall"
#SBATCH --ntasks=128
#SBATCH --mail-type=all
#SBATCH --time=71:59:59

# run compute job

for i in 2 3 4
do

sed -i "s/minLevel.*;/minLevel ${i};/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm
sed -i "s/maxLevel.*;/maxLevel ${i};/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm

# Zero slip;
# Smooth;

# sed -i "s/freeslip.*;/freeslip false;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm
# sed -i "s/mixed.*;/mixed false;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm
# sed -i "s/delta.*;/delta false;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm

# mpirun.openmpi /home/pponkumar/hyteg/hyteg-build/apps/2024-convbench/3DStokesFSConvTest

# Mixed zero and free slip; 
# Smooth;

sed -i "s/freeslip.*;/freeslip false;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm
sed -i "s/mixed.*;/mixed true;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm
sed -i "s/delta.*;/delta false;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm

mpirun.openmpi /home/pponkumar/hyteg/hyteg-build/apps/2024-convbench/3DStokesFSConvTest

# Free slip; 
# Smooth;

sed -i "s/freeslip.*;/freeslip true;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm
sed -i "s/mixed.*;/mixed false;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm
sed -i "s/delta.*;/delta false;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm

mpirun.openmpi /home/pponkumar/hyteg/hyteg-build/apps/2024-convbench/3DStokesFSConvTest

################# DELTA #################

# Zero slip;
# Smooth;

# sed -i "s/freeslip.*;/freeslip false;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm
# sed -i "s/mixed.*;/mixed false;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm
# sed -i "s/delta.*;/delta true;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm

# mpirun.openmpi /home/pponkumar/hyteg/hyteg-build/apps/2024-convbench/3DStokesFSConvTest

# Mixed zero and free slip; 
# Smooth;

sed -i "s/freeslip.*;/freeslip false;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm
sed -i "s/mixed.*;/mixed true;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm
sed -i "s/delta.*;/delta true;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm

mpirun.openmpi /home/pponkumar/hyteg/hyteg-build/apps/2024-convbench/3DStokesFSConvTest

# Free slip; 
# Smooth;

sed -i "s/freeslip.*;/freeslip true;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm
sed -i "s/mixed.*;/mixed false;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm
sed -i "s/delta.*;/delta true;/" /home/pponkumar/hyteg/hyteg/apps/2024-convbench/3DAnnulus.prm

mpirun.openmpi /home/pponkumar/hyteg/hyteg-build/apps/2024-convbench/3DStokesFSConvTest

done

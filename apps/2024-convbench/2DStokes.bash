# run compute job

for i in 2 3 4 5
do

sed -i "s/level.*;/level ${i};/" /home/pponkumar/hyteg/hyteg/apps/annulusAnalytical/2DAnnulus.prm

# Mixed zero and free slip; 
# Smooth;

sed -i "s/freeslip.*;/freeslip false;/" /home/pponkumar/hyteg/hyteg/apps/annulusAnalytical/2DAnnulus.prm
sed -i "s/mixed.*;/mixed true;/" /home/pponkumar/hyteg/hyteg/apps/annulusAnalytical/2DAnnulus.prm
sed -i "s/delta.*;/delta false;/" /home/pponkumar/hyteg/hyteg/apps/annulusAnalytical/2DAnnulus.prm

# mpirun -np 56 /home/pponkumar/hyteg/hyteg-build/apps/annulusAnalytical/2DStokesFSConvTest

# Free slip; 
# Smooth;

sed -i "s/freeslip.*;/freeslip true;/" /home/pponkumar/hyteg/hyteg/apps/annulusAnalytical/2DAnnulus.prm
sed -i "s/mixed.*;/mixed false;/" /home/pponkumar/hyteg/hyteg/apps/annulusAnalytical/2DAnnulus.prm
sed -i "s/delta.*;/delta false;/" /home/pponkumar/hyteg/hyteg/apps/annulusAnalytical/2DAnnulus.prm

mpirun -np 56 /home/pponkumar/hyteg/hyteg-build/apps/annulusAnalytical/2DStokesFSConvTest

################# DELTA #################

# Mixed zero and free slip; 
# Smooth;

sed -i "s/freeslip.*;/freeslip false;/" /home/pponkumar/hyteg/hyteg/apps/annulusAnalytical/2DAnnulus.prm
sed -i "s/mixed.*;/mixed true;/" /home/pponkumar/hyteg/hyteg/apps/annulusAnalytical/2DAnnulus.prm
sed -i "s/delta.*;/delta true;/" /home/pponkumar/hyteg/hyteg/apps/annulusAnalytical/2DAnnulus.prm

# mpirun -np 56 /home/pponkumar/hyteg/hyteg-build/apps/annulusAnalytical/2DStokesFSConvTest

# Free slip; 
# Smooth;

sed -i "s/freeslip.*;/freeslip true;/" /home/pponkumar/hyteg/hyteg/apps/annulusAnalytical/2DAnnulus.prm
sed -i "s/mixed.*;/mixed false;/" /home/pponkumar/hyteg/hyteg/apps/annulusAnalytical/2DAnnulus.prm
sed -i "s/delta.*;/delta true;/" /home/pponkumar/hyteg/hyteg/apps/annulusAnalytical/2DAnnulus.prm

# mpirun -np 56 /home/pponkumar/hyteg/hyteg-build/apps/annulusAnalytical/2DStokesFSConvTest

done

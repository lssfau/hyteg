#!/bin/bash

for i in 2 3 4 5
do

sed -i "s/maxLevel.*;/maxLevel ${i};/" 3DAnnulus.prm

sed -i "s/freeslip.*;/freeslip true;/" 3DAnnulus.prm
sed -i "s/delta.*;/delta true;/" 3DAnnulus.prm

mpiexec -np 128 ./3DStokesFlow

sed -i "s/freeslip.*;/freeslip true;/" 3DAnnulus.prm
sed -i "s/delta.*;/delta false;/" 3DAnnulus.prm

mpiexec -np 128 ./3DStokesFlow

sed -i "s/freeslip.*;/freeslip false;/" 3DAnnulus.prm
sed -i "s/delta.*;/delta true;/" 3DAnnulus.prm

mpiexec -np 128 ./3DStokesFlow

sed -i "s/freeslip.*;/freeslip false;/" 3DAnnulus.prm
sed -i "s/delta.*;/delta false;/" 3DAnnulus.prm

mpiexec -np 128 ./3DStokesFlow

done

# for i in 2 3 4
# do
#     sed -i "s/maxLevel.*;/maxLevel ${i};/" 3DAnnulus.prm
#     mpiexec -np 56 ./3DStokesFlow
# done

#!/bin/bash

Ra=( "7e5" "8e5" "9e5" "1e6" "125e4" "15e5" "175e4" "2e6" "225e4" "25e5" "275e5" "3e6" )
RaPrev=( "6e5" "7e5" "8e5" "9e5" "1e6" "125e4" "15e5" "175e4" "2e6" "225e4" "25e5" "275e5" )

for i in "${!Ra[@]}"; do
    printf "Ra = %s, RaPrev = %s\n" "${Ra[i]}" "${RaPrev[i]}"
    sed -i "s/RayleighNumber .*;/RayleighNumber ${Ra[i]};/" KingTALA.prm
    sed -i "s/outputFilename .*;/outputFilename KingSquare_ALA_Ra_${Ra[i]}_Di_0_25;/" KingTALA.prm
    sed -i "s/cpFilename .*;/cpFilename KingSquare_ALA_Ra_${Ra[i]}_Di_0_25_Continuous.bp;/" KingTALA.prm
    sed -i "s/cpStartFilename .*;/cpStartFilename KingSquare_ALA_Ra_${RaPrev[i]}_Di_0_25_Continuous.bp;/" KingTALA.prm

    mpirun -np 72 ./KingTALA | tee output/gmd_stuff/king/KingALA_Ra_${Ra[i]}_Di_0_25_P1.out
done

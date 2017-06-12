# TinyHHG C++

### Build instructions

    git clone git@i10git.cs.fau.de:walberla/walberla.git
    cd walberla
    git checkout master
    cd ..
    git clone git@i10git.cs.fau.de:terraneo/tinyhhg.git
    cd tinyhhg_cpp
    git checkout master
    cd ..
    mkdir build
    cd build
    cmake ../tinyhhg_cpp -DWALBERLA_DIR=../walberla
    cd apps
    make

### Execution instructions

    cd apps
    ./tinytest_fmg

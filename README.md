# TinyHHG C++

### Build instructions

    git clone git@i10git.cs.fau.de:software/walberla.git
    cd walberla
    git checkout master
    cd ..
    git clone git@i10git.cs.fau.de:drzisga/tinyhhg_cpp.git
    cd tinyhhg_cpp
    git checkout walberla_coupling
    cd ..
    mkdir build
    cd build
    cmake ../tinyhhg_cpp -DWALBERLA_DIR=../walberla
    cd apps
    make

### Execution instructions

    cd apps
    ./tinytest_fmg

#pragma once

#include "tinyhhg_core/types/flags.hpp"
#include <petscvec.h>


namespace hhg {

template<typename FunctionType>
class Vector {
public:

  Vec vec;

  Vector() = delete;

  Vector(uint_t size) {
    VecCreate(walberla::MPIManager::instance()->comm(), &vec);
    VecSetType(vec, VECSTANDARD);
    VecSetSizes(vec, PETSC_DECIDE, size);
    VecSetUp(vec);
  }

  ~Vector() { VecDestroy(&vec); }

  void createVectorFromFunction(FunctionType &src, FunctionType &numerator, uint_t level, DoFType flag = All) {

    src.createVector(numerator, vec, level, flag);

    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);
  }

  void print(const char name[])
  {
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD,name,&viewer);
    PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB );
    VecView(vec,viewer);
    PetscViewerDestroy(&viewer);
  }

};

}
#pragma once

#include "tinyhhg_core/types/flags.hpp"
#include <petscvec.h>


namespace hhg {

template<typename FunctionType>
class PETScVector {
public:

  Vec vec;

  PETScVector() = delete;

  PETScVector(uint_t size) {
    VecCreate(walberla::MPIManager::instance()->comm(), &vec);
    VecSetType(vec, VECSTANDARD);
    VecSetSizes(vec, PETSC_DECIDE, size);
    VecSetUp(vec);
  }

  ~PETScVector() { VecDestroy(&vec); }

  void createVectorFromFunction(FunctionType &src, FunctionType &numerator, uint_t level, DoFType flag = All) {

    src.createVectorFromFunction(numerator, vec, level, flag);

    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);
  }

  void createFunctionFromVector(FunctionType &src, FunctionType &numerator, uint_t level, DoFType flag = All){

    src.createFunctionFromVector(numerator, vec, level, flag);

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
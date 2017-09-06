#pragma once

#include "tinyhhg_core/types/flags.hpp"
#include "PETScWrapper.hpp"

#ifdef HHG_BUILD_WITH_PETSC

namespace hhg {

template<typename FunctionType>
class PETScVector {
protected:

  Vec vec;


public:
  PETScVector() = delete;

  PETScVector(uint_t size,const char name[] = "Vec") {
    VecCreate(walberla::MPIManager::instance()->comm(), &vec);
    VecSetType(vec, VECSTANDARD);
    VecSetSizes(vec, PETSC_DECIDE, size);
    VecSetUp(vec);
    setName(name);
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

  void print(const char filename[])
  {
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
    PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB );
    VecView(vec,viewer);
    PetscViewerDestroy(&viewer);
  }

  inline void setName(const char name[]){ PetscObjectSetName((PetscObject)vec,name); }

  inline Vec& get() { return vec; }


};

}

#endif
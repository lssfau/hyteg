#pragma once

#include "tinyhhg_core/types/flags.hpp"
#include "PETScWrapper.hpp"

#ifdef HHG_BUILD_WITH_PETSC

#include "tinyhhg_core/p1functionspace/P1Petsc.hpp"

namespace hhg {

template<typename ValueType, template <class> class FunctionType>
class PETScVector {
protected:

  Vec vec;


public:
  PETScVector() = delete;

  PETScVector(uint_t localSize, const std::string& name = "Vec") {
    VecCreate(walberla::MPIManager::instance()->comm(), &vec);
    VecSetType(vec, VECSTANDARD);
    VecSetSizes(vec, (PetscInt)localSize, PETSC_DECIDE);
    VecSetUp(vec);
    setName(name.c_str());
  }

  ~PETScVector() { VecDestroy(&vec); }

  void createVectorFromFunction(FunctionType<ValueType> &src, FunctionType<PetscInt> &numerator, uint_t level, DoFType flag = All) {

    hhg::petsc::createVectorFromFunction(src, numerator, vec, level, flag);

    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);
  }

  void createFunctionFromVector(FunctionType<ValueType> &src, FunctionType<PetscInt> &numerator, uint_t level, DoFType flag = All){

    hhg::petsc::createFunctionFromVector(src, numerator, vec, level, flag);

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
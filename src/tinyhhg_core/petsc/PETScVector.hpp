#pragma once

#include "tinyhhg_core/types/flags.hpp"
#include "tinyhhg_core/FunctionProperties.hpp"
#include "PETScWrapper.hpp"

#ifdef HHG_BUILD_WITH_PETSC

#include "tinyhhg_core/p1functionspace/P1Petsc.hpp"
#include "tinyhhg_core/composites/petsc/P1StokesPetsc.hpp"
#include "tinyhhg_core/p2functionspace/P2Petsc.hpp"
#include "tinyhhg_core/composites/petsc/P2P1TaylorHoodPetsc.hpp"
#include "tinyhhg_core/composites/petsc/P2P2StabilizedStokesPetsc.hpp"

namespace hyteg {

template<typename ValueType, template <class> class FunctionType>
class PETScVector {
protected:

  Vec vec;


public:
  PETScVector() = delete;

  PETScVector(const FunctionType< ValueType > & function, const FunctionType<PetscInt> & numerator, const uint_t & level, const DoFType & flag = All, const std::string& name = "Vec" )
    : PETScVector( numberOfLocalDoFs< typename FunctionType< ValueType >::Tag >( *function.getStorage(), level ), name )
  {
    createVectorFromFunction(function, numerator, level, flag);
  }

  PETScVector(uint_t localSize, const std::string& name = "Vec") {
    VecCreate(walberla::MPIManager::instance()->comm(), &vec);
    VecSetType(vec, VECSTANDARD);
    VecSetSizes(vec, (PetscInt)localSize, PETSC_DECIDE);
    VecSetUp(vec);
    setName(name.c_str());
  }

  ~PETScVector() { VecDestroy(&vec); }

  void createVectorFromFunction(const FunctionType<ValueType> &src,const FunctionType<PetscInt> &numerator, uint_t level, DoFType flag = All) {
     hyteg::petsc::createVectorFromFunction(src, numerator, vec, level, flag);

    VecAssemblyBegin(vec);
    VecAssemblyEnd(vec);
  }

  void createFunctionFromVector(const FunctionType<ValueType> &src,const FunctionType<PetscInt> &numerator, uint_t level, DoFType flag = All){
     hyteg::petsc::createFunctionFromVector(src, numerator, vec, level, flag);

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
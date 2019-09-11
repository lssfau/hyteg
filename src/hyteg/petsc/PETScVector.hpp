/*
 * Copyright (c) 2017-2019 Boerge Struempfel, Daniel Drzisga, Dominik Thoennes, Nils Kohl.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include "hyteg/types/flags.hpp"
#include "hyteg/FunctionProperties.hpp"
#include "PETScWrapper.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

#include "hyteg/p1functionspace/P1Petsc.hpp"
#include "hyteg/composites/petsc/P1StokesPetsc.hpp"
#include "hyteg/p2functionspace/P2Petsc.hpp"
#include "hyteg/composites/petsc/P2P1TaylorHoodPetsc.hpp"
#include "hyteg/composites/petsc/P2P2StabilizedStokesPetsc.hpp"

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
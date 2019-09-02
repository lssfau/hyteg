/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/p1functionspace/P1Petsc.hpp"

namespace hyteg {
namespace petsc {

inline void createVectorFromFunction(const P1StokesFunction<PetscScalar> &function,
                                     const P1StokesFunction<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  createVectorFromFunction(function.u, numerator.u, vec, level, flag);
  createVectorFromFunction(function.v, numerator.v, vec, level, flag);
  if ( function.u.getStorage()->hasGlobalCells() )
  {
    createVectorFromFunction(function.w, numerator.w, vec, level, flag);
  }
  createVectorFromFunction(function.p, numerator.p, vec, level, flag);
}

inline void createFunctionFromVector(const P1StokesFunction<PetscScalar> &function,
                                     const P1StokesFunction<PetscInt> &numerator,
                                     Vec &vec,
                                     uint_t level,
                                     DoFType flag) {
  createFunctionFromVector(function.u, numerator.u, vec, level, flag);
  createFunctionFromVector(function.v, numerator.v, vec, level, flag);
  if ( function.u.getStorage()->hasGlobalCells() )
  {
    createFunctionFromVector( function.w, numerator.w, vec, level, flag );
  }
  createFunctionFromVector(function.p, numerator.p, vec, level, flag);
}

inline void applyDirichletBC(const P1StokesFunction<PetscInt> &numerator, std::vector<PetscInt> &mat, uint_t level) {
  applyDirichletBC(numerator.u, mat, level);
  applyDirichletBC(numerator.v, mat, level);
  if ( numerator.u.getStorage()->hasGlobalCells() )
  {
    applyDirichletBC( numerator.w, mat, level );
  }
//  applyDirichletBC(numerator.p, mat, level);
}

template < class OperatorType >
inline void createMatrix( const OperatorType&                 opr,
                          const P1StokesFunction< PetscInt >& src,
                          const P1StokesFunction< PetscInt >& dst,
                          Mat&                                mat,
                          size_t                              level,
                          DoFType                             flag )
{
  createMatrix(opr.A, src.u, dst.u, mat, level, flag);
  createMatrix(opr.divT_x, src.p, dst.u, mat, level, flag);

  createMatrix(opr.A, src.v, dst.v, mat, level, flag);
  createMatrix(opr.divT_y, src.p, dst.v, mat, level, flag);

  if ( src.u.getStorage()->hasGlobalCells() )
  {
    createMatrix( opr.A, src.w, dst.w, mat, level, flag );
    createMatrix( opr.divT_z, src.p, dst.w, mat, level, flag );
  }

  createMatrix(opr.div_x, src.u, dst.p, mat, level, flag | DirichletBoundary);
  createMatrix(opr.div_y, src.v, dst.p, mat, level, flag | DirichletBoundary);
  if ( src.u.getStorage()->hasGlobalCells() )
  {
    createMatrix( opr.div_z, src.w, dst.p, mat, level, flag | DirichletBoundary );
  }
  createMatrix(opr.pspg, src.p, dst.p, mat, level, flag | DirichletBoundary);

}

template <>
inline void createMatrix< P1StokesBlockPreconditioner >( const P1StokesBlockPreconditioner&  opr,
                                                         const P1StokesFunction< PetscInt >& src,
                                                         const P1StokesFunction< PetscInt >& dst,
                                                         Mat&                                mat,
                                                         size_t                              level,
                                                         DoFType                             flag )
{
   createMatrix( opr.A, src.u, dst.u, mat, level, flag );
   createMatrix( opr.A, src.v, dst.v, mat, level, flag );

   if ( src.u.getStorage()->hasGlobalCells() )
   {
      createMatrix( opr.A, src.w, dst.w, mat, level, flag );
   }

   createMatrix( opr.P, src.p, dst.p, mat, level, flag );
}
}
}
/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesBlockPreconditioner.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesBlockPreconditioner.hpp"
#include "hyteg/p1functionspace/P1Petsc.hpp"
#include "hyteg/p2functionspace/P2Petsc.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

namespace hyteg {
namespace petsc {

inline void createVectorFromFunction( const P2P1TaylorHoodFunction< PetscReal >& function,
                                      const P2P1TaylorHoodFunction< idx_t >&     numerator,
                                      const std::shared_ptr< VectorProxy >&      vec,
                                      uint_t                                     level,
                                      DoFType                                    flag )
{
   createVectorFromFunction( function.uvw[0], numerator.uvw[0], vec, level, flag );
   createVectorFromFunction( function.uvw[1], numerator.uvw[1], vec, level, flag );
   if ( function.uvw[0].getStorage()->hasGlobalCells() )
   {
      createVectorFromFunction( function.uvw[2], numerator.uvw[2], vec, level, flag );
   }
   createVectorFromFunction( function.p, numerator.p, vec, level, flag );
}

inline void createFunctionFromVector( const P2P1TaylorHoodFunction< PetscReal >& function,
                                      const P2P1TaylorHoodFunction< idx_t >&     numerator,
                                      const std::shared_ptr< VectorProxy >&      vec,
                                      uint_t                                     level,
                                      DoFType                                    flag )
{
   createFunctionFromVector( function.uvw[0], numerator.uvw[0], vec, level, flag );
   createFunctionFromVector( function.uvw[1], numerator.uvw[1], vec, level, flag );
   if ( function.uvw[0].getStorage()->hasGlobalCells() )
   {
      createFunctionFromVector( function.uvw[2], numerator.uvw[2], vec, level, flag );
   }
   createFunctionFromVector( function.p, numerator.p, vec, level, flag );
}

inline void applyDirichletBC( const P2P1TaylorHoodFunction< idx_t >& numerator, std::vector< PetscInt >& mat, uint_t level )
{
   applyDirichletBC( numerator.uvw[0], mat, level );
   applyDirichletBC( numerator.uvw[1], mat, level );
   if ( numerator.uvw[0].getStorage()->hasGlobalCells() )
   {
      applyDirichletBC( numerator.uvw[2], mat, level );
   }
   //  applyDirichletBC(numerator.p, mat, level);
}

template < class OperatorType >
inline void createMatrix( const OperatorType&                         opr,
                          const P2P1TaylorHoodFunction< idx_t >&      src,
                          const P2P1TaylorHoodFunction< idx_t >&      dst,
                          const std::shared_ptr< SparseMatrixProxy >& mat,
                          size_t                                      level,
                          DoFType                                     flag )
{
   createMatrix( opr.A, src.uvw[0], dst.uvw[0], mat, level, flag );
   createMatrix( opr.divT_x.getVertexToVertexOpr(), src.p, dst.uvw[0].getVertexDoFFunction(), mat, level, flag );
   createMatrix( opr.divT_x.getVertexToEdgeOpr(), src.p, dst.uvw[0].getEdgeDoFFunction(), mat, level, flag );

   createMatrix( opr.A, src.uvw[1], dst.uvw[1], mat, level, flag );
   createMatrix( opr.divT_y.getVertexToVertexOpr(), src.p, dst.uvw[1].getVertexDoFFunction(), mat, level, flag );
   createMatrix( opr.divT_y.getVertexToEdgeOpr(), src.p, dst.uvw[1].getEdgeDoFFunction(), mat, level, flag );

   if ( src.uvw[0].getStorage()->hasGlobalCells() )
   {
      createMatrix( opr.A, src.uvw[2], dst.uvw[2], mat, level, flag );
      createMatrix( opr.divT_z.getVertexToVertexOpr(), src.p, dst.uvw[2].getVertexDoFFunction(), mat, level, flag );
      createMatrix( opr.divT_z.getVertexToEdgeOpr(), src.p, dst.uvw[2].getEdgeDoFFunction(), mat, level, flag );
   }

   createMatrix(
       opr.div_x.getVertexToVertexOpr(), src.uvw[0].getVertexDoFFunction(), dst.p, mat, level, flag | DirichletBoundary );
   createMatrix( opr.div_x.getEdgeToVertexOpr(), src.uvw[0].getEdgeDoFFunction(), dst.p, mat, level, flag | DirichletBoundary );
   createMatrix(
       opr.div_y.getVertexToVertexOpr(), src.uvw[1].getVertexDoFFunction(), dst.p, mat, level, flag | DirichletBoundary );
   createMatrix( opr.div_y.getEdgeToVertexOpr(), src.uvw[1].getEdgeDoFFunction(), dst.p, mat, level, flag | DirichletBoundary );
   if ( src.uvw[0].getStorage()->hasGlobalCells() )
   {
      createMatrix(
          opr.div_z.getVertexToVertexOpr(), src.uvw[2].getVertexDoFFunction(), dst.p, mat, level, flag | DirichletBoundary );
      createMatrix( opr.div_z.getEdgeToVertexOpr(), src.uvw[2].getEdgeDoFFunction(), dst.p, mat, level, flag | DirichletBoundary );
   }
}

template <>
inline void createMatrix< P2P1TaylorHoodStokesBlockPreconditioner >( const P2P1TaylorHoodStokesBlockPreconditioner& opr,
                                                                     const P2P1TaylorHoodFunction< idx_t >&         src,
                                                                     const P2P1TaylorHoodFunction< idx_t >&         dst,
                                                                     const std::shared_ptr< SparseMatrixProxy >&    mat,
                                                                     size_t                                         level,
                                                                     DoFType                                        flag )
{
   createMatrix( opr.A, src.uvw[0], dst.uvw[0], mat, level, flag );
   createMatrix( opr.A, src.uvw[1], dst.uvw[1], mat, level, flag );

   if ( src.uvw[0].getStorage()->hasGlobalCells() )
   {
      createMatrix( opr.A, src.uvw[2], dst.uvw[2], mat, level, flag );
   }

   createMatrix( opr.P, src.p, dst.p, mat, level, flag );
}
} // namespace petsc
} // namespace hyteg

#endif

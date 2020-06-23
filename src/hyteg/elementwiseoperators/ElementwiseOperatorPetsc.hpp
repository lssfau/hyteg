/*
 * Copyright (c) 2017-2020 Marcus Mohr.
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

#include "hyteg/composites/petsc/P2P1TaylorHoodPetsc.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseConstantCoefficientStokesOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

namespace hyteg {
namespace petsc {

/// Version of createMatrix function for P1ElementwiseOperators
template < class ElemOp, class FormType >
inline void createMatrix( const P1ElementwiseOperator< FormType >&    opr,
                          const P1Function< PetscInt >&               src,
                          const P1Function< PetscInt >&               dst,
                          const std::shared_ptr< SparseMatrixProxy >& mat,
                          uint_t                                      level,
                          DoFType                                     flag )
{
   opr.assembleLocalMatrix( mat, src, dst, level, flag );
}

/// Version of createMatrix function for P2ElementwiseOperators
template < class ElemOp, class FormType >
inline void createMatrix( const P2ElementwiseOperator< FormType >&    opr,
                          const P2Function< PetscInt >&               src,
                          const P2Function< PetscInt >&               dst,
                          const std::shared_ptr< SparseMatrixProxy >& mat,
                          uint_t                                      level,
                          DoFType                                     flag )
{
   opr.assembleLocalMatrix( mat, src, dst, level, flag );
}

/// Version of createMatrix function for P2P1ElementwiseConstantCoefficientStokesOperator
template <>
inline void createMatrix( const P2P1ElementwiseConstantCoefficientStokesOperator& opr,
                          const P2P1TaylorHoodFunction< PetscInt >&               src,
                          const P2P1TaylorHoodFunction< PetscInt >&               dst,
                          const std::shared_ptr< SparseMatrixProxy >&             mat,
                          uint_t                                                  level,
                          DoFType                                                 flag )
{
   opr.A.assembleLocalMatrix( mat, src.u, dst.u, level, flag );
   opr.A.assembleLocalMatrix( mat, src.v, dst.v, level, flag );

   opr.divT_x.assembleLocalMatrix( mat, src.p, dst.u, level, flag );
   opr.divT_y.assembleLocalMatrix( mat, src.p, dst.v, level, flag );

   opr.div_x.assembleLocalMatrix( mat, src.u, dst.p, level, flag );
   opr.div_y.assembleLocalMatrix( mat, src.v, dst.p, level, flag );

   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.A.assembleLocalMatrix( mat, src.w, dst.w, level, flag );
      opr.divT_z.assembleLocalMatrix( mat, src.p, dst.w, level, flag );
      opr.div_z.assembleLocalMatrix( mat, src.w, dst.p, level, flag );
   }
}

/// Version of createMatrix function for P2P1ElementwiseBlendingStokesOperator
template <>
inline void createMatrix( const P2P1ElementwiseBlendingStokesOperator& opr,
                          const P2P1TaylorHoodFunction< PetscInt >&    src,
                          const P2P1TaylorHoodFunction< PetscInt >&    dst,
                          const std::shared_ptr< SparseMatrixProxy >&  mat,
                          uint_t                                       level,
                          DoFType                                      flag )
{
   opr.A.assembleLocalMatrix( mat, src.u, dst.u, level, flag );
   opr.A.assembleLocalMatrix( mat, src.v, dst.v, level, flag );

   opr.divT_x.assembleLocalMatrix( mat, src.p, dst.u, level, flag );
   opr.divT_y.assembleLocalMatrix( mat, src.p, dst.v, level, flag );

   opr.div_x.assembleLocalMatrix( mat, src.u, dst.p, level, flag );
   opr.div_y.assembleLocalMatrix( mat, src.v, dst.p, level, flag );

   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.A.assembleLocalMatrix( mat, src.w, dst.w, level, flag );
      opr.divT_z.assembleLocalMatrix( mat, src.p, dst.w, level, flag );
      opr.div_z.assembleLocalMatrix( mat, src.w, dst.p, level, flag );
   }
}

} // namespace petsc
} // namespace hyteg

#endif

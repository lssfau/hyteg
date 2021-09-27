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

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/petsc/P1StokesPetsc.hpp"
#include "hyteg/composites/petsc/P2P1TaylorHoodPetsc.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P1P1ElementwiseAffineEpsilonStokesOperator.hpp"
#include "hyteg/elementwiseoperators/P1P1ElementwiseAffineEpsilonStokesBlockPreconditioner.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesBlockPreconditioner.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseConstantCoefficientStokesOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseAffineEpsilonStokesOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

namespace hyteg {
namespace petsc {

/// Version of createMatrix function for P1ElementwiseOperators
template < class ElemOp, class FormType >
inline void createMatrix( const P1ElementwiseOperator< FormType >&    opr,
                          const P1Function< idx_t >&                  src,
                          const P1Function< idx_t >&                  dst,
                          const std::shared_ptr< SparseMatrixProxy >& mat,
                          uint_t                                      level,
                          DoFType                                     flag )
{
   opr.assembleLocalMatrix( mat, src, dst, level, flag );
}

template <>
inline void createMatrix( const P1P1ElementwiseAffineEpsilonStokesOperator& opr,
                          const P1StokesFunction< idx_t >&                  src,
                          const P1StokesFunction< idx_t >&                  dst,
                          const std::shared_ptr< SparseMatrixProxy >&       mat,
                          uint_t                                            level,
                          DoFType                                           flag )
{
   opr.A_0_0.assembleLocalMatrix( mat, src.uvw[0], dst.uvw[0], level, flag );
   opr.A_0_1.assembleLocalMatrix( mat, src.uvw[1], dst.uvw[0], level, flag );
   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.A_0_2.assembleLocalMatrix( mat, src.uvw[2], dst.uvw[0], level, flag );
   }
   opr.divT_x.assembleLocalMatrix( mat, src.p, dst.uvw[0], level, flag );

   opr.A_1_0.assembleLocalMatrix( mat, src.uvw[0], dst.uvw[1], level, flag );
   opr.A_1_1.assembleLocalMatrix( mat, src.uvw[1], dst.uvw[1], level, flag );
   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.A_1_2.assembleLocalMatrix( mat, src.uvw[2], dst.uvw[1], level, flag );
   }
   opr.divT_y.assembleLocalMatrix( mat, src.p, dst.uvw[1], level, flag );

   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.A_2_0.assembleLocalMatrix( mat, src.uvw[0], dst.uvw[2], level, flag );
      opr.A_2_1.assembleLocalMatrix( mat, src.uvw[1], dst.uvw[2], level, flag );
      opr.A_2_2.assembleLocalMatrix( mat, src.uvw[2], dst.uvw[2], level, flag );

      opr.divT_z.assembleLocalMatrix( mat, src.p, dst.uvw[2], level, flag );
   }

   opr.div_x.assembleLocalMatrix( mat, src.uvw[0], dst.p, level, flag );
   opr.div_y.assembleLocalMatrix( mat, src.uvw[1], dst.p, level, flag );
   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.div_z.assembleLocalMatrix( mat, src.uvw[2], dst.p, level, flag );
   }

   opr.pspg.assembleLocalMatrix( mat, src.p, dst.p, level, flag );
}

template <>
inline void createMatrix( const P1P1ElementwiseAffineEpsilonStokesBlockPreconditioner& opr,
                          const P1StokesFunction< idx_t >&                             src,
                          const P1StokesFunction< idx_t >&                             dst,
                          const std::shared_ptr< SparseMatrixProxy >&                  mat,
                          uint_t                                                       level,
                          DoFType                                                      flag )
{
   opr.A_0_0.assembleLocalMatrix( mat, src.uvw[0], dst.uvw[0], level, flag );
   opr.A_0_1.assembleLocalMatrix( mat, src.uvw[1], dst.uvw[0], level, flag );
   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.A_0_2.assembleLocalMatrix( mat, src.uvw[2], dst.uvw[0], level, flag );
   }

   opr.A_1_0.assembleLocalMatrix( mat, src.uvw[0], dst.uvw[1], level, flag );
   opr.A_1_1.assembleLocalMatrix( mat, src.uvw[1], dst.uvw[1], level, flag );
   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.A_1_2.assembleLocalMatrix( mat, src.uvw[2], dst.uvw[1], level, flag );
   }

   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.A_2_0.assembleLocalMatrix( mat, src.uvw[0], dst.uvw[2], level, flag );
      opr.A_2_1.assembleLocalMatrix( mat, src.uvw[1], dst.uvw[2], level, flag );
      opr.A_2_2.assembleLocalMatrix( mat, src.uvw[2], dst.uvw[2], level, flag );
   }

   createMatrix( opr.P, src.p, dst.p, mat, level, flag );
}

/// Version of createMatrix function for P2ElementwiseOperators
template < class ElemOp, class FormType >
inline void createMatrix( const P2ElementwiseOperator< FormType >&    opr,
                          const P2Function< idx_t >&                  src,
                          const P2Function< idx_t >&                  dst,
                          const std::shared_ptr< SparseMatrixProxy >& mat,
                          uint_t                                      level,
                          DoFType                                     flag )
{
   opr.assembleLocalMatrix( mat, src, dst, level, flag );
}

/// Version of createMatrix function for P2P1ElementwiseConstantCoefficientStokesOperator
template <>
inline void createMatrix( const P2P1ElementwiseConstantCoefficientStokesOperator& opr,
                          const P2P1TaylorHoodFunction< idx_t >&                  src,
                          const P2P1TaylorHoodFunction< idx_t >&                  dst,
                          const std::shared_ptr< SparseMatrixProxy >&             mat,
                          uint_t                                                  level,
                          DoFType                                                 flag )
{
   opr.A.assembleLocalMatrix( mat, src.uvw[0], dst.uvw[0], level, flag );
   opr.A.assembleLocalMatrix( mat, src.uvw[1], dst.uvw[1], level, flag );

   opr.divT_x.assembleLocalMatrix( mat, src.p, dst.uvw[0], level, flag );
   opr.divT_y.assembleLocalMatrix( mat, src.p, dst.uvw[1], level, flag );

   opr.div_x.assembleLocalMatrix( mat, src.uvw[0], dst.p, level, flag );
   opr.div_y.assembleLocalMatrix( mat, src.uvw[1], dst.p, level, flag );

   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.A.assembleLocalMatrix( mat, src.uvw[2], dst.uvw[2], level, flag );
      opr.divT_z.assembleLocalMatrix( mat, src.p, dst.uvw[2], level, flag );
      opr.div_z.assembleLocalMatrix( mat, src.uvw[2], dst.p, level, flag );
   }
}

/// Version of createMatrix function for P2P1ElementwiseBlendingStokesOperator
template <>
inline void createMatrix( const P2P1ElementwiseBlendingStokesOperator& opr,
                          const P2P1TaylorHoodFunction< idx_t >&       src,
                          const P2P1TaylorHoodFunction< idx_t >&       dst,
                          const std::shared_ptr< SparseMatrixProxy >&  mat,
                          uint_t                                       level,
                          DoFType                                      flag )
{
   opr.A.assembleLocalMatrix( mat, src.uvw[0], dst.uvw[0], level, flag );
   opr.A.assembleLocalMatrix( mat, src.uvw[1], dst.uvw[1], level, flag );

   opr.divT_x.assembleLocalMatrix( mat, src.p, dst.uvw[0], level, flag );
   opr.divT_y.assembleLocalMatrix( mat, src.p, dst.uvw[1], level, flag );

   opr.div_x.assembleLocalMatrix( mat, src.uvw[0], dst.p, level, flag );
   opr.div_y.assembleLocalMatrix( mat, src.uvw[1], dst.p, level, flag );

   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.A.assembleLocalMatrix( mat, src.uvw[2], dst.uvw[2], level, flag );
      opr.divT_z.assembleLocalMatrix( mat, src.p, dst.uvw[2], level, flag );
      opr.div_z.assembleLocalMatrix( mat, src.uvw[2], dst.p, level, flag );
   }
}

template <>
inline void
    createMatrix< P2P1ElementwiseBlendingStokesBlockPreconditioner >( const P2P1ElementwiseBlendingStokesBlockPreconditioner& opr,
                                                                      const P2P1TaylorHoodFunction< idx_t >&                  src,
                                                                      const P2P1TaylorHoodFunction< idx_t >&                  dst,
                                                                      const std::shared_ptr< SparseMatrixProxy >&             mat,
                                                                      size_t  level,
                                                                      DoFType flag )
{
   for ( uint_t dim = 0; dim < src.uvw.getDimension(); dim++ )
   {
      // need to help the compiler here to find the correct version of createMatrix
      createMatrix< P2ElementwiseBlendingLaplaceOperator, P2Form_laplace >( opr.A, src.uvw[dim], dst.uvw[dim], mat, level, flag );
   }
   createMatrix( opr.P, src.p, dst.p, mat, level, flag );
}

template <>
inline void createMatrix( const P2P1ElementwiseAffineEpsilonStokesOperator& opr,
                          const P2P1TaylorHoodFunction< idx_t >&            src,
                          const P2P1TaylorHoodFunction< idx_t >&            dst,
                          const std::shared_ptr< SparseMatrixProxy >&       mat,
                          uint_t                                            level,
                          DoFType                                           flag )
{
   opr.A_0_0.assembleLocalMatrix( mat, src.uvw[0], dst.uvw[0], level, flag );
   opr.A_0_1.assembleLocalMatrix( mat, src.uvw[1], dst.uvw[0], level, flag );
   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.A_0_2.assembleLocalMatrix( mat, src.uvw[2], dst.uvw[0], level, flag );
   }
   opr.divT_x.assembleLocalMatrix( mat, src.p, dst.uvw[0], level, flag );

   opr.A_1_0.assembleLocalMatrix( mat, src.uvw[0], dst.uvw[1], level, flag );
   opr.A_1_1.assembleLocalMatrix( mat, src.uvw[1], dst.uvw[1], level, flag );
   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.A_1_2.assembleLocalMatrix( mat, src.uvw[2], dst.uvw[1], level, flag );
   }
   opr.divT_y.assembleLocalMatrix( mat, src.p, dst.uvw[1], level, flag );

   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.A_2_0.assembleLocalMatrix( mat, src.uvw[0], dst.uvw[2], level, flag );
      opr.A_2_1.assembleLocalMatrix( mat, src.uvw[1], dst.uvw[2], level, flag );
      opr.A_2_2.assembleLocalMatrix( mat, src.uvw[2], dst.uvw[2], level, flag );

      opr.divT_z.assembleLocalMatrix( mat, src.p, dst.uvw[2], level, flag );
   }

   opr.div_x.assembleLocalMatrix( mat, src.uvw[0], dst.p, level, flag );
   opr.div_y.assembleLocalMatrix( mat, src.uvw[1], dst.p, level, flag );
   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.div_z.assembleLocalMatrix( mat, src.uvw[2], dst.p, level, flag );
   }
}

template <>
inline void createMatrix( const P2P1ElementwiseAffineEpsilonStokesBlockPreconditioner& opr,
                          const P2P1TaylorHoodFunction< idx_t >&                       src,
                          const P2P1TaylorHoodFunction< idx_t >&                       dst,
                          const std::shared_ptr< SparseMatrixProxy >&                  mat,
                          uint_t                                                       level,
                          DoFType                                                      flag )
{
   opr.A_0_0.assembleLocalMatrix( mat, src.uvw[0], dst.uvw[0], level, flag );
   opr.A_0_1.assembleLocalMatrix( mat, src.uvw[1], dst.uvw[0], level, flag );
   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.A_0_2.assembleLocalMatrix( mat, src.uvw[2], dst.uvw[0], level, flag );
   }

   opr.A_1_0.assembleLocalMatrix( mat, src.uvw[0], dst.uvw[1], level, flag );
   opr.A_1_1.assembleLocalMatrix( mat, src.uvw[1], dst.uvw[1], level, flag );
   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.A_1_2.assembleLocalMatrix( mat, src.uvw[2], dst.uvw[1], level, flag );
   }

   if ( src.getStorage()->hasGlobalCells() )
   {
      opr.A_2_0.assembleLocalMatrix( mat, src.uvw[0], dst.uvw[2], level, flag );
      opr.A_2_1.assembleLocalMatrix( mat, src.uvw[1], dst.uvw[2], level, flag );
      opr.A_2_2.assembleLocalMatrix( mat, src.uvw[2], dst.uvw[2], level, flag );
   }

   createMatrix( opr.P, src.p, dst.p, mat, level, flag );
}

} // namespace petsc
} // namespace hyteg

#endif

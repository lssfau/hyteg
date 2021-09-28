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
#include "hyteg/mixedoperators/P1ScalarToP2VectorOperator.hpp"
#include "hyteg/mixedoperators/P1ToP2Operator.hpp"
#include "hyteg/mixedoperators/P2ToP1Operator.hpp"
#include "hyteg/mixedoperators/P2VectorToP1ScalarOperator.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"

namespace hyteg {

class P2P1TaylorHoodStokesOperator : public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef P2ConstantLaplaceOperator               VelocityOperator_T;
   typedef P2P1TaylorHoodStokesBlockPreconditioner BlockPreconditioner_T;

   P2P1TaylorHoodStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , A( storage, minLevel, maxLevel )
   , Lapl( storage, minLevel, maxLevel )
   , div_x( storage, minLevel, maxLevel )
   , div_y( storage, minLevel, maxLevel )
   , div_z( storage, minLevel, maxLevel )
   , div( storage, minLevel, maxLevel )
   , divT_x( storage, minLevel, maxLevel )
   , divT_y( storage, minLevel, maxLevel )
   , divT_z( storage, minLevel, maxLevel )
   , divT( storage, minLevel, maxLevel )
   , pspg_( storage, minLevel, maxLevel )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               const uint_t                            level,
               const DoFType                           flag ) const
   {
      Lapl.apply( src.uvw, dst.uvw, level, flag, Replace );
      divT.apply( src.p, dst.uvw, level, flag, Add );
      div.apply( src.uvw, dst.p, level, flag, Replace );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P1TaylorHoodFunction< matIdx_t >&   src,
                  const P2P1TaylorHoodFunction< matIdx_t >&   dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      A.toMatrix( mat, src.uvw[0], dst.uvw[0], level, flag );
      divT_x.getVertexToVertexOpr().toMatrix( mat, src.p, dst.uvw[0].getVertexDoFFunction(), level, flag );
      divT_x.getVertexToEdgeOpr().toMatrix( mat, src.p, dst.uvw[0].getEdgeDoFFunction(), level, flag );

      A.toMatrix( mat, src.uvw[1], dst.uvw[1], level, flag );
      divT_y.getVertexToVertexOpr().toMatrix( mat, src.p, dst.uvw[1].getVertexDoFFunction(), level, flag );
      divT_y.getVertexToEdgeOpr().toMatrix( mat, src.p, dst.uvw[1].getEdgeDoFFunction(), level, flag );

      if ( src.uvw[0].getStorage()->hasGlobalCells() )
      {
         A.toMatrix( mat, src.uvw[2], dst.uvw[2], level, flag );
         divT_z.getVertexToVertexOpr().toMatrix( mat, src.p, dst.uvw[2].getVertexDoFFunction(), level, flag );
         divT_z.getVertexToEdgeOpr().toMatrix( mat, src.p, dst.uvw[2].getEdgeDoFFunction(), level, flag );
      }

      div_x.getVertexToVertexOpr().toMatrix( mat, src.uvw[0].getVertexDoFFunction(), dst.p, level, flag | DirichletBoundary );
      div_x.getEdgeToVertexOpr().toMatrix( mat, src.uvw[0].getEdgeDoFFunction(), dst.p, level, flag | DirichletBoundary );

      div_y.getVertexToVertexOpr().toMatrix( mat, src.uvw[1].getVertexDoFFunction(), dst.p, level, flag | DirichletBoundary );
      div_y.getEdgeToVertexOpr().toMatrix( mat, src.uvw[1].getEdgeDoFFunction(), dst.p, level, flag | DirichletBoundary );

      if ( src.uvw[0].getStorage()->hasGlobalCells() )
      {
         div_z.getVertexToVertexOpr().toMatrix( mat, src.uvw[2].getVertexDoFFunction(), dst.p, level, flag | DirichletBoundary );
         div_z.getEdgeToVertexOpr().toMatrix( mat, src.uvw[2].getEdgeDoFFunction(), dst.p, level, flag | DirichletBoundary );
      }
   }

   P2ConstantVectorLaplaceOperator Lapl;
   P2ToP1ConstantDivOperator       div;
   P1ToP2ConstantDivTOperator      divT;

   // currently need these for being able to call createMatrix()
   P2ConstantLaplaceOperator   A;
   P2ToP1ConstantDivxOperator  div_x;
   P2ToP1ConstantDivyOperator  div_y;
   P2ToP1ConstantDivzOperator  div_z;
   P1ToP2ConstantDivTxOperator divT_x;
   P1ToP2ConstantDivTyOperator divT_y;
   P1ToP2ConstantDivTzOperator divT_z;

   /// this operator is need in the uzawa smoother
   P1PSPGOperator        pspg_;
   P1PSPGInvDiagOperator pspg_inv_diag_;
   bool                  hasGlobalCells_;
};

} // namespace hyteg

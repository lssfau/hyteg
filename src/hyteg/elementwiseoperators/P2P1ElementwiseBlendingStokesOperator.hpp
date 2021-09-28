/*
 * Copyright (c) 2017-2020 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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
#include "hyteg/elementwiseoperators/P1ToP2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesBlockPreconditioner.hpp"
#include "hyteg/elementwiseoperators/P2ToP1ElementwiseOperator.hpp"

namespace hyteg {

class P2P1ElementwiseBlendingStokesOperator
: public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef P2ElementwiseBlendingLaplaceOperator VelocityOperator_T;

   // This typedef could be a temporary solution to allow for using
   // the PETScBlockPreconditionedStokesSolver and for compiling
   // the StrongFreeSlipWrapper for the P2P1ElementwiseBlendingStokesOperator
   // typedef P2P1TaylorHoodStokesBlockPreconditioner BlockPreconditioner_T;

   typedef P2P1ElementwiseBlendingStokesBlockPreconditioner BlockPreconditioner_T;

   P2P1ElementwiseBlendingStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , A( storage, minLevel, maxLevel )
   , div_x( storage, minLevel, maxLevel )
   , div_y( storage, minLevel, maxLevel )
   , div_z( storage, minLevel, maxLevel )
   , divT_x( storage, minLevel, maxLevel )
   , divT_y( storage, minLevel, maxLevel )
   , divT_z( storage, minLevel, maxLevel )
   //   , pspg_( storage, minLevel, maxLevel )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void computeAndStoreLocalElementMatrices()
   {
      A.computeAndStoreLocalElementMatrices();

      div_x.computeAndStoreLocalElementMatrices();
      div_y.computeAndStoreLocalElementMatrices();

      divT_x.computeAndStoreLocalElementMatrices();
      divT_y.computeAndStoreLocalElementMatrices();

      if ( hasGlobalCells_ )
      {
         div_z.computeAndStoreLocalElementMatrices();
         divT_z.computeAndStoreLocalElementMatrices();
      }
   }

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               const uint_t                            level,
               const DoFType                           flag ) const
   {
      A.apply( src.uvw[0], dst.uvw[0], level, flag, Replace );
      divT_x.apply( src.p, dst.uvw[0], level, flag, Add );

      A.apply( src.uvw[1], dst.uvw[1], level, flag, Replace );
      divT_y.apply( src.p, dst.uvw[1], level, flag, Add );

      if ( hasGlobalCells_ )
      {
         A.apply( src.uvw[2], dst.uvw[2], level, flag, Replace );
         divT_z.apply( src.p, dst.uvw[2], level, flag, Add );
      }

      div_x.apply( src.uvw[0], dst.p, level, flag, Replace );
      div_y.apply( src.uvw[1], dst.p, level, flag, Add );
      if ( hasGlobalCells_ )
      {
         div_z.apply( src.uvw[2], dst.p, level, flag, Add );
      }
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P1TaylorHoodFunction< matIdx_t >&   src,
                  const P2P1TaylorHoodFunction< matIdx_t >&   dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      A.toMatrix( mat, src.uvw[0], dst.uvw[0], level, flag );
      A.toMatrix( mat, src.uvw[1], dst.uvw[1], level, flag );

      divT_x.toMatrix( mat, src.p, dst.uvw[0], level, flag );
      divT_y.toMatrix( mat, src.p, dst.uvw[1], level, flag );

      div_x.toMatrix( mat, src.uvw[0], dst.p, level, flag );
      div_y.toMatrix( mat, src.uvw[1], dst.p, level, flag );

      if ( src.getStorage()->hasGlobalCells() )
      {
         A.toMatrix( mat, src.uvw[2], dst.uvw[2], level, flag );
         divT_z.toMatrix( mat, src.p, dst.uvw[2], level, flag );
         div_z.toMatrix( mat, src.uvw[2], dst.p, level, flag );
      }
   }

   P2ElementwiseBlendingLaplaceOperator   A;
   P2ToP1ElementwiseBlendingDivxOperator  div_x;
   P2ToP1ElementwiseBlendingDivyOperator  div_y;
   P2ToP1ElementwiseBlendingDivzOperator  div_z;
   P1ToP2ElementwiseBlendingDivTxOperator divT_x;
   P1ToP2ElementwiseBlendingDivTyOperator divT_y;
   P1ToP2ElementwiseBlendingDivTzOperator divT_z;

   /// this operator is need in the uzawa smoother
   //   P1ElementwisePSPGOperator        pspg_;
   P1PSPGInvDiagOperator pspg_inv_diag_;
   bool                  hasGlobalCells_;
};

} // namespace hyteg

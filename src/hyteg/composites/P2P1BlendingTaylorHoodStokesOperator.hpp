/*
 * Copyright (c) 2017-2019 Benjamin Mann.
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
#include "hyteg/mixedoperators/P1ToP2VariableOperator.hpp"
#include "hyteg/mixedoperators/P2ToP1VariableOperator.hpp"
#include "hyteg/p2functionspace/P2VariableOperator.hpp"

namespace hyteg {

class P2P1BlendingTaylorHoodStokesOperator : public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef P2BlendingLaplaceOperator VelocityOperator_T;

   P2P1BlendingTaylorHoodStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , A( storage, minLevel, maxLevel )
   , div_x( storage, minLevel, maxLevel )
   , div_y( storage, minLevel, maxLevel )
   , div_z( storage, minLevel, maxLevel )
   , divT_x( storage, minLevel, maxLevel )
   , divT_y( storage, minLevel, maxLevel )
   , divT_z( storage, minLevel, maxLevel )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               const size_t                            level,
               DoFType                                 flag ) const
   {
      WALBERLA_CHECK( !hasGlobalCells_, "Variable Stokes operator not implemented for 3D." );

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

   P2BlendingLaplaceOperator   A;
   P2ToP1BlendingDivxOperator  div_x;
   P2ToP1BlendingDivyOperator  div_y;
   P2ToP1BlendingDivzOperator  div_z;
   P1ToP2BlendingDivTxOperator divT_x;
   P1ToP2BlendingDivTyOperator divT_y;
   P1ToP2BlendingDivTzOperator divT_z;

   /// this operator is need in the uzawa smoother
   // P1PSPGOperator        pspg_;
   P1PSPGInvDiagOperator pspg_inv_diag_;
   bool                  hasGlobalCells_;
};

} // namespace hyteg

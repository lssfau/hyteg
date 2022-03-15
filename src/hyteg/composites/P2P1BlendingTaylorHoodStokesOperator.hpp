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
#include "hyteg/operators/ScalarToVectorOperator.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/operators/VectorToScalarOperator.hpp"
#include "hyteg/p2functionspace/P2VariableOperator.hpp"

namespace hyteg {

class P2P1BlendingTaylorHoodStokesOperator : public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef P2BlendingLaplaceOperator VelocityOperator_T;

   P2P1BlendingTaylorHoodStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , Lapl( storage, minLevel, maxLevel )
   , div( storage, minLevel, maxLevel )
   , divT( storage, minLevel, maxLevel )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               const size_t                            level,
               DoFType                                 flag ) const
   {
      WALBERLA_CHECK( !hasGlobalCells_, "Variable Stokes operator not implemented for 3D." );

      Lapl.apply( src.uvw(), dst.uvw(), level, flag, Replace );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag, Replace );
   }

   const P2BlendingLaplaceOperator& getA() const
   {
      auto ptr = Lapl.getSubOperator( 0, 0 );
      return dynamic_cast< const P2BlendingLaplaceOperator& >( *ptr );
   }

   P2BlendingVectorLaplaceOperator Lapl;
   P2ToP1VariableDivOperator       div;
   P1ToP2VariableDivTOperator      divT;

   /// this operator is need in the uzawa smoother
   P1PSPGInvDiagOperator pspg_inv_diag_;

   bool hasGlobalCells_;
};

} // namespace hyteg

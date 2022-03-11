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
#include "hyteg/mixedoperators/P1ToP2SurrogateOperator.hpp"
#include "hyteg/mixedoperators/P2ToP1SurrogateOperator.hpp"
#include "hyteg/p2functionspace/P2SurrogateOperator.hpp"

namespace hyteg {

class P2P1SurrogateTaylorHoodStokesOperator
: public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef P2SurrogateLaplaceOperator VelocityOperator_T;

   P2P1SurrogateTaylorHoodStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                          uint_t                                     minLevel,
                                          uint_t                                     maxLevel,
                                          uint_t                                     interpolationLevel )
   : Operator( storage, minLevel, maxLevel )
   , Lapl( storage, minLevel, maxLevel, interpolationLevel )
   , divT( storage, minLevel, maxLevel, interpolationLevel )
   , div( storage, minLevel, maxLevel, interpolationLevel )
   , div_x( storage, minLevel, maxLevel, interpolationLevel )
   , div_y( storage, minLevel, maxLevel, interpolationLevel )
   , div_z( storage, minLevel, maxLevel )
   , divT_x( storage, minLevel, maxLevel, interpolationLevel )
   , divT_y( storage, minLevel, maxLevel, interpolationLevel )
   , divT_z( storage, minLevel, maxLevel )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void interpolateStencils( uint_t polyDegree )
   {
      this->getA().interpolateStencils( polyDegree );
      div_x.interpolateStencils( polyDegree );
      div_y.interpolateStencils( polyDegree );
      divT_x.interpolateStencils( polyDegree );
      divT_y.interpolateStencils( polyDegree );
   }

   void useDegree( uint_t polyDegree )
   {
      this->getA().useDegree( polyDegree );
      div_x.useDegree( polyDegree );
      div_y.useDegree( polyDegree );
      divT_x.useDegree( polyDegree );
      divT_y.useDegree( polyDegree );
   }

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               const size_t                            level,
               DoFType                                 flag ) const
   {
      WALBERLA_CHECK( !hasGlobalCells_, "Surrogate Stokes operator not implemented for 3D." );

      Lapl.apply( src.uvw(), dst.uvw(), level, flag, Replace );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag, Replace );
   }

   const P2SurrogateLaplaceOperator& getA() const
   {
      auto ptr = Lapl.getSubOperator( 0, 0 );
      return dynamic_cast< const P2SurrogateLaplaceOperator& >( *ptr );
   }

   P2SurrogateVectorLaplaceOperator Lapl;
   P1ToP2SurrogateDivTOperator      divT;
   P2ToP1SurrogateDivOperator       div;

   P2ToP1SurrogateDivxOperator  div_x;
   P2ToP1SurrogateDivyOperator  div_y;
   P2ToP1BlendingDivzOperator   div_z;
   P1ToP2SurrogateDivTxOperator divT_x;
   P1ToP2SurrogateDivTyOperator divT_y;
   P1ToP2BlendingDivTzOperator  divT_z;

   /// this operator is need in the uzawa smoother
   // P1PSPGOperator        pspg_;
   P1PSPGInvDiagOperator pspg_inv_diag_;
   bool                  hasGlobalCells_;

private:
   P2SurrogateLaplaceOperator& getA()
   {
      auto ptr = Lapl.getSubOperator( 0, 0 );
      return dynamic_cast< P2SurrogateLaplaceOperator& >( *ptr );
   }
};

} // namespace hyteg

/*
 * Copyright (c) 2017-2021 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseAffineEpsilonStokesBlockPreconditioner.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_epsilonvar_0_0_affine_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_epsilonvar_0_1_affine_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_epsilonvar_0_2_affine_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_epsilonvar_1_0_affine_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_epsilonvar_1_1_affine_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_epsilonvar_1_2_affine_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_epsilonvar_2_0_affine_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_epsilonvar_2_1_affine_q2.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_epsilonvar_2_2_affine_q2.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {

class P2P1ElementwiseAffineEpsilonStokesOperator
: public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef P2P1ElementwiseAffineEpsilonStokesBlockPreconditioner BlockPreconditioner_T;

   P2P1ElementwiseAffineEpsilonStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                               uint_t                                     minLevel,
                                               uint_t                                     maxLevel,
                                               std::function< real_t( const Point3D& ) >  mu )
   : Operator( storage, minLevel, maxLevel )
   , A_0_0( storage, minLevel, maxLevel, forms::p2_epsilonvar_0_0_affine_q2( mu, mu ) )
   , A_0_1( storage, minLevel, maxLevel, forms::p2_epsilonvar_0_1_affine_q2( mu, mu ) )
   , A_0_2( storage, minLevel, maxLevel, forms::p2_epsilonvar_0_2_affine_q2( mu ) )
   , A_1_0( storage, minLevel, maxLevel, forms::p2_epsilonvar_1_0_affine_q2( mu, mu ) )
   , A_1_1( storage, minLevel, maxLevel, forms::p2_epsilonvar_1_1_affine_q2( mu, mu ) )
   , A_1_2( storage, minLevel, maxLevel, forms::p2_epsilonvar_1_2_affine_q2( mu ) )
   , A_2_0( storage, minLevel, maxLevel, forms::p2_epsilonvar_2_0_affine_q2( mu ) )
   , A_2_1( storage, minLevel, maxLevel, forms::p2_epsilonvar_2_1_affine_q2( mu ) )
   , A_2_2( storage, minLevel, maxLevel, forms::p2_epsilonvar_2_2_affine_q2( mu ) )
   , div_x( storage, minLevel, maxLevel )
   , div_y( storage, minLevel, maxLevel )
   , div_z( storage, minLevel, maxLevel )
   , divT_x( storage, minLevel, maxLevel )
   , divT_y( storage, minLevel, maxLevel )
   , divT_z( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   P2P1ElementwiseAffineEpsilonStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                               uint_t                                     minLevel,
                                               uint_t                                     maxLevel )
   : P2P1ElementwiseAffineEpsilonStokesOperator( storage, minLevel, maxLevel, []( const Point3D& ) { return real_c( 1 ); } )
   {}

   void computeAndStoreLocalElementMatrices() { WALBERLA_ABORT( "Not implemented." ) }

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               const uint_t                            level,
               const DoFType                           flag ) const
   {
      A_0_0.apply( src.uvw[0], dst.uvw[0], level, flag );
      A_0_1.apply( src.uvw[1], dst.uvw[0], level, flag, Add );
      if ( hasGlobalCells_ )
      {
         A_0_2.apply( src.uvw[2], dst.uvw[0], level, flag, Add );
      }

      divT_x.apply( src.p, dst.uvw[0], level, flag, Add );

      A_1_0.apply( src.uvw[0], dst.uvw[1], level, flag );
      A_1_1.apply( src.uvw[1], dst.uvw[1], level, flag, Add );
      if ( hasGlobalCells_ )
      {
         A_1_2.apply( src.uvw[2], dst.uvw[1], level, flag, Add );
      }
      divT_y.apply( src.p, dst.uvw[1], level, flag, Add );

      if ( hasGlobalCells_ )
      {
         A_2_0.apply( src.uvw[0], dst.uvw[2], level, flag );
         A_2_1.apply( src.uvw[1], dst.uvw[2], level, flag, Add );
         A_2_2.apply( src.uvw[2], dst.uvw[2], level, flag, Add );
         divT_z.apply( src.p, dst.uvw[2], level, flag, Add );
      }

      div_x.apply( src.uvw[0], dst.p, level, flag );
      div_y.apply( src.uvw[1], dst.p, level, flag, Add );
      if ( hasGlobalCells_ )
      {
         div_z.apply( src.uvw[2], dst.p, level, flag, Add );
      }
   }

   P2ElementwiseOperator< forms::p2_epsilonvar_0_0_affine_q2 > A_0_0;
   P2ElementwiseOperator< forms::p2_epsilonvar_0_1_affine_q2 > A_0_1;
   P2ElementwiseOperator< forms::p2_epsilonvar_0_2_affine_q2 > A_0_2;

   P2ElementwiseOperator< forms::p2_epsilonvar_1_0_affine_q2 > A_1_0;
   P2ElementwiseOperator< forms::p2_epsilonvar_1_1_affine_q2 > A_1_1;
   P2ElementwiseOperator< forms::p2_epsilonvar_1_2_affine_q2 > A_1_2;

   P2ElementwiseOperator< forms::p2_epsilonvar_2_0_affine_q2 > A_2_0;
   P2ElementwiseOperator< forms::p2_epsilonvar_2_1_affine_q2 > A_2_1;
   P2ElementwiseOperator< forms::p2_epsilonvar_2_2_affine_q2 > A_2_2;

   P2ToP1ElementwiseDivxOperator div_x;
   P2ToP1ElementwiseDivyOperator div_y;
   P2ToP1ElementwiseDivzOperator div_z;

   P1ToP2ElementwiseDivTxOperator divT_x;
   P1ToP2ElementwiseDivTyOperator divT_y;
   P1ToP2ElementwiseDivTzOperator divT_z;

   bool hasGlobalCells_;
};

} // namespace hyteg

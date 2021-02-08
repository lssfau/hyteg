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

#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P1P1ElementwiseAffineEpsilonStokesBlockPreconditioner.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/P1EpsilonAffine_0_0.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/P1EpsilonAffine_0_1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/P1EpsilonAffine_0_2.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/P1EpsilonAffine_1_0.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/P1EpsilonAffine_1_1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/P1EpsilonAffine_1_2.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/P1EpsilonAffine_2_0.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/P1EpsilonAffine_2_1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/P1EpsilonAffine_2_2.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {

class P1P1ElementwiseAffineEpsilonStokesOperator : public Operator< P1StokesFunction< real_t >, P1StokesFunction< real_t > >
{
 public:

   typedef P1P1ElementwiseAffineEpsilonStokesBlockPreconditioner BlockPreconditioner_T;

   P1P1ElementwiseAffineEpsilonStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                               uint_t                                     minLevel,
                                               uint_t                                     maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , A_0_0( storage, minLevel, maxLevel )
   , A_0_1( storage, minLevel, maxLevel )
   , A_0_2( storage, minLevel, maxLevel )
   , A_1_0( storage, minLevel, maxLevel )
   , A_1_1( storage, minLevel, maxLevel )
   , A_1_2( storage, minLevel, maxLevel )
   , A_2_0( storage, minLevel, maxLevel )
   , A_2_1( storage, minLevel, maxLevel )
   , A_2_2( storage, minLevel, maxLevel )
   , div_x( storage, minLevel, maxLevel )
   , div_y( storage, minLevel, maxLevel )
   , div_z( storage, minLevel, maxLevel )
   , divT_x( storage, minLevel, maxLevel )
   , divT_y( storage, minLevel, maxLevel )
   , divT_z( storage, minLevel, maxLevel )
   , pspg( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void computeAndStoreLocalElementMatrices() { WALBERLA_ABORT( "Not implemented." ) }

   void apply( const P1StokesFunction< real_t >& src,
               const P1StokesFunction< real_t >& dst,
               const uint_t                      level,
               const DoFType                     flag ) const
   {
      A_0_0.apply( src.uvw.u, dst.uvw.u, level, flag );
      A_0_1.apply( src.uvw.v, dst.uvw.u, level, flag, Add );
      if ( hasGlobalCells_ )
      {
         A_0_2.apply( src.uvw.w, dst.uvw.u, level, flag, Add );
      }

      divT_x.apply( src.p, dst.uvw.u, level, flag, Add );

      A_1_0.apply( src.uvw.u, dst.uvw.v, level, flag );
      A_1_1.apply( src.uvw.v, dst.uvw.v, level, flag, Add );
      if ( hasGlobalCells_ )
      {
         A_1_2.apply( src.uvw.w, dst.uvw.v, level, flag, Add );
      }
      divT_y.apply( src.p, dst.uvw.v, level, flag, Add );

      if ( hasGlobalCells_ )
      {
         A_2_0.apply( src.uvw.u, dst.uvw.w, level, flag );
         A_2_1.apply( src.uvw.v, dst.uvw.w, level, flag, Add );
         A_2_2.apply( src.uvw.w, dst.uvw.w, level, flag, Add );
         divT_z.apply( src.p, dst.uvw.w, level, flag, Add );
      }

      div_x.apply( src.uvw.u, dst.p, level, flag );
      div_y.apply( src.uvw.v, dst.p, level, flag, Add );
      if ( hasGlobalCells_ )
      {
         div_z.apply( src.uvw.w, dst.p, level, flag, Add );
      }

      pspg.apply( src.p, dst.p, level, flag, Add );
   }

   P1ElementwiseOperator< P1EpsilonAffine_0_0 > A_0_0;
   P1ElementwiseOperator< P1EpsilonAffine_0_1 > A_0_1;
   P1ElementwiseOperator< P1EpsilonAffine_0_2 > A_0_2;

   P1ElementwiseOperator< P1EpsilonAffine_1_0 > A_1_0;
   P1ElementwiseOperator< P1EpsilonAffine_1_1 > A_1_1;
   P1ElementwiseOperator< P1EpsilonAffine_1_2 > A_1_2;

   P1ElementwiseOperator< P1EpsilonAffine_2_0 > A_2_0;
   P1ElementwiseOperator< P1EpsilonAffine_2_1 > A_2_1;
   P1ElementwiseOperator< P1EpsilonAffine_2_2 > A_2_2;

   P1ElementwiseDivXOperator div_x;
   P1ElementwiseDivYOperator div_y;
   P1ElementwiseDivZOperator div_z;

   P1ElementwiseDivTXOperator divT_x;
   P1ElementwiseDivTYOperator divT_y;
   P1ElementwiseDivTZOperator divT_z;

   P1ElementwisePSPGOperator pspg;

   bool hasGlobalCells_;
};

} // namespace hyteg

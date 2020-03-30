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
#include "hyteg/elementwiseoperators/P2ToP1ElementwiseOperator.hpp"

namespace hyteg {

class P2P1ElementwiseConstantCoefficientStokesOperator
: public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef P2ElementwiseLaplaceOperator VelocityOperator_T;

   P2P1ElementwiseConstantCoefficientStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                     size_t                                     minLevel,
                                                     size_t                                     maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , A( storage, minLevel, maxLevel )
   , div_x( storage, minLevel, maxLevel )
   , div_y( storage, minLevel, maxLevel )
   , div_z( storage, minLevel, maxLevel )
   , divT_x( storage, minLevel, maxLevel )
   , divT_y( storage, minLevel, maxLevel )
   , divT_z( storage, minLevel, maxLevel )
//   , pspg_( storage, minLevel, maxLevel )
//   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               const uint_t                            level,
               const DoFType                           flag ) const
   {
      A.apply( src.u, dst.u, level, flag, Replace );
      divT_x.apply( src.p, dst.u, level, flag, Add );

      A.apply( src.v, dst.v, level, flag, Replace );
      divT_y.apply( src.p, dst.v, level, flag, Add );

      if ( hasGlobalCells_ )
      {
         A.apply( src.w, dst.w, level, flag, Replace );
         divT_z.apply( src.p, dst.w, level, flag, Add );
      }

      div_x.apply( src.u, dst.p, level, flag, Replace );
      div_y.apply( src.v, dst.p, level, flag, Add );

      if ( hasGlobalCells_ )
      {
         div_z.apply( src.w, dst.p, level, flag, Add );
      }
   }

   P2ElementwiseLaplaceOperator   A;
   P2toP1ElementwiseDivxOperator  div_x;
   P2toP1ElementwiseDivyOperator  div_y;
   P2toP1ElementwiseDivzOperator  div_z;
   P1toP2ElementwiseDivTxOperator divT_x;
   P1toP2ElementwiseDivTyOperator divT_y;
   P1toP2ElementwiseDivTzOperator divT_z;

   /// this operator is need in the uzawa smoother
//   P1ElementwisePSPGOperator        pspg_;
//   P1ElementwisePSPGInvDiagOperator pspg_inv_diag_;
   bool                             hasGlobalCells_;
};

} // namespace hyteg

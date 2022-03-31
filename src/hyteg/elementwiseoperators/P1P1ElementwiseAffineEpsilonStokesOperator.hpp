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
#include "hyteg/forms/form_hyteg_generated/p1/p1_epsilon_all_forms.hpp"
#include "hyteg/p1functionspace/P1EpsilonOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {

class P1P1ElementwiseAffineEpsilonStokesOperator : public Operator< P1StokesFunction< real_t >, P1StokesFunction< real_t > >
{
 public:
   typedef P1P1ElementwiseAffineEpsilonStokesBlockPreconditioner BlockPreconditioner_T;

   P1P1ElementwiseAffineEpsilonStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                               uint_t                                     minLevel,
                                               uint_t                                     maxLevel,
                                               std::function< real_t( const Point3D& ) >  mu )
   : Operator( storage, minLevel, maxLevel )
   , viscOp( storage, minLevel, maxLevel, mu )
   , div( storage, minLevel, maxLevel )
   , divT( storage, minLevel, maxLevel )
   , pspg( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   P1P1ElementwiseAffineEpsilonStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                               uint_t                                     minLevel,
                                               uint_t                                     maxLevel )
   : P1P1ElementwiseAffineEpsilonStokesOperator( storage, minLevel, maxLevel, []( const Point3D& ) { return real_c( 1 ); } )
   {}

   void computeAndStoreLocalElementMatrices() { WALBERLA_ABORT( "Not implemented." ) }

   void apply( const P1StokesFunction< real_t >& src,
               const P1StokesFunction< real_t >& dst,
               const uint_t                      level,
               const DoFType                     flag ) const
   {
      viscOp.apply( src.uvw(), dst.uvw(), level, flag );

      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag, Replace );

      pspg.apply( src.p(), dst.p(), level, flag, Add );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1StokesFunction< idx_t >&            src,
                  const P1StokesFunction< idx_t >&            dst,
                  uint_t                                      level,
                  DoFType                                     flag ) const
   {
      viscOp.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );

      divT.toMatrix( mat, src.p(), dst.uvw(), level, flag );
      div.toMatrix( mat, src.uvw(), dst.p(), level, flag );

      pspg.toMatrix( mat, src.p(), dst.p(), level, flag );
   }

   P1ElementwiseAffineEpsilonOperator viscOp;
   P1ToP1ElementwiseDivOperator  div;
   P1ToP1ElementwiseDivTOperator divT;

   P1ElementwisePSPGOperator pspg;

   bool hasGlobalCells_;
};

} // namespace hyteg

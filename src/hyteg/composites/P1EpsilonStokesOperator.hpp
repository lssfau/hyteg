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

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/StokesOperatorTraits.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"

namespace hyteg {

class P1EpsilonStokesOperator : public Operator< P1StokesFunction< real_t >, P1StokesFunction< real_t > >
{
 public:
   P1EpsilonStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , A_uu( storage, minLevel, maxLevel )
   , A_uv( storage, minLevel, maxLevel )
   , A_vu( storage, minLevel, maxLevel )
   , A_vv( storage, minLevel, maxLevel )
   , div_x( storage, minLevel, maxLevel )
   , div_y( storage, minLevel, maxLevel )
   , divT_x( storage, minLevel, maxLevel )
   , divT_y( storage, minLevel, maxLevel )
   , pspg( storage, minLevel, maxLevel )
   {}

   void apply( const P1StokesFunction< real_t >& src, const P1StokesFunction< real_t >& dst, size_t level, DoFType flag ) const
   {
      WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

      A_uu.apply( src.uvw.u, dst.uvw.u, level, flag, Replace );
      A_uv.apply( src.uvw.v, dst.uvw.u, level, flag, Add );
      divT_x.apply( src.p, dst.uvw.u, level, flag, Add );

      A_vu.apply( src.uvw.u, dst.uvw.v, level, flag, Replace );
      A_vv.apply( src.uvw.v, dst.uvw.v, level, flag, Add );
      divT_y.apply( src.p, dst.uvw.v, level, flag, Add );

      div_x.apply( src.uvw.u, dst.p, level, flag | DirichletBoundary, Replace );
      div_y.apply( src.uvw.v, dst.p, level, flag | DirichletBoundary, Add );
      pspg.apply( src.p, dst.p, level, flag | DirichletBoundary, Add );
   }

   P1ConstantEpsilonOperator_11 A_uu;
   P1ConstantEpsilonOperator_12 A_uv;
   P1ConstantEpsilonOperator_21 A_vu;
   P1ConstantEpsilonOperator_22 A_vv;
   P1DivxOperator               div_x;
   P1DivyOperator               div_y;
   P1DivTxOperator              divT_x;
   P1DivTyOperator              divT_y;
   P1PSPGOperator               pspg;
};

template <>
struct has_pspg_block< P1EpsilonStokesOperator >
{
   static const bool value = true;
};

template <>
struct tensor_variant< P1EpsilonStokesOperator >
{
   static const bool value = true;
};

} // namespace hyteg

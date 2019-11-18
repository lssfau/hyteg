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
#include "hyteg/p1functionspace/P1PolynomialBlendingOperator.hpp"

namespace hyteg {

class P1PolynomialBlendingStokesOperator : public Operator< P1StokesFunction< real_t >, P1StokesFunction< real_t > >
{
 public:
   P1PolynomialBlendingStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                       uint_t                                     minLevel,
                                       uint_t                                     maxLevel,
                                       uint_t                                     interpolationLevel )
   : Operator( storage, minLevel, maxLevel )
   , A_uu( storage, minLevel, maxLevel, interpolationLevel )
   , A_uv( storage, minLevel, maxLevel, interpolationLevel )
   , A_vu( storage, minLevel, maxLevel, interpolationLevel )
   , A_vv( storage, minLevel, maxLevel, interpolationLevel )
   , div_x( storage, minLevel, maxLevel, interpolationLevel )
   , div_y( storage, minLevel, maxLevel, interpolationLevel )
   , divT_x( storage, minLevel, maxLevel, interpolationLevel )
   , divT_y( storage, minLevel, maxLevel, interpolationLevel )
   , pspg( storage, minLevel, maxLevel )
   {}

   void interpolateStencils( uint_t polyDegree )
   {
      A_uu.interpolateStencils( polyDegree );
      A_uv.interpolateStencils( polyDegree );
      A_vu.interpolateStencils( polyDegree );
      A_vv.interpolateStencils( polyDegree );
      div_x.interpolateStencils( polyDegree );
      div_y.interpolateStencils( polyDegree );
      divT_x.interpolateStencils( polyDegree );
      divT_y.interpolateStencils( polyDegree );
   }

   void useDegree( uint_t degree )
   {
      A_uu.useDegree( degree );
      A_uv.useDegree( degree );
      A_vu.useDegree( degree );
      A_vv.useDegree( degree );
      div_x.useDegree( degree );
      div_y.useDegree( degree );
      divT_x.useDegree( degree );
      divT_y.useDegree( degree );
   }

   void apply( const P1StokesFunction< real_t >& src, const P1StokesFunction< real_t >& dst, size_t level, DoFType flag ) const
   {
      WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

      A_uu.apply( src.u, dst.u, level, flag, Replace );
      A_uv.apply( src.v, dst.u, level, flag, Add );
      divT_x.apply( src.p, dst.u, level, flag, Add );

      A_vu.apply( src.u, dst.v, level, flag, Replace );
      A_vv.apply( src.v, dst.v, level, flag, Add );
      divT_y.apply( src.p, dst.v, level, flag, Add );

      div_x.apply( src.u, dst.p, level, flag | DirichletBoundary, Replace );
      div_y.apply( src.v, dst.p, level, flag | DirichletBoundary, Add );
      pspg.apply( src.p, dst.p, level, flag | DirichletBoundary, Add );
   }

   P1PolynomialBlendingEpsilonOperator_11 A_uu;
   P1PolynomialBlendingEpsilonOperator_12 A_uv;
   P1PolynomialBlendingEpsilonOperator_21 A_vu;
   P1PolynomialBlendingEpsilonOperator_22 A_vv;
   P1PolynomialBlendingDivOperator_1      div_x;
   P1PolynomialBlendingDivOperator_2      div_y;
   P1PolynomialBlendingDivTOperator_1     divT_x;
   P1PolynomialBlendingDivTOperator_2     divT_y;
   P1PSPGOperator                         pspg;
};

template <>
struct has_pspg_block< P1PolynomialBlendingStokesOperator >
{
   static const bool value = true;
};

template <>
struct tensor_variant< P1PolynomialBlendingStokesOperator >
{
   static const bool value = true;
};

} // namespace hyteg

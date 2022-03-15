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
#include "hyteg/operators/ScalarToVectorOperator.hpp"
#include "hyteg/operators/VectorToScalarOperator.hpp"
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
   , div( storage, minLevel, maxLevel )
   , divT( storage, minLevel, maxLevel )
   , pspg( storage, minLevel, maxLevel )
   {}

   void apply( const P1StokesFunction< real_t >& src, const P1StokesFunction< real_t >& dst, size_t level, DoFType flag ) const
   {
      WALBERLA_CHECK( !src.uvw()[0].getStorage()->hasGlobalCells(), "P1EpsilonStokesOperator not implemented for 3D." );
      WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

      A_uu.apply( src.uvw()[0], dst.uvw()[0], level, flag, Replace );
      A_uv.apply( src.uvw()[1], dst.uvw()[0], level, flag, Add );
      A_vu.apply( src.uvw()[0], dst.uvw()[1], level, flag, Replace );
      A_vv.apply( src.uvw()[1], dst.uvw()[1], level, flag, Add );

      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag | DirichletBoundary, Replace );

      pspg.apply( src.p(), dst.p(), level, flag | DirichletBoundary, Add );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1StokesFunction< idx_t >&            src,
                  const P1StokesFunction< idx_t >&            dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      WALBERLA_CHECK( !src.uvw()[0].getStorage()->hasGlobalCells(), "P1EpsilonStokesOperator not implemented for 3D." );

      A_uu.toMatrix( mat, src.uvw()[0], dst.uvw()[0], level, flag );
      A_uv.toMatrix( mat, src.uvw()[1], dst.uvw()[0], level, flag );
      A_vu.toMatrix( mat, src.uvw()[0], dst.uvw()[1], level, flag );
      A_vv.toMatrix( mat, src.uvw()[1], dst.uvw()[1], level, flag );

      divT.toMatrix( mat, src.p(), dst.uvw(), level, flag );
      div.toMatrix( mat, src.uvw(), dst.p(), level, flag | DirichletBoundary );

      pspg.toMatrix( mat, src.p(), dst.p(), level, flag | DirichletBoundary );
   }

   P1ConstantEpsilonOperator_11 A_uu;
   P1ConstantEpsilonOperator_12 A_uv;
   P1ConstantEpsilonOperator_21 A_vu;
   P1ConstantEpsilonOperator_22 A_vv;
   P1ConstantDivOperator        div;
   P1ConstantDivTOperator       divT;
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

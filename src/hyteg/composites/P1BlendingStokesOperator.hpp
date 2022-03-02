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

#include "hyteg/composites/P1StokesBlockPreconditioner.hpp"
#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/StokesOperatorTraits.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"

namespace hyteg {

class P1BlendingStokesOperator : public Operator< P1StokesFunction< real_t >, P1StokesFunction< real_t > >
{
 public:
   P1BlendingStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
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

   typedef P1StokesBlockPreconditioner BlockPreconditioner_T;

   void apply( const P1StokesFunction< real_t >& src,
               const P1StokesFunction< real_t >& dst,
               const uint_t                      level,
               const DoFType                     flag ) const
   {
      WALBERLA_CHECK( !src.uvw()[0].getStorage()->hasGlobalCells(), "P1BlendingStokesOperator not implemented for 3D." );

      WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

      A_uu.apply( src.uvw()[0], dst.uvw()[0], level, flag, Replace );
      A_uv.apply( src.uvw()[1], dst.uvw()[0], level, flag, Add );
      divT_x.apply( src.p(), dst.uvw()[0], level, flag, Add );

      A_vu.apply( src.uvw()[0], dst.uvw()[1], level, flag, Replace );
      A_vv.apply( src.uvw()[1], dst.uvw()[1], level, flag, Add );
      divT_y.apply( src.p(), dst.uvw()[1], level, flag, Add );

      div_x.apply( src.uvw()[0], dst.p(), level, flag | DirichletBoundary, Replace );
      div_y.apply( src.uvw()[1], dst.p(), level, flag | DirichletBoundary, Add );
      pspg.apply( src.p(), dst.p(), level, flag | DirichletBoundary, Add );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1StokesFunction< idx_t >&            src,
                  const P1StokesFunction< idx_t >&            dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      WALBERLA_CHECK( !src.uvw()[0].getStorage()->hasGlobalCells(), "P1BlendingStokesOperator not implemented for 3D." );

      A_uu.toMatrix( mat, src.uvw()[0], dst.uvw()[0], level, flag );
      A_uv.toMatrix( mat, src.uvw()[1], dst.uvw()[0], level, flag );
      divT_x.toMatrix( mat, src.p(), dst.uvw()[0], level, flag );

      A_vu.toMatrix( mat, src.uvw()[0], dst.uvw()[1], level, flag );
      A_vv.toMatrix( mat, src.uvw()[1], dst.uvw()[1], level, flag );
      divT_y.toMatrix( mat, src.p(), dst.uvw()[1], level, flag );

      div_x.toMatrix( mat, src.uvw()[0], dst.p(), level, flag | DirichletBoundary );
      div_y.toMatrix( mat, src.uvw()[1], dst.p(), level, flag | DirichletBoundary );

      pspg.toMatrix( mat, src.p(), dst.p(), level, flag | DirichletBoundary );
   }

   P1BlendingEpsilonOperator_11 A_uu;
   P1BlendingEpsilonOperator_12 A_uv;
   P1BlendingEpsilonOperator_21 A_vu;
   P1BlendingEpsilonOperator_22 A_vv;
   P1BlendingDivOperator_1      div_x;
   P1BlendingDivOperator_2      div_y;
   P1BlendingDivTOperator_1     divT_x;
   P1BlendingDivTOperator_2     divT_y;
   P1PSPGOperator               pspg;
};

template <>
struct has_pspg_block< P1BlendingStokesOperator >
{
   static const bool value = true;
};

template <>
struct tensor_variant< P1BlendingStokesOperator >
{
   static const bool value = true;
};

} // namespace hyteg

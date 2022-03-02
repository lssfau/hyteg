/*
 * Copyright (c) 2017-2021 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/composites/P2P2StokesFunction.hpp"
#include "hyteg/composites/StokesOperatorTraits.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"

namespace hyteg {

class P2P2StabilizedStokesOperator : public Operator< P2P2StokesFunction< real_t >, P2P2StokesFunction< real_t > >
{
 public:
   typedef P2ConstantLaplaceOperator VelocityOperator_T;
   typedef P2ConstantLaplaceOperator PressureOperator_T;

   P2P2StabilizedStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , A( storage, minLevel, maxLevel )
   , div_x( storage, minLevel, maxLevel )
   , div_y( storage, minLevel, maxLevel )
   , div_z( storage, minLevel, maxLevel )
   , divT_x( storage, minLevel, maxLevel )
   , divT_y( storage, minLevel, maxLevel )
   , divT_z( storage, minLevel, maxLevel )
   , pspg( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void apply( const P2P2StokesFunction< real_t >& src,
               const P2P2StokesFunction< real_t >& dst,
               const size_t                        level,
               DoFType                             flag ) const
   {
      A.apply( src.uvw()[0], dst.uvw()[0], level, flag, Replace );
      divT_x.apply( src.p(), dst.uvw()[0], level, flag, Add );

      A.apply( src.uvw()[1], dst.uvw()[1], level, flag, Replace );
      divT_y.apply( src.p(), dst.uvw()[1], level, flag, Add );

      if ( hasGlobalCells_ )
      {
         A.apply( src.uvw()[2], dst.uvw()[2], level, flag, Replace );
         divT_z.apply( src.p(), dst.uvw()[2], level, flag, Add );
      }

      div_x.apply( src.uvw()[0], dst.p(), level, flag, Replace );
      div_y.apply( src.uvw()[1], dst.p(), level, flag, Add );

      if ( hasGlobalCells_ )
      {
         div_z.apply( src.uvw()[2], dst.p(), level, flag, Add );
      }

      pspg.apply( src.p(), dst.p(), level, flag, Add );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P2StokesFunction< idx_t >&          src,
                  const P2P2StokesFunction< idx_t >&          dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      A.toMatrix( mat, src.uvw()[0], dst.uvw()[0], level, flag );
      divT_x.toMatrix( mat, src.p(), dst.uvw()[0], level, flag );

      A.toMatrix( mat, src.uvw()[1], dst.uvw()[1], level, flag );
      divT_y.toMatrix( mat, src.p(), dst.uvw()[1], level, flag );

      if ( src.uvw()[0].getStorage()->hasGlobalCells() )
      {
         A.toMatrix( mat, src.uvw()[2], dst.uvw()[2], level, flag );
         divT_z.toMatrix( mat, src.p(), dst.uvw()[2], level, flag );
      }

      div_x.toMatrix( mat, src.uvw()[0], dst.p(), level, flag | DirichletBoundary );
      div_y.toMatrix( mat, src.uvw()[1], dst.p(), level, flag | DirichletBoundary );
      if ( src.uvw()[0].getStorage()->hasGlobalCells() )
      {
         div_z.toMatrix( mat, src.uvw()[2], dst.p(), level, flag | DirichletBoundary );
      }

      pspg.toMatrix( mat, src.p(), dst.p(), level, flag | DirichletBoundary );
   }

   P2ConstantLaplaceOperator A;
   P2ConstantDivxOperator    div_x;
   P2ConstantDivyOperator    div_y;
   P2ConstantDivzOperator    div_z;
   P2ConstantDivTxOperator   divT_x;
   P2ConstantDivTyOperator   divT_y;
   P2ConstantDivTzOperator   divT_z;
   P2ConstantPSPGOperator    pspg;
   bool                      hasGlobalCells_;
};

template <>
struct has_pspg_block< P2P2StabilizedStokesOperator >
{
   static const bool value = true;
};

} // namespace hyteg

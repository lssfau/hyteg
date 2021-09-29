/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1ScalarToP1VectorOperator.hpp"
#include "hyteg/p1functionspace/P1VectorToP1ScalarOperator.hpp"

namespace hyteg {

class P1StokesOperator : public Operator< P1StokesFunction< real_t >, P1StokesFunction< real_t > >
{
 public:
   typedef P1ConstantVectorLaplaceOperator VelocityBlockOperator_T;
   typedef P1ConstantLaplaceOperator       VelocityOperator_T;
   typedef P1ConstantLaplaceOperator       PressureOperator_T;
   typedef P1StokesBlockPreconditioner     BlockPreconditioner_T;

   P1StokesOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , A( storage, minLevel, maxLevel )
   , div_x( storage, minLevel, maxLevel )
   , div_y( storage, minLevel, maxLevel )
   , div_z( storage, minLevel, maxLevel )
   , div( storage, minLevel, maxLevel )
   , divT_x( storage, minLevel, maxLevel )
   , divT_y( storage, minLevel, maxLevel )
   , divT_z( storage, minLevel, maxLevel )
   , divT( storage, minLevel, maxLevel )
   , pspg( storage, minLevel, maxLevel )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void apply( const P1StokesFunction< real_t >& src,
               const P1StokesFunction< real_t >& dst,
               const size_t                      level,
               DoFType                           flag ) const
   {
      WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

      A.apply( src.uvw[0], dst.uvw[0], level, flag, Replace );
      divT_x.apply( src.p, dst.uvw[0], level, flag, Add );

      A.apply( src.uvw[1], dst.uvw[1], level, flag, Replace );
      divT_y.apply( src.p, dst.uvw[1], level, flag, Add );

      if ( hasGlobalCells_ )
      {
         A.apply( src.uvw[2], dst.uvw[2], level, flag, Replace );
         divT_z.apply( src.p, dst.uvw[2], level, flag, Add );
      }

      div_x.apply( src.uvw[0], dst.p, level, flag, Replace );
      div_y.apply( src.uvw[1], dst.p, level, flag, Add );

      if ( hasGlobalCells_ )
      {
         div_z.apply( src.uvw[2], dst.p, level, flag, Add );
      }

      pspg.apply( src.p, dst.p, level, flag, Add );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1StokesFunction< matIdx_t >&         src,
                  const P1StokesFunction< matIdx_t >&         dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      A.toMatrix( mat, src.uvw[0], dst.uvw[0], level, flag );
      divT_x.toMatrix( mat, src.p, dst.uvw[0], level, flag );

      A.toMatrix( mat, src.uvw[1], dst.uvw[1], level, flag );
      divT_y.toMatrix( mat, src.p, dst.uvw[1], level, flag );

      if ( src.uvw[0].getStorage()->hasGlobalCells() )
      {
         A.toMatrix( mat, src.uvw[2], dst.uvw[2], level, flag );
         divT_z.toMatrix( mat, src.p, dst.uvw[2], level, flag );
      }

      div_x.toMatrix( mat, src.uvw[0], dst.p, level, flag | DirichletBoundary );
      div_y.toMatrix( mat, src.uvw[1], dst.p, level, flag | DirichletBoundary );
      if ( src.uvw[0].getStorage()->hasGlobalCells() )
      {
         div_z.toMatrix( mat, src.uvw[2], dst.p, level, flag | DirichletBoundary );
      }

      pspg.toMatrix( mat, src.p, dst.p, level, flag | DirichletBoundary );
   }

   P1ConstantLaplaceOperator A;
   P1DivxOperator            div_x;
   P1DivyOperator            div_y;
   P1DivzOperator            div_z;
   P1DivTxOperator           divT_x;
   P1DivTyOperator           divT_y;
   P1DivTzOperator           divT_z;
   P1PSPGOperator            pspg;
   P1PSPGInvDiagOperator     pspg_inv_diag_;
   bool                      hasGlobalCells_;

   P1ConstantDivOperator  div;
   P1ConstantDivTOperator divT;
};

template <>
struct has_pspg_block< P1StokesOperator >
{
   static const bool value = true;
};

} // namespace hyteg

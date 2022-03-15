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
#include "hyteg/operators/ScalarToVectorOperator.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/operators/VectorToScalarOperator.hpp"

namespace hyteg {

class P1P1StokesOperator : public Operator< P1StokesFunction< real_t >, P1StokesFunction< real_t > >
{
 public:
   typedef P1ConstantVectorLaplaceOperator VelocityBlockOperator_T;
   typedef P1ConstantLaplaceOperator       VelocityOperator_T;
   typedef P1ConstantLaplaceOperator       PressureOperator_T;
   typedef P1StokesBlockPreconditioner     BlockPreconditioner_T;

   P1P1StokesOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , lapl( storage, minLevel, maxLevel )
   , div( storage, minLevel, maxLevel )
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

      lapl.apply( src.uvw(), dst.uvw(), level, flag, Replace );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag, Replace );

      pspg.apply( src.p(), dst.p(), level, flag, Add );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1StokesFunction< idx_t >&            src,
                  const P1StokesFunction< idx_t >&            dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      lapl.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      divT.toMatrix( mat, src.p(), dst.uvw(), level, flag );
      div.toMatrix( mat, src.uvw(), dst.p(), level, flag );

      pspg.toMatrix( mat, src.p(), dst.p(), level, flag | DirichletBoundary );
   }

   const P1ConstantLaplaceOperator& getA() const
   {
      auto ptr = lapl.getSubOperator( 0, 0 );
      return dynamic_cast< const P1ConstantLaplaceOperator& >( *ptr );
   }

   P1ConstantVectorLaplaceOperator lapl;
   P1ConstantDivOperator           div;
   P1ConstantDivTOperator          divT;
   P1PSPGOperator                  pspg;
   P1PSPGInvDiagOperator           pspg_inv_diag_;

   bool hasGlobalCells_;
};

template <>
struct has_pspg_block< P1P1StokesOperator >
{
   static const bool value = true;
};

} // namespace hyteg

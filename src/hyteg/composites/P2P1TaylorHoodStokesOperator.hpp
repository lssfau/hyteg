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

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesBlockPreconditioner.hpp"
#include "hyteg/mixedoperators/P1ToP2Operator.hpp"
#include "hyteg/mixedoperators/P2ToP1Operator.hpp"
#include "hyteg/operators/ScalarToVectorOperator.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/operators/VectorToScalarOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"

namespace hyteg {

class P2P1TaylorHoodStokesOperator : public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef P2ConstantVectorLaplaceOperator         VelocityBlockOperator_T;
   typedef P2ConstantLaplaceOperator               VelocityOperator_T;
   typedef P2P1TaylorHoodStokesBlockPreconditioner BlockPreconditioner_T;

   P2P1TaylorHoodStokesOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , Lapl( storage, minLevel, maxLevel )
   , div( storage, minLevel, maxLevel )
   , divT( storage, minLevel, maxLevel )
   , pspg_( storage, minLevel, maxLevel )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               const uint_t                            level,
               const DoFType                           flag ) const
   {
      Lapl.apply( src.uvw(), dst.uvw(), level, flag, Replace );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag, Replace );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P1TaylorHoodFunction< idx_t >&      src,
                  const P2P1TaylorHoodFunction< idx_t >&      dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      Lapl.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      divT.toMatrix( mat, src.p(), dst.uvw(), level, flag );
      div.toMatrix( mat, src.uvw(), dst.p(), level, flag );
   }

   const P2ConstantLaplaceOperator& getA() const {
     auto ptr = Lapl.getSubOperator( 0, 0 );
     return dynamic_cast< const P2ConstantLaplaceOperator& >( *ptr );
   }

   VelocityBlockOperator_T    Lapl;
   P2ToP1ConstantDivOperator  div;
   P1ToP2ConstantDivTOperator divT;

   // this operator is needed in the uzawa smoother
   P1PSPGOperator        pspg_;
   P1PSPGInvDiagOperator pspg_inv_diag_;
   bool                  hasGlobalCells_;
};

} // namespace hyteg

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
#include "hyteg/operators/ScalarToVectorOperator.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/operators/VectorToScalarOperator.hpp"

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
   , lapl( storage, minLevel, maxLevel )
   , div( storage, minLevel, maxLevel )
   , divT( storage, minLevel, maxLevel )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void computeAndStoreLocalElementMatrices()
   {
      auto scalarA = dynamic_cast< P2ElementwiseLaplaceOperator& >( *lapl.getSubOperator( 0, 0 ) );
      scalarA.computeAndStoreLocalElementMatrices();

      div.getSubOperator< 0 >().computeAndStoreLocalElementMatrices();
      div.getSubOperator< 1 >().computeAndStoreLocalElementMatrices();
      div.getSubOperator< 2 >().computeAndStoreLocalElementMatrices();

      divT.getSubOperator< 0 >().computeAndStoreLocalElementMatrices();
      divT.getSubOperator< 1 >().computeAndStoreLocalElementMatrices();
      divT.getSubOperator< 2 >().computeAndStoreLocalElementMatrices();

      div.getSubOperator< 0 >().computeAndStoreLocalElementMatrices();
      div.getSubOperator< 1 >().computeAndStoreLocalElementMatrices();
      div.getSubOperator< 2 >().computeAndStoreLocalElementMatrices();

      divT.getSubOperator< 0 >().computeAndStoreLocalElementMatrices();
      divT.getSubOperator< 1 >().computeAndStoreLocalElementMatrices();
      divT.getSubOperator< 2 >().computeAndStoreLocalElementMatrices();
   }

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               const uint_t                            level,
               const DoFType                           flag ) const
   {
      WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

      lapl.apply( src.uvw(), dst.uvw(), level, flag, Replace );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag, Replace );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P1TaylorHoodFunction< idx_t >&      src,
                  const P2P1TaylorHoodFunction< idx_t >&      dst,
                  uint_t                                      level,
                  DoFType                                     flag ) const
   {
      lapl.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      divT.toMatrix( mat, src.p(), dst.uvw(), level, flag );
      div.toMatrix( mat, src.uvw(), dst.p(), level, flag );
   }

   const P2ElementwiseLaplaceOperator& getA() const
   {
      auto ptr = lapl.getSubOperator( 0, 0 );
      return dynamic_cast< const P2ElementwiseLaplaceOperator& >( *ptr );
   }

   P2ElementwiseVectorLaplaceOperator lapl;
   P2ToP1ElementwiseDivOperator       div;
   P1ToP2ElementwiseDivTOperator      divT;

   /// this operator is need in the uzawa smoother
   P1PSPGInvDiagOperator pspg_inv_diag_;

   bool hasGlobalCells_;
};

} // namespace hyteg

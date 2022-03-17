/*
 * Copyright (c) 2017-2022 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "hyteg/elementwiseoperators/P1ToP2ElementwiseOperator.hpp"
#include "hyteg/mixedoperators/MixedDummyOperators.hpp"
#include "hyteg/mixedoperators/P1ToP2Operator.hpp"
#include "hyteg/mixedoperators/P1ToP2SurrogateOperator.hpp"
#include "hyteg/mixedoperators/P1ToP2VariableOperator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"

namespace hyteg {

using walberla::real_t;

template < template < typename > class sKind_t, template < typename > class vKind_t, class operX_t, class operY_t, class operZ_t >
class ScalarToVectorOperator : public Operator< sKind_t< real_t >, vKind_t< real_t > >
{
 public:
   template < typename... SpecialCtorArgs >
   ScalarToVectorOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                           size_t                                     minLevel,
                           size_t                                     maxLevel,
                           SpecialCtorArgs... extraArgs )
   : Operator< sKind_t< real_t >, vKind_t< real_t > >( storage, minLevel, maxLevel )
   , operX( storage, minLevel, maxLevel, extraArgs... )
   , operY( storage, minLevel, maxLevel, extraArgs... )
   , operZ( storage, minLevel, maxLevel, extraArgs... )
   {}

   void apply( const sKind_t< real_t >& src,
               const vKind_t< real_t >& dst,
               size_t                   level,
               DoFType                  flag,
               UpdateType               updateType = Replace ) const
   {
      operX.apply( src, dst[0], level, flag, updateType );
      operY.apply( src, dst[1], level, flag, updateType );
      if ( dst.getDimension() == 3 )
         operZ.apply( src, dst[2], level, flag, updateType );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const sKind_t< idx_t >&                     src,
                  const vKind_t< idx_t >&                     dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      operX.toMatrix( mat, src, dst[0], level, flag );
      operY.toMatrix( mat, src, dst[1], level, flag );
      if ( dst.getDimension() == 3 )
         operZ.toMatrix( mat, src, dst[2], level, flag );
   }

   template < uint_t idx >
   const auto& getSubOperator() const
   {
      WALBERLA_ASSERT( idx < 3 );
      if constexpr ( idx == 0 )
      {
         return operX;
      }
      else if constexpr ( idx == 1 )
      {
         return operY;
      }
      else
      {
         return operZ;
      }
   }

   template < uint_t idx >
   auto& getSubOperator()
   {
      WALBERLA_ASSERT( idx < 3 );
      if constexpr ( idx == 0 )
      {
         return operX;
      }
      else if constexpr ( idx == 1 )
      {
         return operY;
      }
      else
      {
         return operZ;
      }
   }

 private:
   operX_t operX;
   operY_t operY;
   operZ_t operZ;

}; // namespace hyteg

// Some operators we might use more often than others
typedef ScalarToVectorOperator< P1Function,
                                P1VectorFunction,
                                P1BlendingDivTOperator_1,
                                P1BlendingDivTOperator_2,
                                P1BlendingDivTOperator_1 > // no 3D implementation exists
    P1ToP1BlendingDivTOperator;

typedef ScalarToVectorOperator< P1Function,
                                P2VectorFunction,
                                P1ToP2ConstantDivTxOperator,
                                P1ToP2ConstantDivTyOperator,
                                P1ToP2ConstantDivTzOperator >
    P1ToP2ConstantDivTOperator;

typedef ScalarToVectorOperator< P1Function,
                                P2VectorFunction,
                                P1ToP2BlendingDivTxOperator,
                                P1ToP2BlendingDivTyOperator,
                                P1ToP2DummyOperator >
    P1ToP2VariableDivTOperator;

typedef ScalarToVectorOperator< P1Function,
                                P2VectorFunction,
                                P1ToP2ElementwiseBlendingDivTxOperator,
                                P1ToP2ElementwiseBlendingDivTyOperator,
                                P1ToP2ElementwiseBlendingDivTzOperator >
    P1ToP2ElementwiseBlendingDivTOperator;

typedef ScalarToVectorOperator< P1Function,
                                P1VectorFunction,
                                P1ElementwiseDivTXOperator,
                                P1ElementwiseDivTYOperator,
                                P1ElementwiseDivTZOperator >
    P1ToP1ElementwiseDivTOperator;

typedef ScalarToVectorOperator< P1Function,
                                P2VectorFunction,
                                P1ToP2ElementwiseDivTxOperator,
                                P1ToP2ElementwiseDivTyOperator,
                                P1ToP2ElementwiseDivTzOperator >
    P1ToP2ElementwiseDivTOperator;

typedef ScalarToVectorOperator< P1Function, P1VectorFunction, P1DivTxOperator, P1DivTyOperator, P1DivTzOperator >
    P1ConstantDivTOperator;

typedef ScalarToVectorOperator< P2Function,
                                P2VectorFunction,
                                P2ConstantDivTxOperator,
                                P2ConstantDivTyOperator,
                                P2ConstantDivTzOperator >
    P2ConstantDivTOperator;

typedef ScalarToVectorOperator< P1Function,
                                P2VectorFunction,
                                P1ToP2SurrogateDivTxOperator,
                                P1ToP2SurrogateDivTyOperator,
                                P1ToP2SurrogateDivTzOperator >
    P1ToP2SurrogateDivTOperator;

} // namespace hyteg

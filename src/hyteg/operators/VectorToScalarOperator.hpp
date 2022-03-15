/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "hyteg/elementwiseoperators/P2ToP1ElementwiseOperator.hpp"
#include "hyteg/mixedoperators/MixedDummyOperators.hpp"
#include "hyteg/mixedoperators/P2ToP1Operator.hpp"
#include "hyteg/mixedoperators/P2ToP1SurrogateOperator.hpp"
#include "hyteg/mixedoperators/P2ToP1VariableOperator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"

namespace hyteg {

using walberla::real_t;

template < template < typename > class vKind_t,
           template < typename >
           class sKind_t,
           typename operX_t,
           typename operY_t,
           typename operZ_t >
class VectorToScalarOperator : public Operator< vKind_t< real_t >, sKind_t< real_t > >
{
 public:
   template < typename... SpecialCtorArgs >
   VectorToScalarOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                           size_t                                     minLevel,
                           size_t                                     maxLevel,
                           SpecialCtorArgs... extraArgs )
   : Operator< vKind_t< real_t >, sKind_t< real_t > >( storage, minLevel, maxLevel )
   , operX( storage, minLevel, maxLevel, extraArgs... )
   , operY( storage, minLevel, maxLevel, extraArgs... )
   , operZ( storage, minLevel, maxLevel, extraArgs... )
   {}

   void apply( const vKind_t< real_t >& src,
               const sKind_t< real_t >& dst,
               size_t                   level,
               DoFType                  flag,
               UpdateType               updateType = Replace ) const
   {
      std::array< UpdateType, 3 > ut = {updateType, Add, Add};
      operX.apply( src[0], dst, level, flag, ut[0] );
      operY.apply( src[1], dst, level, flag, ut[1] );
      if ( src.getDimension() == 3 )
         operZ.apply( src[2], dst, level, flag, ut[2] );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const vKind_t< idx_t >&                     src,
                  const sKind_t< idx_t >&                     dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      operX.toMatrix( mat, src[0], dst, level, flag );
      operY.toMatrix( mat, src[1], dst, level, flag );
      if ( src.getDimension() == 3 )
         operZ.toMatrix( mat, src[2], dst, level, flag );
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
};

// Some operators we might use more often than others
typedef VectorToScalarOperator< P1VectorFunction,
                                P1Function,
                                P1BlendingDivOperator_1,
                                P1BlendingDivOperator_2,
                                P1BlendingDivOperator_1 > // no 3D implementation exists!
    P1ToP1BlendingDivOperator;

typedef VectorToScalarOperator< P2VectorFunction,
                                P1Function,
                                P2ToP1ConstantDivxOperator,
                                P2ToP1ConstantDivyOperator,
                                P2ToP1ConstantDivzOperator >
    P2ToP1ConstantDivOperator;

typedef VectorToScalarOperator< P2VectorFunction,
                                P1Function,
                                P2ToP1BlendingDivxOperator,
                                P2ToP1BlendingDivyOperator,
                                P2ToP1DummyOperator >
    P2ToP1VariableDivOperator;

typedef VectorToScalarOperator< P2VectorFunction,
                                P1Function,
                                P2ToP1ElementwiseBlendingDivxOperator,
                                P2ToP1ElementwiseBlendingDivyOperator,
                                P2ToP1ElementwiseBlendingDivzOperator >
    P2ToP1ElementwiseBlendingDivOperator;

typedef VectorToScalarOperator< P1VectorFunction,
                                P1Function,
                                P1ElementwiseDivXOperator,
                                P1ElementwiseDivYOperator,
                                P1ElementwiseDivZOperator >
    P1ToP1ElementwiseDivOperator;

typedef VectorToScalarOperator< P2VectorFunction,
                                P1Function,
                                P2ToP1ElementwiseDivxOperator,
                                P2ToP1ElementwiseDivyOperator,
                                P2ToP1ElementwiseDivzOperator >
    P2ToP1ElementwiseDivOperator;

typedef VectorToScalarOperator< P1VectorFunction, P1Function, P1DivxOperator, P1DivyOperator, P1DivzOperator >
    P1ConstantDivOperator;

typedef VectorToScalarOperator< P2VectorFunction,
                                P2Function,
                                P2ConstantDivxOperator,
                                P2ConstantDivyOperator,
                                P2ConstantDivzOperator >
    P2ConstantDivOperator;

typedef VectorToScalarOperator< P2VectorFunction,
                                P1Function,
                                P2ToP1SurrogateDivxOperator,
                                P2ToP1SurrogateDivyOperator,
                                P2ToP1SurrogateDivzOperator >
    P2ToP1SurrogateDivOperator;

} // namespace hyteg

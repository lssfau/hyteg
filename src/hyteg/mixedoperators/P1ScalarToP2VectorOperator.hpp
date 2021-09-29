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

#include "hyteg/elementwiseoperators/P1ToP2ElementwiseOperator.hpp"
#include "hyteg/mixedoperators/MixedDummyOperators.hpp"
#include "hyteg/mixedoperators/P1ToP2Operator.hpp"
#include "hyteg/mixedoperators/P1ToP2VariableOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"

namespace hyteg {

using walberla::real_t;

template < class operX_t, class operY_t, class operZ_t >
class P1ScalarToP2VectorOperator : public Operator< P1Function< real_t >, P2VectorFunction< real_t > >
{
 public:
   P1ScalarToP2VectorOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , operX( storage, minLevel, maxLevel )
   , operY( storage, minLevel, maxLevel )
   , operZ( storage, minLevel, maxLevel )
   {}

   void apply( const P1Function< real_t >&       src,
               const P2VectorFunction< real_t >& dst,
               size_t                            level,
               DoFType                           flag,
               UpdateType                        updateType = Replace ) const
   {
      operX.apply( src, dst[0], level, flag, updateType );
      operY.apply( src, dst[1], level, flag, updateType );
      if ( dst.getDimension() == 3 )
         operZ.apply( src, dst[2], level, flag, updateType );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1Function< matIdx_t >&               src,
                  const P2VectorFunction< matIdx_t >&         dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      operX.toMatrix( mat, src, dst[0], level, flag );
      operY.toMatrix( mat, src, dst[1], level, flag );
      if ( dst.getDimension() == 3 )
         operZ.toMatrix( mat, src, dst[2], level, flag );
   }

 private:
   operX_t operX;
   operY_t operY;
   operZ_t operZ;
};

// Some operators we might use more often than others
typedef P1ScalarToP2VectorOperator< P1ToP2ConstantDivTxOperator, P1ToP2ConstantDivTyOperator, P1ToP2ConstantDivTzOperator >
    P1ToP2ConstantDivTOperator;

typedef P1ScalarToP2VectorOperator< P1ToP2BlendingDivTxOperator, P1ToP2BlendingDivTyOperator, P1ToP2DummyOperator >
    P1ToP2VariableDivTOperator;

typedef P1ScalarToP2VectorOperator< P1ToP2ElementwiseBlendingDivTxOperator,
                                    P1ToP2ElementwiseBlendingDivTyOperator,
                                    P1ToP2ElementwiseBlendingDivTzOperator >
    P1ToP2ElementwiseBlendingDivTOperator;

} // namespace hyteg

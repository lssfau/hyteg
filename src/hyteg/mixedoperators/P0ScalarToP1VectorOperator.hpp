/*
 * Copyright (c) 2017-2022 Nils Kohl.
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

#include "hyteg/mixedoperators/MixedDummyOperators.hpp"
#include "hyteg/mixedoperators/P0ToP1Operator.hpp"
#include "hyteg/mixedoperators/P1ToP0Operator.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"

namespace hyteg {

using walberla::real_t;

template < class operX_t, class operY_t, class operZ_t >
class P0ScalarToP1VectorOperator : public Operator< P0Function< real_t >, P1VectorFunction< real_t > >
{
 public:
   P0ScalarToP1VectorOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , operX( storage, minLevel, maxLevel )
   , operY( storage, minLevel, maxLevel )
   , operZ( storage, minLevel, maxLevel )
   {}

   void apply( const P0Function< real_t >&       src,
               const P1VectorFunction< real_t >& dst,
               size_t                            level,
               DoFType                           flag,
               UpdateType                        updateType = Replace ) const
   {
      operX.apply( src, dst[0], level, flag, updateType );
      operY.apply( src, dst[1], level, flag, updateType );
      if ( src.getDimension() == 3 )
         operZ.apply( src, dst[2], level, flag, updateType );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P0Function< idx_t >&                  src,
                  const P1VectorFunction< idx_t >&            dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      operX.toMatrix( mat, src, dst[0], level, flag );
      operY.toMatrix( mat, src, dst[1], level, flag );
      if ( src.getDimension() == 3 )
         operZ.toMatrix( mat, src, dst[2], level, flag );
   }

 private:
   operX_t operX;
   operY_t operY;
   operZ_t operZ;
};

//// Some operators we might use more often than others
typedef P0ScalarToP1VectorOperator< P0ToP1ConstantDivTxOperator, P0ToP1ConstantDivTyOperator, P0ToP1ConstantDivTzOperator >
    P0ToP1ConstantDivTOperator;

typedef P0ScalarToP1VectorOperator< P0ToP1ConstantP1EDGVectorLaplaceXCouplingOperator,
                                    P0ToP1ConstantP1EDGVectorLaplaceYCouplingOperator,
                                    P0ToP1ConstantP1EDGVectorLaplaceZCouplingOperator >
    P0ToP1ConstantP1EDGVectorLaplaceCouplingOperator;


typedef P0ScalarToP1VectorOperator< P0ToP1ConstantP1EDGEpsilonXCouplingOperator,
                                    P0ToP1ConstantP1EDGEpsilonYCouplingOperator,
                                    P0ToP1ConstantP1EDGEpsilonZCouplingOperator >  P0ToP1ConstantP1EDGEpsilonCouplingOperator;


typedef P0ScalarToP1VectorOperator< P0ToP1ConstantP1EDGVectorMassXCouplingOperator,
                                    P0ToP1ConstantP1EDGVectorMassYCouplingOperator,
                                    P0ToP1ConstantP1EDGVectorMassZCouplingOperator >
    P0ToP1ConstantP1EDGVectorMassCouplingOperator;

} // namespace hyteg

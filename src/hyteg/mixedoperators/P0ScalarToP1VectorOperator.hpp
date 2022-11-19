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
   typedef std::tuple< std::shared_ptr< typename operX_t::FormType >,
                       std::shared_ptr< typename operY_t::FormType >,
                       std::shared_ptr< typename operZ_t::FormType > >
       FormTuple;

   P0ScalarToP1VectorOperator(
       const std::shared_ptr< PrimitiveStorage >& storage,
       size_t                                     minLevel,
       size_t                                     maxLevel,
       // hand forms over to xyz operators: enable forms with variable coefficients
       FormTuple forms = std::make_tuple( std::make_shared< typename operX_t::FormType >(),
                                                        std::make_shared< typename operY_t::FormType >(),
                                                        std::make_shared< typename operZ_t::FormType >() ) )
   : Operator( storage, minLevel, maxLevel )
   , operX( storage, minLevel, maxLevel,  std::get<0>(forms)  )
   , operY( storage, minLevel, maxLevel,  std::get<1>(forms)  )
   , operZ( storage, minLevel, maxLevel,  std::get<2>(forms)  )
   {}

   void apply( const P0Function< real_t >&       src,
               const P1VectorFunction< real_t >& dst,
               size_t                            level,
               DoFType                           flag,
               UpdateType                        updateType = Replace ) const
   {
      operX.apply( src, dst[0], level, flag, updateType );
      operY.apply( src, dst[1], level, flag, updateType );
      if ( src.getDimension() == 3 ) {
         operZ.apply( src, dst[2], level, flag, updateType );
      }
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
      {
         operZ.toMatrix( mat, src, dst[2], level, flag );
      }
   }

 private:
   operX_t operX;
   operY_t operY;
   operZ_t operZ;
};

//// Some operators we might use more often than others
typedef P0ScalarToP1VectorOperator< P0ToP1ConstantDivTxOperator, P0ToP1ConstantDivTyOperator, P0ToP1ConstantDivTzOperator >
    P0ToP1ConstantDivTOperator;

typedef P0ScalarToP1VectorOperator< EGVectorLaplaceP0ToP1Coupling_X,
                                    EGVectorLaplaceP0ToP1Coupling_Y,
                                    EGVectorLaplaceP0ToP1Coupling_Z >
    EGVectorLaplaceP0ToP1Coupling;
typedef P0ScalarToP1VectorOperator< EGNIPGVectorLaplaceP0ToP1Coupling_X,
            EGNIPGVectorLaplaceP0ToP1Coupling_Y,
            EGNIPGVectorLaplaceP0ToP1Coupling_Z >
            EGNIPGVectorLaplaceP0ToP1Coupling;
typedef P0ScalarToP1VectorOperator< EGConstantEpsilonP0ToP1Coupling_X,
                                    EGConstantEpsilonP0ToP1Coupling_Y,
                                    EGConstantEpsilonP0ToP1Coupling_Z >
    EGConstantEpsilonP0ToP1Coupling;

typedef P0ScalarToP1VectorOperator< EGEpsilonP0ToP1Coupling_X, EGEpsilonP0ToP1Coupling_Y, EGEpsilonP0ToP1Coupling_Z >
    EGEpsilonP0ToP1Coupling;

typedef P0ScalarToP1VectorOperator< EGMassP0ToP1Coupling_X,
                                    EGMassP0ToP1Coupling_Y, EGMassP0ToP1Coupling_Z >
    EGMassP0toP1Coupling;

} // namespace hyteg

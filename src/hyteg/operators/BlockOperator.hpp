/*
 * Copyright (c) 2021 Marcus Mohr.
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
#include "hyteg/operators/Operator.hpp"
#include "hyteg/operators/OperatorWrapper.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"

namespace hyteg {

using walberla::real_t;

template < class srcBlockFunc_t, class dstBlockFunc_t >
class BlockOperator : public Operator< srcBlockFunc_t, dstBlockFunc_t >
{
 public:
   // temporary, need to add corresponding FunctionTrait?
   typedef typename srcBlockFunc_t::ValueType value_t;

   BlockOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                  size_t                                     minLevel,
                  size_t                                     maxLevel,
                  uint_t                                     nRows,
                  uint_t                                     nCols )
   : Operator< srcBlockFunc_t, dstBlockFunc_t >( storage, minLevel, maxLevel )
   {
      // setup internal 2D array for sub-operators
      nRows_ = nRows;
      nCols_ = nCols;
      subOper_.clear();
      subOper_.resize( nRows_, std::vector< std::shared_ptr< GenericOperator< value_t > > >( nCols_ ) );
   };

   /// Potentially we want to distinguish flag and updateType between components by passing a vector?
   void apply( const srcBlockFunc_t& src,
               const dstBlockFunc_t& dst,
               size_t                level,
               DoFType               flag,
               UpdateType            updateType = Replace ) const
   {
      WALBERLA_ASSERT_EQUAL( nCols_, src.getNumberOfBlocks() );
      WALBERLA_ASSERT_EQUAL( nRows_, dst.getNumberOfBlocks() );

      for ( uint_t i = 0; i < nRows_; i++ )
      {
         UpdateType upType = updateType;
         for ( uint_t j = 0; j < nCols_; j++ )
         {
            if ( subOper_[i][j] != nullptr )
            {
               WALBERLA_LOG_INFO_ON_ROOT( " -> applying sub-operator (" << i << ", " << j << ")" );
               subOper_[i][j]->apply( src[j], dst[i], level, flag, upType );
               upType = Add;
            }
         }
      }
   };

   void setSubOperator( uint_t i, uint_t j, std::shared_ptr< GenericOperator< value_t > > subOp )
   {
      WALBERLA_ASSERT_LESS( i, nRows_ );
      WALBERLA_ASSERT_LESS( j, nCols_ );
      subOper_[i][j] = subOp;
   }

   template < typename OperatorType, typename... ConstructorArguments >
   void createSubOperator( uint_t                                     i,
                           uint_t                                     j,
                           const std::shared_ptr< PrimitiveStorage >& storage,
                           uint_t                                     minLevel,
                           uint_t                                     maxLevel,
                           ConstructorArguments... args )
   {
      const auto op = std::make_shared< OperatorWrapper< OperatorType > >( storage, minLevel, maxLevel, args... );
      setSubOperator(i, j, op);
   }

   const std::shared_ptr< GenericOperator< value_t > > getSubOperator( uint_t i, uint_t j )
   {
      WALBERLA_ASSERT_LESS( i, nRows_ );
      WALBERLA_ASSERT_LESS( j, nCols_ );
      return subOper_[i][j];
   }

 protected:
   std::vector< std::vector< std::shared_ptr< GenericOperator< value_t > > > > subOper_;
   uint_t                                                                      nRows_;
   uint_t                                                                      nCols_;
};

} // namespace hyteg

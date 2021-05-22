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

#include <type_traits>

#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/operators/OperatorWrapper.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/solvers/Smoothables.hpp"

namespace hyteg {

using walberla::real_t;

template < class srcBlockFunc_t, class dstBlockFunc_t >
class BlockOperator : public Operator< srcBlockFunc_t, dstBlockFunc_t >,
                      public GSSmoothable< srcBlockFunc_t >,
                      public GSSmoothable< GenericFunction< typename srcBlockFunc_t::ValueType > >
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
   , tmp_( nullptr )
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

   const std::shared_ptr< GenericOperator< value_t > > getSubOperator( uint_t i, uint_t j )
   {
      WALBERLA_ASSERT_LESS( i, nRows_ );
      WALBERLA_ASSERT_LESS( j, nCols_ );
      return subOper_[i][j];
   }

   void smooth_gs( const srcBlockFunc_t& x, const srcBlockFunc_t& b, size_t level, DoFType flag ) const override
   {
      // check that both the rows and columns of our operator belong to the same type:
      if ( !std::is_same< srcBlockFunc_t, dstBlockFunc_t >::value )
         throw std::runtime_error( "For GaussSeidel src and dst functions must coincide" );

      // importing methods from our parent class:
      using ParentType = Operator< srcBlockFunc_t, dstBlockFunc_t >;

      // for a block Gauss-Seidel we need an intermediate vector
      // otherwise we would have to modify the rhs directly
      if ( !tmp_ )
         // Question: Can functions _always_ be constructed like that?
         tmp_ = std::make_unique< srcBlockFunc_t >(
             "tmp1", ParentType::getStorage(), ParentType::getMinLevel(), ParentType::getMaxLevel() );

      // assign our rhs
      tmp_->interpolate( 0, level, All );

      for ( uint_t rowIdx = 0; rowIdx < nRows_; rowIdx += 1 )
      {
         // calculate off-diagonal contributions: tmp1_[r] = sum_{c != r} A_rc x_c
         for ( uint_t colIdx = 0; colIdx < nCols_; colIdx += 1 )
         {
            // we add everything except the diagonal
            if ( colIdx == rowIdx )
               continue;

            subOper_[rowIdx][colIdx]->apply( x[colIdx], (*tmp_)[rowIdx], level, flag, Add );
         }
         // calculate new rhs: tmp1_[r] = b_r - sum_{c != r} A_rc x_c
         (*tmp_)[rowIdx].assign( { +1, -1 }, { b[rowIdx], (*tmp_)[rowIdx] }, level, flag );

         // diagonal part:
         auto D_rr = subOper_[rowIdx][rowIdx];
         if ( auto* op_diag = dynamic_cast< GSSmoothable< GenericFunction< value_t > >* >( D_rr.get() ) )
         {
            op_diag->smooth_gs( x[rowIdx], (*tmp_)[rowIdx], level, flag );
         }
      }
   }

   void smooth_gs( const GenericFunction< value_t >& dst,
                   const GenericFunction< value_t >& rhs,
                   size_t                            level,
                   DoFType                           flag ) const override
   {
      const auto& dst_unwrapped = dst.template unwrap< srcBlockFunc_t >();
      const auto& rhs_unwrapped = rhs.template unwrap< dstBlockFunc_t >();

      if ( !std::is_same< srcBlockFunc_t, dstBlockFunc_t >::value )
         smooth_gs( dst_unwrapped, rhs_unwrapped, level, flag );
      else
         throw std::runtime_error( "For GaussSeidel src and dst functions must coincide" );
   }

 protected:
   std::vector< std::vector< std::shared_ptr< GenericOperator< value_t > > > > subOper_;
   uint_t                                                                      nRows_;
   uint_t                                                                      nCols_;

   mutable std::unique_ptr< srcBlockFunc_t > tmp_;
};

} // namespace hyteg

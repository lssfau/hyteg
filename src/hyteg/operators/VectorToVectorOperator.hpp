/*
 * Copyright (c) 2017-2020 Marcus Mohr.
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
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"

namespace hyteg {

using walberla::real_t;

template < typename ValueType, template < typename > class SrcVecFuncKind, template < typename > class DstVecFuncKind >
class VectorToVectorOperator : public Operator< SrcVecFuncKind< ValueType >, DstVecFuncKind< ValueType > >
{
 public:
   typedef SrcVecFuncKind< ValueType >                                                                            SrcVecFuncType;
   typedef DstVecFuncKind< ValueType >                                                                            DstVecFuncType;
   typedef Operator< typename SrcVecFuncType::VectorComponentType, typename DstVecFuncType::VectorComponentType > scalarOpType;

   VectorToVectorOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator< SrcVecFuncType, DstVecFuncType >( storage, minLevel, maxLevel )
   {
      // deduce dimension from mesh (not flexible, think manifolds, boundary, interfaces)
      dim_ = storage->hasGlobalCells() ? 3 : 2;

      // setup internal 2D array for scalar sub-operators
      subOper_.clear();
      subOper_.resize( dim_, std::vector< std::shared_ptr< scalarOpType > >( dim_ ) );
   };

   /// Potentially we want to distinguish flag and updateType between components by passing a vector?
   void apply( const SrcVecFuncType& src,
               const DstVecFuncType& dst,
               size_t                level,
               DoFType               flag,
               UpdateType            updateType = Replace ) const
   {
      WALBERLA_ASSERT_EQUAL( dim_, src.getDimension() );
      WALBERLA_ASSERT_EQUAL( dim_, dst.getDimension() );

      for ( uint_t i = 0; i < dim_; i++ )
      {
         UpdateType upType = updateType;
         for ( uint_t j = 0; j < dim_; j++ )
         {
            if ( subOper_[i][j] != nullptr )
            {
               // WALBERLA_LOG_INFO_ON_ROOT( " -> applying sub-operator (" << i << ", " << j << ")" );
               subOper_[i][j]->apply( src[j], dst[i], level, flag, upType );
               upType = Add;
            }
         }
      }
   };

   void setSubOperator( uint_t i, uint_t j, std::shared_ptr< scalarOpType > subOp )
   {
      WALBERLA_ASSERT_LESS( i, dim_ );
      WALBERLA_ASSERT_LESS( j, dim_ );
      subOper_[i][j] = subOp;
   }

   const std::shared_ptr< scalarOpType > getSubOperator( uint_t i, uint_t j )
   {
      WALBERLA_ASSERT_LESS( i, dim_ );
      WALBERLA_ASSERT_LESS( j, dim_ );
      return subOper_[i][j];
   }

 protected:
   std::vector< std::vector< std::shared_ptr< scalarOpType > > > subOper_;
   uint_t                                                        dim_;
};

} // namespace hyteg

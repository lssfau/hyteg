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
#include "hyteg/solvers/Smoothables.hpp"

namespace hyteg {

using walberla::real_t;

template < typename ValueType, template < typename > class SrcVecFuncKind, template < typename > class DstVecFuncKind >
class VectorToVectorOperator : public Operator< SrcVecFuncKind< ValueType >, DstVecFuncKind< ValueType > >
{
 public:
   typedef SrcVecFuncKind< ValueType >                                                                            SrcVecFuncType;
   typedef DstVecFuncKind< ValueType >                                                                            DstVecFuncType;
   typedef Operator< typename SrcVecFuncType::VectorComponentType, typename DstVecFuncType::VectorComponentType > scalarOpType;

   // for compatibility with Operator class
   typedef SrcVecFuncKind< ValueType > srcType;
   typedef DstVecFuncKind< ValueType > dstType;

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

   const std::shared_ptr< scalarOpType > getSubOperator( uint_t i, uint_t j ) const
   {
      WALBERLA_ASSERT_LESS( i, dim_ );
      WALBERLA_ASSERT_LESS( j, dim_ );
      return subOper_[i][j];
   }

   std::shared_ptr< scalarOpType > getSubOperator( uint_t i, uint_t j )
   {
      WALBERLA_ASSERT_LESS( i, dim_ );
      WALBERLA_ASSERT_LESS( j, dim_ );
      return subOper_[i][j];
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const SrcVecFuncKind< idx_t >&              src,
                  const DstVecFuncKind< idx_t >&              dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      for ( uint_t i = 0; i < dim_; i++ )
      {
         for ( uint_t j = 0; j < dim_; j++ )
         {
            if ( subOper_[i][j] != nullptr )
            {
               subOper_[i][j]->toMatrix( mat, src[j], dst[i], level, flag );
            }
         }
      }
   }

   /// Trigger (re)computation of inverse diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeInverseDiagonalOperatorValues()
   {
      // operator must map between the same spaces
      bool consistent = std::is_same< SrcVecFuncType, DstVecFuncType >::value;
      WALBERLA_UNUSED( consistent );
      WALBERLA_ASSERT( consistent );

      using subType = typename SrcVecFuncType::VectorComponentType;

      for ( uint_t i = 0; i < dim_; i++ )
      {
         if ( auto* A_with_inv_diag = dynamic_cast< OperatorWithInverseDiagonal< subType >* >( subOper_[i][i].get() ) )
         {
            A_with_inv_diag->computeInverseDiagonalOperatorValues();
         }
         else
         {
            throw std::runtime_error(
                "VectorToVectorOperator::computeInverseDiagonalOperatorValues() requires sub-operators with the OperatorWithInverseDiagonal interface." );
         }
      }
   }

 protected:
   std::vector< std::vector< std::shared_ptr< scalarOpType > > > subOper_;
   uint_t                                                        dim_;

   std::shared_ptr< SrcVecFuncType > extractInverseDiagonal() const
   {
      // operator must map between the same spaces
      bool consistent = std::is_same< SrcVecFuncType, DstVecFuncType >::value;
      WALBERLA_UNUSED( consistent );
      WALBERLA_ASSERT( consistent );

      using subType = typename SrcVecFuncType::VectorComponentType;
      std::vector< std::shared_ptr< subType > > diags;

      for ( uint_t i = 0; i < dim_; i++ )
      {
         if ( const auto* A_with_inv_diag =
                  dynamic_cast< const OperatorWithInverseDiagonal< subType >* >( subOper_[i][i].get() ) )
         {
            diags.push_back( A_with_inv_diag->getInverseDiagonalValues() );
         }
         else
         {
            throw std::runtime_error(
                "Diagonal extraction from VectorToVectorOperator requires sub-operators with the OperatorWithInverseDiagonal interface." );
         }
      }

      return std::make_shared< SrcVecFuncType >( "diagonal operator part", diags );
   }
};

} // namespace hyteg

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

#include "hyteg/Operator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"

namespace hyteg {

using walberla::real_t;

// template < class FunctionType >
// class ZeroOperator : public Operator< FunctionType, FunctionType >
// {
//  public:
//    void apply( const FunctionType& src, const FunctionType& dst, size_t level, DoFType flag, UpdateType updateType = Replace )
//    {
//       if ( updateType == Add )
//          return;
//       dst.interpolate( FunctionType::ValueType( 0 ), level, flag );
//    }
//    const;
// }

template < class SrcVecFuncType, class DstVecFuncType >
class VectorToVectorOperator : public Operator< SrcVecFuncType, DstVecFuncType >
{
 public:
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
   void apply( const SrcVecFuncType& src, const DstVecFuncType& dst, size_t level, DoFType flag, UpdateType updateType = Replace ) const
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
               WALBERLA_LOG_INFO_ON_ROOT( " -> applying sub-operator (" << i << ", " << j << ")" );
               subOper_[i][j]->apply( src[j], dst[i], level, flag, upType );
               upType = Add;
            }
         }
      }
   };

   void setSubOperator( uint_t i, uint_t j, std::shared_ptr< scalarOpType> subOp ) {
     WALBERLA_ASSERT_LESS_EQUAL( i, dim_ );
     WALBERLA_ASSERT_LESS_EQUAL( j, dim_ );
     subOper_[i][j] = subOp;
   }

 protected:
   std::vector< std::vector< std::shared_ptr< scalarOpType > > > subOper_;
   uint_t                                                        dim_;
};

class P2ConstantVectorLaplaceOperator : public VectorToVectorOperator< P2VectorFunction< real_t >, P2VectorFunction< real_t > >
{
 public:
   P2ConstantVectorLaplaceOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : VectorToVectorOperator< P2VectorFunction< real_t >, P2VectorFunction< real_t > >( storage, minLevel, maxLevel )
   {
      std::shared_ptr< scalarOpType > zero( nullptr );
      std::shared_ptr< scalarOpType > lapl = std::make_shared< P2ConstantLaplaceOperator >( storage, minLevel, maxLevel );

      if ( dim_ == 3 )
      {
         subOper_[0][0] = lapl;
         subOper_[0][1] = zero;
         subOper_[0][2] = zero;

         subOper_[1][0] = zero;
         subOper_[1][1] = lapl;
         subOper_[1][2] = zero;

         subOper_[2][0] = zero;
         subOper_[2][1] = zero;
         subOper_[2][2] = lapl;
      }
      else
      {
         subOper_[0][0] = lapl;
         subOper_[0][1] = zero;

         subOper_[1][0] = zero;
         subOper_[1][1] = lapl;
      }
   };
};


class P1VectorOperator : public VectorToVectorOperator< P1VectorFunction_AltKind< real_t >, P1VectorFunction_AltKind< real_t > >
{
 public:
   P1VectorOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
     : VectorToVectorOperator< P1VectorFunction_AltKind< real_t >, P1VectorFunction_AltKind< real_t > >( storage, minLevel, maxLevel )
   {};
};


// Example using a "factory"
P1VectorOperator generateP1VectorLaplaceOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel ) {

  // setup empty operator
  P1VectorOperator vecOper( storage, minLevel, maxLevel );

  // select a flavour
  bool constant = true;

  using scalarOpType = Operator< P1Function<real_t>, P1Function<real_t> >;
  std::shared_ptr< scalarOpType > zero( nullptr );
  std::shared_ptr< scalarOpType > lapl;

  if( constant ) {
    lapl = std::make_shared< P1ConstantLaplaceOperator >( storage, minLevel, maxLevel );
  }
  else {
    lapl = std::make_shared< P1ElementwiseLaplaceOperator >( storage, minLevel, maxLevel );
  }

  if( storage->hasGlobalCells() ) {
    WALBERLA_LOG_INFO_ON_ROOT( "Creating 3D operator" );
    vecOper.setSubOperator( 0, 0, lapl );
    vecOper.setSubOperator( 0, 1, zero );
    vecOper.setSubOperator( 0, 2, zero );

    vecOper.setSubOperator( 1, 0, zero );
    vecOper.setSubOperator( 1, 1, lapl );
    vecOper.setSubOperator( 1, 2, zero );

    vecOper.setSubOperator( 2, 0, zero );
    vecOper.setSubOperator( 2, 1, zero );
    vecOper.setSubOperator( 2, 2, lapl );
  }
  else {
    WALBERLA_LOG_INFO_ON_ROOT( "Creating 2D operator" );
    vecOper.setSubOperator( 0, 0, lapl );
    vecOper.setSubOperator( 0, 1, zero );

    vecOper.setSubOperator( 1, 0, zero );
    vecOper.setSubOperator( 1, 1, lapl );
  }

  return vecOper;
}

} // namespace hyteg

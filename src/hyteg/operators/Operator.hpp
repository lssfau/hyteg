/*
 * Copyright (c) 2017-2019 Boerge Struempfel, Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include <core/DataTypes.h>
#include <core/timing/TimingTree.h>
#include <memory>

#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

namespace hyteg {

template < typename SourceFunction, typename DestinationFunction >
class Operator
{
 public:
   Operator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   {
      if ( storage->getTimingTree() )
      {
         enableTiming( storage->getTimingTree() );
      }
   }

   typedef SourceFunction      srcType;
   typedef DestinationFunction dstType;

   virtual ~Operator() = default;

   void enableTiming( const std::shared_ptr< walberla::WcTimingTree >& timingTree ) { timingTree_ = timingTree; }

   const std::shared_ptr< PrimitiveStorage > getStorage() const { return storage_; }

   uint_t getMinLevel() const { return minLevel_; }

   uint_t getMaxLevel() const { return maxLevel_; }

   // See the free function of the same name. Method might be revived, if we implement a base class for composite
   // operators
   //
   // uint_t getNumberOfGlobalDoFCouplings( uint_t level ) const
   // {
   //    uint_t nCouplings =
   //        indexing::countLocalDoFCouplings< typename srcType::Tag, typename dstType::Tag >( storage_, level );
   //     walberla::mpi::allReduceInplace( nCouplings, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
   //    return nCouplings;
   // };

   // Needed in VectorToVectorOperator
   virtual void apply( const SourceFunction&      src,
                       const DestinationFunction& dst,
                       size_t                     level,
                       DoFType                    flag,
                       UpdateType                 updateType = Replace ) const
   {
      WALBERLA_ABORT( "Problem with inheritance, this function should have been overridden in the derived class!" );
   };

   /// Implements the BLAS level 2 routine gemv:
   ///
   ///   y ← ɑAx + βy
   ///
   /// where ɑ, β are scalars, A is this operator, x is the "src" function, and y is the "dst" function.
   ///
   /// Note that the apply() method implements
   ///
   ///   y ← Ax + ɣy
   ///
   /// where ɣ ∈ { 0, 1 } for UpdateType Replace and Add respectively. This is just a gemv with ɑ = 1 and β = ɣ.
   virtual void gemv( const typename SourceFunction::valueType&      alpha,
                      const SourceFunction&                          src,
                      const typename DestinationFunction::valueType& beta,
                      const DestinationFunction&                     dst,
                      size_t                                         level,
                      DoFType                                        flag ) const
   {
      WALBERLA_ABORT( "gemv() not implemented in derived class!" );
   };

   virtual void toMatrix( const std::shared_ptr< SparseMatrixProxy >&             mat,
                          const typename srcType::template FunctionType< idx_t >& src,
                          const typename dstType::template FunctionType< idx_t >& dst,
                          size_t                                                  level,
                          DoFType                                                 flag ) const
   {
      WALBERLA_ABORT( "toMatrix() not implemented in derived class!" );
   }

 protected:
   const std::shared_ptr< PrimitiveStorage > storage_;
   const uint_t                              minLevel_;
   const uint_t                              maxLevel_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;

 protected:
   void startTiming( const std::string& timerString ) const
   {
      if ( timingTree_ )
      {
         timingTree_->start( "Operator " + FunctionTrait< SourceFunction >::getTypeName() + " to " +
                             FunctionTrait< DestinationFunction >::getTypeName() );
         timingTree_->start( timerString );
      }
   }

   void stopTiming( const std::string& timerString ) const
   {
      if ( timingTree_ )
      {
         timingTree_->stop( timerString );
         timingTree_->stop( "Operator " + FunctionTrait< SourceFunction >::getTypeName() + " to " +
                            FunctionTrait< DestinationFunction >::getTypeName() );
      }
   }
};

} // namespace hyteg

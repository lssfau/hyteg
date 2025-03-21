/*
 * Copyright (c) 2017-2023 Boerge Struempfel, Daniel Drzisga, Dominik Thoennes,
 * Nils Kohl, Daniel Bauer, Marcus Mohr.
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

#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "hyteg/communication/BufferedCommunication.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {

/// Base class for all HyTeG functions representing scalar fields or
/// vector fields which are not made up of individual scalar component functions.
template < typename FunctionType >
class Function
{
 public:
   typedef typename FunctionTrait< FunctionType >::Tag       Tag;
   typedef typename FunctionTrait< FunctionType >::ValueType valueType;

   Function( std::string name, const std::shared_ptr< PrimitiveStorage >& storage )
   : functionName_( std::move( name ) )
   , storage_( storage )
   , minLevel_( 0 )
   , maxLevel_( 0 )
   {}

   Function( std::string name, const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : functionName_( name )
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   {
      WALBERLA_CHECK_NOT_NULLPTR( storage_, "PrimitiveStorage pointer passed to Function is nullptr." );

      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         communicators_[level]         = std::make_shared< communication::BufferedCommunicator >( storage );
         additiveCommunicators_[level] = std::make_shared< communication::BufferedCommunicator >( storage );
      }
      if ( storage->getTimingTree() )
      {
         enableTiming( storage->getTimingTree() );
      }
   }

   virtual ~Function() = default;

   const std::string& getFunctionName() const { return functionName_; }

   /// Query function object for the dimension of the field it represents
   virtual uint_t getDimension() const = 0;

   /// Query function object for minimal level on which it defined
   uint_t getMinLevel() const { return minLevel_; }

   /// Query function object for maximal level on which it defined
   uint_t getMaxLevel() const { return maxLevel_; }

   /// Set all function DoFs to zero including the ones in the halos
   virtual void setToZero( const uint_t level ) const = 0;

   std::shared_ptr< PrimitiveStorage > getStorage() const { return storage_; }

   void enableTiming( const std::shared_ptr< walberla::WcTimingTree >& timingTree )
   {
      timingTree_ = timingTree;
      for ( auto& communicator : communicators_ )
      {
         communicator.second->enableTiming( timingTree_ );
      }
      for ( auto& communicator : additiveCommunicators_ )
      {
         communicator.second->enableTiming( timingTree_ );
      }
   }

 protected:
   std::string                         functionName_;
   std::shared_ptr< PrimitiveStorage > storage_;
   uint_t                              minLevel_;
   uint_t                              maxLevel_;

   std::map< uint_t, std::shared_ptr< communication::BufferedCommunicator > > communicators_;
   std::map< uint_t, std::shared_ptr< communication::BufferedCommunicator > > additiveCommunicators_;

   std::shared_ptr< walberla::WcTimingTree > timingTree_;

   void startTiming( const std::string& timerString ) const
   {
      if ( timingTree_ )
      {
         timingTree_->start( FunctionTrait< FunctionType >::getTypeName() );
         timingTree_->start( timerString );
      }
   }

   void stopTiming( const std::string& timerString ) const
   {
      if ( timingTree_ )
      {
         timingTree_->stop( timerString );
         timingTree_->stop( FunctionTrait< FunctionType >::getTypeName() );
      }
   }
};

} // namespace hyteg

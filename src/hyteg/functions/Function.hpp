/*
 * Copyright (c) 2017-2022 Boerge Struempfel, Daniel Drzisga, Dominik Thoennes, Nils Kohl, Daniel Bauer.
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
   , isDummy_( true )
   {}

   Function( std::string name, const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : functionName_( name )
   , storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , isDummy_( false )
   {
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         communicators_[level]         = std::make_shared< communication::BufferedCommunicator >( storage );
         additiveCommunicators_[level] = std::make_shared< communication::BufferedCommunicator >( storage );
      }
      if ( storage->getTimingTree() )
      {
         enableTiming( storage->getTimingTree() );
      }

      functionNames_.push_back( name );

      for ( uint_t i = minLevel; i <= maxLevel; ++i )
      {
         levelWiseFunctionCounter_[i]++;
      }
   }

   virtual ~Function() = default;

   const std::string& getFunctionName() const { return functionName_; }

   virtual uint_t getDimension() const = 0;

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

   bool isDummy() const { return isDummy_; }

   static uint_t                     getNumFunctions() { return functionNames_.size(); }
   static std::vector< std::string > getFunctionNames() { return functionNames_; }
   static std::map< uint_t, uint_t > getLevelWiseFunctionCounter() { return levelWiseFunctionCounter_; }

 protected:
   std::string                               functionName_;
   std::shared_ptr< PrimitiveStorage > storage_;
   uint_t                              minLevel_;
   uint_t                              maxLevel_;
   bool                                isDummy_;

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

 private:
   static std::vector< std::string > functionNames_;
   static std::map< uint_t, uint_t > levelWiseFunctionCounter_;
};

template < typename FunctionType >
std::vector< std::string > Function< FunctionType >::functionNames_ = {};

template < typename FunctionType >
std::map< uint_t, uint_t > Function< FunctionType >::levelWiseFunctionCounter_ = {};

} // namespace hyteg

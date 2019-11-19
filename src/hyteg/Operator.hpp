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
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/FunctionTraits.hpp"

#include <memory>

namespace hyteg {

template< typename SourceFunction, typename DestinationFunction >
class Operator
{
public:
  Operator(const std::shared_ptr<PrimitiveStorage> & storage, uint_t minLevel, uint_t maxLevel)
    : storage_(storage), minLevel_(minLevel), maxLevel_(maxLevel)
  {
  if(storage->getTimingTree()){
      enableTiming(storage->getTimingTree());
    }
  }

  typedef SourceFunction srcType;
  typedef DestinationFunction dstType;

  virtual ~Operator() = default;

  void enableTiming( const std::shared_ptr< walberla::WcTimingTree > & timingTree ) { timingTree_ = timingTree; }

  const std::shared_ptr< PrimitiveStorage > getStorage() const { return storage_; }

  uint_t getMinLevel() const { return minLevel_; }

  uint_t getMaxLevel() const { return maxLevel_; }

 protected:

  const std::shared_ptr< PrimitiveStorage > storage_;
  const uint_t minLevel_;
  const uint_t maxLevel_;

  std::shared_ptr< walberla::WcTimingTree > timingTree_;

 protected:

  void startTiming( const std::string & timerString ) const
  {
    if ( timingTree_ )
    {
      timingTree_->start( "Operator " + FunctionTrait< SourceFunction >::getTypeName() + " to " + FunctionTrait< DestinationFunction >::getTypeName() );
      timingTree_->start( timerString );
    }
  }

  void stopTiming ( const std::string & timerString ) const
  {
    if ( timingTree_ )
    {
      timingTree_->stop( timerString );
      timingTree_->stop( "Operator " + FunctionTrait< SourceFunction >::getTypeName() + " to " + FunctionTrait< DestinationFunction >::getTypeName() );
    }
  }
};

}

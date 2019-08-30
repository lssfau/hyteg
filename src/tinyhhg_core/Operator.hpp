
#pragma once

#include <core/DataTypes.h>
#include <core/timing/TimingTree.h>
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/FunctionTraits.hpp"

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

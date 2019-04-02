
#pragma once

#include <utility>
#include <functional>
#include <string>
#include <vector>

#include "core/mpi/Gather.h"

#include "tinyhhg_core/FunctionTraits.hpp"
#include "tinyhhg_core/Operator.hpp"
#include "tinyhhg_core/communication/BufferedCommunication.hpp"
#include "tinyhhg_core/types/flags.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

namespace hhg {


template< typename FunctionType >
class Function {
public:

  typedef typename FunctionTrait< FunctionType >::Tag Tag;
  typedef typename FunctionTrait< FunctionType >::ValueType valueType;

  Function( std::string name, const std::shared_ptr<PrimitiveStorage> & storage ) : functionName_(std::move(name)), storage_( storage ), minLevel_( 0 ), maxLevel_( 0 ), isDummy_( true ) {}

  Function(std::string name, const std::shared_ptr<PrimitiveStorage> & storage, uint_t minLevel, uint_t maxLevel)
      : functionName_(std::move(name))
      , storage_(storage)
      , minLevel_(minLevel)
      , maxLevel_(maxLevel)
      , isDummy_( false )
  {
    for ( uint_t level = minLevel; level <= maxLevel; level++ )
    {
      communicators_[ level ] = std::make_shared< communication::BufferedCommunicator >( storage );
      additiveCommunicators_[ level ] = std::make_shared< communication::BufferedCommunicator >( storage );
    }
    if(storage->getTimingTree())
    {
      enableTiming(storage->getTimingTree());
    }
    for( uint_t i = minLevel; i <= maxLevel; ++i )
    {
      functionCounter_[i]++;
    }

  }

  virtual ~Function() = default;

  const std::string &getFunctionName() const { return functionName_; }

  const std::shared_ptr< PrimitiveStorage > getStorage() const { WALBERLA_ASSERT( !!(storage_.lock()) ); return storage_.lock(); }

  uint_t getMinLevel() const { return minLevel_; }

  uint_t getMaxLevel() const { return maxLevel_; }

  void enableTiming( const std::shared_ptr< walberla::WcTimingTree > & timingTree )
  {
    timingTree_ = timingTree;
    for ( auto & communicator : communicators_ )
    {
      communicator.second->enableTiming( timingTree_ );
    }
    for ( auto & communicator : additiveCommunicators_ )
    {
      communicator.second->enableTiming( timingTree_ );
    }
  }

  bool isDummy() const { return isDummy_; }

  static std::map< uint_t, uint_t > getFunctionCounter() { return functionCounter_; }

protected:

  const std::string functionName_;
  const std::weak_ptr< PrimitiveStorage > storage_;
  const uint_t minLevel_;
  const uint_t maxLevel_;
  const bool   isDummy_;

  std::map< uint_t, std::shared_ptr< communication::BufferedCommunicator > > communicators_;
  std::map< uint_t, std::shared_ptr< communication::BufferedCommunicator > > additiveCommunicators_;

  std::shared_ptr< walberla::WcTimingTree > timingTree_;

  void startTiming( const std::string & timerString ) const
  {
    if ( timingTree_ )
    {
      timingTree_->start( FunctionTrait< FunctionType >::getTypeName() );
      timingTree_->start( timerString );
    }
  }

  void stopTiming ( const std::string & timerString ) const
  {
    if ( timingTree_ )
    {
      timingTree_->stop( timerString );
      timingTree_->stop( FunctionTrait< FunctionType >::getTypeName() );
    }
  }

private:
  static std::map< uint_t, uint_t > functionCounter_;

};

template < typename FunctionType >
std::map< uint_t, uint_t > Function< FunctionType >::functionCounter_ = {};

}



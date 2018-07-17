#pragma once

#include <tinyhhg_core/Operator.hpp>
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/types/flags.hpp"
#include "tinyhhg_core/communication/BufferedCommunication.hpp"
#include <core/mpi/Gather.h>

#include "tinyhhg_core/FunctionTraits.hpp"

#include <string>
#include <functional>
#include <vector>


namespace hhg {

template< typename FunctionType >
class Function {
public:

  typedef typename FunctionTrait< FunctionType >::ValueType ValueType;

  Function( std::string name, const std::shared_ptr<PrimitiveStorage> & storage ) : functionName_( name ), storage_( storage ), minLevel_( 0 ), maxLevel_( 0 ), isDummy_( true ) {}

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
    if(storage->getTimingTree()){
      enableTiming(storage->getTimingTree());
    }
  }

  virtual ~Function() {}

  inline uint_t enumerate( uint_t   level,
                           uint_t & num );


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

protected:

  virtual void
  enumerate_impl( uint_t level, uint_t& num ) = 0;

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
      timingTree_->start( "Function" );
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
      timingTree_->stop( "Function" );
    }
  }

};


template< typename FunctionType >
uint_t Function< FunctionType >::enumerate(size_t level, uint_t& num)
{
  startTiming( "Enumerate" );
  uint_t counter = 0;


  enumerate_impl(level,counter);

  std::vector<uint_t> dofs_per_rank  = walberla::mpi::allGather(counter);

  uint_t start = num;

  for(uint_t i = 0; i<walberla::MPIManager::instance()->rank();++i) {
    start += dofs_per_rank[i];
  }

  for(uint_t i = 0; i<walberla::MPIManager::instance()->numProcesses();++i) {
    num += dofs_per_rank[i];
  }


  enumerate_impl( level, start );

  stopTiming( "Enumerate" );
  return counter;
}

}

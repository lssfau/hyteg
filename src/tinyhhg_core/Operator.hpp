
#pragma once

#include <core/DataTypes.h>
#include <core/timing/TimingTree.h>
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

#include <memory>

namespace hhg
{

template< typename SourceFunction, typename DestinationFunction >
class Operator
{
public:
  Operator(const std::shared_ptr<PrimitiveStorage> & storage, uint_t minLevel, uint_t maxLevel)
    : storage_(storage), minLevel_(minLevel), maxLevel_(maxLevel)
  {
  }

  virtual ~Operator()
  {
  }

  void apply( SourceFunction& src, DestinationFunction& dst, size_t level, DoFType flag, UpdateType updateType = Replace );

  void smooth_gs( DestinationFunction& dst, SourceFunction& rhs, size_t level, DoFType flag );

  void enableTiming( const std::shared_ptr< walberla::WcTimingTree > & timingTree ) { timingTree_ = timingTree; }

 protected:

  virtual void apply_impl( SourceFunction& src, DestinationFunction& dst, size_t level, DoFType flag, UpdateType updateType = Replace ) = 0;
  virtual void smooth_gs_impl( DestinationFunction& dst, SourceFunction& rhs, size_t level, DoFType flag ) = 0;

  const std::shared_ptr< PrimitiveStorage > storage_;
  const uint_t minLevel_;
  const uint_t maxLevel_;

 private:

  std::shared_ptr< walberla::WcTimingTree > timingTree_;

  void startTiming( const std::string & timerString )
  {
    if ( timingTree_ )
    {
      timingTree_->start( "Operator" );
      timingTree_->start( timerString );
    }
  }

  void stopTiming ( const std::string & timerString )
  {
    if ( timingTree_ )
    {
      timingTree_->stop( timerString );
      timingTree_->stop( "Operator" );
    }
  }
};


template< typename SourceFunction, typename DestinationFunction >
void Operator< SourceFunction, DestinationFunction  >::apply( SourceFunction& src, DestinationFunction& dst, size_t level, DoFType flag, UpdateType updateType )
{
  startTiming( "Apply" );

  apply_impl( src, dst, level, flag, updateType );

  stopTiming( "Apply" );
}

template< typename SourceFunction, typename DestinationFunction >
void Operator< SourceFunction, DestinationFunction  >::smooth_gs( DestinationFunction& dst, SourceFunction& rhs, size_t level, DoFType flag )
{
  startTiming( "Smooth GS" );

  smooth_gs_impl( dst, rhs, level, flag );

  stopTiming( "Smooth GS" );
}

}


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
  if(storage->getTimingTree()){
      enableTiming(storage->getTimingTree());
    }
  }

  typedef SourceFunction srcType;
  typedef DestinationFunction dstType;

  virtual ~Operator() = default;

  void apply( SourceFunction& src, DestinationFunction& dst, size_t level, DoFType flag, UpdateType updateType = Replace ) const;

  void smooth_gs( DestinationFunction& dst, SourceFunction& rhs, size_t level, DoFType flag );

  void smooth_sor( DestinationFunction& dst, SourceFunction& rhs, real_t relax, size_t level, DoFType flag );

  void smooth_jac( DestinationFunction& dst, SourceFunction& rhs, DestinationFunction& tmp, size_t level, DoFType flag );

  void enableTiming( const std::shared_ptr< walberla::WcTimingTree > & timingTree ) { timingTree_ = timingTree; }

  const std::shared_ptr< PrimitiveStorage > getStorage() const { return storage_; }

 protected:

  virtual void apply_impl( SourceFunction& src, DestinationFunction& dst, size_t level, DoFType flag, UpdateType updateType = Replace ) const = 0;
  virtual void smooth_gs_impl( DestinationFunction& dst, SourceFunction& rhs, size_t level, DoFType flag ) {
    WALBERLA_ASSERT(false, "Not implemented");
  };

  virtual void smooth_sor_impl( DestinationFunction& dst, SourceFunction& rhs, real_t relax, size_t level, DoFType flag ) {
    WALBERLA_ASSERT(false, "Not implemented");
  };

  virtual void smooth_jac_impl( DestinationFunction& dst, SourceFunction& rhs, DestinationFunction& tmp, size_t level, DoFType flag ) {
    WALBERLA_ASSERT(false, "Not implemented");
  };

  const std::shared_ptr< PrimitiveStorage > storage_;
  const uint_t minLevel_;
  const uint_t maxLevel_;

 private:

  std::shared_ptr< walberla::WcTimingTree > timingTree_;

 protected:

  void startTiming( const std::string & timerString ) const
  {
    if ( timingTree_ )
    {
      timingTree_->start( "Operator" );
      timingTree_->start( timerString );
    }
  }

  void stopTiming ( const std::string & timerString ) const
  {
    if ( timingTree_ )
    {
      timingTree_->stop( timerString );
      timingTree_->stop( "Operator" );
    }
  }
};


template< typename SourceFunction, typename DestinationFunction >
void Operator< SourceFunction, DestinationFunction  >::apply( SourceFunction& src, DestinationFunction& dst, size_t level, DoFType flag, UpdateType updateType ) const
{
  apply_impl( src, dst, level, flag, updateType );
}

template< typename SourceFunction, typename DestinationFunction >
void Operator< SourceFunction, DestinationFunction  >::smooth_gs( DestinationFunction& dst, SourceFunction& rhs, size_t level, DoFType flag )
{
  smooth_gs_impl( dst, rhs, level, flag );
}

template< typename SourceFunction, typename DestinationFunction >
void Operator< SourceFunction, DestinationFunction  >::smooth_sor( DestinationFunction& dst, SourceFunction& rhs, real_t relax, size_t level, DoFType flag )
{
  smooth_sor_impl( dst, rhs, relax, level, flag );
}

template< typename SourceFunction, typename DestinationFunction >
void Operator< SourceFunction, DestinationFunction  >::smooth_jac( DestinationFunction& dst, SourceFunction& rhs, DestinationFunction& tmp, size_t level, DoFType flag )
{
  smooth_jac_impl( dst, rhs, tmp, level, flag );
}

}

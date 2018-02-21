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

  Function(std::string name, const std::shared_ptr<PrimitiveStorage> & storage, uint_t minLevel, uint_t maxLevel)
      : functionName_(std::move(name))
      , storage_(storage)
      , minLevel_(minLevel)
      , maxLevel_(maxLevel)
  {
    for ( uint_t level = minLevel; level <= maxLevel; level++ )
    {
      communicators_[ level ] = std::make_shared< communication::BufferedCommunicator >( storage );
    }
    if(storage->getTimingTree()){
      enableTiming(storage->getTimingTree());
    }
  }

  virtual ~Function() {}

  inline void interpolate( std::function< ValueType( const Point3D & ) > & expr,
                           uint_t                                          level,
                           DoFType                                         flag = All );

  inline void interpolateExtended( std::function< ValueType( const Point3D &, const std::vector<ValueType>& ) > & expr,
                                   const std::vector<FunctionType*>&               srcFunctions,
                                   uint_t                                          level,
                                   DoFType                                         flag = All );

  inline void prolongate( uint_t  level,
                          DoFType flag = All );

  inline void prolongateQuadratic( uint_t  level,
                                   DoFType flag = All );

  inline void restrict( uint_t  level,
                        DoFType flag = All );

  inline uint_t enumerate( uint_t   level,
                           uint_t & num );


  const std::string &getFunctionName() const { return functionName_; }

  const std::shared_ptr< PrimitiveStorage > getStorage() const { WALBERLA_ASSERT( !!(storage_.lock()) ); return storage_.lock(); }

  uint_t getMinLevel() const { return minLevel_; }

  uint_t getMaxLevel() const { return maxLevel_; }

  const std::shared_ptr<communication::BufferedCommunicator> & getCommunicator(uint_t level) const
  {
    WALBERLA_ASSERT(level >= minLevel_ && level <= maxLevel_);
    return communicators_.at( level );
  };

  void enableTiming( const std::shared_ptr< walberla::WcTimingTree > & timingTree )
  {
    timingTree_ = timingTree;
    for ( auto & communicator : communicators_ )
    {
      communicator.second->enableTiming( timingTree_ );
    }
  }

protected:

    virtual void
    interpolate_impl( std::function< ValueType( const Point3D&, const std::vector<ValueType>& ) >& expr,
                      const std::vector<FunctionType*> srcFunctions,
                      uint_t level, DoFType flag = All ) = 0;

    virtual void
    prolongate_impl( uint_t level, DoFType flag = All ) = 0;

    virtual void
    prolongateQuadratic_impl( uint_t level, DoFType flag = All ) = 0;

    virtual void
    restrict_impl( uint_t level, DoFType flag = All ) = 0;

    virtual void
    enumerate_impl( uint_t level, uint_t& num ) = 0;

  const std::string functionName_;
  const std::weak_ptr< PrimitiveStorage > storage_;
  const uint_t minLevel_;
  const uint_t maxLevel_;

  std::map< uint_t, std::shared_ptr< communication::BufferedCommunicator > > communicators_;

  std::shared_ptr< walberla::WcTimingTree > timingTree_;

  void startTiming( const std::string & timerString )
  {
    if ( timingTree_ )
    {
      timingTree_->start( "Function" );
      timingTree_->start( FunctionTrait< FunctionType >::getTypeName() );
      timingTree_->start( timerString );
    }
  }

  void stopTiming ( const std::string & timerString )
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
void Function< FunctionType >::interpolate(std::function< ValueType(const Point3D&)>& expr, uint_t level, DoFType flag)
{
  std::function< ValueType(const Point3D&,const std::vector<ValueType>&)> exprExtended = [&expr](const hhg::Point3D& x, const std::vector<ValueType>&) {
    return expr(x);
  };

  startTiming( "Interpolate" );

  interpolate_impl( exprExtended, {}, level, flag );

  stopTiming( "Interpolate" );
}

template< typename FunctionType >
void Function< FunctionType >::interpolateExtended(std::function< ValueType(const Point3D&,const std::vector<ValueType>&)>& expr,
                                                   const std::vector<FunctionType*>& srcFunctions,
                                                   uint_t level,
                                                   DoFType flag)
{
  startTiming( "InterpolateExt" );

  interpolate_impl( expr, srcFunctions, level, flag );

  stopTiming( "InterpolateExt" );
}

template< typename FunctionType >
void Function< FunctionType >::prolongate(size_t level, DoFType flag)
{
  startTiming( "Prolongate" );

  prolongate_impl( level, flag );

  stopTiming( "Prolongate" );
}

template< typename FunctionType >
void Function< FunctionType >::prolongateQuadratic(size_t level, DoFType flag)
{
  startTiming( "Prolongate Quadratic" );

  prolongateQuadratic_impl( level, flag );

  stopTiming( "Prolongate Quadratic" );
}

template< typename FunctionType >
void Function< FunctionType >::restrict(size_t level, DoFType flag)
{
  startTiming( "Restrict" );

  restrict_impl( level, flag );

  stopTiming( "Restrict" );
}

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

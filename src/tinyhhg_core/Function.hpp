#pragma once

#include <tinyhhg_core/Operator.hpp>
#include "mesh.hpp"
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/types/flags.hpp"
#include "tinyhhg_core/communication/BufferedCommunication.hpp"

#include <string>
#include <functional>

namespace hhg {

template< typename FunctionType >
class Function {
public:
  Function(const std::string& name, const std::shared_ptr<PrimitiveStorage> & storage, size_t minLevel, size_t maxLevel)
      : functionName_(name)
      , storage_(storage)
      , minLevel_(minLevel)
      , maxLevel_(maxLevel)
  {
    for ( uint_t level = minLevel; level <= maxLevel; level++ )
    {
      communicators_[ level ] = std::shared_ptr< communication::BufferedCommunicator >( new communication::BufferedCommunicator( storage ) );
    }
  }

  virtual ~Function()
  {
  }

  inline void interpolate(std::function<real_t(const Point3D&)>& expr, uint_t level, DoFType flag = All);

  inline void assign(const std::vector<walberla::real_t> scalars, const std::vector<FunctionType*> functions, size_t level, DoFType flag = All);

  inline void add(const std::vector<walberla::real_t> scalars, const std::vector<FunctionType*> functions, size_t level, DoFType flag = All);

  inline real_t dot(FunctionType& rhs, size_t level, DoFType flag = All);

  inline void prolongate(size_t level, DoFType flag = All);

  inline void prolongateQuadratic(size_t level, DoFType flag = All);

  inline void restrict(size_t level, DoFType flag = All);


  const std::string &getFunctionName() const { return functionName_; }

  const std::shared_ptr< PrimitiveStorage > getStorage() const { return storage_; }

  uint_t getMinLevel() const { return minLevel_; }

  uint_t getMaxLevel() const { return maxLevel_; }

  std::shared_ptr<communication::BufferedCommunicator>& getCommunicator(uint_t level) {
    WALBERLA_ASSERT(level >= minLevel_ && level <= maxLevel_);
    return communicators_[level];
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

  virtual void interpolate_impl(std::function<real_t(const Point3D&)>& expr, uint_t level, DoFType flag = All) = 0;

  virtual void assign_impl(const std::vector<walberla::real_t> scalars, const std::vector<FunctionType*> functions, size_t level, DoFType flag = All) = 0;

  virtual void add_impl(const std::vector<walberla::real_t> scalars, const std::vector<FunctionType*> functions, size_t level, DoFType flag = All) = 0;

  virtual real_t dot_impl(FunctionType& rhs, size_t level, DoFType flag = All) = 0;

  virtual void prolongate_impl(size_t level, DoFType flag = All) = 0;

  virtual void prolongateQuadratic_impl(size_t level, DoFType flag = All) = 0;

  virtual void restrict_impl(size_t level, DoFType flag = All) = 0;

  const std::string functionName_;
  const std::shared_ptr< PrimitiveStorage > storage_;
  const uint_t minLevel_;
  const uint_t maxLevel_;

  std::map< uint_t, std::shared_ptr< communication::BufferedCommunicator > > communicators_;

  std::shared_ptr< walberla::WcTimingTree > timingTree_;

private:

  void startTiming( const std::string & timerString )
  {
    if ( timingTree_ )
    {
      timingTree_->start( "Function" );
      timingTree_->start( timerString );
    }
  }

  void stopTiming ( const std::string & timerString )
  {
    if ( timingTree_ )
    {
      timingTree_->stop( timerString );
      timingTree_->stop( "Function" );
    }
  }

};


template< typename FunctionType >
void Function< FunctionType >::interpolate(std::function<real_t(const Point3D&)>& expr, uint_t level, DoFType flag)
{
  startTiming( "Interpolate" );

  interpolate_impl( expr, level, flag );

  stopTiming( "Interpolate" );
}

template< typename FunctionType >
void Function< FunctionType >::assign(const std::vector<walberla::real_t> scalars, const std::vector<FunctionType*> functions, size_t level, DoFType flag)
{
  startTiming( "Assign" );

  assign_impl( scalars, functions, level, flag );

  stopTiming( "Assign" );
}

template< typename FunctionType >
void Function< FunctionType >::add(const std::vector<walberla::real_t> scalars, const std::vector<FunctionType*> functions, size_t level, DoFType flag)
{
  startTiming( "Add" );

  add_impl( scalars, functions, level, flag );

  stopTiming( "Add" );
}

template< typename FunctionType >
real_t Function< FunctionType >::dot(FunctionType& rhs, size_t level, DoFType flag)
{
  startTiming( "Dot" );

  real_t dotResult = dot_impl( rhs, level, flag );

  stopTiming( "Dot" );

  return dotResult;
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

}

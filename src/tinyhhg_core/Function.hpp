#pragma once

#include <tinyhhg_core/Operator.hpp>
#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/types/flags.hpp"
#include "tinyhhg_core/communication/BufferedCommunication.hpp"
#include <core/mpi/Gather.h>

#include <string>
#include <functional>
#include <vector>


namespace hhg {

template< typename FunctionType >
class Function {
public:
  Function(const std::string& name, const std::shared_ptr<PrimitiveStorage> & storage, uint_t minLevel, uint_t maxLevel)
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

  inline void assign(const std::vector<walberla::real_t> scalars, const std::vector<FunctionType*> functions, uint_t level, DoFType flag = All);

  inline void add(const std::vector<walberla::real_t> scalars, const std::vector<FunctionType*> functions, uint_t level, DoFType flag = All);

  inline real_t dot(FunctionType& rhs, uint_t level, DoFType flag = All);

  inline void prolongate(uint_t level, DoFType flag = All);

  inline void prolongateQuadratic(uint_t level, DoFType flag = All);

  inline void restrict(uint_t level, DoFType flag = All);

  inline void enumerate(uint_t level, uint_t& num);

  inline void createVectorFromFunction(FunctionType &numerator,Vec &vec, uint_t level, DoFType flag);

  inline void createFunctionFromVector(FunctionType &numerator, Vec &vec, uint_t level, DoFType flag);

  inline void applyDirichletBC(Mat &mat, uint_t level);


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

  virtual void assign_impl(const std::vector<walberla::real_t> scalars, const std::vector<FunctionType*> functions, uint_t level, DoFType flag = All) = 0;

  virtual void add_impl(const std::vector<walberla::real_t> scalars, const std::vector<FunctionType*> functions, uint_t level, DoFType flag = All) = 0;

  virtual real_t dot_impl(FunctionType& rhs, uint_t level, DoFType flag = All) = 0;

  virtual void prolongate_impl(uint_t level, DoFType flag = All) = 0;

  virtual void prolongateQuadratic_impl(uint_t level, DoFType flag = All) = 0;

  virtual void restrict_impl(uint_t level, DoFType flag = All) = 0;

  virtual void enumerate_impl(uint_t level, uint_t& num) = 0;

  virtual void createVectorFromFunction_impl(FunctionType &numerator,Vec &vec, uint_t level,DoFType flag){;} //TODO make this abstract

  virtual void createFunctionFromVector_impl(FunctionType &numerator, Vec &vec, uint_t level, DoFType flag){;}

  virtual void applyDirichletBC_impl(Mat &mat, uint_t level){;}

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

template< typename FunctionType >
void Function< FunctionType >::enumerate(size_t level, uint_t& num)
{
  startTiming( "Enumerate" );
  uint_t counter = 0;


  enumerate_impl(level,counter);

  std::vector<uint_t> dofs_per_rank  = walberla::mpi::allGather(counter);

  uint_t start = num;

  for(int i = 0; i<walberla::MPIManager::instance()->rank();++i)
    start += dofs_per_rank[i];

  for(int i = 0; i<walberla::MPIManager::instance()->numProcesses();++i)
    num += dofs_per_rank[i];


  enumerate_impl( level, start );

  stopTiming( "Enumerate" );
}


template< typename FunctionType >
void Function< FunctionType >::createVectorFromFunction(FunctionType &numerator,Vec &vec, uint_t level,DoFType flag)
{
  startTiming( "createVectorFromFunction" );

  createVectorFromFunction_impl(numerator, vec, level,flag);

  stopTiming( "createVectorFromFunction" );
}


template< typename FunctionType >
void Function< FunctionType >::createFunctionFromVector(FunctionType &numerator,Vec &vec, uint_t level,DoFType flag)
{
  startTiming( "createFunctionFromVector" );

  createFunctionFromVector_impl(numerator, vec, level,flag);

  stopTiming( "createFunctionFromVector" );
}

template< typename FunctionType >
void Function< FunctionType >::applyDirichletBC(Mat &mat, uint_t level)
{
  startTiming( "applyDirichletBC" );

  applyDirichletBC_impl(mat, level);

  stopTiming( "applyDirichletBC" );
}


}

#pragma once

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"

namespace hhg
{

template <typename ValueType>
class P2P1TaylorHoodFunction
{
public:

  P2P1TaylorHoodFunction(const std::string& _name, const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
    : u(_name+"_u", storage, minLevel, maxLevel),
      v(_name+"_v", storage, minLevel, maxLevel),
      p(_name+"_p", storage, minLevel, maxLevel)
  {
  }

  void interpolate(std::function<real_t(const hhg::Point3D&)>& expr, size_t level, DoFType flag = All)
  {
    u.interpolate(expr, level, flag);
    v.interpolate(expr, level, flag);
    p.interpolate(expr, level, flag | DirichletBoundary);
  }

  void assign(const std::vector<walberla::real_t> scalars, const std::vector<P2P1TaylorHoodFunction<ValueType>*> functions, size_t level, DoFType flag = All)
  {
    std::vector< P2Function< ValueType > * > functions_u;
    std::vector< P2Function< ValueType > * > functions_v;
    std::vector< P1Function< ValueType > * > functions_p;

    for (auto& function : functions)
    {
      functions_u.push_back(&function->u);
      functions_v.push_back(&function->v);
      functions_p.push_back(&function->p);
    }

    u.assign(scalars, functions_u, level, flag);
    v.assign(scalars, functions_v, level, flag);
    p.assign(scalars, functions_p, level, flag | DirichletBoundary);
  }

  void add(const std::vector<walberla::real_t> scalars, const std::vector<P2P1TaylorHoodFunction<ValueType>*> functions, size_t level, DoFType flag = All)
  {
    std::vector< P2Function< ValueType > * > functions_u;
    std::vector< P2Function< ValueType > * > functions_v;
    std::vector< P1Function< ValueType > * > functions_p;

    for (auto& function : functions)
    {
      functions_u.push_back(&function->u);
      functions_v.push_back(&function->v);
      functions_p.push_back(&function->p);
    }

    u.add(scalars, functions_u, level, flag);
    v.add(scalars, functions_v, level, flag);
    p.add(scalars, functions_p, level, flag | DirichletBoundary);
  }

  walberla::real_t dotGlobal(P2P1TaylorHoodFunction<ValueType>& rhs, size_t level, DoFType flag = All)
  {
    // this->startTiming( "Dot" );
    walberla::real_t sum = u.dotLocal(rhs.u, level, flag);
    sum += v.dotLocal(rhs.v, level, flag);
    sum += p.dotLocal(rhs.p, level, flag | DirichletBoundary);
    walberla::mpi::allReduceInplace( sum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
    // this->stopTiming( "Dot" );
    return sum;
  }

  void prolongate(size_t level, DoFType flag = All)
  {
    u.prolongate(level, flag);
    v.prolongate(level, flag);
    p.prolongate(level, flag | DirichletBoundary);
  }

  void restrict(size_t level, DoFType flag = All)
  {
    u.restrict(level, flag);
    v.restrict(level, flag);
    p.restrict(level, flag | DirichletBoundary);
  }

  void enableTiming( const std::shared_ptr< walberla::WcTimingTree > & timingTree )
  {
    u.enableTiming(timingTree);
    v.enableTiming(timingTree);
    p.enableTiming(timingTree);
  }

  uint_t enumerate( uint_t level, uint_t& num )
  {
    uint_t counter = 0;
    counter += u.enumerate( level, num );
    counter += v.enumerate( level, num );
    counter += p.enumerate( level, num );
    return counter;
  }


  P2Function< ValueType > u;
  P2Function< ValueType > v;
  P1Function< ValueType > p;
};

}

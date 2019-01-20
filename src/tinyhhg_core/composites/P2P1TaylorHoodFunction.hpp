#pragma once

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/FunctionTraits.hpp"

namespace hhg
{

template <typename ValueType>
class P2P1TaylorHoodFunction 
{
public:

  typedef ValueType valueType;

  template< typename VType >
  using FunctionType = P2P1TaylorHoodFunction< VType >;

  typedef P2Function< ValueType > VelocityFunction_T;
  typedef P1Function< ValueType > PressureFunction_T;
  typedef typename FunctionTrait< P2P1TaylorHoodFunction< ValueType > >::Tag Tag;

  P2P1TaylorHoodFunction( const std::string&                         _name,
                          const std::shared_ptr< PrimitiveStorage >& storage,
                          size_t                                     minLevel,
                          size_t                                     maxLevel )
  : u( _name + "_u", storage, minLevel, maxLevel )
  , v( _name + "_v", storage, minLevel, maxLevel )
  , w( storage->hasGlobalCells() ? P2Function< ValueType >( _name + "_w", storage, minLevel, maxLevel ) :
                                   P2Function< ValueType >( _name + "_w_dummy", storage ) )
  , p( _name + "_p", storage, minLevel, maxLevel, BoundaryCondition::createAllInnerBC() )
  {}

  P2P1TaylorHoodFunction(const std::string& _name, const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel, BoundaryCondition velocityBC)
        : u(_name+"_u", storage, minLevel, maxLevel, velocityBC)
        , v(_name+"_v", storage, minLevel, maxLevel, velocityBC)
        , w( storage->hasGlobalCells() ? P2Function< ValueType >( _name+"_w", storage, minLevel, maxLevel, velocityBC ) :  P2Function< ValueType >( _name+"_w_dummy", storage ))
        , p(_name+"_p", storage, minLevel, maxLevel, BoundaryCondition::createAllInnerBC() )
  {}

  void interpolate(const std::function<real_t(const hhg::Point3D&)>& expr, size_t level, DoFType flag = All) const
  {
    u.interpolate(expr, level, flag);
    v.interpolate(expr, level, flag);
    w.interpolate(expr, level, flag);
    p.interpolate(expr, level, flag);
  }

  void swap( const P2P1TaylorHoodFunction< ValueType > & other,
             const uint_t & level,
             const DoFType & flag = All ) const
  {
    u.swap( other.u, level, flag );
    v.swap( other.v, level, flag );
    w.swap( other.w, level, flag );
    p.swap( other.p, level, flag );
  }

  void assign( const std::vector< walberla::real_t >                                                     scalars,
               const std::vector< std::reference_wrapper< const P2P1TaylorHoodFunction< ValueType > > >& functions,
               size_t                                                                                    level,
               DoFType                                                                                   flag = All ) const
  {
    std::vector< std::reference_wrapper< const VelocityFunction_T > > functions_u;
    std::vector< std::reference_wrapper< const VelocityFunction_T > > functions_v;
    std::vector< std::reference_wrapper< const VelocityFunction_T > > functions_w;
    std::vector< std::reference_wrapper< const PressureFunction_T > > functions_p;

    for (const P2P1TaylorHoodFunction< ValueType > & function : functions)
    {
      functions_u.push_back(function.u);
      functions_v.push_back(function.v);
      functions_w.push_back(function.w);
      functions_p.push_back(function.p);
    }

    u.assign(scalars, functions_u, level, flag);
    v.assign(scalars, functions_v, level, flag);
    w.assign(scalars, functions_w, level, flag);
    p.assign(scalars, functions_p, level, flag);
  }

  void add( const std::vector< walberla::real_t >                                                     scalars,
            const std::vector< std::reference_wrapper< const P2P1TaylorHoodFunction< ValueType > > >& functions,
            size_t                                                                                    level,
            DoFType                                                                                   flag = All ) const
  {
     std::vector< std::reference_wrapper< const VelocityFunction_T > > functions_u;
     std::vector< std::reference_wrapper< const VelocityFunction_T > > functions_v;
     std::vector< std::reference_wrapper< const VelocityFunction_T > > functions_w;
     std::vector< std::reference_wrapper< const PressureFunction_T > > functions_p;

     for( const P2P1TaylorHoodFunction< ValueType >& function : functions )
     {
        functions_u.push_back( function.u );
        functions_v.push_back( function.v );
        functions_w.push_back( function.w );
        functions_p.push_back( function.p );
     }

     u.add( scalars, functions_u, level, flag );
     v.add( scalars, functions_v, level, flag );
     w.add( scalars, functions_w, level, flag );
     p.add( scalars, functions_p, level, flag );
  }

  walberla::real_t dotGlobal(const P2P1TaylorHoodFunction< ValueType >& rhs,const size_t level,const DoFType flag = All ) const
  {
    walberla::real_t sum = u.dotLocal(rhs.u, level, flag);
    sum += v.dotLocal(rhs.v, level, flag);
    sum += w.dotLocal(rhs.w, level, flag);
    sum += p.dotLocal(rhs.p, level, flag | DirichletBoundary);
    walberla::mpi::allReduceInplace( sum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
    return sum;
  }

  void prolongate(const size_t level,const DoFType flag = All) const
  {
    u.prolongate(level, flag);
    v.prolongate(level, flag);
    w.prolongate(level, flag);
    p.prolongate(level, flag);
  }

  void restrict(const size_t level,const DoFType flag = All) const
  {
    u.restrict(level, flag);
    v.restrict(level, flag);
    w.restrict(level, flag);
    p.restrict(level, flag);
  }

  void enableTiming( const std::shared_ptr< walberla::WcTimingTree > & timingTree ) const
  {
    u.enableTiming(timingTree);
    v.enableTiming(timingTree);
    w.enableTiming(timingTree);
    p.enableTiming(timingTree);
  }

  void enumerate( uint_t level ) const
  {
    uint_t counterDoFs = hhg::numberOfLocalDoFs< Tag >( *( u.getStorage() ), level );

    std::vector< uint_t > doFsPerRank = walberla::mpi::allGather( counterDoFs );

    ValueType offset = 0;

    for( uint_t i = 0; i < uint_c( walberla::MPIManager::instance()->rank() ); ++i )
    {
      offset += static_cast< ValueType >( doFsPerRank[i] );
    }

    u.enumerate( level, offset );
    v.enumerate( level, offset );
    w.enumerate( level, offset );
    p.enumerate( level, offset );
  }


  VelocityFunction_T u;
  VelocityFunction_T v;
  VelocityFunction_T w;
  PressureFunction_T p;
};

}

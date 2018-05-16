#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hhg;

int main( int argc, char* argv[] )
{
   /// create enviroment
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   const size_t level = 3;

   /// create timingTree
   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );

   /// read mesh file and create storage
   MeshInfo              meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hhg::loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   /// create operator and functions
   hhg::P2ConstantLaplaceOperator L( storage, level, level );

   hhg::P2Function< real_t > residuumFunction( "residuumFunction", storage, level, level );
   hhg::P2Function< real_t > rightHandSide( "rightHandSide", storage, level, level );
   hhg::P2Function< real_t > u( "u", storage, level, level );
   hhg::P2Function< real_t > u_exact( "u_exact", storage, level, level );
   hhg::P2Function< real_t > error( "error", storage, level, level );
   hhg::P2Function< real_t > npoints_helper( "npoints_helper", storage, level, level );

   std::function< real_t( const hhg::Point3D& ) > exact = []( const hhg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
   std::function< real_t( const hhg::Point3D& ) > ones  = []( const hhg::Point3D& ) { return 1.0; };

   u.interpolate( exact, level, hhg::DirichletBoundary );
   u_exact.interpolate( exact, level );

   real_t residuum = 100;

   uint_t iterations;
   /// use different checks for debug and release since the runtime in debug would be to high otherwise
   WALBERLA_DEBUG_SECTION() { iterations = 100; }
   else { iterations = 5000; }

   /// apply Gauss Seidl smoother i times
   for( uint_t i = 0; i < iterations; ++i )
   {
      L.smooth_gs( u, rightHandSide, level, hhg::Inner );
   }

   /// calculate residuum and check
   L.apply( u, residuumFunction, level, hhg::Inner );
   residuumFunction.add( {-1}, {&rightHandSide}, level, hhg::Inner );
   residuum = std::sqrt( residuumFunction.dotGlobal( residuumFunction, level, hhg::Inner ) );
   WALBERLA_LOG_INFO_ON_ROOT( "residual: = " << residuum );
   WALBERLA_DEBUG_SECTION() { WALBERLA_CHECK_GREATER( 0.2, residuum ); }
   else { WALBERLA_CHECK_GREATER( 1e-14, residuum ); }

   /// calculate and print error
   error.assign( {1.0, -1.0}, {&u, &u_exact}, level );
   npoints_helper.interpolate( ones, level );
   real_t npoints      = npoints_helper.dotGlobal( npoints_helper, level );
   real_t discr_l2_err = std::sqrt( error.dotGlobal( error, level ) / npoints );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );

   //  walberla::WcTimingTree tt = timingTree->getReduced();
   //  WALBERLA_LOG_INFO_ON_ROOT( tt );

   return 0;
}

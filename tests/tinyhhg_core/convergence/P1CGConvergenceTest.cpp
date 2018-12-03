#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"

#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hhg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   const uint_t      level     = 4;
   const std::string meshFile  = "../../data/meshes/quad_8el.msh";

   auto storage = PrimitiveStorage::createFromGmshFile( meshFile );
   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   storage->setTimingTree(timingTree);

   hhg::P1Function< real_t > r( "r", storage, level, level );
   hhg::P1Function< real_t > f( "f", storage, level, level );
   hhg::P1Function< real_t > u( "u", storage, level, level );
   hhg::P1Function< real_t > u_exact( "u_exact", storage, level, level );
   hhg::P1Function< real_t > err( "err", storage, level, level );
   hhg::P1Function< real_t > npoints_helper( "npoints_helper", storage, level, level );

   hhg::P1MassOperator    M( storage, level, level );
   hhg::P1ConstantLaplaceOperator L( storage, level, level );

   std::function< real_t( const hhg::Point3D& ) > exact = []( const hhg::Point3D& x ) {
      return ( 1.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };
   std::function< real_t( const hhg::Point3D& ) > rhs = []( const hhg::Point3D& x ) {
      return ( 3.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };
   std::function< real_t( const hhg::Point3D& ) > ones = []( const hhg::Point3D& ) { return 1.0; };

   u.interpolate( exact, level, hhg::DirichletBoundary );
   u_exact.interpolate( exact, level );
   npoints_helper.interpolate( rhs, level );
   M.apply( npoints_helper, f, level, hhg::All );

   auto solver = hhg::CGSolver< hhg::P1ConstantLaplaceOperator >( storage, level, level );

   solver.solve( L, u, f, level );

   err.assign( {1.0, -1.0}, {u, u_exact}, level );
   npoints_helper.interpolate( ones, level );

   const real_t npoints      = npoints_helper.dotGlobal( npoints_helper, level );
   const real_t discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / npoints );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );
   WALBERLA_CHECK_LESS( discr_l2_err, 1.2e-5 );

   hhg::VTKOutput vtkOutput( "../../output", "P2CGConvergenceTest", storage );
   vtkOutput.add( u );
   vtkOutput.write( level );

   walberla::WcTimingTree tt = timingTree->getReduced();
   WALBERLA_LOG_INFO_ON_ROOT( tt );

   return 0;
}

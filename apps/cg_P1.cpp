#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1Operator.hpp"
#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/petsc/PETScPreconditioner.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
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

#ifdef HHG_BUILD_WITH_PETSC
   PETScManager petscManager;
#endif

   walberla::shared_ptr< walberla::config::Config > cfg( new walberla::config::Config );
   cfg->readParameterFile( "../data/param/cg_P1.prm" );
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );

   size_t      minLevel     = parameters.getParameter< size_t >( "minlevel" );
   size_t      maxLevel     = parameters.getParameter< size_t >( "maxlevel" );
   size_t      maxiter      = parameters.getParameter< size_t >( "maxiter" );
   real_t      tolerance    = parameters.getParameter< real_t >( "tolerance" );
   std::string meshFileName = parameters.getParameter< std::string >( "mesh" );

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hhg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< PrimitiveStorage >       storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   hhg::P1Function< real_t > r( "r", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > f( "f", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > u( "u", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > u_exact( "u_exact", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > err( "err", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > npoints_helper( "npoints_helper", storage, minLevel, maxLevel );

   hhg::P1MassOperator    M( storage, minLevel, maxLevel );
   hhg::P1LaplaceOperator L( storage, minLevel, maxLevel );

   std::function< real_t( const hhg::Point3D& ) > exact = []( const hhg::Point3D& x ) {
      return ( 1.0L / 2.0L ) * sin( 2 * x[0] ) * sinh( x[1] );
   };
   std::function< real_t( const hhg::Point3D& ) > rhs = []( const hhg::Point3D& x ) {
      return ( 3.0L / 2.0L ) * sin( 2 * x[0] ) * sinh( x[1] );
   };
   std::function< real_t( const hhg::Point3D& ) > ones = []( const hhg::Point3D& ) { return 1.0; };

   u.interpolate( exact, maxLevel, hhg::DirichletBoundary );
   u_exact.interpolate( exact, maxLevel );
   npoints_helper.interpolate( rhs, maxLevel );
   M.apply( npoints_helper, f, maxLevel, hhg::All );

#ifdef HHG_BUILD_WITH_PETSC
   typedef hhg::PETScPreconditioner< real_t, hhg::P1Function, hhg::P1LaplaceOperator > PreconditionerType;
   std::shared_ptr< hhg::P1Function< PetscInt > >                                      numerator =
       std::make_shared< hhg::P1Function< PetscInt > >( "numerator", storage, maxLevel, maxLevel );
   uint_t globalSize = 0;
   uint_t localSize  = numerator->enumerate( maxLevel, globalSize );
   auto   prec       = std::make_shared< PreconditionerType >( L, numerator, localSize, globalSize );
#else
   typedef hhg::GaussSeidelPreconditioner< hhg::P1Function< real_t >, hhg::P1LaplaceOperator > PreconditionerType;
   auto prec = std::make_shared< PreconditionerType >( L, 30 );
#endif
   auto solver = hhg::CGSolver< hhg::P1Function< real_t >, hhg::P1LaplaceOperator, PreconditionerType >(
       storage, minLevel, maxLevel, std::numeric_limits< uint_t >::max(), prec );
   walberla::WcTimer timer;
   solver.solve( L, u, f, r, maxLevel, tolerance, maxiter, hhg::Inner, true );
   timer.end();
   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
   err.assign( {1.0, -1.0}, {&u, &u_exact}, maxLevel );

   npoints_helper.interpolate( ones, maxLevel );
   real_t npoints = npoints_helper.dot( npoints_helper, maxLevel );

   real_t discr_l2_err = std::sqrt( err.dot( err, maxLevel ) / npoints );

   //  auto face0data = *storage->beginFaces().operator*().second->getData(u.getFaceDataID());

   //  hhg::P1Face::printFunctionMemory(face0data,maxLevel);

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );

   //hhg::VTKWriter< P1Function< real_t > >({ &u, &u_exact, &f, &r, &err }, maxLevel, "../output", "cg_P1");

   walberla::WcTimingTree tt = timingTree->getReduced();
   WALBERLA_LOG_INFO_ON_ROOT( tt );

   return 0;
}

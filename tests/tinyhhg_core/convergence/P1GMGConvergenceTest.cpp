#include <tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp>
#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/solvers/GeometricMultiGrid.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hhg;

int main( int argc, char* argv[] )
{
  walberla::Environment walberlaEnv( argc, argv );
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  const uint_t      minLevel  = 2;
  const uint_t      maxLevel  = 5;
  const std::string meshFile  = "../../data/meshes/quad_8el.msh";
  const real_t      coarseGridSolverTolerance = 1e-16;
  const uint_t      maxCoarseGridSolverIter   = 10000;
  const uint_t      numVCycles = 10;

  auto storage = PrimitiveStorage::createFromGmshFile( meshFile );

  hhg::P1Function< real_t > r( "r", storage, minLevel, maxLevel );
  hhg::P1Function< real_t > f( "f", storage, minLevel, maxLevel );
  hhg::P1Function< real_t > u( "u", storage, minLevel, maxLevel );
  hhg::P1Function< real_t > Au( "Au", storage, minLevel, maxLevel );
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

  typedef hhg::CGSolver< hhg::P1Function< real_t >, hhg::P1LaplaceOperator > CGSolver;
  typedef hhg::P1toP1LinearRestriction RestrictionOperator;

  RestrictionOperator restrictionOperator;
  auto coarseGridSolver = std::make_shared< CGSolver >( storage, minLevel, minLevel );
  auto gmgSolver = hhg::GMultigridSolver< hhg::P1Function< real_t >, hhg::P1LaplaceOperator, CGSolver, RestrictionOperator >(
    storage, coarseGridSolver, restrictionOperator, minLevel, maxLevel, 3, 3 );

  npoints_helper.interpolate( ones, maxLevel );
  const real_t npoints = npoints_helper.dot( npoints_helper, maxLevel );

  // init residual once
  L.apply(u, Au, maxLevel, hhg::Inner);
  r.assign({1.0, -1.0}, { &f, &Au }, maxLevel, hhg::Inner);

  real_t discr_l2_err;
  real_t discr_l2_res = std::sqrt( r.dot( r, maxLevel, DoFType::Inner ) / npoints );
  real_t discr_l2_res_last_step = discr_l2_res;

  for ( uint_t vCycleCount = 0; vCycleCount < numVCycles; vCycleCount++ )
  {
    gmgSolver.solve( L, u, f, r, maxLevel, coarseGridSolverTolerance, maxCoarseGridSolverIter, DoFType::Inner );

    err.assign( { 1.0, -1.0 }, { &u, &u_exact }, maxLevel );

    L.apply(u, Au, maxLevel, hhg::Inner);
    r.assign({1.0, -1.0}, { &f, &Au }, maxLevel, hhg::Inner);

    discr_l2_err = std::sqrt( err.dot( err, maxLevel, DoFType::All ) / npoints );
    discr_l2_res_last_step = discr_l2_res;
    discr_l2_res = std::sqrt( r.dot( r, maxLevel, DoFType::Inner ) / npoints );

    const real_t convRate = discr_l2_res / discr_l2_res_last_step;

    WALBERLA_LOG_INFO_ON_ROOT( "After " << vCycleCount << " V-cycles: "
                               "discrete L2 error = " << discr_l2_err <<
                               ", discrete L2 residual = " << discr_l2_res <<
                               ", conv rate = " << convRate )
  }
  WALBERLA_CHECK_LESS( discr_l2_err, 2.9e-06 );


  return 0;
}

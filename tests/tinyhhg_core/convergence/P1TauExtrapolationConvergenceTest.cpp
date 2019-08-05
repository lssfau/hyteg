
#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1InjectionRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/solvers/FAS.hpp"
#include "tinyhhg_core/solvers/GeometricMultigridSolver.hpp"
#include "tinyhhg_core/solvers/GaussSeidelSmoother.hpp"

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
  hhg::P1Function< real_t > AIu( "AIu", storage, minLevel, maxLevel );
  hhg::P1Function< real_t > IAu( "IAu", storage, minLevel, maxLevel );
  hhg::P1Function< real_t > u_exact( "u_exact", storage, minLevel, maxLevel );
  hhg::P1Function< real_t > err( "err", storage, minLevel, maxLevel );
  hhg::P1Function< real_t > npoints_helper( "npoints_helper", storage, minLevel, maxLevel );

  hhg::P2Function< real_t > quadraticTmp( "tmp_P2", storage, minLevel, maxLevel - 1 );
  hhg::P2Function< real_t > quadraticRhs( "f_P2", storage, minLevel, maxLevel - 1 );

  hhg::P1ConstantMassOperator    M( storage, minLevel, maxLevel );
  hhg::P2ConstantMassOperator    M_P2( storage, minLevel, maxLevel );
  hhg::P1ConstantLaplaceOperator L( storage, minLevel, maxLevel );

#if 1
  std::function< real_t( const hhg::Point3D& ) > exact = []( const hhg::Point3D& x ) {
      return ( 1.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
  };
  std::function< real_t( const hhg::Point3D& ) > rhs = []( const hhg::Point3D& x ) {
      return ( 3.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
  };
#else
  std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) { return sin(x[0])*sinh(x[1]); };
  std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D& ) { return real_c(0); };
#endif

  std::function< real_t( const hhg::Point3D& ) > ones = []( const hhg::Point3D& ) { return 1.0; };

  u.interpolate( exact, maxLevel, hhg::DirichletBoundary );
  u_exact.interpolate( exact, maxLevel );
  u_exact.interpolate( exact, maxLevel - 1 );
  npoints_helper.interpolate( rhs, maxLevel );
  M.apply( npoints_helper, f, maxLevel, hhg::All );

  auto smoother = std::make_shared< hhg::GaussSeidelSmoother<hhg::P1ConstantLaplaceOperator>  >();
  auto coarseGridSolver = std::make_shared< hhg::CGSolver< hhg::P1ConstantLaplaceOperator > >(
      storage, minLevel, minLevel, maxCoarseGridSolverIter, coarseGridSolverTolerance );
  auto restrictionOperator = std::make_shared< hhg::P1toP1LinearRestriction>();
  auto solutionRestrictionOperator = std::make_shared< hhg::P1toP1InjectionRestriction>();
  auto prolongationOperator = std::make_shared< hhg::P1toP1LinearProlongation >();

  auto gmgSolver = hhg::FASSolver< hhg::P1ConstantLaplaceOperator >(
    storage, smoother, coarseGridSolver, restrictionOperator, solutionRestrictionOperator, prolongationOperator, minLevel, maxLevel, 3, 3 );

  auto gmgSolverTau = hhg::GeometricMultigridSolver< hhg::P1ConstantLaplaceOperator >(
    storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel - 1, 1, 1 );

  npoints_helper.interpolate( ones, maxLevel );
  const real_t npoints = npoints_helper.dotGlobal( npoints_helper, maxLevel );

  npoints_helper.interpolate( ones, maxLevel - 1 );
  const real_t npoints_tau = npoints_helper.dotGlobal( npoints_helper, maxLevel - 1 );

  // init residual once
  L.apply(u, Au, maxLevel, hhg::Inner);
  r.assign({1.0, -1.0}, { f, Au }, maxLevel, hhg::Inner);

  real_t discr_l2_err;
  real_t discr_l2_res = std::sqrt( r.dotGlobal( r, maxLevel, DoFType::Inner ) / npoints );
  real_t discr_l2_res_last_step = discr_l2_res;

  for ( uint_t vCycleCount = 0; vCycleCount < numVCycles; vCycleCount++ )
  {
    gmgSolver.solve( L, u, f, maxLevel);

    err.assign( { 1.0, -1.0 }, { u, u_exact }, maxLevel );

    L.apply(u, Au, maxLevel, hhg::Inner);
    r.assign({1.0, -1.0}, { f, Au }, maxLevel, hhg::Inner);

    discr_l2_err = std::sqrt( err.dotGlobal( err, maxLevel, DoFType::All ) / npoints );
    discr_l2_res_last_step = discr_l2_res;
    discr_l2_res = std::sqrt( r.dotGlobal( r, maxLevel, DoFType::Inner ) / npoints );

    const real_t convRate = discr_l2_res / discr_l2_res_last_step;

    WALBERLA_LOG_INFO_ON_ROOT( "After " << vCycleCount << " V-cycles: "
                                                          "discrete L2 error = " << discr_l2_err <<
                                        ", discrete L2 residual = " << discr_l2_res <<
                                        ", conv rate = " << convRate )
  }

  WALBERLA_CHECK_LESS( discr_l2_err, 2.9e-06 );

  // perform tau-extrapolation

  // I * A * u
  L.apply( u, IAu, maxLevel, hhg::Inner );
  restrictionOperator->restrict( IAu, maxLevel, hhg::All );

  // A * I * u (injection)
  solutionRestrictionOperator->restrict( u, maxLevel, hhg::All );
  L.apply( u, AIu, maxLevel - 1, hhg::All );

  // build RHS (quadratic, then restrict with weighting)
  quadraticTmp.interpolate( rhs, maxLevel - 1, hhg::All );
  M_P2.apply( quadraticTmp, quadraticRhs, maxLevel - 1, hhg::All );
  f.assign( quadraticRhs, maxLevel, hhg::All );

  restrictionOperator->restrict( f, maxLevel, All );
  f.assign( {1.0, -4.0 / 3.0, 4.0 / 3.0}, {f, IAu, AIu}, maxLevel - 1 );

  const uint_t tauMaxLevel = maxLevel - 1;
  
  for ( uint_t vCycleCount = 0; vCycleCount < numVCycles; vCycleCount++ )
  {
    gmgSolverTau.solve( L, u, f, tauMaxLevel );

    err.assign( { 1.0, -1.0 }, { u, u_exact }, tauMaxLevel );

    L.apply(u, Au, tauMaxLevel, hhg::Inner);
    r.assign({1.0, -1.0}, { f, Au }, tauMaxLevel, hhg::Inner);

    discr_l2_err = std::sqrt( err.dotGlobal( err, tauMaxLevel, DoFType::All ) / npoints_tau );
    discr_l2_res_last_step = discr_l2_res;
    discr_l2_res = std::sqrt( r.dotGlobal( r, tauMaxLevel, DoFType::Inner ) / npoints_tau );

    const real_t convRate = discr_l2_res / discr_l2_res_last_step;

    WALBERLA_LOG_INFO_ON_ROOT( "After " << vCycleCount << " V-cycles: "
                               "discrete L2 error = " << discr_l2_err <<
                               ", discrete L2 residual = " << discr_l2_res <<
                               ", conv rate = " << convRate )
  }

  WALBERLA_CHECK_LESS( discr_l2_err, 3.3e-09 );


  return 0;
}

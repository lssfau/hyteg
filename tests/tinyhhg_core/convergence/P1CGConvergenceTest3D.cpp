#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"
#include "core/math/Random.h"

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"
#include "tinyhhg_core/p1functionspace/generated/p1_tet_diffusion.h"
#include "tinyhhg_core/VTKWriter.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hhg;

int main( int argc, char* argv[] )
{
  walberla::Environment walberlaEnv( argc, argv );
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  const uint_t      lowerLevel       = 3;
  const uint_t      higherLevel     = lowerLevel + 1;
  const std::string meshFile        = "../../data/meshes/3D/pyramid_4el.msh";
  const real_t      tolerance       = 1e-17;
  const uint_t      maxIterations   = 10000;
  const bool        writeVTK        = false;
  const bool        enableChecks    = true;

  auto storage = PrimitiveStorage::createFromGmshFile( meshFile );

  WALBERLA_CHECK( storage->hasGlobalCells() );

  writeDomainPartitioningVTK( storage, "../../output", "P1_CG_3D_convergence_partitioning" );

  P1ConstantOperator< p1_tet_diffusion_cell_integral_0_otherwise > laplaceOperator3D( storage, lowerLevel, higherLevel );

  std::function< real_t( const hhg::Point3D& ) > exact = []( const hhg::Point3D & p ) -> real_t
  {
    return sin(p[0]) * sinh(p[1]) * p[2];
  };

  std::function< real_t( const hhg::Point3D& ) > zero = []( const hhg::Point3D & ) -> real_t
  {
    return 0.0;
  };

  std::function< real_t( const hhg::Point3D& ) > one = []( const hhg::Point3D & ) -> real_t
  {
      return 1.0;
  };

  std::function< real_t( const hhg::Point3D& ) > rand = []( const hhg::Point3D & ) -> real_t
  {
    return walberla::math::realRandom( 0.0, 1.0 );
  };

  hhg::P1Function< real_t > res( "r", storage, lowerLevel, higherLevel );
  hhg::P1Function< real_t > f( "f", storage, lowerLevel, higherLevel );
  hhg::P1Function< real_t > u( "u", storage, lowerLevel, higherLevel );
  hhg::P1Function< real_t > uExact( "u_exact", storage, lowerLevel, higherLevel );
  hhg::P1Function< real_t > err( "err", storage, lowerLevel, higherLevel );
  hhg::P1Function< real_t > oneFunction( "oneFunction", storage, lowerLevel, higherLevel );

  u.interpolate( rand, lowerLevel, DoFType::Inner );
  u.interpolate( exact, lowerLevel, DoFType::DirichletBoundary );
  f.interpolate( zero, lowerLevel, DoFType::All );
  res.interpolate( zero, lowerLevel, DoFType::All );
  uExact.interpolate( exact, lowerLevel, DoFType::All );
  oneFunction.interpolate( one, lowerLevel, DoFType::All );

  u.interpolate( rand, higherLevel, DoFType::Inner );
  u.interpolate( exact, higherLevel, DoFType::DirichletBoundary );
  f.interpolate( zero, higherLevel, DoFType::All );
  res.interpolate( zero, higherLevel, DoFType::All );
  uExact.interpolate( exact, higherLevel, DoFType::All );
  oneFunction.interpolate( one, higherLevel, DoFType::All );

  auto solver = hhg::CGSolver< hhg::P1Function< real_t >, P1ConstantOperator< p1_tet_diffusion_cell_integral_0_otherwise > >( storage, lowerLevel, higherLevel );

  WALBERLA_CHECK_LESS( lowerLevel, higherLevel );

  const real_t numPointsLowerLevel  = oneFunction.dot( oneFunction, lowerLevel,  DoFType::Inner );
  const real_t numPointsHigherLevel = oneFunction.dot( oneFunction, higherLevel, DoFType::Inner );

  WALBERLA_ASSERT_EQUAL( walberla::mpi::MPIManager::instance()->numProcesses(), 1 );

  VTKOutput vtkOutput( "../../output", "P1CGConvergenceTest" );
  vtkOutput.set3D();
  vtkOutput.add( &u );
  vtkOutput.add( &err );

  if ( writeVTK )
  {
    vtkOutput.write( higherLevel, 0 );
    vtkOutput.write( lowerLevel, 0 );
  }

  solver.solve( laplaceOperator3D, u, f, res, lowerLevel,  tolerance, maxIterations, hhg::Inner, true );
  solver.solve( laplaceOperator3D, u, f, res, higherLevel, tolerance, maxIterations, hhg::Inner, true );

  err.assign( {1.0, -1.0}, {&u, &uExact}, lowerLevel );
  err.assign( {1.0, -1.0}, {&u, &uExact}, higherLevel );
  laplaceOperator3D.apply( u, res, lowerLevel,  DoFType::Inner );
  laplaceOperator3D.apply( u, res, higherLevel, DoFType::Inner );

  if ( writeVTK )
  {
    vtkOutput.write( higherLevel, 1 );
    vtkOutput.write( lowerLevel, 1 );
  }

  const real_t discrL2ErrLowerLevel  = std::sqrt( err.dot( err, lowerLevel,  DoFType::Inner ) / numPointsLowerLevel );
  const real_t discrL2ErrHigherLevel = std::sqrt( err.dot( err, higherLevel, DoFType::Inner ) / numPointsHigherLevel );

  const real_t discrL2ResLowerLevel  = std::sqrt( res.dot( res, lowerLevel,  DoFType::Inner ) / numPointsLowerLevel );
  const real_t discrL2ResHigherLevel = std::sqrt( res.dot( res, higherLevel, DoFType::Inner ) / numPointsHigherLevel );

  WALBERLA_LOG_INFO( "Residual L2 on level " << lowerLevel  << ": " << std::scientific << discrL2ResLowerLevel  << " | Error L2: " << discrL2ErrLowerLevel );
  WALBERLA_LOG_INFO( "Residual L2 on level " << higherLevel << ": " << std::scientific << discrL2ResHigherLevel << " | Error L2: " << discrL2ErrHigherLevel );

  if ( enableChecks )
  {
    WALBERLA_CHECK_LESS( discrL2ResLowerLevel, 5.7e-17 );
    WALBERLA_CHECK_LESS( discrL2ResHigherLevel, 3.5e-17 );

    // L2 err higher level ~ 0.25 * L2 err lower level
    WALBERLA_CHECK_LESS( discrL2ErrLowerLevel, 8.7e-05 );
    WALBERLA_CHECK_LESS( discrL2ErrHigherLevel, 2.0e-05 );
  }

  return 0;
}

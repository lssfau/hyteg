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

  const uint_t      level           = 4;
  const std::string meshFile        = "../../data/meshes/3D/tet_1el.msh";
  const real_t      tolerance       = 1e-15;
  const uint_t      numIterations   = 300;

  auto storage = PrimitiveStorage::createFromGmshFile( meshFile );

  WALBERLA_CHECK( storage->hasGlobalCells() );

  writeDomainPartitioningVTK( storage, "../../output", "P1_CG_3D_convergence_partitioning" );

  P1ConstantOperator< p1_tet_diffusion_cell_integral_0_otherwise > laplaceOperator3D( storage, level, level );

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

  hhg::P1Function< real_t > res( "r", storage, level, level );
  hhg::P1Function< real_t > f( "f", storage, level, level );
  hhg::P1Function< real_t > u( "u", storage, level, level );
  hhg::P1Function< real_t > uExact( "u_exact", storage, level, level );
  hhg::P1Function< real_t > err( "err", storage, level, level );
  hhg::P1Function< real_t > oneFunction( "oneFunction", storage, level, level );

  u.interpolate( rand, level, DoFType::Inner );
  u.interpolate( one, level, DoFType::DirichletBoundary );
  f.interpolate( zero, level, DoFType::All );
  res.interpolate( zero, level, DoFType::All );
  uExact.interpolate( one, level, DoFType::All );
  oneFunction.interpolate( one, level, DoFType::All );

  auto solver = hhg::CGSolver< hhg::P1Function< real_t >, P1ConstantOperator< p1_tet_diffusion_cell_integral_0_otherwise > >( storage, level, level );

  VTKOutput vtkOutput( "../../output", "P1_CG_3D_convergence_solver", 1 );
  vtkOutput.set3D();

  vtkOutput.add( &u );
  vtkOutput.add( &err );

  real_t discrL2Err;
  real_t discrL2Res;

  const real_t numPoints = oneFunction.dot( oneFunction, level, DoFType::Inner );

  WALBERLA_ASSERT_EQUAL( walberla::mpi::MPIManager::instance()->numProcesses(), 1 );

  for ( uint_t iteration = 0; iteration < numIterations; iteration++ )
  {
    // vtkOutput.write( level, iteration );
    solver.solve( laplaceOperator3D, u, f, res, level, tolerance, 1, hhg::Inner, false );
#if 0
    err.assign( {1.0, -1.0}, {&u, &uExact}, level );
    laplaceOperator3D.apply( u, res, level, DoFType::Inner );
    discrL2Err = std::sqrt( err.dot( err, level, DoFType::Inner) / numPoints );
    discrL2Res = std::sqrt( res.dot( res, level, DoFType::Inner ) / numPoints );
    WALBERLA_LOG_INFO( "Iteration " << std::setw( 10 ) << iteration << ": Residual L2: " << std::scientific << discrL2Res << " | Error L2: " << discrL2Err );
#endif
  }

  err.assign( {1.0, -1.0}, {&u, &uExact}, level );
  laplaceOperator3D.apply( u, res, level, DoFType::Inner );
  discrL2Err = std::sqrt( err.dot( err, level, DoFType::Inner) / numPoints );
  discrL2Res = std::sqrt( res.dot( res, level, DoFType::Inner ) / numPoints );

  WALBERLA_LOG_INFO( "After " << numIterations << " iterations: Residual L2: " << std::scientific << discrL2Res << " | Error L2: " << discrL2Err );
  WALBERLA_CHECK_LESS( discrL2Res, 5.3e-13 );
  WALBERLA_CHECK_LESS( discrL2Err, 1.2e-11 );

  return 0;
}

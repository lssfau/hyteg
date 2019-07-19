#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"
#include "core/math/Random.h"

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"
#include "tinyhhg_core/forms/form_fenics_generated/p1_tet_diffusion.h"
#include "tinyhhg_core/VTKWriter.hpp"

using walberla::real_t;
using walberla::real_c;
using walberla::uint_c;
using walberla::uint_t;

using namespace hhg;

void testLaplace3D( const std::string & meshFile, const uint_t & level )
{
  // Tests (on multiple meshes) if the 3D Laplace operator has the following properties:
  // 1. laplace(u) = 0, if u = const.
  // 2. laplace(u) = 0, if u linear

  const bool   writeVTK   = false;
  const real_t errorLimit = 2.8e-13;

  const auto meshInfo = MeshInfo::fromGmshFile( meshFile );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
  const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

  WALBERLA_CHECK( storage->hasGlobalCells() );

  if ( writeVTK )
    writeDomainPartitioningVTK( storage, "../../output", "P1LaplaceOperatorTest3D_partitioning" );

  P1ConstantLaplaceOperator laplaceOperator3D( storage, level, level );

  std::function< real_t( const hhg::Point3D& ) > zero = []( const hhg::Point3D & ) -> real_t
  {
      return 0.0;
  };

  std::function< real_t( const hhg::Point3D& ) > one = []( const hhg::Point3D & ) -> real_t
  {
      return 1.0;
  };

  std::function< real_t( const hhg::Point3D& ) > linearInX = []( const hhg::Point3D & p ) -> real_t
  {
      return real_c(42) * p[0];
  };

  std::function< real_t( const hhg::Point3D& ) > linearInXYZ = []( const hhg::Point3D & p ) -> real_t
  {
      return real_c(42) * p[0] + p[1] + real_c(1337) * p[2];
  };

  hhg::P1Function< real_t > u          ( "u",           storage, level, level );
  hhg::P1Function< real_t > resultExact( "u_exact",     storage, level, level );
  hhg::P1Function< real_t > result     ( "result",      storage, level, level );
  hhg::P1Function< real_t > err        ( "err",         storage, level, level );
  hhg::P1Function< real_t > oneFunction( "oneFunction", storage, level, level );

  oneFunction.interpolate( one, level, DoFType::All );
  const real_t numPoints  = oneFunction.dotGlobal( oneFunction, level, DoFType::Inner );

  VTKOutput vtkOutput("../../output", "P1LaplaceOperatorTest3D", storage);
  vtkOutput.add( u );
  vtkOutput.add( result );
  vtkOutput.add( resultExact );
  vtkOutput.add( err );

  auto testLaplaceResult = [&]( std::function< real_t( const hhg::Point3D& ) > uFunction,
                                std::function< real_t( const hhg::Point3D& ) > resultExactFunction ) -> real_t
  {
    u.interpolate( uFunction, level, DoFType::All );
    result.interpolate( resultExactFunction, level, DoFType::DirichletBoundary );
    resultExact.interpolate( resultExactFunction, level, DoFType::All );
    err.interpolate( zero, level, DoFType::All );

    laplaceOperator3D.apply( u, result, level, DoFType::Inner );

    err.assign( { 1.0, -1.0 }, { result, resultExact }, level, DoFType::All );
    const real_t discrL2Err = std::sqrt( err.dotGlobal( err, level, DoFType::Inner ) / numPoints );

    return discrL2Err;
  };

  // 1. u = const
  // ------------
  //   a) u = 0
  const real_t errorUZero = testLaplaceResult( zero, zero );
  WALBERLA_LOG_INFO_ON_ROOT( "u = 0: L2 error: " << errorUZero );
  WALBERLA_CHECK_LESS( errorUZero, errorLimit );

  //   b) u = 1
  const real_t errorUOne  = testLaplaceResult( one, zero );
  WALBERLA_LOG_INFO_ON_ROOT( "u = 1: L2 error: " << errorUOne );
  WALBERLA_CHECK_LESS( errorUOne, errorLimit );

  // 2. u linear
  // -----------
  //   a) u linear in x
  const real_t errorULinearInX   = testLaplaceResult( linearInX, zero );
  WALBERLA_LOG_INFO_ON_ROOT( "u linear in x: L2 error: " << errorULinearInX );
  WALBERLA_CHECK_LESS( errorULinearInX, errorLimit );

  //   b) u linear in x, y and z
  const real_t errorULinearInXYZ = testLaplaceResult( linearInXYZ, zero );
  WALBERLA_LOG_INFO_ON_ROOT( "u linear in x, y and z: L2 error: " << errorULinearInXYZ );
  WALBERLA_CHECK_LESS( errorULinearInXYZ, errorLimit );

}

int main( int argc, char* argv[] )
{
  walberla::Environment walberlaEnv( argc, argv );
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  testLaplace3D( "../../data/meshes/3D/tet_1el.msh", 2 );
  testLaplace3D( "../../data/meshes/3D/tet_1el.msh", 3 );
  testLaplace3D( "../../data/meshes/3D/pyramid_2el.msh", 2 );
  testLaplace3D( "../../data/meshes/3D/pyramid_2el.msh", 3 );
  testLaplace3D( "../../data/meshes/3D/pyramid_4el.msh", 3 );
  testLaplace3D( "../../data/meshes/3D/pyramid_tilted_4el.msh", 3 );
  testLaplace3D( "../../data/meshes/3D/regular_octahedron_8el.msh", 3 );

  return 0;
}

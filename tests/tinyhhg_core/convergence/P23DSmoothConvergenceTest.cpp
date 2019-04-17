#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/timing/all.h"
#include "core/math/all.h"

#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/VTKWriter.hpp"

namespace hhg {

static void testP2SmoothConvergence( const uint_t & level, const std::string & meshFile, const uint_t & numIterations, const real_t & expectedL2Error )
{
  MeshInfo mesh  = MeshInfo::fromGmshFile( meshFile );
  SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  P2Function< real_t > x  ( "x",   storage, level, level );
  P2Function< real_t > rhs( "rhs", storage, level, level );
  P2ConstantLaplaceOperator A( storage, level, level );

  VTKOutput vtkOutput("../../output", "P23DSmoothConvergenceTest", storage);
  vtkOutput.add( x );

  std::function< real_t( const Point3D & )> zeros = []( const Point3D & ) { return 0; };
  walberla::math::seedRandomGenerator( 0 );
  std::function< real_t( const Point3D & )> rand  = []( const Point3D & ) { return walberla::math::realRandom( 0.0, 1.0 ); };

  rhs.interpolate( zeros, level, All );
  x.interpolate( zeros, level, DirichletBoundary );
  x.interpolate( rand,  level, Inner );

  real_t discreteL2Norm;

  for ( uint_t step = 0; step < numIterations; step++ )
  {
#if 0
   // if'd out to speed up test
   discreteL2Norm = sqrt( x.dotGlobal( x, level, All ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Iteration " << std::setw(10) << step << " - Discrete L2 Norm: " << std::scientific << discreteL2Norm );
   vtkOutput.write( level, step );
#endif
   A.smooth_gs( x, rhs, level, Inner );
  }

  discreteL2Norm = sqrt( x.dotGlobal( x, level, All ) );
  WALBERLA_LOG_INFO_ON_ROOT( "Discrete L2 norm after " << numIterations << " Gauss-Seidel iterations: " << std::scientific << discreteL2Norm );
  WALBERLA_CHECK_LESS( discreteL2Norm, expectedL2Error );
}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testP2SmoothConvergence( 3, "../../data/meshes/3D/tet_1el.msh", 50, 4.7e-02 );
   hhg::testP2SmoothConvergence( 2, "../../data/meshes/3D/pyramid_2el.msh", 50, 2.1e-06 );
   hhg::testP2SmoothConvergence( 2, "../../data/meshes/3D/pyramid_4el.msh", 50, 2.4e-03 );
   hhg::testP2SmoothConvergence( 2, "../../data/meshes/3D/regular_octahedron_8el.msh", 50, 9.6e-02 );

   return EXIT_SUCCESS;
}

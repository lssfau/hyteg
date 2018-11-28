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

static void testP2SmoothConvergence()
{
  const uint_t level = 2;

  MeshInfo mesh  = MeshInfo::fromGmshFile( "../../data/meshes/quad_16el.msh" );
  SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  auto p2Function      = std::make_shared< P2Function< real_t > >( "p2Function", storage, level, level );
  auto rhs             = std::make_shared< P2Function< real_t > >( "rhs",        storage, level, level );
  auto laplaceOperator = std::make_shared< P2ConstantLaplaceOperator >( storage, level, level );

  VTKOutput vtkOutput("../../output", "P2SmoothConvergenceTest", storage);
  vtkOutput.add( p2Function );

  std::function< real_t( const Point3D & )> zeros = []( const Point3D & ) { return 0; };

  walberla::math::seedRandomGenerator( 0 );
  std::function< real_t( const Point3D & )> rand  = []( const Point3D & ) { return walberla::math::realRandom( 0.0, 1.0 ); };

  rhs->interpolate( zeros, level, All );
  p2Function->interpolate( zeros, level, DirichletBoundary );
  p2Function->interpolate( rand,  level, Inner );

  const uint_t smootherSteps = 2500;
        real_t discreteL2Norm;

  for ( uint_t step = 0; step < smootherSteps; step++ )
  {
#if 0
   // if'd out to speed up test
   discreteL2Norm = sqrt( p2Function->dotGlobal( *p2Function, level, All ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Iteration " << std::setw(10) << step << " - Discrete L2 Norm: " << std::scientific << discreteL2Norm );
   vtkOutput.write( level, step );
#endif
   laplaceOperator->smooth_gs( *p2Function, *rhs, level, Inner );
  }

  discreteL2Norm = sqrt( p2Function->dotGlobal( *p2Function, level, All ) );
  WALBERLA_LOG_INFO_ON_ROOT( "Discrete L2 norm after " << smootherSteps << " Gauss-Seidel iterations: " << std::scientific << discreteL2Norm );

  WALBERLA_CHECK_LESS( discreteL2Norm, 7e-13 );
}

} // namespace hhg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   // walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hhg::testP2SmoothConvergence();

   return EXIT_SUCCESS;
}

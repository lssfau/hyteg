// TODO: Explain
//
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"

using walberla::uint_t;

using namespace hhg;

static real_t xLocMax = 0.0;
static real_t yLocMax = 0.0;
static real_t epsilon = 0.0;

#define TEST_MAX_VALUE 100.0

int main( int argc, char* argv[] )
{
  walberla::debug::enterTestMode();

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  uint_t theLevel = 4;

  // Generate mesh around origin
  MeshInfo meshInfo = MeshInfo::emptyMeshInfo();
  meshInfo = MeshInfo::meshRectangle( Point2D( {0.0, 0.0} ), Point2D( {1.0, 1.0} ), MeshInfo::CRISS, 1, 1 );

  // Generate primitives
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  loadbalancing::roundRobin( setupStorage );
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>( setupStorage );

  // ------------
  //  Test cases
  // ------------

  // Define expressions and functions used for testing

  std::function<real_t(const hhg::Point3D&)> testFunc = []( const hhg::Point3D& x ) {
    real_t distance = std::sqrt( (xLocMax - x[0]) * (xLocMax - x[0]) + (yLocMax - x[1]) * (yLocMax - x[1]) );
    return distance > epsilon ? real_t(0.0) : TEST_MAX_VALUE;
  };

  real_t hmin = std::pow( 2.0, - (real_t)theLevel );
  epsilon = 0.5 * hmin;
  WALBERLA_LOG_INFO_ON_ROOT( "level    = " << theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "h_min    = " << std::scientific << hmin    );
  WALBERLA_LOG_INFO_ON_ROOT( "epsilon  = " << std::scientific << epsilon );

  xLocMax = 0.50;
  yLocMax = 0.50;
  hhg::P1Function< real_t > func1( "testFunc1", storage, theLevel, theLevel );
  func1.interpolate( testFunc, theLevel );

  xLocMax = 0.50;
  yLocMax = 0.25;
  hhg::P1Function< real_t > func2( "testFunc2", storage, theLevel, theLevel );
  func2.interpolate( testFunc, theLevel );

  xLocMax = 0.25;
  yLocMax = 0.25;
  hhg::P1Function< real_t > func3( "testFunc3", storage, theLevel, theLevel );
  func3.interpolate( testFunc, theLevel );

  // Determine maximum
  real_t measure = func1.getMaxValue( theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "Test Case #1: max = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MAX_VALUE );

  measure = func2.getMaxValue( theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "Test Case #2: max = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MAX_VALUE );

  measure = func3.getMaxValue( theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "Test Case #3: max = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MAX_VALUE );

  return EXIT_SUCCESS;
}

// We test finding the maximum, maximum magnitude and minum of a P1 function.
// The corresponding value is selectively placed on a macro vertex, edge and
// face and finally a combined test is done
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
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/dgfunctionspace/DGFunction.hpp"

using walberla::uint_t;

using namespace hhg;

static real_t xLocPos = 0.0;
static real_t yLocPos = 0.0;
static real_t zLocPos = 0.0;
static real_t epsilon = 0.0;
static uint_t counter = 0;

#define TEST_MAX_VALUE  100.0
#define TEST_MIN_VALUE -200.0
#define TEST_MAG_VALUE  200.0

// typedef enum{ FIND_MAX, FIND_MIN, FIND_MAG } findType;

template<class funcType>
void runMagTest( std::string mesg, uint_t theLevel, funcType &dofFunc, std::function<real_t(const hhg::Point3D&)> &testFunc, real_t refVal ) {

  real_t measure = 0.0;

  dofFunc.interpolate( testFunc, theLevel );
  measure = dofFunc.getMaxMagnitude( theLevel );

  WALBERLA_CHECK_FLOAT_EQUAL( measure, refVal );
  WALBERLA_LOG_INFO_ON_ROOT( mesg << std::scientific << measure );
}


int main( int argc, char* argv[] )
{
  walberla::debug::enterTestMode();

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  uint_t theLevel = 4;

  // Generate mesh around origin
  MeshInfo meshInfo = MeshInfo::emptyMeshInfo();
  meshInfo = MeshInfo::meshRectangle( Point2D( {0.0, 0.0} ), Point2D( {1.0, 1.0} ), MeshInfo::CROSS, 1, 1 );

  // Generate primitives
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  loadbalancing::roundRobin( setupStorage );
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>( setupStorage );

  // ------------
  //  Test cases
  // ------------

  // Define expressions and functions used for testing

  std::function<real_t(const hhg::Point3D&)> testFuncMax = []( const hhg::Point3D& x ) {
    real_t distance = std::sqrt( (xLocPos - x[0]) * (xLocPos - x[0]) +
                                 (yLocPos - x[1]) * (yLocPos - x[1]) +
                                 (zLocPos - x[2]) * (zLocPos - x[2]) );
    return distance > epsilon ? real_t(0.0) : TEST_MAX_VALUE;
  };

  std::function<real_t(const hhg::Point3D&)> testFuncMin = []( const hhg::Point3D& x ) {
    real_t distance = std::sqrt( (xLocPos - x[0]) * (xLocPos - x[0]) +
                                 (yLocPos - x[1]) * (yLocPos - x[1]) +
                                 (zLocPos - x[2]) * (zLocPos - x[2]) );
    return distance > epsilon ? real_t(0.0) : TEST_MIN_VALUE;
  };

  std::function<real_t(const hhg::Point3D&)> testFuncCombo = []( const hhg::Point3D& x ) {
    WALBERLA_UNUSED(x);
    real_t retVal = 0.0;
    switch( counter%3 )
      {
      case 0: retVal = -1.0; break;
      case 1: retVal =  2.0; break;
      case 2: retVal =  3.0;
      }
    ++counter;
    return retVal;
  };

  real_t hmin = std::pow( 2.0, - (real_t)theLevel );
  epsilon = 0.5 * hmin;
  WALBERLA_LOG_INFO_ON_ROOT( "level    = " << theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "h_min    = " << std::scientific << hmin    );
  WALBERLA_LOG_INFO_ON_ROOT( "epsilon  = " << std::scientific << epsilon );

  // =============
  //  P1Function 
  // =============

  WALBERLA_LOG_INFO_ON_ROOT( "\n\nP1Function (DoFType=All)\n" );

  real_t measure = 0.0;

  // Special value on macro edge
  xLocPos = 0.50;
  yLocPos = 0.50;

  hhg::P1Function< real_t > funcEdgeMax( "", storage, theLevel, theLevel );
  funcEdgeMax.interpolate( testFuncMax, theLevel );
  measure = funcEdgeMax.getMaxValue( theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "Test #1 (edge     ): maximum   = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MAX_VALUE );

  hhg::P1Function< real_t > funcEdgeMin( "", storage, theLevel, theLevel );
  funcEdgeMin.interpolate( testFuncMin, theLevel );
  measure = funcEdgeMin.getMinValue( theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "                     minumum   = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MIN_VALUE );

  measure = funcEdgeMin.getMaxMagnitude( theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "                     magnitude = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MAG_VALUE );

  // Special value on macro face
  xLocPos = 0.50;
  yLocPos = 0.25;
  hhg::P1Function< real_t > funcFaceMax( "", storage, theLevel, theLevel );
  funcFaceMax.interpolate( testFuncMax, theLevel );
  measure = funcFaceMax.getMaxValue( theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "Test #2 (face     ): maximum   = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MAX_VALUE );

  hhg::P1Function< real_t > funcFaceMin( "", storage, theLevel, theLevel );
  funcFaceMin.interpolate( testFuncMin, theLevel );
  measure = funcFaceMin.getMinValue( theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "                     minumum   = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MIN_VALUE );

  measure = funcFaceMin.getMaxMagnitude( theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "                     magnitude = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MAG_VALUE );

  // Special value on macro vertex
  xLocPos = 0.0;
  yLocPos = 0.0;
  hhg::P1Function< real_t > funcVertMax( "", storage, theLevel, theLevel );
  funcVertMax.interpolate( testFuncMax, theLevel );
  measure = funcVertMax.getMaxValue( theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "Test #3 (vertex   ): maximum   = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MAX_VALUE );

  hhg::P1Function< real_t > funcVertMin( "", storage, theLevel, theLevel );
  funcVertMin.interpolate( testFuncMin, theLevel );
  measure = funcVertMin.getMinValue( theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "                     minimum   = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MIN_VALUE );

  measure = funcVertMin.getMaxMagnitude( theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "                     magnitude = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MAG_VALUE );

  // Combined test
  hhg::P1Function< real_t > comboFunc( "", storage, theLevel, theLevel );
  comboFunc.interpolate( testFuncCombo, theLevel );
  measure = comboFunc.getMaxValue( theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "Test #4 (combo    ): maximum   = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, 3.0 );
  measure = comboFunc.getMinValue( theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "                     minimum   = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, -1.0 );
  measure = comboFunc.getMaxMagnitude( theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "                     magnitude = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, 3.0 );

  WALBERLA_LOG_INFO_ON_ROOT( "\n\n P1Function (DoFType=<varying>)\n" );

  // DoFType test #1
  comboFunc.interpolate( testFuncCombo, theLevel );
  measure = comboFunc.getMaxValue( theLevel, Inner );
  WALBERLA_LOG_INFO_ON_ROOT( "Test #5 (combo    ): maximum   = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, 3.0 );
  measure = comboFunc.getMinValue( theLevel, Inner );
  WALBERLA_LOG_INFO_ON_ROOT( "                     minimum   = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, -1.0 );
  measure = comboFunc.getMaxMagnitude( theLevel, Inner );
  WALBERLA_LOG_INFO_ON_ROOT( "                     magnitude = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, 3.0 );

  // DoFType test #2 (extremum on boundary edge)
  xLocPos = 0.00;
  yLocPos = 0.50;

  comboFunc.interpolate( testFuncMax, theLevel );
  measure = comboFunc.getMaxValue( theLevel, DirichletBoundary );
  WALBERLA_LOG_INFO_ON_ROOT( "Test #6 (bc edge  ): max outer = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MAX_VALUE );

  measure = comboFunc.getMaxValue( theLevel, Inner );
  WALBERLA_LOG_INFO_ON_ROOT( "                     max inner = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, 0.0 );

  comboFunc.interpolate( testFuncMin, theLevel );
  measure = comboFunc.getMinValue( theLevel, DirichletBoundary );
  WALBERLA_LOG_INFO_ON_ROOT( "                     min outer = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MIN_VALUE );

  measure = comboFunc.getMinValue( theLevel, Inner );
  WALBERLA_LOG_INFO_ON_ROOT( "                     min inner = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, 0.0 );

  measure = comboFunc.getMaxMagnitude( theLevel, DirichletBoundary );
  WALBERLA_LOG_INFO_ON_ROOT( "                     mag outer = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MAG_VALUE );

  measure = comboFunc.getMaxMagnitude( theLevel, Inner );
  WALBERLA_LOG_INFO_ON_ROOT( "                     mag outer = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, 0.0 );

  // DoFType test #3 (extremum on boundary vertex)
  xLocPos = 1.00;
  yLocPos = 1.00;

  comboFunc.interpolate( testFuncMax, theLevel );
  measure = comboFunc.getMaxValue( theLevel, DirichletBoundary );
  WALBERLA_LOG_INFO_ON_ROOT( "Test #7 (bc vertex): max outer = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MAX_VALUE );

  measure = comboFunc.getMaxValue( theLevel, Inner );
  WALBERLA_LOG_INFO_ON_ROOT( "                     max inner = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, 0.0 );

  comboFunc.interpolate( testFuncMin, theLevel );
  measure = comboFunc.getMinValue( theLevel, DirichletBoundary );
  WALBERLA_LOG_INFO_ON_ROOT( "                     min outer = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MIN_VALUE );

  measure = comboFunc.getMinValue( theLevel, Inner );
  WALBERLA_LOG_INFO_ON_ROOT( "                     min inner = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, 0.0 );

  measure = comboFunc.getMaxMagnitude( theLevel, DirichletBoundary );
  WALBERLA_LOG_INFO_ON_ROOT( "                     mag outer = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, TEST_MAG_VALUE );

  measure = comboFunc.getMaxMagnitude( theLevel, Inner );
  WALBERLA_LOG_INFO_ON_ROOT( "                     mag inner = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, 0.0 );

  // Test case when DoFType is not present in base mesh
  comboFunc.interpolate( testFuncCombo, theLevel );

  measure = comboFunc.getMaxValue( theLevel, NeumannBoundary );
  WALBERLA_LOG_INFO_ON_ROOT( "Test #8 (Neumann  ): maximum   = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, -std::numeric_limits<real_t>::max() );

  measure = comboFunc.getMinValue( theLevel, NeumannBoundary );
  WALBERLA_LOG_INFO_ON_ROOT( "                     minimum   = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, std::numeric_limits<real_t>::max() );

  measure = comboFunc.getMaxMagnitude( theLevel, NeumannBoundary );
  WALBERLA_LOG_INFO_ON_ROOT( "                     magnitude = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, 0.0 );

  // ============
  //  P2Function 
  // ============

  WALBERLA_LOG_INFO_ON_ROOT( "\n\nP2Function (DoFType=All)\n" );

  theLevel = 2;
  hhg::P2Function< real_t > p2func( "", storage, theLevel, theLevel );

  // Special value on macro edge (vertexdof)
  xLocPos = 0.25;
  yLocPos = 0.25;
  runMagTest( "Test #9 (edge/vert): magnitude = ", theLevel, p2func, testFuncMin, TEST_MAG_VALUE );

  // Special value on macro edge (edgedof)
  xLocPos = 0.125;
  yLocPos = 0.125;
  runMagTest( "Test #A (edge/edge): magnitude = ", theLevel, p2func, testFuncMin, TEST_MAG_VALUE );

  // Special value on macro face (edgedof)
  xLocPos = 0.250;
  yLocPos = 0.125;
  runMagTest( "Test #B (face/edge): magnitude = ", theLevel, p2func, testFuncMin, TEST_MAG_VALUE );

  // ============
  //  DGFunction 
  // ============

  WALBERLA_LOG_INFO_ON_ROOT( "\n\nDGFunction (DoFType=All)\n" );

  theLevel = 2;
  hhg::DGFunction< real_t > dgFunc( "", storage, theLevel, theLevel );
  runMagTest( "Test #C (combo    ): magnitude = ", theLevel, dgFunc, testFuncCombo, 3.0 );

  // ===========
  //  3D Meshes
  // ===========

  WALBERLA_LOG_INFO_ON_ROOT( "\n\n---------------------\n TEST on 3D Meshes" );
  WALBERLA_LOG_INFO_ON_ROOT( "----------------------\n" );

  theLevel = 2;

  // Generate mesh from meshfile (mesh represents unit cube)
  MeshInfo meshInfo3D = MeshInfo::emptyMeshInfo();
  meshInfo3D = MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" );

  // Generate primitives
  SetupPrimitiveStorage setupStorage3D( meshInfo3D,
                                        uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  loadbalancing::roundRobin( setupStorage3D );
  std::shared_ptr<PrimitiveStorage> storage3D = std::make_shared<PrimitiveStorage>( setupStorage3D );

  hhg::P1Function< real_t > p1Func3D( "", storage3D, theLevel, theLevel );
  hhg::P2Function< real_t > p2Func3D( "", storage3D, theLevel, theLevel );
  hhg::DGFunction< real_t > dgFunc3D( "", storage3D, theLevel, theLevel );


  WALBERLA_LOG_INFO_ON_ROOT( "\n\nSingle point test\n" );

  xLocPos = 0.25;
  yLocPos = 0.25;
  zLocPos = 0.25;

  runMagTest( "3D Test P1 function: magnitude =  ", theLevel, p1Func3D, testFuncMin, TEST_MAG_VALUE );
  runMagTest( "3D Test P2 function: magnitude =  ", theLevel, p2Func3D, testFuncMin, TEST_MAG_VALUE );


  WALBERLA_LOG_INFO_ON_ROOT( "\n\nCombo test\n" );

  runMagTest( "3D Test P1 function: magnitude =  ", theLevel, p1Func3D, testFuncCombo, 3.0 );

  measure = p1Func3D.getMaxValue( theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "                     maximum   =  " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, 3.0 );

  measure = p1Func3D.getMinValue( theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "                     minimum   = " << std::scientific << measure );
  WALBERLA_CHECK_FLOAT_EQUAL( measure, -1.0 );

  runMagTest( "3D Test P2 function: magnitude =  ", theLevel, p2Func3D, testFuncCombo, 3.0 );

  return EXIT_SUCCESS;
}

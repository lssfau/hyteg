/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
// We test finding the maximum, maximum magnitude and minum of a P1 function.
// The corresponding value is selectively placed on a macro vertex, edge and
// face and finally a combined test is done
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "hyteg/facedofspace/FaceDoFFunction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::uint_t;

using namespace hyteg;

static real_t xLocPos = 0.0;
static real_t yLocPos = 0.0;
static real_t zLocPos = 0.0;
static real_t epsilon = 0.0;
static uint_t counter = 0;

#define TEST_MAX_VALUE  100.0
#define TEST_MIN_VALUE -200.0
#define TEST_MAG_VALUE  200.0

typedef enum{ FIND_MAX, FIND_MIN, FIND_MAG } findType;

// --------------------------------------------------------------------------------------------------
// Auxilliary function for executing the actual tests
// --------------------------------------------------------------------------------------------------
template<class funcType>
void runFindTest( std::string mesg, findType testType, uint_t theLevel, funcType &dofFunc,
                  std::function<real_t(const hyteg::Point3D&)> &testFunc, real_t refVal,
                  DoFType flag = All ) {

  real_t measure = 0.0;

  dofFunc.interpolate( testFunc, theLevel );

  switch( testType ) {

    case FIND_MAG:
      measure = dofFunc.getMaxMagnitude( theLevel, flag );
      break;

    case FIND_MAX:
      measure = dofFunc.getMaxValue( theLevel, flag );
      break;

    case FIND_MIN:
      measure = dofFunc.getMinValue( theLevel, flag );
      break;
  }

  WALBERLA_CHECK_FLOAT_EQUAL( measure, refVal );
  WALBERLA_LOG_INFO_ON_ROOT( mesg << std::scientific << std::showpos << measure );
}

// --------------------------------------------------------------------------------------------------

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

  std::function<real_t(const hyteg::Point3D&)> testFuncMax = []( const hyteg::Point3D& x ) {
    real_t distance = std::sqrt( (xLocPos - x[0]) * (xLocPos - x[0]) +
                                 (yLocPos - x[1]) * (yLocPos - x[1]) +
                                 (zLocPos - x[2]) * (zLocPos - x[2]) );
    return distance > epsilon ? real_t(0.0) : TEST_MAX_VALUE;
  };

  std::function<real_t(const hyteg::Point3D&)> testFuncMin = []( const hyteg::Point3D& x ) {
    real_t distance = std::sqrt( (xLocPos - x[0]) * (xLocPos - x[0]) +
                                 (yLocPos - x[1]) * (yLocPos - x[1]) +
                                 (zLocPos - x[2]) * (zLocPos - x[2]) );
    return distance > epsilon ? real_t(0.0) : TEST_MIN_VALUE;
  };

  std::function<real_t(const hyteg::Point3D&)> testFuncCombo = []( const hyteg::Point3D& x ) {
    WALBERLA_UNUSED(x);
    real_t retVal = 0.0;
    switch( counter%3 )
      {
      case 0: retVal = -5.0; break;
      case 1: retVal =  2.0; break;
      case 2: retVal =  3.0;
      }
    ++counter;
    return retVal;
  };

  std::function<real_t(const hyteg::Point3D&)> testFuncInvDst = []( const hyteg::Point3D& x ) {
    real_t distance = std::sqrt( (xLocPos - x[0]) * (xLocPos - x[0]) +
                                 (yLocPos - x[1]) * (yLocPos - x[1]) +
                                 (zLocPos - x[2]) * (zLocPos - x[2]) );
    return 2.0 - distance;
    // return std::exp( -distance );
  };

  real_t hmin = std::pow( 2.0, - (real_t)theLevel );
  epsilon = 0.5 * hmin;
  WALBERLA_LOG_INFO_ON_ROOT( "level    = " << theLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "h_min    = " << std::scientific << hmin    );
  WALBERLA_LOG_INFO_ON_ROOT( "epsilon  = " << std::scientific << epsilon );

  // =============
  //  P1Function 
  // =============
  hyteg::P1Function< real_t > p1Func2D( "", storage, theLevel, theLevel );

  WALBERLA_LOG_INFO_ON_ROOT( "\n\nP1Function (DoFType=All)\n" );

  // Special value on macro edge
  xLocPos = 0.50;
  yLocPos = 0.50;

  runFindTest( "Test #1 (edge     ): maximum   = ", FIND_MAX, theLevel, p1Func2D, testFuncMax, TEST_MAX_VALUE );
  runFindTest( "                     minumum   = ", FIND_MIN, theLevel, p1Func2D, testFuncMin, TEST_MIN_VALUE );
  runFindTest( "                     magnitude = ", FIND_MAG, theLevel, p1Func2D, testFuncMin, TEST_MAG_VALUE );

  // Special value on macro face
  xLocPos = 0.50;
  yLocPos = 0.25;

  runFindTest( "Test #2 (face     ): maximum   = ", FIND_MAX, theLevel, p1Func2D, testFuncMax, TEST_MAX_VALUE );
  runFindTest( "                     minumum   = ", FIND_MIN, theLevel, p1Func2D, testFuncMin, TEST_MIN_VALUE );
  runFindTest( "                     magnitude = ", FIND_MAG, theLevel, p1Func2D, testFuncMin, TEST_MAG_VALUE );

  // Special value on macro vertex
  xLocPos = 0.0;
  yLocPos = 0.0;

  runFindTest( "Test #3 (vertex   ): maximum   = ", FIND_MAX, theLevel, p1Func2D, testFuncMax, TEST_MAX_VALUE );
  runFindTest( "                     minumum   = ", FIND_MIN, theLevel, p1Func2D, testFuncMin, TEST_MIN_VALUE );
  runFindTest( "                     magnitude = ", FIND_MAG, theLevel, p1Func2D, testFuncMin, TEST_MAG_VALUE );

  // Combined test
  runFindTest( "Test #4 (combo    ): maximum   = ", FIND_MAX, theLevel, p1Func2D, testFuncMax, TEST_MAX_VALUE );
  runFindTest( "                     minumum   = ", FIND_MIN, theLevel, p1Func2D, testFuncMin, TEST_MIN_VALUE );
  runFindTest( "                     magnitude = ", FIND_MAG, theLevel, p1Func2D, testFuncMin, TEST_MAG_VALUE );

  WALBERLA_LOG_INFO_ON_ROOT( "\n\n P1Function (DoFType=<varying>)\n" );

  // DoFType test #1
  runFindTest( "Test #5 (combo    ): maximum   = ", FIND_MAX, theLevel, p1Func2D, testFuncCombo,  3.0, Inner );
  runFindTest( "                     minumum   = ", FIND_MIN, theLevel, p1Func2D, testFuncCombo, -5.0, Inner );
  runFindTest( "                     magnitude = ", FIND_MAG, theLevel, p1Func2D, testFuncCombo,  5.0, Inner );

  // DoFType test #2 (extremum on boundary edge)
  xLocPos = 0.00;
  yLocPos = 0.50;

  runFindTest( "Test #6 (bc edge  ): max outer = ", FIND_MAX, theLevel, p1Func2D, testFuncMax, TEST_MAX_VALUE, DirichletBoundary );
  runFindTest( "                     max inner = ", FIND_MAX, theLevel, p1Func2D, testFuncMax, 0.0, Inner );
  runFindTest( "                     min outer = ", FIND_MIN, theLevel, p1Func2D, testFuncMin, TEST_MIN_VALUE, DirichletBoundary );
  runFindTest( "                     min inner = ", FIND_MIN, theLevel, p1Func2D, testFuncMin, 0.0, Inner );
  runFindTest( "                     mag outer = ", FIND_MAG, theLevel, p1Func2D, testFuncMin, TEST_MAG_VALUE, DirichletBoundary );
  runFindTest( "                     mag inner = ", FIND_MAG, theLevel, p1Func2D, testFuncMin, 0.0, Inner );

  // DoFType test #3 (extremum on boundary vertex)
  xLocPos = 1.00;
  yLocPos = 1.00;

  runFindTest( "Test #7 (bc vertex): max outer = ", FIND_MAX, theLevel, p1Func2D, testFuncMax, TEST_MAX_VALUE, DirichletBoundary );
  runFindTest( "                     max inner = ", FIND_MAX, theLevel, p1Func2D, testFuncMax, 0.0, Inner );
  runFindTest( "                     min outer = ", FIND_MIN, theLevel, p1Func2D, testFuncMin, TEST_MIN_VALUE, DirichletBoundary );
  runFindTest( "                     min inner = ", FIND_MIN, theLevel, p1Func2D, testFuncMin, 0.0, Inner );
  runFindTest( "                     mag outer = ", FIND_MAG, theLevel, p1Func2D, testFuncMin, TEST_MAG_VALUE, DirichletBoundary );
  runFindTest( "                     mag inner = ", FIND_MAG, theLevel, p1Func2D, testFuncMin, 0.0, Inner );

  // Test case when DoFType is not present in base mesh
  runFindTest( "Test #8 (Neumann  ): maximum   = ", FIND_MAX, theLevel, p1Func2D, testFuncCombo, -std::numeric_limits<real_t>::max(), NeumannBoundary );
  runFindTest( "                     minimum   = ", FIND_MIN, theLevel, p1Func2D, testFuncCombo,  std::numeric_limits<real_t>::max(), NeumannBoundary );
  runFindTest( "                     magnitude = ", FIND_MAG, theLevel, p1Func2D, testFuncCombo,  0.0                               , NeumannBoundary );

  // ============
  //  P2Function 
  // ============

  WALBERLA_LOG_INFO_ON_ROOT( "\n\nP2Function (DoFType=All)\n" );

  theLevel = 2;
  hyteg::P2Function< real_t > p2func( "", storage, theLevel, theLevel );

  // Special value on macro edge (vertexdof)
  xLocPos = 0.25;
  yLocPos = 0.25;
  runFindTest( "Test #9 (edge/vert): maximum   = ", FIND_MAX, theLevel, p2func, testFuncMax, TEST_MAX_VALUE );
  runFindTest( "                     minimum   = ", FIND_MIN, theLevel, p2func, testFuncMin, TEST_MIN_VALUE );
  runFindTest( "                     magnitude = ", FIND_MAG, theLevel, p2func, testFuncMin, TEST_MAG_VALUE );

  // Special value on macro edge (edgedof)
  xLocPos = 0.125;
  yLocPos = 0.125;
  runFindTest( "Test #A (edge/edge): maximum   = ", FIND_MAX, theLevel, p2func, testFuncMax, TEST_MAX_VALUE );
  runFindTest( "                     minimum   = ", FIND_MIN, theLevel, p2func, testFuncMin, TEST_MIN_VALUE );
  runFindTest( "                     magnitude = ", FIND_MAG, theLevel, p2func, testFuncMin, TEST_MAG_VALUE );

  // Special value on macro face (edgedof)
  xLocPos = 0.250;
  yLocPos = 0.125;
  runFindTest( "Test #B (face/edge): maximum   = ", FIND_MAX, theLevel, p2func, testFuncMax, TEST_MAX_VALUE );
  runFindTest( "                     minimum   = ", FIND_MIN, theLevel, p2func, testFuncMin, TEST_MIN_VALUE );
  runFindTest( "                     magnitude = ", FIND_MAG, theLevel, p2func, testFuncMin, TEST_MAG_VALUE );

  // ============
  //  DGFunction 
  // ============

  WALBERLA_LOG_INFO_ON_ROOT( "\n\nDGFunction (DoFType=All)\n" );

  theLevel = 2;
  hyteg::FaceDoFFunction< real_t > dgFunc( "", storage, theLevel, theLevel );
  runFindTest( "Test #C (combo    ): maximum   = ", FIND_MAX, theLevel, dgFunc, testFuncCombo,  3.0 );
  runFindTest( "                     minimum   = ", FIND_MIN, theLevel, dgFunc, testFuncCombo, -5.0 );
  runFindTest( "                     magnitude = ", FIND_MAG, theLevel, dgFunc, testFuncCombo,  5.0 );

  // ===========
  //  3D Meshes
  // ===========

  WALBERLA_LOG_INFO_ON_ROOT( "\n\n------------------\n TESTs on 3D Mesh" );
  WALBERLA_LOG_INFO_ON_ROOT( "------------------\n" );

  theLevel = 2;

  // Generate mesh from meshfile (mesh represents unit cube)
  MeshInfo meshInfo3D = MeshInfo::emptyMeshInfo();
  meshInfo3D = MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_6el.msh" );

  // Generate primitives
  SetupPrimitiveStorage setupStorage3D( meshInfo3D,
                                        uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  loadbalancing::roundRobin( setupStorage3D );
  std::shared_ptr<PrimitiveStorage> storage3D = std::make_shared<PrimitiveStorage>( setupStorage3D );

  hyteg::P1Function< real_t > p1Func3D( "", storage3D, theLevel, theLevel );
  hyteg::P2Function< real_t > p2Func3D( "", storage3D, theLevel, theLevel );
  hyteg::FaceDoFFunction< real_t > dgFunc3D( "", storage3D, theLevel, theLevel );


  // --------------------------------------------------------------------------------------------------
  WALBERLA_LOG_INFO_ON_ROOT( "\nSingle point test (micro-vertex-dof)\n" );

  xLocPos = 0.25;
  yLocPos = 0.25;
  zLocPos = 0.25;

  runFindTest( "3D Test P1 function: magnitude = ", FIND_MAG, theLevel, p1Func3D, testFuncMin, TEST_MAG_VALUE );
  runFindTest( "                     maximum   = ", FIND_MAX, theLevel, p1Func3D, testFuncMax, TEST_MAX_VALUE );
  runFindTest( "                     minimum   = ", FIND_MIN, theLevel, p1Func3D, testFuncMin, TEST_MIN_VALUE );

  runFindTest( "3D Test P2 function: magnitude = ", FIND_MAG, theLevel, p2Func3D, testFuncMin, TEST_MAG_VALUE );
  runFindTest( "                     maximum   = ", FIND_MAX, theLevel, p2Func3D, testFuncMax, TEST_MAX_VALUE );
  runFindTest( "                     minimum   = ", FIND_MIN, theLevel, p2Func3D, testFuncMin, TEST_MIN_VALUE );

  // --------------------------------------------------------------------------------------------------
  WALBERLA_LOG_INFO_ON_ROOT( "\nSingle point test (micro-edge-dof #1)\n" );

  xLocPos = 0.125;
  yLocPos = 0.25;
  zLocPos = 0.25;

  runFindTest( "3D Test P2 function: magnitude = ", FIND_MAG, theLevel, p2Func3D, testFuncMin, TEST_MAG_VALUE );
  runFindTest( "                     maximum   = ", FIND_MAX, theLevel, p2Func3D, testFuncMax, TEST_MAX_VALUE );
  runFindTest( "                     minimum   = ", FIND_MIN, theLevel, p2Func3D, testFuncMin, TEST_MIN_VALUE );

  // --------------------------------------------------------------------------------------------------
  WALBERLA_LOG_INFO_ON_ROOT( "\nSingle point test (micro-edge-dof #2)\n" );

  xLocPos = 0.125;
  yLocPos = 0.125;
  zLocPos = 0.25;

  runFindTest( "3D Test P2 function: magnitude = ", FIND_MAG, theLevel, p2Func3D, testFuncMin, TEST_MAG_VALUE );
  runFindTest( "                     maximum   = ", FIND_MAX, theLevel, p2Func3D, testFuncMax, TEST_MAX_VALUE );
  runFindTest( "                     minimum   = ", FIND_MIN, theLevel, p2Func3D, testFuncMin, TEST_MIN_VALUE );

  // --------------------------------------------------------------------------------------------------
  WALBERLA_LOG_INFO_ON_ROOT( "\nSingle point test (micro-edge-dof #3)\n" );

  xLocPos = 0.125;
  yLocPos = 0.125;
  zLocPos = 0.125;

  runFindTest( "3D Test P2 function: magnitude = ", FIND_MAG, theLevel, p2Func3D, testFuncMin, TEST_MAG_VALUE );
  runFindTest( "                     maximum   = ", FIND_MAX, theLevel, p2Func3D, testFuncMax, TEST_MAX_VALUE );
  runFindTest( "                     minimum   = ", FIND_MIN, theLevel, p2Func3D, testFuncMin, TEST_MIN_VALUE );

  // --------------------------------------------------------------------------------------------------
  WALBERLA_LOG_INFO_ON_ROOT( "\nCombo test\n" );

  runFindTest( "3D Test P1 function: magnitude = ", FIND_MAG, theLevel, p1Func3D, testFuncCombo,  5.0 );
  runFindTest( "                     maximum   = ", FIND_MAX, theLevel, p1Func3D, testFuncCombo,  3.0 );
  runFindTest( "                     minimum   = ", FIND_MIN, theLevel, p1Func3D, testFuncCombo, -5.0 );

  runFindTest( "3D Test P2 function: magnitude = ", FIND_MAG, theLevel, p2Func3D, testFuncCombo,  5.0 );
  runFindTest( "                     maximum   = ", FIND_MAX, theLevel, p2Func3D, testFuncCombo,  3.0 );
  runFindTest( "                     minimum   = ", FIND_MIN, theLevel, p2Func3D, testFuncCombo, -5.0 );

  return EXIT_SUCCESS;
}

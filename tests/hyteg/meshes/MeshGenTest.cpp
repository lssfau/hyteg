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
// Test generation of MeshInfo object for inline mesh generators
//
// Simple test:
// - check that MeshInfo object can be generated
// - use it to setup a PrimitiveStorage
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "core/math/Constants.h"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::math::pi;
using walberla::uint_t;

typedef enum { RECTANGLE, ANNULUS, PARTIAL_ANNULUS, SPHERICAL_SHELL, FACE_CHAIN } testDomainType;

namespace hyteg {

static void testMeshGenerator( testDomainType testDomain, MeshInfo::meshFlavour flavour )
{

  std::shared_ptr< MeshInfo > meshInfo;

  // perform mesh generation
  switch( testDomain )
    {
    case RECTANGLE:
      meshInfo = std::make_shared< MeshInfo >( MeshInfo::meshRectangle( Point2D( {-2.0, 1.0} ), Point2D( {0.0, 3.0} ),
                                                        flavour, 3, 2 ) );
      break;

    case PARTIAL_ANNULUS:
      meshInfo = std::make_shared< MeshInfo >( MeshInfo::meshAnnulus( 1.0, 2.0, 0.25*pi, 0.75*pi, flavour, 4, 2 ) );
      break;

    case ANNULUS:
      meshInfo = std::make_shared< MeshInfo >( MeshInfo::meshAnnulus( 1.0, 2.0, flavour, 4, 2 ) );
      break;

    case SPHERICAL_SHELL:
      meshInfo = std::make_shared< MeshInfo >( MeshInfo::meshSphericalShell( 5, 3, 1.0, 2.0 ) );
      break;

    case FACE_CHAIN:
      meshInfo = std::make_shared< MeshInfo >( MeshInfo::meshFaceChain( 3 ) );

    }

  // test primitive counts
  uint_t numVerts = meshInfo->getVertices().size();
  uint_t numFaces = meshInfo->getFaces().size();
  uint_t numEdges = meshInfo->getEdges().size();
  uint_t numCells = meshInfo->getCells().size();

  if( testDomain == RECTANGLE )
    {
      switch( flavour )
       {
       case hyteg::MeshInfo::CRISS :
       case hyteg::MeshInfo::CROSS :
         WALBERLA_CHECK_EQUAL( numVerts, 12 );
         WALBERLA_CHECK_EQUAL( numFaces, 12 );
         WALBERLA_CHECK_EQUAL( numEdges, 23 );
         WALBERLA_CHECK_EQUAL( numCells,  0 );
         break;

       case hyteg::MeshInfo::CRISSCROSS :
       case hyteg::MeshInfo::DIAMOND :
         WALBERLA_CHECK_EQUAL( numVerts, 18 );
         WALBERLA_CHECK_EQUAL( numFaces, 24 );
         WALBERLA_CHECK_EQUAL( numEdges, 41 );
         WALBERLA_CHECK_EQUAL( numCells,  0 );
         break;
       }
    }
  else if( testDomain == FACE_CHAIN )
    {
      WALBERLA_CHECK_EQUAL( numVerts, 5 );
      WALBERLA_CHECK_EQUAL( numEdges, 7 );
      WALBERLA_CHECK_EQUAL( numFaces, 3 );
      WALBERLA_CHECK_EQUAL( numCells, 0 );
    }
  else if( testDomain == SPHERICAL_SHELL )
    {
      WALBERLA_CHECK_EQUAL( numVerts,  486 );
      WALBERLA_CHECK_EQUAL( numFaces, 4160 );
      WALBERLA_CHECK_EQUAL( numCells, 1920 );
    }

  // use generated MeshInfo
  SetupPrimitiveStorage setupStorage( *meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  loadbalancing::roundRobin( setupStorage );
  WALBERLA_LOG_INFO_ON_ROOT( setupStorage );
  PrimitiveStorage storage( setupStorage );

}

} // namespace hyteg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // Check mesh generator for rectangles
   hyteg::testMeshGenerator( RECTANGLE, hyteg::MeshInfo::CRISS );
   hyteg::testMeshGenerator( RECTANGLE, hyteg::MeshInfo::CROSS );
   hyteg::testMeshGenerator( RECTANGLE, hyteg::MeshInfo::CRISSCROSS );
   hyteg::testMeshGenerator( RECTANGLE, hyteg::MeshInfo::DIAMOND );

   // Check mesh generator for annuli
   hyteg::testMeshGenerator( ANNULUS, hyteg::MeshInfo::CROSS );
   hyteg::testMeshGenerator( ANNULUS, hyteg::MeshInfo::CRISS );
   hyteg::testMeshGenerator( ANNULUS, hyteg::MeshInfo::CRISSCROSS );
   hyteg::testMeshGenerator( ANNULUS, hyteg::MeshInfo::DIAMOND );
   hyteg::testMeshGenerator( PARTIAL_ANNULUS, hyteg::MeshInfo::CRISS );
   hyteg::testMeshGenerator( PARTIAL_ANNULUS, hyteg::MeshInfo::CROSS );
   hyteg::testMeshGenerator( PARTIAL_ANNULUS, hyteg::MeshInfo::CRISSCROSS );
   hyteg::testMeshGenerator( PARTIAL_ANNULUS, hyteg::MeshInfo::DIAMOND );

   // Check mesh generator for face chains
   // (note that flavour argument is meaningless here)
   hyteg::testMeshGenerator( FACE_CHAIN, hyteg::MeshInfo::DIAMOND );

   // Check mesh generator for thick spherical shell
   // (note that flavour argument is meaningless here)
   hyteg::testMeshGenerator( SPHERICAL_SHELL, hyteg::MeshInfo::DIAMOND );

   return EXIT_SUCCESS;
}

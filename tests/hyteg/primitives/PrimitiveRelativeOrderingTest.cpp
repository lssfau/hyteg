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

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/SendBuffer.h"
#include "core/mpi/RecvBuffer.h"
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

static void testPrimitiveRelativeOrdering()
{
  // Checking logical macro-primitive indexing according to documentation

  std::shared_ptr< PrimitiveStorage > storage = PrimitiveStorage::createFromGmshFile( "../../data/meshes/3D/tet_1el.msh" );

  WALBERLA_CHECK_EQUAL( storage->getNumberOfLocalCells(), 1 );
  WALBERLA_CHECK_EQUAL( storage->getNumberOfLocalFaces(), 4 );
  WALBERLA_CHECK_EQUAL( storage->getNumberOfLocalEdges(), 6 );
  WALBERLA_CHECK_EQUAL( storage->getNumberOfLocalVertices(), 4 );

  std::vector< PrimitiveID > cellIDs;
  std::vector< PrimitiveID > faceIDs;
  std::vector< PrimitiveID > edgeIDs;
  std::vector< PrimitiveID > vertexIDs;

  storage->getCellIDs( cellIDs );
  storage->getFaceIDs( faceIDs );
  storage->getEdgeIDs( edgeIDs );
  storage->getVertexIDs( vertexIDs );

  const auto cell = storage->getCell( cellIDs[0] );

  // All cell neighbors

  const auto face0 = storage->getFace( cell->neighborFaces().at( 0 ) );
  const auto face1 = storage->getFace( cell->neighborFaces().at( 1 ) );
  const auto face2 = storage->getFace( cell->neighborFaces().at( 2 ) );
  const auto face3 = storage->getFace( cell->neighborFaces().at( 3 ) );

  const auto edge0 = storage->getEdge( cell->neighborEdges().at( 0 ) );
  const auto edge1 = storage->getEdge( cell->neighborEdges().at( 1 ) );
  const auto edge2 = storage->getEdge( cell->neighborEdges().at( 2 ) );
  const auto edge3 = storage->getEdge( cell->neighborEdges().at( 3 ) );
  const auto edge4 = storage->getEdge( cell->neighborEdges().at( 4 ) );
  const auto edge5 = storage->getEdge( cell->neighborEdges().at( 5 ) );

  // All faces got only one neighbor cell which should be the only cell in our storage

  WALBERLA_CHECK_EQUAL( face0->getNumNeighborCells(), 1 );
  WALBERLA_CHECK_EQUAL( face1->getNumNeighborCells(), 1 );
  WALBERLA_CHECK_EQUAL( face2->getNumNeighborCells(), 1 );
  WALBERLA_CHECK_EQUAL( face3->getNumNeighborCells(), 1 );

  WALBERLA_CHECK_EQUAL( face0->neighborCells().at( 0 ), cell->getID() );
  WALBERLA_CHECK_EQUAL( face1->neighborCells().at( 0 ), cell->getID() );
  WALBERLA_CHECK_EQUAL( face2->neighborCells().at( 0 ), cell->getID() );
  WALBERLA_CHECK_EQUAL( face3->neighborCells().at( 0 ), cell->getID() );

  //////////////////////////////////////////////////////////////////////////////////
  // Checking that the cells neighbor faces are built from the correct vertices
  //

  std::set< PrimitiveID > cellLocalVertices;
  std::set< PrimitiveID > faceLocalVertices;

  // Cell-local face 0

  cellLocalVertices.insert( cell->neighborVertices().at( 0 ) );
  cellLocalVertices.insert( cell->neighborVertices().at( 1 ) );
  cellLocalVertices.insert( cell->neighborVertices().at( 2 ) );

  faceLocalVertices.insert( face0->neighborVertices().at( 0 ) );
  faceLocalVertices.insert( face0->neighborVertices().at( 1 ) );
  faceLocalVertices.insert( face0->neighborVertices().at( 2 ) );

  WALBERLA_CHECK( cellLocalVertices == faceLocalVertices, "Cell-local vertex IDs at face[0] are wrong!" );

  cellLocalVertices.clear();
  faceLocalVertices.clear();

  // Cell-local face 1

  cellLocalVertices.insert( cell->neighborVertices().at( 0 ) );
  cellLocalVertices.insert( cell->neighborVertices().at( 1 ) );
  cellLocalVertices.insert( cell->neighborVertices().at( 3 ) );

  faceLocalVertices.insert( face1->neighborVertices().at( 0 ) );
  faceLocalVertices.insert( face1->neighborVertices().at( 1 ) );
  faceLocalVertices.insert( face1->neighborVertices().at( 2 ) );

  WALBERLA_CHECK( cellLocalVertices == faceLocalVertices, "Cell-local vertex IDs at face[1] are wrong!" );

  cellLocalVertices.clear();
  faceLocalVertices.clear();

  // Cell-local face 2

  cellLocalVertices.insert( cell->neighborVertices().at( 0 ) );
  cellLocalVertices.insert( cell->neighborVertices().at( 2 ) );
  cellLocalVertices.insert( cell->neighborVertices().at( 3 ) );

  faceLocalVertices.insert( face2->neighborVertices().at( 0 ) );
  faceLocalVertices.insert( face2->neighborVertices().at( 1 ) );
  faceLocalVertices.insert( face2->neighborVertices().at( 2 ) );

  WALBERLA_CHECK( cellLocalVertices == faceLocalVertices, "Cell-local vertex IDs at face[2] are wrong!" );

  cellLocalVertices.clear();
  faceLocalVertices.clear();

  // Cell-local face 3

  cellLocalVertices.insert( cell->neighborVertices().at( 1 ) );
  cellLocalVertices.insert( cell->neighborVertices().at( 2 ) );
  cellLocalVertices.insert( cell->neighborVertices().at( 3 ) );

  faceLocalVertices.insert( face3->neighborVertices().at( 0 ) );
  faceLocalVertices.insert( face3->neighborVertices().at( 1 ) );
  faceLocalVertices.insert( face3->neighborVertices().at( 2 ) );

  WALBERLA_CHECK( cellLocalVertices == faceLocalVertices, "Cell-local vertex IDs at face[3] are wrong!" );

  cellLocalVertices.clear();
  faceLocalVertices.clear();

  //////////////////////////////////////////////////////////////////////////////////
  // Checking that the cells neighbor edges are built from the correct vertices
  //

  std::set< PrimitiveID > edgeLocalVertices;

  // Cell-local edge 0

  cellLocalVertices.insert( cell->neighborVertices().at( 0 ) );
  cellLocalVertices.insert( cell->neighborVertices().at( 1 ) );

  edgeLocalVertices.insert( edge0->neighborVertices().at( 0 ) );
  edgeLocalVertices.insert( edge0->neighborVertices().at( 1 ) );

  WALBERLA_CHECK( cellLocalVertices == edgeLocalVertices, "Cell-local vertex IDs at edge[0] are wrong!" );

  cellLocalVertices.clear();
  edgeLocalVertices.clear();

  // Cell-local edge 1

  cellLocalVertices.insert( cell->neighborVertices().at( 0 ) );
  cellLocalVertices.insert( cell->neighborVertices().at( 2 ) );

  edgeLocalVertices.insert( edge1->neighborVertices().at( 0 ) );
  edgeLocalVertices.insert( edge1->neighborVertices().at( 1 ) );

  WALBERLA_CHECK( cellLocalVertices == edgeLocalVertices, "Cell-local vertex IDs at edge[1] are wrong!" );

  cellLocalVertices.clear();
  edgeLocalVertices.clear();

  // Cell-local edge 2

  cellLocalVertices.insert( cell->neighborVertices().at( 1 ) );
  cellLocalVertices.insert( cell->neighborVertices().at( 2 ) );

  edgeLocalVertices.insert( edge2->neighborVertices().at( 0 ) );
  edgeLocalVertices.insert( edge2->neighborVertices().at( 1 ) );

  WALBERLA_CHECK( cellLocalVertices == edgeLocalVertices, "Cell-local vertex IDs at edge[2] are wrong!" );

  cellLocalVertices.clear();
  edgeLocalVertices.clear();

  // Cell-local edge 3

  cellLocalVertices.insert( cell->neighborVertices().at( 0 ) );
  cellLocalVertices.insert( cell->neighborVertices().at( 3 ) );

  edgeLocalVertices.insert( edge3->neighborVertices().at( 0 ) );
  edgeLocalVertices.insert( edge3->neighborVertices().at( 1 ) );

  WALBERLA_CHECK( cellLocalVertices == edgeLocalVertices, "Cell-local vertex IDs at edge[3] are wrong!" );

  cellLocalVertices.clear();
  edgeLocalVertices.clear();

  // Cell-local edge 4

  cellLocalVertices.insert( cell->neighborVertices().at( 1 ) );
  cellLocalVertices.insert( cell->neighborVertices().at( 3 ) );

  edgeLocalVertices.insert( edge4->neighborVertices().at( 0 ) );
  edgeLocalVertices.insert( edge4->neighborVertices().at( 1 ) );

  WALBERLA_CHECK( cellLocalVertices == edgeLocalVertices, "Cell-local vertex IDs at edge[4] are wrong!" );

  cellLocalVertices.clear();
  edgeLocalVertices.clear();

  // Cell-local edge 5

  cellLocalVertices.insert( cell->neighborVertices().at( 2 ) );
  cellLocalVertices.insert( cell->neighborVertices().at( 3 ) );

  edgeLocalVertices.insert( edge5->neighborVertices().at( 0 ) );
  edgeLocalVertices.insert( edge5->neighborVertices().at( 1 ) );

  WALBERLA_CHECK( cellLocalVertices == edgeLocalVertices, "Cell-local vertex IDs at edge[5] are wrong!" );

  cellLocalVertices.clear();
  edgeLocalVertices.clear();

}

} // namespace hyteg


int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();
   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testPrimitiveRelativeOrdering();

   return EXIT_SUCCESS;
}

/*
 * Copyright (c) 2017-2019 Nils Kohl.
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

#include <numeric>

#include "core/DataTypes.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;

namespace hyteg {

void P2Level0IndexingTest()
{
   const uint_t level = 0;

   WALBERLA_CHECK_EQUAL( levelinfo::num_microvertices_per_edge( level ), 2 );
   WALBERLA_CHECK_EQUAL( levelinfo::num_microedges_per_edge( level ), 1 );

   WALBERLA_CHECK_EQUAL( edgedof::macroedge::index( level, 0 ), 0 );
   WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexOnNeighborFace( level, 0, 0, edgedof::EdgeDoFOrientation::XY ), 1 );
   WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexOnNeighborFace( level, 0, 0, edgedof::EdgeDoFOrientation::Y ), 2 );

   WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexOnNeighborFace( level, 0, 1, edgedof::EdgeDoFOrientation::XY ), 3 );
   WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexOnNeighborFace( level, 0, 1, edgedof::EdgeDoFOrientation::Y ), 4 );

   WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexOnNeighborCell( level, 0, 0, 4, edgedof::EdgeDoFOrientation::YZ ), 1 + 2 * 4 );
   WALBERLA_CHECK_EQUAL( edgedof::macroedge::indexOnNeighborCell( level, 0, 1, 4, edgedof::EdgeDoFOrientation::YZ ),
                         1 + 2 * 4 + 1 );
}

void P2Level0InterpolateTest()
{
   const uint_t level = 0;

   MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/regular_octahedron_8el.msh" );
   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   P2Function< real_t > f_interpolation( "f_interpolation", storage, level, level );

   for ( const auto& itEdge : storage->getEdges() )
   {
      auto edge   = itEdge.second;

      auto macroEdgeEdgeDoFDataSize = edge->getData( f_interpolation.getEdgeDoFFunction().getEdgeDataID() )->getSize( level );

      WALBERLA_CHECK_EQUAL( macroEdgeEdgeDoFDataSize, 1 + 2 * edge->getNumNeighborFaces() + 1 * edge->getNumNeighborCells() );
   }

   auto someFunction = []( const Point3D& x ) -> real_t { return real_c( 42.0 + x[0] * 2367. + x[1] * 37. + x[2] * 999. ); };
   f_interpolation.interpolate( someFunction, level, All );

   communication::syncP2FunctionBetweenPrimitives( f_interpolation, level );

   for ( const auto& itVertex : storage->getVertices() )
   {
      auto vertexID = itVertex.first;
      auto vertex   = itVertex.second;

      auto macroVertexVertexDoFDataSize =
          vertex->getData( f_interpolation.getVertexDoFFunction().getVertexDataID() )->getSize( level );
      auto macroVertexVertexDoFData =
          vertex->getData( f_interpolation.getVertexDoFFunction().getVertexDataID() )->getPointer( level );

      auto macroVertexEdgeDoFData =
          vertex->getData( f_interpolation.getEdgeDoFFunction().getVertexDataID() )->getPointer( level );

      WALBERLA_CHECK_EQUAL( macroVertexVertexDoFDataSize, vertex->getNumNeighborEdges() + 1 );

      // neighboring edges
      for ( const auto& edgeID : vertex->neighborEdges() )
      {
         auto edge = storage->getEdge( edgeID );

         auto macroEdgeVertexDoFData =
             edge->getData( f_interpolation.getVertexDoFFunction().getEdgeDataID() )->getPointer( level );

         auto macroEdgeEdgeDoFData = edge->getData( f_interpolation.getEdgeDoFFunction().getEdgeDataID() )->getPointer( level );

         auto vertexLocalEdgeID = vertex->edge_index( edgeID );
         auto edgeLocalVertexID = edge->vertex_index( vertexID );

         if ( edgeLocalVertexID == 0 )
         {
            // vertex DoF data
            WALBERLA_CHECK_FLOAT_EQUAL( macroEdgeVertexDoFData[0], macroVertexVertexDoFData[0] );
            WALBERLA_CHECK_FLOAT_EQUAL( macroEdgeVertexDoFData[1], macroVertexVertexDoFData[vertexLocalEdgeID + 1] );
         }
         else if ( edgeLocalVertexID == 1 )
         {
            // vertex DoF data
            WALBERLA_CHECK_FLOAT_EQUAL( macroEdgeVertexDoFData[1], macroVertexVertexDoFData[0] );
            WALBERLA_CHECK_FLOAT_EQUAL( macroEdgeVertexDoFData[0], macroVertexVertexDoFData[vertexLocalEdgeID + 1] );
         }

         // edge DoF data
         WALBERLA_CHECK_FLOAT_EQUAL( macroEdgeEdgeDoFData[0], macroVertexEdgeDoFData[vertexLocalEdgeID] );
      }

      // opposite edge of vertex in face
      for ( const auto& faceID : vertex->neighborFaces() )
      {
         auto face         = storage->getFace( faceID );
         auto oppositeEdge = storage->getEdge( face->getEdgeOppositeToVertex( vertexID ) );

         auto macroEdgeVertexDoFData =
             oppositeEdge->getData( f_interpolation.getVertexDoFFunction().getEdgeDataID() )->getPointer( level );

         auto macroEdgeEdgeDoFData =
             oppositeEdge->getData( f_interpolation.getEdgeDoFFunction().getEdgeDataID() )->getPointer( level );

         // check that macro-vertex ghost layers are correct
         WALBERLA_CHECK_FLOAT_EQUAL( macroVertexEdgeDoFData[vertex->getNumNeighborEdges() + vertex->face_index( faceID )],
                                     macroEdgeEdgeDoFData[0] );

         // check that macro-edge ghost layers are correct (vertexdofs)
         WALBERLA_CHECK_FLOAT_EQUAL( macroEdgeVertexDoFData[2 + oppositeEdge->face_index( faceID )],
                                     macroVertexVertexDoFData[0] );
      }
   }

   for ( const auto& itEdge : storage->getEdges() )
   {
      auto edgeID = itEdge.first;
      auto edge   = itEdge.second;

      auto macroEdgeEdgeDoFDataSize = edge->getData( f_interpolation.getEdgeDoFFunction().getEdgeDataID() )->getSize( level );
      auto macroEdgeEdgeDoFData     = edge->getData( f_interpolation.getEdgeDoFFunction().getEdgeDataID() )->getPointer( level );

      WALBERLA_CHECK_EQUAL( macroEdgeEdgeDoFDataSize, 1 + 2 * edge->getNumNeighborFaces() + 1 * edge->getNumNeighborCells() );

      for ( const auto & faceID : edge->neighborFaces() )
      {
         auto face = storage->getFace( faceID );

         // check macro-edge ghost layers are correct (edgedofs)
         auto glValue0 = macroEdgeEdgeDoFData[edgedof::macroedge::indexOnNeighborFace(
             level, 0, edge->face_index( faceID ), edgedof::EdgeDoFOrientation::Y )];
         auto glValue1 = macroEdgeEdgeDoFData[edgedof::macroedge::indexOnNeighborFace(
             level, 0, edge->face_index( faceID ), edgedof::EdgeDoFOrientation::XY )];

         // get correct edge
         auto neighborEdge0 = storage->getEdge( face->getEdgeOppositeToVertex( edge->neighborVertices().at( 1 ) ) );
         auto neighborEdge1 = storage->getEdge( face->getEdgeOppositeToVertex( edge->neighborVertices().at( 0 ) ) );

         auto realValue0 = neighborEdge0->getData( f_interpolation.getEdgeDoFFunction().getEdgeDataID() )->getPointer( level )[0];
         auto realValue1 = neighborEdge1->getData( f_interpolation.getEdgeDoFFunction().getEdgeDataID() )->getPointer( level )[0];

         WALBERLA_CHECK_FLOAT_EQUAL( glValue0, realValue0 );
         WALBERLA_CHECK_FLOAT_EQUAL( glValue1, realValue1 );
      }

      for ( const auto & cellID : edge->neighborCells() )
      {
         // now check opposite edge
         auto cell         = storage->getCell( cellID );
         auto oppositeEdge = storage->getEdge( cell->getOppositeEdgeID( edgeID ) );
         auto oppsiteEdgeValue =
             oppositeEdge->getData( f_interpolation.getEdgeDoFFunction().getEdgeDataID() )->getPointer( level )[0];

         auto glValue = macroEdgeEdgeDoFData[edgedof::macroedge::indexOnNeighborCell(
             level, 0, edge->cell_index( cellID ), edge->getNumNeighborFaces(), edgedof::EdgeDoFOrientation::YZ )];

         WALBERLA_CHECK_FLOAT_EQUAL( oppsiteEdgeValue, glValue );

         WALBERLA_LOG_INFO_ON_ROOT(
             "macro-edge: "
             << edgeID << ", opposite edgedof value (in cell): "
             << macroEdgeEdgeDoFData[edgedof::macroedge::indexOnNeighborCell(
                    level, 0, edge->cell_index( cellID ), edge->getNumNeighborFaces(), edgedof::EdgeDoFOrientation::YZ )] )
      }
   }
}

void P2Level0EnumerateTet1elTest()
{
  auto numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );
  auto rank = uint_c( walberla::mpi::MPIManager::instance()->rank() );

  const uint_t level = 0;

  MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" );
  SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
  std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

  P2Function< int > f_interpolation( "f_interpolation", storage, level, level );
  f_interpolation.enumerate( level );

  communication::syncP2FunctionBetweenPrimitives( f_interpolation, level );

  std::map< PrimitiveID, int > realDataLocal;
  std::map< PrimitiveID, int > realDataGlobal;

  // fill real data
  for ( const auto& itVertex : storage->getVertices() )
  {
    auto vertexID = itVertex.first;
    auto vertex = itVertex.second;

    auto macroVertexVertexDoFData =
    vertex->getData( f_interpolation.getVertexDoFFunction().getVertexDataID())->getPointer( level );

    realDataLocal[vertexID] = macroVertexVertexDoFData[0];
  }

  for ( const auto& itEdge : storage->getEdges() )
  {
    auto edgeID = itEdge.first;
    auto edge = itEdge.second;

    auto macroEdgeEdgeDoFData =
    edge->getData( f_interpolation.getEdgeDoFFunction().getEdgeDataID())->getPointer( level );

    realDataLocal[edgeID] = macroEdgeEdgeDoFData[0];
  }

  walberla::mpi::SendBuffer sendBuffer;
  walberla::mpi::RecvBuffer recvBuffer;

  for ( uint_t p = 0; p < numProcesses; p++ )
  {
    if ( p == rank )
    {
      WALBERLA_LOG_INFO( "real data (local):" )
      for ( auto d : realDataLocal )
      {
        WALBERLA_LOG_INFO( "ID " << d.first << ": " << d.second );
        sendBuffer << d.first;
        sendBuffer << d.second;
      }
    }
    WALBERLA_MPI_BARRIER();
  }

  walberla::mpi::allGathervBuffer( sendBuffer, recvBuffer );

  while ( !recvBuffer.isEmpty() )
  {
    PrimitiveID primitiveID;
    int value;

    recvBuffer >> primitiveID;
    recvBuffer >> value;

    realDataGlobal[primitiveID] = value;
  }


  for ( uint_t p = 0; p < numProcesses; p++ )
  {
    if ( p == rank )
    {
      WALBERLA_LOG_INFO( "real data (global):" )
      for ( auto d : realDataGlobal )
      {
        WALBERLA_LOG_INFO( "ID " << d.first << ": " << d.second );
      }
    }
    WALBERLA_MPI_BARRIER();
  }

  WALBERLA_CHECK_EQUAL( realDataGlobal.size(), 10 );

  // fill ghost data maps with all macro-vertex data
  for ( const auto& itVertex : storage->getVertices() )
  {
    std::map< PrimitiveID, int > glData;

    auto vertexID = itVertex.first;
    auto vertex   = itVertex.second;

    WALBERLA_LOG_INFO( "macro-vertex (ID " << vertexID << ") ghost data check..." )

    auto macroVertexVertexDoFDataSize =
    vertex->getData( f_interpolation.getVertexDoFFunction().getVertexDataID() )->getSize( level );
    auto macroVertexVertexDoFData =
    vertex->getData( f_interpolation.getVertexDoFFunction().getVertexDataID() )->getPointer( level );

    auto macroVertexEdgeDoFData =
    vertex->getData( f_interpolation.getEdgeDoFFunction().getVertexDataID() )->getPointer( level );

    WALBERLA_CHECK_EQUAL( macroVertexVertexDoFDataSize, vertex->getNumNeighborEdges() + 1 );

    // "neighboring" vertices and neighboring edges
    for ( const auto& edgeID : vertex->neighborEdges() )
    {
      auto edge = storage->getEdge( edgeID );
      WALBERLA_CHECK_EQUAL( realDataGlobal[edge->get_opposite_vertex( vertexID )], macroVertexVertexDoFData[1 + vertex->edge_index( edgeID )] );
      WALBERLA_CHECK_EQUAL( realDataGlobal[edgeID], macroVertexEdgeDoFData[vertex->edge_index( edgeID )] );
    }

    // opposite edges
    for ( const auto& faceID : vertex->neighborFaces() )
    {
      auto face         = storage->getFace( faceID );
      auto oppositeEdge = storage->getEdge( face->getEdgeOppositeToVertex( vertexID ) );

      WALBERLA_CHECK_NOT_NULLPTR( oppositeEdge, "This test cannot be executed with this number of processes." )

      auto checkLHS = realDataGlobal[oppositeEdge->getID()];
      auto checkRHS = macroVertexEdgeDoFData[vertex->getNumNeighborEdges() + vertex->face_index( faceID )];

      WALBERLA_CHECK_EQUAL( checkLHS, checkRHS );
    }
  }

  WALBERLA_MPI_BARRIER();

  // fill ghost data maps with all macro-edge data
  for ( const auto& itEdge : storage->getEdges() )
  {
    std::map< PrimitiveID, int > glData;

    auto edgeID = itEdge.first;
    auto edge = itEdge.second;

    WALBERLA_LOG_INFO( "macro-edge (ID " << edgeID << ") ghost data check..." )

    auto macroEdgeVertexDoFData =
    edge->getData( f_interpolation.getVertexDoFFunction().getEdgeDataID())->getPointer( level );

    // vertex dofs "on" edge
    WALBERLA_CHECK_EQUAL( realDataGlobal[edge->neighborVertices().at(0)], macroEdgeVertexDoFData[0]);
    WALBERLA_CHECK_EQUAL( realDataGlobal[edge->neighborVertices().at(1)], macroEdgeVertexDoFData[1]);

    // vertex dofs "on" faces
    WALBERLA_CHECK_EQUAL( realDataGlobal[storage->getFace( edge->neighborFaces().at(0) )->get_vertex_opposite_to_edge(edgeID)], macroEdgeVertexDoFData[2] );
    WALBERLA_CHECK_EQUAL( realDataGlobal[storage->getFace( edge->neighborFaces().at(1) )->get_vertex_opposite_to_edge(edgeID)], macroEdgeVertexDoFData[3] );

    auto macroEdgeEdgeDoFData =
    edge->getData( f_interpolation.getEdgeDoFFunction().getEdgeDataID())->getPointer( level );

    // edge dofs "on" faces
    WALBERLA_CHECK_EQUAL( realDataGlobal[storage->getFace( edge->neighborFaces().at(0) )->getEdgeOppositeToVertex(edge->neighborVertices().at(0))], macroEdgeEdgeDoFData[1] ); // XY first
    WALBERLA_CHECK_EQUAL( realDataGlobal[storage->getFace( edge->neighborFaces().at(0) )->getEdgeOppositeToVertex(edge->neighborVertices().at(1))], macroEdgeEdgeDoFData[2] ); // Y second

    WALBERLA_CHECK_EQUAL( realDataGlobal[storage->getFace( edge->neighborFaces().at(1) )->getEdgeOppositeToVertex(edge->neighborVertices().at(0))], macroEdgeEdgeDoFData[3] ); // XY first
    WALBERLA_CHECK_EQUAL( realDataGlobal[storage->getFace( edge->neighborFaces().at(1) )->getEdgeOppositeToVertex(edge->neighborVertices().at(1))], macroEdgeEdgeDoFData[4] ); // Y second

    // edge dofs "in" cells
    WALBERLA_CHECK_EQUAL( realDataGlobal[storage->getCell( edge->neighborCells().at(0) )->getOppositeEdgeID(edgeID)], macroEdgeEdgeDoFData[5] );
  }


}


} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   if ( walberla::mpi::MPIManager::instance()->numProcesses() == 1 )
   {
     hyteg::P2Level0IndexingTest();
     hyteg::P2Level0InterpolateTest();
   }
   hyteg::P2Level0EnumerateTet1elTest();


   return EXIT_SUCCESS;
}

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
#include <hyteg/communication/Syncing.hpp>

#include "core/debug/all.h"
#include "core/mpi/all.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/facedofspace_old/FaceDoFFunction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using namespace hyteg;

using hyteg::vertexdof::VertexDoFFunction;
using walberla::real_t;

void checkComm( const std::string& meshfile, const uint_t maxLevel, bool bufferComm = false )
{
   //MeshInfo meshInfo = MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh");
   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( meshfile );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const uint_t minLevel = 2;
   //const uint_t maxLevel = 4;
   hyteg::FaceDoFFunction_old< int32_t > x( "x", storage, minLevel, maxLevel );
   if ( bufferComm )
   {
      x.setLocalCommunicationMode( communication::BufferedCommunicator::BUFFERED_MPI );
   }

   int32_t num = 1;
   x.enumerate( maxLevel, num );

   uint_t numberOfChecks      = 0;
   uint_t totalExpectedChecks = 0;

   for ( auto& edgeIt : storage->getEdges() )
   {
      if ( edgeIt.second.get()->getNumHigherDimNeighbors() == 1 )
      {
         totalExpectedChecks += 4;
      }
      else if ( edgeIt.second.get()->getNumHigherDimNeighbors() == 2 )
      {
         totalExpectedChecks += 8;
      }
      else
      {
         WALBERLA_CHECK( false )
      }
   }

   //// check vertex to edge comm ////
   for ( auto& edgeIt : storage->getEdges() )
   {
      Edge& edge = *edgeIt.second;
      //BubbleEdge::printFunctionMemory(edge,x.getEdgeDataID(),maxLevel);
      int32_t*                   edgeData = edge.getData( x.getEdgeDataID() )->getPointer( maxLevel );
      std::vector< PrimitiveID > nbrVertices;
      edge.getNeighborVertices( nbrVertices );
      for ( auto& vertexIt : nbrVertices )
      {
         Vertex*  vertex     = storage->getVertex( vertexIt.getID() );
         int32_t* vertexData = vertex->getData( x.getVertexDataID() )->getPointer( maxLevel );
         uint_t   vPerEdge   = levelinfo::num_microvertices_per_edge( maxLevel );
         uint_t   pos        = std::numeric_limits< uint_t >::max();
         if ( edge.vertex_index( vertex->getID() ) == 0 )
         {
            pos = 0;
         }
         else if ( edge.vertex_index( vertex->getID() ) == 1 )
         {
            pos = vPerEdge - 2;
         }
         else
         {
            WALBERLA_CHECK( false, "vertex not on Edge" )
         }
         uint_t index = facedof::macroedge::indexFaceFromVertex( maxLevel, pos, stencilDirection::CELL_GRAY_SE );
         WALBERLA_CHECK_UNEQUAL( 0, edgeData[index] )
         WALBERLA_CHECK_EQUAL( edgeData[index], vertexData[vertex->face_index( edge.neighborFaces()[0] ) * 2] )
         index = facedof::macroedge::indexFaceFromVertex( maxLevel, pos == 0 ? pos + 1 : pos, stencilDirection::CELL_BLUE_SE );
         WALBERLA_CHECK_UNEQUAL( 0, edgeData[index] )
         WALBERLA_CHECK_EQUAL( edgeData[index], vertexData[vertex->face_index( edge.neighborFaces()[0] ) * 2 + 1] )
         numberOfChecks += 2;
         if ( edge.getNumNeighborFaces() == 2 )
         {
            index = facedof::macroedge::indexFaceFromVertex( maxLevel, pos, stencilDirection::CELL_GRAY_NE );
            WALBERLA_CHECK_UNEQUAL( 0, edgeData[index] )
            WALBERLA_CHECK_EQUAL( edgeData[index], vertexData[vertex->face_index( edge.neighborFaces()[1] ) * 2] )
            index = facedof::macroedge::indexFaceFromVertex( maxLevel, pos == 0 ? pos + 1 : pos, stencilDirection::CELL_BLUE_NW );
            WALBERLA_CHECK_UNEQUAL( 0, edgeData[index] )
            WALBERLA_CHECK_EQUAL( edgeData[index], vertexData[vertex->face_index( edge.neighborFaces()[1] ) * 2 + 1] )
            numberOfChecks += 2;
         }
      }
   }

   WALBERLA_CHECK_EQUAL( totalExpectedChecks, numberOfChecks )

   /// check face edge comms ///
   numberOfChecks = 0;
   totalExpectedChecks =
       ( 2 * hyteg::levelinfo::num_microvertices_per_edge( maxLevel ) - 3 ) * 3 * storage->getNumberOfLocalFaces();

   for ( auto& faceIt : storage->getFaces() )
   {
      Face&                      face     = *faceIt.second;
      int32_t*                   faceData = face.getData( x.getFaceDataID() )->getPointer( maxLevel );
      std::vector< PrimitiveID > nbrEdges;
      face.getNeighborEdges( nbrEdges );
      for ( uint_t i = 0; i < nbrEdges.size(); ++i )
      {
         Edge*    edge         = storage->getEdge( nbrEdges[0].getID() );
         int32_t* edgeData     = edge->getData( x.getEdgeDataID() )->getPointer( maxLevel );
         uint_t   idxCounter   = 0;
         uint_t   faceIdOnEdge = edge->face_index( face.getID() );
         //////////////////// GRAY CELL //////////////////////
         idxCounter = 0;
         auto it    = facedof::macroface::indexIterator( face.edge_index( edge->getID() ),
                                                      face.getEdgeOrientation()[face.edge_index( edge->getID() )],
                                                      facedof::macroface::CELL_GRAY,
                                                      maxLevel );
         for ( ; it != facedof::macroface::indexIterator(); ++it )
         {
            if ( faceIdOnEdge == 0 )
            {
               WALBERLA_CHECK_UNEQUAL( 0, faceData[*it] )
               WALBERLA_CHECK_EQUAL(
                   edgeData[facedof::macroedge::indexFaceFromVertex( maxLevel, idxCounter, stencilDirection::CELL_GRAY_SE )],
                   faceData[*it] )
               numberOfChecks++;
            }
            else if ( faceIdOnEdge == 1 )
            {
               WALBERLA_CHECK_UNEQUAL( 0, faceData[*it] )
               WALBERLA_CHECK_EQUAL(
                   edgeData[facedof::macroedge::indexFaceFromVertex( maxLevel, idxCounter, stencilDirection::CELL_GRAY_NE )],
                   faceData[*it] )
               numberOfChecks++;
            }
            else
            {
               WALBERLA_CHECK( false )
            }
            idxCounter++;
         }
         //////////////////// BLUE CELL //////////////////////
         idxCounter = 0;
         it         = facedof::macroface::indexIterator( face.edge_index( edge->getID() ),
                                                 face.getEdgeOrientation()[face.edge_index( edge->getID() )],
                                                 facedof::macroface::CELL_BLUE,
                                                 maxLevel );
         for ( ; it != facedof::macroface::indexIterator(); ++it )
         {
            if ( faceIdOnEdge == 0 )
            {
               WALBERLA_CHECK_EQUAL(
                   edgeData[facedof::macroedge::indexFaceFromVertex( maxLevel, idxCounter + 1, stencilDirection::CELL_BLUE_SE )],
                   faceData[*it] )
               numberOfChecks++;
            }
            else if ( faceIdOnEdge == 1 )
            {
               WALBERLA_CHECK_EQUAL(
                   edgeData[facedof::macroedge::indexFaceFromVertex( maxLevel, idxCounter + 1, stencilDirection::CELL_BLUE_NW )],
                   faceData[*it] )
               numberOfChecks++;
            }
            else
            {
               WALBERLA_CHECK( false )
            }
            idxCounter++;
         }
      }
   }
   WALBERLA_CHECK_EQUAL( totalExpectedChecks, numberOfChecks )
}

void printCellGhostlayer3D( const std::string&                   message,
                            VertexDoFFunction< real_t >&         func,
                            std::shared_ptr< PrimitiveStorage >& storage,
                            const uint_t                         level )
{
   WALBERLA_LOG_INFO_ON_ROOT( message );
   for ( const auto& pid : storage->getCellIDs() )
   {
      const auto cell = storage->getCell( pid );
      //  WALBERLA_LOG_INFO_ON_ROOT( "Iterating cell " << pid );
      for ( const auto& [n, npid] : cell->getIndirectNeighborCellIDsOverFaces() )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Iterating neighbor cell " << n );
         const auto GLVector = cell->getData( func.getCellGLDataID( n ) )->getVector( level );
         for ( const auto& v : GLVector )
         {
            std::cout << v << " ";
         }
         std::cout << std::endl;
      }
   }
}
void printCellGhostlayer2D( const std::string&                   message,
                            VertexDoFFunction< real_t >&         func,
                            std::shared_ptr< PrimitiveStorage >& storage,
                            const uint_t                         level )
{
   WALBERLA_LOG_INFO_ON_ROOT( message );
   for ( const auto& pid : storage->getFaceIDs() )
   {
      const auto face = storage->getFace( pid );
      //  WALBERLA_LOG_INFO_ON_ROOT( "Iterating cell " << pid );
      for ( const auto& [n, npid] : face->getIndirectNeighborFaceIDsOverEdges() )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Iterating neighbor face " << n );
         const auto GLVector = face->getData( func.getFaceGLDataID( n ) )->getVector( level );
         for ( const auto& v : GLVector )
         {
            std::cout << v << " ";
         }
         std::cout << std::endl;
      }
   }
}

void checkVertexDoFCellComm3D( const std::string& meshfile, const uint_t level )
{
   //MeshInfo meshInfo = MeshInfo::fromGmshFile("../../data/meshes/quad_4el.msh");
   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( meshfile );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   VertexDoFFunction< real_t > testFunc( "testFunc", storage, level, level, BoundaryCondition::create0123BC(), 1 );
   testFunc.interpolate( 1, level, All );
   communication::syncFunctionBetweenPrimitives( testFunc, level );

   testFunc.template communicate< Cell, Cell >( level );
   printCellGhostlayer3D( "After communicate:", testFunc, storage, level );
   VTKOutput vtk( "../../output", "DGCommTest_3D", storage );
   vtk.add( testFunc );

   vtk.write( level );
}
void checkVertexDoFCellComm2D( const std::string& meshfile, const uint_t level )
{
   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( meshfile );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   VertexDoFFunction< real_t > testFunc( "testFunc", storage, level, level, BoundaryCondition::create0123BC(), 1 );
   testFunc.interpolate( 1, level, All );
   communication::syncFunctionBetweenPrimitives( testFunc, level );
   testFunc.template communicate< Face, Face >( level );
   printCellGhostlayer2D( "After communicate:", testFunc, storage, level );
   VTKOutput vtk( "../../output", "DGCommTest_2D", storage );
   vtk.add( testFunc );

   vtk.write( level );
}
int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   //  walberla::debug::enterTestMode();
   checkVertexDoFCellComm2D( "../../data/meshes/tri_2el.msh", 1 );
   checkVertexDoFCellComm2D( "../../data/meshes/tri_2el.msh", 2 );

   checkVertexDoFCellComm2D( "../../data/meshes/tri_2el.msh", 3 );

   //checkVertexDoFCellComm3D( "../../data/meshes/3D/pyramid_2el.msh", 1 );
   checkVertexDoFCellComm3D( "../../data/meshes/3D/pyramid_2el.msh", 3 );

   /*
   checkComm( "../../data/meshes/quad_4el.msh", 4, true );

   checkComm( "../../data/meshes/quad_4el.msh", 5, true );

   checkComm( "../../data/meshes/quad_4el.msh", 4, false );

   checkComm( "../../data/meshes/quad_4el.msh", 5, false );

   checkComm( "../../data/meshes/bfs_12el.msh", 3, true );

   checkComm( "../../data/meshes/bfs_12el.msh", 3, false );
*/
   return 0;
}

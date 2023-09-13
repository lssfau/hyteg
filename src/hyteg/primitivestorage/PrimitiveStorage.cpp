/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

#include <algorithm>
#include <map>
#include <vector>

#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/math/DistributedSample.h"
#include "core/mpi/Gatherv.h"
#include "core/mpi/OpenMPBufferSystem.h"

#include "hyteg/Algorithms.hpp"
#include "hyteg/communication/PackageBufferSystem.hpp"
#include "hyteg/primitivedata/PrimitiveDataID.hpp"
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/primitives/Edge.hpp"
#include "hyteg/primitives/Face.hpp"
#include "hyteg/primitives/Primitive.hpp"
#include "hyteg/primitives/Vertex.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

namespace hyteg {

using walberla::uint_t;

std::shared_ptr< PrimitiveStorage > PrimitiveStorage::createFromGmshFile( const std::string& meshFilePath )
{
   const MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFilePath );
   const SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   return std::make_shared< PrimitiveStorage >( setupStorage );
}

PrimitiveStorage::PrimitiveStorage( const SetupPrimitiveStorage&                     setupStorage,
                                    const std::shared_ptr< walberla::WcTimingTree >& timingTree,
                                    const uint_t&                                    additionalHaloDepth )
: primitiveDataHandlers_( 0 )
, modificationStamp_( 0 )
, timingTree_( timingTree )
, hasGlobalCells_( setupStorage.getNumberOfCells() > 0 )
, additionalHaloDepth_( additionalHaloDepth )
{
   // We need to construct at least the maps on the coarsest level.
   vertices_[0];
   edges_[0];
   faces_[0];
   cells_[0];
   neighborVertices_[0];
   neighborEdges_[0];
   neighborFaces_[0];
   neighborCells_[0];

   setupStorage.broadcastPrimitives( vertices_[0],
                                     edges_[0],
                                     faces_[0],
                                     cells_[0],
                                     neighborVertices_[0],
                                     neighborEdges_[0],
                                     neighborFaces_[0],
                                     neighborCells_[0],
                                     neighborRanks_[0]);

//
//   for ( auto it : setupStorage.getVertices() )
//   {
//      if ( uint_c( walberla::mpi::MPIManager::instance()->rank() ) == setupStorage.getTargetRank( it.first ) )
//      {
//         vertices_[0][it.first] = std::make_shared< Vertex >( *( it.second ) );
//      }
//   }
//
//   for ( auto it : setupStorage.getEdges() )
//   {
//      if ( uint_c( walberla::mpi::MPIManager::instance()->rank() ) == setupStorage.getTargetRank( it.first ) )
//      {
//         edges_[0][it.first] = std::make_shared< Edge >( *( it.second ) );
//      }
//   }
//
//   for ( auto it : setupStorage.getFaces() )
//   {
//      if ( uint_c( walberla::mpi::MPIManager::instance()->rank() ) == setupStorage.getTargetRank( it.first ) )
//      {
//         faces_[0][it.first] = std::make_shared< Face >( *( it.second ) );
//      }
//   }
//
//   for ( auto it : setupStorage.getCells() )
//   {
//      if ( uint_c( walberla::mpi::MPIManager::instance()->rank() ) == setupStorage.getTargetRank( it.first ) )
//      {
//         cells_[0][it.first] = std::make_shared< Cell >( *( it.second ) );
//      }
//   }
//
//   // neighbors
//   std::vector< PrimitiveID > vertices{ getVertexIDs() };
//   std::vector< PrimitiveID > edges{ getEdgeIDs() };
//   std::vector< PrimitiveID > faces{ getFaceIDs() };
//   std::vector< PrimitiveID > cells{ getCellIDs() };
//   addDirectNeighbors( setupStorage, vertices, edges, faces, cells );

   // additionally requested neighbors
   for ( uint_t k = 0; k < additionalHaloDepth; k += 1 )
   {
      std::vector< PrimitiveID > additionalVertices;
      getNeighboringVertexIDs( additionalVertices );
      std::vector< PrimitiveID > additionalEdges;
      getNeighboringEdgeIDs( additionalEdges );
      std::vector< PrimitiveID > additionalFaces;
      getNeighboringFaceIDs( additionalFaces );
      std::vector< PrimitiveID > additionalCells;
      getNeighboringCellIDs( additionalCells );
      addDirectNeighbors( setupStorage, additionalVertices, additionalEdges, additionalFaces, additionalCells );
   }

   splitCommunicatorByPrimitiveDistribution();
   updateLeafPrimitiveMaps();

#ifndef NDEBUG
   checkConsistency();
#endif
}

void PrimitiveStorage::addDirectNeighbors( const SetupPrimitiveStorage&      setupStorage,
                                           const std::vector< PrimitiveID >& vertices,
                                           const std::vector< PrimitiveID >& edges,
                                           const std::vector< PrimitiveID >& faces,
                                           const std::vector< PrimitiveID >& cells )
{
   for ( const auto& id : vertices )
   {
      const Vertex* vertex = getVertex( id );
      WALBERLA_ASSERT_NOT_NULLPTR( vertex );

      for ( const auto& neighborVertexID : vertex->neighborVertices() )
      {
         const Vertex* neighborVertex = setupStorage.getVertex( neighborVertexID );
         if ( !vertexExistsLocally( neighborVertexID ) && !vertexExistsInNeighborhood( neighborVertexID ) )
         {
            neighborVertices_[0][neighborVertexID] = std::make_shared< Vertex >( *neighborVertex );
            neighborRanks_[0][neighborVertexID]    = setupStorage.getTargetRank( neighborVertexID );
         }
      }

      for ( const auto& neighborEdgeID : vertex->neighborEdges() )
      {
         const Edge* neighborEdge = setupStorage.getEdge( neighborEdgeID );
         if ( !edgeExistsLocally( neighborEdgeID ) && !edgeExistsInNeighborhood( neighborEdgeID ) )
         {
            neighborEdges_[0][neighborEdgeID] = std::make_shared< Edge >( *neighborEdge );
            neighborRanks_[0][neighborEdgeID] = setupStorage.getTargetRank( neighborEdgeID );
         }
      }

      for ( const auto& neighborFaceID : vertex->neighborFaces() )
      {
         const Face* neighborFace = setupStorage.getFace( neighborFaceID );
         if ( !faceExistsLocally( neighborFaceID ) && !faceExistsInNeighborhood( neighborFaceID ) )
         {
            neighborFaces_[0][neighborFaceID] = std::make_shared< Face >( *neighborFace );
            neighborRanks_[0][neighborFaceID] = setupStorage.getTargetRank( neighborFaceID );
         }
      }

      for ( const auto& neighborCellID : vertex->neighborCells() )
      {
         const Cell* neighborCell = setupStorage.getCell( neighborCellID );
         if ( !cellExistsLocally( neighborCellID ) && !cellExistsInNeighborhood( neighborCellID ) )
         {
            neighborCells_[0][neighborCellID] = std::make_shared< Cell >( *neighborCell );
            neighborRanks_[0][neighborCellID] = setupStorage.getTargetRank( neighborCellID );
         }
      }
   }

   for ( const auto& id : edges )
   {
      const Edge* edge = getEdge( id );
      WALBERLA_ASSERT_NOT_NULLPTR( edge );

      for ( const auto& neighborVertexID : edge->neighborVertices() )
      {
         const Vertex* neighborVertex = setupStorage.getVertex( neighborVertexID );
         if ( !vertexExistsLocally( neighborVertexID ) && !vertexExistsInNeighborhood( neighborVertexID ) )
         {
            neighborVertices_[0][neighborVertexID] = std::make_shared< Vertex >( *neighborVertex );
            neighborRanks_[0][neighborVertexID]    = setupStorage.getTargetRank( neighborVertexID );
         }
      }

      for ( const auto& neighborEdgeID : edge->neighborEdges() )
      {
         const Edge* neighborEdge = setupStorage.getEdge( neighborEdgeID );
         if ( !edgeExistsLocally( neighborEdgeID ) && !edgeExistsInNeighborhood( neighborEdgeID ) )
         {
            neighborEdges_[0][neighborEdgeID] = std::make_shared< Edge >( *neighborEdge );
            neighborRanks_[0][neighborEdgeID] = setupStorage.getTargetRank( neighborEdgeID );
         }
      }

      for ( const auto& neighborFaceID : edge->neighborFaces() )
      {
         const Face* neighborFace = setupStorage.getFace( neighborFaceID );
         if ( !faceExistsLocally( neighborFaceID ) && !faceExistsInNeighborhood( neighborFaceID ) )
         {
            neighborFaces_[0][neighborFaceID] = std::make_shared< Face >( *neighborFace );
            neighborRanks_[0][neighborFaceID] = setupStorage.getTargetRank( neighborFaceID );
         }
      }

      for ( const auto& neighborCellID : edge->neighborCells() )
      {
         const Cell* neighborCell = setupStorage.getCell( neighborCellID );
         if ( !cellExistsLocally( neighborCellID ) && !cellExistsInNeighborhood( neighborCellID ) )
         {
            neighborCells_[0][neighborCellID] = std::make_shared< Cell >( *neighborCell );
            neighborRanks_[0][neighborCellID] = setupStorage.getTargetRank( neighborCellID );
         }
      }
   }

   for ( const auto& id : faces )
   {
      const Face* face = getFace( id );
      WALBERLA_ASSERT_NOT_NULLPTR( face );

      for ( const auto& neighborVertexID : face->neighborVertices() )
      {
         const Vertex* neighborVertex = setupStorage.getVertex( neighborVertexID );
         if ( !vertexExistsLocally( neighborVertexID ) && !vertexExistsInNeighborhood( neighborVertexID ) )
         {
            neighborVertices_[0][neighborVertexID] = std::make_shared< Vertex >( *neighborVertex );
            neighborRanks_[0][neighborVertexID]    = setupStorage.getTargetRank( neighborVertexID );
         }
      }

      for ( const auto& neighborEdgeID : face->neighborEdges() )
      {
         const Edge* neighborEdge = setupStorage.getEdge( neighborEdgeID );
         if ( !edgeExistsLocally( neighborEdgeID ) && !edgeExistsInNeighborhood( neighborEdgeID ) )
         {
            neighborEdges_[0][neighborEdgeID] = std::make_shared< Edge >( *neighborEdge );
            neighborRanks_[0][neighborEdgeID] = setupStorage.getTargetRank( neighborEdgeID );
         }
      }

      for ( const auto& neighborFaceID : face->neighborFaces() )
      {
         const Face* neighborFace = setupStorage.getFace( neighborFaceID );
         if ( !faceExistsLocally( neighborFaceID ) && !faceExistsInNeighborhood( neighborFaceID ) )
         {
            neighborFaces_[0][neighborFaceID] = std::make_shared< Face >( *neighborFace );
            neighborRanks_[0][neighborFaceID] = setupStorage.getTargetRank( neighborFaceID );
         }
      }

      for ( const auto& neighborCellID : face->neighborCells() )
      {
         const Cell* neighborCell = setupStorage.getCell( neighborCellID );
         if ( !cellExistsLocally( neighborCellID ) && !cellExistsInNeighborhood( neighborCellID ) )
         {
            neighborCells_[0][neighborCellID] = std::make_shared< Cell >( *neighborCell );
            neighborRanks_[0][neighborCellID] = setupStorage.getTargetRank( neighborCellID );
         }
      }
   }

   for ( const auto& id : cells )
   {
      const Cell* cell = getCell( id );
      WALBERLA_ASSERT_NOT_NULLPTR( cell );

      for ( const auto& neighborVertexID : cell->neighborVertices() )
      {
         const Vertex* neighborVertex = setupStorage.getVertex( neighborVertexID );
         if ( !vertexExistsLocally( neighborVertexID ) && !vertexExistsInNeighborhood( neighborVertexID ) )
         {
            neighborVertices_[0][neighborVertexID] = std::make_shared< Vertex >( *neighborVertex );
            neighborRanks_[0][neighborVertexID]    = setupStorage.getTargetRank( neighborVertexID );
         }
      }

      for ( const auto& neighborEdgeID : cell->neighborEdges() )
      {
         const Edge* neighborEdge = setupStorage.getEdge( neighborEdgeID );
         if ( !edgeExistsLocally( neighborEdgeID ) && !edgeExistsInNeighborhood( neighborEdgeID ) )
         {
            neighborEdges_[0][neighborEdgeID] = std::make_shared< Edge >( *neighborEdge );
            neighborRanks_[0][neighborEdgeID] = setupStorage.getTargetRank( neighborEdgeID );
         }
      }

      for ( const auto& neighborFaceID : cell->neighborFaces() )
      {
         const Face* neighborFace = setupStorage.getFace( neighborFaceID );
         if ( !faceExistsLocally( neighborFaceID ) && !faceExistsInNeighborhood( neighborFaceID ) )
         {
            neighborFaces_[0][neighborFaceID] = std::make_shared< Face >( *neighborFace );
            neighborRanks_[0][neighborFaceID] = setupStorage.getTargetRank( neighborFaceID );
         }
      }

      for ( const auto& neighborCellID : cell->neighborCells() )
      {
         const Cell* neighborCell = setupStorage.getCell( neighborCellID );
         if ( !cellExistsLocally( neighborCellID ) && !cellExistsInNeighborhood( neighborCellID ) )
         {
            neighborCells_[0][neighborCellID] = std::make_shared< Cell >( *neighborCell );
            neighborRanks_[0][neighborCellID] = setupStorage.getTargetRank( neighborCellID );
         }
      }
   }

   splitCommunicatorByPrimitiveDistribution();

#ifndef NDEBUG
   checkConsistency();
#endif
}

void PrimitiveStorage::addDirectNeighborsDistributed()
{
   // Let's make it simple: send all local primitives to all neighbor ranks.
   // Receive from all neighbor ranks all their local primitives and store ranks as well.
   // Now we could drop information we do not need, but we simply keep it.

   walberla::mpi::BufferSystem bs( walberla::mpi::MPIManager::instance()->comm() );
   std::set< MPIRank >         nranks;
   for ( const auto& it : getNeighboringRanks() )
   {
      nranks.insert( static_cast< MPIRank >( it ) );
   }

   bs.setReceiverInfo( nranks, true );

   for ( auto nbrank : getNeighboringRanks() )
   {
      bs.sendBuffer( nbrank ) << getNumberOfLocalVertices();
      for ( const auto& it : getVertices() )
      {
         bs.sendBuffer( nbrank ) << it.first;
         bs.sendBuffer( nbrank ) << *it.second;
      }

      bs.sendBuffer( nbrank ) << getNumberOfLocalEdges();
      for ( const auto& it : getEdges() )
      {
         bs.sendBuffer( nbrank ) << it.first;
         bs.sendBuffer( nbrank ) << *it.second;
      }

      bs.sendBuffer( nbrank ) << getNumberOfLocalFaces();
      for ( const auto& it : getFaces() )
      {
         bs.sendBuffer( nbrank ) << it.first;
         bs.sendBuffer( nbrank ) << *it.second;
      }

      bs.sendBuffer( nbrank ) << getNumberOfLocalCells();
      for ( const auto& it : getCells() )
      {
         bs.sendBuffer( nbrank ) << it.first;
         bs.sendBuffer( nbrank ) << *it.second;
      }
   }

   bs.sendAll();

   for ( auto msg = bs.begin(); msg != bs.end(); ++msg )
   {
      const auto nbrank = msg.rank();

      uint_t numVertices;
      msg.buffer() >> numVertices;
      for ( uint_t i = 0; i < numVertices; i++ )
      {
         PrimitiveID id;
         msg.buffer() >> id;
         neighborVertices_[0][id] = std::make_shared< Vertex >( msg.buffer() );
         neighborRanks_[0][id]    = uint_c( nbrank );
      }

      uint_t numEdges;
      msg.buffer() >> numEdges;
      for ( uint_t i = 0; i < numEdges; i++ )
      {
         PrimitiveID id;
         msg.buffer() >> id;
         neighborEdges_[0][id] = std::make_shared< Edge >( msg.buffer() );
         neighborRanks_[0][id] = uint_c( nbrank );
      }

      uint_t numFaces;
      msg.buffer() >> numFaces;
      for ( uint_t i = 0; i < numFaces; i++ )
      {
         PrimitiveID id;
         msg.buffer() >> id;
         neighborFaces_[0][id] = std::make_shared< Face >( msg.buffer() );
         neighborRanks_[0][id] = uint_c( nbrank );
      }

      uint_t numCells;
      msg.buffer() >> numCells;
      for ( uint_t i = 0; i < numCells; i++ )
      {
         PrimitiveID id;
         msg.buffer() >> id;
         neighborCells_[0][id] = std::make_shared< Cell >( msg.buffer() );
         neighborRanks_[0][id] = uint_c( nbrank );
      }
   }
}

PrimitiveStorage::PrimitiveStorage( const SetupPrimitiveStorage& setupStorage, const uint_t& additionalHaloDepth )
: PrimitiveStorage( setupStorage, std::make_shared< walberla::WcTimingTree >(), additionalHaloDepth )
{}

PrimitiveStorage::PrimitiveStorage( const VertexMap&      vtxs,
                                    const EdgeMap&        edges,
                                    const FaceMap&        faces,
                                    const CellMap&        cells,
                                    const VertexMap&      nbrvtxs,
                                    const EdgeMap&        nbredges,
                                    const FaceMap&        nbrfaces,
                                    const CellMap&        nbrcells,
                                    const MigrationMap_T& neighborRanks,
                                    const bool&           hasGlobalCells )
: primitiveDataHandlers_( 0 )
, modificationStamp_( 0 )
, timingTree_( std::make_shared< walberla::WcTimingTree >() )
, hasGlobalCells_( hasGlobalCells )
, additionalHaloDepth_( 0 )
{
   vertices_[0] = vtxs;
   edges_[0]    = edges;
   faces_[0]    = faces;
   cells_[0]    = cells;

   neighborVertices_[0] = nbrvtxs;
   neighborEdges_[0]    = nbredges;
   neighborFaces_[0]    = nbrfaces;
   neighborCells_[0]    = nbrcells;

   neighborRanks_[0] = neighborRanks;

   splitCommunicatorByPrimitiveDistribution();
   updateLeafPrimitiveMaps();

#ifndef NDEBUG
   checkConsistency();
#endif
}

uint_t PrimitiveStorage::getCurrentLocalMaxRefinement() const
{
   // We assume here that all primitive maps have the same refinements (which they should have!).
   return vertices_.size() - 1;
}

/// \brief Returns the currently maximum number that a primitive is refined globally (involves global reduction).
uint_t PrimitiveStorage::getCurrentGlobalMaxRefinement() const
{
   return walberla::mpi::allReduce( getCurrentLocalMaxRefinement(), walberla::mpi::MAX );
}

/// \brief Returns the refinement level of the corresponding primitive.
uint_t PrimitiveStorage::getRefinementLevel( const PrimitiveID& pid ) const
{
   return pid.numAncestors();
}

uint_t PrimitiveStorage::getNumberOfLocalVertices() const
{
   uint_t num = 0;
   for ( const auto& [level, primitives] : vertices_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            num++;
         }
         WALBERLA_UNUSED( pid );
      }
      WALBERLA_UNUSED( level );
   }
   return num;
}

uint_t PrimitiveStorage::getNumberOfLocalEdges() const
{
   uint_t num = 0;
   for ( const auto& [level, primitives] : edges_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            num++;
         }
         WALBERLA_UNUSED( pid );
      }
      WALBERLA_UNUSED( level );
   }
   return num;
}

uint_t PrimitiveStorage::getNumberOfLocalFaces() const
{
   uint_t num = 0;
   for ( const auto& [level, primitives] : faces_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            num++;
         }
         WALBERLA_UNUSED( pid );
      }
      WALBERLA_UNUSED( level );
   }
   return num;
}

uint_t PrimitiveStorage::getNumberOfLocalCells() const
{
   uint_t num = 0;
   for ( const auto& [level, primitives] : cells_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            num++;
         }
         WALBERLA_UNUSED( pid );
      }
      WALBERLA_UNUSED( level );
   }
   return num;
}

bool PrimitiveStorage::vertexExistsLocally( const PrimitiveID& id ) const
{
   for ( const auto& [level, primitives] : vertices_ )
   {
      if ( primitives.count( id ) > 0 )
      {
         return true;
      }
      WALBERLA_UNUSED( level );
   }
   return false;
}

bool PrimitiveStorage::edgeExistsLocally( const PrimitiveID& id ) const
{
   for ( const auto& [level, primitives] : edges_ )
   {
      if ( primitives.count( id ) > 0 )
      {
         return true;
      }
      WALBERLA_UNUSED( level );
   }
   return false;
}

bool PrimitiveStorage::faceExistsLocally( const PrimitiveID& id ) const
{
   for ( const auto& [level, primitives] : faces_ )
   {
      if ( primitives.count( id ) > 0 )
      {
         return true;
      }
      WALBERLA_UNUSED( level );
   }
   return false;
}

bool PrimitiveStorage::cellExistsLocally( const PrimitiveID& id ) const
{
   for ( const auto& [level, primitives] : cells_ )
   {
      if ( primitives.count( id ) > 0 )
      {
         return true;
      }
      WALBERLA_UNUSED( level );
   }
   return false;
}

bool PrimitiveStorage::vertexExistsInNeighborhood( const PrimitiveID& id ) const
{
   for ( const auto& [level, primitives] : neighborVertices_ )
   {
      if ( primitives.count( id ) > 0 )
      {
         return true;
      }
      WALBERLA_UNUSED( level );
   }
   return false;
}

bool PrimitiveStorage::edgeExistsInNeighborhood( const PrimitiveID& id ) const
{
   for ( const auto& [level, primitives] : neighborEdges_ )
   {
      if ( primitives.count( id ) > 0 )
      {
         return true;
      }
      WALBERLA_UNUSED( level );
   }
   return false;
}

bool PrimitiveStorage::faceExistsInNeighborhood( const PrimitiveID& id ) const
{
   for ( const auto& [level, primitives] : neighborFaces_ )
   {
      if ( primitives.count( id ) > 0 )
      {
         return true;
      }
      WALBERLA_UNUSED( level );
   }
   return false;
}

bool PrimitiveStorage::cellExistsInNeighborhood( const PrimitiveID& id ) const
{
   for ( const auto& [level, primitives] : neighborCells_ )
   {
      if ( primitives.count( id ) > 0 )
      {
         return true;
      }
      WALBERLA_UNUSED( level );
   }
   return false;
}

std::shared_ptr< PrimitiveStorage > PrimitiveStorage::createCopy() const
{
   WALBERLA_CHECK_EQUAL( getCurrentLocalMaxRefinement(), 0, "Copy not implemented for refined meshes." );

   auto copiedStorage = std::make_shared< PrimitiveStorage >(
       SetupPrimitiveStorage( MeshInfo::emptyMeshInfo(), uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) ) );

   for ( const auto& it : getVertices() )
   {
      copiedStorage->vertices_[0][it.first] = std::make_shared< Vertex >( *it.second );
   }

   for ( const auto& it : getEdges() )
   {
      copiedStorage->edges_[0][it.first] = std::make_shared< Edge >( *it.second );
   }

   for ( const auto& it : getFaces() )
   {
      copiedStorage->faces_[0][it.first] = std::make_shared< Face >( *it.second );
   }

   for ( const auto& it : getCells() )
   {
      copiedStorage->cells_[0][it.first] = std::make_shared< Cell >( *it.second );
   }

   for ( const auto& it : neighborVertices_.at( 0 ) )
   {
      copiedStorage->neighborVertices_[0][it.first] = std::make_shared< Vertex >( *it.second );
   }

   for ( const auto& it : neighborEdges_.at( 0 ) )
   {
      copiedStorage->neighborEdges_[0][it.first] = std::make_shared< Edge >( *it.second );
   }

   for ( const auto& it : neighborFaces_.at( 0 ) )
   {
      copiedStorage->neighborFaces_[0][it.first] = std::make_shared< Face >( *it.second );
   }

   for ( const auto& it : neighborCells_.at( 0 ) )
   {
      copiedStorage->neighborCells_[0][it.first] = std::make_shared< Cell >( *it.second );
   }

   copiedStorage->neighborRanks_  = neighborRanks_;
   copiedStorage->hasGlobalCells_ = hasGlobalCells_;
   copiedStorage->splitComm_      = splitComm_;
   copiedStorage->updateLeafPrimitiveMaps();

   return copiedStorage;
}

void PrimitiveStorage::getPrimitives( PrimitiveMap& primitiveMap ) const
{
   primitiveMap.clear();

   auto vertices = getVertices();
   auto edges    = getEdges();
   auto faces    = getFaces();
   auto cells    = getCells();

   primitiveMap.insert( vertices.begin(), vertices.end() );
   primitiveMap.insert( edges.begin(), edges.end() );
   primitiveMap.insert( faces.begin(), faces.end() );
   primitiveMap.insert( cells.begin(), cells.end() );

   WALBERLA_ASSERT_EQUAL( primitiveMap.size(), vertices.size() + edges.size() + faces.size() + cells.size() );
}

/// returns all vertices without any children
PrimitiveStorage::VertexMap PrimitiveStorage::getVertices() const
{
   return leafVertices_;
}

/// returns all edges without any children
PrimitiveStorage::EdgeMap PrimitiveStorage::getEdges() const
{
   return leafEdges_;
}

/// returns all faces without any children
PrimitiveStorage::FaceMap PrimitiveStorage::getFaces() const
{
   return leafFaces_;
}

/// returns all cells without any children
PrimitiveStorage::CellMap PrimitiveStorage::getCells() const
{
   return leafCells_;
}

PrimitiveStorage::VertexMap PrimitiveStorage::getNeighborVertices() const
{
   PrimitiveStorage::VertexMap pmap;
   for ( const auto& [level, primitives] : neighborVertices_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            pmap[pid] = primitive;
         }
      }
      WALBERLA_UNUSED( level );
   }
   return pmap;
}

PrimitiveStorage::EdgeMap PrimitiveStorage::getNeighborEdges() const
{
   PrimitiveStorage::EdgeMap pmap;
   for ( const auto& [level, primitives] : neighborEdges_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            pmap[pid] = primitive;
         }
      }
      WALBERLA_UNUSED( level );
   }
   return pmap;
}

PrimitiveStorage::FaceMap PrimitiveStorage::getNeighborFaces() const
{
   PrimitiveStorage::FaceMap pmap;
   for ( const auto& [level, primitives] : neighborFaces_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            pmap[pid] = primitive;
         }
      }
      WALBERLA_UNUSED( level );
   }
   return pmap;
}

PrimitiveStorage::CellMap PrimitiveStorage::getNeighborCells() const
{
   PrimitiveStorage::CellMap pmap;
   for ( const auto& [level, primitives] : neighborCells_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            pmap[pid] = primitive;
         }
      }
      WALBERLA_UNUSED( level );
   }
   return pmap;
}

const Primitive* PrimitiveStorage::getPrimitive( const PrimitiveID& id ) const
{
   if ( vertexExistsLocally( id ) || vertexExistsInNeighborhood( id ) )
      return getVertex( id );
   if ( edgeExistsLocally( id ) || edgeExistsInNeighborhood( id ) )
      return getEdge( id );
   if ( faceExistsLocally( id ) || faceExistsInNeighborhood( id ) )
      return getFace( id );
   if ( cellExistsLocally( id ) || cellExistsInNeighborhood( id ) )
      return getCell( id );
   return nullptr;
}

Primitive* PrimitiveStorage::getPrimitive( const PrimitiveID& id )
{
   if ( vertexExistsLocally( id ) || vertexExistsInNeighborhood( id ) )
      return getVertex( id );
   if ( edgeExistsLocally( id ) || edgeExistsInNeighborhood( id ) )
      return getEdge( id );
   if ( faceExistsLocally( id ) || faceExistsInNeighborhood( id ) )
      return getFace( id );
   if ( cellExistsLocally( id ) || cellExistsInNeighborhood( id ) )
      return getCell( id );
   return nullptr;
}

const Vertex* PrimitiveStorage::getVertex( const PrimitiveID& id ) const
{
   if ( vertexExistsLocally( id ) )
   {
      for ( const auto& [level, primitives] : vertices_ )
      {
         if ( primitives.count( id ) > 0 )
         {
            return primitives.at( id ).get();
         }
         WALBERLA_UNUSED( level );
      }
   }
   else if ( vertexExistsInNeighborhood( id ) )
   {
      for ( const auto& [level, primitives] : neighborVertices_ )
      {
         if ( primitives.count( id ) > 0 )
         {
            return primitives.at( id ).get();
         }
         WALBERLA_UNUSED( level );
      }
   }

   return nullptr;
}

Vertex* PrimitiveStorage::getVertex( const PrimitiveID& id )
{
   if ( vertexExistsLocally( id ) )
   {
      for ( auto& [level, primitives] : vertices_ )
      {
         if ( primitives.count( id ) > 0 )
         {
            return primitives[id].get();
         }
         WALBERLA_UNUSED( level );
      }
   }
   else if ( vertexExistsInNeighborhood( id ) )
   {
      for ( auto& [level, primitives] : neighborVertices_ )
      {
         if ( primitives.count( id ) > 0 )
         {
            return primitives[id].get();
         }
         WALBERLA_UNUSED( level );
      }
   }

   return nullptr;
}

const Edge* PrimitiveStorage::getEdge( const PrimitiveID& id ) const
{
   if ( edgeExistsLocally( id ) )
   {
      for ( const auto& [level, primitives] : edges_ )
      {
         if ( primitives.count( id ) > 0 )
         {
            return primitives.at( id ).get();
         }
         WALBERLA_UNUSED( level );
      }
   }
   else if ( edgeExistsInNeighborhood( id ) )
   {
      for ( const auto& [level, primitives] : neighborEdges_ )
      {
         if ( primitives.count( id ) > 0 )
         {
            return primitives.at( id ).get();
         }
         WALBERLA_UNUSED( level );
      }
   }

   return nullptr;
}

Edge* PrimitiveStorage::getEdge( const PrimitiveID& id )
{
   if ( edgeExistsLocally( id ) )
   {
      for ( auto& [level, primitives] : edges_ )
      {
         if ( primitives.count( id ) > 0 )
         {
            return primitives[id].get();
         }
         WALBERLA_UNUSED( level );
      }
   }
   else if ( edgeExistsInNeighborhood( id ) )
   {
      for ( auto& [level, primitives] : neighborEdges_ )
      {
         if ( primitives.count( id ) > 0 )
         {
            return primitives[id].get();
         }
         WALBERLA_UNUSED( level );
      }
   }

   return nullptr;
}

const Face* PrimitiveStorage::getFace( const PrimitiveID& id ) const
{
   if ( faceExistsLocally( id ) )
   {
      for ( const auto& [level, primitives] : faces_ )
      {
         if ( primitives.count( id ) > 0 )
         {
            return primitives.at( id ).get();
         }
         WALBERLA_UNUSED( level );
      }
   }
   else if ( faceExistsInNeighborhood( id ) )
   {
      for ( const auto& [level, primitives] : neighborFaces_ )
      {
         if ( primitives.count( id ) > 0 )
         {
            return primitives.at( id ).get();
         }
         WALBERLA_UNUSED( level );
      }
   }

   return nullptr;
}

Face* PrimitiveStorage::getFace( const PrimitiveID& id )
{
   if ( faceExistsLocally( id ) )
   {
      for ( auto& [level, primitives] : faces_ )
      {
         if ( primitives.count( id ) > 0 )
         {
            return primitives[id].get();
         }
         WALBERLA_UNUSED( level );
      }
   }
   else if ( faceExistsInNeighborhood( id ) )
   {
      for ( auto& [level, primitives] : neighborFaces_ )
      {
         if ( primitives.count( id ) > 0 )
         {
            return primitives[id].get();
         }
         WALBERLA_UNUSED( level );
      }
   }

   return nullptr;
}

const Cell* PrimitiveStorage::getCell( const PrimitiveID& id ) const
{
   if ( cellExistsLocally( id ) )
   {
      for ( const auto& [level, primitives] : cells_ )
      {
         if ( primitives.count( id ) > 0 )
         {
            return primitives.at( id ).get();
         }
         WALBERLA_UNUSED( level );
      }
   }
   else if ( cellExistsInNeighborhood( id ) )
   {
      for ( const auto& [level, primitives] : neighborCells_ )
      {
         if ( primitives.count( id ) > 0 )
         {
            return primitives.at( id ).get();
         }
         WALBERLA_UNUSED( level );
      }
   }

   return nullptr;
}

Cell* PrimitiveStorage::getCell( const PrimitiveID& id )
{
   if ( cellExistsLocally( id ) )
   {
      for ( auto& [level, primitives] : cells_ )
      {
         if ( primitives.count( id ) > 0 )
         {
            return primitives[id].get();
         }
         WALBERLA_UNUSED( level );
      }
   }
   else if ( cellExistsInNeighborhood( id ) )
   {
      for ( auto& [level, primitives] : neighborCells_ )
      {
         if ( primitives.count( id ) > 0 )
         {
            return primitives[id].get();
         }
         WALBERLA_UNUSED( level );
      }
   }

   return nullptr;
}

std::vector< PrimitiveID > PrimitiveStorage::getPrimitiveIDs() const
{
   std::vector< PrimitiveID > ids;
   getPrimitiveIDs( ids );
   return ids;
}

std::vector< PrimitiveID > PrimitiveStorage::getVertexIDs() const
{
   std::vector< PrimitiveID > ids;
   getVertexIDs( ids );
   return ids;
}

std::vector< PrimitiveID > PrimitiveStorage::getEdgeIDs() const
{
   std::vector< PrimitiveID > ids;
   getEdgeIDs( ids );
   return ids;
}

std::vector< PrimitiveID > PrimitiveStorage::getFaceIDs() const
{
   std::vector< PrimitiveID > ids;
   getFaceIDs( ids );
   return ids;
}

std::vector< PrimitiveID > PrimitiveStorage::getCellIDs() const
{
   std::vector< PrimitiveID > ids;
   getCellIDs( ids );
   return ids;
}

std::vector< PrimitiveID > PrimitiveStorage::getVolumeIDs() const
{
   std::vector< PrimitiveID > ids;
   if ( !hasGlobalCells() )
   {
      getFaceIDs( ids );
   }
   else
   {
      getCellIDs( ids );
   }
   return ids;
}

void PrimitiveStorage::getPrimitiveIDs( std::vector< PrimitiveID >& primitiveIDs ) const
{
   primitiveIDs.clear();

   std::vector< PrimitiveID > someIDs;

   getVertexIDs( someIDs );
   primitiveIDs.insert( primitiveIDs.end(), someIDs.begin(), someIDs.end() );

   getEdgeIDs( someIDs );
   primitiveIDs.insert( primitiveIDs.end(), someIDs.begin(), someIDs.end() );

   getFaceIDs( someIDs );
   primitiveIDs.insert( primitiveIDs.end(), someIDs.begin(), someIDs.end() );

   getCellIDs( someIDs );
   primitiveIDs.insert( primitiveIDs.end(), someIDs.begin(), someIDs.end() );
}

void PrimitiveStorage::getNeighboringPrimitiveIDs( std::vector< PrimitiveID >& primitiveIDs ) const
{
   primitiveIDs.clear();

   std::vector< PrimitiveID > someIDs;

   getNeighboringVertexIDs( someIDs );
   primitiveIDs.insert( primitiveIDs.end(), someIDs.begin(), someIDs.end() );

   getNeighboringEdgeIDs( someIDs );
   primitiveIDs.insert( primitiveIDs.end(), someIDs.begin(), someIDs.end() );

   getNeighboringFaceIDs( someIDs );
   primitiveIDs.insert( primitiveIDs.end(), someIDs.begin(), someIDs.end() );

   getNeighboringCellIDs( someIDs );
   primitiveIDs.insert( primitiveIDs.end(), someIDs.begin(), someIDs.end() );
}

void PrimitiveStorage::getVertexIDs( std::vector< PrimitiveID >& vertexIDs ) const
{
   vertexIDs.clear();
   for ( const auto& [level, primitives] : vertices_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            vertexIDs.push_back( pid );
         }
      }
      WALBERLA_UNUSED( level );
   }
}

void PrimitiveStorage::getNeighboringVertexIDs( std::vector< PrimitiveID >& vertexIDs ) const
{
   vertexIDs.clear();
   for ( const auto& [level, primitives] : neighborVertices_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            vertexIDs.push_back( pid );
         }
      }
      WALBERLA_UNUSED( level );
   }
}

void PrimitiveStorage::getEdgeIDs( std::vector< PrimitiveID >& edgeIDs ) const
{
   edgeIDs.clear();
   for ( const auto& [level, primitives] : edges_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            edgeIDs.push_back( pid );
         }
      }
      WALBERLA_UNUSED( level );
   }
}

void PrimitiveStorage::getNeighboringEdgeIDs( std::vector< PrimitiveID >& edgeIDs ) const
{
   edgeIDs.clear();
   for ( const auto& [level, primitives] : neighborEdges_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            edgeIDs.push_back( pid );
         }
      }
      WALBERLA_UNUSED( level );
   }
}

void PrimitiveStorage::getFaceIDs( std::vector< PrimitiveID >& faceIDs ) const
{
   faceIDs.clear();
   for ( const auto& [level, primitives] : faces_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            faceIDs.push_back( pid );
         }
      }
      WALBERLA_UNUSED( level );
   }
}

void PrimitiveStorage::getNeighboringFaceIDs( std::vector< PrimitiveID >& faceIDs ) const
{
   faceIDs.clear();
   for ( const auto& [level, primitives] : neighborFaces_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            faceIDs.push_back( pid );
         }
      }
      WALBERLA_UNUSED( level );
   }
}

void PrimitiveStorage::getCellIDs( std::vector< PrimitiveID >& cellIDs ) const
{
   cellIDs.clear();
   for ( const auto& [level, primitives] : cells_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            cellIDs.push_back( pid );
         }
      }
      WALBERLA_UNUSED( level );
   }
}

void PrimitiveStorage::getNeighboringCellIDs( std::vector< PrimitiveID >& cellIDs ) const
{
   cellIDs.clear();
   for ( const auto& [level, primitives] : neighborCells_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            cellIDs.push_back( pid );
         }
      }
      WALBERLA_UNUSED( level );
   }
}

bool PrimitiveStorage::onBoundary( const PrimitiveID& primitiveID, const bool& highestDimensionAlwaysInner ) const
{
   WALBERLA_CHECK( primitiveExistsLocally( primitiveID ) || primitiveExistsInNeighborhood( primitiveID ),
                   "Cannot check if primitive (ID: " << primitiveID
                                                     << ") "
                                                        "is on boundary since it is not located on this process (rank: "
                                                     << walberla::mpi::MPIManager::instance()->rank() << ")." );

   if ( !hasGlobalCells() )
   {
      // 2D
      if ( highestDimensionAlwaysInner && faceExistsLocally( primitiveID ) )
      {
         return false;
      }
      if ( edgeExistsLocally( primitiveID ) || edgeExistsInNeighborhood( primitiveID ) )
      {
         const auto edge = getEdge( primitiveID );
         WALBERLA_ASSERT_GREATER( edge->getNumNeighborFaces(), 0 );
         WALBERLA_ASSERT_LESS_EQUAL( edge->getNumNeighborFaces(), 2 );
         return edge->getNumNeighborFaces() == 1;
      }
      else
      {
         const auto                 primitive = getPrimitive( primitiveID );
         std::vector< PrimitiveID > neighborEdges;
         primitive->getNeighborEdges( neighborEdges );
         for ( auto it : neighborEdges )
         {
            if ( onBoundary( it ) )
            {
               return true;
            }
         }
         return false;
      }
   }
   else
   {
      // 3D
      if ( highestDimensionAlwaysInner && cellExistsLocally( primitiveID ) )
      {
         return false;
      }
      if ( faceExistsLocally( primitiveID ) || faceExistsInNeighborhood( primitiveID ) )
      {
         const auto face = getFace( primitiveID );
         WALBERLA_ASSERT_GREATER( face->getNumNeighborCells(), 0 );
         WALBERLA_ASSERT_LESS_EQUAL( face->getNumNeighborCells(), 2 );
         return face->getNumNeighborCells() == 1;
      }
      else
      {
         const auto                 primitive = getPrimitive( primitiveID );
         std::vector< PrimitiveID > neighborFaces;
         primitive->getNeighborFaces( neighborFaces );
         for ( auto it : neighborFaces )
         {
            if ( onBoundary( it ) )
            {
               return true;
            }
         }
         return false;
      }
   }
}

uint_t PrimitiveStorage::getPrimitiveRank( const PrimitiveID& id ) const
{
   WALBERLA_ASSERT( primitiveExistsLocally( id ) || primitiveExistsInNeighborhood( id ) );
   if ( primitiveExistsLocally( id ) )
   {
      return uint_c( walberla::mpi::MPIManager::instance()->rank() );
   }
   else
   {
      return getNeighborPrimitiveRank( id );
   }
}

uint_t PrimitiveStorage::getNeighborPrimitiveRank( const PrimitiveID& id ) const
{
   WALBERLA_ASSERT( primitiveExistsInNeighborhood( id ), "Primitive with ID " << id << " does not exist in neighborhood." );
   for ( const auto& [level, nranks] : neighborRanks_ )
   {
      if ( nranks.count( id ) > 0 )
      {
         return nranks.at( id );
      }
      WALBERLA_UNUSED( level );
   }
   WALBERLA_CHECK( false, "Could not determine neighbor rank of neighbor primitive. This is bad." );
   return std::numeric_limits< uint_t >::max();
}

std::map< PrimitiveID, uint_t > PrimitiveStorage::getGlobalPrimitiveRanks() const
{
   std::map< PrimitiveID, uint_t > primitiveRanks;

   uint_t rank = uint_c( walberla::mpi::MPIManager::instance()->rank() );

   SendBuffer migrationInfoSendBuffer;
   RecvBuffer migrationInfoRecvBuffer;

   for ( const auto& primitiveID : getPrimitiveIDs() )
   {
      migrationInfoSendBuffer << primitiveID;
      migrationInfoSendBuffer << rank;
   }

   walberla::mpi::allGathervBuffer(
       migrationInfoSendBuffer, migrationInfoRecvBuffer, walberla::mpi::MPIManager::instance()->comm() );

   while ( !migrationInfoRecvBuffer.isEmpty() )
   {
      PrimitiveID primitiveID;
      uint_t      globalRank;

      migrationInfoRecvBuffer >> primitiveID;
      migrationInfoRecvBuffer >> globalRank;

      primitiveRanks[primitiveID] = globalRank;
   }

   return primitiveRanks;
}

void PrimitiveStorage::migratePrimitives( const MigrationInfo& migrationInfo )
{
   WALBERLA_CHECK_EQUAL( getCurrentLocalMaxRefinement(), 0, "Primitive migration not implemented for refined meshes." );

   WALBERLA_DEBUG_SECTION()
   {
      checkConsistency();
   }

   const auto& primitivesToMigrate = migrationInfo.getMap();

   WALBERLA_CHECK_EQUAL( primitivesToMigrate.size(), getNumberOfLocalPrimitives() );

   WALBERLA_DEBUG_SECTION()
   {
      for ( const auto& it : primitivesToMigrate )
      {
         WALBERLA_CHECK( primitiveExistsLocally( it.first ) );
      }
   }

   const uint_t rank         = uint_c( walberla::mpi::MPIManager::instance()->rank() );
   const uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   // --- Neighborhood rank update ---
   //
   // After the primitive migration is done, we need to update the ranks of all neighboring primitives.
   // To achieve this, each primitive must, after the migration, know the ranks of all its neighbors.
   // We collect this information in multiple steps to avoid expensive global communication:
   //
   //  1. build a map of all local primitives to their ranks after migration
   //  2. send this map to all neighboring processes
   //  3. build a map for each local primitive that maps the primitive's neighbor primitives
   //     to their future ranks
   //  4. migrate this map together with the primitives
   //     after migration, all processes should have for each now local primitive, a map
   //     with all of its neighbors
   //  5. clear the current neighboring ranks map and fill it up again
   //
   // The exact implementation differs a bit but this is the idea.

   // Neighborhood rank update: step 1
   MigrationMap_T localPrimitiveToFutureRankMap;
   for ( const auto& pID : getPrimitiveIDs() )
   {
      if ( primitivesToMigrate.count( pID ) > 0 )
      {
         localPrimitiveToFutureRankMap[pID] = primitivesToMigrate.at( pID );
      }
      else
      {
         localPrimitiveToFutureRankMap[pID] = rank;
      }
   }

   // Neighborhood rank update: step 2
   walberla::mpi::BufferSystem neighborhoodBufferSystem( walberla::mpi::MPIManager::instance()->comm() );
   std::set< uint_t >          neighborRanks;
   std::set< MPIRank >         neighborRanksInt;

   getNeighboringRanks( neighborRanks );
   for ( auto r : neighborRanks )
   {
      neighborRanksInt.insert( MPIRank( r ) );
   }
   neighborhoodBufferSystem.setReceiverInfo( std::begin( neighborRanksInt ), std::end( neighborRanksInt ), true );

   for ( auto r : neighborRanks )
   {
      neighborhoodBufferSystem.sendBuffer( r ) << localPrimitiveToFutureRankMap;
   }

   neighborhoodBufferSystem.sendAll();

   // Neighborhood rank update: step 3, but we have single map that contains all after-migration neighborhood information
   MigrationMap_T neighborhoodPrimitiveToFutureRankMap;
   neighborhoodPrimitiveToFutureRankMap.insert( std::begin( localPrimitiveToFutureRankMap ),
                                                std::end( localPrimitiveToFutureRankMap ) );

   for ( auto it = neighborhoodBufferSystem.begin(); it != neighborhoodBufferSystem.end(); ++it )
   {
      MigrationMap_T singleNeighborMigrationMap;
      it.buffer() >> singleNeighborMigrationMap;
      for ( const auto& itt : singleNeighborMigrationMap )
      {
         // let's do some sanity checks here
         WALBERLA_ASSERT_EQUAL( neighborhoodPrimitiveToFutureRankMap.count( itt.first ), 0 );
         neighborhoodPrimitiveToFutureRankMap[itt.first] = itt.second;
      }
   }

   // Let's start with the actual migration...

   PackageBufferSystem bufferSystem(
       walberla::mpi::MPIManager::instance()->comm(), migrationInfo.getNumReceivingPrimitives(), 8964 );

   ///////////////////////////////////
   // Serialization and sender side //
   ///////////////////////////////////

   WALBERLA_ASSERT_EQUAL( primitivesToMigrate.size(), getNumberOfLocalPrimitives() );

   for ( const auto& primitiveToMigrate : primitivesToMigrate )
   {
      const PrimitiveID primitiveID = primitiveToMigrate.first;
      const uint_t      targetRank  = primitiveToMigrate.second;

      WALBERLA_CHECK( primitiveExistsLocally( primitiveID ), "Cannot migrate non-locally-existent primitives." );
      WALBERLA_CHECK_LESS( targetRank, numProcesses );

      // Serialize primitives:
      // Create one serialization callback per primitive that shall be serialized.
      // for all primitives to be migrated:
      // - true (signals that this is no empty message)
      // - primitive type (Vertex, Edge, ...)
      // - primitive (contains neighborhood IDs and metadata like coordinates etc.)
      // - primitive data
      //   - data that was added to all primitives
      //   - data that was added to the type of primitive
      // - for all neighbors:
      //   - neighbor type
      //   - neighbor primitive
      //   - neighbor primitive rank after redistribution

      auto& sendBuffer = bufferSystem.getPackageSendBuffer( targetRank );

      if ( targetRank == rank )
      {
         // we do not want to send already local primitives
         sendBuffer << false;
         continue;
      }

      const PrimitiveTypeEnum primitiveType = getPrimitiveType( primitiveID );
      const Primitive*        primitive     = getPrimitive( primitiveID );

      WALBERLA_CHECK_NOT_IDENTICAL( primitiveType, PrimitiveTypeEnum::INVALID, "Sending invalid primitive type..." );

      sendBuffer << true;
      sendBuffer << primitiveType;
      sendBuffer << *primitive;

      serializeAllPrimitiveData( sendBuffer, primitiveID );

      // Neighborhood serialization
      // the number of neighbors of the sent primitive is already serialized
      // in its metadata - so we do not need to send it
      std::vector< PrimitiveID > neighborhood;
      primitive->getNeighborPrimitives( neighborhood );
      for ( const auto& neighborID : neighborhood )
      {
         WALBERLA_CHECK( primitiveExistsLocally( neighborID ) || primitiveExistsInNeighborhood( neighborID ) );

         const PrimitiveTypeEnum neighborType      = getPrimitiveType( neighborID );
         const Primitive*        neighborPrimitive = getPrimitive( neighborID );

         sendBuffer << neighborType;
         sendBuffer << *neighborPrimitive;

         // Neighborhood rank update: step 4
         WALBERLA_ASSERT_EQUAL( neighborhoodPrimitiveToFutureRankMap.count( neighborID ), 1 );
         sendBuffer << neighborhoodPrimitiveToFutureRankMap[neighborID];
      }
   }

   bufferSystem.sendAll();

   ///////////////////////////////////////
   // Deserialization and receiver side //
   ///////////////////////////////////////

   while ( !bufferSystem.allPackagesReceived() )
   {
      auto package    = bufferSystem.getNextPackage();
      auto recvBuffer = package.buffer();

      bool hasContent;
      recvBuffer >> hasContent;

      // already local primitives are not sent
      if ( hasContent )
      {
         const PrimitiveID primitiveID = deserializeAndAddPrimitive( recvBuffer, false );

         initializeAndDeserializeAllPrimitiveData( recvBuffer, primitiveID );

         // Neighborhood deserialization
         const Primitive* primitive    = getPrimitive( primitiveID );
         const uint_t     numNeighbors = primitive->getNumNeighborPrimitives();

         for ( uint_t neighborCnt = 0; neighborCnt < numNeighbors; neighborCnt++ )
         {
            PrimitiveID neighborPrimitiveID = deserializeAndAddPrimitive( recvBuffer, true );

            // Here we fill up the neighborhood migration map so that it should contain at least
            // all neighbor-primitive-ranks of staying primitives and
            // all neighbor-primitive-ranks of received primitives,
            // and possibly some more ...
            uint_t rankAfterMigration;
            recvBuffer >> rankAfterMigration;
            neighborhoodPrimitiveToFutureRankMap[neighborPrimitiveID] = rankAfterMigration;
         }
      }
   }

   // Neighborhood rank update: step 5a)
   neighborRanks_.clear();

   /////////////////////////////////////////////////////////////////////////////////////////////////////////
   // Erasing the migrated primitives from the locally allocated but creating copies for the neighborhood //
   /////////////////////////////////////////////////////////////////////////////////////////////////////////

   for ( const auto& it : primitivesToMigrate )
   {
      PrimitiveID idToErase  = it.first;
      uint_t      targetRank = it.second;
      WALBERLA_CHECK( primitiveExistsLocally( idToErase ), "Cannot erase non-locally-existent primitives." );
      if ( targetRank != rank ) // only erase local primitives that were migrated to other ranks than mine
      {
         if ( vertexExistsLocally( idToErase ) )
         {
            neighborVertices_[0][idToErase] = std::make_shared< Vertex >( *getVertex( idToErase ) );
            vertices_[0].erase( idToErase );
         }
         if ( edgeExistsLocally( idToErase ) )
         {
            neighborEdges_[0][idToErase] = std::make_shared< Edge >( *getEdge( idToErase ) );
            edges_[0].erase( idToErase );
         }
         if ( faceExistsLocally( idToErase ) )
         {
            neighborFaces_[0][idToErase] = std::make_shared< Face >( *getFace( idToErase ) );
            faces_[0].erase( idToErase );
         }
         if ( cellExistsLocally( idToErase ) )
         {
            neighborCells_[0][idToErase] = std::make_shared< Cell >( *getCell( idToErase ) );
            cells_[0].erase( idToErase );
         }
      }
   }

   /////////////////////////////////////////////////////////////////////////////
   // Erasing all neighborhood primitives that are also locally allocated now //
   /////////////////////////////////////////////////////////////////////////////

   std::vector< PrimitiveID > localPrimitiveIDs;
   getPrimitiveIDs( localPrimitiveIDs );
   for ( const auto& localID : localPrimitiveIDs )
   {
      if ( vertexExistsInNeighborhood( localID ) )
         neighborVertices_[0].erase( localID );
      if ( edgeExistsInNeighborhood( localID ) )
         neighborEdges_[0].erase( localID );
      if ( faceExistsInNeighborhood( localID ) )
         neighborFaces_[0].erase( localID );
      if ( cellExistsInNeighborhood( localID ) )
         neighborCells_[0].erase( localID );
   }

   /////////////////////////////////////////////////////////////////////////////////////////
   // Erasing all neighbors that are not referenced by local primitives from neighborhood //
   /////////////////////////////////////////////////////////////////////////////////////////

   // enjoy worst possible complexity...

   std::vector< PrimitiveID > neighborhoodIDs;
   for ( const auto& it : neighborVertices_.at( 0 ) )
      neighborhoodIDs.push_back( it.first );
   for ( const auto& it : neighborEdges_.at( 0 ) )
      neighborhoodIDs.push_back( it.first );
   for ( const auto& it : neighborFaces_.at( 0 ) )
      neighborhoodIDs.push_back( it.first );
   for ( const auto& it : neighborCells_.at( 0 ) )
      neighborhoodIDs.push_back( it.first );

   for ( const auto& neighborhoodID : neighborhoodIDs )
   {
      bool referenced = false;
      for ( const auto& localID : localPrimitiveIDs )
      {
         if ( referenced )
            break;

         Primitive* primitive = getPrimitive( localID );

         std::vector< PrimitiveID > neighborIDs;
         primitive->getNeighborPrimitives( neighborIDs );

         for ( const auto& neighborOfLocalID : neighborIDs )
         {
            if ( neighborOfLocalID == neighborhoodID )
            {
               referenced = true;
               break;
            }
         }
      }

      if ( !referenced )
      {
         if ( vertexExistsInNeighborhood( neighborhoodID ) )
            neighborVertices_[0].erase( neighborhoodID );
         if ( edgeExistsInNeighborhood( neighborhoodID ) )
            neighborEdges_[0].erase( neighborhoodID );
         if ( faceExistsInNeighborhood( neighborhoodID ) )
            neighborFaces_[0].erase( neighborhoodID );
         if ( cellExistsInNeighborhood( neighborhoodID ) )
            neighborCells_[0].erase( neighborhoodID );
      }
   }

   /////////////////////////////////
   // Updating neighborhood ranks //
   /////////////////////////////////

   neighborRanks_.clear();
   std::vector< PrimitiveID > neighborPrimitives;
   getNeighboringPrimitiveIDs( neighborPrimitives );

   for ( const auto& np : neighborPrimitives )
   {
      // Neighborhood rank update: step 5b)
      WALBERLA_CHECK_GREATER( neighborhoodPrimitiveToFutureRankMap.count( np ), 0 );
      neighborRanks_[0][np] = neighborhoodPrimitiveToFutureRankMap[np];
   }

   for ( uint_t i = 0; i < additionalHaloDepth_; i++ )
   {
      addDirectNeighborsDistributed();
   }

   splitCommunicatorByPrimitiveDistribution();

   wasModified();
   updateLeafPrimitiveMaps();

   WALBERLA_DEBUG_SECTION()
   {
      checkConsistency();
   }
}

std::set< uint_t > PrimitiveStorage::getNeighboringFaceRanksOfFace( const PrimitiveID& facePrimitiveID ) const
{
   WALBERLA_CHECK( additionalHaloDepth_ > 0,
                   "Additional halo depth must be larger than 0 to access neighbor faces of faces. "
                   "This can be set in the PrimitiveStorage constructor." );
   WALBERLA_CHECK( faceExistsLocally( facePrimitiveID ), "Face " << facePrimitiveID << "does not exist locally." );

   std::set< uint_t > neighboringRanks;
   const auto         face = getFace( facePrimitiveID );
   for ( const auto& npid : face->getIndirectNeighborFaceIDsOverVertices() )
   {
      WALBERLA_CHECK( faceExistsLocally( npid ) || faceExistsInNeighborhood( npid ),
                      "Neighbor face " << npid << " of " << facePrimitiveID << " does not exist locally, nor in neighborhood." );
      if ( faceExistsInNeighborhood( npid ) )
      {
         neighboringRanks.insert( getNeighborPrimitiveRank( npid ) );
      }
   }
   return neighboringRanks;
}

std::set< uint_t > PrimitiveStorage::getNeighboringFaceRanksOfAllFaces() const
{
   std::set< uint_t > neighborRanks;
   for ( const auto& it : faces_.at( 0 ) )
   {
      const auto nr = getNeighboringFaceRanksOfFace( it.first );
      neighborRanks.insert( nr.begin(), nr.end() );
   }
   return neighborRanks;
}

std::set< uint_t > PrimitiveStorage::getNeighboringCellRanksOfCell( const PrimitiveID& cellPrimitiveID ) const
{
   WALBERLA_CHECK_EQUAL( getCurrentLocalMaxRefinement(), 0, "Not implemented for refined meshes." );

   WALBERLA_CHECK( additionalHaloDepth_ > 0,
                   "Additional halo depth must be larger than 0 to access neighbor cells of cells."
                   "This can be set in the PrimitiveStorage constructor." );
   WALBERLA_CHECK( cellExistsLocally( cellPrimitiveID ), "Cell " << cellPrimitiveID << "does not exist locally." );

   std::set< uint_t > neighboringRanks;
   const auto         cell = getCell( cellPrimitiveID );
   for ( const auto& npid : cell->getIndirectNeighborCellIDsOverVertices() )
   {
      WALBERLA_CHECK( cellExistsLocally( npid ) || cellExistsInNeighborhood( npid ),
                      "Neighbor cell " << npid << " of " << cellPrimitiveID << " does not exist locally, nor in neighborhood." );
      if ( cellExistsInNeighborhood( npid ) )
      {
         neighboringRanks.insert( getNeighborPrimitiveRank( npid ) );
      }
   }
   return neighboringRanks;
}

std::set< uint_t > PrimitiveStorage::getNeighboringCellRanksOfAllCells() const
{
   std::set< uint_t > neighborRanks;
   for ( const auto& it : cells_.at( 0 ) )
   {
      const auto nr = getNeighboringCellRanksOfCell( it.first );
      neighborRanks.insert( nr.begin(), nr.end() );
   }
   return neighborRanks;
}

std::set< uint_t > PrimitiveStorage::getNeighboringVolumeRanksOfVolume( const PrimitiveID& volumePrimitiveID ) const
{
   if ( hasGlobalCells() )
   {
      return getNeighboringCellRanksOfCell( volumePrimitiveID );
   }
   else
   {
      return getNeighboringFaceRanksOfFace( volumePrimitiveID );
   }
}

std::set< uint_t > PrimitiveStorage::getNeighboringVolumeRanksOfAllVolumes() const
{
   if ( hasGlobalCells() )
   {
      return getNeighboringCellRanksOfAllCells();
   }
   else
   {
      return getNeighboringFaceRanksOfAllFaces();
   }
}

std::set< uint_t > PrimitiveStorage::getNeighboringRanks() const
{
   std::set< uint_t > neighboringRanks;
   getNeighboringRanks( neighboringRanks );
   return neighboringRanks;
}

void PrimitiveStorage::getNeighboringRanks( std::set< uint_t >& neighboringRanks ) const
{
   neighboringRanks.clear();
   for ( const auto& [level, nranks] : neighborRanks_ )
   {
      for ( const auto& [pid, r] : nranks )
      {
         auto primitive = getPrimitive( pid );
         if ( !primitive->hasChildren() )
         {
            neighboringRanks.insert( r );
         }
      }
      WALBERLA_UNUSED( level );
   }
}

void PrimitiveStorage::getNeighboringRanks( std::set< walberla::mpi::MPIRank >& neighboringRanks ) const
{
   neighboringRanks.clear();

   std::set< uint_t > neighboringRanksUint;
   getNeighboringRanks( neighboringRanksUint );

   for ( const auto& it : neighboringRanksUint )
   {
      neighboringRanks.insert( static_cast< walberla::mpi::MPIRank >( it ) );
   }
}

PrimitiveStorage::PrimitiveTypeEnum PrimitiveStorage::getPrimitiveType( const PrimitiveID& primitiveID ) const
{
   if ( primitiveExistsLocally( primitiveID ) || primitiveExistsInNeighborhood( primitiveID ) )
   {
      return getPrimitive( primitiveID )->getType();
   }
   return PrimitiveTypeEnum::INVALID;
}

PrimitiveID PrimitiveStorage::deserializeAndAddPrimitive( walberla::mpi::RecvBuffer& recvBuffer, const bool& isNeighborPrimitive )
{
   WALBERLA_CHECK_EQUAL( getCurrentLocalMaxRefinement(), 0, "Not implemented for refined meshes." );

   PrimitiveTypeEnum primitiveType;
   PrimitiveID       primitiveID;

   recvBuffer >> primitiveType;

   switch ( primitiveType )
   {
   case PrimitiveTypeEnum::VERTEX: {
      std::shared_ptr< Vertex > vertex = std::make_shared< Vertex >( recvBuffer );
      primitiveID                      = vertex->getID();
      if ( isNeighborPrimitive )
      {
         neighborVertices_[0][primitiveID] = vertex;
      }
      else
      {
         vertices_[0][primitiveID] = vertex;
      }
      break;
   }
   case PrimitiveTypeEnum::EDGE: {
      std::shared_ptr< Edge > edge = std::make_shared< Edge >( recvBuffer );
      primitiveID                  = edge->getID();
      if ( isNeighborPrimitive )
      {
         neighborEdges_[0][primitiveID] = edge;
      }
      else
      {
         edges_[0][primitiveID] = edge;
      }
      break;
   }
   case PrimitiveTypeEnum::FACE: {
      std::shared_ptr< Face > face = std::make_shared< Face >( recvBuffer );
      primitiveID                  = face->getID();
      if ( isNeighborPrimitive )
      {
         neighborFaces_[0][primitiveID] = face;
      }
      else
      {
         faces_[0][primitiveID] = face;
      }
      break;
   }
   case PrimitiveTypeEnum::CELL: {
      std::shared_ptr< Cell > cell = std::make_shared< Cell >( recvBuffer );
      primitiveID                  = cell->getID();
      if ( isNeighborPrimitive )
      {
         neighborCells_[0][primitiveID] = cell;
      }
      else
      {
         cells_[0][primitiveID] = cell;
      }
      break;
   }
   default:
      WALBERLA_ABORT( "Cannot deserialize primitive - unkown primitive type" );
      break;
   }

   return primitiveID;
}

void PrimitiveStorage::serializeAllPrimitiveData( walberla::mpi::SendBuffer& sendBuffer, const PrimitiveID& primitiveID )
{
   WALBERLA_CHECK_EQUAL( getCurrentLocalMaxRefinement(), 0, "Not implemented for refined meshes." );

   WALBERLA_ASSERT( primitiveExistsLocally( primitiveID ) );
   const PrimitiveTypeEnum primitiveType = getPrimitiveType( primitiveID );
   switch ( primitiveType )
   {
   case PrimitiveTypeEnum::VERTEX: {
      WALBERLA_ASSERT( vertexExistsLocally( primitiveID ) );
      auto vertex = vertices_[0][primitiveID];
      for ( const auto& serializationFunction : primitiveDataSerializationFunctions_ )
      {
         serializationFunction.second( vertex, sendBuffer );
      }
      for ( const auto& serializationFunction : vertexDataSerializationFunctions_ )
      {
         serializationFunction.second( vertex, sendBuffer );
      }
      break;
   }
   case PrimitiveTypeEnum::EDGE: {
      WALBERLA_ASSERT( edgeExistsLocally( primitiveID ) );
      auto edge = edges_[0][primitiveID];
      for ( const auto& serializationFunction : primitiveDataSerializationFunctions_ )
      {
         serializationFunction.second( edge, sendBuffer );
      }
      for ( const auto& serializationFunction : edgeDataSerializationFunctions_ )
      {
         serializationFunction.second( edge, sendBuffer );
      }
      break;
   }
   case PrimitiveTypeEnum::FACE: {
      WALBERLA_ASSERT( faceExistsLocally( primitiveID ) );
      auto face = faces_[0][primitiveID];
      for ( const auto& serializationFunction : primitiveDataSerializationFunctions_ )
      {
         serializationFunction.second( face, sendBuffer );
      }
      for ( const auto& serializationFunction : faceDataSerializationFunctions_ )
      {
         serializationFunction.second( face, sendBuffer );
      }
      break;
   }
   case PrimitiveTypeEnum::CELL: {
      WALBERLA_ASSERT( cellExistsLocally( primitiveID ) );
      auto cell = cells_[0][primitiveID];
      for ( const auto& serializationFunction : primitiveDataSerializationFunctions_ )
      {
         serializationFunction.second( cell, sendBuffer );
      }
      for ( const auto& serializationFunction : cellDataSerializationFunctions_ )
      {
         serializationFunction.second( cell, sendBuffer );
      }
      break;
   }
   default:
      WALBERLA_ABORT( "Invalid primitive type during serialization." );
      break;
   }
}

void PrimitiveStorage::initializeAndDeserializeAllPrimitiveData( walberla::mpi::RecvBuffer& recvBuffer,
                                                                 const PrimitiveID&         primitiveID )
{
   WALBERLA_CHECK_EQUAL( getCurrentLocalMaxRefinement(), 0, "Not implemented for refined meshes." );

   WALBERLA_ASSERT( primitiveExistsLocally( primitiveID ) );
   const PrimitiveTypeEnum primitiveType = getPrimitiveType( primitiveID );
   switch ( primitiveType )
   {
   case PrimitiveTypeEnum::VERTEX: {
      WALBERLA_ASSERT( vertexExistsLocally( primitiveID ) );
      auto vertex = vertices_[0][primitiveID];
      for ( const auto& initializationFunction : primitiveDataInitializationFunctions_ )
      {
         initializationFunction.second( vertex );
      }
      for ( const auto& initializationFunction : vertexDataInitializationFunctions_ )
      {
         initializationFunction.second( vertex );
      }
      for ( const auto& deserializationFunction : primitiveDataDeserializationFunctions_ )
      {
         deserializationFunction.second( vertex, recvBuffer );
      }
      for ( const auto& deserializationFunction : vertexDataDeserializationFunctions_ )
      {
         deserializationFunction.second( vertex, recvBuffer );
      }
      break;
   }
   case PrimitiveTypeEnum::EDGE: {
      WALBERLA_ASSERT( edgeExistsLocally( primitiveID ) );
      auto edge = edges_[0][primitiveID];
      for ( const auto& initializationFunction : primitiveDataInitializationFunctions_ )
      {
         initializationFunction.second( edge );
      }
      for ( const auto& initializationFunction : edgeDataInitializationFunctions_ )
      {
         initializationFunction.second( edge );
      }
      for ( const auto& deserializationFunction : primitiveDataDeserializationFunctions_ )
      {
         deserializationFunction.second( edge, recvBuffer );
      }
      for ( const auto& deserializationFunction : edgeDataDeserializationFunctions_ )
      {
         deserializationFunction.second( edge, recvBuffer );
      }
      break;
   }
   case PrimitiveTypeEnum::FACE: {
      WALBERLA_ASSERT( faceExistsLocally( primitiveID ) );
      auto face = faces_[0][primitiveID];
      for ( const auto& initializationFunction : primitiveDataInitializationFunctions_ )
      {
         initializationFunction.second( face );
      }
      for ( const auto& initializationFunction : faceDataInitializationFunctions_ )
      {
         initializationFunction.second( face );
      }
      for ( const auto& deserializationFunction : primitiveDataDeserializationFunctions_ )
      {
         deserializationFunction.second( face, recvBuffer );
      }
      for ( const auto& deserializationFunction : faceDataDeserializationFunctions_ )
      {
         deserializationFunction.second( face, recvBuffer );
      }
      break;
   }
   case PrimitiveTypeEnum::CELL: {
      WALBERLA_ASSERT( cellExistsLocally( primitiveID ) );
      auto cell = cells_[0][primitiveID];
      for ( const auto& initializationFunction : primitiveDataInitializationFunctions_ )
      {
         initializationFunction.second( cell );
      }
      for ( const auto& initializationFunction : cellDataInitializationFunctions_ )
      {
         initializationFunction.second( cell );
      }
      for ( const auto& deserializationFunction : primitiveDataDeserializationFunctions_ )
      {
         deserializationFunction.second( cell, recvBuffer );
      }
      for ( const auto& deserializationFunction : cellDataDeserializationFunctions_ )
      {
         deserializationFunction.second( cell, recvBuffer );
      }
      break;
   }
   default:
      WALBERLA_ABORT( "Invalid primitive type during initialization and deserialization." );
      break;
   }
}

std::string PrimitiveStorage::getGlobalInfo( bool onRootOnly ) const
{
   uint_t globalNumberOfVertices;
   uint_t globalNumberOfEdges;
   uint_t globalNumberOfFaces;
   uint_t globalNumberOfCells;
   uint_t globalNumberOfPrimitives;

   uint_t globalMaxNumberOfVertices;
   uint_t globalMaxNumberOfEdges;
   uint_t globalMaxNumberOfFaces;
   uint_t globalMaxNumberOfCells;
   uint_t globalMaxNumberOfPrimitives;

   uint_t globalMinNumberOfVertices;
   uint_t globalMinNumberOfEdges;
   uint_t globalMinNumberOfFaces;
   uint_t globalMinNumberOfCells;
   uint_t globalMinNumberOfPrimitives;

   if ( onRootOnly )
   {
      globalNumberOfVertices   = walberla::mpi::reduce( getNumberOfLocalVertices(), walberla::mpi::SUM );
      globalNumberOfEdges      = walberla::mpi::reduce( getNumberOfLocalEdges(), walberla::mpi::SUM );
      globalNumberOfFaces      = walberla::mpi::reduce( getNumberOfLocalFaces(), walberla::mpi::SUM );
      globalNumberOfCells      = walberla::mpi::reduce( getNumberOfLocalCells(), walberla::mpi::SUM );
      globalNumberOfPrimitives = walberla::mpi::reduce( getNumberOfLocalPrimitives(), walberla::mpi::SUM );

      globalMaxNumberOfVertices   = walberla::mpi::reduce( getNumberOfLocalVertices(), walberla::mpi::MAX );
      globalMaxNumberOfEdges      = walberla::mpi::reduce( getNumberOfLocalEdges(), walberla::mpi::MAX );
      globalMaxNumberOfFaces      = walberla::mpi::reduce( getNumberOfLocalFaces(), walberla::mpi::MAX );
      globalMaxNumberOfCells      = walberla::mpi::reduce( getNumberOfLocalCells(), walberla::mpi::MAX );
      globalMaxNumberOfPrimitives = walberla::mpi::reduce( getNumberOfLocalPrimitives(), walberla::mpi::MAX );

      globalMinNumberOfVertices   = walberla::mpi::reduce( getNumberOfLocalVertices(), walberla::mpi::MIN );
      globalMinNumberOfEdges      = walberla::mpi::reduce( getNumberOfLocalEdges(), walberla::mpi::MIN );
      globalMinNumberOfFaces      = walberla::mpi::reduce( getNumberOfLocalFaces(), walberla::mpi::MIN );
      globalMinNumberOfCells      = walberla::mpi::reduce( getNumberOfLocalCells(), walberla::mpi::MIN );
      globalMinNumberOfPrimitives = walberla::mpi::reduce( getNumberOfLocalPrimitives(), walberla::mpi::MIN );
   }
   else
   {
      globalNumberOfVertices   = walberla::mpi::allReduce( getNumberOfLocalVertices(), walberla::mpi::SUM );
      globalNumberOfEdges      = walberla::mpi::allReduce( getNumberOfLocalEdges(), walberla::mpi::SUM );
      globalNumberOfFaces      = walberla::mpi::allReduce( getNumberOfLocalFaces(), walberla::mpi::SUM );
      globalNumberOfCells      = walberla::mpi::allReduce( getNumberOfLocalCells(), walberla::mpi::SUM );
      globalNumberOfPrimitives = walberla::mpi::allReduce( getNumberOfLocalPrimitives(), walberla::mpi::SUM );

      globalMaxNumberOfVertices   = walberla::mpi::allReduce( getNumberOfLocalVertices(), walberla::mpi::MAX );
      globalMaxNumberOfEdges      = walberla::mpi::allReduce( getNumberOfLocalEdges(), walberla::mpi::MAX );
      globalMaxNumberOfFaces      = walberla::mpi::allReduce( getNumberOfLocalFaces(), walberla::mpi::MAX );
      globalMaxNumberOfCells      = walberla::mpi::allReduce( getNumberOfLocalCells(), walberla::mpi::MAX );
      globalMaxNumberOfPrimitives = walberla::mpi::allReduce( getNumberOfLocalPrimitives(), walberla::mpi::MAX );

      globalMinNumberOfVertices   = walberla::mpi::allReduce( getNumberOfLocalVertices(), walberla::mpi::MIN );
      globalMinNumberOfEdges      = walberla::mpi::allReduce( getNumberOfLocalEdges(), walberla::mpi::MIN );
      globalMinNumberOfFaces      = walberla::mpi::allReduce( getNumberOfLocalFaces(), walberla::mpi::MIN );
      globalMinNumberOfCells      = walberla::mpi::allReduce( getNumberOfLocalCells(), walberla::mpi::MIN );
      globalMinNumberOfPrimitives = walberla::mpi::allReduce( getNumberOfLocalPrimitives(), walberla::mpi::MIN );
   }

   const uint_t numberOfProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   const double globalAvgNumberOfVertices   = (double) globalNumberOfVertices / (double) numberOfProcesses;
   const double globalAvgNumberOfEdges      = (double) globalNumberOfEdges / (double) numberOfProcesses;
   const double globalAvgNumberOfFaces      = (double) globalNumberOfFaces / (double) numberOfProcesses;
   const double globalAvgNumberOfCells      = (double) globalNumberOfCells / (double) numberOfProcesses;
   const double globalAvgNumberOfPrimitives = (double) globalNumberOfPrimitives / (double) numberOfProcesses;

   walberla::math::DistributedSample neighborhoodSample;
   walberla::math::DistributedSample neighborhoodVolumeSample;
   const auto                        numNeighborProcesses = getNeighboringRanks().size();
   const auto numNeighborVolumeProcesses = additionalHaloDepth_ > 0 ? getNeighboringVolumeRanksOfAllVolumes().size() : 0;
   neighborhoodSample.castToRealAndInsert( numNeighborProcesses );
   neighborhoodVolumeSample.castToRealAndInsert( numNeighborVolumeProcesses );

   if ( onRootOnly )
   {
      neighborhoodSample.mpiGatherRoot();
      neighborhoodVolumeSample.mpiGatherRoot();
   }
   else
   {
      neighborhoodSample.mpiAllGather();
      neighborhoodVolumeSample.mpiAllGather();
   }

   std::stringstream os;
   os << "====================== PrimitiveStorage ======================\n";
   os << " - mesh dimensionality:                              " << ( hasGlobalCells() ? "3D" : "2D" ) << "\n";
   os << " - processes:                                        " << numberOfProcesses << "\n";
   os << " - neighbor processes (min, max, avg):               " << neighborhoodSample.min() << ", " << neighborhoodSample.max()
      << ", " << neighborhoodSample.avg() << "\n";
   if ( additionalHaloDepth_ > 0 )
   {
      os << " - neighbor processes, volumes only (min, max, avg): " << neighborhoodVolumeSample.min() << ", "
         << neighborhoodVolumeSample.max() << ", " << neighborhoodVolumeSample.avg() << "\n";
   }
   os << " - primitive distribution:\n";
   os << "                +-------------------------------------------+\n"
         "                |    total |      min |      max |      avg |\n"
         "   +------------+----------+----------+----------+----------+\n"
         "   | primitives | "
      << std::setw( 8 ) << globalNumberOfPrimitives << " | " << std::setw( 8 ) << globalMinNumberOfPrimitives << " | "
      << std::setw( 8 ) << globalMaxNumberOfPrimitives << " | " << std::fixed << std::setprecision( 1 ) << std::setw( 8 )
      << globalAvgNumberOfPrimitives
      << " |\n"
         "   +------------+----------+----------+----------+----------+\n"
         "   |   vertices | "
      << std::setw( 8 ) << globalNumberOfVertices << " | " << std::setw( 8 ) << globalMinNumberOfVertices << " | "
      << std::setw( 8 ) << globalMaxNumberOfVertices << " | " << std::fixed << std::setprecision( 1 ) << std::setw( 8 )
      << globalAvgNumberOfVertices
      << " |\n"
         "   +------------+----------+----------+----------+----------+\n"
         "   |      edges | "
      << std::setw( 8 ) << globalNumberOfEdges << " | " << std::setw( 8 ) << globalMinNumberOfEdges << " | " << std::setw( 8 )
      << globalMaxNumberOfEdges << " | " << std::fixed << std::setprecision( 1 ) << std::setw( 8 ) << globalAvgNumberOfEdges
      << " |\n"
         "   +------------+----------+----------+----------+----------+\n"
         "   |      faces | "
      << std::setw( 8 ) << globalNumberOfFaces << " | " << std::setw( 8 ) << globalMinNumberOfFaces << " | " << std::setw( 8 )
      << globalMaxNumberOfFaces << " | " << std::fixed << std::setprecision( 1 ) << std::setw( 8 ) << globalAvgNumberOfFaces
      << " |\n"
         "   +------------+----------+----------+----------+----------+\n"
         "   |      cells | "
      << std::setw( 8 ) << globalNumberOfCells << " | " << std::setw( 8 ) << globalMinNumberOfCells << " | " << std::setw( 8 )
      << globalMaxNumberOfCells << " | " << std::fixed << std::setprecision( 1 ) << std::setw( 8 ) << globalAvgNumberOfCells
      << " |\n"
         "   +------------+----------+----------+----------+----------+\n";
   os << "==============================================================\n";

   return os.str();
}

void PrimitiveStorage::checkConsistency()
{
   // 1. Number of data entries less than local counter
   // 2. PrimitiveIDs of maps match IDs of Primitives
   // 3. Neighborhood of Primitives
   for ( auto it = vertices_.at( 0 ).begin(); it != vertices_.at( 0 ).end(); it++ )
   {
      WALBERLA_CHECK_GREATER_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
      WALBERLA_CHECK_EQUAL( it->first, it->second->getID() );
      WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 0 );
   }
   for ( auto it = edges_.at( 0 ).begin(); it != edges_.at( 0 ).end(); it++ )
   {
      WALBERLA_CHECK_GREATER_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
      WALBERLA_CHECK_EQUAL( it->first, it->second->getID() );
      WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 2 );
   }
   for ( auto it = faces_.at( 0 ).begin(); it != faces_.at( 0 ).end(); it++ )
   {
      WALBERLA_CHECK_GREATER_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
      WALBERLA_CHECK_EQUAL( it->first, it->second->getID() );
      WALBERLA_CHECK_EQUAL( it->second->getNumNeighborVertices(), 3 );
      WALBERLA_CHECK_EQUAL( it->second->getNumNeighborEdges(), 3 );
   }
   for ( auto it = cells_.at( 0 ).begin(); it != cells_.at( 0 ).end(); it++ )
   {
      WALBERLA_CHECK_GREATER_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
      WALBERLA_CHECK_EQUAL( it->first, it->second->getID() );
      WALBERLA_CHECK_EQUAL( it->second->getNumNeighborVertices(), 4 );
      WALBERLA_CHECK_EQUAL( it->second->getNumNeighborEdges(), 6 );
      WALBERLA_CHECK_EQUAL( it->second->getNumNeighborFaces(), 4 );
   }

   // 4. Number of data entries of neighbor primitives is zero
   // 5. PrimitiveIDs of neighbor maps match IDs of neighbor Primitives
   // 6. Neighborhood of Primitives
   for ( auto it = neighborVertices_.at( 0 ).begin(); it != neighborVertices_.at( 0 ).end(); it++ )
   {
      WALBERLA_CHECK_GREATER_EQUAL( 0, it->second->getNumberOfDataEntries() );
      WALBERLA_CHECK_EQUAL( it->first, it->second->getID() );
      WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 0 );
   }
   for ( auto it = neighborEdges_.at( 0 ).begin(); it != neighborEdges_.at( 0 ).end(); it++ )
   {
      WALBERLA_CHECK_GREATER_EQUAL( 0, it->second->getNumberOfDataEntries() );
      WALBERLA_CHECK_EQUAL( it->first, it->second->getID() );
      WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 2 );
   }
   for ( auto it = neighborFaces_.at( 0 ).begin(); it != neighborFaces_.at( 0 ).end(); it++ )
   {
      WALBERLA_CHECK_GREATER_EQUAL( 0, it->second->getNumberOfDataEntries() );
      WALBERLA_CHECK_EQUAL( it->first, it->second->getID() );
      WALBERLA_CHECK_EQUAL( it->second->getNumNeighborVertices(), 3 );
      WALBERLA_CHECK_EQUAL( it->second->getNumNeighborEdges(), 3 );
   }
   for ( auto it = neighborCells_.at( 0 ).begin(); it != neighborCells_.at( 0 ).end(); it++ )
   {
      WALBERLA_CHECK_GREATER_EQUAL( 0, it->second->getNumberOfDataEntries() );
      WALBERLA_CHECK_EQUAL( it->first, it->second->getID() );
      WALBERLA_CHECK_EQUAL( it->second->getNumNeighborVertices(), 4 );
      WALBERLA_CHECK_EQUAL( it->second->getNumNeighborEdges(), 6 );
      WALBERLA_CHECK_EQUAL( it->second->getNumNeighborFaces(), 4 );
   }

   // 7. Number of callbacks is less or equal to the data handling counter
   WALBERLA_CHECK_LESS_EQUAL( primitiveDataInitializationFunctions_.size(), primitiveDataHandlers_ );
   WALBERLA_CHECK_LESS_EQUAL( primitiveDataSerializationFunctions_.size(), primitiveDataHandlers_ );
   WALBERLA_CHECK_LESS_EQUAL( primitiveDataDeserializationFunctions_.size(), primitiveDataHandlers_ );
   WALBERLA_CHECK_LESS_EQUAL( vertexDataInitializationFunctions_.size(), primitiveDataHandlers_ );
   WALBERLA_CHECK_LESS_EQUAL( vertexDataSerializationFunctions_.size(), primitiveDataHandlers_ );
   WALBERLA_CHECK_LESS_EQUAL( vertexDataDeserializationFunctions_.size(), primitiveDataHandlers_ );
   WALBERLA_CHECK_LESS_EQUAL( edgeDataInitializationFunctions_.size(), primitiveDataHandlers_ );
   WALBERLA_CHECK_LESS_EQUAL( edgeDataSerializationFunctions_.size(), primitiveDataHandlers_ );
   WALBERLA_CHECK_LESS_EQUAL( edgeDataDeserializationFunctions_.size(), primitiveDataHandlers_ );
   WALBERLA_CHECK_LESS_EQUAL( faceDataInitializationFunctions_.size(), primitiveDataHandlers_ );
   WALBERLA_CHECK_LESS_EQUAL( faceDataSerializationFunctions_.size(), primitiveDataHandlers_ );
   WALBERLA_CHECK_LESS_EQUAL( faceDataDeserializationFunctions_.size(), primitiveDataHandlers_ );

   std::vector< PrimitiveID > primitiveIDs;
   getPrimitiveIDs( primitiveIDs );

   // 8. All neighborIDs that are referenced by locally allocated primitives are allocated in storage
   for ( const auto& id : primitiveIDs )
   {
      std::vector< PrimitiveID > neighborhoodIDs;
      getPrimitive( id )->getNeighborPrimitives( neighborhoodIDs );
      for ( const auto& nID : neighborhoodIDs )
      {
         WALBERLA_CHECK( primitiveExistsLocally( nID ) || primitiveExistsInNeighborhood( nID ),
                         "Neighbor referenced from locally allocated primitive is not allocated in storage. ID: " << nID );
      }
   }

   // 9. As many neighbor ranks as neighbors
   WALBERLA_CHECK_EQUAL( neighborRanks_[0].size(),
                         neighborVertices_[0].size() + neighborEdges_[0].size() + neighborFaces_[0].size() +
                             neighborCells_[0].size() );

   // 10. Local primitives do not exist in neighborhood
   for ( const auto& id : primitiveIDs )
   {
      WALBERLA_CHECK( !primitiveExistsInNeighborhood( id ), "Primitive that exists in neighborhood: " << id );
   }

   // 11. If the mesh is 2D, we require that there are no cells
   if ( !hasGlobalCells() )
   {
      WALBERLA_CHECK_EQUAL( getNumberOfLocalCells(), 0, "Inconsistency regarding the number of cell primitives." );
   }

   // 12. If the mesh is 3D, we require that each vertex, edge and face has at least one neighboring cell.
   if ( hasGlobalCells() )
   {
      for ( auto it : getVertices() )
      {
         WALBERLA_CHECK_GREATER( it.second->getNumNeighborCells(), 0 );
      }
      for ( auto it : getEdges() )
      {
         WALBERLA_CHECK_GREATER( it.second->getNumNeighborCells(), 0 );
      }
      for ( auto it : getFaces() )
      {
         WALBERLA_CHECK_GREATER( it.second->getNumNeighborCells(), 0 );
      }
   }
}

void PrimitiveStorage::splitCommunicatorByPrimitiveDistribution()
{
   MPI_Comm originalComm = walberla::mpi::MPIManager::instance()->comm();

   WALBERLA_NON_MPI_SECTION()
   {
      splitComm_ = originalComm;
      return;
   }

   if ( getNumberOfEmptyProcesses() == 0 )
   {
      splitComm_ = originalComm;
   }
   const int includeInSubComm = getNumberOfLocalPrimitives() > 0 ? 1 : 0;
   int       oldRank;
   MPI_Comm_rank( originalComm, &oldRank );
   MPI_Comm subComm;
   MPI_Comm_split( originalComm, includeInSubComm, oldRank, &subComm );
   splitComm_ = subComm;
}

void PrimitiveStorage::queryRefinementAndCoarseningHanging( const std::vector< PrimitiveID >& volumePIDsRefine,
                                                            const std::vector< PrimitiveID >& volumePIDsCoarsen,
                                                            std::vector< PrimitiveID >&       volumePIDsRefineResult,
                                                            std::vector< PrimitiveID >&       volumePIDsCoarsenResult ) const
{
   volumePIDsRefineResult.clear();
   volumePIDsCoarsenResult.clear();

   //////////////////
   /// Refinement ///
   //////////////////

   WALBERLA_CHECK_EQUAL(
       walberla::mpi::MPIManager::instance()->numProcesses(), 1, "2:1 constraint not implemented for the parallel case." );

   // Building a map with all target levels for all local primitives.
   std::map< PrimitiveID, uint_t > targetLevels;

   for ( const auto& pid : getVolumeIDs() )
   {
      targetLevels[pid] = getRefinementLevel( pid );
   }

   for ( const auto& pid : volumePIDsRefine )
   {
      WALBERLA_CHECK_GREATER( targetLevels.count( pid ), 0 );
      targetLevels[pid]++;
   }

   bool constraintFulfilled = false;

   while ( !constraintFulfilled )
   {
      constraintFulfilled = true;

      if ( !hasGlobalCells() )
      {
         for ( const auto& [faceID, targetLevel] : targetLevels )
         {
            auto face = getFace( faceID );

            for ( const auto& neighborFaceID : face->getIndirectTopLevelNeighborFaceIDsOverVertices() )
            {
               WALBERLA_CHECK( faceID != neighborFaceID );
               auto levelNeighbor = targetLevels[neighborFaceID];

               if ( levelNeighbor + 1 < targetLevel )
               {
                  targetLevels[neighborFaceID]++;
                  constraintFulfilled = false;
               }
            }
         }
      }
      else
      {
         WALBERLA_ABORT( "Not implemented in 3D." );
      }
   }

   for ( const auto& [pid, targetLevel] : targetLevels )
   {
      if ( targetLevel > getRefinementLevel( pid ) )
      {
         WALBERLA_CHECK_EQUAL( targetLevel, getRefinementLevel( pid ) + 1 );
         volumePIDsRefineResult.push_back( pid );
      }
   }

   // TODO: implement coarsening

   WALBERLA_UNUSED( volumePIDsCoarsen );
}

void PrimitiveStorage::refinementAndCoarseningHanging( const std::vector< PrimitiveID >& volumePIDsRefine,
                                                       const std::vector< PrimitiveID >& volumePIDsCoarsen,
                                                       std::vector< PrimitiveID >&       volumePIDsRefineResult,
                                                       std::vector< PrimitiveID >&       volumePIDsCoarsenResult )
{
   WALBERLA_CHECK_GREATER_EQUAL( getAdditionalHaloDepth(),
                                 1,
                                 "For hanging node refinement and coarsening, the number of additional halos must be at least 1. "
                                 "This can be set in the PrimitiveStorage constructor." );

   //////////////////////////////////////////////////////////
   /// Checking what is going to be refined and coarsened ///
   //////////////////////////////////////////////////////////

   queryRefinementAndCoarseningHanging( volumePIDsRefine, volumePIDsCoarsen, volumePIDsRefineResult, volumePIDsCoarsenResult );

   //////////////////
   /// Refinement ///
   //////////////////

   // Since we require all passed faces to be locally allocated we need to communicate that they are flagged for refinement to the
   // neighboring processes. Then during refinement all marked primitives are refined on each process, also the neighboring
   // primitives that are marked.

   // Will contain all PIDs that shall be refined later on.
   // Includes locally allocated primitives as well as those that exist in the direct neighborhood.
   // To be filled below ...
   std::set< PrimitiveID > volumePIDsRefineIncludingNeighborhood;
   volumePIDsRefineIncludingNeighborhood.insert( volumePIDsRefineResult.begin(), volumePIDsRefineResult.end() );

   // Map of neighboring processes to vectors of locally allocated primitives that shall be refined to be sent to these processes.
   std::map< uint_t, std::set< PrimitiveID > > sendData;

   if ( !hasGlobalCells() )
   {
      walberla::mpi::BufferSystem bs( walberla::mpi::MPIManager::instance()->comm() );

      // Need to convert ranks to MPIRank datatype explicitly.
      const auto                         nranksUint = getNeighboringRanks();
      std::set< walberla::mpi::MPIRank > nranks( nranksUint.begin(), nranksUint.end() );

      bs.setReceiverInfo( nranks, true );

      for ( auto r : getNeighboringRanks() )
      {
         bs.sendBuffer( r ) << volumePIDsRefineResult;
      }

      bs.sendAll();

      for ( auto it = bs.begin(); it != bs.end(); ++it )
      {
         while ( !it.buffer().isEmpty() )
         {
            PrimitiveID pid;
            it.buffer() >> pid;
            volumePIDsRefineIncludingNeighborhood.insert( pid );
         }
      }
   }
   else
   {
      WALBERLA_ABORT( "Not implemented." );
   }

   // Loop over primitives that are going to be refined now.
   for ( const auto& coarseVolumeID : volumePIDsRefineIncludingNeighborhood )
   {
      if ( !hasGlobalCells() )
      {
         // 2D
         WALBERLA_CHECK(
             faceExistsLocally( coarseVolumeID ) || faceExistsInNeighborhood( coarseVolumeID ),
             "Face to be refined must exist locally or in neighborhood now that this data should have been communicated." );

         auto coarseFace = getFace( coarseVolumeID );
         auto level      = getRefinementLevel( coarseVolumeID );

         // Creating new primitives.
         //
         // First we need to check if the fine grid primitives already exist.
         // Note that this can only happen for vertices and edges that have parents that are vertices or edges themselves.
         //
         // Also, we are refining both, local and neighborhood primitives.

         // Vertices that are children of coarse vertices
         for ( uint_t i = 0; i < 3; i++ )
         {
            auto coarseVertexID = coarseFace->neighborVertices().at( i );
            auto coarseVertex   = getVertex( coarseVertexID );

            PrimitiveID fineVertexID;

            if ( coarseVertex->hasChildren() )
            {
               // Does exist.
               fineVertexID = coarseVertex->childVertices()[0];
            }
            else
            {
               // Does not exist.
               fineVertexID    = coarseVertexID.createChildren()[0];
               auto fineVertex = std::make_shared< Vertex >( fineVertexID, coarseVertex->getCoordinates() );

               coarseVertex->addChildVertices( { fineVertexID } );
               fineVertex->setParent( coarseVertexID );

               if ( vertexExistsLocally( coarseVertexID ) )
               {
                  vertices_[level + 1][fineVertexID] = fineVertex;
               }
               else
               {
                  neighborVertices_[level + 1][fineVertexID] = fineVertex;
               }

               fineVertex->meshBoundaryFlag_ = coarseVertex->getMeshBoundaryFlag();
            }
         }

         // Vertices that are children of coarse edges
         for ( uint_t i = 0; i < 3; i++ )
         {
            auto coarseEdgeID = coarseFace->neighborEdges().at( i );
            auto coarseEdge   = getEdge( coarseEdgeID );

            PrimitiveID fineVertexID;

            if ( coarseEdge->childVertices().size() > 0 )
            {
               // Does exist.
               fineVertexID = coarseEdge->childVertices()[0];
            }
            else
            {
               // Does not exist.
               fineVertexID    = coarseEdgeID.createChildren()[0];
               auto fineVertex = std::make_shared< Vertex >(
                   fineVertexID, 0.5 * ( coarseEdge->getCoordinates()[0] + coarseEdge->getCoordinates()[1] ) );

               coarseEdge->addChildVertices( { fineVertexID } );
               fineVertex->setParent( coarseEdgeID );

               if ( edgeExistsLocally( coarseEdgeID ) )
               {
                  vertices_[level + 1][fineVertexID] = fineVertex;
               }
               else
               {
                  neighborVertices_[level + 1][fineVertexID] = fineVertex;
               }

               fineVertex->meshBoundaryFlag_ = coarseEdge->getMeshBoundaryFlag();
            }
         }

         // For easy retrieval later, this map shall store the fine edges given the neighboring fine vertices as keys.
         // The vertex ID pair is sorted to allow for unique assignment.
         std::map< std::vector< PrimitiveID >, PrimitiveID > fineEdgeCache;

         // Edges that are children of coarse edges.
         // Locally sorting from bottom to top, then left to right in reference face.
         for ( uint_t i = 0; i < 3; i++ )
         {
            auto coarseEdgeID = coarseFace->neighborEdges().at( i );
            auto coarseEdge   = edges_[level][coarseEdgeID];

            // There are two child edges per edge.
            for ( uint_t ii = 0; ii < 2; ii++ )
            {
               PrimitiveID fineEdgeID;

               if ( coarseEdge->childEdges().size() == 2 )
               {
                  // Does exist.
                  fineEdgeID = coarseEdge->childEdges()[ii];
               }
               else
               {
                  // Does not exist.
                  fineEdgeID = coarseEdgeID.createChildren()[1 + ii]; // adding 1 because a child vertex was created already

                  auto v0 = getVertex( vertices_[level][coarseEdge->neighborVertices().at( ii )]->childVertices()[0] );
                  auto v1 = getVertex( coarseEdge->childVertices()[0] );

                  auto fineEdge =
                      std::make_shared< Edge >( fineEdgeID,
                                                v0->getID(),
                                                v1->getID(),
                                                std::array< Point3D, 2 >( { v0->getCoordinates(), v1->getCoordinates() } ) );

                  coarseEdge->addChildEdges( { fineEdgeID } );
                  fineEdge->setParent( coarseEdgeID );

                  if ( edgeExistsLocally( coarseEdgeID ) )
                  {
                     edges_[level + 1][fineEdgeID] = fineEdge;
                  }
                  else
                  {
                     neighborEdges_[level + 1][fineEdgeID] = fineEdge;
                  }

                  for ( const auto& neighborVertexID : fineEdge->neighborVertices() )
                  {
                     auto vertex = getVertex( neighborVertexID );
                     vertex->addEdge( fineEdge->getID() );
                  }

                  fineEdge->meshBoundaryFlag_ = coarseEdge->getMeshBoundaryFlag();
               }

               // Cache fine edge via neighboring vertex IDs.
               std::vector< PrimitiveID > vertexPair(
                   { getEdge( fineEdgeID )->neighborVertices()[0], getEdge( fineEdgeID )->neighborVertices()[1] } );
               std::sort( vertexPair.begin(), vertexPair.end() );
               fineEdgeCache[vertexPair] = fineEdgeID;
            }
         }

         // All remaining primitives cannot exist yet.
         // Let's begin creating the remaining 3 edges.

         auto edge_x_idx  = coarseVolumeID.createChildren()[0];
         auto edge_y_idx  = coarseVolumeID.createChildren()[1];
         auto edge_xy_idx = coarseVolumeID.createChildren()[2];

         auto coarseEdge0 = edges_[level][coarseFace->neighborEdges()[0]];
         auto coarseEdge1 = edges_[level][coarseFace->neighborEdges()[1]];
         auto coarseEdge2 = edges_[level][coarseFace->neighborEdges()[2]];

         auto fineVertexEdge0 = getVertex( coarseEdge0->childVertices()[0] );
         auto fineVertexEdge1 = getVertex( coarseEdge1->childVertices()[0] );
         auto fineVertexEdge2 = getVertex( coarseEdge2->childVertices()[0] );

         auto edge_x = std::make_shared< Edge >(
             edge_x_idx,
             fineVertexEdge1->getID(),
             fineVertexEdge2->getID(),
             std::array< Point3D, 2 >( { fineVertexEdge1->getCoordinates(), fineVertexEdge2->getCoordinates() } ) );

         auto edge_y = std::make_shared< Edge >(
             edge_y_idx,
             fineVertexEdge0->getID(),
             fineVertexEdge2->getID(),
             std::array< Point3D, 2 >( { fineVertexEdge0->getCoordinates(), fineVertexEdge2->getCoordinates() } ) );

         auto edge_xy = std::make_shared< Edge >(
             edge_xy_idx,
             fineVertexEdge0->getID(),
             fineVertexEdge1->getID(),
             std::array< Point3D, 2 >( { fineVertexEdge0->getCoordinates(), fineVertexEdge1->getCoordinates() } ) );

         if ( faceExistsLocally( coarseVolumeID ) )
         {
            edges_[level + 1][edge_x_idx]  = edge_x;
            edges_[level + 1][edge_y_idx]  = edge_y;
            edges_[level + 1][edge_xy_idx] = edge_xy;
         }
         else
         {
            neighborEdges_[level + 1][edge_x_idx]  = edge_x;
            neighborEdges_[level + 1][edge_y_idx]  = edge_y;
            neighborEdges_[level + 1][edge_xy_idx] = edge_xy;
         }

         for ( const auto& neighborVertexID : edge_x->neighborVertices() )
         {
            auto vertex = getVertex( neighborVertexID );
            vertex->addEdge( edge_x->getID() );
         }

         for ( const auto& neighborVertexID : edge_y->neighborVertices() )
         {
            auto vertex = getVertex( neighborVertexID );
            vertex->addEdge( edge_y->getID() );
         }

         for ( const auto& neighborVertexID : edge_xy->neighborVertices() )
         {
            auto vertex = getVertex( neighborVertexID );
            vertex->addEdge( edge_xy->getID() );
         }

         edge_x->meshBoundaryFlag_  = coarseFace->getMeshBoundaryFlag();
         edge_y->meshBoundaryFlag_  = coarseFace->getMeshBoundaryFlag();
         edge_xy->meshBoundaryFlag_ = coarseFace->getMeshBoundaryFlag();

         // Cache fine edge via neighboring vertex IDs.
         std::vector< PrimitiveID > vertexPair_x(
             { getEdge( edge_x_idx )->neighborVertices()[0], getEdge( edge_x_idx )->neighborVertices()[1] } );
         std::sort( vertexPair_x.begin(), vertexPair_x.end() );
         fineEdgeCache[vertexPair_x] = edge_x_idx;

         std::vector< PrimitiveID > vertexPair_y(
             { getEdge( edge_y_idx )->neighborVertices()[0], getEdge( edge_y_idx )->neighborVertices()[1] } );
         std::sort( vertexPair_y.begin(), vertexPair_y.end() );
         fineEdgeCache[vertexPair_y] = edge_y_idx;

         std::vector< PrimitiveID > vertexPair_xy(
             { getEdge( edge_xy_idx )->neighborVertices()[0], getEdge( edge_xy_idx )->neighborVertices()[1] } );
         std::sort( vertexPair_xy.begin(), vertexPair_xy.end() );
         fineEdgeCache[vertexPair_xy] = edge_xy_idx;

         edge_x->setParent( coarseVolumeID );
         edge_y->setParent( coarseVolumeID );
         edge_xy->setParent( coarseVolumeID );

         coarseFace->addChildEdges( { edge_x_idx, edge_y_idx, edge_xy_idx } );

         // Now the child faces.

         auto face_gray_0_idx = coarseVolumeID.createChildren()[3];
         auto face_gray_1_idx = coarseVolumeID.createChildren()[4];
         auto face_gray_2_idx = coarseVolumeID.createChildren()[5];
         auto face_blue_idx   = coarseVolumeID.createChildren()[6];

         auto coarseVertex0 = vertices_[level][coarseFace->neighborVertices()[0]];
         auto coarseVertex1 = vertices_[level][coarseFace->neighborVertices()[1]];
         auto coarseVertex2 = vertices_[level][coarseFace->neighborVertices()[2]];

         auto fineVertexVertex0 = getVertex( coarseVertex0->childVertices()[0] );
         auto fineVertexVertex1 = getVertex( coarseVertex1->childVertices()[0] );
         auto fineVertexVertex2 = getVertex( coarseVertex2->childVertices()[0] );

         // You better do not change this order - we need to point to specific child faces later!
         std::array< PrimitiveID, 4 > face_face_ids( { face_gray_0_idx, face_gray_1_idx, face_gray_2_idx, face_blue_idx } );

         std::array< std::array< PrimitiveID, 3 >, 4 > fine_face_vertices;

         // You better do not change this order! This is as defined by Bey to ensure proper congruence after successive refinement.
         fine_face_vertices[0] = { { fineVertexVertex0->getID(), fineVertexEdge0->getID(), fineVertexEdge1->getID() } };
         fine_face_vertices[1] = { { fineVertexEdge0->getID(), fineVertexVertex1->getID(), fineVertexEdge2->getID() } };
         fine_face_vertices[2] = { { fineVertexEdge1->getID(), fineVertexEdge2->getID(), fineVertexVertex2->getID() } };
         fine_face_vertices[3] = { { fineVertexEdge0->getID(), fineVertexEdge1->getID(), fineVertexEdge2->getID() } };

         for ( uint_t i = 0; i < 4; i++ )
         {
            auto v0 = fine_face_vertices[i][0];
            auto v1 = fine_face_vertices[i][1];
            auto v2 = fine_face_vertices[i][2];

            std::vector< PrimitiveID > edge_0_vids( { v0, v1 } );
            std::vector< PrimitiveID > edge_1_vids( { v0, v2 } );
            std::vector< PrimitiveID > edge_2_vids( { v1, v2 } );

            std::sort( edge_0_vids.begin(), edge_0_vids.end() );
            std::sort( edge_1_vids.begin(), edge_1_vids.end() );
            std::sort( edge_2_vids.begin(), edge_2_vids.end() );

            auto e0 = fineEdgeCache[edge_0_vids];
            auto e1 = fineEdgeCache[edge_1_vids];
            auto e2 = fineEdgeCache[edge_2_vids];

            auto e0v0 = getEdge( e0 )->neighborVertices()[0];
            auto e0v1 = getEdge( e0 )->neighborVertices()[1];

            auto e1v0 = getEdge( e1 )->neighborVertices()[0];
            auto e1v1 = getEdge( e1 )->neighborVertices()[1];

            auto e2v0 = getEdge( e2 )->neighborVertices()[0];
            auto e2v1 = getEdge( e2 )->neighborVertices()[1];

            std::array< int, 3 > edgeOrientation{};

            if ( e0v0 == v0 )
            {
               WALBERLA_CHECK( e0v1 == v1 );
               edgeOrientation[0] = 1;
            }
            else
            {
               WALBERLA_CHECK( e0v0 == v1 );
               WALBERLA_CHECK( e0v1 == v0 );
               edgeOrientation[0] = -1;
            }

            if ( e1v0 == v0 )
            {
               WALBERLA_CHECK( e1v1 == v2 );
               edgeOrientation[1] = 1;
            }
            else
            {
               WALBERLA_CHECK( e1v0 == v2 );
               WALBERLA_CHECK( e1v1 == v0 );
               edgeOrientation[1] = -1;
            }

            if ( e2v0 == v1 )
            {
               WALBERLA_CHECK( e2v1 == v2 );
               edgeOrientation[2] = 1;
            }
            else
            {
               WALBERLA_CHECK( e2v0 == v2 );
               WALBERLA_CHECK( e2v1 == v1 );
               edgeOrientation[2] = -1;
            }

            auto c0 = getVertex( v0 )->getCoordinates();
            auto c1 = getVertex( v1 )->getCoordinates();
            auto c2 = getVertex( v2 )->getCoordinates();

            auto fineFace = std::make_shared< Face >( face_face_ids[i],
                                                      std::array< PrimitiveID, 3 >( { v0, v1, v2 } ),
                                                      std::array< PrimitiveID, 3 >( { e0, e1, e2 } ),
                                                      edgeOrientation,
                                                      std::array< Point3D, 3 >( { c0, c1, c2 } ) );

            fineFace->setParent( coarseFace->getID() );
            coarseFace->addChildFaces( { fineFace->getID() } );

            if ( faceExistsLocally( coarseVolumeID ) )
            {
               faces_[level + 1][face_face_ids[i]] = fineFace;
            }
            else
            {
               neighborFaces_[level + 1][face_face_ids[i]] = fineFace;
            }

            // Adding this face as a neighbor to lower dim primitives.
            for ( const auto& neighborVertexID : fineFace->neighborVertices() )
            {
               auto vertex = getVertex( neighborVertexID );
               vertex->addFace( fineFace->getID() );
            }

            for ( const auto& neighborEdgeID : fineFace->neighborEdges() )
            {
               auto edge = getEdge( neighborEdgeID );
               edge->addFace( fineFace->getID() );
            }

            fineFace->meshBoundaryFlag_ = coarseFace->getMeshBoundaryFlag();
         }
      }
      else
      {
         // 3D
         WALBERLA_ABORT( "Refinement not implemented in 3D." );

         // TODO: the micro-volumes on the refined macros MUST eventually match the micros on a once refined macro!!!!!!!
         // TODO: also - consult the AMR ordering documentation!
      }

      // TODO: update neighborhood and neighbor processes!!! go via coarse level neighborhood since refinement is process local
      //       no primitives must be created that already exist in neighborhood!!!

      // TODO: cleanup: remove empty refinement levels from maps
   }

   // Clearing and rebuilding indirect volume neighborhood.
   // We only clear the indirect top-level neighborhood since that one changes during refinement.
   for ( const auto& [level, faceMap] : faces_ )
   {
      for ( auto& [faceID, face] : faceMap )
      {
         face->indirectTopLevelNeighborFaceIDsOverVertices_.clear();
         face->indirectTopLevelNeighborFaceIDsOverEdges_.clear();
         WALBERLA_UNUSED( faceID );
      }
      WALBERLA_UNUSED( level );
   }

   // Equal-level neighborhood.
   // This need to be set for all new volume primitives.
   for ( const auto& coarseVolumeID : volumePIDsRefineIncludingNeighborhood )
   {
      auto coarseFace = getFace( coarseVolumeID );
      for ( const auto& fineFaceID : coarseFace->childFaces() )
      {
         auto fineFace = getFace( fineFaceID );

         for ( const auto& neighborVertexID : fineFace->neighborVertices() )
         {
            auto vertex = getVertex( neighborVertexID );
            for ( const auto& neighborFaceID : vertex->neighborFaces() )
            {
               auto nFace = getFace( neighborFaceID );

               if ( neighborFaceID != fineFaceID )
               {
                  if ( !algorithms::contains( fineFace->indirectNeighborFaceIDsOverVertices_, neighborFaceID ) )
                  {
                     fineFace->indirectNeighborFaceIDsOverVertices_.push_back( neighborFaceID );
                     if ( !algorithms::contains( nFace->indirectNeighborFaceIDsOverVertices_, fineFaceID ) )
                     {
                        nFace->indirectNeighborFaceIDsOverVertices_.push_back( fineFaceID );
                     }
                  }
               }
            }
         }

         for ( const auto& neighborEdgeID : fineFace->neighborEdges() )
         {
            auto edge = getEdge( neighborEdgeID );
            for ( const auto& neighborFaceID : edge->neighborFaces() )
            {
               auto nFace = getFace( neighborFaceID );
               if ( neighborFaceID != fineFaceID )
               {
                  WALBERLA_ASSERT( algorithms::contains( nFace->neighborEdges(), neighborEdgeID ) );

                  fineFace->indirectNeighborFaceIDsOverEdges_[fineFace->edge_index( neighborEdgeID )] = neighborFaceID;
                  nFace->indirectNeighborFaceIDsOverEdges_[nFace->edge_index( neighborEdgeID )]       = fineFaceID;
               }
            }
         }
      }
   }

   // Top-level neighborhood.
   // This is fully reset, and we now rebuild it for all top-level volumes.
   for ( const auto& fineFaceID : getFaceIDs() )
   {
      auto fineFace = getFace( fineFaceID );

      WALBERLA_CHECK( !fineFace->hasChildren(), "This face was just refined and therefore should not have any children." );

      // Top-level neighbors over vertices.
      for ( const auto& neighborVertexID : fineFace->neighborVertices() )
      {
         // Go over neighborhood of this vertex, its child and parent if available.
         auto vertexMid = getVertex( neighborVertexID );

         // Finer vertices
         if ( vertexMid->hasChildren() )
         {
            auto vertexFineID = vertexMid->childVertices()[0];
            WALBERLA_CHECK( vertexExistsLocally( vertexFineID ) || vertexExistsInNeighborhood( vertexFineID ) );
            auto vertexFine = getVertex( vertexFineID );

            for ( const auto& neighborFaceID : vertexFine->neighborFaces() )
            {
               WALBERLA_CHECK( neighborFaceID != fineFaceID );
               auto nFace = getFace( neighborFaceID );
               WALBERLA_CHECK( !nFace->hasChildren(), "This would violate 2:1 balance." );

               if ( !algorithms::contains( fineFace->indirectTopLevelNeighborFaceIDsOverVertices_, neighborFaceID ) )
               {
                  fineFace->indirectTopLevelNeighborFaceIDsOverVertices_.push_back( neighborFaceID );
               }

               if ( !algorithms::contains( nFace->indirectTopLevelNeighborFaceIDsOverVertices_, fineFaceID ) )
               {
                  nFace->indirectTopLevelNeighborFaceIDsOverVertices_.push_back( fineFaceID );
               }
            }
         }

         // Equal level
         for ( const auto& neighborFaceID : vertexMid->neighborFaces() )
         {
            auto nFace = getFace( neighborFaceID );
            if ( neighborFaceID != fineFaceID && !nFace->hasChildren() )
            {
               if ( !algorithms::contains( fineFace->indirectTopLevelNeighborFaceIDsOverVertices_, neighborFaceID ) )
               {
                  fineFace->indirectTopLevelNeighborFaceIDsOverVertices_.push_back( neighborFaceID );
               }

               if ( !algorithms::contains( nFace->indirectTopLevelNeighborFaceIDsOverVertices_, fineFaceID ) )
               {
                  nFace->indirectTopLevelNeighborFaceIDsOverVertices_.push_back( fineFaceID );
               }
            }
         }

         // Coarser vertex
         if ( vertexMid->hasParent() )
         {
            auto vertexCoarseID = vertexMid->parent();

            // The parent of the neighbor vertex of the newly created fine face might not be a vertex.
            if ( vertexExistsLocally( vertexCoarseID ) || vertexExistsInNeighborhood( vertexCoarseID ) )
            {
               auto vertexCoarse = getVertex( vertexCoarseID );

               for ( const auto& neighborFaceID : vertexCoarse->neighborFaces() )
               {
                  WALBERLA_CHECK( neighborFaceID != fineFaceID );
                  auto nFace = getFace( neighborFaceID );
                  if ( !nFace->hasChildren() )
                  {
                     if ( !algorithms::contains( fineFace->indirectTopLevelNeighborFaceIDsOverVertices_, neighborFaceID ) )
                     {
                        fineFace->indirectTopLevelNeighborFaceIDsOverVertices_.push_back( neighborFaceID );
                     }

                     if ( !algorithms::contains( nFace->indirectTopLevelNeighborFaceIDsOverVertices_, fineFaceID ) )
                     {
                        nFace->indirectTopLevelNeighborFaceIDsOverVertices_.push_back( fineFaceID );
                     }
                  }
               }
            }
         }
      }

      // Top-level neighbors over edges.

      for ( const auto& neighborEdgeID : fineFace->neighborEdges() )
      {
         // Go over neighborhood of this edge, its child and parent if available.
         auto edgeMid = getEdge( neighborEdgeID );

         auto localEdgeID = fineFace->edge_index( neighborEdgeID );

         // Finer edges
         if ( edgeMid->hasChildren() )
         {
            for ( const auto& childEdgeID : edgeMid->childEdges() )
            {
               auto edgeFine = getEdge( childEdgeID );
               for ( const auto& neighborFaceID : edgeFine->neighborFaces() )
               {
                  WALBERLA_CHECK( neighborFaceID != fineFaceID );
                  auto nFace = getFace( neighborFaceID );
                  WALBERLA_CHECK( !nFace->hasChildren(), "This would violate 2:1 balance." );

                  if ( !algorithms::contains( fineFace->indirectTopLevelNeighborFaceIDsOverEdges_[localEdgeID], neighborFaceID ) )
                  {
                     fineFace->indirectTopLevelNeighborFaceIDsOverEdges_[localEdgeID].push_back( neighborFaceID );
                  }

                  auto nLocalEdgeID = nFace->edge_index( childEdgeID );
                  if ( !algorithms::contains( nFace->indirectTopLevelNeighborFaceIDsOverEdges_[nLocalEdgeID], fineFaceID ) )
                  {
                     nFace->indirectTopLevelNeighborFaceIDsOverEdges_[nLocalEdgeID].push_back( fineFaceID );
                  }
               }
            }
         }

         // Equal level
         for ( const auto& neighborFaceID : edgeMid->neighborFaces() )
         {
            auto nFace = getFace( neighborFaceID );
            if ( neighborFaceID != fineFaceID && !nFace->hasChildren() )
            {
               if ( !algorithms::contains( fineFace->indirectTopLevelNeighborFaceIDsOverEdges_[localEdgeID], neighborFaceID ) )
               {
                  fineFace->indirectTopLevelNeighborFaceIDsOverEdges_[localEdgeID].push_back( neighborFaceID );
               }

               auto nLocalEdgeID = nFace->edge_index( neighborEdgeID );
               if ( !algorithms::contains( nFace->indirectTopLevelNeighborFaceIDsOverEdges_[nLocalEdgeID], fineFaceID ) )
               {
                  nFace->indirectTopLevelNeighborFaceIDsOverEdges_[nLocalEdgeID].push_back( fineFaceID );
               }
            }
         }

         // Coarser edge
         if ( edgeMid->hasParent() )
         {
            auto edgeCoarseID = edgeMid->parent();

            // The parent of the neighbor vertex of the newly created fine face might not be a vertex.
            if ( edgeExistsLocally( edgeCoarseID ) || edgeExistsInNeighborhood( edgeCoarseID ) )
            {
               auto edgeCoarse = getEdge( edgeCoarseID );

               for ( const auto& neighborFaceID : edgeCoarse->neighborFaces() )
               {
                  WALBERLA_CHECK( neighborFaceID != fineFaceID );
                  auto nFace = getFace( neighborFaceID );
                  if ( !nFace->hasChildren() )
                  {
                     if ( !algorithms::contains( fineFace->indirectTopLevelNeighborFaceIDsOverEdges_[localEdgeID],
                                                 neighborFaceID ) )
                     {
                        fineFace->indirectTopLevelNeighborFaceIDsOverEdges_[localEdgeID].push_back( neighborFaceID );
                     }

                     auto nLocalEdgeID = nFace->edge_index( edgeCoarseID );
                     if ( !algorithms::contains( nFace->indirectTopLevelNeighborFaceIDsOverEdges_[nLocalEdgeID], fineFaceID ) )
                     {
                        nFace->indirectTopLevelNeighborFaceIDsOverEdges_[nLocalEdgeID].push_back( fineFaceID );
                     }
                  }
               }
            }
         }
      }
   }

   wasModified();
   updateLeafPrimitiveMaps();
}
void PrimitiveStorage::updateLeafPrimitiveMaps()
{
   leafVertices_.clear();
   for ( const auto& [level, primitives] : vertices_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            leafVertices_[pid] = primitive;
         }
      }
      WALBERLA_UNUSED( level );
   }

   leafEdges_.clear();
   for ( const auto& [level, primitives] : edges_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            leafEdges_[pid] = primitive;
         }
      }
      WALBERLA_UNUSED( level );
   }

   leafFaces_.clear();
   for ( const auto& [level, primitives] : faces_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            leafFaces_[pid] = primitive;
         }
      }
      WALBERLA_UNUSED( level );
   }

   leafCells_.clear();
   for ( const auto& [level, primitives] : cells_ )
   {
      for ( const auto& [pid, primitive] : primitives )
      {
         if ( !primitive->hasChildren() )
         {
            leafCells_[pid] = primitive;
         }
      }
      WALBERLA_UNUSED( level );
   }
}

} // namespace hyteg

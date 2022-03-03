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
#include "core/mpi/Gatherv.h"
#include "core/mpi/OpenMPBufferSystem.h"
#include "core/math/DistributedSample.h"

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
   for ( auto it : setupStorage.getVertices() )
   {
      if ( uint_c( walberla::mpi::MPIManager::instance()->rank() ) == setupStorage.getTargetRank( it.first ) )
      {
         vertices_[it.first] = std::make_shared< Vertex >( *( it.second ) );
      }
   }

   for ( auto it : setupStorage.getEdges() )
   {
      if ( uint_c( walberla::mpi::MPIManager::instance()->rank() ) == setupStorage.getTargetRank( it.first ) )
      {
         edges_[it.first] = std::make_shared< Edge >( *( it.second ) );
      }
   }

   for ( auto it : setupStorage.getFaces() )
   {
      if ( uint_c( walberla::mpi::MPIManager::instance()->rank() ) == setupStorage.getTargetRank( it.first ) )
      {
         faces_[it.first] = std::make_shared< Face >( *( it.second ) );
      }
   }

   for ( auto it : setupStorage.getCells() )
   {
      if ( uint_c( walberla::mpi::MPIManager::instance()->rank() ) == setupStorage.getTargetRank( it.first ) )
      {
         cells_[it.first] = std::make_shared< Cell >( *( it.second ) );
      }
   }

   // neighbors
   std::vector< PrimitiveID > vertices{ getVertexIDs() };
   std::vector< PrimitiveID > edges{ getEdgeIDs() };
   std::vector< PrimitiveID > faces{ getFaceIDs() };
   std::vector< PrimitiveID > cells{ getCellIDs() };
   addDirectNeighbors( setupStorage, vertices, edges, faces, cells );

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
            neighborVertices_[neighborVertexID.getID()] = std::make_shared< Vertex >( *neighborVertex );
            neighborRanks_[neighborVertexID.getID()]    = setupStorage.getTargetRank( neighborVertexID.getID() );
         }
      }

      for ( const auto& neighborEdgeID : vertex->neighborEdges() )
      {
         const Edge* neighborEdge = setupStorage.getEdge( neighborEdgeID );
         if ( !edgeExistsLocally( neighborEdgeID ) && !edgeExistsInNeighborhood( neighborEdgeID ) )
         {
            neighborEdges_[neighborEdgeID.getID()] = std::make_shared< Edge >( *neighborEdge );
            neighborRanks_[neighborEdgeID.getID()] = setupStorage.getTargetRank( neighborEdgeID.getID() );
         }
      }

      for ( const auto& neighborFaceID : vertex->neighborFaces() )
      {
         const Face* neighborFace = setupStorage.getFace( neighborFaceID );
         if ( !faceExistsLocally( neighborFaceID ) && !faceExistsInNeighborhood( neighborFaceID ) )
         {
            neighborFaces_[neighborFaceID.getID()] = std::make_shared< Face >( *neighborFace );
            neighborRanks_[neighborFaceID.getID()] = setupStorage.getTargetRank( neighborFaceID.getID() );
         }
      }

      for ( const auto& neighborCellID : vertex->neighborCells() )
      {
         const Cell* neighborCell = setupStorage.getCell( neighborCellID );
         if ( !cellExistsLocally( neighborCellID ) && !cellExistsInNeighborhood( neighborCellID ) )
         {
            neighborCells_[neighborCellID.getID()] = std::make_shared< Cell >( *neighborCell );
            neighborRanks_[neighborCellID.getID()] = setupStorage.getTargetRank( neighborCellID.getID() );
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
            neighborVertices_[neighborVertexID.getID()] = std::make_shared< Vertex >( *neighborVertex );
            neighborRanks_[neighborVertexID.getID()]    = setupStorage.getTargetRank( neighborVertexID.getID() );
         }
      }

      for ( const auto& neighborEdgeID : edge->neighborEdges() )
      {
         const Edge* neighborEdge = setupStorage.getEdge( neighborEdgeID );
         if ( !edgeExistsLocally( neighborEdgeID ) && !edgeExistsInNeighborhood( neighborEdgeID ) )
         {
            neighborEdges_[neighborEdgeID.getID()] = std::make_shared< Edge >( *neighborEdge );
            neighborRanks_[neighborEdgeID.getID()] = setupStorage.getTargetRank( neighborEdgeID.getID() );
         }
      }

      for ( const auto& neighborFaceID : edge->neighborFaces() )
      {
         const Face* neighborFace = setupStorage.getFace( neighborFaceID );
         if ( !faceExistsLocally( neighborFaceID ) && !faceExistsInNeighborhood( neighborFaceID ) )
         {
            neighborFaces_[neighborFaceID.getID()] = std::make_shared< Face >( *neighborFace );
            neighborRanks_[neighborFaceID.getID()] = setupStorage.getTargetRank( neighborFaceID.getID() );
         }
      }

      for ( const auto& neighborCellID : edge->neighborCells() )
      {
         const Cell* neighborCell = setupStorage.getCell( neighborCellID );
         if ( !cellExistsLocally( neighborCellID ) && !cellExistsInNeighborhood( neighborCellID ) )
         {
            neighborCells_[neighborCellID.getID()] = std::make_shared< Cell >( *neighborCell );
            neighborRanks_[neighborCellID.getID()] = setupStorage.getTargetRank( neighborCellID.getID() );
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
            neighborVertices_[neighborVertexID.getID()] = std::make_shared< Vertex >( *neighborVertex );
            neighborRanks_[neighborVertexID.getID()]    = setupStorage.getTargetRank( neighborVertexID.getID() );
         }
      }

      for ( const auto& neighborEdgeID : face->neighborEdges() )
      {
         const Edge* neighborEdge = setupStorage.getEdge( neighborEdgeID );
         if ( !edgeExistsLocally( neighborEdgeID ) && !edgeExistsInNeighborhood( neighborEdgeID ) )
         {
            neighborEdges_[neighborEdgeID.getID()] = std::make_shared< Edge >( *neighborEdge );
            neighborRanks_[neighborEdgeID.getID()] = setupStorage.getTargetRank( neighborEdgeID.getID() );
         }
      }

      for ( const auto& neighborFaceID : face->neighborFaces() )
      {
         const Face* neighborFace = setupStorage.getFace( neighborFaceID );
         if ( !faceExistsLocally( neighborFaceID ) && !faceExistsInNeighborhood( neighborFaceID ) )
         {
            neighborFaces_[neighborFaceID.getID()] = std::make_shared< Face >( *neighborFace );
            neighborRanks_[neighborFaceID.getID()] = setupStorage.getTargetRank( neighborFaceID.getID() );
         }
      }

      for ( const auto& neighborCellID : face->neighborCells() )
      {
         const Cell* neighborCell = setupStorage.getCell( neighborCellID );
         if ( !cellExistsLocally( neighborCellID ) && !cellExistsInNeighborhood( neighborCellID ) )
         {
            neighborCells_[neighborCellID.getID()] = std::make_shared< Cell >( *neighborCell );
            neighborRanks_[neighborCellID.getID()] = setupStorage.getTargetRank( neighborCellID.getID() );
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
            neighborVertices_[neighborVertexID.getID()] = std::make_shared< Vertex >( *neighborVertex );
            neighborRanks_[neighborVertexID.getID()]    = setupStorage.getTargetRank( neighborVertexID.getID() );
         }
      }

      for ( const auto& neighborEdgeID : cell->neighborEdges() )
      {
         const Edge* neighborEdge = setupStorage.getEdge( neighborEdgeID );
         if ( !edgeExistsLocally( neighborEdgeID ) && !edgeExistsInNeighborhood( neighborEdgeID ) )
         {
            neighborEdges_[neighborEdgeID.getID()] = std::make_shared< Edge >( *neighborEdge );
            neighborRanks_[neighborEdgeID.getID()] = setupStorage.getTargetRank( neighborEdgeID.getID() );
         }
      }

      for ( const auto& neighborFaceID : cell->neighborFaces() )
      {
         const Face* neighborFace = setupStorage.getFace( neighborFaceID );
         if ( !faceExistsLocally( neighborFaceID ) && !faceExistsInNeighborhood( neighborFaceID ) )
         {
            neighborFaces_[neighborFaceID.getID()] = std::make_shared< Face >( *neighborFace );
            neighborRanks_[neighborFaceID.getID()] = setupStorage.getTargetRank( neighborFaceID.getID() );
         }
      }

      for ( const auto& neighborCellID : cell->neighborCells() )
      {
         const Cell* neighborCell = setupStorage.getCell( neighborCellID );
         if ( !cellExistsLocally( neighborCellID ) && !cellExistsInNeighborhood( neighborCellID ) )
         {
            neighborCells_[neighborCellID.getID()] = std::make_shared< Cell >( *neighborCell );
            neighborRanks_[neighborCellID.getID()] = setupStorage.getTargetRank( neighborCellID.getID() );
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
   std::set< MPIRank > nranks;
   for ( const auto & it : getNeighboringRanks() )
   {
      nranks.insert( static_cast< MPIRank >( it ) );
   }

   bs.setReceiverInfo( nranks, true );

   for ( auto nbrank : getNeighboringRanks() )
   {
      bs.sendBuffer( nbrank ) << getNumberOfLocalVertices();
      for ( const auto & it : getVertices() )
      {
         bs.sendBuffer( nbrank ) << it.first;
         bs.sendBuffer( nbrank ) << *it.second;
      }

      bs.sendBuffer( nbrank ) << getNumberOfLocalEdges();
      for ( const auto & it : getEdges() )
      {
         bs.sendBuffer( nbrank ) << it.first;
         bs.sendBuffer( nbrank ) << *it.second;
      }

      bs.sendBuffer( nbrank ) << getNumberOfLocalFaces();
      for ( const auto & it : getFaces() )
      {
         bs.sendBuffer( nbrank ) << it.first;
         bs.sendBuffer( nbrank ) << *it.second;
      }

      bs.sendBuffer( nbrank ) << getNumberOfLocalCells();
      for ( const auto & it : getCells() )
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
         PrimitiveID::IDType id;
         msg.buffer() >> id;
         neighborVertices_[id] = std::make_shared< Vertex >( msg.buffer() );
         neighborRanks_[id] = uint_c( nbrank );
      }

      uint_t numEdges;
      msg.buffer() >> numEdges;
      for ( uint_t i = 0; i < numEdges; i++ )
      {
         PrimitiveID::IDType id;
         msg.buffer() >> id;
         neighborEdges_[id] = std::make_shared< Edge >( msg.buffer() );
         neighborRanks_[id] = uint_c( nbrank );
      }

      uint_t numFaces;
      msg.buffer() >> numFaces;
      for ( uint_t i = 0; i < numFaces; i++ )
      {
         PrimitiveID::IDType id;
         msg.buffer() >> id;
         neighborFaces_[id] = std::make_shared< Face >( msg.buffer() );
         neighborRanks_[id] = uint_c( nbrank );
      }

      uint_t numCells;
      msg.buffer() >> numCells;
      for ( uint_t i = 0; i < numCells; i++ )
      {
         PrimitiveID::IDType id;
         msg.buffer() >> id;
         neighborCells_[id] = std::make_shared< Cell >( msg.buffer() );
         neighborRanks_[id] = uint_c( nbrank );
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
: vertices_( vtxs )
, edges_( edges )
, faces_( faces )
, cells_( cells )
, neighborVertices_( nbrvtxs )
, neighborEdges_( nbredges )
, neighborFaces_( nbrfaces )
, neighborCells_( nbrcells )
, primitiveDataHandlers_( 0 )
, neighborRanks_( neighborRanks )
, modificationStamp_( 0 )
, timingTree_( std::make_shared< walberla::WcTimingTree >() )
, hasGlobalCells_( hasGlobalCells )
, additionalHaloDepth_( 0 )
{
   splitCommunicatorByPrimitiveDistribution();

#ifndef NDEBUG
   checkConsistency();
#endif
}

std::shared_ptr< PrimitiveStorage > PrimitiveStorage::createCopy() const
{
   auto copiedStorage = std::make_shared< PrimitiveStorage >(
       SetupPrimitiveStorage( MeshInfo::emptyMeshInfo(), uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) ) );

   for ( const auto& it : getVertices() )
   {
      copiedStorage->vertices_[it.first] = std::make_shared< Vertex >( *it.second );
   }

   for ( const auto& it : getEdges() )
   {
      copiedStorage->edges_[it.first] = std::make_shared< Edge >( *it.second );
   }

   for ( const auto& it : getFaces() )
   {
      copiedStorage->faces_[it.first] = std::make_shared< Face >( *it.second );
   }

   for ( const auto& it : getCells() )
   {
      copiedStorage->cells_[it.first] = std::make_shared< Cell >( *it.second );
   }

   for ( const auto& it : neighborVertices_ )
   {
      copiedStorage->neighborVertices_[it.first] = std::make_shared< Vertex >( *it.second );
   }

   for ( const auto& it : neighborEdges_ )
   {
      copiedStorage->neighborEdges_[it.first] = std::make_shared< Edge >( *it.second );
   }

   for ( const auto& it : neighborFaces_ )
   {
      copiedStorage->neighborFaces_[it.first] = std::make_shared< Face >( *it.second );
   }

   for ( const auto& it : neighborCells_ )
   {
      copiedStorage->neighborCells_[it.first] = std::make_shared< Cell >( *it.second );
   }

   copiedStorage->neighborRanks_  = neighborRanks_;
   copiedStorage->hasGlobalCells_ = hasGlobalCells_;
   copiedStorage->splitComm_      = splitComm_;

   return copiedStorage;
}

void PrimitiveStorage::getPrimitives( PrimitiveMap& primitiveMap ) const
{
   primitiveMap.clear();

   primitiveMap.insert( vertices_.begin(), vertices_.end() );
   primitiveMap.insert( edges_.begin(), edges_.end() );
   primitiveMap.insert( faces_.begin(), faces_.end() );
   primitiveMap.insert( cells_.begin(), cells_.end() );

   WALBERLA_ASSERT_EQUAL( primitiveMap.size(), vertices_.size() + edges_.size() + faces_.size() + cells_.size() );
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
      return vertices_.at( id.getID() ).get();
   }
   else if ( vertexExistsInNeighborhood( id ) )
   {
      return neighborVertices_.at( id.getID() ).get();
   }
   else
   {
      return nullptr;
   }
}

Vertex* PrimitiveStorage::getVertex( const PrimitiveID& id )
{
   if ( vertexExistsLocally( id ) )
   {
      return vertices_[id.getID()].get();
   }
   else if ( vertexExistsInNeighborhood( id ) )
   {
      return neighborVertices_[id.getID()].get();
   }
   else
   {
      return nullptr;
   }
}

const Edge* PrimitiveStorage::getEdge( const PrimitiveID& id ) const
{
   if ( edgeExistsLocally( id ) )
   {
      return edges_.at( id.getID() ).get();
   }
   else if ( edgeExistsInNeighborhood( id ) )
   {
      return neighborEdges_.at( id.getID() ).get();
   }
   else
   {
      return nullptr;
   }
}

Edge* PrimitiveStorage::getEdge( const PrimitiveID& id )
{
   if ( edgeExistsLocally( id ) )
   {
      return edges_[id.getID()].get();
   }
   else if ( edgeExistsInNeighborhood( id ) )
   {
      return neighborEdges_[id.getID()].get();
   }
   else
   {
      return nullptr;
   }
}

const Face* PrimitiveStorage::getFace( const PrimitiveID& id ) const
{
   if ( faceExistsLocally( id ) )
   {
      return faces_.at( id.getID() ).get();
   }
   else if ( faceExistsInNeighborhood( id ) )
   {
      return neighborFaces_.at( id.getID() ).get();
   }
   else
   {
      return nullptr;
   }
}

Face* PrimitiveStorage::getFace( const PrimitiveID& id )
{
   if ( faceExistsLocally( id ) )
   {
      return faces_[id.getID()].get();
   }
   else if ( faceExistsInNeighborhood( id ) )
   {
      return neighborFaces_[id.getID()].get();
   }
   else
   {
      return nullptr;
   }
}

const Cell* PrimitiveStorage::getCell( const PrimitiveID& id ) const
{
   if ( cellExistsLocally( id ) )
   {
      return cells_.at( id.getID() ).get();
   }
   else if ( cellExistsInNeighborhood( id ) )
   {
      return neighborCells_.at( id.getID() ).get();
   }
   else
   {
      return nullptr;
   }
}

Cell* PrimitiveStorage::getCell( const PrimitiveID& id )
{
   if ( cellExistsLocally( id ) )
   {
      return cells_[id.getID()].get();
   }
   else if ( cellExistsInNeighborhood( id ) )
   {
      return neighborCells_[id.getID()].get();
   }
   else
   {
      return nullptr;
   }
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
   for ( auto const& it : vertices_ )
   {
      vertexIDs.push_back( it.first );
   }
}

void PrimitiveStorage::getNeighboringVertexIDs( std::vector< PrimitiveID >& vertexIDs ) const
{
   vertexIDs.clear();
   for ( auto const& it : neighborVertices_ )
   {
      vertexIDs.push_back( it.first );
   }
}

void PrimitiveStorage::getEdgeIDs( std::vector< PrimitiveID >& edgeIDs ) const
{
   edgeIDs.clear();
   for ( auto const& it : edges_ )
   {
      edgeIDs.push_back( it.first );
   }
}

void PrimitiveStorage::getNeighboringEdgeIDs( std::vector< PrimitiveID >& edgeIDs ) const
{
   edgeIDs.clear();
   for ( auto const& it : neighborEdges_ )
   {
      edgeIDs.push_back( it.first );
   }
}

void PrimitiveStorage::getFaceIDs( std::vector< PrimitiveID >& faceIDs ) const
{
   faceIDs.clear();
   for ( auto const& it : faces_ )
   {
      faceIDs.push_back( it.first );
   }
}

void PrimitiveStorage::getNeighboringFaceIDs( std::vector< PrimitiveID >& faceIDs ) const
{
   faceIDs.clear();
   for ( auto const& it : neighborFaces_ )
   {
      faceIDs.push_back( it.first );
   }
}

void PrimitiveStorage::getCellIDs( std::vector< PrimitiveID >& cellIDs ) const
{
   cellIDs.clear();
   for ( auto const& it : cells_ )
   {
      cellIDs.push_back( it.first );
   }
}

void PrimitiveStorage::getNeighboringCellIDs( std::vector< PrimitiveID >& cellIDs ) const
{
   cellIDs.clear();
   for ( auto const& it : neighborCells_ )
   {
      cellIDs.push_back( it.first );
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
   WALBERLA_DEBUG_SECTION() { checkConsistency(); }

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
      if ( primitivesToMigrate.count( pID.getID() ) > 0 )
      {
         localPrimitiveToFutureRankMap[pID.getID()] = primitivesToMigrate.at( pID.getID() );
      }
      else
      {
         localPrimitiveToFutureRankMap[pID.getID()] = rank;
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

      WALBERLA_CHECK_NOT_IDENTICAL( primitiveType, INVALID, "Sending invalid primitive type..." );

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
         WALBERLA_ASSERT_EQUAL( neighborhoodPrimitiveToFutureRankMap.count( neighborID.getID() ), 1 );
         sendBuffer << neighborhoodPrimitiveToFutureRankMap[neighborID.getID()];
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
            neighborhoodPrimitiveToFutureRankMap[neighborPrimitiveID.getID()] = rankAfterMigration;
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
            neighborVertices_[idToErase.getID()] = std::make_shared< Vertex >( *getVertex( idToErase ) );
            vertices_.erase( idToErase.getID() );
         }
         if ( edgeExistsLocally( idToErase ) )
         {
            neighborEdges_[idToErase.getID()] = std::make_shared< Edge >( *getEdge( idToErase ) );
            edges_.erase( idToErase.getID() );
         }
         if ( faceExistsLocally( idToErase ) )
         {
            neighborFaces_[idToErase.getID()] = std::make_shared< Face >( *getFace( idToErase ) );
            faces_.erase( idToErase.getID() );
         }
         if ( cellExistsLocally( idToErase ) )
         {
            neighborCells_[idToErase.getID()] = std::make_shared< Cell >( *getCell( idToErase ) );
            cells_.erase( idToErase.getID() );
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
         neighborVertices_.erase( localID.getID() );
      if ( edgeExistsInNeighborhood( localID ) )
         neighborEdges_.erase( localID.getID() );
      if ( faceExistsInNeighborhood( localID ) )
         neighborFaces_.erase( localID.getID() );
      if ( cellExistsInNeighborhood( localID ) )
         neighborCells_.erase( localID.getID() );
   }

   /////////////////////////////////////////////////////////////////////////////////////////
   // Erasing all neighbors that are not referenced by local primitives from neighborhood //
   /////////////////////////////////////////////////////////////////////////////////////////

   // enjoy worst possible complexity...

   std::vector< PrimitiveID > neighborhoodIDs;
   for ( const auto& it : neighborVertices_ )
      neighborhoodIDs.push_back( it.first );
   for ( const auto& it : neighborEdges_ )
      neighborhoodIDs.push_back( it.first );
   for ( const auto& it : neighborFaces_ )
      neighborhoodIDs.push_back( it.first );
   for ( const auto& it : neighborCells_ )
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
            neighborVertices_.erase( neighborhoodID.getID() );
         if ( edgeExistsInNeighborhood( neighborhoodID ) )
            neighborEdges_.erase( neighborhoodID.getID() );
         if ( faceExistsInNeighborhood( neighborhoodID ) )
            neighborFaces_.erase( neighborhoodID.getID() );
         if ( cellExistsInNeighborhood( neighborhoodID ) )
            neighborCells_.erase( neighborhoodID.getID() );
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
      WALBERLA_CHECK_GREATER( neighborhoodPrimitiveToFutureRankMap.count( np.getID() ), 0 );
      neighborRanks_[np.getID()] = neighborhoodPrimitiveToFutureRankMap[np.getID()];
   }

   for ( uint_t i = 0; i < additionalHaloDepth_; i++ )
   {
      addDirectNeighborsDistributed();
   }

   splitCommunicatorByPrimitiveDistribution();

   wasModified();

   WALBERLA_DEBUG_SECTION() { checkConsistency(); }
}

std::set< uint_t > PrimitiveStorage::getNeighboringFaceRanksOfFace( const PrimitiveID & facePrimitiveID ) const
{
   WALBERLA_CHECK( additionalHaloDepth_ > 0, "Additional halo depth must be larger than 0 to access neighbor faces of faces." );
   WALBERLA_CHECK( faceExistsLocally( facePrimitiveID ), "Face " << facePrimitiveID << "does not exist locally." );

   std::set< uint_t > neighboringRanks;
   const auto face = getFace( facePrimitiveID );
   for ( const auto & npid : face->getIndirectNeighborFaceIDs() )
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
   for ( const auto & it : faces_ )
   {
      const auto nr = getNeighboringFaceRanksOfFace( it.first );
      neighborRanks.insert( nr.begin(), nr.end() );
   }
   return neighborRanks;
}


std::set< uint_t > PrimitiveStorage::getNeighboringCellRanksOfCell( const PrimitiveID & cellPrimitiveID ) const
{
   WALBERLA_CHECK( additionalHaloDepth_ > 0, "Additional halo depth must be larger than 0 to access neighbor cells of cells." );
   WALBERLA_CHECK( cellExistsLocally( cellPrimitiveID ), "Cell " << cellPrimitiveID << "does not exist locally." );

   std::set< uint_t > neighboringRanks;
   const auto cell = getCell( cellPrimitiveID );
   for ( const auto & npid : cell->getIndirectNeighborCellIDs() )
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
   for ( const auto & it : cells_ )
   {
      const auto nr = getNeighboringCellRanksOfCell( it.first );
      neighborRanks.insert( nr.begin(), nr.end() );
   }
   return neighborRanks;
}

std::set< uint_t > PrimitiveStorage::getNeighboringVolumeRanksOfVolume( const PrimitiveID & volumePrimitiveID ) const
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
   for ( const auto& it : neighborRanks_ )
   {
      const uint_t neighborRank = it.second;
      neighboringRanks.insert( neighborRank );
   }
}

void PrimitiveStorage::getNeighboringRanks( std::set< walberla::mpi::MPIRank >& neighboringRanks ) const
{
   neighboringRanks.clear();
   for ( const auto& it : neighborRanks_ )
   {
      const walberla::mpi::MPIRank neighborRank = static_cast< walberla::mpi::MPIRank >( it.second );
      neighboringRanks.insert( neighborRank );
   }
}

PrimitiveStorage::PrimitiveTypeEnum PrimitiveStorage::getPrimitiveType( const PrimitiveID& primitiveID ) const
{
   if ( vertexExistsLocally( primitiveID ) || vertexExistsInNeighborhood( primitiveID ) )
      return VERTEX;
   if ( edgeExistsLocally( primitiveID ) || edgeExistsInNeighborhood( primitiveID ) )
      return EDGE;
   if ( faceExistsLocally( primitiveID ) || faceExistsInNeighborhood( primitiveID ) )
      return FACE;
   if ( cellExistsLocally( primitiveID ) || cellExistsInNeighborhood( primitiveID ) )
      return CELL;
   return INVALID;
}

PrimitiveID PrimitiveStorage::deserializeAndAddPrimitive( walberla::mpi::RecvBuffer& recvBuffer, const bool& isNeighborPrimitive )
{
   PrimitiveTypeEnum primitiveType;
   PrimitiveID       primitiveID;

   recvBuffer >> primitiveType;

   switch ( primitiveType )
   {
   case VERTEX:
   {
      std::shared_ptr< Vertex > vertex = std::make_shared< Vertex >( recvBuffer );
      primitiveID                      = vertex->getID();
      if ( isNeighborPrimitive )
      {
         neighborVertices_[primitiveID.getID()] = vertex;
      }
      else
      {
         vertices_[primitiveID.getID()] = vertex;
      }
      break;
   }
   case EDGE:
   {
      std::shared_ptr< Edge > edge = std::make_shared< Edge >( recvBuffer );
      primitiveID                  = edge->getID();
      if ( isNeighborPrimitive )
      {
         neighborEdges_[primitiveID.getID()] = edge;
      }
      else
      {
         edges_[primitiveID.getID()] = edge;
      }
      break;
   }
   case FACE:
   {
      std::shared_ptr< Face > face = std::make_shared< Face >( recvBuffer );
      primitiveID                  = face->getID();
      if ( isNeighborPrimitive )
      {
         neighborFaces_[primitiveID.getID()] = face;
      }
      else
      {
         faces_[primitiveID.getID()] = face;
      }
      break;
   }
   case CELL:
   {
      std::shared_ptr< Cell > cell = std::make_shared< Cell >( recvBuffer );
      primitiveID                  = cell->getID();
      if ( isNeighborPrimitive )
      {
         neighborCells_[primitiveID.getID()] = cell;
      }
      else
      {
         cells_[primitiveID.getID()] = cell;
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
   WALBERLA_ASSERT( primitiveExistsLocally( primitiveID ) );
   const PrimitiveTypeEnum primitiveType = getPrimitiveType( primitiveID );
   switch ( primitiveType )
   {
   case VERTEX:
   {
      WALBERLA_ASSERT( vertexExistsLocally( primitiveID ) );
      auto vertex = vertices_[primitiveID.getID()];
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
   case EDGE:
   {
      WALBERLA_ASSERT( edgeExistsLocally( primitiveID ) );
      auto edge = edges_[primitiveID.getID()];
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
   case FACE:
   {
      WALBERLA_ASSERT( faceExistsLocally( primitiveID ) );
      auto face = faces_[primitiveID.getID()];
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
   case CELL:
   {
      WALBERLA_ASSERT( cellExistsLocally( primitiveID ) );
      auto cell = cells_[primitiveID.getID()];
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
   WALBERLA_ASSERT( primitiveExistsLocally( primitiveID ) );
   const PrimitiveTypeEnum primitiveType = getPrimitiveType( primitiveID );
   switch ( primitiveType )
   {
   case VERTEX:
   {
      WALBERLA_ASSERT( vertexExistsLocally( primitiveID ) );
      auto vertex = vertices_[primitiveID.getID()];
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
   case EDGE:
   {
      WALBERLA_ASSERT( edgeExistsLocally( primitiveID ) );
      auto edge = edges_[primitiveID.getID()];
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
   case FACE:
   {
      WALBERLA_ASSERT( faceExistsLocally( primitiveID ) );
      auto face = faces_[primitiveID.getID()];
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
   case CELL:
   {
      WALBERLA_ASSERT( cellExistsLocally( primitiveID ) );
      auto cell = cells_[primitiveID.getID()];
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
   const auto numNeighborProcesses = getNeighboringRanks().size();
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
   os << " - neighbor processes (min, max, avg):               " << neighborhoodSample.min() << ", " << neighborhoodSample.max() << ", " << neighborhoodSample.avg() << "\n";
   if ( additionalHaloDepth_ > 0 )
   {
      os << " - neighbor processes, volumes only (min, max, avg): " << neighborhoodVolumeSample.min() << ", " << neighborhoodVolumeSample.max() << ", " << neighborhoodVolumeSample.avg() << "\n";
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
   for ( auto it = vertices_.begin(); it != vertices_.end(); it++ )
   {
      WALBERLA_CHECK_GREATER_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
      WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
      WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 0 );
   }
   for ( auto it = edges_.begin(); it != edges_.end(); it++ )
   {
      WALBERLA_CHECK_GREATER_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
      WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
      WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 2 );
   }
   for ( auto it = faces_.begin(); it != faces_.end(); it++ )
   {
      WALBERLA_CHECK_GREATER_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
      WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
      WALBERLA_CHECK_EQUAL( it->second->getNumNeighborVertices(), 3 );
      WALBERLA_CHECK_EQUAL( it->second->getNumNeighborEdges(), 3 );
   }
   for ( auto it = cells_.begin(); it != cells_.end(); it++ )
   {
      WALBERLA_CHECK_GREATER_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
      WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
      WALBERLA_CHECK_EQUAL( it->second->getNumNeighborVertices(), 4 );
      WALBERLA_CHECK_EQUAL( it->second->getNumNeighborEdges(), 6 );
      WALBERLA_CHECK_EQUAL( it->second->getNumNeighborFaces(), 4 );
   }

   // 4. Number of data entries of neighbor primitives is zero
   // 5. PrimitiveIDs of neighbor maps match IDs of neighbor Primitives
   // 6. Neighborhood of Primitives
   for ( auto it = neighborVertices_.begin(); it != neighborVertices_.end(); it++ )
   {
      WALBERLA_CHECK_GREATER_EQUAL( 0, it->second->getNumberOfDataEntries() );
      WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
      WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 0 );
   }
   for ( auto it = neighborEdges_.begin(); it != neighborEdges_.end(); it++ )
   {
      WALBERLA_CHECK_GREATER_EQUAL( 0, it->second->getNumberOfDataEntries() );
      WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
      WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 2 );
   }
   for ( auto it = neighborFaces_.begin(); it != neighborFaces_.end(); it++ )
   {
      WALBERLA_CHECK_GREATER_EQUAL( 0, it->second->getNumberOfDataEntries() );
      WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
      WALBERLA_CHECK_EQUAL( it->second->getNumNeighborVertices(), 3 );
      WALBERLA_CHECK_EQUAL( it->second->getNumNeighborEdges(), 3 );
   }
   for ( auto it = neighborCells_.begin(); it != neighborCells_.end(); it++ )
   {
      WALBERLA_CHECK_GREATER_EQUAL( 0, it->second->getNumberOfDataEntries() );
      WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
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
                         "Neighbor referenced from locally allocated primitive is not allocated in storage." );
      }
   }

   // 9. As many neighbor ranks as neighbors
   WALBERLA_CHECK_EQUAL( neighborRanks_.size(),
                         neighborVertices_.size() + neighborEdges_.size() + neighborFaces_.size() + neighborCells_.size() );

   // 10. Local primitives do not exist in neighborhood
   for ( const auto& id : primitiveIDs )
   {
      WALBERLA_CHECK( !primitiveExistsInNeighborhood( id ), "Primitive that exists in neighborhood: " << id.getID() );
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

} // namespace hyteg

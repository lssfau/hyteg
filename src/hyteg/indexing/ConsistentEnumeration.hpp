/*
 * Copyright (c) 2025 Andreas Burkhart.
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

#pragma once
#include <map>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

struct enumeratorBounds
{
   uint_t maxVertex;
   uint_t maxEdge;
   uint_t maxFace;
   uint_t maxCell;
};

std::pair< std::map< uint_t, PrimitiveID >, enumeratorBounds >
    generateEnumeratorMap( const std::shared_ptr< hyteg::PrimitiveStorage >& storage )
{
   enumeratorBounds                enumBounds;
   std::map< uint_t, PrimitiveID > enumeratorMap;
   uint_t                          enumeratorIndex = 0;

   // loop vertices
   {
      std::vector< PrimitiveID > localVertexIds;

      walberla::mpi::SendBuffer sendbuffer;
      walberla::mpi::RecvBuffer recvbuffer;

      for ( auto& [key, value] : ( storage->getVertices() ) )
      {
         localVertexIds.push_back( key );
      }

      std::vector< PrimitiveID > vertexIds;

      sendbuffer << localVertexIds;

      walberla::mpi::allGathervBuffer( sendbuffer, recvbuffer );

      for ( uint_t rank = 0; rank < uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ); rank++ )
      {
         std::vector< PrimitiveID > receiveIds;
         recvbuffer >> receiveIds;

         for ( auto& id : receiveIds )
         {
            vertexIds.push_back( id );
         }
      }

      std::sort( vertexIds.begin(), vertexIds.end() );

      for ( auto& id : vertexIds )
      {
         enumeratorMap[enumeratorIndex] = id;
         // WALBERLA_LOG_INFO_ON_ROOT( "Vertex: " << enumeratorIndex << " id " << id );
         enumeratorIndex++;
      }
   }
   enumBounds.maxVertex = enumeratorIndex;

   // loop edges
   {
      std::vector< PrimitiveID > localEdgeIds;

      walberla::mpi::SendBuffer sendbuffer;
      walberla::mpi::RecvBuffer recvbuffer;

      for ( auto& [key, value] : ( storage->getEdges() ) )
      {
         localEdgeIds.push_back( key );
      }

      std::vector< PrimitiveID > edgeIds;

      sendbuffer << localEdgeIds;

      walberla::mpi::allGathervBuffer( sendbuffer, recvbuffer );

      for ( uint_t rank = 0; rank < uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ); rank++ )
      {
         std::vector< PrimitiveID > receiveIds;
         recvbuffer >> receiveIds;

         for ( auto& id : receiveIds )
         {
            edgeIds.push_back( id );
         }
      }

      std::sort( edgeIds.begin(), edgeIds.end() );

      for ( auto& id : edgeIds )
      {
         enumeratorMap[enumeratorIndex] = id;
         // WALBERLA_LOG_INFO_ON_ROOT( "Edge: " << enumeratorIndex << " id " << id );
         enumeratorIndex++;
      }
   }
   enumBounds.maxEdge = enumeratorIndex;

   // loop faces
   {
      std::vector< PrimitiveID > localFaceIds;

      walberla::mpi::SendBuffer sendbuffer;
      walberla::mpi::RecvBuffer recvbuffer;

      for ( auto& [key, value] : ( storage->getFaces() ) )
      {
         localFaceIds.push_back( key );
      }

      std::vector< PrimitiveID > faceIds;

      sendbuffer << localFaceIds;

      walberla::mpi::allGathervBuffer( sendbuffer, recvbuffer );

      for ( uint_t rank = 0; rank < uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ); rank++ )
      {
         std::vector< PrimitiveID > receiveIds;
         recvbuffer >> receiveIds;

         for ( auto& id : receiveIds )
         {
            faceIds.push_back( id );
         }
      }

      std::sort( faceIds.begin(), faceIds.end() );

      for ( auto& id : faceIds )
      {
         enumeratorMap[enumeratorIndex] = id;
         // WALBERLA_LOG_INFO_ON_ROOT( "Face: " << enumeratorIndex << " id " << id );
         enumeratorIndex++;
      }
   }
   enumBounds.maxFace = enumeratorIndex;

   // loop cells
   {
      std::vector< PrimitiveID > localCellIds;

      walberla::mpi::SendBuffer sendbuffer;
      walberla::mpi::RecvBuffer recvbuffer;

      for ( auto& [key, value] : ( storage->getCells() ) )
      {
         localCellIds.push_back( key );
      }

      std::vector< PrimitiveID > cellIds;

      sendbuffer << localCellIds;

      walberla::mpi::allGathervBuffer( sendbuffer, recvbuffer );

      for ( uint_t rank = 0; rank < uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ); rank++ )
      {
         std::vector< PrimitiveID > receiveIds;
         recvbuffer >> receiveIds;

         for ( auto& id : receiveIds )
         {
            cellIds.push_back( id );
         }
      }

      std::sort( cellIds.begin(), cellIds.end() );

      for ( auto& id : cellIds )
      {
         enumeratorMap[enumeratorIndex] = id;
         // WALBERLA_LOG_INFO_ON_ROOT( "Cell: " << enumeratorIndex << " id " << id );
         enumeratorIndex++;
      }
   }
   enumBounds.maxCell = enumeratorIndex;

   // WALBERLA_LOG_INFO_ON_ROOT( "maxVertex: " << enumBounds.maxVertex );
   // WALBERLA_LOG_INFO_ON_ROOT( "maxEdge: " << enumBounds.maxEdge );
   // WALBERLA_LOG_INFO_ON_ROOT( "maxFace: " << enumBounds.maxFace );
   // WALBERLA_LOG_INFO_ON_ROOT( "maxCell: " << enumBounds.maxCell );

   return { enumeratorMap, enumBounds };
}

std::pair< std::map< uint_t, PrimitiveID >, enumeratorBounds > generateEnumeratorMap( SetupPrimitiveStorage& setupStorage )
{
   enumeratorBounds                enumBounds;
   std::map< uint_t, PrimitiveID > enumeratorMap;
   uint_t                          enumeratorIndex = 0;

   // loop vertices
   for ( auto& [key, value] : ( setupStorage.getVertices() ) )
   {
      enumeratorMap[enumeratorIndex] = key;
      // WALBERLA_LOG_INFO_ON_ROOT( "Vertex: " << enumeratorIndex << " id " << key );
      enumeratorIndex++;
   }
   enumBounds.maxVertex = enumeratorIndex;

   // loop edges
   for ( auto& [key, value] : ( setupStorage.getEdges() ) )
   {
      enumeratorMap[enumeratorIndex] = key;
      // WALBERLA_LOG_INFO_ON_ROOT( "Edge: " << enumeratorIndex << " id " << key );
      enumeratorIndex++;
   }
   enumBounds.maxEdge = enumeratorIndex;

   // loop faces
   for ( auto& [key, value] : ( setupStorage.getFaces() ) )
   {
      enumeratorMap[enumeratorIndex] = key;
      // WALBERLA_LOG_INFO_ON_ROOT( "Face: " << enumeratorIndex << " id " << key );
      enumeratorIndex++;
   }
   enumBounds.maxFace = enumeratorIndex;

   // loop cells
   for ( auto& [key, value] : ( setupStorage.getCells() ) )
   {
      enumeratorMap[enumeratorIndex] = key;
      // WALBERLA_LOG_INFO_ON_ROOT( "Cell: " << enumeratorIndex << " id " << key );
      enumeratorIndex++;
   }
   enumBounds.maxCell = enumeratorIndex;

   // WALBERLA_LOG_INFO_ON_ROOT( "maxVertex: " << enumBounds.maxVertex );
   // WALBERLA_LOG_INFO_ON_ROOT( "maxEdge: " << enumBounds.maxEdge );
   // WALBERLA_LOG_INFO_ON_ROOT( "maxFace: " << enumBounds.maxFace );
   // WALBERLA_LOG_INFO_ON_ROOT( "maxCell: " << enumBounds.maxCell );

   return { enumeratorMap, enumBounds };
}

void enumerateVertexDoFFunctionConsistently( std::map< uint_t, PrimitiveID >&             enumeratorMap,
                                             enumeratorBounds&                            enumBounds,
                                             const vertexdof::VertexDoFFunction< idx_t >& enumerator,
                                             uint_t                                       level,
                                             idx_t&                                       enumIndex )
{
   auto vertexDataID_   = enumerator.getVertexDataID();
   auto edgeDataID_     = enumerator.getEdgeDataID();
   auto faceDataID_     = enumerator.getFaceDataID();
   auto cellDataID_     = enumerator.getCellDataID();
   auto functionStorage = enumerator.getStorage();

   for ( uint_t i = 0; i < enumBounds.maxVertex; i++ )
   {
      auto Id = enumeratorMap[i];

      if ( functionStorage->vertexExistsLocally( Id ) )
      {
         // WALBERLA_LOG_INFO("found vertex i=" << i);
         vertexdof::macrovertex::enumerate( level, *functionStorage->getVertex( Id ), vertexDataID_, enumIndex );
      }
      walberla::mpi::allReduceInplace( enumIndex, walberla::mpi::MAX );
   }

   for ( uint_t i = enumBounds.maxVertex; i < enumBounds.maxEdge; i++ )
   {
      auto Id = enumeratorMap[i];

      if ( functionStorage->edgeExistsLocally( Id ) )
      {
         // WALBERLA_LOG_INFO("found edge i=" << i);
         vertexdof::macroedge::enumerate< idx_t >( level, *functionStorage->getEdge( Id ), edgeDataID_, enumIndex );
      }
      walberla::mpi::allReduceInplace( enumIndex, walberla::mpi::MAX );
   }

   if ( level >= 2 )
   {
      for ( uint_t i = enumBounds.maxEdge; i < enumBounds.maxFace; i++ )
      {
         auto Id = enumeratorMap[i];

         if ( functionStorage->faceExistsLocally( Id ) )
         {
            // WALBERLA_LOG_INFO("found face i=" << i);
            vertexdof::macroface::enumerate< idx_t >( level, *functionStorage->getFace( Id ), faceDataID_, enumIndex );
         }
         walberla::mpi::allReduceInplace( enumIndex, walberla::mpi::MAX );
      }
   }

   if ( level >= 2 )
   {
      for ( uint_t i = enumBounds.maxFace; i < enumBounds.maxCell; i++ )
      {
         auto Id = enumeratorMap[i];

         if ( functionStorage->cellExistsLocally( Id ) )
         {
            // WALBERLA_LOG_INFO("found cell i=" << i);
            vertexdof::macrocell::enumerate< idx_t >( level, *functionStorage->getCell( Id ), cellDataID_, enumIndex );
         }
         walberla::mpi::allReduceInplace( enumIndex, walberla::mpi::MAX );
      }
   }

   // in contrast to other methods in the function class enumerate needs to communicate due to its usage in the PETSc solvers
   communication::syncFunctionBetweenPrimitives( enumerator, level );
}

void enumerateEdgeDoFFunctionConsistently( std::map< uint_t, PrimitiveID >& enumeratorMap,
                                           enumeratorBounds&                enumBounds,
                                           const EdgeDoFFunction< idx_t >&  enumerator,
                                           uint_t                           level,
                                           idx_t&                           enumIndex )
{
   auto edgeDataID_     = enumerator.getEdgeDataID();
   auto faceDataID_     = enumerator.getFaceDataID();
   auto cellDataID_     = enumerator.getCellDataID();
   auto functionStorage = enumerator.getStorage();

   for ( uint_t i = enumBounds.maxVertex; i < enumBounds.maxEdge; i++ )
   {
      auto Id = enumeratorMap[i];

      if ( functionStorage->edgeExistsLocally( Id ) )
      {
         // WALBERLA_LOG_INFO("found edge i=" << i);
         edgedof::macroedge::enumerate< idx_t >( level, *functionStorage->getEdge( Id ), edgeDataID_, enumIndex );
      }
      walberla::mpi::allReduceInplace( enumIndex, walberla::mpi::MAX );
   }

   if ( level >= 1 )
   {
      for ( uint_t i = enumBounds.maxEdge; i < enumBounds.maxFace; i++ )
      {
         auto Id = enumeratorMap[i];

         if ( functionStorage->faceExistsLocally( Id ) )
         {
            // WALBERLA_LOG_INFO("found face i=" << i);
            edgedof::macroface::enumerate< idx_t >( level, *functionStorage->getFace( Id ), faceDataID_, enumIndex );
         }
         walberla::mpi::allReduceInplace( enumIndex, walberla::mpi::MAX );
      }
   }

   if ( level >= 1 )
   {
      for ( uint_t i = enumBounds.maxFace; i < enumBounds.maxCell; i++ )
      {
         auto Id = enumeratorMap[i];

         if ( functionStorage->cellExistsLocally( Id ) )
         {
            // WALBERLA_LOG_INFO("found cell i=" << i);
            edgedof::macrocell::enumerate< idx_t >( level, *functionStorage->getCell( Id ), cellDataID_, enumIndex );
         }
         walberla::mpi::allReduceInplace( enumIndex, walberla::mpi::MAX );
      }
   }

   communication::syncFunctionBetweenPrimitives( enumerator, level );
}

void enumerateConsistently( std::map< uint_t, PrimitiveID >& enumeratorMap,
                            enumeratorBounds&                enumBounds,
                            P2Function< idx_t >&             enumerator,
                            uint_t                           level )
{
   idx_t enumIndex = 0;
   enumerateVertexDoFFunctionConsistently( enumeratorMap, enumBounds, enumerator.getVertexDoFFunction(), level, enumIndex );
   enumerateEdgeDoFFunctionConsistently( enumeratorMap, enumBounds, enumerator.getEdgeDoFFunction(), level, enumIndex );
}

void enumerateConsistently( std::pair< std::map< uint_t, PrimitiveID >, enumeratorBounds >& enumPair,
                            P2Function< idx_t >&                                            enumerator,
                            uint_t                                                          level )
{
   enumerateConsistently( enumPair.first, enumPair.second, enumerator, level );
}

void enumerateConsistently( P2Function< idx_t >& enumerator, uint_t level )
{
   auto enumPair = generateEnumeratorMap( enumerator.getStorage() );
   enumerateConsistently( enumPair.first, enumPair.second, enumerator, level );
}

void enumerateConsistently( std::map< uint_t, PrimitiveID >& enumeratorMap,
                            enumeratorBounds&                enumBounds,
                            P1Function< idx_t >&             enumerator,
                            uint_t                           level )
{
   idx_t enumIndex = 0;
   enumerateVertexDoFFunctionConsistently( enumeratorMap, enumBounds, enumerator, level, enumIndex );
}

void enumerateConsistently( std::pair< std::map< uint_t, PrimitiveID >, enumeratorBounds >& enumPair,
                            P1Function< idx_t >&                                            enumerator,
                            uint_t                                                          level )
{
   enumerateConsistently( enumPair.first, enumPair.second, enumerator, level );
}

void enumerateConsistently( P1Function< idx_t >& enumerator, uint_t level )
{
   auto enumPair = generateEnumeratorMap( enumerator.getStorage() );
   enumerateConsistently( enumPair.first, enumPair.second, enumerator, level );
}

void enumerateConsistently( std::map< uint_t, PrimitiveID >& enumeratorMap,
                            enumeratorBounds&                enumBounds,
                            EdgeDoFFunction< idx_t >&        enumerator,
                            uint_t                           level )
{
   idx_t enumIndex = 0;
   enumerateEdgeDoFFunctionConsistently( enumeratorMap, enumBounds, enumerator, level, enumIndex );
}

void enumerateConsistently( std::pair< std::map< uint_t, PrimitiveID >, enumeratorBounds >& enumPair,
                            EdgeDoFFunction< idx_t >&                                       enumerator,
                            uint_t                                                          level )
{
   enumerateConsistently( enumPair.first, enumPair.second, enumerator, level );
}

void enumerateConsistently( EdgeDoFFunction< idx_t >& enumerator, uint_t level )
{
   auto enumPair = generateEnumeratorMap( enumerator.getStorage() );
   enumerateConsistently( enumPair.first, enumPair.second, enumerator, level );
}

void enumerateConsistently( std::map< uint_t, PrimitiveID >& enumeratorMap,
                            enumeratorBounds&                enumBounds,
                            P1VectorFunction< idx_t >&       enumerator,
                            uint_t                           level )
{
   idx_t enumIndex = 0;
   for ( uint_t i = 0; i < enumerator.getDimension(); i++ )
   {
      enumerateVertexDoFFunctionConsistently( enumeratorMap, enumBounds, enumerator[i], level, enumIndex );
   }
}

void enumerateConsistently( std::pair< std::map< uint_t, PrimitiveID >, enumeratorBounds >& enumPair,
                            P1VectorFunction< idx_t >&                                      enumerator,
                            uint_t                                                          level )
{
   enumerateConsistently( enumPair.first, enumPair.second, enumerator, level );
}

void enumerateConsistently( P1VectorFunction< idx_t >& enumerator, uint_t level )
{
   auto enumPair = generateEnumeratorMap( enumerator.getStorage() );
   enumerateConsistently( enumPair.first, enumPair.second, enumerator, level );
}

void enumerateConsistently( std::map< uint_t, PrimitiveID >& enumeratorMap,
                            enumeratorBounds&                enumBounds,
                            P2VectorFunction< idx_t >&       enumerator,
                            uint_t                           level )
{
   idx_t enumIndex = 0;
   for ( uint_t i = 0; i < enumerator.getDimension(); i++ )
   {
      enumerateVertexDoFFunctionConsistently( enumeratorMap, enumBounds, enumerator[i].getVertexDoFFunction(), level, enumIndex );
      enumerateEdgeDoFFunctionConsistently( enumeratorMap, enumBounds, enumerator[i].getEdgeDoFFunction(), level, enumIndex );
   }
}

void enumerateConsistently( std::pair< std::map< uint_t, PrimitiveID >, enumeratorBounds >& enumPair,
                            P2VectorFunction< idx_t >&                                      enumerator,
                            uint_t                                                          level )
{
   enumerateConsistently( enumPair.first, enumPair.second, enumerator, level );
}

void enumerateConsistently( P2VectorFunction< idx_t >& enumerator, uint_t level )
{
   auto enumPair = generateEnumeratorMap( enumerator.getStorage() );
   enumerateConsistently( enumPair.first, enumPair.second, enumerator, level );
}

void enumerateConsistently( std::map< uint_t, PrimitiveID >& enumeratorMap,
                            enumeratorBounds&                enumBounds,
                            P2P1TaylorHoodFunction< idx_t >& enumerator,
                            uint_t                           level )
{
   idx_t enumIndex = 0;
   for ( uint_t i = 0; i < enumerator.uvw().getDimension(); i++ )
   {
      enumerateVertexDoFFunctionConsistently(
          enumeratorMap, enumBounds, enumerator.uvw()[i].getVertexDoFFunction(), level, enumIndex );
      enumerateEdgeDoFFunctionConsistently(
          enumeratorMap, enumBounds, enumerator.uvw()[i].getEdgeDoFFunction(), level, enumIndex );
   }
   enumerateVertexDoFFunctionConsistently( enumeratorMap, enumBounds, enumerator.p(), level, enumIndex );
}

void enumerateConsistently( std::pair< std::map< uint_t, PrimitiveID >, enumeratorBounds >& enumPair,
                            P2P1TaylorHoodFunction< idx_t >&                                enumerator,
                            uint_t                                                          level )
{
   enumerateConsistently( enumPair.first, enumPair.second, enumerator, level );
}

void enumerateConsistently( P2P1TaylorHoodFunction< idx_t >& enumerator, uint_t level )
{
   auto enumPair = generateEnumeratorMap( enumerator.getStorage() );
   enumerateConsistently( enumPair.first, enumPair.second, enumerator, level );
}

} // namespace hyteg
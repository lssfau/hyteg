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
#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::uint_t;

namespace hyteg {

static void testP1DataHandling()
{
   uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   uint_t minLevel = 1;
   uint_t maxLevel = 10;

   std::string meshFileName = prependHyTeGMeshDir( "tri_2el.msh" );

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
   PrimitiveStorage      storage( setupStorage );

   WALBERLA_LOG_INFO( setupStorage );

   auto vertexP1FunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< real_t >, Vertex > >(
       minLevel, maxLevel, vertexDoFMacroVertexFunctionMemorySize );
   auto edgeP1FunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< real_t >, Edge > >(
       minLevel, maxLevel, vertexDoFMacroEdgeFunctionMemorySize );
   auto faceP1FunctionMemoryDataHandling = std::make_shared< MemoryDataHandling< FunctionMemory< real_t >, Face > >(
       minLevel, maxLevel, vertexDoFMacroFaceFunctionMemorySize );

   PrimitiveDataID< FunctionMemory< real_t >, Vertex > vertexP1FunctionMemoryID;
   storage.addVertexData( vertexP1FunctionMemoryID, vertexP1FunctionMemoryDataHandling, "P1 vertex data" );

   PrimitiveDataID< FunctionMemory< real_t >, Edge > edgeP1FunctionMemoryID;
   storage.addEdgeData( edgeP1FunctionMemoryID, edgeP1FunctionMemoryDataHandling, "P1 edge data" );

   PrimitiveDataID< FunctionMemory< real_t >, Face > faceP1FunctionMemoryID;
   storage.addFaceData( faceP1FunctionMemoryID, faceP1FunctionMemoryDataHandling, "P1 face data" );

   for ( const auto& it : storage.getVertices() )
   {
      PrimitiveID vertexID = it.first;
      auto        vertex   = it.second;

      FunctionMemory< real_t >* mem = vertex->getData( vertexP1FunctionMemoryID );

      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         WALBERLA_LOG_INFO( "Vertex: " << vertexID << " | data size (level " << level << ") : " << mem->getSize( level ) );
      }
   }

   for ( const auto& it : storage.getEdges() )
   {
      PrimitiveID edgeID = it.first;
      auto        edge   = it.second;

      FunctionMemory< real_t >* mem = edge->getData( edgeP1FunctionMemoryID );

      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         WALBERLA_LOG_INFO( "Edge:   " << edgeID << " | data size (level " << level << ") : " << mem->getSize( level ) );
      }
   }

   for ( const auto& it : storage.getFaces() )
   {
      PrimitiveID faceID = it.first;
      auto        face   = it.second;

      FunctionMemory< real_t >* mem = face->getData( faceP1FunctionMemoryID );

      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         WALBERLA_LOG_INFO( "Face:   " << faceID << " | data size (level " << level << ") : " << mem->getSize( level ) );
      }
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::testP1DataHandling();

   return EXIT_SUCCESS;
}

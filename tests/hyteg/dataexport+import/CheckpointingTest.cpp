/*
 * Copyright (c) 2023 Marcus Mohr.
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
#include "core/timing/all.h"

#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointExporter.hpp"
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointHelpers.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using namespace hyteg;

void checkTagGenerationForP2( const P2Function< real_t >& func, uint_t level )
{
   const auto& storage = func.getStorage();

   WALBERLA_LOG_INFO_ON_ROOT( "***********" );
   WALBERLA_LOG_INFO_ON_ROOT( "*  CELLS  *" );
   WALBERLA_LOG_INFO_ON_ROOT( "***********" );

   PrimitiveDataID< FunctionMemory< real_t >, Cell > cellDataIndex = func.getVertexDoFFunction().getCellDataID();
   for ( const auto& macroIterator : storage->getCells() )
   {
      Cell&       cell = *macroIterator.second;
      std::string tag  = adiosCheckpointHelpers::generateVariableName( func.getFunctionName(), cell.getID(), level );
      uint_t      size = cell.getData( cellDataIndex )->getSize( level );
      WALBERLA_LOG_INFO( "Generated tag = '" << tag << "'" );
      WALBERLA_LOG_INFO( "Data buffer size = " << size );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "***********" );
   WALBERLA_LOG_INFO_ON_ROOT( "*  FACES  *" );
   WALBERLA_LOG_INFO_ON_ROOT( "***********" );

   PrimitiveDataID< FunctionMemory< real_t >, Face > faceDataIndex = func.getVertexDoFFunction().getFaceDataID();
   for ( const auto& macroIterator : storage->getFaces() )
   {
      Face&       face = *macroIterator.second;
      std::string tag  = adiosCheckpointHelpers::generateVariableName( func.getFunctionName(), face.getID(), level );
      uint_t      size = face.getData( faceDataIndex )->getSize( level );
      WALBERLA_LOG_INFO( "Generated tag = '" << tag << "'" );
      WALBERLA_LOG_INFO( "Data buffer size = " << size );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "***********" );
   WALBERLA_LOG_INFO_ON_ROOT( "*  EDGES  *" );
   WALBERLA_LOG_INFO_ON_ROOT( "***********" );

   PrimitiveDataID< FunctionMemory< real_t >, Edge > edgeDataIndex = func.getVertexDoFFunction().getEdgeDataID();
   for ( const auto& macroIterator : storage->getEdges() )
   {
      Edge&       edge = *macroIterator.second;
      std::string tag  = adiosCheckpointHelpers::generateVariableName( func.getFunctionName(), edge.getID(), level );
      uint_t      size = edge.getData( edgeDataIndex )->getSize( level );
      WALBERLA_LOG_INFO( "Generated tag = '" << tag << "'" );
      WALBERLA_LOG_INFO( "Data buffer size = " << size );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "***********" );
   WALBERLA_LOG_INFO_ON_ROOT( "*  VERTS  *" );
   WALBERLA_LOG_INFO_ON_ROOT( "***********" );

   PrimitiveDataID< FunctionMemory< real_t >, Vertex > vertexVertexDoFDataIndex = func.getVertexDoFFunction().getVertexDataID();
   PrimitiveDataID< FunctionMemory< real_t >, Vertex > vertexEdgeDoFDataIndex   = func.getEdgeDoFFunction().getVertexDataID();
   for ( const auto& macroIterator : storage->getVertices() )
   {
      Vertex&     vertex = *macroIterator.second;
      std::string tag    = adiosCheckpointHelpers::generateVariableName( func.getFunctionName(), vertex.getID(), level );
      uint_t      size   = vertex.getData( vertexVertexDoFDataIndex )->getSize( level ) +
                    vertex.getData( vertexEdgeDoFDataIndex )->getSize( level );
      WALBERLA_LOG_INFO( "Generated tag = '" << tag << "'" );
      WALBERLA_LOG_INFO( "Data buffer size = " << size );
   }
}

void storeCheckpoint( std::string filePath, std::string fileName )
{
   std::string           meshFileName{ "../../data/meshes/3D/pyramid_2el.msh" };
   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const uint_t minLevel = 2;
   const uint_t maxLevel = 3;

   P1Function< real_t >  funcP1( "P1_Test_Function", storage, minLevel, maxLevel );
   P2Function< real_t >  funcP2( "P2_Test_Function", storage, minLevel, maxLevel );
   P2Function< int64_t > intFunc( "P2_Test_Function<int64>", storage, minLevel, maxLevel );

   P1VectorFunction< real_t > funcP1Vec( "P1_Test_Vector_Function", storage, minLevel, maxLevel );
   P2VectorFunction< real_t > funcP2Vec( "P2_Test_Vector_Function", storage, minLevel, maxLevel );

   P2P1TaylorHoodFunction< real_t > stokesFunc( "Stokes Function", storage, minLevel, maxLevel );

   for ( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
   {
      funcP1Vec.interpolate( { real_c( 1 ), real_c( 2 ), real_c( 3 ) }, lvl );
   }

   AdiosCheckpointExporter checkpointer( "" );
   checkpointer.registerFunction( funcP1, minLevel, maxLevel );
   checkpointer.registerFunction( funcP2, minLevel, maxLevel );
   checkpointer.registerFunction( intFunc, maxLevel, maxLevel );
   checkpointer.registerFunction( funcP1Vec, minLevel, maxLevel );
   checkpointer.registerFunction( funcP2Vec, minLevel, maxLevel );
   checkpointer.registerFunction( stokesFunc, maxLevel, maxLevel );

   checkpointer.storeCheckpoint( filePath, fileName, { "MeshFile" }, { meshFileName } );
}

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   storeCheckpoint( ".", "CheckpointingTest.bp" );
}

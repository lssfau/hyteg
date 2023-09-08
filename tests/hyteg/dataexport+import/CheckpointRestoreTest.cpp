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
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointImporter.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using namespace hyteg;

void exportCheckpoint( const std::string& filePath,
                       const std::string& fileName,
                       const std::string& meshFileName,
                       const uint_t       minLevel,
                       const uint_t       maxLevel )
{
   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   P1Function< real_t > funcP1( "P1_Test_Function", storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > expr = []( const Point3D& p ) -> real_t {
      return real_c( -2 ) * p[0] + real_c( 3 ) * p[1];
   };

   for ( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
   {
      funcP1.interpolate( expr, lvl );
   }

   AdiosCheckpointExporter checkpointer( "" );
   checkpointer.registerFunction( funcP1, minLevel, maxLevel );
   checkpointer.storeCheckpoint( filePath, fileName, { "MeshFile" }, { meshFileName } );
}

void importCheckpoint( const std::string& filePath,
                       const std::string& fileName,
                       const std::string& meshFileName,
                       const uint_t       minLevel,
                       const uint_t       maxLevel,
                       bool               verbose = true )
{
   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   AdiosCheckpointImporter restorer( filePath, fileName, "" );

   if ( verbose )
   {
      restorer.printCheckpointInfo();
   }

   auto&                funcDescr = restorer.getFunctionDetails();
   P1Function< real_t > funcP1( funcDescr[0].name, storage, funcDescr[0].minLevel, funcDescr[0].maxLevel );
   restorer.restoreFunction( funcP1 );

   AdiosWriter adiosWriter( ".", "CheckpointRestoreTestOutput", storage );
   adiosWriter.add( funcP1 );
   adiosWriter.write( maxLevel );
}

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // =====================
   //  Set Test Parameters
   // =====================
   const uint_t minLevel = 2;
   const uint_t maxLevel = 3;

   std::string filePath{ "." };
   std::string fileName{ "CheckpointRestoreTest.bp" };
   std::string meshFile{ "../../data/meshes/LShape_6el.msh" };

   // ===================
   //  Create Checkpoint
   // ===================
   exportCheckpoint( filePath, fileName, meshFile, minLevel, maxLevel );

   // ===================
   //  Import Checkpoint
   // ===================
   importCheckpoint( filePath, fileName, meshFile, minLevel, maxLevel );
}

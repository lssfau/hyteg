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
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

// This program tests export and import of FE functions using ADIOS2
//
// There are two modes:
//
// (1) Just run the executable sequentially or in parallel, in this
//     case FE functions are created, exported and afterwards read
//     in again. The original and restored functions are compared
//     against each other.
//
// (2) Run with cli options "CheckpointRestoreTest.prm -Parameters.onlyImport=true"
//     In this case the program tries to import a function from a previous run of
//     of type (1). We can do this to see that we are not bound to have the same
//     number of MPI processes, or primitive distribution for the export and
//     import phases.

using namespace hyteg;

const std::string markerExport{ "AdiosCheckpointExport" };
const std::string markerImport{ "AdiosCheckpointImport" };

// -------------------------------------------
//  Generate PrimitiveStorage from a meshfile
// -------------------------------------------
std::shared_ptr< PrimitiveStorage > generateStorage( const std::string& meshFileName )
{
   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   loadbalancing::roundRobin( setupStorage );
   return std::make_shared< PrimitiveStorage >( setupStorage );
}

// -------------------------------------------------------
//  Generate PrimitiveStorage for a thick spherical shell
// -------------------------------------------------------
std::shared_ptr< PrimitiveStorage > generateStorage( uint_t nTan, uint_t nRad )
{
   const real_t          rMin = real_c( 1 );
   const real_t          rMax = real_c( 2 );
   MeshInfo              mesh( MeshInfo::meshSphericalShell( nTan, nRad, rMin, rMax ) );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   IcosahedralShellMap::setMap( setupStorage );
   loadbalancing::roundRobin( setupStorage );
   return std::make_shared< PrimitiveStorage >( setupStorage );
}

// --------------
//  runBenchmark
// --------------
void runBenchmark( const std::string&                  filePath,
                   const std::string&                  fileName,
                   std::shared_ptr< PrimitiveStorage > storage,
                   const uint_t                        minLevel,
                   const uint_t                        maxLevel )
{
   // =============
   //  Export Part
   // =============
   {
      // Setup FE functions
      WALBERLA_LOG_INFO_ON_ROOT( " -> generating some FE functions" );
      P2P1TaylorHoodFunction< real_t > funcStokes( "Stokes", storage, minLevel, maxLevel );
      P2Function< real_t >             funcTemperature( "Temperature", storage, minLevel, maxLevel );
      P1Function< real_t >             funcComposition( "Composition", storage, minLevel, maxLevel );

      //  Create Checkpoint
      WALBERLA_LOG_INFO_ON_ROOT( " -> exporting checkpoint" );
      AdiosCheckpointExporter checkpointer( "" );
      checkpointer.registerFunction( funcStokes, minLevel, maxLevel );
      checkpointer.registerFunction( funcTemperature, minLevel, maxLevel );
      checkpointer.registerFunction( funcComposition, minLevel, maxLevel );

      storage->getTimingTree()->start( markerExport );

      checkpointer.storeCheckpoint( filePath, fileName );

      storage->getTimingTree()->stop( markerExport );
   }

   // =============
   //  Import Part
   // =============
   {
      WALBERLA_LOG_INFO_ON_ROOT( " -> importing from checkpoint" );
      P2P1TaylorHoodFunction< real_t > funcStokes( "Stokes", storage, minLevel, maxLevel );
      P2Function< real_t >             funcTemperature( "Temperature", storage, minLevel, maxLevel );
      P1Function< real_t >             funcComposition( "Composition", storage, minLevel, maxLevel );

      AdiosCheckpointImporter restorer( filePath, fileName, "" );
      restorer.printCheckpointInfo();

      storage->getTimingTree()->start( markerImport );

      restorer.restoreFunction( funcStokes.uvw() );
      restorer.restoreFunction( funcStokes.p() );
      restorer.restoreFunction( funcTemperature );
      restorer.restoreFunction( funcComposition );

      storage->getTimingTree()->stop( markerImport );
   }
}

std::string separator{ "------------------------------------------------" };

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // Be talkative
   WALBERLA_LOG_INFO_ON_ROOT( separator );
   WALBERLA_LOG_INFO_ON_ROOT( "--> BENCHMARK CHECKPOINT-RESTORE WITH ADIOS2 <--" );
   WALBERLA_LOG_INFO_ON_ROOT( separator );
   buildinfo::printBuildInfo();
   buildinfo::printGitInfo();
   WALBERLA_LOG_INFO_ON_ROOT( separator );

   // =====================
   //  Set Test Parameters
   // =====================

   // If the name of a config-file was given on the command-line
   // walberlaEnv will have opened it already
   walberla::shared_ptr< walberla::config::Config > cfg;
   if ( argc == 1 )
   {
      std::string defaultParameterFile{ "./BenchmarkCheckpointRestore.prm" };
      WALBERLA_LOG_INFO_ON_ROOT( "No parameter file given on command-line! Using'" << defaultParameterFile << "'" );
      cfg = std::make_shared< walberla::config::Config >();
      cfg->readParameterFile( defaultParameterFile.c_str() );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Reading parameters from '" << argv[1] << "'" );
      cfg = walberlaEnv.config();
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );
   WALBERLA_LOG_INFO_ON_ROOT( parameters );
   WALBERLA_LOG_INFO_ON_ROOT( separator );

   // =================
   //  Mesh Generation
   // =================
   std::string                         meshingType = parameters.getParameter< std::string >( "meshType" );
   std::shared_ptr< PrimitiveStorage > storage     = nullptr;

   if ( meshingType == "SphericalShell" )
   {
      uint_t nTan = parameters.getParameter< uint_t >( "nTan" );
      uint_t nRad = parameters.getParameter< uint_t >( "nRad" );
      storage     = generateStorage( nTan, nRad );
   }
   else if ( meshingType == "fromFile" )
   {
      std::string meshFileName = parameters.getParameter< std::string >( "meshFileName" );
      storage                  = generateStorage( meshFileName );
   }
   else
   {
      WALBERLA_ABORT( "Don't understand meshType = '" << meshingType << "'" );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "Successfully generated mesh" );

   // ===========
   //  Run Tests
   // ===========
   std::string  filePath = parameters.getParameter< std::string >( "filePath" );
   std::string  fileName = parameters.getParameter< std::string >( "fileName" );
   const uint_t minLevel = parameters.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel = parameters.getParameter< uint_t >( "maxLevel" );

   WALBERLA_LOG_INFO_ON_ROOT( "Performing benchmark measurements:" );

   const uint_t numTests = parameters.getParameter< uint_t >( "numTests" );

   for ( uint_t count = 1; count <= numTests; ++count )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "test run " << count << "/" << numTests );
      runBenchmark( filePath, fileName, storage, minLevel, maxLevel );
      WALBERLA_LOG_INFO_ON_ROOT( separator );
   }

   auto timerExport = storage->getTimingTree()->operator[]( markerExport );
   auto timerImport = storage->getTimingTree()->operator[]( markerImport );
   WALBERLA_LOG_INFO_ON_ROOT( "Average time for checkpoint export: " << timerExport.average() );
   WALBERLA_LOG_INFO_ON_ROOT( "Average time for checkpoint import: " << timerImport.average() );
}

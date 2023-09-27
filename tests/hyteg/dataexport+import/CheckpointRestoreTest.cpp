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

std::shared_ptr< PrimitiveStorage > generateStorage( const std::string& meshFileName )
{
   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   return std::make_shared< PrimitiveStorage >( setupStorage );
}

template < template < typename > class func_t, typename value_t >
auto exportCheckpoint( const std::string&                         filePath,
                       const std::string&                         fileName,
                       const std::shared_ptr< PrimitiveStorage >& storage,
                       const uint_t                               minLevel,
                       const uint_t                               maxLevel,
                       bool                                       verbose = false )
{
   WALBERLA_UNUSED( verbose );

   std::string funcName = "Test_" + FunctionTrait< func_t< value_t > >::getTypeName();

   func_t< value_t > feFunc( funcName, storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > expr = []( const Point3D& p ) -> real_t {
      return real_c( -2 ) * p[0] + real_c( 3 ) * p[1];
   };

   for ( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
   {
      feFunc.interpolate( expr, lvl );
   }

   AdiosCheckpointExporter checkpointer( "" );
   checkpointer.registerFunction( feFunc, minLevel, maxLevel );
   checkpointer.storeCheckpoint( filePath, fileName );

   return feFunc;
}

template < template < typename > class func_t, typename value_t >
auto importCheckpoint( const std::string&                         filePath,
                       const std::string&                         fileName,
                       const std::shared_ptr< PrimitiveStorage >& storage,
                       const uint_t                               minLevel,
                       const uint_t                               maxLevel,
                       bool                                       verbose = false )
{
   WALBERLA_UNUSED( verbose );

   AdiosCheckpointImporter restorer( filePath, fileName, "" );

   if ( verbose )
   {
      restorer.printCheckpointInfo();
   }

   auto& funcDescr = restorer.getFunctionDetails();
   WALBERLA_CHECK_EQUAL( funcDescr[0].minLevel, minLevel );
   WALBERLA_CHECK_EQUAL( funcDescr[0].maxLevel, maxLevel );
   func_t< value_t > feFunc( funcDescr[0].name, storage, funcDescr[0].minLevel, funcDescr[0].maxLevel );
   restorer.restoreFunction( feFunc );

   return feFunc;
}

template < template < typename > class func_t, typename value_t >
value_t computeCheckValue( const func_t< value_t >& feFunc, uint_t level )
{
   value_t error = static_cast< value_t >( 0 );

   if constexpr ( std::is_same_v< func_t< value_t >, P2P1TaylorHoodFunction< value_t > > )
   {
      error = feFunc.p().getMaxMagnitude( level, All );
      for ( uint_t k = 0; k < feFunc.getDimension(); ++k )
      {
         error += feFunc.uvw()[k].getMaxMagnitude( level, All );
      }
   }
   else
   {
      for ( uint_t k = 0; k < feFunc.getDimension(); ++k )
      {
         error += feFunc[k].getMaxMagnitude( level, All );
      }
   }
   return error;
}

template < template < typename > class func_t, typename value_t >
void runTestWithIdenticalCommunicator( const std::string& filePath,
                                       const std::string& fileName,
                                       const std::string& meshFileName,
                                       const uint_t       minLevel,
                                       const uint_t       maxLevel,
                                       bool               verbose = false )
{
   bool exportFuncs = false;

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "==============================================================" );
      WALBERLA_LOG_INFO_ON_ROOT( "Testing with the following parameters:" );
      WALBERLA_LOG_INFO_ON_ROOT( " - function of type ... " << FunctionTrait< func_t< value_t > >::getTypeName() );
      WALBERLA_LOG_INFO_ON_ROOT( " - value type ......... " << adiosCheckpointHelpers::valueTypeToString< value_t >() );
      WALBERLA_LOG_INFO_ON_ROOT( " - meshfile ........... '" << meshFileName << "'" );
      WALBERLA_LOG_INFO_ON_ROOT( " - filePath ........... '" << filePath << "'" );
      WALBERLA_LOG_INFO_ON_ROOT( " - fileName ........... '" << fileName << "'" );
      WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------------------" );
   }

   auto storage = generateStorage( meshFileName );

   //  Create Checkpoint
   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * exporting checkpoint" );
   }
   func_t< value_t > funcOriginal = exportCheckpoint< func_t, value_t >( filePath, fileName, storage, minLevel, maxLevel );

   //  Import Checkpoint
   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * importing checkpoint" );
   }
   func_t< value_t > funcRestored = importCheckpoint< func_t, value_t >( filePath, fileName, storage, minLevel, maxLevel );

   func_t< value_t > difference( "Difference", storage, minLevel, maxLevel );
   for ( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
   {
      difference.assign(
          { static_cast< value_t >( 1 ), static_cast< value_t >( -1 ) }, { funcOriginal, funcRestored }, lvl, All );
      value_t error = computeCheckValue( difference, lvl );

      if ( verbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( " * checking differences on refinement level " << lvl << " ..." );
         WALBERLA_CHECK_LESS_EQUAL( error, static_cast< value_t >( 0 ) );
         WALBERLA_LOG_INFO_ON_ROOT( "   ... okay" );
      }
      if ( exportFuncs )
      {
         AdiosWriter adiosWriter( ".", "CheckpointRestoreTestOutput", storage );
         adiosWriter.add( funcOriginal );
         adiosWriter.add( funcRestored );
         adiosWriter.add( difference );
         adiosWriter.write( lvl );
      }
   }
}

template < template < typename > class func_t, typename value_t >
void runTestWithOtherCommunicator( const std::string& filePath,
                                   const std::string& fileName,
                                   const std::string& meshFileName,
                                   const uint_t       minLevel,
                                   const uint_t       maxLevel,
                                   bool               verbose = false )
{
   bool exportFuncs = false;

   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "==============================================================" );
      WALBERLA_LOG_INFO_ON_ROOT( "Testing with the following parameters:" );
      WALBERLA_LOG_INFO_ON_ROOT( " - function of type ... " << FunctionTrait< func_t< value_t > >::getTypeName() );
      WALBERLA_LOG_INFO_ON_ROOT( " - value type ......... " << adiosCheckpointHelpers::valueTypeToString< value_t >() );
      WALBERLA_LOG_INFO_ON_ROOT( " - meshfile ........... '" << meshFileName << "'" );
      WALBERLA_LOG_INFO_ON_ROOT( " - filePath ........... '" << filePath << "'" );
      WALBERLA_LOG_INFO_ON_ROOT( " - fileName ........... '" << fileName << "'" );
      WALBERLA_LOG_INFO_ON_ROOT( "--------------------------------------------------------------" );
   }

   auto storage = generateStorage( meshFileName );

   //  Import Checkpoint
   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * importing checkpoint" );
   }
   func_t< value_t > funcRestored = importCheckpoint< func_t, value_t >( filePath, fileName, storage, minLevel, maxLevel );

   for ( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
   {
      if ( exportFuncs )
      {
         AdiosWriter adiosWriter( ".", "CheckpointRestoreTestOutput", storage );
         adiosWriter.add( funcRestored );
         adiosWriter.write( lvl );
      }

      value_t checkValue = computeCheckValue( funcRestored, maxLevel );
      WALBERLA_CHECK_LESS_EQUAL( std::abs( checkValue - static_cast< value_t >( 3 ) ), static_cast< value_t >( 0 ) );
   }
}

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   if ( argc > 3 )
   {
      WALBERLA_ABORT( "Wrong number of command-line arguments!" );
   }
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // =====================
   //  Set Test Parameters
   // =====================
   const uint_t minLevel = 0;
   const uint_t maxLevel = 3;

   std::string filePath{ "." };
   std::string fileName{ "CheckpointRestoreTest.bp" };
   // std::string meshFile{ "../../data/meshes/LShape_6el.msh" };
   std::string meshFile{ "../../data/meshes/3D/cube_6el.msh" };

   bool onlyImport = false;
   if ( argc > 1 )
   {
      walberla::Config::BlockHandle parameters = walberlaEnv.config()->getOneBlock( "Parameters" );
      onlyImport                               = parameters.getParameter< bool >( "onlyImport" );
   }

   // ===========
   //  Run Tests
   // ===========
   if ( !onlyImport )
   {
      runTestWithIdenticalCommunicator< P1Function, real_t >( filePath, fileName, meshFile, minLevel, maxLevel, true );
      runTestWithIdenticalCommunicator< P1Function, int64_t >( filePath, fileName, meshFile, minLevel, maxLevel, true );

      runTestWithIdenticalCommunicator< P1VectorFunction, real_t >( filePath, fileName, meshFile, minLevel, maxLevel, true );
      runTestWithIdenticalCommunicator< P2Function, int32_t >( filePath, fileName, meshFile, minLevel, maxLevel, true );

      // we are going to reuse this checkpoint in the next part of this pipeline job
      fileName = "CheckpointRestoreTest+reuse.bp";
      runTestWithIdenticalCommunicator< P2Function, real_t >( filePath, fileName, meshFile, minLevel, maxLevel, true );

      // We currently would need to import the two component functions separately; Better than having specialised code here,
      // alter AdiosCheckpoint[Ex|Im]porter
      // runTestWithIdenticalCommunicator< P2P1TaylorHoodFunction, real_t >( filePath, fileName, meshFile, minLevel, maxLevel, true );
   }

   // =====================
   //  Run Other Test Form
   // =====================
   else
   {
      fileName = "CheckpointRestoreTest+reuse.bp";
      runTestWithOtherCommunicator< P2Function, real_t >( filePath, fileName, meshFile, minLevel, maxLevel, true );
      return EXIT_SUCCESS;
   }
}

// ensure speficied interfaces exist by making compiler explicitely instantiate the CRTP "base" class
template class hyteg::CheckpointExporter< hyteg::AdiosCheckpointExporter >;
template class hyteg::CheckpointImporter< hyteg::AdiosCheckpointImporter >;

/*
 * Copyright (c) 2023-2024 Marcus Mohr.
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

#include <utility>

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
// (2) Run with cli options
//     "CheckpointRestoreTest.prm -Parameters.onlyImport=true -Parameters.restoreFromFileWithName=<name of BP file>"
//     In this case the program tries to import a function from a previous run of
//     type (1). We can do this to see that we are not bound to have the same
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
auto exportCheckpoint( const std::string&                                        filePath,
                       const std::string&                                        fileName,
                       const std::shared_ptr< PrimitiveStorage >&                storage,
                       const uint_t                                              minLevel,
                       const uint_t                                              maxLevel,
                       const std::map< std::string, adiosHelpers::adiostype_t >& userAttributes = {},
                       bool                                                      verbose        = false )
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
   checkpointer.storeCheckpoint( filePath, fileName, userAttributes );

   return feFunc;
}

template < template < typename > class func_t, typename value_t >
auto exportCheckpointContinuous( const std::string&                         filePath,
                                 const std::string&                         fileName,
                                 const std::shared_ptr< PrimitiveStorage >& storage,
                                 const uint_t                               minLevel,
                                 const uint_t                               maxLevel,
                                 const uint_t                               nSteps,
                                 bool                                       verbose = false )
{
   WALBERLA_UNUSED( verbose );

   std::vector< std::shared_ptr< func_t< value_t > > > funcsToWrite;
   std::vector< real_t >                               timeVector;

   funcsToWrite.reserve( nSteps );
   timeVector.reserve( nSteps );

   std::string       funcName = "Test_" + FunctionTrait< func_t< value_t > >::getTypeName();
   func_t< value_t > feFunc( funcName, storage, minLevel, maxLevel );

   AdiosCheckpointExporter checkpointer( "" );

   checkpointer.registerFunction( feFunc, minLevel, maxLevel );

   for ( uint_t tStep = 1U; tStep <= nSteps; tStep++ )
   {
      std::function< real_t( const hyteg::Point3D& ) > expr = [tStep]( const Point3D& p ) -> real_t {
         return real_c( -2 ) * p[0] + real_c( 3 ) * p[1] + real_c( 4 ) * real_c( tStep );
      };

      for ( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
      {
         feFunc.interpolate( expr, lvl );
      }

      real_t time = real_c( tStep ) * real_c( 0.1 );

      timeVector.push_back( time );

      auto func = std::make_shared< func_t< value_t > >( funcName + "_" + std::to_string( tStep ), storage, minLevel, maxLevel );

      for ( uint_t l = minLevel; l <= maxLevel; l++ )
      {
         func->assign( { static_cast< value_t >( 1.0 ) }, { feFunc }, l, All );
      }

      funcsToWrite.push_back( func );

      checkpointer.storeCheckpointContinuous( filePath, fileName, time, tStep == nSteps );
   }

   return std::make_pair( funcsToWrite, timeVector );
}

template < template < typename > class func_t, typename value_t >
auto importCheckpoint( const std::string&                                  filePath,
                       const std::string&                                  fileName,
                       const std::shared_ptr< PrimitiveStorage >&          storage,
                       const uint_t                                        minLevel,
                       const uint_t                                        maxLevel,
                       std::map< std::string, adiosHelpers::adiostype_t >& userAttributes,
                       bool                                                verbose = false )
{
   WALBERLA_UNUSED( verbose );

   AdiosCheckpointImporter restorer( filePath, fileName, "" );

   // if ( verbose )
   // {
   restorer.printCheckpointInfo();
   // }

   restorer.readAllUserAttributes( userAttributes );

   auto& funcDescr = restorer.getFunctionDetails();
   WALBERLA_CHECK_EQUAL( funcDescr[0].minLevel, minLevel );
   WALBERLA_CHECK_EQUAL( funcDescr[0].maxLevel, maxLevel );
   func_t< value_t > feFunc( funcDescr[0].name, storage, funcDescr[0].minLevel, funcDescr[0].maxLevel );
   restorer.restoreFunction( feFunc );

   return feFunc;
}

template < template < typename > class func_t, typename value_t >
auto importCheckpointContinuous( const std::string&                         filePath,
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

   std::vector< real_t > timestepInfo = restorer.getTimestepInfo();
   timestepInfo                       = restorer.getTimestepInfo();
   timestepInfo                       = restorer.getTimestepInfo();

   uint_t checkpointSteps = timestepInfo.size();

   std::vector< std::shared_ptr< func_t< value_t > > > funcsToWrite;

   funcsToWrite.resize( checkpointSteps );

   // Let's read in reverse, as we need for adjoints!
   for ( uint_t tStep = checkpointSteps; tStep >= 1; tStep-- )
   {
      restorer.restoreFunction( feFunc, tStep - 1 );

      auto func = std::make_shared< func_t< value_t > >(
          funcDescr[0].name + "_" + std::to_string( tStep - 1 ), storage, minLevel, maxLevel );

      for ( uint_t l = minLevel; l <= maxLevel; l++ )
      {
         func->assign( { static_cast< value_t >( 1.0 ) }, { feFunc }, l, All );
      }

      funcsToWrite[tStep - 1] = func;
   }

   return std::make_pair( funcsToWrite, timestepInfo );
}

template < template < typename > class func_t, typename value_t >
value_t computeCheckValue( const func_t< value_t >& feFunc, uint_t level )
{
   value_t error = static_cast< value_t >( 0 );

   if constexpr ( std::is_same_v< func_t< value_t >, P2P1TaylorHoodFunction< value_t > > )
   {
      error = feFunc.p().getMaxDoFMagnitude( level, All );
      for ( uint_t k = 0; k < feFunc.getDimension(); ++k )
      {
         error += feFunc.uvw()[k].getMaxDoFMagnitude( level, All );
      }
   }
   else
   {
      for ( uint_t k = 0; k < feFunc.getDimension(); ++k )
      {
         error += feFunc[k].getMaxDoFMagnitude( level, All );
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
                                       bool               verbose          = false,
                                       const uint_t       testContinuous   = false,
                                       const uint_t       nStepsContinuous = 0U )
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

   std::map< std::string, adiosHelpers::adiostype_t > userAttributes         = {};
   std::map< std::string, adiosHelpers::adiostype_t > userAttributesToImport = {};

   // test all variant types
   userAttributes["Test0"] = int( -4 );
   userAttributes["Test1"] = long( -12l );
   userAttributes["Test2"] = float( 0.9f );
   userAttributes["Test3"] = double( -1.25 );
   userAttributes["Test4"] = uint_t( 99ul );
   userAttributes["Test5"] = bool( true );
   userAttributes["Test6"] = std::string( "Test string" );
   userAttributes["Test7"]  = std::vector< int >( { 1, -2, 3, -4, 5 } );
   userAttributes["Test8"]  = std::vector< long >( { 5l, -4l, 3l, -2l, 1l } );
   userAttributes["Test9"]  = std::vector< float >( { -0.1f, 0.2f, -0.3f, 0.4f, -0.5f } );
   userAttributes["Test10"] = std::vector< double >( { -0.2, 0.4, -0.6, 0.8, -1.0 } );
   userAttributes["Test11"] = std::vector< uint_t >( { 0ul, 1ul, 10ul, 100ul, 1000ul, 100000ul } );
   userAttributes["Test12"] = std::vector< bool >( { false, false, false, true, false, true } );
   userAttributes["Test13"] = std::vector< std::string >( { "This", "is", "a", "test" } );

   userAttributesToImport["Test0"] = int( 0 );
   userAttributesToImport["Test1"] = long( 0 );
   userAttributesToImport["Test2"] = float( 0.0 );
   userAttributesToImport["Test3"] = double( 0.0 );
   userAttributesToImport["Test4"] = uint_t( 0 );
   userAttributesToImport["Test5"] = bool( false );
   userAttributesToImport["Test6"] = std::string( "" );
   userAttributesToImport["Test7"]  = std::vector< int >( {} );
   userAttributesToImport["Test8"]  = std::vector< long >( {} );
   userAttributesToImport["Test9"]  = std::vector< float >( {} );
   userAttributesToImport["Test10"] = std::vector< double >( {} );
   userAttributesToImport["Test11"] = std::vector< uint_t >( {} );
   userAttributesToImport["Test12"] = std::vector< bool >( {} );
   userAttributesToImport["Test13"] = std::vector< std::string >( {} );   

   //  Create Checkpoint
   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * exporting checkpoint" );
   }
   func_t< value_t > funcOriginal =
       exportCheckpoint< func_t, value_t >( filePath, fileName, storage, minLevel, maxLevel, userAttributes );

   //  Import Checkpoint
   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " * importing checkpoint" );
   }

   func_t< value_t > funcRestored =
       importCheckpoint< func_t, value_t >( filePath, fileName, storage, minLevel, maxLevel, userAttributesToImport );

   // std::cout << std::get< real_t >(userAttributesToImport["Test1"]) << ", " << std::get< uint_t >(userAttributesToImport["Test2"]) << ", " << std::get< bool >(userAttributesToImport["Test3"]) << std::endl;

   // integer datatype for output
   using intData_t = ADIOS2_PARAVIEW_INT_TYPE;

   // check userAttributes after reload
   for ( auto& entry : userAttributes )
   {
      std::visit(
          [&userAttributesToImport, &entry]( const auto& arg ) {
             using T                = std::decay_t< decltype( arg ) >;
             const std::string& key = entry.first;

             if constexpr ( std::is_same_v< T, std::vector< int > > || std::is_same_v< T, std::vector< long > > ||
                            std::is_same_v< T, std::vector< float > > || std::is_same_v< T, std::vector< double > > ||
                            std::is_same_v< T, std::vector< uint_t > > || std::is_same_v< T, std::vector< bool > > ||
                            std::is_same_v< T, std::vector< std::string > > )
             {
                for ( uint_t k = 0; k < arg.size(); k++ )
                {
                   WALBERLA_LOG_INFO_ON_ROOT( key << " ; " << k << " ; " << arg.at( k ) << " ; "
                                                  << std::get< T >( userAttributesToImport[key] ).at( k ) );
                   WALBERLA_CHECK( arg.at( k ) == std::get< T >( userAttributesToImport[key] ).at( k ) );
                }
             }
             else
             {
                WALBERLA_LOG_INFO_ON_ROOT( key << " ; " << arg << " ; " << std::get< T >( userAttributesToImport[key] ) );
                WALBERLA_CHECK( arg == std::get< T >( userAttributesToImport[key] ) );
             }
          },
          entry.second );
   }

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
         // AdiosWriter adiosWriter( ".", "CheckpointRestoreTestOutput", storage );
         // adiosWriter.add( funcOriginal );
         // adiosWriter.add( funcRestored );
         // adiosWriter.add( difference );
         // adiosWriter.write( lvl );
      }
   }

   if ( testContinuous )
   {
      if ( verbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( " * exporting checkpoint continuous" );
      }
      auto funcListTimePairExported = exportCheckpointContinuous< func_t, value_t >(
          filePath, "Continuous" + fileName, storage, minLevel, maxLevel, nStepsContinuous );

      if ( verbose )
      {
         WALBERLA_LOG_INFO_ON_ROOT( " * importing checkpoint continuous" );
      }

      auto funcListTimePairImported =
          importCheckpointContinuous< func_t, value_t >( filePath, "Continuous" + fileName, storage, minLevel, maxLevel );

      WALBERLA_CHECK_EQUAL( nStepsContinuous, funcListTimePairImported.first.size() );
      WALBERLA_CHECK_EQUAL( nStepsContinuous, funcListTimePairImported.second.size() );

      std::string       funcNameExport = "TestWriteExport_" + FunctionTrait< func_t< value_t > >::getTypeName();
      func_t< value_t > feFuncExport( funcNameExport, storage, minLevel, maxLevel );

      std::string       funcNameImport = "TestWriteImport_" + FunctionTrait< func_t< value_t > >::getTypeName();
      func_t< value_t > feFuncImport( funcNameImport, storage, minLevel, maxLevel );

      AdiosWriter adiosWriterContinuous( ".", "ContinuousCheckpointRestoreTestOutput", storage );
      adiosWriterContinuous.add( feFuncExport );
      adiosWriterContinuous.add( feFuncImport );
      adiosWriterContinuous.add( difference );

      for ( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
      {
         if ( verbose )
         {
            WALBERLA_LOG_INFO_ON_ROOT( " * checking differences on refinement level " << lvl << " ..." );
         }
         for ( uint_t tStep = 0U; tStep < nStepsContinuous; tStep++ )
         {
            if ( lvl == minLevel )
            {
               // Check this only once
               WALBERLA_CHECK_LESS_EQUAL( funcListTimePairExported.second[tStep] - funcListTimePairImported.second[tStep],
                                          real_c( 0 ) );
            }

            difference.assign( { static_cast< value_t >( 1 ), static_cast< value_t >( -1 ) },
                               { *( funcListTimePairExported.first[tStep] ), *( funcListTimePairImported.first[tStep] ) },
                               lvl,
                               All );
            value_t error = computeCheckValue( difference, lvl );

            WALBERLA_CHECK_LESS_EQUAL( error, static_cast< value_t >( 0 ) );

            if ( exportFuncs )
            {
               feFuncExport.assign( { static_cast< value_t >( 1 ) }, { *( funcListTimePairExported.first[tStep] ) }, lvl, All );
               feFuncImport.assign( { static_cast< value_t >( 1 ) }, { *( funcListTimePairImported.first[tStep] ) }, lvl, All );
               adiosWriterContinuous.write( lvl, tStep );
            }
         }
         if ( verbose )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "   ... okay" );
         }
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
      WALBERLA_LOG_INFO_ON_ROOT( " * importing checkpoint from '" << filePath << "/" << fileName << "'" );
   }

   std::map< std::string, adiosHelpers::adiostype_t > userAttributes = {};

   func_t< value_t > funcRestored =
       importCheckpoint< func_t, value_t >( filePath, fileName, storage, minLevel, maxLevel, userAttributes );

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

   // either no args, only parameter file, or parameter file + 2 overrides
   if ( argc != 1 && argc != 2 && argc != 4 )
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

   std::string       filePath{ "." };
   std::string       fileName{ "CheckpointRestoreTest" };
   std::stringstream sStr;
   sStr << "-np" << walberla::mpi::MPIManager::instance()->numProcesses() << ".bp";
   fileName += sStr.str();
   std::string meshFile2D{ prependHyTeGMeshDir( "2D/LShape_6el.msh" ) };

   // Some issue with AdiosWriter, so cannot visualise with this mesh!
   std::string meshFile{ prependHyTeGMeshDir( "3D/cube_6el.msh" ) };

   bool        onlyImport = false;
   std::string fileNameForRestore{ "non-existant-file" };
   if ( argc > 1 )
   {
      walberla::Config::BlockHandle parameters = walberlaEnv.config()->getOneBlock( "Parameters" );
      onlyImport                               = parameters.getParameter< bool >( "onlyImport" );
      fileNameForRestore                       = parameters.getParameter< std::string >( "restoreFromFileWithName" );
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

      runTestWithIdenticalCommunicator< P2VectorFunction, real_t >(
          filePath, fileName, meshFile, minLevel, maxLevel, true, true, 4U );
      runTestWithIdenticalCommunicator< P2Function, real_t >( filePath, fileName, meshFile, minLevel, maxLevel, true, true, 4U );
      runTestWithIdenticalCommunicator< P2Function, real_t >(
          filePath, fileName, meshFile2D, minLevel, maxLevel, true, true, 4U );

      // we are going to reuse this checkpoint in the next part of this pipeline job
      runTestWithIdenticalCommunicator< P2Function, real_t >( filePath, fileName, meshFile, minLevel, maxLevel, true );

      // We currently would need to import the two component functions separately; Better than having specialised code here,
      // alter AdiosCheckpoint[Ex|Im]porter meshFile2D
      // runTestWithIdenticalCommunicator< P2P1TaylorHoodFunction, real_t >( filePath, fileName, meshFile, minLevel, maxLevel, true );
   }

   // =====================
   //  Run Other Test Form
   // =====================
   else
   {
      if ( fileNameForRestore == "non-existant-file" )
      {
         WALBERLA_ABORT( "You need to specify 'restoreFromFileWithName' either in the *.prm file or as CLI override" );
      }
      runTestWithOtherCommunicator< P2Function, real_t >( filePath, fileNameForRestore, meshFile, minLevel, maxLevel, true );
      return EXIT_SUCCESS;
   }
}

// ensure speficied interfaces exist by making compiler explicitely instantiate the CRTP "base" class
template class hyteg::CheckpointExporter< hyteg::AdiosCheckpointExporter >;
template class hyteg::CheckpointImporter< hyteg::AdiosCheckpointImporter >;

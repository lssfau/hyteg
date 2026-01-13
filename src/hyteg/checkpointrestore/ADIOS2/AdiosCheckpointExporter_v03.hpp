/*
 * Copyright (c) 2023-2025 Marcus Mohr, Ponsuganth Ilangovan.
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

#include <adios2.h>

#include "core/DataTypes.h"

#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointHelpers.hpp"
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointHelpers_v03.hpp"
#include "hyteg/checkpointrestore/CheckpointExporter.hpp"
#include "hyteg/checkpointrestore/CheckpointImporter.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosHelperFunctions.hpp"
#include "hyteg/functions/FEFunctionRegistry.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

/// Driver class for storing checkpoints with ADIOS2
///
/// This class allows to write checkpoints of FE-functions to files using functionality offered
/// by ADIOS2. Output is using the BP format (currently version BP5). The class currently only supports
/// the following types of functions:
/// - P1Function and P1VectorFunction
/// - P2Function and P2VectorFunction
/// - P2P1TaylorHoodFunction
class AdiosCheckpointExporter_v03 : public CheckpointExporter< AdiosCheckpointExporter_v03 >
{
 public:
#ifdef HYTEG_BUILD_WITH_MPI

   /// \param configFile        name of a file in XML or YAML format with runtime configuration parameters for ADIOS2
   /// \param comm              MPI Communicator, defaults to the HyTeG standard communicator
   AdiosCheckpointExporter_v03( std::string configFile, MPI_Comm comm = walberla::MPIManager::instance()->comm() )
   {
      // setup central ADIOS2 interface object
      if ( configFile.empty() )
      {
         adios_ = adios2::ADIOS( comm );
      }
      else
      {
         adios_ = adios2::ADIOS( configFile, comm );
      }
   }

#else

   /// \param configFile        name of a file in XML or YAML format with runtime configuration parameters for ADIOS2
   AdiosCheckpointExporter_v03( std::string configFile )
   {
      // setup central ADIOS2 interface object
      if ( configFile.empty() )
      {
         adios_ = adios2::ADIOS();
      }
      else
      {
         adios_ = adios2::ADIOS( configFile );
      }
   }

#endif

   /// Register an FE Function to be included into checkpoints
   ///
   /// By calling this method the passed function object we be included into all future checkpoints.
   /// Data will be stored for all levels from minLevel up to maxLevel.
   template < template < typename > class func_t, typename value_t >
   inline void registerFunction( const func_t< value_t >& function, uint_t minLevel, uint_t maxLevel )
   {
      WALBERLA_ASSERT( minLevel <= maxLevel );
      if constexpr ( std::is_same_v< func_t< value_t >, P1Function< value_t > > ||
                     std::is_same_v< func_t< value_t >, P2Function< value_t > > ||
                     std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > ||
                     std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > )
      {
         feFunctionRegistry_.add( function );
         functionMinLevel_[function.getFunctionName()] = minLevel;
         functionMaxLevel_[function.getFunctionName()] = maxLevel;
      }
      else if constexpr ( std::is_same_v< func_t< value_t >, P2P1TaylorHoodFunction< value_t > > )
      {
         std::string uComponent = function.uvw().getFunctionName();
         std::string pComponent = function.p().getFunctionName();

         feFunctionRegistry_.add( function.uvw() );
         functionMinLevel_[uComponent] = minLevel;
         functionMaxLevel_[uComponent] = maxLevel;

         feFunctionRegistry_.add( function.p() );
         functionMinLevel_[pComponent] = minLevel;
         functionMaxLevel_[pComponent] = maxLevel;
      }
      else
      {
         WALBERLA_ABORT( "AdiosCheckpointExporter::registerFunction() called with, as of now, unsupported function type!" );
      }
   }

   /// Deregister an FE Function to be no longer included into checkpoints
   ///
   /// By calling this method the passed function object will be excluded from future checkpoints.
   template < template < typename > class func_t, typename value_t >
   inline void deregisterFunction( const func_t< value_t >& function )
   {
      if constexpr ( std::is_same_v< func_t< value_t >, P1Function< value_t > > ||
                     std::is_same_v< func_t< value_t >, P2Function< value_t > > ||
                     std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > ||
                     std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > )
      {
         feFunctionRegistry_.remove( function );
         size_t numDel = 0;
         numDel        = functionMinLevel_.erase( function.getFunctionName() );
         WALBERLA_ASSERT( numDel == 1 );
         numDel = functionMaxLevel_.erase( function.getFunctionName() );
         WALBERLA_ASSERT( numDel == 1 );
      }
      else if constexpr ( std::is_same_v< func_t< value_t >, P2P1TaylorHoodFunction< value_t > > )
      {
         feFunctionRegistry_.remove( function.uvw() );
         feFunctionRegistry_.remove( function.p() );

         std::string uComponent = function.getFunctionName() + "_uvw";
         std::string pComponent = function.getFunctionName() + "_p";
         size_t      numDel     = 0;

         numDel = functionMinLevel_.erase( uComponent );
         WALBERLA_ASSERT( numDel == 1 );
         numDel = functionMaxLevel_.erase( uComponent );
         WALBERLA_ASSERT( numDel == 1 );
         numDel = functionMinLevel_.erase( pComponent );
         WALBERLA_ASSERT( numDel == 1 );
         numDel = functionMaxLevel_.erase( pComponent );
         WALBERLA_ASSERT( numDel == 1 );
      }
      else
      {
         WALBERLA_ABORT( "AdiosCheckpointExporter::deregisterFunction() called with, as of now, unsupported function type!" );
      }
   }

   /// Trigger storing of a single checkpoint
   ///
   /// \param filePath             path to checkpoint file
   /// \param fileName             name of checkpoint "file" (BP format actually uses a directory)
   inline void storeCheckpoint( std::string filePath, std::string fileName )
   {
      std::vector< std::string >                         userAttributeNames;
      std::vector< std::string >                         userAttributeValues;
      std::map< std::string, adiosHelpers::adiostype_t > userDefinedAttributes;
      storeCheckpoint( filePath, fileName, userDefinedAttributes );
   }

   /// Trigger storing of a single checkpoint
   ///
   /// \param filePath             path to checkpoint file
   /// \param fileName             name of checkpoint "file" (BP format actually uses a directory)
   /// \param userAttributeNames   list of names for additional attributes
   /// \param userAttributeValues  list of strings with data for the additional attributes
   inline void storeCheckpoint( std::string                                     filePath,
                                std::string                                     fileName,
                                const std::vector< std::string >&               userAttributeNames,
                                const std::vector< adiosHelpers::adiostype_t >& userAttributeValues )
   {
      WALBERLA_CHECK( userAttributeNames.size() == userAttributeValues.size() );
      std::map< std::string, adiosHelpers::adiostype_t > userDefinedAttributes;

      for ( uint_t k = 0; k < userAttributeNames.size(); k++ )
      {
         userDefinedAttributes.emplace( userAttributeNames[k], userAttributeValues[k] );
      }

      storeCheckpoint( filePath, fileName, userDefinedAttributes );
   }

   /// Trigger storing of a single checkpoint
   ///
   /// \param filePath              path to checkpoint file
   /// \param fileName              name of checkpoint "file" (BP format actually uses a directory)
   /// \param userDefinedAttributes mapping of string keys to values of adiostype_t
   inline void storeCheckpoint( std::string                                               filePath,
                                std::string                                               fileName,
                                const std::map< std::string, adiosHelpers::adiostype_t >& userDefinedAttributes )
   {
      this->doStoreCheckpoint(
          filePath, fileName, userDefinedAttributes, "AdiosCheckpointSingle", false, true, false, real_c( -1 ) );
   };

   inline void storeCheckpointContinuous( std::string filePath, std::string fileName, real_t time, bool finalCall = false )
   {
      std::vector< std::string >               userAttributeNames;
      std::vector< adiosHelpers::adiostype_t > userAttributeValues;
      storeCheckpointContinuous( filePath, fileName, userAttributeNames, userAttributeValues, time, finalCall );
   }

   inline void storeCheckpointContinuous( std::string                                     filePath,
                                          std::string                                     fileName,
                                          const std::vector< std::string >&               userAttributeNames,
                                          const std::vector< adiosHelpers::adiostype_t >& userAttributeValues,
                                          real_t                                          time,
                                          bool                                            finalCall = false )
   {
      WALBERLA_CHECK( userAttributeNames.size() == userAttributeValues.size() );
      std::map< std::string, adiosHelpers::adiostype_t > userDefinedAttributes;

      for ( uint_t k = 0; k < userAttributeNames.size(); k++ )
      {
         userDefinedAttributes.emplace( userAttributeNames[k], userAttributeValues[k] );
      }

      storeCheckpointContinuous( filePath, fileName, userDefinedAttributes, time, finalCall );
   }

   inline void storeCheckpointContinuous( std::string                                               filePath,
                                          std::string                                               fileName,
                                          const std::map< std::string, adiosHelpers::adiostype_t >& userDefinedAttributes,
                                          real_t                                                    time,
                                          bool                                                      finalCall = false )
   {
      this->doStoreCheckpoint(
          filePath, fileName, userDefinedAttributes, "AdiosCheckpointContinuous", true, finalCall, true, time );
   };

   /// Type of engine to be used for export
   inline static const std::string engineType_{ ADIOS2_CHECKPOINT_FORMAT };

 private:
   /// object that remembers the functions we should export
   FEFunctionRegistry feFunctionRegistry_;

   /// map to remember the minLevel for a registered function
   std::map< std::string, uint_t > functionMinLevel_;

   /// map to remember the maxLevel for a registered function
   std::map< std::string, uint_t > functionMaxLevel_;

   /// central ADIOS2 interface object
   adios2::ADIOS adios_;

   /// central ADIOS2 IO object
   std::shared_ptr< adios2::IO > ptrToIO_;

   /// central ADIOS2 Engine object
   adios2::Engine engine_;

   /// remember if we already had a storeCheckpointContinuous() episode
   bool firstWriteDidHappen_ = false;

   /// auxilliary variable to add management information to checkpoint
   ///@{
   std::vector< std::string > allFunctionNames_;
   std::vector< std::string > allFunctionKinds_;
   std::vector< std::string > allFunctionValueTypes_;

   std::map< std::string, adios2::IO > ioObjectsContainer;
   std::map< std::string, int >        ioObjectsCounter;
   ///@}

   enum class ExportType
   {
      ONLY_DEFINE,
      ONLY_EXPORT,
      DEFINE_AND_EXPORT
   };

   inline void doStoreCheckpoint( std::string                                               filePath,
                                  std::string                                               fileName,
                                  const std::map< std::string, adiosHelpers::adiostype_t >& userDefinedAttributes,
                                  const std::string&                                        engineName,
                                  bool                                                      runContinuous,
                                  bool                                                      finalCall,
                                  bool                                                      storeTime,
                                  real_t                                                    time )
   {
      // create the writer and engine for the export
      std::string cpFileName = filePath + "/" + fileName;

      // our ADIOS2 IO object for this checkpoint export
      std::shared_ptr< adios2::IO > ptrToIO;
      bool                          firstWriteDidHappen = false;
      if ( runContinuous )
      {
         ptrToIO             = ptrToIO_;
         firstWriteDidHappen = firstWriteDidHappen_;
      }

      if ( !firstWriteDidHappen )
      {
         ptrToIO = std::make_shared< adios2::IO >( adios_.DeclareIO( engineName ) );
         if ( runContinuous )
         {
            ptrToIO_ = ptrToIO;
         }

         ptrToIO->SetEngine( engineType_ );
         engine_ = ptrToIO->Open( cpFileName, adios2::Mode::Write );

         // export meta-data
         adiosHelpers::generateSoftwareMetaData( *ptrToIO );
         addVersionInformation( *ptrToIO );

         if ( storeTime )
         {
            ptrToIO->DefineVariable< real_t >( "TIME" );
         }

         // generate variables for export
         defineAndOrExportVariables< P1Function, real_t >( *ptrToIO, engine_, ExportType::ONLY_DEFINE );
         defineAndOrExportVariables< P1Function, int32_t >( *ptrToIO, engine_, ExportType::ONLY_DEFINE );
         defineAndOrExportVariables< P1Function, int64_t >( *ptrToIO, engine_, ExportType::ONLY_DEFINE );

         defineAndOrExportVariables< P1VectorFunction, real_t >( *ptrToIO, engine_, ExportType::ONLY_DEFINE );
         defineAndOrExportVariables< P1VectorFunction, int32_t >( *ptrToIO, engine_, ExportType::ONLY_DEFINE );
         defineAndOrExportVariables< P1VectorFunction, int64_t >( *ptrToIO, engine_, ExportType::ONLY_DEFINE );

         defineAndOrExportVariables< P2Function, real_t >( *ptrToIO, engine_, ExportType::ONLY_DEFINE );
         defineAndOrExportVariables< P2Function, int32_t >( *ptrToIO, engine_, ExportType::ONLY_DEFINE );
         defineAndOrExportVariables< P2Function, int64_t >( *ptrToIO, engine_, ExportType::ONLY_DEFINE );

         defineAndOrExportVariables< P2VectorFunction, real_t >( *ptrToIO, engine_, ExportType::ONLY_DEFINE );
         defineAndOrExportVariables< P2VectorFunction, int32_t >( *ptrToIO, engine_, ExportType::ONLY_DEFINE );
         defineAndOrExportVariables< P2VectorFunction, int64_t >( *ptrToIO, engine_, ExportType::ONLY_DEFINE );

         firstWriteDidHappen_ = true;
      }

      // start the export episode
      engine_.BeginStep();

      // define variables for export
      defineAndOrExportVariables< P1Function, real_t >( *ptrToIO, engine_, ExportType::ONLY_EXPORT );
      defineAndOrExportVariables< P1Function, int32_t >( *ptrToIO, engine_, ExportType::ONLY_EXPORT );
      defineAndOrExportVariables< P1Function, int64_t >( *ptrToIO, engine_, ExportType::ONLY_EXPORT );

      defineAndOrExportVariables< P1VectorFunction, real_t >( *ptrToIO, engine_, ExportType::ONLY_EXPORT );
      defineAndOrExportVariables< P1VectorFunction, int32_t >( *ptrToIO, engine_, ExportType::ONLY_EXPORT );
      defineAndOrExportVariables< P1VectorFunction, int64_t >( *ptrToIO, engine_, ExportType::ONLY_EXPORT );

      defineAndOrExportVariables< P2Function, real_t >( *ptrToIO, engine_, ExportType::ONLY_EXPORT );
      defineAndOrExportVariables< P2Function, int32_t >( *ptrToIO, engine_, ExportType::ONLY_EXPORT );
      defineAndOrExportVariables< P2Function, int64_t >( *ptrToIO, engine_, ExportType::ONLY_EXPORT );

      defineAndOrExportVariables< P2VectorFunction, real_t >( *ptrToIO, engine_, ExportType::ONLY_EXPORT );
      defineAndOrExportVariables< P2VectorFunction, int32_t >( *ptrToIO, engine_, ExportType::ONLY_EXPORT );
      defineAndOrExportVariables< P2VectorFunction, int64_t >( *ptrToIO, engine_, ExportType::ONLY_EXPORT );

      if ( storeTime )
      {
         auto varTimeStep = ptrToIO->InquireVariable< real_t >( "TIME" );
         engine_.Put( varTimeStep, time );
      }

      if ( finalCall )
      {
         // add user defined attributes
         adiosHelpers::writeAllAttributes( *ptrToIO, userDefinedAttributes );

         // add attributes with meta information on functions in the checkpoint
         ptrToIO->DefineAttribute< std::string >( "FunctionNames", allFunctionNames_.data(), allFunctionNames_.size() );
         ptrToIO->DefineAttribute< std::string >( "FunctionKinds", allFunctionKinds_.data(), allFunctionKinds_.size() );
         ptrToIO->DefineAttribute< std::string >(
             "FunctionValueTypes", allFunctionValueTypes_.data(), allFunctionValueTypes_.size() );
         std::vector< uint_t > allMinLevels;
         std::vector< uint_t > allMaxLevels;
         for ( const auto& funcName : allFunctionNames_ )
         {
            allMinLevels.push_back( functionMinLevel_.at( funcName ) );
            allMaxLevels.push_back( functionMaxLevel_.at( funcName ) );
         }
         ptrToIO->DefineAttribute< uint_t >( "FunctionMinLevels", allMinLevels.data(), allMinLevels.size() );
         ptrToIO->DefineAttribute< uint_t >( "FunctionMaxLevels", allMaxLevels.data(), allMaxLevels.size() );
      }

      // actual export performed here (if lazy not overwritten in config file)
      engine_.EndStep();

      if ( finalCall )
      {
         engine_.Close();

         // clean-up for next checkpoint
         allFunctionNames_.clear();
         allFunctionKinds_.clear();
      }
   };

   template < template < typename > class func_t, typename value_t >
   void
       defineAndOrExportVariables( adios2::IO& io, adios2::Engine& engine, ExportType exportType = ExportType::DEFINE_AND_EXPORT )
   {
      if ( exportType != ExportType::DEFINE_AND_EXPORT && exportType != ExportType::ONLY_DEFINE &&
           exportType != ExportType::ONLY_EXPORT )
      {
         WALBERLA_ABORT( "Unsupported enumeration value for exportType argument detected!" );
      }

      // extract all functions of given kind and all value types
      const FunctionMultiStore< func_t >& functionList = feFunctionRegistry_.getFunctions< func_t >();

      // extract functions for desired value type
      const std::vector< func_t< value_t > >& funcs = functionList.template getFunctions< value_t >();

      // for each FE function call specialised free C++ function
      for ( const auto& function : funcs )
      {
         // WALBERLA_LOG_INFO_ON_ROOT( "--> Checkpointing '" << function.getFunctionName() << "'" );
         WALBERLA_ASSERT( functionMinLevel_.at( function.getFunctionName() ) >= 0 );
         WALBERLA_ASSERT( functionMaxLevel_.at( function.getFunctionName() ) >= 0 );

         if ( !firstWriteDidHappen_ )
         {
            // add information on function for later attribute generation
            allFunctionNames_.push_back( function.getFunctionName() );
            allFunctionKinds_.push_back( FunctionTrait< func_t< value_t > >::getTypeName() );
            allFunctionValueTypes_.push_back( adiosCheckpointHelpers::valueTypeToString< value_t >() );
         }

         if constexpr ( std::is_same_v< func_t< value_t >, P1Function< value_t > > ||
                        std::is_same_v< func_t< value_t >, P2Function< value_t > > ||
                        std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > ||
                        std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > )
         {
            if ( exportType == ExportType::DEFINE_AND_EXPORT || exportType == ExportType::ONLY_DEFINE )
            {
               // enumerate primitives, if necessary
               globallyEnumeratePrimitives( function.getStorage() );

               // first define the variable
               for ( uint_t level = functionMinLevel_[function.getFunctionName()];
                     level <= functionMaxLevel_[function.getFunctionName()];
                     ++level )
               {
                  uint_t globalNumberOfHighestDimPrimitives = function.getStorage()->hasGlobalCells() ?
                                                                  numFacesAndCellsGlobally_.second :
                                                                  numFacesAndCellsGlobally_.first;
                  WALBERLA_ASSERT( globalNumberOfHighestDimPrimitives != 0u, "Severe problem detected!" );
                  adiosCheckpointHelpers_v03::generateVariables( io, function, level, globalNumberOfHighestDimPrimitives );
               }
            }

            if ( exportType == ExportType::DEFINE_AND_EXPORT || exportType == ExportType::ONLY_EXPORT )
            {
               // export primitive-to-index map, if necessary
               exportPrimitiveToIndexMap( io, engine, function.getStorage() );

               // as we only store data for the highest dimensional primitives we need to make sure
               // that their halos are up-to-date
               for ( uint_t level = functionMinLevel_[function.getFunctionName()];
                     level <= functionMaxLevel_[function.getFunctionName()];
                     ++level )
               {
                  if constexpr ( std::is_same_v< func_t< value_t >, P1Function< value_t > > ||
                                 std::is_same_v< func_t< value_t >, P2Function< value_t > > )
                  {
                     communication::syncFunctionBetweenPrimitives( function, level, communication::syncDirection_t::LOW2HIGH );
                  }
                  else
                  {
                     communication::syncVectorFunctionBetweenPrimitives(
                         function, level, communication::syncDirection_t::LOW2HIGH );
                  }
               }

               // now schedule the variable for export
               adiosCheckpointHelpers_v03::doSomethingForAFunctionOnAllHighestDimensionalPrimitives(
                   io,
                   engine,
                   function,
                   functionMinLevel_[function.getFunctionName()],
                   functionMaxLevel_[function.getFunctionName()],
                   function.getStorage()->hasGlobalCells() ? rowIndicesInGlobalAdiosArrayForCells_ :
                                                             rowIndicesInGlobalAdiosArrayForFaces_,
                   adiosCheckpointHelpers_v03::exportVariables< func_t, value_t > );
            }
         }
         else
         {
            WALBERLA_ABORT( "Achievement unlocked: 'Detector of the Missing Implementation'!" );
         }
      }
   }

   inline void globallyEnumeratePrimitives( const std::shared_ptr< const PrimitiveStorage >& storage )
   {
      if ( storage->hasGlobalCells() && numFacesAndCellsGlobally_.second == 0u )
      {
         WALBERLA_LOG_PROGRESS_ON_ROOT( "Performing global cell enumeration" );
         auto [localCellIndices, numGlobalCells] = adiosCheckpointHelpers_v03::enumerateCells( storage );
         numFacesAndCellsGlobally_.second        = numGlobalCells;
         rowIndicesInGlobalAdiosArrayForCells_   = localCellIndices;
      }
      else if ( !storage->hasGlobalCells() && numFacesAndCellsGlobally_.first == 0u )
      {
         WALBERLA_LOG_PROGRESS_ON_ROOT( "Performing global face enumeration" );
         auto [localFaceIndices, numGlobalFaces] = adiosCheckpointHelpers_v03::enumerateFaces( storage );
         numFacesAndCellsGlobally_.first         = numGlobalFaces;
         rowIndicesInGlobalAdiosArrayForFaces_   = localFaceIndices;
      }
   }

   /// Provide information on version of checkpoint (currently 0.3)
   inline void addVersionInformation( adios2::IO& io )
   {
      std::vector< std::string > info{ "HyTeG Checkpoint", "0.3" };
      io.DefineAttribute< std::string >( "CheckpointFormat", info.data(), info.size() );
   }

   /// Store information in which array row to find data for a certain PrimitiveID
   inline void exportPrimitiveToIndexMap( adios2::IO&                                      io,
                                          adios2::Engine&                                  engine,
                                          const std::shared_ptr< const PrimitiveStorage >& storage )
   {
      // test whether we need to store a map for faces
      if ( numFacesAndCellsGlobally_.first > 0u )
      {
         if ( !io.InquireVariable< uint8_t >( "RowIndicesOfFaces" ) )
         {
            WALBERLA_LOG_PROGRESS_ON_ROOT( "Need to store faceID to rowIndex map!" );

            uint_t mapEntrySizeInBytes = sizeof( PrimitiveID ) + sizeof( uint_t );

            WALBERLA_LOG_PROGRESS_ON_ROOT( "Going to define joined array 'RowIndicesOfFaces' with\nshape .... {"
                                           << "adios2::JoinedDim, " << mapEntrySizeInBytes << "}\n"
                                           << "start .... {}\n"
                                           << "count .... {" << storage->getNumberOfLocalFaces() << ", " << mapEntrySizeInBytes
                                           << "}" );

            adios2::Variable< uint8_t > varMapData =
                io.DefineVariable< uint8_t >( "RowIndicesOfFaces",
                                              { adios2::JoinedDim, mapEntrySizeInBytes },
                                              {},
                                              { storage->getNumberOfLocalFaces(), mapEntrySizeInBytes } );

            // serialise map
            walberla::mpi::SendBuffer buffer;
            buffer << rowIndicesInGlobalAdiosArrayForFaces_;

            // schedule map data for export
            std::ptrdiff_t offset = buffer.size() - storage->getNumberOfLocalFaces() * mapEntrySizeInBytes;
            engine.Put( varMapData, buffer.ptr() + offset );
         }
      }

      // test whether we need to store a map for cells
      if ( numFacesAndCellsGlobally_.second > 0u )
      {
         std::string varName = "RowIndicesOfCells";
         if ( !io.InquireVariable< uint8_t >( varName ) )
         {
            WALBERLA_LOG_PROGRESS_ON_ROOT( "Need to store cellID to rowIndex map!" );

            uint_t mapEntrySizeInBytes = sizeof( PrimitiveID ) + sizeof( uint_t );

            WALBERLA_LOG_PROGRESS_ON_ROOT( "Going to define joined array '"
                                           << varName << "' with\nshape .... {"
                                           << "adios2::JoinedDim, " << mapEntrySizeInBytes << "}\n"
                                           << "start .... {}\n"
                                           << "count .... {" << storage->getNumberOfLocalCells() << ", " << mapEntrySizeInBytes
                                           << "}" );

            // serialise map
            walberla::mpi::SendBuffer buffer;
            buffer << rowIndicesInGlobalAdiosArrayForCells_;

            adios2::Variable< uint8_t > varMapData =
                io.DefineVariable< uint8_t >( varName,
                                              { adios2::JoinedDim, mapEntrySizeInBytes },
                                              {},
                                              { storage->getNumberOfLocalCells(), mapEntrySizeInBytes } );

            // schedule map data for export
            std::ptrdiff_t offset = buffer.size() - storage->getNumberOfLocalCells() * mapEntrySizeInBytes;
            engine.Put( varMapData, buffer.ptr() + offset );
         }
      }
   }

   /// Global count of faces or cells
   std::pair< uint_t, uint_t > numFacesAndCellsGlobally_{ 0u, 0u };

   /// Local indices of faces within their global enumeration
   std::map< PrimitiveID, uint_t > rowIndicesInGlobalAdiosArrayForFaces_;

   /// Local indices of cells within their global enumeration
   std::map< PrimitiveID, uint_t > rowIndicesInGlobalAdiosArrayForCells_;
};

} // namespace hyteg

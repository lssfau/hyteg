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
#pragma once

#include <adios2.h>

#include "core/DataTypes.h"

#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointHelpers.hpp"
#include "hyteg/checkpointrestore/CheckpointExporter.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosHelperFunctions.hpp"
#include "hyteg/dataexport/FEFunctionRegistry.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

/// Base class for driver classes that store checkpoints
class AdiosCheckpointExporter : CheckpointExporter< AdiosCheckpointExporter >
{
 public:
#ifdef WALBERLA_BUILD_WITH_MPI

   AdiosCheckpointExporter( std::string configFile, MPI_Comm comm = walberla::MPIManager::instance()->comm() )
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

   AdiosCheckpointExporter( std::string configFile )
   {
      // setup central ADIOS2 interface object
      if ( configFile.empty() )
      {
         adios_ = adios2::ADIOS;
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
         std::string uComponent = function.getFunctionName() + "_uvw";
         std::string pComponent = function.getFunctionName() + "_p";

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
      std::vector< std::string > userAttributeNames;
      std::vector< std::string > userAttributeValues;
   }

   /// Trigger storing of a single checkpoint
   ///
   /// \param filePath             path to checkpoint file
   /// \param fileName             name of checkpoint "file" (BP format actually uses a directory)
   /// \param userAttributeNames   list of names for additional attributes
   /// \param userAttributeValues  list of strings with data for the additional attributes
   inline void storeCheckpoint( std::string                       filePath,
                                std::string                       fileName,
                                const std::vector< std::string >& userAttributeNames,
                                const std::vector< std::string >& userAttributeValues )
   {
      // create the writer and engine for the export
      std::string cpFileName = filePath + "/" + fileName;
      adios2::IO  io         = adios_.DeclareIO( "AdiosCheckpointExport" );
      io.SetEngine( "BP5" );
      adios2::Engine engine = io.Open( cpFileName, adios2::Mode::Write );

      // start the export episode
      engine.BeginStep();

      // export meta-data
      adiosHelpers::generateSoftwareMetaData( io );
      addVersionInformation( io );

      // schedule data for export
      defineAndExportVariables< P1Function, real_t >( io, engine );
      defineAndExportVariables< P1Function, int32_t >( io, engine );
      defineAndExportVariables< P1Function, int64_t >( io, engine );

      defineAndExportVariables< P1VectorFunction, real_t >( io, engine );
      defineAndExportVariables< P1VectorFunction, int32_t >( io, engine );
      defineAndExportVariables< P1VectorFunction, int64_t >( io, engine );

      defineAndExportVariables< P2Function, real_t >( io, engine );
      defineAndExportVariables< P2Function, int32_t >( io, engine );
      defineAndExportVariables< P2Function, int64_t >( io, engine );

      defineAndExportVariables< P2VectorFunction, real_t >( io, engine );
      defineAndExportVariables< P2VectorFunction, int32_t >( io, engine );
      defineAndExportVariables< P2VectorFunction, int64_t >( io, engine );

      // add user defined attributes
      WALBERLA_ASSERT( userAttributeNames.size() == userAttributeValues.size() );
      for ( uint_t k = 0; k < userAttributeNames.size(); ++k )
      {
         io.DefineAttribute< std::string >( userAttributeNames[k], userAttributeValues[k] );
      }

      // add attributes with meta information on functions in the checkpoint
      io.DefineAttribute< std::string >( "FunctionNames", allFunctionNames_.data(), allFunctionNames_.size() );
      io.DefineAttribute< std::string >( "FunctionTypes", allFunctionTypes_.data(), allFunctionTypes_.size() );
      std::vector< uint_t > allMinLevels;
      std::vector< uint_t > allMaxLevels;
      for ( const auto& funcName : allFunctionNames_ )
      {
         allMinLevels.push_back( functionMinLevel_.at( funcName ) );
         allMaxLevels.push_back( functionMaxLevel_.at( funcName ) );
      }
      io.DefineAttribute< uint_t >( "FunctionMinLevels", allMinLevels.data(), allMinLevels.size() );
      io.DefineAttribute< uint_t >( "FunctionMaxLevels", allMaxLevels.data(), allMaxLevels.size() );

      // actual export performed here (if lazy not overwritten in config file)
      engine.EndStep();

      // need to close file before Engine and IO objects get destroyed
      engine.Close();

      // clean-up for next checkpoint
      allFunctionNames_.clear();
      allFunctionTypes_.clear();
   };

 private:
   /// object that remembers the functions we should export
   FEFunctionRegistry feFunctionRegistry_;

   /// map to remember the minLevel for a registered function
   std::map< std::string, uint_t > functionMinLevel_;

   /// map to remember the maxLevel for a registered function
   std::map< std::string, uint_t > functionMaxLevel_;

   /// central ADIOS2 interface object
   adios2::ADIOS adios_;

   /// auxilliary variable to add management information to checkpoint
   ///@{
   std::vector< std::string > allFunctionNames_;
   std::vector< std::string > allFunctionTypes_;
   ///@}

   template < template < typename > class func_t, typename value_t >
   void defineAndExportVariables( adios2::IO& io, adios2::Engine& engine )
   {
      // extract all functions of given kind and all value types
      const FunctionMultiStore< func_t >& functionList = feFunctionRegistry_.getFunctions< func_t >();

      // extract functions for desired value type
      const std::vector< func_t< value_t > >& funcs = functionList.template getFunctions< value_t >();

      // for each FE function call specialised free C++ function
      for ( const auto& function : funcs )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "--> Checkpointing '" << function.getFunctionName() << "'" );
         WALBERLA_ASSERT( functionMinLevel_.at( function.getFunctionName() ) >= 0 );
         WALBERLA_ASSERT( functionMaxLevel_.at( function.getFunctionName() ) >= 0 );

         if constexpr ( std::is_same_v< func_t< value_t >, P1Function< value_t > > ||
                        std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > )
         {
            // first define the variable
            adiosCheckpointHelpers::doSomethingForAFunctionOnAllPrimitives(
                io,
                engine,
                function,
                functionMinLevel_[function.getFunctionName()],
                functionMaxLevel_[function.getFunctionName()],
                adiosCheckpointHelpers::generateVariableForP1TypeFunction< func_t, value_t > );

            // now schedule the variable for export
            adiosCheckpointHelpers::doSomethingForAFunctionOnAllPrimitives(
                io,
                engine,
                function,
                functionMinLevel_[function.getFunctionName()],
                functionMaxLevel_[function.getFunctionName()],
                adiosCheckpointHelpers::exportVariableForP1TypeFunction< func_t, value_t > );

            // add information on function for later attribute generation
            allFunctionNames_.push_back( function.getFunctionName() );
            allFunctionTypes_.push_back( FunctionTrait< func_t< value_t > >::getTypeName() );
         }

         else if constexpr ( std::is_same_v< func_t< value_t >, P2Function< value_t > > ||
                             std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > )
         {
            // first define the variable
            adiosCheckpointHelpers::doSomethingForAFunctionOnAllPrimitives(
                io,
                engine,
                function,
                functionMinLevel_[function.getFunctionName()],
                functionMaxLevel_[function.getFunctionName()],
                adiosCheckpointHelpers::generateVariableForP2TypeFunction< func_t, value_t > );

            // now schedule the variable for export
            adiosCheckpointHelpers::doSomethingForAFunctionOnAllPrimitives(
                io,
                engine,
                function,
                functionMinLevel_[function.getFunctionName()],
                functionMaxLevel_[function.getFunctionName()],
                adiosCheckpointHelpers::exportVariableForP2TypeFunction< func_t, value_t > );

            // add information on function for later attribute generation
            allFunctionNames_.push_back( function.getFunctionName() );
            allFunctionTypes_.push_back( FunctionTrait< func_t< value_t > >::getTypeName() );
         }
         else
         {
            WALBERLA_ABORT( "Achievement unlocked: 'Detector of the Missing Implementation'!" );
         }
      }
   }

   /// Provide information on version of checkpoint (currently 0.1)
   inline void addVersionInformation( adios2::IO& io )
   {
      std::vector< std::string > info{ "HyTeG Checkpoint", "0.1" };
      io.DefineAttribute< std::string >( "CheckpointFormat", info.data(), info.size() );
   }
};

} // namespace hyteg

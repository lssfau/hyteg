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
         functionMinLevel[function.getFunctionName()] = minLevel;
         functionMaxLevel[function.getFunctionName()] = maxLevel;
      }
      else if constexpr ( std::is_same_v< func_t< value_t >, P2P1TaylorHoodFunction< value_t > > )
      {
         std::string uComponent = function.getFunctionName() + "_uvw";
         std::string pComponent = function.getFunctionName() + "_p";

         feFunctionRegistry_.add( function.uvw() );
         functionMinLevel[uComponent] = minLevel;
         functionMaxLevel[uComponent] = maxLevel;

         feFunctionRegistry_.add( function.p() );
         functionMinLevel[pComponent] = minLevel;
         functionMaxLevel[pComponent] = maxLevel;
      }
      else
      {
         WALBERLA_ABORT( "AdiosCheckpointExporter::registerFunction() called with, as of now, unsupported function type!" );
      }
   }

#ifdef FE_FUNCTION_REGISTRY_HAS_REMOVE
   /// Deregister an FE Function to be no longer included into checkpoints
   ///
   /// By calling this method the passed function object will be excluded from future checkpoints.
   template < template < typename > class func_t, typename value_t >
   inline void deregisterFunction( const func_t< value_t >& function )
   {
      feFunctionRegistry_.remove( function );
      size_t numDel = 0;
      numDel        = functionMinLevel.erase( function.getFunctionName() );
      WALBERLA_ASSERT( numDel == 1 );
      numDel = functionMaxLevel.erase( function.getFunctionName() );
      WALBERLA_ASSERT( numDel == 1 );
   }
#endif

   /// Trigger storing of a single checkpoint
   inline void storeCheckpoint( std::string filePath, std::string fileName )
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

      // actual export performed here (if lazy not overwritten in config file)
      engine.EndStep();

      // need to close file before Engine and IO objects get destroyed
      engine.Close();
   };

 private:
   /// object that remembers the functions we should export
   FEFunctionRegistry feFunctionRegistry_;

   /// map to remember the minLevel for a registered function
   std::map< std::string, uint_t > functionMinLevel;

   /// map to remember the maxLevel for a registered function
   std::map< std::string, uint_t > functionMaxLevel;

   /// central ADIOS2 interface object
   adios2::ADIOS adios_;

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
         WALBERLA_ASSERT( functionMinLevel.at( function.getFunctionName() ) >= 0 );
         WALBERLA_ASSERT( functionMaxLevel.at( function.getFunctionName() ) >= 0 );

         if constexpr ( std::is_same_v< func_t< value_t >, P1Function< value_t > > ||
                        std::is_same_v< func_t< value_t >, P1VectorFunction< value_t > > )
         {
            // first define the variable
            adiosCheckpointHelpers::doSomethingForAFunctionOnAllPrimitives(
                io,
                engine,
                function,
                functionMinLevel[function.getFunctionName()],
                functionMaxLevel[function.getFunctionName()],
                adiosCheckpointHelpers::generateVariableForP1TypeFunction< func_t, value_t > );

            // now schedule the variable for export
            adiosCheckpointHelpers::doSomethingForAFunctionOnAllPrimitives(
                io,
                engine,
                function,
                functionMinLevel[function.getFunctionName()],
                functionMaxLevel[function.getFunctionName()],
                adiosCheckpointHelpers::exportVariableForP1TypeFunction< func_t, value_t > );
         }

         else if constexpr ( std::is_same_v< func_t< value_t >, P2Function< value_t > > ||
                             std::is_same_v< func_t< value_t >, P2VectorFunction< value_t > > )
         {
            // first define the variable
            adiosCheckpointHelpers::doSomethingForAFunctionOnAllPrimitives(
                io,
                engine,
                function,
                functionMinLevel[function.getFunctionName()],
                functionMaxLevel[function.getFunctionName()],
                adiosCheckpointHelpers::generateVariableForP2TypeFunction< func_t, value_t > );

            // now schedule the variable for export
            adiosCheckpointHelpers::doSomethingForAFunctionOnAllPrimitives(
                io,
                engine,
                function,
                functionMinLevel[function.getFunctionName()],
                functionMaxLevel[function.getFunctionName()],
                adiosCheckpointHelpers::exportVariableForP2TypeFunction< func_t, value_t > );
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
      io.DefineAttribute< std::string >( "Checkpoint_Format", "0.1", "", "" );
   }
};

} // namespace hyteg

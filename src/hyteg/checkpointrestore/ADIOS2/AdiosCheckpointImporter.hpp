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
#include "hyteg/checkpointrestore/CheckpointImporter.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosHelperFunctions.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

/// Driver class for importing function data from checkpoints with ADIOS2
class AdiosCheckpointImporter : CheckpointImporter< AdiosCheckpointImporter >
{
 public:
#ifdef WALBERLA_BUILD_WITH_MPI

   /// \param filePath             path to checkpoint file
   /// \param fileName             name of checkpoint "file" (BP format actually uses a directory)
   /// \param configFile           name of a config file for on-the-fly ADIOS tuning (string is allowed
   ///                             to be empty, in which case it is ignored)
   /// \param comm                 MPI communicator
   AdiosCheckpointImporter( std::string filePath,
                            std::string fileName,
                            std::string configFile,
                            MPI_Comm    comm = walberla::MPIManager::instance()->comm() )
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

      // create sub-objects for ADIOS
      setupImporter( filePath, fileName );
   }

#else

   /// \param filePath             path to checkpoint file
   /// \param fileName             name of checkpoint "file" (BP format actually uses a directory)
   /// \param configFile           name of a config file for on-the-fly ADIOS tuning (string is allowed
   ///                             to be empty, in which case it is ignored)
   AdiosCheckpointImporter( std::string filePath, std::string fileName, std::string configFile )
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

      // create sub-objects for ADIOS
      setupImporter( filePath, fileName );
   }

#endif

   /// Auxilliary class for storing the checkpoint's meta-info on the included FE functions
   struct FunctionDescription
   {
      std::string name;
      std::string kind;
      std::string valueType;
      uint_t      minLevel;
      uint_t      maxLevel;
   };

   const std::vector< FunctionDescription >& getFunctionDetails() const { return funcDescr_; }

   void printCheckpointInfo()
   {
      WALBERLA_ROOT_SECTION()
      {
         std::stringstream sStream;

         adios2::Attribute< std::string > attrFormat = readAttribute< std::string >( "CheckpointFormat" );
         sStream << "Checkpoint Format = " << attrFormat.Data()[0] << ", version " << attrFormat.Data()[1] << '\n';

         writeFunctionDetails( sStream, "" );
         WALBERLA_LOG_INFO( "" << sStream.rdbuf() );
      }
   }

   template < template < typename > class func_t, typename value_t >
   bool restoreFunction( func_t< value_t >& function, bool abortOnError = true )
   {
      return restoreFunction( function, function.getMinLevel(), function.getMaxLevel(), abortOnError );
   }

   template < template < typename > class func_t, typename value_t >
   bool restoreFunction( func_t< value_t >& function, uint_t minLevel, uint_t maxLevel, bool abortOnError = true )
   {
      // check that function is include in checkpoint and that levels make sense
      WALBERLA_ASSERT( minLevel <= maxLevel );
      bool              restorePossible = false;
      std::stringstream msg;
      std::string       funcName = function.getFunctionName();
      for ( const auto& entry : funcDescr_ )
      {
         if ( entry.name == funcName )
         {
            if ( entry.minLevel > minLevel )
            {
               msg << "minLevel = " << minLevel << " not present in checkpoint!";
               return handleError( funcName, msg, abortOnError );
            }
            if ( entry.maxLevel < maxLevel )
            {
               msg << "maxLevel = " << maxLevel << " not present in checkpoint!";
               return handleError( funcName, msg, abortOnError );
            }
            // std::string vType = adiosCheckpointHelpers::valueTypeToString< typename FunctionTrait< func_t >::ValueType >();
            std::string vType = adiosCheckpointHelpers::valueTypeToString< value_t >();
            if ( vType != entry.valueType )
            {
               msg << "requested valueType = " << vType << " not present in checkpoint!";
               return handleError( funcName, msg, abortOnError );
            }
            restorePossible = true;
            break;
         }
      }
      if ( !restorePossible )
      {
         msg << "function not present in checkpoint!";
         return handleError( funcName, msg, abortOnError );
      }

      adiosCheckpointHelpers::doSomethingForAFunctionOnAllPrimitives(
          io_, engine_, function, minLevel, maxLevel, adiosCheckpointHelpers::importVariables< func_t, value_t > );

      // might be more efficient to do this only after ALL functions were schedule for restoration!
      engine_.PerformGets();

      return true;
   }

 private:
   /// central ADIOS2 interface objects
   ///@{
   adios2::ADIOS  adios_;
   adios2::IO     io_;
   adios2::Engine engine_;
   ///@}

   /// meta information on the FE functions stored in the checkpoint
   std::vector< FunctionDescription > funcDescr_;

   /// auxilliary function to read attributes, will abort, if attribute is not found
   template < typename T >
   adios2::Attribute< T > readAttribute( const std::string& name )
   {
      adios2::Attribute< T > attribute = io_.InquireAttribute< T >( name );
      if ( !attribute )
      {
         WALBERLA_ABORT( "Attribute '" << name << "' seems to be missing from checkpoint!" );
      }
      return std::move( attribute );
   }

   /// auxilliary function to avoid code-duplication in c'tors
   void setupImporter( const std::string& filePath, const std::string& fileName )
   {
      // create the reader for the import
      io_ = adios_.DeclareIO( "AdiosCheckpointImport" );
      io_.SetEngine( "BP5" );

      // create the engine for the import
      //
      // at the moment we assume that the file only contains a single checkpoint
      // and do not deal with different steps
      std::string cpFileName = filePath + "/" + fileName;
      engine_ = io_.Open( cpFileName, adios2::Mode::ReadRandomAccess );

      // obtain FE function meta-info
      readFunctionDetailsFromCheckpoint();
   }

   void readFunctionDetailsFromCheckpoint()
   {
      adios2::Attribute< std::string > attrFuncNames      = readAttribute< std::string >( "FunctionNames" );
      adios2::Attribute< std::string > attrFuncKinds      = readAttribute< std::string >( "FunctionKinds" );
      adios2::Attribute< std::string > attrFuncValueTypes = readAttribute< std::string >( "FunctionValueTypes" );
      adios2::Attribute< uint_t >      attrFuncMinLevels  = readAttribute< uint_t >( "FunctionMinLevels" );
      adios2::Attribute< uint_t >      attrFuncMaxLevels  = readAttribute< uint_t >( "FunctionMaxLevels" );

      uint_t numFuncs = attrFuncNames.Data().size();
      WALBERLA_CHECK_EQUAL( numFuncs, attrFuncKinds.Data().size() );
      WALBERLA_CHECK_EQUAL( numFuncs, attrFuncValueTypes.Data().size() );
      WALBERLA_CHECK_EQUAL( numFuncs, attrFuncMinLevels.Data().size() );
      WALBERLA_CHECK_EQUAL( numFuncs, attrFuncMaxLevels.Data().size() );

      for ( uint_t idx = 0; idx < numFuncs; ++idx )
      {
         funcDescr_.emplace_back( FunctionDescription{ attrFuncNames.Data()[idx],
                                                       attrFuncKinds.Data()[idx],
                                                       attrFuncValueTypes.Data()[idx],
                                                       attrFuncMinLevels.Data()[idx],
                                                       attrFuncMaxLevels.Data()[idx] } );
      }
   }

   void writeFunctionDetails( std::ostream& stream, const std::string& prepend )
   {
      uint_t numFuncs = funcDescr_.size();

      uint_t maxLenName{ 0 };
      uint_t maxLenKind{ 0 };
      uint_t maxLenValueType{ 0 };

      for ( uint_t idx = 0; idx < numFuncs; ++idx )
      {
         maxLenName = maxLenName < funcDescr_[idx].name.length() ? funcDescr_[idx].name.length() : maxLenName;
         maxLenKind = maxLenKind < funcDescr_[idx].kind.length() ? funcDescr_[idx].kind.length() : maxLenKind;
         maxLenValueType =
             maxLenValueType < funcDescr_[idx].valueType.length() ? funcDescr_[idx].valueType.length() : maxLenValueType;
      }

      // check against length of "value_t"
      maxLenValueType = maxLenValueType < 7 ? 7 : maxLenValueType;

      uint_t      rowLen = 6u + 5u + ( maxLenName + 2u ) + 17u + ( maxLenKind + 2u ) + ( maxLenValueType + 2u );
      std::string hline( rowLen, '-' );
      std::string rowCol1( 5, '-' );
      std::string rowCol2( maxLenName + 2u, '-' );
      std::string rowCol3 = "--------|--------";
      std::string rowCol4( maxLenKind + 2u, '-' );
      std::string rowCol5( maxLenValueType + 2u, '-' );

      if ( numFuncs == 0 )
      {
         stream << prepend << "NO function inside checkpoint!\n";
         return;
      }
      else if ( numFuncs == 1 )
      {
         stream << prepend << "1 function inside checkpoint:\n";
      }
      else
      {
         stream << prepend << numFuncs << " functions inside checkpoint:\n";
      }
      stream << prepend << hline << '\n';
      stream << prepend << "| idx | " << std::setw( static_cast< int >( maxLenName ) ) << std::left << "name"
             << " | minLvl | maxLvl | " << std::setw( static_cast< int >( maxLenKind ) ) << std::left << "kind"
             << " | " << std::setw( static_cast< int >( maxLenValueType ) ) << "value_t"
             << " |\n";

      for ( uint_t idx = 0; idx < numFuncs; ++idx )
      {
         stream << prepend << "|" << rowCol1 << "|" << rowCol2 << "|" << rowCol3 << "|" << rowCol4 << "|" << rowCol5 << "|\n";

         // clang-format off
         stream << prepend << "| " << std::setw( 3 ) << std::right << idx << " | "
                           << std::setw( static_cast< int >( maxLenName ) ) << std::left << funcDescr_[idx].name << " | "
                           << std::right << std::setw( 6 ) << funcDescr_[idx].minLevel << " | "
                           << std::setw( 6 ) << funcDescr_[idx].maxLevel << " | "
                           << std::setw( static_cast< int >( maxLenKind ) ) << std::left << funcDescr_[idx].kind << " | "
                           << std::setw( static_cast< int >( maxLenValueType ) ) << funcDescr_[idx].valueType << " |\n";
         // clang-format on
      }
      stream << prepend << hline;
   }

   bool handleError( const std::string& functionName, const std::stringstream& msg, bool abortOnError )
   {
      if ( abortOnError )
      {
         WALBERLA_ABORT( "Restoring function '" << functionName << "' not possible; " << msg.str() );
      }
      WALBERLA_LOG_WARNING_ON_ROOT( "Restoring function '" << functionName << "' not possible; " << msg.str() );
      return false;
   }
};

} // namespace hyteg

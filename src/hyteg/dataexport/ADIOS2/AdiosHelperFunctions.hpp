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
#include <variant>

#include "hyteg/BuildInfo.hpp"
#include "hyteg/Git.hpp"
#include "hyteg/functions/FunctionMultiStore.hpp"

using walberla::real_t;
using walberla::uint_t;

namespace hyteg::adiosHelpers {

// Here define valid adios2 data types
typedef std::variant< int, long, std::string, real_t, uint_t, bool > adiostype_t;

// template specialization to check if a std::variant is consistent with adiostype_t
template < typename T, typename VARIANT_T >
struct isAdiosDataType;

template < typename T, typename... ALL_T >
struct isAdiosDataType< T, std::variant< ALL_T... > > : public std::disjunction< std::is_same< T, ALL_T >... >
{};

std::string generateVTKMetaInfo( const std::vector< std::string >& namesOfPointDataFunctions,
                                 const std::vector< std::string >& namesOfCellDataFunctions );

/// Name of scalar variable to be used in exporting time-step information
///
/// We use this variable to ensure consistency between the different places where this
/// variable name is needed, e.g in genewrateVTKMetaInfo and putTimeStepInfo. Note, though,
/// that currently TIME is the only allowed value.
extern const std::string nameOfTimeStepVariable;

/// Schedule information on current time step to be exported
inline void putTimeStepInfo( adios2::IO& io, adios2::Engine& engine, uint_t timestep )
{
   adios2::Variable< real_t > varTimeStep = io.InquireVariable< real_t >( nameOfTimeStepVariable );
   if ( !varTimeStep )
   {
      varTimeStep = io.DefineVariable< real_t >( nameOfTimeStepVariable );
   }
   engine.Put( varTimeStep, real_c( timestep ) );
}

/// Check whether the current process owns primitives with data to export (faces, or cell depending on dimension)
inline bool mpiProcessHasMacrosOfHighestDimension( const std::shared_ptr< PrimitiveStorage >& storage )
{
   bool weStoreRelevantData = false;
   if ( storage->hasGlobalCells() && storage->getCells().size() > 0 )
   {
      return true;
   }
   else if ( storage->getFaces().size() > 0 )
   {
      return true;
   }
   return false;
}

/// Write meta data on the software used for simulation into the given stream
inline std::ostream& printSoftwareMetaData( std::ostream& outStream )
{
   // clang-format off
   outStream << "Data generated with HyTeG (https://i10git.cs.fau.de/hyteg)\n"
             << "git branch         : " << gitBranch()     << '\n'
             << "SHA1 of last commit: " << gitSHA1()       << '\n'
             << "build type         : " << buildType()     << '\n'
             << "compiler           : " << compilerInfo()  << '\n'
             << "compiler flags     : " << compilerFlags() << '\n'
             << "mpi version        : " << mpiVersion()    << std::endl;
   // clang-format on

   return outStream;
}

/// Schedule meta data on the software used for simulation for export
inline void generateSoftwareMetaData( adios2::IO& io )
{
   adios2::Attribute< std::string > infoSoftware = io.InquireAttribute< std::string >( "Software" );
   if ( !infoSoftware )
   {
      std::stringstream oStream;
      printSoftwareMetaData( oStream );
      infoSoftware = io.DefineAttribute( "Software", oStream.str(), "", "" );
   }
}

template < typename T >
inline void writeAttribute( adios2::IO& io, std::string key, T const& val )
{
   adios2::Attribute< T > infoAttribute = io.InquireAttribute< T >( key );
   if ( !infoAttribute )
   {
      infoAttribute = io.DefineAttribute< T >( key, val );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Adios2 file: Attribute " << infoAttribute.Name() << "already written" );
   }
}

inline void writeAllAttributes( adios2::IO& io, std::map< std::string, adiostype_t > additionalAttributes )
{
   // integer datatype for ouptut
   using intData_t = ADIOS2_PARAVIEW_INT_TYPE;

   if ( additionalAttributes.empty() )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "No additional attributes/metadata written to adios2 output" );
      return;
   }

   for ( auto const& [key, val] : additionalAttributes )
   {
      std::visit(
          [&io, key]( auto&& arg ) {
             using T = std::decay_t< decltype( arg ) >;
             if ( isAdiosDataType< T, adiostype_t >::value )
             {
                // For Fortran compatibility, cast unsigned to signed integer
                if constexpr ( std::is_same_v< T, uint_t > )
                   writeAttribute( io, key, static_cast< intData_t >( arg ) );
                else if constexpr ( std::is_same_v< T, bool > )
                {
                   std::string sflag = arg ? "true" : "false";
                   writeAttribute( io, key, sflag );
                }
                else
                   writeAttribute( io, key, arg );
             }
             else
             {
                WALBERLA_LOG_WARNING_ON_ROOT( "The adios2 attribute " << key << " does not hold a valid datatype" );
             }
          },
          val );
   }
}

} // namespace hyteg::adiosHelpers

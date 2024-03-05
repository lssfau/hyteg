/*
 * Copyright (c) 2022 Berta Vilacis, Marcus Mohr.
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

#include <fstream>

#include "core/extern/json.hpp"

#include "terraneo/plates/types.hpp"

namespace terraneo {
namespace io {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

/// Auxilliary function for opening an input file

/// \param fileName  name (including potentially a path) of the input file
/// \return ifstream for the file
inline std::ifstream openFileForReading( std::string& fileName )
{
   std::ifstream infile( fileName, std::ios::in );
   if ( !infile.is_open() )
   {
      std::cerr << "Failed to open file '" << fileName << "' for reading!" << std::endl;
      std::abort();
   }
   return infile;
}

/// Imports a file in JSON format
inline nlohmann::json readJsonFile( std::string filename )
{
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Starting to read datafile '" << filename << "'" );
   std::ifstream  infile{ openFileForReading( filename ) };
   nlohmann::json inobj;
   infile >> inobj;
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Finished reading datafile '" << filename << "'" );
   return inobj;
}

/// Imports data from a text file with rotation information
inline std::vector< terraneo::plates::RotationInfo > readRotationsFile( std::string filename )
{
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Starting to read datafile '" << filename << "'" );
   std::ifstream file{ openFileForReading( filename ) };

   // Determine number of lines in file for reserving space
   std::string line;
   uint_t      numRotations{ 0 };
   while ( std::getline( file, line ) )
   {
      numRotations++;
   }
   file.clear();    // clear EOF
   file.seekg( 0 ); // rewind stream

   WALBERLA_LOG_INFO_ON_ROOT( "Found " << numRotations << " rotations in data-file" );
   std::vector< terraneo::plates::RotationInfo > rotations( numRotations, terraneo::plates::RotationInfo() );

   std::size_t charsUsed;
   uint_t      k{ 0 };

// this was primarily introduced for the regression testing during the
// refactoring of the original coding attempt; the accuracy of the data
// in the input file is smaller than even binary32 ;-)
#ifdef WALBERLA_DOUBLE_ACCURACY
#define PLATES_IO_STR_TO_FP std::stod
#else
#define PLATES_IO_STR_TO_FP std::stof
#endif

   while ( std::getline( file, line ) )
   {
      // position inside line-string
      std::size_t pos{ 0 };

      rotations[k].plateID = std::stoul( line, &charsUsed );

      pos += charsUsed;
      rotations[k].time = real_c( PLATES_IO_STR_TO_FP( line.substr( pos ), &charsUsed ) );

      pos += charsUsed;
      rotations[k].latitude = real_c( PLATES_IO_STR_TO_FP( line.substr( pos ), &charsUsed ) );

      pos += charsUsed;
      rotations[k].longitude = real_c( PLATES_IO_STR_TO_FP( line.substr( pos ), &charsUsed ) );

      pos += charsUsed;
      rotations[k].angle = real_c( PLATES_IO_STR_TO_FP( line.substr( pos ), &charsUsed ) );

      pos += charsUsed;
      rotations[k].conjugateID = std::stoul( line.substr( pos ), &charsUsed );

      // don't forget to increment index into vector
      ++k;
   }

   // Safety check
   if ( k != numRotations )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " Imported " << k << " rotations,\n"
                                              << " but should have been " << numRotations );
      WALBERLA_ABORT( "Data import error" );
   }

   WALBERLA_LOG_PROGRESS_ON_ROOT( "Finished reading datafile '" << filename << "'" );

   return rotations;
}

/**
 * @brief Read profile data into a 2D data_vector.
 *
 * This function reads profile data from a file in either JSON, CSV, or TXT format.
 * If the file has a JSON extension, it expects the data to be in the form of a JSON object
 * with two arrays: "Radius (m)" and the data array (i.e Viscosity or Temperature). 
 * If the file has a CSV or TXT extension, it expects the data to be in two columns: the first column representing the radius 
 * and the second column representing the physical data.
 *
 * @param filename The name of the file to read the profile data from.
 * @param data_vector A 2D array where the radius and profile data will be stored.
 * @param num_columns The expected number of columns in the data. This is used to check if the columns have the same size.
 * @return True if the profile data was successfully read from the file, false otherwise.
 */
inline bool readProfileData( const std::string& filename, std::vector< std::vector< real_t > >& data_vector )
{
   std::ifstream file( filename );

   if ( !file.is_open() )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Failed to open file '" << filename << "' for reading!" );
      return false;
   }

   else
   {
      // Check if the file has a .json extension
      if ( filename.substr( filename.find_last_of( "." ) + 1 ) == "json" )
      {
         // Open the JSON file
         nlohmann::json jsonData = readJsonFile( filename );

         // Check if the JSON file contains an object
         if ( !jsonData.is_object() )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "JSON file does not contain an object." );
            return false;
         }

         // Iterate over the keys in the JSON object
         std::vector< std::string > keys;

         for ( auto it = jsonData.begin(); it != jsonData.end(); ++it )
         {
            keys.push_back( it.key() );
         }

         // Check if the keys are empty (i.e. no string keys found in JSON object)
         if ( keys.empty() )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "No keys found in JSON object." );
            return false;
         }

         // Get the "Radius (m)" array from the JSON file
         auto radiusArray = jsonData[keys[0]];

         // Get the data array from the JSON file
         auto dataArray = jsonData[keys[1]];

         // Check if both data arrays are valid
         if ( !radiusArray.is_array() || !dataArray.is_array() )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "Invalid arrays in JSON object." );
            return false;
         }

         // Check if the arrays have the same size (same number of elements)
         if ( radiusArray.size() != dataArray.size() )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "Arrays in JSON object have different sizes." );
            return false;
         }

         // Iterate over the elements in the arrays and store them in the data_vector
         for ( std::size_t i = 0; i < radiusArray.size(); ++i )
         {
            real_t radius = radiusArray[i].get< real_t >();
            real_t data   = dataArray[i].get< real_t >();

            data_vector.push_back( { radius, data } );
         }
         return true;
      }
      else
      {
         // Must be then either a .csv or .txt file
         std::string line;
         std::size_t numColumns = 0;
         while ( std::getline( file, line ) )
         {
            std::istringstream    sstring( line );
            std::vector< real_t > row;
            real_t                radius, data;
            if ( sstring >> radius >> data )
            {
               row.push_back( radius );
               row.push_back( data );
               data_vector.push_back( row );
               numColumns = row.size();
            }
         }

         // Check if the two data columns have the same size
         for ( const auto& row : data_vector )
         {
            if ( row.size() != numColumns )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "Columns in the file have different sizes." );
               return false;
            }
         }
      }
      file.close();
      return true;
   }
}
} // namespace io
} // namespace terraneo

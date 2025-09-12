/*
 * Copyright (c) 2024-2025 Andreas Burkhart.
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

#include "core/DataTypes.h"
#include "core/config/Config.h"

#include "../Parameters/NondimensionalisationParameters.hpp"
#include "DataTools.hpp"

using walberla::real_t;

namespace MantleConvection {

csvData loadCSV( std::string filename )
{
   csvData data;

   std::ifstream fileIn;
   fileIn.open( filename.c_str(), std::ios::in );
   if ( !fileIn.is_open() )
   {
      WALBERLA_LOG_WARNING( "File " << filename << " not found!" );
      return data;
   };

   bool foundRowsRange          = false;
   bool foundColsClassification = false;

   std::string line;
   std::string input;

   if ( fileIn.good() )
   {
      // read header
      while ( ( ( foundRowsRange == false ) || ( foundColsClassification == false ) ) && std::getline( fileIn, line ) )
      {
         if ( line.at( 0 ) != '#' )
         {
            std::istringstream linestream( line );

            if ( std::getline( linestream, input, ',' ) )
            {
               if ( input == "R" )
               {
                  // rowsRange
                  if ( std::getline( linestream, input, ',' ) )
                  {
                     foundRowsRange = true;
                     // convert string to uint_t
                     std::istringstream convertstream( input );
                     convertstream >> data.rowsRange;
                  }
               }
               else if ( input == "C" )
               {
                  // colsClassification
                  if ( std::getline( linestream, input, ',' ) )
                  {
                     foundColsClassification = true;
                     // convert string to uint_t
                     std::istringstream convertstream( input );
                     convertstream >> data.colsClassification;
                  }
               }
            }
         }
      }

      if ( foundRowsRange == false || foundColsClassification == false )
      {
         WALBERLA_LOG_WARNING( "loadCSV could not find a header in " << filename << "!" );
         return data;
      }

      std::reference_wrapper< std::vector< real_t > > vecRef( data.values.back() );
      uint_t                                          linecounter           = 0;
      uint_t                                          classificationCounter = 0;
      // iterate over all lines
      while ( std::getline( fileIn, line ) )
      {
         // ignore comments
         if ( line.at( 0 ) != '#' )
         {
            std::istringstream linestream( line );

            if ( linecounter < data.rowsRange )
            {
               data.range.push_back( std::vector< real_t >() );
               vecRef = data.range.back();
            }
            else
            {
               data.values.push_back( std::vector< real_t >() );
               vecRef = data.values.back();
            }

            data.classification.push_back( std::vector< MantleConvection::NondimensionalisationType >() );
            classificationCounter = 0;

            // read csv
            while ( std::getline( linestream, input, ',' ) )
            {
               if ( classificationCounter < data.colsClassification )
               {
                  // throws error if type is unknown
                  data.classification.back().push_back( MantleConvection::NondimensionalisationTypeMap.at( input ) );
               }
               else
               {
                  // convert string to real_t
                  std::istringstream convertstream( input );
                  real_t             value;
                  convertstream >> value;
                  vecRef.get().push_back( value );
               }

               classificationCounter++;
            }
         }
         linecounter++;
      }
   }

   return data;
}

void nondimensionaliseCSV( walberla::Config::BlockHandle& parameters,
                           csvData&                       data,
                           bool                           nondimRange  = true,
                           bool                           nondimValues = true )
{
   if ( data.colsClassification >= 1 )
   {
      NondimensionalisationParameters temporaryND( parameters );

      if ( nondimRange )
      {
         for ( uint_t r = 0; r < data.range.size(); r++ )
         {
            std::vector< real_t >&                       vec  = data.range.at( r );
            MantleConvection::NondimensionalisationType& type = data.classification.at( r ).front();

            for ( uint_t c = 0; c < vec.size(); c++ )
            {
               nondimensionaliseValue( temporaryND, type, vec.at( c ) );
            }
         }
      }

      if ( nondimValues )
      {
         for ( uint_t r = 0; r < data.values.size(); r++ )
         {
            std::vector< real_t >&                       vec  = data.values.at( r );
            MantleConvection::NondimensionalisationType& type = data.classification.at( r + data.rowsRange ).front();

            for ( uint_t c = 0; c < vec.size(); c++ )
            {
               nondimensionaliseValue( temporaryND, type, vec.at( c ) );
            }
         }
      }
   }
}

} // namespace MantleConvection
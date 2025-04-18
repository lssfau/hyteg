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

#include <algorithm>

#include "core/extern/json.hpp"

#include "terraneo/dataimport/FileIO.hpp"
#include "terraneo/helpers/conversions.hpp"
#include "terraneo/helpers/typeAliases.hpp"
#include "terraneo/plates/functionsForRotations.hpp"
#include "terraneo/plates/types.hpp"
#include "terraneo/plates/utilities.hpp"

namespace terraneo {
namespace plates {

/// Class for managing plate topology information
class PlateStorage
{
 public:
   /// Information describing a single plate at a certain age
   struct PlateInfo
   {
      /// ID of the plate
      uint_t id{0};

      /// Nodes describing the boundary of the plate
      Polygon boundary;

      /// Textual name of plate
      std::string name;
   };

   using plateVec_t       = std::vector< PlateInfo >;
   using ageToPlatesMap_t = std::map< std::string, plateVec_t >;

   template < typename ImportStrategy >
   PlateStorage( std::string nameOfPlateTopologiesFile, real_t sphereRadius, ImportStrategy readJsonfile )
   : srcFile_( nameOfPlateTopologiesFile )
   {
      // import topology data
      nlohmann::json rootNode = readJsonfile( nameOfPlateTopologiesFile )["ages"];

      // convert required data into our internal format
      extractPlateInfo( rootNode );

      // for testing
      printStatistics();
   }

   /// Report statistics on what an object of this class stores
   void printStatistics()
   {
      uint_t nPlates{0};
      for ( auto entry : ageToPlatesMap_ )
      {
         nPlates += entry.second.size();
      }
      WALBERLA_LOG_INFO_ON_ROOT( "PlateStorage object:\n"
                                 << " - stores plates for " << ageToPlatesMap_.size() << " age stages\n"
                                 << " - stores a total of " << nPlates << " plates\n"
                                 << " - age stages range from " << *listOfPlateStages_.begin() << " to "
                                 << listOfPlateStages_.back() << " Ma\n"
                                 << " - data was obtained from file = '" << srcFile_ << "'" );
   }

   plateVec_t& getPlatesForStage( real_t age )
   {
      auto iter = ageToPlatesMap_.find( ageToKeyStr( age ) );
      if ( iter == ageToPlatesMap_.end() )
      {
         std::cerr << "No plates found for " << ageToKeyStr( age ) << std::endl;
         std::abort();
      }

      return iter->second;
   }

   const plateVec_t& getPlatesForStage( real_t age ) const
   {
      auto iter = ageToPlatesMap_.find( ageToKeyStr( age ) );
      if ( iter == ageToPlatesMap_.end() )
      {
         std::cerr << "No plates found for " << ageToKeyStr( age ) << std::endl;
         std::abort();
      }

      return iter->second;
   }

   const std::vector< real_t >& getListOfPlateStages() const { return listOfPlateStages_; }

 private:
   /// method for data preparation
   ///
   /// this method will convert the imported data into a format more suitable for
   /// using it in the computations; it will convert the plate polygons into
   /// cartesian coordinates, compute the barycenter of each plate and group
   /// plates into vectors depending on their age stage
   void extractPlateInfo( const nlohmann::json& rootNode )
   {
      // we need to fill the list of age stages and do a reservation here
      listOfPlateStages_.clear();
      listOfPlateStages_.reserve( rootNode.size() );

      for ( uint_t idx = 0; idx < rootNode.size(); ++idx )
      {
         // prepare key for map entry and vector to hold plates for this age stage
         std::string ageName = rootNode[idx]["name"];
         uint_t      nPlates = rootNode[idx]["features"].size();
         ageToPlatesMap_.emplace( ageName, plateVec_t( nPlates ) );

         plateVec_t& plates = ageToPlatesMap_[ageName];

         // convert ageName to numeric age value and insert it
         listOfPlateStages_.push_back( keyStrToAge( ageName ) );

         // extract plates for this age stage and put into vector of plates
         for ( uint_t k = 0; k < nPlates; ++k )
         {
            // get ID and name
            plates[k].id   = rootNode[idx]["features"][k]["properties"]["PLATEID1"];
            plates[k].name = static_cast<std::string>(rootNode[idx]["features"][k]["properties"]["NAME"]);

            // extract coordinates of polygon (order: longitude, latitude)
            const nlohmann::json coordinates = rootNode[idx]["features"][k]["geometry"]["coordinates"][0];
            for ( auto& element : coordinates )
            {
               plates[k].boundary.push_back( {element[0], element[1], real_c( 0 )} );
            }
         }
      }

      // sort the list of plate stages (should be sorted in data-file, but hey,
      // better safe than sorry
      std::sort( listOfPlateStages_.begin(), listOfPlateStages_.end() );
   };

   /// name of datafile from which object obtained information
   std::string srcFile_;

   /// map to allow retrieving all plates of a certain age stage by giving that age
   ageToPlatesMap_t ageToPlatesMap_;

   /// a vector containing a sorted list of the ages of the plate stages found in data-file
   std::vector< real_t > listOfPlateStages_;
};

} // namespace plates
} // namespace terraneo

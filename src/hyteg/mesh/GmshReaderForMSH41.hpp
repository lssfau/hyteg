/*
 * Copyright (c) 2024 Marcus Mohr.
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

#include <array>
#include <vector>

#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"

#include "hyteg/mesh/MeshInfo.hpp"

namespace hyteg {

class GmshReaderForMSH41
{
 public:
   GmshReaderForMSH41( const std::string& meshFileName, bool importPhysicalTags )
   : meshFileName_( meshFileName )
   , importPhysicalTags_( importPhysicalTags ){};

   MeshInfo readMesh()
   {
      std::ifstream meshFile = openFileAndCheckFormat();

      std::vector< std::string > sections = getSectionsPresentInFile( meshFile );
      analyseSectionList( sections );

      MeshInfo meshInfo;

      // We are currently not making use of the physical names, if they exist.
      // But potentially this will change in the future.
      if ( importPhysicalTags_ && std::find( sections.begin(), sections.end(), "PhysicalNames" ) != sections.end() )
      {
         readSectionPhysicalNames( meshFile );
      }

      std::map< marker, uint_t > entityToPhysicalTag;
      if ( importPhysicalTags_ && std::find( sections.begin(), sections.end(), "Entities" ) != sections.end() )
      {
         entityToPhysicalTag = readSectionEntities( meshFile );
      }

      readSectionNodes( meshFile, meshInfo, entityToPhysicalTag );
      readSectionElements( meshFile, meshInfo, entityToPhysicalTag );

      meshFile.close();

      WALBERLA_LOG_PROGRESS_ON_ROOT( "----------------------------------------------" );

      return meshInfo;
   }

   /// Gmsh marks stuff by using a pair of geometric dimension and tag
   using marker = std::pair< uint_t, uint_t >;

 private:
   /// name of input file
   const std::string& meshFileName_;

   /// If this flag is false the reader will only import the node and connectivity information from the file
   const bool importPhysicalTags_;

   /// the MSH can contain various kinds of sections, some of which we need, others we cannot work with
   typedef enum
   {
      REQUIRED,                   //< must be present for the reader to work
      REQUIRED_FOR_PHYSICAL_TAGS, //< required, if the user wants us to import physical tags
      OPTIONAL,                   //< if present, we can work with these data
      IGNORABLE,                  //< if present, we can simply ignore it
      PROBABLY_IGNORABLE,         //< if present, we assume we can ignore it
      CRITICAL                    //< when present indicates that the file is special
   } sectionClassification;

   /// classify sections that can be present in MSH 4.1 format
   static inline const std::map< std::string, sectionClassification > sectionType_ = {
       { "MeshFormat", REQUIRED },
       { "PhysicalNames", OPTIONAL },
       { "Entities", REQUIRED_FOR_PHYSICAL_TAGS },
       { "PartitionedEntities", PROBABLY_IGNORABLE },
       { "Nodes", REQUIRED },
       { "Elements", REQUIRED },
       { "Periodic", CRITICAL },
       { "GhostElements", PROBABLY_IGNORABLE },
       { "Parametrizations", CRITICAL },
       { "NodeData", IGNORABLE },
       { "ElementData", IGNORABLE },
       { "ElementNodeData", IGNORABLE },
       { "InterpolationScheme", CRITICAL },
       { "Comments", IGNORABLE } };

   /// Open the mesh-file and verify that its format is MSH4.1
   std::ifstream openFileAndCheckFormat() const
   {
      std::ifstream meshFile;
      meshFile.open( meshFileName_.c_str() );

      WALBERLA_CHECK( !!meshFile, "[Mesh] Error opening file: " << meshFileName_ );
      WALBERLA_LOG_INFO_ON_ROOT( "Reading data from file '" << meshFileName_ << "'" );

      std::string token;
      meshFile >> token; // $MeshFormat

      WALBERLA_CHECK_EQUAL( token, "$MeshFormat", "[Mesh] Missing: $MeshFormat" );

      meshFile >> token; // version

      WALBERLA_CHECK_EQUAL( token, "4.1", "[Mesh] Meshfile version should be 4.1" );

      return std::move( meshFile );
   }

   /// Determine the sections actually present in the MSH file
   std::vector< std::string > getSectionsPresentInFile( std::ifstream& file ) const;

   /// run an analysis on the sections present to see whether our reader can work with the MSH file
   void analyseSectionList( const std::vector< std::string >& sectionList ) const;

   /// locate a specific data section inside the MSH file
   void findSection( std::ifstream& file, const std::string& sectionName ) const
   {
      // rewind the stream
      file.clear();
      file.seekg( 0 );

      // look for section
      std::string line;
      while ( std::getline( file, line ) )
      {
         if ( line[0] == '$' && line.substr( 1, line.size() ) == sectionName )
         {
            return;
         }
      }
      WALBERLA_ABORT( "Could not find section '" << sectionName << "' in file '" << meshFileName_ << "'" );
   }

   /// import information on physical names defined in the MSH file
   std::map< marker, std::string > readSectionPhysicalNames( std::ifstream& file ) const;

   /// import information on entities defined in the MSH file
   [[nodiscard]] std::map< marker, uint_t > readSectionEntities( std::ifstream& file ) const;

   /// import information on the nodes defined in the MSH file
   void readSectionNodes( std::ifstream& file, MeshInfo& meshInfo, const std::map< marker, uint_t >& entityToPhysicalTag ) const;

   /// import information on the elements defined in the MSH file
   void readSectionElements( std::ifstream&                    file,
                             MeshInfo&                         meshInfo,
                             const std::map< marker, uint_t >& entityToPhysicalTag ) const;
};

} // namespace hyteg

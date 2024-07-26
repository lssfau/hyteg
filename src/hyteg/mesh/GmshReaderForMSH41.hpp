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
   GmshReaderForMSH41( const std::string& meshFileName, bool beVerbose = true )
   : meshFileName_( meshFileName )
   , beVerbose_( beVerbose ){};

   MeshInfo readMesh()
   {
      // for testing purposes
      // walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );

      std::ifstream meshFile = openFileAndCheckFormat();

      std::vector< std::string > sections = getSectionsPresentInFile( meshFile );
      analyseSectionList( sections );

      MeshInfo meshInfo;

      readSectionPhysicalNames( meshFile );
      readSectionNodes( meshFile, meshInfo );
      readSectionElements( meshFile, meshInfo );

      meshFile.close();

      WALBERLA_LOG_PROGRESS_ON_ROOT( "----------------------------------------------" );

      return meshInfo;
   }

 private:
   /// if set to true the object generates verbose reports on its actions
   const bool beVerbose_;

   /// name of input file
   const std::string& meshFileName_;

   /// the MSH can contain various kinds of sections, some of which we need, others we cannot work with
   typedef enum
   {
      REQUIRED,           //< must be present for the reader to work
      OPTIONAL,           //< if present, we can work with these data
      IGNORABLE,          //< if present, we can simply ignore it
      PROBABLY_IGNORABLE, //< if present, we assume we can ignore it
      CRITICAL            //< when present indicates that the file is special
   } sectionClassification;

   /// classify sections that can be present in MSH 4.1 format
   static inline const std::map< std::string, sectionClassification > sectionType_ = {
       { "MeshFormat", REQUIRED },
       { "PhysicalNames", OPTIONAL },
       { "Entities", OPTIONAL },
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

      if ( beVerbose_ )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Reading data from file '" << meshFileName_ << "'" );
      }

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
   void readSectionPhysicalNames( std::ifstream& file ) const;

   /// import information on entities defined in the MSH file
   void readSectionEntities( std::ifstream& file ) const { WALBERLA_ABORT( "readSectionEntities() not implemented, yet!" ); };

   /// import information on the nodes defined in the MSH file
   void readSectionNodes( std::ifstream& file, MeshInfo& meshInfo ) const;

   /// import information on the elements defined in the MSH file
   void readSectionElements( std::ifstream& file, MeshInfo& meshInfo ) const;
};

} // namespace hyteg

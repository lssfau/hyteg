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

#include "hyteg/mesh/GmshReaderForMSH41.hpp"

namespace hyteg {

// Determine the sections actually present in the MSH file
std::vector< std::string > GmshReaderForMSH41::getSectionsPresentInFile( std::ifstream& file ) const
{
   std::string                line;
   std::vector< std::string > sections;
   sections.reserve( 14 );

   // rewind the stream
   file.clear();
   file.seekg( 0 );

   while ( std::getline( file, line ) )
   {
      if ( line[0] == '$' )
      {
         if ( line.substr( 1, 3 ) != "End" )
         {
            sections.emplace_back( line.substr( 1, line.size() ) );
         }
      }
   }

   // report
   WALBERLA_LOG_PROGRESS_ON_ROOT( "----------------------------------------------" );
   WALBERLA_LOG_PROGRESS_ON_ROOT( "File contains the following sections:" );
   for ( const auto& section : sections )
   {
      WALBERLA_LOG_PROGRESS_ON_ROOT( " * " << section );
   }

   return sections;
}

// run an analysis on the sections present to see whether our reader can work with the MSH file
void GmshReaderForMSH41::analyseSectionList( const std::vector< std::string >& sectionList ) const
{
   bool   weNeedToAbort       = false;
   uint_t numRequiredSections = 0;

   for ( const auto& section : sectionList )
   {
      if ( sectionType_.count( section ) == 0 )
      {
         WALBERLA_ABORT( "Section '" << section << "' not classifiable! It this really defined for MSH4.1?" );
      }

      switch ( sectionType_.at( section ) )
      {
      case REQUIRED:
         WALBERLA_LOG_PROGRESS_ON_ROOT( "Required section '" << section << "' is present" );
         numRequiredSections++;
         break;
      case OPTIONAL:
         WALBERLA_LOG_PROGRESS_ON_ROOT( "Optional section '" << section << "' is present" );
         break;
      case IGNORABLE:
         WALBERLA_LOG_WARNING_ON_ROOT( "Going to ignore section '" << section << "'" );
         break;
      case PROBABLY_IGNORABLE:
         WALBERLA_LOG_WARNING_ON_ROOT( "Don't know how to handle section '"
                                       << section << "'. Assuming that it can be  ignored. Keep your fingers crossed" );
         break;
      case CRITICAL:
         WALBERLA_LOG_WARNING_ON_ROOT( "Detected section '" << section << "'. Processing this mesh is outside my capabilities!" );
         weNeedToAbort = true;
      }
   }

   if ( weNeedToAbort )
   {
      WALBERLA_ABORT( "Sorry but the reader is not flexible enough to handle the data in file '" << meshFileName_ << "'" );
   }
   else if ( numRequiredSections != 3u )
   {
      WALBERLA_ABORT( "Only " << numRequiredSections << "/3 required sections present!\n"
                              << "We minimally need 'MeshFormat', 'Nodes', and 'Elements'." );
   }
};

// import information on physical names defined in the MSH file
void GmshReaderForMSH41::readSectionPhysicalNames( std::ifstream& file ) const
{
   WALBERLA_LOG_PROGRESS_ON_ROOT( "----------------------------------------------" );
   WALBERLA_LOG_PROGRESS_ON_ROOT( "readSectionPhysicalNames():" );

   // locate section within file
   findSection( file, "PhysicalNames" );

   // read section header
   uint_t numPhysicalNames = 0;
   file >> numPhysicalNames;
   WALBERLA_LOG_PROGRESS_ON_ROOT( "-> File contains " << numPhysicalNames << " 'physical name(s)'" );

   for ( uint_t idx = 0; idx < numPhysicalNames; ++idx )
   {
      uint_t      dimension   = 0;
      uint_t      physicalTag = 0;
      std::string name;

      file >> dimension >> physicalTag >> name;
      // name is allowed to contain blanks, so look for closing "
      while ( name[name.size() - 1] != '"' )
      {
         std::string nextPart;
         file >> nextPart;
         name.append( " " );
         name.append( nextPart );
      }

      WALBERLA_LOG_PROGRESS_ON_ROOT( "" << name << ": (dimension, tag) = (" << dimension << ", " << physicalTag << ")" );
   }
}

// import information on the nodes defined in the MSH file
void GmshReaderForMSH41::readSectionNodes( std::ifstream& file, MeshInfo& meshInfo ) const
{
   WALBERLA_LOG_PROGRESS_ON_ROOT( "----------------------------------------------" );
   WALBERLA_LOG_PROGRESS_ON_ROOT( "readSectionNodes():" );

   // locate section within file
   findSection( file, "Nodes" );

   // read section info
   uint_t numBlocks = 0;
   uint_t numNodes  = 0;
   uint_t minTag    = 0;
   uint_t maxTag    = 0;
   file >> numBlocks >> numNodes >> minTag >> maxTag;

   WALBERLA_LOG_PROGRESS_ON_ROOT( "-> numBlocks = " << numBlocks << ", numNodes = " << numNodes << ", minTag = " << minTag
                                                    << ", maxTag = " << maxTag );

   uint_t nodeCount = 0;

   for ( uint_t block = 0; block < numBlocks; ++block )
   {
      // read block header
      uint_t entityDim       = 0;
      uint_t entityTag       = 0;
      uint_t parametric      = 0;
      uint_t numNodesInBlock = 0;

      file >> entityDim >> entityTag >> parametric >> numNodesInBlock;

      nodeCount += numNodesInBlock;

      if ( parametric != 0 )
      {
         WALBERLA_ABORT( "Cannot handle parametric nodes!" );
      }

      WALBERLA_LOG_PROGRESS_ON_ROOT( "" << numNodesInBlock << " node(s) for entity with (dim, tag) = (" << entityDim << ", "
                                        << entityTag << ")" );

      // read node tags
      std::vector< uint_t > nodeTags;
      nodeTags.reserve( numNodesInBlock );
      for ( uint_t idx = 0; idx < numNodesInBlock; ++idx )
      {
         uint_t tag;
         file >> tag;
         nodeTags.push_back( tag );
      }

      // read node coordinates
      for ( uint_t idx = 0; idx < numNodesInBlock; ++idx )
      {
         real_t xCoord;
         real_t yCoord;
         real_t zCoord;
         file >> xCoord >> yCoord >> zCoord;

         MeshInfo::IDType id    = nodeTags[idx];
         meshInfo.vertices_[id] = MeshInfo::Vertex( id, Point3D( xCoord, yCoord, zCoord ), 0 );
      }
   }

   WALBERLA_CHECK_EQUAL( numNodes, nodeCount, "Inconsistency w.r.t. number of nodes promised and present!" );
}

// import information on the elements defined in the MSH file
void GmshReaderForMSH41::readSectionElements( std::ifstream& file, MeshInfo& meshInfo ) const
{
   WALBERLA_LOG_PROGRESS_ON_ROOT( "----------------------------------------------" );
   WALBERLA_LOG_PROGRESS_ON_ROOT( "readSectionElements():" );

   // locate section within file
   findSection( file, "Elements" );

   // read section header
   uint_t numEntityBlocks = 0;
   uint_t numElements     = 0;
   uint_t minElementTag   = 0;
   uint_t maxElementTag   = 0;

   file >> numEntityBlocks >> numElements >> minElementTag >> maxElementTag;

   WALBERLA_LOG_PROGRESS_ON_ROOT( "-> numEntityBlocks = " << numEntityBlocks << ", numElements = " << numElements
                                                          << ", minElementTag = " << minElementTag
                                                          << ", maxElementTag = " << maxElementTag );

   uint_t elementCount = 0;

   // Gmsh element types we can work with
   const uint_t twoNodeLine         = 1u;
   const uint_t threeNodeTriangle   = 2u;
   const uint_t fourNodeTetrahedron = 4u;

   // We perform a two-pass approach:
   // First we parse the section and store the information on the elements
   // in respective containers. Afterwards we add it to the MeshInfo object
   // going through the elements in increasing dimension.
   MeshInfo::EdgeContainer parsedEdges;
   MeshInfo::FaceContainer parsedFaces;
   MeshInfo::CellContainer parsedCells;

   for ( uint_t block = 0; block < numEntityBlocks; ++block )
   {
      // read block header
      uint_t entityDim          = 0;
      uint_t entityTag          = 0;
      uint_t elementType        = 0;
      uint_t numElementsInBlock = 0;

      file >> entityDim >> entityTag >> elementType >> numElementsInBlock;

      elementCount += numElementsInBlock;

      WALBERLA_LOG_PROGRESS_ON_ROOT( "" << numElementsInBlock << " element(s) for entity with (dim, tag, type) = (" << entityDim
                                        << ", " << entityTag << ", " << elementType << ")" );

      switch ( elementType )
      {
      case twoNodeLine: {
         uint_t                            tag = 0;
         std::array< MeshInfo::IDType, 2 > edgeNodes;

         for ( uint_t idx = 0; idx < numElementsInBlock; ++idx )
         {
            file >> tag >> edgeNodes[0] >> edgeNodes[1];
            parsedEdges[edgeNodes] = MeshInfo::Edge( edgeNodes, 0 );
         }
         break;
      }
      case threeNodeTriangle: {
         uint_t                          tag = 0;
         std::vector< MeshInfo::IDType > triangleNodes( 3 );
         for ( uint_t idx = 0; idx < numElementsInBlock; ++idx )
         {
            file >> tag >> triangleNodes[0] >> triangleNodes[1] >> triangleNodes[2];
            parsedFaces[triangleNodes] = MeshInfo::Face( triangleNodes, 0 );
         }
         break;
      }
      case fourNodeTetrahedron: {
         uint_t                          tag = 0;
         std::vector< MeshInfo::IDType > tetrahedronNodes( 4 );
         for ( uint_t idx = 0; idx < numElementsInBlock; ++idx )
         {
            file >> tag >> tetrahedronNodes[0] >> tetrahedronNodes[1] >> tetrahedronNodes[2] >> tetrahedronNodes[3];
            parsedCells[tetrahedronNodes] = MeshInfo::Cell( tetrahedronNodes, 0 );
         }
         break;
      }
      default: {
         WALBERLA_ABORT( "Detected unsupported element type: " << elementType );
      }
      }
   }

   // second pass: insert primitives into object
   meshInfo.processPrimitivesFromGmshFile( parsedEdges, parsedFaces, parsedCells );
}

} // namespace hyteg

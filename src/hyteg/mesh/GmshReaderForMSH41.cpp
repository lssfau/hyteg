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

using marker = GmshReaderForMSH41::marker;

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
      WALBERLA_UNUSED( section ); // if loglevel is less detailed than progress
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
      case sectionClassification::REQUIRED:
         WALBERLA_LOG_PROGRESS_ON_ROOT( "Required section '" << section << "' is present" );
         numRequiredSections++;
         break;
      case sectionClassification::REQUIRED_FOR_PHYSICAL_TAGS:
         if ( importPhysicalTags_ )
         {
            WALBERLA_LOG_PROGRESS_ON_ROOT( "Required section '" << section << "' is present" );
            numRequiredSections++;
         }
         else
         {
            WALBERLA_LOG_PROGRESS_ON_ROOT( "Optional section '" << section << "' is present" );
         }
         break;
      case sectionClassification::OPTIONAL:
         WALBERLA_LOG_PROGRESS_ON_ROOT( "Optional section '" << section << "' is present" );
         break;
      case sectionClassification::IGNORABLE:
         WALBERLA_LOG_WARNING_ON_ROOT( "Going to ignore section '" << section << "'" );
         break;
      case sectionClassification::PROBABLY_IGNORABLE:
         WALBERLA_LOG_WARNING_ON_ROOT( "Don't know how to handle section '"
                                       << section << "'. Assuming that it can be ignored. Keep your fingers crossed" );
         break;
      case sectionClassification::CRITICAL:
         WALBERLA_LOG_WARNING_ON_ROOT( "Detected section '" << section << "'. Processing this mesh is outside my capabilities!" );
         weNeedToAbort = true;
      }
   }

   if ( weNeedToAbort )
   {
      WALBERLA_ABORT( "Sorry but the reader is not flexible enough to handle the data in file '" << meshFileName_ << "'" );
   }
   else if ( numRequiredSections != 3u + ( importPhysicalTags_ ? 1u : 0u ) )
   {
      if ( !importPhysicalTags_ )
      {
         WALBERLA_ABORT( "Only " << numRequiredSections << "/3 required sections present!\n"
                                 << "We minimally need 'MeshFormat', 'Nodes', and 'Elements'." );
      }
      else
      {
         WALBERLA_ABORT( "Only " << numRequiredSections << "/4 required sections present!\n"
                                 << "We minimally need 'MeshFormat', 'Entities', 'Nodes', and 'Elements'." );
      }
   }
}

// import information on physical names defined in the MSH file
std::map< marker, std::string > GmshReaderForMSH41::readSectionPhysicalNames( std::ifstream& file ) const
{
   WALBERLA_LOG_PROGRESS_ON_ROOT( "----------------------------------------------" );
   WALBERLA_LOG_PROGRESS_ON_ROOT( "readSectionPhysicalNames():" );

   // locate section within file
   findSection( file, "PhysicalNames" );

   // going to store the (dimension,tag) <-> name association in this map
   std::map< marker, std::string > markerToName;

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

      marker onMyMark{ dimension, physicalTag };
      markerToName[onMyMark] = name;
   }

   return markerToName;
}

// import information on physical names defined in the MSH file
std::map< marker, uint_t > GmshReaderForMSH41::readSectionEntities( std::ifstream& file ) const
{
   WALBERLA_LOG_PROGRESS_ON_ROOT( "----------------------------------------------" );
   WALBERLA_LOG_PROGRESS_ON_ROOT( "readSectionEntities():" );

   // locate section within file
   findSection( file, "Entities" );

   // read section header
   uint_t numNodes;
   uint_t numCurves;
   uint_t numSurfaces;
   uint_t numVolumes;

   file >> numNodes >> numCurves >> numSurfaces >> numVolumes;

   uint_t numEntities = numNodes + numCurves + numSurfaces + numVolumes;
   WALBERLA_UNUSED( numEntities );
   WALBERLA_LOG_PROGRESS_ON_ROOT( "-> File contains " << numEntities << " 'entities'" );
   WALBERLA_LOG_PROGRESS_ON_ROOT( "   - " << numNodes << " node(s)" );
   WALBERLA_LOG_PROGRESS_ON_ROOT( "   - " << numCurves << " curve(s)" );
   WALBERLA_LOG_PROGRESS_ON_ROOT( "   - " << numSurfaces << " surface(s)" );
   WALBERLA_LOG_PROGRESS_ON_ROOT( "   - " << numVolumes << " volume(s)" );

   // we will store the entity to physical tags association in this map
   //
   // Note: We do currently not support more than one physical tag per entity!
   //       In case there are more than one physical tags, we will read only
   //       the first one and discard the remaining ones!
   std::map< marker, uint_t > entityToPhysicalTag;

   // count of entities without any physical tag
   uint_t numNonTagged{ 0u };

   // read info on node entities (dim=0)
   for ( uint_t idx = 0; idx < numNodes; ++idx )
   {
      uint_t pointTag;
      double x = -0.1;
      double y = -0.1;
      double z = -0.1;
      uint_t numPhysicalTags;
      uint_t physicalTag;
      file >> pointTag >> x >> y >> z >> numPhysicalTags;

      WALBERLA_LOG_PROGRESS_ON_ROOT( "Processing node with pointTag = " << pointTag
                                                                        << ", numPhysicalTags = " << numPhysicalTags );
      if ( numPhysicalTags >= 1u )
      {
         file >> physicalTag;
         marker entityMarker{ 0u, pointTag };
         entityToPhysicalTag[entityMarker] = physicalTag;

         if ( numPhysicalTags > 1u )
         {
            WALBERLA_LOG_WARNING_ON_ROOT( "Encountered node with " << numPhysicalTags << " physical tags!"
                                                                   << " Only using the first one" );
            for ( uint_t k = 1; k < numPhysicalTags; ++k )
            {
               file >> physicalTag;
            }
         }
      }
      else if ( importPhysicalTags_ )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "Point with tag " << pointTag << " has " << numPhysicalTags << " physical tags!" );
         numNonTagged++;
      }
   }

   // read info on curve entities (dim=1)
   for ( uint_t idx = 0; idx < numCurves; ++idx )
   {
      uint_t curveTag;
      real_t ignoreFloat;
      uint_t numPhysicalTags;
      uint_t physicalTag;

      file >> curveTag;
      // we do not need minX, ..., maxZ
      file >> ignoreFloat >> ignoreFloat >> ignoreFloat >> ignoreFloat >> ignoreFloat >> ignoreFloat;
      file >> numPhysicalTags;

      WALBERLA_LOG_PROGRESS_ON_ROOT( "Processing curve with curveTag = " << curveTag
                                                                         << ", numPhysicalTags = " << numPhysicalTags );
      if ( numPhysicalTags >= 1u )
      {
         file >> physicalTag;
         marker entityMarker{ 1u, curveTag };
         entityToPhysicalTag[entityMarker] = physicalTag;

         if ( numPhysicalTags > 1u )
         {
            WALBERLA_LOG_WARNING_ON_ROOT( "Encountered curve with " << numPhysicalTags << " physical tags!"
                                                                    << " Only using the first one" );
            for ( uint_t k = 1; k < numPhysicalTags; ++k )
            {
               file >> physicalTag;
            }
         }
      }
      else if ( importPhysicalTags_ )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "Curve with tag " << curveTag << " has " << numPhysicalTags << " physical tags!" );
         numNonTagged++;
      }

      // discard remaining curve info (pointTag can be negative)
      uint_t numBoundingPoints;
      int    pointTag;
      file >> numBoundingPoints;
      for ( uint_t k = 0; k < numBoundingPoints; ++k )
      {
         file >> pointTag;
      }
   }

   // read info on surfaces entities (dim=2)
   for ( uint_t idx = 0; idx < numSurfaces; ++idx )
   {
      uint_t surfaceTag;
      real_t ignoreFloat;
      uint_t numPhysicalTags;
      uint_t physicalTag;

      file >> surfaceTag;
      // we do not need minX, ..., maxZ
      file >> ignoreFloat >> ignoreFloat >> ignoreFloat >> ignoreFloat >> ignoreFloat >> ignoreFloat;
      file >> numPhysicalTags;

      WALBERLA_LOG_PROGRESS_ON_ROOT( "Processing surface with surfaceTag = " << surfaceTag
                                                                             << ", numPhysicalTags = " << numPhysicalTags );
      if ( numPhysicalTags >= 1u )
      {
         file >> physicalTag;
         marker entityMarker{ 2u, surfaceTag };
         entityToPhysicalTag[entityMarker] = physicalTag;

         if ( numPhysicalTags > 1u )
         {
            WALBERLA_LOG_WARNING_ON_ROOT( "Encountered surface with " << numPhysicalTags << " physical tags!"
                                                                      << " Only using the first one" );
            for ( uint_t k = 1; k < numPhysicalTags; ++k )
            {
               file >> physicalTag;
            }
         }
      }
      else if ( importPhysicalTags_ )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "Surface with tag " << surfaceTag << " has " << numPhysicalTags << " physical tags!" );
         numNonTagged++;
      }

      // discard remaining surface info (curveTag can be negative)
      uint_t numBoundingCurves;
      int    curveTag;
      file >> numBoundingCurves;
      for ( uint_t k = 0; k < numBoundingCurves; ++k )
      {
         file >> curveTag;
      }
   }

   // read info on volume entities (dim=3)
   for ( uint_t idx = 0; idx < numVolumes; ++idx )
   {
      uint_t volumeTag;
      real_t ignoreFloat;
      uint_t numPhysicalTags;
      uint_t physicalTag;

      file >> volumeTag;
      // we do not need minX, ..., maxZ
      file >> ignoreFloat >> ignoreFloat >> ignoreFloat >> ignoreFloat >> ignoreFloat >> ignoreFloat;
      file >> numPhysicalTags;

      WALBERLA_LOG_PROGRESS_ON_ROOT( "Processing volume with volumeTag = " << volumeTag
                                                                           << ", numPhysicalTags = " << numPhysicalTags );
      if ( numPhysicalTags >= 1u )
      {
         file >> physicalTag;
         marker entityMarker{ 3u, volumeTag };
         entityToPhysicalTag[entityMarker] = physicalTag;

         if ( numPhysicalTags > 1u )
         {
            WALBERLA_LOG_WARNING_ON_ROOT( "Encountered volume with " << numPhysicalTags << " physical tags!"
                                                                     << " Only using the first one" );
            for ( uint_t k = 1; k < numPhysicalTags; ++k )
            {
               file >> physicalTag;
            }
         }
      }
      else if ( importPhysicalTags_ )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "Volume with tag " << volumeTag << " has " << numPhysicalTags << " physical tags!" );
         numNonTagged++;
      }

      // discard remaining surface info (surfaceTag can be negative)
      uint_t numBoundingSurfaces;
      int    surfaceTag;
      file >> numBoundingSurfaces;
      for ( uint_t k = 0; k < numBoundingSurfaces; ++k )
      {
         file >> surfaceTag;
      }
   }

   if ( importPhysicalTags_ && numNonTagged > 0 )
   {
      WALBERLA_LOG_WARNING_ON_ROOT(
          "For importPhysicalTags = true we expect entities to have physical tags!\n"
          << numNonTagged << ( numNonTagged > 1 ? " entities" : " entity" ) << " in the file"
          << ( numNonTagged > 1 ? " do" : " does" ) << " not have a physical tag!\n"
          << "This can only work for redundant entities, such as a point not used in the mesh!\n"
          << "Proceeding ... but keep your fingers crossed!" );
   }

   return entityToPhysicalTag;
}

// import information on the nodes defined in the MSH file
void GmshReaderForMSH41::readSectionNodes( std::ifstream&                    file,
                                           MeshInfo&                         meshInfo,
                                           const std::map< marker, uint_t >& entityToPhysicalTag ) const
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

      // set boundary flag for nodes in this block
      uint_t boundaryFlag = 0u;

      if ( importPhysicalTags_ )
      {
         marker onMyMark{ entityDim, entityTag };
         WALBERLA_CHECK_EQUAL( entityToPhysicalTag.count( onMyMark ), 1, "Entity not present in entity-to-physical-tag map!" );
         boundaryFlag = entityToPhysicalTag.at( onMyMark );
      }

      WALBERLA_LOG_PROGRESS_ON_ROOT( "" << numNodesInBlock << " node(s) for entity with (dim, tag) = (" << entityDim << ", "
                                        << entityTag << ") -> boundaryFlag = " << boundaryFlag );

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
         meshInfo.vertices_[id] = MeshInfo::Vertex( id, Point3D( xCoord, yCoord, zCoord ), boundaryFlag );
      }
   }

   WALBERLA_CHECK_EQUAL( numNodes, nodeCount, "Inconsistency w.r.t. number of nodes promised and present!" );
}

// import information on the elements defined in the MSH file
void GmshReaderForMSH41::readSectionElements( std::ifstream&                    file,
                                              MeshInfo&                         meshInfo,
                                              const std::map< marker, uint_t >& entityToPhysicalTag ) const

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
   const uint_t oneNodePoint        = 15u;
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

      // set boundary flag for elements in this block
      uint_t boundaryFlag = 0u;

      if ( importPhysicalTags_ )
      {
         marker onMyMark{ entityDim, entityTag };
         WALBERLA_CHECK_EQUAL( entityToPhysicalTag.count( onMyMark ), 1, "Entity not present in entity-to-physical-tag map!" );
         boundaryFlag = entityToPhysicalTag.at( onMyMark );
      }

      WALBERLA_LOG_PROGRESS_ON_ROOT( "" << numElementsInBlock << " element(s) of type " << elementType
                                        << " for entity with (dim, tag) = (" << entityDim << ", " << entityTag
                                        << ") -> boundaryFlag = " << boundaryFlag );

      switch ( elementType )
      {
      case oneNodePoint: {
         // we can ignore this, as all nodes from the Nodes section
         // will become macro-vertices
         uint_t ignore;
         file >> ignore >> ignore;
         break;
      }
      case twoNodeLine: {
         uint_t                            tag = 0;
         std::array< MeshInfo::IDType, 2 > edgeNodes;

         for ( uint_t idx = 0; idx < numElementsInBlock; ++idx )
         {
            file >> tag >> edgeNodes[0] >> edgeNodes[1];
            parsedEdges[edgeNodes] = MeshInfo::Edge( edgeNodes, boundaryFlag );
         }
         break;
      }
      case threeNodeTriangle: {
         uint_t                          tag = 0;
         std::vector< MeshInfo::IDType > triangleNodes( 3 );
         for ( uint_t idx = 0; idx < numElementsInBlock; ++idx )
         {
            file >> tag >> triangleNodes[0] >> triangleNodes[1] >> triangleNodes[2];
            parsedFaces[triangleNodes] = MeshInfo::Face( triangleNodes, boundaryFlag );
         }
         break;
      }
      case fourNodeTetrahedron: {
         uint_t                          tag = 0;
         std::vector< MeshInfo::IDType > tetrahedronNodes( 4 );
         for ( uint_t idx = 0; idx < numElementsInBlock; ++idx )
         {
            file >> tag >> tetrahedronNodes[0] >> tetrahedronNodes[1] >> tetrahedronNodes[2] >> tetrahedronNodes[3];
            parsedCells[tetrahedronNodes] = MeshInfo::Cell( tetrahedronNodes, boundaryFlag );
         }
         break;
      }
      default: {
         WALBERLA_ABORT( "Detected unsupported element type: " << elementType );
      }
      }
   }

   WALBERLA_CHECK_EQUAL( numElements, elementCount, "Inconsistency w.r.t. number of elements promised and present!" );

   // --------------------------------------------
   //  second pass: insert primitives into object
   // --------------------------------------------

   WALBERLA_LOG_PROGRESS_ON_ROOT( "----------------------------------------------" );
   WALBERLA_LOG_PROGRESS_ON_ROOT( "Performing second pass" );

   // The MSH file does not contain all elements we need explicitely. E.g. it only contains the triangle
   // elements for surface entities. Hence, we must derive edges from the faces we found and faces from
   // cells.
   //
   // Since all elements we found belong to entities, we can use the parent's boundaryFlag_ for the
   // newly created lower-dimensional ones. As the second pass does not overwrites, this is safe.
   meshInfo.processPrimitivesFromGmshFile( parsedEdges, parsedFaces, parsedCells, importPhysicalTags_ );
}

} // namespace hyteg

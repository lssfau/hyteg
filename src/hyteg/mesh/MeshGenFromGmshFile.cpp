/*
 * Copyright (c) 2017-2024 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/mesh/GmshReaderForMSH22.hpp"
#include "hyteg/mesh/GmshReaderForMSH41.hpp"
#include "hyteg/mesh/MeshInfo.hpp"

namespace hyteg {

MeshInfo MeshInfo::fromGmshFile( const std::string& meshFileName, bool importPhysicalTags )
{
   MeshInfo meshInfo;

   std::ifstream meshFile;
   meshFile.open( meshFileName.c_str() );

   WALBERLA_CHECK( !!meshFile, "[Mesh] Error opening file: " << meshFileName );

   std::string token;
   meshFile >> token; // $MeshFormat

   WALBERLA_CHECK_EQUAL( token, "$MeshFormat", "[Mesh] Missing: $MeshFormat" );

   meshFile >> token; // version

   meshFile.close();

   // delegate work to reader for specific MSH format
   if ( token == "2.2" )
   {
      if ( !importPhysicalTags )
      {
         WALBERLA_ABORT( "MeshReader for MSH 2.2 does not support 'importPhysicalTags=true'.\n"
                         << "It is will always import physical tags (to some extent)!\n Please see documentation for details." );
      }

      meshInfo = GmshReaderForMSH22::readMesh( meshFileName );
   }
   else if ( token == "4.1" )
   {
      GmshReaderForMSH41 mshReader( meshFileName, importPhysicalTags );
      meshInfo = mshReader.readMesh();
   }
   else
   {
      WALBERLA_ABORT( "MSH format " << token << " is not supported.\n\n"
                                    << "If you definitely need it, please open an issue on:\n"
                                    << "https://i10git.cs.fau.de/hyteg/hyteg\n\n"
                                    << "Currently supported formats are 2.2 and 4.1." );
   }

   return meshInfo;
}

void MeshInfo::processPrimitivesFromGmshFile( const EdgeContainer& parsedEdges,
                                              const FaceContainer& parsedFaces,
                                              const CellContainer& parsedCells,
                                              bool                 inheritParentBoundaryFlag )
{
   for ( const auto& it : parsedEdges )
   {
      const MeshInfo::Edge meshInfoEdge = it.second;
      this->addEdge( meshInfoEdge );
   }

   // If the edges associated with the face were not already added, add them
   for ( const auto& it : parsedFaces )
   {
      const std::vector< MeshInfo::IDType > faceCoordinates = it.first;
      const MeshInfo::Face                  meshInfoFace    = it.second;

      WALBERLA_ASSERT_EQUAL( faceCoordinates.size(), 3, "[Mesh] Only triangle faces supported." );

      // Determine type of new edge; by default it will be an Inner edge
      uint_t boundaryFlag = inheritParentBoundaryFlag ? meshInfoFace.getBoundaryFlag() : 0u;

      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { faceCoordinates[0], faceCoordinates[1] } } ), boundaryFlag ) );
      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { faceCoordinates[1], faceCoordinates[2] } } ), boundaryFlag ) );
      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { faceCoordinates[2], faceCoordinates[0] } } ), boundaryFlag ) );

      this->addFace( meshInfoFace );
   }

   // If the edges and faces associated with the cell were not already added, add them
   for ( const auto& it : parsedCells )
   {
      const std::vector< MeshInfo::IDType > cellCoordinates = it.first;
      const MeshInfo::Cell                  meshInfoCell    = it.second;

      WALBERLA_ASSERT_EQUAL( cellCoordinates.size(), 4, "[Mesh] Only tetrahedron cells supported." );

      // Determine type of new edge or face; by default it will be an Inner edge/face
      uint_t boundaryFlag = inheritParentBoundaryFlag ? meshInfoCell.getBoundaryFlag() : 0u;

      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { cellCoordinates[0], cellCoordinates[1] } } ), boundaryFlag ) );
      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { cellCoordinates[0], cellCoordinates[2] } } ), boundaryFlag ) );
      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { cellCoordinates[0], cellCoordinates[3] } } ), boundaryFlag ) );
      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { cellCoordinates[1], cellCoordinates[2] } } ), boundaryFlag ) );
      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { cellCoordinates[1], cellCoordinates[3] } } ), boundaryFlag ) );
      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { cellCoordinates[2], cellCoordinates[3] } } ), boundaryFlag ) );

      this->addFace( Face( std::vector< MeshInfo::IDType >( { { cellCoordinates[0], cellCoordinates[1], cellCoordinates[2] } } ),
                           boundaryFlag ) );
      this->addFace( Face( std::vector< MeshInfo::IDType >( { { cellCoordinates[0], cellCoordinates[1], cellCoordinates[3] } } ),
                           boundaryFlag ) );
      this->addFace( Face( std::vector< MeshInfo::IDType >( { { cellCoordinates[0], cellCoordinates[2], cellCoordinates[3] } } ),
                           boundaryFlag ) );
      this->addFace( Face( std::vector< MeshInfo::IDType >( { { cellCoordinates[1], cellCoordinates[2], cellCoordinates[3] } } ),
                           boundaryFlag ) );

      this->cells_[cellCoordinates] = meshInfoCell;
   }
}

} // namespace hyteg

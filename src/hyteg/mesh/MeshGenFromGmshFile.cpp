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

MeshInfo MeshInfo::fromGmshFile( const std::string& meshFileName )
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
      meshInfo = GmshReaderForMSH22::readMesh( meshFileName );
   }
   else if ( token == "4.1" )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "Support for MSH 4.1 currently under construction" );
      GmshReaderForMSH41 mshReader( meshFileName );
      meshInfo = mshReader.readMesh();
   }
   else
   {
      WALBERLA_ABORT( "MSH format " << token << " is not supported.\n"
                                    << "If you definitely need it, please open an issue on "
                                    << "https://i10git.cs.fau.de/hyteg/hyteg" );
   }

   return meshInfo;
}

void MeshInfo::processPrimitivesFromGmshFile( const EdgeContainer& parsedEdges,
                                              const FaceContainer& parsedFaces,
                                              const CellContainer& parsedCells )
{
   for ( const auto& it : parsedEdges )
   {
      const MeshInfo::Edge meshInfoEdge = it.second;
      this->addEdge( meshInfoEdge );
   }

   for ( const auto& it : parsedFaces )
   {
      const std::vector< MeshInfo::IDType > faceCoordinates = it.first;
      const MeshInfo::Face                  meshInfoFace    = it.second;

      // If the corresponding edge was not already added, add an edge of type Inner
      WALBERLA_ASSERT_EQUAL( faceCoordinates.size(), 3, "[Mesh] Only triangle faces supported." );

      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { faceCoordinates[0], faceCoordinates[1] } } ), 0 ) );
      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { faceCoordinates[1], faceCoordinates[2] } } ), 0 ) );
      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { faceCoordinates[2], faceCoordinates[0] } } ), 0 ) );

      this->addFace( meshInfoFace );
   }

   for ( const auto& it : parsedCells )
   {
      const std::vector< MeshInfo::IDType > cellCoordinates = it.first;
      const MeshInfo::Cell                  meshInfoCell    = it.second;

      // If the corresponding edge was not already added, add an edge of type Inner
      WALBERLA_ASSERT_EQUAL( cellCoordinates.size(), 4, "[Mesh] Only tetrahedron cells supported." );

      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { cellCoordinates[0], cellCoordinates[1] } } ), 0 ) );
      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { cellCoordinates[0], cellCoordinates[2] } } ), 0 ) );
      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { cellCoordinates[0], cellCoordinates[3] } } ), 0 ) );
      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { cellCoordinates[1], cellCoordinates[2] } } ), 0 ) );
      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { cellCoordinates[1], cellCoordinates[3] } } ), 0 ) );
      this->addEdge( Edge( std::array< MeshInfo::IDType, 2 >( { { cellCoordinates[2], cellCoordinates[3] } } ), 0 ) );

      this->addFace(
          Face( std::vector< MeshInfo::IDType >( { { cellCoordinates[0], cellCoordinates[1], cellCoordinates[2] } } ), 0 ) );
      this->addFace(
          Face( std::vector< MeshInfo::IDType >( { { cellCoordinates[0], cellCoordinates[1], cellCoordinates[3] } } ), 0 ) );
      this->addFace(
          Face( std::vector< MeshInfo::IDType >( { { cellCoordinates[0], cellCoordinates[2], cellCoordinates[3] } } ), 0 ) );
      this->addFace(
          Face( std::vector< MeshInfo::IDType >( { { cellCoordinates[1], cellCoordinates[2], cellCoordinates[3] } } ), 0 ) );

      this->cells_[cellCoordinates] = meshInfoCell;
   }
}

} // namespace hyteg

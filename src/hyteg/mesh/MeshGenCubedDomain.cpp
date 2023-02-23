/*
 * Copyright (c) 2017-2012 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

using walberla::int_c;
using walberla::uint_c;
using walberla::uint_t;

static uint_t indexInCube( int x, int y, int z, int xMin, int yMin, int zMin, uint_t nx, uint_t ny )
{
   const auto id = ( z - zMin ) * int( nx * ny ) + ( y - yMin ) * int( nx ) + ( x - xMin );
   WALBERLA_ASSERT_GREATER_EQUAL( id, 0 );
   return uint_c( id );
}

static uint_t symmCubeID( int x, int y, int z, uint_t localVertexID, int xMin, int yMin, int zMin, uint_t numMaxCubesX, uint_t numMaxCubesY, uint_t numMaxCubesZ )
{
   const auto numCubeCornerVerticesX = numMaxCubesX + 1;
   const auto numCubeCornerVerticesY = numMaxCubesY + 1;
   const auto numCubeCornerVerticesZ = numMaxCubesZ + 1;

   const auto numXYFaces = numMaxCubesX * numMaxCubesY * (numMaxCubesZ + 1);
   const auto numXZFaces = numMaxCubesX * numMaxCubesZ * (numMaxCubesY + 1);

   const auto numCubeCornerVertices = numCubeCornerVerticesX * numCubeCornerVerticesY * numCubeCornerVerticesZ;

   uint_t     vertexID;

   switch ( localVertexID )
   {
   case 0:
      vertexID = indexInCube( x, y, z, xMin, yMin, zMin, numCubeCornerVerticesX, numCubeCornerVerticesY );
      break;
   case 1:
      vertexID = indexInCube( x + 1, y, z, xMin, yMin, zMin, numCubeCornerVerticesX, numCubeCornerVerticesY );
      break;
   case 3:
      vertexID = indexInCube( x, y + 1, z, xMin, yMin, zMin, numCubeCornerVerticesX, numCubeCornerVerticesY );
      break;
   case 4:
      vertexID = indexInCube( x + 1, y + 1, z, xMin, yMin, zMin, numCubeCornerVerticesX, numCubeCornerVerticesY );
      break;
   case 9:
      vertexID = indexInCube( x, y, z + 1, xMin, yMin, zMin, numCubeCornerVerticesX, numCubeCornerVerticesY );
      break;
   case 10:
      vertexID = indexInCube( x + 1, y, z + 1, xMin, yMin, zMin, numCubeCornerVerticesX, numCubeCornerVerticesY );
      break;
   case 12:
      vertexID = indexInCube( x, y + 1, z + 1, xMin, yMin, zMin, numCubeCornerVerticesX, numCubeCornerVerticesY );
      break;
   case 13:
      vertexID = indexInCube( x + 1, y + 1, z + 1, xMin, yMin, zMin, numCubeCornerVerticesX, numCubeCornerVerticesY );
      break;

   case 2:
      vertexID = numCubeCornerVertices + indexInCube( x, y, z, xMin, yMin, zMin, numMaxCubesX, numMaxCubesY );
      break;
   case 11:
      vertexID = numCubeCornerVertices + indexInCube( x, y, z + 1, xMin, yMin, zMin, numMaxCubesX, numMaxCubesY );
      break;

   case 5:
      vertexID =
          numCubeCornerVertices + numXYFaces + indexInCube( x, y, z, xMin, yMin, zMin, numMaxCubesX, numMaxCubesY + 1 );
      break;
   case 8:
      vertexID = numCubeCornerVertices + numXYFaces +
                 indexInCube( x, y + 1, z, xMin, yMin, zMin, numMaxCubesX, numMaxCubesY + 1 );
      break;

   case 6:
      vertexID = numCubeCornerVertices + numXYFaces + numXZFaces +
                 indexInCube( x, y, z, xMin, yMin, zMin, numMaxCubesX + 1, numMaxCubesY );
      break;
   case 7:
      vertexID = numCubeCornerVertices + numXYFaces + numXZFaces +
                 indexInCube( x + 1, y, z, xMin, yMin, zMin, numMaxCubesX + 1, numMaxCubesY );
      break;
   default:
      WALBERLA_ABORT( "Invalid local vertex ID." );
      break;
   }
   return vertexID;
}

MeshInfo MeshInfo::meshCubedDomain( const std::set< std::array< int, 3 > >& cubeCoordinates, int cubeType )
{
   if ( cubeType < 0 || cubeType > 1 )
   {
      WALBERLA_ABORT( "Invalid cube type." );
   }

   const real_t cubeSideLength = walberla::real_c( 1.0 );

   MeshInfo meshInfo;

   int xMax = std::numeric_limits< int >::min();
   int yMax = std::numeric_limits< int >::min();
   int zMax = std::numeric_limits< int >::min();
   int xMin = std::numeric_limits< int >::max();
   int yMin = std::numeric_limits< int >::max();
   int zMin = std::numeric_limits< int >::max();
   for ( const auto& it : cubeCoordinates )
   {
      xMax = std::max( it.at( 0 ), xMax );
      yMax = std::max( it.at( 1 ), yMax );
      zMax = std::max( it.at( 2 ), zMax );
      xMin = std::min( it.at( 0 ), xMin );
      yMin = std::min( it.at( 1 ), yMin );
      zMin = std::min( it.at( 2 ), zMin );
   }

   // increment since we need largest "outer" vertices of cubes
   xMax++;
   yMax++;
   zMax++;

   WALBERLA_ASSERT_GREATER_EQUAL( xMax - xMin, 0 );
   WALBERLA_ASSERT_GREATER_EQUAL( yMax - yMin, 0 );
   WALBERLA_ASSERT_GREATER_EQUAL( zMax - zMin, 0 );

   const uint_t numCubeCornerVerticesX = uint_c( xMax - xMin ) + 1;
   const uint_t numCubeCornerVerticesY = uint_c( yMax - yMin ) + 1;
   const uint_t numCubeCornerVerticesZ = uint_c( zMax - zMin ) + 1;

   const uint_t numMaxCubesX = numCubeCornerVerticesX - 1;
   const uint_t numMaxCubesY = numCubeCornerVerticesY - 1;
   const uint_t numMaxCubesZ = numCubeCornerVerticesZ - 1;

   ///// MAPPINGS FOR CUBE TYPE == 0 /////

   // mapping of tetrahedron to vertices of local cell (this is standard
   // hexahedron numbering)
   int tNode[6][4];
   tNode[0][0] = 4;
   tNode[0][1] = 5;
   tNode[0][2] = 7;
   tNode[0][3] = 3;
   tNode[1][0] = 0;
   tNode[1][1] = 3;
   tNode[1][2] = 1;
   tNode[1][3] = 4;
   tNode[2][0] = 4;
   tNode[2][1] = 1;
   tNode[2][2] = 5;
   tNode[2][3] = 3;
   tNode[3][0] = 5;
   tNode[3][1] = 6;
   tNode[3][2] = 7;
   tNode[3][3] = 2;
   tNode[4][0] = 1;
   tNode[4][1] = 3;
   tNode[4][2] = 2;
   tNode[4][3] = 5;
   tNode[5][0] = 7;
   tNode[5][1] = 2;
   tNode[5][2] = 3;
   tNode[5][3] = 5;

   // determine offsets for computing address tuples for the eight
   // vertices of a local cuboid (w.r.t. lower left front vertex)
   int offset[8][3];
   offset[0][0] = 0;
   offset[0][1] = 0;
   offset[0][2] = 0;
   offset[1][0] = 1;
   offset[1][1] = 0;
   offset[1][2] = 0;
   offset[2][0] = 1;
   offset[2][1] = 1;
   offset[2][2] = 0;
   offset[3][0] = 0;
   offset[3][1] = 1;
   offset[3][2] = 0;
   offset[4][0] = 0;
   offset[4][1] = 0;
   offset[4][2] = 1;
   offset[5][0] = 1;
   offset[5][1] = 0;
   offset[5][2] = 1;
   offset[6][0] = 1;
   offset[6][1] = 1;
   offset[6][2] = 1;
   offset[7][0] = 0;
   offset[7][1] = 1;
   offset[7][2] = 1;

   ///// MAPPINGS FOR CUBE TYPE == 1 /////

   /**
    * Vertex numbering for cube type == 1 (symmetric cube with 24 tets)
    *
    * z == 0 | z == 1 | z == 2
    *
    * 3   4      8      12    13
    *   2      6   7       11
    * 0   1      5      9     10
    *
    */

   const Point3D subCubeWidth( {cubeSideLength, cubeSideLength, cubeSideLength} );

   std::vector< Point3D > subCubeVertexCoordinates( 14 );

   subCubeVertexCoordinates[0] = Point3D( 0, 0, 0 );
   subCubeVertexCoordinates[1] = Point3D( subCubeWidth[0], 0, 0 );
   subCubeVertexCoordinates[2] = Point3D( subCubeWidth[0] / 2, subCubeWidth[1] / 2, 0 );
   subCubeVertexCoordinates[3] = Point3D( 0, subCubeWidth[1], 0 );
   subCubeVertexCoordinates[4] = Point3D( subCubeWidth[0], subCubeWidth[1], 0 );

   subCubeVertexCoordinates[5] = Point3D( subCubeWidth[0] / 2, 0, subCubeWidth[2] / 2 );
   subCubeVertexCoordinates[6] = Point3D( 0, subCubeWidth[1] / 2, subCubeWidth[2] / 2 );
   subCubeVertexCoordinates[7] = Point3D( subCubeWidth[0], subCubeWidth[1] / 2, subCubeWidth[2] / 2 );
   subCubeVertexCoordinates[8] = Point3D( subCubeWidth[0] / 2, subCubeWidth[1], subCubeWidth[2] / 2 );

   subCubeVertexCoordinates[9]  = Point3D( 0, 0, subCubeWidth[2] );
   subCubeVertexCoordinates[10] = Point3D( subCubeWidth[0], 0, subCubeWidth[2] );
   subCubeVertexCoordinates[11] = Point3D( subCubeWidth[0] / 2, subCubeWidth[1] / 2, subCubeWidth[2] );
   subCubeVertexCoordinates[12] = Point3D( 0, subCubeWidth[1], subCubeWidth[2] );
   subCubeVertexCoordinates[13] = Point3D( subCubeWidth[0], subCubeWidth[1], subCubeWidth[2] );

   std::vector< std::vector< uint_t > > cellVertices;
   cellVertices.push_back( {{0, 1, 2, 5}} );
   cellVertices.push_back( {{0, 3, 2, 6}} );
   cellVertices.push_back( {{1, 2, 4, 7}} );
   cellVertices.push_back( {{2, 4, 3, 8}} );

   cellVertices.push_back( {{2, 0, 5, 6}} );
   cellVertices.push_back( {{2, 1, 5, 7}} );
   cellVertices.push_back( {{2, 3, 6, 8}} );
   cellVertices.push_back( {{2, 4, 7, 8}} );

   cellVertices.push_back( {{0, 9, 5, 6}} );
   cellVertices.push_back( {{1, 10, 5, 7}} );
   cellVertices.push_back( {{3, 12, 6, 8}} );
   cellVertices.push_back( {{4, 13, 7, 8}} );

   cellVertices.push_back( {{2, 5, 6, 8}} );
   cellVertices.push_back( {{2, 5, 7, 8}} );
   cellVertices.push_back( {{11, 5, 6, 8}} );
   cellVertices.push_back( {{11, 5, 7, 8}} );

   cellVertices.push_back( {{11, 9, 5, 6}} );
   cellVertices.push_back( {{11, 10, 5, 7}} );
   cellVertices.push_back( {{11, 12, 6, 8}} );
   cellVertices.push_back( {{11, 13, 7, 8}} );

   cellVertices.push_back( {{11, 9, 10, 5}} );
   cellVertices.push_back( {{11, 9, 12, 6}} );
   cellVertices.push_back( {{11, 10, 13, 7}} );
   cellVertices.push_back( {{11, 12, 13, 8}} );

   ///////////////////////////////////////

   // compute vertices and insert them
   for ( const auto& it : cubeCoordinates )
   {
      const auto x = it.at( 0 );
      const auto y = it.at( 1 );
      const auto z = it.at( 2 );

      if ( cubeType == 0 )
      {
         for ( int xInc = 0; xInc <= 1; xInc++ )
            for ( int yInc = 0; yInc <= 1; yInc++ )
               for ( int zInc = 0; zInc <= 1; zInc++ )
               {
                  const auto vertexID = indexInCube(
                      x + xInc, y + yInc, z + zInc, xMin, yMin, zMin, numCubeCornerVerticesX, numCubeCornerVerticesY );

                  if ( meshInfo.vertices_.count( vertexID ) == 0 )
                  {
                     const auto    xPos = static_cast< real_t >( x + xInc ) * cubeSideLength;
                     const auto    yPos = static_cast< real_t >( y + yInc ) * cubeSideLength;
                     const auto    zPos = static_cast< real_t >( z + zInc ) * cubeSideLength;
                     const Point3D vertexCoords( {xPos, yPos, zPos} );
                     meshInfo.vertices_[vertexID] = MeshInfo::Vertex( vertexID, vertexCoords, 0 );
                  }
               }
      }
      else if ( cubeType == 1 )
      {
         for ( uint_t localVertexID = 0; localVertexID < 14; localVertexID++ )
         {
            const auto vertexID = symmCubeID( x, y, z, localVertexID, xMin, yMin, zMin, numMaxCubesX, numMaxCubesY, numMaxCubesZ );

            if ( meshInfo.vertices_.count( vertexID ) == 0 )
            {
               const auto    xPos           = static_cast< real_t >( x ) * cubeSideLength;
               const auto    yPos           = static_cast< real_t >( y ) * cubeSideLength;
               const auto    zPos           = static_cast< real_t >( z ) * cubeSideLength;
               const Point3D vertexCoords   = Point3D( xPos, yPos, zPos ) + subCubeVertexCoordinates[localVertexID];
               meshInfo.vertices_[vertexID] = MeshInfo::Vertex( vertexID, vertexCoords, 0 );
            }
         }
      }
   }

   std::vector< uint_t > vertexIDs;
   int                   xid = 0;
   int                   yid = 0;
   int                   zid = 0;

   // loop over cells
   for ( const auto& it : cubeCoordinates )
   {
      const auto ix = it.at( 0 );
      const auto iy = it.at( 1 );
      const auto iz = it.at( 2 );

      if ( cubeType == 0 )
      {
         // split cell into six tets
         for ( uint_t iTet = 0; iTet < 6; ++iTet )
         {
            // determine list of global indices of the four
            // verts making up this tet
            vertexIDs.clear();
            for ( uint_t m = 0; m < 4; ++m )
            {
               xid = ix + offset[tNode[iTet][m]][0];
               yid = iy + offset[tNode[iTet][m]][1];
               zid = iz + offset[tNode[iTet][m]][2];
               vertexIDs.push_back(
                   indexInCube( xid, yid, zid, xMin, yMin, zMin, numCubeCornerVerticesX, numCubeCornerVerticesY ) );
            }

            // insert tet info
            const MeshInfo::Cell cell( vertexIDs, 0 );
            meshInfo.addCellAndAllEdgesAndFaces( cell );
            // std::cout << "Inserted cell with vertices ["
            //           << vertexIDs[0] << ", "
            //           << vertexIDs[1] << ", "
            //           << vertexIDs[2] << ", "
            //           << vertexIDs[3] << "]"
            //           << "\n";
         }
      }
      else if ( cubeType == 1 )
      {
         for ( uint_t i = 0; i < 24; i++ )
         {
            std::vector< uint_t > cv( 4 );
            for ( uint_t ii = 0; ii < 4; ii++ )
               cv[ii] =
                   symmCubeID( ix, iy, iz, cellVertices[i][ii], xMin, yMin, zMin, numMaxCubesX, numMaxCubesY, numMaxCubesZ );

            const MeshInfo::Cell cell( cv, 0 );
            meshInfo.addCellAndAllEdgesAndFaces( cell );
         }
      }
   }

   return meshInfo;
}

} // namespace hyteg

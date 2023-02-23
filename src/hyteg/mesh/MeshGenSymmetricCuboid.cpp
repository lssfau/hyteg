/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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
#include "hyteg/types/PointND.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;

MeshInfo MeshInfo::meshSymmetricCuboid( const Point3D lowerLeftFront,
                                        const Point3D upperRightBack,
                                        uint_t        numCubesX,
                                        uint_t        numCubesY,
                                        uint_t        numCubesZ )
{
   /**
   * Vertex numbering
   *
   * z == 0 | z == 1 | z == 2
   *
   * 3   4      8      12    13
   *   2      6   7       11
   * 0   1      5      9     10
   *
   */

   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();

   if ( numCubesX == 0 || numCubesY == 0 || numCubesZ == 0 )
      return meshInfo;

   Point3D subCubeWidth = upperRightBack - lowerLeftFront;
   subCubeWidth[0] /= real_c( numCubesX );
   subCubeWidth[1] /= real_c( numCubesY );
   subCubeWidth[2] /= real_c( numCubesZ );

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

   std::vector< std::vector< int > > cellVertices;
   cellVertices.push_back( { { 0, 1, 2, 5 } } );
   cellVertices.push_back( { { 0, 3, 2, 6 } } );
   cellVertices.push_back( { { 1, 2, 4, 7 } } );
   cellVertices.push_back( { { 2, 4, 3, 8 } } );

   cellVertices.push_back( { { 2, 0, 5, 6 } } );
   cellVertices.push_back( { { 2, 1, 5, 7 } } );
   cellVertices.push_back( { { 2, 3, 6, 8 } } );
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

   // We collect all 14 vertices for each sub-cube
   // and build a map to the actual vertex IDs.

   int vertexID = 0;

   // map[cubeX][cubeY][cubeZ][localID] = vertexID
   std::map< int, std::map< int, std::map< int, std::map< int, int > > > > vertexMapping;

   for ( int z = 0; z < (int)numCubesZ; z++ )
   {
      for ( int y = 0; y < (int)numCubesY; y++ )
      {
         for ( int x = 0; x < (int)numCubesX; x++ )
         {
            if ( x == 0 && y == 0 && z == 0 )
            {
               for ( int i = 0; i < 14; i++ )
                  vertexMapping[x][y][z][i] = vertexID++;
            }
            else if ( x == 0 && y == 0 )
            {
               vertexMapping[x][y][z][0] = vertexMapping[x][y][z - 1][9];
               vertexMapping[x][y][z][1] = vertexMapping[x][y][z - 1][10];
               vertexMapping[x][y][z][2] = vertexMapping[x][y][z - 1][11];
               vertexMapping[x][y][z][3] = vertexMapping[x][y][z - 1][12];
               vertexMapping[x][y][z][4] = vertexMapping[x][y][z - 1][13];
               for ( int i = 5; i < 14; i++ )
                  vertexMapping[x][y][z][i] = vertexID++;
            }
            else if ( x == 0 && z == 0 )
            {
               vertexMapping[x][y][z][0]  = vertexMapping[x][y - 1][z][3];
               vertexMapping[x][y][z][1]  = vertexMapping[x][y - 1][z][4];
               vertexMapping[x][y][z][5]  = vertexMapping[x][y - 1][z][8];
               vertexMapping[x][y][z][9]  = vertexMapping[x][y - 1][z][12];
               vertexMapping[x][y][z][10] = vertexMapping[x][y - 1][z][13];

               for ( auto i : std::vector< int >( {{2, 3, 4, 6, 7, 8, 11, 12, 13}} ) )
                  vertexMapping[x][y][z][i] = vertexID++;
            }
            else if ( y == 0 && z == 0 )
            {
               vertexMapping[x][y][z][0]  = vertexMapping[x - 1][y][z][1];
               vertexMapping[x][y][z][3]  = vertexMapping[x - 1][y][z][4];
               vertexMapping[x][y][z][6]  = vertexMapping[x - 1][y][z][7];
               vertexMapping[x][y][z][9]  = vertexMapping[x - 1][y][z][10];
               vertexMapping[x][y][z][12] = vertexMapping[x - 1][y][z][13];

               for ( auto i : std::vector< int >( {{1, 2, 4, 5, 7, 8, 10, 11, 13}} ) )
                  vertexMapping[x][y][z][i] = vertexID++;
            }
            else if ( x == 0 )
            {
               vertexMapping[x][y][z][0]  = vertexMapping[x][y][z - 1][9];
               vertexMapping[x][y][z][1]  = vertexMapping[x][y][z - 1][10];
               vertexMapping[x][y][z][2]  = vertexMapping[x][y][z - 1][11];
               vertexMapping[x][y][z][3]  = vertexMapping[x][y][z - 1][12];
               vertexMapping[x][y][z][4]  = vertexMapping[x][y][z - 1][13];
               vertexMapping[x][y][z][5]  = vertexMapping[x][y - 1][z][8];
               vertexMapping[x][y][z][9]  = vertexMapping[x][y - 1][z][12];
               vertexMapping[x][y][z][10] = vertexMapping[x][y - 1][z][13];

               for ( auto i : std::vector< int >( {{6, 7, 8, 11, 12, 13}} ) )
                  vertexMapping[x][y][z][i] = vertexID++;
            }
            else if ( y == 0 )
            {
               vertexMapping[x][y][z][0]  = vertexMapping[x][y][z - 1][9];
               vertexMapping[x][y][z][1]  = vertexMapping[x][y][z - 1][10];
               vertexMapping[x][y][z][2]  = vertexMapping[x][y][z - 1][11];
               vertexMapping[x][y][z][3]  = vertexMapping[x][y][z - 1][12];
               vertexMapping[x][y][z][4]  = vertexMapping[x][y][z - 1][13];
               vertexMapping[x][y][z][6]  = vertexMapping[x - 1][y][z][7];
               vertexMapping[x][y][z][9]  = vertexMapping[x - 1][y][z][10];
               vertexMapping[x][y][z][12] = vertexMapping[x - 1][y][z][13];

               for ( auto i : std::vector< int >( {{5, 7, 8, 10, 11, 13}} ) )
                  vertexMapping[x][y][z][i] = vertexID++;
            }
            else if ( z == 0 )
            {
              vertexMapping[x][y][z][0]  = vertexMapping[x - 1][y][z][1];
              vertexMapping[x][y][z][1]  = vertexMapping[x][y - 1][z][4];
              vertexMapping[x][y][z][3]  = vertexMapping[x - 1][y][z][4];
              vertexMapping[x][y][z][5]  = vertexMapping[x][y - 1][z][8];
              vertexMapping[x][y][z][6]  = vertexMapping[x - 1][y][z][7];
              vertexMapping[x][y][z][9]  = vertexMapping[x - 1][y][z][10];
              vertexMapping[x][y][z][10] = vertexMapping[x][y - 1][z][13];
              vertexMapping[x][y][z][12] = vertexMapping[x - 1][y][z][13];

               for ( auto i : std::vector< int >( {{2, 4, 7, 8, 11, 13}} ) )
                  vertexMapping[x][y][z][i] = vertexID++;
            }
            else
            {
              vertexMapping[x][y][z][0]  = vertexMapping[x][y][z - 1][9];
              vertexMapping[x][y][z][1]  = vertexMapping[x][y][z - 1][10];
              vertexMapping[x][y][z][2]  = vertexMapping[x][y][z - 1][11];
              vertexMapping[x][y][z][3]  = vertexMapping[x][y][z - 1][12];
              vertexMapping[x][y][z][4]  = vertexMapping[x][y][z - 1][13];
              vertexMapping[x][y][z][5]  = vertexMapping[x][y - 1][z][8];
              vertexMapping[x][y][z][6]  = vertexMapping[x - 1][y][z][7];
              vertexMapping[x][y][z][9]  = vertexMapping[x][y - 1][z][12];
              vertexMapping[x][y][z][10] = vertexMapping[x][y - 1][z][13];
              vertexMapping[x][y][z][12] = vertexMapping[x - 1][y][z][13];

              for ( auto i : std::vector< int >( {{7, 8, 11, 13}} ) )
                  vertexMapping[x][y][z][i] = vertexID++;
            }
         }
      }
   }

   // add the actual primitives

   for ( int z = 0; z < (int)numCubesZ; z++ )
   {
      for ( int y = 0; y < (int)numCubesY; y++ )
      {
         for ( int x = 0; x < (int)numCubesX; x++ )
         {
            // vertices

            auto offset = subCubeWidth;
            offset[0] *= real_c( x );
            offset[1] *= real_c( y );
            offset[2] *= real_c( z );

            for ( int i = 0; i < 14; i++ )
            {
               const auto coords      = offset + subCubeVertexCoordinates[uint_c(i)];
               const auto id          = uint_c(vertexMapping[x][y][z][i]);
               meshInfo.vertices_[id] = Vertex( id, coords, 0 );
            }

            // cells

            for ( uint_t i = 0; i < 24; i++ )
            {
               std::vector< uint_t > cv( 4 );
               for ( uint_t ii = 0; ii < 4; ii++ )
                  cv[ii] = uint_c( vertexMapping[x][y][z][cellVertices[i][ii]] );

               const Cell cell( cv, 0 );
               meshInfo.addCellAndAllEdgesAndFaces( cell );
            }
         }
      }
   }

   return meshInfo;
}

} // namespace hyteg
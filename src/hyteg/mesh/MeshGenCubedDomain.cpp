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

#include "hyteg/mesh/MeshInfo.hpp"

#include "core/logging/Logging.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"

#include <array>
#include <vector>

namespace hyteg {

using walberla::uint_t;
using walberla::uint_c;
using walberla::int_c;

static uint_t calculateVertexID( int x, int y, int z, int xMin, int yMin, int zMin, uint_t nx, uint_t ny )
{
   const auto id = (z - zMin) * int( nx * ny ) + (y - yMin) * int( nx ) + (x - xMin);
   WALBERLA_ASSERT_GREATER_EQUAL( id, 0 );
   return uint_c( id );
}

MeshInfo MeshInfo::meshCubedDomain( const std::vector< std::array< int, 3 > >& cubeCoordinates )
{
   const real_t cubeSideLength = 1.0;

   MeshInfo meshInfo;

   int xMax = std::numeric_limits< int >::min();
   int yMax = std::numeric_limits< int >::min();
   int zMax = std::numeric_limits< int >::min();
   int xMin = std::numeric_limits< int >::max();
   int yMin = std::numeric_limits< int >::max();
   int zMin = std::numeric_limits< int >::max();
   for ( const auto& it : cubeCoordinates )
   {
      WALBERLA_LOG_INFO( "contains cube: " << it.at(0) << ", " << it.at(1) << ", " << it.at(2) );

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

   const uint_t nx = uint_c( xMax - xMin ) + 1;
   const uint_t ny = uint_c( yMax - yMin ) + 1;

   WALBERLA_LOG_INFO( "nx " << nx );
   WALBERLA_LOG_INFO( "ny " << ny );
   WALBERLA_LOG_INFO( "xMin " << xMin );
   WALBERLA_LOG_INFO( "yMin " << yMin );
   WALBERLA_LOG_INFO( "zMin " << zMin );

   // compute vertices and insert them
   for ( const auto& it : cubeCoordinates )
   {
      const auto x = it.at( 0 );
      const auto y = it.at( 1 );
      const auto z = it.at( 2 );
      for ( int xInc = 0; xInc <= 1; xInc++ )
         for ( int yInc = 0; yInc <= 1; yInc++ )
            for ( int zInc = 0; zInc <= 1; zInc++ )
            {
               const auto vertexID = calculateVertexID( x + xInc, y + yInc, z + zInc, xMin, yMin, zMin, nx, ny );
               WALBERLA_LOG_INFO( "trying to insert vertex: " << x + xInc << ", " << y + yInc << ", " << z + zInc << "  ... id  " << vertexID );


               if ( meshInfo.vertices_.count( vertexID ) == 0 )
               {
                  const auto xPos              = static_cast< real_t >( x + xInc ) * cubeSideLength;
                  const auto yPos              = static_cast< real_t >( y + yInc ) * cubeSideLength;
                  const auto zPos              = static_cast< real_t >( z + zInc ) * cubeSideLength;
                  const Point3D vertexCoords( {xPos, yPos, zPos} );
                  meshInfo.vertices_[vertexID] = MeshInfo::Vertex( vertexID, vertexCoords, 0 );

                  WALBERLA_LOG_INFO( "v( " << vertexCoords << " ) inserted (id " << vertexID << ")" );
               }
            }
   }

   WALBERLA_LOG_INFO( "all vertices inserted" );


   // mapping of tetrahedron to vertices of local cell (this is standard
   // hexahedron numbering)
   int tNode[6][4];
   tNode[0][0] = 4;  tNode[0][1] = 5;  tNode[0][2] = 7;  tNode[0][3] = 3;
   tNode[1][0] = 0;  tNode[1][1] = 3;  tNode[1][2] = 1;  tNode[1][3] = 4;
   tNode[2][0] = 4;  tNode[2][1] = 1;  tNode[2][2] = 5;  tNode[2][3] = 3;
   tNode[3][0] = 5;  tNode[3][1] = 6;  tNode[3][2] = 7;  tNode[3][3] = 2;
   tNode[4][0] = 1;  tNode[4][1] = 3;  tNode[4][2] = 2;  tNode[4][3] = 5;
   tNode[5][0] = 7;  tNode[5][1] = 2;  tNode[5][2] = 3;  tNode[5][3] = 5;

   // determine offsets for computing address tuples for the eight
   // vertices of a local cuboid (w.r.t. lower left front vertex)
   int offset[8][3];
   offset[0][0] = 0; offset[0][1] = 0; offset[0][2] = 0;
   offset[1][0] = 1; offset[1][1] = 0; offset[1][2] = 0;
   offset[2][0] = 1; offset[2][1] = 1; offset[2][2] = 0;
   offset[3][0] = 0; offset[3][1] = 1; offset[3][2] = 0;
   offset[4][0] = 0; offset[4][1] = 0; offset[4][2] = 1;
   offset[5][0] = 1; offset[5][1] = 0; offset[5][2] = 1;
   offset[6][0] = 1; offset[6][1] = 1; offset[6][2] = 1;
   offset[7][0] = 0; offset[7][1] = 1; offset[7][2] = 1;

   std::vector< uint_t > vertexIDs;
   int                xid = 0;
   int                yid = 0;
   int                zid = 0;

   // loop over cells
   for ( const auto& it : cubeCoordinates )
   {
      const auto ix = it.at( 0 );
      const auto iy = it.at( 1 );
      const auto iz = it.at( 2 );

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
            vertexIDs.push_back( calculateVertexID( xid, yid, zid, xMin, yMin, zMin, nx, ny ) );
         }

         // insert tet info
         meshInfo.cells_[vertexIDs] = MeshInfo::Cell( vertexIDs, 0 );
         // std::cout << "Inserted cell with vertices ["
         //           << vertexIDs[0] << ", "
         //           << vertexIDs[1] << ", "
         //           << vertexIDs[2] << ", "
         //           << vertexIDs[3] << "]"
         //           << "\n";
      }
   }

   // Deduce edges and faces from cell info
   for ( const auto & it : meshInfo.getCells() ) {
      const std::vector< IDType > cellCoordinates = it.first;
      const MeshInfo::Cell        meshInfoCell    = it.second;

      WALBERLA_ASSERT_EQUAL( cellCoordinates.size(), 4, "[Mesh] Only tetrahedron cells supported." );

      meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{ cellCoordinates[0], cellCoordinates[1] }} ), 0 ) );
      meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{ cellCoordinates[0], cellCoordinates[2] }} ), 0 ) );
      meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{ cellCoordinates[0], cellCoordinates[3] }} ), 0 ) );
      meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{ cellCoordinates[1], cellCoordinates[2] }} ), 0 ) );
      meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{ cellCoordinates[1], cellCoordinates[3] }} ), 0 ) );
      meshInfo.addEdge( Edge( std::array< IDType, 2 >( {{ cellCoordinates[2], cellCoordinates[3] }} ), 0 ) );

      meshInfo.addFace( Face( std::vector< IDType >( {{ cellCoordinates[0], cellCoordinates[1], cellCoordinates[2] }} ), 0 ) );
      meshInfo.addFace( Face( std::vector< IDType >( {{ cellCoordinates[0], cellCoordinates[1], cellCoordinates[3] }} ), 0 ) );
      meshInfo.addFace( Face( std::vector< IDType >( {{ cellCoordinates[0], cellCoordinates[2], cellCoordinates[3] }} ), 0 ) );
      meshInfo.addFace( Face( std::vector< IDType >( {{ cellCoordinates[1], cellCoordinates[2], cellCoordinates[3] }} ), 0 ) );
   }

   return meshInfo;
}

}

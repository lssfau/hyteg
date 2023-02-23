/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr.
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

using walberla::real_c;

namespace hyteg {

MeshInfo MeshInfo::meshCuboid( const Point3D lowerLeftFront,
                               const Point3D upperRightBack,
                               uint_t nx,
                               uint_t ny,
                               uint_t nz )
{
  MeshInfo meshInfo;

  // extract corner coordinates
  real_t x0 = lowerLeftFront[ 0 ];
  real_t y0 = lowerLeftFront[ 1 ];
  real_t z0 = lowerLeftFront[ 2 ];

  real_t x1 = upperRightBack[ 0 ];
  real_t y1 = upperRightBack[ 1 ];
  real_t z1 = upperRightBack[ 2 ];

  // check inputs
  WALBERLA_ASSERT_LESS( x0, x1 );
  WALBERLA_ASSERT_LESS( y0, y1 );
  WALBERLA_ASSERT_LESS( z0, z1 );
  WALBERLA_ASSERT_GREATER( nx, 0 );
  WALBERLA_ASSERT_GREATER( ny, 0 );
  WALBERLA_ASSERT_GREATER( nz, 0 );

  // compute mesh spacings
  real_t hx = ( x1 - x0 ) / static_cast<real_t>( nx );
  real_t hy = ( y1 - y0 ) / static_cast<real_t>( ny );
  real_t hz = ( z1 - z0 ) / static_cast<real_t>( nz );

  // compute vertices and insert them
  Point3D vertexCoords;
  real_t xpos = real_c( 0.0 );
  real_t ypos = real_c( 0.0 );
  real_t zpos = real_c( 0.0 );
  uint_t vertexID = 0;

  for ( uint_t iz = 0; iz <= nz; ++iz ) {
    zpos = z0 + static_cast<real_t>( iz ) * hz;
    for ( uint_t iy = 0; iy <= ny; ++iy ) {
      ypos = y0 + static_cast<real_t>( iy ) * hy;
      for ( uint_t ix = 0; ix <= nx; ++ix ) {
         xpos                         = x0 + static_cast< real_t >( ix ) * hx;
         vertexCoords                 = Point3D( xpos, ypos, zpos );
         meshInfo.vertices_[vertexID] = MeshInfo::Vertex( vertexID, vertexCoords, 0 );
         ++vertexID;
      }
    }
  }

  WALBERLA_ASSERT_EQUAL( vertexID, ( nx + 1 ) * ( ny + 1 ) * ( nz + 1 ) );

  // mapping of tetrahedron to vertices of local cell (this is standard
  // hexahedron numbering)
  uint_t tNode[6][4];
  tNode[0][0] = 4;  tNode[0][1] = 5;  tNode[0][2] = 7;  tNode[0][3] = 3;
  tNode[1][0] = 0;  tNode[1][1] = 3;  tNode[1][2] = 1;  tNode[1][3] = 4;
  tNode[2][0] = 4;  tNode[2][1] = 1;  tNode[2][2] = 5;  tNode[2][3] = 3;
  tNode[3][0] = 5;  tNode[3][1] = 6;  tNode[3][2] = 7;  tNode[3][3] = 2;
  tNode[4][0] = 1;  tNode[4][1] = 3;  tNode[4][2] = 2;  tNode[4][3] = 5;
  tNode[5][0] = 7;  tNode[5][1] = 2;  tNode[5][2] = 3;  tNode[5][3] = 5;

  // determine offsets for computing address tuples for the eight
  // vertices of a local cuboid (w.r.t. lower left front vertex)
  uint_t offset[8][3];
  offset[0][0] = 0; offset[0][1] = 0; offset[0][2] = 0;
  offset[1][0] = 1; offset[1][1] = 0; offset[1][2] = 0;
  offset[2][0] = 1; offset[2][1] = 1; offset[2][2] = 0;
  offset[3][0] = 0; offset[3][1] = 1; offset[3][2] = 0;
  offset[4][0] = 0; offset[4][1] = 0; offset[4][2] = 1;
  offset[5][0] = 1; offset[5][1] = 0; offset[5][2] = 1;
  offset[6][0] = 1; offset[6][1] = 1; offset[6][2] = 1;
  offset[7][0] = 0; offset[7][1] = 1; offset[7][2] = 1;

  // map local vertex indices to linear id
  auto cuboidMap = [ nx, ny, nz ] ( uint_t i, uint_t j, uint_t k ) -> IDType {
    IDType id = (k)*(nx+1)*(ny+1)+(j)*(nx+1)+(i);
    WALBERLA_ASSERT_LESS( id, (nx+1)*(ny+1)*(nz+1) );
    WALBERLA_UNUSED( nz );
    return id;
  };

  std::vector< uint_t > vertexIDs;
  uint_t xid = 0;
  uint_t yid = 0;
  uint_t zid = 0;

  // loop over cells
  for ( uint_t iz = 0; iz < nz; ++iz ) {
    for ( uint_t iy = 0; iy < ny; ++iy ) {
      for ( uint_t ix = 0; ix < nx; ++ix ) {

        // split cell into six tets
        for ( uint_t iTet = 0; iTet < 6; ++iTet ) {

          // determine list of global indices of the four
          // verts making up this tet
          vertexIDs.clear();
          for ( uint_t m = 0; m < 4; ++m ) {
            xid = ix + offset[ tNode[iTet][m] ][0];
            yid = iy + offset[ tNode[iTet][m] ][1];
            zid = iz + offset[ tNode[iTet][m] ][2];
            vertexIDs.push_back( cuboidMap( xid, yid, zid ) );
          }

          // insert tet info
          meshInfo.cells_[ vertexIDs ] = MeshInfo::Cell( vertexIDs, 0 );
          // std::cout << "Inserted cell with vertices ["
          //           << vertexIDs[0] << ", "
          //           << vertexIDs[1] << ", "
          //           << vertexIDs[2] << ", "
          //           << vertexIDs[3] << "]"
          //           << "\n";
        }
      }
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

/*
 * Copyright (c) 2017-2021 Nils Kohl.
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
#include "core/math/Constants.h"

#include "hyteg/mesh/MeshInfo.hpp"

namespace hyteg {

using walberla::int_c;
using walberla::real_c;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

static uint_t sliceIndex( uint_t width, uint_t x, uint_t z )
{
   const uint_t rowOffset = z * ( width + 1 ) - ( ( ( z + 1 ) * ( z ) ) / 2 );
   return rowOffset + x;
}

/// Given all 6 vertex IDs of a triangluar prism hull, the following functions return the correct mapping to the three cells.

static std::vector< MeshInfo::Cell > cellsBottomLeftPrism( std::vector< uint_t > vertexIDs_XZY )
{
   std::vector< MeshInfo::Cell > cells;
   cells.push_back( MeshInfo::Cell( { vertexIDs_XZY[0], vertexIDs_XZY[1], vertexIDs_XZY[2], vertexIDs_XZY[3] }, 0 ) );
   cells.push_back( MeshInfo::Cell( { vertexIDs_XZY[3], vertexIDs_XZY[4], vertexIDs_XZY[5], vertexIDs_XZY[1] }, 0 ) );
   cells.push_back( MeshInfo::Cell( { vertexIDs_XZY[1], vertexIDs_XZY[3], vertexIDs_XZY[2], vertexIDs_XZY[5] }, 0 ) );
   return cells;
}

static std::vector< MeshInfo::Cell > cellsTopRightPrism( std::vector< uint_t > vertexIDs_XZY )
{
   std::vector< MeshInfo::Cell > cells;
   cells.push_back( MeshInfo::Cell( { vertexIDs_XZY[0], vertexIDs_XZY[4], vertexIDs_XZY[3], vertexIDs_XZY[2] }, 0 ) );
   cells.push_back( MeshInfo::Cell( { vertexIDs_XZY[0], vertexIDs_XZY[1], vertexIDs_XZY[2], vertexIDs_XZY[4] }, 0 ) );
   cells.push_back( MeshInfo::Cell( { vertexIDs_XZY[2], vertexIDs_XZY[3], vertexIDs_XZY[4], vertexIDs_XZY[5] }, 0 ) );
   return cells;
}

MeshInfo MeshInfo::meshTokamak( uint_t numSlices,
                                uint_t numRadialEdges,
                                real_t innerRadius,
                                real_t outerRadius,
                                real_t radiusZ,
                                uint_t cutSide,
                                uint_t cutTopAndBottom )
{
   MeshInfo meshInfo;

   uint_t id = 0;

   const uint_t width = numRadialEdges + 1;
   const uint_t edges = numRadialEdges;

   // These vectors store the vertex IDs of the slices for the upper and lower triangle respectively.
   std::vector< std::vector< uint_t > > slicesUpper( numSlices );
   std::vector< std::vector< uint_t > > slicesLower( numSlices );

   // We collect vertices in this list and remove them later if necessary.
   std::vector< uint_t > trashVertices;

   // Setting up the vertices first and storing the corresponding IDs
   for ( uint_t slice = 0; slice < numSlices; slice++ )
   {
      for ( uint_t z_idx = 0; z_idx < width; z_idx++ )
      {
         for ( uint_t x_idx = 0; x_idx + z_idx < width; x_idx++ )
         {
            auto radialEdgeLength     = ( outerRadius - innerRadius ) / real_c( edges );
            auto tangentialEdgeLength = radiusZ / real_c( edges );

            // radial angle
            auto phi = slice * ( 2 * pi ) / real_c( numSlices );

            // distance from origin in radial direction
            auto r = innerRadius + real_c( x_idx ) * radialEdgeLength;

            // distance from origin in z-direction
            auto height = tangentialEdgeLength * real_c( z_idx );

            Point3D vertexCoordsUpper( { r * std::cos( phi ), r * std::sin( phi ), height } );
            auto    vertexUpper = MeshInfo::Vertex( id++, vertexCoordsUpper, 0 );
            slicesUpper[slice].push_back( vertexUpper.getID() );
            meshInfo.addVertex( vertexUpper );

            if ( x_idx >= width - cutSide )
            {
               trashVertices.push_back( vertexUpper.getID() );
            }

            // We do not want to construct doubled vertices - so skipping the "middle row" (z == 0) for the bottom part of the slice.
            if ( z_idx > 0 )
            {
               Point3D vertexCoordsLower( { r * std::cos( phi ), r * std::sin( phi ), -height } );
               auto    vertexLower = MeshInfo::Vertex( id++, vertexCoordsLower, 0 );
               slicesLower[slice].push_back( vertexLower.getID() );
               meshInfo.addVertex( vertexLower );

               if ( z_idx >= width - cutTopAndBottom )
               {
                  trashVertices.push_back( vertexUpper.getID() );
                  trashVertices.push_back( vertexLower.getID() );
               }
            }
            else
            {
               // Nevertheless we need to store ID of the vertex also in the lower slice array.
               slicesLower[slice].push_back( vertexUpper.getID() );
            }
         }
      }
   }

   // Connecting the stored IDs to construct the cells
   for ( uint_t slice = 0; slice < numSlices; slice++ )
   {
      auto nextSlice = slice + 1;
      if ( slice == numSlices - 1 )
      {
         nextSlice = 0;
      }

      for ( uint_t z_idx = 0; z_idx < width - 1; z_idx++ )
      {
         for ( uint_t x_idx = 0; x_idx + z_idx < width - 1; x_idx++ )
         {
            if ( x_idx >= width - cutSide - 1 )
            {
               continue;
            }

            if ( z_idx >= width - cutTopAndBottom - 1 )
            {
               continue;
            }

            std::vector< uint_t > cubeIndicesUpper( 8 );

            // reversing the order here from xzy to zxy
            // this is almost to complicated to explain in a comment :/
            cubeIndicesUpper[0] = slicesUpper[slice][sliceIndex( width, x_idx, z_idx )];
            cubeIndicesUpper[1] = slicesUpper[slice][sliceIndex( width, x_idx, z_idx + 1 )];
            cubeIndicesUpper[2] = slicesUpper[slice][sliceIndex( width, x_idx + 1, z_idx )];
            cubeIndicesUpper[3] = slicesUpper[slice][sliceIndex( width, x_idx + 1, z_idx + 1 )];

            cubeIndicesUpper[4] = slicesUpper[nextSlice][sliceIndex( width, x_idx, z_idx )];
            cubeIndicesUpper[5] = slicesUpper[nextSlice][sliceIndex( width, x_idx, z_idx + 1 )];
            cubeIndicesUpper[6] = slicesUpper[nextSlice][sliceIndex( width, x_idx + 1, z_idx )];
            cubeIndicesUpper[7] = slicesUpper[nextSlice][sliceIndex( width, x_idx + 1, z_idx + 1 )];

            auto cellsTopLeft     = cellsBottomLeftPrism( { cubeIndicesUpper[0],
                                                        cubeIndicesUpper[1],
                                                        cubeIndicesUpper[2],
                                                        cubeIndicesUpper[4],
                                                        cubeIndicesUpper[5],
                                                        cubeIndicesUpper[6] } );
            auto cellsBottomRight = cellsTopRightPrism( { cubeIndicesUpper[1],
                                                          cubeIndicesUpper[2],
                                                          cubeIndicesUpper[3],
                                                          cubeIndicesUpper[5],
                                                          cubeIndicesUpper[6],
                                                          cubeIndicesUpper[7] } );

            std::vector< uint_t > cubeIndicesLower( 8 );

            cubeIndicesLower[0] = slicesLower[slice][sliceIndex( width, x_idx, z_idx )];
            cubeIndicesLower[1] = slicesLower[slice][sliceIndex( width, x_idx + 1, z_idx )];
            cubeIndicesLower[2] = slicesLower[slice][sliceIndex( width, x_idx, z_idx + 1 )];
            cubeIndicesLower[3] = slicesLower[slice][sliceIndex( width, x_idx + 1, z_idx + 1 )];

            cubeIndicesLower[4] = slicesLower[nextSlice][sliceIndex( width, x_idx, z_idx )];
            cubeIndicesLower[5] = slicesLower[nextSlice][sliceIndex( width, x_idx + 1, z_idx )];
            cubeIndicesLower[6] = slicesLower[nextSlice][sliceIndex( width, x_idx, z_idx + 1 )];
            cubeIndicesLower[7] = slicesLower[nextSlice][sliceIndex( width, x_idx + 1, z_idx + 1 )];

            auto cellsBottomLeft = cellsBottomLeftPrism( { cubeIndicesLower[0],
                                                           cubeIndicesLower[1],
                                                           cubeIndicesLower[2],
                                                           cubeIndicesLower[4],
                                                           cubeIndicesLower[5],
                                                           cubeIndicesLower[6] } );
            auto cellsTopRight   = cellsTopRightPrism( { cubeIndicesLower[1],
                                                       cubeIndicesLower[2],
                                                       cubeIndicesLower[3],
                                                       cubeIndicesLower[5],
                                                       cubeIndicesLower[6],
                                                       cubeIndicesLower[7] } );

            for ( const auto& cell : cellsTopLeft )
            {
               meshInfo.addCellAndAllEdgesAndFaces( cell );
            }

            for ( const auto& cell : cellsBottomLeft )
            {
               meshInfo.addCellAndAllEdgesAndFaces( cell );
            }

            // The remaining cells are those at the tips of the triangular slice.
            // So we do not add the "full square" but on a single triangle.
            if ( x_idx + z_idx < edges - 1 )
            {
               for ( const auto& cell : cellsBottomRight )
               {
                  meshInfo.addCellAndAllEdgesAndFaces( cell );
               }

               for ( const auto& cell : cellsTopRight )
               {
                  meshInfo.addCellAndAllEdgesAndFaces( cell );
               }
            }
         }
      }
   }

   for ( const auto& trashVertex : trashVertices )
   {
      WALBERLA_ASSERT_GREATER( meshInfo.vertices_.count( trashVertex ), 0 );
      meshInfo.vertices_.erase( trashVertex );
   }

   return meshInfo;
}

} // namespace hyteg

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

#include "hyteg/geometry/Torus.hpp"
#include "hyteg/mesh/MeshInfo.hpp"

namespace hyteg {

using walberla::int_c;
using walberla::real_c;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

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

MeshInfo MeshInfo::meshTorus( uint_t                toroidalResolution,
                              uint_t                poloidalResolution,
                              real_t                radiusOriginToCenterOfTube,
                              std::vector< real_t > tubeLayerRadii,
                              real_t                toroidalStartAngle,
                              real_t                poloidalStartAngle,
                              uint_t                numToroidalSlices )
{
   if ( numToroidalSlices == 0 )
   {
      numToroidalSlices = toroidalResolution;
   }

   MeshInfo meshInfo;

   uint_t id = 0;

   real_t toroidalAngleIncrement = 2 * pi / real_c( toroidalResolution );
   real_t poloidalAngleIncrement = 2 * pi / real_c( poloidalResolution );

   // slices[toroidalSlice][layer][poloidalVertex]
   std::map< uint_t, std::map< uint_t, std::map< uint_t, uint_t > > > slices;
   // std::vector< std::vector< uint_t > > slices( toroidalResolution );

   for ( uint_t toroidalSlice = 0; toroidalSlice < std::min( toroidalResolution, numToroidalSlices + 1 ); toroidalSlice++ )
   {
      real_t toroidalAngle = toroidalStartAngle + real_c( toroidalSlice ) * toroidalAngleIncrement;

      // center vertex
      Point3D coords = torusCoordinates( radiusOriginToCenterOfTube, 0, toroidalAngle, 0 );

      auto vertex                 = MeshInfo::Vertex( id++, coords, 0 );
      slices[toroidalSlice][0][0] = vertex.getID();
      meshInfo.addVertex( vertex );
      // WALBERLA_LOG_DEVEL_ON_ROOT( "Added vertex id " << vertex.getID() << " at: toro " << toroidalSlice << ", layer " << 0
      //                                                << " polo " << 0 );

      for ( uint_t layer = 1; layer <= tubeLayerRadii.size(); layer++ )
      {
         for ( uint_t poloidalSlice = 0; poloidalSlice < poloidalResolution; poloidalSlice++ )
         {
            real_t poloidalAngle = poloidalStartAngle + real_c( poloidalSlice ) * poloidalAngleIncrement;

            auto coordsFirstVertex =
                torusCoordinates( radiusOriginToCenterOfTube, tubeLayerRadii[layer - 1], toroidalAngle, poloidalAngle );
            auto coordsNextVertex = torusCoordinates(
                radiusOriginToCenterOfTube, tubeLayerRadii[layer - 1], toroidalAngle, poloidalAngle + poloidalAngleIncrement );

            auto coordsOffsetVector = ( coordsNextVertex - coordsFirstVertex ) * ( real_c( 1 ) / real_c( layer ) );

            uint_t poloidalVertex = poloidalSlice * layer;

            vertex                                       = MeshInfo::Vertex( id++, coordsFirstVertex, 0 );
            slices[toroidalSlice][layer][poloidalVertex] = vertex.getID();
            meshInfo.addVertex( vertex );
            // WALBERLA_LOG_DEVEL_ON_ROOT( "Added vertex id " << vertex.getID() << " at: toro " << toroidalSlice << ", layer "
            //                                                << layer << " polo " << poloidalVertex );

            for ( uint_t l = 1; l < layer; l++ )
            {
               poloidalVertex                               = poloidalSlice * layer + l;
               coords                                       = coordsFirstVertex + real_c( l ) * coordsOffsetVector;
               vertex                                       = MeshInfo::Vertex( id++, coords, 0 );
               slices[toroidalSlice][layer][poloidalVertex] = vertex.getID();
               meshInfo.addVertex( vertex );
               // WALBERLA_LOG_DEVEL_ON_ROOT( " asdfAdded vertex id " << vertex.getID() << " at: toro " << toroidalSlice
               //                                                     << ", layer " << layer << " polo " << poloidalVertex );
            }
         }
      }
   }

   for ( uint_t toroidalSlice = 0; toroidalSlice < numToroidalSlices; toroidalSlice++ )
   {
      for ( uint_t layer = 0; layer < tubeLayerRadii.size(); layer++ )
      {
         for ( uint_t poloidalSlice = 0; poloidalSlice < poloidalResolution; poloidalSlice++ )
         {
            uint_t prismFronSlice = toroidalSlice;
            uint_t prismBackSlice = toroidalSlice + 1 == toroidalResolution ? 0 : toroidalSlice + 1;

            auto numTriangles = 2 * layer + 1;

            for ( uint_t t = 0; t < numTriangles; t++ )
            {
               std::vector< uint_t > prismVerticesFront;
               std::vector< uint_t > prismVerticesBack;

               if ( t % 2 == 0 )
               {
                  // inner node
                  if ( layer == 0 )
                  {
                     // degenerate case where the triangle is connected to the center
                     prismVerticesFront.push_back( slices[prismFronSlice][layer][0] );
                     prismVerticesBack.push_back( slices[prismBackSlice][layer][0] );
                  }
                  else
                  {
                     if ( poloidalSlice == poloidalResolution - 1 && t == numTriangles - 1 )
                     {
                        prismVerticesFront.push_back( slices[prismFronSlice][layer][0] );
                        prismVerticesBack.push_back( slices[prismBackSlice][layer][0] );
                     }
                     else
                     {
                        prismVerticesFront.push_back( slices[prismFronSlice][layer][poloidalSlice * ( layer ) + ( t / 2 )] );
                        prismVerticesBack.push_back( slices[prismBackSlice][layer][poloidalSlice * ( layer ) + ( t / 2 )] );
                     }
                  }

                  // outer nodes
                  if ( poloidalSlice == poloidalResolution - 1 && t == numTriangles - 1 )
                  {
                     prismVerticesFront.push_back( slices[prismFronSlice][layer + 1][poloidalSlice * ( layer + 1 ) + ( t / 2 )] );
                     prismVerticesBack.push_back( slices[prismBackSlice][layer + 1][poloidalSlice * ( layer + 1 ) + ( t / 2 )] );

                     prismVerticesFront.push_back( slices[prismFronSlice][layer + 1][0] );
                     prismVerticesBack.push_back( slices[prismBackSlice][layer + 1][0] );
                  }
                  else
                  {
                     prismVerticesFront.push_back( slices[prismFronSlice][layer + 1][poloidalSlice * ( layer + 1 ) + ( t / 2 )] );
                     prismVerticesBack.push_back( slices[prismBackSlice][layer + 1][poloidalSlice * ( layer + 1 ) + ( t / 2 )] );

                     prismVerticesFront.push_back(
                         slices[prismFronSlice][layer + 1][poloidalSlice * ( layer + 1 ) + ( t / 2 ) + 1] );
                     prismVerticesBack.push_back(
                         slices[prismBackSlice][layer + 1][poloidalSlice * ( layer + 1 ) + ( t / 2 ) + 1] );
                  }

                  auto cells = cellsBottomLeftPrism( { prismVerticesFront[0],
                                                       prismVerticesFront[1],
                                                       prismVerticesFront[2],
                                                       prismVerticesBack[0],
                                                       prismVerticesBack[1],
                                                       prismVerticesBack[2] } );

                  for ( const auto& cell : cells )
                  {
                     meshInfo.addCellAndAllEdgesAndFaces( cell );
                  }
               }
               else
               {
                  // t % 2 == 1
                  WALBERLA_CHECK_GREATER( layer, 0 );

                  // inner nodes
                  if ( poloidalSlice == poloidalResolution - 1 && t == numTriangles - 2 )
                  {
                     prismVerticesFront.push_back( slices[prismFronSlice][layer][poloidalSlice * ( layer ) + ( t / 2 )] );
                     prismVerticesBack.push_back( slices[prismBackSlice][layer][poloidalSlice * ( layer ) + ( t / 2 )] );

                     prismVerticesFront.push_back( slices[prismFronSlice][layer][0] );
                     prismVerticesBack.push_back( slices[prismBackSlice][layer][0] );
                  }
                  else
                  {
                     prismVerticesFront.push_back( slices[prismFronSlice][layer][poloidalSlice * ( layer ) + ( t / 2 )] );
                     prismVerticesBack.push_back( slices[prismBackSlice][layer][poloidalSlice * ( layer ) + ( t / 2 )] );

                     prismVerticesFront.push_back( slices[prismFronSlice][layer][poloidalSlice * ( layer ) + ( t / 2 ) + 1] );
                     prismVerticesBack.push_back( slices[prismBackSlice][layer][poloidalSlice * ( layer ) + ( t / 2 ) + 1] );
                  }

                  // outer node
                  prismVerticesFront.push_back(
                      slices[prismFronSlice][layer + 1][poloidalSlice * ( layer + 1 ) + ( t / 2 ) + 1] );
                  prismVerticesBack.push_back( slices[prismBackSlice][layer + 1][poloidalSlice * ( layer + 1 ) + ( t / 2 ) + 1] );

                  auto cells = cellsTopRightPrism( { prismVerticesFront[0],
                                                     prismVerticesFront[1],
                                                     prismVerticesFront[2],
                                                     prismVerticesBack[0],
                                                     prismVerticesBack[1],
                                                     prismVerticesBack[2] } );

                  for ( const auto& cell : cells )
                  {
                     meshInfo.addCellAndAllEdgesAndFaces( cell );
                  }
               }
            }
         }
      }
   }

   return meshInfo;
}

} // namespace hyteg

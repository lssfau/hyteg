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

MeshInfo MeshInfo::meshTorus( uint_t                numToroidalSlices,
                              uint_t                numPoloidalSlices,
                              real_t                radiusOriginToCenterOfTube,
                              std::vector< real_t > tubeLayerRadii,
                              real_t                toroidalStartAngle,
                              real_t                poloidalStartAngle )
{
   MeshInfo meshInfo;

   WALBERLA_CHECK_EQUAL( tubeLayerRadii.size(), 1 );

   uint_t id = 0;

   auto toroidalAngleIncrement = 2 * pi / real_c( numToroidalSlices );
   auto poloidalAngleIncrement = 2 * pi / real_c( numPoloidalSlices );

   std::vector< std::vector< uint_t > > slices( numToroidalSlices );

   for ( uint_t toroidalSlice = 0; toroidalSlice < numToroidalSlices; toroidalSlice++ )
   {
      auto toroidalAngle = toroidalStartAngle + real_c( toroidalSlice ) * toroidalAngleIncrement;

      // center vertex
      Point3D coords = torusCoordinates( radiusOriginToCenterOfTube, 0, toroidalAngle, 0 );

      auto vertex = MeshInfo::Vertex( id++, coords, 0 );
      slices[toroidalSlice].push_back( vertex.getID() );
      meshInfo.addVertex( vertex );

      for ( uint_t layer = 0; layer < tubeLayerRadii.size(); layer++ )
      {
         for ( uint_t poloidalSlice = 0; poloidalSlice < numPoloidalSlices; poloidalSlice++ )
         {
            auto poloidalAngle = poloidalStartAngle + real_c( poloidalSlice ) * poloidalAngleIncrement;

            coords = torusCoordinates( radiusOriginToCenterOfTube, tubeLayerRadii[layer], toroidalAngle, poloidalAngle );

            vertex = MeshInfo::Vertex( id++, coords, 0 );
            slices[toroidalSlice].push_back( vertex.getID() );
            meshInfo.addVertex( vertex );
         }
      }
   }

   for ( uint_t toroidalSlice = 0; toroidalSlice < numToroidalSlices; toroidalSlice++ )
   {
      for ( uint_t layer = 0; layer < tubeLayerRadii.size(); layer++ )
      {
         for ( uint_t poloidalSlice = 0; poloidalSlice < numPoloidalSlices; poloidalSlice++ )
         {
            std::vector< uint_t > prismVerticesFront;
            std::vector< uint_t > prismVerticesBack;

            uint_t prismFronSlice = toroidalSlice;
            uint_t prismBackSlice = toroidalSlice + 1 == numToroidalSlices ? 0 : toroidalSlice + 1;

            // center vertex
            prismVerticesFront.push_back( slices[prismFronSlice][0] );
            prismVerticesBack.push_back( slices[prismBackSlice][0] );

            if ( poloidalSlice == numPoloidalSlices - 1 )
            {
               prismVerticesFront.push_back( slices[prismFronSlice][poloidalSlice + 1] );
               prismVerticesFront.push_back( slices[prismFronSlice][1] );

               prismVerticesBack.push_back( slices[prismBackSlice][poloidalSlice + 1] );
               prismVerticesBack.push_back( slices[prismBackSlice][1] );
            }
            else
            {
               prismVerticesFront.push_back( slices[prismFronSlice][poloidalSlice + 1] );
               prismVerticesFront.push_back( slices[prismFronSlice][poloidalSlice + 2] );

               prismVerticesBack.push_back( slices[prismBackSlice][poloidalSlice + 1] );
               prismVerticesBack.push_back( slices[prismBackSlice][poloidalSlice + 2] );
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
      }
   }

   return meshInfo;
}

} // namespace hyteg

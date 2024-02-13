/*
 * Copyright (c) 2017-2024 Dominik Thoennes, Andreas Burkhart.
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

#include <unordered_map>
#include <utility>

#include "MeshInfo.hpp"

namespace hyteg {

Point3D getMidPoint( MeshInfo::IDType v0, MeshInfo::IDType v1, const MeshInfo& originalMeshInfo )
{
   auto& originalVertices = originalMeshInfo.getVertices();
   return walberla::real_c( 0.5 ) * ( originalVertices.at( v0 ).getCoordinates() + originalVertices.at( v1 ).getCoordinates() );
}

std::pair< MeshInfo::IDType, MeshInfo::IDType > sortedPair( MeshInfo::IDType v0, MeshInfo::IDType v1 )
{
   return ( v0 <= v1 ) ? std::pair{ v0, v1 } : std::pair{ v1, v0 };
}

MeshInfo MeshInfo::refinedCoarseMesh( const MeshInfo& originalMesh, uint_t refinementSteps )
{
   // lambda function to do one refinement step
   auto refineOnce = []( const MeshInfo& oldMesh ) {
      MeshInfo newMesh = MeshInfo::emptyMeshInfo();

      // find start for new vertex IDs
      uint_t maxVertexId = 0;
      for ( const auto& vertex : oldMesh.vertices_ )
      {
         if ( vertex.first > maxVertexId )
         {
            maxVertexId = vertex.first;
         }
      }
      uint_t newVertexIdStart = maxVertexId + 1;
      uint_t newVertexIds[6];

      // copy vertices from old mesh
      newMesh.vertices_ = oldMesh.vertices_;
      // map from edge (sorted vertex ids) to vertex at mid-point
      std::map< std::pair< IDType, IDType >, IDType > vertexAtMidPoint;

      for ( const auto& cell : oldMesh.cells_ )
      {
         // create new vertices at the middle points ///
         // numbering corresponds to the edge numbering in the macro tet
         std::vector< IDType > oldVertexIds = cell.second.getVertices();

         std::array< std::pair< IDType, IDType >, 6 > edges = { sortedPair( oldVertexIds[0], oldVertexIds[1] ),
                                                                sortedPair( oldVertexIds[0], oldVertexIds[2] ),
                                                                sortedPair( oldVertexIds[1], oldVertexIds[2] ),
                                                                sortedPair( oldVertexIds[0], oldVertexIds[3] ),
                                                                sortedPair( oldVertexIds[1], oldVertexIds[3] ),
                                                                sortedPair( oldVertexIds[2], oldVertexIds[3] ) };

         std::vector< Point3D > newPoints( 6 );
         for ( uint_t i = 0; i < 6; i++ )
         {
            newPoints[i] = getMidPoint( edges[i].first, edges[i].second, oldMesh );

            // check if newPoint already exists in MeshInfo
            if ( vertexAtMidPoint.count( edges[i] ) != 0 )
            {
               newVertexIds[i] = vertexAtMidPoint[edges[i]];
            }
            else
            {
               newVertexIds[i] = newVertexIdStart;
               newVertexIdStart++;
               MeshInfo::Vertex newVertex{ newVertexIds[i], newPoints[i], 0 };
               newMesh.addVertex( newVertex );
               vertexAtMidPoint[edges[i]] = newVertexIds[i];
            }
         }

         // create new tetrahedrons //
         std::vector< MeshInfo::Cell > newCells( 8 );
         // 4 tets contain the original vertices
         newCells[0] = Cell( { oldVertexIds[0], newVertexIds[0], newVertexIds[1], newVertexIds[3] }, 0 );
         newCells[1] = Cell( { newVertexIds[0], oldVertexIds[1], newVertexIds[2], newVertexIds[4] }, 0 );
         newCells[2] = Cell( { newVertexIds[1], newVertexIds[2], oldVertexIds[2], newVertexIds[5] }, 0 );
         newCells[3] = Cell( { newVertexIds[3], newVertexIds[4], newVertexIds[5], oldVertexIds[3] }, 0 );

         // 4 tets are construced from new points
         newCells[4] = Cell( { newVertexIds[0], newVertexIds[1], newVertexIds[2], newVertexIds[4] }, 0 );
         newCells[5] = Cell( { newVertexIds[0], newVertexIds[1], newVertexIds[3], newVertexIds[4] }, 0 );
         newCells[6] = Cell( { newVertexIds[1], newVertexIds[2], newVertexIds[4], newVertexIds[5] }, 0 );
         newCells[7] = Cell( { newVertexIds[1], newVertexIds[3], newVertexIds[4], newVertexIds[5] }, 0 );

         for ( const auto& newCell : newCells )
         {
            newMesh.addCellAndAllEdgesAndFaces( newCell );
         }
      }
      return newMesh;
   };
   // end of lambda function

   MeshInfo newMesh = originalMesh;
   for ( uint_t steps = 1; steps <= refinementSteps; steps++ )
   {
      newMesh = refineOnce( newMesh );
   }
   return newMesh;
}

MeshInfo MeshInfo::refinedCoarseMesh2D( const MeshInfo& originalMesh, uint_t refinementSteps )
{
   // lambda function to do one refinement step
   auto refineOnce = []( const MeshInfo& oldMesh ) {
      MeshInfo newMesh = MeshInfo::emptyMeshInfo();

      // find start for new vertex IDs
      uint_t maxVertexId = 0;
      for ( const auto& vertex : oldMesh.vertices_ )
      {
         if ( vertex.first > maxVertexId )
         {
            maxVertexId = vertex.first;
         }
      }
      uint_t newVertexIdStart = maxVertexId + 1;
      uint_t newVertexIds[3];

      // copy vertices from old mesh
      newMesh.vertices_ = oldMesh.vertices_;
      // map from edge (sorted vertex ids) to vertex at mid-point
      std::map< std::pair< IDType, IDType >, IDType > vertexAtMidPoint;

      for ( const auto& face : oldMesh.faces_ )
      {
         // create new vertices at the middle points ///
         // numbering corresponds to the edge numbering in the macro triangle
         std::vector< IDType > oldVertexIds = face.second.getVertices();

         std::array< std::pair< IDType, IDType >, 3 > edges = { sortedPair( oldVertexIds[0], oldVertexIds[1] ),
                                                                sortedPair( oldVertexIds[0], oldVertexIds[2] ),
                                                                sortedPair( oldVertexIds[1], oldVertexIds[2] ) };

         std::vector< Point3D > newPoints( 3 );
         for ( uint_t i = 0; i < 3; i++ )
         {
            newPoints[i] = getMidPoint( edges[i].first, edges[i].second, oldMesh );

            // check if newPoint already exists in MeshInfo
            if ( vertexAtMidPoint.count( edges[i] ) != 0 )
            {
               newVertexIds[i] = vertexAtMidPoint[edges[i]];
            }
            else
            {
               newVertexIds[i] = newVertexIdStart;
               newVertexIdStart++;
               MeshInfo::Vertex newVertex{ newVertexIds[i], newPoints[i], 0 };
               newMesh.addVertex( newVertex );
               vertexAtMidPoint[edges[i]] = newVertexIds[i];
            }
         }

         // 3 faces contain the original vertices
         newMesh.addFace( Face( { oldVertexIds[0], newVertexIds[0], newVertexIds[1] }, 0 ) );
         newMesh.addEdge( Edge( std::array< IDType, 2 >( { { oldVertexIds[0], newVertexIds[0] } } ), 0 ) );
         newMesh.addEdge( Edge( std::array< IDType, 2 >( { { newVertexIds[0], newVertexIds[1] } } ), 0 ) );
         newMesh.addEdge( Edge( std::array< IDType, 2 >( { { newVertexIds[1], oldVertexIds[0] } } ), 0 ) );

         newMesh.addFace( Face( { newVertexIds[0], oldVertexIds[1], newVertexIds[2] }, 0 ) );
         newMesh.addEdge( Edge( std::array< IDType, 2 >( { { newVertexIds[0], oldVertexIds[1] } } ), 0 ) );
         newMesh.addEdge( Edge( std::array< IDType, 2 >( { { oldVertexIds[1], newVertexIds[2] } } ), 0 ) );
         newMesh.addEdge( Edge( std::array< IDType, 2 >( { { newVertexIds[2], newVertexIds[0] } } ), 0 ) );

         newMesh.addFace( Face( { newVertexIds[1], newVertexIds[2], oldVertexIds[2] }, 0 ) );
         newMesh.addEdge( Edge( std::array< IDType, 2 >( { { newVertexIds[1], newVertexIds[2] } } ), 0 ) );
         newMesh.addEdge( Edge( std::array< IDType, 2 >( { { newVertexIds[2], oldVertexIds[2] } } ), 0 ) );
         newMesh.addEdge( Edge( std::array< IDType, 2 >( { { oldVertexIds[2], newVertexIds[1] } } ), 0 ) );

         // 1 face is constructed from the new vertices
         newMesh.addFace( Face( { newVertexIds[0], newVertexIds[2], newVertexIds[1] }, 0 ) );
         newMesh.addEdge( Edge( std::array< IDType, 2 >( { { newVertexIds[0], newVertexIds[2] } } ), 0 ) );
         newMesh.addEdge( Edge( std::array< IDType, 2 >( { { newVertexIds[2], newVertexIds[1] } } ), 0 ) );
         newMesh.addEdge( Edge( std::array< IDType, 2 >( { { newVertexIds[1], newVertexIds[0] } } ), 0 ) );
      }
      return newMesh;
   };
   // end of lambda function

   MeshInfo newMesh = originalMesh;
   for ( uint_t steps = 1; steps <= refinementSteps; steps++ )
   {
      newMesh = refineOnce( newMesh );
   }
   return newMesh;
}

} // namespace hyteg
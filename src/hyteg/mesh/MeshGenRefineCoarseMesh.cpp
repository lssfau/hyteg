/*
 * Copyright (c) 2017-2021 Dominik Thoennes.
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

#include "core/math/Vector3.h"

#include "MeshInfo.hpp"

#include <unordered_map>

namespace hyteg {

Point3D getMidPoint( MeshInfo::IDType v0, MeshInfo::IDType v1, const MeshInfo& originalMeshInfo )
{
   auto& originalVertices = originalMeshInfo.getVertices();
   return originalVertices.at( v0 ).getCoordinates() +
          0.5 * ( originalVertices.at( v1 ).getCoordinates() - originalVertices.at( v0 ).getCoordinates() );
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
      // create tmpMap for lookup of vertices:
      std::map< std::array< real_t,3 >, IDType > vertexLookUp;
      for ( auto& it : newMesh.vertices_ )
      {
         const auto coords = it.second.getCoordinates();
         vertexLookUp[{coords[0],coords[1],coords[2]}] = it.first;
      }
      for ( const auto& cell : oldMesh.cells_ )
      {
         // create new vertices at the middle points ///
         // numbering corresponds to the edge numbering in the macro tet
         std::vector< IDType > oldVerticesVector = cell.second.getVertices();

         std::vector< Point3D > newPoints( 6 );
         newPoints[0] = getMidPoint( oldVerticesVector[0], oldVerticesVector[1], oldMesh );
         newPoints[1] = getMidPoint( oldVerticesVector[0], oldVerticesVector[2], oldMesh );
         newPoints[2] = getMidPoint( oldVerticesVector[1], oldVerticesVector[2], oldMesh );
         newPoints[3] = getMidPoint( oldVerticesVector[0], oldVerticesVector[3], oldMesh );
         newPoints[4] = getMidPoint( oldVerticesVector[1], oldVerticesVector[3], oldMesh );
         newPoints[5] = getMidPoint( oldVerticesVector[2], oldVerticesVector[3], oldMesh );
         for ( uint_t i = 0; i < 6; i++ )
         {
            // check if newPoint already exists in MeshInfo
            if ( vertexLookUp.count( {newPoints[i][0], newPoints[i][1], newPoints[i][2] } ) != 0 )
            {
               newVertexIds[i] = vertexLookUp[{newPoints[i][0], newPoints[i][1], newPoints[i][2] }];
            }
            else
            {
               newVertexIds[i] = newVertexIdStart;
               newVertexIdStart++;
               auto newVertex = MeshInfo::Vertex( newVertexIds[i], newPoints[i], 0 );
               newMesh.addVertex( newVertex );
               const auto coords = newVertex.getCoordinates();
               vertexLookUp[{coords[0],coords[1],coords[2]}] = newVertexIds[i];
            }
         }

         // create new tetrahedrons //
         std::vector< MeshInfo::Cell > newCells( 8 );
         // 4 tets contain the original vertices
         newCells[0] = Cell( { oldVerticesVector[0], newVertexIds[0], newVertexIds[1], newVertexIds[3] }, 0 );
         newCells[1] = Cell( { newVertexIds[0], oldVerticesVector[1], newVertexIds[2], newVertexIds[4] }, 0 );
         newCells[2] = Cell( { newVertexIds[1], newVertexIds[2], oldVerticesVector[2], newVertexIds[5] }, 0 );
         newCells[3] = Cell( { newVertexIds[3], newVertexIds[4], newVertexIds[5], oldVerticesVector[3] }, 0 );

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

} // namespace hyteg
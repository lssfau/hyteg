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

#include "MeshInfo.hpp"

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
      uint_t newVertexId = maxVertexId + 1;
      // copy vertices from old mesh
      newMesh.vertices_ = oldMesh.vertices_;
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
            newMesh.addVertex( MeshInfo::Vertex( newVertexId + i, newPoints[i], 0 ) );
         }

         // create new tetrahedrons //
         std::vector< MeshInfo::Cell > newCells( 8 );
         // 4 tets contain the original vertices
         newCells[0] = Cell( { oldVerticesVector[0], newVertexId + 0, newVertexId + 1, newVertexId + 3 }, 0 );
         newCells[1] = Cell( { oldVerticesVector[1], newVertexId + 0, newVertexId + 2, newVertexId + 4 }, 0 );
         newCells[2] = Cell( { oldVerticesVector[2], newVertexId + 1, newVertexId + 2, newVertexId + 5 }, 0 );
         newCells[3] = Cell( { oldVerticesVector[3], newVertexId + 3, newVertexId + 4, newVertexId + 5 }, 0 );

         // 4 tets are construced from new points
         newCells[4] = Cell( { newVertexId + 0, newVertexId + 1, newVertexId + 2, newVertexId + 4 }, 0 );
         newCells[5] = Cell( { newVertexId + 0, newVertexId + 1, newVertexId + 3, newVertexId + 4 }, 0 );
         newCells[6] = Cell( { newVertexId + 1, newVertexId + 2, newVertexId + 4, newVertexId + 5 }, 0 );
         newCells[7] = Cell( { newVertexId + 1, newVertexId + 3, newVertexId + 4, newVertexId + 5 }, 0 );

         for ( const auto& newCell : newCells )
         {
            newMesh.addCellAndAllEdgesAndFaces( newCell );
         }

         newVertexId += 6;
      }
      return newMesh;
   };
   // end of lambda function

   MeshInfo newMesh = originalMesh;
   for ( uint_t steps = 0; steps <= refinementSteps; steps++ )
   {
      newMesh = refineOnce( newMesh );
   }
   return newMesh;
}

} // namespace hyteg
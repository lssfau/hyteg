/*
 * Copyright (c) 2017-2021 Benjamin Mann.
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

#include <unordered_map>

namespace hyteg {

MeshInfo MeshInfo::fromCellData( const std::vector<Point3D>& vertices, const std::vector<std::array<uint_t, 4>>& cells )
{
   MeshInfo mesh = MeshInfo::emptyMeshInfo();

   for (uint_t id = 0; id < vertices.size(); ++id)
   {
      mesh.addVertex(MeshInfo::Vertex(id, vertices[id], 0));
   }

   for (auto& cell: cells)
   {
      mesh.addCellAndAllEdgesAndFaces( MeshInfo::Cell({cell[0], cell[1], cell[2], cell[3]}, 0) );
   }

   return mesh;
}

MeshInfo MeshInfo::fromFaceData( const std::vector<Point3D>& vertices, const std::vector<std::array<uint_t, 3>>& faces )
{
   MeshInfo mesh = MeshInfo::emptyMeshInfo();

   for (uint_t id = 0; id < vertices.size(); ++id)
   {
      mesh.addVertex(MeshInfo::Vertex(id, vertices[id], 0));
   }

   for (auto& face: faces)
   {
      mesh.addEdge( MeshInfo::Edge( { face[0], face[1] } , 0 ) );
      mesh.addEdge( MeshInfo::Edge( { face[1], face[2] } , 0 ) );
      mesh.addEdge( MeshInfo::Edge( { face[2], face[0] } , 0 ) );

      mesh.addFace( MeshInfo::Face( { face[0], face[1], face[2] }, 0 ) );
   }

   return mesh;
}

} // namespace hyteg
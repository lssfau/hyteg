/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "core/logging/Logging.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"

namespace hyteg {

MeshInfo MeshInfo::meshUnitSquare(uint_t level) {
  MeshInfo meshInfo;

  uint_t N = walberla::uint_c(std::pow(2 ,level));
  real_t h = walberla::real_c(1.0)/walberla::real_c(N);

  for (uint_t row = 0; row < N+1; ++row) {
    for (uint_t col = 0; col < N+1; ++col) {
      uint_t id = (N+1)*row + col;
      meshInfo.vertices_[id] = MeshInfo::Vertex( id,
                                                 h * walberla::real_c( col ) * Point3D( 1.0, 0.0, 0.0 ) +
                                                     h * walberla::real_c( row ) * Point3D( 0.0, 1.0, 0.0 ),
                                                 0 );
    }
  }

  std::vector<IDType> triangleVertices(3);

  for (uint_t row = 0; row < N; ++row) {
    for (uint_t col = 0; col < N; ++col) {

      // Up triangle
      triangleVertices[0] = (N+1) * row + col;
      triangleVertices[1] = (N+1) * row + (col+1);
      triangleVertices[2] = (N+1) * (row+1) + col;

      meshInfo.addEdge(MeshInfo::Edge(std::array< IDType, 2 >({triangleVertices[0], triangleVertices[1]}), (row == 0) ? 1 : 0));
      meshInfo.addEdge(MeshInfo::Edge(std::array< IDType, 2 >({triangleVertices[1], triangleVertices[2]}), 0));
      meshInfo.addEdge(MeshInfo::Edge(std::array< IDType, 2 >({triangleVertices[2], triangleVertices[0]}), (col == 0) ? 1 : 0));

      meshInfo.addFace(MeshInfo::Face(triangleVertices, 0));

      // Down  triangle
      triangleVertices[0] = (N+1) * (row+1) + (col+1);
      triangleVertices[1] = (N+1) * (row+1) + col;
      triangleVertices[2] = (N+1) * row + (col+1);

      meshInfo.addEdge(MeshInfo::Edge(std::array< IDType, 2 >({triangleVertices[0], triangleVertices[1]}), (row == N-1) ? 1 : 0));
      meshInfo.addEdge(MeshInfo::Edge(std::array< IDType, 2 >({triangleVertices[1], triangleVertices[2]}), 0));
      meshInfo.addEdge(MeshInfo::Edge(std::array< IDType, 2 >({triangleVertices[2], triangleVertices[0]}), (col == N-1) ? 1 : 0));

      meshInfo.addFace(MeshInfo::Face(triangleVertices, 0));
    }
  }

  return meshInfo;
}

}

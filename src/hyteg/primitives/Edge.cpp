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
#include "Edge.hpp"
#include "Vertex.hpp"

#include <core/mpi/MPIManager.h>

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {

using walberla::uint_c;

uint_t Edge::vertex_index(const PrimitiveID& vertex) const
{
  WALBERLA_ASSERT_EQUAL(getNumNeighborVertices(), 2)

  if (vertex == neighborVertices_[0])
  {
    return 0;
  }
  else if (vertex == neighborVertices_[1])
  {
    return 1;
  }
  else
  {
    WALBERLA_ASSERT(false, "Edge::vertex_index: Vertex does not belong to edge")
    return std::numeric_limits<std::size_t>::max();
  }
}

uint_t Edge::face_index(const PrimitiveID& face) const
{
  WALBERLA_ASSERT_GREATER_EQUAL(getNumNeighborFaces(), 1)

  for ( uint_t localFaceID = 0; localFaceID < getNumNeighborFaces(); localFaceID++ )
  {
    if ( face == neighborFaces_[ localFaceID ] )
    {
      return localFaceID;
    }
  }

  WALBERLA_ASSERT(false, "Edge::face_index: Face does not belong to edge")
  return std::numeric_limits<std::size_t>::max();
}

uint_t Edge::cell_index(const PrimitiveID& cell) const
{
  WALBERLA_ASSERT_GREATER_EQUAL(getNumNeighborCells(), 1)

  for ( uint_t localCellID = 0; localCellID < getNumNeighborCells(); localCellID++ )
  {
    if ( cell == neighborCells_[ localCellID ] )
    {
      return localCellID;
    }
  }

  WALBERLA_ASSERT(false, "Edge::cell_index: Cell does not belong to edge")
  return std::numeric_limits<std::size_t>::max();
}

PrimitiveID Edge::get_opposite_vertex(const PrimitiveID& vertex) const
{
  WALBERLA_ASSERT_EQUAL(getNumNeighborVertices(), 2)

  if (vertex == neighborVertices_[0])
  {
    return neighborVertices_[1];
  }
  else if (vertex == neighborVertices_[1])
  {
    return neighborVertices_[0];
  }

  WALBERLA_ABORT("Edge::get_opposite_vertex: Vertex does not belong to edge")
}

bool Edge::opposite_face_exists(const PrimitiveID& face) const
{
  if (face == neighborFaces_[0]) {
     return getNumNeighborFaces() == 2;
  }

  if (getNumNeighborFaces() == 2 && face == neighborFaces_[1])
  {
    return true;
  }

  WALBERLA_ABORT("Edge::opposite_face_exists: Face does not belong to edge")
}

PrimitiveID Edge::get_opposite_face(const PrimitiveID& face) const
{
  if (face == neighborFaces_[0]) {
    if (getNumNeighborFaces() == 2) {
      return neighborFaces_[1];
    } else {
      WALBERLA_ABORT("Edge::get_opposite_face: Requesting face that does not exist")
    }
  }

  if (getNumNeighborFaces() == 2 && face == neighborFaces_[1])
  {
    return neighborFaces_[0];
  }

  WALBERLA_ABORT("Edge::get_opposite_face: Face does not belong to edge")
}


std::ostream& operator<<(std::ostream &os, const hyteg::Edge &edge)
{
  return os << "Edge { id = " << edge.getID() << "; "
            << "neighborVertices_[0] = " << edge.neighborVertices_[0] << "; "
            << "neighborVertices_[1] = " << edge.neighborVertices_[1] << "; }";
}

void Edge::serializeSubclass ( walberla::mpi::SendBuffer & sendBuffer ) const
{
  sendBuffer << coordinates_[0];
  sendBuffer << coordinates_[1];
  sendBuffer << direction_;
  sendBuffer << length_;
  sendBuffer << tangent_;
  sendBuffer << normal2D_;
}

void Edge::deserializeSubclass ( walberla::mpi::RecvBuffer & recvBuffer )
{
  recvBuffer >> coordinates_[0];
  recvBuffer >> coordinates_[1];
  recvBuffer >> direction_;
  recvBuffer >> length_;
  recvBuffer >> tangent_;
  recvBuffer >> normal2D_;
}

}

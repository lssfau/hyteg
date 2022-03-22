/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "Face.hpp"

#include <core/mpi/MPIManager.h>
#include <cstddef>

#include "hyteg/Math.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/types/types.hpp"

#include "Edge.hpp"
#include "Vertex.hpp"

namespace hyteg {
using walberla::uint_c;

uint_t Face::vertex_index(const PrimitiveID& vertex) const
{
  WALBERLA_ASSERT_EQUAL(getNumNeighborVertices(), 3)

  for (size_t i = 0; i < 3; ++i)
  {
    if (vertex == neighborVertices_[i].getID())
    {
      return i;
    }
  }

  WALBERLA_ASSERT(false, "Face::vertex_index: Vertex does not belong to face")
  return std::numeric_limits<std::size_t>::max();
}

uint_t Face::edge_index(const PrimitiveID& edge) const
{
  WALBERLA_ASSERT_EQUAL(getNumNeighborEdges(), 3)

  for (size_t i = 0; i < 3; ++i)
  {
    if (edge.getID() == neighborEdges_[i].getID())
    {
      return i;
    }
  }

  WALBERLA_ASSERT(false, "Face::edge_index: Edge does not belong to face")
  return std::numeric_limits<std::size_t>::max();
}

uint_t Face::cell_index( const PrimitiveID & cell ) const
{
  WALBERLA_ASSERT_LESS_EQUAL( getNumNeighborCells(), 2 )

  for (size_t i = 0; i < 2; ++i)
  {
    if ( cell.getID() == neighborCells_[i].getID() )
    {
      return i;
    }
  }

  WALBERLA_ASSERT(false, "Face::cell_index: Cell does not belong to face")
  return std::numeric_limits<std::size_t>::max();
}

std::vector<PrimitiveID> Face::adjacent_edges(const PrimitiveID& vertex) const
{
  WALBERLA_ASSERT_EQUAL(getNumNeighborEdges(), 3)

  std::vector<PrimitiveID> e;

  if (vertex_index(vertex) == 0) {
    e.push_back(neighborEdges_[0]);
    e.push_back(neighborEdges_[1]);
  } else if (vertex_index(vertex) == 1) {
    e.push_back(neighborEdges_[0]);
    e.push_back(neighborEdges_[2]);
  } else if (vertex_index(vertex) == 2) {
    e.push_back(neighborEdges_[1]);
    e.push_back(neighborEdges_[2]);
  }

  WALBERLA_ASSERT_EQUAL(e.size(), 2)
  return e;
}

PrimitiveID Face::get_vertex_opposite_to_edge(const PrimitiveID& edge) const
{
  WALBERLA_ASSERT_EQUAL(getNumNeighborVertices(), 3)
  WALBERLA_ASSERT_EQUAL(getNumNeighborEdges(), 3)

  if (edge_index(edge) == 0) {
    return neighborVertices_[2];
  } else if (edge_index(edge) == 1) {
    return neighborVertices_[1];
  } else if (edge_index(edge) == 2) {
    return neighborVertices_[0];
  }

  WALBERLA_ABORT("Face::get_vertex_opposite_to_edge: Edge does not belong to face")
}

PrimitiveID Face::getEdgeOppositeToVertex( const PrimitiveID& vertexID ) const
{
   switch ( vertex_index( vertexID ) )
   {
   case 0:
      return neighborEdges().at( 2 );
   case 1:
      return neighborEdges().at( 1 );
   case 2:
      return neighborEdges().at( 0 );
   default:
      WALBERLA_ABORT( "Invalid neighbor vertex ID." );
   }
}

PrimitiveID Face::get_edge_between_vertices(const PrimitiveID& v0, const PrimitiveID& v1) const
{
  std::vector<PrimitiveID> edges_v0 = adjacent_edges(v0);
  std::vector<PrimitiveID> edges_v1 = adjacent_edges(v1);

  if (edges_v0[0] == edges_v1[0]) {
    return edges_v0[0];
  }

  if (edges_v0[0] == edges_v1[1]) {
    return edges_v0[0];
  }

  if (edges_v0[1] == edges_v1[0]) {
    return edges_v0[1];
  }

  if (edges_v0[1] == edges_v1[1]) {
    return edges_v0[1];
  }

  WALBERLA_ABORT("Face::get_edge_between_vertices: Vertex v1 does not belong to face")
}

std::ostream& operator<<(std::ostream &os, const hyteg::Face &face)
{
  return os << "Face { id = " << face.getID().getID() << "; "
            << "neighborEdges_[0] = " << face.neighborEdges_[0].getID() << "; "
            << "neighborEdges_[1] = " << face.neighborEdges_[1].getID() << "; "
            << "neighborEdges_[2] = " << face.neighborEdges_[2].getID() << "; "
            << "}";
}

void Face::serializeSubclass ( walberla::mpi::SendBuffer & sendBuffer ) const
{
  sendBuffer << area;
  sendBuffer << edge_orientation[0];
  sendBuffer << edge_orientation[1];
  sendBuffer << edge_orientation[2];
  sendBuffer << coords[0];
  sendBuffer << coords[1];
  sendBuffer << coords[2];
  sendBuffer << indirectNeighborFaceIDsOverVertices_;
  sendBuffer << indirectNeighborFaceIDsOverEdges_;
}

void Face::deserializeSubclass ( walberla::mpi::RecvBuffer & recvBuffer )
{
  recvBuffer >> area;
  recvBuffer >> edge_orientation[0];
  recvBuffer >> edge_orientation[1];
  recvBuffer >> edge_orientation[2];
  recvBuffer >> coords[0];
  recvBuffer >> coords[1];
  recvBuffer >> coords[2];
  recvBuffer >> indirectNeighborFaceIDsOverVertices_;
  recvBuffer >> indirectNeighborFaceIDsOverEdges_;
}

}




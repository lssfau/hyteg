/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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
#include "Vertex.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

#include <core/mpi/MPIManager.h>

namespace hyteg {

using walberla::uint_c;

Vertex::Vertex( const PrimitiveID & primitiveID, const Point3D & coordinates )
  : Primitive( primitiveID ), dofType_(Inner), coordinates_( coordinates )
{}

uint_t Vertex::edge_index(const PrimitiveID& edge) const
{
  auto it = std::find(neighborEdges().begin(), neighborEdges().end(), edge);
  return uint_c((it - neighborEdges().begin()));
}

uint_t Vertex::face_index(const PrimitiveID& face) const
{
  auto it = std::find(neighborFaces().begin(), neighborFaces().end(), face);
  return uint_c((it - neighborFaces().begin()));
}

uint_t Vertex::cell_index(const PrimitiveID& cell) const
{
  auto it = std::find(neighborCells().begin(), neighborCells().end(), cell);
  return uint_c((it - neighborCells().begin()));
}

std::ostream& operator<<(std::ostream &os, const hyteg::Vertex &vertex)
{
  return os << "Vertex { id = " << vertex.getID() << "; "
            << "coords = [" << vertex.coordinates_[0] << ", " << vertex.coordinates_[1] << ", " << vertex.coordinates_[2] << "]; "
            << "}";
}

void Vertex::serializeSubclass ( walberla::mpi::SendBuffer & sendBuffer ) const
{
  sendBuffer << dofType_;
  sendBuffer << coordinates_;
}

void Vertex::deserializeSubclass ( walberla::mpi::RecvBuffer & recvBuffer )
{
  recvBuffer >> dofType_;
  recvBuffer >> coordinates_;
}

}

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

#include "hyteg/p1functionspace/P1DataHandling.hpp"

namespace hyteg {


std::shared_ptr< VertexP1LocalMatrixMemory > VertexP1LocalMatrixMemoryDataHandling::initialize( const Vertex * const vertex) const
{
  auto vertexP1LocalMatrixMemory = std::make_shared< VertexP1LocalMatrixMemory >();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    vertexP1LocalMatrixMemory->addlevel( level, vertex->getNumNeighborFaces());
  }
  return vertexP1LocalMatrixMemory;
}

std::shared_ptr< EdgeP1LocalMatrixMemory > EdgeP1LocalMatrixMemoryDataHandling::initialize( const Edge * const ) const
{
  auto edgeP1LocalMatrixMemory = std::make_shared< EdgeP1LocalMatrixMemory >();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    edgeP1LocalMatrixMemory->addlevel( level );
  }
  return edgeP1LocalMatrixMemory;
}

std::shared_ptr< FaceP1LocalMatrixMemory > FaceP1LocalMatrixMemoryDataHandling::initialize( const Face * const ) const
{
  auto faceP1LocalMatrixMemory = std::make_shared< FaceP1LocalMatrixMemory >();
  for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
  {
    faceP1LocalMatrixMemory->addlevel( level );
  }
  return faceP1LocalMatrixMemory;
}

std::shared_ptr< FaceP1PolynomialMemory > FaceP1PolynomialMemoryDataHandling::initialize( const Face * const ) const
{
  return std::make_shared< FaceP1PolynomialMemory >( );
}

}

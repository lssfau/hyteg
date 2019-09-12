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

#pragma once

#include "hyteg/primitivedata/PrimitiveDataHandling.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/FunctionMemory.hpp"

namespace hyteg {

class VertexP1LocalMatrixMemoryDataHandling : public OnlyInitializeDataHandling< VertexP1LocalMatrixMemory, Vertex >
{
public:

VertexP1LocalMatrixMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

std::shared_ptr< VertexP1LocalMatrixMemory > initialize( const Vertex * const vertex ) const;

private:

uint_t minLevel_;
uint_t maxLevel_;

};

class EdgeP1LocalMatrixMemoryDataHandling : public OnlyInitializeDataHandling< EdgeP1LocalMatrixMemory, Edge >
{
public:

  EdgeP1LocalMatrixMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  std::shared_ptr< EdgeP1LocalMatrixMemory > initialize( const Edge * const edge ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

class FaceP1LocalMatrixMemoryDataHandling : public OnlyInitializeDataHandling< FaceP1LocalMatrixMemory, Face >
{
public:

  FaceP1LocalMatrixMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel ) : minLevel_( minLevel ), maxLevel_( maxLevel ) {}

  std::shared_ptr< FaceP1LocalMatrixMemory > initialize( const Face * const face ) const;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

class FaceP1PolynomialMemoryDataHandling : public OnlyInitializeDataHandling< FaceP1PolynomialMemory, Face >
{
public:

  FaceP1PolynomialMemoryDataHandling( const uint_t & maxDegree ) : maxDegree_( maxDegree ) {}

  std::shared_ptr< FaceP1PolynomialMemory > initialize( const Face * const face ) const;

private:
  uint_t maxDegree_;

};


}

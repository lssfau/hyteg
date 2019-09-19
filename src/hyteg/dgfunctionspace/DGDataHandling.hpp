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
#pragma once

#include "hyteg/FunctionMemory.hpp"
#include "hyteg/primitives/all.hpp"
#include "DGMemory.hpp"

namespace hyteg {

template< typename ValueType >
class VertexDGFunctionMemoryDataHandling : public FunctionMemoryDataHandling< FunctionMemory < ValueType >, Vertex >
{
public:

  VertexDGFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel )
    : minLevel_( minLevel ),
      maxLevel_( maxLevel )
  {}

  std::shared_ptr< FunctionMemory< ValueType > > initialize( const Vertex * const vertex ) const override;

private:

uint_t minLevel_;
uint_t maxLevel_;

};


template< typename ValueType >
class EdgeDGFunctionMemoryDataHandling : public FunctionMemoryDataHandling< FunctionMemory < ValueType >, Edge >
{
public:

  EdgeDGFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel )
    : minLevel_( minLevel ),
      maxLevel_( maxLevel )
  {}

  std::shared_ptr< FunctionMemory< ValueType > > initialize( const Edge * const vertex ) const override;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

template< typename ValueType >
class FaceDGFunctionMemoryDataHandling : public FunctionMemoryDataHandling< FunctionMemory < ValueType >, Face >
{
public:

  FaceDGFunctionMemoryDataHandling( const uint_t & minLevel, const uint_t & maxLevel )
    : minLevel_( minLevel ),
      maxLevel_( maxLevel )
  {}

  std::shared_ptr<FunctionMemory<ValueType>> initialize( const Face *const face) const override;

private:

  uint_t minLevel_;
  uint_t maxLevel_;

};

template< typename ValueType >
std::shared_ptr< FunctionMemory< ValueType > > VertexDGFunctionMemoryDataHandling< ValueType >::initialize( const Vertex * const vertex ) const
{
  return std::make_shared< FunctionMemory< ValueType > >( DGVertexFunctionMemorySize, *vertex, minLevel_, maxLevel_ );
}

template< typename ValueType >
std::shared_ptr< FunctionMemory< ValueType > > EdgeDGFunctionMemoryDataHandling< ValueType >::initialize( const Edge * const edge ) const
{
  return std::make_shared< FunctionMemory< ValueType > >( DGEdgeFunctionMemorySize, *edge, minLevel_, maxLevel_ );
}

template< typename ValueType >
std::shared_ptr< FunctionMemory< ValueType > > FaceDGFunctionMemoryDataHandling< ValueType >::initialize( const Face * const face ) const
{
  return std::make_shared< FunctionMemory< ValueType > >( DGFaceFunctionMemorySize, *face, minLevel_, maxLevel_ );
}


}

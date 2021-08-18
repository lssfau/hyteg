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

#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/primitives/all.hpp"

#include "FaceDoFMemory.hpp"

namespace hyteg {

template < typename ValueType >
class VertexFaceDoFMemoryDataHandling : public FunctionMemoryDataHandling< FunctionMemory< ValueType >, Vertex >
{
 public:
   VertexFaceDoFMemoryDataHandling( const uint_t& minLevel, const uint_t& maxLevel )
   : minLevel_( minLevel )
   , maxLevel_( maxLevel )
   {}

   std::shared_ptr< FunctionMemory< ValueType > > initialize( const Vertex* const vertex ) const override;

 private:
   uint_t minLevel_;
   uint_t maxLevel_;
};

template < typename ValueType >
class EdgeFaceDoFMemoryDataHandling : public FunctionMemoryDataHandling< FunctionMemory< ValueType >, Edge >
{
 public:
   EdgeFaceDoFMemoryDataHandling( const uint_t& minLevel, const uint_t& maxLevel )
   : minLevel_( minLevel )
   , maxLevel_( maxLevel )
   {}

   std::shared_ptr< FunctionMemory< ValueType > > initialize( const Edge* const vertex ) const override;

 private:
   uint_t minLevel_;
   uint_t maxLevel_;
};

template < typename ValueType >
class FaceFaceDoFMemoryDataHandling : public FunctionMemoryDataHandling< FunctionMemory< ValueType >, Face >
{
 public:
   FaceFaceDoFMemoryDataHandling( const uint_t& minLevel, const uint_t& maxLevel )
   : minLevel_( minLevel )
   , maxLevel_( maxLevel )
   {}

   std::shared_ptr< FunctionMemory< ValueType > > initialize( const Face* const face ) const override;

 private:
   uint_t minLevel_;
   uint_t maxLevel_;
};

template < typename ValueType >
std::shared_ptr< FunctionMemory< ValueType > >
    VertexFaceDoFMemoryDataHandling< ValueType >::initialize( const Vertex* const vertex ) const
{
   return std::make_shared< FunctionMemory< ValueType > >( faceDoFMacroVertexFunctionMemorySize, *vertex, minLevel_, maxLevel_ );
}

template < typename ValueType >
std::shared_ptr< FunctionMemory< ValueType > >
    EdgeFaceDoFMemoryDataHandling< ValueType >::initialize( const Edge* const edge ) const
{
   return std::make_shared< FunctionMemory< ValueType > >( faceDoFMacroEdgeFunctionMemorySize, *edge, minLevel_, maxLevel_ );
}

template < typename ValueType >
std::shared_ptr< FunctionMemory< ValueType > >
    FaceFaceDoFMemoryDataHandling< ValueType >::initialize( const Face* const face ) const
{
   return std::make_shared< FunctionMemory< ValueType > >( faceDoFMacroFaceFunctionMemorySize, *face, minLevel_, maxLevel_ );
}

} // namespace hyteg

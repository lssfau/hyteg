/*
 * Copyright (c) 2025 Benjamin Mann
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

// Contains various conversions between local and global node numbering for P1 elements.

#pragma once

#include <hyteg/Levelinfo.hpp>
#include <hyteg/polynomial/new/data.hpp>
#include <hyteg/volumedofspace/CellDoFIndexing.hpp>
#include <hyteg/volumedofspace/FaceDoFIndexing.hpp>

namespace hyteg {
namespace p1 {

/* convert local tet indices to local cube indices
   indices are encodes as i = x + 2y + 4z = x + y<<1 + z<<2
   where x,y,z ∈ {0,1}
*/
static constexpr surrogate::ElementTypeWiseData< std::array< uint_t, 4 >, 3 > cubeIndicesFromCellIndices{
    std::array< uint_t, 4 >{ 0b000, 0b001, 0b010, 0b100 }, // WHITE UP
    std::array< uint_t, 4 >{ 0b001, 0b011, 0b010, 0b101 }, // BLUE UP
    std::array< uint_t, 4 >{ 0b001, 0b010, 0b101, 0b100 }, // GREEN UP
    std::array< uint_t, 4 >{ 0b011, 0b111, 0b110, 0b101 }, // WHITE DOWN
    std::array< uint_t, 4 >{ 0b101, 0b110, 0b100, 0b010 }, // BLUE DOWN
    std::array< uint_t, 4 >{ 0b010, 0b011, 0b101, 0b110 }  // GREEN DOWN
};

/**
 * @brief Computes the global vertices of the micro-cube containing the six tets corresponding to the given index
 * @param lvl The refinement level of the grid. Determines the number of micro-vertices per edge.
 * @param microElement The micro-element specified by its x, y, and z indices.
 * @param globalDofIndices Array to store the global indices of the cube's vertices.
 */
static inline void
    getGlobalCubeIndices3D( const uint_t lvl, const indexing::Index& microElement, std::array< uint_t, 8 >& globalDofIndices )
{
   const auto x = microElement.x();
   const auto y = microElement.y();
   const auto z = microElement.z();
   const auto n = levelinfo::num_microvertices_per_edge( lvl );

   // local indices of micro cube
   // indices are encodes as i = dx + 2dy + 4dz = dx + dy<<1 + dz<<2 where dx,dy ∈ {0,1}
   globalDofIndices[0b000] = indexing::macroCellIndex( n, x + 0, y + 0, z + 0 );
   globalDofIndices[0b001] = indexing::macroCellIndex( n, x + 1, y + 0, z + 0 );
   globalDofIndices[0b010] = indexing::macroCellIndex( n, x + 0, y + 1, z + 0 );
   globalDofIndices[0b011] = indexing::macroCellIndex( n, x + 1, y + 1, z + 0 );
   globalDofIndices[0b100] = indexing::macroCellIndex( n, x + 0, y + 0, z + 1 );
   globalDofIndices[0b101] = indexing::macroCellIndex( n, x + 1, y + 0, z + 1 );
   globalDofIndices[0b110] = indexing::macroCellIndex( n, x + 0, y + 1, z + 1 );
   globalDofIndices[0b111] = indexing::macroCellIndex( n, x + 1, y + 1, z + 1 );
}

namespace stencil {
// number of stencil entries
static constexpr inline size_t stencilSize( uint8_t dim )
{
   return ( dim == 3 ) ? 15 : ( dim == 2 ) ? 7 : 3;
}

// container for stencils
template < uint8_t DIM, typename dType = real_t >
using StencilData = std::array< dType, stencilSize( DIM ) >;

// due to the strange ordering of the old stencilDirection enum, we redefine them in
// a way s.th. constant sized arrays with 15 (7 in 2D) entries can be used for stencils.
enum Dir
{
   C,
   W,
   E,
   S,
   SE,
   NW,
   N,
   TC,
   TW,
   TS,
   TSE,
   BC,
   BE,
   BNW,
   BN
};

// convert stencilDirection to stencil::Dir
static constexpr Dir conversion( stencilDirection sd )
{
   switch ( sd )
   {
   case stencilDirection::VERTEX_C:
      return Dir::C;
   case stencilDirection::VERTEX_W:
      return Dir::W;
   case stencilDirection::VERTEX_E:
      return Dir::E;
   case stencilDirection::VERTEX_S:
      return Dir::S;
   case stencilDirection::VERTEX_SE:
      return Dir::SE;
   case stencilDirection::VERTEX_NW:
      return Dir::NW;
   case stencilDirection::VERTEX_N:
      return Dir::N;
   case stencilDirection::VERTEX_TC:
      return Dir::TC;
   case stencilDirection::VERTEX_TW:
      return Dir::TW;
   case stencilDirection::VERTEX_TS:
      return Dir::TS;
   case stencilDirection::VERTEX_TSE:
      return Dir::TSE;
   case stencilDirection::VERTEX_BC:
      return Dir::BC;
   case stencilDirection::VERTEX_BE:
      return Dir::BE;
   case stencilDirection::VERTEX_BNW:
      return Dir::BNW;
   case stencilDirection::VERTEX_BN:
      return Dir::BN;
   default:
      return Dir::C;
   }
}

// convert stencil::Dir to stencilDirection
static constexpr StencilData< 3, stencilDirection > backConversion = {
    stencilDirection::VERTEX_C,
    stencilDirection::VERTEX_W,
    stencilDirection::VERTEX_E,
    stencilDirection::VERTEX_S,
    stencilDirection::VERTEX_SE,
    stencilDirection::VERTEX_NW,
    stencilDirection::VERTEX_N,
    stencilDirection::VERTEX_TC,
    stencilDirection::VERTEX_TW,
    stencilDirection::VERTEX_TS,
    stencilDirection::VERTEX_TSE,
    stencilDirection::VERTEX_BC,
    stencilDirection::VERTEX_BE,
    stencilDirection::VERTEX_BNW,
    stencilDirection::VERTEX_BN //
};

static const StencilData< 3, indexing::Index > offset = {
    indexing::Index{ 0, 0, 0 },   // C
    indexing::Index{ -1, 0, 0 },  // W
    indexing::Index{ 1, 0, 0 },   // E
    indexing::Index{ 0, -1, 0 },  // S
    indexing::Index{ 1, -1, 0 },  // SE
    indexing::Index{ -1, 1, 0 },  // NW
    indexing::Index{ 0, 1, 0 },   // N
    indexing::Index{ 0, 0, 1 },   // TC
    indexing::Index{ -1, 0, 1 },  // TW
    indexing::Index{ 0, -1, 1 },  // TS
    indexing::Index{ 1, -1, 1 },  // TSE
    indexing::Index{ 0, 0, -1 },  // BC
    indexing::Index{ 1, 0, -1 },  // BE
    indexing::Index{ -1, 1, -1 }, // BNW
    indexing::Index{ 0, 1, -1 }   // BN
};

static const StencilData< 3, std::string > dirName =
    { "C", "W", "E", "S", "SE", "NW", "N", "TC", "TW", "TS", "TSE", "BC", "BE", "BNW", "BN" };

// get global indices of adjacent vertices
template < uint8_t DIM >
static void getGlobalIndices( const uint_t& level, const indexing::Index& vtx, StencilData< DIM, uint_t >& globalDofIndices )
{
   const auto n = levelinfo::num_microvertices_per_edge( level );

   for ( uint_t k = 0; k < stencilSize( DIM ); ++k )
   {
      const auto& d = offset[k];
      if constexpr ( DIM == 2 )
      {
         globalDofIndices[k] = indexing::macroFaceIndex( n, vtx.x() + d.x(), vtx.y() + d.y() );
      }
      else
      {
         globalDofIndices[k] = indexing::macroCellIndex( n, vtx.x() + d.x(), vtx.y() + d.y(), vtx.z() + d.z() );
      }
   }
}

} // namespace stencil
} // namespace p1
} // namespace hyteg

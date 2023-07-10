/*
 * Copyright (c) 2017-2022 Dominik Thoennes, Nils Kohl.
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

#include <map>

#include "core/Abort.h"
#include "core/DataTypes.h"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"

namespace hyteg {
namespace celldof {

using indexing::Index;
using walberla::uint_t;

enum class CellType : uint_t
{
   WHITE_UP,
   BLUE_UP,
   GREEN_UP,
   WHITE_DOWN,
   BLUE_DOWN,
   GREEN_DOWN
};

const std::map< CellType, std::string > CellTypeToStr = { { CellType::WHITE_UP, "WHITE_UP" },
                                                          { CellType::BLUE_UP, "BLUE_UP" },
                                                          { CellType::GREEN_UP, "GREEN_UP" },
                                                          { CellType::WHITE_DOWN, "WHITE_DOWN" },
                                                          { CellType::BLUE_DOWN, "BLUE_DOWN" },
                                                          { CellType::GREEN_DOWN, "GREEN_DOWN" } };

const std::array< CellType, 6 > allCellTypes = { { CellType::WHITE_UP,
                                                   CellType::BLUE_UP,
                                                   CellType::GREEN_UP,
                                                   CellType::WHITE_DOWN,
                                                   CellType::BLUE_DOWN,
                                                   CellType::GREEN_DOWN } };

namespace macrocell {

inline constexpr uint_t numCellsPerRowByType( const uint_t& level, const CellType& cellType )
{
   switch ( cellType )
   {
   case CellType::WHITE_UP:
      return levelinfo::num_microedges_per_edge( level );
   case CellType::BLUE_UP:
      return levelinfo::num_microedges_per_edge( level ) - 1;
   case CellType::GREEN_UP:
      return levelinfo::num_microedges_per_edge( level ) - 1;
   case CellType::WHITE_DOWN:
      return levelinfo::num_microedges_per_edge( level ) - 2;
   case CellType::BLUE_DOWN:
      return levelinfo::num_microedges_per_edge( level ) - 1;
   case CellType::GREEN_DOWN:
      return levelinfo::num_microedges_per_edge( level ) - 1;
   default:
      return std::numeric_limits< uint_t >::max();
   }
}

inline constexpr uint_t numCellsPerRowByTypeFromWidth( const uint_t& width, const CellType& cellType )
{
   switch ( cellType )
   {
   case CellType::WHITE_UP:
      return levelinfo::num_microedges_per_edge_from_width( width );
   case CellType::BLUE_UP:
      return levelinfo::num_microedges_per_edge_from_width( width ) - 1;
   case CellType::GREEN_UP:
      return levelinfo::num_microedges_per_edge_from_width( width ) - 1;
   case CellType::WHITE_DOWN:
      return levelinfo::num_microedges_per_edge_from_width( width ) - 2;
   case CellType::BLUE_DOWN:
      return levelinfo::num_microedges_per_edge_from_width( width ) - 1;
   case CellType::GREEN_DOWN:
      return levelinfo::num_microedges_per_edge_from_width( width ) - 1;
   default:
      return std::numeric_limits< uint_t >::max();
   }
}

inline constexpr uint_t numMicroCellsPerMacroCell( const uint_t& level, const CellType& cellType )
{
   return levelinfo::num_microvertices_per_cell_from_width( numCellsPerRowByType( level, cellType ) );
}

inline constexpr uint_t numMicroCellsPerMacroCellTotal( const uint_t& level )
{
   uint_t totalNumCells = 0;

   totalNumCells += numMicroCellsPerMacroCell( level, CellType::WHITE_UP );
   totalNumCells += numMicroCellsPerMacroCell( level, CellType::BLUE_UP );
   totalNumCells += numMicroCellsPerMacroCell( level, CellType::GREEN_UP );
   totalNumCells += numMicroCellsPerMacroCell( level, CellType::WHITE_DOWN );
   totalNumCells += numMicroCellsPerMacroCell( level, CellType::BLUE_DOWN );
   totalNumCells += numMicroCellsPerMacroCell( level, CellType::GREEN_DOWN );

   return totalNumCells;
}

inline constexpr uint_t index( const uint_t& level, const idx_t& x, const idx_t& y, const idx_t& z, const CellType& cellType )
{
   const auto width = numCellsPerRowByType( level, cellType );
   switch ( cellType )
   {
   case CellType::WHITE_UP:
      return indexing::macroCellIndex( width, x, y, z );
   case CellType::BLUE_UP:
      return numMicroCellsPerMacroCell( level, CellType::WHITE_UP ) + indexing::macroCellIndex( width, x, y, z );
   case CellType::GREEN_UP:
      return numMicroCellsPerMacroCell( level, CellType::WHITE_UP ) + numMicroCellsPerMacroCell( level, CellType::BLUE_UP ) +
             indexing::macroCellIndex( width, x, y, z );
   case CellType::WHITE_DOWN:
      return numMicroCellsPerMacroCell( level, CellType::WHITE_UP ) + numMicroCellsPerMacroCell( level, CellType::BLUE_UP ) +
             numMicroCellsPerMacroCell( level, CellType::GREEN_UP ) + indexing::macroCellIndex( width, x, y, z );
   case CellType::BLUE_DOWN:
      return numMicroCellsPerMacroCell( level, CellType::WHITE_UP ) + numMicroCellsPerMacroCell( level, CellType::BLUE_UP ) +
             numMicroCellsPerMacroCell( level, CellType::GREEN_UP ) + numMicroCellsPerMacroCell( level, CellType::WHITE_DOWN ) +
             indexing::macroCellIndex( width, x, y, z );
   case CellType::GREEN_DOWN:
      return numMicroCellsPerMacroCell( level, CellType::WHITE_UP ) + numMicroCellsPerMacroCell( level, CellType::BLUE_UP ) +
             numMicroCellsPerMacroCell( level, CellType::GREEN_UP ) + numMicroCellsPerMacroCell( level, CellType::WHITE_DOWN ) +
             numMicroCellsPerMacroCell( level, CellType::BLUE_DOWN ) + indexing::macroCellIndex( width, x, y, z );
   default:
      return std::numeric_limits< uint_t >::max();
   }
}

/// Returns an array of the four logical micro-vertex-indices that span the micro-cell of the given indices and cell type.
/// \sa hyteg::n1e1::macrocell::getMicroVerticesFromMicroCell
inline std::array< Index, 4 > getMicroVerticesFromMicroCell( const Index& microCellIndex, const CellType& microCellType )
{
   const idx_t cellX = microCellIndex.x();
   const idx_t cellY = microCellIndex.y();
   const idx_t cellZ = microCellIndex.z();

   switch ( microCellType )
   {
   case CellType::WHITE_UP:
      return std::array< Index, 4 >( { { Index( cellX, cellY, cellZ ),
                                         Index( cellX + 1, cellY, cellZ ),
                                         Index( cellX, cellY + 1, cellZ ),
                                         Index( cellX, cellY, cellZ + 1 ) } } );
   case CellType::BLUE_UP:
      return std::array< Index, 4 >( { { Index( cellX + 1, cellY, cellZ ),
                                         Index( cellX + 1, cellY + 1, cellZ ),
                                         Index( cellX, cellY + 1, cellZ ),
                                         Index( cellX + 1, cellY, cellZ + 1 ) } } );
   case CellType::GREEN_UP:
      return std::array< Index, 4 >( { { Index( cellX + 1, cellY, cellZ ),
                                         Index( cellX, cellY + 1, cellZ ),
                                         Index( cellX + 1, cellY, cellZ + 1 ),
                                         Index( cellX, cellY, cellZ + 1 ) } } );
   case CellType::WHITE_DOWN:
      return std::array< Index, 4 >( { { Index( cellX + 1, cellY + 1, cellZ ),
                                         Index( cellX + 1, cellY + 1, cellZ + 1 ),
                                         Index( cellX, cellY + 1, cellZ + 1 ),
                                         Index( cellX + 1, cellY, cellZ + 1 ) } } );
   case CellType::BLUE_DOWN:
      return std::array< Index, 4 >( { { Index( cellX + 1, cellY, cellZ + 1 ),
                                         Index( cellX, cellY + 1, cellZ + 1 ),
                                         Index( cellX, cellY, cellZ + 1 ),
                                         Index( cellX, cellY + 1, cellZ ) } } );
   case CellType::GREEN_DOWN:
      return std::array< Index, 4 >( { { Index( cellX, cellY + 1, cellZ ),
                                         Index( cellX + 1, cellY + 1, cellZ ),
                                         Index( cellX + 1, cellY, cellZ + 1 ),
                                         Index( cellX, cellY + 1, cellZ + 1 ) } } );
   default:
      WALBERLA_ABORT( "Not implemented for this cell type." );
      break;
   }
   return std::array< Index, 4 >();
}

/// \brief Given four micro vertex indices, this function returns the corresponding cell idx and cell type.
///
/// Note that the vertices are sorted internally, so that the order in which they are input does not matter.
///
/// \param vertices micro vertex indices of the micro cell of interest
/// \return pair of micro cell index and type
inline std::pair< Index, CellType > getMicroCellFromMicroVertices( const std::array< Index, 4 >& vertices )
{
   auto v = vertices;
   std::sort( v.begin(), v.end() );

   Index v0( v[0] );
   for ( uint_t i = 0; i < 4; i++ )
   {
      v[i] -= v0;
   }

   std::pair< Index, CellType > ret;

   std::array< Index, 4 > white_up   = { { Index( 0, 0, 0 ), Index( 1, 0, 0 ), Index( 0, 1, 0 ), Index( 0, 0, 1 ) } };
   std::array< Index, 4 > white_down = { { Index( 0, 0, 0 ), Index( 0, -1, 1 ), Index( -1, 0, 1 ), Index( 0, 0, 1 ) } };
   std::array< Index, 4 > blue_up    = { { Index( 0, 0, 0 ), Index( -1, 1, 0 ), Index( 0, 1, 0 ), Index( 0, 0, 1 ) } };
   std::array< Index, 4 > blue_down  = { { Index( 0, 0, 0 ), Index( 0, -1, 1 ), Index( 1, -1, 1 ), Index( 0, 0, 1 ) } };
   std::array< Index, 4 > green_up   = { { Index( 0, 0, 0 ), Index( -1, 1, 0 ), Index( -1, 0, 1 ), Index( 0, 0, 1 ) } };
   std::array< Index, 4 > green_down = { { Index( 0, 0, 0 ), Index( 1, 0, 0 ), Index( 1, -1, 1 ), Index( 0, 0, 1 ) } };

   if ( v == white_up )
   {
      ret = { v0, CellType::WHITE_UP };
   }
   else if ( v == white_down )
   {
      ret = { v0 - Index( 1, 1, 0 ), CellType::WHITE_DOWN };
   }
   else if ( v == blue_up )
   {
      ret = { v0 - Index( 1, 0, 0 ), CellType::BLUE_UP };
   }
   else if ( v == blue_down )
   {
      ret = { v0 - Index( 0, 1, 0 ), CellType::BLUE_DOWN };
   }
   else if ( v == green_up )
   {
      ret = { v0 - Index( 1, 0, 0 ), CellType::GREEN_UP };
   }
   else if ( v == green_down )
   {
      ret = { v0 - Index( 0, 1, 0 ), CellType::GREEN_DOWN };
   }
   else
   {
      WALBERLA_ABORT( "This is not a valid micro tet: \n" << v[0] << "\n" << v[1] << "\n" << v[2] << "\n" << v[3] );
   }

   return ret;
}

// Iterators

class Iterator : public indexing::CellIterator
{
 public:
   Iterator( const uint_t& level, const CellType& cellType, const uint_t& offsetToCenter = 0 )
   : CellIterator( numCellsPerRowByType( level, cellType ), offsetToCenter )
   {}
};

} // namespace macrocell
} // namespace celldof
} // namespace hyteg

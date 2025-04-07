/*
 * Copyright (c) 2022 Daniel Bauer.
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

#include "P1toN1E1Gradient.hpp"

#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/primitives/Edge.hpp"
#include "hyteg/primitives/Face.hpp"
#include "hyteg/primitives/Cell.hpp"

namespace hyteg {
namespace n1e1 {

void gradientMacroEdge( const real_t* src, real_t* dst, const uint_t& lvl );
void gradientMacroFace( const real_t* src, real_t* dst, const uint_t& lvl );
void gradientMacroCell( const real_t* src, real_t* dst, const uint_t& lvl );

void P1toN1E1Gradient( const P1Function< real_t >&         src,
                       const N1E1VectorFunction< real_t >& dst,
                       const uint_t&                       lvl,
                       const DoFType&                      flag )
{
   src.communicate< Vertex, Edge >( lvl );
   src.communicate< Edge, Face >( lvl );
   src.communicate< Face, Cell >( lvl );

   for ( const auto& it : dst.getStorage()->getEdges() )
   {
      const Edge& edge = *it.second;

      if ( testFlag( dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         auto srcData = edge.getData( src.getEdgeDataID() )->getPointer( lvl );
         auto dstData = edge.getData( dst.getDoFs()->getEdgeDataID() )->getPointer( lvl );

         gradientMacroEdge( srcData, dstData, lvl );
      }
   }

   for ( const auto& it : dst.getStorage()->getFaces() )
   {
      const Face& face = *it.second;

      if ( testFlag( dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         auto srcData = face.getData( src.getFaceDataID() )->getPointer( lvl );
         auto dstData = face.getData( dst.getDoFs()->getFaceDataID() )->getPointer( lvl );

         gradientMacroFace( srcData, dstData, lvl );
      }
   }

   for ( const auto& it : dst.getStorage()->getCells() )
   {
      const Cell& cell = *it.second;

      if ( testFlag( dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         auto srcData = cell.getData( src.getCellDataID() )->getPointer( lvl );
         auto dstData = cell.getData( dst.getDoFs()->getCellDataID() )->getPointer( lvl );

         gradientMacroCell( srcData, dstData, lvl );
      }
   }
}

void gradientMacroEdge( const real_t* src, real_t* dst, const uint_t& lvl )
{
   const uint_t n = levelinfo::num_microedges_per_edge( lvl );

   for ( idx_t i = 0; i < idx_t( n ); ++i )
   {
      dst[edgedof::macroedge::index( lvl, i )] =
          src[vertexdof::macroedge::index( lvl, i + 1 )] - src[vertexdof::macroedge::index( lvl, i )];
   }
}

void gradientMacroFace( const real_t* src, real_t* dst, const uint_t& lvl )
{
   for ( auto it : edgedof::macroface::Iterator( lvl ) )
   {
      const idx_t x = it.x();
      const idx_t y = it.y();

      // horizontal edges
      dst[edgedof::macroface::horizontalIndex( lvl, x, y )] =
          src[vertexdof::macroface::index( lvl, x + 1, y )] - src[vertexdof::macroface::index( lvl, x, y )];

      // vertical edges
      dst[edgedof::macroface::verticalIndex( lvl, x, y )] =
          src[vertexdof::macroface::index( lvl, x, y + 1 )] - src[vertexdof::macroface::index( lvl, x, y )];

      // diagonal edges
      dst[edgedof::macroface::diagonalIndex( lvl, x, y )] =
          src[vertexdof::macroface::index( lvl, x, y + 1 )] - src[vertexdof::macroface::index( lvl, x + 1, y )];
   }
}

void gradientMacroCell( const real_t* src, real_t* dst, const uint_t& lvl )
{
   using edgedof::macrocell::xIndex;
   using edgedof::macrocell::xyIndex;
   using edgedof::macrocell::xyzIndex;
   using edgedof::macrocell::xzIndex;
   using edgedof::macrocell::yIndex;
   using edgedof::macrocell::yzIndex;
   using edgedof::macrocell::zIndex;
   using vertexdof::macrocell::index;

   for ( auto it : edgedof::macrocell::Iterator( lvl ) )
   {
      const idx_t x = it.x();
      const idx_t y = it.y();
      const idx_t z = it.z();

      // clang-format off
      dst[ xIndex( lvl, x, y, z )] = src[index( lvl, x + 1, y    , z     )] - src[index( lvl, x    , y    , z )];
      dst[ yIndex( lvl, x, y, z )] = src[index( lvl, x    , y + 1, z     )] - src[index( lvl, x    , y    , z )];
      dst[ zIndex( lvl, x, y, z )] = src[index( lvl, x    , y    , z + 1 )] - src[index( lvl, x    , y    , z )];
      dst[xyIndex( lvl, x, y, z )] = src[index( lvl, x    , y + 1, z     )] - src[index( lvl, x + 1, y    , z )];
      dst[xzIndex( lvl, x, y, z )] = src[index( lvl, x    , y    , z + 1 )] - src[index( lvl, x + 1, y    , z )];
      dst[yzIndex( lvl, x, y, z )] = src[index( lvl, x    , y    , z + 1 )] - src[index( lvl, x    , y + 1, z )];
      // clang-format on
   }

   for ( auto it : edgedof::macrocell::IteratorXYZ( lvl ) )
   {
      const idx_t x = it.x();
      const idx_t y = it.y();
      const idx_t z = it.z();

      dst[xyzIndex( lvl, x, y, z )] = src[index( lvl, x + 1, y, z + 1 )] - src[index( lvl, x, y + 1, z )];
   }
}

} // namespace n1e1
} // namespace hyteg

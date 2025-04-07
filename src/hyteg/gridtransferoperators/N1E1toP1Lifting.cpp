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

#include "N1E1toP1Lifting.hpp"

#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/primitives/all.hpp"

namespace hyteg {
namespace n1e1 {

void liftMacroVertex( const Vertex&                       vertex,
                      const N1E1VectorFunction< real_t >& src,
                      const P1Function< real_t >&         dst,
                      const uint_t                        lvl );
void liftMacroEdge( const real_t* src, real_t* dst, const uint_t lvl );
void liftMacroFace( const real_t* src, real_t* dst, const uint_t lvl );
void liftMacroCell( const real_t* src, real_t* dst, const uint_t lvl );

void N1E1toP1Lifting( const N1E1VectorFunction< real_t >& src,
                      const P1Function< real_t >&         dst,
                      const uint_t                        lvl,
                      const DoFType&                      flag )
{
   src.communicate< Edge, Vertex >( lvl );

   for ( const auto& it : src.getStorage()->getVertices() )
   {
      const Vertex& vertex = *it.second;

      if ( testFlag( src.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         liftMacroVertex( vertex, src, dst, lvl );
      }
   }

   for ( const auto& it : src.getStorage()->getEdges() )
   {
      const Edge& edge = *it.second;

      if ( testFlag( src.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         auto srcData = edge.getData( src.getDoFs()->getEdgeDataID() );
         auto dstData = edge.getData( dst.getEdgeDataID() );

         dstData->setToZero( lvl );
         liftMacroEdge( srcData->getPointer( lvl ), dstData->getPointer( lvl ), lvl );
      }
   }

   for ( const auto& it : src.getStorage()->getFaces() )
   {
      const Face& face = *it.second;

      if ( testFlag( src.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         auto srcData = face.getData( src.getDoFs()->getFaceDataID() );
         auto dstData = face.getData( dst.getFaceDataID() );

         dstData->setToZero( lvl );
         liftMacroFace( srcData->getPointer( lvl ), dstData->getPointer( lvl ), lvl );
      }
   }

   for ( const auto& it : src.getStorage()->getCells() )
   {
      const Cell& cell = *it.second;

      if ( testFlag( src.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         auto srcData = cell.getData( src.getDoFs()->getCellDataID() );
         auto dstData = cell.getData( dst.getCellDataID() );

         dstData->setToZero( lvl );
         liftMacroCell( srcData->getPointer( lvl ), dstData->getPointer( lvl ), lvl );
      }
   }

   // Face → Edge comm must come before Cell → Face comm
   // otherwise the Cell → Edge contributions are propagated twice, via the two faces adjacent to the cell and edge
   dst.communicateAdditively< Face, Edge >( lvl, DoFType::All ^ flag, *src.getStorage(), false );
   dst.communicateAdditively< Cell, Face >( lvl, DoFType::All ^ flag, *src.getStorage(), false );
   dst.communicateAdditively< Cell, Edge >( lvl, DoFType::All ^ flag, *src.getStorage(), false );
}

void liftMacroVertex( const Vertex&                       vertex,
                      const N1E1VectorFunction< real_t >& src,
                      const P1Function< real_t >&         dst,
                      const uint_t                        lvl )
{
   const auto  storage  = dst.getStorage();
   const auto& vertexId = vertex.getID();
   auto        dstData  = vertex.getData( dst.getVertexDataID() )->getPointer( lvl );

   dstData[0] = 0.0;

   for ( const auto& edgeId : vertex.neighborEdges() )
   {
      const auto edge    = storage->getEdge( edgeId );
      const auto srcData = vertex.getData( src.getDoFs()->getVertexDataID() )->getPointer( lvl );

      const real_t sign = ( edge->vertex_index( vertexId ) == 0 ) ? -1.0 : 1.0;
      dstData[0] += sign * srcData[vertex.edge_index( edgeId )];
   }
}

void liftMacroEdge( const real_t* src, real_t* dst, const uint_t lvl )
{
   const uint_t n = levelinfo::num_microedges_per_edge( lvl );

   for ( idx_t i = 0; i < idx_t( n ); ++i )
   {
      dst[vertexdof::macroedge::index( lvl, i + 1 )] += src[edgedof::macroedge::index( lvl, i )];
      dst[vertexdof::macroedge::index( lvl, i )] -= src[edgedof::macroedge::index( lvl, i )];
   }
}

void liftMacroFace( const real_t* src, real_t* dst, const uint_t lvl )
{
   for ( auto it : edgedof::macroface::Iterator( lvl ) )
   {
      const idx_t x = it.x();
      const idx_t y = it.y();

      // horizontal edges
      if ( !edgedof::macroface::isHorizontalEdgeOnBoundary( lvl, it ) )
      {
         dst[vertexdof::macroface::index( lvl, x + 1, y )] += src[edgedof::macroface::horizontalIndex( lvl, x, y )];
         dst[vertexdof::macroface::index( lvl, x, y )] -= src[edgedof::macroface::horizontalIndex( lvl, x, y )];
      }

      // vertical edges
      if ( !edgedof::macroface::isVerticalEdgeOnBoundary( lvl, it ) )
      {
         dst[vertexdof::macroface::index( lvl, x, y + 1 )] += src[edgedof::macroface::verticalIndex( lvl, x, y )];
         dst[vertexdof::macroface::index( lvl, x, y )] -= src[edgedof::macroface::verticalIndex( lvl, x, y )];
      }

      // diagonal edges
      if ( !edgedof::macroface::isDiagonalEdgeOnBoundary( lvl, it ) )
      {
         dst[vertexdof::macroface::index( lvl, x, y + 1 )] += src[edgedof::macroface::diagonalIndex( lvl, x, y )];
         dst[vertexdof::macroface::index( lvl, x + 1, y )] -= src[edgedof::macroface::diagonalIndex( lvl, x, y )];
      }
   }
}

void liftMacroCell( const real_t* src, real_t* dst, const uint_t lvl )
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
      if ( edgedof::macrocell::isInnerXEdgeDoF( lvl, it ) )
      {
         dst[index( lvl, x + 1, y    , z     )] += src[ xIndex( lvl, x, y, z )];
         dst[index( lvl, x    , y    , z     )] -= src[ xIndex( lvl, x, y, z )];
      }

      if ( edgedof::macrocell::isInnerYEdgeDoF( lvl, it ) )
      {
         dst[index( lvl, x    , y + 1, z     )] += src[ yIndex( lvl, x, y, z )];
         dst[index( lvl, x    , y    , z     )] -= src[ yIndex( lvl, x, y, z )];
      }

      if ( edgedof::macrocell::isInnerZEdgeDoF( lvl, it ) )
      {
         dst[index( lvl, x    , y    , z + 1 )] += src[ zIndex( lvl, x, y, z )];
         dst[index( lvl, x    , y    , z     )] -= src[ zIndex( lvl, x, y, z )];
      }

      if ( edgedof::macrocell::isInnerXYEdgeDoF( lvl, it ) )
      {
         dst[index( lvl, x    , y + 1, z     )] += src[xyIndex( lvl, x, y, z )];
         dst[index( lvl, x + 1, y    , z     )] -= src[xyIndex( lvl, x, y, z )];
      }

      if ( edgedof::macrocell::isInnerXZEdgeDoF( lvl, it ) )
      {
         dst[index( lvl, x    , y    , z + 1 )] += src[xzIndex( lvl, x, y, z )];
         dst[index( lvl, x + 1, y    , z     )] -= src[xzIndex( lvl, x, y, z )];
      }

      if ( edgedof::macrocell::isInnerYZEdgeDoF( lvl, it ) )
      {
         dst[index( lvl, x    , y    , z + 1 )] += src[yzIndex( lvl, x, y, z )];
         dst[index( lvl, x    , y + 1, z     )] -= src[yzIndex( lvl, x, y, z )];
      }
      // clang-format on
   }

   for ( auto it : edgedof::macrocell::IteratorXYZ( lvl ) )
   {
      const idx_t x = it.x();
      const idx_t y = it.y();
      const idx_t z = it.z();

      dst[index( lvl, x + 1, y, z + 1 )] += src[xyzIndex( lvl, x, y, z )];
      dst[index( lvl, x, y + 1, z )] -= src[xyzIndex( lvl, x, y, z )];
   }
}

} // namespace n1e1
} // namespace hyteg

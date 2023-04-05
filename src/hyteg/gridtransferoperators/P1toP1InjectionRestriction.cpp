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

#include "hyteg/gridtransferoperators/P1toP1InjectionRestriction.hpp"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"

namespace hyteg {

void P1toP1InjectionRestriction::restrict3D( const P1Function< real_t >& function,
                                             const uint_t&               sourceLevel,
                                             const DoFType&              flag ) const
{
   const uint_t destinationLevel  = sourceLevel - 1;
   const auto   storage           = function.getStorage();
   const auto   boundaryCondition = function.getBoundaryCondition();

   function.communicate< Vertex, Edge >( sourceLevel );
   function.communicate< Edge, Face >( sourceLevel );
   function.communicate< Face, Edge >( sourceLevel );
   function.communicate< Edge, Vertex >( sourceLevel );

   for ( const auto& it : storage->getVertices() )
   {
      const Vertex& vertex  = *it.second;
      const auto    srcData = vertex.getData( function.getVertexDataID() )->getPointer( sourceLevel );
      auto          dstData = vertex.getData( function.getVertexDataID() )->getPointer( destinationLevel );

      if ( testFlag( boundaryCondition.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         restrictMacroVertex( srcData, dstData, sourceLevel );
      }
   }

   for ( const auto& it : storage->getEdges() )
   {
      const Edge& edge    = *it.second;
      const auto  srcData = edge.getData( function.getEdgeDataID() )->getPointer( sourceLevel );
      auto        dstData = edge.getData( function.getEdgeDataID() )->getPointer( destinationLevel );

      if ( testFlag( boundaryCondition.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         restrictMacroEdge( srcData, dstData, sourceLevel );
      }
   }

   for ( const auto& it : storage->getFaces() )
   {
      const Face& face    = *it.second;
      const auto  srcData = face.getData( function.getFaceDataID() )->getPointer( sourceLevel );
      auto        dstData = face.getData( function.getFaceDataID() )->getPointer( destinationLevel );

      if ( testFlag( boundaryCondition.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         restrictMacroFace( srcData, dstData, sourceLevel );
      }
   }

   for ( const auto& it : storage->getCells() )
   {
      const Cell& cell    = *it.second;
      const auto  srcData = cell.getData( function.getCellDataID() )->getPointer( sourceLevel );
      auto        dstData = cell.getData( function.getCellDataID() )->getPointer( destinationLevel );

      if ( testFlag( boundaryCondition.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         restrictMacroCell( srcData, dstData, sourceLevel );
      }
   }
}

void P1toP1InjectionRestriction::restrictMacroVertex( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const
{
   WALBERLA_UNUSED( sourceLevel );
   dst[0] = src[0];
}

void P1toP1InjectionRestriction::restrictMacroEdge( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const
{
   size_t rowsize_c = levelinfo::num_microvertices_per_edge( sourceLevel - 1 );

   idx_t i_c;
   for ( i_c = 1; i_c < idx_t( rowsize_c ) - 1; ++i_c )
   {
      dst[vertexdof::macroedge::indexFromVertex( sourceLevel - 1, i_c, stencilDirection::VERTEX_C )] =
          real_c( 1.0 ) * src[vertexdof::macroedge::indexFromVertex( sourceLevel, 2 * i_c, stencilDirection::VERTEX_C )];
   }
}

void P1toP1InjectionRestriction::restrictMacroFace( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const
{
   uint_t N_c   = levelinfo::num_microvertices_per_edge( sourceLevel - 1 );
   uint_t N_c_i = N_c;

   real_t tmp;

   for ( idx_t j = 1; j < idx_t( N_c ) - 2; ++j )
   {
      for ( idx_t i = 1; i < idx_t( N_c_i ) - 2; ++i )
      {
         tmp = src[vertexdof::macroface::indexFromVertex( sourceLevel, 2 * i, 2 * j, stencilDirection::VERTEX_C )];

         dst[vertexdof::macroface::indexFromVertex( sourceLevel - 1, i, j, stencilDirection::VERTEX_C )] = tmp;
      }

      --N_c_i;
   }
}

void P1toP1InjectionRestriction::restrictMacroCell( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const
{
   for ( const auto& it : vertexdof::macrocell::Iterator( sourceLevel - 1, 1 ) )
   {
      dst[vertexdof::macrocell::index( sourceLevel - 1, it.x(), it.y(), it.z() )] =
          src[vertexdof::macrocell::index( sourceLevel, 2 * it.x(), 2 * it.y(), 2 * it.z() )];
   }
}

} // namespace hyteg

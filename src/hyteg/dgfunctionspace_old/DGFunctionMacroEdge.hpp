/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/facedofspace_old/FaceDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"


namespace hyteg {
namespace dgfunction {
namespace macroedge {

template < typename ValueType >
inline void upwind( const uint_t&                                                                Level,
                    Edge&                                                                        edge,
                    const std::shared_ptr< PrimitiveStorage >&                                   storage,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&                  srcId,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&                  dstId,
                    const std::array< PrimitiveDataID< FunctionMemory< ValueType >, Edge >, 2 >& velocityIds,
                    UpdateType                                                                   updateType )
{
   auto src = edge.getData( srcId )->getPointer( Level );
   auto dst = edge.getData( dstId )->getPointer( Level );
   auto u   = edge.getData( velocityIds[0] )->getPointer( Level );
   auto v   = edge.getData( velocityIds[1] )->getPointer( Level );

   size_t    rowsize = levelinfo::num_microvertices_per_edge( Level );
   ValueType tmp;
   Point2D   u_0, u_1, u_2;
   real_t    un_0, un_1, un_2;
   real_t    c_up_0, c_up_1, c_up_2;

   // first face (south)
   {
      Face*  face        = storage->getFace( edge.neighborFaces()[0] );
      real_t faceArea    = std::pow( 4.0, -walberla::real_c( Level ) ) * face->area;
      real_t faceAreaInv = 1.0 / faceArea;

      auto oppositeVertex = face->get_vertex_opposite_to_edge( edge.getID() );

      uint_t v0 = face->vertex_index( edge.getVertexID0() );
      uint_t v1 = face->vertex_index( edge.getVertexID1() );
      uint_t v2 = face->vertex_index( oppositeVertex );

      // get edge directions
      auto d0 = ( face->coords[v1] - face->coords[v0] ) / walberla::real_c( rowsize - 1 );
      auto d1 = ( face->coords[v2] - face->coords[v1] ) / walberla::real_c( rowsize - 1 );
      auto d2 = ( face->coords[v0] - face->coords[v2] ) / walberla::real_c( rowsize - 1 );

      // compute edge lengths
      real_t d0Length = d0.norm();
      real_t d1Length = d1.norm();
      real_t d2Length = d2.norm();

      // compute normals
      auto n_0 = d0.normal2D() / d0Length;
      auto n_1 = d1.normal2D() / d1Length;
      auto n_2 = d2.normal2D() / d2Length;

      real_t faceOrientation = math::faceOrientation2D( face->coords[v0], face->coords[v1], face->coords[v2] );
      n_0 *= faceOrientation;
      n_1 *= faceOrientation;
      n_2 *= faceOrientation;

      for ( uint_t i = 1; i < rowsize - 2; ++i )
      {
         u_0[0] = 0.5 * ( u[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] +
                          u[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_E )] );
         u_0[1] = 0.5 * ( v[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] +
                          v[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_E )] );

         u_1[0] = 0.5 * ( u[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_E )] +
                          u[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_SE )] );
         u_1[1] = 0.5 * ( v[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_E )] +
                          v[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_SE )] );

         u_2[0] = 0.5 * ( u[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] +
                          u[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_SE )] );
         u_2[1] = 0.5 * ( v[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] +
                          v[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_SE )] );

         // CENTER <-- CELL_GRAY_SE (i)
         // NORTH  <-- CELL_GRAY_NE (i)
         // WEST   <-- CELL_BLUE_SE (i)
         // EAST   <-- CELL_BLUE_SE (i+1)

         un_0 = d0Length * u_0.dot( n_0 );
         un_1 = d1Length * u_1.dot( n_1 );
         un_2 = d2Length * u_2.dot( n_2 );

         if ( un_0 >= 0 )
         {
            c_up_0 = src[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_SE )];
         }
         else
         {
            c_up_0 = src[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_NE )];
         }

         if ( un_1 >= 0 )
         {
            c_up_1 = src[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_SE )];
         }
         else
         {
            c_up_1 = src[facedof::macroedge::indexFaceFromVertex( Level, i + 1, stencilDirection::CELL_BLUE_SE )];
         }

         if ( un_2 >= 0 )
         {
            c_up_2 = src[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_SE )];
         }
         else
         {
            c_up_2 = src[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_BLUE_SE )];
         }

         tmp = un_0 * c_up_0 + un_1 * c_up_1 + un_2 * c_up_2;
         tmp *= faceAreaInv;

         if ( updateType == Replace )
         {
            dst[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_SE )] = tmp;
         }
         else if ( updateType == Add )
         {
            dst[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_SE )] += tmp;
         }
      }
   }

   // second face (north)
   if ( edge.getNumNeighborFaces() == 2 )
   {
      Face*  face        = storage->getFace( edge.neighborFaces()[1] );
      real_t faceArea    = std::pow( 4.0, -walberla::real_c( Level ) ) * face->area;
      real_t faceAreaInv = 1.0 / faceArea;

      auto oppositeVertex = face->get_vertex_opposite_to_edge( edge.getID() );

      uint_t v0 = face->vertex_index( edge.getVertexID0() );
      uint_t v1 = face->vertex_index( edge.getVertexID1() );
      uint_t v2 = face->vertex_index( oppositeVertex );

      // get edge directions
      auto d0 = ( face->coords[v1] - face->coords[v0] ) / walberla::real_c( rowsize - 1 );
      auto d1 = ( face->coords[v2] - face->coords[v1] ) / walberla::real_c( rowsize - 1 );
      auto d2 = ( face->coords[v0] - face->coords[v2] ) / walberla::real_c( rowsize - 1 );

      // compute edge lengths
      real_t d0Length = d0.norm();
      real_t d1Length = d1.norm();
      real_t d2Length = d2.norm();

      // compute normals
      auto n_0 = d0.normal2D() / d0Length;
      auto n_1 = d1.normal2D() / d1Length;
      auto n_2 = d2.normal2D() / d2Length;

      real_t faceOrientation = math::faceOrientation2D( face->coords[v0], face->coords[v1], face->coords[v2] );
      n_0 *= faceOrientation;
      n_1 *= faceOrientation;
      n_2 *= faceOrientation;

      for ( uint_t i = 1; i < rowsize - 2; ++i )
      {
         u_0[0] = 0.5 * ( u[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] +
                          u[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_E )] );
         u_0[1] = 0.5 * ( v[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] +
                          v[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_E )] );

         u_1[0] = 0.5 * ( u[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_E )] +
                          u[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_N )] );
         u_1[1] = 0.5 * ( v[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_E )] +
                          v[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_N )] );

         u_2[0] = 0.5 * ( u[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] +
                          u[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_N )] );
         u_2[1] = 0.5 * ( v[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] +
                          v[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_N )] );

         // CENTER <-- CELL_GRAY_NE (i)
         // SOUTH  <-- CELL_GRAY_SE (i)
         // WEST   <-- CELL_BLUE_NW (i+1)
         // EAST   <-- CELL_BLUE_NW (i)

         un_0 = d0Length * u_0.dot( n_0 );
         un_1 = d1Length * u_1.dot( n_1 );
         un_2 = d2Length * u_2.dot( n_2 );

         if ( un_0 >= 0 )
         {
            c_up_0 = src[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_NE )];
         }
         else
         {
            c_up_0 = src[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_SE )];
         }

         if ( un_1 >= 0 )
         {
            c_up_1 = src[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_NE )];
         }
         else
         {
            c_up_1 = src[facedof::macroedge::indexFaceFromVertex( Level, i + 1, stencilDirection::CELL_BLUE_NW )];
         }

         if ( un_2 >= 0 )
         {
            c_up_2 = src[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_NE )];
         }
         else
         {
            c_up_2 = src[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_BLUE_NW )];
         }

         tmp = un_0 * c_up_0 + un_1 * c_up_1 + un_2 * c_up_2;
         tmp *= faceAreaInv;

         if ( updateType == Replace )
         {
            dst[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_NE )] = tmp;
         }
         else if ( updateType == Add )
         {
            dst[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_NE )] += tmp;
         }
      }
   }
}

template < typename ValueType >
inline void projectP1( const uint_t&                                               Level,
                       Edge&                                                       edge,
                       const std::shared_ptr< PrimitiveStorage >&                  storage,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& srcId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& dstId,
                       UpdateType                                                  updateType )
{
   auto src = edge.getData( srcId )->getPointer( Level );
   auto dst = edge.getData( dstId )->getPointer( Level );

   size_t    rowsize = levelinfo::num_microvertices_per_edge( Level );
   ValueType tmp;

   // first face (south)
   {
      for ( uint_t i = 1; i < rowsize - 2; ++i )
      {
         tmp = 1.0 / 3.0 *
               ( src[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] +
                 src[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_E )] +
                 src[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_SE )] );

         if ( updateType == Replace )
         {
            dst[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_SE )] = tmp;
         }
         else if ( updateType == Add )
         {
            dst[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_SE )] += tmp;
         }
      }
   }

   // second face (north)
   if ( edge.getNumNeighborFaces() == 2 )
   {
      for ( uint_t i = 1; i < rowsize - 2; ++i )
      {
         tmp = 1.0 / 3.0 *
               ( src[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] +
                 src[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_E )] +
                 src[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_N )] );

         if ( updateType == Replace )
         {
            dst[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_NE )] = tmp;
         }
         else if ( updateType == Add )
         {
            dst[facedof::macroedge::indexFaceFromVertex( Level, i, stencilDirection::CELL_GRAY_NE )] += tmp;
         }
      }
   }
}

} // namespace macroedge
} // namespace dgfunction
} //namespace hyteg

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

// #include "DGFaceIndex.hpp"
// #include "hyteg/bubblefunctionspace/BubbleFaceIndex.hpp"

#include "hyteg/facedofspace/FaceDoFIndexing.hpp"

namespace hyteg {
namespace dgfunction {
namespace macroface {

template < typename ValueType >
inline void upwind( const uint_t&                                                                Level,
                    Face&                                                                        face,
                    const std::shared_ptr< PrimitiveStorage >&                                   storage,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                  srcId,
                    const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                  dstId,
                    const std::array< PrimitiveDataID< FunctionMemory< ValueType >, Face >, 2 >& velocityIds,
                    UpdateType                                                                   updateType )
{
   using namespace vertexdof::macroface;

   size_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   size_t inner_rowsize = rowsize;

   // get memories
   auto src = face.getData( srcId )->getPointer( Level );
   auto dst = face.getData( dstId )->getPointer( Level );
   auto u   = face.getData( velocityIds[0] )->getPointer( Level );
   auto v   = face.getData( velocityIds[1] )->getPointer( Level );

   // get edge directions
   auto d0 = ( face.coords[1] - face.coords[0] ) / walberla::real_c( rowsize - 1 );
   auto d1 = ( face.coords[2] - face.coords[1] ) / walberla::real_c( rowsize - 1 );
   auto d2 = ( face.coords[0] - face.coords[2] ) / walberla::real_c( rowsize - 1 );

   // compute edge lengths
   real_t d0Length = d0.norm();
   real_t d1Length = d1.norm();
   real_t d2Length = d2.norm();

   // compute normals
   auto n_0 = d0.normal2D() / d0Length;
   auto n_1 = d1.normal2D() / d1Length;
   auto n_2 = d2.normal2D() / d2Length;

   real_t faceOrientation = math::faceOrientation2D( face.coords[0], face.coords[1], face.coords[2] );

   // correct normals if all normals point in wrong direction
   n_0 *= faceOrientation;
   n_1 *= faceOrientation;
   n_2 *= faceOrientation;

   real_t faceArea    = std::pow( 4.0, -walberla::real_c( Level ) ) * face.area;
   real_t faceAreaInv = 1.0 / faceArea;

   ValueType tmp;

   Point2D u_0, u_1, u_2;
   real_t  un_0, un_1, un_2;
   real_t  c_up_0, c_up_1, c_up_2;

   for ( size_t j = 1; j < rowsize - 2; ++j )
   {
      for ( size_t i = 1; i < inner_rowsize - 3; ++i )
      {
         // evalate velocities
         u_0[0] = 0.5 * ( u[indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] +
                          u[indexFromVertex( Level, i + 1, j, stencilDirection::VERTEX_C )] );
         u_0[1] = 0.5 * ( v[indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] +
                          v[indexFromVertex( Level, i + 1, j, stencilDirection::VERTEX_C )] );

         u_1[0] = 0.5 * ( u[indexFromVertex( Level, i + 1, j, stencilDirection::VERTEX_C )] +
                          u[indexFromVertex( Level, i, j + 1, stencilDirection::VERTEX_C )] );
         u_1[1] = 0.5 * ( v[indexFromVertex( Level, i + 1, j, stencilDirection::VERTEX_C )] +
                          v[indexFromVertex( Level, i, j + 1, stencilDirection::VERTEX_C )] );

         u_2[0] = 0.5 * ( u[indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] +
                          u[indexFromVertex( Level, i, j + 1, stencilDirection::VERTEX_C )] );
         u_2[1] = 0.5 * ( v[indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] +
                          v[indexFromVertex( Level, i, j + 1, stencilDirection::VERTEX_C )] );

         un_0 = d0Length * u_0.dot( n_0 );
         un_1 = d1Length * u_1.dot( n_1 );
         un_2 = d2Length * u_2.dot( n_2 );

         if ( un_0 >= 0 )
         {
            c_up_0 = src[facedof::macroface::indexFaceFromGrayFace( Level, i, j, stencilDirection::CELL_GRAY_C )];
         }
         else
         {
            c_up_0 = src[facedof::macroface::indexFaceFromGrayFace( Level, i, j, stencilDirection::CELL_BLUE_S )];
         }

         if ( un_1 >= 0 )
         {
            c_up_1 = src[facedof::macroface::indexFaceFromGrayFace( Level, i, j, stencilDirection::CELL_GRAY_C )];
         }
         else
         {
            c_up_1 = src[facedof::macroface::indexFaceFromGrayFace( Level, i, j, stencilDirection::CELL_BLUE_E )];
         }

         if ( un_2 >= 0 )
         {
            c_up_2 = src[facedof::macroface::indexFaceFromGrayFace( Level, i, j, stencilDirection::CELL_GRAY_C )];
         }
         else
         {
            c_up_2 = src[facedof::macroface::indexFaceFromGrayFace( Level, i, j, stencilDirection::CELL_BLUE_W )];
         }

         tmp = un_0 * c_up_0 + un_1 * c_up_1 + un_2 * c_up_2;
         tmp *= faceAreaInv;

         if ( updateType == Replace )
         {
            dst[facedof::macroface::indexFaceFromGrayFace( Level, i, j, stencilDirection::CELL_GRAY_C )] = tmp;
         }
         else if ( updateType == Add )
         {
            dst[facedof::macroface::indexFaceFromGrayFace( Level, i, j, stencilDirection::CELL_GRAY_C )] += tmp;
         }
      }
      --inner_rowsize;
   }

   inner_rowsize = rowsize;

   // flip normals
   n_0 *= -1.0;
   n_1 *= -1.0;
   n_2 *= -1.0;

   for ( size_t j = 0; j < rowsize - 2; ++j )
   {
      for ( size_t i = 0; i < inner_rowsize - 2; ++i )
      {
         // evalate velocities
         u_0[0] = 0.5 * ( u[indexFromVertex( Level, i, j + 1, stencilDirection::VERTEX_C )] +
                          u[indexFromVertex( Level, i + 1, j + 1, stencilDirection::VERTEX_C )] );
         u_0[1] = 0.5 * ( v[indexFromVertex( Level, i, j + 1, stencilDirection::VERTEX_C )] +
                          v[indexFromVertex( Level, i + 1, j + 1, stencilDirection::VERTEX_C )] );

         u_1[0] = 0.5 * ( u[indexFromVertex( Level, i, j + 1, stencilDirection::VERTEX_C )] +
                          u[indexFromVertex( Level, i + 1, j, stencilDirection::VERTEX_C )] );
         u_1[1] = 0.5 * ( v[indexFromVertex( Level, i, j + 1, stencilDirection::VERTEX_C )] +
                          v[indexFromVertex( Level, i + 1, j, stencilDirection::VERTEX_C )] );

         u_2[0] = 0.5 * ( u[indexFromVertex( Level, i + 1, j, stencilDirection::VERTEX_C )] +
                          u[indexFromVertex( Level, i + 1, j + 1, stencilDirection::VERTEX_C )] );
         u_2[1] = 0.5 * ( v[indexFromVertex( Level, i + 1, j, stencilDirection::VERTEX_C )] +
                          v[indexFromVertex( Level, i + 1, j + 1, stencilDirection::VERTEX_C )] );

         un_0 = d0Length * u_0.dot( n_0 );
         un_1 = d1Length * u_1.dot( n_1 );
         un_2 = d2Length * u_2.dot( n_2 );

         if ( un_0 >= 0 )
         {
            c_up_0 = src[facedof::macroface::indexFaceFromBlueFace( Level, i, j, stencilDirection::CELL_BLUE_C )];
         }
         else
         {
            c_up_0 = src[facedof::macroface::indexFaceFromBlueFace( Level, i, j, stencilDirection::CELL_GRAY_N )];
         }

         if ( un_1 >= 0 )
         {
            c_up_1 = src[facedof::macroface::indexFaceFromBlueFace( Level, i, j, stencilDirection::CELL_BLUE_C )];
         }
         else
         {
            c_up_1 = src[facedof::macroface::indexFaceFromBlueFace( Level, i, j, stencilDirection::CELL_GRAY_W )];
         }

         if ( un_2 >= 0 )
         {
            c_up_2 = src[facedof::macroface::indexFaceFromBlueFace( Level, i, j, stencilDirection::CELL_BLUE_C )];
         }
         else
         {
            c_up_2 = src[facedof::macroface::indexFaceFromBlueFace( Level, i, j, stencilDirection::CELL_GRAY_E )];
         }

         tmp = un_0 * c_up_0 + un_1 * c_up_1 + un_2 * c_up_2;
         tmp *= faceAreaInv;

         if ( updateType == Replace )
         {
            dst[facedof::macroface::indexFaceFromBlueFace( Level, i, j, stencilDirection::CELL_BLUE_C )] = tmp;
         }
         else if ( updateType == Add )
         {
            dst[facedof::macroface::indexFaceFromBlueFace( Level, i, j, stencilDirection::CELL_BLUE_C )] += tmp;
         }
      }
      --inner_rowsize;
   }
}

template < typename ValueType >
inline void projectP1( const uint_t&                                               Level,
                       Face&                                                       face,
                       const std::shared_ptr< PrimitiveStorage >&                  storage,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Face >& srcId,
                       const PrimitiveDataID< FunctionMemory< ValueType >, Face >& dstId,
                       UpdateType                                                  updateType )
{
   using namespace vertexdof::macroface;

   size_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   size_t inner_rowsize = rowsize;

   // get memories
   auto src = face.getData( srcId )->getPointer( Level );
   auto dst = face.getData( dstId )->getPointer( Level );

   ValueType tmp;

   for ( size_t j = 1; j < rowsize - 2; ++j )
   {
      for ( size_t i = 1; i < inner_rowsize - 3; ++i )
      {
         // evalate velocities
         tmp = 1.0 / 3.0 *
               ( src[indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] +
                 src[indexFromVertex( Level, i + 1, j, stencilDirection::VERTEX_C )] +
                 src[indexFromVertex( Level, i, j + 1, stencilDirection::VERTEX_C )] );

         if ( updateType == Replace )
         {
            dst[facedof::macroface::indexFaceFromGrayFace( Level, i, j, stencilDirection::CELL_GRAY_C )] = tmp;
         }
         else if ( updateType == Add )
         {
            dst[facedof::macroface::indexFaceFromGrayFace( Level, i, j, stencilDirection::CELL_GRAY_C )] += tmp;
         }
      }
      --inner_rowsize;
   }

   inner_rowsize = rowsize;

   for ( size_t j = 0; j < rowsize - 2; ++j )
   {
      for ( size_t i = 0; i < inner_rowsize - 2; ++i )
      {
         // evalate velocities
         tmp = 1.0 / 3.0 *
               ( src[indexFromVertex( Level, i, j + 1, stencilDirection::VERTEX_C )] +
                 src[indexFromVertex( Level, i + 1, j + 1, stencilDirection::VERTEX_C )] +
                 src[indexFromVertex( Level, i + 1, j, stencilDirection::VERTEX_C )] );

         if ( updateType == Replace )
         {
            dst[facedof::macroface::indexFaceFromBlueFace( Level, i, j, stencilDirection::CELL_BLUE_C )] = tmp;
         }
         else if ( updateType == Add )
         {
            dst[facedof::macroface::indexFaceFromBlueFace( Level, i, j, stencilDirection::CELL_BLUE_C )] += tmp;
         }
      }
      --inner_rowsize;
   }
}

} // namespace macroface
} // namespace dgfunction
} //namespace hyteg

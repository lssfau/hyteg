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

#include "hyteg/facedofspace_old/FaceDoFIndexing.hpp"

namespace hyteg {
namespace dgfunction {
namespace macroface {

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

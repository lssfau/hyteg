/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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
#include "AbstractApply.hpp"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"

/// these are copies of the abstract apply function found in HyTeG
/// adjusted a little for better comparison

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

void apply_2d_vertex_to_vertex( double* WALBERLA_RESTRICT dstPtr,
                                double const* WALBERLA_RESTRICT const srcPtr,
                                double const* WALBERLA_RESTRICT const stencilPtr,
                                const uint_t                          level )
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( level );
   for ( idx_t j = 1; j < idx_t( rowsize ) - 2; ++j )
   {
      for ( idx_t i = 1; i < idx_t( rowsize ) -j - 1; ++i )
      {
         auto tmp = real_t( 0 );
         for ( const auto direction : vertexdof::macroface::neighborsWithCenter )
         {
            tmp += stencilPtr[vertexdof::stencilIndexFromVertex( direction )] *
                   srcPtr[vertexdof::macroface::indexFromVertex( level, i, j, direction )];
         }
         dstPtr[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] = tmp;
      }
   }
}
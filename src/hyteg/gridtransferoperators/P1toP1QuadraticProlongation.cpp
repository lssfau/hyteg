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

#include "hyteg/gridtransferoperators/P1toP1QuadraticProlongation.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"

namespace hyteg {

void P1toP1QuadraticProlongation::prolongateMacroVertex( const real_t *src, real_t *dst, const uint_t & ) const
{
  dst[0] = src[0];
}

void P1toP1QuadraticProlongation::prolongateMacroEdge( const real_t *src, real_t *dst, const uint_t & sourceLevel ) const
{
  //TODO: rewrite using index function possible? maybe more generalized notion of Operator between different levels
  // is required.

  size_t rowsize_coarse = levelinfo::num_microvertices_per_edge(sourceLevel);
  size_t i_fine = 1;
  real_t invtemp = 1/8.;
  const real_t s1[3] = {3*invtemp, 6*invtemp, -invtemp};
  const real_t s2[3] = {-invtemp, 6*invtemp, 3*invtemp};

  size_t i_coarse;
  for (i_coarse = 0; i_coarse < rowsize_coarse - 2; ++i_coarse) {
    dst[i_fine] =
    (s1[0]*src[i_coarse] + s1[1]*src[i_coarse + 1] + s1[2]*src[i_coarse + 2]);
    dst[i_fine + 1] = src[i_coarse + 1];
    i_fine += 2;
  }
  i_coarse--;
  dst[i_fine] =
  (s2[0]*src[i_coarse] + s2[1]*src[i_coarse + 1] + s2[2]*src[i_coarse + 2]);
}

void P1toP1QuadraticProlongation::prolongateMacroFace( const real_t *src, real_t *dst, const uint_t & sourceLevel ) const
{
  typedef stencilDirection SD;
  using namespace vertexdof::macroface;

  uint_t N_c = levelinfo::num_microvertices_per_edge(sourceLevel);
  uint_t N_c_i = N_c;

  idx_t i, j;
  real_t linearx, lineary, linearxy, offx, offy, offxy;
  i = 0;
  for ( j = 2; j <= idx_t( N_c ) - 1; j += 2 ) {
// upper triangle inner points
//calculate offsets
    linearx = 0.5 * ( src[indexFromVertex( sourceLevel, i, j - 2, SD::VERTEX_C )] + src[indexFromVertex( sourceLevel, i, j, SD::VERTEX_C )] );
    lineary = 0.5 * ( src[indexFromVertex( sourceLevel, i + 2, j - 2, SD::VERTEX_C )] + src[indexFromVertex( sourceLevel,
                                                                                                       i, j - 2, SD::VERTEX_C )] );
    linearxy = 0.5 * ( src[indexFromVertex( sourceLevel, i + 2, j - 2, SD::VERTEX_C )] + src[indexFromVertex( sourceLevel,
                                                                                                        i, j, SD::VERTEX_C )] );

    offx = src[indexFromVertex( sourceLevel, i, j - 1, SD::VERTEX_C )] - linearx;
    offy = src[indexFromVertex( sourceLevel, i + 1, j - 2, SD::VERTEX_C )] - lineary;
    offxy = src[indexFromVertex( sourceLevel, i + 1, j - 1, SD::VERTEX_C )] - linearxy;

// left bottom corner
    dst[indexFromVertex( sourceLevel + 1, 2 * i + 1, 2 * j - 3, SD::VERTEX_C )] = 0.5 * ( linearx + lineary ) + 0.5 * offx + 0.5 * offy + 0.25 * offxy;
// right bottom corner
    dst[indexFromVertex( sourceLevel + 1, 2 * i + 1, 2 * j - 2, SD::VERTEX_C )] = 0.5 * ( linearx + linearxy ) + 0.5 * offx + 0.25 * offy + 0.5 * offxy;
// top corner
    dst[indexFromVertex( sourceLevel + 1, 2 * i + 2, 2 * j - 3, SD::VERTEX_C )] = 0.5 * ( linearxy + lineary ) + 0.25 * offx + 0.5 * offy + 0.5 * offxy;
  }

  N_c_i -= 1;

  for ( j = 2; j < idx_t( N_c ) - 1; j += 2 ) {
    for ( i = 2; i < idx_t( N_c_i ) - 1; i += 2 ) {

// upper triangle inner points
//calculate offsets
      linearx = 0.5 * ( src[indexFromVertex( sourceLevel, i, j - 2, SD::VERTEX_C )] + src[indexFromVertex( sourceLevel,
                                                                                                     i, j, SD::VERTEX_C )] );
      lineary = 0.5 * ( src[indexFromVertex( sourceLevel, i + 2, j - 2, SD::VERTEX_C )] + src[indexFromVertex( sourceLevel,
                                                                                                         i, j - 2, SD::VERTEX_C )] );
      linearxy = 0.5 * ( src[indexFromVertex( sourceLevel, i + 2, j - 2, SD::VERTEX_C )] + src[indexFromVertex( sourceLevel,
                                                                                                          i, j, SD::VERTEX_C )] );

      offx = src[indexFromVertex( sourceLevel, i, j - 1, SD::VERTEX_C )] - linearx;
      offy = src[indexFromVertex( sourceLevel, i + 1, j - 2, SD::VERTEX_C )] - lineary;
      offxy = src[indexFromVertex( sourceLevel, i + 1, j - 1, SD::VERTEX_C )] - linearxy;
// left bottom corner
      dst[indexFromVertex( sourceLevel + 1, 2 * i + 1, 2 * j - 3, SD::VERTEX_C )] = 0.5 * ( linearx + lineary ) + 0.5 * offx + 0.5 * offy + 0.25 * offxy;
// right bottom corner
      dst[indexFromVertex( sourceLevel + 1, 2 * i + 1, 2 * j - 2, SD::VERTEX_C )] = 0.5 * ( linearx + linearxy ) + 0.5 * offx + 0.25 * offy + 0.5 * offxy;
// top corner
      dst[indexFromVertex( sourceLevel + 1, 2 * i + 2, 2 * j - 3, SD::VERTEX_C )] = 0.5 * ( linearxy + lineary ) + 0.25 * offx + 0.5 * offy + 0.5 * offxy;

// lower triangle all points
//calculate offsets
      lineary = 0.5 * ( src[indexFromVertex( sourceLevel, i - 2, j, SD::VERTEX_C )] + src[indexFromVertex( sourceLevel,
                                                                                                     i, j, SD::VERTEX_C )] );
      linearxy = 0.5 * ( src[indexFromVertex( sourceLevel, i - 2, j, SD::VERTEX_C )] + src[indexFromVertex( sourceLevel,
                                                                                                      i, j - 2, SD::VERTEX_C )] );

      offy = src[indexFromVertex( sourceLevel, i - 1, j, SD::VERTEX_C )] - lineary;
      offxy = src[indexFromVertex( sourceLevel, i - 1, j - 1, SD::VERTEX_C )] - linearxy;
// first inner points
// left bottom corner
      dst[indexFromVertex( sourceLevel + 1, 2 * i - 1, 2 * j - 1, SD::VERTEX_C )] = 0.5 * ( linearx + lineary ) + 0.5 * offx + 0.5 * offy + 0.25 * offxy;
// right bottom corner
      dst[indexFromVertex( sourceLevel + 1, 2 * i - 1, 2 * j - 2, SD::VERTEX_C )] = 0.5 * ( linearx + linearxy ) + 0.5 * offx + 0.25 * offy + 0.5 * offxy;
// top corner
      dst[indexFromVertex( sourceLevel + 1, 2 * i - 2, 2 * j - 1, SD::VERTEX_C )] = 0.5 * ( linearxy + lineary ) + 0.25 * offx + 0.5 * offy + 0.5 * offxy;

// boundary points
// x-direction
      dst[indexFromVertex( sourceLevel + 1, 2 * i, 2 * j - 1, SD::VERTEX_C )] = 0.5 * ( linearx + src[indexFromVertex( sourceLevel,
                                                                                                                 i, j, SD::VERTEX_C )] ) + 0.75 * offx;
      dst[indexFromVertex( sourceLevel + 1, 2 * i, 2 * j - 3, SD::VERTEX_C )] = 0.5 * ( linearx + src[indexFromVertex( sourceLevel,
                                                                                                                 i, j - 2, SD::VERTEX_C )] ) + 0.75 * offx;
//y-direction
      dst[indexFromVertex( sourceLevel + 1, 2 * i - 1, 2 * j, SD::VERTEX_C )] = 0.5 * ( src[indexFromVertex( sourceLevel,

                                                                                                       i, j,
                                                                                                       SD::VERTEX_C )] + lineary ) + 0.75 * offy;
      dst[indexFromVertex( sourceLevel + 1, 2 * i - 3, 2 * j, SD::VERTEX_C )] = 0.5 * ( src[indexFromVertex( sourceLevel,

                                                                                                       i - 2, j,
                                                                                                       SD::VERTEX_C )] + lineary ) + 0.75 * offy;
//xy-direction
      dst[indexFromVertex( sourceLevel + 1, 2 * i - 1, 2 * j - 3, SD::VERTEX_C )] = 0.5 * ( src[indexFromVertex( sourceLevel,
                                                                                                           i, j - 2, SD::VERTEX_C )] + linearxy ) + 0.75 * offxy;
      dst[indexFromVertex( sourceLevel + 1, 2 * i - 3, 2 * j - 1, SD::VERTEX_C )] = 0.5 * ( src[indexFromVertex( sourceLevel,
                                                                                                           i - 2, j, SD::VERTEX_C )] + linearxy ) + 0.75 * offxy;
// coarse points
      dst[indexFromVertex( sourceLevel + 1, 2 * i, 2 * j, SD::VERTEX_C )] = src[indexFromVertex( sourceLevel,
                                                                                           i, j,
                                                                                           SD::VERTEX_C )];
      dst[indexFromVertex( sourceLevel + 1, 2 * i, 2 * j - 2, SD::VERTEX_C )] = src[indexFromVertex( sourceLevel,
                                                                                               i,
                                                                                               j - 1, SD::VERTEX_C )];
      dst[indexFromVertex( sourceLevel + 1, 2 * i - 2, 2 * j, SD::VERTEX_C )] = src[indexFromVertex( sourceLevel,

                                                                                               i - 1, j,
                                                                                               SD::VERTEX_C )];
      dst[indexFromVertex( sourceLevel + 1, 2 * i - 2, 2 * j - 2, SD::VERTEX_C )] = src[indexFromVertex( sourceLevel,

                                                                                                   i - 1, j - 1,
                                                                                                   SD::VERTEX_C )];
    }
    N_c_i -= 2;

  }
}

}
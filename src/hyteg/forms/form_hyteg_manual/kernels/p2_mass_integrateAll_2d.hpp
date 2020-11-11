/*
 * Copyright (c) 2020 Nils Kohl, Marcus Mohr.
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

// Internal Documentation: This is variant V3_star

elMat.setAll( real_c( 0 ) );

const auto coords_0_0 = coords[0][0];
const auto coords_0_1 = coords[0][1];

const auto coords_1_0 = coords[1][0];
const auto coords_1_1 = coords[1][1];

const auto coords_2_0 = coords[2][0];
const auto coords_2_1 = coords[2][1];

// compute Jacobian determinant of inverse pull-back mapping
real_t tmp0 = coords_0_0 - coords_1_0;
real_t tmp1 = coords_0_1 - coords_2_1;
real_t tmp2 = coords_0_0 - coords_2_0;
real_t tmp3 = coords_0_1 - coords_1_1;

real_t aux1 = tmp0 * tmp1 - tmp2 * tmp3;
real_t aux2 = aux1 * aux1;

real_t detDPhiInv = std::abs( aux1 );

// outermost loop is over the quadrature points
for ( uint_t k = 0; k < QUADWEIGHTS.size(); k++ )
{
   // determine barycentric coordinates for current integration point
   real_t L2 = QUADPOINTS[k][0];
   real_t L3 = QUADPOINTS[k][1];
   real_t L1 = 1.0 - L2 - L3;

   // map point to computational element (affine map Phi^{-1} )
   Point3D mappedPt;
   mappedPt[0] = ( coords_1_0 - coords_0_0 ) * L2 + ( coords_2_0 - coords_0_0 ) * L3 + coords_0_0;
   mappedPt[1] = ( coords_1_1 - coords_0_1 ) * L2 + ( coords_2_1 - coords_0_1 ) * L3 + coords_0_1;

   // compute Jacobian determiant of blending map Psi (computational to physical element)
   real_t detDPsi = std::abs( geometryMap_->evalDetDF( mappedPt ) );

   // Evaluate shape functions at current integration point
   std::array< real_t, 6 > sf;
   sf[0] = L1 * ( real_c(2) * L1 - real_c(1) );
   sf[1] = L2 * ( real_c(2) * L2 - real_c(1) );
   sf[2] = L3 * ( real_c(2) * L3 - real_c(1) );
   sf[3] = real_c(4) * L2 * L3;
   sf[4] = real_c(4) * L1 * L3;
   sf[5] = real_c(4) * L1 * L2;

   // compute and add contribution from current integration point to all element matrix entries
   // (computations on upper triangle only)
   for ( uint i = 0; i < 6; i++ )
   {
      for ( uint j = i; j < 6; j++ )
      {
         elMat( i, j ) += QUADWEIGHTS[k] * detDPhiInv * detDPsi * sf[i] * sf[j];
      }
   }
}

// set lower triangular part from symmetry
for ( uint i = 0; i < 6; i++ )
{
   for ( uint j = 0; j < i; j++ )
   {
      elMat( i, j ) = elMat( j, i );
   }
}

/*
 * Copyright (c) 2020 Marcus Mohr.
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

elMat.setAll( real_c( 0 ) );

const auto coords_0_0 = coords[0][0];
const auto coords_0_1 = coords[0][1];
const auto coords_0_2 = coords[0][2];

const auto coords_1_0 = coords[1][0];
const auto coords_1_1 = coords[1][1];
const auto coords_1_2 = coords[1][2];

const auto coords_2_0 = coords[2][0];
const auto coords_2_1 = coords[2][1];
const auto coords_2_2 = coords[2][2];

const auto coords_3_0 = coords[3][0];
const auto coords_3_1 = coords[3][1];
const auto coords_3_2 = coords[3][2];

// compute Jacobian determinant of inverse pull-back mapping
// and location independent auxilliary values
real_t tmp0  = -coords_0_0;
real_t tmp1  = tmp0 + coords_1_0;
real_t tmp2  = -coords_0_1;
real_t tmp3  = tmp2 + coords_2_1;
real_t tmp4  = -coords_0_2;
real_t tmp5  = tmp4 + coords_3_2;
real_t tmp6  = tmp0 + coords_2_0;
real_t tmp7  = tmp2 + coords_3_1;
real_t tmp8  = tmp4 + coords_1_2;
real_t tmp9  = tmp0 + coords_3_0;
real_t tmp10 = tmp2 + coords_1_1;
real_t tmp11 = tmp4 + coords_2_2;

real_t detDPhiInv = -tmp1 * tmp11 * tmp7 + tmp1 * tmp3 * tmp5 + tmp10 * tmp11 * tmp9 - tmp10 * tmp5 * tmp6 - tmp3 * tmp8 * tmp9 +
                    tmp6 * tmp7 * tmp8;
detDPhiInv = std::abs( detDPhiInv );

// outermost loop is over the cubature points
for ( uint_t k = 0; k < CUBAWEIGHTS.size(); k++ )
{
   // determine barycentric coordinates for current integration point
   real_t L2 = CUBAPOINTS[k][0];
   real_t L3 = CUBAPOINTS[k][1];
   real_t L4 = CUBAPOINTS[k][2];
   real_t L1 = 1.0 - L2 - L3 - L4;

   // map point to computational element (affine map Phi^{-1} )
   Point3D mappedPt;
   mappedPt[0] =
       ( coords_1_0 - coords_0_0 ) * L2 + ( coords_2_0 - coords_0_0 ) * L3 + ( coords_3_0 - coords_0_0 ) * L4 + coords_0_0;
   mappedPt[1] =
       ( coords_1_1 - coords_0_1 ) * L2 + ( coords_2_1 - coords_0_1 ) * L3 + ( coords_3_1 - coords_0_1 ) * L4 + coords_0_1;
   mappedPt[2] =
       ( coords_1_2 - coords_0_2 ) * L2 + ( coords_2_2 - coords_0_2 ) * L3 + ( coords_3_2 - coords_0_2 ) * L4 + coords_0_2;

   // compute Jacobian determiant of blending map Psi (computational to physical element)
   Matrix3r dummy;
   real_t detDPsi = std::abs( geometryMap_->evalDF( mappedPt, dummy ) );

   // Evaluate shape functions at current integration point
   std::array< real_t, 10 > sf;

   sf[0] = L1 * ( real_c(2) * L1 - real_c(1) );
   sf[1] = L2 * ( real_c(2) * L2 - real_c(1) );
   sf[2] = L3 * ( real_c(2) * L3 - real_c(1) );
   sf[3] = L4 * ( real_c(2) * L4 - real_c(1) );
   sf[4] = real_c(4) * L3 * L4;
   sf[5] = real_c(4) * L2 * L4;
   sf[6] = real_c(4) * L2 * L3;
   sf[7] = real_c(4) * L1 * L4;
   sf[8] = real_c(4) * L1 * L3;
   sf[9] = real_c(4) * L1 * L2;

   // compute and add contribution from current integration point to all element matrix entries
   // (computations on upper triangle only)
   for ( uint i = 0; i < 10; i++ )
   {
      for ( uint j = i; j < 10; j++ )
      {
        elMat( i, j ) += CUBAWEIGHTS[k] * detDPhiInv * detDPsi * sf[i] * sf[j];        
      }
   }
}

// set lower triangular part from symmetry
for ( uint i = 0; i < 10; i++ )
{
   for ( uint j = 0; j < i; j++ )
   {
      elMat( i, j ) = elMat( j, i );
   }
}

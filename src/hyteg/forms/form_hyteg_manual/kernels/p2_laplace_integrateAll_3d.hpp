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

elMat.setZero();

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
real_t detDPhiFac = std::abs( detDPhiInv ) / ( detDPhiInv * detDPhiInv );

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

   // compute derivative of non-affine blending map Psi (computational to physical element)
   Matrix3r DPsi;
   geometryMap_->evalDF( mappedPt, DPsi );

   const real_t DPsi_local_0_0 = DPsi( 0, 0 );
   const real_t DPsi_local_0_1 = DPsi( 0, 1 );
   const real_t DPsi_local_0_2 = DPsi( 0, 2 );

   const real_t DPsi_local_1_0 = DPsi( 1, 0 );
   const real_t DPsi_local_1_1 = DPsi( 1, 1 );
   const real_t DPsi_local_1_2 = DPsi( 1, 2 );

   const real_t DPsi_local_2_0 = DPsi( 2, 0 );
   const real_t DPsi_local_2_1 = DPsi( 2, 1 );
   const real_t DPsi_local_2_2 = DPsi( 2, 2 );

   real_t tmp12 = DPsi_local_0_1 * DPsi_local_1_0;
   real_t tmp13 = -DPsi_local_1_1 * DPsi_local_0_0 + tmp12;

   real_t tmp24 = DPsi_local_1_2 * DPsi_local_0_0;
   real_t tmp25 = DPsi_local_0_2 * DPsi_local_1_0 - tmp24;

   real_t tmp27 = DPsi_local_0_1 * DPsi_local_1_2;
   real_t tmp28 = DPsi_local_0_2 * DPsi_local_1_1;
   real_t tmp29 = tmp27 - tmp28;

   real_t detDPsi = DPsi_local_0_2 * DPsi_local_2_1 * DPsi_local_1_0 + DPsi_local_1_1 * DPsi_local_2_2 * DPsi_local_0_0 -
                    DPsi_local_2_1 * tmp24 - DPsi_local_2_2 * tmp12 + DPsi_local_2_0 * tmp27 - DPsi_local_2_0 * tmp28;
   real_t detDPsiFac = std::abs( detDPsi ) / ( detDPsi * detDPsi );

   // Evaluate shape function derivatives at current integration point
   std::array< std::array< real_t, 3 >, 10 > sfd;

   sfd[0][0] = L2 + L3 + L4 - 3.0 * L1;
   sfd[0][1] = L2 + L3 + L4 - 3.0 * L1;
   sfd[0][2] = L2 + L3 + L4 - 3.0 * L1;

   sfd[1][0] = 4.0 * L2 - 1.0;
   sfd[1][1] = 0.0;
   sfd[1][2] = 0.0;

   sfd[2][0] = 0.0;
   sfd[2][1] = 4.0 * L3 - 1.0;
   sfd[2][2] = 0.0;

   sfd[3][0] = 0.0;
   sfd[3][1] = 0.0;
   sfd[3][2] = 4.0 * L4 - 1.0;

   sfd[4][0] = 0.0;
   sfd[4][1] = 4.0 * L4;
   sfd[4][2] = 4.0 * L3;

   sfd[5][0] = 4.0 * L4;
   sfd[5][1] = 0.0;
   sfd[5][2] = 4.0 * L2;

   sfd[6][0] = 4.0 * L3;
   sfd[6][1] = 4.0 * L2;
   sfd[6][2] = 0.0;

   sfd[7][0] = -4.0 * L4;
   sfd[7][1] = -4.0 * L4;
   sfd[7][2] = 4.0 * ( L1 - L4 );

   sfd[8][0] = -4.0 * L3;
   sfd[8][1] = 4.0 * ( L1 - L3 );
   sfd[8][2] = -4.0 * L3;

   sfd[9][0] = 4.0 * ( L1 - L2 );
   sfd[9][1] = -4.0 * L2;
   sfd[9][2] = -4.0 * L2;

   // compute and add contribution from current integration point to all element matrix entries
   // (computations on upper triangle only)
   for ( uint_t i = 0; i < 10; i++ )
   {
      real_t tmp14 = sfd[i][0] * coords_0_0;
      real_t tmp15 = sfd[i][0] * coords_2_0;
      real_t tmp16 = sfd[i][0] * coords_3_0;
      real_t tmp17 = sfd[i][1] * coords_0_0;
      real_t tmp18 = sfd[i][1] * coords_1_0;
      real_t tmp19 = sfd[i][1] * coords_3_0;
      real_t tmp20 = sfd[i][2] * coords_0_0;
      real_t tmp21 = sfd[i][2] * coords_1_0;
      real_t tmp22 = sfd[i][2] * coords_2_0;
      real_t tmp23 = tmp14 * coords_2_1 - tmp14 * coords_3_1 - tmp15 * coords_0_1 + tmp15 * coords_3_1 + tmp16 * coords_0_1 -
                     tmp16 * coords_2_1 - tmp17 * coords_1_1 + tmp17 * coords_3_1 + tmp18 * coords_0_1 - tmp18 * coords_3_1 -
                     tmp19 * coords_0_1 + tmp19 * coords_1_1 + tmp20 * coords_1_1 - tmp20 * coords_2_1 - tmp21 * coords_0_1 +
                     tmp21 * coords_2_1 + tmp22 * coords_0_1 - tmp22 * coords_1_1;

      real_t tmp26 = tmp14 * coords_2_2 - tmp14 * coords_3_2 - tmp15 * coords_0_2 + tmp15 * coords_3_2 + tmp16 * coords_0_2 -
                     tmp16 * coords_2_2 - tmp17 * coords_1_2 + tmp17 * coords_3_2 + tmp18 * coords_0_2 - tmp18 * coords_3_2 -
                     tmp19 * coords_0_2 + tmp19 * coords_1_2 + tmp20 * coords_1_2 - tmp20 * coords_2_2 - tmp21 * coords_0_2 +
                     tmp21 * coords_2_2 + tmp22 * coords_0_2 - tmp22 * coords_1_2;

      real_t tmp30 = sfd[i][0] * coords_0_1;
      real_t tmp31 = sfd[i][0] * coords_2_1;
      real_t tmp32 = sfd[i][0] * coords_3_1;
      real_t tmp33 = sfd[i][1] * coords_0_1;
      real_t tmp34 = sfd[i][1] * coords_1_1;
      real_t tmp35 = sfd[i][1] * coords_3_1;
      real_t tmp36 = sfd[i][2] * coords_0_1;
      real_t tmp37 = sfd[i][2] * coords_1_1;
      real_t tmp38 = sfd[i][2] * coords_2_1;

      real_t tmp39 = tmp30 * coords_2_2 - tmp30 * coords_3_2 - tmp31 * coords_0_2 + tmp31 * coords_3_2 + tmp32 * coords_0_2 -
                     tmp32 * coords_2_2 - tmp33 * coords_1_2 + tmp33 * coords_3_2 + tmp34 * coords_0_2 - tmp34 * coords_3_2 -
                     tmp35 * coords_0_2 + tmp35 * coords_1_2 + tmp36 * coords_1_2 - tmp36 * coords_2_2 - tmp37 * coords_0_2 +
                     tmp37 * coords_2_2 + tmp38 * coords_0_2 - tmp38 * coords_1_2;

      for ( uint_t j = i; j < 10; j++ )
      {
         real_t tmp40 = sfd[j][0] * coords_0_0;
         real_t tmp41 = sfd[j][0] * coords_2_0;
         real_t tmp42 = sfd[j][0] * coords_3_0;
         real_t tmp43 = sfd[j][1] * coords_0_0;
         real_t tmp44 = sfd[j][1] * coords_1_0;
         real_t tmp45 = sfd[j][1] * coords_3_0;
         real_t tmp46 = sfd[j][2] * coords_0_0;
         real_t tmp47 = sfd[j][2] * coords_1_0;
         real_t tmp48 = sfd[j][2] * coords_2_0;

         real_t tmp49 = tmp40 * coords_2_1 - tmp40 * coords_3_1 - tmp41 * coords_0_1 + tmp41 * coords_3_1 + tmp42 * coords_0_1 -
                        tmp42 * coords_2_1 - tmp43 * coords_1_1 + tmp43 * coords_3_1 + tmp44 * coords_0_1 - tmp44 * coords_3_1 -
                        tmp45 * coords_0_1 + tmp45 * coords_1_1 + tmp46 * coords_1_1 - tmp46 * coords_2_1 - tmp47 * coords_0_1 +
                        tmp47 * coords_2_1 + tmp48 * coords_0_1 - tmp48 * coords_1_1;

         real_t tmp50 = tmp40 * coords_2_2 - tmp40 * coords_3_2 - tmp41 * coords_0_2 + tmp41 * coords_3_2 + tmp42 * coords_0_2 -
                        tmp42 * coords_2_2 - tmp43 * coords_1_2 + tmp43 * coords_3_2 + tmp44 * coords_0_2 - tmp44 * coords_3_2 -
                        tmp45 * coords_0_2 + tmp45 * coords_1_2 + tmp46 * coords_1_2 - tmp46 * coords_2_2 - tmp47 * coords_0_2 +
                        tmp47 * coords_2_2 + tmp48 * coords_0_2 - tmp48 * coords_1_2;

         real_t tmp51 = sfd[j][0] * coords_0_1;
         real_t tmp52 = sfd[j][0] * coords_2_1;
         real_t tmp53 = sfd[j][0] * coords_3_1;
         real_t tmp54 = sfd[j][1] * coords_0_1;
         real_t tmp55 = sfd[j][1] * coords_1_1;
         real_t tmp56 = sfd[j][1] * coords_3_1;
         real_t tmp57 = sfd[j][2] * coords_0_1;
         real_t tmp58 = sfd[j][2] * coords_1_1;
         real_t tmp59 = sfd[j][2] * coords_2_1;

         real_t tmp60 = tmp51 * coords_2_2 - tmp51 * coords_3_2 - tmp52 * coords_0_2 + tmp52 * coords_3_2 + tmp53 * coords_0_2 -
                        tmp53 * coords_2_2 - tmp54 * coords_1_2 + tmp54 * coords_3_2 + tmp55 * coords_0_2 - tmp55 * coords_3_2 -
                        tmp56 * coords_0_2 + tmp56 * coords_1_2 + tmp57 * coords_1_2 - tmp57 * coords_2_2 - tmp58 * coords_0_2 +
                        tmp58 * coords_2_2 + tmp59 * coords_0_2 - tmp59 * coords_1_2;

         real_t tmp61 = DPsi_local_0_1 * DPsi_local_2_0 - DPsi_local_2_1 * DPsi_local_0_0;
         real_t tmp62 = DPsi_local_0_2 * DPsi_local_2_0 - DPsi_local_2_2 * DPsi_local_0_0;
         real_t tmp63 = DPsi_local_0_1 * DPsi_local_2_2 - DPsi_local_0_2 * DPsi_local_2_1;
         real_t tmp64 = DPsi_local_1_1 * DPsi_local_2_0 - DPsi_local_2_1 * DPsi_local_1_0;
         real_t tmp65 = DPsi_local_1_2 * DPsi_local_2_0 - DPsi_local_2_2 * DPsi_local_1_0;
         real_t tmp66 = DPsi_local_1_1 * DPsi_local_2_2 - DPsi_local_1_2 * DPsi_local_2_1;

         real_t aux = ( tmp13 * tmp23 + tmp25 * tmp26 - tmp29 * tmp39 ) * ( tmp13 * tmp49 + tmp25 * tmp50 - tmp29 * tmp60 ) +
                      ( tmp23 * tmp61 + tmp26 * tmp62 - tmp39 * tmp63 ) * ( tmp49 * tmp61 + tmp50 * tmp62 - tmp60 * tmp63 ) +
                      ( tmp23 * tmp64 + tmp26 * tmp65 - tmp39 * tmp66 ) * ( tmp49 * tmp64 + tmp50 * tmp65 - tmp60 * tmp66 );

         elMat( i, j ) += CUBAWEIGHTS[k] * detDPhiFac * detDPsiFac * aux;
      }
   }
}

// set lower triangular part from symmetry
for ( uint_t i = 0; i < 10; i++ )
{
   for ( uint_t j = 0; j < i; j++ )
   {
      elMat( i, j ) = elMat( j, i );
   }
}

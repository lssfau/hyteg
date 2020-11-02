/*
 * Copyright (c) 2017-2020 Marcus Mohr.
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

#include "hyteg/forms/form_hyteg_base/P2FormHyTeG.hpp"
#include "hyteg/forms/form_hyteg_manual/QuadratureRules.hpp"
#include "hyteg/geometry/GeometryMap.hpp"

// NOTE:
//
// This form is for performance evaluation of different variants.
// Currently it implemtens the following four variants
//
// V_1: Treat shape functions N_i and N_j in innermost loop + computed shape function derivatives always
// V_2: Treat shape functions N_i and N_j in innermost loop + precompute shape function derivatives
// V_3: Treat shape function N_i in loop over i and N_j in loop over j + computed shape function derivatives always
// V_4: Treat shape function N_i in loop over i and N_j in loop over j + precompute shape function derivatives
//
// To run each variant set the following macros
// V_1: none
// V_2: PRECOMUTE_SHAPE_FUNCTION_DERIVATIVES
// V_3: SLOWER
// V_4: SLOWER + PRECOMUTE_SHAPE_FUNCTION_DERIVATIVES
//
#define SLOWER
// #define PRECOMUTE_SHAPE_FUNCTION_DERIVATIVES

#ifdef PRECOMUTE_SHAPE_FUNCTION_DERIVATIVES
#define SFD shapeFunctionDerivatives_[k]
#else
#define SFD sfd
#endif

namespace hyteg {

class P2Form_laplacePimped3D : public P2FormHyTeG
{
#ifdef PRECOMUTE_SHAPE_FUNCTION_DERIVATIVES
 public:
   P2Form_laplacePimped3D()
   {
// Select cubature rule
#define CUBAPOINTS cubature::T3_points
#define CUBAWEIGHTS cubature::T3_weights

      for ( uint_t k = 0; k < CUBAWEIGHTS.size(); k++ )
      {
         // determine barycentric coordinates for current integration point
         real_t L2 = CUBAPOINTS[k][0];
         real_t L3 = CUBAPOINTS[k][1];
         real_t L4 = CUBAPOINTS[k][2];
         real_t L1 = 1.0 - L2 - L3 - L4;

         shapeFunctionDerivatives_[k][0][0] = L2 + L3 + L4 - 3.0 * L1;
         shapeFunctionDerivatives_[k][0][1] = L2 + L3 + L4 - 3.0 * L1;
         shapeFunctionDerivatives_[k][0][2] = L2 + L3 + L4 - 3.0 * L1;

         shapeFunctionDerivatives_[k][1][0] = 4.0 * L2 - 1.0;
         shapeFunctionDerivatives_[k][1][1] = 0.0;
         shapeFunctionDerivatives_[k][1][2] = 0.0;

         shapeFunctionDerivatives_[k][2][0] = 0.0;
         shapeFunctionDerivatives_[k][2][1] = 4.0 * L3 - 1.0;
         shapeFunctionDerivatives_[k][2][2] = 0.0;

         shapeFunctionDerivatives_[k][3][0] = 0.0;
         shapeFunctionDerivatives_[k][3][1] = 0.0;
         shapeFunctionDerivatives_[k][3][2] = 4.0 * L4 - 1.0;

         shapeFunctionDerivatives_[k][4][0] = 0.0;
         shapeFunctionDerivatives_[k][4][1] = 4.0 * L4;
         shapeFunctionDerivatives_[k][4][2] = 4.0 * L3;

         shapeFunctionDerivatives_[k][5][0] = 4.0 * L4;
         shapeFunctionDerivatives_[k][5][1] = 0.0;
         shapeFunctionDerivatives_[k][5][2] = 4.0 * L2;

         shapeFunctionDerivatives_[k][6][0] = 4.0 * L3;
         shapeFunctionDerivatives_[k][6][1] = 4.0 * L2;
         shapeFunctionDerivatives_[k][6][2] = 0.0;

         shapeFunctionDerivatives_[k][7][0] = -4.0 * L4;
         shapeFunctionDerivatives_[k][7][1] = -4.0 * L4;
         shapeFunctionDerivatives_[k][7][2] = 4.0 * ( L1 - L4 );

         shapeFunctionDerivatives_[k][8][0] = -4.0 * L3;
         shapeFunctionDerivatives_[k][8][1] = 4.0 * ( L1 - L3 );
         shapeFunctionDerivatives_[k][8][2] = -4.0 * L3;

         shapeFunctionDerivatives_[k][9][0] = 4.0 * ( L1 - L2 );
         shapeFunctionDerivatives_[k][9][1] = -4.0 * L2;
         shapeFunctionDerivatives_[k][9][2] = -4.0 * L2;
      }
   };

 private:
   std::array< std::array< std::array< real_t, 3 >, 10 >, CUBAWEIGHTS.size() > shapeFunctionDerivatives_;

#endif

 public:
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix6r& elMat ) const final
   {
      WALBERLA_ABORT( "P2Form_laplacePimped3D only works in 3D as the name suggests" );
   };

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix10r& elMat ) const final
   {
// Select cubature rule
#define CUBAPOINTS cubature::T3_points
#define CUBAWEIGHTS cubature::T3_weights

      // initialise element matrix to zero for additive assembly
      elMat.setAll( real_c( 0 ) );

      // compute Jacobian determinant of inverse pull-back mapping
      // and location independent auxilliary values
      real_t tmp0  = -coords[0][0];
      real_t tmp1  = tmp0 + coords[1][0];
      real_t tmp2  = -coords[0][1];
      real_t tmp3  = tmp2 + coords[2][1];
      real_t tmp4  = -coords[0][2];
      real_t tmp5  = tmp4 + coords[3][2];
      real_t tmp6  = tmp0 + coords[2][0];
      real_t tmp7  = tmp2 + coords[3][1];
      real_t tmp8  = tmp4 + coords[1][2];
      real_t tmp9  = tmp0 + coords[3][0];
      real_t tmp10 = tmp2 + coords[1][1];
      real_t tmp11 = tmp4 + coords[2][2];

      real_t detDPhiInv = -tmp1 * tmp11 * tmp7 + tmp1 * tmp3 * tmp5 + tmp10 * tmp11 * tmp9 - tmp10 * tmp5 * tmp6 -
                          tmp3 * tmp8 * tmp9 + tmp6 * tmp7 * tmp8;
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
         mappedPt[0] = ( coords[1][0] - coords[0][0] ) * L2 + ( coords[2][0] - coords[0][0] ) * L3 +
                       ( coords[3][0] - coords[0][0] ) * L4 + coords[0][0];
         mappedPt[1] = ( coords[1][1] - coords[0][1] ) * L2 + ( coords[2][1] - coords[0][1] ) * L3 +
                       ( coords[3][1] - coords[0][1] ) * L4 + coords[0][1];
         mappedPt[2] = ( coords[1][2] - coords[0][2] ) * L2 + ( coords[2][2] - coords[0][2] ) * L3 +
                       ( coords[3][2] - coords[0][2] ) * L4 + coords[0][2];

         // compute derivative of non-affine blending map Psi (computational to physical element)
         Matrix3r DPsi;
         geometryMap_->evalDF( mappedPt, DPsi );

         real_t tmp12 = DPsi( 0, 1 ) * DPsi( 1, 0 );
         real_t tmp13 = -DPsi( 1, 1 ) * DPsi( 0, 0 ) + tmp12;

         real_t tmp24 = DPsi( 1, 2 ) * DPsi( 0, 0 );
         real_t tmp25 = DPsi( 0, 2 ) * DPsi( 1, 0 ) - tmp24;

         real_t tmp27 = DPsi( 0, 1 ) * DPsi( 1, 2 );
         real_t tmp28 = DPsi( 0, 2 ) * DPsi( 1, 1 );
         real_t tmp29 = tmp27 - tmp28;

         real_t detDPsi = DPsi( 0, 2 ) * DPsi( 2, 1 ) * DPsi( 1, 0 ) + DPsi( 1, 1 ) * DPsi( 2, 2 ) * DPsi( 0, 0 ) -
                          DPsi( 2, 1 ) * tmp24 - DPsi( 2, 2 ) * tmp12 + DPsi( 2, 0 ) * tmp27 - DPsi( 2, 0 ) * tmp28;
         real_t detDPsiFac = std::abs( detDPsi ) / ( detDPsi * detDPsi );

#ifndef PRECOMUTE_SHAPE_FUNCTION_DERIVATIVES
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
#endif

#ifdef SLOWER

         // NOTE: This version performs fewer operations, but runs slower in the tests

         // compute and add contribution from current integration point to all element matrix entries
         // (computations on upper triangle only)
         for ( uint i = 0; i < 10; i++ )
         {
            real_t tmp14 = SFD[i][0] * coords[0][0];
            real_t tmp15 = SFD[i][0] * coords[2][0];
            real_t tmp16 = SFD[i][0] * coords[3][0];
            real_t tmp17 = SFD[i][1] * coords[0][0];
            real_t tmp18 = SFD[i][1] * coords[1][0];
            real_t tmp19 = SFD[i][1] * coords[3][0];
            real_t tmp20 = SFD[i][2] * coords[0][0];
            real_t tmp21 = SFD[i][2] * coords[1][0];
            real_t tmp22 = SFD[i][2] * coords[2][0];
            real_t tmp23 = tmp14 * coords[2][1] - tmp14 * coords[3][1] - tmp15 * coords[0][1] + tmp15 * coords[3][1] +
                           tmp16 * coords[0][1] - tmp16 * coords[2][1] - tmp17 * coords[1][1] + tmp17 * coords[3][1] +
                           tmp18 * coords[0][1] - tmp18 * coords[3][1] - tmp19 * coords[0][1] + tmp19 * coords[1][1] +
                           tmp20 * coords[1][1] - tmp20 * coords[2][1] - tmp21 * coords[0][1] + tmp21 * coords[2][1] +
                           tmp22 * coords[0][1] - tmp22 * coords[1][1];

            real_t tmp26 = tmp14 * coords[2][2] - tmp14 * coords[3][2] - tmp15 * coords[0][2] + tmp15 * coords[3][2] +
                           tmp16 * coords[0][2] - tmp16 * coords[2][2] - tmp17 * coords[1][2] + tmp17 * coords[3][2] +
                           tmp18 * coords[0][2] - tmp18 * coords[3][2] - tmp19 * coords[0][2] + tmp19 * coords[1][2] +
                           tmp20 * coords[1][2] - tmp20 * coords[2][2] - tmp21 * coords[0][2] + tmp21 * coords[2][2] +
                           tmp22 * coords[0][2] - tmp22 * coords[1][2];

            real_t tmp30 = SFD[i][0] * coords[0][1];
            real_t tmp31 = SFD[i][0] * coords[2][1];
            real_t tmp32 = SFD[i][0] * coords[3][1];
            real_t tmp33 = SFD[i][1] * coords[0][1];
            real_t tmp34 = SFD[i][1] * coords[1][1];
            real_t tmp35 = SFD[i][1] * coords[3][1];
            real_t tmp36 = SFD[i][2] * coords[0][1];
            real_t tmp37 = SFD[i][2] * coords[1][1];
            real_t tmp38 = SFD[i][2] * coords[2][1];

            real_t tmp39 = tmp30 * coords[2][2] - tmp30 * coords[3][2] - tmp31 * coords[0][2] + tmp31 * coords[3][2] +
                           tmp32 * coords[0][2] - tmp32 * coords[2][2] - tmp33 * coords[1][2] + tmp33 * coords[3][2] +
                           tmp34 * coords[0][2] - tmp34 * coords[3][2] - tmp35 * coords[0][2] + tmp35 * coords[1][2] +
                           tmp36 * coords[1][2] - tmp36 * coords[2][2] - tmp37 * coords[0][2] + tmp37 * coords[2][2] +
                           tmp38 * coords[0][2] - tmp38 * coords[1][2];

            for ( uint j = i; j < 10; j++ )
            {
               real_t tmp40 = SFD[j][0] * coords[0][0];
               real_t tmp41 = SFD[j][0] * coords[2][0];
               real_t tmp42 = SFD[j][0] * coords[3][0];
               real_t tmp43 = SFD[j][1] * coords[0][0];
               real_t tmp44 = SFD[j][1] * coords[1][0];
               real_t tmp45 = SFD[j][1] * coords[3][0];
               real_t tmp46 = SFD[j][2] * coords[0][0];
               real_t tmp47 = SFD[j][2] * coords[1][0];
               real_t tmp48 = SFD[j][2] * coords[2][0];

               real_t tmp49 = tmp40 * coords[2][1] - tmp40 * coords[3][1] - tmp41 * coords[0][1] + tmp41 * coords[3][1] +
                              tmp42 * coords[0][1] - tmp42 * coords[2][1] - tmp43 * coords[1][1] + tmp43 * coords[3][1] +
                              tmp44 * coords[0][1] - tmp44 * coords[3][1] - tmp45 * coords[0][1] + tmp45 * coords[1][1] +
                              tmp46 * coords[1][1] - tmp46 * coords[2][1] - tmp47 * coords[0][1] + tmp47 * coords[2][1] +
                              tmp48 * coords[0][1] - tmp48 * coords[1][1];

               real_t tmp50 = tmp40 * coords[2][2] - tmp40 * coords[3][2] - tmp41 * coords[0][2] + tmp41 * coords[3][2] +
                              tmp42 * coords[0][2] - tmp42 * coords[2][2] - tmp43 * coords[1][2] + tmp43 * coords[3][2] +
                              tmp44 * coords[0][2] - tmp44 * coords[3][2] - tmp45 * coords[0][2] + tmp45 * coords[1][2] +
                              tmp46 * coords[1][2] - tmp46 * coords[2][2] - tmp47 * coords[0][2] + tmp47 * coords[2][2] +
                              tmp48 * coords[0][2] - tmp48 * coords[1][2];

               real_t tmp51 = SFD[j][0] * coords[0][1];
               real_t tmp52 = SFD[j][0] * coords[2][1];
               real_t tmp53 = SFD[j][0] * coords[3][1];
               real_t tmp54 = SFD[j][1] * coords[0][1];
               real_t tmp55 = SFD[j][1] * coords[1][1];
               real_t tmp56 = SFD[j][1] * coords[3][1];
               real_t tmp57 = SFD[j][2] * coords[0][1];
               real_t tmp58 = SFD[j][2] * coords[1][1];
               real_t tmp59 = SFD[j][2] * coords[2][1];

               real_t tmp60 = tmp51 * coords[2][2] - tmp51 * coords[3][2] - tmp52 * coords[0][2] + tmp52 * coords[3][2] +
                              tmp53 * coords[0][2] - tmp53 * coords[2][2] - tmp54 * coords[1][2] + tmp54 * coords[3][2] +
                              tmp55 * coords[0][2] - tmp55 * coords[3][2] - tmp56 * coords[0][2] + tmp56 * coords[1][2] +
                              tmp57 * coords[1][2] - tmp57 * coords[2][2] - tmp58 * coords[0][2] + tmp58 * coords[2][2] +
                              tmp59 * coords[0][2] - tmp59 * coords[1][2];

               real_t tmp61 = DPsi( 0, 1 ) * DPsi( 2, 0 ) - DPsi( 2, 1 ) * DPsi( 0, 0 );
               real_t tmp62 = DPsi( 0, 2 ) * DPsi( 2, 0 ) - DPsi( 2, 2 ) * DPsi( 0, 0 );
               real_t tmp63 = DPsi( 0, 1 ) * DPsi( 2, 2 ) - DPsi( 0, 2 ) * DPsi( 2, 1 );
               real_t tmp64 = DPsi( 1, 1 ) * DPsi( 2, 0 ) - DPsi( 2, 1 ) * DPsi( 1, 0 );
               real_t tmp65 = DPsi( 1, 2 ) * DPsi( 2, 0 ) - DPsi( 2, 2 ) * DPsi( 1, 0 );
               real_t tmp66 = DPsi( 1, 1 ) * DPsi( 2, 2 ) - DPsi( 1, 2 ) * DPsi( 2, 1 );

               real_t aux =
                   ( tmp13 * tmp23 + tmp25 * tmp26 - tmp29 * tmp39 ) * ( tmp13 * tmp49 + tmp25 * tmp50 - tmp29 * tmp60 ) +
                   ( tmp23 * tmp61 + tmp26 * tmp62 - tmp39 * tmp63 ) * ( tmp49 * tmp61 + tmp50 * tmp62 - tmp60 * tmp63 ) +
                   ( tmp23 * tmp64 + tmp26 * tmp65 - tmp39 * tmp66 ) * ( tmp49 * tmp64 + tmp50 * tmp65 - tmp60 * tmp66 );

               elMat( i, j ) += CUBAWEIGHTS[k] * detDPhiFac * detDPsiFac * aux;
            }
         }

#else
         // compute and add contribution from current integration point to all element matrix entries
         // (computations on upper triangle only)
         for ( uint i = 0; i < 10; i++ )
         {
            for ( uint j = i; j < 10; j++ )
            {
               real_t tmp14 = SFD[i][0] * coords[0][0];
               real_t tmp15 = SFD[i][0] * coords[2][0];
               real_t tmp16 = SFD[i][0] * coords[3][0];
               real_t tmp17 = SFD[i][1] * coords[0][0];
               real_t tmp18 = SFD[i][1] * coords[1][0];
               real_t tmp19 = SFD[i][1] * coords[3][0];
               real_t tmp20 = SFD[i][2] * coords[0][0];
               real_t tmp21 = SFD[i][2] * coords[1][0];
               real_t tmp22 = SFD[i][2] * coords[2][0];
               real_t tmp23 = tmp14 * coords[2][1] - tmp14 * coords[3][1] - tmp15 * coords[0][1] + tmp15 * coords[3][1] +
                              tmp16 * coords[0][1] - tmp16 * coords[2][1] - tmp17 * coords[1][1] + tmp17 * coords[3][1] +
                              tmp18 * coords[0][1] - tmp18 * coords[3][1] - tmp19 * coords[0][1] + tmp19 * coords[1][1] +
                              tmp20 * coords[1][1] - tmp20 * coords[2][1] - tmp21 * coords[0][1] + tmp21 * coords[2][1] +
                              tmp22 * coords[0][1] - tmp22 * coords[1][1];

               real_t tmp26 = tmp14 * coords[2][2] - tmp14 * coords[3][2] - tmp15 * coords[0][2] + tmp15 * coords[3][2] +
                              tmp16 * coords[0][2] - tmp16 * coords[2][2] - tmp17 * coords[1][2] + tmp17 * coords[3][2] +
                              tmp18 * coords[0][2] - tmp18 * coords[3][2] - tmp19 * coords[0][2] + tmp19 * coords[1][2] +
                              tmp20 * coords[1][2] - tmp20 * coords[2][2] - tmp21 * coords[0][2] + tmp21 * coords[2][2] +
                              tmp22 * coords[0][2] - tmp22 * coords[1][2];

               real_t tmp30 = SFD[i][0] * coords[0][1];
               real_t tmp31 = SFD[i][0] * coords[2][1];
               real_t tmp32 = SFD[i][0] * coords[3][1];
               real_t tmp33 = SFD[i][1] * coords[0][1];
               real_t tmp34 = SFD[i][1] * coords[1][1];
               real_t tmp35 = SFD[i][1] * coords[3][1];
               real_t tmp36 = SFD[i][2] * coords[0][1];
               real_t tmp37 = SFD[i][2] * coords[1][1];
               real_t tmp38 = SFD[i][2] * coords[2][1];

               real_t tmp39 = tmp30 * coords[2][2] - tmp30 * coords[3][2] - tmp31 * coords[0][2] + tmp31 * coords[3][2] +
                              tmp32 * coords[0][2] - tmp32 * coords[2][2] - tmp33 * coords[1][2] + tmp33 * coords[3][2] +
                              tmp34 * coords[0][2] - tmp34 * coords[3][2] - tmp35 * coords[0][2] + tmp35 * coords[1][2] +
                              tmp36 * coords[1][2] - tmp36 * coords[2][2] - tmp37 * coords[0][2] + tmp37 * coords[2][2] +
                              tmp38 * coords[0][2] - tmp38 * coords[1][2];

               real_t tmp40 = SFD[j][0] * coords[0][0];
               real_t tmp41 = SFD[j][0] * coords[2][0];
               real_t tmp42 = SFD[j][0] * coords[3][0];
               real_t tmp43 = SFD[j][1] * coords[0][0];
               real_t tmp44 = SFD[j][1] * coords[1][0];
               real_t tmp45 = SFD[j][1] * coords[3][0];
               real_t tmp46 = SFD[j][2] * coords[0][0];
               real_t tmp47 = SFD[j][2] * coords[1][0];
               real_t tmp48 = SFD[j][2] * coords[2][0];

               real_t tmp49 = tmp40 * coords[2][1] - tmp40 * coords[3][1] - tmp41 * coords[0][1] + tmp41 * coords[3][1] +
                              tmp42 * coords[0][1] - tmp42 * coords[2][1] - tmp43 * coords[1][1] + tmp43 * coords[3][1] +
                              tmp44 * coords[0][1] - tmp44 * coords[3][1] - tmp45 * coords[0][1] + tmp45 * coords[1][1] +
                              tmp46 * coords[1][1] - tmp46 * coords[2][1] - tmp47 * coords[0][1] + tmp47 * coords[2][1] +
                              tmp48 * coords[0][1] - tmp48 * coords[1][1];

               real_t tmp50 = tmp40 * coords[2][2] - tmp40 * coords[3][2] - tmp41 * coords[0][2] + tmp41 * coords[3][2] +
                              tmp42 * coords[0][2] - tmp42 * coords[2][2] - tmp43 * coords[1][2] + tmp43 * coords[3][2] +
                              tmp44 * coords[0][2] - tmp44 * coords[3][2] - tmp45 * coords[0][2] + tmp45 * coords[1][2] +
                              tmp46 * coords[1][2] - tmp46 * coords[2][2] - tmp47 * coords[0][2] + tmp47 * coords[2][2] +
                              tmp48 * coords[0][2] - tmp48 * coords[1][2];

               real_t tmp51 = SFD[j][0] * coords[0][1];
               real_t tmp52 = SFD[j][0] * coords[2][1];
               real_t tmp53 = SFD[j][0] * coords[3][1];
               real_t tmp54 = SFD[j][1] * coords[0][1];
               real_t tmp55 = SFD[j][1] * coords[1][1];
               real_t tmp56 = SFD[j][1] * coords[3][1];
               real_t tmp57 = SFD[j][2] * coords[0][1];
               real_t tmp58 = SFD[j][2] * coords[1][1];
               real_t tmp59 = SFD[j][2] * coords[2][1];

               real_t tmp60 = tmp51 * coords[2][2] - tmp51 * coords[3][2] - tmp52 * coords[0][2] + tmp52 * coords[3][2] +
                              tmp53 * coords[0][2] - tmp53 * coords[2][2] - tmp54 * coords[1][2] + tmp54 * coords[3][2] +
                              tmp55 * coords[0][2] - tmp55 * coords[3][2] - tmp56 * coords[0][2] + tmp56 * coords[1][2] +
                              tmp57 * coords[1][2] - tmp57 * coords[2][2] - tmp58 * coords[0][2] + tmp58 * coords[2][2] +
                              tmp59 * coords[0][2] - tmp59 * coords[1][2];

               real_t tmp61 = DPsi( 0, 1 ) * DPsi( 2, 0 ) - DPsi( 2, 1 ) * DPsi( 0, 0 );
               real_t tmp62 = DPsi( 0, 2 ) * DPsi( 2, 0 ) - DPsi( 2, 2 ) * DPsi( 0, 0 );
               real_t tmp63 = DPsi( 0, 1 ) * DPsi( 2, 2 ) - DPsi( 0, 2 ) * DPsi( 2, 1 );
               real_t tmp64 = DPsi( 1, 1 ) * DPsi( 2, 0 ) - DPsi( 2, 1 ) * DPsi( 1, 0 );
               real_t tmp65 = DPsi( 1, 2 ) * DPsi( 2, 0 ) - DPsi( 2, 2 ) * DPsi( 1, 0 );
               real_t tmp66 = DPsi( 1, 1 ) * DPsi( 2, 2 ) - DPsi( 1, 2 ) * DPsi( 2, 1 );

               real_t aux =
                   ( tmp13 * tmp23 + tmp25 * tmp26 - tmp29 * tmp39 ) * ( tmp13 * tmp49 + tmp25 * tmp50 - tmp29 * tmp60 ) +
                   ( tmp23 * tmp61 + tmp26 * tmp62 - tmp39 * tmp63 ) * ( tmp49 * tmp61 + tmp50 * tmp62 - tmp60 * tmp63 ) +
                   ( tmp23 * tmp64 + tmp26 * tmp65 - tmp39 * tmp66 ) * ( tmp49 * tmp64 + tmp50 * tmp65 - tmp60 * tmp66 );

               elMat( i, j ) += CUBAWEIGHTS[k] * detDPhiFac * detDPsiFac * aux;
            }
         }
#endif
      }

      // set lower triangular part from symmetry
      for ( uint i = 0; i < 10; i++ )
      {
         for ( uint j = 0; j < i; j++ )
         {
            elMat( i, j ) = elMat( j, i );
         }
      }
   };
};

} // namespace hyteg

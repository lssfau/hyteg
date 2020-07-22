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

namespace hyteg {

class P2Form_laplace : public P2FormHyTeG
{
 public:
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix6r& elMat ) const final
   {

// Derivatives of shape functions on unit triangle
#define DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE
#include "ShapeFunctionMacros.hpp"
#undef DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE

// Select quadrature rule
#define QUADPOINTS quadrature::D5_points
#define QUADWEIGHTS quadrature::D5_weights

// Executing quadrature rule
#define INTEGRATE2D( i, j )                                                                                     \
   elMat( i, j ) = 0.0;                                                                                         \
   for ( uint_t k = 0; k < QUADWEIGHTS.size(); k++ )                                                            \
   {                                                                                                            \
      real_t L2   = QUADPOINTS[k][0];                                                                           \
      real_t L3   = QUADPOINTS[k][1];                                                                           \
      real_t L1   = 1.0 - L2 - L3;                                                                              \
      mappedPt[0] = ( coords[1][0] - coords[0][0] ) * L2 + ( coords[2][0] - coords[0][0] ) * L3 + coords[0][0]; \
      mappedPt[1] = ( coords[1][1] - coords[0][1] ) * L2 + ( coords[2][1] - coords[0][1] ) * L3 + coords[0][1]; \
      Matrix2r DPsi;                                                                                            \
      geometryMap_->evalDF( mappedPt, DPsi );                                                                   \
                                                                                                                \
      real_t tmp4 = DPsi( 0, 1 ) * tmp3 + DPsi( 0, 0 ) * tmp0;                                                  \
      real_t tmp5 = DPsi( 0, 1 ) * tmp1 + DPsi( 0, 0 ) * tmp2;                                                  \
      real_t tmp6 = DPsi( 1, 0 ) * tmp2 + DPsi( 1, 1 ) * tmp1;                                                  \
      real_t tmp7 = DPsi( 1, 0 ) * tmp0 + DPsi( 1, 1 ) * tmp3;                                                  \
                                                                                                                \
      real_t aux3 = DPsi( 0, 1 ) * DPsi( 1, 0 ) - DPsi( 1, 1 ) * DPsi( 0, 0 );                                  \
      real_t aux4 = aux3 * aux3;                                                                                \
                                                                                                                \
      real_t detDPsi = std::abs( aux3 );                                                                        \
                                                                                                                \
      real_t aux5 = ( ( -DxiN##i * tmp5 + DetaN##i * tmp4 ) * ( -DxiN##j * tmp5 + DetaN##j * tmp4 ) +           \
                      ( DxiN##i * tmp6 - DetaN##i * tmp7 ) * ( DxiN##j * tmp6 - DetaN##j * tmp7 ) ) /           \
                    ( aux4 * aux2 );                                                                            \
                                                                                                                \
      elMat( i, j ) += QUADWEIGHTS[k] * detJacPhiInv * detDPsi * aux5;                                          \
   }                                                                                                            \
   elMat( j, i ) = elMat( i, j )

      // compute Jacobian determinant of inverse pull-back mapping
      real_t tmp0 = coords[0][0] - coords[1][0];
      real_t tmp1 = coords[0][1] - coords[2][1];
      real_t tmp2 = coords[0][0] - coords[2][0];
      real_t tmp3 = coords[0][1] - coords[1][1];

      real_t aux1 = tmp0 * tmp1 - tmp2 * tmp3;
      real_t aux2 = aux1 * aux1;

      real_t detJacPhiInv = std::abs( aux1 );

      // Quadrature point mapped to computational triangle
      Point3D mappedPt;

      // ------------
      //  Zeroth row
      // ------------

      INTEGRATE2D( 0, 0 );
      INTEGRATE2D( 0, 1 );
      INTEGRATE2D( 0, 2 );
      INTEGRATE2D( 0, 3 );
      INTEGRATE2D( 0, 4 );
      INTEGRATE2D( 0, 5 );

      // -----------
      //  First row
      // -----------
      INTEGRATE2D( 1, 1 );
      INTEGRATE2D( 1, 2 );
      INTEGRATE2D( 1, 3 );
      INTEGRATE2D( 1, 4 );
      INTEGRATE2D( 1, 5 );

      // ------------
      //  Second row
      // ------------
      INTEGRATE2D( 2, 2 );
      INTEGRATE2D( 2, 3 );
      INTEGRATE2D( 2, 4 );
      INTEGRATE2D( 2, 5 );

      // -----------
      //  Third row
      // -----------
      INTEGRATE2D( 3, 3 );
      INTEGRATE2D( 3, 4 );
      INTEGRATE2D( 3, 5 );

      // -----------
      //  Forth row
      // -----------
      INTEGRATE2D( 4, 4 );
      INTEGRATE2D( 4, 5 );

      // -----------
      //  Fifth row
      // -----------
      INTEGRATE2D( 5, 5 );

#define UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE
#include "ShapeFunctionMacros.hpp"
#undef UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE

#undef INTEGRATE2D
   };

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix10r& elMat ) const final
   {

// Derivatives of shape functions on unit triangle
#define DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET
#include "ShapeFunctionMacros.hpp"
#undef DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET

// Select cubature rule
#define CUBAPOINTS cubature::T3_points
#define CUBAWEIGHTS cubature::T3_weights

// #define P2_FORM_LAPLACE_WO_BLENDING
#ifdef P2_FORM_LAPLACE_WO_BLENDING
// Executing quadrature rule
#define INTEGRATE3D( i, j )                                                                                     \
   elMat( i, j ) = 0.0;                                                                                         \
   for ( uint_t k = 0; k < CUBAWEIGHTS.size(); k++ )                                                            \
   {                                                                                                            \
     real_t L2 = CUBAPOINTS[k][0];                                                                              \
     real_t L3 = CUBAPOINTS[k][1];                                                                              \
     real_t L4 = CUBAPOINTS[k][2];                                                                              \
     real_t L1 = 1.0 - L2 - L3 - L4;                                                                            \
     mappedPt[0] = ( coords[1][0] - coords[0][0] ) * L2 + ( coords[2][0] - coords[0][0] ) * L3 + ( coords[3][0] - coords[0][0] ) * L4 + coords[0][0]; \
     mappedPt[1] = ( coords[1][1] - coords[0][1] ) * L2 + ( coords[2][1] - coords[0][1] ) * L3 + ( coords[3][1] - coords[0][1] ) * L4 + coords[0][1]; \
     mappedPt[2] = ( coords[1][2] - coords[0][2] ) * L2 + ( coords[2][2] - coords[0][2] ) * L3 + ( coords[3][2] - coords[0][2] ) * L4 + coords[0][2]; \
                                                                                                                \
     real_t tmp27 = DX3N ## i*tmp26;                                                                            \
     real_t tmp31 = DX3N ## j*tmp26;                                                                            \
                                                                                                                \
     real_t aux = (tmp13*tmp18*(DX1N ## i*tmp30 + DX2N ## i*tmp25 - tmp27)*(DX1N ## j*tmp30 + DX2N ## j*tmp25 - tmp31) + tmp13*(DX1N ## i*tmp33 - DX2N ## i*tmp32 + tmp23*tmp27)*(DX1N ## j*tmp33 - DX2N ## j*tmp32 + tmp23*tmp31) + (-DX1N ## i*tmp36 + DX2N ## i*tmp35 - tmp27*tmp34)*(-DX1N ## j*tmp36 + DX2N ## j*tmp35 - tmp31*tmp34))/(tmp13*tmp18*pow(tmp24, 2));                                                          \
                                                                                                                \
     real_t detDPsi = 1.0;                                                                                      \
                                                                                                                \
     elMat( i, j ) += CUBAWEIGHTS[k] * detJacPhiInv * detDPsi * aux;                                            \
   }                                                                                                            \
   elMat( j, i ) = elMat( i, j )

      // compute Jacobian determinant of inverse pull-back mapping
      // and location independent auxilliary values
      real_t tmp0 = -coords[0][0];
      real_t tmp1 = tmp0 + coords[1][0];
      real_t tmp2 = -coords[0][1];
      real_t tmp3 = tmp2 + coords[2][1];
      real_t tmp4 = -coords[0][2];
      real_t tmp5 = tmp4 + coords[3][2];
      real_t tmp6 = tmp0 + coords[2][0];
      real_t tmp7 = tmp2 + coords[3][1];
      real_t tmp8 = tmp4 + coords[1][2];
      real_t tmp9 = tmp0 + coords[3][0];
      real_t tmp10 = tmp2 + coords[1][1];
      real_t tmp11 = tmp4 + coords[2][2];
      real_t tmp12 = coords[0][0] - coords[1][0];
      real_t tmp13 = tmp12 * tmp12;
      real_t tmp14 = tmp12*(coords[0][1] - coords[2][1]);
      real_t tmp15 = coords[0][0] - coords[2][0];
      real_t tmp16 = coords[0][1] - coords[1][1];
      real_t tmp17 = tmp14 - tmp15*tmp16;
      real_t tmp18 = tmp17 * tmp17;
      real_t tmp19 = coords[0][0] - coords[3][0];
      real_t tmp20 = coords[0][2] - coords[1][2];
      real_t tmp21 = tmp12*(coords[0][2] - coords[3][2]) - tmp19*tmp20;
      real_t tmp22 = tmp12*(coords[0][1] - coords[3][1]) - tmp16*tmp19;
      real_t tmp23 = tmp12*(coords[0][2] - coords[2][2]) - tmp15*tmp20;
      real_t tmp24 = tmp17*tmp21 - tmp22*tmp23;
      real_t tmp25 = tmp12*tmp22;
      real_t tmp26 = tmp12*tmp17;
      real_t tmp28 = tmp17*tmp19;
      real_t tmp29 = tmp15*tmp22;
      real_t tmp30 = tmp28 - tmp29;
      real_t tmp32 = tmp12*tmp17*tmp21;
      real_t tmp33 = tmp15*tmp24 + tmp23*(-tmp28 + tmp29);
      real_t tmp34 = tmp16*tmp23 - tmp17*tmp20;
      real_t tmp35 = tmp12*(tmp16*tmp24 + tmp22*tmp34);
      real_t tmp36 = tmp14*tmp24 - tmp30*tmp34;
      real_t detJacPhiInv = -tmp1*tmp11*tmp7 + tmp1*tmp3*tmp5 + tmp10*tmp11*tmp9 - tmp10*tmp5*tmp6 - tmp3*tmp8*tmp9 + tmp6*tmp7*tmp8;
      detJacPhiInv = std::abs( detJacPhiInv );

#else

// Executing quadrature rule
#define INTEGRATE3D( i, j )                                             \
      elMat( i, j ) = 0.0;                                              \
      for ( uint_t k = 0; k < CUBAWEIGHTS.size(); k++ )                 \
        {                                                               \
          real_t L2 = CUBAPOINTS[k][0];                                 \
          real_t L3 = CUBAPOINTS[k][1];                                 \
          real_t L4 = CUBAPOINTS[k][2];                                 \
          real_t L1 = 1.0 - L2 - L3 - L4;                               \
          mappedPt[0] = ( coords[1][0] - coords[0][0] ) * L2 + ( coords[2][0] - coords[0][0] ) * L3 + ( coords[3][0] - coords[0][0] ) * L4 + coords[0][0]; \
          mappedPt[1] = ( coords[1][1] - coords[0][1] ) * L2 + ( coords[2][1] - coords[0][1] ) * L3 + ( coords[3][1] - coords[0][1] ) * L4 + coords[0][1]; \
          mappedPt[2] = ( coords[1][2] - coords[0][2] ) * L2 + ( coords[2][2] - coords[0][2] ) * L3 + ( coords[3][2] - coords[0][2] ) * L4 + coords[0][2]; \
          Matrix3r DPsi;                                                \
          geometryMap_->evalDF( mappedPt, DPsi );                       \
                                                                        \
     real_t tmp17 = DPsi(0,1)*DPsi(1,0);                                \
     real_t tmp18 = DPsi(1,1)*DPsi(0,0) - tmp17;                        \
     real_t tmp19 = DPsi(0,2)*DPsi(2,0);                                \
     real_t tmp20 = DPsi(0,1)*DPsi(2,0);                                \
     real_t tmp21 = DPsi(2,1)*DPsi(0,0) - tmp20;                        \
     real_t tmp22 = -DPsi(0,2)*DPsi(1,0) + DPsi(1,2)*DPsi(0,0);         \
     real_t tmp23 = tmp21*tmp22;                                        \
     real_t tmp24 = tmp18*(DPsi(2,2)*DPsi(0,0) - tmp19) - tmp23;        \
     real_t tmp25 = 1.0/tmp24;                                          \
     real_t tmp26 = -tmp10*tmp11 + tmp13;                               \
     real_t tmp27 = tmp12*tmp3 - tmp7*tmp9;                             \
     real_t tmp28 = tmp16*(tmp1*tmp3 - tmp15) - tmp26*tmp27;            \
     real_t tmp29 = 1.0/tmp28;                                          \
     real_t tmp30 = DPsi(0,0)*tmp18*tmp25*tmp29*tmp3;                   \
     real_t tmp31 = tmp27*tmp3;                                         \
     real_t tmp32 = -DPsi(0,1)*tmp22 + DPsi(0,2)*tmp18;                 \
     real_t tmp33 = coords[0][1] - coords[1][1];                        \
     real_t tmp34 = tmp16*(coords[0][2] - coords[1][2]) - tmp27*tmp33;  \
     real_t tmp35 = DPsi(0,0)*tmp22*tmp25*tmp29*tmp31 + tmp16*tmp30 - tmp25*tmp29*tmp32*tmp34; \
     real_t tmp36 = tmp28*tmp3;                                         \
     real_t tmp37 = tmp26*tmp27*tmp3 + tmp36;                           \
     real_t tmp38 = 1.0/tmp16;                                          \
     real_t tmp39 = DPsi(0,0)*tmp22*tmp25*tmp29*tmp38;                  \
     real_t tmp40 = -tmp26*tmp34 + tmp28*tmp33;                         \
     real_t tmp41 = -tmp25*tmp29*tmp32*tmp38*tmp40 - tmp26*tmp30 - tmp37*tmp39; \
     real_t tmp42 = tmp10*tmp16 - tmp26*tmp7;                           \
     real_t tmp43 = 1.0/tmp3;                                           \
     real_t tmp44 = tmp31*tmp42 - tmp36*tmp7;                           \
     real_t tmp45 = tmp28*(tmp16 - tmp33*tmp7) - tmp34*tmp42;           \
     real_t tmp46 = -DPsi(0,0)*tmp18*tmp25*tmp29*tmp42 - tmp25*tmp29*tmp32*tmp38*tmp43*tmp45 - tmp39*tmp43*tmp44; \
     real_t tmp47 = DPsi(0,0)*tmp21*tmp25*tmp29*tmp3;                   \
     real_t tmp48 = DPsi(0,0)*tmp24;                                    \
     real_t tmp49 = DPsi(0,0)*tmp23 + tmp48;                            \
     real_t tmp50 = 1.0/tmp18;                                          \
     real_t tmp51 = tmp25*tmp27*tmp29*tmp3*tmp50;                       \
     real_t tmp52 = DPsi(0,0)*tmp21;                                    \
     real_t tmp53 = -DPsi(0,1)*tmp48 + tmp32*tmp52;                     \
     real_t tmp54 = 1.0/DPsi(0,0);                                      \
     real_t tmp55 = tmp25*tmp29*tmp34*tmp50*tmp54;                      \
     real_t tmp56 = -tmp16*tmp47 - tmp49*tmp51 + tmp53*tmp55;           \
     real_t tmp57 = tmp25*tmp29*tmp37*tmp38*tmp50;                      \
     real_t tmp58 = tmp25*tmp29*tmp38*tmp40*tmp50*tmp54;                \
     real_t tmp59 = tmp26*tmp47 + tmp49*tmp57 + tmp53*tmp58;            \
     real_t tmp60 = tmp25*tmp29*tmp42;                                  \
     real_t tmp61 = tmp25*tmp29*tmp38*tmp43*tmp44*tmp50;                \
     real_t tmp62 = tmp25*tmp29*tmp38*tmp43*tmp45*tmp50*tmp54;          \
     real_t tmp63 = tmp49*tmp61 + tmp52*tmp60 + tmp53*tmp62;            \
     real_t tmp64 = DPsi(1,0)*tmp21 - DPsi(2,0)*tmp18;                  \
     real_t tmp65 = tmp25*tmp29*tmp3*tmp64;                             \
     real_t tmp66 = -DPsi(1,0)*tmp24 - tmp22*tmp64;                     \
     real_t tmp67 = DPsi(1,1)*DPsi(0,0)*tmp24 - tmp32*tmp64;            \
     real_t tmp68 = tmp16*tmp65 - tmp51*tmp66 + tmp55*tmp67;            \
     real_t tmp69 = -tmp26*tmp65 + tmp57*tmp66 + tmp58*tmp67;           \
     real_t tmp70 = -tmp60*tmp64 + tmp61*tmp66 + tmp62*tmp67;           \
                                                                        \
     real_t aux = (DX1N ## i*tmp46 + DX2N ## i*tmp41 + DX3N ## i*tmp35)*(DX1N ## j*tmp46 + DX2N ## j*tmp41 + DX3N ## j*tmp35) + (DX1N ## i*tmp63 + DX2N ## i*tmp59 + DX3N ## i*tmp56)*(DX1N ## j*tmp63 + DX2N ## j*tmp59 + DX3N ## j*tmp56) + (DX1N ## i*tmp70 + DX2N ## i*tmp69 + DX3N ## i*tmp68)*(DX1N ## j*tmp70 + DX2N ## j*tmp69 + DX3N ## j*tmp68); \
                                                                                                                \
     real_t detDPsi = DPsi(0,2)*DPsi(2,1)*DPsi(1,0) + DPsi(1,1)*DPsi(2,2)*DPsi(0,0) - DPsi(1,1)*tmp19 - DPsi(1,2)*DPsi(2,1)*DPsi(0,0) + DPsi(1,2)*tmp20 - DPsi(2,2)*tmp17; \
     detDPsi = std::abs( detDPsi );                                     \
                                                                        \
     elMat( i, j ) += CUBAWEIGHTS[k] * detJacPhiInv * detDPsi * aux;    \
   }                                                                    \
   elMat( j, i ) = elMat( i, j )

      // compute Jacobian determinant of inverse pull-back mapping
      // and location independent auxilliary values
      real_t tmp0 = -coords[0][2];
      real_t tmp1 = tmp0 + coords[3][2];
      real_t tmp2 = -coords[0][0];
      real_t tmp3 = tmp2 + coords[1][0];
      real_t tmp4 = -coords[0][1];
      real_t tmp5 = tmp4 + coords[2][1];
      real_t tmp6 = tmp3*tmp5;
      real_t tmp7 = tmp2 + coords[2][0];
      real_t tmp8 = tmp4 + coords[3][1];
      real_t tmp9 = tmp0 + coords[1][2];
      real_t tmp10 = tmp2 + coords[3][0];
      real_t tmp11 = tmp4 + coords[1][1];
      real_t tmp12 = tmp0 + coords[2][2];
      real_t tmp13 = tmp3*tmp8;
      real_t tmp14 = tmp11*tmp7;
      real_t tmp15 = tmp10*tmp9;
      real_t tmp16 = -tmp14 + tmp6;

      real_t detJacPhiInv = -tmp1*tmp14 + tmp1*tmp6 + tmp10*tmp11*tmp12 - tmp12*tmp13 - tmp15*tmp5 + tmp7*tmp8*tmp9;
      detJacPhiInv = std::abs( detJacPhiInv );
#endif

      // Cubature point mapped to computational tetrahedron
      Point3D mappedPt;

      // dummy matrix for evaluation of Jacobian of 3D map
      Matrix3r dummy;

      // Zeroth row
      INTEGRATE3D( 0, 0 );
      INTEGRATE3D( 0, 1 );
      INTEGRATE3D( 0, 2 );
      INTEGRATE3D( 0, 3 );
      INTEGRATE3D( 0, 4 );
      INTEGRATE3D( 0, 5 );
      INTEGRATE3D( 0, 6 );
      INTEGRATE3D( 0, 7 );
      INTEGRATE3D( 0, 8 );
      INTEGRATE3D( 0, 9 );

      // First row
      INTEGRATE3D( 1, 1 );
      INTEGRATE3D( 1, 2 );
      INTEGRATE3D( 1, 3 );
      INTEGRATE3D( 1, 4 );
      INTEGRATE3D( 1, 5 );
      INTEGRATE3D( 1, 6 );
      INTEGRATE3D( 1, 7 );
      INTEGRATE3D( 1, 8 );
      INTEGRATE3D( 1, 9 );

      // Second row
      INTEGRATE3D( 2, 2 );
      INTEGRATE3D( 2, 3 );
      INTEGRATE3D( 2, 4 );
      INTEGRATE3D( 2, 5 );
      INTEGRATE3D( 2, 6 );
      INTEGRATE3D( 2, 7 );
      INTEGRATE3D( 2, 8 );
      INTEGRATE3D( 2, 9 );

      // Third row
      INTEGRATE3D( 3, 3 );
      INTEGRATE3D( 3, 4 );
      INTEGRATE3D( 3, 5 );
      INTEGRATE3D( 3, 6 );
      INTEGRATE3D( 3, 7 );
      INTEGRATE3D( 3, 8 );
      INTEGRATE3D( 3, 9 );

      // Fourth row
      INTEGRATE3D( 4, 4 );
      INTEGRATE3D( 4, 5 );
      INTEGRATE3D( 4, 6 );
      INTEGRATE3D( 4, 7 );
      INTEGRATE3D( 4, 8 );
      INTEGRATE3D( 4, 9 );

      // Fifth row
      INTEGRATE3D( 5, 5 );
      INTEGRATE3D( 5, 6 );
      INTEGRATE3D( 5, 7 );
      INTEGRATE3D( 5, 8 );
      INTEGRATE3D( 5, 9 );

      // Sixth row
      INTEGRATE3D( 6, 6 );
      INTEGRATE3D( 6, 7 );
      INTEGRATE3D( 6, 8 );
      INTEGRATE3D( 6, 9 );

      // Seventh row
      INTEGRATE3D( 7, 7 );
      INTEGRATE3D( 7, 8 );
      INTEGRATE3D( 7, 9 );

      // Eighth row
      INTEGRATE3D( 8, 8 );
      INTEGRATE3D( 8, 9 );

      // Ninth row
      INTEGRATE3D( 9, 9 );

#define UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET
#include "ShapeFunctionMacros.hpp"
#undef UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET

#undef INTEGRATE3D

   };
};

} // namespace hyteg

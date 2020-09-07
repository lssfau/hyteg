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
#define INTEGRATE3D( i, j )                          \
   elMat( i, j ) = 0.0;                              \
   for ( uint_t k = 0; k < CUBAWEIGHTS.size(); k++ ) \
   {                                                 \
     real_t L2 = CUBAPOINTS[k][0];                   \
     real_t L3 = CUBAPOINTS[k][1];                   \
     real_t L4 = CUBAPOINTS[k][2];                   \
     real_t L1 = 1.0 - L2 - L3 - L4;                 \
                                                     \
     real_t tmp12 = DX1N ## i*coords[0][0];          \
     real_t tmp13 = DX1N ## i*coords[2][0];          \
     real_t tmp14 = DX1N ## i*coords[3][0];          \
     real_t tmp15 = DX2N ## i*coords[0][0];          \
     real_t tmp16 = DX2N ## i*coords[1][0];          \
     real_t tmp17 = DX2N ## i*coords[3][0];          \
     real_t tmp18 = DX3N ## i*coords[0][0];          \
     real_t tmp19 = DX3N ## i*coords[1][0];          \
     real_t tmp20 = DX3N ## i*coords[2][0];          \
     real_t tmp21 = DX1N ## j*coords[0][0];          \
     real_t tmp22 = DX1N ## j*coords[2][0];          \
     real_t tmp23 = DX1N ## j*coords[3][0];          \
     real_t tmp24 = DX2N ## j*coords[0][0];          \
     real_t tmp25 = DX2N ## j*coords[1][0];          \
     real_t tmp26 = DX2N ## j*coords[3][0];          \
     real_t tmp27 = DX3N ## j*coords[0][0];          \
     real_t tmp28 = DX3N ## j*coords[1][0];          \
     real_t tmp29 = DX3N ## j*coords[2][0];          \
     real_t tmp30 = DX1N ## i*coords[0][1];          \
     real_t tmp31 = DX1N ## i*coords[2][1];          \
     real_t tmp32 = DX1N ## i*coords[3][1];          \
     real_t tmp33 = DX2N ## i*coords[0][1];          \
     real_t tmp34 = DX2N ## i*coords[1][1];          \
     real_t tmp35 = DX2N ## i*coords[3][1];          \
     real_t tmp36 = DX3N ## i*coords[0][1];          \
     real_t tmp37 = DX3N ## i*coords[1][1];          \
     real_t tmp38 = DX3N ## i*coords[2][1];          \
     real_t tmp39 = DX1N ## j*coords[0][1];          \
     real_t tmp40 = DX1N ## j*coords[2][1];          \
     real_t tmp41 = DX1N ## j*coords[3][1];          \
     real_t tmp42 = DX2N ## j*coords[0][1];          \
     real_t tmp43 = DX2N ## j*coords[1][1];          \
     real_t tmp44 = DX2N ## j*coords[3][1];          \
     real_t tmp45 = DX3N ## j*coords[0][1];          \
     real_t tmp46 = DX3N ## j*coords[1][1];          \
     real_t tmp47 = DX3N ## j*coords[2][1];          \
                                                     \
     real_t aux = (tmp12*coords[2][1] - tmp12*coords[3][1] - tmp13*coords[0][1] + tmp13*coords[3][1] + tmp14*coords[0][1] - tmp14*coords[2][1] - tmp15*coords[1][1] + tmp15*coords[3][1] + tmp16*coords[0][1] - tmp16*coords[3][1] - tmp17*coords[0][1] + tmp17*coords[1][1] + tmp18*coords[1][1] - tmp18*coords[2][1] - tmp19*coords[0][1] + tmp19*coords[2][1] + tmp20*coords[0][1] - tmp20*coords[1][1])*(tmp21*coords[2][1] - tmp21*coords[3][1] - tmp22*coords[0][1] + tmp22*coords[3][1] + tmp23*coords[0][1] - tmp23*coords[2][1] - tmp24*coords[1][1] + tmp24*coords[3][1] + tmp25*coords[0][1] - tmp25*coords[3][1] - tmp26*coords[0][1] + tmp26*coords[1][1] + tmp27*coords[1][1] - tmp27*coords[2][1] - tmp28*coords[0][1] + tmp28*coords[2][1] + tmp29*coords[0][1] - tmp29*coords[1][1]) + (tmp12*coords[2][2] - tmp12*coords[3][2] - tmp13*coords[0][2] + tmp13*coords[3][2] + tmp14*coords[0][2] - tmp14*coords[2][2] - tmp15*coords[1][2] + tmp15*coords[3][2] + tmp16*coords[0][2] - tmp16*coords[3][2] - tmp17*coords[0][2] + tmp17*coords[1][2] + tmp18*coords[1][2] - tmp18*coords[2][2] - tmp19*coords[0][2] + tmp19*coords[2][2] + tmp20*coords[0][2] - tmp20*coords[1][2])*(tmp21*coords[2][2] - tmp21*coords[3][2] - tmp22*coords[0][2] + tmp22*coords[3][2] + tmp23*coords[0][2] - tmp23*coords[2][2] - tmp24*coords[1][2] + tmp24*coords[3][2] + tmp25*coords[0][2] - tmp25*coords[3][2] - tmp26*coords[0][2] + tmp26*coords[1][2] + tmp27*coords[1][2] - tmp27*coords[2][2] - tmp28*coords[0][2] + tmp28*coords[2][2] + tmp29*coords[0][2] - tmp29*coords[1][2]) + (tmp30*coords[2][2] - tmp30*coords[3][2] - tmp31*coords[0][2] + tmp31*coords[3][2] + tmp32*coords[0][2] - tmp32*coords[2][2] - tmp33*coords[1][2] + tmp33*coords[3][2] + tmp34*coords[0][2] - tmp34*coords[3][2] - tmp35*coords[0][2] + tmp35*coords[1][2] + tmp36*coords[1][2] - tmp36*coords[2][2] - tmp37*coords[0][2] + tmp37*coords[2][2] + tmp38*coords[0][2] - tmp38*coords[1][2])*(tmp39*coords[2][2] - tmp39*coords[3][2] - tmp40*coords[0][2] + tmp40*coords[3][2] + tmp41*coords[0][2] - tmp41*coords[2][2] - tmp42*coords[1][2] + tmp42*coords[3][2] + tmp43*coords[0][2] - tmp43*coords[3][2] - tmp44*coords[0][2] + tmp44*coords[1][2] + tmp45*coords[1][2] - tmp45*coords[2][2] - tmp46*coords[0][2] + tmp46*coords[2][2] + tmp47*coords[0][2] - tmp47*coords[1][2]); \
                                                                        \
     elMat( i, j ) += CUBAWEIGHTS[k] * fac * aux;                       \
   }                                                                    \
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

  real_t detJacPhiInv = -tmp1*tmp11*tmp7 + tmp1*tmp3*tmp5 + tmp10*tmp11*tmp9 - tmp10*tmp5*tmp6 - tmp3*tmp8*tmp9 + tmp6*tmp7*tmp8;
  real_t fac = std::abs( detJacPhiInv ) / ( detJacPhiInv * detJacPhiInv ) ;

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
          real_t tmp12 = DPsi(0,1)*DPsi(1,0);                           \
          real_t tmp13 = -DPsi(1,1)*DPsi(0,0) + tmp12;                  \
          real_t tmp14 = DX1N ## i*coords[0][0];                        \
          real_t tmp15 = DX1N ## i*coords[2][0];                        \
          real_t tmp16 = DX1N ## i*coords[3][0];                        \
          real_t tmp17 = DX2N ## i*coords[0][0];                        \
          real_t tmp18 = DX2N ## i*coords[1][0];                        \
          real_t tmp19 = DX2N ## i*coords[3][0];                        \
          real_t tmp20 = DX3N ## i*coords[0][0];                        \
          real_t tmp21 = DX3N ## i*coords[1][0];                        \
          real_t tmp22 = DX3N ## i*coords[2][0];                        \
          real_t tmp23 = tmp14*coords[2][1] - tmp14*coords[3][1] - tmp15*coords[0][1] + tmp15*coords[3][1] + tmp16*coords[0][1] - tmp16*coords[2][1] - tmp17*coords[1][1] + tmp17*coords[3][1] + tmp18*coords[0][1] - tmp18*coords[3][1] - tmp19*coords[0][1] + tmp19*coords[1][1] + tmp20*coords[1][1] - tmp20*coords[2][1] - tmp21*coords[0][1] + tmp21*coords[2][1] + tmp22*coords[0][1] - tmp22*coords[1][1]; \
          real_t tmp24 = DPsi(1,2)*DPsi(0,0);                           \
          real_t tmp25 = DPsi(0,2)*DPsi(1,0) - tmp24;                   \
          real_t tmp26 = tmp14*coords[2][2] - tmp14*coords[3][2] - tmp15*coords[0][2] + tmp15*coords[3][2] + tmp16*coords[0][2] - tmp16*coords[2][2] - tmp17*coords[1][2] + tmp17*coords[3][2] + tmp18*coords[0][2] - tmp18*coords[3][2] - tmp19*coords[0][2] + tmp19*coords[1][2] + tmp20*coords[1][2] - tmp20*coords[2][2] - tmp21*coords[0][2] + tmp21*coords[2][2] + tmp22*coords[0][2] - tmp22*coords[1][2]; \
          real_t tmp27 = DPsi(0,1)*DPsi(1,2);                           \
          real_t tmp28 = DPsi(0,2)*DPsi(1,1);                           \
          real_t tmp29 = tmp27 - tmp28;                                 \
          real_t tmp30 = DX1N ## i*coords[0][1];                        \
          real_t tmp31 = DX1N ## i*coords[2][1];                        \
          real_t tmp32 = DX1N ## i*coords[3][1];                        \
          real_t tmp33 = DX2N ## i*coords[0][1];                        \
          real_t tmp34 = DX2N ## i*coords[1][1];                        \
          real_t tmp35 = DX2N ## i*coords[3][1];                        \
          real_t tmp36 = DX3N ## i*coords[0][1];                        \
          real_t tmp37 = DX3N ## i*coords[1][1];                        \
          real_t tmp38 = DX3N ## i*coords[2][1];                        \
          real_t tmp39 = tmp30*coords[2][2] - tmp30*coords[3][2] - tmp31*coords[0][2] + tmp31*coords[3][2] + tmp32*coords[0][2] - tmp32*coords[2][2] - tmp33*coords[1][2] + tmp33*coords[3][2] + tmp34*coords[0][2] - tmp34*coords[3][2] - tmp35*coords[0][2] + tmp35*coords[1][2] + tmp36*coords[1][2] - tmp36*coords[2][2] - tmp37*coords[0][2] + tmp37*coords[2][2] + tmp38*coords[0][2] - tmp38*coords[1][2]; \
          real_t tmp40 = DX1N ## j*coords[0][0];                        \
          real_t tmp41 = DX1N ## j*coords[2][0];                        \
          real_t tmp42 = DX1N ## j*coords[3][0];                        \
          real_t tmp43 = DX2N ## j*coords[0][0];                        \
          real_t tmp44 = DX2N ## j*coords[1][0];                        \
          real_t tmp45 = DX2N ## j*coords[3][0];                        \
          real_t tmp46 = DX3N ## j*coords[0][0];                        \
          real_t tmp47 = DX3N ## j*coords[1][0];                        \
          real_t tmp48 = DX3N ## j*coords[2][0];                        \
          real_t tmp49 = tmp40*coords[2][1] - tmp40*coords[3][1] - tmp41*coords[0][1] + tmp41*coords[3][1] + tmp42*coords[0][1] - tmp42*coords[2][1] - tmp43*coords[1][1] + tmp43*coords[3][1] + tmp44*coords[0][1] - tmp44*coords[3][1] - tmp45*coords[0][1] + tmp45*coords[1][1] + tmp46*coords[1][1] - tmp46*coords[2][1] - tmp47*coords[0][1] + tmp47*coords[2][1] + tmp48*coords[0][1] - tmp48*coords[1][1]; \
          real_t tmp50 = tmp40*coords[2][2] - tmp40*coords[3][2] - tmp41*coords[0][2] + tmp41*coords[3][2] + tmp42*coords[0][2] - tmp42*coords[2][2] - tmp43*coords[1][2] + tmp43*coords[3][2] + tmp44*coords[0][2] - tmp44*coords[3][2] - tmp45*coords[0][2] + tmp45*coords[1][2] + tmp46*coords[1][2] - tmp46*coords[2][2] - tmp47*coords[0][2] + tmp47*coords[2][2] + tmp48*coords[0][2] - tmp48*coords[1][2]; \
          real_t tmp51 = DX1N ## j*coords[0][1];                        \
          real_t tmp52 = DX1N ## j*coords[2][1];                        \
          real_t tmp53 = DX1N ## j*coords[3][1];                        \
          real_t tmp54 = DX2N ## j*coords[0][1];                        \
          real_t tmp55 = DX2N ## j*coords[1][1];                        \
          real_t tmp56 = DX2N ## j*coords[3][1];                        \
          real_t tmp57 = DX3N ## j*coords[0][1];                        \
          real_t tmp58 = DX3N ## j*coords[1][1];                        \
          real_t tmp59 = DX3N ## j*coords[2][1];                        \
          real_t tmp60 = tmp51*coords[2][2] - tmp51*coords[3][2] - tmp52*coords[0][2] + tmp52*coords[3][2] + tmp53*coords[0][2] - tmp53*coords[2][2] - tmp54*coords[1][2] + tmp54*coords[3][2] + tmp55*coords[0][2] - tmp55*coords[3][2] - tmp56*coords[0][2] + tmp56*coords[1][2] + tmp57*coords[1][2] - tmp57*coords[2][2] - tmp58*coords[0][2] + tmp58*coords[2][2] + tmp59*coords[0][2] - tmp59*coords[1][2]; \
          real_t tmp61 = DPsi(0,1)*DPsi(2,0) - DPsi(2,1)*DPsi(0,0);     \
          real_t tmp62 = DPsi(0,2)*DPsi(2,0) - DPsi(2,2)*DPsi(0,0);     \
          real_t tmp63 = DPsi(0,1)*DPsi(2,2) - DPsi(0,2)*DPsi(2,1);     \
          real_t tmp64 = DPsi(1,1)*DPsi(2,0) - DPsi(2,1)*DPsi(1,0);     \
          real_t tmp65 = DPsi(1,2)*DPsi(2,0) - DPsi(2,2)*DPsi(1,0);     \
          real_t tmp66 = DPsi(1,1)*DPsi(2,2) - DPsi(1,2)*DPsi(2,1);     \
                                                                        \
     real_t aux = (tmp13*tmp23 + tmp25*tmp26 - tmp29*tmp39)*(tmp13*tmp49 + tmp25*tmp50 - tmp29*tmp60) + (tmp23*tmp61 + tmp26*tmp62 - tmp39*tmp63)*(tmp49*tmp61 + tmp50*tmp62 - tmp60*tmp63) + (tmp23*tmp64 + tmp26*tmp65 - tmp39*tmp66)*(tmp49*tmp64 + tmp50*tmp65 - tmp60*tmp66); \
                                                                        \
     real_t detDPsi = DPsi(0,2)*DPsi(2,1)*DPsi(1,0) + DPsi(1,1)*DPsi(2,2)*DPsi(0,0) - DPsi(2,1)*tmp24 - DPsi(2,2)*tmp12 + DPsi(2,0)*tmp27 - DPsi(2,0)*tmp28; \
     real_t detDPsiFac = std::abs( detDPsi ) / ( detDPsi * detDPsi );   \
                                                                        \
     elMat( i, j ) += CUBAWEIGHTS[k] * detJacPhiFac * detDPsiFac * aux; \
   }                                                                    \
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

      real_t detJacPhiInv = -tmp1*tmp11*tmp7 + tmp1*tmp3*tmp5 + tmp10*tmp11*tmp9 - tmp10*tmp5*tmp6 - tmp3*tmp8*tmp9 + tmp6*tmp7*tmp8;
      real_t detJacPhiFac = std::abs( detJacPhiInv ) / ( detJacPhiInv * detJacPhiInv );
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

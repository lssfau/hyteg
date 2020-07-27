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

#include "hyteg/forms/form_hyteg_base/P2ToP1FormHyTeG.hpp"
#include "hyteg/forms/form_hyteg_manual/QuadratureRules.hpp"
#include "hyteg/geometry/GeometryMap.hpp"

namespace hyteg {

template < uint_t direction >
class P2ToP1Form_div : public P2ToP1FormHyTeG
{
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrixr< 3, 6 >& elMat ) const final
   {
      WALBERLA_ABORT( "P2ToP1Form_div not implemented for requested template value (direction)!" );
   };

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrixr< 4, 10 >& elMat ) const final
   {
      WALBERLA_ABORT( "P2ToP1Form_div not implemented for requested template value (direction)!" );
   };
};


// ========
//  div_x
// ========
template <>
class P2ToP1Form_div< 0 > : public P2ToP1FormHyTeG
{
 public:

   // ---------
   //  2D Case
   // ---------
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrixr< 3, 6 >& elMat ) const final
   {
// Set shape functions and their derivatives
#define DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE
#define DEFINE_P1_SHAPE_FUNCTIONS_TRIANGLE
#include "ShapeFunctionMacros.hpp"
#undef DEFINE_P1_SHAPE_FUNCTIONS_TRIANGLE
#undef DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE

// Select quadrature rule
#define QUADPOINTS quadrature::D5_points
#define QUADWEIGHTS quadrature::D5_weights

// Executing quadrature rule
#define INTEGRATE2D( i, j )                                                                                     \
   elMat( i, j ) = 0.0;                                                                                         \
   for ( uint_t k = 0; k < QUADWEIGHTS.size(); k++ )                                                            \
   {                                                                                                            \
      real_t L2 = QUADPOINTS[k][0];                                                                             \
      real_t L3 = QUADPOINTS[k][1];                                                                             \
      real_t L1 = 1.0 - L2 - L3;                                                                                \
                                                                                                                \
      mappedPt[0] = ( coords[1][0] - coords[0][0] ) * L2 + ( coords[2][0] - coords[0][0] ) * L3 + coords[0][0]; \
      mappedPt[1] = ( coords[1][1] - coords[0][1] ) * L2 + ( coords[2][1] - coords[0][1] ) * L3 + coords[0][1]; \
                                                                                                                \
      Matrix2r DPsi;                                                                                            \
      geometryMap_->evalDF( mappedPt, DPsi );                                                                   \
      real_t detDPsi = DPsi( 0, 0 ) * DPsi( 1, 1 ) - DPsi( 1, 0 ) * DPsi( 0, 1 );                               \
                                                                                                                \
      real_t aux4 = ( tmp3 * DetaN##j - tmp1 * DxiN##j ) * aux3;                                                \
      real_t aux5 = ( tmp2 * DxiN##j - tmp0 * DetaN##j ) * aux3;                                                \
                                                                                                                \
      real_t aux6 = ( DPsi( 1, 1 ) * aux4 - DPsi( 1, 0 ) * aux5 ) / detDPsi;                                    \
                                                                                                                \
      elMat( i, j ) -= QUADWEIGHTS[k] * detJacPhiInv * aux6 * std::abs( detDPsi ) * SF_M##i;                    \
   }

      // compute Jacobian determinant of inverse pull-back mapping
      real_t tmp0 = coords[0][0] - coords[1][0];
      real_t tmp1 = coords[0][1] - coords[2][1];
      real_t tmp2 = coords[0][0] - coords[2][0];
      real_t tmp3 = coords[0][1] - coords[1][1];

      real_t aux1 = tmp0 * tmp1 - tmp2 * tmp3;
      real_t aux3 = 1.0 / aux1;

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
      INTEGRATE2D( 1, 0 );
      INTEGRATE2D( 1, 1 );
      INTEGRATE2D( 1, 2 );
      INTEGRATE2D( 1, 3 );
      INTEGRATE2D( 1, 4 );
      INTEGRATE2D( 1, 5 );

      // ------------
      //  Second row
      // ------------
      INTEGRATE2D( 2, 0 );
      INTEGRATE2D( 2, 1 );
      INTEGRATE2D( 2, 2 );
      INTEGRATE2D( 2, 3 );
      INTEGRATE2D( 2, 4 );
      INTEGRATE2D( 2, 5 );

#define UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE
#define UNDEFINE_P1_SHAPE_FUNCTIONS_TRIANGLE
#include "ShapeFunctionMacros.hpp"
#undef UNDEFINE_P1_SHAPE_FUNCTIONS_TRIANGLE
#undef UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE

#undef INTEGRATE2D
   };


   // ---------
   //  3D Case
   // ---------
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrixr< 4, 10 >& elMat ) const final
   {


// Set shape functions and their derivatives
#define DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET
#define DEFINE_P1_SHAPE_FUNCTIONS_TET
#include "ShapeFunctionMacros.hpp"
#undef DEFINE_P1_SHAPE_FUNCTIONS_TET
#undef DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET

// Select cubature rule
#define CUBAPOINTS cubature::T4_points
#define CUBAWEIGHTS cubature::T4_weights

// Executing quadrature rule
#define INTEGRATE3D( i, j )                          \
   elMat( i, j ) = 0.0;                              \
   for ( uint_t k = 0; k < CUBAWEIGHTS.size(); k++ ) \
   {                                                 \
      real_t L2 = CUBAPOINTS[k][0];                  \
      real_t L3 = CUBAPOINTS[k][1];                  \
      real_t L4 = CUBAPOINTS[k][2];                  \
      real_t L1 = 1.0 - L2 - L3 - L4;                \
                                                     \
      mappedPt[0] = ( coords[1][0] - coords[0][0] ) * L2 + ( coords[2][0] - coords[0][0] ) * L3 + ( coords[3][0] - coords[0][0] ) * L4 + coords[0][0]; \
      mappedPt[1] = ( coords[1][1] - coords[0][1] ) * L2 + ( coords[2][1] - coords[0][1] ) * L3 + ( coords[3][1] - coords[0][1] ) * L4 + coords[0][1]; \
      mappedPt[2] = ( coords[1][2] - coords[0][2] ) * L2 + ( coords[2][2] - coords[0][2] ) * L3 + ( coords[3][2] - coords[0][2] ) * L4 + coords[0][2]; \
                                                     \
      Matrix3r DPsi;                                 \
      geometryMap_->evalDF( mappedPt, DPsi );        \
                                                     \
      real_t tmp12 = DPsi(1,1)*DPsi(2,2);            \
      real_t tmp13 = DPsi(1,2)*DPsi(2,0);            \
      real_t tmp14 = DPsi(1,0)*DPsi(2,1);            \
      real_t tmp15 = DPsi(1,2)*DPsi(2,1);            \
      real_t tmp16 = DPsi(1,0)*DPsi(2,2);            \
      real_t tmp17 = DPsi(1,1)*DPsi(2,0);            \
      real_t tmp18 = tmp14 - tmp17;                  \
      real_t tmp23 = tmp13 - tmp16;                  \
      real_t tmp28 = tmp12 - tmp15;                  \
                                                     \
      real_t diff = DX1N ## j*(tmp18*(tmp19 + tmp20 - tmp21 - tmp22 + coords[2][0]*coords[3][1] - coords[3][0]*coords[2][1]) + tmp23*(tmp24 + tmp25 - tmp26 - tmp27 - coords[2][0]*coords[3][2] + coords[3][0]*coords[2][2]) + tmp28*(tmp29 + tmp30 - tmp31 - coords[2][1]*coords[0][2] + coords[2][1]*coords[3][2] - coords[3][1]*coords[2][2])) + DX2N ## j*(tmp18*(-tmp20 + tmp21 + tmp32 - tmp33 - coords[1][0]*coords[3][1] + coords[3][0]*coords[1][1]) + tmp23*(-tmp24 + tmp27 + tmp34 - tmp35 + coords[1][0]*coords[3][2] - coords[3][0]*coords[1][2]) + tmp28*(-tmp30 + tmp31 + tmp36 - tmp37 - coords[1][1]*coords[3][2] + coords[3][1]*coords[1][2])) + DX3N ## j*(tmp18*(-tmp19 + tmp22 - tmp32 + tmp33 + coords[1][0]*coords[2][1] - coords[2][0]*coords[1][1]) + tmp23*(-tmp25 + tmp26 - tmp34 + tmp35 - coords[1][0]*coords[2][2] + coords[2][0]*coords[1][2]) + tmp28*(-tmp29 - tmp36 + tmp37 + coords[1][1]*coords[2][2] + coords[2][1]*coords[0][2] - coords[2][1]*coords[1][2])); \
                                                     \
      real_t detQ = DPsi(0,0)*tmp12 - DPsi(0,0)*tmp15 + DPsi(0,1)*tmp13 - DPsi(0,1)*tmp16 + DPsi(0,2)*tmp14 - DPsi(0,2)*tmp17; \
      real_t detQfac = std::abs( detQ ) / detQ;      \
      elMat( i, j ) -= CUBAWEIGHTS[k] * detPfac * detQfac * diff * SF_M##i; \
   }

     // Quadrature point mapped to computational triangle
     Point3D mappedPt;

     // pre-compute values including Jacobian determinant of inverse pull-back mapping
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
     real_t tmp19 = coords[0][0]*coords[2][1];
     real_t tmp20 = coords[3][0]*coords[0][1];
     real_t tmp21 = coords[0][0]*coords[3][1];
     real_t tmp22 = coords[2][0]*coords[0][1];
     real_t tmp24 = coords[0][0]*coords[3][2];
     real_t tmp25 = coords[2][0]*coords[0][2];
     real_t tmp26 = coords[0][0]*coords[2][2];
     real_t tmp27 = coords[3][0]*coords[0][2];
     real_t tmp29 = coords[0][1]*coords[2][2];
     real_t tmp30 = coords[3][1]*coords[0][2];
     real_t tmp31 = coords[0][1]*coords[3][2];
     real_t tmp32 = coords[1][0]*coords[0][1];
     real_t tmp33 = coords[0][0]*coords[1][1];
     real_t tmp34 = coords[0][0]*coords[1][2];
     real_t tmp35 = coords[1][0]*coords[0][2];
     real_t tmp36 = coords[1][1]*coords[0][2];
     real_t tmp37 = coords[0][1]*coords[1][2];

     real_t detP = -tmp1*tmp11*tmp7 + tmp1*tmp3*tmp5 + tmp10*tmp11*tmp9 - tmp10*tmp5*tmp6 - tmp3*tmp8*tmp9 + tmp6*tmp7*tmp8;
     real_t detPfac = std::abs( detP ) / detP;

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

     INTEGRATE3D( 1, 0 );
     INTEGRATE3D( 1, 1 );
     INTEGRATE3D( 1, 2 );
     INTEGRATE3D( 1, 3 );
     INTEGRATE3D( 1, 4 );
     INTEGRATE3D( 1, 5 );
     INTEGRATE3D( 1, 6 );
     INTEGRATE3D( 1, 7 );
     INTEGRATE3D( 1, 8 );
     INTEGRATE3D( 1, 9 );

     INTEGRATE3D( 2, 0 );
     INTEGRATE3D( 2, 1 );
     INTEGRATE3D( 2, 2 );
     INTEGRATE3D( 2, 3 );
     INTEGRATE3D( 2, 4 );
     INTEGRATE3D( 2, 5 );
     INTEGRATE3D( 2, 6 );
     INTEGRATE3D( 2, 7 );
     INTEGRATE3D( 2, 8 );
     INTEGRATE3D( 2, 9 );

     INTEGRATE3D( 3, 0 );
     INTEGRATE3D( 3, 1 );
     INTEGRATE3D( 3, 2 );
     INTEGRATE3D( 3, 3 );
     INTEGRATE3D( 3, 4 );
     INTEGRATE3D( 3, 5 );
     INTEGRATE3D( 3, 6 );
     INTEGRATE3D( 3, 7 );
     INTEGRATE3D( 3, 8 );
     INTEGRATE3D( 3, 9 );

// Undefine shape functions
#define UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET
#define UNDEFINE_P1_SHAPE_FUNCTIONS_TET
#include "ShapeFunctionMacros.hpp"
#undef UNDEFINE_P1_SHAPE_FUNCTIONS_TET
#undef UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET

   };
};


// ========
//  div_y
// ========
template <>
class P2ToP1Form_div< 1 > : public P2ToP1FormHyTeG
{
 public:

   // ---------
   //  2D Case
   // ---------
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrixr< 3, 6 >& elMat ) const final
   {

// Derivatives of P2 shape functions on unit triangle
#define DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE
#include "ShapeFunctionMacros.hpp"
#undef DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE

// P1 shape functions on unit triangle
#define DEFINE_P1_SHAPE_FUNCTIONS_TRIANGLE
#include "ShapeFunctionMacros.hpp"
#undef DEFINE_P1_SHAPE_FUNCTIONS_TRIANGLE

// Select quadrature rule
#define QUADPOINTS quadrature::D5_points
#define QUADWEIGHTS quadrature::D5_weights

// Executing quadrature rule
#define INTEGRATE2D( i, j )                                                                                     \
   elMat( i, j ) = 0.0;                                                                                         \
   for ( uint_t k = 0; k < QUADWEIGHTS.size(); k++ )                                                            \
   {                                                                                                            \
      real_t L2 = QUADPOINTS[k][0];                                                                             \
      real_t L3 = QUADPOINTS[k][1];                                                                             \
      real_t L1 = 1.0 - L2 - L3;                                                                                \
                                                                                                                \
      mappedPt[0] = ( coords[1][0] - coords[0][0] ) * L2 + ( coords[2][0] - coords[0][0] ) * L3 + coords[0][0]; \
      mappedPt[1] = ( coords[1][1] - coords[0][1] ) * L2 + ( coords[2][1] - coords[0][1] ) * L3 + coords[0][1]; \
                                                                                                                \
      Matrix2r DPsi;                                                                                            \
      geometryMap_->evalDF( mappedPt, DPsi );                                                                   \
      real_t detDPsi = DPsi( 0, 0 ) * DPsi( 1, 1 ) - DPsi( 1, 0 ) * DPsi( 0, 1 );                               \
                                                                                                                \
      real_t aux4 = ( tmp3 * DetaN##j - tmp1 * DxiN##j ) * aux3;                                                \
      real_t aux5 = ( tmp2 * DxiN##j - tmp0 * DetaN##j ) * aux3;                                                \
                                                                                                                \
      real_t aux6 = ( DPsi( 0, 0 ) * aux5 - DPsi( 0, 1 ) * aux4 ) / detDPsi;                                    \
                                                                                                                \
      elMat( i, j ) -= QUADWEIGHTS[k] * detJacPhiInv * aux6 * std::abs( detDPsi ) * SF_M##i;                    \
   }

      // compute Jacobian determinant of inverse pull-back mapping
      real_t tmp0 = coords[0][0] - coords[1][0];
      real_t tmp1 = coords[0][1] - coords[2][1];
      real_t tmp2 = coords[0][0] - coords[2][0];
      real_t tmp3 = coords[0][1] - coords[1][1];

      real_t aux1 = tmp0 * tmp1 - tmp2 * tmp3;
      real_t aux3 = 1.0 / aux1;

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
      INTEGRATE2D( 1, 0 );
      INTEGRATE2D( 1, 1 );
      INTEGRATE2D( 1, 2 );
      INTEGRATE2D( 1, 3 );
      INTEGRATE2D( 1, 4 );
      INTEGRATE2D( 1, 5 );

      // ------------
      //  Second row
      // ------------
      INTEGRATE2D( 2, 0 );
      INTEGRATE2D( 2, 1 );
      INTEGRATE2D( 2, 2 );
      INTEGRATE2D( 2, 3 );
      INTEGRATE2D( 2, 4 );
      INTEGRATE2D( 2, 5 );

#define UNDEFINE_P1_SHAPE_FUNCTIONS_TRIANGLE
#define UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE
#include "ShapeFunctionMacros.hpp"
#undef UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE
#undef UNDEFINE_P1_SHAPE_FUNCTIONS_TRIANGLE

#undef INTEGRATE2D
   };


   // ---------
   //  3D Case
   // ---------
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrixr< 4, 10 >& elMat ) const final
   {
     // WALBERLA_ABORT( "P2ToP1Form_div<1> not implemented for 3D, yet." );

// Set shape functions and their derivatives
#define DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET
#define DEFINE_P1_SHAPE_FUNCTIONS_TET
#include "ShapeFunctionMacros.hpp"
#undef DEFINE_P1_SHAPE_FUNCTIONS_TET
#undef DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET

// Select cubature rule
#define CUBAPOINTS cubature::T4_points
#define CUBAWEIGHTS cubature::T4_weights

// Executing quadrature rule
#define INTEGRATE3D( i, j )                          \
   elMat( i, j ) = 0.0;                              \
   for ( uint_t k = 0; k < CUBAWEIGHTS.size(); k++ ) \
   {                                                 \
      real_t L2 = CUBAPOINTS[k][0];                  \
      real_t L3 = CUBAPOINTS[k][1];                  \
      real_t L4 = CUBAPOINTS[k][2];                  \
      real_t L1 = 1.0 - L2 - L3 - L4;                \
                                                     \
      mappedPt[0] = ( coords[1][0] - coords[0][0] ) * L2 + ( coords[2][0] - coords[0][0] ) * L3 + ( coords[3][0] - coords[0][0] ) * L4 + coords[0][0]; \
      mappedPt[1] = ( coords[1][1] - coords[0][1] ) * L2 + ( coords[2][1] - coords[0][1] ) * L3 + ( coords[3][1] - coords[0][1] ) * L4 + coords[0][1]; \
      mappedPt[2] = ( coords[1][2] - coords[0][2] ) * L2 + ( coords[2][2] - coords[0][2] ) * L3 + ( coords[3][2] - coords[0][2] ) * L4 + coords[0][2]; \
                                                     \
      Matrix3r DPsi;                                 \
      geometryMap_->evalDF( mappedPt, DPsi );        \
                                                     \
      real_t tmp12 = DPsi(0,0)*DPsi(2,2);            \
      real_t tmp13 = DPsi(0,1)*DPsi(2,0);            \
      real_t tmp14 = DPsi(0,2)*DPsi(2,1);            \
      real_t tmp15 = DPsi(0,0)*DPsi(2,1);            \
      real_t tmp16 = DPsi(0,1)*DPsi(2,2);            \
      real_t tmp17 = DPsi(0,2)*DPsi(2,0);            \
      real_t tmp18 = tmp13 - tmp15;                  \
      real_t tmp23 = tmp12 - tmp17;                  \
      real_t tmp28 = tmp14 - tmp16;                  \
                                                     \
      real_t diff = DX1N ## j*(tmp18*(tmp19 + tmp20 - tmp21 - tmp22 + coords[2][0]*coords[3][1] - coords[3][0]*coords[2][1]) + tmp23*(tmp24 + tmp25 - tmp26 - tmp27 - coords[2][0]*coords[3][2] + coords[3][0]*coords[2][2]) + tmp28*(tmp29 + tmp30 - tmp31 - coords[2][1]*coords[0][2] + coords[2][1]*coords[3][2] - coords[3][1]*coords[2][2])) + DX2N ## j*(tmp18*(-tmp20 + tmp21 + tmp32 - tmp33 - coords[1][0]*coords[3][1] + coords[3][0]*coords[1][1]) + tmp23*(-tmp24 + tmp27 + tmp34 - tmp35 + coords[1][0]*coords[3][2] - coords[3][0]*coords[1][2]) + tmp28*(-tmp30 + tmp31 + tmp36 - tmp37 - coords[1][1]*coords[3][2] + coords[3][1]*coords[1][2])) + DX3N ## j*(tmp18*(-tmp19 + tmp22 - tmp32 + tmp33 + coords[1][0]*coords[2][1] - coords[2][0]*coords[1][1]) + tmp23*(-tmp25 + tmp26 - tmp34 + tmp35 - coords[1][0]*coords[2][2] + coords[2][0]*coords[1][2]) + tmp28*(-tmp29 - tmp36 + tmp37 + coords[1][1]*coords[2][2] + coords[2][1]*coords[0][2] - coords[2][1]*coords[1][2])); \
                                                     \
      real_t detQ = DPsi(1,0)*tmp14 - DPsi(1,0)*tmp16 + DPsi(1,1)*tmp12 - DPsi(1,1)*tmp17 + DPsi(1,2)*tmp13 - DPsi(1,2)*tmp15; \
      real_t detQfac = std::abs( detQ ) / detQ;      \
      elMat( i, j ) -= CUBAWEIGHTS[k] * detPfac * detQfac * diff * SF_M##i; \
   }

     // Quadrature point mapped to computational triangle
     Point3D mappedPt;

     // pre-compute values including Jacobian determinant of inverse pull-back mapping
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
     real_t tmp19 = coords[0][0]*coords[2][1];
     real_t tmp20 = coords[3][0]*coords[0][1];
     real_t tmp21 = coords[0][0]*coords[3][1];
     real_t tmp22 = coords[2][0]*coords[0][1];
     real_t tmp24 = coords[0][0]*coords[3][2];
     real_t tmp25 = coords[2][0]*coords[0][2];
     real_t tmp26 = coords[0][0]*coords[2][2];
     real_t tmp27 = coords[3][0]*coords[0][2];
     real_t tmp29 = coords[0][1]*coords[2][2];
     real_t tmp30 = coords[3][1]*coords[0][2];
     real_t tmp31 = coords[0][1]*coords[3][2];
     real_t tmp32 = coords[1][0]*coords[0][1];
     real_t tmp33 = coords[0][0]*coords[1][1];
     real_t tmp34 = coords[0][0]*coords[1][2];
     real_t tmp35 = coords[1][0]*coords[0][2];
     real_t tmp36 = coords[1][1]*coords[0][2];
     real_t tmp37 = coords[0][1]*coords[1][2];

     real_t detP = -tmp1*tmp11*tmp7 + tmp1*tmp3*tmp5 + tmp10*tmp11*tmp9 - tmp10*tmp5*tmp6 - tmp3*tmp8*tmp9 + tmp6*tmp7*tmp8;
     real_t detPfac = std::abs( detP ) / detP;

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

     INTEGRATE3D( 1, 0 );
     INTEGRATE3D( 1, 1 );
     INTEGRATE3D( 1, 2 );
     INTEGRATE3D( 1, 3 );
     INTEGRATE3D( 1, 4 );
     INTEGRATE3D( 1, 5 );
     INTEGRATE3D( 1, 6 );
     INTEGRATE3D( 1, 7 );
     INTEGRATE3D( 1, 8 );
     INTEGRATE3D( 1, 9 );

     INTEGRATE3D( 2, 0 );
     INTEGRATE3D( 2, 1 );
     INTEGRATE3D( 2, 2 );
     INTEGRATE3D( 2, 3 );
     INTEGRATE3D( 2, 4 );
     INTEGRATE3D( 2, 5 );
     INTEGRATE3D( 2, 6 );
     INTEGRATE3D( 2, 7 );
     INTEGRATE3D( 2, 8 );
     INTEGRATE3D( 2, 9 );

     INTEGRATE3D( 3, 0 );
     INTEGRATE3D( 3, 1 );
     INTEGRATE3D( 3, 2 );
     INTEGRATE3D( 3, 3 );
     INTEGRATE3D( 3, 4 );
     INTEGRATE3D( 3, 5 );
     INTEGRATE3D( 3, 6 );
     INTEGRATE3D( 3, 7 );
     INTEGRATE3D( 3, 8 );
     INTEGRATE3D( 3, 9 );

// Undefine shape functions
#define UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET
#define UNDEFINE_P1_SHAPE_FUNCTIONS_TET
#include "ShapeFunctionMacros.hpp"
#undef UNDEFINE_P1_SHAPE_FUNCTIONS_TET
#undef UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET

   };
};



// ========
//  div_z
// ========
template <>
class P2ToP1Form_div< 2 > : public P2ToP1FormHyTeG
{
 public:

   // ---------
   //  2D Case
   // ---------
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrixr< 3, 6 >& elMat ) const final
   {
     WALBERLA_ABORT( "P2ToP1Form_div< 2 >.integrateAll() not available for 3D!" );
   };


   // ---------
   //  3D Case
   // ---------
   void integrateAll( const std::array< Point3D, 4 >& coords, Matrixr< 4, 10 >& elMat ) const final
   {

// Set shape functions and their derivatives
#define DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET
#define DEFINE_P1_SHAPE_FUNCTIONS_TET
#include "ShapeFunctionMacros.hpp"
#undef DEFINE_P1_SHAPE_FUNCTIONS_TET
#undef DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET

// Select cubature rule
#define CUBAPOINTS cubature::T4_points
#define CUBAWEIGHTS cubature::T4_weights

// Executing quadrature rule
#define INTEGRATE3D( i, j )                          \
   elMat( i, j ) = 0.0;                              \
   for ( uint_t k = 0; k < CUBAWEIGHTS.size(); k++ ) \
   {                                                 \
      real_t L2 = CUBAPOINTS[k][0];                  \
      real_t L3 = CUBAPOINTS[k][1];                  \
      real_t L4 = CUBAPOINTS[k][2];                  \
      real_t L1 = 1.0 - L2 - L3 - L4;                \
                                                     \
      mappedPt[0] = ( coords[1][0] - coords[0][0] ) * L2 + ( coords[2][0] - coords[0][0] ) * L3 + ( coords[3][0] - coords[0][0] ) * L4 + coords[0][0]; \
      mappedPt[1] = ( coords[1][1] - coords[0][1] ) * L2 + ( coords[2][1] - coords[0][1] ) * L3 + ( coords[3][1] - coords[0][1] ) * L4 + coords[0][1]; \
      mappedPt[2] = ( coords[1][2] - coords[0][2] ) * L2 + ( coords[2][2] - coords[0][2] ) * L3 + ( coords[3][2] - coords[0][2] ) * L4 + coords[0][2]; \
                                                     \
      Matrix3r DPsi;                                 \
      geometryMap_->evalDF( mappedPt, DPsi );        \
                                                     \
      real_t tmp12 = DPsi(0,0)*DPsi(1,1);            \
      real_t tmp13 = DPsi(0,1)*DPsi(1,2);            \
      real_t tmp14 = DPsi(0,2)*DPsi(1,0);            \
      real_t tmp15 = DPsi(0,0)*DPsi(1,2);            \
      real_t tmp16 = DPsi(0,1)*DPsi(1,0);            \
      real_t tmp17 = DPsi(0,2)*DPsi(1,1);            \
      real_t tmp18 = tmp12 - tmp16;                  \
      real_t tmp23 = tmp14 - tmp15;                  \
      real_t tmp28 = tmp13 - tmp17;                  \
                                                     \
      real_t diff = DX1N ## j*(tmp18*(tmp19 + tmp20 - tmp21 - tmp22 + coords[2][0]*coords[3][1] - coords[3][0]*coords[2][1]) + tmp23*(tmp24 + tmp25 - tmp26 - tmp27 - coords[2][0]*coords[3][2] + coords[3][0]*coords[2][2]) + tmp28*(tmp29 + tmp30 - tmp31 - coords[2][1]*coords[0][2] + coords[2][1]*coords[3][2] - coords[3][1]*coords[2][2])) + DX2N ## j*(tmp18*(-tmp20 + tmp21 + tmp32 - tmp33 - coords[1][0]*coords[3][1] + coords[3][0]*coords[1][1]) + tmp23*(-tmp24 + tmp27 + tmp34 - tmp35 + coords[1][0]*coords[3][2] - coords[3][0]*coords[1][2]) + tmp28*(-tmp30 + tmp31 + tmp36 - tmp37 - coords[1][1]*coords[3][2] + coords[3][1]*coords[1][2])) + DX3N ## j*(tmp18*(-tmp19 + tmp22 - tmp32 + tmp33 + coords[1][0]*coords[2][1] - coords[2][0]*coords[1][1]) + tmp23*(-tmp25 + tmp26 - tmp34 + tmp35 - coords[1][0]*coords[2][2] + coords[2][0]*coords[1][2]) + tmp28*(-tmp29 - tmp36 + tmp37 + coords[1][1]*coords[2][2] + coords[2][1]*coords[0][2] - coords[2][1]*coords[1][2])); \
                                                     \
      real_t detQ = DPsi(2,0)*tmp13 - DPsi(2,0)*tmp17 + DPsi(2,1)*tmp14 - DPsi(2,1)*tmp15 + DPsi(2,2)*tmp12 - DPsi(2,2)*tmp16; \
      real_t detQfac = std::abs( detQ ) / detQ;      \
      elMat( i, j ) -= CUBAWEIGHTS[k] * detPfac * detQfac * diff * SF_M##i; \
   }

     // Quadrature point mapped to computational triangle
     Point3D mappedPt;

     // pre-compute values including Jacobian determinant of inverse pull-back mapping
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
     real_t tmp19 = coords[0][0]*coords[2][1];
     real_t tmp20 = coords[3][0]*coords[0][1];
     real_t tmp21 = coords[0][0]*coords[3][1];
     real_t tmp22 = coords[2][0]*coords[0][1];
     real_t tmp24 = coords[0][0]*coords[3][2];
     real_t tmp25 = coords[2][0]*coords[0][2];
     real_t tmp26 = coords[0][0]*coords[2][2];
     real_t tmp27 = coords[3][0]*coords[0][2];
     real_t tmp29 = coords[0][1]*coords[2][2];
     real_t tmp30 = coords[3][1]*coords[0][2];
     real_t tmp31 = coords[0][1]*coords[3][2];
     real_t tmp32 = coords[1][0]*coords[0][1];
     real_t tmp33 = coords[0][0]*coords[1][1];
     real_t tmp34 = coords[0][0]*coords[1][2];
     real_t tmp35 = coords[1][0]*coords[0][2];
     real_t tmp36 = coords[1][1]*coords[0][2];
     real_t tmp37 = coords[0][1]*coords[1][2];

     real_t detP = -tmp1*tmp11*tmp7 + tmp1*tmp3*tmp5 + tmp10*tmp11*tmp9 - tmp10*tmp5*tmp6 - tmp3*tmp8*tmp9 + tmp6*tmp7*tmp8;
     real_t detPfac = std::abs( detP ) / detP;

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

     INTEGRATE3D( 1, 0 );
     INTEGRATE3D( 1, 1 );
     INTEGRATE3D( 1, 2 );
     INTEGRATE3D( 1, 3 );
     INTEGRATE3D( 1, 4 );
     INTEGRATE3D( 1, 5 );
     INTEGRATE3D( 1, 6 );
     INTEGRATE3D( 1, 7 );
     INTEGRATE3D( 1, 8 );
     INTEGRATE3D( 1, 9 );

     INTEGRATE3D( 2, 0 );
     INTEGRATE3D( 2, 1 );
     INTEGRATE3D( 2, 2 );
     INTEGRATE3D( 2, 3 );
     INTEGRATE3D( 2, 4 );
     INTEGRATE3D( 2, 5 );
     INTEGRATE3D( 2, 6 );
     INTEGRATE3D( 2, 7 );
     INTEGRATE3D( 2, 8 );
     INTEGRATE3D( 2, 9 );

     INTEGRATE3D( 3, 0 );
     INTEGRATE3D( 3, 1 );
     INTEGRATE3D( 3, 2 );
     INTEGRATE3D( 3, 3 );
     INTEGRATE3D( 3, 4 );
     INTEGRATE3D( 3, 5 );
     INTEGRATE3D( 3, 6 );
     INTEGRATE3D( 3, 7 );
     INTEGRATE3D( 3, 8 );
     INTEGRATE3D( 3, 9 );

// Undefine shape functions
#define UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET
#define UNDEFINE_P1_SHAPE_FUNCTIONS_TET
#include "ShapeFunctionMacros.hpp"
#undef UNDEFINE_P1_SHAPE_FUNCTIONS_TET
#undef UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET

   };
};

} // namespace hyteg

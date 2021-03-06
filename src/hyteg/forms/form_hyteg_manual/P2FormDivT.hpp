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
class P2Form_divt : public P2FormHyTeG
{
 public:
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix6r& elMat ) const final
   {
      WALBERLA_ABORT( "P2FormDivT not implemented for requested template value (direction)!" );
   };

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix10r& elMat ) const final
   {
      WALBERLA_ABORT( "P2FormDivT not implemented for requested template value (direction)!" );
   };
};

// --------
//  div_x
// --------
template <>
class P2Form_divt< 0 > : public P2FormHyTeG
{
 public:
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix6r& elMat ) const final
   {

// Derivatives of P2 shape functions on unit triangle
#define DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE
#include "ShapeFunctionMacros.hpp"
#undef DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE

// P2 shape functions on unit triangle
#define DEFINE_P2_SHAPE_FUNCTIONS_TRIANGLE
#include "ShapeFunctionMacros.hpp"
#undef DEFINE_P2_SHAPE_FUNCTIONS_TRIANGLE

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
      real_t aux4 = ( tmp3 * DetaN##i - tmp1 * DxiN##i ) * aux3;                                                \
      real_t aux5 = ( tmp2 * DxiN##i - tmp0 * DetaN##i ) * aux3;                                                \
                                                                                                                \
      real_t aux6 = ( DPsi( 1, 1 ) * aux4 - DPsi( 1, 0 ) * aux5 ) / detDPsi;                                    \
                                                                                                                \
      elMat( i, j ) -= QUADWEIGHTS[k] * detJacPhiInv * aux6 * std::abs( detDPsi ) * SF_N##j;                    \
   }

      // compute Jacobian determinant of inverse pull-back mapping
      real_t tmp0 = coords[0][0] - coords[1][0];
      real_t tmp1 = coords[0][1] - coords[2][1];
      real_t tmp2 = coords[0][0] - coords[2][0];
      real_t tmp3 = coords[0][1] - coords[1][1];

      real_t aux1 = tmp0 * tmp1 - tmp2 * tmp3;
      real_t aux3 = 1.0 / aux1;

      real_t detJacPhiInv = std::abs( aux1 );

      // quadrature point mapped to computational triangle
      Point3D mappedPt;

      // set matrix entries
      INTEGRATE2D( 0, 0 );
      INTEGRATE2D( 0, 1 );
      INTEGRATE2D( 0, 2 );
      INTEGRATE2D( 0, 3 );
      INTEGRATE2D( 0, 4 );
      INTEGRATE2D( 0, 5 );

      INTEGRATE2D( 1, 0 );
      INTEGRATE2D( 1, 1 );
      INTEGRATE2D( 1, 2 );
      INTEGRATE2D( 1, 3 );
      INTEGRATE2D( 1, 4 );
      INTEGRATE2D( 1, 5 );

      INTEGRATE2D( 2, 0 );
      INTEGRATE2D( 2, 1 );
      INTEGRATE2D( 2, 2 );
      INTEGRATE2D( 2, 3 );
      INTEGRATE2D( 2, 4 );
      INTEGRATE2D( 2, 5 );

      INTEGRATE2D( 3, 0 );
      INTEGRATE2D( 3, 1 );
      INTEGRATE2D( 3, 2 );
      INTEGRATE2D( 3, 3 );
      INTEGRATE2D( 3, 4 );
      INTEGRATE2D( 3, 5 );

      INTEGRATE2D( 4, 0 );
      INTEGRATE2D( 4, 1 );
      INTEGRATE2D( 4, 2 );
      INTEGRATE2D( 4, 3 );
      INTEGRATE2D( 4, 4 );
      INTEGRATE2D( 4, 5 );

      INTEGRATE2D( 5, 0 );
      INTEGRATE2D( 5, 1 );
      INTEGRATE2D( 5, 2 );
      INTEGRATE2D( 5, 3 );
      INTEGRATE2D( 5, 4 );
      INTEGRATE2D( 5, 5 );

#define UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE
#include "ShapeFunctionMacros.hpp"
#undef UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE

#define UNDEFINE_P2_SHAPE_FUNCTIONS_TRIANGLE
#include "ShapeFunctionMacros.hpp"
#undef UNDEFINE_P2_SHAPE_FUNCTIONS_TRIANGLE

#undef INTEGRATE2D
   };

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix10r& elMat ) const final
   {
      WALBERLA_ABORT( "P2FormDivT not implemented for 3D, yet." );
   };
};

// --------
//  div_y
// --------
template <>
class P2Form_divt< 1 > : public P2FormHyTeG
{
 public:
   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix6r& elMat ) const final
   {

// Derivatives of P2 shape functions on unit triangle
#define DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE
#include "ShapeFunctionMacros.hpp"
#undef DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE

// P2 shape functions on unit triangle
#define DEFINE_P2_SHAPE_FUNCTIONS_TRIANGLE
#include "ShapeFunctionMacros.hpp"
#undef DEFINE_P2_SHAPE_FUNCTIONS_TRIANGLE

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
      real_t aux4 = ( tmp3 * DetaN##i - tmp1 * DxiN##i ) * aux3;                                                \
      real_t aux5 = ( tmp2 * DxiN##i - tmp0 * DetaN##i ) * aux3;                                                \
                                                                                                                \
      real_t aux6 = ( DPsi( 0, 0 ) * aux5 - DPsi( 0, 1 ) * aux4 ) / detDPsi;                                    \
                                                                                                                \
      elMat( i, j ) -= QUADWEIGHTS[k] * detJacPhiInv * aux6 * std::abs( detDPsi ) * SF_N##j;                    \
   }

      // compute Jacobian determinant of inverse pull-back mapping
      real_t tmp0 = coords[0][0] - coords[1][0];
      real_t tmp1 = coords[0][1] - coords[2][1];
      real_t tmp2 = coords[0][0] - coords[2][0];
      real_t tmp3 = coords[0][1] - coords[1][1];

      real_t aux1 = tmp0 * tmp1 - tmp2 * tmp3;
      real_t aux3 = 1.0 / aux1;

      real_t detJacPhiInv = std::abs( aux1 );

      // quadrature point mapped to computational triangle
      Point3D mappedPt;

      // set matrix entries
      INTEGRATE2D( 0, 0 );
      INTEGRATE2D( 0, 1 );
      INTEGRATE2D( 0, 2 );
      INTEGRATE2D( 0, 3 );
      INTEGRATE2D( 0, 4 );
      INTEGRATE2D( 0, 5 );

      INTEGRATE2D( 1, 0 );
      INTEGRATE2D( 1, 1 );
      INTEGRATE2D( 1, 2 );
      INTEGRATE2D( 1, 3 );
      INTEGRATE2D( 1, 4 );
      INTEGRATE2D( 1, 5 );

      INTEGRATE2D( 2, 0 );
      INTEGRATE2D( 2, 1 );
      INTEGRATE2D( 2, 2 );
      INTEGRATE2D( 2, 3 );
      INTEGRATE2D( 2, 4 );
      INTEGRATE2D( 2, 5 );

      INTEGRATE2D( 3, 0 );
      INTEGRATE2D( 3, 1 );
      INTEGRATE2D( 3, 2 );
      INTEGRATE2D( 3, 3 );
      INTEGRATE2D( 3, 4 );
      INTEGRATE2D( 3, 5 );

      INTEGRATE2D( 4, 0 );
      INTEGRATE2D( 4, 1 );
      INTEGRATE2D( 4, 2 );
      INTEGRATE2D( 4, 3 );
      INTEGRATE2D( 4, 4 );
      INTEGRATE2D( 4, 5 );

      INTEGRATE2D( 5, 0 );
      INTEGRATE2D( 5, 1 );
      INTEGRATE2D( 5, 2 );
      INTEGRATE2D( 5, 3 );
      INTEGRATE2D( 5, 4 );
      INTEGRATE2D( 5, 5 );

#define UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE
#include "ShapeFunctionMacros.hpp"
#undef UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE

#define UNDEFINE_P2_SHAPE_FUNCTIONS_TRIANGLE
#include "ShapeFunctionMacros.hpp"
#undef UNDEFINE_P2_SHAPE_FUNCTIONS_TRIANGLE

#undef INTEGRATE2D
   };

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix10r& elMat ) const final
   {
      WALBERLA_ABORT( "P2FormDivT not implemented for 3D, yet." );
   };
};

} // namespace hyteg

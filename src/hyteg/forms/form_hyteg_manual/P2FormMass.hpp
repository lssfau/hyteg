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

class P2Form_mass : public P2FormHyTeG {

public:

  void integrateAll( const std::array< Point3D, 3 >& coords, Matrix6r& elMat ) const final {

// Shape functions on unit triangle
#define SF_N0 ( L1 * ( 2.0 * L1 - 1.0 ) )
#define SF_N1 ( L2 * ( 2.0 * L2 - 1.0 ) )
#define SF_N2 ( L3 * ( 2.0 * L3 - 1.0 ) )
#define SF_N3 ( 4.0 * L2 * L3 )
#define SF_N4 ( 4.0 * L1 * L3 )
#define SF_N5 ( 4.0 * L1 * L2 )

// Debugging FPE problem with CircularMap
// #define FPE_DEBUG( cc ) WALBERLA_LOG_INFO_ON_ROOT( "Mapped point: " << cc )
#define FPE_DEBUG( cc )

// Select quadrature rule
// #define QUADPOINTS quadrature::D5_points
// #define QUADWEIGHTS quadrature::D5_weights
#define QUADPOINTS quadrature::D6_points
#define QUADWEIGHTS quadrature::D6_weights

// Executing quadrature rule
#define INTEGRATE2D(i,j)                                                                                      \
  elMat(i,j) = 0.0;                                                                                           \
  for( uint_t k = 0; k < QUADWEIGHTS.size(); k++ ) {                                                          \
    real_t L2 = QUADPOINTS[k][0];                                                                             \
    real_t L3 = QUADPOINTS[k][1];                                                                             \
    real_t L1 = 1.0 - L2 - L3;                                                                                \
    mappedPt[0] = ( coords[1][0] - coords[0][0] ) * L2 + ( coords[2][0] - coords[0][0] ) * L3 + coords[0][0]; \
    mappedPt[1] = ( coords[1][1] - coords[0][1] ) * L2 + ( coords[2][1] - coords[0][1] ) * L3 + coords[0][1]; \
    FPE_DEBUG( mappedPt ); \
    real_t detDPsi = std::abs( geometryMap_->evalDetDF( mappedPt ) );                                         \
    elMat(i,j) += QUADWEIGHTS[k] * detJacPhiInv * detDPsi * SF_N ## i * SF_N ## j;                            \
  }                                                                                                           \
  elMat(j,i) = elMat(i,j)

    // compute Jacobian determinant of inverse pull-back mapping
    real_t detJacPhiInv = ( coords[1][0] - coords[0][0] ) *  ( coords[2][1] - coords[0][1] )
      -  ( coords[2][0] - coords[0][0] ) *  ( coords[1][1] - coords[0][1] );
    detJacPhiInv = std::abs( detJacPhiInv );

    // Quadrature point mapped to computational triangle
    Point3D mappedPt;

    // ------------
    //  Zeroth row
    // ------------

    INTEGRATE2D(0,0);
    INTEGRATE2D(0,1);
    INTEGRATE2D(0,2);
    INTEGRATE2D(0,3);
    INTEGRATE2D(0,4);
    INTEGRATE2D(0,5);

    // -----------
    //  First row
    // -----------
    INTEGRATE2D(1,1);
    INTEGRATE2D(1,2);
    INTEGRATE2D(1,3);
    INTEGRATE2D(1,4);
    INTEGRATE2D(1,5);

    // ------------
    //  Second row
    // ------------
    INTEGRATE2D(2,2);
    INTEGRATE2D(2,3);
    INTEGRATE2D(2,4);
    INTEGRATE2D(2,5);

    // -----------
    //  Third row
    // -----------
    INTEGRATE2D(3,3);
    INTEGRATE2D(3,4);
    INTEGRATE2D(3,5);

    // -----------
    //  Forth row
    // -----------
    INTEGRATE2D(4,4);
    INTEGRATE2D(4,5);

    // -----------
    //  Fifth row
    // -----------
    INTEGRATE2D(5,5);

#undef SF_N0
#undef SF_N1
#undef SF_N2
#undef SF_N3
#undef SF_N4
#undef SF_N5
#undef INTEGRATE2D
  };

  void integrateAll( const std::array< Point3D, 4 >& coords, Matrix10r& elMat ) const final {

// Select quadrature rule
#define CUBAPOINTS cubature::T4_points
#define CUBAWEIGHTS cubature::T4_weights

// Shape functions on unit tetrahedron
#define SF_N0 ( L1 * ( 2.0 * L1 - 1.0 ) )
#define SF_N1 ( L2 * ( 2.0 * L2 - 1.0 ) )
#define SF_N2 ( L3 * ( 2.0 * L3 - 1.0 ) )
#define SF_N3 ( L4 * ( 2.0 * L4 - 1.0 ) )
#define SF_N4 ( 4.0 * L3 * L4 )
#define SF_N5 ( 4.0 * L2 * L4 )
#define SF_N6 ( 4.0 * L2 * L3 )
#define SF_N7 ( 4.0 * L1 * L4 )
#define SF_N8 ( 4.0 * L1 * L3 )
#define SF_N9 ( 4.0 * L1 * L2 )

// Executing quadrature rule
#define INTEGRATE3D(i,j)                                                                                                                             \
  elMat(i,j) = 0.0;                                                                                                                                  \
  for( uint_t k = 0; k < CUBAWEIGHTS.size(); k++ ) {                                                                                                 \
    real_t L2 = CUBAPOINTS[k][0];                                                                                                                    \
    real_t L3 = CUBAPOINTS[k][1];                                                                                                                    \
    real_t L4 = CUBAPOINTS[k][2];                                                                                                                    \
    real_t L1 = 1.0 - L2 - L3 - L4;                                                                                                                  \
    mappedPt[0] = ( coords[1][0] - coords[0][0] ) * L2 + ( coords[2][0] - coords[0][0] ) * L3 + ( coords[3][0] - coords[0][0] ) * L4 + coords[0][0]; \
    mappedPt[1] = ( coords[1][1] - coords[0][1] ) * L2 + ( coords[2][1] - coords[0][1] ) * L3 + ( coords[3][1] - coords[0][1] ) * L4 + coords[0][1]; \
    mappedPt[2] = ( coords[1][2] - coords[0][2] ) * L2 + ( coords[2][2] - coords[0][2] ) * L3 + ( coords[3][2] - coords[0][2] ) * L4 + coords[0][2]; \
    real_t detDPsi = std::abs( geometryMap_->evalDetDF( mappedPt ) );                                                                                \
    elMat(i,j) += CUBAWEIGHTS[k] * detJacPhiInv * detDPsi * SF_N ## i * SF_N ## j;                                                                   \
  }                                                                                                                                                  \
  elMat(j,i) = elMat(i,j)

    // compute Jacobian determinant of inverse pull-back mapping
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
    detJacPhiInv = std::abs( detJacPhiInv );

    // Cubature point mapped to computational tetrahedron
    Point3D mappedPt;

    INTEGRATE3D(0,0);
    INTEGRATE3D(0,1);
    INTEGRATE3D(0,2);
    INTEGRATE3D(0,3);
    INTEGRATE3D(0,4);
    INTEGRATE3D(0,5);
    INTEGRATE3D(0,6);
    INTEGRATE3D(0,7);
    INTEGRATE3D(0,8);
    INTEGRATE3D(0,9);

    INTEGRATE3D(1,1);
    INTEGRATE3D(1,2);
    INTEGRATE3D(1,3);
    INTEGRATE3D(1,4);
    INTEGRATE3D(1,5);
    INTEGRATE3D(1,6);
    INTEGRATE3D(1,7);
    INTEGRATE3D(1,8);
    INTEGRATE3D(1,9);

    INTEGRATE3D(2,2);
    INTEGRATE3D(2,3);
    INTEGRATE3D(2,4);
    INTEGRATE3D(2,5);
    INTEGRATE3D(2,6);
    INTEGRATE3D(2,7);
    INTEGRATE3D(2,8);
    INTEGRATE3D(2,9);

    INTEGRATE3D(3,3);
    INTEGRATE3D(3,4);
    INTEGRATE3D(3,5);
    INTEGRATE3D(3,6);
    INTEGRATE3D(3,7);
    INTEGRATE3D(3,8);
    INTEGRATE3D(3,9);

    INTEGRATE3D(4,4);
    INTEGRATE3D(4,5);
    INTEGRATE3D(4,6);
    INTEGRATE3D(4,7);
    INTEGRATE3D(4,8);
    INTEGRATE3D(4,9);

    INTEGRATE3D(5,5);
    INTEGRATE3D(5,6);
    INTEGRATE3D(5,7);
    INTEGRATE3D(5,8);
    INTEGRATE3D(5,9);

    INTEGRATE3D(6,6);
    INTEGRATE3D(6,7);
    INTEGRATE3D(6,8);
    INTEGRATE3D(6,9);

    INTEGRATE3D(7,7);
    INTEGRATE3D(7,8);
    INTEGRATE3D(7,9);

    INTEGRATE3D(8,8);
    INTEGRATE3D(8,9);

    INTEGRATE3D(9,9);

#undef SF_N0
#undef SF_N1
#undef SF_N2
#undef SF_N3
#undef SF_N4
#undef SF_N5
#undef SF_N6
#undef SF_N7
#undef SF_N8
#undef SF_N9
#undef INTEGRATE3D
#undef CUBAPOINTS
#undef CUBAWEIGHTS
  };
};
}

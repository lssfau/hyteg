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

// Shape functions on unit triangle
#define DEFINE_P2_SHAPE_FUNCTIONS_TRIANGLE
#include "../ShapeFunctionMacros.hpp"
#undef DEFINE_P2_SHAPE_FUNCTIONS_TRIANGLE

// Executing quadrature rule
#define INTEGRATE2D(i,j)                                                                                      \
  elMat(i,j) = 0.0;                                                                                           \
  for( uint_t k = 0; k < QUADWEIGHTS.size(); k++ ) {                                                          \
    real_t L2 = QUADPOINTS[k][0];                                                                             \
    real_t L3 = QUADPOINTS[k][1];                                                                             \
    real_t L1 = 1.0 - L2 - L3;                                                                                \
    mappedPt[0] = ( coords[1][0] - coords[0][0] ) * L2 + ( coords[2][0] - coords[0][0] ) * L3 + coords[0][0]; \
    mappedPt[1] = ( coords[1][1] - coords[0][1] ) * L2 + ( coords[2][1] - coords[0][1] ) * L3 + coords[0][1]; \
                                                                                                              \
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

#define UNDEFINE_P2_SHAPE_FUNCTIONS_TRIANGLE
#include "../ShapeFunctionMacros.hpp"
#undef UNDEFINE_P2_SHAPE_FUNCTIONS_TRIANGLE

#undef INTEGRATE2D

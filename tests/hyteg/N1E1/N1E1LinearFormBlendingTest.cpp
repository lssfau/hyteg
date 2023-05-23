/*
 * Copyright (c) 2023 Daniel Bauer.
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

// The correct values for this test have been obtained by integrating over the
// elements symbolically and without mapping to a reference element.
// Instead, we form the FEM-basis for each affine element and also figure out
// the appropriate integration bounds.
// This way, even if we messed up e.g., the transformation to the reference
// element on paper, we will notice anyway.
//
// We test on the reference tet on refinement level 1.
// We use the reference tet so that we can figure out the integration bounds for
// all micro-cells.
// Since the colored micro-cells on level 1 are rotated, scaled and translated,
// this test is exhaustive.
//
// The integration is implemented in curl-curl-test.py.

#include <memory>

#include "core/debug/TestSubsystem.h"
#include "core/math/Random.h"
#include "core/mpi/Environment.h"

#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_linear_form_affine_q6.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_linear_form_blending_q6.hpp"
#include "hyteg/geometry/AffineMap3D.hpp"

using walberla::real_t;
using namespace hyteg;

void test()
{
   walberla::math::seedRandomGenerator( 42 );
   const uint_t numRandomEvaluations = 100;

   Matrix3r B;
   // clang-format off
   B <<  0.1, 0.2, 0.3,
        -0.8, 0.5, 0.0,
         1.0, 1.0, 0.5;
   // clang-format on
   const Point3D b{ 2.0, 3.0, 4.0 };
   auto          affineMap = std::make_shared< AffineMap3D >( B, b );

   const Point3D                                    a0   = { 1, 2, 3 };
   const Point3D                                    a1   = { 4, 5, 6 };
   const std::function< Point3D( const Point3D& ) > func = [&]( const Point3D& x ) { return Point3D{ a0 + a1.cross( x ) }; };

   forms::n1e1_linear_form_affine_q6   formNoBlending{ func };
   forms::n1e1_linear_form_blending_q6 formBlending{ func };
   formBlending.setGeometryMap( affineMap );

   for ( uint_t r = 0; r < numRandomEvaluations; ++r )
   {
      std::array< Point3D, 4 > coords;
      for ( size_t i = 0; i < 4; ++i )
      {
         for ( int j = 0; j < 3; ++j )
         {
            coords[i][j] = real_c( walberla::math::realRandom( 0.0, 1.0 ) );
         }
      }

      std::array< Point3D, 4 > coordsBlending;
      for ( size_t i = 0; i < 4; ++i )
      {
         affineMap->evalFinv( coords[i], coordsBlending[i] );
      }

      Matrix6r elMatNoBlending;
      Matrix6r elMatBlending;

      formNoBlending.integrateAll( coords, elMatNoBlending );
      formBlending.integrateAll( coordsBlending, elMatBlending );

      for ( int i = 0; i < 6; ++i )
      {
         for ( int j = 0; j < 6; ++j )
         {
            WALBERLA_CHECK_FLOAT_EQUAL( elMatNoBlending( i, j ), elMatBlending( i, j ) )
         }
      }
   }
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   test();
   return EXIT_SUCCESS;
}

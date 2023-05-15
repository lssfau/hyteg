/*
 * Copyright (c) 2022-2023 Daniel Bauer.
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
#include "core/mpi/Environment.h"

#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_curl_curl_affine_qe.hpp"
#include "hyteg/forms/form_hyteg_generated/n1e1/n1e1_curl_curl_blending_q2.hpp"
#include "hyteg/geometry/AffineMap3D.hpp"
#include "hyteg/geometry/IdentityMap.hpp"

using walberla::real_t;
using namespace hyteg;

template < typename Form >
void test( Form form, const std::array< Point3D, 4 >& coords, const Matrix6r& correct )
{
   Matrix6r elMat;

   form.integrateAll( coords, elMat );

   for ( int i = 0; i < 6; ++i )
   {
      for ( int j = 0; j < 6; ++j )
      {
         WALBERLA_CHECK_FLOAT_EQUAL( elMat( i, j ), correct( i, j ) )
      }
   }
}

void test( const std::array< Point3D, 4 >& coords, const Matrix6r& correct )
{
   auto testForms = [&]( auto&&... forms ) { ( test( forms, coords, correct ), ... ); };

   forms::n1e1_curl_curl_blending_q2 formBlendingId;
   formBlendingId.setGeometryMap( std::make_shared< IdentityMap >() );

   testForms( forms::n1e1_curl_curl_affine_qe{}, formBlendingId );

   // blending

   Matrix3r B;
   // clang-format off
   B <<  0.1, 0.2, 0.3,
        -0.8, 0.5, 0.0,
         1.0, 1.0, 0.5;
   // clang-format on
   const Point3D b{ 2.0, 3.0, 4.0 };
   auto          affineMap = std::make_shared< AffineMap3D >( B, b );

   std::array< Point3D, 4 > coordsBlending;
   for ( size_t i = 0; i < 4; ++i )
   {
      affineMap->evalFinv( coords[i], coordsBlending[i] );
   }

   forms::n1e1_curl_curl_blending_q2 formBlendingAf;
   formBlendingAf.setGeometryMap( affineMap );

   test( formBlendingAf, coordsBlending, correct );
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   const Point3D v0{ { 0, 0, 0 } };
   const Point3D v1{ { 0.5, 0, 0 } };
   const Point3D v2{ { 1, 0, 0 } };
   const Point3D v3{ { 0, 0.5, 0 } };
   const Point3D v4{ { 0.5, 0.5, 0 } };
   const Point3D v5{ { 0, 1, 0 } };
   const Point3D v6{ { 0, 0, 0.5 } };
   const Point3D v7{ { 0.5, 0, 0.5 } };
   const Point3D v8{ { 0, 0.5, 0.5 } };
   const Point3D v9{ { 0, 0, 1 } };

   Matrix6r wu, gu, bu, bd, gd;

   // clang-format off
   wu <<  1,  0,  0, -1,  1,  0,
          0,  1,  0, -1,  0,  1,
          0,  0,  1,  0, -1,  1,
         -1, -1,  0,  2, -1, -1,
          1,  0, -1, -1,  2, -1,
          0,  1,  1, -1, -1,  2;

   gu <<  2, -1,  0, -1,  2, -1,
         -1,  2, -1, -1,  0,  1,
          0, -1,  1,  1, -1,  0,
         -1, -1,  1,  2, -2,  0,
          2,  0, -1, -2,  3, -1,
         -1,  1,  0,  0, -1,  1;

   bu <<  2, -1,  0, -1,  2, -1,
         -1,  1,  0,  0, -1,  1,
          0,  0,  1,  0, -1,  1,
         -1,  0,  0,  1, -1,  0,
          2, -1, -1, -1,  3, -2,
         -1,  1,  1,  0, -2,  2;

   bd <<  2, -2,  1,  0,  1, -1,
         -2,  3, -1, -1, -1,  2,
          1, -1,  1,  0,  0,  0,
          0, -1,  0,  1,  0, -1,
          1, -1,  0,  0,  1, -1,
         -1,  2,  0, -1, -1,  2;

   gd <<  1, -1,  0,  0,  1, -1,
         -1,  3, -1, -2,  0,  2,
          0, -1,  1,  1, -1,  0,
          0, -2,  1,  2, -1, -1,
          1,  0, -1, -1,  2, -1,
         -1,  2,  0, -1, -1,  2;
   // clang-format on

   test( { v0, v1, v3, v6 }, 4.0 / 3.0 * wu );
   test( { v1, v2, v4, v7 }, 4.0 / 3.0 * wu );
   test( { v3, v4, v5, v8 }, 4.0 / 3.0 * wu );
   test( { v6, v7, v8, v9 }, 4.0 / 3.0 * wu );
   test( { v1, v3, v6, v7 }, 4.0 / 3.0 * gu );
   test( { v1, v3, v4, v7 }, 4.0 / 3.0 * bu );
   test( { v3, v6, v7, v8 }, 4.0 / 3.0 * bd );
   test( { v3, v4, v7, v8 }, 4.0 / 3.0 * gd );

   return EXIT_SUCCESS;
}

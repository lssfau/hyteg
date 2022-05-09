/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "core/DataTypes.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1dgefunctionspace/P1DGEOperators.hpp"
#include "hyteg/petsc/PETScManager.hpp"

using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

void testDivForm()
{
   hyteg::dg::DGDivFormP1EDG form;

   // we test two triangles with coordinates (p0, p1, p2) and (p0, p2, p3):
   Eigen::Matrix< real_t, 3, 1 > p0{ 0, 0, 0 };
   Eigen::Matrix< real_t, 3, 1 > p1{ 1, 0, 0 };
   Eigen::Matrix< real_t, 3, 1 > p2{ -0.2, 1, 0 };
   Eigen::Matrix< real_t, 3, 1 > p3{ -1, 0.5, 0 };
   Eigen::Matrix< real_t, 3, 1 > n1{ -1, -0.2, 0 };

   // offset to have some shift from the origin in the coordinates
   Eigen::Matrix< real_t, 3, 1 > offset{ 0.05, -0.7, 0 };
   p0 += offset;
   p1 += offset;
   p2 += offset;
   p3 += offset;

   hyteg::dg::DGBasisLinearLagrange_Example basis;

   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > elMat;
   elMat.resize( 3, 1 );

   // checking inner facet integral
   {
      // expected values precalculated in sympy notebook:
      std::vector< real_t > expected = { 0.0849836585598798, 0, 0.0849836585598798 };

      form.integrateFacetInner( 2, { p0, p1, p2 }, { p0, p2 }, p1, n1, basis, basis, 1, 0, elMat );

      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 1, 0 ), expected[1] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 2, 0 ), expected[2] );

      // check a different permutation
      form.integrateFacetInner( 2, { p0, p1, p2 }, { p2, p0 }, p1, n1, basis, basis, 1, 0, elMat );

      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 1, 0 ), expected[1] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 2, 0 ), expected[2] );

      form.integrateFacetInner( 2, { p2, p0, p1 }, { p2, p0 }, p1, n1, basis, basis, 1, 0, elMat );

      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[2] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 1, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 2, 0 ), expected[1] );
   }

   // checking outer facet integral
   {
      // expected values precalculated in sympy notebook:
      std::vector< real_t > expected = { 0.0764852927038918, 0, 0.0764852927038918 };

      form.integrateFacetCoupling( 2, { p0, p1, p2 }, { p0, p2, p3 }, { p0, p2 }, p1, p3, n1, basis, basis, 1, 0, elMat );

      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 1, 0 ), expected[1] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 2, 0 ), expected[2] );

      // check a different permutation
      form.integrateFacetCoupling( 2, { p2, p0, p1 }, { p3, p2, p0 }, { p2, p0 }, p1, p3, n1, basis, basis, 1, 0, elMat );

      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[2] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 1, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 2, 0 ), expected[1] );
   }

   // checking dirichlet integral
   {
      // expected values precalculated in sympy notebook:
      std::vector< real_t > expected = { 2. * 0.0849836585598798, 0, 2. * 0.0849836585598798 };

      form.integrateFacetDirichletBoundary( 2, { p0, p1, p2 }, { p0, p2 }, p1, n1, basis, basis, 1, 0, elMat );

      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 1, 0 ), expected[1] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 2, 0 ), expected[2] );

      // check a different permutation
      form.integrateFacetDirichletBoundary( 2, { p0, p1, p2 }, { p2, p0 }, p1, n1, basis, basis, 1, 0, elMat );

      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 1, 0 ), expected[1] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 2, 0 ), expected[2] );

      form.integrateFacetDirichletBoundary( 2, { p2, p0, p1 }, { p2, p0 }, p1, n1, basis, basis, 1, 0, elMat );

      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[2] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 1, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 2, 0 ), expected[1] );
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager( &argc, &argv );

   hyteg::testDivForm();

   return 0;
}

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
#include "hyteg/p1dgefunctionspace/EDGEpsilonForm.hpp"
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

   // checking inner and outer facet integrals:
   {
      hyteg::dg::DGDivtFormEDGP0 form;

      elMat.resize( 1, 1 );

      // expected values precalculated in sympy notebook:
      std::vector< real_t > expected = { 0.16996731711976 };

      // check inner integral
      form.integrateFacetInner( 2, { p0, p1, p2 }, { p0, p2 }, p1, n1, basis, basis, 0, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );

      // checkout outer integral
      form.integrateFacetCoupling( 2, { p0, p1, p2 }, { p0, p2, p3 }, { p0, p2 }, p1, p3, n1, basis, basis, 0, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );

      // check a different permutation
      form.integrateFacetCoupling( 2, { p0, p1, p2 }, { p0, p2, p3 }, { p2, p0 }, p1, p3, n1, basis, basis, 0, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );

      form.integrateFacetCoupling( 2, { p0, p1, p2 }, { p3, p0, p2 }, { p2, p0 }, p1, p3, n1, basis, basis, 0, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
   }

   // checking inner and outer facet integrals:
   {
      hyteg::dg::DGDivtFormP1P0_0 form;

      elMat.resize( 3, 1 );

      // expected values precalculated in sympy notebook:
      std::vector< real_t > expected = { -0.254950975679639, 0, -0.254950975679639 };

      // check inner integral
      form.integrateFacetInner( 2, { p0, p1, p2 }, { p0, p2 }, p1, n1, basis, basis, 0, 1, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 1, 0 ), expected[1] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 2, 0 ), expected[2] );

      // check a different permutation
      form.integrateFacetCoupling( 2, { p0, p1, p2 }, { p0, p2, p3 }, { p2, p0 }, p1, p3, n1, basis, basis, 0, 1, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 1, 0 ), expected[1] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 2, 0 ), expected[2] );

      form.integrateFacetCoupling( 2, { p2, p0, p1 }, { p0, p2, p3 }, { p2, p0 }, p1, p3, n1, basis, basis, 0, 1, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[2] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 1, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 2, 0 ), expected[1] );
   }

   {
      hyteg::dg::DGDivFormP0EDG form;

      elMat.resize( 1, 1 );

      // expected values precalculated in sympy notebook:
      std::vector< real_t > expected = { 0.152970585407784 };

      // checkout outer integral
      form.integrateFacetCoupling( 2, { p0, p1, p2 }, { p0, p2, p3 }, { p0, p2 }, p1, p3, n1, basis, basis, 0, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );

      // check a different permutation
      form.integrateFacetCoupling( 2, { p0, p1, p2 }, { p0, p2, p3 }, { p2, p0 }, p1, p3, n1, basis, basis, 0, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );

      form.integrateFacetCoupling( 2, { p0, p1, p2 }, { p3, p0, p2 }, { p2, p0 }, p1, p3, n1, basis, basis, 0, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
   }

   {
      hyteg::dg::DGDivFormP0EDG form;

      elMat.resize( 1, 1 );

      // expected values precalculated in sympy notebook:
      std::vector< real_t > expected = { 0.16996731711976 };

      // checkout outer integral
      form.integrateFacetInner( 2, { p0, p1, p2 }, { p0, p2 }, p1, n1, basis, basis, 0, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );

      // check a different permutation
      form.integrateFacetInner( 2, { p0, p1, p2 }, { p2, p0 }, p1, n1, basis, basis, 0, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );

      form.integrateFacetInner( 2, { p2, p0, p1 }, { p2, p0 }, p1, n1, basis, basis, 0, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
   }

   {
      hyteg::dg::DGDivFormP0P1_0 form;

      elMat.resize( 1, 3 );

      // expected values precalculated in sympy notebook:
      std::vector< real_t > expected = { -0.254950975679639, 0, -0.254950975679639 };

      // check inner integral
      form.integrateFacetInner( 2, { p0, p1, p2 }, { p0, p2 }, p1, n1, basis, basis, 1, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 1 ), expected[1] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 2 ), expected[2] );

      // check a different permutation
      form.integrateFacetInner( 2, { p0, p1, p2 }, { p2, p0 }, p1, n1, basis, basis, 1, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 1 ), expected[1] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 2 ), expected[2] );

      form.integrateFacetInner( 2, { p2, p0, p1 }, { p2, p0 }, p1, n1, basis, basis, 1, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[2] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 1 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 2 ), expected[1] );
   }

   {
      hyteg::dg::DGDivFormP0P1_0 form;

      elMat.resize( 1, 3 );

      // expected values precalculated in sympy notebook:
      // we have separated the consistency and symmetry terms from the penalty to allow different penalty values sigma:
      std::vector< real_t > expected = { 0.254950975679639, 0.254950975679639, 0 };

      // check inner integral
      form.integrateFacetCoupling( 2, { p0, p1, p2 }, { p0, p2, p3 }, { p0, p2 }, p1, p3, n1, basis, basis, 1, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 1 ), expected[1] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 2 ), expected[2] );

      // check a different permutation
      form.integrateFacetCoupling( 2, { p0, p1, p2 }, { p0, p2, p3 }, { p2, p0 }, p1, p3, n1, basis, basis, 1, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 1 ), expected[1] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 2 ), expected[2] );

      form.integrateFacetCoupling( 2, { p0, p1, p2 }, { p3, p0, p2 }, { p2, p0 }, p1, p3, n1, basis, basis, 1, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[2] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 1 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 2 ), expected[1] );
   }
}

void testLaplaceForm()
{
   // we test two triangles with coordinates (p0, p1, p2) and (p0, p2, p3):
   Eigen::Matrix< real_t, 3, 1 > p0{ 0, 0, 0 };
   Eigen::Matrix< real_t, 3, 1 > p1{ 1, 0, 0 };
   Eigen::Matrix< real_t, 3, 1 > p2{ -0.2, 1, 0 };
   Eigen::Matrix< real_t, 3, 1 > p3{ -1, 0.5, 0 };
   Eigen::Matrix< real_t, 3, 1 > n1{ -1, -0.2, 0 };

   // offset to have some shift from the origin in the coordinates
   // Eigen::Matrix< real_t, 3, 1 > offset{ 0.05, -0.7, 0 };
   Eigen::Matrix< real_t, 3, 1 > offset{ 0., 0., 0. };
   p0 += offset;
   p1 += offset;
   p2 += offset;
   p3 += offset;

   hyteg::dg::DGBasisLinearLagrange_Example basis;

   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > elMat;
   elMat.resize( 3, 1 );

   const real_t sigma = 6.;

   {
      hyteg::dg::DGVectorLaplaceFormP1P1_00 form;

      Eigen::Matrix< real_t, 3, 3 > expected_consistency_and_symmetry;
      // clang-format off
      expected_consistency_and_symmetry <<
          0.632278419685505, -0.265149014706825, 0.265149014706825,
          -0.265149014706825, 0, -0.265149014706825,
          0.265149014706825, -0.265149014706825, -0.101980390271856;
      // clang-format on

      Eigen::Matrix< real_t, 3, 3 > expected_penalty;
      // clang-format off
      expected_penalty <<
          0.33333333333333337, 0, 0.16666666666666674,
          0, 0, 0,
          0.16666666666666674, 0, 0.3333333333333333;
      // clang-format on

      Eigen::Matrix< real_t, 3, 3 > expected = -expected_consistency_and_symmetry + sigma * expected_penalty;

      elMat.resize( 3, 3 );
      form.integrateFacetInner( 2, { p0, p1, p2 }, { p0, p2 }, p1, n1, basis, basis, 1, 1, elMat );

      for ( uint_t i = 0; i < 3; i += 1 )
         for ( uint_t j = 0; j < 3; j += 1 )
            WALBERLA_CHECK_FLOAT_EQUAL( expected( i, j ), elMat( i, j ) )
   }

   // check inner facet integral EDG-EDG
   {
      hyteg::dg::DGVectorLaplaceFormEDGEDG form;

      elMat.resize( 1, 1 );
      form.integrateFacetInner( 2, { p0, p1, p2 }, { p0, p2 }, p1, n1, basis, basis, 0, 0, elMat );

      real_t expected_consistency_and_symmetry = 0.339934634239519;
      real_t expected_penalty                  = 0.24888888888888888;
      real_t expected                          = -expected_consistency_and_symmetry + sigma * expected_penalty;

      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected );
   }

   // check outer facet integral EDG-EDG
   {
      hyteg::dg::DGVectorLaplaceFormEDGEDG form;

      elMat.resize( 1, 1 );
      form.integrateFacetCoupling( 2, { p0, p1, p2 }, { p0, p2, p3 }, { p0, p2 }, p1, p3, n1, basis, basis, 0, 0, elMat );

      real_t expected_consistency_and_symmetry = 0.322937902527543;
      real_t expected_penalty                  = 0.02333333333333332;
      real_t expected                          = -expected_consistency_and_symmetry + sigma * expected_penalty;

      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected );
   }

   // checking inner facet integral
   {
      hyteg::dg::DGVectorLaplaceFormEDGP1_0 form;

      elMat.resize( 3, 1 );

      // expected values precalculated in sympy notebook:
      // we have separated the consistency and symmetry terms from the penalty to allow different penalty values sigma:
      std::vector< real_t > expected_consistency_and_symmetry = { -0.486786396230991, 0.194442610785005, -0.217558165913292 };

      std::vector< real_t > expected_penalty = { -0.166666666666667, 0, -0.2 };

      std::vector< real_t > expected = {
          -expected_consistency_and_symmetry[0] + sigma * expected_penalty[0],
          -expected_consistency_and_symmetry[1] + sigma * expected_penalty[1],
          -expected_consistency_and_symmetry[2] + sigma * expected_penalty[2],
      };

      form.integrateFacetInner( 2, { p0, p1, p2 }, { p0, p2 }, p1, n1, basis, basis, 1, 0, elMat );

      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 1 ), expected[1] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 2 ), expected[2] );

      // check a different permutation
      form.integrateFacetInner( 2, { p0, p1, p2 }, { p2, p0 }, p1, n1, basis, basis, 1, 0, elMat );

      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 1 ), expected[1] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 2 ), expected[2] );

      form.integrateFacetInner( 2, { p2, p0, p1 }, { p2, p0 }, p1, n1, basis, basis, 1, 0, elMat );

      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[2] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 1 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 2 ), expected[1] );
   }

   // checking outer facet integral
   {
      hyteg::dg::DGVectorLaplaceFormEDGP1_0 form;

      elMat.resize( 3, 1 );

      // expected values precalculated in sympy notebook:
      // we have separated the consistency and symmetry terms from the penalty to allow different penalty values sigma:
      std::vector< real_t > expected_consistency_and_symmetry = { +0.325581838571628, +0.400367458104322, -0.216047345316672 };

      std::vector< real_t > expected_penalty = { 0.16666666666666666, 0.2, 0.0 };

      std::vector< real_t > expected = {
          -expected_consistency_and_symmetry[0] + sigma * expected_penalty[0],
          -expected_consistency_and_symmetry[1] + sigma * expected_penalty[1],
          -expected_consistency_and_symmetry[2] + sigma * expected_penalty[2],
      };

      form.integrateFacetCoupling( 2, { p0, p1, p2 }, { p0, p2, p3 }, { p0, p2 }, p1, p3, n1, basis, basis, 1, 0, elMat );

      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 1 ), expected[1] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 2 ), expected[2] );

      // check a different permutation
      form.integrateFacetCoupling( 2, { p0, p1, p2 }, { p0, p2, p3 }, { p2, p0 }, p1, p3, n1, basis, basis, 1, 0, elMat );

      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 1 ), expected[1] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 2 ), expected[2] );

      form.integrateFacetCoupling( 2, { p0, p1, p2 }, { p3, p0, p2 }, { p2, p0 }, p1, p3, n1, basis, basis, 1, 0, elMat );

      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[2] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 1 ), expected[0] );
      WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 2 ), expected[1] );
   }
}

void testEpsilonForm()
{
   // we test two triangles with coordinates (p0, p1, p2) and (p0, p2, p3):
   Eigen::Matrix< real_t, 3, 1 > p0{ 0, 0, 0 };
   Eigen::Matrix< real_t, 3, 1 > p1{ 1, 0, 0 };
   Eigen::Matrix< real_t, 3, 1 > p2{ -0.2, 1, 0 };
   Eigen::Matrix< real_t, 3, 1 > p3{ -1, 0.5, 0 };
   Eigen::Matrix< real_t, 3, 1 > n1{ -1, -0.2, 0 };

   // offset to have some shift from the origin in the coordinates
   // Eigen::Matrix< real_t, 3, 1 > offset{ 0.05, -0.7, 0 };
   Eigen::Matrix< real_t, 3, 1 > offset{ 0., 0., 0. };
   p0 += offset;
   p1 += offset;
   p2 += offset;
   p3 += offset;

   hyteg::dg::DGBasisLinearLagrange_Example basis;

   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > elMat;
   elMat.resize( 3, 1 );

   const real_t sigma = 6.;

   // volume integrals: (epsilon, epsilon) term
   // edg, edg pairing
   {
      hyteg::dg::EDGEpsilonFormEDGEDG form;

      // (2 * \int_T \epsilon(\phi_{edg}) : \epsilon(\phi_{edg}) dx from sympy notebook
      real_t expected_t1 = 2 * 1.0; // 2*0.96080;
      real_t expected_t2 = 2 * 0.9; // 2*1.88888888888889;

      elMat.resize( 1, 1 );

      // triangle 1
      form.integrateVolume2D( { p0, p1, p2 }, basis, basis, 0, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( expected_t1, elMat( 0, 0 ) )

      // triangle 2
      form.integrateVolume2D( { p0, p2, p3 }, basis, basis, 0, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( expected_t2, elMat( 0, 0 ) )
   }

   // p1, edg pairing
   {
      hyteg::dg::EDGEpsilonFormP1EDG_0 form;

      elMat.resize( 3, 1 );

      // triangle 1
      Eigen::Matrix< real_t, 3, 1 > expected_t1;
      // clang-format off
      expected_t1 <<
      -2*0.50, 2*0.50, 0;// 2: prefactor of volume integral term
      // clang-format on

      form.integrateVolume2D( { p0, p1, p2 }, basis, basis, 0, 1, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( expected_t1( 0, 0 ), elMat( 0, 0 ) )
      WALBERLA_CHECK_FLOAT_EQUAL( expected_t1( 1, 0 ), elMat( 1, 0 ) )
      WALBERLA_CHECK_FLOAT_EQUAL( expected_t1( 2, 0 ), elMat( 2, 0 ) )

      
      // triangle 2
      Eigen::Matrix< real_t, 3, 1 > expected_t2;
      // clang-format off
      expected_t2 <<
          2*0.250, 2*0.250, -2*0.50;
      // clang-format on
  

      form.integrateVolume2D( { p0, p2, p3 }, basis, basis, 0, 1, elMat );
      for ( uint_t i = 0; i < 3; i += 1 )
         WALBERLA_CHECK_FLOAT_EQUAL( expected_t2( i, 0 ), elMat( i, 0 ) )
   }

   // edg, p1 pairing
   {
      hyteg::dg::EDGEpsilonFormEDGP1_0 form;

      elMat.resize( 3, 1 );

      // triangle 1
      Eigen::Matrix< real_t, 3, 1 > expected_t1;
      // clang-format off
      expected_t1 <<
      -2*0.50, 2*0.50, 0;// 2: prefactor of volume integral term
      // clang-format on

      form.integrateVolume2D( { p0, p1, p2 }, basis, basis, 1, 0, elMat );
      WALBERLA_CHECK_FLOAT_EQUAL( expected_t1( 0, 0 ), elMat( 0, 0 ) )
      WALBERLA_CHECK_FLOAT_EQUAL( expected_t1( 1, 0 ), elMat( 1, 0 ) )
      WALBERLA_CHECK_FLOAT_EQUAL( expected_t1( 2, 0 ), elMat( 2, 0 ) )

      
      // triangle 2
      Eigen::Matrix< real_t, 3, 1 > expected_t2;
      // clang-format off
      expected_t2 <<
          2*0.250, 2*0.250, -2*0.50;
      // clang-format on
  

      form.integrateVolume2D( { p0, p2, p3 }, basis, basis, 1, 0, elMat );
      for ( uint_t i = 0; i < 3; i += 1 )
         WALBERLA_CHECK_FLOAT_EQUAL( expected_t2( i, 0 ), elMat( i, 0 ) )
   }
   
   // edge integrals: consistency and symmetry term
   {
      // edg, edg pairing
      {
         hyteg::dg::EDGEpsilonFormEDGEDG form;
         {
            elMat.resize( 1, 1 );

            // inner  side of facet integral
            {
               real_t expected_consistency_symmetry = 0.339934634239519; // consistency, symmetry (penalty does not contain eps)
               real_t expected_penalty              = 0.24888888888888888;
               real_t expected = -expected_consistency_symmetry + sigma*expected_penalty;
               form.integrateFacetInner(2, { p0, p1, p2 }, { p0, p2 }, p1, n1, basis, basis, 0, 0, elMat );
               WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected );
            }

            // outer/coupling side of facet integral
            {
               real_t expected_consistency_symmetry = 0.322937902527543;
               real_t expected_penalty              = 0.02333333333333332;
               real_t expected = -expected_consistency_symmetry + sigma*expected_penalty;
               form.integrateFacetCoupling(2,
                    { p0, p1, p2 }, { p0, p2, p3 }, { p0, p2 }, p1, p3, n1, basis, basis, 0, 0, elMat );

               WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected );
            }
         }
      }

      // edg, p1 pairing
      {
         hyteg::dg::EDGEpsilonFormEDGP1_0 form;

         // inner  side of facet integral
         {
            elMat.resize( 3, 1 );

            std::vector< real_t > expected_consistency_symmetry = {
                -0.413360515235255, 0.182204963952382, -0.278746400076406};
            std::vector< real_t > expected_penalty = { -0.16666666666666674, 0.0, -0.2 };
            std::vector< real_t > expected         = {
                        -expected_consistency_symmetry[0] + sigma * expected_penalty[0],
                        -expected_consistency_symmetry[1] + sigma * expected_penalty[1],
                        -expected_consistency_symmetry[2] + sigma * expected_penalty[2],
            };

            form.integrateFacetInner2D(  { p0, p1, p2 },  { p0, p2 }, p1,  n1, basis, basis, 1, 0, elMat );

            WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
            WALBERLA_CHECK_FLOAT_EQUAL( elMat( 1, 0 ), expected[1] );
            WALBERLA_CHECK_FLOAT_EQUAL( elMat( 2, 0 ), expected[2] );
         }

         // outer/coupling side of facet integral
         {
            std::vector< real_t > expected_consistency_symmetry = {
             0.379971380049951, 0.332380531256419, -0.202449959947091};
            std::vector< real_t > expected_penalty = { 0.16666666666666666, 0.2, 0.0 };
            std::vector< real_t > expected         = {
                        -expected_consistency_symmetry[0] + sigma * expected_penalty[0],
                        -expected_consistency_symmetry[1] + sigma * expected_penalty[1],
                        -expected_consistency_symmetry[2] + sigma * expected_penalty[2],
            };

            form.integrateFacetCoupling( 2, { p0, p1, p2 }, { p0, p2, p3 }, { p2, p0 }, p1, p3, n1, basis, basis, 1, 0, elMat );
            WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
            WALBERLA_CHECK_FLOAT_EQUAL( elMat( 1, 0 ), expected[1] );
            WALBERLA_CHECK_FLOAT_EQUAL( elMat( 2, 0 ), expected[2] );
         }
      }

      
      // p1, edg pairing
      {
         hyteg::dg::EDGEpsilonFormP1EDG_0 form;

         // inner  side of facet integral
         {
            elMat.resize( 3, 1 );

            std::vector< real_t > expected_consistency_symmetry = {
                -0.413360515235255, 0.182204963952382, -0.278746400076406};
            std::vector< real_t > expected_penalty = { -0.16666666666666674, 0.0, -0.2 };
            std::vector< real_t > expected         = {
                        -expected_consistency_symmetry[0] + sigma * expected_penalty[0],
                        -expected_consistency_symmetry[1] + sigma * expected_penalty[1],
                        -expected_consistency_symmetry[2] + sigma * expected_penalty[2],
            };

            form.integrateFacetInner2D(  { p0, p1, p2 },  { p0, p2 }, p1,  n1, basis, basis, 0, 1, elMat );

            WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
            WALBERLA_CHECK_FLOAT_EQUAL( elMat( 1, 0 ), expected[1] );
            WALBERLA_CHECK_FLOAT_EQUAL( elMat( 2, 0 ), expected[2] );
         }

         // outer/coupling side of facet integral
         {
            std::vector< real_t > expected_consistency_symmetry = {
             -0.426278031336357, 0.156029997115939, -0.239653917138861};
            std::vector< real_t > expected_penalty = { -0.16666666666666669, 0.0, -0.1333333333333333};
            std::vector< real_t > expected         = {
                        -expected_consistency_symmetry[0] + sigma * expected_penalty[0],
                        -expected_consistency_symmetry[1] + sigma * expected_penalty[1],
                        -expected_consistency_symmetry[2] + sigma * expected_penalty[2],
            };

            form.integrateFacetCoupling( 2, { p0, p1, p2 }, { p0, p2, p3 }, { p2, p0 }, p1, p3, n1, basis, basis, 0, 1, elMat );
            WALBERLA_CHECK_FLOAT_EQUAL( elMat( 0, 0 ), expected[0] );
            WALBERLA_CHECK_FLOAT_EQUAL( elMat( 1, 0 ), expected[1] );
            WALBERLA_CHECK_FLOAT_EQUAL( elMat( 2, 0 ), expected[2] );
         }
      }
   }

}
} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager( &argc, &argv );

   //hyteg::testDivForm();
   //hyteg::testLaplaceForm();
   hyteg::testEpsilonForm();

   return 0;
}

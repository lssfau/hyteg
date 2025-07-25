/*
* Copyright (c) 2017-2025 Nils Kohl, Marcus Mohr.
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

#include "core/DataTypes.h"

#include "hyteg/dgfunctionspace/DGBasisInfo.hpp"
#include "hyteg/forms/DGForm.hpp"
#include "hyteg/forms/form_hyteg_dg/DGForm2D.hpp"
#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"

#include "Eigen/Eigen"

namespace hyteg {
namespace dg {

class DGStokesP1P0PressureStabFormAffine : public DGForm2D
{
 protected:
   void integrateVolume2D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override
   {
      WALBERLA_ASSERT( trialDegree = 0 );
      WALBERLA_ASSERT( testDegree = 0 );

      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

      elMat.setZero();
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >& coordsElement,
                                       const std::vector< Point3D >& coordsFacet,
                                       const Point3D&                oppositeVertex,
                                       const Point3D&                outwardNormal,
                                       const DGBasisInfo&            trialBasis,
                                       const DGBasisInfo&            testBasis,
                                       int                           trialDegree,
                                       int                           testDegree,
                                       MatrixXr&                     elMat ) const override
   {
      WALBERLA_ASSERT( trialDegree = 0 );
      WALBERLA_ASSERT( testDegree = 0 );

      const auto p_affine_0_0 = coordsElement[0]( 0 );
      const auto p_affine_0_1 = coordsElement[0]( 1 );

      const auto p_affine_1_0 = coordsElement[1]( 0 );
      const auto p_affine_1_1 = coordsElement[1]( 1 );

      const auto p_affine_2_0 = coordsElement[2]( 0 );
      const auto p_affine_2_1 = coordsElement[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t tmp_0 = 0.1 * ( std::abs( std::pow( ( ( -p_affine_6_0 + p_affine_7_0 ) * ( -p_affine_6_0 + p_affine_7_0 ) ) +
                                                     ( ( -p_affine_6_1 + p_affine_7_1 ) * ( -p_affine_6_1 + p_affine_7_1 ) ),
                                                 1.0 / 2.0 ) ) *
                             std::abs( std::pow( ( ( -p_affine_6_0 + p_affine_7_0 ) * ( -p_affine_6_0 + p_affine_7_0 ) ) +
                                                     ( ( -p_affine_6_1 + p_affine_7_1 ) * ( -p_affine_6_1 + p_affine_7_1 ) ),
                                                 1.0 / 2.0 ) ) );
      real_t a_0_0 = 1.0000000000000002 * tmp_0;

      elMat( 0, 0 ) = a_0_0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >& coordsElementInner,
                                          const std::vector< Point3D >& coordsElementOuter,
                                          const std::vector< Point3D >& coordsFacet,
                                          const Point3D&                oppositeVertexInnerElement,
                                          const Point3D&                oppositeVertexOuterElement,
                                          const Point3D&                outwardNormal,
                                          const DGBasisInfo&            trialBasis,
                                          const DGBasisInfo&            testBasis,
                                          int                           trialDegree,
                                          int                           testDegree,
                                          MatrixXr&                     elMat ) const override
   {
      WALBERLA_ASSERT( trialDegree = 0 );
      WALBERLA_ASSERT( testDegree = 0 );

      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coordsElementInner[0]( 0 );
      const auto p_affine_0_1 = coordsElementInner[0]( 1 );

      const auto p_affine_1_0 = coordsElementInner[1]( 0 );
      const auto p_affine_1_1 = coordsElementInner[1]( 1 );

      const auto p_affine_2_0 = coordsElementInner[2]( 0 );
      const auto p_affine_2_1 = coordsElementInner[2]( 1 );

      const auto p_affine_3_0 = coordsElementOuter[0]( 0 );
      const auto p_affine_3_1 = coordsElementOuter[0]( 1 );

      const auto p_affine_4_0 = coordsElementOuter[1]( 0 );
      const auto p_affine_4_1 = coordsElementOuter[1]( 1 );

      const auto p_affine_5_0 = coordsElementOuter[2]( 0 );
      const auto p_affine_5_1 = coordsElementOuter[2]( 1 );

      const auto p_affine_6_0 = coordsFacet[0]( 0 );
      const auto p_affine_6_1 = coordsFacet[0]( 1 );

      const auto p_affine_7_0 = coordsFacet[1]( 0 );
      const auto p_affine_7_1 = coordsFacet[1]( 1 );

      const auto p_affine_8_0 = oppositeVertexInnerElement( 0 );
      const auto p_affine_8_1 = oppositeVertexInnerElement( 1 );

      const auto p_affine_9_0 = oppositeVertexOuterElement( 0 );
      const auto p_affine_9_1 = oppositeVertexOuterElement( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t tmp_0 = 0.1 * ( std::abs( std::pow( ( ( -p_affine_6_0 + p_affine_7_0 ) * ( -p_affine_6_0 + p_affine_7_0 ) ) +
                                                     ( ( -p_affine_6_1 + p_affine_7_1 ) * ( -p_affine_6_1 + p_affine_7_1 ) ),
                                                 1.0 / 2.0 ) ) *
                             std::abs( std::pow( ( ( -p_affine_6_0 + p_affine_7_0 ) * ( -p_affine_6_0 + p_affine_7_0 ) ) +
                                                     ( ( -p_affine_6_1 + p_affine_7_1 ) * ( -p_affine_6_1 + p_affine_7_1 ) ),
                                                 1.0 / 2.0 ) ) );
      real_t a_0_0 = -1.0000000000000002 * tmp_0;

      elMat( 0, 0 ) = a_0_0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                                   const std::vector< Point3D >& coordsFacet,
                                                   const Point3D&                oppositeVertex,
                                                   const Point3D&                outwardNormal,
                                                   const DGBasisInfo&            trialBasis,
                                                   const DGBasisInfo&            testBasis,
                                                   int                           trialDegree,
                                                   int                           testDegree,
                                                   MatrixXr&                     elMat ) const override
   {
      WALBERLA_ASSERT( trialDegree = 0 );
      WALBERLA_ASSERT( testDegree = 0 );

      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );
      elMat.setZero();
   }

   virtual void integrateRHSDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                                 const std::vector< Point3D >& coordsFacet,
                                                 const Point3D&                oppositeVertex,
                                                 const Point3D&                outwardNormal,
                                                 const DGBasisInfo&            basis,
                                                 int                           degree,
                                                 MatrixXr&                     elMat ) const override
   {
      WALBERLA_ASSERT( degree = 0 );

      elMat.resize( basis.numDoFsPerElement( 2, degree ), 1 );
      elMat.setZero();
   }
};

} // namespace dg
} // namespace hyteg

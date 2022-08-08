
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

#pragma once

#include "core/DataTypes.h"

namespace hyteg {
    
namespace dg {
    
namespace eg {

class EGVectorMassFormP1E_0 : public hyteg::dg::DGForm2D
{
 protected:
   void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = p_affine_2_0 + tmp_0;
      real_t tmp_3  = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                               p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_4  = 0.025422453185103409 * tmp_3 * ( -0.27024431884183109 * tmp_1 + 0.54048863768366218 * tmp_2 );
      real_t tmp_5  = 0.058393137863189684 * tmp_3 * ( -0.084046588162422886 * tmp_1 + 0.16809317632484583 * tmp_2 );
      real_t tmp_6  = 0.041425537809186785 * tmp_3 * ( -0.022980882299548921 * tmp_1 - 0.28018828348851638 * tmp_2 );
      real_t tmp_7  = 0.041425537809186785 * tmp_3 * ( 0.30316916578806535 * tmp_1 - 0.022980882299548921 * tmp_2 );
      real_t tmp_8  = 0.025422453185103409 * tmp_3 * ( 0.54048863768366218 * tmp_1 - 0.27024431884183109 * tmp_2 );
      real_t tmp_9  = 0.058393137863189684 * tmp_3 * ( 0.16809317632484583 * tmp_1 - 0.084046588162422886 * tmp_2 );
      real_t tmp_10 = 0.025422453185103409 * tmp_3 * ( -0.27024431884183109 * tmp_1 - 0.27024431884183109 * tmp_2 );
      real_t tmp_11 = 0.058393137863189684 * tmp_3 * ( -0.084046588162422886 * tmp_1 - 0.084046588162422886 * tmp_2 );
      real_t tmp_12 = 0.041425537809186785 * tmp_3 * ( -0.022980882299548921 * tmp_1 + 0.30316916578806535 * tmp_2 );
      real_t tmp_13 = 0.041425537809186785 * tmp_3 * ( -0.28018828348851638 * tmp_1 - 0.022980882299548921 * tmp_2 );
      real_t tmp_14 = 0.041425537809186785 * tmp_3 * ( 0.30316916578806535 * tmp_1 - 0.28018828348851638 * tmp_2 );
      real_t tmp_15 = 0.041425537809186785 * tmp_3 * ( -0.28018828348851638 * tmp_1 + 0.30316916578806535 * tmp_2 );
      real_t a_0_0  = 0.87382197101699566 * tmp_10 + 0.50142650965817914 * tmp_11 + 0.053145049844816938 * tmp_12 +
                     0.63650249912139867 * tmp_13 + 0.31035245103378439 * tmp_14 + 0.31035245103378439 * tmp_15 +
                     0.063089014491502282 * tmp_4 + 0.24928674517091043 * tmp_5 + 0.63650249912139867 * tmp_6 +
                     0.053145049844816938 * tmp_7 + 0.063089014491502227 * tmp_8 + 0.24928674517091043 * tmp_9;
      real_t a_1_0 = 0.063089014491502227 * tmp_10 + 0.24928674517091043 * tmp_11 + 0.31035245103378439 * tmp_12 +
                     0.053145049844816938 * tmp_13 + 0.63650249912139867 * tmp_14 + 0.053145049844816938 * tmp_15 +
                     0.063089014491502227 * tmp_4 + 0.24928674517091043 * tmp_5 + 0.31035245103378439 * tmp_6 +
                     0.63650249912139867 * tmp_7 + 0.87382197101699555 * tmp_8 + 0.50142650965817914 * tmp_9;
      real_t a_2_0 = 0.063089014491502227 * tmp_10 + 0.24928674517091043 * tmp_11 + 0.63650249912139867 * tmp_12 +
                     0.31035245103378439 * tmp_13 + 0.053145049844816938 * tmp_14 + 0.63650249912139867 * tmp_15 +
                     0.87382197101699555 * tmp_4 + 0.50142650965817914 * tmp_5 + 0.053145049844816938 * tmp_6 +
                     0.31035245103378439 * tmp_7 + 0.063089014491502227 * tmp_8 + 0.24928674517091043 * tmp_9;
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
   }

   virtual void integrateFacetInner2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                       const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      elMat( 0, 0 ) = 0;
      elMat( 1, 0 ) = 0;
      elMat( 2, 0 ) = 0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementInner,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementOuter,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexInnerElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexOuterElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
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

      elMat( 0, 0 ) = 0;
      elMat( 1, 0 ) = 0;
      elMat( 2, 0 ) = 0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                   const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      elMat( 0, 0 ) = 0;
      elMat( 1, 0 ) = 0;
      elMat( 2, 0 ) = 0;
   }

   void integrateRHSDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                         const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                         const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                         const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                         const DGBasisInfo&                                       basis,
                                         int                                                      degree,
                                         Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
};

class EGVectorMassFormP1E_1 : public hyteg::dg::DGForm2D
{
 protected:
   void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_1_1 + tmp_0;
      real_t tmp_2  = p_affine_2_1 + tmp_0;
      real_t tmp_3  = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                               p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_4  = 0.025422453185103409 * tmp_3 * ( -0.27024431884183109 * tmp_1 + 0.54048863768366218 * tmp_2 );
      real_t tmp_5  = 0.058393137863189684 * tmp_3 * ( -0.084046588162422886 * tmp_1 + 0.16809317632484583 * tmp_2 );
      real_t tmp_6  = 0.041425537809186785 * tmp_3 * ( -0.022980882299548921 * tmp_1 - 0.28018828348851638 * tmp_2 );
      real_t tmp_7  = 0.041425537809186785 * tmp_3 * ( 0.30316916578806535 * tmp_1 - 0.022980882299548921 * tmp_2 );
      real_t tmp_8  = 0.025422453185103409 * tmp_3 * ( 0.54048863768366218 * tmp_1 - 0.27024431884183109 * tmp_2 );
      real_t tmp_9  = 0.058393137863189684 * tmp_3 * ( 0.16809317632484583 * tmp_1 - 0.084046588162422886 * tmp_2 );
      real_t tmp_10 = 0.025422453185103409 * tmp_3 * ( -0.27024431884183109 * tmp_1 - 0.27024431884183109 * tmp_2 );
      real_t tmp_11 = 0.058393137863189684 * tmp_3 * ( -0.084046588162422886 * tmp_1 - 0.084046588162422886 * tmp_2 );
      real_t tmp_12 = 0.041425537809186785 * tmp_3 * ( -0.022980882299548921 * tmp_1 + 0.30316916578806535 * tmp_2 );
      real_t tmp_13 = 0.041425537809186785 * tmp_3 * ( -0.28018828348851638 * tmp_1 - 0.022980882299548921 * tmp_2 );
      real_t tmp_14 = 0.041425537809186785 * tmp_3 * ( 0.30316916578806535 * tmp_1 - 0.28018828348851638 * tmp_2 );
      real_t tmp_15 = 0.041425537809186785 * tmp_3 * ( -0.28018828348851638 * tmp_1 + 0.30316916578806535 * tmp_2 );
      real_t a_0_0  = 0.87382197101699566 * tmp_10 + 0.50142650965817914 * tmp_11 + 0.053145049844816938 * tmp_12 +
                     0.63650249912139867 * tmp_13 + 0.31035245103378439 * tmp_14 + 0.31035245103378439 * tmp_15 +
                     0.063089014491502282 * tmp_4 + 0.24928674517091043 * tmp_5 + 0.63650249912139867 * tmp_6 +
                     0.053145049844816938 * tmp_7 + 0.063089014491502227 * tmp_8 + 0.24928674517091043 * tmp_9;
      real_t a_1_0 = 0.063089014491502227 * tmp_10 + 0.24928674517091043 * tmp_11 + 0.31035245103378439 * tmp_12 +
                     0.053145049844816938 * tmp_13 + 0.63650249912139867 * tmp_14 + 0.053145049844816938 * tmp_15 +
                     0.063089014491502227 * tmp_4 + 0.24928674517091043 * tmp_5 + 0.31035245103378439 * tmp_6 +
                     0.63650249912139867 * tmp_7 + 0.87382197101699555 * tmp_8 + 0.50142650965817914 * tmp_9;
      real_t a_2_0 = 0.063089014491502227 * tmp_10 + 0.24928674517091043 * tmp_11 + 0.63650249912139867 * tmp_12 +
                     0.31035245103378439 * tmp_13 + 0.053145049844816938 * tmp_14 + 0.63650249912139867 * tmp_15 +
                     0.87382197101699555 * tmp_4 + 0.50142650965817914 * tmp_5 + 0.053145049844816938 * tmp_6 +
                     0.31035245103378439 * tmp_7 + 0.063089014491502227 * tmp_8 + 0.24928674517091043 * tmp_9;
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
   }

   virtual void integrateFacetInner2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                       const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      elMat( 0, 0 ) = 0;
      elMat( 1, 0 ) = 0;
      elMat( 2, 0 ) = 0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementInner,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementOuter,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexInnerElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexOuterElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
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

      elMat( 0, 0 ) = 0;
      elMat( 1, 0 ) = 0;
      elMat( 2, 0 ) = 0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                   const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      elMat( 0, 0 ) = 0;
      elMat( 1, 0 ) = 0;
      elMat( 2, 0 ) = 0;
   }

   void integrateRHSDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                         const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                         const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                         const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                         const DGBasisInfo&                                       basis,
                                         int                                                      degree,
                                         Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
};

class EGVectorMassFormEP1_0 : public hyteg::dg::DGForm2D
{
 protected:
   void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = p_affine_2_0 + tmp_0;
      real_t tmp_3  = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                               p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_4  = 0.025422453185103409 * tmp_3 * ( -0.27024431884183109 * tmp_1 + 0.54048863768366218 * tmp_2 );
      real_t tmp_5  = 0.058393137863189684 * tmp_3 * ( -0.084046588162422886 * tmp_1 + 0.16809317632484583 * tmp_2 );
      real_t tmp_6  = 0.041425537809186785 * tmp_3 * ( -0.022980882299548921 * tmp_1 - 0.28018828348851638 * tmp_2 );
      real_t tmp_7  = 0.041425537809186785 * tmp_3 * ( 0.30316916578806535 * tmp_1 - 0.022980882299548921 * tmp_2 );
      real_t tmp_8  = 0.025422453185103409 * tmp_3 * ( 0.54048863768366218 * tmp_1 - 0.27024431884183109 * tmp_2 );
      real_t tmp_9  = 0.058393137863189684 * tmp_3 * ( 0.16809317632484583 * tmp_1 - 0.084046588162422886 * tmp_2 );
      real_t tmp_10 = 0.025422453185103409 * tmp_3 * ( -0.27024431884183109 * tmp_1 - 0.27024431884183109 * tmp_2 );
      real_t tmp_11 = 0.058393137863189684 * tmp_3 * ( -0.084046588162422886 * tmp_1 - 0.084046588162422886 * tmp_2 );
      real_t tmp_12 = 0.041425537809186785 * tmp_3 * ( -0.022980882299548921 * tmp_1 + 0.30316916578806535 * tmp_2 );
      real_t tmp_13 = 0.041425537809186785 * tmp_3 * ( -0.28018828348851638 * tmp_1 - 0.022980882299548921 * tmp_2 );
      real_t tmp_14 = 0.041425537809186785 * tmp_3 * ( 0.30316916578806535 * tmp_1 - 0.28018828348851638 * tmp_2 );
      real_t tmp_15 = 0.041425537809186785 * tmp_3 * ( -0.28018828348851638 * tmp_1 + 0.30316916578806535 * tmp_2 );
      real_t a_0_0  = 0.87382197101699566 * tmp_10 + 0.50142650965817914 * tmp_11 + 0.053145049844816938 * tmp_12 +
                     0.63650249912139867 * tmp_13 + 0.31035245103378439 * tmp_14 + 0.31035245103378439 * tmp_15 +
                     0.063089014491502282 * tmp_4 + 0.24928674517091043 * tmp_5 + 0.63650249912139867 * tmp_6 +
                     0.053145049844816938 * tmp_7 + 0.063089014491502227 * tmp_8 + 0.24928674517091043 * tmp_9;
      real_t a_0_1 = 0.063089014491502227 * tmp_10 + 0.24928674517091043 * tmp_11 + 0.31035245103378439 * tmp_12 +
                     0.053145049844816938 * tmp_13 + 0.63650249912139867 * tmp_14 + 0.053145049844816938 * tmp_15 +
                     0.063089014491502227 * tmp_4 + 0.24928674517091043 * tmp_5 + 0.31035245103378439 * tmp_6 +
                     0.63650249912139867 * tmp_7 + 0.87382197101699555 * tmp_8 + 0.50142650965817914 * tmp_9;
      real_t a_0_2 = 0.063089014491502227 * tmp_10 + 0.24928674517091043 * tmp_11 + 0.63650249912139867 * tmp_12 +
                     0.31035245103378439 * tmp_13 + 0.053145049844816938 * tmp_14 + 0.63650249912139867 * tmp_15 +
                     0.87382197101699555 * tmp_4 + 0.50142650965817914 * tmp_5 + 0.053145049844816938 * tmp_6 +
                     0.31035245103378439 * tmp_7 + 0.063089014491502227 * tmp_8 + 0.24928674517091043 * tmp_9;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
   }

   virtual void integrateFacetInner2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                       const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      elMat( 0, 0 ) = 0;
      elMat( 0, 1 ) = 0;
      elMat( 0, 2 ) = 0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementInner,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementOuter,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexInnerElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexOuterElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
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

      elMat( 0, 0 ) = 0;
      elMat( 0, 1 ) = 0;
      elMat( 0, 2 ) = 0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                   const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      elMat( 0, 0 ) = 0;
      elMat( 0, 1 ) = 0;
      elMat( 0, 2 ) = 0;
   }

   void integrateRHSDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                         const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                         const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                         const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                         const DGBasisInfo&                                       basis,
                                         int                                                      degree,
                                         Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
};

class EGVectorMassFormEP1_1 : public hyteg::dg::DGForm2D
{
 protected:
   void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_1_1 + tmp_0;
      real_t tmp_2  = p_affine_2_1 + tmp_0;
      real_t tmp_3  = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                               p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_4  = 0.025422453185103409 * tmp_3 * ( -0.27024431884183109 * tmp_1 + 0.54048863768366218 * tmp_2 );
      real_t tmp_5  = 0.058393137863189684 * tmp_3 * ( -0.084046588162422886 * tmp_1 + 0.16809317632484583 * tmp_2 );
      real_t tmp_6  = 0.041425537809186785 * tmp_3 * ( -0.022980882299548921 * tmp_1 - 0.28018828348851638 * tmp_2 );
      real_t tmp_7  = 0.041425537809186785 * tmp_3 * ( 0.30316916578806535 * tmp_1 - 0.022980882299548921 * tmp_2 );
      real_t tmp_8  = 0.025422453185103409 * tmp_3 * ( 0.54048863768366218 * tmp_1 - 0.27024431884183109 * tmp_2 );
      real_t tmp_9  = 0.058393137863189684 * tmp_3 * ( 0.16809317632484583 * tmp_1 - 0.084046588162422886 * tmp_2 );
      real_t tmp_10 = 0.025422453185103409 * tmp_3 * ( -0.27024431884183109 * tmp_1 - 0.27024431884183109 * tmp_2 );
      real_t tmp_11 = 0.058393137863189684 * tmp_3 * ( -0.084046588162422886 * tmp_1 - 0.084046588162422886 * tmp_2 );
      real_t tmp_12 = 0.041425537809186785 * tmp_3 * ( -0.022980882299548921 * tmp_1 + 0.30316916578806535 * tmp_2 );
      real_t tmp_13 = 0.041425537809186785 * tmp_3 * ( -0.28018828348851638 * tmp_1 - 0.022980882299548921 * tmp_2 );
      real_t tmp_14 = 0.041425537809186785 * tmp_3 * ( 0.30316916578806535 * tmp_1 - 0.28018828348851638 * tmp_2 );
      real_t tmp_15 = 0.041425537809186785 * tmp_3 * ( -0.28018828348851638 * tmp_1 + 0.30316916578806535 * tmp_2 );
      real_t a_0_0  = 0.87382197101699566 * tmp_10 + 0.50142650965817914 * tmp_11 + 0.053145049844816938 * tmp_12 +
                     0.63650249912139867 * tmp_13 + 0.31035245103378439 * tmp_14 + 0.31035245103378439 * tmp_15 +
                     0.063089014491502282 * tmp_4 + 0.24928674517091043 * tmp_5 + 0.63650249912139867 * tmp_6 +
                     0.053145049844816938 * tmp_7 + 0.063089014491502227 * tmp_8 + 0.24928674517091043 * tmp_9;
      real_t a_0_1 = 0.063089014491502227 * tmp_10 + 0.24928674517091043 * tmp_11 + 0.31035245103378439 * tmp_12 +
                     0.053145049844816938 * tmp_13 + 0.63650249912139867 * tmp_14 + 0.053145049844816938 * tmp_15 +
                     0.063089014491502227 * tmp_4 + 0.24928674517091043 * tmp_5 + 0.31035245103378439 * tmp_6 +
                     0.63650249912139867 * tmp_7 + 0.87382197101699555 * tmp_8 + 0.50142650965817914 * tmp_9;
      real_t a_0_2 = 0.063089014491502227 * tmp_10 + 0.24928674517091043 * tmp_11 + 0.63650249912139867 * tmp_12 +
                     0.31035245103378439 * tmp_13 + 0.053145049844816938 * tmp_14 + 0.63650249912139867 * tmp_15 +
                     0.87382197101699555 * tmp_4 + 0.50142650965817914 * tmp_5 + 0.053145049844816938 * tmp_6 +
                     0.31035245103378439 * tmp_7 + 0.063089014491502227 * tmp_8 + 0.24928674517091043 * tmp_9;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
   }

   virtual void integrateFacetInner2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                       const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      elMat( 0, 0 ) = 0;
      elMat( 0, 1 ) = 0;
      elMat( 0, 2 ) = 0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementInner,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementOuter,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexInnerElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexOuterElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
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

      elMat( 0, 0 ) = 0;
      elMat( 0, 1 ) = 0;
      elMat( 0, 2 ) = 0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                   const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      elMat( 0, 0 ) = 0;
      elMat( 0, 1 ) = 0;
      elMat( 0, 2 ) = 0;
   }

   void integrateRHSDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                         const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                         const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                         const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                         const DGBasisInfo&                                       basis,
                                         int                                                      degree,
                                         Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
};

class EGVectorMassFormEE : public hyteg::dg::DGForm2D
{
 protected:
   void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t tmp_0  = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                               p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_1  = -0.27024431884183109;
      real_t tmp_2  = -p_affine_0_0;
      real_t tmp_3  = p_affine_1_0 + tmp_2;
      real_t tmp_4  = 0.54048863768366218;
      real_t tmp_5  = p_affine_2_0 + tmp_2;
      real_t tmp_6  = -p_affine_0_1;
      real_t tmp_7  = p_affine_1_1 + tmp_6;
      real_t tmp_8  = p_affine_2_1 + tmp_6;
      real_t tmp_9  = -0.084046588162422886;
      real_t tmp_10 = 0.16809317632484583;
      real_t tmp_11 = -0.022980882299548921;
      real_t tmp_12 = -0.28018828348851638;
      real_t tmp_13 = 0.30316916578806535;
      real_t tmp_14 = -0.022980882299548921;
      real_t tmp_15 = 0.54048863768366218;
      real_t tmp_16 = -0.27024431884183109;
      real_t tmp_17 = 0.16809317632484583;
      real_t tmp_18 = -0.084046588162422886;
      real_t tmp_19 = -0.27024431884183109;
      real_t tmp_20 = -0.27024431884183109;
      real_t tmp_21 = -0.084046588162422886;
      real_t tmp_22 = -0.084046588162422886;
      real_t tmp_23 = -0.022980882299548921;
      real_t tmp_24 = 0.30316916578806535;
      real_t tmp_25 = -0.28018828348851638;
      real_t tmp_26 = -0.022980882299548921;
      real_t tmp_27 = 0.30316916578806535;
      real_t tmp_28 = -0.28018828348851638;
      real_t tmp_29 = -0.28018828348851638;
      real_t tmp_30 = 0.30316916578806535;
      real_t a_0_0  = 0.025422453185103409 * tmp_0 *
                         ( ( ( tmp_1 * tmp_3 + tmp_4 * tmp_5 ) * ( tmp_1 * tmp_3 + tmp_4 * tmp_5 ) ) +
                           ( ( tmp_1 * tmp_7 + tmp_4 * tmp_8 ) * ( tmp_1 * tmp_7 + tmp_4 * tmp_8 ) ) ) +
                     0.058393137863189684 * tmp_0 *
                         ( ( ( tmp_10 * tmp_5 + tmp_3 * tmp_9 ) * ( tmp_10 * tmp_5 + tmp_3 * tmp_9 ) ) +
                           ( ( tmp_10 * tmp_8 + tmp_7 * tmp_9 ) * ( tmp_10 * tmp_8 + tmp_7 * tmp_9 ) ) ) +
                     0.041425537809186785 * tmp_0 *
                         ( ( ( tmp_11 * tmp_3 + tmp_12 * tmp_5 ) * ( tmp_11 * tmp_3 + tmp_12 * tmp_5 ) ) +
                           ( ( tmp_11 * tmp_7 + tmp_12 * tmp_8 ) * ( tmp_11 * tmp_7 + tmp_12 * tmp_8 ) ) ) +
                     0.041425537809186785 * tmp_0 *
                         ( ( ( tmp_13 * tmp_3 + tmp_14 * tmp_5 ) * ( tmp_13 * tmp_3 + tmp_14 * tmp_5 ) ) +
                           ( ( tmp_13 * tmp_7 + tmp_14 * tmp_8 ) * ( tmp_13 * tmp_7 + tmp_14 * tmp_8 ) ) ) +
                     0.025422453185103409 * tmp_0 *
                         ( ( ( tmp_15 * tmp_3 + tmp_16 * tmp_5 ) * ( tmp_15 * tmp_3 + tmp_16 * tmp_5 ) ) +
                           ( ( tmp_15 * tmp_7 + tmp_16 * tmp_8 ) * ( tmp_15 * tmp_7 + tmp_16 * tmp_8 ) ) ) +
                     0.058393137863189684 * tmp_0 *
                         ( ( ( tmp_17 * tmp_3 + tmp_18 * tmp_5 ) * ( tmp_17 * tmp_3 + tmp_18 * tmp_5 ) ) +
                           ( ( tmp_17 * tmp_7 + tmp_18 * tmp_8 ) * ( tmp_17 * tmp_7 + tmp_18 * tmp_8 ) ) ) +
                     0.025422453185103409 * tmp_0 *
                         ( ( ( tmp_19 * tmp_3 + tmp_20 * tmp_5 ) * ( tmp_19 * tmp_3 + tmp_20 * tmp_5 ) ) +
                           ( ( tmp_19 * tmp_7 + tmp_20 * tmp_8 ) * ( tmp_19 * tmp_7 + tmp_20 * tmp_8 ) ) ) +
                     0.058393137863189684 * tmp_0 *
                         ( ( ( tmp_21 * tmp_3 + tmp_22 * tmp_5 ) * ( tmp_21 * tmp_3 + tmp_22 * tmp_5 ) ) +
                           ( ( tmp_21 * tmp_7 + tmp_22 * tmp_8 ) * ( tmp_21 * tmp_7 + tmp_22 * tmp_8 ) ) ) +
                     0.041425537809186785 * tmp_0 *
                         ( ( ( tmp_23 * tmp_3 + tmp_24 * tmp_5 ) * ( tmp_23 * tmp_3 + tmp_24 * tmp_5 ) ) +
                           ( ( tmp_23 * tmp_7 + tmp_24 * tmp_8 ) * ( tmp_23 * tmp_7 + tmp_24 * tmp_8 ) ) ) +
                     0.041425537809186785 * tmp_0 *
                         ( ( ( tmp_25 * tmp_3 + tmp_26 * tmp_5 ) * ( tmp_25 * tmp_3 + tmp_26 * tmp_5 ) ) +
                           ( ( tmp_25 * tmp_7 + tmp_26 * tmp_8 ) * ( tmp_25 * tmp_7 + tmp_26 * tmp_8 ) ) ) +
                     0.041425537809186785 * tmp_0 *
                         ( ( ( tmp_27 * tmp_3 + tmp_28 * tmp_5 ) * ( tmp_27 * tmp_3 + tmp_28 * tmp_5 ) ) +
                           ( ( tmp_27 * tmp_7 + tmp_28 * tmp_8 ) * ( tmp_27 * tmp_7 + tmp_28 * tmp_8 ) ) ) +
                     0.041425537809186785 * tmp_0 *
                         ( ( ( tmp_29 * tmp_3 + tmp_30 * tmp_5 ) * ( tmp_29 * tmp_3 + tmp_30 * tmp_5 ) ) +
                           ( ( tmp_29 * tmp_7 + tmp_30 * tmp_8 ) * ( tmp_29 * tmp_7 + tmp_30 * tmp_8 ) ) );
      elMat( 0, 0 ) = a_0_0;
   }

   virtual void integrateFacetInner2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                       const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      elMat( 0, 0 ) = 0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementInner,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementOuter,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexInnerElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexOuterElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
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

      elMat( 0, 0 ) = 0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                   const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      elMat( 0, 0 ) = 0;
   }

   void integrateRHSDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                         const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                         const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                         const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                         const DGBasisInfo&                                       basis,
                                         int                                                      degree,
                                         Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
};

}
}
}
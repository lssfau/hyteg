
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

# pragma once

# include "core/DataTypes.h"

# include "hyteg/dgfunctionspace/DGBasisInfo.hpp"
# include "hyteg/dgfunctionspace/DGForm.hpp"
# include "hyteg/dgfunctionspace/DGForm2D.hpp"
#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"


# include "Eigen/Eigen"

namespace hyteg {
namespace dg{
namespace eg{

class EGDivFormNitscheBC_P0P1_0 : public hyteg::dg::DGForm
{

 public:
    EGDivFormNitscheBC_P0P1_0()
: callback_Scalar_Variable_Coefficient_3D_g1 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_2D_g1 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_2D_g0 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_3D_g0 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_3D_g2 ([](const Point3D & p) -> real_t { return 0.; })
    {}

void Scalar_Variable_Coefficient_2D_g0( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g0( Point3D( {in_0, in_1, 0} ) );
}
void Scalar_Variable_Coefficient_2D_g1( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g1( Point3D( {in_0, in_1, 0} ) );
}

void Scalar_Variable_Coefficient_3D_g0( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g0( Point3D( {in_0, in_1, in_2} ) );
}
void Scalar_Variable_Coefficient_3D_g1( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g1( Point3D( {in_0, in_1, in_2} ) );
}
void Scalar_Variable_Coefficient_3D_g2( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g2( Point3D( {in_0, in_1, in_2} ) );
}

 protected:
  void integrateVolume2D( const std::vector< Point3D >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t tmp_0 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_1 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_2 = 1.0 / (tmp_0*(-p_affine_0_0 + p_affine_1_0) + tmp_1*(-p_affine_0_0 + p_affine_2_0));
      real_t tmp_3 = tmp_0*tmp_2;
      real_t tmp_4 = tmp_1*tmp_2;
      real_t tmp_5 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_6 = tmp_5*(tmp_3 + tmp_4);
      real_t tmp_7 = tmp_3*tmp_5;
      real_t tmp_8 = tmp_4*tmp_5;
      real_t a_0_0 = 0.5*tmp_6;
      real_t a_0_1 = -0.5*tmp_7;
      real_t a_0_2 = -0.5*tmp_8;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >&      coordsElement,
                                       const std::vector< Point3D >&      coordsFacet,
                                       const Point3D&                     oppositeVertex,
                                       const Point3D&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      real_t tmp_0 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_3 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_4 = 1.0 / (-tmp_0*tmp_3 + tmp_1*tmp_2);
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_9 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_10 = tmp_4*(0.046910077030668018*tmp_8 + tmp_9);
      real_t tmp_11 = tmp_0*tmp_7 + tmp_10*tmp_2;
      real_t tmp_12 = tmp_1*tmp_7 + tmp_10*tmp_3;
      real_t tmp_13 = 0.5*p_affine_10_0*std::abs(std::pow((tmp_5*tmp_5) + (tmp_8*tmp_8), 1.0/2.0));
      real_t tmp_14 = 0.11846344252809471*tmp_13;
      real_t tmp_15 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_16 = tmp_4*(0.23076534494715845*tmp_8 + tmp_9);
      real_t tmp_17 = tmp_0*tmp_15 + tmp_16*tmp_2;
      real_t tmp_18 = tmp_1*tmp_15 + tmp_16*tmp_3;
      real_t tmp_19 = 0.2393143352496831*tmp_13;
      real_t tmp_20 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_21 = tmp_4*(0.5*tmp_8 + tmp_9);
      real_t tmp_22 = tmp_0*tmp_20 + tmp_2*tmp_21;
      real_t tmp_23 = tmp_1*tmp_20 + tmp_21*tmp_3;
      real_t tmp_24 = 0.2844444444444445*tmp_13;
      real_t tmp_25 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_26 = tmp_4*(0.7692346550528415*tmp_8 + tmp_9);
      real_t tmp_27 = tmp_0*tmp_25 + tmp_2*tmp_26;
      real_t tmp_28 = tmp_1*tmp_25 + tmp_26*tmp_3;
      real_t tmp_29 = 0.2393143352496831*tmp_13;
      real_t tmp_30 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_31 = tmp_4*(0.95308992296933193*tmp_8 + tmp_9);
      real_t tmp_32 = tmp_0*tmp_30 + tmp_2*tmp_31;
      real_t tmp_33 = tmp_1*tmp_30 + tmp_3*tmp_31;
      real_t tmp_34 = 0.11846344252809471*tmp_13;
      real_t a_0_0 = tmp_14*(-tmp_11 - tmp_12 + 1) + tmp_19*(-tmp_17 - tmp_18 + 1) + tmp_24*(-tmp_22 - tmp_23 + 1) + tmp_29*(-tmp_27 - tmp_28 + 1) + tmp_34*(-tmp_32 - tmp_33 + 1);
      real_t a_0_1 = tmp_11*tmp_14 + tmp_17*tmp_19 + tmp_22*tmp_24 + tmp_27*tmp_29 + tmp_32*tmp_34;
      real_t a_0_2 = tmp_12*tmp_14 + tmp_18*tmp_19 + tmp_23*tmp_24 + tmp_28*tmp_29 + tmp_33*tmp_34;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >&      coordsElementInner,
                                          const std::vector< Point3D >&      coordsElementOuter,
                                          const std::vector< Point3D >&      coordsFacet,
                                          const Point3D&                     oppositeVertexInnerElement,
                                          const Point3D&                     oppositeVertexOuterElement,
                                          const Point3D&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      real_t tmp_0 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_1 = -p_affine_3_0 + p_affine_4_0;
      real_t tmp_2 = -p_affine_3_1 + p_affine_5_1;
      real_t tmp_3 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_4 = 1.0 / (-tmp_0*tmp_3 + tmp_1*tmp_2);
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = -p_affine_3_1 + p_affine_6_1;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_9 = -p_affine_3_0 + p_affine_6_0;
      real_t tmp_10 = tmp_4*(0.046910077030668018*tmp_8 + tmp_9);
      real_t tmp_11 = tmp_0*tmp_7 + tmp_10*tmp_2;
      real_t tmp_12 = tmp_1*tmp_7 + tmp_10*tmp_3;
      real_t tmp_13 = 0.5*p_affine_10_0*std::abs(std::pow((tmp_5*tmp_5) + (tmp_8*tmp_8), 1.0/2.0));
      real_t tmp_14 = 0.11846344252809471*tmp_13;
      real_t tmp_15 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_16 = tmp_4*(0.23076534494715845*tmp_8 + tmp_9);
      real_t tmp_17 = tmp_0*tmp_15 + tmp_16*tmp_2;
      real_t tmp_18 = tmp_1*tmp_15 + tmp_16*tmp_3;
      real_t tmp_19 = 0.2393143352496831*tmp_13;
      real_t tmp_20 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_21 = tmp_4*(0.5*tmp_8 + tmp_9);
      real_t tmp_22 = tmp_0*tmp_20 + tmp_2*tmp_21;
      real_t tmp_23 = tmp_1*tmp_20 + tmp_21*tmp_3;
      real_t tmp_24 = 0.2844444444444445*tmp_13;
      real_t tmp_25 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_26 = tmp_4*(0.7692346550528415*tmp_8 + tmp_9);
      real_t tmp_27 = tmp_0*tmp_25 + tmp_2*tmp_26;
      real_t tmp_28 = tmp_1*tmp_25 + tmp_26*tmp_3;
      real_t tmp_29 = 0.2393143352496831*tmp_13;
      real_t tmp_30 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_31 = tmp_4*(0.95308992296933193*tmp_8 + tmp_9);
      real_t tmp_32 = tmp_0*tmp_30 + tmp_2*tmp_31;
      real_t tmp_33 = tmp_1*tmp_30 + tmp_3*tmp_31;
      real_t tmp_34 = 0.11846344252809471*tmp_13;
      real_t a_0_0 = -tmp_14*(-tmp_11 - tmp_12 + 1) - tmp_19*(-tmp_17 - tmp_18 + 1) - tmp_24*(-tmp_22 - tmp_23 + 1) - tmp_29*(-tmp_27 - tmp_28 + 1) - tmp_34*(-tmp_32 - tmp_33 + 1);
      real_t a_0_1 = -tmp_11*tmp_14 - tmp_17*tmp_19 - tmp_22*tmp_24 - tmp_27*tmp_29 - tmp_32*tmp_34;
      real_t a_0_2 = -tmp_12*tmp_14 - tmp_18*tmp_19 - tmp_23*tmp_24 - tmp_28*tmp_29 - tmp_33*tmp_34;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                   const std::vector< Point3D >&      coordsFacet,
                                                   const Point3D&                     oppositeVertex,
                                                   const Point3D&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      real_t tmp_0 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_3 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_4 = 1.0 / (-tmp_0*tmp_3 + tmp_1*tmp_2);
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_9 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_10 = tmp_4*(0.046910077030668018*tmp_8 + tmp_9);
      real_t tmp_11 = tmp_0*tmp_7 + tmp_10*tmp_2;
      real_t tmp_12 = tmp_1*tmp_7 + tmp_10*tmp_3;
      real_t tmp_13 = p_affine_10_0*std::abs(std::pow((tmp_5*tmp_5) + (tmp_8*tmp_8), 1.0/2.0));
      real_t tmp_14 = 0.11846344252809471*tmp_13;
      real_t tmp_15 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_16 = tmp_4*(0.23076534494715845*tmp_8 + tmp_9);
      real_t tmp_17 = tmp_0*tmp_15 + tmp_16*tmp_2;
      real_t tmp_18 = tmp_1*tmp_15 + tmp_16*tmp_3;
      real_t tmp_19 = 0.2393143352496831*tmp_13;
      real_t tmp_20 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_21 = tmp_4*(0.5*tmp_8 + tmp_9);
      real_t tmp_22 = tmp_0*tmp_20 + tmp_2*tmp_21;
      real_t tmp_23 = tmp_1*tmp_20 + tmp_21*tmp_3;
      real_t tmp_24 = 0.2844444444444445*tmp_13;
      real_t tmp_25 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_26 = tmp_4*(0.7692346550528415*tmp_8 + tmp_9);
      real_t tmp_27 = tmp_0*tmp_25 + tmp_2*tmp_26;
      real_t tmp_28 = tmp_1*tmp_25 + tmp_26*tmp_3;
      real_t tmp_29 = 0.2393143352496831*tmp_13;
      real_t tmp_30 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_31 = tmp_4*(0.95308992296933193*tmp_8 + tmp_9);
      real_t tmp_32 = tmp_0*tmp_30 + tmp_2*tmp_31;
      real_t tmp_33 = tmp_1*tmp_30 + tmp_3*tmp_31;
      real_t tmp_34 = 0.11846344252809471*tmp_13;
      real_t a_0_0 = tmp_14*(-tmp_11 - tmp_12 + 1) + tmp_19*(-tmp_17 - tmp_18 + 1) + tmp_24*(-tmp_22 - tmp_23 + 1) + tmp_29*(-tmp_27 - tmp_28 + 1) + tmp_34*(-tmp_32 - tmp_33 + 1);
      real_t a_0_1 = tmp_11*tmp_14 + tmp_17*tmp_19 + tmp_22*tmp_24 + tmp_27*tmp_29 + tmp_32*tmp_34;
      real_t a_0_2 = tmp_12*tmp_14 + tmp_18*tmp_19 + tmp_23*tmp_24 + tmp_28*tmp_29 + tmp_33*tmp_34;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

    void integrateRHSDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

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

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_g0_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id9 = 0;
      Scalar_Variable_Coefficient_2D_g0( 0.95308992296933193*p_affine_6_0 + 0.046910077030668018*p_affine_7_0, 0.95308992296933193*p_affine_6_1 + 0.046910077030668018*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id0 );
      Scalar_Variable_Coefficient_2D_g1( 0.95308992296933193*p_affine_6_0 + 0.046910077030668018*p_affine_7_0, 0.95308992296933193*p_affine_6_1 + 0.046910077030668018*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id1 );
      Scalar_Variable_Coefficient_2D_g0( 0.7692346550528415*p_affine_6_0 + 0.23076534494715845*p_affine_7_0, 0.7692346550528415*p_affine_6_1 + 0.23076534494715845*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id2 );
      Scalar_Variable_Coefficient_2D_g1( 0.7692346550528415*p_affine_6_0 + 0.23076534494715845*p_affine_7_0, 0.7692346550528415*p_affine_6_1 + 0.23076534494715845*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id3 );
      Scalar_Variable_Coefficient_2D_g0( 0.5*p_affine_6_0 + 0.5*p_affine_7_0, 0.5*p_affine_6_1 + 0.5*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id4 );
      Scalar_Variable_Coefficient_2D_g1( 0.5*p_affine_6_0 + 0.5*p_affine_7_0, 0.5*p_affine_6_1 + 0.5*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id5 );
      Scalar_Variable_Coefficient_2D_g0( 0.2307653449471585*p_affine_6_0 + 0.7692346550528415*p_affine_7_0, 0.2307653449471585*p_affine_6_1 + 0.7692346550528415*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id6 );
      Scalar_Variable_Coefficient_2D_g1( 0.2307653449471585*p_affine_6_0 + 0.7692346550528415*p_affine_7_0, 0.2307653449471585*p_affine_6_1 + 0.7692346550528415*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id7 );
      Scalar_Variable_Coefficient_2D_g0( 0.046910077030668074*p_affine_6_0 + 0.95308992296933193*p_affine_7_0, 0.046910077030668074*p_affine_6_1 + 0.95308992296933193*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id8 );
      Scalar_Variable_Coefficient_2D_g1( 0.046910077030668074*p_affine_6_0 + 0.95308992296933193*p_affine_7_0, 0.046910077030668074*p_affine_6_1 + 0.95308992296933193*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id9 );
      real_t tmp_0 = std::abs(std::pow(((-p_affine_6_0 + p_affine_7_0)*(-p_affine_6_0 + p_affine_7_0)) + ((-p_affine_6_1 + p_affine_7_1)*(-p_affine_6_1 + p_affine_7_1)), 1.0/2.0));
      real_t a_0_0 = 0.11846344252809471*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id0*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id1*p_affine_10_1) + 0.2393143352496831*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id2*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id3*p_affine_10_1) + 0.2844444444444445*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id4*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id5*p_affine_10_1) + 0.2393143352496831*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id6*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id7*p_affine_10_1) + 0.11846344252809471*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id8*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id9*p_affine_10_1);
      elMat( 0, 0) = a_0_0;
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

      real_t Scalar_Variable_Coefficient_3D_g0_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id13 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id14 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id15 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id16 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id17 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id18 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id19 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id20 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id21 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id22 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id23 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id24 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id25 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id26 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id27 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id28 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id29 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id30 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id31 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id32 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id33 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id34 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id35 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id36 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id37 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id38 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id39 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id40 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id41 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id42 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id43 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id44 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id45 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id46 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id47 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id48 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id49 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id50 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id51 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id52 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id53 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id54 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id55 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id56 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id57 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id58 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id59 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id60 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id61 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id62 = 0;
      Scalar_Variable_Coefficient_3D_g0( 0.93718850182767688*p_affine_10_0 + 0.031405749086161561*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.93718850182767688*p_affine_10_1 + 0.031405749086161561*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.93718850182767688*p_affine_10_2 + 0.031405749086161561*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id0 );
      Scalar_Variable_Coefficient_3D_g1( 0.93718850182767688*p_affine_10_0 + 0.031405749086161561*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.93718850182767688*p_affine_10_1 + 0.031405749086161561*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.93718850182767688*p_affine_10_2 + 0.031405749086161561*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id1 );
      Scalar_Variable_Coefficient_3D_g2( 0.93718850182767688*p_affine_10_0 + 0.031405749086161561*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.93718850182767688*p_affine_10_1 + 0.031405749086161561*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.93718850182767688*p_affine_10_2 + 0.031405749086161561*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id2 );
      Scalar_Variable_Coefficient_3D_g0( 0.60796128279561268*p_affine_10_0 + 0.19601935860219366*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.60796128279561268*p_affine_10_1 + 0.19601935860219366*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.60796128279561268*p_affine_10_2 + 0.19601935860219366*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id3 );
      Scalar_Variable_Coefficient_3D_g1( 0.60796128279561268*p_affine_10_0 + 0.19601935860219366*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.60796128279561268*p_affine_10_1 + 0.19601935860219366*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.60796128279561268*p_affine_10_2 + 0.19601935860219366*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id4 );
      Scalar_Variable_Coefficient_3D_g2( 0.60796128279561268*p_affine_10_0 + 0.19601935860219366*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.60796128279561268*p_affine_10_1 + 0.19601935860219366*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.60796128279561268*p_affine_10_2 + 0.19601935860219366*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id5 );
      Scalar_Variable_Coefficient_3D_g0( 0.039308471900058539*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id6 );
      Scalar_Variable_Coefficient_3D_g1( 0.039308471900058539*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id7 );
      Scalar_Variable_Coefficient_3D_g2( 0.039308471900058539*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id8 );
      Scalar_Variable_Coefficient_3D_g0( 0.1711304259088916*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id9 );
      Scalar_Variable_Coefficient_3D_g1( 0.1711304259088916*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id10 );
      Scalar_Variable_Coefficient_3D_g2( 0.1711304259088916*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id11 );
      Scalar_Variable_Coefficient_3D_g0( 0.37605877282253791*p_affine_10_0 + 0.039308471900058539*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.039308471900058539*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.039308471900058539*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id12 );
      Scalar_Variable_Coefficient_3D_g1( 0.37605877282253791*p_affine_10_0 + 0.039308471900058539*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.039308471900058539*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.039308471900058539*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id13 );
      Scalar_Variable_Coefficient_3D_g2( 0.37605877282253791*p_affine_10_0 + 0.039308471900058539*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.039308471900058539*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.039308471900058539*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id14 );
      Scalar_Variable_Coefficient_3D_g0( 0.78764240869137092*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id15 );
      Scalar_Variable_Coefficient_3D_g1( 0.78764240869137092*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id16 );
      Scalar_Variable_Coefficient_3D_g2( 0.78764240869137092*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id17 );
      Scalar_Variable_Coefficient_3D_g0( 0.58463275527740355*p_affine_10_0 + 0.37605877282253797*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.37605877282253797*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.37605877282253797*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id18 );
      Scalar_Variable_Coefficient_3D_g1( 0.58463275527740355*p_affine_10_0 + 0.37605877282253797*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.37605877282253797*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.37605877282253797*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id19 );
      Scalar_Variable_Coefficient_3D_g2( 0.58463275527740355*p_affine_10_0 + 0.37605877282253797*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.37605877282253797*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.37605877282253797*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id20 );
      Scalar_Variable_Coefficient_3D_g0( 0.041227165399737475*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id21 );
      Scalar_Variable_Coefficient_3D_g1( 0.041227165399737475*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id22 );
      Scalar_Variable_Coefficient_3D_g2( 0.041227165399737475*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id23 );
      Scalar_Variable_Coefficient_3D_g0( 0.039308471900058539*p_affine_10_0 + 0.37605877282253791*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.37605877282253791*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.37605877282253791*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id24 );
      Scalar_Variable_Coefficient_3D_g1( 0.039308471900058539*p_affine_10_0 + 0.37605877282253791*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.37605877282253791*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.37605877282253791*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id25 );
      Scalar_Variable_Coefficient_3D_g2( 0.039308471900058539*p_affine_10_0 + 0.37605877282253791*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.37605877282253791*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.37605877282253791*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id26 );
      Scalar_Variable_Coefficient_3D_g0( 0.78764240869137092*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id27 );
      Scalar_Variable_Coefficient_3D_g1( 0.78764240869137092*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id28 );
      Scalar_Variable_Coefficient_3D_g2( 0.78764240869137092*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id29 );
      Scalar_Variable_Coefficient_3D_g0( 0.58463275527740355*p_affine_10_0 + 0.039308471900058484*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.039308471900058484*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.039308471900058484*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id30 );
      Scalar_Variable_Coefficient_3D_g1( 0.58463275527740355*p_affine_10_0 + 0.039308471900058484*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.039308471900058484*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.039308471900058484*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id31 );
      Scalar_Variable_Coefficient_3D_g2( 0.58463275527740355*p_affine_10_0 + 0.039308471900058484*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.039308471900058484*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.039308471900058484*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id32 );
      Scalar_Variable_Coefficient_3D_g0( 0.1711304259088916*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id33 );
      Scalar_Variable_Coefficient_3D_g1( 0.1711304259088916*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id34 );
      Scalar_Variable_Coefficient_3D_g2( 0.1711304259088916*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id35 );
      Scalar_Variable_Coefficient_3D_g0( 0.19107600050469298*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.19107600050469298*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.19107600050469298*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id36 );
      Scalar_Variable_Coefficient_3D_g1( 0.19107600050469298*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.19107600050469298*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.19107600050469298*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id37 );
      Scalar_Variable_Coefficient_3D_g2( 0.19107600050469298*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.19107600050469298*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.19107600050469298*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id38 );
      Scalar_Variable_Coefficient_3D_g0( 0.37605877282253791*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id39 );
      Scalar_Variable_Coefficient_3D_g1( 0.37605877282253791*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id40 );
      Scalar_Variable_Coefficient_3D_g2( 0.37605877282253791*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id41 );
      Scalar_Variable_Coefficient_3D_g0( 0.031405749086161582*p_affine_10_0 + 0.03140574908616154*p_affine_8_0 + 0.93718850182767688*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.03140574908616154*p_affine_8_1 + 0.93718850182767688*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.03140574908616154*p_affine_8_2 + 0.93718850182767688*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id42 );
      Scalar_Variable_Coefficient_3D_g1( 0.031405749086161582*p_affine_10_0 + 0.03140574908616154*p_affine_8_0 + 0.93718850182767688*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.03140574908616154*p_affine_8_1 + 0.93718850182767688*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.03140574908616154*p_affine_8_2 + 0.93718850182767688*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id43 );
      Scalar_Variable_Coefficient_3D_g2( 0.031405749086161582*p_affine_10_0 + 0.03140574908616154*p_affine_8_0 + 0.93718850182767688*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.03140574908616154*p_affine_8_1 + 0.93718850182767688*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.03140574908616154*p_affine_8_2 + 0.93718850182767688*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id44 );
      Scalar_Variable_Coefficient_3D_g0( 0.19601935860219369*p_affine_10_0 + 0.19601935860219363*p_affine_8_0 + 0.60796128279561268*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.19601935860219363*p_affine_8_1 + 0.60796128279561268*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.19601935860219363*p_affine_8_2 + 0.60796128279561268*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id45 );
      Scalar_Variable_Coefficient_3D_g1( 0.19601935860219369*p_affine_10_0 + 0.19601935860219363*p_affine_8_0 + 0.60796128279561268*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.19601935860219363*p_affine_8_1 + 0.60796128279561268*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.19601935860219363*p_affine_8_2 + 0.60796128279561268*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id46 );
      Scalar_Variable_Coefficient_3D_g2( 0.19601935860219369*p_affine_10_0 + 0.19601935860219363*p_affine_8_0 + 0.60796128279561268*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.19601935860219363*p_affine_8_1 + 0.60796128279561268*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.19601935860219363*p_affine_8_2 + 0.60796128279561268*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id47 );
      Scalar_Variable_Coefficient_3D_g0( 0.40446199974765351*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.19107600050469298*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.19107600050469298*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.19107600050469298*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id48 );
      Scalar_Variable_Coefficient_3D_g1( 0.40446199974765351*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.19107600050469298*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.19107600050469298*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.19107600050469298*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id49 );
      Scalar_Variable_Coefficient_3D_g2( 0.40446199974765351*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.19107600050469298*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.19107600050469298*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.19107600050469298*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id50 );
      Scalar_Variable_Coefficient_3D_g0( 0.031405749086161582*p_affine_10_0 + 0.93718850182767688*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.93718850182767688*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.93718850182767688*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id51 );
      Scalar_Variable_Coefficient_3D_g1( 0.031405749086161582*p_affine_10_0 + 0.93718850182767688*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.93718850182767688*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.93718850182767688*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id52 );
      Scalar_Variable_Coefficient_3D_g2( 0.031405749086161582*p_affine_10_0 + 0.93718850182767688*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.93718850182767688*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.93718850182767688*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id53 );
      Scalar_Variable_Coefficient_3D_g0( 0.19601935860219369*p_affine_10_0 + 0.60796128279561268*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.60796128279561268*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.60796128279561268*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id54 );
      Scalar_Variable_Coefficient_3D_g1( 0.19601935860219369*p_affine_10_0 + 0.60796128279561268*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.60796128279561268*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.60796128279561268*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id55 );
      Scalar_Variable_Coefficient_3D_g2( 0.19601935860219369*p_affine_10_0 + 0.60796128279561268*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.60796128279561268*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.60796128279561268*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id56 );
      Scalar_Variable_Coefficient_3D_g0( 0.40446199974765351*p_affine_10_0 + 0.19107600050469298*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.19107600050469298*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.19107600050469298*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id57 );
      Scalar_Variable_Coefficient_3D_g1( 0.40446199974765351*p_affine_10_0 + 0.19107600050469298*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.19107600050469298*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.19107600050469298*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id58 );
      Scalar_Variable_Coefficient_3D_g2( 0.40446199974765351*p_affine_10_0 + 0.19107600050469298*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.19107600050469298*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.19107600050469298*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id59 );
      Scalar_Variable_Coefficient_3D_g0( 0.041227165399737475*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id60 );
      Scalar_Variable_Coefficient_3D_g1( 0.041227165399737475*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id61 );
      Scalar_Variable_Coefficient_3D_g2( 0.041227165399737475*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id62 );
      real_t tmp_0 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_1 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_2 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_3 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_4 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_5 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_6 = 1.0*std::pow((std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)*std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)) + (std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)*std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)) + (std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)*std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)), 1.0/2.0);
      real_t a_0_0 = 0.0068572537431980923*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id0*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id1*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id2*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id12*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id13*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id14*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id15*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id16*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id17*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id18*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id19*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id20*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id21*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id22*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id23*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id24*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id25*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id26*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id27*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id28*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id29*p_affine_13_2) + 0.037198804536718075*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id3*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id4*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id5*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id30*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id31*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id32*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id33*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id34*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id35*p_affine_13_2) + 0.042507265838595799*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id36*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id37*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id38*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id39*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id40*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id41*p_affine_13_2) + 0.0068572537431980923*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id42*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id43*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id44*p_affine_13_2) + 0.037198804536718075*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id45*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id46*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id47*p_affine_13_2) + 0.042507265838595799*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id48*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id49*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id50*p_affine_13_2) + 0.0068572537431980923*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id51*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id52*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id53*p_affine_13_2) + 0.037198804536718075*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id54*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id55*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id56*p_affine_13_2) + 0.042507265838595799*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id57*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id58*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id59*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id6*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id7*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id8*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id60*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id61*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id62*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id9*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id10*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id11*p_affine_13_2);
      elMat( 0, 0) = a_0_0;
   }
   void integrateVolume3D( const std::vector< Point3D >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coords[0]( 0 );
   const auto p_affine_0_1 = coords[0]( 1 );
   const auto p_affine_0_2 = coords[0]( 2 );

   const auto p_affine_1_0 = coords[1]( 0 );
   const auto p_affine_1_1 = coords[1]( 1 );
   const auto p_affine_1_2 = coords[1]( 2 );

   const auto p_affine_2_0 = coords[2]( 0 );
   const auto p_affine_2_1 = coords[2]( 1 );
   const auto p_affine_2_2 = coords[2]( 2 );

   const auto p_affine_3_0 = coords[3]( 0 );
   const auto p_affine_3_1 = coords[3]( 1 );
   const auto p_affine_3_2 = coords[3]( 2 );

      real_t tmp_0 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_1 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_8 = tmp_3*tmp_7;
      real_t tmp_9 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_10 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_1*tmp_10;
      real_t tmp_13 = tmp_0*tmp_7;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_4*tmp_9 + tmp_11*tmp_2 - tmp_11*tmp_5 - tmp_12*tmp_6 - tmp_13*tmp_9 + tmp_6*tmp_8);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_10*tmp_4 - tmp_13);
      real_t tmp_17 = tmp_14*(-tmp_12 + tmp_8);
      real_t tmp_18 = p_affine_0_0*p_affine_1_1;
      real_t tmp_19 = p_affine_0_0*p_affine_1_2;
      real_t tmp_20 = p_affine_2_1*p_affine_3_2;
      real_t tmp_21 = p_affine_0_1*p_affine_1_0;
      real_t tmp_22 = p_affine_0_1*p_affine_1_2;
      real_t tmp_23 = p_affine_2_2*p_affine_3_0;
      real_t tmp_24 = p_affine_0_2*p_affine_1_0;
      real_t tmp_25 = p_affine_0_2*p_affine_1_1;
      real_t tmp_26 = p_affine_2_0*p_affine_3_1;
      real_t tmp_27 = p_affine_2_2*p_affine_3_1;
      real_t tmp_28 = p_affine_2_0*p_affine_3_2;
      real_t tmp_29 = p_affine_2_1*p_affine_3_0;
      real_t tmp_30 = std::abs(p_affine_0_0*tmp_20 - p_affine_0_0*tmp_27 + p_affine_0_1*tmp_23 - p_affine_0_1*tmp_28 + p_affine_0_2*tmp_26 - p_affine_0_2*tmp_29 - p_affine_1_0*tmp_20 + p_affine_1_0*tmp_27 - p_affine_1_1*tmp_23 + p_affine_1_1*tmp_28 - p_affine_1_2*tmp_26 + p_affine_1_2*tmp_29 + p_affine_2_0*tmp_22 - p_affine_2_0*tmp_25 - p_affine_2_1*tmp_19 + p_affine_2_1*tmp_24 + p_affine_2_2*tmp_18 - p_affine_2_2*tmp_21 - p_affine_3_0*tmp_22 + p_affine_3_0*tmp_25 + p_affine_3_1*tmp_19 - p_affine_3_1*tmp_24 - p_affine_3_2*tmp_18 + p_affine_3_2*tmp_21);
      real_t tmp_31 = tmp_30*(tmp_15 + tmp_16 + tmp_17);
      real_t tmp_32 = tmp_17*tmp_30;
      real_t tmp_33 = tmp_16*tmp_30;
      real_t tmp_34 = tmp_15*tmp_30;
      real_t a_0_0 = 0.1666666666666668*tmp_31;
      real_t a_0_1 = -0.1666666666666668*tmp_32;
      real_t a_0_2 = -0.1666666666666668*tmp_33;
      real_t a_0_3 = -0.1666666666666668*tmp_34;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }



   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                                                     const std::vector< Point3D >& coordsFacet,
                                                     const Point3D&,
                                                     const Point3D&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                               MatrixXr&                            elMat ) const override
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

         real_t tmp_0 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_2 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = tmp_0*tmp_1 - tmp_2*tmp_3;
      real_t tmp_5 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = tmp_3*tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_1*tmp_6;
      real_t tmp_13 = tmp_0*tmp_9;
      real_t tmp_14 = tmp_2*tmp_8;
      real_t tmp_15 = 1.0 / (tmp_0*tmp_1*tmp_8 - tmp_10*tmp_12 + tmp_11*tmp_2 - tmp_13*tmp_5 - tmp_14*tmp_3 + tmp_5*tmp_7);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.031405749086161582*tmp_17 + 0.93718850182767688*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_0*tmp_5 + tmp_10*tmp_2;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.031405749086161582*tmp_23 + 0.93718850182767688*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_1*tmp_10 + tmp_3*tmp_5;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_4 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = -tmp_12 + tmp_2*tmp_9;
      real_t tmp_35 = -tmp_14 + tmp_5*tmp_6;
      real_t tmp_36 = tmp_1*tmp_8 - tmp_5*tmp_9;
      real_t tmp_37 = tmp_20*tmp_34 + tmp_26*tmp_35 + tmp_32*tmp_36;
      real_t tmp_38 = -tmp_13 + tmp_7;
      real_t tmp_39 = tmp_0*tmp_8 - tmp_10*tmp_6;
      real_t tmp_40 = tmp_11 - tmp_3*tmp_8;
      real_t tmp_41 = tmp_20*tmp_38 + tmp_26*tmp_39 + tmp_32*tmp_40;
      real_t tmp_42 = 0.5*p_affine_13_0*std::pow((std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)), 1.0/2.0);
      real_t tmp_43 = 0.0068572537431980923*tmp_42;
      real_t tmp_44 = tmp_15*(0.19601935860219369*tmp_17 + 0.60796128279561268*tmp_18 + tmp_19);
      real_t tmp_45 = tmp_15*(0.19601935860219369*tmp_23 + 0.60796128279561268*tmp_24 + tmp_25);
      real_t tmp_46 = tmp_15*(0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_30 + tmp_31);
      real_t tmp_47 = tmp_21*tmp_45 + tmp_27*tmp_46 + tmp_4*tmp_44;
      real_t tmp_48 = tmp_34*tmp_44 + tmp_35*tmp_45 + tmp_36*tmp_46;
      real_t tmp_49 = tmp_38*tmp_44 + tmp_39*tmp_45 + tmp_40*tmp_46;
      real_t tmp_50 = 0.037198804536718075*tmp_42;
      real_t tmp_51 = tmp_15*(0.37605877282253791*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_52 = tmp_15*(0.37605877282253791*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_53 = tmp_15*(0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_54 = tmp_21*tmp_52 + tmp_27*tmp_53 + tmp_4*tmp_51;
      real_t tmp_55 = tmp_34*tmp_51 + tmp_35*tmp_52 + tmp_36*tmp_53;
      real_t tmp_56 = tmp_38*tmp_51 + tmp_39*tmp_52 + tmp_40*tmp_53;
      real_t tmp_57 = 0.020848748529055869*tmp_42;
      real_t tmp_58 = tmp_15*(0.78764240869137092*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_59 = tmp_15*(0.78764240869137092*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_60 = tmp_15*(0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_61 = tmp_21*tmp_59 + tmp_27*tmp_60 + tmp_4*tmp_58;
      real_t tmp_62 = tmp_34*tmp_58 + tmp_35*tmp_59 + tmp_36*tmp_60;
      real_t tmp_63 = tmp_38*tmp_58 + tmp_39*tmp_59 + tmp_40*tmp_60;
      real_t tmp_64 = 0.019202922745021479*tmp_42;
      real_t tmp_65 = tmp_15*(0.58463275527740355*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_66 = tmp_15*(0.58463275527740355*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_67 = tmp_15*(0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_68 = tmp_21*tmp_66 + tmp_27*tmp_67 + tmp_4*tmp_65;
      real_t tmp_69 = tmp_34*tmp_65 + tmp_35*tmp_66 + tmp_36*tmp_67;
      real_t tmp_70 = tmp_38*tmp_65 + tmp_39*tmp_66 + tmp_40*tmp_67;
      real_t tmp_71 = 0.020848748529055869*tmp_42;
      real_t tmp_72 = tmp_15*(0.041227165399737475*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_73 = tmp_15*(0.041227165399737475*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_74 = tmp_15*(0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_75 = tmp_21*tmp_73 + tmp_27*tmp_74 + tmp_4*tmp_72;
      real_t tmp_76 = tmp_34*tmp_72 + tmp_35*tmp_73 + tmp_36*tmp_74;
      real_t tmp_77 = tmp_38*tmp_72 + tmp_39*tmp_73 + tmp_40*tmp_74;
      real_t tmp_78 = 0.019202922745021479*tmp_42;
      real_t tmp_79 = tmp_15*(0.039308471900058539*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_80 = tmp_15*(0.039308471900058539*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_81 = tmp_15*(0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_82 = tmp_21*tmp_80 + tmp_27*tmp_81 + tmp_4*tmp_79;
      real_t tmp_83 = tmp_34*tmp_79 + tmp_35*tmp_80 + tmp_36*tmp_81;
      real_t tmp_84 = tmp_38*tmp_79 + tmp_39*tmp_80 + tmp_40*tmp_81;
      real_t tmp_85 = 0.020848748529055869*tmp_42;
      real_t tmp_86 = tmp_15*(0.78764240869137092*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_87 = tmp_15*(0.78764240869137092*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_88 = tmp_15*(0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_89 = tmp_21*tmp_87 + tmp_27*tmp_88 + tmp_4*tmp_86;
      real_t tmp_90 = tmp_34*tmp_86 + tmp_35*tmp_87 + tmp_36*tmp_88;
      real_t tmp_91 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88;
      real_t tmp_92 = 0.019202922745021479*tmp_42;
      real_t tmp_93 = tmp_15*(0.58463275527740355*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_94 = tmp_15*(0.58463275527740355*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_95 = tmp_15*(0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_96 = tmp_21*tmp_94 + tmp_27*tmp_95 + tmp_4*tmp_93;
      real_t tmp_97 = tmp_34*tmp_93 + tmp_35*tmp_94 + tmp_36*tmp_95;
      real_t tmp_98 = tmp_38*tmp_93 + tmp_39*tmp_94 + tmp_40*tmp_95;
      real_t tmp_99 = 0.020848748529055869*tmp_42;
      real_t tmp_100 = tmp_15*(0.1711304259088916*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_101 = tmp_15*(0.1711304259088916*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_102 = tmp_15*(0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_103 = tmp_100*tmp_4 + tmp_101*tmp_21 + tmp_102*tmp_27;
      real_t tmp_104 = tmp_100*tmp_34 + tmp_101*tmp_35 + tmp_102*tmp_36;
      real_t tmp_105 = tmp_100*tmp_38 + tmp_101*tmp_39 + tmp_102*tmp_40;
      real_t tmp_106 = 0.019202922745021479*tmp_42;
      real_t tmp_107 = tmp_15*(0.37605877282253791*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_108 = tmp_15*(0.37605877282253791*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_109 = tmp_15*(0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_110 = tmp_107*tmp_4 + tmp_108*tmp_21 + tmp_109*tmp_27;
      real_t tmp_111 = tmp_107*tmp_34 + tmp_108*tmp_35 + tmp_109*tmp_36;
      real_t tmp_112 = tmp_107*tmp_38 + tmp_108*tmp_39 + tmp_109*tmp_40;
      real_t tmp_113 = 0.020848748529055869*tmp_42;
      real_t tmp_114 = tmp_15*(0.041227165399737475*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_115 = tmp_15*(0.041227165399737475*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_116 = tmp_15*(0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_117 = tmp_114*tmp_4 + tmp_115*tmp_21 + tmp_116*tmp_27;
      real_t tmp_118 = tmp_114*tmp_34 + tmp_115*tmp_35 + tmp_116*tmp_36;
      real_t tmp_119 = tmp_114*tmp_38 + tmp_115*tmp_39 + tmp_116*tmp_40;
      real_t tmp_120 = 0.019202922745021479*tmp_42;
      real_t tmp_121 = tmp_15*(0.40446199974765351*tmp_17 + 0.19107600050469298*tmp_18 + tmp_19);
      real_t tmp_122 = tmp_15*(0.40446199974765351*tmp_23 + 0.19107600050469298*tmp_24 + tmp_25);
      real_t tmp_123 = tmp_15*(0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_30 + tmp_31);
      real_t tmp_124 = tmp_121*tmp_4 + tmp_122*tmp_21 + tmp_123*tmp_27;
      real_t tmp_125 = tmp_121*tmp_34 + tmp_122*tmp_35 + tmp_123*tmp_36;
      real_t tmp_126 = tmp_121*tmp_38 + tmp_122*tmp_39 + tmp_123*tmp_40;
      real_t tmp_127 = 0.042507265838595799*tmp_42;
      real_t tmp_128 = tmp_15*(0.039308471900058539*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_129 = tmp_15*(0.039308471900058539*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_130 = tmp_15*(0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_131 = tmp_128*tmp_4 + tmp_129*tmp_21 + tmp_130*tmp_27;
      real_t tmp_132 = tmp_128*tmp_34 + tmp_129*tmp_35 + tmp_130*tmp_36;
      real_t tmp_133 = tmp_128*tmp_38 + tmp_129*tmp_39 + tmp_130*tmp_40;
      real_t tmp_134 = 0.020848748529055869*tmp_42;
      real_t tmp_135 = tmp_15*(0.93718850182767688*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_136 = tmp_15*(0.93718850182767688*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_137 = tmp_15*(0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_138 = tmp_135*tmp_4 + tmp_136*tmp_21 + tmp_137*tmp_27;
      real_t tmp_139 = tmp_135*tmp_34 + tmp_136*tmp_35 + tmp_137*tmp_36;
      real_t tmp_140 = tmp_135*tmp_38 + tmp_136*tmp_39 + tmp_137*tmp_40;
      real_t tmp_141 = 0.0068572537431980923*tmp_42;
      real_t tmp_142 = tmp_15*(0.60796128279561268*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_143 = tmp_15*(0.60796128279561268*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_144 = tmp_15*(0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_145 = tmp_142*tmp_4 + tmp_143*tmp_21 + tmp_144*tmp_27;
      real_t tmp_146 = tmp_142*tmp_34 + tmp_143*tmp_35 + tmp_144*tmp_36;
      real_t tmp_147 = tmp_142*tmp_38 + tmp_143*tmp_39 + tmp_144*tmp_40;
      real_t tmp_148 = 0.037198804536718075*tmp_42;
      real_t tmp_149 = tmp_15*(0.19107600050469298*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_150 = tmp_15*(0.19107600050469298*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_151 = tmp_15*(0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_152 = tmp_149*tmp_4 + tmp_150*tmp_21 + tmp_151*tmp_27;
      real_t tmp_153 = tmp_149*tmp_34 + tmp_150*tmp_35 + tmp_151*tmp_36;
      real_t tmp_154 = tmp_149*tmp_38 + tmp_150*tmp_39 + tmp_151*tmp_40;
      real_t tmp_155 = 0.042507265838595799*tmp_42;
      real_t tmp_156 = tmp_15*(0.031405749086161582*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_157 = tmp_15*(0.031405749086161582*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_158 = tmp_15*(0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_159 = tmp_156*tmp_4 + tmp_157*tmp_21 + tmp_158*tmp_27;
      real_t tmp_160 = tmp_156*tmp_34 + tmp_157*tmp_35 + tmp_158*tmp_36;
      real_t tmp_161 = tmp_156*tmp_38 + tmp_157*tmp_39 + tmp_158*tmp_40;
      real_t tmp_162 = 0.0068572537431980923*tmp_42;
      real_t tmp_163 = tmp_15*(0.19601935860219369*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_164 = tmp_15*(0.19601935860219369*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_165 = tmp_15*(0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_166 = tmp_163*tmp_4 + tmp_164*tmp_21 + tmp_165*tmp_27;
      real_t tmp_167 = tmp_163*tmp_34 + tmp_164*tmp_35 + tmp_165*tmp_36;
      real_t tmp_168 = tmp_163*tmp_38 + tmp_164*tmp_39 + tmp_165*tmp_40;
      real_t tmp_169 = 0.037198804536718075*tmp_42;
      real_t tmp_170 = tmp_15*(0.40446199974765351*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_171 = tmp_15*(0.40446199974765351*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_172 = tmp_15*(0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_173 = tmp_170*tmp_4 + tmp_171*tmp_21 + tmp_172*tmp_27;
      real_t tmp_174 = tmp_170*tmp_34 + tmp_171*tmp_35 + tmp_172*tmp_36;
      real_t tmp_175 = tmp_170*tmp_38 + tmp_171*tmp_39 + tmp_172*tmp_40;
      real_t tmp_176 = 0.042507265838595799*tmp_42;
      real_t tmp_177 = tmp_15*(0.1711304259088916*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_178 = tmp_15*(0.1711304259088916*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_179 = tmp_15*(0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_180 = tmp_177*tmp_4 + tmp_178*tmp_21 + tmp_179*tmp_27;
      real_t tmp_181 = tmp_177*tmp_34 + tmp_178*tmp_35 + tmp_179*tmp_36;
      real_t tmp_182 = tmp_177*tmp_38 + tmp_178*tmp_39 + tmp_179*tmp_40;
      real_t tmp_183 = 0.019202922745021479*tmp_42;
      real_t a_0_0 = tmp_106*(-tmp_103 - tmp_104 - tmp_105 + 1) + tmp_113*(-tmp_110 - tmp_111 - tmp_112 + 1) + tmp_120*(-tmp_117 - tmp_118 - tmp_119 + 1) + tmp_127*(-tmp_124 - tmp_125 - tmp_126 + 1) + tmp_134*(-tmp_131 - tmp_132 - tmp_133 + 1) + tmp_141*(-tmp_138 - tmp_139 - tmp_140 + 1) + tmp_148*(-tmp_145 - tmp_146 - tmp_147 + 1) + tmp_155*(-tmp_152 - tmp_153 - tmp_154 + 1) + tmp_162*(-tmp_159 - tmp_160 - tmp_161 + 1) + tmp_169*(-tmp_166 - tmp_167 - tmp_168 + 1) + tmp_176*(-tmp_173 - tmp_174 - tmp_175 + 1) + tmp_183*(-tmp_180 - tmp_181 - tmp_182 + 1) + tmp_43*(-tmp_33 - tmp_37 - tmp_41 + 1) + tmp_50*(-tmp_47 - tmp_48 - tmp_49 + 1) + tmp_57*(-tmp_54 - tmp_55 - tmp_56 + 1) + tmp_64*(-tmp_61 - tmp_62 - tmp_63 + 1) + tmp_71*(-tmp_68 - tmp_69 - tmp_70 + 1) + tmp_78*(-tmp_75 - tmp_76 - tmp_77 + 1) + tmp_85*(-tmp_82 - tmp_83 - tmp_84 + 1) + tmp_92*(-tmp_89 - tmp_90 - tmp_91 + 1) + tmp_99*(-tmp_96 - tmp_97 - tmp_98 + 1);
      real_t a_0_1 = tmp_103*tmp_106 + tmp_110*tmp_113 + tmp_117*tmp_120 + tmp_124*tmp_127 + tmp_131*tmp_134 + tmp_138*tmp_141 + tmp_145*tmp_148 + tmp_152*tmp_155 + tmp_159*tmp_162 + tmp_166*tmp_169 + tmp_173*tmp_176 + tmp_180*tmp_183 + tmp_33*tmp_43 + tmp_47*tmp_50 + tmp_54*tmp_57 + tmp_61*tmp_64 + tmp_68*tmp_71 + tmp_75*tmp_78 + tmp_82*tmp_85 + tmp_89*tmp_92 + tmp_96*tmp_99;
      real_t a_0_2 = tmp_104*tmp_106 + tmp_111*tmp_113 + tmp_118*tmp_120 + tmp_125*tmp_127 + tmp_132*tmp_134 + tmp_139*tmp_141 + tmp_146*tmp_148 + tmp_153*tmp_155 + tmp_160*tmp_162 + tmp_167*tmp_169 + tmp_174*tmp_176 + tmp_181*tmp_183 + tmp_37*tmp_43 + tmp_48*tmp_50 + tmp_55*tmp_57 + tmp_62*tmp_64 + tmp_69*tmp_71 + tmp_76*tmp_78 + tmp_83*tmp_85 + tmp_90*tmp_92 + tmp_97*tmp_99;
      real_t a_0_3 = tmp_105*tmp_106 + tmp_112*tmp_113 + tmp_119*tmp_120 + tmp_126*tmp_127 + tmp_133*tmp_134 + tmp_140*tmp_141 + tmp_147*tmp_148 + tmp_154*tmp_155 + tmp_161*tmp_162 + tmp_168*tmp_169 + tmp_175*tmp_176 + tmp_182*tmp_183 + tmp_41*tmp_43 + tmp_49*tmp_50 + tmp_56*tmp_57 + tmp_63*tmp_64 + tmp_70*tmp_71 + tmp_77*tmp_78 + tmp_84*tmp_85 + tmp_91*tmp_92 + tmp_98*tmp_99;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }




void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                                        const std::vector< Point3D >& coordsElementOuter,
                                                        const std::vector< Point3D >& coordsFacet,
                                                        const Point3D&,
                                                        const Point3D&,
                                                        const Point3D&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                  MatrixXr&                            elMat ) const override
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElementInner[0]( 0 );
   const auto p_affine_0_1 = coordsElementInner[0]( 1 );
   const auto p_affine_0_2 = coordsElementInner[0]( 2 );

   const auto p_affine_1_0 = coordsElementInner[1]( 0 );
   const auto p_affine_1_1 = coordsElementInner[1]( 1 );
   const auto p_affine_1_2 = coordsElementInner[1]( 2 );

   const auto p_affine_2_0 = coordsElementInner[2]( 0 );
   const auto p_affine_2_1 = coordsElementInner[2]( 1 );
   const auto p_affine_2_2 = coordsElementInner[2]( 2 );

   const auto p_affine_3_0 = coordsElementInner[3]( 0 );
   const auto p_affine_3_1 = coordsElementInner[3]( 1 );
   const auto p_affine_3_2 = coordsElementInner[3]( 2 );

   const auto p_affine_4_0 = coordsElementOuter[0]( 0 );
   const auto p_affine_4_1 = coordsElementOuter[0]( 1 );
   const auto p_affine_4_2 = coordsElementOuter[0]( 2 );

   const auto p_affine_5_0 = coordsElementOuter[1]( 0 );
   const auto p_affine_5_1 = coordsElementOuter[1]( 1 );
   const auto p_affine_5_2 = coordsElementOuter[1]( 2 );

   const auto p_affine_6_0 = coordsElementOuter[2]( 0 );
   const auto p_affine_6_1 = coordsElementOuter[2]( 1 );
   const auto p_affine_6_2 = coordsElementOuter[2]( 2 );

   const auto p_affine_7_0 = coordsElementOuter[3]( 0 );
   const auto p_affine_7_1 = coordsElementOuter[3]( 1 );
   const auto p_affine_7_2 = coordsElementOuter[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t tmp_0 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_1 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_2 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_3 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_4 = tmp_0*tmp_1 - tmp_2*tmp_3;
      real_t tmp_5 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_6 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_7 = tmp_3*tmp_6;
      real_t tmp_8 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_9 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_10 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_1*tmp_6;
      real_t tmp_13 = tmp_0*tmp_9;
      real_t tmp_14 = tmp_2*tmp_8;
      real_t tmp_15 = 1.0 / (tmp_0*tmp_1*tmp_8 - tmp_10*tmp_12 + tmp_11*tmp_2 - tmp_13*tmp_5 - tmp_14*tmp_3 + tmp_5*tmp_7);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.031405749086161582*tmp_17 + 0.93718850182767688*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_0*tmp_5 + tmp_10*tmp_2;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.031405749086161582*tmp_23 + 0.93718850182767688*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_1*tmp_10 + tmp_3*tmp_5;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_4 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = -tmp_12 + tmp_2*tmp_9;
      real_t tmp_35 = -tmp_14 + tmp_5*tmp_6;
      real_t tmp_36 = tmp_1*tmp_8 - tmp_5*tmp_9;
      real_t tmp_37 = tmp_20*tmp_34 + tmp_26*tmp_35 + tmp_32*tmp_36;
      real_t tmp_38 = -tmp_13 + tmp_7;
      real_t tmp_39 = tmp_0*tmp_8 - tmp_10*tmp_6;
      real_t tmp_40 = tmp_11 - tmp_3*tmp_8;
      real_t tmp_41 = tmp_20*tmp_38 + tmp_26*tmp_39 + tmp_32*tmp_40;
      real_t tmp_42 = 0.5*p_affine_13_0*std::pow((std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)), 1.0/2.0);
      real_t tmp_43 = 0.0068572537431980923*tmp_42;
      real_t tmp_44 = tmp_15*(0.19601935860219369*tmp_17 + 0.60796128279561268*tmp_18 + tmp_19);
      real_t tmp_45 = tmp_15*(0.19601935860219369*tmp_23 + 0.60796128279561268*tmp_24 + tmp_25);
      real_t tmp_46 = tmp_15*(0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_30 + tmp_31);
      real_t tmp_47 = tmp_21*tmp_45 + tmp_27*tmp_46 + tmp_4*tmp_44;
      real_t tmp_48 = tmp_34*tmp_44 + tmp_35*tmp_45 + tmp_36*tmp_46;
      real_t tmp_49 = tmp_38*tmp_44 + tmp_39*tmp_45 + tmp_40*tmp_46;
      real_t tmp_50 = 0.037198804536718075*tmp_42;
      real_t tmp_51 = tmp_15*(0.37605877282253791*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_52 = tmp_15*(0.37605877282253791*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_53 = tmp_15*(0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_54 = tmp_21*tmp_52 + tmp_27*tmp_53 + tmp_4*tmp_51;
      real_t tmp_55 = tmp_34*tmp_51 + tmp_35*tmp_52 + tmp_36*tmp_53;
      real_t tmp_56 = tmp_38*tmp_51 + tmp_39*tmp_52 + tmp_40*tmp_53;
      real_t tmp_57 = 0.020848748529055869*tmp_42;
      real_t tmp_58 = tmp_15*(0.78764240869137092*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_59 = tmp_15*(0.78764240869137092*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_60 = tmp_15*(0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_61 = tmp_21*tmp_59 + tmp_27*tmp_60 + tmp_4*tmp_58;
      real_t tmp_62 = tmp_34*tmp_58 + tmp_35*tmp_59 + tmp_36*tmp_60;
      real_t tmp_63 = tmp_38*tmp_58 + tmp_39*tmp_59 + tmp_40*tmp_60;
      real_t tmp_64 = 0.019202922745021479*tmp_42;
      real_t tmp_65 = tmp_15*(0.58463275527740355*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_66 = tmp_15*(0.58463275527740355*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_67 = tmp_15*(0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_68 = tmp_21*tmp_66 + tmp_27*tmp_67 + tmp_4*tmp_65;
      real_t tmp_69 = tmp_34*tmp_65 + tmp_35*tmp_66 + tmp_36*tmp_67;
      real_t tmp_70 = tmp_38*tmp_65 + tmp_39*tmp_66 + tmp_40*tmp_67;
      real_t tmp_71 = 0.020848748529055869*tmp_42;
      real_t tmp_72 = tmp_15*(0.041227165399737475*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_73 = tmp_15*(0.041227165399737475*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_74 = tmp_15*(0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_75 = tmp_21*tmp_73 + tmp_27*tmp_74 + tmp_4*tmp_72;
      real_t tmp_76 = tmp_34*tmp_72 + tmp_35*tmp_73 + tmp_36*tmp_74;
      real_t tmp_77 = tmp_38*tmp_72 + tmp_39*tmp_73 + tmp_40*tmp_74;
      real_t tmp_78 = 0.019202922745021479*tmp_42;
      real_t tmp_79 = tmp_15*(0.039308471900058539*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_80 = tmp_15*(0.039308471900058539*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_81 = tmp_15*(0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_82 = tmp_21*tmp_80 + tmp_27*tmp_81 + tmp_4*tmp_79;
      real_t tmp_83 = tmp_34*tmp_79 + tmp_35*tmp_80 + tmp_36*tmp_81;
      real_t tmp_84 = tmp_38*tmp_79 + tmp_39*tmp_80 + tmp_40*tmp_81;
      real_t tmp_85 = 0.020848748529055869*tmp_42;
      real_t tmp_86 = tmp_15*(0.78764240869137092*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_87 = tmp_15*(0.78764240869137092*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_88 = tmp_15*(0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_89 = tmp_21*tmp_87 + tmp_27*tmp_88 + tmp_4*tmp_86;
      real_t tmp_90 = tmp_34*tmp_86 + tmp_35*tmp_87 + tmp_36*tmp_88;
      real_t tmp_91 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88;
      real_t tmp_92 = 0.019202922745021479*tmp_42;
      real_t tmp_93 = tmp_15*(0.58463275527740355*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_94 = tmp_15*(0.58463275527740355*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_95 = tmp_15*(0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_96 = tmp_21*tmp_94 + tmp_27*tmp_95 + tmp_4*tmp_93;
      real_t tmp_97 = tmp_34*tmp_93 + tmp_35*tmp_94 + tmp_36*tmp_95;
      real_t tmp_98 = tmp_38*tmp_93 + tmp_39*tmp_94 + tmp_40*tmp_95;
      real_t tmp_99 = 0.020848748529055869*tmp_42;
      real_t tmp_100 = tmp_15*(0.1711304259088916*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_101 = tmp_15*(0.1711304259088916*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_102 = tmp_15*(0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_103 = tmp_100*tmp_4 + tmp_101*tmp_21 + tmp_102*tmp_27;
      real_t tmp_104 = tmp_100*tmp_34 + tmp_101*tmp_35 + tmp_102*tmp_36;
      real_t tmp_105 = tmp_100*tmp_38 + tmp_101*tmp_39 + tmp_102*tmp_40;
      real_t tmp_106 = 0.019202922745021479*tmp_42;
      real_t tmp_107 = tmp_15*(0.37605877282253791*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_108 = tmp_15*(0.37605877282253791*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_109 = tmp_15*(0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_110 = tmp_107*tmp_4 + tmp_108*tmp_21 + tmp_109*tmp_27;
      real_t tmp_111 = tmp_107*tmp_34 + tmp_108*tmp_35 + tmp_109*tmp_36;
      real_t tmp_112 = tmp_107*tmp_38 + tmp_108*tmp_39 + tmp_109*tmp_40;
      real_t tmp_113 = 0.020848748529055869*tmp_42;
      real_t tmp_114 = tmp_15*(0.041227165399737475*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_115 = tmp_15*(0.041227165399737475*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_116 = tmp_15*(0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_117 = tmp_114*tmp_4 + tmp_115*tmp_21 + tmp_116*tmp_27;
      real_t tmp_118 = tmp_114*tmp_34 + tmp_115*tmp_35 + tmp_116*tmp_36;
      real_t tmp_119 = tmp_114*tmp_38 + tmp_115*tmp_39 + tmp_116*tmp_40;
      real_t tmp_120 = 0.019202922745021479*tmp_42;
      real_t tmp_121 = tmp_15*(0.40446199974765351*tmp_17 + 0.19107600050469298*tmp_18 + tmp_19);
      real_t tmp_122 = tmp_15*(0.40446199974765351*tmp_23 + 0.19107600050469298*tmp_24 + tmp_25);
      real_t tmp_123 = tmp_15*(0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_30 + tmp_31);
      real_t tmp_124 = tmp_121*tmp_4 + tmp_122*tmp_21 + tmp_123*tmp_27;
      real_t tmp_125 = tmp_121*tmp_34 + tmp_122*tmp_35 + tmp_123*tmp_36;
      real_t tmp_126 = tmp_121*tmp_38 + tmp_122*tmp_39 + tmp_123*tmp_40;
      real_t tmp_127 = 0.042507265838595799*tmp_42;
      real_t tmp_128 = tmp_15*(0.039308471900058539*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_129 = tmp_15*(0.039308471900058539*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_130 = tmp_15*(0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_131 = tmp_128*tmp_4 + tmp_129*tmp_21 + tmp_130*tmp_27;
      real_t tmp_132 = tmp_128*tmp_34 + tmp_129*tmp_35 + tmp_130*tmp_36;
      real_t tmp_133 = tmp_128*tmp_38 + tmp_129*tmp_39 + tmp_130*tmp_40;
      real_t tmp_134 = 0.020848748529055869*tmp_42;
      real_t tmp_135 = tmp_15*(0.93718850182767688*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_136 = tmp_15*(0.93718850182767688*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_137 = tmp_15*(0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_138 = tmp_135*tmp_4 + tmp_136*tmp_21 + tmp_137*tmp_27;
      real_t tmp_139 = tmp_135*tmp_34 + tmp_136*tmp_35 + tmp_137*tmp_36;
      real_t tmp_140 = tmp_135*tmp_38 + tmp_136*tmp_39 + tmp_137*tmp_40;
      real_t tmp_141 = 0.0068572537431980923*tmp_42;
      real_t tmp_142 = tmp_15*(0.60796128279561268*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_143 = tmp_15*(0.60796128279561268*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_144 = tmp_15*(0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_145 = tmp_142*tmp_4 + tmp_143*tmp_21 + tmp_144*tmp_27;
      real_t tmp_146 = tmp_142*tmp_34 + tmp_143*tmp_35 + tmp_144*tmp_36;
      real_t tmp_147 = tmp_142*tmp_38 + tmp_143*tmp_39 + tmp_144*tmp_40;
      real_t tmp_148 = 0.037198804536718075*tmp_42;
      real_t tmp_149 = tmp_15*(0.19107600050469298*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_150 = tmp_15*(0.19107600050469298*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_151 = tmp_15*(0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_152 = tmp_149*tmp_4 + tmp_150*tmp_21 + tmp_151*tmp_27;
      real_t tmp_153 = tmp_149*tmp_34 + tmp_150*tmp_35 + tmp_151*tmp_36;
      real_t tmp_154 = tmp_149*tmp_38 + tmp_150*tmp_39 + tmp_151*tmp_40;
      real_t tmp_155 = 0.042507265838595799*tmp_42;
      real_t tmp_156 = tmp_15*(0.031405749086161582*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_157 = tmp_15*(0.031405749086161582*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_158 = tmp_15*(0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_159 = tmp_156*tmp_4 + tmp_157*tmp_21 + tmp_158*tmp_27;
      real_t tmp_160 = tmp_156*tmp_34 + tmp_157*tmp_35 + tmp_158*tmp_36;
      real_t tmp_161 = tmp_156*tmp_38 + tmp_157*tmp_39 + tmp_158*tmp_40;
      real_t tmp_162 = 0.0068572537431980923*tmp_42;
      real_t tmp_163 = tmp_15*(0.19601935860219369*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_164 = tmp_15*(0.19601935860219369*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_165 = tmp_15*(0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_166 = tmp_163*tmp_4 + tmp_164*tmp_21 + tmp_165*tmp_27;
      real_t tmp_167 = tmp_163*tmp_34 + tmp_164*tmp_35 + tmp_165*tmp_36;
      real_t tmp_168 = tmp_163*tmp_38 + tmp_164*tmp_39 + tmp_165*tmp_40;
      real_t tmp_169 = 0.037198804536718075*tmp_42;
      real_t tmp_170 = tmp_15*(0.40446199974765351*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_171 = tmp_15*(0.40446199974765351*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_172 = tmp_15*(0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_173 = tmp_170*tmp_4 + tmp_171*tmp_21 + tmp_172*tmp_27;
      real_t tmp_174 = tmp_170*tmp_34 + tmp_171*tmp_35 + tmp_172*tmp_36;
      real_t tmp_175 = tmp_170*tmp_38 + tmp_171*tmp_39 + tmp_172*tmp_40;
      real_t tmp_176 = 0.042507265838595799*tmp_42;
      real_t tmp_177 = tmp_15*(0.1711304259088916*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_178 = tmp_15*(0.1711304259088916*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_179 = tmp_15*(0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_180 = tmp_177*tmp_4 + tmp_178*tmp_21 + tmp_179*tmp_27;
      real_t tmp_181 = tmp_177*tmp_34 + tmp_178*tmp_35 + tmp_179*tmp_36;
      real_t tmp_182 = tmp_177*tmp_38 + tmp_178*tmp_39 + tmp_179*tmp_40;
      real_t tmp_183 = 0.019202922745021479*tmp_42;
      real_t a_0_0 = -tmp_106*(-tmp_103 - tmp_104 - tmp_105 + 1) - tmp_113*(-tmp_110 - tmp_111 - tmp_112 + 1) - tmp_120*(-tmp_117 - tmp_118 - tmp_119 + 1) - tmp_127*(-tmp_124 - tmp_125 - tmp_126 + 1) - tmp_134*(-tmp_131 - tmp_132 - tmp_133 + 1) - tmp_141*(-tmp_138 - tmp_139 - tmp_140 + 1) - tmp_148*(-tmp_145 - tmp_146 - tmp_147 + 1) - tmp_155*(-tmp_152 - tmp_153 - tmp_154 + 1) - tmp_162*(-tmp_159 - tmp_160 - tmp_161 + 1) - tmp_169*(-tmp_166 - tmp_167 - tmp_168 + 1) - tmp_176*(-tmp_173 - tmp_174 - tmp_175 + 1) - tmp_183*(-tmp_180 - tmp_181 - tmp_182 + 1) - tmp_43*(-tmp_33 - tmp_37 - tmp_41 + 1) - tmp_50*(-tmp_47 - tmp_48 - tmp_49 + 1) - tmp_57*(-tmp_54 - tmp_55 - tmp_56 + 1) - tmp_64*(-tmp_61 - tmp_62 - tmp_63 + 1) - tmp_71*(-tmp_68 - tmp_69 - tmp_70 + 1) - tmp_78*(-tmp_75 - tmp_76 - tmp_77 + 1) - tmp_85*(-tmp_82 - tmp_83 - tmp_84 + 1) - tmp_92*(-tmp_89 - tmp_90 - tmp_91 + 1) - tmp_99*(-tmp_96 - tmp_97 - tmp_98 + 1);
      real_t a_0_1 = -tmp_103*tmp_106 - tmp_110*tmp_113 - tmp_117*tmp_120 - tmp_124*tmp_127 - tmp_131*tmp_134 - tmp_138*tmp_141 - tmp_145*tmp_148 - tmp_152*tmp_155 - tmp_159*tmp_162 - tmp_166*tmp_169 - tmp_173*tmp_176 - tmp_180*tmp_183 - tmp_33*tmp_43 - tmp_47*tmp_50 - tmp_54*tmp_57 - tmp_61*tmp_64 - tmp_68*tmp_71 - tmp_75*tmp_78 - tmp_82*tmp_85 - tmp_89*tmp_92 - tmp_96*tmp_99;
      real_t a_0_2 = -tmp_104*tmp_106 - tmp_111*tmp_113 - tmp_118*tmp_120 - tmp_125*tmp_127 - tmp_132*tmp_134 - tmp_139*tmp_141 - tmp_146*tmp_148 - tmp_153*tmp_155 - tmp_160*tmp_162 - tmp_167*tmp_169 - tmp_174*tmp_176 - tmp_181*tmp_183 - tmp_37*tmp_43 - tmp_48*tmp_50 - tmp_55*tmp_57 - tmp_62*tmp_64 - tmp_69*tmp_71 - tmp_76*tmp_78 - tmp_83*tmp_85 - tmp_90*tmp_92 - tmp_97*tmp_99;
      real_t a_0_3 = -tmp_105*tmp_106 - tmp_112*tmp_113 - tmp_119*tmp_120 - tmp_126*tmp_127 - tmp_133*tmp_134 - tmp_140*tmp_141 - tmp_147*tmp_148 - tmp_154*tmp_155 - tmp_161*tmp_162 - tmp_168*tmp_169 - tmp_175*tmp_176 - tmp_182*tmp_183 - tmp_41*tmp_43 - tmp_49*tmp_50 - tmp_56*tmp_57 - tmp_63*tmp_64 - tmp_70*tmp_71 - tmp_77*tmp_78 - tmp_84*tmp_85 - tmp_91*tmp_92 - tmp_98*tmp_99;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Point3D >& coordsElement,
    const std::vector< Point3D >& coordsFacet,
    const Point3D&,
    const Point3D&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
                                        MatrixXr&                            elMat ) const override
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t tmp_0 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_2 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = tmp_0*tmp_1 - tmp_2*tmp_3;
      real_t tmp_5 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = tmp_3*tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_1*tmp_6;
      real_t tmp_13 = tmp_0*tmp_9;
      real_t tmp_14 = tmp_2*tmp_8;
      real_t tmp_15 = 1.0 / (tmp_0*tmp_1*tmp_8 - tmp_10*tmp_12 + tmp_11*tmp_2 - tmp_13*tmp_5 - tmp_14*tmp_3 + tmp_5*tmp_7);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.031405749086161582*tmp_17 + 0.93718850182767688*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_0*tmp_5 + tmp_10*tmp_2;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.031405749086161582*tmp_23 + 0.93718850182767688*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_1*tmp_10 + tmp_3*tmp_5;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_4 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = -tmp_12 + tmp_2*tmp_9;
      real_t tmp_35 = -tmp_14 + tmp_5*tmp_6;
      real_t tmp_36 = tmp_1*tmp_8 - tmp_5*tmp_9;
      real_t tmp_37 = tmp_20*tmp_34 + tmp_26*tmp_35 + tmp_32*tmp_36;
      real_t tmp_38 = -tmp_13 + tmp_7;
      real_t tmp_39 = tmp_0*tmp_8 - tmp_10*tmp_6;
      real_t tmp_40 = tmp_11 - tmp_3*tmp_8;
      real_t tmp_41 = tmp_20*tmp_38 + tmp_26*tmp_39 + tmp_32*tmp_40;
      real_t tmp_42 = 1.0*p_affine_13_0*std::pow((std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)), 1.0/2.0);
      real_t tmp_43 = 0.0068572537431980923*tmp_42;
      real_t tmp_44 = tmp_15*(0.19601935860219369*tmp_17 + 0.60796128279561268*tmp_18 + tmp_19);
      real_t tmp_45 = tmp_15*(0.19601935860219369*tmp_23 + 0.60796128279561268*tmp_24 + tmp_25);
      real_t tmp_46 = tmp_15*(0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_30 + tmp_31);
      real_t tmp_47 = tmp_21*tmp_45 + tmp_27*tmp_46 + tmp_4*tmp_44;
      real_t tmp_48 = tmp_34*tmp_44 + tmp_35*tmp_45 + tmp_36*tmp_46;
      real_t tmp_49 = tmp_38*tmp_44 + tmp_39*tmp_45 + tmp_40*tmp_46;
      real_t tmp_50 = 0.037198804536718075*tmp_42;
      real_t tmp_51 = tmp_15*(0.37605877282253791*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_52 = tmp_15*(0.37605877282253791*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_53 = tmp_15*(0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_54 = tmp_21*tmp_52 + tmp_27*tmp_53 + tmp_4*tmp_51;
      real_t tmp_55 = tmp_34*tmp_51 + tmp_35*tmp_52 + tmp_36*tmp_53;
      real_t tmp_56 = tmp_38*tmp_51 + tmp_39*tmp_52 + tmp_40*tmp_53;
      real_t tmp_57 = 0.020848748529055869*tmp_42;
      real_t tmp_58 = tmp_15*(0.78764240869137092*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_59 = tmp_15*(0.78764240869137092*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_60 = tmp_15*(0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_61 = tmp_21*tmp_59 + tmp_27*tmp_60 + tmp_4*tmp_58;
      real_t tmp_62 = tmp_34*tmp_58 + tmp_35*tmp_59 + tmp_36*tmp_60;
      real_t tmp_63 = tmp_38*tmp_58 + tmp_39*tmp_59 + tmp_40*tmp_60;
      real_t tmp_64 = 0.019202922745021479*tmp_42;
      real_t tmp_65 = tmp_15*(0.58463275527740355*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_66 = tmp_15*(0.58463275527740355*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_67 = tmp_15*(0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_68 = tmp_21*tmp_66 + tmp_27*tmp_67 + tmp_4*tmp_65;
      real_t tmp_69 = tmp_34*tmp_65 + tmp_35*tmp_66 + tmp_36*tmp_67;
      real_t tmp_70 = tmp_38*tmp_65 + tmp_39*tmp_66 + tmp_40*tmp_67;
      real_t tmp_71 = 0.020848748529055869*tmp_42;
      real_t tmp_72 = tmp_15*(0.041227165399737475*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_73 = tmp_15*(0.041227165399737475*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_74 = tmp_15*(0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_75 = tmp_21*tmp_73 + tmp_27*tmp_74 + tmp_4*tmp_72;
      real_t tmp_76 = tmp_34*tmp_72 + tmp_35*tmp_73 + tmp_36*tmp_74;
      real_t tmp_77 = tmp_38*tmp_72 + tmp_39*tmp_73 + tmp_40*tmp_74;
      real_t tmp_78 = 0.019202922745021479*tmp_42;
      real_t tmp_79 = tmp_15*(0.039308471900058539*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_80 = tmp_15*(0.039308471900058539*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_81 = tmp_15*(0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_82 = tmp_21*tmp_80 + tmp_27*tmp_81 + tmp_4*tmp_79;
      real_t tmp_83 = tmp_34*tmp_79 + tmp_35*tmp_80 + tmp_36*tmp_81;
      real_t tmp_84 = tmp_38*tmp_79 + tmp_39*tmp_80 + tmp_40*tmp_81;
      real_t tmp_85 = 0.020848748529055869*tmp_42;
      real_t tmp_86 = tmp_15*(0.78764240869137092*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_87 = tmp_15*(0.78764240869137092*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_88 = tmp_15*(0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_89 = tmp_21*tmp_87 + tmp_27*tmp_88 + tmp_4*tmp_86;
      real_t tmp_90 = tmp_34*tmp_86 + tmp_35*tmp_87 + tmp_36*tmp_88;
      real_t tmp_91 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88;
      real_t tmp_92 = 0.019202922745021479*tmp_42;
      real_t tmp_93 = tmp_15*(0.58463275527740355*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_94 = tmp_15*(0.58463275527740355*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_95 = tmp_15*(0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_96 = tmp_21*tmp_94 + tmp_27*tmp_95 + tmp_4*tmp_93;
      real_t tmp_97 = tmp_34*tmp_93 + tmp_35*tmp_94 + tmp_36*tmp_95;
      real_t tmp_98 = tmp_38*tmp_93 + tmp_39*tmp_94 + tmp_40*tmp_95;
      real_t tmp_99 = 0.020848748529055869*tmp_42;
      real_t tmp_100 = tmp_15*(0.1711304259088916*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_101 = tmp_15*(0.1711304259088916*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_102 = tmp_15*(0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_103 = tmp_100*tmp_4 + tmp_101*tmp_21 + tmp_102*tmp_27;
      real_t tmp_104 = tmp_100*tmp_34 + tmp_101*tmp_35 + tmp_102*tmp_36;
      real_t tmp_105 = tmp_100*tmp_38 + tmp_101*tmp_39 + tmp_102*tmp_40;
      real_t tmp_106 = 0.019202922745021479*tmp_42;
      real_t tmp_107 = tmp_15*(0.37605877282253791*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_108 = tmp_15*(0.37605877282253791*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_109 = tmp_15*(0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_110 = tmp_107*tmp_4 + tmp_108*tmp_21 + tmp_109*tmp_27;
      real_t tmp_111 = tmp_107*tmp_34 + tmp_108*tmp_35 + tmp_109*tmp_36;
      real_t tmp_112 = tmp_107*tmp_38 + tmp_108*tmp_39 + tmp_109*tmp_40;
      real_t tmp_113 = 0.020848748529055869*tmp_42;
      real_t tmp_114 = tmp_15*(0.041227165399737475*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_115 = tmp_15*(0.041227165399737475*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_116 = tmp_15*(0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_117 = tmp_114*tmp_4 + tmp_115*tmp_21 + tmp_116*tmp_27;
      real_t tmp_118 = tmp_114*tmp_34 + tmp_115*tmp_35 + tmp_116*tmp_36;
      real_t tmp_119 = tmp_114*tmp_38 + tmp_115*tmp_39 + tmp_116*tmp_40;
      real_t tmp_120 = 0.019202922745021479*tmp_42;
      real_t tmp_121 = tmp_15*(0.40446199974765351*tmp_17 + 0.19107600050469298*tmp_18 + tmp_19);
      real_t tmp_122 = tmp_15*(0.40446199974765351*tmp_23 + 0.19107600050469298*tmp_24 + tmp_25);
      real_t tmp_123 = tmp_15*(0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_30 + tmp_31);
      real_t tmp_124 = tmp_121*tmp_4 + tmp_122*tmp_21 + tmp_123*tmp_27;
      real_t tmp_125 = tmp_121*tmp_34 + tmp_122*tmp_35 + tmp_123*tmp_36;
      real_t tmp_126 = tmp_121*tmp_38 + tmp_122*tmp_39 + tmp_123*tmp_40;
      real_t tmp_127 = 0.042507265838595799*tmp_42;
      real_t tmp_128 = tmp_15*(0.039308471900058539*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_129 = tmp_15*(0.039308471900058539*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_130 = tmp_15*(0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_131 = tmp_128*tmp_4 + tmp_129*tmp_21 + tmp_130*tmp_27;
      real_t tmp_132 = tmp_128*tmp_34 + tmp_129*tmp_35 + tmp_130*tmp_36;
      real_t tmp_133 = tmp_128*tmp_38 + tmp_129*tmp_39 + tmp_130*tmp_40;
      real_t tmp_134 = 0.020848748529055869*tmp_42;
      real_t tmp_135 = tmp_15*(0.93718850182767688*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_136 = tmp_15*(0.93718850182767688*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_137 = tmp_15*(0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_138 = tmp_135*tmp_4 + tmp_136*tmp_21 + tmp_137*tmp_27;
      real_t tmp_139 = tmp_135*tmp_34 + tmp_136*tmp_35 + tmp_137*tmp_36;
      real_t tmp_140 = tmp_135*tmp_38 + tmp_136*tmp_39 + tmp_137*tmp_40;
      real_t tmp_141 = 0.0068572537431980923*tmp_42;
      real_t tmp_142 = tmp_15*(0.60796128279561268*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_143 = tmp_15*(0.60796128279561268*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_144 = tmp_15*(0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_145 = tmp_142*tmp_4 + tmp_143*tmp_21 + tmp_144*tmp_27;
      real_t tmp_146 = tmp_142*tmp_34 + tmp_143*tmp_35 + tmp_144*tmp_36;
      real_t tmp_147 = tmp_142*tmp_38 + tmp_143*tmp_39 + tmp_144*tmp_40;
      real_t tmp_148 = 0.037198804536718075*tmp_42;
      real_t tmp_149 = tmp_15*(0.19107600050469298*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_150 = tmp_15*(0.19107600050469298*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_151 = tmp_15*(0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_152 = tmp_149*tmp_4 + tmp_150*tmp_21 + tmp_151*tmp_27;
      real_t tmp_153 = tmp_149*tmp_34 + tmp_150*tmp_35 + tmp_151*tmp_36;
      real_t tmp_154 = tmp_149*tmp_38 + tmp_150*tmp_39 + tmp_151*tmp_40;
      real_t tmp_155 = 0.042507265838595799*tmp_42;
      real_t tmp_156 = tmp_15*(0.031405749086161582*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_157 = tmp_15*(0.031405749086161582*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_158 = tmp_15*(0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_159 = tmp_156*tmp_4 + tmp_157*tmp_21 + tmp_158*tmp_27;
      real_t tmp_160 = tmp_156*tmp_34 + tmp_157*tmp_35 + tmp_158*tmp_36;
      real_t tmp_161 = tmp_156*tmp_38 + tmp_157*tmp_39 + tmp_158*tmp_40;
      real_t tmp_162 = 0.0068572537431980923*tmp_42;
      real_t tmp_163 = tmp_15*(0.19601935860219369*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_164 = tmp_15*(0.19601935860219369*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_165 = tmp_15*(0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_166 = tmp_163*tmp_4 + tmp_164*tmp_21 + tmp_165*tmp_27;
      real_t tmp_167 = tmp_163*tmp_34 + tmp_164*tmp_35 + tmp_165*tmp_36;
      real_t tmp_168 = tmp_163*tmp_38 + tmp_164*tmp_39 + tmp_165*tmp_40;
      real_t tmp_169 = 0.037198804536718075*tmp_42;
      real_t tmp_170 = tmp_15*(0.40446199974765351*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_171 = tmp_15*(0.40446199974765351*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_172 = tmp_15*(0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_173 = tmp_170*tmp_4 + tmp_171*tmp_21 + tmp_172*tmp_27;
      real_t tmp_174 = tmp_170*tmp_34 + tmp_171*tmp_35 + tmp_172*tmp_36;
      real_t tmp_175 = tmp_170*tmp_38 + tmp_171*tmp_39 + tmp_172*tmp_40;
      real_t tmp_176 = 0.042507265838595799*tmp_42;
      real_t tmp_177 = tmp_15*(0.1711304259088916*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_178 = tmp_15*(0.1711304259088916*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_179 = tmp_15*(0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_180 = tmp_177*tmp_4 + tmp_178*tmp_21 + tmp_179*tmp_27;
      real_t tmp_181 = tmp_177*tmp_34 + tmp_178*tmp_35 + tmp_179*tmp_36;
      real_t tmp_182 = tmp_177*tmp_38 + tmp_178*tmp_39 + tmp_179*tmp_40;
      real_t tmp_183 = 0.019202922745021479*tmp_42;
      real_t a_0_0 = tmp_106*(-tmp_103 - tmp_104 - tmp_105 + 1) + tmp_113*(-tmp_110 - tmp_111 - tmp_112 + 1) + tmp_120*(-tmp_117 - tmp_118 - tmp_119 + 1) + tmp_127*(-tmp_124 - tmp_125 - tmp_126 + 1) + tmp_134*(-tmp_131 - tmp_132 - tmp_133 + 1) + tmp_141*(-tmp_138 - tmp_139 - tmp_140 + 1) + tmp_148*(-tmp_145 - tmp_146 - tmp_147 + 1) + tmp_155*(-tmp_152 - tmp_153 - tmp_154 + 1) + tmp_162*(-tmp_159 - tmp_160 - tmp_161 + 1) + tmp_169*(-tmp_166 - tmp_167 - tmp_168 + 1) + tmp_176*(-tmp_173 - tmp_174 - tmp_175 + 1) + tmp_183*(-tmp_180 - tmp_181 - tmp_182 + 1) + tmp_43*(-tmp_33 - tmp_37 - tmp_41 + 1) + tmp_50*(-tmp_47 - tmp_48 - tmp_49 + 1) + tmp_57*(-tmp_54 - tmp_55 - tmp_56 + 1) + tmp_64*(-tmp_61 - tmp_62 - tmp_63 + 1) + tmp_71*(-tmp_68 - tmp_69 - tmp_70 + 1) + tmp_78*(-tmp_75 - tmp_76 - tmp_77 + 1) + tmp_85*(-tmp_82 - tmp_83 - tmp_84 + 1) + tmp_92*(-tmp_89 - tmp_90 - tmp_91 + 1) + tmp_99*(-tmp_96 - tmp_97 - tmp_98 + 1);
      real_t a_0_1 = tmp_103*tmp_106 + tmp_110*tmp_113 + tmp_117*tmp_120 + tmp_124*tmp_127 + tmp_131*tmp_134 + tmp_138*tmp_141 + tmp_145*tmp_148 + tmp_152*tmp_155 + tmp_159*tmp_162 + tmp_166*tmp_169 + tmp_173*tmp_176 + tmp_180*tmp_183 + tmp_33*tmp_43 + tmp_47*tmp_50 + tmp_54*tmp_57 + tmp_61*tmp_64 + tmp_68*tmp_71 + tmp_75*tmp_78 + tmp_82*tmp_85 + tmp_89*tmp_92 + tmp_96*tmp_99;
      real_t a_0_2 = tmp_104*tmp_106 + tmp_111*tmp_113 + tmp_118*tmp_120 + tmp_125*tmp_127 + tmp_132*tmp_134 + tmp_139*tmp_141 + tmp_146*tmp_148 + tmp_153*tmp_155 + tmp_160*tmp_162 + tmp_167*tmp_169 + tmp_174*tmp_176 + tmp_181*tmp_183 + tmp_37*tmp_43 + tmp_48*tmp_50 + tmp_55*tmp_57 + tmp_62*tmp_64 + tmp_69*tmp_71 + tmp_76*tmp_78 + tmp_83*tmp_85 + tmp_90*tmp_92 + tmp_97*tmp_99;
      real_t a_0_3 = tmp_105*tmp_106 + tmp_112*tmp_113 + tmp_119*tmp_120 + tmp_126*tmp_127 + tmp_133*tmp_134 + tmp_140*tmp_141 + tmp_147*tmp_148 + tmp_154*tmp_155 + tmp_161*tmp_162 + tmp_168*tmp_169 + tmp_175*tmp_176 + tmp_182*tmp_183 + tmp_41*tmp_43 + tmp_49*tmp_50 + tmp_56*tmp_57 + tmp_63*tmp_64 + tmp_70*tmp_71 + tmp_77*tmp_78 + tmp_84*tmp_85 + tmp_91*tmp_92 + tmp_98*tmp_99;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }

public:

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g1;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g1;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g0;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g0;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g2;

};




class EGDivFormNitscheBC_P0P1_1 : public hyteg::dg::DGForm
{

 public:
    EGDivFormNitscheBC_P0P1_1()
: callback_Scalar_Variable_Coefficient_3D_g1 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_2D_g1 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_2D_g0 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_3D_g0 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_3D_g2 ([](const Point3D & p) -> real_t { return 0.; })
    {}

void Scalar_Variable_Coefficient_2D_g0( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g0( Point3D( {in_0, in_1, 0} ) );
}
void Scalar_Variable_Coefficient_2D_g1( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g1( Point3D( {in_0, in_1, 0} ) );
}

void Scalar_Variable_Coefficient_3D_g0( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g0( Point3D( {in_0, in_1, in_2} ) );
}
void Scalar_Variable_Coefficient_3D_g1( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g1( Point3D( {in_0, in_1, in_2} ) );
}
void Scalar_Variable_Coefficient_3D_g2( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g2( Point3D( {in_0, in_1, in_2} ) );
}

 protected:
  void integrateVolume2D( const std::vector< Point3D >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_2 = 1.0 / (tmp_0*(-p_affine_0_1 + p_affine_2_1) + tmp_1*(-p_affine_0_1 + p_affine_1_1));
      real_t tmp_3 = tmp_0*tmp_2;
      real_t tmp_4 = tmp_1*tmp_2;
      real_t tmp_5 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_6 = tmp_5*(tmp_3 + tmp_4);
      real_t tmp_7 = tmp_4*tmp_5;
      real_t tmp_8 = tmp_3*tmp_5;
      real_t a_0_0 = 0.5*tmp_6;
      real_t a_0_1 = -0.5*tmp_7;
      real_t a_0_2 = -0.5*tmp_8;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >&      coordsElement,
                                       const std::vector< Point3D >&      coordsFacet,
                                       const Point3D&                     oppositeVertex,
                                       const Point3D&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      real_t tmp_0 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_3 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_4 = 1.0 / (-tmp_0*tmp_3 + tmp_1*tmp_2);
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_9 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_10 = tmp_4*(0.046910077030668018*tmp_8 + tmp_9);
      real_t tmp_11 = tmp_0*tmp_7 + tmp_10*tmp_2;
      real_t tmp_12 = tmp_1*tmp_7 + tmp_10*tmp_3;
      real_t tmp_13 = 0.5*p_affine_10_1*std::abs(std::pow((tmp_5*tmp_5) + (tmp_8*tmp_8), 1.0/2.0));
      real_t tmp_14 = 0.11846344252809471*tmp_13;
      real_t tmp_15 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_16 = tmp_4*(0.23076534494715845*tmp_8 + tmp_9);
      real_t tmp_17 = tmp_0*tmp_15 + tmp_16*tmp_2;
      real_t tmp_18 = tmp_1*tmp_15 + tmp_16*tmp_3;
      real_t tmp_19 = 0.2393143352496831*tmp_13;
      real_t tmp_20 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_21 = tmp_4*(0.5*tmp_8 + tmp_9);
      real_t tmp_22 = tmp_0*tmp_20 + tmp_2*tmp_21;
      real_t tmp_23 = tmp_1*tmp_20 + tmp_21*tmp_3;
      real_t tmp_24 = 0.2844444444444445*tmp_13;
      real_t tmp_25 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_26 = tmp_4*(0.7692346550528415*tmp_8 + tmp_9);
      real_t tmp_27 = tmp_0*tmp_25 + tmp_2*tmp_26;
      real_t tmp_28 = tmp_1*tmp_25 + tmp_26*tmp_3;
      real_t tmp_29 = 0.2393143352496831*tmp_13;
      real_t tmp_30 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_31 = tmp_4*(0.95308992296933193*tmp_8 + tmp_9);
      real_t tmp_32 = tmp_0*tmp_30 + tmp_2*tmp_31;
      real_t tmp_33 = tmp_1*tmp_30 + tmp_3*tmp_31;
      real_t tmp_34 = 0.11846344252809471*tmp_13;
      real_t a_0_0 = tmp_14*(-tmp_11 - tmp_12 + 1) + tmp_19*(-tmp_17 - tmp_18 + 1) + tmp_24*(-tmp_22 - tmp_23 + 1) + tmp_29*(-tmp_27 - tmp_28 + 1) + tmp_34*(-tmp_32 - tmp_33 + 1);
      real_t a_0_1 = tmp_11*tmp_14 + tmp_17*tmp_19 + tmp_22*tmp_24 + tmp_27*tmp_29 + tmp_32*tmp_34;
      real_t a_0_2 = tmp_12*tmp_14 + tmp_18*tmp_19 + tmp_23*tmp_24 + tmp_28*tmp_29 + tmp_33*tmp_34;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >&      coordsElementInner,
                                          const std::vector< Point3D >&      coordsElementOuter,
                                          const std::vector< Point3D >&      coordsFacet,
                                          const Point3D&                     oppositeVertexInnerElement,
                                          const Point3D&                     oppositeVertexOuterElement,
                                          const Point3D&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      real_t tmp_0 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_1 = -p_affine_3_0 + p_affine_4_0;
      real_t tmp_2 = -p_affine_3_1 + p_affine_5_1;
      real_t tmp_3 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_4 = 1.0 / (-tmp_0*tmp_3 + tmp_1*tmp_2);
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = -p_affine_3_1 + p_affine_6_1;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_9 = -p_affine_3_0 + p_affine_6_0;
      real_t tmp_10 = tmp_4*(0.046910077030668018*tmp_8 + tmp_9);
      real_t tmp_11 = tmp_0*tmp_7 + tmp_10*tmp_2;
      real_t tmp_12 = tmp_1*tmp_7 + tmp_10*tmp_3;
      real_t tmp_13 = 0.5*p_affine_10_1*std::abs(std::pow((tmp_5*tmp_5) + (tmp_8*tmp_8), 1.0/2.0));
      real_t tmp_14 = 0.11846344252809471*tmp_13;
      real_t tmp_15 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_16 = tmp_4*(0.23076534494715845*tmp_8 + tmp_9);
      real_t tmp_17 = tmp_0*tmp_15 + tmp_16*tmp_2;
      real_t tmp_18 = tmp_1*tmp_15 + tmp_16*tmp_3;
      real_t tmp_19 = 0.2393143352496831*tmp_13;
      real_t tmp_20 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_21 = tmp_4*(0.5*tmp_8 + tmp_9);
      real_t tmp_22 = tmp_0*tmp_20 + tmp_2*tmp_21;
      real_t tmp_23 = tmp_1*tmp_20 + tmp_21*tmp_3;
      real_t tmp_24 = 0.2844444444444445*tmp_13;
      real_t tmp_25 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_26 = tmp_4*(0.7692346550528415*tmp_8 + tmp_9);
      real_t tmp_27 = tmp_0*tmp_25 + tmp_2*tmp_26;
      real_t tmp_28 = tmp_1*tmp_25 + tmp_26*tmp_3;
      real_t tmp_29 = 0.2393143352496831*tmp_13;
      real_t tmp_30 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_31 = tmp_4*(0.95308992296933193*tmp_8 + tmp_9);
      real_t tmp_32 = tmp_0*tmp_30 + tmp_2*tmp_31;
      real_t tmp_33 = tmp_1*tmp_30 + tmp_3*tmp_31;
      real_t tmp_34 = 0.11846344252809471*tmp_13;
      real_t a_0_0 = -tmp_14*(-tmp_11 - tmp_12 + 1) - tmp_19*(-tmp_17 - tmp_18 + 1) - tmp_24*(-tmp_22 - tmp_23 + 1) - tmp_29*(-tmp_27 - tmp_28 + 1) - tmp_34*(-tmp_32 - tmp_33 + 1);
      real_t a_0_1 = -tmp_11*tmp_14 - tmp_17*tmp_19 - tmp_22*tmp_24 - tmp_27*tmp_29 - tmp_32*tmp_34;
      real_t a_0_2 = -tmp_12*tmp_14 - tmp_18*tmp_19 - tmp_23*tmp_24 - tmp_28*tmp_29 - tmp_33*tmp_34;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                   const std::vector< Point3D >&      coordsFacet,
                                                   const Point3D&                     oppositeVertex,
                                                   const Point3D&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      real_t tmp_0 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_3 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_4 = 1.0 / (-tmp_0*tmp_3 + tmp_1*tmp_2);
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_9 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_10 = tmp_4*(0.046910077030668018*tmp_8 + tmp_9);
      real_t tmp_11 = tmp_0*tmp_7 + tmp_10*tmp_2;
      real_t tmp_12 = tmp_1*tmp_7 + tmp_10*tmp_3;
      real_t tmp_13 = p_affine_10_1*std::abs(std::pow((tmp_5*tmp_5) + (tmp_8*tmp_8), 1.0/2.0));
      real_t tmp_14 = 0.11846344252809471*tmp_13;
      real_t tmp_15 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_16 = tmp_4*(0.23076534494715845*tmp_8 + tmp_9);
      real_t tmp_17 = tmp_0*tmp_15 + tmp_16*tmp_2;
      real_t tmp_18 = tmp_1*tmp_15 + tmp_16*tmp_3;
      real_t tmp_19 = 0.2393143352496831*tmp_13;
      real_t tmp_20 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_21 = tmp_4*(0.5*tmp_8 + tmp_9);
      real_t tmp_22 = tmp_0*tmp_20 + tmp_2*tmp_21;
      real_t tmp_23 = tmp_1*tmp_20 + tmp_21*tmp_3;
      real_t tmp_24 = 0.2844444444444445*tmp_13;
      real_t tmp_25 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_26 = tmp_4*(0.7692346550528415*tmp_8 + tmp_9);
      real_t tmp_27 = tmp_0*tmp_25 + tmp_2*tmp_26;
      real_t tmp_28 = tmp_1*tmp_25 + tmp_26*tmp_3;
      real_t tmp_29 = 0.2393143352496831*tmp_13;
      real_t tmp_30 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_31 = tmp_4*(0.95308992296933193*tmp_8 + tmp_9);
      real_t tmp_32 = tmp_0*tmp_30 + tmp_2*tmp_31;
      real_t tmp_33 = tmp_1*tmp_30 + tmp_3*tmp_31;
      real_t tmp_34 = 0.11846344252809471*tmp_13;
      real_t a_0_0 = tmp_14*(-tmp_11 - tmp_12 + 1) + tmp_19*(-tmp_17 - tmp_18 + 1) + tmp_24*(-tmp_22 - tmp_23 + 1) + tmp_29*(-tmp_27 - tmp_28 + 1) + tmp_34*(-tmp_32 - tmp_33 + 1);
      real_t a_0_1 = tmp_11*tmp_14 + tmp_17*tmp_19 + tmp_22*tmp_24 + tmp_27*tmp_29 + tmp_32*tmp_34;
      real_t a_0_2 = tmp_12*tmp_14 + tmp_18*tmp_19 + tmp_23*tmp_24 + tmp_28*tmp_29 + tmp_33*tmp_34;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

    void integrateRHSDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

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

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_g0_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id9 = 0;
      Scalar_Variable_Coefficient_2D_g0( 0.95308992296933193*p_affine_6_0 + 0.046910077030668018*p_affine_7_0, 0.95308992296933193*p_affine_6_1 + 0.046910077030668018*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id0 );
      Scalar_Variable_Coefficient_2D_g1( 0.95308992296933193*p_affine_6_0 + 0.046910077030668018*p_affine_7_0, 0.95308992296933193*p_affine_6_1 + 0.046910077030668018*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id1 );
      Scalar_Variable_Coefficient_2D_g0( 0.7692346550528415*p_affine_6_0 + 0.23076534494715845*p_affine_7_0, 0.7692346550528415*p_affine_6_1 + 0.23076534494715845*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id2 );
      Scalar_Variable_Coefficient_2D_g1( 0.7692346550528415*p_affine_6_0 + 0.23076534494715845*p_affine_7_0, 0.7692346550528415*p_affine_6_1 + 0.23076534494715845*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id3 );
      Scalar_Variable_Coefficient_2D_g0( 0.5*p_affine_6_0 + 0.5*p_affine_7_0, 0.5*p_affine_6_1 + 0.5*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id4 );
      Scalar_Variable_Coefficient_2D_g1( 0.5*p_affine_6_0 + 0.5*p_affine_7_0, 0.5*p_affine_6_1 + 0.5*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id5 );
      Scalar_Variable_Coefficient_2D_g0( 0.2307653449471585*p_affine_6_0 + 0.7692346550528415*p_affine_7_0, 0.2307653449471585*p_affine_6_1 + 0.7692346550528415*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id6 );
      Scalar_Variable_Coefficient_2D_g1( 0.2307653449471585*p_affine_6_0 + 0.7692346550528415*p_affine_7_0, 0.2307653449471585*p_affine_6_1 + 0.7692346550528415*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id7 );
      Scalar_Variable_Coefficient_2D_g0( 0.046910077030668074*p_affine_6_0 + 0.95308992296933193*p_affine_7_0, 0.046910077030668074*p_affine_6_1 + 0.95308992296933193*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id8 );
      Scalar_Variable_Coefficient_2D_g1( 0.046910077030668074*p_affine_6_0 + 0.95308992296933193*p_affine_7_0, 0.046910077030668074*p_affine_6_1 + 0.95308992296933193*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id9 );
      real_t tmp_0 = std::abs(std::pow(((-p_affine_6_0 + p_affine_7_0)*(-p_affine_6_0 + p_affine_7_0)) + ((-p_affine_6_1 + p_affine_7_1)*(-p_affine_6_1 + p_affine_7_1)), 1.0/2.0));
      real_t a_0_0 = 0.11846344252809471*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id0*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id1*p_affine_10_1) + 0.2393143352496831*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id2*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id3*p_affine_10_1) + 0.2844444444444445*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id4*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id5*p_affine_10_1) + 0.2393143352496831*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id6*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id7*p_affine_10_1) + 0.11846344252809471*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id8*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id9*p_affine_10_1);
      elMat( 0, 0) = a_0_0;
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

      real_t Scalar_Variable_Coefficient_3D_g0_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id13 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id14 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id15 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id16 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id17 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id18 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id19 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id20 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id21 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id22 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id23 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id24 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id25 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id26 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id27 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id28 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id29 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id30 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id31 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id32 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id33 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id34 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id35 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id36 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id37 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id38 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id39 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id40 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id41 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id42 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id43 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id44 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id45 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id46 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id47 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id48 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id49 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id50 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id51 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id52 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id53 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id54 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id55 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id56 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id57 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id58 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id59 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id60 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id61 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id62 = 0;
      Scalar_Variable_Coefficient_3D_g0( 0.93718850182767688*p_affine_10_0 + 0.031405749086161561*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.93718850182767688*p_affine_10_1 + 0.031405749086161561*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.93718850182767688*p_affine_10_2 + 0.031405749086161561*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id0 );
      Scalar_Variable_Coefficient_3D_g1( 0.93718850182767688*p_affine_10_0 + 0.031405749086161561*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.93718850182767688*p_affine_10_1 + 0.031405749086161561*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.93718850182767688*p_affine_10_2 + 0.031405749086161561*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id1 );
      Scalar_Variable_Coefficient_3D_g2( 0.93718850182767688*p_affine_10_0 + 0.031405749086161561*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.93718850182767688*p_affine_10_1 + 0.031405749086161561*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.93718850182767688*p_affine_10_2 + 0.031405749086161561*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id2 );
      Scalar_Variable_Coefficient_3D_g0( 0.60796128279561268*p_affine_10_0 + 0.19601935860219366*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.60796128279561268*p_affine_10_1 + 0.19601935860219366*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.60796128279561268*p_affine_10_2 + 0.19601935860219366*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id3 );
      Scalar_Variable_Coefficient_3D_g1( 0.60796128279561268*p_affine_10_0 + 0.19601935860219366*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.60796128279561268*p_affine_10_1 + 0.19601935860219366*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.60796128279561268*p_affine_10_2 + 0.19601935860219366*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id4 );
      Scalar_Variable_Coefficient_3D_g2( 0.60796128279561268*p_affine_10_0 + 0.19601935860219366*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.60796128279561268*p_affine_10_1 + 0.19601935860219366*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.60796128279561268*p_affine_10_2 + 0.19601935860219366*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id5 );
      Scalar_Variable_Coefficient_3D_g0( 0.039308471900058539*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id6 );
      Scalar_Variable_Coefficient_3D_g1( 0.039308471900058539*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id7 );
      Scalar_Variable_Coefficient_3D_g2( 0.039308471900058539*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id8 );
      Scalar_Variable_Coefficient_3D_g0( 0.1711304259088916*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id9 );
      Scalar_Variable_Coefficient_3D_g1( 0.1711304259088916*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id10 );
      Scalar_Variable_Coefficient_3D_g2( 0.1711304259088916*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id11 );
      Scalar_Variable_Coefficient_3D_g0( 0.37605877282253791*p_affine_10_0 + 0.039308471900058539*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.039308471900058539*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.039308471900058539*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id12 );
      Scalar_Variable_Coefficient_3D_g1( 0.37605877282253791*p_affine_10_0 + 0.039308471900058539*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.039308471900058539*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.039308471900058539*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id13 );
      Scalar_Variable_Coefficient_3D_g2( 0.37605877282253791*p_affine_10_0 + 0.039308471900058539*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.039308471900058539*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.039308471900058539*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id14 );
      Scalar_Variable_Coefficient_3D_g0( 0.78764240869137092*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id15 );
      Scalar_Variable_Coefficient_3D_g1( 0.78764240869137092*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id16 );
      Scalar_Variable_Coefficient_3D_g2( 0.78764240869137092*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id17 );
      Scalar_Variable_Coefficient_3D_g0( 0.58463275527740355*p_affine_10_0 + 0.37605877282253797*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.37605877282253797*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.37605877282253797*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id18 );
      Scalar_Variable_Coefficient_3D_g1( 0.58463275527740355*p_affine_10_0 + 0.37605877282253797*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.37605877282253797*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.37605877282253797*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id19 );
      Scalar_Variable_Coefficient_3D_g2( 0.58463275527740355*p_affine_10_0 + 0.37605877282253797*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.37605877282253797*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.37605877282253797*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id20 );
      Scalar_Variable_Coefficient_3D_g0( 0.041227165399737475*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id21 );
      Scalar_Variable_Coefficient_3D_g1( 0.041227165399737475*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id22 );
      Scalar_Variable_Coefficient_3D_g2( 0.041227165399737475*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id23 );
      Scalar_Variable_Coefficient_3D_g0( 0.039308471900058539*p_affine_10_0 + 0.37605877282253791*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.37605877282253791*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.37605877282253791*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id24 );
      Scalar_Variable_Coefficient_3D_g1( 0.039308471900058539*p_affine_10_0 + 0.37605877282253791*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.37605877282253791*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.37605877282253791*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id25 );
      Scalar_Variable_Coefficient_3D_g2( 0.039308471900058539*p_affine_10_0 + 0.37605877282253791*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.37605877282253791*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.37605877282253791*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id26 );
      Scalar_Variable_Coefficient_3D_g0( 0.78764240869137092*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id27 );
      Scalar_Variable_Coefficient_3D_g1( 0.78764240869137092*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id28 );
      Scalar_Variable_Coefficient_3D_g2( 0.78764240869137092*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id29 );
      Scalar_Variable_Coefficient_3D_g0( 0.58463275527740355*p_affine_10_0 + 0.039308471900058484*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.039308471900058484*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.039308471900058484*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id30 );
      Scalar_Variable_Coefficient_3D_g1( 0.58463275527740355*p_affine_10_0 + 0.039308471900058484*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.039308471900058484*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.039308471900058484*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id31 );
      Scalar_Variable_Coefficient_3D_g2( 0.58463275527740355*p_affine_10_0 + 0.039308471900058484*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.039308471900058484*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.039308471900058484*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id32 );
      Scalar_Variable_Coefficient_3D_g0( 0.1711304259088916*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id33 );
      Scalar_Variable_Coefficient_3D_g1( 0.1711304259088916*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id34 );
      Scalar_Variable_Coefficient_3D_g2( 0.1711304259088916*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id35 );
      Scalar_Variable_Coefficient_3D_g0( 0.19107600050469298*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.19107600050469298*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.19107600050469298*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id36 );
      Scalar_Variable_Coefficient_3D_g1( 0.19107600050469298*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.19107600050469298*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.19107600050469298*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id37 );
      Scalar_Variable_Coefficient_3D_g2( 0.19107600050469298*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.19107600050469298*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.19107600050469298*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id38 );
      Scalar_Variable_Coefficient_3D_g0( 0.37605877282253791*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id39 );
      Scalar_Variable_Coefficient_3D_g1( 0.37605877282253791*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id40 );
      Scalar_Variable_Coefficient_3D_g2( 0.37605877282253791*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id41 );
      Scalar_Variable_Coefficient_3D_g0( 0.031405749086161582*p_affine_10_0 + 0.03140574908616154*p_affine_8_0 + 0.93718850182767688*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.03140574908616154*p_affine_8_1 + 0.93718850182767688*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.03140574908616154*p_affine_8_2 + 0.93718850182767688*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id42 );
      Scalar_Variable_Coefficient_3D_g1( 0.031405749086161582*p_affine_10_0 + 0.03140574908616154*p_affine_8_0 + 0.93718850182767688*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.03140574908616154*p_affine_8_1 + 0.93718850182767688*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.03140574908616154*p_affine_8_2 + 0.93718850182767688*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id43 );
      Scalar_Variable_Coefficient_3D_g2( 0.031405749086161582*p_affine_10_0 + 0.03140574908616154*p_affine_8_0 + 0.93718850182767688*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.03140574908616154*p_affine_8_1 + 0.93718850182767688*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.03140574908616154*p_affine_8_2 + 0.93718850182767688*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id44 );
      Scalar_Variable_Coefficient_3D_g0( 0.19601935860219369*p_affine_10_0 + 0.19601935860219363*p_affine_8_0 + 0.60796128279561268*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.19601935860219363*p_affine_8_1 + 0.60796128279561268*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.19601935860219363*p_affine_8_2 + 0.60796128279561268*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id45 );
      Scalar_Variable_Coefficient_3D_g1( 0.19601935860219369*p_affine_10_0 + 0.19601935860219363*p_affine_8_0 + 0.60796128279561268*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.19601935860219363*p_affine_8_1 + 0.60796128279561268*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.19601935860219363*p_affine_8_2 + 0.60796128279561268*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id46 );
      Scalar_Variable_Coefficient_3D_g2( 0.19601935860219369*p_affine_10_0 + 0.19601935860219363*p_affine_8_0 + 0.60796128279561268*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.19601935860219363*p_affine_8_1 + 0.60796128279561268*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.19601935860219363*p_affine_8_2 + 0.60796128279561268*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id47 );
      Scalar_Variable_Coefficient_3D_g0( 0.40446199974765351*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.19107600050469298*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.19107600050469298*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.19107600050469298*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id48 );
      Scalar_Variable_Coefficient_3D_g1( 0.40446199974765351*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.19107600050469298*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.19107600050469298*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.19107600050469298*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id49 );
      Scalar_Variable_Coefficient_3D_g2( 0.40446199974765351*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.19107600050469298*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.19107600050469298*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.19107600050469298*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id50 );
      Scalar_Variable_Coefficient_3D_g0( 0.031405749086161582*p_affine_10_0 + 0.93718850182767688*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.93718850182767688*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.93718850182767688*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id51 );
      Scalar_Variable_Coefficient_3D_g1( 0.031405749086161582*p_affine_10_0 + 0.93718850182767688*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.93718850182767688*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.93718850182767688*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id52 );
      Scalar_Variable_Coefficient_3D_g2( 0.031405749086161582*p_affine_10_0 + 0.93718850182767688*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.93718850182767688*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.93718850182767688*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id53 );
      Scalar_Variable_Coefficient_3D_g0( 0.19601935860219369*p_affine_10_0 + 0.60796128279561268*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.60796128279561268*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.60796128279561268*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id54 );
      Scalar_Variable_Coefficient_3D_g1( 0.19601935860219369*p_affine_10_0 + 0.60796128279561268*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.60796128279561268*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.60796128279561268*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id55 );
      Scalar_Variable_Coefficient_3D_g2( 0.19601935860219369*p_affine_10_0 + 0.60796128279561268*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.60796128279561268*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.60796128279561268*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id56 );
      Scalar_Variable_Coefficient_3D_g0( 0.40446199974765351*p_affine_10_0 + 0.19107600050469298*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.19107600050469298*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.19107600050469298*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id57 );
      Scalar_Variable_Coefficient_3D_g1( 0.40446199974765351*p_affine_10_0 + 0.19107600050469298*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.19107600050469298*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.19107600050469298*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id58 );
      Scalar_Variable_Coefficient_3D_g2( 0.40446199974765351*p_affine_10_0 + 0.19107600050469298*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.19107600050469298*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.19107600050469298*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id59 );
      Scalar_Variable_Coefficient_3D_g0( 0.041227165399737475*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id60 );
      Scalar_Variable_Coefficient_3D_g1( 0.041227165399737475*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id61 );
      Scalar_Variable_Coefficient_3D_g2( 0.041227165399737475*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id62 );
      real_t tmp_0 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_1 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_2 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_3 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_4 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_5 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_6 = 1.0*std::pow((std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)*std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)) + (std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)*std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)) + (std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)*std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)), 1.0/2.0);
      real_t a_0_0 = 0.0068572537431980923*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id0*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id1*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id2*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id12*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id13*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id14*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id15*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id16*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id17*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id18*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id19*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id20*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id21*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id22*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id23*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id24*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id25*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id26*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id27*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id28*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id29*p_affine_13_2) + 0.037198804536718075*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id3*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id4*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id5*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id30*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id31*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id32*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id33*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id34*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id35*p_affine_13_2) + 0.042507265838595799*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id36*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id37*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id38*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id39*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id40*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id41*p_affine_13_2) + 0.0068572537431980923*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id42*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id43*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id44*p_affine_13_2) + 0.037198804536718075*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id45*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id46*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id47*p_affine_13_2) + 0.042507265838595799*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id48*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id49*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id50*p_affine_13_2) + 0.0068572537431980923*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id51*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id52*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id53*p_affine_13_2) + 0.037198804536718075*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id54*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id55*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id56*p_affine_13_2) + 0.042507265838595799*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id57*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id58*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id59*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id6*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id7*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id8*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id60*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id61*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id62*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id9*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id10*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id11*p_affine_13_2);
      elMat( 0, 0) = a_0_0;
   }
   void integrateVolume3D( const std::vector< Point3D >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coords[0]( 0 );
   const auto p_affine_0_1 = coords[0]( 1 );
   const auto p_affine_0_2 = coords[0]( 2 );

   const auto p_affine_1_0 = coords[1]( 0 );
   const auto p_affine_1_1 = coords[1]( 1 );
   const auto p_affine_1_2 = coords[1]( 2 );

   const auto p_affine_2_0 = coords[2]( 0 );
   const auto p_affine_2_1 = coords[2]( 1 );
   const auto p_affine_2_2 = coords[2]( 2 );

   const auto p_affine_3_0 = coords[3]( 0 );
   const auto p_affine_3_1 = coords[3]( 1 );
   const auto p_affine_3_2 = coords[3]( 2 );

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_4 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_5 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = tmp_0*tmp_6;
      real_t tmp_8 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_11 = tmp_3*tmp_6;
      real_t tmp_12 = tmp_10*tmp_4;
      real_t tmp_13 = 1.0 / (tmp_1*tmp_10*tmp_9 - tmp_11*tmp_9 - tmp_12*tmp_5 - tmp_2*tmp_8 + tmp_3*tmp_4*tmp_8 + tmp_5*tmp_7);
      real_t tmp_14 = tmp_13*(-tmp_2 + tmp_3*tmp_4);
      real_t tmp_15 = tmp_13*(-tmp_12 + tmp_7);
      real_t tmp_16 = tmp_13*(tmp_1*tmp_10 - tmp_11);
      real_t tmp_17 = p_affine_0_0*p_affine_1_1;
      real_t tmp_18 = p_affine_0_0*p_affine_1_2;
      real_t tmp_19 = p_affine_2_1*p_affine_3_2;
      real_t tmp_20 = p_affine_0_1*p_affine_1_0;
      real_t tmp_21 = p_affine_0_1*p_affine_1_2;
      real_t tmp_22 = p_affine_2_2*p_affine_3_0;
      real_t tmp_23 = p_affine_0_2*p_affine_1_0;
      real_t tmp_24 = p_affine_0_2*p_affine_1_1;
      real_t tmp_25 = p_affine_2_0*p_affine_3_1;
      real_t tmp_26 = p_affine_2_2*p_affine_3_1;
      real_t tmp_27 = p_affine_2_0*p_affine_3_2;
      real_t tmp_28 = p_affine_2_1*p_affine_3_0;
      real_t tmp_29 = std::abs(p_affine_0_0*tmp_19 - p_affine_0_0*tmp_26 + p_affine_0_1*tmp_22 - p_affine_0_1*tmp_27 + p_affine_0_2*tmp_25 - p_affine_0_2*tmp_28 - p_affine_1_0*tmp_19 + p_affine_1_0*tmp_26 - p_affine_1_1*tmp_22 + p_affine_1_1*tmp_27 - p_affine_1_2*tmp_25 + p_affine_1_2*tmp_28 + p_affine_2_0*tmp_21 - p_affine_2_0*tmp_24 - p_affine_2_1*tmp_18 + p_affine_2_1*tmp_23 + p_affine_2_2*tmp_17 - p_affine_2_2*tmp_20 - p_affine_3_0*tmp_21 + p_affine_3_0*tmp_24 + p_affine_3_1*tmp_18 - p_affine_3_1*tmp_23 - p_affine_3_2*tmp_17 + p_affine_3_2*tmp_20);
      real_t tmp_30 = tmp_29*(tmp_14 + tmp_15 + tmp_16);
      real_t tmp_31 = tmp_16*tmp_29;
      real_t tmp_32 = tmp_15*tmp_29;
      real_t tmp_33 = tmp_14*tmp_29;
      real_t a_0_0 = 0.1666666666666668*tmp_30;
      real_t a_0_1 = -0.1666666666666668*tmp_31;
      real_t a_0_2 = -0.1666666666666668*tmp_32;
      real_t a_0_3 = -0.1666666666666668*tmp_33;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }



   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                                                     const std::vector< Point3D >& coordsFacet,
                                                     const Point3D&,
                                                     const Point3D&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                               MatrixXr&                            elMat ) const override
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

         real_t tmp_0 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_2 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = tmp_0*tmp_1 - tmp_2*tmp_3;
      real_t tmp_5 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = tmp_3*tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_1*tmp_6;
      real_t tmp_13 = tmp_0*tmp_9;
      real_t tmp_14 = tmp_2*tmp_8;
      real_t tmp_15 = 1.0 / (tmp_0*tmp_1*tmp_8 - tmp_10*tmp_12 + tmp_11*tmp_2 - tmp_13*tmp_5 - tmp_14*tmp_3 + tmp_5*tmp_7);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.031405749086161582*tmp_17 + 0.93718850182767688*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_0*tmp_5 + tmp_10*tmp_2;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.031405749086161582*tmp_23 + 0.93718850182767688*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_1*tmp_10 + tmp_3*tmp_5;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_4 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = -tmp_12 + tmp_2*tmp_9;
      real_t tmp_35 = -tmp_14 + tmp_5*tmp_6;
      real_t tmp_36 = tmp_1*tmp_8 - tmp_5*tmp_9;
      real_t tmp_37 = tmp_20*tmp_34 + tmp_26*tmp_35 + tmp_32*tmp_36;
      real_t tmp_38 = -tmp_13 + tmp_7;
      real_t tmp_39 = tmp_0*tmp_8 - tmp_10*tmp_6;
      real_t tmp_40 = tmp_11 - tmp_3*tmp_8;
      real_t tmp_41 = tmp_20*tmp_38 + tmp_26*tmp_39 + tmp_32*tmp_40;
      real_t tmp_42 = 0.5*p_affine_13_1*std::pow((std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)), 1.0/2.0);
      real_t tmp_43 = 0.0068572537431980923*tmp_42;
      real_t tmp_44 = tmp_15*(0.19601935860219369*tmp_17 + 0.60796128279561268*tmp_18 + tmp_19);
      real_t tmp_45 = tmp_15*(0.19601935860219369*tmp_23 + 0.60796128279561268*tmp_24 + tmp_25);
      real_t tmp_46 = tmp_15*(0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_30 + tmp_31);
      real_t tmp_47 = tmp_21*tmp_45 + tmp_27*tmp_46 + tmp_4*tmp_44;
      real_t tmp_48 = tmp_34*tmp_44 + tmp_35*tmp_45 + tmp_36*tmp_46;
      real_t tmp_49 = tmp_38*tmp_44 + tmp_39*tmp_45 + tmp_40*tmp_46;
      real_t tmp_50 = 0.037198804536718075*tmp_42;
      real_t tmp_51 = tmp_15*(0.37605877282253791*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_52 = tmp_15*(0.37605877282253791*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_53 = tmp_15*(0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_54 = tmp_21*tmp_52 + tmp_27*tmp_53 + tmp_4*tmp_51;
      real_t tmp_55 = tmp_34*tmp_51 + tmp_35*tmp_52 + tmp_36*tmp_53;
      real_t tmp_56 = tmp_38*tmp_51 + tmp_39*tmp_52 + tmp_40*tmp_53;
      real_t tmp_57 = 0.020848748529055869*tmp_42;
      real_t tmp_58 = tmp_15*(0.78764240869137092*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_59 = tmp_15*(0.78764240869137092*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_60 = tmp_15*(0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_61 = tmp_21*tmp_59 + tmp_27*tmp_60 + tmp_4*tmp_58;
      real_t tmp_62 = tmp_34*tmp_58 + tmp_35*tmp_59 + tmp_36*tmp_60;
      real_t tmp_63 = tmp_38*tmp_58 + tmp_39*tmp_59 + tmp_40*tmp_60;
      real_t tmp_64 = 0.019202922745021479*tmp_42;
      real_t tmp_65 = tmp_15*(0.58463275527740355*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_66 = tmp_15*(0.58463275527740355*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_67 = tmp_15*(0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_68 = tmp_21*tmp_66 + tmp_27*tmp_67 + tmp_4*tmp_65;
      real_t tmp_69 = tmp_34*tmp_65 + tmp_35*tmp_66 + tmp_36*tmp_67;
      real_t tmp_70 = tmp_38*tmp_65 + tmp_39*tmp_66 + tmp_40*tmp_67;
      real_t tmp_71 = 0.020848748529055869*tmp_42;
      real_t tmp_72 = tmp_15*(0.041227165399737475*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_73 = tmp_15*(0.041227165399737475*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_74 = tmp_15*(0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_75 = tmp_21*tmp_73 + tmp_27*tmp_74 + tmp_4*tmp_72;
      real_t tmp_76 = tmp_34*tmp_72 + tmp_35*tmp_73 + tmp_36*tmp_74;
      real_t tmp_77 = tmp_38*tmp_72 + tmp_39*tmp_73 + tmp_40*tmp_74;
      real_t tmp_78 = 0.019202922745021479*tmp_42;
      real_t tmp_79 = tmp_15*(0.039308471900058539*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_80 = tmp_15*(0.039308471900058539*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_81 = tmp_15*(0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_82 = tmp_21*tmp_80 + tmp_27*tmp_81 + tmp_4*tmp_79;
      real_t tmp_83 = tmp_34*tmp_79 + tmp_35*tmp_80 + tmp_36*tmp_81;
      real_t tmp_84 = tmp_38*tmp_79 + tmp_39*tmp_80 + tmp_40*tmp_81;
      real_t tmp_85 = 0.020848748529055869*tmp_42;
      real_t tmp_86 = tmp_15*(0.78764240869137092*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_87 = tmp_15*(0.78764240869137092*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_88 = tmp_15*(0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_89 = tmp_21*tmp_87 + tmp_27*tmp_88 + tmp_4*tmp_86;
      real_t tmp_90 = tmp_34*tmp_86 + tmp_35*tmp_87 + tmp_36*tmp_88;
      real_t tmp_91 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88;
      real_t tmp_92 = 0.019202922745021479*tmp_42;
      real_t tmp_93 = tmp_15*(0.58463275527740355*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_94 = tmp_15*(0.58463275527740355*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_95 = tmp_15*(0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_96 = tmp_21*tmp_94 + tmp_27*tmp_95 + tmp_4*tmp_93;
      real_t tmp_97 = tmp_34*tmp_93 + tmp_35*tmp_94 + tmp_36*tmp_95;
      real_t tmp_98 = tmp_38*tmp_93 + tmp_39*tmp_94 + tmp_40*tmp_95;
      real_t tmp_99 = 0.020848748529055869*tmp_42;
      real_t tmp_100 = tmp_15*(0.1711304259088916*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_101 = tmp_15*(0.1711304259088916*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_102 = tmp_15*(0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_103 = tmp_100*tmp_4 + tmp_101*tmp_21 + tmp_102*tmp_27;
      real_t tmp_104 = tmp_100*tmp_34 + tmp_101*tmp_35 + tmp_102*tmp_36;
      real_t tmp_105 = tmp_100*tmp_38 + tmp_101*tmp_39 + tmp_102*tmp_40;
      real_t tmp_106 = 0.019202922745021479*tmp_42;
      real_t tmp_107 = tmp_15*(0.37605877282253791*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_108 = tmp_15*(0.37605877282253791*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_109 = tmp_15*(0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_110 = tmp_107*tmp_4 + tmp_108*tmp_21 + tmp_109*tmp_27;
      real_t tmp_111 = tmp_107*tmp_34 + tmp_108*tmp_35 + tmp_109*tmp_36;
      real_t tmp_112 = tmp_107*tmp_38 + tmp_108*tmp_39 + tmp_109*tmp_40;
      real_t tmp_113 = 0.020848748529055869*tmp_42;
      real_t tmp_114 = tmp_15*(0.041227165399737475*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_115 = tmp_15*(0.041227165399737475*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_116 = tmp_15*(0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_117 = tmp_114*tmp_4 + tmp_115*tmp_21 + tmp_116*tmp_27;
      real_t tmp_118 = tmp_114*tmp_34 + tmp_115*tmp_35 + tmp_116*tmp_36;
      real_t tmp_119 = tmp_114*tmp_38 + tmp_115*tmp_39 + tmp_116*tmp_40;
      real_t tmp_120 = 0.019202922745021479*tmp_42;
      real_t tmp_121 = tmp_15*(0.40446199974765351*tmp_17 + 0.19107600050469298*tmp_18 + tmp_19);
      real_t tmp_122 = tmp_15*(0.40446199974765351*tmp_23 + 0.19107600050469298*tmp_24 + tmp_25);
      real_t tmp_123 = tmp_15*(0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_30 + tmp_31);
      real_t tmp_124 = tmp_121*tmp_4 + tmp_122*tmp_21 + tmp_123*tmp_27;
      real_t tmp_125 = tmp_121*tmp_34 + tmp_122*tmp_35 + tmp_123*tmp_36;
      real_t tmp_126 = tmp_121*tmp_38 + tmp_122*tmp_39 + tmp_123*tmp_40;
      real_t tmp_127 = 0.042507265838595799*tmp_42;
      real_t tmp_128 = tmp_15*(0.039308471900058539*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_129 = tmp_15*(0.039308471900058539*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_130 = tmp_15*(0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_131 = tmp_128*tmp_4 + tmp_129*tmp_21 + tmp_130*tmp_27;
      real_t tmp_132 = tmp_128*tmp_34 + tmp_129*tmp_35 + tmp_130*tmp_36;
      real_t tmp_133 = tmp_128*tmp_38 + tmp_129*tmp_39 + tmp_130*tmp_40;
      real_t tmp_134 = 0.020848748529055869*tmp_42;
      real_t tmp_135 = tmp_15*(0.93718850182767688*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_136 = tmp_15*(0.93718850182767688*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_137 = tmp_15*(0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_138 = tmp_135*tmp_4 + tmp_136*tmp_21 + tmp_137*tmp_27;
      real_t tmp_139 = tmp_135*tmp_34 + tmp_136*tmp_35 + tmp_137*tmp_36;
      real_t tmp_140 = tmp_135*tmp_38 + tmp_136*tmp_39 + tmp_137*tmp_40;
      real_t tmp_141 = 0.0068572537431980923*tmp_42;
      real_t tmp_142 = tmp_15*(0.60796128279561268*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_143 = tmp_15*(0.60796128279561268*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_144 = tmp_15*(0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_145 = tmp_142*tmp_4 + tmp_143*tmp_21 + tmp_144*tmp_27;
      real_t tmp_146 = tmp_142*tmp_34 + tmp_143*tmp_35 + tmp_144*tmp_36;
      real_t tmp_147 = tmp_142*tmp_38 + tmp_143*tmp_39 + tmp_144*tmp_40;
      real_t tmp_148 = 0.037198804536718075*tmp_42;
      real_t tmp_149 = tmp_15*(0.19107600050469298*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_150 = tmp_15*(0.19107600050469298*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_151 = tmp_15*(0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_152 = tmp_149*tmp_4 + tmp_150*tmp_21 + tmp_151*tmp_27;
      real_t tmp_153 = tmp_149*tmp_34 + tmp_150*tmp_35 + tmp_151*tmp_36;
      real_t tmp_154 = tmp_149*tmp_38 + tmp_150*tmp_39 + tmp_151*tmp_40;
      real_t tmp_155 = 0.042507265838595799*tmp_42;
      real_t tmp_156 = tmp_15*(0.031405749086161582*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_157 = tmp_15*(0.031405749086161582*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_158 = tmp_15*(0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_159 = tmp_156*tmp_4 + tmp_157*tmp_21 + tmp_158*tmp_27;
      real_t tmp_160 = tmp_156*tmp_34 + tmp_157*tmp_35 + tmp_158*tmp_36;
      real_t tmp_161 = tmp_156*tmp_38 + tmp_157*tmp_39 + tmp_158*tmp_40;
      real_t tmp_162 = 0.0068572537431980923*tmp_42;
      real_t tmp_163 = tmp_15*(0.19601935860219369*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_164 = tmp_15*(0.19601935860219369*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_165 = tmp_15*(0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_166 = tmp_163*tmp_4 + tmp_164*tmp_21 + tmp_165*tmp_27;
      real_t tmp_167 = tmp_163*tmp_34 + tmp_164*tmp_35 + tmp_165*tmp_36;
      real_t tmp_168 = tmp_163*tmp_38 + tmp_164*tmp_39 + tmp_165*tmp_40;
      real_t tmp_169 = 0.037198804536718075*tmp_42;
      real_t tmp_170 = tmp_15*(0.40446199974765351*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_171 = tmp_15*(0.40446199974765351*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_172 = tmp_15*(0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_173 = tmp_170*tmp_4 + tmp_171*tmp_21 + tmp_172*tmp_27;
      real_t tmp_174 = tmp_170*tmp_34 + tmp_171*tmp_35 + tmp_172*tmp_36;
      real_t tmp_175 = tmp_170*tmp_38 + tmp_171*tmp_39 + tmp_172*tmp_40;
      real_t tmp_176 = 0.042507265838595799*tmp_42;
      real_t tmp_177 = tmp_15*(0.1711304259088916*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_178 = tmp_15*(0.1711304259088916*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_179 = tmp_15*(0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_180 = tmp_177*tmp_4 + tmp_178*tmp_21 + tmp_179*tmp_27;
      real_t tmp_181 = tmp_177*tmp_34 + tmp_178*tmp_35 + tmp_179*tmp_36;
      real_t tmp_182 = tmp_177*tmp_38 + tmp_178*tmp_39 + tmp_179*tmp_40;
      real_t tmp_183 = 0.019202922745021479*tmp_42;
      real_t a_0_0 = tmp_106*(-tmp_103 - tmp_104 - tmp_105 + 1) + tmp_113*(-tmp_110 - tmp_111 - tmp_112 + 1) + tmp_120*(-tmp_117 - tmp_118 - tmp_119 + 1) + tmp_127*(-tmp_124 - tmp_125 - tmp_126 + 1) + tmp_134*(-tmp_131 - tmp_132 - tmp_133 + 1) + tmp_141*(-tmp_138 - tmp_139 - tmp_140 + 1) + tmp_148*(-tmp_145 - tmp_146 - tmp_147 + 1) + tmp_155*(-tmp_152 - tmp_153 - tmp_154 + 1) + tmp_162*(-tmp_159 - tmp_160 - tmp_161 + 1) + tmp_169*(-tmp_166 - tmp_167 - tmp_168 + 1) + tmp_176*(-tmp_173 - tmp_174 - tmp_175 + 1) + tmp_183*(-tmp_180 - tmp_181 - tmp_182 + 1) + tmp_43*(-tmp_33 - tmp_37 - tmp_41 + 1) + tmp_50*(-tmp_47 - tmp_48 - tmp_49 + 1) + tmp_57*(-tmp_54 - tmp_55 - tmp_56 + 1) + tmp_64*(-tmp_61 - tmp_62 - tmp_63 + 1) + tmp_71*(-tmp_68 - tmp_69 - tmp_70 + 1) + tmp_78*(-tmp_75 - tmp_76 - tmp_77 + 1) + tmp_85*(-tmp_82 - tmp_83 - tmp_84 + 1) + tmp_92*(-tmp_89 - tmp_90 - tmp_91 + 1) + tmp_99*(-tmp_96 - tmp_97 - tmp_98 + 1);
      real_t a_0_1 = tmp_103*tmp_106 + tmp_110*tmp_113 + tmp_117*tmp_120 + tmp_124*tmp_127 + tmp_131*tmp_134 + tmp_138*tmp_141 + tmp_145*tmp_148 + tmp_152*tmp_155 + tmp_159*tmp_162 + tmp_166*tmp_169 + tmp_173*tmp_176 + tmp_180*tmp_183 + tmp_33*tmp_43 + tmp_47*tmp_50 + tmp_54*tmp_57 + tmp_61*tmp_64 + tmp_68*tmp_71 + tmp_75*tmp_78 + tmp_82*tmp_85 + tmp_89*tmp_92 + tmp_96*tmp_99;
      real_t a_0_2 = tmp_104*tmp_106 + tmp_111*tmp_113 + tmp_118*tmp_120 + tmp_125*tmp_127 + tmp_132*tmp_134 + tmp_139*tmp_141 + tmp_146*tmp_148 + tmp_153*tmp_155 + tmp_160*tmp_162 + tmp_167*tmp_169 + tmp_174*tmp_176 + tmp_181*tmp_183 + tmp_37*tmp_43 + tmp_48*tmp_50 + tmp_55*tmp_57 + tmp_62*tmp_64 + tmp_69*tmp_71 + tmp_76*tmp_78 + tmp_83*tmp_85 + tmp_90*tmp_92 + tmp_97*tmp_99;
      real_t a_0_3 = tmp_105*tmp_106 + tmp_112*tmp_113 + tmp_119*tmp_120 + tmp_126*tmp_127 + tmp_133*tmp_134 + tmp_140*tmp_141 + tmp_147*tmp_148 + tmp_154*tmp_155 + tmp_161*tmp_162 + tmp_168*tmp_169 + tmp_175*tmp_176 + tmp_182*tmp_183 + tmp_41*tmp_43 + tmp_49*tmp_50 + tmp_56*tmp_57 + tmp_63*tmp_64 + tmp_70*tmp_71 + tmp_77*tmp_78 + tmp_84*tmp_85 + tmp_91*tmp_92 + tmp_98*tmp_99;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }




void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                                        const std::vector< Point3D >& coordsElementOuter,
                                                        const std::vector< Point3D >& coordsFacet,
                                                        const Point3D&,
                                                        const Point3D&,
                                                        const Point3D&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                  MatrixXr&                            elMat ) const override
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElementInner[0]( 0 );
   const auto p_affine_0_1 = coordsElementInner[0]( 1 );
   const auto p_affine_0_2 = coordsElementInner[0]( 2 );

   const auto p_affine_1_0 = coordsElementInner[1]( 0 );
   const auto p_affine_1_1 = coordsElementInner[1]( 1 );
   const auto p_affine_1_2 = coordsElementInner[1]( 2 );

   const auto p_affine_2_0 = coordsElementInner[2]( 0 );
   const auto p_affine_2_1 = coordsElementInner[2]( 1 );
   const auto p_affine_2_2 = coordsElementInner[2]( 2 );

   const auto p_affine_3_0 = coordsElementInner[3]( 0 );
   const auto p_affine_3_1 = coordsElementInner[3]( 1 );
   const auto p_affine_3_2 = coordsElementInner[3]( 2 );

   const auto p_affine_4_0 = coordsElementOuter[0]( 0 );
   const auto p_affine_4_1 = coordsElementOuter[0]( 1 );
   const auto p_affine_4_2 = coordsElementOuter[0]( 2 );

   const auto p_affine_5_0 = coordsElementOuter[1]( 0 );
   const auto p_affine_5_1 = coordsElementOuter[1]( 1 );
   const auto p_affine_5_2 = coordsElementOuter[1]( 2 );

   const auto p_affine_6_0 = coordsElementOuter[2]( 0 );
   const auto p_affine_6_1 = coordsElementOuter[2]( 1 );
   const auto p_affine_6_2 = coordsElementOuter[2]( 2 );

   const auto p_affine_7_0 = coordsElementOuter[3]( 0 );
   const auto p_affine_7_1 = coordsElementOuter[3]( 1 );
   const auto p_affine_7_2 = coordsElementOuter[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t tmp_0 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_1 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_2 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_3 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_4 = tmp_0*tmp_1 - tmp_2*tmp_3;
      real_t tmp_5 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_6 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_7 = tmp_3*tmp_6;
      real_t tmp_8 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_9 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_10 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_1*tmp_6;
      real_t tmp_13 = tmp_0*tmp_9;
      real_t tmp_14 = tmp_2*tmp_8;
      real_t tmp_15 = 1.0 / (tmp_0*tmp_1*tmp_8 - tmp_10*tmp_12 + tmp_11*tmp_2 - tmp_13*tmp_5 - tmp_14*tmp_3 + tmp_5*tmp_7);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.031405749086161582*tmp_17 + 0.93718850182767688*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_0*tmp_5 + tmp_10*tmp_2;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.031405749086161582*tmp_23 + 0.93718850182767688*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_1*tmp_10 + tmp_3*tmp_5;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_4 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = -tmp_12 + tmp_2*tmp_9;
      real_t tmp_35 = -tmp_14 + tmp_5*tmp_6;
      real_t tmp_36 = tmp_1*tmp_8 - tmp_5*tmp_9;
      real_t tmp_37 = tmp_20*tmp_34 + tmp_26*tmp_35 + tmp_32*tmp_36;
      real_t tmp_38 = -tmp_13 + tmp_7;
      real_t tmp_39 = tmp_0*tmp_8 - tmp_10*tmp_6;
      real_t tmp_40 = tmp_11 - tmp_3*tmp_8;
      real_t tmp_41 = tmp_20*tmp_38 + tmp_26*tmp_39 + tmp_32*tmp_40;
      real_t tmp_42 = 0.5*p_affine_13_1*std::pow((std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)), 1.0/2.0);
      real_t tmp_43 = 0.0068572537431980923*tmp_42;
      real_t tmp_44 = tmp_15*(0.19601935860219369*tmp_17 + 0.60796128279561268*tmp_18 + tmp_19);
      real_t tmp_45 = tmp_15*(0.19601935860219369*tmp_23 + 0.60796128279561268*tmp_24 + tmp_25);
      real_t tmp_46 = tmp_15*(0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_30 + tmp_31);
      real_t tmp_47 = tmp_21*tmp_45 + tmp_27*tmp_46 + tmp_4*tmp_44;
      real_t tmp_48 = tmp_34*tmp_44 + tmp_35*tmp_45 + tmp_36*tmp_46;
      real_t tmp_49 = tmp_38*tmp_44 + tmp_39*tmp_45 + tmp_40*tmp_46;
      real_t tmp_50 = 0.037198804536718075*tmp_42;
      real_t tmp_51 = tmp_15*(0.37605877282253791*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_52 = tmp_15*(0.37605877282253791*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_53 = tmp_15*(0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_54 = tmp_21*tmp_52 + tmp_27*tmp_53 + tmp_4*tmp_51;
      real_t tmp_55 = tmp_34*tmp_51 + tmp_35*tmp_52 + tmp_36*tmp_53;
      real_t tmp_56 = tmp_38*tmp_51 + tmp_39*tmp_52 + tmp_40*tmp_53;
      real_t tmp_57 = 0.020848748529055869*tmp_42;
      real_t tmp_58 = tmp_15*(0.78764240869137092*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_59 = tmp_15*(0.78764240869137092*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_60 = tmp_15*(0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_61 = tmp_21*tmp_59 + tmp_27*tmp_60 + tmp_4*tmp_58;
      real_t tmp_62 = tmp_34*tmp_58 + tmp_35*tmp_59 + tmp_36*tmp_60;
      real_t tmp_63 = tmp_38*tmp_58 + tmp_39*tmp_59 + tmp_40*tmp_60;
      real_t tmp_64 = 0.019202922745021479*tmp_42;
      real_t tmp_65 = tmp_15*(0.58463275527740355*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_66 = tmp_15*(0.58463275527740355*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_67 = tmp_15*(0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_68 = tmp_21*tmp_66 + tmp_27*tmp_67 + tmp_4*tmp_65;
      real_t tmp_69 = tmp_34*tmp_65 + tmp_35*tmp_66 + tmp_36*tmp_67;
      real_t tmp_70 = tmp_38*tmp_65 + tmp_39*tmp_66 + tmp_40*tmp_67;
      real_t tmp_71 = 0.020848748529055869*tmp_42;
      real_t tmp_72 = tmp_15*(0.041227165399737475*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_73 = tmp_15*(0.041227165399737475*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_74 = tmp_15*(0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_75 = tmp_21*tmp_73 + tmp_27*tmp_74 + tmp_4*tmp_72;
      real_t tmp_76 = tmp_34*tmp_72 + tmp_35*tmp_73 + tmp_36*tmp_74;
      real_t tmp_77 = tmp_38*tmp_72 + tmp_39*tmp_73 + tmp_40*tmp_74;
      real_t tmp_78 = 0.019202922745021479*tmp_42;
      real_t tmp_79 = tmp_15*(0.039308471900058539*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_80 = tmp_15*(0.039308471900058539*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_81 = tmp_15*(0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_82 = tmp_21*tmp_80 + tmp_27*tmp_81 + tmp_4*tmp_79;
      real_t tmp_83 = tmp_34*tmp_79 + tmp_35*tmp_80 + tmp_36*tmp_81;
      real_t tmp_84 = tmp_38*tmp_79 + tmp_39*tmp_80 + tmp_40*tmp_81;
      real_t tmp_85 = 0.020848748529055869*tmp_42;
      real_t tmp_86 = tmp_15*(0.78764240869137092*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_87 = tmp_15*(0.78764240869137092*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_88 = tmp_15*(0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_89 = tmp_21*tmp_87 + tmp_27*tmp_88 + tmp_4*tmp_86;
      real_t tmp_90 = tmp_34*tmp_86 + tmp_35*tmp_87 + tmp_36*tmp_88;
      real_t tmp_91 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88;
      real_t tmp_92 = 0.019202922745021479*tmp_42;
      real_t tmp_93 = tmp_15*(0.58463275527740355*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_94 = tmp_15*(0.58463275527740355*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_95 = tmp_15*(0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_96 = tmp_21*tmp_94 + tmp_27*tmp_95 + tmp_4*tmp_93;
      real_t tmp_97 = tmp_34*tmp_93 + tmp_35*tmp_94 + tmp_36*tmp_95;
      real_t tmp_98 = tmp_38*tmp_93 + tmp_39*tmp_94 + tmp_40*tmp_95;
      real_t tmp_99 = 0.020848748529055869*tmp_42;
      real_t tmp_100 = tmp_15*(0.1711304259088916*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_101 = tmp_15*(0.1711304259088916*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_102 = tmp_15*(0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_103 = tmp_100*tmp_4 + tmp_101*tmp_21 + tmp_102*tmp_27;
      real_t tmp_104 = tmp_100*tmp_34 + tmp_101*tmp_35 + tmp_102*tmp_36;
      real_t tmp_105 = tmp_100*tmp_38 + tmp_101*tmp_39 + tmp_102*tmp_40;
      real_t tmp_106 = 0.019202922745021479*tmp_42;
      real_t tmp_107 = tmp_15*(0.37605877282253791*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_108 = tmp_15*(0.37605877282253791*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_109 = tmp_15*(0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_110 = tmp_107*tmp_4 + tmp_108*tmp_21 + tmp_109*tmp_27;
      real_t tmp_111 = tmp_107*tmp_34 + tmp_108*tmp_35 + tmp_109*tmp_36;
      real_t tmp_112 = tmp_107*tmp_38 + tmp_108*tmp_39 + tmp_109*tmp_40;
      real_t tmp_113 = 0.020848748529055869*tmp_42;
      real_t tmp_114 = tmp_15*(0.041227165399737475*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_115 = tmp_15*(0.041227165399737475*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_116 = tmp_15*(0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_117 = tmp_114*tmp_4 + tmp_115*tmp_21 + tmp_116*tmp_27;
      real_t tmp_118 = tmp_114*tmp_34 + tmp_115*tmp_35 + tmp_116*tmp_36;
      real_t tmp_119 = tmp_114*tmp_38 + tmp_115*tmp_39 + tmp_116*tmp_40;
      real_t tmp_120 = 0.019202922745021479*tmp_42;
      real_t tmp_121 = tmp_15*(0.40446199974765351*tmp_17 + 0.19107600050469298*tmp_18 + tmp_19);
      real_t tmp_122 = tmp_15*(0.40446199974765351*tmp_23 + 0.19107600050469298*tmp_24 + tmp_25);
      real_t tmp_123 = tmp_15*(0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_30 + tmp_31);
      real_t tmp_124 = tmp_121*tmp_4 + tmp_122*tmp_21 + tmp_123*tmp_27;
      real_t tmp_125 = tmp_121*tmp_34 + tmp_122*tmp_35 + tmp_123*tmp_36;
      real_t tmp_126 = tmp_121*tmp_38 + tmp_122*tmp_39 + tmp_123*tmp_40;
      real_t tmp_127 = 0.042507265838595799*tmp_42;
      real_t tmp_128 = tmp_15*(0.039308471900058539*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_129 = tmp_15*(0.039308471900058539*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_130 = tmp_15*(0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_131 = tmp_128*tmp_4 + tmp_129*tmp_21 + tmp_130*tmp_27;
      real_t tmp_132 = tmp_128*tmp_34 + tmp_129*tmp_35 + tmp_130*tmp_36;
      real_t tmp_133 = tmp_128*tmp_38 + tmp_129*tmp_39 + tmp_130*tmp_40;
      real_t tmp_134 = 0.020848748529055869*tmp_42;
      real_t tmp_135 = tmp_15*(0.93718850182767688*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_136 = tmp_15*(0.93718850182767688*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_137 = tmp_15*(0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_138 = tmp_135*tmp_4 + tmp_136*tmp_21 + tmp_137*tmp_27;
      real_t tmp_139 = tmp_135*tmp_34 + tmp_136*tmp_35 + tmp_137*tmp_36;
      real_t tmp_140 = tmp_135*tmp_38 + tmp_136*tmp_39 + tmp_137*tmp_40;
      real_t tmp_141 = 0.0068572537431980923*tmp_42;
      real_t tmp_142 = tmp_15*(0.60796128279561268*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_143 = tmp_15*(0.60796128279561268*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_144 = tmp_15*(0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_145 = tmp_142*tmp_4 + tmp_143*tmp_21 + tmp_144*tmp_27;
      real_t tmp_146 = tmp_142*tmp_34 + tmp_143*tmp_35 + tmp_144*tmp_36;
      real_t tmp_147 = tmp_142*tmp_38 + tmp_143*tmp_39 + tmp_144*tmp_40;
      real_t tmp_148 = 0.037198804536718075*tmp_42;
      real_t tmp_149 = tmp_15*(0.19107600050469298*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_150 = tmp_15*(0.19107600050469298*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_151 = tmp_15*(0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_152 = tmp_149*tmp_4 + tmp_150*tmp_21 + tmp_151*tmp_27;
      real_t tmp_153 = tmp_149*tmp_34 + tmp_150*tmp_35 + tmp_151*tmp_36;
      real_t tmp_154 = tmp_149*tmp_38 + tmp_150*tmp_39 + tmp_151*tmp_40;
      real_t tmp_155 = 0.042507265838595799*tmp_42;
      real_t tmp_156 = tmp_15*(0.031405749086161582*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_157 = tmp_15*(0.031405749086161582*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_158 = tmp_15*(0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_159 = tmp_156*tmp_4 + tmp_157*tmp_21 + tmp_158*tmp_27;
      real_t tmp_160 = tmp_156*tmp_34 + tmp_157*tmp_35 + tmp_158*tmp_36;
      real_t tmp_161 = tmp_156*tmp_38 + tmp_157*tmp_39 + tmp_158*tmp_40;
      real_t tmp_162 = 0.0068572537431980923*tmp_42;
      real_t tmp_163 = tmp_15*(0.19601935860219369*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_164 = tmp_15*(0.19601935860219369*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_165 = tmp_15*(0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_166 = tmp_163*tmp_4 + tmp_164*tmp_21 + tmp_165*tmp_27;
      real_t tmp_167 = tmp_163*tmp_34 + tmp_164*tmp_35 + tmp_165*tmp_36;
      real_t tmp_168 = tmp_163*tmp_38 + tmp_164*tmp_39 + tmp_165*tmp_40;
      real_t tmp_169 = 0.037198804536718075*tmp_42;
      real_t tmp_170 = tmp_15*(0.40446199974765351*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_171 = tmp_15*(0.40446199974765351*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_172 = tmp_15*(0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_173 = tmp_170*tmp_4 + tmp_171*tmp_21 + tmp_172*tmp_27;
      real_t tmp_174 = tmp_170*tmp_34 + tmp_171*tmp_35 + tmp_172*tmp_36;
      real_t tmp_175 = tmp_170*tmp_38 + tmp_171*tmp_39 + tmp_172*tmp_40;
      real_t tmp_176 = 0.042507265838595799*tmp_42;
      real_t tmp_177 = tmp_15*(0.1711304259088916*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_178 = tmp_15*(0.1711304259088916*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_179 = tmp_15*(0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_180 = tmp_177*tmp_4 + tmp_178*tmp_21 + tmp_179*tmp_27;
      real_t tmp_181 = tmp_177*tmp_34 + tmp_178*tmp_35 + tmp_179*tmp_36;
      real_t tmp_182 = tmp_177*tmp_38 + tmp_178*tmp_39 + tmp_179*tmp_40;
      real_t tmp_183 = 0.019202922745021479*tmp_42;
      real_t a_0_0 = -tmp_106*(-tmp_103 - tmp_104 - tmp_105 + 1) - tmp_113*(-tmp_110 - tmp_111 - tmp_112 + 1) - tmp_120*(-tmp_117 - tmp_118 - tmp_119 + 1) - tmp_127*(-tmp_124 - tmp_125 - tmp_126 + 1) - tmp_134*(-tmp_131 - tmp_132 - tmp_133 + 1) - tmp_141*(-tmp_138 - tmp_139 - tmp_140 + 1) - tmp_148*(-tmp_145 - tmp_146 - tmp_147 + 1) - tmp_155*(-tmp_152 - tmp_153 - tmp_154 + 1) - tmp_162*(-tmp_159 - tmp_160 - tmp_161 + 1) - tmp_169*(-tmp_166 - tmp_167 - tmp_168 + 1) - tmp_176*(-tmp_173 - tmp_174 - tmp_175 + 1) - tmp_183*(-tmp_180 - tmp_181 - tmp_182 + 1) - tmp_43*(-tmp_33 - tmp_37 - tmp_41 + 1) - tmp_50*(-tmp_47 - tmp_48 - tmp_49 + 1) - tmp_57*(-tmp_54 - tmp_55 - tmp_56 + 1) - tmp_64*(-tmp_61 - tmp_62 - tmp_63 + 1) - tmp_71*(-tmp_68 - tmp_69 - tmp_70 + 1) - tmp_78*(-tmp_75 - tmp_76 - tmp_77 + 1) - tmp_85*(-tmp_82 - tmp_83 - tmp_84 + 1) - tmp_92*(-tmp_89 - tmp_90 - tmp_91 + 1) - tmp_99*(-tmp_96 - tmp_97 - tmp_98 + 1);
      real_t a_0_1 = -tmp_103*tmp_106 - tmp_110*tmp_113 - tmp_117*tmp_120 - tmp_124*tmp_127 - tmp_131*tmp_134 - tmp_138*tmp_141 - tmp_145*tmp_148 - tmp_152*tmp_155 - tmp_159*tmp_162 - tmp_166*tmp_169 - tmp_173*tmp_176 - tmp_180*tmp_183 - tmp_33*tmp_43 - tmp_47*tmp_50 - tmp_54*tmp_57 - tmp_61*tmp_64 - tmp_68*tmp_71 - tmp_75*tmp_78 - tmp_82*tmp_85 - tmp_89*tmp_92 - tmp_96*tmp_99;
      real_t a_0_2 = -tmp_104*tmp_106 - tmp_111*tmp_113 - tmp_118*tmp_120 - tmp_125*tmp_127 - tmp_132*tmp_134 - tmp_139*tmp_141 - tmp_146*tmp_148 - tmp_153*tmp_155 - tmp_160*tmp_162 - tmp_167*tmp_169 - tmp_174*tmp_176 - tmp_181*tmp_183 - tmp_37*tmp_43 - tmp_48*tmp_50 - tmp_55*tmp_57 - tmp_62*tmp_64 - tmp_69*tmp_71 - tmp_76*tmp_78 - tmp_83*tmp_85 - tmp_90*tmp_92 - tmp_97*tmp_99;
      real_t a_0_3 = -tmp_105*tmp_106 - tmp_112*tmp_113 - tmp_119*tmp_120 - tmp_126*tmp_127 - tmp_133*tmp_134 - tmp_140*tmp_141 - tmp_147*tmp_148 - tmp_154*tmp_155 - tmp_161*tmp_162 - tmp_168*tmp_169 - tmp_175*tmp_176 - tmp_182*tmp_183 - tmp_41*tmp_43 - tmp_49*tmp_50 - tmp_56*tmp_57 - tmp_63*tmp_64 - tmp_70*tmp_71 - tmp_77*tmp_78 - tmp_84*tmp_85 - tmp_91*tmp_92 - tmp_98*tmp_99;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Point3D >& coordsElement,
    const std::vector< Point3D >& coordsFacet,
    const Point3D&,
    const Point3D&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
                                        MatrixXr&                            elMat ) const override
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t tmp_0 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_2 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = tmp_0*tmp_1 - tmp_2*tmp_3;
      real_t tmp_5 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = tmp_3*tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_1*tmp_6;
      real_t tmp_13 = tmp_0*tmp_9;
      real_t tmp_14 = tmp_2*tmp_8;
      real_t tmp_15 = 1.0 / (tmp_0*tmp_1*tmp_8 - tmp_10*tmp_12 + tmp_11*tmp_2 - tmp_13*tmp_5 - tmp_14*tmp_3 + tmp_5*tmp_7);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.031405749086161582*tmp_17 + 0.93718850182767688*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_0*tmp_5 + tmp_10*tmp_2;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.031405749086161582*tmp_23 + 0.93718850182767688*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_1*tmp_10 + tmp_3*tmp_5;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_4 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = -tmp_12 + tmp_2*tmp_9;
      real_t tmp_35 = -tmp_14 + tmp_5*tmp_6;
      real_t tmp_36 = tmp_1*tmp_8 - tmp_5*tmp_9;
      real_t tmp_37 = tmp_20*tmp_34 + tmp_26*tmp_35 + tmp_32*tmp_36;
      real_t tmp_38 = -tmp_13 + tmp_7;
      real_t tmp_39 = tmp_0*tmp_8 - tmp_10*tmp_6;
      real_t tmp_40 = tmp_11 - tmp_3*tmp_8;
      real_t tmp_41 = tmp_20*tmp_38 + tmp_26*tmp_39 + tmp_32*tmp_40;
      real_t tmp_42 = 1.0*p_affine_13_1*std::pow((std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)), 1.0/2.0);
      real_t tmp_43 = 0.0068572537431980923*tmp_42;
      real_t tmp_44 = tmp_15*(0.19601935860219369*tmp_17 + 0.60796128279561268*tmp_18 + tmp_19);
      real_t tmp_45 = tmp_15*(0.19601935860219369*tmp_23 + 0.60796128279561268*tmp_24 + tmp_25);
      real_t tmp_46 = tmp_15*(0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_30 + tmp_31);
      real_t tmp_47 = tmp_21*tmp_45 + tmp_27*tmp_46 + tmp_4*tmp_44;
      real_t tmp_48 = tmp_34*tmp_44 + tmp_35*tmp_45 + tmp_36*tmp_46;
      real_t tmp_49 = tmp_38*tmp_44 + tmp_39*tmp_45 + tmp_40*tmp_46;
      real_t tmp_50 = 0.037198804536718075*tmp_42;
      real_t tmp_51 = tmp_15*(0.37605877282253791*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_52 = tmp_15*(0.37605877282253791*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_53 = tmp_15*(0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_54 = tmp_21*tmp_52 + tmp_27*tmp_53 + tmp_4*tmp_51;
      real_t tmp_55 = tmp_34*tmp_51 + tmp_35*tmp_52 + tmp_36*tmp_53;
      real_t tmp_56 = tmp_38*tmp_51 + tmp_39*tmp_52 + tmp_40*tmp_53;
      real_t tmp_57 = 0.020848748529055869*tmp_42;
      real_t tmp_58 = tmp_15*(0.78764240869137092*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_59 = tmp_15*(0.78764240869137092*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_60 = tmp_15*(0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_61 = tmp_21*tmp_59 + tmp_27*tmp_60 + tmp_4*tmp_58;
      real_t tmp_62 = tmp_34*tmp_58 + tmp_35*tmp_59 + tmp_36*tmp_60;
      real_t tmp_63 = tmp_38*tmp_58 + tmp_39*tmp_59 + tmp_40*tmp_60;
      real_t tmp_64 = 0.019202922745021479*tmp_42;
      real_t tmp_65 = tmp_15*(0.58463275527740355*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_66 = tmp_15*(0.58463275527740355*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_67 = tmp_15*(0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_68 = tmp_21*tmp_66 + tmp_27*tmp_67 + tmp_4*tmp_65;
      real_t tmp_69 = tmp_34*tmp_65 + tmp_35*tmp_66 + tmp_36*tmp_67;
      real_t tmp_70 = tmp_38*tmp_65 + tmp_39*tmp_66 + tmp_40*tmp_67;
      real_t tmp_71 = 0.020848748529055869*tmp_42;
      real_t tmp_72 = tmp_15*(0.041227165399737475*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_73 = tmp_15*(0.041227165399737475*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_74 = tmp_15*(0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_75 = tmp_21*tmp_73 + tmp_27*tmp_74 + tmp_4*tmp_72;
      real_t tmp_76 = tmp_34*tmp_72 + tmp_35*tmp_73 + tmp_36*tmp_74;
      real_t tmp_77 = tmp_38*tmp_72 + tmp_39*tmp_73 + tmp_40*tmp_74;
      real_t tmp_78 = 0.019202922745021479*tmp_42;
      real_t tmp_79 = tmp_15*(0.039308471900058539*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_80 = tmp_15*(0.039308471900058539*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_81 = tmp_15*(0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_82 = tmp_21*tmp_80 + tmp_27*tmp_81 + tmp_4*tmp_79;
      real_t tmp_83 = tmp_34*tmp_79 + tmp_35*tmp_80 + tmp_36*tmp_81;
      real_t tmp_84 = tmp_38*tmp_79 + tmp_39*tmp_80 + tmp_40*tmp_81;
      real_t tmp_85 = 0.020848748529055869*tmp_42;
      real_t tmp_86 = tmp_15*(0.78764240869137092*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_87 = tmp_15*(0.78764240869137092*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_88 = tmp_15*(0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_89 = tmp_21*tmp_87 + tmp_27*tmp_88 + tmp_4*tmp_86;
      real_t tmp_90 = tmp_34*tmp_86 + tmp_35*tmp_87 + tmp_36*tmp_88;
      real_t tmp_91 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88;
      real_t tmp_92 = 0.019202922745021479*tmp_42;
      real_t tmp_93 = tmp_15*(0.58463275527740355*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_94 = tmp_15*(0.58463275527740355*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_95 = tmp_15*(0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_96 = tmp_21*tmp_94 + tmp_27*tmp_95 + tmp_4*tmp_93;
      real_t tmp_97 = tmp_34*tmp_93 + tmp_35*tmp_94 + tmp_36*tmp_95;
      real_t tmp_98 = tmp_38*tmp_93 + tmp_39*tmp_94 + tmp_40*tmp_95;
      real_t tmp_99 = 0.020848748529055869*tmp_42;
      real_t tmp_100 = tmp_15*(0.1711304259088916*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_101 = tmp_15*(0.1711304259088916*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_102 = tmp_15*(0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_103 = tmp_100*tmp_4 + tmp_101*tmp_21 + tmp_102*tmp_27;
      real_t tmp_104 = tmp_100*tmp_34 + tmp_101*tmp_35 + tmp_102*tmp_36;
      real_t tmp_105 = tmp_100*tmp_38 + tmp_101*tmp_39 + tmp_102*tmp_40;
      real_t tmp_106 = 0.019202922745021479*tmp_42;
      real_t tmp_107 = tmp_15*(0.37605877282253791*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_108 = tmp_15*(0.37605877282253791*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_109 = tmp_15*(0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_110 = tmp_107*tmp_4 + tmp_108*tmp_21 + tmp_109*tmp_27;
      real_t tmp_111 = tmp_107*tmp_34 + tmp_108*tmp_35 + tmp_109*tmp_36;
      real_t tmp_112 = tmp_107*tmp_38 + tmp_108*tmp_39 + tmp_109*tmp_40;
      real_t tmp_113 = 0.020848748529055869*tmp_42;
      real_t tmp_114 = tmp_15*(0.041227165399737475*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_115 = tmp_15*(0.041227165399737475*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_116 = tmp_15*(0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_117 = tmp_114*tmp_4 + tmp_115*tmp_21 + tmp_116*tmp_27;
      real_t tmp_118 = tmp_114*tmp_34 + tmp_115*tmp_35 + tmp_116*tmp_36;
      real_t tmp_119 = tmp_114*tmp_38 + tmp_115*tmp_39 + tmp_116*tmp_40;
      real_t tmp_120 = 0.019202922745021479*tmp_42;
      real_t tmp_121 = tmp_15*(0.40446199974765351*tmp_17 + 0.19107600050469298*tmp_18 + tmp_19);
      real_t tmp_122 = tmp_15*(0.40446199974765351*tmp_23 + 0.19107600050469298*tmp_24 + tmp_25);
      real_t tmp_123 = tmp_15*(0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_30 + tmp_31);
      real_t tmp_124 = tmp_121*tmp_4 + tmp_122*tmp_21 + tmp_123*tmp_27;
      real_t tmp_125 = tmp_121*tmp_34 + tmp_122*tmp_35 + tmp_123*tmp_36;
      real_t tmp_126 = tmp_121*tmp_38 + tmp_122*tmp_39 + tmp_123*tmp_40;
      real_t tmp_127 = 0.042507265838595799*tmp_42;
      real_t tmp_128 = tmp_15*(0.039308471900058539*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_129 = tmp_15*(0.039308471900058539*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_130 = tmp_15*(0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_131 = tmp_128*tmp_4 + tmp_129*tmp_21 + tmp_130*tmp_27;
      real_t tmp_132 = tmp_128*tmp_34 + tmp_129*tmp_35 + tmp_130*tmp_36;
      real_t tmp_133 = tmp_128*tmp_38 + tmp_129*tmp_39 + tmp_130*tmp_40;
      real_t tmp_134 = 0.020848748529055869*tmp_42;
      real_t tmp_135 = tmp_15*(0.93718850182767688*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_136 = tmp_15*(0.93718850182767688*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_137 = tmp_15*(0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_138 = tmp_135*tmp_4 + tmp_136*tmp_21 + tmp_137*tmp_27;
      real_t tmp_139 = tmp_135*tmp_34 + tmp_136*tmp_35 + tmp_137*tmp_36;
      real_t tmp_140 = tmp_135*tmp_38 + tmp_136*tmp_39 + tmp_137*tmp_40;
      real_t tmp_141 = 0.0068572537431980923*tmp_42;
      real_t tmp_142 = tmp_15*(0.60796128279561268*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_143 = tmp_15*(0.60796128279561268*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_144 = tmp_15*(0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_145 = tmp_142*tmp_4 + tmp_143*tmp_21 + tmp_144*tmp_27;
      real_t tmp_146 = tmp_142*tmp_34 + tmp_143*tmp_35 + tmp_144*tmp_36;
      real_t tmp_147 = tmp_142*tmp_38 + tmp_143*tmp_39 + tmp_144*tmp_40;
      real_t tmp_148 = 0.037198804536718075*tmp_42;
      real_t tmp_149 = tmp_15*(0.19107600050469298*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_150 = tmp_15*(0.19107600050469298*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_151 = tmp_15*(0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_152 = tmp_149*tmp_4 + tmp_150*tmp_21 + tmp_151*tmp_27;
      real_t tmp_153 = tmp_149*tmp_34 + tmp_150*tmp_35 + tmp_151*tmp_36;
      real_t tmp_154 = tmp_149*tmp_38 + tmp_150*tmp_39 + tmp_151*tmp_40;
      real_t tmp_155 = 0.042507265838595799*tmp_42;
      real_t tmp_156 = tmp_15*(0.031405749086161582*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_157 = tmp_15*(0.031405749086161582*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_158 = tmp_15*(0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_159 = tmp_156*tmp_4 + tmp_157*tmp_21 + tmp_158*tmp_27;
      real_t tmp_160 = tmp_156*tmp_34 + tmp_157*tmp_35 + tmp_158*tmp_36;
      real_t tmp_161 = tmp_156*tmp_38 + tmp_157*tmp_39 + tmp_158*tmp_40;
      real_t tmp_162 = 0.0068572537431980923*tmp_42;
      real_t tmp_163 = tmp_15*(0.19601935860219369*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_164 = tmp_15*(0.19601935860219369*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_165 = tmp_15*(0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_166 = tmp_163*tmp_4 + tmp_164*tmp_21 + tmp_165*tmp_27;
      real_t tmp_167 = tmp_163*tmp_34 + tmp_164*tmp_35 + tmp_165*tmp_36;
      real_t tmp_168 = tmp_163*tmp_38 + tmp_164*tmp_39 + tmp_165*tmp_40;
      real_t tmp_169 = 0.037198804536718075*tmp_42;
      real_t tmp_170 = tmp_15*(0.40446199974765351*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_171 = tmp_15*(0.40446199974765351*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_172 = tmp_15*(0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_173 = tmp_170*tmp_4 + tmp_171*tmp_21 + tmp_172*tmp_27;
      real_t tmp_174 = tmp_170*tmp_34 + tmp_171*tmp_35 + tmp_172*tmp_36;
      real_t tmp_175 = tmp_170*tmp_38 + tmp_171*tmp_39 + tmp_172*tmp_40;
      real_t tmp_176 = 0.042507265838595799*tmp_42;
      real_t tmp_177 = tmp_15*(0.1711304259088916*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_178 = tmp_15*(0.1711304259088916*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_179 = tmp_15*(0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_180 = tmp_177*tmp_4 + tmp_178*tmp_21 + tmp_179*tmp_27;
      real_t tmp_181 = tmp_177*tmp_34 + tmp_178*tmp_35 + tmp_179*tmp_36;
      real_t tmp_182 = tmp_177*tmp_38 + tmp_178*tmp_39 + tmp_179*tmp_40;
      real_t tmp_183 = 0.019202922745021479*tmp_42;
      real_t a_0_0 = tmp_106*(-tmp_103 - tmp_104 - tmp_105 + 1) + tmp_113*(-tmp_110 - tmp_111 - tmp_112 + 1) + tmp_120*(-tmp_117 - tmp_118 - tmp_119 + 1) + tmp_127*(-tmp_124 - tmp_125 - tmp_126 + 1) + tmp_134*(-tmp_131 - tmp_132 - tmp_133 + 1) + tmp_141*(-tmp_138 - tmp_139 - tmp_140 + 1) + tmp_148*(-tmp_145 - tmp_146 - tmp_147 + 1) + tmp_155*(-tmp_152 - tmp_153 - tmp_154 + 1) + tmp_162*(-tmp_159 - tmp_160 - tmp_161 + 1) + tmp_169*(-tmp_166 - tmp_167 - tmp_168 + 1) + tmp_176*(-tmp_173 - tmp_174 - tmp_175 + 1) + tmp_183*(-tmp_180 - tmp_181 - tmp_182 + 1) + tmp_43*(-tmp_33 - tmp_37 - tmp_41 + 1) + tmp_50*(-tmp_47 - tmp_48 - tmp_49 + 1) + tmp_57*(-tmp_54 - tmp_55 - tmp_56 + 1) + tmp_64*(-tmp_61 - tmp_62 - tmp_63 + 1) + tmp_71*(-tmp_68 - tmp_69 - tmp_70 + 1) + tmp_78*(-tmp_75 - tmp_76 - tmp_77 + 1) + tmp_85*(-tmp_82 - tmp_83 - tmp_84 + 1) + tmp_92*(-tmp_89 - tmp_90 - tmp_91 + 1) + tmp_99*(-tmp_96 - tmp_97 - tmp_98 + 1);
      real_t a_0_1 = tmp_103*tmp_106 + tmp_110*tmp_113 + tmp_117*tmp_120 + tmp_124*tmp_127 + tmp_131*tmp_134 + tmp_138*tmp_141 + tmp_145*tmp_148 + tmp_152*tmp_155 + tmp_159*tmp_162 + tmp_166*tmp_169 + tmp_173*tmp_176 + tmp_180*tmp_183 + tmp_33*tmp_43 + tmp_47*tmp_50 + tmp_54*tmp_57 + tmp_61*tmp_64 + tmp_68*tmp_71 + tmp_75*tmp_78 + tmp_82*tmp_85 + tmp_89*tmp_92 + tmp_96*tmp_99;
      real_t a_0_2 = tmp_104*tmp_106 + tmp_111*tmp_113 + tmp_118*tmp_120 + tmp_125*tmp_127 + tmp_132*tmp_134 + tmp_139*tmp_141 + tmp_146*tmp_148 + tmp_153*tmp_155 + tmp_160*tmp_162 + tmp_167*tmp_169 + tmp_174*tmp_176 + tmp_181*tmp_183 + tmp_37*tmp_43 + tmp_48*tmp_50 + tmp_55*tmp_57 + tmp_62*tmp_64 + tmp_69*tmp_71 + tmp_76*tmp_78 + tmp_83*tmp_85 + tmp_90*tmp_92 + tmp_97*tmp_99;
      real_t a_0_3 = tmp_105*tmp_106 + tmp_112*tmp_113 + tmp_119*tmp_120 + tmp_126*tmp_127 + tmp_133*tmp_134 + tmp_140*tmp_141 + tmp_147*tmp_148 + tmp_154*tmp_155 + tmp_161*tmp_162 + tmp_168*tmp_169 + tmp_175*tmp_176 + tmp_182*tmp_183 + tmp_41*tmp_43 + tmp_49*tmp_50 + tmp_56*tmp_57 + tmp_63*tmp_64 + tmp_70*tmp_71 + tmp_77*tmp_78 + tmp_84*tmp_85 + tmp_91*tmp_92 + tmp_98*tmp_99;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }

public:

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g1;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g1;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g0;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g0;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g2;

};




class EGDivFormNitscheBC_P0P1_2 : public hyteg::dg::DGForm
{

 public:
    EGDivFormNitscheBC_P0P1_2()
: callback_Scalar_Variable_Coefficient_3D_g1 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_2D_g1 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_2D_g0 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_3D_g0 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_3D_g2 ([](const Point3D & p) -> real_t { return 0.; })
    {}

void Scalar_Variable_Coefficient_2D_g0( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g0( Point3D( {in_0, in_1, 0} ) );
}
void Scalar_Variable_Coefficient_2D_g1( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g1( Point3D( {in_0, in_1, 0} ) );
}

void Scalar_Variable_Coefficient_3D_g0( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g0( Point3D( {in_0, in_1, in_2} ) );
}
void Scalar_Variable_Coefficient_3D_g1( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g1( Point3D( {in_0, in_1, in_2} ) );
}
void Scalar_Variable_Coefficient_3D_g2( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g2( Point3D( {in_0, in_1, in_2} ) );
}

 protected:
  void integrateVolume2D( const std::vector< Point3D >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >&      coordsElement,
                                       const std::vector< Point3D >&      coordsFacet,
                                       const Point3D&                     oppositeVertex,
                                       const Point3D&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >&      coordsElementInner,
                                          const std::vector< Point3D >&      coordsElementOuter,
                                          const std::vector< Point3D >&      coordsFacet,
                                          const Point3D&                     oppositeVertexInnerElement,
                                          const Point3D&                     oppositeVertexOuterElement,
                                          const Point3D&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                   const std::vector< Point3D >&      coordsFacet,
                                                   const Point3D&                     oppositeVertex,
                                                   const Point3D&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
   }

    void integrateRHSDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

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

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_g0_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id9 = 0;
      Scalar_Variable_Coefficient_2D_g0( 0.95308992296933193*p_affine_6_0 + 0.046910077030668018*p_affine_7_0, 0.95308992296933193*p_affine_6_1 + 0.046910077030668018*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id0 );
      Scalar_Variable_Coefficient_2D_g1( 0.95308992296933193*p_affine_6_0 + 0.046910077030668018*p_affine_7_0, 0.95308992296933193*p_affine_6_1 + 0.046910077030668018*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id1 );
      Scalar_Variable_Coefficient_2D_g0( 0.7692346550528415*p_affine_6_0 + 0.23076534494715845*p_affine_7_0, 0.7692346550528415*p_affine_6_1 + 0.23076534494715845*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id2 );
      Scalar_Variable_Coefficient_2D_g1( 0.7692346550528415*p_affine_6_0 + 0.23076534494715845*p_affine_7_0, 0.7692346550528415*p_affine_6_1 + 0.23076534494715845*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id3 );
      Scalar_Variable_Coefficient_2D_g0( 0.5*p_affine_6_0 + 0.5*p_affine_7_0, 0.5*p_affine_6_1 + 0.5*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id4 );
      Scalar_Variable_Coefficient_2D_g1( 0.5*p_affine_6_0 + 0.5*p_affine_7_0, 0.5*p_affine_6_1 + 0.5*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id5 );
      Scalar_Variable_Coefficient_2D_g0( 0.2307653449471585*p_affine_6_0 + 0.7692346550528415*p_affine_7_0, 0.2307653449471585*p_affine_6_1 + 0.7692346550528415*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id6 );
      Scalar_Variable_Coefficient_2D_g1( 0.2307653449471585*p_affine_6_0 + 0.7692346550528415*p_affine_7_0, 0.2307653449471585*p_affine_6_1 + 0.7692346550528415*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id7 );
      Scalar_Variable_Coefficient_2D_g0( 0.046910077030668074*p_affine_6_0 + 0.95308992296933193*p_affine_7_0, 0.046910077030668074*p_affine_6_1 + 0.95308992296933193*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id8 );
      Scalar_Variable_Coefficient_2D_g1( 0.046910077030668074*p_affine_6_0 + 0.95308992296933193*p_affine_7_0, 0.046910077030668074*p_affine_6_1 + 0.95308992296933193*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id9 );
      real_t tmp_0 = std::abs(std::pow(((-p_affine_6_0 + p_affine_7_0)*(-p_affine_6_0 + p_affine_7_0)) + ((-p_affine_6_1 + p_affine_7_1)*(-p_affine_6_1 + p_affine_7_1)), 1.0/2.0));
      real_t a_0_0 = 0.11846344252809471*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id0*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id1*p_affine_10_1) + 0.2393143352496831*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id2*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id3*p_affine_10_1) + 0.2844444444444445*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id4*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id5*p_affine_10_1) + 0.2393143352496831*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id6*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id7*p_affine_10_1) + 0.11846344252809471*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id8*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id9*p_affine_10_1);
      elMat( 0, 0) = a_0_0;
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

      real_t Scalar_Variable_Coefficient_3D_g0_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id13 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id14 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id15 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id16 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id17 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id18 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id19 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id20 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id21 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id22 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id23 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id24 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id25 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id26 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id27 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id28 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id29 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id30 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id31 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id32 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id33 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id34 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id35 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id36 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id37 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id38 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id39 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id40 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id41 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id42 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id43 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id44 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id45 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id46 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id47 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id48 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id49 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id50 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id51 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id52 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id53 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id54 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id55 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id56 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id57 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id58 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id59 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id60 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id61 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id62 = 0;
      Scalar_Variable_Coefficient_3D_g0( 0.93718850182767688*p_affine_10_0 + 0.031405749086161561*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.93718850182767688*p_affine_10_1 + 0.031405749086161561*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.93718850182767688*p_affine_10_2 + 0.031405749086161561*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id0 );
      Scalar_Variable_Coefficient_3D_g1( 0.93718850182767688*p_affine_10_0 + 0.031405749086161561*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.93718850182767688*p_affine_10_1 + 0.031405749086161561*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.93718850182767688*p_affine_10_2 + 0.031405749086161561*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id1 );
      Scalar_Variable_Coefficient_3D_g2( 0.93718850182767688*p_affine_10_0 + 0.031405749086161561*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.93718850182767688*p_affine_10_1 + 0.031405749086161561*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.93718850182767688*p_affine_10_2 + 0.031405749086161561*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id2 );
      Scalar_Variable_Coefficient_3D_g0( 0.60796128279561268*p_affine_10_0 + 0.19601935860219366*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.60796128279561268*p_affine_10_1 + 0.19601935860219366*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.60796128279561268*p_affine_10_2 + 0.19601935860219366*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id3 );
      Scalar_Variable_Coefficient_3D_g1( 0.60796128279561268*p_affine_10_0 + 0.19601935860219366*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.60796128279561268*p_affine_10_1 + 0.19601935860219366*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.60796128279561268*p_affine_10_2 + 0.19601935860219366*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id4 );
      Scalar_Variable_Coefficient_3D_g2( 0.60796128279561268*p_affine_10_0 + 0.19601935860219366*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.60796128279561268*p_affine_10_1 + 0.19601935860219366*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.60796128279561268*p_affine_10_2 + 0.19601935860219366*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id5 );
      Scalar_Variable_Coefficient_3D_g0( 0.039308471900058539*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id6 );
      Scalar_Variable_Coefficient_3D_g1( 0.039308471900058539*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id7 );
      Scalar_Variable_Coefficient_3D_g2( 0.039308471900058539*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id8 );
      Scalar_Variable_Coefficient_3D_g0( 0.1711304259088916*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id9 );
      Scalar_Variable_Coefficient_3D_g1( 0.1711304259088916*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id10 );
      Scalar_Variable_Coefficient_3D_g2( 0.1711304259088916*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id11 );
      Scalar_Variable_Coefficient_3D_g0( 0.37605877282253791*p_affine_10_0 + 0.039308471900058539*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.039308471900058539*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.039308471900058539*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id12 );
      Scalar_Variable_Coefficient_3D_g1( 0.37605877282253791*p_affine_10_0 + 0.039308471900058539*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.039308471900058539*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.039308471900058539*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id13 );
      Scalar_Variable_Coefficient_3D_g2( 0.37605877282253791*p_affine_10_0 + 0.039308471900058539*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.039308471900058539*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.039308471900058539*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id14 );
      Scalar_Variable_Coefficient_3D_g0( 0.78764240869137092*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id15 );
      Scalar_Variable_Coefficient_3D_g1( 0.78764240869137092*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id16 );
      Scalar_Variable_Coefficient_3D_g2( 0.78764240869137092*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id17 );
      Scalar_Variable_Coefficient_3D_g0( 0.58463275527740355*p_affine_10_0 + 0.37605877282253797*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.37605877282253797*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.37605877282253797*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id18 );
      Scalar_Variable_Coefficient_3D_g1( 0.58463275527740355*p_affine_10_0 + 0.37605877282253797*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.37605877282253797*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.37605877282253797*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id19 );
      Scalar_Variable_Coefficient_3D_g2( 0.58463275527740355*p_affine_10_0 + 0.37605877282253797*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.37605877282253797*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.37605877282253797*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id20 );
      Scalar_Variable_Coefficient_3D_g0( 0.041227165399737475*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id21 );
      Scalar_Variable_Coefficient_3D_g1( 0.041227165399737475*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id22 );
      Scalar_Variable_Coefficient_3D_g2( 0.041227165399737475*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id23 );
      Scalar_Variable_Coefficient_3D_g0( 0.039308471900058539*p_affine_10_0 + 0.37605877282253791*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.37605877282253791*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.37605877282253791*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id24 );
      Scalar_Variable_Coefficient_3D_g1( 0.039308471900058539*p_affine_10_0 + 0.37605877282253791*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.37605877282253791*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.37605877282253791*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id25 );
      Scalar_Variable_Coefficient_3D_g2( 0.039308471900058539*p_affine_10_0 + 0.37605877282253791*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.37605877282253791*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.37605877282253791*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id26 );
      Scalar_Variable_Coefficient_3D_g0( 0.78764240869137092*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id27 );
      Scalar_Variable_Coefficient_3D_g1( 0.78764240869137092*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id28 );
      Scalar_Variable_Coefficient_3D_g2( 0.78764240869137092*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id29 );
      Scalar_Variable_Coefficient_3D_g0( 0.58463275527740355*p_affine_10_0 + 0.039308471900058484*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.039308471900058484*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.039308471900058484*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id30 );
      Scalar_Variable_Coefficient_3D_g1( 0.58463275527740355*p_affine_10_0 + 0.039308471900058484*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.039308471900058484*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.039308471900058484*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id31 );
      Scalar_Variable_Coefficient_3D_g2( 0.58463275527740355*p_affine_10_0 + 0.039308471900058484*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.039308471900058484*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.039308471900058484*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id32 );
      Scalar_Variable_Coefficient_3D_g0( 0.1711304259088916*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id33 );
      Scalar_Variable_Coefficient_3D_g1( 0.1711304259088916*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id34 );
      Scalar_Variable_Coefficient_3D_g2( 0.1711304259088916*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id35 );
      Scalar_Variable_Coefficient_3D_g0( 0.19107600050469298*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.19107600050469298*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.19107600050469298*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id36 );
      Scalar_Variable_Coefficient_3D_g1( 0.19107600050469298*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.19107600050469298*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.19107600050469298*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id37 );
      Scalar_Variable_Coefficient_3D_g2( 0.19107600050469298*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.19107600050469298*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.19107600050469298*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id38 );
      Scalar_Variable_Coefficient_3D_g0( 0.37605877282253791*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id39 );
      Scalar_Variable_Coefficient_3D_g1( 0.37605877282253791*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id40 );
      Scalar_Variable_Coefficient_3D_g2( 0.37605877282253791*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id41 );
      Scalar_Variable_Coefficient_3D_g0( 0.031405749086161582*p_affine_10_0 + 0.03140574908616154*p_affine_8_0 + 0.93718850182767688*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.03140574908616154*p_affine_8_1 + 0.93718850182767688*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.03140574908616154*p_affine_8_2 + 0.93718850182767688*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id42 );
      Scalar_Variable_Coefficient_3D_g1( 0.031405749086161582*p_affine_10_0 + 0.03140574908616154*p_affine_8_0 + 0.93718850182767688*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.03140574908616154*p_affine_8_1 + 0.93718850182767688*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.03140574908616154*p_affine_8_2 + 0.93718850182767688*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id43 );
      Scalar_Variable_Coefficient_3D_g2( 0.031405749086161582*p_affine_10_0 + 0.03140574908616154*p_affine_8_0 + 0.93718850182767688*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.03140574908616154*p_affine_8_1 + 0.93718850182767688*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.03140574908616154*p_affine_8_2 + 0.93718850182767688*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id44 );
      Scalar_Variable_Coefficient_3D_g0( 0.19601935860219369*p_affine_10_0 + 0.19601935860219363*p_affine_8_0 + 0.60796128279561268*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.19601935860219363*p_affine_8_1 + 0.60796128279561268*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.19601935860219363*p_affine_8_2 + 0.60796128279561268*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id45 );
      Scalar_Variable_Coefficient_3D_g1( 0.19601935860219369*p_affine_10_0 + 0.19601935860219363*p_affine_8_0 + 0.60796128279561268*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.19601935860219363*p_affine_8_1 + 0.60796128279561268*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.19601935860219363*p_affine_8_2 + 0.60796128279561268*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id46 );
      Scalar_Variable_Coefficient_3D_g2( 0.19601935860219369*p_affine_10_0 + 0.19601935860219363*p_affine_8_0 + 0.60796128279561268*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.19601935860219363*p_affine_8_1 + 0.60796128279561268*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.19601935860219363*p_affine_8_2 + 0.60796128279561268*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id47 );
      Scalar_Variable_Coefficient_3D_g0( 0.40446199974765351*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.19107600050469298*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.19107600050469298*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.19107600050469298*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id48 );
      Scalar_Variable_Coefficient_3D_g1( 0.40446199974765351*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.19107600050469298*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.19107600050469298*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.19107600050469298*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id49 );
      Scalar_Variable_Coefficient_3D_g2( 0.40446199974765351*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.19107600050469298*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.19107600050469298*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.19107600050469298*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id50 );
      Scalar_Variable_Coefficient_3D_g0( 0.031405749086161582*p_affine_10_0 + 0.93718850182767688*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.93718850182767688*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.93718850182767688*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id51 );
      Scalar_Variable_Coefficient_3D_g1( 0.031405749086161582*p_affine_10_0 + 0.93718850182767688*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.93718850182767688*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.93718850182767688*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id52 );
      Scalar_Variable_Coefficient_3D_g2( 0.031405749086161582*p_affine_10_0 + 0.93718850182767688*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.93718850182767688*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.93718850182767688*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id53 );
      Scalar_Variable_Coefficient_3D_g0( 0.19601935860219369*p_affine_10_0 + 0.60796128279561268*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.60796128279561268*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.60796128279561268*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id54 );
      Scalar_Variable_Coefficient_3D_g1( 0.19601935860219369*p_affine_10_0 + 0.60796128279561268*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.60796128279561268*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.60796128279561268*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id55 );
      Scalar_Variable_Coefficient_3D_g2( 0.19601935860219369*p_affine_10_0 + 0.60796128279561268*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.60796128279561268*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.60796128279561268*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id56 );
      Scalar_Variable_Coefficient_3D_g0( 0.40446199974765351*p_affine_10_0 + 0.19107600050469298*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.19107600050469298*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.19107600050469298*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id57 );
      Scalar_Variable_Coefficient_3D_g1( 0.40446199974765351*p_affine_10_0 + 0.19107600050469298*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.19107600050469298*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.19107600050469298*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id58 );
      Scalar_Variable_Coefficient_3D_g2( 0.40446199974765351*p_affine_10_0 + 0.19107600050469298*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.19107600050469298*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.19107600050469298*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id59 );
      Scalar_Variable_Coefficient_3D_g0( 0.041227165399737475*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id60 );
      Scalar_Variable_Coefficient_3D_g1( 0.041227165399737475*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id61 );
      Scalar_Variable_Coefficient_3D_g2( 0.041227165399737475*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id62 );
      real_t tmp_0 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_1 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_2 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_3 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_4 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_5 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_6 = 1.0*std::pow((std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)*std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)) + (std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)*std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)) + (std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)*std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)), 1.0/2.0);
      real_t a_0_0 = 0.0068572537431980923*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id0*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id1*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id2*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id12*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id13*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id14*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id15*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id16*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id17*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id18*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id19*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id20*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id21*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id22*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id23*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id24*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id25*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id26*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id27*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id28*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id29*p_affine_13_2) + 0.037198804536718075*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id3*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id4*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id5*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id30*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id31*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id32*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id33*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id34*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id35*p_affine_13_2) + 0.042507265838595799*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id36*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id37*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id38*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id39*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id40*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id41*p_affine_13_2) + 0.0068572537431980923*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id42*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id43*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id44*p_affine_13_2) + 0.037198804536718075*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id45*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id46*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id47*p_affine_13_2) + 0.042507265838595799*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id48*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id49*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id50*p_affine_13_2) + 0.0068572537431980923*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id51*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id52*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id53*p_affine_13_2) + 0.037198804536718075*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id54*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id55*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id56*p_affine_13_2) + 0.042507265838595799*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id57*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id58*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id59*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id6*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id7*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id8*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id60*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id61*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id62*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id9*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id10*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id11*p_affine_13_2);
      elMat( 0, 0) = a_0_0;
   }
   void integrateVolume3D( const std::vector< Point3D >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coords[0]( 0 );
   const auto p_affine_0_1 = coords[0]( 1 );
   const auto p_affine_0_2 = coords[0]( 2 );

   const auto p_affine_1_0 = coords[1]( 0 );
   const auto p_affine_1_1 = coords[1]( 1 );
   const auto p_affine_1_2 = coords[1]( 2 );

   const auto p_affine_2_0 = coords[2]( 0 );
   const auto p_affine_2_1 = coords[2]( 1 );
   const auto p_affine_2_2 = coords[2]( 2 );

   const auto p_affine_3_0 = coords[3]( 0 );
   const auto p_affine_3_1 = coords[3]( 1 );
   const auto p_affine_3_2 = coords[3]( 2 );

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = tmp_0*tmp_1;
      real_t tmp_3 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_4 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_7 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_8 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_0*tmp_8;
      real_t tmp_13 = tmp_1*tmp_11;
      real_t tmp_14 = 1.0 / (tmp_10*tmp_11*tmp_4 - tmp_10*tmp_12 - tmp_13*tmp_7 + tmp_2*tmp_6 - tmp_5*tmp_6 + tmp_7*tmp_9);
      real_t tmp_15 = tmp_14*(tmp_2 - tmp_5);
      real_t tmp_16 = tmp_14*(tmp_11*tmp_4 - tmp_12);
      real_t tmp_17 = tmp_14*(-tmp_13 + tmp_9);
      real_t tmp_18 = p_affine_0_0*p_affine_1_1;
      real_t tmp_19 = p_affine_0_0*p_affine_1_2;
      real_t tmp_20 = p_affine_2_1*p_affine_3_2;
      real_t tmp_21 = p_affine_0_1*p_affine_1_0;
      real_t tmp_22 = p_affine_0_1*p_affine_1_2;
      real_t tmp_23 = p_affine_2_2*p_affine_3_0;
      real_t tmp_24 = p_affine_0_2*p_affine_1_0;
      real_t tmp_25 = p_affine_0_2*p_affine_1_1;
      real_t tmp_26 = p_affine_2_0*p_affine_3_1;
      real_t tmp_27 = p_affine_2_2*p_affine_3_1;
      real_t tmp_28 = p_affine_2_0*p_affine_3_2;
      real_t tmp_29 = p_affine_2_1*p_affine_3_0;
      real_t tmp_30 = std::abs(p_affine_0_0*tmp_20 - p_affine_0_0*tmp_27 + p_affine_0_1*tmp_23 - p_affine_0_1*tmp_28 + p_affine_0_2*tmp_26 - p_affine_0_2*tmp_29 - p_affine_1_0*tmp_20 + p_affine_1_0*tmp_27 - p_affine_1_1*tmp_23 + p_affine_1_1*tmp_28 - p_affine_1_2*tmp_26 + p_affine_1_2*tmp_29 + p_affine_2_0*tmp_22 - p_affine_2_0*tmp_25 - p_affine_2_1*tmp_19 + p_affine_2_1*tmp_24 + p_affine_2_2*tmp_18 - p_affine_2_2*tmp_21 - p_affine_3_0*tmp_22 + p_affine_3_0*tmp_25 + p_affine_3_1*tmp_19 - p_affine_3_1*tmp_24 - p_affine_3_2*tmp_18 + p_affine_3_2*tmp_21);
      real_t tmp_31 = tmp_30*(tmp_15 + tmp_16 + tmp_17);
      real_t tmp_32 = tmp_17*tmp_30;
      real_t tmp_33 = tmp_16*tmp_30;
      real_t tmp_34 = tmp_15*tmp_30;
      real_t a_0_0 = 0.1666666666666668*tmp_31;
      real_t a_0_1 = -0.1666666666666668*tmp_32;
      real_t a_0_2 = -0.1666666666666668*tmp_33;
      real_t a_0_3 = -0.1666666666666668*tmp_34;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }



   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                                                     const std::vector< Point3D >& coordsFacet,
                                                     const Point3D&,
                                                     const Point3D&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                               MatrixXr&                            elMat ) const override
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

         real_t tmp_0 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_2 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = tmp_0*tmp_1 - tmp_2*tmp_3;
      real_t tmp_5 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = tmp_3*tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_1*tmp_6;
      real_t tmp_13 = tmp_0*tmp_9;
      real_t tmp_14 = tmp_2*tmp_8;
      real_t tmp_15 = 1.0 / (tmp_0*tmp_1*tmp_8 - tmp_10*tmp_12 + tmp_11*tmp_2 - tmp_13*tmp_5 - tmp_14*tmp_3 + tmp_5*tmp_7);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.031405749086161582*tmp_17 + 0.93718850182767688*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_0*tmp_5 + tmp_10*tmp_2;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.031405749086161582*tmp_23 + 0.93718850182767688*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_1*tmp_10 + tmp_3*tmp_5;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_4 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = -tmp_12 + tmp_2*tmp_9;
      real_t tmp_35 = -tmp_14 + tmp_5*tmp_6;
      real_t tmp_36 = tmp_1*tmp_8 - tmp_5*tmp_9;
      real_t tmp_37 = tmp_20*tmp_34 + tmp_26*tmp_35 + tmp_32*tmp_36;
      real_t tmp_38 = -tmp_13 + tmp_7;
      real_t tmp_39 = tmp_0*tmp_8 - tmp_10*tmp_6;
      real_t tmp_40 = tmp_11 - tmp_3*tmp_8;
      real_t tmp_41 = tmp_20*tmp_38 + tmp_26*tmp_39 + tmp_32*tmp_40;
      real_t tmp_42 = 0.5*p_affine_13_2*std::pow((std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)), 1.0/2.0);
      real_t tmp_43 = 0.0068572537431980923*tmp_42;
      real_t tmp_44 = tmp_15*(0.19601935860219369*tmp_17 + 0.60796128279561268*tmp_18 + tmp_19);
      real_t tmp_45 = tmp_15*(0.19601935860219369*tmp_23 + 0.60796128279561268*tmp_24 + tmp_25);
      real_t tmp_46 = tmp_15*(0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_30 + tmp_31);
      real_t tmp_47 = tmp_21*tmp_45 + tmp_27*tmp_46 + tmp_4*tmp_44;
      real_t tmp_48 = tmp_34*tmp_44 + tmp_35*tmp_45 + tmp_36*tmp_46;
      real_t tmp_49 = tmp_38*tmp_44 + tmp_39*tmp_45 + tmp_40*tmp_46;
      real_t tmp_50 = 0.037198804536718075*tmp_42;
      real_t tmp_51 = tmp_15*(0.37605877282253791*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_52 = tmp_15*(0.37605877282253791*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_53 = tmp_15*(0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_54 = tmp_21*tmp_52 + tmp_27*tmp_53 + tmp_4*tmp_51;
      real_t tmp_55 = tmp_34*tmp_51 + tmp_35*tmp_52 + tmp_36*tmp_53;
      real_t tmp_56 = tmp_38*tmp_51 + tmp_39*tmp_52 + tmp_40*tmp_53;
      real_t tmp_57 = 0.020848748529055869*tmp_42;
      real_t tmp_58 = tmp_15*(0.78764240869137092*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_59 = tmp_15*(0.78764240869137092*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_60 = tmp_15*(0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_61 = tmp_21*tmp_59 + tmp_27*tmp_60 + tmp_4*tmp_58;
      real_t tmp_62 = tmp_34*tmp_58 + tmp_35*tmp_59 + tmp_36*tmp_60;
      real_t tmp_63 = tmp_38*tmp_58 + tmp_39*tmp_59 + tmp_40*tmp_60;
      real_t tmp_64 = 0.019202922745021479*tmp_42;
      real_t tmp_65 = tmp_15*(0.58463275527740355*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_66 = tmp_15*(0.58463275527740355*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_67 = tmp_15*(0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_68 = tmp_21*tmp_66 + tmp_27*tmp_67 + tmp_4*tmp_65;
      real_t tmp_69 = tmp_34*tmp_65 + tmp_35*tmp_66 + tmp_36*tmp_67;
      real_t tmp_70 = tmp_38*tmp_65 + tmp_39*tmp_66 + tmp_40*tmp_67;
      real_t tmp_71 = 0.020848748529055869*tmp_42;
      real_t tmp_72 = tmp_15*(0.041227165399737475*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_73 = tmp_15*(0.041227165399737475*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_74 = tmp_15*(0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_75 = tmp_21*tmp_73 + tmp_27*tmp_74 + tmp_4*tmp_72;
      real_t tmp_76 = tmp_34*tmp_72 + tmp_35*tmp_73 + tmp_36*tmp_74;
      real_t tmp_77 = tmp_38*tmp_72 + tmp_39*tmp_73 + tmp_40*tmp_74;
      real_t tmp_78 = 0.019202922745021479*tmp_42;
      real_t tmp_79 = tmp_15*(0.039308471900058539*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_80 = tmp_15*(0.039308471900058539*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_81 = tmp_15*(0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_82 = tmp_21*tmp_80 + tmp_27*tmp_81 + tmp_4*tmp_79;
      real_t tmp_83 = tmp_34*tmp_79 + tmp_35*tmp_80 + tmp_36*tmp_81;
      real_t tmp_84 = tmp_38*tmp_79 + tmp_39*tmp_80 + tmp_40*tmp_81;
      real_t tmp_85 = 0.020848748529055869*tmp_42;
      real_t tmp_86 = tmp_15*(0.78764240869137092*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_87 = tmp_15*(0.78764240869137092*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_88 = tmp_15*(0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_89 = tmp_21*tmp_87 + tmp_27*tmp_88 + tmp_4*tmp_86;
      real_t tmp_90 = tmp_34*tmp_86 + tmp_35*tmp_87 + tmp_36*tmp_88;
      real_t tmp_91 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88;
      real_t tmp_92 = 0.019202922745021479*tmp_42;
      real_t tmp_93 = tmp_15*(0.58463275527740355*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_94 = tmp_15*(0.58463275527740355*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_95 = tmp_15*(0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_96 = tmp_21*tmp_94 + tmp_27*tmp_95 + tmp_4*tmp_93;
      real_t tmp_97 = tmp_34*tmp_93 + tmp_35*tmp_94 + tmp_36*tmp_95;
      real_t tmp_98 = tmp_38*tmp_93 + tmp_39*tmp_94 + tmp_40*tmp_95;
      real_t tmp_99 = 0.020848748529055869*tmp_42;
      real_t tmp_100 = tmp_15*(0.1711304259088916*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_101 = tmp_15*(0.1711304259088916*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_102 = tmp_15*(0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_103 = tmp_100*tmp_4 + tmp_101*tmp_21 + tmp_102*tmp_27;
      real_t tmp_104 = tmp_100*tmp_34 + tmp_101*tmp_35 + tmp_102*tmp_36;
      real_t tmp_105 = tmp_100*tmp_38 + tmp_101*tmp_39 + tmp_102*tmp_40;
      real_t tmp_106 = 0.019202922745021479*tmp_42;
      real_t tmp_107 = tmp_15*(0.37605877282253791*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_108 = tmp_15*(0.37605877282253791*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_109 = tmp_15*(0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_110 = tmp_107*tmp_4 + tmp_108*tmp_21 + tmp_109*tmp_27;
      real_t tmp_111 = tmp_107*tmp_34 + tmp_108*tmp_35 + tmp_109*tmp_36;
      real_t tmp_112 = tmp_107*tmp_38 + tmp_108*tmp_39 + tmp_109*tmp_40;
      real_t tmp_113 = 0.020848748529055869*tmp_42;
      real_t tmp_114 = tmp_15*(0.041227165399737475*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_115 = tmp_15*(0.041227165399737475*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_116 = tmp_15*(0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_117 = tmp_114*tmp_4 + tmp_115*tmp_21 + tmp_116*tmp_27;
      real_t tmp_118 = tmp_114*tmp_34 + tmp_115*tmp_35 + tmp_116*tmp_36;
      real_t tmp_119 = tmp_114*tmp_38 + tmp_115*tmp_39 + tmp_116*tmp_40;
      real_t tmp_120 = 0.019202922745021479*tmp_42;
      real_t tmp_121 = tmp_15*(0.40446199974765351*tmp_17 + 0.19107600050469298*tmp_18 + tmp_19);
      real_t tmp_122 = tmp_15*(0.40446199974765351*tmp_23 + 0.19107600050469298*tmp_24 + tmp_25);
      real_t tmp_123 = tmp_15*(0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_30 + tmp_31);
      real_t tmp_124 = tmp_121*tmp_4 + tmp_122*tmp_21 + tmp_123*tmp_27;
      real_t tmp_125 = tmp_121*tmp_34 + tmp_122*tmp_35 + tmp_123*tmp_36;
      real_t tmp_126 = tmp_121*tmp_38 + tmp_122*tmp_39 + tmp_123*tmp_40;
      real_t tmp_127 = 0.042507265838595799*tmp_42;
      real_t tmp_128 = tmp_15*(0.039308471900058539*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_129 = tmp_15*(0.039308471900058539*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_130 = tmp_15*(0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_131 = tmp_128*tmp_4 + tmp_129*tmp_21 + tmp_130*tmp_27;
      real_t tmp_132 = tmp_128*tmp_34 + tmp_129*tmp_35 + tmp_130*tmp_36;
      real_t tmp_133 = tmp_128*tmp_38 + tmp_129*tmp_39 + tmp_130*tmp_40;
      real_t tmp_134 = 0.020848748529055869*tmp_42;
      real_t tmp_135 = tmp_15*(0.93718850182767688*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_136 = tmp_15*(0.93718850182767688*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_137 = tmp_15*(0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_138 = tmp_135*tmp_4 + tmp_136*tmp_21 + tmp_137*tmp_27;
      real_t tmp_139 = tmp_135*tmp_34 + tmp_136*tmp_35 + tmp_137*tmp_36;
      real_t tmp_140 = tmp_135*tmp_38 + tmp_136*tmp_39 + tmp_137*tmp_40;
      real_t tmp_141 = 0.0068572537431980923*tmp_42;
      real_t tmp_142 = tmp_15*(0.60796128279561268*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_143 = tmp_15*(0.60796128279561268*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_144 = tmp_15*(0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_145 = tmp_142*tmp_4 + tmp_143*tmp_21 + tmp_144*tmp_27;
      real_t tmp_146 = tmp_142*tmp_34 + tmp_143*tmp_35 + tmp_144*tmp_36;
      real_t tmp_147 = tmp_142*tmp_38 + tmp_143*tmp_39 + tmp_144*tmp_40;
      real_t tmp_148 = 0.037198804536718075*tmp_42;
      real_t tmp_149 = tmp_15*(0.19107600050469298*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_150 = tmp_15*(0.19107600050469298*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_151 = tmp_15*(0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_152 = tmp_149*tmp_4 + tmp_150*tmp_21 + tmp_151*tmp_27;
      real_t tmp_153 = tmp_149*tmp_34 + tmp_150*tmp_35 + tmp_151*tmp_36;
      real_t tmp_154 = tmp_149*tmp_38 + tmp_150*tmp_39 + tmp_151*tmp_40;
      real_t tmp_155 = 0.042507265838595799*tmp_42;
      real_t tmp_156 = tmp_15*(0.031405749086161582*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_157 = tmp_15*(0.031405749086161582*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_158 = tmp_15*(0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_159 = tmp_156*tmp_4 + tmp_157*tmp_21 + tmp_158*tmp_27;
      real_t tmp_160 = tmp_156*tmp_34 + tmp_157*tmp_35 + tmp_158*tmp_36;
      real_t tmp_161 = tmp_156*tmp_38 + tmp_157*tmp_39 + tmp_158*tmp_40;
      real_t tmp_162 = 0.0068572537431980923*tmp_42;
      real_t tmp_163 = tmp_15*(0.19601935860219369*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_164 = tmp_15*(0.19601935860219369*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_165 = tmp_15*(0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_166 = tmp_163*tmp_4 + tmp_164*tmp_21 + tmp_165*tmp_27;
      real_t tmp_167 = tmp_163*tmp_34 + tmp_164*tmp_35 + tmp_165*tmp_36;
      real_t tmp_168 = tmp_163*tmp_38 + tmp_164*tmp_39 + tmp_165*tmp_40;
      real_t tmp_169 = 0.037198804536718075*tmp_42;
      real_t tmp_170 = tmp_15*(0.40446199974765351*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_171 = tmp_15*(0.40446199974765351*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_172 = tmp_15*(0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_173 = tmp_170*tmp_4 + tmp_171*tmp_21 + tmp_172*tmp_27;
      real_t tmp_174 = tmp_170*tmp_34 + tmp_171*tmp_35 + tmp_172*tmp_36;
      real_t tmp_175 = tmp_170*tmp_38 + tmp_171*tmp_39 + tmp_172*tmp_40;
      real_t tmp_176 = 0.042507265838595799*tmp_42;
      real_t tmp_177 = tmp_15*(0.1711304259088916*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_178 = tmp_15*(0.1711304259088916*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_179 = tmp_15*(0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_180 = tmp_177*tmp_4 + tmp_178*tmp_21 + tmp_179*tmp_27;
      real_t tmp_181 = tmp_177*tmp_34 + tmp_178*tmp_35 + tmp_179*tmp_36;
      real_t tmp_182 = tmp_177*tmp_38 + tmp_178*tmp_39 + tmp_179*tmp_40;
      real_t tmp_183 = 0.019202922745021479*tmp_42;
      real_t a_0_0 = tmp_106*(-tmp_103 - tmp_104 - tmp_105 + 1) + tmp_113*(-tmp_110 - tmp_111 - tmp_112 + 1) + tmp_120*(-tmp_117 - tmp_118 - tmp_119 + 1) + tmp_127*(-tmp_124 - tmp_125 - tmp_126 + 1) + tmp_134*(-tmp_131 - tmp_132 - tmp_133 + 1) + tmp_141*(-tmp_138 - tmp_139 - tmp_140 + 1) + tmp_148*(-tmp_145 - tmp_146 - tmp_147 + 1) + tmp_155*(-tmp_152 - tmp_153 - tmp_154 + 1) + tmp_162*(-tmp_159 - tmp_160 - tmp_161 + 1) + tmp_169*(-tmp_166 - tmp_167 - tmp_168 + 1) + tmp_176*(-tmp_173 - tmp_174 - tmp_175 + 1) + tmp_183*(-tmp_180 - tmp_181 - tmp_182 + 1) + tmp_43*(-tmp_33 - tmp_37 - tmp_41 + 1) + tmp_50*(-tmp_47 - tmp_48 - tmp_49 + 1) + tmp_57*(-tmp_54 - tmp_55 - tmp_56 + 1) + tmp_64*(-tmp_61 - tmp_62 - tmp_63 + 1) + tmp_71*(-tmp_68 - tmp_69 - tmp_70 + 1) + tmp_78*(-tmp_75 - tmp_76 - tmp_77 + 1) + tmp_85*(-tmp_82 - tmp_83 - tmp_84 + 1) + tmp_92*(-tmp_89 - tmp_90 - tmp_91 + 1) + tmp_99*(-tmp_96 - tmp_97 - tmp_98 + 1);
      real_t a_0_1 = tmp_103*tmp_106 + tmp_110*tmp_113 + tmp_117*tmp_120 + tmp_124*tmp_127 + tmp_131*tmp_134 + tmp_138*tmp_141 + tmp_145*tmp_148 + tmp_152*tmp_155 + tmp_159*tmp_162 + tmp_166*tmp_169 + tmp_173*tmp_176 + tmp_180*tmp_183 + tmp_33*tmp_43 + tmp_47*tmp_50 + tmp_54*tmp_57 + tmp_61*tmp_64 + tmp_68*tmp_71 + tmp_75*tmp_78 + tmp_82*tmp_85 + tmp_89*tmp_92 + tmp_96*tmp_99;
      real_t a_0_2 = tmp_104*tmp_106 + tmp_111*tmp_113 + tmp_118*tmp_120 + tmp_125*tmp_127 + tmp_132*tmp_134 + tmp_139*tmp_141 + tmp_146*tmp_148 + tmp_153*tmp_155 + tmp_160*tmp_162 + tmp_167*tmp_169 + tmp_174*tmp_176 + tmp_181*tmp_183 + tmp_37*tmp_43 + tmp_48*tmp_50 + tmp_55*tmp_57 + tmp_62*tmp_64 + tmp_69*tmp_71 + tmp_76*tmp_78 + tmp_83*tmp_85 + tmp_90*tmp_92 + tmp_97*tmp_99;
      real_t a_0_3 = tmp_105*tmp_106 + tmp_112*tmp_113 + tmp_119*tmp_120 + tmp_126*tmp_127 + tmp_133*tmp_134 + tmp_140*tmp_141 + tmp_147*tmp_148 + tmp_154*tmp_155 + tmp_161*tmp_162 + tmp_168*tmp_169 + tmp_175*tmp_176 + tmp_182*tmp_183 + tmp_41*tmp_43 + tmp_49*tmp_50 + tmp_56*tmp_57 + tmp_63*tmp_64 + tmp_70*tmp_71 + tmp_77*tmp_78 + tmp_84*tmp_85 + tmp_91*tmp_92 + tmp_98*tmp_99;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }




void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                                        const std::vector< Point3D >& coordsElementOuter,
                                                        const std::vector< Point3D >& coordsFacet,
                                                        const Point3D&,
                                                        const Point3D&,
                                                        const Point3D&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                  MatrixXr&                            elMat ) const override
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElementInner[0]( 0 );
   const auto p_affine_0_1 = coordsElementInner[0]( 1 );
   const auto p_affine_0_2 = coordsElementInner[0]( 2 );

   const auto p_affine_1_0 = coordsElementInner[1]( 0 );
   const auto p_affine_1_1 = coordsElementInner[1]( 1 );
   const auto p_affine_1_2 = coordsElementInner[1]( 2 );

   const auto p_affine_2_0 = coordsElementInner[2]( 0 );
   const auto p_affine_2_1 = coordsElementInner[2]( 1 );
   const auto p_affine_2_2 = coordsElementInner[2]( 2 );

   const auto p_affine_3_0 = coordsElementInner[3]( 0 );
   const auto p_affine_3_1 = coordsElementInner[3]( 1 );
   const auto p_affine_3_2 = coordsElementInner[3]( 2 );

   const auto p_affine_4_0 = coordsElementOuter[0]( 0 );
   const auto p_affine_4_1 = coordsElementOuter[0]( 1 );
   const auto p_affine_4_2 = coordsElementOuter[0]( 2 );

   const auto p_affine_5_0 = coordsElementOuter[1]( 0 );
   const auto p_affine_5_1 = coordsElementOuter[1]( 1 );
   const auto p_affine_5_2 = coordsElementOuter[1]( 2 );

   const auto p_affine_6_0 = coordsElementOuter[2]( 0 );
   const auto p_affine_6_1 = coordsElementOuter[2]( 1 );
   const auto p_affine_6_2 = coordsElementOuter[2]( 2 );

   const auto p_affine_7_0 = coordsElementOuter[3]( 0 );
   const auto p_affine_7_1 = coordsElementOuter[3]( 1 );
   const auto p_affine_7_2 = coordsElementOuter[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t tmp_0 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_1 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_2 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_3 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_4 = tmp_0*tmp_1 - tmp_2*tmp_3;
      real_t tmp_5 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_6 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_7 = tmp_3*tmp_6;
      real_t tmp_8 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_9 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_10 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_1*tmp_6;
      real_t tmp_13 = tmp_0*tmp_9;
      real_t tmp_14 = tmp_2*tmp_8;
      real_t tmp_15 = 1.0 / (tmp_0*tmp_1*tmp_8 - tmp_10*tmp_12 + tmp_11*tmp_2 - tmp_13*tmp_5 - tmp_14*tmp_3 + tmp_5*tmp_7);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.031405749086161582*tmp_17 + 0.93718850182767688*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_0*tmp_5 + tmp_10*tmp_2;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.031405749086161582*tmp_23 + 0.93718850182767688*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_1*tmp_10 + tmp_3*tmp_5;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_4 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = -tmp_12 + tmp_2*tmp_9;
      real_t tmp_35 = -tmp_14 + tmp_5*tmp_6;
      real_t tmp_36 = tmp_1*tmp_8 - tmp_5*tmp_9;
      real_t tmp_37 = tmp_20*tmp_34 + tmp_26*tmp_35 + tmp_32*tmp_36;
      real_t tmp_38 = -tmp_13 + tmp_7;
      real_t tmp_39 = tmp_0*tmp_8 - tmp_10*tmp_6;
      real_t tmp_40 = tmp_11 - tmp_3*tmp_8;
      real_t tmp_41 = tmp_20*tmp_38 + tmp_26*tmp_39 + tmp_32*tmp_40;
      real_t tmp_42 = 0.5*p_affine_13_2*std::pow((std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)), 1.0/2.0);
      real_t tmp_43 = 0.0068572537431980923*tmp_42;
      real_t tmp_44 = tmp_15*(0.19601935860219369*tmp_17 + 0.60796128279561268*tmp_18 + tmp_19);
      real_t tmp_45 = tmp_15*(0.19601935860219369*tmp_23 + 0.60796128279561268*tmp_24 + tmp_25);
      real_t tmp_46 = tmp_15*(0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_30 + tmp_31);
      real_t tmp_47 = tmp_21*tmp_45 + tmp_27*tmp_46 + tmp_4*tmp_44;
      real_t tmp_48 = tmp_34*tmp_44 + tmp_35*tmp_45 + tmp_36*tmp_46;
      real_t tmp_49 = tmp_38*tmp_44 + tmp_39*tmp_45 + tmp_40*tmp_46;
      real_t tmp_50 = 0.037198804536718075*tmp_42;
      real_t tmp_51 = tmp_15*(0.37605877282253791*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_52 = tmp_15*(0.37605877282253791*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_53 = tmp_15*(0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_54 = tmp_21*tmp_52 + tmp_27*tmp_53 + tmp_4*tmp_51;
      real_t tmp_55 = tmp_34*tmp_51 + tmp_35*tmp_52 + tmp_36*tmp_53;
      real_t tmp_56 = tmp_38*tmp_51 + tmp_39*tmp_52 + tmp_40*tmp_53;
      real_t tmp_57 = 0.020848748529055869*tmp_42;
      real_t tmp_58 = tmp_15*(0.78764240869137092*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_59 = tmp_15*(0.78764240869137092*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_60 = tmp_15*(0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_61 = tmp_21*tmp_59 + tmp_27*tmp_60 + tmp_4*tmp_58;
      real_t tmp_62 = tmp_34*tmp_58 + tmp_35*tmp_59 + tmp_36*tmp_60;
      real_t tmp_63 = tmp_38*tmp_58 + tmp_39*tmp_59 + tmp_40*tmp_60;
      real_t tmp_64 = 0.019202922745021479*tmp_42;
      real_t tmp_65 = tmp_15*(0.58463275527740355*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_66 = tmp_15*(0.58463275527740355*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_67 = tmp_15*(0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_68 = tmp_21*tmp_66 + tmp_27*tmp_67 + tmp_4*tmp_65;
      real_t tmp_69 = tmp_34*tmp_65 + tmp_35*tmp_66 + tmp_36*tmp_67;
      real_t tmp_70 = tmp_38*tmp_65 + tmp_39*tmp_66 + tmp_40*tmp_67;
      real_t tmp_71 = 0.020848748529055869*tmp_42;
      real_t tmp_72 = tmp_15*(0.041227165399737475*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_73 = tmp_15*(0.041227165399737475*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_74 = tmp_15*(0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_75 = tmp_21*tmp_73 + tmp_27*tmp_74 + tmp_4*tmp_72;
      real_t tmp_76 = tmp_34*tmp_72 + tmp_35*tmp_73 + tmp_36*tmp_74;
      real_t tmp_77 = tmp_38*tmp_72 + tmp_39*tmp_73 + tmp_40*tmp_74;
      real_t tmp_78 = 0.019202922745021479*tmp_42;
      real_t tmp_79 = tmp_15*(0.039308471900058539*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_80 = tmp_15*(0.039308471900058539*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_81 = tmp_15*(0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_82 = tmp_21*tmp_80 + tmp_27*tmp_81 + tmp_4*tmp_79;
      real_t tmp_83 = tmp_34*tmp_79 + tmp_35*tmp_80 + tmp_36*tmp_81;
      real_t tmp_84 = tmp_38*tmp_79 + tmp_39*tmp_80 + tmp_40*tmp_81;
      real_t tmp_85 = 0.020848748529055869*tmp_42;
      real_t tmp_86 = tmp_15*(0.78764240869137092*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_87 = tmp_15*(0.78764240869137092*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_88 = tmp_15*(0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_89 = tmp_21*tmp_87 + tmp_27*tmp_88 + tmp_4*tmp_86;
      real_t tmp_90 = tmp_34*tmp_86 + tmp_35*tmp_87 + tmp_36*tmp_88;
      real_t tmp_91 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88;
      real_t tmp_92 = 0.019202922745021479*tmp_42;
      real_t tmp_93 = tmp_15*(0.58463275527740355*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_94 = tmp_15*(0.58463275527740355*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_95 = tmp_15*(0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_96 = tmp_21*tmp_94 + tmp_27*tmp_95 + tmp_4*tmp_93;
      real_t tmp_97 = tmp_34*tmp_93 + tmp_35*tmp_94 + tmp_36*tmp_95;
      real_t tmp_98 = tmp_38*tmp_93 + tmp_39*tmp_94 + tmp_40*tmp_95;
      real_t tmp_99 = 0.020848748529055869*tmp_42;
      real_t tmp_100 = tmp_15*(0.1711304259088916*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_101 = tmp_15*(0.1711304259088916*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_102 = tmp_15*(0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_103 = tmp_100*tmp_4 + tmp_101*tmp_21 + tmp_102*tmp_27;
      real_t tmp_104 = tmp_100*tmp_34 + tmp_101*tmp_35 + tmp_102*tmp_36;
      real_t tmp_105 = tmp_100*tmp_38 + tmp_101*tmp_39 + tmp_102*tmp_40;
      real_t tmp_106 = 0.019202922745021479*tmp_42;
      real_t tmp_107 = tmp_15*(0.37605877282253791*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_108 = tmp_15*(0.37605877282253791*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_109 = tmp_15*(0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_110 = tmp_107*tmp_4 + tmp_108*tmp_21 + tmp_109*tmp_27;
      real_t tmp_111 = tmp_107*tmp_34 + tmp_108*tmp_35 + tmp_109*tmp_36;
      real_t tmp_112 = tmp_107*tmp_38 + tmp_108*tmp_39 + tmp_109*tmp_40;
      real_t tmp_113 = 0.020848748529055869*tmp_42;
      real_t tmp_114 = tmp_15*(0.041227165399737475*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_115 = tmp_15*(0.041227165399737475*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_116 = tmp_15*(0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_117 = tmp_114*tmp_4 + tmp_115*tmp_21 + tmp_116*tmp_27;
      real_t tmp_118 = tmp_114*tmp_34 + tmp_115*tmp_35 + tmp_116*tmp_36;
      real_t tmp_119 = tmp_114*tmp_38 + tmp_115*tmp_39 + tmp_116*tmp_40;
      real_t tmp_120 = 0.019202922745021479*tmp_42;
      real_t tmp_121 = tmp_15*(0.40446199974765351*tmp_17 + 0.19107600050469298*tmp_18 + tmp_19);
      real_t tmp_122 = tmp_15*(0.40446199974765351*tmp_23 + 0.19107600050469298*tmp_24 + tmp_25);
      real_t tmp_123 = tmp_15*(0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_30 + tmp_31);
      real_t tmp_124 = tmp_121*tmp_4 + tmp_122*tmp_21 + tmp_123*tmp_27;
      real_t tmp_125 = tmp_121*tmp_34 + tmp_122*tmp_35 + tmp_123*tmp_36;
      real_t tmp_126 = tmp_121*tmp_38 + tmp_122*tmp_39 + tmp_123*tmp_40;
      real_t tmp_127 = 0.042507265838595799*tmp_42;
      real_t tmp_128 = tmp_15*(0.039308471900058539*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_129 = tmp_15*(0.039308471900058539*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_130 = tmp_15*(0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_131 = tmp_128*tmp_4 + tmp_129*tmp_21 + tmp_130*tmp_27;
      real_t tmp_132 = tmp_128*tmp_34 + tmp_129*tmp_35 + tmp_130*tmp_36;
      real_t tmp_133 = tmp_128*tmp_38 + tmp_129*tmp_39 + tmp_130*tmp_40;
      real_t tmp_134 = 0.020848748529055869*tmp_42;
      real_t tmp_135 = tmp_15*(0.93718850182767688*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_136 = tmp_15*(0.93718850182767688*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_137 = tmp_15*(0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_138 = tmp_135*tmp_4 + tmp_136*tmp_21 + tmp_137*tmp_27;
      real_t tmp_139 = tmp_135*tmp_34 + tmp_136*tmp_35 + tmp_137*tmp_36;
      real_t tmp_140 = tmp_135*tmp_38 + tmp_136*tmp_39 + tmp_137*tmp_40;
      real_t tmp_141 = 0.0068572537431980923*tmp_42;
      real_t tmp_142 = tmp_15*(0.60796128279561268*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_143 = tmp_15*(0.60796128279561268*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_144 = tmp_15*(0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_145 = tmp_142*tmp_4 + tmp_143*tmp_21 + tmp_144*tmp_27;
      real_t tmp_146 = tmp_142*tmp_34 + tmp_143*tmp_35 + tmp_144*tmp_36;
      real_t tmp_147 = tmp_142*tmp_38 + tmp_143*tmp_39 + tmp_144*tmp_40;
      real_t tmp_148 = 0.037198804536718075*tmp_42;
      real_t tmp_149 = tmp_15*(0.19107600050469298*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_150 = tmp_15*(0.19107600050469298*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_151 = tmp_15*(0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_152 = tmp_149*tmp_4 + tmp_150*tmp_21 + tmp_151*tmp_27;
      real_t tmp_153 = tmp_149*tmp_34 + tmp_150*tmp_35 + tmp_151*tmp_36;
      real_t tmp_154 = tmp_149*tmp_38 + tmp_150*tmp_39 + tmp_151*tmp_40;
      real_t tmp_155 = 0.042507265838595799*tmp_42;
      real_t tmp_156 = tmp_15*(0.031405749086161582*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_157 = tmp_15*(0.031405749086161582*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_158 = tmp_15*(0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_159 = tmp_156*tmp_4 + tmp_157*tmp_21 + tmp_158*tmp_27;
      real_t tmp_160 = tmp_156*tmp_34 + tmp_157*tmp_35 + tmp_158*tmp_36;
      real_t tmp_161 = tmp_156*tmp_38 + tmp_157*tmp_39 + tmp_158*tmp_40;
      real_t tmp_162 = 0.0068572537431980923*tmp_42;
      real_t tmp_163 = tmp_15*(0.19601935860219369*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_164 = tmp_15*(0.19601935860219369*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_165 = tmp_15*(0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_166 = tmp_163*tmp_4 + tmp_164*tmp_21 + tmp_165*tmp_27;
      real_t tmp_167 = tmp_163*tmp_34 + tmp_164*tmp_35 + tmp_165*tmp_36;
      real_t tmp_168 = tmp_163*tmp_38 + tmp_164*tmp_39 + tmp_165*tmp_40;
      real_t tmp_169 = 0.037198804536718075*tmp_42;
      real_t tmp_170 = tmp_15*(0.40446199974765351*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_171 = tmp_15*(0.40446199974765351*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_172 = tmp_15*(0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_173 = tmp_170*tmp_4 + tmp_171*tmp_21 + tmp_172*tmp_27;
      real_t tmp_174 = tmp_170*tmp_34 + tmp_171*tmp_35 + tmp_172*tmp_36;
      real_t tmp_175 = tmp_170*tmp_38 + tmp_171*tmp_39 + tmp_172*tmp_40;
      real_t tmp_176 = 0.042507265838595799*tmp_42;
      real_t tmp_177 = tmp_15*(0.1711304259088916*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_178 = tmp_15*(0.1711304259088916*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_179 = tmp_15*(0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_180 = tmp_177*tmp_4 + tmp_178*tmp_21 + tmp_179*tmp_27;
      real_t tmp_181 = tmp_177*tmp_34 + tmp_178*tmp_35 + tmp_179*tmp_36;
      real_t tmp_182 = tmp_177*tmp_38 + tmp_178*tmp_39 + tmp_179*tmp_40;
      real_t tmp_183 = 0.019202922745021479*tmp_42;
      real_t a_0_0 = -tmp_106*(-tmp_103 - tmp_104 - tmp_105 + 1) - tmp_113*(-tmp_110 - tmp_111 - tmp_112 + 1) - tmp_120*(-tmp_117 - tmp_118 - tmp_119 + 1) - tmp_127*(-tmp_124 - tmp_125 - tmp_126 + 1) - tmp_134*(-tmp_131 - tmp_132 - tmp_133 + 1) - tmp_141*(-tmp_138 - tmp_139 - tmp_140 + 1) - tmp_148*(-tmp_145 - tmp_146 - tmp_147 + 1) - tmp_155*(-tmp_152 - tmp_153 - tmp_154 + 1) - tmp_162*(-tmp_159 - tmp_160 - tmp_161 + 1) - tmp_169*(-tmp_166 - tmp_167 - tmp_168 + 1) - tmp_176*(-tmp_173 - tmp_174 - tmp_175 + 1) - tmp_183*(-tmp_180 - tmp_181 - tmp_182 + 1) - tmp_43*(-tmp_33 - tmp_37 - tmp_41 + 1) - tmp_50*(-tmp_47 - tmp_48 - tmp_49 + 1) - tmp_57*(-tmp_54 - tmp_55 - tmp_56 + 1) - tmp_64*(-tmp_61 - tmp_62 - tmp_63 + 1) - tmp_71*(-tmp_68 - tmp_69 - tmp_70 + 1) - tmp_78*(-tmp_75 - tmp_76 - tmp_77 + 1) - tmp_85*(-tmp_82 - tmp_83 - tmp_84 + 1) - tmp_92*(-tmp_89 - tmp_90 - tmp_91 + 1) - tmp_99*(-tmp_96 - tmp_97 - tmp_98 + 1);
      real_t a_0_1 = -tmp_103*tmp_106 - tmp_110*tmp_113 - tmp_117*tmp_120 - tmp_124*tmp_127 - tmp_131*tmp_134 - tmp_138*tmp_141 - tmp_145*tmp_148 - tmp_152*tmp_155 - tmp_159*tmp_162 - tmp_166*tmp_169 - tmp_173*tmp_176 - tmp_180*tmp_183 - tmp_33*tmp_43 - tmp_47*tmp_50 - tmp_54*tmp_57 - tmp_61*tmp_64 - tmp_68*tmp_71 - tmp_75*tmp_78 - tmp_82*tmp_85 - tmp_89*tmp_92 - tmp_96*tmp_99;
      real_t a_0_2 = -tmp_104*tmp_106 - tmp_111*tmp_113 - tmp_118*tmp_120 - tmp_125*tmp_127 - tmp_132*tmp_134 - tmp_139*tmp_141 - tmp_146*tmp_148 - tmp_153*tmp_155 - tmp_160*tmp_162 - tmp_167*tmp_169 - tmp_174*tmp_176 - tmp_181*tmp_183 - tmp_37*tmp_43 - tmp_48*tmp_50 - tmp_55*tmp_57 - tmp_62*tmp_64 - tmp_69*tmp_71 - tmp_76*tmp_78 - tmp_83*tmp_85 - tmp_90*tmp_92 - tmp_97*tmp_99;
      real_t a_0_3 = -tmp_105*tmp_106 - tmp_112*tmp_113 - tmp_119*tmp_120 - tmp_126*tmp_127 - tmp_133*tmp_134 - tmp_140*tmp_141 - tmp_147*tmp_148 - tmp_154*tmp_155 - tmp_161*tmp_162 - tmp_168*tmp_169 - tmp_175*tmp_176 - tmp_182*tmp_183 - tmp_41*tmp_43 - tmp_49*tmp_50 - tmp_56*tmp_57 - tmp_63*tmp_64 - tmp_70*tmp_71 - tmp_77*tmp_78 - tmp_84*tmp_85 - tmp_91*tmp_92 - tmp_98*tmp_99;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Point3D >& coordsElement,
    const std::vector< Point3D >& coordsFacet,
    const Point3D&,
    const Point3D&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
                                        MatrixXr&                            elMat ) const override
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t tmp_0 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_2 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_3 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_4 = tmp_0*tmp_1 - tmp_2*tmp_3;
      real_t tmp_5 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_6 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_7 = tmp_3*tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_9 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_10 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = tmp_1*tmp_6;
      real_t tmp_13 = tmp_0*tmp_9;
      real_t tmp_14 = tmp_2*tmp_8;
      real_t tmp_15 = 1.0 / (tmp_0*tmp_1*tmp_8 - tmp_10*tmp_12 + tmp_11*tmp_2 - tmp_13*tmp_5 - tmp_14*tmp_3 + tmp_5*tmp_7);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.031405749086161582*tmp_17 + 0.93718850182767688*tmp_18 + tmp_19);
      real_t tmp_21 = -tmp_0*tmp_5 + tmp_10*tmp_2;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.031405749086161582*tmp_23 + 0.93718850182767688*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_1*tmp_10 + tmp_3*tmp_5;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_4 + tmp_21*tmp_26 + tmp_27*tmp_32;
      real_t tmp_34 = -tmp_12 + tmp_2*tmp_9;
      real_t tmp_35 = -tmp_14 + tmp_5*tmp_6;
      real_t tmp_36 = tmp_1*tmp_8 - tmp_5*tmp_9;
      real_t tmp_37 = tmp_20*tmp_34 + tmp_26*tmp_35 + tmp_32*tmp_36;
      real_t tmp_38 = -tmp_13 + tmp_7;
      real_t tmp_39 = tmp_0*tmp_8 - tmp_10*tmp_6;
      real_t tmp_40 = tmp_11 - tmp_3*tmp_8;
      real_t tmp_41 = tmp_20*tmp_38 + tmp_26*tmp_39 + tmp_32*tmp_40;
      real_t tmp_42 = 1.0*p_affine_13_2*std::pow((std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)), 1.0/2.0);
      real_t tmp_43 = 0.0068572537431980923*tmp_42;
      real_t tmp_44 = tmp_15*(0.19601935860219369*tmp_17 + 0.60796128279561268*tmp_18 + tmp_19);
      real_t tmp_45 = tmp_15*(0.19601935860219369*tmp_23 + 0.60796128279561268*tmp_24 + tmp_25);
      real_t tmp_46 = tmp_15*(0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_30 + tmp_31);
      real_t tmp_47 = tmp_21*tmp_45 + tmp_27*tmp_46 + tmp_4*tmp_44;
      real_t tmp_48 = tmp_34*tmp_44 + tmp_35*tmp_45 + tmp_36*tmp_46;
      real_t tmp_49 = tmp_38*tmp_44 + tmp_39*tmp_45 + tmp_40*tmp_46;
      real_t tmp_50 = 0.037198804536718075*tmp_42;
      real_t tmp_51 = tmp_15*(0.37605877282253791*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_52 = tmp_15*(0.37605877282253791*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_53 = tmp_15*(0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_54 = tmp_21*tmp_52 + tmp_27*tmp_53 + tmp_4*tmp_51;
      real_t tmp_55 = tmp_34*tmp_51 + tmp_35*tmp_52 + tmp_36*tmp_53;
      real_t tmp_56 = tmp_38*tmp_51 + tmp_39*tmp_52 + tmp_40*tmp_53;
      real_t tmp_57 = 0.020848748529055869*tmp_42;
      real_t tmp_58 = tmp_15*(0.78764240869137092*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_59 = tmp_15*(0.78764240869137092*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_60 = tmp_15*(0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_61 = tmp_21*tmp_59 + tmp_27*tmp_60 + tmp_4*tmp_58;
      real_t tmp_62 = tmp_34*tmp_58 + tmp_35*tmp_59 + tmp_36*tmp_60;
      real_t tmp_63 = tmp_38*tmp_58 + tmp_39*tmp_59 + tmp_40*tmp_60;
      real_t tmp_64 = 0.019202922745021479*tmp_42;
      real_t tmp_65 = tmp_15*(0.58463275527740355*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_66 = tmp_15*(0.58463275527740355*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_67 = tmp_15*(0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_68 = tmp_21*tmp_66 + tmp_27*tmp_67 + tmp_4*tmp_65;
      real_t tmp_69 = tmp_34*tmp_65 + tmp_35*tmp_66 + tmp_36*tmp_67;
      real_t tmp_70 = tmp_38*tmp_65 + tmp_39*tmp_66 + tmp_40*tmp_67;
      real_t tmp_71 = 0.020848748529055869*tmp_42;
      real_t tmp_72 = tmp_15*(0.041227165399737475*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_73 = tmp_15*(0.041227165399737475*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_74 = tmp_15*(0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_75 = tmp_21*tmp_73 + tmp_27*tmp_74 + tmp_4*tmp_72;
      real_t tmp_76 = tmp_34*tmp_72 + tmp_35*tmp_73 + tmp_36*tmp_74;
      real_t tmp_77 = tmp_38*tmp_72 + tmp_39*tmp_73 + tmp_40*tmp_74;
      real_t tmp_78 = 0.019202922745021479*tmp_42;
      real_t tmp_79 = tmp_15*(0.039308471900058539*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_80 = tmp_15*(0.039308471900058539*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_81 = tmp_15*(0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_82 = tmp_21*tmp_80 + tmp_27*tmp_81 + tmp_4*tmp_79;
      real_t tmp_83 = tmp_34*tmp_79 + tmp_35*tmp_80 + tmp_36*tmp_81;
      real_t tmp_84 = tmp_38*tmp_79 + tmp_39*tmp_80 + tmp_40*tmp_81;
      real_t tmp_85 = 0.020848748529055869*tmp_42;
      real_t tmp_86 = tmp_15*(0.78764240869137092*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_87 = tmp_15*(0.78764240869137092*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_88 = tmp_15*(0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_89 = tmp_21*tmp_87 + tmp_27*tmp_88 + tmp_4*tmp_86;
      real_t tmp_90 = tmp_34*tmp_86 + tmp_35*tmp_87 + tmp_36*tmp_88;
      real_t tmp_91 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88;
      real_t tmp_92 = 0.019202922745021479*tmp_42;
      real_t tmp_93 = tmp_15*(0.58463275527740355*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_94 = tmp_15*(0.58463275527740355*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_95 = tmp_15*(0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_96 = tmp_21*tmp_94 + tmp_27*tmp_95 + tmp_4*tmp_93;
      real_t tmp_97 = tmp_34*tmp_93 + tmp_35*tmp_94 + tmp_36*tmp_95;
      real_t tmp_98 = tmp_38*tmp_93 + tmp_39*tmp_94 + tmp_40*tmp_95;
      real_t tmp_99 = 0.020848748529055869*tmp_42;
      real_t tmp_100 = tmp_15*(0.1711304259088916*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_101 = tmp_15*(0.1711304259088916*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_102 = tmp_15*(0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_103 = tmp_100*tmp_4 + tmp_101*tmp_21 + tmp_102*tmp_27;
      real_t tmp_104 = tmp_100*tmp_34 + tmp_101*tmp_35 + tmp_102*tmp_36;
      real_t tmp_105 = tmp_100*tmp_38 + tmp_101*tmp_39 + tmp_102*tmp_40;
      real_t tmp_106 = 0.019202922745021479*tmp_42;
      real_t tmp_107 = tmp_15*(0.37605877282253791*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_108 = tmp_15*(0.37605877282253791*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_109 = tmp_15*(0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_110 = tmp_107*tmp_4 + tmp_108*tmp_21 + tmp_109*tmp_27;
      real_t tmp_111 = tmp_107*tmp_34 + tmp_108*tmp_35 + tmp_109*tmp_36;
      real_t tmp_112 = tmp_107*tmp_38 + tmp_108*tmp_39 + tmp_109*tmp_40;
      real_t tmp_113 = 0.020848748529055869*tmp_42;
      real_t tmp_114 = tmp_15*(0.041227165399737475*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_115 = tmp_15*(0.041227165399737475*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_116 = tmp_15*(0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_117 = tmp_114*tmp_4 + tmp_115*tmp_21 + tmp_116*tmp_27;
      real_t tmp_118 = tmp_114*tmp_34 + tmp_115*tmp_35 + tmp_116*tmp_36;
      real_t tmp_119 = tmp_114*tmp_38 + tmp_115*tmp_39 + tmp_116*tmp_40;
      real_t tmp_120 = 0.019202922745021479*tmp_42;
      real_t tmp_121 = tmp_15*(0.40446199974765351*tmp_17 + 0.19107600050469298*tmp_18 + tmp_19);
      real_t tmp_122 = tmp_15*(0.40446199974765351*tmp_23 + 0.19107600050469298*tmp_24 + tmp_25);
      real_t tmp_123 = tmp_15*(0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_30 + tmp_31);
      real_t tmp_124 = tmp_121*tmp_4 + tmp_122*tmp_21 + tmp_123*tmp_27;
      real_t tmp_125 = tmp_121*tmp_34 + tmp_122*tmp_35 + tmp_123*tmp_36;
      real_t tmp_126 = tmp_121*tmp_38 + tmp_122*tmp_39 + tmp_123*tmp_40;
      real_t tmp_127 = 0.042507265838595799*tmp_42;
      real_t tmp_128 = tmp_15*(0.039308471900058539*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_129 = tmp_15*(0.039308471900058539*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_130 = tmp_15*(0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_131 = tmp_128*tmp_4 + tmp_129*tmp_21 + tmp_130*tmp_27;
      real_t tmp_132 = tmp_128*tmp_34 + tmp_129*tmp_35 + tmp_130*tmp_36;
      real_t tmp_133 = tmp_128*tmp_38 + tmp_129*tmp_39 + tmp_130*tmp_40;
      real_t tmp_134 = 0.020848748529055869*tmp_42;
      real_t tmp_135 = tmp_15*(0.93718850182767688*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_136 = tmp_15*(0.93718850182767688*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_137 = tmp_15*(0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_138 = tmp_135*tmp_4 + tmp_136*tmp_21 + tmp_137*tmp_27;
      real_t tmp_139 = tmp_135*tmp_34 + tmp_136*tmp_35 + tmp_137*tmp_36;
      real_t tmp_140 = tmp_135*tmp_38 + tmp_136*tmp_39 + tmp_137*tmp_40;
      real_t tmp_141 = 0.0068572537431980923*tmp_42;
      real_t tmp_142 = tmp_15*(0.60796128279561268*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_143 = tmp_15*(0.60796128279561268*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_144 = tmp_15*(0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_145 = tmp_142*tmp_4 + tmp_143*tmp_21 + tmp_144*tmp_27;
      real_t tmp_146 = tmp_142*tmp_34 + tmp_143*tmp_35 + tmp_144*tmp_36;
      real_t tmp_147 = tmp_142*tmp_38 + tmp_143*tmp_39 + tmp_144*tmp_40;
      real_t tmp_148 = 0.037198804536718075*tmp_42;
      real_t tmp_149 = tmp_15*(0.19107600050469298*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_150 = tmp_15*(0.19107600050469298*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_151 = tmp_15*(0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_152 = tmp_149*tmp_4 + tmp_150*tmp_21 + tmp_151*tmp_27;
      real_t tmp_153 = tmp_149*tmp_34 + tmp_150*tmp_35 + tmp_151*tmp_36;
      real_t tmp_154 = tmp_149*tmp_38 + tmp_150*tmp_39 + tmp_151*tmp_40;
      real_t tmp_155 = 0.042507265838595799*tmp_42;
      real_t tmp_156 = tmp_15*(0.031405749086161582*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_157 = tmp_15*(0.031405749086161582*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_158 = tmp_15*(0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_159 = tmp_156*tmp_4 + tmp_157*tmp_21 + tmp_158*tmp_27;
      real_t tmp_160 = tmp_156*tmp_34 + tmp_157*tmp_35 + tmp_158*tmp_36;
      real_t tmp_161 = tmp_156*tmp_38 + tmp_157*tmp_39 + tmp_158*tmp_40;
      real_t tmp_162 = 0.0068572537431980923*tmp_42;
      real_t tmp_163 = tmp_15*(0.19601935860219369*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_164 = tmp_15*(0.19601935860219369*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_165 = tmp_15*(0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_166 = tmp_163*tmp_4 + tmp_164*tmp_21 + tmp_165*tmp_27;
      real_t tmp_167 = tmp_163*tmp_34 + tmp_164*tmp_35 + tmp_165*tmp_36;
      real_t tmp_168 = tmp_163*tmp_38 + tmp_164*tmp_39 + tmp_165*tmp_40;
      real_t tmp_169 = 0.037198804536718075*tmp_42;
      real_t tmp_170 = tmp_15*(0.40446199974765351*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_171 = tmp_15*(0.40446199974765351*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_172 = tmp_15*(0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_173 = tmp_170*tmp_4 + tmp_171*tmp_21 + tmp_172*tmp_27;
      real_t tmp_174 = tmp_170*tmp_34 + tmp_171*tmp_35 + tmp_172*tmp_36;
      real_t tmp_175 = tmp_170*tmp_38 + tmp_171*tmp_39 + tmp_172*tmp_40;
      real_t tmp_176 = 0.042507265838595799*tmp_42;
      real_t tmp_177 = tmp_15*(0.1711304259088916*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_178 = tmp_15*(0.1711304259088916*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_179 = tmp_15*(0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_180 = tmp_177*tmp_4 + tmp_178*tmp_21 + tmp_179*tmp_27;
      real_t tmp_181 = tmp_177*tmp_34 + tmp_178*tmp_35 + tmp_179*tmp_36;
      real_t tmp_182 = tmp_177*tmp_38 + tmp_178*tmp_39 + tmp_179*tmp_40;
      real_t tmp_183 = 0.019202922745021479*tmp_42;
      real_t a_0_0 = tmp_106*(-tmp_103 - tmp_104 - tmp_105 + 1) + tmp_113*(-tmp_110 - tmp_111 - tmp_112 + 1) + tmp_120*(-tmp_117 - tmp_118 - tmp_119 + 1) + tmp_127*(-tmp_124 - tmp_125 - tmp_126 + 1) + tmp_134*(-tmp_131 - tmp_132 - tmp_133 + 1) + tmp_141*(-tmp_138 - tmp_139 - tmp_140 + 1) + tmp_148*(-tmp_145 - tmp_146 - tmp_147 + 1) + tmp_155*(-tmp_152 - tmp_153 - tmp_154 + 1) + tmp_162*(-tmp_159 - tmp_160 - tmp_161 + 1) + tmp_169*(-tmp_166 - tmp_167 - tmp_168 + 1) + tmp_176*(-tmp_173 - tmp_174 - tmp_175 + 1) + tmp_183*(-tmp_180 - tmp_181 - tmp_182 + 1) + tmp_43*(-tmp_33 - tmp_37 - tmp_41 + 1) + tmp_50*(-tmp_47 - tmp_48 - tmp_49 + 1) + tmp_57*(-tmp_54 - tmp_55 - tmp_56 + 1) + tmp_64*(-tmp_61 - tmp_62 - tmp_63 + 1) + tmp_71*(-tmp_68 - tmp_69 - tmp_70 + 1) + tmp_78*(-tmp_75 - tmp_76 - tmp_77 + 1) + tmp_85*(-tmp_82 - tmp_83 - tmp_84 + 1) + tmp_92*(-tmp_89 - tmp_90 - tmp_91 + 1) + tmp_99*(-tmp_96 - tmp_97 - tmp_98 + 1);
      real_t a_0_1 = tmp_103*tmp_106 + tmp_110*tmp_113 + tmp_117*tmp_120 + tmp_124*tmp_127 + tmp_131*tmp_134 + tmp_138*tmp_141 + tmp_145*tmp_148 + tmp_152*tmp_155 + tmp_159*tmp_162 + tmp_166*tmp_169 + tmp_173*tmp_176 + tmp_180*tmp_183 + tmp_33*tmp_43 + tmp_47*tmp_50 + tmp_54*tmp_57 + tmp_61*tmp_64 + tmp_68*tmp_71 + tmp_75*tmp_78 + tmp_82*tmp_85 + tmp_89*tmp_92 + tmp_96*tmp_99;
      real_t a_0_2 = tmp_104*tmp_106 + tmp_111*tmp_113 + tmp_118*tmp_120 + tmp_125*tmp_127 + tmp_132*tmp_134 + tmp_139*tmp_141 + tmp_146*tmp_148 + tmp_153*tmp_155 + tmp_160*tmp_162 + tmp_167*tmp_169 + tmp_174*tmp_176 + tmp_181*tmp_183 + tmp_37*tmp_43 + tmp_48*tmp_50 + tmp_55*tmp_57 + tmp_62*tmp_64 + tmp_69*tmp_71 + tmp_76*tmp_78 + tmp_83*tmp_85 + tmp_90*tmp_92 + tmp_97*tmp_99;
      real_t a_0_3 = tmp_105*tmp_106 + tmp_112*tmp_113 + tmp_119*tmp_120 + tmp_126*tmp_127 + tmp_133*tmp_134 + tmp_140*tmp_141 + tmp_147*tmp_148 + tmp_154*tmp_155 + tmp_161*tmp_162 + tmp_168*tmp_169 + tmp_175*tmp_176 + tmp_182*tmp_183 + tmp_41*tmp_43 + tmp_49*tmp_50 + tmp_56*tmp_57 + tmp_63*tmp_64 + tmp_70*tmp_71 + tmp_77*tmp_78 + tmp_84*tmp_85 + tmp_91*tmp_92 + tmp_98*tmp_99;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }

public:

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g1;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g1;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g0;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g0;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g2;

};




class EGDivFormNitscheBC_P0E : public hyteg::dg::DGForm
{

 public:
    EGDivFormNitscheBC_P0E()
: callback_Scalar_Variable_Coefficient_3D_g1 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_2D_g1 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_2D_g0 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_3D_g0 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_3D_g2 ([](const Point3D & p) -> real_t { return 0.; })
    {}

void Scalar_Variable_Coefficient_2D_g0( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g0( Point3D( {in_0, in_1, 0} ) );
}
void Scalar_Variable_Coefficient_2D_g1( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g1( Point3D( {in_0, in_1, 0} ) );
}

void Scalar_Variable_Coefficient_3D_g0( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g0( Point3D( {in_0, in_1, in_2} ) );
}
void Scalar_Variable_Coefficient_3D_g1( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g1( Point3D( {in_0, in_1, in_2} ) );
}
void Scalar_Variable_Coefficient_3D_g2( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_3D_g2( Point3D( {in_0, in_1, in_2} ) );
}

 protected:
  void integrateVolume2D( const std::vector< Point3D >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t tmp_0 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_1 = -tmp_0;
      real_t tmp_2 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_3 = (-p_affine_0_0 + p_affine_1_0)*(-p_affine_0_1 + p_affine_2_1);
      real_t tmp_4 = -tmp_2;
      real_t tmp_5 = 1.0 / (-tmp_1*tmp_4 + tmp_3);
      real_t tmp_6 = (-tmp_0*tmp_4*tmp_5 - tmp_1*tmp_2*tmp_5 - 2*tmp_3*tmp_5)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = 0.5*tmp_6;
      elMat( 0, 0) = a_0_0;
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >&      coordsElement,
                                       const std::vector< Point3D >&      coordsFacet,
                                       const Point3D&                     oppositeVertex,
                                       const Point3D&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_4 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_5 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_6 = -tmp_4;
      real_t tmp_7 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_8 = -tmp_7;
      real_t tmp_9 = 1.0 / (tmp_3*tmp_5 - tmp_6*tmp_8);
      real_t tmp_10 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_11 = tmp_9*(0.046910077030668018*tmp_1 + tmp_10);
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = tmp_9*(0.046910077030668018*tmp_0 + tmp_12);
      real_t tmp_14 = tmp_11*tmp_4 + tmp_13*tmp_5 - 1.0/3.0;
      real_t tmp_15 = tmp_11*tmp_3 + tmp_13*tmp_7 - 1.0/3.0;
      real_t tmp_16 = 0.5*p_affine_10_0;
      real_t tmp_17 = 0.5*p_affine_10_1;
      real_t tmp_18 = tmp_9*(0.23076534494715845*tmp_1 + tmp_10);
      real_t tmp_19 = tmp_9*(0.23076534494715845*tmp_0 + tmp_12);
      real_t tmp_20 = tmp_18*tmp_4 + tmp_19*tmp_5 - 1.0/3.0;
      real_t tmp_21 = tmp_18*tmp_3 + tmp_19*tmp_7 - 1.0/3.0;
      real_t tmp_22 = tmp_9*(0.5*tmp_1 + tmp_10);
      real_t tmp_23 = tmp_9*(0.5*tmp_0 + tmp_12);
      real_t tmp_24 = tmp_22*tmp_4 + tmp_23*tmp_5 - 1.0/3.0;
      real_t tmp_25 = tmp_22*tmp_3 + tmp_23*tmp_7 - 1.0/3.0;
      real_t tmp_26 = tmp_9*(0.7692346550528415*tmp_1 + tmp_10);
      real_t tmp_27 = tmp_9*(0.7692346550528415*tmp_0 + tmp_12);
      real_t tmp_28 = tmp_26*tmp_4 + tmp_27*tmp_5 - 1.0/3.0;
      real_t tmp_29 = tmp_26*tmp_3 + tmp_27*tmp_7 - 1.0/3.0;
      real_t tmp_30 = tmp_9*(0.95308992296933193*tmp_1 + tmp_10);
      real_t tmp_31 = tmp_9*(0.95308992296933193*tmp_0 + tmp_12);
      real_t tmp_32 = tmp_30*tmp_4 + tmp_31*tmp_5 - 1.0/3.0;
      real_t tmp_33 = tmp_3*tmp_30 + tmp_31*tmp_7 - 1.0/3.0;
      real_t a_0_0 = 0.11846344252809471*tmp_2*(tmp_16*(tmp_14*tmp_3 + tmp_15*tmp_6) + tmp_17*(tmp_14*tmp_8 + tmp_15*tmp_5)) + 0.2393143352496831*tmp_2*(tmp_16*(tmp_20*tmp_3 + tmp_21*tmp_6) + tmp_17*(tmp_20*tmp_8 + tmp_21*tmp_5)) + 0.2844444444444445*tmp_2*(tmp_16*(tmp_24*tmp_3 + tmp_25*tmp_6) + tmp_17*(tmp_24*tmp_8 + tmp_25*tmp_5)) + 0.2393143352496831*tmp_2*(tmp_16*(tmp_28*tmp_3 + tmp_29*tmp_6) + tmp_17*(tmp_28*tmp_8 + tmp_29*tmp_5)) + 0.11846344252809471*tmp_2*(tmp_16*(tmp_3*tmp_32 + tmp_33*tmp_6) + tmp_17*(tmp_32*tmp_8 + tmp_33*tmp_5));
      elMat( 0, 0) = a_0_0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >&      coordsElementInner,
                                          const std::vector< Point3D >&      coordsElementOuter,
                                          const std::vector< Point3D >&      coordsFacet,
                                          const Point3D&                     oppositeVertexInnerElement,
                                          const Point3D&                     oppositeVertexOuterElement,
                                          const Point3D&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = -p_affine_3_0 + p_affine_4_0;
      real_t tmp_4 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_5 = -p_affine_3_1 + p_affine_5_1;
      real_t tmp_6 = -tmp_4;
      real_t tmp_7 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_8 = -tmp_7;
      real_t tmp_9 = 1.0 / (tmp_3*tmp_5 - tmp_6*tmp_8);
      real_t tmp_10 = -p_affine_3_1 + p_affine_6_1;
      real_t tmp_11 = tmp_9*(0.046910077030668018*tmp_1 + tmp_10);
      real_t tmp_12 = -p_affine_3_0 + p_affine_6_0;
      real_t tmp_13 = tmp_9*(0.046910077030668018*tmp_0 + tmp_12);
      real_t tmp_14 = tmp_11*tmp_4 + tmp_13*tmp_5 - 1.0/3.0;
      real_t tmp_15 = tmp_11*tmp_3 + tmp_13*tmp_7 - 1.0/3.0;
      real_t tmp_16 = 0.5*p_affine_10_0;
      real_t tmp_17 = 0.5*p_affine_10_1;
      real_t tmp_18 = tmp_9*(0.23076534494715845*tmp_1 + tmp_10);
      real_t tmp_19 = tmp_9*(0.23076534494715845*tmp_0 + tmp_12);
      real_t tmp_20 = tmp_18*tmp_4 + tmp_19*tmp_5 - 1.0/3.0;
      real_t tmp_21 = tmp_18*tmp_3 + tmp_19*tmp_7 - 1.0/3.0;
      real_t tmp_22 = tmp_9*(0.5*tmp_1 + tmp_10);
      real_t tmp_23 = tmp_9*(0.5*tmp_0 + tmp_12);
      real_t tmp_24 = tmp_22*tmp_4 + tmp_23*tmp_5 - 1.0/3.0;
      real_t tmp_25 = tmp_22*tmp_3 + tmp_23*tmp_7 - 1.0/3.0;
      real_t tmp_26 = tmp_9*(0.7692346550528415*tmp_1 + tmp_10);
      real_t tmp_27 = tmp_9*(0.7692346550528415*tmp_0 + tmp_12);
      real_t tmp_28 = tmp_26*tmp_4 + tmp_27*tmp_5 - 1.0/3.0;
      real_t tmp_29 = tmp_26*tmp_3 + tmp_27*tmp_7 - 1.0/3.0;
      real_t tmp_30 = tmp_9*(0.95308992296933193*tmp_1 + tmp_10);
      real_t tmp_31 = tmp_9*(0.95308992296933193*tmp_0 + tmp_12);
      real_t tmp_32 = tmp_30*tmp_4 + tmp_31*tmp_5 - 1.0/3.0;
      real_t tmp_33 = tmp_3*tmp_30 + tmp_31*tmp_7 - 1.0/3.0;
      real_t a_0_0 = 0.11846344252809471*tmp_2*(-tmp_16*(tmp_14*tmp_3 + tmp_15*tmp_6) - tmp_17*(tmp_14*tmp_8 + tmp_15*tmp_5)) + 0.2393143352496831*tmp_2*(-tmp_16*(tmp_20*tmp_3 + tmp_21*tmp_6) - tmp_17*(tmp_20*tmp_8 + tmp_21*tmp_5)) + 0.2844444444444445*tmp_2*(-tmp_16*(tmp_24*tmp_3 + tmp_25*tmp_6) - tmp_17*(tmp_24*tmp_8 + tmp_25*tmp_5)) + 0.2393143352496831*tmp_2*(-tmp_16*(tmp_28*tmp_3 + tmp_29*tmp_6) - tmp_17*(tmp_28*tmp_8 + tmp_29*tmp_5)) + 0.11846344252809471*tmp_2*(-tmp_16*(tmp_3*tmp_32 + tmp_33*tmp_6) - tmp_17*(tmp_32*tmp_8 + tmp_33*tmp_5));
      elMat( 0, 0) = a_0_0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                   const std::vector< Point3D >&      coordsFacet,
                                                   const Point3D&                     oppositeVertex,
                                                   const Point3D&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   MatrixXr&                                           elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_4 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_5 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_6 = -tmp_4;
      real_t tmp_7 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_8 = -tmp_7;
      real_t tmp_9 = 1.0 / (tmp_3*tmp_5 - tmp_6*tmp_8);
      real_t tmp_10 = -p_affine_0_1 + p_affine_6_1;
      real_t tmp_11 = tmp_9*(0.046910077030668018*tmp_1 + tmp_10);
      real_t tmp_12 = -p_affine_0_0 + p_affine_6_0;
      real_t tmp_13 = tmp_9*(0.046910077030668018*tmp_0 + tmp_12);
      real_t tmp_14 = tmp_11*tmp_4 + tmp_13*tmp_5 - 1.0/3.0;
      real_t tmp_15 = tmp_11*tmp_3 + tmp_13*tmp_7 - 1.0/3.0;
      real_t tmp_16 = tmp_9*(0.23076534494715845*tmp_1 + tmp_10);
      real_t tmp_17 = tmp_9*(0.23076534494715845*tmp_0 + tmp_12);
      real_t tmp_18 = tmp_16*tmp_4 + tmp_17*tmp_5 - 1.0/3.0;
      real_t tmp_19 = tmp_16*tmp_3 + tmp_17*tmp_7 - 1.0/3.0;
      real_t tmp_20 = tmp_9*(0.5*tmp_1 + tmp_10);
      real_t tmp_21 = tmp_9*(0.5*tmp_0 + tmp_12);
      real_t tmp_22 = tmp_20*tmp_4 + tmp_21*tmp_5 - 1.0/3.0;
      real_t tmp_23 = tmp_20*tmp_3 + tmp_21*tmp_7 - 1.0/3.0;
      real_t tmp_24 = tmp_9*(0.7692346550528415*tmp_1 + tmp_10);
      real_t tmp_25 = tmp_9*(0.7692346550528415*tmp_0 + tmp_12);
      real_t tmp_26 = tmp_24*tmp_4 + tmp_25*tmp_5 - 1.0/3.0;
      real_t tmp_27 = tmp_24*tmp_3 + tmp_25*tmp_7 - 1.0/3.0;
      real_t tmp_28 = tmp_9*(0.95308992296933193*tmp_1 + tmp_10);
      real_t tmp_29 = tmp_9*(0.95308992296933193*tmp_0 + tmp_12);
      real_t tmp_30 = tmp_28*tmp_4 + tmp_29*tmp_5 - 1.0/3.0;
      real_t tmp_31 = tmp_28*tmp_3 + tmp_29*tmp_7 - 1.0/3.0;
      real_t a_0_0 = 0.11846344252809471*tmp_2*(p_affine_10_0*(tmp_14*tmp_3 + tmp_15*tmp_6) + p_affine_10_1*(tmp_14*tmp_8 + tmp_15*tmp_5)) + 0.2393143352496831*tmp_2*(p_affine_10_0*(tmp_18*tmp_3 + tmp_19*tmp_6) + p_affine_10_1*(tmp_18*tmp_8 + tmp_19*tmp_5)) + 0.2844444444444445*tmp_2*(p_affine_10_0*(tmp_22*tmp_3 + tmp_23*tmp_6) + p_affine_10_1*(tmp_22*tmp_8 + tmp_23*tmp_5)) + 0.2393143352496831*tmp_2*(p_affine_10_0*(tmp_26*tmp_3 + tmp_27*tmp_6) + p_affine_10_1*(tmp_26*tmp_8 + tmp_27*tmp_5)) + 0.11846344252809471*tmp_2*(p_affine_10_0*(tmp_3*tmp_30 + tmp_31*tmp_6) + p_affine_10_1*(tmp_30*tmp_8 + tmp_31*tmp_5));
      elMat( 0, 0) = a_0_0;
   }

    void integrateRHSDirichletBoundary2D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

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

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_g0_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id9 = 0;
      Scalar_Variable_Coefficient_2D_g0( 0.95308992296933193*p_affine_6_0 + 0.046910077030668018*p_affine_7_0, 0.95308992296933193*p_affine_6_1 + 0.046910077030668018*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id0 );
      Scalar_Variable_Coefficient_2D_g1( 0.95308992296933193*p_affine_6_0 + 0.046910077030668018*p_affine_7_0, 0.95308992296933193*p_affine_6_1 + 0.046910077030668018*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id1 );
      Scalar_Variable_Coefficient_2D_g0( 0.7692346550528415*p_affine_6_0 + 0.23076534494715845*p_affine_7_0, 0.7692346550528415*p_affine_6_1 + 0.23076534494715845*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id2 );
      Scalar_Variable_Coefficient_2D_g1( 0.7692346550528415*p_affine_6_0 + 0.23076534494715845*p_affine_7_0, 0.7692346550528415*p_affine_6_1 + 0.23076534494715845*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id3 );
      Scalar_Variable_Coefficient_2D_g0( 0.5*p_affine_6_0 + 0.5*p_affine_7_0, 0.5*p_affine_6_1 + 0.5*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id4 );
      Scalar_Variable_Coefficient_2D_g1( 0.5*p_affine_6_0 + 0.5*p_affine_7_0, 0.5*p_affine_6_1 + 0.5*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id5 );
      Scalar_Variable_Coefficient_2D_g0( 0.2307653449471585*p_affine_6_0 + 0.7692346550528415*p_affine_7_0, 0.2307653449471585*p_affine_6_1 + 0.7692346550528415*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id6 );
      Scalar_Variable_Coefficient_2D_g1( 0.2307653449471585*p_affine_6_0 + 0.7692346550528415*p_affine_7_0, 0.2307653449471585*p_affine_6_1 + 0.7692346550528415*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id7 );
      Scalar_Variable_Coefficient_2D_g0( 0.046910077030668074*p_affine_6_0 + 0.95308992296933193*p_affine_7_0, 0.046910077030668074*p_affine_6_1 + 0.95308992296933193*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id8 );
      Scalar_Variable_Coefficient_2D_g1( 0.046910077030668074*p_affine_6_0 + 0.95308992296933193*p_affine_7_0, 0.046910077030668074*p_affine_6_1 + 0.95308992296933193*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id9 );
      real_t tmp_0 = std::abs(std::pow(((-p_affine_6_0 + p_affine_7_0)*(-p_affine_6_0 + p_affine_7_0)) + ((-p_affine_6_1 + p_affine_7_1)*(-p_affine_6_1 + p_affine_7_1)), 1.0/2.0));
      real_t a_0_0 = 0.11846344252809471*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id0*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id1*p_affine_10_1) + 0.2393143352496831*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id2*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id3*p_affine_10_1) + 0.2844444444444445*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id4*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id5*p_affine_10_1) + 0.2393143352496831*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id6*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id7*p_affine_10_1) + 0.11846344252809471*tmp_0*(Scalar_Variable_Coefficient_2D_g0_out0_id8*p_affine_10_0 + Scalar_Variable_Coefficient_2D_g1_out0_id9*p_affine_10_1);
      elMat( 0, 0) = a_0_0;
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >&      coordsElement,
                                                 const std::vector< Point3D >&      coordsFacet,
                                                 const Point3D&                     oppositeVertex,
                                                 const Point3D&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                         MatrixXr&                                           elMat ) const override
   {
     elMat.resize( Eigen::Index( basis.numDoFsPerElement( 3, walberla::uint_c( degree ) ) ), 1 );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

      real_t Scalar_Variable_Coefficient_3D_g0_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id13 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id14 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id15 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id16 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id17 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id18 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id19 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id20 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id21 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id22 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id23 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id24 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id25 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id26 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id27 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id28 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id29 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id30 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id31 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id32 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id33 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id34 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id35 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id36 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id37 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id38 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id39 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id40 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id41 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id42 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id43 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id44 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id45 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id46 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id47 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id48 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id49 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id50 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id51 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id52 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id53 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id54 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id55 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id56 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id57 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id58 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id59 = 0;
      real_t Scalar_Variable_Coefficient_3D_g0_out0_id60 = 0;
      real_t Scalar_Variable_Coefficient_3D_g1_out0_id61 = 0;
      real_t Scalar_Variable_Coefficient_3D_g2_out0_id62 = 0;
      Scalar_Variable_Coefficient_3D_g0( 0.93718850182767688*p_affine_10_0 + 0.031405749086161561*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.93718850182767688*p_affine_10_1 + 0.031405749086161561*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.93718850182767688*p_affine_10_2 + 0.031405749086161561*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id0 );
      Scalar_Variable_Coefficient_3D_g1( 0.93718850182767688*p_affine_10_0 + 0.031405749086161561*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.93718850182767688*p_affine_10_1 + 0.031405749086161561*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.93718850182767688*p_affine_10_2 + 0.031405749086161561*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id1 );
      Scalar_Variable_Coefficient_3D_g2( 0.93718850182767688*p_affine_10_0 + 0.031405749086161561*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.93718850182767688*p_affine_10_1 + 0.031405749086161561*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.93718850182767688*p_affine_10_2 + 0.031405749086161561*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id2 );
      Scalar_Variable_Coefficient_3D_g0( 0.60796128279561268*p_affine_10_0 + 0.19601935860219366*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.60796128279561268*p_affine_10_1 + 0.19601935860219366*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.60796128279561268*p_affine_10_2 + 0.19601935860219366*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id3 );
      Scalar_Variable_Coefficient_3D_g1( 0.60796128279561268*p_affine_10_0 + 0.19601935860219366*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.60796128279561268*p_affine_10_1 + 0.19601935860219366*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.60796128279561268*p_affine_10_2 + 0.19601935860219366*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id4 );
      Scalar_Variable_Coefficient_3D_g2( 0.60796128279561268*p_affine_10_0 + 0.19601935860219366*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.60796128279561268*p_affine_10_1 + 0.19601935860219366*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.60796128279561268*p_affine_10_2 + 0.19601935860219366*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id5 );
      Scalar_Variable_Coefficient_3D_g0( 0.039308471900058539*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id6 );
      Scalar_Variable_Coefficient_3D_g1( 0.039308471900058539*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id7 );
      Scalar_Variable_Coefficient_3D_g2( 0.039308471900058539*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id8 );
      Scalar_Variable_Coefficient_3D_g0( 0.1711304259088916*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id9 );
      Scalar_Variable_Coefficient_3D_g1( 0.1711304259088916*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id10 );
      Scalar_Variable_Coefficient_3D_g2( 0.1711304259088916*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id11 );
      Scalar_Variable_Coefficient_3D_g0( 0.37605877282253791*p_affine_10_0 + 0.039308471900058539*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.039308471900058539*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.039308471900058539*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id12 );
      Scalar_Variable_Coefficient_3D_g1( 0.37605877282253791*p_affine_10_0 + 0.039308471900058539*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.039308471900058539*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.039308471900058539*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id13 );
      Scalar_Variable_Coefficient_3D_g2( 0.37605877282253791*p_affine_10_0 + 0.039308471900058539*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.039308471900058539*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.039308471900058539*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id14 );
      Scalar_Variable_Coefficient_3D_g0( 0.78764240869137092*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id15 );
      Scalar_Variable_Coefficient_3D_g1( 0.78764240869137092*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id16 );
      Scalar_Variable_Coefficient_3D_g2( 0.78764240869137092*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id17 );
      Scalar_Variable_Coefficient_3D_g0( 0.58463275527740355*p_affine_10_0 + 0.37605877282253797*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.37605877282253797*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.37605877282253797*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id18 );
      Scalar_Variable_Coefficient_3D_g1( 0.58463275527740355*p_affine_10_0 + 0.37605877282253797*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.37605877282253797*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.37605877282253797*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id19 );
      Scalar_Variable_Coefficient_3D_g2( 0.58463275527740355*p_affine_10_0 + 0.37605877282253797*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.37605877282253797*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.37605877282253797*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id20 );
      Scalar_Variable_Coefficient_3D_g0( 0.041227165399737475*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id21 );
      Scalar_Variable_Coefficient_3D_g1( 0.041227165399737475*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id22 );
      Scalar_Variable_Coefficient_3D_g2( 0.041227165399737475*p_affine_10_0 + 0.1711304259088916*p_affine_8_0 + 0.78764240869137092*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.1711304259088916*p_affine_8_1 + 0.78764240869137092*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.1711304259088916*p_affine_8_2 + 0.78764240869137092*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id23 );
      Scalar_Variable_Coefficient_3D_g0( 0.039308471900058539*p_affine_10_0 + 0.37605877282253791*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.37605877282253791*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.37605877282253791*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id24 );
      Scalar_Variable_Coefficient_3D_g1( 0.039308471900058539*p_affine_10_0 + 0.37605877282253791*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.37605877282253791*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.37605877282253791*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id25 );
      Scalar_Variable_Coefficient_3D_g2( 0.039308471900058539*p_affine_10_0 + 0.37605877282253791*p_affine_8_0 + 0.58463275527740355*p_affine_9_0, 0.039308471900058539*p_affine_10_1 + 0.37605877282253791*p_affine_8_1 + 0.58463275527740355*p_affine_9_1, 0.039308471900058539*p_affine_10_2 + 0.37605877282253791*p_affine_8_2 + 0.58463275527740355*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id26 );
      Scalar_Variable_Coefficient_3D_g0( 0.78764240869137092*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id27 );
      Scalar_Variable_Coefficient_3D_g1( 0.78764240869137092*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id28 );
      Scalar_Variable_Coefficient_3D_g2( 0.78764240869137092*p_affine_10_0 + 0.041227165399737475*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.78764240869137092*p_affine_10_1 + 0.041227165399737475*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.78764240869137092*p_affine_10_2 + 0.041227165399737475*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id29 );
      Scalar_Variable_Coefficient_3D_g0( 0.58463275527740355*p_affine_10_0 + 0.039308471900058484*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.039308471900058484*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.039308471900058484*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id30 );
      Scalar_Variable_Coefficient_3D_g1( 0.58463275527740355*p_affine_10_0 + 0.039308471900058484*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.039308471900058484*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.039308471900058484*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id31 );
      Scalar_Variable_Coefficient_3D_g2( 0.58463275527740355*p_affine_10_0 + 0.039308471900058484*p_affine_8_0 + 0.37605877282253791*p_affine_9_0, 0.58463275527740355*p_affine_10_1 + 0.039308471900058484*p_affine_8_1 + 0.37605877282253791*p_affine_9_1, 0.58463275527740355*p_affine_10_2 + 0.039308471900058484*p_affine_8_2 + 0.37605877282253791*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id32 );
      Scalar_Variable_Coefficient_3D_g0( 0.1711304259088916*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id33 );
      Scalar_Variable_Coefficient_3D_g1( 0.1711304259088916*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id34 );
      Scalar_Variable_Coefficient_3D_g2( 0.1711304259088916*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.041227165399737475*p_affine_9_0, 0.1711304259088916*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.041227165399737475*p_affine_9_1, 0.1711304259088916*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.041227165399737475*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id35 );
      Scalar_Variable_Coefficient_3D_g0( 0.19107600050469298*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.19107600050469298*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.19107600050469298*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id36 );
      Scalar_Variable_Coefficient_3D_g1( 0.19107600050469298*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.19107600050469298*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.19107600050469298*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id37 );
      Scalar_Variable_Coefficient_3D_g2( 0.19107600050469298*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.19107600050469298*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.19107600050469298*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id38 );
      Scalar_Variable_Coefficient_3D_g0( 0.37605877282253791*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id39 );
      Scalar_Variable_Coefficient_3D_g1( 0.37605877282253791*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id40 );
      Scalar_Variable_Coefficient_3D_g2( 0.37605877282253791*p_affine_10_0 + 0.58463275527740355*p_affine_8_0 + 0.039308471900058539*p_affine_9_0, 0.37605877282253791*p_affine_10_1 + 0.58463275527740355*p_affine_8_1 + 0.039308471900058539*p_affine_9_1, 0.37605877282253791*p_affine_10_2 + 0.58463275527740355*p_affine_8_2 + 0.039308471900058539*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id41 );
      Scalar_Variable_Coefficient_3D_g0( 0.031405749086161582*p_affine_10_0 + 0.03140574908616154*p_affine_8_0 + 0.93718850182767688*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.03140574908616154*p_affine_8_1 + 0.93718850182767688*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.03140574908616154*p_affine_8_2 + 0.93718850182767688*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id42 );
      Scalar_Variable_Coefficient_3D_g1( 0.031405749086161582*p_affine_10_0 + 0.03140574908616154*p_affine_8_0 + 0.93718850182767688*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.03140574908616154*p_affine_8_1 + 0.93718850182767688*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.03140574908616154*p_affine_8_2 + 0.93718850182767688*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id43 );
      Scalar_Variable_Coefficient_3D_g2( 0.031405749086161582*p_affine_10_0 + 0.03140574908616154*p_affine_8_0 + 0.93718850182767688*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.03140574908616154*p_affine_8_1 + 0.93718850182767688*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.03140574908616154*p_affine_8_2 + 0.93718850182767688*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id44 );
      Scalar_Variable_Coefficient_3D_g0( 0.19601935860219369*p_affine_10_0 + 0.19601935860219363*p_affine_8_0 + 0.60796128279561268*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.19601935860219363*p_affine_8_1 + 0.60796128279561268*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.19601935860219363*p_affine_8_2 + 0.60796128279561268*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id45 );
      Scalar_Variable_Coefficient_3D_g1( 0.19601935860219369*p_affine_10_0 + 0.19601935860219363*p_affine_8_0 + 0.60796128279561268*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.19601935860219363*p_affine_8_1 + 0.60796128279561268*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.19601935860219363*p_affine_8_2 + 0.60796128279561268*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id46 );
      Scalar_Variable_Coefficient_3D_g2( 0.19601935860219369*p_affine_10_0 + 0.19601935860219363*p_affine_8_0 + 0.60796128279561268*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.19601935860219363*p_affine_8_1 + 0.60796128279561268*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.19601935860219363*p_affine_8_2 + 0.60796128279561268*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id47 );
      Scalar_Variable_Coefficient_3D_g0( 0.40446199974765351*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.19107600050469298*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.19107600050469298*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.19107600050469298*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id48 );
      Scalar_Variable_Coefficient_3D_g1( 0.40446199974765351*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.19107600050469298*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.19107600050469298*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.19107600050469298*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id49 );
      Scalar_Variable_Coefficient_3D_g2( 0.40446199974765351*p_affine_10_0 + 0.40446199974765351*p_affine_8_0 + 0.19107600050469298*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.40446199974765351*p_affine_8_1 + 0.19107600050469298*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.40446199974765351*p_affine_8_2 + 0.19107600050469298*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id50 );
      Scalar_Variable_Coefficient_3D_g0( 0.031405749086161582*p_affine_10_0 + 0.93718850182767688*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.93718850182767688*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.93718850182767688*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id51 );
      Scalar_Variable_Coefficient_3D_g1( 0.031405749086161582*p_affine_10_0 + 0.93718850182767688*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.93718850182767688*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.93718850182767688*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id52 );
      Scalar_Variable_Coefficient_3D_g2( 0.031405749086161582*p_affine_10_0 + 0.93718850182767688*p_affine_8_0 + 0.031405749086161582*p_affine_9_0, 0.031405749086161582*p_affine_10_1 + 0.93718850182767688*p_affine_8_1 + 0.031405749086161582*p_affine_9_1, 0.031405749086161582*p_affine_10_2 + 0.93718850182767688*p_affine_8_2 + 0.031405749086161582*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id53 );
      Scalar_Variable_Coefficient_3D_g0( 0.19601935860219369*p_affine_10_0 + 0.60796128279561268*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.60796128279561268*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.60796128279561268*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id54 );
      Scalar_Variable_Coefficient_3D_g1( 0.19601935860219369*p_affine_10_0 + 0.60796128279561268*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.60796128279561268*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.60796128279561268*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id55 );
      Scalar_Variable_Coefficient_3D_g2( 0.19601935860219369*p_affine_10_0 + 0.60796128279561268*p_affine_8_0 + 0.19601935860219369*p_affine_9_0, 0.19601935860219369*p_affine_10_1 + 0.60796128279561268*p_affine_8_1 + 0.19601935860219369*p_affine_9_1, 0.19601935860219369*p_affine_10_2 + 0.60796128279561268*p_affine_8_2 + 0.19601935860219369*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id56 );
      Scalar_Variable_Coefficient_3D_g0( 0.40446199974765351*p_affine_10_0 + 0.19107600050469298*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.19107600050469298*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.19107600050469298*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id57 );
      Scalar_Variable_Coefficient_3D_g1( 0.40446199974765351*p_affine_10_0 + 0.19107600050469298*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.19107600050469298*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.19107600050469298*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id58 );
      Scalar_Variable_Coefficient_3D_g2( 0.40446199974765351*p_affine_10_0 + 0.19107600050469298*p_affine_8_0 + 0.40446199974765351*p_affine_9_0, 0.40446199974765351*p_affine_10_1 + 0.19107600050469298*p_affine_8_1 + 0.40446199974765351*p_affine_9_1, 0.40446199974765351*p_affine_10_2 + 0.19107600050469298*p_affine_8_2 + 0.40446199974765351*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id59 );
      Scalar_Variable_Coefficient_3D_g0( 0.041227165399737475*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g0_out0_id60 );
      Scalar_Variable_Coefficient_3D_g1( 0.041227165399737475*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g1_out0_id61 );
      Scalar_Variable_Coefficient_3D_g2( 0.041227165399737475*p_affine_10_0 + 0.78764240869137092*p_affine_8_0 + 0.1711304259088916*p_affine_9_0, 0.041227165399737475*p_affine_10_1 + 0.78764240869137092*p_affine_8_1 + 0.1711304259088916*p_affine_9_1, 0.041227165399737475*p_affine_10_2 + 0.78764240869137092*p_affine_8_2 + 0.1711304259088916*p_affine_9_2, &Scalar_Variable_Coefficient_3D_g2_out0_id62 );
      real_t tmp_0 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_1 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_2 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_3 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_4 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_5 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_6 = 1.0*std::pow((std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)*std::abs(tmp_0*tmp_1 - tmp_2*tmp_3)) + (std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)*std::abs(tmp_0*tmp_4 - tmp_3*tmp_5)) + (std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)*std::abs(tmp_1*tmp_5 - tmp_2*tmp_4)), 1.0/2.0);
      real_t a_0_0 = 0.0068572537431980923*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id0*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id1*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id2*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id12*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id13*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id14*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id15*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id16*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id17*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id18*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id19*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id20*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id21*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id22*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id23*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id24*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id25*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id26*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id27*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id28*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id29*p_affine_13_2) + 0.037198804536718075*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id3*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id4*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id5*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id30*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id31*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id32*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id33*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id34*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id35*p_affine_13_2) + 0.042507265838595799*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id36*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id37*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id38*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id39*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id40*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id41*p_affine_13_2) + 0.0068572537431980923*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id42*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id43*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id44*p_affine_13_2) + 0.037198804536718075*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id45*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id46*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id47*p_affine_13_2) + 0.042507265838595799*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id48*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id49*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id50*p_affine_13_2) + 0.0068572537431980923*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id51*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id52*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id53*p_affine_13_2) + 0.037198804536718075*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id54*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id55*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id56*p_affine_13_2) + 0.042507265838595799*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id57*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id58*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id59*p_affine_13_2) + 0.020848748529055869*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id6*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id7*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id8*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id60*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id61*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id62*p_affine_13_2) + 0.019202922745021479*tmp_6*(Scalar_Variable_Coefficient_3D_g0_out0_id9*p_affine_13_0 + Scalar_Variable_Coefficient_3D_g1_out0_id10*p_affine_13_1 + Scalar_Variable_Coefficient_3D_g2_out0_id11*p_affine_13_2);
      elMat( 0, 0) = a_0_0;
   }
   void integrateVolume3D( const std::vector< Point3D >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                           MatrixXr&                                           elMat ) const override
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coords[0]( 0 );
   const auto p_affine_0_1 = coords[0]( 1 );
   const auto p_affine_0_2 = coords[0]( 2 );

   const auto p_affine_1_0 = coords[1]( 0 );
   const auto p_affine_1_1 = coords[1]( 1 );
   const auto p_affine_1_2 = coords[1]( 2 );

   const auto p_affine_2_0 = coords[2]( 0 );
   const auto p_affine_2_1 = coords[2]( 1 );
   const auto p_affine_2_2 = coords[2]( 2 );

   const auto p_affine_3_0 = coords[3]( 0 );
   const auto p_affine_3_1 = coords[3]( 1 );
   const auto p_affine_3_2 = coords[3]( 2 );

      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_2 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_3 = tmp_1*tmp_2;
      real_t tmp_4 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_5 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_8 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_9 = tmp_4*tmp_8;
      real_t tmp_10 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_11 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_12 = tmp_2*tmp_8;
      real_t tmp_13 = tmp_1*tmp_11;
      real_t tmp_14 = 1.0 / (tmp_0*tmp_3 - tmp_0*tmp_6 + tmp_10*tmp_11*tmp_5 - tmp_10*tmp_12 - tmp_13*tmp_7 + tmp_7*tmp_9);
      real_t tmp_15 = p_affine_0_0*p_affine_1_1;
      real_t tmp_16 = p_affine_0_0*p_affine_1_2;
      real_t tmp_17 = p_affine_2_1*p_affine_3_2;
      real_t tmp_18 = p_affine_0_1*p_affine_1_0;
      real_t tmp_19 = p_affine_0_1*p_affine_1_2;
      real_t tmp_20 = p_affine_2_2*p_affine_3_0;
      real_t tmp_21 = p_affine_0_2*p_affine_1_0;
      real_t tmp_22 = p_affine_0_2*p_affine_1_1;
      real_t tmp_23 = p_affine_2_0*p_affine_3_1;
      real_t tmp_24 = p_affine_2_2*p_affine_3_1;
      real_t tmp_25 = p_affine_2_0*p_affine_3_2;
      real_t tmp_26 = p_affine_2_1*p_affine_3_0;
      real_t tmp_27 = (-tmp_0*tmp_14*(tmp_3 - tmp_6) - tmp_1*tmp_14*(tmp_0*tmp_2 - tmp_11*tmp_7) - tmp_10*tmp_14*(tmp_11*tmp_5 - tmp_12) - tmp_11*tmp_14*(-tmp_1*tmp_7 + tmp_10*tmp_5) - tmp_14*tmp_2*(tmp_0*tmp_1 - tmp_10*tmp_8) - tmp_14*tmp_4*(-tmp_0*tmp_5 + tmp_7*tmp_8) - tmp_14*tmp_5*(-tmp_0*tmp_4 + tmp_10*tmp_11) - tmp_14*tmp_7*(-tmp_13 + tmp_9) - tmp_14*tmp_8*(-tmp_10*tmp_2 + tmp_4*tmp_7))*std::abs(p_affine_0_0*tmp_17 - p_affine_0_0*tmp_24 + p_affine_0_1*tmp_20 - p_affine_0_1*tmp_25 + p_affine_0_2*tmp_23 - p_affine_0_2*tmp_26 - p_affine_1_0*tmp_17 + p_affine_1_0*tmp_24 - p_affine_1_1*tmp_20 + p_affine_1_1*tmp_25 - p_affine_1_2*tmp_23 + p_affine_1_2*tmp_26 + p_affine_2_0*tmp_19 - p_affine_2_0*tmp_22 - p_affine_2_1*tmp_16 + p_affine_2_1*tmp_21 + p_affine_2_2*tmp_15 - p_affine_2_2*tmp_18 - p_affine_3_0*tmp_19 + p_affine_3_0*tmp_22 + p_affine_3_1*tmp_16 - p_affine_3_1*tmp_21 - p_affine_3_2*tmp_15 + p_affine_3_2*tmp_18);
      real_t a_0_0 = 0.1666666666666668*tmp_27;
      elMat( 0, 0) = a_0_0;
   }



   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                                                     const std::vector< Point3D >& coordsFacet,
                                                     const Point3D&,
                                                     const Point3D&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                               MatrixXr&                            elMat ) const override
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );

         real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = tmp_1*tmp_2;
      real_t tmp_4 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_5 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = tmp_3 - tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_9 = tmp_5*tmp_8;
      real_t tmp_10 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_11 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_12 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_13 = tmp_12*tmp_2;
      real_t tmp_14 = tmp_1*tmp_8;
      real_t tmp_15 = 1.0 / (-tmp_0*tmp_13 + tmp_0*tmp_9 + tmp_10*tmp_3 - tmp_10*tmp_6 + tmp_11*tmp_12*tmp_4 - tmp_11*tmp_14);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.031405749086161582*tmp_17 + 0.93718850182767688*tmp_18 + tmp_19);
      real_t tmp_21 = tmp_12*tmp_4 - tmp_14;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.031405749086161582*tmp_23 + 0.93718850182767688*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_13 + tmp_9;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_7 + tmp_21*tmp_26 + tmp_27*tmp_32 - 1.0/4.0;
      real_t tmp_34 = -tmp_0*tmp_2 + tmp_11*tmp_4;
      real_t tmp_35 = tmp_0*tmp_8 - tmp_10*tmp_4;
      real_t tmp_36 = tmp_10*tmp_2 - tmp_11*tmp_8;
      real_t tmp_37 = tmp_20*tmp_34 + tmp_26*tmp_35 + tmp_32*tmp_36 - 1.0/4.0;
      real_t tmp_38 = tmp_0*tmp_5 - tmp_1*tmp_11;
      real_t tmp_39 = -tmp_0*tmp_12 + tmp_1*tmp_10;
      real_t tmp_40 = -tmp_10*tmp_5 + tmp_11*tmp_12;
      real_t tmp_41 = tmp_20*tmp_38 + tmp_26*tmp_39 + tmp_32*tmp_40 - 1.0/4.0;
      real_t tmp_42 = 0.5*p_affine_13_0;
      real_t tmp_43 = 0.5*p_affine_13_1;
      real_t tmp_44 = 0.5*p_affine_13_2;
      real_t tmp_45 = 1.0*std::pow((std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)), 1.0/2.0);
      real_t tmp_46 = tmp_15*(0.19601935860219369*tmp_17 + 0.60796128279561268*tmp_18 + tmp_19);
      real_t tmp_47 = tmp_15*(0.19601935860219369*tmp_23 + 0.60796128279561268*tmp_24 + tmp_25);
      real_t tmp_48 = tmp_15*(0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_30 + tmp_31);
      real_t tmp_49 = tmp_21*tmp_47 + tmp_27*tmp_48 + tmp_46*tmp_7 - 1.0/4.0;
      real_t tmp_50 = tmp_34*tmp_46 + tmp_35*tmp_47 + tmp_36*tmp_48 - 1.0/4.0;
      real_t tmp_51 = tmp_38*tmp_46 + tmp_39*tmp_47 + tmp_40*tmp_48 - 1.0/4.0;
      real_t tmp_52 = tmp_15*(0.37605877282253791*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_53 = tmp_15*(0.37605877282253791*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_54 = tmp_15*(0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_55 = tmp_21*tmp_53 + tmp_27*tmp_54 + tmp_52*tmp_7 - 1.0/4.0;
      real_t tmp_56 = tmp_34*tmp_52 + tmp_35*tmp_53 + tmp_36*tmp_54 - 1.0/4.0;
      real_t tmp_57 = tmp_38*tmp_52 + tmp_39*tmp_53 + tmp_40*tmp_54 - 1.0/4.0;
      real_t tmp_58 = tmp_15*(0.78764240869137092*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_59 = tmp_15*(0.78764240869137092*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_60 = tmp_15*(0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_61 = tmp_21*tmp_59 + tmp_27*tmp_60 + tmp_58*tmp_7 - 1.0/4.0;
      real_t tmp_62 = tmp_34*tmp_58 + tmp_35*tmp_59 + tmp_36*tmp_60 - 1.0/4.0;
      real_t tmp_63 = tmp_38*tmp_58 + tmp_39*tmp_59 + tmp_40*tmp_60 - 1.0/4.0;
      real_t tmp_64 = tmp_15*(0.58463275527740355*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_65 = tmp_15*(0.58463275527740355*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_66 = tmp_15*(0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_67 = tmp_21*tmp_65 + tmp_27*tmp_66 + tmp_64*tmp_7 - 1.0/4.0;
      real_t tmp_68 = tmp_34*tmp_64 + tmp_35*tmp_65 + tmp_36*tmp_66 - 1.0/4.0;
      real_t tmp_69 = tmp_38*tmp_64 + tmp_39*tmp_65 + tmp_40*tmp_66 - 1.0/4.0;
      real_t tmp_70 = tmp_15*(0.041227165399737475*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_71 = tmp_15*(0.041227165399737475*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_72 = tmp_15*(0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_73 = tmp_21*tmp_71 + tmp_27*tmp_72 + tmp_7*tmp_70 - 1.0/4.0;
      real_t tmp_74 = tmp_34*tmp_70 + tmp_35*tmp_71 + tmp_36*tmp_72 - 1.0/4.0;
      real_t tmp_75 = tmp_38*tmp_70 + tmp_39*tmp_71 + tmp_40*tmp_72 - 1.0/4.0;
      real_t tmp_76 = tmp_15*(0.039308471900058539*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_77 = tmp_15*(0.039308471900058539*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_78 = tmp_15*(0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_79 = tmp_21*tmp_77 + tmp_27*tmp_78 + tmp_7*tmp_76 - 1.0/4.0;
      real_t tmp_80 = tmp_34*tmp_76 + tmp_35*tmp_77 + tmp_36*tmp_78 - 1.0/4.0;
      real_t tmp_81 = tmp_38*tmp_76 + tmp_39*tmp_77 + tmp_40*tmp_78 - 1.0/4.0;
      real_t tmp_82 = tmp_15*(0.78764240869137092*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_83 = tmp_15*(0.78764240869137092*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_84 = tmp_15*(0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_85 = tmp_21*tmp_83 + tmp_27*tmp_84 + tmp_7*tmp_82 - 1.0/4.0;
      real_t tmp_86 = tmp_34*tmp_82 + tmp_35*tmp_83 + tmp_36*tmp_84 - 1.0/4.0;
      real_t tmp_87 = tmp_38*tmp_82 + tmp_39*tmp_83 + tmp_40*tmp_84 - 1.0/4.0;
      real_t tmp_88 = tmp_15*(0.58463275527740355*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_89 = tmp_15*(0.58463275527740355*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_90 = tmp_15*(0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_91 = tmp_21*tmp_89 + tmp_27*tmp_90 + tmp_7*tmp_88 - 1.0/4.0;
      real_t tmp_92 = tmp_34*tmp_88 + tmp_35*tmp_89 + tmp_36*tmp_90 - 1.0/4.0;
      real_t tmp_93 = tmp_38*tmp_88 + tmp_39*tmp_89 + tmp_40*tmp_90 - 1.0/4.0;
      real_t tmp_94 = tmp_15*(0.1711304259088916*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_95 = tmp_15*(0.1711304259088916*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_96 = tmp_15*(0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_97 = tmp_21*tmp_95 + tmp_27*tmp_96 + tmp_7*tmp_94 - 1.0/4.0;
      real_t tmp_98 = tmp_34*tmp_94 + tmp_35*tmp_95 + tmp_36*tmp_96 - 1.0/4.0;
      real_t tmp_99 = tmp_38*tmp_94 + tmp_39*tmp_95 + tmp_40*tmp_96 - 1.0/4.0;
      real_t tmp_100 = tmp_15*(0.37605877282253791*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_101 = tmp_15*(0.37605877282253791*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_102 = tmp_15*(0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_103 = tmp_100*tmp_7 + tmp_101*tmp_21 + tmp_102*tmp_27 - 1.0/4.0;
      real_t tmp_104 = tmp_100*tmp_34 + tmp_101*tmp_35 + tmp_102*tmp_36 - 1.0/4.0;
      real_t tmp_105 = tmp_100*tmp_38 + tmp_101*tmp_39 + tmp_102*tmp_40 - 1.0/4.0;
      real_t tmp_106 = tmp_15*(0.041227165399737475*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_107 = tmp_15*(0.041227165399737475*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_108 = tmp_15*(0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_109 = tmp_106*tmp_7 + tmp_107*tmp_21 + tmp_108*tmp_27 - 1.0/4.0;
      real_t tmp_110 = tmp_106*tmp_34 + tmp_107*tmp_35 + tmp_108*tmp_36 - 1.0/4.0;
      real_t tmp_111 = tmp_106*tmp_38 + tmp_107*tmp_39 + tmp_108*tmp_40 - 1.0/4.0;
      real_t tmp_112 = tmp_15*(0.40446199974765351*tmp_17 + 0.19107600050469298*tmp_18 + tmp_19);
      real_t tmp_113 = tmp_15*(0.40446199974765351*tmp_23 + 0.19107600050469298*tmp_24 + tmp_25);
      real_t tmp_114 = tmp_15*(0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_30 + tmp_31);
      real_t tmp_115 = tmp_112*tmp_7 + tmp_113*tmp_21 + tmp_114*tmp_27 - 1.0/4.0;
      real_t tmp_116 = tmp_112*tmp_34 + tmp_113*tmp_35 + tmp_114*tmp_36 - 1.0/4.0;
      real_t tmp_117 = tmp_112*tmp_38 + tmp_113*tmp_39 + tmp_114*tmp_40 - 1.0/4.0;
      real_t tmp_118 = tmp_15*(0.039308471900058539*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_119 = tmp_15*(0.039308471900058539*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_120 = tmp_15*(0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_121 = tmp_118*tmp_7 + tmp_119*tmp_21 + tmp_120*tmp_27 - 1.0/4.0;
      real_t tmp_122 = tmp_118*tmp_34 + tmp_119*tmp_35 + tmp_120*tmp_36 - 1.0/4.0;
      real_t tmp_123 = tmp_118*tmp_38 + tmp_119*tmp_39 + tmp_120*tmp_40 - 1.0/4.0;
      real_t tmp_124 = tmp_15*(0.93718850182767688*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_125 = tmp_15*(0.93718850182767688*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_126 = tmp_15*(0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_127 = tmp_124*tmp_7 + tmp_125*tmp_21 + tmp_126*tmp_27 - 1.0/4.0;
      real_t tmp_128 = tmp_124*tmp_34 + tmp_125*tmp_35 + tmp_126*tmp_36 - 1.0/4.0;
      real_t tmp_129 = tmp_124*tmp_38 + tmp_125*tmp_39 + tmp_126*tmp_40 - 1.0/4.0;
      real_t tmp_130 = tmp_15*(0.60796128279561268*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_131 = tmp_15*(0.60796128279561268*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_132 = tmp_15*(0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_133 = tmp_130*tmp_7 + tmp_131*tmp_21 + tmp_132*tmp_27 - 1.0/4.0;
      real_t tmp_134 = tmp_130*tmp_34 + tmp_131*tmp_35 + tmp_132*tmp_36 - 1.0/4.0;
      real_t tmp_135 = tmp_130*tmp_38 + tmp_131*tmp_39 + tmp_132*tmp_40 - 1.0/4.0;
      real_t tmp_136 = tmp_15*(0.19107600050469298*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_137 = tmp_15*(0.19107600050469298*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_138 = tmp_15*(0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_139 = tmp_136*tmp_7 + tmp_137*tmp_21 + tmp_138*tmp_27 - 1.0/4.0;
      real_t tmp_140 = tmp_136*tmp_34 + tmp_137*tmp_35 + tmp_138*tmp_36 - 1.0/4.0;
      real_t tmp_141 = tmp_136*tmp_38 + tmp_137*tmp_39 + tmp_138*tmp_40 - 1.0/4.0;
      real_t tmp_142 = tmp_15*(0.031405749086161582*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_143 = tmp_15*(0.031405749086161582*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_144 = tmp_15*(0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_145 = tmp_142*tmp_7 + tmp_143*tmp_21 + tmp_144*tmp_27 - 1.0/4.0;
      real_t tmp_146 = tmp_142*tmp_34 + tmp_143*tmp_35 + tmp_144*tmp_36 - 1.0/4.0;
      real_t tmp_147 = tmp_142*tmp_38 + tmp_143*tmp_39 + tmp_144*tmp_40 - 1.0/4.0;
      real_t tmp_148 = tmp_15*(0.19601935860219369*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_149 = tmp_15*(0.19601935860219369*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_150 = tmp_15*(0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_151 = tmp_148*tmp_7 + tmp_149*tmp_21 + tmp_150*tmp_27 - 1.0/4.0;
      real_t tmp_152 = tmp_148*tmp_34 + tmp_149*tmp_35 + tmp_150*tmp_36 - 1.0/4.0;
      real_t tmp_153 = tmp_148*tmp_38 + tmp_149*tmp_39 + tmp_150*tmp_40 - 1.0/4.0;
      real_t tmp_154 = tmp_15*(0.40446199974765351*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_155 = tmp_15*(0.40446199974765351*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_156 = tmp_15*(0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_157 = tmp_154*tmp_7 + tmp_155*tmp_21 + tmp_156*tmp_27 - 1.0/4.0;
      real_t tmp_158 = tmp_154*tmp_34 + tmp_155*tmp_35 + tmp_156*tmp_36 - 1.0/4.0;
      real_t tmp_159 = tmp_154*tmp_38 + tmp_155*tmp_39 + tmp_156*tmp_40 - 1.0/4.0;
      real_t tmp_160 = tmp_15*(0.1711304259088916*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_161 = tmp_15*(0.1711304259088916*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_162 = tmp_15*(0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_163 = tmp_160*tmp_7 + tmp_161*tmp_21 + tmp_162*tmp_27 - 1.0/4.0;
      real_t tmp_164 = tmp_160*tmp_34 + tmp_161*tmp_35 + tmp_162*tmp_36 - 1.0/4.0;
      real_t tmp_165 = tmp_160*tmp_38 + tmp_161*tmp_39 + tmp_162*tmp_40 - 1.0/4.0;
      real_t a_0_0 = 0.020848748529055869*tmp_45*(tmp_42*(tmp_0*tmp_103 + tmp_1*tmp_104 + tmp_105*tmp_4) + tmp_43*(tmp_103*tmp_11 + tmp_104*tmp_5 + tmp_105*tmp_2) + tmp_44*(tmp_10*tmp_103 + tmp_104*tmp_12 + tmp_105*tmp_8)) + 0.019202922745021479*tmp_45*(tmp_42*(tmp_0*tmp_109 + tmp_1*tmp_110 + tmp_111*tmp_4) + tmp_43*(tmp_109*tmp_11 + tmp_110*tmp_5 + tmp_111*tmp_2) + tmp_44*(tmp_10*tmp_109 + tmp_110*tmp_12 + tmp_111*tmp_8)) + 0.042507265838595799*tmp_45*(tmp_42*(tmp_0*tmp_115 + tmp_1*tmp_116 + tmp_117*tmp_4) + tmp_43*(tmp_11*tmp_115 + tmp_116*tmp_5 + tmp_117*tmp_2) + tmp_44*(tmp_10*tmp_115 + tmp_116*tmp_12 + tmp_117*tmp_8)) + 0.020848748529055869*tmp_45*(tmp_42*(tmp_0*tmp_121 + tmp_1*tmp_122 + tmp_123*tmp_4) + tmp_43*(tmp_11*tmp_121 + tmp_122*tmp_5 + tmp_123*tmp_2) + tmp_44*(tmp_10*tmp_121 + tmp_12*tmp_122 + tmp_123*tmp_8)) + 0.0068572537431980923*tmp_45*(tmp_42*(tmp_0*tmp_127 + tmp_1*tmp_128 + tmp_129*tmp_4) + tmp_43*(tmp_11*tmp_127 + tmp_128*tmp_5 + tmp_129*tmp_2) + tmp_44*(tmp_10*tmp_127 + tmp_12*tmp_128 + tmp_129*tmp_8)) + 0.037198804536718075*tmp_45*(tmp_42*(tmp_0*tmp_133 + tmp_1*tmp_134 + tmp_135*tmp_4) + tmp_43*(tmp_11*tmp_133 + tmp_134*tmp_5 + tmp_135*tmp_2) + tmp_44*(tmp_10*tmp_133 + tmp_12*tmp_134 + tmp_135*tmp_8)) + 0.042507265838595799*tmp_45*(tmp_42*(tmp_0*tmp_139 + tmp_1*tmp_140 + tmp_141*tmp_4) + tmp_43*(tmp_11*tmp_139 + tmp_140*tmp_5 + tmp_141*tmp_2) + tmp_44*(tmp_10*tmp_139 + tmp_12*tmp_140 + tmp_141*tmp_8)) + 0.0068572537431980923*tmp_45*(tmp_42*(tmp_0*tmp_145 + tmp_1*tmp_146 + tmp_147*tmp_4) + tmp_43*(tmp_11*tmp_145 + tmp_146*tmp_5 + tmp_147*tmp_2) + tmp_44*(tmp_10*tmp_145 + tmp_12*tmp_146 + tmp_147*tmp_8)) + 0.037198804536718075*tmp_45*(tmp_42*(tmp_0*tmp_151 + tmp_1*tmp_152 + tmp_153*tmp_4) + tmp_43*(tmp_11*tmp_151 + tmp_152*tmp_5 + tmp_153*tmp_2) + tmp_44*(tmp_10*tmp_151 + tmp_12*tmp_152 + tmp_153*tmp_8)) + 0.042507265838595799*tmp_45*(tmp_42*(tmp_0*tmp_157 + tmp_1*tmp_158 + tmp_159*tmp_4) + tmp_43*(tmp_11*tmp_157 + tmp_158*tmp_5 + tmp_159*tmp_2) + tmp_44*(tmp_10*tmp_157 + tmp_12*tmp_158 + tmp_159*tmp_8)) + 0.019202922745021479*tmp_45*(tmp_42*(tmp_0*tmp_163 + tmp_1*tmp_164 + tmp_165*tmp_4) + tmp_43*(tmp_11*tmp_163 + tmp_164*tmp_5 + tmp_165*tmp_2) + tmp_44*(tmp_10*tmp_163 + tmp_12*tmp_164 + tmp_165*tmp_8)) + 0.0068572537431980923*tmp_45*(tmp_42*(tmp_0*tmp_33 + tmp_1*tmp_37 + tmp_4*tmp_41) + tmp_43*(tmp_11*tmp_33 + tmp_2*tmp_41 + tmp_37*tmp_5) + tmp_44*(tmp_10*tmp_33 + tmp_12*tmp_37 + tmp_41*tmp_8)) + 0.037198804536718075*tmp_45*(tmp_42*(tmp_0*tmp_49 + tmp_1*tmp_50 + tmp_4*tmp_51) + tmp_43*(tmp_11*tmp_49 + tmp_2*tmp_51 + tmp_5*tmp_50) + tmp_44*(tmp_10*tmp_49 + tmp_12*tmp_50 + tmp_51*tmp_8)) + 0.020848748529055869*tmp_45*(tmp_42*(tmp_0*tmp_55 + tmp_1*tmp_56 + tmp_4*tmp_57) + tmp_43*(tmp_11*tmp_55 + tmp_2*tmp_57 + tmp_5*tmp_56) + tmp_44*(tmp_10*tmp_55 + tmp_12*tmp_56 + tmp_57*tmp_8)) + 0.019202922745021479*tmp_45*(tmp_42*(tmp_0*tmp_61 + tmp_1*tmp_62 + tmp_4*tmp_63) + tmp_43*(tmp_11*tmp_61 + tmp_2*tmp_63 + tmp_5*tmp_62) + tmp_44*(tmp_10*tmp_61 + tmp_12*tmp_62 + tmp_63*tmp_8)) + 0.020848748529055869*tmp_45*(tmp_42*(tmp_0*tmp_67 + tmp_1*tmp_68 + tmp_4*tmp_69) + tmp_43*(tmp_11*tmp_67 + tmp_2*tmp_69 + tmp_5*tmp_68) + tmp_44*(tmp_10*tmp_67 + tmp_12*tmp_68 + tmp_69*tmp_8)) + 0.019202922745021479*tmp_45*(tmp_42*(tmp_0*tmp_73 + tmp_1*tmp_74 + tmp_4*tmp_75) + tmp_43*(tmp_11*tmp_73 + tmp_2*tmp_75 + tmp_5*tmp_74) + tmp_44*(tmp_10*tmp_73 + tmp_12*tmp_74 + tmp_75*tmp_8)) + 0.020848748529055869*tmp_45*(tmp_42*(tmp_0*tmp_79 + tmp_1*tmp_80 + tmp_4*tmp_81) + tmp_43*(tmp_11*tmp_79 + tmp_2*tmp_81 + tmp_5*tmp_80) + tmp_44*(tmp_10*tmp_79 + tmp_12*tmp_80 + tmp_8*tmp_81)) + 0.019202922745021479*tmp_45*(tmp_42*(tmp_0*tmp_85 + tmp_1*tmp_86 + tmp_4*tmp_87) + tmp_43*(tmp_11*tmp_85 + tmp_2*tmp_87 + tmp_5*tmp_86) + tmp_44*(tmp_10*tmp_85 + tmp_12*tmp_86 + tmp_8*tmp_87)) + 0.020848748529055869*tmp_45*(tmp_42*(tmp_0*tmp_91 + tmp_1*tmp_92 + tmp_4*tmp_93) + tmp_43*(tmp_11*tmp_91 + tmp_2*tmp_93 + tmp_5*tmp_92) + tmp_44*(tmp_10*tmp_91 + tmp_12*tmp_92 + tmp_8*tmp_93)) + 0.019202922745021479*tmp_45*(tmp_42*(tmp_0*tmp_97 + tmp_1*tmp_98 + tmp_4*tmp_99) + tmp_43*(tmp_11*tmp_97 + tmp_2*tmp_99 + tmp_5*tmp_98) + tmp_44*(tmp_10*tmp_97 + tmp_12*tmp_98 + tmp_8*tmp_99));
      elMat( 0, 0) = a_0_0;
   }




void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                                        const std::vector< Point3D >& coordsElementOuter,
                                                        const std::vector< Point3D >& coordsFacet,
                                                        const Point3D&,
                                                        const Point3D&,
                                                        const Point3D&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                  MatrixXr&                            elMat ) const override
   {
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElementInner[0]( 0 );
   const auto p_affine_0_1 = coordsElementInner[0]( 1 );
   const auto p_affine_0_2 = coordsElementInner[0]( 2 );

   const auto p_affine_1_0 = coordsElementInner[1]( 0 );
   const auto p_affine_1_1 = coordsElementInner[1]( 1 );
   const auto p_affine_1_2 = coordsElementInner[1]( 2 );

   const auto p_affine_2_0 = coordsElementInner[2]( 0 );
   const auto p_affine_2_1 = coordsElementInner[2]( 1 );
   const auto p_affine_2_2 = coordsElementInner[2]( 2 );

   const auto p_affine_3_0 = coordsElementInner[3]( 0 );
   const auto p_affine_3_1 = coordsElementInner[3]( 1 );
   const auto p_affine_3_2 = coordsElementInner[3]( 2 );

   const auto p_affine_4_0 = coordsElementOuter[0]( 0 );
   const auto p_affine_4_1 = coordsElementOuter[0]( 1 );
   const auto p_affine_4_2 = coordsElementOuter[0]( 2 );

   const auto p_affine_5_0 = coordsElementOuter[1]( 0 );
   const auto p_affine_5_1 = coordsElementOuter[1]( 1 );
   const auto p_affine_5_2 = coordsElementOuter[1]( 2 );

   const auto p_affine_6_0 = coordsElementOuter[2]( 0 );
   const auto p_affine_6_1 = coordsElementOuter[2]( 1 );
   const auto p_affine_6_2 = coordsElementOuter[2]( 2 );

   const auto p_affine_7_0 = coordsElementOuter[3]( 0 );
   const auto p_affine_7_1 = coordsElementOuter[3]( 1 );
   const auto p_affine_7_2 = coordsElementOuter[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t tmp_0 = -p_affine_4_0 + p_affine_5_0;
      real_t tmp_1 = -p_affine_4_0 + p_affine_6_0;
      real_t tmp_2 = -p_affine_4_1 + p_affine_7_1;
      real_t tmp_3 = tmp_1*tmp_2;
      real_t tmp_4 = -p_affine_4_0 + p_affine_7_0;
      real_t tmp_5 = -p_affine_4_1 + p_affine_6_1;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = tmp_3 - tmp_6;
      real_t tmp_8 = -p_affine_4_2 + p_affine_7_2;
      real_t tmp_9 = tmp_5*tmp_8;
      real_t tmp_10 = -p_affine_4_2 + p_affine_5_2;
      real_t tmp_11 = -p_affine_4_1 + p_affine_5_1;
      real_t tmp_12 = -p_affine_4_2 + p_affine_6_2;
      real_t tmp_13 = tmp_12*tmp_2;
      real_t tmp_14 = tmp_1*tmp_8;
      real_t tmp_15 = 1.0 / (-tmp_0*tmp_13 + tmp_0*tmp_9 + tmp_10*tmp_3 - tmp_10*tmp_6 + tmp_11*tmp_12*tmp_4 - tmp_11*tmp_14);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_4_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.031405749086161582*tmp_17 + 0.93718850182767688*tmp_18 + tmp_19);
      real_t tmp_21 = tmp_12*tmp_4 - tmp_14;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_4_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.031405749086161582*tmp_23 + 0.93718850182767688*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_13 + tmp_9;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_4_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_7 + tmp_21*tmp_26 + tmp_27*tmp_32 - 1.0/4.0;
      real_t tmp_34 = -tmp_0*tmp_2 + tmp_11*tmp_4;
      real_t tmp_35 = tmp_0*tmp_8 - tmp_10*tmp_4;
      real_t tmp_36 = tmp_10*tmp_2 - tmp_11*tmp_8;
      real_t tmp_37 = tmp_20*tmp_34 + tmp_26*tmp_35 + tmp_32*tmp_36 - 1.0/4.0;
      real_t tmp_38 = tmp_0*tmp_5 - tmp_1*tmp_11;
      real_t tmp_39 = -tmp_0*tmp_12 + tmp_1*tmp_10;
      real_t tmp_40 = -tmp_10*tmp_5 + tmp_11*tmp_12;
      real_t tmp_41 = tmp_20*tmp_38 + tmp_26*tmp_39 + tmp_32*tmp_40 - 1.0/4.0;
      real_t tmp_42 = 0.5*p_affine_13_0;
      real_t tmp_43 = 0.5*p_affine_13_1;
      real_t tmp_44 = 0.5*p_affine_13_2;
      real_t tmp_45 = 1.0*std::pow((std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)), 1.0/2.0);
      real_t tmp_46 = tmp_15*(0.19601935860219369*tmp_17 + 0.60796128279561268*tmp_18 + tmp_19);
      real_t tmp_47 = tmp_15*(0.19601935860219369*tmp_23 + 0.60796128279561268*tmp_24 + tmp_25);
      real_t tmp_48 = tmp_15*(0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_30 + tmp_31);
      real_t tmp_49 = tmp_21*tmp_47 + tmp_27*tmp_48 + tmp_46*tmp_7 - 1.0/4.0;
      real_t tmp_50 = tmp_34*tmp_46 + tmp_35*tmp_47 + tmp_36*tmp_48 - 1.0/4.0;
      real_t tmp_51 = tmp_38*tmp_46 + tmp_39*tmp_47 + tmp_40*tmp_48 - 1.0/4.0;
      real_t tmp_52 = tmp_15*(0.37605877282253791*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_53 = tmp_15*(0.37605877282253791*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_54 = tmp_15*(0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_55 = tmp_21*tmp_53 + tmp_27*tmp_54 + tmp_52*tmp_7 - 1.0/4.0;
      real_t tmp_56 = tmp_34*tmp_52 + tmp_35*tmp_53 + tmp_36*tmp_54 - 1.0/4.0;
      real_t tmp_57 = tmp_38*tmp_52 + tmp_39*tmp_53 + tmp_40*tmp_54 - 1.0/4.0;
      real_t tmp_58 = tmp_15*(0.78764240869137092*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_59 = tmp_15*(0.78764240869137092*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_60 = tmp_15*(0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_61 = tmp_21*tmp_59 + tmp_27*tmp_60 + tmp_58*tmp_7 - 1.0/4.0;
      real_t tmp_62 = tmp_34*tmp_58 + tmp_35*tmp_59 + tmp_36*tmp_60 - 1.0/4.0;
      real_t tmp_63 = tmp_38*tmp_58 + tmp_39*tmp_59 + tmp_40*tmp_60 - 1.0/4.0;
      real_t tmp_64 = tmp_15*(0.58463275527740355*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_65 = tmp_15*(0.58463275527740355*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_66 = tmp_15*(0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_67 = tmp_21*tmp_65 + tmp_27*tmp_66 + tmp_64*tmp_7 - 1.0/4.0;
      real_t tmp_68 = tmp_34*tmp_64 + tmp_35*tmp_65 + tmp_36*tmp_66 - 1.0/4.0;
      real_t tmp_69 = tmp_38*tmp_64 + tmp_39*tmp_65 + tmp_40*tmp_66 - 1.0/4.0;
      real_t tmp_70 = tmp_15*(0.041227165399737475*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_71 = tmp_15*(0.041227165399737475*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_72 = tmp_15*(0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_73 = tmp_21*tmp_71 + tmp_27*tmp_72 + tmp_7*tmp_70 - 1.0/4.0;
      real_t tmp_74 = tmp_34*tmp_70 + tmp_35*tmp_71 + tmp_36*tmp_72 - 1.0/4.0;
      real_t tmp_75 = tmp_38*tmp_70 + tmp_39*tmp_71 + tmp_40*tmp_72 - 1.0/4.0;
      real_t tmp_76 = tmp_15*(0.039308471900058539*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_77 = tmp_15*(0.039308471900058539*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_78 = tmp_15*(0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_79 = tmp_21*tmp_77 + tmp_27*tmp_78 + tmp_7*tmp_76 - 1.0/4.0;
      real_t tmp_80 = tmp_34*tmp_76 + tmp_35*tmp_77 + tmp_36*tmp_78 - 1.0/4.0;
      real_t tmp_81 = tmp_38*tmp_76 + tmp_39*tmp_77 + tmp_40*tmp_78 - 1.0/4.0;
      real_t tmp_82 = tmp_15*(0.78764240869137092*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_83 = tmp_15*(0.78764240869137092*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_84 = tmp_15*(0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_85 = tmp_21*tmp_83 + tmp_27*tmp_84 + tmp_7*tmp_82 - 1.0/4.0;
      real_t tmp_86 = tmp_34*tmp_82 + tmp_35*tmp_83 + tmp_36*tmp_84 - 1.0/4.0;
      real_t tmp_87 = tmp_38*tmp_82 + tmp_39*tmp_83 + tmp_40*tmp_84 - 1.0/4.0;
      real_t tmp_88 = tmp_15*(0.58463275527740355*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_89 = tmp_15*(0.58463275527740355*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_90 = tmp_15*(0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_91 = tmp_21*tmp_89 + tmp_27*tmp_90 + tmp_7*tmp_88 - 1.0/4.0;
      real_t tmp_92 = tmp_34*tmp_88 + tmp_35*tmp_89 + tmp_36*tmp_90 - 1.0/4.0;
      real_t tmp_93 = tmp_38*tmp_88 + tmp_39*tmp_89 + tmp_40*tmp_90 - 1.0/4.0;
      real_t tmp_94 = tmp_15*(0.1711304259088916*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_95 = tmp_15*(0.1711304259088916*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_96 = tmp_15*(0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_97 = tmp_21*tmp_95 + tmp_27*tmp_96 + tmp_7*tmp_94 - 1.0/4.0;
      real_t tmp_98 = tmp_34*tmp_94 + tmp_35*tmp_95 + tmp_36*tmp_96 - 1.0/4.0;
      real_t tmp_99 = tmp_38*tmp_94 + tmp_39*tmp_95 + tmp_40*tmp_96 - 1.0/4.0;
      real_t tmp_100 = tmp_15*(0.37605877282253791*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_101 = tmp_15*(0.37605877282253791*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_102 = tmp_15*(0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_103 = tmp_100*tmp_7 + tmp_101*tmp_21 + tmp_102*tmp_27 - 1.0/4.0;
      real_t tmp_104 = tmp_100*tmp_34 + tmp_101*tmp_35 + tmp_102*tmp_36 - 1.0/4.0;
      real_t tmp_105 = tmp_100*tmp_38 + tmp_101*tmp_39 + tmp_102*tmp_40 - 1.0/4.0;
      real_t tmp_106 = tmp_15*(0.041227165399737475*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_107 = tmp_15*(0.041227165399737475*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_108 = tmp_15*(0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_109 = tmp_106*tmp_7 + tmp_107*tmp_21 + tmp_108*tmp_27 - 1.0/4.0;
      real_t tmp_110 = tmp_106*tmp_34 + tmp_107*tmp_35 + tmp_108*tmp_36 - 1.0/4.0;
      real_t tmp_111 = tmp_106*tmp_38 + tmp_107*tmp_39 + tmp_108*tmp_40 - 1.0/4.0;
      real_t tmp_112 = tmp_15*(0.40446199974765351*tmp_17 + 0.19107600050469298*tmp_18 + tmp_19);
      real_t tmp_113 = tmp_15*(0.40446199974765351*tmp_23 + 0.19107600050469298*tmp_24 + tmp_25);
      real_t tmp_114 = tmp_15*(0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_30 + tmp_31);
      real_t tmp_115 = tmp_112*tmp_7 + tmp_113*tmp_21 + tmp_114*tmp_27 - 1.0/4.0;
      real_t tmp_116 = tmp_112*tmp_34 + tmp_113*tmp_35 + tmp_114*tmp_36 - 1.0/4.0;
      real_t tmp_117 = tmp_112*tmp_38 + tmp_113*tmp_39 + tmp_114*tmp_40 - 1.0/4.0;
      real_t tmp_118 = tmp_15*(0.039308471900058539*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_119 = tmp_15*(0.039308471900058539*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_120 = tmp_15*(0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_121 = tmp_118*tmp_7 + tmp_119*tmp_21 + tmp_120*tmp_27 - 1.0/4.0;
      real_t tmp_122 = tmp_118*tmp_34 + tmp_119*tmp_35 + tmp_120*tmp_36 - 1.0/4.0;
      real_t tmp_123 = tmp_118*tmp_38 + tmp_119*tmp_39 + tmp_120*tmp_40 - 1.0/4.0;
      real_t tmp_124 = tmp_15*(0.93718850182767688*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_125 = tmp_15*(0.93718850182767688*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_126 = tmp_15*(0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_127 = tmp_124*tmp_7 + tmp_125*tmp_21 + tmp_126*tmp_27 - 1.0/4.0;
      real_t tmp_128 = tmp_124*tmp_34 + tmp_125*tmp_35 + tmp_126*tmp_36 - 1.0/4.0;
      real_t tmp_129 = tmp_124*tmp_38 + tmp_125*tmp_39 + tmp_126*tmp_40 - 1.0/4.0;
      real_t tmp_130 = tmp_15*(0.60796128279561268*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_131 = tmp_15*(0.60796128279561268*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_132 = tmp_15*(0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_133 = tmp_130*tmp_7 + tmp_131*tmp_21 + tmp_132*tmp_27 - 1.0/4.0;
      real_t tmp_134 = tmp_130*tmp_34 + tmp_131*tmp_35 + tmp_132*tmp_36 - 1.0/4.0;
      real_t tmp_135 = tmp_130*tmp_38 + tmp_131*tmp_39 + tmp_132*tmp_40 - 1.0/4.0;
      real_t tmp_136 = tmp_15*(0.19107600050469298*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_137 = tmp_15*(0.19107600050469298*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_138 = tmp_15*(0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_139 = tmp_136*tmp_7 + tmp_137*tmp_21 + tmp_138*tmp_27 - 1.0/4.0;
      real_t tmp_140 = tmp_136*tmp_34 + tmp_137*tmp_35 + tmp_138*tmp_36 - 1.0/4.0;
      real_t tmp_141 = tmp_136*tmp_38 + tmp_137*tmp_39 + tmp_138*tmp_40 - 1.0/4.0;
      real_t tmp_142 = tmp_15*(0.031405749086161582*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_143 = tmp_15*(0.031405749086161582*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_144 = tmp_15*(0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_145 = tmp_142*tmp_7 + tmp_143*tmp_21 + tmp_144*tmp_27 - 1.0/4.0;
      real_t tmp_146 = tmp_142*tmp_34 + tmp_143*tmp_35 + tmp_144*tmp_36 - 1.0/4.0;
      real_t tmp_147 = tmp_142*tmp_38 + tmp_143*tmp_39 + tmp_144*tmp_40 - 1.0/4.0;
      real_t tmp_148 = tmp_15*(0.19601935860219369*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_149 = tmp_15*(0.19601935860219369*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_150 = tmp_15*(0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_151 = tmp_148*tmp_7 + tmp_149*tmp_21 + tmp_150*tmp_27 - 1.0/4.0;
      real_t tmp_152 = tmp_148*tmp_34 + tmp_149*tmp_35 + tmp_150*tmp_36 - 1.0/4.0;
      real_t tmp_153 = tmp_148*tmp_38 + tmp_149*tmp_39 + tmp_150*tmp_40 - 1.0/4.0;
      real_t tmp_154 = tmp_15*(0.40446199974765351*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_155 = tmp_15*(0.40446199974765351*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_156 = tmp_15*(0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_157 = tmp_154*tmp_7 + tmp_155*tmp_21 + tmp_156*tmp_27 - 1.0/4.0;
      real_t tmp_158 = tmp_154*tmp_34 + tmp_155*tmp_35 + tmp_156*tmp_36 - 1.0/4.0;
      real_t tmp_159 = tmp_154*tmp_38 + tmp_155*tmp_39 + tmp_156*tmp_40 - 1.0/4.0;
      real_t tmp_160 = tmp_15*(0.1711304259088916*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_161 = tmp_15*(0.1711304259088916*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_162 = tmp_15*(0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_163 = tmp_160*tmp_7 + tmp_161*tmp_21 + tmp_162*tmp_27 - 1.0/4.0;
      real_t tmp_164 = tmp_160*tmp_34 + tmp_161*tmp_35 + tmp_162*tmp_36 - 1.0/4.0;
      real_t tmp_165 = tmp_160*tmp_38 + tmp_161*tmp_39 + tmp_162*tmp_40 - 1.0/4.0;
      real_t a_0_0 = 0.020848748529055869*tmp_45*(-tmp_42*(tmp_0*tmp_103 + tmp_1*tmp_104 + tmp_105*tmp_4) - tmp_43*(tmp_103*tmp_11 + tmp_104*tmp_5 + tmp_105*tmp_2) - tmp_44*(tmp_10*tmp_103 + tmp_104*tmp_12 + tmp_105*tmp_8)) + 0.019202922745021479*tmp_45*(-tmp_42*(tmp_0*tmp_109 + tmp_1*tmp_110 + tmp_111*tmp_4) - tmp_43*(tmp_109*tmp_11 + tmp_110*tmp_5 + tmp_111*tmp_2) - tmp_44*(tmp_10*tmp_109 + tmp_110*tmp_12 + tmp_111*tmp_8)) + 0.042507265838595799*tmp_45*(-tmp_42*(tmp_0*tmp_115 + tmp_1*tmp_116 + tmp_117*tmp_4) - tmp_43*(tmp_11*tmp_115 + tmp_116*tmp_5 + tmp_117*tmp_2) - tmp_44*(tmp_10*tmp_115 + tmp_116*tmp_12 + tmp_117*tmp_8)) + 0.020848748529055869*tmp_45*(-tmp_42*(tmp_0*tmp_121 + tmp_1*tmp_122 + tmp_123*tmp_4) - tmp_43*(tmp_11*tmp_121 + tmp_122*tmp_5 + tmp_123*tmp_2) - tmp_44*(tmp_10*tmp_121 + tmp_12*tmp_122 + tmp_123*tmp_8)) + 0.0068572537431980923*tmp_45*(-tmp_42*(tmp_0*tmp_127 + tmp_1*tmp_128 + tmp_129*tmp_4) - tmp_43*(tmp_11*tmp_127 + tmp_128*tmp_5 + tmp_129*tmp_2) - tmp_44*(tmp_10*tmp_127 + tmp_12*tmp_128 + tmp_129*tmp_8)) + 0.037198804536718075*tmp_45*(-tmp_42*(tmp_0*tmp_133 + tmp_1*tmp_134 + tmp_135*tmp_4) - tmp_43*(tmp_11*tmp_133 + tmp_134*tmp_5 + tmp_135*tmp_2) - tmp_44*(tmp_10*tmp_133 + tmp_12*tmp_134 + tmp_135*tmp_8)) + 0.042507265838595799*tmp_45*(-tmp_42*(tmp_0*tmp_139 + tmp_1*tmp_140 + tmp_141*tmp_4) - tmp_43*(tmp_11*tmp_139 + tmp_140*tmp_5 + tmp_141*tmp_2) - tmp_44*(tmp_10*tmp_139 + tmp_12*tmp_140 + tmp_141*tmp_8)) + 0.0068572537431980923*tmp_45*(-tmp_42*(tmp_0*tmp_145 + tmp_1*tmp_146 + tmp_147*tmp_4) - tmp_43*(tmp_11*tmp_145 + tmp_146*tmp_5 + tmp_147*tmp_2) - tmp_44*(tmp_10*tmp_145 + tmp_12*tmp_146 + tmp_147*tmp_8)) + 0.037198804536718075*tmp_45*(-tmp_42*(tmp_0*tmp_151 + tmp_1*tmp_152 + tmp_153*tmp_4) - tmp_43*(tmp_11*tmp_151 + tmp_152*tmp_5 + tmp_153*tmp_2) - tmp_44*(tmp_10*tmp_151 + tmp_12*tmp_152 + tmp_153*tmp_8)) + 0.042507265838595799*tmp_45*(-tmp_42*(tmp_0*tmp_157 + tmp_1*tmp_158 + tmp_159*tmp_4) - tmp_43*(tmp_11*tmp_157 + tmp_158*tmp_5 + tmp_159*tmp_2) - tmp_44*(tmp_10*tmp_157 + tmp_12*tmp_158 + tmp_159*tmp_8)) + 0.019202922745021479*tmp_45*(-tmp_42*(tmp_0*tmp_163 + tmp_1*tmp_164 + tmp_165*tmp_4) - tmp_43*(tmp_11*tmp_163 + tmp_164*tmp_5 + tmp_165*tmp_2) - tmp_44*(tmp_10*tmp_163 + tmp_12*tmp_164 + tmp_165*tmp_8)) + 0.0068572537431980923*tmp_45*(-tmp_42*(tmp_0*tmp_33 + tmp_1*tmp_37 + tmp_4*tmp_41) - tmp_43*(tmp_11*tmp_33 + tmp_2*tmp_41 + tmp_37*tmp_5) - tmp_44*(tmp_10*tmp_33 + tmp_12*tmp_37 + tmp_41*tmp_8)) + 0.037198804536718075*tmp_45*(-tmp_42*(tmp_0*tmp_49 + tmp_1*tmp_50 + tmp_4*tmp_51) - tmp_43*(tmp_11*tmp_49 + tmp_2*tmp_51 + tmp_5*tmp_50) - tmp_44*(tmp_10*tmp_49 + tmp_12*tmp_50 + tmp_51*tmp_8)) + 0.020848748529055869*tmp_45*(-tmp_42*(tmp_0*tmp_55 + tmp_1*tmp_56 + tmp_4*tmp_57) - tmp_43*(tmp_11*tmp_55 + tmp_2*tmp_57 + tmp_5*tmp_56) - tmp_44*(tmp_10*tmp_55 + tmp_12*tmp_56 + tmp_57*tmp_8)) + 0.019202922745021479*tmp_45*(-tmp_42*(tmp_0*tmp_61 + tmp_1*tmp_62 + tmp_4*tmp_63) - tmp_43*(tmp_11*tmp_61 + tmp_2*tmp_63 + tmp_5*tmp_62) - tmp_44*(tmp_10*tmp_61 + tmp_12*tmp_62 + tmp_63*tmp_8)) + 0.020848748529055869*tmp_45*(-tmp_42*(tmp_0*tmp_67 + tmp_1*tmp_68 + tmp_4*tmp_69) - tmp_43*(tmp_11*tmp_67 + tmp_2*tmp_69 + tmp_5*tmp_68) - tmp_44*(tmp_10*tmp_67 + tmp_12*tmp_68 + tmp_69*tmp_8)) + 0.019202922745021479*tmp_45*(-tmp_42*(tmp_0*tmp_73 + tmp_1*tmp_74 + tmp_4*tmp_75) - tmp_43*(tmp_11*tmp_73 + tmp_2*tmp_75 + tmp_5*tmp_74) - tmp_44*(tmp_10*tmp_73 + tmp_12*tmp_74 + tmp_75*tmp_8)) + 0.020848748529055869*tmp_45*(-tmp_42*(tmp_0*tmp_79 + tmp_1*tmp_80 + tmp_4*tmp_81) - tmp_43*(tmp_11*tmp_79 + tmp_2*tmp_81 + tmp_5*tmp_80) - tmp_44*(tmp_10*tmp_79 + tmp_12*tmp_80 + tmp_8*tmp_81)) + 0.019202922745021479*tmp_45*(-tmp_42*(tmp_0*tmp_85 + tmp_1*tmp_86 + tmp_4*tmp_87) - tmp_43*(tmp_11*tmp_85 + tmp_2*tmp_87 + tmp_5*tmp_86) - tmp_44*(tmp_10*tmp_85 + tmp_12*tmp_86 + tmp_8*tmp_87)) + 0.020848748529055869*tmp_45*(-tmp_42*(tmp_0*tmp_91 + tmp_1*tmp_92 + tmp_4*tmp_93) - tmp_43*(tmp_11*tmp_91 + tmp_2*tmp_93 + tmp_5*tmp_92) - tmp_44*(tmp_10*tmp_91 + tmp_12*tmp_92 + tmp_8*tmp_93)) + 0.019202922745021479*tmp_45*(-tmp_42*(tmp_0*tmp_97 + tmp_1*tmp_98 + tmp_4*tmp_99) - tmp_43*(tmp_11*tmp_97 + tmp_2*tmp_99 + tmp_5*tmp_98) - tmp_44*(tmp_10*tmp_97 + tmp_12*tmp_98 + tmp_8*tmp_99));
      elMat( 0, 0) = a_0_0;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Point3D >& coordsElement,
    const std::vector< Point3D >& coordsFacet,
    const Point3D&,
    const Point3D&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
                                        MatrixXr&                            elMat ) const override
{
   elMat.resize( Eigen::Index( testBasis.numDoFsPerElement( 3, uint_c( testDegree ) ) ),
                 Eigen::Index( trialBasis.numDoFsPerElement( 3, uint_c( trialDegree ) ) ) );

   const auto p_affine_0_0 = coordsElement[0]( 0 );
   const auto p_affine_0_1 = coordsElement[0]( 1 );
   const auto p_affine_0_2 = coordsElement[0]( 2 );

   const auto p_affine_1_0 = coordsElement[1]( 0 );
   const auto p_affine_1_1 = coordsElement[1]( 1 );
   const auto p_affine_1_2 = coordsElement[1]( 2 );

   const auto p_affine_2_0 = coordsElement[2]( 0 );
   const auto p_affine_2_1 = coordsElement[2]( 1 );
   const auto p_affine_2_2 = coordsElement[2]( 2 );

   const auto p_affine_3_0 = coordsElement[3]( 0 );
   const auto p_affine_3_1 = coordsElement[3]( 1 );
   const auto p_affine_3_2 = coordsElement[3]( 2 );

   const auto p_affine_8_0 = coordsFacet[0]( 0 );
   const auto p_affine_8_1 = coordsFacet[0]( 1 );
   const auto p_affine_8_2 = coordsFacet[0]( 2 );

   const auto p_affine_9_0 = coordsFacet[1]( 0 );
   const auto p_affine_9_1 = coordsFacet[1]( 1 );
   const auto p_affine_9_2 = coordsFacet[1]( 2 );

   const auto p_affine_10_0 = coordsFacet[2]( 0 );
   const auto p_affine_10_1 = coordsFacet[2]( 1 );
   const auto p_affine_10_2 = coordsFacet[2]( 2 );

   const auto p_affine_13_0 = outwardNormal( 0 );
   const auto p_affine_13_1 = outwardNormal( 1 );
   const auto p_affine_13_2 = outwardNormal( 2 );


      real_t tmp_0 = -p_affine_0_0 + p_affine_1_0;
      real_t tmp_1 = -p_affine_0_0 + p_affine_2_0;
      real_t tmp_2 = -p_affine_0_1 + p_affine_3_1;
      real_t tmp_3 = tmp_1*tmp_2;
      real_t tmp_4 = -p_affine_0_0 + p_affine_3_0;
      real_t tmp_5 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = tmp_3 - tmp_6;
      real_t tmp_8 = -p_affine_0_2 + p_affine_3_2;
      real_t tmp_9 = tmp_5*tmp_8;
      real_t tmp_10 = -p_affine_0_2 + p_affine_1_2;
      real_t tmp_11 = -p_affine_0_1 + p_affine_1_1;
      real_t tmp_12 = -p_affine_0_2 + p_affine_2_2;
      real_t tmp_13 = tmp_12*tmp_2;
      real_t tmp_14 = tmp_1*tmp_8;
      real_t tmp_15 = 1.0 / (-tmp_0*tmp_13 + tmp_0*tmp_9 + tmp_10*tmp_3 - tmp_10*tmp_6 + tmp_11*tmp_12*tmp_4 - tmp_11*tmp_14);
      real_t tmp_16 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_17 = -tmp_16;
      real_t tmp_18 = p_affine_10_2 - p_affine_8_2;
      real_t tmp_19 = -p_affine_0_2 + p_affine_8_2;
      real_t tmp_20 = tmp_15*(0.031405749086161582*tmp_17 + 0.93718850182767688*tmp_18 + tmp_19);
      real_t tmp_21 = tmp_12*tmp_4 - tmp_14;
      real_t tmp_22 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_23 = -tmp_22;
      real_t tmp_24 = p_affine_10_1 - p_affine_8_1;
      real_t tmp_25 = -p_affine_0_1 + p_affine_8_1;
      real_t tmp_26 = tmp_15*(0.031405749086161582*tmp_23 + 0.93718850182767688*tmp_24 + tmp_25);
      real_t tmp_27 = -tmp_13 + tmp_9;
      real_t tmp_28 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = p_affine_10_0 - p_affine_8_0;
      real_t tmp_31 = -p_affine_0_0 + p_affine_8_0;
      real_t tmp_32 = tmp_15*(0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_30 + tmp_31);
      real_t tmp_33 = tmp_20*tmp_7 + tmp_21*tmp_26 + tmp_27*tmp_32 - 1.0/4.0;
      real_t tmp_34 = -tmp_0*tmp_2 + tmp_11*tmp_4;
      real_t tmp_35 = tmp_0*tmp_8 - tmp_10*tmp_4;
      real_t tmp_36 = tmp_10*tmp_2 - tmp_11*tmp_8;
      real_t tmp_37 = tmp_20*tmp_34 + tmp_26*tmp_35 + tmp_32*tmp_36 - 1.0/4.0;
      real_t tmp_38 = tmp_0*tmp_5 - tmp_1*tmp_11;
      real_t tmp_39 = -tmp_0*tmp_12 + tmp_1*tmp_10;
      real_t tmp_40 = -tmp_10*tmp_5 + tmp_11*tmp_12;
      real_t tmp_41 = tmp_20*tmp_38 + tmp_26*tmp_39 + tmp_32*tmp_40 - 1.0/4.0;
      real_t tmp_42 = 1.0*std::pow((std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)*std::abs(tmp_16*tmp_24 - tmp_18*tmp_22)) + (std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)*std::abs(tmp_16*tmp_30 - tmp_18*tmp_28)) + (std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)*std::abs(tmp_22*tmp_30 - tmp_24*tmp_28)), 1.0/2.0);
      real_t tmp_43 = tmp_15*(0.19601935860219369*tmp_17 + 0.60796128279561268*tmp_18 + tmp_19);
      real_t tmp_44 = tmp_15*(0.19601935860219369*tmp_23 + 0.60796128279561268*tmp_24 + tmp_25);
      real_t tmp_45 = tmp_15*(0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_30 + tmp_31);
      real_t tmp_46 = tmp_21*tmp_44 + tmp_27*tmp_45 + tmp_43*tmp_7 - 1.0/4.0;
      real_t tmp_47 = tmp_34*tmp_43 + tmp_35*tmp_44 + tmp_36*tmp_45 - 1.0/4.0;
      real_t tmp_48 = tmp_38*tmp_43 + tmp_39*tmp_44 + tmp_40*tmp_45 - 1.0/4.0;
      real_t tmp_49 = tmp_15*(0.37605877282253791*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_50 = tmp_15*(0.37605877282253791*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_51 = tmp_15*(0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_52 = tmp_21*tmp_50 + tmp_27*tmp_51 + tmp_49*tmp_7 - 1.0/4.0;
      real_t tmp_53 = tmp_34*tmp_49 + tmp_35*tmp_50 + tmp_36*tmp_51 - 1.0/4.0;
      real_t tmp_54 = tmp_38*tmp_49 + tmp_39*tmp_50 + tmp_40*tmp_51 - 1.0/4.0;
      real_t tmp_55 = tmp_15*(0.78764240869137092*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_56 = tmp_15*(0.78764240869137092*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_57 = tmp_15*(0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_58 = tmp_21*tmp_56 + tmp_27*tmp_57 + tmp_55*tmp_7 - 1.0/4.0;
      real_t tmp_59 = tmp_34*tmp_55 + tmp_35*tmp_56 + tmp_36*tmp_57 - 1.0/4.0;
      real_t tmp_60 = tmp_38*tmp_55 + tmp_39*tmp_56 + tmp_40*tmp_57 - 1.0/4.0;
      real_t tmp_61 = tmp_15*(0.58463275527740355*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_62 = tmp_15*(0.58463275527740355*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_63 = tmp_15*(0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_64 = tmp_21*tmp_62 + tmp_27*tmp_63 + tmp_61*tmp_7 - 1.0/4.0;
      real_t tmp_65 = tmp_34*tmp_61 + tmp_35*tmp_62 + tmp_36*tmp_63 - 1.0/4.0;
      real_t tmp_66 = tmp_38*tmp_61 + tmp_39*tmp_62 + tmp_40*tmp_63 - 1.0/4.0;
      real_t tmp_67 = tmp_15*(0.041227165399737475*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_68 = tmp_15*(0.041227165399737475*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_69 = tmp_15*(0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_70 = tmp_21*tmp_68 + tmp_27*tmp_69 + tmp_67*tmp_7 - 1.0/4.0;
      real_t tmp_71 = tmp_34*tmp_67 + tmp_35*tmp_68 + tmp_36*tmp_69 - 1.0/4.0;
      real_t tmp_72 = tmp_38*tmp_67 + tmp_39*tmp_68 + tmp_40*tmp_69 - 1.0/4.0;
      real_t tmp_73 = tmp_15*(0.039308471900058539*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_74 = tmp_15*(0.039308471900058539*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_75 = tmp_15*(0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_76 = tmp_21*tmp_74 + tmp_27*tmp_75 + tmp_7*tmp_73 - 1.0/4.0;
      real_t tmp_77 = tmp_34*tmp_73 + tmp_35*tmp_74 + tmp_36*tmp_75 - 1.0/4.0;
      real_t tmp_78 = tmp_38*tmp_73 + tmp_39*tmp_74 + tmp_40*tmp_75 - 1.0/4.0;
      real_t tmp_79 = tmp_15*(0.78764240869137092*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_80 = tmp_15*(0.78764240869137092*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_81 = tmp_15*(0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_82 = tmp_21*tmp_80 + tmp_27*tmp_81 + tmp_7*tmp_79 - 1.0/4.0;
      real_t tmp_83 = tmp_34*tmp_79 + tmp_35*tmp_80 + tmp_36*tmp_81 - 1.0/4.0;
      real_t tmp_84 = tmp_38*tmp_79 + tmp_39*tmp_80 + tmp_40*tmp_81 - 1.0/4.0;
      real_t tmp_85 = tmp_15*(0.58463275527740355*tmp_17 + 0.039308471900058539*tmp_18 + tmp_19);
      real_t tmp_86 = tmp_15*(0.58463275527740355*tmp_23 + 0.039308471900058539*tmp_24 + tmp_25);
      real_t tmp_87 = tmp_15*(0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_30 + tmp_31);
      real_t tmp_88 = tmp_21*tmp_86 + tmp_27*tmp_87 + tmp_7*tmp_85 - 1.0/4.0;
      real_t tmp_89 = tmp_34*tmp_85 + tmp_35*tmp_86 + tmp_36*tmp_87 - 1.0/4.0;
      real_t tmp_90 = tmp_38*tmp_85 + tmp_39*tmp_86 + tmp_40*tmp_87 - 1.0/4.0;
      real_t tmp_91 = tmp_15*(0.1711304259088916*tmp_17 + 0.78764240869137092*tmp_18 + tmp_19);
      real_t tmp_92 = tmp_15*(0.1711304259088916*tmp_23 + 0.78764240869137092*tmp_24 + tmp_25);
      real_t tmp_93 = tmp_15*(0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_30 + tmp_31);
      real_t tmp_94 = tmp_21*tmp_92 + tmp_27*tmp_93 + tmp_7*tmp_91 - 1.0/4.0;
      real_t tmp_95 = tmp_34*tmp_91 + tmp_35*tmp_92 + tmp_36*tmp_93 - 1.0/4.0;
      real_t tmp_96 = tmp_38*tmp_91 + tmp_39*tmp_92 + tmp_40*tmp_93 - 1.0/4.0;
      real_t tmp_97 = tmp_15*(0.37605877282253791*tmp_17 + 0.58463275527740355*tmp_18 + tmp_19);
      real_t tmp_98 = tmp_15*(0.37605877282253791*tmp_23 + 0.58463275527740355*tmp_24 + tmp_25);
      real_t tmp_99 = tmp_15*(0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_30 + tmp_31);
      real_t tmp_100 = tmp_21*tmp_98 + tmp_27*tmp_99 + tmp_7*tmp_97 - 1.0/4.0;
      real_t tmp_101 = tmp_34*tmp_97 + tmp_35*tmp_98 + tmp_36*tmp_99 - 1.0/4.0;
      real_t tmp_102 = tmp_38*tmp_97 + tmp_39*tmp_98 + tmp_40*tmp_99 - 1.0/4.0;
      real_t tmp_103 = tmp_15*(0.041227165399737475*tmp_17 + 0.1711304259088916*tmp_18 + tmp_19);
      real_t tmp_104 = tmp_15*(0.041227165399737475*tmp_23 + 0.1711304259088916*tmp_24 + tmp_25);
      real_t tmp_105 = tmp_15*(0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_30 + tmp_31);
      real_t tmp_106 = tmp_103*tmp_7 + tmp_104*tmp_21 + tmp_105*tmp_27 - 1.0/4.0;
      real_t tmp_107 = tmp_103*tmp_34 + tmp_104*tmp_35 + tmp_105*tmp_36 - 1.0/4.0;
      real_t tmp_108 = tmp_103*tmp_38 + tmp_104*tmp_39 + tmp_105*tmp_40 - 1.0/4.0;
      real_t tmp_109 = tmp_15*(0.40446199974765351*tmp_17 + 0.19107600050469298*tmp_18 + tmp_19);
      real_t tmp_110 = tmp_15*(0.40446199974765351*tmp_23 + 0.19107600050469298*tmp_24 + tmp_25);
      real_t tmp_111 = tmp_15*(0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_30 + tmp_31);
      real_t tmp_112 = tmp_109*tmp_7 + tmp_110*tmp_21 + tmp_111*tmp_27 - 1.0/4.0;
      real_t tmp_113 = tmp_109*tmp_34 + tmp_110*tmp_35 + tmp_111*tmp_36 - 1.0/4.0;
      real_t tmp_114 = tmp_109*tmp_38 + tmp_110*tmp_39 + tmp_111*tmp_40 - 1.0/4.0;
      real_t tmp_115 = tmp_15*(0.039308471900058539*tmp_17 + 0.37605877282253791*tmp_18 + tmp_19);
      real_t tmp_116 = tmp_15*(0.039308471900058539*tmp_23 + 0.37605877282253791*tmp_24 + tmp_25);
      real_t tmp_117 = tmp_15*(0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_30 + tmp_31);
      real_t tmp_118 = tmp_115*tmp_7 + tmp_116*tmp_21 + tmp_117*tmp_27 - 1.0/4.0;
      real_t tmp_119 = tmp_115*tmp_34 + tmp_116*tmp_35 + tmp_117*tmp_36 - 1.0/4.0;
      real_t tmp_120 = tmp_115*tmp_38 + tmp_116*tmp_39 + tmp_117*tmp_40 - 1.0/4.0;
      real_t tmp_121 = tmp_15*(0.93718850182767688*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_122 = tmp_15*(0.93718850182767688*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_123 = tmp_15*(0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_124 = tmp_121*tmp_7 + tmp_122*tmp_21 + tmp_123*tmp_27 - 1.0/4.0;
      real_t tmp_125 = tmp_121*tmp_34 + tmp_122*tmp_35 + tmp_123*tmp_36 - 1.0/4.0;
      real_t tmp_126 = tmp_121*tmp_38 + tmp_122*tmp_39 + tmp_123*tmp_40 - 1.0/4.0;
      real_t tmp_127 = tmp_15*(0.60796128279561268*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_128 = tmp_15*(0.60796128279561268*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_129 = tmp_15*(0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_130 = tmp_127*tmp_7 + tmp_128*tmp_21 + tmp_129*tmp_27 - 1.0/4.0;
      real_t tmp_131 = tmp_127*tmp_34 + tmp_128*tmp_35 + tmp_129*tmp_36 - 1.0/4.0;
      real_t tmp_132 = tmp_127*tmp_38 + tmp_128*tmp_39 + tmp_129*tmp_40 - 1.0/4.0;
      real_t tmp_133 = tmp_15*(0.19107600050469298*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_134 = tmp_15*(0.19107600050469298*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_135 = tmp_15*(0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_136 = tmp_133*tmp_7 + tmp_134*tmp_21 + tmp_135*tmp_27 - 1.0/4.0;
      real_t tmp_137 = tmp_133*tmp_34 + tmp_134*tmp_35 + tmp_135*tmp_36 - 1.0/4.0;
      real_t tmp_138 = tmp_133*tmp_38 + tmp_134*tmp_39 + tmp_135*tmp_40 - 1.0/4.0;
      real_t tmp_139 = tmp_15*(0.031405749086161582*tmp_17 + 0.031405749086161582*tmp_18 + tmp_19);
      real_t tmp_140 = tmp_15*(0.031405749086161582*tmp_23 + 0.031405749086161582*tmp_24 + tmp_25);
      real_t tmp_141 = tmp_15*(0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_30 + tmp_31);
      real_t tmp_142 = tmp_139*tmp_7 + tmp_140*tmp_21 + tmp_141*tmp_27 - 1.0/4.0;
      real_t tmp_143 = tmp_139*tmp_34 + tmp_140*tmp_35 + tmp_141*tmp_36 - 1.0/4.0;
      real_t tmp_144 = tmp_139*tmp_38 + tmp_140*tmp_39 + tmp_141*tmp_40 - 1.0/4.0;
      real_t tmp_145 = tmp_15*(0.19601935860219369*tmp_17 + 0.19601935860219369*tmp_18 + tmp_19);
      real_t tmp_146 = tmp_15*(0.19601935860219369*tmp_23 + 0.19601935860219369*tmp_24 + tmp_25);
      real_t tmp_147 = tmp_15*(0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_30 + tmp_31);
      real_t tmp_148 = tmp_145*tmp_7 + tmp_146*tmp_21 + tmp_147*tmp_27 - 1.0/4.0;
      real_t tmp_149 = tmp_145*tmp_34 + tmp_146*tmp_35 + tmp_147*tmp_36 - 1.0/4.0;
      real_t tmp_150 = tmp_145*tmp_38 + tmp_146*tmp_39 + tmp_147*tmp_40 - 1.0/4.0;
      real_t tmp_151 = tmp_15*(0.40446199974765351*tmp_17 + 0.40446199974765351*tmp_18 + tmp_19);
      real_t tmp_152 = tmp_15*(0.40446199974765351*tmp_23 + 0.40446199974765351*tmp_24 + tmp_25);
      real_t tmp_153 = tmp_15*(0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_30 + tmp_31);
      real_t tmp_154 = tmp_151*tmp_7 + tmp_152*tmp_21 + tmp_153*tmp_27 - 1.0/4.0;
      real_t tmp_155 = tmp_151*tmp_34 + tmp_152*tmp_35 + tmp_153*tmp_36 - 1.0/4.0;
      real_t tmp_156 = tmp_151*tmp_38 + tmp_152*tmp_39 + tmp_153*tmp_40 - 1.0/4.0;
      real_t tmp_157 = tmp_15*(0.1711304259088916*tmp_17 + 0.041227165399737475*tmp_18 + tmp_19);
      real_t tmp_158 = tmp_15*(0.1711304259088916*tmp_23 + 0.041227165399737475*tmp_24 + tmp_25);
      real_t tmp_159 = tmp_15*(0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_30 + tmp_31);
      real_t tmp_160 = tmp_157*tmp_7 + tmp_158*tmp_21 + tmp_159*tmp_27 - 1.0/4.0;
      real_t tmp_161 = tmp_157*tmp_34 + tmp_158*tmp_35 + tmp_159*tmp_36 - 1.0/4.0;
      real_t tmp_162 = tmp_157*tmp_38 + tmp_158*tmp_39 + tmp_159*tmp_40 - 1.0/4.0;
      real_t a_0_0 = 0.020848748529055869*tmp_42*(p_affine_13_0*(tmp_0*tmp_100 + tmp_1*tmp_101 + tmp_102*tmp_4) + p_affine_13_1*(tmp_100*tmp_11 + tmp_101*tmp_5 + tmp_102*tmp_2) + p_affine_13_2*(tmp_10*tmp_100 + tmp_101*tmp_12 + tmp_102*tmp_8)) + 0.019202922745021479*tmp_42*(p_affine_13_0*(tmp_0*tmp_106 + tmp_1*tmp_107 + tmp_108*tmp_4) + p_affine_13_1*(tmp_106*tmp_11 + tmp_107*tmp_5 + tmp_108*tmp_2) + p_affine_13_2*(tmp_10*tmp_106 + tmp_107*tmp_12 + tmp_108*tmp_8)) + 0.042507265838595799*tmp_42*(p_affine_13_0*(tmp_0*tmp_112 + tmp_1*tmp_113 + tmp_114*tmp_4) + p_affine_13_1*(tmp_11*tmp_112 + tmp_113*tmp_5 + tmp_114*tmp_2) + p_affine_13_2*(tmp_10*tmp_112 + tmp_113*tmp_12 + tmp_114*tmp_8)) + 0.020848748529055869*tmp_42*(p_affine_13_0*(tmp_0*tmp_118 + tmp_1*tmp_119 + tmp_120*tmp_4) + p_affine_13_1*(tmp_11*tmp_118 + tmp_119*tmp_5 + tmp_120*tmp_2) + p_affine_13_2*(tmp_10*tmp_118 + tmp_119*tmp_12 + tmp_120*tmp_8)) + 0.0068572537431980923*tmp_42*(p_affine_13_0*(tmp_0*tmp_124 + tmp_1*tmp_125 + tmp_126*tmp_4) + p_affine_13_1*(tmp_11*tmp_124 + tmp_125*tmp_5 + tmp_126*tmp_2) + p_affine_13_2*(tmp_10*tmp_124 + tmp_12*tmp_125 + tmp_126*tmp_8)) + 0.037198804536718075*tmp_42*(p_affine_13_0*(tmp_0*tmp_130 + tmp_1*tmp_131 + tmp_132*tmp_4) + p_affine_13_1*(tmp_11*tmp_130 + tmp_131*tmp_5 + tmp_132*tmp_2) + p_affine_13_2*(tmp_10*tmp_130 + tmp_12*tmp_131 + tmp_132*tmp_8)) + 0.042507265838595799*tmp_42*(p_affine_13_0*(tmp_0*tmp_136 + tmp_1*tmp_137 + tmp_138*tmp_4) + p_affine_13_1*(tmp_11*tmp_136 + tmp_137*tmp_5 + tmp_138*tmp_2) + p_affine_13_2*(tmp_10*tmp_136 + tmp_12*tmp_137 + tmp_138*tmp_8)) + 0.0068572537431980923*tmp_42*(p_affine_13_0*(tmp_0*tmp_142 + tmp_1*tmp_143 + tmp_144*tmp_4) + p_affine_13_1*(tmp_11*tmp_142 + tmp_143*tmp_5 + tmp_144*tmp_2) + p_affine_13_2*(tmp_10*tmp_142 + tmp_12*tmp_143 + tmp_144*tmp_8)) + 0.037198804536718075*tmp_42*(p_affine_13_0*(tmp_0*tmp_148 + tmp_1*tmp_149 + tmp_150*tmp_4) + p_affine_13_1*(tmp_11*tmp_148 + tmp_149*tmp_5 + tmp_150*tmp_2) + p_affine_13_2*(tmp_10*tmp_148 + tmp_12*tmp_149 + tmp_150*tmp_8)) + 0.042507265838595799*tmp_42*(p_affine_13_0*(tmp_0*tmp_154 + tmp_1*tmp_155 + tmp_156*tmp_4) + p_affine_13_1*(tmp_11*tmp_154 + tmp_155*tmp_5 + tmp_156*tmp_2) + p_affine_13_2*(tmp_10*tmp_154 + tmp_12*tmp_155 + tmp_156*tmp_8)) + 0.019202922745021479*tmp_42*(p_affine_13_0*(tmp_0*tmp_160 + tmp_1*tmp_161 + tmp_162*tmp_4) + p_affine_13_1*(tmp_11*tmp_160 + tmp_161*tmp_5 + tmp_162*tmp_2) + p_affine_13_2*(tmp_10*tmp_160 + tmp_12*tmp_161 + tmp_162*tmp_8)) + 0.0068572537431980923*tmp_42*(p_affine_13_0*(tmp_0*tmp_33 + tmp_1*tmp_37 + tmp_4*tmp_41) + p_affine_13_1*(tmp_11*tmp_33 + tmp_2*tmp_41 + tmp_37*tmp_5) + p_affine_13_2*(tmp_10*tmp_33 + tmp_12*tmp_37 + tmp_41*tmp_8)) + 0.037198804536718075*tmp_42*(p_affine_13_0*(tmp_0*tmp_46 + tmp_1*tmp_47 + tmp_4*tmp_48) + p_affine_13_1*(tmp_11*tmp_46 + tmp_2*tmp_48 + tmp_47*tmp_5) + p_affine_13_2*(tmp_10*tmp_46 + tmp_12*tmp_47 + tmp_48*tmp_8)) + 0.020848748529055869*tmp_42*(p_affine_13_0*(tmp_0*tmp_52 + tmp_1*tmp_53 + tmp_4*tmp_54) + p_affine_13_1*(tmp_11*tmp_52 + tmp_2*tmp_54 + tmp_5*tmp_53) + p_affine_13_2*(tmp_10*tmp_52 + tmp_12*tmp_53 + tmp_54*tmp_8)) + 0.019202922745021479*tmp_42*(p_affine_13_0*(tmp_0*tmp_58 + tmp_1*tmp_59 + tmp_4*tmp_60) + p_affine_13_1*(tmp_11*tmp_58 + tmp_2*tmp_60 + tmp_5*tmp_59) + p_affine_13_2*(tmp_10*tmp_58 + tmp_12*tmp_59 + tmp_60*tmp_8)) + 0.020848748529055869*tmp_42*(p_affine_13_0*(tmp_0*tmp_64 + tmp_1*tmp_65 + tmp_4*tmp_66) + p_affine_13_1*(tmp_11*tmp_64 + tmp_2*tmp_66 + tmp_5*tmp_65) + p_affine_13_2*(tmp_10*tmp_64 + tmp_12*tmp_65 + tmp_66*tmp_8)) + 0.019202922745021479*tmp_42*(p_affine_13_0*(tmp_0*tmp_70 + tmp_1*tmp_71 + tmp_4*tmp_72) + p_affine_13_1*(tmp_11*tmp_70 + tmp_2*tmp_72 + tmp_5*tmp_71) + p_affine_13_2*(tmp_10*tmp_70 + tmp_12*tmp_71 + tmp_72*tmp_8)) + 0.020848748529055869*tmp_42*(p_affine_13_0*(tmp_0*tmp_76 + tmp_1*tmp_77 + tmp_4*tmp_78) + p_affine_13_1*(tmp_11*tmp_76 + tmp_2*tmp_78 + tmp_5*tmp_77) + p_affine_13_2*(tmp_10*tmp_76 + tmp_12*tmp_77 + tmp_78*tmp_8)) + 0.019202922745021479*tmp_42*(p_affine_13_0*(tmp_0*tmp_82 + tmp_1*tmp_83 + tmp_4*tmp_84) + p_affine_13_1*(tmp_11*tmp_82 + tmp_2*tmp_84 + tmp_5*tmp_83) + p_affine_13_2*(tmp_10*tmp_82 + tmp_12*tmp_83 + tmp_8*tmp_84)) + 0.020848748529055869*tmp_42*(p_affine_13_0*(tmp_0*tmp_88 + tmp_1*tmp_89 + tmp_4*tmp_90) + p_affine_13_1*(tmp_11*tmp_88 + tmp_2*tmp_90 + tmp_5*tmp_89) + p_affine_13_2*(tmp_10*tmp_88 + tmp_12*tmp_89 + tmp_8*tmp_90)) + 0.019202922745021479*tmp_42*(p_affine_13_0*(tmp_0*tmp_94 + tmp_1*tmp_95 + tmp_4*tmp_96) + p_affine_13_1*(tmp_11*tmp_94 + tmp_2*tmp_96 + tmp_5*tmp_95) + p_affine_13_2*(tmp_10*tmp_94 + tmp_12*tmp_95 + tmp_8*tmp_96));
      elMat( 0, 0) = a_0_0;
   }

public:

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g1;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g1;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g0;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g0;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_3D_g2;

};


} //eg
} // dg
} // hyteg

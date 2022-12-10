
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

#include "hyteg/dgfunctionspace/DGBasisInfo.hpp"
#include "hyteg/dgfunctionspace/DGForm.hpp"
#include "hyteg/dgfunctionspace/DGForm2D.hpp"
#include "hyteg/types/matrix.hpp"
#include "hyteg/types/pointnd.hpp"

#include "Eigen/Eigen"

namespace hyteg {
namespace dg{    

class EGDivFormP0P1_new_0 : public hyteg::dg::DGForm2D
{
 public:
    EGDivFormP0P1_new_0()
: callback_Scalar_Variable_Coefficient_2D_g0 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_2D_g1 ([](const Point3D & p) -> real_t { return 0.; })
    {}

void Scalar_Variable_Coefficient_2D_g0( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g0( Point3D( {in_0, in_1, 0} ) );
}
void Scalar_Variable_Coefficient_2D_g1( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g1( Point3D( {in_0, in_1, 0} ) );
}

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

      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_0;
      real_t tmp_3 = 1.0 / (tmp_1*(p_affine_1_0 + tmp_2) - (p_affine_1_1 + tmp_0)*(p_affine_2_0 + tmp_2));
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = tmp_3*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_6 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_7 = tmp_6*(tmp_4 + tmp_5);
      real_t tmp_8 = tmp_4*tmp_6;
      real_t tmp_9 = tmp_5*tmp_6;
      real_t a_0_0 = 0.5*tmp_7;
      real_t a_0_1 = -0.5*tmp_8;
      real_t a_0_2 = -0.5*tmp_9;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = p_affine_6_1 + tmp_2;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.046910077030668018*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = 0.5*p_affine_10_0*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_18 = 0.11846344252809471*tmp_17;
      real_t tmp_19 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_20 = tmp_1*tmp_19;
      real_t tmp_21 = tmp_19*tmp_9;
      real_t tmp_22 = tmp_4*(0.23076534494715845*tmp_11 + tmp_12);
      real_t tmp_23 = tmp_22*tmp_3;
      real_t tmp_24 = tmp_15*tmp_22;
      real_t tmp_25 = 0.2393143352496831*tmp_17;
      real_t tmp_26 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_27 = tmp_1*tmp_26;
      real_t tmp_28 = tmp_26*tmp_9;
      real_t tmp_29 = tmp_4*(0.5*tmp_11 + tmp_12);
      real_t tmp_30 = tmp_29*tmp_3;
      real_t tmp_31 = tmp_15*tmp_29;
      real_t tmp_32 = 0.2844444444444445*tmp_17;
      real_t tmp_33 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_34 = tmp_1*tmp_33;
      real_t tmp_35 = tmp_33*tmp_9;
      real_t tmp_36 = tmp_4*(0.7692346550528415*tmp_11 + tmp_12);
      real_t tmp_37 = tmp_3*tmp_36;
      real_t tmp_38 = tmp_15*tmp_36;
      real_t tmp_39 = 0.2393143352496831*tmp_17;
      real_t tmp_40 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_41 = tmp_1*tmp_40;
      real_t tmp_42 = tmp_40*tmp_9;
      real_t tmp_43 = tmp_4*(0.95308992296933193*tmp_11 + tmp_12);
      real_t tmp_44 = tmp_3*tmp_43;
      real_t tmp_45 = tmp_15*tmp_43;
      real_t tmp_46 = 0.11846344252809471*tmp_17;
      real_t a_0_0 = tmp_18*(-tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1) + tmp_25*(-tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1) + tmp_32*(-tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1) + tmp_39*(-tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1) + tmp_46*(-tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1);
      real_t a_0_1 = tmp_18*(tmp_10 + tmp_14) + tmp_25*(tmp_21 + tmp_23) + tmp_32*(tmp_28 + tmp_30) + tmp_39*(tmp_35 + tmp_37) + tmp_46*(tmp_42 + tmp_44);
      real_t a_0_2 = tmp_18*(tmp_16 + tmp_8) + tmp_25*(tmp_20 + tmp_24) + tmp_32*(tmp_27 + tmp_31) + tmp_39*(tmp_34 + tmp_38) + tmp_46*(tmp_41 + tmp_45);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
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

      real_t tmp_0 = -p_affine_3_0;
      real_t tmp_1 = p_affine_4_0 + tmp_0;
      real_t tmp_2 = -p_affine_3_1;
      real_t tmp_3 = p_affine_5_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_4_1 + tmp_2)*(p_affine_5_0 + tmp_0));
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = p_affine_6_1 + tmp_2;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.046910077030668018*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = 0.5*p_affine_10_0*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_18 = 0.11846344252809471*tmp_17;
      real_t tmp_19 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_20 = tmp_1*tmp_19;
      real_t tmp_21 = tmp_19*tmp_9;
      real_t tmp_22 = tmp_4*(0.23076534494715845*tmp_11 + tmp_12);
      real_t tmp_23 = tmp_22*tmp_3;
      real_t tmp_24 = tmp_15*tmp_22;
      real_t tmp_25 = 0.2393143352496831*tmp_17;
      real_t tmp_26 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_27 = tmp_1*tmp_26;
      real_t tmp_28 = tmp_26*tmp_9;
      real_t tmp_29 = tmp_4*(0.5*tmp_11 + tmp_12);
      real_t tmp_30 = tmp_29*tmp_3;
      real_t tmp_31 = tmp_15*tmp_29;
      real_t tmp_32 = 0.2844444444444445*tmp_17;
      real_t tmp_33 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_34 = tmp_1*tmp_33;
      real_t tmp_35 = tmp_33*tmp_9;
      real_t tmp_36 = tmp_4*(0.7692346550528415*tmp_11 + tmp_12);
      real_t tmp_37 = tmp_3*tmp_36;
      real_t tmp_38 = tmp_15*tmp_36;
      real_t tmp_39 = 0.2393143352496831*tmp_17;
      real_t tmp_40 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_41 = tmp_1*tmp_40;
      real_t tmp_42 = tmp_40*tmp_9;
      real_t tmp_43 = tmp_4*(0.95308992296933193*tmp_11 + tmp_12);
      real_t tmp_44 = tmp_3*tmp_43;
      real_t tmp_45 = tmp_15*tmp_43;
      real_t tmp_46 = 0.11846344252809471*tmp_17;
      real_t a_0_0 = -tmp_18*(-tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1) - tmp_25*(-tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1) - tmp_32*(-tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1) - tmp_39*(-tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1) - tmp_46*(-tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1);
      real_t a_0_1 = -tmp_18*(tmp_10 + tmp_14) - tmp_25*(tmp_21 + tmp_23) - tmp_32*(tmp_28 + tmp_30) - tmp_39*(tmp_35 + tmp_37) - tmp_46*(tmp_42 + tmp_44);
      real_t a_0_2 = -tmp_18*(tmp_16 + tmp_8) - tmp_25*(tmp_20 + tmp_24) - tmp_32*(tmp_27 + tmp_31) - tmp_39*(tmp_34 + tmp_38) - tmp_46*(tmp_41 + tmp_45);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = p_affine_6_1 + tmp_2;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.046910077030668018*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = p_affine_10_0*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_18 = 0.11846344252809471*tmp_17;
      real_t tmp_19 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_20 = tmp_1*tmp_19;
      real_t tmp_21 = tmp_19*tmp_9;
      real_t tmp_22 = tmp_4*(0.23076534494715845*tmp_11 + tmp_12);
      real_t tmp_23 = tmp_22*tmp_3;
      real_t tmp_24 = tmp_15*tmp_22;
      real_t tmp_25 = 0.2393143352496831*tmp_17;
      real_t tmp_26 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_27 = tmp_1*tmp_26;
      real_t tmp_28 = tmp_26*tmp_9;
      real_t tmp_29 = tmp_4*(0.5*tmp_11 + tmp_12);
      real_t tmp_30 = tmp_29*tmp_3;
      real_t tmp_31 = tmp_15*tmp_29;
      real_t tmp_32 = 0.2844444444444445*tmp_17;
      real_t tmp_33 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_34 = tmp_1*tmp_33;
      real_t tmp_35 = tmp_33*tmp_9;
      real_t tmp_36 = tmp_4*(0.7692346550528415*tmp_11 + tmp_12);
      real_t tmp_37 = tmp_3*tmp_36;
      real_t tmp_38 = tmp_15*tmp_36;
      real_t tmp_39 = 0.2393143352496831*tmp_17;
      real_t tmp_40 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_41 = tmp_1*tmp_40;
      real_t tmp_42 = tmp_40*tmp_9;
      real_t tmp_43 = tmp_4*(0.95308992296933193*tmp_11 + tmp_12);
      real_t tmp_44 = tmp_3*tmp_43;
      real_t tmp_45 = tmp_15*tmp_43;
      real_t tmp_46 = 0.11846344252809471*tmp_17;
      real_t a_0_0 = tmp_18*(-tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1) + tmp_25*(-tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1) + tmp_32*(-tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1) + tmp_39*(-tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1) + tmp_46*(-tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1);
      real_t a_0_1 = tmp_18*(tmp_10 + tmp_14) + tmp_25*(tmp_21 + tmp_23) + tmp_32*(tmp_28 + tmp_30) + tmp_39*(tmp_35 + tmp_37) + tmp_46*(tmp_42 + tmp_44);
      real_t a_0_2 = tmp_18*(tmp_16 + tmp_8) + tmp_25*(tmp_20 + tmp_24) + tmp_32*(tmp_27 + tmp_31) + tmp_39*(tmp_34 + tmp_38) + tmp_46*(tmp_41 + tmp_45);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
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


      // Does nothing.
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

public:

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g0;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g1;
};




class EGDivFormP0P1_new_1 : public hyteg::dg::DGForm2D
{
 public:
    EGDivFormP0P1_new_1()
: callback_Scalar_Variable_Coefficient_2D_g0 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_2D_g1 ([](const Point3D & p) -> real_t { return 0.; })
    {}

void Scalar_Variable_Coefficient_2D_g0( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g0( Point3D( {in_0, in_1, 0} ) );
}
void Scalar_Variable_Coefficient_2D_g1( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g1( Point3D( {in_0, in_1, 0} ) );
}

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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = 1.0 / (tmp_1*(p_affine_2_1 + tmp_2) - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = tmp_3*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_6 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_7 = tmp_6*(tmp_4 + tmp_5);
      real_t tmp_8 = tmp_5*tmp_6;
      real_t tmp_9 = tmp_4*tmp_6;
      real_t a_0_0 = 0.5*tmp_7;
      real_t a_0_1 = -0.5*tmp_8;
      real_t a_0_2 = -0.5*tmp_9;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = p_affine_6_1 + tmp_2;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.046910077030668018*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = 0.5*p_affine_10_1*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_18 = 0.11846344252809471*tmp_17;
      real_t tmp_19 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_20 = tmp_1*tmp_19;
      real_t tmp_21 = tmp_19*tmp_9;
      real_t tmp_22 = tmp_4*(0.23076534494715845*tmp_11 + tmp_12);
      real_t tmp_23 = tmp_22*tmp_3;
      real_t tmp_24 = tmp_15*tmp_22;
      real_t tmp_25 = 0.2393143352496831*tmp_17;
      real_t tmp_26 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_27 = tmp_1*tmp_26;
      real_t tmp_28 = tmp_26*tmp_9;
      real_t tmp_29 = tmp_4*(0.5*tmp_11 + tmp_12);
      real_t tmp_30 = tmp_29*tmp_3;
      real_t tmp_31 = tmp_15*tmp_29;
      real_t tmp_32 = 0.2844444444444445*tmp_17;
      real_t tmp_33 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_34 = tmp_1*tmp_33;
      real_t tmp_35 = tmp_33*tmp_9;
      real_t tmp_36 = tmp_4*(0.7692346550528415*tmp_11 + tmp_12);
      real_t tmp_37 = tmp_3*tmp_36;
      real_t tmp_38 = tmp_15*tmp_36;
      real_t tmp_39 = 0.2393143352496831*tmp_17;
      real_t tmp_40 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_41 = tmp_1*tmp_40;
      real_t tmp_42 = tmp_40*tmp_9;
      real_t tmp_43 = tmp_4*(0.95308992296933193*tmp_11 + tmp_12);
      real_t tmp_44 = tmp_3*tmp_43;
      real_t tmp_45 = tmp_15*tmp_43;
      real_t tmp_46 = 0.11846344252809471*tmp_17;
      real_t a_0_0 = tmp_18*(-tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1) + tmp_25*(-tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1) + tmp_32*(-tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1) + tmp_39*(-tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1) + tmp_46*(-tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1);
      real_t a_0_1 = tmp_18*(tmp_10 + tmp_14) + tmp_25*(tmp_21 + tmp_23) + tmp_32*(tmp_28 + tmp_30) + tmp_39*(tmp_35 + tmp_37) + tmp_46*(tmp_42 + tmp_44);
      real_t a_0_2 = tmp_18*(tmp_16 + tmp_8) + tmp_25*(tmp_20 + tmp_24) + tmp_32*(tmp_27 + tmp_31) + tmp_39*(tmp_34 + tmp_38) + tmp_46*(tmp_41 + tmp_45);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
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

      real_t tmp_0 = -p_affine_3_0;
      real_t tmp_1 = p_affine_4_0 + tmp_0;
      real_t tmp_2 = -p_affine_3_1;
      real_t tmp_3 = p_affine_5_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_4_1 + tmp_2)*(p_affine_5_0 + tmp_0));
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = p_affine_6_1 + tmp_2;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.046910077030668018*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = 0.5*p_affine_10_1*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_18 = 0.11846344252809471*tmp_17;
      real_t tmp_19 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_20 = tmp_1*tmp_19;
      real_t tmp_21 = tmp_19*tmp_9;
      real_t tmp_22 = tmp_4*(0.23076534494715845*tmp_11 + tmp_12);
      real_t tmp_23 = tmp_22*tmp_3;
      real_t tmp_24 = tmp_15*tmp_22;
      real_t tmp_25 = 0.2393143352496831*tmp_17;
      real_t tmp_26 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_27 = tmp_1*tmp_26;
      real_t tmp_28 = tmp_26*tmp_9;
      real_t tmp_29 = tmp_4*(0.5*tmp_11 + tmp_12);
      real_t tmp_30 = tmp_29*tmp_3;
      real_t tmp_31 = tmp_15*tmp_29;
      real_t tmp_32 = 0.2844444444444445*tmp_17;
      real_t tmp_33 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_34 = tmp_1*tmp_33;
      real_t tmp_35 = tmp_33*tmp_9;
      real_t tmp_36 = tmp_4*(0.7692346550528415*tmp_11 + tmp_12);
      real_t tmp_37 = tmp_3*tmp_36;
      real_t tmp_38 = tmp_15*tmp_36;
      real_t tmp_39 = 0.2393143352496831*tmp_17;
      real_t tmp_40 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_41 = tmp_1*tmp_40;
      real_t tmp_42 = tmp_40*tmp_9;
      real_t tmp_43 = tmp_4*(0.95308992296933193*tmp_11 + tmp_12);
      real_t tmp_44 = tmp_3*tmp_43;
      real_t tmp_45 = tmp_15*tmp_43;
      real_t tmp_46 = 0.11846344252809471*tmp_17;
      real_t a_0_0 = -tmp_18*(-tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1) - tmp_25*(-tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1) - tmp_32*(-tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1) - tmp_39*(-tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1) - tmp_46*(-tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1);
      real_t a_0_1 = -tmp_18*(tmp_10 + tmp_14) - tmp_25*(tmp_21 + tmp_23) - tmp_32*(tmp_28 + tmp_30) - tmp_39*(tmp_35 + tmp_37) - tmp_46*(tmp_42 + tmp_44);
      real_t a_0_2 = -tmp_18*(tmp_16 + tmp_8) - tmp_25*(tmp_20 + tmp_24) - tmp_32*(tmp_27 + tmp_31) - tmp_39*(tmp_34 + tmp_38) - tmp_46*(tmp_41 + tmp_45);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = p_affine_6_1 + tmp_2;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.046910077030668018*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = p_affine_10_1*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_18 = 0.11846344252809471*tmp_17;
      real_t tmp_19 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_20 = tmp_1*tmp_19;
      real_t tmp_21 = tmp_19*tmp_9;
      real_t tmp_22 = tmp_4*(0.23076534494715845*tmp_11 + tmp_12);
      real_t tmp_23 = tmp_22*tmp_3;
      real_t tmp_24 = tmp_15*tmp_22;
      real_t tmp_25 = 0.2393143352496831*tmp_17;
      real_t tmp_26 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_27 = tmp_1*tmp_26;
      real_t tmp_28 = tmp_26*tmp_9;
      real_t tmp_29 = tmp_4*(0.5*tmp_11 + tmp_12);
      real_t tmp_30 = tmp_29*tmp_3;
      real_t tmp_31 = tmp_15*tmp_29;
      real_t tmp_32 = 0.2844444444444445*tmp_17;
      real_t tmp_33 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_34 = tmp_1*tmp_33;
      real_t tmp_35 = tmp_33*tmp_9;
      real_t tmp_36 = tmp_4*(0.7692346550528415*tmp_11 + tmp_12);
      real_t tmp_37 = tmp_3*tmp_36;
      real_t tmp_38 = tmp_15*tmp_36;
      real_t tmp_39 = 0.2393143352496831*tmp_17;
      real_t tmp_40 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_41 = tmp_1*tmp_40;
      real_t tmp_42 = tmp_40*tmp_9;
      real_t tmp_43 = tmp_4*(0.95308992296933193*tmp_11 + tmp_12);
      real_t tmp_44 = tmp_3*tmp_43;
      real_t tmp_45 = tmp_15*tmp_43;
      real_t tmp_46 = 0.11846344252809471*tmp_17;
      real_t a_0_0 = tmp_18*(-tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1) + tmp_25*(-tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1) + tmp_32*(-tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1) + tmp_39*(-tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1) + tmp_46*(-tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1);
      real_t a_0_1 = tmp_18*(tmp_10 + tmp_14) + tmp_25*(tmp_21 + tmp_23) + tmp_32*(tmp_28 + tmp_30) + tmp_39*(tmp_35 + tmp_37) + tmp_46*(tmp_42 + tmp_44);
      real_t a_0_2 = tmp_18*(tmp_16 + tmp_8) + tmp_25*(tmp_20 + tmp_24) + tmp_32*(tmp_27 + tmp_31) + tmp_39*(tmp_34 + tmp_38) + tmp_46*(tmp_41 + tmp_45);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
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


      // Does nothing.
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

public:

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g0;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g1;
};




class EGDivFormP0EDG_new : public hyteg::dg::DGForm2D
{
 public:
    EGDivFormP0EDG_new()
: callback_Scalar_Variable_Coefficient_2D_g0 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_2D_g1 ([](const Point3D & p) -> real_t { return 0.; })
    {}

void Scalar_Variable_Coefficient_2D_g0( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g0( Point3D( {in_0, in_1, 0} ) );
}
void Scalar_Variable_Coefficient_2D_g1( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g1( Point3D( {in_0, in_1, 0} ) );
}

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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = -p_affine_0_1;
      real_t tmp_2 = (p_affine_1_0 + tmp_0)*(p_affine_2_1 + tmp_1);
      real_t tmp_3 = p_affine_2_0 + tmp_0;
      real_t tmp_4 = p_affine_1_1 + tmp_1;
      real_t tmp_5 = 1.0 / (tmp_2 - tmp_3*tmp_4);
      real_t tmp_6 = (-2*tmp_2*tmp_5 - tmp_3*tmp_5*(p_affine_0_1 - p_affine_1_1) - tmp_4*tmp_5*(p_affine_0_0 - p_affine_2_0))*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = 0.5*tmp_6;
      elMat( 0, 0) = a_0_0;
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

      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = p_affine_1_0 + tmp_3;
      real_t tmp_5 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_6 = -p_affine_0_1;
      real_t tmp_7 = p_affine_2_1 + tmp_6;
      real_t tmp_8 = p_affine_2_0 + tmp_3;
      real_t tmp_9 = p_affine_1_1 + tmp_6;
      real_t tmp_10 = 1.0 / (tmp_4*tmp_7 - tmp_8*tmp_9);
      real_t tmp_11 = p_affine_6_1 + tmp_6;
      real_t tmp_12 = tmp_10*(0.046910077030668018*tmp_1 + tmp_11);
      real_t tmp_13 = p_affine_6_0 + tmp_3;
      real_t tmp_14 = tmp_10*(0.046910077030668018*tmp_0 + tmp_13);
      real_t tmp_15 = tmp_12*tmp_5 + tmp_14*tmp_7 - 1.0/3.0;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_12*tmp_4 + tmp_14*tmp_16 - 1.0/3.0;
      real_t tmp_18 = 0.5*p_affine_10_0;
      real_t tmp_19 = 0.5*p_affine_10_1;
      real_t tmp_20 = tmp_10*(0.23076534494715845*tmp_1 + tmp_11);
      real_t tmp_21 = tmp_10*(0.23076534494715845*tmp_0 + tmp_13);
      real_t tmp_22 = tmp_20*tmp_5 + tmp_21*tmp_7 - 1.0/3.0;
      real_t tmp_23 = tmp_16*tmp_21 + tmp_20*tmp_4 - 1.0/3.0;
      real_t tmp_24 = tmp_10*(0.5*tmp_1 + tmp_11);
      real_t tmp_25 = tmp_10*(0.5*tmp_0 + tmp_13);
      real_t tmp_26 = tmp_24*tmp_5 + tmp_25*tmp_7 - 1.0/3.0;
      real_t tmp_27 = tmp_16*tmp_25 + tmp_24*tmp_4 - 1.0/3.0;
      real_t tmp_28 = tmp_10*(0.7692346550528415*tmp_1 + tmp_11);
      real_t tmp_29 = tmp_10*(0.7692346550528415*tmp_0 + tmp_13);
      real_t tmp_30 = tmp_28*tmp_5 + tmp_29*tmp_7 - 1.0/3.0;
      real_t tmp_31 = tmp_16*tmp_29 + tmp_28*tmp_4 - 1.0/3.0;
      real_t tmp_32 = tmp_10*(0.95308992296933193*tmp_1 + tmp_11);
      real_t tmp_33 = tmp_10*(0.95308992296933193*tmp_0 + tmp_13);
      real_t tmp_34 = tmp_32*tmp_5 + tmp_33*tmp_7 - 1.0/3.0;
      real_t tmp_35 = tmp_16*tmp_33 + tmp_32*tmp_4 - 1.0/3.0;
      real_t a_0_0 = 0.11846344252809471*tmp_2*(tmp_18*(tmp_15*tmp_4 + tmp_17*tmp_8) + tmp_19*(tmp_15*tmp_9 + tmp_17*tmp_7)) + 0.2393143352496831*tmp_2*(tmp_18*(tmp_22*tmp_4 + tmp_23*tmp_8) + tmp_19*(tmp_22*tmp_9 + tmp_23*tmp_7)) + 0.2844444444444445*tmp_2*(tmp_18*(tmp_26*tmp_4 + tmp_27*tmp_8) + tmp_19*(tmp_26*tmp_9 + tmp_27*tmp_7)) + 0.2393143352496831*tmp_2*(tmp_18*(tmp_30*tmp_4 + tmp_31*tmp_8) + tmp_19*(tmp_30*tmp_9 + tmp_31*tmp_7)) + 0.11846344252809471*tmp_2*(tmp_18*(tmp_34*tmp_4 + tmp_35*tmp_8) + tmp_19*(tmp_34*tmp_9 + tmp_35*tmp_7));
      elMat( 0, 0) = a_0_0;
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

      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = -p_affine_3_0;
      real_t tmp_4 = p_affine_4_0 + tmp_3;
      real_t tmp_5 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_6 = -p_affine_3_1;
      real_t tmp_7 = p_affine_5_1 + tmp_6;
      real_t tmp_8 = p_affine_5_0 + tmp_3;
      real_t tmp_9 = p_affine_4_1 + tmp_6;
      real_t tmp_10 = 1.0 / (tmp_4*tmp_7 - tmp_8*tmp_9);
      real_t tmp_11 = p_affine_6_1 + tmp_6;
      real_t tmp_12 = tmp_10*(0.046910077030668018*tmp_1 + tmp_11);
      real_t tmp_13 = p_affine_6_0 + tmp_3;
      real_t tmp_14 = tmp_10*(0.046910077030668018*tmp_0 + tmp_13);
      real_t tmp_15 = tmp_12*tmp_5 + tmp_14*tmp_7 - 1.0/3.0;
      real_t tmp_16 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_17 = tmp_12*tmp_4 + tmp_14*tmp_16 - 1.0/3.0;
      real_t tmp_18 = 0.5*p_affine_10_0;
      real_t tmp_19 = 0.5*p_affine_10_1;
      real_t tmp_20 = tmp_10*(0.23076534494715845*tmp_1 + tmp_11);
      real_t tmp_21 = tmp_10*(0.23076534494715845*tmp_0 + tmp_13);
      real_t tmp_22 = tmp_20*tmp_5 + tmp_21*tmp_7 - 1.0/3.0;
      real_t tmp_23 = tmp_16*tmp_21 + tmp_20*tmp_4 - 1.0/3.0;
      real_t tmp_24 = tmp_10*(0.5*tmp_1 + tmp_11);
      real_t tmp_25 = tmp_10*(0.5*tmp_0 + tmp_13);
      real_t tmp_26 = tmp_24*tmp_5 + tmp_25*tmp_7 - 1.0/3.0;
      real_t tmp_27 = tmp_16*tmp_25 + tmp_24*tmp_4 - 1.0/3.0;
      real_t tmp_28 = tmp_10*(0.7692346550528415*tmp_1 + tmp_11);
      real_t tmp_29 = tmp_10*(0.7692346550528415*tmp_0 + tmp_13);
      real_t tmp_30 = tmp_28*tmp_5 + tmp_29*tmp_7 - 1.0/3.0;
      real_t tmp_31 = tmp_16*tmp_29 + tmp_28*tmp_4 - 1.0/3.0;
      real_t tmp_32 = tmp_10*(0.95308992296933193*tmp_1 + tmp_11);
      real_t tmp_33 = tmp_10*(0.95308992296933193*tmp_0 + tmp_13);
      real_t tmp_34 = tmp_32*tmp_5 + tmp_33*tmp_7 - 1.0/3.0;
      real_t tmp_35 = tmp_16*tmp_33 + tmp_32*tmp_4 - 1.0/3.0;
      real_t a_0_0 = 0.11846344252809471*tmp_2*(-tmp_18*(tmp_15*tmp_4 + tmp_17*tmp_8) - tmp_19*(tmp_15*tmp_9 + tmp_17*tmp_7)) + 0.2393143352496831*tmp_2*(-tmp_18*(tmp_22*tmp_4 + tmp_23*tmp_8) - tmp_19*(tmp_22*tmp_9 + tmp_23*tmp_7)) + 0.2844444444444445*tmp_2*(-tmp_18*(tmp_26*tmp_4 + tmp_27*tmp_8) - tmp_19*(tmp_26*tmp_9 + tmp_27*tmp_7)) + 0.2393143352496831*tmp_2*(-tmp_18*(tmp_30*tmp_4 + tmp_31*tmp_8) - tmp_19*(tmp_30*tmp_9 + tmp_31*tmp_7)) + 0.11846344252809471*tmp_2*(-tmp_18*(tmp_34*tmp_4 + tmp_35*tmp_8) - tmp_19*(tmp_34*tmp_9 + tmp_35*tmp_7));
      elMat( 0, 0) = a_0_0;
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

      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = p_affine_1_0 + tmp_3;
      real_t tmp_5 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_6 = -p_affine_0_1;
      real_t tmp_7 = p_affine_2_1 + tmp_6;
      real_t tmp_8 = p_affine_2_0 + tmp_3;
      real_t tmp_9 = p_affine_1_1 + tmp_6;
      real_t tmp_10 = 1.0 / (tmp_4*tmp_7 - tmp_8*tmp_9);
      real_t tmp_11 = p_affine_6_1 + tmp_6;
      real_t tmp_12 = tmp_10*(0.046910077030668018*tmp_1 + tmp_11);
      real_t tmp_13 = p_affine_6_0 + tmp_3;
      real_t tmp_14 = tmp_10*(0.046910077030668018*tmp_0 + tmp_13);
      real_t tmp_15 = tmp_12*tmp_5 + tmp_14*tmp_7 - 1.0/3.0;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_12*tmp_4 + tmp_14*tmp_16 - 1.0/3.0;
      real_t tmp_18 = tmp_10*(0.23076534494715845*tmp_1 + tmp_11);
      real_t tmp_19 = tmp_10*(0.23076534494715845*tmp_0 + tmp_13);
      real_t tmp_20 = tmp_18*tmp_5 + tmp_19*tmp_7 - 1.0/3.0;
      real_t tmp_21 = tmp_16*tmp_19 + tmp_18*tmp_4 - 1.0/3.0;
      real_t tmp_22 = tmp_10*(0.5*tmp_1 + tmp_11);
      real_t tmp_23 = tmp_10*(0.5*tmp_0 + tmp_13);
      real_t tmp_24 = tmp_22*tmp_5 + tmp_23*tmp_7 - 1.0/3.0;
      real_t tmp_25 = tmp_16*tmp_23 + tmp_22*tmp_4 - 1.0/3.0;
      real_t tmp_26 = tmp_10*(0.7692346550528415*tmp_1 + tmp_11);
      real_t tmp_27 = tmp_10*(0.7692346550528415*tmp_0 + tmp_13);
      real_t tmp_28 = tmp_26*tmp_5 + tmp_27*tmp_7 - 1.0/3.0;
      real_t tmp_29 = tmp_16*tmp_27 + tmp_26*tmp_4 - 1.0/3.0;
      real_t tmp_30 = tmp_10*(0.95308992296933193*tmp_1 + tmp_11);
      real_t tmp_31 = tmp_10*(0.95308992296933193*tmp_0 + tmp_13);
      real_t tmp_32 = tmp_30*tmp_5 + tmp_31*tmp_7 - 1.0/3.0;
      real_t tmp_33 = tmp_16*tmp_31 + tmp_30*tmp_4 - 1.0/3.0;
      real_t a_0_0 = 0.11846344252809471*tmp_2*(p_affine_10_0*(tmp_15*tmp_4 + tmp_17*tmp_8) + p_affine_10_1*(tmp_15*tmp_9 + tmp_17*tmp_7)) + 0.2393143352496831*tmp_2*(p_affine_10_0*(tmp_20*tmp_4 + tmp_21*tmp_8) + p_affine_10_1*(tmp_20*tmp_9 + tmp_21*tmp_7)) + 0.2844444444444445*tmp_2*(p_affine_10_0*(tmp_24*tmp_4 + tmp_25*tmp_8) + p_affine_10_1*(tmp_24*tmp_9 + tmp_25*tmp_7)) + 0.2393143352496831*tmp_2*(p_affine_10_0*(tmp_28*tmp_4 + tmp_29*tmp_8) + p_affine_10_1*(tmp_28*tmp_9 + tmp_29*tmp_7)) + 0.11846344252809471*tmp_2*(p_affine_10_0*(tmp_32*tmp_4 + tmp_33*tmp_8) + p_affine_10_1*(tmp_32*tmp_9 + tmp_33*tmp_7));
      elMat( 0, 0) = a_0_0;
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


      // Does nothing.
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

public:

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g0;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g1;
};




class EGDivtFormP1P0_new_0 : public hyteg::dg::DGForm2D
{
 public:
    EGDivtFormP1P0_new_0()

    {}



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

      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_0;
      real_t tmp_3 = 1.0 / (tmp_1*(p_affine_1_0 + tmp_2) - (p_affine_1_1 + tmp_0)*(p_affine_2_0 + tmp_2));
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = tmp_3*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_6 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_7 = tmp_6*(tmp_4 + tmp_5);
      real_t tmp_8 = tmp_4*tmp_6;
      real_t tmp_9 = tmp_5*tmp_6;
      real_t a_0_0 = 0.5*tmp_7;
      real_t a_1_0 = -0.5*tmp_8;
      real_t a_2_0 = -0.5*tmp_9;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = p_affine_6_1 + tmp_2;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.046910077030668018*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = 0.5*p_affine_10_0*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_18 = 0.11846344252809471*tmp_17;
      real_t tmp_19 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_20 = tmp_1*tmp_19;
      real_t tmp_21 = tmp_19*tmp_9;
      real_t tmp_22 = tmp_4*(0.23076534494715845*tmp_11 + tmp_12);
      real_t tmp_23 = tmp_22*tmp_3;
      real_t tmp_24 = tmp_15*tmp_22;
      real_t tmp_25 = 0.2393143352496831*tmp_17;
      real_t tmp_26 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_27 = tmp_1*tmp_26;
      real_t tmp_28 = tmp_26*tmp_9;
      real_t tmp_29 = tmp_4*(0.5*tmp_11 + tmp_12);
      real_t tmp_30 = tmp_29*tmp_3;
      real_t tmp_31 = tmp_15*tmp_29;
      real_t tmp_32 = 0.2844444444444445*tmp_17;
      real_t tmp_33 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_34 = tmp_1*tmp_33;
      real_t tmp_35 = tmp_33*tmp_9;
      real_t tmp_36 = tmp_4*(0.7692346550528415*tmp_11 + tmp_12);
      real_t tmp_37 = tmp_3*tmp_36;
      real_t tmp_38 = tmp_15*tmp_36;
      real_t tmp_39 = 0.2393143352496831*tmp_17;
      real_t tmp_40 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_41 = tmp_1*tmp_40;
      real_t tmp_42 = tmp_40*tmp_9;
      real_t tmp_43 = tmp_4*(0.95308992296933193*tmp_11 + tmp_12);
      real_t tmp_44 = tmp_3*tmp_43;
      real_t tmp_45 = tmp_15*tmp_43;
      real_t tmp_46 = 0.11846344252809471*tmp_17;
      real_t a_0_0 = tmp_18*(-tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1) + tmp_25*(-tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1) + tmp_32*(-tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1) + tmp_39*(-tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1) + tmp_46*(-tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1);
      real_t a_1_0 = tmp_18*(tmp_10 + tmp_14) + tmp_25*(tmp_21 + tmp_23) + tmp_32*(tmp_28 + tmp_30) + tmp_39*(tmp_35 + tmp_37) + tmp_46*(tmp_42 + tmp_44);
      real_t a_2_0 = tmp_18*(tmp_16 + tmp_8) + tmp_25*(tmp_20 + tmp_24) + tmp_32*(tmp_27 + tmp_31) + tmp_39*(tmp_34 + tmp_38) + tmp_46*(tmp_41 + tmp_45);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = p_affine_6_1 + tmp_2;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.046910077030668018*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = 0.5*p_affine_10_0*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_18 = 0.11846344252809471*tmp_17;
      real_t tmp_19 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_20 = tmp_1*tmp_19;
      real_t tmp_21 = tmp_19*tmp_9;
      real_t tmp_22 = tmp_4*(0.23076534494715845*tmp_11 + tmp_12);
      real_t tmp_23 = tmp_22*tmp_3;
      real_t tmp_24 = tmp_15*tmp_22;
      real_t tmp_25 = 0.2393143352496831*tmp_17;
      real_t tmp_26 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_27 = tmp_1*tmp_26;
      real_t tmp_28 = tmp_26*tmp_9;
      real_t tmp_29 = tmp_4*(0.5*tmp_11 + tmp_12);
      real_t tmp_30 = tmp_29*tmp_3;
      real_t tmp_31 = tmp_15*tmp_29;
      real_t tmp_32 = 0.2844444444444445*tmp_17;
      real_t tmp_33 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_34 = tmp_1*tmp_33;
      real_t tmp_35 = tmp_33*tmp_9;
      real_t tmp_36 = tmp_4*(0.7692346550528415*tmp_11 + tmp_12);
      real_t tmp_37 = tmp_3*tmp_36;
      real_t tmp_38 = tmp_15*tmp_36;
      real_t tmp_39 = 0.2393143352496831*tmp_17;
      real_t tmp_40 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_41 = tmp_1*tmp_40;
      real_t tmp_42 = tmp_40*tmp_9;
      real_t tmp_43 = tmp_4*(0.95308992296933193*tmp_11 + tmp_12);
      real_t tmp_44 = tmp_3*tmp_43;
      real_t tmp_45 = tmp_15*tmp_43;
      real_t tmp_46 = 0.11846344252809471*tmp_17;
      real_t a_0_0 = tmp_18*(-tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1) + tmp_25*(-tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1) + tmp_32*(-tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1) + tmp_39*(-tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1) + tmp_46*(-tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1);
      real_t a_1_0 = tmp_18*(tmp_10 + tmp_14) + tmp_25*(tmp_21 + tmp_23) + tmp_32*(tmp_28 + tmp_30) + tmp_39*(tmp_35 + tmp_37) + tmp_46*(tmp_42 + tmp_44);
      real_t a_2_0 = tmp_18*(tmp_16 + tmp_8) + tmp_25*(tmp_20 + tmp_24) + tmp_32*(tmp_27 + tmp_31) + tmp_39*(tmp_34 + tmp_38) + tmp_46*(tmp_41 + tmp_45);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = p_affine_6_1 + tmp_2;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.046910077030668018*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = p_affine_10_0*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_18 = 0.11846344252809471*tmp_17;
      real_t tmp_19 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_20 = tmp_1*tmp_19;
      real_t tmp_21 = tmp_19*tmp_9;
      real_t tmp_22 = tmp_4*(0.23076534494715845*tmp_11 + tmp_12);
      real_t tmp_23 = tmp_22*tmp_3;
      real_t tmp_24 = tmp_15*tmp_22;
      real_t tmp_25 = 0.2393143352496831*tmp_17;
      real_t tmp_26 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_27 = tmp_1*tmp_26;
      real_t tmp_28 = tmp_26*tmp_9;
      real_t tmp_29 = tmp_4*(0.5*tmp_11 + tmp_12);
      real_t tmp_30 = tmp_29*tmp_3;
      real_t tmp_31 = tmp_15*tmp_29;
      real_t tmp_32 = 0.2844444444444445*tmp_17;
      real_t tmp_33 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_34 = tmp_1*tmp_33;
      real_t tmp_35 = tmp_33*tmp_9;
      real_t tmp_36 = tmp_4*(0.7692346550528415*tmp_11 + tmp_12);
      real_t tmp_37 = tmp_3*tmp_36;
      real_t tmp_38 = tmp_15*tmp_36;
      real_t tmp_39 = 0.2393143352496831*tmp_17;
      real_t tmp_40 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_41 = tmp_1*tmp_40;
      real_t tmp_42 = tmp_40*tmp_9;
      real_t tmp_43 = tmp_4*(0.95308992296933193*tmp_11 + tmp_12);
      real_t tmp_44 = tmp_3*tmp_43;
      real_t tmp_45 = tmp_15*tmp_43;
      real_t tmp_46 = 0.11846344252809471*tmp_17;
      real_t a_0_0 = tmp_18*(-tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1) + tmp_25*(-tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1) + tmp_32*(-tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1) + tmp_39*(-tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1) + tmp_46*(-tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1);
      real_t a_1_0 = tmp_18*(tmp_10 + tmp_14) + tmp_25*(tmp_21 + tmp_23) + tmp_32*(tmp_28 + tmp_30) + tmp_39*(tmp_35 + tmp_37) + tmp_46*(tmp_42 + tmp_44);
      real_t a_2_0 = tmp_18*(tmp_16 + tmp_8) + tmp_25*(tmp_20 + tmp_24) + tmp_32*(tmp_27 + tmp_31) + tmp_39*(tmp_34 + tmp_38) + tmp_46*(tmp_41 + tmp_45);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
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


      // Does nothing.
      elMat( 0, 0) = 0;
      elMat( 1, 0) = 0;
      elMat( 2, 0) = 0;
   }

public:


};




class EGDivtFormP1P0_new_1 : public hyteg::dg::DGForm2D
{
 public:
    EGDivtFormP1P0_new_1()

    {}



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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = 1.0 / (tmp_1*(p_affine_2_1 + tmp_2) - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = tmp_3*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_6 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_7 = tmp_6*(tmp_4 + tmp_5);
      real_t tmp_8 = tmp_5*tmp_6;
      real_t tmp_9 = tmp_4*tmp_6;
      real_t a_0_0 = 0.5*tmp_7;
      real_t a_1_0 = -0.5*tmp_8;
      real_t a_2_0 = -0.5*tmp_9;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = p_affine_6_1 + tmp_2;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.046910077030668018*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = 0.5*p_affine_10_1*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_18 = 0.11846344252809471*tmp_17;
      real_t tmp_19 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_20 = tmp_1*tmp_19;
      real_t tmp_21 = tmp_19*tmp_9;
      real_t tmp_22 = tmp_4*(0.23076534494715845*tmp_11 + tmp_12);
      real_t tmp_23 = tmp_22*tmp_3;
      real_t tmp_24 = tmp_15*tmp_22;
      real_t tmp_25 = 0.2393143352496831*tmp_17;
      real_t tmp_26 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_27 = tmp_1*tmp_26;
      real_t tmp_28 = tmp_26*tmp_9;
      real_t tmp_29 = tmp_4*(0.5*tmp_11 + tmp_12);
      real_t tmp_30 = tmp_29*tmp_3;
      real_t tmp_31 = tmp_15*tmp_29;
      real_t tmp_32 = 0.2844444444444445*tmp_17;
      real_t tmp_33 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_34 = tmp_1*tmp_33;
      real_t tmp_35 = tmp_33*tmp_9;
      real_t tmp_36 = tmp_4*(0.7692346550528415*tmp_11 + tmp_12);
      real_t tmp_37 = tmp_3*tmp_36;
      real_t tmp_38 = tmp_15*tmp_36;
      real_t tmp_39 = 0.2393143352496831*tmp_17;
      real_t tmp_40 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_41 = tmp_1*tmp_40;
      real_t tmp_42 = tmp_40*tmp_9;
      real_t tmp_43 = tmp_4*(0.95308992296933193*tmp_11 + tmp_12);
      real_t tmp_44 = tmp_3*tmp_43;
      real_t tmp_45 = tmp_15*tmp_43;
      real_t tmp_46 = 0.11846344252809471*tmp_17;
      real_t a_0_0 = tmp_18*(-tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1) + tmp_25*(-tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1) + tmp_32*(-tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1) + tmp_39*(-tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1) + tmp_46*(-tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1);
      real_t a_1_0 = tmp_18*(tmp_10 + tmp_14) + tmp_25*(tmp_21 + tmp_23) + tmp_32*(tmp_28 + tmp_30) + tmp_39*(tmp_35 + tmp_37) + tmp_46*(tmp_42 + tmp_44);
      real_t a_2_0 = tmp_18*(tmp_16 + tmp_8) + tmp_25*(tmp_20 + tmp_24) + tmp_32*(tmp_27 + tmp_31) + tmp_39*(tmp_34 + tmp_38) + tmp_46*(tmp_41 + tmp_45);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = p_affine_6_1 + tmp_2;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.046910077030668018*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = 0.5*p_affine_10_1*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_18 = 0.11846344252809471*tmp_17;
      real_t tmp_19 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_20 = tmp_1*tmp_19;
      real_t tmp_21 = tmp_19*tmp_9;
      real_t tmp_22 = tmp_4*(0.23076534494715845*tmp_11 + tmp_12);
      real_t tmp_23 = tmp_22*tmp_3;
      real_t tmp_24 = tmp_15*tmp_22;
      real_t tmp_25 = 0.2393143352496831*tmp_17;
      real_t tmp_26 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_27 = tmp_1*tmp_26;
      real_t tmp_28 = tmp_26*tmp_9;
      real_t tmp_29 = tmp_4*(0.5*tmp_11 + tmp_12);
      real_t tmp_30 = tmp_29*tmp_3;
      real_t tmp_31 = tmp_15*tmp_29;
      real_t tmp_32 = 0.2844444444444445*tmp_17;
      real_t tmp_33 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_34 = tmp_1*tmp_33;
      real_t tmp_35 = tmp_33*tmp_9;
      real_t tmp_36 = tmp_4*(0.7692346550528415*tmp_11 + tmp_12);
      real_t tmp_37 = tmp_3*tmp_36;
      real_t tmp_38 = tmp_15*tmp_36;
      real_t tmp_39 = 0.2393143352496831*tmp_17;
      real_t tmp_40 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_41 = tmp_1*tmp_40;
      real_t tmp_42 = tmp_40*tmp_9;
      real_t tmp_43 = tmp_4*(0.95308992296933193*tmp_11 + tmp_12);
      real_t tmp_44 = tmp_3*tmp_43;
      real_t tmp_45 = tmp_15*tmp_43;
      real_t tmp_46 = 0.11846344252809471*tmp_17;
      real_t a_0_0 = tmp_18*(-tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1) + tmp_25*(-tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1) + tmp_32*(-tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1) + tmp_39*(-tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1) + tmp_46*(-tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1);
      real_t a_1_0 = tmp_18*(tmp_10 + tmp_14) + tmp_25*(tmp_21 + tmp_23) + tmp_32*(tmp_28 + tmp_30) + tmp_39*(tmp_35 + tmp_37) + tmp_46*(tmp_42 + tmp_44);
      real_t a_2_0 = tmp_18*(tmp_16 + tmp_8) + tmp_25*(tmp_20 + tmp_24) + tmp_32*(tmp_27 + tmp_31) + tmp_39*(tmp_34 + tmp_38) + tmp_46*(tmp_41 + tmp_45);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6 = p_affine_6_1 + tmp_2;
      real_t tmp_7 = tmp_4*(0.046910077030668018*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.046910077030668018*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = p_affine_10_1*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_18 = 0.11846344252809471*tmp_17;
      real_t tmp_19 = tmp_4*(0.23076534494715845*tmp_5 + tmp_6);
      real_t tmp_20 = tmp_1*tmp_19;
      real_t tmp_21 = tmp_19*tmp_9;
      real_t tmp_22 = tmp_4*(0.23076534494715845*tmp_11 + tmp_12);
      real_t tmp_23 = tmp_22*tmp_3;
      real_t tmp_24 = tmp_15*tmp_22;
      real_t tmp_25 = 0.2393143352496831*tmp_17;
      real_t tmp_26 = tmp_4*(0.5*tmp_5 + tmp_6);
      real_t tmp_27 = tmp_1*tmp_26;
      real_t tmp_28 = tmp_26*tmp_9;
      real_t tmp_29 = tmp_4*(0.5*tmp_11 + tmp_12);
      real_t tmp_30 = tmp_29*tmp_3;
      real_t tmp_31 = tmp_15*tmp_29;
      real_t tmp_32 = 0.2844444444444445*tmp_17;
      real_t tmp_33 = tmp_4*(0.7692346550528415*tmp_5 + tmp_6);
      real_t tmp_34 = tmp_1*tmp_33;
      real_t tmp_35 = tmp_33*tmp_9;
      real_t tmp_36 = tmp_4*(0.7692346550528415*tmp_11 + tmp_12);
      real_t tmp_37 = tmp_3*tmp_36;
      real_t tmp_38 = tmp_15*tmp_36;
      real_t tmp_39 = 0.2393143352496831*tmp_17;
      real_t tmp_40 = tmp_4*(0.95308992296933193*tmp_5 + tmp_6);
      real_t tmp_41 = tmp_1*tmp_40;
      real_t tmp_42 = tmp_40*tmp_9;
      real_t tmp_43 = tmp_4*(0.95308992296933193*tmp_11 + tmp_12);
      real_t tmp_44 = tmp_3*tmp_43;
      real_t tmp_45 = tmp_15*tmp_43;
      real_t tmp_46 = 0.11846344252809471*tmp_17;
      real_t a_0_0 = tmp_18*(-tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1) + tmp_25*(-tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1) + tmp_32*(-tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1) + tmp_39*(-tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1) + tmp_46*(-tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1);
      real_t a_1_0 = tmp_18*(tmp_10 + tmp_14) + tmp_25*(tmp_21 + tmp_23) + tmp_32*(tmp_28 + tmp_30) + tmp_39*(tmp_35 + tmp_37) + tmp_46*(tmp_42 + tmp_44);
      real_t a_2_0 = tmp_18*(tmp_16 + tmp_8) + tmp_25*(tmp_20 + tmp_24) + tmp_32*(tmp_27 + tmp_31) + tmp_39*(tmp_34 + tmp_38) + tmp_46*(tmp_41 + tmp_45);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
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


      // Does nothing.
      elMat( 0, 0) = 0;
      elMat( 1, 0) = 0;
      elMat( 2, 0) = 0;
   }

public:


};




class EGDivtFormEDGP0_new : public hyteg::dg::DGForm2D
{
 public:
    EGDivtFormEDGP0_new()

    {}



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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = -p_affine_0_1;
      real_t tmp_2 = (p_affine_1_0 + tmp_0)*(p_affine_2_1 + tmp_1);
      real_t tmp_3 = p_affine_2_0 + tmp_0;
      real_t tmp_4 = p_affine_1_1 + tmp_1;
      real_t tmp_5 = 1.0 / (tmp_2 - tmp_3*tmp_4);
      real_t tmp_6 = (-2*tmp_2*tmp_5 - tmp_3*tmp_5*(p_affine_0_1 - p_affine_1_1) - tmp_4*tmp_5*(p_affine_0_0 - p_affine_2_0))*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = 0.5*tmp_6;
      elMat( 0, 0) = a_0_0;
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

      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = p_affine_1_0 + tmp_3;
      real_t tmp_5 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_6 = -p_affine_0_1;
      real_t tmp_7 = p_affine_2_1 + tmp_6;
      real_t tmp_8 = p_affine_2_0 + tmp_3;
      real_t tmp_9 = p_affine_1_1 + tmp_6;
      real_t tmp_10 = 1.0 / (tmp_4*tmp_7 - tmp_8*tmp_9);
      real_t tmp_11 = p_affine_6_1 + tmp_6;
      real_t tmp_12 = tmp_10*(0.046910077030668018*tmp_1 + tmp_11);
      real_t tmp_13 = p_affine_6_0 + tmp_3;
      real_t tmp_14 = tmp_10*(0.046910077030668018*tmp_0 + tmp_13);
      real_t tmp_15 = tmp_12*tmp_5 + tmp_14*tmp_7 - 1.0/3.0;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_12*tmp_4 + tmp_14*tmp_16 - 1.0/3.0;
      real_t tmp_18 = 0.5*p_affine_10_0;
      real_t tmp_19 = 0.5*p_affine_10_1;
      real_t tmp_20 = tmp_10*(0.23076534494715845*tmp_1 + tmp_11);
      real_t tmp_21 = tmp_10*(0.23076534494715845*tmp_0 + tmp_13);
      real_t tmp_22 = tmp_20*tmp_5 + tmp_21*tmp_7 - 1.0/3.0;
      real_t tmp_23 = tmp_16*tmp_21 + tmp_20*tmp_4 - 1.0/3.0;
      real_t tmp_24 = tmp_10*(0.5*tmp_1 + tmp_11);
      real_t tmp_25 = tmp_10*(0.5*tmp_0 + tmp_13);
      real_t tmp_26 = tmp_24*tmp_5 + tmp_25*tmp_7 - 1.0/3.0;
      real_t tmp_27 = tmp_16*tmp_25 + tmp_24*tmp_4 - 1.0/3.0;
      real_t tmp_28 = tmp_10*(0.7692346550528415*tmp_1 + tmp_11);
      real_t tmp_29 = tmp_10*(0.7692346550528415*tmp_0 + tmp_13);
      real_t tmp_30 = tmp_28*tmp_5 + tmp_29*tmp_7 - 1.0/3.0;
      real_t tmp_31 = tmp_16*tmp_29 + tmp_28*tmp_4 - 1.0/3.0;
      real_t tmp_32 = tmp_10*(0.95308992296933193*tmp_1 + tmp_11);
      real_t tmp_33 = tmp_10*(0.95308992296933193*tmp_0 + tmp_13);
      real_t tmp_34 = tmp_32*tmp_5 + tmp_33*tmp_7 - 1.0/3.0;
      real_t tmp_35 = tmp_16*tmp_33 + tmp_32*tmp_4 - 1.0/3.0;
      real_t a_0_0 = 0.11846344252809471*tmp_2*(tmp_18*(tmp_15*tmp_4 + tmp_17*tmp_8) + tmp_19*(tmp_15*tmp_9 + tmp_17*tmp_7)) + 0.2393143352496831*tmp_2*(tmp_18*(tmp_22*tmp_4 + tmp_23*tmp_8) + tmp_19*(tmp_22*tmp_9 + tmp_23*tmp_7)) + 0.2844444444444445*tmp_2*(tmp_18*(tmp_26*tmp_4 + tmp_27*tmp_8) + tmp_19*(tmp_26*tmp_9 + tmp_27*tmp_7)) + 0.2393143352496831*tmp_2*(tmp_18*(tmp_30*tmp_4 + tmp_31*tmp_8) + tmp_19*(tmp_30*tmp_9 + tmp_31*tmp_7)) + 0.11846344252809471*tmp_2*(tmp_18*(tmp_34*tmp_4 + tmp_35*tmp_8) + tmp_19*(tmp_34*tmp_9 + tmp_35*tmp_7));
      elMat( 0, 0) = a_0_0;
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

      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = p_affine_1_0 + tmp_3;
      real_t tmp_5 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_6 = -p_affine_0_1;
      real_t tmp_7 = p_affine_2_1 + tmp_6;
      real_t tmp_8 = p_affine_2_0 + tmp_3;
      real_t tmp_9 = p_affine_1_1 + tmp_6;
      real_t tmp_10 = 1.0 / (tmp_4*tmp_7 - tmp_8*tmp_9);
      real_t tmp_11 = p_affine_6_1 + tmp_6;
      real_t tmp_12 = tmp_10*(0.046910077030668018*tmp_1 + tmp_11);
      real_t tmp_13 = p_affine_6_0 + tmp_3;
      real_t tmp_14 = tmp_10*(0.046910077030668018*tmp_0 + tmp_13);
      real_t tmp_15 = tmp_12*tmp_5 + tmp_14*tmp_7 - 1.0/3.0;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_12*tmp_4 + tmp_14*tmp_16 - 1.0/3.0;
      real_t tmp_18 = 0.5*p_affine_10_0;
      real_t tmp_19 = 0.5*p_affine_10_1;
      real_t tmp_20 = tmp_10*(0.23076534494715845*tmp_1 + tmp_11);
      real_t tmp_21 = tmp_10*(0.23076534494715845*tmp_0 + tmp_13);
      real_t tmp_22 = tmp_20*tmp_5 + tmp_21*tmp_7 - 1.0/3.0;
      real_t tmp_23 = tmp_16*tmp_21 + tmp_20*tmp_4 - 1.0/3.0;
      real_t tmp_24 = tmp_10*(0.5*tmp_1 + tmp_11);
      real_t tmp_25 = tmp_10*(0.5*tmp_0 + tmp_13);
      real_t tmp_26 = tmp_24*tmp_5 + tmp_25*tmp_7 - 1.0/3.0;
      real_t tmp_27 = tmp_16*tmp_25 + tmp_24*tmp_4 - 1.0/3.0;
      real_t tmp_28 = tmp_10*(0.7692346550528415*tmp_1 + tmp_11);
      real_t tmp_29 = tmp_10*(0.7692346550528415*tmp_0 + tmp_13);
      real_t tmp_30 = tmp_28*tmp_5 + tmp_29*tmp_7 - 1.0/3.0;
      real_t tmp_31 = tmp_16*tmp_29 + tmp_28*tmp_4 - 1.0/3.0;
      real_t tmp_32 = tmp_10*(0.95308992296933193*tmp_1 + tmp_11);
      real_t tmp_33 = tmp_10*(0.95308992296933193*tmp_0 + tmp_13);
      real_t tmp_34 = tmp_32*tmp_5 + tmp_33*tmp_7 - 1.0/3.0;
      real_t tmp_35 = tmp_16*tmp_33 + tmp_32*tmp_4 - 1.0/3.0;
      real_t a_0_0 = 0.11846344252809471*tmp_2*(tmp_18*(tmp_15*tmp_4 + tmp_17*tmp_8) + tmp_19*(tmp_15*tmp_9 + tmp_17*tmp_7)) + 0.2393143352496831*tmp_2*(tmp_18*(tmp_22*tmp_4 + tmp_23*tmp_8) + tmp_19*(tmp_22*tmp_9 + tmp_23*tmp_7)) + 0.2844444444444445*tmp_2*(tmp_18*(tmp_26*tmp_4 + tmp_27*tmp_8) + tmp_19*(tmp_26*tmp_9 + tmp_27*tmp_7)) + 0.2393143352496831*tmp_2*(tmp_18*(tmp_30*tmp_4 + tmp_31*tmp_8) + tmp_19*(tmp_30*tmp_9 + tmp_31*tmp_7)) + 0.11846344252809471*tmp_2*(tmp_18*(tmp_34*tmp_4 + tmp_35*tmp_8) + tmp_19*(tmp_34*tmp_9 + tmp_35*tmp_7));
      elMat( 0, 0) = a_0_0;
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

      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = p_affine_1_0 + tmp_3;
      real_t tmp_5 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_6 = -p_affine_0_1;
      real_t tmp_7 = p_affine_2_1 + tmp_6;
      real_t tmp_8 = p_affine_2_0 + tmp_3;
      real_t tmp_9 = p_affine_1_1 + tmp_6;
      real_t tmp_10 = 1.0 / (tmp_4*tmp_7 - tmp_8*tmp_9);
      real_t tmp_11 = p_affine_6_1 + tmp_6;
      real_t tmp_12 = tmp_10*(0.046910077030668018*tmp_1 + tmp_11);
      real_t tmp_13 = p_affine_6_0 + tmp_3;
      real_t tmp_14 = tmp_10*(0.046910077030668018*tmp_0 + tmp_13);
      real_t tmp_15 = tmp_12*tmp_5 + tmp_14*tmp_7 - 1.0/3.0;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_12*tmp_4 + tmp_14*tmp_16 - 1.0/3.0;
      real_t tmp_18 = tmp_10*(0.23076534494715845*tmp_1 + tmp_11);
      real_t tmp_19 = tmp_10*(0.23076534494715845*tmp_0 + tmp_13);
      real_t tmp_20 = tmp_18*tmp_5 + tmp_19*tmp_7 - 1.0/3.0;
      real_t tmp_21 = tmp_16*tmp_19 + tmp_18*tmp_4 - 1.0/3.0;
      real_t tmp_22 = tmp_10*(0.5*tmp_1 + tmp_11);
      real_t tmp_23 = tmp_10*(0.5*tmp_0 + tmp_13);
      real_t tmp_24 = tmp_22*tmp_5 + tmp_23*tmp_7 - 1.0/3.0;
      real_t tmp_25 = tmp_16*tmp_23 + tmp_22*tmp_4 - 1.0/3.0;
      real_t tmp_26 = tmp_10*(0.7692346550528415*tmp_1 + tmp_11);
      real_t tmp_27 = tmp_10*(0.7692346550528415*tmp_0 + tmp_13);
      real_t tmp_28 = tmp_26*tmp_5 + tmp_27*tmp_7 - 1.0/3.0;
      real_t tmp_29 = tmp_16*tmp_27 + tmp_26*tmp_4 - 1.0/3.0;
      real_t tmp_30 = tmp_10*(0.95308992296933193*tmp_1 + tmp_11);
      real_t tmp_31 = tmp_10*(0.95308992296933193*tmp_0 + tmp_13);
      real_t tmp_32 = tmp_30*tmp_5 + tmp_31*tmp_7 - 1.0/3.0;
      real_t tmp_33 = tmp_16*tmp_31 + tmp_30*tmp_4 - 1.0/3.0;
      real_t a_0_0 = 0.11846344252809471*tmp_2*(p_affine_10_0*(tmp_15*tmp_4 + tmp_17*tmp_8) + p_affine_10_1*(tmp_15*tmp_9 + tmp_17*tmp_7)) + 0.2393143352496831*tmp_2*(p_affine_10_0*(tmp_20*tmp_4 + tmp_21*tmp_8) + p_affine_10_1*(tmp_20*tmp_9 + tmp_21*tmp_7)) + 0.2844444444444445*tmp_2*(p_affine_10_0*(tmp_24*tmp_4 + tmp_25*tmp_8) + p_affine_10_1*(tmp_24*tmp_9 + tmp_25*tmp_7)) + 0.2393143352496831*tmp_2*(p_affine_10_0*(tmp_28*tmp_4 + tmp_29*tmp_8) + p_affine_10_1*(tmp_28*tmp_9 + tmp_29*tmp_7)) + 0.11846344252809471*tmp_2*(p_affine_10_0*(tmp_32*tmp_4 + tmp_33*tmp_8) + p_affine_10_1*(tmp_32*tmp_9 + tmp_33*tmp_7));
      elMat( 0, 0) = a_0_0;
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


      // Does nothing.
      elMat( 0, 0) = 0;
   }

public:


};


} // dg
} // hyteg

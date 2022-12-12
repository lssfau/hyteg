
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
# include "hyteg/types/matrix.hpp"
# include "hyteg/types/pointnd.hpp"

# include "Eigen/Eigen"

namespace hyteg {
namespace dg{
namespace eg{

class EGDivtForm_P1P0_0 : public hyteg::dg::DGForm
{

 public:
    EGDivtForm_P1P0_0()

    {}





 protected:
  void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      elMat( 0, 0) = 0;
      elMat( 1, 0) = 0;
      elMat( 2, 0) = 0;
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
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

      elMat( 0, 0) = 0;
      elMat( 1, 0) = 0;
      elMat( 2, 0) = 0;
      elMat( 3, 0) = 0;
   }
   void integrateVolume3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_1_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_2;
      real_t tmp_3 = p_affine_2_2 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_1 + tmp_0;
      real_t tmp_6 = p_affine_1_2 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = -p_affine_0_0;
      real_t tmp_9 = p_affine_1_0 + tmp_8;
      real_t tmp_10 = p_affine_3_2 + tmp_2;
      real_t tmp_11 = tmp_10*tmp_5;
      real_t tmp_12 = p_affine_2_0 + tmp_8;
      real_t tmp_13 = p_affine_3_1 + tmp_0;
      real_t tmp_14 = tmp_13*tmp_6;
      real_t tmp_15 = p_affine_3_0 + tmp_8;
      real_t tmp_16 = tmp_13*tmp_3;
      real_t tmp_17 = tmp_1*tmp_10;
      real_t tmp_18 = 1.0 / (tmp_11*tmp_9 + tmp_12*tmp_14 - tmp_12*tmp_17 + tmp_15*tmp_4 - tmp_15*tmp_7 - tmp_16*tmp_9);
      real_t tmp_19 = tmp_18*(tmp_4 - tmp_7);
      real_t tmp_20 = tmp_18*(tmp_14 - tmp_17);
      real_t tmp_21 = tmp_18*(tmp_11 - tmp_16);
      real_t tmp_22 = p_affine_0_0*p_affine_1_1;
      real_t tmp_23 = p_affine_0_0*p_affine_1_2;
      real_t tmp_24 = p_affine_2_1*p_affine_3_2;
      real_t tmp_25 = p_affine_0_1*p_affine_1_0;
      real_t tmp_26 = p_affine_0_1*p_affine_1_2;
      real_t tmp_27 = p_affine_2_2*p_affine_3_0;
      real_t tmp_28 = p_affine_0_2*p_affine_1_0;
      real_t tmp_29 = p_affine_0_2*p_affine_1_1;
      real_t tmp_30 = p_affine_2_0*p_affine_3_1;
      real_t tmp_31 = p_affine_2_2*p_affine_3_1;
      real_t tmp_32 = p_affine_2_0*p_affine_3_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_0;
      real_t tmp_34 = std::abs(p_affine_0_0*tmp_24 - p_affine_0_0*tmp_31 + p_affine_0_1*tmp_27 - p_affine_0_1*tmp_32 + p_affine_0_2*tmp_30 - p_affine_0_2*tmp_33 - p_affine_1_0*tmp_24 + p_affine_1_0*tmp_31 - p_affine_1_1*tmp_27 + p_affine_1_1*tmp_32 - p_affine_1_2*tmp_30 + p_affine_1_2*tmp_33 + p_affine_2_0*tmp_26 - p_affine_2_0*tmp_29 - p_affine_2_1*tmp_23 + p_affine_2_1*tmp_28 + p_affine_2_2*tmp_22 - p_affine_2_2*tmp_25 - p_affine_3_0*tmp_26 + p_affine_3_0*tmp_29 + p_affine_3_1*tmp_23 - p_affine_3_1*tmp_28 - p_affine_3_2*tmp_22 + p_affine_3_2*tmp_25);
      real_t tmp_35 = tmp_34*(tmp_19 + tmp_20 + tmp_21);
      real_t tmp_36 = tmp_21*tmp_34;
      real_t tmp_37 = tmp_20*tmp_34;
      real_t tmp_38 = tmp_19*tmp_34;
      real_t a_0_0 = 0.1666666666666668*tmp_35;
      real_t a_1_0 = -0.1666666666666668*tmp_36;
      real_t a_2_0 = -0.1666666666666668*tmp_37;
      real_t a_3_0 = -0.1666666666666668*tmp_38;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }



   void integrateFacetInner3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
                                                     const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                     const Eigen::Matrix< real_t, 3, 1 >&,
                                                     const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                                                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

         real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = p_affine_1_2 + tmp_9;
      real_t tmp_13 = tmp_12*tmp_5;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = p_affine_2_2 + tmp_9;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_1*tmp_11;
      real_t tmp_18 = tmp_12*tmp_14;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_4 - tmp_10*tmp_7 + tmp_11*tmp_13 + tmp_14*tmp_16 - tmp_15*tmp_17 - tmp_18*tmp_3);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_9;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_8;
      real_t tmp_26 = tmp_14*tmp_6 - tmp_17;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_15 + tmp_13;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_10 - tmp_18;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_11*tmp_5 - tmp_14*tmp_3;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_12*tmp_3 + tmp_16;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = tmp_10*tmp_3 - tmp_11*tmp_15;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_53 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_54 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_55 = 0.5*p_affine_13_0*std::pow((std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)*std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)) + (std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)*std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)) + (std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)*std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)), 1.0/2.0);
      real_t tmp_56 = 0.0068572537431980923*tmp_55;
      real_t tmp_57 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_58 = tmp_57*tmp_8;
      real_t tmp_59 = tmp_26*tmp_57;
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_61 = tmp_28*tmp_60;
      real_t tmp_62 = tmp_35*tmp_60;
      real_t tmp_63 = tmp_37*tmp_57;
      real_t tmp_64 = tmp_39*tmp_60;
      real_t tmp_65 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_66 = tmp_41*tmp_65;
      real_t tmp_67 = tmp_48*tmp_65;
      real_t tmp_68 = tmp_50*tmp_65;
      real_t tmp_69 = 0.037198804536718075*tmp_55;
      real_t tmp_70 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_71 = tmp_70*tmp_8;
      real_t tmp_72 = tmp_26*tmp_70;
      real_t tmp_73 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_74 = tmp_28*tmp_73;
      real_t tmp_75 = tmp_35*tmp_73;
      real_t tmp_76 = tmp_37*tmp_70;
      real_t tmp_77 = tmp_39*tmp_73;
      real_t tmp_78 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_79 = tmp_41*tmp_78;
      real_t tmp_80 = tmp_48*tmp_78;
      real_t tmp_81 = tmp_50*tmp_78;
      real_t tmp_82 = 0.020848748529055869*tmp_55;
      real_t tmp_83 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_84 = tmp_8*tmp_83;
      real_t tmp_85 = tmp_26*tmp_83;
      real_t tmp_86 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_87 = tmp_28*tmp_86;
      real_t tmp_88 = tmp_35*tmp_86;
      real_t tmp_89 = tmp_37*tmp_83;
      real_t tmp_90 = tmp_39*tmp_86;
      real_t tmp_91 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_92 = tmp_41*tmp_91;
      real_t tmp_93 = tmp_48*tmp_91;
      real_t tmp_94 = tmp_50*tmp_91;
      real_t tmp_95 = 0.019202922745021479*tmp_55;
      real_t tmp_96 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_97 = tmp_8*tmp_96;
      real_t tmp_98 = tmp_26*tmp_96;
      real_t tmp_99 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_100 = tmp_28*tmp_99;
      real_t tmp_101 = tmp_35*tmp_99;
      real_t tmp_102 = tmp_37*tmp_96;
      real_t tmp_103 = tmp_39*tmp_99;
      real_t tmp_104 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_105 = tmp_104*tmp_41;
      real_t tmp_106 = tmp_104*tmp_48;
      real_t tmp_107 = tmp_104*tmp_50;
      real_t tmp_108 = 0.020848748529055869*tmp_55;
      real_t tmp_109 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_110 = tmp_109*tmp_8;
      real_t tmp_111 = tmp_109*tmp_26;
      real_t tmp_112 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_113 = tmp_112*tmp_28;
      real_t tmp_114 = tmp_112*tmp_35;
      real_t tmp_115 = tmp_109*tmp_37;
      real_t tmp_116 = tmp_112*tmp_39;
      real_t tmp_117 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_118 = tmp_117*tmp_41;
      real_t tmp_119 = tmp_117*tmp_48;
      real_t tmp_120 = tmp_117*tmp_50;
      real_t tmp_121 = 0.019202922745021479*tmp_55;
      real_t tmp_122 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_123 = tmp_122*tmp_8;
      real_t tmp_124 = tmp_122*tmp_26;
      real_t tmp_125 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_126 = tmp_125*tmp_28;
      real_t tmp_127 = tmp_125*tmp_35;
      real_t tmp_128 = tmp_122*tmp_37;
      real_t tmp_129 = tmp_125*tmp_39;
      real_t tmp_130 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_131 = tmp_130*tmp_41;
      real_t tmp_132 = tmp_130*tmp_48;
      real_t tmp_133 = tmp_130*tmp_50;
      real_t tmp_134 = 0.020848748529055869*tmp_55;
      real_t tmp_135 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_136 = tmp_135*tmp_8;
      real_t tmp_137 = tmp_135*tmp_26;
      real_t tmp_138 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_139 = tmp_138*tmp_28;
      real_t tmp_140 = tmp_138*tmp_35;
      real_t tmp_141 = tmp_135*tmp_37;
      real_t tmp_142 = tmp_138*tmp_39;
      real_t tmp_143 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_144 = tmp_143*tmp_41;
      real_t tmp_145 = tmp_143*tmp_48;
      real_t tmp_146 = tmp_143*tmp_50;
      real_t tmp_147 = 0.019202922745021479*tmp_55;
      real_t tmp_148 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_149 = tmp_148*tmp_8;
      real_t tmp_150 = tmp_148*tmp_26;
      real_t tmp_151 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_152 = tmp_151*tmp_28;
      real_t tmp_153 = tmp_151*tmp_35;
      real_t tmp_154 = tmp_148*tmp_37;
      real_t tmp_155 = tmp_151*tmp_39;
      real_t tmp_156 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_157 = tmp_156*tmp_41;
      real_t tmp_158 = tmp_156*tmp_48;
      real_t tmp_159 = tmp_156*tmp_50;
      real_t tmp_160 = 0.020848748529055869*tmp_55;
      real_t tmp_161 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_162 = tmp_161*tmp_8;
      real_t tmp_163 = tmp_161*tmp_26;
      real_t tmp_164 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_165 = tmp_164*tmp_28;
      real_t tmp_166 = tmp_164*tmp_35;
      real_t tmp_167 = tmp_161*tmp_37;
      real_t tmp_168 = tmp_164*tmp_39;
      real_t tmp_169 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_170 = tmp_169*tmp_41;
      real_t tmp_171 = tmp_169*tmp_48;
      real_t tmp_172 = tmp_169*tmp_50;
      real_t tmp_173 = 0.019202922745021479*tmp_55;
      real_t tmp_174 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_175 = tmp_174*tmp_8;
      real_t tmp_176 = tmp_174*tmp_26;
      real_t tmp_177 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_178 = tmp_177*tmp_28;
      real_t tmp_179 = tmp_177*tmp_35;
      real_t tmp_180 = tmp_174*tmp_37;
      real_t tmp_181 = tmp_177*tmp_39;
      real_t tmp_182 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_183 = tmp_182*tmp_41;
      real_t tmp_184 = tmp_182*tmp_48;
      real_t tmp_185 = tmp_182*tmp_50;
      real_t tmp_186 = 0.020848748529055869*tmp_55;
      real_t tmp_187 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_188 = tmp_187*tmp_8;
      real_t tmp_189 = tmp_187*tmp_26;
      real_t tmp_190 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_191 = tmp_190*tmp_28;
      real_t tmp_192 = tmp_190*tmp_35;
      real_t tmp_193 = tmp_187*tmp_37;
      real_t tmp_194 = tmp_190*tmp_39;
      real_t tmp_195 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_196 = tmp_195*tmp_41;
      real_t tmp_197 = tmp_195*tmp_48;
      real_t tmp_198 = tmp_195*tmp_50;
      real_t tmp_199 = 0.019202922745021479*tmp_55;
      real_t tmp_200 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_201 = tmp_200*tmp_8;
      real_t tmp_202 = tmp_200*tmp_26;
      real_t tmp_203 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_204 = tmp_203*tmp_28;
      real_t tmp_205 = tmp_203*tmp_35;
      real_t tmp_206 = tmp_200*tmp_37;
      real_t tmp_207 = tmp_203*tmp_39;
      real_t tmp_208 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_209 = tmp_208*tmp_41;
      real_t tmp_210 = tmp_208*tmp_48;
      real_t tmp_211 = tmp_208*tmp_50;
      real_t tmp_212 = 0.042507265838595799*tmp_55;
      real_t tmp_213 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_214 = tmp_213*tmp_8;
      real_t tmp_215 = tmp_213*tmp_26;
      real_t tmp_216 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_217 = tmp_216*tmp_28;
      real_t tmp_218 = tmp_216*tmp_35;
      real_t tmp_219 = tmp_213*tmp_37;
      real_t tmp_220 = tmp_216*tmp_39;
      real_t tmp_221 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_222 = tmp_221*tmp_41;
      real_t tmp_223 = tmp_221*tmp_48;
      real_t tmp_224 = tmp_221*tmp_50;
      real_t tmp_225 = 0.020848748529055869*tmp_55;
      real_t tmp_226 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_227 = tmp_226*tmp_8;
      real_t tmp_228 = tmp_226*tmp_26;
      real_t tmp_229 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_230 = tmp_229*tmp_28;
      real_t tmp_231 = tmp_229*tmp_35;
      real_t tmp_232 = tmp_226*tmp_37;
      real_t tmp_233 = tmp_229*tmp_39;
      real_t tmp_234 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_235 = tmp_234*tmp_41;
      real_t tmp_236 = tmp_234*tmp_48;
      real_t tmp_237 = tmp_234*tmp_50;
      real_t tmp_238 = 0.0068572537431980923*tmp_55;
      real_t tmp_239 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_240 = tmp_239*tmp_8;
      real_t tmp_241 = tmp_239*tmp_26;
      real_t tmp_242 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_243 = tmp_242*tmp_28;
      real_t tmp_244 = tmp_242*tmp_35;
      real_t tmp_245 = tmp_239*tmp_37;
      real_t tmp_246 = tmp_242*tmp_39;
      real_t tmp_247 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_248 = tmp_247*tmp_41;
      real_t tmp_249 = tmp_247*tmp_48;
      real_t tmp_250 = tmp_247*tmp_50;
      real_t tmp_251 = 0.037198804536718075*tmp_55;
      real_t tmp_252 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_8;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = 0.042507265838595799*tmp_55;
      real_t tmp_265 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_266 = tmp_265*tmp_8;
      real_t tmp_267 = tmp_26*tmp_265;
      real_t tmp_268 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_269 = tmp_268*tmp_28;
      real_t tmp_270 = tmp_268*tmp_35;
      real_t tmp_271 = tmp_265*tmp_37;
      real_t tmp_272 = tmp_268*tmp_39;
      real_t tmp_273 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_274 = tmp_273*tmp_41;
      real_t tmp_275 = tmp_273*tmp_48;
      real_t tmp_276 = tmp_273*tmp_50;
      real_t tmp_277 = 0.0068572537431980923*tmp_55;
      real_t tmp_278 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_279 = tmp_278*tmp_8;
      real_t tmp_280 = tmp_26*tmp_278;
      real_t tmp_281 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_282 = tmp_28*tmp_281;
      real_t tmp_283 = tmp_281*tmp_35;
      real_t tmp_284 = tmp_278*tmp_37;
      real_t tmp_285 = tmp_281*tmp_39;
      real_t tmp_286 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_287 = tmp_286*tmp_41;
      real_t tmp_288 = tmp_286*tmp_48;
      real_t tmp_289 = tmp_286*tmp_50;
      real_t tmp_290 = 0.037198804536718075*tmp_55;
      real_t tmp_291 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_292 = tmp_291*tmp_8;
      real_t tmp_293 = tmp_26*tmp_291;
      real_t tmp_294 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_295 = tmp_28*tmp_294;
      real_t tmp_296 = tmp_294*tmp_35;
      real_t tmp_297 = tmp_291*tmp_37;
      real_t tmp_298 = tmp_294*tmp_39;
      real_t tmp_299 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_300 = tmp_299*tmp_41;
      real_t tmp_301 = tmp_299*tmp_48;
      real_t tmp_302 = tmp_299*tmp_50;
      real_t tmp_303 = 0.042507265838595799*tmp_55;
      real_t tmp_304 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_305 = tmp_304*tmp_8;
      real_t tmp_306 = tmp_26*tmp_304;
      real_t tmp_307 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_308 = tmp_28*tmp_307;
      real_t tmp_309 = tmp_307*tmp_35;
      real_t tmp_310 = tmp_304*tmp_37;
      real_t tmp_311 = tmp_307*tmp_39;
      real_t tmp_312 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_313 = tmp_312*tmp_41;
      real_t tmp_314 = tmp_312*tmp_48;
      real_t tmp_315 = tmp_312*tmp_50;
      real_t tmp_316 = 0.019202922745021479*tmp_55;
      real_t a_0_0 = tmp_108*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_105 - tmp_106 - tmp_107 - tmp_97 - tmp_98 + 1) + tmp_121*(-tmp_110 - tmp_111 - tmp_113 - tmp_114 - tmp_115 - tmp_116 - tmp_118 - tmp_119 - tmp_120 + 1) + tmp_134*(-tmp_123 - tmp_124 - tmp_126 - tmp_127 - tmp_128 - tmp_129 - tmp_131 - tmp_132 - tmp_133 + 1) + tmp_147*(-tmp_136 - tmp_137 - tmp_139 - tmp_140 - tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 + 1) + tmp_160*(-tmp_149 - tmp_150 - tmp_152 - tmp_153 - tmp_154 - tmp_155 - tmp_157 - tmp_158 - tmp_159 + 1) + tmp_173*(-tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 - tmp_168 - tmp_170 - tmp_171 - tmp_172 + 1) + tmp_186*(-tmp_175 - tmp_176 - tmp_178 - tmp_179 - tmp_180 - tmp_181 - tmp_183 - tmp_184 - tmp_185 + 1) + tmp_199*(-tmp_188 - tmp_189 - tmp_191 - tmp_192 - tmp_193 - tmp_194 - tmp_196 - tmp_197 - tmp_198 + 1) + tmp_212*(-tmp_201 - tmp_202 - tmp_204 - tmp_205 - tmp_206 - tmp_207 - tmp_209 - tmp_210 - tmp_211 + 1) + tmp_225*(-tmp_214 - tmp_215 - tmp_217 - tmp_218 - tmp_219 - tmp_220 - tmp_222 - tmp_223 - tmp_224 + 1) + tmp_238*(-tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 - tmp_233 - tmp_235 - tmp_236 - tmp_237 + 1) + tmp_251*(-tmp_240 - tmp_241 - tmp_243 - tmp_244 - tmp_245 - tmp_246 - tmp_248 - tmp_249 - tmp_250 + 1) + tmp_264*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_277*(-tmp_266 - tmp_267 - tmp_269 - tmp_270 - tmp_271 - tmp_272 - tmp_274 - tmp_275 - tmp_276 + 1) + tmp_290*(-tmp_279 - tmp_280 - tmp_282 - tmp_283 - tmp_284 - tmp_285 - tmp_287 - tmp_288 - tmp_289 + 1) + tmp_303*(-tmp_292 - tmp_293 - tmp_295 - tmp_296 - tmp_297 - tmp_298 - tmp_300 - tmp_301 - tmp_302 + 1) + tmp_316*(-tmp_305 - tmp_306 - tmp_308 - tmp_309 - tmp_310 - tmp_311 - tmp_313 - tmp_314 - tmp_315 + 1) + tmp_56*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_69*(-tmp_58 - tmp_59 - tmp_61 - tmp_62 - tmp_63 - tmp_64 - tmp_66 - tmp_67 - tmp_68 + 1) + tmp_82*(-tmp_71 - tmp_72 - tmp_74 - tmp_75 - tmp_76 - tmp_77 - tmp_79 - tmp_80 - tmp_81 + 1) + tmp_95*(-tmp_84 - tmp_85 - tmp_87 - tmp_88 - tmp_89 - tmp_90 - tmp_92 - tmp_93 - tmp_94 + 1);
      real_t a_1_0 = tmp_108*(tmp_102 + tmp_103 + tmp_107) + tmp_121*(tmp_115 + tmp_116 + tmp_120) + tmp_134*(tmp_128 + tmp_129 + tmp_133) + tmp_147*(tmp_141 + tmp_142 + tmp_146) + tmp_160*(tmp_154 + tmp_155 + tmp_159) + tmp_173*(tmp_167 + tmp_168 + tmp_172) + tmp_186*(tmp_180 + tmp_181 + tmp_185) + tmp_199*(tmp_193 + tmp_194 + tmp_198) + tmp_212*(tmp_206 + tmp_207 + tmp_211) + tmp_225*(tmp_219 + tmp_220 + tmp_224) + tmp_238*(tmp_232 + tmp_233 + tmp_237) + tmp_251*(tmp_245 + tmp_246 + tmp_250) + tmp_264*(tmp_258 + tmp_259 + tmp_263) + tmp_277*(tmp_271 + tmp_272 + tmp_276) + tmp_290*(tmp_284 + tmp_285 + tmp_289) + tmp_303*(tmp_297 + tmp_298 + tmp_302) + tmp_316*(tmp_310 + tmp_311 + tmp_315) + tmp_56*(tmp_38 + tmp_40 + tmp_51) + tmp_69*(tmp_63 + tmp_64 + tmp_68) + tmp_82*(tmp_76 + tmp_77 + tmp_81) + tmp_95*(tmp_89 + tmp_90 + tmp_94);
      real_t a_2_0 = tmp_108*(tmp_101 + tmp_106 + tmp_98) + tmp_121*(tmp_111 + tmp_114 + tmp_119) + tmp_134*(tmp_124 + tmp_127 + tmp_132) + tmp_147*(tmp_137 + tmp_140 + tmp_145) + tmp_160*(tmp_150 + tmp_153 + tmp_158) + tmp_173*(tmp_163 + tmp_166 + tmp_171) + tmp_186*(tmp_176 + tmp_179 + tmp_184) + tmp_199*(tmp_189 + tmp_192 + tmp_197) + tmp_212*(tmp_202 + tmp_205 + tmp_210) + tmp_225*(tmp_215 + tmp_218 + tmp_223) + tmp_238*(tmp_228 + tmp_231 + tmp_236) + tmp_251*(tmp_241 + tmp_244 + tmp_249) + tmp_264*(tmp_254 + tmp_257 + tmp_262) + tmp_277*(tmp_267 + tmp_270 + tmp_275) + tmp_290*(tmp_280 + tmp_283 + tmp_288) + tmp_303*(tmp_293 + tmp_296 + tmp_301) + tmp_316*(tmp_306 + tmp_309 + tmp_314) + tmp_56*(tmp_27 + tmp_36 + tmp_49) + tmp_69*(tmp_59 + tmp_62 + tmp_67) + tmp_82*(tmp_72 + tmp_75 + tmp_80) + tmp_95*(tmp_85 + tmp_88 + tmp_93);
      real_t a_3_0 = tmp_108*(tmp_100 + tmp_105 + tmp_97) + tmp_121*(tmp_110 + tmp_113 + tmp_118) + tmp_134*(tmp_123 + tmp_126 + tmp_131) + tmp_147*(tmp_136 + tmp_139 + tmp_144) + tmp_160*(tmp_149 + tmp_152 + tmp_157) + tmp_173*(tmp_162 + tmp_165 + tmp_170) + tmp_186*(tmp_175 + tmp_178 + tmp_183) + tmp_199*(tmp_188 + tmp_191 + tmp_196) + tmp_212*(tmp_201 + tmp_204 + tmp_209) + tmp_225*(tmp_214 + tmp_217 + tmp_222) + tmp_238*(tmp_227 + tmp_230 + tmp_235) + tmp_251*(tmp_240 + tmp_243 + tmp_248) + tmp_264*(tmp_253 + tmp_256 + tmp_261) + tmp_277*(tmp_266 + tmp_269 + tmp_274) + tmp_290*(tmp_279 + tmp_282 + tmp_287) + tmp_303*(tmp_292 + tmp_295 + tmp_300) + tmp_316*(tmp_305 + tmp_308 + tmp_313) + tmp_56*(tmp_25 + tmp_34 + tmp_47) + tmp_69*(tmp_58 + tmp_61 + tmp_66) + tmp_82*(tmp_71 + tmp_74 + tmp_79) + tmp_95*(tmp_84 + tmp_87 + tmp_92);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }




void integrateFacetCoupling3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementInner,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementOuter,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                                        Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = p_affine_1_2 + tmp_9;
      real_t tmp_13 = tmp_12*tmp_5;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = p_affine_2_2 + tmp_9;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_1*tmp_11;
      real_t tmp_18 = tmp_12*tmp_14;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_4 - tmp_10*tmp_7 + tmp_11*tmp_13 + tmp_14*tmp_16 - tmp_15*tmp_17 - tmp_18*tmp_3);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_9;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_8;
      real_t tmp_26 = tmp_14*tmp_6 - tmp_17;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_15 + tmp_13;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_10 - tmp_18;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_11*tmp_5 - tmp_14*tmp_3;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_12*tmp_3 + tmp_16;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = tmp_10*tmp_3 - tmp_11*tmp_15;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_53 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_54 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_55 = 0.5*p_affine_13_0*std::pow((std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)*std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)) + (std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)*std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)) + (std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)*std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)), 1.0/2.0);
      real_t tmp_56 = 0.0068572537431980923*tmp_55;
      real_t tmp_57 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_58 = tmp_57*tmp_8;
      real_t tmp_59 = tmp_26*tmp_57;
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_61 = tmp_28*tmp_60;
      real_t tmp_62 = tmp_35*tmp_60;
      real_t tmp_63 = tmp_37*tmp_57;
      real_t tmp_64 = tmp_39*tmp_60;
      real_t tmp_65 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_66 = tmp_41*tmp_65;
      real_t tmp_67 = tmp_48*tmp_65;
      real_t tmp_68 = tmp_50*tmp_65;
      real_t tmp_69 = 0.037198804536718075*tmp_55;
      real_t tmp_70 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_71 = tmp_70*tmp_8;
      real_t tmp_72 = tmp_26*tmp_70;
      real_t tmp_73 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_74 = tmp_28*tmp_73;
      real_t tmp_75 = tmp_35*tmp_73;
      real_t tmp_76 = tmp_37*tmp_70;
      real_t tmp_77 = tmp_39*tmp_73;
      real_t tmp_78 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_79 = tmp_41*tmp_78;
      real_t tmp_80 = tmp_48*tmp_78;
      real_t tmp_81 = tmp_50*tmp_78;
      real_t tmp_82 = 0.020848748529055869*tmp_55;
      real_t tmp_83 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_84 = tmp_8*tmp_83;
      real_t tmp_85 = tmp_26*tmp_83;
      real_t tmp_86 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_87 = tmp_28*tmp_86;
      real_t tmp_88 = tmp_35*tmp_86;
      real_t tmp_89 = tmp_37*tmp_83;
      real_t tmp_90 = tmp_39*tmp_86;
      real_t tmp_91 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_92 = tmp_41*tmp_91;
      real_t tmp_93 = tmp_48*tmp_91;
      real_t tmp_94 = tmp_50*tmp_91;
      real_t tmp_95 = 0.019202922745021479*tmp_55;
      real_t tmp_96 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_97 = tmp_8*tmp_96;
      real_t tmp_98 = tmp_26*tmp_96;
      real_t tmp_99 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_100 = tmp_28*tmp_99;
      real_t tmp_101 = tmp_35*tmp_99;
      real_t tmp_102 = tmp_37*tmp_96;
      real_t tmp_103 = tmp_39*tmp_99;
      real_t tmp_104 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_105 = tmp_104*tmp_41;
      real_t tmp_106 = tmp_104*tmp_48;
      real_t tmp_107 = tmp_104*tmp_50;
      real_t tmp_108 = 0.020848748529055869*tmp_55;
      real_t tmp_109 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_110 = tmp_109*tmp_8;
      real_t tmp_111 = tmp_109*tmp_26;
      real_t tmp_112 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_113 = tmp_112*tmp_28;
      real_t tmp_114 = tmp_112*tmp_35;
      real_t tmp_115 = tmp_109*tmp_37;
      real_t tmp_116 = tmp_112*tmp_39;
      real_t tmp_117 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_118 = tmp_117*tmp_41;
      real_t tmp_119 = tmp_117*tmp_48;
      real_t tmp_120 = tmp_117*tmp_50;
      real_t tmp_121 = 0.019202922745021479*tmp_55;
      real_t tmp_122 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_123 = tmp_122*tmp_8;
      real_t tmp_124 = tmp_122*tmp_26;
      real_t tmp_125 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_126 = tmp_125*tmp_28;
      real_t tmp_127 = tmp_125*tmp_35;
      real_t tmp_128 = tmp_122*tmp_37;
      real_t tmp_129 = tmp_125*tmp_39;
      real_t tmp_130 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_131 = tmp_130*tmp_41;
      real_t tmp_132 = tmp_130*tmp_48;
      real_t tmp_133 = tmp_130*tmp_50;
      real_t tmp_134 = 0.020848748529055869*tmp_55;
      real_t tmp_135 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_136 = tmp_135*tmp_8;
      real_t tmp_137 = tmp_135*tmp_26;
      real_t tmp_138 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_139 = tmp_138*tmp_28;
      real_t tmp_140 = tmp_138*tmp_35;
      real_t tmp_141 = tmp_135*tmp_37;
      real_t tmp_142 = tmp_138*tmp_39;
      real_t tmp_143 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_144 = tmp_143*tmp_41;
      real_t tmp_145 = tmp_143*tmp_48;
      real_t tmp_146 = tmp_143*tmp_50;
      real_t tmp_147 = 0.019202922745021479*tmp_55;
      real_t tmp_148 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_149 = tmp_148*tmp_8;
      real_t tmp_150 = tmp_148*tmp_26;
      real_t tmp_151 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_152 = tmp_151*tmp_28;
      real_t tmp_153 = tmp_151*tmp_35;
      real_t tmp_154 = tmp_148*tmp_37;
      real_t tmp_155 = tmp_151*tmp_39;
      real_t tmp_156 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_157 = tmp_156*tmp_41;
      real_t tmp_158 = tmp_156*tmp_48;
      real_t tmp_159 = tmp_156*tmp_50;
      real_t tmp_160 = 0.020848748529055869*tmp_55;
      real_t tmp_161 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_162 = tmp_161*tmp_8;
      real_t tmp_163 = tmp_161*tmp_26;
      real_t tmp_164 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_165 = tmp_164*tmp_28;
      real_t tmp_166 = tmp_164*tmp_35;
      real_t tmp_167 = tmp_161*tmp_37;
      real_t tmp_168 = tmp_164*tmp_39;
      real_t tmp_169 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_170 = tmp_169*tmp_41;
      real_t tmp_171 = tmp_169*tmp_48;
      real_t tmp_172 = tmp_169*tmp_50;
      real_t tmp_173 = 0.019202922745021479*tmp_55;
      real_t tmp_174 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_175 = tmp_174*tmp_8;
      real_t tmp_176 = tmp_174*tmp_26;
      real_t tmp_177 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_178 = tmp_177*tmp_28;
      real_t tmp_179 = tmp_177*tmp_35;
      real_t tmp_180 = tmp_174*tmp_37;
      real_t tmp_181 = tmp_177*tmp_39;
      real_t tmp_182 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_183 = tmp_182*tmp_41;
      real_t tmp_184 = tmp_182*tmp_48;
      real_t tmp_185 = tmp_182*tmp_50;
      real_t tmp_186 = 0.020848748529055869*tmp_55;
      real_t tmp_187 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_188 = tmp_187*tmp_8;
      real_t tmp_189 = tmp_187*tmp_26;
      real_t tmp_190 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_191 = tmp_190*tmp_28;
      real_t tmp_192 = tmp_190*tmp_35;
      real_t tmp_193 = tmp_187*tmp_37;
      real_t tmp_194 = tmp_190*tmp_39;
      real_t tmp_195 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_196 = tmp_195*tmp_41;
      real_t tmp_197 = tmp_195*tmp_48;
      real_t tmp_198 = tmp_195*tmp_50;
      real_t tmp_199 = 0.019202922745021479*tmp_55;
      real_t tmp_200 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_201 = tmp_200*tmp_8;
      real_t tmp_202 = tmp_200*tmp_26;
      real_t tmp_203 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_204 = tmp_203*tmp_28;
      real_t tmp_205 = tmp_203*tmp_35;
      real_t tmp_206 = tmp_200*tmp_37;
      real_t tmp_207 = tmp_203*tmp_39;
      real_t tmp_208 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_209 = tmp_208*tmp_41;
      real_t tmp_210 = tmp_208*tmp_48;
      real_t tmp_211 = tmp_208*tmp_50;
      real_t tmp_212 = 0.042507265838595799*tmp_55;
      real_t tmp_213 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_214 = tmp_213*tmp_8;
      real_t tmp_215 = tmp_213*tmp_26;
      real_t tmp_216 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_217 = tmp_216*tmp_28;
      real_t tmp_218 = tmp_216*tmp_35;
      real_t tmp_219 = tmp_213*tmp_37;
      real_t tmp_220 = tmp_216*tmp_39;
      real_t tmp_221 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_222 = tmp_221*tmp_41;
      real_t tmp_223 = tmp_221*tmp_48;
      real_t tmp_224 = tmp_221*tmp_50;
      real_t tmp_225 = 0.020848748529055869*tmp_55;
      real_t tmp_226 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_227 = tmp_226*tmp_8;
      real_t tmp_228 = tmp_226*tmp_26;
      real_t tmp_229 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_230 = tmp_229*tmp_28;
      real_t tmp_231 = tmp_229*tmp_35;
      real_t tmp_232 = tmp_226*tmp_37;
      real_t tmp_233 = tmp_229*tmp_39;
      real_t tmp_234 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_235 = tmp_234*tmp_41;
      real_t tmp_236 = tmp_234*tmp_48;
      real_t tmp_237 = tmp_234*tmp_50;
      real_t tmp_238 = 0.0068572537431980923*tmp_55;
      real_t tmp_239 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_240 = tmp_239*tmp_8;
      real_t tmp_241 = tmp_239*tmp_26;
      real_t tmp_242 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_243 = tmp_242*tmp_28;
      real_t tmp_244 = tmp_242*tmp_35;
      real_t tmp_245 = tmp_239*tmp_37;
      real_t tmp_246 = tmp_242*tmp_39;
      real_t tmp_247 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_248 = tmp_247*tmp_41;
      real_t tmp_249 = tmp_247*tmp_48;
      real_t tmp_250 = tmp_247*tmp_50;
      real_t tmp_251 = 0.037198804536718075*tmp_55;
      real_t tmp_252 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_8;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = 0.042507265838595799*tmp_55;
      real_t tmp_265 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_266 = tmp_265*tmp_8;
      real_t tmp_267 = tmp_26*tmp_265;
      real_t tmp_268 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_269 = tmp_268*tmp_28;
      real_t tmp_270 = tmp_268*tmp_35;
      real_t tmp_271 = tmp_265*tmp_37;
      real_t tmp_272 = tmp_268*tmp_39;
      real_t tmp_273 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_274 = tmp_273*tmp_41;
      real_t tmp_275 = tmp_273*tmp_48;
      real_t tmp_276 = tmp_273*tmp_50;
      real_t tmp_277 = 0.0068572537431980923*tmp_55;
      real_t tmp_278 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_279 = tmp_278*tmp_8;
      real_t tmp_280 = tmp_26*tmp_278;
      real_t tmp_281 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_282 = tmp_28*tmp_281;
      real_t tmp_283 = tmp_281*tmp_35;
      real_t tmp_284 = tmp_278*tmp_37;
      real_t tmp_285 = tmp_281*tmp_39;
      real_t tmp_286 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_287 = tmp_286*tmp_41;
      real_t tmp_288 = tmp_286*tmp_48;
      real_t tmp_289 = tmp_286*tmp_50;
      real_t tmp_290 = 0.037198804536718075*tmp_55;
      real_t tmp_291 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_292 = tmp_291*tmp_8;
      real_t tmp_293 = tmp_26*tmp_291;
      real_t tmp_294 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_295 = tmp_28*tmp_294;
      real_t tmp_296 = tmp_294*tmp_35;
      real_t tmp_297 = tmp_291*tmp_37;
      real_t tmp_298 = tmp_294*tmp_39;
      real_t tmp_299 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_300 = tmp_299*tmp_41;
      real_t tmp_301 = tmp_299*tmp_48;
      real_t tmp_302 = tmp_299*tmp_50;
      real_t tmp_303 = 0.042507265838595799*tmp_55;
      real_t tmp_304 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_305 = tmp_304*tmp_8;
      real_t tmp_306 = tmp_26*tmp_304;
      real_t tmp_307 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_308 = tmp_28*tmp_307;
      real_t tmp_309 = tmp_307*tmp_35;
      real_t tmp_310 = tmp_304*tmp_37;
      real_t tmp_311 = tmp_307*tmp_39;
      real_t tmp_312 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_313 = tmp_312*tmp_41;
      real_t tmp_314 = tmp_312*tmp_48;
      real_t tmp_315 = tmp_312*tmp_50;
      real_t tmp_316 = 0.019202922745021479*tmp_55;
      real_t a_0_0 = tmp_108*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_105 - tmp_106 - tmp_107 - tmp_97 - tmp_98 + 1) + tmp_121*(-tmp_110 - tmp_111 - tmp_113 - tmp_114 - tmp_115 - tmp_116 - tmp_118 - tmp_119 - tmp_120 + 1) + tmp_134*(-tmp_123 - tmp_124 - tmp_126 - tmp_127 - tmp_128 - tmp_129 - tmp_131 - tmp_132 - tmp_133 + 1) + tmp_147*(-tmp_136 - tmp_137 - tmp_139 - tmp_140 - tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 + 1) + tmp_160*(-tmp_149 - tmp_150 - tmp_152 - tmp_153 - tmp_154 - tmp_155 - tmp_157 - tmp_158 - tmp_159 + 1) + tmp_173*(-tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 - tmp_168 - tmp_170 - tmp_171 - tmp_172 + 1) + tmp_186*(-tmp_175 - tmp_176 - tmp_178 - tmp_179 - tmp_180 - tmp_181 - tmp_183 - tmp_184 - tmp_185 + 1) + tmp_199*(-tmp_188 - tmp_189 - tmp_191 - tmp_192 - tmp_193 - tmp_194 - tmp_196 - tmp_197 - tmp_198 + 1) + tmp_212*(-tmp_201 - tmp_202 - tmp_204 - tmp_205 - tmp_206 - tmp_207 - tmp_209 - tmp_210 - tmp_211 + 1) + tmp_225*(-tmp_214 - tmp_215 - tmp_217 - tmp_218 - tmp_219 - tmp_220 - tmp_222 - tmp_223 - tmp_224 + 1) + tmp_238*(-tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 - tmp_233 - tmp_235 - tmp_236 - tmp_237 + 1) + tmp_251*(-tmp_240 - tmp_241 - tmp_243 - tmp_244 - tmp_245 - tmp_246 - tmp_248 - tmp_249 - tmp_250 + 1) + tmp_264*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_277*(-tmp_266 - tmp_267 - tmp_269 - tmp_270 - tmp_271 - tmp_272 - tmp_274 - tmp_275 - tmp_276 + 1) + tmp_290*(-tmp_279 - tmp_280 - tmp_282 - tmp_283 - tmp_284 - tmp_285 - tmp_287 - tmp_288 - tmp_289 + 1) + tmp_303*(-tmp_292 - tmp_293 - tmp_295 - tmp_296 - tmp_297 - tmp_298 - tmp_300 - tmp_301 - tmp_302 + 1) + tmp_316*(-tmp_305 - tmp_306 - tmp_308 - tmp_309 - tmp_310 - tmp_311 - tmp_313 - tmp_314 - tmp_315 + 1) + tmp_56*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_69*(-tmp_58 - tmp_59 - tmp_61 - tmp_62 - tmp_63 - tmp_64 - tmp_66 - tmp_67 - tmp_68 + 1) + tmp_82*(-tmp_71 - tmp_72 - tmp_74 - tmp_75 - tmp_76 - tmp_77 - tmp_79 - tmp_80 - tmp_81 + 1) + tmp_95*(-tmp_84 - tmp_85 - tmp_87 - tmp_88 - tmp_89 - tmp_90 - tmp_92 - tmp_93 - tmp_94 + 1);
      real_t a_1_0 = tmp_108*(tmp_102 + tmp_103 + tmp_107) + tmp_121*(tmp_115 + tmp_116 + tmp_120) + tmp_134*(tmp_128 + tmp_129 + tmp_133) + tmp_147*(tmp_141 + tmp_142 + tmp_146) + tmp_160*(tmp_154 + tmp_155 + tmp_159) + tmp_173*(tmp_167 + tmp_168 + tmp_172) + tmp_186*(tmp_180 + tmp_181 + tmp_185) + tmp_199*(tmp_193 + tmp_194 + tmp_198) + tmp_212*(tmp_206 + tmp_207 + tmp_211) + tmp_225*(tmp_219 + tmp_220 + tmp_224) + tmp_238*(tmp_232 + tmp_233 + tmp_237) + tmp_251*(tmp_245 + tmp_246 + tmp_250) + tmp_264*(tmp_258 + tmp_259 + tmp_263) + tmp_277*(tmp_271 + tmp_272 + tmp_276) + tmp_290*(tmp_284 + tmp_285 + tmp_289) + tmp_303*(tmp_297 + tmp_298 + tmp_302) + tmp_316*(tmp_310 + tmp_311 + tmp_315) + tmp_56*(tmp_38 + tmp_40 + tmp_51) + tmp_69*(tmp_63 + tmp_64 + tmp_68) + tmp_82*(tmp_76 + tmp_77 + tmp_81) + tmp_95*(tmp_89 + tmp_90 + tmp_94);
      real_t a_2_0 = tmp_108*(tmp_101 + tmp_106 + tmp_98) + tmp_121*(tmp_111 + tmp_114 + tmp_119) + tmp_134*(tmp_124 + tmp_127 + tmp_132) + tmp_147*(tmp_137 + tmp_140 + tmp_145) + tmp_160*(tmp_150 + tmp_153 + tmp_158) + tmp_173*(tmp_163 + tmp_166 + tmp_171) + tmp_186*(tmp_176 + tmp_179 + tmp_184) + tmp_199*(tmp_189 + tmp_192 + tmp_197) + tmp_212*(tmp_202 + tmp_205 + tmp_210) + tmp_225*(tmp_215 + tmp_218 + tmp_223) + tmp_238*(tmp_228 + tmp_231 + tmp_236) + tmp_251*(tmp_241 + tmp_244 + tmp_249) + tmp_264*(tmp_254 + tmp_257 + tmp_262) + tmp_277*(tmp_267 + tmp_270 + tmp_275) + tmp_290*(tmp_280 + tmp_283 + tmp_288) + tmp_303*(tmp_293 + tmp_296 + tmp_301) + tmp_316*(tmp_306 + tmp_309 + tmp_314) + tmp_56*(tmp_27 + tmp_36 + tmp_49) + tmp_69*(tmp_59 + tmp_62 + tmp_67) + tmp_82*(tmp_72 + tmp_75 + tmp_80) + tmp_95*(tmp_85 + tmp_88 + tmp_93);
      real_t a_3_0 = tmp_108*(tmp_100 + tmp_105 + tmp_97) + tmp_121*(tmp_110 + tmp_113 + tmp_118) + tmp_134*(tmp_123 + tmp_126 + tmp_131) + tmp_147*(tmp_136 + tmp_139 + tmp_144) + tmp_160*(tmp_149 + tmp_152 + tmp_157) + tmp_173*(tmp_162 + tmp_165 + tmp_170) + tmp_186*(tmp_175 + tmp_178 + tmp_183) + tmp_199*(tmp_188 + tmp_191 + tmp_196) + tmp_212*(tmp_201 + tmp_204 + tmp_209) + tmp_225*(tmp_214 + tmp_217 + tmp_222) + tmp_238*(tmp_227 + tmp_230 + tmp_235) + tmp_251*(tmp_240 + tmp_243 + tmp_248) + tmp_264*(tmp_253 + tmp_256 + tmp_261) + tmp_277*(tmp_266 + tmp_269 + tmp_274) + tmp_290*(tmp_279 + tmp_282 + tmp_287) + tmp_303*(tmp_292 + tmp_295 + tmp_300) + tmp_316*(tmp_305 + tmp_308 + tmp_313) + tmp_56*(tmp_25 + tmp_34 + tmp_47) + tmp_69*(tmp_58 + tmp_61 + tmp_66) + tmp_82*(tmp_71 + tmp_74 + tmp_79) + tmp_95*(tmp_84 + tmp_87 + tmp_92);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
    const Eigen::Matrix< real_t, 3, 1 >&,
    const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
    Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = p_affine_1_2 + tmp_9;
      real_t tmp_13 = tmp_12*tmp_5;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = p_affine_2_2 + tmp_9;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_1*tmp_11;
      real_t tmp_18 = tmp_12*tmp_14;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_4 - tmp_10*tmp_7 + tmp_11*tmp_13 + tmp_14*tmp_16 - tmp_15*tmp_17 - tmp_18*tmp_3);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_9;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_8;
      real_t tmp_26 = tmp_14*tmp_6 - tmp_17;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_15 + tmp_13;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_10 - tmp_18;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_11*tmp_5 - tmp_14*tmp_3;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_12*tmp_3 + tmp_16;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = tmp_10*tmp_3 - tmp_11*tmp_15;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_53 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_54 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_55 = 1.0*p_affine_13_0*std::pow((std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)*std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)) + (std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)*std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)) + (std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)*std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)), 1.0/2.0);
      real_t tmp_56 = 0.0068572537431980923*tmp_55;
      real_t tmp_57 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_58 = tmp_57*tmp_8;
      real_t tmp_59 = tmp_26*tmp_57;
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_61 = tmp_28*tmp_60;
      real_t tmp_62 = tmp_35*tmp_60;
      real_t tmp_63 = tmp_37*tmp_57;
      real_t tmp_64 = tmp_39*tmp_60;
      real_t tmp_65 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_66 = tmp_41*tmp_65;
      real_t tmp_67 = tmp_48*tmp_65;
      real_t tmp_68 = tmp_50*tmp_65;
      real_t tmp_69 = 0.037198804536718075*tmp_55;
      real_t tmp_70 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_71 = tmp_70*tmp_8;
      real_t tmp_72 = tmp_26*tmp_70;
      real_t tmp_73 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_74 = tmp_28*tmp_73;
      real_t tmp_75 = tmp_35*tmp_73;
      real_t tmp_76 = tmp_37*tmp_70;
      real_t tmp_77 = tmp_39*tmp_73;
      real_t tmp_78 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_79 = tmp_41*tmp_78;
      real_t tmp_80 = tmp_48*tmp_78;
      real_t tmp_81 = tmp_50*tmp_78;
      real_t tmp_82 = 0.020848748529055869*tmp_55;
      real_t tmp_83 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_84 = tmp_8*tmp_83;
      real_t tmp_85 = tmp_26*tmp_83;
      real_t tmp_86 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_87 = tmp_28*tmp_86;
      real_t tmp_88 = tmp_35*tmp_86;
      real_t tmp_89 = tmp_37*tmp_83;
      real_t tmp_90 = tmp_39*tmp_86;
      real_t tmp_91 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_92 = tmp_41*tmp_91;
      real_t tmp_93 = tmp_48*tmp_91;
      real_t tmp_94 = tmp_50*tmp_91;
      real_t tmp_95 = 0.019202922745021479*tmp_55;
      real_t tmp_96 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_97 = tmp_8*tmp_96;
      real_t tmp_98 = tmp_26*tmp_96;
      real_t tmp_99 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_100 = tmp_28*tmp_99;
      real_t tmp_101 = tmp_35*tmp_99;
      real_t tmp_102 = tmp_37*tmp_96;
      real_t tmp_103 = tmp_39*tmp_99;
      real_t tmp_104 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_105 = tmp_104*tmp_41;
      real_t tmp_106 = tmp_104*tmp_48;
      real_t tmp_107 = tmp_104*tmp_50;
      real_t tmp_108 = 0.020848748529055869*tmp_55;
      real_t tmp_109 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_110 = tmp_109*tmp_8;
      real_t tmp_111 = tmp_109*tmp_26;
      real_t tmp_112 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_113 = tmp_112*tmp_28;
      real_t tmp_114 = tmp_112*tmp_35;
      real_t tmp_115 = tmp_109*tmp_37;
      real_t tmp_116 = tmp_112*tmp_39;
      real_t tmp_117 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_118 = tmp_117*tmp_41;
      real_t tmp_119 = tmp_117*tmp_48;
      real_t tmp_120 = tmp_117*tmp_50;
      real_t tmp_121 = 0.019202922745021479*tmp_55;
      real_t tmp_122 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_123 = tmp_122*tmp_8;
      real_t tmp_124 = tmp_122*tmp_26;
      real_t tmp_125 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_126 = tmp_125*tmp_28;
      real_t tmp_127 = tmp_125*tmp_35;
      real_t tmp_128 = tmp_122*tmp_37;
      real_t tmp_129 = tmp_125*tmp_39;
      real_t tmp_130 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_131 = tmp_130*tmp_41;
      real_t tmp_132 = tmp_130*tmp_48;
      real_t tmp_133 = tmp_130*tmp_50;
      real_t tmp_134 = 0.020848748529055869*tmp_55;
      real_t tmp_135 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_136 = tmp_135*tmp_8;
      real_t tmp_137 = tmp_135*tmp_26;
      real_t tmp_138 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_139 = tmp_138*tmp_28;
      real_t tmp_140 = tmp_138*tmp_35;
      real_t tmp_141 = tmp_135*tmp_37;
      real_t tmp_142 = tmp_138*tmp_39;
      real_t tmp_143 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_144 = tmp_143*tmp_41;
      real_t tmp_145 = tmp_143*tmp_48;
      real_t tmp_146 = tmp_143*tmp_50;
      real_t tmp_147 = 0.019202922745021479*tmp_55;
      real_t tmp_148 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_149 = tmp_148*tmp_8;
      real_t tmp_150 = tmp_148*tmp_26;
      real_t tmp_151 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_152 = tmp_151*tmp_28;
      real_t tmp_153 = tmp_151*tmp_35;
      real_t tmp_154 = tmp_148*tmp_37;
      real_t tmp_155 = tmp_151*tmp_39;
      real_t tmp_156 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_157 = tmp_156*tmp_41;
      real_t tmp_158 = tmp_156*tmp_48;
      real_t tmp_159 = tmp_156*tmp_50;
      real_t tmp_160 = 0.020848748529055869*tmp_55;
      real_t tmp_161 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_162 = tmp_161*tmp_8;
      real_t tmp_163 = tmp_161*tmp_26;
      real_t tmp_164 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_165 = tmp_164*tmp_28;
      real_t tmp_166 = tmp_164*tmp_35;
      real_t tmp_167 = tmp_161*tmp_37;
      real_t tmp_168 = tmp_164*tmp_39;
      real_t tmp_169 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_170 = tmp_169*tmp_41;
      real_t tmp_171 = tmp_169*tmp_48;
      real_t tmp_172 = tmp_169*tmp_50;
      real_t tmp_173 = 0.019202922745021479*tmp_55;
      real_t tmp_174 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_175 = tmp_174*tmp_8;
      real_t tmp_176 = tmp_174*tmp_26;
      real_t tmp_177 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_178 = tmp_177*tmp_28;
      real_t tmp_179 = tmp_177*tmp_35;
      real_t tmp_180 = tmp_174*tmp_37;
      real_t tmp_181 = tmp_177*tmp_39;
      real_t tmp_182 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_183 = tmp_182*tmp_41;
      real_t tmp_184 = tmp_182*tmp_48;
      real_t tmp_185 = tmp_182*tmp_50;
      real_t tmp_186 = 0.020848748529055869*tmp_55;
      real_t tmp_187 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_188 = tmp_187*tmp_8;
      real_t tmp_189 = tmp_187*tmp_26;
      real_t tmp_190 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_191 = tmp_190*tmp_28;
      real_t tmp_192 = tmp_190*tmp_35;
      real_t tmp_193 = tmp_187*tmp_37;
      real_t tmp_194 = tmp_190*tmp_39;
      real_t tmp_195 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_196 = tmp_195*tmp_41;
      real_t tmp_197 = tmp_195*tmp_48;
      real_t tmp_198 = tmp_195*tmp_50;
      real_t tmp_199 = 0.019202922745021479*tmp_55;
      real_t tmp_200 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_201 = tmp_200*tmp_8;
      real_t tmp_202 = tmp_200*tmp_26;
      real_t tmp_203 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_204 = tmp_203*tmp_28;
      real_t tmp_205 = tmp_203*tmp_35;
      real_t tmp_206 = tmp_200*tmp_37;
      real_t tmp_207 = tmp_203*tmp_39;
      real_t tmp_208 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_209 = tmp_208*tmp_41;
      real_t tmp_210 = tmp_208*tmp_48;
      real_t tmp_211 = tmp_208*tmp_50;
      real_t tmp_212 = 0.042507265838595799*tmp_55;
      real_t tmp_213 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_214 = tmp_213*tmp_8;
      real_t tmp_215 = tmp_213*tmp_26;
      real_t tmp_216 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_217 = tmp_216*tmp_28;
      real_t tmp_218 = tmp_216*tmp_35;
      real_t tmp_219 = tmp_213*tmp_37;
      real_t tmp_220 = tmp_216*tmp_39;
      real_t tmp_221 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_222 = tmp_221*tmp_41;
      real_t tmp_223 = tmp_221*tmp_48;
      real_t tmp_224 = tmp_221*tmp_50;
      real_t tmp_225 = 0.020848748529055869*tmp_55;
      real_t tmp_226 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_227 = tmp_226*tmp_8;
      real_t tmp_228 = tmp_226*tmp_26;
      real_t tmp_229 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_230 = tmp_229*tmp_28;
      real_t tmp_231 = tmp_229*tmp_35;
      real_t tmp_232 = tmp_226*tmp_37;
      real_t tmp_233 = tmp_229*tmp_39;
      real_t tmp_234 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_235 = tmp_234*tmp_41;
      real_t tmp_236 = tmp_234*tmp_48;
      real_t tmp_237 = tmp_234*tmp_50;
      real_t tmp_238 = 0.0068572537431980923*tmp_55;
      real_t tmp_239 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_240 = tmp_239*tmp_8;
      real_t tmp_241 = tmp_239*tmp_26;
      real_t tmp_242 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_243 = tmp_242*tmp_28;
      real_t tmp_244 = tmp_242*tmp_35;
      real_t tmp_245 = tmp_239*tmp_37;
      real_t tmp_246 = tmp_242*tmp_39;
      real_t tmp_247 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_248 = tmp_247*tmp_41;
      real_t tmp_249 = tmp_247*tmp_48;
      real_t tmp_250 = tmp_247*tmp_50;
      real_t tmp_251 = 0.037198804536718075*tmp_55;
      real_t tmp_252 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_8;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = 0.042507265838595799*tmp_55;
      real_t tmp_265 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_266 = tmp_265*tmp_8;
      real_t tmp_267 = tmp_26*tmp_265;
      real_t tmp_268 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_269 = tmp_268*tmp_28;
      real_t tmp_270 = tmp_268*tmp_35;
      real_t tmp_271 = tmp_265*tmp_37;
      real_t tmp_272 = tmp_268*tmp_39;
      real_t tmp_273 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_274 = tmp_273*tmp_41;
      real_t tmp_275 = tmp_273*tmp_48;
      real_t tmp_276 = tmp_273*tmp_50;
      real_t tmp_277 = 0.0068572537431980923*tmp_55;
      real_t tmp_278 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_279 = tmp_278*tmp_8;
      real_t tmp_280 = tmp_26*tmp_278;
      real_t tmp_281 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_282 = tmp_28*tmp_281;
      real_t tmp_283 = tmp_281*tmp_35;
      real_t tmp_284 = tmp_278*tmp_37;
      real_t tmp_285 = tmp_281*tmp_39;
      real_t tmp_286 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_287 = tmp_286*tmp_41;
      real_t tmp_288 = tmp_286*tmp_48;
      real_t tmp_289 = tmp_286*tmp_50;
      real_t tmp_290 = 0.037198804536718075*tmp_55;
      real_t tmp_291 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_292 = tmp_291*tmp_8;
      real_t tmp_293 = tmp_26*tmp_291;
      real_t tmp_294 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_295 = tmp_28*tmp_294;
      real_t tmp_296 = tmp_294*tmp_35;
      real_t tmp_297 = tmp_291*tmp_37;
      real_t tmp_298 = tmp_294*tmp_39;
      real_t tmp_299 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_300 = tmp_299*tmp_41;
      real_t tmp_301 = tmp_299*tmp_48;
      real_t tmp_302 = tmp_299*tmp_50;
      real_t tmp_303 = 0.042507265838595799*tmp_55;
      real_t tmp_304 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_305 = tmp_304*tmp_8;
      real_t tmp_306 = tmp_26*tmp_304;
      real_t tmp_307 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_308 = tmp_28*tmp_307;
      real_t tmp_309 = tmp_307*tmp_35;
      real_t tmp_310 = tmp_304*tmp_37;
      real_t tmp_311 = tmp_307*tmp_39;
      real_t tmp_312 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_313 = tmp_312*tmp_41;
      real_t tmp_314 = tmp_312*tmp_48;
      real_t tmp_315 = tmp_312*tmp_50;
      real_t tmp_316 = 0.019202922745021479*tmp_55;
      real_t a_0_0 = tmp_108*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_105 - tmp_106 - tmp_107 - tmp_97 - tmp_98 + 1) + tmp_121*(-tmp_110 - tmp_111 - tmp_113 - tmp_114 - tmp_115 - tmp_116 - tmp_118 - tmp_119 - tmp_120 + 1) + tmp_134*(-tmp_123 - tmp_124 - tmp_126 - tmp_127 - tmp_128 - tmp_129 - tmp_131 - tmp_132 - tmp_133 + 1) + tmp_147*(-tmp_136 - tmp_137 - tmp_139 - tmp_140 - tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 + 1) + tmp_160*(-tmp_149 - tmp_150 - tmp_152 - tmp_153 - tmp_154 - tmp_155 - tmp_157 - tmp_158 - tmp_159 + 1) + tmp_173*(-tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 - tmp_168 - tmp_170 - tmp_171 - tmp_172 + 1) + tmp_186*(-tmp_175 - tmp_176 - tmp_178 - tmp_179 - tmp_180 - tmp_181 - tmp_183 - tmp_184 - tmp_185 + 1) + tmp_199*(-tmp_188 - tmp_189 - tmp_191 - tmp_192 - tmp_193 - tmp_194 - tmp_196 - tmp_197 - tmp_198 + 1) + tmp_212*(-tmp_201 - tmp_202 - tmp_204 - tmp_205 - tmp_206 - tmp_207 - tmp_209 - tmp_210 - tmp_211 + 1) + tmp_225*(-tmp_214 - tmp_215 - tmp_217 - tmp_218 - tmp_219 - tmp_220 - tmp_222 - tmp_223 - tmp_224 + 1) + tmp_238*(-tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 - tmp_233 - tmp_235 - tmp_236 - tmp_237 + 1) + tmp_251*(-tmp_240 - tmp_241 - tmp_243 - tmp_244 - tmp_245 - tmp_246 - tmp_248 - tmp_249 - tmp_250 + 1) + tmp_264*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_277*(-tmp_266 - tmp_267 - tmp_269 - tmp_270 - tmp_271 - tmp_272 - tmp_274 - tmp_275 - tmp_276 + 1) + tmp_290*(-tmp_279 - tmp_280 - tmp_282 - tmp_283 - tmp_284 - tmp_285 - tmp_287 - tmp_288 - tmp_289 + 1) + tmp_303*(-tmp_292 - tmp_293 - tmp_295 - tmp_296 - tmp_297 - tmp_298 - tmp_300 - tmp_301 - tmp_302 + 1) + tmp_316*(-tmp_305 - tmp_306 - tmp_308 - tmp_309 - tmp_310 - tmp_311 - tmp_313 - tmp_314 - tmp_315 + 1) + tmp_56*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_69*(-tmp_58 - tmp_59 - tmp_61 - tmp_62 - tmp_63 - tmp_64 - tmp_66 - tmp_67 - tmp_68 + 1) + tmp_82*(-tmp_71 - tmp_72 - tmp_74 - tmp_75 - tmp_76 - tmp_77 - tmp_79 - tmp_80 - tmp_81 + 1) + tmp_95*(-tmp_84 - tmp_85 - tmp_87 - tmp_88 - tmp_89 - tmp_90 - tmp_92 - tmp_93 - tmp_94 + 1);
      real_t a_1_0 = tmp_108*(tmp_102 + tmp_103 + tmp_107) + tmp_121*(tmp_115 + tmp_116 + tmp_120) + tmp_134*(tmp_128 + tmp_129 + tmp_133) + tmp_147*(tmp_141 + tmp_142 + tmp_146) + tmp_160*(tmp_154 + tmp_155 + tmp_159) + tmp_173*(tmp_167 + tmp_168 + tmp_172) + tmp_186*(tmp_180 + tmp_181 + tmp_185) + tmp_199*(tmp_193 + tmp_194 + tmp_198) + tmp_212*(tmp_206 + tmp_207 + tmp_211) + tmp_225*(tmp_219 + tmp_220 + tmp_224) + tmp_238*(tmp_232 + tmp_233 + tmp_237) + tmp_251*(tmp_245 + tmp_246 + tmp_250) + tmp_264*(tmp_258 + tmp_259 + tmp_263) + tmp_277*(tmp_271 + tmp_272 + tmp_276) + tmp_290*(tmp_284 + tmp_285 + tmp_289) + tmp_303*(tmp_297 + tmp_298 + tmp_302) + tmp_316*(tmp_310 + tmp_311 + tmp_315) + tmp_56*(tmp_38 + tmp_40 + tmp_51) + tmp_69*(tmp_63 + tmp_64 + tmp_68) + tmp_82*(tmp_76 + tmp_77 + tmp_81) + tmp_95*(tmp_89 + tmp_90 + tmp_94);
      real_t a_2_0 = tmp_108*(tmp_101 + tmp_106 + tmp_98) + tmp_121*(tmp_111 + tmp_114 + tmp_119) + tmp_134*(tmp_124 + tmp_127 + tmp_132) + tmp_147*(tmp_137 + tmp_140 + tmp_145) + tmp_160*(tmp_150 + tmp_153 + tmp_158) + tmp_173*(tmp_163 + tmp_166 + tmp_171) + tmp_186*(tmp_176 + tmp_179 + tmp_184) + tmp_199*(tmp_189 + tmp_192 + tmp_197) + tmp_212*(tmp_202 + tmp_205 + tmp_210) + tmp_225*(tmp_215 + tmp_218 + tmp_223) + tmp_238*(tmp_228 + tmp_231 + tmp_236) + tmp_251*(tmp_241 + tmp_244 + tmp_249) + tmp_264*(tmp_254 + tmp_257 + tmp_262) + tmp_277*(tmp_267 + tmp_270 + tmp_275) + tmp_290*(tmp_280 + tmp_283 + tmp_288) + tmp_303*(tmp_293 + tmp_296 + tmp_301) + tmp_316*(tmp_306 + tmp_309 + tmp_314) + tmp_56*(tmp_27 + tmp_36 + tmp_49) + tmp_69*(tmp_59 + tmp_62 + tmp_67) + tmp_82*(tmp_72 + tmp_75 + tmp_80) + tmp_95*(tmp_85 + tmp_88 + tmp_93);
      real_t a_3_0 = tmp_108*(tmp_100 + tmp_105 + tmp_97) + tmp_121*(tmp_110 + tmp_113 + tmp_118) + tmp_134*(tmp_123 + tmp_126 + tmp_131) + tmp_147*(tmp_136 + tmp_139 + tmp_144) + tmp_160*(tmp_149 + tmp_152 + tmp_157) + tmp_173*(tmp_162 + tmp_165 + tmp_170) + tmp_186*(tmp_175 + tmp_178 + tmp_183) + tmp_199*(tmp_188 + tmp_191 + tmp_196) + tmp_212*(tmp_201 + tmp_204 + tmp_209) + tmp_225*(tmp_214 + tmp_217 + tmp_222) + tmp_238*(tmp_227 + tmp_230 + tmp_235) + tmp_251*(tmp_240 + tmp_243 + tmp_248) + tmp_264*(tmp_253 + tmp_256 + tmp_261) + tmp_277*(tmp_266 + tmp_269 + tmp_274) + tmp_290*(tmp_279 + tmp_282 + tmp_287) + tmp_303*(tmp_292 + tmp_295 + tmp_300) + tmp_316*(tmp_305 + tmp_308 + tmp_313) + tmp_56*(tmp_25 + tmp_34 + tmp_47) + tmp_69*(tmp_58 + tmp_61 + tmp_66) + tmp_82*(tmp_71 + tmp_74 + tmp_79) + tmp_95*(tmp_84 + tmp_87 + tmp_92);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }

public:



};




class EGDivtForm_P1P0_1 : public hyteg::dg::DGForm
{

 public:
    EGDivtForm_P1P0_1()

    {}





 protected:
  void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      elMat( 0, 0) = 0;
      elMat( 1, 0) = 0;
      elMat( 2, 0) = 0;
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
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

      elMat( 0, 0) = 0;
      elMat( 1, 0) = 0;
      elMat( 2, 0) = 0;
      elMat( 3, 0) = 0;
   }
   void integrateVolume3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_2_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_2;
      real_t tmp_3 = p_affine_1_2 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_1_0 + tmp_0;
      real_t tmp_6 = p_affine_2_2 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = -p_affine_0_1;
      real_t tmp_9 = p_affine_2_1 + tmp_8;
      real_t tmp_10 = p_affine_3_2 + tmp_2;
      real_t tmp_11 = tmp_10*tmp_5;
      real_t tmp_12 = p_affine_3_1 + tmp_8;
      real_t tmp_13 = p_affine_1_1 + tmp_8;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = tmp_14*tmp_6;
      real_t tmp_16 = tmp_1*tmp_10;
      real_t tmp_17 = tmp_14*tmp_3;
      real_t tmp_18 = 1.0 / (tmp_11*tmp_9 + tmp_12*tmp_4 - tmp_12*tmp_7 + tmp_13*tmp_15 - tmp_13*tmp_16 - tmp_17*tmp_9);
      real_t tmp_19 = tmp_18*(tmp_4 - tmp_7);
      real_t tmp_20 = tmp_18*(tmp_11 - tmp_17);
      real_t tmp_21 = tmp_18*(tmp_15 - tmp_16);
      real_t tmp_22 = p_affine_0_0*p_affine_1_1;
      real_t tmp_23 = p_affine_0_0*p_affine_1_2;
      real_t tmp_24 = p_affine_2_1*p_affine_3_2;
      real_t tmp_25 = p_affine_0_1*p_affine_1_0;
      real_t tmp_26 = p_affine_0_1*p_affine_1_2;
      real_t tmp_27 = p_affine_2_2*p_affine_3_0;
      real_t tmp_28 = p_affine_0_2*p_affine_1_0;
      real_t tmp_29 = p_affine_0_2*p_affine_1_1;
      real_t tmp_30 = p_affine_2_0*p_affine_3_1;
      real_t tmp_31 = p_affine_2_2*p_affine_3_1;
      real_t tmp_32 = p_affine_2_0*p_affine_3_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_0;
      real_t tmp_34 = std::abs(p_affine_0_0*tmp_24 - p_affine_0_0*tmp_31 + p_affine_0_1*tmp_27 - p_affine_0_1*tmp_32 + p_affine_0_2*tmp_30 - p_affine_0_2*tmp_33 - p_affine_1_0*tmp_24 + p_affine_1_0*tmp_31 - p_affine_1_1*tmp_27 + p_affine_1_1*tmp_32 - p_affine_1_2*tmp_30 + p_affine_1_2*tmp_33 + p_affine_2_0*tmp_26 - p_affine_2_0*tmp_29 - p_affine_2_1*tmp_23 + p_affine_2_1*tmp_28 + p_affine_2_2*tmp_22 - p_affine_2_2*tmp_25 - p_affine_3_0*tmp_26 + p_affine_3_0*tmp_29 + p_affine_3_1*tmp_23 - p_affine_3_1*tmp_28 - p_affine_3_2*tmp_22 + p_affine_3_2*tmp_25);
      real_t tmp_35 = tmp_34*(tmp_19 + tmp_20 + tmp_21);
      real_t tmp_36 = tmp_21*tmp_34;
      real_t tmp_37 = tmp_20*tmp_34;
      real_t tmp_38 = tmp_19*tmp_34;
      real_t a_0_0 = 0.1666666666666668*tmp_35;
      real_t a_1_0 = -0.1666666666666668*tmp_36;
      real_t a_2_0 = -0.1666666666666668*tmp_37;
      real_t a_3_0 = -0.1666666666666668*tmp_38;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }



   void integrateFacetInner3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
                                                     const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                     const Eigen::Matrix< real_t, 3, 1 >&,
                                                     const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                                                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

         real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = p_affine_1_2 + tmp_9;
      real_t tmp_13 = tmp_12*tmp_5;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = p_affine_2_2 + tmp_9;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_1*tmp_11;
      real_t tmp_18 = tmp_12*tmp_14;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_4 - tmp_10*tmp_7 + tmp_11*tmp_13 + tmp_14*tmp_16 - tmp_15*tmp_17 - tmp_18*tmp_3);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_9;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_8;
      real_t tmp_26 = tmp_14*tmp_6 - tmp_17;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_15 + tmp_13;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_10 - tmp_18;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_11*tmp_5 - tmp_14*tmp_3;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_12*tmp_3 + tmp_16;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = tmp_10*tmp_3 - tmp_11*tmp_15;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_53 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_54 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_55 = 0.5*p_affine_13_1*std::pow((std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)*std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)) + (std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)*std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)) + (std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)*std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)), 1.0/2.0);
      real_t tmp_56 = 0.0068572537431980923*tmp_55;
      real_t tmp_57 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_58 = tmp_57*tmp_8;
      real_t tmp_59 = tmp_26*tmp_57;
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_61 = tmp_28*tmp_60;
      real_t tmp_62 = tmp_35*tmp_60;
      real_t tmp_63 = tmp_37*tmp_57;
      real_t tmp_64 = tmp_39*tmp_60;
      real_t tmp_65 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_66 = tmp_41*tmp_65;
      real_t tmp_67 = tmp_48*tmp_65;
      real_t tmp_68 = tmp_50*tmp_65;
      real_t tmp_69 = 0.037198804536718075*tmp_55;
      real_t tmp_70 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_71 = tmp_70*tmp_8;
      real_t tmp_72 = tmp_26*tmp_70;
      real_t tmp_73 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_74 = tmp_28*tmp_73;
      real_t tmp_75 = tmp_35*tmp_73;
      real_t tmp_76 = tmp_37*tmp_70;
      real_t tmp_77 = tmp_39*tmp_73;
      real_t tmp_78 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_79 = tmp_41*tmp_78;
      real_t tmp_80 = tmp_48*tmp_78;
      real_t tmp_81 = tmp_50*tmp_78;
      real_t tmp_82 = 0.020848748529055869*tmp_55;
      real_t tmp_83 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_84 = tmp_8*tmp_83;
      real_t tmp_85 = tmp_26*tmp_83;
      real_t tmp_86 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_87 = tmp_28*tmp_86;
      real_t tmp_88 = tmp_35*tmp_86;
      real_t tmp_89 = tmp_37*tmp_83;
      real_t tmp_90 = tmp_39*tmp_86;
      real_t tmp_91 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_92 = tmp_41*tmp_91;
      real_t tmp_93 = tmp_48*tmp_91;
      real_t tmp_94 = tmp_50*tmp_91;
      real_t tmp_95 = 0.019202922745021479*tmp_55;
      real_t tmp_96 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_97 = tmp_8*tmp_96;
      real_t tmp_98 = tmp_26*tmp_96;
      real_t tmp_99 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_100 = tmp_28*tmp_99;
      real_t tmp_101 = tmp_35*tmp_99;
      real_t tmp_102 = tmp_37*tmp_96;
      real_t tmp_103 = tmp_39*tmp_99;
      real_t tmp_104 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_105 = tmp_104*tmp_41;
      real_t tmp_106 = tmp_104*tmp_48;
      real_t tmp_107 = tmp_104*tmp_50;
      real_t tmp_108 = 0.020848748529055869*tmp_55;
      real_t tmp_109 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_110 = tmp_109*tmp_8;
      real_t tmp_111 = tmp_109*tmp_26;
      real_t tmp_112 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_113 = tmp_112*tmp_28;
      real_t tmp_114 = tmp_112*tmp_35;
      real_t tmp_115 = tmp_109*tmp_37;
      real_t tmp_116 = tmp_112*tmp_39;
      real_t tmp_117 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_118 = tmp_117*tmp_41;
      real_t tmp_119 = tmp_117*tmp_48;
      real_t tmp_120 = tmp_117*tmp_50;
      real_t tmp_121 = 0.019202922745021479*tmp_55;
      real_t tmp_122 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_123 = tmp_122*tmp_8;
      real_t tmp_124 = tmp_122*tmp_26;
      real_t tmp_125 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_126 = tmp_125*tmp_28;
      real_t tmp_127 = tmp_125*tmp_35;
      real_t tmp_128 = tmp_122*tmp_37;
      real_t tmp_129 = tmp_125*tmp_39;
      real_t tmp_130 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_131 = tmp_130*tmp_41;
      real_t tmp_132 = tmp_130*tmp_48;
      real_t tmp_133 = tmp_130*tmp_50;
      real_t tmp_134 = 0.020848748529055869*tmp_55;
      real_t tmp_135 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_136 = tmp_135*tmp_8;
      real_t tmp_137 = tmp_135*tmp_26;
      real_t tmp_138 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_139 = tmp_138*tmp_28;
      real_t tmp_140 = tmp_138*tmp_35;
      real_t tmp_141 = tmp_135*tmp_37;
      real_t tmp_142 = tmp_138*tmp_39;
      real_t tmp_143 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_144 = tmp_143*tmp_41;
      real_t tmp_145 = tmp_143*tmp_48;
      real_t tmp_146 = tmp_143*tmp_50;
      real_t tmp_147 = 0.019202922745021479*tmp_55;
      real_t tmp_148 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_149 = tmp_148*tmp_8;
      real_t tmp_150 = tmp_148*tmp_26;
      real_t tmp_151 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_152 = tmp_151*tmp_28;
      real_t tmp_153 = tmp_151*tmp_35;
      real_t tmp_154 = tmp_148*tmp_37;
      real_t tmp_155 = tmp_151*tmp_39;
      real_t tmp_156 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_157 = tmp_156*tmp_41;
      real_t tmp_158 = tmp_156*tmp_48;
      real_t tmp_159 = tmp_156*tmp_50;
      real_t tmp_160 = 0.020848748529055869*tmp_55;
      real_t tmp_161 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_162 = tmp_161*tmp_8;
      real_t tmp_163 = tmp_161*tmp_26;
      real_t tmp_164 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_165 = tmp_164*tmp_28;
      real_t tmp_166 = tmp_164*tmp_35;
      real_t tmp_167 = tmp_161*tmp_37;
      real_t tmp_168 = tmp_164*tmp_39;
      real_t tmp_169 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_170 = tmp_169*tmp_41;
      real_t tmp_171 = tmp_169*tmp_48;
      real_t tmp_172 = tmp_169*tmp_50;
      real_t tmp_173 = 0.019202922745021479*tmp_55;
      real_t tmp_174 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_175 = tmp_174*tmp_8;
      real_t tmp_176 = tmp_174*tmp_26;
      real_t tmp_177 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_178 = tmp_177*tmp_28;
      real_t tmp_179 = tmp_177*tmp_35;
      real_t tmp_180 = tmp_174*tmp_37;
      real_t tmp_181 = tmp_177*tmp_39;
      real_t tmp_182 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_183 = tmp_182*tmp_41;
      real_t tmp_184 = tmp_182*tmp_48;
      real_t tmp_185 = tmp_182*tmp_50;
      real_t tmp_186 = 0.020848748529055869*tmp_55;
      real_t tmp_187 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_188 = tmp_187*tmp_8;
      real_t tmp_189 = tmp_187*tmp_26;
      real_t tmp_190 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_191 = tmp_190*tmp_28;
      real_t tmp_192 = tmp_190*tmp_35;
      real_t tmp_193 = tmp_187*tmp_37;
      real_t tmp_194 = tmp_190*tmp_39;
      real_t tmp_195 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_196 = tmp_195*tmp_41;
      real_t tmp_197 = tmp_195*tmp_48;
      real_t tmp_198 = tmp_195*tmp_50;
      real_t tmp_199 = 0.019202922745021479*tmp_55;
      real_t tmp_200 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_201 = tmp_200*tmp_8;
      real_t tmp_202 = tmp_200*tmp_26;
      real_t tmp_203 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_204 = tmp_203*tmp_28;
      real_t tmp_205 = tmp_203*tmp_35;
      real_t tmp_206 = tmp_200*tmp_37;
      real_t tmp_207 = tmp_203*tmp_39;
      real_t tmp_208 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_209 = tmp_208*tmp_41;
      real_t tmp_210 = tmp_208*tmp_48;
      real_t tmp_211 = tmp_208*tmp_50;
      real_t tmp_212 = 0.042507265838595799*tmp_55;
      real_t tmp_213 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_214 = tmp_213*tmp_8;
      real_t tmp_215 = tmp_213*tmp_26;
      real_t tmp_216 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_217 = tmp_216*tmp_28;
      real_t tmp_218 = tmp_216*tmp_35;
      real_t tmp_219 = tmp_213*tmp_37;
      real_t tmp_220 = tmp_216*tmp_39;
      real_t tmp_221 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_222 = tmp_221*tmp_41;
      real_t tmp_223 = tmp_221*tmp_48;
      real_t tmp_224 = tmp_221*tmp_50;
      real_t tmp_225 = 0.020848748529055869*tmp_55;
      real_t tmp_226 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_227 = tmp_226*tmp_8;
      real_t tmp_228 = tmp_226*tmp_26;
      real_t tmp_229 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_230 = tmp_229*tmp_28;
      real_t tmp_231 = tmp_229*tmp_35;
      real_t tmp_232 = tmp_226*tmp_37;
      real_t tmp_233 = tmp_229*tmp_39;
      real_t tmp_234 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_235 = tmp_234*tmp_41;
      real_t tmp_236 = tmp_234*tmp_48;
      real_t tmp_237 = tmp_234*tmp_50;
      real_t tmp_238 = 0.0068572537431980923*tmp_55;
      real_t tmp_239 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_240 = tmp_239*tmp_8;
      real_t tmp_241 = tmp_239*tmp_26;
      real_t tmp_242 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_243 = tmp_242*tmp_28;
      real_t tmp_244 = tmp_242*tmp_35;
      real_t tmp_245 = tmp_239*tmp_37;
      real_t tmp_246 = tmp_242*tmp_39;
      real_t tmp_247 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_248 = tmp_247*tmp_41;
      real_t tmp_249 = tmp_247*tmp_48;
      real_t tmp_250 = tmp_247*tmp_50;
      real_t tmp_251 = 0.037198804536718075*tmp_55;
      real_t tmp_252 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_8;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = 0.042507265838595799*tmp_55;
      real_t tmp_265 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_266 = tmp_265*tmp_8;
      real_t tmp_267 = tmp_26*tmp_265;
      real_t tmp_268 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_269 = tmp_268*tmp_28;
      real_t tmp_270 = tmp_268*tmp_35;
      real_t tmp_271 = tmp_265*tmp_37;
      real_t tmp_272 = tmp_268*tmp_39;
      real_t tmp_273 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_274 = tmp_273*tmp_41;
      real_t tmp_275 = tmp_273*tmp_48;
      real_t tmp_276 = tmp_273*tmp_50;
      real_t tmp_277 = 0.0068572537431980923*tmp_55;
      real_t tmp_278 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_279 = tmp_278*tmp_8;
      real_t tmp_280 = tmp_26*tmp_278;
      real_t tmp_281 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_282 = tmp_28*tmp_281;
      real_t tmp_283 = tmp_281*tmp_35;
      real_t tmp_284 = tmp_278*tmp_37;
      real_t tmp_285 = tmp_281*tmp_39;
      real_t tmp_286 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_287 = tmp_286*tmp_41;
      real_t tmp_288 = tmp_286*tmp_48;
      real_t tmp_289 = tmp_286*tmp_50;
      real_t tmp_290 = 0.037198804536718075*tmp_55;
      real_t tmp_291 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_292 = tmp_291*tmp_8;
      real_t tmp_293 = tmp_26*tmp_291;
      real_t tmp_294 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_295 = tmp_28*tmp_294;
      real_t tmp_296 = tmp_294*tmp_35;
      real_t tmp_297 = tmp_291*tmp_37;
      real_t tmp_298 = tmp_294*tmp_39;
      real_t tmp_299 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_300 = tmp_299*tmp_41;
      real_t tmp_301 = tmp_299*tmp_48;
      real_t tmp_302 = tmp_299*tmp_50;
      real_t tmp_303 = 0.042507265838595799*tmp_55;
      real_t tmp_304 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_305 = tmp_304*tmp_8;
      real_t tmp_306 = tmp_26*tmp_304;
      real_t tmp_307 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_308 = tmp_28*tmp_307;
      real_t tmp_309 = tmp_307*tmp_35;
      real_t tmp_310 = tmp_304*tmp_37;
      real_t tmp_311 = tmp_307*tmp_39;
      real_t tmp_312 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_313 = tmp_312*tmp_41;
      real_t tmp_314 = tmp_312*tmp_48;
      real_t tmp_315 = tmp_312*tmp_50;
      real_t tmp_316 = 0.019202922745021479*tmp_55;
      real_t a_0_0 = tmp_108*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_105 - tmp_106 - tmp_107 - tmp_97 - tmp_98 + 1) + tmp_121*(-tmp_110 - tmp_111 - tmp_113 - tmp_114 - tmp_115 - tmp_116 - tmp_118 - tmp_119 - tmp_120 + 1) + tmp_134*(-tmp_123 - tmp_124 - tmp_126 - tmp_127 - tmp_128 - tmp_129 - tmp_131 - tmp_132 - tmp_133 + 1) + tmp_147*(-tmp_136 - tmp_137 - tmp_139 - tmp_140 - tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 + 1) + tmp_160*(-tmp_149 - tmp_150 - tmp_152 - tmp_153 - tmp_154 - tmp_155 - tmp_157 - tmp_158 - tmp_159 + 1) + tmp_173*(-tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 - tmp_168 - tmp_170 - tmp_171 - tmp_172 + 1) + tmp_186*(-tmp_175 - tmp_176 - tmp_178 - tmp_179 - tmp_180 - tmp_181 - tmp_183 - tmp_184 - tmp_185 + 1) + tmp_199*(-tmp_188 - tmp_189 - tmp_191 - tmp_192 - tmp_193 - tmp_194 - tmp_196 - tmp_197 - tmp_198 + 1) + tmp_212*(-tmp_201 - tmp_202 - tmp_204 - tmp_205 - tmp_206 - tmp_207 - tmp_209 - tmp_210 - tmp_211 + 1) + tmp_225*(-tmp_214 - tmp_215 - tmp_217 - tmp_218 - tmp_219 - tmp_220 - tmp_222 - tmp_223 - tmp_224 + 1) + tmp_238*(-tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 - tmp_233 - tmp_235 - tmp_236 - tmp_237 + 1) + tmp_251*(-tmp_240 - tmp_241 - tmp_243 - tmp_244 - tmp_245 - tmp_246 - tmp_248 - tmp_249 - tmp_250 + 1) + tmp_264*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_277*(-tmp_266 - tmp_267 - tmp_269 - tmp_270 - tmp_271 - tmp_272 - tmp_274 - tmp_275 - tmp_276 + 1) + tmp_290*(-tmp_279 - tmp_280 - tmp_282 - tmp_283 - tmp_284 - tmp_285 - tmp_287 - tmp_288 - tmp_289 + 1) + tmp_303*(-tmp_292 - tmp_293 - tmp_295 - tmp_296 - tmp_297 - tmp_298 - tmp_300 - tmp_301 - tmp_302 + 1) + tmp_316*(-tmp_305 - tmp_306 - tmp_308 - tmp_309 - tmp_310 - tmp_311 - tmp_313 - tmp_314 - tmp_315 + 1) + tmp_56*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_69*(-tmp_58 - tmp_59 - tmp_61 - tmp_62 - tmp_63 - tmp_64 - tmp_66 - tmp_67 - tmp_68 + 1) + tmp_82*(-tmp_71 - tmp_72 - tmp_74 - tmp_75 - tmp_76 - tmp_77 - tmp_79 - tmp_80 - tmp_81 + 1) + tmp_95*(-tmp_84 - tmp_85 - tmp_87 - tmp_88 - tmp_89 - tmp_90 - tmp_92 - tmp_93 - tmp_94 + 1);
      real_t a_1_0 = tmp_108*(tmp_102 + tmp_103 + tmp_107) + tmp_121*(tmp_115 + tmp_116 + tmp_120) + tmp_134*(tmp_128 + tmp_129 + tmp_133) + tmp_147*(tmp_141 + tmp_142 + tmp_146) + tmp_160*(tmp_154 + tmp_155 + tmp_159) + tmp_173*(tmp_167 + tmp_168 + tmp_172) + tmp_186*(tmp_180 + tmp_181 + tmp_185) + tmp_199*(tmp_193 + tmp_194 + tmp_198) + tmp_212*(tmp_206 + tmp_207 + tmp_211) + tmp_225*(tmp_219 + tmp_220 + tmp_224) + tmp_238*(tmp_232 + tmp_233 + tmp_237) + tmp_251*(tmp_245 + tmp_246 + tmp_250) + tmp_264*(tmp_258 + tmp_259 + tmp_263) + tmp_277*(tmp_271 + tmp_272 + tmp_276) + tmp_290*(tmp_284 + tmp_285 + tmp_289) + tmp_303*(tmp_297 + tmp_298 + tmp_302) + tmp_316*(tmp_310 + tmp_311 + tmp_315) + tmp_56*(tmp_38 + tmp_40 + tmp_51) + tmp_69*(tmp_63 + tmp_64 + tmp_68) + tmp_82*(tmp_76 + tmp_77 + tmp_81) + tmp_95*(tmp_89 + tmp_90 + tmp_94);
      real_t a_2_0 = tmp_108*(tmp_101 + tmp_106 + tmp_98) + tmp_121*(tmp_111 + tmp_114 + tmp_119) + tmp_134*(tmp_124 + tmp_127 + tmp_132) + tmp_147*(tmp_137 + tmp_140 + tmp_145) + tmp_160*(tmp_150 + tmp_153 + tmp_158) + tmp_173*(tmp_163 + tmp_166 + tmp_171) + tmp_186*(tmp_176 + tmp_179 + tmp_184) + tmp_199*(tmp_189 + tmp_192 + tmp_197) + tmp_212*(tmp_202 + tmp_205 + tmp_210) + tmp_225*(tmp_215 + tmp_218 + tmp_223) + tmp_238*(tmp_228 + tmp_231 + tmp_236) + tmp_251*(tmp_241 + tmp_244 + tmp_249) + tmp_264*(tmp_254 + tmp_257 + tmp_262) + tmp_277*(tmp_267 + tmp_270 + tmp_275) + tmp_290*(tmp_280 + tmp_283 + tmp_288) + tmp_303*(tmp_293 + tmp_296 + tmp_301) + tmp_316*(tmp_306 + tmp_309 + tmp_314) + tmp_56*(tmp_27 + tmp_36 + tmp_49) + tmp_69*(tmp_59 + tmp_62 + tmp_67) + tmp_82*(tmp_72 + tmp_75 + tmp_80) + tmp_95*(tmp_85 + tmp_88 + tmp_93);
      real_t a_3_0 = tmp_108*(tmp_100 + tmp_105 + tmp_97) + tmp_121*(tmp_110 + tmp_113 + tmp_118) + tmp_134*(tmp_123 + tmp_126 + tmp_131) + tmp_147*(tmp_136 + tmp_139 + tmp_144) + tmp_160*(tmp_149 + tmp_152 + tmp_157) + tmp_173*(tmp_162 + tmp_165 + tmp_170) + tmp_186*(tmp_175 + tmp_178 + tmp_183) + tmp_199*(tmp_188 + tmp_191 + tmp_196) + tmp_212*(tmp_201 + tmp_204 + tmp_209) + tmp_225*(tmp_214 + tmp_217 + tmp_222) + tmp_238*(tmp_227 + tmp_230 + tmp_235) + tmp_251*(tmp_240 + tmp_243 + tmp_248) + tmp_264*(tmp_253 + tmp_256 + tmp_261) + tmp_277*(tmp_266 + tmp_269 + tmp_274) + tmp_290*(tmp_279 + tmp_282 + tmp_287) + tmp_303*(tmp_292 + tmp_295 + tmp_300) + tmp_316*(tmp_305 + tmp_308 + tmp_313) + tmp_56*(tmp_25 + tmp_34 + tmp_47) + tmp_69*(tmp_58 + tmp_61 + tmp_66) + tmp_82*(tmp_71 + tmp_74 + tmp_79) + tmp_95*(tmp_84 + tmp_87 + tmp_92);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }




void integrateFacetCoupling3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementInner,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementOuter,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                                        Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = p_affine_1_2 + tmp_9;
      real_t tmp_13 = tmp_12*tmp_5;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = p_affine_2_2 + tmp_9;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_1*tmp_11;
      real_t tmp_18 = tmp_12*tmp_14;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_4 - tmp_10*tmp_7 + tmp_11*tmp_13 + tmp_14*tmp_16 - tmp_15*tmp_17 - tmp_18*tmp_3);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_9;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_8;
      real_t tmp_26 = tmp_14*tmp_6 - tmp_17;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_15 + tmp_13;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_10 - tmp_18;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_11*tmp_5 - tmp_14*tmp_3;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_12*tmp_3 + tmp_16;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = tmp_10*tmp_3 - tmp_11*tmp_15;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_53 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_54 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_55 = 0.5*p_affine_13_1*std::pow((std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)*std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)) + (std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)*std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)) + (std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)*std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)), 1.0/2.0);
      real_t tmp_56 = 0.0068572537431980923*tmp_55;
      real_t tmp_57 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_58 = tmp_57*tmp_8;
      real_t tmp_59 = tmp_26*tmp_57;
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_61 = tmp_28*tmp_60;
      real_t tmp_62 = tmp_35*tmp_60;
      real_t tmp_63 = tmp_37*tmp_57;
      real_t tmp_64 = tmp_39*tmp_60;
      real_t tmp_65 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_66 = tmp_41*tmp_65;
      real_t tmp_67 = tmp_48*tmp_65;
      real_t tmp_68 = tmp_50*tmp_65;
      real_t tmp_69 = 0.037198804536718075*tmp_55;
      real_t tmp_70 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_71 = tmp_70*tmp_8;
      real_t tmp_72 = tmp_26*tmp_70;
      real_t tmp_73 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_74 = tmp_28*tmp_73;
      real_t tmp_75 = tmp_35*tmp_73;
      real_t tmp_76 = tmp_37*tmp_70;
      real_t tmp_77 = tmp_39*tmp_73;
      real_t tmp_78 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_79 = tmp_41*tmp_78;
      real_t tmp_80 = tmp_48*tmp_78;
      real_t tmp_81 = tmp_50*tmp_78;
      real_t tmp_82 = 0.020848748529055869*tmp_55;
      real_t tmp_83 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_84 = tmp_8*tmp_83;
      real_t tmp_85 = tmp_26*tmp_83;
      real_t tmp_86 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_87 = tmp_28*tmp_86;
      real_t tmp_88 = tmp_35*tmp_86;
      real_t tmp_89 = tmp_37*tmp_83;
      real_t tmp_90 = tmp_39*tmp_86;
      real_t tmp_91 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_92 = tmp_41*tmp_91;
      real_t tmp_93 = tmp_48*tmp_91;
      real_t tmp_94 = tmp_50*tmp_91;
      real_t tmp_95 = 0.019202922745021479*tmp_55;
      real_t tmp_96 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_97 = tmp_8*tmp_96;
      real_t tmp_98 = tmp_26*tmp_96;
      real_t tmp_99 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_100 = tmp_28*tmp_99;
      real_t tmp_101 = tmp_35*tmp_99;
      real_t tmp_102 = tmp_37*tmp_96;
      real_t tmp_103 = tmp_39*tmp_99;
      real_t tmp_104 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_105 = tmp_104*tmp_41;
      real_t tmp_106 = tmp_104*tmp_48;
      real_t tmp_107 = tmp_104*tmp_50;
      real_t tmp_108 = 0.020848748529055869*tmp_55;
      real_t tmp_109 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_110 = tmp_109*tmp_8;
      real_t tmp_111 = tmp_109*tmp_26;
      real_t tmp_112 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_113 = tmp_112*tmp_28;
      real_t tmp_114 = tmp_112*tmp_35;
      real_t tmp_115 = tmp_109*tmp_37;
      real_t tmp_116 = tmp_112*tmp_39;
      real_t tmp_117 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_118 = tmp_117*tmp_41;
      real_t tmp_119 = tmp_117*tmp_48;
      real_t tmp_120 = tmp_117*tmp_50;
      real_t tmp_121 = 0.019202922745021479*tmp_55;
      real_t tmp_122 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_123 = tmp_122*tmp_8;
      real_t tmp_124 = tmp_122*tmp_26;
      real_t tmp_125 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_126 = tmp_125*tmp_28;
      real_t tmp_127 = tmp_125*tmp_35;
      real_t tmp_128 = tmp_122*tmp_37;
      real_t tmp_129 = tmp_125*tmp_39;
      real_t tmp_130 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_131 = tmp_130*tmp_41;
      real_t tmp_132 = tmp_130*tmp_48;
      real_t tmp_133 = tmp_130*tmp_50;
      real_t tmp_134 = 0.020848748529055869*tmp_55;
      real_t tmp_135 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_136 = tmp_135*tmp_8;
      real_t tmp_137 = tmp_135*tmp_26;
      real_t tmp_138 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_139 = tmp_138*tmp_28;
      real_t tmp_140 = tmp_138*tmp_35;
      real_t tmp_141 = tmp_135*tmp_37;
      real_t tmp_142 = tmp_138*tmp_39;
      real_t tmp_143 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_144 = tmp_143*tmp_41;
      real_t tmp_145 = tmp_143*tmp_48;
      real_t tmp_146 = tmp_143*tmp_50;
      real_t tmp_147 = 0.019202922745021479*tmp_55;
      real_t tmp_148 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_149 = tmp_148*tmp_8;
      real_t tmp_150 = tmp_148*tmp_26;
      real_t tmp_151 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_152 = tmp_151*tmp_28;
      real_t tmp_153 = tmp_151*tmp_35;
      real_t tmp_154 = tmp_148*tmp_37;
      real_t tmp_155 = tmp_151*tmp_39;
      real_t tmp_156 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_157 = tmp_156*tmp_41;
      real_t tmp_158 = tmp_156*tmp_48;
      real_t tmp_159 = tmp_156*tmp_50;
      real_t tmp_160 = 0.020848748529055869*tmp_55;
      real_t tmp_161 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_162 = tmp_161*tmp_8;
      real_t tmp_163 = tmp_161*tmp_26;
      real_t tmp_164 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_165 = tmp_164*tmp_28;
      real_t tmp_166 = tmp_164*tmp_35;
      real_t tmp_167 = tmp_161*tmp_37;
      real_t tmp_168 = tmp_164*tmp_39;
      real_t tmp_169 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_170 = tmp_169*tmp_41;
      real_t tmp_171 = tmp_169*tmp_48;
      real_t tmp_172 = tmp_169*tmp_50;
      real_t tmp_173 = 0.019202922745021479*tmp_55;
      real_t tmp_174 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_175 = tmp_174*tmp_8;
      real_t tmp_176 = tmp_174*tmp_26;
      real_t tmp_177 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_178 = tmp_177*tmp_28;
      real_t tmp_179 = tmp_177*tmp_35;
      real_t tmp_180 = tmp_174*tmp_37;
      real_t tmp_181 = tmp_177*tmp_39;
      real_t tmp_182 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_183 = tmp_182*tmp_41;
      real_t tmp_184 = tmp_182*tmp_48;
      real_t tmp_185 = tmp_182*tmp_50;
      real_t tmp_186 = 0.020848748529055869*tmp_55;
      real_t tmp_187 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_188 = tmp_187*tmp_8;
      real_t tmp_189 = tmp_187*tmp_26;
      real_t tmp_190 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_191 = tmp_190*tmp_28;
      real_t tmp_192 = tmp_190*tmp_35;
      real_t tmp_193 = tmp_187*tmp_37;
      real_t tmp_194 = tmp_190*tmp_39;
      real_t tmp_195 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_196 = tmp_195*tmp_41;
      real_t tmp_197 = tmp_195*tmp_48;
      real_t tmp_198 = tmp_195*tmp_50;
      real_t tmp_199 = 0.019202922745021479*tmp_55;
      real_t tmp_200 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_201 = tmp_200*tmp_8;
      real_t tmp_202 = tmp_200*tmp_26;
      real_t tmp_203 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_204 = tmp_203*tmp_28;
      real_t tmp_205 = tmp_203*tmp_35;
      real_t tmp_206 = tmp_200*tmp_37;
      real_t tmp_207 = tmp_203*tmp_39;
      real_t tmp_208 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_209 = tmp_208*tmp_41;
      real_t tmp_210 = tmp_208*tmp_48;
      real_t tmp_211 = tmp_208*tmp_50;
      real_t tmp_212 = 0.042507265838595799*tmp_55;
      real_t tmp_213 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_214 = tmp_213*tmp_8;
      real_t tmp_215 = tmp_213*tmp_26;
      real_t tmp_216 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_217 = tmp_216*tmp_28;
      real_t tmp_218 = tmp_216*tmp_35;
      real_t tmp_219 = tmp_213*tmp_37;
      real_t tmp_220 = tmp_216*tmp_39;
      real_t tmp_221 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_222 = tmp_221*tmp_41;
      real_t tmp_223 = tmp_221*tmp_48;
      real_t tmp_224 = tmp_221*tmp_50;
      real_t tmp_225 = 0.020848748529055869*tmp_55;
      real_t tmp_226 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_227 = tmp_226*tmp_8;
      real_t tmp_228 = tmp_226*tmp_26;
      real_t tmp_229 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_230 = tmp_229*tmp_28;
      real_t tmp_231 = tmp_229*tmp_35;
      real_t tmp_232 = tmp_226*tmp_37;
      real_t tmp_233 = tmp_229*tmp_39;
      real_t tmp_234 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_235 = tmp_234*tmp_41;
      real_t tmp_236 = tmp_234*tmp_48;
      real_t tmp_237 = tmp_234*tmp_50;
      real_t tmp_238 = 0.0068572537431980923*tmp_55;
      real_t tmp_239 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_240 = tmp_239*tmp_8;
      real_t tmp_241 = tmp_239*tmp_26;
      real_t tmp_242 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_243 = tmp_242*tmp_28;
      real_t tmp_244 = tmp_242*tmp_35;
      real_t tmp_245 = tmp_239*tmp_37;
      real_t tmp_246 = tmp_242*tmp_39;
      real_t tmp_247 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_248 = tmp_247*tmp_41;
      real_t tmp_249 = tmp_247*tmp_48;
      real_t tmp_250 = tmp_247*tmp_50;
      real_t tmp_251 = 0.037198804536718075*tmp_55;
      real_t tmp_252 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_8;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = 0.042507265838595799*tmp_55;
      real_t tmp_265 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_266 = tmp_265*tmp_8;
      real_t tmp_267 = tmp_26*tmp_265;
      real_t tmp_268 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_269 = tmp_268*tmp_28;
      real_t tmp_270 = tmp_268*tmp_35;
      real_t tmp_271 = tmp_265*tmp_37;
      real_t tmp_272 = tmp_268*tmp_39;
      real_t tmp_273 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_274 = tmp_273*tmp_41;
      real_t tmp_275 = tmp_273*tmp_48;
      real_t tmp_276 = tmp_273*tmp_50;
      real_t tmp_277 = 0.0068572537431980923*tmp_55;
      real_t tmp_278 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_279 = tmp_278*tmp_8;
      real_t tmp_280 = tmp_26*tmp_278;
      real_t tmp_281 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_282 = tmp_28*tmp_281;
      real_t tmp_283 = tmp_281*tmp_35;
      real_t tmp_284 = tmp_278*tmp_37;
      real_t tmp_285 = tmp_281*tmp_39;
      real_t tmp_286 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_287 = tmp_286*tmp_41;
      real_t tmp_288 = tmp_286*tmp_48;
      real_t tmp_289 = tmp_286*tmp_50;
      real_t tmp_290 = 0.037198804536718075*tmp_55;
      real_t tmp_291 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_292 = tmp_291*tmp_8;
      real_t tmp_293 = tmp_26*tmp_291;
      real_t tmp_294 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_295 = tmp_28*tmp_294;
      real_t tmp_296 = tmp_294*tmp_35;
      real_t tmp_297 = tmp_291*tmp_37;
      real_t tmp_298 = tmp_294*tmp_39;
      real_t tmp_299 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_300 = tmp_299*tmp_41;
      real_t tmp_301 = tmp_299*tmp_48;
      real_t tmp_302 = tmp_299*tmp_50;
      real_t tmp_303 = 0.042507265838595799*tmp_55;
      real_t tmp_304 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_305 = tmp_304*tmp_8;
      real_t tmp_306 = tmp_26*tmp_304;
      real_t tmp_307 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_308 = tmp_28*tmp_307;
      real_t tmp_309 = tmp_307*tmp_35;
      real_t tmp_310 = tmp_304*tmp_37;
      real_t tmp_311 = tmp_307*tmp_39;
      real_t tmp_312 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_313 = tmp_312*tmp_41;
      real_t tmp_314 = tmp_312*tmp_48;
      real_t tmp_315 = tmp_312*tmp_50;
      real_t tmp_316 = 0.019202922745021479*tmp_55;
      real_t a_0_0 = tmp_108*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_105 - tmp_106 - tmp_107 - tmp_97 - tmp_98 + 1) + tmp_121*(-tmp_110 - tmp_111 - tmp_113 - tmp_114 - tmp_115 - tmp_116 - tmp_118 - tmp_119 - tmp_120 + 1) + tmp_134*(-tmp_123 - tmp_124 - tmp_126 - tmp_127 - tmp_128 - tmp_129 - tmp_131 - tmp_132 - tmp_133 + 1) + tmp_147*(-tmp_136 - tmp_137 - tmp_139 - tmp_140 - tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 + 1) + tmp_160*(-tmp_149 - tmp_150 - tmp_152 - tmp_153 - tmp_154 - tmp_155 - tmp_157 - tmp_158 - tmp_159 + 1) + tmp_173*(-tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 - tmp_168 - tmp_170 - tmp_171 - tmp_172 + 1) + tmp_186*(-tmp_175 - tmp_176 - tmp_178 - tmp_179 - tmp_180 - tmp_181 - tmp_183 - tmp_184 - tmp_185 + 1) + tmp_199*(-tmp_188 - tmp_189 - tmp_191 - tmp_192 - tmp_193 - tmp_194 - tmp_196 - tmp_197 - tmp_198 + 1) + tmp_212*(-tmp_201 - tmp_202 - tmp_204 - tmp_205 - tmp_206 - tmp_207 - tmp_209 - tmp_210 - tmp_211 + 1) + tmp_225*(-tmp_214 - tmp_215 - tmp_217 - tmp_218 - tmp_219 - tmp_220 - tmp_222 - tmp_223 - tmp_224 + 1) + tmp_238*(-tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 - tmp_233 - tmp_235 - tmp_236 - tmp_237 + 1) + tmp_251*(-tmp_240 - tmp_241 - tmp_243 - tmp_244 - tmp_245 - tmp_246 - tmp_248 - tmp_249 - tmp_250 + 1) + tmp_264*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_277*(-tmp_266 - tmp_267 - tmp_269 - tmp_270 - tmp_271 - tmp_272 - tmp_274 - tmp_275 - tmp_276 + 1) + tmp_290*(-tmp_279 - tmp_280 - tmp_282 - tmp_283 - tmp_284 - tmp_285 - tmp_287 - tmp_288 - tmp_289 + 1) + tmp_303*(-tmp_292 - tmp_293 - tmp_295 - tmp_296 - tmp_297 - tmp_298 - tmp_300 - tmp_301 - tmp_302 + 1) + tmp_316*(-tmp_305 - tmp_306 - tmp_308 - tmp_309 - tmp_310 - tmp_311 - tmp_313 - tmp_314 - tmp_315 + 1) + tmp_56*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_69*(-tmp_58 - tmp_59 - tmp_61 - tmp_62 - tmp_63 - tmp_64 - tmp_66 - tmp_67 - tmp_68 + 1) + tmp_82*(-tmp_71 - tmp_72 - tmp_74 - tmp_75 - tmp_76 - tmp_77 - tmp_79 - tmp_80 - tmp_81 + 1) + tmp_95*(-tmp_84 - tmp_85 - tmp_87 - tmp_88 - tmp_89 - tmp_90 - tmp_92 - tmp_93 - tmp_94 + 1);
      real_t a_1_0 = tmp_108*(tmp_102 + tmp_103 + tmp_107) + tmp_121*(tmp_115 + tmp_116 + tmp_120) + tmp_134*(tmp_128 + tmp_129 + tmp_133) + tmp_147*(tmp_141 + tmp_142 + tmp_146) + tmp_160*(tmp_154 + tmp_155 + tmp_159) + tmp_173*(tmp_167 + tmp_168 + tmp_172) + tmp_186*(tmp_180 + tmp_181 + tmp_185) + tmp_199*(tmp_193 + tmp_194 + tmp_198) + tmp_212*(tmp_206 + tmp_207 + tmp_211) + tmp_225*(tmp_219 + tmp_220 + tmp_224) + tmp_238*(tmp_232 + tmp_233 + tmp_237) + tmp_251*(tmp_245 + tmp_246 + tmp_250) + tmp_264*(tmp_258 + tmp_259 + tmp_263) + tmp_277*(tmp_271 + tmp_272 + tmp_276) + tmp_290*(tmp_284 + tmp_285 + tmp_289) + tmp_303*(tmp_297 + tmp_298 + tmp_302) + tmp_316*(tmp_310 + tmp_311 + tmp_315) + tmp_56*(tmp_38 + tmp_40 + tmp_51) + tmp_69*(tmp_63 + tmp_64 + tmp_68) + tmp_82*(tmp_76 + tmp_77 + tmp_81) + tmp_95*(tmp_89 + tmp_90 + tmp_94);
      real_t a_2_0 = tmp_108*(tmp_101 + tmp_106 + tmp_98) + tmp_121*(tmp_111 + tmp_114 + tmp_119) + tmp_134*(tmp_124 + tmp_127 + tmp_132) + tmp_147*(tmp_137 + tmp_140 + tmp_145) + tmp_160*(tmp_150 + tmp_153 + tmp_158) + tmp_173*(tmp_163 + tmp_166 + tmp_171) + tmp_186*(tmp_176 + tmp_179 + tmp_184) + tmp_199*(tmp_189 + tmp_192 + tmp_197) + tmp_212*(tmp_202 + tmp_205 + tmp_210) + tmp_225*(tmp_215 + tmp_218 + tmp_223) + tmp_238*(tmp_228 + tmp_231 + tmp_236) + tmp_251*(tmp_241 + tmp_244 + tmp_249) + tmp_264*(tmp_254 + tmp_257 + tmp_262) + tmp_277*(tmp_267 + tmp_270 + tmp_275) + tmp_290*(tmp_280 + tmp_283 + tmp_288) + tmp_303*(tmp_293 + tmp_296 + tmp_301) + tmp_316*(tmp_306 + tmp_309 + tmp_314) + tmp_56*(tmp_27 + tmp_36 + tmp_49) + tmp_69*(tmp_59 + tmp_62 + tmp_67) + tmp_82*(tmp_72 + tmp_75 + tmp_80) + tmp_95*(tmp_85 + tmp_88 + tmp_93);
      real_t a_3_0 = tmp_108*(tmp_100 + tmp_105 + tmp_97) + tmp_121*(tmp_110 + tmp_113 + tmp_118) + tmp_134*(tmp_123 + tmp_126 + tmp_131) + tmp_147*(tmp_136 + tmp_139 + tmp_144) + tmp_160*(tmp_149 + tmp_152 + tmp_157) + tmp_173*(tmp_162 + tmp_165 + tmp_170) + tmp_186*(tmp_175 + tmp_178 + tmp_183) + tmp_199*(tmp_188 + tmp_191 + tmp_196) + tmp_212*(tmp_201 + tmp_204 + tmp_209) + tmp_225*(tmp_214 + tmp_217 + tmp_222) + tmp_238*(tmp_227 + tmp_230 + tmp_235) + tmp_251*(tmp_240 + tmp_243 + tmp_248) + tmp_264*(tmp_253 + tmp_256 + tmp_261) + tmp_277*(tmp_266 + tmp_269 + tmp_274) + tmp_290*(tmp_279 + tmp_282 + tmp_287) + tmp_303*(tmp_292 + tmp_295 + tmp_300) + tmp_316*(tmp_305 + tmp_308 + tmp_313) + tmp_56*(tmp_25 + tmp_34 + tmp_47) + tmp_69*(tmp_58 + tmp_61 + tmp_66) + tmp_82*(tmp_71 + tmp_74 + tmp_79) + tmp_95*(tmp_84 + tmp_87 + tmp_92);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
    const Eigen::Matrix< real_t, 3, 1 >&,
    const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
    Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = p_affine_1_2 + tmp_9;
      real_t tmp_13 = tmp_12*tmp_5;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = p_affine_2_2 + tmp_9;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_1*tmp_11;
      real_t tmp_18 = tmp_12*tmp_14;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_4 - tmp_10*tmp_7 + tmp_11*tmp_13 + tmp_14*tmp_16 - tmp_15*tmp_17 - tmp_18*tmp_3);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_9;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_8;
      real_t tmp_26 = tmp_14*tmp_6 - tmp_17;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_15 + tmp_13;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_10 - tmp_18;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_11*tmp_5 - tmp_14*tmp_3;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_12*tmp_3 + tmp_16;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = tmp_10*tmp_3 - tmp_11*tmp_15;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_53 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_54 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_55 = 1.0*p_affine_13_1*std::pow((std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)*std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)) + (std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)*std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)) + (std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)*std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)), 1.0/2.0);
      real_t tmp_56 = 0.0068572537431980923*tmp_55;
      real_t tmp_57 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_58 = tmp_57*tmp_8;
      real_t tmp_59 = tmp_26*tmp_57;
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_61 = tmp_28*tmp_60;
      real_t tmp_62 = tmp_35*tmp_60;
      real_t tmp_63 = tmp_37*tmp_57;
      real_t tmp_64 = tmp_39*tmp_60;
      real_t tmp_65 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_66 = tmp_41*tmp_65;
      real_t tmp_67 = tmp_48*tmp_65;
      real_t tmp_68 = tmp_50*tmp_65;
      real_t tmp_69 = 0.037198804536718075*tmp_55;
      real_t tmp_70 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_71 = tmp_70*tmp_8;
      real_t tmp_72 = tmp_26*tmp_70;
      real_t tmp_73 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_74 = tmp_28*tmp_73;
      real_t tmp_75 = tmp_35*tmp_73;
      real_t tmp_76 = tmp_37*tmp_70;
      real_t tmp_77 = tmp_39*tmp_73;
      real_t tmp_78 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_79 = tmp_41*tmp_78;
      real_t tmp_80 = tmp_48*tmp_78;
      real_t tmp_81 = tmp_50*tmp_78;
      real_t tmp_82 = 0.020848748529055869*tmp_55;
      real_t tmp_83 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_84 = tmp_8*tmp_83;
      real_t tmp_85 = tmp_26*tmp_83;
      real_t tmp_86 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_87 = tmp_28*tmp_86;
      real_t tmp_88 = tmp_35*tmp_86;
      real_t tmp_89 = tmp_37*tmp_83;
      real_t tmp_90 = tmp_39*tmp_86;
      real_t tmp_91 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_92 = tmp_41*tmp_91;
      real_t tmp_93 = tmp_48*tmp_91;
      real_t tmp_94 = tmp_50*tmp_91;
      real_t tmp_95 = 0.019202922745021479*tmp_55;
      real_t tmp_96 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_97 = tmp_8*tmp_96;
      real_t tmp_98 = tmp_26*tmp_96;
      real_t tmp_99 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_100 = tmp_28*tmp_99;
      real_t tmp_101 = tmp_35*tmp_99;
      real_t tmp_102 = tmp_37*tmp_96;
      real_t tmp_103 = tmp_39*tmp_99;
      real_t tmp_104 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_105 = tmp_104*tmp_41;
      real_t tmp_106 = tmp_104*tmp_48;
      real_t tmp_107 = tmp_104*tmp_50;
      real_t tmp_108 = 0.020848748529055869*tmp_55;
      real_t tmp_109 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_110 = tmp_109*tmp_8;
      real_t tmp_111 = tmp_109*tmp_26;
      real_t tmp_112 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_113 = tmp_112*tmp_28;
      real_t tmp_114 = tmp_112*tmp_35;
      real_t tmp_115 = tmp_109*tmp_37;
      real_t tmp_116 = tmp_112*tmp_39;
      real_t tmp_117 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_118 = tmp_117*tmp_41;
      real_t tmp_119 = tmp_117*tmp_48;
      real_t tmp_120 = tmp_117*tmp_50;
      real_t tmp_121 = 0.019202922745021479*tmp_55;
      real_t tmp_122 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_123 = tmp_122*tmp_8;
      real_t tmp_124 = tmp_122*tmp_26;
      real_t tmp_125 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_126 = tmp_125*tmp_28;
      real_t tmp_127 = tmp_125*tmp_35;
      real_t tmp_128 = tmp_122*tmp_37;
      real_t tmp_129 = tmp_125*tmp_39;
      real_t tmp_130 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_131 = tmp_130*tmp_41;
      real_t tmp_132 = tmp_130*tmp_48;
      real_t tmp_133 = tmp_130*tmp_50;
      real_t tmp_134 = 0.020848748529055869*tmp_55;
      real_t tmp_135 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_136 = tmp_135*tmp_8;
      real_t tmp_137 = tmp_135*tmp_26;
      real_t tmp_138 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_139 = tmp_138*tmp_28;
      real_t tmp_140 = tmp_138*tmp_35;
      real_t tmp_141 = tmp_135*tmp_37;
      real_t tmp_142 = tmp_138*tmp_39;
      real_t tmp_143 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_144 = tmp_143*tmp_41;
      real_t tmp_145 = tmp_143*tmp_48;
      real_t tmp_146 = tmp_143*tmp_50;
      real_t tmp_147 = 0.019202922745021479*tmp_55;
      real_t tmp_148 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_149 = tmp_148*tmp_8;
      real_t tmp_150 = tmp_148*tmp_26;
      real_t tmp_151 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_152 = tmp_151*tmp_28;
      real_t tmp_153 = tmp_151*tmp_35;
      real_t tmp_154 = tmp_148*tmp_37;
      real_t tmp_155 = tmp_151*tmp_39;
      real_t tmp_156 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_157 = tmp_156*tmp_41;
      real_t tmp_158 = tmp_156*tmp_48;
      real_t tmp_159 = tmp_156*tmp_50;
      real_t tmp_160 = 0.020848748529055869*tmp_55;
      real_t tmp_161 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_162 = tmp_161*tmp_8;
      real_t tmp_163 = tmp_161*tmp_26;
      real_t tmp_164 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_165 = tmp_164*tmp_28;
      real_t tmp_166 = tmp_164*tmp_35;
      real_t tmp_167 = tmp_161*tmp_37;
      real_t tmp_168 = tmp_164*tmp_39;
      real_t tmp_169 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_170 = tmp_169*tmp_41;
      real_t tmp_171 = tmp_169*tmp_48;
      real_t tmp_172 = tmp_169*tmp_50;
      real_t tmp_173 = 0.019202922745021479*tmp_55;
      real_t tmp_174 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_175 = tmp_174*tmp_8;
      real_t tmp_176 = tmp_174*tmp_26;
      real_t tmp_177 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_178 = tmp_177*tmp_28;
      real_t tmp_179 = tmp_177*tmp_35;
      real_t tmp_180 = tmp_174*tmp_37;
      real_t tmp_181 = tmp_177*tmp_39;
      real_t tmp_182 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_183 = tmp_182*tmp_41;
      real_t tmp_184 = tmp_182*tmp_48;
      real_t tmp_185 = tmp_182*tmp_50;
      real_t tmp_186 = 0.020848748529055869*tmp_55;
      real_t tmp_187 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_188 = tmp_187*tmp_8;
      real_t tmp_189 = tmp_187*tmp_26;
      real_t tmp_190 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_191 = tmp_190*tmp_28;
      real_t tmp_192 = tmp_190*tmp_35;
      real_t tmp_193 = tmp_187*tmp_37;
      real_t tmp_194 = tmp_190*tmp_39;
      real_t tmp_195 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_196 = tmp_195*tmp_41;
      real_t tmp_197 = tmp_195*tmp_48;
      real_t tmp_198 = tmp_195*tmp_50;
      real_t tmp_199 = 0.019202922745021479*tmp_55;
      real_t tmp_200 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_201 = tmp_200*tmp_8;
      real_t tmp_202 = tmp_200*tmp_26;
      real_t tmp_203 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_204 = tmp_203*tmp_28;
      real_t tmp_205 = tmp_203*tmp_35;
      real_t tmp_206 = tmp_200*tmp_37;
      real_t tmp_207 = tmp_203*tmp_39;
      real_t tmp_208 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_209 = tmp_208*tmp_41;
      real_t tmp_210 = tmp_208*tmp_48;
      real_t tmp_211 = tmp_208*tmp_50;
      real_t tmp_212 = 0.042507265838595799*tmp_55;
      real_t tmp_213 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_214 = tmp_213*tmp_8;
      real_t tmp_215 = tmp_213*tmp_26;
      real_t tmp_216 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_217 = tmp_216*tmp_28;
      real_t tmp_218 = tmp_216*tmp_35;
      real_t tmp_219 = tmp_213*tmp_37;
      real_t tmp_220 = tmp_216*tmp_39;
      real_t tmp_221 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_222 = tmp_221*tmp_41;
      real_t tmp_223 = tmp_221*tmp_48;
      real_t tmp_224 = tmp_221*tmp_50;
      real_t tmp_225 = 0.020848748529055869*tmp_55;
      real_t tmp_226 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_227 = tmp_226*tmp_8;
      real_t tmp_228 = tmp_226*tmp_26;
      real_t tmp_229 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_230 = tmp_229*tmp_28;
      real_t tmp_231 = tmp_229*tmp_35;
      real_t tmp_232 = tmp_226*tmp_37;
      real_t tmp_233 = tmp_229*tmp_39;
      real_t tmp_234 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_235 = tmp_234*tmp_41;
      real_t tmp_236 = tmp_234*tmp_48;
      real_t tmp_237 = tmp_234*tmp_50;
      real_t tmp_238 = 0.0068572537431980923*tmp_55;
      real_t tmp_239 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_240 = tmp_239*tmp_8;
      real_t tmp_241 = tmp_239*tmp_26;
      real_t tmp_242 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_243 = tmp_242*tmp_28;
      real_t tmp_244 = tmp_242*tmp_35;
      real_t tmp_245 = tmp_239*tmp_37;
      real_t tmp_246 = tmp_242*tmp_39;
      real_t tmp_247 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_248 = tmp_247*tmp_41;
      real_t tmp_249 = tmp_247*tmp_48;
      real_t tmp_250 = tmp_247*tmp_50;
      real_t tmp_251 = 0.037198804536718075*tmp_55;
      real_t tmp_252 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_8;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = 0.042507265838595799*tmp_55;
      real_t tmp_265 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_266 = tmp_265*tmp_8;
      real_t tmp_267 = tmp_26*tmp_265;
      real_t tmp_268 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_269 = tmp_268*tmp_28;
      real_t tmp_270 = tmp_268*tmp_35;
      real_t tmp_271 = tmp_265*tmp_37;
      real_t tmp_272 = tmp_268*tmp_39;
      real_t tmp_273 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_274 = tmp_273*tmp_41;
      real_t tmp_275 = tmp_273*tmp_48;
      real_t tmp_276 = tmp_273*tmp_50;
      real_t tmp_277 = 0.0068572537431980923*tmp_55;
      real_t tmp_278 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_279 = tmp_278*tmp_8;
      real_t tmp_280 = tmp_26*tmp_278;
      real_t tmp_281 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_282 = tmp_28*tmp_281;
      real_t tmp_283 = tmp_281*tmp_35;
      real_t tmp_284 = tmp_278*tmp_37;
      real_t tmp_285 = tmp_281*tmp_39;
      real_t tmp_286 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_287 = tmp_286*tmp_41;
      real_t tmp_288 = tmp_286*tmp_48;
      real_t tmp_289 = tmp_286*tmp_50;
      real_t tmp_290 = 0.037198804536718075*tmp_55;
      real_t tmp_291 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_292 = tmp_291*tmp_8;
      real_t tmp_293 = tmp_26*tmp_291;
      real_t tmp_294 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_295 = tmp_28*tmp_294;
      real_t tmp_296 = tmp_294*tmp_35;
      real_t tmp_297 = tmp_291*tmp_37;
      real_t tmp_298 = tmp_294*tmp_39;
      real_t tmp_299 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_300 = tmp_299*tmp_41;
      real_t tmp_301 = tmp_299*tmp_48;
      real_t tmp_302 = tmp_299*tmp_50;
      real_t tmp_303 = 0.042507265838595799*tmp_55;
      real_t tmp_304 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_305 = tmp_304*tmp_8;
      real_t tmp_306 = tmp_26*tmp_304;
      real_t tmp_307 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_308 = tmp_28*tmp_307;
      real_t tmp_309 = tmp_307*tmp_35;
      real_t tmp_310 = tmp_304*tmp_37;
      real_t tmp_311 = tmp_307*tmp_39;
      real_t tmp_312 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_313 = tmp_312*tmp_41;
      real_t tmp_314 = tmp_312*tmp_48;
      real_t tmp_315 = tmp_312*tmp_50;
      real_t tmp_316 = 0.019202922745021479*tmp_55;
      real_t a_0_0 = tmp_108*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_105 - tmp_106 - tmp_107 - tmp_97 - tmp_98 + 1) + tmp_121*(-tmp_110 - tmp_111 - tmp_113 - tmp_114 - tmp_115 - tmp_116 - tmp_118 - tmp_119 - tmp_120 + 1) + tmp_134*(-tmp_123 - tmp_124 - tmp_126 - tmp_127 - tmp_128 - tmp_129 - tmp_131 - tmp_132 - tmp_133 + 1) + tmp_147*(-tmp_136 - tmp_137 - tmp_139 - tmp_140 - tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 + 1) + tmp_160*(-tmp_149 - tmp_150 - tmp_152 - tmp_153 - tmp_154 - tmp_155 - tmp_157 - tmp_158 - tmp_159 + 1) + tmp_173*(-tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 - tmp_168 - tmp_170 - tmp_171 - tmp_172 + 1) + tmp_186*(-tmp_175 - tmp_176 - tmp_178 - tmp_179 - tmp_180 - tmp_181 - tmp_183 - tmp_184 - tmp_185 + 1) + tmp_199*(-tmp_188 - tmp_189 - tmp_191 - tmp_192 - tmp_193 - tmp_194 - tmp_196 - tmp_197 - tmp_198 + 1) + tmp_212*(-tmp_201 - tmp_202 - tmp_204 - tmp_205 - tmp_206 - tmp_207 - tmp_209 - tmp_210 - tmp_211 + 1) + tmp_225*(-tmp_214 - tmp_215 - tmp_217 - tmp_218 - tmp_219 - tmp_220 - tmp_222 - tmp_223 - tmp_224 + 1) + tmp_238*(-tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 - tmp_233 - tmp_235 - tmp_236 - tmp_237 + 1) + tmp_251*(-tmp_240 - tmp_241 - tmp_243 - tmp_244 - tmp_245 - tmp_246 - tmp_248 - tmp_249 - tmp_250 + 1) + tmp_264*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_277*(-tmp_266 - tmp_267 - tmp_269 - tmp_270 - tmp_271 - tmp_272 - tmp_274 - tmp_275 - tmp_276 + 1) + tmp_290*(-tmp_279 - tmp_280 - tmp_282 - tmp_283 - tmp_284 - tmp_285 - tmp_287 - tmp_288 - tmp_289 + 1) + tmp_303*(-tmp_292 - tmp_293 - tmp_295 - tmp_296 - tmp_297 - tmp_298 - tmp_300 - tmp_301 - tmp_302 + 1) + tmp_316*(-tmp_305 - tmp_306 - tmp_308 - tmp_309 - tmp_310 - tmp_311 - tmp_313 - tmp_314 - tmp_315 + 1) + tmp_56*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_69*(-tmp_58 - tmp_59 - tmp_61 - tmp_62 - tmp_63 - tmp_64 - tmp_66 - tmp_67 - tmp_68 + 1) + tmp_82*(-tmp_71 - tmp_72 - tmp_74 - tmp_75 - tmp_76 - tmp_77 - tmp_79 - tmp_80 - tmp_81 + 1) + tmp_95*(-tmp_84 - tmp_85 - tmp_87 - tmp_88 - tmp_89 - tmp_90 - tmp_92 - tmp_93 - tmp_94 + 1);
      real_t a_1_0 = tmp_108*(tmp_102 + tmp_103 + tmp_107) + tmp_121*(tmp_115 + tmp_116 + tmp_120) + tmp_134*(tmp_128 + tmp_129 + tmp_133) + tmp_147*(tmp_141 + tmp_142 + tmp_146) + tmp_160*(tmp_154 + tmp_155 + tmp_159) + tmp_173*(tmp_167 + tmp_168 + tmp_172) + tmp_186*(tmp_180 + tmp_181 + tmp_185) + tmp_199*(tmp_193 + tmp_194 + tmp_198) + tmp_212*(tmp_206 + tmp_207 + tmp_211) + tmp_225*(tmp_219 + tmp_220 + tmp_224) + tmp_238*(tmp_232 + tmp_233 + tmp_237) + tmp_251*(tmp_245 + tmp_246 + tmp_250) + tmp_264*(tmp_258 + tmp_259 + tmp_263) + tmp_277*(tmp_271 + tmp_272 + tmp_276) + tmp_290*(tmp_284 + tmp_285 + tmp_289) + tmp_303*(tmp_297 + tmp_298 + tmp_302) + tmp_316*(tmp_310 + tmp_311 + tmp_315) + tmp_56*(tmp_38 + tmp_40 + tmp_51) + tmp_69*(tmp_63 + tmp_64 + tmp_68) + tmp_82*(tmp_76 + tmp_77 + tmp_81) + tmp_95*(tmp_89 + tmp_90 + tmp_94);
      real_t a_2_0 = tmp_108*(tmp_101 + tmp_106 + tmp_98) + tmp_121*(tmp_111 + tmp_114 + tmp_119) + tmp_134*(tmp_124 + tmp_127 + tmp_132) + tmp_147*(tmp_137 + tmp_140 + tmp_145) + tmp_160*(tmp_150 + tmp_153 + tmp_158) + tmp_173*(tmp_163 + tmp_166 + tmp_171) + tmp_186*(tmp_176 + tmp_179 + tmp_184) + tmp_199*(tmp_189 + tmp_192 + tmp_197) + tmp_212*(tmp_202 + tmp_205 + tmp_210) + tmp_225*(tmp_215 + tmp_218 + tmp_223) + tmp_238*(tmp_228 + tmp_231 + tmp_236) + tmp_251*(tmp_241 + tmp_244 + tmp_249) + tmp_264*(tmp_254 + tmp_257 + tmp_262) + tmp_277*(tmp_267 + tmp_270 + tmp_275) + tmp_290*(tmp_280 + tmp_283 + tmp_288) + tmp_303*(tmp_293 + tmp_296 + tmp_301) + tmp_316*(tmp_306 + tmp_309 + tmp_314) + tmp_56*(tmp_27 + tmp_36 + tmp_49) + tmp_69*(tmp_59 + tmp_62 + tmp_67) + tmp_82*(tmp_72 + tmp_75 + tmp_80) + tmp_95*(tmp_85 + tmp_88 + tmp_93);
      real_t a_3_0 = tmp_108*(tmp_100 + tmp_105 + tmp_97) + tmp_121*(tmp_110 + tmp_113 + tmp_118) + tmp_134*(tmp_123 + tmp_126 + tmp_131) + tmp_147*(tmp_136 + tmp_139 + tmp_144) + tmp_160*(tmp_149 + tmp_152 + tmp_157) + tmp_173*(tmp_162 + tmp_165 + tmp_170) + tmp_186*(tmp_175 + tmp_178 + tmp_183) + tmp_199*(tmp_188 + tmp_191 + tmp_196) + tmp_212*(tmp_201 + tmp_204 + tmp_209) + tmp_225*(tmp_214 + tmp_217 + tmp_222) + tmp_238*(tmp_227 + tmp_230 + tmp_235) + tmp_251*(tmp_240 + tmp_243 + tmp_248) + tmp_264*(tmp_253 + tmp_256 + tmp_261) + tmp_277*(tmp_266 + tmp_269 + tmp_274) + tmp_290*(tmp_279 + tmp_282 + tmp_287) + tmp_303*(tmp_292 + tmp_295 + tmp_300) + tmp_316*(tmp_305 + tmp_308 + tmp_313) + tmp_56*(tmp_25 + tmp_34 + tmp_47) + tmp_69*(tmp_58 + tmp_61 + tmp_66) + tmp_82*(tmp_71 + tmp_74 + tmp_79) + tmp_95*(tmp_84 + tmp_87 + tmp_92);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }

public:



};




class EGDivtForm_P1P0_2 : public hyteg::dg::DGForm
{

 public:
    EGDivtForm_P1P0_2()

    {}





 protected:
  void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
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
      real_t a_1_0 = 0;
      real_t a_2_0 = 0;
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
      real_t a_1_0 = 0;
      real_t a_2_0 = 0;
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
      real_t a_1_0 = 0;
      real_t a_2_0 = 0;
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
      real_t a_1_0 = 0;
      real_t a_2_0 = 0;
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

      elMat( 0, 0) = 0;
      elMat( 1, 0) = 0;
      elMat( 2, 0) = 0;
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
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

      elMat( 0, 0) = 0;
      elMat( 1, 0) = 0;
      elMat( 2, 0) = 0;
      elMat( 3, 0) = 0;
   }
   void integrateVolume3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = -p_affine_0_2;
      real_t tmp_9 = p_affine_3_2 + tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_8;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = p_affine_2_2 + tmp_8;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = tmp_14*tmp_6;
      real_t tmp_16 = tmp_1*tmp_11;
      real_t tmp_17 = tmp_14*tmp_3;
      real_t tmp_18 = 1.0 / (tmp_10*tmp_12 - tmp_10*tmp_17 + tmp_13*tmp_15 - tmp_13*tmp_16 + tmp_4*tmp_9 - tmp_7*tmp_9);
      real_t tmp_19 = tmp_18*(tmp_4 - tmp_7);
      real_t tmp_20 = tmp_18*(tmp_15 - tmp_16);
      real_t tmp_21 = tmp_18*(tmp_12 - tmp_17);
      real_t tmp_22 = p_affine_0_0*p_affine_1_1;
      real_t tmp_23 = p_affine_0_0*p_affine_1_2;
      real_t tmp_24 = p_affine_2_1*p_affine_3_2;
      real_t tmp_25 = p_affine_0_1*p_affine_1_0;
      real_t tmp_26 = p_affine_0_1*p_affine_1_2;
      real_t tmp_27 = p_affine_2_2*p_affine_3_0;
      real_t tmp_28 = p_affine_0_2*p_affine_1_0;
      real_t tmp_29 = p_affine_0_2*p_affine_1_1;
      real_t tmp_30 = p_affine_2_0*p_affine_3_1;
      real_t tmp_31 = p_affine_2_2*p_affine_3_1;
      real_t tmp_32 = p_affine_2_0*p_affine_3_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_0;
      real_t tmp_34 = std::abs(p_affine_0_0*tmp_24 - p_affine_0_0*tmp_31 + p_affine_0_1*tmp_27 - p_affine_0_1*tmp_32 + p_affine_0_2*tmp_30 - p_affine_0_2*tmp_33 - p_affine_1_0*tmp_24 + p_affine_1_0*tmp_31 - p_affine_1_1*tmp_27 + p_affine_1_1*tmp_32 - p_affine_1_2*tmp_30 + p_affine_1_2*tmp_33 + p_affine_2_0*tmp_26 - p_affine_2_0*tmp_29 - p_affine_2_1*tmp_23 + p_affine_2_1*tmp_28 + p_affine_2_2*tmp_22 - p_affine_2_2*tmp_25 - p_affine_3_0*tmp_26 + p_affine_3_0*tmp_29 + p_affine_3_1*tmp_23 - p_affine_3_1*tmp_28 - p_affine_3_2*tmp_22 + p_affine_3_2*tmp_25);
      real_t tmp_35 = tmp_34*(tmp_19 + tmp_20 + tmp_21);
      real_t tmp_36 = tmp_21*tmp_34;
      real_t tmp_37 = tmp_20*tmp_34;
      real_t tmp_38 = tmp_19*tmp_34;
      real_t a_0_0 = 0.1666666666666668*tmp_35;
      real_t a_1_0 = -0.1666666666666668*tmp_36;
      real_t a_2_0 = -0.1666666666666668*tmp_37;
      real_t a_3_0 = -0.1666666666666668*tmp_38;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }



   void integrateFacetInner3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
                                                     const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                     const Eigen::Matrix< real_t, 3, 1 >&,
                                                     const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                                                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

         real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = p_affine_1_2 + tmp_9;
      real_t tmp_13 = tmp_12*tmp_5;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = p_affine_2_2 + tmp_9;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_1*tmp_11;
      real_t tmp_18 = tmp_12*tmp_14;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_4 - tmp_10*tmp_7 + tmp_11*tmp_13 + tmp_14*tmp_16 - tmp_15*tmp_17 - tmp_18*tmp_3);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_9;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_8;
      real_t tmp_26 = tmp_14*tmp_6 - tmp_17;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_15 + tmp_13;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_10 - tmp_18;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_11*tmp_5 - tmp_14*tmp_3;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_12*tmp_3 + tmp_16;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = tmp_10*tmp_3 - tmp_11*tmp_15;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_53 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_54 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_55 = 0.5*p_affine_13_2*std::pow((std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)*std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)) + (std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)*std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)) + (std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)*std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)), 1.0/2.0);
      real_t tmp_56 = 0.0068572537431980923*tmp_55;
      real_t tmp_57 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_58 = tmp_57*tmp_8;
      real_t tmp_59 = tmp_26*tmp_57;
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_61 = tmp_28*tmp_60;
      real_t tmp_62 = tmp_35*tmp_60;
      real_t tmp_63 = tmp_37*tmp_57;
      real_t tmp_64 = tmp_39*tmp_60;
      real_t tmp_65 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_66 = tmp_41*tmp_65;
      real_t tmp_67 = tmp_48*tmp_65;
      real_t tmp_68 = tmp_50*tmp_65;
      real_t tmp_69 = 0.037198804536718075*tmp_55;
      real_t tmp_70 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_71 = tmp_70*tmp_8;
      real_t tmp_72 = tmp_26*tmp_70;
      real_t tmp_73 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_74 = tmp_28*tmp_73;
      real_t tmp_75 = tmp_35*tmp_73;
      real_t tmp_76 = tmp_37*tmp_70;
      real_t tmp_77 = tmp_39*tmp_73;
      real_t tmp_78 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_79 = tmp_41*tmp_78;
      real_t tmp_80 = tmp_48*tmp_78;
      real_t tmp_81 = tmp_50*tmp_78;
      real_t tmp_82 = 0.020848748529055869*tmp_55;
      real_t tmp_83 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_84 = tmp_8*tmp_83;
      real_t tmp_85 = tmp_26*tmp_83;
      real_t tmp_86 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_87 = tmp_28*tmp_86;
      real_t tmp_88 = tmp_35*tmp_86;
      real_t tmp_89 = tmp_37*tmp_83;
      real_t tmp_90 = tmp_39*tmp_86;
      real_t tmp_91 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_92 = tmp_41*tmp_91;
      real_t tmp_93 = tmp_48*tmp_91;
      real_t tmp_94 = tmp_50*tmp_91;
      real_t tmp_95 = 0.019202922745021479*tmp_55;
      real_t tmp_96 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_97 = tmp_8*tmp_96;
      real_t tmp_98 = tmp_26*tmp_96;
      real_t tmp_99 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_100 = tmp_28*tmp_99;
      real_t tmp_101 = tmp_35*tmp_99;
      real_t tmp_102 = tmp_37*tmp_96;
      real_t tmp_103 = tmp_39*tmp_99;
      real_t tmp_104 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_105 = tmp_104*tmp_41;
      real_t tmp_106 = tmp_104*tmp_48;
      real_t tmp_107 = tmp_104*tmp_50;
      real_t tmp_108 = 0.020848748529055869*tmp_55;
      real_t tmp_109 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_110 = tmp_109*tmp_8;
      real_t tmp_111 = tmp_109*tmp_26;
      real_t tmp_112 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_113 = tmp_112*tmp_28;
      real_t tmp_114 = tmp_112*tmp_35;
      real_t tmp_115 = tmp_109*tmp_37;
      real_t tmp_116 = tmp_112*tmp_39;
      real_t tmp_117 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_118 = tmp_117*tmp_41;
      real_t tmp_119 = tmp_117*tmp_48;
      real_t tmp_120 = tmp_117*tmp_50;
      real_t tmp_121 = 0.019202922745021479*tmp_55;
      real_t tmp_122 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_123 = tmp_122*tmp_8;
      real_t tmp_124 = tmp_122*tmp_26;
      real_t tmp_125 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_126 = tmp_125*tmp_28;
      real_t tmp_127 = tmp_125*tmp_35;
      real_t tmp_128 = tmp_122*tmp_37;
      real_t tmp_129 = tmp_125*tmp_39;
      real_t tmp_130 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_131 = tmp_130*tmp_41;
      real_t tmp_132 = tmp_130*tmp_48;
      real_t tmp_133 = tmp_130*tmp_50;
      real_t tmp_134 = 0.020848748529055869*tmp_55;
      real_t tmp_135 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_136 = tmp_135*tmp_8;
      real_t tmp_137 = tmp_135*tmp_26;
      real_t tmp_138 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_139 = tmp_138*tmp_28;
      real_t tmp_140 = tmp_138*tmp_35;
      real_t tmp_141 = tmp_135*tmp_37;
      real_t tmp_142 = tmp_138*tmp_39;
      real_t tmp_143 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_144 = tmp_143*tmp_41;
      real_t tmp_145 = tmp_143*tmp_48;
      real_t tmp_146 = tmp_143*tmp_50;
      real_t tmp_147 = 0.019202922745021479*tmp_55;
      real_t tmp_148 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_149 = tmp_148*tmp_8;
      real_t tmp_150 = tmp_148*tmp_26;
      real_t tmp_151 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_152 = tmp_151*tmp_28;
      real_t tmp_153 = tmp_151*tmp_35;
      real_t tmp_154 = tmp_148*tmp_37;
      real_t tmp_155 = tmp_151*tmp_39;
      real_t tmp_156 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_157 = tmp_156*tmp_41;
      real_t tmp_158 = tmp_156*tmp_48;
      real_t tmp_159 = tmp_156*tmp_50;
      real_t tmp_160 = 0.020848748529055869*tmp_55;
      real_t tmp_161 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_162 = tmp_161*tmp_8;
      real_t tmp_163 = tmp_161*tmp_26;
      real_t tmp_164 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_165 = tmp_164*tmp_28;
      real_t tmp_166 = tmp_164*tmp_35;
      real_t tmp_167 = tmp_161*tmp_37;
      real_t tmp_168 = tmp_164*tmp_39;
      real_t tmp_169 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_170 = tmp_169*tmp_41;
      real_t tmp_171 = tmp_169*tmp_48;
      real_t tmp_172 = tmp_169*tmp_50;
      real_t tmp_173 = 0.019202922745021479*tmp_55;
      real_t tmp_174 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_175 = tmp_174*tmp_8;
      real_t tmp_176 = tmp_174*tmp_26;
      real_t tmp_177 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_178 = tmp_177*tmp_28;
      real_t tmp_179 = tmp_177*tmp_35;
      real_t tmp_180 = tmp_174*tmp_37;
      real_t tmp_181 = tmp_177*tmp_39;
      real_t tmp_182 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_183 = tmp_182*tmp_41;
      real_t tmp_184 = tmp_182*tmp_48;
      real_t tmp_185 = tmp_182*tmp_50;
      real_t tmp_186 = 0.020848748529055869*tmp_55;
      real_t tmp_187 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_188 = tmp_187*tmp_8;
      real_t tmp_189 = tmp_187*tmp_26;
      real_t tmp_190 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_191 = tmp_190*tmp_28;
      real_t tmp_192 = tmp_190*tmp_35;
      real_t tmp_193 = tmp_187*tmp_37;
      real_t tmp_194 = tmp_190*tmp_39;
      real_t tmp_195 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_196 = tmp_195*tmp_41;
      real_t tmp_197 = tmp_195*tmp_48;
      real_t tmp_198 = tmp_195*tmp_50;
      real_t tmp_199 = 0.019202922745021479*tmp_55;
      real_t tmp_200 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_201 = tmp_200*tmp_8;
      real_t tmp_202 = tmp_200*tmp_26;
      real_t tmp_203 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_204 = tmp_203*tmp_28;
      real_t tmp_205 = tmp_203*tmp_35;
      real_t tmp_206 = tmp_200*tmp_37;
      real_t tmp_207 = tmp_203*tmp_39;
      real_t tmp_208 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_209 = tmp_208*tmp_41;
      real_t tmp_210 = tmp_208*tmp_48;
      real_t tmp_211 = tmp_208*tmp_50;
      real_t tmp_212 = 0.042507265838595799*tmp_55;
      real_t tmp_213 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_214 = tmp_213*tmp_8;
      real_t tmp_215 = tmp_213*tmp_26;
      real_t tmp_216 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_217 = tmp_216*tmp_28;
      real_t tmp_218 = tmp_216*tmp_35;
      real_t tmp_219 = tmp_213*tmp_37;
      real_t tmp_220 = tmp_216*tmp_39;
      real_t tmp_221 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_222 = tmp_221*tmp_41;
      real_t tmp_223 = tmp_221*tmp_48;
      real_t tmp_224 = tmp_221*tmp_50;
      real_t tmp_225 = 0.020848748529055869*tmp_55;
      real_t tmp_226 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_227 = tmp_226*tmp_8;
      real_t tmp_228 = tmp_226*tmp_26;
      real_t tmp_229 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_230 = tmp_229*tmp_28;
      real_t tmp_231 = tmp_229*tmp_35;
      real_t tmp_232 = tmp_226*tmp_37;
      real_t tmp_233 = tmp_229*tmp_39;
      real_t tmp_234 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_235 = tmp_234*tmp_41;
      real_t tmp_236 = tmp_234*tmp_48;
      real_t tmp_237 = tmp_234*tmp_50;
      real_t tmp_238 = 0.0068572537431980923*tmp_55;
      real_t tmp_239 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_240 = tmp_239*tmp_8;
      real_t tmp_241 = tmp_239*tmp_26;
      real_t tmp_242 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_243 = tmp_242*tmp_28;
      real_t tmp_244 = tmp_242*tmp_35;
      real_t tmp_245 = tmp_239*tmp_37;
      real_t tmp_246 = tmp_242*tmp_39;
      real_t tmp_247 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_248 = tmp_247*tmp_41;
      real_t tmp_249 = tmp_247*tmp_48;
      real_t tmp_250 = tmp_247*tmp_50;
      real_t tmp_251 = 0.037198804536718075*tmp_55;
      real_t tmp_252 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_8;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = 0.042507265838595799*tmp_55;
      real_t tmp_265 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_266 = tmp_265*tmp_8;
      real_t tmp_267 = tmp_26*tmp_265;
      real_t tmp_268 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_269 = tmp_268*tmp_28;
      real_t tmp_270 = tmp_268*tmp_35;
      real_t tmp_271 = tmp_265*tmp_37;
      real_t tmp_272 = tmp_268*tmp_39;
      real_t tmp_273 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_274 = tmp_273*tmp_41;
      real_t tmp_275 = tmp_273*tmp_48;
      real_t tmp_276 = tmp_273*tmp_50;
      real_t tmp_277 = 0.0068572537431980923*tmp_55;
      real_t tmp_278 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_279 = tmp_278*tmp_8;
      real_t tmp_280 = tmp_26*tmp_278;
      real_t tmp_281 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_282 = tmp_28*tmp_281;
      real_t tmp_283 = tmp_281*tmp_35;
      real_t tmp_284 = tmp_278*tmp_37;
      real_t tmp_285 = tmp_281*tmp_39;
      real_t tmp_286 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_287 = tmp_286*tmp_41;
      real_t tmp_288 = tmp_286*tmp_48;
      real_t tmp_289 = tmp_286*tmp_50;
      real_t tmp_290 = 0.037198804536718075*tmp_55;
      real_t tmp_291 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_292 = tmp_291*tmp_8;
      real_t tmp_293 = tmp_26*tmp_291;
      real_t tmp_294 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_295 = tmp_28*tmp_294;
      real_t tmp_296 = tmp_294*tmp_35;
      real_t tmp_297 = tmp_291*tmp_37;
      real_t tmp_298 = tmp_294*tmp_39;
      real_t tmp_299 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_300 = tmp_299*tmp_41;
      real_t tmp_301 = tmp_299*tmp_48;
      real_t tmp_302 = tmp_299*tmp_50;
      real_t tmp_303 = 0.042507265838595799*tmp_55;
      real_t tmp_304 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_305 = tmp_304*tmp_8;
      real_t tmp_306 = tmp_26*tmp_304;
      real_t tmp_307 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_308 = tmp_28*tmp_307;
      real_t tmp_309 = tmp_307*tmp_35;
      real_t tmp_310 = tmp_304*tmp_37;
      real_t tmp_311 = tmp_307*tmp_39;
      real_t tmp_312 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_313 = tmp_312*tmp_41;
      real_t tmp_314 = tmp_312*tmp_48;
      real_t tmp_315 = tmp_312*tmp_50;
      real_t tmp_316 = 0.019202922745021479*tmp_55;
      real_t a_0_0 = tmp_108*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_105 - tmp_106 - tmp_107 - tmp_97 - tmp_98 + 1) + tmp_121*(-tmp_110 - tmp_111 - tmp_113 - tmp_114 - tmp_115 - tmp_116 - tmp_118 - tmp_119 - tmp_120 + 1) + tmp_134*(-tmp_123 - tmp_124 - tmp_126 - tmp_127 - tmp_128 - tmp_129 - tmp_131 - tmp_132 - tmp_133 + 1) + tmp_147*(-tmp_136 - tmp_137 - tmp_139 - tmp_140 - tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 + 1) + tmp_160*(-tmp_149 - tmp_150 - tmp_152 - tmp_153 - tmp_154 - tmp_155 - tmp_157 - tmp_158 - tmp_159 + 1) + tmp_173*(-tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 - tmp_168 - tmp_170 - tmp_171 - tmp_172 + 1) + tmp_186*(-tmp_175 - tmp_176 - tmp_178 - tmp_179 - tmp_180 - tmp_181 - tmp_183 - tmp_184 - tmp_185 + 1) + tmp_199*(-tmp_188 - tmp_189 - tmp_191 - tmp_192 - tmp_193 - tmp_194 - tmp_196 - tmp_197 - tmp_198 + 1) + tmp_212*(-tmp_201 - tmp_202 - tmp_204 - tmp_205 - tmp_206 - tmp_207 - tmp_209 - tmp_210 - tmp_211 + 1) + tmp_225*(-tmp_214 - tmp_215 - tmp_217 - tmp_218 - tmp_219 - tmp_220 - tmp_222 - tmp_223 - tmp_224 + 1) + tmp_238*(-tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 - tmp_233 - tmp_235 - tmp_236 - tmp_237 + 1) + tmp_251*(-tmp_240 - tmp_241 - tmp_243 - tmp_244 - tmp_245 - tmp_246 - tmp_248 - tmp_249 - tmp_250 + 1) + tmp_264*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_277*(-tmp_266 - tmp_267 - tmp_269 - tmp_270 - tmp_271 - tmp_272 - tmp_274 - tmp_275 - tmp_276 + 1) + tmp_290*(-tmp_279 - tmp_280 - tmp_282 - tmp_283 - tmp_284 - tmp_285 - tmp_287 - tmp_288 - tmp_289 + 1) + tmp_303*(-tmp_292 - tmp_293 - tmp_295 - tmp_296 - tmp_297 - tmp_298 - tmp_300 - tmp_301 - tmp_302 + 1) + tmp_316*(-tmp_305 - tmp_306 - tmp_308 - tmp_309 - tmp_310 - tmp_311 - tmp_313 - tmp_314 - tmp_315 + 1) + tmp_56*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_69*(-tmp_58 - tmp_59 - tmp_61 - tmp_62 - tmp_63 - tmp_64 - tmp_66 - tmp_67 - tmp_68 + 1) + tmp_82*(-tmp_71 - tmp_72 - tmp_74 - tmp_75 - tmp_76 - tmp_77 - tmp_79 - tmp_80 - tmp_81 + 1) + tmp_95*(-tmp_84 - tmp_85 - tmp_87 - tmp_88 - tmp_89 - tmp_90 - tmp_92 - tmp_93 - tmp_94 + 1);
      real_t a_1_0 = tmp_108*(tmp_102 + tmp_103 + tmp_107) + tmp_121*(tmp_115 + tmp_116 + tmp_120) + tmp_134*(tmp_128 + tmp_129 + tmp_133) + tmp_147*(tmp_141 + tmp_142 + tmp_146) + tmp_160*(tmp_154 + tmp_155 + tmp_159) + tmp_173*(tmp_167 + tmp_168 + tmp_172) + tmp_186*(tmp_180 + tmp_181 + tmp_185) + tmp_199*(tmp_193 + tmp_194 + tmp_198) + tmp_212*(tmp_206 + tmp_207 + tmp_211) + tmp_225*(tmp_219 + tmp_220 + tmp_224) + tmp_238*(tmp_232 + tmp_233 + tmp_237) + tmp_251*(tmp_245 + tmp_246 + tmp_250) + tmp_264*(tmp_258 + tmp_259 + tmp_263) + tmp_277*(tmp_271 + tmp_272 + tmp_276) + tmp_290*(tmp_284 + tmp_285 + tmp_289) + tmp_303*(tmp_297 + tmp_298 + tmp_302) + tmp_316*(tmp_310 + tmp_311 + tmp_315) + tmp_56*(tmp_38 + tmp_40 + tmp_51) + tmp_69*(tmp_63 + tmp_64 + tmp_68) + tmp_82*(tmp_76 + tmp_77 + tmp_81) + tmp_95*(tmp_89 + tmp_90 + tmp_94);
      real_t a_2_0 = tmp_108*(tmp_101 + tmp_106 + tmp_98) + tmp_121*(tmp_111 + tmp_114 + tmp_119) + tmp_134*(tmp_124 + tmp_127 + tmp_132) + tmp_147*(tmp_137 + tmp_140 + tmp_145) + tmp_160*(tmp_150 + tmp_153 + tmp_158) + tmp_173*(tmp_163 + tmp_166 + tmp_171) + tmp_186*(tmp_176 + tmp_179 + tmp_184) + tmp_199*(tmp_189 + tmp_192 + tmp_197) + tmp_212*(tmp_202 + tmp_205 + tmp_210) + tmp_225*(tmp_215 + tmp_218 + tmp_223) + tmp_238*(tmp_228 + tmp_231 + tmp_236) + tmp_251*(tmp_241 + tmp_244 + tmp_249) + tmp_264*(tmp_254 + tmp_257 + tmp_262) + tmp_277*(tmp_267 + tmp_270 + tmp_275) + tmp_290*(tmp_280 + tmp_283 + tmp_288) + tmp_303*(tmp_293 + tmp_296 + tmp_301) + tmp_316*(tmp_306 + tmp_309 + tmp_314) + tmp_56*(tmp_27 + tmp_36 + tmp_49) + tmp_69*(tmp_59 + tmp_62 + tmp_67) + tmp_82*(tmp_72 + tmp_75 + tmp_80) + tmp_95*(tmp_85 + tmp_88 + tmp_93);
      real_t a_3_0 = tmp_108*(tmp_100 + tmp_105 + tmp_97) + tmp_121*(tmp_110 + tmp_113 + tmp_118) + tmp_134*(tmp_123 + tmp_126 + tmp_131) + tmp_147*(tmp_136 + tmp_139 + tmp_144) + tmp_160*(tmp_149 + tmp_152 + tmp_157) + tmp_173*(tmp_162 + tmp_165 + tmp_170) + tmp_186*(tmp_175 + tmp_178 + tmp_183) + tmp_199*(tmp_188 + tmp_191 + tmp_196) + tmp_212*(tmp_201 + tmp_204 + tmp_209) + tmp_225*(tmp_214 + tmp_217 + tmp_222) + tmp_238*(tmp_227 + tmp_230 + tmp_235) + tmp_251*(tmp_240 + tmp_243 + tmp_248) + tmp_264*(tmp_253 + tmp_256 + tmp_261) + tmp_277*(tmp_266 + tmp_269 + tmp_274) + tmp_290*(tmp_279 + tmp_282 + tmp_287) + tmp_303*(tmp_292 + tmp_295 + tmp_300) + tmp_316*(tmp_305 + tmp_308 + tmp_313) + tmp_56*(tmp_25 + tmp_34 + tmp_47) + tmp_69*(tmp_58 + tmp_61 + tmp_66) + tmp_82*(tmp_71 + tmp_74 + tmp_79) + tmp_95*(tmp_84 + tmp_87 + tmp_92);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }




void integrateFacetCoupling3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementInner,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementOuter,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                                        Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = p_affine_1_2 + tmp_9;
      real_t tmp_13 = tmp_12*tmp_5;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = p_affine_2_2 + tmp_9;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_1*tmp_11;
      real_t tmp_18 = tmp_12*tmp_14;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_4 - tmp_10*tmp_7 + tmp_11*tmp_13 + tmp_14*tmp_16 - tmp_15*tmp_17 - tmp_18*tmp_3);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_9;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_8;
      real_t tmp_26 = tmp_14*tmp_6 - tmp_17;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_15 + tmp_13;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_10 - tmp_18;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_11*tmp_5 - tmp_14*tmp_3;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_12*tmp_3 + tmp_16;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = tmp_10*tmp_3 - tmp_11*tmp_15;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_53 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_54 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_55 = 0.5*p_affine_13_2*std::pow((std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)*std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)) + (std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)*std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)) + (std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)*std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)), 1.0/2.0);
      real_t tmp_56 = 0.0068572537431980923*tmp_55;
      real_t tmp_57 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_58 = tmp_57*tmp_8;
      real_t tmp_59 = tmp_26*tmp_57;
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_61 = tmp_28*tmp_60;
      real_t tmp_62 = tmp_35*tmp_60;
      real_t tmp_63 = tmp_37*tmp_57;
      real_t tmp_64 = tmp_39*tmp_60;
      real_t tmp_65 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_66 = tmp_41*tmp_65;
      real_t tmp_67 = tmp_48*tmp_65;
      real_t tmp_68 = tmp_50*tmp_65;
      real_t tmp_69 = 0.037198804536718075*tmp_55;
      real_t tmp_70 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_71 = tmp_70*tmp_8;
      real_t tmp_72 = tmp_26*tmp_70;
      real_t tmp_73 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_74 = tmp_28*tmp_73;
      real_t tmp_75 = tmp_35*tmp_73;
      real_t tmp_76 = tmp_37*tmp_70;
      real_t tmp_77 = tmp_39*tmp_73;
      real_t tmp_78 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_79 = tmp_41*tmp_78;
      real_t tmp_80 = tmp_48*tmp_78;
      real_t tmp_81 = tmp_50*tmp_78;
      real_t tmp_82 = 0.020848748529055869*tmp_55;
      real_t tmp_83 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_84 = tmp_8*tmp_83;
      real_t tmp_85 = tmp_26*tmp_83;
      real_t tmp_86 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_87 = tmp_28*tmp_86;
      real_t tmp_88 = tmp_35*tmp_86;
      real_t tmp_89 = tmp_37*tmp_83;
      real_t tmp_90 = tmp_39*tmp_86;
      real_t tmp_91 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_92 = tmp_41*tmp_91;
      real_t tmp_93 = tmp_48*tmp_91;
      real_t tmp_94 = tmp_50*tmp_91;
      real_t tmp_95 = 0.019202922745021479*tmp_55;
      real_t tmp_96 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_97 = tmp_8*tmp_96;
      real_t tmp_98 = tmp_26*tmp_96;
      real_t tmp_99 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_100 = tmp_28*tmp_99;
      real_t tmp_101 = tmp_35*tmp_99;
      real_t tmp_102 = tmp_37*tmp_96;
      real_t tmp_103 = tmp_39*tmp_99;
      real_t tmp_104 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_105 = tmp_104*tmp_41;
      real_t tmp_106 = tmp_104*tmp_48;
      real_t tmp_107 = tmp_104*tmp_50;
      real_t tmp_108 = 0.020848748529055869*tmp_55;
      real_t tmp_109 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_110 = tmp_109*tmp_8;
      real_t tmp_111 = tmp_109*tmp_26;
      real_t tmp_112 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_113 = tmp_112*tmp_28;
      real_t tmp_114 = tmp_112*tmp_35;
      real_t tmp_115 = tmp_109*tmp_37;
      real_t tmp_116 = tmp_112*tmp_39;
      real_t tmp_117 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_118 = tmp_117*tmp_41;
      real_t tmp_119 = tmp_117*tmp_48;
      real_t tmp_120 = tmp_117*tmp_50;
      real_t tmp_121 = 0.019202922745021479*tmp_55;
      real_t tmp_122 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_123 = tmp_122*tmp_8;
      real_t tmp_124 = tmp_122*tmp_26;
      real_t tmp_125 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_126 = tmp_125*tmp_28;
      real_t tmp_127 = tmp_125*tmp_35;
      real_t tmp_128 = tmp_122*tmp_37;
      real_t tmp_129 = tmp_125*tmp_39;
      real_t tmp_130 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_131 = tmp_130*tmp_41;
      real_t tmp_132 = tmp_130*tmp_48;
      real_t tmp_133 = tmp_130*tmp_50;
      real_t tmp_134 = 0.020848748529055869*tmp_55;
      real_t tmp_135 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_136 = tmp_135*tmp_8;
      real_t tmp_137 = tmp_135*tmp_26;
      real_t tmp_138 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_139 = tmp_138*tmp_28;
      real_t tmp_140 = tmp_138*tmp_35;
      real_t tmp_141 = tmp_135*tmp_37;
      real_t tmp_142 = tmp_138*tmp_39;
      real_t tmp_143 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_144 = tmp_143*tmp_41;
      real_t tmp_145 = tmp_143*tmp_48;
      real_t tmp_146 = tmp_143*tmp_50;
      real_t tmp_147 = 0.019202922745021479*tmp_55;
      real_t tmp_148 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_149 = tmp_148*tmp_8;
      real_t tmp_150 = tmp_148*tmp_26;
      real_t tmp_151 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_152 = tmp_151*tmp_28;
      real_t tmp_153 = tmp_151*tmp_35;
      real_t tmp_154 = tmp_148*tmp_37;
      real_t tmp_155 = tmp_151*tmp_39;
      real_t tmp_156 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_157 = tmp_156*tmp_41;
      real_t tmp_158 = tmp_156*tmp_48;
      real_t tmp_159 = tmp_156*tmp_50;
      real_t tmp_160 = 0.020848748529055869*tmp_55;
      real_t tmp_161 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_162 = tmp_161*tmp_8;
      real_t tmp_163 = tmp_161*tmp_26;
      real_t tmp_164 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_165 = tmp_164*tmp_28;
      real_t tmp_166 = tmp_164*tmp_35;
      real_t tmp_167 = tmp_161*tmp_37;
      real_t tmp_168 = tmp_164*tmp_39;
      real_t tmp_169 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_170 = tmp_169*tmp_41;
      real_t tmp_171 = tmp_169*tmp_48;
      real_t tmp_172 = tmp_169*tmp_50;
      real_t tmp_173 = 0.019202922745021479*tmp_55;
      real_t tmp_174 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_175 = tmp_174*tmp_8;
      real_t tmp_176 = tmp_174*tmp_26;
      real_t tmp_177 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_178 = tmp_177*tmp_28;
      real_t tmp_179 = tmp_177*tmp_35;
      real_t tmp_180 = tmp_174*tmp_37;
      real_t tmp_181 = tmp_177*tmp_39;
      real_t tmp_182 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_183 = tmp_182*tmp_41;
      real_t tmp_184 = tmp_182*tmp_48;
      real_t tmp_185 = tmp_182*tmp_50;
      real_t tmp_186 = 0.020848748529055869*tmp_55;
      real_t tmp_187 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_188 = tmp_187*tmp_8;
      real_t tmp_189 = tmp_187*tmp_26;
      real_t tmp_190 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_191 = tmp_190*tmp_28;
      real_t tmp_192 = tmp_190*tmp_35;
      real_t tmp_193 = tmp_187*tmp_37;
      real_t tmp_194 = tmp_190*tmp_39;
      real_t tmp_195 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_196 = tmp_195*tmp_41;
      real_t tmp_197 = tmp_195*tmp_48;
      real_t tmp_198 = tmp_195*tmp_50;
      real_t tmp_199 = 0.019202922745021479*tmp_55;
      real_t tmp_200 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_201 = tmp_200*tmp_8;
      real_t tmp_202 = tmp_200*tmp_26;
      real_t tmp_203 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_204 = tmp_203*tmp_28;
      real_t tmp_205 = tmp_203*tmp_35;
      real_t tmp_206 = tmp_200*tmp_37;
      real_t tmp_207 = tmp_203*tmp_39;
      real_t tmp_208 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_209 = tmp_208*tmp_41;
      real_t tmp_210 = tmp_208*tmp_48;
      real_t tmp_211 = tmp_208*tmp_50;
      real_t tmp_212 = 0.042507265838595799*tmp_55;
      real_t tmp_213 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_214 = tmp_213*tmp_8;
      real_t tmp_215 = tmp_213*tmp_26;
      real_t tmp_216 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_217 = tmp_216*tmp_28;
      real_t tmp_218 = tmp_216*tmp_35;
      real_t tmp_219 = tmp_213*tmp_37;
      real_t tmp_220 = tmp_216*tmp_39;
      real_t tmp_221 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_222 = tmp_221*tmp_41;
      real_t tmp_223 = tmp_221*tmp_48;
      real_t tmp_224 = tmp_221*tmp_50;
      real_t tmp_225 = 0.020848748529055869*tmp_55;
      real_t tmp_226 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_227 = tmp_226*tmp_8;
      real_t tmp_228 = tmp_226*tmp_26;
      real_t tmp_229 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_230 = tmp_229*tmp_28;
      real_t tmp_231 = tmp_229*tmp_35;
      real_t tmp_232 = tmp_226*tmp_37;
      real_t tmp_233 = tmp_229*tmp_39;
      real_t tmp_234 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_235 = tmp_234*tmp_41;
      real_t tmp_236 = tmp_234*tmp_48;
      real_t tmp_237 = tmp_234*tmp_50;
      real_t tmp_238 = 0.0068572537431980923*tmp_55;
      real_t tmp_239 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_240 = tmp_239*tmp_8;
      real_t tmp_241 = tmp_239*tmp_26;
      real_t tmp_242 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_243 = tmp_242*tmp_28;
      real_t tmp_244 = tmp_242*tmp_35;
      real_t tmp_245 = tmp_239*tmp_37;
      real_t tmp_246 = tmp_242*tmp_39;
      real_t tmp_247 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_248 = tmp_247*tmp_41;
      real_t tmp_249 = tmp_247*tmp_48;
      real_t tmp_250 = tmp_247*tmp_50;
      real_t tmp_251 = 0.037198804536718075*tmp_55;
      real_t tmp_252 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_8;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = 0.042507265838595799*tmp_55;
      real_t tmp_265 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_266 = tmp_265*tmp_8;
      real_t tmp_267 = tmp_26*tmp_265;
      real_t tmp_268 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_269 = tmp_268*tmp_28;
      real_t tmp_270 = tmp_268*tmp_35;
      real_t tmp_271 = tmp_265*tmp_37;
      real_t tmp_272 = tmp_268*tmp_39;
      real_t tmp_273 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_274 = tmp_273*tmp_41;
      real_t tmp_275 = tmp_273*tmp_48;
      real_t tmp_276 = tmp_273*tmp_50;
      real_t tmp_277 = 0.0068572537431980923*tmp_55;
      real_t tmp_278 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_279 = tmp_278*tmp_8;
      real_t tmp_280 = tmp_26*tmp_278;
      real_t tmp_281 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_282 = tmp_28*tmp_281;
      real_t tmp_283 = tmp_281*tmp_35;
      real_t tmp_284 = tmp_278*tmp_37;
      real_t tmp_285 = tmp_281*tmp_39;
      real_t tmp_286 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_287 = tmp_286*tmp_41;
      real_t tmp_288 = tmp_286*tmp_48;
      real_t tmp_289 = tmp_286*tmp_50;
      real_t tmp_290 = 0.037198804536718075*tmp_55;
      real_t tmp_291 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_292 = tmp_291*tmp_8;
      real_t tmp_293 = tmp_26*tmp_291;
      real_t tmp_294 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_295 = tmp_28*tmp_294;
      real_t tmp_296 = tmp_294*tmp_35;
      real_t tmp_297 = tmp_291*tmp_37;
      real_t tmp_298 = tmp_294*tmp_39;
      real_t tmp_299 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_300 = tmp_299*tmp_41;
      real_t tmp_301 = tmp_299*tmp_48;
      real_t tmp_302 = tmp_299*tmp_50;
      real_t tmp_303 = 0.042507265838595799*tmp_55;
      real_t tmp_304 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_305 = tmp_304*tmp_8;
      real_t tmp_306 = tmp_26*tmp_304;
      real_t tmp_307 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_308 = tmp_28*tmp_307;
      real_t tmp_309 = tmp_307*tmp_35;
      real_t tmp_310 = tmp_304*tmp_37;
      real_t tmp_311 = tmp_307*tmp_39;
      real_t tmp_312 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_313 = tmp_312*tmp_41;
      real_t tmp_314 = tmp_312*tmp_48;
      real_t tmp_315 = tmp_312*tmp_50;
      real_t tmp_316 = 0.019202922745021479*tmp_55;
      real_t a_0_0 = tmp_108*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_105 - tmp_106 - tmp_107 - tmp_97 - tmp_98 + 1) + tmp_121*(-tmp_110 - tmp_111 - tmp_113 - tmp_114 - tmp_115 - tmp_116 - tmp_118 - tmp_119 - tmp_120 + 1) + tmp_134*(-tmp_123 - tmp_124 - tmp_126 - tmp_127 - tmp_128 - tmp_129 - tmp_131 - tmp_132 - tmp_133 + 1) + tmp_147*(-tmp_136 - tmp_137 - tmp_139 - tmp_140 - tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 + 1) + tmp_160*(-tmp_149 - tmp_150 - tmp_152 - tmp_153 - tmp_154 - tmp_155 - tmp_157 - tmp_158 - tmp_159 + 1) + tmp_173*(-tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 - tmp_168 - tmp_170 - tmp_171 - tmp_172 + 1) + tmp_186*(-tmp_175 - tmp_176 - tmp_178 - tmp_179 - tmp_180 - tmp_181 - tmp_183 - tmp_184 - tmp_185 + 1) + tmp_199*(-tmp_188 - tmp_189 - tmp_191 - tmp_192 - tmp_193 - tmp_194 - tmp_196 - tmp_197 - tmp_198 + 1) + tmp_212*(-tmp_201 - tmp_202 - tmp_204 - tmp_205 - tmp_206 - tmp_207 - tmp_209 - tmp_210 - tmp_211 + 1) + tmp_225*(-tmp_214 - tmp_215 - tmp_217 - tmp_218 - tmp_219 - tmp_220 - tmp_222 - tmp_223 - tmp_224 + 1) + tmp_238*(-tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 - tmp_233 - tmp_235 - tmp_236 - tmp_237 + 1) + tmp_251*(-tmp_240 - tmp_241 - tmp_243 - tmp_244 - tmp_245 - tmp_246 - tmp_248 - tmp_249 - tmp_250 + 1) + tmp_264*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_277*(-tmp_266 - tmp_267 - tmp_269 - tmp_270 - tmp_271 - tmp_272 - tmp_274 - tmp_275 - tmp_276 + 1) + tmp_290*(-tmp_279 - tmp_280 - tmp_282 - tmp_283 - tmp_284 - tmp_285 - tmp_287 - tmp_288 - tmp_289 + 1) + tmp_303*(-tmp_292 - tmp_293 - tmp_295 - tmp_296 - tmp_297 - tmp_298 - tmp_300 - tmp_301 - tmp_302 + 1) + tmp_316*(-tmp_305 - tmp_306 - tmp_308 - tmp_309 - tmp_310 - tmp_311 - tmp_313 - tmp_314 - tmp_315 + 1) + tmp_56*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_69*(-tmp_58 - tmp_59 - tmp_61 - tmp_62 - tmp_63 - tmp_64 - tmp_66 - tmp_67 - tmp_68 + 1) + tmp_82*(-tmp_71 - tmp_72 - tmp_74 - tmp_75 - tmp_76 - tmp_77 - tmp_79 - tmp_80 - tmp_81 + 1) + tmp_95*(-tmp_84 - tmp_85 - tmp_87 - tmp_88 - tmp_89 - tmp_90 - tmp_92 - tmp_93 - tmp_94 + 1);
      real_t a_1_0 = tmp_108*(tmp_102 + tmp_103 + tmp_107) + tmp_121*(tmp_115 + tmp_116 + tmp_120) + tmp_134*(tmp_128 + tmp_129 + tmp_133) + tmp_147*(tmp_141 + tmp_142 + tmp_146) + tmp_160*(tmp_154 + tmp_155 + tmp_159) + tmp_173*(tmp_167 + tmp_168 + tmp_172) + tmp_186*(tmp_180 + tmp_181 + tmp_185) + tmp_199*(tmp_193 + tmp_194 + tmp_198) + tmp_212*(tmp_206 + tmp_207 + tmp_211) + tmp_225*(tmp_219 + tmp_220 + tmp_224) + tmp_238*(tmp_232 + tmp_233 + tmp_237) + tmp_251*(tmp_245 + tmp_246 + tmp_250) + tmp_264*(tmp_258 + tmp_259 + tmp_263) + tmp_277*(tmp_271 + tmp_272 + tmp_276) + tmp_290*(tmp_284 + tmp_285 + tmp_289) + tmp_303*(tmp_297 + tmp_298 + tmp_302) + tmp_316*(tmp_310 + tmp_311 + tmp_315) + tmp_56*(tmp_38 + tmp_40 + tmp_51) + tmp_69*(tmp_63 + tmp_64 + tmp_68) + tmp_82*(tmp_76 + tmp_77 + tmp_81) + tmp_95*(tmp_89 + tmp_90 + tmp_94);
      real_t a_2_0 = tmp_108*(tmp_101 + tmp_106 + tmp_98) + tmp_121*(tmp_111 + tmp_114 + tmp_119) + tmp_134*(tmp_124 + tmp_127 + tmp_132) + tmp_147*(tmp_137 + tmp_140 + tmp_145) + tmp_160*(tmp_150 + tmp_153 + tmp_158) + tmp_173*(tmp_163 + tmp_166 + tmp_171) + tmp_186*(tmp_176 + tmp_179 + tmp_184) + tmp_199*(tmp_189 + tmp_192 + tmp_197) + tmp_212*(tmp_202 + tmp_205 + tmp_210) + tmp_225*(tmp_215 + tmp_218 + tmp_223) + tmp_238*(tmp_228 + tmp_231 + tmp_236) + tmp_251*(tmp_241 + tmp_244 + tmp_249) + tmp_264*(tmp_254 + tmp_257 + tmp_262) + tmp_277*(tmp_267 + tmp_270 + tmp_275) + tmp_290*(tmp_280 + tmp_283 + tmp_288) + tmp_303*(tmp_293 + tmp_296 + tmp_301) + tmp_316*(tmp_306 + tmp_309 + tmp_314) + tmp_56*(tmp_27 + tmp_36 + tmp_49) + tmp_69*(tmp_59 + tmp_62 + tmp_67) + tmp_82*(tmp_72 + tmp_75 + tmp_80) + tmp_95*(tmp_85 + tmp_88 + tmp_93);
      real_t a_3_0 = tmp_108*(tmp_100 + tmp_105 + tmp_97) + tmp_121*(tmp_110 + tmp_113 + tmp_118) + tmp_134*(tmp_123 + tmp_126 + tmp_131) + tmp_147*(tmp_136 + tmp_139 + tmp_144) + tmp_160*(tmp_149 + tmp_152 + tmp_157) + tmp_173*(tmp_162 + tmp_165 + tmp_170) + tmp_186*(tmp_175 + tmp_178 + tmp_183) + tmp_199*(tmp_188 + tmp_191 + tmp_196) + tmp_212*(tmp_201 + tmp_204 + tmp_209) + tmp_225*(tmp_214 + tmp_217 + tmp_222) + tmp_238*(tmp_227 + tmp_230 + tmp_235) + tmp_251*(tmp_240 + tmp_243 + tmp_248) + tmp_264*(tmp_253 + tmp_256 + tmp_261) + tmp_277*(tmp_266 + tmp_269 + tmp_274) + tmp_290*(tmp_279 + tmp_282 + tmp_287) + tmp_303*(tmp_292 + tmp_295 + tmp_300) + tmp_316*(tmp_305 + tmp_308 + tmp_313) + tmp_56*(tmp_25 + tmp_34 + tmp_47) + tmp_69*(tmp_58 + tmp_61 + tmp_66) + tmp_82*(tmp_71 + tmp_74 + tmp_79) + tmp_95*(tmp_84 + tmp_87 + tmp_92);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
    const Eigen::Matrix< real_t, 3, 1 >&,
    const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
    Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = p_affine_1_2 + tmp_9;
      real_t tmp_13 = tmp_12*tmp_5;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = p_affine_2_2 + tmp_9;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_1*tmp_11;
      real_t tmp_18 = tmp_12*tmp_14;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_4 - tmp_10*tmp_7 + tmp_11*tmp_13 + tmp_14*tmp_16 - tmp_15*tmp_17 - tmp_18*tmp_3);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_9;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_8;
      real_t tmp_26 = tmp_14*tmp_6 - tmp_17;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_15 + tmp_13;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_10 - tmp_18;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_11*tmp_5 - tmp_14*tmp_3;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_12*tmp_3 + tmp_16;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = tmp_10*tmp_3 - tmp_11*tmp_15;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_53 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_54 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_55 = 1.0*p_affine_13_2*std::pow((std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)*std::abs(tmp_22*tmp_52 - tmp_31*tmp_54)) + (std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)*std::abs(tmp_22*tmp_53 - tmp_44*tmp_54)) + (std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)*std::abs(tmp_31*tmp_53 - tmp_44*tmp_52)), 1.0/2.0);
      real_t tmp_56 = 0.0068572537431980923*tmp_55;
      real_t tmp_57 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_58 = tmp_57*tmp_8;
      real_t tmp_59 = tmp_26*tmp_57;
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_61 = tmp_28*tmp_60;
      real_t tmp_62 = tmp_35*tmp_60;
      real_t tmp_63 = tmp_37*tmp_57;
      real_t tmp_64 = tmp_39*tmp_60;
      real_t tmp_65 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_66 = tmp_41*tmp_65;
      real_t tmp_67 = tmp_48*tmp_65;
      real_t tmp_68 = tmp_50*tmp_65;
      real_t tmp_69 = 0.037198804536718075*tmp_55;
      real_t tmp_70 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_71 = tmp_70*tmp_8;
      real_t tmp_72 = tmp_26*tmp_70;
      real_t tmp_73 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_74 = tmp_28*tmp_73;
      real_t tmp_75 = tmp_35*tmp_73;
      real_t tmp_76 = tmp_37*tmp_70;
      real_t tmp_77 = tmp_39*tmp_73;
      real_t tmp_78 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_79 = tmp_41*tmp_78;
      real_t tmp_80 = tmp_48*tmp_78;
      real_t tmp_81 = tmp_50*tmp_78;
      real_t tmp_82 = 0.020848748529055869*tmp_55;
      real_t tmp_83 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_84 = tmp_8*tmp_83;
      real_t tmp_85 = tmp_26*tmp_83;
      real_t tmp_86 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_87 = tmp_28*tmp_86;
      real_t tmp_88 = tmp_35*tmp_86;
      real_t tmp_89 = tmp_37*tmp_83;
      real_t tmp_90 = tmp_39*tmp_86;
      real_t tmp_91 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_92 = tmp_41*tmp_91;
      real_t tmp_93 = tmp_48*tmp_91;
      real_t tmp_94 = tmp_50*tmp_91;
      real_t tmp_95 = 0.019202922745021479*tmp_55;
      real_t tmp_96 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_97 = tmp_8*tmp_96;
      real_t tmp_98 = tmp_26*tmp_96;
      real_t tmp_99 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_100 = tmp_28*tmp_99;
      real_t tmp_101 = tmp_35*tmp_99;
      real_t tmp_102 = tmp_37*tmp_96;
      real_t tmp_103 = tmp_39*tmp_99;
      real_t tmp_104 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_105 = tmp_104*tmp_41;
      real_t tmp_106 = tmp_104*tmp_48;
      real_t tmp_107 = tmp_104*tmp_50;
      real_t tmp_108 = 0.020848748529055869*tmp_55;
      real_t tmp_109 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_110 = tmp_109*tmp_8;
      real_t tmp_111 = tmp_109*tmp_26;
      real_t tmp_112 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_113 = tmp_112*tmp_28;
      real_t tmp_114 = tmp_112*tmp_35;
      real_t tmp_115 = tmp_109*tmp_37;
      real_t tmp_116 = tmp_112*tmp_39;
      real_t tmp_117 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_118 = tmp_117*tmp_41;
      real_t tmp_119 = tmp_117*tmp_48;
      real_t tmp_120 = tmp_117*tmp_50;
      real_t tmp_121 = 0.019202922745021479*tmp_55;
      real_t tmp_122 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_123 = tmp_122*tmp_8;
      real_t tmp_124 = tmp_122*tmp_26;
      real_t tmp_125 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_126 = tmp_125*tmp_28;
      real_t tmp_127 = tmp_125*tmp_35;
      real_t tmp_128 = tmp_122*tmp_37;
      real_t tmp_129 = tmp_125*tmp_39;
      real_t tmp_130 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_131 = tmp_130*tmp_41;
      real_t tmp_132 = tmp_130*tmp_48;
      real_t tmp_133 = tmp_130*tmp_50;
      real_t tmp_134 = 0.020848748529055869*tmp_55;
      real_t tmp_135 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_136 = tmp_135*tmp_8;
      real_t tmp_137 = tmp_135*tmp_26;
      real_t tmp_138 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_139 = tmp_138*tmp_28;
      real_t tmp_140 = tmp_138*tmp_35;
      real_t tmp_141 = tmp_135*tmp_37;
      real_t tmp_142 = tmp_138*tmp_39;
      real_t tmp_143 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_144 = tmp_143*tmp_41;
      real_t tmp_145 = tmp_143*tmp_48;
      real_t tmp_146 = tmp_143*tmp_50;
      real_t tmp_147 = 0.019202922745021479*tmp_55;
      real_t tmp_148 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_149 = tmp_148*tmp_8;
      real_t tmp_150 = tmp_148*tmp_26;
      real_t tmp_151 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_152 = tmp_151*tmp_28;
      real_t tmp_153 = tmp_151*tmp_35;
      real_t tmp_154 = tmp_148*tmp_37;
      real_t tmp_155 = tmp_151*tmp_39;
      real_t tmp_156 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_157 = tmp_156*tmp_41;
      real_t tmp_158 = tmp_156*tmp_48;
      real_t tmp_159 = tmp_156*tmp_50;
      real_t tmp_160 = 0.020848748529055869*tmp_55;
      real_t tmp_161 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_162 = tmp_161*tmp_8;
      real_t tmp_163 = tmp_161*tmp_26;
      real_t tmp_164 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_165 = tmp_164*tmp_28;
      real_t tmp_166 = tmp_164*tmp_35;
      real_t tmp_167 = tmp_161*tmp_37;
      real_t tmp_168 = tmp_164*tmp_39;
      real_t tmp_169 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_170 = tmp_169*tmp_41;
      real_t tmp_171 = tmp_169*tmp_48;
      real_t tmp_172 = tmp_169*tmp_50;
      real_t tmp_173 = 0.019202922745021479*tmp_55;
      real_t tmp_174 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_175 = tmp_174*tmp_8;
      real_t tmp_176 = tmp_174*tmp_26;
      real_t tmp_177 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_178 = tmp_177*tmp_28;
      real_t tmp_179 = tmp_177*tmp_35;
      real_t tmp_180 = tmp_174*tmp_37;
      real_t tmp_181 = tmp_177*tmp_39;
      real_t tmp_182 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_183 = tmp_182*tmp_41;
      real_t tmp_184 = tmp_182*tmp_48;
      real_t tmp_185 = tmp_182*tmp_50;
      real_t tmp_186 = 0.020848748529055869*tmp_55;
      real_t tmp_187 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_188 = tmp_187*tmp_8;
      real_t tmp_189 = tmp_187*tmp_26;
      real_t tmp_190 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_191 = tmp_190*tmp_28;
      real_t tmp_192 = tmp_190*tmp_35;
      real_t tmp_193 = tmp_187*tmp_37;
      real_t tmp_194 = tmp_190*tmp_39;
      real_t tmp_195 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_196 = tmp_195*tmp_41;
      real_t tmp_197 = tmp_195*tmp_48;
      real_t tmp_198 = tmp_195*tmp_50;
      real_t tmp_199 = 0.019202922745021479*tmp_55;
      real_t tmp_200 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_201 = tmp_200*tmp_8;
      real_t tmp_202 = tmp_200*tmp_26;
      real_t tmp_203 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_204 = tmp_203*tmp_28;
      real_t tmp_205 = tmp_203*tmp_35;
      real_t tmp_206 = tmp_200*tmp_37;
      real_t tmp_207 = tmp_203*tmp_39;
      real_t tmp_208 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_209 = tmp_208*tmp_41;
      real_t tmp_210 = tmp_208*tmp_48;
      real_t tmp_211 = tmp_208*tmp_50;
      real_t tmp_212 = 0.042507265838595799*tmp_55;
      real_t tmp_213 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_214 = tmp_213*tmp_8;
      real_t tmp_215 = tmp_213*tmp_26;
      real_t tmp_216 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_217 = tmp_216*tmp_28;
      real_t tmp_218 = tmp_216*tmp_35;
      real_t tmp_219 = tmp_213*tmp_37;
      real_t tmp_220 = tmp_216*tmp_39;
      real_t tmp_221 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_222 = tmp_221*tmp_41;
      real_t tmp_223 = tmp_221*tmp_48;
      real_t tmp_224 = tmp_221*tmp_50;
      real_t tmp_225 = 0.020848748529055869*tmp_55;
      real_t tmp_226 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_227 = tmp_226*tmp_8;
      real_t tmp_228 = tmp_226*tmp_26;
      real_t tmp_229 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_230 = tmp_229*tmp_28;
      real_t tmp_231 = tmp_229*tmp_35;
      real_t tmp_232 = tmp_226*tmp_37;
      real_t tmp_233 = tmp_229*tmp_39;
      real_t tmp_234 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_235 = tmp_234*tmp_41;
      real_t tmp_236 = tmp_234*tmp_48;
      real_t tmp_237 = tmp_234*tmp_50;
      real_t tmp_238 = 0.0068572537431980923*tmp_55;
      real_t tmp_239 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_240 = tmp_239*tmp_8;
      real_t tmp_241 = tmp_239*tmp_26;
      real_t tmp_242 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_243 = tmp_242*tmp_28;
      real_t tmp_244 = tmp_242*tmp_35;
      real_t tmp_245 = tmp_239*tmp_37;
      real_t tmp_246 = tmp_242*tmp_39;
      real_t tmp_247 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_248 = tmp_247*tmp_41;
      real_t tmp_249 = tmp_247*tmp_48;
      real_t tmp_250 = tmp_247*tmp_50;
      real_t tmp_251 = 0.037198804536718075*tmp_55;
      real_t tmp_252 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_8;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = 0.042507265838595799*tmp_55;
      real_t tmp_265 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_266 = tmp_265*tmp_8;
      real_t tmp_267 = tmp_26*tmp_265;
      real_t tmp_268 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_269 = tmp_268*tmp_28;
      real_t tmp_270 = tmp_268*tmp_35;
      real_t tmp_271 = tmp_265*tmp_37;
      real_t tmp_272 = tmp_268*tmp_39;
      real_t tmp_273 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_274 = tmp_273*tmp_41;
      real_t tmp_275 = tmp_273*tmp_48;
      real_t tmp_276 = tmp_273*tmp_50;
      real_t tmp_277 = 0.0068572537431980923*tmp_55;
      real_t tmp_278 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_279 = tmp_278*tmp_8;
      real_t tmp_280 = tmp_26*tmp_278;
      real_t tmp_281 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_282 = tmp_28*tmp_281;
      real_t tmp_283 = tmp_281*tmp_35;
      real_t tmp_284 = tmp_278*tmp_37;
      real_t tmp_285 = tmp_281*tmp_39;
      real_t tmp_286 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_287 = tmp_286*tmp_41;
      real_t tmp_288 = tmp_286*tmp_48;
      real_t tmp_289 = tmp_286*tmp_50;
      real_t tmp_290 = 0.037198804536718075*tmp_55;
      real_t tmp_291 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_292 = tmp_291*tmp_8;
      real_t tmp_293 = tmp_26*tmp_291;
      real_t tmp_294 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_295 = tmp_28*tmp_294;
      real_t tmp_296 = tmp_294*tmp_35;
      real_t tmp_297 = tmp_291*tmp_37;
      real_t tmp_298 = tmp_294*tmp_39;
      real_t tmp_299 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_300 = tmp_299*tmp_41;
      real_t tmp_301 = tmp_299*tmp_48;
      real_t tmp_302 = tmp_299*tmp_50;
      real_t tmp_303 = 0.042507265838595799*tmp_55;
      real_t tmp_304 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_305 = tmp_304*tmp_8;
      real_t tmp_306 = tmp_26*tmp_304;
      real_t tmp_307 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_308 = tmp_28*tmp_307;
      real_t tmp_309 = tmp_307*tmp_35;
      real_t tmp_310 = tmp_304*tmp_37;
      real_t tmp_311 = tmp_307*tmp_39;
      real_t tmp_312 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_313 = tmp_312*tmp_41;
      real_t tmp_314 = tmp_312*tmp_48;
      real_t tmp_315 = tmp_312*tmp_50;
      real_t tmp_316 = 0.019202922745021479*tmp_55;
      real_t a_0_0 = tmp_108*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_105 - tmp_106 - tmp_107 - tmp_97 - tmp_98 + 1) + tmp_121*(-tmp_110 - tmp_111 - tmp_113 - tmp_114 - tmp_115 - tmp_116 - tmp_118 - tmp_119 - tmp_120 + 1) + tmp_134*(-tmp_123 - tmp_124 - tmp_126 - tmp_127 - tmp_128 - tmp_129 - tmp_131 - tmp_132 - tmp_133 + 1) + tmp_147*(-tmp_136 - tmp_137 - tmp_139 - tmp_140 - tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 + 1) + tmp_160*(-tmp_149 - tmp_150 - tmp_152 - tmp_153 - tmp_154 - tmp_155 - tmp_157 - tmp_158 - tmp_159 + 1) + tmp_173*(-tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 - tmp_168 - tmp_170 - tmp_171 - tmp_172 + 1) + tmp_186*(-tmp_175 - tmp_176 - tmp_178 - tmp_179 - tmp_180 - tmp_181 - tmp_183 - tmp_184 - tmp_185 + 1) + tmp_199*(-tmp_188 - tmp_189 - tmp_191 - tmp_192 - tmp_193 - tmp_194 - tmp_196 - tmp_197 - tmp_198 + 1) + tmp_212*(-tmp_201 - tmp_202 - tmp_204 - tmp_205 - tmp_206 - tmp_207 - tmp_209 - tmp_210 - tmp_211 + 1) + tmp_225*(-tmp_214 - tmp_215 - tmp_217 - tmp_218 - tmp_219 - tmp_220 - tmp_222 - tmp_223 - tmp_224 + 1) + tmp_238*(-tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 - tmp_233 - tmp_235 - tmp_236 - tmp_237 + 1) + tmp_251*(-tmp_240 - tmp_241 - tmp_243 - tmp_244 - tmp_245 - tmp_246 - tmp_248 - tmp_249 - tmp_250 + 1) + tmp_264*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_277*(-tmp_266 - tmp_267 - tmp_269 - tmp_270 - tmp_271 - tmp_272 - tmp_274 - tmp_275 - tmp_276 + 1) + tmp_290*(-tmp_279 - tmp_280 - tmp_282 - tmp_283 - tmp_284 - tmp_285 - tmp_287 - tmp_288 - tmp_289 + 1) + tmp_303*(-tmp_292 - tmp_293 - tmp_295 - tmp_296 - tmp_297 - tmp_298 - tmp_300 - tmp_301 - tmp_302 + 1) + tmp_316*(-tmp_305 - tmp_306 - tmp_308 - tmp_309 - tmp_310 - tmp_311 - tmp_313 - tmp_314 - tmp_315 + 1) + tmp_56*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_69*(-tmp_58 - tmp_59 - tmp_61 - tmp_62 - tmp_63 - tmp_64 - tmp_66 - tmp_67 - tmp_68 + 1) + tmp_82*(-tmp_71 - tmp_72 - tmp_74 - tmp_75 - tmp_76 - tmp_77 - tmp_79 - tmp_80 - tmp_81 + 1) + tmp_95*(-tmp_84 - tmp_85 - tmp_87 - tmp_88 - tmp_89 - tmp_90 - tmp_92 - tmp_93 - tmp_94 + 1);
      real_t a_1_0 = tmp_108*(tmp_102 + tmp_103 + tmp_107) + tmp_121*(tmp_115 + tmp_116 + tmp_120) + tmp_134*(tmp_128 + tmp_129 + tmp_133) + tmp_147*(tmp_141 + tmp_142 + tmp_146) + tmp_160*(tmp_154 + tmp_155 + tmp_159) + tmp_173*(tmp_167 + tmp_168 + tmp_172) + tmp_186*(tmp_180 + tmp_181 + tmp_185) + tmp_199*(tmp_193 + tmp_194 + tmp_198) + tmp_212*(tmp_206 + tmp_207 + tmp_211) + tmp_225*(tmp_219 + tmp_220 + tmp_224) + tmp_238*(tmp_232 + tmp_233 + tmp_237) + tmp_251*(tmp_245 + tmp_246 + tmp_250) + tmp_264*(tmp_258 + tmp_259 + tmp_263) + tmp_277*(tmp_271 + tmp_272 + tmp_276) + tmp_290*(tmp_284 + tmp_285 + tmp_289) + tmp_303*(tmp_297 + tmp_298 + tmp_302) + tmp_316*(tmp_310 + tmp_311 + tmp_315) + tmp_56*(tmp_38 + tmp_40 + tmp_51) + tmp_69*(tmp_63 + tmp_64 + tmp_68) + tmp_82*(tmp_76 + tmp_77 + tmp_81) + tmp_95*(tmp_89 + tmp_90 + tmp_94);
      real_t a_2_0 = tmp_108*(tmp_101 + tmp_106 + tmp_98) + tmp_121*(tmp_111 + tmp_114 + tmp_119) + tmp_134*(tmp_124 + tmp_127 + tmp_132) + tmp_147*(tmp_137 + tmp_140 + tmp_145) + tmp_160*(tmp_150 + tmp_153 + tmp_158) + tmp_173*(tmp_163 + tmp_166 + tmp_171) + tmp_186*(tmp_176 + tmp_179 + tmp_184) + tmp_199*(tmp_189 + tmp_192 + tmp_197) + tmp_212*(tmp_202 + tmp_205 + tmp_210) + tmp_225*(tmp_215 + tmp_218 + tmp_223) + tmp_238*(tmp_228 + tmp_231 + tmp_236) + tmp_251*(tmp_241 + tmp_244 + tmp_249) + tmp_264*(tmp_254 + tmp_257 + tmp_262) + tmp_277*(tmp_267 + tmp_270 + tmp_275) + tmp_290*(tmp_280 + tmp_283 + tmp_288) + tmp_303*(tmp_293 + tmp_296 + tmp_301) + tmp_316*(tmp_306 + tmp_309 + tmp_314) + tmp_56*(tmp_27 + tmp_36 + tmp_49) + tmp_69*(tmp_59 + tmp_62 + tmp_67) + tmp_82*(tmp_72 + tmp_75 + tmp_80) + tmp_95*(tmp_85 + tmp_88 + tmp_93);
      real_t a_3_0 = tmp_108*(tmp_100 + tmp_105 + tmp_97) + tmp_121*(tmp_110 + tmp_113 + tmp_118) + tmp_134*(tmp_123 + tmp_126 + tmp_131) + tmp_147*(tmp_136 + tmp_139 + tmp_144) + tmp_160*(tmp_149 + tmp_152 + tmp_157) + tmp_173*(tmp_162 + tmp_165 + tmp_170) + tmp_186*(tmp_175 + tmp_178 + tmp_183) + tmp_199*(tmp_188 + tmp_191 + tmp_196) + tmp_212*(tmp_201 + tmp_204 + tmp_209) + tmp_225*(tmp_214 + tmp_217 + tmp_222) + tmp_238*(tmp_227 + tmp_230 + tmp_235) + tmp_251*(tmp_240 + tmp_243 + tmp_248) + tmp_264*(tmp_253 + tmp_256 + tmp_261) + tmp_277*(tmp_266 + tmp_269 + tmp_274) + tmp_290*(tmp_279 + tmp_282 + tmp_287) + tmp_303*(tmp_292 + tmp_295 + tmp_300) + tmp_316*(tmp_305 + tmp_308 + tmp_313) + tmp_56*(tmp_25 + tmp_34 + tmp_47) + tmp_69*(tmp_58 + tmp_61 + tmp_66) + tmp_82*(tmp_71 + tmp_74 + tmp_79) + tmp_95*(tmp_84 + tmp_87 + tmp_92);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }

public:



};




class EGDivtForm_EP0 : public hyteg::dg::DGForm
{

 public:
    EGDivtForm_EP0()

    {}





 protected:
  void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ),
                   trialBasis.numDoFsPerElement( 2, trialDegree ) );

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

      elMat( 0, 0) = 0;
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
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

      elMat( 0, 0) = 0;
   }
   void integrateVolume3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = -p_affine_0_2;
      real_t tmp_5 = p_affine_3_2 + tmp_4;
      real_t tmp_6 = tmp_3*tmp_5;
      real_t tmp_7 = p_affine_3_1 + tmp_2;
      real_t tmp_8 = p_affine_2_2 + tmp_4;
      real_t tmp_9 = tmp_7*tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_4;
      real_t tmp_11 = p_affine_2_0 + tmp_0;
      real_t tmp_12 = tmp_11*tmp_7;
      real_t tmp_13 = p_affine_1_1 + tmp_2;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = tmp_14*tmp_8;
      real_t tmp_16 = tmp_11*tmp_5;
      real_t tmp_17 = tmp_14*tmp_3;
      real_t tmp_18 = 1.0 / (tmp_1*tmp_6 - tmp_1*tmp_9 + tmp_10*tmp_12 - tmp_10*tmp_17 + tmp_13*tmp_15 - tmp_13*tmp_16);
      real_t tmp_19 = p_affine_0_0*p_affine_1_1;
      real_t tmp_20 = p_affine_0_0*p_affine_1_2;
      real_t tmp_21 = p_affine_2_1*p_affine_3_2;
      real_t tmp_22 = p_affine_0_1*p_affine_1_0;
      real_t tmp_23 = p_affine_0_1*p_affine_1_2;
      real_t tmp_24 = p_affine_2_2*p_affine_3_0;
      real_t tmp_25 = p_affine_0_2*p_affine_1_0;
      real_t tmp_26 = p_affine_0_2*p_affine_1_1;
      real_t tmp_27 = p_affine_2_0*p_affine_3_1;
      real_t tmp_28 = p_affine_2_2*p_affine_3_1;
      real_t tmp_29 = p_affine_2_0*p_affine_3_2;
      real_t tmp_30 = p_affine_2_1*p_affine_3_0;
      real_t tmp_31 = (-tmp_1*tmp_18*(tmp_6 - tmp_9) - tmp_10*tmp_18*(tmp_12 - tmp_17) - tmp_11*tmp_18*(tmp_10*tmp_7 - tmp_13*tmp_5) - tmp_13*tmp_18*(tmp_15 - tmp_16) - tmp_14*tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_8) - tmp_18*tmp_3*(tmp_1*tmp_5 - tmp_10*tmp_14) - tmp_18*tmp_5*(tmp_1*tmp_3 - tmp_11*tmp_13) - tmp_18*tmp_7*(-tmp_1*tmp_8 + tmp_10*tmp_11) - tmp_18*tmp_8*(-tmp_1*tmp_7 + tmp_13*tmp_14))*std::abs(p_affine_0_0*tmp_21 - p_affine_0_0*tmp_28 + p_affine_0_1*tmp_24 - p_affine_0_1*tmp_29 + p_affine_0_2*tmp_27 - p_affine_0_2*tmp_30 - p_affine_1_0*tmp_21 + p_affine_1_0*tmp_28 - p_affine_1_1*tmp_24 + p_affine_1_1*tmp_29 - p_affine_1_2*tmp_27 + p_affine_1_2*tmp_30 + p_affine_2_0*tmp_23 - p_affine_2_0*tmp_26 - p_affine_2_1*tmp_20 + p_affine_2_1*tmp_25 + p_affine_2_2*tmp_19 - p_affine_2_2*tmp_22 - p_affine_3_0*tmp_23 + p_affine_3_0*tmp_26 + p_affine_3_1*tmp_20 - p_affine_3_1*tmp_25 - p_affine_3_2*tmp_19 + p_affine_3_2*tmp_22);
      real_t a_0_0 = 0.1666666666666668*tmp_31;
      elMat( 0, 0) = a_0_0;
   }



   void integrateFacetInner3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
                                                     const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                     const Eigen::Matrix< real_t, 3, 1 >&,
                                                     const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                     const DGBasisInfo&                                       trialBasis,
                                                     const DGBasisInfo&                                       testBasis,
                                                     int                                                      trialDegree,
                                                     int                                                      testDegree,
                                                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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

         real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = p_affine_2_0 + tmp_0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_3_1 + tmp_3;
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = p_affine_3_0 + tmp_0;
      real_t tmp_7 = p_affine_2_1 + tmp_3;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = tmp_5 - tmp_8;
      real_t tmp_10 = -p_affine_0_2;
      real_t tmp_11 = p_affine_3_2 + tmp_10;
      real_t tmp_12 = tmp_11*tmp_7;
      real_t tmp_13 = p_affine_1_2 + tmp_10;
      real_t tmp_14 = p_affine_1_1 + tmp_3;
      real_t tmp_15 = p_affine_2_2 + tmp_10;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_15*tmp_4;
      real_t tmp_18 = tmp_11*tmp_2;
      real_t tmp_19 = 1.0 / (tmp_1*tmp_12 - tmp_1*tmp_17 + tmp_13*tmp_5 - tmp_13*tmp_8 + tmp_14*tmp_16 - tmp_14*tmp_18);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_10;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_16 - tmp_18;
      real_t tmp_26 = -p_affine_8_1;
      real_t tmp_27 = p_affine_9_1 + tmp_26;
      real_t tmp_28 = p_affine_10_1 + tmp_26;
      real_t tmp_29 = p_affine_8_1 + tmp_3;
      real_t tmp_30 = tmp_19*(0.031405749086161582*tmp_27 + 0.93718850182767688*tmp_28 + tmp_29);
      real_t tmp_31 = tmp_12 - tmp_17;
      real_t tmp_32 = -p_affine_8_0;
      real_t tmp_33 = p_affine_9_0 + tmp_32;
      real_t tmp_34 = p_affine_10_0 + tmp_32;
      real_t tmp_35 = p_affine_8_0 + tmp_0;
      real_t tmp_36 = tmp_19*(0.031405749086161582*tmp_33 + 0.93718850182767688*tmp_34 + tmp_35);
      real_t tmp_37 = tmp_24*tmp_9 + tmp_25*tmp_30 + tmp_31*tmp_36 - 1.0/4.0;
      real_t tmp_38 = -tmp_1*tmp_4 + tmp_14*tmp_6;
      real_t tmp_39 = tmp_1*tmp_11 - tmp_13*tmp_6;
      real_t tmp_40 = -tmp_11*tmp_14 + tmp_13*tmp_4;
      real_t tmp_41 = tmp_24*tmp_38 + tmp_30*tmp_39 + tmp_36*tmp_40 - 1.0/4.0;
      real_t tmp_42 = tmp_1*tmp_7 - tmp_14*tmp_2;
      real_t tmp_43 = -tmp_1*tmp_15 + tmp_13*tmp_2;
      real_t tmp_44 = -tmp_13*tmp_7 + tmp_14*tmp_15;
      real_t tmp_45 = tmp_24*tmp_42 + tmp_30*tmp_43 + tmp_36*tmp_44 - 1.0/4.0;
      real_t tmp_46 = 0.5*p_affine_13_0;
      real_t tmp_47 = 0.5*p_affine_13_1;
      real_t tmp_48 = 0.5*p_affine_13_2;
      real_t tmp_49 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_50 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_51 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_52 = 1.0*std::pow((std::abs(tmp_22*tmp_49 - tmp_28*tmp_51)*std::abs(tmp_22*tmp_49 - tmp_28*tmp_51)) + (std::abs(tmp_22*tmp_50 - tmp_34*tmp_51)*std::abs(tmp_22*tmp_50 - tmp_34*tmp_51)) + (std::abs(tmp_28*tmp_50 - tmp_34*tmp_49)*std::abs(tmp_28*tmp_50 - tmp_34*tmp_49)), 1.0/2.0);
      real_t tmp_53 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_54 = tmp_19*(0.19601935860219369*tmp_27 + 0.60796128279561268*tmp_28 + tmp_29);
      real_t tmp_55 = tmp_19*(0.19601935860219369*tmp_33 + 0.60796128279561268*tmp_34 + tmp_35);
      real_t tmp_56 = tmp_25*tmp_54 + tmp_31*tmp_55 + tmp_53*tmp_9 - 1.0/4.0;
      real_t tmp_57 = tmp_38*tmp_53 + tmp_39*tmp_54 + tmp_40*tmp_55 - 1.0/4.0;
      real_t tmp_58 = tmp_42*tmp_53 + tmp_43*tmp_54 + tmp_44*tmp_55 - 1.0/4.0;
      real_t tmp_59 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_60 = tmp_19*(0.37605877282253791*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29);
      real_t tmp_61 = tmp_19*(0.37605877282253791*tmp_33 + 0.039308471900058539*tmp_34 + tmp_35);
      real_t tmp_62 = tmp_25*tmp_60 + tmp_31*tmp_61 + tmp_59*tmp_9 - 1.0/4.0;
      real_t tmp_63 = tmp_38*tmp_59 + tmp_39*tmp_60 + tmp_40*tmp_61 - 1.0/4.0;
      real_t tmp_64 = tmp_42*tmp_59 + tmp_43*tmp_60 + tmp_44*tmp_61 - 1.0/4.0;
      real_t tmp_65 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_66 = tmp_19*(0.78764240869137092*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29);
      real_t tmp_67 = tmp_19*(0.78764240869137092*tmp_33 + 0.1711304259088916*tmp_34 + tmp_35);
      real_t tmp_68 = tmp_25*tmp_66 + tmp_31*tmp_67 + tmp_65*tmp_9 - 1.0/4.0;
      real_t tmp_69 = tmp_38*tmp_65 + tmp_39*tmp_66 + tmp_40*tmp_67 - 1.0/4.0;
      real_t tmp_70 = tmp_42*tmp_65 + tmp_43*tmp_66 + tmp_44*tmp_67 - 1.0/4.0;
      real_t tmp_71 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_72 = tmp_19*(0.58463275527740355*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29);
      real_t tmp_73 = tmp_19*(0.58463275527740355*tmp_33 + 0.37605877282253791*tmp_34 + tmp_35);
      real_t tmp_74 = tmp_25*tmp_72 + tmp_31*tmp_73 + tmp_71*tmp_9 - 1.0/4.0;
      real_t tmp_75 = tmp_38*tmp_71 + tmp_39*tmp_72 + tmp_40*tmp_73 - 1.0/4.0;
      real_t tmp_76 = tmp_42*tmp_71 + tmp_43*tmp_72 + tmp_44*tmp_73 - 1.0/4.0;
      real_t tmp_77 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_78 = tmp_19*(0.041227165399737475*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29);
      real_t tmp_79 = tmp_19*(0.041227165399737475*tmp_33 + 0.78764240869137092*tmp_34 + tmp_35);
      real_t tmp_80 = tmp_25*tmp_78 + tmp_31*tmp_79 + tmp_77*tmp_9 - 1.0/4.0;
      real_t tmp_81 = tmp_38*tmp_77 + tmp_39*tmp_78 + tmp_40*tmp_79 - 1.0/4.0;
      real_t tmp_82 = tmp_42*tmp_77 + tmp_43*tmp_78 + tmp_44*tmp_79 - 1.0/4.0;
      real_t tmp_83 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_84 = tmp_19*(0.039308471900058539*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29);
      real_t tmp_85 = tmp_19*(0.039308471900058539*tmp_33 + 0.58463275527740355*tmp_34 + tmp_35);
      real_t tmp_86 = tmp_25*tmp_84 + tmp_31*tmp_85 + tmp_83*tmp_9 - 1.0/4.0;
      real_t tmp_87 = tmp_38*tmp_83 + tmp_39*tmp_84 + tmp_40*tmp_85 - 1.0/4.0;
      real_t tmp_88 = tmp_42*tmp_83 + tmp_43*tmp_84 + tmp_44*tmp_85 - 1.0/4.0;
      real_t tmp_89 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_90 = tmp_19*(0.78764240869137092*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29);
      real_t tmp_91 = tmp_19*(0.78764240869137092*tmp_33 + 0.041227165399737475*tmp_34 + tmp_35);
      real_t tmp_92 = tmp_25*tmp_90 + tmp_31*tmp_91 + tmp_89*tmp_9 - 1.0/4.0;
      real_t tmp_93 = tmp_38*tmp_89 + tmp_39*tmp_90 + tmp_40*tmp_91 - 1.0/4.0;
      real_t tmp_94 = tmp_42*tmp_89 + tmp_43*tmp_90 + tmp_44*tmp_91 - 1.0/4.0;
      real_t tmp_95 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_96 = tmp_19*(0.58463275527740355*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29);
      real_t tmp_97 = tmp_19*(0.58463275527740355*tmp_33 + 0.039308471900058539*tmp_34 + tmp_35);
      real_t tmp_98 = tmp_25*tmp_96 + tmp_31*tmp_97 + tmp_9*tmp_95 - 1.0/4.0;
      real_t tmp_99 = tmp_38*tmp_95 + tmp_39*tmp_96 + tmp_40*tmp_97 - 1.0/4.0;
      real_t tmp_100 = tmp_42*tmp_95 + tmp_43*tmp_96 + tmp_44*tmp_97 - 1.0/4.0;
      real_t tmp_101 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_102 = tmp_19*(0.1711304259088916*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29);
      real_t tmp_103 = tmp_19*(0.1711304259088916*tmp_33 + 0.78764240869137092*tmp_34 + tmp_35);
      real_t tmp_104 = tmp_101*tmp_9 + tmp_102*tmp_25 + tmp_103*tmp_31 - 1.0/4.0;
      real_t tmp_105 = tmp_101*tmp_38 + tmp_102*tmp_39 + tmp_103*tmp_40 - 1.0/4.0;
      real_t tmp_106 = tmp_101*tmp_42 + tmp_102*tmp_43 + tmp_103*tmp_44 - 1.0/4.0;
      real_t tmp_107 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_108 = tmp_19*(0.37605877282253791*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29);
      real_t tmp_109 = tmp_19*(0.37605877282253791*tmp_33 + 0.58463275527740355*tmp_34 + tmp_35);
      real_t tmp_110 = tmp_107*tmp_9 + tmp_108*tmp_25 + tmp_109*tmp_31 - 1.0/4.0;
      real_t tmp_111 = tmp_107*tmp_38 + tmp_108*tmp_39 + tmp_109*tmp_40 - 1.0/4.0;
      real_t tmp_112 = tmp_107*tmp_42 + tmp_108*tmp_43 + tmp_109*tmp_44 - 1.0/4.0;
      real_t tmp_113 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_114 = tmp_19*(0.041227165399737475*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29);
      real_t tmp_115 = tmp_19*(0.041227165399737475*tmp_33 + 0.1711304259088916*tmp_34 + tmp_35);
      real_t tmp_116 = tmp_113*tmp_9 + tmp_114*tmp_25 + tmp_115*tmp_31 - 1.0/4.0;
      real_t tmp_117 = tmp_113*tmp_38 + tmp_114*tmp_39 + tmp_115*tmp_40 - 1.0/4.0;
      real_t tmp_118 = tmp_113*tmp_42 + tmp_114*tmp_43 + tmp_115*tmp_44 - 1.0/4.0;
      real_t tmp_119 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_120 = tmp_19*(0.40446199974765351*tmp_27 + 0.19107600050469298*tmp_28 + tmp_29);
      real_t tmp_121 = tmp_19*(0.40446199974765351*tmp_33 + 0.19107600050469298*tmp_34 + tmp_35);
      real_t tmp_122 = tmp_119*tmp_9 + tmp_120*tmp_25 + tmp_121*tmp_31 - 1.0/4.0;
      real_t tmp_123 = tmp_119*tmp_38 + tmp_120*tmp_39 + tmp_121*tmp_40 - 1.0/4.0;
      real_t tmp_124 = tmp_119*tmp_42 + tmp_120*tmp_43 + tmp_121*tmp_44 - 1.0/4.0;
      real_t tmp_125 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_126 = tmp_19*(0.039308471900058539*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29);
      real_t tmp_127 = tmp_19*(0.039308471900058539*tmp_33 + 0.37605877282253791*tmp_34 + tmp_35);
      real_t tmp_128 = tmp_125*tmp_9 + tmp_126*tmp_25 + tmp_127*tmp_31 - 1.0/4.0;
      real_t tmp_129 = tmp_125*tmp_38 + tmp_126*tmp_39 + tmp_127*tmp_40 - 1.0/4.0;
      real_t tmp_130 = tmp_125*tmp_42 + tmp_126*tmp_43 + tmp_127*tmp_44 - 1.0/4.0;
      real_t tmp_131 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_132 = tmp_19*(0.93718850182767688*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29);
      real_t tmp_133 = tmp_19*(0.93718850182767688*tmp_33 + 0.031405749086161582*tmp_34 + tmp_35);
      real_t tmp_134 = tmp_131*tmp_9 + tmp_132*tmp_25 + tmp_133*tmp_31 - 1.0/4.0;
      real_t tmp_135 = tmp_131*tmp_38 + tmp_132*tmp_39 + tmp_133*tmp_40 - 1.0/4.0;
      real_t tmp_136 = tmp_131*tmp_42 + tmp_132*tmp_43 + tmp_133*tmp_44 - 1.0/4.0;
      real_t tmp_137 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_138 = tmp_19*(0.60796128279561268*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29);
      real_t tmp_139 = tmp_19*(0.60796128279561268*tmp_33 + 0.19601935860219369*tmp_34 + tmp_35);
      real_t tmp_140 = tmp_137*tmp_9 + tmp_138*tmp_25 + tmp_139*tmp_31 - 1.0/4.0;
      real_t tmp_141 = tmp_137*tmp_38 + tmp_138*tmp_39 + tmp_139*tmp_40 - 1.0/4.0;
      real_t tmp_142 = tmp_137*tmp_42 + tmp_138*tmp_43 + tmp_139*tmp_44 - 1.0/4.0;
      real_t tmp_143 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_144 = tmp_19*(0.19107600050469298*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29);
      real_t tmp_145 = tmp_19*(0.19107600050469298*tmp_33 + 0.40446199974765351*tmp_34 + tmp_35);
      real_t tmp_146 = tmp_143*tmp_9 + tmp_144*tmp_25 + tmp_145*tmp_31 - 1.0/4.0;
      real_t tmp_147 = tmp_143*tmp_38 + tmp_144*tmp_39 + tmp_145*tmp_40 - 1.0/4.0;
      real_t tmp_148 = tmp_143*tmp_42 + tmp_144*tmp_43 + tmp_145*tmp_44 - 1.0/4.0;
      real_t tmp_149 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_150 = tmp_19*(0.031405749086161582*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29);
      real_t tmp_151 = tmp_19*(0.031405749086161582*tmp_33 + 0.031405749086161582*tmp_34 + tmp_35);
      real_t tmp_152 = tmp_149*tmp_9 + tmp_150*tmp_25 + tmp_151*tmp_31 - 1.0/4.0;
      real_t tmp_153 = tmp_149*tmp_38 + tmp_150*tmp_39 + tmp_151*tmp_40 - 1.0/4.0;
      real_t tmp_154 = tmp_149*tmp_42 + tmp_150*tmp_43 + tmp_151*tmp_44 - 1.0/4.0;
      real_t tmp_155 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_156 = tmp_19*(0.19601935860219369*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29);
      real_t tmp_157 = tmp_19*(0.19601935860219369*tmp_33 + 0.19601935860219369*tmp_34 + tmp_35);
      real_t tmp_158 = tmp_155*tmp_9 + tmp_156*tmp_25 + tmp_157*tmp_31 - 1.0/4.0;
      real_t tmp_159 = tmp_155*tmp_38 + tmp_156*tmp_39 + tmp_157*tmp_40 - 1.0/4.0;
      real_t tmp_160 = tmp_155*tmp_42 + tmp_156*tmp_43 + tmp_157*tmp_44 - 1.0/4.0;
      real_t tmp_161 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_162 = tmp_19*(0.40446199974765351*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29);
      real_t tmp_163 = tmp_19*(0.40446199974765351*tmp_33 + 0.40446199974765351*tmp_34 + tmp_35);
      real_t tmp_164 = tmp_161*tmp_9 + tmp_162*tmp_25 + tmp_163*tmp_31 - 1.0/4.0;
      real_t tmp_165 = tmp_161*tmp_38 + tmp_162*tmp_39 + tmp_163*tmp_40 - 1.0/4.0;
      real_t tmp_166 = tmp_161*tmp_42 + tmp_162*tmp_43 + tmp_163*tmp_44 - 1.0/4.0;
      real_t tmp_167 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_168 = tmp_19*(0.1711304259088916*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29);
      real_t tmp_169 = tmp_19*(0.1711304259088916*tmp_33 + 0.041227165399737475*tmp_34 + tmp_35);
      real_t tmp_170 = tmp_167*tmp_9 + tmp_168*tmp_25 + tmp_169*tmp_31 - 1.0/4.0;
      real_t tmp_171 = tmp_167*tmp_38 + tmp_168*tmp_39 + tmp_169*tmp_40 - 1.0/4.0;
      real_t tmp_172 = tmp_167*tmp_42 + tmp_168*tmp_43 + tmp_169*tmp_44 - 1.0/4.0;
      real_t a_0_0 = 0.019202922745021479*tmp_52*(tmp_46*(tmp_1*tmp_104 + tmp_105*tmp_2 + tmp_106*tmp_6) + tmp_47*(tmp_104*tmp_14 + tmp_105*tmp_7 + tmp_106*tmp_4) + tmp_48*(tmp_104*tmp_13 + tmp_105*tmp_15 + tmp_106*tmp_11)) + 0.020848748529055869*tmp_52*(tmp_46*(tmp_1*tmp_110 + tmp_111*tmp_2 + tmp_112*tmp_6) + tmp_47*(tmp_110*tmp_14 + tmp_111*tmp_7 + tmp_112*tmp_4) + tmp_48*(tmp_11*tmp_112 + tmp_110*tmp_13 + tmp_111*tmp_15)) + 0.019202922745021479*tmp_52*(tmp_46*(tmp_1*tmp_116 + tmp_117*tmp_2 + tmp_118*tmp_6) + tmp_47*(tmp_116*tmp_14 + tmp_117*tmp_7 + tmp_118*tmp_4) + tmp_48*(tmp_11*tmp_118 + tmp_116*tmp_13 + tmp_117*tmp_15)) + 0.042507265838595799*tmp_52*(tmp_46*(tmp_1*tmp_122 + tmp_123*tmp_2 + tmp_124*tmp_6) + tmp_47*(tmp_122*tmp_14 + tmp_123*tmp_7 + tmp_124*tmp_4) + tmp_48*(tmp_11*tmp_124 + tmp_122*tmp_13 + tmp_123*tmp_15)) + 0.020848748529055869*tmp_52*(tmp_46*(tmp_1*tmp_128 + tmp_129*tmp_2 + tmp_130*tmp_6) + tmp_47*(tmp_128*tmp_14 + tmp_129*tmp_7 + tmp_130*tmp_4) + tmp_48*(tmp_11*tmp_130 + tmp_128*tmp_13 + tmp_129*tmp_15)) + 0.0068572537431980923*tmp_52*(tmp_46*(tmp_1*tmp_134 + tmp_135*tmp_2 + tmp_136*tmp_6) + tmp_47*(tmp_134*tmp_14 + tmp_135*tmp_7 + tmp_136*tmp_4) + tmp_48*(tmp_11*tmp_136 + tmp_13*tmp_134 + tmp_135*tmp_15)) + 0.037198804536718075*tmp_52*(tmp_46*(tmp_1*tmp_140 + tmp_141*tmp_2 + tmp_142*tmp_6) + tmp_47*(tmp_14*tmp_140 + tmp_141*tmp_7 + tmp_142*tmp_4) + tmp_48*(tmp_11*tmp_142 + tmp_13*tmp_140 + tmp_141*tmp_15)) + 0.042507265838595799*tmp_52*(tmp_46*(tmp_1*tmp_146 + tmp_147*tmp_2 + tmp_148*tmp_6) + tmp_47*(tmp_14*tmp_146 + tmp_147*tmp_7 + tmp_148*tmp_4) + tmp_48*(tmp_11*tmp_148 + tmp_13*tmp_146 + tmp_147*tmp_15)) + 0.0068572537431980923*tmp_52*(tmp_46*(tmp_1*tmp_152 + tmp_153*tmp_2 + tmp_154*tmp_6) + tmp_47*(tmp_14*tmp_152 + tmp_153*tmp_7 + tmp_154*tmp_4) + tmp_48*(tmp_11*tmp_154 + tmp_13*tmp_152 + tmp_15*tmp_153)) + 0.037198804536718075*tmp_52*(tmp_46*(tmp_1*tmp_158 + tmp_159*tmp_2 + tmp_160*tmp_6) + tmp_47*(tmp_14*tmp_158 + tmp_159*tmp_7 + tmp_160*tmp_4) + tmp_48*(tmp_11*tmp_160 + tmp_13*tmp_158 + tmp_15*tmp_159)) + 0.042507265838595799*tmp_52*(tmp_46*(tmp_1*tmp_164 + tmp_165*tmp_2 + tmp_166*tmp_6) + tmp_47*(tmp_14*tmp_164 + tmp_165*tmp_7 + tmp_166*tmp_4) + tmp_48*(tmp_11*tmp_166 + tmp_13*tmp_164 + tmp_15*tmp_165)) + 0.019202922745021479*tmp_52*(tmp_46*(tmp_1*tmp_170 + tmp_171*tmp_2 + tmp_172*tmp_6) + tmp_47*(tmp_14*tmp_170 + tmp_171*tmp_7 + tmp_172*tmp_4) + tmp_48*(tmp_11*tmp_172 + tmp_13*tmp_170 + tmp_15*tmp_171)) + 0.0068572537431980923*tmp_52*(tmp_46*(tmp_1*tmp_37 + tmp_2*tmp_41 + tmp_45*tmp_6) + tmp_47*(tmp_14*tmp_37 + tmp_4*tmp_45 + tmp_41*tmp_7) + tmp_48*(tmp_11*tmp_45 + tmp_13*tmp_37 + tmp_15*tmp_41)) + 0.037198804536718075*tmp_52*(tmp_46*(tmp_1*tmp_56 + tmp_2*tmp_57 + tmp_58*tmp_6) + tmp_47*(tmp_14*tmp_56 + tmp_4*tmp_58 + tmp_57*tmp_7) + tmp_48*(tmp_11*tmp_58 + tmp_13*tmp_56 + tmp_15*tmp_57)) + 0.020848748529055869*tmp_52*(tmp_46*(tmp_1*tmp_62 + tmp_2*tmp_63 + tmp_6*tmp_64) + tmp_47*(tmp_14*tmp_62 + tmp_4*tmp_64 + tmp_63*tmp_7) + tmp_48*(tmp_11*tmp_64 + tmp_13*tmp_62 + tmp_15*tmp_63)) + 0.019202922745021479*tmp_52*(tmp_46*(tmp_1*tmp_68 + tmp_2*tmp_69 + tmp_6*tmp_70) + tmp_47*(tmp_14*tmp_68 + tmp_4*tmp_70 + tmp_69*tmp_7) + tmp_48*(tmp_11*tmp_70 + tmp_13*tmp_68 + tmp_15*tmp_69)) + 0.020848748529055869*tmp_52*(tmp_46*(tmp_1*tmp_74 + tmp_2*tmp_75 + tmp_6*tmp_76) + tmp_47*(tmp_14*tmp_74 + tmp_4*tmp_76 + tmp_7*tmp_75) + tmp_48*(tmp_11*tmp_76 + tmp_13*tmp_74 + tmp_15*tmp_75)) + 0.019202922745021479*tmp_52*(tmp_46*(tmp_1*tmp_80 + tmp_2*tmp_81 + tmp_6*tmp_82) + tmp_47*(tmp_14*tmp_80 + tmp_4*tmp_82 + tmp_7*tmp_81) + tmp_48*(tmp_11*tmp_82 + tmp_13*tmp_80 + tmp_15*tmp_81)) + 0.020848748529055869*tmp_52*(tmp_46*(tmp_1*tmp_86 + tmp_2*tmp_87 + tmp_6*tmp_88) + tmp_47*(tmp_14*tmp_86 + tmp_4*tmp_88 + tmp_7*tmp_87) + tmp_48*(tmp_11*tmp_88 + tmp_13*tmp_86 + tmp_15*tmp_87)) + 0.019202922745021479*tmp_52*(tmp_46*(tmp_1*tmp_92 + tmp_2*tmp_93 + tmp_6*tmp_94) + tmp_47*(tmp_14*tmp_92 + tmp_4*tmp_94 + tmp_7*tmp_93) + tmp_48*(tmp_11*tmp_94 + tmp_13*tmp_92 + tmp_15*tmp_93)) + 0.020848748529055869*tmp_52*(tmp_46*(tmp_1*tmp_98 + tmp_100*tmp_6 + tmp_2*tmp_99) + tmp_47*(tmp_100*tmp_4 + tmp_14*tmp_98 + tmp_7*tmp_99) + tmp_48*(tmp_100*tmp_11 + tmp_13*tmp_98 + tmp_15*tmp_99));
      elMat( 0, 0) = a_0_0;
   }




void integrateFacetCoupling3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementInner,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElementOuter,
                                                        const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&,
                                                        const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                        const DGBasisInfo&                                       trialBasis,
                                                        const DGBasisInfo&                                       testBasis,
                                                        int                                                      trialDegree,
                                                        int                                                      testDegree,
                                                        Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = p_affine_2_0 + tmp_0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_3_1 + tmp_3;
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = p_affine_3_0 + tmp_0;
      real_t tmp_7 = p_affine_2_1 + tmp_3;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = tmp_5 - tmp_8;
      real_t tmp_10 = -p_affine_0_2;
      real_t tmp_11 = p_affine_3_2 + tmp_10;
      real_t tmp_12 = tmp_11*tmp_7;
      real_t tmp_13 = p_affine_1_2 + tmp_10;
      real_t tmp_14 = p_affine_1_1 + tmp_3;
      real_t tmp_15 = p_affine_2_2 + tmp_10;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_15*tmp_4;
      real_t tmp_18 = tmp_11*tmp_2;
      real_t tmp_19 = 1.0 / (tmp_1*tmp_12 - tmp_1*tmp_17 + tmp_13*tmp_5 - tmp_13*tmp_8 + tmp_14*tmp_16 - tmp_14*tmp_18);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_10;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_16 - tmp_18;
      real_t tmp_26 = -p_affine_8_1;
      real_t tmp_27 = p_affine_9_1 + tmp_26;
      real_t tmp_28 = p_affine_10_1 + tmp_26;
      real_t tmp_29 = p_affine_8_1 + tmp_3;
      real_t tmp_30 = tmp_19*(0.031405749086161582*tmp_27 + 0.93718850182767688*tmp_28 + tmp_29);
      real_t tmp_31 = tmp_12 - tmp_17;
      real_t tmp_32 = -p_affine_8_0;
      real_t tmp_33 = p_affine_9_0 + tmp_32;
      real_t tmp_34 = p_affine_10_0 + tmp_32;
      real_t tmp_35 = p_affine_8_0 + tmp_0;
      real_t tmp_36 = tmp_19*(0.031405749086161582*tmp_33 + 0.93718850182767688*tmp_34 + tmp_35);
      real_t tmp_37 = tmp_24*tmp_9 + tmp_25*tmp_30 + tmp_31*tmp_36 - 1.0/4.0;
      real_t tmp_38 = -tmp_1*tmp_4 + tmp_14*tmp_6;
      real_t tmp_39 = tmp_1*tmp_11 - tmp_13*tmp_6;
      real_t tmp_40 = -tmp_11*tmp_14 + tmp_13*tmp_4;
      real_t tmp_41 = tmp_24*tmp_38 + tmp_30*tmp_39 + tmp_36*tmp_40 - 1.0/4.0;
      real_t tmp_42 = tmp_1*tmp_7 - tmp_14*tmp_2;
      real_t tmp_43 = -tmp_1*tmp_15 + tmp_13*tmp_2;
      real_t tmp_44 = -tmp_13*tmp_7 + tmp_14*tmp_15;
      real_t tmp_45 = tmp_24*tmp_42 + tmp_30*tmp_43 + tmp_36*tmp_44 - 1.0/4.0;
      real_t tmp_46 = 0.5*p_affine_13_0;
      real_t tmp_47 = 0.5*p_affine_13_1;
      real_t tmp_48 = 0.5*p_affine_13_2;
      real_t tmp_49 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_50 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_51 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_52 = 1.0*std::pow((std::abs(tmp_22*tmp_49 - tmp_28*tmp_51)*std::abs(tmp_22*tmp_49 - tmp_28*tmp_51)) + (std::abs(tmp_22*tmp_50 - tmp_34*tmp_51)*std::abs(tmp_22*tmp_50 - tmp_34*tmp_51)) + (std::abs(tmp_28*tmp_50 - tmp_34*tmp_49)*std::abs(tmp_28*tmp_50 - tmp_34*tmp_49)), 1.0/2.0);
      real_t tmp_53 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_54 = tmp_19*(0.19601935860219369*tmp_27 + 0.60796128279561268*tmp_28 + tmp_29);
      real_t tmp_55 = tmp_19*(0.19601935860219369*tmp_33 + 0.60796128279561268*tmp_34 + tmp_35);
      real_t tmp_56 = tmp_25*tmp_54 + tmp_31*tmp_55 + tmp_53*tmp_9 - 1.0/4.0;
      real_t tmp_57 = tmp_38*tmp_53 + tmp_39*tmp_54 + tmp_40*tmp_55 - 1.0/4.0;
      real_t tmp_58 = tmp_42*tmp_53 + tmp_43*tmp_54 + tmp_44*tmp_55 - 1.0/4.0;
      real_t tmp_59 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_60 = tmp_19*(0.37605877282253791*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29);
      real_t tmp_61 = tmp_19*(0.37605877282253791*tmp_33 + 0.039308471900058539*tmp_34 + tmp_35);
      real_t tmp_62 = tmp_25*tmp_60 + tmp_31*tmp_61 + tmp_59*tmp_9 - 1.0/4.0;
      real_t tmp_63 = tmp_38*tmp_59 + tmp_39*tmp_60 + tmp_40*tmp_61 - 1.0/4.0;
      real_t tmp_64 = tmp_42*tmp_59 + tmp_43*tmp_60 + tmp_44*tmp_61 - 1.0/4.0;
      real_t tmp_65 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_66 = tmp_19*(0.78764240869137092*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29);
      real_t tmp_67 = tmp_19*(0.78764240869137092*tmp_33 + 0.1711304259088916*tmp_34 + tmp_35);
      real_t tmp_68 = tmp_25*tmp_66 + tmp_31*tmp_67 + tmp_65*tmp_9 - 1.0/4.0;
      real_t tmp_69 = tmp_38*tmp_65 + tmp_39*tmp_66 + tmp_40*tmp_67 - 1.0/4.0;
      real_t tmp_70 = tmp_42*tmp_65 + tmp_43*tmp_66 + tmp_44*tmp_67 - 1.0/4.0;
      real_t tmp_71 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_72 = tmp_19*(0.58463275527740355*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29);
      real_t tmp_73 = tmp_19*(0.58463275527740355*tmp_33 + 0.37605877282253791*tmp_34 + tmp_35);
      real_t tmp_74 = tmp_25*tmp_72 + tmp_31*tmp_73 + tmp_71*tmp_9 - 1.0/4.0;
      real_t tmp_75 = tmp_38*tmp_71 + tmp_39*tmp_72 + tmp_40*tmp_73 - 1.0/4.0;
      real_t tmp_76 = tmp_42*tmp_71 + tmp_43*tmp_72 + tmp_44*tmp_73 - 1.0/4.0;
      real_t tmp_77 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_78 = tmp_19*(0.041227165399737475*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29);
      real_t tmp_79 = tmp_19*(0.041227165399737475*tmp_33 + 0.78764240869137092*tmp_34 + tmp_35);
      real_t tmp_80 = tmp_25*tmp_78 + tmp_31*tmp_79 + tmp_77*tmp_9 - 1.0/4.0;
      real_t tmp_81 = tmp_38*tmp_77 + tmp_39*tmp_78 + tmp_40*tmp_79 - 1.0/4.0;
      real_t tmp_82 = tmp_42*tmp_77 + tmp_43*tmp_78 + tmp_44*tmp_79 - 1.0/4.0;
      real_t tmp_83 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_84 = tmp_19*(0.039308471900058539*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29);
      real_t tmp_85 = tmp_19*(0.039308471900058539*tmp_33 + 0.58463275527740355*tmp_34 + tmp_35);
      real_t tmp_86 = tmp_25*tmp_84 + tmp_31*tmp_85 + tmp_83*tmp_9 - 1.0/4.0;
      real_t tmp_87 = tmp_38*tmp_83 + tmp_39*tmp_84 + tmp_40*tmp_85 - 1.0/4.0;
      real_t tmp_88 = tmp_42*tmp_83 + tmp_43*tmp_84 + tmp_44*tmp_85 - 1.0/4.0;
      real_t tmp_89 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_90 = tmp_19*(0.78764240869137092*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29);
      real_t tmp_91 = tmp_19*(0.78764240869137092*tmp_33 + 0.041227165399737475*tmp_34 + tmp_35);
      real_t tmp_92 = tmp_25*tmp_90 + tmp_31*tmp_91 + tmp_89*tmp_9 - 1.0/4.0;
      real_t tmp_93 = tmp_38*tmp_89 + tmp_39*tmp_90 + tmp_40*tmp_91 - 1.0/4.0;
      real_t tmp_94 = tmp_42*tmp_89 + tmp_43*tmp_90 + tmp_44*tmp_91 - 1.0/4.0;
      real_t tmp_95 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_96 = tmp_19*(0.58463275527740355*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29);
      real_t tmp_97 = tmp_19*(0.58463275527740355*tmp_33 + 0.039308471900058539*tmp_34 + tmp_35);
      real_t tmp_98 = tmp_25*tmp_96 + tmp_31*tmp_97 + tmp_9*tmp_95 - 1.0/4.0;
      real_t tmp_99 = tmp_38*tmp_95 + tmp_39*tmp_96 + tmp_40*tmp_97 - 1.0/4.0;
      real_t tmp_100 = tmp_42*tmp_95 + tmp_43*tmp_96 + tmp_44*tmp_97 - 1.0/4.0;
      real_t tmp_101 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_102 = tmp_19*(0.1711304259088916*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29);
      real_t tmp_103 = tmp_19*(0.1711304259088916*tmp_33 + 0.78764240869137092*tmp_34 + tmp_35);
      real_t tmp_104 = tmp_101*tmp_9 + tmp_102*tmp_25 + tmp_103*tmp_31 - 1.0/4.0;
      real_t tmp_105 = tmp_101*tmp_38 + tmp_102*tmp_39 + tmp_103*tmp_40 - 1.0/4.0;
      real_t tmp_106 = tmp_101*tmp_42 + tmp_102*tmp_43 + tmp_103*tmp_44 - 1.0/4.0;
      real_t tmp_107 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_108 = tmp_19*(0.37605877282253791*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29);
      real_t tmp_109 = tmp_19*(0.37605877282253791*tmp_33 + 0.58463275527740355*tmp_34 + tmp_35);
      real_t tmp_110 = tmp_107*tmp_9 + tmp_108*tmp_25 + tmp_109*tmp_31 - 1.0/4.0;
      real_t tmp_111 = tmp_107*tmp_38 + tmp_108*tmp_39 + tmp_109*tmp_40 - 1.0/4.0;
      real_t tmp_112 = tmp_107*tmp_42 + tmp_108*tmp_43 + tmp_109*tmp_44 - 1.0/4.0;
      real_t tmp_113 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_114 = tmp_19*(0.041227165399737475*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29);
      real_t tmp_115 = tmp_19*(0.041227165399737475*tmp_33 + 0.1711304259088916*tmp_34 + tmp_35);
      real_t tmp_116 = tmp_113*tmp_9 + tmp_114*tmp_25 + tmp_115*tmp_31 - 1.0/4.0;
      real_t tmp_117 = tmp_113*tmp_38 + tmp_114*tmp_39 + tmp_115*tmp_40 - 1.0/4.0;
      real_t tmp_118 = tmp_113*tmp_42 + tmp_114*tmp_43 + tmp_115*tmp_44 - 1.0/4.0;
      real_t tmp_119 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_120 = tmp_19*(0.40446199974765351*tmp_27 + 0.19107600050469298*tmp_28 + tmp_29);
      real_t tmp_121 = tmp_19*(0.40446199974765351*tmp_33 + 0.19107600050469298*tmp_34 + tmp_35);
      real_t tmp_122 = tmp_119*tmp_9 + tmp_120*tmp_25 + tmp_121*tmp_31 - 1.0/4.0;
      real_t tmp_123 = tmp_119*tmp_38 + tmp_120*tmp_39 + tmp_121*tmp_40 - 1.0/4.0;
      real_t tmp_124 = tmp_119*tmp_42 + tmp_120*tmp_43 + tmp_121*tmp_44 - 1.0/4.0;
      real_t tmp_125 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_126 = tmp_19*(0.039308471900058539*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29);
      real_t tmp_127 = tmp_19*(0.039308471900058539*tmp_33 + 0.37605877282253791*tmp_34 + tmp_35);
      real_t tmp_128 = tmp_125*tmp_9 + tmp_126*tmp_25 + tmp_127*tmp_31 - 1.0/4.0;
      real_t tmp_129 = tmp_125*tmp_38 + tmp_126*tmp_39 + tmp_127*tmp_40 - 1.0/4.0;
      real_t tmp_130 = tmp_125*tmp_42 + tmp_126*tmp_43 + tmp_127*tmp_44 - 1.0/4.0;
      real_t tmp_131 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_132 = tmp_19*(0.93718850182767688*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29);
      real_t tmp_133 = tmp_19*(0.93718850182767688*tmp_33 + 0.031405749086161582*tmp_34 + tmp_35);
      real_t tmp_134 = tmp_131*tmp_9 + tmp_132*tmp_25 + tmp_133*tmp_31 - 1.0/4.0;
      real_t tmp_135 = tmp_131*tmp_38 + tmp_132*tmp_39 + tmp_133*tmp_40 - 1.0/4.0;
      real_t tmp_136 = tmp_131*tmp_42 + tmp_132*tmp_43 + tmp_133*tmp_44 - 1.0/4.0;
      real_t tmp_137 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_138 = tmp_19*(0.60796128279561268*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29);
      real_t tmp_139 = tmp_19*(0.60796128279561268*tmp_33 + 0.19601935860219369*tmp_34 + tmp_35);
      real_t tmp_140 = tmp_137*tmp_9 + tmp_138*tmp_25 + tmp_139*tmp_31 - 1.0/4.0;
      real_t tmp_141 = tmp_137*tmp_38 + tmp_138*tmp_39 + tmp_139*tmp_40 - 1.0/4.0;
      real_t tmp_142 = tmp_137*tmp_42 + tmp_138*tmp_43 + tmp_139*tmp_44 - 1.0/4.0;
      real_t tmp_143 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_144 = tmp_19*(0.19107600050469298*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29);
      real_t tmp_145 = tmp_19*(0.19107600050469298*tmp_33 + 0.40446199974765351*tmp_34 + tmp_35);
      real_t tmp_146 = tmp_143*tmp_9 + tmp_144*tmp_25 + tmp_145*tmp_31 - 1.0/4.0;
      real_t tmp_147 = tmp_143*tmp_38 + tmp_144*tmp_39 + tmp_145*tmp_40 - 1.0/4.0;
      real_t tmp_148 = tmp_143*tmp_42 + tmp_144*tmp_43 + tmp_145*tmp_44 - 1.0/4.0;
      real_t tmp_149 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_150 = tmp_19*(0.031405749086161582*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29);
      real_t tmp_151 = tmp_19*(0.031405749086161582*tmp_33 + 0.031405749086161582*tmp_34 + tmp_35);
      real_t tmp_152 = tmp_149*tmp_9 + tmp_150*tmp_25 + tmp_151*tmp_31 - 1.0/4.0;
      real_t tmp_153 = tmp_149*tmp_38 + tmp_150*tmp_39 + tmp_151*tmp_40 - 1.0/4.0;
      real_t tmp_154 = tmp_149*tmp_42 + tmp_150*tmp_43 + tmp_151*tmp_44 - 1.0/4.0;
      real_t tmp_155 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_156 = tmp_19*(0.19601935860219369*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29);
      real_t tmp_157 = tmp_19*(0.19601935860219369*tmp_33 + 0.19601935860219369*tmp_34 + tmp_35);
      real_t tmp_158 = tmp_155*tmp_9 + tmp_156*tmp_25 + tmp_157*tmp_31 - 1.0/4.0;
      real_t tmp_159 = tmp_155*tmp_38 + tmp_156*tmp_39 + tmp_157*tmp_40 - 1.0/4.0;
      real_t tmp_160 = tmp_155*tmp_42 + tmp_156*tmp_43 + tmp_157*tmp_44 - 1.0/4.0;
      real_t tmp_161 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_162 = tmp_19*(0.40446199974765351*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29);
      real_t tmp_163 = tmp_19*(0.40446199974765351*tmp_33 + 0.40446199974765351*tmp_34 + tmp_35);
      real_t tmp_164 = tmp_161*tmp_9 + tmp_162*tmp_25 + tmp_163*tmp_31 - 1.0/4.0;
      real_t tmp_165 = tmp_161*tmp_38 + tmp_162*tmp_39 + tmp_163*tmp_40 - 1.0/4.0;
      real_t tmp_166 = tmp_161*tmp_42 + tmp_162*tmp_43 + tmp_163*tmp_44 - 1.0/4.0;
      real_t tmp_167 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_168 = tmp_19*(0.1711304259088916*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29);
      real_t tmp_169 = tmp_19*(0.1711304259088916*tmp_33 + 0.041227165399737475*tmp_34 + tmp_35);
      real_t tmp_170 = tmp_167*tmp_9 + tmp_168*tmp_25 + tmp_169*tmp_31 - 1.0/4.0;
      real_t tmp_171 = tmp_167*tmp_38 + tmp_168*tmp_39 + tmp_169*tmp_40 - 1.0/4.0;
      real_t tmp_172 = tmp_167*tmp_42 + tmp_168*tmp_43 + tmp_169*tmp_44 - 1.0/4.0;
      real_t a_0_0 = 0.019202922745021479*tmp_52*(tmp_46*(tmp_1*tmp_104 + tmp_105*tmp_2 + tmp_106*tmp_6) + tmp_47*(tmp_104*tmp_14 + tmp_105*tmp_7 + tmp_106*tmp_4) + tmp_48*(tmp_104*tmp_13 + tmp_105*tmp_15 + tmp_106*tmp_11)) + 0.020848748529055869*tmp_52*(tmp_46*(tmp_1*tmp_110 + tmp_111*tmp_2 + tmp_112*tmp_6) + tmp_47*(tmp_110*tmp_14 + tmp_111*tmp_7 + tmp_112*tmp_4) + tmp_48*(tmp_11*tmp_112 + tmp_110*tmp_13 + tmp_111*tmp_15)) + 0.019202922745021479*tmp_52*(tmp_46*(tmp_1*tmp_116 + tmp_117*tmp_2 + tmp_118*tmp_6) + tmp_47*(tmp_116*tmp_14 + tmp_117*tmp_7 + tmp_118*tmp_4) + tmp_48*(tmp_11*tmp_118 + tmp_116*tmp_13 + tmp_117*tmp_15)) + 0.042507265838595799*tmp_52*(tmp_46*(tmp_1*tmp_122 + tmp_123*tmp_2 + tmp_124*tmp_6) + tmp_47*(tmp_122*tmp_14 + tmp_123*tmp_7 + tmp_124*tmp_4) + tmp_48*(tmp_11*tmp_124 + tmp_122*tmp_13 + tmp_123*tmp_15)) + 0.020848748529055869*tmp_52*(tmp_46*(tmp_1*tmp_128 + tmp_129*tmp_2 + tmp_130*tmp_6) + tmp_47*(tmp_128*tmp_14 + tmp_129*tmp_7 + tmp_130*tmp_4) + tmp_48*(tmp_11*tmp_130 + tmp_128*tmp_13 + tmp_129*tmp_15)) + 0.0068572537431980923*tmp_52*(tmp_46*(tmp_1*tmp_134 + tmp_135*tmp_2 + tmp_136*tmp_6) + tmp_47*(tmp_134*tmp_14 + tmp_135*tmp_7 + tmp_136*tmp_4) + tmp_48*(tmp_11*tmp_136 + tmp_13*tmp_134 + tmp_135*tmp_15)) + 0.037198804536718075*tmp_52*(tmp_46*(tmp_1*tmp_140 + tmp_141*tmp_2 + tmp_142*tmp_6) + tmp_47*(tmp_14*tmp_140 + tmp_141*tmp_7 + tmp_142*tmp_4) + tmp_48*(tmp_11*tmp_142 + tmp_13*tmp_140 + tmp_141*tmp_15)) + 0.042507265838595799*tmp_52*(tmp_46*(tmp_1*tmp_146 + tmp_147*tmp_2 + tmp_148*tmp_6) + tmp_47*(tmp_14*tmp_146 + tmp_147*tmp_7 + tmp_148*tmp_4) + tmp_48*(tmp_11*tmp_148 + tmp_13*tmp_146 + tmp_147*tmp_15)) + 0.0068572537431980923*tmp_52*(tmp_46*(tmp_1*tmp_152 + tmp_153*tmp_2 + tmp_154*tmp_6) + tmp_47*(tmp_14*tmp_152 + tmp_153*tmp_7 + tmp_154*tmp_4) + tmp_48*(tmp_11*tmp_154 + tmp_13*tmp_152 + tmp_15*tmp_153)) + 0.037198804536718075*tmp_52*(tmp_46*(tmp_1*tmp_158 + tmp_159*tmp_2 + tmp_160*tmp_6) + tmp_47*(tmp_14*tmp_158 + tmp_159*tmp_7 + tmp_160*tmp_4) + tmp_48*(tmp_11*tmp_160 + tmp_13*tmp_158 + tmp_15*tmp_159)) + 0.042507265838595799*tmp_52*(tmp_46*(tmp_1*tmp_164 + tmp_165*tmp_2 + tmp_166*tmp_6) + tmp_47*(tmp_14*tmp_164 + tmp_165*tmp_7 + tmp_166*tmp_4) + tmp_48*(tmp_11*tmp_166 + tmp_13*tmp_164 + tmp_15*tmp_165)) + 0.019202922745021479*tmp_52*(tmp_46*(tmp_1*tmp_170 + tmp_171*tmp_2 + tmp_172*tmp_6) + tmp_47*(tmp_14*tmp_170 + tmp_171*tmp_7 + tmp_172*tmp_4) + tmp_48*(tmp_11*tmp_172 + tmp_13*tmp_170 + tmp_15*tmp_171)) + 0.0068572537431980923*tmp_52*(tmp_46*(tmp_1*tmp_37 + tmp_2*tmp_41 + tmp_45*tmp_6) + tmp_47*(tmp_14*tmp_37 + tmp_4*tmp_45 + tmp_41*tmp_7) + tmp_48*(tmp_11*tmp_45 + tmp_13*tmp_37 + tmp_15*tmp_41)) + 0.037198804536718075*tmp_52*(tmp_46*(tmp_1*tmp_56 + tmp_2*tmp_57 + tmp_58*tmp_6) + tmp_47*(tmp_14*tmp_56 + tmp_4*tmp_58 + tmp_57*tmp_7) + tmp_48*(tmp_11*tmp_58 + tmp_13*tmp_56 + tmp_15*tmp_57)) + 0.020848748529055869*tmp_52*(tmp_46*(tmp_1*tmp_62 + tmp_2*tmp_63 + tmp_6*tmp_64) + tmp_47*(tmp_14*tmp_62 + tmp_4*tmp_64 + tmp_63*tmp_7) + tmp_48*(tmp_11*tmp_64 + tmp_13*tmp_62 + tmp_15*tmp_63)) + 0.019202922745021479*tmp_52*(tmp_46*(tmp_1*tmp_68 + tmp_2*tmp_69 + tmp_6*tmp_70) + tmp_47*(tmp_14*tmp_68 + tmp_4*tmp_70 + tmp_69*tmp_7) + tmp_48*(tmp_11*tmp_70 + tmp_13*tmp_68 + tmp_15*tmp_69)) + 0.020848748529055869*tmp_52*(tmp_46*(tmp_1*tmp_74 + tmp_2*tmp_75 + tmp_6*tmp_76) + tmp_47*(tmp_14*tmp_74 + tmp_4*tmp_76 + tmp_7*tmp_75) + tmp_48*(tmp_11*tmp_76 + tmp_13*tmp_74 + tmp_15*tmp_75)) + 0.019202922745021479*tmp_52*(tmp_46*(tmp_1*tmp_80 + tmp_2*tmp_81 + tmp_6*tmp_82) + tmp_47*(tmp_14*tmp_80 + tmp_4*tmp_82 + tmp_7*tmp_81) + tmp_48*(tmp_11*tmp_82 + tmp_13*tmp_80 + tmp_15*tmp_81)) + 0.020848748529055869*tmp_52*(tmp_46*(tmp_1*tmp_86 + tmp_2*tmp_87 + tmp_6*tmp_88) + tmp_47*(tmp_14*tmp_86 + tmp_4*tmp_88 + tmp_7*tmp_87) + tmp_48*(tmp_11*tmp_88 + tmp_13*tmp_86 + tmp_15*tmp_87)) + 0.019202922745021479*tmp_52*(tmp_46*(tmp_1*tmp_92 + tmp_2*tmp_93 + tmp_6*tmp_94) + tmp_47*(tmp_14*tmp_92 + tmp_4*tmp_94 + tmp_7*tmp_93) + tmp_48*(tmp_11*tmp_94 + tmp_13*tmp_92 + tmp_15*tmp_93)) + 0.020848748529055869*tmp_52*(tmp_46*(tmp_1*tmp_98 + tmp_100*tmp_6 + tmp_2*tmp_99) + tmp_47*(tmp_100*tmp_4 + tmp_14*tmp_98 + tmp_7*tmp_99) + tmp_48*(tmp_100*tmp_11 + tmp_13*tmp_98 + tmp_15*tmp_99));
      elMat( 0, 0) = a_0_0;
}



void integrateFacetDirichletBoundary3D(
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsElement,
    const std::vector< Eigen::Matrix< real_t, 3, 1 > >& coordsFacet,
    const Eigen::Matrix< real_t, 3, 1 >&,
    const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
    const DGBasisInfo&                                       trialBasis,
    const DGBasisInfo&                                       testBasis,
    int                                                      trialDegree,
    int                                                      testDegree,
    Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
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


      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = p_affine_2_0 + tmp_0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_3_1 + tmp_3;
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = p_affine_3_0 + tmp_0;
      real_t tmp_7 = p_affine_2_1 + tmp_3;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = tmp_5 - tmp_8;
      real_t tmp_10 = -p_affine_0_2;
      real_t tmp_11 = p_affine_3_2 + tmp_10;
      real_t tmp_12 = tmp_11*tmp_7;
      real_t tmp_13 = p_affine_1_2 + tmp_10;
      real_t tmp_14 = p_affine_1_1 + tmp_3;
      real_t tmp_15 = p_affine_2_2 + tmp_10;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_15*tmp_4;
      real_t tmp_18 = tmp_11*tmp_2;
      real_t tmp_19 = 1.0 / (tmp_1*tmp_12 - tmp_1*tmp_17 + tmp_13*tmp_5 - tmp_13*tmp_8 + tmp_14*tmp_16 - tmp_14*tmp_18);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_10;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_16 - tmp_18;
      real_t tmp_26 = -p_affine_8_1;
      real_t tmp_27 = p_affine_9_1 + tmp_26;
      real_t tmp_28 = p_affine_10_1 + tmp_26;
      real_t tmp_29 = p_affine_8_1 + tmp_3;
      real_t tmp_30 = tmp_19*(0.031405749086161582*tmp_27 + 0.93718850182767688*tmp_28 + tmp_29);
      real_t tmp_31 = tmp_12 - tmp_17;
      real_t tmp_32 = -p_affine_8_0;
      real_t tmp_33 = p_affine_9_0 + tmp_32;
      real_t tmp_34 = p_affine_10_0 + tmp_32;
      real_t tmp_35 = p_affine_8_0 + tmp_0;
      real_t tmp_36 = tmp_19*(0.031405749086161582*tmp_33 + 0.93718850182767688*tmp_34 + tmp_35);
      real_t tmp_37 = tmp_24*tmp_9 + tmp_25*tmp_30 + tmp_31*tmp_36 - 1.0/4.0;
      real_t tmp_38 = -tmp_1*tmp_4 + tmp_14*tmp_6;
      real_t tmp_39 = tmp_1*tmp_11 - tmp_13*tmp_6;
      real_t tmp_40 = -tmp_11*tmp_14 + tmp_13*tmp_4;
      real_t tmp_41 = tmp_24*tmp_38 + tmp_30*tmp_39 + tmp_36*tmp_40 - 1.0/4.0;
      real_t tmp_42 = tmp_1*tmp_7 - tmp_14*tmp_2;
      real_t tmp_43 = -tmp_1*tmp_15 + tmp_13*tmp_2;
      real_t tmp_44 = -tmp_13*tmp_7 + tmp_14*tmp_15;
      real_t tmp_45 = tmp_24*tmp_42 + tmp_30*tmp_43 + tmp_36*tmp_44 - 1.0/4.0;
      real_t tmp_46 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_47 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_48 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_49 = 1.0*std::pow((std::abs(tmp_22*tmp_46 - tmp_28*tmp_48)*std::abs(tmp_22*tmp_46 - tmp_28*tmp_48)) + (std::abs(tmp_22*tmp_47 - tmp_34*tmp_48)*std::abs(tmp_22*tmp_47 - tmp_34*tmp_48)) + (std::abs(tmp_28*tmp_47 - tmp_34*tmp_46)*std::abs(tmp_28*tmp_47 - tmp_34*tmp_46)), 1.0/2.0);
      real_t tmp_50 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_51 = tmp_19*(0.19601935860219369*tmp_27 + 0.60796128279561268*tmp_28 + tmp_29);
      real_t tmp_52 = tmp_19*(0.19601935860219369*tmp_33 + 0.60796128279561268*tmp_34 + tmp_35);
      real_t tmp_53 = tmp_25*tmp_51 + tmp_31*tmp_52 + tmp_50*tmp_9 - 1.0/4.0;
      real_t tmp_54 = tmp_38*tmp_50 + tmp_39*tmp_51 + tmp_40*tmp_52 - 1.0/4.0;
      real_t tmp_55 = tmp_42*tmp_50 + tmp_43*tmp_51 + tmp_44*tmp_52 - 1.0/4.0;
      real_t tmp_56 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_57 = tmp_19*(0.37605877282253791*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29);
      real_t tmp_58 = tmp_19*(0.37605877282253791*tmp_33 + 0.039308471900058539*tmp_34 + tmp_35);
      real_t tmp_59 = tmp_25*tmp_57 + tmp_31*tmp_58 + tmp_56*tmp_9 - 1.0/4.0;
      real_t tmp_60 = tmp_38*tmp_56 + tmp_39*tmp_57 + tmp_40*tmp_58 - 1.0/4.0;
      real_t tmp_61 = tmp_42*tmp_56 + tmp_43*tmp_57 + tmp_44*tmp_58 - 1.0/4.0;
      real_t tmp_62 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_63 = tmp_19*(0.78764240869137092*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29);
      real_t tmp_64 = tmp_19*(0.78764240869137092*tmp_33 + 0.1711304259088916*tmp_34 + tmp_35);
      real_t tmp_65 = tmp_25*tmp_63 + tmp_31*tmp_64 + tmp_62*tmp_9 - 1.0/4.0;
      real_t tmp_66 = tmp_38*tmp_62 + tmp_39*tmp_63 + tmp_40*tmp_64 - 1.0/4.0;
      real_t tmp_67 = tmp_42*tmp_62 + tmp_43*tmp_63 + tmp_44*tmp_64 - 1.0/4.0;
      real_t tmp_68 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_69 = tmp_19*(0.58463275527740355*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29);
      real_t tmp_70 = tmp_19*(0.58463275527740355*tmp_33 + 0.37605877282253791*tmp_34 + tmp_35);
      real_t tmp_71 = tmp_25*tmp_69 + tmp_31*tmp_70 + tmp_68*tmp_9 - 1.0/4.0;
      real_t tmp_72 = tmp_38*tmp_68 + tmp_39*tmp_69 + tmp_40*tmp_70 - 1.0/4.0;
      real_t tmp_73 = tmp_42*tmp_68 + tmp_43*tmp_69 + tmp_44*tmp_70 - 1.0/4.0;
      real_t tmp_74 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_75 = tmp_19*(0.041227165399737475*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29);
      real_t tmp_76 = tmp_19*(0.041227165399737475*tmp_33 + 0.78764240869137092*tmp_34 + tmp_35);
      real_t tmp_77 = tmp_25*tmp_75 + tmp_31*tmp_76 + tmp_74*tmp_9 - 1.0/4.0;
      real_t tmp_78 = tmp_38*tmp_74 + tmp_39*tmp_75 + tmp_40*tmp_76 - 1.0/4.0;
      real_t tmp_79 = tmp_42*tmp_74 + tmp_43*tmp_75 + tmp_44*tmp_76 - 1.0/4.0;
      real_t tmp_80 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_81 = tmp_19*(0.039308471900058539*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29);
      real_t tmp_82 = tmp_19*(0.039308471900058539*tmp_33 + 0.58463275527740355*tmp_34 + tmp_35);
      real_t tmp_83 = tmp_25*tmp_81 + tmp_31*tmp_82 + tmp_80*tmp_9 - 1.0/4.0;
      real_t tmp_84 = tmp_38*tmp_80 + tmp_39*tmp_81 + tmp_40*tmp_82 - 1.0/4.0;
      real_t tmp_85 = tmp_42*tmp_80 + tmp_43*tmp_81 + tmp_44*tmp_82 - 1.0/4.0;
      real_t tmp_86 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_87 = tmp_19*(0.78764240869137092*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29);
      real_t tmp_88 = tmp_19*(0.78764240869137092*tmp_33 + 0.041227165399737475*tmp_34 + tmp_35);
      real_t tmp_89 = tmp_25*tmp_87 + tmp_31*tmp_88 + tmp_86*tmp_9 - 1.0/4.0;
      real_t tmp_90 = tmp_38*tmp_86 + tmp_39*tmp_87 + tmp_40*tmp_88 - 1.0/4.0;
      real_t tmp_91 = tmp_42*tmp_86 + tmp_43*tmp_87 + tmp_44*tmp_88 - 1.0/4.0;
      real_t tmp_92 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_93 = tmp_19*(0.58463275527740355*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29);
      real_t tmp_94 = tmp_19*(0.58463275527740355*tmp_33 + 0.039308471900058539*tmp_34 + tmp_35);
      real_t tmp_95 = tmp_25*tmp_93 + tmp_31*tmp_94 + tmp_9*tmp_92 - 1.0/4.0;
      real_t tmp_96 = tmp_38*tmp_92 + tmp_39*tmp_93 + tmp_40*tmp_94 - 1.0/4.0;
      real_t tmp_97 = tmp_42*tmp_92 + tmp_43*tmp_93 + tmp_44*tmp_94 - 1.0/4.0;
      real_t tmp_98 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_99 = tmp_19*(0.1711304259088916*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29);
      real_t tmp_100 = tmp_19*(0.1711304259088916*tmp_33 + 0.78764240869137092*tmp_34 + tmp_35);
      real_t tmp_101 = tmp_100*tmp_31 + tmp_25*tmp_99 + tmp_9*tmp_98 - 1.0/4.0;
      real_t tmp_102 = tmp_100*tmp_40 + tmp_38*tmp_98 + tmp_39*tmp_99 - 1.0/4.0;
      real_t tmp_103 = tmp_100*tmp_44 + tmp_42*tmp_98 + tmp_43*tmp_99 - 1.0/4.0;
      real_t tmp_104 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_105 = tmp_19*(0.37605877282253791*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29);
      real_t tmp_106 = tmp_19*(0.37605877282253791*tmp_33 + 0.58463275527740355*tmp_34 + tmp_35);
      real_t tmp_107 = tmp_104*tmp_9 + tmp_105*tmp_25 + tmp_106*tmp_31 - 1.0/4.0;
      real_t tmp_108 = tmp_104*tmp_38 + tmp_105*tmp_39 + tmp_106*tmp_40 - 1.0/4.0;
      real_t tmp_109 = tmp_104*tmp_42 + tmp_105*tmp_43 + tmp_106*tmp_44 - 1.0/4.0;
      real_t tmp_110 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_111 = tmp_19*(0.041227165399737475*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29);
      real_t tmp_112 = tmp_19*(0.041227165399737475*tmp_33 + 0.1711304259088916*tmp_34 + tmp_35);
      real_t tmp_113 = tmp_110*tmp_9 + tmp_111*tmp_25 + tmp_112*tmp_31 - 1.0/4.0;
      real_t tmp_114 = tmp_110*tmp_38 + tmp_111*tmp_39 + tmp_112*tmp_40 - 1.0/4.0;
      real_t tmp_115 = tmp_110*tmp_42 + tmp_111*tmp_43 + tmp_112*tmp_44 - 1.0/4.0;
      real_t tmp_116 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_117 = tmp_19*(0.40446199974765351*tmp_27 + 0.19107600050469298*tmp_28 + tmp_29);
      real_t tmp_118 = tmp_19*(0.40446199974765351*tmp_33 + 0.19107600050469298*tmp_34 + tmp_35);
      real_t tmp_119 = tmp_116*tmp_9 + tmp_117*tmp_25 + tmp_118*tmp_31 - 1.0/4.0;
      real_t tmp_120 = tmp_116*tmp_38 + tmp_117*tmp_39 + tmp_118*tmp_40 - 1.0/4.0;
      real_t tmp_121 = tmp_116*tmp_42 + tmp_117*tmp_43 + tmp_118*tmp_44 - 1.0/4.0;
      real_t tmp_122 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_123 = tmp_19*(0.039308471900058539*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29);
      real_t tmp_124 = tmp_19*(0.039308471900058539*tmp_33 + 0.37605877282253791*tmp_34 + tmp_35);
      real_t tmp_125 = tmp_122*tmp_9 + tmp_123*tmp_25 + tmp_124*tmp_31 - 1.0/4.0;
      real_t tmp_126 = tmp_122*tmp_38 + tmp_123*tmp_39 + tmp_124*tmp_40 - 1.0/4.0;
      real_t tmp_127 = tmp_122*tmp_42 + tmp_123*tmp_43 + tmp_124*tmp_44 - 1.0/4.0;
      real_t tmp_128 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_129 = tmp_19*(0.93718850182767688*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29);
      real_t tmp_130 = tmp_19*(0.93718850182767688*tmp_33 + 0.031405749086161582*tmp_34 + tmp_35);
      real_t tmp_131 = tmp_128*tmp_9 + tmp_129*tmp_25 + tmp_130*tmp_31 - 1.0/4.0;
      real_t tmp_132 = tmp_128*tmp_38 + tmp_129*tmp_39 + tmp_130*tmp_40 - 1.0/4.0;
      real_t tmp_133 = tmp_128*tmp_42 + tmp_129*tmp_43 + tmp_130*tmp_44 - 1.0/4.0;
      real_t tmp_134 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_135 = tmp_19*(0.60796128279561268*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29);
      real_t tmp_136 = tmp_19*(0.60796128279561268*tmp_33 + 0.19601935860219369*tmp_34 + tmp_35);
      real_t tmp_137 = tmp_134*tmp_9 + tmp_135*tmp_25 + tmp_136*tmp_31 - 1.0/4.0;
      real_t tmp_138 = tmp_134*tmp_38 + tmp_135*tmp_39 + tmp_136*tmp_40 - 1.0/4.0;
      real_t tmp_139 = tmp_134*tmp_42 + tmp_135*tmp_43 + tmp_136*tmp_44 - 1.0/4.0;
      real_t tmp_140 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_141 = tmp_19*(0.19107600050469298*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29);
      real_t tmp_142 = tmp_19*(0.19107600050469298*tmp_33 + 0.40446199974765351*tmp_34 + tmp_35);
      real_t tmp_143 = tmp_140*tmp_9 + tmp_141*tmp_25 + tmp_142*tmp_31 - 1.0/4.0;
      real_t tmp_144 = tmp_140*tmp_38 + tmp_141*tmp_39 + tmp_142*tmp_40 - 1.0/4.0;
      real_t tmp_145 = tmp_140*tmp_42 + tmp_141*tmp_43 + tmp_142*tmp_44 - 1.0/4.0;
      real_t tmp_146 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_147 = tmp_19*(0.031405749086161582*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29);
      real_t tmp_148 = tmp_19*(0.031405749086161582*tmp_33 + 0.031405749086161582*tmp_34 + tmp_35);
      real_t tmp_149 = tmp_146*tmp_9 + tmp_147*tmp_25 + tmp_148*tmp_31 - 1.0/4.0;
      real_t tmp_150 = tmp_146*tmp_38 + tmp_147*tmp_39 + tmp_148*tmp_40 - 1.0/4.0;
      real_t tmp_151 = tmp_146*tmp_42 + tmp_147*tmp_43 + tmp_148*tmp_44 - 1.0/4.0;
      real_t tmp_152 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_153 = tmp_19*(0.19601935860219369*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29);
      real_t tmp_154 = tmp_19*(0.19601935860219369*tmp_33 + 0.19601935860219369*tmp_34 + tmp_35);
      real_t tmp_155 = tmp_152*tmp_9 + tmp_153*tmp_25 + tmp_154*tmp_31 - 1.0/4.0;
      real_t tmp_156 = tmp_152*tmp_38 + tmp_153*tmp_39 + tmp_154*tmp_40 - 1.0/4.0;
      real_t tmp_157 = tmp_152*tmp_42 + tmp_153*tmp_43 + tmp_154*tmp_44 - 1.0/4.0;
      real_t tmp_158 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_159 = tmp_19*(0.40446199974765351*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29);
      real_t tmp_160 = tmp_19*(0.40446199974765351*tmp_33 + 0.40446199974765351*tmp_34 + tmp_35);
      real_t tmp_161 = tmp_158*tmp_9 + tmp_159*tmp_25 + tmp_160*tmp_31 - 1.0/4.0;
      real_t tmp_162 = tmp_158*tmp_38 + tmp_159*tmp_39 + tmp_160*tmp_40 - 1.0/4.0;
      real_t tmp_163 = tmp_158*tmp_42 + tmp_159*tmp_43 + tmp_160*tmp_44 - 1.0/4.0;
      real_t tmp_164 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_165 = tmp_19*(0.1711304259088916*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29);
      real_t tmp_166 = tmp_19*(0.1711304259088916*tmp_33 + 0.041227165399737475*tmp_34 + tmp_35);
      real_t tmp_167 = tmp_164*tmp_9 + tmp_165*tmp_25 + tmp_166*tmp_31 - 1.0/4.0;
      real_t tmp_168 = tmp_164*tmp_38 + tmp_165*tmp_39 + tmp_166*tmp_40 - 1.0/4.0;
      real_t tmp_169 = tmp_164*tmp_42 + tmp_165*tmp_43 + tmp_166*tmp_44 - 1.0/4.0;
      real_t a_0_0 = 0.019202922745021479*tmp_49*(p_affine_13_0*(tmp_1*tmp_101 + tmp_102*tmp_2 + tmp_103*tmp_6) + p_affine_13_1*(tmp_101*tmp_14 + tmp_102*tmp_7 + tmp_103*tmp_4) + p_affine_13_2*(tmp_101*tmp_13 + tmp_102*tmp_15 + tmp_103*tmp_11)) + 0.020848748529055869*tmp_49*(p_affine_13_0*(tmp_1*tmp_107 + tmp_108*tmp_2 + tmp_109*tmp_6) + p_affine_13_1*(tmp_107*tmp_14 + tmp_108*tmp_7 + tmp_109*tmp_4) + p_affine_13_2*(tmp_107*tmp_13 + tmp_108*tmp_15 + tmp_109*tmp_11)) + 0.019202922745021479*tmp_49*(p_affine_13_0*(tmp_1*tmp_113 + tmp_114*tmp_2 + tmp_115*tmp_6) + p_affine_13_1*(tmp_113*tmp_14 + tmp_114*tmp_7 + tmp_115*tmp_4) + p_affine_13_2*(tmp_11*tmp_115 + tmp_113*tmp_13 + tmp_114*tmp_15)) + 0.042507265838595799*tmp_49*(p_affine_13_0*(tmp_1*tmp_119 + tmp_120*tmp_2 + tmp_121*tmp_6) + p_affine_13_1*(tmp_119*tmp_14 + tmp_120*tmp_7 + tmp_121*tmp_4) + p_affine_13_2*(tmp_11*tmp_121 + tmp_119*tmp_13 + tmp_120*tmp_15)) + 0.020848748529055869*tmp_49*(p_affine_13_0*(tmp_1*tmp_125 + tmp_126*tmp_2 + tmp_127*tmp_6) + p_affine_13_1*(tmp_125*tmp_14 + tmp_126*tmp_7 + tmp_127*tmp_4) + p_affine_13_2*(tmp_11*tmp_127 + tmp_125*tmp_13 + tmp_126*tmp_15)) + 0.0068572537431980923*tmp_49*(p_affine_13_0*(tmp_1*tmp_131 + tmp_132*tmp_2 + tmp_133*tmp_6) + p_affine_13_1*(tmp_131*tmp_14 + tmp_132*tmp_7 + tmp_133*tmp_4) + p_affine_13_2*(tmp_11*tmp_133 + tmp_13*tmp_131 + tmp_132*tmp_15)) + 0.037198804536718075*tmp_49*(p_affine_13_0*(tmp_1*tmp_137 + tmp_138*tmp_2 + tmp_139*tmp_6) + p_affine_13_1*(tmp_137*tmp_14 + tmp_138*tmp_7 + tmp_139*tmp_4) + p_affine_13_2*(tmp_11*tmp_139 + tmp_13*tmp_137 + tmp_138*tmp_15)) + 0.042507265838595799*tmp_49*(p_affine_13_0*(tmp_1*tmp_143 + tmp_144*tmp_2 + tmp_145*tmp_6) + p_affine_13_1*(tmp_14*tmp_143 + tmp_144*tmp_7 + tmp_145*tmp_4) + p_affine_13_2*(tmp_11*tmp_145 + tmp_13*tmp_143 + tmp_144*tmp_15)) + 0.0068572537431980923*tmp_49*(p_affine_13_0*(tmp_1*tmp_149 + tmp_150*tmp_2 + tmp_151*tmp_6) + p_affine_13_1*(tmp_14*tmp_149 + tmp_150*tmp_7 + tmp_151*tmp_4) + p_affine_13_2*(tmp_11*tmp_151 + tmp_13*tmp_149 + tmp_15*tmp_150)) + 0.037198804536718075*tmp_49*(p_affine_13_0*(tmp_1*tmp_155 + tmp_156*tmp_2 + tmp_157*tmp_6) + p_affine_13_1*(tmp_14*tmp_155 + tmp_156*tmp_7 + tmp_157*tmp_4) + p_affine_13_2*(tmp_11*tmp_157 + tmp_13*tmp_155 + tmp_15*tmp_156)) + 0.042507265838595799*tmp_49*(p_affine_13_0*(tmp_1*tmp_161 + tmp_162*tmp_2 + tmp_163*tmp_6) + p_affine_13_1*(tmp_14*tmp_161 + tmp_162*tmp_7 + tmp_163*tmp_4) + p_affine_13_2*(tmp_11*tmp_163 + tmp_13*tmp_161 + tmp_15*tmp_162)) + 0.019202922745021479*tmp_49*(p_affine_13_0*(tmp_1*tmp_167 + tmp_168*tmp_2 + tmp_169*tmp_6) + p_affine_13_1*(tmp_14*tmp_167 + tmp_168*tmp_7 + tmp_169*tmp_4) + p_affine_13_2*(tmp_11*tmp_169 + tmp_13*tmp_167 + tmp_15*tmp_168)) + 0.0068572537431980923*tmp_49*(p_affine_13_0*(tmp_1*tmp_37 + tmp_2*tmp_41 + tmp_45*tmp_6) + p_affine_13_1*(tmp_14*tmp_37 + tmp_4*tmp_45 + tmp_41*tmp_7) + p_affine_13_2*(tmp_11*tmp_45 + tmp_13*tmp_37 + tmp_15*tmp_41)) + 0.037198804536718075*tmp_49*(p_affine_13_0*(tmp_1*tmp_53 + tmp_2*tmp_54 + tmp_55*tmp_6) + p_affine_13_1*(tmp_14*tmp_53 + tmp_4*tmp_55 + tmp_54*tmp_7) + p_affine_13_2*(tmp_11*tmp_55 + tmp_13*tmp_53 + tmp_15*tmp_54)) + 0.020848748529055869*tmp_49*(p_affine_13_0*(tmp_1*tmp_59 + tmp_2*tmp_60 + tmp_6*tmp_61) + p_affine_13_1*(tmp_14*tmp_59 + tmp_4*tmp_61 + tmp_60*tmp_7) + p_affine_13_2*(tmp_11*tmp_61 + tmp_13*tmp_59 + tmp_15*tmp_60)) + 0.019202922745021479*tmp_49*(p_affine_13_0*(tmp_1*tmp_65 + tmp_2*tmp_66 + tmp_6*tmp_67) + p_affine_13_1*(tmp_14*tmp_65 + tmp_4*tmp_67 + tmp_66*tmp_7) + p_affine_13_2*(tmp_11*tmp_67 + tmp_13*tmp_65 + tmp_15*tmp_66)) + 0.020848748529055869*tmp_49*(p_affine_13_0*(tmp_1*tmp_71 + tmp_2*tmp_72 + tmp_6*tmp_73) + p_affine_13_1*(tmp_14*tmp_71 + tmp_4*tmp_73 + tmp_7*tmp_72) + p_affine_13_2*(tmp_11*tmp_73 + tmp_13*tmp_71 + tmp_15*tmp_72)) + 0.019202922745021479*tmp_49*(p_affine_13_0*(tmp_1*tmp_77 + tmp_2*tmp_78 + tmp_6*tmp_79) + p_affine_13_1*(tmp_14*tmp_77 + tmp_4*tmp_79 + tmp_7*tmp_78) + p_affine_13_2*(tmp_11*tmp_79 + tmp_13*tmp_77 + tmp_15*tmp_78)) + 0.020848748529055869*tmp_49*(p_affine_13_0*(tmp_1*tmp_83 + tmp_2*tmp_84 + tmp_6*tmp_85) + p_affine_13_1*(tmp_14*tmp_83 + tmp_4*tmp_85 + tmp_7*tmp_84) + p_affine_13_2*(tmp_11*tmp_85 + tmp_13*tmp_83 + tmp_15*tmp_84)) + 0.019202922745021479*tmp_49*(p_affine_13_0*(tmp_1*tmp_89 + tmp_2*tmp_90 + tmp_6*tmp_91) + p_affine_13_1*(tmp_14*tmp_89 + tmp_4*tmp_91 + tmp_7*tmp_90) + p_affine_13_2*(tmp_11*tmp_91 + tmp_13*tmp_89 + tmp_15*tmp_90)) + 0.020848748529055869*tmp_49*(p_affine_13_0*(tmp_1*tmp_95 + tmp_2*tmp_96 + tmp_6*tmp_97) + p_affine_13_1*(tmp_14*tmp_95 + tmp_4*tmp_97 + tmp_7*tmp_96) + p_affine_13_2*(tmp_11*tmp_97 + tmp_13*tmp_95 + tmp_15*tmp_96));
      elMat( 0, 0) = a_0_0;
   }

public:



};


} //eg
} // dg
} // hyteg

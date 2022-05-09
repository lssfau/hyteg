
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

class DGDivFormP1P1_0 : public hyteg::dg::DGForm2D
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

      real_t tmp_0 = 0.091576213509770743;
      real_t tmp_1 = -p_affine_0_1;
      real_t tmp_2 = p_affine_2_1 + tmp_1;
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = 1.0 / (tmp_2*(p_affine_1_0 + tmp_3) - (p_affine_1_1 + tmp_1)*(p_affine_2_0 + tmp_3));
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = tmp_4*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_7 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_8 = tmp_7*(tmp_5 + tmp_6);
      real_t tmp_9 = 0.054975871827660928*tmp_8;
      real_t tmp_10 = 0.44594849091596489;
      real_t tmp_11 = 0.11169079483900572*tmp_8;
      real_t tmp_12 = 0.091576213509770743;
      real_t tmp_13 = 0.054975871827660928*tmp_8;
      real_t tmp_14 = 0.44594849091596489;
      real_t tmp_15 = 0.11169079483900572*tmp_8;
      real_t tmp_16 = 0.81684757298045851;
      real_t tmp_17 = 0.054975871827660928*tmp_8;
      real_t tmp_18 = 0.10810301816807022;
      real_t tmp_19 = 0.11169079483900572*tmp_8;
      real_t tmp_20 = tmp_5*tmp_7;
      real_t tmp_21 = 0.054975871827660928*tmp_0;
      real_t tmp_22 = 0.11169079483900572*tmp_10;
      real_t tmp_23 = 0.054975871827660928*tmp_12;
      real_t tmp_24 = 0.11169079483900572*tmp_14;
      real_t tmp_25 = 0.054975871827660928*tmp_16;
      real_t tmp_26 = 0.11169079483900572*tmp_18;
      real_t tmp_27 = tmp_6*tmp_7;
      real_t tmp_28 = 0.0050344821763756674;
      real_t tmp_29 = 0.049808341407659239;
      real_t tmp_30 = 0.044906907474909594;
      real_t tmp_31 = 0.012074112023687239;
      real_t tmp_32 = 0.0050344821763756674;
      real_t tmp_33 = 0.049808341407659239;
      real_t tmp_34 = 0.044906907474909594;
      real_t tmp_35 = 0.012074112023687239;
      real_t tmp_36 = 0.0050344821763756674;
      real_t tmp_37 = 0.049808341407659239;
      real_t tmp_38 = 0.0050344821763756674;
      real_t tmp_39 = 0.049808341407659239;
      real_t a_0_0 = tmp_0*tmp_9 + tmp_10*tmp_11 + tmp_12*tmp_13 + tmp_14*tmp_15 + tmp_16*tmp_17 + tmp_18*tmp_19;
      real_t a_0_1 = -tmp_20*tmp_21 - tmp_20*tmp_22 - tmp_20*tmp_23 - tmp_20*tmp_24 - tmp_20*tmp_25 - tmp_20*tmp_26;
      real_t a_0_2 = -tmp_21*tmp_27 - tmp_22*tmp_27 - tmp_23*tmp_27 - tmp_24*tmp_27 - tmp_25*tmp_27 - tmp_26*tmp_27;
      real_t a_1_0 = 0.44594849091596489*tmp_11 + 0.81684757298045851*tmp_13 + 0.10810301816807022*tmp_15 + 0.091576213509770743*tmp_17 + 0.44594849091596489*tmp_19 + 0.091576213509770743*tmp_9;
      real_t a_1_1 = -tmp_20*tmp_28 - tmp_20*tmp_29 - tmp_20*tmp_30 - tmp_20*tmp_31 - tmp_20*tmp_32 - tmp_20*tmp_33;
      real_t a_1_2 = -tmp_27*tmp_28 - tmp_27*tmp_29 - tmp_27*tmp_30 - tmp_27*tmp_31 - tmp_27*tmp_32 - tmp_27*tmp_33;
      real_t a_2_0 = 0.10810301816807022*tmp_11 + 0.091576213509770743*tmp_13 + 0.44594849091596489*tmp_15 + 0.091576213509770743*tmp_17 + 0.44594849091596489*tmp_19 + 0.81684757298045851*tmp_9;
      real_t a_2_1 = -tmp_20*tmp_34 - tmp_20*tmp_35 - tmp_20*tmp_36 - tmp_20*tmp_37 - tmp_20*tmp_38 - tmp_20*tmp_39;
      real_t a_2_2 = -tmp_27*tmp_34 - tmp_27*tmp_35 - tmp_27*tmp_36 - tmp_27*tmp_37 - tmp_27*tmp_38 - tmp_27*tmp_39;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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
      real_t tmp_7 = tmp_4*(0.21132486540518713*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.21132486540518713*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = 0.5*p_affine_10_0*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_19 = 0.5*tmp_18;
      real_t tmp_20 = tmp_4*(0.78867513459481287*tmp_5 + tmp_6);
      real_t tmp_21 = tmp_1*tmp_20;
      real_t tmp_22 = tmp_20*tmp_9;
      real_t tmp_23 = tmp_4*(0.78867513459481287*tmp_11 + tmp_12);
      real_t tmp_24 = tmp_23*tmp_3;
      real_t tmp_25 = tmp_15*tmp_23;
      real_t tmp_26 = -tmp_21 - tmp_22 - tmp_24 - tmp_25 + 1;
      real_t tmp_27 = 0.5*tmp_18;
      real_t tmp_28 = tmp_10 + tmp_14;
      real_t tmp_29 = tmp_17*tmp_19;
      real_t tmp_30 = tmp_22 + tmp_24;
      real_t tmp_31 = tmp_26*tmp_27;
      real_t tmp_32 = tmp_28*tmp_29 + tmp_30*tmp_31;
      real_t tmp_33 = tmp_16 + tmp_8;
      real_t tmp_34 = tmp_21 + tmp_25;
      real_t tmp_35 = tmp_29*tmp_33 + tmp_31*tmp_34;
      real_t tmp_36 = tmp_19*tmp_28*tmp_33 + tmp_27*tmp_30*tmp_34;
      real_t a_0_0 = (tmp_17*tmp_17)*tmp_19 + (tmp_26*tmp_26)*tmp_27;
      real_t a_0_1 = tmp_32;
      real_t a_0_2 = tmp_35;
      real_t a_1_0 = tmp_32;
      real_t a_1_1 = tmp_19*(tmp_28*tmp_28) + tmp_27*(tmp_30*tmp_30);
      real_t a_1_2 = tmp_36;
      real_t a_2_0 = tmp_35;
      real_t a_2_1 = tmp_36;
      real_t a_2_2 = tmp_19*(tmp_33*tmp_33) + tmp_27*(tmp_34*tmp_34);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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
      real_t tmp_6 = p_affine_6_1 + 0.21132486540518713*tmp_5;
      real_t tmp_7 = tmp_4*(tmp_2 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.21132486540518713*tmp_11;
      real_t tmp_13 = tmp_4*(tmp_0 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = -p_affine_0_0;
      real_t tmp_19 = p_affine_1_0 + tmp_18;
      real_t tmp_20 = -p_affine_0_1;
      real_t tmp_21 = p_affine_2_1 + tmp_20;
      real_t tmp_22 = 1.0 / (tmp_19*tmp_21 - (p_affine_1_1 + tmp_20)*(p_affine_2_0 + tmp_18));
      real_t tmp_23 = tmp_22*(tmp_20 + tmp_6);
      real_t tmp_24 = tmp_19*tmp_23;
      real_t tmp_25 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_26 = tmp_23*tmp_25;
      real_t tmp_27 = tmp_22*(tmp_12 + tmp_18);
      real_t tmp_28 = tmp_21*tmp_27;
      real_t tmp_29 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = 0.5*p_affine_10_0*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_32 = 0.5*tmp_31;
      real_t tmp_33 = tmp_32*(-tmp_24 - tmp_26 - tmp_28 - tmp_30 + 1);
      real_t tmp_34 = p_affine_6_1 + 0.78867513459481287*tmp_5;
      real_t tmp_35 = tmp_4*(tmp_2 + tmp_34);
      real_t tmp_36 = tmp_1*tmp_35;
      real_t tmp_37 = tmp_35*tmp_9;
      real_t tmp_38 = p_affine_6_0 + 0.78867513459481287*tmp_11;
      real_t tmp_39 = tmp_4*(tmp_0 + tmp_38);
      real_t tmp_40 = tmp_3*tmp_39;
      real_t tmp_41 = tmp_15*tmp_39;
      real_t tmp_42 = -tmp_36 - tmp_37 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = tmp_22*(tmp_20 + tmp_34);
      real_t tmp_44 = tmp_19*tmp_43;
      real_t tmp_45 = tmp_25*tmp_43;
      real_t tmp_46 = tmp_22*(tmp_18 + tmp_38);
      real_t tmp_47 = tmp_21*tmp_46;
      real_t tmp_48 = tmp_29*tmp_46;
      real_t tmp_49 = 0.5*tmp_31;
      real_t tmp_50 = tmp_49*(-tmp_44 - tmp_45 - tmp_47 - tmp_48 + 1);
      real_t tmp_51 = tmp_10 + tmp_14;
      real_t tmp_52 = tmp_37 + tmp_40;
      real_t tmp_53 = tmp_16 + tmp_8;
      real_t tmp_54 = tmp_36 + tmp_41;
      real_t tmp_55 = tmp_32*(tmp_26 + tmp_28);
      real_t tmp_56 = tmp_49*(tmp_45 + tmp_47);
      real_t tmp_57 = tmp_32*(tmp_24 + tmp_30);
      real_t tmp_58 = tmp_49*(tmp_44 + tmp_48);
      real_t a_0_0 = -tmp_17*tmp_33 - tmp_42*tmp_50;
      real_t a_0_1 = -tmp_33*tmp_51 - tmp_50*tmp_52;
      real_t a_0_2 = -tmp_33*tmp_53 - tmp_50*tmp_54;
      real_t a_1_0 = -tmp_17*tmp_55 - tmp_42*tmp_56;
      real_t a_1_1 = -tmp_51*tmp_55 - tmp_52*tmp_56;
      real_t a_1_2 = -tmp_53*tmp_55 - tmp_54*tmp_56;
      real_t a_2_0 = -tmp_17*tmp_57 - tmp_42*tmp_58;
      real_t a_2_1 = -tmp_51*tmp_57 - tmp_52*tmp_58;
      real_t a_2_2 = -tmp_53*tmp_57 - tmp_54*tmp_58;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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
      real_t tmp_7 = tmp_4*(0.21132486540518713*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.21132486540518713*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = p_affine_10_0*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_19 = 0.5*tmp_18;
      real_t tmp_20 = tmp_4*(0.78867513459481287*tmp_5 + tmp_6);
      real_t tmp_21 = tmp_1*tmp_20;
      real_t tmp_22 = tmp_20*tmp_9;
      real_t tmp_23 = tmp_4*(0.78867513459481287*tmp_11 + tmp_12);
      real_t tmp_24 = tmp_23*tmp_3;
      real_t tmp_25 = tmp_15*tmp_23;
      real_t tmp_26 = -tmp_21 - tmp_22 - tmp_24 - tmp_25 + 1;
      real_t tmp_27 = 0.5*tmp_18;
      real_t tmp_28 = tmp_10 + tmp_14;
      real_t tmp_29 = tmp_17*tmp_19;
      real_t tmp_30 = tmp_22 + tmp_24;
      real_t tmp_31 = tmp_26*tmp_27;
      real_t tmp_32 = tmp_28*tmp_29 + tmp_30*tmp_31;
      real_t tmp_33 = tmp_16 + tmp_8;
      real_t tmp_34 = tmp_21 + tmp_25;
      real_t tmp_35 = tmp_29*tmp_33 + tmp_31*tmp_34;
      real_t tmp_36 = tmp_19*tmp_28*tmp_33 + tmp_27*tmp_30*tmp_34;
      real_t a_0_0 = (tmp_17*tmp_17)*tmp_19 + (tmp_26*tmp_26)*tmp_27;
      real_t a_0_1 = tmp_32;
      real_t a_0_2 = tmp_35;
      real_t a_1_0 = tmp_32;
      real_t a_1_1 = tmp_19*(tmp_28*tmp_28) + tmp_27*(tmp_30*tmp_30);
      real_t a_1_2 = tmp_36;
      real_t a_2_0 = tmp_35;
      real_t a_2_1 = tmp_36;
      real_t a_2_2 = tmp_19*(tmp_33*tmp_33) + tmp_27*(tmp_34*tmp_34);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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




class DGDivFormP1P1_1 : public hyteg::dg::DGForm2D
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

      real_t tmp_0 = 0.091576213509770743;
      real_t tmp_1 = -p_affine_0_0;
      real_t tmp_2 = p_affine_1_0 + tmp_1;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = 1.0 / (tmp_2*(p_affine_2_1 + tmp_3) - (p_affine_1_1 + tmp_3)*(p_affine_2_0 + tmp_1));
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = tmp_4*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_7 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_8 = tmp_7*(tmp_5 + tmp_6);
      real_t tmp_9 = 0.054975871827660928*tmp_8;
      real_t tmp_10 = 0.44594849091596489;
      real_t tmp_11 = 0.11169079483900572*tmp_8;
      real_t tmp_12 = 0.091576213509770743;
      real_t tmp_13 = 0.054975871827660928*tmp_8;
      real_t tmp_14 = 0.44594849091596489;
      real_t tmp_15 = 0.11169079483900572*tmp_8;
      real_t tmp_16 = 0.81684757298045851;
      real_t tmp_17 = 0.054975871827660928*tmp_8;
      real_t tmp_18 = 0.10810301816807022;
      real_t tmp_19 = 0.11169079483900572*tmp_8;
      real_t tmp_20 = tmp_6*tmp_7;
      real_t tmp_21 = 0.054975871827660928*tmp_0;
      real_t tmp_22 = 0.11169079483900572*tmp_10;
      real_t tmp_23 = 0.054975871827660928*tmp_12;
      real_t tmp_24 = 0.11169079483900572*tmp_14;
      real_t tmp_25 = 0.054975871827660928*tmp_16;
      real_t tmp_26 = 0.11169079483900572*tmp_18;
      real_t tmp_27 = tmp_5*tmp_7;
      real_t tmp_28 = 0.0050344821763756674;
      real_t tmp_29 = 0.049808341407659239;
      real_t tmp_30 = 0.044906907474909594;
      real_t tmp_31 = 0.012074112023687239;
      real_t tmp_32 = 0.0050344821763756674;
      real_t tmp_33 = 0.049808341407659239;
      real_t tmp_34 = 0.044906907474909594;
      real_t tmp_35 = 0.012074112023687239;
      real_t tmp_36 = 0.0050344821763756674;
      real_t tmp_37 = 0.049808341407659239;
      real_t tmp_38 = 0.0050344821763756674;
      real_t tmp_39 = 0.049808341407659239;
      real_t a_0_0 = tmp_0*tmp_9 + tmp_10*tmp_11 + tmp_12*tmp_13 + tmp_14*tmp_15 + tmp_16*tmp_17 + tmp_18*tmp_19;
      real_t a_0_1 = -tmp_20*tmp_21 - tmp_20*tmp_22 - tmp_20*tmp_23 - tmp_20*tmp_24 - tmp_20*tmp_25 - tmp_20*tmp_26;
      real_t a_0_2 = -tmp_21*tmp_27 - tmp_22*tmp_27 - tmp_23*tmp_27 - tmp_24*tmp_27 - tmp_25*tmp_27 - tmp_26*tmp_27;
      real_t a_1_0 = 0.44594849091596489*tmp_11 + 0.81684757298045851*tmp_13 + 0.10810301816807022*tmp_15 + 0.091576213509770743*tmp_17 + 0.44594849091596489*tmp_19 + 0.091576213509770743*tmp_9;
      real_t a_1_1 = -tmp_20*tmp_28 - tmp_20*tmp_29 - tmp_20*tmp_30 - tmp_20*tmp_31 - tmp_20*tmp_32 - tmp_20*tmp_33;
      real_t a_1_2 = -tmp_27*tmp_28 - tmp_27*tmp_29 - tmp_27*tmp_30 - tmp_27*tmp_31 - tmp_27*tmp_32 - tmp_27*tmp_33;
      real_t a_2_0 = 0.10810301816807022*tmp_11 + 0.091576213509770743*tmp_13 + 0.44594849091596489*tmp_15 + 0.091576213509770743*tmp_17 + 0.44594849091596489*tmp_19 + 0.81684757298045851*tmp_9;
      real_t a_2_1 = -tmp_20*tmp_34 - tmp_20*tmp_35 - tmp_20*tmp_36 - tmp_20*tmp_37 - tmp_20*tmp_38 - tmp_20*tmp_39;
      real_t a_2_2 = -tmp_27*tmp_34 - tmp_27*tmp_35 - tmp_27*tmp_36 - tmp_27*tmp_37 - tmp_27*tmp_38 - tmp_27*tmp_39;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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
      real_t tmp_7 = tmp_4*(0.21132486540518713*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.21132486540518713*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = 0.5*p_affine_10_1*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_19 = 0.5*tmp_18;
      real_t tmp_20 = tmp_4*(0.78867513459481287*tmp_5 + tmp_6);
      real_t tmp_21 = tmp_1*tmp_20;
      real_t tmp_22 = tmp_20*tmp_9;
      real_t tmp_23 = tmp_4*(0.78867513459481287*tmp_11 + tmp_12);
      real_t tmp_24 = tmp_23*tmp_3;
      real_t tmp_25 = tmp_15*tmp_23;
      real_t tmp_26 = -tmp_21 - tmp_22 - tmp_24 - tmp_25 + 1;
      real_t tmp_27 = 0.5*tmp_18;
      real_t tmp_28 = tmp_10 + tmp_14;
      real_t tmp_29 = tmp_17*tmp_19;
      real_t tmp_30 = tmp_22 + tmp_24;
      real_t tmp_31 = tmp_26*tmp_27;
      real_t tmp_32 = tmp_28*tmp_29 + tmp_30*tmp_31;
      real_t tmp_33 = tmp_16 + tmp_8;
      real_t tmp_34 = tmp_21 + tmp_25;
      real_t tmp_35 = tmp_29*tmp_33 + tmp_31*tmp_34;
      real_t tmp_36 = tmp_19*tmp_28*tmp_33 + tmp_27*tmp_30*tmp_34;
      real_t a_0_0 = (tmp_17*tmp_17)*tmp_19 + (tmp_26*tmp_26)*tmp_27;
      real_t a_0_1 = tmp_32;
      real_t a_0_2 = tmp_35;
      real_t a_1_0 = tmp_32;
      real_t a_1_1 = tmp_19*(tmp_28*tmp_28) + tmp_27*(tmp_30*tmp_30);
      real_t a_1_2 = tmp_36;
      real_t a_2_0 = tmp_35;
      real_t a_2_1 = tmp_36;
      real_t a_2_2 = tmp_19*(tmp_33*tmp_33) + tmp_27*(tmp_34*tmp_34);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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
      real_t tmp_6 = p_affine_6_1 + 0.21132486540518713*tmp_5;
      real_t tmp_7 = tmp_4*(tmp_2 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.21132486540518713*tmp_11;
      real_t tmp_13 = tmp_4*(tmp_0 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = -p_affine_0_0;
      real_t tmp_19 = p_affine_1_0 + tmp_18;
      real_t tmp_20 = -p_affine_0_1;
      real_t tmp_21 = p_affine_2_1 + tmp_20;
      real_t tmp_22 = 1.0 / (tmp_19*tmp_21 - (p_affine_1_1 + tmp_20)*(p_affine_2_0 + tmp_18));
      real_t tmp_23 = tmp_22*(tmp_20 + tmp_6);
      real_t tmp_24 = tmp_19*tmp_23;
      real_t tmp_25 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_26 = tmp_23*tmp_25;
      real_t tmp_27 = tmp_22*(tmp_12 + tmp_18);
      real_t tmp_28 = tmp_21*tmp_27;
      real_t tmp_29 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = 0.5*p_affine_10_1*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_32 = 0.5*tmp_31;
      real_t tmp_33 = tmp_32*(-tmp_24 - tmp_26 - tmp_28 - tmp_30 + 1);
      real_t tmp_34 = p_affine_6_1 + 0.78867513459481287*tmp_5;
      real_t tmp_35 = tmp_4*(tmp_2 + tmp_34);
      real_t tmp_36 = tmp_1*tmp_35;
      real_t tmp_37 = tmp_35*tmp_9;
      real_t tmp_38 = p_affine_6_0 + 0.78867513459481287*tmp_11;
      real_t tmp_39 = tmp_4*(tmp_0 + tmp_38);
      real_t tmp_40 = tmp_3*tmp_39;
      real_t tmp_41 = tmp_15*tmp_39;
      real_t tmp_42 = -tmp_36 - tmp_37 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = tmp_22*(tmp_20 + tmp_34);
      real_t tmp_44 = tmp_19*tmp_43;
      real_t tmp_45 = tmp_25*tmp_43;
      real_t tmp_46 = tmp_22*(tmp_18 + tmp_38);
      real_t tmp_47 = tmp_21*tmp_46;
      real_t tmp_48 = tmp_29*tmp_46;
      real_t tmp_49 = 0.5*tmp_31;
      real_t tmp_50 = tmp_49*(-tmp_44 - tmp_45 - tmp_47 - tmp_48 + 1);
      real_t tmp_51 = tmp_10 + tmp_14;
      real_t tmp_52 = tmp_37 + tmp_40;
      real_t tmp_53 = tmp_16 + tmp_8;
      real_t tmp_54 = tmp_36 + tmp_41;
      real_t tmp_55 = tmp_32*(tmp_26 + tmp_28);
      real_t tmp_56 = tmp_49*(tmp_45 + tmp_47);
      real_t tmp_57 = tmp_32*(tmp_24 + tmp_30);
      real_t tmp_58 = tmp_49*(tmp_44 + tmp_48);
      real_t a_0_0 = -tmp_17*tmp_33 - tmp_42*tmp_50;
      real_t a_0_1 = -tmp_33*tmp_51 - tmp_50*tmp_52;
      real_t a_0_2 = -tmp_33*tmp_53 - tmp_50*tmp_54;
      real_t a_1_0 = -tmp_17*tmp_55 - tmp_42*tmp_56;
      real_t a_1_1 = -tmp_51*tmp_55 - tmp_52*tmp_56;
      real_t a_1_2 = -tmp_53*tmp_55 - tmp_54*tmp_56;
      real_t a_2_0 = -tmp_17*tmp_57 - tmp_42*tmp_58;
      real_t a_2_1 = -tmp_51*tmp_57 - tmp_52*tmp_58;
      real_t a_2_2 = -tmp_53*tmp_57 - tmp_54*tmp_58;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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
      real_t tmp_7 = tmp_4*(0.21132486540518713*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.21132486540518713*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = p_affine_10_1*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_19 = 0.5*tmp_18;
      real_t tmp_20 = tmp_4*(0.78867513459481287*tmp_5 + tmp_6);
      real_t tmp_21 = tmp_1*tmp_20;
      real_t tmp_22 = tmp_20*tmp_9;
      real_t tmp_23 = tmp_4*(0.78867513459481287*tmp_11 + tmp_12);
      real_t tmp_24 = tmp_23*tmp_3;
      real_t tmp_25 = tmp_15*tmp_23;
      real_t tmp_26 = -tmp_21 - tmp_22 - tmp_24 - tmp_25 + 1;
      real_t tmp_27 = 0.5*tmp_18;
      real_t tmp_28 = tmp_10 + tmp_14;
      real_t tmp_29 = tmp_17*tmp_19;
      real_t tmp_30 = tmp_22 + tmp_24;
      real_t tmp_31 = tmp_26*tmp_27;
      real_t tmp_32 = tmp_28*tmp_29 + tmp_30*tmp_31;
      real_t tmp_33 = tmp_16 + tmp_8;
      real_t tmp_34 = tmp_21 + tmp_25;
      real_t tmp_35 = tmp_29*tmp_33 + tmp_31*tmp_34;
      real_t tmp_36 = tmp_19*tmp_28*tmp_33 + tmp_27*tmp_30*tmp_34;
      real_t a_0_0 = (tmp_17*tmp_17)*tmp_19 + (tmp_26*tmp_26)*tmp_27;
      real_t a_0_1 = tmp_32;
      real_t a_0_2 = tmp_35;
      real_t a_1_0 = tmp_32;
      real_t a_1_1 = tmp_19*(tmp_28*tmp_28) + tmp_27*(tmp_30*tmp_30);
      real_t a_1_2 = tmp_36;
      real_t a_2_0 = tmp_35;
      real_t a_2_1 = tmp_36;
      real_t a_2_2 = tmp_19*(tmp_33*tmp_33) + tmp_27*(tmp_34*tmp_34);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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




class DGDivFormP1EDG : public hyteg::dg::DGForm2D
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = (-tmp_1*tmp_4 - tmp_3*tmp_4)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_6 = 0.054975871827660928*tmp_5;
      real_t tmp_7 = 0.11169079483900572*tmp_5;
      real_t tmp_8 = 0.054975871827660928*tmp_5;
      real_t tmp_9 = 0.11169079483900572*tmp_5;
      real_t tmp_10 = 0.054975871827660928*tmp_5;
      real_t tmp_11 = 0.11169079483900572*tmp_5;
      real_t a_0_0 = 0.81684757298045851*tmp_10 + 0.10810301816807022*tmp_11 + 0.091576213509770743*tmp_6 + 0.44594849091596489*tmp_7 + 0.091576213509770743*tmp_8 + 0.44594849091596489*tmp_9;
      real_t a_1_0 = 0.091576213509770743*tmp_10 + 0.44594849091596489*tmp_11 + 0.091576213509770743*tmp_6 + 0.44594849091596489*tmp_7 + 0.81684757298045851*tmp_8 + 0.10810301816807022*tmp_9;
      real_t a_2_0 = 0.091576213509770743*tmp_10 + 0.44594849091596489*tmp_11 + 0.81684757298045851*tmp_6 + 0.10810301816807022*tmp_7 + 0.091576213509770743*tmp_8 + 0.44594849091596489*tmp_9;
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
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = 1.0 / (tmp_1*tmp_3 - tmp_4*tmp_5);
      real_t tmp_7 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_8 = p_affine_6_1 + tmp_2;
      real_t tmp_9 = tmp_6*(0.21132486540518713*tmp_7 + tmp_8);
      real_t tmp_10 = tmp_1*tmp_9;
      real_t tmp_11 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_12 = tmp_11*tmp_9;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_0;
      real_t tmp_15 = tmp_6*(0.21132486540518713*tmp_13 + tmp_14);
      real_t tmp_16 = tmp_15*tmp_3;
      real_t tmp_17 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_18 = tmp_15*tmp_17;
      real_t tmp_19 = -0.5*tmp_10 - 0.5*tmp_12 - 0.5*tmp_16 - 0.5*tmp_18 + 0.5;
      real_t tmp_20 = tmp_12 + tmp_16;
      real_t tmp_21 = tmp_20 - 1.0/3.0;
      real_t tmp_22 = tmp_10 + tmp_18;
      real_t tmp_23 = tmp_22 - 1.0/3.0;
      real_t tmp_24 = p_affine_10_0*(tmp_1*tmp_21 + tmp_23*tmp_4);
      real_t tmp_25 = p_affine_10_1*(tmp_21*tmp_5 + tmp_23*tmp_3);
      real_t tmp_26 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_7*tmp_7), 1.0/2.0));
      real_t tmp_27 = 0.5*tmp_26;
      real_t tmp_28 = tmp_6*(0.78867513459481287*tmp_7 + tmp_8);
      real_t tmp_29 = tmp_1*tmp_28;
      real_t tmp_30 = tmp_11*tmp_28;
      real_t tmp_31 = tmp_6*(0.78867513459481287*tmp_13 + tmp_14);
      real_t tmp_32 = tmp_3*tmp_31;
      real_t tmp_33 = tmp_17*tmp_31;
      real_t tmp_34 = -0.5*tmp_29 - 0.5*tmp_30 - 0.5*tmp_32 - 0.5*tmp_33 + 0.5;
      real_t tmp_35 = tmp_30 + tmp_32;
      real_t tmp_36 = tmp_35 - 1.0/3.0;
      real_t tmp_37 = tmp_29 + tmp_33;
      real_t tmp_38 = tmp_37 - 1.0/3.0;
      real_t tmp_39 = p_affine_10_0*(tmp_1*tmp_36 + tmp_38*tmp_4);
      real_t tmp_40 = p_affine_10_1*(tmp_3*tmp_38 + tmp_36*tmp_5);
      real_t tmp_41 = 0.5*tmp_26;
      real_t tmp_42 = 0.5*tmp_20;
      real_t tmp_43 = 0.5*tmp_35;
      real_t tmp_44 = 0.5*tmp_22;
      real_t tmp_45 = 0.5*tmp_37;
      real_t a_0_0 = tmp_27*(tmp_19*tmp_24 + tmp_19*tmp_25) + tmp_41*(tmp_34*tmp_39 + tmp_34*tmp_40);
      real_t a_1_0 = tmp_27*(tmp_24*tmp_42 + tmp_25*tmp_42) + tmp_41*(tmp_39*tmp_43 + tmp_40*tmp_43);
      real_t a_2_0 = tmp_27*(tmp_24*tmp_44 + tmp_25*tmp_44) + tmp_41*(tmp_39*tmp_45 + tmp_40*tmp_45);
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
      real_t tmp_6 = p_affine_6_1 + 0.21132486540518713*tmp_5;
      real_t tmp_7 = tmp_4*(tmp_2 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.21132486540518713*tmp_11;
      real_t tmp_13 = tmp_4*(tmp_0 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -0.5*tmp_10 - 0.5*tmp_14 - 0.5*tmp_16 - 0.5*tmp_8 + 0.5;
      real_t tmp_18 = -p_affine_3_0;
      real_t tmp_19 = p_affine_4_0 + tmp_18;
      real_t tmp_20 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_21 = -p_affine_3_1;
      real_t tmp_22 = p_affine_5_1 + tmp_21;
      real_t tmp_23 = p_affine_5_0 + tmp_18;
      real_t tmp_24 = p_affine_4_1 + tmp_21;
      real_t tmp_25 = 1.0 / (tmp_19*tmp_22 - tmp_23*tmp_24);
      real_t tmp_26 = tmp_25*(tmp_21 + tmp_6);
      real_t tmp_27 = tmp_25*(tmp_12 + tmp_18);
      real_t tmp_28 = tmp_20*tmp_26 + tmp_22*tmp_27 - 1.0/3.0;
      real_t tmp_29 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_30 = tmp_19*tmp_26 + tmp_27*tmp_29 - 1.0/3.0;
      real_t tmp_31 = p_affine_10_0*(tmp_19*tmp_28 + tmp_23*tmp_30);
      real_t tmp_32 = p_affine_10_1*(tmp_22*tmp_30 + tmp_24*tmp_28);
      real_t tmp_33 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_34 = 0.5*tmp_33;
      real_t tmp_35 = p_affine_6_1 + 0.78867513459481287*tmp_5;
      real_t tmp_36 = tmp_4*(tmp_2 + tmp_35);
      real_t tmp_37 = tmp_1*tmp_36;
      real_t tmp_38 = tmp_36*tmp_9;
      real_t tmp_39 = p_affine_6_0 + 0.78867513459481287*tmp_11;
      real_t tmp_40 = tmp_4*(tmp_0 + tmp_39);
      real_t tmp_41 = tmp_3*tmp_40;
      real_t tmp_42 = tmp_15*tmp_40;
      real_t tmp_43 = -0.5*tmp_37 - 0.5*tmp_38 - 0.5*tmp_41 - 0.5*tmp_42 + 0.5;
      real_t tmp_44 = tmp_25*(tmp_21 + tmp_35);
      real_t tmp_45 = tmp_25*(tmp_18 + tmp_39);
      real_t tmp_46 = tmp_20*tmp_44 + tmp_22*tmp_45 - 1.0/3.0;
      real_t tmp_47 = tmp_19*tmp_44 + tmp_29*tmp_45 - 1.0/3.0;
      real_t tmp_48 = p_affine_10_0*(tmp_19*tmp_46 + tmp_23*tmp_47);
      real_t tmp_49 = p_affine_10_1*(tmp_22*tmp_47 + tmp_24*tmp_46);
      real_t tmp_50 = 0.5*tmp_33;
      real_t tmp_51 = 0.5*tmp_10 + 0.5*tmp_14;
      real_t tmp_52 = 0.5*tmp_38 + 0.5*tmp_41;
      real_t tmp_53 = 0.5*tmp_16 + 0.5*tmp_8;
      real_t tmp_54 = 0.5*tmp_37 + 0.5*tmp_42;
      real_t a_0_0 = tmp_34*(-tmp_17*tmp_31 - tmp_17*tmp_32) + tmp_50*(-tmp_43*tmp_48 - tmp_43*tmp_49);
      real_t a_1_0 = tmp_34*(-tmp_31*tmp_51 - tmp_32*tmp_51) + tmp_50*(-tmp_48*tmp_52 - tmp_49*tmp_52);
      real_t a_2_0 = tmp_34*(-tmp_31*tmp_53 - tmp_32*tmp_53) + tmp_50*(-tmp_48*tmp_54 - tmp_49*tmp_54);
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
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = 1.0 / (tmp_1*tmp_3 - tmp_4*tmp_5);
      real_t tmp_7 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_8 = p_affine_6_1 + tmp_2;
      real_t tmp_9 = tmp_6*(0.21132486540518713*tmp_7 + tmp_8);
      real_t tmp_10 = tmp_1*tmp_9;
      real_t tmp_11 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_12 = tmp_11*tmp_9;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_0;
      real_t tmp_15 = tmp_6*(0.21132486540518713*tmp_13 + tmp_14);
      real_t tmp_16 = tmp_15*tmp_3;
      real_t tmp_17 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_18 = tmp_15*tmp_17;
      real_t tmp_19 = -tmp_10 - tmp_12 - tmp_16 - tmp_18 + 1;
      real_t tmp_20 = tmp_12 + tmp_16;
      real_t tmp_21 = tmp_20 - 1.0/3.0;
      real_t tmp_22 = tmp_10 + tmp_18;
      real_t tmp_23 = tmp_22 - 1.0/3.0;
      real_t tmp_24 = p_affine_10_0*(tmp_1*tmp_21 + tmp_23*tmp_4);
      real_t tmp_25 = p_affine_10_1*(tmp_21*tmp_5 + tmp_23*tmp_3);
      real_t tmp_26 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_7*tmp_7), 1.0/2.0));
      real_t tmp_27 = 0.5*tmp_26;
      real_t tmp_28 = tmp_6*(0.78867513459481287*tmp_7 + tmp_8);
      real_t tmp_29 = tmp_1*tmp_28;
      real_t tmp_30 = tmp_11*tmp_28;
      real_t tmp_31 = tmp_6*(0.78867513459481287*tmp_13 + tmp_14);
      real_t tmp_32 = tmp_3*tmp_31;
      real_t tmp_33 = tmp_17*tmp_31;
      real_t tmp_34 = -tmp_29 - tmp_30 - tmp_32 - tmp_33 + 1;
      real_t tmp_35 = tmp_30 + tmp_32;
      real_t tmp_36 = tmp_35 - 1.0/3.0;
      real_t tmp_37 = tmp_29 + tmp_33;
      real_t tmp_38 = tmp_37 - 1.0/3.0;
      real_t tmp_39 = p_affine_10_0*(tmp_1*tmp_36 + tmp_38*tmp_4);
      real_t tmp_40 = p_affine_10_1*(tmp_3*tmp_38 + tmp_36*tmp_5);
      real_t tmp_41 = 0.5*tmp_26;
      real_t a_0_0 = tmp_27*(tmp_19*tmp_24 + tmp_19*tmp_25) + tmp_41*(tmp_34*tmp_39 + tmp_34*tmp_40);
      real_t a_1_0 = tmp_27*(tmp_20*tmp_24 + tmp_20*tmp_25) + tmp_41*(tmp_35*tmp_39 + tmp_35*tmp_40);
      real_t a_2_0 = tmp_27*(tmp_22*tmp_24 + tmp_22*tmp_25) + tmp_41*(tmp_37*tmp_39 + tmp_37*tmp_40);
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

      // Does nothing.
   }
};




class DGDivtFormP1P1_0 : public hyteg::dg::DGForm2D
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

      real_t tmp_0 = 0.091576213509770743;
      real_t tmp_1 = -p_affine_0_1;
      real_t tmp_2 = p_affine_2_1 + tmp_1;
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = 1.0 / (tmp_2*(p_affine_1_0 + tmp_3) - (p_affine_1_1 + tmp_1)*(p_affine_2_0 + tmp_3));
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = tmp_4*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_7 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_8 = tmp_7*(tmp_5 + tmp_6);
      real_t tmp_9 = 0.054975871827660928*tmp_8;
      real_t tmp_10 = 0.44594849091596489;
      real_t tmp_11 = 0.11169079483900572*tmp_8;
      real_t tmp_12 = 0.091576213509770743;
      real_t tmp_13 = 0.054975871827660928*tmp_8;
      real_t tmp_14 = 0.44594849091596489;
      real_t tmp_15 = 0.11169079483900572*tmp_8;
      real_t tmp_16 = 0.81684757298045851;
      real_t tmp_17 = 0.054975871827660928*tmp_8;
      real_t tmp_18 = 0.10810301816807022;
      real_t tmp_19 = 0.11169079483900572*tmp_8;
      real_t tmp_20 = tmp_5*tmp_7;
      real_t tmp_21 = 0.054975871827660928*tmp_20;
      real_t tmp_22 = 0.11169079483900572*tmp_20;
      real_t tmp_23 = 0.054975871827660928*tmp_20;
      real_t tmp_24 = 0.11169079483900572*tmp_20;
      real_t tmp_25 = 0.054975871827660928*tmp_20;
      real_t tmp_26 = 0.11169079483900572*tmp_20;
      real_t tmp_27 = tmp_6*tmp_7;
      real_t tmp_28 = 0.054975871827660928*tmp_27;
      real_t tmp_29 = 0.11169079483900572*tmp_27;
      real_t tmp_30 = 0.054975871827660928*tmp_27;
      real_t tmp_31 = 0.11169079483900572*tmp_27;
      real_t tmp_32 = 0.054975871827660928*tmp_27;
      real_t tmp_33 = 0.11169079483900572*tmp_27;
      real_t a_0_0 = tmp_0*tmp_9 + tmp_10*tmp_11 + tmp_12*tmp_13 + tmp_14*tmp_15 + tmp_16*tmp_17 + tmp_18*tmp_19;
      real_t a_0_1 = 0.44594849091596489*tmp_11 + 0.81684757298045851*tmp_13 + 0.10810301816807022*tmp_15 + 0.091576213509770743*tmp_17 + 0.44594849091596489*tmp_19 + 0.091576213509770743*tmp_9;
      real_t a_0_2 = 0.10810301816807022*tmp_11 + 0.091576213509770743*tmp_13 + 0.44594849091596489*tmp_15 + 0.091576213509770743*tmp_17 + 0.44594849091596489*tmp_19 + 0.81684757298045851*tmp_9;
      real_t a_1_0 = -tmp_0*tmp_21 - tmp_10*tmp_22 - tmp_12*tmp_23 - tmp_14*tmp_24 - tmp_16*tmp_25 - tmp_18*tmp_26;
      real_t a_1_1 = -0.091576213509770743*tmp_21 - 0.44594849091596489*tmp_22 - 0.81684757298045851*tmp_23 - 0.10810301816807022*tmp_24 - 0.091576213509770743*tmp_25 - 0.44594849091596489*tmp_26;
      real_t a_1_2 = -0.81684757298045851*tmp_21 - 0.10810301816807022*tmp_22 - 0.091576213509770743*tmp_23 - 0.44594849091596489*tmp_24 - 0.091576213509770743*tmp_25 - 0.44594849091596489*tmp_26;
      real_t a_2_0 = -tmp_0*tmp_28 - tmp_10*tmp_29 - tmp_12*tmp_30 - tmp_14*tmp_31 - tmp_16*tmp_32 - tmp_18*tmp_33;
      real_t a_2_1 = -0.091576213509770743*tmp_28 - 0.44594849091596489*tmp_29 - 0.81684757298045851*tmp_30 - 0.10810301816807022*tmp_31 - 0.091576213509770743*tmp_32 - 0.44594849091596489*tmp_33;
      real_t a_2_2 = -0.81684757298045851*tmp_28 - 0.10810301816807022*tmp_29 - 0.091576213509770743*tmp_30 - 0.44594849091596489*tmp_31 - 0.091576213509770743*tmp_32 - 0.44594849091596489*tmp_33;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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
      real_t tmp_7 = tmp_4*(0.21132486540518713*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.21132486540518713*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = 0.5*p_affine_10_0*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_19 = 0.5*tmp_18;
      real_t tmp_20 = tmp_4*(0.78867513459481287*tmp_5 + tmp_6);
      real_t tmp_21 = tmp_1*tmp_20;
      real_t tmp_22 = tmp_20*tmp_9;
      real_t tmp_23 = tmp_4*(0.78867513459481287*tmp_11 + tmp_12);
      real_t tmp_24 = tmp_23*tmp_3;
      real_t tmp_25 = tmp_15*tmp_23;
      real_t tmp_26 = -tmp_21 - tmp_22 - tmp_24 - tmp_25 + 1;
      real_t tmp_27 = 0.5*tmp_18;
      real_t tmp_28 = tmp_10 + tmp_14;
      real_t tmp_29 = tmp_17*tmp_19;
      real_t tmp_30 = tmp_22 + tmp_24;
      real_t tmp_31 = tmp_26*tmp_27;
      real_t tmp_32 = tmp_28*tmp_29 + tmp_30*tmp_31;
      real_t tmp_33 = tmp_16 + tmp_8;
      real_t tmp_34 = tmp_21 + tmp_25;
      real_t tmp_35 = tmp_29*tmp_33 + tmp_31*tmp_34;
      real_t tmp_36 = tmp_19*tmp_28*tmp_33 + tmp_27*tmp_30*tmp_34;
      real_t a_0_0 = (tmp_17*tmp_17)*tmp_19 + (tmp_26*tmp_26)*tmp_27;
      real_t a_0_1 = tmp_32;
      real_t a_0_2 = tmp_35;
      real_t a_1_0 = tmp_32;
      real_t a_1_1 = tmp_19*(tmp_28*tmp_28) + tmp_27*(tmp_30*tmp_30);
      real_t a_1_2 = tmp_36;
      real_t a_2_0 = tmp_35;
      real_t a_2_1 = tmp_36;
      real_t a_2_2 = tmp_19*(tmp_33*tmp_33) + tmp_27*(tmp_34*tmp_34);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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
      real_t tmp_6 = p_affine_6_1 + 0.21132486540518713*tmp_5;
      real_t tmp_7 = tmp_4*(tmp_2 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.21132486540518713*tmp_11;
      real_t tmp_13 = tmp_4*(tmp_0 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = -p_affine_0_0;
      real_t tmp_19 = p_affine_1_0 + tmp_18;
      real_t tmp_20 = -p_affine_0_1;
      real_t tmp_21 = p_affine_2_1 + tmp_20;
      real_t tmp_22 = 1.0 / (tmp_19*tmp_21 - (p_affine_1_1 + tmp_20)*(p_affine_2_0 + tmp_18));
      real_t tmp_23 = tmp_22*(tmp_20 + tmp_6);
      real_t tmp_24 = tmp_19*tmp_23;
      real_t tmp_25 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_26 = tmp_23*tmp_25;
      real_t tmp_27 = tmp_22*(tmp_12 + tmp_18);
      real_t tmp_28 = tmp_21*tmp_27;
      real_t tmp_29 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = 0.5*p_affine_10_0*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_32 = 0.5*tmp_31;
      real_t tmp_33 = tmp_32*(-tmp_24 - tmp_26 - tmp_28 - tmp_30 + 1);
      real_t tmp_34 = p_affine_6_1 + 0.78867513459481287*tmp_5;
      real_t tmp_35 = tmp_4*(tmp_2 + tmp_34);
      real_t tmp_36 = tmp_1*tmp_35;
      real_t tmp_37 = tmp_35*tmp_9;
      real_t tmp_38 = p_affine_6_0 + 0.78867513459481287*tmp_11;
      real_t tmp_39 = tmp_4*(tmp_0 + tmp_38);
      real_t tmp_40 = tmp_3*tmp_39;
      real_t tmp_41 = tmp_15*tmp_39;
      real_t tmp_42 = -tmp_36 - tmp_37 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = tmp_22*(tmp_20 + tmp_34);
      real_t tmp_44 = tmp_19*tmp_43;
      real_t tmp_45 = tmp_25*tmp_43;
      real_t tmp_46 = tmp_22*(tmp_18 + tmp_38);
      real_t tmp_47 = tmp_21*tmp_46;
      real_t tmp_48 = tmp_29*tmp_46;
      real_t tmp_49 = 0.5*tmp_31;
      real_t tmp_50 = tmp_49*(-tmp_44 - tmp_45 - tmp_47 - tmp_48 + 1);
      real_t tmp_51 = tmp_10 + tmp_14;
      real_t tmp_52 = tmp_37 + tmp_40;
      real_t tmp_53 = tmp_16 + tmp_8;
      real_t tmp_54 = tmp_36 + tmp_41;
      real_t tmp_55 = tmp_32*(tmp_26 + tmp_28);
      real_t tmp_56 = tmp_49*(tmp_45 + tmp_47);
      real_t tmp_57 = tmp_32*(tmp_24 + tmp_30);
      real_t tmp_58 = tmp_49*(tmp_44 + tmp_48);
      real_t a_0_0 = tmp_17*tmp_33 + tmp_42*tmp_50;
      real_t a_0_1 = tmp_33*tmp_51 + tmp_50*tmp_52;
      real_t a_0_2 = tmp_33*tmp_53 + tmp_50*tmp_54;
      real_t a_1_0 = tmp_17*tmp_55 + tmp_42*tmp_56;
      real_t a_1_1 = tmp_51*tmp_55 + tmp_52*tmp_56;
      real_t a_1_2 = tmp_53*tmp_55 + tmp_54*tmp_56;
      real_t a_2_0 = tmp_17*tmp_57 + tmp_42*tmp_58;
      real_t a_2_1 = tmp_51*tmp_57 + tmp_52*tmp_58;
      real_t a_2_2 = tmp_53*tmp_57 + tmp_54*tmp_58;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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
      real_t tmp_7 = tmp_4*(0.21132486540518713*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.21132486540518713*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = p_affine_10_0*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_19 = 0.5*tmp_18;
      real_t tmp_20 = tmp_4*(0.78867513459481287*tmp_5 + tmp_6);
      real_t tmp_21 = tmp_1*tmp_20;
      real_t tmp_22 = tmp_20*tmp_9;
      real_t tmp_23 = tmp_4*(0.78867513459481287*tmp_11 + tmp_12);
      real_t tmp_24 = tmp_23*tmp_3;
      real_t tmp_25 = tmp_15*tmp_23;
      real_t tmp_26 = -tmp_21 - tmp_22 - tmp_24 - tmp_25 + 1;
      real_t tmp_27 = 0.5*tmp_18;
      real_t tmp_28 = tmp_10 + tmp_14;
      real_t tmp_29 = tmp_17*tmp_19;
      real_t tmp_30 = tmp_22 + tmp_24;
      real_t tmp_31 = tmp_26*tmp_27;
      real_t tmp_32 = tmp_28*tmp_29 + tmp_30*tmp_31;
      real_t tmp_33 = tmp_16 + tmp_8;
      real_t tmp_34 = tmp_21 + tmp_25;
      real_t tmp_35 = tmp_29*tmp_33 + tmp_31*tmp_34;
      real_t tmp_36 = tmp_19*tmp_28*tmp_33 + tmp_27*tmp_30*tmp_34;
      real_t a_0_0 = (tmp_17*tmp_17)*tmp_19 + (tmp_26*tmp_26)*tmp_27;
      real_t a_0_1 = tmp_32;
      real_t a_0_2 = tmp_35;
      real_t a_1_0 = tmp_32;
      real_t a_1_1 = tmp_19*(tmp_28*tmp_28) + tmp_27*(tmp_30*tmp_30);
      real_t a_1_2 = tmp_36;
      real_t a_2_0 = tmp_35;
      real_t a_2_1 = tmp_36;
      real_t a_2_2 = tmp_19*(tmp_33*tmp_33) + tmp_27*(tmp_34*tmp_34);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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




class DGDivtFormP1P1_1 : public hyteg::dg::DGForm2D
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

      real_t tmp_0 = 0.091576213509770743;
      real_t tmp_1 = -p_affine_0_0;
      real_t tmp_2 = p_affine_1_0 + tmp_1;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = 1.0 / (tmp_2*(p_affine_2_1 + tmp_3) - (p_affine_1_1 + tmp_3)*(p_affine_2_0 + tmp_1));
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = tmp_4*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_7 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_8 = tmp_7*(tmp_5 + tmp_6);
      real_t tmp_9 = 0.054975871827660928*tmp_8;
      real_t tmp_10 = 0.44594849091596489;
      real_t tmp_11 = 0.11169079483900572*tmp_8;
      real_t tmp_12 = 0.091576213509770743;
      real_t tmp_13 = 0.054975871827660928*tmp_8;
      real_t tmp_14 = 0.44594849091596489;
      real_t tmp_15 = 0.11169079483900572*tmp_8;
      real_t tmp_16 = 0.81684757298045851;
      real_t tmp_17 = 0.054975871827660928*tmp_8;
      real_t tmp_18 = 0.10810301816807022;
      real_t tmp_19 = 0.11169079483900572*tmp_8;
      real_t tmp_20 = tmp_6*tmp_7;
      real_t tmp_21 = 0.054975871827660928*tmp_20;
      real_t tmp_22 = 0.11169079483900572*tmp_20;
      real_t tmp_23 = 0.054975871827660928*tmp_20;
      real_t tmp_24 = 0.11169079483900572*tmp_20;
      real_t tmp_25 = 0.054975871827660928*tmp_20;
      real_t tmp_26 = 0.11169079483900572*tmp_20;
      real_t tmp_27 = tmp_5*tmp_7;
      real_t tmp_28 = 0.054975871827660928*tmp_27;
      real_t tmp_29 = 0.11169079483900572*tmp_27;
      real_t tmp_30 = 0.054975871827660928*tmp_27;
      real_t tmp_31 = 0.11169079483900572*tmp_27;
      real_t tmp_32 = 0.054975871827660928*tmp_27;
      real_t tmp_33 = 0.11169079483900572*tmp_27;
      real_t a_0_0 = tmp_0*tmp_9 + tmp_10*tmp_11 + tmp_12*tmp_13 + tmp_14*tmp_15 + tmp_16*tmp_17 + tmp_18*tmp_19;
      real_t a_0_1 = 0.44594849091596489*tmp_11 + 0.81684757298045851*tmp_13 + 0.10810301816807022*tmp_15 + 0.091576213509770743*tmp_17 + 0.44594849091596489*tmp_19 + 0.091576213509770743*tmp_9;
      real_t a_0_2 = 0.10810301816807022*tmp_11 + 0.091576213509770743*tmp_13 + 0.44594849091596489*tmp_15 + 0.091576213509770743*tmp_17 + 0.44594849091596489*tmp_19 + 0.81684757298045851*tmp_9;
      real_t a_1_0 = -tmp_0*tmp_21 - tmp_10*tmp_22 - tmp_12*tmp_23 - tmp_14*tmp_24 - tmp_16*tmp_25 - tmp_18*tmp_26;
      real_t a_1_1 = -0.091576213509770743*tmp_21 - 0.44594849091596489*tmp_22 - 0.81684757298045851*tmp_23 - 0.10810301816807022*tmp_24 - 0.091576213509770743*tmp_25 - 0.44594849091596489*tmp_26;
      real_t a_1_2 = -0.81684757298045851*tmp_21 - 0.10810301816807022*tmp_22 - 0.091576213509770743*tmp_23 - 0.44594849091596489*tmp_24 - 0.091576213509770743*tmp_25 - 0.44594849091596489*tmp_26;
      real_t a_2_0 = -tmp_0*tmp_28 - tmp_10*tmp_29 - tmp_12*tmp_30 - tmp_14*tmp_31 - tmp_16*tmp_32 - tmp_18*tmp_33;
      real_t a_2_1 = -0.091576213509770743*tmp_28 - 0.44594849091596489*tmp_29 - 0.81684757298045851*tmp_30 - 0.10810301816807022*tmp_31 - 0.091576213509770743*tmp_32 - 0.44594849091596489*tmp_33;
      real_t a_2_2 = -0.81684757298045851*tmp_28 - 0.10810301816807022*tmp_29 - 0.091576213509770743*tmp_30 - 0.44594849091596489*tmp_31 - 0.091576213509770743*tmp_32 - 0.44594849091596489*tmp_33;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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
      real_t tmp_7 = tmp_4*(0.21132486540518713*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.21132486540518713*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = 0.5*p_affine_10_1*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_19 = 0.5*tmp_18;
      real_t tmp_20 = tmp_4*(0.78867513459481287*tmp_5 + tmp_6);
      real_t tmp_21 = tmp_1*tmp_20;
      real_t tmp_22 = tmp_20*tmp_9;
      real_t tmp_23 = tmp_4*(0.78867513459481287*tmp_11 + tmp_12);
      real_t tmp_24 = tmp_23*tmp_3;
      real_t tmp_25 = tmp_15*tmp_23;
      real_t tmp_26 = -tmp_21 - tmp_22 - tmp_24 - tmp_25 + 1;
      real_t tmp_27 = 0.5*tmp_18;
      real_t tmp_28 = tmp_10 + tmp_14;
      real_t tmp_29 = tmp_17*tmp_19;
      real_t tmp_30 = tmp_22 + tmp_24;
      real_t tmp_31 = tmp_26*tmp_27;
      real_t tmp_32 = tmp_28*tmp_29 + tmp_30*tmp_31;
      real_t tmp_33 = tmp_16 + tmp_8;
      real_t tmp_34 = tmp_21 + tmp_25;
      real_t tmp_35 = tmp_29*tmp_33 + tmp_31*tmp_34;
      real_t tmp_36 = tmp_19*tmp_28*tmp_33 + tmp_27*tmp_30*tmp_34;
      real_t a_0_0 = (tmp_17*tmp_17)*tmp_19 + (tmp_26*tmp_26)*tmp_27;
      real_t a_0_1 = tmp_32;
      real_t a_0_2 = tmp_35;
      real_t a_1_0 = tmp_32;
      real_t a_1_1 = tmp_19*(tmp_28*tmp_28) + tmp_27*(tmp_30*tmp_30);
      real_t a_1_2 = tmp_36;
      real_t a_2_0 = tmp_35;
      real_t a_2_1 = tmp_36;
      real_t a_2_2 = tmp_19*(tmp_33*tmp_33) + tmp_27*(tmp_34*tmp_34);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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
      real_t tmp_6 = p_affine_6_1 + 0.21132486540518713*tmp_5;
      real_t tmp_7 = tmp_4*(tmp_2 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.21132486540518713*tmp_11;
      real_t tmp_13 = tmp_4*(tmp_0 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = -p_affine_0_0;
      real_t tmp_19 = p_affine_1_0 + tmp_18;
      real_t tmp_20 = -p_affine_0_1;
      real_t tmp_21 = p_affine_2_1 + tmp_20;
      real_t tmp_22 = 1.0 / (tmp_19*tmp_21 - (p_affine_1_1 + tmp_20)*(p_affine_2_0 + tmp_18));
      real_t tmp_23 = tmp_22*(tmp_20 + tmp_6);
      real_t tmp_24 = tmp_19*tmp_23;
      real_t tmp_25 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_26 = tmp_23*tmp_25;
      real_t tmp_27 = tmp_22*(tmp_12 + tmp_18);
      real_t tmp_28 = tmp_21*tmp_27;
      real_t tmp_29 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = 0.5*p_affine_10_1*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_32 = 0.5*tmp_31;
      real_t tmp_33 = tmp_32*(-tmp_24 - tmp_26 - tmp_28 - tmp_30 + 1);
      real_t tmp_34 = p_affine_6_1 + 0.78867513459481287*tmp_5;
      real_t tmp_35 = tmp_4*(tmp_2 + tmp_34);
      real_t tmp_36 = tmp_1*tmp_35;
      real_t tmp_37 = tmp_35*tmp_9;
      real_t tmp_38 = p_affine_6_0 + 0.78867513459481287*tmp_11;
      real_t tmp_39 = tmp_4*(tmp_0 + tmp_38);
      real_t tmp_40 = tmp_3*tmp_39;
      real_t tmp_41 = tmp_15*tmp_39;
      real_t tmp_42 = -tmp_36 - tmp_37 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = tmp_22*(tmp_20 + tmp_34);
      real_t tmp_44 = tmp_19*tmp_43;
      real_t tmp_45 = tmp_25*tmp_43;
      real_t tmp_46 = tmp_22*(tmp_18 + tmp_38);
      real_t tmp_47 = tmp_21*tmp_46;
      real_t tmp_48 = tmp_29*tmp_46;
      real_t tmp_49 = 0.5*tmp_31;
      real_t tmp_50 = tmp_49*(-tmp_44 - tmp_45 - tmp_47 - tmp_48 + 1);
      real_t tmp_51 = tmp_10 + tmp_14;
      real_t tmp_52 = tmp_37 + tmp_40;
      real_t tmp_53 = tmp_16 + tmp_8;
      real_t tmp_54 = tmp_36 + tmp_41;
      real_t tmp_55 = tmp_32*(tmp_26 + tmp_28);
      real_t tmp_56 = tmp_49*(tmp_45 + tmp_47);
      real_t tmp_57 = tmp_32*(tmp_24 + tmp_30);
      real_t tmp_58 = tmp_49*(tmp_44 + tmp_48);
      real_t a_0_0 = tmp_17*tmp_33 + tmp_42*tmp_50;
      real_t a_0_1 = tmp_33*tmp_51 + tmp_50*tmp_52;
      real_t a_0_2 = tmp_33*tmp_53 + tmp_50*tmp_54;
      real_t a_1_0 = tmp_17*tmp_55 + tmp_42*tmp_56;
      real_t a_1_1 = tmp_51*tmp_55 + tmp_52*tmp_56;
      real_t a_1_2 = tmp_53*tmp_55 + tmp_54*tmp_56;
      real_t a_2_0 = tmp_17*tmp_57 + tmp_42*tmp_58;
      real_t a_2_1 = tmp_51*tmp_57 + tmp_52*tmp_58;
      real_t a_2_2 = tmp_53*tmp_57 + tmp_54*tmp_58;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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
      real_t tmp_7 = tmp_4*(0.21132486540518713*tmp_5 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4*(0.21132486540518713*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = p_affine_10_1*std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_19 = 0.5*tmp_18;
      real_t tmp_20 = tmp_4*(0.78867513459481287*tmp_5 + tmp_6);
      real_t tmp_21 = tmp_1*tmp_20;
      real_t tmp_22 = tmp_20*tmp_9;
      real_t tmp_23 = tmp_4*(0.78867513459481287*tmp_11 + tmp_12);
      real_t tmp_24 = tmp_23*tmp_3;
      real_t tmp_25 = tmp_15*tmp_23;
      real_t tmp_26 = -tmp_21 - tmp_22 - tmp_24 - tmp_25 + 1;
      real_t tmp_27 = 0.5*tmp_18;
      real_t tmp_28 = tmp_10 + tmp_14;
      real_t tmp_29 = tmp_17*tmp_19;
      real_t tmp_30 = tmp_22 + tmp_24;
      real_t tmp_31 = tmp_26*tmp_27;
      real_t tmp_32 = tmp_28*tmp_29 + tmp_30*tmp_31;
      real_t tmp_33 = tmp_16 + tmp_8;
      real_t tmp_34 = tmp_21 + tmp_25;
      real_t tmp_35 = tmp_29*tmp_33 + tmp_31*tmp_34;
      real_t tmp_36 = tmp_19*tmp_28*tmp_33 + tmp_27*tmp_30*tmp_34;
      real_t a_0_0 = (tmp_17*tmp_17)*tmp_19 + (tmp_26*tmp_26)*tmp_27;
      real_t a_0_1 = tmp_32;
      real_t a_0_2 = tmp_35;
      real_t a_1_0 = tmp_32;
      real_t a_1_1 = tmp_19*(tmp_28*tmp_28) + tmp_27*(tmp_30*tmp_30);
      real_t a_1_2 = tmp_36;
      real_t a_2_0 = tmp_35;
      real_t a_2_1 = tmp_36;
      real_t a_2_2 = tmp_19*(tmp_33*tmp_33) + tmp_27*(tmp_34*tmp_34);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 1, 0) = a_1_0;
      elMat( 1, 1) = a_1_1;
      elMat( 1, 2) = a_1_2;
      elMat( 2, 0) = a_2_0;
      elMat( 2, 1) = a_2_1;
      elMat( 2, 2) = a_2_2;
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




class DGDivtFormEDGP1 : public hyteg::dg::DGForm2D
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = (-tmp_1*tmp_4 - tmp_3*tmp_4)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_6 = 0.054975871827660928*tmp_5;
      real_t tmp_7 = 0.11169079483900572*tmp_5;
      real_t tmp_8 = 0.054975871827660928*tmp_5;
      real_t tmp_9 = 0.11169079483900572*tmp_5;
      real_t tmp_10 = 0.054975871827660928*tmp_5;
      real_t tmp_11 = 0.11169079483900572*tmp_5;
      real_t a_0_0 = 0.81684757298045851*tmp_10 + 0.10810301816807022*tmp_11 + 0.091576213509770743*tmp_6 + 0.44594849091596489*tmp_7 + 0.091576213509770743*tmp_8 + 0.44594849091596489*tmp_9;
      real_t a_0_1 = 0.091576213509770743*tmp_10 + 0.44594849091596489*tmp_11 + 0.091576213509770743*tmp_6 + 0.44594849091596489*tmp_7 + 0.81684757298045851*tmp_8 + 0.10810301816807022*tmp_9;
      real_t a_0_2 = 0.091576213509770743*tmp_10 + 0.44594849091596489*tmp_11 + 0.81684757298045851*tmp_6 + 0.10810301816807022*tmp_7 + 0.091576213509770743*tmp_8 + 0.44594849091596489*tmp_9;
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
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = 1.0 / (tmp_1*tmp_3 - tmp_4*tmp_5);
      real_t tmp_7 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_8 = p_affine_6_1 + tmp_2;
      real_t tmp_9 = tmp_6*(0.21132486540518713*tmp_7 + tmp_8);
      real_t tmp_10 = tmp_1*tmp_9;
      real_t tmp_11 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_12 = tmp_11*tmp_9;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_0;
      real_t tmp_15 = tmp_6*(0.21132486540518713*tmp_13 + tmp_14);
      real_t tmp_16 = tmp_15*tmp_3;
      real_t tmp_17 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_18 = tmp_15*tmp_17;
      real_t tmp_19 = -0.5*tmp_10 - 0.5*tmp_12 - 0.5*tmp_16 - 0.5*tmp_18 + 0.5;
      real_t tmp_20 = tmp_12 + tmp_16;
      real_t tmp_21 = tmp_20 - 1.0/3.0;
      real_t tmp_22 = tmp_10 + tmp_18;
      real_t tmp_23 = tmp_22 - 1.0/3.0;
      real_t tmp_24 = p_affine_10_0*(tmp_1*tmp_21 + tmp_23*tmp_4);
      real_t tmp_25 = p_affine_10_1*(tmp_21*tmp_5 + tmp_23*tmp_3);
      real_t tmp_26 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_7*tmp_7), 1.0/2.0));
      real_t tmp_27 = 0.5*tmp_26;
      real_t tmp_28 = tmp_6*(0.78867513459481287*tmp_7 + tmp_8);
      real_t tmp_29 = tmp_1*tmp_28;
      real_t tmp_30 = tmp_11*tmp_28;
      real_t tmp_31 = tmp_6*(0.78867513459481287*tmp_13 + tmp_14);
      real_t tmp_32 = tmp_3*tmp_31;
      real_t tmp_33 = tmp_17*tmp_31;
      real_t tmp_34 = -0.5*tmp_29 - 0.5*tmp_30 - 0.5*tmp_32 - 0.5*tmp_33 + 0.5;
      real_t tmp_35 = tmp_30 + tmp_32;
      real_t tmp_36 = tmp_35 - 1.0/3.0;
      real_t tmp_37 = tmp_29 + tmp_33;
      real_t tmp_38 = tmp_37 - 1.0/3.0;
      real_t tmp_39 = p_affine_10_0*(tmp_1*tmp_36 + tmp_38*tmp_4);
      real_t tmp_40 = p_affine_10_1*(tmp_3*tmp_38 + tmp_36*tmp_5);
      real_t tmp_41 = 0.5*tmp_26;
      real_t tmp_42 = 0.5*tmp_20;
      real_t tmp_43 = 0.5*tmp_35;
      real_t tmp_44 = 0.5*tmp_22;
      real_t tmp_45 = 0.5*tmp_37;
      real_t a_0_0 = tmp_27*(tmp_19*tmp_24 + tmp_19*tmp_25) + tmp_41*(tmp_34*tmp_39 + tmp_34*tmp_40);
      real_t a_0_1 = tmp_27*(tmp_24*tmp_42 + tmp_25*tmp_42) + tmp_41*(tmp_39*tmp_43 + tmp_40*tmp_43);
      real_t a_0_2 = tmp_27*(tmp_24*tmp_44 + tmp_25*tmp_44) + tmp_41*(tmp_39*tmp_45 + tmp_40*tmp_45);
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
      real_t tmp_6 = p_affine_6_1 + 0.21132486540518713*tmp_5;
      real_t tmp_7 = tmp_4*(tmp_2 + tmp_6);
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_10 = tmp_7*tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.21132486540518713*tmp_11;
      real_t tmp_13 = tmp_4*(tmp_0 + tmp_12);
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_16 = tmp_13*tmp_15;
      real_t tmp_17 = -0.5*tmp_10 - 0.5*tmp_14 - 0.5*tmp_16 - 0.5*tmp_8 + 0.5;
      real_t tmp_18 = -p_affine_0_0;
      real_t tmp_19 = p_affine_1_0 + tmp_18;
      real_t tmp_20 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_21 = -p_affine_0_1;
      real_t tmp_22 = p_affine_2_1 + tmp_21;
      real_t tmp_23 = p_affine_2_0 + tmp_18;
      real_t tmp_24 = p_affine_1_1 + tmp_21;
      real_t tmp_25 = 1.0 / (tmp_19*tmp_22 - tmp_23*tmp_24);
      real_t tmp_26 = tmp_25*(tmp_21 + tmp_6);
      real_t tmp_27 = tmp_25*(tmp_12 + tmp_18);
      real_t tmp_28 = tmp_20*tmp_26 + tmp_22*tmp_27 - 1.0/3.0;
      real_t tmp_29 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_30 = tmp_19*tmp_26 + tmp_27*tmp_29 - 1.0/3.0;
      real_t tmp_31 = p_affine_10_0*(tmp_19*tmp_28 + tmp_23*tmp_30);
      real_t tmp_32 = p_affine_10_1*(tmp_22*tmp_30 + tmp_24*tmp_28);
      real_t tmp_33 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_5*tmp_5), 1.0/2.0));
      real_t tmp_34 = 0.5*tmp_33;
      real_t tmp_35 = p_affine_6_1 + 0.78867513459481287*tmp_5;
      real_t tmp_36 = tmp_4*(tmp_2 + tmp_35);
      real_t tmp_37 = tmp_1*tmp_36;
      real_t tmp_38 = tmp_36*tmp_9;
      real_t tmp_39 = p_affine_6_0 + 0.78867513459481287*tmp_11;
      real_t tmp_40 = tmp_4*(tmp_0 + tmp_39);
      real_t tmp_41 = tmp_3*tmp_40;
      real_t tmp_42 = tmp_15*tmp_40;
      real_t tmp_43 = -0.5*tmp_37 - 0.5*tmp_38 - 0.5*tmp_41 - 0.5*tmp_42 + 0.5;
      real_t tmp_44 = tmp_25*(tmp_21 + tmp_35);
      real_t tmp_45 = tmp_25*(tmp_18 + tmp_39);
      real_t tmp_46 = tmp_20*tmp_44 + tmp_22*tmp_45 - 1.0/3.0;
      real_t tmp_47 = tmp_19*tmp_44 + tmp_29*tmp_45 - 1.0/3.0;
      real_t tmp_48 = p_affine_10_0*(tmp_19*tmp_46 + tmp_23*tmp_47);
      real_t tmp_49 = p_affine_10_1*(tmp_22*tmp_47 + tmp_24*tmp_46);
      real_t tmp_50 = 0.5*tmp_33;
      real_t tmp_51 = 0.5*tmp_10 + 0.5*tmp_14;
      real_t tmp_52 = 0.5*tmp_38 + 0.5*tmp_41;
      real_t tmp_53 = 0.5*tmp_16 + 0.5*tmp_8;
      real_t tmp_54 = 0.5*tmp_37 + 0.5*tmp_42;
      real_t a_0_0 = tmp_34*(tmp_17*tmp_31 + tmp_17*tmp_32) + tmp_50*(tmp_43*tmp_48 + tmp_43*tmp_49);
      real_t a_0_1 = tmp_34*(tmp_31*tmp_51 + tmp_32*tmp_51) + tmp_50*(tmp_48*tmp_52 + tmp_49*tmp_52);
      real_t a_0_2 = tmp_34*(tmp_31*tmp_53 + tmp_32*tmp_53) + tmp_50*(tmp_48*tmp_54 + tmp_49*tmp_54);
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
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = 1.0 / (tmp_1*tmp_3 - tmp_4*tmp_5);
      real_t tmp_7 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_8 = p_affine_6_1 + tmp_2;
      real_t tmp_9 = tmp_6*(0.21132486540518713*tmp_7 + tmp_8);
      real_t tmp_10 = tmp_1*tmp_9;
      real_t tmp_11 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_12 = tmp_11*tmp_9;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_0;
      real_t tmp_15 = tmp_6*(0.21132486540518713*tmp_13 + tmp_14);
      real_t tmp_16 = tmp_15*tmp_3;
      real_t tmp_17 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_18 = tmp_15*tmp_17;
      real_t tmp_19 = -tmp_10 - tmp_12 - tmp_16 - tmp_18 + 1;
      real_t tmp_20 = tmp_12 + tmp_16;
      real_t tmp_21 = tmp_20 - 1.0/3.0;
      real_t tmp_22 = tmp_10 + tmp_18;
      real_t tmp_23 = tmp_22 - 1.0/3.0;
      real_t tmp_24 = p_affine_10_0*(tmp_1*tmp_21 + tmp_23*tmp_4);
      real_t tmp_25 = p_affine_10_1*(tmp_21*tmp_5 + tmp_23*tmp_3);
      real_t tmp_26 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_7*tmp_7), 1.0/2.0));
      real_t tmp_27 = 0.5*tmp_26;
      real_t tmp_28 = tmp_6*(0.78867513459481287*tmp_7 + tmp_8);
      real_t tmp_29 = tmp_1*tmp_28;
      real_t tmp_30 = tmp_11*tmp_28;
      real_t tmp_31 = tmp_6*(0.78867513459481287*tmp_13 + tmp_14);
      real_t tmp_32 = tmp_3*tmp_31;
      real_t tmp_33 = tmp_17*tmp_31;
      real_t tmp_34 = -tmp_29 - tmp_30 - tmp_32 - tmp_33 + 1;
      real_t tmp_35 = tmp_30 + tmp_32;
      real_t tmp_36 = tmp_35 - 1.0/3.0;
      real_t tmp_37 = tmp_29 + tmp_33;
      real_t tmp_38 = tmp_37 - 1.0/3.0;
      real_t tmp_39 = p_affine_10_0*(tmp_1*tmp_36 + tmp_38*tmp_4);
      real_t tmp_40 = p_affine_10_1*(tmp_3*tmp_38 + tmp_36*tmp_5);
      real_t tmp_41 = 0.5*tmp_26;
      real_t a_0_0 = tmp_27*(tmp_19*tmp_24 + tmp_19*tmp_25) + tmp_41*(tmp_34*tmp_39 + tmp_34*tmp_40);
      real_t a_0_1 = tmp_27*(tmp_20*tmp_24 + tmp_20*tmp_25) + tmp_41*(tmp_35*tmp_39 + tmp_35*tmp_40);
      real_t a_0_2 = tmp_27*(tmp_22*tmp_24 + tmp_22*tmp_25) + tmp_41*(tmp_37*tmp_39 + tmp_37*tmp_40);
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

      // Does nothing.
   }
};


} // dg
} // hyteg

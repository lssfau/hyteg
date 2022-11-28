
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

class EGVectorLaplaceFormP1P1_new_00 : public hyteg::dg::DGForm2D
{
 public:
    EGVectorLaplaceFormP1P1_new_00()
: callback_Scalar_Variable_Coefficient_2D_g0 ([](const Point3D & p) -> real_t { return 0.; })
    {}

void Scalar_Variable_Coefficient_2D_g0( real_t in_0, real_t in_1, real_t * out_0 ) const
{
   *out_0 = callback_Scalar_Variable_Coefficient_2D_g0( Point3D( {in_0, in_1, 0} ) );
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
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0);
      real_t tmp_5 = 1.0 / (tmp_4);
      real_t tmp_6 = tmp_1*tmp_5;
      real_t tmp_7 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_8 = tmp_5*tmp_7;
      real_t tmp_9 = -tmp_6 - tmp_8;
      real_t tmp_10 = tmp_3*tmp_5;
      real_t tmp_11 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = -tmp_10 - tmp_12;
      real_t tmp_14 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_15 = tmp_14*((tmp_13*tmp_13) + (tmp_9*tmp_9));
      real_t tmp_16 = tmp_14*(tmp_10*tmp_13 + tmp_8*tmp_9);
      real_t tmp_17 = 0.5*tmp_16;
      real_t tmp_18 = tmp_14*(tmp_12*tmp_13 + tmp_6*tmp_9);
      real_t tmp_19 = 0.5*tmp_18;
      real_t tmp_20 = 1.0 / (tmp_4*tmp_4);
      real_t tmp_21 = tmp_14*(tmp_20*(tmp_3*tmp_3) + tmp_20*(tmp_7*tmp_7));
      real_t tmp_22 = tmp_14*(tmp_1*tmp_20*tmp_7 + tmp_11*tmp_20*tmp_3);
      real_t tmp_23 = 0.5*tmp_22;
      real_t tmp_24 = tmp_14*((tmp_1*tmp_1)*tmp_20 + (tmp_11*tmp_11)*tmp_20);
      real_t a_0_0 = 0.5*tmp_15;
      real_t a_0_1 = tmp_17;
      real_t a_0_2 = tmp_19;
      real_t a_1_0 = tmp_17;
      real_t a_1_1 = 0.5*tmp_21;
      real_t a_1_2 = tmp_23;
      real_t a_2_0 = tmp_19;
      real_t a_2_1 = tmp_23;
      real_t a_2_2 = 0.5*tmp_24;
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

      real_t tmp_0 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_1 = -p_affine_0_1;
      real_t tmp_2 = p_affine_6_1 + tmp_1;
      real_t tmp_3 = 0.046910077030668018*tmp_0 + tmp_2;
      real_t tmp_4 = -p_affine_0_0;
      real_t tmp_5 = p_affine_1_0 + tmp_4;
      real_t tmp_6 = p_affine_2_1 + tmp_1;
      real_t tmp_7 = 1.0 / (tmp_5*tmp_6 - (p_affine_1_1 + tmp_1)*(p_affine_2_0 + tmp_4));
      real_t tmp_8 = tmp_5*tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = tmp_7*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_11 = tmp_10*tmp_3;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_4;
      real_t tmp_14 = 0.046910077030668018*tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6*tmp_7;
      real_t tmp_16 = tmp_14*tmp_15;
      real_t tmp_17 = tmp_7*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_18 = tmp_14*tmp_17;
      real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_12*tmp_12), 1.0/2.0));
      real_t tmp_21 = 6/tmp_20;
      real_t tmp_22 = p_affine_10_0*(-tmp_15 - tmp_17) + p_affine_10_1*(-tmp_10 - tmp_8);
      real_t tmp_23 = 1.0*tmp_22;
      real_t tmp_24 = 0.11846344252809471*tmp_20;
      real_t tmp_25 = 0.23076534494715845*tmp_0 + tmp_2;
      real_t tmp_26 = tmp_25*tmp_8;
      real_t tmp_27 = tmp_10*tmp_25;
      real_t tmp_28 = 0.23076534494715845*tmp_12 + tmp_13;
      real_t tmp_29 = tmp_15*tmp_28;
      real_t tmp_30 = tmp_17*tmp_28;
      real_t tmp_31 = -tmp_26 - tmp_27 - tmp_29 - tmp_30 + 1;
      real_t tmp_32 = 0.2393143352496831*tmp_20;
      real_t tmp_33 = 0.5*tmp_0 + tmp_2;
      real_t tmp_34 = tmp_33*tmp_8;
      real_t tmp_35 = tmp_10*tmp_33;
      real_t tmp_36 = 0.5*tmp_12 + tmp_13;
      real_t tmp_37 = tmp_15*tmp_36;
      real_t tmp_38 = tmp_17*tmp_36;
      real_t tmp_39 = -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1;
      real_t tmp_40 = 0.2844444444444445*tmp_20;
      real_t tmp_41 = 0.7692346550528415*tmp_0 + tmp_2;
      real_t tmp_42 = tmp_41*tmp_8;
      real_t tmp_43 = tmp_10*tmp_41;
      real_t tmp_44 = 0.7692346550528415*tmp_12 + tmp_13;
      real_t tmp_45 = tmp_15*tmp_44;
      real_t tmp_46 = tmp_17*tmp_44;
      real_t tmp_47 = -tmp_42 - tmp_43 - tmp_45 - tmp_46 + 1;
      real_t tmp_48 = 0.2393143352496831*tmp_20;
      real_t tmp_49 = 0.95308992296933193*tmp_0 + tmp_2;
      real_t tmp_50 = tmp_49*tmp_8;
      real_t tmp_51 = tmp_10*tmp_49;
      real_t tmp_52 = 0.95308992296933193*tmp_12 + tmp_13;
      real_t tmp_53 = tmp_15*tmp_52;
      real_t tmp_54 = tmp_17*tmp_52;
      real_t tmp_55 = -tmp_50 - tmp_51 - tmp_53 - tmp_54 + 1;
      real_t tmp_56 = 0.11846344252809471*tmp_20;
      real_t tmp_57 = tmp_11 + tmp_16;
      real_t tmp_58 = 0.5*tmp_22;
      real_t tmp_59 = p_affine_10_0*tmp_15 + p_affine_10_1*tmp_10;
      real_t tmp_60 = 0.5*tmp_59;
      real_t tmp_61 = tmp_19*tmp_21;
      real_t tmp_62 = tmp_27 + tmp_29;
      real_t tmp_63 = tmp_21*tmp_31;
      real_t tmp_64 = tmp_35 + tmp_37;
      real_t tmp_65 = tmp_21*tmp_39;
      real_t tmp_66 = tmp_43 + tmp_45;
      real_t tmp_67 = tmp_21*tmp_47;
      real_t tmp_68 = tmp_51 + tmp_53;
      real_t tmp_69 = tmp_21*tmp_55;
      real_t tmp_70 = tmp_24*(-tmp_19*tmp_60 - tmp_57*tmp_58 + tmp_57*tmp_61) + tmp_32*(-tmp_31*tmp_60 - tmp_58*tmp_62 + tmp_62*tmp_63) + tmp_40*(-tmp_39*tmp_60 - tmp_58*tmp_64 + tmp_64*tmp_65) + tmp_48*(-tmp_47*tmp_60 - tmp_58*tmp_66 + tmp_66*tmp_67) + tmp_56*(-tmp_55*tmp_60 - tmp_58*tmp_68 + tmp_68*tmp_69);
      real_t tmp_71 = tmp_18 + tmp_9;
      real_t tmp_72 = p_affine_10_0*tmp_17 + p_affine_10_1*tmp_8;
      real_t tmp_73 = 0.5*tmp_72;
      real_t tmp_74 = tmp_26 + tmp_30;
      real_t tmp_75 = tmp_34 + tmp_38;
      real_t tmp_76 = tmp_42 + tmp_46;
      real_t tmp_77 = tmp_50 + tmp_54;
      real_t tmp_78 = tmp_24*(-tmp_19*tmp_73 - tmp_58*tmp_71 + tmp_61*tmp_71) + tmp_32*(-tmp_31*tmp_73 - tmp_58*tmp_74 + tmp_63*tmp_74) + tmp_40*(-tmp_39*tmp_73 - tmp_58*tmp_75 + tmp_65*tmp_75) + tmp_48*(-tmp_47*tmp_73 - tmp_58*tmp_76 + tmp_67*tmp_76) + tmp_56*(-tmp_55*tmp_73 - tmp_58*tmp_77 + tmp_69*tmp_77);
      real_t tmp_79 = 1.0*tmp_59;
      real_t tmp_80 = tmp_24*(tmp_21*tmp_57*tmp_71 - tmp_57*tmp_73 - tmp_60*tmp_71) + tmp_32*(tmp_21*tmp_62*tmp_74 - tmp_60*tmp_74 - tmp_62*tmp_73) + tmp_40*(tmp_21*tmp_64*tmp_75 - tmp_60*tmp_75 - tmp_64*tmp_73) + tmp_48*(tmp_21*tmp_66*tmp_76 - tmp_60*tmp_76 - tmp_66*tmp_73) + tmp_56*(tmp_21*tmp_68*tmp_77 - tmp_60*tmp_77 - tmp_68*tmp_73);
      real_t tmp_81 = 1.0*tmp_72;
      real_t a_0_0 = tmp_24*((tmp_19*tmp_19)*tmp_21 - tmp_19*tmp_23) + tmp_32*(tmp_21*(tmp_31*tmp_31) - tmp_23*tmp_31) + tmp_40*(tmp_21*(tmp_39*tmp_39) - tmp_23*tmp_39) + tmp_48*(tmp_21*(tmp_47*tmp_47) - tmp_23*tmp_47) + tmp_56*(tmp_21*(tmp_55*tmp_55) - tmp_23*tmp_55);
      real_t a_0_1 = tmp_70;
      real_t a_0_2 = tmp_78;
      real_t a_1_0 = tmp_70;
      real_t a_1_1 = tmp_24*(tmp_21*(tmp_57*tmp_57) - tmp_57*tmp_79) + tmp_32*(tmp_21*(tmp_62*tmp_62) - tmp_62*tmp_79) + tmp_40*(tmp_21*(tmp_64*tmp_64) - tmp_64*tmp_79) + tmp_48*(tmp_21*(tmp_66*tmp_66) - tmp_66*tmp_79) + tmp_56*(tmp_21*(tmp_68*tmp_68) - tmp_68*tmp_79);
      real_t a_1_2 = tmp_80;
      real_t a_2_0 = tmp_78;
      real_t a_2_1 = tmp_80;
      real_t a_2_2 = tmp_24*(tmp_21*(tmp_71*tmp_71) - tmp_71*tmp_81) + tmp_32*(tmp_21*(tmp_74*tmp_74) - tmp_74*tmp_81) + tmp_40*(tmp_21*(tmp_75*tmp_75) - tmp_75*tmp_81) + tmp_48*(tmp_21*(tmp_76*tmp_76) - tmp_76*tmp_81) + tmp_56*(tmp_21*(tmp_77*tmp_77) - tmp_77*tmp_81);
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

      real_t tmp_0 = -p_affine_3_1;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = p_affine_6_1 + 0.046910077030668018*tmp_1;
      real_t tmp_3 = tmp_0 + tmp_2;
      real_t tmp_4 = -p_affine_3_0;
      real_t tmp_5 = p_affine_4_0 + tmp_4;
      real_t tmp_6 = p_affine_5_1 + tmp_0;
      real_t tmp_7 = 1.0 / (tmp_5*tmp_6 - (p_affine_4_1 + tmp_0)*(p_affine_5_0 + tmp_4));
      real_t tmp_8 = tmp_5*tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = tmp_7*(p_affine_3_0 - p_affine_5_0);
      real_t tmp_11 = tmp_10*tmp_3;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + 0.046910077030668018*tmp_12;
      real_t tmp_14 = tmp_13 + tmp_4;
      real_t tmp_15 = tmp_6*tmp_7;
      real_t tmp_16 = tmp_14*tmp_15;
      real_t tmp_17 = tmp_7*(p_affine_3_1 - p_affine_4_1);
      real_t tmp_18 = tmp_14*tmp_17;
      real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20 = -p_affine_0_1;
      real_t tmp_21 = p_affine_2_1 + tmp_20;
      real_t tmp_22 = -p_affine_0_0;
      real_t tmp_23 = p_affine_1_0 + tmp_22;
      real_t tmp_24 = 1.0 / (tmp_21*tmp_23 - (p_affine_1_1 + tmp_20)*(p_affine_2_0 + tmp_22));
      real_t tmp_25 = tmp_21*tmp_24;
      real_t tmp_26 = tmp_24*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_27 = tmp_23*tmp_24;
      real_t tmp_28 = tmp_24*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_29 = 0.5*p_affine_10_0*(-tmp_25 - tmp_26) + 0.5*p_affine_10_1*(-tmp_27 - tmp_28);
      real_t tmp_30 = tmp_2 + tmp_20;
      real_t tmp_31 = tmp_27*tmp_30;
      real_t tmp_32 = tmp_28*tmp_30;
      real_t tmp_33 = tmp_13 + tmp_22;
      real_t tmp_34 = tmp_25*tmp_33;
      real_t tmp_35 = tmp_26*tmp_33;
      real_t tmp_36 = -tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1;
      real_t tmp_37 = 0.5*p_affine_10_0*(-tmp_15 - tmp_17) + 0.5*p_affine_10_1*(-tmp_10 - tmp_8);
      real_t tmp_38 = std::abs(std::pow((tmp_1*tmp_1) + (tmp_12*tmp_12), 1.0/2.0));
      real_t tmp_39 = 6/tmp_38;
      real_t tmp_40 = tmp_36*tmp_39;
      real_t tmp_41 = 0.11846344252809471*tmp_38;
      real_t tmp_42 = p_affine_6_1 + 0.23076534494715845*tmp_1;
      real_t tmp_43 = tmp_0 + tmp_42;
      real_t tmp_44 = tmp_43*tmp_8;
      real_t tmp_45 = tmp_10*tmp_43;
      real_t tmp_46 = p_affine_6_0 + 0.23076534494715845*tmp_12;
      real_t tmp_47 = tmp_4 + tmp_46;
      real_t tmp_48 = tmp_15*tmp_47;
      real_t tmp_49 = tmp_17*tmp_47;
      real_t tmp_50 = -tmp_44 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_51 = tmp_20 + tmp_42;
      real_t tmp_52 = tmp_27*tmp_51;
      real_t tmp_53 = tmp_28*tmp_51;
      real_t tmp_54 = tmp_22 + tmp_46;
      real_t tmp_55 = tmp_25*tmp_54;
      real_t tmp_56 = tmp_26*tmp_54;
      real_t tmp_57 = -tmp_52 - tmp_53 - tmp_55 - tmp_56 + 1;
      real_t tmp_58 = tmp_39*tmp_57;
      real_t tmp_59 = 0.2393143352496831*tmp_38;
      real_t tmp_60 = p_affine_6_1 + 0.5*tmp_1;
      real_t tmp_61 = tmp_0 + tmp_60;
      real_t tmp_62 = tmp_61*tmp_8;
      real_t tmp_63 = tmp_10*tmp_61;
      real_t tmp_64 = p_affine_6_0 + 0.5*tmp_12;
      real_t tmp_65 = tmp_4 + tmp_64;
      real_t tmp_66 = tmp_15*tmp_65;
      real_t tmp_67 = tmp_17*tmp_65;
      real_t tmp_68 = -tmp_62 - tmp_63 - tmp_66 - tmp_67 + 1;
      real_t tmp_69 = tmp_20 + tmp_60;
      real_t tmp_70 = tmp_27*tmp_69;
      real_t tmp_71 = tmp_28*tmp_69;
      real_t tmp_72 = tmp_22 + tmp_64;
      real_t tmp_73 = tmp_25*tmp_72;
      real_t tmp_74 = tmp_26*tmp_72;
      real_t tmp_75 = -tmp_70 - tmp_71 - tmp_73 - tmp_74 + 1;
      real_t tmp_76 = tmp_39*tmp_75;
      real_t tmp_77 = 0.2844444444444445*tmp_38;
      real_t tmp_78 = p_affine_6_1 + 0.7692346550528415*tmp_1;
      real_t tmp_79 = tmp_0 + tmp_78;
      real_t tmp_80 = tmp_79*tmp_8;
      real_t tmp_81 = tmp_10*tmp_79;
      real_t tmp_82 = p_affine_6_0 + 0.7692346550528415*tmp_12;
      real_t tmp_83 = tmp_4 + tmp_82;
      real_t tmp_84 = tmp_15*tmp_83;
      real_t tmp_85 = tmp_17*tmp_83;
      real_t tmp_86 = -tmp_80 - tmp_81 - tmp_84 - tmp_85 + 1;
      real_t tmp_87 = tmp_20 + tmp_78;
      real_t tmp_88 = tmp_27*tmp_87;
      real_t tmp_89 = tmp_28*tmp_87;
      real_t tmp_90 = tmp_22 + tmp_82;
      real_t tmp_91 = tmp_25*tmp_90;
      real_t tmp_92 = tmp_26*tmp_90;
      real_t tmp_93 = -tmp_88 - tmp_89 - tmp_91 - tmp_92 + 1;
      real_t tmp_94 = tmp_39*tmp_93;
      real_t tmp_95 = 0.2393143352496831*tmp_38;
      real_t tmp_96 = p_affine_6_1 + 0.95308992296933193*tmp_1;
      real_t tmp_97 = tmp_0 + tmp_96;
      real_t tmp_98 = tmp_8*tmp_97;
      real_t tmp_99 = tmp_10*tmp_97;
      real_t tmp_100 = p_affine_6_0 + 0.95308992296933193*tmp_12;
      real_t tmp_101 = tmp_100 + tmp_4;
      real_t tmp_102 = tmp_101*tmp_15;
      real_t tmp_103 = tmp_101*tmp_17;
      real_t tmp_104 = -tmp_102 - tmp_103 - tmp_98 - tmp_99 + 1;
      real_t tmp_105 = tmp_20 + tmp_96;
      real_t tmp_106 = tmp_105*tmp_27;
      real_t tmp_107 = tmp_105*tmp_28;
      real_t tmp_108 = tmp_100 + tmp_22;
      real_t tmp_109 = tmp_108*tmp_25;
      real_t tmp_110 = tmp_108*tmp_26;
      real_t tmp_111 = -tmp_106 - tmp_107 - tmp_109 - tmp_110 + 1;
      real_t tmp_112 = tmp_111*tmp_39;
      real_t tmp_113 = 0.11846344252809471*tmp_38;
      real_t tmp_114 = tmp_11 + tmp_16;
      real_t tmp_115 = 0.5*p_affine_10_0*tmp_15 + 0.5*p_affine_10_1*tmp_10;
      real_t tmp_116 = tmp_45 + tmp_48;
      real_t tmp_117 = tmp_63 + tmp_66;
      real_t tmp_118 = tmp_81 + tmp_84;
      real_t tmp_119 = tmp_102 + tmp_99;
      real_t tmp_120 = tmp_18 + tmp_9;
      real_t tmp_121 = 0.5*p_affine_10_0*tmp_17 + 0.5*p_affine_10_1*tmp_8;
      real_t tmp_122 = tmp_44 + tmp_49;
      real_t tmp_123 = tmp_62 + tmp_67;
      real_t tmp_124 = tmp_80 + tmp_85;
      real_t tmp_125 = tmp_103 + tmp_98;
      real_t tmp_126 = tmp_32 + tmp_34;
      real_t tmp_127 = 0.5*p_affine_10_0*tmp_25 + 0.5*p_affine_10_1*tmp_28;
      real_t tmp_128 = tmp_126*tmp_39;
      real_t tmp_129 = tmp_53 + tmp_55;
      real_t tmp_130 = tmp_129*tmp_39;
      real_t tmp_131 = tmp_71 + tmp_73;
      real_t tmp_132 = tmp_131*tmp_39;
      real_t tmp_133 = tmp_89 + tmp_91;
      real_t tmp_134 = tmp_133*tmp_39;
      real_t tmp_135 = tmp_107 + tmp_109;
      real_t tmp_136 = tmp_135*tmp_39;
      real_t tmp_137 = tmp_31 + tmp_35;
      real_t tmp_138 = 0.5*p_affine_10_0*tmp_26 + 0.5*p_affine_10_1*tmp_27;
      real_t tmp_139 = tmp_137*tmp_39;
      real_t tmp_140 = tmp_52 + tmp_56;
      real_t tmp_141 = tmp_140*tmp_39;
      real_t tmp_142 = tmp_70 + tmp_74;
      real_t tmp_143 = tmp_142*tmp_39;
      real_t tmp_144 = tmp_88 + tmp_92;
      real_t tmp_145 = tmp_144*tmp_39;
      real_t tmp_146 = tmp_106 + tmp_110;
      real_t tmp_147 = tmp_146*tmp_39;
      real_t a_0_0 = tmp_113*(-tmp_104*tmp_112 + tmp_104*tmp_29 - tmp_111*tmp_37) + tmp_41*(tmp_19*tmp_29 - tmp_19*tmp_40 - tmp_36*tmp_37) + tmp_59*(tmp_29*tmp_50 - tmp_37*tmp_57 - tmp_50*tmp_58) + tmp_77*(tmp_29*tmp_68 - tmp_37*tmp_75 - tmp_68*tmp_76) + tmp_95*(tmp_29*tmp_86 - tmp_37*tmp_93 - tmp_86*tmp_94);
      real_t a_0_1 = tmp_113*(-tmp_111*tmp_115 - tmp_112*tmp_119 + tmp_119*tmp_29) + tmp_41*(tmp_114*tmp_29 - tmp_114*tmp_40 - tmp_115*tmp_36) + tmp_59*(-tmp_115*tmp_57 + tmp_116*tmp_29 - tmp_116*tmp_58) + tmp_77*(-tmp_115*tmp_75 + tmp_117*tmp_29 - tmp_117*tmp_76) + tmp_95*(-tmp_115*tmp_93 + tmp_118*tmp_29 - tmp_118*tmp_94);
      real_t a_0_2 = tmp_113*(-tmp_111*tmp_121 - tmp_112*tmp_125 + tmp_125*tmp_29) + tmp_41*(tmp_120*tmp_29 - tmp_120*tmp_40 - tmp_121*tmp_36) + tmp_59*(-tmp_121*tmp_57 + tmp_122*tmp_29 - tmp_122*tmp_58) + tmp_77*(-tmp_121*tmp_75 + tmp_123*tmp_29 - tmp_123*tmp_76) + tmp_95*(-tmp_121*tmp_93 + tmp_124*tmp_29 - tmp_124*tmp_94);
      real_t a_1_0 = tmp_113*(tmp_104*tmp_127 - tmp_104*tmp_136 - tmp_135*tmp_37) + tmp_41*(-tmp_126*tmp_37 + tmp_127*tmp_19 - tmp_128*tmp_19) + tmp_59*(tmp_127*tmp_50 - tmp_129*tmp_37 - tmp_130*tmp_50) + tmp_77*(tmp_127*tmp_68 - tmp_131*tmp_37 - tmp_132*tmp_68) + tmp_95*(tmp_127*tmp_86 - tmp_133*tmp_37 - tmp_134*tmp_86);
      real_t a_1_1 = tmp_113*(-tmp_115*tmp_135 + tmp_119*tmp_127 - tmp_119*tmp_136) + tmp_41*(tmp_114*tmp_127 - tmp_114*tmp_128 - tmp_115*tmp_126) + tmp_59*(-tmp_115*tmp_129 + tmp_116*tmp_127 - tmp_116*tmp_130) + tmp_77*(-tmp_115*tmp_131 + tmp_117*tmp_127 - tmp_117*tmp_132) + tmp_95*(-tmp_115*tmp_133 + tmp_118*tmp_127 - tmp_118*tmp_134);
      real_t a_1_2 = tmp_113*(-tmp_121*tmp_135 + tmp_125*tmp_127 - tmp_125*tmp_136) + tmp_41*(tmp_120*tmp_127 - tmp_120*tmp_128 - tmp_121*tmp_126) + tmp_59*(-tmp_121*tmp_129 + tmp_122*tmp_127 - tmp_122*tmp_130) + tmp_77*(-tmp_121*tmp_131 + tmp_123*tmp_127 - tmp_123*tmp_132) + tmp_95*(-tmp_121*tmp_133 + tmp_124*tmp_127 - tmp_124*tmp_134);
      real_t a_2_0 = tmp_113*(tmp_104*tmp_138 - tmp_104*tmp_147 - tmp_146*tmp_37) + tmp_41*(-tmp_137*tmp_37 + tmp_138*tmp_19 - tmp_139*tmp_19) + tmp_59*(tmp_138*tmp_50 - tmp_140*tmp_37 - tmp_141*tmp_50) + tmp_77*(tmp_138*tmp_68 - tmp_142*tmp_37 - tmp_143*tmp_68) + tmp_95*(tmp_138*tmp_86 - tmp_144*tmp_37 - tmp_145*tmp_86);
      real_t a_2_1 = tmp_113*(-tmp_115*tmp_146 + tmp_119*tmp_138 - tmp_119*tmp_147) + tmp_41*(tmp_114*tmp_138 - tmp_114*tmp_139 - tmp_115*tmp_137) + tmp_59*(-tmp_115*tmp_140 + tmp_116*tmp_138 - tmp_116*tmp_141) + tmp_77*(-tmp_115*tmp_142 + tmp_117*tmp_138 - tmp_117*tmp_143) + tmp_95*(-tmp_115*tmp_144 + tmp_118*tmp_138 - tmp_118*tmp_145);
      real_t a_2_2 = tmp_113*(-tmp_121*tmp_146 + tmp_125*tmp_138 - tmp_125*tmp_147) + tmp_41*(tmp_120*tmp_138 - tmp_120*tmp_139 - tmp_121*tmp_137) + tmp_59*(-tmp_121*tmp_140 + tmp_122*tmp_138 - tmp_122*tmp_141) + tmp_77*(-tmp_121*tmp_142 + tmp_123*tmp_138 - tmp_123*tmp_143) + tmp_95*(-tmp_121*tmp_144 + tmp_124*tmp_138 - tmp_124*tmp_145);
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

      real_t tmp_0 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_1 = -p_affine_0_1;
      real_t tmp_2 = p_affine_6_1 + tmp_1;
      real_t tmp_3 = 0.046910077030668018*tmp_0 + tmp_2;
      real_t tmp_4 = -p_affine_0_0;
      real_t tmp_5 = p_affine_1_0 + tmp_4;
      real_t tmp_6 = p_affine_2_1 + tmp_1;
      real_t tmp_7 = 1.0 / (tmp_5*tmp_6 - (p_affine_1_1 + tmp_1)*(p_affine_2_0 + tmp_4));
      real_t tmp_8 = tmp_5*tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = tmp_7*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_11 = tmp_10*tmp_3;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_4;
      real_t tmp_14 = 0.046910077030668018*tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6*tmp_7;
      real_t tmp_16 = tmp_14*tmp_15;
      real_t tmp_17 = tmp_7*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_18 = tmp_14*tmp_17;
      real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_12*tmp_12), 1.0/2.0));
      real_t tmp_21 = 6/tmp_20;
      real_t tmp_22 = p_affine_10_0*(-tmp_15 - tmp_17) + p_affine_10_1*(-tmp_10 - tmp_8);
      real_t tmp_23 = 2*tmp_22;
      real_t tmp_24 = 0.11846344252809471*tmp_20;
      real_t tmp_25 = 0.23076534494715845*tmp_0 + tmp_2;
      real_t tmp_26 = tmp_25*tmp_8;
      real_t tmp_27 = tmp_10*tmp_25;
      real_t tmp_28 = 0.23076534494715845*tmp_12 + tmp_13;
      real_t tmp_29 = tmp_15*tmp_28;
      real_t tmp_30 = tmp_17*tmp_28;
      real_t tmp_31 = -tmp_26 - tmp_27 - tmp_29 - tmp_30 + 1;
      real_t tmp_32 = 0.2393143352496831*tmp_20;
      real_t tmp_33 = 0.5*tmp_0 + tmp_2;
      real_t tmp_34 = tmp_33*tmp_8;
      real_t tmp_35 = tmp_10*tmp_33;
      real_t tmp_36 = 0.5*tmp_12 + tmp_13;
      real_t tmp_37 = tmp_15*tmp_36;
      real_t tmp_38 = tmp_17*tmp_36;
      real_t tmp_39 = -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1;
      real_t tmp_40 = 0.2844444444444445*tmp_20;
      real_t tmp_41 = 0.7692346550528415*tmp_0 + tmp_2;
      real_t tmp_42 = tmp_41*tmp_8;
      real_t tmp_43 = tmp_10*tmp_41;
      real_t tmp_44 = 0.7692346550528415*tmp_12 + tmp_13;
      real_t tmp_45 = tmp_15*tmp_44;
      real_t tmp_46 = tmp_17*tmp_44;
      real_t tmp_47 = -tmp_42 - tmp_43 - tmp_45 - tmp_46 + 1;
      real_t tmp_48 = 0.2393143352496831*tmp_20;
      real_t tmp_49 = 0.95308992296933193*tmp_0 + tmp_2;
      real_t tmp_50 = tmp_49*tmp_8;
      real_t tmp_51 = tmp_10*tmp_49;
      real_t tmp_52 = 0.95308992296933193*tmp_12 + tmp_13;
      real_t tmp_53 = tmp_15*tmp_52;
      real_t tmp_54 = tmp_17*tmp_52;
      real_t tmp_55 = -tmp_50 - tmp_51 - tmp_53 - tmp_54 + 1;
      real_t tmp_56 = 0.11846344252809471*tmp_20;
      real_t tmp_57 = tmp_11 + tmp_16;
      real_t tmp_58 = p_affine_10_0*tmp_15 + p_affine_10_1*tmp_10;
      real_t tmp_59 = tmp_19*tmp_21;
      real_t tmp_60 = tmp_27 + tmp_29;
      real_t tmp_61 = tmp_21*tmp_31;
      real_t tmp_62 = tmp_35 + tmp_37;
      real_t tmp_63 = tmp_21*tmp_39;
      real_t tmp_64 = tmp_43 + tmp_45;
      real_t tmp_65 = tmp_21*tmp_47;
      real_t tmp_66 = tmp_51 + tmp_53;
      real_t tmp_67 = tmp_21*tmp_55;
      real_t tmp_68 = tmp_24*(-tmp_19*tmp_58 - tmp_22*tmp_57 + tmp_57*tmp_59) + tmp_32*(-tmp_22*tmp_60 - tmp_31*tmp_58 + tmp_60*tmp_61) + tmp_40*(-tmp_22*tmp_62 - tmp_39*tmp_58 + tmp_62*tmp_63) + tmp_48*(-tmp_22*tmp_64 - tmp_47*tmp_58 + tmp_64*tmp_65) + tmp_56*(-tmp_22*tmp_66 - tmp_55*tmp_58 + tmp_66*tmp_67);
      real_t tmp_69 = tmp_18 + tmp_9;
      real_t tmp_70 = p_affine_10_0*tmp_17 + p_affine_10_1*tmp_8;
      real_t tmp_71 = tmp_26 + tmp_30;
      real_t tmp_72 = tmp_34 + tmp_38;
      real_t tmp_73 = tmp_42 + tmp_46;
      real_t tmp_74 = tmp_50 + tmp_54;
      real_t tmp_75 = tmp_24*(-tmp_19*tmp_70 - tmp_22*tmp_69 + tmp_59*tmp_69) + tmp_32*(-tmp_22*tmp_71 - tmp_31*tmp_70 + tmp_61*tmp_71) + tmp_40*(-tmp_22*tmp_72 - tmp_39*tmp_70 + tmp_63*tmp_72) + tmp_48*(-tmp_22*tmp_73 - tmp_47*tmp_70 + tmp_65*tmp_73) + tmp_56*(-tmp_22*tmp_74 - tmp_55*tmp_70 + tmp_67*tmp_74);
      real_t tmp_76 = 2*tmp_58;
      real_t tmp_77 = tmp_24*(tmp_21*tmp_57*tmp_69 - tmp_57*tmp_70 - tmp_58*tmp_69) + tmp_32*(tmp_21*tmp_60*tmp_71 - tmp_58*tmp_71 - tmp_60*tmp_70) + tmp_40*(tmp_21*tmp_62*tmp_72 - tmp_58*tmp_72 - tmp_62*tmp_70) + tmp_48*(tmp_21*tmp_64*tmp_73 - tmp_58*tmp_73 - tmp_64*tmp_70) + tmp_56*(tmp_21*tmp_66*tmp_74 - tmp_58*tmp_74 - tmp_66*tmp_70);
      real_t tmp_78 = 2*tmp_70;
      real_t a_0_0 = tmp_24*((tmp_19*tmp_19)*tmp_21 - tmp_19*tmp_23) + tmp_32*(tmp_21*(tmp_31*tmp_31) - tmp_23*tmp_31) + tmp_40*(tmp_21*(tmp_39*tmp_39) - tmp_23*tmp_39) + tmp_48*(tmp_21*(tmp_47*tmp_47) - tmp_23*tmp_47) + tmp_56*(tmp_21*(tmp_55*tmp_55) - tmp_23*tmp_55);
      real_t a_0_1 = tmp_68;
      real_t a_0_2 = tmp_75;
      real_t a_1_0 = tmp_68;
      real_t a_1_1 = tmp_24*(tmp_21*(tmp_57*tmp_57) - tmp_57*tmp_76) + tmp_32*(tmp_21*(tmp_60*tmp_60) - tmp_60*tmp_76) + tmp_40*(tmp_21*(tmp_62*tmp_62) - tmp_62*tmp_76) + tmp_48*(tmp_21*(tmp_64*tmp_64) - tmp_64*tmp_76) + tmp_56*(tmp_21*(tmp_66*tmp_66) - tmp_66*tmp_76);
      real_t a_1_2 = tmp_77;
      real_t a_2_0 = tmp_75;
      real_t a_2_1 = tmp_77;
      real_t a_2_2 = tmp_24*(tmp_21*(tmp_69*tmp_69) - tmp_69*tmp_78) + tmp_32*(tmp_21*(tmp_71*tmp_71) - tmp_71*tmp_78) + tmp_40*(tmp_21*(tmp_72*tmp_72) - tmp_72*tmp_78) + tmp_48*(tmp_21*(tmp_73*tmp_73) - tmp_73*tmp_78) + tmp_56*(tmp_21*(tmp_74*tmp_74) - tmp_74*tmp_78);
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
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_g0_out0_id4 = 0;
      Scalar_Variable_Coefficient_2D_g0( 0.95308992296933193*p_affine_6_0 + 0.046910077030668018*p_affine_7_0, 0.95308992296933193*p_affine_6_1 + 0.046910077030668018*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id0 );
      Scalar_Variable_Coefficient_2D_g0( 0.7692346550528415*p_affine_6_0 + 0.23076534494715845*p_affine_7_0, 0.7692346550528415*p_affine_6_1 + 0.23076534494715845*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id1 );
      Scalar_Variable_Coefficient_2D_g0( 0.5*p_affine_6_0 + 0.5*p_affine_7_0, 0.5*p_affine_6_1 + 0.5*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id2 );
      Scalar_Variable_Coefficient_2D_g0( 0.2307653449471585*p_affine_6_0 + 0.7692346550528415*p_affine_7_0, 0.2307653449471585*p_affine_6_1 + 0.7692346550528415*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id3 );
      Scalar_Variable_Coefficient_2D_g0( 0.046910077030668074*p_affine_6_0 + 0.95308992296933193*p_affine_7_0, 0.046910077030668074*p_affine_6_1 + 0.95308992296933193*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g0_out0_id4 );
      real_t tmp_0 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_1 = -p_affine_0_1;
      real_t tmp_2 = p_affine_6_1 + tmp_1;
      real_t tmp_3 = 0.046910077030668018*tmp_0 + tmp_2;
      real_t tmp_4 = -p_affine_0_0;
      real_t tmp_5 = p_affine_1_0 + tmp_4;
      real_t tmp_6 = p_affine_2_1 + tmp_1;
      real_t tmp_7 = 1.0 / (tmp_5*tmp_6 - (p_affine_1_1 + tmp_1)*(p_affine_2_0 + tmp_4));
      real_t tmp_8 = tmp_5*tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = tmp_7*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_11 = tmp_10*tmp_3;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_4;
      real_t tmp_14 = 0.046910077030668018*tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6*tmp_7;
      real_t tmp_16 = tmp_14*tmp_15;
      real_t tmp_17 = tmp_7*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_18 = tmp_14*tmp_17;
      real_t tmp_19 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_12*tmp_12), 1.0/2.0));
      real_t tmp_20 = 6/tmp_19;
      real_t tmp_21 = -p_affine_10_0*(-tmp_15 - tmp_17) - p_affine_10_1*(-tmp_10 - tmp_8);
      real_t tmp_22 = 0.11846344252809471*Scalar_Variable_Coefficient_2D_g0_out0_id0*tmp_19;
      real_t tmp_23 = 0.23076534494715845*tmp_0 + tmp_2;
      real_t tmp_24 = tmp_23*tmp_8;
      real_t tmp_25 = tmp_10*tmp_23;
      real_t tmp_26 = 0.23076534494715845*tmp_12 + tmp_13;
      real_t tmp_27 = tmp_15*tmp_26;
      real_t tmp_28 = tmp_17*tmp_26;
      real_t tmp_29 = 0.2393143352496831*Scalar_Variable_Coefficient_2D_g0_out0_id1*tmp_19;
      real_t tmp_30 = 0.5*tmp_0 + tmp_2;
      real_t tmp_31 = tmp_30*tmp_8;
      real_t tmp_32 = tmp_10*tmp_30;
      real_t tmp_33 = 0.5*tmp_12 + tmp_13;
      real_t tmp_34 = tmp_15*tmp_33;
      real_t tmp_35 = tmp_17*tmp_33;
      real_t tmp_36 = 0.2844444444444445*Scalar_Variable_Coefficient_2D_g0_out0_id2*tmp_19;
      real_t tmp_37 = 0.7692346550528415*tmp_0 + tmp_2;
      real_t tmp_38 = tmp_37*tmp_8;
      real_t tmp_39 = tmp_10*tmp_37;
      real_t tmp_40 = 0.7692346550528415*tmp_12 + tmp_13;
      real_t tmp_41 = tmp_15*tmp_40;
      real_t tmp_42 = tmp_17*tmp_40;
      real_t tmp_43 = 0.2393143352496831*Scalar_Variable_Coefficient_2D_g0_out0_id3*tmp_19;
      real_t tmp_44 = 0.95308992296933193*tmp_0 + tmp_2;
      real_t tmp_45 = tmp_44*tmp_8;
      real_t tmp_46 = tmp_10*tmp_44;
      real_t tmp_47 = 0.95308992296933193*tmp_12 + tmp_13;
      real_t tmp_48 = tmp_15*tmp_47;
      real_t tmp_49 = tmp_17*tmp_47;
      real_t tmp_50 = 0.11846344252809471*Scalar_Variable_Coefficient_2D_g0_out0_id4*tmp_19;
      real_t tmp_51 = -p_affine_10_0*tmp_15 - p_affine_10_1*tmp_10;
      real_t tmp_52 = -p_affine_10_0*tmp_17 - p_affine_10_1*tmp_8;
      real_t a_0_0 = tmp_22*(tmp_20*(-tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1) + tmp_21) + tmp_29*(tmp_20*(-tmp_24 - tmp_25 - tmp_27 - tmp_28 + 1) + tmp_21) + tmp_36*(tmp_20*(-tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1) + tmp_21) + tmp_43*(tmp_20*(-tmp_38 - tmp_39 - tmp_41 - tmp_42 + 1) + tmp_21) + tmp_50*(tmp_20*(-tmp_45 - tmp_46 - tmp_48 - tmp_49 + 1) + tmp_21);
      real_t a_1_0 = tmp_22*(tmp_20*(tmp_11 + tmp_16) + tmp_51) + tmp_29*(tmp_20*(tmp_25 + tmp_27) + tmp_51) + tmp_36*(tmp_20*(tmp_32 + tmp_34) + tmp_51) + tmp_43*(tmp_20*(tmp_39 + tmp_41) + tmp_51) + tmp_50*(tmp_20*(tmp_46 + tmp_48) + tmp_51);
      real_t a_2_0 = tmp_22*(tmp_20*(tmp_18 + tmp_9) + tmp_52) + tmp_29*(tmp_20*(tmp_24 + tmp_28) + tmp_52) + tmp_36*(tmp_20*(tmp_31 + tmp_35) + tmp_52) + tmp_43*(tmp_20*(tmp_38 + tmp_42) + tmp_52) + tmp_50*(tmp_20*(tmp_45 + tmp_49) + tmp_52);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   }

public:

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g0;
};




class EGVectorLaplaceFormP1P1_new_10 : public hyteg::dg::DGForm2D
{
 public:
    EGVectorLaplaceFormP1P1_new_10()

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

      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      real_t a_1_0 = 0;
      real_t a_1_1 = 0;
      real_t a_1_2 = 0;
      real_t a_2_0 = 0;
      real_t a_2_1 = 0;
      real_t a_2_2 = 0;
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

      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      real_t a_1_0 = 0;
      real_t a_1_1 = 0;
      real_t a_1_2 = 0;
      real_t a_2_0 = 0;
      real_t a_2_1 = 0;
      real_t a_2_2 = 0;
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

      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      real_t a_1_0 = 0;
      real_t a_1_1 = 0;
      real_t a_1_2 = 0;
      real_t a_2_0 = 0;
      real_t a_2_1 = 0;
      real_t a_2_2 = 0;
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

      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      real_t a_1_0 = 0;
      real_t a_1_1 = 0;
      real_t a_1_2 = 0;
      real_t a_2_0 = 0;
      real_t a_2_1 = 0;
      real_t a_2_2 = 0;
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




class EGVectorLaplaceFormP1P1_new_01 : public hyteg::dg::DGForm2D
{
 public:
    EGVectorLaplaceFormP1P1_new_01()

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

      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      real_t a_1_0 = 0;
      real_t a_1_1 = 0;
      real_t a_1_2 = 0;
      real_t a_2_0 = 0;
      real_t a_2_1 = 0;
      real_t a_2_2 = 0;
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

      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      real_t a_1_0 = 0;
      real_t a_1_1 = 0;
      real_t a_1_2 = 0;
      real_t a_2_0 = 0;
      real_t a_2_1 = 0;
      real_t a_2_2 = 0;
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

      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      real_t a_1_0 = 0;
      real_t a_1_1 = 0;
      real_t a_1_2 = 0;
      real_t a_2_0 = 0;
      real_t a_2_1 = 0;
      real_t a_2_2 = 0;
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

      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      real_t a_1_0 = 0;
      real_t a_1_1 = 0;
      real_t a_1_2 = 0;
      real_t a_2_0 = 0;
      real_t a_2_1 = 0;
      real_t a_2_2 = 0;
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




class EGVectorLaplaceFormP1P1_new_11 : public hyteg::dg::DGForm2D
{
 public:
    EGVectorLaplaceFormP1P1_new_11()
: callback_Scalar_Variable_Coefficient_2D_g1 ([](const Point3D & p) -> real_t { return 0.; })
    {}

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
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0);
      real_t tmp_5 = 1.0 / (tmp_4);
      real_t tmp_6 = tmp_1*tmp_5;
      real_t tmp_7 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_8 = tmp_5*tmp_7;
      real_t tmp_9 = -tmp_6 - tmp_8;
      real_t tmp_10 = tmp_3*tmp_5;
      real_t tmp_11 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = -tmp_10 - tmp_12;
      real_t tmp_14 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_15 = tmp_14*((tmp_13*tmp_13) + (tmp_9*tmp_9));
      real_t tmp_16 = tmp_14*(tmp_10*tmp_13 + tmp_8*tmp_9);
      real_t tmp_17 = 0.5*tmp_16;
      real_t tmp_18 = tmp_14*(tmp_12*tmp_13 + tmp_6*tmp_9);
      real_t tmp_19 = 0.5*tmp_18;
      real_t tmp_20 = 1.0 / (tmp_4*tmp_4);
      real_t tmp_21 = tmp_14*(tmp_20*(tmp_3*tmp_3) + tmp_20*(tmp_7*tmp_7));
      real_t tmp_22 = tmp_14*(tmp_1*tmp_20*tmp_7 + tmp_11*tmp_20*tmp_3);
      real_t tmp_23 = 0.5*tmp_22;
      real_t tmp_24 = tmp_14*((tmp_1*tmp_1)*tmp_20 + (tmp_11*tmp_11)*tmp_20);
      real_t a_0_0 = 0.5*tmp_15;
      real_t a_0_1 = tmp_17;
      real_t a_0_2 = tmp_19;
      real_t a_1_0 = tmp_17;
      real_t a_1_1 = 0.5*tmp_21;
      real_t a_1_2 = tmp_23;
      real_t a_2_0 = tmp_19;
      real_t a_2_1 = tmp_23;
      real_t a_2_2 = 0.5*tmp_24;
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

      real_t tmp_0 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_1 = -p_affine_0_1;
      real_t tmp_2 = p_affine_6_1 + tmp_1;
      real_t tmp_3 = 0.046910077030668018*tmp_0 + tmp_2;
      real_t tmp_4 = -p_affine_0_0;
      real_t tmp_5 = p_affine_1_0 + tmp_4;
      real_t tmp_6 = p_affine_2_1 + tmp_1;
      real_t tmp_7 = 1.0 / (tmp_5*tmp_6 - (p_affine_1_1 + tmp_1)*(p_affine_2_0 + tmp_4));
      real_t tmp_8 = tmp_5*tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = tmp_7*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_11 = tmp_10*tmp_3;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_4;
      real_t tmp_14 = 0.046910077030668018*tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6*tmp_7;
      real_t tmp_16 = tmp_14*tmp_15;
      real_t tmp_17 = tmp_7*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_18 = tmp_14*tmp_17;
      real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_12*tmp_12), 1.0/2.0));
      real_t tmp_21 = 6/tmp_20;
      real_t tmp_22 = p_affine_10_0*(-tmp_15 - tmp_17) + p_affine_10_1*(-tmp_10 - tmp_8);
      real_t tmp_23 = 1.0*tmp_22;
      real_t tmp_24 = 0.11846344252809471*tmp_20;
      real_t tmp_25 = 0.23076534494715845*tmp_0 + tmp_2;
      real_t tmp_26 = tmp_25*tmp_8;
      real_t tmp_27 = tmp_10*tmp_25;
      real_t tmp_28 = 0.23076534494715845*tmp_12 + tmp_13;
      real_t tmp_29 = tmp_15*tmp_28;
      real_t tmp_30 = tmp_17*tmp_28;
      real_t tmp_31 = -tmp_26 - tmp_27 - tmp_29 - tmp_30 + 1;
      real_t tmp_32 = 0.2393143352496831*tmp_20;
      real_t tmp_33 = 0.5*tmp_0 + tmp_2;
      real_t tmp_34 = tmp_33*tmp_8;
      real_t tmp_35 = tmp_10*tmp_33;
      real_t tmp_36 = 0.5*tmp_12 + tmp_13;
      real_t tmp_37 = tmp_15*tmp_36;
      real_t tmp_38 = tmp_17*tmp_36;
      real_t tmp_39 = -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1;
      real_t tmp_40 = 0.2844444444444445*tmp_20;
      real_t tmp_41 = 0.7692346550528415*tmp_0 + tmp_2;
      real_t tmp_42 = tmp_41*tmp_8;
      real_t tmp_43 = tmp_10*tmp_41;
      real_t tmp_44 = 0.7692346550528415*tmp_12 + tmp_13;
      real_t tmp_45 = tmp_15*tmp_44;
      real_t tmp_46 = tmp_17*tmp_44;
      real_t tmp_47 = -tmp_42 - tmp_43 - tmp_45 - tmp_46 + 1;
      real_t tmp_48 = 0.2393143352496831*tmp_20;
      real_t tmp_49 = 0.95308992296933193*tmp_0 + tmp_2;
      real_t tmp_50 = tmp_49*tmp_8;
      real_t tmp_51 = tmp_10*tmp_49;
      real_t tmp_52 = 0.95308992296933193*tmp_12 + tmp_13;
      real_t tmp_53 = tmp_15*tmp_52;
      real_t tmp_54 = tmp_17*tmp_52;
      real_t tmp_55 = -tmp_50 - tmp_51 - tmp_53 - tmp_54 + 1;
      real_t tmp_56 = 0.11846344252809471*tmp_20;
      real_t tmp_57 = tmp_11 + tmp_16;
      real_t tmp_58 = 0.5*tmp_22;
      real_t tmp_59 = p_affine_10_0*tmp_15 + p_affine_10_1*tmp_10;
      real_t tmp_60 = 0.5*tmp_59;
      real_t tmp_61 = tmp_19*tmp_21;
      real_t tmp_62 = tmp_27 + tmp_29;
      real_t tmp_63 = tmp_21*tmp_31;
      real_t tmp_64 = tmp_35 + tmp_37;
      real_t tmp_65 = tmp_21*tmp_39;
      real_t tmp_66 = tmp_43 + tmp_45;
      real_t tmp_67 = tmp_21*tmp_47;
      real_t tmp_68 = tmp_51 + tmp_53;
      real_t tmp_69 = tmp_21*tmp_55;
      real_t tmp_70 = tmp_24*(-tmp_19*tmp_60 - tmp_57*tmp_58 + tmp_57*tmp_61) + tmp_32*(-tmp_31*tmp_60 - tmp_58*tmp_62 + tmp_62*tmp_63) + tmp_40*(-tmp_39*tmp_60 - tmp_58*tmp_64 + tmp_64*tmp_65) + tmp_48*(-tmp_47*tmp_60 - tmp_58*tmp_66 + tmp_66*tmp_67) + tmp_56*(-tmp_55*tmp_60 - tmp_58*tmp_68 + tmp_68*tmp_69);
      real_t tmp_71 = tmp_18 + tmp_9;
      real_t tmp_72 = p_affine_10_0*tmp_17 + p_affine_10_1*tmp_8;
      real_t tmp_73 = 0.5*tmp_72;
      real_t tmp_74 = tmp_26 + tmp_30;
      real_t tmp_75 = tmp_34 + tmp_38;
      real_t tmp_76 = tmp_42 + tmp_46;
      real_t tmp_77 = tmp_50 + tmp_54;
      real_t tmp_78 = tmp_24*(-tmp_19*tmp_73 - tmp_58*tmp_71 + tmp_61*tmp_71) + tmp_32*(-tmp_31*tmp_73 - tmp_58*tmp_74 + tmp_63*tmp_74) + tmp_40*(-tmp_39*tmp_73 - tmp_58*tmp_75 + tmp_65*tmp_75) + tmp_48*(-tmp_47*tmp_73 - tmp_58*tmp_76 + tmp_67*tmp_76) + tmp_56*(-tmp_55*tmp_73 - tmp_58*tmp_77 + tmp_69*tmp_77);
      real_t tmp_79 = 1.0*tmp_59;
      real_t tmp_80 = tmp_24*(tmp_21*tmp_57*tmp_71 - tmp_57*tmp_73 - tmp_60*tmp_71) + tmp_32*(tmp_21*tmp_62*tmp_74 - tmp_60*tmp_74 - tmp_62*tmp_73) + tmp_40*(tmp_21*tmp_64*tmp_75 - tmp_60*tmp_75 - tmp_64*tmp_73) + tmp_48*(tmp_21*tmp_66*tmp_76 - tmp_60*tmp_76 - tmp_66*tmp_73) + tmp_56*(tmp_21*tmp_68*tmp_77 - tmp_60*tmp_77 - tmp_68*tmp_73);
      real_t tmp_81 = 1.0*tmp_72;
      real_t a_0_0 = tmp_24*((tmp_19*tmp_19)*tmp_21 - tmp_19*tmp_23) + tmp_32*(tmp_21*(tmp_31*tmp_31) - tmp_23*tmp_31) + tmp_40*(tmp_21*(tmp_39*tmp_39) - tmp_23*tmp_39) + tmp_48*(tmp_21*(tmp_47*tmp_47) - tmp_23*tmp_47) + tmp_56*(tmp_21*(tmp_55*tmp_55) - tmp_23*tmp_55);
      real_t a_0_1 = tmp_70;
      real_t a_0_2 = tmp_78;
      real_t a_1_0 = tmp_70;
      real_t a_1_1 = tmp_24*(tmp_21*(tmp_57*tmp_57) - tmp_57*tmp_79) + tmp_32*(tmp_21*(tmp_62*tmp_62) - tmp_62*tmp_79) + tmp_40*(tmp_21*(tmp_64*tmp_64) - tmp_64*tmp_79) + tmp_48*(tmp_21*(tmp_66*tmp_66) - tmp_66*tmp_79) + tmp_56*(tmp_21*(tmp_68*tmp_68) - tmp_68*tmp_79);
      real_t a_1_2 = tmp_80;
      real_t a_2_0 = tmp_78;
      real_t a_2_1 = tmp_80;
      real_t a_2_2 = tmp_24*(tmp_21*(tmp_71*tmp_71) - tmp_71*tmp_81) + tmp_32*(tmp_21*(tmp_74*tmp_74) - tmp_74*tmp_81) + tmp_40*(tmp_21*(tmp_75*tmp_75) - tmp_75*tmp_81) + tmp_48*(tmp_21*(tmp_76*tmp_76) - tmp_76*tmp_81) + tmp_56*(tmp_21*(tmp_77*tmp_77) - tmp_77*tmp_81);
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

      real_t tmp_0 = -p_affine_3_1;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = p_affine_6_1 + 0.046910077030668018*tmp_1;
      real_t tmp_3 = tmp_0 + tmp_2;
      real_t tmp_4 = -p_affine_3_0;
      real_t tmp_5 = p_affine_4_0 + tmp_4;
      real_t tmp_6 = p_affine_5_1 + tmp_0;
      real_t tmp_7 = 1.0 / (tmp_5*tmp_6 - (p_affine_4_1 + tmp_0)*(p_affine_5_0 + tmp_4));
      real_t tmp_8 = tmp_5*tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = tmp_7*(p_affine_3_0 - p_affine_5_0);
      real_t tmp_11 = tmp_10*tmp_3;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + 0.046910077030668018*tmp_12;
      real_t tmp_14 = tmp_13 + tmp_4;
      real_t tmp_15 = tmp_6*tmp_7;
      real_t tmp_16 = tmp_14*tmp_15;
      real_t tmp_17 = tmp_7*(p_affine_3_1 - p_affine_4_1);
      real_t tmp_18 = tmp_14*tmp_17;
      real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20 = -p_affine_0_1;
      real_t tmp_21 = p_affine_2_1 + tmp_20;
      real_t tmp_22 = -p_affine_0_0;
      real_t tmp_23 = p_affine_1_0 + tmp_22;
      real_t tmp_24 = 1.0 / (tmp_21*tmp_23 - (p_affine_1_1 + tmp_20)*(p_affine_2_0 + tmp_22));
      real_t tmp_25 = tmp_21*tmp_24;
      real_t tmp_26 = tmp_24*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_27 = tmp_23*tmp_24;
      real_t tmp_28 = tmp_24*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_29 = 0.5*p_affine_10_0*(-tmp_25 - tmp_26) + 0.5*p_affine_10_1*(-tmp_27 - tmp_28);
      real_t tmp_30 = tmp_2 + tmp_20;
      real_t tmp_31 = tmp_27*tmp_30;
      real_t tmp_32 = tmp_28*tmp_30;
      real_t tmp_33 = tmp_13 + tmp_22;
      real_t tmp_34 = tmp_25*tmp_33;
      real_t tmp_35 = tmp_26*tmp_33;
      real_t tmp_36 = -tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1;
      real_t tmp_37 = 0.5*p_affine_10_0*(-tmp_15 - tmp_17) + 0.5*p_affine_10_1*(-tmp_10 - tmp_8);
      real_t tmp_38 = std::abs(std::pow((tmp_1*tmp_1) + (tmp_12*tmp_12), 1.0/2.0));
      real_t tmp_39 = 6/tmp_38;
      real_t tmp_40 = tmp_36*tmp_39;
      real_t tmp_41 = 0.11846344252809471*tmp_38;
      real_t tmp_42 = p_affine_6_1 + 0.23076534494715845*tmp_1;
      real_t tmp_43 = tmp_0 + tmp_42;
      real_t tmp_44 = tmp_43*tmp_8;
      real_t tmp_45 = tmp_10*tmp_43;
      real_t tmp_46 = p_affine_6_0 + 0.23076534494715845*tmp_12;
      real_t tmp_47 = tmp_4 + tmp_46;
      real_t tmp_48 = tmp_15*tmp_47;
      real_t tmp_49 = tmp_17*tmp_47;
      real_t tmp_50 = -tmp_44 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_51 = tmp_20 + tmp_42;
      real_t tmp_52 = tmp_27*tmp_51;
      real_t tmp_53 = tmp_28*tmp_51;
      real_t tmp_54 = tmp_22 + tmp_46;
      real_t tmp_55 = tmp_25*tmp_54;
      real_t tmp_56 = tmp_26*tmp_54;
      real_t tmp_57 = -tmp_52 - tmp_53 - tmp_55 - tmp_56 + 1;
      real_t tmp_58 = tmp_39*tmp_57;
      real_t tmp_59 = 0.2393143352496831*tmp_38;
      real_t tmp_60 = p_affine_6_1 + 0.5*tmp_1;
      real_t tmp_61 = tmp_0 + tmp_60;
      real_t tmp_62 = tmp_61*tmp_8;
      real_t tmp_63 = tmp_10*tmp_61;
      real_t tmp_64 = p_affine_6_0 + 0.5*tmp_12;
      real_t tmp_65 = tmp_4 + tmp_64;
      real_t tmp_66 = tmp_15*tmp_65;
      real_t tmp_67 = tmp_17*tmp_65;
      real_t tmp_68 = -tmp_62 - tmp_63 - tmp_66 - tmp_67 + 1;
      real_t tmp_69 = tmp_20 + tmp_60;
      real_t tmp_70 = tmp_27*tmp_69;
      real_t tmp_71 = tmp_28*tmp_69;
      real_t tmp_72 = tmp_22 + tmp_64;
      real_t tmp_73 = tmp_25*tmp_72;
      real_t tmp_74 = tmp_26*tmp_72;
      real_t tmp_75 = -tmp_70 - tmp_71 - tmp_73 - tmp_74 + 1;
      real_t tmp_76 = tmp_39*tmp_75;
      real_t tmp_77 = 0.2844444444444445*tmp_38;
      real_t tmp_78 = p_affine_6_1 + 0.7692346550528415*tmp_1;
      real_t tmp_79 = tmp_0 + tmp_78;
      real_t tmp_80 = tmp_79*tmp_8;
      real_t tmp_81 = tmp_10*tmp_79;
      real_t tmp_82 = p_affine_6_0 + 0.7692346550528415*tmp_12;
      real_t tmp_83 = tmp_4 + tmp_82;
      real_t tmp_84 = tmp_15*tmp_83;
      real_t tmp_85 = tmp_17*tmp_83;
      real_t tmp_86 = -tmp_80 - tmp_81 - tmp_84 - tmp_85 + 1;
      real_t tmp_87 = tmp_20 + tmp_78;
      real_t tmp_88 = tmp_27*tmp_87;
      real_t tmp_89 = tmp_28*tmp_87;
      real_t tmp_90 = tmp_22 + tmp_82;
      real_t tmp_91 = tmp_25*tmp_90;
      real_t tmp_92 = tmp_26*tmp_90;
      real_t tmp_93 = -tmp_88 - tmp_89 - tmp_91 - tmp_92 + 1;
      real_t tmp_94 = tmp_39*tmp_93;
      real_t tmp_95 = 0.2393143352496831*tmp_38;
      real_t tmp_96 = p_affine_6_1 + 0.95308992296933193*tmp_1;
      real_t tmp_97 = tmp_0 + tmp_96;
      real_t tmp_98 = tmp_8*tmp_97;
      real_t tmp_99 = tmp_10*tmp_97;
      real_t tmp_100 = p_affine_6_0 + 0.95308992296933193*tmp_12;
      real_t tmp_101 = tmp_100 + tmp_4;
      real_t tmp_102 = tmp_101*tmp_15;
      real_t tmp_103 = tmp_101*tmp_17;
      real_t tmp_104 = -tmp_102 - tmp_103 - tmp_98 - tmp_99 + 1;
      real_t tmp_105 = tmp_20 + tmp_96;
      real_t tmp_106 = tmp_105*tmp_27;
      real_t tmp_107 = tmp_105*tmp_28;
      real_t tmp_108 = tmp_100 + tmp_22;
      real_t tmp_109 = tmp_108*tmp_25;
      real_t tmp_110 = tmp_108*tmp_26;
      real_t tmp_111 = -tmp_106 - tmp_107 - tmp_109 - tmp_110 + 1;
      real_t tmp_112 = tmp_111*tmp_39;
      real_t tmp_113 = 0.11846344252809471*tmp_38;
      real_t tmp_114 = tmp_11 + tmp_16;
      real_t tmp_115 = 0.5*p_affine_10_0*tmp_15 + 0.5*p_affine_10_1*tmp_10;
      real_t tmp_116 = tmp_45 + tmp_48;
      real_t tmp_117 = tmp_63 + tmp_66;
      real_t tmp_118 = tmp_81 + tmp_84;
      real_t tmp_119 = tmp_102 + tmp_99;
      real_t tmp_120 = tmp_18 + tmp_9;
      real_t tmp_121 = 0.5*p_affine_10_0*tmp_17 + 0.5*p_affine_10_1*tmp_8;
      real_t tmp_122 = tmp_44 + tmp_49;
      real_t tmp_123 = tmp_62 + tmp_67;
      real_t tmp_124 = tmp_80 + tmp_85;
      real_t tmp_125 = tmp_103 + tmp_98;
      real_t tmp_126 = tmp_32 + tmp_34;
      real_t tmp_127 = 0.5*p_affine_10_0*tmp_25 + 0.5*p_affine_10_1*tmp_28;
      real_t tmp_128 = tmp_126*tmp_39;
      real_t tmp_129 = tmp_53 + tmp_55;
      real_t tmp_130 = tmp_129*tmp_39;
      real_t tmp_131 = tmp_71 + tmp_73;
      real_t tmp_132 = tmp_131*tmp_39;
      real_t tmp_133 = tmp_89 + tmp_91;
      real_t tmp_134 = tmp_133*tmp_39;
      real_t tmp_135 = tmp_107 + tmp_109;
      real_t tmp_136 = tmp_135*tmp_39;
      real_t tmp_137 = tmp_31 + tmp_35;
      real_t tmp_138 = 0.5*p_affine_10_0*tmp_26 + 0.5*p_affine_10_1*tmp_27;
      real_t tmp_139 = tmp_137*tmp_39;
      real_t tmp_140 = tmp_52 + tmp_56;
      real_t tmp_141 = tmp_140*tmp_39;
      real_t tmp_142 = tmp_70 + tmp_74;
      real_t tmp_143 = tmp_142*tmp_39;
      real_t tmp_144 = tmp_88 + tmp_92;
      real_t tmp_145 = tmp_144*tmp_39;
      real_t tmp_146 = tmp_106 + tmp_110;
      real_t tmp_147 = tmp_146*tmp_39;
      real_t a_0_0 = tmp_113*(-tmp_104*tmp_112 + tmp_104*tmp_29 - tmp_111*tmp_37) + tmp_41*(tmp_19*tmp_29 - tmp_19*tmp_40 - tmp_36*tmp_37) + tmp_59*(tmp_29*tmp_50 - tmp_37*tmp_57 - tmp_50*tmp_58) + tmp_77*(tmp_29*tmp_68 - tmp_37*tmp_75 - tmp_68*tmp_76) + tmp_95*(tmp_29*tmp_86 - tmp_37*tmp_93 - tmp_86*tmp_94);
      real_t a_0_1 = tmp_113*(-tmp_111*tmp_115 - tmp_112*tmp_119 + tmp_119*tmp_29) + tmp_41*(tmp_114*tmp_29 - tmp_114*tmp_40 - tmp_115*tmp_36) + tmp_59*(-tmp_115*tmp_57 + tmp_116*tmp_29 - tmp_116*tmp_58) + tmp_77*(-tmp_115*tmp_75 + tmp_117*tmp_29 - tmp_117*tmp_76) + tmp_95*(-tmp_115*tmp_93 + tmp_118*tmp_29 - tmp_118*tmp_94);
      real_t a_0_2 = tmp_113*(-tmp_111*tmp_121 - tmp_112*tmp_125 + tmp_125*tmp_29) + tmp_41*(tmp_120*tmp_29 - tmp_120*tmp_40 - tmp_121*tmp_36) + tmp_59*(-tmp_121*tmp_57 + tmp_122*tmp_29 - tmp_122*tmp_58) + tmp_77*(-tmp_121*tmp_75 + tmp_123*tmp_29 - tmp_123*tmp_76) + tmp_95*(-tmp_121*tmp_93 + tmp_124*tmp_29 - tmp_124*tmp_94);
      real_t a_1_0 = tmp_113*(tmp_104*tmp_127 - tmp_104*tmp_136 - tmp_135*tmp_37) + tmp_41*(-tmp_126*tmp_37 + tmp_127*tmp_19 - tmp_128*tmp_19) + tmp_59*(tmp_127*tmp_50 - tmp_129*tmp_37 - tmp_130*tmp_50) + tmp_77*(tmp_127*tmp_68 - tmp_131*tmp_37 - tmp_132*tmp_68) + tmp_95*(tmp_127*tmp_86 - tmp_133*tmp_37 - tmp_134*tmp_86);
      real_t a_1_1 = tmp_113*(-tmp_115*tmp_135 + tmp_119*tmp_127 - tmp_119*tmp_136) + tmp_41*(tmp_114*tmp_127 - tmp_114*tmp_128 - tmp_115*tmp_126) + tmp_59*(-tmp_115*tmp_129 + tmp_116*tmp_127 - tmp_116*tmp_130) + tmp_77*(-tmp_115*tmp_131 + tmp_117*tmp_127 - tmp_117*tmp_132) + tmp_95*(-tmp_115*tmp_133 + tmp_118*tmp_127 - tmp_118*tmp_134);
      real_t a_1_2 = tmp_113*(-tmp_121*tmp_135 + tmp_125*tmp_127 - tmp_125*tmp_136) + tmp_41*(tmp_120*tmp_127 - tmp_120*tmp_128 - tmp_121*tmp_126) + tmp_59*(-tmp_121*tmp_129 + tmp_122*tmp_127 - tmp_122*tmp_130) + tmp_77*(-tmp_121*tmp_131 + tmp_123*tmp_127 - tmp_123*tmp_132) + tmp_95*(-tmp_121*tmp_133 + tmp_124*tmp_127 - tmp_124*tmp_134);
      real_t a_2_0 = tmp_113*(tmp_104*tmp_138 - tmp_104*tmp_147 - tmp_146*tmp_37) + tmp_41*(-tmp_137*tmp_37 + tmp_138*tmp_19 - tmp_139*tmp_19) + tmp_59*(tmp_138*tmp_50 - tmp_140*tmp_37 - tmp_141*tmp_50) + tmp_77*(tmp_138*tmp_68 - tmp_142*tmp_37 - tmp_143*tmp_68) + tmp_95*(tmp_138*tmp_86 - tmp_144*tmp_37 - tmp_145*tmp_86);
      real_t a_2_1 = tmp_113*(-tmp_115*tmp_146 + tmp_119*tmp_138 - tmp_119*tmp_147) + tmp_41*(tmp_114*tmp_138 - tmp_114*tmp_139 - tmp_115*tmp_137) + tmp_59*(-tmp_115*tmp_140 + tmp_116*tmp_138 - tmp_116*tmp_141) + tmp_77*(-tmp_115*tmp_142 + tmp_117*tmp_138 - tmp_117*tmp_143) + tmp_95*(-tmp_115*tmp_144 + tmp_118*tmp_138 - tmp_118*tmp_145);
      real_t a_2_2 = tmp_113*(-tmp_121*tmp_146 + tmp_125*tmp_138 - tmp_125*tmp_147) + tmp_41*(tmp_120*tmp_138 - tmp_120*tmp_139 - tmp_121*tmp_137) + tmp_59*(-tmp_121*tmp_140 + tmp_122*tmp_138 - tmp_122*tmp_141) + tmp_77*(-tmp_121*tmp_142 + tmp_123*tmp_138 - tmp_123*tmp_143) + tmp_95*(-tmp_121*tmp_144 + tmp_124*tmp_138 - tmp_124*tmp_145);
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

      real_t tmp_0 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_1 = -p_affine_0_1;
      real_t tmp_2 = p_affine_6_1 + tmp_1;
      real_t tmp_3 = 0.046910077030668018*tmp_0 + tmp_2;
      real_t tmp_4 = -p_affine_0_0;
      real_t tmp_5 = p_affine_1_0 + tmp_4;
      real_t tmp_6 = p_affine_2_1 + tmp_1;
      real_t tmp_7 = 1.0 / (tmp_5*tmp_6 - (p_affine_1_1 + tmp_1)*(p_affine_2_0 + tmp_4));
      real_t tmp_8 = tmp_5*tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = tmp_7*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_11 = tmp_10*tmp_3;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_4;
      real_t tmp_14 = 0.046910077030668018*tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6*tmp_7;
      real_t tmp_16 = tmp_14*tmp_15;
      real_t tmp_17 = tmp_7*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_18 = tmp_14*tmp_17;
      real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_12*tmp_12), 1.0/2.0));
      real_t tmp_21 = 6/tmp_20;
      real_t tmp_22 = p_affine_10_0*(-tmp_15 - tmp_17) + p_affine_10_1*(-tmp_10 - tmp_8);
      real_t tmp_23 = 2*tmp_22;
      real_t tmp_24 = 0.11846344252809471*tmp_20;
      real_t tmp_25 = 0.23076534494715845*tmp_0 + tmp_2;
      real_t tmp_26 = tmp_25*tmp_8;
      real_t tmp_27 = tmp_10*tmp_25;
      real_t tmp_28 = 0.23076534494715845*tmp_12 + tmp_13;
      real_t tmp_29 = tmp_15*tmp_28;
      real_t tmp_30 = tmp_17*tmp_28;
      real_t tmp_31 = -tmp_26 - tmp_27 - tmp_29 - tmp_30 + 1;
      real_t tmp_32 = 0.2393143352496831*tmp_20;
      real_t tmp_33 = 0.5*tmp_0 + tmp_2;
      real_t tmp_34 = tmp_33*tmp_8;
      real_t tmp_35 = tmp_10*tmp_33;
      real_t tmp_36 = 0.5*tmp_12 + tmp_13;
      real_t tmp_37 = tmp_15*tmp_36;
      real_t tmp_38 = tmp_17*tmp_36;
      real_t tmp_39 = -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1;
      real_t tmp_40 = 0.2844444444444445*tmp_20;
      real_t tmp_41 = 0.7692346550528415*tmp_0 + tmp_2;
      real_t tmp_42 = tmp_41*tmp_8;
      real_t tmp_43 = tmp_10*tmp_41;
      real_t tmp_44 = 0.7692346550528415*tmp_12 + tmp_13;
      real_t tmp_45 = tmp_15*tmp_44;
      real_t tmp_46 = tmp_17*tmp_44;
      real_t tmp_47 = -tmp_42 - tmp_43 - tmp_45 - tmp_46 + 1;
      real_t tmp_48 = 0.2393143352496831*tmp_20;
      real_t tmp_49 = 0.95308992296933193*tmp_0 + tmp_2;
      real_t tmp_50 = tmp_49*tmp_8;
      real_t tmp_51 = tmp_10*tmp_49;
      real_t tmp_52 = 0.95308992296933193*tmp_12 + tmp_13;
      real_t tmp_53 = tmp_15*tmp_52;
      real_t tmp_54 = tmp_17*tmp_52;
      real_t tmp_55 = -tmp_50 - tmp_51 - tmp_53 - tmp_54 + 1;
      real_t tmp_56 = 0.11846344252809471*tmp_20;
      real_t tmp_57 = tmp_11 + tmp_16;
      real_t tmp_58 = p_affine_10_0*tmp_15 + p_affine_10_1*tmp_10;
      real_t tmp_59 = tmp_19*tmp_21;
      real_t tmp_60 = tmp_27 + tmp_29;
      real_t tmp_61 = tmp_21*tmp_31;
      real_t tmp_62 = tmp_35 + tmp_37;
      real_t tmp_63 = tmp_21*tmp_39;
      real_t tmp_64 = tmp_43 + tmp_45;
      real_t tmp_65 = tmp_21*tmp_47;
      real_t tmp_66 = tmp_51 + tmp_53;
      real_t tmp_67 = tmp_21*tmp_55;
      real_t tmp_68 = tmp_24*(-tmp_19*tmp_58 - tmp_22*tmp_57 + tmp_57*tmp_59) + tmp_32*(-tmp_22*tmp_60 - tmp_31*tmp_58 + tmp_60*tmp_61) + tmp_40*(-tmp_22*tmp_62 - tmp_39*tmp_58 + tmp_62*tmp_63) + tmp_48*(-tmp_22*tmp_64 - tmp_47*tmp_58 + tmp_64*tmp_65) + tmp_56*(-tmp_22*tmp_66 - tmp_55*tmp_58 + tmp_66*tmp_67);
      real_t tmp_69 = tmp_18 + tmp_9;
      real_t tmp_70 = p_affine_10_0*tmp_17 + p_affine_10_1*tmp_8;
      real_t tmp_71 = tmp_26 + tmp_30;
      real_t tmp_72 = tmp_34 + tmp_38;
      real_t tmp_73 = tmp_42 + tmp_46;
      real_t tmp_74 = tmp_50 + tmp_54;
      real_t tmp_75 = tmp_24*(-tmp_19*tmp_70 - tmp_22*tmp_69 + tmp_59*tmp_69) + tmp_32*(-tmp_22*tmp_71 - tmp_31*tmp_70 + tmp_61*tmp_71) + tmp_40*(-tmp_22*tmp_72 - tmp_39*tmp_70 + tmp_63*tmp_72) + tmp_48*(-tmp_22*tmp_73 - tmp_47*tmp_70 + tmp_65*tmp_73) + tmp_56*(-tmp_22*tmp_74 - tmp_55*tmp_70 + tmp_67*tmp_74);
      real_t tmp_76 = 2*tmp_58;
      real_t tmp_77 = tmp_24*(tmp_21*tmp_57*tmp_69 - tmp_57*tmp_70 - tmp_58*tmp_69) + tmp_32*(tmp_21*tmp_60*tmp_71 - tmp_58*tmp_71 - tmp_60*tmp_70) + tmp_40*(tmp_21*tmp_62*tmp_72 - tmp_58*tmp_72 - tmp_62*tmp_70) + tmp_48*(tmp_21*tmp_64*tmp_73 - tmp_58*tmp_73 - tmp_64*tmp_70) + tmp_56*(tmp_21*tmp_66*tmp_74 - tmp_58*tmp_74 - tmp_66*tmp_70);
      real_t tmp_78 = 2*tmp_70;
      real_t a_0_0 = tmp_24*((tmp_19*tmp_19)*tmp_21 - tmp_19*tmp_23) + tmp_32*(tmp_21*(tmp_31*tmp_31) - tmp_23*tmp_31) + tmp_40*(tmp_21*(tmp_39*tmp_39) - tmp_23*tmp_39) + tmp_48*(tmp_21*(tmp_47*tmp_47) - tmp_23*tmp_47) + tmp_56*(tmp_21*(tmp_55*tmp_55) - tmp_23*tmp_55);
      real_t a_0_1 = tmp_68;
      real_t a_0_2 = tmp_75;
      real_t a_1_0 = tmp_68;
      real_t a_1_1 = tmp_24*(tmp_21*(tmp_57*tmp_57) - tmp_57*tmp_76) + tmp_32*(tmp_21*(tmp_60*tmp_60) - tmp_60*tmp_76) + tmp_40*(tmp_21*(tmp_62*tmp_62) - tmp_62*tmp_76) + tmp_48*(tmp_21*(tmp_64*tmp_64) - tmp_64*tmp_76) + tmp_56*(tmp_21*(tmp_66*tmp_66) - tmp_66*tmp_76);
      real_t a_1_2 = tmp_77;
      real_t a_2_0 = tmp_75;
      real_t a_2_1 = tmp_77;
      real_t a_2_2 = tmp_24*(tmp_21*(tmp_69*tmp_69) - tmp_69*tmp_78) + tmp_32*(tmp_21*(tmp_71*tmp_71) - tmp_71*tmp_78) + tmp_40*(tmp_21*(tmp_72*tmp_72) - tmp_72*tmp_78) + tmp_48*(tmp_21*(tmp_73*tmp_73) - tmp_73*tmp_78) + tmp_56*(tmp_21*(tmp_74*tmp_74) - tmp_74*tmp_78);
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
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_g1_out0_id4 = 0;
      Scalar_Variable_Coefficient_2D_g1( 0.95308992296933193*p_affine_6_0 + 0.046910077030668018*p_affine_7_0, 0.95308992296933193*p_affine_6_1 + 0.046910077030668018*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id0 );
      Scalar_Variable_Coefficient_2D_g1( 0.7692346550528415*p_affine_6_0 + 0.23076534494715845*p_affine_7_0, 0.7692346550528415*p_affine_6_1 + 0.23076534494715845*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id1 );
      Scalar_Variable_Coefficient_2D_g1( 0.5*p_affine_6_0 + 0.5*p_affine_7_0, 0.5*p_affine_6_1 + 0.5*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id2 );
      Scalar_Variable_Coefficient_2D_g1( 0.2307653449471585*p_affine_6_0 + 0.7692346550528415*p_affine_7_0, 0.2307653449471585*p_affine_6_1 + 0.7692346550528415*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id3 );
      Scalar_Variable_Coefficient_2D_g1( 0.046910077030668074*p_affine_6_0 + 0.95308992296933193*p_affine_7_0, 0.046910077030668074*p_affine_6_1 + 0.95308992296933193*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g1_out0_id4 );
      real_t tmp_0 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_1 = -p_affine_0_1;
      real_t tmp_2 = p_affine_6_1 + tmp_1;
      real_t tmp_3 = 0.046910077030668018*tmp_0 + tmp_2;
      real_t tmp_4 = -p_affine_0_0;
      real_t tmp_5 = p_affine_1_0 + tmp_4;
      real_t tmp_6 = p_affine_2_1 + tmp_1;
      real_t tmp_7 = 1.0 / (tmp_5*tmp_6 - (p_affine_1_1 + tmp_1)*(p_affine_2_0 + tmp_4));
      real_t tmp_8 = tmp_5*tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = tmp_7*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_11 = tmp_10*tmp_3;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_4;
      real_t tmp_14 = 0.046910077030668018*tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6*tmp_7;
      real_t tmp_16 = tmp_14*tmp_15;
      real_t tmp_17 = tmp_7*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_18 = tmp_14*tmp_17;
      real_t tmp_19 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_12*tmp_12), 1.0/2.0));
      real_t tmp_20 = 6/tmp_19;
      real_t tmp_21 = -p_affine_10_0*(-tmp_15 - tmp_17) - p_affine_10_1*(-tmp_10 - tmp_8);
      real_t tmp_22 = 0.11846344252809471*Scalar_Variable_Coefficient_2D_g1_out0_id0*tmp_19;
      real_t tmp_23 = 0.23076534494715845*tmp_0 + tmp_2;
      real_t tmp_24 = tmp_23*tmp_8;
      real_t tmp_25 = tmp_10*tmp_23;
      real_t tmp_26 = 0.23076534494715845*tmp_12 + tmp_13;
      real_t tmp_27 = tmp_15*tmp_26;
      real_t tmp_28 = tmp_17*tmp_26;
      real_t tmp_29 = 0.2393143352496831*Scalar_Variable_Coefficient_2D_g1_out0_id1*tmp_19;
      real_t tmp_30 = 0.5*tmp_0 + tmp_2;
      real_t tmp_31 = tmp_30*tmp_8;
      real_t tmp_32 = tmp_10*tmp_30;
      real_t tmp_33 = 0.5*tmp_12 + tmp_13;
      real_t tmp_34 = tmp_15*tmp_33;
      real_t tmp_35 = tmp_17*tmp_33;
      real_t tmp_36 = 0.2844444444444445*Scalar_Variable_Coefficient_2D_g1_out0_id2*tmp_19;
      real_t tmp_37 = 0.7692346550528415*tmp_0 + tmp_2;
      real_t tmp_38 = tmp_37*tmp_8;
      real_t tmp_39 = tmp_10*tmp_37;
      real_t tmp_40 = 0.7692346550528415*tmp_12 + tmp_13;
      real_t tmp_41 = tmp_15*tmp_40;
      real_t tmp_42 = tmp_17*tmp_40;
      real_t tmp_43 = 0.2393143352496831*Scalar_Variable_Coefficient_2D_g1_out0_id3*tmp_19;
      real_t tmp_44 = 0.95308992296933193*tmp_0 + tmp_2;
      real_t tmp_45 = tmp_44*tmp_8;
      real_t tmp_46 = tmp_10*tmp_44;
      real_t tmp_47 = 0.95308992296933193*tmp_12 + tmp_13;
      real_t tmp_48 = tmp_15*tmp_47;
      real_t tmp_49 = tmp_17*tmp_47;
      real_t tmp_50 = 0.11846344252809471*Scalar_Variable_Coefficient_2D_g1_out0_id4*tmp_19;
      real_t tmp_51 = -p_affine_10_0*tmp_15 - p_affine_10_1*tmp_10;
      real_t tmp_52 = -p_affine_10_0*tmp_17 - p_affine_10_1*tmp_8;
      real_t a_0_0 = tmp_22*(tmp_20*(-tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1) + tmp_21) + tmp_29*(tmp_20*(-tmp_24 - tmp_25 - tmp_27 - tmp_28 + 1) + tmp_21) + tmp_36*(tmp_20*(-tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1) + tmp_21) + tmp_43*(tmp_20*(-tmp_38 - tmp_39 - tmp_41 - tmp_42 + 1) + tmp_21) + tmp_50*(tmp_20*(-tmp_45 - tmp_46 - tmp_48 - tmp_49 + 1) + tmp_21);
      real_t a_1_0 = tmp_22*(tmp_20*(tmp_11 + tmp_16) + tmp_51) + tmp_29*(tmp_20*(tmp_25 + tmp_27) + tmp_51) + tmp_36*(tmp_20*(tmp_32 + tmp_34) + tmp_51) + tmp_43*(tmp_20*(tmp_39 + tmp_41) + tmp_51) + tmp_50*(tmp_20*(tmp_46 + tmp_48) + tmp_51);
      real_t a_2_0 = tmp_22*(tmp_20*(tmp_18 + tmp_9) + tmp_52) + tmp_29*(tmp_20*(tmp_24 + tmp_28) + tmp_52) + tmp_36*(tmp_20*(tmp_31 + tmp_35) + tmp_52) + tmp_43*(tmp_20*(tmp_38 + tmp_42) + tmp_52) + tmp_50*(tmp_20*(tmp_45 + tmp_49) + tmp_52);
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
   }

public:

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g1;
};




class EGVectorLaplaceFormP1EDG_new_0 : public hyteg::dg::DGForm2D
{
 public:
    EGVectorLaplaceFormP1EDG_new_0()

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
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = 1.0 / (tmp_4 - tmp_5*(p_affine_1_1 + tmp_2));
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = tmp_6*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_9 = tmp_1*tmp_8 + tmp_5*tmp_7;
      real_t tmp_10 = tmp_3*tmp_6;
      real_t tmp_11 = tmp_6*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_12 = tmp_11*tmp_5 + tmp_4*tmp_6;
      real_t tmp_13 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_14 = tmp_13*(tmp_12*(-tmp_10 - tmp_11) + tmp_9*(-tmp_7 - tmp_8));
      real_t tmp_15 = tmp_13*(tmp_10*tmp_12 + tmp_8*tmp_9);
      real_t tmp_16 = tmp_13*(tmp_11*tmp_12 + tmp_7*tmp_9);
      real_t a_0_0 = 0.5*tmp_14;
      real_t a_1_0 = 0.5*tmp_15;
      real_t a_2_0 = 0.5*tmp_16;
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
      real_t tmp_2 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_6_1 + tmp_3;
      real_t tmp_5 = 0.046910077030668018*tmp_2 + tmp_4;
      real_t tmp_6 = p_affine_2_1 + tmp_3;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = p_affine_2_0 + tmp_0;
      real_t tmp_9 = 1.0 / (tmp_7 - tmp_8*(p_affine_1_1 + tmp_3));
      real_t tmp_10 = tmp_9*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_11 = tmp_10*tmp_5;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_0;
      real_t tmp_14 = 0.046910077030668018*tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6*tmp_9;
      real_t tmp_16 = tmp_14*tmp_15;
      real_t tmp_17 = tmp_11 + tmp_16;
      real_t tmp_18 = tmp_1*tmp_9;
      real_t tmp_19 = tmp_18*tmp_5;
      real_t tmp_20 = tmp_9*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_21 = tmp_14*tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_1*(tmp_17 - 1.0/3.0) + tmp_8*(tmp_22 - 1.0/3.0);
      real_t tmp_24 = 0.5*p_affine_10_0*(-tmp_15 - tmp_20) + 0.5*p_affine_10_1*(-tmp_10 - tmp_18);
      real_t tmp_25 = -tmp_11 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_26 = 0.5*p_affine_10_0*(tmp_20*tmp_8 + tmp_7*tmp_9) + 0.5*p_affine_10_1*(tmp_1*tmp_10 + tmp_18*tmp_8);
      real_t tmp_27 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_2*tmp_2), 1.0/2.0));
      real_t tmp_28 = 6/tmp_27;
      real_t tmp_29 = tmp_23*tmp_28;
      real_t tmp_30 = 0.11846344252809471*tmp_27;
      real_t tmp_31 = 0.23076534494715845*tmp_2 + tmp_4;
      real_t tmp_32 = tmp_10*tmp_31;
      real_t tmp_33 = 0.23076534494715845*tmp_12 + tmp_13;
      real_t tmp_34 = tmp_15*tmp_33;
      real_t tmp_35 = tmp_32 + tmp_34;
      real_t tmp_36 = tmp_18*tmp_31;
      real_t tmp_37 = tmp_20*tmp_33;
      real_t tmp_38 = tmp_36 + tmp_37;
      real_t tmp_39 = tmp_1*(tmp_35 - 1.0/3.0) + tmp_8*(tmp_38 - 1.0/3.0);
      real_t tmp_40 = -tmp_32 - tmp_34 - tmp_36 - tmp_37 + 1;
      real_t tmp_41 = tmp_28*tmp_39;
      real_t tmp_42 = 0.2393143352496831*tmp_27;
      real_t tmp_43 = 0.5*tmp_2 + tmp_4;
      real_t tmp_44 = tmp_10*tmp_43;
      real_t tmp_45 = 0.5*tmp_12 + tmp_13;
      real_t tmp_46 = tmp_15*tmp_45;
      real_t tmp_47 = tmp_44 + tmp_46;
      real_t tmp_48 = tmp_18*tmp_43;
      real_t tmp_49 = tmp_20*tmp_45;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_1*(tmp_47 - 1.0/3.0) + tmp_8*(tmp_50 - 1.0/3.0);
      real_t tmp_52 = -tmp_44 - tmp_46 - tmp_48 - tmp_49 + 1;
      real_t tmp_53 = tmp_28*tmp_51;
      real_t tmp_54 = 0.2844444444444445*tmp_27;
      real_t tmp_55 = 0.7692346550528415*tmp_2 + tmp_4;
      real_t tmp_56 = tmp_10*tmp_55;
      real_t tmp_57 = 0.7692346550528415*tmp_12 + tmp_13;
      real_t tmp_58 = tmp_15*tmp_57;
      real_t tmp_59 = tmp_56 + tmp_58;
      real_t tmp_60 = tmp_18*tmp_55;
      real_t tmp_61 = tmp_20*tmp_57;
      real_t tmp_62 = tmp_60 + tmp_61;
      real_t tmp_63 = tmp_1*(tmp_59 - 1.0/3.0) + tmp_8*(tmp_62 - 1.0/3.0);
      real_t tmp_64 = -tmp_56 - tmp_58 - tmp_60 - tmp_61 + 1;
      real_t tmp_65 = tmp_28*tmp_63;
      real_t tmp_66 = 0.2393143352496831*tmp_27;
      real_t tmp_67 = 0.95308992296933193*tmp_2 + tmp_4;
      real_t tmp_68 = tmp_10*tmp_67;
      real_t tmp_69 = 0.95308992296933193*tmp_12 + tmp_13;
      real_t tmp_70 = tmp_15*tmp_69;
      real_t tmp_71 = tmp_68 + tmp_70;
      real_t tmp_72 = tmp_18*tmp_67;
      real_t tmp_73 = tmp_20*tmp_69;
      real_t tmp_74 = tmp_72 + tmp_73;
      real_t tmp_75 = tmp_1*(tmp_71 - 1.0/3.0) + tmp_8*(tmp_74 - 1.0/3.0);
      real_t tmp_76 = -tmp_68 - tmp_70 - tmp_72 - tmp_73 + 1;
      real_t tmp_77 = tmp_28*tmp_75;
      real_t tmp_78 = 0.11846344252809471*tmp_27;
      real_t tmp_79 = 0.5*p_affine_10_0*tmp_15 + 0.5*p_affine_10_1*tmp_10;
      real_t tmp_80 = 0.5*p_affine_10_0*tmp_20 + 0.5*p_affine_10_1*tmp_18;
      real_t a_0_0 = tmp_30*(-tmp_23*tmp_24 - tmp_25*tmp_26 + tmp_25*tmp_29) + tmp_42*(-tmp_24*tmp_39 - tmp_26*tmp_40 + tmp_40*tmp_41) + tmp_54*(-tmp_24*tmp_51 - tmp_26*tmp_52 + tmp_52*tmp_53) + tmp_66*(-tmp_24*tmp_63 - tmp_26*tmp_64 + tmp_64*tmp_65) + tmp_78*(-tmp_24*tmp_75 - tmp_26*tmp_76 + tmp_76*tmp_77);
      real_t a_1_0 = tmp_30*(-tmp_17*tmp_26 + tmp_17*tmp_29 - tmp_23*tmp_79) + tmp_42*(-tmp_26*tmp_35 + tmp_35*tmp_41 - tmp_39*tmp_79) + tmp_54*(-tmp_26*tmp_47 + tmp_47*tmp_53 - tmp_51*tmp_79) + tmp_66*(-tmp_26*tmp_59 + tmp_59*tmp_65 - tmp_63*tmp_79) + tmp_78*(-tmp_26*tmp_71 + tmp_71*tmp_77 - tmp_75*tmp_79);
      real_t a_2_0 = tmp_30*(-tmp_22*tmp_26 + tmp_22*tmp_29 - tmp_23*tmp_80) + tmp_42*(-tmp_26*tmp_38 + tmp_38*tmp_41 - tmp_39*tmp_80) + tmp_54*(-tmp_26*tmp_50 + tmp_50*tmp_53 - tmp_51*tmp_80) + tmp_66*(-tmp_26*tmp_62 + tmp_62*tmp_65 - tmp_63*tmp_80) + tmp_78*(-tmp_26*tmp_74 + tmp_74*tmp_77 - tmp_75*tmp_80);
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

      real_t tmp_0 = -p_affine_3_0;
      real_t tmp_1 = p_affine_4_0 + tmp_0;
      real_t tmp_2 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_3 = -p_affine_3_1;
      real_t tmp_4 = p_affine_5_1 + tmp_3;
      real_t tmp_5 = tmp_1*tmp_4;
      real_t tmp_6 = p_affine_5_0 + tmp_0;
      real_t tmp_7 = 1.0 / (tmp_5 - tmp_6*(p_affine_4_1 + tmp_3));
      real_t tmp_8 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9 = p_affine_6_1 + 0.046910077030668018*tmp_8;
      real_t tmp_10 = tmp_7*(tmp_3 + tmp_9);
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.046910077030668018*tmp_11;
      real_t tmp_13 = tmp_7*(tmp_0 + tmp_12);
      real_t tmp_14 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_15 = tmp_1*(tmp_10*tmp_2 + tmp_13*tmp_4 - 1.0/3.0) + tmp_6*(tmp_1*tmp_10 + tmp_13*tmp_14 - 1.0/3.0);
      real_t tmp_16 = -p_affine_0_1;
      real_t tmp_17 = p_affine_2_1 + tmp_16;
      real_t tmp_18 = -p_affine_0_0;
      real_t tmp_19 = p_affine_1_0 + tmp_18;
      real_t tmp_20 = 1.0 / (tmp_17*tmp_19 - (p_affine_1_1 + tmp_16)*(p_affine_2_0 + tmp_18));
      real_t tmp_21 = tmp_17*tmp_20;
      real_t tmp_22 = tmp_20*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_23 = tmp_19*tmp_20;
      real_t tmp_24 = tmp_20*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_25 = 0.5*p_affine_10_0*(-tmp_21 - tmp_22) + 0.5*p_affine_10_1*(-tmp_23 - tmp_24);
      real_t tmp_26 = tmp_16 + tmp_9;
      real_t tmp_27 = tmp_23*tmp_26;
      real_t tmp_28 = tmp_24*tmp_26;
      real_t tmp_29 = tmp_12 + tmp_18;
      real_t tmp_30 = tmp_21*tmp_29;
      real_t tmp_31 = tmp_22*tmp_29;
      real_t tmp_32 = -tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1;
      real_t tmp_33 = tmp_6*tmp_7;
      real_t tmp_34 = tmp_2*tmp_7;
      real_t tmp_35 = 0.5*p_affine_10_0*(tmp_14*tmp_33 + tmp_5*tmp_7) + 0.5*p_affine_10_1*(tmp_1*tmp_33 + tmp_1*tmp_34);
      real_t tmp_36 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_8*tmp_8), 1.0/2.0));
      real_t tmp_37 = 6/tmp_36;
      real_t tmp_38 = tmp_15*tmp_37;
      real_t tmp_39 = 0.11846344252809471*tmp_36;
      real_t tmp_40 = p_affine_6_1 + 0.23076534494715845*tmp_8;
      real_t tmp_41 = tmp_3 + tmp_40;
      real_t tmp_42 = p_affine_6_0 + 0.23076534494715845*tmp_11;
      real_t tmp_43 = tmp_7*(tmp_0 + tmp_42);
      real_t tmp_44 = tmp_1*tmp_7;
      real_t tmp_45 = tmp_1*(tmp_34*tmp_41 + tmp_4*tmp_43 - 1.0/3.0) + tmp_6*(tmp_14*tmp_43 + tmp_41*tmp_44 - 1.0/3.0);
      real_t tmp_46 = tmp_16 + tmp_40;
      real_t tmp_47 = tmp_23*tmp_46;
      real_t tmp_48 = tmp_24*tmp_46;
      real_t tmp_49 = tmp_18 + tmp_42;
      real_t tmp_50 = tmp_21*tmp_49;
      real_t tmp_51 = tmp_22*tmp_49;
      real_t tmp_52 = -tmp_47 - tmp_48 - tmp_50 - tmp_51 + 1;
      real_t tmp_53 = tmp_37*tmp_45;
      real_t tmp_54 = 0.2393143352496831*tmp_36;
      real_t tmp_55 = p_affine_6_1 + 0.5*tmp_8;
      real_t tmp_56 = tmp_3 + tmp_55;
      real_t tmp_57 = p_affine_6_0 + 0.5*tmp_11;
      real_t tmp_58 = tmp_7*(tmp_0 + tmp_57);
      real_t tmp_59 = tmp_1*(tmp_34*tmp_56 + tmp_4*tmp_58 - 1.0/3.0) + tmp_6*(tmp_14*tmp_58 + tmp_44*tmp_56 - 1.0/3.0);
      real_t tmp_60 = tmp_16 + tmp_55;
      real_t tmp_61 = tmp_23*tmp_60;
      real_t tmp_62 = tmp_24*tmp_60;
      real_t tmp_63 = tmp_18 + tmp_57;
      real_t tmp_64 = tmp_21*tmp_63;
      real_t tmp_65 = tmp_22*tmp_63;
      real_t tmp_66 = -tmp_61 - tmp_62 - tmp_64 - tmp_65 + 1;
      real_t tmp_67 = tmp_37*tmp_59;
      real_t tmp_68 = 0.2844444444444445*tmp_36;
      real_t tmp_69 = p_affine_6_1 + 0.7692346550528415*tmp_8;
      real_t tmp_70 = tmp_3 + tmp_69;
      real_t tmp_71 = p_affine_6_0 + 0.7692346550528415*tmp_11;
      real_t tmp_72 = tmp_7*(tmp_0 + tmp_71);
      real_t tmp_73 = tmp_1*(tmp_34*tmp_70 + tmp_4*tmp_72 - 1.0/3.0) + tmp_6*(tmp_14*tmp_72 + tmp_44*tmp_70 - 1.0/3.0);
      real_t tmp_74 = tmp_16 + tmp_69;
      real_t tmp_75 = tmp_23*tmp_74;
      real_t tmp_76 = tmp_24*tmp_74;
      real_t tmp_77 = tmp_18 + tmp_71;
      real_t tmp_78 = tmp_21*tmp_77;
      real_t tmp_79 = tmp_22*tmp_77;
      real_t tmp_80 = -tmp_75 - tmp_76 - tmp_78 - tmp_79 + 1;
      real_t tmp_81 = tmp_37*tmp_73;
      real_t tmp_82 = 0.2393143352496831*tmp_36;
      real_t tmp_83 = p_affine_6_1 + 0.95308992296933193*tmp_8;
      real_t tmp_84 = tmp_3 + tmp_83;
      real_t tmp_85 = p_affine_6_0 + 0.95308992296933193*tmp_11;
      real_t tmp_86 = tmp_7*(tmp_0 + tmp_85);
      real_t tmp_87 = tmp_1*(tmp_34*tmp_84 + tmp_4*tmp_86 - 1.0/3.0) + tmp_6*(tmp_14*tmp_86 + tmp_44*tmp_84 - 1.0/3.0);
      real_t tmp_88 = tmp_16 + tmp_83;
      real_t tmp_89 = tmp_23*tmp_88;
      real_t tmp_90 = tmp_24*tmp_88;
      real_t tmp_91 = tmp_18 + tmp_85;
      real_t tmp_92 = tmp_21*tmp_91;
      real_t tmp_93 = tmp_22*tmp_91;
      real_t tmp_94 = -tmp_89 - tmp_90 - tmp_92 - tmp_93 + 1;
      real_t tmp_95 = tmp_37*tmp_87;
      real_t tmp_96 = 0.11846344252809471*tmp_36;
      real_t tmp_97 = tmp_28 + tmp_30;
      real_t tmp_98 = 0.5*p_affine_10_0*tmp_21 + 0.5*p_affine_10_1*tmp_24;
      real_t tmp_99 = tmp_48 + tmp_50;
      real_t tmp_100 = tmp_62 + tmp_64;
      real_t tmp_101 = tmp_76 + tmp_78;
      real_t tmp_102 = tmp_90 + tmp_92;
      real_t tmp_103 = tmp_27 + tmp_31;
      real_t tmp_104 = 0.5*p_affine_10_0*tmp_22 + 0.5*p_affine_10_1*tmp_23;
      real_t tmp_105 = tmp_47 + tmp_51;
      real_t tmp_106 = tmp_61 + tmp_65;
      real_t tmp_107 = tmp_75 + tmp_79;
      real_t tmp_108 = tmp_89 + tmp_93;
      real_t a_0_0 = tmp_39*(tmp_15*tmp_25 - tmp_32*tmp_35 - tmp_32*tmp_38) + tmp_54*(tmp_25*tmp_45 - tmp_35*tmp_52 - tmp_52*tmp_53) + tmp_68*(tmp_25*tmp_59 - tmp_35*tmp_66 - tmp_66*tmp_67) + tmp_82*(tmp_25*tmp_73 - tmp_35*tmp_80 - tmp_80*tmp_81) + tmp_96*(tmp_25*tmp_87 - tmp_35*tmp_94 - tmp_94*tmp_95);
      real_t a_1_0 = tmp_39*(tmp_15*tmp_98 - tmp_35*tmp_97 - tmp_38*tmp_97) + tmp_54*(-tmp_35*tmp_99 + tmp_45*tmp_98 - tmp_53*tmp_99) + tmp_68*(-tmp_100*tmp_35 - tmp_100*tmp_67 + tmp_59*tmp_98) + tmp_82*(-tmp_101*tmp_35 - tmp_101*tmp_81 + tmp_73*tmp_98) + tmp_96*(-tmp_102*tmp_35 - tmp_102*tmp_95 + tmp_87*tmp_98);
      real_t a_2_0 = tmp_39*(-tmp_103*tmp_35 - tmp_103*tmp_38 + tmp_104*tmp_15) + tmp_54*(tmp_104*tmp_45 - tmp_105*tmp_35 - tmp_105*tmp_53) + tmp_68*(tmp_104*tmp_59 - tmp_106*tmp_35 - tmp_106*tmp_67) + tmp_82*(tmp_104*tmp_73 - tmp_107*tmp_35 - tmp_107*tmp_81) + tmp_96*(tmp_104*tmp_87 - tmp_108*tmp_35 - tmp_108*tmp_95);
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

      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_0;
      real_t tmp_3 = p_affine_1_0 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_2;
      real_t tmp_6 = 1.0 / (tmp_4 - tmp_5*(p_affine_1_1 + tmp_0));
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = tmp_6*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_9 = tmp_3*tmp_6;
      real_t tmp_10 = tmp_6*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_11 = p_affine_10_0*(-tmp_7 - tmp_8) + p_affine_10_1*(-tmp_10 - tmp_9);
      real_t tmp_12 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_13 = p_affine_6_1 + tmp_0;
      real_t tmp_14 = 0.046910077030668018*tmp_12 + tmp_13;
      real_t tmp_15 = tmp_10*tmp_14;
      real_t tmp_16 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_17 = p_affine_6_0 + tmp_2;
      real_t tmp_18 = 0.046910077030668018*tmp_16 + tmp_17;
      real_t tmp_19 = tmp_18*tmp_7;
      real_t tmp_20 = tmp_15 + tmp_19;
      real_t tmp_21 = tmp_14*tmp_9;
      real_t tmp_22 = tmp_18*tmp_8;
      real_t tmp_23 = tmp_21 + tmp_22;
      real_t tmp_24 = tmp_3*(tmp_20 - 1.0/3.0) + tmp_5*(tmp_23 - 1.0/3.0);
      real_t tmp_25 = p_affine_10_0*(tmp_4*tmp_6 + tmp_5*tmp_8) + p_affine_10_1*(tmp_10*tmp_3 + tmp_5*tmp_9);
      real_t tmp_26 = -tmp_15 - tmp_19 - tmp_21 - tmp_22 + 1;
      real_t tmp_27 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_16*tmp_16), 1.0/2.0));
      real_t tmp_28 = 6/tmp_27;
      real_t tmp_29 = tmp_24*tmp_28;
      real_t tmp_30 = 0.11846344252809471*tmp_27;
      real_t tmp_31 = 0.23076534494715845*tmp_12 + tmp_13;
      real_t tmp_32 = tmp_10*tmp_31;
      real_t tmp_33 = 0.23076534494715845*tmp_16 + tmp_17;
      real_t tmp_34 = tmp_33*tmp_7;
      real_t tmp_35 = tmp_32 + tmp_34;
      real_t tmp_36 = tmp_31*tmp_9;
      real_t tmp_37 = tmp_33*tmp_8;
      real_t tmp_38 = tmp_36 + tmp_37;
      real_t tmp_39 = tmp_3*(tmp_35 - 1.0/3.0) + tmp_5*(tmp_38 - 1.0/3.0);
      real_t tmp_40 = -tmp_32 - tmp_34 - tmp_36 - tmp_37 + 1;
      real_t tmp_41 = tmp_28*tmp_39;
      real_t tmp_42 = 0.2393143352496831*tmp_27;
      real_t tmp_43 = 0.5*tmp_12 + tmp_13;
      real_t tmp_44 = tmp_10*tmp_43;
      real_t tmp_45 = 0.5*tmp_16 + tmp_17;
      real_t tmp_46 = tmp_45*tmp_7;
      real_t tmp_47 = tmp_44 + tmp_46;
      real_t tmp_48 = tmp_43*tmp_9;
      real_t tmp_49 = tmp_45*tmp_8;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_3*(tmp_47 - 1.0/3.0) + tmp_5*(tmp_50 - 1.0/3.0);
      real_t tmp_52 = -tmp_44 - tmp_46 - tmp_48 - tmp_49 + 1;
      real_t tmp_53 = tmp_28*tmp_51;
      real_t tmp_54 = 0.2844444444444445*tmp_27;
      real_t tmp_55 = 0.7692346550528415*tmp_12 + tmp_13;
      real_t tmp_56 = tmp_10*tmp_55;
      real_t tmp_57 = 0.7692346550528415*tmp_16 + tmp_17;
      real_t tmp_58 = tmp_57*tmp_7;
      real_t tmp_59 = tmp_56 + tmp_58;
      real_t tmp_60 = tmp_55*tmp_9;
      real_t tmp_61 = tmp_57*tmp_8;
      real_t tmp_62 = tmp_60 + tmp_61;
      real_t tmp_63 = tmp_3*(tmp_59 - 1.0/3.0) + tmp_5*(tmp_62 - 1.0/3.0);
      real_t tmp_64 = -tmp_56 - tmp_58 - tmp_60 - tmp_61 + 1;
      real_t tmp_65 = tmp_28*tmp_63;
      real_t tmp_66 = 0.2393143352496831*tmp_27;
      real_t tmp_67 = 0.95308992296933193*tmp_12 + tmp_13;
      real_t tmp_68 = tmp_10*tmp_67;
      real_t tmp_69 = 0.95308992296933193*tmp_16 + tmp_17;
      real_t tmp_70 = tmp_69*tmp_7;
      real_t tmp_71 = tmp_68 + tmp_70;
      real_t tmp_72 = tmp_67*tmp_9;
      real_t tmp_73 = tmp_69*tmp_8;
      real_t tmp_74 = tmp_72 + tmp_73;
      real_t tmp_75 = tmp_3*(tmp_71 - 1.0/3.0) + tmp_5*(tmp_74 - 1.0/3.0);
      real_t tmp_76 = -tmp_68 - tmp_70 - tmp_72 - tmp_73 + 1;
      real_t tmp_77 = tmp_28*tmp_75;
      real_t tmp_78 = 0.11846344252809471*tmp_27;
      real_t tmp_79 = p_affine_10_0*tmp_7 + p_affine_10_1*tmp_10;
      real_t tmp_80 = p_affine_10_0*tmp_8 + p_affine_10_1*tmp_9;
      real_t a_0_0 = tmp_30*(-tmp_11*tmp_24 - tmp_25*tmp_26 + tmp_26*tmp_29) + tmp_42*(-tmp_11*tmp_39 - tmp_25*tmp_40 + tmp_40*tmp_41) + tmp_54*(-tmp_11*tmp_51 - tmp_25*tmp_52 + tmp_52*tmp_53) + tmp_66*(-tmp_11*tmp_63 - tmp_25*tmp_64 + tmp_64*tmp_65) + tmp_78*(-tmp_11*tmp_75 - tmp_25*tmp_76 + tmp_76*tmp_77);
      real_t a_1_0 = tmp_30*(-tmp_20*tmp_25 + tmp_20*tmp_29 - tmp_24*tmp_79) + tmp_42*(-tmp_25*tmp_35 + tmp_35*tmp_41 - tmp_39*tmp_79) + tmp_54*(-tmp_25*tmp_47 + tmp_47*tmp_53 - tmp_51*tmp_79) + tmp_66*(-tmp_25*tmp_59 + tmp_59*tmp_65 - tmp_63*tmp_79) + tmp_78*(-tmp_25*tmp_71 + tmp_71*tmp_77 - tmp_75*tmp_79);
      real_t a_2_0 = tmp_30*(-tmp_23*tmp_25 + tmp_23*tmp_29 - tmp_24*tmp_80) + tmp_42*(-tmp_25*tmp_38 + tmp_38*tmp_41 - tmp_39*tmp_80) + tmp_54*(-tmp_25*tmp_50 + tmp_50*tmp_53 - tmp_51*tmp_80) + tmp_66*(-tmp_25*tmp_62 + tmp_62*tmp_65 - tmp_63*tmp_80) + tmp_78*(-tmp_25*tmp_74 + tmp_74*tmp_77 - tmp_75*tmp_80);
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




class EGVectorLaplaceFormP1EDG_new_1 : public hyteg::dg::DGForm2D
{
 public:
    EGVectorLaplaceFormP1EDG_new_1()

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
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = 1.0 / (tmp_4 - tmp_5*(p_affine_2_0 + tmp_0));
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = tmp_6*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_9 = tmp_4*tmp_6 + tmp_5*tmp_8;
      real_t tmp_10 = tmp_3*tmp_6;
      real_t tmp_11 = tmp_6*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_12 = tmp_10*tmp_5 + tmp_11*tmp_3;
      real_t tmp_13 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_14 = tmp_13*(tmp_12*(-tmp_10 - tmp_11) + tmp_9*(-tmp_7 - tmp_8));
      real_t tmp_15 = tmp_13*(tmp_10*tmp_12 + tmp_8*tmp_9);
      real_t tmp_16 = tmp_13*(tmp_11*tmp_12 + tmp_7*tmp_9);
      real_t a_0_0 = 0.5*tmp_14;
      real_t a_1_0 = 0.5*tmp_15;
      real_t a_2_0 = 0.5*tmp_16;
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

      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_1_1 + tmp_0;
      real_t tmp_2 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_3 = p_affine_6_1 + tmp_0;
      real_t tmp_4 = 0.046910077030668018*tmp_2 + tmp_3;
      real_t tmp_5 = -p_affine_0_0;
      real_t tmp_6 = p_affine_1_0 + tmp_5;
      real_t tmp_7 = p_affine_2_1 + tmp_0;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = 1.0 / (-tmp_1*(p_affine_2_0 + tmp_5) + tmp_8);
      real_t tmp_10 = tmp_9*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_11 = tmp_10*tmp_4;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_5;
      real_t tmp_14 = 0.046910077030668018*tmp_12 + tmp_13;
      real_t tmp_15 = tmp_7*tmp_9;
      real_t tmp_16 = tmp_14*tmp_15;
      real_t tmp_17 = tmp_11 + tmp_16;
      real_t tmp_18 = tmp_6*tmp_9;
      real_t tmp_19 = tmp_18*tmp_4;
      real_t tmp_20 = tmp_9*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_21 = tmp_14*tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_1*(tmp_17 - 1.0/3.0) + tmp_7*(tmp_22 - 1.0/3.0);
      real_t tmp_24 = 0.5*p_affine_10_0*(-tmp_15 - tmp_20) + 0.5*p_affine_10_1*(-tmp_10 - tmp_18);
      real_t tmp_25 = -tmp_11 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_26 = 0.5*p_affine_10_0*(tmp_1*tmp_15 + tmp_20*tmp_7) + 0.5*p_affine_10_1*(tmp_1*tmp_10 + tmp_8*tmp_9);
      real_t tmp_27 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_2*tmp_2), 1.0/2.0));
      real_t tmp_28 = 6/tmp_27;
      real_t tmp_29 = tmp_23*tmp_28;
      real_t tmp_30 = 0.11846344252809471*tmp_27;
      real_t tmp_31 = 0.23076534494715845*tmp_2 + tmp_3;
      real_t tmp_32 = tmp_10*tmp_31;
      real_t tmp_33 = 0.23076534494715845*tmp_12 + tmp_13;
      real_t tmp_34 = tmp_15*tmp_33;
      real_t tmp_35 = tmp_32 + tmp_34;
      real_t tmp_36 = tmp_18*tmp_31;
      real_t tmp_37 = tmp_20*tmp_33;
      real_t tmp_38 = tmp_36 + tmp_37;
      real_t tmp_39 = tmp_1*(tmp_35 - 1.0/3.0) + tmp_7*(tmp_38 - 1.0/3.0);
      real_t tmp_40 = -tmp_32 - tmp_34 - tmp_36 - tmp_37 + 1;
      real_t tmp_41 = tmp_28*tmp_39;
      real_t tmp_42 = 0.2393143352496831*tmp_27;
      real_t tmp_43 = 0.5*tmp_2 + tmp_3;
      real_t tmp_44 = tmp_10*tmp_43;
      real_t tmp_45 = 0.5*tmp_12 + tmp_13;
      real_t tmp_46 = tmp_15*tmp_45;
      real_t tmp_47 = tmp_44 + tmp_46;
      real_t tmp_48 = tmp_18*tmp_43;
      real_t tmp_49 = tmp_20*tmp_45;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_1*(tmp_47 - 1.0/3.0) + tmp_7*(tmp_50 - 1.0/3.0);
      real_t tmp_52 = -tmp_44 - tmp_46 - tmp_48 - tmp_49 + 1;
      real_t tmp_53 = tmp_28*tmp_51;
      real_t tmp_54 = 0.2844444444444445*tmp_27;
      real_t tmp_55 = 0.7692346550528415*tmp_2 + tmp_3;
      real_t tmp_56 = tmp_10*tmp_55;
      real_t tmp_57 = 0.7692346550528415*tmp_12 + tmp_13;
      real_t tmp_58 = tmp_15*tmp_57;
      real_t tmp_59 = tmp_56 + tmp_58;
      real_t tmp_60 = tmp_18*tmp_55;
      real_t tmp_61 = tmp_20*tmp_57;
      real_t tmp_62 = tmp_60 + tmp_61;
      real_t tmp_63 = tmp_1*(tmp_59 - 1.0/3.0) + tmp_7*(tmp_62 - 1.0/3.0);
      real_t tmp_64 = -tmp_56 - tmp_58 - tmp_60 - tmp_61 + 1;
      real_t tmp_65 = tmp_28*tmp_63;
      real_t tmp_66 = 0.2393143352496831*tmp_27;
      real_t tmp_67 = 0.95308992296933193*tmp_2 + tmp_3;
      real_t tmp_68 = tmp_10*tmp_67;
      real_t tmp_69 = 0.95308992296933193*tmp_12 + tmp_13;
      real_t tmp_70 = tmp_15*tmp_69;
      real_t tmp_71 = tmp_68 + tmp_70;
      real_t tmp_72 = tmp_18*tmp_67;
      real_t tmp_73 = tmp_20*tmp_69;
      real_t tmp_74 = tmp_72 + tmp_73;
      real_t tmp_75 = tmp_1*(tmp_71 - 1.0/3.0) + tmp_7*(tmp_74 - 1.0/3.0);
      real_t tmp_76 = -tmp_68 - tmp_70 - tmp_72 - tmp_73 + 1;
      real_t tmp_77 = tmp_28*tmp_75;
      real_t tmp_78 = 0.11846344252809471*tmp_27;
      real_t tmp_79 = 0.5*p_affine_10_0*tmp_15 + 0.5*p_affine_10_1*tmp_10;
      real_t tmp_80 = 0.5*p_affine_10_0*tmp_20 + 0.5*p_affine_10_1*tmp_18;
      real_t a_0_0 = tmp_30*(-tmp_23*tmp_24 - tmp_25*tmp_26 + tmp_25*tmp_29) + tmp_42*(-tmp_24*tmp_39 - tmp_26*tmp_40 + tmp_40*tmp_41) + tmp_54*(-tmp_24*tmp_51 - tmp_26*tmp_52 + tmp_52*tmp_53) + tmp_66*(-tmp_24*tmp_63 - tmp_26*tmp_64 + tmp_64*tmp_65) + tmp_78*(-tmp_24*tmp_75 - tmp_26*tmp_76 + tmp_76*tmp_77);
      real_t a_1_0 = tmp_30*(-tmp_17*tmp_26 + tmp_17*tmp_29 - tmp_23*tmp_79) + tmp_42*(-tmp_26*tmp_35 + tmp_35*tmp_41 - tmp_39*tmp_79) + tmp_54*(-tmp_26*tmp_47 + tmp_47*tmp_53 - tmp_51*tmp_79) + tmp_66*(-tmp_26*tmp_59 + tmp_59*tmp_65 - tmp_63*tmp_79) + tmp_78*(-tmp_26*tmp_71 + tmp_71*tmp_77 - tmp_75*tmp_79);
      real_t a_2_0 = tmp_30*(-tmp_22*tmp_26 + tmp_22*tmp_29 - tmp_23*tmp_80) + tmp_42*(-tmp_26*tmp_38 + tmp_38*tmp_41 - tmp_39*tmp_80) + tmp_54*(-tmp_26*tmp_50 + tmp_50*tmp_53 - tmp_51*tmp_80) + tmp_66*(-tmp_26*tmp_62 + tmp_62*tmp_65 - tmp_63*tmp_80) + tmp_78*(-tmp_26*tmp_74 + tmp_74*tmp_77 - tmp_75*tmp_80);
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

      real_t tmp_0 = -p_affine_3_1;
      real_t tmp_1 = p_affine_4_1 + tmp_0;
      real_t tmp_2 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_3 = -p_affine_3_0;
      real_t tmp_4 = p_affine_4_0 + tmp_3;
      real_t tmp_5 = p_affine_5_1 + tmp_0;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = 1.0 / (-tmp_1*(p_affine_5_0 + tmp_3) + tmp_6);
      real_t tmp_8 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9 = p_affine_6_1 + 0.046910077030668018*tmp_8;
      real_t tmp_10 = tmp_7*(tmp_0 + tmp_9);
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.046910077030668018*tmp_11;
      real_t tmp_13 = tmp_7*(tmp_12 + tmp_3);
      real_t tmp_14 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_15 = tmp_1*(tmp_10*tmp_2 + tmp_13*tmp_5 - 1.0/3.0) + tmp_5*(tmp_10*tmp_4 + tmp_13*tmp_14 - 1.0/3.0);
      real_t tmp_16 = -p_affine_0_1;
      real_t tmp_17 = p_affine_2_1 + tmp_16;
      real_t tmp_18 = -p_affine_0_0;
      real_t tmp_19 = p_affine_1_0 + tmp_18;
      real_t tmp_20 = 1.0 / (tmp_17*tmp_19 - (p_affine_1_1 + tmp_16)*(p_affine_2_0 + tmp_18));
      real_t tmp_21 = tmp_17*tmp_20;
      real_t tmp_22 = tmp_20*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_23 = tmp_19*tmp_20;
      real_t tmp_24 = tmp_20*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_25 = 0.5*p_affine_10_0*(-tmp_21 - tmp_22) + 0.5*p_affine_10_1*(-tmp_23 - tmp_24);
      real_t tmp_26 = tmp_16 + tmp_9;
      real_t tmp_27 = tmp_23*tmp_26;
      real_t tmp_28 = tmp_24*tmp_26;
      real_t tmp_29 = tmp_12 + tmp_18;
      real_t tmp_30 = tmp_21*tmp_29;
      real_t tmp_31 = tmp_22*tmp_29;
      real_t tmp_32 = -tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1;
      real_t tmp_33 = tmp_5*tmp_7;
      real_t tmp_34 = tmp_2*tmp_7;
      real_t tmp_35 = 0.5*p_affine_10_0*(tmp_1*tmp_33 + tmp_14*tmp_33) + 0.5*p_affine_10_1*(tmp_1*tmp_34 + tmp_6*tmp_7);
      real_t tmp_36 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_8*tmp_8), 1.0/2.0));
      real_t tmp_37 = 6/tmp_36;
      real_t tmp_38 = tmp_15*tmp_37;
      real_t tmp_39 = 0.11846344252809471*tmp_36;
      real_t tmp_40 = p_affine_6_1 + 0.23076534494715845*tmp_8;
      real_t tmp_41 = tmp_0 + tmp_40;
      real_t tmp_42 = p_affine_6_0 + 0.23076534494715845*tmp_11;
      real_t tmp_43 = tmp_3 + tmp_42;
      real_t tmp_44 = tmp_4*tmp_7;
      real_t tmp_45 = tmp_14*tmp_7;
      real_t tmp_46 = tmp_1*(tmp_33*tmp_43 + tmp_34*tmp_41 - 1.0/3.0) + tmp_5*(tmp_41*tmp_44 + tmp_43*tmp_45 - 1.0/3.0);
      real_t tmp_47 = tmp_16 + tmp_40;
      real_t tmp_48 = tmp_23*tmp_47;
      real_t tmp_49 = tmp_24*tmp_47;
      real_t tmp_50 = tmp_18 + tmp_42;
      real_t tmp_51 = tmp_21*tmp_50;
      real_t tmp_52 = tmp_22*tmp_50;
      real_t tmp_53 = -tmp_48 - tmp_49 - tmp_51 - tmp_52 + 1;
      real_t tmp_54 = tmp_37*tmp_46;
      real_t tmp_55 = 0.2393143352496831*tmp_36;
      real_t tmp_56 = p_affine_6_1 + 0.5*tmp_8;
      real_t tmp_57 = tmp_0 + tmp_56;
      real_t tmp_58 = p_affine_6_0 + 0.5*tmp_11;
      real_t tmp_59 = tmp_3 + tmp_58;
      real_t tmp_60 = tmp_1*(tmp_33*tmp_59 + tmp_34*tmp_57 - 1.0/3.0) + tmp_5*(tmp_44*tmp_57 + tmp_45*tmp_59 - 1.0/3.0);
      real_t tmp_61 = tmp_16 + tmp_56;
      real_t tmp_62 = tmp_23*tmp_61;
      real_t tmp_63 = tmp_24*tmp_61;
      real_t tmp_64 = tmp_18 + tmp_58;
      real_t tmp_65 = tmp_21*tmp_64;
      real_t tmp_66 = tmp_22*tmp_64;
      real_t tmp_67 = -tmp_62 - tmp_63 - tmp_65 - tmp_66 + 1;
      real_t tmp_68 = tmp_37*tmp_60;
      real_t tmp_69 = 0.2844444444444445*tmp_36;
      real_t tmp_70 = p_affine_6_1 + 0.7692346550528415*tmp_8;
      real_t tmp_71 = tmp_0 + tmp_70;
      real_t tmp_72 = p_affine_6_0 + 0.7692346550528415*tmp_11;
      real_t tmp_73 = tmp_3 + tmp_72;
      real_t tmp_74 = tmp_1*(tmp_33*tmp_73 + tmp_34*tmp_71 - 1.0/3.0) + tmp_5*(tmp_44*tmp_71 + tmp_45*tmp_73 - 1.0/3.0);
      real_t tmp_75 = tmp_16 + tmp_70;
      real_t tmp_76 = tmp_23*tmp_75;
      real_t tmp_77 = tmp_24*tmp_75;
      real_t tmp_78 = tmp_18 + tmp_72;
      real_t tmp_79 = tmp_21*tmp_78;
      real_t tmp_80 = tmp_22*tmp_78;
      real_t tmp_81 = -tmp_76 - tmp_77 - tmp_79 - tmp_80 + 1;
      real_t tmp_82 = tmp_37*tmp_74;
      real_t tmp_83 = 0.2393143352496831*tmp_36;
      real_t tmp_84 = p_affine_6_1 + 0.95308992296933193*tmp_8;
      real_t tmp_85 = tmp_0 + tmp_84;
      real_t tmp_86 = p_affine_6_0 + 0.95308992296933193*tmp_11;
      real_t tmp_87 = tmp_3 + tmp_86;
      real_t tmp_88 = tmp_1*(tmp_33*tmp_87 + tmp_34*tmp_85 - 1.0/3.0) + tmp_5*(tmp_44*tmp_85 + tmp_45*tmp_87 - 1.0/3.0);
      real_t tmp_89 = tmp_16 + tmp_84;
      real_t tmp_90 = tmp_23*tmp_89;
      real_t tmp_91 = tmp_24*tmp_89;
      real_t tmp_92 = tmp_18 + tmp_86;
      real_t tmp_93 = tmp_21*tmp_92;
      real_t tmp_94 = tmp_22*tmp_92;
      real_t tmp_95 = -tmp_90 - tmp_91 - tmp_93 - tmp_94 + 1;
      real_t tmp_96 = tmp_37*tmp_88;
      real_t tmp_97 = 0.11846344252809471*tmp_36;
      real_t tmp_98 = tmp_28 + tmp_30;
      real_t tmp_99 = 0.5*p_affine_10_0*tmp_21 + 0.5*p_affine_10_1*tmp_24;
      real_t tmp_100 = tmp_49 + tmp_51;
      real_t tmp_101 = tmp_63 + tmp_65;
      real_t tmp_102 = tmp_77 + tmp_79;
      real_t tmp_103 = tmp_91 + tmp_93;
      real_t tmp_104 = tmp_27 + tmp_31;
      real_t tmp_105 = 0.5*p_affine_10_0*tmp_22 + 0.5*p_affine_10_1*tmp_23;
      real_t tmp_106 = tmp_48 + tmp_52;
      real_t tmp_107 = tmp_62 + tmp_66;
      real_t tmp_108 = tmp_76 + tmp_80;
      real_t tmp_109 = tmp_90 + tmp_94;
      real_t a_0_0 = tmp_39*(tmp_15*tmp_25 - tmp_32*tmp_35 - tmp_32*tmp_38) + tmp_55*(tmp_25*tmp_46 - tmp_35*tmp_53 - tmp_53*tmp_54) + tmp_69*(tmp_25*tmp_60 - tmp_35*tmp_67 - tmp_67*tmp_68) + tmp_83*(tmp_25*tmp_74 - tmp_35*tmp_81 - tmp_81*tmp_82) + tmp_97*(tmp_25*tmp_88 - tmp_35*tmp_95 - tmp_95*tmp_96);
      real_t a_1_0 = tmp_39*(tmp_15*tmp_99 - tmp_35*tmp_98 - tmp_38*tmp_98) + tmp_55*(-tmp_100*tmp_35 - tmp_100*tmp_54 + tmp_46*tmp_99) + tmp_69*(-tmp_101*tmp_35 - tmp_101*tmp_68 + tmp_60*tmp_99) + tmp_83*(-tmp_102*tmp_35 - tmp_102*tmp_82 + tmp_74*tmp_99) + tmp_97*(-tmp_103*tmp_35 - tmp_103*tmp_96 + tmp_88*tmp_99);
      real_t a_2_0 = tmp_39*(-tmp_104*tmp_35 - tmp_104*tmp_38 + tmp_105*tmp_15) + tmp_55*(tmp_105*tmp_46 - tmp_106*tmp_35 - tmp_106*tmp_54) + tmp_69*(tmp_105*tmp_60 - tmp_107*tmp_35 - tmp_107*tmp_68) + tmp_83*(tmp_105*tmp_74 - tmp_108*tmp_35 - tmp_108*tmp_82) + tmp_97*(tmp_105*tmp_88 - tmp_109*tmp_35 - tmp_109*tmp_96);
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

      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_0;
      real_t tmp_3 = p_affine_1_0 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_1_1 + tmp_0;
      real_t tmp_6 = 1.0 / (tmp_4 - tmp_5*(p_affine_2_0 + tmp_2));
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = tmp_6*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_9 = tmp_3*tmp_6;
      real_t tmp_10 = tmp_6*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_11 = p_affine_10_0*(-tmp_7 - tmp_8) + p_affine_10_1*(-tmp_10 - tmp_9);
      real_t tmp_12 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_13 = p_affine_6_1 + tmp_0;
      real_t tmp_14 = 0.046910077030668018*tmp_12 + tmp_13;
      real_t tmp_15 = tmp_10*tmp_14;
      real_t tmp_16 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_17 = p_affine_6_0 + tmp_2;
      real_t tmp_18 = 0.046910077030668018*tmp_16 + tmp_17;
      real_t tmp_19 = tmp_18*tmp_7;
      real_t tmp_20 = tmp_15 + tmp_19;
      real_t tmp_21 = tmp_14*tmp_9;
      real_t tmp_22 = tmp_18*tmp_8;
      real_t tmp_23 = tmp_21 + tmp_22;
      real_t tmp_24 = tmp_1*(tmp_23 - 1.0/3.0) + tmp_5*(tmp_20 - 1.0/3.0);
      real_t tmp_25 = p_affine_10_0*(tmp_1*tmp_8 + tmp_5*tmp_7) + p_affine_10_1*(tmp_10*tmp_5 + tmp_4*tmp_6);
      real_t tmp_26 = -tmp_15 - tmp_19 - tmp_21 - tmp_22 + 1;
      real_t tmp_27 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_16*tmp_16), 1.0/2.0));
      real_t tmp_28 = 6/tmp_27;
      real_t tmp_29 = tmp_24*tmp_28;
      real_t tmp_30 = 0.11846344252809471*tmp_27;
      real_t tmp_31 = 0.23076534494715845*tmp_12 + tmp_13;
      real_t tmp_32 = tmp_10*tmp_31;
      real_t tmp_33 = 0.23076534494715845*tmp_16 + tmp_17;
      real_t tmp_34 = tmp_33*tmp_7;
      real_t tmp_35 = tmp_32 + tmp_34;
      real_t tmp_36 = tmp_31*tmp_9;
      real_t tmp_37 = tmp_33*tmp_8;
      real_t tmp_38 = tmp_36 + tmp_37;
      real_t tmp_39 = tmp_1*(tmp_38 - 1.0/3.0) + tmp_5*(tmp_35 - 1.0/3.0);
      real_t tmp_40 = -tmp_32 - tmp_34 - tmp_36 - tmp_37 + 1;
      real_t tmp_41 = tmp_28*tmp_39;
      real_t tmp_42 = 0.2393143352496831*tmp_27;
      real_t tmp_43 = 0.5*tmp_12 + tmp_13;
      real_t tmp_44 = tmp_10*tmp_43;
      real_t tmp_45 = 0.5*tmp_16 + tmp_17;
      real_t tmp_46 = tmp_45*tmp_7;
      real_t tmp_47 = tmp_44 + tmp_46;
      real_t tmp_48 = tmp_43*tmp_9;
      real_t tmp_49 = tmp_45*tmp_8;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_1*(tmp_50 - 1.0/3.0) + tmp_5*(tmp_47 - 1.0/3.0);
      real_t tmp_52 = -tmp_44 - tmp_46 - tmp_48 - tmp_49 + 1;
      real_t tmp_53 = tmp_28*tmp_51;
      real_t tmp_54 = 0.2844444444444445*tmp_27;
      real_t tmp_55 = 0.7692346550528415*tmp_12 + tmp_13;
      real_t tmp_56 = tmp_10*tmp_55;
      real_t tmp_57 = 0.7692346550528415*tmp_16 + tmp_17;
      real_t tmp_58 = tmp_57*tmp_7;
      real_t tmp_59 = tmp_56 + tmp_58;
      real_t tmp_60 = tmp_55*tmp_9;
      real_t tmp_61 = tmp_57*tmp_8;
      real_t tmp_62 = tmp_60 + tmp_61;
      real_t tmp_63 = tmp_1*(tmp_62 - 1.0/3.0) + tmp_5*(tmp_59 - 1.0/3.0);
      real_t tmp_64 = -tmp_56 - tmp_58 - tmp_60 - tmp_61 + 1;
      real_t tmp_65 = tmp_28*tmp_63;
      real_t tmp_66 = 0.2393143352496831*tmp_27;
      real_t tmp_67 = 0.95308992296933193*tmp_12 + tmp_13;
      real_t tmp_68 = tmp_10*tmp_67;
      real_t tmp_69 = 0.95308992296933193*tmp_16 + tmp_17;
      real_t tmp_70 = tmp_69*tmp_7;
      real_t tmp_71 = tmp_68 + tmp_70;
      real_t tmp_72 = tmp_67*tmp_9;
      real_t tmp_73 = tmp_69*tmp_8;
      real_t tmp_74 = tmp_72 + tmp_73;
      real_t tmp_75 = tmp_1*(tmp_74 - 1.0/3.0) + tmp_5*(tmp_71 - 1.0/3.0);
      real_t tmp_76 = -tmp_68 - tmp_70 - tmp_72 - tmp_73 + 1;
      real_t tmp_77 = tmp_28*tmp_75;
      real_t tmp_78 = 0.11846344252809471*tmp_27;
      real_t tmp_79 = p_affine_10_0*tmp_7 + p_affine_10_1*tmp_10;
      real_t tmp_80 = p_affine_10_0*tmp_8 + p_affine_10_1*tmp_9;
      real_t a_0_0 = tmp_30*(-tmp_11*tmp_24 - tmp_25*tmp_26 + tmp_26*tmp_29) + tmp_42*(-tmp_11*tmp_39 - tmp_25*tmp_40 + tmp_40*tmp_41) + tmp_54*(-tmp_11*tmp_51 - tmp_25*tmp_52 + tmp_52*tmp_53) + tmp_66*(-tmp_11*tmp_63 - tmp_25*tmp_64 + tmp_64*tmp_65) + tmp_78*(-tmp_11*tmp_75 - tmp_25*tmp_76 + tmp_76*tmp_77);
      real_t a_1_0 = tmp_30*(-tmp_20*tmp_25 + tmp_20*tmp_29 - tmp_24*tmp_79) + tmp_42*(-tmp_25*tmp_35 + tmp_35*tmp_41 - tmp_39*tmp_79) + tmp_54*(-tmp_25*tmp_47 + tmp_47*tmp_53 - tmp_51*tmp_79) + tmp_66*(-tmp_25*tmp_59 + tmp_59*tmp_65 - tmp_63*tmp_79) + tmp_78*(-tmp_25*tmp_71 + tmp_71*tmp_77 - tmp_75*tmp_79);
      real_t a_2_0 = tmp_30*(-tmp_23*tmp_25 + tmp_23*tmp_29 - tmp_24*tmp_80) + tmp_42*(-tmp_25*tmp_38 + tmp_38*tmp_41 - tmp_39*tmp_80) + tmp_54*(-tmp_25*tmp_50 + tmp_50*tmp_53 - tmp_51*tmp_80) + tmp_66*(-tmp_25*tmp_62 + tmp_62*tmp_65 - tmp_63*tmp_80) + tmp_78*(-tmp_25*tmp_74 + tmp_74*tmp_77 - tmp_75*tmp_80);
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




class EGVectorLaplaceFormEDGP1_new_0 : public hyteg::dg::DGForm2D
{
 public:
    EGVectorLaplaceFormEDGP1_new_0()

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
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = 1.0 / (tmp_4 - tmp_5*(p_affine_1_1 + tmp_2));
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = tmp_6*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_9 = tmp_1*tmp_8 + tmp_5*tmp_7;
      real_t tmp_10 = tmp_3*tmp_6;
      real_t tmp_11 = tmp_6*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_12 = tmp_11*tmp_5 + tmp_4*tmp_6;
      real_t tmp_13 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_14 = tmp_13*(tmp_12*(-tmp_10 - tmp_11) + tmp_9*(-tmp_7 - tmp_8));
      real_t tmp_15 = tmp_13*(tmp_10*tmp_12 + tmp_8*tmp_9);
      real_t tmp_16 = tmp_13*(tmp_11*tmp_12 + tmp_7*tmp_9);
      real_t a_0_0 = 0.5*tmp_14;
      real_t a_0_1 = 0.5*tmp_15;
      real_t a_0_2 = 0.5*tmp_16;
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
      real_t tmp_2 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_6_1 + tmp_3;
      real_t tmp_5 = 0.046910077030668018*tmp_2 + tmp_4;
      real_t tmp_6 = p_affine_2_1 + tmp_3;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = p_affine_2_0 + tmp_0;
      real_t tmp_9 = 1.0 / (tmp_7 - tmp_8*(p_affine_1_1 + tmp_3));
      real_t tmp_10 = tmp_9*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_11 = tmp_10*tmp_5;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_0;
      real_t tmp_14 = 0.046910077030668018*tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6*tmp_9;
      real_t tmp_16 = tmp_14*tmp_15;
      real_t tmp_17 = tmp_11 + tmp_16;
      real_t tmp_18 = tmp_1*tmp_9;
      real_t tmp_19 = tmp_18*tmp_5;
      real_t tmp_20 = tmp_9*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_21 = tmp_14*tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_1*(tmp_17 - 1.0/3.0) + tmp_8*(tmp_22 - 1.0/3.0);
      real_t tmp_24 = 0.5*p_affine_10_0*(-tmp_15 - tmp_20) + 0.5*p_affine_10_1*(-tmp_10 - tmp_18);
      real_t tmp_25 = -tmp_11 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_26 = 0.5*p_affine_10_0*(tmp_20*tmp_8 + tmp_7*tmp_9) + 0.5*p_affine_10_1*(tmp_1*tmp_10 + tmp_18*tmp_8);
      real_t tmp_27 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_2*tmp_2), 1.0/2.0));
      real_t tmp_28 = 6/tmp_27;
      real_t tmp_29 = tmp_23*tmp_28;
      real_t tmp_30 = 0.11846344252809471*tmp_27;
      real_t tmp_31 = 0.23076534494715845*tmp_2 + tmp_4;
      real_t tmp_32 = tmp_10*tmp_31;
      real_t tmp_33 = 0.23076534494715845*tmp_12 + tmp_13;
      real_t tmp_34 = tmp_15*tmp_33;
      real_t tmp_35 = tmp_32 + tmp_34;
      real_t tmp_36 = tmp_18*tmp_31;
      real_t tmp_37 = tmp_20*tmp_33;
      real_t tmp_38 = tmp_36 + tmp_37;
      real_t tmp_39 = tmp_1*(tmp_35 - 1.0/3.0) + tmp_8*(tmp_38 - 1.0/3.0);
      real_t tmp_40 = -tmp_32 - tmp_34 - tmp_36 - tmp_37 + 1;
      real_t tmp_41 = tmp_28*tmp_39;
      real_t tmp_42 = 0.2393143352496831*tmp_27;
      real_t tmp_43 = 0.5*tmp_2 + tmp_4;
      real_t tmp_44 = tmp_10*tmp_43;
      real_t tmp_45 = 0.5*tmp_12 + tmp_13;
      real_t tmp_46 = tmp_15*tmp_45;
      real_t tmp_47 = tmp_44 + tmp_46;
      real_t tmp_48 = tmp_18*tmp_43;
      real_t tmp_49 = tmp_20*tmp_45;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_1*(tmp_47 - 1.0/3.0) + tmp_8*(tmp_50 - 1.0/3.0);
      real_t tmp_52 = -tmp_44 - tmp_46 - tmp_48 - tmp_49 + 1;
      real_t tmp_53 = tmp_28*tmp_51;
      real_t tmp_54 = 0.2844444444444445*tmp_27;
      real_t tmp_55 = 0.7692346550528415*tmp_2 + tmp_4;
      real_t tmp_56 = tmp_10*tmp_55;
      real_t tmp_57 = 0.7692346550528415*tmp_12 + tmp_13;
      real_t tmp_58 = tmp_15*tmp_57;
      real_t tmp_59 = tmp_56 + tmp_58;
      real_t tmp_60 = tmp_18*tmp_55;
      real_t tmp_61 = tmp_20*tmp_57;
      real_t tmp_62 = tmp_60 + tmp_61;
      real_t tmp_63 = tmp_1*(tmp_59 - 1.0/3.0) + tmp_8*(tmp_62 - 1.0/3.0);
      real_t tmp_64 = -tmp_56 - tmp_58 - tmp_60 - tmp_61 + 1;
      real_t tmp_65 = tmp_28*tmp_63;
      real_t tmp_66 = 0.2393143352496831*tmp_27;
      real_t tmp_67 = 0.95308992296933193*tmp_2 + tmp_4;
      real_t tmp_68 = tmp_10*tmp_67;
      real_t tmp_69 = 0.95308992296933193*tmp_12 + tmp_13;
      real_t tmp_70 = tmp_15*tmp_69;
      real_t tmp_71 = tmp_68 + tmp_70;
      real_t tmp_72 = tmp_18*tmp_67;
      real_t tmp_73 = tmp_20*tmp_69;
      real_t tmp_74 = tmp_72 + tmp_73;
      real_t tmp_75 = tmp_1*(tmp_71 - 1.0/3.0) + tmp_8*(tmp_74 - 1.0/3.0);
      real_t tmp_76 = -tmp_68 - tmp_70 - tmp_72 - tmp_73 + 1;
      real_t tmp_77 = tmp_28*tmp_75;
      real_t tmp_78 = 0.11846344252809471*tmp_27;
      real_t tmp_79 = 0.5*p_affine_10_0*tmp_15 + 0.5*p_affine_10_1*tmp_10;
      real_t tmp_80 = 0.5*p_affine_10_0*tmp_20 + 0.5*p_affine_10_1*tmp_18;
      real_t a_0_0 = tmp_30*(-tmp_23*tmp_24 - tmp_25*tmp_26 + tmp_25*tmp_29) + tmp_42*(-tmp_24*tmp_39 - tmp_26*tmp_40 + tmp_40*tmp_41) + tmp_54*(-tmp_24*tmp_51 - tmp_26*tmp_52 + tmp_52*tmp_53) + tmp_66*(-tmp_24*tmp_63 - tmp_26*tmp_64 + tmp_64*tmp_65) + tmp_78*(-tmp_24*tmp_75 - tmp_26*tmp_76 + tmp_76*tmp_77);
      real_t a_0_1 = tmp_30*(-tmp_17*tmp_26 + tmp_17*tmp_29 - tmp_23*tmp_79) + tmp_42*(-tmp_26*tmp_35 + tmp_35*tmp_41 - tmp_39*tmp_79) + tmp_54*(-tmp_26*tmp_47 + tmp_47*tmp_53 - tmp_51*tmp_79) + tmp_66*(-tmp_26*tmp_59 + tmp_59*tmp_65 - tmp_63*tmp_79) + tmp_78*(-tmp_26*tmp_71 + tmp_71*tmp_77 - tmp_75*tmp_79);
      real_t a_0_2 = tmp_30*(-tmp_22*tmp_26 + tmp_22*tmp_29 - tmp_23*tmp_80) + tmp_42*(-tmp_26*tmp_38 + tmp_38*tmp_41 - tmp_39*tmp_80) + tmp_54*(-tmp_26*tmp_50 + tmp_50*tmp_53 - tmp_51*tmp_80) + tmp_66*(-tmp_26*tmp_62 + tmp_62*tmp_65 - tmp_63*tmp_80) + tmp_78*(-tmp_26*tmp_74 + tmp_74*tmp_77 - tmp_75*tmp_80);
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_4 = p_affine_6_1 + 0.046910077030668018*tmp_3;
      real_t tmp_5 = tmp_2 + tmp_4;
      real_t tmp_6 = p_affine_2_1 + tmp_2;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = p_affine_2_0 + tmp_0;
      real_t tmp_9 = 1.0 / (tmp_7 - tmp_8*(p_affine_1_1 + tmp_2));
      real_t tmp_10 = tmp_9*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.046910077030668018*tmp_11;
      real_t tmp_13 = tmp_9*(tmp_0 + tmp_12);
      real_t tmp_14 = tmp_1*tmp_9;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_1*(tmp_10*tmp_5 + tmp_13*tmp_6 - 1.0/3.0) + tmp_8*(tmp_13*tmp_15 + tmp_14*tmp_5 - 1.0/3.0);
      real_t tmp_17 = -p_affine_3_1;
      real_t tmp_18 = p_affine_5_1 + tmp_17;
      real_t tmp_19 = -p_affine_3_0;
      real_t tmp_20 = p_affine_4_0 + tmp_19;
      real_t tmp_21 = 1.0 / (tmp_18*tmp_20 - (p_affine_4_1 + tmp_17)*(p_affine_5_0 + tmp_19));
      real_t tmp_22 = tmp_18*tmp_21;
      real_t tmp_23 = tmp_21*(p_affine_3_1 - p_affine_4_1);
      real_t tmp_24 = tmp_20*tmp_21;
      real_t tmp_25 = tmp_21*(p_affine_3_0 - p_affine_5_0);
      real_t tmp_26 = 0.5*p_affine_10_0*(-tmp_22 - tmp_23) + 0.5*p_affine_10_1*(-tmp_24 - tmp_25);
      real_t tmp_27 = tmp_17 + tmp_4;
      real_t tmp_28 = tmp_24*tmp_27;
      real_t tmp_29 = tmp_25*tmp_27;
      real_t tmp_30 = tmp_12 + tmp_19;
      real_t tmp_31 = tmp_22*tmp_30;
      real_t tmp_32 = tmp_23*tmp_30;
      real_t tmp_33 = -tmp_28 - tmp_29 - tmp_31 - tmp_32 + 1;
      real_t tmp_34 = tmp_8*tmp_9;
      real_t tmp_35 = 0.5*p_affine_10_0*(tmp_15*tmp_34 + tmp_7*tmp_9) + 0.5*p_affine_10_1*(tmp_1*tmp_10 + tmp_1*tmp_34);
      real_t tmp_36 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_3*tmp_3), 1.0/2.0));
      real_t tmp_37 = 6/tmp_36;
      real_t tmp_38 = tmp_16*tmp_37;
      real_t tmp_39 = 0.11846344252809471*tmp_36;
      real_t tmp_40 = p_affine_6_1 + 0.23076534494715845*tmp_3;
      real_t tmp_41 = tmp_2 + tmp_40;
      real_t tmp_42 = p_affine_6_0 + 0.23076534494715845*tmp_11;
      real_t tmp_43 = tmp_9*(tmp_0 + tmp_42);
      real_t tmp_44 = tmp_1*(tmp_10*tmp_41 + tmp_43*tmp_6 - 1.0/3.0) + tmp_8*(tmp_14*tmp_41 + tmp_15*tmp_43 - 1.0/3.0);
      real_t tmp_45 = tmp_17 + tmp_40;
      real_t tmp_46 = tmp_24*tmp_45;
      real_t tmp_47 = tmp_25*tmp_45;
      real_t tmp_48 = tmp_19 + tmp_42;
      real_t tmp_49 = tmp_22*tmp_48;
      real_t tmp_50 = tmp_23*tmp_48;
      real_t tmp_51 = -tmp_46 - tmp_47 - tmp_49 - tmp_50 + 1;
      real_t tmp_52 = tmp_37*tmp_44;
      real_t tmp_53 = 0.2393143352496831*tmp_36;
      real_t tmp_54 = p_affine_6_1 + 0.5*tmp_3;
      real_t tmp_55 = tmp_2 + tmp_54;
      real_t tmp_56 = p_affine_6_0 + 0.5*tmp_11;
      real_t tmp_57 = tmp_9*(tmp_0 + tmp_56);
      real_t tmp_58 = tmp_1*(tmp_10*tmp_55 + tmp_57*tmp_6 - 1.0/3.0) + tmp_8*(tmp_14*tmp_55 + tmp_15*tmp_57 - 1.0/3.0);
      real_t tmp_59 = tmp_17 + tmp_54;
      real_t tmp_60 = tmp_24*tmp_59;
      real_t tmp_61 = tmp_25*tmp_59;
      real_t tmp_62 = tmp_19 + tmp_56;
      real_t tmp_63 = tmp_22*tmp_62;
      real_t tmp_64 = tmp_23*tmp_62;
      real_t tmp_65 = -tmp_60 - tmp_61 - tmp_63 - tmp_64 + 1;
      real_t tmp_66 = tmp_37*tmp_58;
      real_t tmp_67 = 0.2844444444444445*tmp_36;
      real_t tmp_68 = p_affine_6_1 + 0.7692346550528415*tmp_3;
      real_t tmp_69 = tmp_2 + tmp_68;
      real_t tmp_70 = p_affine_6_0 + 0.7692346550528415*tmp_11;
      real_t tmp_71 = tmp_9*(tmp_0 + tmp_70);
      real_t tmp_72 = tmp_1*(tmp_10*tmp_69 + tmp_6*tmp_71 - 1.0/3.0) + tmp_8*(tmp_14*tmp_69 + tmp_15*tmp_71 - 1.0/3.0);
      real_t tmp_73 = tmp_17 + tmp_68;
      real_t tmp_74 = tmp_24*tmp_73;
      real_t tmp_75 = tmp_25*tmp_73;
      real_t tmp_76 = tmp_19 + tmp_70;
      real_t tmp_77 = tmp_22*tmp_76;
      real_t tmp_78 = tmp_23*tmp_76;
      real_t tmp_79 = -tmp_74 - tmp_75 - tmp_77 - tmp_78 + 1;
      real_t tmp_80 = tmp_37*tmp_72;
      real_t tmp_81 = 0.2393143352496831*tmp_36;
      real_t tmp_82 = p_affine_6_1 + 0.95308992296933193*tmp_3;
      real_t tmp_83 = tmp_2 + tmp_82;
      real_t tmp_84 = p_affine_6_0 + 0.95308992296933193*tmp_11;
      real_t tmp_85 = tmp_9*(tmp_0 + tmp_84);
      real_t tmp_86 = tmp_1*(tmp_10*tmp_83 + tmp_6*tmp_85 - 1.0/3.0) + tmp_8*(tmp_14*tmp_83 + tmp_15*tmp_85 - 1.0/3.0);
      real_t tmp_87 = tmp_17 + tmp_82;
      real_t tmp_88 = tmp_24*tmp_87;
      real_t tmp_89 = tmp_25*tmp_87;
      real_t tmp_90 = tmp_19 + tmp_84;
      real_t tmp_91 = tmp_22*tmp_90;
      real_t tmp_92 = tmp_23*tmp_90;
      real_t tmp_93 = -tmp_88 - tmp_89 - tmp_91 - tmp_92 + 1;
      real_t tmp_94 = tmp_37*tmp_86;
      real_t tmp_95 = 0.11846344252809471*tmp_36;
      real_t tmp_96 = tmp_29 + tmp_31;
      real_t tmp_97 = 0.5*p_affine_10_0*tmp_22 + 0.5*p_affine_10_1*tmp_25;
      real_t tmp_98 = tmp_47 + tmp_49;
      real_t tmp_99 = tmp_61 + tmp_63;
      real_t tmp_100 = tmp_75 + tmp_77;
      real_t tmp_101 = tmp_89 + tmp_91;
      real_t tmp_102 = tmp_28 + tmp_32;
      real_t tmp_103 = 0.5*p_affine_10_0*tmp_23 + 0.5*p_affine_10_1*tmp_24;
      real_t tmp_104 = tmp_46 + tmp_50;
      real_t tmp_105 = tmp_60 + tmp_64;
      real_t tmp_106 = tmp_74 + tmp_78;
      real_t tmp_107 = tmp_88 + tmp_92;
      real_t a_0_0 = tmp_39*(-tmp_16*tmp_26 + tmp_33*tmp_35 - tmp_33*tmp_38) + tmp_53*(-tmp_26*tmp_44 + tmp_35*tmp_51 - tmp_51*tmp_52) + tmp_67*(-tmp_26*tmp_58 + tmp_35*tmp_65 - tmp_65*tmp_66) + tmp_81*(-tmp_26*tmp_72 + tmp_35*tmp_79 - tmp_79*tmp_80) + tmp_95*(-tmp_26*tmp_86 + tmp_35*tmp_93 - tmp_93*tmp_94);
      real_t a_0_1 = tmp_39*(-tmp_16*tmp_97 + tmp_35*tmp_96 - tmp_38*tmp_96) + tmp_53*(tmp_35*tmp_98 - tmp_44*tmp_97 - tmp_52*tmp_98) + tmp_67*(tmp_35*tmp_99 - tmp_58*tmp_97 - tmp_66*tmp_99) + tmp_81*(tmp_100*tmp_35 - tmp_100*tmp_80 - tmp_72*tmp_97) + tmp_95*(tmp_101*tmp_35 - tmp_101*tmp_94 - tmp_86*tmp_97);
      real_t a_0_2 = tmp_39*(tmp_102*tmp_35 - tmp_102*tmp_38 - tmp_103*tmp_16) + tmp_53*(-tmp_103*tmp_44 + tmp_104*tmp_35 - tmp_104*tmp_52) + tmp_67*(-tmp_103*tmp_58 + tmp_105*tmp_35 - tmp_105*tmp_66) + tmp_81*(-tmp_103*tmp_72 + tmp_106*tmp_35 - tmp_106*tmp_80) + tmp_95*(-tmp_103*tmp_86 + tmp_107*tmp_35 - tmp_107*tmp_94);
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

      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_0;
      real_t tmp_3 = p_affine_1_0 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_2;
      real_t tmp_6 = 1.0 / (tmp_4 - tmp_5*(p_affine_1_1 + tmp_0));
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = tmp_6*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_9 = tmp_3*tmp_6;
      real_t tmp_10 = tmp_6*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_11 = p_affine_10_0*(-tmp_7 - tmp_8) + p_affine_10_1*(-tmp_10 - tmp_9);
      real_t tmp_12 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_13 = p_affine_6_1 + tmp_0;
      real_t tmp_14 = 0.046910077030668018*tmp_12 + tmp_13;
      real_t tmp_15 = tmp_10*tmp_14;
      real_t tmp_16 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_17 = p_affine_6_0 + tmp_2;
      real_t tmp_18 = 0.046910077030668018*tmp_16 + tmp_17;
      real_t tmp_19 = tmp_18*tmp_7;
      real_t tmp_20 = tmp_15 + tmp_19;
      real_t tmp_21 = tmp_14*tmp_9;
      real_t tmp_22 = tmp_18*tmp_8;
      real_t tmp_23 = tmp_21 + tmp_22;
      real_t tmp_24 = tmp_3*(tmp_20 - 1.0/3.0) + tmp_5*(tmp_23 - 1.0/3.0);
      real_t tmp_25 = p_affine_10_0*(tmp_4*tmp_6 + tmp_5*tmp_8) + p_affine_10_1*(tmp_10*tmp_3 + tmp_5*tmp_9);
      real_t tmp_26 = -tmp_15 - tmp_19 - tmp_21 - tmp_22 + 1;
      real_t tmp_27 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_16*tmp_16), 1.0/2.0));
      real_t tmp_28 = 6/tmp_27;
      real_t tmp_29 = tmp_24*tmp_28;
      real_t tmp_30 = 0.11846344252809471*tmp_27;
      real_t tmp_31 = 0.23076534494715845*tmp_12 + tmp_13;
      real_t tmp_32 = tmp_10*tmp_31;
      real_t tmp_33 = 0.23076534494715845*tmp_16 + tmp_17;
      real_t tmp_34 = tmp_33*tmp_7;
      real_t tmp_35 = tmp_32 + tmp_34;
      real_t tmp_36 = tmp_31*tmp_9;
      real_t tmp_37 = tmp_33*tmp_8;
      real_t tmp_38 = tmp_36 + tmp_37;
      real_t tmp_39 = tmp_3*(tmp_35 - 1.0/3.0) + tmp_5*(tmp_38 - 1.0/3.0);
      real_t tmp_40 = -tmp_32 - tmp_34 - tmp_36 - tmp_37 + 1;
      real_t tmp_41 = tmp_28*tmp_39;
      real_t tmp_42 = 0.2393143352496831*tmp_27;
      real_t tmp_43 = 0.5*tmp_12 + tmp_13;
      real_t tmp_44 = tmp_10*tmp_43;
      real_t tmp_45 = 0.5*tmp_16 + tmp_17;
      real_t tmp_46 = tmp_45*tmp_7;
      real_t tmp_47 = tmp_44 + tmp_46;
      real_t tmp_48 = tmp_43*tmp_9;
      real_t tmp_49 = tmp_45*tmp_8;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_3*(tmp_47 - 1.0/3.0) + tmp_5*(tmp_50 - 1.0/3.0);
      real_t tmp_52 = -tmp_44 - tmp_46 - tmp_48 - tmp_49 + 1;
      real_t tmp_53 = tmp_28*tmp_51;
      real_t tmp_54 = 0.2844444444444445*tmp_27;
      real_t tmp_55 = 0.7692346550528415*tmp_12 + tmp_13;
      real_t tmp_56 = tmp_10*tmp_55;
      real_t tmp_57 = 0.7692346550528415*tmp_16 + tmp_17;
      real_t tmp_58 = tmp_57*tmp_7;
      real_t tmp_59 = tmp_56 + tmp_58;
      real_t tmp_60 = tmp_55*tmp_9;
      real_t tmp_61 = tmp_57*tmp_8;
      real_t tmp_62 = tmp_60 + tmp_61;
      real_t tmp_63 = tmp_3*(tmp_59 - 1.0/3.0) + tmp_5*(tmp_62 - 1.0/3.0);
      real_t tmp_64 = -tmp_56 - tmp_58 - tmp_60 - tmp_61 + 1;
      real_t tmp_65 = tmp_28*tmp_63;
      real_t tmp_66 = 0.2393143352496831*tmp_27;
      real_t tmp_67 = 0.95308992296933193*tmp_12 + tmp_13;
      real_t tmp_68 = tmp_10*tmp_67;
      real_t tmp_69 = 0.95308992296933193*tmp_16 + tmp_17;
      real_t tmp_70 = tmp_69*tmp_7;
      real_t tmp_71 = tmp_68 + tmp_70;
      real_t tmp_72 = tmp_67*tmp_9;
      real_t tmp_73 = tmp_69*tmp_8;
      real_t tmp_74 = tmp_72 + tmp_73;
      real_t tmp_75 = tmp_3*(tmp_71 - 1.0/3.0) + tmp_5*(tmp_74 - 1.0/3.0);
      real_t tmp_76 = -tmp_68 - tmp_70 - tmp_72 - tmp_73 + 1;
      real_t tmp_77 = tmp_28*tmp_75;
      real_t tmp_78 = 0.11846344252809471*tmp_27;
      real_t tmp_79 = p_affine_10_0*tmp_7 + p_affine_10_1*tmp_10;
      real_t tmp_80 = p_affine_10_0*tmp_8 + p_affine_10_1*tmp_9;
      real_t a_0_0 = tmp_30*(-tmp_11*tmp_24 - tmp_25*tmp_26 + tmp_26*tmp_29) + tmp_42*(-tmp_11*tmp_39 - tmp_25*tmp_40 + tmp_40*tmp_41) + tmp_54*(-tmp_11*tmp_51 - tmp_25*tmp_52 + tmp_52*tmp_53) + tmp_66*(-tmp_11*tmp_63 - tmp_25*tmp_64 + tmp_64*tmp_65) + tmp_78*(-tmp_11*tmp_75 - tmp_25*tmp_76 + tmp_76*tmp_77);
      real_t a_0_1 = tmp_30*(-tmp_20*tmp_25 + tmp_20*tmp_29 - tmp_24*tmp_79) + tmp_42*(-tmp_25*tmp_35 + tmp_35*tmp_41 - tmp_39*tmp_79) + tmp_54*(-tmp_25*tmp_47 + tmp_47*tmp_53 - tmp_51*tmp_79) + tmp_66*(-tmp_25*tmp_59 + tmp_59*tmp_65 - tmp_63*tmp_79) + tmp_78*(-tmp_25*tmp_71 + tmp_71*tmp_77 - tmp_75*tmp_79);
      real_t a_0_2 = tmp_30*(-tmp_23*tmp_25 + tmp_23*tmp_29 - tmp_24*tmp_80) + tmp_42*(-tmp_25*tmp_38 + tmp_38*tmp_41 - tmp_39*tmp_80) + tmp_54*(-tmp_25*tmp_50 + tmp_50*tmp_53 - tmp_51*tmp_80) + tmp_66*(-tmp_25*tmp_62 + tmp_62*tmp_65 - tmp_63*tmp_80) + tmp_78*(-tmp_25*tmp_74 + tmp_74*tmp_77 - tmp_75*tmp_80);
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
      elMat( 0, 0) = 0;
   }

public:


};




class EGVectorLaplaceFormEDGP1_new_1 : public hyteg::dg::DGForm2D
{
 public:
    EGVectorLaplaceFormEDGP1_new_1()

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
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = 1.0 / (tmp_4 - tmp_5*(p_affine_2_0 + tmp_0));
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = tmp_6*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_9 = tmp_4*tmp_6 + tmp_5*tmp_8;
      real_t tmp_10 = tmp_3*tmp_6;
      real_t tmp_11 = tmp_6*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_12 = tmp_10*tmp_5 + tmp_11*tmp_3;
      real_t tmp_13 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_14 = tmp_13*(tmp_12*(-tmp_10 - tmp_11) + tmp_9*(-tmp_7 - tmp_8));
      real_t tmp_15 = tmp_13*(tmp_10*tmp_12 + tmp_8*tmp_9);
      real_t tmp_16 = tmp_13*(tmp_11*tmp_12 + tmp_7*tmp_9);
      real_t a_0_0 = 0.5*tmp_14;
      real_t a_0_1 = 0.5*tmp_15;
      real_t a_0_2 = 0.5*tmp_16;
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

      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_1_1 + tmp_0;
      real_t tmp_2 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_3 = p_affine_6_1 + tmp_0;
      real_t tmp_4 = 0.046910077030668018*tmp_2 + tmp_3;
      real_t tmp_5 = -p_affine_0_0;
      real_t tmp_6 = p_affine_1_0 + tmp_5;
      real_t tmp_7 = p_affine_2_1 + tmp_0;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = 1.0 / (-tmp_1*(p_affine_2_0 + tmp_5) + tmp_8);
      real_t tmp_10 = tmp_9*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_11 = tmp_10*tmp_4;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_5;
      real_t tmp_14 = 0.046910077030668018*tmp_12 + tmp_13;
      real_t tmp_15 = tmp_7*tmp_9;
      real_t tmp_16 = tmp_14*tmp_15;
      real_t tmp_17 = tmp_11 + tmp_16;
      real_t tmp_18 = tmp_6*tmp_9;
      real_t tmp_19 = tmp_18*tmp_4;
      real_t tmp_20 = tmp_9*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_21 = tmp_14*tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_1*(tmp_17 - 1.0/3.0) + tmp_7*(tmp_22 - 1.0/3.0);
      real_t tmp_24 = 0.5*p_affine_10_0*(-tmp_15 - tmp_20) + 0.5*p_affine_10_1*(-tmp_10 - tmp_18);
      real_t tmp_25 = -tmp_11 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_26 = 0.5*p_affine_10_0*(tmp_1*tmp_15 + tmp_20*tmp_7) + 0.5*p_affine_10_1*(tmp_1*tmp_10 + tmp_8*tmp_9);
      real_t tmp_27 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_2*tmp_2), 1.0/2.0));
      real_t tmp_28 = 6/tmp_27;
      real_t tmp_29 = tmp_23*tmp_28;
      real_t tmp_30 = 0.11846344252809471*tmp_27;
      real_t tmp_31 = 0.23076534494715845*tmp_2 + tmp_3;
      real_t tmp_32 = tmp_10*tmp_31;
      real_t tmp_33 = 0.23076534494715845*tmp_12 + tmp_13;
      real_t tmp_34 = tmp_15*tmp_33;
      real_t tmp_35 = tmp_32 + tmp_34;
      real_t tmp_36 = tmp_18*tmp_31;
      real_t tmp_37 = tmp_20*tmp_33;
      real_t tmp_38 = tmp_36 + tmp_37;
      real_t tmp_39 = tmp_1*(tmp_35 - 1.0/3.0) + tmp_7*(tmp_38 - 1.0/3.0);
      real_t tmp_40 = -tmp_32 - tmp_34 - tmp_36 - tmp_37 + 1;
      real_t tmp_41 = tmp_28*tmp_39;
      real_t tmp_42 = 0.2393143352496831*tmp_27;
      real_t tmp_43 = 0.5*tmp_2 + tmp_3;
      real_t tmp_44 = tmp_10*tmp_43;
      real_t tmp_45 = 0.5*tmp_12 + tmp_13;
      real_t tmp_46 = tmp_15*tmp_45;
      real_t tmp_47 = tmp_44 + tmp_46;
      real_t tmp_48 = tmp_18*tmp_43;
      real_t tmp_49 = tmp_20*tmp_45;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_1*(tmp_47 - 1.0/3.0) + tmp_7*(tmp_50 - 1.0/3.0);
      real_t tmp_52 = -tmp_44 - tmp_46 - tmp_48 - tmp_49 + 1;
      real_t tmp_53 = tmp_28*tmp_51;
      real_t tmp_54 = 0.2844444444444445*tmp_27;
      real_t tmp_55 = 0.7692346550528415*tmp_2 + tmp_3;
      real_t tmp_56 = tmp_10*tmp_55;
      real_t tmp_57 = 0.7692346550528415*tmp_12 + tmp_13;
      real_t tmp_58 = tmp_15*tmp_57;
      real_t tmp_59 = tmp_56 + tmp_58;
      real_t tmp_60 = tmp_18*tmp_55;
      real_t tmp_61 = tmp_20*tmp_57;
      real_t tmp_62 = tmp_60 + tmp_61;
      real_t tmp_63 = tmp_1*(tmp_59 - 1.0/3.0) + tmp_7*(tmp_62 - 1.0/3.0);
      real_t tmp_64 = -tmp_56 - tmp_58 - tmp_60 - tmp_61 + 1;
      real_t tmp_65 = tmp_28*tmp_63;
      real_t tmp_66 = 0.2393143352496831*tmp_27;
      real_t tmp_67 = 0.95308992296933193*tmp_2 + tmp_3;
      real_t tmp_68 = tmp_10*tmp_67;
      real_t tmp_69 = 0.95308992296933193*tmp_12 + tmp_13;
      real_t tmp_70 = tmp_15*tmp_69;
      real_t tmp_71 = tmp_68 + tmp_70;
      real_t tmp_72 = tmp_18*tmp_67;
      real_t tmp_73 = tmp_20*tmp_69;
      real_t tmp_74 = tmp_72 + tmp_73;
      real_t tmp_75 = tmp_1*(tmp_71 - 1.0/3.0) + tmp_7*(tmp_74 - 1.0/3.0);
      real_t tmp_76 = -tmp_68 - tmp_70 - tmp_72 - tmp_73 + 1;
      real_t tmp_77 = tmp_28*tmp_75;
      real_t tmp_78 = 0.11846344252809471*tmp_27;
      real_t tmp_79 = 0.5*p_affine_10_0*tmp_15 + 0.5*p_affine_10_1*tmp_10;
      real_t tmp_80 = 0.5*p_affine_10_0*tmp_20 + 0.5*p_affine_10_1*tmp_18;
      real_t a_0_0 = tmp_30*(-tmp_23*tmp_24 - tmp_25*tmp_26 + tmp_25*tmp_29) + tmp_42*(-tmp_24*tmp_39 - tmp_26*tmp_40 + tmp_40*tmp_41) + tmp_54*(-tmp_24*tmp_51 - tmp_26*tmp_52 + tmp_52*tmp_53) + tmp_66*(-tmp_24*tmp_63 - tmp_26*tmp_64 + tmp_64*tmp_65) + tmp_78*(-tmp_24*tmp_75 - tmp_26*tmp_76 + tmp_76*tmp_77);
      real_t a_0_1 = tmp_30*(-tmp_17*tmp_26 + tmp_17*tmp_29 - tmp_23*tmp_79) + tmp_42*(-tmp_26*tmp_35 + tmp_35*tmp_41 - tmp_39*tmp_79) + tmp_54*(-tmp_26*tmp_47 + tmp_47*tmp_53 - tmp_51*tmp_79) + tmp_66*(-tmp_26*tmp_59 + tmp_59*tmp_65 - tmp_63*tmp_79) + tmp_78*(-tmp_26*tmp_71 + tmp_71*tmp_77 - tmp_75*tmp_79);
      real_t a_0_2 = tmp_30*(-tmp_22*tmp_26 + tmp_22*tmp_29 - tmp_23*tmp_80) + tmp_42*(-tmp_26*tmp_38 + tmp_38*tmp_41 - tmp_39*tmp_80) + tmp_54*(-tmp_26*tmp_50 + tmp_50*tmp_53 - tmp_51*tmp_80) + tmp_66*(-tmp_26*tmp_62 + tmp_62*tmp_65 - tmp_63*tmp_80) + tmp_78*(-tmp_26*tmp_74 + tmp_74*tmp_77 - tmp_75*tmp_80);
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

      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_1_1 + tmp_0;
      real_t tmp_2 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_3 = p_affine_6_1 + 0.046910077030668018*tmp_2;
      real_t tmp_4 = tmp_0 + tmp_3;
      real_t tmp_5 = -p_affine_0_0;
      real_t tmp_6 = p_affine_1_0 + tmp_5;
      real_t tmp_7 = p_affine_2_1 + tmp_0;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = 1.0 / (-tmp_1*(p_affine_2_0 + tmp_5) + tmp_8);
      real_t tmp_10 = tmp_9*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.046910077030668018*tmp_11;
      real_t tmp_13 = tmp_12 + tmp_5;
      real_t tmp_14 = tmp_7*tmp_9;
      real_t tmp_15 = tmp_6*tmp_9;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_16*tmp_9;
      real_t tmp_18 = tmp_1*(tmp_10*tmp_4 + tmp_13*tmp_14 - 1.0/3.0) + tmp_7*(tmp_13*tmp_17 + tmp_15*tmp_4 - 1.0/3.0);
      real_t tmp_19 = -p_affine_3_1;
      real_t tmp_20 = p_affine_5_1 + tmp_19;
      real_t tmp_21 = -p_affine_3_0;
      real_t tmp_22 = p_affine_4_0 + tmp_21;
      real_t tmp_23 = 1.0 / (tmp_20*tmp_22 - (p_affine_4_1 + tmp_19)*(p_affine_5_0 + tmp_21));
      real_t tmp_24 = tmp_20*tmp_23;
      real_t tmp_25 = tmp_23*(p_affine_3_1 - p_affine_4_1);
      real_t tmp_26 = tmp_22*tmp_23;
      real_t tmp_27 = tmp_23*(p_affine_3_0 - p_affine_5_0);
      real_t tmp_28 = 0.5*p_affine_10_0*(-tmp_24 - tmp_25) + 0.5*p_affine_10_1*(-tmp_26 - tmp_27);
      real_t tmp_29 = tmp_19 + tmp_3;
      real_t tmp_30 = tmp_26*tmp_29;
      real_t tmp_31 = tmp_27*tmp_29;
      real_t tmp_32 = tmp_12 + tmp_21;
      real_t tmp_33 = tmp_24*tmp_32;
      real_t tmp_34 = tmp_25*tmp_32;
      real_t tmp_35 = -tmp_30 - tmp_31 - tmp_33 - tmp_34 + 1;
      real_t tmp_36 = 0.5*p_affine_10_0*(tmp_1*tmp_14 + tmp_14*tmp_16) + 0.5*p_affine_10_1*(tmp_1*tmp_10 + tmp_8*tmp_9);
      real_t tmp_37 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_2*tmp_2), 1.0/2.0));
      real_t tmp_38 = 6/tmp_37;
      real_t tmp_39 = tmp_18*tmp_38;
      real_t tmp_40 = 0.11846344252809471*tmp_37;
      real_t tmp_41 = p_affine_6_1 + 0.23076534494715845*tmp_2;
      real_t tmp_42 = tmp_0 + tmp_41;
      real_t tmp_43 = p_affine_6_0 + 0.23076534494715845*tmp_11;
      real_t tmp_44 = tmp_43 + tmp_5;
      real_t tmp_45 = tmp_1*(tmp_10*tmp_42 + tmp_14*tmp_44 - 1.0/3.0) + tmp_7*(tmp_15*tmp_42 + tmp_17*tmp_44 - 1.0/3.0);
      real_t tmp_46 = tmp_19 + tmp_41;
      real_t tmp_47 = tmp_26*tmp_46;
      real_t tmp_48 = tmp_27*tmp_46;
      real_t tmp_49 = tmp_21 + tmp_43;
      real_t tmp_50 = tmp_24*tmp_49;
      real_t tmp_51 = tmp_25*tmp_49;
      real_t tmp_52 = -tmp_47 - tmp_48 - tmp_50 - tmp_51 + 1;
      real_t tmp_53 = tmp_38*tmp_45;
      real_t tmp_54 = 0.2393143352496831*tmp_37;
      real_t tmp_55 = p_affine_6_1 + 0.5*tmp_2;
      real_t tmp_56 = tmp_0 + tmp_55;
      real_t tmp_57 = p_affine_6_0 + 0.5*tmp_11;
      real_t tmp_58 = tmp_5 + tmp_57;
      real_t tmp_59 = tmp_1*(tmp_10*tmp_56 + tmp_14*tmp_58 - 1.0/3.0) + tmp_7*(tmp_15*tmp_56 + tmp_17*tmp_58 - 1.0/3.0);
      real_t tmp_60 = tmp_19 + tmp_55;
      real_t tmp_61 = tmp_26*tmp_60;
      real_t tmp_62 = tmp_27*tmp_60;
      real_t tmp_63 = tmp_21 + tmp_57;
      real_t tmp_64 = tmp_24*tmp_63;
      real_t tmp_65 = tmp_25*tmp_63;
      real_t tmp_66 = -tmp_61 - tmp_62 - tmp_64 - tmp_65 + 1;
      real_t tmp_67 = tmp_38*tmp_59;
      real_t tmp_68 = 0.2844444444444445*tmp_37;
      real_t tmp_69 = p_affine_6_1 + 0.7692346550528415*tmp_2;
      real_t tmp_70 = tmp_0 + tmp_69;
      real_t tmp_71 = p_affine_6_0 + 0.7692346550528415*tmp_11;
      real_t tmp_72 = tmp_5 + tmp_71;
      real_t tmp_73 = tmp_1*(tmp_10*tmp_70 + tmp_14*tmp_72 - 1.0/3.0) + tmp_7*(tmp_15*tmp_70 + tmp_17*tmp_72 - 1.0/3.0);
      real_t tmp_74 = tmp_19 + tmp_69;
      real_t tmp_75 = tmp_26*tmp_74;
      real_t tmp_76 = tmp_27*tmp_74;
      real_t tmp_77 = tmp_21 + tmp_71;
      real_t tmp_78 = tmp_24*tmp_77;
      real_t tmp_79 = tmp_25*tmp_77;
      real_t tmp_80 = -tmp_75 - tmp_76 - tmp_78 - tmp_79 + 1;
      real_t tmp_81 = tmp_38*tmp_73;
      real_t tmp_82 = 0.2393143352496831*tmp_37;
      real_t tmp_83 = p_affine_6_1 + 0.95308992296933193*tmp_2;
      real_t tmp_84 = tmp_0 + tmp_83;
      real_t tmp_85 = p_affine_6_0 + 0.95308992296933193*tmp_11;
      real_t tmp_86 = tmp_5 + tmp_85;
      real_t tmp_87 = tmp_1*(tmp_10*tmp_84 + tmp_14*tmp_86 - 1.0/3.0) + tmp_7*(tmp_15*tmp_84 + tmp_17*tmp_86 - 1.0/3.0);
      real_t tmp_88 = tmp_19 + tmp_83;
      real_t tmp_89 = tmp_26*tmp_88;
      real_t tmp_90 = tmp_27*tmp_88;
      real_t tmp_91 = tmp_21 + tmp_85;
      real_t tmp_92 = tmp_24*tmp_91;
      real_t tmp_93 = tmp_25*tmp_91;
      real_t tmp_94 = -tmp_89 - tmp_90 - tmp_92 - tmp_93 + 1;
      real_t tmp_95 = tmp_38*tmp_87;
      real_t tmp_96 = 0.11846344252809471*tmp_37;
      real_t tmp_97 = tmp_31 + tmp_33;
      real_t tmp_98 = 0.5*p_affine_10_0*tmp_24 + 0.5*p_affine_10_1*tmp_27;
      real_t tmp_99 = tmp_48 + tmp_50;
      real_t tmp_100 = tmp_62 + tmp_64;
      real_t tmp_101 = tmp_76 + tmp_78;
      real_t tmp_102 = tmp_90 + tmp_92;
      real_t tmp_103 = tmp_30 + tmp_34;
      real_t tmp_104 = 0.5*p_affine_10_0*tmp_25 + 0.5*p_affine_10_1*tmp_26;
      real_t tmp_105 = tmp_47 + tmp_51;
      real_t tmp_106 = tmp_61 + tmp_65;
      real_t tmp_107 = tmp_75 + tmp_79;
      real_t tmp_108 = tmp_89 + tmp_93;
      real_t a_0_0 = tmp_40*(-tmp_18*tmp_28 + tmp_35*tmp_36 - tmp_35*tmp_39) + tmp_54*(-tmp_28*tmp_45 + tmp_36*tmp_52 - tmp_52*tmp_53) + tmp_68*(-tmp_28*tmp_59 + tmp_36*tmp_66 - tmp_66*tmp_67) + tmp_82*(-tmp_28*tmp_73 + tmp_36*tmp_80 - tmp_80*tmp_81) + tmp_96*(-tmp_28*tmp_87 + tmp_36*tmp_94 - tmp_94*tmp_95);
      real_t a_0_1 = tmp_40*(-tmp_18*tmp_98 + tmp_36*tmp_97 - tmp_39*tmp_97) + tmp_54*(tmp_36*tmp_99 - tmp_45*tmp_98 - tmp_53*tmp_99) + tmp_68*(tmp_100*tmp_36 - tmp_100*tmp_67 - tmp_59*tmp_98) + tmp_82*(tmp_101*tmp_36 - tmp_101*tmp_81 - tmp_73*tmp_98) + tmp_96*(tmp_102*tmp_36 - tmp_102*tmp_95 - tmp_87*tmp_98);
      real_t a_0_2 = tmp_40*(tmp_103*tmp_36 - tmp_103*tmp_39 - tmp_104*tmp_18) + tmp_54*(-tmp_104*tmp_45 + tmp_105*tmp_36 - tmp_105*tmp_53) + tmp_68*(-tmp_104*tmp_59 + tmp_106*tmp_36 - tmp_106*tmp_67) + tmp_82*(-tmp_104*tmp_73 + tmp_107*tmp_36 - tmp_107*tmp_81) + tmp_96*(-tmp_104*tmp_87 + tmp_108*tmp_36 - tmp_108*tmp_95);
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

      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_0;
      real_t tmp_3 = p_affine_1_0 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_1_1 + tmp_0;
      real_t tmp_6 = 1.0 / (tmp_4 - tmp_5*(p_affine_2_0 + tmp_2));
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = tmp_6*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_9 = tmp_3*tmp_6;
      real_t tmp_10 = tmp_6*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_11 = p_affine_10_0*(-tmp_7 - tmp_8) + p_affine_10_1*(-tmp_10 - tmp_9);
      real_t tmp_12 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_13 = p_affine_6_1 + tmp_0;
      real_t tmp_14 = 0.046910077030668018*tmp_12 + tmp_13;
      real_t tmp_15 = tmp_10*tmp_14;
      real_t tmp_16 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_17 = p_affine_6_0 + tmp_2;
      real_t tmp_18 = 0.046910077030668018*tmp_16 + tmp_17;
      real_t tmp_19 = tmp_18*tmp_7;
      real_t tmp_20 = tmp_15 + tmp_19;
      real_t tmp_21 = tmp_14*tmp_9;
      real_t tmp_22 = tmp_18*tmp_8;
      real_t tmp_23 = tmp_21 + tmp_22;
      real_t tmp_24 = tmp_1*(tmp_23 - 1.0/3.0) + tmp_5*(tmp_20 - 1.0/3.0);
      real_t tmp_25 = p_affine_10_0*(tmp_1*tmp_8 + tmp_5*tmp_7) + p_affine_10_1*(tmp_10*tmp_5 + tmp_4*tmp_6);
      real_t tmp_26 = -tmp_15 - tmp_19 - tmp_21 - tmp_22 + 1;
      real_t tmp_27 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_16*tmp_16), 1.0/2.0));
      real_t tmp_28 = 6/tmp_27;
      real_t tmp_29 = tmp_24*tmp_28;
      real_t tmp_30 = 0.11846344252809471*tmp_27;
      real_t tmp_31 = 0.23076534494715845*tmp_12 + tmp_13;
      real_t tmp_32 = tmp_10*tmp_31;
      real_t tmp_33 = 0.23076534494715845*tmp_16 + tmp_17;
      real_t tmp_34 = tmp_33*tmp_7;
      real_t tmp_35 = tmp_32 + tmp_34;
      real_t tmp_36 = tmp_31*tmp_9;
      real_t tmp_37 = tmp_33*tmp_8;
      real_t tmp_38 = tmp_36 + tmp_37;
      real_t tmp_39 = tmp_1*(tmp_38 - 1.0/3.0) + tmp_5*(tmp_35 - 1.0/3.0);
      real_t tmp_40 = -tmp_32 - tmp_34 - tmp_36 - tmp_37 + 1;
      real_t tmp_41 = tmp_28*tmp_39;
      real_t tmp_42 = 0.2393143352496831*tmp_27;
      real_t tmp_43 = 0.5*tmp_12 + tmp_13;
      real_t tmp_44 = tmp_10*tmp_43;
      real_t tmp_45 = 0.5*tmp_16 + tmp_17;
      real_t tmp_46 = tmp_45*tmp_7;
      real_t tmp_47 = tmp_44 + tmp_46;
      real_t tmp_48 = tmp_43*tmp_9;
      real_t tmp_49 = tmp_45*tmp_8;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_1*(tmp_50 - 1.0/3.0) + tmp_5*(tmp_47 - 1.0/3.0);
      real_t tmp_52 = -tmp_44 - tmp_46 - tmp_48 - tmp_49 + 1;
      real_t tmp_53 = tmp_28*tmp_51;
      real_t tmp_54 = 0.2844444444444445*tmp_27;
      real_t tmp_55 = 0.7692346550528415*tmp_12 + tmp_13;
      real_t tmp_56 = tmp_10*tmp_55;
      real_t tmp_57 = 0.7692346550528415*tmp_16 + tmp_17;
      real_t tmp_58 = tmp_57*tmp_7;
      real_t tmp_59 = tmp_56 + tmp_58;
      real_t tmp_60 = tmp_55*tmp_9;
      real_t tmp_61 = tmp_57*tmp_8;
      real_t tmp_62 = tmp_60 + tmp_61;
      real_t tmp_63 = tmp_1*(tmp_62 - 1.0/3.0) + tmp_5*(tmp_59 - 1.0/3.0);
      real_t tmp_64 = -tmp_56 - tmp_58 - tmp_60 - tmp_61 + 1;
      real_t tmp_65 = tmp_28*tmp_63;
      real_t tmp_66 = 0.2393143352496831*tmp_27;
      real_t tmp_67 = 0.95308992296933193*tmp_12 + tmp_13;
      real_t tmp_68 = tmp_10*tmp_67;
      real_t tmp_69 = 0.95308992296933193*tmp_16 + tmp_17;
      real_t tmp_70 = tmp_69*tmp_7;
      real_t tmp_71 = tmp_68 + tmp_70;
      real_t tmp_72 = tmp_67*tmp_9;
      real_t tmp_73 = tmp_69*tmp_8;
      real_t tmp_74 = tmp_72 + tmp_73;
      real_t tmp_75 = tmp_1*(tmp_74 - 1.0/3.0) + tmp_5*(tmp_71 - 1.0/3.0);
      real_t tmp_76 = -tmp_68 - tmp_70 - tmp_72 - tmp_73 + 1;
      real_t tmp_77 = tmp_28*tmp_75;
      real_t tmp_78 = 0.11846344252809471*tmp_27;
      real_t tmp_79 = p_affine_10_0*tmp_7 + p_affine_10_1*tmp_10;
      real_t tmp_80 = p_affine_10_0*tmp_8 + p_affine_10_1*tmp_9;
      real_t a_0_0 = tmp_30*(-tmp_11*tmp_24 - tmp_25*tmp_26 + tmp_26*tmp_29) + tmp_42*(-tmp_11*tmp_39 - tmp_25*tmp_40 + tmp_40*tmp_41) + tmp_54*(-tmp_11*tmp_51 - tmp_25*tmp_52 + tmp_52*tmp_53) + tmp_66*(-tmp_11*tmp_63 - tmp_25*tmp_64 + tmp_64*tmp_65) + tmp_78*(-tmp_11*tmp_75 - tmp_25*tmp_76 + tmp_76*tmp_77);
      real_t a_0_1 = tmp_30*(-tmp_20*tmp_25 + tmp_20*tmp_29 - tmp_24*tmp_79) + tmp_42*(-tmp_25*tmp_35 + tmp_35*tmp_41 - tmp_39*tmp_79) + tmp_54*(-tmp_25*tmp_47 + tmp_47*tmp_53 - tmp_51*tmp_79) + tmp_66*(-tmp_25*tmp_59 + tmp_59*tmp_65 - tmp_63*tmp_79) + tmp_78*(-tmp_25*tmp_71 + tmp_71*tmp_77 - tmp_75*tmp_79);
      real_t a_0_2 = tmp_30*(-tmp_23*tmp_25 + tmp_23*tmp_29 - tmp_24*tmp_80) + tmp_42*(-tmp_25*tmp_38 + tmp_38*tmp_41 - tmp_39*tmp_80) + tmp_54*(-tmp_25*tmp_50 + tmp_50*tmp_53 - tmp_51*tmp_80) + tmp_66*(-tmp_25*tmp_62 + tmp_62*tmp_65 - tmp_63*tmp_80) + tmp_78*(-tmp_25*tmp_74 + tmp_74*tmp_77 - tmp_75*tmp_80);
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
      elMat( 0, 0) = 0;
   }

public:


};




class EGVectorLaplaceFormEDGEDG_new : public hyteg::dg::DGForm2D
{
 public:
    EGVectorLaplaceFormEDGEDG_new()
: callback_Scalar_Variable_Coefficient_2D_g1 ([](const Point3D & p) -> real_t { return 0.; })
, callback_Scalar_Variable_Coefficient_2D_g0 ([](const Point3D & p) -> real_t { return 0.; })
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
      real_t tmp_1 = p_affine_2_0 + tmp_0;
      real_t tmp_2 = p_affine_1_0 + tmp_0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = p_affine_1_1 + tmp_3;
      real_t tmp_7 = 1.0 / (-tmp_1*tmp_6 + tmp_5);
      real_t tmp_8 = tmp_2*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_5*tmp_7;
      real_t tmp_11 = tmp_7*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_12 = tmp_6*tmp_7;
      real_t tmp_13 = (((tmp_10 + tmp_12*tmp_9)*(tmp_10 + tmp_12*tmp_9)) + ((tmp_1*tmp_11 + tmp_10)*(tmp_1*tmp_11 + tmp_10)) + ((tmp_1*tmp_8 + tmp_8*tmp_9)*(tmp_1*tmp_8 + tmp_8*tmp_9)) + ((tmp_11*tmp_4 + tmp_12*tmp_4)*(tmp_11*tmp_4 + tmp_12*tmp_4)))*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = 0.5*tmp_13;
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
      real_t tmp_8 = tmp_4*tmp_7;
      real_t tmp_9 = p_affine_2_0 + tmp_3;
      real_t tmp_10 = p_affine_1_1 + tmp_6;
      real_t tmp_11 = 1.0 / (-tmp_10*tmp_9 + tmp_8);
      real_t tmp_12 = p_affine_6_1 + tmp_6;
      real_t tmp_13 = tmp_11*(0.046910077030668018*tmp_1 + tmp_12);
      real_t tmp_14 = p_affine_6_0 + tmp_3;
      real_t tmp_15 = tmp_11*(0.046910077030668018*tmp_0 + tmp_14);
      real_t tmp_16 = tmp_13*tmp_5 + tmp_15*tmp_7 - 1.0/3.0;
      real_t tmp_17 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_18 = tmp_13*tmp_4 + tmp_15*tmp_17 - 1.0/3.0;
      real_t tmp_19 = tmp_16*tmp_4 + tmp_18*tmp_9;
      real_t tmp_20 = tmp_11*tmp_8;
      real_t tmp_21 = tmp_11*tmp_9;
      real_t tmp_22 = tmp_11*tmp_5;
      real_t tmp_23 = 1.0*p_affine_10_0*(tmp_17*tmp_21 + tmp_20) + 1.0*p_affine_10_1*(tmp_21*tmp_4 + tmp_22*tmp_4);
      real_t tmp_24 = tmp_10*tmp_16 + tmp_18*tmp_7;
      real_t tmp_25 = tmp_11*tmp_7;
      real_t tmp_26 = 1.0*p_affine_10_0*(tmp_10*tmp_25 + tmp_17*tmp_25) + 1.0*p_affine_10_1*(tmp_10*tmp_22 + tmp_20);
      real_t tmp_27 = 6/tmp_2;
      real_t tmp_28 = 0.23076534494715845*tmp_1 + tmp_12;
      real_t tmp_29 = 0.23076534494715845*tmp_0 + tmp_14;
      real_t tmp_30 = tmp_22*tmp_28 + tmp_25*tmp_29 - 1.0/3.0;
      real_t tmp_31 = tmp_11*tmp_4;
      real_t tmp_32 = tmp_11*tmp_17;
      real_t tmp_33 = tmp_28*tmp_31 + tmp_29*tmp_32 - 1.0/3.0;
      real_t tmp_34 = tmp_30*tmp_4 + tmp_33*tmp_9;
      real_t tmp_35 = tmp_10*tmp_30 + tmp_33*tmp_7;
      real_t tmp_36 = 0.5*tmp_1 + tmp_12;
      real_t tmp_37 = 0.5*tmp_0 + tmp_14;
      real_t tmp_38 = tmp_22*tmp_36 + tmp_25*tmp_37 - 1.0/3.0;
      real_t tmp_39 = tmp_31*tmp_36 + tmp_32*tmp_37 - 1.0/3.0;
      real_t tmp_40 = tmp_38*tmp_4 + tmp_39*tmp_9;
      real_t tmp_41 = tmp_10*tmp_38 + tmp_39*tmp_7;
      real_t tmp_42 = 0.7692346550528415*tmp_1 + tmp_12;
      real_t tmp_43 = 0.7692346550528415*tmp_0 + tmp_14;
      real_t tmp_44 = tmp_22*tmp_42 + tmp_25*tmp_43 - 1.0/3.0;
      real_t tmp_45 = tmp_31*tmp_42 + tmp_32*tmp_43 - 1.0/3.0;
      real_t tmp_46 = tmp_4*tmp_44 + tmp_45*tmp_9;
      real_t tmp_47 = tmp_10*tmp_44 + tmp_45*tmp_7;
      real_t tmp_48 = 0.95308992296933193*tmp_1 + tmp_12;
      real_t tmp_49 = 0.95308992296933193*tmp_0 + tmp_14;
      real_t tmp_50 = tmp_22*tmp_48 + tmp_25*tmp_49 - 1.0/3.0;
      real_t tmp_51 = tmp_31*tmp_48 + tmp_32*tmp_49 - 1.0/3.0;
      real_t tmp_52 = tmp_4*tmp_50 + tmp_51*tmp_9;
      real_t tmp_53 = tmp_10*tmp_50 + tmp_51*tmp_7;
      real_t a_0_0 = 0.11846344252809471*tmp_2*(-tmp_19*tmp_23 - tmp_24*tmp_26 + tmp_27*((tmp_19*tmp_19) + (tmp_24*tmp_24))) + 0.2393143352496831*tmp_2*(-tmp_23*tmp_34 - tmp_26*tmp_35 + tmp_27*((tmp_34*tmp_34) + (tmp_35*tmp_35))) + 0.2844444444444445*tmp_2*(-tmp_23*tmp_40 - tmp_26*tmp_41 + tmp_27*((tmp_40*tmp_40) + (tmp_41*tmp_41))) + 0.2393143352496831*tmp_2*(-tmp_23*tmp_46 - tmp_26*tmp_47 + tmp_27*((tmp_46*tmp_46) + (tmp_47*tmp_47))) + 0.11846344252809471*tmp_2*(-tmp_23*tmp_52 - tmp_26*tmp_53 + tmp_27*((tmp_52*tmp_52) + (tmp_53*tmp_53)));
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
      real_t tmp_8 = tmp_4*tmp_7;
      real_t tmp_9 = p_affine_5_0 + tmp_3;
      real_t tmp_10 = p_affine_4_1 + tmp_6;
      real_t tmp_11 = 1.0 / (-tmp_10*tmp_9 + tmp_8);
      real_t tmp_12 = p_affine_6_1 + 0.046910077030668018*tmp_1;
      real_t tmp_13 = tmp_11*(tmp_12 + tmp_6);
      real_t tmp_14 = p_affine_6_0 + 0.046910077030668018*tmp_0;
      real_t tmp_15 = tmp_11*(tmp_14 + tmp_3);
      real_t tmp_16 = tmp_13*tmp_5 + tmp_15*tmp_7 - 1.0/3.0;
      real_t tmp_17 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_18 = tmp_13*tmp_4 + tmp_15*tmp_17 - 1.0/3.0;
      real_t tmp_19 = tmp_16*tmp_4 + tmp_18*tmp_9;
      real_t tmp_20 = -p_affine_0_0;
      real_t tmp_21 = p_affine_1_0 + tmp_20;
      real_t tmp_22 = -p_affine_0_1;
      real_t tmp_23 = p_affine_2_1 + tmp_22;
      real_t tmp_24 = tmp_21*tmp_23;
      real_t tmp_25 = p_affine_2_0 + tmp_20;
      real_t tmp_26 = p_affine_1_1 + tmp_22;
      real_t tmp_27 = 1.0 / (tmp_24 - tmp_25*tmp_26);
      real_t tmp_28 = tmp_24*tmp_27;
      real_t tmp_29 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_30 = tmp_25*tmp_27;
      real_t tmp_31 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = 0.5*p_affine_10_0*(tmp_28 + tmp_29*tmp_30) + 0.5*p_affine_10_1*(tmp_21*tmp_30 + tmp_21*tmp_32);
      real_t tmp_34 = tmp_10*tmp_16 + tmp_18*tmp_7;
      real_t tmp_35 = tmp_23*tmp_27;
      real_t tmp_36 = 0.5*p_affine_10_0*(tmp_26*tmp_35 + tmp_29*tmp_35) + 0.5*p_affine_10_1*(tmp_26*tmp_32 + tmp_28);
      real_t tmp_37 = tmp_27*(tmp_12 + tmp_22);
      real_t tmp_38 = tmp_27*(tmp_14 + tmp_20);
      real_t tmp_39 = tmp_23*tmp_38 + tmp_31*tmp_37 - 1.0/3.0;
      real_t tmp_40 = tmp_21*tmp_37 + tmp_29*tmp_38 - 1.0/3.0;
      real_t tmp_41 = tmp_21*tmp_39 + tmp_25*tmp_40;
      real_t tmp_42 = tmp_11*tmp_8;
      real_t tmp_43 = tmp_11*tmp_9;
      real_t tmp_44 = tmp_11*tmp_5;
      real_t tmp_45 = 0.5*p_affine_10_0*(tmp_17*tmp_43 + tmp_42) + 0.5*p_affine_10_1*(tmp_4*tmp_43 + tmp_4*tmp_44);
      real_t tmp_46 = tmp_23*tmp_40 + tmp_26*tmp_39;
      real_t tmp_47 = tmp_11*tmp_7;
      real_t tmp_48 = 0.5*p_affine_10_0*(tmp_10*tmp_47 + tmp_17*tmp_47) + 0.5*p_affine_10_1*(tmp_10*tmp_44 + tmp_42);
      real_t tmp_49 = 6/tmp_2;
      real_t tmp_50 = p_affine_6_1 + 0.23076534494715845*tmp_1;
      real_t tmp_51 = tmp_50 + tmp_6;
      real_t tmp_52 = p_affine_6_0 + 0.23076534494715845*tmp_0;
      real_t tmp_53 = tmp_3 + tmp_52;
      real_t tmp_54 = tmp_44*tmp_51 + tmp_47*tmp_53 - 1.0/3.0;
      real_t tmp_55 = tmp_11*tmp_4;
      real_t tmp_56 = tmp_11*tmp_17;
      real_t tmp_57 = tmp_51*tmp_55 + tmp_53*tmp_56 - 1.0/3.0;
      real_t tmp_58 = tmp_4*tmp_54 + tmp_57*tmp_9;
      real_t tmp_59 = tmp_10*tmp_54 + tmp_57*tmp_7;
      real_t tmp_60 = tmp_22 + tmp_50;
      real_t tmp_61 = tmp_20 + tmp_52;
      real_t tmp_62 = tmp_32*tmp_60 + tmp_35*tmp_61 - 1.0/3.0;
      real_t tmp_63 = tmp_21*tmp_27;
      real_t tmp_64 = tmp_27*tmp_29;
      real_t tmp_65 = tmp_60*tmp_63 + tmp_61*tmp_64 - 1.0/3.0;
      real_t tmp_66 = tmp_21*tmp_62 + tmp_25*tmp_65;
      real_t tmp_67 = tmp_23*tmp_65 + tmp_26*tmp_62;
      real_t tmp_68 = p_affine_6_1 + 0.5*tmp_1;
      real_t tmp_69 = tmp_6 + tmp_68;
      real_t tmp_70 = p_affine_6_0 + 0.5*tmp_0;
      real_t tmp_71 = tmp_3 + tmp_70;
      real_t tmp_72 = tmp_44*tmp_69 + tmp_47*tmp_71 - 1.0/3.0;
      real_t tmp_73 = tmp_55*tmp_69 + tmp_56*tmp_71 - 1.0/3.0;
      real_t tmp_74 = tmp_4*tmp_72 + tmp_73*tmp_9;
      real_t tmp_75 = tmp_10*tmp_72 + tmp_7*tmp_73;
      real_t tmp_76 = tmp_22 + tmp_68;
      real_t tmp_77 = tmp_20 + tmp_70;
      real_t tmp_78 = tmp_32*tmp_76 + tmp_35*tmp_77 - 1.0/3.0;
      real_t tmp_79 = tmp_63*tmp_76 + tmp_64*tmp_77 - 1.0/3.0;
      real_t tmp_80 = tmp_21*tmp_78 + tmp_25*tmp_79;
      real_t tmp_81 = tmp_23*tmp_79 + tmp_26*tmp_78;
      real_t tmp_82 = p_affine_6_1 + 0.7692346550528415*tmp_1;
      real_t tmp_83 = tmp_6 + tmp_82;
      real_t tmp_84 = p_affine_6_0 + 0.7692346550528415*tmp_0;
      real_t tmp_85 = tmp_3 + tmp_84;
      real_t tmp_86 = tmp_44*tmp_83 + tmp_47*tmp_85 - 1.0/3.0;
      real_t tmp_87 = tmp_55*tmp_83 + tmp_56*tmp_85 - 1.0/3.0;
      real_t tmp_88 = tmp_4*tmp_86 + tmp_87*tmp_9;
      real_t tmp_89 = tmp_10*tmp_86 + tmp_7*tmp_87;
      real_t tmp_90 = tmp_22 + tmp_82;
      real_t tmp_91 = tmp_20 + tmp_84;
      real_t tmp_92 = tmp_32*tmp_90 + tmp_35*tmp_91 - 1.0/3.0;
      real_t tmp_93 = tmp_63*tmp_90 + tmp_64*tmp_91 - 1.0/3.0;
      real_t tmp_94 = tmp_21*tmp_92 + tmp_25*tmp_93;
      real_t tmp_95 = tmp_23*tmp_93 + tmp_26*tmp_92;
      real_t tmp_96 = p_affine_6_1 + 0.95308992296933193*tmp_1;
      real_t tmp_97 = tmp_6 + tmp_96;
      real_t tmp_98 = p_affine_6_0 + 0.95308992296933193*tmp_0;
      real_t tmp_99 = tmp_3 + tmp_98;
      real_t tmp_100 = tmp_44*tmp_97 + tmp_47*tmp_99 - 1.0/3.0;
      real_t tmp_101 = tmp_55*tmp_97 + tmp_56*tmp_99 - 1.0/3.0;
      real_t tmp_102 = tmp_100*tmp_4 + tmp_101*tmp_9;
      real_t tmp_103 = tmp_10*tmp_100 + tmp_101*tmp_7;
      real_t tmp_104 = tmp_22 + tmp_96;
      real_t tmp_105 = tmp_20 + tmp_98;
      real_t tmp_106 = tmp_104*tmp_32 + tmp_105*tmp_35 - 1.0/3.0;
      real_t tmp_107 = tmp_104*tmp_63 + tmp_105*tmp_64 - 1.0/3.0;
      real_t tmp_108 = tmp_106*tmp_21 + tmp_107*tmp_25;
      real_t tmp_109 = tmp_106*tmp_26 + tmp_107*tmp_23;
      real_t a_0_0 = 0.11846344252809471*tmp_2*(tmp_102*tmp_33 + tmp_103*tmp_36 - tmp_108*tmp_45 - tmp_109*tmp_48 - tmp_49*(tmp_102*tmp_108 + tmp_103*tmp_109)) + 0.11846344252809471*tmp_2*(tmp_19*tmp_33 + tmp_34*tmp_36 - tmp_41*tmp_45 - tmp_46*tmp_48 - tmp_49*(tmp_19*tmp_41 + tmp_34*tmp_46)) + 0.2393143352496831*tmp_2*(tmp_33*tmp_58 + tmp_36*tmp_59 - tmp_45*tmp_66 - tmp_48*tmp_67 - tmp_49*(tmp_58*tmp_66 + tmp_59*tmp_67)) + 0.2844444444444445*tmp_2*(tmp_33*tmp_74 + tmp_36*tmp_75 - tmp_45*tmp_80 - tmp_48*tmp_81 - tmp_49*(tmp_74*tmp_80 + tmp_75*tmp_81)) + 0.2393143352496831*tmp_2*(tmp_33*tmp_88 + tmp_36*tmp_89 - tmp_45*tmp_94 - tmp_48*tmp_95 - tmp_49*(tmp_88*tmp_94 + tmp_89*tmp_95));
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
      real_t tmp_5 = -p_affine_0_1;
      real_t tmp_6 = p_affine_6_1 + tmp_5;
      real_t tmp_7 = 0.046910077030668018*tmp_1 + tmp_6;
      real_t tmp_8 = p_affine_2_1 + tmp_5;
      real_t tmp_9 = tmp_4*tmp_8;
      real_t tmp_10 = p_affine_2_0 + tmp_3;
      real_t tmp_11 = p_affine_1_1 + tmp_5;
      real_t tmp_12 = 1.0 / (-tmp_10*tmp_11 + tmp_9);
      real_t tmp_13 = tmp_12*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_14 = p_affine_6_0 + tmp_3;
      real_t tmp_15 = tmp_12*(0.046910077030668018*tmp_0 + tmp_14);
      real_t tmp_16 = tmp_13*tmp_7 + tmp_15*tmp_8 - 1.0/3.0;
      real_t tmp_17 = tmp_12*tmp_4;
      real_t tmp_18 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_19 = tmp_15*tmp_18 + tmp_17*tmp_7 - 1.0/3.0;
      real_t tmp_20 = tmp_10*tmp_19 + tmp_16*tmp_4;
      real_t tmp_21 = tmp_12*tmp_9;
      real_t tmp_22 = tmp_10*tmp_12;
      real_t tmp_23 = 2*p_affine_10_0*(tmp_18*tmp_22 + tmp_21) + 2*p_affine_10_1*(tmp_13*tmp_4 + tmp_22*tmp_4);
      real_t tmp_24 = tmp_11*tmp_16 + tmp_19*tmp_8;
      real_t tmp_25 = tmp_12*tmp_8;
      real_t tmp_26 = 2*p_affine_10_0*(tmp_11*tmp_25 + tmp_18*tmp_25) + 2*p_affine_10_1*(tmp_11*tmp_13 + tmp_21);
      real_t tmp_27 = 6/tmp_2;
      real_t tmp_28 = 0.23076534494715845*tmp_1 + tmp_6;
      real_t tmp_29 = 0.23076534494715845*tmp_0 + tmp_14;
      real_t tmp_30 = tmp_13*tmp_28 + tmp_25*tmp_29 - 1.0/3.0;
      real_t tmp_31 = tmp_12*tmp_18;
      real_t tmp_32 = tmp_17*tmp_28 + tmp_29*tmp_31 - 1.0/3.0;
      real_t tmp_33 = tmp_10*tmp_32 + tmp_30*tmp_4;
      real_t tmp_34 = tmp_11*tmp_30 + tmp_32*tmp_8;
      real_t tmp_35 = 0.5*tmp_1 + tmp_6;
      real_t tmp_36 = 0.5*tmp_0 + tmp_14;
      real_t tmp_37 = tmp_13*tmp_35 + tmp_25*tmp_36 - 1.0/3.0;
      real_t tmp_38 = tmp_17*tmp_35 + tmp_31*tmp_36 - 1.0/3.0;
      real_t tmp_39 = tmp_10*tmp_38 + tmp_37*tmp_4;
      real_t tmp_40 = tmp_11*tmp_37 + tmp_38*tmp_8;
      real_t tmp_41 = 0.7692346550528415*tmp_1 + tmp_6;
      real_t tmp_42 = 0.7692346550528415*tmp_0 + tmp_14;
      real_t tmp_43 = tmp_13*tmp_41 + tmp_25*tmp_42 - 1.0/3.0;
      real_t tmp_44 = tmp_17*tmp_41 + tmp_31*tmp_42 - 1.0/3.0;
      real_t tmp_45 = tmp_10*tmp_44 + tmp_4*tmp_43;
      real_t tmp_46 = tmp_11*tmp_43 + tmp_44*tmp_8;
      real_t tmp_47 = 0.95308992296933193*tmp_1 + tmp_6;
      real_t tmp_48 = 0.95308992296933193*tmp_0 + tmp_14;
      real_t tmp_49 = tmp_13*tmp_47 + tmp_25*tmp_48 - 1.0/3.0;
      real_t tmp_50 = tmp_17*tmp_47 + tmp_31*tmp_48 - 1.0/3.0;
      real_t tmp_51 = tmp_10*tmp_50 + tmp_4*tmp_49;
      real_t tmp_52 = tmp_11*tmp_49 + tmp_50*tmp_8;
      real_t a_0_0 = 0.11846344252809471*tmp_2*(-tmp_20*tmp_23 - tmp_24*tmp_26 + tmp_27*((tmp_20*tmp_20) + (tmp_24*tmp_24))) + 0.2393143352496831*tmp_2*(-tmp_23*tmp_33 - tmp_26*tmp_34 + tmp_27*((tmp_33*tmp_33) + (tmp_34*tmp_34))) + 0.2844444444444445*tmp_2*(-tmp_23*tmp_39 - tmp_26*tmp_40 + tmp_27*((tmp_39*tmp_39) + (tmp_40*tmp_40))) + 0.2393143352496831*tmp_2*(-tmp_23*tmp_45 - tmp_26*tmp_46 + tmp_27*((tmp_45*tmp_45) + (tmp_46*tmp_46))) + 0.11846344252809471*tmp_2*(-tmp_23*tmp_51 - tmp_26*tmp_52 + tmp_27*((tmp_51*tmp_51) + (tmp_52*tmp_52)));
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
      real_t tmp_0 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_1*tmp_1), 1.0/2.0));
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = p_affine_1_0 + tmp_3;
      real_t tmp_5 = -p_affine_0_1;
      real_t tmp_6 = p_affine_6_1 + tmp_5;
      real_t tmp_7 = 0.046910077030668018*tmp_1 + tmp_6;
      real_t tmp_8 = p_affine_2_1 + tmp_5;
      real_t tmp_9 = tmp_4*tmp_8;
      real_t tmp_10 = p_affine_2_0 + tmp_3;
      real_t tmp_11 = p_affine_1_1 + tmp_5;
      real_t tmp_12 = 1.0 / (-tmp_10*tmp_11 + tmp_9);
      real_t tmp_13 = tmp_12*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_14 = p_affine_6_0 + tmp_3;
      real_t tmp_15 = tmp_12*(0.046910077030668018*tmp_0 + tmp_14);
      real_t tmp_16 = tmp_13*tmp_7 + tmp_15*tmp_8 - 1.0/3.0;
      real_t tmp_17 = tmp_12*tmp_4;
      real_t tmp_18 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_19 = tmp_15*tmp_18 + tmp_17*tmp_7 - 1.0/3.0;
      real_t tmp_20 = 6/tmp_2;
      real_t tmp_21 = tmp_12*tmp_9;
      real_t tmp_22 = tmp_10*tmp_12;
      real_t tmp_23 = -p_affine_10_0*(tmp_18*tmp_22 + tmp_21) - p_affine_10_1*(tmp_13*tmp_4 + tmp_22*tmp_4);
      real_t tmp_24 = tmp_12*tmp_8;
      real_t tmp_25 = -p_affine_10_0*(tmp_11*tmp_24 + tmp_18*tmp_24) - p_affine_10_1*(tmp_11*tmp_13 + tmp_21);
      real_t tmp_26 = 0.23076534494715845*tmp_1 + tmp_6;
      real_t tmp_27 = 0.23076534494715845*tmp_0 + tmp_14;
      real_t tmp_28 = tmp_13*tmp_26 + tmp_24*tmp_27 - 1.0/3.0;
      real_t tmp_29 = tmp_12*tmp_18;
      real_t tmp_30 = tmp_17*tmp_26 + tmp_27*tmp_29 - 1.0/3.0;
      real_t tmp_31 = 0.5*tmp_1 + tmp_6;
      real_t tmp_32 = 0.5*tmp_0 + tmp_14;
      real_t tmp_33 = tmp_13*tmp_31 + tmp_24*tmp_32 - 1.0/3.0;
      real_t tmp_34 = tmp_17*tmp_31 + tmp_29*tmp_32 - 1.0/3.0;
      real_t tmp_35 = 0.7692346550528415*tmp_1 + tmp_6;
      real_t tmp_36 = 0.7692346550528415*tmp_0 + tmp_14;
      real_t tmp_37 = tmp_13*tmp_35 + tmp_24*tmp_36 - 1.0/3.0;
      real_t tmp_38 = tmp_17*tmp_35 + tmp_29*tmp_36 - 1.0/3.0;
      real_t tmp_39 = 0.95308992296933193*tmp_1 + tmp_6;
      real_t tmp_40 = 0.95308992296933193*tmp_0 + tmp_14;
      real_t tmp_41 = tmp_13*tmp_39 + tmp_24*tmp_40 - 1.0/3.0;
      real_t tmp_42 = tmp_17*tmp_39 + tmp_29*tmp_40 - 1.0/3.0;
      real_t a_0_0 = 0.11846344252809471*tmp_2*(Scalar_Variable_Coefficient_2D_g0_out0_id0*(tmp_20*(tmp_10*tmp_19 + tmp_16*tmp_4) + tmp_23) + Scalar_Variable_Coefficient_2D_g1_out0_id1*(tmp_20*(tmp_11*tmp_16 + tmp_19*tmp_8) + tmp_25)) + 0.2393143352496831*tmp_2*(Scalar_Variable_Coefficient_2D_g0_out0_id2*(tmp_20*(tmp_10*tmp_30 + tmp_28*tmp_4) + tmp_23) + Scalar_Variable_Coefficient_2D_g1_out0_id3*(tmp_20*(tmp_11*tmp_28 + tmp_30*tmp_8) + tmp_25)) + 0.2844444444444445*tmp_2*(Scalar_Variable_Coefficient_2D_g0_out0_id4*(tmp_20*(tmp_10*tmp_34 + tmp_33*tmp_4) + tmp_23) + Scalar_Variable_Coefficient_2D_g1_out0_id5*(tmp_20*(tmp_11*tmp_33 + tmp_34*tmp_8) + tmp_25)) + 0.2393143352496831*tmp_2*(Scalar_Variable_Coefficient_2D_g0_out0_id6*(tmp_20*(tmp_10*tmp_38 + tmp_37*tmp_4) + tmp_23) + Scalar_Variable_Coefficient_2D_g1_out0_id7*(tmp_20*(tmp_11*tmp_37 + tmp_38*tmp_8) + tmp_25)) + 0.11846344252809471*tmp_2*(Scalar_Variable_Coefficient_2D_g0_out0_id8*(tmp_20*(tmp_10*tmp_42 + tmp_4*tmp_41) + tmp_23) + Scalar_Variable_Coefficient_2D_g1_out0_id9*(tmp_20*(tmp_11*tmp_41 + tmp_42*tmp_8) + tmp_25));
      elMat( 0, 0) = a_0_0;
   }

public:

std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g0;
std::function< real_t ( const Point3D & ) > callback_Scalar_Variable_Coefficient_2D_g1;
};


} // dg
} // hyteg

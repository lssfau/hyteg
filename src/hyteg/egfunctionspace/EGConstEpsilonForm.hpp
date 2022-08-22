
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

namespace eg{    
class EGConstEpsilonFormEP1_0 : public hyteg::dg::DGForm2D
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

      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_0;
      real_t tmp_3 = p_affine_1_0 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_2;
      real_t tmp_6 = p_affine_1_1 + tmp_0;
      real_t tmp_7 = 1.0 / (tmp_4 - tmp_5*tmp_6);
      real_t tmp_8 = 2.0*tmp_7;
      real_t tmp_9 = tmp_1*tmp_8;
      real_t tmp_10 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_11 = tmp_10*tmp_8;
      real_t tmp_12 = 1.0*tmp_7;
      real_t tmp_13 = tmp_10*tmp_12*tmp_5 + tmp_12*tmp_4;
      real_t tmp_14 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_15 = 0.5*tmp_7;
      real_t tmp_16 = tmp_15*tmp_3;
      real_t tmp_17 = tmp_1*tmp_15;
      real_t tmp_18 = tmp_10*tmp_17 + tmp_14*tmp_16 + tmp_16*tmp_5 + tmp_17*tmp_6;
      real_t tmp_19 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_20 = tmp_19*(tmp_13*(-tmp_11 - tmp_9) + 2*tmp_18*(-tmp_12*tmp_14 - tmp_12*tmp_3));
      real_t tmp_21 = tmp_18*tmp_8;
      real_t tmp_22 = tmp_19*(tmp_13*tmp_9 + tmp_14*tmp_21);
      real_t tmp_23 = tmp_19*(tmp_11*tmp_13 + tmp_21*tmp_3);
      real_t a_0_0 = 0.5*tmp_20;
      real_t a_0_1 = 0.5*tmp_22;
      real_t a_0_2 = 0.5*tmp_23;
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
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = p_affine_1_0 + tmp_3;
      real_t tmp_5 = p_affine_2_1 + tmp_0;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = p_affine_2_0 + tmp_3;
      real_t tmp_8 = 1.0 / (-tmp_1*tmp_7 + tmp_6);
      real_t tmp_9 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + tmp_0;
      real_t tmp_11 = tmp_8*(tmp_10 + 0.046910077030668018*tmp_9);
      real_t tmp_12 = tmp_11*tmp_2;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_3;
      real_t tmp_15 = tmp_8*(0.046910077030668018*tmp_13 + tmp_14);
      real_t tmp_16 = tmp_15*tmp_5;
      real_t tmp_17 = tmp_12 + tmp_16;
      real_t tmp_18 = tmp_17 - 1.0/3.0;
      real_t tmp_19 = tmp_11*tmp_4;
      real_t tmp_20 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_21 = tmp_15*tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_22 - 1.0/3.0;
      real_t tmp_24 = p_affine_10_0*(tmp_1*tmp_18 + tmp_23*tmp_5);
      real_t tmp_25 = 0.5*tmp_8;
      real_t tmp_26 = tmp_25*tmp_4;
      real_t tmp_27 = tmp_2*tmp_25;
      real_t tmp_28 = -tmp_26 - tmp_27;
      real_t tmp_29 = 1.0*tmp_28;
      real_t tmp_30 = tmp_18*tmp_4 + tmp_23*tmp_7;
      real_t tmp_31 = 1.0*tmp_8;
      real_t tmp_32 = tmp_31*tmp_5;
      real_t tmp_33 = tmp_20*tmp_31;
      real_t tmp_34 = 1.0*p_affine_10_0*(-tmp_32 - tmp_33) + 1.0*p_affine_10_1*tmp_28;
      real_t tmp_35 = -tmp_12 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_36 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_37 = 12/tmp_36;
      real_t tmp_38 = tmp_30*tmp_37;
      real_t tmp_39 = tmp_25*tmp_5;
      real_t tmp_40 = 1.0*p_affine_10_0*(tmp_31*tmp_6 + tmp_33*tmp_7) + 1.0*p_affine_10_1*(tmp_1*tmp_39 + tmp_20*tmp_39 + tmp_26*tmp_7 + tmp_27*tmp_4);
      real_t tmp_41 = 0.11846344252809471*tmp_36;
      real_t tmp_42 = tmp_8*(tmp_10 + 0.23076534494715845*tmp_9);
      real_t tmp_43 = tmp_2*tmp_42;
      real_t tmp_44 = tmp_8*(0.23076534494715845*tmp_13 + tmp_14);
      real_t tmp_45 = tmp_44*tmp_5;
      real_t tmp_46 = tmp_43 + tmp_45;
      real_t tmp_47 = tmp_46 - 1.0/3.0;
      real_t tmp_48 = tmp_4*tmp_42;
      real_t tmp_49 = tmp_20*tmp_44;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_50 - 1.0/3.0;
      real_t tmp_52 = tmp_1*tmp_47 + tmp_5*tmp_51;
      real_t tmp_53 = p_affine_10_0*tmp_29;
      real_t tmp_54 = tmp_4*tmp_47 + tmp_51*tmp_7;
      real_t tmp_55 = -tmp_43 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_56 = tmp_37*tmp_54;
      real_t tmp_57 = 0.2393143352496831*tmp_36;
      real_t tmp_58 = tmp_8*(tmp_10 + 0.5*tmp_9);
      real_t tmp_59 = tmp_2*tmp_58;
      real_t tmp_60 = tmp_8*(0.5*tmp_13 + tmp_14);
      real_t tmp_61 = tmp_5*tmp_60;
      real_t tmp_62 = tmp_59 + tmp_61;
      real_t tmp_63 = tmp_62 - 1.0/3.0;
      real_t tmp_64 = tmp_4*tmp_58;
      real_t tmp_65 = tmp_20*tmp_60;
      real_t tmp_66 = tmp_64 + tmp_65;
      real_t tmp_67 = tmp_66 - 1.0/3.0;
      real_t tmp_68 = tmp_1*tmp_63 + tmp_5*tmp_67;
      real_t tmp_69 = tmp_4*tmp_63 + tmp_67*tmp_7;
      real_t tmp_70 = -tmp_59 - tmp_61 - tmp_64 - tmp_65 + 1;
      real_t tmp_71 = tmp_37*tmp_69;
      real_t tmp_72 = 0.2844444444444445*tmp_36;
      real_t tmp_73 = tmp_8*(tmp_10 + 0.7692346550528415*tmp_9);
      real_t tmp_74 = tmp_2*tmp_73;
      real_t tmp_75 = tmp_8*(0.7692346550528415*tmp_13 + tmp_14);
      real_t tmp_76 = tmp_5*tmp_75;
      real_t tmp_77 = tmp_74 + tmp_76;
      real_t tmp_78 = tmp_77 - 1.0/3.0;
      real_t tmp_79 = tmp_4*tmp_73;
      real_t tmp_80 = tmp_20*tmp_75;
      real_t tmp_81 = tmp_79 + tmp_80;
      real_t tmp_82 = tmp_81 - 1.0/3.0;
      real_t tmp_83 = tmp_1*tmp_78 + tmp_5*tmp_82;
      real_t tmp_84 = tmp_4*tmp_78 + tmp_7*tmp_82;
      real_t tmp_85 = -tmp_74 - tmp_76 - tmp_79 - tmp_80 + 1;
      real_t tmp_86 = tmp_37*tmp_84;
      real_t tmp_87 = 0.2393143352496831*tmp_36;
      real_t tmp_88 = tmp_8*(tmp_10 + 0.95308992296933193*tmp_9);
      real_t tmp_89 = tmp_2*tmp_88;
      real_t tmp_90 = tmp_8*(0.95308992296933193*tmp_13 + tmp_14);
      real_t tmp_91 = tmp_5*tmp_90;
      real_t tmp_92 = tmp_89 + tmp_91;
      real_t tmp_93 = tmp_92 - 1.0/3.0;
      real_t tmp_94 = tmp_4*tmp_88;
      real_t tmp_95 = tmp_20*tmp_90;
      real_t tmp_96 = tmp_94 + tmp_95;
      real_t tmp_97 = tmp_96 - 1.0/3.0;
      real_t tmp_98 = tmp_1*tmp_93 + tmp_5*tmp_97;
      real_t tmp_99 = tmp_4*tmp_93 + tmp_7*tmp_97;
      real_t tmp_100 = -tmp_89 - tmp_91 - tmp_94 - tmp_95 + 1;
      real_t tmp_101 = tmp_37*tmp_99;
      real_t tmp_102 = 0.11846344252809471*tmp_36;
      real_t tmp_103 = 1.0*p_affine_10_0*tmp_32 + 1.0*p_affine_10_1*tmp_27;
      real_t tmp_104 = p_affine_10_0*tmp_27;
      real_t tmp_105 = 1.0*p_affine_10_0*tmp_33 + 1.0*p_affine_10_1*tmp_26;
      real_t tmp_106 = p_affine_10_0*tmp_26;
      real_t a_0_0 = tmp_102*(tmp_100*tmp_101 - tmp_100*tmp_40 - tmp_34*tmp_99 - tmp_53*tmp_98) + tmp_41*(-tmp_24*tmp_29 - tmp_30*tmp_34 + tmp_35*tmp_38 - tmp_35*tmp_40) + tmp_57*(-tmp_34*tmp_54 - tmp_40*tmp_55 - tmp_52*tmp_53 + tmp_55*tmp_56) + tmp_72*(-tmp_34*tmp_69 - tmp_40*tmp_70 - tmp_53*tmp_68 + tmp_70*tmp_71) + tmp_87*(-tmp_34*tmp_84 - tmp_40*tmp_85 - tmp_53*tmp_83 + tmp_85*tmp_86);
      real_t a_0_1 = tmp_102*(tmp_101*tmp_92 - tmp_103*tmp_99 - tmp_104*tmp_98 - tmp_40*tmp_92) + tmp_41*(-tmp_103*tmp_30 + tmp_17*tmp_38 - tmp_17*tmp_40 - tmp_24*tmp_27) + tmp_57*(-tmp_103*tmp_54 - tmp_104*tmp_52 - tmp_40*tmp_46 + tmp_46*tmp_56) + tmp_72*(-tmp_103*tmp_69 - tmp_104*tmp_68 - tmp_40*tmp_62 + tmp_62*tmp_71) + tmp_87*(-tmp_103*tmp_84 - tmp_104*tmp_83 - tmp_40*tmp_77 + tmp_77*tmp_86);
      real_t a_0_2 = tmp_102*(tmp_101*tmp_96 - tmp_105*tmp_99 - tmp_106*tmp_98 - tmp_40*tmp_96) + tmp_41*(-tmp_105*tmp_30 + tmp_22*tmp_38 - tmp_22*tmp_40 - tmp_24*tmp_26) + tmp_57*(-tmp_105*tmp_54 - tmp_106*tmp_52 - tmp_40*tmp_50 + tmp_50*tmp_56) + tmp_72*(-tmp_105*tmp_69 - tmp_106*tmp_68 - tmp_40*tmp_66 + tmp_66*tmp_71) + tmp_87*(-tmp_105*tmp_84 - tmp_106*tmp_83 - tmp_40*tmp_81 + tmp_81*tmp_86);
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
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = p_affine_1_0 + tmp_3;
      real_t tmp_5 = p_affine_2_1 + tmp_0;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = p_affine_2_0 + tmp_3;
      real_t tmp_8 = 1.0 / (-tmp_1*tmp_7 + tmp_6);
      real_t tmp_9 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + 0.046910077030668018*tmp_9;
      real_t tmp_11 = tmp_8*(tmp_0 + tmp_10);
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + 0.046910077030668018*tmp_12;
      real_t tmp_14 = tmp_8*(tmp_13 + tmp_3);
      real_t tmp_15 = tmp_11*tmp_2 + tmp_14*tmp_5 - 1.0/3.0;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_11*tmp_4 + tmp_14*tmp_16 - 1.0/3.0;
      real_t tmp_18 = p_affine_10_0*(tmp_1*tmp_15 + tmp_17*tmp_5);
      real_t tmp_19 = -p_affine_3_0;
      real_t tmp_20 = p_affine_4_0 + tmp_19;
      real_t tmp_21 = -p_affine_3_1;
      real_t tmp_22 = p_affine_5_1 + tmp_21;
      real_t tmp_23 = 1.0 / (tmp_20*tmp_22 - (p_affine_4_1 + tmp_21)*(p_affine_5_0 + tmp_19));
      real_t tmp_24 = 0.5*tmp_23;
      real_t tmp_25 = tmp_20*tmp_24;
      real_t tmp_26 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_25 - tmp_27;
      real_t tmp_29 = 1.0*tmp_28;
      real_t tmp_30 = tmp_15*tmp_4 + tmp_17*tmp_7;
      real_t tmp_31 = 1.0*tmp_23;
      real_t tmp_32 = tmp_22*tmp_31;
      real_t tmp_33 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_34 = tmp_31*tmp_33;
      real_t tmp_35 = 1.0*p_affine_10_0*(-tmp_32 - tmp_34) + 1.0*p_affine_10_1*tmp_28;
      real_t tmp_36 = tmp_23*(tmp_10 + tmp_21);
      real_t tmp_37 = tmp_20*tmp_36;
      real_t tmp_38 = tmp_26*tmp_36;
      real_t tmp_39 = tmp_23*(tmp_13 + tmp_19);
      real_t tmp_40 = tmp_22*tmp_39;
      real_t tmp_41 = tmp_33*tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_44 = 12/tmp_43;
      real_t tmp_45 = tmp_30*tmp_44;
      real_t tmp_46 = 1.0*tmp_8;
      real_t tmp_47 = 0.5*tmp_8;
      real_t tmp_48 = tmp_4*tmp_47;
      real_t tmp_49 = tmp_47*tmp_5;
      real_t tmp_50 = p_affine_10_0*(tmp_16*tmp_46*tmp_7 + tmp_46*tmp_6) + p_affine_10_1*(tmp_1*tmp_49 + tmp_16*tmp_49 + tmp_2*tmp_48 + tmp_48*tmp_7);
      real_t tmp_51 = 0.11846344252809471*tmp_43;
      real_t tmp_52 = p_affine_6_1 + 0.23076534494715845*tmp_9;
      real_t tmp_53 = tmp_8*(tmp_0 + tmp_52);
      real_t tmp_54 = p_affine_6_0 + 0.23076534494715845*tmp_12;
      real_t tmp_55 = tmp_8*(tmp_3 + tmp_54);
      real_t tmp_56 = tmp_2*tmp_53 + tmp_5*tmp_55 - 1.0/3.0;
      real_t tmp_57 = tmp_16*tmp_55 + tmp_4*tmp_53 - 1.0/3.0;
      real_t tmp_58 = tmp_1*tmp_56 + tmp_5*tmp_57;
      real_t tmp_59 = p_affine_10_0*tmp_29;
      real_t tmp_60 = tmp_4*tmp_56 + tmp_57*tmp_7;
      real_t tmp_61 = tmp_23*(tmp_21 + tmp_52);
      real_t tmp_62 = tmp_20*tmp_61;
      real_t tmp_63 = tmp_26*tmp_61;
      real_t tmp_64 = tmp_23*(tmp_19 + tmp_54);
      real_t tmp_65 = tmp_22*tmp_64;
      real_t tmp_66 = tmp_33*tmp_64;
      real_t tmp_67 = -tmp_62 - tmp_63 - tmp_65 - tmp_66 + 1;
      real_t tmp_68 = tmp_44*tmp_60;
      real_t tmp_69 = 0.2393143352496831*tmp_43;
      real_t tmp_70 = p_affine_6_1 + 0.5*tmp_9;
      real_t tmp_71 = tmp_8*(tmp_0 + tmp_70);
      real_t tmp_72 = p_affine_6_0 + 0.5*tmp_12;
      real_t tmp_73 = tmp_8*(tmp_3 + tmp_72);
      real_t tmp_74 = tmp_2*tmp_71 + tmp_5*tmp_73 - 1.0/3.0;
      real_t tmp_75 = tmp_16*tmp_73 + tmp_4*tmp_71 - 1.0/3.0;
      real_t tmp_76 = tmp_1*tmp_74 + tmp_5*tmp_75;
      real_t tmp_77 = tmp_4*tmp_74 + tmp_7*tmp_75;
      real_t tmp_78 = tmp_23*(tmp_21 + tmp_70);
      real_t tmp_79 = tmp_20*tmp_78;
      real_t tmp_80 = tmp_26*tmp_78;
      real_t tmp_81 = tmp_23*(tmp_19 + tmp_72);
      real_t tmp_82 = tmp_22*tmp_81;
      real_t tmp_83 = tmp_33*tmp_81;
      real_t tmp_84 = -tmp_79 - tmp_80 - tmp_82 - tmp_83 + 1;
      real_t tmp_85 = tmp_44*tmp_77;
      real_t tmp_86 = 0.2844444444444445*tmp_43;
      real_t tmp_87 = p_affine_6_1 + 0.7692346550528415*tmp_9;
      real_t tmp_88 = tmp_8*(tmp_0 + tmp_87);
      real_t tmp_89 = p_affine_6_0 + 0.7692346550528415*tmp_12;
      real_t tmp_90 = tmp_8*(tmp_3 + tmp_89);
      real_t tmp_91 = tmp_2*tmp_88 + tmp_5*tmp_90 - 1.0/3.0;
      real_t tmp_92 = tmp_16*tmp_90 + tmp_4*tmp_88 - 1.0/3.0;
      real_t tmp_93 = tmp_1*tmp_91 + tmp_5*tmp_92;
      real_t tmp_94 = tmp_4*tmp_91 + tmp_7*tmp_92;
      real_t tmp_95 = tmp_23*(tmp_21 + tmp_87);
      real_t tmp_96 = tmp_20*tmp_95;
      real_t tmp_97 = tmp_26*tmp_95;
      real_t tmp_98 = tmp_23*(tmp_19 + tmp_89);
      real_t tmp_99 = tmp_22*tmp_98;
      real_t tmp_100 = tmp_33*tmp_98;
      real_t tmp_101 = -tmp_100 - tmp_96 - tmp_97 - tmp_99 + 1;
      real_t tmp_102 = tmp_44*tmp_94;
      real_t tmp_103 = 0.2393143352496831*tmp_43;
      real_t tmp_104 = p_affine_6_1 + 0.95308992296933193*tmp_9;
      real_t tmp_105 = tmp_8*(tmp_0 + tmp_104);
      real_t tmp_106 = p_affine_6_0 + 0.95308992296933193*tmp_12;
      real_t tmp_107 = tmp_8*(tmp_106 + tmp_3);
      real_t tmp_108 = tmp_105*tmp_2 + tmp_107*tmp_5 - 1.0/3.0;
      real_t tmp_109 = tmp_105*tmp_4 + tmp_107*tmp_16 - 1.0/3.0;
      real_t tmp_110 = tmp_1*tmp_108 + tmp_109*tmp_5;
      real_t tmp_111 = tmp_108*tmp_4 + tmp_109*tmp_7;
      real_t tmp_112 = tmp_23*(tmp_104 + tmp_21);
      real_t tmp_113 = tmp_112*tmp_20;
      real_t tmp_114 = tmp_112*tmp_26;
      real_t tmp_115 = tmp_23*(tmp_106 + tmp_19);
      real_t tmp_116 = tmp_115*tmp_22;
      real_t tmp_117 = tmp_115*tmp_33;
      real_t tmp_118 = -tmp_113 - tmp_114 - tmp_116 - tmp_117 + 1;
      real_t tmp_119 = tmp_111*tmp_44;
      real_t tmp_120 = 0.11846344252809471*tmp_43;
      real_t tmp_121 = 1.0*p_affine_10_0*tmp_32 + 1.0*p_affine_10_1*tmp_27;
      real_t tmp_122 = tmp_38 + tmp_40;
      real_t tmp_123 = p_affine_10_0*tmp_27;
      real_t tmp_124 = tmp_63 + tmp_65;
      real_t tmp_125 = tmp_80 + tmp_82;
      real_t tmp_126 = tmp_97 + tmp_99;
      real_t tmp_127 = tmp_114 + tmp_116;
      real_t tmp_128 = 1.0*p_affine_10_0*tmp_34 + 1.0*p_affine_10_1*tmp_25;
      real_t tmp_129 = tmp_37 + tmp_41;
      real_t tmp_130 = p_affine_10_0*tmp_25;
      real_t tmp_131 = tmp_62 + tmp_66;
      real_t tmp_132 = tmp_79 + tmp_83;
      real_t tmp_133 = tmp_100 + tmp_96;
      real_t tmp_134 = tmp_113 + tmp_117;
      real_t a_0_0 = tmp_103*(-tmp_101*tmp_102 + tmp_101*tmp_50 - tmp_35*tmp_94 - tmp_59*tmp_93) + tmp_120*(-tmp_110*tmp_59 - tmp_111*tmp_35 - tmp_118*tmp_119 + tmp_118*tmp_50) + tmp_51*(-tmp_18*tmp_29 - tmp_30*tmp_35 - tmp_42*tmp_45 + tmp_42*tmp_50) + tmp_69*(-tmp_35*tmp_60 + tmp_50*tmp_67 - tmp_58*tmp_59 - tmp_67*tmp_68) + tmp_86*(-tmp_35*tmp_77 + tmp_50*tmp_84 - tmp_59*tmp_76 - tmp_84*tmp_85);
      real_t a_0_1 = tmp_103*(-tmp_102*tmp_126 - tmp_121*tmp_94 - tmp_123*tmp_93 + tmp_126*tmp_50) + tmp_120*(-tmp_110*tmp_123 - tmp_111*tmp_121 - tmp_119*tmp_127 + tmp_127*tmp_50) + tmp_51*(-tmp_121*tmp_30 - tmp_122*tmp_45 + tmp_122*tmp_50 - tmp_18*tmp_27) + tmp_69*(-tmp_121*tmp_60 - tmp_123*tmp_58 + tmp_124*tmp_50 - tmp_124*tmp_68) + tmp_86*(-tmp_121*tmp_77 - tmp_123*tmp_76 + tmp_125*tmp_50 - tmp_125*tmp_85);
      real_t a_0_2 = tmp_103*(-tmp_102*tmp_133 - tmp_128*tmp_94 - tmp_130*tmp_93 + tmp_133*tmp_50) + tmp_120*(-tmp_110*tmp_130 - tmp_111*tmp_128 - tmp_119*tmp_134 + tmp_134*tmp_50) + tmp_51*(-tmp_128*tmp_30 - tmp_129*tmp_45 + tmp_129*tmp_50 - tmp_18*tmp_25) + tmp_69*(-tmp_128*tmp_60 - tmp_130*tmp_58 + tmp_131*tmp_50 - tmp_131*tmp_68) + tmp_86*(-tmp_128*tmp_77 - tmp_130*tmp_76 + tmp_132*tmp_50 - tmp_132*tmp_85);
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
      real_t tmp_1 = p_affine_1_1 + tmp_0;
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = p_affine_1_0 + tmp_3;
      real_t tmp_5 = p_affine_2_1 + tmp_0;
      real_t tmp_6 = p_affine_2_0 + tmp_3;
      real_t tmp_7 = 1.0 / (-tmp_1*tmp_6 + tmp_4*tmp_5);
      real_t tmp_8 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9 = p_affine_6_1 + tmp_0;
      real_t tmp_10 = tmp_7*(0.046910077030668018*tmp_8 + tmp_9);
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_3;
      real_t tmp_13 = tmp_7*(0.046910077030668018*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_10*tmp_2 + tmp_13*tmp_5 - 1.0/3.0;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_10*tmp_4 + tmp_13*tmp_15 - 1.0/3.0;
      real_t tmp_17 = p_affine_10_0*(tmp_1*tmp_14 + tmp_16*tmp_5);
      real_t tmp_18 = 0.5*tmp_7;
      real_t tmp_19 = tmp_18*tmp_4;
      real_t tmp_20 = tmp_18*tmp_2;
      real_t tmp_21 = -tmp_19 - tmp_20;
      real_t tmp_22 = 2*tmp_21;
      real_t tmp_23 = tmp_14*tmp_4 + tmp_16*tmp_6;
      real_t tmp_24 = 1.0*tmp_7;
      real_t tmp_25 = tmp_24*tmp_5;
      real_t tmp_26 = tmp_15*tmp_24;
      real_t tmp_27 = 2*p_affine_10_0*(-tmp_25 - tmp_26) + 2*p_affine_10_1*tmp_21;
      real_t tmp_28 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_8*tmp_8), 1.0/2.0));
      real_t tmp_29 = 0.11846344252809471*tmp_28;
      real_t tmp_30 = tmp_7*(0.23076534494715845*tmp_8 + tmp_9);
      real_t tmp_31 = tmp_7*(0.23076534494715845*tmp_11 + tmp_12);
      real_t tmp_32 = tmp_2*tmp_30 + tmp_31*tmp_5 - 1.0/3.0;
      real_t tmp_33 = tmp_15*tmp_31 + tmp_30*tmp_4 - 1.0/3.0;
      real_t tmp_34 = tmp_1*tmp_32 + tmp_33*tmp_5;
      real_t tmp_35 = p_affine_10_0*tmp_22;
      real_t tmp_36 = tmp_32*tmp_4 + tmp_33*tmp_6;
      real_t tmp_37 = 0.2393143352496831*tmp_28;
      real_t tmp_38 = tmp_7*(0.5*tmp_8 + tmp_9);
      real_t tmp_39 = tmp_7*(0.5*tmp_11 + tmp_12);
      real_t tmp_40 = tmp_2*tmp_38 + tmp_39*tmp_5 - 1.0/3.0;
      real_t tmp_41 = tmp_15*tmp_39 + tmp_38*tmp_4 - 1.0/3.0;
      real_t tmp_42 = tmp_1*tmp_40 + tmp_41*tmp_5;
      real_t tmp_43 = tmp_4*tmp_40 + tmp_41*tmp_6;
      real_t tmp_44 = 0.2844444444444445*tmp_28;
      real_t tmp_45 = tmp_7*(0.7692346550528415*tmp_8 + tmp_9);
      real_t tmp_46 = tmp_7*(0.7692346550528415*tmp_11 + tmp_12);
      real_t tmp_47 = tmp_2*tmp_45 + tmp_46*tmp_5 - 1.0/3.0;
      real_t tmp_48 = tmp_15*tmp_46 + tmp_4*tmp_45 - 1.0/3.0;
      real_t tmp_49 = tmp_1*tmp_47 + tmp_48*tmp_5;
      real_t tmp_50 = tmp_4*tmp_47 + tmp_48*tmp_6;
      real_t tmp_51 = 0.2393143352496831*tmp_28;
      real_t tmp_52 = tmp_7*(0.95308992296933193*tmp_8 + tmp_9);
      real_t tmp_53 = tmp_7*(0.95308992296933193*tmp_11 + tmp_12);
      real_t tmp_54 = tmp_2*tmp_52 + tmp_5*tmp_53 - 1.0/3.0;
      real_t tmp_55 = tmp_15*tmp_53 + tmp_4*tmp_52 - 1.0/3.0;
      real_t tmp_56 = tmp_1*tmp_54 + tmp_5*tmp_55;
      real_t tmp_57 = tmp_4*tmp_54 + tmp_55*tmp_6;
      real_t tmp_58 = 0.11846344252809471*tmp_28;
      real_t tmp_59 = tmp_2*tmp_24;
      real_t tmp_60 = 2*p_affine_10_0*tmp_25 + 2*p_affine_10_1*tmp_20;
      real_t tmp_61 = p_affine_10_0*tmp_59;
      real_t tmp_62 = tmp_24*tmp_4;
      real_t tmp_63 = 2*p_affine_10_0*tmp_26 + 2*p_affine_10_1*tmp_19;
      real_t tmp_64 = p_affine_10_0*tmp_62;
      real_t a_0_0 = tmp_29*(-tmp_17*tmp_22 - tmp_23*tmp_27) + tmp_37*(-tmp_27*tmp_36 - tmp_34*tmp_35) + tmp_44*(-tmp_27*tmp_43 - tmp_35*tmp_42) + tmp_51*(-tmp_27*tmp_50 - tmp_35*tmp_49) + tmp_58*(-tmp_27*tmp_57 - tmp_35*tmp_56);
      real_t a_0_1 = tmp_29*(-tmp_17*tmp_59 - tmp_23*tmp_60) + tmp_37*(-tmp_34*tmp_61 - tmp_36*tmp_60) + tmp_44*(-tmp_42*tmp_61 - tmp_43*tmp_60) + tmp_51*(-tmp_49*tmp_61 - tmp_50*tmp_60) + tmp_58*(-tmp_56*tmp_61 - tmp_57*tmp_60);
      real_t a_0_2 = tmp_29*(-tmp_17*tmp_62 - tmp_23*tmp_63) + tmp_37*(-tmp_34*tmp_64 - tmp_36*tmp_63) + tmp_44*(-tmp_42*tmp_64 - tmp_43*tmp_63) + tmp_51*(-tmp_49*tmp_64 - tmp_50*tmp_63) + tmp_58*(-tmp_56*tmp_64 - tmp_57*tmp_63);
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




class EGConstEpsilonFormP1E_0 : public hyteg::dg::DGForm2D
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

      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_0;
      real_t tmp_3 = p_affine_1_0 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_2;
      real_t tmp_6 = p_affine_1_1 + tmp_0;
      real_t tmp_7 = 1.0 / (tmp_4 - tmp_5*tmp_6);
      real_t tmp_8 = 1.0*tmp_7;
      real_t tmp_9 = tmp_1*tmp_8;
      real_t tmp_10 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_11 = tmp_10*tmp_8;
      real_t tmp_12 = 2.0*tmp_7;
      real_t tmp_13 = tmp_10*tmp_12*tmp_5 + tmp_12*tmp_4;
      real_t tmp_14 = 0.5*tmp_7;
      real_t tmp_15 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_16 = tmp_3*tmp_7;
      real_t tmp_17 = tmp_1*tmp_7;
      real_t tmp_18 = tmp_10*tmp_17 + tmp_15*tmp_16 + tmp_16*tmp_5 + tmp_17*tmp_6;
      real_t tmp_19 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_20 = tmp_19*(tmp_13*(-tmp_11 - tmp_9) + 2*tmp_18*(-tmp_14*tmp_15 - tmp_14*tmp_3));
      real_t tmp_21 = tmp_18*tmp_8;
      real_t tmp_22 = tmp_19*(tmp_13*tmp_9 + tmp_15*tmp_21);
      real_t tmp_23 = tmp_19*(tmp_11*tmp_13 + tmp_21*tmp_3);
      real_t a_0_0 = 0.5*tmp_20;
      real_t a_1_0 = 0.5*tmp_22;
      real_t a_2_0 = 0.5*tmp_23;
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
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = p_affine_1_0 + tmp_3;
      real_t tmp_5 = p_affine_2_1 + tmp_0;
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = p_affine_2_0 + tmp_3;
      real_t tmp_8 = 1.0 / (-tmp_1*tmp_7 + tmp_6);
      real_t tmp_9 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + tmp_0;
      real_t tmp_11 = tmp_8*(tmp_10 + 0.046910077030668018*tmp_9);
      real_t tmp_12 = tmp_11*tmp_2;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_3;
      real_t tmp_15 = tmp_8*(0.046910077030668018*tmp_13 + tmp_14);
      real_t tmp_16 = tmp_15*tmp_5;
      real_t tmp_17 = tmp_12 + tmp_16;
      real_t tmp_18 = tmp_17 - 1.0/3.0;
      real_t tmp_19 = tmp_11*tmp_4;
      real_t tmp_20 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_21 = tmp_15*tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_22 - 1.0/3.0;
      real_t tmp_24 = p_affine_10_0*(tmp_1*tmp_18 + tmp_23*tmp_5);
      real_t tmp_25 = 0.5*tmp_8;
      real_t tmp_26 = tmp_25*tmp_4;
      real_t tmp_27 = tmp_2*tmp_25;
      real_t tmp_28 = -tmp_26 - tmp_27;
      real_t tmp_29 = 1.0*tmp_28;
      real_t tmp_30 = tmp_18*tmp_4 + tmp_23*tmp_7;
      real_t tmp_31 = 1.0*tmp_8;
      real_t tmp_32 = tmp_31*tmp_5;
      real_t tmp_33 = tmp_20*tmp_31;
      real_t tmp_34 = 1.0*p_affine_10_0*(-tmp_32 - tmp_33) + 1.0*p_affine_10_1*tmp_28;
      real_t tmp_35 = -tmp_12 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_36 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_37 = 12/tmp_36;
      real_t tmp_38 = tmp_30*tmp_37;
      real_t tmp_39 = tmp_25*tmp_5;
      real_t tmp_40 = 1.0*p_affine_10_0*(tmp_31*tmp_6 + tmp_33*tmp_7) + 1.0*p_affine_10_1*(tmp_1*tmp_39 + tmp_20*tmp_39 + tmp_26*tmp_7 + tmp_27*tmp_4);
      real_t tmp_41 = 0.11846344252809471*tmp_36;
      real_t tmp_42 = tmp_8*(tmp_10 + 0.23076534494715845*tmp_9);
      real_t tmp_43 = tmp_2*tmp_42;
      real_t tmp_44 = tmp_8*(0.23076534494715845*tmp_13 + tmp_14);
      real_t tmp_45 = tmp_44*tmp_5;
      real_t tmp_46 = tmp_43 + tmp_45;
      real_t tmp_47 = tmp_46 - 1.0/3.0;
      real_t tmp_48 = tmp_4*tmp_42;
      real_t tmp_49 = tmp_20*tmp_44;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_50 - 1.0/3.0;
      real_t tmp_52 = tmp_1*tmp_47 + tmp_5*tmp_51;
      real_t tmp_53 = p_affine_10_0*tmp_29;
      real_t tmp_54 = tmp_4*tmp_47 + tmp_51*tmp_7;
      real_t tmp_55 = -tmp_43 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_56 = tmp_37*tmp_54;
      real_t tmp_57 = 0.2393143352496831*tmp_36;
      real_t tmp_58 = tmp_8*(tmp_10 + 0.5*tmp_9);
      real_t tmp_59 = tmp_2*tmp_58;
      real_t tmp_60 = tmp_8*(0.5*tmp_13 + tmp_14);
      real_t tmp_61 = tmp_5*tmp_60;
      real_t tmp_62 = tmp_59 + tmp_61;
      real_t tmp_63 = tmp_62 - 1.0/3.0;
      real_t tmp_64 = tmp_4*tmp_58;
      real_t tmp_65 = tmp_20*tmp_60;
      real_t tmp_66 = tmp_64 + tmp_65;
      real_t tmp_67 = tmp_66 - 1.0/3.0;
      real_t tmp_68 = tmp_1*tmp_63 + tmp_5*tmp_67;
      real_t tmp_69 = tmp_4*tmp_63 + tmp_67*tmp_7;
      real_t tmp_70 = -tmp_59 - tmp_61 - tmp_64 - tmp_65 + 1;
      real_t tmp_71 = tmp_37*tmp_69;
      real_t tmp_72 = 0.2844444444444445*tmp_36;
      real_t tmp_73 = tmp_8*(tmp_10 + 0.7692346550528415*tmp_9);
      real_t tmp_74 = tmp_2*tmp_73;
      real_t tmp_75 = tmp_8*(0.7692346550528415*tmp_13 + tmp_14);
      real_t tmp_76 = tmp_5*tmp_75;
      real_t tmp_77 = tmp_74 + tmp_76;
      real_t tmp_78 = tmp_77 - 1.0/3.0;
      real_t tmp_79 = tmp_4*tmp_73;
      real_t tmp_80 = tmp_20*tmp_75;
      real_t tmp_81 = tmp_79 + tmp_80;
      real_t tmp_82 = tmp_81 - 1.0/3.0;
      real_t tmp_83 = tmp_1*tmp_78 + tmp_5*tmp_82;
      real_t tmp_84 = tmp_4*tmp_78 + tmp_7*tmp_82;
      real_t tmp_85 = -tmp_74 - tmp_76 - tmp_79 - tmp_80 + 1;
      real_t tmp_86 = tmp_37*tmp_84;
      real_t tmp_87 = 0.2393143352496831*tmp_36;
      real_t tmp_88 = tmp_8*(tmp_10 + 0.95308992296933193*tmp_9);
      real_t tmp_89 = tmp_2*tmp_88;
      real_t tmp_90 = tmp_8*(0.95308992296933193*tmp_13 + tmp_14);
      real_t tmp_91 = tmp_5*tmp_90;
      real_t tmp_92 = tmp_89 + tmp_91;
      real_t tmp_93 = tmp_92 - 1.0/3.0;
      real_t tmp_94 = tmp_4*tmp_88;
      real_t tmp_95 = tmp_20*tmp_90;
      real_t tmp_96 = tmp_94 + tmp_95;
      real_t tmp_97 = tmp_96 - 1.0/3.0;
      real_t tmp_98 = tmp_1*tmp_93 + tmp_5*tmp_97;
      real_t tmp_99 = tmp_4*tmp_93 + tmp_7*tmp_97;
      real_t tmp_100 = -tmp_89 - tmp_91 - tmp_94 - tmp_95 + 1;
      real_t tmp_101 = tmp_37*tmp_99;
      real_t tmp_102 = 0.11846344252809471*tmp_36;
      real_t tmp_103 = 1.0*p_affine_10_0*tmp_32 + 1.0*p_affine_10_1*tmp_27;
      real_t tmp_104 = p_affine_10_0*tmp_27;
      real_t tmp_105 = 1.0*p_affine_10_0*tmp_33 + 1.0*p_affine_10_1*tmp_26;
      real_t tmp_106 = p_affine_10_0*tmp_26;
      real_t a_0_0 = tmp_102*(tmp_100*tmp_101 - tmp_100*tmp_40 - tmp_34*tmp_99 - tmp_53*tmp_98) + tmp_41*(-tmp_24*tmp_29 - tmp_30*tmp_34 + tmp_35*tmp_38 - tmp_35*tmp_40) + tmp_57*(-tmp_34*tmp_54 - tmp_40*tmp_55 - tmp_52*tmp_53 + tmp_55*tmp_56) + tmp_72*(-tmp_34*tmp_69 - tmp_40*tmp_70 - tmp_53*tmp_68 + tmp_70*tmp_71) + tmp_87*(-tmp_34*tmp_84 - tmp_40*tmp_85 - tmp_53*tmp_83 + tmp_85*tmp_86);
      real_t a_1_0 = tmp_102*(tmp_101*tmp_92 - tmp_103*tmp_99 - tmp_104*tmp_98 - tmp_40*tmp_92) + tmp_41*(-tmp_103*tmp_30 + tmp_17*tmp_38 - tmp_17*tmp_40 - tmp_24*tmp_27) + tmp_57*(-tmp_103*tmp_54 - tmp_104*tmp_52 - tmp_40*tmp_46 + tmp_46*tmp_56) + tmp_72*(-tmp_103*tmp_69 - tmp_104*tmp_68 - tmp_40*tmp_62 + tmp_62*tmp_71) + tmp_87*(-tmp_103*tmp_84 - tmp_104*tmp_83 - tmp_40*tmp_77 + tmp_77*tmp_86);
      real_t a_2_0 = tmp_102*(tmp_101*tmp_96 - tmp_105*tmp_99 - tmp_106*tmp_98 - tmp_40*tmp_96) + tmp_41*(-tmp_105*tmp_30 + tmp_22*tmp_38 - tmp_22*tmp_40 - tmp_24*tmp_26) + tmp_57*(-tmp_105*tmp_54 - tmp_106*tmp_52 - tmp_40*tmp_50 + tmp_50*tmp_56) + tmp_72*(-tmp_105*tmp_69 - tmp_106*tmp_68 - tmp_40*tmp_66 + tmp_66*tmp_71) + tmp_87*(-tmp_105*tmp_84 - tmp_106*tmp_83 - tmp_40*tmp_81 + tmp_81*tmp_86);
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
      real_t tmp_7 = p_affine_5_0 + tmp_3;
      real_t tmp_8 = 1.0 / (-tmp_1*tmp_7 + tmp_6);
      real_t tmp_9 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + 0.046910077030668018*tmp_9;
      real_t tmp_11 = tmp_8*(tmp_0 + tmp_10);
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + 0.046910077030668018*tmp_12;
      real_t tmp_14 = tmp_8*(tmp_13 + tmp_3);
      real_t tmp_15 = tmp_11*tmp_2 + tmp_14*tmp_5 - 1.0/3.0;
      real_t tmp_16 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_17 = tmp_11*tmp_4 + tmp_14*tmp_16 - 1.0/3.0;
      real_t tmp_18 = tmp_1*tmp_15 + tmp_17*tmp_5;
      real_t tmp_19 = -p_affine_0_0;
      real_t tmp_20 = p_affine_1_0 + tmp_19;
      real_t tmp_21 = -p_affine_0_1;
      real_t tmp_22 = p_affine_2_1 + tmp_21;
      real_t tmp_23 = 1.0 / (tmp_20*tmp_22 - (p_affine_1_1 + tmp_21)*(p_affine_2_0 + tmp_19));
      real_t tmp_24 = 0.5*tmp_23;
      real_t tmp_25 = tmp_20*tmp_24;
      real_t tmp_26 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_25 - tmp_27;
      real_t tmp_29 = p_affine_10_0*tmp_28;
      real_t tmp_30 = 1.0*tmp_23;
      real_t tmp_31 = tmp_22*tmp_30;
      real_t tmp_32 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_33 = tmp_30*tmp_32;
      real_t tmp_34 = p_affine_10_0*(-tmp_31 - tmp_33) + p_affine_10_1*tmp_28;
      real_t tmp_35 = tmp_15*tmp_4 + tmp_17*tmp_7;
      real_t tmp_36 = tmp_23*(tmp_10 + tmp_21);
      real_t tmp_37 = tmp_20*tmp_36;
      real_t tmp_38 = tmp_26*tmp_36;
      real_t tmp_39 = tmp_23*(tmp_13 + tmp_19);
      real_t tmp_40 = tmp_22*tmp_39;
      real_t tmp_41 = tmp_32*tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_44 = 12/tmp_43;
      real_t tmp_45 = tmp_35*tmp_44;
      real_t tmp_46 = 1.0*tmp_8;
      real_t tmp_47 = 0.5*tmp_8;
      real_t tmp_48 = tmp_4*tmp_47;
      real_t tmp_49 = tmp_47*tmp_5;
      real_t tmp_50 = 1.0*p_affine_10_0*(tmp_16*tmp_46*tmp_7 + tmp_46*tmp_6) + 1.0*p_affine_10_1*(tmp_1*tmp_49 + tmp_16*tmp_49 + tmp_2*tmp_48 + tmp_48*tmp_7);
      real_t tmp_51 = 0.11846344252809471*tmp_43;
      real_t tmp_52 = p_affine_6_1 + 0.23076534494715845*tmp_9;
      real_t tmp_53 = tmp_8*(tmp_0 + tmp_52);
      real_t tmp_54 = p_affine_6_0 + 0.23076534494715845*tmp_12;
      real_t tmp_55 = tmp_8*(tmp_3 + tmp_54);
      real_t tmp_56 = tmp_2*tmp_53 + tmp_5*tmp_55 - 1.0/3.0;
      real_t tmp_57 = tmp_16*tmp_55 + tmp_4*tmp_53 - 1.0/3.0;
      real_t tmp_58 = tmp_1*tmp_56 + tmp_5*tmp_57;
      real_t tmp_59 = tmp_4*tmp_56 + tmp_57*tmp_7;
      real_t tmp_60 = tmp_23*(tmp_21 + tmp_52);
      real_t tmp_61 = tmp_20*tmp_60;
      real_t tmp_62 = tmp_26*tmp_60;
      real_t tmp_63 = tmp_23*(tmp_19 + tmp_54);
      real_t tmp_64 = tmp_22*tmp_63;
      real_t tmp_65 = tmp_32*tmp_63;
      real_t tmp_66 = -tmp_61 - tmp_62 - tmp_64 - tmp_65 + 1;
      real_t tmp_67 = tmp_44*tmp_59;
      real_t tmp_68 = 0.2393143352496831*tmp_43;
      real_t tmp_69 = p_affine_6_1 + 0.5*tmp_9;
      real_t tmp_70 = tmp_8*(tmp_0 + tmp_69);
      real_t tmp_71 = p_affine_6_0 + 0.5*tmp_12;
      real_t tmp_72 = tmp_8*(tmp_3 + tmp_71);
      real_t tmp_73 = tmp_2*tmp_70 + tmp_5*tmp_72 - 1.0/3.0;
      real_t tmp_74 = tmp_16*tmp_72 + tmp_4*tmp_70 - 1.0/3.0;
      real_t tmp_75 = tmp_1*tmp_73 + tmp_5*tmp_74;
      real_t tmp_76 = tmp_4*tmp_73 + tmp_7*tmp_74;
      real_t tmp_77 = tmp_23*(tmp_21 + tmp_69);
      real_t tmp_78 = tmp_20*tmp_77;
      real_t tmp_79 = tmp_26*tmp_77;
      real_t tmp_80 = tmp_23*(tmp_19 + tmp_71);
      real_t tmp_81 = tmp_22*tmp_80;
      real_t tmp_82 = tmp_32*tmp_80;
      real_t tmp_83 = -tmp_78 - tmp_79 - tmp_81 - tmp_82 + 1;
      real_t tmp_84 = tmp_44*tmp_76;
      real_t tmp_85 = 0.2844444444444445*tmp_43;
      real_t tmp_86 = p_affine_6_1 + 0.7692346550528415*tmp_9;
      real_t tmp_87 = tmp_8*(tmp_0 + tmp_86);
      real_t tmp_88 = p_affine_6_0 + 0.7692346550528415*tmp_12;
      real_t tmp_89 = tmp_8*(tmp_3 + tmp_88);
      real_t tmp_90 = tmp_2*tmp_87 + tmp_5*tmp_89 - 1.0/3.0;
      real_t tmp_91 = tmp_16*tmp_89 + tmp_4*tmp_87 - 1.0/3.0;
      real_t tmp_92 = tmp_1*tmp_90 + tmp_5*tmp_91;
      real_t tmp_93 = tmp_4*tmp_90 + tmp_7*tmp_91;
      real_t tmp_94 = tmp_23*(tmp_21 + tmp_86);
      real_t tmp_95 = tmp_20*tmp_94;
      real_t tmp_96 = tmp_26*tmp_94;
      real_t tmp_97 = tmp_23*(tmp_19 + tmp_88);
      real_t tmp_98 = tmp_22*tmp_97;
      real_t tmp_99 = tmp_32*tmp_97;
      real_t tmp_100 = -tmp_95 - tmp_96 - tmp_98 - tmp_99 + 1;
      real_t tmp_101 = tmp_44*tmp_93;
      real_t tmp_102 = 0.2393143352496831*tmp_43;
      real_t tmp_103 = p_affine_6_1 + 0.95308992296933193*tmp_9;
      real_t tmp_104 = tmp_8*(tmp_0 + tmp_103);
      real_t tmp_105 = p_affine_6_0 + 0.95308992296933193*tmp_12;
      real_t tmp_106 = tmp_8*(tmp_105 + tmp_3);
      real_t tmp_107 = tmp_104*tmp_2 + tmp_106*tmp_5 - 1.0/3.0;
      real_t tmp_108 = tmp_104*tmp_4 + tmp_106*tmp_16 - 1.0/3.0;
      real_t tmp_109 = tmp_1*tmp_107 + tmp_108*tmp_5;
      real_t tmp_110 = tmp_107*tmp_4 + tmp_108*tmp_7;
      real_t tmp_111 = tmp_23*(tmp_103 + tmp_21);
      real_t tmp_112 = tmp_111*tmp_20;
      real_t tmp_113 = tmp_111*tmp_26;
      real_t tmp_114 = tmp_23*(tmp_105 + tmp_19);
      real_t tmp_115 = tmp_114*tmp_22;
      real_t tmp_116 = tmp_114*tmp_32;
      real_t tmp_117 = -tmp_112 - tmp_113 - tmp_115 - tmp_116 + 1;
      real_t tmp_118 = tmp_110*tmp_44;
      real_t tmp_119 = 0.11846344252809471*tmp_43;
      real_t tmp_120 = p_affine_10_0*tmp_27;
      real_t tmp_121 = p_affine_10_0*tmp_31 + p_affine_10_1*tmp_27;
      real_t tmp_122 = tmp_38 + tmp_40;
      real_t tmp_123 = tmp_62 + tmp_64;
      real_t tmp_124 = tmp_79 + tmp_81;
      real_t tmp_125 = tmp_96 + tmp_98;
      real_t tmp_126 = tmp_113 + tmp_115;
      real_t tmp_127 = p_affine_10_0*tmp_25;
      real_t tmp_128 = p_affine_10_0*tmp_33 + p_affine_10_1*tmp_25;
      real_t tmp_129 = tmp_37 + tmp_41;
      real_t tmp_130 = tmp_61 + tmp_65;
      real_t tmp_131 = tmp_78 + tmp_82;
      real_t tmp_132 = tmp_95 + tmp_99;
      real_t tmp_133 = tmp_112 + tmp_116;
      real_t a_0_0 = tmp_102*(-tmp_100*tmp_101 - tmp_100*tmp_50 + tmp_29*tmp_92 + tmp_34*tmp_93) + tmp_119*(tmp_109*tmp_29 + tmp_110*tmp_34 - tmp_117*tmp_118 - tmp_117*tmp_50) + tmp_51*(tmp_18*tmp_29 + tmp_34*tmp_35 - tmp_42*tmp_45 - tmp_42*tmp_50) + tmp_68*(tmp_29*tmp_58 + tmp_34*tmp_59 - tmp_50*tmp_66 - tmp_66*tmp_67) + tmp_85*(tmp_29*tmp_75 + tmp_34*tmp_76 - tmp_50*tmp_83 - tmp_83*tmp_84);
      real_t a_1_0 = tmp_102*(-tmp_101*tmp_125 + tmp_120*tmp_92 + tmp_121*tmp_93 - tmp_125*tmp_50) + tmp_119*(tmp_109*tmp_120 + tmp_110*tmp_121 - tmp_118*tmp_126 - tmp_126*tmp_50) + tmp_51*(tmp_120*tmp_18 + tmp_121*tmp_35 - tmp_122*tmp_45 - tmp_122*tmp_50) + tmp_68*(tmp_120*tmp_58 + tmp_121*tmp_59 - tmp_123*tmp_50 - tmp_123*tmp_67) + tmp_85*(tmp_120*tmp_75 + tmp_121*tmp_76 - tmp_124*tmp_50 - tmp_124*tmp_84);
      real_t a_2_0 = tmp_102*(-tmp_101*tmp_132 + tmp_127*tmp_92 + tmp_128*tmp_93 - tmp_132*tmp_50) + tmp_119*(tmp_109*tmp_127 + tmp_110*tmp_128 - tmp_118*tmp_133 - tmp_133*tmp_50) + tmp_51*(tmp_127*tmp_18 + tmp_128*tmp_35 - tmp_129*tmp_45 - tmp_129*tmp_50) + tmp_68*(tmp_127*tmp_58 + tmp_128*tmp_59 - tmp_130*tmp_50 - tmp_130*tmp_67) + tmp_85*(tmp_127*tmp_75 + tmp_128*tmp_76 - tmp_131*tmp_50 - tmp_131*tmp_84);
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
      real_t tmp_5 = 1.0 / (tmp_1*tmp_3 - tmp_4*(p_affine_1_1 + tmp_2));
      real_t tmp_6 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_7 = p_affine_6_1 + tmp_2;
      real_t tmp_8 = tmp_5*(0.046910077030668018*tmp_6 + tmp_7);
      real_t tmp_9 = tmp_1*tmp_8;
      real_t tmp_10 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_11 = tmp_10*tmp_8;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_0;
      real_t tmp_14 = tmp_5*(0.046910077030668018*tmp_12 + tmp_13);
      real_t tmp_15 = tmp_14*tmp_3;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_14*tmp_16;
      real_t tmp_18 = tmp_11 + tmp_15;
      real_t tmp_19 = tmp_17 + tmp_9;
      real_t tmp_20 = 1.4215613103371365*tmp_1*(tmp_18 - 1.0/3.0) + 1.4215613103371365*tmp_4*(tmp_19 - 1.0/3.0);
      real_t tmp_21 = tmp_5*(0.23076534494715845*tmp_6 + tmp_7);
      real_t tmp_22 = tmp_1*tmp_21;
      real_t tmp_23 = tmp_10*tmp_21;
      real_t tmp_24 = tmp_5*(0.23076534494715845*tmp_12 + tmp_13);
      real_t tmp_25 = tmp_24*tmp_3;
      real_t tmp_26 = tmp_16*tmp_24;
      real_t tmp_27 = tmp_23 + tmp_25;
      real_t tmp_28 = tmp_22 + tmp_26;
      real_t tmp_29 = 2.8717720229961969*tmp_1*(tmp_27 - 1.0/3.0) + 2.8717720229961969*tmp_4*(tmp_28 - 1.0/3.0);
      real_t tmp_30 = tmp_5*(0.5*tmp_6 + tmp_7);
      real_t tmp_31 = tmp_1*tmp_30;
      real_t tmp_32 = tmp_10*tmp_30;
      real_t tmp_33 = tmp_5*(0.5*tmp_12 + tmp_13);
      real_t tmp_34 = tmp_3*tmp_33;
      real_t tmp_35 = tmp_16*tmp_33;
      real_t tmp_36 = tmp_32 + tmp_34;
      real_t tmp_37 = tmp_31 + tmp_35;
      real_t tmp_38 = 3.413333333333334*tmp_1*(tmp_36 - 1.0/3.0) + 3.413333333333334*tmp_4*(tmp_37 - 1.0/3.0);
      real_t tmp_39 = tmp_5*(0.7692346550528415*tmp_6 + tmp_7);
      real_t tmp_40 = tmp_1*tmp_39;
      real_t tmp_41 = tmp_10*tmp_39;
      real_t tmp_42 = tmp_5*(0.7692346550528415*tmp_12 + tmp_13);
      real_t tmp_43 = tmp_3*tmp_42;
      real_t tmp_44 = tmp_16*tmp_42;
      real_t tmp_45 = tmp_41 + tmp_43;
      real_t tmp_46 = tmp_40 + tmp_44;
      real_t tmp_47 = 2.8717720229961969*tmp_1*(tmp_45 - 1.0/3.0) + 2.8717720229961969*tmp_4*(tmp_46 - 1.0/3.0);
      real_t tmp_48 = tmp_5*(0.95308992296933193*tmp_6 + tmp_7);
      real_t tmp_49 = tmp_1*tmp_48;
      real_t tmp_50 = tmp_10*tmp_48;
      real_t tmp_51 = tmp_5*(0.95308992296933193*tmp_12 + tmp_13);
      real_t tmp_52 = tmp_3*tmp_51;
      real_t tmp_53 = tmp_16*tmp_51;
      real_t tmp_54 = tmp_50 + tmp_52;
      real_t tmp_55 = tmp_49 + tmp_53;
      real_t tmp_56 = 1.4215613103371365*tmp_1*(tmp_54 - 1.0/3.0) + 1.4215613103371365*tmp_4*(tmp_55 - 1.0/3.0);
      real_t a_0_0 = tmp_20*(-tmp_11 - tmp_15 - tmp_17 - tmp_9 + 1) + tmp_29*(-tmp_22 - tmp_23 - tmp_25 - tmp_26 + 1) + tmp_38*(-tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1) + tmp_47*(-tmp_40 - tmp_41 - tmp_43 - tmp_44 + 1) + tmp_56*(-tmp_49 - tmp_50 - tmp_52 - tmp_53 + 1);
      real_t a_1_0 = tmp_18*tmp_20 + tmp_27*tmp_29 + tmp_36*tmp_38 + tmp_45*tmp_47 + tmp_54*tmp_56;
      real_t a_2_0 = tmp_19*tmp_20 + tmp_28*tmp_29 + tmp_37*tmp_38 + tmp_46*tmp_47 + tmp_55*tmp_56;
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




class EGConstEpsilonFormEP1_1 : public hyteg::dg::DGForm2D
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
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = 1.0 / (tmp_4 - tmp_5*tmp_6);
      real_t tmp_8 = 2.0*tmp_7;
      real_t tmp_9 = tmp_1*tmp_8;
      real_t tmp_10 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_11 = tmp_10*tmp_8;
      real_t tmp_12 = 1.0*tmp_7;
      real_t tmp_13 = tmp_10*tmp_12*tmp_6 + tmp_12*tmp_4;
      real_t tmp_14 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_15 = 0.5*tmp_7;
      real_t tmp_16 = tmp_1*tmp_15;
      real_t tmp_17 = tmp_15*tmp_3;
      real_t tmp_18 = tmp_10*tmp_16 + tmp_14*tmp_17 + tmp_16*tmp_5 + tmp_17*tmp_6;
      real_t tmp_19 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_20 = tmp_19*(tmp_13*(-tmp_11 - tmp_9) + 2*tmp_18*(-tmp_12*tmp_14 - tmp_12*tmp_3));
      real_t tmp_21 = tmp_18*tmp_8;
      real_t tmp_22 = tmp_19*(tmp_11*tmp_13 + tmp_21*tmp_3);
      real_t tmp_23 = tmp_19*(tmp_13*tmp_9 + tmp_14*tmp_21);
      real_t a_0_0 = 0.5*tmp_20;
      real_t a_0_1 = 0.5*tmp_22;
      real_t a_0_2 = 0.5*tmp_23;
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
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = tmp_1*tmp_4;
      real_t tmp_6 = p_affine_2_0 + tmp_0;
      real_t tmp_7 = p_affine_1_1 + tmp_3;
      real_t tmp_8 = 1.0 / (tmp_5 - tmp_6*tmp_7);
      real_t tmp_9 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + tmp_3;
      real_t tmp_11 = tmp_8*(tmp_10 + 0.046910077030668018*tmp_9);
      real_t tmp_12 = tmp_11*tmp_2;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_0;
      real_t tmp_15 = tmp_8*(0.046910077030668018*tmp_13 + tmp_14);
      real_t tmp_16 = tmp_15*tmp_4;
      real_t tmp_17 = tmp_12 + tmp_16;
      real_t tmp_18 = tmp_17 - 1.0/3.0;
      real_t tmp_19 = tmp_1*tmp_11;
      real_t tmp_20 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_21 = tmp_15*tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_22 - 1.0/3.0;
      real_t tmp_24 = p_affine_10_1*(tmp_1*tmp_18 + tmp_23*tmp_6);
      real_t tmp_25 = 0.5*tmp_8;
      real_t tmp_26 = tmp_25*tmp_4;
      real_t tmp_27 = tmp_20*tmp_25;
      real_t tmp_28 = -tmp_26 - tmp_27;
      real_t tmp_29 = 1.0*tmp_28;
      real_t tmp_30 = tmp_18*tmp_7 + tmp_23*tmp_4;
      real_t tmp_31 = 1.0*tmp_8;
      real_t tmp_32 = tmp_1*tmp_31;
      real_t tmp_33 = tmp_2*tmp_31;
      real_t tmp_34 = 1.0*p_affine_10_0*tmp_28 + 1.0*p_affine_10_1*(-tmp_32 - tmp_33);
      real_t tmp_35 = -tmp_12 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_36 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_37 = 12/tmp_36;
      real_t tmp_38 = tmp_30*tmp_37;
      real_t tmp_39 = tmp_1*tmp_25;
      real_t tmp_40 = 1.0*p_affine_10_0*(tmp_2*tmp_39 + tmp_26*tmp_7 + tmp_27*tmp_4 + tmp_39*tmp_6) + 1.0*p_affine_10_1*(tmp_31*tmp_5 + tmp_33*tmp_7);
      real_t tmp_41 = 0.11846344252809471*tmp_36;
      real_t tmp_42 = tmp_8*(tmp_10 + 0.23076534494715845*tmp_9);
      real_t tmp_43 = tmp_2*tmp_42;
      real_t tmp_44 = tmp_8*(0.23076534494715845*tmp_13 + tmp_14);
      real_t tmp_45 = tmp_4*tmp_44;
      real_t tmp_46 = tmp_43 + tmp_45;
      real_t tmp_47 = tmp_46 - 1.0/3.0;
      real_t tmp_48 = tmp_1*tmp_42;
      real_t tmp_49 = tmp_20*tmp_44;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_50 - 1.0/3.0;
      real_t tmp_52 = tmp_1*tmp_47 + tmp_51*tmp_6;
      real_t tmp_53 = p_affine_10_1*tmp_29;
      real_t tmp_54 = tmp_4*tmp_51 + tmp_47*tmp_7;
      real_t tmp_55 = -tmp_43 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_56 = tmp_37*tmp_54;
      real_t tmp_57 = 0.2393143352496831*tmp_36;
      real_t tmp_58 = tmp_8*(tmp_10 + 0.5*tmp_9);
      real_t tmp_59 = tmp_2*tmp_58;
      real_t tmp_60 = tmp_8*(0.5*tmp_13 + tmp_14);
      real_t tmp_61 = tmp_4*tmp_60;
      real_t tmp_62 = tmp_59 + tmp_61;
      real_t tmp_63 = tmp_62 - 1.0/3.0;
      real_t tmp_64 = tmp_1*tmp_58;
      real_t tmp_65 = tmp_20*tmp_60;
      real_t tmp_66 = tmp_64 + tmp_65;
      real_t tmp_67 = tmp_66 - 1.0/3.0;
      real_t tmp_68 = tmp_1*tmp_63 + tmp_6*tmp_67;
      real_t tmp_69 = tmp_4*tmp_67 + tmp_63*tmp_7;
      real_t tmp_70 = -tmp_59 - tmp_61 - tmp_64 - tmp_65 + 1;
      real_t tmp_71 = tmp_37*tmp_69;
      real_t tmp_72 = 0.2844444444444445*tmp_36;
      real_t tmp_73 = tmp_8*(tmp_10 + 0.7692346550528415*tmp_9);
      real_t tmp_74 = tmp_2*tmp_73;
      real_t tmp_75 = tmp_8*(0.7692346550528415*tmp_13 + tmp_14);
      real_t tmp_76 = tmp_4*tmp_75;
      real_t tmp_77 = tmp_74 + tmp_76;
      real_t tmp_78 = tmp_77 - 1.0/3.0;
      real_t tmp_79 = tmp_1*tmp_73;
      real_t tmp_80 = tmp_20*tmp_75;
      real_t tmp_81 = tmp_79 + tmp_80;
      real_t tmp_82 = tmp_81 - 1.0/3.0;
      real_t tmp_83 = tmp_1*tmp_78 + tmp_6*tmp_82;
      real_t tmp_84 = tmp_4*tmp_82 + tmp_7*tmp_78;
      real_t tmp_85 = -tmp_74 - tmp_76 - tmp_79 - tmp_80 + 1;
      real_t tmp_86 = tmp_37*tmp_84;
      real_t tmp_87 = 0.2393143352496831*tmp_36;
      real_t tmp_88 = tmp_8*(tmp_10 + 0.95308992296933193*tmp_9);
      real_t tmp_89 = tmp_2*tmp_88;
      real_t tmp_90 = tmp_8*(0.95308992296933193*tmp_13 + tmp_14);
      real_t tmp_91 = tmp_4*tmp_90;
      real_t tmp_92 = tmp_89 + tmp_91;
      real_t tmp_93 = tmp_92 - 1.0/3.0;
      real_t tmp_94 = tmp_1*tmp_88;
      real_t tmp_95 = tmp_20*tmp_90;
      real_t tmp_96 = tmp_94 + tmp_95;
      real_t tmp_97 = tmp_96 - 1.0/3.0;
      real_t tmp_98 = tmp_1*tmp_93 + tmp_6*tmp_97;
      real_t tmp_99 = tmp_4*tmp_97 + tmp_7*tmp_93;
      real_t tmp_100 = -tmp_89 - tmp_91 - tmp_94 - tmp_95 + 1;
      real_t tmp_101 = tmp_37*tmp_99;
      real_t tmp_102 = 0.11846344252809471*tmp_36;
      real_t tmp_103 = 1.0*p_affine_10_0*tmp_26 + 1.0*p_affine_10_1*tmp_33;
      real_t tmp_104 = p_affine_10_1*tmp_26;
      real_t tmp_105 = 1.0*p_affine_10_0*tmp_27 + 1.0*p_affine_10_1*tmp_32;
      real_t tmp_106 = p_affine_10_1*tmp_27;
      real_t a_0_0 = tmp_102*(tmp_100*tmp_101 - tmp_100*tmp_40 - tmp_34*tmp_99 - tmp_53*tmp_98) + tmp_41*(-tmp_24*tmp_29 - tmp_30*tmp_34 + tmp_35*tmp_38 - tmp_35*tmp_40) + tmp_57*(-tmp_34*tmp_54 - tmp_40*tmp_55 - tmp_52*tmp_53 + tmp_55*tmp_56) + tmp_72*(-tmp_34*tmp_69 - tmp_40*tmp_70 - tmp_53*tmp_68 + tmp_70*tmp_71) + tmp_87*(-tmp_34*tmp_84 - tmp_40*tmp_85 - tmp_53*tmp_83 + tmp_85*tmp_86);
      real_t a_0_1 = tmp_102*(tmp_101*tmp_92 - tmp_103*tmp_99 - tmp_104*tmp_98 - tmp_40*tmp_92) + tmp_41*(-tmp_103*tmp_30 + tmp_17*tmp_38 - tmp_17*tmp_40 - tmp_24*tmp_26) + tmp_57*(-tmp_103*tmp_54 - tmp_104*tmp_52 - tmp_40*tmp_46 + tmp_46*tmp_56) + tmp_72*(-tmp_103*tmp_69 - tmp_104*tmp_68 - tmp_40*tmp_62 + tmp_62*tmp_71) + tmp_87*(-tmp_103*tmp_84 - tmp_104*tmp_83 - tmp_40*tmp_77 + tmp_77*tmp_86);
      real_t a_0_2 = tmp_102*(tmp_101*tmp_96 - tmp_105*tmp_99 - tmp_106*tmp_98 - tmp_40*tmp_96) + tmp_41*(-tmp_105*tmp_30 + tmp_22*tmp_38 - tmp_22*tmp_40 - tmp_24*tmp_27) + tmp_57*(-tmp_105*tmp_54 - tmp_106*tmp_52 - tmp_40*tmp_50 + tmp_50*tmp_56) + tmp_72*(-tmp_105*tmp_69 - tmp_106*tmp_68 - tmp_40*tmp_66 + tmp_66*tmp_71) + tmp_87*(-tmp_105*tmp_84 - tmp_106*tmp_83 - tmp_40*tmp_81 + tmp_81*tmp_86);
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
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = tmp_1*tmp_4;
      real_t tmp_6 = p_affine_2_0 + tmp_0;
      real_t tmp_7 = p_affine_1_1 + tmp_3;
      real_t tmp_8 = 1.0 / (tmp_5 - tmp_6*tmp_7);
      real_t tmp_9 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + 0.046910077030668018*tmp_9;
      real_t tmp_11 = tmp_8*(tmp_10 + tmp_3);
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + 0.046910077030668018*tmp_12;
      real_t tmp_14 = tmp_8*(tmp_0 + tmp_13);
      real_t tmp_15 = tmp_11*tmp_2 + tmp_14*tmp_4 - 1.0/3.0;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_1*tmp_11 + tmp_14*tmp_16 - 1.0/3.0;
      real_t tmp_18 = p_affine_10_1*(tmp_1*tmp_15 + tmp_17*tmp_6);
      real_t tmp_19 = -p_affine_3_1;
      real_t tmp_20 = p_affine_5_1 + tmp_19;
      real_t tmp_21 = -p_affine_3_0;
      real_t tmp_22 = p_affine_4_0 + tmp_21;
      real_t tmp_23 = 1.0 / (tmp_20*tmp_22 - (p_affine_4_1 + tmp_19)*(p_affine_5_0 + tmp_21));
      real_t tmp_24 = 0.5*tmp_23;
      real_t tmp_25 = tmp_20*tmp_24;
      real_t tmp_26 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_25 - tmp_27;
      real_t tmp_29 = 1.0*tmp_28;
      real_t tmp_30 = tmp_15*tmp_7 + tmp_17*tmp_4;
      real_t tmp_31 = 1.0*tmp_23;
      real_t tmp_32 = tmp_22*tmp_31;
      real_t tmp_33 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_34 = tmp_31*tmp_33;
      real_t tmp_35 = 1.0*p_affine_10_0*tmp_28 + 1.0*p_affine_10_1*(-tmp_32 - tmp_34);
      real_t tmp_36 = tmp_23*(tmp_10 + tmp_19);
      real_t tmp_37 = tmp_22*tmp_36;
      real_t tmp_38 = tmp_33*tmp_36;
      real_t tmp_39 = tmp_23*(tmp_13 + tmp_21);
      real_t tmp_40 = tmp_20*tmp_39;
      real_t tmp_41 = tmp_26*tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_44 = 12/tmp_43;
      real_t tmp_45 = tmp_30*tmp_44;
      real_t tmp_46 = 1.0*tmp_8;
      real_t tmp_47 = 0.5*tmp_8;
      real_t tmp_48 = tmp_1*tmp_47;
      real_t tmp_49 = tmp_4*tmp_47;
      real_t tmp_50 = p_affine_10_0*(tmp_16*tmp_49 + tmp_2*tmp_48 + tmp_48*tmp_6 + tmp_49*tmp_7) + p_affine_10_1*(tmp_2*tmp_46*tmp_7 + tmp_46*tmp_5);
      real_t tmp_51 = 0.11846344252809471*tmp_43;
      real_t tmp_52 = p_affine_6_1 + 0.23076534494715845*tmp_9;
      real_t tmp_53 = tmp_8*(tmp_3 + tmp_52);
      real_t tmp_54 = p_affine_6_0 + 0.23076534494715845*tmp_12;
      real_t tmp_55 = tmp_8*(tmp_0 + tmp_54);
      real_t tmp_56 = tmp_2*tmp_53 + tmp_4*tmp_55 - 1.0/3.0;
      real_t tmp_57 = tmp_1*tmp_53 + tmp_16*tmp_55 - 1.0/3.0;
      real_t tmp_58 = tmp_1*tmp_56 + tmp_57*tmp_6;
      real_t tmp_59 = p_affine_10_1*tmp_29;
      real_t tmp_60 = tmp_4*tmp_57 + tmp_56*tmp_7;
      real_t tmp_61 = tmp_23*(tmp_19 + tmp_52);
      real_t tmp_62 = tmp_22*tmp_61;
      real_t tmp_63 = tmp_33*tmp_61;
      real_t tmp_64 = tmp_23*(tmp_21 + tmp_54);
      real_t tmp_65 = tmp_20*tmp_64;
      real_t tmp_66 = tmp_26*tmp_64;
      real_t tmp_67 = -tmp_62 - tmp_63 - tmp_65 - tmp_66 + 1;
      real_t tmp_68 = tmp_44*tmp_60;
      real_t tmp_69 = 0.2393143352496831*tmp_43;
      real_t tmp_70 = p_affine_6_1 + 0.5*tmp_9;
      real_t tmp_71 = tmp_8*(tmp_3 + tmp_70);
      real_t tmp_72 = p_affine_6_0 + 0.5*tmp_12;
      real_t tmp_73 = tmp_8*(tmp_0 + tmp_72);
      real_t tmp_74 = tmp_2*tmp_71 + tmp_4*tmp_73 - 1.0/3.0;
      real_t tmp_75 = tmp_1*tmp_71 + tmp_16*tmp_73 - 1.0/3.0;
      real_t tmp_76 = tmp_1*tmp_74 + tmp_6*tmp_75;
      real_t tmp_77 = tmp_4*tmp_75 + tmp_7*tmp_74;
      real_t tmp_78 = tmp_23*(tmp_19 + tmp_70);
      real_t tmp_79 = tmp_22*tmp_78;
      real_t tmp_80 = tmp_33*tmp_78;
      real_t tmp_81 = tmp_23*(tmp_21 + tmp_72);
      real_t tmp_82 = tmp_20*tmp_81;
      real_t tmp_83 = tmp_26*tmp_81;
      real_t tmp_84 = -tmp_79 - tmp_80 - tmp_82 - tmp_83 + 1;
      real_t tmp_85 = tmp_44*tmp_77;
      real_t tmp_86 = 0.2844444444444445*tmp_43;
      real_t tmp_87 = p_affine_6_1 + 0.7692346550528415*tmp_9;
      real_t tmp_88 = tmp_8*(tmp_3 + tmp_87);
      real_t tmp_89 = p_affine_6_0 + 0.7692346550528415*tmp_12;
      real_t tmp_90 = tmp_8*(tmp_0 + tmp_89);
      real_t tmp_91 = tmp_2*tmp_88 + tmp_4*tmp_90 - 1.0/3.0;
      real_t tmp_92 = tmp_1*tmp_88 + tmp_16*tmp_90 - 1.0/3.0;
      real_t tmp_93 = tmp_1*tmp_91 + tmp_6*tmp_92;
      real_t tmp_94 = tmp_4*tmp_92 + tmp_7*tmp_91;
      real_t tmp_95 = tmp_23*(tmp_19 + tmp_87);
      real_t tmp_96 = tmp_22*tmp_95;
      real_t tmp_97 = tmp_33*tmp_95;
      real_t tmp_98 = tmp_23*(tmp_21 + tmp_89);
      real_t tmp_99 = tmp_20*tmp_98;
      real_t tmp_100 = tmp_26*tmp_98;
      real_t tmp_101 = -tmp_100 - tmp_96 - tmp_97 - tmp_99 + 1;
      real_t tmp_102 = tmp_44*tmp_94;
      real_t tmp_103 = 0.2393143352496831*tmp_43;
      real_t tmp_104 = p_affine_6_1 + 0.95308992296933193*tmp_9;
      real_t tmp_105 = tmp_8*(tmp_104 + tmp_3);
      real_t tmp_106 = p_affine_6_0 + 0.95308992296933193*tmp_12;
      real_t tmp_107 = tmp_8*(tmp_0 + tmp_106);
      real_t tmp_108 = tmp_105*tmp_2 + tmp_107*tmp_4 - 1.0/3.0;
      real_t tmp_109 = tmp_1*tmp_105 + tmp_107*tmp_16 - 1.0/3.0;
      real_t tmp_110 = tmp_1*tmp_108 + tmp_109*tmp_6;
      real_t tmp_111 = tmp_108*tmp_7 + tmp_109*tmp_4;
      real_t tmp_112 = tmp_23*(tmp_104 + tmp_19);
      real_t tmp_113 = tmp_112*tmp_22;
      real_t tmp_114 = tmp_112*tmp_33;
      real_t tmp_115 = tmp_23*(tmp_106 + tmp_21);
      real_t tmp_116 = tmp_115*tmp_20;
      real_t tmp_117 = tmp_115*tmp_26;
      real_t tmp_118 = -tmp_113 - tmp_114 - tmp_116 - tmp_117 + 1;
      real_t tmp_119 = tmp_111*tmp_44;
      real_t tmp_120 = 0.11846344252809471*tmp_43;
      real_t tmp_121 = 1.0*p_affine_10_0*tmp_25 + 1.0*p_affine_10_1*tmp_34;
      real_t tmp_122 = tmp_38 + tmp_40;
      real_t tmp_123 = p_affine_10_1*tmp_25;
      real_t tmp_124 = tmp_63 + tmp_65;
      real_t tmp_125 = tmp_80 + tmp_82;
      real_t tmp_126 = tmp_97 + tmp_99;
      real_t tmp_127 = tmp_114 + tmp_116;
      real_t tmp_128 = 1.0*p_affine_10_0*tmp_27 + 1.0*p_affine_10_1*tmp_32;
      real_t tmp_129 = tmp_37 + tmp_41;
      real_t tmp_130 = p_affine_10_1*tmp_27;
      real_t tmp_131 = tmp_62 + tmp_66;
      real_t tmp_132 = tmp_79 + tmp_83;
      real_t tmp_133 = tmp_100 + tmp_96;
      real_t tmp_134 = tmp_113 + tmp_117;
      real_t a_0_0 = tmp_103*(-tmp_101*tmp_102 + tmp_101*tmp_50 - tmp_35*tmp_94 - tmp_59*tmp_93) + tmp_120*(-tmp_110*tmp_59 - tmp_111*tmp_35 - tmp_118*tmp_119 + tmp_118*tmp_50) + tmp_51*(-tmp_18*tmp_29 - tmp_30*tmp_35 - tmp_42*tmp_45 + tmp_42*tmp_50) + tmp_69*(-tmp_35*tmp_60 + tmp_50*tmp_67 - tmp_58*tmp_59 - tmp_67*tmp_68) + tmp_86*(-tmp_35*tmp_77 + tmp_50*tmp_84 - tmp_59*tmp_76 - tmp_84*tmp_85);
      real_t a_0_1 = tmp_103*(-tmp_102*tmp_126 - tmp_121*tmp_94 - tmp_123*tmp_93 + tmp_126*tmp_50) + tmp_120*(-tmp_110*tmp_123 - tmp_111*tmp_121 - tmp_119*tmp_127 + tmp_127*tmp_50) + tmp_51*(-tmp_121*tmp_30 - tmp_122*tmp_45 + tmp_122*tmp_50 - tmp_18*tmp_25) + tmp_69*(-tmp_121*tmp_60 - tmp_123*tmp_58 + tmp_124*tmp_50 - tmp_124*tmp_68) + tmp_86*(-tmp_121*tmp_77 - tmp_123*tmp_76 + tmp_125*tmp_50 - tmp_125*tmp_85);
      real_t a_0_2 = tmp_103*(-tmp_102*tmp_133 - tmp_128*tmp_94 - tmp_130*tmp_93 + tmp_133*tmp_50) + tmp_120*(-tmp_110*tmp_130 - tmp_111*tmp_128 - tmp_119*tmp_134 + tmp_134*tmp_50) + tmp_51*(-tmp_128*tmp_30 - tmp_129*tmp_45 + tmp_129*tmp_50 - tmp_18*tmp_27) + tmp_69*(-tmp_128*tmp_60 - tmp_130*tmp_58 + tmp_131*tmp_50 - tmp_131*tmp_68) + tmp_86*(-tmp_128*tmp_77 - tmp_130*tmp_76 + tmp_132*tmp_50 - tmp_132*tmp_85);
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
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_3;
      real_t tmp_7 = 1.0 / (tmp_1*tmp_4 - tmp_5*tmp_6);
      real_t tmp_8 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9 = p_affine_6_1 + tmp_3;
      real_t tmp_10 = tmp_7*(0.046910077030668018*tmp_8 + tmp_9);
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_7*(0.046910077030668018*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_10*tmp_2 + tmp_13*tmp_4 - 1.0/3.0;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_1*tmp_10 + tmp_13*tmp_15 - 1.0/3.0;
      real_t tmp_17 = p_affine_10_1*(tmp_1*tmp_14 + tmp_16*tmp_5);
      real_t tmp_18 = 0.5*tmp_7;
      real_t tmp_19 = tmp_18*tmp_4;
      real_t tmp_20 = tmp_15*tmp_18;
      real_t tmp_21 = -tmp_19 - tmp_20;
      real_t tmp_22 = 2*tmp_21;
      real_t tmp_23 = tmp_14*tmp_6 + tmp_16*tmp_4;
      real_t tmp_24 = 1.0*tmp_7;
      real_t tmp_25 = tmp_1*tmp_24;
      real_t tmp_26 = tmp_2*tmp_24;
      real_t tmp_27 = 2*p_affine_10_0*tmp_21 + 2*p_affine_10_1*(-tmp_25 - tmp_26);
      real_t tmp_28 = std::abs(std::pow((tmp_11*tmp_11) + (tmp_8*tmp_8), 1.0/2.0));
      real_t tmp_29 = 0.11846344252809471*tmp_28;
      real_t tmp_30 = tmp_7*(0.23076534494715845*tmp_8 + tmp_9);
      real_t tmp_31 = tmp_7*(0.23076534494715845*tmp_11 + tmp_12);
      real_t tmp_32 = tmp_2*tmp_30 + tmp_31*tmp_4 - 1.0/3.0;
      real_t tmp_33 = tmp_1*tmp_30 + tmp_15*tmp_31 - 1.0/3.0;
      real_t tmp_34 = tmp_1*tmp_32 + tmp_33*tmp_5;
      real_t tmp_35 = p_affine_10_1*tmp_22;
      real_t tmp_36 = tmp_32*tmp_6 + tmp_33*tmp_4;
      real_t tmp_37 = 0.2393143352496831*tmp_28;
      real_t tmp_38 = tmp_7*(0.5*tmp_8 + tmp_9);
      real_t tmp_39 = tmp_7*(0.5*tmp_11 + tmp_12);
      real_t tmp_40 = tmp_2*tmp_38 + tmp_39*tmp_4 - 1.0/3.0;
      real_t tmp_41 = tmp_1*tmp_38 + tmp_15*tmp_39 - 1.0/3.0;
      real_t tmp_42 = tmp_1*tmp_40 + tmp_41*tmp_5;
      real_t tmp_43 = tmp_4*tmp_41 + tmp_40*tmp_6;
      real_t tmp_44 = 0.2844444444444445*tmp_28;
      real_t tmp_45 = tmp_7*(0.7692346550528415*tmp_8 + tmp_9);
      real_t tmp_46 = tmp_7*(0.7692346550528415*tmp_11 + tmp_12);
      real_t tmp_47 = tmp_2*tmp_45 + tmp_4*tmp_46 - 1.0/3.0;
      real_t tmp_48 = tmp_1*tmp_45 + tmp_15*tmp_46 - 1.0/3.0;
      real_t tmp_49 = tmp_1*tmp_47 + tmp_48*tmp_5;
      real_t tmp_50 = tmp_4*tmp_48 + tmp_47*tmp_6;
      real_t tmp_51 = 0.2393143352496831*tmp_28;
      real_t tmp_52 = tmp_7*(0.95308992296933193*tmp_8 + tmp_9);
      real_t tmp_53 = tmp_7*(0.95308992296933193*tmp_11 + tmp_12);
      real_t tmp_54 = tmp_2*tmp_52 + tmp_4*tmp_53 - 1.0/3.0;
      real_t tmp_55 = tmp_1*tmp_52 + tmp_15*tmp_53 - 1.0/3.0;
      real_t tmp_56 = tmp_1*tmp_54 + tmp_5*tmp_55;
      real_t tmp_57 = tmp_4*tmp_55 + tmp_54*tmp_6;
      real_t tmp_58 = 0.11846344252809471*tmp_28;
      real_t tmp_59 = tmp_24*tmp_4;
      real_t tmp_60 = 2*p_affine_10_0*tmp_19 + 2*p_affine_10_1*tmp_26;
      real_t tmp_61 = p_affine_10_1*tmp_59;
      real_t tmp_62 = tmp_15*tmp_24;
      real_t tmp_63 = 2*p_affine_10_0*tmp_20 + 2*p_affine_10_1*tmp_25;
      real_t tmp_64 = p_affine_10_1*tmp_62;
      real_t a_0_0 = tmp_29*(-tmp_17*tmp_22 - tmp_23*tmp_27) + tmp_37*(-tmp_27*tmp_36 - tmp_34*tmp_35) + tmp_44*(-tmp_27*tmp_43 - tmp_35*tmp_42) + tmp_51*(-tmp_27*tmp_50 - tmp_35*tmp_49) + tmp_58*(-tmp_27*tmp_57 - tmp_35*tmp_56);
      real_t a_0_1 = tmp_29*(-tmp_17*tmp_59 - tmp_23*tmp_60) + tmp_37*(-tmp_34*tmp_61 - tmp_36*tmp_60) + tmp_44*(-tmp_42*tmp_61 - tmp_43*tmp_60) + tmp_51*(-tmp_49*tmp_61 - tmp_50*tmp_60) + tmp_58*(-tmp_56*tmp_61 - tmp_57*tmp_60);
      real_t a_0_2 = tmp_29*(-tmp_17*tmp_62 - tmp_23*tmp_63) + tmp_37*(-tmp_34*tmp_64 - tmp_36*tmp_63) + tmp_44*(-tmp_42*tmp_64 - tmp_43*tmp_63) + tmp_51*(-tmp_49*tmp_64 - tmp_50*tmp_63) + tmp_58*(-tmp_56*tmp_64 - tmp_57*tmp_63);
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




class EGConstEpsilonFormP1E_1 : public hyteg::dg::DGForm2D
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
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = 1.0 / (tmp_4 - tmp_5*tmp_6);
      real_t tmp_8 = 1.0*tmp_7;
      real_t tmp_9 = tmp_1*tmp_8;
      real_t tmp_10 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_11 = tmp_10*tmp_8;
      real_t tmp_12 = 2.0*tmp_7;
      real_t tmp_13 = tmp_10*tmp_12*tmp_6 + tmp_12*tmp_4;
      real_t tmp_14 = 0.5*tmp_7;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_1*tmp_7;
      real_t tmp_17 = tmp_3*tmp_7;
      real_t tmp_18 = tmp_10*tmp_16 + tmp_15*tmp_17 + tmp_16*tmp_5 + tmp_17*tmp_6;
      real_t tmp_19 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_20 = tmp_19*(tmp_13*(-tmp_11 - tmp_9) + 2*tmp_18*(-tmp_14*tmp_15 - tmp_14*tmp_3));
      real_t tmp_21 = tmp_18*tmp_8;
      real_t tmp_22 = tmp_19*(tmp_11*tmp_13 + tmp_21*tmp_3);
      real_t tmp_23 = tmp_19*(tmp_13*tmp_9 + tmp_15*tmp_21);
      real_t a_0_0 = 0.5*tmp_20;
      real_t a_1_0 = 0.5*tmp_22;
      real_t a_2_0 = 0.5*tmp_23;
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
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = tmp_1*tmp_4;
      real_t tmp_6 = p_affine_2_0 + tmp_0;
      real_t tmp_7 = p_affine_1_1 + tmp_3;
      real_t tmp_8 = 1.0 / (tmp_5 - tmp_6*tmp_7);
      real_t tmp_9 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + tmp_3;
      real_t tmp_11 = tmp_8*(tmp_10 + 0.046910077030668018*tmp_9);
      real_t tmp_12 = tmp_11*tmp_2;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_0;
      real_t tmp_15 = tmp_8*(0.046910077030668018*tmp_13 + tmp_14);
      real_t tmp_16 = tmp_15*tmp_4;
      real_t tmp_17 = tmp_12 + tmp_16;
      real_t tmp_18 = tmp_17 - 1.0/3.0;
      real_t tmp_19 = tmp_1*tmp_11;
      real_t tmp_20 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_21 = tmp_15*tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_22 - 1.0/3.0;
      real_t tmp_24 = p_affine_10_1*(tmp_1*tmp_18 + tmp_23*tmp_6);
      real_t tmp_25 = 0.5*tmp_8;
      real_t tmp_26 = tmp_25*tmp_4;
      real_t tmp_27 = tmp_20*tmp_25;
      real_t tmp_28 = -tmp_26 - tmp_27;
      real_t tmp_29 = 1.0*tmp_28;
      real_t tmp_30 = tmp_18*tmp_7 + tmp_23*tmp_4;
      real_t tmp_31 = 1.0*tmp_8;
      real_t tmp_32 = tmp_1*tmp_31;
      real_t tmp_33 = tmp_2*tmp_31;
      real_t tmp_34 = 1.0*p_affine_10_0*tmp_28 + 1.0*p_affine_10_1*(-tmp_32 - tmp_33);
      real_t tmp_35 = -tmp_12 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_36 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_37 = 12/tmp_36;
      real_t tmp_38 = tmp_30*tmp_37;
      real_t tmp_39 = tmp_1*tmp_25;
      real_t tmp_40 = 1.0*p_affine_10_0*(tmp_2*tmp_39 + tmp_26*tmp_7 + tmp_27*tmp_4 + tmp_39*tmp_6) + 1.0*p_affine_10_1*(tmp_31*tmp_5 + tmp_33*tmp_7);
      real_t tmp_41 = 0.11846344252809471*tmp_36;
      real_t tmp_42 = tmp_8*(tmp_10 + 0.23076534494715845*tmp_9);
      real_t tmp_43 = tmp_2*tmp_42;
      real_t tmp_44 = tmp_8*(0.23076534494715845*tmp_13 + tmp_14);
      real_t tmp_45 = tmp_4*tmp_44;
      real_t tmp_46 = tmp_43 + tmp_45;
      real_t tmp_47 = tmp_46 - 1.0/3.0;
      real_t tmp_48 = tmp_1*tmp_42;
      real_t tmp_49 = tmp_20*tmp_44;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_50 - 1.0/3.0;
      real_t tmp_52 = tmp_1*tmp_47 + tmp_51*tmp_6;
      real_t tmp_53 = p_affine_10_1*tmp_29;
      real_t tmp_54 = tmp_4*tmp_51 + tmp_47*tmp_7;
      real_t tmp_55 = -tmp_43 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_56 = tmp_37*tmp_54;
      real_t tmp_57 = 0.2393143352496831*tmp_36;
      real_t tmp_58 = tmp_8*(tmp_10 + 0.5*tmp_9);
      real_t tmp_59 = tmp_2*tmp_58;
      real_t tmp_60 = tmp_8*(0.5*tmp_13 + tmp_14);
      real_t tmp_61 = tmp_4*tmp_60;
      real_t tmp_62 = tmp_59 + tmp_61;
      real_t tmp_63 = tmp_62 - 1.0/3.0;
      real_t tmp_64 = tmp_1*tmp_58;
      real_t tmp_65 = tmp_20*tmp_60;
      real_t tmp_66 = tmp_64 + tmp_65;
      real_t tmp_67 = tmp_66 - 1.0/3.0;
      real_t tmp_68 = tmp_1*tmp_63 + tmp_6*tmp_67;
      real_t tmp_69 = tmp_4*tmp_67 + tmp_63*tmp_7;
      real_t tmp_70 = -tmp_59 - tmp_61 - tmp_64 - tmp_65 + 1;
      real_t tmp_71 = tmp_37*tmp_69;
      real_t tmp_72 = 0.2844444444444445*tmp_36;
      real_t tmp_73 = tmp_8*(tmp_10 + 0.7692346550528415*tmp_9);
      real_t tmp_74 = tmp_2*tmp_73;
      real_t tmp_75 = tmp_8*(0.7692346550528415*tmp_13 + tmp_14);
      real_t tmp_76 = tmp_4*tmp_75;
      real_t tmp_77 = tmp_74 + tmp_76;
      real_t tmp_78 = tmp_77 - 1.0/3.0;
      real_t tmp_79 = tmp_1*tmp_73;
      real_t tmp_80 = tmp_20*tmp_75;
      real_t tmp_81 = tmp_79 + tmp_80;
      real_t tmp_82 = tmp_81 - 1.0/3.0;
      real_t tmp_83 = tmp_1*tmp_78 + tmp_6*tmp_82;
      real_t tmp_84 = tmp_4*tmp_82 + tmp_7*tmp_78;
      real_t tmp_85 = -tmp_74 - tmp_76 - tmp_79 - tmp_80 + 1;
      real_t tmp_86 = tmp_37*tmp_84;
      real_t tmp_87 = 0.2393143352496831*tmp_36;
      real_t tmp_88 = tmp_8*(tmp_10 + 0.95308992296933193*tmp_9);
      real_t tmp_89 = tmp_2*tmp_88;
      real_t tmp_90 = tmp_8*(0.95308992296933193*tmp_13 + tmp_14);
      real_t tmp_91 = tmp_4*tmp_90;
      real_t tmp_92 = tmp_89 + tmp_91;
      real_t tmp_93 = tmp_92 - 1.0/3.0;
      real_t tmp_94 = tmp_1*tmp_88;
      real_t tmp_95 = tmp_20*tmp_90;
      real_t tmp_96 = tmp_94 + tmp_95;
      real_t tmp_97 = tmp_96 - 1.0/3.0;
      real_t tmp_98 = tmp_1*tmp_93 + tmp_6*tmp_97;
      real_t tmp_99 = tmp_4*tmp_97 + tmp_7*tmp_93;
      real_t tmp_100 = -tmp_89 - tmp_91 - tmp_94 - tmp_95 + 1;
      real_t tmp_101 = tmp_37*tmp_99;
      real_t tmp_102 = 0.11846344252809471*tmp_36;
      real_t tmp_103 = 1.0*p_affine_10_0*tmp_26 + 1.0*p_affine_10_1*tmp_33;
      real_t tmp_104 = p_affine_10_1*tmp_26;
      real_t tmp_105 = 1.0*p_affine_10_0*tmp_27 + 1.0*p_affine_10_1*tmp_32;
      real_t tmp_106 = p_affine_10_1*tmp_27;
      real_t a_0_0 = tmp_102*(tmp_100*tmp_101 - tmp_100*tmp_40 - tmp_34*tmp_99 - tmp_53*tmp_98) + tmp_41*(-tmp_24*tmp_29 - tmp_30*tmp_34 + tmp_35*tmp_38 - tmp_35*tmp_40) + tmp_57*(-tmp_34*tmp_54 - tmp_40*tmp_55 - tmp_52*tmp_53 + tmp_55*tmp_56) + tmp_72*(-tmp_34*tmp_69 - tmp_40*tmp_70 - tmp_53*tmp_68 + tmp_70*tmp_71) + tmp_87*(-tmp_34*tmp_84 - tmp_40*tmp_85 - tmp_53*tmp_83 + tmp_85*tmp_86);
      real_t a_1_0 = tmp_102*(tmp_101*tmp_92 - tmp_103*tmp_99 - tmp_104*tmp_98 - tmp_40*tmp_92) + tmp_41*(-tmp_103*tmp_30 + tmp_17*tmp_38 - tmp_17*tmp_40 - tmp_24*tmp_26) + tmp_57*(-tmp_103*tmp_54 - tmp_104*tmp_52 - tmp_40*tmp_46 + tmp_46*tmp_56) + tmp_72*(-tmp_103*tmp_69 - tmp_104*tmp_68 - tmp_40*tmp_62 + tmp_62*tmp_71) + tmp_87*(-tmp_103*tmp_84 - tmp_104*tmp_83 - tmp_40*tmp_77 + tmp_77*tmp_86);
      real_t a_2_0 = tmp_102*(tmp_101*tmp_96 - tmp_105*tmp_99 - tmp_106*tmp_98 - tmp_40*tmp_96) + tmp_41*(-tmp_105*tmp_30 + tmp_22*tmp_38 - tmp_22*tmp_40 - tmp_24*tmp_27) + tmp_57*(-tmp_105*tmp_54 - tmp_106*tmp_52 - tmp_40*tmp_50 + tmp_50*tmp_56) + tmp_72*(-tmp_105*tmp_69 - tmp_106*tmp_68 - tmp_40*tmp_66 + tmp_66*tmp_71) + tmp_87*(-tmp_105*tmp_84 - tmp_106*tmp_83 - tmp_40*tmp_81 + tmp_81*tmp_86);
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
      real_t tmp_7 = p_affine_4_1 + tmp_3;
      real_t tmp_8 = 1.0 / (tmp_5 - tmp_6*tmp_7);
      real_t tmp_9 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + 0.046910077030668018*tmp_9;
      real_t tmp_11 = tmp_8*(tmp_10 + tmp_3);
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + 0.046910077030668018*tmp_12;
      real_t tmp_14 = tmp_8*(tmp_0 + tmp_13);
      real_t tmp_15 = tmp_11*tmp_2 + tmp_14*tmp_4 - 1.0/3.0;
      real_t tmp_16 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_17 = tmp_1*tmp_11 + tmp_14*tmp_16 - 1.0/3.0;
      real_t tmp_18 = tmp_1*tmp_15 + tmp_17*tmp_6;
      real_t tmp_19 = -p_affine_0_1;
      real_t tmp_20 = p_affine_2_1 + tmp_19;
      real_t tmp_21 = -p_affine_0_0;
      real_t tmp_22 = p_affine_1_0 + tmp_21;
      real_t tmp_23 = 1.0 / (tmp_20*tmp_22 - (p_affine_1_1 + tmp_19)*(p_affine_2_0 + tmp_21));
      real_t tmp_24 = 0.5*tmp_23;
      real_t tmp_25 = tmp_20*tmp_24;
      real_t tmp_26 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_25 - tmp_27;
      real_t tmp_29 = p_affine_10_1*tmp_28;
      real_t tmp_30 = 1.0*tmp_23;
      real_t tmp_31 = tmp_22*tmp_30;
      real_t tmp_32 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_33 = tmp_30*tmp_32;
      real_t tmp_34 = p_affine_10_0*tmp_28 + p_affine_10_1*(-tmp_31 - tmp_33);
      real_t tmp_35 = tmp_15*tmp_7 + tmp_17*tmp_4;
      real_t tmp_36 = tmp_23*(tmp_10 + tmp_19);
      real_t tmp_37 = tmp_22*tmp_36;
      real_t tmp_38 = tmp_32*tmp_36;
      real_t tmp_39 = tmp_23*(tmp_13 + tmp_21);
      real_t tmp_40 = tmp_20*tmp_39;
      real_t tmp_41 = tmp_26*tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_44 = 12/tmp_43;
      real_t tmp_45 = tmp_35*tmp_44;
      real_t tmp_46 = 1.0*tmp_8;
      real_t tmp_47 = 0.5*tmp_8;
      real_t tmp_48 = tmp_1*tmp_47;
      real_t tmp_49 = tmp_4*tmp_47;
      real_t tmp_50 = 1.0*p_affine_10_0*(tmp_16*tmp_49 + tmp_2*tmp_48 + tmp_48*tmp_6 + tmp_49*tmp_7) + 1.0*p_affine_10_1*(tmp_2*tmp_46*tmp_7 + tmp_46*tmp_5);
      real_t tmp_51 = 0.11846344252809471*tmp_43;
      real_t tmp_52 = p_affine_6_1 + 0.23076534494715845*tmp_9;
      real_t tmp_53 = tmp_8*(tmp_3 + tmp_52);
      real_t tmp_54 = p_affine_6_0 + 0.23076534494715845*tmp_12;
      real_t tmp_55 = tmp_8*(tmp_0 + tmp_54);
      real_t tmp_56 = tmp_2*tmp_53 + tmp_4*tmp_55 - 1.0/3.0;
      real_t tmp_57 = tmp_1*tmp_53 + tmp_16*tmp_55 - 1.0/3.0;
      real_t tmp_58 = tmp_1*tmp_56 + tmp_57*tmp_6;
      real_t tmp_59 = tmp_4*tmp_57 + tmp_56*tmp_7;
      real_t tmp_60 = tmp_23*(tmp_19 + tmp_52);
      real_t tmp_61 = tmp_22*tmp_60;
      real_t tmp_62 = tmp_32*tmp_60;
      real_t tmp_63 = tmp_23*(tmp_21 + tmp_54);
      real_t tmp_64 = tmp_20*tmp_63;
      real_t tmp_65 = tmp_26*tmp_63;
      real_t tmp_66 = -tmp_61 - tmp_62 - tmp_64 - tmp_65 + 1;
      real_t tmp_67 = tmp_44*tmp_59;
      real_t tmp_68 = 0.2393143352496831*tmp_43;
      real_t tmp_69 = p_affine_6_1 + 0.5*tmp_9;
      real_t tmp_70 = tmp_8*(tmp_3 + tmp_69);
      real_t tmp_71 = p_affine_6_0 + 0.5*tmp_12;
      real_t tmp_72 = tmp_8*(tmp_0 + tmp_71);
      real_t tmp_73 = tmp_2*tmp_70 + tmp_4*tmp_72 - 1.0/3.0;
      real_t tmp_74 = tmp_1*tmp_70 + tmp_16*tmp_72 - 1.0/3.0;
      real_t tmp_75 = tmp_1*tmp_73 + tmp_6*tmp_74;
      real_t tmp_76 = tmp_4*tmp_74 + tmp_7*tmp_73;
      real_t tmp_77 = tmp_23*(tmp_19 + tmp_69);
      real_t tmp_78 = tmp_22*tmp_77;
      real_t tmp_79 = tmp_32*tmp_77;
      real_t tmp_80 = tmp_23*(tmp_21 + tmp_71);
      real_t tmp_81 = tmp_20*tmp_80;
      real_t tmp_82 = tmp_26*tmp_80;
      real_t tmp_83 = -tmp_78 - tmp_79 - tmp_81 - tmp_82 + 1;
      real_t tmp_84 = tmp_44*tmp_76;
      real_t tmp_85 = 0.2844444444444445*tmp_43;
      real_t tmp_86 = p_affine_6_1 + 0.7692346550528415*tmp_9;
      real_t tmp_87 = tmp_8*(tmp_3 + tmp_86);
      real_t tmp_88 = p_affine_6_0 + 0.7692346550528415*tmp_12;
      real_t tmp_89 = tmp_8*(tmp_0 + tmp_88);
      real_t tmp_90 = tmp_2*tmp_87 + tmp_4*tmp_89 - 1.0/3.0;
      real_t tmp_91 = tmp_1*tmp_87 + tmp_16*tmp_89 - 1.0/3.0;
      real_t tmp_92 = tmp_1*tmp_90 + tmp_6*tmp_91;
      real_t tmp_93 = tmp_4*tmp_91 + tmp_7*tmp_90;
      real_t tmp_94 = tmp_23*(tmp_19 + tmp_86);
      real_t tmp_95 = tmp_22*tmp_94;
      real_t tmp_96 = tmp_32*tmp_94;
      real_t tmp_97 = tmp_23*(tmp_21 + tmp_88);
      real_t tmp_98 = tmp_20*tmp_97;
      real_t tmp_99 = tmp_26*tmp_97;
      real_t tmp_100 = -tmp_95 - tmp_96 - tmp_98 - tmp_99 + 1;
      real_t tmp_101 = tmp_44*tmp_93;
      real_t tmp_102 = 0.2393143352496831*tmp_43;
      real_t tmp_103 = p_affine_6_1 + 0.95308992296933193*tmp_9;
      real_t tmp_104 = tmp_8*(tmp_103 + tmp_3);
      real_t tmp_105 = p_affine_6_0 + 0.95308992296933193*tmp_12;
      real_t tmp_106 = tmp_8*(tmp_0 + tmp_105);
      real_t tmp_107 = tmp_104*tmp_2 + tmp_106*tmp_4 - 1.0/3.0;
      real_t tmp_108 = tmp_1*tmp_104 + tmp_106*tmp_16 - 1.0/3.0;
      real_t tmp_109 = tmp_1*tmp_107 + tmp_108*tmp_6;
      real_t tmp_110 = tmp_107*tmp_7 + tmp_108*tmp_4;
      real_t tmp_111 = tmp_23*(tmp_103 + tmp_19);
      real_t tmp_112 = tmp_111*tmp_22;
      real_t tmp_113 = tmp_111*tmp_32;
      real_t tmp_114 = tmp_23*(tmp_105 + tmp_21);
      real_t tmp_115 = tmp_114*tmp_20;
      real_t tmp_116 = tmp_114*tmp_26;
      real_t tmp_117 = -tmp_112 - tmp_113 - tmp_115 - tmp_116 + 1;
      real_t tmp_118 = tmp_110*tmp_44;
      real_t tmp_119 = 0.11846344252809471*tmp_43;
      real_t tmp_120 = p_affine_10_1*tmp_25;
      real_t tmp_121 = p_affine_10_0*tmp_25 + p_affine_10_1*tmp_33;
      real_t tmp_122 = tmp_38 + tmp_40;
      real_t tmp_123 = tmp_62 + tmp_64;
      real_t tmp_124 = tmp_79 + tmp_81;
      real_t tmp_125 = tmp_96 + tmp_98;
      real_t tmp_126 = tmp_113 + tmp_115;
      real_t tmp_127 = p_affine_10_1*tmp_27;
      real_t tmp_128 = p_affine_10_0*tmp_27 + p_affine_10_1*tmp_31;
      real_t tmp_129 = tmp_37 + tmp_41;
      real_t tmp_130 = tmp_61 + tmp_65;
      real_t tmp_131 = tmp_78 + tmp_82;
      real_t tmp_132 = tmp_95 + tmp_99;
      real_t tmp_133 = tmp_112 + tmp_116;
      real_t a_0_0 = tmp_102*(-tmp_100*tmp_101 - tmp_100*tmp_50 + tmp_29*tmp_92 + tmp_34*tmp_93) + tmp_119*(tmp_109*tmp_29 + tmp_110*tmp_34 - tmp_117*tmp_118 - tmp_117*tmp_50) + tmp_51*(tmp_18*tmp_29 + tmp_34*tmp_35 - tmp_42*tmp_45 - tmp_42*tmp_50) + tmp_68*(tmp_29*tmp_58 + tmp_34*tmp_59 - tmp_50*tmp_66 - tmp_66*tmp_67) + tmp_85*(tmp_29*tmp_75 + tmp_34*tmp_76 - tmp_50*tmp_83 - tmp_83*tmp_84);
      real_t a_1_0 = tmp_102*(-tmp_101*tmp_125 + tmp_120*tmp_92 + tmp_121*tmp_93 - tmp_125*tmp_50) + tmp_119*(tmp_109*tmp_120 + tmp_110*tmp_121 - tmp_118*tmp_126 - tmp_126*tmp_50) + tmp_51*(tmp_120*tmp_18 + tmp_121*tmp_35 - tmp_122*tmp_45 - tmp_122*tmp_50) + tmp_68*(tmp_120*tmp_58 + tmp_121*tmp_59 - tmp_123*tmp_50 - tmp_123*tmp_67) + tmp_85*(tmp_120*tmp_75 + tmp_121*tmp_76 - tmp_124*tmp_50 - tmp_124*tmp_84);
      real_t a_2_0 = tmp_102*(-tmp_101*tmp_132 + tmp_127*tmp_92 + tmp_128*tmp_93 - tmp_132*tmp_50) + tmp_119*(tmp_109*tmp_127 + tmp_110*tmp_128 - tmp_118*tmp_133 - tmp_133*tmp_50) + tmp_51*(tmp_127*tmp_18 + tmp_128*tmp_35 - tmp_129*tmp_45 - tmp_129*tmp_50) + tmp_68*(tmp_127*tmp_58 + tmp_128*tmp_59 - tmp_130*tmp_50 - tmp_130*tmp_67) + tmp_85*(tmp_127*tmp_75 + tmp_128*tmp_76 - tmp_131*tmp_50 - tmp_131*tmp_84);
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
      real_t tmp_4 = p_affine_1_1 + tmp_2;
      real_t tmp_5 = 1.0 / (tmp_1*tmp_3 - tmp_4*(p_affine_2_0 + tmp_0));
      real_t tmp_6 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_7 = p_affine_6_1 + tmp_2;
      real_t tmp_8 = tmp_5*(0.046910077030668018*tmp_6 + tmp_7);
      real_t tmp_9 = tmp_1*tmp_8;
      real_t tmp_10 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_11 = tmp_10*tmp_8;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_0;
      real_t tmp_14 = tmp_5*(0.046910077030668018*tmp_12 + tmp_13);
      real_t tmp_15 = tmp_14*tmp_3;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_14*tmp_16;
      real_t tmp_18 = tmp_11 + tmp_15;
      real_t tmp_19 = tmp_17 + tmp_9;
      real_t tmp_20 = 1.4215613103371365*tmp_3*(tmp_19 - 1.0/3.0) + 1.4215613103371365*tmp_4*(tmp_18 - 1.0/3.0);
      real_t tmp_21 = tmp_5*(0.23076534494715845*tmp_6 + tmp_7);
      real_t tmp_22 = tmp_1*tmp_21;
      real_t tmp_23 = tmp_10*tmp_21;
      real_t tmp_24 = tmp_5*(0.23076534494715845*tmp_12 + tmp_13);
      real_t tmp_25 = tmp_24*tmp_3;
      real_t tmp_26 = tmp_16*tmp_24;
      real_t tmp_27 = tmp_23 + tmp_25;
      real_t tmp_28 = tmp_22 + tmp_26;
      real_t tmp_29 = 2.8717720229961969*tmp_3*(tmp_28 - 1.0/3.0) + 2.8717720229961969*tmp_4*(tmp_27 - 1.0/3.0);
      real_t tmp_30 = tmp_5*(0.5*tmp_6 + tmp_7);
      real_t tmp_31 = tmp_1*tmp_30;
      real_t tmp_32 = tmp_10*tmp_30;
      real_t tmp_33 = tmp_5*(0.5*tmp_12 + tmp_13);
      real_t tmp_34 = tmp_3*tmp_33;
      real_t tmp_35 = tmp_16*tmp_33;
      real_t tmp_36 = tmp_32 + tmp_34;
      real_t tmp_37 = tmp_31 + tmp_35;
      real_t tmp_38 = 3.413333333333334*tmp_3*(tmp_37 - 1.0/3.0) + 3.413333333333334*tmp_4*(tmp_36 - 1.0/3.0);
      real_t tmp_39 = tmp_5*(0.7692346550528415*tmp_6 + tmp_7);
      real_t tmp_40 = tmp_1*tmp_39;
      real_t tmp_41 = tmp_10*tmp_39;
      real_t tmp_42 = tmp_5*(0.7692346550528415*tmp_12 + tmp_13);
      real_t tmp_43 = tmp_3*tmp_42;
      real_t tmp_44 = tmp_16*tmp_42;
      real_t tmp_45 = tmp_41 + tmp_43;
      real_t tmp_46 = tmp_40 + tmp_44;
      real_t tmp_47 = 2.8717720229961969*tmp_3*(tmp_46 - 1.0/3.0) + 2.8717720229961969*tmp_4*(tmp_45 - 1.0/3.0);
      real_t tmp_48 = tmp_5*(0.95308992296933193*tmp_6 + tmp_7);
      real_t tmp_49 = tmp_1*tmp_48;
      real_t tmp_50 = tmp_10*tmp_48;
      real_t tmp_51 = tmp_5*(0.95308992296933193*tmp_12 + tmp_13);
      real_t tmp_52 = tmp_3*tmp_51;
      real_t tmp_53 = tmp_16*tmp_51;
      real_t tmp_54 = tmp_50 + tmp_52;
      real_t tmp_55 = tmp_49 + tmp_53;
      real_t tmp_56 = 1.4215613103371365*tmp_3*(tmp_55 - 1.0/3.0) + 1.4215613103371365*tmp_4*(tmp_54 - 1.0/3.0);
      real_t a_0_0 = tmp_20*(-tmp_11 - tmp_15 - tmp_17 - tmp_9 + 1) + tmp_29*(-tmp_22 - tmp_23 - tmp_25 - tmp_26 + 1) + tmp_38*(-tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1) + tmp_47*(-tmp_40 - tmp_41 - tmp_43 - tmp_44 + 1) + tmp_56*(-tmp_49 - tmp_50 - tmp_52 - tmp_53 + 1);
      real_t a_1_0 = tmp_18*tmp_20 + tmp_27*tmp_29 + tmp_36*tmp_38 + tmp_45*tmp_47 + tmp_54*tmp_56;
      real_t a_2_0 = tmp_19*tmp_20 + tmp_28*tmp_29 + tmp_37*tmp_38 + tmp_46*tmp_47 + tmp_55*tmp_56;
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




class EGConstEpsilonFormEE : public hyteg::dg::DGForm2D
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
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = 1.0 / (tmp_4 - tmp_5*tmp_6);
      real_t tmp_8 = 1.0*tmp_7;
      real_t tmp_9 = tmp_4*tmp_8;
      real_t tmp_10 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_11 = tmp_10*tmp_5;
      real_t tmp_12 = 2.0*tmp_7;
      real_t tmp_13 = tmp_12*tmp_4;
      real_t tmp_14 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_15 = tmp_14*tmp_6;
      real_t tmp_16 = tmp_1*tmp_7;
      real_t tmp_17 = tmp_16*tmp_5;
      real_t tmp_18 = tmp_14*tmp_16;
      real_t tmp_19 = tmp_3*tmp_7;
      real_t tmp_20 = tmp_19*tmp_6;
      real_t tmp_21 = tmp_10*tmp_19;
      real_t tmp_22 = ((tmp_11*tmp_12 + tmp_13)*(tmp_11*tmp_8 + tmp_9) + (tmp_12*tmp_15 + tmp_13)*(tmp_15*tmp_8 + tmp_9) + 2*(0.5*tmp_17 + 0.5*tmp_18 + 0.5*tmp_20 + 0.5*tmp_21)*(tmp_17 + tmp_18 + tmp_20 + tmp_21))*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = 0.5*tmp_22;
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
      real_t tmp_20 = tmp_10*tmp_16 + tmp_18*tmp_7;
      real_t tmp_21 = 12/tmp_2;
      real_t tmp_22 = 1.0*tmp_11;
      real_t tmp_23 = tmp_22*tmp_8;
      real_t tmp_24 = 0.5*tmp_11;
      real_t tmp_25 = tmp_24*tmp_4;
      real_t tmp_26 = tmp_24*tmp_7;
      real_t tmp_27 = tmp_10*tmp_26 + tmp_17*tmp_26 + tmp_25*tmp_5 + tmp_25*tmp_9;
      real_t tmp_28 = 2.0*p_affine_10_0*(tmp_17*tmp_22*tmp_9 + tmp_23) + 2.0*p_affine_10_1*tmp_27;
      real_t tmp_29 = 2.0*p_affine_10_0*tmp_27 + 2.0*p_affine_10_1*(tmp_10*tmp_22*tmp_5 + tmp_23);
      real_t tmp_30 = tmp_11*(0.23076534494715845*tmp_1 + tmp_12);
      real_t tmp_31 = tmp_11*(0.23076534494715845*tmp_0 + tmp_14);
      real_t tmp_32 = tmp_30*tmp_5 + tmp_31*tmp_7 - 1.0/3.0;
      real_t tmp_33 = tmp_17*tmp_31 + tmp_30*tmp_4 - 1.0/3.0;
      real_t tmp_34 = tmp_32*tmp_4 + tmp_33*tmp_9;
      real_t tmp_35 = tmp_10*tmp_32 + tmp_33*tmp_7;
      real_t tmp_36 = tmp_11*(0.5*tmp_1 + tmp_12);
      real_t tmp_37 = tmp_11*(0.5*tmp_0 + tmp_14);
      real_t tmp_38 = tmp_36*tmp_5 + tmp_37*tmp_7 - 1.0/3.0;
      real_t tmp_39 = tmp_17*tmp_37 + tmp_36*tmp_4 - 1.0/3.0;
      real_t tmp_40 = tmp_38*tmp_4 + tmp_39*tmp_9;
      real_t tmp_41 = tmp_10*tmp_38 + tmp_39*tmp_7;
      real_t tmp_42 = tmp_11*(0.7692346550528415*tmp_1 + tmp_12);
      real_t tmp_43 = tmp_11*(0.7692346550528415*tmp_0 + tmp_14);
      real_t tmp_44 = tmp_42*tmp_5 + tmp_43*tmp_7 - 1.0/3.0;
      real_t tmp_45 = tmp_17*tmp_43 + tmp_4*tmp_42 - 1.0/3.0;
      real_t tmp_46 = tmp_4*tmp_44 + tmp_45*tmp_9;
      real_t tmp_47 = tmp_10*tmp_44 + tmp_45*tmp_7;
      real_t tmp_48 = tmp_11*(0.95308992296933193*tmp_1 + tmp_12);
      real_t tmp_49 = tmp_11*(0.95308992296933193*tmp_0 + tmp_14);
      real_t tmp_50 = tmp_48*tmp_5 + tmp_49*tmp_7 - 1.0/3.0;
      real_t tmp_51 = tmp_17*tmp_49 + tmp_4*tmp_48 - 1.0/3.0;
      real_t tmp_52 = tmp_4*tmp_50 + tmp_51*tmp_9;
      real_t tmp_53 = tmp_10*tmp_50 + tmp_51*tmp_7;
      real_t a_0_0 = 0.11846344252809471*tmp_2*(-tmp_19*tmp_28 - tmp_20*tmp_29 + tmp_21*((tmp_19*tmp_19) + (tmp_20*tmp_20))) + 0.2393143352496831*tmp_2*(tmp_21*((tmp_34*tmp_34) + (tmp_35*tmp_35)) - tmp_28*tmp_34 - tmp_29*tmp_35) + 0.2844444444444445*tmp_2*(tmp_21*((tmp_40*tmp_40) + (tmp_41*tmp_41)) - tmp_28*tmp_40 - tmp_29*tmp_41) + 0.2393143352496831*tmp_2*(tmp_21*((tmp_46*tmp_46) + (tmp_47*tmp_47)) - tmp_28*tmp_46 - tmp_29*tmp_47) + 0.11846344252809471*tmp_2*(tmp_21*((tmp_52*tmp_52) + (tmp_53*tmp_53)) - tmp_28*tmp_52 - tmp_29*tmp_53);
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
      real_t tmp_28 = 1.0*tmp_27;
      real_t tmp_29 = tmp_24*tmp_28;
      real_t tmp_30 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_31 = 0.5*tmp_27;
      real_t tmp_32 = tmp_21*tmp_31;
      real_t tmp_33 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_34 = tmp_23*tmp_31;
      real_t tmp_35 = tmp_25*tmp_32 + tmp_26*tmp_34 + tmp_30*tmp_34 + tmp_32*tmp_33;
      real_t tmp_36 = p_affine_10_0*(tmp_25*tmp_28*tmp_30 + tmp_29) + p_affine_10_1*tmp_35;
      real_t tmp_37 = tmp_10*tmp_16 + tmp_18*tmp_7;
      real_t tmp_38 = p_affine_10_0*tmp_35 + p_affine_10_1*(tmp_26*tmp_28*tmp_33 + tmp_29);
      real_t tmp_39 = tmp_27*(tmp_12 + tmp_22);
      real_t tmp_40 = tmp_27*(tmp_14 + tmp_20);
      real_t tmp_41 = tmp_23*tmp_40 + tmp_33*tmp_39 - 1.0/3.0;
      real_t tmp_42 = tmp_21*tmp_39 + tmp_30*tmp_40 - 1.0/3.0;
      real_t tmp_43 = tmp_21*tmp_41 + tmp_25*tmp_42;
      real_t tmp_44 = 1.0*tmp_11;
      real_t tmp_45 = tmp_44*tmp_8;
      real_t tmp_46 = 0.5*tmp_11;
      real_t tmp_47 = tmp_4*tmp_46;
      real_t tmp_48 = tmp_46*tmp_7;
      real_t tmp_49 = tmp_10*tmp_48 + tmp_17*tmp_48 + tmp_47*tmp_5 + tmp_47*tmp_9;
      real_t tmp_50 = 1.0*p_affine_10_0*(tmp_17*tmp_44*tmp_9 + tmp_45) + 1.0*p_affine_10_1*tmp_49;
      real_t tmp_51 = tmp_23*tmp_42 + tmp_26*tmp_41;
      real_t tmp_52 = 1.0*p_affine_10_0*tmp_49 + 1.0*p_affine_10_1*(tmp_10*tmp_44*tmp_5 + tmp_45);
      real_t tmp_53 = 12/tmp_2;
      real_t tmp_54 = p_affine_6_1 + 0.23076534494715845*tmp_1;
      real_t tmp_55 = tmp_11*(tmp_54 + tmp_6);
      real_t tmp_56 = p_affine_6_0 + 0.23076534494715845*tmp_0;
      real_t tmp_57 = tmp_11*(tmp_3 + tmp_56);
      real_t tmp_58 = tmp_5*tmp_55 + tmp_57*tmp_7 - 1.0/3.0;
      real_t tmp_59 = tmp_17*tmp_57 + tmp_4*tmp_55 - 1.0/3.0;
      real_t tmp_60 = tmp_4*tmp_58 + tmp_59*tmp_9;
      real_t tmp_61 = tmp_10*tmp_58 + tmp_59*tmp_7;
      real_t tmp_62 = tmp_27*(tmp_22 + tmp_54);
      real_t tmp_63 = tmp_27*(tmp_20 + tmp_56);
      real_t tmp_64 = tmp_23*tmp_63 + tmp_33*tmp_62 - 1.0/3.0;
      real_t tmp_65 = tmp_21*tmp_62 + tmp_30*tmp_63 - 1.0/3.0;
      real_t tmp_66 = tmp_21*tmp_64 + tmp_25*tmp_65;
      real_t tmp_67 = tmp_23*tmp_65 + tmp_26*tmp_64;
      real_t tmp_68 = p_affine_6_1 + 0.5*tmp_1;
      real_t tmp_69 = tmp_11*(tmp_6 + tmp_68);
      real_t tmp_70 = p_affine_6_0 + 0.5*tmp_0;
      real_t tmp_71 = tmp_11*(tmp_3 + tmp_70);
      real_t tmp_72 = tmp_5*tmp_69 + tmp_7*tmp_71 - 1.0/3.0;
      real_t tmp_73 = tmp_17*tmp_71 + tmp_4*tmp_69 - 1.0/3.0;
      real_t tmp_74 = tmp_4*tmp_72 + tmp_73*tmp_9;
      real_t tmp_75 = tmp_10*tmp_72 + tmp_7*tmp_73;
      real_t tmp_76 = tmp_27*(tmp_22 + tmp_68);
      real_t tmp_77 = tmp_27*(tmp_20 + tmp_70);
      real_t tmp_78 = tmp_23*tmp_77 + tmp_33*tmp_76 - 1.0/3.0;
      real_t tmp_79 = tmp_21*tmp_76 + tmp_30*tmp_77 - 1.0/3.0;
      real_t tmp_80 = tmp_21*tmp_78 + tmp_25*tmp_79;
      real_t tmp_81 = tmp_23*tmp_79 + tmp_26*tmp_78;
      real_t tmp_82 = p_affine_6_1 + 0.7692346550528415*tmp_1;
      real_t tmp_83 = tmp_11*(tmp_6 + tmp_82);
      real_t tmp_84 = p_affine_6_0 + 0.7692346550528415*tmp_0;
      real_t tmp_85 = tmp_11*(tmp_3 + tmp_84);
      real_t tmp_86 = tmp_5*tmp_83 + tmp_7*tmp_85 - 1.0/3.0;
      real_t tmp_87 = tmp_17*tmp_85 + tmp_4*tmp_83 - 1.0/3.0;
      real_t tmp_88 = tmp_4*tmp_86 + tmp_87*tmp_9;
      real_t tmp_89 = tmp_10*tmp_86 + tmp_7*tmp_87;
      real_t tmp_90 = tmp_27*(tmp_22 + tmp_82);
      real_t tmp_91 = tmp_27*(tmp_20 + tmp_84);
      real_t tmp_92 = tmp_23*tmp_91 + tmp_33*tmp_90 - 1.0/3.0;
      real_t tmp_93 = tmp_21*tmp_90 + tmp_30*tmp_91 - 1.0/3.0;
      real_t tmp_94 = tmp_21*tmp_92 + tmp_25*tmp_93;
      real_t tmp_95 = tmp_23*tmp_93 + tmp_26*tmp_92;
      real_t tmp_96 = p_affine_6_1 + 0.95308992296933193*tmp_1;
      real_t tmp_97 = tmp_11*(tmp_6 + tmp_96);
      real_t tmp_98 = p_affine_6_0 + 0.95308992296933193*tmp_0;
      real_t tmp_99 = tmp_11*(tmp_3 + tmp_98);
      real_t tmp_100 = tmp_5*tmp_97 + tmp_7*tmp_99 - 1.0/3.0;
      real_t tmp_101 = tmp_17*tmp_99 + tmp_4*tmp_97 - 1.0/3.0;
      real_t tmp_102 = tmp_100*tmp_4 + tmp_101*tmp_9;
      real_t tmp_103 = tmp_10*tmp_100 + tmp_101*tmp_7;
      real_t tmp_104 = tmp_27*(tmp_22 + tmp_96);
      real_t tmp_105 = tmp_27*(tmp_20 + tmp_98);
      real_t tmp_106 = tmp_104*tmp_33 + tmp_105*tmp_23 - 1.0/3.0;
      real_t tmp_107 = tmp_104*tmp_21 + tmp_105*tmp_30 - 1.0/3.0;
      real_t tmp_108 = tmp_106*tmp_21 + tmp_107*tmp_25;
      real_t tmp_109 = tmp_106*tmp_26 + tmp_107*tmp_23;
      real_t a_0_0 = 0.11846344252809471*tmp_2*(tmp_102*tmp_36 + tmp_103*tmp_38 - tmp_108*tmp_50 - tmp_109*tmp_52 - tmp_53*(tmp_102*tmp_108 + tmp_103*tmp_109)) + 0.11846344252809471*tmp_2*(tmp_19*tmp_36 + tmp_37*tmp_38 - tmp_43*tmp_50 - tmp_51*tmp_52 - tmp_53*(tmp_19*tmp_43 + tmp_37*tmp_51)) + 0.2393143352496831*tmp_2*(tmp_36*tmp_60 + tmp_38*tmp_61 - tmp_50*tmp_66 - tmp_52*tmp_67 - tmp_53*(tmp_60*tmp_66 + tmp_61*tmp_67)) + 0.2844444444444445*tmp_2*(tmp_36*tmp_74 + tmp_38*tmp_75 - tmp_50*tmp_80 - tmp_52*tmp_81 - tmp_53*(tmp_74*tmp_80 + tmp_75*tmp_81)) + 0.2393143352496831*tmp_2*(tmp_36*tmp_88 + tmp_38*tmp_89 - tmp_50*tmp_94 - tmp_52*tmp_95 - tmp_53*(tmp_88*tmp_94 + tmp_89*tmp_95));
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

      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_3;
      real_t tmp_7 = 1.0 / (tmp_1*tmp_4 - tmp_5*tmp_6);
      real_t tmp_8 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9 = p_affine_6_1 + tmp_3;
      real_t tmp_10 = tmp_7*(0.046910077030668018*tmp_8 + tmp_9);
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_7*(0.046910077030668018*tmp_11 + tmp_12);
      real_t tmp_14 = tmp_10*tmp_2 + tmp_13*tmp_4 - 1.0/3.0;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_1*tmp_10 + tmp_13*tmp_15 - 1.0/3.0;
      real_t tmp_17 = tmp_7*(0.23076534494715845*tmp_8 + tmp_9);
      real_t tmp_18 = tmp_7*(0.23076534494715845*tmp_11 + tmp_12);
      real_t tmp_19 = tmp_17*tmp_2 + tmp_18*tmp_4 - 1.0/3.0;
      real_t tmp_20 = tmp_1*tmp_17 + tmp_15*tmp_18 - 1.0/3.0;
      real_t tmp_21 = tmp_7*(0.5*tmp_8 + tmp_9);
      real_t tmp_22 = tmp_7*(0.5*tmp_11 + tmp_12);
      real_t tmp_23 = tmp_2*tmp_21 + tmp_22*tmp_4 - 1.0/3.0;
      real_t tmp_24 = tmp_1*tmp_21 + tmp_15*tmp_22 - 1.0/3.0;
      real_t tmp_25 = tmp_7*(0.7692346550528415*tmp_8 + tmp_9);
      real_t tmp_26 = tmp_7*(0.7692346550528415*tmp_11 + tmp_12);
      real_t tmp_27 = tmp_2*tmp_25 + tmp_26*tmp_4 - 1.0/3.0;
      real_t tmp_28 = tmp_1*tmp_25 + tmp_15*tmp_26 - 1.0/3.0;
      real_t tmp_29 = tmp_7*(0.95308992296933193*tmp_8 + tmp_9);
      real_t tmp_30 = tmp_7*(0.95308992296933193*tmp_11 + tmp_12);
      real_t tmp_31 = tmp_2*tmp_29 + tmp_30*tmp_4 - 1.0/3.0;
      real_t tmp_32 = tmp_1*tmp_29 + tmp_15*tmp_30 - 1.0/3.0;
      real_t a_0_0 = 1.4215613103371365*((tmp_1*tmp_14 + tmp_16*tmp_5)*(tmp_1*tmp_14 + tmp_16*tmp_5)) + 2.8717720229961969*((tmp_1*tmp_19 + tmp_20*tmp_5)*(tmp_1*tmp_19 + tmp_20*tmp_5)) + 3.413333333333334*((tmp_1*tmp_23 + tmp_24*tmp_5)*(tmp_1*tmp_23 + tmp_24*tmp_5)) + 2.8717720229961969*((tmp_1*tmp_27 + tmp_28*tmp_5)*(tmp_1*tmp_27 + tmp_28*tmp_5)) + 1.4215613103371365*((tmp_1*tmp_31 + tmp_32*tmp_5)*(tmp_1*tmp_31 + tmp_32*tmp_5)) + 1.4215613103371365*((tmp_14*tmp_6 + tmp_16*tmp_4)*(tmp_14*tmp_6 + tmp_16*tmp_4)) + 2.8717720229961969*((tmp_19*tmp_6 + tmp_20*tmp_4)*(tmp_19*tmp_6 + tmp_20*tmp_4)) + 3.413333333333334*((tmp_23*tmp_6 + tmp_24*tmp_4)*(tmp_23*tmp_6 + tmp_24*tmp_4)) + 2.8717720229961969*((tmp_27*tmp_6 + tmp_28*tmp_4)*(tmp_27*tmp_6 + tmp_28*tmp_4)) + 1.4215613103371365*((tmp_31*tmp_6 + tmp_32*tmp_4)*(tmp_31*tmp_6 + tmp_32*tmp_4));
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

      // Does nothing.
   }
};

} // eg
} // dg
} // hyteg


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

class EDGEpsilonFormEDGP1_0 : public hyteg::dg::DGForm2D
{
 public:
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
      real_t tmp_11 = tmp_8*(tmp_10 + 0.21132486540518713*tmp_9);
      real_t tmp_12 = tmp_11*tmp_2;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_3;
      real_t tmp_15 = tmp_8*(0.21132486540518713*tmp_13 + tmp_14);
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
      real_t tmp_29 = 0.5*tmp_28;
      real_t tmp_30 = tmp_18*tmp_4 + tmp_23*tmp_7;
      real_t tmp_31 = 1.0*tmp_8;
      real_t tmp_32 = tmp_31*tmp_5;
      real_t tmp_33 = tmp_20*tmp_31;
      real_t tmp_34 = 0.5*p_affine_10_0*(-tmp_32 - tmp_33) + 0.5*p_affine_10_1*tmp_28;
      real_t tmp_35 = -tmp_12 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_36 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_37 = 6/tmp_36;
      real_t tmp_38 = tmp_30*tmp_37;
      real_t tmp_39 = tmp_25*tmp_5;
      real_t tmp_40 = 0.5*p_affine_10_0*(tmp_31*tmp_6 + tmp_33*tmp_7) + 0.5*p_affine_10_1*(tmp_1*tmp_39 + tmp_20*tmp_39 + tmp_26*tmp_7 + tmp_27*tmp_4);
      real_t tmp_41 = 0.5*tmp_36;
      real_t tmp_42 = tmp_8*(tmp_10 + 0.78867513459481287*tmp_9);
      real_t tmp_43 = tmp_2*tmp_42;
      real_t tmp_44 = tmp_8*(0.78867513459481287*tmp_13 + tmp_14);
      real_t tmp_45 = tmp_44*tmp_5;
      real_t tmp_46 = tmp_43 + tmp_45;
      real_t tmp_47 = tmp_46 - 1.0/3.0;
      real_t tmp_48 = tmp_4*tmp_42;
      real_t tmp_49 = tmp_20*tmp_44;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_50 - 1.0/3.0;
      real_t tmp_52 = p_affine_10_0*(tmp_1*tmp_47 + tmp_5*tmp_51);
      real_t tmp_53 = tmp_4*tmp_47 + tmp_51*tmp_7;
      real_t tmp_54 = -tmp_43 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_55 = tmp_37*tmp_53;
      real_t tmp_56 = 0.5*tmp_36;
      real_t tmp_57 = 0.25*tmp_8;
      real_t tmp_58 = tmp_2*tmp_57;
      real_t tmp_59 = 0.5*p_affine_10_0*tmp_32 + 0.5*p_affine_10_1*tmp_27;
      real_t tmp_60 = tmp_4*tmp_57;
      real_t tmp_61 = 0.5*p_affine_10_0*tmp_33 + 0.5*p_affine_10_1*tmp_26;
      real_t a_0_0 = tmp_41*(-tmp_24*tmp_29 - tmp_30*tmp_34 + tmp_35*tmp_38 - tmp_35*tmp_40) + tmp_56*(-tmp_29*tmp_52 - tmp_34*tmp_53 - tmp_40*tmp_54 + tmp_54*tmp_55);
      real_t a_1_0 = tmp_41*(tmp_17*tmp_38 - tmp_17*tmp_40 - tmp_24*tmp_58 - tmp_30*tmp_59) + tmp_56*(-tmp_40*tmp_46 + tmp_46*tmp_55 - tmp_52*tmp_58 - tmp_53*tmp_59);
      real_t a_2_0 = tmp_41*(tmp_22*tmp_38 - tmp_22*tmp_40 - tmp_24*tmp_60 - tmp_30*tmp_61) + tmp_56*(-tmp_40*tmp_50 + tmp_50*tmp_55 - tmp_52*tmp_60 - tmp_53*tmp_61);
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
      real_t tmp_10 = p_affine_6_1 + 0.21132486540518713*tmp_9;
      real_t tmp_11 = tmp_8*(tmp_0 + tmp_10);
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + 0.21132486540518713*tmp_12;
      real_t tmp_14 = tmp_8*(tmp_13 + tmp_3);
      real_t tmp_15 = tmp_11*tmp_2 + tmp_14*tmp_5 - 1.0/3.0;
      real_t tmp_16 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_17 = tmp_11*tmp_4 + tmp_14*tmp_16 - 1.0/3.0;
      real_t tmp_18 = p_affine_10_0*(tmp_1*tmp_15 + tmp_17*tmp_5);
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
      real_t tmp_29 = 0.5*tmp_28;
      real_t tmp_30 = tmp_15*tmp_4 + tmp_17*tmp_7;
      real_t tmp_31 = 1.0*tmp_23;
      real_t tmp_32 = tmp_22*tmp_31;
      real_t tmp_33 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_34 = tmp_31*tmp_33;
      real_t tmp_35 = 0.5*p_affine_10_0*(-tmp_32 - tmp_34) + 0.5*p_affine_10_1*tmp_28;
      real_t tmp_36 = tmp_23*(tmp_10 + tmp_21);
      real_t tmp_37 = tmp_20*tmp_36;
      real_t tmp_38 = tmp_26*tmp_36;
      real_t tmp_39 = tmp_23*(tmp_13 + tmp_19);
      real_t tmp_40 = tmp_22*tmp_39;
      real_t tmp_41 = tmp_33*tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_44 = 6/tmp_43;
      real_t tmp_45 = tmp_30*tmp_44;
      real_t tmp_46 = 1.0*tmp_8;
      real_t tmp_47 = 0.5*tmp_8;
      real_t tmp_48 = tmp_4*tmp_47;
      real_t tmp_49 = tmp_47*tmp_5;
      real_t tmp_50 = 0.5*p_affine_10_0*(tmp_16*tmp_46*tmp_7 + tmp_46*tmp_6) + 0.5*p_affine_10_1*(tmp_1*tmp_49 + tmp_16*tmp_49 + tmp_2*tmp_48 + tmp_48*tmp_7);
      real_t tmp_51 = 0.5*tmp_43;
      real_t tmp_52 = p_affine_6_1 + 0.78867513459481287*tmp_9;
      real_t tmp_53 = tmp_8*(tmp_0 + tmp_52);
      real_t tmp_54 = p_affine_6_0 + 0.78867513459481287*tmp_12;
      real_t tmp_55 = tmp_8*(tmp_3 + tmp_54);
      real_t tmp_56 = tmp_2*tmp_53 + tmp_5*tmp_55 - 1.0/3.0;
      real_t tmp_57 = tmp_16*tmp_55 + tmp_4*tmp_53 - 1.0/3.0;
      real_t tmp_58 = p_affine_10_0*(tmp_1*tmp_56 + tmp_5*tmp_57);
      real_t tmp_59 = tmp_4*tmp_56 + tmp_57*tmp_7;
      real_t tmp_60 = tmp_23*(tmp_21 + tmp_52);
      real_t tmp_61 = tmp_20*tmp_60;
      real_t tmp_62 = tmp_26*tmp_60;
      real_t tmp_63 = tmp_23*(tmp_19 + tmp_54);
      real_t tmp_64 = tmp_22*tmp_63;
      real_t tmp_65 = tmp_33*tmp_63;
      real_t tmp_66 = -tmp_61 - tmp_62 - tmp_64 - tmp_65 + 1;
      real_t tmp_67 = tmp_44*tmp_59;
      real_t tmp_68 = 0.5*tmp_43;
      real_t tmp_69 = 0.25*tmp_23;
      real_t tmp_70 = tmp_26*tmp_69;
      real_t tmp_71 = 0.5*p_affine_10_0*tmp_32 + 0.5*p_affine_10_1*tmp_27;
      real_t tmp_72 = tmp_38 + tmp_40;
      real_t tmp_73 = tmp_62 + tmp_64;
      real_t tmp_74 = tmp_20*tmp_69;
      real_t tmp_75 = 0.5*p_affine_10_0*tmp_34 + 0.5*p_affine_10_1*tmp_25;
      real_t tmp_76 = tmp_37 + tmp_41;
      real_t tmp_77 = tmp_61 + tmp_65;
      real_t a_0_0 = tmp_51*(tmp_18*tmp_29 + tmp_30*tmp_35 - tmp_42*tmp_45 - tmp_42*tmp_50) + tmp_68*(tmp_29*tmp_58 + tmp_35*tmp_59 - tmp_50*tmp_66 - tmp_66*tmp_67);
      real_t a_1_0 = tmp_51*(tmp_18*tmp_70 + tmp_30*tmp_71 - tmp_45*tmp_72 - tmp_50*tmp_72) + tmp_68*(-tmp_50*tmp_73 + tmp_58*tmp_70 + tmp_59*tmp_71 - tmp_67*tmp_73);
      real_t a_2_0 = tmp_51*(tmp_18*tmp_74 + tmp_30*tmp_75 - tmp_45*tmp_76 - tmp_50*tmp_76) + tmp_68*(-tmp_50*tmp_77 + tmp_58*tmp_74 + tmp_59*tmp_75 - tmp_67*tmp_77);
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
      real_t tmp_11 = tmp_8*(tmp_10 + 0.21132486540518713*tmp_9);
      real_t tmp_12 = tmp_11*tmp_2;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_3;
      real_t tmp_15 = tmp_8*(0.21132486540518713*tmp_13 + tmp_14);
      real_t tmp_16 = tmp_15*tmp_5;
      real_t tmp_17 = tmp_12 + tmp_16;
      real_t tmp_18 = tmp_17 - 1.0/3.0;
      real_t tmp_19 = tmp_11*tmp_4;
      real_t tmp_20 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_21 = tmp_15*tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_22 - 1.0/3.0;
      real_t tmp_24 = tmp_1*tmp_18 + tmp_23*tmp_5;
      real_t tmp_25 = 0.5*tmp_8;
      real_t tmp_26 = tmp_25*tmp_4;
      real_t tmp_27 = tmp_2*tmp_25;
      real_t tmp_28 = -tmp_26 - tmp_27;
      real_t tmp_29 = p_affine_10_0*tmp_28;
      real_t tmp_30 = 1.0*tmp_8;
      real_t tmp_31 = tmp_30*tmp_5;
      real_t tmp_32 = tmp_20*tmp_30;
      real_t tmp_33 = p_affine_10_0*(-tmp_31 - tmp_32) + p_affine_10_1*tmp_28;
      real_t tmp_34 = tmp_18*tmp_4 + tmp_23*tmp_7;
      real_t tmp_35 = -tmp_12 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_36 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_37 = 24/tmp_36;
      real_t tmp_38 = tmp_34*tmp_37;
      real_t tmp_39 = tmp_25*tmp_5;
      real_t tmp_40 = p_affine_10_0*(tmp_30*tmp_6 + tmp_32*tmp_7) + p_affine_10_1*(tmp_1*tmp_39 + tmp_20*tmp_39 + tmp_26*tmp_7 + tmp_27*tmp_4);
      real_t tmp_41 = 0.5*tmp_36;
      real_t tmp_42 = tmp_8*(tmp_10 + 0.78867513459481287*tmp_9);
      real_t tmp_43 = tmp_2*tmp_42;
      real_t tmp_44 = tmp_8*(0.78867513459481287*tmp_13 + tmp_14);
      real_t tmp_45 = tmp_44*tmp_5;
      real_t tmp_46 = tmp_43 + tmp_45;
      real_t tmp_47 = tmp_46 - 1.0/3.0;
      real_t tmp_48 = tmp_4*tmp_42;
      real_t tmp_49 = tmp_20*tmp_44;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_50 - 1.0/3.0;
      real_t tmp_52 = tmp_1*tmp_47 + tmp_5*tmp_51;
      real_t tmp_53 = tmp_4*tmp_47 + tmp_51*tmp_7;
      real_t tmp_54 = -tmp_43 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_55 = tmp_37*tmp_53;
      real_t tmp_56 = 0.5*tmp_36;
      real_t tmp_57 = p_affine_10_0*tmp_27;
      real_t tmp_58 = p_affine_10_0*tmp_31 + p_affine_10_1*tmp_27;
      real_t tmp_59 = p_affine_10_0*tmp_26;
      real_t tmp_60 = p_affine_10_0*tmp_32 + p_affine_10_1*tmp_26;
      real_t a_0_0 = tmp_41*(-tmp_24*tmp_29 - tmp_33*tmp_34 + tmp_35*tmp_38 - tmp_35*tmp_40) + tmp_56*(-tmp_29*tmp_52 - tmp_33*tmp_53 - tmp_40*tmp_54 + tmp_54*tmp_55);
      real_t a_1_0 = tmp_41*(tmp_17*tmp_38 - tmp_17*tmp_40 - tmp_24*tmp_57 - tmp_34*tmp_58) + tmp_56*(-tmp_40*tmp_46 + tmp_46*tmp_55 - tmp_52*tmp_57 - tmp_53*tmp_58);
      real_t a_2_0 = tmp_41*(tmp_22*tmp_38 - tmp_22*tmp_40 - tmp_24*tmp_59 - tmp_34*tmp_60) + tmp_56*(-tmp_40*tmp_50 + tmp_50*tmp_55 - tmp_52*tmp_59 - tmp_53*tmp_60);
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




class EDGEpsilonFormP1EDG_0 : public hyteg::dg::DGForm2D
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
      real_t tmp_11 = tmp_8*(tmp_10 + 0.21132486540518713*tmp_9);
      real_t tmp_12 = tmp_11*tmp_2;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_3;
      real_t tmp_15 = tmp_8*(0.21132486540518713*tmp_13 + tmp_14);
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
      real_t tmp_29 = 0.5*tmp_28;
      real_t tmp_30 = tmp_18*tmp_4 + tmp_23*tmp_7;
      real_t tmp_31 = 1.0*tmp_8;
      real_t tmp_32 = tmp_31*tmp_5;
      real_t tmp_33 = tmp_20*tmp_31;
      real_t tmp_34 = 0.5*p_affine_10_0*(-tmp_32 - tmp_33) + 0.5*p_affine_10_1*tmp_28;
      real_t tmp_35 = -tmp_12 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_36 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_37 = 6/tmp_36;
      real_t tmp_38 = tmp_30*tmp_37;
      real_t tmp_39 = tmp_25*tmp_5;
      real_t tmp_40 = 0.5*p_affine_10_0*(tmp_31*tmp_6 + tmp_33*tmp_7) + 0.5*p_affine_10_1*(tmp_1*tmp_39 + tmp_20*tmp_39 + tmp_26*tmp_7 + tmp_27*tmp_4);
      real_t tmp_41 = 0.5*tmp_36;
      real_t tmp_42 = tmp_8*(tmp_10 + 0.78867513459481287*tmp_9);
      real_t tmp_43 = tmp_2*tmp_42;
      real_t tmp_44 = tmp_8*(0.78867513459481287*tmp_13 + tmp_14);
      real_t tmp_45 = tmp_44*tmp_5;
      real_t tmp_46 = tmp_43 + tmp_45;
      real_t tmp_47 = tmp_46 - 1.0/3.0;
      real_t tmp_48 = tmp_4*tmp_42;
      real_t tmp_49 = tmp_20*tmp_44;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_50 - 1.0/3.0;
      real_t tmp_52 = p_affine_10_0*(tmp_1*tmp_47 + tmp_5*tmp_51);
      real_t tmp_53 = tmp_4*tmp_47 + tmp_51*tmp_7;
      real_t tmp_54 = -tmp_43 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_55 = tmp_37*tmp_53;
      real_t tmp_56 = 0.5*tmp_36;
      real_t tmp_57 = 0.25*tmp_8;
      real_t tmp_58 = tmp_2*tmp_57;
      real_t tmp_59 = 0.5*p_affine_10_0*tmp_32 + 0.5*p_affine_10_1*tmp_27;
      real_t tmp_60 = tmp_4*tmp_57;
      real_t tmp_61 = 0.5*p_affine_10_0*tmp_33 + 0.5*p_affine_10_1*tmp_26;
      real_t a_0_0 = tmp_41*(-tmp_24*tmp_29 - tmp_30*tmp_34 + tmp_35*tmp_38 - tmp_35*tmp_40) + tmp_56*(-tmp_29*tmp_52 - tmp_34*tmp_53 - tmp_40*tmp_54 + tmp_54*tmp_55);
      real_t a_0_1 = tmp_41*(tmp_17*tmp_38 - tmp_17*tmp_40 - tmp_24*tmp_58 - tmp_30*tmp_59) + tmp_56*(-tmp_40*tmp_46 + tmp_46*tmp_55 - tmp_52*tmp_58 - tmp_53*tmp_59);
      real_t a_0_2 = tmp_41*(tmp_22*tmp_38 - tmp_22*tmp_40 - tmp_24*tmp_60 - tmp_30*tmp_61) + tmp_56*(-tmp_40*tmp_50 + tmp_50*tmp_55 - tmp_52*tmp_60 - tmp_53*tmp_61);
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
      real_t tmp_10 = p_affine_6_1 + 0.21132486540518713*tmp_9;
      real_t tmp_11 = tmp_8*(tmp_0 + tmp_10);
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + 0.21132486540518713*tmp_12;
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
      real_t tmp_29 = 0.5*tmp_28;
      real_t tmp_30 = tmp_15*tmp_4 + tmp_17*tmp_7;
      real_t tmp_31 = 1.0*tmp_23;
      real_t tmp_32 = tmp_22*tmp_31;
      real_t tmp_33 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_34 = tmp_31*tmp_33;
      real_t tmp_35 = 0.5*p_affine_10_0*(-tmp_32 - tmp_34) + 0.5*p_affine_10_1*tmp_28;
      real_t tmp_36 = tmp_23*(tmp_10 + tmp_21);
      real_t tmp_37 = tmp_20*tmp_36;
      real_t tmp_38 = tmp_26*tmp_36;
      real_t tmp_39 = tmp_23*(tmp_13 + tmp_19);
      real_t tmp_40 = tmp_22*tmp_39;
      real_t tmp_41 = tmp_33*tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_44 = 6/tmp_43;
      real_t tmp_45 = tmp_30*tmp_44;
      real_t tmp_46 = 1.0*tmp_8;
      real_t tmp_47 = 0.5*tmp_8;
      real_t tmp_48 = tmp_4*tmp_47;
      real_t tmp_49 = tmp_47*tmp_5;
      real_t tmp_50 = 0.5*p_affine_10_0*(tmp_16*tmp_46*tmp_7 + tmp_46*tmp_6) + 0.5*p_affine_10_1*(tmp_1*tmp_49 + tmp_16*tmp_49 + tmp_2*tmp_48 + tmp_48*tmp_7);
      real_t tmp_51 = 0.5*tmp_43;
      real_t tmp_52 = p_affine_6_1 + 0.78867513459481287*tmp_9;
      real_t tmp_53 = tmp_8*(tmp_0 + tmp_52);
      real_t tmp_54 = p_affine_6_0 + 0.78867513459481287*tmp_12;
      real_t tmp_55 = tmp_8*(tmp_3 + tmp_54);
      real_t tmp_56 = tmp_2*tmp_53 + tmp_5*tmp_55 - 1.0/3.0;
      real_t tmp_57 = tmp_16*tmp_55 + tmp_4*tmp_53 - 1.0/3.0;
      real_t tmp_58 = p_affine_10_0*(tmp_1*tmp_56 + tmp_5*tmp_57);
      real_t tmp_59 = tmp_4*tmp_56 + tmp_57*tmp_7;
      real_t tmp_60 = tmp_23*(tmp_21 + tmp_52);
      real_t tmp_61 = tmp_20*tmp_60;
      real_t tmp_62 = tmp_26*tmp_60;
      real_t tmp_63 = tmp_23*(tmp_19 + tmp_54);
      real_t tmp_64 = tmp_22*tmp_63;
      real_t tmp_65 = tmp_33*tmp_63;
      real_t tmp_66 = -tmp_61 - tmp_62 - tmp_64 - tmp_65 + 1;
      real_t tmp_67 = tmp_44*tmp_59;
      real_t tmp_68 = 0.5*tmp_43;
      real_t tmp_69 = 0.25*tmp_23;
      real_t tmp_70 = tmp_26*tmp_69;
      real_t tmp_71 = 0.5*p_affine_10_0*tmp_32 + 0.5*p_affine_10_1*tmp_27;
      real_t tmp_72 = tmp_38 + tmp_40;
      real_t tmp_73 = tmp_62 + tmp_64;
      real_t tmp_74 = tmp_20*tmp_69;
      real_t tmp_75 = 0.5*p_affine_10_0*tmp_34 + 0.5*p_affine_10_1*tmp_25;
      real_t tmp_76 = tmp_37 + tmp_41;
      real_t tmp_77 = tmp_61 + tmp_65;
      real_t a_0_0 = tmp_51*(-tmp_18*tmp_29 - tmp_30*tmp_35 - tmp_42*tmp_45 + tmp_42*tmp_50) + tmp_68*(-tmp_29*tmp_58 - tmp_35*tmp_59 + tmp_50*tmp_66 - tmp_66*tmp_67);
      real_t a_0_1 = tmp_51*(-tmp_18*tmp_70 - tmp_30*tmp_71 - tmp_45*tmp_72 + tmp_50*tmp_72) + tmp_68*(tmp_50*tmp_73 - tmp_58*tmp_70 - tmp_59*tmp_71 - tmp_67*tmp_73);
      real_t a_0_2 = tmp_51*(-tmp_18*tmp_74 - tmp_30*tmp_75 - tmp_45*tmp_76 + tmp_50*tmp_76) + tmp_68*(tmp_50*tmp_77 - tmp_58*tmp_74 - tmp_59*tmp_75 - tmp_67*tmp_77);
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
      real_t tmp_6 = tmp_4*tmp_5;
      real_t tmp_7 = p_affine_2_0 + tmp_3;
      real_t tmp_8 = 1.0 / (-tmp_1*tmp_7 + tmp_6);
      real_t tmp_9 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + tmp_0;
      real_t tmp_11 = tmp_8*(tmp_10 + 0.21132486540518713*tmp_9);
      real_t tmp_12 = tmp_11*tmp_2;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_3;
      real_t tmp_15 = tmp_8*(0.21132486540518713*tmp_13 + tmp_14);
      real_t tmp_16 = tmp_15*tmp_5;
      real_t tmp_17 = tmp_12 + tmp_16;
      real_t tmp_18 = tmp_17 - 1.0/3.0;
      real_t tmp_19 = tmp_11*tmp_4;
      real_t tmp_20 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_21 = tmp_15*tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_22 - 1.0/3.0;
      real_t tmp_24 = tmp_1*tmp_18 + tmp_23*tmp_5;
      real_t tmp_25 = 0.5*tmp_8;
      real_t tmp_26 = tmp_25*tmp_4;
      real_t tmp_27 = tmp_2*tmp_25;
      real_t tmp_28 = -tmp_26 - tmp_27;
      real_t tmp_29 = p_affine_10_0*tmp_28;
      real_t tmp_30 = 1.0*tmp_8;
      real_t tmp_31 = tmp_30*tmp_5;
      real_t tmp_32 = tmp_20*tmp_30;
      real_t tmp_33 = p_affine_10_0*(-tmp_31 - tmp_32) + p_affine_10_1*tmp_28;
      real_t tmp_34 = tmp_18*tmp_4 + tmp_23*tmp_7;
      real_t tmp_35 = -tmp_12 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_36 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_37 = 24/tmp_36;
      real_t tmp_38 = tmp_34*tmp_37;
      real_t tmp_39 = tmp_25*tmp_5;
      real_t tmp_40 = p_affine_10_0*(tmp_30*tmp_6 + tmp_32*tmp_7) + p_affine_10_1*(tmp_1*tmp_39 + tmp_20*tmp_39 + tmp_26*tmp_7 + tmp_27*tmp_4);
      real_t tmp_41 = 0.5*tmp_36;
      real_t tmp_42 = tmp_8*(tmp_10 + 0.78867513459481287*tmp_9);
      real_t tmp_43 = tmp_2*tmp_42;
      real_t tmp_44 = tmp_8*(0.78867513459481287*tmp_13 + tmp_14);
      real_t tmp_45 = tmp_44*tmp_5;
      real_t tmp_46 = tmp_43 + tmp_45;
      real_t tmp_47 = tmp_46 - 1.0/3.0;
      real_t tmp_48 = tmp_4*tmp_42;
      real_t tmp_49 = tmp_20*tmp_44;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_50 - 1.0/3.0;
      real_t tmp_52 = tmp_1*tmp_47 + tmp_5*tmp_51;
      real_t tmp_53 = tmp_4*tmp_47 + tmp_51*tmp_7;
      real_t tmp_54 = -tmp_43 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_55 = tmp_37*tmp_53;
      real_t tmp_56 = 0.5*tmp_36;
      real_t tmp_57 = p_affine_10_0*tmp_27;
      real_t tmp_58 = p_affine_10_0*tmp_31 + p_affine_10_1*tmp_27;
      real_t tmp_59 = p_affine_10_0*tmp_26;
      real_t tmp_60 = p_affine_10_0*tmp_32 + p_affine_10_1*tmp_26;
      real_t a_0_0 = tmp_41*(-tmp_24*tmp_29 - tmp_33*tmp_34 + tmp_35*tmp_38 - tmp_35*tmp_40) + tmp_56*(-tmp_29*tmp_52 - tmp_33*tmp_53 - tmp_40*tmp_54 + tmp_54*tmp_55);
      real_t a_0_1 = tmp_41*(tmp_17*tmp_38 - tmp_17*tmp_40 - tmp_24*tmp_57 - tmp_34*tmp_58) + tmp_56*(-tmp_40*tmp_46 + tmp_46*tmp_55 - tmp_52*tmp_57 - tmp_53*tmp_58);
      real_t a_0_2 = tmp_41*(tmp_22*tmp_38 - tmp_22*tmp_40 - tmp_24*tmp_59 - tmp_34*tmp_60) + tmp_56*(-tmp_40*tmp_50 + tmp_50*tmp_55 - tmp_52*tmp_59 - tmp_53*tmp_60);
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




class EDGEpsilonFormEDGP1_1 : public hyteg::dg::DGForm2D
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
      real_t tmp_11 = tmp_8*(tmp_10 + 0.21132486540518713*tmp_9);
      real_t tmp_12 = tmp_11*tmp_2;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_0;
      real_t tmp_15 = tmp_8*(0.21132486540518713*tmp_13 + tmp_14);
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
      real_t tmp_29 = 0.5*tmp_28;
      real_t tmp_30 = tmp_18*tmp_7 + tmp_23*tmp_4;
      real_t tmp_31 = 1.0*tmp_8;
      real_t tmp_32 = tmp_1*tmp_31;
      real_t tmp_33 = tmp_2*tmp_31;
      real_t tmp_34 = 0.5*p_affine_10_0*tmp_28 + 0.5*p_affine_10_1*(-tmp_32 - tmp_33);
      real_t tmp_35 = -tmp_12 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_36 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_37 = 6/tmp_36;
      real_t tmp_38 = tmp_30*tmp_37;
      real_t tmp_39 = tmp_1*tmp_25;
      real_t tmp_40 = 0.5*p_affine_10_0*(tmp_2*tmp_39 + tmp_26*tmp_7 + tmp_27*tmp_4 + tmp_39*tmp_6) + 0.5*p_affine_10_1*(tmp_31*tmp_5 + tmp_33*tmp_7);
      real_t tmp_41 = 0.5*tmp_36;
      real_t tmp_42 = tmp_8*(tmp_10 + 0.78867513459481287*tmp_9);
      real_t tmp_43 = tmp_2*tmp_42;
      real_t tmp_44 = tmp_8*(0.78867513459481287*tmp_13 + tmp_14);
      real_t tmp_45 = tmp_4*tmp_44;
      real_t tmp_46 = tmp_43 + tmp_45;
      real_t tmp_47 = tmp_46 - 1.0/3.0;
      real_t tmp_48 = tmp_1*tmp_42;
      real_t tmp_49 = tmp_20*tmp_44;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_50 - 1.0/3.0;
      real_t tmp_52 = p_affine_10_1*(tmp_1*tmp_47 + tmp_51*tmp_6);
      real_t tmp_53 = tmp_4*tmp_51 + tmp_47*tmp_7;
      real_t tmp_54 = -tmp_43 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_55 = tmp_37*tmp_53;
      real_t tmp_56 = 0.5*tmp_36;
      real_t tmp_57 = 0.25*tmp_8;
      real_t tmp_58 = tmp_4*tmp_57;
      real_t tmp_59 = 0.5*p_affine_10_0*tmp_26 + 0.5*p_affine_10_1*tmp_33;
      real_t tmp_60 = tmp_20*tmp_57;
      real_t tmp_61 = 0.5*p_affine_10_0*tmp_27 + 0.5*p_affine_10_1*tmp_32;
      real_t a_0_0 = tmp_41*(-tmp_24*tmp_29 - tmp_30*tmp_34 + tmp_35*tmp_38 - tmp_35*tmp_40) + tmp_56*(-tmp_29*tmp_52 - tmp_34*tmp_53 - tmp_40*tmp_54 + tmp_54*tmp_55);
      real_t a_1_0 = tmp_41*(tmp_17*tmp_38 - tmp_17*tmp_40 - tmp_24*tmp_58 - tmp_30*tmp_59) + tmp_56*(-tmp_40*tmp_46 + tmp_46*tmp_55 - tmp_52*tmp_58 - tmp_53*tmp_59);
      real_t a_2_0 = tmp_41*(tmp_22*tmp_38 - tmp_22*tmp_40 - tmp_24*tmp_60 - tmp_30*tmp_61) + tmp_56*(-tmp_40*tmp_50 + tmp_50*tmp_55 - tmp_52*tmp_60 - tmp_53*tmp_61);
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
      real_t tmp_10 = p_affine_6_1 + 0.21132486540518713*tmp_9;
      real_t tmp_11 = tmp_8*(tmp_10 + tmp_3);
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + 0.21132486540518713*tmp_12;
      real_t tmp_14 = tmp_8*(tmp_0 + tmp_13);
      real_t tmp_15 = tmp_11*tmp_2 + tmp_14*tmp_4 - 1.0/3.0;
      real_t tmp_16 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_17 = tmp_1*tmp_11 + tmp_14*tmp_16 - 1.0/3.0;
      real_t tmp_18 = p_affine_10_1*(tmp_1*tmp_15 + tmp_17*tmp_6);
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
      real_t tmp_29 = 0.5*tmp_28;
      real_t tmp_30 = tmp_15*tmp_7 + tmp_17*tmp_4;
      real_t tmp_31 = 1.0*tmp_23;
      real_t tmp_32 = tmp_22*tmp_31;
      real_t tmp_33 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_34 = tmp_31*tmp_33;
      real_t tmp_35 = 0.5*p_affine_10_0*tmp_28 + 0.5*p_affine_10_1*(-tmp_32 - tmp_34);
      real_t tmp_36 = tmp_23*(tmp_10 + tmp_19);
      real_t tmp_37 = tmp_22*tmp_36;
      real_t tmp_38 = tmp_33*tmp_36;
      real_t tmp_39 = tmp_23*(tmp_13 + tmp_21);
      real_t tmp_40 = tmp_20*tmp_39;
      real_t tmp_41 = tmp_26*tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_44 = 6/tmp_43;
      real_t tmp_45 = tmp_30*tmp_44;
      real_t tmp_46 = 1.0*tmp_8;
      real_t tmp_47 = 0.5*tmp_8;
      real_t tmp_48 = tmp_1*tmp_47;
      real_t tmp_49 = tmp_4*tmp_47;
      real_t tmp_50 = 0.5*p_affine_10_0*(tmp_16*tmp_49 + tmp_2*tmp_48 + tmp_48*tmp_6 + tmp_49*tmp_7) + 0.5*p_affine_10_1*(tmp_2*tmp_46*tmp_7 + tmp_46*tmp_5);
      real_t tmp_51 = 0.5*tmp_43;
      real_t tmp_52 = p_affine_6_1 + 0.78867513459481287*tmp_9;
      real_t tmp_53 = tmp_8*(tmp_3 + tmp_52);
      real_t tmp_54 = p_affine_6_0 + 0.78867513459481287*tmp_12;
      real_t tmp_55 = tmp_8*(tmp_0 + tmp_54);
      real_t tmp_56 = tmp_2*tmp_53 + tmp_4*tmp_55 - 1.0/3.0;
      real_t tmp_57 = tmp_1*tmp_53 + tmp_16*tmp_55 - 1.0/3.0;
      real_t tmp_58 = p_affine_10_1*(tmp_1*tmp_56 + tmp_57*tmp_6);
      real_t tmp_59 = tmp_4*tmp_57 + tmp_56*tmp_7;
      real_t tmp_60 = tmp_23*(tmp_19 + tmp_52);
      real_t tmp_61 = tmp_22*tmp_60;
      real_t tmp_62 = tmp_33*tmp_60;
      real_t tmp_63 = tmp_23*(tmp_21 + tmp_54);
      real_t tmp_64 = tmp_20*tmp_63;
      real_t tmp_65 = tmp_26*tmp_63;
      real_t tmp_66 = -tmp_61 - tmp_62 - tmp_64 - tmp_65 + 1;
      real_t tmp_67 = tmp_44*tmp_59;
      real_t tmp_68 = 0.5*tmp_43;
      real_t tmp_69 = 0.25*tmp_23;
      real_t tmp_70 = tmp_20*tmp_69;
      real_t tmp_71 = 0.5*p_affine_10_0*tmp_25 + 0.5*p_affine_10_1*tmp_34;
      real_t tmp_72 = tmp_38 + tmp_40;
      real_t tmp_73 = tmp_62 + tmp_64;
      real_t tmp_74 = tmp_26*tmp_69;
      real_t tmp_75 = 0.5*p_affine_10_0*tmp_27 + 0.5*p_affine_10_1*tmp_32;
      real_t tmp_76 = tmp_37 + tmp_41;
      real_t tmp_77 = tmp_61 + tmp_65;
      real_t a_0_0 = tmp_51*(tmp_18*tmp_29 + tmp_30*tmp_35 - tmp_42*tmp_45 - tmp_42*tmp_50) + tmp_68*(tmp_29*tmp_58 + tmp_35*tmp_59 - tmp_50*tmp_66 - tmp_66*tmp_67);
      real_t a_1_0 = tmp_51*(tmp_18*tmp_70 + tmp_30*tmp_71 - tmp_45*tmp_72 - tmp_50*tmp_72) + tmp_68*(-tmp_50*tmp_73 + tmp_58*tmp_70 + tmp_59*tmp_71 - tmp_67*tmp_73);
      real_t a_2_0 = tmp_51*(tmp_18*tmp_74 + tmp_30*tmp_75 - tmp_45*tmp_76 - tmp_50*tmp_76) + tmp_68*(-tmp_50*tmp_77 + tmp_58*tmp_74 + tmp_59*tmp_75 - tmp_67*tmp_77);
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
      real_t tmp_2 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = tmp_1*tmp_4;
      real_t tmp_6 = p_affine_2_0 + tmp_0;
      real_t tmp_7 = p_affine_1_1 + tmp_3;
      real_t tmp_8 = 1.0 / (tmp_5 - tmp_6*tmp_7);
      real_t tmp_9 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + tmp_3;
      real_t tmp_11 = tmp_8*(tmp_10 + 0.21132486540518713*tmp_9);
      real_t tmp_12 = tmp_11*tmp_2;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_0;
      real_t tmp_15 = tmp_8*(0.21132486540518713*tmp_13 + tmp_14);
      real_t tmp_16 = tmp_15*tmp_4;
      real_t tmp_17 = tmp_12 + tmp_16;
      real_t tmp_18 = tmp_17 - 1.0/3.0;
      real_t tmp_19 = tmp_1*tmp_11;
      real_t tmp_20 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_21 = tmp_15*tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_22 - 1.0/3.0;
      real_t tmp_24 = tmp_1*tmp_18 + tmp_23*tmp_6;
      real_t tmp_25 = 0.5*tmp_8;
      real_t tmp_26 = tmp_25*tmp_4;
      real_t tmp_27 = tmp_20*tmp_25;
      real_t tmp_28 = -tmp_26 - tmp_27;
      real_t tmp_29 = p_affine_10_1*tmp_28;
      real_t tmp_30 = 1.0*tmp_8;
      real_t tmp_31 = tmp_1*tmp_30;
      real_t tmp_32 = tmp_2*tmp_30;
      real_t tmp_33 = p_affine_10_0*tmp_28 + p_affine_10_1*(-tmp_31 - tmp_32);
      real_t tmp_34 = tmp_18*tmp_7 + tmp_23*tmp_4;
      real_t tmp_35 = -tmp_12 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_36 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_37 = 24/tmp_36;
      real_t tmp_38 = tmp_34*tmp_37;
      real_t tmp_39 = tmp_1*tmp_25;
      real_t tmp_40 = p_affine_10_0*(tmp_2*tmp_39 + tmp_26*tmp_7 + tmp_27*tmp_4 + tmp_39*tmp_6) + p_affine_10_1*(tmp_30*tmp_5 + tmp_32*tmp_7);
      real_t tmp_41 = 0.5*tmp_36;
      real_t tmp_42 = tmp_8*(tmp_10 + 0.78867513459481287*tmp_9);
      real_t tmp_43 = tmp_2*tmp_42;
      real_t tmp_44 = tmp_8*(0.78867513459481287*tmp_13 + tmp_14);
      real_t tmp_45 = tmp_4*tmp_44;
      real_t tmp_46 = tmp_43 + tmp_45;
      real_t tmp_47 = tmp_46 - 1.0/3.0;
      real_t tmp_48 = tmp_1*tmp_42;
      real_t tmp_49 = tmp_20*tmp_44;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_50 - 1.0/3.0;
      real_t tmp_52 = tmp_1*tmp_47 + tmp_51*tmp_6;
      real_t tmp_53 = tmp_4*tmp_51 + tmp_47*tmp_7;
      real_t tmp_54 = -tmp_43 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_55 = tmp_37*tmp_53;
      real_t tmp_56 = 0.5*tmp_36;
      real_t tmp_57 = p_affine_10_1*tmp_26;
      real_t tmp_58 = p_affine_10_0*tmp_26 + p_affine_10_1*tmp_32;
      real_t tmp_59 = p_affine_10_1*tmp_27;
      real_t tmp_60 = p_affine_10_0*tmp_27 + p_affine_10_1*tmp_31;
      real_t a_0_0 = tmp_41*(-tmp_24*tmp_29 - tmp_33*tmp_34 + tmp_35*tmp_38 - tmp_35*tmp_40) + tmp_56*(-tmp_29*tmp_52 - tmp_33*tmp_53 - tmp_40*tmp_54 + tmp_54*tmp_55);
      real_t a_1_0 = tmp_41*(tmp_17*tmp_38 - tmp_17*tmp_40 - tmp_24*tmp_57 - tmp_34*tmp_58) + tmp_56*(-tmp_40*tmp_46 + tmp_46*tmp_55 - tmp_52*tmp_57 - tmp_53*tmp_58);
      real_t a_2_0 = tmp_41*(tmp_22*tmp_38 - tmp_22*tmp_40 - tmp_24*tmp_59 - tmp_34*tmp_60) + tmp_56*(-tmp_40*tmp_50 + tmp_50*tmp_55 - tmp_52*tmp_59 - tmp_53*tmp_60);
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




class EDGEpsilonFormP1EDG_1 : public hyteg::dg::DGForm2D
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
      real_t tmp_11 = tmp_8*(tmp_10 + 0.21132486540518713*tmp_9);
      real_t tmp_12 = tmp_11*tmp_2;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_0;
      real_t tmp_15 = tmp_8*(0.21132486540518713*tmp_13 + tmp_14);
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
      real_t tmp_29 = 0.5*tmp_28;
      real_t tmp_30 = tmp_18*tmp_7 + tmp_23*tmp_4;
      real_t tmp_31 = 1.0*tmp_8;
      real_t tmp_32 = tmp_1*tmp_31;
      real_t tmp_33 = tmp_2*tmp_31;
      real_t tmp_34 = 0.5*p_affine_10_0*tmp_28 + 0.5*p_affine_10_1*(-tmp_32 - tmp_33);
      real_t tmp_35 = -tmp_12 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_36 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_37 = 6/tmp_36;
      real_t tmp_38 = tmp_30*tmp_37;
      real_t tmp_39 = tmp_1*tmp_25;
      real_t tmp_40 = 0.5*p_affine_10_0*(tmp_2*tmp_39 + tmp_26*tmp_7 + tmp_27*tmp_4 + tmp_39*tmp_6) + 0.5*p_affine_10_1*(tmp_31*tmp_5 + tmp_33*tmp_7);
      real_t tmp_41 = 0.5*tmp_36;
      real_t tmp_42 = tmp_8*(tmp_10 + 0.78867513459481287*tmp_9);
      real_t tmp_43 = tmp_2*tmp_42;
      real_t tmp_44 = tmp_8*(0.78867513459481287*tmp_13 + tmp_14);
      real_t tmp_45 = tmp_4*tmp_44;
      real_t tmp_46 = tmp_43 + tmp_45;
      real_t tmp_47 = tmp_46 - 1.0/3.0;
      real_t tmp_48 = tmp_1*tmp_42;
      real_t tmp_49 = tmp_20*tmp_44;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_50 - 1.0/3.0;
      real_t tmp_52 = p_affine_10_1*(tmp_1*tmp_47 + tmp_51*tmp_6);
      real_t tmp_53 = tmp_4*tmp_51 + tmp_47*tmp_7;
      real_t tmp_54 = -tmp_43 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_55 = tmp_37*tmp_53;
      real_t tmp_56 = 0.5*tmp_36;
      real_t tmp_57 = 0.25*tmp_8;
      real_t tmp_58 = tmp_4*tmp_57;
      real_t tmp_59 = 0.5*p_affine_10_0*tmp_26 + 0.5*p_affine_10_1*tmp_33;
      real_t tmp_60 = tmp_20*tmp_57;
      real_t tmp_61 = 0.5*p_affine_10_0*tmp_27 + 0.5*p_affine_10_1*tmp_32;
      real_t a_0_0 = tmp_41*(-tmp_24*tmp_29 - tmp_30*tmp_34 + tmp_35*tmp_38 - tmp_35*tmp_40) + tmp_56*(-tmp_29*tmp_52 - tmp_34*tmp_53 - tmp_40*tmp_54 + tmp_54*tmp_55);
      real_t a_0_1 = tmp_41*(tmp_17*tmp_38 - tmp_17*tmp_40 - tmp_24*tmp_58 - tmp_30*tmp_59) + tmp_56*(-tmp_40*tmp_46 + tmp_46*tmp_55 - tmp_52*tmp_58 - tmp_53*tmp_59);
      real_t a_0_2 = tmp_41*(tmp_22*tmp_38 - tmp_22*tmp_40 - tmp_24*tmp_60 - tmp_30*tmp_61) + tmp_56*(-tmp_40*tmp_50 + tmp_50*tmp_55 - tmp_52*tmp_60 - tmp_53*tmp_61);
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
      real_t tmp_10 = p_affine_6_1 + 0.21132486540518713*tmp_9;
      real_t tmp_11 = tmp_8*(tmp_10 + tmp_3);
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + 0.21132486540518713*tmp_12;
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
      real_t tmp_29 = 0.5*tmp_28;
      real_t tmp_30 = tmp_15*tmp_7 + tmp_17*tmp_4;
      real_t tmp_31 = 1.0*tmp_23;
      real_t tmp_32 = tmp_22*tmp_31;
      real_t tmp_33 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_34 = tmp_31*tmp_33;
      real_t tmp_35 = 0.5*p_affine_10_0*tmp_28 + 0.5*p_affine_10_1*(-tmp_32 - tmp_34);
      real_t tmp_36 = tmp_23*(tmp_10 + tmp_19);
      real_t tmp_37 = tmp_22*tmp_36;
      real_t tmp_38 = tmp_33*tmp_36;
      real_t tmp_39 = tmp_23*(tmp_13 + tmp_21);
      real_t tmp_40 = tmp_20*tmp_39;
      real_t tmp_41 = tmp_26*tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = std::abs(std::pow((tmp_12*tmp_12) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_44 = 6/tmp_43;
      real_t tmp_45 = tmp_30*tmp_44;
      real_t tmp_46 = 1.0*tmp_8;
      real_t tmp_47 = 0.5*tmp_8;
      real_t tmp_48 = tmp_1*tmp_47;
      real_t tmp_49 = tmp_4*tmp_47;
      real_t tmp_50 = 0.5*p_affine_10_0*(tmp_16*tmp_49 + tmp_2*tmp_48 + tmp_48*tmp_6 + tmp_49*tmp_7) + 0.5*p_affine_10_1*(tmp_2*tmp_46*tmp_7 + tmp_46*tmp_5);
      real_t tmp_51 = 0.5*tmp_43;
      real_t tmp_52 = p_affine_6_1 + 0.78867513459481287*tmp_9;
      real_t tmp_53 = tmp_8*(tmp_3 + tmp_52);
      real_t tmp_54 = p_affine_6_0 + 0.78867513459481287*tmp_12;
      real_t tmp_55 = tmp_8*(tmp_0 + tmp_54);
      real_t tmp_56 = tmp_2*tmp_53 + tmp_4*tmp_55 - 1.0/3.0;
      real_t tmp_57 = tmp_1*tmp_53 + tmp_16*tmp_55 - 1.0/3.0;
      real_t tmp_58 = p_affine_10_1*(tmp_1*tmp_56 + tmp_57*tmp_6);
      real_t tmp_59 = tmp_4*tmp_57 + tmp_56*tmp_7;
      real_t tmp_60 = tmp_23*(tmp_19 + tmp_52);
      real_t tmp_61 = tmp_22*tmp_60;
      real_t tmp_62 = tmp_33*tmp_60;
      real_t tmp_63 = tmp_23*(tmp_21 + tmp_54);
      real_t tmp_64 = tmp_20*tmp_63;
      real_t tmp_65 = tmp_26*tmp_63;
      real_t tmp_66 = -tmp_61 - tmp_62 - tmp_64 - tmp_65 + 1;
      real_t tmp_67 = tmp_44*tmp_59;
      real_t tmp_68 = 0.5*tmp_43;
      real_t tmp_69 = 0.25*tmp_23;
      real_t tmp_70 = tmp_20*tmp_69;
      real_t tmp_71 = 0.5*p_affine_10_0*tmp_25 + 0.5*p_affine_10_1*tmp_34;
      real_t tmp_72 = tmp_38 + tmp_40;
      real_t tmp_73 = tmp_62 + tmp_64;
      real_t tmp_74 = tmp_26*tmp_69;
      real_t tmp_75 = 0.5*p_affine_10_0*tmp_27 + 0.5*p_affine_10_1*tmp_32;
      real_t tmp_76 = tmp_37 + tmp_41;
      real_t tmp_77 = tmp_61 + tmp_65;
      real_t a_0_0 = tmp_51*(-tmp_18*tmp_29 - tmp_30*tmp_35 - tmp_42*tmp_45 + tmp_42*tmp_50) + tmp_68*(-tmp_29*tmp_58 - tmp_35*tmp_59 + tmp_50*tmp_66 - tmp_66*tmp_67);
      real_t a_0_1 = tmp_51*(-tmp_18*tmp_70 - tmp_30*tmp_71 - tmp_45*tmp_72 + tmp_50*tmp_72) + tmp_68*(tmp_50*tmp_73 - tmp_58*tmp_70 - tmp_59*tmp_71 - tmp_67*tmp_73);
      real_t a_0_2 = tmp_51*(-tmp_18*tmp_74 - tmp_30*tmp_75 - tmp_45*tmp_76 + tmp_50*tmp_76) + tmp_68*(tmp_50*tmp_77 - tmp_58*tmp_74 - tmp_59*tmp_75 - tmp_67*tmp_77);
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
      real_t tmp_5 = tmp_1*tmp_4;
      real_t tmp_6 = p_affine_2_0 + tmp_0;
      real_t tmp_7 = p_affine_1_1 + tmp_3;
      real_t tmp_8 = 1.0 / (tmp_5 - tmp_6*tmp_7);
      real_t tmp_9 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + tmp_3;
      real_t tmp_11 = tmp_8*(tmp_10 + 0.21132486540518713*tmp_9);
      real_t tmp_12 = tmp_11*tmp_2;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_0;
      real_t tmp_15 = tmp_8*(0.21132486540518713*tmp_13 + tmp_14);
      real_t tmp_16 = tmp_15*tmp_4;
      real_t tmp_17 = tmp_12 + tmp_16;
      real_t tmp_18 = tmp_17 - 1.0/3.0;
      real_t tmp_19 = tmp_1*tmp_11;
      real_t tmp_20 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_21 = tmp_15*tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_22 - 1.0/3.0;
      real_t tmp_24 = tmp_1*tmp_18 + tmp_23*tmp_6;
      real_t tmp_25 = 0.5*tmp_8;
      real_t tmp_26 = tmp_25*tmp_4;
      real_t tmp_27 = tmp_20*tmp_25;
      real_t tmp_28 = -tmp_26 - tmp_27;
      real_t tmp_29 = p_affine_10_1*tmp_28;
      real_t tmp_30 = 1.0*tmp_8;
      real_t tmp_31 = tmp_1*tmp_30;
      real_t tmp_32 = tmp_2*tmp_30;
      real_t tmp_33 = p_affine_10_0*tmp_28 + p_affine_10_1*(-tmp_31 - tmp_32);
      real_t tmp_34 = tmp_18*tmp_7 + tmp_23*tmp_4;
      real_t tmp_35 = -tmp_12 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_36 = std::abs(std::pow((tmp_13*tmp_13) + (tmp_9*tmp_9), 1.0/2.0));
      real_t tmp_37 = 24/tmp_36;
      real_t tmp_38 = tmp_34*tmp_37;
      real_t tmp_39 = tmp_1*tmp_25;
      real_t tmp_40 = p_affine_10_0*(tmp_2*tmp_39 + tmp_26*tmp_7 + tmp_27*tmp_4 + tmp_39*tmp_6) + p_affine_10_1*(tmp_30*tmp_5 + tmp_32*tmp_7);
      real_t tmp_41 = 0.5*tmp_36;
      real_t tmp_42 = tmp_8*(tmp_10 + 0.78867513459481287*tmp_9);
      real_t tmp_43 = tmp_2*tmp_42;
      real_t tmp_44 = tmp_8*(0.78867513459481287*tmp_13 + tmp_14);
      real_t tmp_45 = tmp_4*tmp_44;
      real_t tmp_46 = tmp_43 + tmp_45;
      real_t tmp_47 = tmp_46 - 1.0/3.0;
      real_t tmp_48 = tmp_1*tmp_42;
      real_t tmp_49 = tmp_20*tmp_44;
      real_t tmp_50 = tmp_48 + tmp_49;
      real_t tmp_51 = tmp_50 - 1.0/3.0;
      real_t tmp_52 = tmp_1*tmp_47 + tmp_51*tmp_6;
      real_t tmp_53 = tmp_4*tmp_51 + tmp_47*tmp_7;
      real_t tmp_54 = -tmp_43 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_55 = tmp_37*tmp_53;
      real_t tmp_56 = 0.5*tmp_36;
      real_t tmp_57 = p_affine_10_1*tmp_26;
      real_t tmp_58 = p_affine_10_0*tmp_26 + p_affine_10_1*tmp_32;
      real_t tmp_59 = p_affine_10_1*tmp_27;
      real_t tmp_60 = p_affine_10_0*tmp_27 + p_affine_10_1*tmp_31;
      real_t a_0_0 = tmp_41*(-tmp_24*tmp_29 - tmp_33*tmp_34 + tmp_35*tmp_38 - tmp_35*tmp_40) + tmp_56*(-tmp_29*tmp_52 - tmp_33*tmp_53 - tmp_40*tmp_54 + tmp_54*tmp_55);
      real_t a_0_1 = tmp_41*(tmp_17*tmp_38 - tmp_17*tmp_40 - tmp_24*tmp_57 - tmp_34*tmp_58) + tmp_56*(-tmp_40*tmp_46 + tmp_46*tmp_55 - tmp_52*tmp_57 - tmp_53*tmp_58);
      real_t a_0_2 = tmp_41*(tmp_22*tmp_38 - tmp_22*tmp_40 - tmp_24*tmp_59 - tmp_34*tmp_60) + tmp_56*(-tmp_40*tmp_50 + tmp_50*tmp_55 - tmp_52*tmp_59 - tmp_53*tmp_60);
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




class EDGEpsilonFormEDGEDG : public hyteg::dg::DGForm2D
{
 public:
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
      real_t tmp_13 = tmp_11*(0.21132486540518713*tmp_1 + tmp_12);
      real_t tmp_14 = p_affine_6_0 + tmp_3;
      real_t tmp_15 = tmp_11*(0.21132486540518713*tmp_0 + tmp_14);
      real_t tmp_16 = tmp_13*tmp_5 + tmp_15*tmp_7 - 1.0/3.0;
      real_t tmp_17 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_18 = tmp_13*tmp_4 + tmp_15*tmp_17 - 1.0/3.0;
      real_t tmp_19 = tmp_16*tmp_4 + tmp_18*tmp_9;
      real_t tmp_20 = tmp_10*tmp_16 + tmp_18*tmp_7;
      real_t tmp_21 = 6/tmp_2;
      real_t tmp_22 = 1.0*tmp_11;
      real_t tmp_23 = tmp_22*tmp_8;
      real_t tmp_24 = 0.5*tmp_11;
      real_t tmp_25 = tmp_24*tmp_4;
      real_t tmp_26 = tmp_24*tmp_7;
      real_t tmp_27 = tmp_10*tmp_26 + tmp_17*tmp_26 + tmp_25*tmp_5 + tmp_25*tmp_9;
      real_t tmp_28 = 1.0*p_affine_10_0*(tmp_17*tmp_22*tmp_9 + tmp_23) + 1.0*p_affine_10_1*tmp_27;
      real_t tmp_29 = 1.0*p_affine_10_0*tmp_27 + 1.0*p_affine_10_1*(tmp_10*tmp_22*tmp_5 + tmp_23);
      real_t tmp_30 = tmp_11*(0.78867513459481287*tmp_1 + tmp_12);
      real_t tmp_31 = tmp_11*(0.78867513459481287*tmp_0 + tmp_14);
      real_t tmp_32 = tmp_30*tmp_5 + tmp_31*tmp_7 - 1.0/3.0;
      real_t tmp_33 = tmp_17*tmp_31 + tmp_30*tmp_4 - 1.0/3.0;
      real_t tmp_34 = tmp_32*tmp_4 + tmp_33*tmp_9;
      real_t tmp_35 = tmp_10*tmp_32 + tmp_33*tmp_7;
      real_t a_0_0 = 0.5*tmp_2*(-tmp_19*tmp_28 - tmp_20*tmp_29 + tmp_21*((tmp_19*tmp_19) + (tmp_20*tmp_20))) + 0.5*tmp_2*(tmp_21*((tmp_34*tmp_34) + (tmp_35*tmp_35)) - tmp_28*tmp_34 - tmp_29*tmp_35);
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
      real_t tmp_12 = p_affine_6_1 + 0.21132486540518713*tmp_1;
      real_t tmp_13 = tmp_11*(tmp_12 + tmp_6);
      real_t tmp_14 = p_affine_6_0 + 0.21132486540518713*tmp_0;
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
      real_t tmp_36 = 0.5*p_affine_10_0*(tmp_25*tmp_28*tmp_30 + tmp_29) + 0.5*p_affine_10_1*tmp_35;
      real_t tmp_37 = tmp_27*(tmp_12 + tmp_22);
      real_t tmp_38 = tmp_27*(tmp_14 + tmp_20);
      real_t tmp_39 = tmp_23*tmp_38 + tmp_33*tmp_37 - 1.0/3.0;
      real_t tmp_40 = tmp_21*tmp_37 + tmp_30*tmp_38 - 1.0/3.0;
      real_t tmp_41 = tmp_21*tmp_39 + tmp_25*tmp_40;
      real_t tmp_42 = 1.0*tmp_11;
      real_t tmp_43 = tmp_42*tmp_8;
      real_t tmp_44 = 0.5*tmp_11;
      real_t tmp_45 = tmp_4*tmp_44;
      real_t tmp_46 = tmp_44*tmp_7;
      real_t tmp_47 = tmp_10*tmp_46 + tmp_17*tmp_46 + tmp_45*tmp_5 + tmp_45*tmp_9;
      real_t tmp_48 = 0.5*p_affine_10_0*(tmp_17*tmp_42*tmp_9 + tmp_43) + 0.5*p_affine_10_1*tmp_47;
      real_t tmp_49 = tmp_10*tmp_16 + tmp_18*tmp_7;
      real_t tmp_50 = 0.5*p_affine_10_0*tmp_35 + 0.5*p_affine_10_1*(tmp_26*tmp_28*tmp_33 + tmp_29);
      real_t tmp_51 = tmp_23*tmp_40 + tmp_26*tmp_39;
      real_t tmp_52 = 0.5*p_affine_10_0*tmp_47 + 0.5*p_affine_10_1*(tmp_10*tmp_42*tmp_5 + tmp_43);
      real_t tmp_53 = 6/tmp_2;
      real_t tmp_54 = p_affine_6_1 + 0.78867513459481287*tmp_1;
      real_t tmp_55 = tmp_11*(tmp_54 + tmp_6);
      real_t tmp_56 = p_affine_6_0 + 0.78867513459481287*tmp_0;
      real_t tmp_57 = tmp_11*(tmp_3 + tmp_56);
      real_t tmp_58 = tmp_5*tmp_55 + tmp_57*tmp_7 - 1.0/3.0;
      real_t tmp_59 = tmp_17*tmp_57 + tmp_4*tmp_55 - 1.0/3.0;
      real_t tmp_60 = tmp_4*tmp_58 + tmp_59*tmp_9;
      real_t tmp_61 = tmp_27*(tmp_22 + tmp_54);
      real_t tmp_62 = tmp_27*(tmp_20 + tmp_56);
      real_t tmp_63 = tmp_23*tmp_62 + tmp_33*tmp_61 - 1.0/3.0;
      real_t tmp_64 = tmp_21*tmp_61 + tmp_30*tmp_62 - 1.0/3.0;
      real_t tmp_65 = tmp_21*tmp_63 + tmp_25*tmp_64;
      real_t tmp_66 = tmp_10*tmp_58 + tmp_59*tmp_7;
      real_t tmp_67 = tmp_23*tmp_64 + tmp_26*tmp_63;
      real_t a_0_0 = 0.5*tmp_2*(tmp_19*tmp_36 - tmp_41*tmp_48 + tmp_49*tmp_50 - tmp_51*tmp_52 - tmp_53*(tmp_19*tmp_41 + tmp_49*tmp_51)) + 0.5*tmp_2*(tmp_36*tmp_60 - tmp_48*tmp_65 + tmp_50*tmp_66 - tmp_52*tmp_67 - tmp_53*(tmp_60*tmp_65 + tmp_66*tmp_67));
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
      real_t tmp_8 = tmp_4*tmp_7;
      real_t tmp_9 = p_affine_2_0 + tmp_3;
      real_t tmp_10 = p_affine_1_1 + tmp_6;
      real_t tmp_11 = 1.0 / (-tmp_10*tmp_9 + tmp_8);
      real_t tmp_12 = p_affine_6_1 + tmp_6;
      real_t tmp_13 = tmp_11*(0.21132486540518713*tmp_1 + tmp_12);
      real_t tmp_14 = p_affine_6_0 + tmp_3;
      real_t tmp_15 = tmp_11*(0.21132486540518713*tmp_0 + tmp_14);
      real_t tmp_16 = tmp_13*tmp_5 + tmp_15*tmp_7 - 1.0/3.0;
      real_t tmp_17 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_18 = tmp_13*tmp_4 + tmp_15*tmp_17 - 1.0/3.0;
      real_t tmp_19 = tmp_16*tmp_4 + tmp_18*tmp_9;
      real_t tmp_20 = tmp_10*tmp_16 + tmp_18*tmp_7;
      real_t tmp_21 = 24/tmp_2;
      real_t tmp_22 = 1.0*tmp_11;
      real_t tmp_23 = tmp_22*tmp_8;
      real_t tmp_24 = 0.5*tmp_11;
      real_t tmp_25 = tmp_24*tmp_4;
      real_t tmp_26 = tmp_24*tmp_7;
      real_t tmp_27 = tmp_10*tmp_26 + tmp_17*tmp_26 + tmp_25*tmp_5 + tmp_25*tmp_9;
      real_t tmp_28 = 2*p_affine_10_0*(tmp_17*tmp_22*tmp_9 + tmp_23) + 2*p_affine_10_1*tmp_27;
      real_t tmp_29 = 2*p_affine_10_0*tmp_27 + 2*p_affine_10_1*(tmp_10*tmp_22*tmp_5 + tmp_23);
      real_t tmp_30 = tmp_11*(0.78867513459481287*tmp_1 + tmp_12);
      real_t tmp_31 = tmp_11*(0.78867513459481287*tmp_0 + tmp_14);
      real_t tmp_32 = tmp_30*tmp_5 + tmp_31*tmp_7 - 1.0/3.0;
      real_t tmp_33 = tmp_17*tmp_31 + tmp_30*tmp_4 - 1.0/3.0;
      real_t tmp_34 = tmp_32*tmp_4 + tmp_33*tmp_9;
      real_t tmp_35 = tmp_10*tmp_32 + tmp_33*tmp_7;
      real_t a_0_0 = 0.5*tmp_2*(-tmp_19*tmp_28 - tmp_20*tmp_29 + tmp_21*((tmp_19*tmp_19) + (tmp_20*tmp_20))) + 0.5*tmp_2*(tmp_21*((tmp_34*tmp_34) + (tmp_35*tmp_35)) - tmp_28*tmp_34 - tmp_29*tmp_35);
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


} // dg
} // hyteg

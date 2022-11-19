
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

class EGNIPGVectorLaplaceFormP1E_0 : public hyteg::dg::DGForm
{
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
      real_t tmp_28 = 5/tmp_27;
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
      real_t tmp_37 = 5/tmp_36;
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
      real_t tmp_20 = 0.59231721264047354*tmp_1*(tmp_18 - 1.0/3.0) + 0.59231721264047354*tmp_4*(tmp_19 - 1.0/3.0);
      real_t tmp_21 = tmp_5*(0.23076534494715845*tmp_6 + tmp_7);
      real_t tmp_22 = tmp_1*tmp_21;
      real_t tmp_23 = tmp_10*tmp_21;
      real_t tmp_24 = tmp_5*(0.23076534494715845*tmp_12 + tmp_13);
      real_t tmp_25 = tmp_24*tmp_3;
      real_t tmp_26 = tmp_16*tmp_24;
      real_t tmp_27 = tmp_23 + tmp_25;
      real_t tmp_28 = tmp_22 + tmp_26;
      real_t tmp_29 = 1.1965716762484155*tmp_1*(tmp_27 - 1.0/3.0) + 1.1965716762484155*tmp_4*(tmp_28 - 1.0/3.0);
      real_t tmp_30 = tmp_5*(0.5*tmp_6 + tmp_7);
      real_t tmp_31 = tmp_1*tmp_30;
      real_t tmp_32 = tmp_10*tmp_30;
      real_t tmp_33 = tmp_5*(0.5*tmp_12 + tmp_13);
      real_t tmp_34 = tmp_3*tmp_33;
      real_t tmp_35 = tmp_16*tmp_33;
      real_t tmp_36 = tmp_32 + tmp_34;
      real_t tmp_37 = tmp_31 + tmp_35;
      real_t tmp_38 = 1.4222222222222225*tmp_1*(tmp_36 - 1.0/3.0) + 1.4222222222222225*tmp_4*(tmp_37 - 1.0/3.0);
      real_t tmp_39 = tmp_5*(0.7692346550528415*tmp_6 + tmp_7);
      real_t tmp_40 = tmp_1*tmp_39;
      real_t tmp_41 = tmp_10*tmp_39;
      real_t tmp_42 = tmp_5*(0.7692346550528415*tmp_12 + tmp_13);
      real_t tmp_43 = tmp_3*tmp_42;
      real_t tmp_44 = tmp_16*tmp_42;
      real_t tmp_45 = tmp_41 + tmp_43;
      real_t tmp_46 = tmp_40 + tmp_44;
      real_t tmp_47 = 1.1965716762484155*tmp_1*(tmp_45 - 1.0/3.0) + 1.1965716762484155*tmp_4*(tmp_46 - 1.0/3.0);
      real_t tmp_48 = tmp_5*(0.95308992296933193*tmp_6 + tmp_7);
      real_t tmp_49 = tmp_1*tmp_48;
      real_t tmp_50 = tmp_10*tmp_48;
      real_t tmp_51 = tmp_5*(0.95308992296933193*tmp_12 + tmp_13);
      real_t tmp_52 = tmp_3*tmp_51;
      real_t tmp_53 = tmp_16*tmp_51;
      real_t tmp_54 = tmp_50 + tmp_52;
      real_t tmp_55 = tmp_49 + tmp_53;
      real_t tmp_56 = 0.59231721264047354*tmp_1*(tmp_54 - 1.0/3.0) + 0.59231721264047354*tmp_4*(tmp_55 - 1.0/3.0);
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
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
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
      real_t tmp_22 = tmp_1*tmp_21 + tmp_14*tmp_19 + tmp_20*tmp_5;
      real_t tmp_23 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_5);
      real_t tmp_24 = tmp_18*(tmp_1*tmp_9 - tmp_10*tmp_14);
      real_t tmp_25 = tmp_18*(tmp_13*tmp_14 - tmp_5*tmp_9);
      real_t tmp_26 = tmp_1*tmp_25 + tmp_14*tmp_23 + tmp_24*tmp_5;
      real_t tmp_27 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_6);
      real_t tmp_28 = tmp_18*(tmp_10*tmp_11 - tmp_6*tmp_9);
      real_t tmp_29 = tmp_18*(-tmp_11*tmp_13 + tmp_3*tmp_9);
      real_t tmp_30 = tmp_1*tmp_29 + tmp_14*tmp_27 + tmp_28*tmp_5;
      real_t tmp_31 = p_affine_0_0*p_affine_1_1;
      real_t tmp_32 = p_affine_0_0*p_affine_1_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_2;
      real_t tmp_34 = p_affine_0_1*p_affine_1_0;
      real_t tmp_35 = p_affine_0_1*p_affine_1_2;
      real_t tmp_36 = p_affine_2_2*p_affine_3_0;
      real_t tmp_37 = p_affine_0_2*p_affine_1_0;
      real_t tmp_38 = p_affine_0_2*p_affine_1_1;
      real_t tmp_39 = p_affine_2_0*p_affine_3_1;
      real_t tmp_40 = p_affine_2_2*p_affine_3_1;
      real_t tmp_41 = p_affine_2_0*p_affine_3_2;
      real_t tmp_42 = p_affine_2_1*p_affine_3_0;
      real_t tmp_43 = std::abs(p_affine_0_0*tmp_33 - p_affine_0_0*tmp_40 + p_affine_0_1*tmp_36 - p_affine_0_1*tmp_41 + p_affine_0_2*tmp_39 - p_affine_0_2*tmp_42 - p_affine_1_0*tmp_33 + p_affine_1_0*tmp_40 - p_affine_1_1*tmp_36 + p_affine_1_1*tmp_41 - p_affine_1_2*tmp_39 + p_affine_1_2*tmp_42 + p_affine_2_0*tmp_35 - p_affine_2_0*tmp_38 - p_affine_2_1*tmp_32 + p_affine_2_1*tmp_37 + p_affine_2_2*tmp_31 - p_affine_2_2*tmp_34 - p_affine_3_0*tmp_35 + p_affine_3_0*tmp_38 + p_affine_3_1*tmp_32 - p_affine_3_1*tmp_37 - p_affine_3_2*tmp_31 + p_affine_3_2*tmp_34);
      real_t tmp_44 = tmp_43*(tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25) + tmp_30*(-tmp_27 - tmp_28 - tmp_29));
      real_t tmp_45 = tmp_43*(tmp_21*tmp_22 + tmp_25*tmp_26 + tmp_29*tmp_30);
      real_t tmp_46 = tmp_43*(tmp_20*tmp_22 + tmp_24*tmp_26 + tmp_28*tmp_30);
      real_t tmp_47 = tmp_43*(tmp_19*tmp_22 + tmp_23*tmp_26 + tmp_27*tmp_30);
      real_t a_0_0 = 0.1666666666666668*tmp_44;
      real_t a_1_0 = 0.1666666666666668*tmp_45;
      real_t a_2_0 = 0.1666666666666668*tmp_46;
      real_t a_3_0 = 0.1666666666666668*tmp_47;
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
      real_t tmp_2 = -p_affine_8_2;
      real_t tmp_3 = p_affine_9_2 + tmp_2;
      real_t tmp_4 = p_affine_10_2 + tmp_2;
      real_t tmp_5 = -p_affine_0_2;
      real_t tmp_6 = p_affine_8_2 + tmp_5;
      real_t tmp_7 = 0.031405749086161582*tmp_3 + 0.93718850182767688*tmp_4 + tmp_6;
      real_t tmp_8 = p_affine_2_0 + tmp_0;
      real_t tmp_9 = -p_affine_0_1;
      real_t tmp_10 = p_affine_3_1 + tmp_9;
      real_t tmp_11 = p_affine_3_0 + tmp_0;
      real_t tmp_12 = p_affine_2_1 + tmp_9;
      real_t tmp_13 = p_affine_3_2 + tmp_5;
      real_t tmp_14 = tmp_12*tmp_13;
      real_t tmp_15 = p_affine_1_2 + tmp_5;
      real_t tmp_16 = tmp_10*tmp_15;
      real_t tmp_17 = p_affine_1_1 + tmp_9;
      real_t tmp_18 = p_affine_2_2 + tmp_5;
      real_t tmp_19 = tmp_17*tmp_18;
      real_t tmp_20 = tmp_10*tmp_18;
      real_t tmp_21 = tmp_13*tmp_17;
      real_t tmp_22 = tmp_12*tmp_15;
      real_t tmp_23 = 1.0 / (tmp_1*tmp_14 - tmp_1*tmp_20 + tmp_11*tmp_19 - tmp_11*tmp_22 + tmp_16*tmp_8 - tmp_21*tmp_8);
      real_t tmp_24 = tmp_23*(tmp_10*tmp_8 - tmp_11*tmp_12);
      real_t tmp_25 = tmp_24*tmp_7;
      real_t tmp_26 = -p_affine_8_1;
      real_t tmp_27 = p_affine_9_1 + tmp_26;
      real_t tmp_28 = p_affine_10_1 + tmp_26;
      real_t tmp_29 = p_affine_8_1 + tmp_9;
      real_t tmp_30 = 0.031405749086161582*tmp_27 + 0.93718850182767688*tmp_28 + tmp_29;
      real_t tmp_31 = tmp_23*(tmp_11*tmp_18 - tmp_13*tmp_8);
      real_t tmp_32 = tmp_30*tmp_31;
      real_t tmp_33 = -p_affine_8_0;
      real_t tmp_34 = p_affine_9_0 + tmp_33;
      real_t tmp_35 = p_affine_10_0 + tmp_33;
      real_t tmp_36 = p_affine_8_0 + tmp_0;
      real_t tmp_37 = 0.031405749086161582*tmp_34 + 0.93718850182767688*tmp_35 + tmp_36;
      real_t tmp_38 = tmp_23*(tmp_14 - tmp_20);
      real_t tmp_39 = tmp_37*tmp_38;
      real_t tmp_40 = tmp_25 + tmp_32 + tmp_39;
      real_t tmp_41 = tmp_23*(-tmp_1*tmp_10 + tmp_11*tmp_17);
      real_t tmp_42 = tmp_41*tmp_7;
      real_t tmp_43 = tmp_23*(tmp_1*tmp_13 - tmp_11*tmp_15);
      real_t tmp_44 = tmp_30*tmp_43;
      real_t tmp_45 = tmp_23*(tmp_16 - tmp_21);
      real_t tmp_46 = tmp_37*tmp_45;
      real_t tmp_47 = tmp_42 + tmp_44 + tmp_46;
      real_t tmp_48 = tmp_23*(tmp_1*tmp_12 - tmp_17*tmp_8);
      real_t tmp_49 = tmp_48*tmp_7;
      real_t tmp_50 = tmp_23*(-tmp_1*tmp_18 + tmp_15*tmp_8);
      real_t tmp_51 = tmp_30*tmp_50;
      real_t tmp_52 = tmp_23*(tmp_19 - tmp_22);
      real_t tmp_53 = tmp_37*tmp_52;
      real_t tmp_54 = tmp_49 + tmp_51 + tmp_53;
      real_t tmp_55 = tmp_1*(tmp_40 - 1.0/4.0) + tmp_11*(tmp_54 - 1.0/4.0) + tmp_8*(tmp_47 - 1.0/4.0);
      real_t tmp_56 = 0.5*p_affine_13_0*(-tmp_38 - tmp_45 - tmp_52) + 0.5*p_affine_13_1*(-tmp_31 - tmp_43 - tmp_50) + 0.5*p_affine_13_2*(-tmp_24 - tmp_41 - tmp_48);
      real_t tmp_57 = -tmp_25 - tmp_32 - tmp_39 - tmp_42 - tmp_44 - tmp_46 - tmp_49 - tmp_51 - tmp_53 + 1;
      real_t tmp_58 = 0.5*p_affine_13_0*(tmp_1*tmp_38 + tmp_11*tmp_52 + tmp_45*tmp_8) + 0.5*p_affine_13_1*(tmp_1*tmp_31 + tmp_11*tmp_50 + tmp_43*tmp_8) + 0.5*p_affine_13_2*(tmp_1*tmp_24 + tmp_11*tmp_48 + tmp_41*tmp_8);
      real_t tmp_59 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_60 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_61 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_62 = (std::abs(tmp_28*tmp_60 - tmp_35*tmp_59)*std::abs(tmp_28*tmp_60 - tmp_35*tmp_59)) + (std::abs(tmp_28*tmp_61 - tmp_4*tmp_59)*std::abs(tmp_28*tmp_61 - tmp_4*tmp_59)) + (std::abs(tmp_35*tmp_61 - tmp_4*tmp_60)*std::abs(tmp_35*tmp_61 - tmp_4*tmp_60));
      real_t tmp_63 = 5.0*std::pow(tmp_62, -0.25);
      real_t tmp_64 = tmp_55*tmp_63;
      real_t tmp_65 = 1.0*std::pow(tmp_62, 1.0/2.0);
      real_t tmp_66 = 0.0068572537431980923*tmp_65;
      real_t tmp_67 = 0.19601935860219369*tmp_3 + 0.60796128279561268*tmp_4 + tmp_6;
      real_t tmp_68 = tmp_24*tmp_67;
      real_t tmp_69 = 0.19601935860219369*tmp_27 + 0.60796128279561268*tmp_28 + tmp_29;
      real_t tmp_70 = tmp_31*tmp_69;
      real_t tmp_71 = 0.19601935860219369*tmp_34 + 0.60796128279561268*tmp_35 + tmp_36;
      real_t tmp_72 = tmp_38*tmp_71;
      real_t tmp_73 = tmp_68 + tmp_70 + tmp_72;
      real_t tmp_74 = tmp_41*tmp_67;
      real_t tmp_75 = tmp_43*tmp_69;
      real_t tmp_76 = tmp_45*tmp_71;
      real_t tmp_77 = tmp_74 + tmp_75 + tmp_76;
      real_t tmp_78 = tmp_48*tmp_67;
      real_t tmp_79 = tmp_50*tmp_69;
      real_t tmp_80 = tmp_52*tmp_71;
      real_t tmp_81 = tmp_78 + tmp_79 + tmp_80;
      real_t tmp_82 = tmp_1*(tmp_73 - 1.0/4.0) + tmp_11*(tmp_81 - 1.0/4.0) + tmp_8*(tmp_77 - 1.0/4.0);
      real_t tmp_83 = -tmp_68 - tmp_70 - tmp_72 - tmp_74 - tmp_75 - tmp_76 - tmp_78 - tmp_79 - tmp_80 + 1;
      real_t tmp_84 = tmp_63*tmp_82;
      real_t tmp_85 = 0.037198804536718075*tmp_65;
      real_t tmp_86 = 0.37605877282253791*tmp_3 + 0.039308471900058539*tmp_4 + tmp_6;
      real_t tmp_87 = tmp_24*tmp_86;
      real_t tmp_88 = 0.37605877282253791*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29;
      real_t tmp_89 = tmp_31*tmp_88;
      real_t tmp_90 = 0.37605877282253791*tmp_34 + 0.039308471900058539*tmp_35 + tmp_36;
      real_t tmp_91 = tmp_38*tmp_90;
      real_t tmp_92 = tmp_87 + tmp_89 + tmp_91;
      real_t tmp_93 = tmp_41*tmp_86;
      real_t tmp_94 = tmp_43*tmp_88;
      real_t tmp_95 = tmp_45*tmp_90;
      real_t tmp_96 = tmp_93 + tmp_94 + tmp_95;
      real_t tmp_97 = tmp_48*tmp_86;
      real_t tmp_98 = tmp_50*tmp_88;
      real_t tmp_99 = tmp_52*tmp_90;
      real_t tmp_100 = tmp_97 + tmp_98 + tmp_99;
      real_t tmp_101 = tmp_1*(tmp_92 - 1.0/4.0) + tmp_11*(tmp_100 - 1.0/4.0) + tmp_8*(tmp_96 - 1.0/4.0);
      real_t tmp_102 = -tmp_87 - tmp_89 - tmp_91 - tmp_93 - tmp_94 - tmp_95 - tmp_97 - tmp_98 - tmp_99 + 1;
      real_t tmp_103 = tmp_101*tmp_63;
      real_t tmp_104 = 0.020848748529055869*tmp_65;
      real_t tmp_105 = 0.78764240869137092*tmp_3 + 0.1711304259088916*tmp_4 + tmp_6;
      real_t tmp_106 = tmp_105*tmp_24;
      real_t tmp_107 = 0.78764240869137092*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29;
      real_t tmp_108 = tmp_107*tmp_31;
      real_t tmp_109 = 0.78764240869137092*tmp_34 + 0.1711304259088916*tmp_35 + tmp_36;
      real_t tmp_110 = tmp_109*tmp_38;
      real_t tmp_111 = tmp_106 + tmp_108 + tmp_110;
      real_t tmp_112 = tmp_105*tmp_41;
      real_t tmp_113 = tmp_107*tmp_43;
      real_t tmp_114 = tmp_109*tmp_45;
      real_t tmp_115 = tmp_112 + tmp_113 + tmp_114;
      real_t tmp_116 = tmp_105*tmp_48;
      real_t tmp_117 = tmp_107*tmp_50;
      real_t tmp_118 = tmp_109*tmp_52;
      real_t tmp_119 = tmp_116 + tmp_117 + tmp_118;
      real_t tmp_120 = tmp_1*(tmp_111 - 1.0/4.0) + tmp_11*(tmp_119 - 1.0/4.0) + tmp_8*(tmp_115 - 1.0/4.0);
      real_t tmp_121 = -tmp_106 - tmp_108 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_116 - tmp_117 - tmp_118 + 1;
      real_t tmp_122 = tmp_120*tmp_63;
      real_t tmp_123 = 0.019202922745021479*tmp_65;
      real_t tmp_124 = 0.58463275527740355*tmp_3 + 0.37605877282253791*tmp_4 + tmp_6;
      real_t tmp_125 = tmp_124*tmp_24;
      real_t tmp_126 = 0.58463275527740355*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29;
      real_t tmp_127 = tmp_126*tmp_31;
      real_t tmp_128 = 0.58463275527740355*tmp_34 + 0.37605877282253791*tmp_35 + tmp_36;
      real_t tmp_129 = tmp_128*tmp_38;
      real_t tmp_130 = tmp_125 + tmp_127 + tmp_129;
      real_t tmp_131 = tmp_124*tmp_41;
      real_t tmp_132 = tmp_126*tmp_43;
      real_t tmp_133 = tmp_128*tmp_45;
      real_t tmp_134 = tmp_131 + tmp_132 + tmp_133;
      real_t tmp_135 = tmp_124*tmp_48;
      real_t tmp_136 = tmp_126*tmp_50;
      real_t tmp_137 = tmp_128*tmp_52;
      real_t tmp_138 = tmp_135 + tmp_136 + tmp_137;
      real_t tmp_139 = tmp_1*(tmp_130 - 1.0/4.0) + tmp_11*(tmp_138 - 1.0/4.0) + tmp_8*(tmp_134 - 1.0/4.0);
      real_t tmp_140 = -tmp_125 - tmp_127 - tmp_129 - tmp_131 - tmp_132 - tmp_133 - tmp_135 - tmp_136 - tmp_137 + 1;
      real_t tmp_141 = tmp_139*tmp_63;
      real_t tmp_142 = 0.020848748529055869*tmp_65;
      real_t tmp_143 = 0.041227165399737475*tmp_3 + 0.78764240869137092*tmp_4 + tmp_6;
      real_t tmp_144 = tmp_143*tmp_24;
      real_t tmp_145 = 0.041227165399737475*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29;
      real_t tmp_146 = tmp_145*tmp_31;
      real_t tmp_147 = 0.041227165399737475*tmp_34 + 0.78764240869137092*tmp_35 + tmp_36;
      real_t tmp_148 = tmp_147*tmp_38;
      real_t tmp_149 = tmp_144 + tmp_146 + tmp_148;
      real_t tmp_150 = tmp_143*tmp_41;
      real_t tmp_151 = tmp_145*tmp_43;
      real_t tmp_152 = tmp_147*tmp_45;
      real_t tmp_153 = tmp_150 + tmp_151 + tmp_152;
      real_t tmp_154 = tmp_143*tmp_48;
      real_t tmp_155 = tmp_145*tmp_50;
      real_t tmp_156 = tmp_147*tmp_52;
      real_t tmp_157 = tmp_154 + tmp_155 + tmp_156;
      real_t tmp_158 = tmp_1*(tmp_149 - 1.0/4.0) + tmp_11*(tmp_157 - 1.0/4.0) + tmp_8*(tmp_153 - 1.0/4.0);
      real_t tmp_159 = -tmp_144 - tmp_146 - tmp_148 - tmp_150 - tmp_151 - tmp_152 - tmp_154 - tmp_155 - tmp_156 + 1;
      real_t tmp_160 = tmp_158*tmp_63;
      real_t tmp_161 = 0.019202922745021479*tmp_65;
      real_t tmp_162 = 0.039308471900058539*tmp_3 + 0.58463275527740355*tmp_4 + tmp_6;
      real_t tmp_163 = tmp_162*tmp_24;
      real_t tmp_164 = 0.039308471900058539*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29;
      real_t tmp_165 = tmp_164*tmp_31;
      real_t tmp_166 = 0.039308471900058539*tmp_34 + 0.58463275527740355*tmp_35 + tmp_36;
      real_t tmp_167 = tmp_166*tmp_38;
      real_t tmp_168 = tmp_163 + tmp_165 + tmp_167;
      real_t tmp_169 = tmp_162*tmp_41;
      real_t tmp_170 = tmp_164*tmp_43;
      real_t tmp_171 = tmp_166*tmp_45;
      real_t tmp_172 = tmp_169 + tmp_170 + tmp_171;
      real_t tmp_173 = tmp_162*tmp_48;
      real_t tmp_174 = tmp_164*tmp_50;
      real_t tmp_175 = tmp_166*tmp_52;
      real_t tmp_176 = tmp_173 + tmp_174 + tmp_175;
      real_t tmp_177 = tmp_1*(tmp_168 - 1.0/4.0) + tmp_11*(tmp_176 - 1.0/4.0) + tmp_8*(tmp_172 - 1.0/4.0);
      real_t tmp_178 = -tmp_163 - tmp_165 - tmp_167 - tmp_169 - tmp_170 - tmp_171 - tmp_173 - tmp_174 - tmp_175 + 1;
      real_t tmp_179 = tmp_177*tmp_63;
      real_t tmp_180 = 0.020848748529055869*tmp_65;
      real_t tmp_181 = 0.78764240869137092*tmp_3 + 0.041227165399737475*tmp_4 + tmp_6;
      real_t tmp_182 = tmp_181*tmp_24;
      real_t tmp_183 = 0.78764240869137092*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29;
      real_t tmp_184 = tmp_183*tmp_31;
      real_t tmp_185 = 0.78764240869137092*tmp_34 + 0.041227165399737475*tmp_35 + tmp_36;
      real_t tmp_186 = tmp_185*tmp_38;
      real_t tmp_187 = tmp_182 + tmp_184 + tmp_186;
      real_t tmp_188 = tmp_181*tmp_41;
      real_t tmp_189 = tmp_183*tmp_43;
      real_t tmp_190 = tmp_185*tmp_45;
      real_t tmp_191 = tmp_188 + tmp_189 + tmp_190;
      real_t tmp_192 = tmp_181*tmp_48;
      real_t tmp_193 = tmp_183*tmp_50;
      real_t tmp_194 = tmp_185*tmp_52;
      real_t tmp_195 = tmp_192 + tmp_193 + tmp_194;
      real_t tmp_196 = tmp_1*(tmp_187 - 1.0/4.0) + tmp_11*(tmp_195 - 1.0/4.0) + tmp_8*(tmp_191 - 1.0/4.0);
      real_t tmp_197 = -tmp_182 - tmp_184 - tmp_186 - tmp_188 - tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 + 1;
      real_t tmp_198 = tmp_196*tmp_63;
      real_t tmp_199 = 0.019202922745021479*tmp_65;
      real_t tmp_200 = 0.58463275527740355*tmp_3 + 0.039308471900058539*tmp_4 + tmp_6;
      real_t tmp_201 = tmp_200*tmp_24;
      real_t tmp_202 = 0.58463275527740355*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29;
      real_t tmp_203 = tmp_202*tmp_31;
      real_t tmp_204 = 0.58463275527740355*tmp_34 + 0.039308471900058539*tmp_35 + tmp_36;
      real_t tmp_205 = tmp_204*tmp_38;
      real_t tmp_206 = tmp_201 + tmp_203 + tmp_205;
      real_t tmp_207 = tmp_200*tmp_41;
      real_t tmp_208 = tmp_202*tmp_43;
      real_t tmp_209 = tmp_204*tmp_45;
      real_t tmp_210 = tmp_207 + tmp_208 + tmp_209;
      real_t tmp_211 = tmp_200*tmp_48;
      real_t tmp_212 = tmp_202*tmp_50;
      real_t tmp_213 = tmp_204*tmp_52;
      real_t tmp_214 = tmp_211 + tmp_212 + tmp_213;
      real_t tmp_215 = tmp_1*(tmp_206 - 1.0/4.0) + tmp_11*(tmp_214 - 1.0/4.0) + tmp_8*(tmp_210 - 1.0/4.0);
      real_t tmp_216 = -tmp_201 - tmp_203 - tmp_205 - tmp_207 - tmp_208 - tmp_209 - tmp_211 - tmp_212 - tmp_213 + 1;
      real_t tmp_217 = tmp_215*tmp_63;
      real_t tmp_218 = 0.020848748529055869*tmp_65;
      real_t tmp_219 = 0.1711304259088916*tmp_3 + 0.78764240869137092*tmp_4 + tmp_6;
      real_t tmp_220 = tmp_219*tmp_24;
      real_t tmp_221 = 0.1711304259088916*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29;
      real_t tmp_222 = tmp_221*tmp_31;
      real_t tmp_223 = 0.1711304259088916*tmp_34 + 0.78764240869137092*tmp_35 + tmp_36;
      real_t tmp_224 = tmp_223*tmp_38;
      real_t tmp_225 = tmp_220 + tmp_222 + tmp_224;
      real_t tmp_226 = tmp_219*tmp_41;
      real_t tmp_227 = tmp_221*tmp_43;
      real_t tmp_228 = tmp_223*tmp_45;
      real_t tmp_229 = tmp_226 + tmp_227 + tmp_228;
      real_t tmp_230 = tmp_219*tmp_48;
      real_t tmp_231 = tmp_221*tmp_50;
      real_t tmp_232 = tmp_223*tmp_52;
      real_t tmp_233 = tmp_230 + tmp_231 + tmp_232;
      real_t tmp_234 = tmp_1*(tmp_225 - 1.0/4.0) + tmp_11*(tmp_233 - 1.0/4.0) + tmp_8*(tmp_229 - 1.0/4.0);
      real_t tmp_235 = -tmp_220 - tmp_222 - tmp_224 - tmp_226 - tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 + 1;
      real_t tmp_236 = tmp_234*tmp_63;
      real_t tmp_237 = 0.019202922745021479*tmp_65;
      real_t tmp_238 = 0.37605877282253791*tmp_3 + 0.58463275527740355*tmp_4 + tmp_6;
      real_t tmp_239 = tmp_238*tmp_24;
      real_t tmp_240 = 0.37605877282253791*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29;
      real_t tmp_241 = tmp_240*tmp_31;
      real_t tmp_242 = 0.37605877282253791*tmp_34 + 0.58463275527740355*tmp_35 + tmp_36;
      real_t tmp_243 = tmp_242*tmp_38;
      real_t tmp_244 = tmp_239 + tmp_241 + tmp_243;
      real_t tmp_245 = tmp_238*tmp_41;
      real_t tmp_246 = tmp_240*tmp_43;
      real_t tmp_247 = tmp_242*tmp_45;
      real_t tmp_248 = tmp_245 + tmp_246 + tmp_247;
      real_t tmp_249 = tmp_238*tmp_48;
      real_t tmp_250 = tmp_240*tmp_50;
      real_t tmp_251 = tmp_242*tmp_52;
      real_t tmp_252 = tmp_249 + tmp_250 + tmp_251;
      real_t tmp_253 = tmp_1*(tmp_244 - 1.0/4.0) + tmp_11*(tmp_252 - 1.0/4.0) + tmp_8*(tmp_248 - 1.0/4.0);
      real_t tmp_254 = -tmp_239 - tmp_241 - tmp_243 - tmp_245 - tmp_246 - tmp_247 - tmp_249 - tmp_250 - tmp_251 + 1;
      real_t tmp_255 = tmp_253*tmp_63;
      real_t tmp_256 = 0.020848748529055869*tmp_65;
      real_t tmp_257 = 0.041227165399737475*tmp_3 + 0.1711304259088916*tmp_4 + tmp_6;
      real_t tmp_258 = tmp_24*tmp_257;
      real_t tmp_259 = 0.041227165399737475*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29;
      real_t tmp_260 = tmp_259*tmp_31;
      real_t tmp_261 = 0.041227165399737475*tmp_34 + 0.1711304259088916*tmp_35 + tmp_36;
      real_t tmp_262 = tmp_261*tmp_38;
      real_t tmp_263 = tmp_258 + tmp_260 + tmp_262;
      real_t tmp_264 = tmp_257*tmp_41;
      real_t tmp_265 = tmp_259*tmp_43;
      real_t tmp_266 = tmp_261*tmp_45;
      real_t tmp_267 = tmp_264 + tmp_265 + tmp_266;
      real_t tmp_268 = tmp_257*tmp_48;
      real_t tmp_269 = tmp_259*tmp_50;
      real_t tmp_270 = tmp_261*tmp_52;
      real_t tmp_271 = tmp_268 + tmp_269 + tmp_270;
      real_t tmp_272 = tmp_1*(tmp_263 - 1.0/4.0) + tmp_11*(tmp_271 - 1.0/4.0) + tmp_8*(tmp_267 - 1.0/4.0);
      real_t tmp_273 = -tmp_258 - tmp_260 - tmp_262 - tmp_264 - tmp_265 - tmp_266 - tmp_268 - tmp_269 - tmp_270 + 1;
      real_t tmp_274 = tmp_272*tmp_63;
      real_t tmp_275 = 0.019202922745021479*tmp_65;
      real_t tmp_276 = 0.40446199974765351*tmp_3 + 0.19107600050469298*tmp_4 + tmp_6;
      real_t tmp_277 = tmp_24*tmp_276;
      real_t tmp_278 = 0.40446199974765351*tmp_27 + 0.19107600050469298*tmp_28 + tmp_29;
      real_t tmp_279 = tmp_278*tmp_31;
      real_t tmp_280 = 0.40446199974765351*tmp_34 + 0.19107600050469298*tmp_35 + tmp_36;
      real_t tmp_281 = tmp_280*tmp_38;
      real_t tmp_282 = tmp_277 + tmp_279 + tmp_281;
      real_t tmp_283 = tmp_276*tmp_41;
      real_t tmp_284 = tmp_278*tmp_43;
      real_t tmp_285 = tmp_280*tmp_45;
      real_t tmp_286 = tmp_283 + tmp_284 + tmp_285;
      real_t tmp_287 = tmp_276*tmp_48;
      real_t tmp_288 = tmp_278*tmp_50;
      real_t tmp_289 = tmp_280*tmp_52;
      real_t tmp_290 = tmp_287 + tmp_288 + tmp_289;
      real_t tmp_291 = tmp_1*(tmp_282 - 1.0/4.0) + tmp_11*(tmp_290 - 1.0/4.0) + tmp_8*(tmp_286 - 1.0/4.0);
      real_t tmp_292 = -tmp_277 - tmp_279 - tmp_281 - tmp_283 - tmp_284 - tmp_285 - tmp_287 - tmp_288 - tmp_289 + 1;
      real_t tmp_293 = tmp_291*tmp_63;
      real_t tmp_294 = 0.042507265838595799*tmp_65;
      real_t tmp_295 = 0.039308471900058539*tmp_3 + 0.37605877282253791*tmp_4 + tmp_6;
      real_t tmp_296 = tmp_24*tmp_295;
      real_t tmp_297 = 0.039308471900058539*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29;
      real_t tmp_298 = tmp_297*tmp_31;
      real_t tmp_299 = 0.039308471900058539*tmp_34 + 0.37605877282253791*tmp_35 + tmp_36;
      real_t tmp_300 = tmp_299*tmp_38;
      real_t tmp_301 = tmp_296 + tmp_298 + tmp_300;
      real_t tmp_302 = tmp_295*tmp_41;
      real_t tmp_303 = tmp_297*tmp_43;
      real_t tmp_304 = tmp_299*tmp_45;
      real_t tmp_305 = tmp_302 + tmp_303 + tmp_304;
      real_t tmp_306 = tmp_295*tmp_48;
      real_t tmp_307 = tmp_297*tmp_50;
      real_t tmp_308 = tmp_299*tmp_52;
      real_t tmp_309 = tmp_306 + tmp_307 + tmp_308;
      real_t tmp_310 = tmp_1*(tmp_301 - 1.0/4.0) + tmp_11*(tmp_309 - 1.0/4.0) + tmp_8*(tmp_305 - 1.0/4.0);
      real_t tmp_311 = -tmp_296 - tmp_298 - tmp_300 - tmp_302 - tmp_303 - tmp_304 - tmp_306 - tmp_307 - tmp_308 + 1;
      real_t tmp_312 = tmp_310*tmp_63;
      real_t tmp_313 = 0.020848748529055869*tmp_65;
      real_t tmp_314 = 0.93718850182767688*tmp_3 + 0.031405749086161582*tmp_4 + tmp_6;
      real_t tmp_315 = tmp_24*tmp_314;
      real_t tmp_316 = 0.93718850182767688*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29;
      real_t tmp_317 = tmp_31*tmp_316;
      real_t tmp_318 = 0.93718850182767688*tmp_34 + 0.031405749086161582*tmp_35 + tmp_36;
      real_t tmp_319 = tmp_318*tmp_38;
      real_t tmp_320 = tmp_315 + tmp_317 + tmp_319;
      real_t tmp_321 = tmp_314*tmp_41;
      real_t tmp_322 = tmp_316*tmp_43;
      real_t tmp_323 = tmp_318*tmp_45;
      real_t tmp_324 = tmp_321 + tmp_322 + tmp_323;
      real_t tmp_325 = tmp_314*tmp_48;
      real_t tmp_326 = tmp_316*tmp_50;
      real_t tmp_327 = tmp_318*tmp_52;
      real_t tmp_328 = tmp_325 + tmp_326 + tmp_327;
      real_t tmp_329 = tmp_1*(tmp_320 - 1.0/4.0) + tmp_11*(tmp_328 - 1.0/4.0) + tmp_8*(tmp_324 - 1.0/4.0);
      real_t tmp_330 = -tmp_315 - tmp_317 - tmp_319 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1;
      real_t tmp_331 = tmp_329*tmp_63;
      real_t tmp_332 = 0.0068572537431980923*tmp_65;
      real_t tmp_333 = 0.60796128279561268*tmp_3 + 0.19601935860219369*tmp_4 + tmp_6;
      real_t tmp_334 = tmp_24*tmp_333;
      real_t tmp_335 = 0.60796128279561268*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29;
      real_t tmp_336 = tmp_31*tmp_335;
      real_t tmp_337 = 0.60796128279561268*tmp_34 + 0.19601935860219369*tmp_35 + tmp_36;
      real_t tmp_338 = tmp_337*tmp_38;
      real_t tmp_339 = tmp_334 + tmp_336 + tmp_338;
      real_t tmp_340 = tmp_333*tmp_41;
      real_t tmp_341 = tmp_335*tmp_43;
      real_t tmp_342 = tmp_337*tmp_45;
      real_t tmp_343 = tmp_340 + tmp_341 + tmp_342;
      real_t tmp_344 = tmp_333*tmp_48;
      real_t tmp_345 = tmp_335*tmp_50;
      real_t tmp_346 = tmp_337*tmp_52;
      real_t tmp_347 = tmp_344 + tmp_345 + tmp_346;
      real_t tmp_348 = tmp_1*(tmp_339 - 1.0/4.0) + tmp_11*(tmp_347 - 1.0/4.0) + tmp_8*(tmp_343 - 1.0/4.0);
      real_t tmp_349 = -tmp_334 - tmp_336 - tmp_338 - tmp_340 - tmp_341 - tmp_342 - tmp_344 - tmp_345 - tmp_346 + 1;
      real_t tmp_350 = tmp_348*tmp_63;
      real_t tmp_351 = 0.037198804536718075*tmp_65;
      real_t tmp_352 = 0.19107600050469298*tmp_3 + 0.40446199974765351*tmp_4 + tmp_6;
      real_t tmp_353 = tmp_24*tmp_352;
      real_t tmp_354 = 0.19107600050469298*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29;
      real_t tmp_355 = tmp_31*tmp_354;
      real_t tmp_356 = 0.19107600050469298*tmp_34 + 0.40446199974765351*tmp_35 + tmp_36;
      real_t tmp_357 = tmp_356*tmp_38;
      real_t tmp_358 = tmp_353 + tmp_355 + tmp_357;
      real_t tmp_359 = tmp_352*tmp_41;
      real_t tmp_360 = tmp_354*tmp_43;
      real_t tmp_361 = tmp_356*tmp_45;
      real_t tmp_362 = tmp_359 + tmp_360 + tmp_361;
      real_t tmp_363 = tmp_352*tmp_48;
      real_t tmp_364 = tmp_354*tmp_50;
      real_t tmp_365 = tmp_356*tmp_52;
      real_t tmp_366 = tmp_363 + tmp_364 + tmp_365;
      real_t tmp_367 = tmp_1*(tmp_358 - 1.0/4.0) + tmp_11*(tmp_366 - 1.0/4.0) + tmp_8*(tmp_362 - 1.0/4.0);
      real_t tmp_368 = -tmp_353 - tmp_355 - tmp_357 - tmp_359 - tmp_360 - tmp_361 - tmp_363 - tmp_364 - tmp_365 + 1;
      real_t tmp_369 = tmp_367*tmp_63;
      real_t tmp_370 = 0.042507265838595799*tmp_65;
      real_t tmp_371 = 0.031405749086161582*tmp_3 + 0.031405749086161582*tmp_4 + tmp_6;
      real_t tmp_372 = tmp_24*tmp_371;
      real_t tmp_373 = 0.031405749086161582*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29;
      real_t tmp_374 = tmp_31*tmp_373;
      real_t tmp_375 = 0.031405749086161582*tmp_34 + 0.031405749086161582*tmp_35 + tmp_36;
      real_t tmp_376 = tmp_375*tmp_38;
      real_t tmp_377 = tmp_372 + tmp_374 + tmp_376;
      real_t tmp_378 = tmp_371*tmp_41;
      real_t tmp_379 = tmp_373*tmp_43;
      real_t tmp_380 = tmp_375*tmp_45;
      real_t tmp_381 = tmp_378 + tmp_379 + tmp_380;
      real_t tmp_382 = tmp_371*tmp_48;
      real_t tmp_383 = tmp_373*tmp_50;
      real_t tmp_384 = tmp_375*tmp_52;
      real_t tmp_385 = tmp_382 + tmp_383 + tmp_384;
      real_t tmp_386 = tmp_1*(tmp_377 - 1.0/4.0) + tmp_11*(tmp_385 - 1.0/4.0) + tmp_8*(tmp_381 - 1.0/4.0);
      real_t tmp_387 = -tmp_372 - tmp_374 - tmp_376 - tmp_378 - tmp_379 - tmp_380 - tmp_382 - tmp_383 - tmp_384 + 1;
      real_t tmp_388 = tmp_386*tmp_63;
      real_t tmp_389 = 0.0068572537431980923*tmp_65;
      real_t tmp_390 = 0.19601935860219369*tmp_3 + 0.19601935860219369*tmp_4 + tmp_6;
      real_t tmp_391 = tmp_24*tmp_390;
      real_t tmp_392 = 0.19601935860219369*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29;
      real_t tmp_393 = tmp_31*tmp_392;
      real_t tmp_394 = 0.19601935860219369*tmp_34 + 0.19601935860219369*tmp_35 + tmp_36;
      real_t tmp_395 = tmp_38*tmp_394;
      real_t tmp_396 = tmp_391 + tmp_393 + tmp_395;
      real_t tmp_397 = tmp_390*tmp_41;
      real_t tmp_398 = tmp_392*tmp_43;
      real_t tmp_399 = tmp_394*tmp_45;
      real_t tmp_400 = tmp_397 + tmp_398 + tmp_399;
      real_t tmp_401 = tmp_390*tmp_48;
      real_t tmp_402 = tmp_392*tmp_50;
      real_t tmp_403 = tmp_394*tmp_52;
      real_t tmp_404 = tmp_401 + tmp_402 + tmp_403;
      real_t tmp_405 = tmp_1*(tmp_396 - 1.0/4.0) + tmp_11*(tmp_404 - 1.0/4.0) + tmp_8*(tmp_400 - 1.0/4.0);
      real_t tmp_406 = -tmp_391 - tmp_393 - tmp_395 - tmp_397 - tmp_398 - tmp_399 - tmp_401 - tmp_402 - tmp_403 + 1;
      real_t tmp_407 = tmp_405*tmp_63;
      real_t tmp_408 = 0.037198804536718075*tmp_65;
      real_t tmp_409 = 0.40446199974765351*tmp_3 + 0.40446199974765351*tmp_4 + tmp_6;
      real_t tmp_410 = tmp_24*tmp_409;
      real_t tmp_411 = 0.40446199974765351*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29;
      real_t tmp_412 = tmp_31*tmp_411;
      real_t tmp_413 = 0.40446199974765351*tmp_34 + 0.40446199974765351*tmp_35 + tmp_36;
      real_t tmp_414 = tmp_38*tmp_413;
      real_t tmp_415 = tmp_410 + tmp_412 + tmp_414;
      real_t tmp_416 = tmp_409*tmp_41;
      real_t tmp_417 = tmp_411*tmp_43;
      real_t tmp_418 = tmp_413*tmp_45;
      real_t tmp_419 = tmp_416 + tmp_417 + tmp_418;
      real_t tmp_420 = tmp_409*tmp_48;
      real_t tmp_421 = tmp_411*tmp_50;
      real_t tmp_422 = tmp_413*tmp_52;
      real_t tmp_423 = tmp_420 + tmp_421 + tmp_422;
      real_t tmp_424 = tmp_1*(tmp_415 - 1.0/4.0) + tmp_11*(tmp_423 - 1.0/4.0) + tmp_8*(tmp_419 - 1.0/4.0);
      real_t tmp_425 = -tmp_410 - tmp_412 - tmp_414 - tmp_416 - tmp_417 - tmp_418 - tmp_420 - tmp_421 - tmp_422 + 1;
      real_t tmp_426 = tmp_424*tmp_63;
      real_t tmp_427 = 0.042507265838595799*tmp_65;
      real_t tmp_428 = 0.1711304259088916*tmp_3 + 0.041227165399737475*tmp_4 + tmp_6;
      real_t tmp_429 = tmp_24*tmp_428;
      real_t tmp_430 = 0.1711304259088916*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29;
      real_t tmp_431 = tmp_31*tmp_430;
      real_t tmp_432 = 0.1711304259088916*tmp_34 + 0.041227165399737475*tmp_35 + tmp_36;
      real_t tmp_433 = tmp_38*tmp_432;
      real_t tmp_434 = tmp_429 + tmp_431 + tmp_433;
      real_t tmp_435 = tmp_41*tmp_428;
      real_t tmp_436 = tmp_43*tmp_430;
      real_t tmp_437 = tmp_432*tmp_45;
      real_t tmp_438 = tmp_435 + tmp_436 + tmp_437;
      real_t tmp_439 = tmp_428*tmp_48;
      real_t tmp_440 = tmp_430*tmp_50;
      real_t tmp_441 = tmp_432*tmp_52;
      real_t tmp_442 = tmp_439 + tmp_440 + tmp_441;
      real_t tmp_443 = tmp_1*(tmp_434 - 1.0/4.0) + tmp_11*(tmp_442 - 1.0/4.0) + tmp_8*(tmp_438 - 1.0/4.0);
      real_t tmp_444 = -tmp_429 - tmp_431 - tmp_433 - tmp_435 - tmp_436 - tmp_437 - tmp_439 - tmp_440 - tmp_441 + 1;
      real_t tmp_445 = tmp_443*tmp_63;
      real_t tmp_446 = 0.019202922745021479*tmp_65;
      real_t tmp_447 = 0.5*p_affine_13_0*tmp_38 + 0.5*p_affine_13_1*tmp_31 + 0.5*p_affine_13_2*tmp_24;
      real_t tmp_448 = 0.5*p_affine_13_0*tmp_45 + 0.5*p_affine_13_1*tmp_43 + 0.5*p_affine_13_2*tmp_41;
      real_t tmp_449 = 0.5*p_affine_13_0*tmp_52 + 0.5*p_affine_13_1*tmp_50 + 0.5*p_affine_13_2*tmp_48;
      real_t a_0_0 = tmp_104*(-tmp_101*tmp_56 + tmp_102*tmp_103 - tmp_102*tmp_58) + tmp_123*(-tmp_120*tmp_56 + tmp_121*tmp_122 - tmp_121*tmp_58) + tmp_142*(-tmp_139*tmp_56 + tmp_140*tmp_141 - tmp_140*tmp_58) + tmp_161*(-tmp_158*tmp_56 + tmp_159*tmp_160 - tmp_159*tmp_58) + tmp_180*(-tmp_177*tmp_56 + tmp_178*tmp_179 - tmp_178*tmp_58) + tmp_199*(-tmp_196*tmp_56 + tmp_197*tmp_198 - tmp_197*tmp_58) + tmp_218*(-tmp_215*tmp_56 + tmp_216*tmp_217 - tmp_216*tmp_58) + tmp_237*(-tmp_234*tmp_56 + tmp_235*tmp_236 - tmp_235*tmp_58) + tmp_256*(-tmp_253*tmp_56 + tmp_254*tmp_255 - tmp_254*tmp_58) + tmp_275*(-tmp_272*tmp_56 + tmp_273*tmp_274 - tmp_273*tmp_58) + tmp_294*(-tmp_291*tmp_56 + tmp_292*tmp_293 - tmp_292*tmp_58) + tmp_313*(-tmp_310*tmp_56 + tmp_311*tmp_312 - tmp_311*tmp_58) + tmp_332*(-tmp_329*tmp_56 + tmp_330*tmp_331 - tmp_330*tmp_58) + tmp_351*(-tmp_348*tmp_56 + tmp_349*tmp_350 - tmp_349*tmp_58) + tmp_370*(-tmp_367*tmp_56 + tmp_368*tmp_369 - tmp_368*tmp_58) + tmp_389*(-tmp_386*tmp_56 + tmp_387*tmp_388 - tmp_387*tmp_58) + tmp_408*(-tmp_405*tmp_56 + tmp_406*tmp_407 - tmp_406*tmp_58) + tmp_427*(-tmp_424*tmp_56 + tmp_425*tmp_426 - tmp_425*tmp_58) + tmp_446*(-tmp_443*tmp_56 + tmp_444*tmp_445 - tmp_444*tmp_58) + tmp_66*(-tmp_55*tmp_56 - tmp_57*tmp_58 + tmp_57*tmp_64) + tmp_85*(-tmp_56*tmp_82 - tmp_58*tmp_83 + tmp_83*tmp_84);
      real_t a_1_0 = tmp_104*(-tmp_101*tmp_447 + tmp_103*tmp_92 - tmp_58*tmp_92) + tmp_123*(tmp_111*tmp_122 - tmp_111*tmp_58 - tmp_120*tmp_447) + tmp_142*(tmp_130*tmp_141 - tmp_130*tmp_58 - tmp_139*tmp_447) + tmp_161*(tmp_149*tmp_160 - tmp_149*tmp_58 - tmp_158*tmp_447) + tmp_180*(tmp_168*tmp_179 - tmp_168*tmp_58 - tmp_177*tmp_447) + tmp_199*(tmp_187*tmp_198 - tmp_187*tmp_58 - tmp_196*tmp_447) + tmp_218*(tmp_206*tmp_217 - tmp_206*tmp_58 - tmp_215*tmp_447) + tmp_237*(tmp_225*tmp_236 - tmp_225*tmp_58 - tmp_234*tmp_447) + tmp_256*(tmp_244*tmp_255 - tmp_244*tmp_58 - tmp_253*tmp_447) + tmp_275*(tmp_263*tmp_274 - tmp_263*tmp_58 - tmp_272*tmp_447) + tmp_294*(tmp_282*tmp_293 - tmp_282*tmp_58 - tmp_291*tmp_447) + tmp_313*(tmp_301*tmp_312 - tmp_301*tmp_58 - tmp_310*tmp_447) + tmp_332*(tmp_320*tmp_331 - tmp_320*tmp_58 - tmp_329*tmp_447) + tmp_351*(tmp_339*tmp_350 - tmp_339*tmp_58 - tmp_348*tmp_447) + tmp_370*(tmp_358*tmp_369 - tmp_358*tmp_58 - tmp_367*tmp_447) + tmp_389*(tmp_377*tmp_388 - tmp_377*tmp_58 - tmp_386*tmp_447) + tmp_408*(tmp_396*tmp_407 - tmp_396*tmp_58 - tmp_405*tmp_447) + tmp_427*(tmp_415*tmp_426 - tmp_415*tmp_58 - tmp_424*tmp_447) + tmp_446*(tmp_434*tmp_445 - tmp_434*tmp_58 - tmp_443*tmp_447) + tmp_66*(-tmp_40*tmp_58 + tmp_40*tmp_64 - tmp_447*tmp_55) + tmp_85*(-tmp_447*tmp_82 - tmp_58*tmp_73 + tmp_73*tmp_84);
      real_t a_2_0 = tmp_104*(-tmp_101*tmp_448 + tmp_103*tmp_96 - tmp_58*tmp_96) + tmp_123*(tmp_115*tmp_122 - tmp_115*tmp_58 - tmp_120*tmp_448) + tmp_142*(tmp_134*tmp_141 - tmp_134*tmp_58 - tmp_139*tmp_448) + tmp_161*(tmp_153*tmp_160 - tmp_153*tmp_58 - tmp_158*tmp_448) + tmp_180*(tmp_172*tmp_179 - tmp_172*tmp_58 - tmp_177*tmp_448) + tmp_199*(tmp_191*tmp_198 - tmp_191*tmp_58 - tmp_196*tmp_448) + tmp_218*(tmp_210*tmp_217 - tmp_210*tmp_58 - tmp_215*tmp_448) + tmp_237*(tmp_229*tmp_236 - tmp_229*tmp_58 - tmp_234*tmp_448) + tmp_256*(tmp_248*tmp_255 - tmp_248*tmp_58 - tmp_253*tmp_448) + tmp_275*(tmp_267*tmp_274 - tmp_267*tmp_58 - tmp_272*tmp_448) + tmp_294*(tmp_286*tmp_293 - tmp_286*tmp_58 - tmp_291*tmp_448) + tmp_313*(tmp_305*tmp_312 - tmp_305*tmp_58 - tmp_310*tmp_448) + tmp_332*(tmp_324*tmp_331 - tmp_324*tmp_58 - tmp_329*tmp_448) + tmp_351*(tmp_343*tmp_350 - tmp_343*tmp_58 - tmp_348*tmp_448) + tmp_370*(tmp_362*tmp_369 - tmp_362*tmp_58 - tmp_367*tmp_448) + tmp_389*(tmp_381*tmp_388 - tmp_381*tmp_58 - tmp_386*tmp_448) + tmp_408*(tmp_400*tmp_407 - tmp_400*tmp_58 - tmp_405*tmp_448) + tmp_427*(tmp_419*tmp_426 - tmp_419*tmp_58 - tmp_424*tmp_448) + tmp_446*(tmp_438*tmp_445 - tmp_438*tmp_58 - tmp_443*tmp_448) + tmp_66*(-tmp_448*tmp_55 - tmp_47*tmp_58 + tmp_47*tmp_64) + tmp_85*(-tmp_448*tmp_82 - tmp_58*tmp_77 + tmp_77*tmp_84);
      real_t a_3_0 = tmp_104*(tmp_100*tmp_103 - tmp_100*tmp_58 - tmp_101*tmp_449) + tmp_123*(tmp_119*tmp_122 - tmp_119*tmp_58 - tmp_120*tmp_449) + tmp_142*(tmp_138*tmp_141 - tmp_138*tmp_58 - tmp_139*tmp_449) + tmp_161*(tmp_157*tmp_160 - tmp_157*tmp_58 - tmp_158*tmp_449) + tmp_180*(tmp_176*tmp_179 - tmp_176*tmp_58 - tmp_177*tmp_449) + tmp_199*(tmp_195*tmp_198 - tmp_195*tmp_58 - tmp_196*tmp_449) + tmp_218*(tmp_214*tmp_217 - tmp_214*tmp_58 - tmp_215*tmp_449) + tmp_237*(tmp_233*tmp_236 - tmp_233*tmp_58 - tmp_234*tmp_449) + tmp_256*(tmp_252*tmp_255 - tmp_252*tmp_58 - tmp_253*tmp_449) + tmp_275*(tmp_271*tmp_274 - tmp_271*tmp_58 - tmp_272*tmp_449) + tmp_294*(tmp_290*tmp_293 - tmp_290*tmp_58 - tmp_291*tmp_449) + tmp_313*(tmp_309*tmp_312 - tmp_309*tmp_58 - tmp_310*tmp_449) + tmp_332*(tmp_328*tmp_331 - tmp_328*tmp_58 - tmp_329*tmp_449) + tmp_351*(tmp_347*tmp_350 - tmp_347*tmp_58 - tmp_348*tmp_449) + tmp_370*(tmp_366*tmp_369 - tmp_366*tmp_58 - tmp_367*tmp_449) + tmp_389*(tmp_385*tmp_388 - tmp_385*tmp_58 - tmp_386*tmp_449) + tmp_408*(tmp_404*tmp_407 - tmp_404*tmp_58 - tmp_405*tmp_449) + tmp_427*(tmp_423*tmp_426 - tmp_423*tmp_58 - tmp_424*tmp_449) + tmp_446*(tmp_442*tmp_445 - tmp_442*tmp_58 - tmp_443*tmp_449) + tmp_66*(-tmp_449*tmp_55 - tmp_54*tmp_58 + tmp_54*tmp_64) + tmp_85*(-tmp_449*tmp_82 - tmp_58*tmp_81 + tmp_81*tmp_84);
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


      real_t tmp_0 = -p_affine_4_0;
      real_t tmp_1 = p_affine_5_0 + tmp_0;
      real_t tmp_2 = p_affine_6_0 + tmp_0;
      real_t tmp_3 = -p_affine_4_1;
      real_t tmp_4 = p_affine_7_1 + tmp_3;
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = p_affine_7_0 + tmp_0;
      real_t tmp_7 = p_affine_6_1 + tmp_3;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = tmp_5 - tmp_8;
      real_t tmp_10 = -p_affine_4_2;
      real_t tmp_11 = p_affine_7_2 + tmp_10;
      real_t tmp_12 = tmp_11*tmp_7;
      real_t tmp_13 = p_affine_5_2 + tmp_10;
      real_t tmp_14 = p_affine_5_1 + tmp_3;
      real_t tmp_15 = p_affine_6_2 + tmp_10;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_15*tmp_4;
      real_t tmp_18 = tmp_11*tmp_2;
      real_t tmp_19 = 1.0 / (tmp_1*tmp_12 - tmp_1*tmp_17 + tmp_13*tmp_5 - tmp_13*tmp_8 + tmp_14*tmp_16 - tmp_14*tmp_18);
      real_t tmp_20 = p_affine_8_2 + tmp_10;
      real_t tmp_21 = -p_affine_8_2;
      real_t tmp_22 = p_affine_9_2 + tmp_21;
      real_t tmp_23 = p_affine_10_2 + tmp_21;
      real_t tmp_24 = 0.031405749086161582*tmp_22 + 0.93718850182767688*tmp_23;
      real_t tmp_25 = tmp_19*(tmp_20 + tmp_24);
      real_t tmp_26 = tmp_16 - tmp_18;
      real_t tmp_27 = p_affine_8_1 + tmp_3;
      real_t tmp_28 = -p_affine_8_1;
      real_t tmp_29 = p_affine_9_1 + tmp_28;
      real_t tmp_30 = p_affine_10_1 + tmp_28;
      real_t tmp_31 = 0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_30;
      real_t tmp_32 = tmp_19*(tmp_27 + tmp_31);
      real_t tmp_33 = tmp_12 - tmp_17;
      real_t tmp_34 = p_affine_8_0 + tmp_0;
      real_t tmp_35 = -p_affine_8_0;
      real_t tmp_36 = p_affine_9_0 + tmp_35;
      real_t tmp_37 = p_affine_10_0 + tmp_35;
      real_t tmp_38 = 0.031405749086161582*tmp_36 + 0.93718850182767688*tmp_37;
      real_t tmp_39 = tmp_19*(tmp_34 + tmp_38);
      real_t tmp_40 = -tmp_1*tmp_4 + tmp_14*tmp_6;
      real_t tmp_41 = tmp_1*tmp_11 - tmp_13*tmp_6;
      real_t tmp_42 = -tmp_11*tmp_14 + tmp_13*tmp_4;
      real_t tmp_43 = tmp_1*tmp_7 - tmp_14*tmp_2;
      real_t tmp_44 = -tmp_1*tmp_15 + tmp_13*tmp_2;
      real_t tmp_45 = -tmp_13*tmp_7 + tmp_14*tmp_15;
      real_t tmp_46 = tmp_1*(tmp_25*tmp_9 + tmp_26*tmp_32 + tmp_33*tmp_39 - 1.0/4.0) + tmp_2*(tmp_25*tmp_40 + tmp_32*tmp_41 + tmp_39*tmp_42 - 1.0/4.0) + tmp_6*(tmp_25*tmp_43 + tmp_32*tmp_44 + tmp_39*tmp_45 - 1.0/4.0);
      real_t tmp_47 = -p_affine_0_1;
      real_t tmp_48 = p_affine_1_1 + tmp_47;
      real_t tmp_49 = -p_affine_0_2;
      real_t tmp_50 = p_affine_2_2 + tmp_49;
      real_t tmp_51 = tmp_48*tmp_50;
      real_t tmp_52 = p_affine_2_1 + tmp_47;
      real_t tmp_53 = p_affine_1_2 + tmp_49;
      real_t tmp_54 = tmp_52*tmp_53;
      real_t tmp_55 = -p_affine_0_0;
      real_t tmp_56 = p_affine_1_0 + tmp_55;
      real_t tmp_57 = p_affine_3_2 + tmp_49;
      real_t tmp_58 = tmp_52*tmp_57;
      real_t tmp_59 = p_affine_2_0 + tmp_55;
      real_t tmp_60 = p_affine_3_1 + tmp_47;
      real_t tmp_61 = tmp_53*tmp_60;
      real_t tmp_62 = p_affine_3_0 + tmp_55;
      real_t tmp_63 = tmp_50*tmp_60;
      real_t tmp_64 = tmp_48*tmp_57;
      real_t tmp_65 = 1.0 / (tmp_51*tmp_62 - tmp_54*tmp_62 + tmp_56*tmp_58 - tmp_56*tmp_63 + tmp_59*tmp_61 - tmp_59*tmp_64);
      real_t tmp_66 = tmp_65*(tmp_51 - tmp_54);
      real_t tmp_67 = tmp_65*(tmp_61 - tmp_64);
      real_t tmp_68 = tmp_65*(tmp_58 - tmp_63);
      real_t tmp_69 = tmp_65*(-tmp_50*tmp_56 + tmp_53*tmp_59);
      real_t tmp_70 = tmp_65*(-tmp_53*tmp_62 + tmp_56*tmp_57);
      real_t tmp_71 = tmp_65*(tmp_50*tmp_62 - tmp_57*tmp_59);
      real_t tmp_72 = tmp_65*(-tmp_48*tmp_59 + tmp_52*tmp_56);
      real_t tmp_73 = tmp_65*(tmp_48*tmp_62 - tmp_56*tmp_60);
      real_t tmp_74 = tmp_65*(-tmp_52*tmp_62 + tmp_59*tmp_60);
      real_t tmp_75 = 0.5*p_affine_13_0*(-tmp_66 - tmp_67 - tmp_68) + 0.5*p_affine_13_1*(-tmp_69 - tmp_70 - tmp_71) + 0.5*p_affine_13_2*(-tmp_72 - tmp_73 - tmp_74);
      real_t tmp_76 = p_affine_8_2 + tmp_49;
      real_t tmp_77 = tmp_24 + tmp_76;
      real_t tmp_78 = tmp_72*tmp_77;
      real_t tmp_79 = tmp_73*tmp_77;
      real_t tmp_80 = p_affine_8_1 + tmp_47;
      real_t tmp_81 = tmp_31 + tmp_80;
      real_t tmp_82 = tmp_69*tmp_81;
      real_t tmp_83 = tmp_70*tmp_81;
      real_t tmp_84 = tmp_74*tmp_77;
      real_t tmp_85 = tmp_71*tmp_81;
      real_t tmp_86 = p_affine_8_0 + tmp_55;
      real_t tmp_87 = tmp_38 + tmp_86;
      real_t tmp_88 = tmp_66*tmp_87;
      real_t tmp_89 = tmp_67*tmp_87;
      real_t tmp_90 = tmp_68*tmp_87;
      real_t tmp_91 = -tmp_78 - tmp_79 - tmp_82 - tmp_83 - tmp_84 - tmp_85 - tmp_88 - tmp_89 - tmp_90 + 1;
      real_t tmp_92 = tmp_1*tmp_19;
      real_t tmp_93 = tmp_19*tmp_2;
      real_t tmp_94 = tmp_19*tmp_6;
      real_t tmp_95 = 0.5*p_affine_13_0*(tmp_33*tmp_92 + tmp_42*tmp_93 + tmp_45*tmp_94) + 0.5*p_affine_13_1*(tmp_26*tmp_92 + tmp_41*tmp_93 + tmp_44*tmp_94) + 0.5*p_affine_13_2*(tmp_40*tmp_93 + tmp_43*tmp_94 + tmp_9*tmp_92);
      real_t tmp_96 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_97 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_98 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_99 = (std::abs(tmp_23*tmp_96 - tmp_30*tmp_98)*std::abs(tmp_23*tmp_96 - tmp_30*tmp_98)) + (std::abs(tmp_23*tmp_97 - tmp_37*tmp_98)*std::abs(tmp_23*tmp_97 - tmp_37*tmp_98)) + (std::abs(tmp_30*tmp_97 - tmp_37*tmp_96)*std::abs(tmp_30*tmp_97 - tmp_37*tmp_96));
      real_t tmp_100 = 5.0*std::pow(tmp_99, -0.25);
      real_t tmp_101 = tmp_100*tmp_46;
      real_t tmp_102 = 1.0*std::pow(tmp_99, 1.0/2.0);
      real_t tmp_103 = 0.0068572537431980923*tmp_102;
      real_t tmp_104 = 0.19601935860219369*tmp_22 + 0.60796128279561268*tmp_23;
      real_t tmp_105 = tmp_19*(tmp_104 + tmp_20);
      real_t tmp_106 = 0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_30;
      real_t tmp_107 = tmp_19*(tmp_106 + tmp_27);
      real_t tmp_108 = 0.19601935860219369*tmp_36 + 0.60796128279561268*tmp_37;
      real_t tmp_109 = tmp_19*(tmp_108 + tmp_34);
      real_t tmp_110 = tmp_1*(tmp_105*tmp_9 + tmp_107*tmp_26 + tmp_109*tmp_33 - 1.0/4.0) + tmp_2*(tmp_105*tmp_40 + tmp_107*tmp_41 + tmp_109*tmp_42 - 1.0/4.0) + tmp_6*(tmp_105*tmp_43 + tmp_107*tmp_44 + tmp_109*tmp_45 - 1.0/4.0);
      real_t tmp_111 = tmp_104 + tmp_76;
      real_t tmp_112 = tmp_111*tmp_72;
      real_t tmp_113 = tmp_111*tmp_73;
      real_t tmp_114 = tmp_106 + tmp_80;
      real_t tmp_115 = tmp_114*tmp_69;
      real_t tmp_116 = tmp_114*tmp_70;
      real_t tmp_117 = tmp_111*tmp_74;
      real_t tmp_118 = tmp_114*tmp_71;
      real_t tmp_119 = tmp_108 + tmp_86;
      real_t tmp_120 = tmp_119*tmp_66;
      real_t tmp_121 = tmp_119*tmp_67;
      real_t tmp_122 = tmp_119*tmp_68;
      real_t tmp_123 = -tmp_112 - tmp_113 - tmp_115 - tmp_116 - tmp_117 - tmp_118 - tmp_120 - tmp_121 - tmp_122 + 1;
      real_t tmp_124 = tmp_100*tmp_110;
      real_t tmp_125 = 0.037198804536718075*tmp_102;
      real_t tmp_126 = 0.37605877282253791*tmp_22 + 0.039308471900058539*tmp_23;
      real_t tmp_127 = tmp_19*(tmp_126 + tmp_20);
      real_t tmp_128 = 0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_30;
      real_t tmp_129 = tmp_19*(tmp_128 + tmp_27);
      real_t tmp_130 = 0.37605877282253791*tmp_36 + 0.039308471900058539*tmp_37;
      real_t tmp_131 = tmp_19*(tmp_130 + tmp_34);
      real_t tmp_132 = tmp_1*(tmp_127*tmp_9 + tmp_129*tmp_26 + tmp_131*tmp_33 - 1.0/4.0) + tmp_2*(tmp_127*tmp_40 + tmp_129*tmp_41 + tmp_131*tmp_42 - 1.0/4.0) + tmp_6*(tmp_127*tmp_43 + tmp_129*tmp_44 + tmp_131*tmp_45 - 1.0/4.0);
      real_t tmp_133 = tmp_126 + tmp_76;
      real_t tmp_134 = tmp_133*tmp_72;
      real_t tmp_135 = tmp_133*tmp_73;
      real_t tmp_136 = tmp_128 + tmp_80;
      real_t tmp_137 = tmp_136*tmp_69;
      real_t tmp_138 = tmp_136*tmp_70;
      real_t tmp_139 = tmp_133*tmp_74;
      real_t tmp_140 = tmp_136*tmp_71;
      real_t tmp_141 = tmp_130 + tmp_86;
      real_t tmp_142 = tmp_141*tmp_66;
      real_t tmp_143 = tmp_141*tmp_67;
      real_t tmp_144 = tmp_141*tmp_68;
      real_t tmp_145 = -tmp_134 - tmp_135 - tmp_137 - tmp_138 - tmp_139 - tmp_140 - tmp_142 - tmp_143 - tmp_144 + 1;
      real_t tmp_146 = tmp_100*tmp_132;
      real_t tmp_147 = 0.020848748529055869*tmp_102;
      real_t tmp_148 = 0.78764240869137092*tmp_22 + 0.1711304259088916*tmp_23;
      real_t tmp_149 = tmp_19*(tmp_148 + tmp_20);
      real_t tmp_150 = 0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_30;
      real_t tmp_151 = tmp_19*(tmp_150 + tmp_27);
      real_t tmp_152 = 0.78764240869137092*tmp_36 + 0.1711304259088916*tmp_37;
      real_t tmp_153 = tmp_19*(tmp_152 + tmp_34);
      real_t tmp_154 = tmp_1*(tmp_149*tmp_9 + tmp_151*tmp_26 + tmp_153*tmp_33 - 1.0/4.0) + tmp_2*(tmp_149*tmp_40 + tmp_151*tmp_41 + tmp_153*tmp_42 - 1.0/4.0) + tmp_6*(tmp_149*tmp_43 + tmp_151*tmp_44 + tmp_153*tmp_45 - 1.0/4.0);
      real_t tmp_155 = tmp_148 + tmp_76;
      real_t tmp_156 = tmp_155*tmp_72;
      real_t tmp_157 = tmp_155*tmp_73;
      real_t tmp_158 = tmp_150 + tmp_80;
      real_t tmp_159 = tmp_158*tmp_69;
      real_t tmp_160 = tmp_158*tmp_70;
      real_t tmp_161 = tmp_155*tmp_74;
      real_t tmp_162 = tmp_158*tmp_71;
      real_t tmp_163 = tmp_152 + tmp_86;
      real_t tmp_164 = tmp_163*tmp_66;
      real_t tmp_165 = tmp_163*tmp_67;
      real_t tmp_166 = tmp_163*tmp_68;
      real_t tmp_167 = -tmp_156 - tmp_157 - tmp_159 - tmp_160 - tmp_161 - tmp_162 - tmp_164 - tmp_165 - tmp_166 + 1;
      real_t tmp_168 = tmp_100*tmp_154;
      real_t tmp_169 = 0.019202922745021479*tmp_102;
      real_t tmp_170 = 0.58463275527740355*tmp_22 + 0.37605877282253791*tmp_23;
      real_t tmp_171 = tmp_19*(tmp_170 + tmp_20);
      real_t tmp_172 = 0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_30;
      real_t tmp_173 = tmp_19*(tmp_172 + tmp_27);
      real_t tmp_174 = 0.58463275527740355*tmp_36 + 0.37605877282253791*tmp_37;
      real_t tmp_175 = tmp_19*(tmp_174 + tmp_34);
      real_t tmp_176 = tmp_1*(tmp_171*tmp_9 + tmp_173*tmp_26 + tmp_175*tmp_33 - 1.0/4.0) + tmp_2*(tmp_171*tmp_40 + tmp_173*tmp_41 + tmp_175*tmp_42 - 1.0/4.0) + tmp_6*(tmp_171*tmp_43 + tmp_173*tmp_44 + tmp_175*tmp_45 - 1.0/4.0);
      real_t tmp_177 = tmp_170 + tmp_76;
      real_t tmp_178 = tmp_177*tmp_72;
      real_t tmp_179 = tmp_177*tmp_73;
      real_t tmp_180 = tmp_172 + tmp_80;
      real_t tmp_181 = tmp_180*tmp_69;
      real_t tmp_182 = tmp_180*tmp_70;
      real_t tmp_183 = tmp_177*tmp_74;
      real_t tmp_184 = tmp_180*tmp_71;
      real_t tmp_185 = tmp_174 + tmp_86;
      real_t tmp_186 = tmp_185*tmp_66;
      real_t tmp_187 = tmp_185*tmp_67;
      real_t tmp_188 = tmp_185*tmp_68;
      real_t tmp_189 = -tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 - tmp_184 - tmp_186 - tmp_187 - tmp_188 + 1;
      real_t tmp_190 = tmp_100*tmp_176;
      real_t tmp_191 = 0.020848748529055869*tmp_102;
      real_t tmp_192 = 0.041227165399737475*tmp_22 + 0.78764240869137092*tmp_23;
      real_t tmp_193 = tmp_19*(tmp_192 + tmp_20);
      real_t tmp_194 = 0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_30;
      real_t tmp_195 = tmp_19*(tmp_194 + tmp_27);
      real_t tmp_196 = 0.041227165399737475*tmp_36 + 0.78764240869137092*tmp_37;
      real_t tmp_197 = tmp_19*(tmp_196 + tmp_34);
      real_t tmp_198 = tmp_1*(tmp_193*tmp_9 + tmp_195*tmp_26 + tmp_197*tmp_33 - 1.0/4.0) + tmp_2*(tmp_193*tmp_40 + tmp_195*tmp_41 + tmp_197*tmp_42 - 1.0/4.0) + tmp_6*(tmp_193*tmp_43 + tmp_195*tmp_44 + tmp_197*tmp_45 - 1.0/4.0);
      real_t tmp_199 = tmp_192 + tmp_76;
      real_t tmp_200 = tmp_199*tmp_72;
      real_t tmp_201 = tmp_199*tmp_73;
      real_t tmp_202 = tmp_194 + tmp_80;
      real_t tmp_203 = tmp_202*tmp_69;
      real_t tmp_204 = tmp_202*tmp_70;
      real_t tmp_205 = tmp_199*tmp_74;
      real_t tmp_206 = tmp_202*tmp_71;
      real_t tmp_207 = tmp_196 + tmp_86;
      real_t tmp_208 = tmp_207*tmp_66;
      real_t tmp_209 = tmp_207*tmp_67;
      real_t tmp_210 = tmp_207*tmp_68;
      real_t tmp_211 = -tmp_200 - tmp_201 - tmp_203 - tmp_204 - tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 + 1;
      real_t tmp_212 = tmp_100*tmp_198;
      real_t tmp_213 = 0.019202922745021479*tmp_102;
      real_t tmp_214 = 0.039308471900058539*tmp_22 + 0.58463275527740355*tmp_23;
      real_t tmp_215 = tmp_19*(tmp_20 + tmp_214);
      real_t tmp_216 = 0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_30;
      real_t tmp_217 = tmp_19*(tmp_216 + tmp_27);
      real_t tmp_218 = 0.039308471900058539*tmp_36 + 0.58463275527740355*tmp_37;
      real_t tmp_219 = tmp_19*(tmp_218 + tmp_34);
      real_t tmp_220 = tmp_1*(tmp_215*tmp_9 + tmp_217*tmp_26 + tmp_219*tmp_33 - 1.0/4.0) + tmp_2*(tmp_215*tmp_40 + tmp_217*tmp_41 + tmp_219*tmp_42 - 1.0/4.0) + tmp_6*(tmp_215*tmp_43 + tmp_217*tmp_44 + tmp_219*tmp_45 - 1.0/4.0);
      real_t tmp_221 = tmp_214 + tmp_76;
      real_t tmp_222 = tmp_221*tmp_72;
      real_t tmp_223 = tmp_221*tmp_73;
      real_t tmp_224 = tmp_216 + tmp_80;
      real_t tmp_225 = tmp_224*tmp_69;
      real_t tmp_226 = tmp_224*tmp_70;
      real_t tmp_227 = tmp_221*tmp_74;
      real_t tmp_228 = tmp_224*tmp_71;
      real_t tmp_229 = tmp_218 + tmp_86;
      real_t tmp_230 = tmp_229*tmp_66;
      real_t tmp_231 = tmp_229*tmp_67;
      real_t tmp_232 = tmp_229*tmp_68;
      real_t tmp_233 = -tmp_222 - tmp_223 - tmp_225 - tmp_226 - tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 + 1;
      real_t tmp_234 = tmp_100*tmp_220;
      real_t tmp_235 = 0.020848748529055869*tmp_102;
      real_t tmp_236 = 0.78764240869137092*tmp_22 + 0.041227165399737475*tmp_23;
      real_t tmp_237 = tmp_19*(tmp_20 + tmp_236);
      real_t tmp_238 = 0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_30;
      real_t tmp_239 = tmp_19*(tmp_238 + tmp_27);
      real_t tmp_240 = 0.78764240869137092*tmp_36 + 0.041227165399737475*tmp_37;
      real_t tmp_241 = tmp_19*(tmp_240 + tmp_34);
      real_t tmp_242 = tmp_1*(tmp_237*tmp_9 + tmp_239*tmp_26 + tmp_241*tmp_33 - 1.0/4.0) + tmp_2*(tmp_237*tmp_40 + tmp_239*tmp_41 + tmp_241*tmp_42 - 1.0/4.0) + tmp_6*(tmp_237*tmp_43 + tmp_239*tmp_44 + tmp_241*tmp_45 - 1.0/4.0);
      real_t tmp_243 = tmp_236 + tmp_76;
      real_t tmp_244 = tmp_243*tmp_72;
      real_t tmp_245 = tmp_243*tmp_73;
      real_t tmp_246 = tmp_238 + tmp_80;
      real_t tmp_247 = tmp_246*tmp_69;
      real_t tmp_248 = tmp_246*tmp_70;
      real_t tmp_249 = tmp_243*tmp_74;
      real_t tmp_250 = tmp_246*tmp_71;
      real_t tmp_251 = tmp_240 + tmp_86;
      real_t tmp_252 = tmp_251*tmp_66;
      real_t tmp_253 = tmp_251*tmp_67;
      real_t tmp_254 = tmp_251*tmp_68;
      real_t tmp_255 = -tmp_244 - tmp_245 - tmp_247 - tmp_248 - tmp_249 - tmp_250 - tmp_252 - tmp_253 - tmp_254 + 1;
      real_t tmp_256 = tmp_100*tmp_242;
      real_t tmp_257 = 0.019202922745021479*tmp_102;
      real_t tmp_258 = 0.58463275527740355*tmp_22 + 0.039308471900058539*tmp_23;
      real_t tmp_259 = tmp_19*(tmp_20 + tmp_258);
      real_t tmp_260 = 0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_30;
      real_t tmp_261 = tmp_19*(tmp_260 + tmp_27);
      real_t tmp_262 = 0.58463275527740355*tmp_36 + 0.039308471900058539*tmp_37;
      real_t tmp_263 = tmp_19*(tmp_262 + tmp_34);
      real_t tmp_264 = tmp_1*(tmp_259*tmp_9 + tmp_26*tmp_261 + tmp_263*tmp_33 - 1.0/4.0) + tmp_2*(tmp_259*tmp_40 + tmp_261*tmp_41 + tmp_263*tmp_42 - 1.0/4.0) + tmp_6*(tmp_259*tmp_43 + tmp_261*tmp_44 + tmp_263*tmp_45 - 1.0/4.0);
      real_t tmp_265 = tmp_258 + tmp_76;
      real_t tmp_266 = tmp_265*tmp_72;
      real_t tmp_267 = tmp_265*tmp_73;
      real_t tmp_268 = tmp_260 + tmp_80;
      real_t tmp_269 = tmp_268*tmp_69;
      real_t tmp_270 = tmp_268*tmp_70;
      real_t tmp_271 = tmp_265*tmp_74;
      real_t tmp_272 = tmp_268*tmp_71;
      real_t tmp_273 = tmp_262 + tmp_86;
      real_t tmp_274 = tmp_273*tmp_66;
      real_t tmp_275 = tmp_273*tmp_67;
      real_t tmp_276 = tmp_273*tmp_68;
      real_t tmp_277 = -tmp_266 - tmp_267 - tmp_269 - tmp_270 - tmp_271 - tmp_272 - tmp_274 - tmp_275 - tmp_276 + 1;
      real_t tmp_278 = tmp_100*tmp_264;
      real_t tmp_279 = 0.020848748529055869*tmp_102;
      real_t tmp_280 = 0.1711304259088916*tmp_22 + 0.78764240869137092*tmp_23;
      real_t tmp_281 = tmp_19*(tmp_20 + tmp_280);
      real_t tmp_282 = 0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_30;
      real_t tmp_283 = tmp_19*(tmp_27 + tmp_282);
      real_t tmp_284 = 0.1711304259088916*tmp_36 + 0.78764240869137092*tmp_37;
      real_t tmp_285 = tmp_19*(tmp_284 + tmp_34);
      real_t tmp_286 = tmp_1*(tmp_26*tmp_283 + tmp_281*tmp_9 + tmp_285*tmp_33 - 1.0/4.0) + tmp_2*(tmp_281*tmp_40 + tmp_283*tmp_41 + tmp_285*tmp_42 - 1.0/4.0) + tmp_6*(tmp_281*tmp_43 + tmp_283*tmp_44 + tmp_285*tmp_45 - 1.0/4.0);
      real_t tmp_287 = tmp_280 + tmp_76;
      real_t tmp_288 = tmp_287*tmp_72;
      real_t tmp_289 = tmp_287*tmp_73;
      real_t tmp_290 = tmp_282 + tmp_80;
      real_t tmp_291 = tmp_290*tmp_69;
      real_t tmp_292 = tmp_290*tmp_70;
      real_t tmp_293 = tmp_287*tmp_74;
      real_t tmp_294 = tmp_290*tmp_71;
      real_t tmp_295 = tmp_284 + tmp_86;
      real_t tmp_296 = tmp_295*tmp_66;
      real_t tmp_297 = tmp_295*tmp_67;
      real_t tmp_298 = tmp_295*tmp_68;
      real_t tmp_299 = -tmp_288 - tmp_289 - tmp_291 - tmp_292 - tmp_293 - tmp_294 - tmp_296 - tmp_297 - tmp_298 + 1;
      real_t tmp_300 = tmp_100*tmp_286;
      real_t tmp_301 = 0.019202922745021479*tmp_102;
      real_t tmp_302 = 0.37605877282253791*tmp_22 + 0.58463275527740355*tmp_23;
      real_t tmp_303 = tmp_19*(tmp_20 + tmp_302);
      real_t tmp_304 = 0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_30;
      real_t tmp_305 = tmp_19*(tmp_27 + tmp_304);
      real_t tmp_306 = 0.37605877282253791*tmp_36 + 0.58463275527740355*tmp_37;
      real_t tmp_307 = tmp_19*(tmp_306 + tmp_34);
      real_t tmp_308 = tmp_1*(tmp_26*tmp_305 + tmp_303*tmp_9 + tmp_307*tmp_33 - 1.0/4.0) + tmp_2*(tmp_303*tmp_40 + tmp_305*tmp_41 + tmp_307*tmp_42 - 1.0/4.0) + tmp_6*(tmp_303*tmp_43 + tmp_305*tmp_44 + tmp_307*tmp_45 - 1.0/4.0);
      real_t tmp_309 = tmp_302 + tmp_76;
      real_t tmp_310 = tmp_309*tmp_72;
      real_t tmp_311 = tmp_309*tmp_73;
      real_t tmp_312 = tmp_304 + tmp_80;
      real_t tmp_313 = tmp_312*tmp_69;
      real_t tmp_314 = tmp_312*tmp_70;
      real_t tmp_315 = tmp_309*tmp_74;
      real_t tmp_316 = tmp_312*tmp_71;
      real_t tmp_317 = tmp_306 + tmp_86;
      real_t tmp_318 = tmp_317*tmp_66;
      real_t tmp_319 = tmp_317*tmp_67;
      real_t tmp_320 = tmp_317*tmp_68;
      real_t tmp_321 = -tmp_310 - tmp_311 - tmp_313 - tmp_314 - tmp_315 - tmp_316 - tmp_318 - tmp_319 - tmp_320 + 1;
      real_t tmp_322 = tmp_100*tmp_308;
      real_t tmp_323 = 0.020848748529055869*tmp_102;
      real_t tmp_324 = 0.041227165399737475*tmp_22 + 0.1711304259088916*tmp_23;
      real_t tmp_325 = tmp_19*(tmp_20 + tmp_324);
      real_t tmp_326 = 0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_30;
      real_t tmp_327 = tmp_19*(tmp_27 + tmp_326);
      real_t tmp_328 = 0.041227165399737475*tmp_36 + 0.1711304259088916*tmp_37;
      real_t tmp_329 = tmp_19*(tmp_328 + tmp_34);
      real_t tmp_330 = tmp_1*(tmp_26*tmp_327 + tmp_325*tmp_9 + tmp_329*tmp_33 - 1.0/4.0) + tmp_2*(tmp_325*tmp_40 + tmp_327*tmp_41 + tmp_329*tmp_42 - 1.0/4.0) + tmp_6*(tmp_325*tmp_43 + tmp_327*tmp_44 + tmp_329*tmp_45 - 1.0/4.0);
      real_t tmp_331 = tmp_324 + tmp_76;
      real_t tmp_332 = tmp_331*tmp_72;
      real_t tmp_333 = tmp_331*tmp_73;
      real_t tmp_334 = tmp_326 + tmp_80;
      real_t tmp_335 = tmp_334*tmp_69;
      real_t tmp_336 = tmp_334*tmp_70;
      real_t tmp_337 = tmp_331*tmp_74;
      real_t tmp_338 = tmp_334*tmp_71;
      real_t tmp_339 = tmp_328 + tmp_86;
      real_t tmp_340 = tmp_339*tmp_66;
      real_t tmp_341 = tmp_339*tmp_67;
      real_t tmp_342 = tmp_339*tmp_68;
      real_t tmp_343 = -tmp_332 - tmp_333 - tmp_335 - tmp_336 - tmp_337 - tmp_338 - tmp_340 - tmp_341 - tmp_342 + 1;
      real_t tmp_344 = tmp_100*tmp_330;
      real_t tmp_345 = 0.019202922745021479*tmp_102;
      real_t tmp_346 = 0.40446199974765351*tmp_22 + 0.19107600050469298*tmp_23;
      real_t tmp_347 = tmp_19*(tmp_20 + tmp_346);
      real_t tmp_348 = 0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_30;
      real_t tmp_349 = tmp_19*(tmp_27 + tmp_348);
      real_t tmp_350 = 0.40446199974765351*tmp_36 + 0.19107600050469298*tmp_37;
      real_t tmp_351 = tmp_19*(tmp_34 + tmp_350);
      real_t tmp_352 = tmp_1*(tmp_26*tmp_349 + tmp_33*tmp_351 + tmp_347*tmp_9 - 1.0/4.0) + tmp_2*(tmp_347*tmp_40 + tmp_349*tmp_41 + tmp_351*tmp_42 - 1.0/4.0) + tmp_6*(tmp_347*tmp_43 + tmp_349*tmp_44 + tmp_351*tmp_45 - 1.0/4.0);
      real_t tmp_353 = tmp_346 + tmp_76;
      real_t tmp_354 = tmp_353*tmp_72;
      real_t tmp_355 = tmp_353*tmp_73;
      real_t tmp_356 = tmp_348 + tmp_80;
      real_t tmp_357 = tmp_356*tmp_69;
      real_t tmp_358 = tmp_356*tmp_70;
      real_t tmp_359 = tmp_353*tmp_74;
      real_t tmp_360 = tmp_356*tmp_71;
      real_t tmp_361 = tmp_350 + tmp_86;
      real_t tmp_362 = tmp_361*tmp_66;
      real_t tmp_363 = tmp_361*tmp_67;
      real_t tmp_364 = tmp_361*tmp_68;
      real_t tmp_365 = -tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 - tmp_360 - tmp_362 - tmp_363 - tmp_364 + 1;
      real_t tmp_366 = tmp_100*tmp_352;
      real_t tmp_367 = 0.042507265838595799*tmp_102;
      real_t tmp_368 = 0.039308471900058539*tmp_22 + 0.37605877282253791*tmp_23;
      real_t tmp_369 = tmp_19*(tmp_20 + tmp_368);
      real_t tmp_370 = 0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_30;
      real_t tmp_371 = tmp_19*(tmp_27 + tmp_370);
      real_t tmp_372 = 0.039308471900058539*tmp_36 + 0.37605877282253791*tmp_37;
      real_t tmp_373 = tmp_19*(tmp_34 + tmp_372);
      real_t tmp_374 = tmp_1*(tmp_26*tmp_371 + tmp_33*tmp_373 + tmp_369*tmp_9 - 1.0/4.0) + tmp_2*(tmp_369*tmp_40 + tmp_371*tmp_41 + tmp_373*tmp_42 - 1.0/4.0) + tmp_6*(tmp_369*tmp_43 + tmp_371*tmp_44 + tmp_373*tmp_45 - 1.0/4.0);
      real_t tmp_375 = tmp_368 + tmp_76;
      real_t tmp_376 = tmp_375*tmp_72;
      real_t tmp_377 = tmp_375*tmp_73;
      real_t tmp_378 = tmp_370 + tmp_80;
      real_t tmp_379 = tmp_378*tmp_69;
      real_t tmp_380 = tmp_378*tmp_70;
      real_t tmp_381 = tmp_375*tmp_74;
      real_t tmp_382 = tmp_378*tmp_71;
      real_t tmp_383 = tmp_372 + tmp_86;
      real_t tmp_384 = tmp_383*tmp_66;
      real_t tmp_385 = tmp_383*tmp_67;
      real_t tmp_386 = tmp_383*tmp_68;
      real_t tmp_387 = -tmp_376 - tmp_377 - tmp_379 - tmp_380 - tmp_381 - tmp_382 - tmp_384 - tmp_385 - tmp_386 + 1;
      real_t tmp_388 = tmp_100*tmp_374;
      real_t tmp_389 = 0.020848748529055869*tmp_102;
      real_t tmp_390 = 0.93718850182767688*tmp_22 + 0.031405749086161582*tmp_23;
      real_t tmp_391 = tmp_19*(tmp_20 + tmp_390);
      real_t tmp_392 = 0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_30;
      real_t tmp_393 = tmp_19*(tmp_27 + tmp_392);
      real_t tmp_394 = 0.93718850182767688*tmp_36 + 0.031405749086161582*tmp_37;
      real_t tmp_395 = tmp_19*(tmp_34 + tmp_394);
      real_t tmp_396 = tmp_1*(tmp_26*tmp_393 + tmp_33*tmp_395 + tmp_391*tmp_9 - 1.0/4.0) + tmp_2*(tmp_391*tmp_40 + tmp_393*tmp_41 + tmp_395*tmp_42 - 1.0/4.0) + tmp_6*(tmp_391*tmp_43 + tmp_393*tmp_44 + tmp_395*tmp_45 - 1.0/4.0);
      real_t tmp_397 = tmp_390 + tmp_76;
      real_t tmp_398 = tmp_397*tmp_72;
      real_t tmp_399 = tmp_397*tmp_73;
      real_t tmp_400 = tmp_392 + tmp_80;
      real_t tmp_401 = tmp_400*tmp_69;
      real_t tmp_402 = tmp_400*tmp_70;
      real_t tmp_403 = tmp_397*tmp_74;
      real_t tmp_404 = tmp_400*tmp_71;
      real_t tmp_405 = tmp_394 + tmp_86;
      real_t tmp_406 = tmp_405*tmp_66;
      real_t tmp_407 = tmp_405*tmp_67;
      real_t tmp_408 = tmp_405*tmp_68;
      real_t tmp_409 = -tmp_398 - tmp_399 - tmp_401 - tmp_402 - tmp_403 - tmp_404 - tmp_406 - tmp_407 - tmp_408 + 1;
      real_t tmp_410 = tmp_100*tmp_396;
      real_t tmp_411 = 0.0068572537431980923*tmp_102;
      real_t tmp_412 = 0.60796128279561268*tmp_22 + 0.19601935860219369*tmp_23;
      real_t tmp_413 = tmp_19*(tmp_20 + tmp_412);
      real_t tmp_414 = 0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_30;
      real_t tmp_415 = tmp_19*(tmp_27 + tmp_414);
      real_t tmp_416 = 0.60796128279561268*tmp_36 + 0.19601935860219369*tmp_37;
      real_t tmp_417 = tmp_19*(tmp_34 + tmp_416);
      real_t tmp_418 = tmp_1*(tmp_26*tmp_415 + tmp_33*tmp_417 + tmp_413*tmp_9 - 1.0/4.0) + tmp_2*(tmp_40*tmp_413 + tmp_41*tmp_415 + tmp_417*tmp_42 - 1.0/4.0) + tmp_6*(tmp_413*tmp_43 + tmp_415*tmp_44 + tmp_417*tmp_45 - 1.0/4.0);
      real_t tmp_419 = tmp_412 + tmp_76;
      real_t tmp_420 = tmp_419*tmp_72;
      real_t tmp_421 = tmp_419*tmp_73;
      real_t tmp_422 = tmp_414 + tmp_80;
      real_t tmp_423 = tmp_422*tmp_69;
      real_t tmp_424 = tmp_422*tmp_70;
      real_t tmp_425 = tmp_419*tmp_74;
      real_t tmp_426 = tmp_422*tmp_71;
      real_t tmp_427 = tmp_416 + tmp_86;
      real_t tmp_428 = tmp_427*tmp_66;
      real_t tmp_429 = tmp_427*tmp_67;
      real_t tmp_430 = tmp_427*tmp_68;
      real_t tmp_431 = -tmp_420 - tmp_421 - tmp_423 - tmp_424 - tmp_425 - tmp_426 - tmp_428 - tmp_429 - tmp_430 + 1;
      real_t tmp_432 = tmp_100*tmp_418;
      real_t tmp_433 = 0.037198804536718075*tmp_102;
      real_t tmp_434 = 0.19107600050469298*tmp_22 + 0.40446199974765351*tmp_23;
      real_t tmp_435 = tmp_19*(tmp_20 + tmp_434);
      real_t tmp_436 = 0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_30;
      real_t tmp_437 = tmp_19*(tmp_27 + tmp_436);
      real_t tmp_438 = 0.19107600050469298*tmp_36 + 0.40446199974765351*tmp_37;
      real_t tmp_439 = tmp_19*(tmp_34 + tmp_438);
      real_t tmp_440 = tmp_1*(tmp_26*tmp_437 + tmp_33*tmp_439 + tmp_435*tmp_9 - 1.0/4.0) + tmp_2*(tmp_40*tmp_435 + tmp_41*tmp_437 + tmp_42*tmp_439 - 1.0/4.0) + tmp_6*(tmp_43*tmp_435 + tmp_437*tmp_44 + tmp_439*tmp_45 - 1.0/4.0);
      real_t tmp_441 = tmp_434 + tmp_76;
      real_t tmp_442 = tmp_441*tmp_72;
      real_t tmp_443 = tmp_441*tmp_73;
      real_t tmp_444 = tmp_436 + tmp_80;
      real_t tmp_445 = tmp_444*tmp_69;
      real_t tmp_446 = tmp_444*tmp_70;
      real_t tmp_447 = tmp_441*tmp_74;
      real_t tmp_448 = tmp_444*tmp_71;
      real_t tmp_449 = tmp_438 + tmp_86;
      real_t tmp_450 = tmp_449*tmp_66;
      real_t tmp_451 = tmp_449*tmp_67;
      real_t tmp_452 = tmp_449*tmp_68;
      real_t tmp_453 = -tmp_442 - tmp_443 - tmp_445 - tmp_446 - tmp_447 - tmp_448 - tmp_450 - tmp_451 - tmp_452 + 1;
      real_t tmp_454 = tmp_100*tmp_440;
      real_t tmp_455 = 0.042507265838595799*tmp_102;
      real_t tmp_456 = 0.031405749086161582*tmp_22 + 0.031405749086161582*tmp_23;
      real_t tmp_457 = tmp_19*(tmp_20 + tmp_456);
      real_t tmp_458 = 0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_30;
      real_t tmp_459 = tmp_19*(tmp_27 + tmp_458);
      real_t tmp_460 = 0.031405749086161582*tmp_36 + 0.031405749086161582*tmp_37;
      real_t tmp_461 = tmp_19*(tmp_34 + tmp_460);
      real_t tmp_462 = tmp_1*(tmp_26*tmp_459 + tmp_33*tmp_461 + tmp_457*tmp_9 - 1.0/4.0) + tmp_2*(tmp_40*tmp_457 + tmp_41*tmp_459 + tmp_42*tmp_461 - 1.0/4.0) + tmp_6*(tmp_43*tmp_457 + tmp_44*tmp_459 + tmp_45*tmp_461 - 1.0/4.0);
      real_t tmp_463 = tmp_456 + tmp_76;
      real_t tmp_464 = tmp_463*tmp_72;
      real_t tmp_465 = tmp_463*tmp_73;
      real_t tmp_466 = tmp_458 + tmp_80;
      real_t tmp_467 = tmp_466*tmp_69;
      real_t tmp_468 = tmp_466*tmp_70;
      real_t tmp_469 = tmp_463*tmp_74;
      real_t tmp_470 = tmp_466*tmp_71;
      real_t tmp_471 = tmp_460 + tmp_86;
      real_t tmp_472 = tmp_471*tmp_66;
      real_t tmp_473 = tmp_471*tmp_67;
      real_t tmp_474 = tmp_471*tmp_68;
      real_t tmp_475 = -tmp_464 - tmp_465 - tmp_467 - tmp_468 - tmp_469 - tmp_470 - tmp_472 - tmp_473 - tmp_474 + 1;
      real_t tmp_476 = tmp_100*tmp_462;
      real_t tmp_477 = 0.0068572537431980923*tmp_102;
      real_t tmp_478 = 0.19601935860219369*tmp_22 + 0.19601935860219369*tmp_23;
      real_t tmp_479 = tmp_19*(tmp_20 + tmp_478);
      real_t tmp_480 = 0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_30;
      real_t tmp_481 = tmp_19*(tmp_27 + tmp_480);
      real_t tmp_482 = 0.19601935860219369*tmp_36 + 0.19601935860219369*tmp_37;
      real_t tmp_483 = tmp_19*(tmp_34 + tmp_482);
      real_t tmp_484 = tmp_1*(tmp_26*tmp_481 + tmp_33*tmp_483 + tmp_479*tmp_9 - 1.0/4.0) + tmp_2*(tmp_40*tmp_479 + tmp_41*tmp_481 + tmp_42*tmp_483 - 1.0/4.0) + tmp_6*(tmp_43*tmp_479 + tmp_44*tmp_481 + tmp_45*tmp_483 - 1.0/4.0);
      real_t tmp_485 = tmp_478 + tmp_76;
      real_t tmp_486 = tmp_485*tmp_72;
      real_t tmp_487 = tmp_485*tmp_73;
      real_t tmp_488 = tmp_480 + tmp_80;
      real_t tmp_489 = tmp_488*tmp_69;
      real_t tmp_490 = tmp_488*tmp_70;
      real_t tmp_491 = tmp_485*tmp_74;
      real_t tmp_492 = tmp_488*tmp_71;
      real_t tmp_493 = tmp_482 + tmp_86;
      real_t tmp_494 = tmp_493*tmp_66;
      real_t tmp_495 = tmp_493*tmp_67;
      real_t tmp_496 = tmp_493*tmp_68;
      real_t tmp_497 = -tmp_486 - tmp_487 - tmp_489 - tmp_490 - tmp_491 - tmp_492 - tmp_494 - tmp_495 - tmp_496 + 1;
      real_t tmp_498 = tmp_100*tmp_484;
      real_t tmp_499 = 0.037198804536718075*tmp_102;
      real_t tmp_500 = 0.40446199974765351*tmp_22 + 0.40446199974765351*tmp_23;
      real_t tmp_501 = tmp_19*(tmp_20 + tmp_500);
      real_t tmp_502 = 0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_30;
      real_t tmp_503 = tmp_19*(tmp_27 + tmp_502);
      real_t tmp_504 = 0.40446199974765351*tmp_36 + 0.40446199974765351*tmp_37;
      real_t tmp_505 = tmp_19*(tmp_34 + tmp_504);
      real_t tmp_506 = tmp_1*(tmp_26*tmp_503 + tmp_33*tmp_505 + tmp_501*tmp_9 - 1.0/4.0) + tmp_2*(tmp_40*tmp_501 + tmp_41*tmp_503 + tmp_42*tmp_505 - 1.0/4.0) + tmp_6*(tmp_43*tmp_501 + tmp_44*tmp_503 + tmp_45*tmp_505 - 1.0/4.0);
      real_t tmp_507 = tmp_500 + tmp_76;
      real_t tmp_508 = tmp_507*tmp_72;
      real_t tmp_509 = tmp_507*tmp_73;
      real_t tmp_510 = tmp_502 + tmp_80;
      real_t tmp_511 = tmp_510*tmp_69;
      real_t tmp_512 = tmp_510*tmp_70;
      real_t tmp_513 = tmp_507*tmp_74;
      real_t tmp_514 = tmp_510*tmp_71;
      real_t tmp_515 = tmp_504 + tmp_86;
      real_t tmp_516 = tmp_515*tmp_66;
      real_t tmp_517 = tmp_515*tmp_67;
      real_t tmp_518 = tmp_515*tmp_68;
      real_t tmp_519 = -tmp_508 - tmp_509 - tmp_511 - tmp_512 - tmp_513 - tmp_514 - tmp_516 - tmp_517 - tmp_518 + 1;
      real_t tmp_520 = tmp_100*tmp_506;
      real_t tmp_521 = 0.042507265838595799*tmp_102;
      real_t tmp_522 = 0.1711304259088916*tmp_22 + 0.041227165399737475*tmp_23;
      real_t tmp_523 = tmp_19*(tmp_20 + tmp_522);
      real_t tmp_524 = 0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_30;
      real_t tmp_525 = tmp_19*(tmp_27 + tmp_524);
      real_t tmp_526 = 0.1711304259088916*tmp_36 + 0.041227165399737475*tmp_37;
      real_t tmp_527 = tmp_19*(tmp_34 + tmp_526);
      real_t tmp_528 = tmp_1*(tmp_26*tmp_525 + tmp_33*tmp_527 + tmp_523*tmp_9 - 1.0/4.0) + tmp_2*(tmp_40*tmp_523 + tmp_41*tmp_525 + tmp_42*tmp_527 - 1.0/4.0) + tmp_6*(tmp_43*tmp_523 + tmp_44*tmp_525 + tmp_45*tmp_527 - 1.0/4.0);
      real_t tmp_529 = tmp_522 + tmp_76;
      real_t tmp_530 = tmp_529*tmp_72;
      real_t tmp_531 = tmp_529*tmp_73;
      real_t tmp_532 = tmp_524 + tmp_80;
      real_t tmp_533 = tmp_532*tmp_69;
      real_t tmp_534 = tmp_532*tmp_70;
      real_t tmp_535 = tmp_529*tmp_74;
      real_t tmp_536 = tmp_532*tmp_71;
      real_t tmp_537 = tmp_526 + tmp_86;
      real_t tmp_538 = tmp_537*tmp_66;
      real_t tmp_539 = tmp_537*tmp_67;
      real_t tmp_540 = tmp_537*tmp_68;
      real_t tmp_541 = -tmp_530 - tmp_531 - tmp_533 - tmp_534 - tmp_535 - tmp_536 - tmp_538 - tmp_539 - tmp_540 + 1;
      real_t tmp_542 = tmp_100*tmp_528;
      real_t tmp_543 = 0.019202922745021479*tmp_102;
      real_t tmp_544 = tmp_84 + tmp_85 + tmp_90;
      real_t tmp_545 = 0.5*p_affine_13_0*tmp_68 + 0.5*p_affine_13_1*tmp_71 + 0.5*p_affine_13_2*tmp_74;
      real_t tmp_546 = tmp_117 + tmp_118 + tmp_122;
      real_t tmp_547 = tmp_139 + tmp_140 + tmp_144;
      real_t tmp_548 = tmp_161 + tmp_162 + tmp_166;
      real_t tmp_549 = tmp_183 + tmp_184 + tmp_188;
      real_t tmp_550 = tmp_205 + tmp_206 + tmp_210;
      real_t tmp_551 = tmp_227 + tmp_228 + tmp_232;
      real_t tmp_552 = tmp_249 + tmp_250 + tmp_254;
      real_t tmp_553 = tmp_271 + tmp_272 + tmp_276;
      real_t tmp_554 = tmp_293 + tmp_294 + tmp_298;
      real_t tmp_555 = tmp_315 + tmp_316 + tmp_320;
      real_t tmp_556 = tmp_337 + tmp_338 + tmp_342;
      real_t tmp_557 = tmp_359 + tmp_360 + tmp_364;
      real_t tmp_558 = tmp_381 + tmp_382 + tmp_386;
      real_t tmp_559 = tmp_403 + tmp_404 + tmp_408;
      real_t tmp_560 = tmp_425 + tmp_426 + tmp_430;
      real_t tmp_561 = tmp_447 + tmp_448 + tmp_452;
      real_t tmp_562 = tmp_469 + tmp_470 + tmp_474;
      real_t tmp_563 = tmp_491 + tmp_492 + tmp_496;
      real_t tmp_564 = tmp_513 + tmp_514 + tmp_518;
      real_t tmp_565 = tmp_535 + tmp_536 + tmp_540;
      real_t tmp_566 = tmp_79 + tmp_83 + tmp_89;
      real_t tmp_567 = 0.5*p_affine_13_0*tmp_67 + 0.5*p_affine_13_1*tmp_70 + 0.5*p_affine_13_2*tmp_73;
      real_t tmp_568 = tmp_113 + tmp_116 + tmp_121;
      real_t tmp_569 = tmp_135 + tmp_138 + tmp_143;
      real_t tmp_570 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_571 = tmp_179 + tmp_182 + tmp_187;
      real_t tmp_572 = tmp_201 + tmp_204 + tmp_209;
      real_t tmp_573 = tmp_223 + tmp_226 + tmp_231;
      real_t tmp_574 = tmp_245 + tmp_248 + tmp_253;
      real_t tmp_575 = tmp_267 + tmp_270 + tmp_275;
      real_t tmp_576 = tmp_289 + tmp_292 + tmp_297;
      real_t tmp_577 = tmp_311 + tmp_314 + tmp_319;
      real_t tmp_578 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_579 = tmp_355 + tmp_358 + tmp_363;
      real_t tmp_580 = tmp_377 + tmp_380 + tmp_385;
      real_t tmp_581 = tmp_399 + tmp_402 + tmp_407;
      real_t tmp_582 = tmp_421 + tmp_424 + tmp_429;
      real_t tmp_583 = tmp_443 + tmp_446 + tmp_451;
      real_t tmp_584 = tmp_465 + tmp_468 + tmp_473;
      real_t tmp_585 = tmp_487 + tmp_490 + tmp_495;
      real_t tmp_586 = tmp_509 + tmp_512 + tmp_517;
      real_t tmp_587 = tmp_531 + tmp_534 + tmp_539;
      real_t tmp_588 = tmp_78 + tmp_82 + tmp_88;
      real_t tmp_589 = 0.5*p_affine_13_0*tmp_66 + 0.5*p_affine_13_1*tmp_69 + 0.5*p_affine_13_2*tmp_72;
      real_t tmp_590 = tmp_112 + tmp_115 + tmp_120;
      real_t tmp_591 = tmp_134 + tmp_137 + tmp_142;
      real_t tmp_592 = tmp_156 + tmp_159 + tmp_164;
      real_t tmp_593 = tmp_178 + tmp_181 + tmp_186;
      real_t tmp_594 = tmp_200 + tmp_203 + tmp_208;
      real_t tmp_595 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_596 = tmp_244 + tmp_247 + tmp_252;
      real_t tmp_597 = tmp_266 + tmp_269 + tmp_274;
      real_t tmp_598 = tmp_288 + tmp_291 + tmp_296;
      real_t tmp_599 = tmp_310 + tmp_313 + tmp_318;
      real_t tmp_600 = tmp_332 + tmp_335 + tmp_340;
      real_t tmp_601 = tmp_354 + tmp_357 + tmp_362;
      real_t tmp_602 = tmp_376 + tmp_379 + tmp_384;
      real_t tmp_603 = tmp_398 + tmp_401 + tmp_406;
      real_t tmp_604 = tmp_420 + tmp_423 + tmp_428;
      real_t tmp_605 = tmp_442 + tmp_445 + tmp_450;
      real_t tmp_606 = tmp_464 + tmp_467 + tmp_472;
      real_t tmp_607 = tmp_486 + tmp_489 + tmp_494;
      real_t tmp_608 = tmp_508 + tmp_511 + tmp_516;
      real_t tmp_609 = tmp_530 + tmp_533 + tmp_538;
      real_t a_0_0 = tmp_103*(-tmp_101*tmp_91 + tmp_46*tmp_75 - tmp_91*tmp_95) + tmp_125*(tmp_110*tmp_75 - tmp_123*tmp_124 - tmp_123*tmp_95) + tmp_147*(tmp_132*tmp_75 - tmp_145*tmp_146 - tmp_145*tmp_95) + tmp_169*(tmp_154*tmp_75 - tmp_167*tmp_168 - tmp_167*tmp_95) + tmp_191*(tmp_176*tmp_75 - tmp_189*tmp_190 - tmp_189*tmp_95) + tmp_213*(tmp_198*tmp_75 - tmp_211*tmp_212 - tmp_211*tmp_95) + tmp_235*(tmp_220*tmp_75 - tmp_233*tmp_234 - tmp_233*tmp_95) + tmp_257*(tmp_242*tmp_75 - tmp_255*tmp_256 - tmp_255*tmp_95) + tmp_279*(tmp_264*tmp_75 - tmp_277*tmp_278 - tmp_277*tmp_95) + tmp_301*(tmp_286*tmp_75 - tmp_299*tmp_300 - tmp_299*tmp_95) + tmp_323*(tmp_308*tmp_75 - tmp_321*tmp_322 - tmp_321*tmp_95) + tmp_345*(tmp_330*tmp_75 - tmp_343*tmp_344 - tmp_343*tmp_95) + tmp_367*(tmp_352*tmp_75 - tmp_365*tmp_366 - tmp_365*tmp_95) + tmp_389*(tmp_374*tmp_75 - tmp_387*tmp_388 - tmp_387*tmp_95) + tmp_411*(tmp_396*tmp_75 - tmp_409*tmp_410 - tmp_409*tmp_95) + tmp_433*(tmp_418*tmp_75 - tmp_431*tmp_432 - tmp_431*tmp_95) + tmp_455*(tmp_440*tmp_75 - tmp_453*tmp_454 - tmp_453*tmp_95) + tmp_477*(tmp_462*tmp_75 - tmp_475*tmp_476 - tmp_475*tmp_95) + tmp_499*(tmp_484*tmp_75 - tmp_497*tmp_498 - tmp_497*tmp_95) + tmp_521*(tmp_506*tmp_75 - tmp_519*tmp_520 - tmp_519*tmp_95) + tmp_543*(tmp_528*tmp_75 - tmp_541*tmp_542 - tmp_541*tmp_95);
      real_t a_1_0 = tmp_103*(-tmp_101*tmp_544 + tmp_46*tmp_545 - tmp_544*tmp_95) + tmp_125*(tmp_110*tmp_545 - tmp_124*tmp_546 - tmp_546*tmp_95) + tmp_147*(tmp_132*tmp_545 - tmp_146*tmp_547 - tmp_547*tmp_95) + tmp_169*(tmp_154*tmp_545 - tmp_168*tmp_548 - tmp_548*tmp_95) + tmp_191*(tmp_176*tmp_545 - tmp_190*tmp_549 - tmp_549*tmp_95) + tmp_213*(tmp_198*tmp_545 - tmp_212*tmp_550 - tmp_550*tmp_95) + tmp_235*(tmp_220*tmp_545 - tmp_234*tmp_551 - tmp_551*tmp_95) + tmp_257*(tmp_242*tmp_545 - tmp_256*tmp_552 - tmp_552*tmp_95) + tmp_279*(tmp_264*tmp_545 - tmp_278*tmp_553 - tmp_553*tmp_95) + tmp_301*(tmp_286*tmp_545 - tmp_300*tmp_554 - tmp_554*tmp_95) + tmp_323*(tmp_308*tmp_545 - tmp_322*tmp_555 - tmp_555*tmp_95) + tmp_345*(tmp_330*tmp_545 - tmp_344*tmp_556 - tmp_556*tmp_95) + tmp_367*(tmp_352*tmp_545 - tmp_366*tmp_557 - tmp_557*tmp_95) + tmp_389*(tmp_374*tmp_545 - tmp_388*tmp_558 - tmp_558*tmp_95) + tmp_411*(tmp_396*tmp_545 - tmp_410*tmp_559 - tmp_559*tmp_95) + tmp_433*(tmp_418*tmp_545 - tmp_432*tmp_560 - tmp_560*tmp_95) + tmp_455*(tmp_440*tmp_545 - tmp_454*tmp_561 - tmp_561*tmp_95) + tmp_477*(tmp_462*tmp_545 - tmp_476*tmp_562 - tmp_562*tmp_95) + tmp_499*(tmp_484*tmp_545 - tmp_498*tmp_563 - tmp_563*tmp_95) + tmp_521*(tmp_506*tmp_545 - tmp_520*tmp_564 - tmp_564*tmp_95) + tmp_543*(tmp_528*tmp_545 - tmp_542*tmp_565 - tmp_565*tmp_95);
      real_t a_2_0 = tmp_103*(-tmp_101*tmp_566 + tmp_46*tmp_567 - tmp_566*tmp_95) + tmp_125*(tmp_110*tmp_567 - tmp_124*tmp_568 - tmp_568*tmp_95) + tmp_147*(tmp_132*tmp_567 - tmp_146*tmp_569 - tmp_569*tmp_95) + tmp_169*(tmp_154*tmp_567 - tmp_168*tmp_570 - tmp_570*tmp_95) + tmp_191*(tmp_176*tmp_567 - tmp_190*tmp_571 - tmp_571*tmp_95) + tmp_213*(tmp_198*tmp_567 - tmp_212*tmp_572 - tmp_572*tmp_95) + tmp_235*(tmp_220*tmp_567 - tmp_234*tmp_573 - tmp_573*tmp_95) + tmp_257*(tmp_242*tmp_567 - tmp_256*tmp_574 - tmp_574*tmp_95) + tmp_279*(tmp_264*tmp_567 - tmp_278*tmp_575 - tmp_575*tmp_95) + tmp_301*(tmp_286*tmp_567 - tmp_300*tmp_576 - tmp_576*tmp_95) + tmp_323*(tmp_308*tmp_567 - tmp_322*tmp_577 - tmp_577*tmp_95) + tmp_345*(tmp_330*tmp_567 - tmp_344*tmp_578 - tmp_578*tmp_95) + tmp_367*(tmp_352*tmp_567 - tmp_366*tmp_579 - tmp_579*tmp_95) + tmp_389*(tmp_374*tmp_567 - tmp_388*tmp_580 - tmp_580*tmp_95) + tmp_411*(tmp_396*tmp_567 - tmp_410*tmp_581 - tmp_581*tmp_95) + tmp_433*(tmp_418*tmp_567 - tmp_432*tmp_582 - tmp_582*tmp_95) + tmp_455*(tmp_440*tmp_567 - tmp_454*tmp_583 - tmp_583*tmp_95) + tmp_477*(tmp_462*tmp_567 - tmp_476*tmp_584 - tmp_584*tmp_95) + tmp_499*(tmp_484*tmp_567 - tmp_498*tmp_585 - tmp_585*tmp_95) + tmp_521*(tmp_506*tmp_567 - tmp_520*tmp_586 - tmp_586*tmp_95) + tmp_543*(tmp_528*tmp_567 - tmp_542*tmp_587 - tmp_587*tmp_95);
      real_t a_3_0 = tmp_103*(-tmp_101*tmp_588 + tmp_46*tmp_589 - tmp_588*tmp_95) + tmp_125*(tmp_110*tmp_589 - tmp_124*tmp_590 - tmp_590*tmp_95) + tmp_147*(tmp_132*tmp_589 - tmp_146*tmp_591 - tmp_591*tmp_95) + tmp_169*(tmp_154*tmp_589 - tmp_168*tmp_592 - tmp_592*tmp_95) + tmp_191*(tmp_176*tmp_589 - tmp_190*tmp_593 - tmp_593*tmp_95) + tmp_213*(tmp_198*tmp_589 - tmp_212*tmp_594 - tmp_594*tmp_95) + tmp_235*(tmp_220*tmp_589 - tmp_234*tmp_595 - tmp_595*tmp_95) + tmp_257*(tmp_242*tmp_589 - tmp_256*tmp_596 - tmp_596*tmp_95) + tmp_279*(tmp_264*tmp_589 - tmp_278*tmp_597 - tmp_597*tmp_95) + tmp_301*(tmp_286*tmp_589 - tmp_300*tmp_598 - tmp_598*tmp_95) + tmp_323*(tmp_308*tmp_589 - tmp_322*tmp_599 - tmp_599*tmp_95) + tmp_345*(tmp_330*tmp_589 - tmp_344*tmp_600 - tmp_600*tmp_95) + tmp_367*(tmp_352*tmp_589 - tmp_366*tmp_601 - tmp_601*tmp_95) + tmp_389*(tmp_374*tmp_589 - tmp_388*tmp_602 - tmp_602*tmp_95) + tmp_411*(tmp_396*tmp_589 - tmp_410*tmp_603 - tmp_603*tmp_95) + tmp_433*(tmp_418*tmp_589 - tmp_432*tmp_604 - tmp_604*tmp_95) + tmp_455*(tmp_440*tmp_589 - tmp_454*tmp_605 - tmp_605*tmp_95) + tmp_477*(tmp_462*tmp_589 - tmp_476*tmp_606 - tmp_606*tmp_95) + tmp_499*(tmp_484*tmp_589 - tmp_498*tmp_607 - tmp_607*tmp_95) + tmp_521*(tmp_506*tmp_589 - tmp_520*tmp_608 - tmp_608*tmp_95) + tmp_543*(tmp_528*tmp_589 - tmp_542*tmp_609 - tmp_609*tmp_95);
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
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_2;
      real_t tmp_8 = p_affine_3_2 + tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_7;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_4;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_11*tmp_14;
      real_t tmp_17 = tmp_4*tmp_8;
      real_t tmp_18 = tmp_13*tmp_3;
      real_t tmp_19 = 1.0 / (-tmp_1*tmp_16 + tmp_1*tmp_9 + tmp_10*tmp_12 - tmp_10*tmp_18 + tmp_15*tmp_5 - tmp_17*tmp_5);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_6;
      real_t tmp_26 = -tmp_1*tmp_11 + tmp_13*tmp_5;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_14 + tmp_10*tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_8 - tmp_10*tmp_13;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_12 - tmp_18;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = tmp_15 - tmp_17;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_3 + tmp_14*tmp_5;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = tmp_10*tmp_11 - tmp_5*tmp_8;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_53 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_54 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_55 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_58 = 5.0*std::pow((std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)*std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)) + (std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)*std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)) + (std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)*std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)), 0.25);
      real_t tmp_59 = 0.0068572537431980923*tmp_58*(tmp_1*(tmp_52 - 1.0/4.0) + tmp_13*(tmp_54 - 1.0/4.0) + tmp_4*(tmp_53 - 1.0/4.0));
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_61 = tmp_6*tmp_60;
      real_t tmp_62 = tmp_26*tmp_60;
      real_t tmp_63 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_64 = tmp_28*tmp_63;
      real_t tmp_65 = tmp_35*tmp_63;
      real_t tmp_66 = tmp_37*tmp_60;
      real_t tmp_67 = tmp_39*tmp_63;
      real_t tmp_68 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_69 = tmp_41*tmp_68;
      real_t tmp_70 = tmp_48*tmp_68;
      real_t tmp_71 = tmp_50*tmp_68;
      real_t tmp_72 = tmp_66 + tmp_67 + tmp_71;
      real_t tmp_73 = tmp_62 + tmp_65 + tmp_70;
      real_t tmp_74 = tmp_61 + tmp_64 + tmp_69;
      real_t tmp_75 = 0.037198804536718075*tmp_58*(tmp_1*(tmp_72 - 1.0/4.0) + tmp_13*(tmp_74 - 1.0/4.0) + tmp_4*(tmp_73 - 1.0/4.0));
      real_t tmp_76 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_77 = tmp_6*tmp_76;
      real_t tmp_78 = tmp_26*tmp_76;
      real_t tmp_79 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_80 = tmp_28*tmp_79;
      real_t tmp_81 = tmp_35*tmp_79;
      real_t tmp_82 = tmp_37*tmp_76;
      real_t tmp_83 = tmp_39*tmp_79;
      real_t tmp_84 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_85 = tmp_41*tmp_84;
      real_t tmp_86 = tmp_48*tmp_84;
      real_t tmp_87 = tmp_50*tmp_84;
      real_t tmp_88 = tmp_82 + tmp_83 + tmp_87;
      real_t tmp_89 = tmp_78 + tmp_81 + tmp_86;
      real_t tmp_90 = tmp_77 + tmp_80 + tmp_85;
      real_t tmp_91 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_88 - 1.0/4.0) + tmp_13*(tmp_90 - 1.0/4.0) + tmp_4*(tmp_89 - 1.0/4.0));
      real_t tmp_92 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_93 = tmp_6*tmp_92;
      real_t tmp_94 = tmp_26*tmp_92;
      real_t tmp_95 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_96 = tmp_28*tmp_95;
      real_t tmp_97 = tmp_35*tmp_95;
      real_t tmp_98 = tmp_37*tmp_92;
      real_t tmp_99 = tmp_39*tmp_95;
      real_t tmp_100 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_101 = tmp_100*tmp_41;
      real_t tmp_102 = tmp_100*tmp_48;
      real_t tmp_103 = tmp_100*tmp_50;
      real_t tmp_104 = tmp_103 + tmp_98 + tmp_99;
      real_t tmp_105 = tmp_102 + tmp_94 + tmp_97;
      real_t tmp_106 = tmp_101 + tmp_93 + tmp_96;
      real_t tmp_107 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_104 - 1.0/4.0) + tmp_13*(tmp_106 - 1.0/4.0) + tmp_4*(tmp_105 - 1.0/4.0));
      real_t tmp_108 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_109 = tmp_108*tmp_6;
      real_t tmp_110 = tmp_108*tmp_26;
      real_t tmp_111 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_112 = tmp_111*tmp_28;
      real_t tmp_113 = tmp_111*tmp_35;
      real_t tmp_114 = tmp_108*tmp_37;
      real_t tmp_115 = tmp_111*tmp_39;
      real_t tmp_116 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_117 = tmp_116*tmp_41;
      real_t tmp_118 = tmp_116*tmp_48;
      real_t tmp_119 = tmp_116*tmp_50;
      real_t tmp_120 = tmp_114 + tmp_115 + tmp_119;
      real_t tmp_121 = tmp_110 + tmp_113 + tmp_118;
      real_t tmp_122 = tmp_109 + tmp_112 + tmp_117;
      real_t tmp_123 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_120 - 1.0/4.0) + tmp_13*(tmp_122 - 1.0/4.0) + tmp_4*(tmp_121 - 1.0/4.0));
      real_t tmp_124 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_125 = tmp_124*tmp_6;
      real_t tmp_126 = tmp_124*tmp_26;
      real_t tmp_127 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_128 = tmp_127*tmp_28;
      real_t tmp_129 = tmp_127*tmp_35;
      real_t tmp_130 = tmp_124*tmp_37;
      real_t tmp_131 = tmp_127*tmp_39;
      real_t tmp_132 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_133 = tmp_132*tmp_41;
      real_t tmp_134 = tmp_132*tmp_48;
      real_t tmp_135 = tmp_132*tmp_50;
      real_t tmp_136 = tmp_130 + tmp_131 + tmp_135;
      real_t tmp_137 = tmp_126 + tmp_129 + tmp_134;
      real_t tmp_138 = tmp_125 + tmp_128 + tmp_133;
      real_t tmp_139 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_136 - 1.0/4.0) + tmp_13*(tmp_138 - 1.0/4.0) + tmp_4*(tmp_137 - 1.0/4.0));
      real_t tmp_140 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_141 = tmp_140*tmp_6;
      real_t tmp_142 = tmp_140*tmp_26;
      real_t tmp_143 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_144 = tmp_143*tmp_28;
      real_t tmp_145 = tmp_143*tmp_35;
      real_t tmp_146 = tmp_140*tmp_37;
      real_t tmp_147 = tmp_143*tmp_39;
      real_t tmp_148 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_149 = tmp_148*tmp_41;
      real_t tmp_150 = tmp_148*tmp_48;
      real_t tmp_151 = tmp_148*tmp_50;
      real_t tmp_152 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_153 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_154 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_155 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_152 - 1.0/4.0) + tmp_13*(tmp_154 - 1.0/4.0) + tmp_4*(tmp_153 - 1.0/4.0));
      real_t tmp_156 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_157 = tmp_156*tmp_6;
      real_t tmp_158 = tmp_156*tmp_26;
      real_t tmp_159 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_160 = tmp_159*tmp_28;
      real_t tmp_161 = tmp_159*tmp_35;
      real_t tmp_162 = tmp_156*tmp_37;
      real_t tmp_163 = tmp_159*tmp_39;
      real_t tmp_164 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_165 = tmp_164*tmp_41;
      real_t tmp_166 = tmp_164*tmp_48;
      real_t tmp_167 = tmp_164*tmp_50;
      real_t tmp_168 = tmp_162 + tmp_163 + tmp_167;
      real_t tmp_169 = tmp_158 + tmp_161 + tmp_166;
      real_t tmp_170 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_171 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_168 - 1.0/4.0) + tmp_13*(tmp_170 - 1.0/4.0) + tmp_4*(tmp_169 - 1.0/4.0));
      real_t tmp_172 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_173 = tmp_172*tmp_6;
      real_t tmp_174 = tmp_172*tmp_26;
      real_t tmp_175 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_176 = tmp_175*tmp_28;
      real_t tmp_177 = tmp_175*tmp_35;
      real_t tmp_178 = tmp_172*tmp_37;
      real_t tmp_179 = tmp_175*tmp_39;
      real_t tmp_180 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_181 = tmp_180*tmp_41;
      real_t tmp_182 = tmp_180*tmp_48;
      real_t tmp_183 = tmp_180*tmp_50;
      real_t tmp_184 = tmp_178 + tmp_179 + tmp_183;
      real_t tmp_185 = tmp_174 + tmp_177 + tmp_182;
      real_t tmp_186 = tmp_173 + tmp_176 + tmp_181;
      real_t tmp_187 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_184 - 1.0/4.0) + tmp_13*(tmp_186 - 1.0/4.0) + tmp_4*(tmp_185 - 1.0/4.0));
      real_t tmp_188 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_189 = tmp_188*tmp_6;
      real_t tmp_190 = tmp_188*tmp_26;
      real_t tmp_191 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_192 = tmp_191*tmp_28;
      real_t tmp_193 = tmp_191*tmp_35;
      real_t tmp_194 = tmp_188*tmp_37;
      real_t tmp_195 = tmp_191*tmp_39;
      real_t tmp_196 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_197 = tmp_196*tmp_41;
      real_t tmp_198 = tmp_196*tmp_48;
      real_t tmp_199 = tmp_196*tmp_50;
      real_t tmp_200 = tmp_194 + tmp_195 + tmp_199;
      real_t tmp_201 = tmp_190 + tmp_193 + tmp_198;
      real_t tmp_202 = tmp_189 + tmp_192 + tmp_197;
      real_t tmp_203 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_200 - 1.0/4.0) + tmp_13*(tmp_202 - 1.0/4.0) + tmp_4*(tmp_201 - 1.0/4.0));
      real_t tmp_204 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_205 = tmp_204*tmp_6;
      real_t tmp_206 = tmp_204*tmp_26;
      real_t tmp_207 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_208 = tmp_207*tmp_28;
      real_t tmp_209 = tmp_207*tmp_35;
      real_t tmp_210 = tmp_204*tmp_37;
      real_t tmp_211 = tmp_207*tmp_39;
      real_t tmp_212 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_213 = tmp_212*tmp_41;
      real_t tmp_214 = tmp_212*tmp_48;
      real_t tmp_215 = tmp_212*tmp_50;
      real_t tmp_216 = tmp_210 + tmp_211 + tmp_215;
      real_t tmp_217 = tmp_206 + tmp_209 + tmp_214;
      real_t tmp_218 = tmp_205 + tmp_208 + tmp_213;
      real_t tmp_219 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_216 - 1.0/4.0) + tmp_13*(tmp_218 - 1.0/4.0) + tmp_4*(tmp_217 - 1.0/4.0));
      real_t tmp_220 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_221 = tmp_220*tmp_6;
      real_t tmp_222 = tmp_220*tmp_26;
      real_t tmp_223 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_224 = tmp_223*tmp_28;
      real_t tmp_225 = tmp_223*tmp_35;
      real_t tmp_226 = tmp_220*tmp_37;
      real_t tmp_227 = tmp_223*tmp_39;
      real_t tmp_228 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_229 = tmp_228*tmp_41;
      real_t tmp_230 = tmp_228*tmp_48;
      real_t tmp_231 = tmp_228*tmp_50;
      real_t tmp_232 = tmp_226 + tmp_227 + tmp_231;
      real_t tmp_233 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_234 = tmp_221 + tmp_224 + tmp_229;
      real_t tmp_235 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_232 - 1.0/4.0) + tmp_13*(tmp_234 - 1.0/4.0) + tmp_4*(tmp_233 - 1.0/4.0));
      real_t tmp_236 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_237 = tmp_236*tmp_6;
      real_t tmp_238 = tmp_236*tmp_26;
      real_t tmp_239 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_240 = tmp_239*tmp_28;
      real_t tmp_241 = tmp_239*tmp_35;
      real_t tmp_242 = tmp_236*tmp_37;
      real_t tmp_243 = tmp_239*tmp_39;
      real_t tmp_244 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_245 = tmp_244*tmp_41;
      real_t tmp_246 = tmp_244*tmp_48;
      real_t tmp_247 = tmp_244*tmp_50;
      real_t tmp_248 = tmp_242 + tmp_243 + tmp_247;
      real_t tmp_249 = tmp_238 + tmp_241 + tmp_246;
      real_t tmp_250 = tmp_237 + tmp_240 + tmp_245;
      real_t tmp_251 = 0.042507265838595799*tmp_58*(tmp_1*(tmp_248 - 1.0/4.0) + tmp_13*(tmp_250 - 1.0/4.0) + tmp_4*(tmp_249 - 1.0/4.0));
      real_t tmp_252 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_6;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = tmp_258 + tmp_259 + tmp_263;
      real_t tmp_265 = tmp_254 + tmp_257 + tmp_262;
      real_t tmp_266 = tmp_253 + tmp_256 + tmp_261;
      real_t tmp_267 = 0.020848748529055869*tmp_58*(tmp_1*(tmp_264 - 1.0/4.0) + tmp_13*(tmp_266 - 1.0/4.0) + tmp_4*(tmp_265 - 1.0/4.0));
      real_t tmp_268 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_269 = tmp_268*tmp_6;
      real_t tmp_270 = tmp_26*tmp_268;
      real_t tmp_271 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_272 = tmp_271*tmp_28;
      real_t tmp_273 = tmp_271*tmp_35;
      real_t tmp_274 = tmp_268*tmp_37;
      real_t tmp_275 = tmp_271*tmp_39;
      real_t tmp_276 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_277 = tmp_276*tmp_41;
      real_t tmp_278 = tmp_276*tmp_48;
      real_t tmp_279 = tmp_276*tmp_50;
      real_t tmp_280 = tmp_274 + tmp_275 + tmp_279;
      real_t tmp_281 = tmp_270 + tmp_273 + tmp_278;
      real_t tmp_282 = tmp_269 + tmp_272 + tmp_277;
      real_t tmp_283 = 0.0068572537431980923*tmp_58*(tmp_1*(tmp_280 - 1.0/4.0) + tmp_13*(tmp_282 - 1.0/4.0) + tmp_4*(tmp_281 - 1.0/4.0));
      real_t tmp_284 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_285 = tmp_284*tmp_6;
      real_t tmp_286 = tmp_26*tmp_284;
      real_t tmp_287 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_288 = tmp_28*tmp_287;
      real_t tmp_289 = tmp_287*tmp_35;
      real_t tmp_290 = tmp_284*tmp_37;
      real_t tmp_291 = tmp_287*tmp_39;
      real_t tmp_292 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_293 = tmp_292*tmp_41;
      real_t tmp_294 = tmp_292*tmp_48;
      real_t tmp_295 = tmp_292*tmp_50;
      real_t tmp_296 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_297 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_298 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_299 = 0.037198804536718075*tmp_58*(tmp_1*(tmp_296 - 1.0/4.0) + tmp_13*(tmp_298 - 1.0/4.0) + tmp_4*(tmp_297 - 1.0/4.0));
      real_t tmp_300 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_301 = tmp_300*tmp_6;
      real_t tmp_302 = tmp_26*tmp_300;
      real_t tmp_303 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_304 = tmp_28*tmp_303;
      real_t tmp_305 = tmp_303*tmp_35;
      real_t tmp_306 = tmp_300*tmp_37;
      real_t tmp_307 = tmp_303*tmp_39;
      real_t tmp_308 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_309 = tmp_308*tmp_41;
      real_t tmp_310 = tmp_308*tmp_48;
      real_t tmp_311 = tmp_308*tmp_50;
      real_t tmp_312 = tmp_306 + tmp_307 + tmp_311;
      real_t tmp_313 = tmp_302 + tmp_305 + tmp_310;
      real_t tmp_314 = tmp_301 + tmp_304 + tmp_309;
      real_t tmp_315 = 0.042507265838595799*tmp_58*(tmp_1*(tmp_312 - 1.0/4.0) + tmp_13*(tmp_314 - 1.0/4.0) + tmp_4*(tmp_313 - 1.0/4.0));
      real_t tmp_316 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_317 = tmp_316*tmp_6;
      real_t tmp_318 = tmp_26*tmp_316;
      real_t tmp_319 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_320 = tmp_28*tmp_319;
      real_t tmp_321 = tmp_319*tmp_35;
      real_t tmp_322 = tmp_316*tmp_37;
      real_t tmp_323 = tmp_319*tmp_39;
      real_t tmp_324 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_325 = tmp_324*tmp_41;
      real_t tmp_326 = tmp_324*tmp_48;
      real_t tmp_327 = tmp_324*tmp_50;
      real_t tmp_328 = tmp_322 + tmp_323 + tmp_327;
      real_t tmp_329 = tmp_318 + tmp_321 + tmp_326;
      real_t tmp_330 = tmp_317 + tmp_320 + tmp_325;
      real_t tmp_331 = 0.0068572537431980923*tmp_58*(tmp_1*(tmp_328 - 1.0/4.0) + tmp_13*(tmp_330 - 1.0/4.0) + tmp_4*(tmp_329 - 1.0/4.0));
      real_t tmp_332 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_333 = tmp_332*tmp_6;
      real_t tmp_334 = tmp_26*tmp_332;
      real_t tmp_335 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_336 = tmp_28*tmp_335;
      real_t tmp_337 = tmp_335*tmp_35;
      real_t tmp_338 = tmp_332*tmp_37;
      real_t tmp_339 = tmp_335*tmp_39;
      real_t tmp_340 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_341 = tmp_340*tmp_41;
      real_t tmp_342 = tmp_340*tmp_48;
      real_t tmp_343 = tmp_340*tmp_50;
      real_t tmp_344 = tmp_338 + tmp_339 + tmp_343;
      real_t tmp_345 = tmp_334 + tmp_337 + tmp_342;
      real_t tmp_346 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_347 = 0.037198804536718075*tmp_58*(tmp_1*(tmp_344 - 1.0/4.0) + tmp_13*(tmp_346 - 1.0/4.0) + tmp_4*(tmp_345 - 1.0/4.0));
      real_t tmp_348 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_349 = tmp_348*tmp_6;
      real_t tmp_350 = tmp_26*tmp_348;
      real_t tmp_351 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_352 = tmp_28*tmp_351;
      real_t tmp_353 = tmp_35*tmp_351;
      real_t tmp_354 = tmp_348*tmp_37;
      real_t tmp_355 = tmp_351*tmp_39;
      real_t tmp_356 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_357 = tmp_356*tmp_41;
      real_t tmp_358 = tmp_356*tmp_48;
      real_t tmp_359 = tmp_356*tmp_50;
      real_t tmp_360 = tmp_354 + tmp_355 + tmp_359;
      real_t tmp_361 = tmp_350 + tmp_353 + tmp_358;
      real_t tmp_362 = tmp_349 + tmp_352 + tmp_357;
      real_t tmp_363 = 0.042507265838595799*tmp_58*(tmp_1*(tmp_360 - 1.0/4.0) + tmp_13*(tmp_362 - 1.0/4.0) + tmp_4*(tmp_361 - 1.0/4.0));
      real_t tmp_364 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_365 = tmp_364*tmp_6;
      real_t tmp_366 = tmp_26*tmp_364;
      real_t tmp_367 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_368 = tmp_28*tmp_367;
      real_t tmp_369 = tmp_35*tmp_367;
      real_t tmp_370 = tmp_364*tmp_37;
      real_t tmp_371 = tmp_367*tmp_39;
      real_t tmp_372 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_373 = tmp_372*tmp_41;
      real_t tmp_374 = tmp_372*tmp_48;
      real_t tmp_375 = tmp_372*tmp_50;
      real_t tmp_376 = tmp_370 + tmp_371 + tmp_375;
      real_t tmp_377 = tmp_366 + tmp_369 + tmp_374;
      real_t tmp_378 = tmp_365 + tmp_368 + tmp_373;
      real_t tmp_379 = 0.019202922745021479*tmp_58*(tmp_1*(tmp_376 - 1.0/4.0) + tmp_13*(tmp_378 - 1.0/4.0) + tmp_4*(tmp_377 - 1.0/4.0));
      real_t a_0_0 = tmp_107*(-tmp_101 - tmp_102 - tmp_103 - tmp_93 - tmp_94 - tmp_96 - tmp_97 - tmp_98 - tmp_99 + 1) + tmp_123*(-tmp_109 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_115 - tmp_117 - tmp_118 - tmp_119 + 1) + tmp_139*(-tmp_125 - tmp_126 - tmp_128 - tmp_129 - tmp_130 - tmp_131 - tmp_133 - tmp_134 - tmp_135 + 1) + tmp_155*(-tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1) + tmp_171*(-tmp_157 - tmp_158 - tmp_160 - tmp_161 - tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 + 1) + tmp_187*(-tmp_173 - tmp_174 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 + 1) + tmp_203*(-tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_197 - tmp_198 - tmp_199 + 1) + tmp_219*(-tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 - tmp_211 - tmp_213 - tmp_214 - tmp_215 + 1) + tmp_235*(-tmp_221 - tmp_222 - tmp_224 - tmp_225 - tmp_226 - tmp_227 - tmp_229 - tmp_230 - tmp_231 + 1) + tmp_251*(-tmp_237 - tmp_238 - tmp_240 - tmp_241 - tmp_242 - tmp_243 - tmp_245 - tmp_246 - tmp_247 + 1) + tmp_267*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_283*(-tmp_269 - tmp_270 - tmp_272 - tmp_273 - tmp_274 - tmp_275 - tmp_277 - tmp_278 - tmp_279 + 1) + tmp_299*(-tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1) + tmp_315*(-tmp_301 - tmp_302 - tmp_304 - tmp_305 - tmp_306 - tmp_307 - tmp_309 - tmp_310 - tmp_311 + 1) + tmp_331*(-tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1) + tmp_347*(-tmp_333 - tmp_334 - tmp_336 - tmp_337 - tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 + 1) + tmp_363*(-tmp_349 - tmp_350 - tmp_352 - tmp_353 - tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 + 1) + tmp_379*(-tmp_365 - tmp_366 - tmp_368 - tmp_369 - tmp_370 - tmp_371 - tmp_373 - tmp_374 - tmp_375 + 1) + tmp_59*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_75*(-tmp_61 - tmp_62 - tmp_64 - tmp_65 - tmp_66 - tmp_67 - tmp_69 - tmp_70 - tmp_71 + 1) + tmp_91*(-tmp_77 - tmp_78 - tmp_80 - tmp_81 - tmp_82 - tmp_83 - tmp_85 - tmp_86 - tmp_87 + 1);
      real_t a_1_0 = tmp_104*tmp_107 + tmp_120*tmp_123 + tmp_136*tmp_139 + tmp_152*tmp_155 + tmp_168*tmp_171 + tmp_184*tmp_187 + tmp_200*tmp_203 + tmp_216*tmp_219 + tmp_232*tmp_235 + tmp_248*tmp_251 + tmp_264*tmp_267 + tmp_280*tmp_283 + tmp_296*tmp_299 + tmp_312*tmp_315 + tmp_328*tmp_331 + tmp_344*tmp_347 + tmp_360*tmp_363 + tmp_376*tmp_379 + tmp_52*tmp_59 + tmp_72*tmp_75 + tmp_88*tmp_91;
      real_t a_2_0 = tmp_105*tmp_107 + tmp_121*tmp_123 + tmp_137*tmp_139 + tmp_153*tmp_155 + tmp_169*tmp_171 + tmp_185*tmp_187 + tmp_201*tmp_203 + tmp_217*tmp_219 + tmp_233*tmp_235 + tmp_249*tmp_251 + tmp_265*tmp_267 + tmp_281*tmp_283 + tmp_297*tmp_299 + tmp_313*tmp_315 + tmp_329*tmp_331 + tmp_345*tmp_347 + tmp_361*tmp_363 + tmp_377*tmp_379 + tmp_53*tmp_59 + tmp_73*tmp_75 + tmp_89*tmp_91;
      real_t a_3_0 = tmp_106*tmp_107 + tmp_122*tmp_123 + tmp_138*tmp_139 + tmp_154*tmp_155 + tmp_170*tmp_171 + tmp_186*tmp_187 + tmp_202*tmp_203 + tmp_218*tmp_219 + tmp_234*tmp_235 + tmp_250*tmp_251 + tmp_266*tmp_267 + tmp_282*tmp_283 + tmp_298*tmp_299 + tmp_314*tmp_315 + tmp_330*tmp_331 + tmp_346*tmp_347 + tmp_362*tmp_363 + tmp_378*tmp_379 + tmp_54*tmp_59 + tmp_74*tmp_75 + tmp_90*tmp_91;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }


};




class EGNIPGVectorLaplaceFormP1E_1 : public hyteg::dg::DGForm
{
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
      real_t tmp_28 = 5/tmp_27;
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
      real_t tmp_37 = 5/tmp_36;
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
      real_t tmp_20 = 0.59231721264047354*tmp_3*(tmp_19 - 1.0/3.0) + 0.59231721264047354*tmp_4*(tmp_18 - 1.0/3.0);
      real_t tmp_21 = tmp_5*(0.23076534494715845*tmp_6 + tmp_7);
      real_t tmp_22 = tmp_1*tmp_21;
      real_t tmp_23 = tmp_10*tmp_21;
      real_t tmp_24 = tmp_5*(0.23076534494715845*tmp_12 + tmp_13);
      real_t tmp_25 = tmp_24*tmp_3;
      real_t tmp_26 = tmp_16*tmp_24;
      real_t tmp_27 = tmp_23 + tmp_25;
      real_t tmp_28 = tmp_22 + tmp_26;
      real_t tmp_29 = 1.1965716762484155*tmp_3*(tmp_28 - 1.0/3.0) + 1.1965716762484155*tmp_4*(tmp_27 - 1.0/3.0);
      real_t tmp_30 = tmp_5*(0.5*tmp_6 + tmp_7);
      real_t tmp_31 = tmp_1*tmp_30;
      real_t tmp_32 = tmp_10*tmp_30;
      real_t tmp_33 = tmp_5*(0.5*tmp_12 + tmp_13);
      real_t tmp_34 = tmp_3*tmp_33;
      real_t tmp_35 = tmp_16*tmp_33;
      real_t tmp_36 = tmp_32 + tmp_34;
      real_t tmp_37 = tmp_31 + tmp_35;
      real_t tmp_38 = 1.4222222222222225*tmp_3*(tmp_37 - 1.0/3.0) + 1.4222222222222225*tmp_4*(tmp_36 - 1.0/3.0);
      real_t tmp_39 = tmp_5*(0.7692346550528415*tmp_6 + tmp_7);
      real_t tmp_40 = tmp_1*tmp_39;
      real_t tmp_41 = tmp_10*tmp_39;
      real_t tmp_42 = tmp_5*(0.7692346550528415*tmp_12 + tmp_13);
      real_t tmp_43 = tmp_3*tmp_42;
      real_t tmp_44 = tmp_16*tmp_42;
      real_t tmp_45 = tmp_41 + tmp_43;
      real_t tmp_46 = tmp_40 + tmp_44;
      real_t tmp_47 = 1.1965716762484155*tmp_3*(tmp_46 - 1.0/3.0) + 1.1965716762484155*tmp_4*(tmp_45 - 1.0/3.0);
      real_t tmp_48 = tmp_5*(0.95308992296933193*tmp_6 + tmp_7);
      real_t tmp_49 = tmp_1*tmp_48;
      real_t tmp_50 = tmp_10*tmp_48;
      real_t tmp_51 = tmp_5*(0.95308992296933193*tmp_12 + tmp_13);
      real_t tmp_52 = tmp_3*tmp_51;
      real_t tmp_53 = tmp_16*tmp_51;
      real_t tmp_54 = tmp_50 + tmp_52;
      real_t tmp_55 = tmp_49 + tmp_53;
      real_t tmp_56 = 0.59231721264047354*tmp_3*(tmp_55 - 1.0/3.0) + 0.59231721264047354*tmp_4*(tmp_54 - 1.0/3.0);
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
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
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
      real_t tmp_22 = tmp_11*tmp_19 + tmp_20*tmp_3 + tmp_21*tmp_6;
      real_t tmp_23 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_5);
      real_t tmp_24 = tmp_18*(tmp_1*tmp_9 - tmp_10*tmp_14);
      real_t tmp_25 = tmp_18*(tmp_13*tmp_14 - tmp_5*tmp_9);
      real_t tmp_26 = tmp_11*tmp_23 + tmp_24*tmp_3 + tmp_25*tmp_6;
      real_t tmp_27 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_6);
      real_t tmp_28 = tmp_18*(tmp_10*tmp_11 - tmp_6*tmp_9);
      real_t tmp_29 = tmp_18*(-tmp_11*tmp_13 + tmp_3*tmp_9);
      real_t tmp_30 = tmp_11*tmp_27 + tmp_28*tmp_3 + tmp_29*tmp_6;
      real_t tmp_31 = p_affine_0_0*p_affine_1_1;
      real_t tmp_32 = p_affine_0_0*p_affine_1_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_2;
      real_t tmp_34 = p_affine_0_1*p_affine_1_0;
      real_t tmp_35 = p_affine_0_1*p_affine_1_2;
      real_t tmp_36 = p_affine_2_2*p_affine_3_0;
      real_t tmp_37 = p_affine_0_2*p_affine_1_0;
      real_t tmp_38 = p_affine_0_2*p_affine_1_1;
      real_t tmp_39 = p_affine_2_0*p_affine_3_1;
      real_t tmp_40 = p_affine_2_2*p_affine_3_1;
      real_t tmp_41 = p_affine_2_0*p_affine_3_2;
      real_t tmp_42 = p_affine_2_1*p_affine_3_0;
      real_t tmp_43 = std::abs(p_affine_0_0*tmp_33 - p_affine_0_0*tmp_40 + p_affine_0_1*tmp_36 - p_affine_0_1*tmp_41 + p_affine_0_2*tmp_39 - p_affine_0_2*tmp_42 - p_affine_1_0*tmp_33 + p_affine_1_0*tmp_40 - p_affine_1_1*tmp_36 + p_affine_1_1*tmp_41 - p_affine_1_2*tmp_39 + p_affine_1_2*tmp_42 + p_affine_2_0*tmp_35 - p_affine_2_0*tmp_38 - p_affine_2_1*tmp_32 + p_affine_2_1*tmp_37 + p_affine_2_2*tmp_31 - p_affine_2_2*tmp_34 - p_affine_3_0*tmp_35 + p_affine_3_0*tmp_38 + p_affine_3_1*tmp_32 - p_affine_3_1*tmp_37 - p_affine_3_2*tmp_31 + p_affine_3_2*tmp_34);
      real_t tmp_44 = tmp_43*(tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25) + tmp_30*(-tmp_27 - tmp_28 - tmp_29));
      real_t tmp_45 = tmp_43*(tmp_21*tmp_22 + tmp_25*tmp_26 + tmp_29*tmp_30);
      real_t tmp_46 = tmp_43*(tmp_20*tmp_22 + tmp_24*tmp_26 + tmp_28*tmp_30);
      real_t tmp_47 = tmp_43*(tmp_19*tmp_22 + tmp_23*tmp_26 + tmp_27*tmp_30);
      real_t a_0_0 = 0.1666666666666668*tmp_44;
      real_t a_1_0 = 0.1666666666666668*tmp_45;
      real_t a_2_0 = 0.1666666666666668*tmp_46;
      real_t a_3_0 = 0.1666666666666668*tmp_47;
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

         real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_1_1 + tmp_0;
      real_t tmp_2 = -p_affine_8_2;
      real_t tmp_3 = p_affine_9_2 + tmp_2;
      real_t tmp_4 = p_affine_10_2 + tmp_2;
      real_t tmp_5 = -p_affine_0_2;
      real_t tmp_6 = p_affine_8_2 + tmp_5;
      real_t tmp_7 = 0.031405749086161582*tmp_3 + 0.93718850182767688*tmp_4 + tmp_6;
      real_t tmp_8 = -p_affine_0_0;
      real_t tmp_9 = p_affine_2_0 + tmp_8;
      real_t tmp_10 = p_affine_3_1 + tmp_0;
      real_t tmp_11 = p_affine_3_0 + tmp_8;
      real_t tmp_12 = p_affine_2_1 + tmp_0;
      real_t tmp_13 = p_affine_1_0 + tmp_8;
      real_t tmp_14 = p_affine_3_2 + tmp_5;
      real_t tmp_15 = tmp_12*tmp_14;
      real_t tmp_16 = p_affine_1_2 + tmp_5;
      real_t tmp_17 = tmp_10*tmp_16;
      real_t tmp_18 = p_affine_2_2 + tmp_5;
      real_t tmp_19 = tmp_1*tmp_18;
      real_t tmp_20 = tmp_10*tmp_18;
      real_t tmp_21 = tmp_1*tmp_14;
      real_t tmp_22 = tmp_12*tmp_16;
      real_t tmp_23 = 1.0 / (tmp_11*tmp_19 - tmp_11*tmp_22 + tmp_13*tmp_15 - tmp_13*tmp_20 + tmp_17*tmp_9 - tmp_21*tmp_9);
      real_t tmp_24 = tmp_23*(tmp_10*tmp_9 - tmp_11*tmp_12);
      real_t tmp_25 = tmp_24*tmp_7;
      real_t tmp_26 = -p_affine_8_1;
      real_t tmp_27 = p_affine_9_1 + tmp_26;
      real_t tmp_28 = p_affine_10_1 + tmp_26;
      real_t tmp_29 = p_affine_8_1 + tmp_0;
      real_t tmp_30 = 0.031405749086161582*tmp_27 + 0.93718850182767688*tmp_28 + tmp_29;
      real_t tmp_31 = tmp_23*(tmp_11*tmp_18 - tmp_14*tmp_9);
      real_t tmp_32 = tmp_30*tmp_31;
      real_t tmp_33 = -p_affine_8_0;
      real_t tmp_34 = p_affine_9_0 + tmp_33;
      real_t tmp_35 = p_affine_10_0 + tmp_33;
      real_t tmp_36 = p_affine_8_0 + tmp_8;
      real_t tmp_37 = 0.031405749086161582*tmp_34 + 0.93718850182767688*tmp_35 + tmp_36;
      real_t tmp_38 = tmp_23*(tmp_15 - tmp_20);
      real_t tmp_39 = tmp_37*tmp_38;
      real_t tmp_40 = tmp_25 + tmp_32 + tmp_39;
      real_t tmp_41 = tmp_23*(tmp_1*tmp_11 - tmp_10*tmp_13);
      real_t tmp_42 = tmp_41*tmp_7;
      real_t tmp_43 = tmp_23*(-tmp_11*tmp_16 + tmp_13*tmp_14);
      real_t tmp_44 = tmp_30*tmp_43;
      real_t tmp_45 = tmp_23*(tmp_17 - tmp_21);
      real_t tmp_46 = tmp_37*tmp_45;
      real_t tmp_47 = tmp_42 + tmp_44 + tmp_46;
      real_t tmp_48 = tmp_23*(-tmp_1*tmp_9 + tmp_12*tmp_13);
      real_t tmp_49 = tmp_48*tmp_7;
      real_t tmp_50 = tmp_23*(-tmp_13*tmp_18 + tmp_16*tmp_9);
      real_t tmp_51 = tmp_30*tmp_50;
      real_t tmp_52 = tmp_23*(tmp_19 - tmp_22);
      real_t tmp_53 = tmp_37*tmp_52;
      real_t tmp_54 = tmp_49 + tmp_51 + tmp_53;
      real_t tmp_55 = tmp_1*(tmp_40 - 1.0/4.0) + tmp_10*(tmp_54 - 1.0/4.0) + tmp_12*(tmp_47 - 1.0/4.0);
      real_t tmp_56 = 0.5*p_affine_13_0*(-tmp_38 - tmp_45 - tmp_52) + 0.5*p_affine_13_1*(-tmp_31 - tmp_43 - tmp_50) + 0.5*p_affine_13_2*(-tmp_24 - tmp_41 - tmp_48);
      real_t tmp_57 = -tmp_25 - tmp_32 - tmp_39 - tmp_42 - tmp_44 - tmp_46 - tmp_49 - tmp_51 - tmp_53 + 1;
      real_t tmp_58 = 0.5*p_affine_13_0*(tmp_1*tmp_38 + tmp_10*tmp_52 + tmp_12*tmp_45) + 0.5*p_affine_13_1*(tmp_1*tmp_31 + tmp_10*tmp_50 + tmp_12*tmp_43) + 0.5*p_affine_13_2*(tmp_1*tmp_24 + tmp_10*tmp_48 + tmp_12*tmp_41);
      real_t tmp_59 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_60 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_61 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_62 = (std::abs(tmp_28*tmp_60 - tmp_35*tmp_59)*std::abs(tmp_28*tmp_60 - tmp_35*tmp_59)) + (std::abs(tmp_28*tmp_61 - tmp_4*tmp_59)*std::abs(tmp_28*tmp_61 - tmp_4*tmp_59)) + (std::abs(tmp_35*tmp_61 - tmp_4*tmp_60)*std::abs(tmp_35*tmp_61 - tmp_4*tmp_60));
      real_t tmp_63 = 5.0*std::pow(tmp_62, -0.25);
      real_t tmp_64 = tmp_55*tmp_63;
      real_t tmp_65 = 1.0*std::pow(tmp_62, 1.0/2.0);
      real_t tmp_66 = 0.0068572537431980923*tmp_65;
      real_t tmp_67 = 0.19601935860219369*tmp_3 + 0.60796128279561268*tmp_4 + tmp_6;
      real_t tmp_68 = tmp_24*tmp_67;
      real_t tmp_69 = 0.19601935860219369*tmp_27 + 0.60796128279561268*tmp_28 + tmp_29;
      real_t tmp_70 = tmp_31*tmp_69;
      real_t tmp_71 = 0.19601935860219369*tmp_34 + 0.60796128279561268*tmp_35 + tmp_36;
      real_t tmp_72 = tmp_38*tmp_71;
      real_t tmp_73 = tmp_68 + tmp_70 + tmp_72;
      real_t tmp_74 = tmp_41*tmp_67;
      real_t tmp_75 = tmp_43*tmp_69;
      real_t tmp_76 = tmp_45*tmp_71;
      real_t tmp_77 = tmp_74 + tmp_75 + tmp_76;
      real_t tmp_78 = tmp_48*tmp_67;
      real_t tmp_79 = tmp_50*tmp_69;
      real_t tmp_80 = tmp_52*tmp_71;
      real_t tmp_81 = tmp_78 + tmp_79 + tmp_80;
      real_t tmp_82 = tmp_1*(tmp_73 - 1.0/4.0) + tmp_10*(tmp_81 - 1.0/4.0) + tmp_12*(tmp_77 - 1.0/4.0);
      real_t tmp_83 = -tmp_68 - tmp_70 - tmp_72 - tmp_74 - tmp_75 - tmp_76 - tmp_78 - tmp_79 - tmp_80 + 1;
      real_t tmp_84 = tmp_63*tmp_82;
      real_t tmp_85 = 0.037198804536718075*tmp_65;
      real_t tmp_86 = 0.37605877282253791*tmp_3 + 0.039308471900058539*tmp_4 + tmp_6;
      real_t tmp_87 = tmp_24*tmp_86;
      real_t tmp_88 = 0.37605877282253791*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29;
      real_t tmp_89 = tmp_31*tmp_88;
      real_t tmp_90 = 0.37605877282253791*tmp_34 + 0.039308471900058539*tmp_35 + tmp_36;
      real_t tmp_91 = tmp_38*tmp_90;
      real_t tmp_92 = tmp_87 + tmp_89 + tmp_91;
      real_t tmp_93 = tmp_41*tmp_86;
      real_t tmp_94 = tmp_43*tmp_88;
      real_t tmp_95 = tmp_45*tmp_90;
      real_t tmp_96 = tmp_93 + tmp_94 + tmp_95;
      real_t tmp_97 = tmp_48*tmp_86;
      real_t tmp_98 = tmp_50*tmp_88;
      real_t tmp_99 = tmp_52*tmp_90;
      real_t tmp_100 = tmp_97 + tmp_98 + tmp_99;
      real_t tmp_101 = tmp_1*(tmp_92 - 1.0/4.0) + tmp_10*(tmp_100 - 1.0/4.0) + tmp_12*(tmp_96 - 1.0/4.0);
      real_t tmp_102 = -tmp_87 - tmp_89 - tmp_91 - tmp_93 - tmp_94 - tmp_95 - tmp_97 - tmp_98 - tmp_99 + 1;
      real_t tmp_103 = tmp_101*tmp_63;
      real_t tmp_104 = 0.020848748529055869*tmp_65;
      real_t tmp_105 = 0.78764240869137092*tmp_3 + 0.1711304259088916*tmp_4 + tmp_6;
      real_t tmp_106 = tmp_105*tmp_24;
      real_t tmp_107 = 0.78764240869137092*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29;
      real_t tmp_108 = tmp_107*tmp_31;
      real_t tmp_109 = 0.78764240869137092*tmp_34 + 0.1711304259088916*tmp_35 + tmp_36;
      real_t tmp_110 = tmp_109*tmp_38;
      real_t tmp_111 = tmp_106 + tmp_108 + tmp_110;
      real_t tmp_112 = tmp_105*tmp_41;
      real_t tmp_113 = tmp_107*tmp_43;
      real_t tmp_114 = tmp_109*tmp_45;
      real_t tmp_115 = tmp_112 + tmp_113 + tmp_114;
      real_t tmp_116 = tmp_105*tmp_48;
      real_t tmp_117 = tmp_107*tmp_50;
      real_t tmp_118 = tmp_109*tmp_52;
      real_t tmp_119 = tmp_116 + tmp_117 + tmp_118;
      real_t tmp_120 = tmp_1*(tmp_111 - 1.0/4.0) + tmp_10*(tmp_119 - 1.0/4.0) + tmp_12*(tmp_115 - 1.0/4.0);
      real_t tmp_121 = -tmp_106 - tmp_108 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_116 - tmp_117 - tmp_118 + 1;
      real_t tmp_122 = tmp_120*tmp_63;
      real_t tmp_123 = 0.019202922745021479*tmp_65;
      real_t tmp_124 = 0.58463275527740355*tmp_3 + 0.37605877282253791*tmp_4 + tmp_6;
      real_t tmp_125 = tmp_124*tmp_24;
      real_t tmp_126 = 0.58463275527740355*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29;
      real_t tmp_127 = tmp_126*tmp_31;
      real_t tmp_128 = 0.58463275527740355*tmp_34 + 0.37605877282253791*tmp_35 + tmp_36;
      real_t tmp_129 = tmp_128*tmp_38;
      real_t tmp_130 = tmp_125 + tmp_127 + tmp_129;
      real_t tmp_131 = tmp_124*tmp_41;
      real_t tmp_132 = tmp_126*tmp_43;
      real_t tmp_133 = tmp_128*tmp_45;
      real_t tmp_134 = tmp_131 + tmp_132 + tmp_133;
      real_t tmp_135 = tmp_124*tmp_48;
      real_t tmp_136 = tmp_126*tmp_50;
      real_t tmp_137 = tmp_128*tmp_52;
      real_t tmp_138 = tmp_135 + tmp_136 + tmp_137;
      real_t tmp_139 = tmp_1*(tmp_130 - 1.0/4.0) + tmp_10*(tmp_138 - 1.0/4.0) + tmp_12*(tmp_134 - 1.0/4.0);
      real_t tmp_140 = -tmp_125 - tmp_127 - tmp_129 - tmp_131 - tmp_132 - tmp_133 - tmp_135 - tmp_136 - tmp_137 + 1;
      real_t tmp_141 = tmp_139*tmp_63;
      real_t tmp_142 = 0.020848748529055869*tmp_65;
      real_t tmp_143 = 0.041227165399737475*tmp_3 + 0.78764240869137092*tmp_4 + tmp_6;
      real_t tmp_144 = tmp_143*tmp_24;
      real_t tmp_145 = 0.041227165399737475*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29;
      real_t tmp_146 = tmp_145*tmp_31;
      real_t tmp_147 = 0.041227165399737475*tmp_34 + 0.78764240869137092*tmp_35 + tmp_36;
      real_t tmp_148 = tmp_147*tmp_38;
      real_t tmp_149 = tmp_144 + tmp_146 + tmp_148;
      real_t tmp_150 = tmp_143*tmp_41;
      real_t tmp_151 = tmp_145*tmp_43;
      real_t tmp_152 = tmp_147*tmp_45;
      real_t tmp_153 = tmp_150 + tmp_151 + tmp_152;
      real_t tmp_154 = tmp_143*tmp_48;
      real_t tmp_155 = tmp_145*tmp_50;
      real_t tmp_156 = tmp_147*tmp_52;
      real_t tmp_157 = tmp_154 + tmp_155 + tmp_156;
      real_t tmp_158 = tmp_1*(tmp_149 - 1.0/4.0) + tmp_10*(tmp_157 - 1.0/4.0) + tmp_12*(tmp_153 - 1.0/4.0);
      real_t tmp_159 = -tmp_144 - tmp_146 - tmp_148 - tmp_150 - tmp_151 - tmp_152 - tmp_154 - tmp_155 - tmp_156 + 1;
      real_t tmp_160 = tmp_158*tmp_63;
      real_t tmp_161 = 0.019202922745021479*tmp_65;
      real_t tmp_162 = 0.039308471900058539*tmp_3 + 0.58463275527740355*tmp_4 + tmp_6;
      real_t tmp_163 = tmp_162*tmp_24;
      real_t tmp_164 = 0.039308471900058539*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29;
      real_t tmp_165 = tmp_164*tmp_31;
      real_t tmp_166 = 0.039308471900058539*tmp_34 + 0.58463275527740355*tmp_35 + tmp_36;
      real_t tmp_167 = tmp_166*tmp_38;
      real_t tmp_168 = tmp_163 + tmp_165 + tmp_167;
      real_t tmp_169 = tmp_162*tmp_41;
      real_t tmp_170 = tmp_164*tmp_43;
      real_t tmp_171 = tmp_166*tmp_45;
      real_t tmp_172 = tmp_169 + tmp_170 + tmp_171;
      real_t tmp_173 = tmp_162*tmp_48;
      real_t tmp_174 = tmp_164*tmp_50;
      real_t tmp_175 = tmp_166*tmp_52;
      real_t tmp_176 = tmp_173 + tmp_174 + tmp_175;
      real_t tmp_177 = tmp_1*(tmp_168 - 1.0/4.0) + tmp_10*(tmp_176 - 1.0/4.0) + tmp_12*(tmp_172 - 1.0/4.0);
      real_t tmp_178 = -tmp_163 - tmp_165 - tmp_167 - tmp_169 - tmp_170 - tmp_171 - tmp_173 - tmp_174 - tmp_175 + 1;
      real_t tmp_179 = tmp_177*tmp_63;
      real_t tmp_180 = 0.020848748529055869*tmp_65;
      real_t tmp_181 = 0.78764240869137092*tmp_3 + 0.041227165399737475*tmp_4 + tmp_6;
      real_t tmp_182 = tmp_181*tmp_24;
      real_t tmp_183 = 0.78764240869137092*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29;
      real_t tmp_184 = tmp_183*tmp_31;
      real_t tmp_185 = 0.78764240869137092*tmp_34 + 0.041227165399737475*tmp_35 + tmp_36;
      real_t tmp_186 = tmp_185*tmp_38;
      real_t tmp_187 = tmp_182 + tmp_184 + tmp_186;
      real_t tmp_188 = tmp_181*tmp_41;
      real_t tmp_189 = tmp_183*tmp_43;
      real_t tmp_190 = tmp_185*tmp_45;
      real_t tmp_191 = tmp_188 + tmp_189 + tmp_190;
      real_t tmp_192 = tmp_181*tmp_48;
      real_t tmp_193 = tmp_183*tmp_50;
      real_t tmp_194 = tmp_185*tmp_52;
      real_t tmp_195 = tmp_192 + tmp_193 + tmp_194;
      real_t tmp_196 = tmp_1*(tmp_187 - 1.0/4.0) + tmp_10*(tmp_195 - 1.0/4.0) + tmp_12*(tmp_191 - 1.0/4.0);
      real_t tmp_197 = -tmp_182 - tmp_184 - tmp_186 - tmp_188 - tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 + 1;
      real_t tmp_198 = tmp_196*tmp_63;
      real_t tmp_199 = 0.019202922745021479*tmp_65;
      real_t tmp_200 = 0.58463275527740355*tmp_3 + 0.039308471900058539*tmp_4 + tmp_6;
      real_t tmp_201 = tmp_200*tmp_24;
      real_t tmp_202 = 0.58463275527740355*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29;
      real_t tmp_203 = tmp_202*tmp_31;
      real_t tmp_204 = 0.58463275527740355*tmp_34 + 0.039308471900058539*tmp_35 + tmp_36;
      real_t tmp_205 = tmp_204*tmp_38;
      real_t tmp_206 = tmp_201 + tmp_203 + tmp_205;
      real_t tmp_207 = tmp_200*tmp_41;
      real_t tmp_208 = tmp_202*tmp_43;
      real_t tmp_209 = tmp_204*tmp_45;
      real_t tmp_210 = tmp_207 + tmp_208 + tmp_209;
      real_t tmp_211 = tmp_200*tmp_48;
      real_t tmp_212 = tmp_202*tmp_50;
      real_t tmp_213 = tmp_204*tmp_52;
      real_t tmp_214 = tmp_211 + tmp_212 + tmp_213;
      real_t tmp_215 = tmp_1*(tmp_206 - 1.0/4.0) + tmp_10*(tmp_214 - 1.0/4.0) + tmp_12*(tmp_210 - 1.0/4.0);
      real_t tmp_216 = -tmp_201 - tmp_203 - tmp_205 - tmp_207 - tmp_208 - tmp_209 - tmp_211 - tmp_212 - tmp_213 + 1;
      real_t tmp_217 = tmp_215*tmp_63;
      real_t tmp_218 = 0.020848748529055869*tmp_65;
      real_t tmp_219 = 0.1711304259088916*tmp_3 + 0.78764240869137092*tmp_4 + tmp_6;
      real_t tmp_220 = tmp_219*tmp_24;
      real_t tmp_221 = 0.1711304259088916*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29;
      real_t tmp_222 = tmp_221*tmp_31;
      real_t tmp_223 = 0.1711304259088916*tmp_34 + 0.78764240869137092*tmp_35 + tmp_36;
      real_t tmp_224 = tmp_223*tmp_38;
      real_t tmp_225 = tmp_220 + tmp_222 + tmp_224;
      real_t tmp_226 = tmp_219*tmp_41;
      real_t tmp_227 = tmp_221*tmp_43;
      real_t tmp_228 = tmp_223*tmp_45;
      real_t tmp_229 = tmp_226 + tmp_227 + tmp_228;
      real_t tmp_230 = tmp_219*tmp_48;
      real_t tmp_231 = tmp_221*tmp_50;
      real_t tmp_232 = tmp_223*tmp_52;
      real_t tmp_233 = tmp_230 + tmp_231 + tmp_232;
      real_t tmp_234 = tmp_1*(tmp_225 - 1.0/4.0) + tmp_10*(tmp_233 - 1.0/4.0) + tmp_12*(tmp_229 - 1.0/4.0);
      real_t tmp_235 = -tmp_220 - tmp_222 - tmp_224 - tmp_226 - tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 + 1;
      real_t tmp_236 = tmp_234*tmp_63;
      real_t tmp_237 = 0.019202922745021479*tmp_65;
      real_t tmp_238 = 0.37605877282253791*tmp_3 + 0.58463275527740355*tmp_4 + tmp_6;
      real_t tmp_239 = tmp_238*tmp_24;
      real_t tmp_240 = 0.37605877282253791*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29;
      real_t tmp_241 = tmp_240*tmp_31;
      real_t tmp_242 = 0.37605877282253791*tmp_34 + 0.58463275527740355*tmp_35 + tmp_36;
      real_t tmp_243 = tmp_242*tmp_38;
      real_t tmp_244 = tmp_239 + tmp_241 + tmp_243;
      real_t tmp_245 = tmp_238*tmp_41;
      real_t tmp_246 = tmp_240*tmp_43;
      real_t tmp_247 = tmp_242*tmp_45;
      real_t tmp_248 = tmp_245 + tmp_246 + tmp_247;
      real_t tmp_249 = tmp_238*tmp_48;
      real_t tmp_250 = tmp_240*tmp_50;
      real_t tmp_251 = tmp_242*tmp_52;
      real_t tmp_252 = tmp_249 + tmp_250 + tmp_251;
      real_t tmp_253 = tmp_1*(tmp_244 - 1.0/4.0) + tmp_10*(tmp_252 - 1.0/4.0) + tmp_12*(tmp_248 - 1.0/4.0);
      real_t tmp_254 = -tmp_239 - tmp_241 - tmp_243 - tmp_245 - tmp_246 - tmp_247 - tmp_249 - tmp_250 - tmp_251 + 1;
      real_t tmp_255 = tmp_253*tmp_63;
      real_t tmp_256 = 0.020848748529055869*tmp_65;
      real_t tmp_257 = 0.041227165399737475*tmp_3 + 0.1711304259088916*tmp_4 + tmp_6;
      real_t tmp_258 = tmp_24*tmp_257;
      real_t tmp_259 = 0.041227165399737475*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29;
      real_t tmp_260 = tmp_259*tmp_31;
      real_t tmp_261 = 0.041227165399737475*tmp_34 + 0.1711304259088916*tmp_35 + tmp_36;
      real_t tmp_262 = tmp_261*tmp_38;
      real_t tmp_263 = tmp_258 + tmp_260 + tmp_262;
      real_t tmp_264 = tmp_257*tmp_41;
      real_t tmp_265 = tmp_259*tmp_43;
      real_t tmp_266 = tmp_261*tmp_45;
      real_t tmp_267 = tmp_264 + tmp_265 + tmp_266;
      real_t tmp_268 = tmp_257*tmp_48;
      real_t tmp_269 = tmp_259*tmp_50;
      real_t tmp_270 = tmp_261*tmp_52;
      real_t tmp_271 = tmp_268 + tmp_269 + tmp_270;
      real_t tmp_272 = tmp_1*(tmp_263 - 1.0/4.0) + tmp_10*(tmp_271 - 1.0/4.0) + tmp_12*(tmp_267 - 1.0/4.0);
      real_t tmp_273 = -tmp_258 - tmp_260 - tmp_262 - tmp_264 - tmp_265 - tmp_266 - tmp_268 - tmp_269 - tmp_270 + 1;
      real_t tmp_274 = tmp_272*tmp_63;
      real_t tmp_275 = 0.019202922745021479*tmp_65;
      real_t tmp_276 = 0.40446199974765351*tmp_3 + 0.19107600050469298*tmp_4 + tmp_6;
      real_t tmp_277 = tmp_24*tmp_276;
      real_t tmp_278 = 0.40446199974765351*tmp_27 + 0.19107600050469298*tmp_28 + tmp_29;
      real_t tmp_279 = tmp_278*tmp_31;
      real_t tmp_280 = 0.40446199974765351*tmp_34 + 0.19107600050469298*tmp_35 + tmp_36;
      real_t tmp_281 = tmp_280*tmp_38;
      real_t tmp_282 = tmp_277 + tmp_279 + tmp_281;
      real_t tmp_283 = tmp_276*tmp_41;
      real_t tmp_284 = tmp_278*tmp_43;
      real_t tmp_285 = tmp_280*tmp_45;
      real_t tmp_286 = tmp_283 + tmp_284 + tmp_285;
      real_t tmp_287 = tmp_276*tmp_48;
      real_t tmp_288 = tmp_278*tmp_50;
      real_t tmp_289 = tmp_280*tmp_52;
      real_t tmp_290 = tmp_287 + tmp_288 + tmp_289;
      real_t tmp_291 = tmp_1*(tmp_282 - 1.0/4.0) + tmp_10*(tmp_290 - 1.0/4.0) + tmp_12*(tmp_286 - 1.0/4.0);
      real_t tmp_292 = -tmp_277 - tmp_279 - tmp_281 - tmp_283 - tmp_284 - tmp_285 - tmp_287 - tmp_288 - tmp_289 + 1;
      real_t tmp_293 = tmp_291*tmp_63;
      real_t tmp_294 = 0.042507265838595799*tmp_65;
      real_t tmp_295 = 0.039308471900058539*tmp_3 + 0.37605877282253791*tmp_4 + tmp_6;
      real_t tmp_296 = tmp_24*tmp_295;
      real_t tmp_297 = 0.039308471900058539*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29;
      real_t tmp_298 = tmp_297*tmp_31;
      real_t tmp_299 = 0.039308471900058539*tmp_34 + 0.37605877282253791*tmp_35 + tmp_36;
      real_t tmp_300 = tmp_299*tmp_38;
      real_t tmp_301 = tmp_296 + tmp_298 + tmp_300;
      real_t tmp_302 = tmp_295*tmp_41;
      real_t tmp_303 = tmp_297*tmp_43;
      real_t tmp_304 = tmp_299*tmp_45;
      real_t tmp_305 = tmp_302 + tmp_303 + tmp_304;
      real_t tmp_306 = tmp_295*tmp_48;
      real_t tmp_307 = tmp_297*tmp_50;
      real_t tmp_308 = tmp_299*tmp_52;
      real_t tmp_309 = tmp_306 + tmp_307 + tmp_308;
      real_t tmp_310 = tmp_1*(tmp_301 - 1.0/4.0) + tmp_10*(tmp_309 - 1.0/4.0) + tmp_12*(tmp_305 - 1.0/4.0);
      real_t tmp_311 = -tmp_296 - tmp_298 - tmp_300 - tmp_302 - tmp_303 - tmp_304 - tmp_306 - tmp_307 - tmp_308 + 1;
      real_t tmp_312 = tmp_310*tmp_63;
      real_t tmp_313 = 0.020848748529055869*tmp_65;
      real_t tmp_314 = 0.93718850182767688*tmp_3 + 0.031405749086161582*tmp_4 + tmp_6;
      real_t tmp_315 = tmp_24*tmp_314;
      real_t tmp_316 = 0.93718850182767688*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29;
      real_t tmp_317 = tmp_31*tmp_316;
      real_t tmp_318 = 0.93718850182767688*tmp_34 + 0.031405749086161582*tmp_35 + tmp_36;
      real_t tmp_319 = tmp_318*tmp_38;
      real_t tmp_320 = tmp_315 + tmp_317 + tmp_319;
      real_t tmp_321 = tmp_314*tmp_41;
      real_t tmp_322 = tmp_316*tmp_43;
      real_t tmp_323 = tmp_318*tmp_45;
      real_t tmp_324 = tmp_321 + tmp_322 + tmp_323;
      real_t tmp_325 = tmp_314*tmp_48;
      real_t tmp_326 = tmp_316*tmp_50;
      real_t tmp_327 = tmp_318*tmp_52;
      real_t tmp_328 = tmp_325 + tmp_326 + tmp_327;
      real_t tmp_329 = tmp_1*(tmp_320 - 1.0/4.0) + tmp_10*(tmp_328 - 1.0/4.0) + tmp_12*(tmp_324 - 1.0/4.0);
      real_t tmp_330 = -tmp_315 - tmp_317 - tmp_319 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1;
      real_t tmp_331 = tmp_329*tmp_63;
      real_t tmp_332 = 0.0068572537431980923*tmp_65;
      real_t tmp_333 = 0.60796128279561268*tmp_3 + 0.19601935860219369*tmp_4 + tmp_6;
      real_t tmp_334 = tmp_24*tmp_333;
      real_t tmp_335 = 0.60796128279561268*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29;
      real_t tmp_336 = tmp_31*tmp_335;
      real_t tmp_337 = 0.60796128279561268*tmp_34 + 0.19601935860219369*tmp_35 + tmp_36;
      real_t tmp_338 = tmp_337*tmp_38;
      real_t tmp_339 = tmp_334 + tmp_336 + tmp_338;
      real_t tmp_340 = tmp_333*tmp_41;
      real_t tmp_341 = tmp_335*tmp_43;
      real_t tmp_342 = tmp_337*tmp_45;
      real_t tmp_343 = tmp_340 + tmp_341 + tmp_342;
      real_t tmp_344 = tmp_333*tmp_48;
      real_t tmp_345 = tmp_335*tmp_50;
      real_t tmp_346 = tmp_337*tmp_52;
      real_t tmp_347 = tmp_344 + tmp_345 + tmp_346;
      real_t tmp_348 = tmp_1*(tmp_339 - 1.0/4.0) + tmp_10*(tmp_347 - 1.0/4.0) + tmp_12*(tmp_343 - 1.0/4.0);
      real_t tmp_349 = -tmp_334 - tmp_336 - tmp_338 - tmp_340 - tmp_341 - tmp_342 - tmp_344 - tmp_345 - tmp_346 + 1;
      real_t tmp_350 = tmp_348*tmp_63;
      real_t tmp_351 = 0.037198804536718075*tmp_65;
      real_t tmp_352 = 0.19107600050469298*tmp_3 + 0.40446199974765351*tmp_4 + tmp_6;
      real_t tmp_353 = tmp_24*tmp_352;
      real_t tmp_354 = 0.19107600050469298*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29;
      real_t tmp_355 = tmp_31*tmp_354;
      real_t tmp_356 = 0.19107600050469298*tmp_34 + 0.40446199974765351*tmp_35 + tmp_36;
      real_t tmp_357 = tmp_356*tmp_38;
      real_t tmp_358 = tmp_353 + tmp_355 + tmp_357;
      real_t tmp_359 = tmp_352*tmp_41;
      real_t tmp_360 = tmp_354*tmp_43;
      real_t tmp_361 = tmp_356*tmp_45;
      real_t tmp_362 = tmp_359 + tmp_360 + tmp_361;
      real_t tmp_363 = tmp_352*tmp_48;
      real_t tmp_364 = tmp_354*tmp_50;
      real_t tmp_365 = tmp_356*tmp_52;
      real_t tmp_366 = tmp_363 + tmp_364 + tmp_365;
      real_t tmp_367 = tmp_1*(tmp_358 - 1.0/4.0) + tmp_10*(tmp_366 - 1.0/4.0) + tmp_12*(tmp_362 - 1.0/4.0);
      real_t tmp_368 = -tmp_353 - tmp_355 - tmp_357 - tmp_359 - tmp_360 - tmp_361 - tmp_363 - tmp_364 - tmp_365 + 1;
      real_t tmp_369 = tmp_367*tmp_63;
      real_t tmp_370 = 0.042507265838595799*tmp_65;
      real_t tmp_371 = 0.031405749086161582*tmp_3 + 0.031405749086161582*tmp_4 + tmp_6;
      real_t tmp_372 = tmp_24*tmp_371;
      real_t tmp_373 = 0.031405749086161582*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29;
      real_t tmp_374 = tmp_31*tmp_373;
      real_t tmp_375 = 0.031405749086161582*tmp_34 + 0.031405749086161582*tmp_35 + tmp_36;
      real_t tmp_376 = tmp_375*tmp_38;
      real_t tmp_377 = tmp_372 + tmp_374 + tmp_376;
      real_t tmp_378 = tmp_371*tmp_41;
      real_t tmp_379 = tmp_373*tmp_43;
      real_t tmp_380 = tmp_375*tmp_45;
      real_t tmp_381 = tmp_378 + tmp_379 + tmp_380;
      real_t tmp_382 = tmp_371*tmp_48;
      real_t tmp_383 = tmp_373*tmp_50;
      real_t tmp_384 = tmp_375*tmp_52;
      real_t tmp_385 = tmp_382 + tmp_383 + tmp_384;
      real_t tmp_386 = tmp_1*(tmp_377 - 1.0/4.0) + tmp_10*(tmp_385 - 1.0/4.0) + tmp_12*(tmp_381 - 1.0/4.0);
      real_t tmp_387 = -tmp_372 - tmp_374 - tmp_376 - tmp_378 - tmp_379 - tmp_380 - tmp_382 - tmp_383 - tmp_384 + 1;
      real_t tmp_388 = tmp_386*tmp_63;
      real_t tmp_389 = 0.0068572537431980923*tmp_65;
      real_t tmp_390 = 0.19601935860219369*tmp_3 + 0.19601935860219369*tmp_4 + tmp_6;
      real_t tmp_391 = tmp_24*tmp_390;
      real_t tmp_392 = 0.19601935860219369*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29;
      real_t tmp_393 = tmp_31*tmp_392;
      real_t tmp_394 = 0.19601935860219369*tmp_34 + 0.19601935860219369*tmp_35 + tmp_36;
      real_t tmp_395 = tmp_38*tmp_394;
      real_t tmp_396 = tmp_391 + tmp_393 + tmp_395;
      real_t tmp_397 = tmp_390*tmp_41;
      real_t tmp_398 = tmp_392*tmp_43;
      real_t tmp_399 = tmp_394*tmp_45;
      real_t tmp_400 = tmp_397 + tmp_398 + tmp_399;
      real_t tmp_401 = tmp_390*tmp_48;
      real_t tmp_402 = tmp_392*tmp_50;
      real_t tmp_403 = tmp_394*tmp_52;
      real_t tmp_404 = tmp_401 + tmp_402 + tmp_403;
      real_t tmp_405 = tmp_1*(tmp_396 - 1.0/4.0) + tmp_10*(tmp_404 - 1.0/4.0) + tmp_12*(tmp_400 - 1.0/4.0);
      real_t tmp_406 = -tmp_391 - tmp_393 - tmp_395 - tmp_397 - tmp_398 - tmp_399 - tmp_401 - tmp_402 - tmp_403 + 1;
      real_t tmp_407 = tmp_405*tmp_63;
      real_t tmp_408 = 0.037198804536718075*tmp_65;
      real_t tmp_409 = 0.40446199974765351*tmp_3 + 0.40446199974765351*tmp_4 + tmp_6;
      real_t tmp_410 = tmp_24*tmp_409;
      real_t tmp_411 = 0.40446199974765351*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29;
      real_t tmp_412 = tmp_31*tmp_411;
      real_t tmp_413 = 0.40446199974765351*tmp_34 + 0.40446199974765351*tmp_35 + tmp_36;
      real_t tmp_414 = tmp_38*tmp_413;
      real_t tmp_415 = tmp_410 + tmp_412 + tmp_414;
      real_t tmp_416 = tmp_409*tmp_41;
      real_t tmp_417 = tmp_411*tmp_43;
      real_t tmp_418 = tmp_413*tmp_45;
      real_t tmp_419 = tmp_416 + tmp_417 + tmp_418;
      real_t tmp_420 = tmp_409*tmp_48;
      real_t tmp_421 = tmp_411*tmp_50;
      real_t tmp_422 = tmp_413*tmp_52;
      real_t tmp_423 = tmp_420 + tmp_421 + tmp_422;
      real_t tmp_424 = tmp_1*(tmp_415 - 1.0/4.0) + tmp_10*(tmp_423 - 1.0/4.0) + tmp_12*(tmp_419 - 1.0/4.0);
      real_t tmp_425 = -tmp_410 - tmp_412 - tmp_414 - tmp_416 - tmp_417 - tmp_418 - tmp_420 - tmp_421 - tmp_422 + 1;
      real_t tmp_426 = tmp_424*tmp_63;
      real_t tmp_427 = 0.042507265838595799*tmp_65;
      real_t tmp_428 = 0.1711304259088916*tmp_3 + 0.041227165399737475*tmp_4 + tmp_6;
      real_t tmp_429 = tmp_24*tmp_428;
      real_t tmp_430 = 0.1711304259088916*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29;
      real_t tmp_431 = tmp_31*tmp_430;
      real_t tmp_432 = 0.1711304259088916*tmp_34 + 0.041227165399737475*tmp_35 + tmp_36;
      real_t tmp_433 = tmp_38*tmp_432;
      real_t tmp_434 = tmp_429 + tmp_431 + tmp_433;
      real_t tmp_435 = tmp_41*tmp_428;
      real_t tmp_436 = tmp_43*tmp_430;
      real_t tmp_437 = tmp_432*tmp_45;
      real_t tmp_438 = tmp_435 + tmp_436 + tmp_437;
      real_t tmp_439 = tmp_428*tmp_48;
      real_t tmp_440 = tmp_430*tmp_50;
      real_t tmp_441 = tmp_432*tmp_52;
      real_t tmp_442 = tmp_439 + tmp_440 + tmp_441;
      real_t tmp_443 = tmp_1*(tmp_434 - 1.0/4.0) + tmp_10*(tmp_442 - 1.0/4.0) + tmp_12*(tmp_438 - 1.0/4.0);
      real_t tmp_444 = -tmp_429 - tmp_431 - tmp_433 - tmp_435 - tmp_436 - tmp_437 - tmp_439 - tmp_440 - tmp_441 + 1;
      real_t tmp_445 = tmp_443*tmp_63;
      real_t tmp_446 = 0.019202922745021479*tmp_65;
      real_t tmp_447 = 0.5*p_affine_13_0*tmp_38 + 0.5*p_affine_13_1*tmp_31 + 0.5*p_affine_13_2*tmp_24;
      real_t tmp_448 = 0.5*p_affine_13_0*tmp_45 + 0.5*p_affine_13_1*tmp_43 + 0.5*p_affine_13_2*tmp_41;
      real_t tmp_449 = 0.5*p_affine_13_0*tmp_52 + 0.5*p_affine_13_1*tmp_50 + 0.5*p_affine_13_2*tmp_48;
      real_t a_0_0 = tmp_104*(-tmp_101*tmp_56 + tmp_102*tmp_103 - tmp_102*tmp_58) + tmp_123*(-tmp_120*tmp_56 + tmp_121*tmp_122 - tmp_121*tmp_58) + tmp_142*(-tmp_139*tmp_56 + tmp_140*tmp_141 - tmp_140*tmp_58) + tmp_161*(-tmp_158*tmp_56 + tmp_159*tmp_160 - tmp_159*tmp_58) + tmp_180*(-tmp_177*tmp_56 + tmp_178*tmp_179 - tmp_178*tmp_58) + tmp_199*(-tmp_196*tmp_56 + tmp_197*tmp_198 - tmp_197*tmp_58) + tmp_218*(-tmp_215*tmp_56 + tmp_216*tmp_217 - tmp_216*tmp_58) + tmp_237*(-tmp_234*tmp_56 + tmp_235*tmp_236 - tmp_235*tmp_58) + tmp_256*(-tmp_253*tmp_56 + tmp_254*tmp_255 - tmp_254*tmp_58) + tmp_275*(-tmp_272*tmp_56 + tmp_273*tmp_274 - tmp_273*tmp_58) + tmp_294*(-tmp_291*tmp_56 + tmp_292*tmp_293 - tmp_292*tmp_58) + tmp_313*(-tmp_310*tmp_56 + tmp_311*tmp_312 - tmp_311*tmp_58) + tmp_332*(-tmp_329*tmp_56 + tmp_330*tmp_331 - tmp_330*tmp_58) + tmp_351*(-tmp_348*tmp_56 + tmp_349*tmp_350 - tmp_349*tmp_58) + tmp_370*(-tmp_367*tmp_56 + tmp_368*tmp_369 - tmp_368*tmp_58) + tmp_389*(-tmp_386*tmp_56 + tmp_387*tmp_388 - tmp_387*tmp_58) + tmp_408*(-tmp_405*tmp_56 + tmp_406*tmp_407 - tmp_406*tmp_58) + tmp_427*(-tmp_424*tmp_56 + tmp_425*tmp_426 - tmp_425*tmp_58) + tmp_446*(-tmp_443*tmp_56 + tmp_444*tmp_445 - tmp_444*tmp_58) + tmp_66*(-tmp_55*tmp_56 - tmp_57*tmp_58 + tmp_57*tmp_64) + tmp_85*(-tmp_56*tmp_82 - tmp_58*tmp_83 + tmp_83*tmp_84);
      real_t a_1_0 = tmp_104*(-tmp_101*tmp_447 + tmp_103*tmp_92 - tmp_58*tmp_92) + tmp_123*(tmp_111*tmp_122 - tmp_111*tmp_58 - tmp_120*tmp_447) + tmp_142*(tmp_130*tmp_141 - tmp_130*tmp_58 - tmp_139*tmp_447) + tmp_161*(tmp_149*tmp_160 - tmp_149*tmp_58 - tmp_158*tmp_447) + tmp_180*(tmp_168*tmp_179 - tmp_168*tmp_58 - tmp_177*tmp_447) + tmp_199*(tmp_187*tmp_198 - tmp_187*tmp_58 - tmp_196*tmp_447) + tmp_218*(tmp_206*tmp_217 - tmp_206*tmp_58 - tmp_215*tmp_447) + tmp_237*(tmp_225*tmp_236 - tmp_225*tmp_58 - tmp_234*tmp_447) + tmp_256*(tmp_244*tmp_255 - tmp_244*tmp_58 - tmp_253*tmp_447) + tmp_275*(tmp_263*tmp_274 - tmp_263*tmp_58 - tmp_272*tmp_447) + tmp_294*(tmp_282*tmp_293 - tmp_282*tmp_58 - tmp_291*tmp_447) + tmp_313*(tmp_301*tmp_312 - tmp_301*tmp_58 - tmp_310*tmp_447) + tmp_332*(tmp_320*tmp_331 - tmp_320*tmp_58 - tmp_329*tmp_447) + tmp_351*(tmp_339*tmp_350 - tmp_339*tmp_58 - tmp_348*tmp_447) + tmp_370*(tmp_358*tmp_369 - tmp_358*tmp_58 - tmp_367*tmp_447) + tmp_389*(tmp_377*tmp_388 - tmp_377*tmp_58 - tmp_386*tmp_447) + tmp_408*(tmp_396*tmp_407 - tmp_396*tmp_58 - tmp_405*tmp_447) + tmp_427*(tmp_415*tmp_426 - tmp_415*tmp_58 - tmp_424*tmp_447) + tmp_446*(tmp_434*tmp_445 - tmp_434*tmp_58 - tmp_443*tmp_447) + tmp_66*(-tmp_40*tmp_58 + tmp_40*tmp_64 - tmp_447*tmp_55) + tmp_85*(-tmp_447*tmp_82 - tmp_58*tmp_73 + tmp_73*tmp_84);
      real_t a_2_0 = tmp_104*(-tmp_101*tmp_448 + tmp_103*tmp_96 - tmp_58*tmp_96) + tmp_123*(tmp_115*tmp_122 - tmp_115*tmp_58 - tmp_120*tmp_448) + tmp_142*(tmp_134*tmp_141 - tmp_134*tmp_58 - tmp_139*tmp_448) + tmp_161*(tmp_153*tmp_160 - tmp_153*tmp_58 - tmp_158*tmp_448) + tmp_180*(tmp_172*tmp_179 - tmp_172*tmp_58 - tmp_177*tmp_448) + tmp_199*(tmp_191*tmp_198 - tmp_191*tmp_58 - tmp_196*tmp_448) + tmp_218*(tmp_210*tmp_217 - tmp_210*tmp_58 - tmp_215*tmp_448) + tmp_237*(tmp_229*tmp_236 - tmp_229*tmp_58 - tmp_234*tmp_448) + tmp_256*(tmp_248*tmp_255 - tmp_248*tmp_58 - tmp_253*tmp_448) + tmp_275*(tmp_267*tmp_274 - tmp_267*tmp_58 - tmp_272*tmp_448) + tmp_294*(tmp_286*tmp_293 - tmp_286*tmp_58 - tmp_291*tmp_448) + tmp_313*(tmp_305*tmp_312 - tmp_305*tmp_58 - tmp_310*tmp_448) + tmp_332*(tmp_324*tmp_331 - tmp_324*tmp_58 - tmp_329*tmp_448) + tmp_351*(tmp_343*tmp_350 - tmp_343*tmp_58 - tmp_348*tmp_448) + tmp_370*(tmp_362*tmp_369 - tmp_362*tmp_58 - tmp_367*tmp_448) + tmp_389*(tmp_381*tmp_388 - tmp_381*tmp_58 - tmp_386*tmp_448) + tmp_408*(tmp_400*tmp_407 - tmp_400*tmp_58 - tmp_405*tmp_448) + tmp_427*(tmp_419*tmp_426 - tmp_419*tmp_58 - tmp_424*tmp_448) + tmp_446*(tmp_438*tmp_445 - tmp_438*tmp_58 - tmp_443*tmp_448) + tmp_66*(-tmp_448*tmp_55 - tmp_47*tmp_58 + tmp_47*tmp_64) + tmp_85*(-tmp_448*tmp_82 - tmp_58*tmp_77 + tmp_77*tmp_84);
      real_t a_3_0 = tmp_104*(tmp_100*tmp_103 - tmp_100*tmp_58 - tmp_101*tmp_449) + tmp_123*(tmp_119*tmp_122 - tmp_119*tmp_58 - tmp_120*tmp_449) + tmp_142*(tmp_138*tmp_141 - tmp_138*tmp_58 - tmp_139*tmp_449) + tmp_161*(tmp_157*tmp_160 - tmp_157*tmp_58 - tmp_158*tmp_449) + tmp_180*(tmp_176*tmp_179 - tmp_176*tmp_58 - tmp_177*tmp_449) + tmp_199*(tmp_195*tmp_198 - tmp_195*tmp_58 - tmp_196*tmp_449) + tmp_218*(tmp_214*tmp_217 - tmp_214*tmp_58 - tmp_215*tmp_449) + tmp_237*(tmp_233*tmp_236 - tmp_233*tmp_58 - tmp_234*tmp_449) + tmp_256*(tmp_252*tmp_255 - tmp_252*tmp_58 - tmp_253*tmp_449) + tmp_275*(tmp_271*tmp_274 - tmp_271*tmp_58 - tmp_272*tmp_449) + tmp_294*(tmp_290*tmp_293 - tmp_290*tmp_58 - tmp_291*tmp_449) + tmp_313*(tmp_309*tmp_312 - tmp_309*tmp_58 - tmp_310*tmp_449) + tmp_332*(tmp_328*tmp_331 - tmp_328*tmp_58 - tmp_329*tmp_449) + tmp_351*(tmp_347*tmp_350 - tmp_347*tmp_58 - tmp_348*tmp_449) + tmp_370*(tmp_366*tmp_369 - tmp_366*tmp_58 - tmp_367*tmp_449) + tmp_389*(tmp_385*tmp_388 - tmp_385*tmp_58 - tmp_386*tmp_449) + tmp_408*(tmp_404*tmp_407 - tmp_404*tmp_58 - tmp_405*tmp_449) + tmp_427*(tmp_423*tmp_426 - tmp_423*tmp_58 - tmp_424*tmp_449) + tmp_446*(tmp_442*tmp_445 - tmp_442*tmp_58 - tmp_443*tmp_449) + tmp_66*(-tmp_449*tmp_55 - tmp_54*tmp_58 + tmp_54*tmp_64) + tmp_85*(-tmp_449*tmp_82 - tmp_58*tmp_81 + tmp_81*tmp_84);
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


      real_t tmp_0 = -p_affine_4_1;
      real_t tmp_1 = p_affine_5_1 + tmp_0;
      real_t tmp_2 = -p_affine_4_0;
      real_t tmp_3 = p_affine_6_0 + tmp_2;
      real_t tmp_4 = p_affine_7_1 + tmp_0;
      real_t tmp_5 = tmp_3*tmp_4;
      real_t tmp_6 = p_affine_7_0 + tmp_2;
      real_t tmp_7 = p_affine_6_1 + tmp_0;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = tmp_5 - tmp_8;
      real_t tmp_10 = p_affine_5_0 + tmp_2;
      real_t tmp_11 = -p_affine_4_2;
      real_t tmp_12 = p_affine_7_2 + tmp_11;
      real_t tmp_13 = tmp_12*tmp_7;
      real_t tmp_14 = p_affine_5_2 + tmp_11;
      real_t tmp_15 = p_affine_6_2 + tmp_11;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_15*tmp_4;
      real_t tmp_18 = tmp_12*tmp_3;
      real_t tmp_19 = 1.0 / (tmp_1*tmp_16 - tmp_1*tmp_18 + tmp_10*tmp_13 - tmp_10*tmp_17 + tmp_14*tmp_5 - tmp_14*tmp_8);
      real_t tmp_20 = p_affine_8_2 + tmp_11;
      real_t tmp_21 = -p_affine_8_2;
      real_t tmp_22 = p_affine_9_2 + tmp_21;
      real_t tmp_23 = p_affine_10_2 + tmp_21;
      real_t tmp_24 = 0.031405749086161582*tmp_22 + 0.93718850182767688*tmp_23;
      real_t tmp_25 = tmp_19*(tmp_20 + tmp_24);
      real_t tmp_26 = tmp_16 - tmp_18;
      real_t tmp_27 = p_affine_8_1 + tmp_0;
      real_t tmp_28 = -p_affine_8_1;
      real_t tmp_29 = p_affine_9_1 + tmp_28;
      real_t tmp_30 = p_affine_10_1 + tmp_28;
      real_t tmp_31 = 0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_30;
      real_t tmp_32 = tmp_19*(tmp_27 + tmp_31);
      real_t tmp_33 = tmp_13 - tmp_17;
      real_t tmp_34 = p_affine_8_0 + tmp_2;
      real_t tmp_35 = -p_affine_8_0;
      real_t tmp_36 = p_affine_9_0 + tmp_35;
      real_t tmp_37 = p_affine_10_0 + tmp_35;
      real_t tmp_38 = 0.031405749086161582*tmp_36 + 0.93718850182767688*tmp_37;
      real_t tmp_39 = tmp_19*(tmp_34 + tmp_38);
      real_t tmp_40 = tmp_1*tmp_6 - tmp_10*tmp_4;
      real_t tmp_41 = tmp_10*tmp_12 - tmp_14*tmp_6;
      real_t tmp_42 = -tmp_1*tmp_12 + tmp_14*tmp_4;
      real_t tmp_43 = -tmp_1*tmp_3 + tmp_10*tmp_7;
      real_t tmp_44 = -tmp_10*tmp_15 + tmp_14*tmp_3;
      real_t tmp_45 = tmp_1*tmp_15 - tmp_14*tmp_7;
      real_t tmp_46 = tmp_1*(tmp_25*tmp_9 + tmp_26*tmp_32 + tmp_33*tmp_39 - 1.0/4.0) + tmp_4*(tmp_25*tmp_43 + tmp_32*tmp_44 + tmp_39*tmp_45 - 1.0/4.0) + tmp_7*(tmp_25*tmp_40 + tmp_32*tmp_41 + tmp_39*tmp_42 - 1.0/4.0);
      real_t tmp_47 = -p_affine_0_1;
      real_t tmp_48 = p_affine_1_1 + tmp_47;
      real_t tmp_49 = -p_affine_0_2;
      real_t tmp_50 = p_affine_2_2 + tmp_49;
      real_t tmp_51 = tmp_48*tmp_50;
      real_t tmp_52 = p_affine_2_1 + tmp_47;
      real_t tmp_53 = p_affine_1_2 + tmp_49;
      real_t tmp_54 = tmp_52*tmp_53;
      real_t tmp_55 = -p_affine_0_0;
      real_t tmp_56 = p_affine_1_0 + tmp_55;
      real_t tmp_57 = p_affine_3_2 + tmp_49;
      real_t tmp_58 = tmp_52*tmp_57;
      real_t tmp_59 = p_affine_2_0 + tmp_55;
      real_t tmp_60 = p_affine_3_1 + tmp_47;
      real_t tmp_61 = tmp_53*tmp_60;
      real_t tmp_62 = p_affine_3_0 + tmp_55;
      real_t tmp_63 = tmp_50*tmp_60;
      real_t tmp_64 = tmp_48*tmp_57;
      real_t tmp_65 = 1.0 / (tmp_51*tmp_62 - tmp_54*tmp_62 + tmp_56*tmp_58 - tmp_56*tmp_63 + tmp_59*tmp_61 - tmp_59*tmp_64);
      real_t tmp_66 = tmp_65*(tmp_51 - tmp_54);
      real_t tmp_67 = tmp_65*(tmp_61 - tmp_64);
      real_t tmp_68 = tmp_65*(tmp_58 - tmp_63);
      real_t tmp_69 = tmp_65*(-tmp_50*tmp_56 + tmp_53*tmp_59);
      real_t tmp_70 = tmp_65*(-tmp_53*tmp_62 + tmp_56*tmp_57);
      real_t tmp_71 = tmp_65*(tmp_50*tmp_62 - tmp_57*tmp_59);
      real_t tmp_72 = tmp_65*(-tmp_48*tmp_59 + tmp_52*tmp_56);
      real_t tmp_73 = tmp_65*(tmp_48*tmp_62 - tmp_56*tmp_60);
      real_t tmp_74 = tmp_65*(-tmp_52*tmp_62 + tmp_59*tmp_60);
      real_t tmp_75 = 0.5*p_affine_13_0*(-tmp_66 - tmp_67 - tmp_68) + 0.5*p_affine_13_1*(-tmp_69 - tmp_70 - tmp_71) + 0.5*p_affine_13_2*(-tmp_72 - tmp_73 - tmp_74);
      real_t tmp_76 = p_affine_8_2 + tmp_49;
      real_t tmp_77 = tmp_24 + tmp_76;
      real_t tmp_78 = tmp_72*tmp_77;
      real_t tmp_79 = tmp_73*tmp_77;
      real_t tmp_80 = p_affine_8_1 + tmp_47;
      real_t tmp_81 = tmp_31 + tmp_80;
      real_t tmp_82 = tmp_69*tmp_81;
      real_t tmp_83 = tmp_70*tmp_81;
      real_t tmp_84 = tmp_74*tmp_77;
      real_t tmp_85 = tmp_71*tmp_81;
      real_t tmp_86 = p_affine_8_0 + tmp_55;
      real_t tmp_87 = tmp_38 + tmp_86;
      real_t tmp_88 = tmp_66*tmp_87;
      real_t tmp_89 = tmp_67*tmp_87;
      real_t tmp_90 = tmp_68*tmp_87;
      real_t tmp_91 = -tmp_78 - tmp_79 - tmp_82 - tmp_83 - tmp_84 - tmp_85 - tmp_88 - tmp_89 - tmp_90 + 1;
      real_t tmp_92 = tmp_1*tmp_19;
      real_t tmp_93 = tmp_19*tmp_7;
      real_t tmp_94 = tmp_19*tmp_4;
      real_t tmp_95 = 0.5*p_affine_13_0*(tmp_33*tmp_92 + tmp_42*tmp_93 + tmp_45*tmp_94) + 0.5*p_affine_13_1*(tmp_26*tmp_92 + tmp_41*tmp_93 + tmp_44*tmp_94) + 0.5*p_affine_13_2*(tmp_40*tmp_93 + tmp_43*tmp_94 + tmp_9*tmp_92);
      real_t tmp_96 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_97 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_98 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_99 = (std::abs(tmp_23*tmp_96 - tmp_30*tmp_98)*std::abs(tmp_23*tmp_96 - tmp_30*tmp_98)) + (std::abs(tmp_23*tmp_97 - tmp_37*tmp_98)*std::abs(tmp_23*tmp_97 - tmp_37*tmp_98)) + (std::abs(tmp_30*tmp_97 - tmp_37*tmp_96)*std::abs(tmp_30*tmp_97 - tmp_37*tmp_96));
      real_t tmp_100 = 5.0*std::pow(tmp_99, -0.25);
      real_t tmp_101 = tmp_100*tmp_46;
      real_t tmp_102 = 1.0*std::pow(tmp_99, 1.0/2.0);
      real_t tmp_103 = 0.0068572537431980923*tmp_102;
      real_t tmp_104 = 0.19601935860219369*tmp_22 + 0.60796128279561268*tmp_23;
      real_t tmp_105 = tmp_19*(tmp_104 + tmp_20);
      real_t tmp_106 = 0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_30;
      real_t tmp_107 = tmp_19*(tmp_106 + tmp_27);
      real_t tmp_108 = 0.19601935860219369*tmp_36 + 0.60796128279561268*tmp_37;
      real_t tmp_109 = tmp_19*(tmp_108 + tmp_34);
      real_t tmp_110 = tmp_1*(tmp_105*tmp_9 + tmp_107*tmp_26 + tmp_109*tmp_33 - 1.0/4.0) + tmp_4*(tmp_105*tmp_43 + tmp_107*tmp_44 + tmp_109*tmp_45 - 1.0/4.0) + tmp_7*(tmp_105*tmp_40 + tmp_107*tmp_41 + tmp_109*tmp_42 - 1.0/4.0);
      real_t tmp_111 = tmp_104 + tmp_76;
      real_t tmp_112 = tmp_111*tmp_72;
      real_t tmp_113 = tmp_111*tmp_73;
      real_t tmp_114 = tmp_106 + tmp_80;
      real_t tmp_115 = tmp_114*tmp_69;
      real_t tmp_116 = tmp_114*tmp_70;
      real_t tmp_117 = tmp_111*tmp_74;
      real_t tmp_118 = tmp_114*tmp_71;
      real_t tmp_119 = tmp_108 + tmp_86;
      real_t tmp_120 = tmp_119*tmp_66;
      real_t tmp_121 = tmp_119*tmp_67;
      real_t tmp_122 = tmp_119*tmp_68;
      real_t tmp_123 = -tmp_112 - tmp_113 - tmp_115 - tmp_116 - tmp_117 - tmp_118 - tmp_120 - tmp_121 - tmp_122 + 1;
      real_t tmp_124 = tmp_100*tmp_110;
      real_t tmp_125 = 0.037198804536718075*tmp_102;
      real_t tmp_126 = 0.37605877282253791*tmp_22 + 0.039308471900058539*tmp_23;
      real_t tmp_127 = tmp_19*(tmp_126 + tmp_20);
      real_t tmp_128 = 0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_30;
      real_t tmp_129 = tmp_19*(tmp_128 + tmp_27);
      real_t tmp_130 = 0.37605877282253791*tmp_36 + 0.039308471900058539*tmp_37;
      real_t tmp_131 = tmp_19*(tmp_130 + tmp_34);
      real_t tmp_132 = tmp_1*(tmp_127*tmp_9 + tmp_129*tmp_26 + tmp_131*tmp_33 - 1.0/4.0) + tmp_4*(tmp_127*tmp_43 + tmp_129*tmp_44 + tmp_131*tmp_45 - 1.0/4.0) + tmp_7*(tmp_127*tmp_40 + tmp_129*tmp_41 + tmp_131*tmp_42 - 1.0/4.0);
      real_t tmp_133 = tmp_126 + tmp_76;
      real_t tmp_134 = tmp_133*tmp_72;
      real_t tmp_135 = tmp_133*tmp_73;
      real_t tmp_136 = tmp_128 + tmp_80;
      real_t tmp_137 = tmp_136*tmp_69;
      real_t tmp_138 = tmp_136*tmp_70;
      real_t tmp_139 = tmp_133*tmp_74;
      real_t tmp_140 = tmp_136*tmp_71;
      real_t tmp_141 = tmp_130 + tmp_86;
      real_t tmp_142 = tmp_141*tmp_66;
      real_t tmp_143 = tmp_141*tmp_67;
      real_t tmp_144 = tmp_141*tmp_68;
      real_t tmp_145 = -tmp_134 - tmp_135 - tmp_137 - tmp_138 - tmp_139 - tmp_140 - tmp_142 - tmp_143 - tmp_144 + 1;
      real_t tmp_146 = tmp_100*tmp_132;
      real_t tmp_147 = 0.020848748529055869*tmp_102;
      real_t tmp_148 = 0.78764240869137092*tmp_22 + 0.1711304259088916*tmp_23;
      real_t tmp_149 = tmp_19*(tmp_148 + tmp_20);
      real_t tmp_150 = 0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_30;
      real_t tmp_151 = tmp_19*(tmp_150 + tmp_27);
      real_t tmp_152 = 0.78764240869137092*tmp_36 + 0.1711304259088916*tmp_37;
      real_t tmp_153 = tmp_19*(tmp_152 + tmp_34);
      real_t tmp_154 = tmp_1*(tmp_149*tmp_9 + tmp_151*tmp_26 + tmp_153*tmp_33 - 1.0/4.0) + tmp_4*(tmp_149*tmp_43 + tmp_151*tmp_44 + tmp_153*tmp_45 - 1.0/4.0) + tmp_7*(tmp_149*tmp_40 + tmp_151*tmp_41 + tmp_153*tmp_42 - 1.0/4.0);
      real_t tmp_155 = tmp_148 + tmp_76;
      real_t tmp_156 = tmp_155*tmp_72;
      real_t tmp_157 = tmp_155*tmp_73;
      real_t tmp_158 = tmp_150 + tmp_80;
      real_t tmp_159 = tmp_158*tmp_69;
      real_t tmp_160 = tmp_158*tmp_70;
      real_t tmp_161 = tmp_155*tmp_74;
      real_t tmp_162 = tmp_158*tmp_71;
      real_t tmp_163 = tmp_152 + tmp_86;
      real_t tmp_164 = tmp_163*tmp_66;
      real_t tmp_165 = tmp_163*tmp_67;
      real_t tmp_166 = tmp_163*tmp_68;
      real_t tmp_167 = -tmp_156 - tmp_157 - tmp_159 - tmp_160 - tmp_161 - tmp_162 - tmp_164 - tmp_165 - tmp_166 + 1;
      real_t tmp_168 = tmp_100*tmp_154;
      real_t tmp_169 = 0.019202922745021479*tmp_102;
      real_t tmp_170 = 0.58463275527740355*tmp_22 + 0.37605877282253791*tmp_23;
      real_t tmp_171 = tmp_19*(tmp_170 + tmp_20);
      real_t tmp_172 = 0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_30;
      real_t tmp_173 = tmp_19*(tmp_172 + tmp_27);
      real_t tmp_174 = 0.58463275527740355*tmp_36 + 0.37605877282253791*tmp_37;
      real_t tmp_175 = tmp_19*(tmp_174 + tmp_34);
      real_t tmp_176 = tmp_1*(tmp_171*tmp_9 + tmp_173*tmp_26 + tmp_175*tmp_33 - 1.0/4.0) + tmp_4*(tmp_171*tmp_43 + tmp_173*tmp_44 + tmp_175*tmp_45 - 1.0/4.0) + tmp_7*(tmp_171*tmp_40 + tmp_173*tmp_41 + tmp_175*tmp_42 - 1.0/4.0);
      real_t tmp_177 = tmp_170 + tmp_76;
      real_t tmp_178 = tmp_177*tmp_72;
      real_t tmp_179 = tmp_177*tmp_73;
      real_t tmp_180 = tmp_172 + tmp_80;
      real_t tmp_181 = tmp_180*tmp_69;
      real_t tmp_182 = tmp_180*tmp_70;
      real_t tmp_183 = tmp_177*tmp_74;
      real_t tmp_184 = tmp_180*tmp_71;
      real_t tmp_185 = tmp_174 + tmp_86;
      real_t tmp_186 = tmp_185*tmp_66;
      real_t tmp_187 = tmp_185*tmp_67;
      real_t tmp_188 = tmp_185*tmp_68;
      real_t tmp_189 = -tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 - tmp_184 - tmp_186 - tmp_187 - tmp_188 + 1;
      real_t tmp_190 = tmp_100*tmp_176;
      real_t tmp_191 = 0.020848748529055869*tmp_102;
      real_t tmp_192 = 0.041227165399737475*tmp_22 + 0.78764240869137092*tmp_23;
      real_t tmp_193 = tmp_19*(tmp_192 + tmp_20);
      real_t tmp_194 = 0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_30;
      real_t tmp_195 = tmp_19*(tmp_194 + tmp_27);
      real_t tmp_196 = 0.041227165399737475*tmp_36 + 0.78764240869137092*tmp_37;
      real_t tmp_197 = tmp_19*(tmp_196 + tmp_34);
      real_t tmp_198 = tmp_1*(tmp_193*tmp_9 + tmp_195*tmp_26 + tmp_197*tmp_33 - 1.0/4.0) + tmp_4*(tmp_193*tmp_43 + tmp_195*tmp_44 + tmp_197*tmp_45 - 1.0/4.0) + tmp_7*(tmp_193*tmp_40 + tmp_195*tmp_41 + tmp_197*tmp_42 - 1.0/4.0);
      real_t tmp_199 = tmp_192 + tmp_76;
      real_t tmp_200 = tmp_199*tmp_72;
      real_t tmp_201 = tmp_199*tmp_73;
      real_t tmp_202 = tmp_194 + tmp_80;
      real_t tmp_203 = tmp_202*tmp_69;
      real_t tmp_204 = tmp_202*tmp_70;
      real_t tmp_205 = tmp_199*tmp_74;
      real_t tmp_206 = tmp_202*tmp_71;
      real_t tmp_207 = tmp_196 + tmp_86;
      real_t tmp_208 = tmp_207*tmp_66;
      real_t tmp_209 = tmp_207*tmp_67;
      real_t tmp_210 = tmp_207*tmp_68;
      real_t tmp_211 = -tmp_200 - tmp_201 - tmp_203 - tmp_204 - tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 + 1;
      real_t tmp_212 = tmp_100*tmp_198;
      real_t tmp_213 = 0.019202922745021479*tmp_102;
      real_t tmp_214 = 0.039308471900058539*tmp_22 + 0.58463275527740355*tmp_23;
      real_t tmp_215 = tmp_19*(tmp_20 + tmp_214);
      real_t tmp_216 = 0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_30;
      real_t tmp_217 = tmp_19*(tmp_216 + tmp_27);
      real_t tmp_218 = 0.039308471900058539*tmp_36 + 0.58463275527740355*tmp_37;
      real_t tmp_219 = tmp_19*(tmp_218 + tmp_34);
      real_t tmp_220 = tmp_1*(tmp_215*tmp_9 + tmp_217*tmp_26 + tmp_219*tmp_33 - 1.0/4.0) + tmp_4*(tmp_215*tmp_43 + tmp_217*tmp_44 + tmp_219*tmp_45 - 1.0/4.0) + tmp_7*(tmp_215*tmp_40 + tmp_217*tmp_41 + tmp_219*tmp_42 - 1.0/4.0);
      real_t tmp_221 = tmp_214 + tmp_76;
      real_t tmp_222 = tmp_221*tmp_72;
      real_t tmp_223 = tmp_221*tmp_73;
      real_t tmp_224 = tmp_216 + tmp_80;
      real_t tmp_225 = tmp_224*tmp_69;
      real_t tmp_226 = tmp_224*tmp_70;
      real_t tmp_227 = tmp_221*tmp_74;
      real_t tmp_228 = tmp_224*tmp_71;
      real_t tmp_229 = tmp_218 + tmp_86;
      real_t tmp_230 = tmp_229*tmp_66;
      real_t tmp_231 = tmp_229*tmp_67;
      real_t tmp_232 = tmp_229*tmp_68;
      real_t tmp_233 = -tmp_222 - tmp_223 - tmp_225 - tmp_226 - tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 + 1;
      real_t tmp_234 = tmp_100*tmp_220;
      real_t tmp_235 = 0.020848748529055869*tmp_102;
      real_t tmp_236 = 0.78764240869137092*tmp_22 + 0.041227165399737475*tmp_23;
      real_t tmp_237 = tmp_19*(tmp_20 + tmp_236);
      real_t tmp_238 = 0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_30;
      real_t tmp_239 = tmp_19*(tmp_238 + tmp_27);
      real_t tmp_240 = 0.78764240869137092*tmp_36 + 0.041227165399737475*tmp_37;
      real_t tmp_241 = tmp_19*(tmp_240 + tmp_34);
      real_t tmp_242 = tmp_1*(tmp_237*tmp_9 + tmp_239*tmp_26 + tmp_241*tmp_33 - 1.0/4.0) + tmp_4*(tmp_237*tmp_43 + tmp_239*tmp_44 + tmp_241*tmp_45 - 1.0/4.0) + tmp_7*(tmp_237*tmp_40 + tmp_239*tmp_41 + tmp_241*tmp_42 - 1.0/4.0);
      real_t tmp_243 = tmp_236 + tmp_76;
      real_t tmp_244 = tmp_243*tmp_72;
      real_t tmp_245 = tmp_243*tmp_73;
      real_t tmp_246 = tmp_238 + tmp_80;
      real_t tmp_247 = tmp_246*tmp_69;
      real_t tmp_248 = tmp_246*tmp_70;
      real_t tmp_249 = tmp_243*tmp_74;
      real_t tmp_250 = tmp_246*tmp_71;
      real_t tmp_251 = tmp_240 + tmp_86;
      real_t tmp_252 = tmp_251*tmp_66;
      real_t tmp_253 = tmp_251*tmp_67;
      real_t tmp_254 = tmp_251*tmp_68;
      real_t tmp_255 = -tmp_244 - tmp_245 - tmp_247 - tmp_248 - tmp_249 - tmp_250 - tmp_252 - tmp_253 - tmp_254 + 1;
      real_t tmp_256 = tmp_100*tmp_242;
      real_t tmp_257 = 0.019202922745021479*tmp_102;
      real_t tmp_258 = 0.58463275527740355*tmp_22 + 0.039308471900058539*tmp_23;
      real_t tmp_259 = tmp_19*(tmp_20 + tmp_258);
      real_t tmp_260 = 0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_30;
      real_t tmp_261 = tmp_19*(tmp_260 + tmp_27);
      real_t tmp_262 = 0.58463275527740355*tmp_36 + 0.039308471900058539*tmp_37;
      real_t tmp_263 = tmp_19*(tmp_262 + tmp_34);
      real_t tmp_264 = tmp_1*(tmp_259*tmp_9 + tmp_26*tmp_261 + tmp_263*tmp_33 - 1.0/4.0) + tmp_4*(tmp_259*tmp_43 + tmp_261*tmp_44 + tmp_263*tmp_45 - 1.0/4.0) + tmp_7*(tmp_259*tmp_40 + tmp_261*tmp_41 + tmp_263*tmp_42 - 1.0/4.0);
      real_t tmp_265 = tmp_258 + tmp_76;
      real_t tmp_266 = tmp_265*tmp_72;
      real_t tmp_267 = tmp_265*tmp_73;
      real_t tmp_268 = tmp_260 + tmp_80;
      real_t tmp_269 = tmp_268*tmp_69;
      real_t tmp_270 = tmp_268*tmp_70;
      real_t tmp_271 = tmp_265*tmp_74;
      real_t tmp_272 = tmp_268*tmp_71;
      real_t tmp_273 = tmp_262 + tmp_86;
      real_t tmp_274 = tmp_273*tmp_66;
      real_t tmp_275 = tmp_273*tmp_67;
      real_t tmp_276 = tmp_273*tmp_68;
      real_t tmp_277 = -tmp_266 - tmp_267 - tmp_269 - tmp_270 - tmp_271 - tmp_272 - tmp_274 - tmp_275 - tmp_276 + 1;
      real_t tmp_278 = tmp_100*tmp_264;
      real_t tmp_279 = 0.020848748529055869*tmp_102;
      real_t tmp_280 = 0.1711304259088916*tmp_22 + 0.78764240869137092*tmp_23;
      real_t tmp_281 = tmp_19*(tmp_20 + tmp_280);
      real_t tmp_282 = 0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_30;
      real_t tmp_283 = tmp_19*(tmp_27 + tmp_282);
      real_t tmp_284 = 0.1711304259088916*tmp_36 + 0.78764240869137092*tmp_37;
      real_t tmp_285 = tmp_19*(tmp_284 + tmp_34);
      real_t tmp_286 = tmp_1*(tmp_26*tmp_283 + tmp_281*tmp_9 + tmp_285*tmp_33 - 1.0/4.0) + tmp_4*(tmp_281*tmp_43 + tmp_283*tmp_44 + tmp_285*tmp_45 - 1.0/4.0) + tmp_7*(tmp_281*tmp_40 + tmp_283*tmp_41 + tmp_285*tmp_42 - 1.0/4.0);
      real_t tmp_287 = tmp_280 + tmp_76;
      real_t tmp_288 = tmp_287*tmp_72;
      real_t tmp_289 = tmp_287*tmp_73;
      real_t tmp_290 = tmp_282 + tmp_80;
      real_t tmp_291 = tmp_290*tmp_69;
      real_t tmp_292 = tmp_290*tmp_70;
      real_t tmp_293 = tmp_287*tmp_74;
      real_t tmp_294 = tmp_290*tmp_71;
      real_t tmp_295 = tmp_284 + tmp_86;
      real_t tmp_296 = tmp_295*tmp_66;
      real_t tmp_297 = tmp_295*tmp_67;
      real_t tmp_298 = tmp_295*tmp_68;
      real_t tmp_299 = -tmp_288 - tmp_289 - tmp_291 - tmp_292 - tmp_293 - tmp_294 - tmp_296 - tmp_297 - tmp_298 + 1;
      real_t tmp_300 = tmp_100*tmp_286;
      real_t tmp_301 = 0.019202922745021479*tmp_102;
      real_t tmp_302 = 0.37605877282253791*tmp_22 + 0.58463275527740355*tmp_23;
      real_t tmp_303 = tmp_19*(tmp_20 + tmp_302);
      real_t tmp_304 = 0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_30;
      real_t tmp_305 = tmp_19*(tmp_27 + tmp_304);
      real_t tmp_306 = 0.37605877282253791*tmp_36 + 0.58463275527740355*tmp_37;
      real_t tmp_307 = tmp_19*(tmp_306 + tmp_34);
      real_t tmp_308 = tmp_1*(tmp_26*tmp_305 + tmp_303*tmp_9 + tmp_307*tmp_33 - 1.0/4.0) + tmp_4*(tmp_303*tmp_43 + tmp_305*tmp_44 + tmp_307*tmp_45 - 1.0/4.0) + tmp_7*(tmp_303*tmp_40 + tmp_305*tmp_41 + tmp_307*tmp_42 - 1.0/4.0);
      real_t tmp_309 = tmp_302 + tmp_76;
      real_t tmp_310 = tmp_309*tmp_72;
      real_t tmp_311 = tmp_309*tmp_73;
      real_t tmp_312 = tmp_304 + tmp_80;
      real_t tmp_313 = tmp_312*tmp_69;
      real_t tmp_314 = tmp_312*tmp_70;
      real_t tmp_315 = tmp_309*tmp_74;
      real_t tmp_316 = tmp_312*tmp_71;
      real_t tmp_317 = tmp_306 + tmp_86;
      real_t tmp_318 = tmp_317*tmp_66;
      real_t tmp_319 = tmp_317*tmp_67;
      real_t tmp_320 = tmp_317*tmp_68;
      real_t tmp_321 = -tmp_310 - tmp_311 - tmp_313 - tmp_314 - tmp_315 - tmp_316 - tmp_318 - tmp_319 - tmp_320 + 1;
      real_t tmp_322 = tmp_100*tmp_308;
      real_t tmp_323 = 0.020848748529055869*tmp_102;
      real_t tmp_324 = 0.041227165399737475*tmp_22 + 0.1711304259088916*tmp_23;
      real_t tmp_325 = tmp_19*(tmp_20 + tmp_324);
      real_t tmp_326 = 0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_30;
      real_t tmp_327 = tmp_19*(tmp_27 + tmp_326);
      real_t tmp_328 = 0.041227165399737475*tmp_36 + 0.1711304259088916*tmp_37;
      real_t tmp_329 = tmp_19*(tmp_328 + tmp_34);
      real_t tmp_330 = tmp_1*(tmp_26*tmp_327 + tmp_325*tmp_9 + tmp_329*tmp_33 - 1.0/4.0) + tmp_4*(tmp_325*tmp_43 + tmp_327*tmp_44 + tmp_329*tmp_45 - 1.0/4.0) + tmp_7*(tmp_325*tmp_40 + tmp_327*tmp_41 + tmp_329*tmp_42 - 1.0/4.0);
      real_t tmp_331 = tmp_324 + tmp_76;
      real_t tmp_332 = tmp_331*tmp_72;
      real_t tmp_333 = tmp_331*tmp_73;
      real_t tmp_334 = tmp_326 + tmp_80;
      real_t tmp_335 = tmp_334*tmp_69;
      real_t tmp_336 = tmp_334*tmp_70;
      real_t tmp_337 = tmp_331*tmp_74;
      real_t tmp_338 = tmp_334*tmp_71;
      real_t tmp_339 = tmp_328 + tmp_86;
      real_t tmp_340 = tmp_339*tmp_66;
      real_t tmp_341 = tmp_339*tmp_67;
      real_t tmp_342 = tmp_339*tmp_68;
      real_t tmp_343 = -tmp_332 - tmp_333 - tmp_335 - tmp_336 - tmp_337 - tmp_338 - tmp_340 - tmp_341 - tmp_342 + 1;
      real_t tmp_344 = tmp_100*tmp_330;
      real_t tmp_345 = 0.019202922745021479*tmp_102;
      real_t tmp_346 = 0.40446199974765351*tmp_22 + 0.19107600050469298*tmp_23;
      real_t tmp_347 = tmp_19*(tmp_20 + tmp_346);
      real_t tmp_348 = 0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_30;
      real_t tmp_349 = tmp_19*(tmp_27 + tmp_348);
      real_t tmp_350 = 0.40446199974765351*tmp_36 + 0.19107600050469298*tmp_37;
      real_t tmp_351 = tmp_19*(tmp_34 + tmp_350);
      real_t tmp_352 = tmp_1*(tmp_26*tmp_349 + tmp_33*tmp_351 + tmp_347*tmp_9 - 1.0/4.0) + tmp_4*(tmp_347*tmp_43 + tmp_349*tmp_44 + tmp_351*tmp_45 - 1.0/4.0) + tmp_7*(tmp_347*tmp_40 + tmp_349*tmp_41 + tmp_351*tmp_42 - 1.0/4.0);
      real_t tmp_353 = tmp_346 + tmp_76;
      real_t tmp_354 = tmp_353*tmp_72;
      real_t tmp_355 = tmp_353*tmp_73;
      real_t tmp_356 = tmp_348 + tmp_80;
      real_t tmp_357 = tmp_356*tmp_69;
      real_t tmp_358 = tmp_356*tmp_70;
      real_t tmp_359 = tmp_353*tmp_74;
      real_t tmp_360 = tmp_356*tmp_71;
      real_t tmp_361 = tmp_350 + tmp_86;
      real_t tmp_362 = tmp_361*tmp_66;
      real_t tmp_363 = tmp_361*tmp_67;
      real_t tmp_364 = tmp_361*tmp_68;
      real_t tmp_365 = -tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 - tmp_360 - tmp_362 - tmp_363 - tmp_364 + 1;
      real_t tmp_366 = tmp_100*tmp_352;
      real_t tmp_367 = 0.042507265838595799*tmp_102;
      real_t tmp_368 = 0.039308471900058539*tmp_22 + 0.37605877282253791*tmp_23;
      real_t tmp_369 = tmp_19*(tmp_20 + tmp_368);
      real_t tmp_370 = 0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_30;
      real_t tmp_371 = tmp_19*(tmp_27 + tmp_370);
      real_t tmp_372 = 0.039308471900058539*tmp_36 + 0.37605877282253791*tmp_37;
      real_t tmp_373 = tmp_19*(tmp_34 + tmp_372);
      real_t tmp_374 = tmp_1*(tmp_26*tmp_371 + tmp_33*tmp_373 + tmp_369*tmp_9 - 1.0/4.0) + tmp_4*(tmp_369*tmp_43 + tmp_371*tmp_44 + tmp_373*tmp_45 - 1.0/4.0) + tmp_7*(tmp_369*tmp_40 + tmp_371*tmp_41 + tmp_373*tmp_42 - 1.0/4.0);
      real_t tmp_375 = tmp_368 + tmp_76;
      real_t tmp_376 = tmp_375*tmp_72;
      real_t tmp_377 = tmp_375*tmp_73;
      real_t tmp_378 = tmp_370 + tmp_80;
      real_t tmp_379 = tmp_378*tmp_69;
      real_t tmp_380 = tmp_378*tmp_70;
      real_t tmp_381 = tmp_375*tmp_74;
      real_t tmp_382 = tmp_378*tmp_71;
      real_t tmp_383 = tmp_372 + tmp_86;
      real_t tmp_384 = tmp_383*tmp_66;
      real_t tmp_385 = tmp_383*tmp_67;
      real_t tmp_386 = tmp_383*tmp_68;
      real_t tmp_387 = -tmp_376 - tmp_377 - tmp_379 - tmp_380 - tmp_381 - tmp_382 - tmp_384 - tmp_385 - tmp_386 + 1;
      real_t tmp_388 = tmp_100*tmp_374;
      real_t tmp_389 = 0.020848748529055869*tmp_102;
      real_t tmp_390 = 0.93718850182767688*tmp_22 + 0.031405749086161582*tmp_23;
      real_t tmp_391 = tmp_19*(tmp_20 + tmp_390);
      real_t tmp_392 = 0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_30;
      real_t tmp_393 = tmp_19*(tmp_27 + tmp_392);
      real_t tmp_394 = 0.93718850182767688*tmp_36 + 0.031405749086161582*tmp_37;
      real_t tmp_395 = tmp_19*(tmp_34 + tmp_394);
      real_t tmp_396 = tmp_1*(tmp_26*tmp_393 + tmp_33*tmp_395 + tmp_391*tmp_9 - 1.0/4.0) + tmp_4*(tmp_391*tmp_43 + tmp_393*tmp_44 + tmp_395*tmp_45 - 1.0/4.0) + tmp_7*(tmp_391*tmp_40 + tmp_393*tmp_41 + tmp_395*tmp_42 - 1.0/4.0);
      real_t tmp_397 = tmp_390 + tmp_76;
      real_t tmp_398 = tmp_397*tmp_72;
      real_t tmp_399 = tmp_397*tmp_73;
      real_t tmp_400 = tmp_392 + tmp_80;
      real_t tmp_401 = tmp_400*tmp_69;
      real_t tmp_402 = tmp_400*tmp_70;
      real_t tmp_403 = tmp_397*tmp_74;
      real_t tmp_404 = tmp_400*tmp_71;
      real_t tmp_405 = tmp_394 + tmp_86;
      real_t tmp_406 = tmp_405*tmp_66;
      real_t tmp_407 = tmp_405*tmp_67;
      real_t tmp_408 = tmp_405*tmp_68;
      real_t tmp_409 = -tmp_398 - tmp_399 - tmp_401 - tmp_402 - tmp_403 - tmp_404 - tmp_406 - tmp_407 - tmp_408 + 1;
      real_t tmp_410 = tmp_100*tmp_396;
      real_t tmp_411 = 0.0068572537431980923*tmp_102;
      real_t tmp_412 = 0.60796128279561268*tmp_22 + 0.19601935860219369*tmp_23;
      real_t tmp_413 = tmp_19*(tmp_20 + tmp_412);
      real_t tmp_414 = 0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_30;
      real_t tmp_415 = tmp_19*(tmp_27 + tmp_414);
      real_t tmp_416 = 0.60796128279561268*tmp_36 + 0.19601935860219369*tmp_37;
      real_t tmp_417 = tmp_19*(tmp_34 + tmp_416);
      real_t tmp_418 = tmp_1*(tmp_26*tmp_415 + tmp_33*tmp_417 + tmp_413*tmp_9 - 1.0/4.0) + tmp_4*(tmp_413*tmp_43 + tmp_415*tmp_44 + tmp_417*tmp_45 - 1.0/4.0) + tmp_7*(tmp_40*tmp_413 + tmp_41*tmp_415 + tmp_417*tmp_42 - 1.0/4.0);
      real_t tmp_419 = tmp_412 + tmp_76;
      real_t tmp_420 = tmp_419*tmp_72;
      real_t tmp_421 = tmp_419*tmp_73;
      real_t tmp_422 = tmp_414 + tmp_80;
      real_t tmp_423 = tmp_422*tmp_69;
      real_t tmp_424 = tmp_422*tmp_70;
      real_t tmp_425 = tmp_419*tmp_74;
      real_t tmp_426 = tmp_422*tmp_71;
      real_t tmp_427 = tmp_416 + tmp_86;
      real_t tmp_428 = tmp_427*tmp_66;
      real_t tmp_429 = tmp_427*tmp_67;
      real_t tmp_430 = tmp_427*tmp_68;
      real_t tmp_431 = -tmp_420 - tmp_421 - tmp_423 - tmp_424 - tmp_425 - tmp_426 - tmp_428 - tmp_429 - tmp_430 + 1;
      real_t tmp_432 = tmp_100*tmp_418;
      real_t tmp_433 = 0.037198804536718075*tmp_102;
      real_t tmp_434 = 0.19107600050469298*tmp_22 + 0.40446199974765351*tmp_23;
      real_t tmp_435 = tmp_19*(tmp_20 + tmp_434);
      real_t tmp_436 = 0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_30;
      real_t tmp_437 = tmp_19*(tmp_27 + tmp_436);
      real_t tmp_438 = 0.19107600050469298*tmp_36 + 0.40446199974765351*tmp_37;
      real_t tmp_439 = tmp_19*(tmp_34 + tmp_438);
      real_t tmp_440 = tmp_1*(tmp_26*tmp_437 + tmp_33*tmp_439 + tmp_435*tmp_9 - 1.0/4.0) + tmp_4*(tmp_43*tmp_435 + tmp_437*tmp_44 + tmp_439*tmp_45 - 1.0/4.0) + tmp_7*(tmp_40*tmp_435 + tmp_41*tmp_437 + tmp_42*tmp_439 - 1.0/4.0);
      real_t tmp_441 = tmp_434 + tmp_76;
      real_t tmp_442 = tmp_441*tmp_72;
      real_t tmp_443 = tmp_441*tmp_73;
      real_t tmp_444 = tmp_436 + tmp_80;
      real_t tmp_445 = tmp_444*tmp_69;
      real_t tmp_446 = tmp_444*tmp_70;
      real_t tmp_447 = tmp_441*tmp_74;
      real_t tmp_448 = tmp_444*tmp_71;
      real_t tmp_449 = tmp_438 + tmp_86;
      real_t tmp_450 = tmp_449*tmp_66;
      real_t tmp_451 = tmp_449*tmp_67;
      real_t tmp_452 = tmp_449*tmp_68;
      real_t tmp_453 = -tmp_442 - tmp_443 - tmp_445 - tmp_446 - tmp_447 - tmp_448 - tmp_450 - tmp_451 - tmp_452 + 1;
      real_t tmp_454 = tmp_100*tmp_440;
      real_t tmp_455 = 0.042507265838595799*tmp_102;
      real_t tmp_456 = 0.031405749086161582*tmp_22 + 0.031405749086161582*tmp_23;
      real_t tmp_457 = tmp_19*(tmp_20 + tmp_456);
      real_t tmp_458 = 0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_30;
      real_t tmp_459 = tmp_19*(tmp_27 + tmp_458);
      real_t tmp_460 = 0.031405749086161582*tmp_36 + 0.031405749086161582*tmp_37;
      real_t tmp_461 = tmp_19*(tmp_34 + tmp_460);
      real_t tmp_462 = tmp_1*(tmp_26*tmp_459 + tmp_33*tmp_461 + tmp_457*tmp_9 - 1.0/4.0) + tmp_4*(tmp_43*tmp_457 + tmp_44*tmp_459 + tmp_45*tmp_461 - 1.0/4.0) + tmp_7*(tmp_40*tmp_457 + tmp_41*tmp_459 + tmp_42*tmp_461 - 1.0/4.0);
      real_t tmp_463 = tmp_456 + tmp_76;
      real_t tmp_464 = tmp_463*tmp_72;
      real_t tmp_465 = tmp_463*tmp_73;
      real_t tmp_466 = tmp_458 + tmp_80;
      real_t tmp_467 = tmp_466*tmp_69;
      real_t tmp_468 = tmp_466*tmp_70;
      real_t tmp_469 = tmp_463*tmp_74;
      real_t tmp_470 = tmp_466*tmp_71;
      real_t tmp_471 = tmp_460 + tmp_86;
      real_t tmp_472 = tmp_471*tmp_66;
      real_t tmp_473 = tmp_471*tmp_67;
      real_t tmp_474 = tmp_471*tmp_68;
      real_t tmp_475 = -tmp_464 - tmp_465 - tmp_467 - tmp_468 - tmp_469 - tmp_470 - tmp_472 - tmp_473 - tmp_474 + 1;
      real_t tmp_476 = tmp_100*tmp_462;
      real_t tmp_477 = 0.0068572537431980923*tmp_102;
      real_t tmp_478 = 0.19601935860219369*tmp_22 + 0.19601935860219369*tmp_23;
      real_t tmp_479 = tmp_19*(tmp_20 + tmp_478);
      real_t tmp_480 = 0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_30;
      real_t tmp_481 = tmp_19*(tmp_27 + tmp_480);
      real_t tmp_482 = 0.19601935860219369*tmp_36 + 0.19601935860219369*tmp_37;
      real_t tmp_483 = tmp_19*(tmp_34 + tmp_482);
      real_t tmp_484 = tmp_1*(tmp_26*tmp_481 + tmp_33*tmp_483 + tmp_479*tmp_9 - 1.0/4.0) + tmp_4*(tmp_43*tmp_479 + tmp_44*tmp_481 + tmp_45*tmp_483 - 1.0/4.0) + tmp_7*(tmp_40*tmp_479 + tmp_41*tmp_481 + tmp_42*tmp_483 - 1.0/4.0);
      real_t tmp_485 = tmp_478 + tmp_76;
      real_t tmp_486 = tmp_485*tmp_72;
      real_t tmp_487 = tmp_485*tmp_73;
      real_t tmp_488 = tmp_480 + tmp_80;
      real_t tmp_489 = tmp_488*tmp_69;
      real_t tmp_490 = tmp_488*tmp_70;
      real_t tmp_491 = tmp_485*tmp_74;
      real_t tmp_492 = tmp_488*tmp_71;
      real_t tmp_493 = tmp_482 + tmp_86;
      real_t tmp_494 = tmp_493*tmp_66;
      real_t tmp_495 = tmp_493*tmp_67;
      real_t tmp_496 = tmp_493*tmp_68;
      real_t tmp_497 = -tmp_486 - tmp_487 - tmp_489 - tmp_490 - tmp_491 - tmp_492 - tmp_494 - tmp_495 - tmp_496 + 1;
      real_t tmp_498 = tmp_100*tmp_484;
      real_t tmp_499 = 0.037198804536718075*tmp_102;
      real_t tmp_500 = 0.40446199974765351*tmp_22 + 0.40446199974765351*tmp_23;
      real_t tmp_501 = tmp_19*(tmp_20 + tmp_500);
      real_t tmp_502 = 0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_30;
      real_t tmp_503 = tmp_19*(tmp_27 + tmp_502);
      real_t tmp_504 = 0.40446199974765351*tmp_36 + 0.40446199974765351*tmp_37;
      real_t tmp_505 = tmp_19*(tmp_34 + tmp_504);
      real_t tmp_506 = tmp_1*(tmp_26*tmp_503 + tmp_33*tmp_505 + tmp_501*tmp_9 - 1.0/4.0) + tmp_4*(tmp_43*tmp_501 + tmp_44*tmp_503 + tmp_45*tmp_505 - 1.0/4.0) + tmp_7*(tmp_40*tmp_501 + tmp_41*tmp_503 + tmp_42*tmp_505 - 1.0/4.0);
      real_t tmp_507 = tmp_500 + tmp_76;
      real_t tmp_508 = tmp_507*tmp_72;
      real_t tmp_509 = tmp_507*tmp_73;
      real_t tmp_510 = tmp_502 + tmp_80;
      real_t tmp_511 = tmp_510*tmp_69;
      real_t tmp_512 = tmp_510*tmp_70;
      real_t tmp_513 = tmp_507*tmp_74;
      real_t tmp_514 = tmp_510*tmp_71;
      real_t tmp_515 = tmp_504 + tmp_86;
      real_t tmp_516 = tmp_515*tmp_66;
      real_t tmp_517 = tmp_515*tmp_67;
      real_t tmp_518 = tmp_515*tmp_68;
      real_t tmp_519 = -tmp_508 - tmp_509 - tmp_511 - tmp_512 - tmp_513 - tmp_514 - tmp_516 - tmp_517 - tmp_518 + 1;
      real_t tmp_520 = tmp_100*tmp_506;
      real_t tmp_521 = 0.042507265838595799*tmp_102;
      real_t tmp_522 = 0.1711304259088916*tmp_22 + 0.041227165399737475*tmp_23;
      real_t tmp_523 = tmp_19*(tmp_20 + tmp_522);
      real_t tmp_524 = 0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_30;
      real_t tmp_525 = tmp_19*(tmp_27 + tmp_524);
      real_t tmp_526 = 0.1711304259088916*tmp_36 + 0.041227165399737475*tmp_37;
      real_t tmp_527 = tmp_19*(tmp_34 + tmp_526);
      real_t tmp_528 = tmp_1*(tmp_26*tmp_525 + tmp_33*tmp_527 + tmp_523*tmp_9 - 1.0/4.0) + tmp_4*(tmp_43*tmp_523 + tmp_44*tmp_525 + tmp_45*tmp_527 - 1.0/4.0) + tmp_7*(tmp_40*tmp_523 + tmp_41*tmp_525 + tmp_42*tmp_527 - 1.0/4.0);
      real_t tmp_529 = tmp_522 + tmp_76;
      real_t tmp_530 = tmp_529*tmp_72;
      real_t tmp_531 = tmp_529*tmp_73;
      real_t tmp_532 = tmp_524 + tmp_80;
      real_t tmp_533 = tmp_532*tmp_69;
      real_t tmp_534 = tmp_532*tmp_70;
      real_t tmp_535 = tmp_529*tmp_74;
      real_t tmp_536 = tmp_532*tmp_71;
      real_t tmp_537 = tmp_526 + tmp_86;
      real_t tmp_538 = tmp_537*tmp_66;
      real_t tmp_539 = tmp_537*tmp_67;
      real_t tmp_540 = tmp_537*tmp_68;
      real_t tmp_541 = -tmp_530 - tmp_531 - tmp_533 - tmp_534 - tmp_535 - tmp_536 - tmp_538 - tmp_539 - tmp_540 + 1;
      real_t tmp_542 = tmp_100*tmp_528;
      real_t tmp_543 = 0.019202922745021479*tmp_102;
      real_t tmp_544 = tmp_84 + tmp_85 + tmp_90;
      real_t tmp_545 = 0.5*p_affine_13_0*tmp_68 + 0.5*p_affine_13_1*tmp_71 + 0.5*p_affine_13_2*tmp_74;
      real_t tmp_546 = tmp_117 + tmp_118 + tmp_122;
      real_t tmp_547 = tmp_139 + tmp_140 + tmp_144;
      real_t tmp_548 = tmp_161 + tmp_162 + tmp_166;
      real_t tmp_549 = tmp_183 + tmp_184 + tmp_188;
      real_t tmp_550 = tmp_205 + tmp_206 + tmp_210;
      real_t tmp_551 = tmp_227 + tmp_228 + tmp_232;
      real_t tmp_552 = tmp_249 + tmp_250 + tmp_254;
      real_t tmp_553 = tmp_271 + tmp_272 + tmp_276;
      real_t tmp_554 = tmp_293 + tmp_294 + tmp_298;
      real_t tmp_555 = tmp_315 + tmp_316 + tmp_320;
      real_t tmp_556 = tmp_337 + tmp_338 + tmp_342;
      real_t tmp_557 = tmp_359 + tmp_360 + tmp_364;
      real_t tmp_558 = tmp_381 + tmp_382 + tmp_386;
      real_t tmp_559 = tmp_403 + tmp_404 + tmp_408;
      real_t tmp_560 = tmp_425 + tmp_426 + tmp_430;
      real_t tmp_561 = tmp_447 + tmp_448 + tmp_452;
      real_t tmp_562 = tmp_469 + tmp_470 + tmp_474;
      real_t tmp_563 = tmp_491 + tmp_492 + tmp_496;
      real_t tmp_564 = tmp_513 + tmp_514 + tmp_518;
      real_t tmp_565 = tmp_535 + tmp_536 + tmp_540;
      real_t tmp_566 = tmp_79 + tmp_83 + tmp_89;
      real_t tmp_567 = 0.5*p_affine_13_0*tmp_67 + 0.5*p_affine_13_1*tmp_70 + 0.5*p_affine_13_2*tmp_73;
      real_t tmp_568 = tmp_113 + tmp_116 + tmp_121;
      real_t tmp_569 = tmp_135 + tmp_138 + tmp_143;
      real_t tmp_570 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_571 = tmp_179 + tmp_182 + tmp_187;
      real_t tmp_572 = tmp_201 + tmp_204 + tmp_209;
      real_t tmp_573 = tmp_223 + tmp_226 + tmp_231;
      real_t tmp_574 = tmp_245 + tmp_248 + tmp_253;
      real_t tmp_575 = tmp_267 + tmp_270 + tmp_275;
      real_t tmp_576 = tmp_289 + tmp_292 + tmp_297;
      real_t tmp_577 = tmp_311 + tmp_314 + tmp_319;
      real_t tmp_578 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_579 = tmp_355 + tmp_358 + tmp_363;
      real_t tmp_580 = tmp_377 + tmp_380 + tmp_385;
      real_t tmp_581 = tmp_399 + tmp_402 + tmp_407;
      real_t tmp_582 = tmp_421 + tmp_424 + tmp_429;
      real_t tmp_583 = tmp_443 + tmp_446 + tmp_451;
      real_t tmp_584 = tmp_465 + tmp_468 + tmp_473;
      real_t tmp_585 = tmp_487 + tmp_490 + tmp_495;
      real_t tmp_586 = tmp_509 + tmp_512 + tmp_517;
      real_t tmp_587 = tmp_531 + tmp_534 + tmp_539;
      real_t tmp_588 = tmp_78 + tmp_82 + tmp_88;
      real_t tmp_589 = 0.5*p_affine_13_0*tmp_66 + 0.5*p_affine_13_1*tmp_69 + 0.5*p_affine_13_2*tmp_72;
      real_t tmp_590 = tmp_112 + tmp_115 + tmp_120;
      real_t tmp_591 = tmp_134 + tmp_137 + tmp_142;
      real_t tmp_592 = tmp_156 + tmp_159 + tmp_164;
      real_t tmp_593 = tmp_178 + tmp_181 + tmp_186;
      real_t tmp_594 = tmp_200 + tmp_203 + tmp_208;
      real_t tmp_595 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_596 = tmp_244 + tmp_247 + tmp_252;
      real_t tmp_597 = tmp_266 + tmp_269 + tmp_274;
      real_t tmp_598 = tmp_288 + tmp_291 + tmp_296;
      real_t tmp_599 = tmp_310 + tmp_313 + tmp_318;
      real_t tmp_600 = tmp_332 + tmp_335 + tmp_340;
      real_t tmp_601 = tmp_354 + tmp_357 + tmp_362;
      real_t tmp_602 = tmp_376 + tmp_379 + tmp_384;
      real_t tmp_603 = tmp_398 + tmp_401 + tmp_406;
      real_t tmp_604 = tmp_420 + tmp_423 + tmp_428;
      real_t tmp_605 = tmp_442 + tmp_445 + tmp_450;
      real_t tmp_606 = tmp_464 + tmp_467 + tmp_472;
      real_t tmp_607 = tmp_486 + tmp_489 + tmp_494;
      real_t tmp_608 = tmp_508 + tmp_511 + tmp_516;
      real_t tmp_609 = tmp_530 + tmp_533 + tmp_538;
      real_t a_0_0 = tmp_103*(-tmp_101*tmp_91 + tmp_46*tmp_75 - tmp_91*tmp_95) + tmp_125*(tmp_110*tmp_75 - tmp_123*tmp_124 - tmp_123*tmp_95) + tmp_147*(tmp_132*tmp_75 - tmp_145*tmp_146 - tmp_145*tmp_95) + tmp_169*(tmp_154*tmp_75 - tmp_167*tmp_168 - tmp_167*tmp_95) + tmp_191*(tmp_176*tmp_75 - tmp_189*tmp_190 - tmp_189*tmp_95) + tmp_213*(tmp_198*tmp_75 - tmp_211*tmp_212 - tmp_211*tmp_95) + tmp_235*(tmp_220*tmp_75 - tmp_233*tmp_234 - tmp_233*tmp_95) + tmp_257*(tmp_242*tmp_75 - tmp_255*tmp_256 - tmp_255*tmp_95) + tmp_279*(tmp_264*tmp_75 - tmp_277*tmp_278 - tmp_277*tmp_95) + tmp_301*(tmp_286*tmp_75 - tmp_299*tmp_300 - tmp_299*tmp_95) + tmp_323*(tmp_308*tmp_75 - tmp_321*tmp_322 - tmp_321*tmp_95) + tmp_345*(tmp_330*tmp_75 - tmp_343*tmp_344 - tmp_343*tmp_95) + tmp_367*(tmp_352*tmp_75 - tmp_365*tmp_366 - tmp_365*tmp_95) + tmp_389*(tmp_374*tmp_75 - tmp_387*tmp_388 - tmp_387*tmp_95) + tmp_411*(tmp_396*tmp_75 - tmp_409*tmp_410 - tmp_409*tmp_95) + tmp_433*(tmp_418*tmp_75 - tmp_431*tmp_432 - tmp_431*tmp_95) + tmp_455*(tmp_440*tmp_75 - tmp_453*tmp_454 - tmp_453*tmp_95) + tmp_477*(tmp_462*tmp_75 - tmp_475*tmp_476 - tmp_475*tmp_95) + tmp_499*(tmp_484*tmp_75 - tmp_497*tmp_498 - tmp_497*tmp_95) + tmp_521*(tmp_506*tmp_75 - tmp_519*tmp_520 - tmp_519*tmp_95) + tmp_543*(tmp_528*tmp_75 - tmp_541*tmp_542 - tmp_541*tmp_95);
      real_t a_1_0 = tmp_103*(-tmp_101*tmp_544 + tmp_46*tmp_545 - tmp_544*tmp_95) + tmp_125*(tmp_110*tmp_545 - tmp_124*tmp_546 - tmp_546*tmp_95) + tmp_147*(tmp_132*tmp_545 - tmp_146*tmp_547 - tmp_547*tmp_95) + tmp_169*(tmp_154*tmp_545 - tmp_168*tmp_548 - tmp_548*tmp_95) + tmp_191*(tmp_176*tmp_545 - tmp_190*tmp_549 - tmp_549*tmp_95) + tmp_213*(tmp_198*tmp_545 - tmp_212*tmp_550 - tmp_550*tmp_95) + tmp_235*(tmp_220*tmp_545 - tmp_234*tmp_551 - tmp_551*tmp_95) + tmp_257*(tmp_242*tmp_545 - tmp_256*tmp_552 - tmp_552*tmp_95) + tmp_279*(tmp_264*tmp_545 - tmp_278*tmp_553 - tmp_553*tmp_95) + tmp_301*(tmp_286*tmp_545 - tmp_300*tmp_554 - tmp_554*tmp_95) + tmp_323*(tmp_308*tmp_545 - tmp_322*tmp_555 - tmp_555*tmp_95) + tmp_345*(tmp_330*tmp_545 - tmp_344*tmp_556 - tmp_556*tmp_95) + tmp_367*(tmp_352*tmp_545 - tmp_366*tmp_557 - tmp_557*tmp_95) + tmp_389*(tmp_374*tmp_545 - tmp_388*tmp_558 - tmp_558*tmp_95) + tmp_411*(tmp_396*tmp_545 - tmp_410*tmp_559 - tmp_559*tmp_95) + tmp_433*(tmp_418*tmp_545 - tmp_432*tmp_560 - tmp_560*tmp_95) + tmp_455*(tmp_440*tmp_545 - tmp_454*tmp_561 - tmp_561*tmp_95) + tmp_477*(tmp_462*tmp_545 - tmp_476*tmp_562 - tmp_562*tmp_95) + tmp_499*(tmp_484*tmp_545 - tmp_498*tmp_563 - tmp_563*tmp_95) + tmp_521*(tmp_506*tmp_545 - tmp_520*tmp_564 - tmp_564*tmp_95) + tmp_543*(tmp_528*tmp_545 - tmp_542*tmp_565 - tmp_565*tmp_95);
      real_t a_2_0 = tmp_103*(-tmp_101*tmp_566 + tmp_46*tmp_567 - tmp_566*tmp_95) + tmp_125*(tmp_110*tmp_567 - tmp_124*tmp_568 - tmp_568*tmp_95) + tmp_147*(tmp_132*tmp_567 - tmp_146*tmp_569 - tmp_569*tmp_95) + tmp_169*(tmp_154*tmp_567 - tmp_168*tmp_570 - tmp_570*tmp_95) + tmp_191*(tmp_176*tmp_567 - tmp_190*tmp_571 - tmp_571*tmp_95) + tmp_213*(tmp_198*tmp_567 - tmp_212*tmp_572 - tmp_572*tmp_95) + tmp_235*(tmp_220*tmp_567 - tmp_234*tmp_573 - tmp_573*tmp_95) + tmp_257*(tmp_242*tmp_567 - tmp_256*tmp_574 - tmp_574*tmp_95) + tmp_279*(tmp_264*tmp_567 - tmp_278*tmp_575 - tmp_575*tmp_95) + tmp_301*(tmp_286*tmp_567 - tmp_300*tmp_576 - tmp_576*tmp_95) + tmp_323*(tmp_308*tmp_567 - tmp_322*tmp_577 - tmp_577*tmp_95) + tmp_345*(tmp_330*tmp_567 - tmp_344*tmp_578 - tmp_578*tmp_95) + tmp_367*(tmp_352*tmp_567 - tmp_366*tmp_579 - tmp_579*tmp_95) + tmp_389*(tmp_374*tmp_567 - tmp_388*tmp_580 - tmp_580*tmp_95) + tmp_411*(tmp_396*tmp_567 - tmp_410*tmp_581 - tmp_581*tmp_95) + tmp_433*(tmp_418*tmp_567 - tmp_432*tmp_582 - tmp_582*tmp_95) + tmp_455*(tmp_440*tmp_567 - tmp_454*tmp_583 - tmp_583*tmp_95) + tmp_477*(tmp_462*tmp_567 - tmp_476*tmp_584 - tmp_584*tmp_95) + tmp_499*(tmp_484*tmp_567 - tmp_498*tmp_585 - tmp_585*tmp_95) + tmp_521*(tmp_506*tmp_567 - tmp_520*tmp_586 - tmp_586*tmp_95) + tmp_543*(tmp_528*tmp_567 - tmp_542*tmp_587 - tmp_587*tmp_95);
      real_t a_3_0 = tmp_103*(-tmp_101*tmp_588 + tmp_46*tmp_589 - tmp_588*tmp_95) + tmp_125*(tmp_110*tmp_589 - tmp_124*tmp_590 - tmp_590*tmp_95) + tmp_147*(tmp_132*tmp_589 - tmp_146*tmp_591 - tmp_591*tmp_95) + tmp_169*(tmp_154*tmp_589 - tmp_168*tmp_592 - tmp_592*tmp_95) + tmp_191*(tmp_176*tmp_589 - tmp_190*tmp_593 - tmp_593*tmp_95) + tmp_213*(tmp_198*tmp_589 - tmp_212*tmp_594 - tmp_594*tmp_95) + tmp_235*(tmp_220*tmp_589 - tmp_234*tmp_595 - tmp_595*tmp_95) + tmp_257*(tmp_242*tmp_589 - tmp_256*tmp_596 - tmp_596*tmp_95) + tmp_279*(tmp_264*tmp_589 - tmp_278*tmp_597 - tmp_597*tmp_95) + tmp_301*(tmp_286*tmp_589 - tmp_300*tmp_598 - tmp_598*tmp_95) + tmp_323*(tmp_308*tmp_589 - tmp_322*tmp_599 - tmp_599*tmp_95) + tmp_345*(tmp_330*tmp_589 - tmp_344*tmp_600 - tmp_600*tmp_95) + tmp_367*(tmp_352*tmp_589 - tmp_366*tmp_601 - tmp_601*tmp_95) + tmp_389*(tmp_374*tmp_589 - tmp_388*tmp_602 - tmp_602*tmp_95) + tmp_411*(tmp_396*tmp_589 - tmp_410*tmp_603 - tmp_603*tmp_95) + tmp_433*(tmp_418*tmp_589 - tmp_432*tmp_604 - tmp_604*tmp_95) + tmp_455*(tmp_440*tmp_589 - tmp_454*tmp_605 - tmp_605*tmp_95) + tmp_477*(tmp_462*tmp_589 - tmp_476*tmp_606 - tmp_606*tmp_95) + tmp_499*(tmp_484*tmp_589 - tmp_498*tmp_607 - tmp_607*tmp_95) + tmp_521*(tmp_506*tmp_589 - tmp_520*tmp_608 - tmp_608*tmp_95) + tmp_543*(tmp_528*tmp_589 - tmp_542*tmp_609 - tmp_609*tmp_95);
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
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_2;
      real_t tmp_8 = p_affine_3_2 + tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_7;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_4;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_11*tmp_14;
      real_t tmp_17 = tmp_4*tmp_8;
      real_t tmp_18 = tmp_13*tmp_3;
      real_t tmp_19 = 1.0 / (-tmp_1*tmp_16 + tmp_1*tmp_9 + tmp_10*tmp_12 - tmp_10*tmp_18 + tmp_15*tmp_5 - tmp_17*tmp_5);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_6;
      real_t tmp_26 = -tmp_1*tmp_11 + tmp_13*tmp_5;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_14 + tmp_10*tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_8 - tmp_10*tmp_13;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_12 - tmp_18;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = tmp_15 - tmp_17;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_3 + tmp_14*tmp_5;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = tmp_10*tmp_11 - tmp_5*tmp_8;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_53 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_54 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_55 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_58 = 5.0*std::pow((std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)*std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)) + (std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)*std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)) + (std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)*std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)), 0.25);
      real_t tmp_59 = 0.0068572537431980923*tmp_58*(tmp_11*(tmp_54 - 1.0/4.0) + tmp_3*(tmp_53 - 1.0/4.0) + tmp_5*(tmp_52 - 1.0/4.0));
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_61 = tmp_6*tmp_60;
      real_t tmp_62 = tmp_26*tmp_60;
      real_t tmp_63 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_64 = tmp_28*tmp_63;
      real_t tmp_65 = tmp_35*tmp_63;
      real_t tmp_66 = tmp_37*tmp_60;
      real_t tmp_67 = tmp_39*tmp_63;
      real_t tmp_68 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_69 = tmp_41*tmp_68;
      real_t tmp_70 = tmp_48*tmp_68;
      real_t tmp_71 = tmp_50*tmp_68;
      real_t tmp_72 = tmp_66 + tmp_67 + tmp_71;
      real_t tmp_73 = tmp_62 + tmp_65 + tmp_70;
      real_t tmp_74 = tmp_61 + tmp_64 + tmp_69;
      real_t tmp_75 = 0.037198804536718075*tmp_58*(tmp_11*(tmp_74 - 1.0/4.0) + tmp_3*(tmp_73 - 1.0/4.0) + tmp_5*(tmp_72 - 1.0/4.0));
      real_t tmp_76 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_77 = tmp_6*tmp_76;
      real_t tmp_78 = tmp_26*tmp_76;
      real_t tmp_79 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_80 = tmp_28*tmp_79;
      real_t tmp_81 = tmp_35*tmp_79;
      real_t tmp_82 = tmp_37*tmp_76;
      real_t tmp_83 = tmp_39*tmp_79;
      real_t tmp_84 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_85 = tmp_41*tmp_84;
      real_t tmp_86 = tmp_48*tmp_84;
      real_t tmp_87 = tmp_50*tmp_84;
      real_t tmp_88 = tmp_82 + tmp_83 + tmp_87;
      real_t tmp_89 = tmp_78 + tmp_81 + tmp_86;
      real_t tmp_90 = tmp_77 + tmp_80 + tmp_85;
      real_t tmp_91 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_90 - 1.0/4.0) + tmp_3*(tmp_89 - 1.0/4.0) + tmp_5*(tmp_88 - 1.0/4.0));
      real_t tmp_92 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_93 = tmp_6*tmp_92;
      real_t tmp_94 = tmp_26*tmp_92;
      real_t tmp_95 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_96 = tmp_28*tmp_95;
      real_t tmp_97 = tmp_35*tmp_95;
      real_t tmp_98 = tmp_37*tmp_92;
      real_t tmp_99 = tmp_39*tmp_95;
      real_t tmp_100 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_101 = tmp_100*tmp_41;
      real_t tmp_102 = tmp_100*tmp_48;
      real_t tmp_103 = tmp_100*tmp_50;
      real_t tmp_104 = tmp_103 + tmp_98 + tmp_99;
      real_t tmp_105 = tmp_102 + tmp_94 + tmp_97;
      real_t tmp_106 = tmp_101 + tmp_93 + tmp_96;
      real_t tmp_107 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_106 - 1.0/4.0) + tmp_3*(tmp_105 - 1.0/4.0) + tmp_5*(tmp_104 - 1.0/4.0));
      real_t tmp_108 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_109 = tmp_108*tmp_6;
      real_t tmp_110 = tmp_108*tmp_26;
      real_t tmp_111 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_112 = tmp_111*tmp_28;
      real_t tmp_113 = tmp_111*tmp_35;
      real_t tmp_114 = tmp_108*tmp_37;
      real_t tmp_115 = tmp_111*tmp_39;
      real_t tmp_116 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_117 = tmp_116*tmp_41;
      real_t tmp_118 = tmp_116*tmp_48;
      real_t tmp_119 = tmp_116*tmp_50;
      real_t tmp_120 = tmp_114 + tmp_115 + tmp_119;
      real_t tmp_121 = tmp_110 + tmp_113 + tmp_118;
      real_t tmp_122 = tmp_109 + tmp_112 + tmp_117;
      real_t tmp_123 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_122 - 1.0/4.0) + tmp_3*(tmp_121 - 1.0/4.0) + tmp_5*(tmp_120 - 1.0/4.0));
      real_t tmp_124 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_125 = tmp_124*tmp_6;
      real_t tmp_126 = tmp_124*tmp_26;
      real_t tmp_127 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_128 = tmp_127*tmp_28;
      real_t tmp_129 = tmp_127*tmp_35;
      real_t tmp_130 = tmp_124*tmp_37;
      real_t tmp_131 = tmp_127*tmp_39;
      real_t tmp_132 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_133 = tmp_132*tmp_41;
      real_t tmp_134 = tmp_132*tmp_48;
      real_t tmp_135 = tmp_132*tmp_50;
      real_t tmp_136 = tmp_130 + tmp_131 + tmp_135;
      real_t tmp_137 = tmp_126 + tmp_129 + tmp_134;
      real_t tmp_138 = tmp_125 + tmp_128 + tmp_133;
      real_t tmp_139 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_138 - 1.0/4.0) + tmp_3*(tmp_137 - 1.0/4.0) + tmp_5*(tmp_136 - 1.0/4.0));
      real_t tmp_140 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_141 = tmp_140*tmp_6;
      real_t tmp_142 = tmp_140*tmp_26;
      real_t tmp_143 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_144 = tmp_143*tmp_28;
      real_t tmp_145 = tmp_143*tmp_35;
      real_t tmp_146 = tmp_140*tmp_37;
      real_t tmp_147 = tmp_143*tmp_39;
      real_t tmp_148 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_149 = tmp_148*tmp_41;
      real_t tmp_150 = tmp_148*tmp_48;
      real_t tmp_151 = tmp_148*tmp_50;
      real_t tmp_152 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_153 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_154 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_155 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_154 - 1.0/4.0) + tmp_3*(tmp_153 - 1.0/4.0) + tmp_5*(tmp_152 - 1.0/4.0));
      real_t tmp_156 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_157 = tmp_156*tmp_6;
      real_t tmp_158 = tmp_156*tmp_26;
      real_t tmp_159 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_160 = tmp_159*tmp_28;
      real_t tmp_161 = tmp_159*tmp_35;
      real_t tmp_162 = tmp_156*tmp_37;
      real_t tmp_163 = tmp_159*tmp_39;
      real_t tmp_164 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_165 = tmp_164*tmp_41;
      real_t tmp_166 = tmp_164*tmp_48;
      real_t tmp_167 = tmp_164*tmp_50;
      real_t tmp_168 = tmp_162 + tmp_163 + tmp_167;
      real_t tmp_169 = tmp_158 + tmp_161 + tmp_166;
      real_t tmp_170 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_171 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_170 - 1.0/4.0) + tmp_3*(tmp_169 - 1.0/4.0) + tmp_5*(tmp_168 - 1.0/4.0));
      real_t tmp_172 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_173 = tmp_172*tmp_6;
      real_t tmp_174 = tmp_172*tmp_26;
      real_t tmp_175 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_176 = tmp_175*tmp_28;
      real_t tmp_177 = tmp_175*tmp_35;
      real_t tmp_178 = tmp_172*tmp_37;
      real_t tmp_179 = tmp_175*tmp_39;
      real_t tmp_180 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_181 = tmp_180*tmp_41;
      real_t tmp_182 = tmp_180*tmp_48;
      real_t tmp_183 = tmp_180*tmp_50;
      real_t tmp_184 = tmp_178 + tmp_179 + tmp_183;
      real_t tmp_185 = tmp_174 + tmp_177 + tmp_182;
      real_t tmp_186 = tmp_173 + tmp_176 + tmp_181;
      real_t tmp_187 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_186 - 1.0/4.0) + tmp_3*(tmp_185 - 1.0/4.0) + tmp_5*(tmp_184 - 1.0/4.0));
      real_t tmp_188 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_189 = tmp_188*tmp_6;
      real_t tmp_190 = tmp_188*tmp_26;
      real_t tmp_191 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_192 = tmp_191*tmp_28;
      real_t tmp_193 = tmp_191*tmp_35;
      real_t tmp_194 = tmp_188*tmp_37;
      real_t tmp_195 = tmp_191*tmp_39;
      real_t tmp_196 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_197 = tmp_196*tmp_41;
      real_t tmp_198 = tmp_196*tmp_48;
      real_t tmp_199 = tmp_196*tmp_50;
      real_t tmp_200 = tmp_194 + tmp_195 + tmp_199;
      real_t tmp_201 = tmp_190 + tmp_193 + tmp_198;
      real_t tmp_202 = tmp_189 + tmp_192 + tmp_197;
      real_t tmp_203 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_202 - 1.0/4.0) + tmp_3*(tmp_201 - 1.0/4.0) + tmp_5*(tmp_200 - 1.0/4.0));
      real_t tmp_204 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_205 = tmp_204*tmp_6;
      real_t tmp_206 = tmp_204*tmp_26;
      real_t tmp_207 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_208 = tmp_207*tmp_28;
      real_t tmp_209 = tmp_207*tmp_35;
      real_t tmp_210 = tmp_204*tmp_37;
      real_t tmp_211 = tmp_207*tmp_39;
      real_t tmp_212 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_213 = tmp_212*tmp_41;
      real_t tmp_214 = tmp_212*tmp_48;
      real_t tmp_215 = tmp_212*tmp_50;
      real_t tmp_216 = tmp_210 + tmp_211 + tmp_215;
      real_t tmp_217 = tmp_206 + tmp_209 + tmp_214;
      real_t tmp_218 = tmp_205 + tmp_208 + tmp_213;
      real_t tmp_219 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_218 - 1.0/4.0) + tmp_3*(tmp_217 - 1.0/4.0) + tmp_5*(tmp_216 - 1.0/4.0));
      real_t tmp_220 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_221 = tmp_220*tmp_6;
      real_t tmp_222 = tmp_220*tmp_26;
      real_t tmp_223 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_224 = tmp_223*tmp_28;
      real_t tmp_225 = tmp_223*tmp_35;
      real_t tmp_226 = tmp_220*tmp_37;
      real_t tmp_227 = tmp_223*tmp_39;
      real_t tmp_228 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_229 = tmp_228*tmp_41;
      real_t tmp_230 = tmp_228*tmp_48;
      real_t tmp_231 = tmp_228*tmp_50;
      real_t tmp_232 = tmp_226 + tmp_227 + tmp_231;
      real_t tmp_233 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_234 = tmp_221 + tmp_224 + tmp_229;
      real_t tmp_235 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_234 - 1.0/4.0) + tmp_3*(tmp_233 - 1.0/4.0) + tmp_5*(tmp_232 - 1.0/4.0));
      real_t tmp_236 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_237 = tmp_236*tmp_6;
      real_t tmp_238 = tmp_236*tmp_26;
      real_t tmp_239 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_240 = tmp_239*tmp_28;
      real_t tmp_241 = tmp_239*tmp_35;
      real_t tmp_242 = tmp_236*tmp_37;
      real_t tmp_243 = tmp_239*tmp_39;
      real_t tmp_244 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_245 = tmp_244*tmp_41;
      real_t tmp_246 = tmp_244*tmp_48;
      real_t tmp_247 = tmp_244*tmp_50;
      real_t tmp_248 = tmp_242 + tmp_243 + tmp_247;
      real_t tmp_249 = tmp_238 + tmp_241 + tmp_246;
      real_t tmp_250 = tmp_237 + tmp_240 + tmp_245;
      real_t tmp_251 = 0.042507265838595799*tmp_58*(tmp_11*(tmp_250 - 1.0/4.0) + tmp_3*(tmp_249 - 1.0/4.0) + tmp_5*(tmp_248 - 1.0/4.0));
      real_t tmp_252 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_6;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = tmp_258 + tmp_259 + tmp_263;
      real_t tmp_265 = tmp_254 + tmp_257 + tmp_262;
      real_t tmp_266 = tmp_253 + tmp_256 + tmp_261;
      real_t tmp_267 = 0.020848748529055869*tmp_58*(tmp_11*(tmp_266 - 1.0/4.0) + tmp_3*(tmp_265 - 1.0/4.0) + tmp_5*(tmp_264 - 1.0/4.0));
      real_t tmp_268 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_269 = tmp_268*tmp_6;
      real_t tmp_270 = tmp_26*tmp_268;
      real_t tmp_271 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_272 = tmp_271*tmp_28;
      real_t tmp_273 = tmp_271*tmp_35;
      real_t tmp_274 = tmp_268*tmp_37;
      real_t tmp_275 = tmp_271*tmp_39;
      real_t tmp_276 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_277 = tmp_276*tmp_41;
      real_t tmp_278 = tmp_276*tmp_48;
      real_t tmp_279 = tmp_276*tmp_50;
      real_t tmp_280 = tmp_274 + tmp_275 + tmp_279;
      real_t tmp_281 = tmp_270 + tmp_273 + tmp_278;
      real_t tmp_282 = tmp_269 + tmp_272 + tmp_277;
      real_t tmp_283 = 0.0068572537431980923*tmp_58*(tmp_11*(tmp_282 - 1.0/4.0) + tmp_3*(tmp_281 - 1.0/4.0) + tmp_5*(tmp_280 - 1.0/4.0));
      real_t tmp_284 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_285 = tmp_284*tmp_6;
      real_t tmp_286 = tmp_26*tmp_284;
      real_t tmp_287 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_288 = tmp_28*tmp_287;
      real_t tmp_289 = tmp_287*tmp_35;
      real_t tmp_290 = tmp_284*tmp_37;
      real_t tmp_291 = tmp_287*tmp_39;
      real_t tmp_292 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_293 = tmp_292*tmp_41;
      real_t tmp_294 = tmp_292*tmp_48;
      real_t tmp_295 = tmp_292*tmp_50;
      real_t tmp_296 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_297 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_298 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_299 = 0.037198804536718075*tmp_58*(tmp_11*(tmp_298 - 1.0/4.0) + tmp_3*(tmp_297 - 1.0/4.0) + tmp_5*(tmp_296 - 1.0/4.0));
      real_t tmp_300 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_301 = tmp_300*tmp_6;
      real_t tmp_302 = tmp_26*tmp_300;
      real_t tmp_303 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_304 = tmp_28*tmp_303;
      real_t tmp_305 = tmp_303*tmp_35;
      real_t tmp_306 = tmp_300*tmp_37;
      real_t tmp_307 = tmp_303*tmp_39;
      real_t tmp_308 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_309 = tmp_308*tmp_41;
      real_t tmp_310 = tmp_308*tmp_48;
      real_t tmp_311 = tmp_308*tmp_50;
      real_t tmp_312 = tmp_306 + tmp_307 + tmp_311;
      real_t tmp_313 = tmp_302 + tmp_305 + tmp_310;
      real_t tmp_314 = tmp_301 + tmp_304 + tmp_309;
      real_t tmp_315 = 0.042507265838595799*tmp_58*(tmp_11*(tmp_314 - 1.0/4.0) + tmp_3*(tmp_313 - 1.0/4.0) + tmp_5*(tmp_312 - 1.0/4.0));
      real_t tmp_316 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_317 = tmp_316*tmp_6;
      real_t tmp_318 = tmp_26*tmp_316;
      real_t tmp_319 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_320 = tmp_28*tmp_319;
      real_t tmp_321 = tmp_319*tmp_35;
      real_t tmp_322 = tmp_316*tmp_37;
      real_t tmp_323 = tmp_319*tmp_39;
      real_t tmp_324 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_325 = tmp_324*tmp_41;
      real_t tmp_326 = tmp_324*tmp_48;
      real_t tmp_327 = tmp_324*tmp_50;
      real_t tmp_328 = tmp_322 + tmp_323 + tmp_327;
      real_t tmp_329 = tmp_318 + tmp_321 + tmp_326;
      real_t tmp_330 = tmp_317 + tmp_320 + tmp_325;
      real_t tmp_331 = 0.0068572537431980923*tmp_58*(tmp_11*(tmp_330 - 1.0/4.0) + tmp_3*(tmp_329 - 1.0/4.0) + tmp_5*(tmp_328 - 1.0/4.0));
      real_t tmp_332 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_333 = tmp_332*tmp_6;
      real_t tmp_334 = tmp_26*tmp_332;
      real_t tmp_335 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_336 = tmp_28*tmp_335;
      real_t tmp_337 = tmp_335*tmp_35;
      real_t tmp_338 = tmp_332*tmp_37;
      real_t tmp_339 = tmp_335*tmp_39;
      real_t tmp_340 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_341 = tmp_340*tmp_41;
      real_t tmp_342 = tmp_340*tmp_48;
      real_t tmp_343 = tmp_340*tmp_50;
      real_t tmp_344 = tmp_338 + tmp_339 + tmp_343;
      real_t tmp_345 = tmp_334 + tmp_337 + tmp_342;
      real_t tmp_346 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_347 = 0.037198804536718075*tmp_58*(tmp_11*(tmp_346 - 1.0/4.0) + tmp_3*(tmp_345 - 1.0/4.0) + tmp_5*(tmp_344 - 1.0/4.0));
      real_t tmp_348 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_349 = tmp_348*tmp_6;
      real_t tmp_350 = tmp_26*tmp_348;
      real_t tmp_351 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_352 = tmp_28*tmp_351;
      real_t tmp_353 = tmp_35*tmp_351;
      real_t tmp_354 = tmp_348*tmp_37;
      real_t tmp_355 = tmp_351*tmp_39;
      real_t tmp_356 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_357 = tmp_356*tmp_41;
      real_t tmp_358 = tmp_356*tmp_48;
      real_t tmp_359 = tmp_356*tmp_50;
      real_t tmp_360 = tmp_354 + tmp_355 + tmp_359;
      real_t tmp_361 = tmp_350 + tmp_353 + tmp_358;
      real_t tmp_362 = tmp_349 + tmp_352 + tmp_357;
      real_t tmp_363 = 0.042507265838595799*tmp_58*(tmp_11*(tmp_362 - 1.0/4.0) + tmp_3*(tmp_361 - 1.0/4.0) + tmp_5*(tmp_360 - 1.0/4.0));
      real_t tmp_364 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_365 = tmp_364*tmp_6;
      real_t tmp_366 = tmp_26*tmp_364;
      real_t tmp_367 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_368 = tmp_28*tmp_367;
      real_t tmp_369 = tmp_35*tmp_367;
      real_t tmp_370 = tmp_364*tmp_37;
      real_t tmp_371 = tmp_367*tmp_39;
      real_t tmp_372 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_373 = tmp_372*tmp_41;
      real_t tmp_374 = tmp_372*tmp_48;
      real_t tmp_375 = tmp_372*tmp_50;
      real_t tmp_376 = tmp_370 + tmp_371 + tmp_375;
      real_t tmp_377 = tmp_366 + tmp_369 + tmp_374;
      real_t tmp_378 = tmp_365 + tmp_368 + tmp_373;
      real_t tmp_379 = 0.019202922745021479*tmp_58*(tmp_11*(tmp_378 - 1.0/4.0) + tmp_3*(tmp_377 - 1.0/4.0) + tmp_5*(tmp_376 - 1.0/4.0));
      real_t a_0_0 = tmp_107*(-tmp_101 - tmp_102 - tmp_103 - tmp_93 - tmp_94 - tmp_96 - tmp_97 - tmp_98 - tmp_99 + 1) + tmp_123*(-tmp_109 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_115 - tmp_117 - tmp_118 - tmp_119 + 1) + tmp_139*(-tmp_125 - tmp_126 - tmp_128 - tmp_129 - tmp_130 - tmp_131 - tmp_133 - tmp_134 - tmp_135 + 1) + tmp_155*(-tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1) + tmp_171*(-tmp_157 - tmp_158 - tmp_160 - tmp_161 - tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 + 1) + tmp_187*(-tmp_173 - tmp_174 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 + 1) + tmp_203*(-tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_197 - tmp_198 - tmp_199 + 1) + tmp_219*(-tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 - tmp_211 - tmp_213 - tmp_214 - tmp_215 + 1) + tmp_235*(-tmp_221 - tmp_222 - tmp_224 - tmp_225 - tmp_226 - tmp_227 - tmp_229 - tmp_230 - tmp_231 + 1) + tmp_251*(-tmp_237 - tmp_238 - tmp_240 - tmp_241 - tmp_242 - tmp_243 - tmp_245 - tmp_246 - tmp_247 + 1) + tmp_267*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_283*(-tmp_269 - tmp_270 - tmp_272 - tmp_273 - tmp_274 - tmp_275 - tmp_277 - tmp_278 - tmp_279 + 1) + tmp_299*(-tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1) + tmp_315*(-tmp_301 - tmp_302 - tmp_304 - tmp_305 - tmp_306 - tmp_307 - tmp_309 - tmp_310 - tmp_311 + 1) + tmp_331*(-tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1) + tmp_347*(-tmp_333 - tmp_334 - tmp_336 - tmp_337 - tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 + 1) + tmp_363*(-tmp_349 - tmp_350 - tmp_352 - tmp_353 - tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 + 1) + tmp_379*(-tmp_365 - tmp_366 - tmp_368 - tmp_369 - tmp_370 - tmp_371 - tmp_373 - tmp_374 - tmp_375 + 1) + tmp_59*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_75*(-tmp_61 - tmp_62 - tmp_64 - tmp_65 - tmp_66 - tmp_67 - tmp_69 - tmp_70 - tmp_71 + 1) + tmp_91*(-tmp_77 - tmp_78 - tmp_80 - tmp_81 - tmp_82 - tmp_83 - tmp_85 - tmp_86 - tmp_87 + 1);
      real_t a_1_0 = tmp_104*tmp_107 + tmp_120*tmp_123 + tmp_136*tmp_139 + tmp_152*tmp_155 + tmp_168*tmp_171 + tmp_184*tmp_187 + tmp_200*tmp_203 + tmp_216*tmp_219 + tmp_232*tmp_235 + tmp_248*tmp_251 + tmp_264*tmp_267 + tmp_280*tmp_283 + tmp_296*tmp_299 + tmp_312*tmp_315 + tmp_328*tmp_331 + tmp_344*tmp_347 + tmp_360*tmp_363 + tmp_376*tmp_379 + tmp_52*tmp_59 + tmp_72*tmp_75 + tmp_88*tmp_91;
      real_t a_2_0 = tmp_105*tmp_107 + tmp_121*tmp_123 + tmp_137*tmp_139 + tmp_153*tmp_155 + tmp_169*tmp_171 + tmp_185*tmp_187 + tmp_201*tmp_203 + tmp_217*tmp_219 + tmp_233*tmp_235 + tmp_249*tmp_251 + tmp_265*tmp_267 + tmp_281*tmp_283 + tmp_297*tmp_299 + tmp_313*tmp_315 + tmp_329*tmp_331 + tmp_345*tmp_347 + tmp_361*tmp_363 + tmp_377*tmp_379 + tmp_53*tmp_59 + tmp_73*tmp_75 + tmp_89*tmp_91;
      real_t a_3_0 = tmp_106*tmp_107 + tmp_122*tmp_123 + tmp_138*tmp_139 + tmp_154*tmp_155 + tmp_170*tmp_171 + tmp_186*tmp_187 + tmp_202*tmp_203 + tmp_218*tmp_219 + tmp_234*tmp_235 + tmp_250*tmp_251 + tmp_266*tmp_267 + tmp_282*tmp_283 + tmp_298*tmp_299 + tmp_314*tmp_315 + tmp_330*tmp_331 + tmp_346*tmp_347 + tmp_362*tmp_363 + tmp_378*tmp_379 + tmp_54*tmp_59 + tmp_74*tmp_75 + tmp_90*tmp_91;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }


};




class EGNIPGVectorLaplaceFormP1E_2 : public hyteg::dg::DGForm
{
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
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
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
      real_t tmp_22 = tmp_10*tmp_21 + tmp_13*tmp_20 + tmp_19*tmp_9;
      real_t tmp_23 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_5);
      real_t tmp_24 = tmp_18*(tmp_1*tmp_9 - tmp_10*tmp_14);
      real_t tmp_25 = tmp_18*(tmp_13*tmp_14 - tmp_5*tmp_9);
      real_t tmp_26 = tmp_10*tmp_25 + tmp_13*tmp_24 + tmp_23*tmp_9;
      real_t tmp_27 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_6);
      real_t tmp_28 = tmp_18*(tmp_10*tmp_11 - tmp_6*tmp_9);
      real_t tmp_29 = tmp_18*(-tmp_11*tmp_13 + tmp_3*tmp_9);
      real_t tmp_30 = tmp_10*tmp_29 + tmp_13*tmp_28 + tmp_27*tmp_9;
      real_t tmp_31 = p_affine_0_0*p_affine_1_1;
      real_t tmp_32 = p_affine_0_0*p_affine_1_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_2;
      real_t tmp_34 = p_affine_0_1*p_affine_1_0;
      real_t tmp_35 = p_affine_0_1*p_affine_1_2;
      real_t tmp_36 = p_affine_2_2*p_affine_3_0;
      real_t tmp_37 = p_affine_0_2*p_affine_1_0;
      real_t tmp_38 = p_affine_0_2*p_affine_1_1;
      real_t tmp_39 = p_affine_2_0*p_affine_3_1;
      real_t tmp_40 = p_affine_2_2*p_affine_3_1;
      real_t tmp_41 = p_affine_2_0*p_affine_3_2;
      real_t tmp_42 = p_affine_2_1*p_affine_3_0;
      real_t tmp_43 = std::abs(p_affine_0_0*tmp_33 - p_affine_0_0*tmp_40 + p_affine_0_1*tmp_36 - p_affine_0_1*tmp_41 + p_affine_0_2*tmp_39 - p_affine_0_2*tmp_42 - p_affine_1_0*tmp_33 + p_affine_1_0*tmp_40 - p_affine_1_1*tmp_36 + p_affine_1_1*tmp_41 - p_affine_1_2*tmp_39 + p_affine_1_2*tmp_42 + p_affine_2_0*tmp_35 - p_affine_2_0*tmp_38 - p_affine_2_1*tmp_32 + p_affine_2_1*tmp_37 + p_affine_2_2*tmp_31 - p_affine_2_2*tmp_34 - p_affine_3_0*tmp_35 + p_affine_3_0*tmp_38 + p_affine_3_1*tmp_32 - p_affine_3_1*tmp_37 - p_affine_3_2*tmp_31 + p_affine_3_2*tmp_34);
      real_t tmp_44 = tmp_43*(tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25) + tmp_30*(-tmp_27 - tmp_28 - tmp_29));
      real_t tmp_45 = tmp_43*(tmp_21*tmp_22 + tmp_25*tmp_26 + tmp_29*tmp_30);
      real_t tmp_46 = tmp_43*(tmp_20*tmp_22 + tmp_24*tmp_26 + tmp_28*tmp_30);
      real_t tmp_47 = tmp_43*(tmp_19*tmp_22 + tmp_23*tmp_26 + tmp_27*tmp_30);
      real_t a_0_0 = 0.1666666666666668*tmp_44;
      real_t a_1_0 = 0.1666666666666668*tmp_45;
      real_t a_2_0 = 0.1666666666666668*tmp_46;
      real_t a_3_0 = 0.1666666666666668*tmp_47;
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

         real_t tmp_0 = -p_affine_0_2;
      real_t tmp_1 = p_affine_1_2 + tmp_0;
      real_t tmp_2 = -p_affine_8_2;
      real_t tmp_3 = p_affine_9_2 + tmp_2;
      real_t tmp_4 = p_affine_10_2 + tmp_2;
      real_t tmp_5 = p_affine_8_2 + tmp_0;
      real_t tmp_6 = 0.031405749086161582*tmp_3 + 0.93718850182767688*tmp_4 + tmp_5;
      real_t tmp_7 = -p_affine_0_0;
      real_t tmp_8 = p_affine_2_0 + tmp_7;
      real_t tmp_9 = -p_affine_0_1;
      real_t tmp_10 = p_affine_3_1 + tmp_9;
      real_t tmp_11 = p_affine_3_0 + tmp_7;
      real_t tmp_12 = p_affine_2_1 + tmp_9;
      real_t tmp_13 = p_affine_1_0 + tmp_7;
      real_t tmp_14 = p_affine_3_2 + tmp_0;
      real_t tmp_15 = tmp_12*tmp_14;
      real_t tmp_16 = tmp_1*tmp_10;
      real_t tmp_17 = p_affine_1_1 + tmp_9;
      real_t tmp_18 = p_affine_2_2 + tmp_0;
      real_t tmp_19 = tmp_17*tmp_18;
      real_t tmp_20 = tmp_10*tmp_18;
      real_t tmp_21 = tmp_14*tmp_17;
      real_t tmp_22 = tmp_1*tmp_12;
      real_t tmp_23 = 1.0 / (tmp_11*tmp_19 - tmp_11*tmp_22 + tmp_13*tmp_15 - tmp_13*tmp_20 + tmp_16*tmp_8 - tmp_21*tmp_8);
      real_t tmp_24 = tmp_23*(tmp_10*tmp_8 - tmp_11*tmp_12);
      real_t tmp_25 = tmp_24*tmp_6;
      real_t tmp_26 = -p_affine_8_1;
      real_t tmp_27 = p_affine_9_1 + tmp_26;
      real_t tmp_28 = p_affine_10_1 + tmp_26;
      real_t tmp_29 = p_affine_8_1 + tmp_9;
      real_t tmp_30 = 0.031405749086161582*tmp_27 + 0.93718850182767688*tmp_28 + tmp_29;
      real_t tmp_31 = tmp_23*(tmp_11*tmp_18 - tmp_14*tmp_8);
      real_t tmp_32 = tmp_30*tmp_31;
      real_t tmp_33 = -p_affine_8_0;
      real_t tmp_34 = p_affine_9_0 + tmp_33;
      real_t tmp_35 = p_affine_10_0 + tmp_33;
      real_t tmp_36 = p_affine_8_0 + tmp_7;
      real_t tmp_37 = 0.031405749086161582*tmp_34 + 0.93718850182767688*tmp_35 + tmp_36;
      real_t tmp_38 = tmp_23*(tmp_15 - tmp_20);
      real_t tmp_39 = tmp_37*tmp_38;
      real_t tmp_40 = tmp_25 + tmp_32 + tmp_39;
      real_t tmp_41 = tmp_23*(-tmp_10*tmp_13 + tmp_11*tmp_17);
      real_t tmp_42 = tmp_41*tmp_6;
      real_t tmp_43 = tmp_23*(-tmp_1*tmp_11 + tmp_13*tmp_14);
      real_t tmp_44 = tmp_30*tmp_43;
      real_t tmp_45 = tmp_23*(tmp_16 - tmp_21);
      real_t tmp_46 = tmp_37*tmp_45;
      real_t tmp_47 = tmp_42 + tmp_44 + tmp_46;
      real_t tmp_48 = tmp_23*(tmp_12*tmp_13 - tmp_17*tmp_8);
      real_t tmp_49 = tmp_48*tmp_6;
      real_t tmp_50 = tmp_23*(tmp_1*tmp_8 - tmp_13*tmp_18);
      real_t tmp_51 = tmp_30*tmp_50;
      real_t tmp_52 = tmp_23*(tmp_19 - tmp_22);
      real_t tmp_53 = tmp_37*tmp_52;
      real_t tmp_54 = tmp_49 + tmp_51 + tmp_53;
      real_t tmp_55 = tmp_1*(tmp_40 - 1.0/4.0) + tmp_14*(tmp_54 - 1.0/4.0) + tmp_18*(tmp_47 - 1.0/4.0);
      real_t tmp_56 = 0.5*p_affine_13_0*(-tmp_38 - tmp_45 - tmp_52) + 0.5*p_affine_13_1*(-tmp_31 - tmp_43 - tmp_50) + 0.5*p_affine_13_2*(-tmp_24 - tmp_41 - tmp_48);
      real_t tmp_57 = -tmp_25 - tmp_32 - tmp_39 - tmp_42 - tmp_44 - tmp_46 - tmp_49 - tmp_51 - tmp_53 + 1;
      real_t tmp_58 = 0.5*p_affine_13_0*(tmp_1*tmp_38 + tmp_14*tmp_52 + tmp_18*tmp_45) + 0.5*p_affine_13_1*(tmp_1*tmp_31 + tmp_14*tmp_50 + tmp_18*tmp_43) + 0.5*p_affine_13_2*(tmp_1*tmp_24 + tmp_14*tmp_48 + tmp_18*tmp_41);
      real_t tmp_59 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_60 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_61 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_62 = (std::abs(tmp_28*tmp_60 - tmp_35*tmp_59)*std::abs(tmp_28*tmp_60 - tmp_35*tmp_59)) + (std::abs(tmp_28*tmp_61 - tmp_4*tmp_59)*std::abs(tmp_28*tmp_61 - tmp_4*tmp_59)) + (std::abs(tmp_35*tmp_61 - tmp_4*tmp_60)*std::abs(tmp_35*tmp_61 - tmp_4*tmp_60));
      real_t tmp_63 = 5.0*std::pow(tmp_62, -0.25);
      real_t tmp_64 = tmp_55*tmp_63;
      real_t tmp_65 = 1.0*std::pow(tmp_62, 1.0/2.0);
      real_t tmp_66 = 0.0068572537431980923*tmp_65;
      real_t tmp_67 = 0.19601935860219369*tmp_3 + 0.60796128279561268*tmp_4 + tmp_5;
      real_t tmp_68 = tmp_24*tmp_67;
      real_t tmp_69 = 0.19601935860219369*tmp_27 + 0.60796128279561268*tmp_28 + tmp_29;
      real_t tmp_70 = tmp_31*tmp_69;
      real_t tmp_71 = 0.19601935860219369*tmp_34 + 0.60796128279561268*tmp_35 + tmp_36;
      real_t tmp_72 = tmp_38*tmp_71;
      real_t tmp_73 = tmp_68 + tmp_70 + tmp_72;
      real_t tmp_74 = tmp_41*tmp_67;
      real_t tmp_75 = tmp_43*tmp_69;
      real_t tmp_76 = tmp_45*tmp_71;
      real_t tmp_77 = tmp_74 + tmp_75 + tmp_76;
      real_t tmp_78 = tmp_48*tmp_67;
      real_t tmp_79 = tmp_50*tmp_69;
      real_t tmp_80 = tmp_52*tmp_71;
      real_t tmp_81 = tmp_78 + tmp_79 + tmp_80;
      real_t tmp_82 = tmp_1*(tmp_73 - 1.0/4.0) + tmp_14*(tmp_81 - 1.0/4.0) + tmp_18*(tmp_77 - 1.0/4.0);
      real_t tmp_83 = -tmp_68 - tmp_70 - tmp_72 - tmp_74 - tmp_75 - tmp_76 - tmp_78 - tmp_79 - tmp_80 + 1;
      real_t tmp_84 = tmp_63*tmp_82;
      real_t tmp_85 = 0.037198804536718075*tmp_65;
      real_t tmp_86 = 0.37605877282253791*tmp_3 + 0.039308471900058539*tmp_4 + tmp_5;
      real_t tmp_87 = tmp_24*tmp_86;
      real_t tmp_88 = 0.37605877282253791*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29;
      real_t tmp_89 = tmp_31*tmp_88;
      real_t tmp_90 = 0.37605877282253791*tmp_34 + 0.039308471900058539*tmp_35 + tmp_36;
      real_t tmp_91 = tmp_38*tmp_90;
      real_t tmp_92 = tmp_87 + tmp_89 + tmp_91;
      real_t tmp_93 = tmp_41*tmp_86;
      real_t tmp_94 = tmp_43*tmp_88;
      real_t tmp_95 = tmp_45*tmp_90;
      real_t tmp_96 = tmp_93 + tmp_94 + tmp_95;
      real_t tmp_97 = tmp_48*tmp_86;
      real_t tmp_98 = tmp_50*tmp_88;
      real_t tmp_99 = tmp_52*tmp_90;
      real_t tmp_100 = tmp_97 + tmp_98 + tmp_99;
      real_t tmp_101 = tmp_1*(tmp_92 - 1.0/4.0) + tmp_14*(tmp_100 - 1.0/4.0) + tmp_18*(tmp_96 - 1.0/4.0);
      real_t tmp_102 = -tmp_87 - tmp_89 - tmp_91 - tmp_93 - tmp_94 - tmp_95 - tmp_97 - tmp_98 - tmp_99 + 1;
      real_t tmp_103 = tmp_101*tmp_63;
      real_t tmp_104 = 0.020848748529055869*tmp_65;
      real_t tmp_105 = 0.78764240869137092*tmp_3 + 0.1711304259088916*tmp_4 + tmp_5;
      real_t tmp_106 = tmp_105*tmp_24;
      real_t tmp_107 = 0.78764240869137092*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29;
      real_t tmp_108 = tmp_107*tmp_31;
      real_t tmp_109 = 0.78764240869137092*tmp_34 + 0.1711304259088916*tmp_35 + tmp_36;
      real_t tmp_110 = tmp_109*tmp_38;
      real_t tmp_111 = tmp_106 + tmp_108 + tmp_110;
      real_t tmp_112 = tmp_105*tmp_41;
      real_t tmp_113 = tmp_107*tmp_43;
      real_t tmp_114 = tmp_109*tmp_45;
      real_t tmp_115 = tmp_112 + tmp_113 + tmp_114;
      real_t tmp_116 = tmp_105*tmp_48;
      real_t tmp_117 = tmp_107*tmp_50;
      real_t tmp_118 = tmp_109*tmp_52;
      real_t tmp_119 = tmp_116 + tmp_117 + tmp_118;
      real_t tmp_120 = tmp_1*(tmp_111 - 1.0/4.0) + tmp_14*(tmp_119 - 1.0/4.0) + tmp_18*(tmp_115 - 1.0/4.0);
      real_t tmp_121 = -tmp_106 - tmp_108 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_116 - tmp_117 - tmp_118 + 1;
      real_t tmp_122 = tmp_120*tmp_63;
      real_t tmp_123 = 0.019202922745021479*tmp_65;
      real_t tmp_124 = 0.58463275527740355*tmp_3 + 0.37605877282253791*tmp_4 + tmp_5;
      real_t tmp_125 = tmp_124*tmp_24;
      real_t tmp_126 = 0.58463275527740355*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29;
      real_t tmp_127 = tmp_126*tmp_31;
      real_t tmp_128 = 0.58463275527740355*tmp_34 + 0.37605877282253791*tmp_35 + tmp_36;
      real_t tmp_129 = tmp_128*tmp_38;
      real_t tmp_130 = tmp_125 + tmp_127 + tmp_129;
      real_t tmp_131 = tmp_124*tmp_41;
      real_t tmp_132 = tmp_126*tmp_43;
      real_t tmp_133 = tmp_128*tmp_45;
      real_t tmp_134 = tmp_131 + tmp_132 + tmp_133;
      real_t tmp_135 = tmp_124*tmp_48;
      real_t tmp_136 = tmp_126*tmp_50;
      real_t tmp_137 = tmp_128*tmp_52;
      real_t tmp_138 = tmp_135 + tmp_136 + tmp_137;
      real_t tmp_139 = tmp_1*(tmp_130 - 1.0/4.0) + tmp_14*(tmp_138 - 1.0/4.0) + tmp_18*(tmp_134 - 1.0/4.0);
      real_t tmp_140 = -tmp_125 - tmp_127 - tmp_129 - tmp_131 - tmp_132 - tmp_133 - tmp_135 - tmp_136 - tmp_137 + 1;
      real_t tmp_141 = tmp_139*tmp_63;
      real_t tmp_142 = 0.020848748529055869*tmp_65;
      real_t tmp_143 = 0.041227165399737475*tmp_3 + 0.78764240869137092*tmp_4 + tmp_5;
      real_t tmp_144 = tmp_143*tmp_24;
      real_t tmp_145 = 0.041227165399737475*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29;
      real_t tmp_146 = tmp_145*tmp_31;
      real_t tmp_147 = 0.041227165399737475*tmp_34 + 0.78764240869137092*tmp_35 + tmp_36;
      real_t tmp_148 = tmp_147*tmp_38;
      real_t tmp_149 = tmp_144 + tmp_146 + tmp_148;
      real_t tmp_150 = tmp_143*tmp_41;
      real_t tmp_151 = tmp_145*tmp_43;
      real_t tmp_152 = tmp_147*tmp_45;
      real_t tmp_153 = tmp_150 + tmp_151 + tmp_152;
      real_t tmp_154 = tmp_143*tmp_48;
      real_t tmp_155 = tmp_145*tmp_50;
      real_t tmp_156 = tmp_147*tmp_52;
      real_t tmp_157 = tmp_154 + tmp_155 + tmp_156;
      real_t tmp_158 = tmp_1*(tmp_149 - 1.0/4.0) + tmp_14*(tmp_157 - 1.0/4.0) + tmp_18*(tmp_153 - 1.0/4.0);
      real_t tmp_159 = -tmp_144 - tmp_146 - tmp_148 - tmp_150 - tmp_151 - tmp_152 - tmp_154 - tmp_155 - tmp_156 + 1;
      real_t tmp_160 = tmp_158*tmp_63;
      real_t tmp_161 = 0.019202922745021479*tmp_65;
      real_t tmp_162 = 0.039308471900058539*tmp_3 + 0.58463275527740355*tmp_4 + tmp_5;
      real_t tmp_163 = tmp_162*tmp_24;
      real_t tmp_164 = 0.039308471900058539*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29;
      real_t tmp_165 = tmp_164*tmp_31;
      real_t tmp_166 = 0.039308471900058539*tmp_34 + 0.58463275527740355*tmp_35 + tmp_36;
      real_t tmp_167 = tmp_166*tmp_38;
      real_t tmp_168 = tmp_163 + tmp_165 + tmp_167;
      real_t tmp_169 = tmp_162*tmp_41;
      real_t tmp_170 = tmp_164*tmp_43;
      real_t tmp_171 = tmp_166*tmp_45;
      real_t tmp_172 = tmp_169 + tmp_170 + tmp_171;
      real_t tmp_173 = tmp_162*tmp_48;
      real_t tmp_174 = tmp_164*tmp_50;
      real_t tmp_175 = tmp_166*tmp_52;
      real_t tmp_176 = tmp_173 + tmp_174 + tmp_175;
      real_t tmp_177 = tmp_1*(tmp_168 - 1.0/4.0) + tmp_14*(tmp_176 - 1.0/4.0) + tmp_18*(tmp_172 - 1.0/4.0);
      real_t tmp_178 = -tmp_163 - tmp_165 - tmp_167 - tmp_169 - tmp_170 - tmp_171 - tmp_173 - tmp_174 - tmp_175 + 1;
      real_t tmp_179 = tmp_177*tmp_63;
      real_t tmp_180 = 0.020848748529055869*tmp_65;
      real_t tmp_181 = 0.78764240869137092*tmp_3 + 0.041227165399737475*tmp_4 + tmp_5;
      real_t tmp_182 = tmp_181*tmp_24;
      real_t tmp_183 = 0.78764240869137092*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29;
      real_t tmp_184 = tmp_183*tmp_31;
      real_t tmp_185 = 0.78764240869137092*tmp_34 + 0.041227165399737475*tmp_35 + tmp_36;
      real_t tmp_186 = tmp_185*tmp_38;
      real_t tmp_187 = tmp_182 + tmp_184 + tmp_186;
      real_t tmp_188 = tmp_181*tmp_41;
      real_t tmp_189 = tmp_183*tmp_43;
      real_t tmp_190 = tmp_185*tmp_45;
      real_t tmp_191 = tmp_188 + tmp_189 + tmp_190;
      real_t tmp_192 = tmp_181*tmp_48;
      real_t tmp_193 = tmp_183*tmp_50;
      real_t tmp_194 = tmp_185*tmp_52;
      real_t tmp_195 = tmp_192 + tmp_193 + tmp_194;
      real_t tmp_196 = tmp_1*(tmp_187 - 1.0/4.0) + tmp_14*(tmp_195 - 1.0/4.0) + tmp_18*(tmp_191 - 1.0/4.0);
      real_t tmp_197 = -tmp_182 - tmp_184 - tmp_186 - tmp_188 - tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 + 1;
      real_t tmp_198 = tmp_196*tmp_63;
      real_t tmp_199 = 0.019202922745021479*tmp_65;
      real_t tmp_200 = 0.58463275527740355*tmp_3 + 0.039308471900058539*tmp_4 + tmp_5;
      real_t tmp_201 = tmp_200*tmp_24;
      real_t tmp_202 = 0.58463275527740355*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29;
      real_t tmp_203 = tmp_202*tmp_31;
      real_t tmp_204 = 0.58463275527740355*tmp_34 + 0.039308471900058539*tmp_35 + tmp_36;
      real_t tmp_205 = tmp_204*tmp_38;
      real_t tmp_206 = tmp_201 + tmp_203 + tmp_205;
      real_t tmp_207 = tmp_200*tmp_41;
      real_t tmp_208 = tmp_202*tmp_43;
      real_t tmp_209 = tmp_204*tmp_45;
      real_t tmp_210 = tmp_207 + tmp_208 + tmp_209;
      real_t tmp_211 = tmp_200*tmp_48;
      real_t tmp_212 = tmp_202*tmp_50;
      real_t tmp_213 = tmp_204*tmp_52;
      real_t tmp_214 = tmp_211 + tmp_212 + tmp_213;
      real_t tmp_215 = tmp_1*(tmp_206 - 1.0/4.0) + tmp_14*(tmp_214 - 1.0/4.0) + tmp_18*(tmp_210 - 1.0/4.0);
      real_t tmp_216 = -tmp_201 - tmp_203 - tmp_205 - tmp_207 - tmp_208 - tmp_209 - tmp_211 - tmp_212 - tmp_213 + 1;
      real_t tmp_217 = tmp_215*tmp_63;
      real_t tmp_218 = 0.020848748529055869*tmp_65;
      real_t tmp_219 = 0.1711304259088916*tmp_3 + 0.78764240869137092*tmp_4 + tmp_5;
      real_t tmp_220 = tmp_219*tmp_24;
      real_t tmp_221 = 0.1711304259088916*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29;
      real_t tmp_222 = tmp_221*tmp_31;
      real_t tmp_223 = 0.1711304259088916*tmp_34 + 0.78764240869137092*tmp_35 + tmp_36;
      real_t tmp_224 = tmp_223*tmp_38;
      real_t tmp_225 = tmp_220 + tmp_222 + tmp_224;
      real_t tmp_226 = tmp_219*tmp_41;
      real_t tmp_227 = tmp_221*tmp_43;
      real_t tmp_228 = tmp_223*tmp_45;
      real_t tmp_229 = tmp_226 + tmp_227 + tmp_228;
      real_t tmp_230 = tmp_219*tmp_48;
      real_t tmp_231 = tmp_221*tmp_50;
      real_t tmp_232 = tmp_223*tmp_52;
      real_t tmp_233 = tmp_230 + tmp_231 + tmp_232;
      real_t tmp_234 = tmp_1*(tmp_225 - 1.0/4.0) + tmp_14*(tmp_233 - 1.0/4.0) + tmp_18*(tmp_229 - 1.0/4.0);
      real_t tmp_235 = -tmp_220 - tmp_222 - tmp_224 - tmp_226 - tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 + 1;
      real_t tmp_236 = tmp_234*tmp_63;
      real_t tmp_237 = 0.019202922745021479*tmp_65;
      real_t tmp_238 = 0.37605877282253791*tmp_3 + 0.58463275527740355*tmp_4 + tmp_5;
      real_t tmp_239 = tmp_238*tmp_24;
      real_t tmp_240 = 0.37605877282253791*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29;
      real_t tmp_241 = tmp_240*tmp_31;
      real_t tmp_242 = 0.37605877282253791*tmp_34 + 0.58463275527740355*tmp_35 + tmp_36;
      real_t tmp_243 = tmp_242*tmp_38;
      real_t tmp_244 = tmp_239 + tmp_241 + tmp_243;
      real_t tmp_245 = tmp_238*tmp_41;
      real_t tmp_246 = tmp_240*tmp_43;
      real_t tmp_247 = tmp_242*tmp_45;
      real_t tmp_248 = tmp_245 + tmp_246 + tmp_247;
      real_t tmp_249 = tmp_238*tmp_48;
      real_t tmp_250 = tmp_240*tmp_50;
      real_t tmp_251 = tmp_242*tmp_52;
      real_t tmp_252 = tmp_249 + tmp_250 + tmp_251;
      real_t tmp_253 = tmp_1*(tmp_244 - 1.0/4.0) + tmp_14*(tmp_252 - 1.0/4.0) + tmp_18*(tmp_248 - 1.0/4.0);
      real_t tmp_254 = -tmp_239 - tmp_241 - tmp_243 - tmp_245 - tmp_246 - tmp_247 - tmp_249 - tmp_250 - tmp_251 + 1;
      real_t tmp_255 = tmp_253*tmp_63;
      real_t tmp_256 = 0.020848748529055869*tmp_65;
      real_t tmp_257 = 0.041227165399737475*tmp_3 + 0.1711304259088916*tmp_4 + tmp_5;
      real_t tmp_258 = tmp_24*tmp_257;
      real_t tmp_259 = 0.041227165399737475*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29;
      real_t tmp_260 = tmp_259*tmp_31;
      real_t tmp_261 = 0.041227165399737475*tmp_34 + 0.1711304259088916*tmp_35 + tmp_36;
      real_t tmp_262 = tmp_261*tmp_38;
      real_t tmp_263 = tmp_258 + tmp_260 + tmp_262;
      real_t tmp_264 = tmp_257*tmp_41;
      real_t tmp_265 = tmp_259*tmp_43;
      real_t tmp_266 = tmp_261*tmp_45;
      real_t tmp_267 = tmp_264 + tmp_265 + tmp_266;
      real_t tmp_268 = tmp_257*tmp_48;
      real_t tmp_269 = tmp_259*tmp_50;
      real_t tmp_270 = tmp_261*tmp_52;
      real_t tmp_271 = tmp_268 + tmp_269 + tmp_270;
      real_t tmp_272 = tmp_1*(tmp_263 - 1.0/4.0) + tmp_14*(tmp_271 - 1.0/4.0) + tmp_18*(tmp_267 - 1.0/4.0);
      real_t tmp_273 = -tmp_258 - tmp_260 - tmp_262 - tmp_264 - tmp_265 - tmp_266 - tmp_268 - tmp_269 - tmp_270 + 1;
      real_t tmp_274 = tmp_272*tmp_63;
      real_t tmp_275 = 0.019202922745021479*tmp_65;
      real_t tmp_276 = 0.40446199974765351*tmp_3 + 0.19107600050469298*tmp_4 + tmp_5;
      real_t tmp_277 = tmp_24*tmp_276;
      real_t tmp_278 = 0.40446199974765351*tmp_27 + 0.19107600050469298*tmp_28 + tmp_29;
      real_t tmp_279 = tmp_278*tmp_31;
      real_t tmp_280 = 0.40446199974765351*tmp_34 + 0.19107600050469298*tmp_35 + tmp_36;
      real_t tmp_281 = tmp_280*tmp_38;
      real_t tmp_282 = tmp_277 + tmp_279 + tmp_281;
      real_t tmp_283 = tmp_276*tmp_41;
      real_t tmp_284 = tmp_278*tmp_43;
      real_t tmp_285 = tmp_280*tmp_45;
      real_t tmp_286 = tmp_283 + tmp_284 + tmp_285;
      real_t tmp_287 = tmp_276*tmp_48;
      real_t tmp_288 = tmp_278*tmp_50;
      real_t tmp_289 = tmp_280*tmp_52;
      real_t tmp_290 = tmp_287 + tmp_288 + tmp_289;
      real_t tmp_291 = tmp_1*(tmp_282 - 1.0/4.0) + tmp_14*(tmp_290 - 1.0/4.0) + tmp_18*(tmp_286 - 1.0/4.0);
      real_t tmp_292 = -tmp_277 - tmp_279 - tmp_281 - tmp_283 - tmp_284 - tmp_285 - tmp_287 - tmp_288 - tmp_289 + 1;
      real_t tmp_293 = tmp_291*tmp_63;
      real_t tmp_294 = 0.042507265838595799*tmp_65;
      real_t tmp_295 = 0.039308471900058539*tmp_3 + 0.37605877282253791*tmp_4 + tmp_5;
      real_t tmp_296 = tmp_24*tmp_295;
      real_t tmp_297 = 0.039308471900058539*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29;
      real_t tmp_298 = tmp_297*tmp_31;
      real_t tmp_299 = 0.039308471900058539*tmp_34 + 0.37605877282253791*tmp_35 + tmp_36;
      real_t tmp_300 = tmp_299*tmp_38;
      real_t tmp_301 = tmp_296 + tmp_298 + tmp_300;
      real_t tmp_302 = tmp_295*tmp_41;
      real_t tmp_303 = tmp_297*tmp_43;
      real_t tmp_304 = tmp_299*tmp_45;
      real_t tmp_305 = tmp_302 + tmp_303 + tmp_304;
      real_t tmp_306 = tmp_295*tmp_48;
      real_t tmp_307 = tmp_297*tmp_50;
      real_t tmp_308 = tmp_299*tmp_52;
      real_t tmp_309 = tmp_306 + tmp_307 + tmp_308;
      real_t tmp_310 = tmp_1*(tmp_301 - 1.0/4.0) + tmp_14*(tmp_309 - 1.0/4.0) + tmp_18*(tmp_305 - 1.0/4.0);
      real_t tmp_311 = -tmp_296 - tmp_298 - tmp_300 - tmp_302 - tmp_303 - tmp_304 - tmp_306 - tmp_307 - tmp_308 + 1;
      real_t tmp_312 = tmp_310*tmp_63;
      real_t tmp_313 = 0.020848748529055869*tmp_65;
      real_t tmp_314 = 0.93718850182767688*tmp_3 + 0.031405749086161582*tmp_4 + tmp_5;
      real_t tmp_315 = tmp_24*tmp_314;
      real_t tmp_316 = 0.93718850182767688*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29;
      real_t tmp_317 = tmp_31*tmp_316;
      real_t tmp_318 = 0.93718850182767688*tmp_34 + 0.031405749086161582*tmp_35 + tmp_36;
      real_t tmp_319 = tmp_318*tmp_38;
      real_t tmp_320 = tmp_315 + tmp_317 + tmp_319;
      real_t tmp_321 = tmp_314*tmp_41;
      real_t tmp_322 = tmp_316*tmp_43;
      real_t tmp_323 = tmp_318*tmp_45;
      real_t tmp_324 = tmp_321 + tmp_322 + tmp_323;
      real_t tmp_325 = tmp_314*tmp_48;
      real_t tmp_326 = tmp_316*tmp_50;
      real_t tmp_327 = tmp_318*tmp_52;
      real_t tmp_328 = tmp_325 + tmp_326 + tmp_327;
      real_t tmp_329 = tmp_1*(tmp_320 - 1.0/4.0) + tmp_14*(tmp_328 - 1.0/4.0) + tmp_18*(tmp_324 - 1.0/4.0);
      real_t tmp_330 = -tmp_315 - tmp_317 - tmp_319 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1;
      real_t tmp_331 = tmp_329*tmp_63;
      real_t tmp_332 = 0.0068572537431980923*tmp_65;
      real_t tmp_333 = 0.60796128279561268*tmp_3 + 0.19601935860219369*tmp_4 + tmp_5;
      real_t tmp_334 = tmp_24*tmp_333;
      real_t tmp_335 = 0.60796128279561268*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29;
      real_t tmp_336 = tmp_31*tmp_335;
      real_t tmp_337 = 0.60796128279561268*tmp_34 + 0.19601935860219369*tmp_35 + tmp_36;
      real_t tmp_338 = tmp_337*tmp_38;
      real_t tmp_339 = tmp_334 + tmp_336 + tmp_338;
      real_t tmp_340 = tmp_333*tmp_41;
      real_t tmp_341 = tmp_335*tmp_43;
      real_t tmp_342 = tmp_337*tmp_45;
      real_t tmp_343 = tmp_340 + tmp_341 + tmp_342;
      real_t tmp_344 = tmp_333*tmp_48;
      real_t tmp_345 = tmp_335*tmp_50;
      real_t tmp_346 = tmp_337*tmp_52;
      real_t tmp_347 = tmp_344 + tmp_345 + tmp_346;
      real_t tmp_348 = tmp_1*(tmp_339 - 1.0/4.0) + tmp_14*(tmp_347 - 1.0/4.0) + tmp_18*(tmp_343 - 1.0/4.0);
      real_t tmp_349 = -tmp_334 - tmp_336 - tmp_338 - tmp_340 - tmp_341 - tmp_342 - tmp_344 - tmp_345 - tmp_346 + 1;
      real_t tmp_350 = tmp_348*tmp_63;
      real_t tmp_351 = 0.037198804536718075*tmp_65;
      real_t tmp_352 = 0.19107600050469298*tmp_3 + 0.40446199974765351*tmp_4 + tmp_5;
      real_t tmp_353 = tmp_24*tmp_352;
      real_t tmp_354 = 0.19107600050469298*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29;
      real_t tmp_355 = tmp_31*tmp_354;
      real_t tmp_356 = 0.19107600050469298*tmp_34 + 0.40446199974765351*tmp_35 + tmp_36;
      real_t tmp_357 = tmp_356*tmp_38;
      real_t tmp_358 = tmp_353 + tmp_355 + tmp_357;
      real_t tmp_359 = tmp_352*tmp_41;
      real_t tmp_360 = tmp_354*tmp_43;
      real_t tmp_361 = tmp_356*tmp_45;
      real_t tmp_362 = tmp_359 + tmp_360 + tmp_361;
      real_t tmp_363 = tmp_352*tmp_48;
      real_t tmp_364 = tmp_354*tmp_50;
      real_t tmp_365 = tmp_356*tmp_52;
      real_t tmp_366 = tmp_363 + tmp_364 + tmp_365;
      real_t tmp_367 = tmp_1*(tmp_358 - 1.0/4.0) + tmp_14*(tmp_366 - 1.0/4.0) + tmp_18*(tmp_362 - 1.0/4.0);
      real_t tmp_368 = -tmp_353 - tmp_355 - tmp_357 - tmp_359 - tmp_360 - tmp_361 - tmp_363 - tmp_364 - tmp_365 + 1;
      real_t tmp_369 = tmp_367*tmp_63;
      real_t tmp_370 = 0.042507265838595799*tmp_65;
      real_t tmp_371 = 0.031405749086161582*tmp_3 + 0.031405749086161582*tmp_4 + tmp_5;
      real_t tmp_372 = tmp_24*tmp_371;
      real_t tmp_373 = 0.031405749086161582*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29;
      real_t tmp_374 = tmp_31*tmp_373;
      real_t tmp_375 = 0.031405749086161582*tmp_34 + 0.031405749086161582*tmp_35 + tmp_36;
      real_t tmp_376 = tmp_375*tmp_38;
      real_t tmp_377 = tmp_372 + tmp_374 + tmp_376;
      real_t tmp_378 = tmp_371*tmp_41;
      real_t tmp_379 = tmp_373*tmp_43;
      real_t tmp_380 = tmp_375*tmp_45;
      real_t tmp_381 = tmp_378 + tmp_379 + tmp_380;
      real_t tmp_382 = tmp_371*tmp_48;
      real_t tmp_383 = tmp_373*tmp_50;
      real_t tmp_384 = tmp_375*tmp_52;
      real_t tmp_385 = tmp_382 + tmp_383 + tmp_384;
      real_t tmp_386 = tmp_1*(tmp_377 - 1.0/4.0) + tmp_14*(tmp_385 - 1.0/4.0) + tmp_18*(tmp_381 - 1.0/4.0);
      real_t tmp_387 = -tmp_372 - tmp_374 - tmp_376 - tmp_378 - tmp_379 - tmp_380 - tmp_382 - tmp_383 - tmp_384 + 1;
      real_t tmp_388 = tmp_386*tmp_63;
      real_t tmp_389 = 0.0068572537431980923*tmp_65;
      real_t tmp_390 = 0.19601935860219369*tmp_3 + 0.19601935860219369*tmp_4 + tmp_5;
      real_t tmp_391 = tmp_24*tmp_390;
      real_t tmp_392 = 0.19601935860219369*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29;
      real_t tmp_393 = tmp_31*tmp_392;
      real_t tmp_394 = 0.19601935860219369*tmp_34 + 0.19601935860219369*tmp_35 + tmp_36;
      real_t tmp_395 = tmp_38*tmp_394;
      real_t tmp_396 = tmp_391 + tmp_393 + tmp_395;
      real_t tmp_397 = tmp_390*tmp_41;
      real_t tmp_398 = tmp_392*tmp_43;
      real_t tmp_399 = tmp_394*tmp_45;
      real_t tmp_400 = tmp_397 + tmp_398 + tmp_399;
      real_t tmp_401 = tmp_390*tmp_48;
      real_t tmp_402 = tmp_392*tmp_50;
      real_t tmp_403 = tmp_394*tmp_52;
      real_t tmp_404 = tmp_401 + tmp_402 + tmp_403;
      real_t tmp_405 = tmp_1*(tmp_396 - 1.0/4.0) + tmp_14*(tmp_404 - 1.0/4.0) + tmp_18*(tmp_400 - 1.0/4.0);
      real_t tmp_406 = -tmp_391 - tmp_393 - tmp_395 - tmp_397 - tmp_398 - tmp_399 - tmp_401 - tmp_402 - tmp_403 + 1;
      real_t tmp_407 = tmp_405*tmp_63;
      real_t tmp_408 = 0.037198804536718075*tmp_65;
      real_t tmp_409 = 0.40446199974765351*tmp_3 + 0.40446199974765351*tmp_4 + tmp_5;
      real_t tmp_410 = tmp_24*tmp_409;
      real_t tmp_411 = 0.40446199974765351*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29;
      real_t tmp_412 = tmp_31*tmp_411;
      real_t tmp_413 = 0.40446199974765351*tmp_34 + 0.40446199974765351*tmp_35 + tmp_36;
      real_t tmp_414 = tmp_38*tmp_413;
      real_t tmp_415 = tmp_410 + tmp_412 + tmp_414;
      real_t tmp_416 = tmp_409*tmp_41;
      real_t tmp_417 = tmp_411*tmp_43;
      real_t tmp_418 = tmp_413*tmp_45;
      real_t tmp_419 = tmp_416 + tmp_417 + tmp_418;
      real_t tmp_420 = tmp_409*tmp_48;
      real_t tmp_421 = tmp_411*tmp_50;
      real_t tmp_422 = tmp_413*tmp_52;
      real_t tmp_423 = tmp_420 + tmp_421 + tmp_422;
      real_t tmp_424 = tmp_1*(tmp_415 - 1.0/4.0) + tmp_14*(tmp_423 - 1.0/4.0) + tmp_18*(tmp_419 - 1.0/4.0);
      real_t tmp_425 = -tmp_410 - tmp_412 - tmp_414 - tmp_416 - tmp_417 - tmp_418 - tmp_420 - tmp_421 - tmp_422 + 1;
      real_t tmp_426 = tmp_424*tmp_63;
      real_t tmp_427 = 0.042507265838595799*tmp_65;
      real_t tmp_428 = 0.1711304259088916*tmp_3 + 0.041227165399737475*tmp_4 + tmp_5;
      real_t tmp_429 = tmp_24*tmp_428;
      real_t tmp_430 = 0.1711304259088916*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29;
      real_t tmp_431 = tmp_31*tmp_430;
      real_t tmp_432 = 0.1711304259088916*tmp_34 + 0.041227165399737475*tmp_35 + tmp_36;
      real_t tmp_433 = tmp_38*tmp_432;
      real_t tmp_434 = tmp_429 + tmp_431 + tmp_433;
      real_t tmp_435 = tmp_41*tmp_428;
      real_t tmp_436 = tmp_43*tmp_430;
      real_t tmp_437 = tmp_432*tmp_45;
      real_t tmp_438 = tmp_435 + tmp_436 + tmp_437;
      real_t tmp_439 = tmp_428*tmp_48;
      real_t tmp_440 = tmp_430*tmp_50;
      real_t tmp_441 = tmp_432*tmp_52;
      real_t tmp_442 = tmp_439 + tmp_440 + tmp_441;
      real_t tmp_443 = tmp_1*(tmp_434 - 1.0/4.0) + tmp_14*(tmp_442 - 1.0/4.0) + tmp_18*(tmp_438 - 1.0/4.0);
      real_t tmp_444 = -tmp_429 - tmp_431 - tmp_433 - tmp_435 - tmp_436 - tmp_437 - tmp_439 - tmp_440 - tmp_441 + 1;
      real_t tmp_445 = tmp_443*tmp_63;
      real_t tmp_446 = 0.019202922745021479*tmp_65;
      real_t tmp_447 = 0.5*p_affine_13_0*tmp_38 + 0.5*p_affine_13_1*tmp_31 + 0.5*p_affine_13_2*tmp_24;
      real_t tmp_448 = 0.5*p_affine_13_0*tmp_45 + 0.5*p_affine_13_1*tmp_43 + 0.5*p_affine_13_2*tmp_41;
      real_t tmp_449 = 0.5*p_affine_13_0*tmp_52 + 0.5*p_affine_13_1*tmp_50 + 0.5*p_affine_13_2*tmp_48;
      real_t a_0_0 = tmp_104*(-tmp_101*tmp_56 + tmp_102*tmp_103 - tmp_102*tmp_58) + tmp_123*(-tmp_120*tmp_56 + tmp_121*tmp_122 - tmp_121*tmp_58) + tmp_142*(-tmp_139*tmp_56 + tmp_140*tmp_141 - tmp_140*tmp_58) + tmp_161*(-tmp_158*tmp_56 + tmp_159*tmp_160 - tmp_159*tmp_58) + tmp_180*(-tmp_177*tmp_56 + tmp_178*tmp_179 - tmp_178*tmp_58) + tmp_199*(-tmp_196*tmp_56 + tmp_197*tmp_198 - tmp_197*tmp_58) + tmp_218*(-tmp_215*tmp_56 + tmp_216*tmp_217 - tmp_216*tmp_58) + tmp_237*(-tmp_234*tmp_56 + tmp_235*tmp_236 - tmp_235*tmp_58) + tmp_256*(-tmp_253*tmp_56 + tmp_254*tmp_255 - tmp_254*tmp_58) + tmp_275*(-tmp_272*tmp_56 + tmp_273*tmp_274 - tmp_273*tmp_58) + tmp_294*(-tmp_291*tmp_56 + tmp_292*tmp_293 - tmp_292*tmp_58) + tmp_313*(-tmp_310*tmp_56 + tmp_311*tmp_312 - tmp_311*tmp_58) + tmp_332*(-tmp_329*tmp_56 + tmp_330*tmp_331 - tmp_330*tmp_58) + tmp_351*(-tmp_348*tmp_56 + tmp_349*tmp_350 - tmp_349*tmp_58) + tmp_370*(-tmp_367*tmp_56 + tmp_368*tmp_369 - tmp_368*tmp_58) + tmp_389*(-tmp_386*tmp_56 + tmp_387*tmp_388 - tmp_387*tmp_58) + tmp_408*(-tmp_405*tmp_56 + tmp_406*tmp_407 - tmp_406*tmp_58) + tmp_427*(-tmp_424*tmp_56 + tmp_425*tmp_426 - tmp_425*tmp_58) + tmp_446*(-tmp_443*tmp_56 + tmp_444*tmp_445 - tmp_444*tmp_58) + tmp_66*(-tmp_55*tmp_56 - tmp_57*tmp_58 + tmp_57*tmp_64) + tmp_85*(-tmp_56*tmp_82 - tmp_58*tmp_83 + tmp_83*tmp_84);
      real_t a_1_0 = tmp_104*(-tmp_101*tmp_447 + tmp_103*tmp_92 - tmp_58*tmp_92) + tmp_123*(tmp_111*tmp_122 - tmp_111*tmp_58 - tmp_120*tmp_447) + tmp_142*(tmp_130*tmp_141 - tmp_130*tmp_58 - tmp_139*tmp_447) + tmp_161*(tmp_149*tmp_160 - tmp_149*tmp_58 - tmp_158*tmp_447) + tmp_180*(tmp_168*tmp_179 - tmp_168*tmp_58 - tmp_177*tmp_447) + tmp_199*(tmp_187*tmp_198 - tmp_187*tmp_58 - tmp_196*tmp_447) + tmp_218*(tmp_206*tmp_217 - tmp_206*tmp_58 - tmp_215*tmp_447) + tmp_237*(tmp_225*tmp_236 - tmp_225*tmp_58 - tmp_234*tmp_447) + tmp_256*(tmp_244*tmp_255 - tmp_244*tmp_58 - tmp_253*tmp_447) + tmp_275*(tmp_263*tmp_274 - tmp_263*tmp_58 - tmp_272*tmp_447) + tmp_294*(tmp_282*tmp_293 - tmp_282*tmp_58 - tmp_291*tmp_447) + tmp_313*(tmp_301*tmp_312 - tmp_301*tmp_58 - tmp_310*tmp_447) + tmp_332*(tmp_320*tmp_331 - tmp_320*tmp_58 - tmp_329*tmp_447) + tmp_351*(tmp_339*tmp_350 - tmp_339*tmp_58 - tmp_348*tmp_447) + tmp_370*(tmp_358*tmp_369 - tmp_358*tmp_58 - tmp_367*tmp_447) + tmp_389*(tmp_377*tmp_388 - tmp_377*tmp_58 - tmp_386*tmp_447) + tmp_408*(tmp_396*tmp_407 - tmp_396*tmp_58 - tmp_405*tmp_447) + tmp_427*(tmp_415*tmp_426 - tmp_415*tmp_58 - tmp_424*tmp_447) + tmp_446*(tmp_434*tmp_445 - tmp_434*tmp_58 - tmp_443*tmp_447) + tmp_66*(-tmp_40*tmp_58 + tmp_40*tmp_64 - tmp_447*tmp_55) + tmp_85*(-tmp_447*tmp_82 - tmp_58*tmp_73 + tmp_73*tmp_84);
      real_t a_2_0 = tmp_104*(-tmp_101*tmp_448 + tmp_103*tmp_96 - tmp_58*tmp_96) + tmp_123*(tmp_115*tmp_122 - tmp_115*tmp_58 - tmp_120*tmp_448) + tmp_142*(tmp_134*tmp_141 - tmp_134*tmp_58 - tmp_139*tmp_448) + tmp_161*(tmp_153*tmp_160 - tmp_153*tmp_58 - tmp_158*tmp_448) + tmp_180*(tmp_172*tmp_179 - tmp_172*tmp_58 - tmp_177*tmp_448) + tmp_199*(tmp_191*tmp_198 - tmp_191*tmp_58 - tmp_196*tmp_448) + tmp_218*(tmp_210*tmp_217 - tmp_210*tmp_58 - tmp_215*tmp_448) + tmp_237*(tmp_229*tmp_236 - tmp_229*tmp_58 - tmp_234*tmp_448) + tmp_256*(tmp_248*tmp_255 - tmp_248*tmp_58 - tmp_253*tmp_448) + tmp_275*(tmp_267*tmp_274 - tmp_267*tmp_58 - tmp_272*tmp_448) + tmp_294*(tmp_286*tmp_293 - tmp_286*tmp_58 - tmp_291*tmp_448) + tmp_313*(tmp_305*tmp_312 - tmp_305*tmp_58 - tmp_310*tmp_448) + tmp_332*(tmp_324*tmp_331 - tmp_324*tmp_58 - tmp_329*tmp_448) + tmp_351*(tmp_343*tmp_350 - tmp_343*tmp_58 - tmp_348*tmp_448) + tmp_370*(tmp_362*tmp_369 - tmp_362*tmp_58 - tmp_367*tmp_448) + tmp_389*(tmp_381*tmp_388 - tmp_381*tmp_58 - tmp_386*tmp_448) + tmp_408*(tmp_400*tmp_407 - tmp_400*tmp_58 - tmp_405*tmp_448) + tmp_427*(tmp_419*tmp_426 - tmp_419*tmp_58 - tmp_424*tmp_448) + tmp_446*(tmp_438*tmp_445 - tmp_438*tmp_58 - tmp_443*tmp_448) + tmp_66*(-tmp_448*tmp_55 - tmp_47*tmp_58 + tmp_47*tmp_64) + tmp_85*(-tmp_448*tmp_82 - tmp_58*tmp_77 + tmp_77*tmp_84);
      real_t a_3_0 = tmp_104*(tmp_100*tmp_103 - tmp_100*tmp_58 - tmp_101*tmp_449) + tmp_123*(tmp_119*tmp_122 - tmp_119*tmp_58 - tmp_120*tmp_449) + tmp_142*(tmp_138*tmp_141 - tmp_138*tmp_58 - tmp_139*tmp_449) + tmp_161*(tmp_157*tmp_160 - tmp_157*tmp_58 - tmp_158*tmp_449) + tmp_180*(tmp_176*tmp_179 - tmp_176*tmp_58 - tmp_177*tmp_449) + tmp_199*(tmp_195*tmp_198 - tmp_195*tmp_58 - tmp_196*tmp_449) + tmp_218*(tmp_214*tmp_217 - tmp_214*tmp_58 - tmp_215*tmp_449) + tmp_237*(tmp_233*tmp_236 - tmp_233*tmp_58 - tmp_234*tmp_449) + tmp_256*(tmp_252*tmp_255 - tmp_252*tmp_58 - tmp_253*tmp_449) + tmp_275*(tmp_271*tmp_274 - tmp_271*tmp_58 - tmp_272*tmp_449) + tmp_294*(tmp_290*tmp_293 - tmp_290*tmp_58 - tmp_291*tmp_449) + tmp_313*(tmp_309*tmp_312 - tmp_309*tmp_58 - tmp_310*tmp_449) + tmp_332*(tmp_328*tmp_331 - tmp_328*tmp_58 - tmp_329*tmp_449) + tmp_351*(tmp_347*tmp_350 - tmp_347*tmp_58 - tmp_348*tmp_449) + tmp_370*(tmp_366*tmp_369 - tmp_366*tmp_58 - tmp_367*tmp_449) + tmp_389*(tmp_385*tmp_388 - tmp_385*tmp_58 - tmp_386*tmp_449) + tmp_408*(tmp_404*tmp_407 - tmp_404*tmp_58 - tmp_405*tmp_449) + tmp_427*(tmp_423*tmp_426 - tmp_423*tmp_58 - tmp_424*tmp_449) + tmp_446*(tmp_442*tmp_445 - tmp_442*tmp_58 - tmp_443*tmp_449) + tmp_66*(-tmp_449*tmp_55 - tmp_54*tmp_58 + tmp_54*tmp_64) + tmp_85*(-tmp_449*tmp_82 - tmp_58*tmp_81 + tmp_81*tmp_84);
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


      real_t tmp_0 = -p_affine_4_2;
      real_t tmp_1 = p_affine_5_2 + tmp_0;
      real_t tmp_2 = -p_affine_4_0;
      real_t tmp_3 = p_affine_6_0 + tmp_2;
      real_t tmp_4 = -p_affine_4_1;
      real_t tmp_5 = p_affine_7_1 + tmp_4;
      real_t tmp_6 = tmp_3*tmp_5;
      real_t tmp_7 = p_affine_7_0 + tmp_2;
      real_t tmp_8 = p_affine_6_1 + tmp_4;
      real_t tmp_9 = tmp_7*tmp_8;
      real_t tmp_10 = tmp_6 - tmp_9;
      real_t tmp_11 = p_affine_5_0 + tmp_2;
      real_t tmp_12 = p_affine_7_2 + tmp_0;
      real_t tmp_13 = tmp_12*tmp_8;
      real_t tmp_14 = p_affine_5_1 + tmp_4;
      real_t tmp_15 = p_affine_6_2 + tmp_0;
      real_t tmp_16 = tmp_15*tmp_7;
      real_t tmp_17 = tmp_15*tmp_5;
      real_t tmp_18 = tmp_12*tmp_3;
      real_t tmp_19 = 1.0 / (tmp_1*tmp_6 - tmp_1*tmp_9 + tmp_11*tmp_13 - tmp_11*tmp_17 + tmp_14*tmp_16 - tmp_14*tmp_18);
      real_t tmp_20 = p_affine_8_2 + tmp_0;
      real_t tmp_21 = -p_affine_8_2;
      real_t tmp_22 = p_affine_9_2 + tmp_21;
      real_t tmp_23 = p_affine_10_2 + tmp_21;
      real_t tmp_24 = 0.031405749086161582*tmp_22 + 0.93718850182767688*tmp_23;
      real_t tmp_25 = tmp_19*(tmp_20 + tmp_24);
      real_t tmp_26 = tmp_16 - tmp_18;
      real_t tmp_27 = p_affine_8_1 + tmp_4;
      real_t tmp_28 = -p_affine_8_1;
      real_t tmp_29 = p_affine_9_1 + tmp_28;
      real_t tmp_30 = p_affine_10_1 + tmp_28;
      real_t tmp_31 = 0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_30;
      real_t tmp_32 = tmp_19*(tmp_27 + tmp_31);
      real_t tmp_33 = tmp_13 - tmp_17;
      real_t tmp_34 = p_affine_8_0 + tmp_2;
      real_t tmp_35 = -p_affine_8_0;
      real_t tmp_36 = p_affine_9_0 + tmp_35;
      real_t tmp_37 = p_affine_10_0 + tmp_35;
      real_t tmp_38 = 0.031405749086161582*tmp_36 + 0.93718850182767688*tmp_37;
      real_t tmp_39 = tmp_19*(tmp_34 + tmp_38);
      real_t tmp_40 = -tmp_11*tmp_5 + tmp_14*tmp_7;
      real_t tmp_41 = -tmp_1*tmp_7 + tmp_11*tmp_12;
      real_t tmp_42 = tmp_1*tmp_5 - tmp_12*tmp_14;
      real_t tmp_43 = tmp_11*tmp_8 - tmp_14*tmp_3;
      real_t tmp_44 = tmp_1*tmp_3 - tmp_11*tmp_15;
      real_t tmp_45 = -tmp_1*tmp_8 + tmp_14*tmp_15;
      real_t tmp_46 = tmp_1*(tmp_10*tmp_25 + tmp_26*tmp_32 + tmp_33*tmp_39 - 1.0/4.0) + tmp_12*(tmp_25*tmp_43 + tmp_32*tmp_44 + tmp_39*tmp_45 - 1.0/4.0) + tmp_15*(tmp_25*tmp_40 + tmp_32*tmp_41 + tmp_39*tmp_42 - 1.0/4.0);
      real_t tmp_47 = -p_affine_0_1;
      real_t tmp_48 = p_affine_1_1 + tmp_47;
      real_t tmp_49 = -p_affine_0_2;
      real_t tmp_50 = p_affine_2_2 + tmp_49;
      real_t tmp_51 = tmp_48*tmp_50;
      real_t tmp_52 = p_affine_2_1 + tmp_47;
      real_t tmp_53 = p_affine_1_2 + tmp_49;
      real_t tmp_54 = tmp_52*tmp_53;
      real_t tmp_55 = -p_affine_0_0;
      real_t tmp_56 = p_affine_1_0 + tmp_55;
      real_t tmp_57 = p_affine_3_2 + tmp_49;
      real_t tmp_58 = tmp_52*tmp_57;
      real_t tmp_59 = p_affine_2_0 + tmp_55;
      real_t tmp_60 = p_affine_3_1 + tmp_47;
      real_t tmp_61 = tmp_53*tmp_60;
      real_t tmp_62 = p_affine_3_0 + tmp_55;
      real_t tmp_63 = tmp_50*tmp_60;
      real_t tmp_64 = tmp_48*tmp_57;
      real_t tmp_65 = 1.0 / (tmp_51*tmp_62 - tmp_54*tmp_62 + tmp_56*tmp_58 - tmp_56*tmp_63 + tmp_59*tmp_61 - tmp_59*tmp_64);
      real_t tmp_66 = tmp_65*(tmp_51 - tmp_54);
      real_t tmp_67 = tmp_65*(tmp_61 - tmp_64);
      real_t tmp_68 = tmp_65*(tmp_58 - tmp_63);
      real_t tmp_69 = tmp_65*(-tmp_50*tmp_56 + tmp_53*tmp_59);
      real_t tmp_70 = tmp_65*(-tmp_53*tmp_62 + tmp_56*tmp_57);
      real_t tmp_71 = tmp_65*(tmp_50*tmp_62 - tmp_57*tmp_59);
      real_t tmp_72 = tmp_65*(-tmp_48*tmp_59 + tmp_52*tmp_56);
      real_t tmp_73 = tmp_65*(tmp_48*tmp_62 - tmp_56*tmp_60);
      real_t tmp_74 = tmp_65*(-tmp_52*tmp_62 + tmp_59*tmp_60);
      real_t tmp_75 = 0.5*p_affine_13_0*(-tmp_66 - tmp_67 - tmp_68) + 0.5*p_affine_13_1*(-tmp_69 - tmp_70 - tmp_71) + 0.5*p_affine_13_2*(-tmp_72 - tmp_73 - tmp_74);
      real_t tmp_76 = p_affine_8_2 + tmp_49;
      real_t tmp_77 = tmp_24 + tmp_76;
      real_t tmp_78 = tmp_72*tmp_77;
      real_t tmp_79 = tmp_73*tmp_77;
      real_t tmp_80 = p_affine_8_1 + tmp_47;
      real_t tmp_81 = tmp_31 + tmp_80;
      real_t tmp_82 = tmp_69*tmp_81;
      real_t tmp_83 = tmp_70*tmp_81;
      real_t tmp_84 = tmp_74*tmp_77;
      real_t tmp_85 = tmp_71*tmp_81;
      real_t tmp_86 = p_affine_8_0 + tmp_55;
      real_t tmp_87 = tmp_38 + tmp_86;
      real_t tmp_88 = tmp_66*tmp_87;
      real_t tmp_89 = tmp_67*tmp_87;
      real_t tmp_90 = tmp_68*tmp_87;
      real_t tmp_91 = -tmp_78 - tmp_79 - tmp_82 - tmp_83 - tmp_84 - tmp_85 - tmp_88 - tmp_89 - tmp_90 + 1;
      real_t tmp_92 = tmp_1*tmp_19;
      real_t tmp_93 = tmp_15*tmp_19;
      real_t tmp_94 = tmp_12*tmp_19;
      real_t tmp_95 = 0.5*p_affine_13_0*(tmp_33*tmp_92 + tmp_42*tmp_93 + tmp_45*tmp_94) + 0.5*p_affine_13_1*(tmp_26*tmp_92 + tmp_41*tmp_93 + tmp_44*tmp_94) + 0.5*p_affine_13_2*(tmp_10*tmp_92 + tmp_40*tmp_93 + tmp_43*tmp_94);
      real_t tmp_96 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_97 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_98 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_99 = (std::abs(tmp_23*tmp_96 - tmp_30*tmp_98)*std::abs(tmp_23*tmp_96 - tmp_30*tmp_98)) + (std::abs(tmp_23*tmp_97 - tmp_37*tmp_98)*std::abs(tmp_23*tmp_97 - tmp_37*tmp_98)) + (std::abs(tmp_30*tmp_97 - tmp_37*tmp_96)*std::abs(tmp_30*tmp_97 - tmp_37*tmp_96));
      real_t tmp_100 = 5.0*std::pow(tmp_99, -0.25);
      real_t tmp_101 = tmp_100*tmp_46;
      real_t tmp_102 = 1.0*std::pow(tmp_99, 1.0/2.0);
      real_t tmp_103 = 0.0068572537431980923*tmp_102;
      real_t tmp_104 = 0.19601935860219369*tmp_22 + 0.60796128279561268*tmp_23;
      real_t tmp_105 = tmp_19*(tmp_104 + tmp_20);
      real_t tmp_106 = 0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_30;
      real_t tmp_107 = tmp_19*(tmp_106 + tmp_27);
      real_t tmp_108 = 0.19601935860219369*tmp_36 + 0.60796128279561268*tmp_37;
      real_t tmp_109 = tmp_19*(tmp_108 + tmp_34);
      real_t tmp_110 = tmp_1*(tmp_10*tmp_105 + tmp_107*tmp_26 + tmp_109*tmp_33 - 1.0/4.0) + tmp_12*(tmp_105*tmp_43 + tmp_107*tmp_44 + tmp_109*tmp_45 - 1.0/4.0) + tmp_15*(tmp_105*tmp_40 + tmp_107*tmp_41 + tmp_109*tmp_42 - 1.0/4.0);
      real_t tmp_111 = tmp_104 + tmp_76;
      real_t tmp_112 = tmp_111*tmp_72;
      real_t tmp_113 = tmp_111*tmp_73;
      real_t tmp_114 = tmp_106 + tmp_80;
      real_t tmp_115 = tmp_114*tmp_69;
      real_t tmp_116 = tmp_114*tmp_70;
      real_t tmp_117 = tmp_111*tmp_74;
      real_t tmp_118 = tmp_114*tmp_71;
      real_t tmp_119 = tmp_108 + tmp_86;
      real_t tmp_120 = tmp_119*tmp_66;
      real_t tmp_121 = tmp_119*tmp_67;
      real_t tmp_122 = tmp_119*tmp_68;
      real_t tmp_123 = -tmp_112 - tmp_113 - tmp_115 - tmp_116 - tmp_117 - tmp_118 - tmp_120 - tmp_121 - tmp_122 + 1;
      real_t tmp_124 = tmp_100*tmp_110;
      real_t tmp_125 = 0.037198804536718075*tmp_102;
      real_t tmp_126 = 0.37605877282253791*tmp_22 + 0.039308471900058539*tmp_23;
      real_t tmp_127 = tmp_19*(tmp_126 + tmp_20);
      real_t tmp_128 = 0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_30;
      real_t tmp_129 = tmp_19*(tmp_128 + tmp_27);
      real_t tmp_130 = 0.37605877282253791*tmp_36 + 0.039308471900058539*tmp_37;
      real_t tmp_131 = tmp_19*(tmp_130 + tmp_34);
      real_t tmp_132 = tmp_1*(tmp_10*tmp_127 + tmp_129*tmp_26 + tmp_131*tmp_33 - 1.0/4.0) + tmp_12*(tmp_127*tmp_43 + tmp_129*tmp_44 + tmp_131*tmp_45 - 1.0/4.0) + tmp_15*(tmp_127*tmp_40 + tmp_129*tmp_41 + tmp_131*tmp_42 - 1.0/4.0);
      real_t tmp_133 = tmp_126 + tmp_76;
      real_t tmp_134 = tmp_133*tmp_72;
      real_t tmp_135 = tmp_133*tmp_73;
      real_t tmp_136 = tmp_128 + tmp_80;
      real_t tmp_137 = tmp_136*tmp_69;
      real_t tmp_138 = tmp_136*tmp_70;
      real_t tmp_139 = tmp_133*tmp_74;
      real_t tmp_140 = tmp_136*tmp_71;
      real_t tmp_141 = tmp_130 + tmp_86;
      real_t tmp_142 = tmp_141*tmp_66;
      real_t tmp_143 = tmp_141*tmp_67;
      real_t tmp_144 = tmp_141*tmp_68;
      real_t tmp_145 = -tmp_134 - tmp_135 - tmp_137 - tmp_138 - tmp_139 - tmp_140 - tmp_142 - tmp_143 - tmp_144 + 1;
      real_t tmp_146 = tmp_100*tmp_132;
      real_t tmp_147 = 0.020848748529055869*tmp_102;
      real_t tmp_148 = 0.78764240869137092*tmp_22 + 0.1711304259088916*tmp_23;
      real_t tmp_149 = tmp_19*(tmp_148 + tmp_20);
      real_t tmp_150 = 0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_30;
      real_t tmp_151 = tmp_19*(tmp_150 + tmp_27);
      real_t tmp_152 = 0.78764240869137092*tmp_36 + 0.1711304259088916*tmp_37;
      real_t tmp_153 = tmp_19*(tmp_152 + tmp_34);
      real_t tmp_154 = tmp_1*(tmp_10*tmp_149 + tmp_151*tmp_26 + tmp_153*tmp_33 - 1.0/4.0) + tmp_12*(tmp_149*tmp_43 + tmp_151*tmp_44 + tmp_153*tmp_45 - 1.0/4.0) + tmp_15*(tmp_149*tmp_40 + tmp_151*tmp_41 + tmp_153*tmp_42 - 1.0/4.0);
      real_t tmp_155 = tmp_148 + tmp_76;
      real_t tmp_156 = tmp_155*tmp_72;
      real_t tmp_157 = tmp_155*tmp_73;
      real_t tmp_158 = tmp_150 + tmp_80;
      real_t tmp_159 = tmp_158*tmp_69;
      real_t tmp_160 = tmp_158*tmp_70;
      real_t tmp_161 = tmp_155*tmp_74;
      real_t tmp_162 = tmp_158*tmp_71;
      real_t tmp_163 = tmp_152 + tmp_86;
      real_t tmp_164 = tmp_163*tmp_66;
      real_t tmp_165 = tmp_163*tmp_67;
      real_t tmp_166 = tmp_163*tmp_68;
      real_t tmp_167 = -tmp_156 - tmp_157 - tmp_159 - tmp_160 - tmp_161 - tmp_162 - tmp_164 - tmp_165 - tmp_166 + 1;
      real_t tmp_168 = tmp_100*tmp_154;
      real_t tmp_169 = 0.019202922745021479*tmp_102;
      real_t tmp_170 = 0.58463275527740355*tmp_22 + 0.37605877282253791*tmp_23;
      real_t tmp_171 = tmp_19*(tmp_170 + tmp_20);
      real_t tmp_172 = 0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_30;
      real_t tmp_173 = tmp_19*(tmp_172 + tmp_27);
      real_t tmp_174 = 0.58463275527740355*tmp_36 + 0.37605877282253791*tmp_37;
      real_t tmp_175 = tmp_19*(tmp_174 + tmp_34);
      real_t tmp_176 = tmp_1*(tmp_10*tmp_171 + tmp_173*tmp_26 + tmp_175*tmp_33 - 1.0/4.0) + tmp_12*(tmp_171*tmp_43 + tmp_173*tmp_44 + tmp_175*tmp_45 - 1.0/4.0) + tmp_15*(tmp_171*tmp_40 + tmp_173*tmp_41 + tmp_175*tmp_42 - 1.0/4.0);
      real_t tmp_177 = tmp_170 + tmp_76;
      real_t tmp_178 = tmp_177*tmp_72;
      real_t tmp_179 = tmp_177*tmp_73;
      real_t tmp_180 = tmp_172 + tmp_80;
      real_t tmp_181 = tmp_180*tmp_69;
      real_t tmp_182 = tmp_180*tmp_70;
      real_t tmp_183 = tmp_177*tmp_74;
      real_t tmp_184 = tmp_180*tmp_71;
      real_t tmp_185 = tmp_174 + tmp_86;
      real_t tmp_186 = tmp_185*tmp_66;
      real_t tmp_187 = tmp_185*tmp_67;
      real_t tmp_188 = tmp_185*tmp_68;
      real_t tmp_189 = -tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 - tmp_184 - tmp_186 - tmp_187 - tmp_188 + 1;
      real_t tmp_190 = tmp_100*tmp_176;
      real_t tmp_191 = 0.020848748529055869*tmp_102;
      real_t tmp_192 = 0.041227165399737475*tmp_22 + 0.78764240869137092*tmp_23;
      real_t tmp_193 = tmp_19*(tmp_192 + tmp_20);
      real_t tmp_194 = 0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_30;
      real_t tmp_195 = tmp_19*(tmp_194 + tmp_27);
      real_t tmp_196 = 0.041227165399737475*tmp_36 + 0.78764240869137092*tmp_37;
      real_t tmp_197 = tmp_19*(tmp_196 + tmp_34);
      real_t tmp_198 = tmp_1*(tmp_10*tmp_193 + tmp_195*tmp_26 + tmp_197*tmp_33 - 1.0/4.0) + tmp_12*(tmp_193*tmp_43 + tmp_195*tmp_44 + tmp_197*tmp_45 - 1.0/4.0) + tmp_15*(tmp_193*tmp_40 + tmp_195*tmp_41 + tmp_197*tmp_42 - 1.0/4.0);
      real_t tmp_199 = tmp_192 + tmp_76;
      real_t tmp_200 = tmp_199*tmp_72;
      real_t tmp_201 = tmp_199*tmp_73;
      real_t tmp_202 = tmp_194 + tmp_80;
      real_t tmp_203 = tmp_202*tmp_69;
      real_t tmp_204 = tmp_202*tmp_70;
      real_t tmp_205 = tmp_199*tmp_74;
      real_t tmp_206 = tmp_202*tmp_71;
      real_t tmp_207 = tmp_196 + tmp_86;
      real_t tmp_208 = tmp_207*tmp_66;
      real_t tmp_209 = tmp_207*tmp_67;
      real_t tmp_210 = tmp_207*tmp_68;
      real_t tmp_211 = -tmp_200 - tmp_201 - tmp_203 - tmp_204 - tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 + 1;
      real_t tmp_212 = tmp_100*tmp_198;
      real_t tmp_213 = 0.019202922745021479*tmp_102;
      real_t tmp_214 = 0.039308471900058539*tmp_22 + 0.58463275527740355*tmp_23;
      real_t tmp_215 = tmp_19*(tmp_20 + tmp_214);
      real_t tmp_216 = 0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_30;
      real_t tmp_217 = tmp_19*(tmp_216 + tmp_27);
      real_t tmp_218 = 0.039308471900058539*tmp_36 + 0.58463275527740355*tmp_37;
      real_t tmp_219 = tmp_19*(tmp_218 + tmp_34);
      real_t tmp_220 = tmp_1*(tmp_10*tmp_215 + tmp_217*tmp_26 + tmp_219*tmp_33 - 1.0/4.0) + tmp_12*(tmp_215*tmp_43 + tmp_217*tmp_44 + tmp_219*tmp_45 - 1.0/4.0) + tmp_15*(tmp_215*tmp_40 + tmp_217*tmp_41 + tmp_219*tmp_42 - 1.0/4.0);
      real_t tmp_221 = tmp_214 + tmp_76;
      real_t tmp_222 = tmp_221*tmp_72;
      real_t tmp_223 = tmp_221*tmp_73;
      real_t tmp_224 = tmp_216 + tmp_80;
      real_t tmp_225 = tmp_224*tmp_69;
      real_t tmp_226 = tmp_224*tmp_70;
      real_t tmp_227 = tmp_221*tmp_74;
      real_t tmp_228 = tmp_224*tmp_71;
      real_t tmp_229 = tmp_218 + tmp_86;
      real_t tmp_230 = tmp_229*tmp_66;
      real_t tmp_231 = tmp_229*tmp_67;
      real_t tmp_232 = tmp_229*tmp_68;
      real_t tmp_233 = -tmp_222 - tmp_223 - tmp_225 - tmp_226 - tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 + 1;
      real_t tmp_234 = tmp_100*tmp_220;
      real_t tmp_235 = 0.020848748529055869*tmp_102;
      real_t tmp_236 = 0.78764240869137092*tmp_22 + 0.041227165399737475*tmp_23;
      real_t tmp_237 = tmp_19*(tmp_20 + tmp_236);
      real_t tmp_238 = 0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_30;
      real_t tmp_239 = tmp_19*(tmp_238 + tmp_27);
      real_t tmp_240 = 0.78764240869137092*tmp_36 + 0.041227165399737475*tmp_37;
      real_t tmp_241 = tmp_19*(tmp_240 + tmp_34);
      real_t tmp_242 = tmp_1*(tmp_10*tmp_237 + tmp_239*tmp_26 + tmp_241*tmp_33 - 1.0/4.0) + tmp_12*(tmp_237*tmp_43 + tmp_239*tmp_44 + tmp_241*tmp_45 - 1.0/4.0) + tmp_15*(tmp_237*tmp_40 + tmp_239*tmp_41 + tmp_241*tmp_42 - 1.0/4.0);
      real_t tmp_243 = tmp_236 + tmp_76;
      real_t tmp_244 = tmp_243*tmp_72;
      real_t tmp_245 = tmp_243*tmp_73;
      real_t tmp_246 = tmp_238 + tmp_80;
      real_t tmp_247 = tmp_246*tmp_69;
      real_t tmp_248 = tmp_246*tmp_70;
      real_t tmp_249 = tmp_243*tmp_74;
      real_t tmp_250 = tmp_246*tmp_71;
      real_t tmp_251 = tmp_240 + tmp_86;
      real_t tmp_252 = tmp_251*tmp_66;
      real_t tmp_253 = tmp_251*tmp_67;
      real_t tmp_254 = tmp_251*tmp_68;
      real_t tmp_255 = -tmp_244 - tmp_245 - tmp_247 - tmp_248 - tmp_249 - tmp_250 - tmp_252 - tmp_253 - tmp_254 + 1;
      real_t tmp_256 = tmp_100*tmp_242;
      real_t tmp_257 = 0.019202922745021479*tmp_102;
      real_t tmp_258 = 0.58463275527740355*tmp_22 + 0.039308471900058539*tmp_23;
      real_t tmp_259 = tmp_19*(tmp_20 + tmp_258);
      real_t tmp_260 = 0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_30;
      real_t tmp_261 = tmp_19*(tmp_260 + tmp_27);
      real_t tmp_262 = 0.58463275527740355*tmp_36 + 0.039308471900058539*tmp_37;
      real_t tmp_263 = tmp_19*(tmp_262 + tmp_34);
      real_t tmp_264 = tmp_1*(tmp_10*tmp_259 + tmp_26*tmp_261 + tmp_263*tmp_33 - 1.0/4.0) + tmp_12*(tmp_259*tmp_43 + tmp_261*tmp_44 + tmp_263*tmp_45 - 1.0/4.0) + tmp_15*(tmp_259*tmp_40 + tmp_261*tmp_41 + tmp_263*tmp_42 - 1.0/4.0);
      real_t tmp_265 = tmp_258 + tmp_76;
      real_t tmp_266 = tmp_265*tmp_72;
      real_t tmp_267 = tmp_265*tmp_73;
      real_t tmp_268 = tmp_260 + tmp_80;
      real_t tmp_269 = tmp_268*tmp_69;
      real_t tmp_270 = tmp_268*tmp_70;
      real_t tmp_271 = tmp_265*tmp_74;
      real_t tmp_272 = tmp_268*tmp_71;
      real_t tmp_273 = tmp_262 + tmp_86;
      real_t tmp_274 = tmp_273*tmp_66;
      real_t tmp_275 = tmp_273*tmp_67;
      real_t tmp_276 = tmp_273*tmp_68;
      real_t tmp_277 = -tmp_266 - tmp_267 - tmp_269 - tmp_270 - tmp_271 - tmp_272 - tmp_274 - tmp_275 - tmp_276 + 1;
      real_t tmp_278 = tmp_100*tmp_264;
      real_t tmp_279 = 0.020848748529055869*tmp_102;
      real_t tmp_280 = 0.1711304259088916*tmp_22 + 0.78764240869137092*tmp_23;
      real_t tmp_281 = tmp_19*(tmp_20 + tmp_280);
      real_t tmp_282 = 0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_30;
      real_t tmp_283 = tmp_19*(tmp_27 + tmp_282);
      real_t tmp_284 = 0.1711304259088916*tmp_36 + 0.78764240869137092*tmp_37;
      real_t tmp_285 = tmp_19*(tmp_284 + tmp_34);
      real_t tmp_286 = tmp_1*(tmp_10*tmp_281 + tmp_26*tmp_283 + tmp_285*tmp_33 - 1.0/4.0) + tmp_12*(tmp_281*tmp_43 + tmp_283*tmp_44 + tmp_285*tmp_45 - 1.0/4.0) + tmp_15*(tmp_281*tmp_40 + tmp_283*tmp_41 + tmp_285*tmp_42 - 1.0/4.0);
      real_t tmp_287 = tmp_280 + tmp_76;
      real_t tmp_288 = tmp_287*tmp_72;
      real_t tmp_289 = tmp_287*tmp_73;
      real_t tmp_290 = tmp_282 + tmp_80;
      real_t tmp_291 = tmp_290*tmp_69;
      real_t tmp_292 = tmp_290*tmp_70;
      real_t tmp_293 = tmp_287*tmp_74;
      real_t tmp_294 = tmp_290*tmp_71;
      real_t tmp_295 = tmp_284 + tmp_86;
      real_t tmp_296 = tmp_295*tmp_66;
      real_t tmp_297 = tmp_295*tmp_67;
      real_t tmp_298 = tmp_295*tmp_68;
      real_t tmp_299 = -tmp_288 - tmp_289 - tmp_291 - tmp_292 - tmp_293 - tmp_294 - tmp_296 - tmp_297 - tmp_298 + 1;
      real_t tmp_300 = tmp_100*tmp_286;
      real_t tmp_301 = 0.019202922745021479*tmp_102;
      real_t tmp_302 = 0.37605877282253791*tmp_22 + 0.58463275527740355*tmp_23;
      real_t tmp_303 = tmp_19*(tmp_20 + tmp_302);
      real_t tmp_304 = 0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_30;
      real_t tmp_305 = tmp_19*(tmp_27 + tmp_304);
      real_t tmp_306 = 0.37605877282253791*tmp_36 + 0.58463275527740355*tmp_37;
      real_t tmp_307 = tmp_19*(tmp_306 + tmp_34);
      real_t tmp_308 = tmp_1*(tmp_10*tmp_303 + tmp_26*tmp_305 + tmp_307*tmp_33 - 1.0/4.0) + tmp_12*(tmp_303*tmp_43 + tmp_305*tmp_44 + tmp_307*tmp_45 - 1.0/4.0) + tmp_15*(tmp_303*tmp_40 + tmp_305*tmp_41 + tmp_307*tmp_42 - 1.0/4.0);
      real_t tmp_309 = tmp_302 + tmp_76;
      real_t tmp_310 = tmp_309*tmp_72;
      real_t tmp_311 = tmp_309*tmp_73;
      real_t tmp_312 = tmp_304 + tmp_80;
      real_t tmp_313 = tmp_312*tmp_69;
      real_t tmp_314 = tmp_312*tmp_70;
      real_t tmp_315 = tmp_309*tmp_74;
      real_t tmp_316 = tmp_312*tmp_71;
      real_t tmp_317 = tmp_306 + tmp_86;
      real_t tmp_318 = tmp_317*tmp_66;
      real_t tmp_319 = tmp_317*tmp_67;
      real_t tmp_320 = tmp_317*tmp_68;
      real_t tmp_321 = -tmp_310 - tmp_311 - tmp_313 - tmp_314 - tmp_315 - tmp_316 - tmp_318 - tmp_319 - tmp_320 + 1;
      real_t tmp_322 = tmp_100*tmp_308;
      real_t tmp_323 = 0.020848748529055869*tmp_102;
      real_t tmp_324 = 0.041227165399737475*tmp_22 + 0.1711304259088916*tmp_23;
      real_t tmp_325 = tmp_19*(tmp_20 + tmp_324);
      real_t tmp_326 = 0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_30;
      real_t tmp_327 = tmp_19*(tmp_27 + tmp_326);
      real_t tmp_328 = 0.041227165399737475*tmp_36 + 0.1711304259088916*tmp_37;
      real_t tmp_329 = tmp_19*(tmp_328 + tmp_34);
      real_t tmp_330 = tmp_1*(tmp_10*tmp_325 + tmp_26*tmp_327 + tmp_329*tmp_33 - 1.0/4.0) + tmp_12*(tmp_325*tmp_43 + tmp_327*tmp_44 + tmp_329*tmp_45 - 1.0/4.0) + tmp_15*(tmp_325*tmp_40 + tmp_327*tmp_41 + tmp_329*tmp_42 - 1.0/4.0);
      real_t tmp_331 = tmp_324 + tmp_76;
      real_t tmp_332 = tmp_331*tmp_72;
      real_t tmp_333 = tmp_331*tmp_73;
      real_t tmp_334 = tmp_326 + tmp_80;
      real_t tmp_335 = tmp_334*tmp_69;
      real_t tmp_336 = tmp_334*tmp_70;
      real_t tmp_337 = tmp_331*tmp_74;
      real_t tmp_338 = tmp_334*tmp_71;
      real_t tmp_339 = tmp_328 + tmp_86;
      real_t tmp_340 = tmp_339*tmp_66;
      real_t tmp_341 = tmp_339*tmp_67;
      real_t tmp_342 = tmp_339*tmp_68;
      real_t tmp_343 = -tmp_332 - tmp_333 - tmp_335 - tmp_336 - tmp_337 - tmp_338 - tmp_340 - tmp_341 - tmp_342 + 1;
      real_t tmp_344 = tmp_100*tmp_330;
      real_t tmp_345 = 0.019202922745021479*tmp_102;
      real_t tmp_346 = 0.40446199974765351*tmp_22 + 0.19107600050469298*tmp_23;
      real_t tmp_347 = tmp_19*(tmp_20 + tmp_346);
      real_t tmp_348 = 0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_30;
      real_t tmp_349 = tmp_19*(tmp_27 + tmp_348);
      real_t tmp_350 = 0.40446199974765351*tmp_36 + 0.19107600050469298*tmp_37;
      real_t tmp_351 = tmp_19*(tmp_34 + tmp_350);
      real_t tmp_352 = tmp_1*(tmp_10*tmp_347 + tmp_26*tmp_349 + tmp_33*tmp_351 - 1.0/4.0) + tmp_12*(tmp_347*tmp_43 + tmp_349*tmp_44 + tmp_351*tmp_45 - 1.0/4.0) + tmp_15*(tmp_347*tmp_40 + tmp_349*tmp_41 + tmp_351*tmp_42 - 1.0/4.0);
      real_t tmp_353 = tmp_346 + tmp_76;
      real_t tmp_354 = tmp_353*tmp_72;
      real_t tmp_355 = tmp_353*tmp_73;
      real_t tmp_356 = tmp_348 + tmp_80;
      real_t tmp_357 = tmp_356*tmp_69;
      real_t tmp_358 = tmp_356*tmp_70;
      real_t tmp_359 = tmp_353*tmp_74;
      real_t tmp_360 = tmp_356*tmp_71;
      real_t tmp_361 = tmp_350 + tmp_86;
      real_t tmp_362 = tmp_361*tmp_66;
      real_t tmp_363 = tmp_361*tmp_67;
      real_t tmp_364 = tmp_361*tmp_68;
      real_t tmp_365 = -tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 - tmp_360 - tmp_362 - tmp_363 - tmp_364 + 1;
      real_t tmp_366 = tmp_100*tmp_352;
      real_t tmp_367 = 0.042507265838595799*tmp_102;
      real_t tmp_368 = 0.039308471900058539*tmp_22 + 0.37605877282253791*tmp_23;
      real_t tmp_369 = tmp_19*(tmp_20 + tmp_368);
      real_t tmp_370 = 0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_30;
      real_t tmp_371 = tmp_19*(tmp_27 + tmp_370);
      real_t tmp_372 = 0.039308471900058539*tmp_36 + 0.37605877282253791*tmp_37;
      real_t tmp_373 = tmp_19*(tmp_34 + tmp_372);
      real_t tmp_374 = tmp_1*(tmp_10*tmp_369 + tmp_26*tmp_371 + tmp_33*tmp_373 - 1.0/4.0) + tmp_12*(tmp_369*tmp_43 + tmp_371*tmp_44 + tmp_373*tmp_45 - 1.0/4.0) + tmp_15*(tmp_369*tmp_40 + tmp_371*tmp_41 + tmp_373*tmp_42 - 1.0/4.0);
      real_t tmp_375 = tmp_368 + tmp_76;
      real_t tmp_376 = tmp_375*tmp_72;
      real_t tmp_377 = tmp_375*tmp_73;
      real_t tmp_378 = tmp_370 + tmp_80;
      real_t tmp_379 = tmp_378*tmp_69;
      real_t tmp_380 = tmp_378*tmp_70;
      real_t tmp_381 = tmp_375*tmp_74;
      real_t tmp_382 = tmp_378*tmp_71;
      real_t tmp_383 = tmp_372 + tmp_86;
      real_t tmp_384 = tmp_383*tmp_66;
      real_t tmp_385 = tmp_383*tmp_67;
      real_t tmp_386 = tmp_383*tmp_68;
      real_t tmp_387 = -tmp_376 - tmp_377 - tmp_379 - tmp_380 - tmp_381 - tmp_382 - tmp_384 - tmp_385 - tmp_386 + 1;
      real_t tmp_388 = tmp_100*tmp_374;
      real_t tmp_389 = 0.020848748529055869*tmp_102;
      real_t tmp_390 = 0.93718850182767688*tmp_22 + 0.031405749086161582*tmp_23;
      real_t tmp_391 = tmp_19*(tmp_20 + tmp_390);
      real_t tmp_392 = 0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_30;
      real_t tmp_393 = tmp_19*(tmp_27 + tmp_392);
      real_t tmp_394 = 0.93718850182767688*tmp_36 + 0.031405749086161582*tmp_37;
      real_t tmp_395 = tmp_19*(tmp_34 + tmp_394);
      real_t tmp_396 = tmp_1*(tmp_10*tmp_391 + tmp_26*tmp_393 + tmp_33*tmp_395 - 1.0/4.0) + tmp_12*(tmp_391*tmp_43 + tmp_393*tmp_44 + tmp_395*tmp_45 - 1.0/4.0) + tmp_15*(tmp_391*tmp_40 + tmp_393*tmp_41 + tmp_395*tmp_42 - 1.0/4.0);
      real_t tmp_397 = tmp_390 + tmp_76;
      real_t tmp_398 = tmp_397*tmp_72;
      real_t tmp_399 = tmp_397*tmp_73;
      real_t tmp_400 = tmp_392 + tmp_80;
      real_t tmp_401 = tmp_400*tmp_69;
      real_t tmp_402 = tmp_400*tmp_70;
      real_t tmp_403 = tmp_397*tmp_74;
      real_t tmp_404 = tmp_400*tmp_71;
      real_t tmp_405 = tmp_394 + tmp_86;
      real_t tmp_406 = tmp_405*tmp_66;
      real_t tmp_407 = tmp_405*tmp_67;
      real_t tmp_408 = tmp_405*tmp_68;
      real_t tmp_409 = -tmp_398 - tmp_399 - tmp_401 - tmp_402 - tmp_403 - tmp_404 - tmp_406 - tmp_407 - tmp_408 + 1;
      real_t tmp_410 = tmp_100*tmp_396;
      real_t tmp_411 = 0.0068572537431980923*tmp_102;
      real_t tmp_412 = 0.60796128279561268*tmp_22 + 0.19601935860219369*tmp_23;
      real_t tmp_413 = tmp_19*(tmp_20 + tmp_412);
      real_t tmp_414 = 0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_30;
      real_t tmp_415 = tmp_19*(tmp_27 + tmp_414);
      real_t tmp_416 = 0.60796128279561268*tmp_36 + 0.19601935860219369*tmp_37;
      real_t tmp_417 = tmp_19*(tmp_34 + tmp_416);
      real_t tmp_418 = tmp_1*(tmp_10*tmp_413 + tmp_26*tmp_415 + tmp_33*tmp_417 - 1.0/4.0) + tmp_12*(tmp_413*tmp_43 + tmp_415*tmp_44 + tmp_417*tmp_45 - 1.0/4.0) + tmp_15*(tmp_40*tmp_413 + tmp_41*tmp_415 + tmp_417*tmp_42 - 1.0/4.0);
      real_t tmp_419 = tmp_412 + tmp_76;
      real_t tmp_420 = tmp_419*tmp_72;
      real_t tmp_421 = tmp_419*tmp_73;
      real_t tmp_422 = tmp_414 + tmp_80;
      real_t tmp_423 = tmp_422*tmp_69;
      real_t tmp_424 = tmp_422*tmp_70;
      real_t tmp_425 = tmp_419*tmp_74;
      real_t tmp_426 = tmp_422*tmp_71;
      real_t tmp_427 = tmp_416 + tmp_86;
      real_t tmp_428 = tmp_427*tmp_66;
      real_t tmp_429 = tmp_427*tmp_67;
      real_t tmp_430 = tmp_427*tmp_68;
      real_t tmp_431 = -tmp_420 - tmp_421 - tmp_423 - tmp_424 - tmp_425 - tmp_426 - tmp_428 - tmp_429 - tmp_430 + 1;
      real_t tmp_432 = tmp_100*tmp_418;
      real_t tmp_433 = 0.037198804536718075*tmp_102;
      real_t tmp_434 = 0.19107600050469298*tmp_22 + 0.40446199974765351*tmp_23;
      real_t tmp_435 = tmp_19*(tmp_20 + tmp_434);
      real_t tmp_436 = 0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_30;
      real_t tmp_437 = tmp_19*(tmp_27 + tmp_436);
      real_t tmp_438 = 0.19107600050469298*tmp_36 + 0.40446199974765351*tmp_37;
      real_t tmp_439 = tmp_19*(tmp_34 + tmp_438);
      real_t tmp_440 = tmp_1*(tmp_10*tmp_435 + tmp_26*tmp_437 + tmp_33*tmp_439 - 1.0/4.0) + tmp_12*(tmp_43*tmp_435 + tmp_437*tmp_44 + tmp_439*tmp_45 - 1.0/4.0) + tmp_15*(tmp_40*tmp_435 + tmp_41*tmp_437 + tmp_42*tmp_439 - 1.0/4.0);
      real_t tmp_441 = tmp_434 + tmp_76;
      real_t tmp_442 = tmp_441*tmp_72;
      real_t tmp_443 = tmp_441*tmp_73;
      real_t tmp_444 = tmp_436 + tmp_80;
      real_t tmp_445 = tmp_444*tmp_69;
      real_t tmp_446 = tmp_444*tmp_70;
      real_t tmp_447 = tmp_441*tmp_74;
      real_t tmp_448 = tmp_444*tmp_71;
      real_t tmp_449 = tmp_438 + tmp_86;
      real_t tmp_450 = tmp_449*tmp_66;
      real_t tmp_451 = tmp_449*tmp_67;
      real_t tmp_452 = tmp_449*tmp_68;
      real_t tmp_453 = -tmp_442 - tmp_443 - tmp_445 - tmp_446 - tmp_447 - tmp_448 - tmp_450 - tmp_451 - tmp_452 + 1;
      real_t tmp_454 = tmp_100*tmp_440;
      real_t tmp_455 = 0.042507265838595799*tmp_102;
      real_t tmp_456 = 0.031405749086161582*tmp_22 + 0.031405749086161582*tmp_23;
      real_t tmp_457 = tmp_19*(tmp_20 + tmp_456);
      real_t tmp_458 = 0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_30;
      real_t tmp_459 = tmp_19*(tmp_27 + tmp_458);
      real_t tmp_460 = 0.031405749086161582*tmp_36 + 0.031405749086161582*tmp_37;
      real_t tmp_461 = tmp_19*(tmp_34 + tmp_460);
      real_t tmp_462 = tmp_1*(tmp_10*tmp_457 + tmp_26*tmp_459 + tmp_33*tmp_461 - 1.0/4.0) + tmp_12*(tmp_43*tmp_457 + tmp_44*tmp_459 + tmp_45*tmp_461 - 1.0/4.0) + tmp_15*(tmp_40*tmp_457 + tmp_41*tmp_459 + tmp_42*tmp_461 - 1.0/4.0);
      real_t tmp_463 = tmp_456 + tmp_76;
      real_t tmp_464 = tmp_463*tmp_72;
      real_t tmp_465 = tmp_463*tmp_73;
      real_t tmp_466 = tmp_458 + tmp_80;
      real_t tmp_467 = tmp_466*tmp_69;
      real_t tmp_468 = tmp_466*tmp_70;
      real_t tmp_469 = tmp_463*tmp_74;
      real_t tmp_470 = tmp_466*tmp_71;
      real_t tmp_471 = tmp_460 + tmp_86;
      real_t tmp_472 = tmp_471*tmp_66;
      real_t tmp_473 = tmp_471*tmp_67;
      real_t tmp_474 = tmp_471*tmp_68;
      real_t tmp_475 = -tmp_464 - tmp_465 - tmp_467 - tmp_468 - tmp_469 - tmp_470 - tmp_472 - tmp_473 - tmp_474 + 1;
      real_t tmp_476 = tmp_100*tmp_462;
      real_t tmp_477 = 0.0068572537431980923*tmp_102;
      real_t tmp_478 = 0.19601935860219369*tmp_22 + 0.19601935860219369*tmp_23;
      real_t tmp_479 = tmp_19*(tmp_20 + tmp_478);
      real_t tmp_480 = 0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_30;
      real_t tmp_481 = tmp_19*(tmp_27 + tmp_480);
      real_t tmp_482 = 0.19601935860219369*tmp_36 + 0.19601935860219369*tmp_37;
      real_t tmp_483 = tmp_19*(tmp_34 + tmp_482);
      real_t tmp_484 = tmp_1*(tmp_10*tmp_479 + tmp_26*tmp_481 + tmp_33*tmp_483 - 1.0/4.0) + tmp_12*(tmp_43*tmp_479 + tmp_44*tmp_481 + tmp_45*tmp_483 - 1.0/4.0) + tmp_15*(tmp_40*tmp_479 + tmp_41*tmp_481 + tmp_42*tmp_483 - 1.0/4.0);
      real_t tmp_485 = tmp_478 + tmp_76;
      real_t tmp_486 = tmp_485*tmp_72;
      real_t tmp_487 = tmp_485*tmp_73;
      real_t tmp_488 = tmp_480 + tmp_80;
      real_t tmp_489 = tmp_488*tmp_69;
      real_t tmp_490 = tmp_488*tmp_70;
      real_t tmp_491 = tmp_485*tmp_74;
      real_t tmp_492 = tmp_488*tmp_71;
      real_t tmp_493 = tmp_482 + tmp_86;
      real_t tmp_494 = tmp_493*tmp_66;
      real_t tmp_495 = tmp_493*tmp_67;
      real_t tmp_496 = tmp_493*tmp_68;
      real_t tmp_497 = -tmp_486 - tmp_487 - tmp_489 - tmp_490 - tmp_491 - tmp_492 - tmp_494 - tmp_495 - tmp_496 + 1;
      real_t tmp_498 = tmp_100*tmp_484;
      real_t tmp_499 = 0.037198804536718075*tmp_102;
      real_t tmp_500 = 0.40446199974765351*tmp_22 + 0.40446199974765351*tmp_23;
      real_t tmp_501 = tmp_19*(tmp_20 + tmp_500);
      real_t tmp_502 = 0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_30;
      real_t tmp_503 = tmp_19*(tmp_27 + tmp_502);
      real_t tmp_504 = 0.40446199974765351*tmp_36 + 0.40446199974765351*tmp_37;
      real_t tmp_505 = tmp_19*(tmp_34 + tmp_504);
      real_t tmp_506 = tmp_1*(tmp_10*tmp_501 + tmp_26*tmp_503 + tmp_33*tmp_505 - 1.0/4.0) + tmp_12*(tmp_43*tmp_501 + tmp_44*tmp_503 + tmp_45*tmp_505 - 1.0/4.0) + tmp_15*(tmp_40*tmp_501 + tmp_41*tmp_503 + tmp_42*tmp_505 - 1.0/4.0);
      real_t tmp_507 = tmp_500 + tmp_76;
      real_t tmp_508 = tmp_507*tmp_72;
      real_t tmp_509 = tmp_507*tmp_73;
      real_t tmp_510 = tmp_502 + tmp_80;
      real_t tmp_511 = tmp_510*tmp_69;
      real_t tmp_512 = tmp_510*tmp_70;
      real_t tmp_513 = tmp_507*tmp_74;
      real_t tmp_514 = tmp_510*tmp_71;
      real_t tmp_515 = tmp_504 + tmp_86;
      real_t tmp_516 = tmp_515*tmp_66;
      real_t tmp_517 = tmp_515*tmp_67;
      real_t tmp_518 = tmp_515*tmp_68;
      real_t tmp_519 = -tmp_508 - tmp_509 - tmp_511 - tmp_512 - tmp_513 - tmp_514 - tmp_516 - tmp_517 - tmp_518 + 1;
      real_t tmp_520 = tmp_100*tmp_506;
      real_t tmp_521 = 0.042507265838595799*tmp_102;
      real_t tmp_522 = 0.1711304259088916*tmp_22 + 0.041227165399737475*tmp_23;
      real_t tmp_523 = tmp_19*(tmp_20 + tmp_522);
      real_t tmp_524 = 0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_30;
      real_t tmp_525 = tmp_19*(tmp_27 + tmp_524);
      real_t tmp_526 = 0.1711304259088916*tmp_36 + 0.041227165399737475*tmp_37;
      real_t tmp_527 = tmp_19*(tmp_34 + tmp_526);
      real_t tmp_528 = tmp_1*(tmp_10*tmp_523 + tmp_26*tmp_525 + tmp_33*tmp_527 - 1.0/4.0) + tmp_12*(tmp_43*tmp_523 + tmp_44*tmp_525 + tmp_45*tmp_527 - 1.0/4.0) + tmp_15*(tmp_40*tmp_523 + tmp_41*tmp_525 + tmp_42*tmp_527 - 1.0/4.0);
      real_t tmp_529 = tmp_522 + tmp_76;
      real_t tmp_530 = tmp_529*tmp_72;
      real_t tmp_531 = tmp_529*tmp_73;
      real_t tmp_532 = tmp_524 + tmp_80;
      real_t tmp_533 = tmp_532*tmp_69;
      real_t tmp_534 = tmp_532*tmp_70;
      real_t tmp_535 = tmp_529*tmp_74;
      real_t tmp_536 = tmp_532*tmp_71;
      real_t tmp_537 = tmp_526 + tmp_86;
      real_t tmp_538 = tmp_537*tmp_66;
      real_t tmp_539 = tmp_537*tmp_67;
      real_t tmp_540 = tmp_537*tmp_68;
      real_t tmp_541 = -tmp_530 - tmp_531 - tmp_533 - tmp_534 - tmp_535 - tmp_536 - tmp_538 - tmp_539 - tmp_540 + 1;
      real_t tmp_542 = tmp_100*tmp_528;
      real_t tmp_543 = 0.019202922745021479*tmp_102;
      real_t tmp_544 = tmp_84 + tmp_85 + tmp_90;
      real_t tmp_545 = 0.5*p_affine_13_0*tmp_68 + 0.5*p_affine_13_1*tmp_71 + 0.5*p_affine_13_2*tmp_74;
      real_t tmp_546 = tmp_117 + tmp_118 + tmp_122;
      real_t tmp_547 = tmp_139 + tmp_140 + tmp_144;
      real_t tmp_548 = tmp_161 + tmp_162 + tmp_166;
      real_t tmp_549 = tmp_183 + tmp_184 + tmp_188;
      real_t tmp_550 = tmp_205 + tmp_206 + tmp_210;
      real_t tmp_551 = tmp_227 + tmp_228 + tmp_232;
      real_t tmp_552 = tmp_249 + tmp_250 + tmp_254;
      real_t tmp_553 = tmp_271 + tmp_272 + tmp_276;
      real_t tmp_554 = tmp_293 + tmp_294 + tmp_298;
      real_t tmp_555 = tmp_315 + tmp_316 + tmp_320;
      real_t tmp_556 = tmp_337 + tmp_338 + tmp_342;
      real_t tmp_557 = tmp_359 + tmp_360 + tmp_364;
      real_t tmp_558 = tmp_381 + tmp_382 + tmp_386;
      real_t tmp_559 = tmp_403 + tmp_404 + tmp_408;
      real_t tmp_560 = tmp_425 + tmp_426 + tmp_430;
      real_t tmp_561 = tmp_447 + tmp_448 + tmp_452;
      real_t tmp_562 = tmp_469 + tmp_470 + tmp_474;
      real_t tmp_563 = tmp_491 + tmp_492 + tmp_496;
      real_t tmp_564 = tmp_513 + tmp_514 + tmp_518;
      real_t tmp_565 = tmp_535 + tmp_536 + tmp_540;
      real_t tmp_566 = tmp_79 + tmp_83 + tmp_89;
      real_t tmp_567 = 0.5*p_affine_13_0*tmp_67 + 0.5*p_affine_13_1*tmp_70 + 0.5*p_affine_13_2*tmp_73;
      real_t tmp_568 = tmp_113 + tmp_116 + tmp_121;
      real_t tmp_569 = tmp_135 + tmp_138 + tmp_143;
      real_t tmp_570 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_571 = tmp_179 + tmp_182 + tmp_187;
      real_t tmp_572 = tmp_201 + tmp_204 + tmp_209;
      real_t tmp_573 = tmp_223 + tmp_226 + tmp_231;
      real_t tmp_574 = tmp_245 + tmp_248 + tmp_253;
      real_t tmp_575 = tmp_267 + tmp_270 + tmp_275;
      real_t tmp_576 = tmp_289 + tmp_292 + tmp_297;
      real_t tmp_577 = tmp_311 + tmp_314 + tmp_319;
      real_t tmp_578 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_579 = tmp_355 + tmp_358 + tmp_363;
      real_t tmp_580 = tmp_377 + tmp_380 + tmp_385;
      real_t tmp_581 = tmp_399 + tmp_402 + tmp_407;
      real_t tmp_582 = tmp_421 + tmp_424 + tmp_429;
      real_t tmp_583 = tmp_443 + tmp_446 + tmp_451;
      real_t tmp_584 = tmp_465 + tmp_468 + tmp_473;
      real_t tmp_585 = tmp_487 + tmp_490 + tmp_495;
      real_t tmp_586 = tmp_509 + tmp_512 + tmp_517;
      real_t tmp_587 = tmp_531 + tmp_534 + tmp_539;
      real_t tmp_588 = tmp_78 + tmp_82 + tmp_88;
      real_t tmp_589 = 0.5*p_affine_13_0*tmp_66 + 0.5*p_affine_13_1*tmp_69 + 0.5*p_affine_13_2*tmp_72;
      real_t tmp_590 = tmp_112 + tmp_115 + tmp_120;
      real_t tmp_591 = tmp_134 + tmp_137 + tmp_142;
      real_t tmp_592 = tmp_156 + tmp_159 + tmp_164;
      real_t tmp_593 = tmp_178 + tmp_181 + tmp_186;
      real_t tmp_594 = tmp_200 + tmp_203 + tmp_208;
      real_t tmp_595 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_596 = tmp_244 + tmp_247 + tmp_252;
      real_t tmp_597 = tmp_266 + tmp_269 + tmp_274;
      real_t tmp_598 = tmp_288 + tmp_291 + tmp_296;
      real_t tmp_599 = tmp_310 + tmp_313 + tmp_318;
      real_t tmp_600 = tmp_332 + tmp_335 + tmp_340;
      real_t tmp_601 = tmp_354 + tmp_357 + tmp_362;
      real_t tmp_602 = tmp_376 + tmp_379 + tmp_384;
      real_t tmp_603 = tmp_398 + tmp_401 + tmp_406;
      real_t tmp_604 = tmp_420 + tmp_423 + tmp_428;
      real_t tmp_605 = tmp_442 + tmp_445 + tmp_450;
      real_t tmp_606 = tmp_464 + tmp_467 + tmp_472;
      real_t tmp_607 = tmp_486 + tmp_489 + tmp_494;
      real_t tmp_608 = tmp_508 + tmp_511 + tmp_516;
      real_t tmp_609 = tmp_530 + tmp_533 + tmp_538;
      real_t a_0_0 = tmp_103*(-tmp_101*tmp_91 + tmp_46*tmp_75 - tmp_91*tmp_95) + tmp_125*(tmp_110*tmp_75 - tmp_123*tmp_124 - tmp_123*tmp_95) + tmp_147*(tmp_132*tmp_75 - tmp_145*tmp_146 - tmp_145*tmp_95) + tmp_169*(tmp_154*tmp_75 - tmp_167*tmp_168 - tmp_167*tmp_95) + tmp_191*(tmp_176*tmp_75 - tmp_189*tmp_190 - tmp_189*tmp_95) + tmp_213*(tmp_198*tmp_75 - tmp_211*tmp_212 - tmp_211*tmp_95) + tmp_235*(tmp_220*tmp_75 - tmp_233*tmp_234 - tmp_233*tmp_95) + tmp_257*(tmp_242*tmp_75 - tmp_255*tmp_256 - tmp_255*tmp_95) + tmp_279*(tmp_264*tmp_75 - tmp_277*tmp_278 - tmp_277*tmp_95) + tmp_301*(tmp_286*tmp_75 - tmp_299*tmp_300 - tmp_299*tmp_95) + tmp_323*(tmp_308*tmp_75 - tmp_321*tmp_322 - tmp_321*tmp_95) + tmp_345*(tmp_330*tmp_75 - tmp_343*tmp_344 - tmp_343*tmp_95) + tmp_367*(tmp_352*tmp_75 - tmp_365*tmp_366 - tmp_365*tmp_95) + tmp_389*(tmp_374*tmp_75 - tmp_387*tmp_388 - tmp_387*tmp_95) + tmp_411*(tmp_396*tmp_75 - tmp_409*tmp_410 - tmp_409*tmp_95) + tmp_433*(tmp_418*tmp_75 - tmp_431*tmp_432 - tmp_431*tmp_95) + tmp_455*(tmp_440*tmp_75 - tmp_453*tmp_454 - tmp_453*tmp_95) + tmp_477*(tmp_462*tmp_75 - tmp_475*tmp_476 - tmp_475*tmp_95) + tmp_499*(tmp_484*tmp_75 - tmp_497*tmp_498 - tmp_497*tmp_95) + tmp_521*(tmp_506*tmp_75 - tmp_519*tmp_520 - tmp_519*tmp_95) + tmp_543*(tmp_528*tmp_75 - tmp_541*tmp_542 - tmp_541*tmp_95);
      real_t a_1_0 = tmp_103*(-tmp_101*tmp_544 + tmp_46*tmp_545 - tmp_544*tmp_95) + tmp_125*(tmp_110*tmp_545 - tmp_124*tmp_546 - tmp_546*tmp_95) + tmp_147*(tmp_132*tmp_545 - tmp_146*tmp_547 - tmp_547*tmp_95) + tmp_169*(tmp_154*tmp_545 - tmp_168*tmp_548 - tmp_548*tmp_95) + tmp_191*(tmp_176*tmp_545 - tmp_190*tmp_549 - tmp_549*tmp_95) + tmp_213*(tmp_198*tmp_545 - tmp_212*tmp_550 - tmp_550*tmp_95) + tmp_235*(tmp_220*tmp_545 - tmp_234*tmp_551 - tmp_551*tmp_95) + tmp_257*(tmp_242*tmp_545 - tmp_256*tmp_552 - tmp_552*tmp_95) + tmp_279*(tmp_264*tmp_545 - tmp_278*tmp_553 - tmp_553*tmp_95) + tmp_301*(tmp_286*tmp_545 - tmp_300*tmp_554 - tmp_554*tmp_95) + tmp_323*(tmp_308*tmp_545 - tmp_322*tmp_555 - tmp_555*tmp_95) + tmp_345*(tmp_330*tmp_545 - tmp_344*tmp_556 - tmp_556*tmp_95) + tmp_367*(tmp_352*tmp_545 - tmp_366*tmp_557 - tmp_557*tmp_95) + tmp_389*(tmp_374*tmp_545 - tmp_388*tmp_558 - tmp_558*tmp_95) + tmp_411*(tmp_396*tmp_545 - tmp_410*tmp_559 - tmp_559*tmp_95) + tmp_433*(tmp_418*tmp_545 - tmp_432*tmp_560 - tmp_560*tmp_95) + tmp_455*(tmp_440*tmp_545 - tmp_454*tmp_561 - tmp_561*tmp_95) + tmp_477*(tmp_462*tmp_545 - tmp_476*tmp_562 - tmp_562*tmp_95) + tmp_499*(tmp_484*tmp_545 - tmp_498*tmp_563 - tmp_563*tmp_95) + tmp_521*(tmp_506*tmp_545 - tmp_520*tmp_564 - tmp_564*tmp_95) + tmp_543*(tmp_528*tmp_545 - tmp_542*tmp_565 - tmp_565*tmp_95);
      real_t a_2_0 = tmp_103*(-tmp_101*tmp_566 + tmp_46*tmp_567 - tmp_566*tmp_95) + tmp_125*(tmp_110*tmp_567 - tmp_124*tmp_568 - tmp_568*tmp_95) + tmp_147*(tmp_132*tmp_567 - tmp_146*tmp_569 - tmp_569*tmp_95) + tmp_169*(tmp_154*tmp_567 - tmp_168*tmp_570 - tmp_570*tmp_95) + tmp_191*(tmp_176*tmp_567 - tmp_190*tmp_571 - tmp_571*tmp_95) + tmp_213*(tmp_198*tmp_567 - tmp_212*tmp_572 - tmp_572*tmp_95) + tmp_235*(tmp_220*tmp_567 - tmp_234*tmp_573 - tmp_573*tmp_95) + tmp_257*(tmp_242*tmp_567 - tmp_256*tmp_574 - tmp_574*tmp_95) + tmp_279*(tmp_264*tmp_567 - tmp_278*tmp_575 - tmp_575*tmp_95) + tmp_301*(tmp_286*tmp_567 - tmp_300*tmp_576 - tmp_576*tmp_95) + tmp_323*(tmp_308*tmp_567 - tmp_322*tmp_577 - tmp_577*tmp_95) + tmp_345*(tmp_330*tmp_567 - tmp_344*tmp_578 - tmp_578*tmp_95) + tmp_367*(tmp_352*tmp_567 - tmp_366*tmp_579 - tmp_579*tmp_95) + tmp_389*(tmp_374*tmp_567 - tmp_388*tmp_580 - tmp_580*tmp_95) + tmp_411*(tmp_396*tmp_567 - tmp_410*tmp_581 - tmp_581*tmp_95) + tmp_433*(tmp_418*tmp_567 - tmp_432*tmp_582 - tmp_582*tmp_95) + tmp_455*(tmp_440*tmp_567 - tmp_454*tmp_583 - tmp_583*tmp_95) + tmp_477*(tmp_462*tmp_567 - tmp_476*tmp_584 - tmp_584*tmp_95) + tmp_499*(tmp_484*tmp_567 - tmp_498*tmp_585 - tmp_585*tmp_95) + tmp_521*(tmp_506*tmp_567 - tmp_520*tmp_586 - tmp_586*tmp_95) + tmp_543*(tmp_528*tmp_567 - tmp_542*tmp_587 - tmp_587*tmp_95);
      real_t a_3_0 = tmp_103*(-tmp_101*tmp_588 + tmp_46*tmp_589 - tmp_588*tmp_95) + tmp_125*(tmp_110*tmp_589 - tmp_124*tmp_590 - tmp_590*tmp_95) + tmp_147*(tmp_132*tmp_589 - tmp_146*tmp_591 - tmp_591*tmp_95) + tmp_169*(tmp_154*tmp_589 - tmp_168*tmp_592 - tmp_592*tmp_95) + tmp_191*(tmp_176*tmp_589 - tmp_190*tmp_593 - tmp_593*tmp_95) + tmp_213*(tmp_198*tmp_589 - tmp_212*tmp_594 - tmp_594*tmp_95) + tmp_235*(tmp_220*tmp_589 - tmp_234*tmp_595 - tmp_595*tmp_95) + tmp_257*(tmp_242*tmp_589 - tmp_256*tmp_596 - tmp_596*tmp_95) + tmp_279*(tmp_264*tmp_589 - tmp_278*tmp_597 - tmp_597*tmp_95) + tmp_301*(tmp_286*tmp_589 - tmp_300*tmp_598 - tmp_598*tmp_95) + tmp_323*(tmp_308*tmp_589 - tmp_322*tmp_599 - tmp_599*tmp_95) + tmp_345*(tmp_330*tmp_589 - tmp_344*tmp_600 - tmp_600*tmp_95) + tmp_367*(tmp_352*tmp_589 - tmp_366*tmp_601 - tmp_601*tmp_95) + tmp_389*(tmp_374*tmp_589 - tmp_388*tmp_602 - tmp_602*tmp_95) + tmp_411*(tmp_396*tmp_589 - tmp_410*tmp_603 - tmp_603*tmp_95) + tmp_433*(tmp_418*tmp_589 - tmp_432*tmp_604 - tmp_604*tmp_95) + tmp_455*(tmp_440*tmp_589 - tmp_454*tmp_605 - tmp_605*tmp_95) + tmp_477*(tmp_462*tmp_589 - tmp_476*tmp_606 - tmp_606*tmp_95) + tmp_499*(tmp_484*tmp_589 - tmp_498*tmp_607 - tmp_607*tmp_95) + tmp_521*(tmp_506*tmp_589 - tmp_520*tmp_608 - tmp_608*tmp_95) + tmp_543*(tmp_528*tmp_589 - tmp_542*tmp_609 - tmp_609*tmp_95);
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
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_2;
      real_t tmp_8 = p_affine_3_2 + tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_7;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_4;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_11*tmp_14;
      real_t tmp_17 = tmp_4*tmp_8;
      real_t tmp_18 = tmp_13*tmp_3;
      real_t tmp_19 = 1.0 / (-tmp_1*tmp_16 + tmp_1*tmp_9 + tmp_10*tmp_12 - tmp_10*tmp_18 + tmp_15*tmp_5 - tmp_17*tmp_5);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19*(0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22 + tmp_23);
      real_t tmp_25 = tmp_24*tmp_6;
      real_t tmp_26 = -tmp_1*tmp_11 + tmp_13*tmp_5;
      real_t tmp_27 = tmp_24*tmp_26;
      real_t tmp_28 = -tmp_1*tmp_14 + tmp_10*tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19*(0.031405749086161582*tmp_30 + 0.93718850182767688*tmp_31 + tmp_32);
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_8 - tmp_10*tmp_13;
      real_t tmp_36 = tmp_33*tmp_35;
      real_t tmp_37 = tmp_12 - tmp_18;
      real_t tmp_38 = tmp_24*tmp_37;
      real_t tmp_39 = tmp_15 - tmp_17;
      real_t tmp_40 = tmp_33*tmp_39;
      real_t tmp_41 = -tmp_10*tmp_3 + tmp_14*tmp_5;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19*(0.031405749086161582*tmp_43 + 0.93718850182767688*tmp_44 + tmp_45);
      real_t tmp_47 = tmp_41*tmp_46;
      real_t tmp_48 = tmp_10*tmp_11 - tmp_5*tmp_8;
      real_t tmp_49 = tmp_46*tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46*tmp_50;
      real_t tmp_52 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_53 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_54 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_55 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_58 = 5.0*std::pow((std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)*std::abs(tmp_22*tmp_55 - tmp_31*tmp_57)) + (std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)*std::abs(tmp_22*tmp_56 - tmp_44*tmp_57)) + (std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)*std::abs(tmp_31*tmp_56 - tmp_44*tmp_55)), 0.25);
      real_t tmp_59 = 0.0068572537431980923*tmp_58*(tmp_10*(tmp_52 - 1.0/4.0) + tmp_14*(tmp_53 - 1.0/4.0) + tmp_8*(tmp_54 - 1.0/4.0));
      real_t tmp_60 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_61 = tmp_6*tmp_60;
      real_t tmp_62 = tmp_26*tmp_60;
      real_t tmp_63 = tmp_19*(0.19601935860219369*tmp_30 + 0.60796128279561268*tmp_31 + tmp_32);
      real_t tmp_64 = tmp_28*tmp_63;
      real_t tmp_65 = tmp_35*tmp_63;
      real_t tmp_66 = tmp_37*tmp_60;
      real_t tmp_67 = tmp_39*tmp_63;
      real_t tmp_68 = tmp_19*(0.19601935860219369*tmp_43 + 0.60796128279561268*tmp_44 + tmp_45);
      real_t tmp_69 = tmp_41*tmp_68;
      real_t tmp_70 = tmp_48*tmp_68;
      real_t tmp_71 = tmp_50*tmp_68;
      real_t tmp_72 = tmp_66 + tmp_67 + tmp_71;
      real_t tmp_73 = tmp_62 + tmp_65 + tmp_70;
      real_t tmp_74 = tmp_61 + tmp_64 + tmp_69;
      real_t tmp_75 = 0.037198804536718075*tmp_58*(tmp_10*(tmp_72 - 1.0/4.0) + tmp_14*(tmp_73 - 1.0/4.0) + tmp_8*(tmp_74 - 1.0/4.0));
      real_t tmp_76 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_77 = tmp_6*tmp_76;
      real_t tmp_78 = tmp_26*tmp_76;
      real_t tmp_79 = tmp_19*(0.37605877282253791*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_80 = tmp_28*tmp_79;
      real_t tmp_81 = tmp_35*tmp_79;
      real_t tmp_82 = tmp_37*tmp_76;
      real_t tmp_83 = tmp_39*tmp_79;
      real_t tmp_84 = tmp_19*(0.37605877282253791*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_85 = tmp_41*tmp_84;
      real_t tmp_86 = tmp_48*tmp_84;
      real_t tmp_87 = tmp_50*tmp_84;
      real_t tmp_88 = tmp_82 + tmp_83 + tmp_87;
      real_t tmp_89 = tmp_78 + tmp_81 + tmp_86;
      real_t tmp_90 = tmp_77 + tmp_80 + tmp_85;
      real_t tmp_91 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_88 - 1.0/4.0) + tmp_14*(tmp_89 - 1.0/4.0) + tmp_8*(tmp_90 - 1.0/4.0));
      real_t tmp_92 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_93 = tmp_6*tmp_92;
      real_t tmp_94 = tmp_26*tmp_92;
      real_t tmp_95 = tmp_19*(0.78764240869137092*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_96 = tmp_28*tmp_95;
      real_t tmp_97 = tmp_35*tmp_95;
      real_t tmp_98 = tmp_37*tmp_92;
      real_t tmp_99 = tmp_39*tmp_95;
      real_t tmp_100 = tmp_19*(0.78764240869137092*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_101 = tmp_100*tmp_41;
      real_t tmp_102 = tmp_100*tmp_48;
      real_t tmp_103 = tmp_100*tmp_50;
      real_t tmp_104 = tmp_103 + tmp_98 + tmp_99;
      real_t tmp_105 = tmp_102 + tmp_94 + tmp_97;
      real_t tmp_106 = tmp_101 + tmp_93 + tmp_96;
      real_t tmp_107 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_104 - 1.0/4.0) + tmp_14*(tmp_105 - 1.0/4.0) + tmp_8*(tmp_106 - 1.0/4.0));
      real_t tmp_108 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_109 = tmp_108*tmp_6;
      real_t tmp_110 = tmp_108*tmp_26;
      real_t tmp_111 = tmp_19*(0.58463275527740355*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_112 = tmp_111*tmp_28;
      real_t tmp_113 = tmp_111*tmp_35;
      real_t tmp_114 = tmp_108*tmp_37;
      real_t tmp_115 = tmp_111*tmp_39;
      real_t tmp_116 = tmp_19*(0.58463275527740355*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_117 = tmp_116*tmp_41;
      real_t tmp_118 = tmp_116*tmp_48;
      real_t tmp_119 = tmp_116*tmp_50;
      real_t tmp_120 = tmp_114 + tmp_115 + tmp_119;
      real_t tmp_121 = tmp_110 + tmp_113 + tmp_118;
      real_t tmp_122 = tmp_109 + tmp_112 + tmp_117;
      real_t tmp_123 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_120 - 1.0/4.0) + tmp_14*(tmp_121 - 1.0/4.0) + tmp_8*(tmp_122 - 1.0/4.0));
      real_t tmp_124 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_125 = tmp_124*tmp_6;
      real_t tmp_126 = tmp_124*tmp_26;
      real_t tmp_127 = tmp_19*(0.041227165399737475*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_128 = tmp_127*tmp_28;
      real_t tmp_129 = tmp_127*tmp_35;
      real_t tmp_130 = tmp_124*tmp_37;
      real_t tmp_131 = tmp_127*tmp_39;
      real_t tmp_132 = tmp_19*(0.041227165399737475*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_133 = tmp_132*tmp_41;
      real_t tmp_134 = tmp_132*tmp_48;
      real_t tmp_135 = tmp_132*tmp_50;
      real_t tmp_136 = tmp_130 + tmp_131 + tmp_135;
      real_t tmp_137 = tmp_126 + tmp_129 + tmp_134;
      real_t tmp_138 = tmp_125 + tmp_128 + tmp_133;
      real_t tmp_139 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_136 - 1.0/4.0) + tmp_14*(tmp_137 - 1.0/4.0) + tmp_8*(tmp_138 - 1.0/4.0));
      real_t tmp_140 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_141 = tmp_140*tmp_6;
      real_t tmp_142 = tmp_140*tmp_26;
      real_t tmp_143 = tmp_19*(0.039308471900058539*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_144 = tmp_143*tmp_28;
      real_t tmp_145 = tmp_143*tmp_35;
      real_t tmp_146 = tmp_140*tmp_37;
      real_t tmp_147 = tmp_143*tmp_39;
      real_t tmp_148 = tmp_19*(0.039308471900058539*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_149 = tmp_148*tmp_41;
      real_t tmp_150 = tmp_148*tmp_48;
      real_t tmp_151 = tmp_148*tmp_50;
      real_t tmp_152 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_153 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_154 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_155 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_152 - 1.0/4.0) + tmp_14*(tmp_153 - 1.0/4.0) + tmp_8*(tmp_154 - 1.0/4.0));
      real_t tmp_156 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_157 = tmp_156*tmp_6;
      real_t tmp_158 = tmp_156*tmp_26;
      real_t tmp_159 = tmp_19*(0.78764240869137092*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_160 = tmp_159*tmp_28;
      real_t tmp_161 = tmp_159*tmp_35;
      real_t tmp_162 = tmp_156*tmp_37;
      real_t tmp_163 = tmp_159*tmp_39;
      real_t tmp_164 = tmp_19*(0.78764240869137092*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_165 = tmp_164*tmp_41;
      real_t tmp_166 = tmp_164*tmp_48;
      real_t tmp_167 = tmp_164*tmp_50;
      real_t tmp_168 = tmp_162 + tmp_163 + tmp_167;
      real_t tmp_169 = tmp_158 + tmp_161 + tmp_166;
      real_t tmp_170 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_171 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_168 - 1.0/4.0) + tmp_14*(tmp_169 - 1.0/4.0) + tmp_8*(tmp_170 - 1.0/4.0));
      real_t tmp_172 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_173 = tmp_172*tmp_6;
      real_t tmp_174 = tmp_172*tmp_26;
      real_t tmp_175 = tmp_19*(0.58463275527740355*tmp_30 + 0.039308471900058539*tmp_31 + tmp_32);
      real_t tmp_176 = tmp_175*tmp_28;
      real_t tmp_177 = tmp_175*tmp_35;
      real_t tmp_178 = tmp_172*tmp_37;
      real_t tmp_179 = tmp_175*tmp_39;
      real_t tmp_180 = tmp_19*(0.58463275527740355*tmp_43 + 0.039308471900058539*tmp_44 + tmp_45);
      real_t tmp_181 = tmp_180*tmp_41;
      real_t tmp_182 = tmp_180*tmp_48;
      real_t tmp_183 = tmp_180*tmp_50;
      real_t tmp_184 = tmp_178 + tmp_179 + tmp_183;
      real_t tmp_185 = tmp_174 + tmp_177 + tmp_182;
      real_t tmp_186 = tmp_173 + tmp_176 + tmp_181;
      real_t tmp_187 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_184 - 1.0/4.0) + tmp_14*(tmp_185 - 1.0/4.0) + tmp_8*(tmp_186 - 1.0/4.0));
      real_t tmp_188 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_189 = tmp_188*tmp_6;
      real_t tmp_190 = tmp_188*tmp_26;
      real_t tmp_191 = tmp_19*(0.1711304259088916*tmp_30 + 0.78764240869137092*tmp_31 + tmp_32);
      real_t tmp_192 = tmp_191*tmp_28;
      real_t tmp_193 = tmp_191*tmp_35;
      real_t tmp_194 = tmp_188*tmp_37;
      real_t tmp_195 = tmp_191*tmp_39;
      real_t tmp_196 = tmp_19*(0.1711304259088916*tmp_43 + 0.78764240869137092*tmp_44 + tmp_45);
      real_t tmp_197 = tmp_196*tmp_41;
      real_t tmp_198 = tmp_196*tmp_48;
      real_t tmp_199 = tmp_196*tmp_50;
      real_t tmp_200 = tmp_194 + tmp_195 + tmp_199;
      real_t tmp_201 = tmp_190 + tmp_193 + tmp_198;
      real_t tmp_202 = tmp_189 + tmp_192 + tmp_197;
      real_t tmp_203 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_200 - 1.0/4.0) + tmp_14*(tmp_201 - 1.0/4.0) + tmp_8*(tmp_202 - 1.0/4.0));
      real_t tmp_204 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_205 = tmp_204*tmp_6;
      real_t tmp_206 = tmp_204*tmp_26;
      real_t tmp_207 = tmp_19*(0.37605877282253791*tmp_30 + 0.58463275527740355*tmp_31 + tmp_32);
      real_t tmp_208 = tmp_207*tmp_28;
      real_t tmp_209 = tmp_207*tmp_35;
      real_t tmp_210 = tmp_204*tmp_37;
      real_t tmp_211 = tmp_207*tmp_39;
      real_t tmp_212 = tmp_19*(0.37605877282253791*tmp_43 + 0.58463275527740355*tmp_44 + tmp_45);
      real_t tmp_213 = tmp_212*tmp_41;
      real_t tmp_214 = tmp_212*tmp_48;
      real_t tmp_215 = tmp_212*tmp_50;
      real_t tmp_216 = tmp_210 + tmp_211 + tmp_215;
      real_t tmp_217 = tmp_206 + tmp_209 + tmp_214;
      real_t tmp_218 = tmp_205 + tmp_208 + tmp_213;
      real_t tmp_219 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_216 - 1.0/4.0) + tmp_14*(tmp_217 - 1.0/4.0) + tmp_8*(tmp_218 - 1.0/4.0));
      real_t tmp_220 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_221 = tmp_220*tmp_6;
      real_t tmp_222 = tmp_220*tmp_26;
      real_t tmp_223 = tmp_19*(0.041227165399737475*tmp_30 + 0.1711304259088916*tmp_31 + tmp_32);
      real_t tmp_224 = tmp_223*tmp_28;
      real_t tmp_225 = tmp_223*tmp_35;
      real_t tmp_226 = tmp_220*tmp_37;
      real_t tmp_227 = tmp_223*tmp_39;
      real_t tmp_228 = tmp_19*(0.041227165399737475*tmp_43 + 0.1711304259088916*tmp_44 + tmp_45);
      real_t tmp_229 = tmp_228*tmp_41;
      real_t tmp_230 = tmp_228*tmp_48;
      real_t tmp_231 = tmp_228*tmp_50;
      real_t tmp_232 = tmp_226 + tmp_227 + tmp_231;
      real_t tmp_233 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_234 = tmp_221 + tmp_224 + tmp_229;
      real_t tmp_235 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_232 - 1.0/4.0) + tmp_14*(tmp_233 - 1.0/4.0) + tmp_8*(tmp_234 - 1.0/4.0));
      real_t tmp_236 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_237 = tmp_236*tmp_6;
      real_t tmp_238 = tmp_236*tmp_26;
      real_t tmp_239 = tmp_19*(0.40446199974765351*tmp_30 + 0.19107600050469298*tmp_31 + tmp_32);
      real_t tmp_240 = tmp_239*tmp_28;
      real_t tmp_241 = tmp_239*tmp_35;
      real_t tmp_242 = tmp_236*tmp_37;
      real_t tmp_243 = tmp_239*tmp_39;
      real_t tmp_244 = tmp_19*(0.40446199974765351*tmp_43 + 0.19107600050469298*tmp_44 + tmp_45);
      real_t tmp_245 = tmp_244*tmp_41;
      real_t tmp_246 = tmp_244*tmp_48;
      real_t tmp_247 = tmp_244*tmp_50;
      real_t tmp_248 = tmp_242 + tmp_243 + tmp_247;
      real_t tmp_249 = tmp_238 + tmp_241 + tmp_246;
      real_t tmp_250 = tmp_237 + tmp_240 + tmp_245;
      real_t tmp_251 = 0.042507265838595799*tmp_58*(tmp_10*(tmp_248 - 1.0/4.0) + tmp_14*(tmp_249 - 1.0/4.0) + tmp_8*(tmp_250 - 1.0/4.0));
      real_t tmp_252 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_253 = tmp_252*tmp_6;
      real_t tmp_254 = tmp_252*tmp_26;
      real_t tmp_255 = tmp_19*(0.039308471900058539*tmp_30 + 0.37605877282253791*tmp_31 + tmp_32);
      real_t tmp_256 = tmp_255*tmp_28;
      real_t tmp_257 = tmp_255*tmp_35;
      real_t tmp_258 = tmp_252*tmp_37;
      real_t tmp_259 = tmp_255*tmp_39;
      real_t tmp_260 = tmp_19*(0.039308471900058539*tmp_43 + 0.37605877282253791*tmp_44 + tmp_45);
      real_t tmp_261 = tmp_260*tmp_41;
      real_t tmp_262 = tmp_260*tmp_48;
      real_t tmp_263 = tmp_260*tmp_50;
      real_t tmp_264 = tmp_258 + tmp_259 + tmp_263;
      real_t tmp_265 = tmp_254 + tmp_257 + tmp_262;
      real_t tmp_266 = tmp_253 + tmp_256 + tmp_261;
      real_t tmp_267 = 0.020848748529055869*tmp_58*(tmp_10*(tmp_264 - 1.0/4.0) + tmp_14*(tmp_265 - 1.0/4.0) + tmp_8*(tmp_266 - 1.0/4.0));
      real_t tmp_268 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_269 = tmp_268*tmp_6;
      real_t tmp_270 = tmp_26*tmp_268;
      real_t tmp_271 = tmp_19*(0.93718850182767688*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_272 = tmp_271*tmp_28;
      real_t tmp_273 = tmp_271*tmp_35;
      real_t tmp_274 = tmp_268*tmp_37;
      real_t tmp_275 = tmp_271*tmp_39;
      real_t tmp_276 = tmp_19*(0.93718850182767688*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_277 = tmp_276*tmp_41;
      real_t tmp_278 = tmp_276*tmp_48;
      real_t tmp_279 = tmp_276*tmp_50;
      real_t tmp_280 = tmp_274 + tmp_275 + tmp_279;
      real_t tmp_281 = tmp_270 + tmp_273 + tmp_278;
      real_t tmp_282 = tmp_269 + tmp_272 + tmp_277;
      real_t tmp_283 = 0.0068572537431980923*tmp_58*(tmp_10*(tmp_280 - 1.0/4.0) + tmp_14*(tmp_281 - 1.0/4.0) + tmp_8*(tmp_282 - 1.0/4.0));
      real_t tmp_284 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_285 = tmp_284*tmp_6;
      real_t tmp_286 = tmp_26*tmp_284;
      real_t tmp_287 = tmp_19*(0.60796128279561268*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_288 = tmp_28*tmp_287;
      real_t tmp_289 = tmp_287*tmp_35;
      real_t tmp_290 = tmp_284*tmp_37;
      real_t tmp_291 = tmp_287*tmp_39;
      real_t tmp_292 = tmp_19*(0.60796128279561268*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_293 = tmp_292*tmp_41;
      real_t tmp_294 = tmp_292*tmp_48;
      real_t tmp_295 = tmp_292*tmp_50;
      real_t tmp_296 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_297 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_298 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_299 = 0.037198804536718075*tmp_58*(tmp_10*(tmp_296 - 1.0/4.0) + tmp_14*(tmp_297 - 1.0/4.0) + tmp_8*(tmp_298 - 1.0/4.0));
      real_t tmp_300 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_301 = tmp_300*tmp_6;
      real_t tmp_302 = tmp_26*tmp_300;
      real_t tmp_303 = tmp_19*(0.19107600050469298*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_304 = tmp_28*tmp_303;
      real_t tmp_305 = tmp_303*tmp_35;
      real_t tmp_306 = tmp_300*tmp_37;
      real_t tmp_307 = tmp_303*tmp_39;
      real_t tmp_308 = tmp_19*(0.19107600050469298*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_309 = tmp_308*tmp_41;
      real_t tmp_310 = tmp_308*tmp_48;
      real_t tmp_311 = tmp_308*tmp_50;
      real_t tmp_312 = tmp_306 + tmp_307 + tmp_311;
      real_t tmp_313 = tmp_302 + tmp_305 + tmp_310;
      real_t tmp_314 = tmp_301 + tmp_304 + tmp_309;
      real_t tmp_315 = 0.042507265838595799*tmp_58*(tmp_10*(tmp_312 - 1.0/4.0) + tmp_14*(tmp_313 - 1.0/4.0) + tmp_8*(tmp_314 - 1.0/4.0));
      real_t tmp_316 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_317 = tmp_316*tmp_6;
      real_t tmp_318 = tmp_26*tmp_316;
      real_t tmp_319 = tmp_19*(0.031405749086161582*tmp_30 + 0.031405749086161582*tmp_31 + tmp_32);
      real_t tmp_320 = tmp_28*tmp_319;
      real_t tmp_321 = tmp_319*tmp_35;
      real_t tmp_322 = tmp_316*tmp_37;
      real_t tmp_323 = tmp_319*tmp_39;
      real_t tmp_324 = tmp_19*(0.031405749086161582*tmp_43 + 0.031405749086161582*tmp_44 + tmp_45);
      real_t tmp_325 = tmp_324*tmp_41;
      real_t tmp_326 = tmp_324*tmp_48;
      real_t tmp_327 = tmp_324*tmp_50;
      real_t tmp_328 = tmp_322 + tmp_323 + tmp_327;
      real_t tmp_329 = tmp_318 + tmp_321 + tmp_326;
      real_t tmp_330 = tmp_317 + tmp_320 + tmp_325;
      real_t tmp_331 = 0.0068572537431980923*tmp_58*(tmp_10*(tmp_328 - 1.0/4.0) + tmp_14*(tmp_329 - 1.0/4.0) + tmp_8*(tmp_330 - 1.0/4.0));
      real_t tmp_332 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_333 = tmp_332*tmp_6;
      real_t tmp_334 = tmp_26*tmp_332;
      real_t tmp_335 = tmp_19*(0.19601935860219369*tmp_30 + 0.19601935860219369*tmp_31 + tmp_32);
      real_t tmp_336 = tmp_28*tmp_335;
      real_t tmp_337 = tmp_335*tmp_35;
      real_t tmp_338 = tmp_332*tmp_37;
      real_t tmp_339 = tmp_335*tmp_39;
      real_t tmp_340 = tmp_19*(0.19601935860219369*tmp_43 + 0.19601935860219369*tmp_44 + tmp_45);
      real_t tmp_341 = tmp_340*tmp_41;
      real_t tmp_342 = tmp_340*tmp_48;
      real_t tmp_343 = tmp_340*tmp_50;
      real_t tmp_344 = tmp_338 + tmp_339 + tmp_343;
      real_t tmp_345 = tmp_334 + tmp_337 + tmp_342;
      real_t tmp_346 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_347 = 0.037198804536718075*tmp_58*(tmp_10*(tmp_344 - 1.0/4.0) + tmp_14*(tmp_345 - 1.0/4.0) + tmp_8*(tmp_346 - 1.0/4.0));
      real_t tmp_348 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_349 = tmp_348*tmp_6;
      real_t tmp_350 = tmp_26*tmp_348;
      real_t tmp_351 = tmp_19*(0.40446199974765351*tmp_30 + 0.40446199974765351*tmp_31 + tmp_32);
      real_t tmp_352 = tmp_28*tmp_351;
      real_t tmp_353 = tmp_35*tmp_351;
      real_t tmp_354 = tmp_348*tmp_37;
      real_t tmp_355 = tmp_351*tmp_39;
      real_t tmp_356 = tmp_19*(0.40446199974765351*tmp_43 + 0.40446199974765351*tmp_44 + tmp_45);
      real_t tmp_357 = tmp_356*tmp_41;
      real_t tmp_358 = tmp_356*tmp_48;
      real_t tmp_359 = tmp_356*tmp_50;
      real_t tmp_360 = tmp_354 + tmp_355 + tmp_359;
      real_t tmp_361 = tmp_350 + tmp_353 + tmp_358;
      real_t tmp_362 = tmp_349 + tmp_352 + tmp_357;
      real_t tmp_363 = 0.042507265838595799*tmp_58*(tmp_10*(tmp_360 - 1.0/4.0) + tmp_14*(tmp_361 - 1.0/4.0) + tmp_8*(tmp_362 - 1.0/4.0));
      real_t tmp_364 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_365 = tmp_364*tmp_6;
      real_t tmp_366 = tmp_26*tmp_364;
      real_t tmp_367 = tmp_19*(0.1711304259088916*tmp_30 + 0.041227165399737475*tmp_31 + tmp_32);
      real_t tmp_368 = tmp_28*tmp_367;
      real_t tmp_369 = tmp_35*tmp_367;
      real_t tmp_370 = tmp_364*tmp_37;
      real_t tmp_371 = tmp_367*tmp_39;
      real_t tmp_372 = tmp_19*(0.1711304259088916*tmp_43 + 0.041227165399737475*tmp_44 + tmp_45);
      real_t tmp_373 = tmp_372*tmp_41;
      real_t tmp_374 = tmp_372*tmp_48;
      real_t tmp_375 = tmp_372*tmp_50;
      real_t tmp_376 = tmp_370 + tmp_371 + tmp_375;
      real_t tmp_377 = tmp_366 + tmp_369 + tmp_374;
      real_t tmp_378 = tmp_365 + tmp_368 + tmp_373;
      real_t tmp_379 = 0.019202922745021479*tmp_58*(tmp_10*(tmp_376 - 1.0/4.0) + tmp_14*(tmp_377 - 1.0/4.0) + tmp_8*(tmp_378 - 1.0/4.0));
      real_t a_0_0 = tmp_107*(-tmp_101 - tmp_102 - tmp_103 - tmp_93 - tmp_94 - tmp_96 - tmp_97 - tmp_98 - tmp_99 + 1) + tmp_123*(-tmp_109 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_115 - tmp_117 - tmp_118 - tmp_119 + 1) + tmp_139*(-tmp_125 - tmp_126 - tmp_128 - tmp_129 - tmp_130 - tmp_131 - tmp_133 - tmp_134 - tmp_135 + 1) + tmp_155*(-tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1) + tmp_171*(-tmp_157 - tmp_158 - tmp_160 - tmp_161 - tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 + 1) + tmp_187*(-tmp_173 - tmp_174 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 + 1) + tmp_203*(-tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_197 - tmp_198 - tmp_199 + 1) + tmp_219*(-tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 - tmp_211 - tmp_213 - tmp_214 - tmp_215 + 1) + tmp_235*(-tmp_221 - tmp_222 - tmp_224 - tmp_225 - tmp_226 - tmp_227 - tmp_229 - tmp_230 - tmp_231 + 1) + tmp_251*(-tmp_237 - tmp_238 - tmp_240 - tmp_241 - tmp_242 - tmp_243 - tmp_245 - tmp_246 - tmp_247 + 1) + tmp_267*(-tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1) + tmp_283*(-tmp_269 - tmp_270 - tmp_272 - tmp_273 - tmp_274 - tmp_275 - tmp_277 - tmp_278 - tmp_279 + 1) + tmp_299*(-tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1) + tmp_315*(-tmp_301 - tmp_302 - tmp_304 - tmp_305 - tmp_306 - tmp_307 - tmp_309 - tmp_310 - tmp_311 + 1) + tmp_331*(-tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1) + tmp_347*(-tmp_333 - tmp_334 - tmp_336 - tmp_337 - tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 + 1) + tmp_363*(-tmp_349 - tmp_350 - tmp_352 - tmp_353 - tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 + 1) + tmp_379*(-tmp_365 - tmp_366 - tmp_368 - tmp_369 - tmp_370 - tmp_371 - tmp_373 - tmp_374 - tmp_375 + 1) + tmp_59*(-tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1) + tmp_75*(-tmp_61 - tmp_62 - tmp_64 - tmp_65 - tmp_66 - tmp_67 - tmp_69 - tmp_70 - tmp_71 + 1) + tmp_91*(-tmp_77 - tmp_78 - tmp_80 - tmp_81 - tmp_82 - tmp_83 - tmp_85 - tmp_86 - tmp_87 + 1);
      real_t a_1_0 = tmp_104*tmp_107 + tmp_120*tmp_123 + tmp_136*tmp_139 + tmp_152*tmp_155 + tmp_168*tmp_171 + tmp_184*tmp_187 + tmp_200*tmp_203 + tmp_216*tmp_219 + tmp_232*tmp_235 + tmp_248*tmp_251 + tmp_264*tmp_267 + tmp_280*tmp_283 + tmp_296*tmp_299 + tmp_312*tmp_315 + tmp_328*tmp_331 + tmp_344*tmp_347 + tmp_360*tmp_363 + tmp_376*tmp_379 + tmp_52*tmp_59 + tmp_72*tmp_75 + tmp_88*tmp_91;
      real_t a_2_0 = tmp_105*tmp_107 + tmp_121*tmp_123 + tmp_137*tmp_139 + tmp_153*tmp_155 + tmp_169*tmp_171 + tmp_185*tmp_187 + tmp_201*tmp_203 + tmp_217*tmp_219 + tmp_233*tmp_235 + tmp_249*tmp_251 + tmp_265*tmp_267 + tmp_281*tmp_283 + tmp_297*tmp_299 + tmp_313*tmp_315 + tmp_329*tmp_331 + tmp_345*tmp_347 + tmp_361*tmp_363 + tmp_377*tmp_379 + tmp_53*tmp_59 + tmp_73*tmp_75 + tmp_89*tmp_91;
      real_t a_3_0 = tmp_106*tmp_107 + tmp_122*tmp_123 + tmp_138*tmp_139 + tmp_154*tmp_155 + tmp_170*tmp_171 + tmp_186*tmp_187 + tmp_202*tmp_203 + tmp_218*tmp_219 + tmp_234*tmp_235 + tmp_250*tmp_251 + tmp_266*tmp_267 + tmp_282*tmp_283 + tmp_298*tmp_299 + tmp_314*tmp_315 + tmp_330*tmp_331 + tmp_346*tmp_347 + tmp_362*tmp_363 + tmp_378*tmp_379 + tmp_54*tmp_59 + tmp_74*tmp_75 + tmp_90*tmp_91;
      elMat( 0, 0) = a_0_0;
      elMat( 1, 0) = a_1_0;
      elMat( 2, 0) = a_2_0;
      elMat( 3, 0) = a_3_0;
   }


};




class EGNIPGVectorLaplaceFormEP1_0 : public hyteg::dg::DGForm
{
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
      real_t tmp_28 = 5/tmp_27;
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
      real_t tmp_37 = 5/tmp_36;
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

      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_0;
      real_t tmp_3 = p_affine_1_0 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_2;
      real_t tmp_5 = 1.0 / (tmp_1*tmp_3 - tmp_4*(p_affine_1_1 + tmp_0));
      real_t tmp_6 = tmp_1*tmp_5;
      real_t tmp_7 = tmp_5*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_8 = tmp_3*tmp_5;
      real_t tmp_9 = tmp_5*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_10 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_11 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_12 = std::abs(std::pow((tmp_10*tmp_10) + (tmp_11*tmp_11), 1.0/2.0));
      real_t tmp_13 = tmp_12*(p_affine_10_0*(-tmp_6 - tmp_7) + p_affine_10_1*(-tmp_8 - tmp_9));
      real_t tmp_14 = p_affine_6_1 + tmp_0;
      real_t tmp_15 = 0.046910077030668018*tmp_11 + tmp_14;
      real_t tmp_16 = p_affine_6_0 + tmp_2;
      real_t tmp_17 = 0.046910077030668018*tmp_10 + tmp_16;
      real_t tmp_18 = 0.11846344252809471*tmp_3*(tmp_15*tmp_9 + tmp_17*tmp_6 - 1.0/3.0) + 0.11846344252809471*tmp_4*(tmp_15*tmp_8 + tmp_17*tmp_7 - 1.0/3.0);
      real_t tmp_19 = 0.23076534494715845*tmp_11 + tmp_14;
      real_t tmp_20 = 0.23076534494715845*tmp_10 + tmp_16;
      real_t tmp_21 = 0.2393143352496831*tmp_3*(tmp_19*tmp_9 + tmp_20*tmp_6 - 1.0/3.0) + 0.2393143352496831*tmp_4*(tmp_19*tmp_8 + tmp_20*tmp_7 - 1.0/3.0);
      real_t tmp_22 = 0.5*tmp_11 + tmp_14;
      real_t tmp_23 = 0.5*tmp_10 + tmp_16;
      real_t tmp_24 = 0.2844444444444445*tmp_3*(tmp_22*tmp_9 + tmp_23*tmp_6 - 1.0/3.0) + 0.2844444444444445*tmp_4*(tmp_22*tmp_8 + tmp_23*tmp_7 - 1.0/3.0);
      real_t tmp_25 = 0.7692346550528415*tmp_11 + tmp_14;
      real_t tmp_26 = 0.7692346550528415*tmp_10 + tmp_16;
      real_t tmp_27 = 0.2393143352496831*tmp_3*(tmp_25*tmp_9 + tmp_26*tmp_6 - 1.0/3.0) + 0.2393143352496831*tmp_4*(tmp_25*tmp_8 + tmp_26*tmp_7 - 1.0/3.0);
      real_t tmp_28 = 0.95308992296933193*tmp_11 + tmp_14;
      real_t tmp_29 = 0.95308992296933193*tmp_10 + tmp_16;
      real_t tmp_30 = 0.11846344252809471*tmp_3*(tmp_28*tmp_9 + tmp_29*tmp_6 - 1.0/3.0) + 0.11846344252809471*tmp_4*(tmp_28*tmp_8 + tmp_29*tmp_7 - 1.0/3.0);
      real_t tmp_31 = tmp_12*(p_affine_10_0*tmp_6 + p_affine_10_1*tmp_9);
      real_t tmp_32 = tmp_12*(p_affine_10_0*tmp_7 + p_affine_10_1*tmp_8);
      real_t a_0_0 = -tmp_13*tmp_18 - tmp_13*tmp_21 - tmp_13*tmp_24 - tmp_13*tmp_27 - tmp_13*tmp_30;
      real_t a_0_1 = -tmp_18*tmp_31 - tmp_21*tmp_31 - tmp_24*tmp_31 - tmp_27*tmp_31 - tmp_30*tmp_31;
      real_t a_0_2 = -tmp_18*tmp_32 - tmp_21*tmp_32 - tmp_24*tmp_32 - tmp_27*tmp_32 - tmp_30*tmp_32;
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
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
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
      real_t tmp_22 = tmp_1*tmp_21 + tmp_14*tmp_19 + tmp_20*tmp_5;
      real_t tmp_23 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_5);
      real_t tmp_24 = tmp_18*(tmp_1*tmp_9 - tmp_10*tmp_14);
      real_t tmp_25 = tmp_18*(tmp_13*tmp_14 - tmp_5*tmp_9);
      real_t tmp_26 = tmp_1*tmp_25 + tmp_14*tmp_23 + tmp_24*tmp_5;
      real_t tmp_27 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_6);
      real_t tmp_28 = tmp_18*(tmp_10*tmp_11 - tmp_6*tmp_9);
      real_t tmp_29 = tmp_18*(-tmp_11*tmp_13 + tmp_3*tmp_9);
      real_t tmp_30 = tmp_1*tmp_29 + tmp_14*tmp_27 + tmp_28*tmp_5;
      real_t tmp_31 = p_affine_0_0*p_affine_1_1;
      real_t tmp_32 = p_affine_0_0*p_affine_1_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_2;
      real_t tmp_34 = p_affine_0_1*p_affine_1_0;
      real_t tmp_35 = p_affine_0_1*p_affine_1_2;
      real_t tmp_36 = p_affine_2_2*p_affine_3_0;
      real_t tmp_37 = p_affine_0_2*p_affine_1_0;
      real_t tmp_38 = p_affine_0_2*p_affine_1_1;
      real_t tmp_39 = p_affine_2_0*p_affine_3_1;
      real_t tmp_40 = p_affine_2_2*p_affine_3_1;
      real_t tmp_41 = p_affine_2_0*p_affine_3_2;
      real_t tmp_42 = p_affine_2_1*p_affine_3_0;
      real_t tmp_43 = std::abs(p_affine_0_0*tmp_33 - p_affine_0_0*tmp_40 + p_affine_0_1*tmp_36 - p_affine_0_1*tmp_41 + p_affine_0_2*tmp_39 - p_affine_0_2*tmp_42 - p_affine_1_0*tmp_33 + p_affine_1_0*tmp_40 - p_affine_1_1*tmp_36 + p_affine_1_1*tmp_41 - p_affine_1_2*tmp_39 + p_affine_1_2*tmp_42 + p_affine_2_0*tmp_35 - p_affine_2_0*tmp_38 - p_affine_2_1*tmp_32 + p_affine_2_1*tmp_37 + p_affine_2_2*tmp_31 - p_affine_2_2*tmp_34 - p_affine_3_0*tmp_35 + p_affine_3_0*tmp_38 + p_affine_3_1*tmp_32 - p_affine_3_1*tmp_37 - p_affine_3_2*tmp_31 + p_affine_3_2*tmp_34);
      real_t tmp_44 = tmp_43*(tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25) + tmp_30*(-tmp_27 - tmp_28 - tmp_29));
      real_t tmp_45 = tmp_43*(tmp_21*tmp_22 + tmp_25*tmp_26 + tmp_29*tmp_30);
      real_t tmp_46 = tmp_43*(tmp_20*tmp_22 + tmp_24*tmp_26 + tmp_28*tmp_30);
      real_t tmp_47 = tmp_43*(tmp_19*tmp_22 + tmp_23*tmp_26 + tmp_27*tmp_30);
      real_t a_0_0 = 0.1666666666666668*tmp_44;
      real_t a_0_1 = 0.1666666666666668*tmp_45;
      real_t a_0_2 = 0.1666666666666668*tmp_46;
      real_t a_0_3 = 0.1666666666666668*tmp_47;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
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
      real_t tmp_2 = -p_affine_8_2;
      real_t tmp_3 = p_affine_9_2 + tmp_2;
      real_t tmp_4 = p_affine_10_2 + tmp_2;
      real_t tmp_5 = -p_affine_0_2;
      real_t tmp_6 = p_affine_8_2 + tmp_5;
      real_t tmp_7 = 0.031405749086161582*tmp_3 + 0.93718850182767688*tmp_4 + tmp_6;
      real_t tmp_8 = p_affine_2_0 + tmp_0;
      real_t tmp_9 = -p_affine_0_1;
      real_t tmp_10 = p_affine_3_1 + tmp_9;
      real_t tmp_11 = p_affine_3_0 + tmp_0;
      real_t tmp_12 = p_affine_2_1 + tmp_9;
      real_t tmp_13 = p_affine_3_2 + tmp_5;
      real_t tmp_14 = tmp_12*tmp_13;
      real_t tmp_15 = p_affine_1_2 + tmp_5;
      real_t tmp_16 = tmp_10*tmp_15;
      real_t tmp_17 = p_affine_1_1 + tmp_9;
      real_t tmp_18 = p_affine_2_2 + tmp_5;
      real_t tmp_19 = tmp_17*tmp_18;
      real_t tmp_20 = tmp_10*tmp_18;
      real_t tmp_21 = tmp_13*tmp_17;
      real_t tmp_22 = tmp_12*tmp_15;
      real_t tmp_23 = 1.0 / (tmp_1*tmp_14 - tmp_1*tmp_20 + tmp_11*tmp_19 - tmp_11*tmp_22 + tmp_16*tmp_8 - tmp_21*tmp_8);
      real_t tmp_24 = tmp_23*(tmp_10*tmp_8 - tmp_11*tmp_12);
      real_t tmp_25 = tmp_24*tmp_7;
      real_t tmp_26 = -p_affine_8_1;
      real_t tmp_27 = p_affine_9_1 + tmp_26;
      real_t tmp_28 = p_affine_10_1 + tmp_26;
      real_t tmp_29 = p_affine_8_1 + tmp_9;
      real_t tmp_30 = 0.031405749086161582*tmp_27 + 0.93718850182767688*tmp_28 + tmp_29;
      real_t tmp_31 = tmp_23*(tmp_11*tmp_18 - tmp_13*tmp_8);
      real_t tmp_32 = tmp_30*tmp_31;
      real_t tmp_33 = -p_affine_8_0;
      real_t tmp_34 = p_affine_9_0 + tmp_33;
      real_t tmp_35 = p_affine_10_0 + tmp_33;
      real_t tmp_36 = p_affine_8_0 + tmp_0;
      real_t tmp_37 = 0.031405749086161582*tmp_34 + 0.93718850182767688*tmp_35 + tmp_36;
      real_t tmp_38 = tmp_23*(tmp_14 - tmp_20);
      real_t tmp_39 = tmp_37*tmp_38;
      real_t tmp_40 = tmp_25 + tmp_32 + tmp_39;
      real_t tmp_41 = tmp_23*(-tmp_1*tmp_10 + tmp_11*tmp_17);
      real_t tmp_42 = tmp_41*tmp_7;
      real_t tmp_43 = tmp_23*(tmp_1*tmp_13 - tmp_11*tmp_15);
      real_t tmp_44 = tmp_30*tmp_43;
      real_t tmp_45 = tmp_23*(tmp_16 - tmp_21);
      real_t tmp_46 = tmp_37*tmp_45;
      real_t tmp_47 = tmp_42 + tmp_44 + tmp_46;
      real_t tmp_48 = tmp_23*(tmp_1*tmp_12 - tmp_17*tmp_8);
      real_t tmp_49 = tmp_48*tmp_7;
      real_t tmp_50 = tmp_23*(-tmp_1*tmp_18 + tmp_15*tmp_8);
      real_t tmp_51 = tmp_30*tmp_50;
      real_t tmp_52 = tmp_23*(tmp_19 - tmp_22);
      real_t tmp_53 = tmp_37*tmp_52;
      real_t tmp_54 = tmp_49 + tmp_51 + tmp_53;
      real_t tmp_55 = tmp_1*(tmp_40 - 1.0/4.0) + tmp_11*(tmp_54 - 1.0/4.0) + tmp_8*(tmp_47 - 1.0/4.0);
      real_t tmp_56 = 0.5*p_affine_13_0*(-tmp_38 - tmp_45 - tmp_52) + 0.5*p_affine_13_1*(-tmp_31 - tmp_43 - tmp_50) + 0.5*p_affine_13_2*(-tmp_24 - tmp_41 - tmp_48);
      real_t tmp_57 = -tmp_25 - tmp_32 - tmp_39 - tmp_42 - tmp_44 - tmp_46 - tmp_49 - tmp_51 - tmp_53 + 1;
      real_t tmp_58 = 0.5*p_affine_13_0*(tmp_1*tmp_38 + tmp_11*tmp_52 + tmp_45*tmp_8) + 0.5*p_affine_13_1*(tmp_1*tmp_31 + tmp_11*tmp_50 + tmp_43*tmp_8) + 0.5*p_affine_13_2*(tmp_1*tmp_24 + tmp_11*tmp_48 + tmp_41*tmp_8);
      real_t tmp_59 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_60 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_61 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_62 = (std::abs(tmp_28*tmp_60 - tmp_35*tmp_59)*std::abs(tmp_28*tmp_60 - tmp_35*tmp_59)) + (std::abs(tmp_28*tmp_61 - tmp_4*tmp_59)*std::abs(tmp_28*tmp_61 - tmp_4*tmp_59)) + (std::abs(tmp_35*tmp_61 - tmp_4*tmp_60)*std::abs(tmp_35*tmp_61 - tmp_4*tmp_60));
      real_t tmp_63 = 5.0*std::pow(tmp_62, -0.25);
      real_t tmp_64 = tmp_55*tmp_63;
      real_t tmp_65 = 1.0*std::pow(tmp_62, 1.0/2.0);
      real_t tmp_66 = 0.0068572537431980923*tmp_65;
      real_t tmp_67 = 0.19601935860219369*tmp_3 + 0.60796128279561268*tmp_4 + tmp_6;
      real_t tmp_68 = tmp_24*tmp_67;
      real_t tmp_69 = 0.19601935860219369*tmp_27 + 0.60796128279561268*tmp_28 + tmp_29;
      real_t tmp_70 = tmp_31*tmp_69;
      real_t tmp_71 = 0.19601935860219369*tmp_34 + 0.60796128279561268*tmp_35 + tmp_36;
      real_t tmp_72 = tmp_38*tmp_71;
      real_t tmp_73 = tmp_68 + tmp_70 + tmp_72;
      real_t tmp_74 = tmp_41*tmp_67;
      real_t tmp_75 = tmp_43*tmp_69;
      real_t tmp_76 = tmp_45*tmp_71;
      real_t tmp_77 = tmp_74 + tmp_75 + tmp_76;
      real_t tmp_78 = tmp_48*tmp_67;
      real_t tmp_79 = tmp_50*tmp_69;
      real_t tmp_80 = tmp_52*tmp_71;
      real_t tmp_81 = tmp_78 + tmp_79 + tmp_80;
      real_t tmp_82 = tmp_1*(tmp_73 - 1.0/4.0) + tmp_11*(tmp_81 - 1.0/4.0) + tmp_8*(tmp_77 - 1.0/4.0);
      real_t tmp_83 = -tmp_68 - tmp_70 - tmp_72 - tmp_74 - tmp_75 - tmp_76 - tmp_78 - tmp_79 - tmp_80 + 1;
      real_t tmp_84 = tmp_63*tmp_82;
      real_t tmp_85 = 0.037198804536718075*tmp_65;
      real_t tmp_86 = 0.37605877282253791*tmp_3 + 0.039308471900058539*tmp_4 + tmp_6;
      real_t tmp_87 = tmp_24*tmp_86;
      real_t tmp_88 = 0.37605877282253791*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29;
      real_t tmp_89 = tmp_31*tmp_88;
      real_t tmp_90 = 0.37605877282253791*tmp_34 + 0.039308471900058539*tmp_35 + tmp_36;
      real_t tmp_91 = tmp_38*tmp_90;
      real_t tmp_92 = tmp_87 + tmp_89 + tmp_91;
      real_t tmp_93 = tmp_41*tmp_86;
      real_t tmp_94 = tmp_43*tmp_88;
      real_t tmp_95 = tmp_45*tmp_90;
      real_t tmp_96 = tmp_93 + tmp_94 + tmp_95;
      real_t tmp_97 = tmp_48*tmp_86;
      real_t tmp_98 = tmp_50*tmp_88;
      real_t tmp_99 = tmp_52*tmp_90;
      real_t tmp_100 = tmp_97 + tmp_98 + tmp_99;
      real_t tmp_101 = tmp_1*(tmp_92 - 1.0/4.0) + tmp_11*(tmp_100 - 1.0/4.0) + tmp_8*(tmp_96 - 1.0/4.0);
      real_t tmp_102 = -tmp_87 - tmp_89 - tmp_91 - tmp_93 - tmp_94 - tmp_95 - tmp_97 - tmp_98 - tmp_99 + 1;
      real_t tmp_103 = tmp_101*tmp_63;
      real_t tmp_104 = 0.020848748529055869*tmp_65;
      real_t tmp_105 = 0.78764240869137092*tmp_3 + 0.1711304259088916*tmp_4 + tmp_6;
      real_t tmp_106 = tmp_105*tmp_24;
      real_t tmp_107 = 0.78764240869137092*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29;
      real_t tmp_108 = tmp_107*tmp_31;
      real_t tmp_109 = 0.78764240869137092*tmp_34 + 0.1711304259088916*tmp_35 + tmp_36;
      real_t tmp_110 = tmp_109*tmp_38;
      real_t tmp_111 = tmp_106 + tmp_108 + tmp_110;
      real_t tmp_112 = tmp_105*tmp_41;
      real_t tmp_113 = tmp_107*tmp_43;
      real_t tmp_114 = tmp_109*tmp_45;
      real_t tmp_115 = tmp_112 + tmp_113 + tmp_114;
      real_t tmp_116 = tmp_105*tmp_48;
      real_t tmp_117 = tmp_107*tmp_50;
      real_t tmp_118 = tmp_109*tmp_52;
      real_t tmp_119 = tmp_116 + tmp_117 + tmp_118;
      real_t tmp_120 = tmp_1*(tmp_111 - 1.0/4.0) + tmp_11*(tmp_119 - 1.0/4.0) + tmp_8*(tmp_115 - 1.0/4.0);
      real_t tmp_121 = -tmp_106 - tmp_108 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_116 - tmp_117 - tmp_118 + 1;
      real_t tmp_122 = tmp_120*tmp_63;
      real_t tmp_123 = 0.019202922745021479*tmp_65;
      real_t tmp_124 = 0.58463275527740355*tmp_3 + 0.37605877282253791*tmp_4 + tmp_6;
      real_t tmp_125 = tmp_124*tmp_24;
      real_t tmp_126 = 0.58463275527740355*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29;
      real_t tmp_127 = tmp_126*tmp_31;
      real_t tmp_128 = 0.58463275527740355*tmp_34 + 0.37605877282253791*tmp_35 + tmp_36;
      real_t tmp_129 = tmp_128*tmp_38;
      real_t tmp_130 = tmp_125 + tmp_127 + tmp_129;
      real_t tmp_131 = tmp_124*tmp_41;
      real_t tmp_132 = tmp_126*tmp_43;
      real_t tmp_133 = tmp_128*tmp_45;
      real_t tmp_134 = tmp_131 + tmp_132 + tmp_133;
      real_t tmp_135 = tmp_124*tmp_48;
      real_t tmp_136 = tmp_126*tmp_50;
      real_t tmp_137 = tmp_128*tmp_52;
      real_t tmp_138 = tmp_135 + tmp_136 + tmp_137;
      real_t tmp_139 = tmp_1*(tmp_130 - 1.0/4.0) + tmp_11*(tmp_138 - 1.0/4.0) + tmp_8*(tmp_134 - 1.0/4.0);
      real_t tmp_140 = -tmp_125 - tmp_127 - tmp_129 - tmp_131 - tmp_132 - tmp_133 - tmp_135 - tmp_136 - tmp_137 + 1;
      real_t tmp_141 = tmp_139*tmp_63;
      real_t tmp_142 = 0.020848748529055869*tmp_65;
      real_t tmp_143 = 0.041227165399737475*tmp_3 + 0.78764240869137092*tmp_4 + tmp_6;
      real_t tmp_144 = tmp_143*tmp_24;
      real_t tmp_145 = 0.041227165399737475*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29;
      real_t tmp_146 = tmp_145*tmp_31;
      real_t tmp_147 = 0.041227165399737475*tmp_34 + 0.78764240869137092*tmp_35 + tmp_36;
      real_t tmp_148 = tmp_147*tmp_38;
      real_t tmp_149 = tmp_144 + tmp_146 + tmp_148;
      real_t tmp_150 = tmp_143*tmp_41;
      real_t tmp_151 = tmp_145*tmp_43;
      real_t tmp_152 = tmp_147*tmp_45;
      real_t tmp_153 = tmp_150 + tmp_151 + tmp_152;
      real_t tmp_154 = tmp_143*tmp_48;
      real_t tmp_155 = tmp_145*tmp_50;
      real_t tmp_156 = tmp_147*tmp_52;
      real_t tmp_157 = tmp_154 + tmp_155 + tmp_156;
      real_t tmp_158 = tmp_1*(tmp_149 - 1.0/4.0) + tmp_11*(tmp_157 - 1.0/4.0) + tmp_8*(tmp_153 - 1.0/4.0);
      real_t tmp_159 = -tmp_144 - tmp_146 - tmp_148 - tmp_150 - tmp_151 - tmp_152 - tmp_154 - tmp_155 - tmp_156 + 1;
      real_t tmp_160 = tmp_158*tmp_63;
      real_t tmp_161 = 0.019202922745021479*tmp_65;
      real_t tmp_162 = 0.039308471900058539*tmp_3 + 0.58463275527740355*tmp_4 + tmp_6;
      real_t tmp_163 = tmp_162*tmp_24;
      real_t tmp_164 = 0.039308471900058539*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29;
      real_t tmp_165 = tmp_164*tmp_31;
      real_t tmp_166 = 0.039308471900058539*tmp_34 + 0.58463275527740355*tmp_35 + tmp_36;
      real_t tmp_167 = tmp_166*tmp_38;
      real_t tmp_168 = tmp_163 + tmp_165 + tmp_167;
      real_t tmp_169 = tmp_162*tmp_41;
      real_t tmp_170 = tmp_164*tmp_43;
      real_t tmp_171 = tmp_166*tmp_45;
      real_t tmp_172 = tmp_169 + tmp_170 + tmp_171;
      real_t tmp_173 = tmp_162*tmp_48;
      real_t tmp_174 = tmp_164*tmp_50;
      real_t tmp_175 = tmp_166*tmp_52;
      real_t tmp_176 = tmp_173 + tmp_174 + tmp_175;
      real_t tmp_177 = tmp_1*(tmp_168 - 1.0/4.0) + tmp_11*(tmp_176 - 1.0/4.0) + tmp_8*(tmp_172 - 1.0/4.0);
      real_t tmp_178 = -tmp_163 - tmp_165 - tmp_167 - tmp_169 - tmp_170 - tmp_171 - tmp_173 - tmp_174 - tmp_175 + 1;
      real_t tmp_179 = tmp_177*tmp_63;
      real_t tmp_180 = 0.020848748529055869*tmp_65;
      real_t tmp_181 = 0.78764240869137092*tmp_3 + 0.041227165399737475*tmp_4 + tmp_6;
      real_t tmp_182 = tmp_181*tmp_24;
      real_t tmp_183 = 0.78764240869137092*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29;
      real_t tmp_184 = tmp_183*tmp_31;
      real_t tmp_185 = 0.78764240869137092*tmp_34 + 0.041227165399737475*tmp_35 + tmp_36;
      real_t tmp_186 = tmp_185*tmp_38;
      real_t tmp_187 = tmp_182 + tmp_184 + tmp_186;
      real_t tmp_188 = tmp_181*tmp_41;
      real_t tmp_189 = tmp_183*tmp_43;
      real_t tmp_190 = tmp_185*tmp_45;
      real_t tmp_191 = tmp_188 + tmp_189 + tmp_190;
      real_t tmp_192 = tmp_181*tmp_48;
      real_t tmp_193 = tmp_183*tmp_50;
      real_t tmp_194 = tmp_185*tmp_52;
      real_t tmp_195 = tmp_192 + tmp_193 + tmp_194;
      real_t tmp_196 = tmp_1*(tmp_187 - 1.0/4.0) + tmp_11*(tmp_195 - 1.0/4.0) + tmp_8*(tmp_191 - 1.0/4.0);
      real_t tmp_197 = -tmp_182 - tmp_184 - tmp_186 - tmp_188 - tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 + 1;
      real_t tmp_198 = tmp_196*tmp_63;
      real_t tmp_199 = 0.019202922745021479*tmp_65;
      real_t tmp_200 = 0.58463275527740355*tmp_3 + 0.039308471900058539*tmp_4 + tmp_6;
      real_t tmp_201 = tmp_200*tmp_24;
      real_t tmp_202 = 0.58463275527740355*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29;
      real_t tmp_203 = tmp_202*tmp_31;
      real_t tmp_204 = 0.58463275527740355*tmp_34 + 0.039308471900058539*tmp_35 + tmp_36;
      real_t tmp_205 = tmp_204*tmp_38;
      real_t tmp_206 = tmp_201 + tmp_203 + tmp_205;
      real_t tmp_207 = tmp_200*tmp_41;
      real_t tmp_208 = tmp_202*tmp_43;
      real_t tmp_209 = tmp_204*tmp_45;
      real_t tmp_210 = tmp_207 + tmp_208 + tmp_209;
      real_t tmp_211 = tmp_200*tmp_48;
      real_t tmp_212 = tmp_202*tmp_50;
      real_t tmp_213 = tmp_204*tmp_52;
      real_t tmp_214 = tmp_211 + tmp_212 + tmp_213;
      real_t tmp_215 = tmp_1*(tmp_206 - 1.0/4.0) + tmp_11*(tmp_214 - 1.0/4.0) + tmp_8*(tmp_210 - 1.0/4.0);
      real_t tmp_216 = -tmp_201 - tmp_203 - tmp_205 - tmp_207 - tmp_208 - tmp_209 - tmp_211 - tmp_212 - tmp_213 + 1;
      real_t tmp_217 = tmp_215*tmp_63;
      real_t tmp_218 = 0.020848748529055869*tmp_65;
      real_t tmp_219 = 0.1711304259088916*tmp_3 + 0.78764240869137092*tmp_4 + tmp_6;
      real_t tmp_220 = tmp_219*tmp_24;
      real_t tmp_221 = 0.1711304259088916*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29;
      real_t tmp_222 = tmp_221*tmp_31;
      real_t tmp_223 = 0.1711304259088916*tmp_34 + 0.78764240869137092*tmp_35 + tmp_36;
      real_t tmp_224 = tmp_223*tmp_38;
      real_t tmp_225 = tmp_220 + tmp_222 + tmp_224;
      real_t tmp_226 = tmp_219*tmp_41;
      real_t tmp_227 = tmp_221*tmp_43;
      real_t tmp_228 = tmp_223*tmp_45;
      real_t tmp_229 = tmp_226 + tmp_227 + tmp_228;
      real_t tmp_230 = tmp_219*tmp_48;
      real_t tmp_231 = tmp_221*tmp_50;
      real_t tmp_232 = tmp_223*tmp_52;
      real_t tmp_233 = tmp_230 + tmp_231 + tmp_232;
      real_t tmp_234 = tmp_1*(tmp_225 - 1.0/4.0) + tmp_11*(tmp_233 - 1.0/4.0) + tmp_8*(tmp_229 - 1.0/4.0);
      real_t tmp_235 = -tmp_220 - tmp_222 - tmp_224 - tmp_226 - tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 + 1;
      real_t tmp_236 = tmp_234*tmp_63;
      real_t tmp_237 = 0.019202922745021479*tmp_65;
      real_t tmp_238 = 0.37605877282253791*tmp_3 + 0.58463275527740355*tmp_4 + tmp_6;
      real_t tmp_239 = tmp_238*tmp_24;
      real_t tmp_240 = 0.37605877282253791*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29;
      real_t tmp_241 = tmp_240*tmp_31;
      real_t tmp_242 = 0.37605877282253791*tmp_34 + 0.58463275527740355*tmp_35 + tmp_36;
      real_t tmp_243 = tmp_242*tmp_38;
      real_t tmp_244 = tmp_239 + tmp_241 + tmp_243;
      real_t tmp_245 = tmp_238*tmp_41;
      real_t tmp_246 = tmp_240*tmp_43;
      real_t tmp_247 = tmp_242*tmp_45;
      real_t tmp_248 = tmp_245 + tmp_246 + tmp_247;
      real_t tmp_249 = tmp_238*tmp_48;
      real_t tmp_250 = tmp_240*tmp_50;
      real_t tmp_251 = tmp_242*tmp_52;
      real_t tmp_252 = tmp_249 + tmp_250 + tmp_251;
      real_t tmp_253 = tmp_1*(tmp_244 - 1.0/4.0) + tmp_11*(tmp_252 - 1.0/4.0) + tmp_8*(tmp_248 - 1.0/4.0);
      real_t tmp_254 = -tmp_239 - tmp_241 - tmp_243 - tmp_245 - tmp_246 - tmp_247 - tmp_249 - tmp_250 - tmp_251 + 1;
      real_t tmp_255 = tmp_253*tmp_63;
      real_t tmp_256 = 0.020848748529055869*tmp_65;
      real_t tmp_257 = 0.041227165399737475*tmp_3 + 0.1711304259088916*tmp_4 + tmp_6;
      real_t tmp_258 = tmp_24*tmp_257;
      real_t tmp_259 = 0.041227165399737475*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29;
      real_t tmp_260 = tmp_259*tmp_31;
      real_t tmp_261 = 0.041227165399737475*tmp_34 + 0.1711304259088916*tmp_35 + tmp_36;
      real_t tmp_262 = tmp_261*tmp_38;
      real_t tmp_263 = tmp_258 + tmp_260 + tmp_262;
      real_t tmp_264 = tmp_257*tmp_41;
      real_t tmp_265 = tmp_259*tmp_43;
      real_t tmp_266 = tmp_261*tmp_45;
      real_t tmp_267 = tmp_264 + tmp_265 + tmp_266;
      real_t tmp_268 = tmp_257*tmp_48;
      real_t tmp_269 = tmp_259*tmp_50;
      real_t tmp_270 = tmp_261*tmp_52;
      real_t tmp_271 = tmp_268 + tmp_269 + tmp_270;
      real_t tmp_272 = tmp_1*(tmp_263 - 1.0/4.0) + tmp_11*(tmp_271 - 1.0/4.0) + tmp_8*(tmp_267 - 1.0/4.0);
      real_t tmp_273 = -tmp_258 - tmp_260 - tmp_262 - tmp_264 - tmp_265 - tmp_266 - tmp_268 - tmp_269 - tmp_270 + 1;
      real_t tmp_274 = tmp_272*tmp_63;
      real_t tmp_275 = 0.019202922745021479*tmp_65;
      real_t tmp_276 = 0.40446199974765351*tmp_3 + 0.19107600050469298*tmp_4 + tmp_6;
      real_t tmp_277 = tmp_24*tmp_276;
      real_t tmp_278 = 0.40446199974765351*tmp_27 + 0.19107600050469298*tmp_28 + tmp_29;
      real_t tmp_279 = tmp_278*tmp_31;
      real_t tmp_280 = 0.40446199974765351*tmp_34 + 0.19107600050469298*tmp_35 + tmp_36;
      real_t tmp_281 = tmp_280*tmp_38;
      real_t tmp_282 = tmp_277 + tmp_279 + tmp_281;
      real_t tmp_283 = tmp_276*tmp_41;
      real_t tmp_284 = tmp_278*tmp_43;
      real_t tmp_285 = tmp_280*tmp_45;
      real_t tmp_286 = tmp_283 + tmp_284 + tmp_285;
      real_t tmp_287 = tmp_276*tmp_48;
      real_t tmp_288 = tmp_278*tmp_50;
      real_t tmp_289 = tmp_280*tmp_52;
      real_t tmp_290 = tmp_287 + tmp_288 + tmp_289;
      real_t tmp_291 = tmp_1*(tmp_282 - 1.0/4.0) + tmp_11*(tmp_290 - 1.0/4.0) + tmp_8*(tmp_286 - 1.0/4.0);
      real_t tmp_292 = -tmp_277 - tmp_279 - tmp_281 - tmp_283 - tmp_284 - tmp_285 - tmp_287 - tmp_288 - tmp_289 + 1;
      real_t tmp_293 = tmp_291*tmp_63;
      real_t tmp_294 = 0.042507265838595799*tmp_65;
      real_t tmp_295 = 0.039308471900058539*tmp_3 + 0.37605877282253791*tmp_4 + tmp_6;
      real_t tmp_296 = tmp_24*tmp_295;
      real_t tmp_297 = 0.039308471900058539*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29;
      real_t tmp_298 = tmp_297*tmp_31;
      real_t tmp_299 = 0.039308471900058539*tmp_34 + 0.37605877282253791*tmp_35 + tmp_36;
      real_t tmp_300 = tmp_299*tmp_38;
      real_t tmp_301 = tmp_296 + tmp_298 + tmp_300;
      real_t tmp_302 = tmp_295*tmp_41;
      real_t tmp_303 = tmp_297*tmp_43;
      real_t tmp_304 = tmp_299*tmp_45;
      real_t tmp_305 = tmp_302 + tmp_303 + tmp_304;
      real_t tmp_306 = tmp_295*tmp_48;
      real_t tmp_307 = tmp_297*tmp_50;
      real_t tmp_308 = tmp_299*tmp_52;
      real_t tmp_309 = tmp_306 + tmp_307 + tmp_308;
      real_t tmp_310 = tmp_1*(tmp_301 - 1.0/4.0) + tmp_11*(tmp_309 - 1.0/4.0) + tmp_8*(tmp_305 - 1.0/4.0);
      real_t tmp_311 = -tmp_296 - tmp_298 - tmp_300 - tmp_302 - tmp_303 - tmp_304 - tmp_306 - tmp_307 - tmp_308 + 1;
      real_t tmp_312 = tmp_310*tmp_63;
      real_t tmp_313 = 0.020848748529055869*tmp_65;
      real_t tmp_314 = 0.93718850182767688*tmp_3 + 0.031405749086161582*tmp_4 + tmp_6;
      real_t tmp_315 = tmp_24*tmp_314;
      real_t tmp_316 = 0.93718850182767688*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29;
      real_t tmp_317 = tmp_31*tmp_316;
      real_t tmp_318 = 0.93718850182767688*tmp_34 + 0.031405749086161582*tmp_35 + tmp_36;
      real_t tmp_319 = tmp_318*tmp_38;
      real_t tmp_320 = tmp_315 + tmp_317 + tmp_319;
      real_t tmp_321 = tmp_314*tmp_41;
      real_t tmp_322 = tmp_316*tmp_43;
      real_t tmp_323 = tmp_318*tmp_45;
      real_t tmp_324 = tmp_321 + tmp_322 + tmp_323;
      real_t tmp_325 = tmp_314*tmp_48;
      real_t tmp_326 = tmp_316*tmp_50;
      real_t tmp_327 = tmp_318*tmp_52;
      real_t tmp_328 = tmp_325 + tmp_326 + tmp_327;
      real_t tmp_329 = tmp_1*(tmp_320 - 1.0/4.0) + tmp_11*(tmp_328 - 1.0/4.0) + tmp_8*(tmp_324 - 1.0/4.0);
      real_t tmp_330 = -tmp_315 - tmp_317 - tmp_319 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1;
      real_t tmp_331 = tmp_329*tmp_63;
      real_t tmp_332 = 0.0068572537431980923*tmp_65;
      real_t tmp_333 = 0.60796128279561268*tmp_3 + 0.19601935860219369*tmp_4 + tmp_6;
      real_t tmp_334 = tmp_24*tmp_333;
      real_t tmp_335 = 0.60796128279561268*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29;
      real_t tmp_336 = tmp_31*tmp_335;
      real_t tmp_337 = 0.60796128279561268*tmp_34 + 0.19601935860219369*tmp_35 + tmp_36;
      real_t tmp_338 = tmp_337*tmp_38;
      real_t tmp_339 = tmp_334 + tmp_336 + tmp_338;
      real_t tmp_340 = tmp_333*tmp_41;
      real_t tmp_341 = tmp_335*tmp_43;
      real_t tmp_342 = tmp_337*tmp_45;
      real_t tmp_343 = tmp_340 + tmp_341 + tmp_342;
      real_t tmp_344 = tmp_333*tmp_48;
      real_t tmp_345 = tmp_335*tmp_50;
      real_t tmp_346 = tmp_337*tmp_52;
      real_t tmp_347 = tmp_344 + tmp_345 + tmp_346;
      real_t tmp_348 = tmp_1*(tmp_339 - 1.0/4.0) + tmp_11*(tmp_347 - 1.0/4.0) + tmp_8*(tmp_343 - 1.0/4.0);
      real_t tmp_349 = -tmp_334 - tmp_336 - tmp_338 - tmp_340 - tmp_341 - tmp_342 - tmp_344 - tmp_345 - tmp_346 + 1;
      real_t tmp_350 = tmp_348*tmp_63;
      real_t tmp_351 = 0.037198804536718075*tmp_65;
      real_t tmp_352 = 0.19107600050469298*tmp_3 + 0.40446199974765351*tmp_4 + tmp_6;
      real_t tmp_353 = tmp_24*tmp_352;
      real_t tmp_354 = 0.19107600050469298*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29;
      real_t tmp_355 = tmp_31*tmp_354;
      real_t tmp_356 = 0.19107600050469298*tmp_34 + 0.40446199974765351*tmp_35 + tmp_36;
      real_t tmp_357 = tmp_356*tmp_38;
      real_t tmp_358 = tmp_353 + tmp_355 + tmp_357;
      real_t tmp_359 = tmp_352*tmp_41;
      real_t tmp_360 = tmp_354*tmp_43;
      real_t tmp_361 = tmp_356*tmp_45;
      real_t tmp_362 = tmp_359 + tmp_360 + tmp_361;
      real_t tmp_363 = tmp_352*tmp_48;
      real_t tmp_364 = tmp_354*tmp_50;
      real_t tmp_365 = tmp_356*tmp_52;
      real_t tmp_366 = tmp_363 + tmp_364 + tmp_365;
      real_t tmp_367 = tmp_1*(tmp_358 - 1.0/4.0) + tmp_11*(tmp_366 - 1.0/4.0) + tmp_8*(tmp_362 - 1.0/4.0);
      real_t tmp_368 = -tmp_353 - tmp_355 - tmp_357 - tmp_359 - tmp_360 - tmp_361 - tmp_363 - tmp_364 - tmp_365 + 1;
      real_t tmp_369 = tmp_367*tmp_63;
      real_t tmp_370 = 0.042507265838595799*tmp_65;
      real_t tmp_371 = 0.031405749086161582*tmp_3 + 0.031405749086161582*tmp_4 + tmp_6;
      real_t tmp_372 = tmp_24*tmp_371;
      real_t tmp_373 = 0.031405749086161582*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29;
      real_t tmp_374 = tmp_31*tmp_373;
      real_t tmp_375 = 0.031405749086161582*tmp_34 + 0.031405749086161582*tmp_35 + tmp_36;
      real_t tmp_376 = tmp_375*tmp_38;
      real_t tmp_377 = tmp_372 + tmp_374 + tmp_376;
      real_t tmp_378 = tmp_371*tmp_41;
      real_t tmp_379 = tmp_373*tmp_43;
      real_t tmp_380 = tmp_375*tmp_45;
      real_t tmp_381 = tmp_378 + tmp_379 + tmp_380;
      real_t tmp_382 = tmp_371*tmp_48;
      real_t tmp_383 = tmp_373*tmp_50;
      real_t tmp_384 = tmp_375*tmp_52;
      real_t tmp_385 = tmp_382 + tmp_383 + tmp_384;
      real_t tmp_386 = tmp_1*(tmp_377 - 1.0/4.0) + tmp_11*(tmp_385 - 1.0/4.0) + tmp_8*(tmp_381 - 1.0/4.0);
      real_t tmp_387 = -tmp_372 - tmp_374 - tmp_376 - tmp_378 - tmp_379 - tmp_380 - tmp_382 - tmp_383 - tmp_384 + 1;
      real_t tmp_388 = tmp_386*tmp_63;
      real_t tmp_389 = 0.0068572537431980923*tmp_65;
      real_t tmp_390 = 0.19601935860219369*tmp_3 + 0.19601935860219369*tmp_4 + tmp_6;
      real_t tmp_391 = tmp_24*tmp_390;
      real_t tmp_392 = 0.19601935860219369*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29;
      real_t tmp_393 = tmp_31*tmp_392;
      real_t tmp_394 = 0.19601935860219369*tmp_34 + 0.19601935860219369*tmp_35 + tmp_36;
      real_t tmp_395 = tmp_38*tmp_394;
      real_t tmp_396 = tmp_391 + tmp_393 + tmp_395;
      real_t tmp_397 = tmp_390*tmp_41;
      real_t tmp_398 = tmp_392*tmp_43;
      real_t tmp_399 = tmp_394*tmp_45;
      real_t tmp_400 = tmp_397 + tmp_398 + tmp_399;
      real_t tmp_401 = tmp_390*tmp_48;
      real_t tmp_402 = tmp_392*tmp_50;
      real_t tmp_403 = tmp_394*tmp_52;
      real_t tmp_404 = tmp_401 + tmp_402 + tmp_403;
      real_t tmp_405 = tmp_1*(tmp_396 - 1.0/4.0) + tmp_11*(tmp_404 - 1.0/4.0) + tmp_8*(tmp_400 - 1.0/4.0);
      real_t tmp_406 = -tmp_391 - tmp_393 - tmp_395 - tmp_397 - tmp_398 - tmp_399 - tmp_401 - tmp_402 - tmp_403 + 1;
      real_t tmp_407 = tmp_405*tmp_63;
      real_t tmp_408 = 0.037198804536718075*tmp_65;
      real_t tmp_409 = 0.40446199974765351*tmp_3 + 0.40446199974765351*tmp_4 + tmp_6;
      real_t tmp_410 = tmp_24*tmp_409;
      real_t tmp_411 = 0.40446199974765351*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29;
      real_t tmp_412 = tmp_31*tmp_411;
      real_t tmp_413 = 0.40446199974765351*tmp_34 + 0.40446199974765351*tmp_35 + tmp_36;
      real_t tmp_414 = tmp_38*tmp_413;
      real_t tmp_415 = tmp_410 + tmp_412 + tmp_414;
      real_t tmp_416 = tmp_409*tmp_41;
      real_t tmp_417 = tmp_411*tmp_43;
      real_t tmp_418 = tmp_413*tmp_45;
      real_t tmp_419 = tmp_416 + tmp_417 + tmp_418;
      real_t tmp_420 = tmp_409*tmp_48;
      real_t tmp_421 = tmp_411*tmp_50;
      real_t tmp_422 = tmp_413*tmp_52;
      real_t tmp_423 = tmp_420 + tmp_421 + tmp_422;
      real_t tmp_424 = tmp_1*(tmp_415 - 1.0/4.0) + tmp_11*(tmp_423 - 1.0/4.0) + tmp_8*(tmp_419 - 1.0/4.0);
      real_t tmp_425 = -tmp_410 - tmp_412 - tmp_414 - tmp_416 - tmp_417 - tmp_418 - tmp_420 - tmp_421 - tmp_422 + 1;
      real_t tmp_426 = tmp_424*tmp_63;
      real_t tmp_427 = 0.042507265838595799*tmp_65;
      real_t tmp_428 = 0.1711304259088916*tmp_3 + 0.041227165399737475*tmp_4 + tmp_6;
      real_t tmp_429 = tmp_24*tmp_428;
      real_t tmp_430 = 0.1711304259088916*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29;
      real_t tmp_431 = tmp_31*tmp_430;
      real_t tmp_432 = 0.1711304259088916*tmp_34 + 0.041227165399737475*tmp_35 + tmp_36;
      real_t tmp_433 = tmp_38*tmp_432;
      real_t tmp_434 = tmp_429 + tmp_431 + tmp_433;
      real_t tmp_435 = tmp_41*tmp_428;
      real_t tmp_436 = tmp_43*tmp_430;
      real_t tmp_437 = tmp_432*tmp_45;
      real_t tmp_438 = tmp_435 + tmp_436 + tmp_437;
      real_t tmp_439 = tmp_428*tmp_48;
      real_t tmp_440 = tmp_430*tmp_50;
      real_t tmp_441 = tmp_432*tmp_52;
      real_t tmp_442 = tmp_439 + tmp_440 + tmp_441;
      real_t tmp_443 = tmp_1*(tmp_434 - 1.0/4.0) + tmp_11*(tmp_442 - 1.0/4.0) + tmp_8*(tmp_438 - 1.0/4.0);
      real_t tmp_444 = -tmp_429 - tmp_431 - tmp_433 - tmp_435 - tmp_436 - tmp_437 - tmp_439 - tmp_440 - tmp_441 + 1;
      real_t tmp_445 = tmp_443*tmp_63;
      real_t tmp_446 = 0.019202922745021479*tmp_65;
      real_t tmp_447 = 0.5*p_affine_13_0*tmp_38 + 0.5*p_affine_13_1*tmp_31 + 0.5*p_affine_13_2*tmp_24;
      real_t tmp_448 = 0.5*p_affine_13_0*tmp_45 + 0.5*p_affine_13_1*tmp_43 + 0.5*p_affine_13_2*tmp_41;
      real_t tmp_449 = 0.5*p_affine_13_0*tmp_52 + 0.5*p_affine_13_1*tmp_50 + 0.5*p_affine_13_2*tmp_48;
      real_t a_0_0 = tmp_104*(-tmp_101*tmp_56 + tmp_102*tmp_103 - tmp_102*tmp_58) + tmp_123*(-tmp_120*tmp_56 + tmp_121*tmp_122 - tmp_121*tmp_58) + tmp_142*(-tmp_139*tmp_56 + tmp_140*tmp_141 - tmp_140*tmp_58) + tmp_161*(-tmp_158*tmp_56 + tmp_159*tmp_160 - tmp_159*tmp_58) + tmp_180*(-tmp_177*tmp_56 + tmp_178*tmp_179 - tmp_178*tmp_58) + tmp_199*(-tmp_196*tmp_56 + tmp_197*tmp_198 - tmp_197*tmp_58) + tmp_218*(-tmp_215*tmp_56 + tmp_216*tmp_217 - tmp_216*tmp_58) + tmp_237*(-tmp_234*tmp_56 + tmp_235*tmp_236 - tmp_235*tmp_58) + tmp_256*(-tmp_253*tmp_56 + tmp_254*tmp_255 - tmp_254*tmp_58) + tmp_275*(-tmp_272*tmp_56 + tmp_273*tmp_274 - tmp_273*tmp_58) + tmp_294*(-tmp_291*tmp_56 + tmp_292*tmp_293 - tmp_292*tmp_58) + tmp_313*(-tmp_310*tmp_56 + tmp_311*tmp_312 - tmp_311*tmp_58) + tmp_332*(-tmp_329*tmp_56 + tmp_330*tmp_331 - tmp_330*tmp_58) + tmp_351*(-tmp_348*tmp_56 + tmp_349*tmp_350 - tmp_349*tmp_58) + tmp_370*(-tmp_367*tmp_56 + tmp_368*tmp_369 - tmp_368*tmp_58) + tmp_389*(-tmp_386*tmp_56 + tmp_387*tmp_388 - tmp_387*tmp_58) + tmp_408*(-tmp_405*tmp_56 + tmp_406*tmp_407 - tmp_406*tmp_58) + tmp_427*(-tmp_424*tmp_56 + tmp_425*tmp_426 - tmp_425*tmp_58) + tmp_446*(-tmp_443*tmp_56 + tmp_444*tmp_445 - tmp_444*tmp_58) + tmp_66*(-tmp_55*tmp_56 - tmp_57*tmp_58 + tmp_57*tmp_64) + tmp_85*(-tmp_56*tmp_82 - tmp_58*tmp_83 + tmp_83*tmp_84);
      real_t a_0_1 = tmp_104*(-tmp_101*tmp_447 + tmp_103*tmp_92 - tmp_58*tmp_92) + tmp_123*(tmp_111*tmp_122 - tmp_111*tmp_58 - tmp_120*tmp_447) + tmp_142*(tmp_130*tmp_141 - tmp_130*tmp_58 - tmp_139*tmp_447) + tmp_161*(tmp_149*tmp_160 - tmp_149*tmp_58 - tmp_158*tmp_447) + tmp_180*(tmp_168*tmp_179 - tmp_168*tmp_58 - tmp_177*tmp_447) + tmp_199*(tmp_187*tmp_198 - tmp_187*tmp_58 - tmp_196*tmp_447) + tmp_218*(tmp_206*tmp_217 - tmp_206*tmp_58 - tmp_215*tmp_447) + tmp_237*(tmp_225*tmp_236 - tmp_225*tmp_58 - tmp_234*tmp_447) + tmp_256*(tmp_244*tmp_255 - tmp_244*tmp_58 - tmp_253*tmp_447) + tmp_275*(tmp_263*tmp_274 - tmp_263*tmp_58 - tmp_272*tmp_447) + tmp_294*(tmp_282*tmp_293 - tmp_282*tmp_58 - tmp_291*tmp_447) + tmp_313*(tmp_301*tmp_312 - tmp_301*tmp_58 - tmp_310*tmp_447) + tmp_332*(tmp_320*tmp_331 - tmp_320*tmp_58 - tmp_329*tmp_447) + tmp_351*(tmp_339*tmp_350 - tmp_339*tmp_58 - tmp_348*tmp_447) + tmp_370*(tmp_358*tmp_369 - tmp_358*tmp_58 - tmp_367*tmp_447) + tmp_389*(tmp_377*tmp_388 - tmp_377*tmp_58 - tmp_386*tmp_447) + tmp_408*(tmp_396*tmp_407 - tmp_396*tmp_58 - tmp_405*tmp_447) + tmp_427*(tmp_415*tmp_426 - tmp_415*tmp_58 - tmp_424*tmp_447) + tmp_446*(tmp_434*tmp_445 - tmp_434*tmp_58 - tmp_443*tmp_447) + tmp_66*(-tmp_40*tmp_58 + tmp_40*tmp_64 - tmp_447*tmp_55) + tmp_85*(-tmp_447*tmp_82 - tmp_58*tmp_73 + tmp_73*tmp_84);
      real_t a_0_2 = tmp_104*(-tmp_101*tmp_448 + tmp_103*tmp_96 - tmp_58*tmp_96) + tmp_123*(tmp_115*tmp_122 - tmp_115*tmp_58 - tmp_120*tmp_448) + tmp_142*(tmp_134*tmp_141 - tmp_134*tmp_58 - tmp_139*tmp_448) + tmp_161*(tmp_153*tmp_160 - tmp_153*tmp_58 - tmp_158*tmp_448) + tmp_180*(tmp_172*tmp_179 - tmp_172*tmp_58 - tmp_177*tmp_448) + tmp_199*(tmp_191*tmp_198 - tmp_191*tmp_58 - tmp_196*tmp_448) + tmp_218*(tmp_210*tmp_217 - tmp_210*tmp_58 - tmp_215*tmp_448) + tmp_237*(tmp_229*tmp_236 - tmp_229*tmp_58 - tmp_234*tmp_448) + tmp_256*(tmp_248*tmp_255 - tmp_248*tmp_58 - tmp_253*tmp_448) + tmp_275*(tmp_267*tmp_274 - tmp_267*tmp_58 - tmp_272*tmp_448) + tmp_294*(tmp_286*tmp_293 - tmp_286*tmp_58 - tmp_291*tmp_448) + tmp_313*(tmp_305*tmp_312 - tmp_305*tmp_58 - tmp_310*tmp_448) + tmp_332*(tmp_324*tmp_331 - tmp_324*tmp_58 - tmp_329*tmp_448) + tmp_351*(tmp_343*tmp_350 - tmp_343*tmp_58 - tmp_348*tmp_448) + tmp_370*(tmp_362*tmp_369 - tmp_362*tmp_58 - tmp_367*tmp_448) + tmp_389*(tmp_381*tmp_388 - tmp_381*tmp_58 - tmp_386*tmp_448) + tmp_408*(tmp_400*tmp_407 - tmp_400*tmp_58 - tmp_405*tmp_448) + tmp_427*(tmp_419*tmp_426 - tmp_419*tmp_58 - tmp_424*tmp_448) + tmp_446*(tmp_438*tmp_445 - tmp_438*tmp_58 - tmp_443*tmp_448) + tmp_66*(-tmp_448*tmp_55 - tmp_47*tmp_58 + tmp_47*tmp_64) + tmp_85*(-tmp_448*tmp_82 - tmp_58*tmp_77 + tmp_77*tmp_84);
      real_t a_0_3 = tmp_104*(tmp_100*tmp_103 - tmp_100*tmp_58 - tmp_101*tmp_449) + tmp_123*(tmp_119*tmp_122 - tmp_119*tmp_58 - tmp_120*tmp_449) + tmp_142*(tmp_138*tmp_141 - tmp_138*tmp_58 - tmp_139*tmp_449) + tmp_161*(tmp_157*tmp_160 - tmp_157*tmp_58 - tmp_158*tmp_449) + tmp_180*(tmp_176*tmp_179 - tmp_176*tmp_58 - tmp_177*tmp_449) + tmp_199*(tmp_195*tmp_198 - tmp_195*tmp_58 - tmp_196*tmp_449) + tmp_218*(tmp_214*tmp_217 - tmp_214*tmp_58 - tmp_215*tmp_449) + tmp_237*(tmp_233*tmp_236 - tmp_233*tmp_58 - tmp_234*tmp_449) + tmp_256*(tmp_252*tmp_255 - tmp_252*tmp_58 - tmp_253*tmp_449) + tmp_275*(tmp_271*tmp_274 - tmp_271*tmp_58 - tmp_272*tmp_449) + tmp_294*(tmp_290*tmp_293 - tmp_290*tmp_58 - tmp_291*tmp_449) + tmp_313*(tmp_309*tmp_312 - tmp_309*tmp_58 - tmp_310*tmp_449) + tmp_332*(tmp_328*tmp_331 - tmp_328*tmp_58 - tmp_329*tmp_449) + tmp_351*(tmp_347*tmp_350 - tmp_347*tmp_58 - tmp_348*tmp_449) + tmp_370*(tmp_366*tmp_369 - tmp_366*tmp_58 - tmp_367*tmp_449) + tmp_389*(tmp_385*tmp_388 - tmp_385*tmp_58 - tmp_386*tmp_449) + tmp_408*(tmp_404*tmp_407 - tmp_404*tmp_58 - tmp_405*tmp_449) + tmp_427*(tmp_423*tmp_426 - tmp_423*tmp_58 - tmp_424*tmp_449) + tmp_446*(tmp_442*tmp_445 - tmp_442*tmp_58 - tmp_443*tmp_449) + tmp_66*(-tmp_449*tmp_55 - tmp_54*tmp_58 + tmp_54*tmp_64) + tmp_85*(-tmp_449*tmp_82 - tmp_58*tmp_81 + tmp_81*tmp_84);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
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
      real_t tmp_5 = p_affine_3_0 + tmp_0;
      real_t tmp_6 = p_affine_2_1 + tmp_3;
      real_t tmp_7 = tmp_2*tmp_4 - tmp_5*tmp_6;
      real_t tmp_8 = -p_affine_0_2;
      real_t tmp_9 = p_affine_3_2 + tmp_8;
      real_t tmp_10 = tmp_6*tmp_9;
      real_t tmp_11 = p_affine_1_2 + tmp_8;
      real_t tmp_12 = tmp_11*tmp_4;
      real_t tmp_13 = p_affine_1_1 + tmp_3;
      real_t tmp_14 = p_affine_2_2 + tmp_8;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_14*tmp_4;
      real_t tmp_17 = tmp_13*tmp_9;
      real_t tmp_18 = tmp_11*tmp_6;
      real_t tmp_19 = 1.0 / (tmp_1*tmp_10 - tmp_1*tmp_16 + tmp_12*tmp_2 + tmp_15*tmp_5 - tmp_17*tmp_2 - tmp_18*tmp_5);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = 0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22;
      real_t tmp_24 = p_affine_8_2 + tmp_8;
      real_t tmp_25 = tmp_19*(tmp_23 + tmp_24);
      real_t tmp_26 = tmp_14*tmp_5 - tmp_2*tmp_9;
      real_t tmp_27 = -p_affine_8_1;
      real_t tmp_28 = p_affine_9_1 + tmp_27;
      real_t tmp_29 = p_affine_10_1 + tmp_27;
      real_t tmp_30 = 0.031405749086161582*tmp_28 + 0.93718850182767688*tmp_29;
      real_t tmp_31 = p_affine_8_1 + tmp_3;
      real_t tmp_32 = tmp_19*(tmp_30 + tmp_31);
      real_t tmp_33 = tmp_10 - tmp_16;
      real_t tmp_34 = -p_affine_8_0;
      real_t tmp_35 = p_affine_9_0 + tmp_34;
      real_t tmp_36 = p_affine_10_0 + tmp_34;
      real_t tmp_37 = 0.031405749086161582*tmp_35 + 0.93718850182767688*tmp_36;
      real_t tmp_38 = p_affine_8_0 + tmp_0;
      real_t tmp_39 = tmp_19*(tmp_37 + tmp_38);
      real_t tmp_40 = -tmp_1*tmp_4 + tmp_13*tmp_5;
      real_t tmp_41 = tmp_1*tmp_9 - tmp_11*tmp_5;
      real_t tmp_42 = tmp_12 - tmp_17;
      real_t tmp_43 = tmp_1*tmp_6 - tmp_13*tmp_2;
      real_t tmp_44 = -tmp_1*tmp_14 + tmp_11*tmp_2;
      real_t tmp_45 = tmp_15 - tmp_18;
      real_t tmp_46 = tmp_1*(tmp_25*tmp_7 + tmp_26*tmp_32 + tmp_33*tmp_39 - 1.0/4.0) + tmp_2*(tmp_25*tmp_40 + tmp_32*tmp_41 + tmp_39*tmp_42 - 1.0/4.0) + tmp_5*(tmp_25*tmp_43 + tmp_32*tmp_44 + tmp_39*tmp_45 - 1.0/4.0);
      real_t tmp_47 = -p_affine_4_1;
      real_t tmp_48 = p_affine_5_1 + tmp_47;
      real_t tmp_49 = -p_affine_4_2;
      real_t tmp_50 = p_affine_6_2 + tmp_49;
      real_t tmp_51 = tmp_48*tmp_50;
      real_t tmp_52 = p_affine_6_1 + tmp_47;
      real_t tmp_53 = p_affine_5_2 + tmp_49;
      real_t tmp_54 = p_affine_7_2 + tmp_49;
      real_t tmp_55 = -p_affine_4_0;
      real_t tmp_56 = p_affine_5_0 + tmp_55;
      real_t tmp_57 = tmp_52*tmp_56;
      real_t tmp_58 = p_affine_7_1 + tmp_47;
      real_t tmp_59 = p_affine_6_0 + tmp_55;
      real_t tmp_60 = tmp_53*tmp_59;
      real_t tmp_61 = p_affine_7_0 + tmp_55;
      real_t tmp_62 = tmp_56*tmp_58;
      real_t tmp_63 = tmp_48*tmp_59;
      real_t tmp_64 = tmp_53*tmp_61;
      real_t tmp_65 = 1.0 / (-tmp_50*tmp_62 + tmp_51*tmp_61 - tmp_52*tmp_64 + tmp_54*tmp_57 - tmp_54*tmp_63 + tmp_58*tmp_60);
      real_t tmp_66 = tmp_65*(tmp_51 - tmp_52*tmp_53);
      real_t tmp_67 = tmp_65*(-tmp_48*tmp_54 + tmp_53*tmp_58);
      real_t tmp_68 = tmp_65*(-tmp_50*tmp_58 + tmp_52*tmp_54);
      real_t tmp_69 = tmp_65*(-tmp_50*tmp_56 + tmp_60);
      real_t tmp_70 = tmp_65*(tmp_54*tmp_56 - tmp_64);
      real_t tmp_71 = tmp_65*(tmp_50*tmp_61 - tmp_54*tmp_59);
      real_t tmp_72 = tmp_65*(tmp_57 - tmp_63);
      real_t tmp_73 = tmp_65*(tmp_48*tmp_61 - tmp_62);
      real_t tmp_74 = tmp_65*(-tmp_52*tmp_61 + tmp_58*tmp_59);
      real_t tmp_75 = 0.5*p_affine_13_0*(-tmp_66 - tmp_67 - tmp_68) + 0.5*p_affine_13_1*(-tmp_69 - tmp_70 - tmp_71) + 0.5*p_affine_13_2*(-tmp_72 - tmp_73 - tmp_74);
      real_t tmp_76 = p_affine_8_2 + tmp_49;
      real_t tmp_77 = tmp_23 + tmp_76;
      real_t tmp_78 = tmp_72*tmp_77;
      real_t tmp_79 = tmp_73*tmp_77;
      real_t tmp_80 = p_affine_8_1 + tmp_47;
      real_t tmp_81 = tmp_30 + tmp_80;
      real_t tmp_82 = tmp_69*tmp_81;
      real_t tmp_83 = tmp_70*tmp_81;
      real_t tmp_84 = tmp_74*tmp_77;
      real_t tmp_85 = tmp_71*tmp_81;
      real_t tmp_86 = p_affine_8_0 + tmp_55;
      real_t tmp_87 = tmp_37 + tmp_86;
      real_t tmp_88 = tmp_66*tmp_87;
      real_t tmp_89 = tmp_67*tmp_87;
      real_t tmp_90 = tmp_68*tmp_87;
      real_t tmp_91 = -tmp_78 - tmp_79 - tmp_82 - tmp_83 - tmp_84 - tmp_85 - tmp_88 - tmp_89 - tmp_90 + 1;
      real_t tmp_92 = tmp_1*tmp_19;
      real_t tmp_93 = tmp_19*tmp_2;
      real_t tmp_94 = tmp_19*tmp_5;
      real_t tmp_95 = 0.5*p_affine_13_0*(tmp_33*tmp_92 + tmp_42*tmp_93 + tmp_45*tmp_94) + 0.5*p_affine_13_1*(tmp_26*tmp_92 + tmp_41*tmp_93 + tmp_44*tmp_94) + 0.5*p_affine_13_2*(tmp_40*tmp_93 + tmp_43*tmp_94 + tmp_7*tmp_92);
      real_t tmp_96 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_97 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_98 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_99 = (std::abs(tmp_22*tmp_96 - tmp_29*tmp_98)*std::abs(tmp_22*tmp_96 - tmp_29*tmp_98)) + (std::abs(tmp_22*tmp_97 - tmp_36*tmp_98)*std::abs(tmp_22*tmp_97 - tmp_36*tmp_98)) + (std::abs(tmp_29*tmp_97 - tmp_36*tmp_96)*std::abs(tmp_29*tmp_97 - tmp_36*tmp_96));
      real_t tmp_100 = 5.0*std::pow(tmp_99, -0.25);
      real_t tmp_101 = tmp_100*tmp_46;
      real_t tmp_102 = 1.0*std::pow(tmp_99, 1.0/2.0);
      real_t tmp_103 = 0.0068572537431980923*tmp_102;
      real_t tmp_104 = 0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22;
      real_t tmp_105 = tmp_19*(tmp_104 + tmp_24);
      real_t tmp_106 = 0.19601935860219369*tmp_28 + 0.60796128279561268*tmp_29;
      real_t tmp_107 = tmp_19*(tmp_106 + tmp_31);
      real_t tmp_108 = 0.19601935860219369*tmp_35 + 0.60796128279561268*tmp_36;
      real_t tmp_109 = tmp_19*(tmp_108 + tmp_38);
      real_t tmp_110 = tmp_1*(tmp_105*tmp_7 + tmp_107*tmp_26 + tmp_109*tmp_33 - 1.0/4.0) + tmp_2*(tmp_105*tmp_40 + tmp_107*tmp_41 + tmp_109*tmp_42 - 1.0/4.0) + tmp_5*(tmp_105*tmp_43 + tmp_107*tmp_44 + tmp_109*tmp_45 - 1.0/4.0);
      real_t tmp_111 = tmp_104 + tmp_76;
      real_t tmp_112 = tmp_111*tmp_72;
      real_t tmp_113 = tmp_111*tmp_73;
      real_t tmp_114 = tmp_106 + tmp_80;
      real_t tmp_115 = tmp_114*tmp_69;
      real_t tmp_116 = tmp_114*tmp_70;
      real_t tmp_117 = tmp_111*tmp_74;
      real_t tmp_118 = tmp_114*tmp_71;
      real_t tmp_119 = tmp_108 + tmp_86;
      real_t tmp_120 = tmp_119*tmp_66;
      real_t tmp_121 = tmp_119*tmp_67;
      real_t tmp_122 = tmp_119*tmp_68;
      real_t tmp_123 = -tmp_112 - tmp_113 - tmp_115 - tmp_116 - tmp_117 - tmp_118 - tmp_120 - tmp_121 - tmp_122 + 1;
      real_t tmp_124 = tmp_100*tmp_110;
      real_t tmp_125 = 0.037198804536718075*tmp_102;
      real_t tmp_126 = 0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_127 = tmp_19*(tmp_126 + tmp_24);
      real_t tmp_128 = 0.37605877282253791*tmp_28 + 0.039308471900058539*tmp_29;
      real_t tmp_129 = tmp_19*(tmp_128 + tmp_31);
      real_t tmp_130 = 0.37605877282253791*tmp_35 + 0.039308471900058539*tmp_36;
      real_t tmp_131 = tmp_19*(tmp_130 + tmp_38);
      real_t tmp_132 = tmp_1*(tmp_127*tmp_7 + tmp_129*tmp_26 + tmp_131*tmp_33 - 1.0/4.0) + tmp_2*(tmp_127*tmp_40 + tmp_129*tmp_41 + tmp_131*tmp_42 - 1.0/4.0) + tmp_5*(tmp_127*tmp_43 + tmp_129*tmp_44 + tmp_131*tmp_45 - 1.0/4.0);
      real_t tmp_133 = tmp_126 + tmp_76;
      real_t tmp_134 = tmp_133*tmp_72;
      real_t tmp_135 = tmp_133*tmp_73;
      real_t tmp_136 = tmp_128 + tmp_80;
      real_t tmp_137 = tmp_136*tmp_69;
      real_t tmp_138 = tmp_136*tmp_70;
      real_t tmp_139 = tmp_133*tmp_74;
      real_t tmp_140 = tmp_136*tmp_71;
      real_t tmp_141 = tmp_130 + tmp_86;
      real_t tmp_142 = tmp_141*tmp_66;
      real_t tmp_143 = tmp_141*tmp_67;
      real_t tmp_144 = tmp_141*tmp_68;
      real_t tmp_145 = -tmp_134 - tmp_135 - tmp_137 - tmp_138 - tmp_139 - tmp_140 - tmp_142 - tmp_143 - tmp_144 + 1;
      real_t tmp_146 = tmp_100*tmp_132;
      real_t tmp_147 = 0.020848748529055869*tmp_102;
      real_t tmp_148 = 0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_149 = tmp_19*(tmp_148 + tmp_24);
      real_t tmp_150 = 0.78764240869137092*tmp_28 + 0.1711304259088916*tmp_29;
      real_t tmp_151 = tmp_19*(tmp_150 + tmp_31);
      real_t tmp_152 = 0.78764240869137092*tmp_35 + 0.1711304259088916*tmp_36;
      real_t tmp_153 = tmp_19*(tmp_152 + tmp_38);
      real_t tmp_154 = tmp_1*(tmp_149*tmp_7 + tmp_151*tmp_26 + tmp_153*tmp_33 - 1.0/4.0) + tmp_2*(tmp_149*tmp_40 + tmp_151*tmp_41 + tmp_153*tmp_42 - 1.0/4.0) + tmp_5*(tmp_149*tmp_43 + tmp_151*tmp_44 + tmp_153*tmp_45 - 1.0/4.0);
      real_t tmp_155 = tmp_148 + tmp_76;
      real_t tmp_156 = tmp_155*tmp_72;
      real_t tmp_157 = tmp_155*tmp_73;
      real_t tmp_158 = tmp_150 + tmp_80;
      real_t tmp_159 = tmp_158*tmp_69;
      real_t tmp_160 = tmp_158*tmp_70;
      real_t tmp_161 = tmp_155*tmp_74;
      real_t tmp_162 = tmp_158*tmp_71;
      real_t tmp_163 = tmp_152 + tmp_86;
      real_t tmp_164 = tmp_163*tmp_66;
      real_t tmp_165 = tmp_163*tmp_67;
      real_t tmp_166 = tmp_163*tmp_68;
      real_t tmp_167 = -tmp_156 - tmp_157 - tmp_159 - tmp_160 - tmp_161 - tmp_162 - tmp_164 - tmp_165 - tmp_166 + 1;
      real_t tmp_168 = tmp_100*tmp_154;
      real_t tmp_169 = 0.019202922745021479*tmp_102;
      real_t tmp_170 = 0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_171 = tmp_19*(tmp_170 + tmp_24);
      real_t tmp_172 = 0.58463275527740355*tmp_28 + 0.37605877282253791*tmp_29;
      real_t tmp_173 = tmp_19*(tmp_172 + tmp_31);
      real_t tmp_174 = 0.58463275527740355*tmp_35 + 0.37605877282253791*tmp_36;
      real_t tmp_175 = tmp_19*(tmp_174 + tmp_38);
      real_t tmp_176 = tmp_1*(tmp_171*tmp_7 + tmp_173*tmp_26 + tmp_175*tmp_33 - 1.0/4.0) + tmp_2*(tmp_171*tmp_40 + tmp_173*tmp_41 + tmp_175*tmp_42 - 1.0/4.0) + tmp_5*(tmp_171*tmp_43 + tmp_173*tmp_44 + tmp_175*tmp_45 - 1.0/4.0);
      real_t tmp_177 = tmp_170 + tmp_76;
      real_t tmp_178 = tmp_177*tmp_72;
      real_t tmp_179 = tmp_177*tmp_73;
      real_t tmp_180 = tmp_172 + tmp_80;
      real_t tmp_181 = tmp_180*tmp_69;
      real_t tmp_182 = tmp_180*tmp_70;
      real_t tmp_183 = tmp_177*tmp_74;
      real_t tmp_184 = tmp_180*tmp_71;
      real_t tmp_185 = tmp_174 + tmp_86;
      real_t tmp_186 = tmp_185*tmp_66;
      real_t tmp_187 = tmp_185*tmp_67;
      real_t tmp_188 = tmp_185*tmp_68;
      real_t tmp_189 = -tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 - tmp_184 - tmp_186 - tmp_187 - tmp_188 + 1;
      real_t tmp_190 = tmp_100*tmp_176;
      real_t tmp_191 = 0.020848748529055869*tmp_102;
      real_t tmp_192 = 0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_193 = tmp_19*(tmp_192 + tmp_24);
      real_t tmp_194 = 0.041227165399737475*tmp_28 + 0.78764240869137092*tmp_29;
      real_t tmp_195 = tmp_19*(tmp_194 + tmp_31);
      real_t tmp_196 = 0.041227165399737475*tmp_35 + 0.78764240869137092*tmp_36;
      real_t tmp_197 = tmp_19*(tmp_196 + tmp_38);
      real_t tmp_198 = tmp_1*(tmp_193*tmp_7 + tmp_195*tmp_26 + tmp_197*tmp_33 - 1.0/4.0) + tmp_2*(tmp_193*tmp_40 + tmp_195*tmp_41 + tmp_197*tmp_42 - 1.0/4.0) + tmp_5*(tmp_193*tmp_43 + tmp_195*tmp_44 + tmp_197*tmp_45 - 1.0/4.0);
      real_t tmp_199 = tmp_192 + tmp_76;
      real_t tmp_200 = tmp_199*tmp_72;
      real_t tmp_201 = tmp_199*tmp_73;
      real_t tmp_202 = tmp_194 + tmp_80;
      real_t tmp_203 = tmp_202*tmp_69;
      real_t tmp_204 = tmp_202*tmp_70;
      real_t tmp_205 = tmp_199*tmp_74;
      real_t tmp_206 = tmp_202*tmp_71;
      real_t tmp_207 = tmp_196 + tmp_86;
      real_t tmp_208 = tmp_207*tmp_66;
      real_t tmp_209 = tmp_207*tmp_67;
      real_t tmp_210 = tmp_207*tmp_68;
      real_t tmp_211 = -tmp_200 - tmp_201 - tmp_203 - tmp_204 - tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 + 1;
      real_t tmp_212 = tmp_100*tmp_198;
      real_t tmp_213 = 0.019202922745021479*tmp_102;
      real_t tmp_214 = 0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_215 = tmp_19*(tmp_214 + tmp_24);
      real_t tmp_216 = 0.039308471900058539*tmp_28 + 0.58463275527740355*tmp_29;
      real_t tmp_217 = tmp_19*(tmp_216 + tmp_31);
      real_t tmp_218 = 0.039308471900058539*tmp_35 + 0.58463275527740355*tmp_36;
      real_t tmp_219 = tmp_19*(tmp_218 + tmp_38);
      real_t tmp_220 = tmp_1*(tmp_215*tmp_7 + tmp_217*tmp_26 + tmp_219*tmp_33 - 1.0/4.0) + tmp_2*(tmp_215*tmp_40 + tmp_217*tmp_41 + tmp_219*tmp_42 - 1.0/4.0) + tmp_5*(tmp_215*tmp_43 + tmp_217*tmp_44 + tmp_219*tmp_45 - 1.0/4.0);
      real_t tmp_221 = tmp_214 + tmp_76;
      real_t tmp_222 = tmp_221*tmp_72;
      real_t tmp_223 = tmp_221*tmp_73;
      real_t tmp_224 = tmp_216 + tmp_80;
      real_t tmp_225 = tmp_224*tmp_69;
      real_t tmp_226 = tmp_224*tmp_70;
      real_t tmp_227 = tmp_221*tmp_74;
      real_t tmp_228 = tmp_224*tmp_71;
      real_t tmp_229 = tmp_218 + tmp_86;
      real_t tmp_230 = tmp_229*tmp_66;
      real_t tmp_231 = tmp_229*tmp_67;
      real_t tmp_232 = tmp_229*tmp_68;
      real_t tmp_233 = -tmp_222 - tmp_223 - tmp_225 - tmp_226 - tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 + 1;
      real_t tmp_234 = tmp_100*tmp_220;
      real_t tmp_235 = 0.020848748529055869*tmp_102;
      real_t tmp_236 = 0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_237 = tmp_19*(tmp_236 + tmp_24);
      real_t tmp_238 = 0.78764240869137092*tmp_28 + 0.041227165399737475*tmp_29;
      real_t tmp_239 = tmp_19*(tmp_238 + tmp_31);
      real_t tmp_240 = 0.78764240869137092*tmp_35 + 0.041227165399737475*tmp_36;
      real_t tmp_241 = tmp_19*(tmp_240 + tmp_38);
      real_t tmp_242 = tmp_1*(tmp_237*tmp_7 + tmp_239*tmp_26 + tmp_241*tmp_33 - 1.0/4.0) + tmp_2*(tmp_237*tmp_40 + tmp_239*tmp_41 + tmp_241*tmp_42 - 1.0/4.0) + tmp_5*(tmp_237*tmp_43 + tmp_239*tmp_44 + tmp_241*tmp_45 - 1.0/4.0);
      real_t tmp_243 = tmp_236 + tmp_76;
      real_t tmp_244 = tmp_243*tmp_72;
      real_t tmp_245 = tmp_243*tmp_73;
      real_t tmp_246 = tmp_238 + tmp_80;
      real_t tmp_247 = tmp_246*tmp_69;
      real_t tmp_248 = tmp_246*tmp_70;
      real_t tmp_249 = tmp_243*tmp_74;
      real_t tmp_250 = tmp_246*tmp_71;
      real_t tmp_251 = tmp_240 + tmp_86;
      real_t tmp_252 = tmp_251*tmp_66;
      real_t tmp_253 = tmp_251*tmp_67;
      real_t tmp_254 = tmp_251*tmp_68;
      real_t tmp_255 = -tmp_244 - tmp_245 - tmp_247 - tmp_248 - tmp_249 - tmp_250 - tmp_252 - tmp_253 - tmp_254 + 1;
      real_t tmp_256 = tmp_100*tmp_242;
      real_t tmp_257 = 0.019202922745021479*tmp_102;
      real_t tmp_258 = 0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_259 = tmp_19*(tmp_24 + tmp_258);
      real_t tmp_260 = 0.58463275527740355*tmp_28 + 0.039308471900058539*tmp_29;
      real_t tmp_261 = tmp_19*(tmp_260 + tmp_31);
      real_t tmp_262 = 0.58463275527740355*tmp_35 + 0.039308471900058539*tmp_36;
      real_t tmp_263 = tmp_19*(tmp_262 + tmp_38);
      real_t tmp_264 = tmp_1*(tmp_259*tmp_7 + tmp_26*tmp_261 + tmp_263*tmp_33 - 1.0/4.0) + tmp_2*(tmp_259*tmp_40 + tmp_261*tmp_41 + tmp_263*tmp_42 - 1.0/4.0) + tmp_5*(tmp_259*tmp_43 + tmp_261*tmp_44 + tmp_263*tmp_45 - 1.0/4.0);
      real_t tmp_265 = tmp_258 + tmp_76;
      real_t tmp_266 = tmp_265*tmp_72;
      real_t tmp_267 = tmp_265*tmp_73;
      real_t tmp_268 = tmp_260 + tmp_80;
      real_t tmp_269 = tmp_268*tmp_69;
      real_t tmp_270 = tmp_268*tmp_70;
      real_t tmp_271 = tmp_265*tmp_74;
      real_t tmp_272 = tmp_268*tmp_71;
      real_t tmp_273 = tmp_262 + tmp_86;
      real_t tmp_274 = tmp_273*tmp_66;
      real_t tmp_275 = tmp_273*tmp_67;
      real_t tmp_276 = tmp_273*tmp_68;
      real_t tmp_277 = -tmp_266 - tmp_267 - tmp_269 - tmp_270 - tmp_271 - tmp_272 - tmp_274 - tmp_275 - tmp_276 + 1;
      real_t tmp_278 = tmp_100*tmp_264;
      real_t tmp_279 = 0.020848748529055869*tmp_102;
      real_t tmp_280 = 0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_281 = tmp_19*(tmp_24 + tmp_280);
      real_t tmp_282 = 0.1711304259088916*tmp_28 + 0.78764240869137092*tmp_29;
      real_t tmp_283 = tmp_19*(tmp_282 + tmp_31);
      real_t tmp_284 = 0.1711304259088916*tmp_35 + 0.78764240869137092*tmp_36;
      real_t tmp_285 = tmp_19*(tmp_284 + tmp_38);
      real_t tmp_286 = tmp_1*(tmp_26*tmp_283 + tmp_281*tmp_7 + tmp_285*tmp_33 - 1.0/4.0) + tmp_2*(tmp_281*tmp_40 + tmp_283*tmp_41 + tmp_285*tmp_42 - 1.0/4.0) + tmp_5*(tmp_281*tmp_43 + tmp_283*tmp_44 + tmp_285*tmp_45 - 1.0/4.0);
      real_t tmp_287 = tmp_280 + tmp_76;
      real_t tmp_288 = tmp_287*tmp_72;
      real_t tmp_289 = tmp_287*tmp_73;
      real_t tmp_290 = tmp_282 + tmp_80;
      real_t tmp_291 = tmp_290*tmp_69;
      real_t tmp_292 = tmp_290*tmp_70;
      real_t tmp_293 = tmp_287*tmp_74;
      real_t tmp_294 = tmp_290*tmp_71;
      real_t tmp_295 = tmp_284 + tmp_86;
      real_t tmp_296 = tmp_295*tmp_66;
      real_t tmp_297 = tmp_295*tmp_67;
      real_t tmp_298 = tmp_295*tmp_68;
      real_t tmp_299 = -tmp_288 - tmp_289 - tmp_291 - tmp_292 - tmp_293 - tmp_294 - tmp_296 - tmp_297 - tmp_298 + 1;
      real_t tmp_300 = tmp_100*tmp_286;
      real_t tmp_301 = 0.019202922745021479*tmp_102;
      real_t tmp_302 = 0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_303 = tmp_19*(tmp_24 + tmp_302);
      real_t tmp_304 = 0.37605877282253791*tmp_28 + 0.58463275527740355*tmp_29;
      real_t tmp_305 = tmp_19*(tmp_304 + tmp_31);
      real_t tmp_306 = 0.37605877282253791*tmp_35 + 0.58463275527740355*tmp_36;
      real_t tmp_307 = tmp_19*(tmp_306 + tmp_38);
      real_t tmp_308 = tmp_1*(tmp_26*tmp_305 + tmp_303*tmp_7 + tmp_307*tmp_33 - 1.0/4.0) + tmp_2*(tmp_303*tmp_40 + tmp_305*tmp_41 + tmp_307*tmp_42 - 1.0/4.0) + tmp_5*(tmp_303*tmp_43 + tmp_305*tmp_44 + tmp_307*tmp_45 - 1.0/4.0);
      real_t tmp_309 = tmp_302 + tmp_76;
      real_t tmp_310 = tmp_309*tmp_72;
      real_t tmp_311 = tmp_309*tmp_73;
      real_t tmp_312 = tmp_304 + tmp_80;
      real_t tmp_313 = tmp_312*tmp_69;
      real_t tmp_314 = tmp_312*tmp_70;
      real_t tmp_315 = tmp_309*tmp_74;
      real_t tmp_316 = tmp_312*tmp_71;
      real_t tmp_317 = tmp_306 + tmp_86;
      real_t tmp_318 = tmp_317*tmp_66;
      real_t tmp_319 = tmp_317*tmp_67;
      real_t tmp_320 = tmp_317*tmp_68;
      real_t tmp_321 = -tmp_310 - tmp_311 - tmp_313 - tmp_314 - tmp_315 - tmp_316 - tmp_318 - tmp_319 - tmp_320 + 1;
      real_t tmp_322 = tmp_100*tmp_308;
      real_t tmp_323 = 0.020848748529055869*tmp_102;
      real_t tmp_324 = 0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_325 = tmp_19*(tmp_24 + tmp_324);
      real_t tmp_326 = 0.041227165399737475*tmp_28 + 0.1711304259088916*tmp_29;
      real_t tmp_327 = tmp_19*(tmp_31 + tmp_326);
      real_t tmp_328 = 0.041227165399737475*tmp_35 + 0.1711304259088916*tmp_36;
      real_t tmp_329 = tmp_19*(tmp_328 + tmp_38);
      real_t tmp_330 = tmp_1*(tmp_26*tmp_327 + tmp_325*tmp_7 + tmp_329*tmp_33 - 1.0/4.0) + tmp_2*(tmp_325*tmp_40 + tmp_327*tmp_41 + tmp_329*tmp_42 - 1.0/4.0) + tmp_5*(tmp_325*tmp_43 + tmp_327*tmp_44 + tmp_329*tmp_45 - 1.0/4.0);
      real_t tmp_331 = tmp_324 + tmp_76;
      real_t tmp_332 = tmp_331*tmp_72;
      real_t tmp_333 = tmp_331*tmp_73;
      real_t tmp_334 = tmp_326 + tmp_80;
      real_t tmp_335 = tmp_334*tmp_69;
      real_t tmp_336 = tmp_334*tmp_70;
      real_t tmp_337 = tmp_331*tmp_74;
      real_t tmp_338 = tmp_334*tmp_71;
      real_t tmp_339 = tmp_328 + tmp_86;
      real_t tmp_340 = tmp_339*tmp_66;
      real_t tmp_341 = tmp_339*tmp_67;
      real_t tmp_342 = tmp_339*tmp_68;
      real_t tmp_343 = -tmp_332 - tmp_333 - tmp_335 - tmp_336 - tmp_337 - tmp_338 - tmp_340 - tmp_341 - tmp_342 + 1;
      real_t tmp_344 = tmp_100*tmp_330;
      real_t tmp_345 = 0.019202922745021479*tmp_102;
      real_t tmp_346 = 0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22;
      real_t tmp_347 = tmp_19*(tmp_24 + tmp_346);
      real_t tmp_348 = 0.40446199974765351*tmp_28 + 0.19107600050469298*tmp_29;
      real_t tmp_349 = tmp_19*(tmp_31 + tmp_348);
      real_t tmp_350 = 0.40446199974765351*tmp_35 + 0.19107600050469298*tmp_36;
      real_t tmp_351 = tmp_19*(tmp_350 + tmp_38);
      real_t tmp_352 = tmp_1*(tmp_26*tmp_349 + tmp_33*tmp_351 + tmp_347*tmp_7 - 1.0/4.0) + tmp_2*(tmp_347*tmp_40 + tmp_349*tmp_41 + tmp_351*tmp_42 - 1.0/4.0) + tmp_5*(tmp_347*tmp_43 + tmp_349*tmp_44 + tmp_351*tmp_45 - 1.0/4.0);
      real_t tmp_353 = tmp_346 + tmp_76;
      real_t tmp_354 = tmp_353*tmp_72;
      real_t tmp_355 = tmp_353*tmp_73;
      real_t tmp_356 = tmp_348 + tmp_80;
      real_t tmp_357 = tmp_356*tmp_69;
      real_t tmp_358 = tmp_356*tmp_70;
      real_t tmp_359 = tmp_353*tmp_74;
      real_t tmp_360 = tmp_356*tmp_71;
      real_t tmp_361 = tmp_350 + tmp_86;
      real_t tmp_362 = tmp_361*tmp_66;
      real_t tmp_363 = tmp_361*tmp_67;
      real_t tmp_364 = tmp_361*tmp_68;
      real_t tmp_365 = -tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 - tmp_360 - tmp_362 - tmp_363 - tmp_364 + 1;
      real_t tmp_366 = tmp_100*tmp_352;
      real_t tmp_367 = 0.042507265838595799*tmp_102;
      real_t tmp_368 = 0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_369 = tmp_19*(tmp_24 + tmp_368);
      real_t tmp_370 = 0.039308471900058539*tmp_28 + 0.37605877282253791*tmp_29;
      real_t tmp_371 = tmp_19*(tmp_31 + tmp_370);
      real_t tmp_372 = 0.039308471900058539*tmp_35 + 0.37605877282253791*tmp_36;
      real_t tmp_373 = tmp_19*(tmp_372 + tmp_38);
      real_t tmp_374 = tmp_1*(tmp_26*tmp_371 + tmp_33*tmp_373 + tmp_369*tmp_7 - 1.0/4.0) + tmp_2*(tmp_369*tmp_40 + tmp_371*tmp_41 + tmp_373*tmp_42 - 1.0/4.0) + tmp_5*(tmp_369*tmp_43 + tmp_371*tmp_44 + tmp_373*tmp_45 - 1.0/4.0);
      real_t tmp_375 = tmp_368 + tmp_76;
      real_t tmp_376 = tmp_375*tmp_72;
      real_t tmp_377 = tmp_375*tmp_73;
      real_t tmp_378 = tmp_370 + tmp_80;
      real_t tmp_379 = tmp_378*tmp_69;
      real_t tmp_380 = tmp_378*tmp_70;
      real_t tmp_381 = tmp_375*tmp_74;
      real_t tmp_382 = tmp_378*tmp_71;
      real_t tmp_383 = tmp_372 + tmp_86;
      real_t tmp_384 = tmp_383*tmp_66;
      real_t tmp_385 = tmp_383*tmp_67;
      real_t tmp_386 = tmp_383*tmp_68;
      real_t tmp_387 = -tmp_376 - tmp_377 - tmp_379 - tmp_380 - tmp_381 - tmp_382 - tmp_384 - tmp_385 - tmp_386 + 1;
      real_t tmp_388 = tmp_100*tmp_374;
      real_t tmp_389 = 0.020848748529055869*tmp_102;
      real_t tmp_390 = 0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_391 = tmp_19*(tmp_24 + tmp_390);
      real_t tmp_392 = 0.93718850182767688*tmp_28 + 0.031405749086161582*tmp_29;
      real_t tmp_393 = tmp_19*(tmp_31 + tmp_392);
      real_t tmp_394 = 0.93718850182767688*tmp_35 + 0.031405749086161582*tmp_36;
      real_t tmp_395 = tmp_19*(tmp_38 + tmp_394);
      real_t tmp_396 = tmp_1*(tmp_26*tmp_393 + tmp_33*tmp_395 + tmp_391*tmp_7 - 1.0/4.0) + tmp_2*(tmp_391*tmp_40 + tmp_393*tmp_41 + tmp_395*tmp_42 - 1.0/4.0) + tmp_5*(tmp_391*tmp_43 + tmp_393*tmp_44 + tmp_395*tmp_45 - 1.0/4.0);
      real_t tmp_397 = tmp_390 + tmp_76;
      real_t tmp_398 = tmp_397*tmp_72;
      real_t tmp_399 = tmp_397*tmp_73;
      real_t tmp_400 = tmp_392 + tmp_80;
      real_t tmp_401 = tmp_400*tmp_69;
      real_t tmp_402 = tmp_400*tmp_70;
      real_t tmp_403 = tmp_397*tmp_74;
      real_t tmp_404 = tmp_400*tmp_71;
      real_t tmp_405 = tmp_394 + tmp_86;
      real_t tmp_406 = tmp_405*tmp_66;
      real_t tmp_407 = tmp_405*tmp_67;
      real_t tmp_408 = tmp_405*tmp_68;
      real_t tmp_409 = -tmp_398 - tmp_399 - tmp_401 - tmp_402 - tmp_403 - tmp_404 - tmp_406 - tmp_407 - tmp_408 + 1;
      real_t tmp_410 = tmp_100*tmp_396;
      real_t tmp_411 = 0.0068572537431980923*tmp_102;
      real_t tmp_412 = 0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_413 = tmp_19*(tmp_24 + tmp_412);
      real_t tmp_414 = 0.60796128279561268*tmp_28 + 0.19601935860219369*tmp_29;
      real_t tmp_415 = tmp_19*(tmp_31 + tmp_414);
      real_t tmp_416 = 0.60796128279561268*tmp_35 + 0.19601935860219369*tmp_36;
      real_t tmp_417 = tmp_19*(tmp_38 + tmp_416);
      real_t tmp_418 = tmp_1*(tmp_26*tmp_415 + tmp_33*tmp_417 + tmp_413*tmp_7 - 1.0/4.0) + tmp_2*(tmp_40*tmp_413 + tmp_41*tmp_415 + tmp_417*tmp_42 - 1.0/4.0) + tmp_5*(tmp_413*tmp_43 + tmp_415*tmp_44 + tmp_417*tmp_45 - 1.0/4.0);
      real_t tmp_419 = tmp_412 + tmp_76;
      real_t tmp_420 = tmp_419*tmp_72;
      real_t tmp_421 = tmp_419*tmp_73;
      real_t tmp_422 = tmp_414 + tmp_80;
      real_t tmp_423 = tmp_422*tmp_69;
      real_t tmp_424 = tmp_422*tmp_70;
      real_t tmp_425 = tmp_419*tmp_74;
      real_t tmp_426 = tmp_422*tmp_71;
      real_t tmp_427 = tmp_416 + tmp_86;
      real_t tmp_428 = tmp_427*tmp_66;
      real_t tmp_429 = tmp_427*tmp_67;
      real_t tmp_430 = tmp_427*tmp_68;
      real_t tmp_431 = -tmp_420 - tmp_421 - tmp_423 - tmp_424 - tmp_425 - tmp_426 - tmp_428 - tmp_429 - tmp_430 + 1;
      real_t tmp_432 = tmp_100*tmp_418;
      real_t tmp_433 = 0.037198804536718075*tmp_102;
      real_t tmp_434 = 0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_435 = tmp_19*(tmp_24 + tmp_434);
      real_t tmp_436 = 0.19107600050469298*tmp_28 + 0.40446199974765351*tmp_29;
      real_t tmp_437 = tmp_19*(tmp_31 + tmp_436);
      real_t tmp_438 = 0.19107600050469298*tmp_35 + 0.40446199974765351*tmp_36;
      real_t tmp_439 = tmp_19*(tmp_38 + tmp_438);
      real_t tmp_440 = tmp_1*(tmp_26*tmp_437 + tmp_33*tmp_439 + tmp_435*tmp_7 - 1.0/4.0) + tmp_2*(tmp_40*tmp_435 + tmp_41*tmp_437 + tmp_42*tmp_439 - 1.0/4.0) + tmp_5*(tmp_43*tmp_435 + tmp_437*tmp_44 + tmp_439*tmp_45 - 1.0/4.0);
      real_t tmp_441 = tmp_434 + tmp_76;
      real_t tmp_442 = tmp_441*tmp_72;
      real_t tmp_443 = tmp_441*tmp_73;
      real_t tmp_444 = tmp_436 + tmp_80;
      real_t tmp_445 = tmp_444*tmp_69;
      real_t tmp_446 = tmp_444*tmp_70;
      real_t tmp_447 = tmp_441*tmp_74;
      real_t tmp_448 = tmp_444*tmp_71;
      real_t tmp_449 = tmp_438 + tmp_86;
      real_t tmp_450 = tmp_449*tmp_66;
      real_t tmp_451 = tmp_449*tmp_67;
      real_t tmp_452 = tmp_449*tmp_68;
      real_t tmp_453 = -tmp_442 - tmp_443 - tmp_445 - tmp_446 - tmp_447 - tmp_448 - tmp_450 - tmp_451 - tmp_452 + 1;
      real_t tmp_454 = tmp_100*tmp_440;
      real_t tmp_455 = 0.042507265838595799*tmp_102;
      real_t tmp_456 = 0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_457 = tmp_19*(tmp_24 + tmp_456);
      real_t tmp_458 = 0.031405749086161582*tmp_28 + 0.031405749086161582*tmp_29;
      real_t tmp_459 = tmp_19*(tmp_31 + tmp_458);
      real_t tmp_460 = 0.031405749086161582*tmp_35 + 0.031405749086161582*tmp_36;
      real_t tmp_461 = tmp_19*(tmp_38 + tmp_460);
      real_t tmp_462 = tmp_1*(tmp_26*tmp_459 + tmp_33*tmp_461 + tmp_457*tmp_7 - 1.0/4.0) + tmp_2*(tmp_40*tmp_457 + tmp_41*tmp_459 + tmp_42*tmp_461 - 1.0/4.0) + tmp_5*(tmp_43*tmp_457 + tmp_44*tmp_459 + tmp_45*tmp_461 - 1.0/4.0);
      real_t tmp_463 = tmp_456 + tmp_76;
      real_t tmp_464 = tmp_463*tmp_72;
      real_t tmp_465 = tmp_463*tmp_73;
      real_t tmp_466 = tmp_458 + tmp_80;
      real_t tmp_467 = tmp_466*tmp_69;
      real_t tmp_468 = tmp_466*tmp_70;
      real_t tmp_469 = tmp_463*tmp_74;
      real_t tmp_470 = tmp_466*tmp_71;
      real_t tmp_471 = tmp_460 + tmp_86;
      real_t tmp_472 = tmp_471*tmp_66;
      real_t tmp_473 = tmp_471*tmp_67;
      real_t tmp_474 = tmp_471*tmp_68;
      real_t tmp_475 = -tmp_464 - tmp_465 - tmp_467 - tmp_468 - tmp_469 - tmp_470 - tmp_472 - tmp_473 - tmp_474 + 1;
      real_t tmp_476 = tmp_100*tmp_462;
      real_t tmp_477 = 0.0068572537431980923*tmp_102;
      real_t tmp_478 = 0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_479 = tmp_19*(tmp_24 + tmp_478);
      real_t tmp_480 = 0.19601935860219369*tmp_28 + 0.19601935860219369*tmp_29;
      real_t tmp_481 = tmp_19*(tmp_31 + tmp_480);
      real_t tmp_482 = 0.19601935860219369*tmp_35 + 0.19601935860219369*tmp_36;
      real_t tmp_483 = tmp_19*(tmp_38 + tmp_482);
      real_t tmp_484 = tmp_1*(tmp_26*tmp_481 + tmp_33*tmp_483 + tmp_479*tmp_7 - 1.0/4.0) + tmp_2*(tmp_40*tmp_479 + tmp_41*tmp_481 + tmp_42*tmp_483 - 1.0/4.0) + tmp_5*(tmp_43*tmp_479 + tmp_44*tmp_481 + tmp_45*tmp_483 - 1.0/4.0);
      real_t tmp_485 = tmp_478 + tmp_76;
      real_t tmp_486 = tmp_485*tmp_72;
      real_t tmp_487 = tmp_485*tmp_73;
      real_t tmp_488 = tmp_480 + tmp_80;
      real_t tmp_489 = tmp_488*tmp_69;
      real_t tmp_490 = tmp_488*tmp_70;
      real_t tmp_491 = tmp_485*tmp_74;
      real_t tmp_492 = tmp_488*tmp_71;
      real_t tmp_493 = tmp_482 + tmp_86;
      real_t tmp_494 = tmp_493*tmp_66;
      real_t tmp_495 = tmp_493*tmp_67;
      real_t tmp_496 = tmp_493*tmp_68;
      real_t tmp_497 = -tmp_486 - tmp_487 - tmp_489 - tmp_490 - tmp_491 - tmp_492 - tmp_494 - tmp_495 - tmp_496 + 1;
      real_t tmp_498 = tmp_100*tmp_484;
      real_t tmp_499 = 0.037198804536718075*tmp_102;
      real_t tmp_500 = 0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_501 = tmp_19*(tmp_24 + tmp_500);
      real_t tmp_502 = 0.40446199974765351*tmp_28 + 0.40446199974765351*tmp_29;
      real_t tmp_503 = tmp_19*(tmp_31 + tmp_502);
      real_t tmp_504 = 0.40446199974765351*tmp_35 + 0.40446199974765351*tmp_36;
      real_t tmp_505 = tmp_19*(tmp_38 + tmp_504);
      real_t tmp_506 = tmp_1*(tmp_26*tmp_503 + tmp_33*tmp_505 + tmp_501*tmp_7 - 1.0/4.0) + tmp_2*(tmp_40*tmp_501 + tmp_41*tmp_503 + tmp_42*tmp_505 - 1.0/4.0) + tmp_5*(tmp_43*tmp_501 + tmp_44*tmp_503 + tmp_45*tmp_505 - 1.0/4.0);
      real_t tmp_507 = tmp_500 + tmp_76;
      real_t tmp_508 = tmp_507*tmp_72;
      real_t tmp_509 = tmp_507*tmp_73;
      real_t tmp_510 = tmp_502 + tmp_80;
      real_t tmp_511 = tmp_510*tmp_69;
      real_t tmp_512 = tmp_510*tmp_70;
      real_t tmp_513 = tmp_507*tmp_74;
      real_t tmp_514 = tmp_510*tmp_71;
      real_t tmp_515 = tmp_504 + tmp_86;
      real_t tmp_516 = tmp_515*tmp_66;
      real_t tmp_517 = tmp_515*tmp_67;
      real_t tmp_518 = tmp_515*tmp_68;
      real_t tmp_519 = -tmp_508 - tmp_509 - tmp_511 - tmp_512 - tmp_513 - tmp_514 - tmp_516 - tmp_517 - tmp_518 + 1;
      real_t tmp_520 = tmp_100*tmp_506;
      real_t tmp_521 = 0.042507265838595799*tmp_102;
      real_t tmp_522 = 0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_523 = tmp_19*(tmp_24 + tmp_522);
      real_t tmp_524 = 0.1711304259088916*tmp_28 + 0.041227165399737475*tmp_29;
      real_t tmp_525 = tmp_19*(tmp_31 + tmp_524);
      real_t tmp_526 = 0.1711304259088916*tmp_35 + 0.041227165399737475*tmp_36;
      real_t tmp_527 = tmp_19*(tmp_38 + tmp_526);
      real_t tmp_528 = tmp_1*(tmp_26*tmp_525 + tmp_33*tmp_527 + tmp_523*tmp_7 - 1.0/4.0) + tmp_2*(tmp_40*tmp_523 + tmp_41*tmp_525 + tmp_42*tmp_527 - 1.0/4.0) + tmp_5*(tmp_43*tmp_523 + tmp_44*tmp_525 + tmp_45*tmp_527 - 1.0/4.0);
      real_t tmp_529 = tmp_522 + tmp_76;
      real_t tmp_530 = tmp_529*tmp_72;
      real_t tmp_531 = tmp_529*tmp_73;
      real_t tmp_532 = tmp_524 + tmp_80;
      real_t tmp_533 = tmp_532*tmp_69;
      real_t tmp_534 = tmp_532*tmp_70;
      real_t tmp_535 = tmp_529*tmp_74;
      real_t tmp_536 = tmp_532*tmp_71;
      real_t tmp_537 = tmp_526 + tmp_86;
      real_t tmp_538 = tmp_537*tmp_66;
      real_t tmp_539 = tmp_537*tmp_67;
      real_t tmp_540 = tmp_537*tmp_68;
      real_t tmp_541 = -tmp_530 - tmp_531 - tmp_533 - tmp_534 - tmp_535 - tmp_536 - tmp_538 - tmp_539 - tmp_540 + 1;
      real_t tmp_542 = tmp_100*tmp_528;
      real_t tmp_543 = 0.019202922745021479*tmp_102;
      real_t tmp_544 = tmp_84 + tmp_85 + tmp_90;
      real_t tmp_545 = 0.5*p_affine_13_0*tmp_68 + 0.5*p_affine_13_1*tmp_71 + 0.5*p_affine_13_2*tmp_74;
      real_t tmp_546 = tmp_117 + tmp_118 + tmp_122;
      real_t tmp_547 = tmp_139 + tmp_140 + tmp_144;
      real_t tmp_548 = tmp_161 + tmp_162 + tmp_166;
      real_t tmp_549 = tmp_183 + tmp_184 + tmp_188;
      real_t tmp_550 = tmp_205 + tmp_206 + tmp_210;
      real_t tmp_551 = tmp_227 + tmp_228 + tmp_232;
      real_t tmp_552 = tmp_249 + tmp_250 + tmp_254;
      real_t tmp_553 = tmp_271 + tmp_272 + tmp_276;
      real_t tmp_554 = tmp_293 + tmp_294 + tmp_298;
      real_t tmp_555 = tmp_315 + tmp_316 + tmp_320;
      real_t tmp_556 = tmp_337 + tmp_338 + tmp_342;
      real_t tmp_557 = tmp_359 + tmp_360 + tmp_364;
      real_t tmp_558 = tmp_381 + tmp_382 + tmp_386;
      real_t tmp_559 = tmp_403 + tmp_404 + tmp_408;
      real_t tmp_560 = tmp_425 + tmp_426 + tmp_430;
      real_t tmp_561 = tmp_447 + tmp_448 + tmp_452;
      real_t tmp_562 = tmp_469 + tmp_470 + tmp_474;
      real_t tmp_563 = tmp_491 + tmp_492 + tmp_496;
      real_t tmp_564 = tmp_513 + tmp_514 + tmp_518;
      real_t tmp_565 = tmp_535 + tmp_536 + tmp_540;
      real_t tmp_566 = tmp_79 + tmp_83 + tmp_89;
      real_t tmp_567 = 0.5*p_affine_13_0*tmp_67 + 0.5*p_affine_13_1*tmp_70 + 0.5*p_affine_13_2*tmp_73;
      real_t tmp_568 = tmp_113 + tmp_116 + tmp_121;
      real_t tmp_569 = tmp_135 + tmp_138 + tmp_143;
      real_t tmp_570 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_571 = tmp_179 + tmp_182 + tmp_187;
      real_t tmp_572 = tmp_201 + tmp_204 + tmp_209;
      real_t tmp_573 = tmp_223 + tmp_226 + tmp_231;
      real_t tmp_574 = tmp_245 + tmp_248 + tmp_253;
      real_t tmp_575 = tmp_267 + tmp_270 + tmp_275;
      real_t tmp_576 = tmp_289 + tmp_292 + tmp_297;
      real_t tmp_577 = tmp_311 + tmp_314 + tmp_319;
      real_t tmp_578 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_579 = tmp_355 + tmp_358 + tmp_363;
      real_t tmp_580 = tmp_377 + tmp_380 + tmp_385;
      real_t tmp_581 = tmp_399 + tmp_402 + tmp_407;
      real_t tmp_582 = tmp_421 + tmp_424 + tmp_429;
      real_t tmp_583 = tmp_443 + tmp_446 + tmp_451;
      real_t tmp_584 = tmp_465 + tmp_468 + tmp_473;
      real_t tmp_585 = tmp_487 + tmp_490 + tmp_495;
      real_t tmp_586 = tmp_509 + tmp_512 + tmp_517;
      real_t tmp_587 = tmp_531 + tmp_534 + tmp_539;
      real_t tmp_588 = tmp_78 + tmp_82 + tmp_88;
      real_t tmp_589 = 0.5*p_affine_13_0*tmp_66 + 0.5*p_affine_13_1*tmp_69 + 0.5*p_affine_13_2*tmp_72;
      real_t tmp_590 = tmp_112 + tmp_115 + tmp_120;
      real_t tmp_591 = tmp_134 + tmp_137 + tmp_142;
      real_t tmp_592 = tmp_156 + tmp_159 + tmp_164;
      real_t tmp_593 = tmp_178 + tmp_181 + tmp_186;
      real_t tmp_594 = tmp_200 + tmp_203 + tmp_208;
      real_t tmp_595 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_596 = tmp_244 + tmp_247 + tmp_252;
      real_t tmp_597 = tmp_266 + tmp_269 + tmp_274;
      real_t tmp_598 = tmp_288 + tmp_291 + tmp_296;
      real_t tmp_599 = tmp_310 + tmp_313 + tmp_318;
      real_t tmp_600 = tmp_332 + tmp_335 + tmp_340;
      real_t tmp_601 = tmp_354 + tmp_357 + tmp_362;
      real_t tmp_602 = tmp_376 + tmp_379 + tmp_384;
      real_t tmp_603 = tmp_398 + tmp_401 + tmp_406;
      real_t tmp_604 = tmp_420 + tmp_423 + tmp_428;
      real_t tmp_605 = tmp_442 + tmp_445 + tmp_450;
      real_t tmp_606 = tmp_464 + tmp_467 + tmp_472;
      real_t tmp_607 = tmp_486 + tmp_489 + tmp_494;
      real_t tmp_608 = tmp_508 + tmp_511 + tmp_516;
      real_t tmp_609 = tmp_530 + tmp_533 + tmp_538;
      real_t a_0_0 = tmp_103*(-tmp_101*tmp_91 - tmp_46*tmp_75 + tmp_91*tmp_95) + tmp_125*(-tmp_110*tmp_75 - tmp_123*tmp_124 + tmp_123*tmp_95) + tmp_147*(-tmp_132*tmp_75 - tmp_145*tmp_146 + tmp_145*tmp_95) + tmp_169*(-tmp_154*tmp_75 - tmp_167*tmp_168 + tmp_167*tmp_95) + tmp_191*(-tmp_176*tmp_75 - tmp_189*tmp_190 + tmp_189*tmp_95) + tmp_213*(-tmp_198*tmp_75 - tmp_211*tmp_212 + tmp_211*tmp_95) + tmp_235*(-tmp_220*tmp_75 - tmp_233*tmp_234 + tmp_233*tmp_95) + tmp_257*(-tmp_242*tmp_75 - tmp_255*tmp_256 + tmp_255*tmp_95) + tmp_279*(-tmp_264*tmp_75 - tmp_277*tmp_278 + tmp_277*tmp_95) + tmp_301*(-tmp_286*tmp_75 - tmp_299*tmp_300 + tmp_299*tmp_95) + tmp_323*(-tmp_308*tmp_75 - tmp_321*tmp_322 + tmp_321*tmp_95) + tmp_345*(-tmp_330*tmp_75 - tmp_343*tmp_344 + tmp_343*tmp_95) + tmp_367*(-tmp_352*tmp_75 - tmp_365*tmp_366 + tmp_365*tmp_95) + tmp_389*(-tmp_374*tmp_75 - tmp_387*tmp_388 + tmp_387*tmp_95) + tmp_411*(-tmp_396*tmp_75 - tmp_409*tmp_410 + tmp_409*tmp_95) + tmp_433*(-tmp_418*tmp_75 - tmp_431*tmp_432 + tmp_431*tmp_95) + tmp_455*(-tmp_440*tmp_75 - tmp_453*tmp_454 + tmp_453*tmp_95) + tmp_477*(-tmp_462*tmp_75 - tmp_475*tmp_476 + tmp_475*tmp_95) + tmp_499*(-tmp_484*tmp_75 - tmp_497*tmp_498 + tmp_497*tmp_95) + tmp_521*(-tmp_506*tmp_75 - tmp_519*tmp_520 + tmp_519*tmp_95) + tmp_543*(-tmp_528*tmp_75 - tmp_541*tmp_542 + tmp_541*tmp_95);
      real_t a_0_1 = tmp_103*(-tmp_101*tmp_544 - tmp_46*tmp_545 + tmp_544*tmp_95) + tmp_125*(-tmp_110*tmp_545 - tmp_124*tmp_546 + tmp_546*tmp_95) + tmp_147*(-tmp_132*tmp_545 - tmp_146*tmp_547 + tmp_547*tmp_95) + tmp_169*(-tmp_154*tmp_545 - tmp_168*tmp_548 + tmp_548*tmp_95) + tmp_191*(-tmp_176*tmp_545 - tmp_190*tmp_549 + tmp_549*tmp_95) + tmp_213*(-tmp_198*tmp_545 - tmp_212*tmp_550 + tmp_550*tmp_95) + tmp_235*(-tmp_220*tmp_545 - tmp_234*tmp_551 + tmp_551*tmp_95) + tmp_257*(-tmp_242*tmp_545 - tmp_256*tmp_552 + tmp_552*tmp_95) + tmp_279*(-tmp_264*tmp_545 - tmp_278*tmp_553 + tmp_553*tmp_95) + tmp_301*(-tmp_286*tmp_545 - tmp_300*tmp_554 + tmp_554*tmp_95) + tmp_323*(-tmp_308*tmp_545 - tmp_322*tmp_555 + tmp_555*tmp_95) + tmp_345*(-tmp_330*tmp_545 - tmp_344*tmp_556 + tmp_556*tmp_95) + tmp_367*(-tmp_352*tmp_545 - tmp_366*tmp_557 + tmp_557*tmp_95) + tmp_389*(-tmp_374*tmp_545 - tmp_388*tmp_558 + tmp_558*tmp_95) + tmp_411*(-tmp_396*tmp_545 - tmp_410*tmp_559 + tmp_559*tmp_95) + tmp_433*(-tmp_418*tmp_545 - tmp_432*tmp_560 + tmp_560*tmp_95) + tmp_455*(-tmp_440*tmp_545 - tmp_454*tmp_561 + tmp_561*tmp_95) + tmp_477*(-tmp_462*tmp_545 - tmp_476*tmp_562 + tmp_562*tmp_95) + tmp_499*(-tmp_484*tmp_545 - tmp_498*tmp_563 + tmp_563*tmp_95) + tmp_521*(-tmp_506*tmp_545 - tmp_520*tmp_564 + tmp_564*tmp_95) + tmp_543*(-tmp_528*tmp_545 - tmp_542*tmp_565 + tmp_565*tmp_95);
      real_t a_0_2 = tmp_103*(-tmp_101*tmp_566 - tmp_46*tmp_567 + tmp_566*tmp_95) + tmp_125*(-tmp_110*tmp_567 - tmp_124*tmp_568 + tmp_568*tmp_95) + tmp_147*(-tmp_132*tmp_567 - tmp_146*tmp_569 + tmp_569*tmp_95) + tmp_169*(-tmp_154*tmp_567 - tmp_168*tmp_570 + tmp_570*tmp_95) + tmp_191*(-tmp_176*tmp_567 - tmp_190*tmp_571 + tmp_571*tmp_95) + tmp_213*(-tmp_198*tmp_567 - tmp_212*tmp_572 + tmp_572*tmp_95) + tmp_235*(-tmp_220*tmp_567 - tmp_234*tmp_573 + tmp_573*tmp_95) + tmp_257*(-tmp_242*tmp_567 - tmp_256*tmp_574 + tmp_574*tmp_95) + tmp_279*(-tmp_264*tmp_567 - tmp_278*tmp_575 + tmp_575*tmp_95) + tmp_301*(-tmp_286*tmp_567 - tmp_300*tmp_576 + tmp_576*tmp_95) + tmp_323*(-tmp_308*tmp_567 - tmp_322*tmp_577 + tmp_577*tmp_95) + tmp_345*(-tmp_330*tmp_567 - tmp_344*tmp_578 + tmp_578*tmp_95) + tmp_367*(-tmp_352*tmp_567 - tmp_366*tmp_579 + tmp_579*tmp_95) + tmp_389*(-tmp_374*tmp_567 - tmp_388*tmp_580 + tmp_580*tmp_95) + tmp_411*(-tmp_396*tmp_567 - tmp_410*tmp_581 + tmp_581*tmp_95) + tmp_433*(-tmp_418*tmp_567 - tmp_432*tmp_582 + tmp_582*tmp_95) + tmp_455*(-tmp_440*tmp_567 - tmp_454*tmp_583 + tmp_583*tmp_95) + tmp_477*(-tmp_462*tmp_567 - tmp_476*tmp_584 + tmp_584*tmp_95) + tmp_499*(-tmp_484*tmp_567 - tmp_498*tmp_585 + tmp_585*tmp_95) + tmp_521*(-tmp_506*tmp_567 - tmp_520*tmp_586 + tmp_586*tmp_95) + tmp_543*(-tmp_528*tmp_567 - tmp_542*tmp_587 + tmp_587*tmp_95);
      real_t a_0_3 = tmp_103*(-tmp_101*tmp_588 - tmp_46*tmp_589 + tmp_588*tmp_95) + tmp_125*(-tmp_110*tmp_589 - tmp_124*tmp_590 + tmp_590*tmp_95) + tmp_147*(-tmp_132*tmp_589 - tmp_146*tmp_591 + tmp_591*tmp_95) + tmp_169*(-tmp_154*tmp_589 - tmp_168*tmp_592 + tmp_592*tmp_95) + tmp_191*(-tmp_176*tmp_589 - tmp_190*tmp_593 + tmp_593*tmp_95) + tmp_213*(-tmp_198*tmp_589 - tmp_212*tmp_594 + tmp_594*tmp_95) + tmp_235*(-tmp_220*tmp_589 - tmp_234*tmp_595 + tmp_595*tmp_95) + tmp_257*(-tmp_242*tmp_589 - tmp_256*tmp_596 + tmp_596*tmp_95) + tmp_279*(-tmp_264*tmp_589 - tmp_278*tmp_597 + tmp_597*tmp_95) + tmp_301*(-tmp_286*tmp_589 - tmp_300*tmp_598 + tmp_598*tmp_95) + tmp_323*(-tmp_308*tmp_589 - tmp_322*tmp_599 + tmp_599*tmp_95) + tmp_345*(-tmp_330*tmp_589 - tmp_344*tmp_600 + tmp_600*tmp_95) + tmp_367*(-tmp_352*tmp_589 - tmp_366*tmp_601 + tmp_601*tmp_95) + tmp_389*(-tmp_374*tmp_589 - tmp_388*tmp_602 + tmp_602*tmp_95) + tmp_411*(-tmp_396*tmp_589 - tmp_410*tmp_603 + tmp_603*tmp_95) + tmp_433*(-tmp_418*tmp_589 - tmp_432*tmp_604 + tmp_604*tmp_95) + tmp_455*(-tmp_440*tmp_589 - tmp_454*tmp_605 + tmp_605*tmp_95) + tmp_477*(-tmp_462*tmp_589 - tmp_476*tmp_606 + tmp_606*tmp_95) + tmp_499*(-tmp_484*tmp_589 - tmp_498*tmp_607 + tmp_607*tmp_95) + tmp_521*(-tmp_506*tmp_589 - tmp_520*tmp_608 + tmp_608*tmp_95) + tmp_543*(-tmp_528*tmp_589 - tmp_542*tmp_609 + tmp_609*tmp_95);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
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
      real_t tmp_22 = tmp_18*(tmp_12*tmp_6 - tmp_3*tmp_9);
      real_t tmp_23 = tmp_18*(tmp_10*tmp_9 - tmp_15*tmp_6);
      real_t tmp_24 = tmp_18*(-tmp_10*tmp_12 + tmp_15*tmp_3);
      real_t tmp_25 = tmp_18*(-tmp_1*tmp_12 + tmp_5*tmp_9);
      real_t tmp_26 = tmp_18*(tmp_1*tmp_15 - tmp_13*tmp_9);
      real_t tmp_27 = tmp_18*(tmp_12*tmp_13 - tmp_15*tmp_5);
      real_t tmp_28 = -p_affine_8_0;
      real_t tmp_29 = p_affine_10_0 + tmp_28;
      real_t tmp_30 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_31 = -p_affine_8_1;
      real_t tmp_32 = p_affine_10_1 + tmp_31;
      real_t tmp_33 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_34 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_35 = -p_affine_8_2;
      real_t tmp_36 = p_affine_10_2 + tmp_35;
      real_t tmp_37 = 1.0*std::pow((std::abs(tmp_29*tmp_30 - tmp_32*tmp_33)*std::abs(tmp_29*tmp_30 - tmp_32*tmp_33)) + (std::abs(tmp_29*tmp_34 - tmp_33*tmp_36)*std::abs(tmp_29*tmp_34 - tmp_33*tmp_36)) + (std::abs(tmp_30*tmp_36 - tmp_32*tmp_34)*std::abs(tmp_30*tmp_36 - tmp_32*tmp_34)), 1.0/2.0);
      real_t tmp_38 = tmp_37*(p_affine_13_0*(-tmp_19 - tmp_20 - tmp_21) + p_affine_13_1*(-tmp_22 - tmp_23 - tmp_24) + p_affine_13_2*(-tmp_25 - tmp_26 - tmp_27));
      real_t tmp_39 = p_affine_9_2 + tmp_35;
      real_t tmp_40 = p_affine_8_2 + tmp_2;
      real_t tmp_41 = 0.93718850182767688*tmp_36 + 0.031405749086161582*tmp_39 + tmp_40;
      real_t tmp_42 = p_affine_9_1 + tmp_31;
      real_t tmp_43 = p_affine_8_1 + tmp_0;
      real_t tmp_44 = 0.93718850182767688*tmp_32 + 0.031405749086161582*tmp_42 + tmp_43;
      real_t tmp_45 = p_affine_9_0 + tmp_28;
      real_t tmp_46 = p_affine_8_0 + tmp_8;
      real_t tmp_47 = 0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_45 + tmp_46;
      real_t tmp_48 = 0.0068572537431980923*tmp_12*(tmp_20*tmp_47 + tmp_23*tmp_44 + tmp_26*tmp_41 - 1.0/4.0) + 0.0068572537431980923*tmp_15*(tmp_19*tmp_47 + tmp_22*tmp_44 + tmp_25*tmp_41 - 1.0/4.0) + 0.0068572537431980923*tmp_9*(tmp_21*tmp_47 + tmp_24*tmp_44 + tmp_27*tmp_41 - 1.0/4.0);
      real_t tmp_49 = 0.60796128279561268*tmp_36 + 0.19601935860219369*tmp_39 + tmp_40;
      real_t tmp_50 = 0.60796128279561268*tmp_32 + 0.19601935860219369*tmp_42 + tmp_43;
      real_t tmp_51 = 0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_45 + tmp_46;
      real_t tmp_52 = 0.037198804536718075*tmp_12*(tmp_20*tmp_51 + tmp_23*tmp_50 + tmp_26*tmp_49 - 1.0/4.0) + 0.037198804536718075*tmp_15*(tmp_19*tmp_51 + tmp_22*tmp_50 + tmp_25*tmp_49 - 1.0/4.0) + 0.037198804536718075*tmp_9*(tmp_21*tmp_51 + tmp_24*tmp_50 + tmp_27*tmp_49 - 1.0/4.0);
      real_t tmp_53 = 0.039308471900058539*tmp_36 + 0.37605877282253791*tmp_39 + tmp_40;
      real_t tmp_54 = 0.039308471900058539*tmp_32 + 0.37605877282253791*tmp_42 + tmp_43;
      real_t tmp_55 = 0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_45 + tmp_46;
      real_t tmp_56 = 0.020848748529055869*tmp_12*(tmp_20*tmp_55 + tmp_23*tmp_54 + tmp_26*tmp_53 - 1.0/4.0) + 0.020848748529055869*tmp_15*(tmp_19*tmp_55 + tmp_22*tmp_54 + tmp_25*tmp_53 - 1.0/4.0) + 0.020848748529055869*tmp_9*(tmp_21*tmp_55 + tmp_24*tmp_54 + tmp_27*tmp_53 - 1.0/4.0);
      real_t tmp_57 = 0.1711304259088916*tmp_36 + 0.78764240869137092*tmp_39 + tmp_40;
      real_t tmp_58 = 0.1711304259088916*tmp_32 + 0.78764240869137092*tmp_42 + tmp_43;
      real_t tmp_59 = 0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_45 + tmp_46;
      real_t tmp_60 = 0.019202922745021479*tmp_12*(tmp_20*tmp_59 + tmp_23*tmp_58 + tmp_26*tmp_57 - 1.0/4.0) + 0.019202922745021479*tmp_15*(tmp_19*tmp_59 + tmp_22*tmp_58 + tmp_25*tmp_57 - 1.0/4.0) + 0.019202922745021479*tmp_9*(tmp_21*tmp_59 + tmp_24*tmp_58 + tmp_27*tmp_57 - 1.0/4.0);
      real_t tmp_61 = 0.37605877282253791*tmp_36 + 0.58463275527740355*tmp_39 + tmp_40;
      real_t tmp_62 = 0.37605877282253791*tmp_32 + 0.58463275527740355*tmp_42 + tmp_43;
      real_t tmp_63 = 0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_45 + tmp_46;
      real_t tmp_64 = 0.020848748529055869*tmp_12*(tmp_20*tmp_63 + tmp_23*tmp_62 + tmp_26*tmp_61 - 1.0/4.0) + 0.020848748529055869*tmp_15*(tmp_19*tmp_63 + tmp_22*tmp_62 + tmp_25*tmp_61 - 1.0/4.0) + 0.020848748529055869*tmp_9*(tmp_21*tmp_63 + tmp_24*tmp_62 + tmp_27*tmp_61 - 1.0/4.0);
      real_t tmp_65 = 0.78764240869137092*tmp_36 + 0.041227165399737475*tmp_39 + tmp_40;
      real_t tmp_66 = 0.78764240869137092*tmp_32 + 0.041227165399737475*tmp_42 + tmp_43;
      real_t tmp_67 = 0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_45 + tmp_46;
      real_t tmp_68 = 0.019202922745021479*tmp_12*(tmp_20*tmp_67 + tmp_23*tmp_66 + tmp_26*tmp_65 - 1.0/4.0) + 0.019202922745021479*tmp_15*(tmp_19*tmp_67 + tmp_22*tmp_66 + tmp_25*tmp_65 - 1.0/4.0) + 0.019202922745021479*tmp_9*(tmp_21*tmp_67 + tmp_24*tmp_66 + tmp_27*tmp_65 - 1.0/4.0);
      real_t tmp_69 = 0.58463275527740355*tmp_36 + 0.039308471900058539*tmp_39 + tmp_40;
      real_t tmp_70 = 0.58463275527740355*tmp_32 + 0.039308471900058539*tmp_42 + tmp_43;
      real_t tmp_71 = 0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_45 + tmp_46;
      real_t tmp_72 = 0.020848748529055869*tmp_12*(tmp_20*tmp_71 + tmp_23*tmp_70 + tmp_26*tmp_69 - 1.0/4.0) + 0.020848748529055869*tmp_15*(tmp_19*tmp_71 + tmp_22*tmp_70 + tmp_25*tmp_69 - 1.0/4.0) + 0.020848748529055869*tmp_9*(tmp_21*tmp_71 + tmp_24*tmp_70 + tmp_27*tmp_69 - 1.0/4.0);
      real_t tmp_73 = 0.041227165399737475*tmp_36 + 0.78764240869137092*tmp_39 + tmp_40;
      real_t tmp_74 = 0.041227165399737475*tmp_32 + 0.78764240869137092*tmp_42 + tmp_43;
      real_t tmp_75 = 0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_45 + tmp_46;
      real_t tmp_76 = 0.019202922745021479*tmp_12*(tmp_20*tmp_75 + tmp_23*tmp_74 + tmp_26*tmp_73 - 1.0/4.0) + 0.019202922745021479*tmp_15*(tmp_19*tmp_75 + tmp_22*tmp_74 + tmp_25*tmp_73 - 1.0/4.0) + 0.019202922745021479*tmp_9*(tmp_21*tmp_75 + tmp_24*tmp_74 + tmp_27*tmp_73 - 1.0/4.0);
      real_t tmp_77 = 0.039308471900058539*tmp_36 + 0.58463275527740355*tmp_39 + tmp_40;
      real_t tmp_78 = 0.039308471900058539*tmp_32 + 0.58463275527740355*tmp_42 + tmp_43;
      real_t tmp_79 = 0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_45 + tmp_46;
      real_t tmp_80 = 0.020848748529055869*tmp_12*(tmp_20*tmp_79 + tmp_23*tmp_78 + tmp_26*tmp_77 - 1.0/4.0) + 0.020848748529055869*tmp_15*(tmp_19*tmp_79 + tmp_22*tmp_78 + tmp_25*tmp_77 - 1.0/4.0) + 0.020848748529055869*tmp_9*(tmp_21*tmp_79 + tmp_24*tmp_78 + tmp_27*tmp_77 - 1.0/4.0);
      real_t tmp_81 = 0.78764240869137092*tmp_36 + 0.1711304259088916*tmp_39 + tmp_40;
      real_t tmp_82 = 0.78764240869137092*tmp_32 + 0.1711304259088916*tmp_42 + tmp_43;
      real_t tmp_83 = 0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_45 + tmp_46;
      real_t tmp_84 = 0.019202922745021479*tmp_12*(tmp_20*tmp_83 + tmp_23*tmp_82 + tmp_26*tmp_81 - 1.0/4.0) + 0.019202922745021479*tmp_15*(tmp_19*tmp_83 + tmp_22*tmp_82 + tmp_25*tmp_81 - 1.0/4.0) + 0.019202922745021479*tmp_9*(tmp_21*tmp_83 + tmp_24*tmp_82 + tmp_27*tmp_81 - 1.0/4.0);
      real_t tmp_85 = 0.58463275527740355*tmp_36 + 0.37605877282253791*tmp_39 + tmp_40;
      real_t tmp_86 = 0.58463275527740355*tmp_32 + 0.37605877282253791*tmp_42 + tmp_43;
      real_t tmp_87 = 0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_45 + tmp_46;
      real_t tmp_88 = 0.020848748529055869*tmp_12*(tmp_20*tmp_87 + tmp_23*tmp_86 + tmp_26*tmp_85 - 1.0/4.0) + 0.020848748529055869*tmp_15*(tmp_19*tmp_87 + tmp_22*tmp_86 + tmp_25*tmp_85 - 1.0/4.0) + 0.020848748529055869*tmp_9*(tmp_21*tmp_87 + tmp_24*tmp_86 + tmp_27*tmp_85 - 1.0/4.0);
      real_t tmp_89 = 0.1711304259088916*tmp_36 + 0.041227165399737475*tmp_39 + tmp_40;
      real_t tmp_90 = 0.1711304259088916*tmp_32 + 0.041227165399737475*tmp_42 + tmp_43;
      real_t tmp_91 = 0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_45 + tmp_46;
      real_t tmp_92 = 0.019202922745021479*tmp_12*(tmp_20*tmp_91 + tmp_23*tmp_90 + tmp_26*tmp_89 - 1.0/4.0) + 0.019202922745021479*tmp_15*(tmp_19*tmp_91 + tmp_22*tmp_90 + tmp_25*tmp_89 - 1.0/4.0) + 0.019202922745021479*tmp_9*(tmp_21*tmp_91 + tmp_24*tmp_90 + tmp_27*tmp_89 - 1.0/4.0);
      real_t tmp_93 = 0.19107600050469298*tmp_36 + 0.40446199974765351*tmp_39 + tmp_40;
      real_t tmp_94 = 0.19107600050469298*tmp_32 + 0.40446199974765351*tmp_42 + tmp_43;
      real_t tmp_95 = 0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_45 + tmp_46;
      real_t tmp_96 = 0.042507265838595799*tmp_12*(tmp_20*tmp_95 + tmp_23*tmp_94 + tmp_26*tmp_93 - 1.0/4.0) + 0.042507265838595799*tmp_15*(tmp_19*tmp_95 + tmp_22*tmp_94 + tmp_25*tmp_93 - 1.0/4.0) + 0.042507265838595799*tmp_9*(tmp_21*tmp_95 + tmp_24*tmp_94 + tmp_27*tmp_93 - 1.0/4.0);
      real_t tmp_97 = 0.37605877282253791*tmp_36 + 0.039308471900058539*tmp_39 + tmp_40;
      real_t tmp_98 = 0.37605877282253791*tmp_32 + 0.039308471900058539*tmp_42 + tmp_43;
      real_t tmp_99 = 0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_45 + tmp_46;
      real_t tmp_100 = 0.020848748529055869*tmp_12*(tmp_20*tmp_99 + tmp_23*tmp_98 + tmp_26*tmp_97 - 1.0/4.0) + 0.020848748529055869*tmp_15*(tmp_19*tmp_99 + tmp_22*tmp_98 + tmp_25*tmp_97 - 1.0/4.0) + 0.020848748529055869*tmp_9*(tmp_21*tmp_99 + tmp_24*tmp_98 + tmp_27*tmp_97 - 1.0/4.0);
      real_t tmp_101 = 0.031405749086161582*tmp_36 + 0.93718850182767688*tmp_39 + tmp_40;
      real_t tmp_102 = 0.031405749086161582*tmp_32 + 0.93718850182767688*tmp_42 + tmp_43;
      real_t tmp_103 = 0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_45 + tmp_46;
      real_t tmp_104 = 0.0068572537431980923*tmp_12*(tmp_101*tmp_26 + tmp_102*tmp_23 + tmp_103*tmp_20 - 1.0/4.0) + 0.0068572537431980923*tmp_15*(tmp_101*tmp_25 + tmp_102*tmp_22 + tmp_103*tmp_19 - 1.0/4.0) + 0.0068572537431980923*tmp_9*(tmp_101*tmp_27 + tmp_102*tmp_24 + tmp_103*tmp_21 - 1.0/4.0);
      real_t tmp_105 = 0.19601935860219369*tmp_36 + 0.60796128279561268*tmp_39 + tmp_40;
      real_t tmp_106 = 0.19601935860219369*tmp_32 + 0.60796128279561268*tmp_42 + tmp_43;
      real_t tmp_107 = 0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_45 + tmp_46;
      real_t tmp_108 = 0.037198804536718075*tmp_12*(tmp_105*tmp_26 + tmp_106*tmp_23 + tmp_107*tmp_20 - 1.0/4.0) + 0.037198804536718075*tmp_15*(tmp_105*tmp_25 + tmp_106*tmp_22 + tmp_107*tmp_19 - 1.0/4.0) + 0.037198804536718075*tmp_9*(tmp_105*tmp_27 + tmp_106*tmp_24 + tmp_107*tmp_21 - 1.0/4.0);
      real_t tmp_109 = 0.40446199974765351*tmp_36 + 0.19107600050469298*tmp_39 + tmp_40;
      real_t tmp_110 = 0.40446199974765351*tmp_32 + 0.19107600050469298*tmp_42 + tmp_43;
      real_t tmp_111 = 0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_45 + tmp_46;
      real_t tmp_112 = 0.042507265838595799*tmp_12*(tmp_109*tmp_26 + tmp_110*tmp_23 + tmp_111*tmp_20 - 1.0/4.0) + 0.042507265838595799*tmp_15*(tmp_109*tmp_25 + tmp_110*tmp_22 + tmp_111*tmp_19 - 1.0/4.0) + 0.042507265838595799*tmp_9*(tmp_109*tmp_27 + tmp_110*tmp_24 + tmp_111*tmp_21 - 1.0/4.0);
      real_t tmp_113 = 0.031405749086161582*tmp_36 + 0.031405749086161582*tmp_39 + tmp_40;
      real_t tmp_114 = 0.031405749086161582*tmp_32 + 0.031405749086161582*tmp_42 + tmp_43;
      real_t tmp_115 = 0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_45 + tmp_46;
      real_t tmp_116 = 0.0068572537431980923*tmp_12*(tmp_113*tmp_26 + tmp_114*tmp_23 + tmp_115*tmp_20 - 1.0/4.0) + 0.0068572537431980923*tmp_15*(tmp_113*tmp_25 + tmp_114*tmp_22 + tmp_115*tmp_19 - 1.0/4.0) + 0.0068572537431980923*tmp_9*(tmp_113*tmp_27 + tmp_114*tmp_24 + tmp_115*tmp_21 - 1.0/4.0);
      real_t tmp_117 = 0.19601935860219369*tmp_36 + 0.19601935860219369*tmp_39 + tmp_40;
      real_t tmp_118 = 0.19601935860219369*tmp_32 + 0.19601935860219369*tmp_42 + tmp_43;
      real_t tmp_119 = 0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_45 + tmp_46;
      real_t tmp_120 = 0.037198804536718075*tmp_12*(tmp_117*tmp_26 + tmp_118*tmp_23 + tmp_119*tmp_20 - 1.0/4.0) + 0.037198804536718075*tmp_15*(tmp_117*tmp_25 + tmp_118*tmp_22 + tmp_119*tmp_19 - 1.0/4.0) + 0.037198804536718075*tmp_9*(tmp_117*tmp_27 + tmp_118*tmp_24 + tmp_119*tmp_21 - 1.0/4.0);
      real_t tmp_121 = 0.40446199974765351*tmp_36 + 0.40446199974765351*tmp_39 + tmp_40;
      real_t tmp_122 = 0.40446199974765351*tmp_32 + 0.40446199974765351*tmp_42 + tmp_43;
      real_t tmp_123 = 0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_45 + tmp_46;
      real_t tmp_124 = 0.042507265838595799*tmp_12*(tmp_121*tmp_26 + tmp_122*tmp_23 + tmp_123*tmp_20 - 1.0/4.0) + 0.042507265838595799*tmp_15*(tmp_121*tmp_25 + tmp_122*tmp_22 + tmp_123*tmp_19 - 1.0/4.0) + 0.042507265838595799*tmp_9*(tmp_121*tmp_27 + tmp_122*tmp_24 + tmp_123*tmp_21 - 1.0/4.0);
      real_t tmp_125 = 0.041227165399737475*tmp_36 + 0.1711304259088916*tmp_39 + tmp_40;
      real_t tmp_126 = 0.041227165399737475*tmp_32 + 0.1711304259088916*tmp_42 + tmp_43;
      real_t tmp_127 = 0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_45 + tmp_46;
      real_t tmp_128 = 0.019202922745021479*tmp_12*(tmp_125*tmp_26 + tmp_126*tmp_23 + tmp_127*tmp_20 - 1.0/4.0) + 0.019202922745021479*tmp_15*(tmp_125*tmp_25 + tmp_126*tmp_22 + tmp_127*tmp_19 - 1.0/4.0) + 0.019202922745021479*tmp_9*(tmp_125*tmp_27 + tmp_126*tmp_24 + tmp_127*tmp_21 - 1.0/4.0);
      real_t tmp_129 = tmp_37*(p_affine_13_0*tmp_21 + p_affine_13_1*tmp_24 + p_affine_13_2*tmp_27);
      real_t tmp_130 = tmp_37*(p_affine_13_0*tmp_20 + p_affine_13_1*tmp_23 + p_affine_13_2*tmp_26);
      real_t tmp_131 = tmp_37*(p_affine_13_0*tmp_19 + p_affine_13_1*tmp_22 + p_affine_13_2*tmp_25);
      real_t a_0_0 = -tmp_100*tmp_38 - tmp_104*tmp_38 - tmp_108*tmp_38 - tmp_112*tmp_38 - tmp_116*tmp_38 - tmp_120*tmp_38 - tmp_124*tmp_38 - tmp_128*tmp_38 - tmp_38*tmp_48 - tmp_38*tmp_52 - tmp_38*tmp_56 - tmp_38*tmp_60 - tmp_38*tmp_64 - tmp_38*tmp_68 - tmp_38*tmp_72 - tmp_38*tmp_76 - tmp_38*tmp_80 - tmp_38*tmp_84 - tmp_38*tmp_88 - tmp_38*tmp_92 - tmp_38*tmp_96;
      real_t a_0_1 = -tmp_100*tmp_129 - tmp_104*tmp_129 - tmp_108*tmp_129 - tmp_112*tmp_129 - tmp_116*tmp_129 - tmp_120*tmp_129 - tmp_124*tmp_129 - tmp_128*tmp_129 - tmp_129*tmp_48 - tmp_129*tmp_52 - tmp_129*tmp_56 - tmp_129*tmp_60 - tmp_129*tmp_64 - tmp_129*tmp_68 - tmp_129*tmp_72 - tmp_129*tmp_76 - tmp_129*tmp_80 - tmp_129*tmp_84 - tmp_129*tmp_88 - tmp_129*tmp_92 - tmp_129*tmp_96;
      real_t a_0_2 = -tmp_100*tmp_130 - tmp_104*tmp_130 - tmp_108*tmp_130 - tmp_112*tmp_130 - tmp_116*tmp_130 - tmp_120*tmp_130 - tmp_124*tmp_130 - tmp_128*tmp_130 - tmp_130*tmp_48 - tmp_130*tmp_52 - tmp_130*tmp_56 - tmp_130*tmp_60 - tmp_130*tmp_64 - tmp_130*tmp_68 - tmp_130*tmp_72 - tmp_130*tmp_76 - tmp_130*tmp_80 - tmp_130*tmp_84 - tmp_130*tmp_88 - tmp_130*tmp_92 - tmp_130*tmp_96;
      real_t a_0_3 = -tmp_100*tmp_131 - tmp_104*tmp_131 - tmp_108*tmp_131 - tmp_112*tmp_131 - tmp_116*tmp_131 - tmp_120*tmp_131 - tmp_124*tmp_131 - tmp_128*tmp_131 - tmp_131*tmp_48 - tmp_131*tmp_52 - tmp_131*tmp_56 - tmp_131*tmp_60 - tmp_131*tmp_64 - tmp_131*tmp_68 - tmp_131*tmp_72 - tmp_131*tmp_76 - tmp_131*tmp_80 - tmp_131*tmp_84 - tmp_131*tmp_88 - tmp_131*tmp_92 - tmp_131*tmp_96;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }


};




class EGNIPGVectorLaplaceFormEP1_1 : public hyteg::dg::DGForm
{
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
      real_t tmp_28 = 5/tmp_27;
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
      real_t tmp_38 = 5/tmp_37;
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

      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_0;
      real_t tmp_3 = p_affine_1_0 + tmp_2;
      real_t tmp_4 = p_affine_1_1 + tmp_0;
      real_t tmp_5 = 1.0 / (tmp_1*tmp_3 - tmp_4*(p_affine_2_0 + tmp_2));
      real_t tmp_6 = tmp_1*tmp_5;
      real_t tmp_7 = tmp_5*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_8 = tmp_3*tmp_5;
      real_t tmp_9 = tmp_5*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_10 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_11 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_12 = std::abs(std::pow((tmp_10*tmp_10) + (tmp_11*tmp_11), 1.0/2.0));
      real_t tmp_13 = tmp_12*(p_affine_10_0*(-tmp_6 - tmp_7) + p_affine_10_1*(-tmp_8 - tmp_9));
      real_t tmp_14 = p_affine_6_1 + tmp_0;
      real_t tmp_15 = 0.046910077030668018*tmp_11 + tmp_14;
      real_t tmp_16 = p_affine_6_0 + tmp_2;
      real_t tmp_17 = 0.046910077030668018*tmp_10 + tmp_16;
      real_t tmp_18 = 0.11846344252809471*tmp_1*(tmp_15*tmp_8 + tmp_17*tmp_7 - 1.0/3.0) + 0.11846344252809471*tmp_4*(tmp_15*tmp_9 + tmp_17*tmp_6 - 1.0/3.0);
      real_t tmp_19 = 0.23076534494715845*tmp_11 + tmp_14;
      real_t tmp_20 = 0.23076534494715845*tmp_10 + tmp_16;
      real_t tmp_21 = 0.2393143352496831*tmp_1*(tmp_19*tmp_8 + tmp_20*tmp_7 - 1.0/3.0) + 0.2393143352496831*tmp_4*(tmp_19*tmp_9 + tmp_20*tmp_6 - 1.0/3.0);
      real_t tmp_22 = 0.5*tmp_11 + tmp_14;
      real_t tmp_23 = 0.5*tmp_10 + tmp_16;
      real_t tmp_24 = 0.2844444444444445*tmp_1*(tmp_22*tmp_8 + tmp_23*tmp_7 - 1.0/3.0) + 0.2844444444444445*tmp_4*(tmp_22*tmp_9 + tmp_23*tmp_6 - 1.0/3.0);
      real_t tmp_25 = 0.7692346550528415*tmp_11 + tmp_14;
      real_t tmp_26 = 0.7692346550528415*tmp_10 + tmp_16;
      real_t tmp_27 = 0.2393143352496831*tmp_1*(tmp_25*tmp_8 + tmp_26*tmp_7 - 1.0/3.0) + 0.2393143352496831*tmp_4*(tmp_25*tmp_9 + tmp_26*tmp_6 - 1.0/3.0);
      real_t tmp_28 = 0.95308992296933193*tmp_11 + tmp_14;
      real_t tmp_29 = 0.95308992296933193*tmp_10 + tmp_16;
      real_t tmp_30 = 0.11846344252809471*tmp_1*(tmp_28*tmp_8 + tmp_29*tmp_7 - 1.0/3.0) + 0.11846344252809471*tmp_4*(tmp_28*tmp_9 + tmp_29*tmp_6 - 1.0/3.0);
      real_t tmp_31 = tmp_12*(p_affine_10_0*tmp_6 + p_affine_10_1*tmp_9);
      real_t tmp_32 = tmp_12*(p_affine_10_0*tmp_7 + p_affine_10_1*tmp_8);
      real_t a_0_0 = -tmp_13*tmp_18 - tmp_13*tmp_21 - tmp_13*tmp_24 - tmp_13*tmp_27 - tmp_13*tmp_30;
      real_t a_0_1 = -tmp_18*tmp_31 - tmp_21*tmp_31 - tmp_24*tmp_31 - tmp_27*tmp_31 - tmp_30*tmp_31;
      real_t a_0_2 = -tmp_18*tmp_32 - tmp_21*tmp_32 - tmp_24*tmp_32 - tmp_27*tmp_32 - tmp_30*tmp_32;
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
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
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
      real_t tmp_22 = tmp_11*tmp_19 + tmp_20*tmp_3 + tmp_21*tmp_6;
      real_t tmp_23 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_5);
      real_t tmp_24 = tmp_18*(tmp_1*tmp_9 - tmp_10*tmp_14);
      real_t tmp_25 = tmp_18*(tmp_13*tmp_14 - tmp_5*tmp_9);
      real_t tmp_26 = tmp_11*tmp_23 + tmp_24*tmp_3 + tmp_25*tmp_6;
      real_t tmp_27 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_6);
      real_t tmp_28 = tmp_18*(tmp_10*tmp_11 - tmp_6*tmp_9);
      real_t tmp_29 = tmp_18*(-tmp_11*tmp_13 + tmp_3*tmp_9);
      real_t tmp_30 = tmp_11*tmp_27 + tmp_28*tmp_3 + tmp_29*tmp_6;
      real_t tmp_31 = p_affine_0_0*p_affine_1_1;
      real_t tmp_32 = p_affine_0_0*p_affine_1_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_2;
      real_t tmp_34 = p_affine_0_1*p_affine_1_0;
      real_t tmp_35 = p_affine_0_1*p_affine_1_2;
      real_t tmp_36 = p_affine_2_2*p_affine_3_0;
      real_t tmp_37 = p_affine_0_2*p_affine_1_0;
      real_t tmp_38 = p_affine_0_2*p_affine_1_1;
      real_t tmp_39 = p_affine_2_0*p_affine_3_1;
      real_t tmp_40 = p_affine_2_2*p_affine_3_1;
      real_t tmp_41 = p_affine_2_0*p_affine_3_2;
      real_t tmp_42 = p_affine_2_1*p_affine_3_0;
      real_t tmp_43 = std::abs(p_affine_0_0*tmp_33 - p_affine_0_0*tmp_40 + p_affine_0_1*tmp_36 - p_affine_0_1*tmp_41 + p_affine_0_2*tmp_39 - p_affine_0_2*tmp_42 - p_affine_1_0*tmp_33 + p_affine_1_0*tmp_40 - p_affine_1_1*tmp_36 + p_affine_1_1*tmp_41 - p_affine_1_2*tmp_39 + p_affine_1_2*tmp_42 + p_affine_2_0*tmp_35 - p_affine_2_0*tmp_38 - p_affine_2_1*tmp_32 + p_affine_2_1*tmp_37 + p_affine_2_2*tmp_31 - p_affine_2_2*tmp_34 - p_affine_3_0*tmp_35 + p_affine_3_0*tmp_38 + p_affine_3_1*tmp_32 - p_affine_3_1*tmp_37 - p_affine_3_2*tmp_31 + p_affine_3_2*tmp_34);
      real_t tmp_44 = tmp_43*(tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25) + tmp_30*(-tmp_27 - tmp_28 - tmp_29));
      real_t tmp_45 = tmp_43*(tmp_21*tmp_22 + tmp_25*tmp_26 + tmp_29*tmp_30);
      real_t tmp_46 = tmp_43*(tmp_20*tmp_22 + tmp_24*tmp_26 + tmp_28*tmp_30);
      real_t tmp_47 = tmp_43*(tmp_19*tmp_22 + tmp_23*tmp_26 + tmp_27*tmp_30);
      real_t a_0_0 = 0.1666666666666668*tmp_44;
      real_t a_0_1 = 0.1666666666666668*tmp_45;
      real_t a_0_2 = 0.1666666666666668*tmp_46;
      real_t a_0_3 = 0.1666666666666668*tmp_47;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
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

         real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_1_1 + tmp_0;
      real_t tmp_2 = -p_affine_8_2;
      real_t tmp_3 = p_affine_9_2 + tmp_2;
      real_t tmp_4 = p_affine_10_2 + tmp_2;
      real_t tmp_5 = -p_affine_0_2;
      real_t tmp_6 = p_affine_8_2 + tmp_5;
      real_t tmp_7 = 0.031405749086161582*tmp_3 + 0.93718850182767688*tmp_4 + tmp_6;
      real_t tmp_8 = -p_affine_0_0;
      real_t tmp_9 = p_affine_2_0 + tmp_8;
      real_t tmp_10 = p_affine_3_1 + tmp_0;
      real_t tmp_11 = p_affine_3_0 + tmp_8;
      real_t tmp_12 = p_affine_2_1 + tmp_0;
      real_t tmp_13 = p_affine_1_0 + tmp_8;
      real_t tmp_14 = p_affine_3_2 + tmp_5;
      real_t tmp_15 = tmp_12*tmp_14;
      real_t tmp_16 = p_affine_1_2 + tmp_5;
      real_t tmp_17 = tmp_10*tmp_16;
      real_t tmp_18 = p_affine_2_2 + tmp_5;
      real_t tmp_19 = tmp_1*tmp_18;
      real_t tmp_20 = tmp_10*tmp_18;
      real_t tmp_21 = tmp_1*tmp_14;
      real_t tmp_22 = tmp_12*tmp_16;
      real_t tmp_23 = 1.0 / (tmp_11*tmp_19 - tmp_11*tmp_22 + tmp_13*tmp_15 - tmp_13*tmp_20 + tmp_17*tmp_9 - tmp_21*tmp_9);
      real_t tmp_24 = tmp_23*(tmp_10*tmp_9 - tmp_11*tmp_12);
      real_t tmp_25 = tmp_24*tmp_7;
      real_t tmp_26 = -p_affine_8_1;
      real_t tmp_27 = p_affine_9_1 + tmp_26;
      real_t tmp_28 = p_affine_10_1 + tmp_26;
      real_t tmp_29 = p_affine_8_1 + tmp_0;
      real_t tmp_30 = 0.031405749086161582*tmp_27 + 0.93718850182767688*tmp_28 + tmp_29;
      real_t tmp_31 = tmp_23*(tmp_11*tmp_18 - tmp_14*tmp_9);
      real_t tmp_32 = tmp_30*tmp_31;
      real_t tmp_33 = -p_affine_8_0;
      real_t tmp_34 = p_affine_9_0 + tmp_33;
      real_t tmp_35 = p_affine_10_0 + tmp_33;
      real_t tmp_36 = p_affine_8_0 + tmp_8;
      real_t tmp_37 = 0.031405749086161582*tmp_34 + 0.93718850182767688*tmp_35 + tmp_36;
      real_t tmp_38 = tmp_23*(tmp_15 - tmp_20);
      real_t tmp_39 = tmp_37*tmp_38;
      real_t tmp_40 = tmp_25 + tmp_32 + tmp_39;
      real_t tmp_41 = tmp_23*(tmp_1*tmp_11 - tmp_10*tmp_13);
      real_t tmp_42 = tmp_41*tmp_7;
      real_t tmp_43 = tmp_23*(-tmp_11*tmp_16 + tmp_13*tmp_14);
      real_t tmp_44 = tmp_30*tmp_43;
      real_t tmp_45 = tmp_23*(tmp_17 - tmp_21);
      real_t tmp_46 = tmp_37*tmp_45;
      real_t tmp_47 = tmp_42 + tmp_44 + tmp_46;
      real_t tmp_48 = tmp_23*(-tmp_1*tmp_9 + tmp_12*tmp_13);
      real_t tmp_49 = tmp_48*tmp_7;
      real_t tmp_50 = tmp_23*(-tmp_13*tmp_18 + tmp_16*tmp_9);
      real_t tmp_51 = tmp_30*tmp_50;
      real_t tmp_52 = tmp_23*(tmp_19 - tmp_22);
      real_t tmp_53 = tmp_37*tmp_52;
      real_t tmp_54 = tmp_49 + tmp_51 + tmp_53;
      real_t tmp_55 = tmp_1*(tmp_40 - 1.0/4.0) + tmp_10*(tmp_54 - 1.0/4.0) + tmp_12*(tmp_47 - 1.0/4.0);
      real_t tmp_56 = 0.5*p_affine_13_0*(-tmp_38 - tmp_45 - tmp_52) + 0.5*p_affine_13_1*(-tmp_31 - tmp_43 - tmp_50) + 0.5*p_affine_13_2*(-tmp_24 - tmp_41 - tmp_48);
      real_t tmp_57 = -tmp_25 - tmp_32 - tmp_39 - tmp_42 - tmp_44 - tmp_46 - tmp_49 - tmp_51 - tmp_53 + 1;
      real_t tmp_58 = 0.5*p_affine_13_0*(tmp_1*tmp_38 + tmp_10*tmp_52 + tmp_12*tmp_45) + 0.5*p_affine_13_1*(tmp_1*tmp_31 + tmp_10*tmp_50 + tmp_12*tmp_43) + 0.5*p_affine_13_2*(tmp_1*tmp_24 + tmp_10*tmp_48 + tmp_12*tmp_41);
      real_t tmp_59 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_60 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_61 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_62 = (std::abs(tmp_28*tmp_60 - tmp_35*tmp_59)*std::abs(tmp_28*tmp_60 - tmp_35*tmp_59)) + (std::abs(tmp_28*tmp_61 - tmp_4*tmp_59)*std::abs(tmp_28*tmp_61 - tmp_4*tmp_59)) + (std::abs(tmp_35*tmp_61 - tmp_4*tmp_60)*std::abs(tmp_35*tmp_61 - tmp_4*tmp_60));
      real_t tmp_63 = 5.0*std::pow(tmp_62, -0.25);
      real_t tmp_64 = tmp_55*tmp_63;
      real_t tmp_65 = 1.0*std::pow(tmp_62, 1.0/2.0);
      real_t tmp_66 = 0.0068572537431980923*tmp_65;
      real_t tmp_67 = 0.19601935860219369*tmp_3 + 0.60796128279561268*tmp_4 + tmp_6;
      real_t tmp_68 = tmp_24*tmp_67;
      real_t tmp_69 = 0.19601935860219369*tmp_27 + 0.60796128279561268*tmp_28 + tmp_29;
      real_t tmp_70 = tmp_31*tmp_69;
      real_t tmp_71 = 0.19601935860219369*tmp_34 + 0.60796128279561268*tmp_35 + tmp_36;
      real_t tmp_72 = tmp_38*tmp_71;
      real_t tmp_73 = tmp_68 + tmp_70 + tmp_72;
      real_t tmp_74 = tmp_41*tmp_67;
      real_t tmp_75 = tmp_43*tmp_69;
      real_t tmp_76 = tmp_45*tmp_71;
      real_t tmp_77 = tmp_74 + tmp_75 + tmp_76;
      real_t tmp_78 = tmp_48*tmp_67;
      real_t tmp_79 = tmp_50*tmp_69;
      real_t tmp_80 = tmp_52*tmp_71;
      real_t tmp_81 = tmp_78 + tmp_79 + tmp_80;
      real_t tmp_82 = tmp_1*(tmp_73 - 1.0/4.0) + tmp_10*(tmp_81 - 1.0/4.0) + tmp_12*(tmp_77 - 1.0/4.0);
      real_t tmp_83 = -tmp_68 - tmp_70 - tmp_72 - tmp_74 - tmp_75 - tmp_76 - tmp_78 - tmp_79 - tmp_80 + 1;
      real_t tmp_84 = tmp_63*tmp_82;
      real_t tmp_85 = 0.037198804536718075*tmp_65;
      real_t tmp_86 = 0.37605877282253791*tmp_3 + 0.039308471900058539*tmp_4 + tmp_6;
      real_t tmp_87 = tmp_24*tmp_86;
      real_t tmp_88 = 0.37605877282253791*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29;
      real_t tmp_89 = tmp_31*tmp_88;
      real_t tmp_90 = 0.37605877282253791*tmp_34 + 0.039308471900058539*tmp_35 + tmp_36;
      real_t tmp_91 = tmp_38*tmp_90;
      real_t tmp_92 = tmp_87 + tmp_89 + tmp_91;
      real_t tmp_93 = tmp_41*tmp_86;
      real_t tmp_94 = tmp_43*tmp_88;
      real_t tmp_95 = tmp_45*tmp_90;
      real_t tmp_96 = tmp_93 + tmp_94 + tmp_95;
      real_t tmp_97 = tmp_48*tmp_86;
      real_t tmp_98 = tmp_50*tmp_88;
      real_t tmp_99 = tmp_52*tmp_90;
      real_t tmp_100 = tmp_97 + tmp_98 + tmp_99;
      real_t tmp_101 = tmp_1*(tmp_92 - 1.0/4.0) + tmp_10*(tmp_100 - 1.0/4.0) + tmp_12*(tmp_96 - 1.0/4.0);
      real_t tmp_102 = -tmp_87 - tmp_89 - tmp_91 - tmp_93 - tmp_94 - tmp_95 - tmp_97 - tmp_98 - tmp_99 + 1;
      real_t tmp_103 = tmp_101*tmp_63;
      real_t tmp_104 = 0.020848748529055869*tmp_65;
      real_t tmp_105 = 0.78764240869137092*tmp_3 + 0.1711304259088916*tmp_4 + tmp_6;
      real_t tmp_106 = tmp_105*tmp_24;
      real_t tmp_107 = 0.78764240869137092*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29;
      real_t tmp_108 = tmp_107*tmp_31;
      real_t tmp_109 = 0.78764240869137092*tmp_34 + 0.1711304259088916*tmp_35 + tmp_36;
      real_t tmp_110 = tmp_109*tmp_38;
      real_t tmp_111 = tmp_106 + tmp_108 + tmp_110;
      real_t tmp_112 = tmp_105*tmp_41;
      real_t tmp_113 = tmp_107*tmp_43;
      real_t tmp_114 = tmp_109*tmp_45;
      real_t tmp_115 = tmp_112 + tmp_113 + tmp_114;
      real_t tmp_116 = tmp_105*tmp_48;
      real_t tmp_117 = tmp_107*tmp_50;
      real_t tmp_118 = tmp_109*tmp_52;
      real_t tmp_119 = tmp_116 + tmp_117 + tmp_118;
      real_t tmp_120 = tmp_1*(tmp_111 - 1.0/4.0) + tmp_10*(tmp_119 - 1.0/4.0) + tmp_12*(tmp_115 - 1.0/4.0);
      real_t tmp_121 = -tmp_106 - tmp_108 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_116 - tmp_117 - tmp_118 + 1;
      real_t tmp_122 = tmp_120*tmp_63;
      real_t tmp_123 = 0.019202922745021479*tmp_65;
      real_t tmp_124 = 0.58463275527740355*tmp_3 + 0.37605877282253791*tmp_4 + tmp_6;
      real_t tmp_125 = tmp_124*tmp_24;
      real_t tmp_126 = 0.58463275527740355*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29;
      real_t tmp_127 = tmp_126*tmp_31;
      real_t tmp_128 = 0.58463275527740355*tmp_34 + 0.37605877282253791*tmp_35 + tmp_36;
      real_t tmp_129 = tmp_128*tmp_38;
      real_t tmp_130 = tmp_125 + tmp_127 + tmp_129;
      real_t tmp_131 = tmp_124*tmp_41;
      real_t tmp_132 = tmp_126*tmp_43;
      real_t tmp_133 = tmp_128*tmp_45;
      real_t tmp_134 = tmp_131 + tmp_132 + tmp_133;
      real_t tmp_135 = tmp_124*tmp_48;
      real_t tmp_136 = tmp_126*tmp_50;
      real_t tmp_137 = tmp_128*tmp_52;
      real_t tmp_138 = tmp_135 + tmp_136 + tmp_137;
      real_t tmp_139 = tmp_1*(tmp_130 - 1.0/4.0) + tmp_10*(tmp_138 - 1.0/4.0) + tmp_12*(tmp_134 - 1.0/4.0);
      real_t tmp_140 = -tmp_125 - tmp_127 - tmp_129 - tmp_131 - tmp_132 - tmp_133 - tmp_135 - tmp_136 - tmp_137 + 1;
      real_t tmp_141 = tmp_139*tmp_63;
      real_t tmp_142 = 0.020848748529055869*tmp_65;
      real_t tmp_143 = 0.041227165399737475*tmp_3 + 0.78764240869137092*tmp_4 + tmp_6;
      real_t tmp_144 = tmp_143*tmp_24;
      real_t tmp_145 = 0.041227165399737475*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29;
      real_t tmp_146 = tmp_145*tmp_31;
      real_t tmp_147 = 0.041227165399737475*tmp_34 + 0.78764240869137092*tmp_35 + tmp_36;
      real_t tmp_148 = tmp_147*tmp_38;
      real_t tmp_149 = tmp_144 + tmp_146 + tmp_148;
      real_t tmp_150 = tmp_143*tmp_41;
      real_t tmp_151 = tmp_145*tmp_43;
      real_t tmp_152 = tmp_147*tmp_45;
      real_t tmp_153 = tmp_150 + tmp_151 + tmp_152;
      real_t tmp_154 = tmp_143*tmp_48;
      real_t tmp_155 = tmp_145*tmp_50;
      real_t tmp_156 = tmp_147*tmp_52;
      real_t tmp_157 = tmp_154 + tmp_155 + tmp_156;
      real_t tmp_158 = tmp_1*(tmp_149 - 1.0/4.0) + tmp_10*(tmp_157 - 1.0/4.0) + tmp_12*(tmp_153 - 1.0/4.0);
      real_t tmp_159 = -tmp_144 - tmp_146 - tmp_148 - tmp_150 - tmp_151 - tmp_152 - tmp_154 - tmp_155 - tmp_156 + 1;
      real_t tmp_160 = tmp_158*tmp_63;
      real_t tmp_161 = 0.019202922745021479*tmp_65;
      real_t tmp_162 = 0.039308471900058539*tmp_3 + 0.58463275527740355*tmp_4 + tmp_6;
      real_t tmp_163 = tmp_162*tmp_24;
      real_t tmp_164 = 0.039308471900058539*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29;
      real_t tmp_165 = tmp_164*tmp_31;
      real_t tmp_166 = 0.039308471900058539*tmp_34 + 0.58463275527740355*tmp_35 + tmp_36;
      real_t tmp_167 = tmp_166*tmp_38;
      real_t tmp_168 = tmp_163 + tmp_165 + tmp_167;
      real_t tmp_169 = tmp_162*tmp_41;
      real_t tmp_170 = tmp_164*tmp_43;
      real_t tmp_171 = tmp_166*tmp_45;
      real_t tmp_172 = tmp_169 + tmp_170 + tmp_171;
      real_t tmp_173 = tmp_162*tmp_48;
      real_t tmp_174 = tmp_164*tmp_50;
      real_t tmp_175 = tmp_166*tmp_52;
      real_t tmp_176 = tmp_173 + tmp_174 + tmp_175;
      real_t tmp_177 = tmp_1*(tmp_168 - 1.0/4.0) + tmp_10*(tmp_176 - 1.0/4.0) + tmp_12*(tmp_172 - 1.0/4.0);
      real_t tmp_178 = -tmp_163 - tmp_165 - tmp_167 - tmp_169 - tmp_170 - tmp_171 - tmp_173 - tmp_174 - tmp_175 + 1;
      real_t tmp_179 = tmp_177*tmp_63;
      real_t tmp_180 = 0.020848748529055869*tmp_65;
      real_t tmp_181 = 0.78764240869137092*tmp_3 + 0.041227165399737475*tmp_4 + tmp_6;
      real_t tmp_182 = tmp_181*tmp_24;
      real_t tmp_183 = 0.78764240869137092*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29;
      real_t tmp_184 = tmp_183*tmp_31;
      real_t tmp_185 = 0.78764240869137092*tmp_34 + 0.041227165399737475*tmp_35 + tmp_36;
      real_t tmp_186 = tmp_185*tmp_38;
      real_t tmp_187 = tmp_182 + tmp_184 + tmp_186;
      real_t tmp_188 = tmp_181*tmp_41;
      real_t tmp_189 = tmp_183*tmp_43;
      real_t tmp_190 = tmp_185*tmp_45;
      real_t tmp_191 = tmp_188 + tmp_189 + tmp_190;
      real_t tmp_192 = tmp_181*tmp_48;
      real_t tmp_193 = tmp_183*tmp_50;
      real_t tmp_194 = tmp_185*tmp_52;
      real_t tmp_195 = tmp_192 + tmp_193 + tmp_194;
      real_t tmp_196 = tmp_1*(tmp_187 - 1.0/4.0) + tmp_10*(tmp_195 - 1.0/4.0) + tmp_12*(tmp_191 - 1.0/4.0);
      real_t tmp_197 = -tmp_182 - tmp_184 - tmp_186 - tmp_188 - tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 + 1;
      real_t tmp_198 = tmp_196*tmp_63;
      real_t tmp_199 = 0.019202922745021479*tmp_65;
      real_t tmp_200 = 0.58463275527740355*tmp_3 + 0.039308471900058539*tmp_4 + tmp_6;
      real_t tmp_201 = tmp_200*tmp_24;
      real_t tmp_202 = 0.58463275527740355*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29;
      real_t tmp_203 = tmp_202*tmp_31;
      real_t tmp_204 = 0.58463275527740355*tmp_34 + 0.039308471900058539*tmp_35 + tmp_36;
      real_t tmp_205 = tmp_204*tmp_38;
      real_t tmp_206 = tmp_201 + tmp_203 + tmp_205;
      real_t tmp_207 = tmp_200*tmp_41;
      real_t tmp_208 = tmp_202*tmp_43;
      real_t tmp_209 = tmp_204*tmp_45;
      real_t tmp_210 = tmp_207 + tmp_208 + tmp_209;
      real_t tmp_211 = tmp_200*tmp_48;
      real_t tmp_212 = tmp_202*tmp_50;
      real_t tmp_213 = tmp_204*tmp_52;
      real_t tmp_214 = tmp_211 + tmp_212 + tmp_213;
      real_t tmp_215 = tmp_1*(tmp_206 - 1.0/4.0) + tmp_10*(tmp_214 - 1.0/4.0) + tmp_12*(tmp_210 - 1.0/4.0);
      real_t tmp_216 = -tmp_201 - tmp_203 - tmp_205 - tmp_207 - tmp_208 - tmp_209 - tmp_211 - tmp_212 - tmp_213 + 1;
      real_t tmp_217 = tmp_215*tmp_63;
      real_t tmp_218 = 0.020848748529055869*tmp_65;
      real_t tmp_219 = 0.1711304259088916*tmp_3 + 0.78764240869137092*tmp_4 + tmp_6;
      real_t tmp_220 = tmp_219*tmp_24;
      real_t tmp_221 = 0.1711304259088916*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29;
      real_t tmp_222 = tmp_221*tmp_31;
      real_t tmp_223 = 0.1711304259088916*tmp_34 + 0.78764240869137092*tmp_35 + tmp_36;
      real_t tmp_224 = tmp_223*tmp_38;
      real_t tmp_225 = tmp_220 + tmp_222 + tmp_224;
      real_t tmp_226 = tmp_219*tmp_41;
      real_t tmp_227 = tmp_221*tmp_43;
      real_t tmp_228 = tmp_223*tmp_45;
      real_t tmp_229 = tmp_226 + tmp_227 + tmp_228;
      real_t tmp_230 = tmp_219*tmp_48;
      real_t tmp_231 = tmp_221*tmp_50;
      real_t tmp_232 = tmp_223*tmp_52;
      real_t tmp_233 = tmp_230 + tmp_231 + tmp_232;
      real_t tmp_234 = tmp_1*(tmp_225 - 1.0/4.0) + tmp_10*(tmp_233 - 1.0/4.0) + tmp_12*(tmp_229 - 1.0/4.0);
      real_t tmp_235 = -tmp_220 - tmp_222 - tmp_224 - tmp_226 - tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 + 1;
      real_t tmp_236 = tmp_234*tmp_63;
      real_t tmp_237 = 0.019202922745021479*tmp_65;
      real_t tmp_238 = 0.37605877282253791*tmp_3 + 0.58463275527740355*tmp_4 + tmp_6;
      real_t tmp_239 = tmp_238*tmp_24;
      real_t tmp_240 = 0.37605877282253791*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29;
      real_t tmp_241 = tmp_240*tmp_31;
      real_t tmp_242 = 0.37605877282253791*tmp_34 + 0.58463275527740355*tmp_35 + tmp_36;
      real_t tmp_243 = tmp_242*tmp_38;
      real_t tmp_244 = tmp_239 + tmp_241 + tmp_243;
      real_t tmp_245 = tmp_238*tmp_41;
      real_t tmp_246 = tmp_240*tmp_43;
      real_t tmp_247 = tmp_242*tmp_45;
      real_t tmp_248 = tmp_245 + tmp_246 + tmp_247;
      real_t tmp_249 = tmp_238*tmp_48;
      real_t tmp_250 = tmp_240*tmp_50;
      real_t tmp_251 = tmp_242*tmp_52;
      real_t tmp_252 = tmp_249 + tmp_250 + tmp_251;
      real_t tmp_253 = tmp_1*(tmp_244 - 1.0/4.0) + tmp_10*(tmp_252 - 1.0/4.0) + tmp_12*(tmp_248 - 1.0/4.0);
      real_t tmp_254 = -tmp_239 - tmp_241 - tmp_243 - tmp_245 - tmp_246 - tmp_247 - tmp_249 - tmp_250 - tmp_251 + 1;
      real_t tmp_255 = tmp_253*tmp_63;
      real_t tmp_256 = 0.020848748529055869*tmp_65;
      real_t tmp_257 = 0.041227165399737475*tmp_3 + 0.1711304259088916*tmp_4 + tmp_6;
      real_t tmp_258 = tmp_24*tmp_257;
      real_t tmp_259 = 0.041227165399737475*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29;
      real_t tmp_260 = tmp_259*tmp_31;
      real_t tmp_261 = 0.041227165399737475*tmp_34 + 0.1711304259088916*tmp_35 + tmp_36;
      real_t tmp_262 = tmp_261*tmp_38;
      real_t tmp_263 = tmp_258 + tmp_260 + tmp_262;
      real_t tmp_264 = tmp_257*tmp_41;
      real_t tmp_265 = tmp_259*tmp_43;
      real_t tmp_266 = tmp_261*tmp_45;
      real_t tmp_267 = tmp_264 + tmp_265 + tmp_266;
      real_t tmp_268 = tmp_257*tmp_48;
      real_t tmp_269 = tmp_259*tmp_50;
      real_t tmp_270 = tmp_261*tmp_52;
      real_t tmp_271 = tmp_268 + tmp_269 + tmp_270;
      real_t tmp_272 = tmp_1*(tmp_263 - 1.0/4.0) + tmp_10*(tmp_271 - 1.0/4.0) + tmp_12*(tmp_267 - 1.0/4.0);
      real_t tmp_273 = -tmp_258 - tmp_260 - tmp_262 - tmp_264 - tmp_265 - tmp_266 - tmp_268 - tmp_269 - tmp_270 + 1;
      real_t tmp_274 = tmp_272*tmp_63;
      real_t tmp_275 = 0.019202922745021479*tmp_65;
      real_t tmp_276 = 0.40446199974765351*tmp_3 + 0.19107600050469298*tmp_4 + tmp_6;
      real_t tmp_277 = tmp_24*tmp_276;
      real_t tmp_278 = 0.40446199974765351*tmp_27 + 0.19107600050469298*tmp_28 + tmp_29;
      real_t tmp_279 = tmp_278*tmp_31;
      real_t tmp_280 = 0.40446199974765351*tmp_34 + 0.19107600050469298*tmp_35 + tmp_36;
      real_t tmp_281 = tmp_280*tmp_38;
      real_t tmp_282 = tmp_277 + tmp_279 + tmp_281;
      real_t tmp_283 = tmp_276*tmp_41;
      real_t tmp_284 = tmp_278*tmp_43;
      real_t tmp_285 = tmp_280*tmp_45;
      real_t tmp_286 = tmp_283 + tmp_284 + tmp_285;
      real_t tmp_287 = tmp_276*tmp_48;
      real_t tmp_288 = tmp_278*tmp_50;
      real_t tmp_289 = tmp_280*tmp_52;
      real_t tmp_290 = tmp_287 + tmp_288 + tmp_289;
      real_t tmp_291 = tmp_1*(tmp_282 - 1.0/4.0) + tmp_10*(tmp_290 - 1.0/4.0) + tmp_12*(tmp_286 - 1.0/4.0);
      real_t tmp_292 = -tmp_277 - tmp_279 - tmp_281 - tmp_283 - tmp_284 - tmp_285 - tmp_287 - tmp_288 - tmp_289 + 1;
      real_t tmp_293 = tmp_291*tmp_63;
      real_t tmp_294 = 0.042507265838595799*tmp_65;
      real_t tmp_295 = 0.039308471900058539*tmp_3 + 0.37605877282253791*tmp_4 + tmp_6;
      real_t tmp_296 = tmp_24*tmp_295;
      real_t tmp_297 = 0.039308471900058539*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29;
      real_t tmp_298 = tmp_297*tmp_31;
      real_t tmp_299 = 0.039308471900058539*tmp_34 + 0.37605877282253791*tmp_35 + tmp_36;
      real_t tmp_300 = tmp_299*tmp_38;
      real_t tmp_301 = tmp_296 + tmp_298 + tmp_300;
      real_t tmp_302 = tmp_295*tmp_41;
      real_t tmp_303 = tmp_297*tmp_43;
      real_t tmp_304 = tmp_299*tmp_45;
      real_t tmp_305 = tmp_302 + tmp_303 + tmp_304;
      real_t tmp_306 = tmp_295*tmp_48;
      real_t tmp_307 = tmp_297*tmp_50;
      real_t tmp_308 = tmp_299*tmp_52;
      real_t tmp_309 = tmp_306 + tmp_307 + tmp_308;
      real_t tmp_310 = tmp_1*(tmp_301 - 1.0/4.0) + tmp_10*(tmp_309 - 1.0/4.0) + tmp_12*(tmp_305 - 1.0/4.0);
      real_t tmp_311 = -tmp_296 - tmp_298 - tmp_300 - tmp_302 - tmp_303 - tmp_304 - tmp_306 - tmp_307 - tmp_308 + 1;
      real_t tmp_312 = tmp_310*tmp_63;
      real_t tmp_313 = 0.020848748529055869*tmp_65;
      real_t tmp_314 = 0.93718850182767688*tmp_3 + 0.031405749086161582*tmp_4 + tmp_6;
      real_t tmp_315 = tmp_24*tmp_314;
      real_t tmp_316 = 0.93718850182767688*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29;
      real_t tmp_317 = tmp_31*tmp_316;
      real_t tmp_318 = 0.93718850182767688*tmp_34 + 0.031405749086161582*tmp_35 + tmp_36;
      real_t tmp_319 = tmp_318*tmp_38;
      real_t tmp_320 = tmp_315 + tmp_317 + tmp_319;
      real_t tmp_321 = tmp_314*tmp_41;
      real_t tmp_322 = tmp_316*tmp_43;
      real_t tmp_323 = tmp_318*tmp_45;
      real_t tmp_324 = tmp_321 + tmp_322 + tmp_323;
      real_t tmp_325 = tmp_314*tmp_48;
      real_t tmp_326 = tmp_316*tmp_50;
      real_t tmp_327 = tmp_318*tmp_52;
      real_t tmp_328 = tmp_325 + tmp_326 + tmp_327;
      real_t tmp_329 = tmp_1*(tmp_320 - 1.0/4.0) + tmp_10*(tmp_328 - 1.0/4.0) + tmp_12*(tmp_324 - 1.0/4.0);
      real_t tmp_330 = -tmp_315 - tmp_317 - tmp_319 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1;
      real_t tmp_331 = tmp_329*tmp_63;
      real_t tmp_332 = 0.0068572537431980923*tmp_65;
      real_t tmp_333 = 0.60796128279561268*tmp_3 + 0.19601935860219369*tmp_4 + tmp_6;
      real_t tmp_334 = tmp_24*tmp_333;
      real_t tmp_335 = 0.60796128279561268*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29;
      real_t tmp_336 = tmp_31*tmp_335;
      real_t tmp_337 = 0.60796128279561268*tmp_34 + 0.19601935860219369*tmp_35 + tmp_36;
      real_t tmp_338 = tmp_337*tmp_38;
      real_t tmp_339 = tmp_334 + tmp_336 + tmp_338;
      real_t tmp_340 = tmp_333*tmp_41;
      real_t tmp_341 = tmp_335*tmp_43;
      real_t tmp_342 = tmp_337*tmp_45;
      real_t tmp_343 = tmp_340 + tmp_341 + tmp_342;
      real_t tmp_344 = tmp_333*tmp_48;
      real_t tmp_345 = tmp_335*tmp_50;
      real_t tmp_346 = tmp_337*tmp_52;
      real_t tmp_347 = tmp_344 + tmp_345 + tmp_346;
      real_t tmp_348 = tmp_1*(tmp_339 - 1.0/4.0) + tmp_10*(tmp_347 - 1.0/4.0) + tmp_12*(tmp_343 - 1.0/4.0);
      real_t tmp_349 = -tmp_334 - tmp_336 - tmp_338 - tmp_340 - tmp_341 - tmp_342 - tmp_344 - tmp_345 - tmp_346 + 1;
      real_t tmp_350 = tmp_348*tmp_63;
      real_t tmp_351 = 0.037198804536718075*tmp_65;
      real_t tmp_352 = 0.19107600050469298*tmp_3 + 0.40446199974765351*tmp_4 + tmp_6;
      real_t tmp_353 = tmp_24*tmp_352;
      real_t tmp_354 = 0.19107600050469298*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29;
      real_t tmp_355 = tmp_31*tmp_354;
      real_t tmp_356 = 0.19107600050469298*tmp_34 + 0.40446199974765351*tmp_35 + tmp_36;
      real_t tmp_357 = tmp_356*tmp_38;
      real_t tmp_358 = tmp_353 + tmp_355 + tmp_357;
      real_t tmp_359 = tmp_352*tmp_41;
      real_t tmp_360 = tmp_354*tmp_43;
      real_t tmp_361 = tmp_356*tmp_45;
      real_t tmp_362 = tmp_359 + tmp_360 + tmp_361;
      real_t tmp_363 = tmp_352*tmp_48;
      real_t tmp_364 = tmp_354*tmp_50;
      real_t tmp_365 = tmp_356*tmp_52;
      real_t tmp_366 = tmp_363 + tmp_364 + tmp_365;
      real_t tmp_367 = tmp_1*(tmp_358 - 1.0/4.0) + tmp_10*(tmp_366 - 1.0/4.0) + tmp_12*(tmp_362 - 1.0/4.0);
      real_t tmp_368 = -tmp_353 - tmp_355 - tmp_357 - tmp_359 - tmp_360 - tmp_361 - tmp_363 - tmp_364 - tmp_365 + 1;
      real_t tmp_369 = tmp_367*tmp_63;
      real_t tmp_370 = 0.042507265838595799*tmp_65;
      real_t tmp_371 = 0.031405749086161582*tmp_3 + 0.031405749086161582*tmp_4 + tmp_6;
      real_t tmp_372 = tmp_24*tmp_371;
      real_t tmp_373 = 0.031405749086161582*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29;
      real_t tmp_374 = tmp_31*tmp_373;
      real_t tmp_375 = 0.031405749086161582*tmp_34 + 0.031405749086161582*tmp_35 + tmp_36;
      real_t tmp_376 = tmp_375*tmp_38;
      real_t tmp_377 = tmp_372 + tmp_374 + tmp_376;
      real_t tmp_378 = tmp_371*tmp_41;
      real_t tmp_379 = tmp_373*tmp_43;
      real_t tmp_380 = tmp_375*tmp_45;
      real_t tmp_381 = tmp_378 + tmp_379 + tmp_380;
      real_t tmp_382 = tmp_371*tmp_48;
      real_t tmp_383 = tmp_373*tmp_50;
      real_t tmp_384 = tmp_375*tmp_52;
      real_t tmp_385 = tmp_382 + tmp_383 + tmp_384;
      real_t tmp_386 = tmp_1*(tmp_377 - 1.0/4.0) + tmp_10*(tmp_385 - 1.0/4.0) + tmp_12*(tmp_381 - 1.0/4.0);
      real_t tmp_387 = -tmp_372 - tmp_374 - tmp_376 - tmp_378 - tmp_379 - tmp_380 - tmp_382 - tmp_383 - tmp_384 + 1;
      real_t tmp_388 = tmp_386*tmp_63;
      real_t tmp_389 = 0.0068572537431980923*tmp_65;
      real_t tmp_390 = 0.19601935860219369*tmp_3 + 0.19601935860219369*tmp_4 + tmp_6;
      real_t tmp_391 = tmp_24*tmp_390;
      real_t tmp_392 = 0.19601935860219369*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29;
      real_t tmp_393 = tmp_31*tmp_392;
      real_t tmp_394 = 0.19601935860219369*tmp_34 + 0.19601935860219369*tmp_35 + tmp_36;
      real_t tmp_395 = tmp_38*tmp_394;
      real_t tmp_396 = tmp_391 + tmp_393 + tmp_395;
      real_t tmp_397 = tmp_390*tmp_41;
      real_t tmp_398 = tmp_392*tmp_43;
      real_t tmp_399 = tmp_394*tmp_45;
      real_t tmp_400 = tmp_397 + tmp_398 + tmp_399;
      real_t tmp_401 = tmp_390*tmp_48;
      real_t tmp_402 = tmp_392*tmp_50;
      real_t tmp_403 = tmp_394*tmp_52;
      real_t tmp_404 = tmp_401 + tmp_402 + tmp_403;
      real_t tmp_405 = tmp_1*(tmp_396 - 1.0/4.0) + tmp_10*(tmp_404 - 1.0/4.0) + tmp_12*(tmp_400 - 1.0/4.0);
      real_t tmp_406 = -tmp_391 - tmp_393 - tmp_395 - tmp_397 - tmp_398 - tmp_399 - tmp_401 - tmp_402 - tmp_403 + 1;
      real_t tmp_407 = tmp_405*tmp_63;
      real_t tmp_408 = 0.037198804536718075*tmp_65;
      real_t tmp_409 = 0.40446199974765351*tmp_3 + 0.40446199974765351*tmp_4 + tmp_6;
      real_t tmp_410 = tmp_24*tmp_409;
      real_t tmp_411 = 0.40446199974765351*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29;
      real_t tmp_412 = tmp_31*tmp_411;
      real_t tmp_413 = 0.40446199974765351*tmp_34 + 0.40446199974765351*tmp_35 + tmp_36;
      real_t tmp_414 = tmp_38*tmp_413;
      real_t tmp_415 = tmp_410 + tmp_412 + tmp_414;
      real_t tmp_416 = tmp_409*tmp_41;
      real_t tmp_417 = tmp_411*tmp_43;
      real_t tmp_418 = tmp_413*tmp_45;
      real_t tmp_419 = tmp_416 + tmp_417 + tmp_418;
      real_t tmp_420 = tmp_409*tmp_48;
      real_t tmp_421 = tmp_411*tmp_50;
      real_t tmp_422 = tmp_413*tmp_52;
      real_t tmp_423 = tmp_420 + tmp_421 + tmp_422;
      real_t tmp_424 = tmp_1*(tmp_415 - 1.0/4.0) + tmp_10*(tmp_423 - 1.0/4.0) + tmp_12*(tmp_419 - 1.0/4.0);
      real_t tmp_425 = -tmp_410 - tmp_412 - tmp_414 - tmp_416 - tmp_417 - tmp_418 - tmp_420 - tmp_421 - tmp_422 + 1;
      real_t tmp_426 = tmp_424*tmp_63;
      real_t tmp_427 = 0.042507265838595799*tmp_65;
      real_t tmp_428 = 0.1711304259088916*tmp_3 + 0.041227165399737475*tmp_4 + tmp_6;
      real_t tmp_429 = tmp_24*tmp_428;
      real_t tmp_430 = 0.1711304259088916*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29;
      real_t tmp_431 = tmp_31*tmp_430;
      real_t tmp_432 = 0.1711304259088916*tmp_34 + 0.041227165399737475*tmp_35 + tmp_36;
      real_t tmp_433 = tmp_38*tmp_432;
      real_t tmp_434 = tmp_429 + tmp_431 + tmp_433;
      real_t tmp_435 = tmp_41*tmp_428;
      real_t tmp_436 = tmp_43*tmp_430;
      real_t tmp_437 = tmp_432*tmp_45;
      real_t tmp_438 = tmp_435 + tmp_436 + tmp_437;
      real_t tmp_439 = tmp_428*tmp_48;
      real_t tmp_440 = tmp_430*tmp_50;
      real_t tmp_441 = tmp_432*tmp_52;
      real_t tmp_442 = tmp_439 + tmp_440 + tmp_441;
      real_t tmp_443 = tmp_1*(tmp_434 - 1.0/4.0) + tmp_10*(tmp_442 - 1.0/4.0) + tmp_12*(tmp_438 - 1.0/4.0);
      real_t tmp_444 = -tmp_429 - tmp_431 - tmp_433 - tmp_435 - tmp_436 - tmp_437 - tmp_439 - tmp_440 - tmp_441 + 1;
      real_t tmp_445 = tmp_443*tmp_63;
      real_t tmp_446 = 0.019202922745021479*tmp_65;
      real_t tmp_447 = 0.5*p_affine_13_0*tmp_38 + 0.5*p_affine_13_1*tmp_31 + 0.5*p_affine_13_2*tmp_24;
      real_t tmp_448 = 0.5*p_affine_13_0*tmp_45 + 0.5*p_affine_13_1*tmp_43 + 0.5*p_affine_13_2*tmp_41;
      real_t tmp_449 = 0.5*p_affine_13_0*tmp_52 + 0.5*p_affine_13_1*tmp_50 + 0.5*p_affine_13_2*tmp_48;
      real_t a_0_0 = tmp_104*(-tmp_101*tmp_56 + tmp_102*tmp_103 - tmp_102*tmp_58) + tmp_123*(-tmp_120*tmp_56 + tmp_121*tmp_122 - tmp_121*tmp_58) + tmp_142*(-tmp_139*tmp_56 + tmp_140*tmp_141 - tmp_140*tmp_58) + tmp_161*(-tmp_158*tmp_56 + tmp_159*tmp_160 - tmp_159*tmp_58) + tmp_180*(-tmp_177*tmp_56 + tmp_178*tmp_179 - tmp_178*tmp_58) + tmp_199*(-tmp_196*tmp_56 + tmp_197*tmp_198 - tmp_197*tmp_58) + tmp_218*(-tmp_215*tmp_56 + tmp_216*tmp_217 - tmp_216*tmp_58) + tmp_237*(-tmp_234*tmp_56 + tmp_235*tmp_236 - tmp_235*tmp_58) + tmp_256*(-tmp_253*tmp_56 + tmp_254*tmp_255 - tmp_254*tmp_58) + tmp_275*(-tmp_272*tmp_56 + tmp_273*tmp_274 - tmp_273*tmp_58) + tmp_294*(-tmp_291*tmp_56 + tmp_292*tmp_293 - tmp_292*tmp_58) + tmp_313*(-tmp_310*tmp_56 + tmp_311*tmp_312 - tmp_311*tmp_58) + tmp_332*(-tmp_329*tmp_56 + tmp_330*tmp_331 - tmp_330*tmp_58) + tmp_351*(-tmp_348*tmp_56 + tmp_349*tmp_350 - tmp_349*tmp_58) + tmp_370*(-tmp_367*tmp_56 + tmp_368*tmp_369 - tmp_368*tmp_58) + tmp_389*(-tmp_386*tmp_56 + tmp_387*tmp_388 - tmp_387*tmp_58) + tmp_408*(-tmp_405*tmp_56 + tmp_406*tmp_407 - tmp_406*tmp_58) + tmp_427*(-tmp_424*tmp_56 + tmp_425*tmp_426 - tmp_425*tmp_58) + tmp_446*(-tmp_443*tmp_56 + tmp_444*tmp_445 - tmp_444*tmp_58) + tmp_66*(-tmp_55*tmp_56 - tmp_57*tmp_58 + tmp_57*tmp_64) + tmp_85*(-tmp_56*tmp_82 - tmp_58*tmp_83 + tmp_83*tmp_84);
      real_t a_0_1 = tmp_104*(-tmp_101*tmp_447 + tmp_103*tmp_92 - tmp_58*tmp_92) + tmp_123*(tmp_111*tmp_122 - tmp_111*tmp_58 - tmp_120*tmp_447) + tmp_142*(tmp_130*tmp_141 - tmp_130*tmp_58 - tmp_139*tmp_447) + tmp_161*(tmp_149*tmp_160 - tmp_149*tmp_58 - tmp_158*tmp_447) + tmp_180*(tmp_168*tmp_179 - tmp_168*tmp_58 - tmp_177*tmp_447) + tmp_199*(tmp_187*tmp_198 - tmp_187*tmp_58 - tmp_196*tmp_447) + tmp_218*(tmp_206*tmp_217 - tmp_206*tmp_58 - tmp_215*tmp_447) + tmp_237*(tmp_225*tmp_236 - tmp_225*tmp_58 - tmp_234*tmp_447) + tmp_256*(tmp_244*tmp_255 - tmp_244*tmp_58 - tmp_253*tmp_447) + tmp_275*(tmp_263*tmp_274 - tmp_263*tmp_58 - tmp_272*tmp_447) + tmp_294*(tmp_282*tmp_293 - tmp_282*tmp_58 - tmp_291*tmp_447) + tmp_313*(tmp_301*tmp_312 - tmp_301*tmp_58 - tmp_310*tmp_447) + tmp_332*(tmp_320*tmp_331 - tmp_320*tmp_58 - tmp_329*tmp_447) + tmp_351*(tmp_339*tmp_350 - tmp_339*tmp_58 - tmp_348*tmp_447) + tmp_370*(tmp_358*tmp_369 - tmp_358*tmp_58 - tmp_367*tmp_447) + tmp_389*(tmp_377*tmp_388 - tmp_377*tmp_58 - tmp_386*tmp_447) + tmp_408*(tmp_396*tmp_407 - tmp_396*tmp_58 - tmp_405*tmp_447) + tmp_427*(tmp_415*tmp_426 - tmp_415*tmp_58 - tmp_424*tmp_447) + tmp_446*(tmp_434*tmp_445 - tmp_434*tmp_58 - tmp_443*tmp_447) + tmp_66*(-tmp_40*tmp_58 + tmp_40*tmp_64 - tmp_447*tmp_55) + tmp_85*(-tmp_447*tmp_82 - tmp_58*tmp_73 + tmp_73*tmp_84);
      real_t a_0_2 = tmp_104*(-tmp_101*tmp_448 + tmp_103*tmp_96 - tmp_58*tmp_96) + tmp_123*(tmp_115*tmp_122 - tmp_115*tmp_58 - tmp_120*tmp_448) + tmp_142*(tmp_134*tmp_141 - tmp_134*tmp_58 - tmp_139*tmp_448) + tmp_161*(tmp_153*tmp_160 - tmp_153*tmp_58 - tmp_158*tmp_448) + tmp_180*(tmp_172*tmp_179 - tmp_172*tmp_58 - tmp_177*tmp_448) + tmp_199*(tmp_191*tmp_198 - tmp_191*tmp_58 - tmp_196*tmp_448) + tmp_218*(tmp_210*tmp_217 - tmp_210*tmp_58 - tmp_215*tmp_448) + tmp_237*(tmp_229*tmp_236 - tmp_229*tmp_58 - tmp_234*tmp_448) + tmp_256*(tmp_248*tmp_255 - tmp_248*tmp_58 - tmp_253*tmp_448) + tmp_275*(tmp_267*tmp_274 - tmp_267*tmp_58 - tmp_272*tmp_448) + tmp_294*(tmp_286*tmp_293 - tmp_286*tmp_58 - tmp_291*tmp_448) + tmp_313*(tmp_305*tmp_312 - tmp_305*tmp_58 - tmp_310*tmp_448) + tmp_332*(tmp_324*tmp_331 - tmp_324*tmp_58 - tmp_329*tmp_448) + tmp_351*(tmp_343*tmp_350 - tmp_343*tmp_58 - tmp_348*tmp_448) + tmp_370*(tmp_362*tmp_369 - tmp_362*tmp_58 - tmp_367*tmp_448) + tmp_389*(tmp_381*tmp_388 - tmp_381*tmp_58 - tmp_386*tmp_448) + tmp_408*(tmp_400*tmp_407 - tmp_400*tmp_58 - tmp_405*tmp_448) + tmp_427*(tmp_419*tmp_426 - tmp_419*tmp_58 - tmp_424*tmp_448) + tmp_446*(tmp_438*tmp_445 - tmp_438*tmp_58 - tmp_443*tmp_448) + tmp_66*(-tmp_448*tmp_55 - tmp_47*tmp_58 + tmp_47*tmp_64) + tmp_85*(-tmp_448*tmp_82 - tmp_58*tmp_77 + tmp_77*tmp_84);
      real_t a_0_3 = tmp_104*(tmp_100*tmp_103 - tmp_100*tmp_58 - tmp_101*tmp_449) + tmp_123*(tmp_119*tmp_122 - tmp_119*tmp_58 - tmp_120*tmp_449) + tmp_142*(tmp_138*tmp_141 - tmp_138*tmp_58 - tmp_139*tmp_449) + tmp_161*(tmp_157*tmp_160 - tmp_157*tmp_58 - tmp_158*tmp_449) + tmp_180*(tmp_176*tmp_179 - tmp_176*tmp_58 - tmp_177*tmp_449) + tmp_199*(tmp_195*tmp_198 - tmp_195*tmp_58 - tmp_196*tmp_449) + tmp_218*(tmp_214*tmp_217 - tmp_214*tmp_58 - tmp_215*tmp_449) + tmp_237*(tmp_233*tmp_236 - tmp_233*tmp_58 - tmp_234*tmp_449) + tmp_256*(tmp_252*tmp_255 - tmp_252*tmp_58 - tmp_253*tmp_449) + tmp_275*(tmp_271*tmp_274 - tmp_271*tmp_58 - tmp_272*tmp_449) + tmp_294*(tmp_290*tmp_293 - tmp_290*tmp_58 - tmp_291*tmp_449) + tmp_313*(tmp_309*tmp_312 - tmp_309*tmp_58 - tmp_310*tmp_449) + tmp_332*(tmp_328*tmp_331 - tmp_328*tmp_58 - tmp_329*tmp_449) + tmp_351*(tmp_347*tmp_350 - tmp_347*tmp_58 - tmp_348*tmp_449) + tmp_370*(tmp_366*tmp_369 - tmp_366*tmp_58 - tmp_367*tmp_449) + tmp_389*(tmp_385*tmp_388 - tmp_385*tmp_58 - tmp_386*tmp_449) + tmp_408*(tmp_404*tmp_407 - tmp_404*tmp_58 - tmp_405*tmp_449) + tmp_427*(tmp_423*tmp_426 - tmp_423*tmp_58 - tmp_424*tmp_449) + tmp_446*(tmp_442*tmp_445 - tmp_442*tmp_58 - tmp_443*tmp_449) + tmp_66*(-tmp_449*tmp_55 - tmp_54*tmp_58 + tmp_54*tmp_64) + tmp_85*(-tmp_449*tmp_82 - tmp_58*tmp_81 + tmp_81*tmp_84);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
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


      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_1_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_0;
      real_t tmp_3 = p_affine_2_0 + tmp_2;
      real_t tmp_4 = p_affine_3_1 + tmp_0;
      real_t tmp_5 = p_affine_3_0 + tmp_2;
      real_t tmp_6 = p_affine_2_1 + tmp_0;
      real_t tmp_7 = tmp_3*tmp_4 - tmp_5*tmp_6;
      real_t tmp_8 = p_affine_1_0 + tmp_2;
      real_t tmp_9 = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = tmp_10*tmp_6;
      real_t tmp_12 = p_affine_1_2 + tmp_9;
      real_t tmp_13 = tmp_12*tmp_4;
      real_t tmp_14 = p_affine_2_2 + tmp_9;
      real_t tmp_15 = tmp_1*tmp_14;
      real_t tmp_16 = tmp_14*tmp_4;
      real_t tmp_17 = tmp_1*tmp_10;
      real_t tmp_18 = tmp_12*tmp_6;
      real_t tmp_19 = 1.0 / (tmp_11*tmp_8 + tmp_13*tmp_3 + tmp_15*tmp_5 - tmp_16*tmp_8 - tmp_17*tmp_3 - tmp_18*tmp_5);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = 0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22;
      real_t tmp_24 = p_affine_8_2 + tmp_9;
      real_t tmp_25 = tmp_19*(tmp_23 + tmp_24);
      real_t tmp_26 = -tmp_10*tmp_3 + tmp_14*tmp_5;
      real_t tmp_27 = -p_affine_8_1;
      real_t tmp_28 = p_affine_9_1 + tmp_27;
      real_t tmp_29 = p_affine_10_1 + tmp_27;
      real_t tmp_30 = 0.031405749086161582*tmp_28 + 0.93718850182767688*tmp_29;
      real_t tmp_31 = p_affine_8_1 + tmp_0;
      real_t tmp_32 = tmp_19*(tmp_30 + tmp_31);
      real_t tmp_33 = tmp_11 - tmp_16;
      real_t tmp_34 = -p_affine_8_0;
      real_t tmp_35 = p_affine_9_0 + tmp_34;
      real_t tmp_36 = p_affine_10_0 + tmp_34;
      real_t tmp_37 = 0.031405749086161582*tmp_35 + 0.93718850182767688*tmp_36;
      real_t tmp_38 = p_affine_8_0 + tmp_2;
      real_t tmp_39 = tmp_19*(tmp_37 + tmp_38);
      real_t tmp_40 = tmp_1*tmp_5 - tmp_4*tmp_8;
      real_t tmp_41 = tmp_10*tmp_8 - tmp_12*tmp_5;
      real_t tmp_42 = tmp_13 - tmp_17;
      real_t tmp_43 = -tmp_1*tmp_3 + tmp_6*tmp_8;
      real_t tmp_44 = tmp_12*tmp_3 - tmp_14*tmp_8;
      real_t tmp_45 = tmp_15 - tmp_18;
      real_t tmp_46 = tmp_1*(tmp_25*tmp_7 + tmp_26*tmp_32 + tmp_33*tmp_39 - 1.0/4.0) + tmp_4*(tmp_25*tmp_43 + tmp_32*tmp_44 + tmp_39*tmp_45 - 1.0/4.0) + tmp_6*(tmp_25*tmp_40 + tmp_32*tmp_41 + tmp_39*tmp_42 - 1.0/4.0);
      real_t tmp_47 = -p_affine_4_1;
      real_t tmp_48 = p_affine_5_1 + tmp_47;
      real_t tmp_49 = -p_affine_4_2;
      real_t tmp_50 = p_affine_6_2 + tmp_49;
      real_t tmp_51 = tmp_48*tmp_50;
      real_t tmp_52 = p_affine_6_1 + tmp_47;
      real_t tmp_53 = p_affine_5_2 + tmp_49;
      real_t tmp_54 = p_affine_7_2 + tmp_49;
      real_t tmp_55 = -p_affine_4_0;
      real_t tmp_56 = p_affine_5_0 + tmp_55;
      real_t tmp_57 = tmp_52*tmp_56;
      real_t tmp_58 = p_affine_7_1 + tmp_47;
      real_t tmp_59 = p_affine_6_0 + tmp_55;
      real_t tmp_60 = tmp_53*tmp_59;
      real_t tmp_61 = p_affine_7_0 + tmp_55;
      real_t tmp_62 = tmp_56*tmp_58;
      real_t tmp_63 = tmp_48*tmp_59;
      real_t tmp_64 = tmp_53*tmp_61;
      real_t tmp_65 = 1.0 / (-tmp_50*tmp_62 + tmp_51*tmp_61 - tmp_52*tmp_64 + tmp_54*tmp_57 - tmp_54*tmp_63 + tmp_58*tmp_60);
      real_t tmp_66 = tmp_65*(tmp_51 - tmp_52*tmp_53);
      real_t tmp_67 = tmp_65*(-tmp_48*tmp_54 + tmp_53*tmp_58);
      real_t tmp_68 = tmp_65*(-tmp_50*tmp_58 + tmp_52*tmp_54);
      real_t tmp_69 = tmp_65*(-tmp_50*tmp_56 + tmp_60);
      real_t tmp_70 = tmp_65*(tmp_54*tmp_56 - tmp_64);
      real_t tmp_71 = tmp_65*(tmp_50*tmp_61 - tmp_54*tmp_59);
      real_t tmp_72 = tmp_65*(tmp_57 - tmp_63);
      real_t tmp_73 = tmp_65*(tmp_48*tmp_61 - tmp_62);
      real_t tmp_74 = tmp_65*(-tmp_52*tmp_61 + tmp_58*tmp_59);
      real_t tmp_75 = 0.5*p_affine_13_0*(-tmp_66 - tmp_67 - tmp_68) + 0.5*p_affine_13_1*(-tmp_69 - tmp_70 - tmp_71) + 0.5*p_affine_13_2*(-tmp_72 - tmp_73 - tmp_74);
      real_t tmp_76 = p_affine_8_2 + tmp_49;
      real_t tmp_77 = tmp_23 + tmp_76;
      real_t tmp_78 = tmp_72*tmp_77;
      real_t tmp_79 = tmp_73*tmp_77;
      real_t tmp_80 = p_affine_8_1 + tmp_47;
      real_t tmp_81 = tmp_30 + tmp_80;
      real_t tmp_82 = tmp_69*tmp_81;
      real_t tmp_83 = tmp_70*tmp_81;
      real_t tmp_84 = tmp_74*tmp_77;
      real_t tmp_85 = tmp_71*tmp_81;
      real_t tmp_86 = p_affine_8_0 + tmp_55;
      real_t tmp_87 = tmp_37 + tmp_86;
      real_t tmp_88 = tmp_66*tmp_87;
      real_t tmp_89 = tmp_67*tmp_87;
      real_t tmp_90 = tmp_68*tmp_87;
      real_t tmp_91 = -tmp_78 - tmp_79 - tmp_82 - tmp_83 - tmp_84 - tmp_85 - tmp_88 - tmp_89 - tmp_90 + 1;
      real_t tmp_92 = tmp_1*tmp_19;
      real_t tmp_93 = tmp_19*tmp_6;
      real_t tmp_94 = tmp_19*tmp_4;
      real_t tmp_95 = 0.5*p_affine_13_0*(tmp_33*tmp_92 + tmp_42*tmp_93 + tmp_45*tmp_94) + 0.5*p_affine_13_1*(tmp_26*tmp_92 + tmp_41*tmp_93 + tmp_44*tmp_94) + 0.5*p_affine_13_2*(tmp_40*tmp_93 + tmp_43*tmp_94 + tmp_7*tmp_92);
      real_t tmp_96 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_97 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_98 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_99 = (std::abs(tmp_22*tmp_96 - tmp_29*tmp_98)*std::abs(tmp_22*tmp_96 - tmp_29*tmp_98)) + (std::abs(tmp_22*tmp_97 - tmp_36*tmp_98)*std::abs(tmp_22*tmp_97 - tmp_36*tmp_98)) + (std::abs(tmp_29*tmp_97 - tmp_36*tmp_96)*std::abs(tmp_29*tmp_97 - tmp_36*tmp_96));
      real_t tmp_100 = 5.0*std::pow(tmp_99, -0.25);
      real_t tmp_101 = tmp_100*tmp_46;
      real_t tmp_102 = 1.0*std::pow(tmp_99, 1.0/2.0);
      real_t tmp_103 = 0.0068572537431980923*tmp_102;
      real_t tmp_104 = 0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22;
      real_t tmp_105 = tmp_19*(tmp_104 + tmp_24);
      real_t tmp_106 = 0.19601935860219369*tmp_28 + 0.60796128279561268*tmp_29;
      real_t tmp_107 = tmp_19*(tmp_106 + tmp_31);
      real_t tmp_108 = 0.19601935860219369*tmp_35 + 0.60796128279561268*tmp_36;
      real_t tmp_109 = tmp_19*(tmp_108 + tmp_38);
      real_t tmp_110 = tmp_1*(tmp_105*tmp_7 + tmp_107*tmp_26 + tmp_109*tmp_33 - 1.0/4.0) + tmp_4*(tmp_105*tmp_43 + tmp_107*tmp_44 + tmp_109*tmp_45 - 1.0/4.0) + tmp_6*(tmp_105*tmp_40 + tmp_107*tmp_41 + tmp_109*tmp_42 - 1.0/4.0);
      real_t tmp_111 = tmp_104 + tmp_76;
      real_t tmp_112 = tmp_111*tmp_72;
      real_t tmp_113 = tmp_111*tmp_73;
      real_t tmp_114 = tmp_106 + tmp_80;
      real_t tmp_115 = tmp_114*tmp_69;
      real_t tmp_116 = tmp_114*tmp_70;
      real_t tmp_117 = tmp_111*tmp_74;
      real_t tmp_118 = tmp_114*tmp_71;
      real_t tmp_119 = tmp_108 + tmp_86;
      real_t tmp_120 = tmp_119*tmp_66;
      real_t tmp_121 = tmp_119*tmp_67;
      real_t tmp_122 = tmp_119*tmp_68;
      real_t tmp_123 = -tmp_112 - tmp_113 - tmp_115 - tmp_116 - tmp_117 - tmp_118 - tmp_120 - tmp_121 - tmp_122 + 1;
      real_t tmp_124 = tmp_100*tmp_110;
      real_t tmp_125 = 0.037198804536718075*tmp_102;
      real_t tmp_126 = 0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_127 = tmp_19*(tmp_126 + tmp_24);
      real_t tmp_128 = 0.37605877282253791*tmp_28 + 0.039308471900058539*tmp_29;
      real_t tmp_129 = tmp_19*(tmp_128 + tmp_31);
      real_t tmp_130 = 0.37605877282253791*tmp_35 + 0.039308471900058539*tmp_36;
      real_t tmp_131 = tmp_19*(tmp_130 + tmp_38);
      real_t tmp_132 = tmp_1*(tmp_127*tmp_7 + tmp_129*tmp_26 + tmp_131*tmp_33 - 1.0/4.0) + tmp_4*(tmp_127*tmp_43 + tmp_129*tmp_44 + tmp_131*tmp_45 - 1.0/4.0) + tmp_6*(tmp_127*tmp_40 + tmp_129*tmp_41 + tmp_131*tmp_42 - 1.0/4.0);
      real_t tmp_133 = tmp_126 + tmp_76;
      real_t tmp_134 = tmp_133*tmp_72;
      real_t tmp_135 = tmp_133*tmp_73;
      real_t tmp_136 = tmp_128 + tmp_80;
      real_t tmp_137 = tmp_136*tmp_69;
      real_t tmp_138 = tmp_136*tmp_70;
      real_t tmp_139 = tmp_133*tmp_74;
      real_t tmp_140 = tmp_136*tmp_71;
      real_t tmp_141 = tmp_130 + tmp_86;
      real_t tmp_142 = tmp_141*tmp_66;
      real_t tmp_143 = tmp_141*tmp_67;
      real_t tmp_144 = tmp_141*tmp_68;
      real_t tmp_145 = -tmp_134 - tmp_135 - tmp_137 - tmp_138 - tmp_139 - tmp_140 - tmp_142 - tmp_143 - tmp_144 + 1;
      real_t tmp_146 = tmp_100*tmp_132;
      real_t tmp_147 = 0.020848748529055869*tmp_102;
      real_t tmp_148 = 0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_149 = tmp_19*(tmp_148 + tmp_24);
      real_t tmp_150 = 0.78764240869137092*tmp_28 + 0.1711304259088916*tmp_29;
      real_t tmp_151 = tmp_19*(tmp_150 + tmp_31);
      real_t tmp_152 = 0.78764240869137092*tmp_35 + 0.1711304259088916*tmp_36;
      real_t tmp_153 = tmp_19*(tmp_152 + tmp_38);
      real_t tmp_154 = tmp_1*(tmp_149*tmp_7 + tmp_151*tmp_26 + tmp_153*tmp_33 - 1.0/4.0) + tmp_4*(tmp_149*tmp_43 + tmp_151*tmp_44 + tmp_153*tmp_45 - 1.0/4.0) + tmp_6*(tmp_149*tmp_40 + tmp_151*tmp_41 + tmp_153*tmp_42 - 1.0/4.0);
      real_t tmp_155 = tmp_148 + tmp_76;
      real_t tmp_156 = tmp_155*tmp_72;
      real_t tmp_157 = tmp_155*tmp_73;
      real_t tmp_158 = tmp_150 + tmp_80;
      real_t tmp_159 = tmp_158*tmp_69;
      real_t tmp_160 = tmp_158*tmp_70;
      real_t tmp_161 = tmp_155*tmp_74;
      real_t tmp_162 = tmp_158*tmp_71;
      real_t tmp_163 = tmp_152 + tmp_86;
      real_t tmp_164 = tmp_163*tmp_66;
      real_t tmp_165 = tmp_163*tmp_67;
      real_t tmp_166 = tmp_163*tmp_68;
      real_t tmp_167 = -tmp_156 - tmp_157 - tmp_159 - tmp_160 - tmp_161 - tmp_162 - tmp_164 - tmp_165 - tmp_166 + 1;
      real_t tmp_168 = tmp_100*tmp_154;
      real_t tmp_169 = 0.019202922745021479*tmp_102;
      real_t tmp_170 = 0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_171 = tmp_19*(tmp_170 + tmp_24);
      real_t tmp_172 = 0.58463275527740355*tmp_28 + 0.37605877282253791*tmp_29;
      real_t tmp_173 = tmp_19*(tmp_172 + tmp_31);
      real_t tmp_174 = 0.58463275527740355*tmp_35 + 0.37605877282253791*tmp_36;
      real_t tmp_175 = tmp_19*(tmp_174 + tmp_38);
      real_t tmp_176 = tmp_1*(tmp_171*tmp_7 + tmp_173*tmp_26 + tmp_175*tmp_33 - 1.0/4.0) + tmp_4*(tmp_171*tmp_43 + tmp_173*tmp_44 + tmp_175*tmp_45 - 1.0/4.0) + tmp_6*(tmp_171*tmp_40 + tmp_173*tmp_41 + tmp_175*tmp_42 - 1.0/4.0);
      real_t tmp_177 = tmp_170 + tmp_76;
      real_t tmp_178 = tmp_177*tmp_72;
      real_t tmp_179 = tmp_177*tmp_73;
      real_t tmp_180 = tmp_172 + tmp_80;
      real_t tmp_181 = tmp_180*tmp_69;
      real_t tmp_182 = tmp_180*tmp_70;
      real_t tmp_183 = tmp_177*tmp_74;
      real_t tmp_184 = tmp_180*tmp_71;
      real_t tmp_185 = tmp_174 + tmp_86;
      real_t tmp_186 = tmp_185*tmp_66;
      real_t tmp_187 = tmp_185*tmp_67;
      real_t tmp_188 = tmp_185*tmp_68;
      real_t tmp_189 = -tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 - tmp_184 - tmp_186 - tmp_187 - tmp_188 + 1;
      real_t tmp_190 = tmp_100*tmp_176;
      real_t tmp_191 = 0.020848748529055869*tmp_102;
      real_t tmp_192 = 0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_193 = tmp_19*(tmp_192 + tmp_24);
      real_t tmp_194 = 0.041227165399737475*tmp_28 + 0.78764240869137092*tmp_29;
      real_t tmp_195 = tmp_19*(tmp_194 + tmp_31);
      real_t tmp_196 = 0.041227165399737475*tmp_35 + 0.78764240869137092*tmp_36;
      real_t tmp_197 = tmp_19*(tmp_196 + tmp_38);
      real_t tmp_198 = tmp_1*(tmp_193*tmp_7 + tmp_195*tmp_26 + tmp_197*tmp_33 - 1.0/4.0) + tmp_4*(tmp_193*tmp_43 + tmp_195*tmp_44 + tmp_197*tmp_45 - 1.0/4.0) + tmp_6*(tmp_193*tmp_40 + tmp_195*tmp_41 + tmp_197*tmp_42 - 1.0/4.0);
      real_t tmp_199 = tmp_192 + tmp_76;
      real_t tmp_200 = tmp_199*tmp_72;
      real_t tmp_201 = tmp_199*tmp_73;
      real_t tmp_202 = tmp_194 + tmp_80;
      real_t tmp_203 = tmp_202*tmp_69;
      real_t tmp_204 = tmp_202*tmp_70;
      real_t tmp_205 = tmp_199*tmp_74;
      real_t tmp_206 = tmp_202*tmp_71;
      real_t tmp_207 = tmp_196 + tmp_86;
      real_t tmp_208 = tmp_207*tmp_66;
      real_t tmp_209 = tmp_207*tmp_67;
      real_t tmp_210 = tmp_207*tmp_68;
      real_t tmp_211 = -tmp_200 - tmp_201 - tmp_203 - tmp_204 - tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 + 1;
      real_t tmp_212 = tmp_100*tmp_198;
      real_t tmp_213 = 0.019202922745021479*tmp_102;
      real_t tmp_214 = 0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_215 = tmp_19*(tmp_214 + tmp_24);
      real_t tmp_216 = 0.039308471900058539*tmp_28 + 0.58463275527740355*tmp_29;
      real_t tmp_217 = tmp_19*(tmp_216 + tmp_31);
      real_t tmp_218 = 0.039308471900058539*tmp_35 + 0.58463275527740355*tmp_36;
      real_t tmp_219 = tmp_19*(tmp_218 + tmp_38);
      real_t tmp_220 = tmp_1*(tmp_215*tmp_7 + tmp_217*tmp_26 + tmp_219*tmp_33 - 1.0/4.0) + tmp_4*(tmp_215*tmp_43 + tmp_217*tmp_44 + tmp_219*tmp_45 - 1.0/4.0) + tmp_6*(tmp_215*tmp_40 + tmp_217*tmp_41 + tmp_219*tmp_42 - 1.0/4.0);
      real_t tmp_221 = tmp_214 + tmp_76;
      real_t tmp_222 = tmp_221*tmp_72;
      real_t tmp_223 = tmp_221*tmp_73;
      real_t tmp_224 = tmp_216 + tmp_80;
      real_t tmp_225 = tmp_224*tmp_69;
      real_t tmp_226 = tmp_224*tmp_70;
      real_t tmp_227 = tmp_221*tmp_74;
      real_t tmp_228 = tmp_224*tmp_71;
      real_t tmp_229 = tmp_218 + tmp_86;
      real_t tmp_230 = tmp_229*tmp_66;
      real_t tmp_231 = tmp_229*tmp_67;
      real_t tmp_232 = tmp_229*tmp_68;
      real_t tmp_233 = -tmp_222 - tmp_223 - tmp_225 - tmp_226 - tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 + 1;
      real_t tmp_234 = tmp_100*tmp_220;
      real_t tmp_235 = 0.020848748529055869*tmp_102;
      real_t tmp_236 = 0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_237 = tmp_19*(tmp_236 + tmp_24);
      real_t tmp_238 = 0.78764240869137092*tmp_28 + 0.041227165399737475*tmp_29;
      real_t tmp_239 = tmp_19*(tmp_238 + tmp_31);
      real_t tmp_240 = 0.78764240869137092*tmp_35 + 0.041227165399737475*tmp_36;
      real_t tmp_241 = tmp_19*(tmp_240 + tmp_38);
      real_t tmp_242 = tmp_1*(tmp_237*tmp_7 + tmp_239*tmp_26 + tmp_241*tmp_33 - 1.0/4.0) + tmp_4*(tmp_237*tmp_43 + tmp_239*tmp_44 + tmp_241*tmp_45 - 1.0/4.0) + tmp_6*(tmp_237*tmp_40 + tmp_239*tmp_41 + tmp_241*tmp_42 - 1.0/4.0);
      real_t tmp_243 = tmp_236 + tmp_76;
      real_t tmp_244 = tmp_243*tmp_72;
      real_t tmp_245 = tmp_243*tmp_73;
      real_t tmp_246 = tmp_238 + tmp_80;
      real_t tmp_247 = tmp_246*tmp_69;
      real_t tmp_248 = tmp_246*tmp_70;
      real_t tmp_249 = tmp_243*tmp_74;
      real_t tmp_250 = tmp_246*tmp_71;
      real_t tmp_251 = tmp_240 + tmp_86;
      real_t tmp_252 = tmp_251*tmp_66;
      real_t tmp_253 = tmp_251*tmp_67;
      real_t tmp_254 = tmp_251*tmp_68;
      real_t tmp_255 = -tmp_244 - tmp_245 - tmp_247 - tmp_248 - tmp_249 - tmp_250 - tmp_252 - tmp_253 - tmp_254 + 1;
      real_t tmp_256 = tmp_100*tmp_242;
      real_t tmp_257 = 0.019202922745021479*tmp_102;
      real_t tmp_258 = 0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_259 = tmp_19*(tmp_24 + tmp_258);
      real_t tmp_260 = 0.58463275527740355*tmp_28 + 0.039308471900058539*tmp_29;
      real_t tmp_261 = tmp_19*(tmp_260 + tmp_31);
      real_t tmp_262 = 0.58463275527740355*tmp_35 + 0.039308471900058539*tmp_36;
      real_t tmp_263 = tmp_19*(tmp_262 + tmp_38);
      real_t tmp_264 = tmp_1*(tmp_259*tmp_7 + tmp_26*tmp_261 + tmp_263*tmp_33 - 1.0/4.0) + tmp_4*(tmp_259*tmp_43 + tmp_261*tmp_44 + tmp_263*tmp_45 - 1.0/4.0) + tmp_6*(tmp_259*tmp_40 + tmp_261*tmp_41 + tmp_263*tmp_42 - 1.0/4.0);
      real_t tmp_265 = tmp_258 + tmp_76;
      real_t tmp_266 = tmp_265*tmp_72;
      real_t tmp_267 = tmp_265*tmp_73;
      real_t tmp_268 = tmp_260 + tmp_80;
      real_t tmp_269 = tmp_268*tmp_69;
      real_t tmp_270 = tmp_268*tmp_70;
      real_t tmp_271 = tmp_265*tmp_74;
      real_t tmp_272 = tmp_268*tmp_71;
      real_t tmp_273 = tmp_262 + tmp_86;
      real_t tmp_274 = tmp_273*tmp_66;
      real_t tmp_275 = tmp_273*tmp_67;
      real_t tmp_276 = tmp_273*tmp_68;
      real_t tmp_277 = -tmp_266 - tmp_267 - tmp_269 - tmp_270 - tmp_271 - tmp_272 - tmp_274 - tmp_275 - tmp_276 + 1;
      real_t tmp_278 = tmp_100*tmp_264;
      real_t tmp_279 = 0.020848748529055869*tmp_102;
      real_t tmp_280 = 0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_281 = tmp_19*(tmp_24 + tmp_280);
      real_t tmp_282 = 0.1711304259088916*tmp_28 + 0.78764240869137092*tmp_29;
      real_t tmp_283 = tmp_19*(tmp_282 + tmp_31);
      real_t tmp_284 = 0.1711304259088916*tmp_35 + 0.78764240869137092*tmp_36;
      real_t tmp_285 = tmp_19*(tmp_284 + tmp_38);
      real_t tmp_286 = tmp_1*(tmp_26*tmp_283 + tmp_281*tmp_7 + tmp_285*tmp_33 - 1.0/4.0) + tmp_4*(tmp_281*tmp_43 + tmp_283*tmp_44 + tmp_285*tmp_45 - 1.0/4.0) + tmp_6*(tmp_281*tmp_40 + tmp_283*tmp_41 + tmp_285*tmp_42 - 1.0/4.0);
      real_t tmp_287 = tmp_280 + tmp_76;
      real_t tmp_288 = tmp_287*tmp_72;
      real_t tmp_289 = tmp_287*tmp_73;
      real_t tmp_290 = tmp_282 + tmp_80;
      real_t tmp_291 = tmp_290*tmp_69;
      real_t tmp_292 = tmp_290*tmp_70;
      real_t tmp_293 = tmp_287*tmp_74;
      real_t tmp_294 = tmp_290*tmp_71;
      real_t tmp_295 = tmp_284 + tmp_86;
      real_t tmp_296 = tmp_295*tmp_66;
      real_t tmp_297 = tmp_295*tmp_67;
      real_t tmp_298 = tmp_295*tmp_68;
      real_t tmp_299 = -tmp_288 - tmp_289 - tmp_291 - tmp_292 - tmp_293 - tmp_294 - tmp_296 - tmp_297 - tmp_298 + 1;
      real_t tmp_300 = tmp_100*tmp_286;
      real_t tmp_301 = 0.019202922745021479*tmp_102;
      real_t tmp_302 = 0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_303 = tmp_19*(tmp_24 + tmp_302);
      real_t tmp_304 = 0.37605877282253791*tmp_28 + 0.58463275527740355*tmp_29;
      real_t tmp_305 = tmp_19*(tmp_304 + tmp_31);
      real_t tmp_306 = 0.37605877282253791*tmp_35 + 0.58463275527740355*tmp_36;
      real_t tmp_307 = tmp_19*(tmp_306 + tmp_38);
      real_t tmp_308 = tmp_1*(tmp_26*tmp_305 + tmp_303*tmp_7 + tmp_307*tmp_33 - 1.0/4.0) + tmp_4*(tmp_303*tmp_43 + tmp_305*tmp_44 + tmp_307*tmp_45 - 1.0/4.0) + tmp_6*(tmp_303*tmp_40 + tmp_305*tmp_41 + tmp_307*tmp_42 - 1.0/4.0);
      real_t tmp_309 = tmp_302 + tmp_76;
      real_t tmp_310 = tmp_309*tmp_72;
      real_t tmp_311 = tmp_309*tmp_73;
      real_t tmp_312 = tmp_304 + tmp_80;
      real_t tmp_313 = tmp_312*tmp_69;
      real_t tmp_314 = tmp_312*tmp_70;
      real_t tmp_315 = tmp_309*tmp_74;
      real_t tmp_316 = tmp_312*tmp_71;
      real_t tmp_317 = tmp_306 + tmp_86;
      real_t tmp_318 = tmp_317*tmp_66;
      real_t tmp_319 = tmp_317*tmp_67;
      real_t tmp_320 = tmp_317*tmp_68;
      real_t tmp_321 = -tmp_310 - tmp_311 - tmp_313 - tmp_314 - tmp_315 - tmp_316 - tmp_318 - tmp_319 - tmp_320 + 1;
      real_t tmp_322 = tmp_100*tmp_308;
      real_t tmp_323 = 0.020848748529055869*tmp_102;
      real_t tmp_324 = 0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_325 = tmp_19*(tmp_24 + tmp_324);
      real_t tmp_326 = 0.041227165399737475*tmp_28 + 0.1711304259088916*tmp_29;
      real_t tmp_327 = tmp_19*(tmp_31 + tmp_326);
      real_t tmp_328 = 0.041227165399737475*tmp_35 + 0.1711304259088916*tmp_36;
      real_t tmp_329 = tmp_19*(tmp_328 + tmp_38);
      real_t tmp_330 = tmp_1*(tmp_26*tmp_327 + tmp_325*tmp_7 + tmp_329*tmp_33 - 1.0/4.0) + tmp_4*(tmp_325*tmp_43 + tmp_327*tmp_44 + tmp_329*tmp_45 - 1.0/4.0) + tmp_6*(tmp_325*tmp_40 + tmp_327*tmp_41 + tmp_329*tmp_42 - 1.0/4.0);
      real_t tmp_331 = tmp_324 + tmp_76;
      real_t tmp_332 = tmp_331*tmp_72;
      real_t tmp_333 = tmp_331*tmp_73;
      real_t tmp_334 = tmp_326 + tmp_80;
      real_t tmp_335 = tmp_334*tmp_69;
      real_t tmp_336 = tmp_334*tmp_70;
      real_t tmp_337 = tmp_331*tmp_74;
      real_t tmp_338 = tmp_334*tmp_71;
      real_t tmp_339 = tmp_328 + tmp_86;
      real_t tmp_340 = tmp_339*tmp_66;
      real_t tmp_341 = tmp_339*tmp_67;
      real_t tmp_342 = tmp_339*tmp_68;
      real_t tmp_343 = -tmp_332 - tmp_333 - tmp_335 - tmp_336 - tmp_337 - tmp_338 - tmp_340 - tmp_341 - tmp_342 + 1;
      real_t tmp_344 = tmp_100*tmp_330;
      real_t tmp_345 = 0.019202922745021479*tmp_102;
      real_t tmp_346 = 0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22;
      real_t tmp_347 = tmp_19*(tmp_24 + tmp_346);
      real_t tmp_348 = 0.40446199974765351*tmp_28 + 0.19107600050469298*tmp_29;
      real_t tmp_349 = tmp_19*(tmp_31 + tmp_348);
      real_t tmp_350 = 0.40446199974765351*tmp_35 + 0.19107600050469298*tmp_36;
      real_t tmp_351 = tmp_19*(tmp_350 + tmp_38);
      real_t tmp_352 = tmp_1*(tmp_26*tmp_349 + tmp_33*tmp_351 + tmp_347*tmp_7 - 1.0/4.0) + tmp_4*(tmp_347*tmp_43 + tmp_349*tmp_44 + tmp_351*tmp_45 - 1.0/4.0) + tmp_6*(tmp_347*tmp_40 + tmp_349*tmp_41 + tmp_351*tmp_42 - 1.0/4.0);
      real_t tmp_353 = tmp_346 + tmp_76;
      real_t tmp_354 = tmp_353*tmp_72;
      real_t tmp_355 = tmp_353*tmp_73;
      real_t tmp_356 = tmp_348 + tmp_80;
      real_t tmp_357 = tmp_356*tmp_69;
      real_t tmp_358 = tmp_356*tmp_70;
      real_t tmp_359 = tmp_353*tmp_74;
      real_t tmp_360 = tmp_356*tmp_71;
      real_t tmp_361 = tmp_350 + tmp_86;
      real_t tmp_362 = tmp_361*tmp_66;
      real_t tmp_363 = tmp_361*tmp_67;
      real_t tmp_364 = tmp_361*tmp_68;
      real_t tmp_365 = -tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 - tmp_360 - tmp_362 - tmp_363 - tmp_364 + 1;
      real_t tmp_366 = tmp_100*tmp_352;
      real_t tmp_367 = 0.042507265838595799*tmp_102;
      real_t tmp_368 = 0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_369 = tmp_19*(tmp_24 + tmp_368);
      real_t tmp_370 = 0.039308471900058539*tmp_28 + 0.37605877282253791*tmp_29;
      real_t tmp_371 = tmp_19*(tmp_31 + tmp_370);
      real_t tmp_372 = 0.039308471900058539*tmp_35 + 0.37605877282253791*tmp_36;
      real_t tmp_373 = tmp_19*(tmp_372 + tmp_38);
      real_t tmp_374 = tmp_1*(tmp_26*tmp_371 + tmp_33*tmp_373 + tmp_369*tmp_7 - 1.0/4.0) + tmp_4*(tmp_369*tmp_43 + tmp_371*tmp_44 + tmp_373*tmp_45 - 1.0/4.0) + tmp_6*(tmp_369*tmp_40 + tmp_371*tmp_41 + tmp_373*tmp_42 - 1.0/4.0);
      real_t tmp_375 = tmp_368 + tmp_76;
      real_t tmp_376 = tmp_375*tmp_72;
      real_t tmp_377 = tmp_375*tmp_73;
      real_t tmp_378 = tmp_370 + tmp_80;
      real_t tmp_379 = tmp_378*tmp_69;
      real_t tmp_380 = tmp_378*tmp_70;
      real_t tmp_381 = tmp_375*tmp_74;
      real_t tmp_382 = tmp_378*tmp_71;
      real_t tmp_383 = tmp_372 + tmp_86;
      real_t tmp_384 = tmp_383*tmp_66;
      real_t tmp_385 = tmp_383*tmp_67;
      real_t tmp_386 = tmp_383*tmp_68;
      real_t tmp_387 = -tmp_376 - tmp_377 - tmp_379 - tmp_380 - tmp_381 - tmp_382 - tmp_384 - tmp_385 - tmp_386 + 1;
      real_t tmp_388 = tmp_100*tmp_374;
      real_t tmp_389 = 0.020848748529055869*tmp_102;
      real_t tmp_390 = 0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_391 = tmp_19*(tmp_24 + tmp_390);
      real_t tmp_392 = 0.93718850182767688*tmp_28 + 0.031405749086161582*tmp_29;
      real_t tmp_393 = tmp_19*(tmp_31 + tmp_392);
      real_t tmp_394 = 0.93718850182767688*tmp_35 + 0.031405749086161582*tmp_36;
      real_t tmp_395 = tmp_19*(tmp_38 + tmp_394);
      real_t tmp_396 = tmp_1*(tmp_26*tmp_393 + tmp_33*tmp_395 + tmp_391*tmp_7 - 1.0/4.0) + tmp_4*(tmp_391*tmp_43 + tmp_393*tmp_44 + tmp_395*tmp_45 - 1.0/4.0) + tmp_6*(tmp_391*tmp_40 + tmp_393*tmp_41 + tmp_395*tmp_42 - 1.0/4.0);
      real_t tmp_397 = tmp_390 + tmp_76;
      real_t tmp_398 = tmp_397*tmp_72;
      real_t tmp_399 = tmp_397*tmp_73;
      real_t tmp_400 = tmp_392 + tmp_80;
      real_t tmp_401 = tmp_400*tmp_69;
      real_t tmp_402 = tmp_400*tmp_70;
      real_t tmp_403 = tmp_397*tmp_74;
      real_t tmp_404 = tmp_400*tmp_71;
      real_t tmp_405 = tmp_394 + tmp_86;
      real_t tmp_406 = tmp_405*tmp_66;
      real_t tmp_407 = tmp_405*tmp_67;
      real_t tmp_408 = tmp_405*tmp_68;
      real_t tmp_409 = -tmp_398 - tmp_399 - tmp_401 - tmp_402 - tmp_403 - tmp_404 - tmp_406 - tmp_407 - tmp_408 + 1;
      real_t tmp_410 = tmp_100*tmp_396;
      real_t tmp_411 = 0.0068572537431980923*tmp_102;
      real_t tmp_412 = 0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_413 = tmp_19*(tmp_24 + tmp_412);
      real_t tmp_414 = 0.60796128279561268*tmp_28 + 0.19601935860219369*tmp_29;
      real_t tmp_415 = tmp_19*(tmp_31 + tmp_414);
      real_t tmp_416 = 0.60796128279561268*tmp_35 + 0.19601935860219369*tmp_36;
      real_t tmp_417 = tmp_19*(tmp_38 + tmp_416);
      real_t tmp_418 = tmp_1*(tmp_26*tmp_415 + tmp_33*tmp_417 + tmp_413*tmp_7 - 1.0/4.0) + tmp_4*(tmp_413*tmp_43 + tmp_415*tmp_44 + tmp_417*tmp_45 - 1.0/4.0) + tmp_6*(tmp_40*tmp_413 + tmp_41*tmp_415 + tmp_417*tmp_42 - 1.0/4.0);
      real_t tmp_419 = tmp_412 + tmp_76;
      real_t tmp_420 = tmp_419*tmp_72;
      real_t tmp_421 = tmp_419*tmp_73;
      real_t tmp_422 = tmp_414 + tmp_80;
      real_t tmp_423 = tmp_422*tmp_69;
      real_t tmp_424 = tmp_422*tmp_70;
      real_t tmp_425 = tmp_419*tmp_74;
      real_t tmp_426 = tmp_422*tmp_71;
      real_t tmp_427 = tmp_416 + tmp_86;
      real_t tmp_428 = tmp_427*tmp_66;
      real_t tmp_429 = tmp_427*tmp_67;
      real_t tmp_430 = tmp_427*tmp_68;
      real_t tmp_431 = -tmp_420 - tmp_421 - tmp_423 - tmp_424 - tmp_425 - tmp_426 - tmp_428 - tmp_429 - tmp_430 + 1;
      real_t tmp_432 = tmp_100*tmp_418;
      real_t tmp_433 = 0.037198804536718075*tmp_102;
      real_t tmp_434 = 0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_435 = tmp_19*(tmp_24 + tmp_434);
      real_t tmp_436 = 0.19107600050469298*tmp_28 + 0.40446199974765351*tmp_29;
      real_t tmp_437 = tmp_19*(tmp_31 + tmp_436);
      real_t tmp_438 = 0.19107600050469298*tmp_35 + 0.40446199974765351*tmp_36;
      real_t tmp_439 = tmp_19*(tmp_38 + tmp_438);
      real_t tmp_440 = tmp_1*(tmp_26*tmp_437 + tmp_33*tmp_439 + tmp_435*tmp_7 - 1.0/4.0) + tmp_4*(tmp_43*tmp_435 + tmp_437*tmp_44 + tmp_439*tmp_45 - 1.0/4.0) + tmp_6*(tmp_40*tmp_435 + tmp_41*tmp_437 + tmp_42*tmp_439 - 1.0/4.0);
      real_t tmp_441 = tmp_434 + tmp_76;
      real_t tmp_442 = tmp_441*tmp_72;
      real_t tmp_443 = tmp_441*tmp_73;
      real_t tmp_444 = tmp_436 + tmp_80;
      real_t tmp_445 = tmp_444*tmp_69;
      real_t tmp_446 = tmp_444*tmp_70;
      real_t tmp_447 = tmp_441*tmp_74;
      real_t tmp_448 = tmp_444*tmp_71;
      real_t tmp_449 = tmp_438 + tmp_86;
      real_t tmp_450 = tmp_449*tmp_66;
      real_t tmp_451 = tmp_449*tmp_67;
      real_t tmp_452 = tmp_449*tmp_68;
      real_t tmp_453 = -tmp_442 - tmp_443 - tmp_445 - tmp_446 - tmp_447 - tmp_448 - tmp_450 - tmp_451 - tmp_452 + 1;
      real_t tmp_454 = tmp_100*tmp_440;
      real_t tmp_455 = 0.042507265838595799*tmp_102;
      real_t tmp_456 = 0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_457 = tmp_19*(tmp_24 + tmp_456);
      real_t tmp_458 = 0.031405749086161582*tmp_28 + 0.031405749086161582*tmp_29;
      real_t tmp_459 = tmp_19*(tmp_31 + tmp_458);
      real_t tmp_460 = 0.031405749086161582*tmp_35 + 0.031405749086161582*tmp_36;
      real_t tmp_461 = tmp_19*(tmp_38 + tmp_460);
      real_t tmp_462 = tmp_1*(tmp_26*tmp_459 + tmp_33*tmp_461 + tmp_457*tmp_7 - 1.0/4.0) + tmp_4*(tmp_43*tmp_457 + tmp_44*tmp_459 + tmp_45*tmp_461 - 1.0/4.0) + tmp_6*(tmp_40*tmp_457 + tmp_41*tmp_459 + tmp_42*tmp_461 - 1.0/4.0);
      real_t tmp_463 = tmp_456 + tmp_76;
      real_t tmp_464 = tmp_463*tmp_72;
      real_t tmp_465 = tmp_463*tmp_73;
      real_t tmp_466 = tmp_458 + tmp_80;
      real_t tmp_467 = tmp_466*tmp_69;
      real_t tmp_468 = tmp_466*tmp_70;
      real_t tmp_469 = tmp_463*tmp_74;
      real_t tmp_470 = tmp_466*tmp_71;
      real_t tmp_471 = tmp_460 + tmp_86;
      real_t tmp_472 = tmp_471*tmp_66;
      real_t tmp_473 = tmp_471*tmp_67;
      real_t tmp_474 = tmp_471*tmp_68;
      real_t tmp_475 = -tmp_464 - tmp_465 - tmp_467 - tmp_468 - tmp_469 - tmp_470 - tmp_472 - tmp_473 - tmp_474 + 1;
      real_t tmp_476 = tmp_100*tmp_462;
      real_t tmp_477 = 0.0068572537431980923*tmp_102;
      real_t tmp_478 = 0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_479 = tmp_19*(tmp_24 + tmp_478);
      real_t tmp_480 = 0.19601935860219369*tmp_28 + 0.19601935860219369*tmp_29;
      real_t tmp_481 = tmp_19*(tmp_31 + tmp_480);
      real_t tmp_482 = 0.19601935860219369*tmp_35 + 0.19601935860219369*tmp_36;
      real_t tmp_483 = tmp_19*(tmp_38 + tmp_482);
      real_t tmp_484 = tmp_1*(tmp_26*tmp_481 + tmp_33*tmp_483 + tmp_479*tmp_7 - 1.0/4.0) + tmp_4*(tmp_43*tmp_479 + tmp_44*tmp_481 + tmp_45*tmp_483 - 1.0/4.0) + tmp_6*(tmp_40*tmp_479 + tmp_41*tmp_481 + tmp_42*tmp_483 - 1.0/4.0);
      real_t tmp_485 = tmp_478 + tmp_76;
      real_t tmp_486 = tmp_485*tmp_72;
      real_t tmp_487 = tmp_485*tmp_73;
      real_t tmp_488 = tmp_480 + tmp_80;
      real_t tmp_489 = tmp_488*tmp_69;
      real_t tmp_490 = tmp_488*tmp_70;
      real_t tmp_491 = tmp_485*tmp_74;
      real_t tmp_492 = tmp_488*tmp_71;
      real_t tmp_493 = tmp_482 + tmp_86;
      real_t tmp_494 = tmp_493*tmp_66;
      real_t tmp_495 = tmp_493*tmp_67;
      real_t tmp_496 = tmp_493*tmp_68;
      real_t tmp_497 = -tmp_486 - tmp_487 - tmp_489 - tmp_490 - tmp_491 - tmp_492 - tmp_494 - tmp_495 - tmp_496 + 1;
      real_t tmp_498 = tmp_100*tmp_484;
      real_t tmp_499 = 0.037198804536718075*tmp_102;
      real_t tmp_500 = 0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_501 = tmp_19*(tmp_24 + tmp_500);
      real_t tmp_502 = 0.40446199974765351*tmp_28 + 0.40446199974765351*tmp_29;
      real_t tmp_503 = tmp_19*(tmp_31 + tmp_502);
      real_t tmp_504 = 0.40446199974765351*tmp_35 + 0.40446199974765351*tmp_36;
      real_t tmp_505 = tmp_19*(tmp_38 + tmp_504);
      real_t tmp_506 = tmp_1*(tmp_26*tmp_503 + tmp_33*tmp_505 + tmp_501*tmp_7 - 1.0/4.0) + tmp_4*(tmp_43*tmp_501 + tmp_44*tmp_503 + tmp_45*tmp_505 - 1.0/4.0) + tmp_6*(tmp_40*tmp_501 + tmp_41*tmp_503 + tmp_42*tmp_505 - 1.0/4.0);
      real_t tmp_507 = tmp_500 + tmp_76;
      real_t tmp_508 = tmp_507*tmp_72;
      real_t tmp_509 = tmp_507*tmp_73;
      real_t tmp_510 = tmp_502 + tmp_80;
      real_t tmp_511 = tmp_510*tmp_69;
      real_t tmp_512 = tmp_510*tmp_70;
      real_t tmp_513 = tmp_507*tmp_74;
      real_t tmp_514 = tmp_510*tmp_71;
      real_t tmp_515 = tmp_504 + tmp_86;
      real_t tmp_516 = tmp_515*tmp_66;
      real_t tmp_517 = tmp_515*tmp_67;
      real_t tmp_518 = tmp_515*tmp_68;
      real_t tmp_519 = -tmp_508 - tmp_509 - tmp_511 - tmp_512 - tmp_513 - tmp_514 - tmp_516 - tmp_517 - tmp_518 + 1;
      real_t tmp_520 = tmp_100*tmp_506;
      real_t tmp_521 = 0.042507265838595799*tmp_102;
      real_t tmp_522 = 0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_523 = tmp_19*(tmp_24 + tmp_522);
      real_t tmp_524 = 0.1711304259088916*tmp_28 + 0.041227165399737475*tmp_29;
      real_t tmp_525 = tmp_19*(tmp_31 + tmp_524);
      real_t tmp_526 = 0.1711304259088916*tmp_35 + 0.041227165399737475*tmp_36;
      real_t tmp_527 = tmp_19*(tmp_38 + tmp_526);
      real_t tmp_528 = tmp_1*(tmp_26*tmp_525 + tmp_33*tmp_527 + tmp_523*tmp_7 - 1.0/4.0) + tmp_4*(tmp_43*tmp_523 + tmp_44*tmp_525 + tmp_45*tmp_527 - 1.0/4.0) + tmp_6*(tmp_40*tmp_523 + tmp_41*tmp_525 + tmp_42*tmp_527 - 1.0/4.0);
      real_t tmp_529 = tmp_522 + tmp_76;
      real_t tmp_530 = tmp_529*tmp_72;
      real_t tmp_531 = tmp_529*tmp_73;
      real_t tmp_532 = tmp_524 + tmp_80;
      real_t tmp_533 = tmp_532*tmp_69;
      real_t tmp_534 = tmp_532*tmp_70;
      real_t tmp_535 = tmp_529*tmp_74;
      real_t tmp_536 = tmp_532*tmp_71;
      real_t tmp_537 = tmp_526 + tmp_86;
      real_t tmp_538 = tmp_537*tmp_66;
      real_t tmp_539 = tmp_537*tmp_67;
      real_t tmp_540 = tmp_537*tmp_68;
      real_t tmp_541 = -tmp_530 - tmp_531 - tmp_533 - tmp_534 - tmp_535 - tmp_536 - tmp_538 - tmp_539 - tmp_540 + 1;
      real_t tmp_542 = tmp_100*tmp_528;
      real_t tmp_543 = 0.019202922745021479*tmp_102;
      real_t tmp_544 = tmp_84 + tmp_85 + tmp_90;
      real_t tmp_545 = 0.5*p_affine_13_0*tmp_68 + 0.5*p_affine_13_1*tmp_71 + 0.5*p_affine_13_2*tmp_74;
      real_t tmp_546 = tmp_117 + tmp_118 + tmp_122;
      real_t tmp_547 = tmp_139 + tmp_140 + tmp_144;
      real_t tmp_548 = tmp_161 + tmp_162 + tmp_166;
      real_t tmp_549 = tmp_183 + tmp_184 + tmp_188;
      real_t tmp_550 = tmp_205 + tmp_206 + tmp_210;
      real_t tmp_551 = tmp_227 + tmp_228 + tmp_232;
      real_t tmp_552 = tmp_249 + tmp_250 + tmp_254;
      real_t tmp_553 = tmp_271 + tmp_272 + tmp_276;
      real_t tmp_554 = tmp_293 + tmp_294 + tmp_298;
      real_t tmp_555 = tmp_315 + tmp_316 + tmp_320;
      real_t tmp_556 = tmp_337 + tmp_338 + tmp_342;
      real_t tmp_557 = tmp_359 + tmp_360 + tmp_364;
      real_t tmp_558 = tmp_381 + tmp_382 + tmp_386;
      real_t tmp_559 = tmp_403 + tmp_404 + tmp_408;
      real_t tmp_560 = tmp_425 + tmp_426 + tmp_430;
      real_t tmp_561 = tmp_447 + tmp_448 + tmp_452;
      real_t tmp_562 = tmp_469 + tmp_470 + tmp_474;
      real_t tmp_563 = tmp_491 + tmp_492 + tmp_496;
      real_t tmp_564 = tmp_513 + tmp_514 + tmp_518;
      real_t tmp_565 = tmp_535 + tmp_536 + tmp_540;
      real_t tmp_566 = tmp_79 + tmp_83 + tmp_89;
      real_t tmp_567 = 0.5*p_affine_13_0*tmp_67 + 0.5*p_affine_13_1*tmp_70 + 0.5*p_affine_13_2*tmp_73;
      real_t tmp_568 = tmp_113 + tmp_116 + tmp_121;
      real_t tmp_569 = tmp_135 + tmp_138 + tmp_143;
      real_t tmp_570 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_571 = tmp_179 + tmp_182 + tmp_187;
      real_t tmp_572 = tmp_201 + tmp_204 + tmp_209;
      real_t tmp_573 = tmp_223 + tmp_226 + tmp_231;
      real_t tmp_574 = tmp_245 + tmp_248 + tmp_253;
      real_t tmp_575 = tmp_267 + tmp_270 + tmp_275;
      real_t tmp_576 = tmp_289 + tmp_292 + tmp_297;
      real_t tmp_577 = tmp_311 + tmp_314 + tmp_319;
      real_t tmp_578 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_579 = tmp_355 + tmp_358 + tmp_363;
      real_t tmp_580 = tmp_377 + tmp_380 + tmp_385;
      real_t tmp_581 = tmp_399 + tmp_402 + tmp_407;
      real_t tmp_582 = tmp_421 + tmp_424 + tmp_429;
      real_t tmp_583 = tmp_443 + tmp_446 + tmp_451;
      real_t tmp_584 = tmp_465 + tmp_468 + tmp_473;
      real_t tmp_585 = tmp_487 + tmp_490 + tmp_495;
      real_t tmp_586 = tmp_509 + tmp_512 + tmp_517;
      real_t tmp_587 = tmp_531 + tmp_534 + tmp_539;
      real_t tmp_588 = tmp_78 + tmp_82 + tmp_88;
      real_t tmp_589 = 0.5*p_affine_13_0*tmp_66 + 0.5*p_affine_13_1*tmp_69 + 0.5*p_affine_13_2*tmp_72;
      real_t tmp_590 = tmp_112 + tmp_115 + tmp_120;
      real_t tmp_591 = tmp_134 + tmp_137 + tmp_142;
      real_t tmp_592 = tmp_156 + tmp_159 + tmp_164;
      real_t tmp_593 = tmp_178 + tmp_181 + tmp_186;
      real_t tmp_594 = tmp_200 + tmp_203 + tmp_208;
      real_t tmp_595 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_596 = tmp_244 + tmp_247 + tmp_252;
      real_t tmp_597 = tmp_266 + tmp_269 + tmp_274;
      real_t tmp_598 = tmp_288 + tmp_291 + tmp_296;
      real_t tmp_599 = tmp_310 + tmp_313 + tmp_318;
      real_t tmp_600 = tmp_332 + tmp_335 + tmp_340;
      real_t tmp_601 = tmp_354 + tmp_357 + tmp_362;
      real_t tmp_602 = tmp_376 + tmp_379 + tmp_384;
      real_t tmp_603 = tmp_398 + tmp_401 + tmp_406;
      real_t tmp_604 = tmp_420 + tmp_423 + tmp_428;
      real_t tmp_605 = tmp_442 + tmp_445 + tmp_450;
      real_t tmp_606 = tmp_464 + tmp_467 + tmp_472;
      real_t tmp_607 = tmp_486 + tmp_489 + tmp_494;
      real_t tmp_608 = tmp_508 + tmp_511 + tmp_516;
      real_t tmp_609 = tmp_530 + tmp_533 + tmp_538;
      real_t a_0_0 = tmp_103*(-tmp_101*tmp_91 - tmp_46*tmp_75 + tmp_91*tmp_95) + tmp_125*(-tmp_110*tmp_75 - tmp_123*tmp_124 + tmp_123*tmp_95) + tmp_147*(-tmp_132*tmp_75 - tmp_145*tmp_146 + tmp_145*tmp_95) + tmp_169*(-tmp_154*tmp_75 - tmp_167*tmp_168 + tmp_167*tmp_95) + tmp_191*(-tmp_176*tmp_75 - tmp_189*tmp_190 + tmp_189*tmp_95) + tmp_213*(-tmp_198*tmp_75 - tmp_211*tmp_212 + tmp_211*tmp_95) + tmp_235*(-tmp_220*tmp_75 - tmp_233*tmp_234 + tmp_233*tmp_95) + tmp_257*(-tmp_242*tmp_75 - tmp_255*tmp_256 + tmp_255*tmp_95) + tmp_279*(-tmp_264*tmp_75 - tmp_277*tmp_278 + tmp_277*tmp_95) + tmp_301*(-tmp_286*tmp_75 - tmp_299*tmp_300 + tmp_299*tmp_95) + tmp_323*(-tmp_308*tmp_75 - tmp_321*tmp_322 + tmp_321*tmp_95) + tmp_345*(-tmp_330*tmp_75 - tmp_343*tmp_344 + tmp_343*tmp_95) + tmp_367*(-tmp_352*tmp_75 - tmp_365*tmp_366 + tmp_365*tmp_95) + tmp_389*(-tmp_374*tmp_75 - tmp_387*tmp_388 + tmp_387*tmp_95) + tmp_411*(-tmp_396*tmp_75 - tmp_409*tmp_410 + tmp_409*tmp_95) + tmp_433*(-tmp_418*tmp_75 - tmp_431*tmp_432 + tmp_431*tmp_95) + tmp_455*(-tmp_440*tmp_75 - tmp_453*tmp_454 + tmp_453*tmp_95) + tmp_477*(-tmp_462*tmp_75 - tmp_475*tmp_476 + tmp_475*tmp_95) + tmp_499*(-tmp_484*tmp_75 - tmp_497*tmp_498 + tmp_497*tmp_95) + tmp_521*(-tmp_506*tmp_75 - tmp_519*tmp_520 + tmp_519*tmp_95) + tmp_543*(-tmp_528*tmp_75 - tmp_541*tmp_542 + tmp_541*tmp_95);
      real_t a_0_1 = tmp_103*(-tmp_101*tmp_544 - tmp_46*tmp_545 + tmp_544*tmp_95) + tmp_125*(-tmp_110*tmp_545 - tmp_124*tmp_546 + tmp_546*tmp_95) + tmp_147*(-tmp_132*tmp_545 - tmp_146*tmp_547 + tmp_547*tmp_95) + tmp_169*(-tmp_154*tmp_545 - tmp_168*tmp_548 + tmp_548*tmp_95) + tmp_191*(-tmp_176*tmp_545 - tmp_190*tmp_549 + tmp_549*tmp_95) + tmp_213*(-tmp_198*tmp_545 - tmp_212*tmp_550 + tmp_550*tmp_95) + tmp_235*(-tmp_220*tmp_545 - tmp_234*tmp_551 + tmp_551*tmp_95) + tmp_257*(-tmp_242*tmp_545 - tmp_256*tmp_552 + tmp_552*tmp_95) + tmp_279*(-tmp_264*tmp_545 - tmp_278*tmp_553 + tmp_553*tmp_95) + tmp_301*(-tmp_286*tmp_545 - tmp_300*tmp_554 + tmp_554*tmp_95) + tmp_323*(-tmp_308*tmp_545 - tmp_322*tmp_555 + tmp_555*tmp_95) + tmp_345*(-tmp_330*tmp_545 - tmp_344*tmp_556 + tmp_556*tmp_95) + tmp_367*(-tmp_352*tmp_545 - tmp_366*tmp_557 + tmp_557*tmp_95) + tmp_389*(-tmp_374*tmp_545 - tmp_388*tmp_558 + tmp_558*tmp_95) + tmp_411*(-tmp_396*tmp_545 - tmp_410*tmp_559 + tmp_559*tmp_95) + tmp_433*(-tmp_418*tmp_545 - tmp_432*tmp_560 + tmp_560*tmp_95) + tmp_455*(-tmp_440*tmp_545 - tmp_454*tmp_561 + tmp_561*tmp_95) + tmp_477*(-tmp_462*tmp_545 - tmp_476*tmp_562 + tmp_562*tmp_95) + tmp_499*(-tmp_484*tmp_545 - tmp_498*tmp_563 + tmp_563*tmp_95) + tmp_521*(-tmp_506*tmp_545 - tmp_520*tmp_564 + tmp_564*tmp_95) + tmp_543*(-tmp_528*tmp_545 - tmp_542*tmp_565 + tmp_565*tmp_95);
      real_t a_0_2 = tmp_103*(-tmp_101*tmp_566 - tmp_46*tmp_567 + tmp_566*tmp_95) + tmp_125*(-tmp_110*tmp_567 - tmp_124*tmp_568 + tmp_568*tmp_95) + tmp_147*(-tmp_132*tmp_567 - tmp_146*tmp_569 + tmp_569*tmp_95) + tmp_169*(-tmp_154*tmp_567 - tmp_168*tmp_570 + tmp_570*tmp_95) + tmp_191*(-tmp_176*tmp_567 - tmp_190*tmp_571 + tmp_571*tmp_95) + tmp_213*(-tmp_198*tmp_567 - tmp_212*tmp_572 + tmp_572*tmp_95) + tmp_235*(-tmp_220*tmp_567 - tmp_234*tmp_573 + tmp_573*tmp_95) + tmp_257*(-tmp_242*tmp_567 - tmp_256*tmp_574 + tmp_574*tmp_95) + tmp_279*(-tmp_264*tmp_567 - tmp_278*tmp_575 + tmp_575*tmp_95) + tmp_301*(-tmp_286*tmp_567 - tmp_300*tmp_576 + tmp_576*tmp_95) + tmp_323*(-tmp_308*tmp_567 - tmp_322*tmp_577 + tmp_577*tmp_95) + tmp_345*(-tmp_330*tmp_567 - tmp_344*tmp_578 + tmp_578*tmp_95) + tmp_367*(-tmp_352*tmp_567 - tmp_366*tmp_579 + tmp_579*tmp_95) + tmp_389*(-tmp_374*tmp_567 - tmp_388*tmp_580 + tmp_580*tmp_95) + tmp_411*(-tmp_396*tmp_567 - tmp_410*tmp_581 + tmp_581*tmp_95) + tmp_433*(-tmp_418*tmp_567 - tmp_432*tmp_582 + tmp_582*tmp_95) + tmp_455*(-tmp_440*tmp_567 - tmp_454*tmp_583 + tmp_583*tmp_95) + tmp_477*(-tmp_462*tmp_567 - tmp_476*tmp_584 + tmp_584*tmp_95) + tmp_499*(-tmp_484*tmp_567 - tmp_498*tmp_585 + tmp_585*tmp_95) + tmp_521*(-tmp_506*tmp_567 - tmp_520*tmp_586 + tmp_586*tmp_95) + tmp_543*(-tmp_528*tmp_567 - tmp_542*tmp_587 + tmp_587*tmp_95);
      real_t a_0_3 = tmp_103*(-tmp_101*tmp_588 - tmp_46*tmp_589 + tmp_588*tmp_95) + tmp_125*(-tmp_110*tmp_589 - tmp_124*tmp_590 + tmp_590*tmp_95) + tmp_147*(-tmp_132*tmp_589 - tmp_146*tmp_591 + tmp_591*tmp_95) + tmp_169*(-tmp_154*tmp_589 - tmp_168*tmp_592 + tmp_592*tmp_95) + tmp_191*(-tmp_176*tmp_589 - tmp_190*tmp_593 + tmp_593*tmp_95) + tmp_213*(-tmp_198*tmp_589 - tmp_212*tmp_594 + tmp_594*tmp_95) + tmp_235*(-tmp_220*tmp_589 - tmp_234*tmp_595 + tmp_595*tmp_95) + tmp_257*(-tmp_242*tmp_589 - tmp_256*tmp_596 + tmp_596*tmp_95) + tmp_279*(-tmp_264*tmp_589 - tmp_278*tmp_597 + tmp_597*tmp_95) + tmp_301*(-tmp_286*tmp_589 - tmp_300*tmp_598 + tmp_598*tmp_95) + tmp_323*(-tmp_308*tmp_589 - tmp_322*tmp_599 + tmp_599*tmp_95) + tmp_345*(-tmp_330*tmp_589 - tmp_344*tmp_600 + tmp_600*tmp_95) + tmp_367*(-tmp_352*tmp_589 - tmp_366*tmp_601 + tmp_601*tmp_95) + tmp_389*(-tmp_374*tmp_589 - tmp_388*tmp_602 + tmp_602*tmp_95) + tmp_411*(-tmp_396*tmp_589 - tmp_410*tmp_603 + tmp_603*tmp_95) + tmp_433*(-tmp_418*tmp_589 - tmp_432*tmp_604 + tmp_604*tmp_95) + tmp_455*(-tmp_440*tmp_589 - tmp_454*tmp_605 + tmp_605*tmp_95) + tmp_477*(-tmp_462*tmp_589 - tmp_476*tmp_606 + tmp_606*tmp_95) + tmp_499*(-tmp_484*tmp_589 - tmp_498*tmp_607 + tmp_607*tmp_95) + tmp_521*(-tmp_506*tmp_589 - tmp_520*tmp_608 + tmp_608*tmp_95) + tmp_543*(-tmp_528*tmp_589 - tmp_542*tmp_609 + tmp_609*tmp_95);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
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
      real_t tmp_22 = tmp_18*(tmp_12*tmp_6 - tmp_3*tmp_9);
      real_t tmp_23 = tmp_18*(tmp_10*tmp_9 - tmp_15*tmp_6);
      real_t tmp_24 = tmp_18*(-tmp_10*tmp_12 + tmp_15*tmp_3);
      real_t tmp_25 = tmp_18*(-tmp_1*tmp_12 + tmp_5*tmp_9);
      real_t tmp_26 = tmp_18*(tmp_1*tmp_15 - tmp_13*tmp_9);
      real_t tmp_27 = tmp_18*(tmp_12*tmp_13 - tmp_15*tmp_5);
      real_t tmp_28 = -p_affine_8_0;
      real_t tmp_29 = p_affine_10_0 + tmp_28;
      real_t tmp_30 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_31 = -p_affine_8_1;
      real_t tmp_32 = p_affine_10_1 + tmp_31;
      real_t tmp_33 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_34 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_35 = -p_affine_8_2;
      real_t tmp_36 = p_affine_10_2 + tmp_35;
      real_t tmp_37 = 1.0*std::pow((std::abs(tmp_29*tmp_30 - tmp_32*tmp_33)*std::abs(tmp_29*tmp_30 - tmp_32*tmp_33)) + (std::abs(tmp_29*tmp_34 - tmp_33*tmp_36)*std::abs(tmp_29*tmp_34 - tmp_33*tmp_36)) + (std::abs(tmp_30*tmp_36 - tmp_32*tmp_34)*std::abs(tmp_30*tmp_36 - tmp_32*tmp_34)), 1.0/2.0);
      real_t tmp_38 = tmp_37*(p_affine_13_0*(-tmp_19 - tmp_20 - tmp_21) + p_affine_13_1*(-tmp_22 - tmp_23 - tmp_24) + p_affine_13_2*(-tmp_25 - tmp_26 - tmp_27));
      real_t tmp_39 = p_affine_9_2 + tmp_35;
      real_t tmp_40 = p_affine_8_2 + tmp_2;
      real_t tmp_41 = 0.93718850182767688*tmp_36 + 0.031405749086161582*tmp_39 + tmp_40;
      real_t tmp_42 = p_affine_9_1 + tmp_31;
      real_t tmp_43 = p_affine_8_1 + tmp_0;
      real_t tmp_44 = 0.93718850182767688*tmp_32 + 0.031405749086161582*tmp_42 + tmp_43;
      real_t tmp_45 = p_affine_9_0 + tmp_28;
      real_t tmp_46 = p_affine_8_0 + tmp_8;
      real_t tmp_47 = 0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_45 + tmp_46;
      real_t tmp_48 = 0.0068572537431980923*tmp_1*(tmp_21*tmp_47 + tmp_24*tmp_44 + tmp_27*tmp_41 - 1.0/4.0) + 0.0068572537431980923*tmp_13*(tmp_19*tmp_47 + tmp_22*tmp_44 + tmp_25*tmp_41 - 1.0/4.0) + 0.0068572537431980923*tmp_5*(tmp_20*tmp_47 + tmp_23*tmp_44 + tmp_26*tmp_41 - 1.0/4.0);
      real_t tmp_49 = 0.60796128279561268*tmp_36 + 0.19601935860219369*tmp_39 + tmp_40;
      real_t tmp_50 = 0.60796128279561268*tmp_32 + 0.19601935860219369*tmp_42 + tmp_43;
      real_t tmp_51 = 0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_45 + tmp_46;
      real_t tmp_52 = 0.037198804536718075*tmp_1*(tmp_21*tmp_51 + tmp_24*tmp_50 + tmp_27*tmp_49 - 1.0/4.0) + 0.037198804536718075*tmp_13*(tmp_19*tmp_51 + tmp_22*tmp_50 + tmp_25*tmp_49 - 1.0/4.0) + 0.037198804536718075*tmp_5*(tmp_20*tmp_51 + tmp_23*tmp_50 + tmp_26*tmp_49 - 1.0/4.0);
      real_t tmp_53 = 0.039308471900058539*tmp_36 + 0.37605877282253791*tmp_39 + tmp_40;
      real_t tmp_54 = 0.039308471900058539*tmp_32 + 0.37605877282253791*tmp_42 + tmp_43;
      real_t tmp_55 = 0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_45 + tmp_46;
      real_t tmp_56 = 0.020848748529055869*tmp_1*(tmp_21*tmp_55 + tmp_24*tmp_54 + tmp_27*tmp_53 - 1.0/4.0) + 0.020848748529055869*tmp_13*(tmp_19*tmp_55 + tmp_22*tmp_54 + tmp_25*tmp_53 - 1.0/4.0) + 0.020848748529055869*tmp_5*(tmp_20*tmp_55 + tmp_23*tmp_54 + tmp_26*tmp_53 - 1.0/4.0);
      real_t tmp_57 = 0.1711304259088916*tmp_36 + 0.78764240869137092*tmp_39 + tmp_40;
      real_t tmp_58 = 0.1711304259088916*tmp_32 + 0.78764240869137092*tmp_42 + tmp_43;
      real_t tmp_59 = 0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_45 + tmp_46;
      real_t tmp_60 = 0.019202922745021479*tmp_1*(tmp_21*tmp_59 + tmp_24*tmp_58 + tmp_27*tmp_57 - 1.0/4.0) + 0.019202922745021479*tmp_13*(tmp_19*tmp_59 + tmp_22*tmp_58 + tmp_25*tmp_57 - 1.0/4.0) + 0.019202922745021479*tmp_5*(tmp_20*tmp_59 + tmp_23*tmp_58 + tmp_26*tmp_57 - 1.0/4.0);
      real_t tmp_61 = 0.37605877282253791*tmp_36 + 0.58463275527740355*tmp_39 + tmp_40;
      real_t tmp_62 = 0.37605877282253791*tmp_32 + 0.58463275527740355*tmp_42 + tmp_43;
      real_t tmp_63 = 0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_45 + tmp_46;
      real_t tmp_64 = 0.020848748529055869*tmp_1*(tmp_21*tmp_63 + tmp_24*tmp_62 + tmp_27*tmp_61 - 1.0/4.0) + 0.020848748529055869*tmp_13*(tmp_19*tmp_63 + tmp_22*tmp_62 + tmp_25*tmp_61 - 1.0/4.0) + 0.020848748529055869*tmp_5*(tmp_20*tmp_63 + tmp_23*tmp_62 + tmp_26*tmp_61 - 1.0/4.0);
      real_t tmp_65 = 0.78764240869137092*tmp_36 + 0.041227165399737475*tmp_39 + tmp_40;
      real_t tmp_66 = 0.78764240869137092*tmp_32 + 0.041227165399737475*tmp_42 + tmp_43;
      real_t tmp_67 = 0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_45 + tmp_46;
      real_t tmp_68 = 0.019202922745021479*tmp_1*(tmp_21*tmp_67 + tmp_24*tmp_66 + tmp_27*tmp_65 - 1.0/4.0) + 0.019202922745021479*tmp_13*(tmp_19*tmp_67 + tmp_22*tmp_66 + tmp_25*tmp_65 - 1.0/4.0) + 0.019202922745021479*tmp_5*(tmp_20*tmp_67 + tmp_23*tmp_66 + tmp_26*tmp_65 - 1.0/4.0);
      real_t tmp_69 = 0.58463275527740355*tmp_36 + 0.039308471900058539*tmp_39 + tmp_40;
      real_t tmp_70 = 0.58463275527740355*tmp_32 + 0.039308471900058539*tmp_42 + tmp_43;
      real_t tmp_71 = 0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_45 + tmp_46;
      real_t tmp_72 = 0.020848748529055869*tmp_1*(tmp_21*tmp_71 + tmp_24*tmp_70 + tmp_27*tmp_69 - 1.0/4.0) + 0.020848748529055869*tmp_13*(tmp_19*tmp_71 + tmp_22*tmp_70 + tmp_25*tmp_69 - 1.0/4.0) + 0.020848748529055869*tmp_5*(tmp_20*tmp_71 + tmp_23*tmp_70 + tmp_26*tmp_69 - 1.0/4.0);
      real_t tmp_73 = 0.041227165399737475*tmp_36 + 0.78764240869137092*tmp_39 + tmp_40;
      real_t tmp_74 = 0.041227165399737475*tmp_32 + 0.78764240869137092*tmp_42 + tmp_43;
      real_t tmp_75 = 0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_45 + tmp_46;
      real_t tmp_76 = 0.019202922745021479*tmp_1*(tmp_21*tmp_75 + tmp_24*tmp_74 + tmp_27*tmp_73 - 1.0/4.0) + 0.019202922745021479*tmp_13*(tmp_19*tmp_75 + tmp_22*tmp_74 + tmp_25*tmp_73 - 1.0/4.0) + 0.019202922745021479*tmp_5*(tmp_20*tmp_75 + tmp_23*tmp_74 + tmp_26*tmp_73 - 1.0/4.0);
      real_t tmp_77 = 0.039308471900058539*tmp_36 + 0.58463275527740355*tmp_39 + tmp_40;
      real_t tmp_78 = 0.039308471900058539*tmp_32 + 0.58463275527740355*tmp_42 + tmp_43;
      real_t tmp_79 = 0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_45 + tmp_46;
      real_t tmp_80 = 0.020848748529055869*tmp_1*(tmp_21*tmp_79 + tmp_24*tmp_78 + tmp_27*tmp_77 - 1.0/4.0) + 0.020848748529055869*tmp_13*(tmp_19*tmp_79 + tmp_22*tmp_78 + tmp_25*tmp_77 - 1.0/4.0) + 0.020848748529055869*tmp_5*(tmp_20*tmp_79 + tmp_23*tmp_78 + tmp_26*tmp_77 - 1.0/4.0);
      real_t tmp_81 = 0.78764240869137092*tmp_36 + 0.1711304259088916*tmp_39 + tmp_40;
      real_t tmp_82 = 0.78764240869137092*tmp_32 + 0.1711304259088916*tmp_42 + tmp_43;
      real_t tmp_83 = 0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_45 + tmp_46;
      real_t tmp_84 = 0.019202922745021479*tmp_1*(tmp_21*tmp_83 + tmp_24*tmp_82 + tmp_27*tmp_81 - 1.0/4.0) + 0.019202922745021479*tmp_13*(tmp_19*tmp_83 + tmp_22*tmp_82 + tmp_25*tmp_81 - 1.0/4.0) + 0.019202922745021479*tmp_5*(tmp_20*tmp_83 + tmp_23*tmp_82 + tmp_26*tmp_81 - 1.0/4.0);
      real_t tmp_85 = 0.58463275527740355*tmp_36 + 0.37605877282253791*tmp_39 + tmp_40;
      real_t tmp_86 = 0.58463275527740355*tmp_32 + 0.37605877282253791*tmp_42 + tmp_43;
      real_t tmp_87 = 0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_45 + tmp_46;
      real_t tmp_88 = 0.020848748529055869*tmp_1*(tmp_21*tmp_87 + tmp_24*tmp_86 + tmp_27*tmp_85 - 1.0/4.0) + 0.020848748529055869*tmp_13*(tmp_19*tmp_87 + tmp_22*tmp_86 + tmp_25*tmp_85 - 1.0/4.0) + 0.020848748529055869*tmp_5*(tmp_20*tmp_87 + tmp_23*tmp_86 + tmp_26*tmp_85 - 1.0/4.0);
      real_t tmp_89 = 0.1711304259088916*tmp_36 + 0.041227165399737475*tmp_39 + tmp_40;
      real_t tmp_90 = 0.1711304259088916*tmp_32 + 0.041227165399737475*tmp_42 + tmp_43;
      real_t tmp_91 = 0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_45 + tmp_46;
      real_t tmp_92 = 0.019202922745021479*tmp_1*(tmp_21*tmp_91 + tmp_24*tmp_90 + tmp_27*tmp_89 - 1.0/4.0) + 0.019202922745021479*tmp_13*(tmp_19*tmp_91 + tmp_22*tmp_90 + tmp_25*tmp_89 - 1.0/4.0) + 0.019202922745021479*tmp_5*(tmp_20*tmp_91 + tmp_23*tmp_90 + tmp_26*tmp_89 - 1.0/4.0);
      real_t tmp_93 = 0.19107600050469298*tmp_36 + 0.40446199974765351*tmp_39 + tmp_40;
      real_t tmp_94 = 0.19107600050469298*tmp_32 + 0.40446199974765351*tmp_42 + tmp_43;
      real_t tmp_95 = 0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_45 + tmp_46;
      real_t tmp_96 = 0.042507265838595799*tmp_1*(tmp_21*tmp_95 + tmp_24*tmp_94 + tmp_27*tmp_93 - 1.0/4.0) + 0.042507265838595799*tmp_13*(tmp_19*tmp_95 + tmp_22*tmp_94 + tmp_25*tmp_93 - 1.0/4.0) + 0.042507265838595799*tmp_5*(tmp_20*tmp_95 + tmp_23*tmp_94 + tmp_26*tmp_93 - 1.0/4.0);
      real_t tmp_97 = 0.37605877282253791*tmp_36 + 0.039308471900058539*tmp_39 + tmp_40;
      real_t tmp_98 = 0.37605877282253791*tmp_32 + 0.039308471900058539*tmp_42 + tmp_43;
      real_t tmp_99 = 0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_45 + tmp_46;
      real_t tmp_100 = 0.020848748529055869*tmp_1*(tmp_21*tmp_99 + tmp_24*tmp_98 + tmp_27*tmp_97 - 1.0/4.0) + 0.020848748529055869*tmp_13*(tmp_19*tmp_99 + tmp_22*tmp_98 + tmp_25*tmp_97 - 1.0/4.0) + 0.020848748529055869*tmp_5*(tmp_20*tmp_99 + tmp_23*tmp_98 + tmp_26*tmp_97 - 1.0/4.0);
      real_t tmp_101 = 0.031405749086161582*tmp_36 + 0.93718850182767688*tmp_39 + tmp_40;
      real_t tmp_102 = 0.031405749086161582*tmp_32 + 0.93718850182767688*tmp_42 + tmp_43;
      real_t tmp_103 = 0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_45 + tmp_46;
      real_t tmp_104 = 0.0068572537431980923*tmp_1*(tmp_101*tmp_27 + tmp_102*tmp_24 + tmp_103*tmp_21 - 1.0/4.0) + 0.0068572537431980923*tmp_13*(tmp_101*tmp_25 + tmp_102*tmp_22 + tmp_103*tmp_19 - 1.0/4.0) + 0.0068572537431980923*tmp_5*(tmp_101*tmp_26 + tmp_102*tmp_23 + tmp_103*tmp_20 - 1.0/4.0);
      real_t tmp_105 = 0.19601935860219369*tmp_36 + 0.60796128279561268*tmp_39 + tmp_40;
      real_t tmp_106 = 0.19601935860219369*tmp_32 + 0.60796128279561268*tmp_42 + tmp_43;
      real_t tmp_107 = 0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_45 + tmp_46;
      real_t tmp_108 = 0.037198804536718075*tmp_1*(tmp_105*tmp_27 + tmp_106*tmp_24 + tmp_107*tmp_21 - 1.0/4.0) + 0.037198804536718075*tmp_13*(tmp_105*tmp_25 + tmp_106*tmp_22 + tmp_107*tmp_19 - 1.0/4.0) + 0.037198804536718075*tmp_5*(tmp_105*tmp_26 + tmp_106*tmp_23 + tmp_107*tmp_20 - 1.0/4.0);
      real_t tmp_109 = 0.40446199974765351*tmp_36 + 0.19107600050469298*tmp_39 + tmp_40;
      real_t tmp_110 = 0.40446199974765351*tmp_32 + 0.19107600050469298*tmp_42 + tmp_43;
      real_t tmp_111 = 0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_45 + tmp_46;
      real_t tmp_112 = 0.042507265838595799*tmp_1*(tmp_109*tmp_27 + tmp_110*tmp_24 + tmp_111*tmp_21 - 1.0/4.0) + 0.042507265838595799*tmp_13*(tmp_109*tmp_25 + tmp_110*tmp_22 + tmp_111*tmp_19 - 1.0/4.0) + 0.042507265838595799*tmp_5*(tmp_109*tmp_26 + tmp_110*tmp_23 + tmp_111*tmp_20 - 1.0/4.0);
      real_t tmp_113 = 0.031405749086161582*tmp_36 + 0.031405749086161582*tmp_39 + tmp_40;
      real_t tmp_114 = 0.031405749086161582*tmp_32 + 0.031405749086161582*tmp_42 + tmp_43;
      real_t tmp_115 = 0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_45 + tmp_46;
      real_t tmp_116 = 0.0068572537431980923*tmp_1*(tmp_113*tmp_27 + tmp_114*tmp_24 + tmp_115*tmp_21 - 1.0/4.0) + 0.0068572537431980923*tmp_13*(tmp_113*tmp_25 + tmp_114*tmp_22 + tmp_115*tmp_19 - 1.0/4.0) + 0.0068572537431980923*tmp_5*(tmp_113*tmp_26 + tmp_114*tmp_23 + tmp_115*tmp_20 - 1.0/4.0);
      real_t tmp_117 = 0.19601935860219369*tmp_36 + 0.19601935860219369*tmp_39 + tmp_40;
      real_t tmp_118 = 0.19601935860219369*tmp_32 + 0.19601935860219369*tmp_42 + tmp_43;
      real_t tmp_119 = 0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_45 + tmp_46;
      real_t tmp_120 = 0.037198804536718075*tmp_1*(tmp_117*tmp_27 + tmp_118*tmp_24 + tmp_119*tmp_21 - 1.0/4.0) + 0.037198804536718075*tmp_13*(tmp_117*tmp_25 + tmp_118*tmp_22 + tmp_119*tmp_19 - 1.0/4.0) + 0.037198804536718075*tmp_5*(tmp_117*tmp_26 + tmp_118*tmp_23 + tmp_119*tmp_20 - 1.0/4.0);
      real_t tmp_121 = 0.40446199974765351*tmp_36 + 0.40446199974765351*tmp_39 + tmp_40;
      real_t tmp_122 = 0.40446199974765351*tmp_32 + 0.40446199974765351*tmp_42 + tmp_43;
      real_t tmp_123 = 0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_45 + tmp_46;
      real_t tmp_124 = 0.042507265838595799*tmp_1*(tmp_121*tmp_27 + tmp_122*tmp_24 + tmp_123*tmp_21 - 1.0/4.0) + 0.042507265838595799*tmp_13*(tmp_121*tmp_25 + tmp_122*tmp_22 + tmp_123*tmp_19 - 1.0/4.0) + 0.042507265838595799*tmp_5*(tmp_121*tmp_26 + tmp_122*tmp_23 + tmp_123*tmp_20 - 1.0/4.0);
      real_t tmp_125 = 0.041227165399737475*tmp_36 + 0.1711304259088916*tmp_39 + tmp_40;
      real_t tmp_126 = 0.041227165399737475*tmp_32 + 0.1711304259088916*tmp_42 + tmp_43;
      real_t tmp_127 = 0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_45 + tmp_46;
      real_t tmp_128 = 0.019202922745021479*tmp_1*(tmp_125*tmp_27 + tmp_126*tmp_24 + tmp_127*tmp_21 - 1.0/4.0) + 0.019202922745021479*tmp_13*(tmp_125*tmp_25 + tmp_126*tmp_22 + tmp_127*tmp_19 - 1.0/4.0) + 0.019202922745021479*tmp_5*(tmp_125*tmp_26 + tmp_126*tmp_23 + tmp_127*tmp_20 - 1.0/4.0);
      real_t tmp_129 = tmp_37*(p_affine_13_0*tmp_21 + p_affine_13_1*tmp_24 + p_affine_13_2*tmp_27);
      real_t tmp_130 = tmp_37*(p_affine_13_0*tmp_20 + p_affine_13_1*tmp_23 + p_affine_13_2*tmp_26);
      real_t tmp_131 = tmp_37*(p_affine_13_0*tmp_19 + p_affine_13_1*tmp_22 + p_affine_13_2*tmp_25);
      real_t a_0_0 = -tmp_100*tmp_38 - tmp_104*tmp_38 - tmp_108*tmp_38 - tmp_112*tmp_38 - tmp_116*tmp_38 - tmp_120*tmp_38 - tmp_124*tmp_38 - tmp_128*tmp_38 - tmp_38*tmp_48 - tmp_38*tmp_52 - tmp_38*tmp_56 - tmp_38*tmp_60 - tmp_38*tmp_64 - tmp_38*tmp_68 - tmp_38*tmp_72 - tmp_38*tmp_76 - tmp_38*tmp_80 - tmp_38*tmp_84 - tmp_38*tmp_88 - tmp_38*tmp_92 - tmp_38*tmp_96;
      real_t a_0_1 = -tmp_100*tmp_129 - tmp_104*tmp_129 - tmp_108*tmp_129 - tmp_112*tmp_129 - tmp_116*tmp_129 - tmp_120*tmp_129 - tmp_124*tmp_129 - tmp_128*tmp_129 - tmp_129*tmp_48 - tmp_129*tmp_52 - tmp_129*tmp_56 - tmp_129*tmp_60 - tmp_129*tmp_64 - tmp_129*tmp_68 - tmp_129*tmp_72 - tmp_129*tmp_76 - tmp_129*tmp_80 - tmp_129*tmp_84 - tmp_129*tmp_88 - tmp_129*tmp_92 - tmp_129*tmp_96;
      real_t a_0_2 = -tmp_100*tmp_130 - tmp_104*tmp_130 - tmp_108*tmp_130 - tmp_112*tmp_130 - tmp_116*tmp_130 - tmp_120*tmp_130 - tmp_124*tmp_130 - tmp_128*tmp_130 - tmp_130*tmp_48 - tmp_130*tmp_52 - tmp_130*tmp_56 - tmp_130*tmp_60 - tmp_130*tmp_64 - tmp_130*tmp_68 - tmp_130*tmp_72 - tmp_130*tmp_76 - tmp_130*tmp_80 - tmp_130*tmp_84 - tmp_130*tmp_88 - tmp_130*tmp_92 - tmp_130*tmp_96;
      real_t a_0_3 = -tmp_100*tmp_131 - tmp_104*tmp_131 - tmp_108*tmp_131 - tmp_112*tmp_131 - tmp_116*tmp_131 - tmp_120*tmp_131 - tmp_124*tmp_131 - tmp_128*tmp_131 - tmp_131*tmp_48 - tmp_131*tmp_52 - tmp_131*tmp_56 - tmp_131*tmp_60 - tmp_131*tmp_64 - tmp_131*tmp_68 - tmp_131*tmp_72 - tmp_131*tmp_76 - tmp_131*tmp_80 - tmp_131*tmp_84 - tmp_131*tmp_88 - tmp_131*tmp_92 - tmp_131*tmp_96;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }


};




class EGNIPGVectorLaplaceFormEP1_2 : public hyteg::dg::DGForm
{
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
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
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
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
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
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
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
      real_t tmp_22 = tmp_10*tmp_21 + tmp_13*tmp_20 + tmp_19*tmp_9;
      real_t tmp_23 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_5);
      real_t tmp_24 = tmp_18*(tmp_1*tmp_9 - tmp_10*tmp_14);
      real_t tmp_25 = tmp_18*(tmp_13*tmp_14 - tmp_5*tmp_9);
      real_t tmp_26 = tmp_10*tmp_25 + tmp_13*tmp_24 + tmp_23*tmp_9;
      real_t tmp_27 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_6);
      real_t tmp_28 = tmp_18*(tmp_10*tmp_11 - tmp_6*tmp_9);
      real_t tmp_29 = tmp_18*(-tmp_11*tmp_13 + tmp_3*tmp_9);
      real_t tmp_30 = tmp_10*tmp_29 + tmp_13*tmp_28 + tmp_27*tmp_9;
      real_t tmp_31 = p_affine_0_0*p_affine_1_1;
      real_t tmp_32 = p_affine_0_0*p_affine_1_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_2;
      real_t tmp_34 = p_affine_0_1*p_affine_1_0;
      real_t tmp_35 = p_affine_0_1*p_affine_1_2;
      real_t tmp_36 = p_affine_2_2*p_affine_3_0;
      real_t tmp_37 = p_affine_0_2*p_affine_1_0;
      real_t tmp_38 = p_affine_0_2*p_affine_1_1;
      real_t tmp_39 = p_affine_2_0*p_affine_3_1;
      real_t tmp_40 = p_affine_2_2*p_affine_3_1;
      real_t tmp_41 = p_affine_2_0*p_affine_3_2;
      real_t tmp_42 = p_affine_2_1*p_affine_3_0;
      real_t tmp_43 = std::abs(p_affine_0_0*tmp_33 - p_affine_0_0*tmp_40 + p_affine_0_1*tmp_36 - p_affine_0_1*tmp_41 + p_affine_0_2*tmp_39 - p_affine_0_2*tmp_42 - p_affine_1_0*tmp_33 + p_affine_1_0*tmp_40 - p_affine_1_1*tmp_36 + p_affine_1_1*tmp_41 - p_affine_1_2*tmp_39 + p_affine_1_2*tmp_42 + p_affine_2_0*tmp_35 - p_affine_2_0*tmp_38 - p_affine_2_1*tmp_32 + p_affine_2_1*tmp_37 + p_affine_2_2*tmp_31 - p_affine_2_2*tmp_34 - p_affine_3_0*tmp_35 + p_affine_3_0*tmp_38 + p_affine_3_1*tmp_32 - p_affine_3_1*tmp_37 - p_affine_3_2*tmp_31 + p_affine_3_2*tmp_34);
      real_t tmp_44 = tmp_43*(tmp_22*(-tmp_19 - tmp_20 - tmp_21) + tmp_26*(-tmp_23 - tmp_24 - tmp_25) + tmp_30*(-tmp_27 - tmp_28 - tmp_29));
      real_t tmp_45 = tmp_43*(tmp_21*tmp_22 + tmp_25*tmp_26 + tmp_29*tmp_30);
      real_t tmp_46 = tmp_43*(tmp_20*tmp_22 + tmp_24*tmp_26 + tmp_28*tmp_30);
      real_t tmp_47 = tmp_43*(tmp_19*tmp_22 + tmp_23*tmp_26 + tmp_27*tmp_30);
      real_t a_0_0 = 0.1666666666666668*tmp_44;
      real_t a_0_1 = 0.1666666666666668*tmp_45;
      real_t a_0_2 = 0.1666666666666668*tmp_46;
      real_t a_0_3 = 0.1666666666666668*tmp_47;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
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

         real_t tmp_0 = -p_affine_0_2;
      real_t tmp_1 = p_affine_1_2 + tmp_0;
      real_t tmp_2 = -p_affine_8_2;
      real_t tmp_3 = p_affine_9_2 + tmp_2;
      real_t tmp_4 = p_affine_10_2 + tmp_2;
      real_t tmp_5 = p_affine_8_2 + tmp_0;
      real_t tmp_6 = 0.031405749086161582*tmp_3 + 0.93718850182767688*tmp_4 + tmp_5;
      real_t tmp_7 = -p_affine_0_0;
      real_t tmp_8 = p_affine_2_0 + tmp_7;
      real_t tmp_9 = -p_affine_0_1;
      real_t tmp_10 = p_affine_3_1 + tmp_9;
      real_t tmp_11 = p_affine_3_0 + tmp_7;
      real_t tmp_12 = p_affine_2_1 + tmp_9;
      real_t tmp_13 = p_affine_1_0 + tmp_7;
      real_t tmp_14 = p_affine_3_2 + tmp_0;
      real_t tmp_15 = tmp_12*tmp_14;
      real_t tmp_16 = tmp_1*tmp_10;
      real_t tmp_17 = p_affine_1_1 + tmp_9;
      real_t tmp_18 = p_affine_2_2 + tmp_0;
      real_t tmp_19 = tmp_17*tmp_18;
      real_t tmp_20 = tmp_10*tmp_18;
      real_t tmp_21 = tmp_14*tmp_17;
      real_t tmp_22 = tmp_1*tmp_12;
      real_t tmp_23 = 1.0 / (tmp_11*tmp_19 - tmp_11*tmp_22 + tmp_13*tmp_15 - tmp_13*tmp_20 + tmp_16*tmp_8 - tmp_21*tmp_8);
      real_t tmp_24 = tmp_23*(tmp_10*tmp_8 - tmp_11*tmp_12);
      real_t tmp_25 = tmp_24*tmp_6;
      real_t tmp_26 = -p_affine_8_1;
      real_t tmp_27 = p_affine_9_1 + tmp_26;
      real_t tmp_28 = p_affine_10_1 + tmp_26;
      real_t tmp_29 = p_affine_8_1 + tmp_9;
      real_t tmp_30 = 0.031405749086161582*tmp_27 + 0.93718850182767688*tmp_28 + tmp_29;
      real_t tmp_31 = tmp_23*(tmp_11*tmp_18 - tmp_14*tmp_8);
      real_t tmp_32 = tmp_30*tmp_31;
      real_t tmp_33 = -p_affine_8_0;
      real_t tmp_34 = p_affine_9_0 + tmp_33;
      real_t tmp_35 = p_affine_10_0 + tmp_33;
      real_t tmp_36 = p_affine_8_0 + tmp_7;
      real_t tmp_37 = 0.031405749086161582*tmp_34 + 0.93718850182767688*tmp_35 + tmp_36;
      real_t tmp_38 = tmp_23*(tmp_15 - tmp_20);
      real_t tmp_39 = tmp_37*tmp_38;
      real_t tmp_40 = tmp_25 + tmp_32 + tmp_39;
      real_t tmp_41 = tmp_23*(-tmp_10*tmp_13 + tmp_11*tmp_17);
      real_t tmp_42 = tmp_41*tmp_6;
      real_t tmp_43 = tmp_23*(-tmp_1*tmp_11 + tmp_13*tmp_14);
      real_t tmp_44 = tmp_30*tmp_43;
      real_t tmp_45 = tmp_23*(tmp_16 - tmp_21);
      real_t tmp_46 = tmp_37*tmp_45;
      real_t tmp_47 = tmp_42 + tmp_44 + tmp_46;
      real_t tmp_48 = tmp_23*(tmp_12*tmp_13 - tmp_17*tmp_8);
      real_t tmp_49 = tmp_48*tmp_6;
      real_t tmp_50 = tmp_23*(tmp_1*tmp_8 - tmp_13*tmp_18);
      real_t tmp_51 = tmp_30*tmp_50;
      real_t tmp_52 = tmp_23*(tmp_19 - tmp_22);
      real_t tmp_53 = tmp_37*tmp_52;
      real_t tmp_54 = tmp_49 + tmp_51 + tmp_53;
      real_t tmp_55 = tmp_1*(tmp_40 - 1.0/4.0) + tmp_14*(tmp_54 - 1.0/4.0) + tmp_18*(tmp_47 - 1.0/4.0);
      real_t tmp_56 = 0.5*p_affine_13_0*(-tmp_38 - tmp_45 - tmp_52) + 0.5*p_affine_13_1*(-tmp_31 - tmp_43 - tmp_50) + 0.5*p_affine_13_2*(-tmp_24 - tmp_41 - tmp_48);
      real_t tmp_57 = -tmp_25 - tmp_32 - tmp_39 - tmp_42 - tmp_44 - tmp_46 - tmp_49 - tmp_51 - tmp_53 + 1;
      real_t tmp_58 = 0.5*p_affine_13_0*(tmp_1*tmp_38 + tmp_14*tmp_52 + tmp_18*tmp_45) + 0.5*p_affine_13_1*(tmp_1*tmp_31 + tmp_14*tmp_50 + tmp_18*tmp_43) + 0.5*p_affine_13_2*(tmp_1*tmp_24 + tmp_14*tmp_48 + tmp_18*tmp_41);
      real_t tmp_59 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_60 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_61 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_62 = (std::abs(tmp_28*tmp_60 - tmp_35*tmp_59)*std::abs(tmp_28*tmp_60 - tmp_35*tmp_59)) + (std::abs(tmp_28*tmp_61 - tmp_4*tmp_59)*std::abs(tmp_28*tmp_61 - tmp_4*tmp_59)) + (std::abs(tmp_35*tmp_61 - tmp_4*tmp_60)*std::abs(tmp_35*tmp_61 - tmp_4*tmp_60));
      real_t tmp_63 = 5.0*std::pow(tmp_62, -0.25);
      real_t tmp_64 = tmp_55*tmp_63;
      real_t tmp_65 = 1.0*std::pow(tmp_62, 1.0/2.0);
      real_t tmp_66 = 0.0068572537431980923*tmp_65;
      real_t tmp_67 = 0.19601935860219369*tmp_3 + 0.60796128279561268*tmp_4 + tmp_5;
      real_t tmp_68 = tmp_24*tmp_67;
      real_t tmp_69 = 0.19601935860219369*tmp_27 + 0.60796128279561268*tmp_28 + tmp_29;
      real_t tmp_70 = tmp_31*tmp_69;
      real_t tmp_71 = 0.19601935860219369*tmp_34 + 0.60796128279561268*tmp_35 + tmp_36;
      real_t tmp_72 = tmp_38*tmp_71;
      real_t tmp_73 = tmp_68 + tmp_70 + tmp_72;
      real_t tmp_74 = tmp_41*tmp_67;
      real_t tmp_75 = tmp_43*tmp_69;
      real_t tmp_76 = tmp_45*tmp_71;
      real_t tmp_77 = tmp_74 + tmp_75 + tmp_76;
      real_t tmp_78 = tmp_48*tmp_67;
      real_t tmp_79 = tmp_50*tmp_69;
      real_t tmp_80 = tmp_52*tmp_71;
      real_t tmp_81 = tmp_78 + tmp_79 + tmp_80;
      real_t tmp_82 = tmp_1*(tmp_73 - 1.0/4.0) + tmp_14*(tmp_81 - 1.0/4.0) + tmp_18*(tmp_77 - 1.0/4.0);
      real_t tmp_83 = -tmp_68 - tmp_70 - tmp_72 - tmp_74 - tmp_75 - tmp_76 - tmp_78 - tmp_79 - tmp_80 + 1;
      real_t tmp_84 = tmp_63*tmp_82;
      real_t tmp_85 = 0.037198804536718075*tmp_65;
      real_t tmp_86 = 0.37605877282253791*tmp_3 + 0.039308471900058539*tmp_4 + tmp_5;
      real_t tmp_87 = tmp_24*tmp_86;
      real_t tmp_88 = 0.37605877282253791*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29;
      real_t tmp_89 = tmp_31*tmp_88;
      real_t tmp_90 = 0.37605877282253791*tmp_34 + 0.039308471900058539*tmp_35 + tmp_36;
      real_t tmp_91 = tmp_38*tmp_90;
      real_t tmp_92 = tmp_87 + tmp_89 + tmp_91;
      real_t tmp_93 = tmp_41*tmp_86;
      real_t tmp_94 = tmp_43*tmp_88;
      real_t tmp_95 = tmp_45*tmp_90;
      real_t tmp_96 = tmp_93 + tmp_94 + tmp_95;
      real_t tmp_97 = tmp_48*tmp_86;
      real_t tmp_98 = tmp_50*tmp_88;
      real_t tmp_99 = tmp_52*tmp_90;
      real_t tmp_100 = tmp_97 + tmp_98 + tmp_99;
      real_t tmp_101 = tmp_1*(tmp_92 - 1.0/4.0) + tmp_14*(tmp_100 - 1.0/4.0) + tmp_18*(tmp_96 - 1.0/4.0);
      real_t tmp_102 = -tmp_87 - tmp_89 - tmp_91 - tmp_93 - tmp_94 - tmp_95 - tmp_97 - tmp_98 - tmp_99 + 1;
      real_t tmp_103 = tmp_101*tmp_63;
      real_t tmp_104 = 0.020848748529055869*tmp_65;
      real_t tmp_105 = 0.78764240869137092*tmp_3 + 0.1711304259088916*tmp_4 + tmp_5;
      real_t tmp_106 = tmp_105*tmp_24;
      real_t tmp_107 = 0.78764240869137092*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29;
      real_t tmp_108 = tmp_107*tmp_31;
      real_t tmp_109 = 0.78764240869137092*tmp_34 + 0.1711304259088916*tmp_35 + tmp_36;
      real_t tmp_110 = tmp_109*tmp_38;
      real_t tmp_111 = tmp_106 + tmp_108 + tmp_110;
      real_t tmp_112 = tmp_105*tmp_41;
      real_t tmp_113 = tmp_107*tmp_43;
      real_t tmp_114 = tmp_109*tmp_45;
      real_t tmp_115 = tmp_112 + tmp_113 + tmp_114;
      real_t tmp_116 = tmp_105*tmp_48;
      real_t tmp_117 = tmp_107*tmp_50;
      real_t tmp_118 = tmp_109*tmp_52;
      real_t tmp_119 = tmp_116 + tmp_117 + tmp_118;
      real_t tmp_120 = tmp_1*(tmp_111 - 1.0/4.0) + tmp_14*(tmp_119 - 1.0/4.0) + tmp_18*(tmp_115 - 1.0/4.0);
      real_t tmp_121 = -tmp_106 - tmp_108 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_116 - tmp_117 - tmp_118 + 1;
      real_t tmp_122 = tmp_120*tmp_63;
      real_t tmp_123 = 0.019202922745021479*tmp_65;
      real_t tmp_124 = 0.58463275527740355*tmp_3 + 0.37605877282253791*tmp_4 + tmp_5;
      real_t tmp_125 = tmp_124*tmp_24;
      real_t tmp_126 = 0.58463275527740355*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29;
      real_t tmp_127 = tmp_126*tmp_31;
      real_t tmp_128 = 0.58463275527740355*tmp_34 + 0.37605877282253791*tmp_35 + tmp_36;
      real_t tmp_129 = tmp_128*tmp_38;
      real_t tmp_130 = tmp_125 + tmp_127 + tmp_129;
      real_t tmp_131 = tmp_124*tmp_41;
      real_t tmp_132 = tmp_126*tmp_43;
      real_t tmp_133 = tmp_128*tmp_45;
      real_t tmp_134 = tmp_131 + tmp_132 + tmp_133;
      real_t tmp_135 = tmp_124*tmp_48;
      real_t tmp_136 = tmp_126*tmp_50;
      real_t tmp_137 = tmp_128*tmp_52;
      real_t tmp_138 = tmp_135 + tmp_136 + tmp_137;
      real_t tmp_139 = tmp_1*(tmp_130 - 1.0/4.0) + tmp_14*(tmp_138 - 1.0/4.0) + tmp_18*(tmp_134 - 1.0/4.0);
      real_t tmp_140 = -tmp_125 - tmp_127 - tmp_129 - tmp_131 - tmp_132 - tmp_133 - tmp_135 - tmp_136 - tmp_137 + 1;
      real_t tmp_141 = tmp_139*tmp_63;
      real_t tmp_142 = 0.020848748529055869*tmp_65;
      real_t tmp_143 = 0.041227165399737475*tmp_3 + 0.78764240869137092*tmp_4 + tmp_5;
      real_t tmp_144 = tmp_143*tmp_24;
      real_t tmp_145 = 0.041227165399737475*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29;
      real_t tmp_146 = tmp_145*tmp_31;
      real_t tmp_147 = 0.041227165399737475*tmp_34 + 0.78764240869137092*tmp_35 + tmp_36;
      real_t tmp_148 = tmp_147*tmp_38;
      real_t tmp_149 = tmp_144 + tmp_146 + tmp_148;
      real_t tmp_150 = tmp_143*tmp_41;
      real_t tmp_151 = tmp_145*tmp_43;
      real_t tmp_152 = tmp_147*tmp_45;
      real_t tmp_153 = tmp_150 + tmp_151 + tmp_152;
      real_t tmp_154 = tmp_143*tmp_48;
      real_t tmp_155 = tmp_145*tmp_50;
      real_t tmp_156 = tmp_147*tmp_52;
      real_t tmp_157 = tmp_154 + tmp_155 + tmp_156;
      real_t tmp_158 = tmp_1*(tmp_149 - 1.0/4.0) + tmp_14*(tmp_157 - 1.0/4.0) + tmp_18*(tmp_153 - 1.0/4.0);
      real_t tmp_159 = -tmp_144 - tmp_146 - tmp_148 - tmp_150 - tmp_151 - tmp_152 - tmp_154 - tmp_155 - tmp_156 + 1;
      real_t tmp_160 = tmp_158*tmp_63;
      real_t tmp_161 = 0.019202922745021479*tmp_65;
      real_t tmp_162 = 0.039308471900058539*tmp_3 + 0.58463275527740355*tmp_4 + tmp_5;
      real_t tmp_163 = tmp_162*tmp_24;
      real_t tmp_164 = 0.039308471900058539*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29;
      real_t tmp_165 = tmp_164*tmp_31;
      real_t tmp_166 = 0.039308471900058539*tmp_34 + 0.58463275527740355*tmp_35 + tmp_36;
      real_t tmp_167 = tmp_166*tmp_38;
      real_t tmp_168 = tmp_163 + tmp_165 + tmp_167;
      real_t tmp_169 = tmp_162*tmp_41;
      real_t tmp_170 = tmp_164*tmp_43;
      real_t tmp_171 = tmp_166*tmp_45;
      real_t tmp_172 = tmp_169 + tmp_170 + tmp_171;
      real_t tmp_173 = tmp_162*tmp_48;
      real_t tmp_174 = tmp_164*tmp_50;
      real_t tmp_175 = tmp_166*tmp_52;
      real_t tmp_176 = tmp_173 + tmp_174 + tmp_175;
      real_t tmp_177 = tmp_1*(tmp_168 - 1.0/4.0) + tmp_14*(tmp_176 - 1.0/4.0) + tmp_18*(tmp_172 - 1.0/4.0);
      real_t tmp_178 = -tmp_163 - tmp_165 - tmp_167 - tmp_169 - tmp_170 - tmp_171 - tmp_173 - tmp_174 - tmp_175 + 1;
      real_t tmp_179 = tmp_177*tmp_63;
      real_t tmp_180 = 0.020848748529055869*tmp_65;
      real_t tmp_181 = 0.78764240869137092*tmp_3 + 0.041227165399737475*tmp_4 + tmp_5;
      real_t tmp_182 = tmp_181*tmp_24;
      real_t tmp_183 = 0.78764240869137092*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29;
      real_t tmp_184 = tmp_183*tmp_31;
      real_t tmp_185 = 0.78764240869137092*tmp_34 + 0.041227165399737475*tmp_35 + tmp_36;
      real_t tmp_186 = tmp_185*tmp_38;
      real_t tmp_187 = tmp_182 + tmp_184 + tmp_186;
      real_t tmp_188 = tmp_181*tmp_41;
      real_t tmp_189 = tmp_183*tmp_43;
      real_t tmp_190 = tmp_185*tmp_45;
      real_t tmp_191 = tmp_188 + tmp_189 + tmp_190;
      real_t tmp_192 = tmp_181*tmp_48;
      real_t tmp_193 = tmp_183*tmp_50;
      real_t tmp_194 = tmp_185*tmp_52;
      real_t tmp_195 = tmp_192 + tmp_193 + tmp_194;
      real_t tmp_196 = tmp_1*(tmp_187 - 1.0/4.0) + tmp_14*(tmp_195 - 1.0/4.0) + tmp_18*(tmp_191 - 1.0/4.0);
      real_t tmp_197 = -tmp_182 - tmp_184 - tmp_186 - tmp_188 - tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 + 1;
      real_t tmp_198 = tmp_196*tmp_63;
      real_t tmp_199 = 0.019202922745021479*tmp_65;
      real_t tmp_200 = 0.58463275527740355*tmp_3 + 0.039308471900058539*tmp_4 + tmp_5;
      real_t tmp_201 = tmp_200*tmp_24;
      real_t tmp_202 = 0.58463275527740355*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29;
      real_t tmp_203 = tmp_202*tmp_31;
      real_t tmp_204 = 0.58463275527740355*tmp_34 + 0.039308471900058539*tmp_35 + tmp_36;
      real_t tmp_205 = tmp_204*tmp_38;
      real_t tmp_206 = tmp_201 + tmp_203 + tmp_205;
      real_t tmp_207 = tmp_200*tmp_41;
      real_t tmp_208 = tmp_202*tmp_43;
      real_t tmp_209 = tmp_204*tmp_45;
      real_t tmp_210 = tmp_207 + tmp_208 + tmp_209;
      real_t tmp_211 = tmp_200*tmp_48;
      real_t tmp_212 = tmp_202*tmp_50;
      real_t tmp_213 = tmp_204*tmp_52;
      real_t tmp_214 = tmp_211 + tmp_212 + tmp_213;
      real_t tmp_215 = tmp_1*(tmp_206 - 1.0/4.0) + tmp_14*(tmp_214 - 1.0/4.0) + tmp_18*(tmp_210 - 1.0/4.0);
      real_t tmp_216 = -tmp_201 - tmp_203 - tmp_205 - tmp_207 - tmp_208 - tmp_209 - tmp_211 - tmp_212 - tmp_213 + 1;
      real_t tmp_217 = tmp_215*tmp_63;
      real_t tmp_218 = 0.020848748529055869*tmp_65;
      real_t tmp_219 = 0.1711304259088916*tmp_3 + 0.78764240869137092*tmp_4 + tmp_5;
      real_t tmp_220 = tmp_219*tmp_24;
      real_t tmp_221 = 0.1711304259088916*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29;
      real_t tmp_222 = tmp_221*tmp_31;
      real_t tmp_223 = 0.1711304259088916*tmp_34 + 0.78764240869137092*tmp_35 + tmp_36;
      real_t tmp_224 = tmp_223*tmp_38;
      real_t tmp_225 = tmp_220 + tmp_222 + tmp_224;
      real_t tmp_226 = tmp_219*tmp_41;
      real_t tmp_227 = tmp_221*tmp_43;
      real_t tmp_228 = tmp_223*tmp_45;
      real_t tmp_229 = tmp_226 + tmp_227 + tmp_228;
      real_t tmp_230 = tmp_219*tmp_48;
      real_t tmp_231 = tmp_221*tmp_50;
      real_t tmp_232 = tmp_223*tmp_52;
      real_t tmp_233 = tmp_230 + tmp_231 + tmp_232;
      real_t tmp_234 = tmp_1*(tmp_225 - 1.0/4.0) + tmp_14*(tmp_233 - 1.0/4.0) + tmp_18*(tmp_229 - 1.0/4.0);
      real_t tmp_235 = -tmp_220 - tmp_222 - tmp_224 - tmp_226 - tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 + 1;
      real_t tmp_236 = tmp_234*tmp_63;
      real_t tmp_237 = 0.019202922745021479*tmp_65;
      real_t tmp_238 = 0.37605877282253791*tmp_3 + 0.58463275527740355*tmp_4 + tmp_5;
      real_t tmp_239 = tmp_238*tmp_24;
      real_t tmp_240 = 0.37605877282253791*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29;
      real_t tmp_241 = tmp_240*tmp_31;
      real_t tmp_242 = 0.37605877282253791*tmp_34 + 0.58463275527740355*tmp_35 + tmp_36;
      real_t tmp_243 = tmp_242*tmp_38;
      real_t tmp_244 = tmp_239 + tmp_241 + tmp_243;
      real_t tmp_245 = tmp_238*tmp_41;
      real_t tmp_246 = tmp_240*tmp_43;
      real_t tmp_247 = tmp_242*tmp_45;
      real_t tmp_248 = tmp_245 + tmp_246 + tmp_247;
      real_t tmp_249 = tmp_238*tmp_48;
      real_t tmp_250 = tmp_240*tmp_50;
      real_t tmp_251 = tmp_242*tmp_52;
      real_t tmp_252 = tmp_249 + tmp_250 + tmp_251;
      real_t tmp_253 = tmp_1*(tmp_244 - 1.0/4.0) + tmp_14*(tmp_252 - 1.0/4.0) + tmp_18*(tmp_248 - 1.0/4.0);
      real_t tmp_254 = -tmp_239 - tmp_241 - tmp_243 - tmp_245 - tmp_246 - tmp_247 - tmp_249 - tmp_250 - tmp_251 + 1;
      real_t tmp_255 = tmp_253*tmp_63;
      real_t tmp_256 = 0.020848748529055869*tmp_65;
      real_t tmp_257 = 0.041227165399737475*tmp_3 + 0.1711304259088916*tmp_4 + tmp_5;
      real_t tmp_258 = tmp_24*tmp_257;
      real_t tmp_259 = 0.041227165399737475*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29;
      real_t tmp_260 = tmp_259*tmp_31;
      real_t tmp_261 = 0.041227165399737475*tmp_34 + 0.1711304259088916*tmp_35 + tmp_36;
      real_t tmp_262 = tmp_261*tmp_38;
      real_t tmp_263 = tmp_258 + tmp_260 + tmp_262;
      real_t tmp_264 = tmp_257*tmp_41;
      real_t tmp_265 = tmp_259*tmp_43;
      real_t tmp_266 = tmp_261*tmp_45;
      real_t tmp_267 = tmp_264 + tmp_265 + tmp_266;
      real_t tmp_268 = tmp_257*tmp_48;
      real_t tmp_269 = tmp_259*tmp_50;
      real_t tmp_270 = tmp_261*tmp_52;
      real_t tmp_271 = tmp_268 + tmp_269 + tmp_270;
      real_t tmp_272 = tmp_1*(tmp_263 - 1.0/4.0) + tmp_14*(tmp_271 - 1.0/4.0) + tmp_18*(tmp_267 - 1.0/4.0);
      real_t tmp_273 = -tmp_258 - tmp_260 - tmp_262 - tmp_264 - tmp_265 - tmp_266 - tmp_268 - tmp_269 - tmp_270 + 1;
      real_t tmp_274 = tmp_272*tmp_63;
      real_t tmp_275 = 0.019202922745021479*tmp_65;
      real_t tmp_276 = 0.40446199974765351*tmp_3 + 0.19107600050469298*tmp_4 + tmp_5;
      real_t tmp_277 = tmp_24*tmp_276;
      real_t tmp_278 = 0.40446199974765351*tmp_27 + 0.19107600050469298*tmp_28 + tmp_29;
      real_t tmp_279 = tmp_278*tmp_31;
      real_t tmp_280 = 0.40446199974765351*tmp_34 + 0.19107600050469298*tmp_35 + tmp_36;
      real_t tmp_281 = tmp_280*tmp_38;
      real_t tmp_282 = tmp_277 + tmp_279 + tmp_281;
      real_t tmp_283 = tmp_276*tmp_41;
      real_t tmp_284 = tmp_278*tmp_43;
      real_t tmp_285 = tmp_280*tmp_45;
      real_t tmp_286 = tmp_283 + tmp_284 + tmp_285;
      real_t tmp_287 = tmp_276*tmp_48;
      real_t tmp_288 = tmp_278*tmp_50;
      real_t tmp_289 = tmp_280*tmp_52;
      real_t tmp_290 = tmp_287 + tmp_288 + tmp_289;
      real_t tmp_291 = tmp_1*(tmp_282 - 1.0/4.0) + tmp_14*(tmp_290 - 1.0/4.0) + tmp_18*(tmp_286 - 1.0/4.0);
      real_t tmp_292 = -tmp_277 - tmp_279 - tmp_281 - tmp_283 - tmp_284 - tmp_285 - tmp_287 - tmp_288 - tmp_289 + 1;
      real_t tmp_293 = tmp_291*tmp_63;
      real_t tmp_294 = 0.042507265838595799*tmp_65;
      real_t tmp_295 = 0.039308471900058539*tmp_3 + 0.37605877282253791*tmp_4 + tmp_5;
      real_t tmp_296 = tmp_24*tmp_295;
      real_t tmp_297 = 0.039308471900058539*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29;
      real_t tmp_298 = tmp_297*tmp_31;
      real_t tmp_299 = 0.039308471900058539*tmp_34 + 0.37605877282253791*tmp_35 + tmp_36;
      real_t tmp_300 = tmp_299*tmp_38;
      real_t tmp_301 = tmp_296 + tmp_298 + tmp_300;
      real_t tmp_302 = tmp_295*tmp_41;
      real_t tmp_303 = tmp_297*tmp_43;
      real_t tmp_304 = tmp_299*tmp_45;
      real_t tmp_305 = tmp_302 + tmp_303 + tmp_304;
      real_t tmp_306 = tmp_295*tmp_48;
      real_t tmp_307 = tmp_297*tmp_50;
      real_t tmp_308 = tmp_299*tmp_52;
      real_t tmp_309 = tmp_306 + tmp_307 + tmp_308;
      real_t tmp_310 = tmp_1*(tmp_301 - 1.0/4.0) + tmp_14*(tmp_309 - 1.0/4.0) + tmp_18*(tmp_305 - 1.0/4.0);
      real_t tmp_311 = -tmp_296 - tmp_298 - tmp_300 - tmp_302 - tmp_303 - tmp_304 - tmp_306 - tmp_307 - tmp_308 + 1;
      real_t tmp_312 = tmp_310*tmp_63;
      real_t tmp_313 = 0.020848748529055869*tmp_65;
      real_t tmp_314 = 0.93718850182767688*tmp_3 + 0.031405749086161582*tmp_4 + tmp_5;
      real_t tmp_315 = tmp_24*tmp_314;
      real_t tmp_316 = 0.93718850182767688*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29;
      real_t tmp_317 = tmp_31*tmp_316;
      real_t tmp_318 = 0.93718850182767688*tmp_34 + 0.031405749086161582*tmp_35 + tmp_36;
      real_t tmp_319 = tmp_318*tmp_38;
      real_t tmp_320 = tmp_315 + tmp_317 + tmp_319;
      real_t tmp_321 = tmp_314*tmp_41;
      real_t tmp_322 = tmp_316*tmp_43;
      real_t tmp_323 = tmp_318*tmp_45;
      real_t tmp_324 = tmp_321 + tmp_322 + tmp_323;
      real_t tmp_325 = tmp_314*tmp_48;
      real_t tmp_326 = tmp_316*tmp_50;
      real_t tmp_327 = tmp_318*tmp_52;
      real_t tmp_328 = tmp_325 + tmp_326 + tmp_327;
      real_t tmp_329 = tmp_1*(tmp_320 - 1.0/4.0) + tmp_14*(tmp_328 - 1.0/4.0) + tmp_18*(tmp_324 - 1.0/4.0);
      real_t tmp_330 = -tmp_315 - tmp_317 - tmp_319 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1;
      real_t tmp_331 = tmp_329*tmp_63;
      real_t tmp_332 = 0.0068572537431980923*tmp_65;
      real_t tmp_333 = 0.60796128279561268*tmp_3 + 0.19601935860219369*tmp_4 + tmp_5;
      real_t tmp_334 = tmp_24*tmp_333;
      real_t tmp_335 = 0.60796128279561268*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29;
      real_t tmp_336 = tmp_31*tmp_335;
      real_t tmp_337 = 0.60796128279561268*tmp_34 + 0.19601935860219369*tmp_35 + tmp_36;
      real_t tmp_338 = tmp_337*tmp_38;
      real_t tmp_339 = tmp_334 + tmp_336 + tmp_338;
      real_t tmp_340 = tmp_333*tmp_41;
      real_t tmp_341 = tmp_335*tmp_43;
      real_t tmp_342 = tmp_337*tmp_45;
      real_t tmp_343 = tmp_340 + tmp_341 + tmp_342;
      real_t tmp_344 = tmp_333*tmp_48;
      real_t tmp_345 = tmp_335*tmp_50;
      real_t tmp_346 = tmp_337*tmp_52;
      real_t tmp_347 = tmp_344 + tmp_345 + tmp_346;
      real_t tmp_348 = tmp_1*(tmp_339 - 1.0/4.0) + tmp_14*(tmp_347 - 1.0/4.0) + tmp_18*(tmp_343 - 1.0/4.0);
      real_t tmp_349 = -tmp_334 - tmp_336 - tmp_338 - tmp_340 - tmp_341 - tmp_342 - tmp_344 - tmp_345 - tmp_346 + 1;
      real_t tmp_350 = tmp_348*tmp_63;
      real_t tmp_351 = 0.037198804536718075*tmp_65;
      real_t tmp_352 = 0.19107600050469298*tmp_3 + 0.40446199974765351*tmp_4 + tmp_5;
      real_t tmp_353 = tmp_24*tmp_352;
      real_t tmp_354 = 0.19107600050469298*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29;
      real_t tmp_355 = tmp_31*tmp_354;
      real_t tmp_356 = 0.19107600050469298*tmp_34 + 0.40446199974765351*tmp_35 + tmp_36;
      real_t tmp_357 = tmp_356*tmp_38;
      real_t tmp_358 = tmp_353 + tmp_355 + tmp_357;
      real_t tmp_359 = tmp_352*tmp_41;
      real_t tmp_360 = tmp_354*tmp_43;
      real_t tmp_361 = tmp_356*tmp_45;
      real_t tmp_362 = tmp_359 + tmp_360 + tmp_361;
      real_t tmp_363 = tmp_352*tmp_48;
      real_t tmp_364 = tmp_354*tmp_50;
      real_t tmp_365 = tmp_356*tmp_52;
      real_t tmp_366 = tmp_363 + tmp_364 + tmp_365;
      real_t tmp_367 = tmp_1*(tmp_358 - 1.0/4.0) + tmp_14*(tmp_366 - 1.0/4.0) + tmp_18*(tmp_362 - 1.0/4.0);
      real_t tmp_368 = -tmp_353 - tmp_355 - tmp_357 - tmp_359 - tmp_360 - tmp_361 - tmp_363 - tmp_364 - tmp_365 + 1;
      real_t tmp_369 = tmp_367*tmp_63;
      real_t tmp_370 = 0.042507265838595799*tmp_65;
      real_t tmp_371 = 0.031405749086161582*tmp_3 + 0.031405749086161582*tmp_4 + tmp_5;
      real_t tmp_372 = tmp_24*tmp_371;
      real_t tmp_373 = 0.031405749086161582*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29;
      real_t tmp_374 = tmp_31*tmp_373;
      real_t tmp_375 = 0.031405749086161582*tmp_34 + 0.031405749086161582*tmp_35 + tmp_36;
      real_t tmp_376 = tmp_375*tmp_38;
      real_t tmp_377 = tmp_372 + tmp_374 + tmp_376;
      real_t tmp_378 = tmp_371*tmp_41;
      real_t tmp_379 = tmp_373*tmp_43;
      real_t tmp_380 = tmp_375*tmp_45;
      real_t tmp_381 = tmp_378 + tmp_379 + tmp_380;
      real_t tmp_382 = tmp_371*tmp_48;
      real_t tmp_383 = tmp_373*tmp_50;
      real_t tmp_384 = tmp_375*tmp_52;
      real_t tmp_385 = tmp_382 + tmp_383 + tmp_384;
      real_t tmp_386 = tmp_1*(tmp_377 - 1.0/4.0) + tmp_14*(tmp_385 - 1.0/4.0) + tmp_18*(tmp_381 - 1.0/4.0);
      real_t tmp_387 = -tmp_372 - tmp_374 - tmp_376 - tmp_378 - tmp_379 - tmp_380 - tmp_382 - tmp_383 - tmp_384 + 1;
      real_t tmp_388 = tmp_386*tmp_63;
      real_t tmp_389 = 0.0068572537431980923*tmp_65;
      real_t tmp_390 = 0.19601935860219369*tmp_3 + 0.19601935860219369*tmp_4 + tmp_5;
      real_t tmp_391 = tmp_24*tmp_390;
      real_t tmp_392 = 0.19601935860219369*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29;
      real_t tmp_393 = tmp_31*tmp_392;
      real_t tmp_394 = 0.19601935860219369*tmp_34 + 0.19601935860219369*tmp_35 + tmp_36;
      real_t tmp_395 = tmp_38*tmp_394;
      real_t tmp_396 = tmp_391 + tmp_393 + tmp_395;
      real_t tmp_397 = tmp_390*tmp_41;
      real_t tmp_398 = tmp_392*tmp_43;
      real_t tmp_399 = tmp_394*tmp_45;
      real_t tmp_400 = tmp_397 + tmp_398 + tmp_399;
      real_t tmp_401 = tmp_390*tmp_48;
      real_t tmp_402 = tmp_392*tmp_50;
      real_t tmp_403 = tmp_394*tmp_52;
      real_t tmp_404 = tmp_401 + tmp_402 + tmp_403;
      real_t tmp_405 = tmp_1*(tmp_396 - 1.0/4.0) + tmp_14*(tmp_404 - 1.0/4.0) + tmp_18*(tmp_400 - 1.0/4.0);
      real_t tmp_406 = -tmp_391 - tmp_393 - tmp_395 - tmp_397 - tmp_398 - tmp_399 - tmp_401 - tmp_402 - tmp_403 + 1;
      real_t tmp_407 = tmp_405*tmp_63;
      real_t tmp_408 = 0.037198804536718075*tmp_65;
      real_t tmp_409 = 0.40446199974765351*tmp_3 + 0.40446199974765351*tmp_4 + tmp_5;
      real_t tmp_410 = tmp_24*tmp_409;
      real_t tmp_411 = 0.40446199974765351*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29;
      real_t tmp_412 = tmp_31*tmp_411;
      real_t tmp_413 = 0.40446199974765351*tmp_34 + 0.40446199974765351*tmp_35 + tmp_36;
      real_t tmp_414 = tmp_38*tmp_413;
      real_t tmp_415 = tmp_410 + tmp_412 + tmp_414;
      real_t tmp_416 = tmp_409*tmp_41;
      real_t tmp_417 = tmp_411*tmp_43;
      real_t tmp_418 = tmp_413*tmp_45;
      real_t tmp_419 = tmp_416 + tmp_417 + tmp_418;
      real_t tmp_420 = tmp_409*tmp_48;
      real_t tmp_421 = tmp_411*tmp_50;
      real_t tmp_422 = tmp_413*tmp_52;
      real_t tmp_423 = tmp_420 + tmp_421 + tmp_422;
      real_t tmp_424 = tmp_1*(tmp_415 - 1.0/4.0) + tmp_14*(tmp_423 - 1.0/4.0) + tmp_18*(tmp_419 - 1.0/4.0);
      real_t tmp_425 = -tmp_410 - tmp_412 - tmp_414 - tmp_416 - tmp_417 - tmp_418 - tmp_420 - tmp_421 - tmp_422 + 1;
      real_t tmp_426 = tmp_424*tmp_63;
      real_t tmp_427 = 0.042507265838595799*tmp_65;
      real_t tmp_428 = 0.1711304259088916*tmp_3 + 0.041227165399737475*tmp_4 + tmp_5;
      real_t tmp_429 = tmp_24*tmp_428;
      real_t tmp_430 = 0.1711304259088916*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29;
      real_t tmp_431 = tmp_31*tmp_430;
      real_t tmp_432 = 0.1711304259088916*tmp_34 + 0.041227165399737475*tmp_35 + tmp_36;
      real_t tmp_433 = tmp_38*tmp_432;
      real_t tmp_434 = tmp_429 + tmp_431 + tmp_433;
      real_t tmp_435 = tmp_41*tmp_428;
      real_t tmp_436 = tmp_43*tmp_430;
      real_t tmp_437 = tmp_432*tmp_45;
      real_t tmp_438 = tmp_435 + tmp_436 + tmp_437;
      real_t tmp_439 = tmp_428*tmp_48;
      real_t tmp_440 = tmp_430*tmp_50;
      real_t tmp_441 = tmp_432*tmp_52;
      real_t tmp_442 = tmp_439 + tmp_440 + tmp_441;
      real_t tmp_443 = tmp_1*(tmp_434 - 1.0/4.0) + tmp_14*(tmp_442 - 1.0/4.0) + tmp_18*(tmp_438 - 1.0/4.0);
      real_t tmp_444 = -tmp_429 - tmp_431 - tmp_433 - tmp_435 - tmp_436 - tmp_437 - tmp_439 - tmp_440 - tmp_441 + 1;
      real_t tmp_445 = tmp_443*tmp_63;
      real_t tmp_446 = 0.019202922745021479*tmp_65;
      real_t tmp_447 = 0.5*p_affine_13_0*tmp_38 + 0.5*p_affine_13_1*tmp_31 + 0.5*p_affine_13_2*tmp_24;
      real_t tmp_448 = 0.5*p_affine_13_0*tmp_45 + 0.5*p_affine_13_1*tmp_43 + 0.5*p_affine_13_2*tmp_41;
      real_t tmp_449 = 0.5*p_affine_13_0*tmp_52 + 0.5*p_affine_13_1*tmp_50 + 0.5*p_affine_13_2*tmp_48;
      real_t a_0_0 = tmp_104*(-tmp_101*tmp_56 + tmp_102*tmp_103 - tmp_102*tmp_58) + tmp_123*(-tmp_120*tmp_56 + tmp_121*tmp_122 - tmp_121*tmp_58) + tmp_142*(-tmp_139*tmp_56 + tmp_140*tmp_141 - tmp_140*tmp_58) + tmp_161*(-tmp_158*tmp_56 + tmp_159*tmp_160 - tmp_159*tmp_58) + tmp_180*(-tmp_177*tmp_56 + tmp_178*tmp_179 - tmp_178*tmp_58) + tmp_199*(-tmp_196*tmp_56 + tmp_197*tmp_198 - tmp_197*tmp_58) + tmp_218*(-tmp_215*tmp_56 + tmp_216*tmp_217 - tmp_216*tmp_58) + tmp_237*(-tmp_234*tmp_56 + tmp_235*tmp_236 - tmp_235*tmp_58) + tmp_256*(-tmp_253*tmp_56 + tmp_254*tmp_255 - tmp_254*tmp_58) + tmp_275*(-tmp_272*tmp_56 + tmp_273*tmp_274 - tmp_273*tmp_58) + tmp_294*(-tmp_291*tmp_56 + tmp_292*tmp_293 - tmp_292*tmp_58) + tmp_313*(-tmp_310*tmp_56 + tmp_311*tmp_312 - tmp_311*tmp_58) + tmp_332*(-tmp_329*tmp_56 + tmp_330*tmp_331 - tmp_330*tmp_58) + tmp_351*(-tmp_348*tmp_56 + tmp_349*tmp_350 - tmp_349*tmp_58) + tmp_370*(-tmp_367*tmp_56 + tmp_368*tmp_369 - tmp_368*tmp_58) + tmp_389*(-tmp_386*tmp_56 + tmp_387*tmp_388 - tmp_387*tmp_58) + tmp_408*(-tmp_405*tmp_56 + tmp_406*tmp_407 - tmp_406*tmp_58) + tmp_427*(-tmp_424*tmp_56 + tmp_425*tmp_426 - tmp_425*tmp_58) + tmp_446*(-tmp_443*tmp_56 + tmp_444*tmp_445 - tmp_444*tmp_58) + tmp_66*(-tmp_55*tmp_56 - tmp_57*tmp_58 + tmp_57*tmp_64) + tmp_85*(-tmp_56*tmp_82 - tmp_58*tmp_83 + tmp_83*tmp_84);
      real_t a_0_1 = tmp_104*(-tmp_101*tmp_447 + tmp_103*tmp_92 - tmp_58*tmp_92) + tmp_123*(tmp_111*tmp_122 - tmp_111*tmp_58 - tmp_120*tmp_447) + tmp_142*(tmp_130*tmp_141 - tmp_130*tmp_58 - tmp_139*tmp_447) + tmp_161*(tmp_149*tmp_160 - tmp_149*tmp_58 - tmp_158*tmp_447) + tmp_180*(tmp_168*tmp_179 - tmp_168*tmp_58 - tmp_177*tmp_447) + tmp_199*(tmp_187*tmp_198 - tmp_187*tmp_58 - tmp_196*tmp_447) + tmp_218*(tmp_206*tmp_217 - tmp_206*tmp_58 - tmp_215*tmp_447) + tmp_237*(tmp_225*tmp_236 - tmp_225*tmp_58 - tmp_234*tmp_447) + tmp_256*(tmp_244*tmp_255 - tmp_244*tmp_58 - tmp_253*tmp_447) + tmp_275*(tmp_263*tmp_274 - tmp_263*tmp_58 - tmp_272*tmp_447) + tmp_294*(tmp_282*tmp_293 - tmp_282*tmp_58 - tmp_291*tmp_447) + tmp_313*(tmp_301*tmp_312 - tmp_301*tmp_58 - tmp_310*tmp_447) + tmp_332*(tmp_320*tmp_331 - tmp_320*tmp_58 - tmp_329*tmp_447) + tmp_351*(tmp_339*tmp_350 - tmp_339*tmp_58 - tmp_348*tmp_447) + tmp_370*(tmp_358*tmp_369 - tmp_358*tmp_58 - tmp_367*tmp_447) + tmp_389*(tmp_377*tmp_388 - tmp_377*tmp_58 - tmp_386*tmp_447) + tmp_408*(tmp_396*tmp_407 - tmp_396*tmp_58 - tmp_405*tmp_447) + tmp_427*(tmp_415*tmp_426 - tmp_415*tmp_58 - tmp_424*tmp_447) + tmp_446*(tmp_434*tmp_445 - tmp_434*tmp_58 - tmp_443*tmp_447) + tmp_66*(-tmp_40*tmp_58 + tmp_40*tmp_64 - tmp_447*tmp_55) + tmp_85*(-tmp_447*tmp_82 - tmp_58*tmp_73 + tmp_73*tmp_84);
      real_t a_0_2 = tmp_104*(-tmp_101*tmp_448 + tmp_103*tmp_96 - tmp_58*tmp_96) + tmp_123*(tmp_115*tmp_122 - tmp_115*tmp_58 - tmp_120*tmp_448) + tmp_142*(tmp_134*tmp_141 - tmp_134*tmp_58 - tmp_139*tmp_448) + tmp_161*(tmp_153*tmp_160 - tmp_153*tmp_58 - tmp_158*tmp_448) + tmp_180*(tmp_172*tmp_179 - tmp_172*tmp_58 - tmp_177*tmp_448) + tmp_199*(tmp_191*tmp_198 - tmp_191*tmp_58 - tmp_196*tmp_448) + tmp_218*(tmp_210*tmp_217 - tmp_210*tmp_58 - tmp_215*tmp_448) + tmp_237*(tmp_229*tmp_236 - tmp_229*tmp_58 - tmp_234*tmp_448) + tmp_256*(tmp_248*tmp_255 - tmp_248*tmp_58 - tmp_253*tmp_448) + tmp_275*(tmp_267*tmp_274 - tmp_267*tmp_58 - tmp_272*tmp_448) + tmp_294*(tmp_286*tmp_293 - tmp_286*tmp_58 - tmp_291*tmp_448) + tmp_313*(tmp_305*tmp_312 - tmp_305*tmp_58 - tmp_310*tmp_448) + tmp_332*(tmp_324*tmp_331 - tmp_324*tmp_58 - tmp_329*tmp_448) + tmp_351*(tmp_343*tmp_350 - tmp_343*tmp_58 - tmp_348*tmp_448) + tmp_370*(tmp_362*tmp_369 - tmp_362*tmp_58 - tmp_367*tmp_448) + tmp_389*(tmp_381*tmp_388 - tmp_381*tmp_58 - tmp_386*tmp_448) + tmp_408*(tmp_400*tmp_407 - tmp_400*tmp_58 - tmp_405*tmp_448) + tmp_427*(tmp_419*tmp_426 - tmp_419*tmp_58 - tmp_424*tmp_448) + tmp_446*(tmp_438*tmp_445 - tmp_438*tmp_58 - tmp_443*tmp_448) + tmp_66*(-tmp_448*tmp_55 - tmp_47*tmp_58 + tmp_47*tmp_64) + tmp_85*(-tmp_448*tmp_82 - tmp_58*tmp_77 + tmp_77*tmp_84);
      real_t a_0_3 = tmp_104*(tmp_100*tmp_103 - tmp_100*tmp_58 - tmp_101*tmp_449) + tmp_123*(tmp_119*tmp_122 - tmp_119*tmp_58 - tmp_120*tmp_449) + tmp_142*(tmp_138*tmp_141 - tmp_138*tmp_58 - tmp_139*tmp_449) + tmp_161*(tmp_157*tmp_160 - tmp_157*tmp_58 - tmp_158*tmp_449) + tmp_180*(tmp_176*tmp_179 - tmp_176*tmp_58 - tmp_177*tmp_449) + tmp_199*(tmp_195*tmp_198 - tmp_195*tmp_58 - tmp_196*tmp_449) + tmp_218*(tmp_214*tmp_217 - tmp_214*tmp_58 - tmp_215*tmp_449) + tmp_237*(tmp_233*tmp_236 - tmp_233*tmp_58 - tmp_234*tmp_449) + tmp_256*(tmp_252*tmp_255 - tmp_252*tmp_58 - tmp_253*tmp_449) + tmp_275*(tmp_271*tmp_274 - tmp_271*tmp_58 - tmp_272*tmp_449) + tmp_294*(tmp_290*tmp_293 - tmp_290*tmp_58 - tmp_291*tmp_449) + tmp_313*(tmp_309*tmp_312 - tmp_309*tmp_58 - tmp_310*tmp_449) + tmp_332*(tmp_328*tmp_331 - tmp_328*tmp_58 - tmp_329*tmp_449) + tmp_351*(tmp_347*tmp_350 - tmp_347*tmp_58 - tmp_348*tmp_449) + tmp_370*(tmp_366*tmp_369 - tmp_366*tmp_58 - tmp_367*tmp_449) + tmp_389*(tmp_385*tmp_388 - tmp_385*tmp_58 - tmp_386*tmp_449) + tmp_408*(tmp_404*tmp_407 - tmp_404*tmp_58 - tmp_405*tmp_449) + tmp_427*(tmp_423*tmp_426 - tmp_423*tmp_58 - tmp_424*tmp_449) + tmp_446*(tmp_442*tmp_445 - tmp_442*tmp_58 - tmp_443*tmp_449) + tmp_66*(-tmp_449*tmp_55 - tmp_54*tmp_58 + tmp_54*tmp_64) + tmp_85*(-tmp_449*tmp_82 - tmp_58*tmp_81 + tmp_81*tmp_84);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
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


      real_t tmp_0 = -p_affine_0_2;
      real_t tmp_1 = p_affine_1_2 + tmp_0;
      real_t tmp_2 = -p_affine_0_0;
      real_t tmp_3 = p_affine_2_0 + tmp_2;
      real_t tmp_4 = -p_affine_0_1;
      real_t tmp_5 = p_affine_3_1 + tmp_4;
      real_t tmp_6 = p_affine_3_0 + tmp_2;
      real_t tmp_7 = p_affine_2_1 + tmp_4;
      real_t tmp_8 = tmp_3*tmp_5 - tmp_6*tmp_7;
      real_t tmp_9 = p_affine_1_0 + tmp_2;
      real_t tmp_10 = p_affine_3_2 + tmp_0;
      real_t tmp_11 = tmp_10*tmp_7;
      real_t tmp_12 = tmp_1*tmp_5;
      real_t tmp_13 = p_affine_1_1 + tmp_4;
      real_t tmp_14 = p_affine_2_2 + tmp_0;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_14*tmp_5;
      real_t tmp_17 = tmp_10*tmp_13;
      real_t tmp_18 = tmp_1*tmp_7;
      real_t tmp_19 = 1.0 / (tmp_11*tmp_9 + tmp_12*tmp_3 + tmp_15*tmp_6 - tmp_16*tmp_9 - tmp_17*tmp_3 - tmp_18*tmp_6);
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = 0.031405749086161582*tmp_21 + 0.93718850182767688*tmp_22;
      real_t tmp_24 = p_affine_8_2 + tmp_0;
      real_t tmp_25 = tmp_19*(tmp_23 + tmp_24);
      real_t tmp_26 = -tmp_10*tmp_3 + tmp_14*tmp_6;
      real_t tmp_27 = -p_affine_8_1;
      real_t tmp_28 = p_affine_9_1 + tmp_27;
      real_t tmp_29 = p_affine_10_1 + tmp_27;
      real_t tmp_30 = 0.031405749086161582*tmp_28 + 0.93718850182767688*tmp_29;
      real_t tmp_31 = p_affine_8_1 + tmp_4;
      real_t tmp_32 = tmp_19*(tmp_30 + tmp_31);
      real_t tmp_33 = tmp_11 - tmp_16;
      real_t tmp_34 = -p_affine_8_0;
      real_t tmp_35 = p_affine_9_0 + tmp_34;
      real_t tmp_36 = p_affine_10_0 + tmp_34;
      real_t tmp_37 = 0.031405749086161582*tmp_35 + 0.93718850182767688*tmp_36;
      real_t tmp_38 = p_affine_8_0 + tmp_2;
      real_t tmp_39 = tmp_19*(tmp_37 + tmp_38);
      real_t tmp_40 = tmp_13*tmp_6 - tmp_5*tmp_9;
      real_t tmp_41 = -tmp_1*tmp_6 + tmp_10*tmp_9;
      real_t tmp_42 = tmp_12 - tmp_17;
      real_t tmp_43 = -tmp_13*tmp_3 + tmp_7*tmp_9;
      real_t tmp_44 = tmp_1*tmp_3 - tmp_14*tmp_9;
      real_t tmp_45 = tmp_15 - tmp_18;
      real_t tmp_46 = tmp_1*(tmp_25*tmp_8 + tmp_26*tmp_32 + tmp_33*tmp_39 - 1.0/4.0) + tmp_10*(tmp_25*tmp_43 + tmp_32*tmp_44 + tmp_39*tmp_45 - 1.0/4.0) + tmp_14*(tmp_25*tmp_40 + tmp_32*tmp_41 + tmp_39*tmp_42 - 1.0/4.0);
      real_t tmp_47 = -p_affine_4_1;
      real_t tmp_48 = p_affine_5_1 + tmp_47;
      real_t tmp_49 = -p_affine_4_2;
      real_t tmp_50 = p_affine_6_2 + tmp_49;
      real_t tmp_51 = tmp_48*tmp_50;
      real_t tmp_52 = p_affine_6_1 + tmp_47;
      real_t tmp_53 = p_affine_5_2 + tmp_49;
      real_t tmp_54 = p_affine_7_2 + tmp_49;
      real_t tmp_55 = -p_affine_4_0;
      real_t tmp_56 = p_affine_5_0 + tmp_55;
      real_t tmp_57 = tmp_52*tmp_56;
      real_t tmp_58 = p_affine_7_1 + tmp_47;
      real_t tmp_59 = p_affine_6_0 + tmp_55;
      real_t tmp_60 = tmp_53*tmp_59;
      real_t tmp_61 = p_affine_7_0 + tmp_55;
      real_t tmp_62 = tmp_56*tmp_58;
      real_t tmp_63 = tmp_48*tmp_59;
      real_t tmp_64 = tmp_53*tmp_61;
      real_t tmp_65 = 1.0 / (-tmp_50*tmp_62 + tmp_51*tmp_61 - tmp_52*tmp_64 + tmp_54*tmp_57 - tmp_54*tmp_63 + tmp_58*tmp_60);
      real_t tmp_66 = tmp_65*(tmp_51 - tmp_52*tmp_53);
      real_t tmp_67 = tmp_65*(-tmp_48*tmp_54 + tmp_53*tmp_58);
      real_t tmp_68 = tmp_65*(-tmp_50*tmp_58 + tmp_52*tmp_54);
      real_t tmp_69 = tmp_65*(-tmp_50*tmp_56 + tmp_60);
      real_t tmp_70 = tmp_65*(tmp_54*tmp_56 - tmp_64);
      real_t tmp_71 = tmp_65*(tmp_50*tmp_61 - tmp_54*tmp_59);
      real_t tmp_72 = tmp_65*(tmp_57 - tmp_63);
      real_t tmp_73 = tmp_65*(tmp_48*tmp_61 - tmp_62);
      real_t tmp_74 = tmp_65*(-tmp_52*tmp_61 + tmp_58*tmp_59);
      real_t tmp_75 = 0.5*p_affine_13_0*(-tmp_66 - tmp_67 - tmp_68) + 0.5*p_affine_13_1*(-tmp_69 - tmp_70 - tmp_71) + 0.5*p_affine_13_2*(-tmp_72 - tmp_73 - tmp_74);
      real_t tmp_76 = p_affine_8_2 + tmp_49;
      real_t tmp_77 = tmp_23 + tmp_76;
      real_t tmp_78 = tmp_72*tmp_77;
      real_t tmp_79 = tmp_73*tmp_77;
      real_t tmp_80 = p_affine_8_1 + tmp_47;
      real_t tmp_81 = tmp_30 + tmp_80;
      real_t tmp_82 = tmp_69*tmp_81;
      real_t tmp_83 = tmp_70*tmp_81;
      real_t tmp_84 = tmp_74*tmp_77;
      real_t tmp_85 = tmp_71*tmp_81;
      real_t tmp_86 = p_affine_8_0 + tmp_55;
      real_t tmp_87 = tmp_37 + tmp_86;
      real_t tmp_88 = tmp_66*tmp_87;
      real_t tmp_89 = tmp_67*tmp_87;
      real_t tmp_90 = tmp_68*tmp_87;
      real_t tmp_91 = -tmp_78 - tmp_79 - tmp_82 - tmp_83 - tmp_84 - tmp_85 - tmp_88 - tmp_89 - tmp_90 + 1;
      real_t tmp_92 = tmp_1*tmp_19;
      real_t tmp_93 = tmp_14*tmp_19;
      real_t tmp_94 = tmp_10*tmp_19;
      real_t tmp_95 = 0.5*p_affine_13_0*(tmp_33*tmp_92 + tmp_42*tmp_93 + tmp_45*tmp_94) + 0.5*p_affine_13_1*(tmp_26*tmp_92 + tmp_41*tmp_93 + tmp_44*tmp_94) + 0.5*p_affine_13_2*(tmp_40*tmp_93 + tmp_43*tmp_94 + tmp_8*tmp_92);
      real_t tmp_96 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_97 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_98 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_99 = (std::abs(tmp_22*tmp_96 - tmp_29*tmp_98)*std::abs(tmp_22*tmp_96 - tmp_29*tmp_98)) + (std::abs(tmp_22*tmp_97 - tmp_36*tmp_98)*std::abs(tmp_22*tmp_97 - tmp_36*tmp_98)) + (std::abs(tmp_29*tmp_97 - tmp_36*tmp_96)*std::abs(tmp_29*tmp_97 - tmp_36*tmp_96));
      real_t tmp_100 = 5.0*std::pow(tmp_99, -0.25);
      real_t tmp_101 = tmp_100*tmp_46;
      real_t tmp_102 = 1.0*std::pow(tmp_99, 1.0/2.0);
      real_t tmp_103 = 0.0068572537431980923*tmp_102;
      real_t tmp_104 = 0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22;
      real_t tmp_105 = tmp_19*(tmp_104 + tmp_24);
      real_t tmp_106 = 0.19601935860219369*tmp_28 + 0.60796128279561268*tmp_29;
      real_t tmp_107 = tmp_19*(tmp_106 + tmp_31);
      real_t tmp_108 = 0.19601935860219369*tmp_35 + 0.60796128279561268*tmp_36;
      real_t tmp_109 = tmp_19*(tmp_108 + tmp_38);
      real_t tmp_110 = tmp_1*(tmp_105*tmp_8 + tmp_107*tmp_26 + tmp_109*tmp_33 - 1.0/4.0) + tmp_10*(tmp_105*tmp_43 + tmp_107*tmp_44 + tmp_109*tmp_45 - 1.0/4.0) + tmp_14*(tmp_105*tmp_40 + tmp_107*tmp_41 + tmp_109*tmp_42 - 1.0/4.0);
      real_t tmp_111 = tmp_104 + tmp_76;
      real_t tmp_112 = tmp_111*tmp_72;
      real_t tmp_113 = tmp_111*tmp_73;
      real_t tmp_114 = tmp_106 + tmp_80;
      real_t tmp_115 = tmp_114*tmp_69;
      real_t tmp_116 = tmp_114*tmp_70;
      real_t tmp_117 = tmp_111*tmp_74;
      real_t tmp_118 = tmp_114*tmp_71;
      real_t tmp_119 = tmp_108 + tmp_86;
      real_t tmp_120 = tmp_119*tmp_66;
      real_t tmp_121 = tmp_119*tmp_67;
      real_t tmp_122 = tmp_119*tmp_68;
      real_t tmp_123 = -tmp_112 - tmp_113 - tmp_115 - tmp_116 - tmp_117 - tmp_118 - tmp_120 - tmp_121 - tmp_122 + 1;
      real_t tmp_124 = tmp_100*tmp_110;
      real_t tmp_125 = 0.037198804536718075*tmp_102;
      real_t tmp_126 = 0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_127 = tmp_19*(tmp_126 + tmp_24);
      real_t tmp_128 = 0.37605877282253791*tmp_28 + 0.039308471900058539*tmp_29;
      real_t tmp_129 = tmp_19*(tmp_128 + tmp_31);
      real_t tmp_130 = 0.37605877282253791*tmp_35 + 0.039308471900058539*tmp_36;
      real_t tmp_131 = tmp_19*(tmp_130 + tmp_38);
      real_t tmp_132 = tmp_1*(tmp_127*tmp_8 + tmp_129*tmp_26 + tmp_131*tmp_33 - 1.0/4.0) + tmp_10*(tmp_127*tmp_43 + tmp_129*tmp_44 + tmp_131*tmp_45 - 1.0/4.0) + tmp_14*(tmp_127*tmp_40 + tmp_129*tmp_41 + tmp_131*tmp_42 - 1.0/4.0);
      real_t tmp_133 = tmp_126 + tmp_76;
      real_t tmp_134 = tmp_133*tmp_72;
      real_t tmp_135 = tmp_133*tmp_73;
      real_t tmp_136 = tmp_128 + tmp_80;
      real_t tmp_137 = tmp_136*tmp_69;
      real_t tmp_138 = tmp_136*tmp_70;
      real_t tmp_139 = tmp_133*tmp_74;
      real_t tmp_140 = tmp_136*tmp_71;
      real_t tmp_141 = tmp_130 + tmp_86;
      real_t tmp_142 = tmp_141*tmp_66;
      real_t tmp_143 = tmp_141*tmp_67;
      real_t tmp_144 = tmp_141*tmp_68;
      real_t tmp_145 = -tmp_134 - tmp_135 - tmp_137 - tmp_138 - tmp_139 - tmp_140 - tmp_142 - tmp_143 - tmp_144 + 1;
      real_t tmp_146 = tmp_100*tmp_132;
      real_t tmp_147 = 0.020848748529055869*tmp_102;
      real_t tmp_148 = 0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_149 = tmp_19*(tmp_148 + tmp_24);
      real_t tmp_150 = 0.78764240869137092*tmp_28 + 0.1711304259088916*tmp_29;
      real_t tmp_151 = tmp_19*(tmp_150 + tmp_31);
      real_t tmp_152 = 0.78764240869137092*tmp_35 + 0.1711304259088916*tmp_36;
      real_t tmp_153 = tmp_19*(tmp_152 + tmp_38);
      real_t tmp_154 = tmp_1*(tmp_149*tmp_8 + tmp_151*tmp_26 + tmp_153*tmp_33 - 1.0/4.0) + tmp_10*(tmp_149*tmp_43 + tmp_151*tmp_44 + tmp_153*tmp_45 - 1.0/4.0) + tmp_14*(tmp_149*tmp_40 + tmp_151*tmp_41 + tmp_153*tmp_42 - 1.0/4.0);
      real_t tmp_155 = tmp_148 + tmp_76;
      real_t tmp_156 = tmp_155*tmp_72;
      real_t tmp_157 = tmp_155*tmp_73;
      real_t tmp_158 = tmp_150 + tmp_80;
      real_t tmp_159 = tmp_158*tmp_69;
      real_t tmp_160 = tmp_158*tmp_70;
      real_t tmp_161 = tmp_155*tmp_74;
      real_t tmp_162 = tmp_158*tmp_71;
      real_t tmp_163 = tmp_152 + tmp_86;
      real_t tmp_164 = tmp_163*tmp_66;
      real_t tmp_165 = tmp_163*tmp_67;
      real_t tmp_166 = tmp_163*tmp_68;
      real_t tmp_167 = -tmp_156 - tmp_157 - tmp_159 - tmp_160 - tmp_161 - tmp_162 - tmp_164 - tmp_165 - tmp_166 + 1;
      real_t tmp_168 = tmp_100*tmp_154;
      real_t tmp_169 = 0.019202922745021479*tmp_102;
      real_t tmp_170 = 0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_171 = tmp_19*(tmp_170 + tmp_24);
      real_t tmp_172 = 0.58463275527740355*tmp_28 + 0.37605877282253791*tmp_29;
      real_t tmp_173 = tmp_19*(tmp_172 + tmp_31);
      real_t tmp_174 = 0.58463275527740355*tmp_35 + 0.37605877282253791*tmp_36;
      real_t tmp_175 = tmp_19*(tmp_174 + tmp_38);
      real_t tmp_176 = tmp_1*(tmp_171*tmp_8 + tmp_173*tmp_26 + tmp_175*tmp_33 - 1.0/4.0) + tmp_10*(tmp_171*tmp_43 + tmp_173*tmp_44 + tmp_175*tmp_45 - 1.0/4.0) + tmp_14*(tmp_171*tmp_40 + tmp_173*tmp_41 + tmp_175*tmp_42 - 1.0/4.0);
      real_t tmp_177 = tmp_170 + tmp_76;
      real_t tmp_178 = tmp_177*tmp_72;
      real_t tmp_179 = tmp_177*tmp_73;
      real_t tmp_180 = tmp_172 + tmp_80;
      real_t tmp_181 = tmp_180*tmp_69;
      real_t tmp_182 = tmp_180*tmp_70;
      real_t tmp_183 = tmp_177*tmp_74;
      real_t tmp_184 = tmp_180*tmp_71;
      real_t tmp_185 = tmp_174 + tmp_86;
      real_t tmp_186 = tmp_185*tmp_66;
      real_t tmp_187 = tmp_185*tmp_67;
      real_t tmp_188 = tmp_185*tmp_68;
      real_t tmp_189 = -tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 - tmp_184 - tmp_186 - tmp_187 - tmp_188 + 1;
      real_t tmp_190 = tmp_100*tmp_176;
      real_t tmp_191 = 0.020848748529055869*tmp_102;
      real_t tmp_192 = 0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_193 = tmp_19*(tmp_192 + tmp_24);
      real_t tmp_194 = 0.041227165399737475*tmp_28 + 0.78764240869137092*tmp_29;
      real_t tmp_195 = tmp_19*(tmp_194 + tmp_31);
      real_t tmp_196 = 0.041227165399737475*tmp_35 + 0.78764240869137092*tmp_36;
      real_t tmp_197 = tmp_19*(tmp_196 + tmp_38);
      real_t tmp_198 = tmp_1*(tmp_193*tmp_8 + tmp_195*tmp_26 + tmp_197*tmp_33 - 1.0/4.0) + tmp_10*(tmp_193*tmp_43 + tmp_195*tmp_44 + tmp_197*tmp_45 - 1.0/4.0) + tmp_14*(tmp_193*tmp_40 + tmp_195*tmp_41 + tmp_197*tmp_42 - 1.0/4.0);
      real_t tmp_199 = tmp_192 + tmp_76;
      real_t tmp_200 = tmp_199*tmp_72;
      real_t tmp_201 = tmp_199*tmp_73;
      real_t tmp_202 = tmp_194 + tmp_80;
      real_t tmp_203 = tmp_202*tmp_69;
      real_t tmp_204 = tmp_202*tmp_70;
      real_t tmp_205 = tmp_199*tmp_74;
      real_t tmp_206 = tmp_202*tmp_71;
      real_t tmp_207 = tmp_196 + tmp_86;
      real_t tmp_208 = tmp_207*tmp_66;
      real_t tmp_209 = tmp_207*tmp_67;
      real_t tmp_210 = tmp_207*tmp_68;
      real_t tmp_211 = -tmp_200 - tmp_201 - tmp_203 - tmp_204 - tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 + 1;
      real_t tmp_212 = tmp_100*tmp_198;
      real_t tmp_213 = 0.019202922745021479*tmp_102;
      real_t tmp_214 = 0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_215 = tmp_19*(tmp_214 + tmp_24);
      real_t tmp_216 = 0.039308471900058539*tmp_28 + 0.58463275527740355*tmp_29;
      real_t tmp_217 = tmp_19*(tmp_216 + tmp_31);
      real_t tmp_218 = 0.039308471900058539*tmp_35 + 0.58463275527740355*tmp_36;
      real_t tmp_219 = tmp_19*(tmp_218 + tmp_38);
      real_t tmp_220 = tmp_1*(tmp_215*tmp_8 + tmp_217*tmp_26 + tmp_219*tmp_33 - 1.0/4.0) + tmp_10*(tmp_215*tmp_43 + tmp_217*tmp_44 + tmp_219*tmp_45 - 1.0/4.0) + tmp_14*(tmp_215*tmp_40 + tmp_217*tmp_41 + tmp_219*tmp_42 - 1.0/4.0);
      real_t tmp_221 = tmp_214 + tmp_76;
      real_t tmp_222 = tmp_221*tmp_72;
      real_t tmp_223 = tmp_221*tmp_73;
      real_t tmp_224 = tmp_216 + tmp_80;
      real_t tmp_225 = tmp_224*tmp_69;
      real_t tmp_226 = tmp_224*tmp_70;
      real_t tmp_227 = tmp_221*tmp_74;
      real_t tmp_228 = tmp_224*tmp_71;
      real_t tmp_229 = tmp_218 + tmp_86;
      real_t tmp_230 = tmp_229*tmp_66;
      real_t tmp_231 = tmp_229*tmp_67;
      real_t tmp_232 = tmp_229*tmp_68;
      real_t tmp_233 = -tmp_222 - tmp_223 - tmp_225 - tmp_226 - tmp_227 - tmp_228 - tmp_230 - tmp_231 - tmp_232 + 1;
      real_t tmp_234 = tmp_100*tmp_220;
      real_t tmp_235 = 0.020848748529055869*tmp_102;
      real_t tmp_236 = 0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_237 = tmp_19*(tmp_236 + tmp_24);
      real_t tmp_238 = 0.78764240869137092*tmp_28 + 0.041227165399737475*tmp_29;
      real_t tmp_239 = tmp_19*(tmp_238 + tmp_31);
      real_t tmp_240 = 0.78764240869137092*tmp_35 + 0.041227165399737475*tmp_36;
      real_t tmp_241 = tmp_19*(tmp_240 + tmp_38);
      real_t tmp_242 = tmp_1*(tmp_237*tmp_8 + tmp_239*tmp_26 + tmp_241*tmp_33 - 1.0/4.0) + tmp_10*(tmp_237*tmp_43 + tmp_239*tmp_44 + tmp_241*tmp_45 - 1.0/4.0) + tmp_14*(tmp_237*tmp_40 + tmp_239*tmp_41 + tmp_241*tmp_42 - 1.0/4.0);
      real_t tmp_243 = tmp_236 + tmp_76;
      real_t tmp_244 = tmp_243*tmp_72;
      real_t tmp_245 = tmp_243*tmp_73;
      real_t tmp_246 = tmp_238 + tmp_80;
      real_t tmp_247 = tmp_246*tmp_69;
      real_t tmp_248 = tmp_246*tmp_70;
      real_t tmp_249 = tmp_243*tmp_74;
      real_t tmp_250 = tmp_246*tmp_71;
      real_t tmp_251 = tmp_240 + tmp_86;
      real_t tmp_252 = tmp_251*tmp_66;
      real_t tmp_253 = tmp_251*tmp_67;
      real_t tmp_254 = tmp_251*tmp_68;
      real_t tmp_255 = -tmp_244 - tmp_245 - tmp_247 - tmp_248 - tmp_249 - tmp_250 - tmp_252 - tmp_253 - tmp_254 + 1;
      real_t tmp_256 = tmp_100*tmp_242;
      real_t tmp_257 = 0.019202922745021479*tmp_102;
      real_t tmp_258 = 0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22;
      real_t tmp_259 = tmp_19*(tmp_24 + tmp_258);
      real_t tmp_260 = 0.58463275527740355*tmp_28 + 0.039308471900058539*tmp_29;
      real_t tmp_261 = tmp_19*(tmp_260 + tmp_31);
      real_t tmp_262 = 0.58463275527740355*tmp_35 + 0.039308471900058539*tmp_36;
      real_t tmp_263 = tmp_19*(tmp_262 + tmp_38);
      real_t tmp_264 = tmp_1*(tmp_259*tmp_8 + tmp_26*tmp_261 + tmp_263*tmp_33 - 1.0/4.0) + tmp_10*(tmp_259*tmp_43 + tmp_261*tmp_44 + tmp_263*tmp_45 - 1.0/4.0) + tmp_14*(tmp_259*tmp_40 + tmp_261*tmp_41 + tmp_263*tmp_42 - 1.0/4.0);
      real_t tmp_265 = tmp_258 + tmp_76;
      real_t tmp_266 = tmp_265*tmp_72;
      real_t tmp_267 = tmp_265*tmp_73;
      real_t tmp_268 = tmp_260 + tmp_80;
      real_t tmp_269 = tmp_268*tmp_69;
      real_t tmp_270 = tmp_268*tmp_70;
      real_t tmp_271 = tmp_265*tmp_74;
      real_t tmp_272 = tmp_268*tmp_71;
      real_t tmp_273 = tmp_262 + tmp_86;
      real_t tmp_274 = tmp_273*tmp_66;
      real_t tmp_275 = tmp_273*tmp_67;
      real_t tmp_276 = tmp_273*tmp_68;
      real_t tmp_277 = -tmp_266 - tmp_267 - tmp_269 - tmp_270 - tmp_271 - tmp_272 - tmp_274 - tmp_275 - tmp_276 + 1;
      real_t tmp_278 = tmp_100*tmp_264;
      real_t tmp_279 = 0.020848748529055869*tmp_102;
      real_t tmp_280 = 0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22;
      real_t tmp_281 = tmp_19*(tmp_24 + tmp_280);
      real_t tmp_282 = 0.1711304259088916*tmp_28 + 0.78764240869137092*tmp_29;
      real_t tmp_283 = tmp_19*(tmp_282 + tmp_31);
      real_t tmp_284 = 0.1711304259088916*tmp_35 + 0.78764240869137092*tmp_36;
      real_t tmp_285 = tmp_19*(tmp_284 + tmp_38);
      real_t tmp_286 = tmp_1*(tmp_26*tmp_283 + tmp_281*tmp_8 + tmp_285*tmp_33 - 1.0/4.0) + tmp_10*(tmp_281*tmp_43 + tmp_283*tmp_44 + tmp_285*tmp_45 - 1.0/4.0) + tmp_14*(tmp_281*tmp_40 + tmp_283*tmp_41 + tmp_285*tmp_42 - 1.0/4.0);
      real_t tmp_287 = tmp_280 + tmp_76;
      real_t tmp_288 = tmp_287*tmp_72;
      real_t tmp_289 = tmp_287*tmp_73;
      real_t tmp_290 = tmp_282 + tmp_80;
      real_t tmp_291 = tmp_290*tmp_69;
      real_t tmp_292 = tmp_290*tmp_70;
      real_t tmp_293 = tmp_287*tmp_74;
      real_t tmp_294 = tmp_290*tmp_71;
      real_t tmp_295 = tmp_284 + tmp_86;
      real_t tmp_296 = tmp_295*tmp_66;
      real_t tmp_297 = tmp_295*tmp_67;
      real_t tmp_298 = tmp_295*tmp_68;
      real_t tmp_299 = -tmp_288 - tmp_289 - tmp_291 - tmp_292 - tmp_293 - tmp_294 - tmp_296 - tmp_297 - tmp_298 + 1;
      real_t tmp_300 = tmp_100*tmp_286;
      real_t tmp_301 = 0.019202922745021479*tmp_102;
      real_t tmp_302 = 0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22;
      real_t tmp_303 = tmp_19*(tmp_24 + tmp_302);
      real_t tmp_304 = 0.37605877282253791*tmp_28 + 0.58463275527740355*tmp_29;
      real_t tmp_305 = tmp_19*(tmp_304 + tmp_31);
      real_t tmp_306 = 0.37605877282253791*tmp_35 + 0.58463275527740355*tmp_36;
      real_t tmp_307 = tmp_19*(tmp_306 + tmp_38);
      real_t tmp_308 = tmp_1*(tmp_26*tmp_305 + tmp_303*tmp_8 + tmp_307*tmp_33 - 1.0/4.0) + tmp_10*(tmp_303*tmp_43 + tmp_305*tmp_44 + tmp_307*tmp_45 - 1.0/4.0) + tmp_14*(tmp_303*tmp_40 + tmp_305*tmp_41 + tmp_307*tmp_42 - 1.0/4.0);
      real_t tmp_309 = tmp_302 + tmp_76;
      real_t tmp_310 = tmp_309*tmp_72;
      real_t tmp_311 = tmp_309*tmp_73;
      real_t tmp_312 = tmp_304 + tmp_80;
      real_t tmp_313 = tmp_312*tmp_69;
      real_t tmp_314 = tmp_312*tmp_70;
      real_t tmp_315 = tmp_309*tmp_74;
      real_t tmp_316 = tmp_312*tmp_71;
      real_t tmp_317 = tmp_306 + tmp_86;
      real_t tmp_318 = tmp_317*tmp_66;
      real_t tmp_319 = tmp_317*tmp_67;
      real_t tmp_320 = tmp_317*tmp_68;
      real_t tmp_321 = -tmp_310 - tmp_311 - tmp_313 - tmp_314 - tmp_315 - tmp_316 - tmp_318 - tmp_319 - tmp_320 + 1;
      real_t tmp_322 = tmp_100*tmp_308;
      real_t tmp_323 = 0.020848748529055869*tmp_102;
      real_t tmp_324 = 0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22;
      real_t tmp_325 = tmp_19*(tmp_24 + tmp_324);
      real_t tmp_326 = 0.041227165399737475*tmp_28 + 0.1711304259088916*tmp_29;
      real_t tmp_327 = tmp_19*(tmp_31 + tmp_326);
      real_t tmp_328 = 0.041227165399737475*tmp_35 + 0.1711304259088916*tmp_36;
      real_t tmp_329 = tmp_19*(tmp_328 + tmp_38);
      real_t tmp_330 = tmp_1*(tmp_26*tmp_327 + tmp_325*tmp_8 + tmp_329*tmp_33 - 1.0/4.0) + tmp_10*(tmp_325*tmp_43 + tmp_327*tmp_44 + tmp_329*tmp_45 - 1.0/4.0) + tmp_14*(tmp_325*tmp_40 + tmp_327*tmp_41 + tmp_329*tmp_42 - 1.0/4.0);
      real_t tmp_331 = tmp_324 + tmp_76;
      real_t tmp_332 = tmp_331*tmp_72;
      real_t tmp_333 = tmp_331*tmp_73;
      real_t tmp_334 = tmp_326 + tmp_80;
      real_t tmp_335 = tmp_334*tmp_69;
      real_t tmp_336 = tmp_334*tmp_70;
      real_t tmp_337 = tmp_331*tmp_74;
      real_t tmp_338 = tmp_334*tmp_71;
      real_t tmp_339 = tmp_328 + tmp_86;
      real_t tmp_340 = tmp_339*tmp_66;
      real_t tmp_341 = tmp_339*tmp_67;
      real_t tmp_342 = tmp_339*tmp_68;
      real_t tmp_343 = -tmp_332 - tmp_333 - tmp_335 - tmp_336 - tmp_337 - tmp_338 - tmp_340 - tmp_341 - tmp_342 + 1;
      real_t tmp_344 = tmp_100*tmp_330;
      real_t tmp_345 = 0.019202922745021479*tmp_102;
      real_t tmp_346 = 0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22;
      real_t tmp_347 = tmp_19*(tmp_24 + tmp_346);
      real_t tmp_348 = 0.40446199974765351*tmp_28 + 0.19107600050469298*tmp_29;
      real_t tmp_349 = tmp_19*(tmp_31 + tmp_348);
      real_t tmp_350 = 0.40446199974765351*tmp_35 + 0.19107600050469298*tmp_36;
      real_t tmp_351 = tmp_19*(tmp_350 + tmp_38);
      real_t tmp_352 = tmp_1*(tmp_26*tmp_349 + tmp_33*tmp_351 + tmp_347*tmp_8 - 1.0/4.0) + tmp_10*(tmp_347*tmp_43 + tmp_349*tmp_44 + tmp_351*tmp_45 - 1.0/4.0) + tmp_14*(tmp_347*tmp_40 + tmp_349*tmp_41 + tmp_351*tmp_42 - 1.0/4.0);
      real_t tmp_353 = tmp_346 + tmp_76;
      real_t tmp_354 = tmp_353*tmp_72;
      real_t tmp_355 = tmp_353*tmp_73;
      real_t tmp_356 = tmp_348 + tmp_80;
      real_t tmp_357 = tmp_356*tmp_69;
      real_t tmp_358 = tmp_356*tmp_70;
      real_t tmp_359 = tmp_353*tmp_74;
      real_t tmp_360 = tmp_356*tmp_71;
      real_t tmp_361 = tmp_350 + tmp_86;
      real_t tmp_362 = tmp_361*tmp_66;
      real_t tmp_363 = tmp_361*tmp_67;
      real_t tmp_364 = tmp_361*tmp_68;
      real_t tmp_365 = -tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 - tmp_360 - tmp_362 - tmp_363 - tmp_364 + 1;
      real_t tmp_366 = tmp_100*tmp_352;
      real_t tmp_367 = 0.042507265838595799*tmp_102;
      real_t tmp_368 = 0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22;
      real_t tmp_369 = tmp_19*(tmp_24 + tmp_368);
      real_t tmp_370 = 0.039308471900058539*tmp_28 + 0.37605877282253791*tmp_29;
      real_t tmp_371 = tmp_19*(tmp_31 + tmp_370);
      real_t tmp_372 = 0.039308471900058539*tmp_35 + 0.37605877282253791*tmp_36;
      real_t tmp_373 = tmp_19*(tmp_372 + tmp_38);
      real_t tmp_374 = tmp_1*(tmp_26*tmp_371 + tmp_33*tmp_373 + tmp_369*tmp_8 - 1.0/4.0) + tmp_10*(tmp_369*tmp_43 + tmp_371*tmp_44 + tmp_373*tmp_45 - 1.0/4.0) + tmp_14*(tmp_369*tmp_40 + tmp_371*tmp_41 + tmp_373*tmp_42 - 1.0/4.0);
      real_t tmp_375 = tmp_368 + tmp_76;
      real_t tmp_376 = tmp_375*tmp_72;
      real_t tmp_377 = tmp_375*tmp_73;
      real_t tmp_378 = tmp_370 + tmp_80;
      real_t tmp_379 = tmp_378*tmp_69;
      real_t tmp_380 = tmp_378*tmp_70;
      real_t tmp_381 = tmp_375*tmp_74;
      real_t tmp_382 = tmp_378*tmp_71;
      real_t tmp_383 = tmp_372 + tmp_86;
      real_t tmp_384 = tmp_383*tmp_66;
      real_t tmp_385 = tmp_383*tmp_67;
      real_t tmp_386 = tmp_383*tmp_68;
      real_t tmp_387 = -tmp_376 - tmp_377 - tmp_379 - tmp_380 - tmp_381 - tmp_382 - tmp_384 - tmp_385 - tmp_386 + 1;
      real_t tmp_388 = tmp_100*tmp_374;
      real_t tmp_389 = 0.020848748529055869*tmp_102;
      real_t tmp_390 = 0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_391 = tmp_19*(tmp_24 + tmp_390);
      real_t tmp_392 = 0.93718850182767688*tmp_28 + 0.031405749086161582*tmp_29;
      real_t tmp_393 = tmp_19*(tmp_31 + tmp_392);
      real_t tmp_394 = 0.93718850182767688*tmp_35 + 0.031405749086161582*tmp_36;
      real_t tmp_395 = tmp_19*(tmp_38 + tmp_394);
      real_t tmp_396 = tmp_1*(tmp_26*tmp_393 + tmp_33*tmp_395 + tmp_391*tmp_8 - 1.0/4.0) + tmp_10*(tmp_391*tmp_43 + tmp_393*tmp_44 + tmp_395*tmp_45 - 1.0/4.0) + tmp_14*(tmp_391*tmp_40 + tmp_393*tmp_41 + tmp_395*tmp_42 - 1.0/4.0);
      real_t tmp_397 = tmp_390 + tmp_76;
      real_t tmp_398 = tmp_397*tmp_72;
      real_t tmp_399 = tmp_397*tmp_73;
      real_t tmp_400 = tmp_392 + tmp_80;
      real_t tmp_401 = tmp_400*tmp_69;
      real_t tmp_402 = tmp_400*tmp_70;
      real_t tmp_403 = tmp_397*tmp_74;
      real_t tmp_404 = tmp_400*tmp_71;
      real_t tmp_405 = tmp_394 + tmp_86;
      real_t tmp_406 = tmp_405*tmp_66;
      real_t tmp_407 = tmp_405*tmp_67;
      real_t tmp_408 = tmp_405*tmp_68;
      real_t tmp_409 = -tmp_398 - tmp_399 - tmp_401 - tmp_402 - tmp_403 - tmp_404 - tmp_406 - tmp_407 - tmp_408 + 1;
      real_t tmp_410 = tmp_100*tmp_396;
      real_t tmp_411 = 0.0068572537431980923*tmp_102;
      real_t tmp_412 = 0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_413 = tmp_19*(tmp_24 + tmp_412);
      real_t tmp_414 = 0.60796128279561268*tmp_28 + 0.19601935860219369*tmp_29;
      real_t tmp_415 = tmp_19*(tmp_31 + tmp_414);
      real_t tmp_416 = 0.60796128279561268*tmp_35 + 0.19601935860219369*tmp_36;
      real_t tmp_417 = tmp_19*(tmp_38 + tmp_416);
      real_t tmp_418 = tmp_1*(tmp_26*tmp_415 + tmp_33*tmp_417 + tmp_413*tmp_8 - 1.0/4.0) + tmp_10*(tmp_413*tmp_43 + tmp_415*tmp_44 + tmp_417*tmp_45 - 1.0/4.0) + tmp_14*(tmp_40*tmp_413 + tmp_41*tmp_415 + tmp_417*tmp_42 - 1.0/4.0);
      real_t tmp_419 = tmp_412 + tmp_76;
      real_t tmp_420 = tmp_419*tmp_72;
      real_t tmp_421 = tmp_419*tmp_73;
      real_t tmp_422 = tmp_414 + tmp_80;
      real_t tmp_423 = tmp_422*tmp_69;
      real_t tmp_424 = tmp_422*tmp_70;
      real_t tmp_425 = tmp_419*tmp_74;
      real_t tmp_426 = tmp_422*tmp_71;
      real_t tmp_427 = tmp_416 + tmp_86;
      real_t tmp_428 = tmp_427*tmp_66;
      real_t tmp_429 = tmp_427*tmp_67;
      real_t tmp_430 = tmp_427*tmp_68;
      real_t tmp_431 = -tmp_420 - tmp_421 - tmp_423 - tmp_424 - tmp_425 - tmp_426 - tmp_428 - tmp_429 - tmp_430 + 1;
      real_t tmp_432 = tmp_100*tmp_418;
      real_t tmp_433 = 0.037198804536718075*tmp_102;
      real_t tmp_434 = 0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_435 = tmp_19*(tmp_24 + tmp_434);
      real_t tmp_436 = 0.19107600050469298*tmp_28 + 0.40446199974765351*tmp_29;
      real_t tmp_437 = tmp_19*(tmp_31 + tmp_436);
      real_t tmp_438 = 0.19107600050469298*tmp_35 + 0.40446199974765351*tmp_36;
      real_t tmp_439 = tmp_19*(tmp_38 + tmp_438);
      real_t tmp_440 = tmp_1*(tmp_26*tmp_437 + tmp_33*tmp_439 + tmp_435*tmp_8 - 1.0/4.0) + tmp_10*(tmp_43*tmp_435 + tmp_437*tmp_44 + tmp_439*tmp_45 - 1.0/4.0) + tmp_14*(tmp_40*tmp_435 + tmp_41*tmp_437 + tmp_42*tmp_439 - 1.0/4.0);
      real_t tmp_441 = tmp_434 + tmp_76;
      real_t tmp_442 = tmp_441*tmp_72;
      real_t tmp_443 = tmp_441*tmp_73;
      real_t tmp_444 = tmp_436 + tmp_80;
      real_t tmp_445 = tmp_444*tmp_69;
      real_t tmp_446 = tmp_444*tmp_70;
      real_t tmp_447 = tmp_441*tmp_74;
      real_t tmp_448 = tmp_444*tmp_71;
      real_t tmp_449 = tmp_438 + tmp_86;
      real_t tmp_450 = tmp_449*tmp_66;
      real_t tmp_451 = tmp_449*tmp_67;
      real_t tmp_452 = tmp_449*tmp_68;
      real_t tmp_453 = -tmp_442 - tmp_443 - tmp_445 - tmp_446 - tmp_447 - tmp_448 - tmp_450 - tmp_451 - tmp_452 + 1;
      real_t tmp_454 = tmp_100*tmp_440;
      real_t tmp_455 = 0.042507265838595799*tmp_102;
      real_t tmp_456 = 0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22;
      real_t tmp_457 = tmp_19*(tmp_24 + tmp_456);
      real_t tmp_458 = 0.031405749086161582*tmp_28 + 0.031405749086161582*tmp_29;
      real_t tmp_459 = tmp_19*(tmp_31 + tmp_458);
      real_t tmp_460 = 0.031405749086161582*tmp_35 + 0.031405749086161582*tmp_36;
      real_t tmp_461 = tmp_19*(tmp_38 + tmp_460);
      real_t tmp_462 = tmp_1*(tmp_26*tmp_459 + tmp_33*tmp_461 + tmp_457*tmp_8 - 1.0/4.0) + tmp_10*(tmp_43*tmp_457 + tmp_44*tmp_459 + tmp_45*tmp_461 - 1.0/4.0) + tmp_14*(tmp_40*tmp_457 + tmp_41*tmp_459 + tmp_42*tmp_461 - 1.0/4.0);
      real_t tmp_463 = tmp_456 + tmp_76;
      real_t tmp_464 = tmp_463*tmp_72;
      real_t tmp_465 = tmp_463*tmp_73;
      real_t tmp_466 = tmp_458 + tmp_80;
      real_t tmp_467 = tmp_466*tmp_69;
      real_t tmp_468 = tmp_466*tmp_70;
      real_t tmp_469 = tmp_463*tmp_74;
      real_t tmp_470 = tmp_466*tmp_71;
      real_t tmp_471 = tmp_460 + tmp_86;
      real_t tmp_472 = tmp_471*tmp_66;
      real_t tmp_473 = tmp_471*tmp_67;
      real_t tmp_474 = tmp_471*tmp_68;
      real_t tmp_475 = -tmp_464 - tmp_465 - tmp_467 - tmp_468 - tmp_469 - tmp_470 - tmp_472 - tmp_473 - tmp_474 + 1;
      real_t tmp_476 = tmp_100*tmp_462;
      real_t tmp_477 = 0.0068572537431980923*tmp_102;
      real_t tmp_478 = 0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22;
      real_t tmp_479 = tmp_19*(tmp_24 + tmp_478);
      real_t tmp_480 = 0.19601935860219369*tmp_28 + 0.19601935860219369*tmp_29;
      real_t tmp_481 = tmp_19*(tmp_31 + tmp_480);
      real_t tmp_482 = 0.19601935860219369*tmp_35 + 0.19601935860219369*tmp_36;
      real_t tmp_483 = tmp_19*(tmp_38 + tmp_482);
      real_t tmp_484 = tmp_1*(tmp_26*tmp_481 + tmp_33*tmp_483 + tmp_479*tmp_8 - 1.0/4.0) + tmp_10*(tmp_43*tmp_479 + tmp_44*tmp_481 + tmp_45*tmp_483 - 1.0/4.0) + tmp_14*(tmp_40*tmp_479 + tmp_41*tmp_481 + tmp_42*tmp_483 - 1.0/4.0);
      real_t tmp_485 = tmp_478 + tmp_76;
      real_t tmp_486 = tmp_485*tmp_72;
      real_t tmp_487 = tmp_485*tmp_73;
      real_t tmp_488 = tmp_480 + tmp_80;
      real_t tmp_489 = tmp_488*tmp_69;
      real_t tmp_490 = tmp_488*tmp_70;
      real_t tmp_491 = tmp_485*tmp_74;
      real_t tmp_492 = tmp_488*tmp_71;
      real_t tmp_493 = tmp_482 + tmp_86;
      real_t tmp_494 = tmp_493*tmp_66;
      real_t tmp_495 = tmp_493*tmp_67;
      real_t tmp_496 = tmp_493*tmp_68;
      real_t tmp_497 = -tmp_486 - tmp_487 - tmp_489 - tmp_490 - tmp_491 - tmp_492 - tmp_494 - tmp_495 - tmp_496 + 1;
      real_t tmp_498 = tmp_100*tmp_484;
      real_t tmp_499 = 0.037198804536718075*tmp_102;
      real_t tmp_500 = 0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22;
      real_t tmp_501 = tmp_19*(tmp_24 + tmp_500);
      real_t tmp_502 = 0.40446199974765351*tmp_28 + 0.40446199974765351*tmp_29;
      real_t tmp_503 = tmp_19*(tmp_31 + tmp_502);
      real_t tmp_504 = 0.40446199974765351*tmp_35 + 0.40446199974765351*tmp_36;
      real_t tmp_505 = tmp_19*(tmp_38 + tmp_504);
      real_t tmp_506 = tmp_1*(tmp_26*tmp_503 + tmp_33*tmp_505 + tmp_501*tmp_8 - 1.0/4.0) + tmp_10*(tmp_43*tmp_501 + tmp_44*tmp_503 + tmp_45*tmp_505 - 1.0/4.0) + tmp_14*(tmp_40*tmp_501 + tmp_41*tmp_503 + tmp_42*tmp_505 - 1.0/4.0);
      real_t tmp_507 = tmp_500 + tmp_76;
      real_t tmp_508 = tmp_507*tmp_72;
      real_t tmp_509 = tmp_507*tmp_73;
      real_t tmp_510 = tmp_502 + tmp_80;
      real_t tmp_511 = tmp_510*tmp_69;
      real_t tmp_512 = tmp_510*tmp_70;
      real_t tmp_513 = tmp_507*tmp_74;
      real_t tmp_514 = tmp_510*tmp_71;
      real_t tmp_515 = tmp_504 + tmp_86;
      real_t tmp_516 = tmp_515*tmp_66;
      real_t tmp_517 = tmp_515*tmp_67;
      real_t tmp_518 = tmp_515*tmp_68;
      real_t tmp_519 = -tmp_508 - tmp_509 - tmp_511 - tmp_512 - tmp_513 - tmp_514 - tmp_516 - tmp_517 - tmp_518 + 1;
      real_t tmp_520 = tmp_100*tmp_506;
      real_t tmp_521 = 0.042507265838595799*tmp_102;
      real_t tmp_522 = 0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22;
      real_t tmp_523 = tmp_19*(tmp_24 + tmp_522);
      real_t tmp_524 = 0.1711304259088916*tmp_28 + 0.041227165399737475*tmp_29;
      real_t tmp_525 = tmp_19*(tmp_31 + tmp_524);
      real_t tmp_526 = 0.1711304259088916*tmp_35 + 0.041227165399737475*tmp_36;
      real_t tmp_527 = tmp_19*(tmp_38 + tmp_526);
      real_t tmp_528 = tmp_1*(tmp_26*tmp_525 + tmp_33*tmp_527 + tmp_523*tmp_8 - 1.0/4.0) + tmp_10*(tmp_43*tmp_523 + tmp_44*tmp_525 + tmp_45*tmp_527 - 1.0/4.0) + tmp_14*(tmp_40*tmp_523 + tmp_41*tmp_525 + tmp_42*tmp_527 - 1.0/4.0);
      real_t tmp_529 = tmp_522 + tmp_76;
      real_t tmp_530 = tmp_529*tmp_72;
      real_t tmp_531 = tmp_529*tmp_73;
      real_t tmp_532 = tmp_524 + tmp_80;
      real_t tmp_533 = tmp_532*tmp_69;
      real_t tmp_534 = tmp_532*tmp_70;
      real_t tmp_535 = tmp_529*tmp_74;
      real_t tmp_536 = tmp_532*tmp_71;
      real_t tmp_537 = tmp_526 + tmp_86;
      real_t tmp_538 = tmp_537*tmp_66;
      real_t tmp_539 = tmp_537*tmp_67;
      real_t tmp_540 = tmp_537*tmp_68;
      real_t tmp_541 = -tmp_530 - tmp_531 - tmp_533 - tmp_534 - tmp_535 - tmp_536 - tmp_538 - tmp_539 - tmp_540 + 1;
      real_t tmp_542 = tmp_100*tmp_528;
      real_t tmp_543 = 0.019202922745021479*tmp_102;
      real_t tmp_544 = tmp_84 + tmp_85 + tmp_90;
      real_t tmp_545 = 0.5*p_affine_13_0*tmp_68 + 0.5*p_affine_13_1*tmp_71 + 0.5*p_affine_13_2*tmp_74;
      real_t tmp_546 = tmp_117 + tmp_118 + tmp_122;
      real_t tmp_547 = tmp_139 + tmp_140 + tmp_144;
      real_t tmp_548 = tmp_161 + tmp_162 + tmp_166;
      real_t tmp_549 = tmp_183 + tmp_184 + tmp_188;
      real_t tmp_550 = tmp_205 + tmp_206 + tmp_210;
      real_t tmp_551 = tmp_227 + tmp_228 + tmp_232;
      real_t tmp_552 = tmp_249 + tmp_250 + tmp_254;
      real_t tmp_553 = tmp_271 + tmp_272 + tmp_276;
      real_t tmp_554 = tmp_293 + tmp_294 + tmp_298;
      real_t tmp_555 = tmp_315 + tmp_316 + tmp_320;
      real_t tmp_556 = tmp_337 + tmp_338 + tmp_342;
      real_t tmp_557 = tmp_359 + tmp_360 + tmp_364;
      real_t tmp_558 = tmp_381 + tmp_382 + tmp_386;
      real_t tmp_559 = tmp_403 + tmp_404 + tmp_408;
      real_t tmp_560 = tmp_425 + tmp_426 + tmp_430;
      real_t tmp_561 = tmp_447 + tmp_448 + tmp_452;
      real_t tmp_562 = tmp_469 + tmp_470 + tmp_474;
      real_t tmp_563 = tmp_491 + tmp_492 + tmp_496;
      real_t tmp_564 = tmp_513 + tmp_514 + tmp_518;
      real_t tmp_565 = tmp_535 + tmp_536 + tmp_540;
      real_t tmp_566 = tmp_79 + tmp_83 + tmp_89;
      real_t tmp_567 = 0.5*p_affine_13_0*tmp_67 + 0.5*p_affine_13_1*tmp_70 + 0.5*p_affine_13_2*tmp_73;
      real_t tmp_568 = tmp_113 + tmp_116 + tmp_121;
      real_t tmp_569 = tmp_135 + tmp_138 + tmp_143;
      real_t tmp_570 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_571 = tmp_179 + tmp_182 + tmp_187;
      real_t tmp_572 = tmp_201 + tmp_204 + tmp_209;
      real_t tmp_573 = tmp_223 + tmp_226 + tmp_231;
      real_t tmp_574 = tmp_245 + tmp_248 + tmp_253;
      real_t tmp_575 = tmp_267 + tmp_270 + tmp_275;
      real_t tmp_576 = tmp_289 + tmp_292 + tmp_297;
      real_t tmp_577 = tmp_311 + tmp_314 + tmp_319;
      real_t tmp_578 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_579 = tmp_355 + tmp_358 + tmp_363;
      real_t tmp_580 = tmp_377 + tmp_380 + tmp_385;
      real_t tmp_581 = tmp_399 + tmp_402 + tmp_407;
      real_t tmp_582 = tmp_421 + tmp_424 + tmp_429;
      real_t tmp_583 = tmp_443 + tmp_446 + tmp_451;
      real_t tmp_584 = tmp_465 + tmp_468 + tmp_473;
      real_t tmp_585 = tmp_487 + tmp_490 + tmp_495;
      real_t tmp_586 = tmp_509 + tmp_512 + tmp_517;
      real_t tmp_587 = tmp_531 + tmp_534 + tmp_539;
      real_t tmp_588 = tmp_78 + tmp_82 + tmp_88;
      real_t tmp_589 = 0.5*p_affine_13_0*tmp_66 + 0.5*p_affine_13_1*tmp_69 + 0.5*p_affine_13_2*tmp_72;
      real_t tmp_590 = tmp_112 + tmp_115 + tmp_120;
      real_t tmp_591 = tmp_134 + tmp_137 + tmp_142;
      real_t tmp_592 = tmp_156 + tmp_159 + tmp_164;
      real_t tmp_593 = tmp_178 + tmp_181 + tmp_186;
      real_t tmp_594 = tmp_200 + tmp_203 + tmp_208;
      real_t tmp_595 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_596 = tmp_244 + tmp_247 + tmp_252;
      real_t tmp_597 = tmp_266 + tmp_269 + tmp_274;
      real_t tmp_598 = tmp_288 + tmp_291 + tmp_296;
      real_t tmp_599 = tmp_310 + tmp_313 + tmp_318;
      real_t tmp_600 = tmp_332 + tmp_335 + tmp_340;
      real_t tmp_601 = tmp_354 + tmp_357 + tmp_362;
      real_t tmp_602 = tmp_376 + tmp_379 + tmp_384;
      real_t tmp_603 = tmp_398 + tmp_401 + tmp_406;
      real_t tmp_604 = tmp_420 + tmp_423 + tmp_428;
      real_t tmp_605 = tmp_442 + tmp_445 + tmp_450;
      real_t tmp_606 = tmp_464 + tmp_467 + tmp_472;
      real_t tmp_607 = tmp_486 + tmp_489 + tmp_494;
      real_t tmp_608 = tmp_508 + tmp_511 + tmp_516;
      real_t tmp_609 = tmp_530 + tmp_533 + tmp_538;
      real_t a_0_0 = tmp_103*(-tmp_101*tmp_91 - tmp_46*tmp_75 + tmp_91*tmp_95) + tmp_125*(-tmp_110*tmp_75 - tmp_123*tmp_124 + tmp_123*tmp_95) + tmp_147*(-tmp_132*tmp_75 - tmp_145*tmp_146 + tmp_145*tmp_95) + tmp_169*(-tmp_154*tmp_75 - tmp_167*tmp_168 + tmp_167*tmp_95) + tmp_191*(-tmp_176*tmp_75 - tmp_189*tmp_190 + tmp_189*tmp_95) + tmp_213*(-tmp_198*tmp_75 - tmp_211*tmp_212 + tmp_211*tmp_95) + tmp_235*(-tmp_220*tmp_75 - tmp_233*tmp_234 + tmp_233*tmp_95) + tmp_257*(-tmp_242*tmp_75 - tmp_255*tmp_256 + tmp_255*tmp_95) + tmp_279*(-tmp_264*tmp_75 - tmp_277*tmp_278 + tmp_277*tmp_95) + tmp_301*(-tmp_286*tmp_75 - tmp_299*tmp_300 + tmp_299*tmp_95) + tmp_323*(-tmp_308*tmp_75 - tmp_321*tmp_322 + tmp_321*tmp_95) + tmp_345*(-tmp_330*tmp_75 - tmp_343*tmp_344 + tmp_343*tmp_95) + tmp_367*(-tmp_352*tmp_75 - tmp_365*tmp_366 + tmp_365*tmp_95) + tmp_389*(-tmp_374*tmp_75 - tmp_387*tmp_388 + tmp_387*tmp_95) + tmp_411*(-tmp_396*tmp_75 - tmp_409*tmp_410 + tmp_409*tmp_95) + tmp_433*(-tmp_418*tmp_75 - tmp_431*tmp_432 + tmp_431*tmp_95) + tmp_455*(-tmp_440*tmp_75 - tmp_453*tmp_454 + tmp_453*tmp_95) + tmp_477*(-tmp_462*tmp_75 - tmp_475*tmp_476 + tmp_475*tmp_95) + tmp_499*(-tmp_484*tmp_75 - tmp_497*tmp_498 + tmp_497*tmp_95) + tmp_521*(-tmp_506*tmp_75 - tmp_519*tmp_520 + tmp_519*tmp_95) + tmp_543*(-tmp_528*tmp_75 - tmp_541*tmp_542 + tmp_541*tmp_95);
      real_t a_0_1 = tmp_103*(-tmp_101*tmp_544 - tmp_46*tmp_545 + tmp_544*tmp_95) + tmp_125*(-tmp_110*tmp_545 - tmp_124*tmp_546 + tmp_546*tmp_95) + tmp_147*(-tmp_132*tmp_545 - tmp_146*tmp_547 + tmp_547*tmp_95) + tmp_169*(-tmp_154*tmp_545 - tmp_168*tmp_548 + tmp_548*tmp_95) + tmp_191*(-tmp_176*tmp_545 - tmp_190*tmp_549 + tmp_549*tmp_95) + tmp_213*(-tmp_198*tmp_545 - tmp_212*tmp_550 + tmp_550*tmp_95) + tmp_235*(-tmp_220*tmp_545 - tmp_234*tmp_551 + tmp_551*tmp_95) + tmp_257*(-tmp_242*tmp_545 - tmp_256*tmp_552 + tmp_552*tmp_95) + tmp_279*(-tmp_264*tmp_545 - tmp_278*tmp_553 + tmp_553*tmp_95) + tmp_301*(-tmp_286*tmp_545 - tmp_300*tmp_554 + tmp_554*tmp_95) + tmp_323*(-tmp_308*tmp_545 - tmp_322*tmp_555 + tmp_555*tmp_95) + tmp_345*(-tmp_330*tmp_545 - tmp_344*tmp_556 + tmp_556*tmp_95) + tmp_367*(-tmp_352*tmp_545 - tmp_366*tmp_557 + tmp_557*tmp_95) + tmp_389*(-tmp_374*tmp_545 - tmp_388*tmp_558 + tmp_558*tmp_95) + tmp_411*(-tmp_396*tmp_545 - tmp_410*tmp_559 + tmp_559*tmp_95) + tmp_433*(-tmp_418*tmp_545 - tmp_432*tmp_560 + tmp_560*tmp_95) + tmp_455*(-tmp_440*tmp_545 - tmp_454*tmp_561 + tmp_561*tmp_95) + tmp_477*(-tmp_462*tmp_545 - tmp_476*tmp_562 + tmp_562*tmp_95) + tmp_499*(-tmp_484*tmp_545 - tmp_498*tmp_563 + tmp_563*tmp_95) + tmp_521*(-tmp_506*tmp_545 - tmp_520*tmp_564 + tmp_564*tmp_95) + tmp_543*(-tmp_528*tmp_545 - tmp_542*tmp_565 + tmp_565*tmp_95);
      real_t a_0_2 = tmp_103*(-tmp_101*tmp_566 - tmp_46*tmp_567 + tmp_566*tmp_95) + tmp_125*(-tmp_110*tmp_567 - tmp_124*tmp_568 + tmp_568*tmp_95) + tmp_147*(-tmp_132*tmp_567 - tmp_146*tmp_569 + tmp_569*tmp_95) + tmp_169*(-tmp_154*tmp_567 - tmp_168*tmp_570 + tmp_570*tmp_95) + tmp_191*(-tmp_176*tmp_567 - tmp_190*tmp_571 + tmp_571*tmp_95) + tmp_213*(-tmp_198*tmp_567 - tmp_212*tmp_572 + tmp_572*tmp_95) + tmp_235*(-tmp_220*tmp_567 - tmp_234*tmp_573 + tmp_573*tmp_95) + tmp_257*(-tmp_242*tmp_567 - tmp_256*tmp_574 + tmp_574*tmp_95) + tmp_279*(-tmp_264*tmp_567 - tmp_278*tmp_575 + tmp_575*tmp_95) + tmp_301*(-tmp_286*tmp_567 - tmp_300*tmp_576 + tmp_576*tmp_95) + tmp_323*(-tmp_308*tmp_567 - tmp_322*tmp_577 + tmp_577*tmp_95) + tmp_345*(-tmp_330*tmp_567 - tmp_344*tmp_578 + tmp_578*tmp_95) + tmp_367*(-tmp_352*tmp_567 - tmp_366*tmp_579 + tmp_579*tmp_95) + tmp_389*(-tmp_374*tmp_567 - tmp_388*tmp_580 + tmp_580*tmp_95) + tmp_411*(-tmp_396*tmp_567 - tmp_410*tmp_581 + tmp_581*tmp_95) + tmp_433*(-tmp_418*tmp_567 - tmp_432*tmp_582 + tmp_582*tmp_95) + tmp_455*(-tmp_440*tmp_567 - tmp_454*tmp_583 + tmp_583*tmp_95) + tmp_477*(-tmp_462*tmp_567 - tmp_476*tmp_584 + tmp_584*tmp_95) + tmp_499*(-tmp_484*tmp_567 - tmp_498*tmp_585 + tmp_585*tmp_95) + tmp_521*(-tmp_506*tmp_567 - tmp_520*tmp_586 + tmp_586*tmp_95) + tmp_543*(-tmp_528*tmp_567 - tmp_542*tmp_587 + tmp_587*tmp_95);
      real_t a_0_3 = tmp_103*(-tmp_101*tmp_588 - tmp_46*tmp_589 + tmp_588*tmp_95) + tmp_125*(-tmp_110*tmp_589 - tmp_124*tmp_590 + tmp_590*tmp_95) + tmp_147*(-tmp_132*tmp_589 - tmp_146*tmp_591 + tmp_591*tmp_95) + tmp_169*(-tmp_154*tmp_589 - tmp_168*tmp_592 + tmp_592*tmp_95) + tmp_191*(-tmp_176*tmp_589 - tmp_190*tmp_593 + tmp_593*tmp_95) + tmp_213*(-tmp_198*tmp_589 - tmp_212*tmp_594 + tmp_594*tmp_95) + tmp_235*(-tmp_220*tmp_589 - tmp_234*tmp_595 + tmp_595*tmp_95) + tmp_257*(-tmp_242*tmp_589 - tmp_256*tmp_596 + tmp_596*tmp_95) + tmp_279*(-tmp_264*tmp_589 - tmp_278*tmp_597 + tmp_597*tmp_95) + tmp_301*(-tmp_286*tmp_589 - tmp_300*tmp_598 + tmp_598*tmp_95) + tmp_323*(-tmp_308*tmp_589 - tmp_322*tmp_599 + tmp_599*tmp_95) + tmp_345*(-tmp_330*tmp_589 - tmp_344*tmp_600 + tmp_600*tmp_95) + tmp_367*(-tmp_352*tmp_589 - tmp_366*tmp_601 + tmp_601*tmp_95) + tmp_389*(-tmp_374*tmp_589 - tmp_388*tmp_602 + tmp_602*tmp_95) + tmp_411*(-tmp_396*tmp_589 - tmp_410*tmp_603 + tmp_603*tmp_95) + tmp_433*(-tmp_418*tmp_589 - tmp_432*tmp_604 + tmp_604*tmp_95) + tmp_455*(-tmp_440*tmp_589 - tmp_454*tmp_605 + tmp_605*tmp_95) + tmp_477*(-tmp_462*tmp_589 - tmp_476*tmp_606 + tmp_606*tmp_95) + tmp_499*(-tmp_484*tmp_589 - tmp_498*tmp_607 + tmp_607*tmp_95) + tmp_521*(-tmp_506*tmp_589 - tmp_520*tmp_608 + tmp_608*tmp_95) + tmp_543*(-tmp_528*tmp_589 - tmp_542*tmp_609 + tmp_609*tmp_95);
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
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
      real_t tmp_22 = tmp_18*(tmp_12*tmp_6 - tmp_3*tmp_9);
      real_t tmp_23 = tmp_18*(tmp_10*tmp_9 - tmp_15*tmp_6);
      real_t tmp_24 = tmp_18*(-tmp_10*tmp_12 + tmp_15*tmp_3);
      real_t tmp_25 = tmp_18*(-tmp_1*tmp_12 + tmp_5*tmp_9);
      real_t tmp_26 = tmp_18*(tmp_1*tmp_15 - tmp_13*tmp_9);
      real_t tmp_27 = tmp_18*(tmp_12*tmp_13 - tmp_15*tmp_5);
      real_t tmp_28 = -p_affine_8_0;
      real_t tmp_29 = p_affine_10_0 + tmp_28;
      real_t tmp_30 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_31 = -p_affine_8_1;
      real_t tmp_32 = p_affine_10_1 + tmp_31;
      real_t tmp_33 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_34 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_35 = -p_affine_8_2;
      real_t tmp_36 = p_affine_10_2 + tmp_35;
      real_t tmp_37 = 1.0*std::pow((std::abs(tmp_29*tmp_30 - tmp_32*tmp_33)*std::abs(tmp_29*tmp_30 - tmp_32*tmp_33)) + (std::abs(tmp_29*tmp_34 - tmp_33*tmp_36)*std::abs(tmp_29*tmp_34 - tmp_33*tmp_36)) + (std::abs(tmp_30*tmp_36 - tmp_32*tmp_34)*std::abs(tmp_30*tmp_36 - tmp_32*tmp_34)), 1.0/2.0);
      real_t tmp_38 = tmp_37*(p_affine_13_0*(-tmp_19 - tmp_20 - tmp_21) + p_affine_13_1*(-tmp_22 - tmp_23 - tmp_24) + p_affine_13_2*(-tmp_25 - tmp_26 - tmp_27));
      real_t tmp_39 = p_affine_9_2 + tmp_35;
      real_t tmp_40 = p_affine_8_2 + tmp_2;
      real_t tmp_41 = 0.93718850182767688*tmp_36 + 0.031405749086161582*tmp_39 + tmp_40;
      real_t tmp_42 = p_affine_9_1 + tmp_31;
      real_t tmp_43 = p_affine_8_1 + tmp_0;
      real_t tmp_44 = 0.93718850182767688*tmp_32 + 0.031405749086161582*tmp_42 + tmp_43;
      real_t tmp_45 = p_affine_9_0 + tmp_28;
      real_t tmp_46 = p_affine_8_0 + tmp_8;
      real_t tmp_47 = 0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_45 + tmp_46;
      real_t tmp_48 = 0.0068572537431980923*tmp_10*(tmp_19*tmp_47 + tmp_22*tmp_44 + tmp_25*tmp_41 - 1.0/4.0) + 0.0068572537431980923*tmp_3*(tmp_20*tmp_47 + tmp_23*tmp_44 + tmp_26*tmp_41 - 1.0/4.0) + 0.0068572537431980923*tmp_6*(tmp_21*tmp_47 + tmp_24*tmp_44 + tmp_27*tmp_41 - 1.0/4.0);
      real_t tmp_49 = 0.60796128279561268*tmp_36 + 0.19601935860219369*tmp_39 + tmp_40;
      real_t tmp_50 = 0.60796128279561268*tmp_32 + 0.19601935860219369*tmp_42 + tmp_43;
      real_t tmp_51 = 0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_45 + tmp_46;
      real_t tmp_52 = 0.037198804536718075*tmp_10*(tmp_19*tmp_51 + tmp_22*tmp_50 + tmp_25*tmp_49 - 1.0/4.0) + 0.037198804536718075*tmp_3*(tmp_20*tmp_51 + tmp_23*tmp_50 + tmp_26*tmp_49 - 1.0/4.0) + 0.037198804536718075*tmp_6*(tmp_21*tmp_51 + tmp_24*tmp_50 + tmp_27*tmp_49 - 1.0/4.0);
      real_t tmp_53 = 0.039308471900058539*tmp_36 + 0.37605877282253791*tmp_39 + tmp_40;
      real_t tmp_54 = 0.039308471900058539*tmp_32 + 0.37605877282253791*tmp_42 + tmp_43;
      real_t tmp_55 = 0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_45 + tmp_46;
      real_t tmp_56 = 0.020848748529055869*tmp_10*(tmp_19*tmp_55 + tmp_22*tmp_54 + tmp_25*tmp_53 - 1.0/4.0) + 0.020848748529055869*tmp_3*(tmp_20*tmp_55 + tmp_23*tmp_54 + tmp_26*tmp_53 - 1.0/4.0) + 0.020848748529055869*tmp_6*(tmp_21*tmp_55 + tmp_24*tmp_54 + tmp_27*tmp_53 - 1.0/4.0);
      real_t tmp_57 = 0.1711304259088916*tmp_36 + 0.78764240869137092*tmp_39 + tmp_40;
      real_t tmp_58 = 0.1711304259088916*tmp_32 + 0.78764240869137092*tmp_42 + tmp_43;
      real_t tmp_59 = 0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_45 + tmp_46;
      real_t tmp_60 = 0.019202922745021479*tmp_10*(tmp_19*tmp_59 + tmp_22*tmp_58 + tmp_25*tmp_57 - 1.0/4.0) + 0.019202922745021479*tmp_3*(tmp_20*tmp_59 + tmp_23*tmp_58 + tmp_26*tmp_57 - 1.0/4.0) + 0.019202922745021479*tmp_6*(tmp_21*tmp_59 + tmp_24*tmp_58 + tmp_27*tmp_57 - 1.0/4.0);
      real_t tmp_61 = 0.37605877282253791*tmp_36 + 0.58463275527740355*tmp_39 + tmp_40;
      real_t tmp_62 = 0.37605877282253791*tmp_32 + 0.58463275527740355*tmp_42 + tmp_43;
      real_t tmp_63 = 0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_45 + tmp_46;
      real_t tmp_64 = 0.020848748529055869*tmp_10*(tmp_19*tmp_63 + tmp_22*tmp_62 + tmp_25*tmp_61 - 1.0/4.0) + 0.020848748529055869*tmp_3*(tmp_20*tmp_63 + tmp_23*tmp_62 + tmp_26*tmp_61 - 1.0/4.0) + 0.020848748529055869*tmp_6*(tmp_21*tmp_63 + tmp_24*tmp_62 + tmp_27*tmp_61 - 1.0/4.0);
      real_t tmp_65 = 0.78764240869137092*tmp_36 + 0.041227165399737475*tmp_39 + tmp_40;
      real_t tmp_66 = 0.78764240869137092*tmp_32 + 0.041227165399737475*tmp_42 + tmp_43;
      real_t tmp_67 = 0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_45 + tmp_46;
      real_t tmp_68 = 0.019202922745021479*tmp_10*(tmp_19*tmp_67 + tmp_22*tmp_66 + tmp_25*tmp_65 - 1.0/4.0) + 0.019202922745021479*tmp_3*(tmp_20*tmp_67 + tmp_23*tmp_66 + tmp_26*tmp_65 - 1.0/4.0) + 0.019202922745021479*tmp_6*(tmp_21*tmp_67 + tmp_24*tmp_66 + tmp_27*tmp_65 - 1.0/4.0);
      real_t tmp_69 = 0.58463275527740355*tmp_36 + 0.039308471900058539*tmp_39 + tmp_40;
      real_t tmp_70 = 0.58463275527740355*tmp_32 + 0.039308471900058539*tmp_42 + tmp_43;
      real_t tmp_71 = 0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_45 + tmp_46;
      real_t tmp_72 = 0.020848748529055869*tmp_10*(tmp_19*tmp_71 + tmp_22*tmp_70 + tmp_25*tmp_69 - 1.0/4.0) + 0.020848748529055869*tmp_3*(tmp_20*tmp_71 + tmp_23*tmp_70 + tmp_26*tmp_69 - 1.0/4.0) + 0.020848748529055869*tmp_6*(tmp_21*tmp_71 + tmp_24*tmp_70 + tmp_27*tmp_69 - 1.0/4.0);
      real_t tmp_73 = 0.041227165399737475*tmp_36 + 0.78764240869137092*tmp_39 + tmp_40;
      real_t tmp_74 = 0.041227165399737475*tmp_32 + 0.78764240869137092*tmp_42 + tmp_43;
      real_t tmp_75 = 0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_45 + tmp_46;
      real_t tmp_76 = 0.019202922745021479*tmp_10*(tmp_19*tmp_75 + tmp_22*tmp_74 + tmp_25*tmp_73 - 1.0/4.0) + 0.019202922745021479*tmp_3*(tmp_20*tmp_75 + tmp_23*tmp_74 + tmp_26*tmp_73 - 1.0/4.0) + 0.019202922745021479*tmp_6*(tmp_21*tmp_75 + tmp_24*tmp_74 + tmp_27*tmp_73 - 1.0/4.0);
      real_t tmp_77 = 0.039308471900058539*tmp_36 + 0.58463275527740355*tmp_39 + tmp_40;
      real_t tmp_78 = 0.039308471900058539*tmp_32 + 0.58463275527740355*tmp_42 + tmp_43;
      real_t tmp_79 = 0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_45 + tmp_46;
      real_t tmp_80 = 0.020848748529055869*tmp_10*(tmp_19*tmp_79 + tmp_22*tmp_78 + tmp_25*tmp_77 - 1.0/4.0) + 0.020848748529055869*tmp_3*(tmp_20*tmp_79 + tmp_23*tmp_78 + tmp_26*tmp_77 - 1.0/4.0) + 0.020848748529055869*tmp_6*(tmp_21*tmp_79 + tmp_24*tmp_78 + tmp_27*tmp_77 - 1.0/4.0);
      real_t tmp_81 = 0.78764240869137092*tmp_36 + 0.1711304259088916*tmp_39 + tmp_40;
      real_t tmp_82 = 0.78764240869137092*tmp_32 + 0.1711304259088916*tmp_42 + tmp_43;
      real_t tmp_83 = 0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_45 + tmp_46;
      real_t tmp_84 = 0.019202922745021479*tmp_10*(tmp_19*tmp_83 + tmp_22*tmp_82 + tmp_25*tmp_81 - 1.0/4.0) + 0.019202922745021479*tmp_3*(tmp_20*tmp_83 + tmp_23*tmp_82 + tmp_26*tmp_81 - 1.0/4.0) + 0.019202922745021479*tmp_6*(tmp_21*tmp_83 + tmp_24*tmp_82 + tmp_27*tmp_81 - 1.0/4.0);
      real_t tmp_85 = 0.58463275527740355*tmp_36 + 0.37605877282253791*tmp_39 + tmp_40;
      real_t tmp_86 = 0.58463275527740355*tmp_32 + 0.37605877282253791*tmp_42 + tmp_43;
      real_t tmp_87 = 0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_45 + tmp_46;
      real_t tmp_88 = 0.020848748529055869*tmp_10*(tmp_19*tmp_87 + tmp_22*tmp_86 + tmp_25*tmp_85 - 1.0/4.0) + 0.020848748529055869*tmp_3*(tmp_20*tmp_87 + tmp_23*tmp_86 + tmp_26*tmp_85 - 1.0/4.0) + 0.020848748529055869*tmp_6*(tmp_21*tmp_87 + tmp_24*tmp_86 + tmp_27*tmp_85 - 1.0/4.0);
      real_t tmp_89 = 0.1711304259088916*tmp_36 + 0.041227165399737475*tmp_39 + tmp_40;
      real_t tmp_90 = 0.1711304259088916*tmp_32 + 0.041227165399737475*tmp_42 + tmp_43;
      real_t tmp_91 = 0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_45 + tmp_46;
      real_t tmp_92 = 0.019202922745021479*tmp_10*(tmp_19*tmp_91 + tmp_22*tmp_90 + tmp_25*tmp_89 - 1.0/4.0) + 0.019202922745021479*tmp_3*(tmp_20*tmp_91 + tmp_23*tmp_90 + tmp_26*tmp_89 - 1.0/4.0) + 0.019202922745021479*tmp_6*(tmp_21*tmp_91 + tmp_24*tmp_90 + tmp_27*tmp_89 - 1.0/4.0);
      real_t tmp_93 = 0.19107600050469298*tmp_36 + 0.40446199974765351*tmp_39 + tmp_40;
      real_t tmp_94 = 0.19107600050469298*tmp_32 + 0.40446199974765351*tmp_42 + tmp_43;
      real_t tmp_95 = 0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_45 + tmp_46;
      real_t tmp_96 = 0.042507265838595799*tmp_10*(tmp_19*tmp_95 + tmp_22*tmp_94 + tmp_25*tmp_93 - 1.0/4.0) + 0.042507265838595799*tmp_3*(tmp_20*tmp_95 + tmp_23*tmp_94 + tmp_26*tmp_93 - 1.0/4.0) + 0.042507265838595799*tmp_6*(tmp_21*tmp_95 + tmp_24*tmp_94 + tmp_27*tmp_93 - 1.0/4.0);
      real_t tmp_97 = 0.37605877282253791*tmp_36 + 0.039308471900058539*tmp_39 + tmp_40;
      real_t tmp_98 = 0.37605877282253791*tmp_32 + 0.039308471900058539*tmp_42 + tmp_43;
      real_t tmp_99 = 0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_45 + tmp_46;
      real_t tmp_100 = 0.020848748529055869*tmp_10*(tmp_19*tmp_99 + tmp_22*tmp_98 + tmp_25*tmp_97 - 1.0/4.0) + 0.020848748529055869*tmp_3*(tmp_20*tmp_99 + tmp_23*tmp_98 + tmp_26*tmp_97 - 1.0/4.0) + 0.020848748529055869*tmp_6*(tmp_21*tmp_99 + tmp_24*tmp_98 + tmp_27*tmp_97 - 1.0/4.0);
      real_t tmp_101 = 0.031405749086161582*tmp_36 + 0.93718850182767688*tmp_39 + tmp_40;
      real_t tmp_102 = 0.031405749086161582*tmp_32 + 0.93718850182767688*tmp_42 + tmp_43;
      real_t tmp_103 = 0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_45 + tmp_46;
      real_t tmp_104 = 0.0068572537431980923*tmp_10*(tmp_101*tmp_25 + tmp_102*tmp_22 + tmp_103*tmp_19 - 1.0/4.0) + 0.0068572537431980923*tmp_3*(tmp_101*tmp_26 + tmp_102*tmp_23 + tmp_103*tmp_20 - 1.0/4.0) + 0.0068572537431980923*tmp_6*(tmp_101*tmp_27 + tmp_102*tmp_24 + tmp_103*tmp_21 - 1.0/4.0);
      real_t tmp_105 = 0.19601935860219369*tmp_36 + 0.60796128279561268*tmp_39 + tmp_40;
      real_t tmp_106 = 0.19601935860219369*tmp_32 + 0.60796128279561268*tmp_42 + tmp_43;
      real_t tmp_107 = 0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_45 + tmp_46;
      real_t tmp_108 = 0.037198804536718075*tmp_10*(tmp_105*tmp_25 + tmp_106*tmp_22 + tmp_107*tmp_19 - 1.0/4.0) + 0.037198804536718075*tmp_3*(tmp_105*tmp_26 + tmp_106*tmp_23 + tmp_107*tmp_20 - 1.0/4.0) + 0.037198804536718075*tmp_6*(tmp_105*tmp_27 + tmp_106*tmp_24 + tmp_107*tmp_21 - 1.0/4.0);
      real_t tmp_109 = 0.40446199974765351*tmp_36 + 0.19107600050469298*tmp_39 + tmp_40;
      real_t tmp_110 = 0.40446199974765351*tmp_32 + 0.19107600050469298*tmp_42 + tmp_43;
      real_t tmp_111 = 0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_45 + tmp_46;
      real_t tmp_112 = 0.042507265838595799*tmp_10*(tmp_109*tmp_25 + tmp_110*tmp_22 + tmp_111*tmp_19 - 1.0/4.0) + 0.042507265838595799*tmp_3*(tmp_109*tmp_26 + tmp_110*tmp_23 + tmp_111*tmp_20 - 1.0/4.0) + 0.042507265838595799*tmp_6*(tmp_109*tmp_27 + tmp_110*tmp_24 + tmp_111*tmp_21 - 1.0/4.0);
      real_t tmp_113 = 0.031405749086161582*tmp_36 + 0.031405749086161582*tmp_39 + tmp_40;
      real_t tmp_114 = 0.031405749086161582*tmp_32 + 0.031405749086161582*tmp_42 + tmp_43;
      real_t tmp_115 = 0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_45 + tmp_46;
      real_t tmp_116 = 0.0068572537431980923*tmp_10*(tmp_113*tmp_25 + tmp_114*tmp_22 + tmp_115*tmp_19 - 1.0/4.0) + 0.0068572537431980923*tmp_3*(tmp_113*tmp_26 + tmp_114*tmp_23 + tmp_115*tmp_20 - 1.0/4.0) + 0.0068572537431980923*tmp_6*(tmp_113*tmp_27 + tmp_114*tmp_24 + tmp_115*tmp_21 - 1.0/4.0);
      real_t tmp_117 = 0.19601935860219369*tmp_36 + 0.19601935860219369*tmp_39 + tmp_40;
      real_t tmp_118 = 0.19601935860219369*tmp_32 + 0.19601935860219369*tmp_42 + tmp_43;
      real_t tmp_119 = 0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_45 + tmp_46;
      real_t tmp_120 = 0.037198804536718075*tmp_10*(tmp_117*tmp_25 + tmp_118*tmp_22 + tmp_119*tmp_19 - 1.0/4.0) + 0.037198804536718075*tmp_3*(tmp_117*tmp_26 + tmp_118*tmp_23 + tmp_119*tmp_20 - 1.0/4.0) + 0.037198804536718075*tmp_6*(tmp_117*tmp_27 + tmp_118*tmp_24 + tmp_119*tmp_21 - 1.0/4.0);
      real_t tmp_121 = 0.40446199974765351*tmp_36 + 0.40446199974765351*tmp_39 + tmp_40;
      real_t tmp_122 = 0.40446199974765351*tmp_32 + 0.40446199974765351*tmp_42 + tmp_43;
      real_t tmp_123 = 0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_45 + tmp_46;
      real_t tmp_124 = 0.042507265838595799*tmp_10*(tmp_121*tmp_25 + tmp_122*tmp_22 + tmp_123*tmp_19 - 1.0/4.0) + 0.042507265838595799*tmp_3*(tmp_121*tmp_26 + tmp_122*tmp_23 + tmp_123*tmp_20 - 1.0/4.0) + 0.042507265838595799*tmp_6*(tmp_121*tmp_27 + tmp_122*tmp_24 + tmp_123*tmp_21 - 1.0/4.0);
      real_t tmp_125 = 0.041227165399737475*tmp_36 + 0.1711304259088916*tmp_39 + tmp_40;
      real_t tmp_126 = 0.041227165399737475*tmp_32 + 0.1711304259088916*tmp_42 + tmp_43;
      real_t tmp_127 = 0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_45 + tmp_46;
      real_t tmp_128 = 0.019202922745021479*tmp_10*(tmp_125*tmp_25 + tmp_126*tmp_22 + tmp_127*tmp_19 - 1.0/4.0) + 0.019202922745021479*tmp_3*(tmp_125*tmp_26 + tmp_126*tmp_23 + tmp_127*tmp_20 - 1.0/4.0) + 0.019202922745021479*tmp_6*(tmp_125*tmp_27 + tmp_126*tmp_24 + tmp_127*tmp_21 - 1.0/4.0);
      real_t tmp_129 = tmp_37*(p_affine_13_0*tmp_21 + p_affine_13_1*tmp_24 + p_affine_13_2*tmp_27);
      real_t tmp_130 = tmp_37*(p_affine_13_0*tmp_20 + p_affine_13_1*tmp_23 + p_affine_13_2*tmp_26);
      real_t tmp_131 = tmp_37*(p_affine_13_0*tmp_19 + p_affine_13_1*tmp_22 + p_affine_13_2*tmp_25);
      real_t a_0_0 = -tmp_100*tmp_38 - tmp_104*tmp_38 - tmp_108*tmp_38 - tmp_112*tmp_38 - tmp_116*tmp_38 - tmp_120*tmp_38 - tmp_124*tmp_38 - tmp_128*tmp_38 - tmp_38*tmp_48 - tmp_38*tmp_52 - tmp_38*tmp_56 - tmp_38*tmp_60 - tmp_38*tmp_64 - tmp_38*tmp_68 - tmp_38*tmp_72 - tmp_38*tmp_76 - tmp_38*tmp_80 - tmp_38*tmp_84 - tmp_38*tmp_88 - tmp_38*tmp_92 - tmp_38*tmp_96;
      real_t a_0_1 = -tmp_100*tmp_129 - tmp_104*tmp_129 - tmp_108*tmp_129 - tmp_112*tmp_129 - tmp_116*tmp_129 - tmp_120*tmp_129 - tmp_124*tmp_129 - tmp_128*tmp_129 - tmp_129*tmp_48 - tmp_129*tmp_52 - tmp_129*tmp_56 - tmp_129*tmp_60 - tmp_129*tmp_64 - tmp_129*tmp_68 - tmp_129*tmp_72 - tmp_129*tmp_76 - tmp_129*tmp_80 - tmp_129*tmp_84 - tmp_129*tmp_88 - tmp_129*tmp_92 - tmp_129*tmp_96;
      real_t a_0_2 = -tmp_100*tmp_130 - tmp_104*tmp_130 - tmp_108*tmp_130 - tmp_112*tmp_130 - tmp_116*tmp_130 - tmp_120*tmp_130 - tmp_124*tmp_130 - tmp_128*tmp_130 - tmp_130*tmp_48 - tmp_130*tmp_52 - tmp_130*tmp_56 - tmp_130*tmp_60 - tmp_130*tmp_64 - tmp_130*tmp_68 - tmp_130*tmp_72 - tmp_130*tmp_76 - tmp_130*tmp_80 - tmp_130*tmp_84 - tmp_130*tmp_88 - tmp_130*tmp_92 - tmp_130*tmp_96;
      real_t a_0_3 = -tmp_100*tmp_131 - tmp_104*tmp_131 - tmp_108*tmp_131 - tmp_112*tmp_131 - tmp_116*tmp_131 - tmp_120*tmp_131 - tmp_124*tmp_131 - tmp_128*tmp_131 - tmp_131*tmp_48 - tmp_131*tmp_52 - tmp_131*tmp_56 - tmp_131*tmp_60 - tmp_131*tmp_64 - tmp_131*tmp_68 - tmp_131*tmp_72 - tmp_131*tmp_76 - tmp_131*tmp_80 - tmp_131*tmp_84 - tmp_131*tmp_88 - tmp_131*tmp_92 - tmp_131*tmp_96;
      elMat( 0, 0) = a_0_0;
      elMat( 0, 1) = a_0_1;
      elMat( 0, 2) = a_0_2;
      elMat( 0, 3) = a_0_3;
   }


};




class EGNIPGVectorLaplaceFormEE : public hyteg::dg::DGForm
{
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
      real_t tmp_27 = 5/tmp_2;
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
      real_t tmp_49 = 5/tmp_2;
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
      real_t a_0_0 = 0.59231721264047354*((tmp_1*tmp_14 + tmp_16*tmp_5)*(tmp_1*tmp_14 + tmp_16*tmp_5)) + 1.1965716762484155*((tmp_1*tmp_19 + tmp_20*tmp_5)*(tmp_1*tmp_19 + tmp_20*tmp_5)) + 1.4222222222222225*((tmp_1*tmp_23 + tmp_24*tmp_5)*(tmp_1*tmp_23 + tmp_24*tmp_5)) + 1.1965716762484155*((tmp_1*tmp_27 + tmp_28*tmp_5)*(tmp_1*tmp_27 + tmp_28*tmp_5)) + 0.59231721264047354*((tmp_1*tmp_31 + tmp_32*tmp_5)*(tmp_1*tmp_31 + tmp_32*tmp_5)) + 0.59231721264047354*((tmp_14*tmp_6 + tmp_16*tmp_4)*(tmp_14*tmp_6 + tmp_16*tmp_4)) + 1.1965716762484155*((tmp_19*tmp_6 + tmp_20*tmp_4)*(tmp_19*tmp_6 + tmp_20*tmp_4)) + 1.4222222222222225*((tmp_23*tmp_6 + tmp_24*tmp_4)*(tmp_23*tmp_6 + tmp_24*tmp_4)) + 1.1965716762484155*((tmp_27*tmp_6 + tmp_28*tmp_4)*(tmp_27*tmp_6 + tmp_28*tmp_4)) + 0.59231721264047354*((tmp_31*tmp_6 + tmp_32*tmp_4)*(tmp_31*tmp_6 + tmp_32*tmp_4));
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
   void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
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
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_3_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_3_0 + tmp_0;
      real_t tmp_6 = p_affine_2_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = p_affine_1_0 + tmp_0;
      real_t tmp_10 = -p_affine_0_2;
      real_t tmp_11 = p_affine_3_2 + tmp_10;
      real_t tmp_12 = tmp_6*tmp_9;
      real_t tmp_13 = p_affine_1_2 + tmp_10;
      real_t tmp_14 = p_affine_2_2 + tmp_10;
      real_t tmp_15 = p_affine_1_1 + tmp_2;
      real_t tmp_16 = tmp_15*tmp_5;
      real_t tmp_17 = tmp_3*tmp_9;
      real_t tmp_18 = tmp_1*tmp_15;
      real_t tmp_19 = 1.0 / (tmp_11*tmp_12 - tmp_11*tmp_18 + tmp_13*tmp_4 - tmp_13*tmp_7 + tmp_14*tmp_16 - tmp_14*tmp_17);
      real_t tmp_20 = tmp_19*tmp_9;
      real_t tmp_21 = tmp_16 - tmp_17;
      real_t tmp_22 = tmp_1*tmp_19;
      real_t tmp_23 = tmp_12 - tmp_18;
      real_t tmp_24 = tmp_19*tmp_5;
      real_t tmp_25 = -tmp_1*tmp_11 + tmp_14*tmp_5;
      real_t tmp_26 = tmp_11*tmp_9 - tmp_13*tmp_5;
      real_t tmp_27 = tmp_1*tmp_13 - tmp_14*tmp_9;
      real_t tmp_28 = tmp_11*tmp_6 - tmp_14*tmp_3;
      real_t tmp_29 = -tmp_11*tmp_15 + tmp_13*tmp_3;
      real_t tmp_30 = -tmp_13*tmp_6 + tmp_14*tmp_15;
      real_t tmp_31 = tmp_15*tmp_19;
      real_t tmp_32 = tmp_19*tmp_6;
      real_t tmp_33 = tmp_19*tmp_3;
      real_t tmp_34 = tmp_13*tmp_19;
      real_t tmp_35 = tmp_14*tmp_19;
      real_t tmp_36 = tmp_11*tmp_19;
      real_t tmp_37 = p_affine_0_0*p_affine_1_1;
      real_t tmp_38 = p_affine_0_0*p_affine_1_2;
      real_t tmp_39 = p_affine_2_1*p_affine_3_2;
      real_t tmp_40 = p_affine_0_1*p_affine_1_0;
      real_t tmp_41 = p_affine_0_1*p_affine_1_2;
      real_t tmp_42 = p_affine_2_2*p_affine_3_0;
      real_t tmp_43 = p_affine_0_2*p_affine_1_0;
      real_t tmp_44 = p_affine_0_2*p_affine_1_1;
      real_t tmp_45 = p_affine_2_0*p_affine_3_1;
      real_t tmp_46 = p_affine_2_2*p_affine_3_1;
      real_t tmp_47 = p_affine_2_0*p_affine_3_2;
      real_t tmp_48 = p_affine_2_1*p_affine_3_0;
      real_t tmp_49 = (((tmp_20*tmp_25 + tmp_22*tmp_26 + tmp_24*tmp_27)*(tmp_20*tmp_25 + tmp_22*tmp_26 + tmp_24*tmp_27)) + ((tmp_20*tmp_28 + tmp_22*tmp_29 + tmp_24*tmp_30)*(tmp_20*tmp_28 + tmp_22*tmp_29 + tmp_24*tmp_30)) + ((tmp_20*tmp_8 + tmp_21*tmp_22 + tmp_23*tmp_24)*(tmp_20*tmp_8 + tmp_21*tmp_22 + tmp_23*tmp_24)) + ((tmp_21*tmp_32 + tmp_23*tmp_33 + tmp_31*tmp_8)*(tmp_21*tmp_32 + tmp_23*tmp_33 + tmp_31*tmp_8)) + ((tmp_21*tmp_35 + tmp_23*tmp_36 + tmp_34*tmp_8)*(tmp_21*tmp_35 + tmp_23*tmp_36 + tmp_34*tmp_8)) + ((tmp_25*tmp_31 + tmp_26*tmp_32 + tmp_27*tmp_33)*(tmp_25*tmp_31 + tmp_26*tmp_32 + tmp_27*tmp_33)) + ((tmp_25*tmp_34 + tmp_26*tmp_35 + tmp_27*tmp_36)*(tmp_25*tmp_34 + tmp_26*tmp_35 + tmp_27*tmp_36)) + ((tmp_28*tmp_31 + tmp_29*tmp_32 + tmp_30*tmp_33)*(tmp_28*tmp_31 + tmp_29*tmp_32 + tmp_30*tmp_33)) + ((tmp_28*tmp_34 + tmp_29*tmp_35 + tmp_30*tmp_36)*(tmp_28*tmp_34 + tmp_29*tmp_35 + tmp_30*tmp_36)))*std::abs(p_affine_0_0*tmp_39 - p_affine_0_0*tmp_46 + p_affine_0_1*tmp_42 - p_affine_0_1*tmp_47 + p_affine_0_2*tmp_45 - p_affine_0_2*tmp_48 - p_affine_1_0*tmp_39 + p_affine_1_0*tmp_46 - p_affine_1_1*tmp_42 + p_affine_1_1*tmp_47 - p_affine_1_2*tmp_45 + p_affine_1_2*tmp_48 + p_affine_2_0*tmp_41 - p_affine_2_0*tmp_44 - p_affine_2_1*tmp_38 + p_affine_2_1*tmp_43 + p_affine_2_2*tmp_37 - p_affine_2_2*tmp_40 - p_affine_3_0*tmp_41 + p_affine_3_0*tmp_44 + p_affine_3_1*tmp_38 - p_affine_3_1*tmp_43 - p_affine_3_2*tmp_37 + p_affine_3_2*tmp_40);
      real_t a_0_0 = 0.1666666666666668*tmp_49;
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
      real_t tmp_46 = tmp_1*tmp_37 + tmp_2*tmp_41 + tmp_45*tmp_6;
      real_t tmp_47 = tmp_1*tmp_19;
      real_t tmp_48 = tmp_19*tmp_2;
      real_t tmp_49 = tmp_19*tmp_6;
      real_t tmp_50 = 1.0*p_affine_13_0*(tmp_31*tmp_47 + tmp_40*tmp_48 + tmp_44*tmp_49) + 1.0*p_affine_13_1*(tmp_25*tmp_47 + tmp_39*tmp_48 + tmp_43*tmp_49) + 1.0*p_affine_13_2*(tmp_38*tmp_48 + tmp_42*tmp_49 + tmp_47*tmp_9);
      real_t tmp_51 = tmp_14*tmp_37 + tmp_4*tmp_45 + tmp_41*tmp_7;
      real_t tmp_52 = tmp_14*tmp_19;
      real_t tmp_53 = tmp_19*tmp_7;
      real_t tmp_54 = tmp_19*tmp_4;
      real_t tmp_55 = 1.0*p_affine_13_0*(tmp_31*tmp_52 + tmp_40*tmp_53 + tmp_44*tmp_54) + 1.0*p_affine_13_1*(tmp_25*tmp_52 + tmp_39*tmp_53 + tmp_43*tmp_54) + 1.0*p_affine_13_2*(tmp_38*tmp_53 + tmp_42*tmp_54 + tmp_52*tmp_9);
      real_t tmp_56 = tmp_11*tmp_45 + tmp_13*tmp_37 + tmp_15*tmp_41;
      real_t tmp_57 = tmp_13*tmp_19;
      real_t tmp_58 = tmp_15*tmp_19;
      real_t tmp_59 = tmp_11*tmp_19;
      real_t tmp_60 = 1.0*p_affine_13_0*(tmp_31*tmp_57 + tmp_40*tmp_58 + tmp_44*tmp_59) + 1.0*p_affine_13_1*(tmp_25*tmp_57 + tmp_39*tmp_58 + tmp_43*tmp_59) + 1.0*p_affine_13_2*(tmp_38*tmp_58 + tmp_42*tmp_59 + tmp_57*tmp_9);
      real_t tmp_61 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_62 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_63 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_64 = (std::abs(tmp_22*tmp_61 - tmp_28*tmp_63)*std::abs(tmp_22*tmp_61 - tmp_28*tmp_63)) + (std::abs(tmp_22*tmp_62 - tmp_34*tmp_63)*std::abs(tmp_22*tmp_62 - tmp_34*tmp_63)) + (std::abs(tmp_28*tmp_62 - tmp_34*tmp_61)*std::abs(tmp_28*tmp_62 - tmp_34*tmp_61));
      real_t tmp_65 = 5.0*std::pow(tmp_64, -0.25);
      real_t tmp_66 = 1.0*std::pow(tmp_64, 1.0/2.0);
      real_t tmp_67 = tmp_19*(0.19601935860219369*tmp_21 + 0.60796128279561268*tmp_22 + tmp_23);
      real_t tmp_68 = tmp_19*(0.19601935860219369*tmp_27 + 0.60796128279561268*tmp_28 + tmp_29);
      real_t tmp_69 = tmp_19*(0.19601935860219369*tmp_33 + 0.60796128279561268*tmp_34 + tmp_35);
      real_t tmp_70 = tmp_25*tmp_68 + tmp_31*tmp_69 + tmp_67*tmp_9 - 1.0/4.0;
      real_t tmp_71 = tmp_38*tmp_67 + tmp_39*tmp_68 + tmp_40*tmp_69 - 1.0/4.0;
      real_t tmp_72 = tmp_42*tmp_67 + tmp_43*tmp_68 + tmp_44*tmp_69 - 1.0/4.0;
      real_t tmp_73 = tmp_1*tmp_70 + tmp_2*tmp_71 + tmp_6*tmp_72;
      real_t tmp_74 = tmp_14*tmp_70 + tmp_4*tmp_72 + tmp_7*tmp_71;
      real_t tmp_75 = tmp_11*tmp_72 + tmp_13*tmp_70 + tmp_15*tmp_71;
      real_t tmp_76 = tmp_19*(0.37605877282253791*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_77 = tmp_19*(0.37605877282253791*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29);
      real_t tmp_78 = tmp_19*(0.37605877282253791*tmp_33 + 0.039308471900058539*tmp_34 + tmp_35);
      real_t tmp_79 = tmp_25*tmp_77 + tmp_31*tmp_78 + tmp_76*tmp_9 - 1.0/4.0;
      real_t tmp_80 = tmp_38*tmp_76 + tmp_39*tmp_77 + tmp_40*tmp_78 - 1.0/4.0;
      real_t tmp_81 = tmp_42*tmp_76 + tmp_43*tmp_77 + tmp_44*tmp_78 - 1.0/4.0;
      real_t tmp_82 = tmp_1*tmp_79 + tmp_2*tmp_80 + tmp_6*tmp_81;
      real_t tmp_83 = tmp_14*tmp_79 + tmp_4*tmp_81 + tmp_7*tmp_80;
      real_t tmp_84 = tmp_11*tmp_81 + tmp_13*tmp_79 + tmp_15*tmp_80;
      real_t tmp_85 = tmp_19*(0.78764240869137092*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_86 = tmp_19*(0.78764240869137092*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29);
      real_t tmp_87 = tmp_19*(0.78764240869137092*tmp_33 + 0.1711304259088916*tmp_34 + tmp_35);
      real_t tmp_88 = tmp_25*tmp_86 + tmp_31*tmp_87 + tmp_85*tmp_9 - 1.0/4.0;
      real_t tmp_89 = tmp_38*tmp_85 + tmp_39*tmp_86 + tmp_40*tmp_87 - 1.0/4.0;
      real_t tmp_90 = tmp_42*tmp_85 + tmp_43*tmp_86 + tmp_44*tmp_87 - 1.0/4.0;
      real_t tmp_91 = tmp_1*tmp_88 + tmp_2*tmp_89 + tmp_6*tmp_90;
      real_t tmp_92 = tmp_14*tmp_88 + tmp_4*tmp_90 + tmp_7*tmp_89;
      real_t tmp_93 = tmp_11*tmp_90 + tmp_13*tmp_88 + tmp_15*tmp_89;
      real_t tmp_94 = tmp_19*(0.58463275527740355*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_95 = tmp_19*(0.58463275527740355*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29);
      real_t tmp_96 = tmp_19*(0.58463275527740355*tmp_33 + 0.37605877282253791*tmp_34 + tmp_35);
      real_t tmp_97 = tmp_25*tmp_95 + tmp_31*tmp_96 + tmp_9*tmp_94 - 1.0/4.0;
      real_t tmp_98 = tmp_38*tmp_94 + tmp_39*tmp_95 + tmp_40*tmp_96 - 1.0/4.0;
      real_t tmp_99 = tmp_42*tmp_94 + tmp_43*tmp_95 + tmp_44*tmp_96 - 1.0/4.0;
      real_t tmp_100 = tmp_1*tmp_97 + tmp_2*tmp_98 + tmp_6*tmp_99;
      real_t tmp_101 = tmp_14*tmp_97 + tmp_4*tmp_99 + tmp_7*tmp_98;
      real_t tmp_102 = tmp_11*tmp_99 + tmp_13*tmp_97 + tmp_15*tmp_98;
      real_t tmp_103 = tmp_19*(0.041227165399737475*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_104 = tmp_19*(0.041227165399737475*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29);
      real_t tmp_105 = tmp_19*(0.041227165399737475*tmp_33 + 0.78764240869137092*tmp_34 + tmp_35);
      real_t tmp_106 = tmp_103*tmp_9 + tmp_104*tmp_25 + tmp_105*tmp_31 - 1.0/4.0;
      real_t tmp_107 = tmp_103*tmp_38 + tmp_104*tmp_39 + tmp_105*tmp_40 - 1.0/4.0;
      real_t tmp_108 = tmp_103*tmp_42 + tmp_104*tmp_43 + tmp_105*tmp_44 - 1.0/4.0;
      real_t tmp_109 = tmp_1*tmp_106 + tmp_107*tmp_2 + tmp_108*tmp_6;
      real_t tmp_110 = tmp_106*tmp_14 + tmp_107*tmp_7 + tmp_108*tmp_4;
      real_t tmp_111 = tmp_106*tmp_13 + tmp_107*tmp_15 + tmp_108*tmp_11;
      real_t tmp_112 = tmp_19*(0.039308471900058539*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_113 = tmp_19*(0.039308471900058539*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29);
      real_t tmp_114 = tmp_19*(0.039308471900058539*tmp_33 + 0.58463275527740355*tmp_34 + tmp_35);
      real_t tmp_115 = tmp_112*tmp_9 + tmp_113*tmp_25 + tmp_114*tmp_31 - 1.0/4.0;
      real_t tmp_116 = tmp_112*tmp_38 + tmp_113*tmp_39 + tmp_114*tmp_40 - 1.0/4.0;
      real_t tmp_117 = tmp_112*tmp_42 + tmp_113*tmp_43 + tmp_114*tmp_44 - 1.0/4.0;
      real_t tmp_118 = tmp_1*tmp_115 + tmp_116*tmp_2 + tmp_117*tmp_6;
      real_t tmp_119 = tmp_115*tmp_14 + tmp_116*tmp_7 + tmp_117*tmp_4;
      real_t tmp_120 = tmp_11*tmp_117 + tmp_115*tmp_13 + tmp_116*tmp_15;
      real_t tmp_121 = tmp_19*(0.78764240869137092*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_122 = tmp_19*(0.78764240869137092*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29);
      real_t tmp_123 = tmp_19*(0.78764240869137092*tmp_33 + 0.041227165399737475*tmp_34 + tmp_35);
      real_t tmp_124 = tmp_121*tmp_9 + tmp_122*tmp_25 + tmp_123*tmp_31 - 1.0/4.0;
      real_t tmp_125 = tmp_121*tmp_38 + tmp_122*tmp_39 + tmp_123*tmp_40 - 1.0/4.0;
      real_t tmp_126 = tmp_121*tmp_42 + tmp_122*tmp_43 + tmp_123*tmp_44 - 1.0/4.0;
      real_t tmp_127 = tmp_1*tmp_124 + tmp_125*tmp_2 + tmp_126*tmp_6;
      real_t tmp_128 = tmp_124*tmp_14 + tmp_125*tmp_7 + tmp_126*tmp_4;
      real_t tmp_129 = tmp_11*tmp_126 + tmp_124*tmp_13 + tmp_125*tmp_15;
      real_t tmp_130 = tmp_19*(0.58463275527740355*tmp_21 + 0.039308471900058539*tmp_22 + tmp_23);
      real_t tmp_131 = tmp_19*(0.58463275527740355*tmp_27 + 0.039308471900058539*tmp_28 + tmp_29);
      real_t tmp_132 = tmp_19*(0.58463275527740355*tmp_33 + 0.039308471900058539*tmp_34 + tmp_35);
      real_t tmp_133 = tmp_130*tmp_9 + tmp_131*tmp_25 + tmp_132*tmp_31 - 1.0/4.0;
      real_t tmp_134 = tmp_130*tmp_38 + tmp_131*tmp_39 + tmp_132*tmp_40 - 1.0/4.0;
      real_t tmp_135 = tmp_130*tmp_42 + tmp_131*tmp_43 + tmp_132*tmp_44 - 1.0/4.0;
      real_t tmp_136 = tmp_1*tmp_133 + tmp_134*tmp_2 + tmp_135*tmp_6;
      real_t tmp_137 = tmp_133*tmp_14 + tmp_134*tmp_7 + tmp_135*tmp_4;
      real_t tmp_138 = tmp_11*tmp_135 + tmp_13*tmp_133 + tmp_134*tmp_15;
      real_t tmp_139 = tmp_19*(0.1711304259088916*tmp_21 + 0.78764240869137092*tmp_22 + tmp_23);
      real_t tmp_140 = tmp_19*(0.1711304259088916*tmp_27 + 0.78764240869137092*tmp_28 + tmp_29);
      real_t tmp_141 = tmp_19*(0.1711304259088916*tmp_33 + 0.78764240869137092*tmp_34 + tmp_35);
      real_t tmp_142 = tmp_139*tmp_9 + tmp_140*tmp_25 + tmp_141*tmp_31 - 1.0/4.0;
      real_t tmp_143 = tmp_139*tmp_38 + tmp_140*tmp_39 + tmp_141*tmp_40 - 1.0/4.0;
      real_t tmp_144 = tmp_139*tmp_42 + tmp_140*tmp_43 + tmp_141*tmp_44 - 1.0/4.0;
      real_t tmp_145 = tmp_1*tmp_142 + tmp_143*tmp_2 + tmp_144*tmp_6;
      real_t tmp_146 = tmp_14*tmp_142 + tmp_143*tmp_7 + tmp_144*tmp_4;
      real_t tmp_147 = tmp_11*tmp_144 + tmp_13*tmp_142 + tmp_143*tmp_15;
      real_t tmp_148 = tmp_19*(0.37605877282253791*tmp_21 + 0.58463275527740355*tmp_22 + tmp_23);
      real_t tmp_149 = tmp_19*(0.37605877282253791*tmp_27 + 0.58463275527740355*tmp_28 + tmp_29);
      real_t tmp_150 = tmp_19*(0.37605877282253791*tmp_33 + 0.58463275527740355*tmp_34 + tmp_35);
      real_t tmp_151 = tmp_148*tmp_9 + tmp_149*tmp_25 + tmp_150*tmp_31 - 1.0/4.0;
      real_t tmp_152 = tmp_148*tmp_38 + tmp_149*tmp_39 + tmp_150*tmp_40 - 1.0/4.0;
      real_t tmp_153 = tmp_148*tmp_42 + tmp_149*tmp_43 + tmp_150*tmp_44 - 1.0/4.0;
      real_t tmp_154 = tmp_1*tmp_151 + tmp_152*tmp_2 + tmp_153*tmp_6;
      real_t tmp_155 = tmp_14*tmp_151 + tmp_152*tmp_7 + tmp_153*tmp_4;
      real_t tmp_156 = tmp_11*tmp_153 + tmp_13*tmp_151 + tmp_15*tmp_152;
      real_t tmp_157 = tmp_19*(0.041227165399737475*tmp_21 + 0.1711304259088916*tmp_22 + tmp_23);
      real_t tmp_158 = tmp_19*(0.041227165399737475*tmp_27 + 0.1711304259088916*tmp_28 + tmp_29);
      real_t tmp_159 = tmp_19*(0.041227165399737475*tmp_33 + 0.1711304259088916*tmp_34 + tmp_35);
      real_t tmp_160 = tmp_157*tmp_9 + tmp_158*tmp_25 + tmp_159*tmp_31 - 1.0/4.0;
      real_t tmp_161 = tmp_157*tmp_38 + tmp_158*tmp_39 + tmp_159*tmp_40 - 1.0/4.0;
      real_t tmp_162 = tmp_157*tmp_42 + tmp_158*tmp_43 + tmp_159*tmp_44 - 1.0/4.0;
      real_t tmp_163 = tmp_1*tmp_160 + tmp_161*tmp_2 + tmp_162*tmp_6;
      real_t tmp_164 = tmp_14*tmp_160 + tmp_161*tmp_7 + tmp_162*tmp_4;
      real_t tmp_165 = tmp_11*tmp_162 + tmp_13*tmp_160 + tmp_15*tmp_161;
      real_t tmp_166 = tmp_19*(0.40446199974765351*tmp_21 + 0.19107600050469298*tmp_22 + tmp_23);
      real_t tmp_167 = tmp_19*(0.40446199974765351*tmp_27 + 0.19107600050469298*tmp_28 + tmp_29);
      real_t tmp_168 = tmp_19*(0.40446199974765351*tmp_33 + 0.19107600050469298*tmp_34 + tmp_35);
      real_t tmp_169 = tmp_166*tmp_9 + tmp_167*tmp_25 + tmp_168*tmp_31 - 1.0/4.0;
      real_t tmp_170 = tmp_166*tmp_38 + tmp_167*tmp_39 + tmp_168*tmp_40 - 1.0/4.0;
      real_t tmp_171 = tmp_166*tmp_42 + tmp_167*tmp_43 + tmp_168*tmp_44 - 1.0/4.0;
      real_t tmp_172 = tmp_1*tmp_169 + tmp_170*tmp_2 + tmp_171*tmp_6;
      real_t tmp_173 = tmp_14*tmp_169 + tmp_170*tmp_7 + tmp_171*tmp_4;
      real_t tmp_174 = tmp_11*tmp_171 + tmp_13*tmp_169 + tmp_15*tmp_170;
      real_t tmp_175 = tmp_19*(0.039308471900058539*tmp_21 + 0.37605877282253791*tmp_22 + tmp_23);
      real_t tmp_176 = tmp_19*(0.039308471900058539*tmp_27 + 0.37605877282253791*tmp_28 + tmp_29);
      real_t tmp_177 = tmp_19*(0.039308471900058539*tmp_33 + 0.37605877282253791*tmp_34 + tmp_35);
      real_t tmp_178 = tmp_175*tmp_9 + tmp_176*tmp_25 + tmp_177*tmp_31 - 1.0/4.0;
      real_t tmp_179 = tmp_175*tmp_38 + tmp_176*tmp_39 + tmp_177*tmp_40 - 1.0/4.0;
      real_t tmp_180 = tmp_175*tmp_42 + tmp_176*tmp_43 + tmp_177*tmp_44 - 1.0/4.0;
      real_t tmp_181 = tmp_1*tmp_178 + tmp_179*tmp_2 + tmp_180*tmp_6;
      real_t tmp_182 = tmp_14*tmp_178 + tmp_179*tmp_7 + tmp_180*tmp_4;
      real_t tmp_183 = tmp_11*tmp_180 + tmp_13*tmp_178 + tmp_15*tmp_179;
      real_t tmp_184 = tmp_19*(0.93718850182767688*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_185 = tmp_19*(0.93718850182767688*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29);
      real_t tmp_186 = tmp_19*(0.93718850182767688*tmp_33 + 0.031405749086161582*tmp_34 + tmp_35);
      real_t tmp_187 = tmp_184*tmp_9 + tmp_185*tmp_25 + tmp_186*tmp_31 - 1.0/4.0;
      real_t tmp_188 = tmp_184*tmp_38 + tmp_185*tmp_39 + tmp_186*tmp_40 - 1.0/4.0;
      real_t tmp_189 = tmp_184*tmp_42 + tmp_185*tmp_43 + tmp_186*tmp_44 - 1.0/4.0;
      real_t tmp_190 = tmp_1*tmp_187 + tmp_188*tmp_2 + tmp_189*tmp_6;
      real_t tmp_191 = tmp_14*tmp_187 + tmp_188*tmp_7 + tmp_189*tmp_4;
      real_t tmp_192 = tmp_11*tmp_189 + tmp_13*tmp_187 + tmp_15*tmp_188;
      real_t tmp_193 = tmp_19*(0.60796128279561268*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_194 = tmp_19*(0.60796128279561268*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29);
      real_t tmp_195 = tmp_19*(0.60796128279561268*tmp_33 + 0.19601935860219369*tmp_34 + tmp_35);
      real_t tmp_196 = tmp_193*tmp_9 + tmp_194*tmp_25 + tmp_195*tmp_31 - 1.0/4.0;
      real_t tmp_197 = tmp_193*tmp_38 + tmp_194*tmp_39 + tmp_195*tmp_40 - 1.0/4.0;
      real_t tmp_198 = tmp_193*tmp_42 + tmp_194*tmp_43 + tmp_195*tmp_44 - 1.0/4.0;
      real_t tmp_199 = tmp_1*tmp_196 + tmp_197*tmp_2 + tmp_198*tmp_6;
      real_t tmp_200 = tmp_14*tmp_196 + tmp_197*tmp_7 + tmp_198*tmp_4;
      real_t tmp_201 = tmp_11*tmp_198 + tmp_13*tmp_196 + tmp_15*tmp_197;
      real_t tmp_202 = tmp_19*(0.19107600050469298*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_203 = tmp_19*(0.19107600050469298*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29);
      real_t tmp_204 = tmp_19*(0.19107600050469298*tmp_33 + 0.40446199974765351*tmp_34 + tmp_35);
      real_t tmp_205 = tmp_202*tmp_9 + tmp_203*tmp_25 + tmp_204*tmp_31 - 1.0/4.0;
      real_t tmp_206 = tmp_202*tmp_38 + tmp_203*tmp_39 + tmp_204*tmp_40 - 1.0/4.0;
      real_t tmp_207 = tmp_202*tmp_42 + tmp_203*tmp_43 + tmp_204*tmp_44 - 1.0/4.0;
      real_t tmp_208 = tmp_1*tmp_205 + tmp_2*tmp_206 + tmp_207*tmp_6;
      real_t tmp_209 = tmp_14*tmp_205 + tmp_206*tmp_7 + tmp_207*tmp_4;
      real_t tmp_210 = tmp_11*tmp_207 + tmp_13*tmp_205 + tmp_15*tmp_206;
      real_t tmp_211 = tmp_19*(0.031405749086161582*tmp_21 + 0.031405749086161582*tmp_22 + tmp_23);
      real_t tmp_212 = tmp_19*(0.031405749086161582*tmp_27 + 0.031405749086161582*tmp_28 + tmp_29);
      real_t tmp_213 = tmp_19*(0.031405749086161582*tmp_33 + 0.031405749086161582*tmp_34 + tmp_35);
      real_t tmp_214 = tmp_211*tmp_9 + tmp_212*tmp_25 + tmp_213*tmp_31 - 1.0/4.0;
      real_t tmp_215 = tmp_211*tmp_38 + tmp_212*tmp_39 + tmp_213*tmp_40 - 1.0/4.0;
      real_t tmp_216 = tmp_211*tmp_42 + tmp_212*tmp_43 + tmp_213*tmp_44 - 1.0/4.0;
      real_t tmp_217 = tmp_1*tmp_214 + tmp_2*tmp_215 + tmp_216*tmp_6;
      real_t tmp_218 = tmp_14*tmp_214 + tmp_215*tmp_7 + tmp_216*tmp_4;
      real_t tmp_219 = tmp_11*tmp_216 + tmp_13*tmp_214 + tmp_15*tmp_215;
      real_t tmp_220 = tmp_19*(0.19601935860219369*tmp_21 + 0.19601935860219369*tmp_22 + tmp_23);
      real_t tmp_221 = tmp_19*(0.19601935860219369*tmp_27 + 0.19601935860219369*tmp_28 + tmp_29);
      real_t tmp_222 = tmp_19*(0.19601935860219369*tmp_33 + 0.19601935860219369*tmp_34 + tmp_35);
      real_t tmp_223 = tmp_220*tmp_9 + tmp_221*tmp_25 + tmp_222*tmp_31 - 1.0/4.0;
      real_t tmp_224 = tmp_220*tmp_38 + tmp_221*tmp_39 + tmp_222*tmp_40 - 1.0/4.0;
      real_t tmp_225 = tmp_220*tmp_42 + tmp_221*tmp_43 + tmp_222*tmp_44 - 1.0/4.0;
      real_t tmp_226 = tmp_1*tmp_223 + tmp_2*tmp_224 + tmp_225*tmp_6;
      real_t tmp_227 = tmp_14*tmp_223 + tmp_224*tmp_7 + tmp_225*tmp_4;
      real_t tmp_228 = tmp_11*tmp_225 + tmp_13*tmp_223 + tmp_15*tmp_224;
      real_t tmp_229 = tmp_19*(0.40446199974765351*tmp_21 + 0.40446199974765351*tmp_22 + tmp_23);
      real_t tmp_230 = tmp_19*(0.40446199974765351*tmp_27 + 0.40446199974765351*tmp_28 + tmp_29);
      real_t tmp_231 = tmp_19*(0.40446199974765351*tmp_33 + 0.40446199974765351*tmp_34 + tmp_35);
      real_t tmp_232 = tmp_229*tmp_9 + tmp_230*tmp_25 + tmp_231*tmp_31 - 1.0/4.0;
      real_t tmp_233 = tmp_229*tmp_38 + tmp_230*tmp_39 + tmp_231*tmp_40 - 1.0/4.0;
      real_t tmp_234 = tmp_229*tmp_42 + tmp_230*tmp_43 + tmp_231*tmp_44 - 1.0/4.0;
      real_t tmp_235 = tmp_1*tmp_232 + tmp_2*tmp_233 + tmp_234*tmp_6;
      real_t tmp_236 = tmp_14*tmp_232 + tmp_233*tmp_7 + tmp_234*tmp_4;
      real_t tmp_237 = tmp_11*tmp_234 + tmp_13*tmp_232 + tmp_15*tmp_233;
      real_t tmp_238 = tmp_19*(0.1711304259088916*tmp_21 + 0.041227165399737475*tmp_22 + tmp_23);
      real_t tmp_239 = tmp_19*(0.1711304259088916*tmp_27 + 0.041227165399737475*tmp_28 + tmp_29);
      real_t tmp_240 = tmp_19*(0.1711304259088916*tmp_33 + 0.041227165399737475*tmp_34 + tmp_35);
      real_t tmp_241 = tmp_238*tmp_9 + tmp_239*tmp_25 + tmp_240*tmp_31 - 1.0/4.0;
      real_t tmp_242 = tmp_238*tmp_38 + tmp_239*tmp_39 + tmp_240*tmp_40 - 1.0/4.0;
      real_t tmp_243 = tmp_238*tmp_42 + tmp_239*tmp_43 + tmp_240*tmp_44 - 1.0/4.0;
      real_t tmp_244 = tmp_1*tmp_241 + tmp_2*tmp_242 + tmp_243*tmp_6;
      real_t tmp_245 = tmp_14*tmp_241 + tmp_242*tmp_7 + tmp_243*tmp_4;
      real_t tmp_246 = tmp_11*tmp_243 + tmp_13*tmp_241 + tmp_15*tmp_242;
      real_t a_0_0 = 0.020848748529055869*tmp_66*(-tmp_100*tmp_50 - tmp_101*tmp_55 - tmp_102*tmp_60 + tmp_65*((tmp_100*tmp_100) + (tmp_101*tmp_101) + (tmp_102*tmp_102))) + 0.019202922745021479*tmp_66*(-tmp_109*tmp_50 - tmp_110*tmp_55 - tmp_111*tmp_60 + tmp_65*((tmp_109*tmp_109) + (tmp_110*tmp_110) + (tmp_111*tmp_111))) + 0.020848748529055869*tmp_66*(-tmp_118*tmp_50 - tmp_119*tmp_55 - tmp_120*tmp_60 + tmp_65*((tmp_118*tmp_118) + (tmp_119*tmp_119) + (tmp_120*tmp_120))) + 0.019202922745021479*tmp_66*(-tmp_127*tmp_50 - tmp_128*tmp_55 - tmp_129*tmp_60 + tmp_65*((tmp_127*tmp_127) + (tmp_128*tmp_128) + (tmp_129*tmp_129))) + 0.020848748529055869*tmp_66*(-tmp_136*tmp_50 - tmp_137*tmp_55 - tmp_138*tmp_60 + tmp_65*((tmp_136*tmp_136) + (tmp_137*tmp_137) + (tmp_138*tmp_138))) + 0.019202922745021479*tmp_66*(-tmp_145*tmp_50 - tmp_146*tmp_55 - tmp_147*tmp_60 + tmp_65*((tmp_145*tmp_145) + (tmp_146*tmp_146) + (tmp_147*tmp_147))) + 0.020848748529055869*tmp_66*(-tmp_154*tmp_50 - tmp_155*tmp_55 - tmp_156*tmp_60 + tmp_65*((tmp_154*tmp_154) + (tmp_155*tmp_155) + (tmp_156*tmp_156))) + 0.019202922745021479*tmp_66*(-tmp_163*tmp_50 - tmp_164*tmp_55 - tmp_165*tmp_60 + tmp_65*((tmp_163*tmp_163) + (tmp_164*tmp_164) + (tmp_165*tmp_165))) + 0.042507265838595799*tmp_66*(-tmp_172*tmp_50 - tmp_173*tmp_55 - tmp_174*tmp_60 + tmp_65*((tmp_172*tmp_172) + (tmp_173*tmp_173) + (tmp_174*tmp_174))) + 0.020848748529055869*tmp_66*(-tmp_181*tmp_50 - tmp_182*tmp_55 - tmp_183*tmp_60 + tmp_65*((tmp_181*tmp_181) + (tmp_182*tmp_182) + (tmp_183*tmp_183))) + 0.0068572537431980923*tmp_66*(-tmp_190*tmp_50 - tmp_191*tmp_55 - tmp_192*tmp_60 + tmp_65*((tmp_190*tmp_190) + (tmp_191*tmp_191) + (tmp_192*tmp_192))) + 0.037198804536718075*tmp_66*(-tmp_199*tmp_50 - tmp_200*tmp_55 - tmp_201*tmp_60 + tmp_65*((tmp_199*tmp_199) + (tmp_200*tmp_200) + (tmp_201*tmp_201))) + 0.042507265838595799*tmp_66*(-tmp_208*tmp_50 - tmp_209*tmp_55 - tmp_210*tmp_60 + tmp_65*((tmp_208*tmp_208) + (tmp_209*tmp_209) + (tmp_210*tmp_210))) + 0.0068572537431980923*tmp_66*(-tmp_217*tmp_50 - tmp_218*tmp_55 - tmp_219*tmp_60 + tmp_65*((tmp_217*tmp_217) + (tmp_218*tmp_218) + (tmp_219*tmp_219))) + 0.037198804536718075*tmp_66*(-tmp_226*tmp_50 - tmp_227*tmp_55 - tmp_228*tmp_60 + tmp_65*((tmp_226*tmp_226) + (tmp_227*tmp_227) + (tmp_228*tmp_228))) + 0.042507265838595799*tmp_66*(-tmp_235*tmp_50 - tmp_236*tmp_55 - tmp_237*tmp_60 + tmp_65*((tmp_235*tmp_235) + (tmp_236*tmp_236) + (tmp_237*tmp_237))) + 0.019202922745021479*tmp_66*(-tmp_244*tmp_50 - tmp_245*tmp_55 - tmp_246*tmp_60 + tmp_65*((tmp_244*tmp_244) + (tmp_245*tmp_245) + (tmp_246*tmp_246))) + 0.0068572537431980923*tmp_66*(-tmp_46*tmp_50 - tmp_51*tmp_55 - tmp_56*tmp_60 + tmp_65*((tmp_46*tmp_46) + (tmp_51*tmp_51) + (tmp_56*tmp_56))) + 0.037198804536718075*tmp_66*(-tmp_50*tmp_73 - tmp_55*tmp_74 - tmp_60*tmp_75 + tmp_65*((tmp_73*tmp_73) + (tmp_74*tmp_74) + (tmp_75*tmp_75))) + 0.020848748529055869*tmp_66*(-tmp_50*tmp_82 - tmp_55*tmp_83 - tmp_60*tmp_84 + tmp_65*((tmp_82*tmp_82) + (tmp_83*tmp_83) + (tmp_84*tmp_84))) + 0.019202922745021479*tmp_66*(-tmp_50*tmp_91 - tmp_55*tmp_92 - tmp_60*tmp_93 + tmp_65*((tmp_91*tmp_91) + (tmp_92*tmp_92) + (tmp_93*tmp_93)));
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


      real_t tmp_0 = -p_affine_4_0;
      real_t tmp_1 = p_affine_5_0 + tmp_0;
      real_t tmp_2 = p_affine_6_0 + tmp_0;
      real_t tmp_3 = -p_affine_4_1;
      real_t tmp_4 = p_affine_7_1 + tmp_3;
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = p_affine_7_0 + tmp_0;
      real_t tmp_7 = p_affine_6_1 + tmp_3;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = tmp_5 - tmp_8;
      real_t tmp_10 = -p_affine_4_2;
      real_t tmp_11 = p_affine_7_2 + tmp_10;
      real_t tmp_12 = tmp_11*tmp_7;
      real_t tmp_13 = p_affine_5_2 + tmp_10;
      real_t tmp_14 = p_affine_5_1 + tmp_3;
      real_t tmp_15 = p_affine_6_2 + tmp_10;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_15*tmp_4;
      real_t tmp_18 = tmp_11*tmp_2;
      real_t tmp_19 = 1.0 / (tmp_1*tmp_12 - tmp_1*tmp_17 + tmp_13*tmp_5 - tmp_13*tmp_8 + tmp_14*tmp_16 - tmp_14*tmp_18);
      real_t tmp_20 = p_affine_8_2 + tmp_10;
      real_t tmp_21 = -p_affine_8_2;
      real_t tmp_22 = p_affine_9_2 + tmp_21;
      real_t tmp_23 = p_affine_10_2 + tmp_21;
      real_t tmp_24 = 0.031405749086161582*tmp_22 + 0.93718850182767688*tmp_23;
      real_t tmp_25 = tmp_19*(tmp_20 + tmp_24);
      real_t tmp_26 = tmp_16 - tmp_18;
      real_t tmp_27 = p_affine_8_1 + tmp_3;
      real_t tmp_28 = -p_affine_8_1;
      real_t tmp_29 = p_affine_9_1 + tmp_28;
      real_t tmp_30 = p_affine_10_1 + tmp_28;
      real_t tmp_31 = 0.031405749086161582*tmp_29 + 0.93718850182767688*tmp_30;
      real_t tmp_32 = tmp_19*(tmp_27 + tmp_31);
      real_t tmp_33 = tmp_12 - tmp_17;
      real_t tmp_34 = p_affine_8_0 + tmp_0;
      real_t tmp_35 = -p_affine_8_0;
      real_t tmp_36 = p_affine_9_0 + tmp_35;
      real_t tmp_37 = p_affine_10_0 + tmp_35;
      real_t tmp_38 = 0.031405749086161582*tmp_36 + 0.93718850182767688*tmp_37;
      real_t tmp_39 = tmp_19*(tmp_34 + tmp_38);
      real_t tmp_40 = tmp_25*tmp_9 + tmp_26*tmp_32 + tmp_33*tmp_39 - 1.0/4.0;
      real_t tmp_41 = -tmp_1*tmp_4 + tmp_14*tmp_6;
      real_t tmp_42 = tmp_1*tmp_11 - tmp_13*tmp_6;
      real_t tmp_43 = -tmp_11*tmp_14 + tmp_13*tmp_4;
      real_t tmp_44 = tmp_25*tmp_41 + tmp_32*tmp_42 + tmp_39*tmp_43 - 1.0/4.0;
      real_t tmp_45 = tmp_1*tmp_7 - tmp_14*tmp_2;
      real_t tmp_46 = -tmp_1*tmp_15 + tmp_13*tmp_2;
      real_t tmp_47 = -tmp_13*tmp_7 + tmp_14*tmp_15;
      real_t tmp_48 = tmp_25*tmp_45 + tmp_32*tmp_46 + tmp_39*tmp_47 - 1.0/4.0;
      real_t tmp_49 = tmp_1*tmp_40 + tmp_2*tmp_44 + tmp_48*tmp_6;
      real_t tmp_50 = -p_affine_0_1;
      real_t tmp_51 = p_affine_2_1 + tmp_50;
      real_t tmp_52 = -p_affine_0_2;
      real_t tmp_53 = p_affine_3_2 + tmp_52;
      real_t tmp_54 = tmp_51*tmp_53;
      real_t tmp_55 = p_affine_3_1 + tmp_50;
      real_t tmp_56 = p_affine_2_2 + tmp_52;
      real_t tmp_57 = tmp_55*tmp_56;
      real_t tmp_58 = tmp_54 - tmp_57;
      real_t tmp_59 = -p_affine_0_0;
      real_t tmp_60 = p_affine_1_0 + tmp_59;
      real_t tmp_61 = p_affine_2_0 + tmp_59;
      real_t tmp_62 = p_affine_1_2 + tmp_52;
      real_t tmp_63 = tmp_55*tmp_62;
      real_t tmp_64 = p_affine_3_0 + tmp_59;
      real_t tmp_65 = p_affine_1_1 + tmp_50;
      real_t tmp_66 = tmp_56*tmp_65;
      real_t tmp_67 = tmp_53*tmp_65;
      real_t tmp_68 = tmp_51*tmp_62;
      real_t tmp_69 = 1.0 / (tmp_54*tmp_60 - tmp_57*tmp_60 + tmp_61*tmp_63 - tmp_61*tmp_67 + tmp_64*tmp_66 - tmp_64*tmp_68);
      real_t tmp_70 = tmp_60*tmp_69;
      real_t tmp_71 = tmp_63 - tmp_67;
      real_t tmp_72 = tmp_61*tmp_69;
      real_t tmp_73 = tmp_66 - tmp_68;
      real_t tmp_74 = tmp_64*tmp_69;
      real_t tmp_75 = -tmp_53*tmp_61 + tmp_56*tmp_64;
      real_t tmp_76 = tmp_53*tmp_60 - tmp_62*tmp_64;
      real_t tmp_77 = -tmp_56*tmp_60 + tmp_61*tmp_62;
      real_t tmp_78 = -tmp_51*tmp_64 + tmp_55*tmp_61;
      real_t tmp_79 = -tmp_55*tmp_60 + tmp_64*tmp_65;
      real_t tmp_80 = tmp_51*tmp_60 - tmp_61*tmp_65;
      real_t tmp_81 = 0.5*p_affine_13_0*(tmp_58*tmp_70 + tmp_71*tmp_72 + tmp_73*tmp_74) + 0.5*p_affine_13_1*(tmp_70*tmp_75 + tmp_72*tmp_76 + tmp_74*tmp_77) + 0.5*p_affine_13_2*(tmp_70*tmp_78 + tmp_72*tmp_79 + tmp_74*tmp_80);
      real_t tmp_82 = tmp_14*tmp_40 + tmp_4*tmp_48 + tmp_44*tmp_7;
      real_t tmp_83 = tmp_65*tmp_69;
      real_t tmp_84 = tmp_51*tmp_69;
      real_t tmp_85 = tmp_55*tmp_69;
      real_t tmp_86 = 0.5*p_affine_13_0*(tmp_58*tmp_83 + tmp_71*tmp_84 + tmp_73*tmp_85) + 0.5*p_affine_13_1*(tmp_75*tmp_83 + tmp_76*tmp_84 + tmp_77*tmp_85) + 0.5*p_affine_13_2*(tmp_78*tmp_83 + tmp_79*tmp_84 + tmp_80*tmp_85);
      real_t tmp_87 = tmp_11*tmp_48 + tmp_13*tmp_40 + tmp_15*tmp_44;
      real_t tmp_88 = tmp_62*tmp_69;
      real_t tmp_89 = tmp_56*tmp_69;
      real_t tmp_90 = tmp_53*tmp_69;
      real_t tmp_91 = 0.5*p_affine_13_0*(tmp_58*tmp_88 + tmp_71*tmp_89 + tmp_73*tmp_90) + 0.5*p_affine_13_1*(tmp_75*tmp_88 + tmp_76*tmp_89 + tmp_77*tmp_90) + 0.5*p_affine_13_2*(tmp_78*tmp_88 + tmp_79*tmp_89 + tmp_80*tmp_90);
      real_t tmp_92 = p_affine_8_2 + tmp_52;
      real_t tmp_93 = tmp_69*(tmp_24 + tmp_92);
      real_t tmp_94 = p_affine_8_1 + tmp_50;
      real_t tmp_95 = tmp_69*(tmp_31 + tmp_94);
      real_t tmp_96 = p_affine_8_0 + tmp_59;
      real_t tmp_97 = tmp_69*(tmp_38 + tmp_96);
      real_t tmp_98 = tmp_58*tmp_97 + tmp_75*tmp_95 + tmp_78*tmp_93 - 1.0/4.0;
      real_t tmp_99 = tmp_71*tmp_97 + tmp_76*tmp_95 + tmp_79*tmp_93 - 1.0/4.0;
      real_t tmp_100 = tmp_73*tmp_97 + tmp_77*tmp_95 + tmp_80*tmp_93 - 1.0/4.0;
      real_t tmp_101 = tmp_100*tmp_64 + tmp_60*tmp_98 + tmp_61*tmp_99;
      real_t tmp_102 = tmp_1*tmp_19;
      real_t tmp_103 = tmp_19*tmp_2;
      real_t tmp_104 = tmp_19*tmp_6;
      real_t tmp_105 = 0.5*p_affine_13_0*(tmp_102*tmp_33 + tmp_103*tmp_43 + tmp_104*tmp_47) + 0.5*p_affine_13_1*(tmp_102*tmp_26 + tmp_103*tmp_42 + tmp_104*tmp_46) + 0.5*p_affine_13_2*(tmp_102*tmp_9 + tmp_103*tmp_41 + tmp_104*tmp_45);
      real_t tmp_106 = tmp_100*tmp_55 + tmp_51*tmp_99 + tmp_65*tmp_98;
      real_t tmp_107 = tmp_14*tmp_19;
      real_t tmp_108 = tmp_19*tmp_7;
      real_t tmp_109 = tmp_19*tmp_4;
      real_t tmp_110 = 0.5*p_affine_13_0*(tmp_107*tmp_33 + tmp_108*tmp_43 + tmp_109*tmp_47) + 0.5*p_affine_13_1*(tmp_107*tmp_26 + tmp_108*tmp_42 + tmp_109*tmp_46) + 0.5*p_affine_13_2*(tmp_107*tmp_9 + tmp_108*tmp_41 + tmp_109*tmp_45);
      real_t tmp_111 = tmp_100*tmp_53 + tmp_56*tmp_99 + tmp_62*tmp_98;
      real_t tmp_112 = tmp_13*tmp_19;
      real_t tmp_113 = tmp_15*tmp_19;
      real_t tmp_114 = tmp_11*tmp_19;
      real_t tmp_115 = 0.5*p_affine_13_0*(tmp_112*tmp_33 + tmp_113*tmp_43 + tmp_114*tmp_47) + 0.5*p_affine_13_1*(tmp_112*tmp_26 + tmp_113*tmp_42 + tmp_114*tmp_46) + 0.5*p_affine_13_2*(tmp_112*tmp_9 + tmp_113*tmp_41 + tmp_114*tmp_45);
      real_t tmp_116 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_117 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_118 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_119 = (std::abs(tmp_116*tmp_23 - tmp_118*tmp_30)*std::abs(tmp_116*tmp_23 - tmp_118*tmp_30)) + (std::abs(tmp_116*tmp_37 - tmp_117*tmp_30)*std::abs(tmp_116*tmp_37 - tmp_117*tmp_30)) + (std::abs(tmp_117*tmp_23 - tmp_118*tmp_37)*std::abs(tmp_117*tmp_23 - tmp_118*tmp_37));
      real_t tmp_120 = 5.0*std::pow(tmp_119, -0.25);
      real_t tmp_121 = 1.0*std::pow(tmp_119, 1.0/2.0);
      real_t tmp_122 = 0.19601935860219369*tmp_22 + 0.60796128279561268*tmp_23;
      real_t tmp_123 = tmp_19*(tmp_122 + tmp_20);
      real_t tmp_124 = 0.19601935860219369*tmp_29 + 0.60796128279561268*tmp_30;
      real_t tmp_125 = tmp_19*(tmp_124 + tmp_27);
      real_t tmp_126 = 0.19601935860219369*tmp_36 + 0.60796128279561268*tmp_37;
      real_t tmp_127 = tmp_19*(tmp_126 + tmp_34);
      real_t tmp_128 = tmp_123*tmp_9 + tmp_125*tmp_26 + tmp_127*tmp_33 - 1.0/4.0;
      real_t tmp_129 = tmp_123*tmp_41 + tmp_125*tmp_42 + tmp_127*tmp_43 - 1.0/4.0;
      real_t tmp_130 = tmp_123*tmp_45 + tmp_125*tmp_46 + tmp_127*tmp_47 - 1.0/4.0;
      real_t tmp_131 = tmp_1*tmp_128 + tmp_129*tmp_2 + tmp_130*tmp_6;
      real_t tmp_132 = tmp_128*tmp_14 + tmp_129*tmp_7 + tmp_130*tmp_4;
      real_t tmp_133 = tmp_11*tmp_130 + tmp_128*tmp_13 + tmp_129*tmp_15;
      real_t tmp_134 = tmp_69*(tmp_122 + tmp_92);
      real_t tmp_135 = tmp_69*(tmp_124 + tmp_94);
      real_t tmp_136 = tmp_69*(tmp_126 + tmp_96);
      real_t tmp_137 = tmp_134*tmp_78 + tmp_135*tmp_75 + tmp_136*tmp_58 - 1.0/4.0;
      real_t tmp_138 = tmp_134*tmp_79 + tmp_135*tmp_76 + tmp_136*tmp_71 - 1.0/4.0;
      real_t tmp_139 = tmp_134*tmp_80 + tmp_135*tmp_77 + tmp_136*tmp_73 - 1.0/4.0;
      real_t tmp_140 = tmp_137*tmp_60 + tmp_138*tmp_61 + tmp_139*tmp_64;
      real_t tmp_141 = tmp_137*tmp_65 + tmp_138*tmp_51 + tmp_139*tmp_55;
      real_t tmp_142 = tmp_137*tmp_62 + tmp_138*tmp_56 + tmp_139*tmp_53;
      real_t tmp_143 = 0.37605877282253791*tmp_22 + 0.039308471900058539*tmp_23;
      real_t tmp_144 = tmp_19*(tmp_143 + tmp_20);
      real_t tmp_145 = 0.37605877282253791*tmp_29 + 0.039308471900058539*tmp_30;
      real_t tmp_146 = tmp_19*(tmp_145 + tmp_27);
      real_t tmp_147 = 0.37605877282253791*tmp_36 + 0.039308471900058539*tmp_37;
      real_t tmp_148 = tmp_19*(tmp_147 + tmp_34);
      real_t tmp_149 = tmp_144*tmp_9 + tmp_146*tmp_26 + tmp_148*tmp_33 - 1.0/4.0;
      real_t tmp_150 = tmp_144*tmp_41 + tmp_146*tmp_42 + tmp_148*tmp_43 - 1.0/4.0;
      real_t tmp_151 = tmp_144*tmp_45 + tmp_146*tmp_46 + tmp_148*tmp_47 - 1.0/4.0;
      real_t tmp_152 = tmp_1*tmp_149 + tmp_150*tmp_2 + tmp_151*tmp_6;
      real_t tmp_153 = tmp_14*tmp_149 + tmp_150*tmp_7 + tmp_151*tmp_4;
      real_t tmp_154 = tmp_11*tmp_151 + tmp_13*tmp_149 + tmp_15*tmp_150;
      real_t tmp_155 = tmp_69*(tmp_143 + tmp_92);
      real_t tmp_156 = tmp_69*(tmp_145 + tmp_94);
      real_t tmp_157 = tmp_69*(tmp_147 + tmp_96);
      real_t tmp_158 = tmp_155*tmp_78 + tmp_156*tmp_75 + tmp_157*tmp_58 - 1.0/4.0;
      real_t tmp_159 = tmp_155*tmp_79 + tmp_156*tmp_76 + tmp_157*tmp_71 - 1.0/4.0;
      real_t tmp_160 = tmp_155*tmp_80 + tmp_156*tmp_77 + tmp_157*tmp_73 - 1.0/4.0;
      real_t tmp_161 = tmp_158*tmp_60 + tmp_159*tmp_61 + tmp_160*tmp_64;
      real_t tmp_162 = tmp_158*tmp_65 + tmp_159*tmp_51 + tmp_160*tmp_55;
      real_t tmp_163 = tmp_158*tmp_62 + tmp_159*tmp_56 + tmp_160*tmp_53;
      real_t tmp_164 = 0.78764240869137092*tmp_22 + 0.1711304259088916*tmp_23;
      real_t tmp_165 = tmp_19*(tmp_164 + tmp_20);
      real_t tmp_166 = 0.78764240869137092*tmp_29 + 0.1711304259088916*tmp_30;
      real_t tmp_167 = tmp_19*(tmp_166 + tmp_27);
      real_t tmp_168 = 0.78764240869137092*tmp_36 + 0.1711304259088916*tmp_37;
      real_t tmp_169 = tmp_19*(tmp_168 + tmp_34);
      real_t tmp_170 = tmp_165*tmp_9 + tmp_167*tmp_26 + tmp_169*tmp_33 - 1.0/4.0;
      real_t tmp_171 = tmp_165*tmp_41 + tmp_167*tmp_42 + tmp_169*tmp_43 - 1.0/4.0;
      real_t tmp_172 = tmp_165*tmp_45 + tmp_167*tmp_46 + tmp_169*tmp_47 - 1.0/4.0;
      real_t tmp_173 = tmp_1*tmp_170 + tmp_171*tmp_2 + tmp_172*tmp_6;
      real_t tmp_174 = tmp_14*tmp_170 + tmp_171*tmp_7 + tmp_172*tmp_4;
      real_t tmp_175 = tmp_11*tmp_172 + tmp_13*tmp_170 + tmp_15*tmp_171;
      real_t tmp_176 = tmp_69*(tmp_164 + tmp_92);
      real_t tmp_177 = tmp_69*(tmp_166 + tmp_94);
      real_t tmp_178 = tmp_69*(tmp_168 + tmp_96);
      real_t tmp_179 = tmp_176*tmp_78 + tmp_177*tmp_75 + tmp_178*tmp_58 - 1.0/4.0;
      real_t tmp_180 = tmp_176*tmp_79 + tmp_177*tmp_76 + tmp_178*tmp_71 - 1.0/4.0;
      real_t tmp_181 = tmp_176*tmp_80 + tmp_177*tmp_77 + tmp_178*tmp_73 - 1.0/4.0;
      real_t tmp_182 = tmp_179*tmp_60 + tmp_180*tmp_61 + tmp_181*tmp_64;
      real_t tmp_183 = tmp_179*tmp_65 + tmp_180*tmp_51 + tmp_181*tmp_55;
      real_t tmp_184 = tmp_179*tmp_62 + tmp_180*tmp_56 + tmp_181*tmp_53;
      real_t tmp_185 = 0.58463275527740355*tmp_22 + 0.37605877282253791*tmp_23;
      real_t tmp_186 = tmp_19*(tmp_185 + tmp_20);
      real_t tmp_187 = 0.58463275527740355*tmp_29 + 0.37605877282253791*tmp_30;
      real_t tmp_188 = tmp_19*(tmp_187 + tmp_27);
      real_t tmp_189 = 0.58463275527740355*tmp_36 + 0.37605877282253791*tmp_37;
      real_t tmp_190 = tmp_19*(tmp_189 + tmp_34);
      real_t tmp_191 = tmp_186*tmp_9 + tmp_188*tmp_26 + tmp_190*tmp_33 - 1.0/4.0;
      real_t tmp_192 = tmp_186*tmp_41 + tmp_188*tmp_42 + tmp_190*tmp_43 - 1.0/4.0;
      real_t tmp_193 = tmp_186*tmp_45 + tmp_188*tmp_46 + tmp_190*tmp_47 - 1.0/4.0;
      real_t tmp_194 = tmp_1*tmp_191 + tmp_192*tmp_2 + tmp_193*tmp_6;
      real_t tmp_195 = tmp_14*tmp_191 + tmp_192*tmp_7 + tmp_193*tmp_4;
      real_t tmp_196 = tmp_11*tmp_193 + tmp_13*tmp_191 + tmp_15*tmp_192;
      real_t tmp_197 = tmp_69*(tmp_185 + tmp_92);
      real_t tmp_198 = tmp_69*(tmp_187 + tmp_94);
      real_t tmp_199 = tmp_69*(tmp_189 + tmp_96);
      real_t tmp_200 = tmp_197*tmp_78 + tmp_198*tmp_75 + tmp_199*tmp_58 - 1.0/4.0;
      real_t tmp_201 = tmp_197*tmp_79 + tmp_198*tmp_76 + tmp_199*tmp_71 - 1.0/4.0;
      real_t tmp_202 = tmp_197*tmp_80 + tmp_198*tmp_77 + tmp_199*tmp_73 - 1.0/4.0;
      real_t tmp_203 = tmp_200*tmp_60 + tmp_201*tmp_61 + tmp_202*tmp_64;
      real_t tmp_204 = tmp_200*tmp_65 + tmp_201*tmp_51 + tmp_202*tmp_55;
      real_t tmp_205 = tmp_200*tmp_62 + tmp_201*tmp_56 + tmp_202*tmp_53;
      real_t tmp_206 = 0.041227165399737475*tmp_22 + 0.78764240869137092*tmp_23;
      real_t tmp_207 = tmp_19*(tmp_20 + tmp_206);
      real_t tmp_208 = 0.041227165399737475*tmp_29 + 0.78764240869137092*tmp_30;
      real_t tmp_209 = tmp_19*(tmp_208 + tmp_27);
      real_t tmp_210 = 0.041227165399737475*tmp_36 + 0.78764240869137092*tmp_37;
      real_t tmp_211 = tmp_19*(tmp_210 + tmp_34);
      real_t tmp_212 = tmp_207*tmp_9 + tmp_209*tmp_26 + tmp_211*tmp_33 - 1.0/4.0;
      real_t tmp_213 = tmp_207*tmp_41 + tmp_209*tmp_42 + tmp_211*tmp_43 - 1.0/4.0;
      real_t tmp_214 = tmp_207*tmp_45 + tmp_209*tmp_46 + tmp_211*tmp_47 - 1.0/4.0;
      real_t tmp_215 = tmp_1*tmp_212 + tmp_2*tmp_213 + tmp_214*tmp_6;
      real_t tmp_216 = tmp_14*tmp_212 + tmp_213*tmp_7 + tmp_214*tmp_4;
      real_t tmp_217 = tmp_11*tmp_214 + tmp_13*tmp_212 + tmp_15*tmp_213;
      real_t tmp_218 = tmp_69*(tmp_206 + tmp_92);
      real_t tmp_219 = tmp_69*(tmp_208 + tmp_94);
      real_t tmp_220 = tmp_69*(tmp_210 + tmp_96);
      real_t tmp_221 = tmp_218*tmp_78 + tmp_219*tmp_75 + tmp_220*tmp_58 - 1.0/4.0;
      real_t tmp_222 = tmp_218*tmp_79 + tmp_219*tmp_76 + tmp_220*tmp_71 - 1.0/4.0;
      real_t tmp_223 = tmp_218*tmp_80 + tmp_219*tmp_77 + tmp_220*tmp_73 - 1.0/4.0;
      real_t tmp_224 = tmp_221*tmp_60 + tmp_222*tmp_61 + tmp_223*tmp_64;
      real_t tmp_225 = tmp_221*tmp_65 + tmp_222*tmp_51 + tmp_223*tmp_55;
      real_t tmp_226 = tmp_221*tmp_62 + tmp_222*tmp_56 + tmp_223*tmp_53;
      real_t tmp_227 = 0.039308471900058539*tmp_22 + 0.58463275527740355*tmp_23;
      real_t tmp_228 = tmp_19*(tmp_20 + tmp_227);
      real_t tmp_229 = 0.039308471900058539*tmp_29 + 0.58463275527740355*tmp_30;
      real_t tmp_230 = tmp_19*(tmp_229 + tmp_27);
      real_t tmp_231 = 0.039308471900058539*tmp_36 + 0.58463275527740355*tmp_37;
      real_t tmp_232 = tmp_19*(tmp_231 + tmp_34);
      real_t tmp_233 = tmp_228*tmp_9 + tmp_230*tmp_26 + tmp_232*tmp_33 - 1.0/4.0;
      real_t tmp_234 = tmp_228*tmp_41 + tmp_230*tmp_42 + tmp_232*tmp_43 - 1.0/4.0;
      real_t tmp_235 = tmp_228*tmp_45 + tmp_230*tmp_46 + tmp_232*tmp_47 - 1.0/4.0;
      real_t tmp_236 = tmp_1*tmp_233 + tmp_2*tmp_234 + tmp_235*tmp_6;
      real_t tmp_237 = tmp_14*tmp_233 + tmp_234*tmp_7 + tmp_235*tmp_4;
      real_t tmp_238 = tmp_11*tmp_235 + tmp_13*tmp_233 + tmp_15*tmp_234;
      real_t tmp_239 = tmp_69*(tmp_227 + tmp_92);
      real_t tmp_240 = tmp_69*(tmp_229 + tmp_94);
      real_t tmp_241 = tmp_69*(tmp_231 + tmp_96);
      real_t tmp_242 = tmp_239*tmp_78 + tmp_240*tmp_75 + tmp_241*tmp_58 - 1.0/4.0;
      real_t tmp_243 = tmp_239*tmp_79 + tmp_240*tmp_76 + tmp_241*tmp_71 - 1.0/4.0;
      real_t tmp_244 = tmp_239*tmp_80 + tmp_240*tmp_77 + tmp_241*tmp_73 - 1.0/4.0;
      real_t tmp_245 = tmp_242*tmp_60 + tmp_243*tmp_61 + tmp_244*tmp_64;
      real_t tmp_246 = tmp_242*tmp_65 + tmp_243*tmp_51 + tmp_244*tmp_55;
      real_t tmp_247 = tmp_242*tmp_62 + tmp_243*tmp_56 + tmp_244*tmp_53;
      real_t tmp_248 = 0.78764240869137092*tmp_22 + 0.041227165399737475*tmp_23;
      real_t tmp_249 = tmp_19*(tmp_20 + tmp_248);
      real_t tmp_250 = 0.78764240869137092*tmp_29 + 0.041227165399737475*tmp_30;
      real_t tmp_251 = tmp_19*(tmp_250 + tmp_27);
      real_t tmp_252 = 0.78764240869137092*tmp_36 + 0.041227165399737475*tmp_37;
      real_t tmp_253 = tmp_19*(tmp_252 + tmp_34);
      real_t tmp_254 = tmp_249*tmp_9 + tmp_251*tmp_26 + tmp_253*tmp_33 - 1.0/4.0;
      real_t tmp_255 = tmp_249*tmp_41 + tmp_251*tmp_42 + tmp_253*tmp_43 - 1.0/4.0;
      real_t tmp_256 = tmp_249*tmp_45 + tmp_251*tmp_46 + tmp_253*tmp_47 - 1.0/4.0;
      real_t tmp_257 = tmp_1*tmp_254 + tmp_2*tmp_255 + tmp_256*tmp_6;
      real_t tmp_258 = tmp_14*tmp_254 + tmp_255*tmp_7 + tmp_256*tmp_4;
      real_t tmp_259 = tmp_11*tmp_256 + tmp_13*tmp_254 + tmp_15*tmp_255;
      real_t tmp_260 = tmp_69*(tmp_248 + tmp_92);
      real_t tmp_261 = tmp_69*(tmp_250 + tmp_94);
      real_t tmp_262 = tmp_69*(tmp_252 + tmp_96);
      real_t tmp_263 = tmp_260*tmp_78 + tmp_261*tmp_75 + tmp_262*tmp_58 - 1.0/4.0;
      real_t tmp_264 = tmp_260*tmp_79 + tmp_261*tmp_76 + tmp_262*tmp_71 - 1.0/4.0;
      real_t tmp_265 = tmp_260*tmp_80 + tmp_261*tmp_77 + tmp_262*tmp_73 - 1.0/4.0;
      real_t tmp_266 = tmp_263*tmp_60 + tmp_264*tmp_61 + tmp_265*tmp_64;
      real_t tmp_267 = tmp_263*tmp_65 + tmp_264*tmp_51 + tmp_265*tmp_55;
      real_t tmp_268 = tmp_263*tmp_62 + tmp_264*tmp_56 + tmp_265*tmp_53;
      real_t tmp_269 = 0.58463275527740355*tmp_22 + 0.039308471900058539*tmp_23;
      real_t tmp_270 = tmp_19*(tmp_20 + tmp_269);
      real_t tmp_271 = 0.58463275527740355*tmp_29 + 0.039308471900058539*tmp_30;
      real_t tmp_272 = tmp_19*(tmp_27 + tmp_271);
      real_t tmp_273 = 0.58463275527740355*tmp_36 + 0.039308471900058539*tmp_37;
      real_t tmp_274 = tmp_19*(tmp_273 + tmp_34);
      real_t tmp_275 = tmp_26*tmp_272 + tmp_270*tmp_9 + tmp_274*tmp_33 - 1.0/4.0;
      real_t tmp_276 = tmp_270*tmp_41 + tmp_272*tmp_42 + tmp_274*tmp_43 - 1.0/4.0;
      real_t tmp_277 = tmp_270*tmp_45 + tmp_272*tmp_46 + tmp_274*tmp_47 - 1.0/4.0;
      real_t tmp_278 = tmp_1*tmp_275 + tmp_2*tmp_276 + tmp_277*tmp_6;
      real_t tmp_279 = tmp_14*tmp_275 + tmp_276*tmp_7 + tmp_277*tmp_4;
      real_t tmp_280 = tmp_11*tmp_277 + tmp_13*tmp_275 + tmp_15*tmp_276;
      real_t tmp_281 = tmp_69*(tmp_269 + tmp_92);
      real_t tmp_282 = tmp_69*(tmp_271 + tmp_94);
      real_t tmp_283 = tmp_69*(tmp_273 + tmp_96);
      real_t tmp_284 = tmp_281*tmp_78 + tmp_282*tmp_75 + tmp_283*tmp_58 - 1.0/4.0;
      real_t tmp_285 = tmp_281*tmp_79 + tmp_282*tmp_76 + tmp_283*tmp_71 - 1.0/4.0;
      real_t tmp_286 = tmp_281*tmp_80 + tmp_282*tmp_77 + tmp_283*tmp_73 - 1.0/4.0;
      real_t tmp_287 = tmp_284*tmp_60 + tmp_285*tmp_61 + tmp_286*tmp_64;
      real_t tmp_288 = tmp_284*tmp_65 + tmp_285*tmp_51 + tmp_286*tmp_55;
      real_t tmp_289 = tmp_284*tmp_62 + tmp_285*tmp_56 + tmp_286*tmp_53;
      real_t tmp_290 = 0.1711304259088916*tmp_22 + 0.78764240869137092*tmp_23;
      real_t tmp_291 = tmp_19*(tmp_20 + tmp_290);
      real_t tmp_292 = 0.1711304259088916*tmp_29 + 0.78764240869137092*tmp_30;
      real_t tmp_293 = tmp_19*(tmp_27 + tmp_292);
      real_t tmp_294 = 0.1711304259088916*tmp_36 + 0.78764240869137092*tmp_37;
      real_t tmp_295 = tmp_19*(tmp_294 + tmp_34);
      real_t tmp_296 = tmp_26*tmp_293 + tmp_291*tmp_9 + tmp_295*tmp_33 - 1.0/4.0;
      real_t tmp_297 = tmp_291*tmp_41 + tmp_293*tmp_42 + tmp_295*tmp_43 - 1.0/4.0;
      real_t tmp_298 = tmp_291*tmp_45 + tmp_293*tmp_46 + tmp_295*tmp_47 - 1.0/4.0;
      real_t tmp_299 = tmp_1*tmp_296 + tmp_2*tmp_297 + tmp_298*tmp_6;
      real_t tmp_300 = tmp_14*tmp_296 + tmp_297*tmp_7 + tmp_298*tmp_4;
      real_t tmp_301 = tmp_11*tmp_298 + tmp_13*tmp_296 + tmp_15*tmp_297;
      real_t tmp_302 = tmp_69*(tmp_290 + tmp_92);
      real_t tmp_303 = tmp_69*(tmp_292 + tmp_94);
      real_t tmp_304 = tmp_69*(tmp_294 + tmp_96);
      real_t tmp_305 = tmp_302*tmp_78 + tmp_303*tmp_75 + tmp_304*tmp_58 - 1.0/4.0;
      real_t tmp_306 = tmp_302*tmp_79 + tmp_303*tmp_76 + tmp_304*tmp_71 - 1.0/4.0;
      real_t tmp_307 = tmp_302*tmp_80 + tmp_303*tmp_77 + tmp_304*tmp_73 - 1.0/4.0;
      real_t tmp_308 = tmp_305*tmp_60 + tmp_306*tmp_61 + tmp_307*tmp_64;
      real_t tmp_309 = tmp_305*tmp_65 + tmp_306*tmp_51 + tmp_307*tmp_55;
      real_t tmp_310 = tmp_305*tmp_62 + tmp_306*tmp_56 + tmp_307*tmp_53;
      real_t tmp_311 = 0.37605877282253791*tmp_22 + 0.58463275527740355*tmp_23;
      real_t tmp_312 = tmp_19*(tmp_20 + tmp_311);
      real_t tmp_313 = 0.37605877282253791*tmp_29 + 0.58463275527740355*tmp_30;
      real_t tmp_314 = tmp_19*(tmp_27 + tmp_313);
      real_t tmp_315 = 0.37605877282253791*tmp_36 + 0.58463275527740355*tmp_37;
      real_t tmp_316 = tmp_19*(tmp_315 + tmp_34);
      real_t tmp_317 = tmp_26*tmp_314 + tmp_312*tmp_9 + tmp_316*tmp_33 - 1.0/4.0;
      real_t tmp_318 = tmp_312*tmp_41 + tmp_314*tmp_42 + tmp_316*tmp_43 - 1.0/4.0;
      real_t tmp_319 = tmp_312*tmp_45 + tmp_314*tmp_46 + tmp_316*tmp_47 - 1.0/4.0;
      real_t tmp_320 = tmp_1*tmp_317 + tmp_2*tmp_318 + tmp_319*tmp_6;
      real_t tmp_321 = tmp_14*tmp_317 + tmp_318*tmp_7 + tmp_319*tmp_4;
      real_t tmp_322 = tmp_11*tmp_319 + tmp_13*tmp_317 + tmp_15*tmp_318;
      real_t tmp_323 = tmp_69*(tmp_311 + tmp_92);
      real_t tmp_324 = tmp_69*(tmp_313 + tmp_94);
      real_t tmp_325 = tmp_69*(tmp_315 + tmp_96);
      real_t tmp_326 = tmp_323*tmp_78 + tmp_324*tmp_75 + tmp_325*tmp_58 - 1.0/4.0;
      real_t tmp_327 = tmp_323*tmp_79 + tmp_324*tmp_76 + tmp_325*tmp_71 - 1.0/4.0;
      real_t tmp_328 = tmp_323*tmp_80 + tmp_324*tmp_77 + tmp_325*tmp_73 - 1.0/4.0;
      real_t tmp_329 = tmp_326*tmp_60 + tmp_327*tmp_61 + tmp_328*tmp_64;
      real_t tmp_330 = tmp_326*tmp_65 + tmp_327*tmp_51 + tmp_328*tmp_55;
      real_t tmp_331 = tmp_326*tmp_62 + tmp_327*tmp_56 + tmp_328*tmp_53;
      real_t tmp_332 = 0.041227165399737475*tmp_22 + 0.1711304259088916*tmp_23;
      real_t tmp_333 = tmp_19*(tmp_20 + tmp_332);
      real_t tmp_334 = 0.041227165399737475*tmp_29 + 0.1711304259088916*tmp_30;
      real_t tmp_335 = tmp_19*(tmp_27 + tmp_334);
      real_t tmp_336 = 0.041227165399737475*tmp_36 + 0.1711304259088916*tmp_37;
      real_t tmp_337 = tmp_19*(tmp_336 + tmp_34);
      real_t tmp_338 = tmp_26*tmp_335 + tmp_33*tmp_337 + tmp_333*tmp_9 - 1.0/4.0;
      real_t tmp_339 = tmp_333*tmp_41 + tmp_335*tmp_42 + tmp_337*tmp_43 - 1.0/4.0;
      real_t tmp_340 = tmp_333*tmp_45 + tmp_335*tmp_46 + tmp_337*tmp_47 - 1.0/4.0;
      real_t tmp_341 = tmp_1*tmp_338 + tmp_2*tmp_339 + tmp_340*tmp_6;
      real_t tmp_342 = tmp_14*tmp_338 + tmp_339*tmp_7 + tmp_340*tmp_4;
      real_t tmp_343 = tmp_11*tmp_340 + tmp_13*tmp_338 + tmp_15*tmp_339;
      real_t tmp_344 = tmp_69*(tmp_332 + tmp_92);
      real_t tmp_345 = tmp_69*(tmp_334 + tmp_94);
      real_t tmp_346 = tmp_69*(tmp_336 + tmp_96);
      real_t tmp_347 = tmp_344*tmp_78 + tmp_345*tmp_75 + tmp_346*tmp_58 - 1.0/4.0;
      real_t tmp_348 = tmp_344*tmp_79 + tmp_345*tmp_76 + tmp_346*tmp_71 - 1.0/4.0;
      real_t tmp_349 = tmp_344*tmp_80 + tmp_345*tmp_77 + tmp_346*tmp_73 - 1.0/4.0;
      real_t tmp_350 = tmp_347*tmp_60 + tmp_348*tmp_61 + tmp_349*tmp_64;
      real_t tmp_351 = tmp_347*tmp_65 + tmp_348*tmp_51 + tmp_349*tmp_55;
      real_t tmp_352 = tmp_347*tmp_62 + tmp_348*tmp_56 + tmp_349*tmp_53;
      real_t tmp_353 = 0.40446199974765351*tmp_22 + 0.19107600050469298*tmp_23;
      real_t tmp_354 = tmp_19*(tmp_20 + tmp_353);
      real_t tmp_355 = 0.40446199974765351*tmp_29 + 0.19107600050469298*tmp_30;
      real_t tmp_356 = tmp_19*(tmp_27 + tmp_355);
      real_t tmp_357 = 0.40446199974765351*tmp_36 + 0.19107600050469298*tmp_37;
      real_t tmp_358 = tmp_19*(tmp_34 + tmp_357);
      real_t tmp_359 = tmp_26*tmp_356 + tmp_33*tmp_358 + tmp_354*tmp_9 - 1.0/4.0;
      real_t tmp_360 = tmp_354*tmp_41 + tmp_356*tmp_42 + tmp_358*tmp_43 - 1.0/4.0;
      real_t tmp_361 = tmp_354*tmp_45 + tmp_356*tmp_46 + tmp_358*tmp_47 - 1.0/4.0;
      real_t tmp_362 = tmp_1*tmp_359 + tmp_2*tmp_360 + tmp_361*tmp_6;
      real_t tmp_363 = tmp_14*tmp_359 + tmp_360*tmp_7 + tmp_361*tmp_4;
      real_t tmp_364 = tmp_11*tmp_361 + tmp_13*tmp_359 + tmp_15*tmp_360;
      real_t tmp_365 = tmp_69*(tmp_353 + tmp_92);
      real_t tmp_366 = tmp_69*(tmp_355 + tmp_94);
      real_t tmp_367 = tmp_69*(tmp_357 + tmp_96);
      real_t tmp_368 = tmp_365*tmp_78 + tmp_366*tmp_75 + tmp_367*tmp_58 - 1.0/4.0;
      real_t tmp_369 = tmp_365*tmp_79 + tmp_366*tmp_76 + tmp_367*tmp_71 - 1.0/4.0;
      real_t tmp_370 = tmp_365*tmp_80 + tmp_366*tmp_77 + tmp_367*tmp_73 - 1.0/4.0;
      real_t tmp_371 = tmp_368*tmp_60 + tmp_369*tmp_61 + tmp_370*tmp_64;
      real_t tmp_372 = tmp_368*tmp_65 + tmp_369*tmp_51 + tmp_370*tmp_55;
      real_t tmp_373 = tmp_368*tmp_62 + tmp_369*tmp_56 + tmp_370*tmp_53;
      real_t tmp_374 = 0.039308471900058539*tmp_22 + 0.37605877282253791*tmp_23;
      real_t tmp_375 = tmp_19*(tmp_20 + tmp_374);
      real_t tmp_376 = 0.039308471900058539*tmp_29 + 0.37605877282253791*tmp_30;
      real_t tmp_377 = tmp_19*(tmp_27 + tmp_376);
      real_t tmp_378 = 0.039308471900058539*tmp_36 + 0.37605877282253791*tmp_37;
      real_t tmp_379 = tmp_19*(tmp_34 + tmp_378);
      real_t tmp_380 = tmp_26*tmp_377 + tmp_33*tmp_379 + tmp_375*tmp_9 - 1.0/4.0;
      real_t tmp_381 = tmp_375*tmp_41 + tmp_377*tmp_42 + tmp_379*tmp_43 - 1.0/4.0;
      real_t tmp_382 = tmp_375*tmp_45 + tmp_377*tmp_46 + tmp_379*tmp_47 - 1.0/4.0;
      real_t tmp_383 = tmp_1*tmp_380 + tmp_2*tmp_381 + tmp_382*tmp_6;
      real_t tmp_384 = tmp_14*tmp_380 + tmp_381*tmp_7 + tmp_382*tmp_4;
      real_t tmp_385 = tmp_11*tmp_382 + tmp_13*tmp_380 + tmp_15*tmp_381;
      real_t tmp_386 = tmp_69*(tmp_374 + tmp_92);
      real_t tmp_387 = tmp_69*(tmp_376 + tmp_94);
      real_t tmp_388 = tmp_69*(tmp_378 + tmp_96);
      real_t tmp_389 = tmp_386*tmp_78 + tmp_387*tmp_75 + tmp_388*tmp_58 - 1.0/4.0;
      real_t tmp_390 = tmp_386*tmp_79 + tmp_387*tmp_76 + tmp_388*tmp_71 - 1.0/4.0;
      real_t tmp_391 = tmp_386*tmp_80 + tmp_387*tmp_77 + tmp_388*tmp_73 - 1.0/4.0;
      real_t tmp_392 = tmp_389*tmp_60 + tmp_390*tmp_61 + tmp_391*tmp_64;
      real_t tmp_393 = tmp_389*tmp_65 + tmp_390*tmp_51 + tmp_391*tmp_55;
      real_t tmp_394 = tmp_389*tmp_62 + tmp_390*tmp_56 + tmp_391*tmp_53;
      real_t tmp_395 = 0.93718850182767688*tmp_22 + 0.031405749086161582*tmp_23;
      real_t tmp_396 = tmp_19*(tmp_20 + tmp_395);
      real_t tmp_397 = 0.93718850182767688*tmp_29 + 0.031405749086161582*tmp_30;
      real_t tmp_398 = tmp_19*(tmp_27 + tmp_397);
      real_t tmp_399 = 0.93718850182767688*tmp_36 + 0.031405749086161582*tmp_37;
      real_t tmp_400 = tmp_19*(tmp_34 + tmp_399);
      real_t tmp_401 = tmp_26*tmp_398 + tmp_33*tmp_400 + tmp_396*tmp_9 - 1.0/4.0;
      real_t tmp_402 = tmp_396*tmp_41 + tmp_398*tmp_42 + tmp_400*tmp_43 - 1.0/4.0;
      real_t tmp_403 = tmp_396*tmp_45 + tmp_398*tmp_46 + tmp_400*tmp_47 - 1.0/4.0;
      real_t tmp_404 = tmp_1*tmp_401 + tmp_2*tmp_402 + tmp_403*tmp_6;
      real_t tmp_405 = tmp_14*tmp_401 + tmp_4*tmp_403 + tmp_402*tmp_7;
      real_t tmp_406 = tmp_11*tmp_403 + tmp_13*tmp_401 + tmp_15*tmp_402;
      real_t tmp_407 = tmp_69*(tmp_395 + tmp_92);
      real_t tmp_408 = tmp_69*(tmp_397 + tmp_94);
      real_t tmp_409 = tmp_69*(tmp_399 + tmp_96);
      real_t tmp_410 = tmp_407*tmp_78 + tmp_408*tmp_75 + tmp_409*tmp_58 - 1.0/4.0;
      real_t tmp_411 = tmp_407*tmp_79 + tmp_408*tmp_76 + tmp_409*tmp_71 - 1.0/4.0;
      real_t tmp_412 = tmp_407*tmp_80 + tmp_408*tmp_77 + tmp_409*tmp_73 - 1.0/4.0;
      real_t tmp_413 = tmp_410*tmp_60 + tmp_411*tmp_61 + tmp_412*tmp_64;
      real_t tmp_414 = tmp_410*tmp_65 + tmp_411*tmp_51 + tmp_412*tmp_55;
      real_t tmp_415 = tmp_410*tmp_62 + tmp_411*tmp_56 + tmp_412*tmp_53;
      real_t tmp_416 = 0.60796128279561268*tmp_22 + 0.19601935860219369*tmp_23;
      real_t tmp_417 = tmp_19*(tmp_20 + tmp_416);
      real_t tmp_418 = 0.60796128279561268*tmp_29 + 0.19601935860219369*tmp_30;
      real_t tmp_419 = tmp_19*(tmp_27 + tmp_418);
      real_t tmp_420 = 0.60796128279561268*tmp_36 + 0.19601935860219369*tmp_37;
      real_t tmp_421 = tmp_19*(tmp_34 + tmp_420);
      real_t tmp_422 = tmp_26*tmp_419 + tmp_33*tmp_421 + tmp_417*tmp_9 - 1.0/4.0;
      real_t tmp_423 = tmp_41*tmp_417 + tmp_419*tmp_42 + tmp_421*tmp_43 - 1.0/4.0;
      real_t tmp_424 = tmp_417*tmp_45 + tmp_419*tmp_46 + tmp_421*tmp_47 - 1.0/4.0;
      real_t tmp_425 = tmp_1*tmp_422 + tmp_2*tmp_423 + tmp_424*tmp_6;
      real_t tmp_426 = tmp_14*tmp_422 + tmp_4*tmp_424 + tmp_423*tmp_7;
      real_t tmp_427 = tmp_11*tmp_424 + tmp_13*tmp_422 + tmp_15*tmp_423;
      real_t tmp_428 = tmp_69*(tmp_416 + tmp_92);
      real_t tmp_429 = tmp_69*(tmp_418 + tmp_94);
      real_t tmp_430 = tmp_69*(tmp_420 + tmp_96);
      real_t tmp_431 = tmp_428*tmp_78 + tmp_429*tmp_75 + tmp_430*tmp_58 - 1.0/4.0;
      real_t tmp_432 = tmp_428*tmp_79 + tmp_429*tmp_76 + tmp_430*tmp_71 - 1.0/4.0;
      real_t tmp_433 = tmp_428*tmp_80 + tmp_429*tmp_77 + tmp_430*tmp_73 - 1.0/4.0;
      real_t tmp_434 = tmp_431*tmp_60 + tmp_432*tmp_61 + tmp_433*tmp_64;
      real_t tmp_435 = tmp_431*tmp_65 + tmp_432*tmp_51 + tmp_433*tmp_55;
      real_t tmp_436 = tmp_431*tmp_62 + tmp_432*tmp_56 + tmp_433*tmp_53;
      real_t tmp_437 = 0.19107600050469298*tmp_22 + 0.40446199974765351*tmp_23;
      real_t tmp_438 = tmp_19*(tmp_20 + tmp_437);
      real_t tmp_439 = 0.19107600050469298*tmp_29 + 0.40446199974765351*tmp_30;
      real_t tmp_440 = tmp_19*(tmp_27 + tmp_439);
      real_t tmp_441 = 0.19107600050469298*tmp_36 + 0.40446199974765351*tmp_37;
      real_t tmp_442 = tmp_19*(tmp_34 + tmp_441);
      real_t tmp_443 = tmp_26*tmp_440 + tmp_33*tmp_442 + tmp_438*tmp_9 - 1.0/4.0;
      real_t tmp_444 = tmp_41*tmp_438 + tmp_42*tmp_440 + tmp_43*tmp_442 - 1.0/4.0;
      real_t tmp_445 = tmp_438*tmp_45 + tmp_440*tmp_46 + tmp_442*tmp_47 - 1.0/4.0;
      real_t tmp_446 = tmp_1*tmp_443 + tmp_2*tmp_444 + tmp_445*tmp_6;
      real_t tmp_447 = tmp_14*tmp_443 + tmp_4*tmp_445 + tmp_444*tmp_7;
      real_t tmp_448 = tmp_11*tmp_445 + tmp_13*tmp_443 + tmp_15*tmp_444;
      real_t tmp_449 = tmp_69*(tmp_437 + tmp_92);
      real_t tmp_450 = tmp_69*(tmp_439 + tmp_94);
      real_t tmp_451 = tmp_69*(tmp_441 + tmp_96);
      real_t tmp_452 = tmp_449*tmp_78 + tmp_450*tmp_75 + tmp_451*tmp_58 - 1.0/4.0;
      real_t tmp_453 = tmp_449*tmp_79 + tmp_450*tmp_76 + tmp_451*tmp_71 - 1.0/4.0;
      real_t tmp_454 = tmp_449*tmp_80 + tmp_450*tmp_77 + tmp_451*tmp_73 - 1.0/4.0;
      real_t tmp_455 = tmp_452*tmp_60 + tmp_453*tmp_61 + tmp_454*tmp_64;
      real_t tmp_456 = tmp_452*tmp_65 + tmp_453*tmp_51 + tmp_454*tmp_55;
      real_t tmp_457 = tmp_452*tmp_62 + tmp_453*tmp_56 + tmp_454*tmp_53;
      real_t tmp_458 = 0.031405749086161582*tmp_22 + 0.031405749086161582*tmp_23;
      real_t tmp_459 = tmp_19*(tmp_20 + tmp_458);
      real_t tmp_460 = 0.031405749086161582*tmp_29 + 0.031405749086161582*tmp_30;
      real_t tmp_461 = tmp_19*(tmp_27 + tmp_460);
      real_t tmp_462 = 0.031405749086161582*tmp_36 + 0.031405749086161582*tmp_37;
      real_t tmp_463 = tmp_19*(tmp_34 + tmp_462);
      real_t tmp_464 = tmp_26*tmp_461 + tmp_33*tmp_463 + tmp_459*tmp_9 - 1.0/4.0;
      real_t tmp_465 = tmp_41*tmp_459 + tmp_42*tmp_461 + tmp_43*tmp_463 - 1.0/4.0;
      real_t tmp_466 = tmp_45*tmp_459 + tmp_46*tmp_461 + tmp_463*tmp_47 - 1.0/4.0;
      real_t tmp_467 = tmp_1*tmp_464 + tmp_2*tmp_465 + tmp_466*tmp_6;
      real_t tmp_468 = tmp_14*tmp_464 + tmp_4*tmp_466 + tmp_465*tmp_7;
      real_t tmp_469 = tmp_11*tmp_466 + tmp_13*tmp_464 + tmp_15*tmp_465;
      real_t tmp_470 = tmp_69*(tmp_458 + tmp_92);
      real_t tmp_471 = tmp_69*(tmp_460 + tmp_94);
      real_t tmp_472 = tmp_69*(tmp_462 + tmp_96);
      real_t tmp_473 = tmp_470*tmp_78 + tmp_471*tmp_75 + tmp_472*tmp_58 - 1.0/4.0;
      real_t tmp_474 = tmp_470*tmp_79 + tmp_471*tmp_76 + tmp_472*tmp_71 - 1.0/4.0;
      real_t tmp_475 = tmp_470*tmp_80 + tmp_471*tmp_77 + tmp_472*tmp_73 - 1.0/4.0;
      real_t tmp_476 = tmp_473*tmp_60 + tmp_474*tmp_61 + tmp_475*tmp_64;
      real_t tmp_477 = tmp_473*tmp_65 + tmp_474*tmp_51 + tmp_475*tmp_55;
      real_t tmp_478 = tmp_473*tmp_62 + tmp_474*tmp_56 + tmp_475*tmp_53;
      real_t tmp_479 = 0.19601935860219369*tmp_22 + 0.19601935860219369*tmp_23;
      real_t tmp_480 = tmp_19*(tmp_20 + tmp_479);
      real_t tmp_481 = 0.19601935860219369*tmp_29 + 0.19601935860219369*tmp_30;
      real_t tmp_482 = tmp_19*(tmp_27 + tmp_481);
      real_t tmp_483 = 0.19601935860219369*tmp_36 + 0.19601935860219369*tmp_37;
      real_t tmp_484 = tmp_19*(tmp_34 + tmp_483);
      real_t tmp_485 = tmp_26*tmp_482 + tmp_33*tmp_484 + tmp_480*tmp_9 - 1.0/4.0;
      real_t tmp_486 = tmp_41*tmp_480 + tmp_42*tmp_482 + tmp_43*tmp_484 - 1.0/4.0;
      real_t tmp_487 = tmp_45*tmp_480 + tmp_46*tmp_482 + tmp_47*tmp_484 - 1.0/4.0;
      real_t tmp_488 = tmp_1*tmp_485 + tmp_2*tmp_486 + tmp_487*tmp_6;
      real_t tmp_489 = tmp_14*tmp_485 + tmp_4*tmp_487 + tmp_486*tmp_7;
      real_t tmp_490 = tmp_11*tmp_487 + tmp_13*tmp_485 + tmp_15*tmp_486;
      real_t tmp_491 = tmp_69*(tmp_479 + tmp_92);
      real_t tmp_492 = tmp_69*(tmp_481 + tmp_94);
      real_t tmp_493 = tmp_69*(tmp_483 + tmp_96);
      real_t tmp_494 = tmp_491*tmp_78 + tmp_492*tmp_75 + tmp_493*tmp_58 - 1.0/4.0;
      real_t tmp_495 = tmp_491*tmp_79 + tmp_492*tmp_76 + tmp_493*tmp_71 - 1.0/4.0;
      real_t tmp_496 = tmp_491*tmp_80 + tmp_492*tmp_77 + tmp_493*tmp_73 - 1.0/4.0;
      real_t tmp_497 = tmp_494*tmp_60 + tmp_495*tmp_61 + tmp_496*tmp_64;
      real_t tmp_498 = tmp_494*tmp_65 + tmp_495*tmp_51 + tmp_496*tmp_55;
      real_t tmp_499 = tmp_494*tmp_62 + tmp_495*tmp_56 + tmp_496*tmp_53;
      real_t tmp_500 = 0.40446199974765351*tmp_22 + 0.40446199974765351*tmp_23;
      real_t tmp_501 = tmp_19*(tmp_20 + tmp_500);
      real_t tmp_502 = 0.40446199974765351*tmp_29 + 0.40446199974765351*tmp_30;
      real_t tmp_503 = tmp_19*(tmp_27 + tmp_502);
      real_t tmp_504 = 0.40446199974765351*tmp_36 + 0.40446199974765351*tmp_37;
      real_t tmp_505 = tmp_19*(tmp_34 + tmp_504);
      real_t tmp_506 = tmp_26*tmp_503 + tmp_33*tmp_505 + tmp_501*tmp_9 - 1.0/4.0;
      real_t tmp_507 = tmp_41*tmp_501 + tmp_42*tmp_503 + tmp_43*tmp_505 - 1.0/4.0;
      real_t tmp_508 = tmp_45*tmp_501 + tmp_46*tmp_503 + tmp_47*tmp_505 - 1.0/4.0;
      real_t tmp_509 = tmp_1*tmp_506 + tmp_2*tmp_507 + tmp_508*tmp_6;
      real_t tmp_510 = tmp_14*tmp_506 + tmp_4*tmp_508 + tmp_507*tmp_7;
      real_t tmp_511 = tmp_11*tmp_508 + tmp_13*tmp_506 + tmp_15*tmp_507;
      real_t tmp_512 = tmp_69*(tmp_500 + tmp_92);
      real_t tmp_513 = tmp_69*(tmp_502 + tmp_94);
      real_t tmp_514 = tmp_69*(tmp_504 + tmp_96);
      real_t tmp_515 = tmp_512*tmp_78 + tmp_513*tmp_75 + tmp_514*tmp_58 - 1.0/4.0;
      real_t tmp_516 = tmp_512*tmp_79 + tmp_513*tmp_76 + tmp_514*tmp_71 - 1.0/4.0;
      real_t tmp_517 = tmp_512*tmp_80 + tmp_513*tmp_77 + tmp_514*tmp_73 - 1.0/4.0;
      real_t tmp_518 = tmp_515*tmp_60 + tmp_516*tmp_61 + tmp_517*tmp_64;
      real_t tmp_519 = tmp_51*tmp_516 + tmp_515*tmp_65 + tmp_517*tmp_55;
      real_t tmp_520 = tmp_515*tmp_62 + tmp_516*tmp_56 + tmp_517*tmp_53;
      real_t tmp_521 = 0.1711304259088916*tmp_22 + 0.041227165399737475*tmp_23;
      real_t tmp_522 = tmp_19*(tmp_20 + tmp_521);
      real_t tmp_523 = 0.1711304259088916*tmp_29 + 0.041227165399737475*tmp_30;
      real_t tmp_524 = tmp_19*(tmp_27 + tmp_523);
      real_t tmp_525 = 0.1711304259088916*tmp_36 + 0.041227165399737475*tmp_37;
      real_t tmp_526 = tmp_19*(tmp_34 + tmp_525);
      real_t tmp_527 = tmp_26*tmp_524 + tmp_33*tmp_526 + tmp_522*tmp_9 - 1.0/4.0;
      real_t tmp_528 = tmp_41*tmp_522 + tmp_42*tmp_524 + tmp_43*tmp_526 - 1.0/4.0;
      real_t tmp_529 = tmp_45*tmp_522 + tmp_46*tmp_524 + tmp_47*tmp_526 - 1.0/4.0;
      real_t tmp_530 = tmp_1*tmp_527 + tmp_2*tmp_528 + tmp_529*tmp_6;
      real_t tmp_531 = tmp_14*tmp_527 + tmp_4*tmp_529 + tmp_528*tmp_7;
      real_t tmp_532 = tmp_11*tmp_529 + tmp_13*tmp_527 + tmp_15*tmp_528;
      real_t tmp_533 = tmp_69*(tmp_521 + tmp_92);
      real_t tmp_534 = tmp_69*(tmp_523 + tmp_94);
      real_t tmp_535 = tmp_69*(tmp_525 + tmp_96);
      real_t tmp_536 = tmp_533*tmp_78 + tmp_534*tmp_75 + tmp_535*tmp_58 - 1.0/4.0;
      real_t tmp_537 = tmp_533*tmp_79 + tmp_534*tmp_76 + tmp_535*tmp_71 - 1.0/4.0;
      real_t tmp_538 = tmp_533*tmp_80 + tmp_534*tmp_77 + tmp_535*tmp_73 - 1.0/4.0;
      real_t tmp_539 = tmp_536*tmp_60 + tmp_537*tmp_61 + tmp_538*tmp_64;
      real_t tmp_540 = tmp_51*tmp_537 + tmp_536*tmp_65 + tmp_538*tmp_55;
      real_t tmp_541 = tmp_53*tmp_538 + tmp_536*tmp_62 + tmp_537*tmp_56;
      real_t a_0_0 = 0.0068572537431980923*tmp_121*(-tmp_101*tmp_105 - tmp_106*tmp_110 - tmp_111*tmp_115 - tmp_120*(tmp_101*tmp_49 + tmp_106*tmp_82 + tmp_111*tmp_87) + tmp_49*tmp_81 + tmp_82*tmp_86 + tmp_87*tmp_91) + 0.037198804536718075*tmp_121*(-tmp_105*tmp_140 - tmp_110*tmp_141 - tmp_115*tmp_142 - tmp_120*(tmp_131*tmp_140 + tmp_132*tmp_141 + tmp_133*tmp_142) + tmp_131*tmp_81 + tmp_132*tmp_86 + tmp_133*tmp_91) + 0.020848748529055869*tmp_121*(-tmp_105*tmp_161 - tmp_110*tmp_162 - tmp_115*tmp_163 - tmp_120*(tmp_152*tmp_161 + tmp_153*tmp_162 + tmp_154*tmp_163) + tmp_152*tmp_81 + tmp_153*tmp_86 + tmp_154*tmp_91) + 0.019202922745021479*tmp_121*(-tmp_105*tmp_182 - tmp_110*tmp_183 - tmp_115*tmp_184 - tmp_120*(tmp_173*tmp_182 + tmp_174*tmp_183 + tmp_175*tmp_184) + tmp_173*tmp_81 + tmp_174*tmp_86 + tmp_175*tmp_91) + 0.020848748529055869*tmp_121*(-tmp_105*tmp_203 - tmp_110*tmp_204 - tmp_115*tmp_205 - tmp_120*(tmp_194*tmp_203 + tmp_195*tmp_204 + tmp_196*tmp_205) + tmp_194*tmp_81 + tmp_195*tmp_86 + tmp_196*tmp_91) + 0.019202922745021479*tmp_121*(-tmp_105*tmp_224 - tmp_110*tmp_225 - tmp_115*tmp_226 - tmp_120*(tmp_215*tmp_224 + tmp_216*tmp_225 + tmp_217*tmp_226) + tmp_215*tmp_81 + tmp_216*tmp_86 + tmp_217*tmp_91) + 0.020848748529055869*tmp_121*(-tmp_105*tmp_245 - tmp_110*tmp_246 - tmp_115*tmp_247 - tmp_120*(tmp_236*tmp_245 + tmp_237*tmp_246 + tmp_238*tmp_247) + tmp_236*tmp_81 + tmp_237*tmp_86 + tmp_238*tmp_91) + 0.019202922745021479*tmp_121*(-tmp_105*tmp_266 - tmp_110*tmp_267 - tmp_115*tmp_268 - tmp_120*(tmp_257*tmp_266 + tmp_258*tmp_267 + tmp_259*tmp_268) + tmp_257*tmp_81 + tmp_258*tmp_86 + tmp_259*tmp_91) + 0.020848748529055869*tmp_121*(-tmp_105*tmp_287 - tmp_110*tmp_288 - tmp_115*tmp_289 - tmp_120*(tmp_278*tmp_287 + tmp_279*tmp_288 + tmp_280*tmp_289) + tmp_278*tmp_81 + tmp_279*tmp_86 + tmp_280*tmp_91) + 0.019202922745021479*tmp_121*(-tmp_105*tmp_308 - tmp_110*tmp_309 - tmp_115*tmp_310 - tmp_120*(tmp_299*tmp_308 + tmp_300*tmp_309 + tmp_301*tmp_310) + tmp_299*tmp_81 + tmp_300*tmp_86 + tmp_301*tmp_91) + 0.020848748529055869*tmp_121*(-tmp_105*tmp_329 - tmp_110*tmp_330 - tmp_115*tmp_331 - tmp_120*(tmp_320*tmp_329 + tmp_321*tmp_330 + tmp_322*tmp_331) + tmp_320*tmp_81 + tmp_321*tmp_86 + tmp_322*tmp_91) + 0.019202922745021479*tmp_121*(-tmp_105*tmp_350 - tmp_110*tmp_351 - tmp_115*tmp_352 - tmp_120*(tmp_341*tmp_350 + tmp_342*tmp_351 + tmp_343*tmp_352) + tmp_341*tmp_81 + tmp_342*tmp_86 + tmp_343*tmp_91) + 0.042507265838595799*tmp_121*(-tmp_105*tmp_371 - tmp_110*tmp_372 - tmp_115*tmp_373 - tmp_120*(tmp_362*tmp_371 + tmp_363*tmp_372 + tmp_364*tmp_373) + tmp_362*tmp_81 + tmp_363*tmp_86 + tmp_364*tmp_91) + 0.020848748529055869*tmp_121*(-tmp_105*tmp_392 - tmp_110*tmp_393 - tmp_115*tmp_394 - tmp_120*(tmp_383*tmp_392 + tmp_384*tmp_393 + tmp_385*tmp_394) + tmp_383*tmp_81 + tmp_384*tmp_86 + tmp_385*tmp_91) + 0.0068572537431980923*tmp_121*(-tmp_105*tmp_413 - tmp_110*tmp_414 - tmp_115*tmp_415 - tmp_120*(tmp_404*tmp_413 + tmp_405*tmp_414 + tmp_406*tmp_415) + tmp_404*tmp_81 + tmp_405*tmp_86 + tmp_406*tmp_91) + 0.037198804536718075*tmp_121*(-tmp_105*tmp_434 - tmp_110*tmp_435 - tmp_115*tmp_436 - tmp_120*(tmp_425*tmp_434 + tmp_426*tmp_435 + tmp_427*tmp_436) + tmp_425*tmp_81 + tmp_426*tmp_86 + tmp_427*tmp_91) + 0.042507265838595799*tmp_121*(-tmp_105*tmp_455 - tmp_110*tmp_456 - tmp_115*tmp_457 - tmp_120*(tmp_446*tmp_455 + tmp_447*tmp_456 + tmp_448*tmp_457) + tmp_446*tmp_81 + tmp_447*tmp_86 + tmp_448*tmp_91) + 0.0068572537431980923*tmp_121*(-tmp_105*tmp_476 - tmp_110*tmp_477 - tmp_115*tmp_478 - tmp_120*(tmp_467*tmp_476 + tmp_468*tmp_477 + tmp_469*tmp_478) + tmp_467*tmp_81 + tmp_468*tmp_86 + tmp_469*tmp_91) + 0.037198804536718075*tmp_121*(-tmp_105*tmp_497 - tmp_110*tmp_498 - tmp_115*tmp_499 - tmp_120*(tmp_488*tmp_497 + tmp_489*tmp_498 + tmp_490*tmp_499) + tmp_488*tmp_81 + tmp_489*tmp_86 + tmp_490*tmp_91) + 0.042507265838595799*tmp_121*(-tmp_105*tmp_518 - tmp_110*tmp_519 - tmp_115*tmp_520 - tmp_120*(tmp_509*tmp_518 + tmp_510*tmp_519 + tmp_511*tmp_520) + tmp_509*tmp_81 + tmp_510*tmp_86 + tmp_511*tmp_91) + 0.019202922745021479*tmp_121*(-tmp_105*tmp_539 - tmp_110*tmp_540 - tmp_115*tmp_541 - tmp_120*(tmp_530*tmp_539 + tmp_531*tmp_540 + tmp_532*tmp_541) + tmp_530*tmp_81 + tmp_531*tmp_86 + tmp_532*tmp_91);
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
      real_t tmp_49 = 5.0*std::pow((std::abs(tmp_22*tmp_46 - tmp_28*tmp_48)*std::abs(tmp_22*tmp_46 - tmp_28*tmp_48)) + (std::abs(tmp_22*tmp_47 - tmp_34*tmp_48)*std::abs(tmp_22*tmp_47 - tmp_34*tmp_48)) + (std::abs(tmp_28*tmp_47 - tmp_34*tmp_46)*std::abs(tmp_28*tmp_47 - tmp_34*tmp_46)), 0.25);
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
      real_t a_0_0 = 0.019202922745021479*tmp_49*(((tmp_1*tmp_101 + tmp_102*tmp_2 + tmp_103*tmp_6)*(tmp_1*tmp_101 + tmp_102*tmp_2 + tmp_103*tmp_6)) + ((tmp_101*tmp_13 + tmp_102*tmp_15 + tmp_103*tmp_11)*(tmp_101*tmp_13 + tmp_102*tmp_15 + tmp_103*tmp_11)) + ((tmp_101*tmp_14 + tmp_102*tmp_7 + tmp_103*tmp_4)*(tmp_101*tmp_14 + tmp_102*tmp_7 + tmp_103*tmp_4))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_107 + tmp_108*tmp_2 + tmp_109*tmp_6)*(tmp_1*tmp_107 + tmp_108*tmp_2 + tmp_109*tmp_6)) + ((tmp_107*tmp_13 + tmp_108*tmp_15 + tmp_109*tmp_11)*(tmp_107*tmp_13 + tmp_108*tmp_15 + tmp_109*tmp_11)) + ((tmp_107*tmp_14 + tmp_108*tmp_7 + tmp_109*tmp_4)*(tmp_107*tmp_14 + tmp_108*tmp_7 + tmp_109*tmp_4))) + 0.019202922745021479*tmp_49*(((tmp_1*tmp_113 + tmp_114*tmp_2 + tmp_115*tmp_6)*(tmp_1*tmp_113 + tmp_114*tmp_2 + tmp_115*tmp_6)) + ((tmp_11*tmp_115 + tmp_113*tmp_13 + tmp_114*tmp_15)*(tmp_11*tmp_115 + tmp_113*tmp_13 + tmp_114*tmp_15)) + ((tmp_113*tmp_14 + tmp_114*tmp_7 + tmp_115*tmp_4)*(tmp_113*tmp_14 + tmp_114*tmp_7 + tmp_115*tmp_4))) + 0.042507265838595799*tmp_49*(((tmp_1*tmp_119 + tmp_120*tmp_2 + tmp_121*tmp_6)*(tmp_1*tmp_119 + tmp_120*tmp_2 + tmp_121*tmp_6)) + ((tmp_11*tmp_121 + tmp_119*tmp_13 + tmp_120*tmp_15)*(tmp_11*tmp_121 + tmp_119*tmp_13 + tmp_120*tmp_15)) + ((tmp_119*tmp_14 + tmp_120*tmp_7 + tmp_121*tmp_4)*(tmp_119*tmp_14 + tmp_120*tmp_7 + tmp_121*tmp_4))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_125 + tmp_126*tmp_2 + tmp_127*tmp_6)*(tmp_1*tmp_125 + tmp_126*tmp_2 + tmp_127*tmp_6)) + ((tmp_11*tmp_127 + tmp_125*tmp_13 + tmp_126*tmp_15)*(tmp_11*tmp_127 + tmp_125*tmp_13 + tmp_126*tmp_15)) + ((tmp_125*tmp_14 + tmp_126*tmp_7 + tmp_127*tmp_4)*(tmp_125*tmp_14 + tmp_126*tmp_7 + tmp_127*tmp_4))) + 0.0068572537431980923*tmp_49*(((tmp_1*tmp_131 + tmp_132*tmp_2 + tmp_133*tmp_6)*(tmp_1*tmp_131 + tmp_132*tmp_2 + tmp_133*tmp_6)) + ((tmp_11*tmp_133 + tmp_13*tmp_131 + tmp_132*tmp_15)*(tmp_11*tmp_133 + tmp_13*tmp_131 + tmp_132*tmp_15)) + ((tmp_131*tmp_14 + tmp_132*tmp_7 + tmp_133*tmp_4)*(tmp_131*tmp_14 + tmp_132*tmp_7 + tmp_133*tmp_4))) + 0.037198804536718075*tmp_49*(((tmp_1*tmp_137 + tmp_138*tmp_2 + tmp_139*tmp_6)*(tmp_1*tmp_137 + tmp_138*tmp_2 + tmp_139*tmp_6)) + ((tmp_11*tmp_139 + tmp_13*tmp_137 + tmp_138*tmp_15)*(tmp_11*tmp_139 + tmp_13*tmp_137 + tmp_138*tmp_15)) + ((tmp_137*tmp_14 + tmp_138*tmp_7 + tmp_139*tmp_4)*(tmp_137*tmp_14 + tmp_138*tmp_7 + tmp_139*tmp_4))) + 0.042507265838595799*tmp_49*(((tmp_1*tmp_143 + tmp_144*tmp_2 + tmp_145*tmp_6)*(tmp_1*tmp_143 + tmp_144*tmp_2 + tmp_145*tmp_6)) + ((tmp_11*tmp_145 + tmp_13*tmp_143 + tmp_144*tmp_15)*(tmp_11*tmp_145 + tmp_13*tmp_143 + tmp_144*tmp_15)) + ((tmp_14*tmp_143 + tmp_144*tmp_7 + tmp_145*tmp_4)*(tmp_14*tmp_143 + tmp_144*tmp_7 + tmp_145*tmp_4))) + 0.0068572537431980923*tmp_49*(((tmp_1*tmp_149 + tmp_150*tmp_2 + tmp_151*tmp_6)*(tmp_1*tmp_149 + tmp_150*tmp_2 + tmp_151*tmp_6)) + ((tmp_11*tmp_151 + tmp_13*tmp_149 + tmp_15*tmp_150)*(tmp_11*tmp_151 + tmp_13*tmp_149 + tmp_15*tmp_150)) + ((tmp_14*tmp_149 + tmp_150*tmp_7 + tmp_151*tmp_4)*(tmp_14*tmp_149 + tmp_150*tmp_7 + tmp_151*tmp_4))) + 0.037198804536718075*tmp_49*(((tmp_1*tmp_155 + tmp_156*tmp_2 + tmp_157*tmp_6)*(tmp_1*tmp_155 + tmp_156*tmp_2 + tmp_157*tmp_6)) + ((tmp_11*tmp_157 + tmp_13*tmp_155 + tmp_15*tmp_156)*(tmp_11*tmp_157 + tmp_13*tmp_155 + tmp_15*tmp_156)) + ((tmp_14*tmp_155 + tmp_156*tmp_7 + tmp_157*tmp_4)*(tmp_14*tmp_155 + tmp_156*tmp_7 + tmp_157*tmp_4))) + 0.042507265838595799*tmp_49*(((tmp_1*tmp_161 + tmp_162*tmp_2 + tmp_163*tmp_6)*(tmp_1*tmp_161 + tmp_162*tmp_2 + tmp_163*tmp_6)) + ((tmp_11*tmp_163 + tmp_13*tmp_161 + tmp_15*tmp_162)*(tmp_11*tmp_163 + tmp_13*tmp_161 + tmp_15*tmp_162)) + ((tmp_14*tmp_161 + tmp_162*tmp_7 + tmp_163*tmp_4)*(tmp_14*tmp_161 + tmp_162*tmp_7 + tmp_163*tmp_4))) + 0.019202922745021479*tmp_49*(((tmp_1*tmp_167 + tmp_168*tmp_2 + tmp_169*tmp_6)*(tmp_1*tmp_167 + tmp_168*tmp_2 + tmp_169*tmp_6)) + ((tmp_11*tmp_169 + tmp_13*tmp_167 + tmp_15*tmp_168)*(tmp_11*tmp_169 + tmp_13*tmp_167 + tmp_15*tmp_168)) + ((tmp_14*tmp_167 + tmp_168*tmp_7 + tmp_169*tmp_4)*(tmp_14*tmp_167 + tmp_168*tmp_7 + tmp_169*tmp_4))) + 0.0068572537431980923*tmp_49*(((tmp_1*tmp_37 + tmp_2*tmp_41 + tmp_45*tmp_6)*(tmp_1*tmp_37 + tmp_2*tmp_41 + tmp_45*tmp_6)) + ((tmp_11*tmp_45 + tmp_13*tmp_37 + tmp_15*tmp_41)*(tmp_11*tmp_45 + tmp_13*tmp_37 + tmp_15*tmp_41)) + ((tmp_14*tmp_37 + tmp_4*tmp_45 + tmp_41*tmp_7)*(tmp_14*tmp_37 + tmp_4*tmp_45 + tmp_41*tmp_7))) + 0.037198804536718075*tmp_49*(((tmp_1*tmp_53 + tmp_2*tmp_54 + tmp_55*tmp_6)*(tmp_1*tmp_53 + tmp_2*tmp_54 + tmp_55*tmp_6)) + ((tmp_11*tmp_55 + tmp_13*tmp_53 + tmp_15*tmp_54)*(tmp_11*tmp_55 + tmp_13*tmp_53 + tmp_15*tmp_54)) + ((tmp_14*tmp_53 + tmp_4*tmp_55 + tmp_54*tmp_7)*(tmp_14*tmp_53 + tmp_4*tmp_55 + tmp_54*tmp_7))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_59 + tmp_2*tmp_60 + tmp_6*tmp_61)*(tmp_1*tmp_59 + tmp_2*tmp_60 + tmp_6*tmp_61)) + ((tmp_11*tmp_61 + tmp_13*tmp_59 + tmp_15*tmp_60)*(tmp_11*tmp_61 + tmp_13*tmp_59 + tmp_15*tmp_60)) + ((tmp_14*tmp_59 + tmp_4*tmp_61 + tmp_60*tmp_7)*(tmp_14*tmp_59 + tmp_4*tmp_61 + tmp_60*tmp_7))) + 0.019202922745021479*tmp_49*(((tmp_1*tmp_65 + tmp_2*tmp_66 + tmp_6*tmp_67)*(tmp_1*tmp_65 + tmp_2*tmp_66 + tmp_6*tmp_67)) + ((tmp_11*tmp_67 + tmp_13*tmp_65 + tmp_15*tmp_66)*(tmp_11*tmp_67 + tmp_13*tmp_65 + tmp_15*tmp_66)) + ((tmp_14*tmp_65 + tmp_4*tmp_67 + tmp_66*tmp_7)*(tmp_14*tmp_65 + tmp_4*tmp_67 + tmp_66*tmp_7))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_71 + tmp_2*tmp_72 + tmp_6*tmp_73)*(tmp_1*tmp_71 + tmp_2*tmp_72 + tmp_6*tmp_73)) + ((tmp_11*tmp_73 + tmp_13*tmp_71 + tmp_15*tmp_72)*(tmp_11*tmp_73 + tmp_13*tmp_71 + tmp_15*tmp_72)) + ((tmp_14*tmp_71 + tmp_4*tmp_73 + tmp_7*tmp_72)*(tmp_14*tmp_71 + tmp_4*tmp_73 + tmp_7*tmp_72))) + 0.019202922745021479*tmp_49*(((tmp_1*tmp_77 + tmp_2*tmp_78 + tmp_6*tmp_79)*(tmp_1*tmp_77 + tmp_2*tmp_78 + tmp_6*tmp_79)) + ((tmp_11*tmp_79 + tmp_13*tmp_77 + tmp_15*tmp_78)*(tmp_11*tmp_79 + tmp_13*tmp_77 + tmp_15*tmp_78)) + ((tmp_14*tmp_77 + tmp_4*tmp_79 + tmp_7*tmp_78)*(tmp_14*tmp_77 + tmp_4*tmp_79 + tmp_7*tmp_78))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_83 + tmp_2*tmp_84 + tmp_6*tmp_85)*(tmp_1*tmp_83 + tmp_2*tmp_84 + tmp_6*tmp_85)) + ((tmp_11*tmp_85 + tmp_13*tmp_83 + tmp_15*tmp_84)*(tmp_11*tmp_85 + tmp_13*tmp_83 + tmp_15*tmp_84)) + ((tmp_14*tmp_83 + tmp_4*tmp_85 + tmp_7*tmp_84)*(tmp_14*tmp_83 + tmp_4*tmp_85 + tmp_7*tmp_84))) + 0.019202922745021479*tmp_49*(((tmp_1*tmp_89 + tmp_2*tmp_90 + tmp_6*tmp_91)*(tmp_1*tmp_89 + tmp_2*tmp_90 + tmp_6*tmp_91)) + ((tmp_11*tmp_91 + tmp_13*tmp_89 + tmp_15*tmp_90)*(tmp_11*tmp_91 + tmp_13*tmp_89 + tmp_15*tmp_90)) + ((tmp_14*tmp_89 + tmp_4*tmp_91 + tmp_7*tmp_90)*(tmp_14*tmp_89 + tmp_4*tmp_91 + tmp_7*tmp_90))) + 0.020848748529055869*tmp_49*(((tmp_1*tmp_95 + tmp_2*tmp_96 + tmp_6*tmp_97)*(tmp_1*tmp_95 + tmp_2*tmp_96 + tmp_6*tmp_97)) + ((tmp_11*tmp_97 + tmp_13*tmp_95 + tmp_15*tmp_96)*(tmp_11*tmp_97 + tmp_13*tmp_95 + tmp_15*tmp_96)) + ((tmp_14*tmp_95 + tmp_4*tmp_97 + tmp_7*tmp_96)*(tmp_14*tmp_95 + tmp_4*tmp_97 + tmp_7*tmp_96)));
      elMat( 0, 0) = a_0_0;
   }


};


} //eg
} // dg
} // hyteg

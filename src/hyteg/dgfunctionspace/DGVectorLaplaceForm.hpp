
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
namespace dg {

class DGVectorLaplaceFormP1P1_00 : public hyteg::dg::DGForm2D
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
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 );
      real_t tmp_5  = 1.0 / ( tmp_4 );
      real_t tmp_6  = tmp_1 * tmp_5;
      real_t tmp_7  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_8  = tmp_5 * tmp_7;
      real_t tmp_9  = -tmp_6 - tmp_8;
      real_t tmp_10 = tmp_3 * tmp_5;
      real_t tmp_11 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_12 = tmp_11 * tmp_5;
      real_t tmp_13 = -tmp_10 - tmp_12;
      real_t tmp_14 = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                                p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_15 = tmp_14 * ( ( tmp_13 * tmp_13 ) + ( tmp_9 * tmp_9 ) );
      real_t tmp_16 = tmp_14 * ( tmp_10 * tmp_13 + tmp_8 * tmp_9 );
      real_t tmp_17 = 0.5 * tmp_16;
      real_t tmp_18 = tmp_14 * ( tmp_12 * tmp_13 + tmp_6 * tmp_9 );
      real_t tmp_19 = 0.5 * tmp_18;
      real_t tmp_20 = 1.0 / ( tmp_4 * tmp_4 );
      real_t tmp_21 = tmp_14 * ( tmp_20 * ( tmp_3 * tmp_3 ) + tmp_20 * ( tmp_7 * tmp_7 ) );
      real_t tmp_22 = tmp_14 * ( tmp_1 * tmp_20 * tmp_7 + tmp_11 * tmp_20 * tmp_3 );
      real_t tmp_23 = 0.5 * tmp_22;
      real_t tmp_24 = tmp_14 * ( ( tmp_1 * tmp_1 ) * tmp_20 + ( tmp_11 * tmp_11 ) * tmp_20 );
      real_t a_0_0  = 0.5 * tmp_15;
      real_t a_0_1  = tmp_17;
      real_t a_0_2  = tmp_19;
      real_t a_1_0  = tmp_17;
      real_t a_1_1  = 0.5 * tmp_21;
      real_t a_1_2  = tmp_23;
      real_t a_2_0  = tmp_19;
      real_t a_2_1  = tmp_23;
      real_t a_2_2  = 0.5 * tmp_24;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;
      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
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

      real_t tmp_0  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_1  = -p_affine_0_1;
      real_t tmp_2  = p_affine_6_1 + tmp_1;
      real_t tmp_3  = 0.046910077030668018 * tmp_0 + tmp_2;
      real_t tmp_4  = -p_affine_0_0;
      real_t tmp_5  = p_affine_1_0 + tmp_4;
      real_t tmp_6  = p_affine_2_1 + tmp_1;
      real_t tmp_7  = 1.0 / ( tmp_5 * tmp_6 - ( p_affine_1_1 + tmp_1 ) * ( p_affine_2_0 + tmp_4 ) );
      real_t tmp_8  = tmp_5 * tmp_7;
      real_t tmp_9  = tmp_3 * tmp_8;
      real_t tmp_10 = tmp_7 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_11 = tmp_10 * tmp_3;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_4;
      real_t tmp_14 = 0.046910077030668018 * tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6 * tmp_7;
      real_t tmp_16 = tmp_14 * tmp_15;
      real_t tmp_17 = tmp_7 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_18 = tmp_14 * tmp_17;
      real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20 = std::abs( std::pow( ( tmp_0 * tmp_0 ) + ( tmp_12 * tmp_12 ), 1.0 / 2.0 ) );
      real_t tmp_21 = 6 / tmp_20;
      real_t tmp_22 = p_affine_10_0 * ( -tmp_15 - tmp_17 ) + p_affine_10_1 * ( -tmp_10 - tmp_8 );
      real_t tmp_23 = 1.0 * tmp_22;
      real_t tmp_24 = 0.11846344252809471 * tmp_20;
      real_t tmp_25 = 0.23076534494715845 * tmp_0 + tmp_2;
      real_t tmp_26 = tmp_25 * tmp_8;
      real_t tmp_27 = tmp_10 * tmp_25;
      real_t tmp_28 = 0.23076534494715845 * tmp_12 + tmp_13;
      real_t tmp_29 = tmp_15 * tmp_28;
      real_t tmp_30 = tmp_17 * tmp_28;
      real_t tmp_31 = -tmp_26 - tmp_27 - tmp_29 - tmp_30 + 1;
      real_t tmp_32 = 0.2393143352496831 * tmp_20;
      real_t tmp_33 = 0.5 * tmp_0 + tmp_2;
      real_t tmp_34 = tmp_33 * tmp_8;
      real_t tmp_35 = tmp_10 * tmp_33;
      real_t tmp_36 = 0.5 * tmp_12 + tmp_13;
      real_t tmp_37 = tmp_15 * tmp_36;
      real_t tmp_38 = tmp_17 * tmp_36;
      real_t tmp_39 = -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1;
      real_t tmp_40 = 0.2844444444444445 * tmp_20;
      real_t tmp_41 = 0.7692346550528415 * tmp_0 + tmp_2;
      real_t tmp_42 = tmp_41 * tmp_8;
      real_t tmp_43 = tmp_10 * tmp_41;
      real_t tmp_44 = 0.7692346550528415 * tmp_12 + tmp_13;
      real_t tmp_45 = tmp_15 * tmp_44;
      real_t tmp_46 = tmp_17 * tmp_44;
      real_t tmp_47 = -tmp_42 - tmp_43 - tmp_45 - tmp_46 + 1;
      real_t tmp_48 = 0.2393143352496831 * tmp_20;
      real_t tmp_49 = 0.95308992296933193 * tmp_0 + tmp_2;
      real_t tmp_50 = tmp_49 * tmp_8;
      real_t tmp_51 = tmp_10 * tmp_49;
      real_t tmp_52 = 0.95308992296933193 * tmp_12 + tmp_13;
      real_t tmp_53 = tmp_15 * tmp_52;
      real_t tmp_54 = tmp_17 * tmp_52;
      real_t tmp_55 = -tmp_50 - tmp_51 - tmp_53 - tmp_54 + 1;
      real_t tmp_56 = 0.11846344252809471 * tmp_20;
      real_t tmp_57 = tmp_11 + tmp_16;
      real_t tmp_58 = 0.5 * tmp_22;
      real_t tmp_59 = p_affine_10_0 * tmp_15 + p_affine_10_1 * tmp_10;
      real_t tmp_60 = 0.5 * tmp_59;
      real_t tmp_61 = tmp_19 * tmp_21;
      real_t tmp_62 = tmp_27 + tmp_29;
      real_t tmp_63 = tmp_21 * tmp_31;
      real_t tmp_64 = tmp_35 + tmp_37;
      real_t tmp_65 = tmp_21 * tmp_39;
      real_t tmp_66 = tmp_43 + tmp_45;
      real_t tmp_67 = tmp_21 * tmp_47;
      real_t tmp_68 = tmp_51 + tmp_53;
      real_t tmp_69 = tmp_21 * tmp_55;
      real_t tmp_70 = tmp_24 * ( -tmp_19 * tmp_60 - tmp_57 * tmp_58 + tmp_57 * tmp_61 ) +
                      tmp_32 * ( -tmp_31 * tmp_60 - tmp_58 * tmp_62 + tmp_62 * tmp_63 ) +
                      tmp_40 * ( -tmp_39 * tmp_60 - tmp_58 * tmp_64 + tmp_64 * tmp_65 ) +
                      tmp_48 * ( -tmp_47 * tmp_60 - tmp_58 * tmp_66 + tmp_66 * tmp_67 ) +
                      tmp_56 * ( -tmp_55 * tmp_60 - tmp_58 * tmp_68 + tmp_68 * tmp_69 );
      real_t tmp_71 = tmp_18 + tmp_9;
      real_t tmp_72 = p_affine_10_0 * tmp_17 + p_affine_10_1 * tmp_8;
      real_t tmp_73 = 0.5 * tmp_72;
      real_t tmp_74 = tmp_26 + tmp_30;
      real_t tmp_75 = tmp_34 + tmp_38;
      real_t tmp_76 = tmp_42 + tmp_46;
      real_t tmp_77 = tmp_50 + tmp_54;
      real_t tmp_78 = tmp_24 * ( -tmp_19 * tmp_73 - tmp_58 * tmp_71 + tmp_61 * tmp_71 ) +
                      tmp_32 * ( -tmp_31 * tmp_73 - tmp_58 * tmp_74 + tmp_63 * tmp_74 ) +
                      tmp_40 * ( -tmp_39 * tmp_73 - tmp_58 * tmp_75 + tmp_65 * tmp_75 ) +
                      tmp_48 * ( -tmp_47 * tmp_73 - tmp_58 * tmp_76 + tmp_67 * tmp_76 ) +
                      tmp_56 * ( -tmp_55 * tmp_73 - tmp_58 * tmp_77 + tmp_69 * tmp_77 );
      real_t tmp_79 = 1.0 * tmp_59;
      real_t tmp_80 = tmp_24 * ( tmp_21 * tmp_57 * tmp_71 - tmp_57 * tmp_73 - tmp_60 * tmp_71 ) +
                      tmp_32 * ( tmp_21 * tmp_62 * tmp_74 - tmp_60 * tmp_74 - tmp_62 * tmp_73 ) +
                      tmp_40 * ( tmp_21 * tmp_64 * tmp_75 - tmp_60 * tmp_75 - tmp_64 * tmp_73 ) +
                      tmp_48 * ( tmp_21 * tmp_66 * tmp_76 - tmp_60 * tmp_76 - tmp_66 * tmp_73 ) +
                      tmp_56 * ( tmp_21 * tmp_68 * tmp_77 - tmp_60 * tmp_77 - tmp_68 * tmp_73 );
      real_t tmp_81 = 1.0 * tmp_72;
      real_t a_0_0  = tmp_24 * ( ( tmp_19 * tmp_19 ) * tmp_21 - tmp_19 * tmp_23 ) +
                     tmp_32 * ( tmp_21 * ( tmp_31 * tmp_31 ) - tmp_23 * tmp_31 ) +
                     tmp_40 * ( tmp_21 * ( tmp_39 * tmp_39 ) - tmp_23 * tmp_39 ) +
                     tmp_48 * ( tmp_21 * ( tmp_47 * tmp_47 ) - tmp_23 * tmp_47 ) +
                     tmp_56 * ( tmp_21 * ( tmp_55 * tmp_55 ) - tmp_23 * tmp_55 );
      real_t a_0_1 = tmp_70;
      real_t a_0_2 = tmp_78;
      real_t a_1_0 = tmp_70;
      real_t a_1_1 = tmp_24 * ( tmp_21 * ( tmp_57 * tmp_57 ) - tmp_57 * tmp_79 ) +
                     tmp_32 * ( tmp_21 * ( tmp_62 * tmp_62 ) - tmp_62 * tmp_79 ) +
                     tmp_40 * ( tmp_21 * ( tmp_64 * tmp_64 ) - tmp_64 * tmp_79 ) +
                     tmp_48 * ( tmp_21 * ( tmp_66 * tmp_66 ) - tmp_66 * tmp_79 ) +
                     tmp_56 * ( tmp_21 * ( tmp_68 * tmp_68 ) - tmp_68 * tmp_79 );
      real_t a_1_2 = tmp_80;
      real_t a_2_0 = tmp_78;
      real_t a_2_1 = tmp_80;
      real_t a_2_2 = tmp_24 * ( tmp_21 * ( tmp_71 * tmp_71 ) - tmp_71 * tmp_81 ) +
                     tmp_32 * ( tmp_21 * ( tmp_74 * tmp_74 ) - tmp_74 * tmp_81 ) +
                     tmp_40 * ( tmp_21 * ( tmp_75 * tmp_75 ) - tmp_75 * tmp_81 ) +
                     tmp_48 * ( tmp_21 * ( tmp_76 * tmp_76 ) - tmp_76 * tmp_81 ) +
                     tmp_56 * ( tmp_21 * ( tmp_77 * tmp_77 ) - tmp_77 * tmp_81 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;
      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
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

      real_t tmp_0   = -p_affine_3_1;
      real_t tmp_1   = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2   = p_affine_6_1 + 0.046910077030668018 * tmp_1;
      real_t tmp_3   = tmp_0 + tmp_2;
      real_t tmp_4   = -p_affine_3_0;
      real_t tmp_5   = p_affine_4_0 + tmp_4;
      real_t tmp_6   = p_affine_5_1 + tmp_0;
      real_t tmp_7   = 1.0 / ( tmp_5 * tmp_6 - ( p_affine_4_1 + tmp_0 ) * ( p_affine_5_0 + tmp_4 ) );
      real_t tmp_8   = tmp_5 * tmp_7;
      real_t tmp_9   = tmp_3 * tmp_8;
      real_t tmp_10  = tmp_7 * ( p_affine_3_0 - p_affine_5_0 );
      real_t tmp_11  = tmp_10 * tmp_3;
      real_t tmp_12  = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13  = p_affine_6_0 + 0.046910077030668018 * tmp_12;
      real_t tmp_14  = tmp_13 + tmp_4;
      real_t tmp_15  = tmp_6 * tmp_7;
      real_t tmp_16  = tmp_14 * tmp_15;
      real_t tmp_17  = tmp_7 * ( p_affine_3_1 - p_affine_4_1 );
      real_t tmp_18  = tmp_14 * tmp_17;
      real_t tmp_19  = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20  = -p_affine_0_1;
      real_t tmp_21  = p_affine_2_1 + tmp_20;
      real_t tmp_22  = -p_affine_0_0;
      real_t tmp_23  = p_affine_1_0 + tmp_22;
      real_t tmp_24  = 1.0 / ( tmp_21 * tmp_23 - ( p_affine_1_1 + tmp_20 ) * ( p_affine_2_0 + tmp_22 ) );
      real_t tmp_25  = tmp_21 * tmp_24;
      real_t tmp_26  = tmp_24 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_27  = tmp_23 * tmp_24;
      real_t tmp_28  = tmp_24 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_29  = 0.5 * p_affine_10_0 * ( -tmp_25 - tmp_26 ) + 0.5 * p_affine_10_1 * ( -tmp_27 - tmp_28 );
      real_t tmp_30  = tmp_2 + tmp_20;
      real_t tmp_31  = tmp_27 * tmp_30;
      real_t tmp_32  = tmp_28 * tmp_30;
      real_t tmp_33  = tmp_13 + tmp_22;
      real_t tmp_34  = tmp_25 * tmp_33;
      real_t tmp_35  = tmp_26 * tmp_33;
      real_t tmp_36  = -tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1;
      real_t tmp_37  = 0.5 * p_affine_10_0 * ( -tmp_15 - tmp_17 ) + 0.5 * p_affine_10_1 * ( -tmp_10 - tmp_8 );
      real_t tmp_38  = std::abs( std::pow( ( tmp_1 * tmp_1 ) + ( tmp_12 * tmp_12 ), 1.0 / 2.0 ) );
      real_t tmp_39  = 6 / tmp_38;
      real_t tmp_40  = tmp_36 * tmp_39;
      real_t tmp_41  = 0.11846344252809471 * tmp_38;
      real_t tmp_42  = p_affine_6_1 + 0.23076534494715845 * tmp_1;
      real_t tmp_43  = tmp_0 + tmp_42;
      real_t tmp_44  = tmp_43 * tmp_8;
      real_t tmp_45  = tmp_10 * tmp_43;
      real_t tmp_46  = p_affine_6_0 + 0.23076534494715845 * tmp_12;
      real_t tmp_47  = tmp_4 + tmp_46;
      real_t tmp_48  = tmp_15 * tmp_47;
      real_t tmp_49  = tmp_17 * tmp_47;
      real_t tmp_50  = -tmp_44 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_51  = tmp_20 + tmp_42;
      real_t tmp_52  = tmp_27 * tmp_51;
      real_t tmp_53  = tmp_28 * tmp_51;
      real_t tmp_54  = tmp_22 + tmp_46;
      real_t tmp_55  = tmp_25 * tmp_54;
      real_t tmp_56  = tmp_26 * tmp_54;
      real_t tmp_57  = -tmp_52 - tmp_53 - tmp_55 - tmp_56 + 1;
      real_t tmp_58  = tmp_39 * tmp_57;
      real_t tmp_59  = 0.2393143352496831 * tmp_38;
      real_t tmp_60  = p_affine_6_1 + 0.5 * tmp_1;
      real_t tmp_61  = tmp_0 + tmp_60;
      real_t tmp_62  = tmp_61 * tmp_8;
      real_t tmp_63  = tmp_10 * tmp_61;
      real_t tmp_64  = p_affine_6_0 + 0.5 * tmp_12;
      real_t tmp_65  = tmp_4 + tmp_64;
      real_t tmp_66  = tmp_15 * tmp_65;
      real_t tmp_67  = tmp_17 * tmp_65;
      real_t tmp_68  = -tmp_62 - tmp_63 - tmp_66 - tmp_67 + 1;
      real_t tmp_69  = tmp_20 + tmp_60;
      real_t tmp_70  = tmp_27 * tmp_69;
      real_t tmp_71  = tmp_28 * tmp_69;
      real_t tmp_72  = tmp_22 + tmp_64;
      real_t tmp_73  = tmp_25 * tmp_72;
      real_t tmp_74  = tmp_26 * tmp_72;
      real_t tmp_75  = -tmp_70 - tmp_71 - tmp_73 - tmp_74 + 1;
      real_t tmp_76  = tmp_39 * tmp_75;
      real_t tmp_77  = 0.2844444444444445 * tmp_38;
      real_t tmp_78  = p_affine_6_1 + 0.7692346550528415 * tmp_1;
      real_t tmp_79  = tmp_0 + tmp_78;
      real_t tmp_80  = tmp_79 * tmp_8;
      real_t tmp_81  = tmp_10 * tmp_79;
      real_t tmp_82  = p_affine_6_0 + 0.7692346550528415 * tmp_12;
      real_t tmp_83  = tmp_4 + tmp_82;
      real_t tmp_84  = tmp_15 * tmp_83;
      real_t tmp_85  = tmp_17 * tmp_83;
      real_t tmp_86  = -tmp_80 - tmp_81 - tmp_84 - tmp_85 + 1;
      real_t tmp_87  = tmp_20 + tmp_78;
      real_t tmp_88  = tmp_27 * tmp_87;
      real_t tmp_89  = tmp_28 * tmp_87;
      real_t tmp_90  = tmp_22 + tmp_82;
      real_t tmp_91  = tmp_25 * tmp_90;
      real_t tmp_92  = tmp_26 * tmp_90;
      real_t tmp_93  = -tmp_88 - tmp_89 - tmp_91 - tmp_92 + 1;
      real_t tmp_94  = tmp_39 * tmp_93;
      real_t tmp_95  = 0.2393143352496831 * tmp_38;
      real_t tmp_96  = p_affine_6_1 + 0.95308992296933193 * tmp_1;
      real_t tmp_97  = tmp_0 + tmp_96;
      real_t tmp_98  = tmp_8 * tmp_97;
      real_t tmp_99  = tmp_10 * tmp_97;
      real_t tmp_100 = p_affine_6_0 + 0.95308992296933193 * tmp_12;
      real_t tmp_101 = tmp_100 + tmp_4;
      real_t tmp_102 = tmp_101 * tmp_15;
      real_t tmp_103 = tmp_101 * tmp_17;
      real_t tmp_104 = -tmp_102 - tmp_103 - tmp_98 - tmp_99 + 1;
      real_t tmp_105 = tmp_20 + tmp_96;
      real_t tmp_106 = tmp_105 * tmp_27;
      real_t tmp_107 = tmp_105 * tmp_28;
      real_t tmp_108 = tmp_100 + tmp_22;
      real_t tmp_109 = tmp_108 * tmp_25;
      real_t tmp_110 = tmp_108 * tmp_26;
      real_t tmp_111 = -tmp_106 - tmp_107 - tmp_109 - tmp_110 + 1;
      real_t tmp_112 = tmp_111 * tmp_39;
      real_t tmp_113 = 0.11846344252809471 * tmp_38;
      real_t tmp_114 = tmp_11 + tmp_16;
      real_t tmp_115 = 0.5 * p_affine_10_0 * tmp_15 + 0.5 * p_affine_10_1 * tmp_10;
      real_t tmp_116 = tmp_45 + tmp_48;
      real_t tmp_117 = tmp_63 + tmp_66;
      real_t tmp_118 = tmp_81 + tmp_84;
      real_t tmp_119 = tmp_102 + tmp_99;
      real_t tmp_120 = tmp_18 + tmp_9;
      real_t tmp_121 = 0.5 * p_affine_10_0 * tmp_17 + 0.5 * p_affine_10_1 * tmp_8;
      real_t tmp_122 = tmp_44 + tmp_49;
      real_t tmp_123 = tmp_62 + tmp_67;
      real_t tmp_124 = tmp_80 + tmp_85;
      real_t tmp_125 = tmp_103 + tmp_98;
      real_t tmp_126 = tmp_32 + tmp_34;
      real_t tmp_127 = 0.5 * p_affine_10_0 * tmp_25 + 0.5 * p_affine_10_1 * tmp_28;
      real_t tmp_128 = tmp_126 * tmp_39;
      real_t tmp_129 = tmp_53 + tmp_55;
      real_t tmp_130 = tmp_129 * tmp_39;
      real_t tmp_131 = tmp_71 + tmp_73;
      real_t tmp_132 = tmp_131 * tmp_39;
      real_t tmp_133 = tmp_89 + tmp_91;
      real_t tmp_134 = tmp_133 * tmp_39;
      real_t tmp_135 = tmp_107 + tmp_109;
      real_t tmp_136 = tmp_135 * tmp_39;
      real_t tmp_137 = tmp_31 + tmp_35;
      real_t tmp_138 = 0.5 * p_affine_10_0 * tmp_26 + 0.5 * p_affine_10_1 * tmp_27;
      real_t tmp_139 = tmp_137 * tmp_39;
      real_t tmp_140 = tmp_52 + tmp_56;
      real_t tmp_141 = tmp_140 * tmp_39;
      real_t tmp_142 = tmp_70 + tmp_74;
      real_t tmp_143 = tmp_142 * tmp_39;
      real_t tmp_144 = tmp_88 + tmp_92;
      real_t tmp_145 = tmp_144 * tmp_39;
      real_t tmp_146 = tmp_106 + tmp_110;
      real_t tmp_147 = tmp_146 * tmp_39;
      real_t a_0_0   = tmp_113 * ( -tmp_104 * tmp_112 + tmp_104 * tmp_29 - tmp_111 * tmp_37 ) +
                     tmp_41 * ( tmp_19 * tmp_29 - tmp_19 * tmp_40 - tmp_36 * tmp_37 ) +
                     tmp_59 * ( tmp_29 * tmp_50 - tmp_37 * tmp_57 - tmp_50 * tmp_58 ) +
                     tmp_77 * ( tmp_29 * tmp_68 - tmp_37 * tmp_75 - tmp_68 * tmp_76 ) +
                     tmp_95 * ( tmp_29 * tmp_86 - tmp_37 * tmp_93 - tmp_86 * tmp_94 );
      real_t a_0_1 = tmp_113 * ( -tmp_111 * tmp_115 - tmp_112 * tmp_119 + tmp_119 * tmp_29 ) +
                     tmp_41 * ( tmp_114 * tmp_29 - tmp_114 * tmp_40 - tmp_115 * tmp_36 ) +
                     tmp_59 * ( -tmp_115 * tmp_57 + tmp_116 * tmp_29 - tmp_116 * tmp_58 ) +
                     tmp_77 * ( -tmp_115 * tmp_75 + tmp_117 * tmp_29 - tmp_117 * tmp_76 ) +
                     tmp_95 * ( -tmp_115 * tmp_93 + tmp_118 * tmp_29 - tmp_118 * tmp_94 );
      real_t a_0_2 = tmp_113 * ( -tmp_111 * tmp_121 - tmp_112 * tmp_125 + tmp_125 * tmp_29 ) +
                     tmp_41 * ( tmp_120 * tmp_29 - tmp_120 * tmp_40 - tmp_121 * tmp_36 ) +
                     tmp_59 * ( -tmp_121 * tmp_57 + tmp_122 * tmp_29 - tmp_122 * tmp_58 ) +
                     tmp_77 * ( -tmp_121 * tmp_75 + tmp_123 * tmp_29 - tmp_123 * tmp_76 ) +
                     tmp_95 * ( -tmp_121 * tmp_93 + tmp_124 * tmp_29 - tmp_124 * tmp_94 );
      real_t a_1_0 = tmp_113 * ( tmp_104 * tmp_127 - tmp_104 * tmp_136 - tmp_135 * tmp_37 ) +
                     tmp_41 * ( -tmp_126 * tmp_37 + tmp_127 * tmp_19 - tmp_128 * tmp_19 ) +
                     tmp_59 * ( tmp_127 * tmp_50 - tmp_129 * tmp_37 - tmp_130 * tmp_50 ) +
                     tmp_77 * ( tmp_127 * tmp_68 - tmp_131 * tmp_37 - tmp_132 * tmp_68 ) +
                     tmp_95 * ( tmp_127 * tmp_86 - tmp_133 * tmp_37 - tmp_134 * tmp_86 );
      real_t a_1_1 = tmp_113 * ( -tmp_115 * tmp_135 + tmp_119 * tmp_127 - tmp_119 * tmp_136 ) +
                     tmp_41 * ( tmp_114 * tmp_127 - tmp_114 * tmp_128 - tmp_115 * tmp_126 ) +
                     tmp_59 * ( -tmp_115 * tmp_129 + tmp_116 * tmp_127 - tmp_116 * tmp_130 ) +
                     tmp_77 * ( -tmp_115 * tmp_131 + tmp_117 * tmp_127 - tmp_117 * tmp_132 ) +
                     tmp_95 * ( -tmp_115 * tmp_133 + tmp_118 * tmp_127 - tmp_118 * tmp_134 );
      real_t a_1_2 = tmp_113 * ( -tmp_121 * tmp_135 + tmp_125 * tmp_127 - tmp_125 * tmp_136 ) +
                     tmp_41 * ( tmp_120 * tmp_127 - tmp_120 * tmp_128 - tmp_121 * tmp_126 ) +
                     tmp_59 * ( -tmp_121 * tmp_129 + tmp_122 * tmp_127 - tmp_122 * tmp_130 ) +
                     tmp_77 * ( -tmp_121 * tmp_131 + tmp_123 * tmp_127 - tmp_123 * tmp_132 ) +
                     tmp_95 * ( -tmp_121 * tmp_133 + tmp_124 * tmp_127 - tmp_124 * tmp_134 );
      real_t a_2_0 = tmp_113 * ( tmp_104 * tmp_138 - tmp_104 * tmp_147 - tmp_146 * tmp_37 ) +
                     tmp_41 * ( -tmp_137 * tmp_37 + tmp_138 * tmp_19 - tmp_139 * tmp_19 ) +
                     tmp_59 * ( tmp_138 * tmp_50 - tmp_140 * tmp_37 - tmp_141 * tmp_50 ) +
                     tmp_77 * ( tmp_138 * tmp_68 - tmp_142 * tmp_37 - tmp_143 * tmp_68 ) +
                     tmp_95 * ( tmp_138 * tmp_86 - tmp_144 * tmp_37 - tmp_145 * tmp_86 );
      real_t a_2_1 = tmp_113 * ( -tmp_115 * tmp_146 + tmp_119 * tmp_138 - tmp_119 * tmp_147 ) +
                     tmp_41 * ( tmp_114 * tmp_138 - tmp_114 * tmp_139 - tmp_115 * tmp_137 ) +
                     tmp_59 * ( -tmp_115 * tmp_140 + tmp_116 * tmp_138 - tmp_116 * tmp_141 ) +
                     tmp_77 * ( -tmp_115 * tmp_142 + tmp_117 * tmp_138 - tmp_117 * tmp_143 ) +
                     tmp_95 * ( -tmp_115 * tmp_144 + tmp_118 * tmp_138 - tmp_118 * tmp_145 );
      real_t a_2_2 = tmp_113 * ( -tmp_121 * tmp_146 + tmp_125 * tmp_138 - tmp_125 * tmp_147 ) +
                     tmp_41 * ( tmp_120 * tmp_138 - tmp_120 * tmp_139 - tmp_121 * tmp_137 ) +
                     tmp_59 * ( -tmp_121 * tmp_140 + tmp_122 * tmp_138 - tmp_122 * tmp_141 ) +
                     tmp_77 * ( -tmp_121 * tmp_142 + tmp_123 * tmp_138 - tmp_123 * tmp_143 ) +
                     tmp_95 * ( -tmp_121 * tmp_144 + tmp_124 * tmp_138 - tmp_124 * tmp_145 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;
      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
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

      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_2_1 + tmp_0;
      real_t tmp_2  = -p_affine_0_0;
      real_t tmp_3  = p_affine_1_0 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_0 ) * ( p_affine_2_0 + tmp_2 ) );
      real_t tmp_5  = tmp_1 * tmp_4;
      real_t tmp_6  = tmp_4 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_7  = tmp_3 * tmp_4;
      real_t tmp_8  = tmp_4 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_9  = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_10 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_11 = std::abs( std::pow( ( tmp_10 * tmp_10 ) + ( tmp_9 * tmp_9 ), 1.0 / 2.0 ) );
      real_t tmp_12 = tmp_11 * ( p_affine_10_0 * ( -tmp_5 - tmp_6 ) + p_affine_10_1 * ( -tmp_7 - tmp_8 ) );
      real_t tmp_13 = p_affine_6_1 + tmp_0;
      real_t tmp_14 = 0.046910077030668018 * tmp_10 + tmp_13;
      real_t tmp_15 = tmp_14 * tmp_7;
      real_t tmp_16 = tmp_14 * tmp_8;
      real_t tmp_17 = p_affine_6_0 + tmp_2;
      real_t tmp_18 = tmp_17 + 0.046910077030668018 * tmp_9;
      real_t tmp_19 = tmp_18 * tmp_5;
      real_t tmp_20 = tmp_18 * tmp_6;
      real_t tmp_21 = -0.11846344252809471 * tmp_15 - 0.11846344252809471 * tmp_16 - 0.11846344252809471 * tmp_19 -
                      0.11846344252809471 * tmp_20 + 0.11846344252809471;
      real_t tmp_22 = 0.23076534494715845 * tmp_10 + tmp_13;
      real_t tmp_23 = tmp_22 * tmp_7;
      real_t tmp_24 = tmp_22 * tmp_8;
      real_t tmp_25 = tmp_17 + 0.23076534494715845 * tmp_9;
      real_t tmp_26 = tmp_25 * tmp_5;
      real_t tmp_27 = tmp_25 * tmp_6;
      real_t tmp_28 = -0.2393143352496831 * tmp_23 - 0.2393143352496831 * tmp_24 - 0.2393143352496831 * tmp_26 -
                      0.2393143352496831 * tmp_27 + 0.2393143352496831;
      real_t tmp_29 = 0.5 * tmp_10 + tmp_13;
      real_t tmp_30 = tmp_29 * tmp_7;
      real_t tmp_31 = tmp_29 * tmp_8;
      real_t tmp_32 = tmp_17 + 0.5 * tmp_9;
      real_t tmp_33 = tmp_32 * tmp_5;
      real_t tmp_34 = tmp_32 * tmp_6;
      real_t tmp_35 = -0.2844444444444445 * tmp_30 - 0.2844444444444445 * tmp_31 - 0.2844444444444445 * tmp_33 -
                      0.2844444444444445 * tmp_34 + 0.2844444444444445;
      real_t tmp_36 = 0.7692346550528415 * tmp_10 + tmp_13;
      real_t tmp_37 = tmp_36 * tmp_7;
      real_t tmp_38 = tmp_36 * tmp_8;
      real_t tmp_39 = tmp_17 + 0.7692346550528415 * tmp_9;
      real_t tmp_40 = tmp_39 * tmp_5;
      real_t tmp_41 = tmp_39 * tmp_6;
      real_t tmp_42 = -0.2393143352496831 * tmp_37 - 0.2393143352496831 * tmp_38 - 0.2393143352496831 * tmp_40 -
                      0.2393143352496831 * tmp_41 + 0.2393143352496831;
      real_t tmp_43 = 0.95308992296933193 * tmp_10 + tmp_13;
      real_t tmp_44 = tmp_43 * tmp_7;
      real_t tmp_45 = tmp_43 * tmp_8;
      real_t tmp_46 = tmp_17 + 0.95308992296933193 * tmp_9;
      real_t tmp_47 = tmp_46 * tmp_5;
      real_t tmp_48 = tmp_46 * tmp_6;
      real_t tmp_49 = -0.11846344252809471 * tmp_44 - 0.11846344252809471 * tmp_45 - 0.11846344252809471 * tmp_47 -
                      0.11846344252809471 * tmp_48 + 0.11846344252809471;
      real_t tmp_50 = tmp_11 * ( p_affine_10_0 * tmp_5 + p_affine_10_1 * tmp_8 );
      real_t tmp_51 = tmp_11 * ( p_affine_10_0 * tmp_6 + p_affine_10_1 * tmp_7 );
      real_t tmp_52 = 0.11846344252809471 * tmp_16 + 0.11846344252809471 * tmp_19;
      real_t tmp_53 = 0.2393143352496831 * tmp_24 + 0.2393143352496831 * tmp_26;
      real_t tmp_54 = 0.2844444444444445 * tmp_31 + 0.2844444444444445 * tmp_33;
      real_t tmp_55 = 0.2393143352496831 * tmp_38 + 0.2393143352496831 * tmp_40;
      real_t tmp_56 = 0.11846344252809471 * tmp_45 + 0.11846344252809471 * tmp_47;
      real_t tmp_57 = 0.11846344252809471 * tmp_15 + 0.11846344252809471 * tmp_20;
      real_t tmp_58 = 0.2393143352496831 * tmp_23 + 0.2393143352496831 * tmp_27;
      real_t tmp_59 = 0.2844444444444445 * tmp_30 + 0.2844444444444445 * tmp_34;
      real_t tmp_60 = 0.2393143352496831 * tmp_37 + 0.2393143352496831 * tmp_41;
      real_t tmp_61 = 0.11846344252809471 * tmp_44 + 0.11846344252809471 * tmp_48;
      real_t a_0_0  = -tmp_12 * tmp_21 - tmp_12 * tmp_28 - tmp_12 * tmp_35 - tmp_12 * tmp_42 - tmp_12 * tmp_49;
      real_t a_0_1  = -tmp_21 * tmp_50 - tmp_28 * tmp_50 - tmp_35 * tmp_50 - tmp_42 * tmp_50 - tmp_49 * tmp_50;
      real_t a_0_2  = -tmp_21 * tmp_51 - tmp_28 * tmp_51 - tmp_35 * tmp_51 - tmp_42 * tmp_51 - tmp_49 * tmp_51;
      real_t a_1_0  = -tmp_12 * tmp_52 - tmp_12 * tmp_53 - tmp_12 * tmp_54 - tmp_12 * tmp_55 - tmp_12 * tmp_56;
      real_t a_1_1  = -tmp_50 * tmp_52 - tmp_50 * tmp_53 - tmp_50 * tmp_54 - tmp_50 * tmp_55 - tmp_50 * tmp_56;
      real_t a_1_2  = -tmp_51 * tmp_52 - tmp_51 * tmp_53 - tmp_51 * tmp_54 - tmp_51 * tmp_55 - tmp_51 * tmp_56;
      real_t a_2_0  = -tmp_12 * tmp_57 - tmp_12 * tmp_58 - tmp_12 * tmp_59 - tmp_12 * tmp_60 - tmp_12 * tmp_61;
      real_t a_2_1  = -tmp_50 * tmp_57 - tmp_50 * tmp_58 - tmp_50 * tmp_59 - tmp_50 * tmp_60 - tmp_50 * tmp_61;
      real_t a_2_2  = -tmp_51 * tmp_57 - tmp_51 * tmp_58 - tmp_51 * tmp_59 - tmp_51 * tmp_60 - tmp_51 * tmp_61;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;
      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
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

class DGVectorLaplaceFormP1P1_10 : public hyteg::dg::DGForm2D
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

      real_t a_0_0  = 0;
      real_t a_0_1  = 0;
      real_t a_0_2  = 0;
      real_t a_1_0  = 0;
      real_t a_1_1  = 0;
      real_t a_1_2  = 0;
      real_t a_2_0  = 0;
      real_t a_2_1  = 0;
      real_t a_2_2  = 0;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;
      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
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

      real_t a_0_0  = 0;
      real_t a_0_1  = 0;
      real_t a_0_2  = 0;
      real_t a_1_0  = 0;
      real_t a_1_1  = 0;
      real_t a_1_2  = 0;
      real_t a_2_0  = 0;
      real_t a_2_1  = 0;
      real_t a_2_2  = 0;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;
      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
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

      real_t a_0_0  = 0;
      real_t a_0_1  = 0;
      real_t a_0_2  = 0;
      real_t a_1_0  = 0;
      real_t a_1_1  = 0;
      real_t a_1_2  = 0;
      real_t a_2_0  = 0;
      real_t a_2_1  = 0;
      real_t a_2_2  = 0;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;
      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
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

      real_t a_0_0  = 0;
      real_t a_0_1  = 0;
      real_t a_0_2  = 0;
      real_t a_1_0  = 0;
      real_t a_1_1  = 0;
      real_t a_1_2  = 0;
      real_t a_2_0  = 0;
      real_t a_2_1  = 0;
      real_t a_2_2  = 0;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;
      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
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

class DGVectorLaplaceFormP1P1_01 : public hyteg::dg::DGForm2D
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

      real_t a_0_0  = 0;
      real_t a_0_1  = 0;
      real_t a_0_2  = 0;
      real_t a_1_0  = 0;
      real_t a_1_1  = 0;
      real_t a_1_2  = 0;
      real_t a_2_0  = 0;
      real_t a_2_1  = 0;
      real_t a_2_2  = 0;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;
      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
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

      real_t a_0_0  = 0;
      real_t a_0_1  = 0;
      real_t a_0_2  = 0;
      real_t a_1_0  = 0;
      real_t a_1_1  = 0;
      real_t a_1_2  = 0;
      real_t a_2_0  = 0;
      real_t a_2_1  = 0;
      real_t a_2_2  = 0;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;
      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
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

      real_t a_0_0  = 0;
      real_t a_0_1  = 0;
      real_t a_0_2  = 0;
      real_t a_1_0  = 0;
      real_t a_1_1  = 0;
      real_t a_1_2  = 0;
      real_t a_2_0  = 0;
      real_t a_2_1  = 0;
      real_t a_2_2  = 0;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;
      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
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

      real_t a_0_0  = 0;
      real_t a_0_1  = 0;
      real_t a_0_2  = 0;
      real_t a_1_0  = 0;
      real_t a_1_1  = 0;
      real_t a_1_2  = 0;
      real_t a_2_0  = 0;
      real_t a_2_1  = 0;
      real_t a_2_2  = 0;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;
      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
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

class DGVectorLaplaceFormP1P1_11 : public hyteg::dg::DGForm2D
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
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 );
      real_t tmp_5  = 1.0 / ( tmp_4 );
      real_t tmp_6  = tmp_1 * tmp_5;
      real_t tmp_7  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_8  = tmp_5 * tmp_7;
      real_t tmp_9  = -tmp_6 - tmp_8;
      real_t tmp_10 = tmp_3 * tmp_5;
      real_t tmp_11 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_12 = tmp_11 * tmp_5;
      real_t tmp_13 = -tmp_10 - tmp_12;
      real_t tmp_14 = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                                p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_15 = tmp_14 * ( ( tmp_13 * tmp_13 ) + ( tmp_9 * tmp_9 ) );
      real_t tmp_16 = tmp_14 * ( tmp_10 * tmp_13 + tmp_8 * tmp_9 );
      real_t tmp_17 = 0.5 * tmp_16;
      real_t tmp_18 = tmp_14 * ( tmp_12 * tmp_13 + tmp_6 * tmp_9 );
      real_t tmp_19 = 0.5 * tmp_18;
      real_t tmp_20 = 1.0 / ( tmp_4 * tmp_4 );
      real_t tmp_21 = tmp_14 * ( tmp_20 * ( tmp_3 * tmp_3 ) + tmp_20 * ( tmp_7 * tmp_7 ) );
      real_t tmp_22 = tmp_14 * ( tmp_1 * tmp_20 * tmp_7 + tmp_11 * tmp_20 * tmp_3 );
      real_t tmp_23 = 0.5 * tmp_22;
      real_t tmp_24 = tmp_14 * ( ( tmp_1 * tmp_1 ) * tmp_20 + ( tmp_11 * tmp_11 ) * tmp_20 );
      real_t a_0_0  = 0.5 * tmp_15;
      real_t a_0_1  = tmp_17;
      real_t a_0_2  = tmp_19;
      real_t a_1_0  = tmp_17;
      real_t a_1_1  = 0.5 * tmp_21;
      real_t a_1_2  = tmp_23;
      real_t a_2_0  = tmp_19;
      real_t a_2_1  = tmp_23;
      real_t a_2_2  = 0.5 * tmp_24;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;
      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
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

      real_t tmp_0  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_1  = -p_affine_0_1;
      real_t tmp_2  = p_affine_6_1 + tmp_1;
      real_t tmp_3  = 0.046910077030668018 * tmp_0 + tmp_2;
      real_t tmp_4  = -p_affine_0_0;
      real_t tmp_5  = p_affine_1_0 + tmp_4;
      real_t tmp_6  = p_affine_2_1 + tmp_1;
      real_t tmp_7  = 1.0 / ( tmp_5 * tmp_6 - ( p_affine_1_1 + tmp_1 ) * ( p_affine_2_0 + tmp_4 ) );
      real_t tmp_8  = tmp_5 * tmp_7;
      real_t tmp_9  = tmp_3 * tmp_8;
      real_t tmp_10 = tmp_7 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_11 = tmp_10 * tmp_3;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_4;
      real_t tmp_14 = 0.046910077030668018 * tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6 * tmp_7;
      real_t tmp_16 = tmp_14 * tmp_15;
      real_t tmp_17 = tmp_7 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_18 = tmp_14 * tmp_17;
      real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20 = std::abs( std::pow( ( tmp_0 * tmp_0 ) + ( tmp_12 * tmp_12 ), 1.0 / 2.0 ) );
      real_t tmp_21 = 6 / tmp_20;
      real_t tmp_22 = p_affine_10_0 * ( -tmp_15 - tmp_17 ) + p_affine_10_1 * ( -tmp_10 - tmp_8 );
      real_t tmp_23 = 1.0 * tmp_22;
      real_t tmp_24 = 0.11846344252809471 * tmp_20;
      real_t tmp_25 = 0.23076534494715845 * tmp_0 + tmp_2;
      real_t tmp_26 = tmp_25 * tmp_8;
      real_t tmp_27 = tmp_10 * tmp_25;
      real_t tmp_28 = 0.23076534494715845 * tmp_12 + tmp_13;
      real_t tmp_29 = tmp_15 * tmp_28;
      real_t tmp_30 = tmp_17 * tmp_28;
      real_t tmp_31 = -tmp_26 - tmp_27 - tmp_29 - tmp_30 + 1;
      real_t tmp_32 = 0.2393143352496831 * tmp_20;
      real_t tmp_33 = 0.5 * tmp_0 + tmp_2;
      real_t tmp_34 = tmp_33 * tmp_8;
      real_t tmp_35 = tmp_10 * tmp_33;
      real_t tmp_36 = 0.5 * tmp_12 + tmp_13;
      real_t tmp_37 = tmp_15 * tmp_36;
      real_t tmp_38 = tmp_17 * tmp_36;
      real_t tmp_39 = -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1;
      real_t tmp_40 = 0.2844444444444445 * tmp_20;
      real_t tmp_41 = 0.7692346550528415 * tmp_0 + tmp_2;
      real_t tmp_42 = tmp_41 * tmp_8;
      real_t tmp_43 = tmp_10 * tmp_41;
      real_t tmp_44 = 0.7692346550528415 * tmp_12 + tmp_13;
      real_t tmp_45 = tmp_15 * tmp_44;
      real_t tmp_46 = tmp_17 * tmp_44;
      real_t tmp_47 = -tmp_42 - tmp_43 - tmp_45 - tmp_46 + 1;
      real_t tmp_48 = 0.2393143352496831 * tmp_20;
      real_t tmp_49 = 0.95308992296933193 * tmp_0 + tmp_2;
      real_t tmp_50 = tmp_49 * tmp_8;
      real_t tmp_51 = tmp_10 * tmp_49;
      real_t tmp_52 = 0.95308992296933193 * tmp_12 + tmp_13;
      real_t tmp_53 = tmp_15 * tmp_52;
      real_t tmp_54 = tmp_17 * tmp_52;
      real_t tmp_55 = -tmp_50 - tmp_51 - tmp_53 - tmp_54 + 1;
      real_t tmp_56 = 0.11846344252809471 * tmp_20;
      real_t tmp_57 = tmp_11 + tmp_16;
      real_t tmp_58 = 0.5 * tmp_22;
      real_t tmp_59 = p_affine_10_0 * tmp_15 + p_affine_10_1 * tmp_10;
      real_t tmp_60 = 0.5 * tmp_59;
      real_t tmp_61 = tmp_19 * tmp_21;
      real_t tmp_62 = tmp_27 + tmp_29;
      real_t tmp_63 = tmp_21 * tmp_31;
      real_t tmp_64 = tmp_35 + tmp_37;
      real_t tmp_65 = tmp_21 * tmp_39;
      real_t tmp_66 = tmp_43 + tmp_45;
      real_t tmp_67 = tmp_21 * tmp_47;
      real_t tmp_68 = tmp_51 + tmp_53;
      real_t tmp_69 = tmp_21 * tmp_55;
      real_t tmp_70 = tmp_24 * ( -tmp_19 * tmp_60 - tmp_57 * tmp_58 + tmp_57 * tmp_61 ) +
                      tmp_32 * ( -tmp_31 * tmp_60 - tmp_58 * tmp_62 + tmp_62 * tmp_63 ) +
                      tmp_40 * ( -tmp_39 * tmp_60 - tmp_58 * tmp_64 + tmp_64 * tmp_65 ) +
                      tmp_48 * ( -tmp_47 * tmp_60 - tmp_58 * tmp_66 + tmp_66 * tmp_67 ) +
                      tmp_56 * ( -tmp_55 * tmp_60 - tmp_58 * tmp_68 + tmp_68 * tmp_69 );
      real_t tmp_71 = tmp_18 + tmp_9;
      real_t tmp_72 = p_affine_10_0 * tmp_17 + p_affine_10_1 * tmp_8;
      real_t tmp_73 = 0.5 * tmp_72;
      real_t tmp_74 = tmp_26 + tmp_30;
      real_t tmp_75 = tmp_34 + tmp_38;
      real_t tmp_76 = tmp_42 + tmp_46;
      real_t tmp_77 = tmp_50 + tmp_54;
      real_t tmp_78 = tmp_24 * ( -tmp_19 * tmp_73 - tmp_58 * tmp_71 + tmp_61 * tmp_71 ) +
                      tmp_32 * ( -tmp_31 * tmp_73 - tmp_58 * tmp_74 + tmp_63 * tmp_74 ) +
                      tmp_40 * ( -tmp_39 * tmp_73 - tmp_58 * tmp_75 + tmp_65 * tmp_75 ) +
                      tmp_48 * ( -tmp_47 * tmp_73 - tmp_58 * tmp_76 + tmp_67 * tmp_76 ) +
                      tmp_56 * ( -tmp_55 * tmp_73 - tmp_58 * tmp_77 + tmp_69 * tmp_77 );
      real_t tmp_79 = 1.0 * tmp_59;
      real_t tmp_80 = tmp_24 * ( tmp_21 * tmp_57 * tmp_71 - tmp_57 * tmp_73 - tmp_60 * tmp_71 ) +
                      tmp_32 * ( tmp_21 * tmp_62 * tmp_74 - tmp_60 * tmp_74 - tmp_62 * tmp_73 ) +
                      tmp_40 * ( tmp_21 * tmp_64 * tmp_75 - tmp_60 * tmp_75 - tmp_64 * tmp_73 ) +
                      tmp_48 * ( tmp_21 * tmp_66 * tmp_76 - tmp_60 * tmp_76 - tmp_66 * tmp_73 ) +
                      tmp_56 * ( tmp_21 * tmp_68 * tmp_77 - tmp_60 * tmp_77 - tmp_68 * tmp_73 );
      real_t tmp_81 = 1.0 * tmp_72;
      real_t a_0_0  = tmp_24 * ( ( tmp_19 * tmp_19 ) * tmp_21 - tmp_19 * tmp_23 ) +
                     tmp_32 * ( tmp_21 * ( tmp_31 * tmp_31 ) - tmp_23 * tmp_31 ) +
                     tmp_40 * ( tmp_21 * ( tmp_39 * tmp_39 ) - tmp_23 * tmp_39 ) +
                     tmp_48 * ( tmp_21 * ( tmp_47 * tmp_47 ) - tmp_23 * tmp_47 ) +
                     tmp_56 * ( tmp_21 * ( tmp_55 * tmp_55 ) - tmp_23 * tmp_55 );
      real_t a_0_1 = tmp_70;
      real_t a_0_2 = tmp_78;
      real_t a_1_0 = tmp_70;
      real_t a_1_1 = tmp_24 * ( tmp_21 * ( tmp_57 * tmp_57 ) - tmp_57 * tmp_79 ) +
                     tmp_32 * ( tmp_21 * ( tmp_62 * tmp_62 ) - tmp_62 * tmp_79 ) +
                     tmp_40 * ( tmp_21 * ( tmp_64 * tmp_64 ) - tmp_64 * tmp_79 ) +
                     tmp_48 * ( tmp_21 * ( tmp_66 * tmp_66 ) - tmp_66 * tmp_79 ) +
                     tmp_56 * ( tmp_21 * ( tmp_68 * tmp_68 ) - tmp_68 * tmp_79 );
      real_t a_1_2 = tmp_80;
      real_t a_2_0 = tmp_78;
      real_t a_2_1 = tmp_80;
      real_t a_2_2 = tmp_24 * ( tmp_21 * ( tmp_71 * tmp_71 ) - tmp_71 * tmp_81 ) +
                     tmp_32 * ( tmp_21 * ( tmp_74 * tmp_74 ) - tmp_74 * tmp_81 ) +
                     tmp_40 * ( tmp_21 * ( tmp_75 * tmp_75 ) - tmp_75 * tmp_81 ) +
                     tmp_48 * ( tmp_21 * ( tmp_76 * tmp_76 ) - tmp_76 * tmp_81 ) +
                     tmp_56 * ( tmp_21 * ( tmp_77 * tmp_77 ) - tmp_77 * tmp_81 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;
      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
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

      real_t tmp_0   = -p_affine_3_1;
      real_t tmp_1   = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2   = p_affine_6_1 + 0.046910077030668018 * tmp_1;
      real_t tmp_3   = tmp_0 + tmp_2;
      real_t tmp_4   = -p_affine_3_0;
      real_t tmp_5   = p_affine_4_0 + tmp_4;
      real_t tmp_6   = p_affine_5_1 + tmp_0;
      real_t tmp_7   = 1.0 / ( tmp_5 * tmp_6 - ( p_affine_4_1 + tmp_0 ) * ( p_affine_5_0 + tmp_4 ) );
      real_t tmp_8   = tmp_5 * tmp_7;
      real_t tmp_9   = tmp_3 * tmp_8;
      real_t tmp_10  = tmp_7 * ( p_affine_3_0 - p_affine_5_0 );
      real_t tmp_11  = tmp_10 * tmp_3;
      real_t tmp_12  = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13  = p_affine_6_0 + 0.046910077030668018 * tmp_12;
      real_t tmp_14  = tmp_13 + tmp_4;
      real_t tmp_15  = tmp_6 * tmp_7;
      real_t tmp_16  = tmp_14 * tmp_15;
      real_t tmp_17  = tmp_7 * ( p_affine_3_1 - p_affine_4_1 );
      real_t tmp_18  = tmp_14 * tmp_17;
      real_t tmp_19  = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20  = -p_affine_0_1;
      real_t tmp_21  = p_affine_2_1 + tmp_20;
      real_t tmp_22  = -p_affine_0_0;
      real_t tmp_23  = p_affine_1_0 + tmp_22;
      real_t tmp_24  = 1.0 / ( tmp_21 * tmp_23 - ( p_affine_1_1 + tmp_20 ) * ( p_affine_2_0 + tmp_22 ) );
      real_t tmp_25  = tmp_21 * tmp_24;
      real_t tmp_26  = tmp_24 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_27  = tmp_23 * tmp_24;
      real_t tmp_28  = tmp_24 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_29  = 0.5 * p_affine_10_0 * ( -tmp_25 - tmp_26 ) + 0.5 * p_affine_10_1 * ( -tmp_27 - tmp_28 );
      real_t tmp_30  = tmp_2 + tmp_20;
      real_t tmp_31  = tmp_27 * tmp_30;
      real_t tmp_32  = tmp_28 * tmp_30;
      real_t tmp_33  = tmp_13 + tmp_22;
      real_t tmp_34  = tmp_25 * tmp_33;
      real_t tmp_35  = tmp_26 * tmp_33;
      real_t tmp_36  = -tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1;
      real_t tmp_37  = 0.5 * p_affine_10_0 * ( -tmp_15 - tmp_17 ) + 0.5 * p_affine_10_1 * ( -tmp_10 - tmp_8 );
      real_t tmp_38  = std::abs( std::pow( ( tmp_1 * tmp_1 ) + ( tmp_12 * tmp_12 ), 1.0 / 2.0 ) );
      real_t tmp_39  = 6 / tmp_38;
      real_t tmp_40  = tmp_36 * tmp_39;
      real_t tmp_41  = 0.11846344252809471 * tmp_38;
      real_t tmp_42  = p_affine_6_1 + 0.23076534494715845 * tmp_1;
      real_t tmp_43  = tmp_0 + tmp_42;
      real_t tmp_44  = tmp_43 * tmp_8;
      real_t tmp_45  = tmp_10 * tmp_43;
      real_t tmp_46  = p_affine_6_0 + 0.23076534494715845 * tmp_12;
      real_t tmp_47  = tmp_4 + tmp_46;
      real_t tmp_48  = tmp_15 * tmp_47;
      real_t tmp_49  = tmp_17 * tmp_47;
      real_t tmp_50  = -tmp_44 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_51  = tmp_20 + tmp_42;
      real_t tmp_52  = tmp_27 * tmp_51;
      real_t tmp_53  = tmp_28 * tmp_51;
      real_t tmp_54  = tmp_22 + tmp_46;
      real_t tmp_55  = tmp_25 * tmp_54;
      real_t tmp_56  = tmp_26 * tmp_54;
      real_t tmp_57  = -tmp_52 - tmp_53 - tmp_55 - tmp_56 + 1;
      real_t tmp_58  = tmp_39 * tmp_57;
      real_t tmp_59  = 0.2393143352496831 * tmp_38;
      real_t tmp_60  = p_affine_6_1 + 0.5 * tmp_1;
      real_t tmp_61  = tmp_0 + tmp_60;
      real_t tmp_62  = tmp_61 * tmp_8;
      real_t tmp_63  = tmp_10 * tmp_61;
      real_t tmp_64  = p_affine_6_0 + 0.5 * tmp_12;
      real_t tmp_65  = tmp_4 + tmp_64;
      real_t tmp_66  = tmp_15 * tmp_65;
      real_t tmp_67  = tmp_17 * tmp_65;
      real_t tmp_68  = -tmp_62 - tmp_63 - tmp_66 - tmp_67 + 1;
      real_t tmp_69  = tmp_20 + tmp_60;
      real_t tmp_70  = tmp_27 * tmp_69;
      real_t tmp_71  = tmp_28 * tmp_69;
      real_t tmp_72  = tmp_22 + tmp_64;
      real_t tmp_73  = tmp_25 * tmp_72;
      real_t tmp_74  = tmp_26 * tmp_72;
      real_t tmp_75  = -tmp_70 - tmp_71 - tmp_73 - tmp_74 + 1;
      real_t tmp_76  = tmp_39 * tmp_75;
      real_t tmp_77  = 0.2844444444444445 * tmp_38;
      real_t tmp_78  = p_affine_6_1 + 0.7692346550528415 * tmp_1;
      real_t tmp_79  = tmp_0 + tmp_78;
      real_t tmp_80  = tmp_79 * tmp_8;
      real_t tmp_81  = tmp_10 * tmp_79;
      real_t tmp_82  = p_affine_6_0 + 0.7692346550528415 * tmp_12;
      real_t tmp_83  = tmp_4 + tmp_82;
      real_t tmp_84  = tmp_15 * tmp_83;
      real_t tmp_85  = tmp_17 * tmp_83;
      real_t tmp_86  = -tmp_80 - tmp_81 - tmp_84 - tmp_85 + 1;
      real_t tmp_87  = tmp_20 + tmp_78;
      real_t tmp_88  = tmp_27 * tmp_87;
      real_t tmp_89  = tmp_28 * tmp_87;
      real_t tmp_90  = tmp_22 + tmp_82;
      real_t tmp_91  = tmp_25 * tmp_90;
      real_t tmp_92  = tmp_26 * tmp_90;
      real_t tmp_93  = -tmp_88 - tmp_89 - tmp_91 - tmp_92 + 1;
      real_t tmp_94  = tmp_39 * tmp_93;
      real_t tmp_95  = 0.2393143352496831 * tmp_38;
      real_t tmp_96  = p_affine_6_1 + 0.95308992296933193 * tmp_1;
      real_t tmp_97  = tmp_0 + tmp_96;
      real_t tmp_98  = tmp_8 * tmp_97;
      real_t tmp_99  = tmp_10 * tmp_97;
      real_t tmp_100 = p_affine_6_0 + 0.95308992296933193 * tmp_12;
      real_t tmp_101 = tmp_100 + tmp_4;
      real_t tmp_102 = tmp_101 * tmp_15;
      real_t tmp_103 = tmp_101 * tmp_17;
      real_t tmp_104 = -tmp_102 - tmp_103 - tmp_98 - tmp_99 + 1;
      real_t tmp_105 = tmp_20 + tmp_96;
      real_t tmp_106 = tmp_105 * tmp_27;
      real_t tmp_107 = tmp_105 * tmp_28;
      real_t tmp_108 = tmp_100 + tmp_22;
      real_t tmp_109 = tmp_108 * tmp_25;
      real_t tmp_110 = tmp_108 * tmp_26;
      real_t tmp_111 = -tmp_106 - tmp_107 - tmp_109 - tmp_110 + 1;
      real_t tmp_112 = tmp_111 * tmp_39;
      real_t tmp_113 = 0.11846344252809471 * tmp_38;
      real_t tmp_114 = tmp_11 + tmp_16;
      real_t tmp_115 = 0.5 * p_affine_10_0 * tmp_15 + 0.5 * p_affine_10_1 * tmp_10;
      real_t tmp_116 = tmp_45 + tmp_48;
      real_t tmp_117 = tmp_63 + tmp_66;
      real_t tmp_118 = tmp_81 + tmp_84;
      real_t tmp_119 = tmp_102 + tmp_99;
      real_t tmp_120 = tmp_18 + tmp_9;
      real_t tmp_121 = 0.5 * p_affine_10_0 * tmp_17 + 0.5 * p_affine_10_1 * tmp_8;
      real_t tmp_122 = tmp_44 + tmp_49;
      real_t tmp_123 = tmp_62 + tmp_67;
      real_t tmp_124 = tmp_80 + tmp_85;
      real_t tmp_125 = tmp_103 + tmp_98;
      real_t tmp_126 = tmp_32 + tmp_34;
      real_t tmp_127 = 0.5 * p_affine_10_0 * tmp_25 + 0.5 * p_affine_10_1 * tmp_28;
      real_t tmp_128 = tmp_126 * tmp_39;
      real_t tmp_129 = tmp_53 + tmp_55;
      real_t tmp_130 = tmp_129 * tmp_39;
      real_t tmp_131 = tmp_71 + tmp_73;
      real_t tmp_132 = tmp_131 * tmp_39;
      real_t tmp_133 = tmp_89 + tmp_91;
      real_t tmp_134 = tmp_133 * tmp_39;
      real_t tmp_135 = tmp_107 + tmp_109;
      real_t tmp_136 = tmp_135 * tmp_39;
      real_t tmp_137 = tmp_31 + tmp_35;
      real_t tmp_138 = 0.5 * p_affine_10_0 * tmp_26 + 0.5 * p_affine_10_1 * tmp_27;
      real_t tmp_139 = tmp_137 * tmp_39;
      real_t tmp_140 = tmp_52 + tmp_56;
      real_t tmp_141 = tmp_140 * tmp_39;
      real_t tmp_142 = tmp_70 + tmp_74;
      real_t tmp_143 = tmp_142 * tmp_39;
      real_t tmp_144 = tmp_88 + tmp_92;
      real_t tmp_145 = tmp_144 * tmp_39;
      real_t tmp_146 = tmp_106 + tmp_110;
      real_t tmp_147 = tmp_146 * tmp_39;
      real_t a_0_0   = tmp_113 * ( -tmp_104 * tmp_112 + tmp_104 * tmp_29 - tmp_111 * tmp_37 ) +
                     tmp_41 * ( tmp_19 * tmp_29 - tmp_19 * tmp_40 - tmp_36 * tmp_37 ) +
                     tmp_59 * ( tmp_29 * tmp_50 - tmp_37 * tmp_57 - tmp_50 * tmp_58 ) +
                     tmp_77 * ( tmp_29 * tmp_68 - tmp_37 * tmp_75 - tmp_68 * tmp_76 ) +
                     tmp_95 * ( tmp_29 * tmp_86 - tmp_37 * tmp_93 - tmp_86 * tmp_94 );
      real_t a_0_1 = tmp_113 * ( -tmp_111 * tmp_115 - tmp_112 * tmp_119 + tmp_119 * tmp_29 ) +
                     tmp_41 * ( tmp_114 * tmp_29 - tmp_114 * tmp_40 - tmp_115 * tmp_36 ) +
                     tmp_59 * ( -tmp_115 * tmp_57 + tmp_116 * tmp_29 - tmp_116 * tmp_58 ) +
                     tmp_77 * ( -tmp_115 * tmp_75 + tmp_117 * tmp_29 - tmp_117 * tmp_76 ) +
                     tmp_95 * ( -tmp_115 * tmp_93 + tmp_118 * tmp_29 - tmp_118 * tmp_94 );
      real_t a_0_2 = tmp_113 * ( -tmp_111 * tmp_121 - tmp_112 * tmp_125 + tmp_125 * tmp_29 ) +
                     tmp_41 * ( tmp_120 * tmp_29 - tmp_120 * tmp_40 - tmp_121 * tmp_36 ) +
                     tmp_59 * ( -tmp_121 * tmp_57 + tmp_122 * tmp_29 - tmp_122 * tmp_58 ) +
                     tmp_77 * ( -tmp_121 * tmp_75 + tmp_123 * tmp_29 - tmp_123 * tmp_76 ) +
                     tmp_95 * ( -tmp_121 * tmp_93 + tmp_124 * tmp_29 - tmp_124 * tmp_94 );
      real_t a_1_0 = tmp_113 * ( tmp_104 * tmp_127 - tmp_104 * tmp_136 - tmp_135 * tmp_37 ) +
                     tmp_41 * ( -tmp_126 * tmp_37 + tmp_127 * tmp_19 - tmp_128 * tmp_19 ) +
                     tmp_59 * ( tmp_127 * tmp_50 - tmp_129 * tmp_37 - tmp_130 * tmp_50 ) +
                     tmp_77 * ( tmp_127 * tmp_68 - tmp_131 * tmp_37 - tmp_132 * tmp_68 ) +
                     tmp_95 * ( tmp_127 * tmp_86 - tmp_133 * tmp_37 - tmp_134 * tmp_86 );
      real_t a_1_1 = tmp_113 * ( -tmp_115 * tmp_135 + tmp_119 * tmp_127 - tmp_119 * tmp_136 ) +
                     tmp_41 * ( tmp_114 * tmp_127 - tmp_114 * tmp_128 - tmp_115 * tmp_126 ) +
                     tmp_59 * ( -tmp_115 * tmp_129 + tmp_116 * tmp_127 - tmp_116 * tmp_130 ) +
                     tmp_77 * ( -tmp_115 * tmp_131 + tmp_117 * tmp_127 - tmp_117 * tmp_132 ) +
                     tmp_95 * ( -tmp_115 * tmp_133 + tmp_118 * tmp_127 - tmp_118 * tmp_134 );
      real_t a_1_2 = tmp_113 * ( -tmp_121 * tmp_135 + tmp_125 * tmp_127 - tmp_125 * tmp_136 ) +
                     tmp_41 * ( tmp_120 * tmp_127 - tmp_120 * tmp_128 - tmp_121 * tmp_126 ) +
                     tmp_59 * ( -tmp_121 * tmp_129 + tmp_122 * tmp_127 - tmp_122 * tmp_130 ) +
                     tmp_77 * ( -tmp_121 * tmp_131 + tmp_123 * tmp_127 - tmp_123 * tmp_132 ) +
                     tmp_95 * ( -tmp_121 * tmp_133 + tmp_124 * tmp_127 - tmp_124 * tmp_134 );
      real_t a_2_0 = tmp_113 * ( tmp_104 * tmp_138 - tmp_104 * tmp_147 - tmp_146 * tmp_37 ) +
                     tmp_41 * ( -tmp_137 * tmp_37 + tmp_138 * tmp_19 - tmp_139 * tmp_19 ) +
                     tmp_59 * ( tmp_138 * tmp_50 - tmp_140 * tmp_37 - tmp_141 * tmp_50 ) +
                     tmp_77 * ( tmp_138 * tmp_68 - tmp_142 * tmp_37 - tmp_143 * tmp_68 ) +
                     tmp_95 * ( tmp_138 * tmp_86 - tmp_144 * tmp_37 - tmp_145 * tmp_86 );
      real_t a_2_1 = tmp_113 * ( -tmp_115 * tmp_146 + tmp_119 * tmp_138 - tmp_119 * tmp_147 ) +
                     tmp_41 * ( tmp_114 * tmp_138 - tmp_114 * tmp_139 - tmp_115 * tmp_137 ) +
                     tmp_59 * ( -tmp_115 * tmp_140 + tmp_116 * tmp_138 - tmp_116 * tmp_141 ) +
                     tmp_77 * ( -tmp_115 * tmp_142 + tmp_117 * tmp_138 - tmp_117 * tmp_143 ) +
                     tmp_95 * ( -tmp_115 * tmp_144 + tmp_118 * tmp_138 - tmp_118 * tmp_145 );
      real_t a_2_2 = tmp_113 * ( -tmp_121 * tmp_146 + tmp_125 * tmp_138 - tmp_125 * tmp_147 ) +
                     tmp_41 * ( tmp_120 * tmp_138 - tmp_120 * tmp_139 - tmp_121 * tmp_137 ) +
                     tmp_59 * ( -tmp_121 * tmp_140 + tmp_122 * tmp_138 - tmp_122 * tmp_141 ) +
                     tmp_77 * ( -tmp_121 * tmp_142 + tmp_123 * tmp_138 - tmp_123 * tmp_143 ) +
                     tmp_95 * ( -tmp_121 * tmp_144 + tmp_124 * tmp_138 - tmp_124 * tmp_145 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;
      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
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

      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_2_1 + tmp_0;
      real_t tmp_2  = -p_affine_0_0;
      real_t tmp_3  = p_affine_1_0 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_0 ) * ( p_affine_2_0 + tmp_2 ) );
      real_t tmp_5  = tmp_1 * tmp_4;
      real_t tmp_6  = tmp_4 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_7  = tmp_3 * tmp_4;
      real_t tmp_8  = tmp_4 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_9  = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_10 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_11 = std::abs( std::pow( ( tmp_10 * tmp_10 ) + ( tmp_9 * tmp_9 ), 1.0 / 2.0 ) );
      real_t tmp_12 = tmp_11 * ( p_affine_10_0 * ( -tmp_5 - tmp_6 ) + p_affine_10_1 * ( -tmp_7 - tmp_8 ) );
      real_t tmp_13 = p_affine_6_1 + tmp_0;
      real_t tmp_14 = 0.046910077030668018 * tmp_10 + tmp_13;
      real_t tmp_15 = tmp_14 * tmp_7;
      real_t tmp_16 = tmp_14 * tmp_8;
      real_t tmp_17 = p_affine_6_0 + tmp_2;
      real_t tmp_18 = tmp_17 + 0.046910077030668018 * tmp_9;
      real_t tmp_19 = tmp_18 * tmp_5;
      real_t tmp_20 = tmp_18 * tmp_6;
      real_t tmp_21 = -0.11846344252809471 * tmp_15 - 0.11846344252809471 * tmp_16 - 0.11846344252809471 * tmp_19 -
                      0.11846344252809471 * tmp_20 + 0.11846344252809471;
      real_t tmp_22 = 0.23076534494715845 * tmp_10 + tmp_13;
      real_t tmp_23 = tmp_22 * tmp_7;
      real_t tmp_24 = tmp_22 * tmp_8;
      real_t tmp_25 = tmp_17 + 0.23076534494715845 * tmp_9;
      real_t tmp_26 = tmp_25 * tmp_5;
      real_t tmp_27 = tmp_25 * tmp_6;
      real_t tmp_28 = -0.2393143352496831 * tmp_23 - 0.2393143352496831 * tmp_24 - 0.2393143352496831 * tmp_26 -
                      0.2393143352496831 * tmp_27 + 0.2393143352496831;
      real_t tmp_29 = 0.5 * tmp_10 + tmp_13;
      real_t tmp_30 = tmp_29 * tmp_7;
      real_t tmp_31 = tmp_29 * tmp_8;
      real_t tmp_32 = tmp_17 + 0.5 * tmp_9;
      real_t tmp_33 = tmp_32 * tmp_5;
      real_t tmp_34 = tmp_32 * tmp_6;
      real_t tmp_35 = -0.2844444444444445 * tmp_30 - 0.2844444444444445 * tmp_31 - 0.2844444444444445 * tmp_33 -
                      0.2844444444444445 * tmp_34 + 0.2844444444444445;
      real_t tmp_36 = 0.7692346550528415 * tmp_10 + tmp_13;
      real_t tmp_37 = tmp_36 * tmp_7;
      real_t tmp_38 = tmp_36 * tmp_8;
      real_t tmp_39 = tmp_17 + 0.7692346550528415 * tmp_9;
      real_t tmp_40 = tmp_39 * tmp_5;
      real_t tmp_41 = tmp_39 * tmp_6;
      real_t tmp_42 = -0.2393143352496831 * tmp_37 - 0.2393143352496831 * tmp_38 - 0.2393143352496831 * tmp_40 -
                      0.2393143352496831 * tmp_41 + 0.2393143352496831;
      real_t tmp_43 = 0.95308992296933193 * tmp_10 + tmp_13;
      real_t tmp_44 = tmp_43 * tmp_7;
      real_t tmp_45 = tmp_43 * tmp_8;
      real_t tmp_46 = tmp_17 + 0.95308992296933193 * tmp_9;
      real_t tmp_47 = tmp_46 * tmp_5;
      real_t tmp_48 = tmp_46 * tmp_6;
      real_t tmp_49 = -0.11846344252809471 * tmp_44 - 0.11846344252809471 * tmp_45 - 0.11846344252809471 * tmp_47 -
                      0.11846344252809471 * tmp_48 + 0.11846344252809471;
      real_t tmp_50 = tmp_11 * ( p_affine_10_0 * tmp_5 + p_affine_10_1 * tmp_8 );
      real_t tmp_51 = tmp_11 * ( p_affine_10_0 * tmp_6 + p_affine_10_1 * tmp_7 );
      real_t tmp_52 = 0.11846344252809471 * tmp_16 + 0.11846344252809471 * tmp_19;
      real_t tmp_53 = 0.2393143352496831 * tmp_24 + 0.2393143352496831 * tmp_26;
      real_t tmp_54 = 0.2844444444444445 * tmp_31 + 0.2844444444444445 * tmp_33;
      real_t tmp_55 = 0.2393143352496831 * tmp_38 + 0.2393143352496831 * tmp_40;
      real_t tmp_56 = 0.11846344252809471 * tmp_45 + 0.11846344252809471 * tmp_47;
      real_t tmp_57 = 0.11846344252809471 * tmp_15 + 0.11846344252809471 * tmp_20;
      real_t tmp_58 = 0.2393143352496831 * tmp_23 + 0.2393143352496831 * tmp_27;
      real_t tmp_59 = 0.2844444444444445 * tmp_30 + 0.2844444444444445 * tmp_34;
      real_t tmp_60 = 0.2393143352496831 * tmp_37 + 0.2393143352496831 * tmp_41;
      real_t tmp_61 = 0.11846344252809471 * tmp_44 + 0.11846344252809471 * tmp_48;
      real_t a_0_0  = -tmp_12 * tmp_21 - tmp_12 * tmp_28 - tmp_12 * tmp_35 - tmp_12 * tmp_42 - tmp_12 * tmp_49;
      real_t a_0_1  = -tmp_21 * tmp_50 - tmp_28 * tmp_50 - tmp_35 * tmp_50 - tmp_42 * tmp_50 - tmp_49 * tmp_50;
      real_t a_0_2  = -tmp_21 * tmp_51 - tmp_28 * tmp_51 - tmp_35 * tmp_51 - tmp_42 * tmp_51 - tmp_49 * tmp_51;
      real_t a_1_0  = -tmp_12 * tmp_52 - tmp_12 * tmp_53 - tmp_12 * tmp_54 - tmp_12 * tmp_55 - tmp_12 * tmp_56;
      real_t a_1_1  = -tmp_50 * tmp_52 - tmp_50 * tmp_53 - tmp_50 * tmp_54 - tmp_50 * tmp_55 - tmp_50 * tmp_56;
      real_t a_1_2  = -tmp_51 * tmp_52 - tmp_51 * tmp_53 - tmp_51 * tmp_54 - tmp_51 * tmp_55 - tmp_51 * tmp_56;
      real_t a_2_0  = -tmp_12 * tmp_57 - tmp_12 * tmp_58 - tmp_12 * tmp_59 - tmp_12 * tmp_60 - tmp_12 * tmp_61;
      real_t a_2_1  = -tmp_50 * tmp_57 - tmp_50 * tmp_58 - tmp_50 * tmp_59 - tmp_50 * tmp_60 - tmp_50 * tmp_61;
      real_t a_2_2  = -tmp_51 * tmp_57 - tmp_51 * tmp_58 - tmp_51 * tmp_59 - tmp_51 * tmp_60 - tmp_51 * tmp_61;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 1, 0 ) = a_1_0;
      elMat( 1, 1 ) = a_1_1;
      elMat( 1, 2 ) = a_1_2;
      elMat( 2, 0 ) = a_2_0;
      elMat( 2, 1 ) = a_2_1;
      elMat( 2, 2 ) = a_2_2;
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


} // namespace dg
} // namespace hyteg

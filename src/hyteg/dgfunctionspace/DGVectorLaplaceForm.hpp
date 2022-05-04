
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
      real_t tmp_3  = 0.21132486540518713 * tmp_0 + tmp_2;
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
      real_t tmp_14 = 0.21132486540518713 * tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6 * tmp_7;
      real_t tmp_16 = tmp_14 * tmp_15;
      real_t tmp_17 = tmp_7 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_18 = tmp_14 * tmp_17;
      real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20 = std::abs( std::pow( ( tmp_0 * tmp_0 ) + ( tmp_12 * tmp_12 ), 1.0 / 2.0 ) );
      real_t tmp_21 = 6 / tmp_20;
      real_t tmp_22 = p_affine_10_0 * ( -tmp_15 - tmp_17 ) + p_affine_10_1 * ( -tmp_10 - tmp_8 );
      real_t tmp_23 = 1.0 * tmp_22;
      real_t tmp_24 = 0.5 * tmp_20;
      real_t tmp_25 = 0.78867513459481287 * tmp_0 + tmp_2;
      real_t tmp_26 = tmp_25 * tmp_8;
      real_t tmp_27 = tmp_10 * tmp_25;
      real_t tmp_28 = 0.78867513459481287 * tmp_12 + tmp_13;
      real_t tmp_29 = tmp_15 * tmp_28;
      real_t tmp_30 = tmp_17 * tmp_28;
      real_t tmp_31 = -tmp_26 - tmp_27 - tmp_29 - tmp_30 + 1;
      real_t tmp_32 = 0.5 * tmp_20;
      real_t tmp_33 = tmp_11 + tmp_16;
      real_t tmp_34 = 0.5 * tmp_22;
      real_t tmp_35 = p_affine_10_0 * tmp_15 + p_affine_10_1 * tmp_10;
      real_t tmp_36 = 0.5 * tmp_35;
      real_t tmp_37 = tmp_19 * tmp_21;
      real_t tmp_38 = tmp_27 + tmp_29;
      real_t tmp_39 = tmp_21 * tmp_31;
      real_t tmp_40 = tmp_24 * ( -tmp_19 * tmp_36 - tmp_33 * tmp_34 + tmp_33 * tmp_37 ) +
                      tmp_32 * ( -tmp_31 * tmp_36 - tmp_34 * tmp_38 + tmp_38 * tmp_39 );
      real_t tmp_41 = tmp_18 + tmp_9;
      real_t tmp_42 = p_affine_10_0 * tmp_17 + p_affine_10_1 * tmp_8;
      real_t tmp_43 = 0.5 * tmp_42;
      real_t tmp_44 = tmp_26 + tmp_30;
      real_t tmp_45 = tmp_24 * ( -tmp_19 * tmp_43 - tmp_34 * tmp_41 + tmp_37 * tmp_41 ) +
                      tmp_32 * ( -tmp_31 * tmp_43 - tmp_34 * tmp_44 + tmp_39 * tmp_44 );
      real_t tmp_46 = 1.0 * tmp_35;
      real_t tmp_47 = tmp_24 * ( tmp_21 * tmp_33 * tmp_41 - tmp_33 * tmp_43 - tmp_36 * tmp_41 ) +
                      tmp_32 * ( tmp_21 * tmp_38 * tmp_44 - tmp_36 * tmp_44 - tmp_38 * tmp_43 );
      real_t tmp_48 = 1.0 * tmp_42;
      real_t a_0_0  = tmp_24 * ( ( tmp_19 * tmp_19 ) * tmp_21 - tmp_19 * tmp_23 ) +
                     tmp_32 * ( tmp_21 * ( tmp_31 * tmp_31 ) - tmp_23 * tmp_31 );
      real_t a_0_1 = tmp_40;
      real_t a_0_2 = tmp_45;
      real_t a_1_0 = tmp_40;
      real_t a_1_1 = tmp_24 * ( tmp_21 * ( tmp_33 * tmp_33 ) - tmp_33 * tmp_46 ) +
                     tmp_32 * ( tmp_21 * ( tmp_38 * tmp_38 ) - tmp_38 * tmp_46 );
      real_t a_1_2 = tmp_47;
      real_t a_2_0 = tmp_45;
      real_t a_2_1 = tmp_47;
      real_t a_2_2 = tmp_24 * ( tmp_21 * ( tmp_41 * tmp_41 ) - tmp_41 * tmp_48 ) +
                     tmp_32 * ( tmp_21 * ( tmp_44 * tmp_44 ) - tmp_44 * tmp_48 );
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

      real_t tmp_0  = -p_affine_3_1;
      real_t tmp_1  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2  = p_affine_6_1 + 0.21132486540518713 * tmp_1;
      real_t tmp_3  = tmp_0 + tmp_2;
      real_t tmp_4  = -p_affine_3_0;
      real_t tmp_5  = p_affine_4_0 + tmp_4;
      real_t tmp_6  = p_affine_5_1 + tmp_0;
      real_t tmp_7  = 1.0 / ( tmp_5 * tmp_6 - ( p_affine_4_1 + tmp_0 ) * ( p_affine_5_0 + tmp_4 ) );
      real_t tmp_8  = tmp_5 * tmp_7;
      real_t tmp_9  = tmp_3 * tmp_8;
      real_t tmp_10 = tmp_7 * ( p_affine_3_0 - p_affine_5_0 );
      real_t tmp_11 = tmp_10 * tmp_3;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + 0.21132486540518713 * tmp_12;
      real_t tmp_14 = tmp_13 + tmp_4;
      real_t tmp_15 = tmp_6 * tmp_7;
      real_t tmp_16 = tmp_14 * tmp_15;
      real_t tmp_17 = tmp_7 * ( p_affine_3_1 - p_affine_4_1 );
      real_t tmp_18 = tmp_14 * tmp_17;
      real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20 = -p_affine_0_1;
      real_t tmp_21 = p_affine_2_1 + tmp_20;
      real_t tmp_22 = -p_affine_0_0;
      real_t tmp_23 = p_affine_1_0 + tmp_22;
      real_t tmp_24 = 1.0 / ( tmp_21 * tmp_23 - ( p_affine_1_1 + tmp_20 ) * ( p_affine_2_0 + tmp_22 ) );
      real_t tmp_25 = tmp_21 * tmp_24;
      real_t tmp_26 = tmp_24 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_27 = tmp_23 * tmp_24;
      real_t tmp_28 = tmp_24 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_29 = 0.5 * p_affine_10_0 * ( -tmp_25 - tmp_26 ) + 0.5 * p_affine_10_1 * ( -tmp_27 - tmp_28 );
      real_t tmp_30 = tmp_2 + tmp_20;
      real_t tmp_31 = tmp_27 * tmp_30;
      real_t tmp_32 = tmp_28 * tmp_30;
      real_t tmp_33 = tmp_13 + tmp_22;
      real_t tmp_34 = tmp_25 * tmp_33;
      real_t tmp_35 = tmp_26 * tmp_33;
      real_t tmp_36 = -tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1;
      real_t tmp_37 = 0.5 * p_affine_10_0 * ( -tmp_15 - tmp_17 ) + 0.5 * p_affine_10_1 * ( -tmp_10 - tmp_8 );
      real_t tmp_38 = std::abs( std::pow( ( tmp_1 * tmp_1 ) + ( tmp_12 * tmp_12 ), 1.0 / 2.0 ) );
      real_t tmp_39 = 6 / tmp_38;
      real_t tmp_40 = tmp_36 * tmp_39;
      real_t tmp_41 = 0.5 * tmp_38;
      real_t tmp_42 = p_affine_6_1 + 0.78867513459481287 * tmp_1;
      real_t tmp_43 = tmp_0 + tmp_42;
      real_t tmp_44 = tmp_43 * tmp_8;
      real_t tmp_45 = tmp_10 * tmp_43;
      real_t tmp_46 = p_affine_6_0 + 0.78867513459481287 * tmp_12;
      real_t tmp_47 = tmp_4 + tmp_46;
      real_t tmp_48 = tmp_15 * tmp_47;
      real_t tmp_49 = tmp_17 * tmp_47;
      real_t tmp_50 = -tmp_44 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_51 = tmp_20 + tmp_42;
      real_t tmp_52 = tmp_27 * tmp_51;
      real_t tmp_53 = tmp_28 * tmp_51;
      real_t tmp_54 = tmp_22 + tmp_46;
      real_t tmp_55 = tmp_25 * tmp_54;
      real_t tmp_56 = tmp_26 * tmp_54;
      real_t tmp_57 = -tmp_52 - tmp_53 - tmp_55 - tmp_56 + 1;
      real_t tmp_58 = tmp_39 * tmp_57;
      real_t tmp_59 = 0.5 * tmp_38;
      real_t tmp_60 = tmp_11 + tmp_16;
      real_t tmp_61 = 0.5 * p_affine_10_0 * tmp_15 + 0.5 * p_affine_10_1 * tmp_10;
      real_t tmp_62 = tmp_45 + tmp_48;
      real_t tmp_63 = tmp_18 + tmp_9;
      real_t tmp_64 = 0.5 * p_affine_10_0 * tmp_17 + 0.5 * p_affine_10_1 * tmp_8;
      real_t tmp_65 = tmp_44 + tmp_49;
      real_t tmp_66 = tmp_32 + tmp_34;
      real_t tmp_67 = 0.5 * p_affine_10_0 * tmp_25 + 0.5 * p_affine_10_1 * tmp_28;
      real_t tmp_68 = tmp_39 * tmp_66;
      real_t tmp_69 = tmp_53 + tmp_55;
      real_t tmp_70 = tmp_39 * tmp_69;
      real_t tmp_71 = tmp_31 + tmp_35;
      real_t tmp_72 = 0.5 * p_affine_10_0 * tmp_26 + 0.5 * p_affine_10_1 * tmp_27;
      real_t tmp_73 = tmp_39 * tmp_71;
      real_t tmp_74 = tmp_52 + tmp_56;
      real_t tmp_75 = tmp_39 * tmp_74;
      real_t a_0_0  = tmp_41 * ( tmp_19 * tmp_29 - tmp_19 * tmp_40 - tmp_36 * tmp_37 ) +
                     tmp_59 * ( tmp_29 * tmp_50 - tmp_37 * tmp_57 - tmp_50 * tmp_58 );
      real_t a_0_1 = tmp_41 * ( tmp_29 * tmp_60 - tmp_36 * tmp_61 - tmp_40 * tmp_60 ) +
                     tmp_59 * ( tmp_29 * tmp_62 - tmp_57 * tmp_61 - tmp_58 * tmp_62 );
      real_t a_0_2 = tmp_41 * ( tmp_29 * tmp_63 - tmp_36 * tmp_64 - tmp_40 * tmp_63 ) +
                     tmp_59 * ( tmp_29 * tmp_65 - tmp_57 * tmp_64 - tmp_58 * tmp_65 );
      real_t a_1_0 = tmp_41 * ( tmp_19 * tmp_67 - tmp_19 * tmp_68 - tmp_37 * tmp_66 ) +
                     tmp_59 * ( -tmp_37 * tmp_69 + tmp_50 * tmp_67 - tmp_50 * tmp_70 );
      real_t a_1_1 = tmp_41 * ( tmp_60 * tmp_67 - tmp_60 * tmp_68 - tmp_61 * tmp_66 ) +
                     tmp_59 * ( -tmp_61 * tmp_69 + tmp_62 * tmp_67 - tmp_62 * tmp_70 );
      real_t a_1_2 = tmp_41 * ( tmp_63 * tmp_67 - tmp_63 * tmp_68 - tmp_64 * tmp_66 ) +
                     tmp_59 * ( -tmp_64 * tmp_69 + tmp_65 * tmp_67 - tmp_65 * tmp_70 );
      real_t a_2_0 = tmp_41 * ( tmp_19 * tmp_72 - tmp_19 * tmp_73 - tmp_37 * tmp_71 ) +
                     tmp_59 * ( -tmp_37 * tmp_74 + tmp_50 * tmp_72 - tmp_50 * tmp_75 );
      real_t a_2_1 = tmp_41 * ( tmp_60 * tmp_72 - tmp_60 * tmp_73 - tmp_61 * tmp_71 ) +
                     tmp_59 * ( -tmp_61 * tmp_74 + tmp_62 * tmp_72 - tmp_62 * tmp_75 );
      real_t a_2_2 = tmp_41 * ( tmp_63 * tmp_72 - tmp_63 * tmp_73 - tmp_64 * tmp_71 ) +
                     tmp_59 * ( -tmp_64 * tmp_74 + tmp_65 * tmp_72 - tmp_65 * tmp_75 );
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

      real_t tmp_0  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_1  = -p_affine_0_1;
      real_t tmp_2  = p_affine_6_1 + tmp_1;
      real_t tmp_3  = 0.21132486540518713 * tmp_0 + tmp_2;
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
      real_t tmp_14 = 0.21132486540518713 * tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6 * tmp_7;
      real_t tmp_16 = tmp_14 * tmp_15;
      real_t tmp_17 = tmp_7 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_18 = tmp_14 * tmp_17;
      real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20 = std::abs( std::pow( ( tmp_0 * tmp_0 ) + ( tmp_12 * tmp_12 ), 1.0 / 2.0 ) );
      real_t tmp_21 = 24 / tmp_20;
      real_t tmp_22 = p_affine_10_0 * ( -tmp_15 - tmp_17 ) + p_affine_10_1 * ( -tmp_10 - tmp_8 );
      real_t tmp_23 = 2 * tmp_22;
      real_t tmp_24 = 0.5 * tmp_20;
      real_t tmp_25 = 0.78867513459481287 * tmp_0 + tmp_2;
      real_t tmp_26 = tmp_25 * tmp_8;
      real_t tmp_27 = tmp_10 * tmp_25;
      real_t tmp_28 = 0.78867513459481287 * tmp_12 + tmp_13;
      real_t tmp_29 = tmp_15 * tmp_28;
      real_t tmp_30 = tmp_17 * tmp_28;
      real_t tmp_31 = -tmp_26 - tmp_27 - tmp_29 - tmp_30 + 1;
      real_t tmp_32 = 0.5 * tmp_20;
      real_t tmp_33 = tmp_11 + tmp_16;
      real_t tmp_34 = p_affine_10_0 * tmp_15 + p_affine_10_1 * tmp_10;
      real_t tmp_35 = tmp_19 * tmp_21;
      real_t tmp_36 = tmp_27 + tmp_29;
      real_t tmp_37 = tmp_21 * tmp_31;
      real_t tmp_38 = tmp_24 * ( -tmp_19 * tmp_34 - tmp_22 * tmp_33 + tmp_33 * tmp_35 ) +
                      tmp_32 * ( -tmp_22 * tmp_36 - tmp_31 * tmp_34 + tmp_36 * tmp_37 );
      real_t tmp_39 = tmp_18 + tmp_9;
      real_t tmp_40 = p_affine_10_0 * tmp_17 + p_affine_10_1 * tmp_8;
      real_t tmp_41 = tmp_26 + tmp_30;
      real_t tmp_42 = tmp_24 * ( -tmp_19 * tmp_40 - tmp_22 * tmp_39 + tmp_35 * tmp_39 ) +
                      tmp_32 * ( -tmp_22 * tmp_41 - tmp_31 * tmp_40 + tmp_37 * tmp_41 );
      real_t tmp_43 = 2 * tmp_34;
      real_t tmp_44 = tmp_24 * ( tmp_21 * tmp_33 * tmp_39 - tmp_33 * tmp_40 - tmp_34 * tmp_39 ) +
                      tmp_32 * ( tmp_21 * tmp_36 * tmp_41 - tmp_34 * tmp_41 - tmp_36 * tmp_40 );
      real_t tmp_45 = 2 * tmp_40;
      real_t a_0_0  = tmp_24 * ( ( tmp_19 * tmp_19 ) * tmp_21 - tmp_19 * tmp_23 ) +
                     tmp_32 * ( tmp_21 * ( tmp_31 * tmp_31 ) - tmp_23 * tmp_31 );
      real_t a_0_1 = tmp_38;
      real_t a_0_2 = tmp_42;
      real_t a_1_0 = tmp_38;
      real_t a_1_1 = tmp_24 * ( tmp_21 * ( tmp_33 * tmp_33 ) - tmp_33 * tmp_43 ) +
                     tmp_32 * ( tmp_21 * ( tmp_36 * tmp_36 ) - tmp_36 * tmp_43 );
      real_t a_1_2 = tmp_44;
      real_t a_2_0 = tmp_42;
      real_t a_2_1 = tmp_44;
      real_t a_2_2 = tmp_24 * ( tmp_21 * ( tmp_39 * tmp_39 ) - tmp_39 * tmp_45 ) +
                     tmp_32 * ( tmp_21 * ( tmp_41 * tmp_41 ) - tmp_41 * tmp_45 );
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
      real_t tmp_3  = 0.21132486540518713 * tmp_0 + tmp_2;
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
      real_t tmp_14 = 0.21132486540518713 * tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6 * tmp_7;
      real_t tmp_16 = tmp_14 * tmp_15;
      real_t tmp_17 = tmp_7 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_18 = tmp_14 * tmp_17;
      real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20 = std::abs( std::pow( ( tmp_0 * tmp_0 ) + ( tmp_12 * tmp_12 ), 1.0 / 2.0 ) );
      real_t tmp_21 = 6 / tmp_20;
      real_t tmp_22 = p_affine_10_0 * ( -tmp_15 - tmp_17 ) + p_affine_10_1 * ( -tmp_10 - tmp_8 );
      real_t tmp_23 = 1.0 * tmp_22;
      real_t tmp_24 = 0.5 * tmp_20;
      real_t tmp_25 = 0.78867513459481287 * tmp_0 + tmp_2;
      real_t tmp_26 = tmp_25 * tmp_8;
      real_t tmp_27 = tmp_10 * tmp_25;
      real_t tmp_28 = 0.78867513459481287 * tmp_12 + tmp_13;
      real_t tmp_29 = tmp_15 * tmp_28;
      real_t tmp_30 = tmp_17 * tmp_28;
      real_t tmp_31 = -tmp_26 - tmp_27 - tmp_29 - tmp_30 + 1;
      real_t tmp_32 = 0.5 * tmp_20;
      real_t tmp_33 = tmp_11 + tmp_16;
      real_t tmp_34 = 0.5 * tmp_22;
      real_t tmp_35 = p_affine_10_0 * tmp_15 + p_affine_10_1 * tmp_10;
      real_t tmp_36 = 0.5 * tmp_35;
      real_t tmp_37 = tmp_19 * tmp_21;
      real_t tmp_38 = tmp_27 + tmp_29;
      real_t tmp_39 = tmp_21 * tmp_31;
      real_t tmp_40 = tmp_24 * ( -tmp_19 * tmp_36 - tmp_33 * tmp_34 + tmp_33 * tmp_37 ) +
                      tmp_32 * ( -tmp_31 * tmp_36 - tmp_34 * tmp_38 + tmp_38 * tmp_39 );
      real_t tmp_41 = tmp_18 + tmp_9;
      real_t tmp_42 = p_affine_10_0 * tmp_17 + p_affine_10_1 * tmp_8;
      real_t tmp_43 = 0.5 * tmp_42;
      real_t tmp_44 = tmp_26 + tmp_30;
      real_t tmp_45 = tmp_24 * ( -tmp_19 * tmp_43 - tmp_34 * tmp_41 + tmp_37 * tmp_41 ) +
                      tmp_32 * ( -tmp_31 * tmp_43 - tmp_34 * tmp_44 + tmp_39 * tmp_44 );
      real_t tmp_46 = 1.0 * tmp_35;
      real_t tmp_47 = tmp_24 * ( tmp_21 * tmp_33 * tmp_41 - tmp_33 * tmp_43 - tmp_36 * tmp_41 ) +
                      tmp_32 * ( tmp_21 * tmp_38 * tmp_44 - tmp_36 * tmp_44 - tmp_38 * tmp_43 );
      real_t tmp_48 = 1.0 * tmp_42;
      real_t a_0_0  = tmp_24 * ( ( tmp_19 * tmp_19 ) * tmp_21 - tmp_19 * tmp_23 ) +
                     tmp_32 * ( tmp_21 * ( tmp_31 * tmp_31 ) - tmp_23 * tmp_31 );
      real_t a_0_1 = tmp_40;
      real_t a_0_2 = tmp_45;
      real_t a_1_0 = tmp_40;
      real_t a_1_1 = tmp_24 * ( tmp_21 * ( tmp_33 * tmp_33 ) - tmp_33 * tmp_46 ) +
                     tmp_32 * ( tmp_21 * ( tmp_38 * tmp_38 ) - tmp_38 * tmp_46 );
      real_t a_1_2 = tmp_47;
      real_t a_2_0 = tmp_45;
      real_t a_2_1 = tmp_47;
      real_t a_2_2 = tmp_24 * ( tmp_21 * ( tmp_41 * tmp_41 ) - tmp_41 * tmp_48 ) +
                     tmp_32 * ( tmp_21 * ( tmp_44 * tmp_44 ) - tmp_44 * tmp_48 );
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

      real_t tmp_0  = -p_affine_3_1;
      real_t tmp_1  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2  = p_affine_6_1 + 0.21132486540518713 * tmp_1;
      real_t tmp_3  = tmp_0 + tmp_2;
      real_t tmp_4  = -p_affine_3_0;
      real_t tmp_5  = p_affine_4_0 + tmp_4;
      real_t tmp_6  = p_affine_5_1 + tmp_0;
      real_t tmp_7  = 1.0 / ( tmp_5 * tmp_6 - ( p_affine_4_1 + tmp_0 ) * ( p_affine_5_0 + tmp_4 ) );
      real_t tmp_8  = tmp_5 * tmp_7;
      real_t tmp_9  = tmp_3 * tmp_8;
      real_t tmp_10 = tmp_7 * ( p_affine_3_0 - p_affine_5_0 );
      real_t tmp_11 = tmp_10 * tmp_3;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + 0.21132486540518713 * tmp_12;
      real_t tmp_14 = tmp_13 + tmp_4;
      real_t tmp_15 = tmp_6 * tmp_7;
      real_t tmp_16 = tmp_14 * tmp_15;
      real_t tmp_17 = tmp_7 * ( p_affine_3_1 - p_affine_4_1 );
      real_t tmp_18 = tmp_14 * tmp_17;
      real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20 = -p_affine_0_1;
      real_t tmp_21 = p_affine_2_1 + tmp_20;
      real_t tmp_22 = -p_affine_0_0;
      real_t tmp_23 = p_affine_1_0 + tmp_22;
      real_t tmp_24 = 1.0 / ( tmp_21 * tmp_23 - ( p_affine_1_1 + tmp_20 ) * ( p_affine_2_0 + tmp_22 ) );
      real_t tmp_25 = tmp_21 * tmp_24;
      real_t tmp_26 = tmp_24 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_27 = tmp_23 * tmp_24;
      real_t tmp_28 = tmp_24 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_29 = 0.5 * p_affine_10_0 * ( -tmp_25 - tmp_26 ) + 0.5 * p_affine_10_1 * ( -tmp_27 - tmp_28 );
      real_t tmp_30 = tmp_2 + tmp_20;
      real_t tmp_31 = tmp_27 * tmp_30;
      real_t tmp_32 = tmp_28 * tmp_30;
      real_t tmp_33 = tmp_13 + tmp_22;
      real_t tmp_34 = tmp_25 * tmp_33;
      real_t tmp_35 = tmp_26 * tmp_33;
      real_t tmp_36 = -tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1;
      real_t tmp_37 = 0.5 * p_affine_10_0 * ( -tmp_15 - tmp_17 ) + 0.5 * p_affine_10_1 * ( -tmp_10 - tmp_8 );
      real_t tmp_38 = std::abs( std::pow( ( tmp_1 * tmp_1 ) + ( tmp_12 * tmp_12 ), 1.0 / 2.0 ) );
      real_t tmp_39 = 6 / tmp_38;
      real_t tmp_40 = tmp_36 * tmp_39;
      real_t tmp_41 = 0.5 * tmp_38;
      real_t tmp_42 = p_affine_6_1 + 0.78867513459481287 * tmp_1;
      real_t tmp_43 = tmp_0 + tmp_42;
      real_t tmp_44 = tmp_43 * tmp_8;
      real_t tmp_45 = tmp_10 * tmp_43;
      real_t tmp_46 = p_affine_6_0 + 0.78867513459481287 * tmp_12;
      real_t tmp_47 = tmp_4 + tmp_46;
      real_t tmp_48 = tmp_15 * tmp_47;
      real_t tmp_49 = tmp_17 * tmp_47;
      real_t tmp_50 = -tmp_44 - tmp_45 - tmp_48 - tmp_49 + 1;
      real_t tmp_51 = tmp_20 + tmp_42;
      real_t tmp_52 = tmp_27 * tmp_51;
      real_t tmp_53 = tmp_28 * tmp_51;
      real_t tmp_54 = tmp_22 + tmp_46;
      real_t tmp_55 = tmp_25 * tmp_54;
      real_t tmp_56 = tmp_26 * tmp_54;
      real_t tmp_57 = -tmp_52 - tmp_53 - tmp_55 - tmp_56 + 1;
      real_t tmp_58 = tmp_39 * tmp_57;
      real_t tmp_59 = 0.5 * tmp_38;
      real_t tmp_60 = tmp_11 + tmp_16;
      real_t tmp_61 = 0.5 * p_affine_10_0 * tmp_15 + 0.5 * p_affine_10_1 * tmp_10;
      real_t tmp_62 = tmp_45 + tmp_48;
      real_t tmp_63 = tmp_18 + tmp_9;
      real_t tmp_64 = 0.5 * p_affine_10_0 * tmp_17 + 0.5 * p_affine_10_1 * tmp_8;
      real_t tmp_65 = tmp_44 + tmp_49;
      real_t tmp_66 = tmp_32 + tmp_34;
      real_t tmp_67 = 0.5 * p_affine_10_0 * tmp_25 + 0.5 * p_affine_10_1 * tmp_28;
      real_t tmp_68 = tmp_39 * tmp_66;
      real_t tmp_69 = tmp_53 + tmp_55;
      real_t tmp_70 = tmp_39 * tmp_69;
      real_t tmp_71 = tmp_31 + tmp_35;
      real_t tmp_72 = 0.5 * p_affine_10_0 * tmp_26 + 0.5 * p_affine_10_1 * tmp_27;
      real_t tmp_73 = tmp_39 * tmp_71;
      real_t tmp_74 = tmp_52 + tmp_56;
      real_t tmp_75 = tmp_39 * tmp_74;
      real_t a_0_0  = tmp_41 * ( tmp_19 * tmp_29 - tmp_19 * tmp_40 - tmp_36 * tmp_37 ) +
                     tmp_59 * ( tmp_29 * tmp_50 - tmp_37 * tmp_57 - tmp_50 * tmp_58 );
      real_t a_0_1 = tmp_41 * ( tmp_29 * tmp_60 - tmp_36 * tmp_61 - tmp_40 * tmp_60 ) +
                     tmp_59 * ( tmp_29 * tmp_62 - tmp_57 * tmp_61 - tmp_58 * tmp_62 );
      real_t a_0_2 = tmp_41 * ( tmp_29 * tmp_63 - tmp_36 * tmp_64 - tmp_40 * tmp_63 ) +
                     tmp_59 * ( tmp_29 * tmp_65 - tmp_57 * tmp_64 - tmp_58 * tmp_65 );
      real_t a_1_0 = tmp_41 * ( tmp_19 * tmp_67 - tmp_19 * tmp_68 - tmp_37 * tmp_66 ) +
                     tmp_59 * ( -tmp_37 * tmp_69 + tmp_50 * tmp_67 - tmp_50 * tmp_70 );
      real_t a_1_1 = tmp_41 * ( tmp_60 * tmp_67 - tmp_60 * tmp_68 - tmp_61 * tmp_66 ) +
                     tmp_59 * ( -tmp_61 * tmp_69 + tmp_62 * tmp_67 - tmp_62 * tmp_70 );
      real_t a_1_2 = tmp_41 * ( tmp_63 * tmp_67 - tmp_63 * tmp_68 - tmp_64 * tmp_66 ) +
                     tmp_59 * ( -tmp_64 * tmp_69 + tmp_65 * tmp_67 - tmp_65 * tmp_70 );
      real_t a_2_0 = tmp_41 * ( tmp_19 * tmp_72 - tmp_19 * tmp_73 - tmp_37 * tmp_71 ) +
                     tmp_59 * ( -tmp_37 * tmp_74 + tmp_50 * tmp_72 - tmp_50 * tmp_75 );
      real_t a_2_1 = tmp_41 * ( tmp_60 * tmp_72 - tmp_60 * tmp_73 - tmp_61 * tmp_71 ) +
                     tmp_59 * ( -tmp_61 * tmp_74 + tmp_62 * tmp_72 - tmp_62 * tmp_75 );
      real_t a_2_2 = tmp_41 * ( tmp_63 * tmp_72 - tmp_63 * tmp_73 - tmp_64 * tmp_71 ) +
                     tmp_59 * ( -tmp_64 * tmp_74 + tmp_65 * tmp_72 - tmp_65 * tmp_75 );
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

      real_t tmp_0  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_1  = -p_affine_0_1;
      real_t tmp_2  = p_affine_6_1 + tmp_1;
      real_t tmp_3  = 0.21132486540518713 * tmp_0 + tmp_2;
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
      real_t tmp_14 = 0.21132486540518713 * tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6 * tmp_7;
      real_t tmp_16 = tmp_14 * tmp_15;
      real_t tmp_17 = tmp_7 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_18 = tmp_14 * tmp_17;
      real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20 = std::abs( std::pow( ( tmp_0 * tmp_0 ) + ( tmp_12 * tmp_12 ), 1.0 / 2.0 ) );
      real_t tmp_21 = 24 / tmp_20;
      real_t tmp_22 = p_affine_10_0 * ( -tmp_15 - tmp_17 ) + p_affine_10_1 * ( -tmp_10 - tmp_8 );
      real_t tmp_23 = 2 * tmp_22;
      real_t tmp_24 = 0.5 * tmp_20;
      real_t tmp_25 = 0.78867513459481287 * tmp_0 + tmp_2;
      real_t tmp_26 = tmp_25 * tmp_8;
      real_t tmp_27 = tmp_10 * tmp_25;
      real_t tmp_28 = 0.78867513459481287 * tmp_12 + tmp_13;
      real_t tmp_29 = tmp_15 * tmp_28;
      real_t tmp_30 = tmp_17 * tmp_28;
      real_t tmp_31 = -tmp_26 - tmp_27 - tmp_29 - tmp_30 + 1;
      real_t tmp_32 = 0.5 * tmp_20;
      real_t tmp_33 = tmp_11 + tmp_16;
      real_t tmp_34 = p_affine_10_0 * tmp_15 + p_affine_10_1 * tmp_10;
      real_t tmp_35 = tmp_19 * tmp_21;
      real_t tmp_36 = tmp_27 + tmp_29;
      real_t tmp_37 = tmp_21 * tmp_31;
      real_t tmp_38 = tmp_24 * ( -tmp_19 * tmp_34 - tmp_22 * tmp_33 + tmp_33 * tmp_35 ) +
                      tmp_32 * ( -tmp_22 * tmp_36 - tmp_31 * tmp_34 + tmp_36 * tmp_37 );
      real_t tmp_39 = tmp_18 + tmp_9;
      real_t tmp_40 = p_affine_10_0 * tmp_17 + p_affine_10_1 * tmp_8;
      real_t tmp_41 = tmp_26 + tmp_30;
      real_t tmp_42 = tmp_24 * ( -tmp_19 * tmp_40 - tmp_22 * tmp_39 + tmp_35 * tmp_39 ) +
                      tmp_32 * ( -tmp_22 * tmp_41 - tmp_31 * tmp_40 + tmp_37 * tmp_41 );
      real_t tmp_43 = 2 * tmp_34;
      real_t tmp_44 = tmp_24 * ( tmp_21 * tmp_33 * tmp_39 - tmp_33 * tmp_40 - tmp_34 * tmp_39 ) +
                      tmp_32 * ( tmp_21 * tmp_36 * tmp_41 - tmp_34 * tmp_41 - tmp_36 * tmp_40 );
      real_t tmp_45 = 2 * tmp_40;
      real_t a_0_0  = tmp_24 * ( ( tmp_19 * tmp_19 ) * tmp_21 - tmp_19 * tmp_23 ) +
                     tmp_32 * ( tmp_21 * ( tmp_31 * tmp_31 ) - tmp_23 * tmp_31 );
      real_t a_0_1 = tmp_38;
      real_t a_0_2 = tmp_42;
      real_t a_1_0 = tmp_38;
      real_t a_1_1 = tmp_24 * ( tmp_21 * ( tmp_33 * tmp_33 ) - tmp_33 * tmp_43 ) +
                     tmp_32 * ( tmp_21 * ( tmp_36 * tmp_36 ) - tmp_36 * tmp_43 );
      real_t a_1_2 = tmp_44;
      real_t a_2_0 = tmp_42;
      real_t a_2_1 = tmp_44;
      real_t a_2_2 = tmp_24 * ( tmp_21 * ( tmp_39 * tmp_39 ) - tmp_39 * tmp_45 ) +
                     tmp_32 * ( tmp_21 * ( tmp_41 * tmp_41 ) - tmp_41 * tmp_45 );
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

class DGVectorLaplaceFormP1EDG_0 : public hyteg::dg::DGForm2D
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
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_0;
      real_t tmp_6  = 1.0 / ( tmp_4 - tmp_5 * ( p_affine_1_1 + tmp_2 ) );
      real_t tmp_7  = tmp_1 * tmp_6;
      real_t tmp_8  = tmp_6 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_9  = tmp_1 * tmp_8 + tmp_5 * tmp_7;
      real_t tmp_10 = tmp_3 * tmp_6;
      real_t tmp_11 = tmp_6 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_12 = tmp_11 * tmp_5 + tmp_4 * tmp_6;
      real_t tmp_13 = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                                p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_14 = tmp_13 * ( tmp_12 * ( -tmp_10 - tmp_11 ) + tmp_9 * ( -tmp_7 - tmp_8 ) );
      real_t tmp_15 = tmp_13 * ( tmp_10 * tmp_12 + tmp_8 * tmp_9 );
      real_t tmp_16 = tmp_13 * ( tmp_11 * tmp_12 + tmp_7 * tmp_9 );
      real_t a_0_0  = 0.5 * tmp_14;
      real_t a_1_0  = 0.5 * tmp_15;
      real_t a_2_0  = 0.5 * tmp_16;
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_3  = -p_affine_0_1;
      real_t tmp_4  = p_affine_6_1 + tmp_3;
      real_t tmp_5  = 0.21132486540518713 * tmp_2 + tmp_4;
      real_t tmp_6  = p_affine_2_1 + tmp_3;
      real_t tmp_7  = tmp_1 * tmp_6;
      real_t tmp_8  = p_affine_2_0 + tmp_0;
      real_t tmp_9  = 1.0 / ( tmp_7 - tmp_8 * ( p_affine_1_1 + tmp_3 ) );
      real_t tmp_10 = tmp_9 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_11 = tmp_10 * tmp_5;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_0;
      real_t tmp_14 = 0.21132486540518713 * tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6 * tmp_9;
      real_t tmp_16 = tmp_14 * tmp_15;
      real_t tmp_17 = tmp_11 + tmp_16;
      real_t tmp_18 = tmp_1 * tmp_9;
      real_t tmp_19 = tmp_18 * tmp_5;
      real_t tmp_20 = tmp_9 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_21 = tmp_14 * tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_1 * ( tmp_17 - 1.0 / 3.0 ) + tmp_8 * ( tmp_22 - 1.0 / 3.0 );
      real_t tmp_24 = 0.5 * p_affine_10_0 * ( -tmp_15 - tmp_20 ) + 0.5 * p_affine_10_1 * ( -tmp_10 - tmp_18 );
      real_t tmp_25 = -tmp_11 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_26 =
          0.5 * p_affine_10_0 * ( tmp_20 * tmp_8 + tmp_7 * tmp_9 ) + 0.5 * p_affine_10_1 * ( tmp_1 * tmp_10 + tmp_18 * tmp_8 );
      real_t tmp_27 = std::abs( std::pow( ( tmp_12 * tmp_12 ) + ( tmp_2 * tmp_2 ), 1.0 / 2.0 ) );
      real_t tmp_28 = 6 / tmp_27;
      real_t tmp_29 = tmp_23 * tmp_28;
      real_t tmp_30 = 0.5 * tmp_27;
      real_t tmp_31 = 0.78867513459481287 * tmp_2 + tmp_4;
      real_t tmp_32 = tmp_10 * tmp_31;
      real_t tmp_33 = 0.78867513459481287 * tmp_12 + tmp_13;
      real_t tmp_34 = tmp_15 * tmp_33;
      real_t tmp_35 = tmp_32 + tmp_34;
      real_t tmp_36 = tmp_18 * tmp_31;
      real_t tmp_37 = tmp_20 * tmp_33;
      real_t tmp_38 = tmp_36 + tmp_37;
      real_t tmp_39 = tmp_1 * ( tmp_35 - 1.0 / 3.0 ) + tmp_8 * ( tmp_38 - 1.0 / 3.0 );
      real_t tmp_40 = -tmp_32 - tmp_34 - tmp_36 - tmp_37 + 1;
      real_t tmp_41 = tmp_28 * tmp_39;
      real_t tmp_42 = 0.5 * tmp_27;
      real_t tmp_43 = 0.5 * p_affine_10_0 * tmp_15 + 0.5 * p_affine_10_1 * tmp_10;
      real_t tmp_44 = 0.5 * p_affine_10_0 * tmp_20 + 0.5 * p_affine_10_1 * tmp_18;
      real_t a_0_0  = tmp_30 * ( -tmp_23 * tmp_24 - tmp_25 * tmp_26 + tmp_25 * tmp_29 ) +
                     tmp_42 * ( -tmp_24 * tmp_39 - tmp_26 * tmp_40 + tmp_40 * tmp_41 );
      real_t a_1_0 = tmp_30 * ( -tmp_17 * tmp_26 + tmp_17 * tmp_29 - tmp_23 * tmp_43 ) +
                     tmp_42 * ( -tmp_26 * tmp_35 + tmp_35 * tmp_41 - tmp_39 * tmp_43 );
      real_t a_2_0 = tmp_30 * ( -tmp_22 * tmp_26 + tmp_22 * tmp_29 - tmp_23 * tmp_44 ) +
                     tmp_42 * ( -tmp_26 * tmp_38 + tmp_38 * tmp_41 - tmp_39 * tmp_44 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
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

      real_t tmp_0  = -p_affine_3_0;
      real_t tmp_1  = p_affine_4_0 + tmp_0;
      real_t tmp_2  = p_affine_3_0 - p_affine_5_0;
      real_t tmp_3  = -p_affine_3_1;
      real_t tmp_4  = p_affine_5_1 + tmp_3;
      real_t tmp_5  = tmp_1 * tmp_4;
      real_t tmp_6  = p_affine_5_0 + tmp_0;
      real_t tmp_7  = 1.0 / ( tmp_5 - tmp_6 * ( p_affine_4_1 + tmp_3 ) );
      real_t tmp_8  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9  = p_affine_6_1 + 0.21132486540518713 * tmp_8;
      real_t tmp_10 = tmp_7 * ( tmp_3 + tmp_9 );
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.21132486540518713 * tmp_11;
      real_t tmp_13 = tmp_7 * ( tmp_0 + tmp_12 );
      real_t tmp_14 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_15 =
          tmp_1 * ( tmp_10 * tmp_2 + tmp_13 * tmp_4 - 1.0 / 3.0 ) + tmp_6 * ( tmp_1 * tmp_10 + tmp_13 * tmp_14 - 1.0 / 3.0 );
      real_t tmp_16 = -p_affine_0_1;
      real_t tmp_17 = p_affine_2_1 + tmp_16;
      real_t tmp_18 = -p_affine_0_0;
      real_t tmp_19 = p_affine_1_0 + tmp_18;
      real_t tmp_20 = 1.0 / ( tmp_17 * tmp_19 - ( p_affine_1_1 + tmp_16 ) * ( p_affine_2_0 + tmp_18 ) );
      real_t tmp_21 = tmp_17 * tmp_20;
      real_t tmp_22 = tmp_20 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_23 = tmp_19 * tmp_20;
      real_t tmp_24 = tmp_20 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_25 = 0.5 * p_affine_10_0 * ( -tmp_21 - tmp_22 ) + 0.5 * p_affine_10_1 * ( -tmp_23 - tmp_24 );
      real_t tmp_26 = tmp_16 + tmp_9;
      real_t tmp_27 = tmp_23 * tmp_26;
      real_t tmp_28 = tmp_24 * tmp_26;
      real_t tmp_29 = tmp_12 + tmp_18;
      real_t tmp_30 = tmp_21 * tmp_29;
      real_t tmp_31 = tmp_22 * tmp_29;
      real_t tmp_32 = -tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1;
      real_t tmp_33 = tmp_6 * tmp_7;
      real_t tmp_34 = tmp_2 * tmp_7;
      real_t tmp_35 =
          0.5 * p_affine_10_0 * ( tmp_14 * tmp_33 + tmp_5 * tmp_7 ) + 0.5 * p_affine_10_1 * ( tmp_1 * tmp_33 + tmp_1 * tmp_34 );
      real_t tmp_36 = std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_8 * tmp_8 ), 1.0 / 2.0 ) );
      real_t tmp_37 = 6 / tmp_36;
      real_t tmp_38 = tmp_15 * tmp_37;
      real_t tmp_39 = 0.5 * tmp_36;
      real_t tmp_40 = p_affine_6_1 + 0.78867513459481287 * tmp_8;
      real_t tmp_41 = tmp_3 + tmp_40;
      real_t tmp_42 = p_affine_6_0 + 0.78867513459481287 * tmp_11;
      real_t tmp_43 = tmp_7 * ( tmp_0 + tmp_42 );
      real_t tmp_44 = tmp_1 * ( tmp_34 * tmp_41 + tmp_4 * tmp_43 - 1.0 / 3.0 ) +
                      tmp_6 * ( tmp_1 * tmp_41 * tmp_7 + tmp_14 * tmp_43 - 1.0 / 3.0 );
      real_t tmp_45 = tmp_16 + tmp_40;
      real_t tmp_46 = tmp_23 * tmp_45;
      real_t tmp_47 = tmp_24 * tmp_45;
      real_t tmp_48 = tmp_18 + tmp_42;
      real_t tmp_49 = tmp_21 * tmp_48;
      real_t tmp_50 = tmp_22 * tmp_48;
      real_t tmp_51 = -tmp_46 - tmp_47 - tmp_49 - tmp_50 + 1;
      real_t tmp_52 = tmp_37 * tmp_44;
      real_t tmp_53 = 0.5 * tmp_36;
      real_t tmp_54 = tmp_28 + tmp_30;
      real_t tmp_55 = 0.5 * p_affine_10_0 * tmp_21 + 0.5 * p_affine_10_1 * tmp_24;
      real_t tmp_56 = tmp_47 + tmp_49;
      real_t tmp_57 = tmp_27 + tmp_31;
      real_t tmp_58 = 0.5 * p_affine_10_0 * tmp_22 + 0.5 * p_affine_10_1 * tmp_23;
      real_t tmp_59 = tmp_46 + tmp_50;
      real_t a_0_0  = tmp_39 * ( tmp_15 * tmp_25 - tmp_32 * tmp_35 - tmp_32 * tmp_38 ) +
                     tmp_53 * ( tmp_25 * tmp_44 - tmp_35 * tmp_51 - tmp_51 * tmp_52 );
      real_t a_1_0 = tmp_39 * ( tmp_15 * tmp_55 - tmp_35 * tmp_54 - tmp_38 * tmp_54 ) +
                     tmp_53 * ( -tmp_35 * tmp_56 + tmp_44 * tmp_55 - tmp_52 * tmp_56 );
      real_t a_2_0 = tmp_39 * ( tmp_15 * tmp_58 - tmp_35 * tmp_57 - tmp_38 * tmp_57 ) +
                     tmp_53 * ( -tmp_35 * tmp_59 + tmp_44 * tmp_58 - tmp_52 * tmp_59 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
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
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_2;
      real_t tmp_6  = 1.0 / ( tmp_4 - tmp_5 * ( p_affine_1_1 + tmp_0 ) );
      real_t tmp_7  = tmp_1 * tmp_6;
      real_t tmp_8  = tmp_6 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_9  = tmp_3 * tmp_6;
      real_t tmp_10 = tmp_6 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_11 = p_affine_10_0 * ( -tmp_7 - tmp_8 ) + p_affine_10_1 * ( -tmp_10 - tmp_9 );
      real_t tmp_12 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_13 = p_affine_6_1 + tmp_0;
      real_t tmp_14 = 0.21132486540518713 * tmp_12 + tmp_13;
      real_t tmp_15 = tmp_10 * tmp_14;
      real_t tmp_16 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_17 = p_affine_6_0 + tmp_2;
      real_t tmp_18 = 0.21132486540518713 * tmp_16 + tmp_17;
      real_t tmp_19 = tmp_18 * tmp_7;
      real_t tmp_20 = tmp_15 + tmp_19;
      real_t tmp_21 = tmp_14 * tmp_9;
      real_t tmp_22 = tmp_18 * tmp_8;
      real_t tmp_23 = tmp_21 + tmp_22;
      real_t tmp_24 = tmp_3 * ( tmp_20 - 1.0 / 3.0 ) + tmp_5 * ( tmp_23 - 1.0 / 3.0 );
      real_t tmp_25 = p_affine_10_0 * ( tmp_4 * tmp_6 + tmp_5 * tmp_8 ) + p_affine_10_1 * ( tmp_10 * tmp_3 + tmp_5 * tmp_9 );
      real_t tmp_26 = -tmp_15 - tmp_19 - tmp_21 - tmp_22 + 1;
      real_t tmp_27 = std::abs( std::pow( ( tmp_12 * tmp_12 ) + ( tmp_16 * tmp_16 ), 1.0 / 2.0 ) );
      real_t tmp_28 = 24 / tmp_27;
      real_t tmp_29 = tmp_24 * tmp_28;
      real_t tmp_30 = 0.5 * tmp_27;
      real_t tmp_31 = 0.78867513459481287 * tmp_12 + tmp_13;
      real_t tmp_32 = tmp_10 * tmp_31;
      real_t tmp_33 = 0.78867513459481287 * tmp_16 + tmp_17;
      real_t tmp_34 = tmp_33 * tmp_7;
      real_t tmp_35 = tmp_32 + tmp_34;
      real_t tmp_36 = tmp_31 * tmp_9;
      real_t tmp_37 = tmp_33 * tmp_8;
      real_t tmp_38 = tmp_36 + tmp_37;
      real_t tmp_39 = tmp_3 * ( tmp_35 - 1.0 / 3.0 ) + tmp_5 * ( tmp_38 - 1.0 / 3.0 );
      real_t tmp_40 = -tmp_32 - tmp_34 - tmp_36 - tmp_37 + 1;
      real_t tmp_41 = tmp_28 * tmp_39;
      real_t tmp_42 = 0.5 * tmp_27;
      real_t tmp_43 = p_affine_10_0 * tmp_7 + p_affine_10_1 * tmp_10;
      real_t tmp_44 = p_affine_10_0 * tmp_8 + p_affine_10_1 * tmp_9;
      real_t a_0_0  = tmp_30 * ( -tmp_11 * tmp_24 - tmp_25 * tmp_26 + tmp_26 * tmp_29 ) +
                     tmp_42 * ( -tmp_11 * tmp_39 - tmp_25 * tmp_40 + tmp_40 * tmp_41 );
      real_t a_1_0 = tmp_30 * ( -tmp_20 * tmp_25 + tmp_20 * tmp_29 - tmp_24 * tmp_43 ) +
                     tmp_42 * ( -tmp_25 * tmp_35 + tmp_35 * tmp_41 - tmp_39 * tmp_43 );
      real_t a_2_0 = tmp_30 * ( -tmp_23 * tmp_25 + tmp_23 * tmp_29 - tmp_24 * tmp_44 ) +
                     tmp_42 * ( -tmp_25 * tmp_38 + tmp_38 * tmp_41 - tmp_39 * tmp_44 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
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

class DGVectorLaplaceFormP1EDG_1 : public hyteg::dg::DGForm2D
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
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_1_1 + tmp_2;
      real_t tmp_6  = 1.0 / ( tmp_4 - tmp_5 * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_7  = tmp_1 * tmp_6;
      real_t tmp_8  = tmp_6 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_9  = tmp_4 * tmp_6 + tmp_5 * tmp_8;
      real_t tmp_10 = tmp_3 * tmp_6;
      real_t tmp_11 = tmp_6 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_12 = tmp_10 * tmp_5 + tmp_11 * tmp_3;
      real_t tmp_13 = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                                p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_14 = tmp_13 * ( tmp_12 * ( -tmp_10 - tmp_11 ) + tmp_9 * ( -tmp_7 - tmp_8 ) );
      real_t tmp_15 = tmp_13 * ( tmp_10 * tmp_12 + tmp_8 * tmp_9 );
      real_t tmp_16 = tmp_13 * ( tmp_11 * tmp_12 + tmp_7 * tmp_9 );
      real_t a_0_0  = 0.5 * tmp_14;
      real_t a_1_0  = 0.5 * tmp_15;
      real_t a_2_0  = 0.5 * tmp_16;
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

      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_1_1 + tmp_0;
      real_t tmp_2  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_3  = p_affine_6_1 + tmp_0;
      real_t tmp_4  = 0.21132486540518713 * tmp_2 + tmp_3;
      real_t tmp_5  = -p_affine_0_0;
      real_t tmp_6  = p_affine_1_0 + tmp_5;
      real_t tmp_7  = p_affine_2_1 + tmp_0;
      real_t tmp_8  = tmp_6 * tmp_7;
      real_t tmp_9  = 1.0 / ( -tmp_1 * ( p_affine_2_0 + tmp_5 ) + tmp_8 );
      real_t tmp_10 = tmp_9 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_11 = tmp_10 * tmp_4;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_5;
      real_t tmp_14 = 0.21132486540518713 * tmp_12 + tmp_13;
      real_t tmp_15 = tmp_7 * tmp_9;
      real_t tmp_16 = tmp_14 * tmp_15;
      real_t tmp_17 = tmp_11 + tmp_16;
      real_t tmp_18 = tmp_6 * tmp_9;
      real_t tmp_19 = tmp_18 * tmp_4;
      real_t tmp_20 = tmp_9 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_21 = tmp_14 * tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_1 * ( tmp_17 - 1.0 / 3.0 ) + tmp_7 * ( tmp_22 - 1.0 / 3.0 );
      real_t tmp_24 = 0.5 * p_affine_10_0 * ( -tmp_15 - tmp_20 ) + 0.5 * p_affine_10_1 * ( -tmp_10 - tmp_18 );
      real_t tmp_25 = -tmp_11 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_26 =
          0.5 * p_affine_10_0 * ( tmp_1 * tmp_15 + tmp_20 * tmp_7 ) + 0.5 * p_affine_10_1 * ( tmp_1 * tmp_10 + tmp_8 * tmp_9 );
      real_t tmp_27 = std::abs( std::pow( ( tmp_12 * tmp_12 ) + ( tmp_2 * tmp_2 ), 1.0 / 2.0 ) );
      real_t tmp_28 = 6 / tmp_27;
      real_t tmp_29 = tmp_23 * tmp_28;
      real_t tmp_30 = 0.5 * tmp_27;
      real_t tmp_31 = 0.78867513459481287 * tmp_2 + tmp_3;
      real_t tmp_32 = tmp_10 * tmp_31;
      real_t tmp_33 = 0.78867513459481287 * tmp_12 + tmp_13;
      real_t tmp_34 = tmp_15 * tmp_33;
      real_t tmp_35 = tmp_32 + tmp_34;
      real_t tmp_36 = tmp_18 * tmp_31;
      real_t tmp_37 = tmp_20 * tmp_33;
      real_t tmp_38 = tmp_36 + tmp_37;
      real_t tmp_39 = tmp_1 * ( tmp_35 - 1.0 / 3.0 ) + tmp_7 * ( tmp_38 - 1.0 / 3.0 );
      real_t tmp_40 = -tmp_32 - tmp_34 - tmp_36 - tmp_37 + 1;
      real_t tmp_41 = tmp_28 * tmp_39;
      real_t tmp_42 = 0.5 * tmp_27;
      real_t tmp_43 = 0.5 * p_affine_10_0 * tmp_15 + 0.5 * p_affine_10_1 * tmp_10;
      real_t tmp_44 = 0.5 * p_affine_10_0 * tmp_20 + 0.5 * p_affine_10_1 * tmp_18;
      real_t a_0_0  = tmp_30 * ( -tmp_23 * tmp_24 - tmp_25 * tmp_26 + tmp_25 * tmp_29 ) +
                     tmp_42 * ( -tmp_24 * tmp_39 - tmp_26 * tmp_40 + tmp_40 * tmp_41 );
      real_t a_1_0 = tmp_30 * ( -tmp_17 * tmp_26 + tmp_17 * tmp_29 - tmp_23 * tmp_43 ) +
                     tmp_42 * ( -tmp_26 * tmp_35 + tmp_35 * tmp_41 - tmp_39 * tmp_43 );
      real_t a_2_0 = tmp_30 * ( -tmp_22 * tmp_26 + tmp_22 * tmp_29 - tmp_23 * tmp_44 ) +
                     tmp_42 * ( -tmp_26 * tmp_38 + tmp_38 * tmp_41 - tmp_39 * tmp_44 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
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

      real_t tmp_0  = -p_affine_3_1;
      real_t tmp_1  = p_affine_4_1 + tmp_0;
      real_t tmp_2  = p_affine_3_0 - p_affine_5_0;
      real_t tmp_3  = -p_affine_3_0;
      real_t tmp_4  = p_affine_4_0 + tmp_3;
      real_t tmp_5  = p_affine_5_1 + tmp_0;
      real_t tmp_6  = tmp_4 * tmp_5;
      real_t tmp_7  = 1.0 / ( -tmp_1 * ( p_affine_5_0 + tmp_3 ) + tmp_6 );
      real_t tmp_8  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9  = p_affine_6_1 + 0.21132486540518713 * tmp_8;
      real_t tmp_10 = tmp_7 * ( tmp_0 + tmp_9 );
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.21132486540518713 * tmp_11;
      real_t tmp_13 = tmp_7 * ( tmp_12 + tmp_3 );
      real_t tmp_14 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_15 =
          tmp_1 * ( tmp_10 * tmp_2 + tmp_13 * tmp_5 - 1.0 / 3.0 ) + tmp_5 * ( tmp_10 * tmp_4 + tmp_13 * tmp_14 - 1.0 / 3.0 );
      real_t tmp_16 = -p_affine_0_1;
      real_t tmp_17 = p_affine_2_1 + tmp_16;
      real_t tmp_18 = -p_affine_0_0;
      real_t tmp_19 = p_affine_1_0 + tmp_18;
      real_t tmp_20 = 1.0 / ( tmp_17 * tmp_19 - ( p_affine_1_1 + tmp_16 ) * ( p_affine_2_0 + tmp_18 ) );
      real_t tmp_21 = tmp_17 * tmp_20;
      real_t tmp_22 = tmp_20 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_23 = tmp_19 * tmp_20;
      real_t tmp_24 = tmp_20 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_25 = 0.5 * p_affine_10_0 * ( -tmp_21 - tmp_22 ) + 0.5 * p_affine_10_1 * ( -tmp_23 - tmp_24 );
      real_t tmp_26 = tmp_16 + tmp_9;
      real_t tmp_27 = tmp_23 * tmp_26;
      real_t tmp_28 = tmp_24 * tmp_26;
      real_t tmp_29 = tmp_12 + tmp_18;
      real_t tmp_30 = tmp_21 * tmp_29;
      real_t tmp_31 = tmp_22 * tmp_29;
      real_t tmp_32 = -tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1;
      real_t tmp_33 = tmp_5 * tmp_7;
      real_t tmp_34 = tmp_2 * tmp_7;
      real_t tmp_35 =
          0.5 * p_affine_10_0 * ( tmp_1 * tmp_33 + tmp_14 * tmp_33 ) + 0.5 * p_affine_10_1 * ( tmp_1 * tmp_34 + tmp_6 * tmp_7 );
      real_t tmp_36 = std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_8 * tmp_8 ), 1.0 / 2.0 ) );
      real_t tmp_37 = 6 / tmp_36;
      real_t tmp_38 = tmp_15 * tmp_37;
      real_t tmp_39 = 0.5 * tmp_36;
      real_t tmp_40 = p_affine_6_1 + 0.78867513459481287 * tmp_8;
      real_t tmp_41 = tmp_0 + tmp_40;
      real_t tmp_42 = p_affine_6_0 + 0.78867513459481287 * tmp_11;
      real_t tmp_43 = tmp_3 + tmp_42;
      real_t tmp_44 = tmp_1 * ( tmp_33 * tmp_43 + tmp_34 * tmp_41 - 1.0 / 3.0 ) +
                      tmp_5 * ( tmp_14 * tmp_43 * tmp_7 + tmp_4 * tmp_41 * tmp_7 - 1.0 / 3.0 );
      real_t tmp_45 = tmp_16 + tmp_40;
      real_t tmp_46 = tmp_23 * tmp_45;
      real_t tmp_47 = tmp_24 * tmp_45;
      real_t tmp_48 = tmp_18 + tmp_42;
      real_t tmp_49 = tmp_21 * tmp_48;
      real_t tmp_50 = tmp_22 * tmp_48;
      real_t tmp_51 = -tmp_46 - tmp_47 - tmp_49 - tmp_50 + 1;
      real_t tmp_52 = tmp_37 * tmp_44;
      real_t tmp_53 = 0.5 * tmp_36;
      real_t tmp_54 = tmp_28 + tmp_30;
      real_t tmp_55 = 0.5 * p_affine_10_0 * tmp_21 + 0.5 * p_affine_10_1 * tmp_24;
      real_t tmp_56 = tmp_47 + tmp_49;
      real_t tmp_57 = tmp_27 + tmp_31;
      real_t tmp_58 = 0.5 * p_affine_10_0 * tmp_22 + 0.5 * p_affine_10_1 * tmp_23;
      real_t tmp_59 = tmp_46 + tmp_50;
      real_t a_0_0  = tmp_39 * ( tmp_15 * tmp_25 - tmp_32 * tmp_35 - tmp_32 * tmp_38 ) +
                     tmp_53 * ( tmp_25 * tmp_44 - tmp_35 * tmp_51 - tmp_51 * tmp_52 );
      real_t a_1_0 = tmp_39 * ( tmp_15 * tmp_55 - tmp_35 * tmp_54 - tmp_38 * tmp_54 ) +
                     tmp_53 * ( -tmp_35 * tmp_56 + tmp_44 * tmp_55 - tmp_52 * tmp_56 );
      real_t a_2_0 = tmp_39 * ( tmp_15 * tmp_58 - tmp_35 * tmp_57 - tmp_38 * tmp_57 ) +
                     tmp_53 * ( -tmp_35 * tmp_59 + tmp_44 * tmp_58 - tmp_52 * tmp_59 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
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
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_1_1 + tmp_0;
      real_t tmp_6  = 1.0 / ( tmp_4 - tmp_5 * ( p_affine_2_0 + tmp_2 ) );
      real_t tmp_7  = tmp_1 * tmp_6;
      real_t tmp_8  = tmp_6 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_9  = tmp_3 * tmp_6;
      real_t tmp_10 = tmp_6 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_11 = p_affine_10_0 * ( -tmp_7 - tmp_8 ) + p_affine_10_1 * ( -tmp_10 - tmp_9 );
      real_t tmp_12 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_13 = p_affine_6_1 + tmp_0;
      real_t tmp_14 = 0.21132486540518713 * tmp_12 + tmp_13;
      real_t tmp_15 = tmp_10 * tmp_14;
      real_t tmp_16 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_17 = p_affine_6_0 + tmp_2;
      real_t tmp_18 = 0.21132486540518713 * tmp_16 + tmp_17;
      real_t tmp_19 = tmp_18 * tmp_7;
      real_t tmp_20 = tmp_15 + tmp_19;
      real_t tmp_21 = tmp_14 * tmp_9;
      real_t tmp_22 = tmp_18 * tmp_8;
      real_t tmp_23 = tmp_21 + tmp_22;
      real_t tmp_24 = tmp_1 * ( tmp_23 - 1.0 / 3.0 ) + tmp_5 * ( tmp_20 - 1.0 / 3.0 );
      real_t tmp_25 = p_affine_10_0 * ( tmp_1 * tmp_8 + tmp_5 * tmp_7 ) + p_affine_10_1 * ( tmp_10 * tmp_5 + tmp_4 * tmp_6 );
      real_t tmp_26 = -tmp_15 - tmp_19 - tmp_21 - tmp_22 + 1;
      real_t tmp_27 = std::abs( std::pow( ( tmp_12 * tmp_12 ) + ( tmp_16 * tmp_16 ), 1.0 / 2.0 ) );
      real_t tmp_28 = 24 / tmp_27;
      real_t tmp_29 = tmp_24 * tmp_28;
      real_t tmp_30 = 0.5 * tmp_27;
      real_t tmp_31 = 0.78867513459481287 * tmp_12 + tmp_13;
      real_t tmp_32 = tmp_10 * tmp_31;
      real_t tmp_33 = 0.78867513459481287 * tmp_16 + tmp_17;
      real_t tmp_34 = tmp_33 * tmp_7;
      real_t tmp_35 = tmp_32 + tmp_34;
      real_t tmp_36 = tmp_31 * tmp_9;
      real_t tmp_37 = tmp_33 * tmp_8;
      real_t tmp_38 = tmp_36 + tmp_37;
      real_t tmp_39 = tmp_1 * ( tmp_38 - 1.0 / 3.0 ) + tmp_5 * ( tmp_35 - 1.0 / 3.0 );
      real_t tmp_40 = -tmp_32 - tmp_34 - tmp_36 - tmp_37 + 1;
      real_t tmp_41 = tmp_28 * tmp_39;
      real_t tmp_42 = 0.5 * tmp_27;
      real_t tmp_43 = p_affine_10_0 * tmp_7 + p_affine_10_1 * tmp_10;
      real_t tmp_44 = p_affine_10_0 * tmp_8 + p_affine_10_1 * tmp_9;
      real_t a_0_0  = tmp_30 * ( -tmp_11 * tmp_24 - tmp_25 * tmp_26 + tmp_26 * tmp_29 ) +
                     tmp_42 * ( -tmp_11 * tmp_39 - tmp_25 * tmp_40 + tmp_40 * tmp_41 );
      real_t a_1_0 = tmp_30 * ( -tmp_20 * tmp_25 + tmp_20 * tmp_29 - tmp_24 * tmp_43 ) +
                     tmp_42 * ( -tmp_25 * tmp_35 + tmp_35 * tmp_41 - tmp_39 * tmp_43 );
      real_t a_2_0 = tmp_30 * ( -tmp_23 * tmp_25 + tmp_23 * tmp_29 - tmp_24 * tmp_44 ) +
                     tmp_42 * ( -tmp_25 * tmp_38 + tmp_38 * tmp_41 - tmp_39 * tmp_44 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
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

class DGVectorLaplaceFormEDGP1_0 : public hyteg::dg::DGForm2D
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
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_0;
      real_t tmp_6  = 1.0 / ( tmp_4 - tmp_5 * ( p_affine_1_1 + tmp_2 ) );
      real_t tmp_7  = tmp_1 * tmp_6;
      real_t tmp_8  = tmp_6 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_9  = tmp_1 * tmp_8 + tmp_5 * tmp_7;
      real_t tmp_10 = tmp_3 * tmp_6;
      real_t tmp_11 = tmp_6 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_12 = tmp_11 * tmp_5 + tmp_4 * tmp_6;
      real_t tmp_13 = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                                p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_14 = tmp_13 * ( tmp_12 * ( -tmp_10 - tmp_11 ) + tmp_9 * ( -tmp_7 - tmp_8 ) );
      real_t tmp_15 = tmp_13 * ( tmp_10 * tmp_12 + tmp_8 * tmp_9 );
      real_t tmp_16 = tmp_13 * ( tmp_11 * tmp_12 + tmp_7 * tmp_9 );
      real_t a_0_0  = 0.5 * tmp_14;
      real_t a_0_1  = 0.5 * tmp_15;
      real_t a_0_2  = 0.5 * tmp_16;
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_3  = -p_affine_0_1;
      real_t tmp_4  = p_affine_6_1 + tmp_3;
      real_t tmp_5  = 0.21132486540518713 * tmp_2 + tmp_4;
      real_t tmp_6  = p_affine_2_1 + tmp_3;
      real_t tmp_7  = tmp_1 * tmp_6;
      real_t tmp_8  = p_affine_2_0 + tmp_0;
      real_t tmp_9  = 1.0 / ( tmp_7 - tmp_8 * ( p_affine_1_1 + tmp_3 ) );
      real_t tmp_10 = tmp_9 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_11 = tmp_10 * tmp_5;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_0;
      real_t tmp_14 = 0.21132486540518713 * tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6 * tmp_9;
      real_t tmp_16 = tmp_14 * tmp_15;
      real_t tmp_17 = tmp_11 + tmp_16;
      real_t tmp_18 = tmp_1 * tmp_9;
      real_t tmp_19 = tmp_18 * tmp_5;
      real_t tmp_20 = tmp_9 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_21 = tmp_14 * tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_1 * ( tmp_17 - 1.0 / 3.0 ) + tmp_8 * ( tmp_22 - 1.0 / 3.0 );
      real_t tmp_24 = 0.5 * p_affine_10_0 * ( -tmp_15 - tmp_20 ) + 0.5 * p_affine_10_1 * ( -tmp_10 - tmp_18 );
      real_t tmp_25 = -tmp_11 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_26 =
          0.5 * p_affine_10_0 * ( tmp_20 * tmp_8 + tmp_7 * tmp_9 ) + 0.5 * p_affine_10_1 * ( tmp_1 * tmp_10 + tmp_18 * tmp_8 );
      real_t tmp_27 = std::abs( std::pow( ( tmp_12 * tmp_12 ) + ( tmp_2 * tmp_2 ), 1.0 / 2.0 ) );
      real_t tmp_28 = 6 / tmp_27;
      real_t tmp_29 = tmp_23 * tmp_28;
      real_t tmp_30 = 0.5 * tmp_27;
      real_t tmp_31 = 0.78867513459481287 * tmp_2 + tmp_4;
      real_t tmp_32 = tmp_10 * tmp_31;
      real_t tmp_33 = 0.78867513459481287 * tmp_12 + tmp_13;
      real_t tmp_34 = tmp_15 * tmp_33;
      real_t tmp_35 = tmp_32 + tmp_34;
      real_t tmp_36 = tmp_18 * tmp_31;
      real_t tmp_37 = tmp_20 * tmp_33;
      real_t tmp_38 = tmp_36 + tmp_37;
      real_t tmp_39 = tmp_1 * ( tmp_35 - 1.0 / 3.0 ) + tmp_8 * ( tmp_38 - 1.0 / 3.0 );
      real_t tmp_40 = -tmp_32 - tmp_34 - tmp_36 - tmp_37 + 1;
      real_t tmp_41 = tmp_28 * tmp_39;
      real_t tmp_42 = 0.5 * tmp_27;
      real_t tmp_43 = 0.5 * p_affine_10_0 * tmp_15 + 0.5 * p_affine_10_1 * tmp_10;
      real_t tmp_44 = 0.5 * p_affine_10_0 * tmp_20 + 0.5 * p_affine_10_1 * tmp_18;
      real_t a_0_0  = tmp_30 * ( -tmp_23 * tmp_24 - tmp_25 * tmp_26 + tmp_25 * tmp_29 ) +
                     tmp_42 * ( -tmp_24 * tmp_39 - tmp_26 * tmp_40 + tmp_40 * tmp_41 );
      real_t a_0_1 = tmp_30 * ( -tmp_17 * tmp_26 + tmp_17 * tmp_29 - tmp_23 * tmp_43 ) +
                     tmp_42 * ( -tmp_26 * tmp_35 + tmp_35 * tmp_41 - tmp_39 * tmp_43 );
      real_t a_0_2 = tmp_30 * ( -tmp_22 * tmp_26 + tmp_22 * tmp_29 - tmp_23 * tmp_44 ) +
                     tmp_42 * ( -tmp_26 * tmp_38 + tmp_38 * tmp_41 - tmp_39 * tmp_44 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_4  = p_affine_6_1 + 0.21132486540518713 * tmp_3;
      real_t tmp_5  = tmp_2 + tmp_4;
      real_t tmp_6  = p_affine_2_1 + tmp_2;
      real_t tmp_7  = tmp_1 * tmp_6;
      real_t tmp_8  = p_affine_2_0 + tmp_0;
      real_t tmp_9  = 1.0 / ( tmp_7 - tmp_8 * ( p_affine_1_1 + tmp_2 ) );
      real_t tmp_10 = tmp_9 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.21132486540518713 * tmp_11;
      real_t tmp_13 = tmp_9 * ( tmp_0 + tmp_12 );
      real_t tmp_14 = tmp_1 * tmp_9;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 =
          tmp_1 * ( tmp_10 * tmp_5 + tmp_13 * tmp_6 - 1.0 / 3.0 ) + tmp_8 * ( tmp_13 * tmp_15 + tmp_14 * tmp_5 - 1.0 / 3.0 );
      real_t tmp_17 = -p_affine_3_1;
      real_t tmp_18 = p_affine_5_1 + tmp_17;
      real_t tmp_19 = -p_affine_3_0;
      real_t tmp_20 = p_affine_4_0 + tmp_19;
      real_t tmp_21 = 1.0 / ( tmp_18 * tmp_20 - ( p_affine_4_1 + tmp_17 ) * ( p_affine_5_0 + tmp_19 ) );
      real_t tmp_22 = tmp_18 * tmp_21;
      real_t tmp_23 = tmp_21 * ( p_affine_3_1 - p_affine_4_1 );
      real_t tmp_24 = tmp_20 * tmp_21;
      real_t tmp_25 = tmp_21 * ( p_affine_3_0 - p_affine_5_0 );
      real_t tmp_26 = 0.5 * p_affine_10_0 * ( -tmp_22 - tmp_23 ) + 0.5 * p_affine_10_1 * ( -tmp_24 - tmp_25 );
      real_t tmp_27 = tmp_17 + tmp_4;
      real_t tmp_28 = tmp_24 * tmp_27;
      real_t tmp_29 = tmp_25 * tmp_27;
      real_t tmp_30 = tmp_12 + tmp_19;
      real_t tmp_31 = tmp_22 * tmp_30;
      real_t tmp_32 = tmp_23 * tmp_30;
      real_t tmp_33 = -tmp_28 - tmp_29 - tmp_31 - tmp_32 + 1;
      real_t tmp_34 = tmp_8 * tmp_9;
      real_t tmp_35 =
          0.5 * p_affine_10_0 * ( tmp_15 * tmp_34 + tmp_7 * tmp_9 ) + 0.5 * p_affine_10_1 * ( tmp_1 * tmp_10 + tmp_1 * tmp_34 );
      real_t tmp_36 = std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_3 * tmp_3 ), 1.0 / 2.0 ) );
      real_t tmp_37 = 6 / tmp_36;
      real_t tmp_38 = tmp_16 * tmp_37;
      real_t tmp_39 = 0.5 * tmp_36;
      real_t tmp_40 = p_affine_6_1 + 0.78867513459481287 * tmp_3;
      real_t tmp_41 = tmp_2 + tmp_40;
      real_t tmp_42 = p_affine_6_0 + 0.78867513459481287 * tmp_11;
      real_t tmp_43 = tmp_9 * ( tmp_0 + tmp_42 );
      real_t tmp_44 =
          tmp_1 * ( tmp_10 * tmp_41 + tmp_43 * tmp_6 - 1.0 / 3.0 ) + tmp_8 * ( tmp_14 * tmp_41 + tmp_15 * tmp_43 - 1.0 / 3.0 );
      real_t tmp_45 = tmp_17 + tmp_40;
      real_t tmp_46 = tmp_24 * tmp_45;
      real_t tmp_47 = tmp_25 * tmp_45;
      real_t tmp_48 = tmp_19 + tmp_42;
      real_t tmp_49 = tmp_22 * tmp_48;
      real_t tmp_50 = tmp_23 * tmp_48;
      real_t tmp_51 = -tmp_46 - tmp_47 - tmp_49 - tmp_50 + 1;
      real_t tmp_52 = tmp_37 * tmp_44;
      real_t tmp_53 = 0.5 * tmp_36;
      real_t tmp_54 = tmp_29 + tmp_31;
      real_t tmp_55 = 0.5 * p_affine_10_0 * tmp_22 + 0.5 * p_affine_10_1 * tmp_25;
      real_t tmp_56 = tmp_47 + tmp_49;
      real_t tmp_57 = tmp_28 + tmp_32;
      real_t tmp_58 = 0.5 * p_affine_10_0 * tmp_23 + 0.5 * p_affine_10_1 * tmp_24;
      real_t tmp_59 = tmp_46 + tmp_50;
      real_t a_0_0  = tmp_39 * ( -tmp_16 * tmp_26 + tmp_33 * tmp_35 - tmp_33 * tmp_38 ) +
                     tmp_53 * ( -tmp_26 * tmp_44 + tmp_35 * tmp_51 - tmp_51 * tmp_52 );
      real_t a_0_1 = tmp_39 * ( -tmp_16 * tmp_55 + tmp_35 * tmp_54 - tmp_38 * tmp_54 ) +
                     tmp_53 * ( tmp_35 * tmp_56 - tmp_44 * tmp_55 - tmp_52 * tmp_56 );
      real_t a_0_2 = tmp_39 * ( -tmp_16 * tmp_58 + tmp_35 * tmp_57 - tmp_38 * tmp_57 ) +
                     tmp_53 * ( tmp_35 * tmp_59 - tmp_44 * tmp_58 - tmp_52 * tmp_59 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
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
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_2;
      real_t tmp_6  = 1.0 / ( tmp_4 - tmp_5 * ( p_affine_1_1 + tmp_0 ) );
      real_t tmp_7  = tmp_1 * tmp_6;
      real_t tmp_8  = tmp_6 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_9  = tmp_3 * tmp_6;
      real_t tmp_10 = tmp_6 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_11 = p_affine_10_0 * ( -tmp_7 - tmp_8 ) + p_affine_10_1 * ( -tmp_10 - tmp_9 );
      real_t tmp_12 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_13 = p_affine_6_1 + tmp_0;
      real_t tmp_14 = 0.21132486540518713 * tmp_12 + tmp_13;
      real_t tmp_15 = tmp_10 * tmp_14;
      real_t tmp_16 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_17 = p_affine_6_0 + tmp_2;
      real_t tmp_18 = 0.21132486540518713 * tmp_16 + tmp_17;
      real_t tmp_19 = tmp_18 * tmp_7;
      real_t tmp_20 = tmp_15 + tmp_19;
      real_t tmp_21 = tmp_14 * tmp_9;
      real_t tmp_22 = tmp_18 * tmp_8;
      real_t tmp_23 = tmp_21 + tmp_22;
      real_t tmp_24 = tmp_3 * ( tmp_20 - 1.0 / 3.0 ) + tmp_5 * ( tmp_23 - 1.0 / 3.0 );
      real_t tmp_25 = p_affine_10_0 * ( tmp_4 * tmp_6 + tmp_5 * tmp_8 ) + p_affine_10_1 * ( tmp_10 * tmp_3 + tmp_5 * tmp_9 );
      real_t tmp_26 = -tmp_15 - tmp_19 - tmp_21 - tmp_22 + 1;
      real_t tmp_27 = std::abs( std::pow( ( tmp_12 * tmp_12 ) + ( tmp_16 * tmp_16 ), 1.0 / 2.0 ) );
      real_t tmp_28 = 24 / tmp_27;
      real_t tmp_29 = tmp_24 * tmp_28;
      real_t tmp_30 = 0.5 * tmp_27;
      real_t tmp_31 = 0.78867513459481287 * tmp_12 + tmp_13;
      real_t tmp_32 = tmp_10 * tmp_31;
      real_t tmp_33 = 0.78867513459481287 * tmp_16 + tmp_17;
      real_t tmp_34 = tmp_33 * tmp_7;
      real_t tmp_35 = tmp_32 + tmp_34;
      real_t tmp_36 = tmp_31 * tmp_9;
      real_t tmp_37 = tmp_33 * tmp_8;
      real_t tmp_38 = tmp_36 + tmp_37;
      real_t tmp_39 = tmp_3 * ( tmp_35 - 1.0 / 3.0 ) + tmp_5 * ( tmp_38 - 1.0 / 3.0 );
      real_t tmp_40 = -tmp_32 - tmp_34 - tmp_36 - tmp_37 + 1;
      real_t tmp_41 = tmp_28 * tmp_39;
      real_t tmp_42 = 0.5 * tmp_27;
      real_t tmp_43 = p_affine_10_0 * tmp_7 + p_affine_10_1 * tmp_10;
      real_t tmp_44 = p_affine_10_0 * tmp_8 + p_affine_10_1 * tmp_9;
      real_t a_0_0  = tmp_30 * ( -tmp_11 * tmp_24 - tmp_25 * tmp_26 + tmp_26 * tmp_29 ) +
                     tmp_42 * ( -tmp_11 * tmp_39 - tmp_25 * tmp_40 + tmp_40 * tmp_41 );
      real_t a_0_1 = tmp_30 * ( -tmp_20 * tmp_25 + tmp_20 * tmp_29 - tmp_24 * tmp_43 ) +
                     tmp_42 * ( -tmp_25 * tmp_35 + tmp_35 * tmp_41 - tmp_39 * tmp_43 );
      real_t a_0_2 = tmp_30 * ( -tmp_23 * tmp_25 + tmp_23 * tmp_29 - tmp_24 * tmp_44 ) +
                     tmp_42 * ( -tmp_25 * tmp_38 + tmp_38 * tmp_41 - tmp_39 * tmp_44 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
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

class DGVectorLaplaceFormEDGP1_1 : public hyteg::dg::DGForm2D
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
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_1_1 + tmp_2;
      real_t tmp_6  = 1.0 / ( tmp_4 - tmp_5 * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_7  = tmp_1 * tmp_6;
      real_t tmp_8  = tmp_6 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_9  = tmp_4 * tmp_6 + tmp_5 * tmp_8;
      real_t tmp_10 = tmp_3 * tmp_6;
      real_t tmp_11 = tmp_6 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_12 = tmp_10 * tmp_5 + tmp_11 * tmp_3;
      real_t tmp_13 = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                                p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_14 = tmp_13 * ( tmp_12 * ( -tmp_10 - tmp_11 ) + tmp_9 * ( -tmp_7 - tmp_8 ) );
      real_t tmp_15 = tmp_13 * ( tmp_10 * tmp_12 + tmp_8 * tmp_9 );
      real_t tmp_16 = tmp_13 * ( tmp_11 * tmp_12 + tmp_7 * tmp_9 );
      real_t a_0_0  = 0.5 * tmp_14;
      real_t a_0_1  = 0.5 * tmp_15;
      real_t a_0_2  = 0.5 * tmp_16;
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

      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_1_1 + tmp_0;
      real_t tmp_2  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_3  = p_affine_6_1 + tmp_0;
      real_t tmp_4  = 0.21132486540518713 * tmp_2 + tmp_3;
      real_t tmp_5  = -p_affine_0_0;
      real_t tmp_6  = p_affine_1_0 + tmp_5;
      real_t tmp_7  = p_affine_2_1 + tmp_0;
      real_t tmp_8  = tmp_6 * tmp_7;
      real_t tmp_9  = 1.0 / ( -tmp_1 * ( p_affine_2_0 + tmp_5 ) + tmp_8 );
      real_t tmp_10 = tmp_9 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_11 = tmp_10 * tmp_4;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_5;
      real_t tmp_14 = 0.21132486540518713 * tmp_12 + tmp_13;
      real_t tmp_15 = tmp_7 * tmp_9;
      real_t tmp_16 = tmp_14 * tmp_15;
      real_t tmp_17 = tmp_11 + tmp_16;
      real_t tmp_18 = tmp_6 * tmp_9;
      real_t tmp_19 = tmp_18 * tmp_4;
      real_t tmp_20 = tmp_9 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_21 = tmp_14 * tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_1 * ( tmp_17 - 1.0 / 3.0 ) + tmp_7 * ( tmp_22 - 1.0 / 3.0 );
      real_t tmp_24 = 0.5 * p_affine_10_0 * ( -tmp_15 - tmp_20 ) + 0.5 * p_affine_10_1 * ( -tmp_10 - tmp_18 );
      real_t tmp_25 = -tmp_11 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_26 =
          0.5 * p_affine_10_0 * ( tmp_1 * tmp_15 + tmp_20 * tmp_7 ) + 0.5 * p_affine_10_1 * ( tmp_1 * tmp_10 + tmp_8 * tmp_9 );
      real_t tmp_27 = std::abs( std::pow( ( tmp_12 * tmp_12 ) + ( tmp_2 * tmp_2 ), 1.0 / 2.0 ) );
      real_t tmp_28 = 6 / tmp_27;
      real_t tmp_29 = tmp_23 * tmp_28;
      real_t tmp_30 = 0.5 * tmp_27;
      real_t tmp_31 = 0.78867513459481287 * tmp_2 + tmp_3;
      real_t tmp_32 = tmp_10 * tmp_31;
      real_t tmp_33 = 0.78867513459481287 * tmp_12 + tmp_13;
      real_t tmp_34 = tmp_15 * tmp_33;
      real_t tmp_35 = tmp_32 + tmp_34;
      real_t tmp_36 = tmp_18 * tmp_31;
      real_t tmp_37 = tmp_20 * tmp_33;
      real_t tmp_38 = tmp_36 + tmp_37;
      real_t tmp_39 = tmp_1 * ( tmp_35 - 1.0 / 3.0 ) + tmp_7 * ( tmp_38 - 1.0 / 3.0 );
      real_t tmp_40 = -tmp_32 - tmp_34 - tmp_36 - tmp_37 + 1;
      real_t tmp_41 = tmp_28 * tmp_39;
      real_t tmp_42 = 0.5 * tmp_27;
      real_t tmp_43 = 0.5 * p_affine_10_0 * tmp_15 + 0.5 * p_affine_10_1 * tmp_10;
      real_t tmp_44 = 0.5 * p_affine_10_0 * tmp_20 + 0.5 * p_affine_10_1 * tmp_18;
      real_t a_0_0  = tmp_30 * ( -tmp_23 * tmp_24 - tmp_25 * tmp_26 + tmp_25 * tmp_29 ) +
                     tmp_42 * ( -tmp_24 * tmp_39 - tmp_26 * tmp_40 + tmp_40 * tmp_41 );
      real_t a_0_1 = tmp_30 * ( -tmp_17 * tmp_26 + tmp_17 * tmp_29 - tmp_23 * tmp_43 ) +
                     tmp_42 * ( -tmp_26 * tmp_35 + tmp_35 * tmp_41 - tmp_39 * tmp_43 );
      real_t a_0_2 = tmp_30 * ( -tmp_22 * tmp_26 + tmp_22 * tmp_29 - tmp_23 * tmp_44 ) +
                     tmp_42 * ( -tmp_26 * tmp_38 + tmp_38 * tmp_41 - tmp_39 * tmp_44 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
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

      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_1_1 + tmp_0;
      real_t tmp_2  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_3  = p_affine_6_1 + 0.21132486540518713 * tmp_2;
      real_t tmp_4  = tmp_0 + tmp_3;
      real_t tmp_5  = -p_affine_0_0;
      real_t tmp_6  = p_affine_1_0 + tmp_5;
      real_t tmp_7  = p_affine_2_1 + tmp_0;
      real_t tmp_8  = tmp_6 * tmp_7;
      real_t tmp_9  = 1.0 / ( -tmp_1 * ( p_affine_2_0 + tmp_5 ) + tmp_8 );
      real_t tmp_10 = tmp_9 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.21132486540518713 * tmp_11;
      real_t tmp_13 = tmp_12 + tmp_5;
      real_t tmp_14 = tmp_7 * tmp_9;
      real_t tmp_15 = tmp_6 * tmp_9;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_16 * tmp_9;
      real_t tmp_18 =
          tmp_1 * ( tmp_10 * tmp_4 + tmp_13 * tmp_14 - 1.0 / 3.0 ) + tmp_7 * ( tmp_13 * tmp_17 + tmp_15 * tmp_4 - 1.0 / 3.0 );
      real_t tmp_19 = -p_affine_3_1;
      real_t tmp_20 = p_affine_5_1 + tmp_19;
      real_t tmp_21 = -p_affine_3_0;
      real_t tmp_22 = p_affine_4_0 + tmp_21;
      real_t tmp_23 = 1.0 / ( tmp_20 * tmp_22 - ( p_affine_4_1 + tmp_19 ) * ( p_affine_5_0 + tmp_21 ) );
      real_t tmp_24 = tmp_20 * tmp_23;
      real_t tmp_25 = tmp_23 * ( p_affine_3_1 - p_affine_4_1 );
      real_t tmp_26 = tmp_22 * tmp_23;
      real_t tmp_27 = tmp_23 * ( p_affine_3_0 - p_affine_5_0 );
      real_t tmp_28 = 0.5 * p_affine_10_0 * ( -tmp_24 - tmp_25 ) + 0.5 * p_affine_10_1 * ( -tmp_26 - tmp_27 );
      real_t tmp_29 = tmp_19 + tmp_3;
      real_t tmp_30 = tmp_26 * tmp_29;
      real_t tmp_31 = tmp_27 * tmp_29;
      real_t tmp_32 = tmp_12 + tmp_21;
      real_t tmp_33 = tmp_24 * tmp_32;
      real_t tmp_34 = tmp_25 * tmp_32;
      real_t tmp_35 = -tmp_30 - tmp_31 - tmp_33 - tmp_34 + 1;
      real_t tmp_36 =
          0.5 * p_affine_10_0 * ( tmp_1 * tmp_14 + tmp_14 * tmp_16 ) + 0.5 * p_affine_10_1 * ( tmp_1 * tmp_10 + tmp_8 * tmp_9 );
      real_t tmp_37 = std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_2 * tmp_2 ), 1.0 / 2.0 ) );
      real_t tmp_38 = 6 / tmp_37;
      real_t tmp_39 = tmp_18 * tmp_38;
      real_t tmp_40 = 0.5 * tmp_37;
      real_t tmp_41 = p_affine_6_1 + 0.78867513459481287 * tmp_2;
      real_t tmp_42 = tmp_0 + tmp_41;
      real_t tmp_43 = p_affine_6_0 + 0.78867513459481287 * tmp_11;
      real_t tmp_44 = tmp_43 + tmp_5;
      real_t tmp_45 =
          tmp_1 * ( tmp_10 * tmp_42 + tmp_14 * tmp_44 - 1.0 / 3.0 ) + tmp_7 * ( tmp_15 * tmp_42 + tmp_17 * tmp_44 - 1.0 / 3.0 );
      real_t tmp_46 = tmp_19 + tmp_41;
      real_t tmp_47 = tmp_26 * tmp_46;
      real_t tmp_48 = tmp_27 * tmp_46;
      real_t tmp_49 = tmp_21 + tmp_43;
      real_t tmp_50 = tmp_24 * tmp_49;
      real_t tmp_51 = tmp_25 * tmp_49;
      real_t tmp_52 = -tmp_47 - tmp_48 - tmp_50 - tmp_51 + 1;
      real_t tmp_53 = tmp_38 * tmp_45;
      real_t tmp_54 = 0.5 * tmp_37;
      real_t tmp_55 = tmp_31 + tmp_33;
      real_t tmp_56 = 0.5 * p_affine_10_0 * tmp_24 + 0.5 * p_affine_10_1 * tmp_27;
      real_t tmp_57 = tmp_48 + tmp_50;
      real_t tmp_58 = tmp_30 + tmp_34;
      real_t tmp_59 = 0.5 * p_affine_10_0 * tmp_25 + 0.5 * p_affine_10_1 * tmp_26;
      real_t tmp_60 = tmp_47 + tmp_51;
      real_t a_0_0  = tmp_40 * ( -tmp_18 * tmp_28 + tmp_35 * tmp_36 - tmp_35 * tmp_39 ) +
                     tmp_54 * ( -tmp_28 * tmp_45 + tmp_36 * tmp_52 - tmp_52 * tmp_53 );
      real_t a_0_1 = tmp_40 * ( -tmp_18 * tmp_56 + tmp_36 * tmp_55 - tmp_39 * tmp_55 ) +
                     tmp_54 * ( tmp_36 * tmp_57 - tmp_45 * tmp_56 - tmp_53 * tmp_57 );
      real_t a_0_2 = tmp_40 * ( -tmp_18 * tmp_59 + tmp_36 * tmp_58 - tmp_39 * tmp_58 ) +
                     tmp_54 * ( tmp_36 * tmp_60 - tmp_45 * tmp_59 - tmp_53 * tmp_60 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
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
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_1_1 + tmp_0;
      real_t tmp_6  = 1.0 / ( tmp_4 - tmp_5 * ( p_affine_2_0 + tmp_2 ) );
      real_t tmp_7  = tmp_1 * tmp_6;
      real_t tmp_8  = tmp_6 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_9  = tmp_3 * tmp_6;
      real_t tmp_10 = tmp_6 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_11 = p_affine_10_0 * ( -tmp_7 - tmp_8 ) + p_affine_10_1 * ( -tmp_10 - tmp_9 );
      real_t tmp_12 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_13 = p_affine_6_1 + tmp_0;
      real_t tmp_14 = 0.21132486540518713 * tmp_12 + tmp_13;
      real_t tmp_15 = tmp_10 * tmp_14;
      real_t tmp_16 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_17 = p_affine_6_0 + tmp_2;
      real_t tmp_18 = 0.21132486540518713 * tmp_16 + tmp_17;
      real_t tmp_19 = tmp_18 * tmp_7;
      real_t tmp_20 = tmp_15 + tmp_19;
      real_t tmp_21 = tmp_14 * tmp_9;
      real_t tmp_22 = tmp_18 * tmp_8;
      real_t tmp_23 = tmp_21 + tmp_22;
      real_t tmp_24 = tmp_1 * ( tmp_23 - 1.0 / 3.0 ) + tmp_5 * ( tmp_20 - 1.0 / 3.0 );
      real_t tmp_25 = p_affine_10_0 * ( tmp_1 * tmp_8 + tmp_5 * tmp_7 ) + p_affine_10_1 * ( tmp_10 * tmp_5 + tmp_4 * tmp_6 );
      real_t tmp_26 = -tmp_15 - tmp_19 - tmp_21 - tmp_22 + 1;
      real_t tmp_27 = std::abs( std::pow( ( tmp_12 * tmp_12 ) + ( tmp_16 * tmp_16 ), 1.0 / 2.0 ) );
      real_t tmp_28 = 24 / tmp_27;
      real_t tmp_29 = tmp_24 * tmp_28;
      real_t tmp_30 = 0.5 * tmp_27;
      real_t tmp_31 = 0.78867513459481287 * tmp_12 + tmp_13;
      real_t tmp_32 = tmp_10 * tmp_31;
      real_t tmp_33 = 0.78867513459481287 * tmp_16 + tmp_17;
      real_t tmp_34 = tmp_33 * tmp_7;
      real_t tmp_35 = tmp_32 + tmp_34;
      real_t tmp_36 = tmp_31 * tmp_9;
      real_t tmp_37 = tmp_33 * tmp_8;
      real_t tmp_38 = tmp_36 + tmp_37;
      real_t tmp_39 = tmp_1 * ( tmp_38 - 1.0 / 3.0 ) + tmp_5 * ( tmp_35 - 1.0 / 3.0 );
      real_t tmp_40 = -tmp_32 - tmp_34 - tmp_36 - tmp_37 + 1;
      real_t tmp_41 = tmp_28 * tmp_39;
      real_t tmp_42 = 0.5 * tmp_27;
      real_t tmp_43 = p_affine_10_0 * tmp_7 + p_affine_10_1 * tmp_10;
      real_t tmp_44 = p_affine_10_0 * tmp_8 + p_affine_10_1 * tmp_9;
      real_t a_0_0  = tmp_30 * ( -tmp_11 * tmp_24 - tmp_25 * tmp_26 + tmp_26 * tmp_29 ) +
                     tmp_42 * ( -tmp_11 * tmp_39 - tmp_25 * tmp_40 + tmp_40 * tmp_41 );
      real_t a_0_1 = tmp_30 * ( -tmp_20 * tmp_25 + tmp_20 * tmp_29 - tmp_24 * tmp_43 ) +
                     tmp_42 * ( -tmp_25 * tmp_35 + tmp_35 * tmp_41 - tmp_39 * tmp_43 );
      real_t a_0_2 = tmp_30 * ( -tmp_23 * tmp_25 + tmp_23 * tmp_29 - tmp_24 * tmp_44 ) +
                     tmp_42 * ( -tmp_25 * tmp_38 + tmp_38 * tmp_41 - tmp_39 * tmp_44 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
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

class DGVectorLaplaceFormEDGEDG : public hyteg::dg::DGForm2D
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
      real_t tmp_1  = p_affine_2_0 + tmp_0;
      real_t tmp_2  = p_affine_1_0 + tmp_0;
      real_t tmp_3  = -p_affine_0_1;
      real_t tmp_4  = p_affine_2_1 + tmp_3;
      real_t tmp_5  = tmp_2 * tmp_4;
      real_t tmp_6  = p_affine_1_1 + tmp_3;
      real_t tmp_7  = 1.0 / ( -tmp_1 * tmp_6 + tmp_5 );
      real_t tmp_8  = tmp_2 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_5 * tmp_7;
      real_t tmp_11 = tmp_7 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_12 = tmp_6 * tmp_7;
      real_t tmp_13 = ( ( ( tmp_10 + tmp_12 * tmp_9 ) * ( tmp_10 + tmp_12 * tmp_9 ) ) +
                        ( ( tmp_1 * tmp_11 + tmp_10 ) * ( tmp_1 * tmp_11 + tmp_10 ) ) +
                        ( ( tmp_1 * tmp_8 + tmp_8 * tmp_9 ) * ( tmp_1 * tmp_8 + tmp_8 * tmp_9 ) ) +
                        ( ( tmp_11 * tmp_4 + tmp_12 * tmp_4 ) * ( tmp_11 * tmp_4 + tmp_12 * tmp_4 ) ) ) *
                      std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                                p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t a_0_0  = 0.5 * tmp_13;
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

      real_t tmp_0  = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2  = std::abs( std::pow( ( tmp_0 * tmp_0 ) + ( tmp_1 * tmp_1 ), 1.0 / 2.0 ) );
      real_t tmp_3  = -p_affine_0_0;
      real_t tmp_4  = p_affine_1_0 + tmp_3;
      real_t tmp_5  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_6  = -p_affine_0_1;
      real_t tmp_7  = p_affine_2_1 + tmp_6;
      real_t tmp_8  = tmp_4 * tmp_7;
      real_t tmp_9  = p_affine_2_0 + tmp_3;
      real_t tmp_10 = p_affine_1_1 + tmp_6;
      real_t tmp_11 = 1.0 / ( -tmp_10 * tmp_9 + tmp_8 );
      real_t tmp_12 = p_affine_6_1 + tmp_6;
      real_t tmp_13 = tmp_11 * ( 0.21132486540518713 * tmp_1 + tmp_12 );
      real_t tmp_14 = p_affine_6_0 + tmp_3;
      real_t tmp_15 = tmp_11 * ( 0.21132486540518713 * tmp_0 + tmp_14 );
      real_t tmp_16 = tmp_13 * tmp_5 + tmp_15 * tmp_7 - 1.0 / 3.0;
      real_t tmp_17 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_18 = tmp_13 * tmp_4 + tmp_15 * tmp_17 - 1.0 / 3.0;
      real_t tmp_19 = tmp_16 * tmp_4 + tmp_18 * tmp_9;
      real_t tmp_20 = tmp_11 * tmp_8;
      real_t tmp_21 = tmp_11 * tmp_9;
      real_t tmp_22 = tmp_11 * tmp_5;
      real_t tmp_23 =
          1.0 * p_affine_10_0 * ( tmp_17 * tmp_21 + tmp_20 ) + 1.0 * p_affine_10_1 * ( tmp_21 * tmp_4 + tmp_22 * tmp_4 );
      real_t tmp_24 = tmp_10 * tmp_16 + tmp_18 * tmp_7;
      real_t tmp_25 = tmp_11 * tmp_7;
      real_t tmp_26 =
          1.0 * p_affine_10_0 * ( tmp_10 * tmp_25 + tmp_17 * tmp_25 ) + 1.0 * p_affine_10_1 * ( tmp_10 * tmp_22 + tmp_20 );
      real_t tmp_27 = 6 / tmp_2;
      real_t tmp_28 = 0.78867513459481287 * tmp_1 + tmp_12;
      real_t tmp_29 = 0.78867513459481287 * tmp_0 + tmp_14;
      real_t tmp_30 = tmp_22 * tmp_28 + tmp_25 * tmp_29 - 1.0 / 3.0;
      real_t tmp_31 = tmp_11 * tmp_17 * tmp_29 + tmp_11 * tmp_28 * tmp_4 - 1.0 / 3.0;
      real_t tmp_32 = tmp_30 * tmp_4 + tmp_31 * tmp_9;
      real_t tmp_33 = tmp_10 * tmp_30 + tmp_31 * tmp_7;
      real_t a_0_0 =
          0.5 * tmp_2 * ( -tmp_19 * tmp_23 - tmp_24 * tmp_26 + tmp_27 * ( ( tmp_19 * tmp_19 ) + ( tmp_24 * tmp_24 ) ) ) +
          0.5 * tmp_2 * ( -tmp_23 * tmp_32 - tmp_26 * tmp_33 + tmp_27 * ( ( tmp_32 * tmp_32 ) + ( tmp_33 * tmp_33 ) ) );
      elMat( 0, 0 ) = a_0_0;
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

      real_t tmp_0  = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2  = std::abs( std::pow( ( tmp_0 * tmp_0 ) + ( tmp_1 * tmp_1 ), 1.0 / 2.0 ) );
      real_t tmp_3  = -p_affine_3_0;
      real_t tmp_4  = p_affine_4_0 + tmp_3;
      real_t tmp_5  = p_affine_3_0 - p_affine_5_0;
      real_t tmp_6  = -p_affine_3_1;
      real_t tmp_7  = p_affine_5_1 + tmp_6;
      real_t tmp_8  = tmp_4 * tmp_7;
      real_t tmp_9  = p_affine_5_0 + tmp_3;
      real_t tmp_10 = p_affine_4_1 + tmp_6;
      real_t tmp_11 = 1.0 / ( -tmp_10 * tmp_9 + tmp_8 );
      real_t tmp_12 = p_affine_6_1 + 0.21132486540518713 * tmp_1;
      real_t tmp_13 = tmp_11 * ( tmp_12 + tmp_6 );
      real_t tmp_14 = p_affine_6_0 + 0.21132486540518713 * tmp_0;
      real_t tmp_15 = tmp_11 * ( tmp_14 + tmp_3 );
      real_t tmp_16 = tmp_13 * tmp_5 + tmp_15 * tmp_7 - 1.0 / 3.0;
      real_t tmp_17 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_18 = tmp_13 * tmp_4 + tmp_15 * tmp_17 - 1.0 / 3.0;
      real_t tmp_19 = tmp_16 * tmp_4 + tmp_18 * tmp_9;
      real_t tmp_20 = -p_affine_0_0;
      real_t tmp_21 = p_affine_1_0 + tmp_20;
      real_t tmp_22 = -p_affine_0_1;
      real_t tmp_23 = p_affine_2_1 + tmp_22;
      real_t tmp_24 = tmp_21 * tmp_23;
      real_t tmp_25 = p_affine_2_0 + tmp_20;
      real_t tmp_26 = p_affine_1_1 + tmp_22;
      real_t tmp_27 = 1.0 / ( tmp_24 - tmp_25 * tmp_26 );
      real_t tmp_28 = tmp_24 * tmp_27;
      real_t tmp_29 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_30 = tmp_25 * tmp_27;
      real_t tmp_31 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_32 = tmp_27 * tmp_31;
      real_t tmp_33 =
          0.5 * p_affine_10_0 * ( tmp_28 + tmp_29 * tmp_30 ) + 0.5 * p_affine_10_1 * ( tmp_21 * tmp_30 + tmp_21 * tmp_32 );
      real_t tmp_34 = tmp_10 * tmp_16 + tmp_18 * tmp_7;
      real_t tmp_35 = tmp_23 * tmp_27;
      real_t tmp_36 =
          0.5 * p_affine_10_0 * ( tmp_26 * tmp_35 + tmp_29 * tmp_35 ) + 0.5 * p_affine_10_1 * ( tmp_26 * tmp_32 + tmp_28 );
      real_t tmp_37 = tmp_27 * ( tmp_12 + tmp_22 );
      real_t tmp_38 = tmp_27 * ( tmp_14 + tmp_20 );
      real_t tmp_39 = tmp_23 * tmp_38 + tmp_31 * tmp_37 - 1.0 / 3.0;
      real_t tmp_40 = tmp_21 * tmp_37 + tmp_29 * tmp_38 - 1.0 / 3.0;
      real_t tmp_41 = tmp_21 * tmp_39 + tmp_25 * tmp_40;
      real_t tmp_42 = tmp_11 * tmp_8;
      real_t tmp_43 = tmp_11 * tmp_9;
      real_t tmp_44 = tmp_11 * tmp_5;
      real_t tmp_45 =
          0.5 * p_affine_10_0 * ( tmp_17 * tmp_43 + tmp_42 ) + 0.5 * p_affine_10_1 * ( tmp_4 * tmp_43 + tmp_4 * tmp_44 );
      real_t tmp_46 = tmp_23 * tmp_40 + tmp_26 * tmp_39;
      real_t tmp_47 = tmp_11 * tmp_7;
      real_t tmp_48 =
          0.5 * p_affine_10_0 * ( tmp_10 * tmp_47 + tmp_17 * tmp_47 ) + 0.5 * p_affine_10_1 * ( tmp_10 * tmp_44 + tmp_42 );
      real_t tmp_49 = 6 / tmp_2;
      real_t tmp_50 = p_affine_6_1 + 0.78867513459481287 * tmp_1;
      real_t tmp_51 = tmp_50 + tmp_6;
      real_t tmp_52 = p_affine_6_0 + 0.78867513459481287 * tmp_0;
      real_t tmp_53 = tmp_3 + tmp_52;
      real_t tmp_54 = tmp_44 * tmp_51 + tmp_47 * tmp_53 - 1.0 / 3.0;
      real_t tmp_55 = tmp_11 * tmp_17 * tmp_53 + tmp_11 * tmp_4 * tmp_51 - 1.0 / 3.0;
      real_t tmp_56 = tmp_4 * tmp_54 + tmp_55 * tmp_9;
      real_t tmp_57 = tmp_10 * tmp_54 + tmp_55 * tmp_7;
      real_t tmp_58 = tmp_22 + tmp_50;
      real_t tmp_59 = tmp_20 + tmp_52;
      real_t tmp_60 = tmp_32 * tmp_58 + tmp_35 * tmp_59 - 1.0 / 3.0;
      real_t tmp_61 = tmp_21 * tmp_27 * tmp_58 + tmp_27 * tmp_29 * tmp_59 - 1.0 / 3.0;
      real_t tmp_62 = tmp_21 * tmp_60 + tmp_25 * tmp_61;
      real_t tmp_63 = tmp_23 * tmp_61 + tmp_26 * tmp_60;
      real_t a_0_0  = 0.5 * tmp_2 *
                         ( tmp_19 * tmp_33 + tmp_34 * tmp_36 - tmp_41 * tmp_45 - tmp_46 * tmp_48 -
                           tmp_49 * ( tmp_19 * tmp_41 + tmp_34 * tmp_46 ) ) +
                     0.5 * tmp_2 *
                         ( tmp_33 * tmp_56 + tmp_36 * tmp_57 - tmp_45 * tmp_62 - tmp_48 * tmp_63 -
                           tmp_49 * ( tmp_56 * tmp_62 + tmp_57 * tmp_63 ) );
      elMat( 0, 0 ) = a_0_0;
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

      real_t tmp_0  = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_1  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_2  = std::abs( std::pow( ( tmp_0 * tmp_0 ) + ( tmp_1 * tmp_1 ), 1.0 / 2.0 ) );
      real_t tmp_3  = -p_affine_0_0;
      real_t tmp_4  = p_affine_1_0 + tmp_3;
      real_t tmp_5  = -p_affine_0_1;
      real_t tmp_6  = p_affine_6_1 + tmp_5;
      real_t tmp_7  = 0.21132486540518713 * tmp_1 + tmp_6;
      real_t tmp_8  = p_affine_2_1 + tmp_5;
      real_t tmp_9  = tmp_4 * tmp_8;
      real_t tmp_10 = p_affine_2_0 + tmp_3;
      real_t tmp_11 = p_affine_1_1 + tmp_5;
      real_t tmp_12 = 1.0 / ( -tmp_10 * tmp_11 + tmp_9 );
      real_t tmp_13 = tmp_12 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_14 = p_affine_6_0 + tmp_3;
      real_t tmp_15 = tmp_12 * ( 0.21132486540518713 * tmp_0 + tmp_14 );
      real_t tmp_16 = tmp_13 * tmp_7 + tmp_15 * tmp_8 - 1.0 / 3.0;
      real_t tmp_17 = tmp_12 * tmp_4;
      real_t tmp_18 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_19 = tmp_15 * tmp_18 + tmp_17 * tmp_7 - 1.0 / 3.0;
      real_t tmp_20 = tmp_10 * tmp_19 + tmp_16 * tmp_4;
      real_t tmp_21 = tmp_12 * tmp_9;
      real_t tmp_22 = tmp_10 * tmp_12;
      real_t tmp_23 = 2 * p_affine_10_0 * ( tmp_18 * tmp_22 + tmp_21 ) + 2 * p_affine_10_1 * ( tmp_13 * tmp_4 + tmp_22 * tmp_4 );
      real_t tmp_24 = tmp_11 * tmp_16 + tmp_19 * tmp_8;
      real_t tmp_25 = tmp_12 * tmp_8;
      real_t tmp_26 =
          2 * p_affine_10_0 * ( tmp_11 * tmp_25 + tmp_18 * tmp_25 ) + 2 * p_affine_10_1 * ( tmp_11 * tmp_13 + tmp_21 );
      real_t tmp_27 = 24 / tmp_2;
      real_t tmp_28 = 0.78867513459481287 * tmp_1 + tmp_6;
      real_t tmp_29 = 0.78867513459481287 * tmp_0 + tmp_14;
      real_t tmp_30 = tmp_13 * tmp_28 + tmp_25 * tmp_29 - 1.0 / 3.0;
      real_t tmp_31 = tmp_12 * tmp_18 * tmp_29 + tmp_17 * tmp_28 - 1.0 / 3.0;
      real_t tmp_32 = tmp_10 * tmp_31 + tmp_30 * tmp_4;
      real_t tmp_33 = tmp_11 * tmp_30 + tmp_31 * tmp_8;
      real_t a_0_0 =
          0.5 * tmp_2 * ( -tmp_20 * tmp_23 - tmp_24 * tmp_26 + tmp_27 * ( ( tmp_20 * tmp_20 ) + ( tmp_24 * tmp_24 ) ) ) +
          0.5 * tmp_2 * ( -tmp_23 * tmp_32 - tmp_26 * tmp_33 + tmp_27 * ( ( tmp_32 * tmp_32 ) + ( tmp_33 * tmp_33 ) ) );
      elMat( 0, 0 ) = a_0_0;
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

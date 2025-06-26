/*
* Copyright (c) 2017-2025 Nils Kohl, Marcus Mohr.
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
#include "hyteg/forms/DGForm.hpp"
#include "hyteg/forms/form_hyteg_dg/DGForm2D.hpp"
#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"

#include "Eigen/Eigen"

namespace hyteg {
namespace dg {
namespace eg {

class EGIIPGVectorLaplaceFormP1E_0 : public hyteg::dg::DGForm
{
 protected:
   void integrateVolume2D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override
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

   virtual void integrateFacetInner2D( const std::vector< Point3D >& coordsElement,
                                       const std::vector< Point3D >& coordsFacet,
                                       const Point3D&                oppositeVertex,
                                       const Point3D&                outwardNormal,
                                       const DGBasisInfo&            trialBasis,
                                       const DGBasisInfo&            testBasis,
                                       int                           trialDegree,
                                       int                           testDegree,
                                       MatrixXr&                     elMat ) const override
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
      real_t tmp_7  = tmp_5 * tmp_6;
      real_t tmp_8  = p_affine_2_0 + tmp_4;
      real_t tmp_9  = 1.0 / ( tmp_7 - tmp_8 * ( p_affine_1_1 + tmp_1 ) );
      real_t tmp_10 = tmp_5 * tmp_9;
      real_t tmp_11 = tmp_10 * tmp_3;
      real_t tmp_12 = tmp_9 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_13 = tmp_12 * tmp_3;
      real_t tmp_14 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_15 = p_affine_6_0 + tmp_4;
      real_t tmp_16 = tmp_9 * ( 0.046910077030668018 * tmp_14 + tmp_15 );
      real_t tmp_17 = tmp_16 * tmp_6;
      real_t tmp_18 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_19 = tmp_16 * tmp_18;
      real_t tmp_20 = -tmp_11 - tmp_13 - tmp_17 - tmp_19 + 1;
      real_t tmp_21 = tmp_8 * tmp_9;
      real_t tmp_22 =
          0.5 * p_affine_10_0 * ( tmp_18 * tmp_21 + tmp_7 * tmp_9 ) + 0.5 * p_affine_10_1 * ( tmp_12 * tmp_5 + tmp_21 * tmp_5 );
      real_t tmp_23 = std::abs( std::pow( ( tmp_0 * tmp_0 ) + ( tmp_14 * tmp_14 ), 1.0 / 2.0 ) );
      real_t tmp_24 = 1.0 / ( tmp_23 );
      real_t tmp_25 = tmp_13 + tmp_17;
      real_t tmp_26 = tmp_11 + tmp_19;
      real_t tmp_27 = tmp_24 * ( tmp_5 * ( tmp_25 - 1.0 / 3.0 ) + tmp_8 * ( tmp_26 - 1.0 / 3.0 ) );
      real_t tmp_28 = 0.11846344252809471 * tmp_23;
      real_t tmp_29 = 0.23076534494715845 * tmp_0 + tmp_2;
      real_t tmp_30 = tmp_10 * tmp_29;
      real_t tmp_31 = tmp_12 * tmp_29;
      real_t tmp_32 = tmp_9 * ( 0.23076534494715845 * tmp_14 + tmp_15 );
      real_t tmp_33 = tmp_32 * tmp_6;
      real_t tmp_34 = tmp_18 * tmp_32;
      real_t tmp_35 = -tmp_30 - tmp_31 - tmp_33 - tmp_34 + 1;
      real_t tmp_36 = tmp_31 + tmp_33;
      real_t tmp_37 = tmp_30 + tmp_34;
      real_t tmp_38 = tmp_24 * ( tmp_5 * ( tmp_36 - 1.0 / 3.0 ) + tmp_8 * ( tmp_37 - 1.0 / 3.0 ) );
      real_t tmp_39 = 0.2393143352496831 * tmp_23;
      real_t tmp_40 = 0.5 * tmp_0 + tmp_2;
      real_t tmp_41 = tmp_10 * tmp_40;
      real_t tmp_42 = tmp_12 * tmp_40;
      real_t tmp_43 = tmp_9 * ( 0.5 * tmp_14 + tmp_15 );
      real_t tmp_44 = tmp_43 * tmp_6;
      real_t tmp_45 = tmp_18 * tmp_43;
      real_t tmp_46 = -tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1;
      real_t tmp_47 = tmp_42 + tmp_44;
      real_t tmp_48 = tmp_41 + tmp_45;
      real_t tmp_49 = tmp_24 * ( tmp_5 * ( tmp_47 - 1.0 / 3.0 ) + tmp_8 * ( tmp_48 - 1.0 / 3.0 ) );
      real_t tmp_50 = 0.2844444444444445 * tmp_23;
      real_t tmp_51 = 0.7692346550528415 * tmp_0 + tmp_2;
      real_t tmp_52 = tmp_10 * tmp_51;
      real_t tmp_53 = tmp_12 * tmp_51;
      real_t tmp_54 = tmp_9 * ( 0.7692346550528415 * tmp_14 + tmp_15 );
      real_t tmp_55 = tmp_54 * tmp_6;
      real_t tmp_56 = tmp_18 * tmp_54;
      real_t tmp_57 = -tmp_52 - tmp_53 - tmp_55 - tmp_56 + 1;
      real_t tmp_58 = tmp_53 + tmp_55;
      real_t tmp_59 = tmp_52 + tmp_56;
      real_t tmp_60 = tmp_24 * ( tmp_5 * ( tmp_58 - 1.0 / 3.0 ) + tmp_8 * ( tmp_59 - 1.0 / 3.0 ) );
      real_t tmp_61 = 0.2393143352496831 * tmp_23;
      real_t tmp_62 = 0.95308992296933193 * tmp_0 + tmp_2;
      real_t tmp_63 = tmp_10 * tmp_62;
      real_t tmp_64 = tmp_12 * tmp_62;
      real_t tmp_65 = tmp_9 * ( 0.95308992296933193 * tmp_14 + tmp_15 );
      real_t tmp_66 = tmp_6 * tmp_65;
      real_t tmp_67 = tmp_18 * tmp_65;
      real_t tmp_68 = -tmp_63 - tmp_64 - tmp_66 - tmp_67 + 1;
      real_t tmp_69 = tmp_64 + tmp_66;
      real_t tmp_70 = tmp_63 + tmp_67;
      real_t tmp_71 = tmp_24 * ( tmp_5 * ( tmp_69 - 1.0 / 3.0 ) + tmp_8 * ( tmp_70 - 1.0 / 3.0 ) );
      real_t tmp_72 = 0.11846344252809471 * tmp_23;
      real_t a_0_0  = tmp_28 * ( -tmp_20 * tmp_22 + tmp_20 * tmp_27 ) + tmp_39 * ( -tmp_22 * tmp_35 + tmp_35 * tmp_38 ) +
                     tmp_50 * ( -tmp_22 * tmp_46 + tmp_46 * tmp_49 ) + tmp_61 * ( -tmp_22 * tmp_57 + tmp_57 * tmp_60 ) +
                     tmp_72 * ( -tmp_22 * tmp_68 + tmp_68 * tmp_71 );
      real_t a_1_0 = tmp_28 * ( -tmp_22 * tmp_25 + tmp_25 * tmp_27 ) + tmp_39 * ( -tmp_22 * tmp_36 + tmp_36 * tmp_38 ) +
                     tmp_50 * ( -tmp_22 * tmp_47 + tmp_47 * tmp_49 ) + tmp_61 * ( -tmp_22 * tmp_58 + tmp_58 * tmp_60 ) +
                     tmp_72 * ( -tmp_22 * tmp_69 + tmp_69 * tmp_71 );
      real_t a_2_0 = tmp_28 * ( -tmp_22 * tmp_26 + tmp_26 * tmp_27 ) + tmp_39 * ( -tmp_22 * tmp_37 + tmp_37 * tmp_38 ) +
                     tmp_50 * ( -tmp_22 * tmp_48 + tmp_48 * tmp_49 ) + tmp_61 * ( -tmp_22 * tmp_59 + tmp_59 * tmp_60 ) +
                     tmp_72 * ( -tmp_22 * tmp_70 + tmp_70 * tmp_71 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >& coordsElementInner,
                                          const std::vector< Point3D >& coordsElementOuter,
                                          const std::vector< Point3D >& coordsFacet,
                                          const Point3D&                oppositeVertexInnerElement,
                                          const Point3D&                oppositeVertexOuterElement,
                                          const Point3D&                outwardNormal,
                                          const DGBasisInfo&            trialBasis,
                                          const DGBasisInfo&            testBasis,
                                          int                           trialDegree,
                                          int                           testDegree,
                                          MatrixXr&                     elMat ) const override
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
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + 0.046910077030668018 * tmp_5;
      real_t tmp_7  = tmp_4 * ( tmp_2 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.046910077030668018 * tmp_11;
      real_t tmp_13 = tmp_4 * ( tmp_0 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = -p_affine_3_0;
      real_t tmp_19 = p_affine_4_0 + tmp_18;
      real_t tmp_20 = -p_affine_3_1;
      real_t tmp_21 = p_affine_5_1 + tmp_20;
      real_t tmp_22 = tmp_19 * tmp_21;
      real_t tmp_23 = p_affine_5_0 + tmp_18;
      real_t tmp_24 = 1.0 / ( tmp_22 - tmp_23 * ( p_affine_4_1 + tmp_20 ) );
      real_t tmp_25 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_26 = tmp_23 * tmp_24;
      real_t tmp_27 = tmp_24 * ( p_affine_3_0 - p_affine_5_0 );
      real_t tmp_28 = 0.5 * p_affine_10_0 * ( tmp_22 * tmp_24 + tmp_25 * tmp_26 ) +
                      0.5 * p_affine_10_1 * ( tmp_19 * tmp_26 + tmp_19 * tmp_27 );
      real_t tmp_29 = std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_30 = 1.0 / ( tmp_29 );
      real_t tmp_31 = tmp_20 + tmp_6;
      real_t tmp_32 = tmp_24 * ( tmp_12 + tmp_18 );
      real_t tmp_33 = tmp_19 * tmp_24;
      real_t tmp_34 = tmp_30 * ( tmp_19 * ( tmp_21 * tmp_32 + tmp_27 * tmp_31 - 1.0 / 3.0 ) +
                                 tmp_23 * ( tmp_25 * tmp_32 + tmp_31 * tmp_33 - 1.0 / 3.0 ) );
      real_t tmp_35 = 0.11846344252809471 * tmp_29;
      real_t tmp_36 = p_affine_6_1 + 0.23076534494715845 * tmp_5;
      real_t tmp_37 = tmp_4 * ( tmp_2 + tmp_36 );
      real_t tmp_38 = tmp_1 * tmp_37;
      real_t tmp_39 = tmp_37 * tmp_9;
      real_t tmp_40 = p_affine_6_0 + 0.23076534494715845 * tmp_11;
      real_t tmp_41 = tmp_4 * ( tmp_0 + tmp_40 );
      real_t tmp_42 = tmp_3 * tmp_41;
      real_t tmp_43 = tmp_15 * tmp_41;
      real_t tmp_44 = -tmp_38 - tmp_39 - tmp_42 - tmp_43 + 1;
      real_t tmp_45 = tmp_20 + tmp_36;
      real_t tmp_46 = tmp_24 * ( tmp_18 + tmp_40 );
      real_t tmp_47 = tmp_30 * ( tmp_19 * ( tmp_21 * tmp_46 + tmp_27 * tmp_45 - 1.0 / 3.0 ) +
                                 tmp_23 * ( tmp_25 * tmp_46 + tmp_33 * tmp_45 - 1.0 / 3.0 ) );
      real_t tmp_48 = 0.2393143352496831 * tmp_29;
      real_t tmp_49 = p_affine_6_1 + 0.5 * tmp_5;
      real_t tmp_50 = tmp_4 * ( tmp_2 + tmp_49 );
      real_t tmp_51 = tmp_1 * tmp_50;
      real_t tmp_52 = tmp_50 * tmp_9;
      real_t tmp_53 = p_affine_6_0 + 0.5 * tmp_11;
      real_t tmp_54 = tmp_4 * ( tmp_0 + tmp_53 );
      real_t tmp_55 = tmp_3 * tmp_54;
      real_t tmp_56 = tmp_15 * tmp_54;
      real_t tmp_57 = -tmp_51 - tmp_52 - tmp_55 - tmp_56 + 1;
      real_t tmp_58 = tmp_20 + tmp_49;
      real_t tmp_59 = tmp_24 * ( tmp_18 + tmp_53 );
      real_t tmp_60 = tmp_30 * ( tmp_19 * ( tmp_21 * tmp_59 + tmp_27 * tmp_58 - 1.0 / 3.0 ) +
                                 tmp_23 * ( tmp_25 * tmp_59 + tmp_33 * tmp_58 - 1.0 / 3.0 ) );
      real_t tmp_61 = 0.2844444444444445 * tmp_29;
      real_t tmp_62 = p_affine_6_1 + 0.7692346550528415 * tmp_5;
      real_t tmp_63 = tmp_4 * ( tmp_2 + tmp_62 );
      real_t tmp_64 = tmp_1 * tmp_63;
      real_t tmp_65 = tmp_63 * tmp_9;
      real_t tmp_66 = p_affine_6_0 + 0.7692346550528415 * tmp_11;
      real_t tmp_67 = tmp_4 * ( tmp_0 + tmp_66 );
      real_t tmp_68 = tmp_3 * tmp_67;
      real_t tmp_69 = tmp_15 * tmp_67;
      real_t tmp_70 = -tmp_64 - tmp_65 - tmp_68 - tmp_69 + 1;
      real_t tmp_71 = tmp_20 + tmp_62;
      real_t tmp_72 = tmp_24 * ( tmp_18 + tmp_66 );
      real_t tmp_73 = tmp_30 * ( tmp_19 * ( tmp_21 * tmp_72 + tmp_27 * tmp_71 - 1.0 / 3.0 ) +
                                 tmp_23 * ( tmp_25 * tmp_72 + tmp_33 * tmp_71 - 1.0 / 3.0 ) );
      real_t tmp_74 = 0.2393143352496831 * tmp_29;
      real_t tmp_75 = p_affine_6_1 + 0.95308992296933193 * tmp_5;
      real_t tmp_76 = tmp_4 * ( tmp_2 + tmp_75 );
      real_t tmp_77 = tmp_1 * tmp_76;
      real_t tmp_78 = tmp_76 * tmp_9;
      real_t tmp_79 = p_affine_6_0 + 0.95308992296933193 * tmp_11;
      real_t tmp_80 = tmp_4 * ( tmp_0 + tmp_79 );
      real_t tmp_81 = tmp_3 * tmp_80;
      real_t tmp_82 = tmp_15 * tmp_80;
      real_t tmp_83 = -tmp_77 - tmp_78 - tmp_81 - tmp_82 + 1;
      real_t tmp_84 = tmp_20 + tmp_75;
      real_t tmp_85 = tmp_24 * ( tmp_18 + tmp_79 );
      real_t tmp_86 = tmp_30 * ( tmp_19 * ( tmp_21 * tmp_85 + tmp_27 * tmp_84 - 1.0 / 3.0 ) +
                                 tmp_23 * ( tmp_25 * tmp_85 + tmp_33 * tmp_84 - 1.0 / 3.0 ) );
      real_t tmp_87 = 0.11846344252809471 * tmp_29;
      real_t tmp_88 = tmp_10 + tmp_14;
      real_t tmp_89 = tmp_39 + tmp_42;
      real_t tmp_90 = tmp_52 + tmp_55;
      real_t tmp_91 = tmp_65 + tmp_68;
      real_t tmp_92 = tmp_78 + tmp_81;
      real_t tmp_93 = tmp_16 + tmp_8;
      real_t tmp_94 = tmp_38 + tmp_43;
      real_t tmp_95 = tmp_51 + tmp_56;
      real_t tmp_96 = tmp_64 + tmp_69;
      real_t tmp_97 = tmp_77 + tmp_82;
      real_t a_0_0  = tmp_35 * ( -tmp_17 * tmp_28 - tmp_17 * tmp_34 ) + tmp_48 * ( -tmp_28 * tmp_44 - tmp_44 * tmp_47 ) +
                     tmp_61 * ( -tmp_28 * tmp_57 - tmp_57 * tmp_60 ) + tmp_74 * ( -tmp_28 * tmp_70 - tmp_70 * tmp_73 ) +
                     tmp_87 * ( -tmp_28 * tmp_83 - tmp_83 * tmp_86 );
      real_t a_1_0 = tmp_35 * ( -tmp_28 * tmp_88 - tmp_34 * tmp_88 ) + tmp_48 * ( -tmp_28 * tmp_89 - tmp_47 * tmp_89 ) +
                     tmp_61 * ( -tmp_28 * tmp_90 - tmp_60 * tmp_90 ) + tmp_74 * ( -tmp_28 * tmp_91 - tmp_73 * tmp_91 ) +
                     tmp_87 * ( -tmp_28 * tmp_92 - tmp_86 * tmp_92 );
      real_t a_2_0 = tmp_35 * ( -tmp_28 * tmp_93 - tmp_34 * tmp_93 ) + tmp_48 * ( -tmp_28 * tmp_94 - tmp_47 * tmp_94 ) +
                     tmp_61 * ( -tmp_28 * tmp_95 - tmp_60 * tmp_95 ) + tmp_74 * ( -tmp_28 * tmp_96 - tmp_73 * tmp_96 ) +
                     tmp_87 * ( -tmp_28 * tmp_97 - tmp_86 * tmp_97 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                                   const std::vector< Point3D >& coordsFacet,
                                                   const Point3D&                oppositeVertex,
                                                   const Point3D&                outwardNormal,
                                                   const DGBasisInfo&            trialBasis,
                                                   const DGBasisInfo&            testBasis,
                                                   int                           trialDegree,
                                                   int                           testDegree,
                                                   MatrixXr&                     elMat ) const override
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
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = p_affine_2_0 + tmp_0;
      real_t tmp_5  = 1.0 / ( tmp_1 * tmp_3 - tmp_4 * ( p_affine_1_1 + tmp_2 ) );
      real_t tmp_6  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_7  = p_affine_6_1 + tmp_2;
      real_t tmp_8  = tmp_5 * ( 0.046910077030668018 * tmp_6 + tmp_7 );
      real_t tmp_9  = tmp_1 * tmp_8;
      real_t tmp_10 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_11 = tmp_10 * tmp_8;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_0;
      real_t tmp_14 = tmp_5 * ( 0.046910077030668018 * tmp_12 + tmp_13 );
      real_t tmp_15 = tmp_14 * tmp_3;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_14 * tmp_16;
      real_t tmp_18 = tmp_11 + tmp_15;
      real_t tmp_19 = tmp_17 + tmp_9;
      real_t tmp_20 = 0.11846344252809471 * tmp_1 * ( tmp_18 - 1.0 / 3.0 ) + 0.11846344252809471 * tmp_4 * ( tmp_19 - 1.0 / 3.0 );
      real_t tmp_21 = tmp_5 * ( 0.23076534494715845 * tmp_6 + tmp_7 );
      real_t tmp_22 = tmp_1 * tmp_21;
      real_t tmp_23 = tmp_10 * tmp_21;
      real_t tmp_24 = tmp_5 * ( 0.23076534494715845 * tmp_12 + tmp_13 );
      real_t tmp_25 = tmp_24 * tmp_3;
      real_t tmp_26 = tmp_16 * tmp_24;
      real_t tmp_27 = tmp_23 + tmp_25;
      real_t tmp_28 = tmp_22 + tmp_26;
      real_t tmp_29 = 0.2393143352496831 * tmp_1 * ( tmp_27 - 1.0 / 3.0 ) + 0.2393143352496831 * tmp_4 * ( tmp_28 - 1.0 / 3.0 );
      real_t tmp_30 = tmp_5 * ( 0.5 * tmp_6 + tmp_7 );
      real_t tmp_31 = tmp_1 * tmp_30;
      real_t tmp_32 = tmp_10 * tmp_30;
      real_t tmp_33 = tmp_5 * ( 0.5 * tmp_12 + tmp_13 );
      real_t tmp_34 = tmp_3 * tmp_33;
      real_t tmp_35 = tmp_16 * tmp_33;
      real_t tmp_36 = tmp_32 + tmp_34;
      real_t tmp_37 = tmp_31 + tmp_35;
      real_t tmp_38 = 0.2844444444444445 * tmp_1 * ( tmp_36 - 1.0 / 3.0 ) + 0.2844444444444445 * tmp_4 * ( tmp_37 - 1.0 / 3.0 );
      real_t tmp_39 = tmp_5 * ( 0.7692346550528415 * tmp_6 + tmp_7 );
      real_t tmp_40 = tmp_1 * tmp_39;
      real_t tmp_41 = tmp_10 * tmp_39;
      real_t tmp_42 = tmp_5 * ( 0.7692346550528415 * tmp_12 + tmp_13 );
      real_t tmp_43 = tmp_3 * tmp_42;
      real_t tmp_44 = tmp_16 * tmp_42;
      real_t tmp_45 = tmp_41 + tmp_43;
      real_t tmp_46 = tmp_40 + tmp_44;
      real_t tmp_47 = 0.2393143352496831 * tmp_1 * ( tmp_45 - 1.0 / 3.0 ) + 0.2393143352496831 * tmp_4 * ( tmp_46 - 1.0 / 3.0 );
      real_t tmp_48 = tmp_5 * ( 0.95308992296933193 * tmp_6 + tmp_7 );
      real_t tmp_49 = tmp_1 * tmp_48;
      real_t tmp_50 = tmp_10 * tmp_48;
      real_t tmp_51 = tmp_5 * ( 0.95308992296933193 * tmp_12 + tmp_13 );
      real_t tmp_52 = tmp_3 * tmp_51;
      real_t tmp_53 = tmp_16 * tmp_51;
      real_t tmp_54 = tmp_50 + tmp_52;
      real_t tmp_55 = tmp_49 + tmp_53;
      real_t tmp_56 = 0.11846344252809471 * tmp_1 * ( tmp_54 - 1.0 / 3.0 ) + 0.11846344252809471 * tmp_4 * ( tmp_55 - 1.0 / 3.0 );
      real_t a_0_0  = tmp_20 * ( -tmp_11 - tmp_15 - tmp_17 - tmp_9 + 1 ) + tmp_29 * ( -tmp_22 - tmp_23 - tmp_25 - tmp_26 + 1 ) +
                     tmp_38 * ( -tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1 ) + tmp_47 * ( -tmp_40 - tmp_41 - tmp_43 - tmp_44 + 1 ) +
                     tmp_56 * ( -tmp_49 - tmp_50 - tmp_52 - tmp_53 + 1 );
      real_t a_1_0  = tmp_18 * tmp_20 + tmp_27 * tmp_29 + tmp_36 * tmp_38 + tmp_45 * tmp_47 + tmp_54 * tmp_56;
      real_t a_2_0  = tmp_19 * tmp_20 + tmp_28 * tmp_29 + tmp_37 * tmp_38 + tmp_46 * tmp_47 + tmp_55 * tmp_56;
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
   }

   void integrateRHSDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                         const std::vector< Point3D >& coordsFacet,
                                         const Point3D&                oppositeVertex,
                                         const Point3D&                outwardNormal,
                                         const DGBasisInfo&            basis,
                                         int                           degree,
                                         MatrixXr&                     elMat ) const override
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
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                         const std::vector< Point3D >& coordsFacet,
                                         const Point3D&                oppositeVertex,
                                         const Point3D&                outwardNormal,
                                         const DGBasisInfo&            basis,
                                         int                           degree,
                                         MatrixXr&                     elMat ) const override
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
   void integrateVolume3D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_0;
      real_t tmp_6  = p_affine_1_1 + tmp_2;
      real_t tmp_7  = tmp_5 * tmp_6;
      real_t tmp_8  = -p_affine_0_2;
      real_t tmp_9  = p_affine_3_2 + tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_8;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11 * tmp_5;
      real_t tmp_13 = p_affine_2_2 + tmp_8;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = tmp_14 * tmp_6;
      real_t tmp_16 = tmp_1 * tmp_11;
      real_t tmp_17 = tmp_14 * tmp_3;
      real_t tmp_18 =
          1.0 / ( tmp_10 * tmp_12 - tmp_10 * tmp_17 + tmp_13 * tmp_15 - tmp_13 * tmp_16 + tmp_4 * tmp_9 - tmp_7 * tmp_9 );
      real_t tmp_19 = tmp_18 * ( tmp_4 - tmp_7 );
      real_t tmp_20 = tmp_18 * ( tmp_15 - tmp_16 );
      real_t tmp_21 = tmp_18 * ( tmp_12 - tmp_17 );
      real_t tmp_22 = tmp_1 * tmp_21 + tmp_14 * tmp_19 + tmp_20 * tmp_5;
      real_t tmp_23 = tmp_18 * ( -tmp_1 * tmp_13 + tmp_10 * tmp_5 );
      real_t tmp_24 = tmp_18 * ( tmp_1 * tmp_9 - tmp_10 * tmp_14 );
      real_t tmp_25 = tmp_18 * ( tmp_13 * tmp_14 - tmp_5 * tmp_9 );
      real_t tmp_26 = tmp_1 * tmp_25 + tmp_14 * tmp_23 + tmp_24 * tmp_5;
      real_t tmp_27 = tmp_18 * ( -tmp_10 * tmp_3 + tmp_13 * tmp_6 );
      real_t tmp_28 = tmp_18 * ( tmp_10 * tmp_11 - tmp_6 * tmp_9 );
      real_t tmp_29 = tmp_18 * ( -tmp_11 * tmp_13 + tmp_3 * tmp_9 );
      real_t tmp_30 = tmp_1 * tmp_29 + tmp_14 * tmp_27 + tmp_28 * tmp_5;
      real_t tmp_31 = p_affine_0_0 * p_affine_1_1;
      real_t tmp_32 = p_affine_0_0 * p_affine_1_2;
      real_t tmp_33 = p_affine_2_1 * p_affine_3_2;
      real_t tmp_34 = p_affine_0_1 * p_affine_1_0;
      real_t tmp_35 = p_affine_0_1 * p_affine_1_2;
      real_t tmp_36 = p_affine_2_2 * p_affine_3_0;
      real_t tmp_37 = p_affine_0_2 * p_affine_1_0;
      real_t tmp_38 = p_affine_0_2 * p_affine_1_1;
      real_t tmp_39 = p_affine_2_0 * p_affine_3_1;
      real_t tmp_40 = p_affine_2_2 * p_affine_3_1;
      real_t tmp_41 = p_affine_2_0 * p_affine_3_2;
      real_t tmp_42 = p_affine_2_1 * p_affine_3_0;
      real_t tmp_43 = std::abs( p_affine_0_0 * tmp_33 - p_affine_0_0 * tmp_40 + p_affine_0_1 * tmp_36 - p_affine_0_1 * tmp_41 +
                                p_affine_0_2 * tmp_39 - p_affine_0_2 * tmp_42 - p_affine_1_0 * tmp_33 + p_affine_1_0 * tmp_40 -
                                p_affine_1_1 * tmp_36 + p_affine_1_1 * tmp_41 - p_affine_1_2 * tmp_39 + p_affine_1_2 * tmp_42 +
                                p_affine_2_0 * tmp_35 - p_affine_2_0 * tmp_38 - p_affine_2_1 * tmp_32 + p_affine_2_1 * tmp_37 +
                                p_affine_2_2 * tmp_31 - p_affine_2_2 * tmp_34 - p_affine_3_0 * tmp_35 + p_affine_3_0 * tmp_38 +
                                p_affine_3_1 * tmp_32 - p_affine_3_1 * tmp_37 - p_affine_3_2 * tmp_31 + p_affine_3_2 * tmp_34 );
      real_t tmp_44 = tmp_43 * ( tmp_22 * ( -tmp_19 - tmp_20 - tmp_21 ) + tmp_26 * ( -tmp_23 - tmp_24 - tmp_25 ) +
                                 tmp_30 * ( -tmp_27 - tmp_28 - tmp_29 ) );
      real_t tmp_45 = tmp_43 * ( tmp_21 * tmp_22 + tmp_25 * tmp_26 + tmp_29 * tmp_30 );
      real_t tmp_46 = tmp_43 * ( tmp_20 * tmp_22 + tmp_24 * tmp_26 + tmp_28 * tmp_30 );
      real_t tmp_47 = tmp_43 * ( tmp_19 * tmp_22 + tmp_23 * tmp_26 + tmp_27 * tmp_30 );
      real_t a_0_0  = 0.1666666666666668 * tmp_44;
      real_t a_1_0  = 0.1666666666666668 * tmp_45;
      real_t a_2_0  = 0.1666666666666668 * tmp_46;
      real_t a_3_0  = 0.1666666666666668 * tmp_47;
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
      elMat( 3, 0 ) = a_3_0;
   }

   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                               const std::vector< Point3D >& coordsFacet,
                               const Point3D&,
                               const Point3D&     outwardNormal,
                               const DGBasisInfo& trialBasis,
                               const DGBasisInfo& testBasis,
                               int                trialDegree,
                               int                testDegree,
                               MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = p_affine_2_0 + tmp_0;
      real_t tmp_5  = p_affine_1_1 + tmp_2;
      real_t tmp_6  = tmp_1 * tmp_3 - tmp_4 * tmp_5;
      real_t tmp_7  = -p_affine_0_2;
      real_t tmp_8  = p_affine_3_2 + tmp_7;
      real_t tmp_9  = tmp_3 * tmp_8;
      real_t tmp_10 = p_affine_3_1 + tmp_2;
      real_t tmp_11 = p_affine_1_2 + tmp_7;
      real_t tmp_12 = tmp_10 * tmp_11;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_14 * tmp_5;
      real_t tmp_16 = tmp_10 * tmp_14;
      real_t tmp_17 = tmp_5 * tmp_8;
      real_t tmp_18 = tmp_11 * tmp_3;
      real_t tmp_19 =
          1.0 / ( -tmp_1 * tmp_16 + tmp_1 * tmp_9 + tmp_12 * tmp_4 + tmp_13 * tmp_15 - tmp_13 * tmp_18 - tmp_17 * tmp_4 );
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19 * ( 0.031405749086161582 * tmp_21 + 0.93718850182767688 * tmp_22 + tmp_23 );
      real_t tmp_25 = tmp_24 * tmp_6;
      real_t tmp_26 = -tmp_1 * tmp_10 + tmp_13 * tmp_5;
      real_t tmp_27 = tmp_24 * tmp_26;
      real_t tmp_28 = -tmp_1 * tmp_14 + tmp_11 * tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19 * ( 0.031405749086161582 * tmp_30 + 0.93718850182767688 * tmp_31 + tmp_32 );
      real_t tmp_34 = tmp_28 * tmp_33;
      real_t tmp_35 = tmp_1 * tmp_8 - tmp_11 * tmp_13;
      real_t tmp_36 = tmp_33 * tmp_35;
      real_t tmp_37 = tmp_10 * tmp_4 - tmp_13 * tmp_3;
      real_t tmp_38 = tmp_24 * tmp_37;
      real_t tmp_39 = tmp_13 * tmp_14 - tmp_4 * tmp_8;
      real_t tmp_40 = tmp_33 * tmp_39;
      real_t tmp_41 = tmp_15 - tmp_18;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19 * ( 0.031405749086161582 * tmp_43 + 0.93718850182767688 * tmp_44 + tmp_45 );
      real_t tmp_47 = tmp_41 * tmp_46;
      real_t tmp_48 = tmp_12 - tmp_17;
      real_t tmp_49 = tmp_46 * tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46 * tmp_50;
      real_t tmp_52 = -tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1;
      real_t tmp_53 = tmp_1 * tmp_19;
      real_t tmp_54 = tmp_19 * tmp_4;
      real_t tmp_55 = tmp_13 * tmp_19;
      real_t tmp_56 = 0.5 * p_affine_13_0 * ( tmp_41 * tmp_55 + tmp_48 * tmp_54 + tmp_50 * tmp_53 ) +
                      0.5 * p_affine_13_1 * ( tmp_28 * tmp_55 + tmp_35 * tmp_54 + tmp_39 * tmp_53 ) +
                      0.5 * p_affine_13_2 * ( tmp_26 * tmp_54 + tmp_37 * tmp_53 + tmp_55 * tmp_6 );
      real_t tmp_57 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_58 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_59 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_60 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_61 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_62 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_63 = ( std::abs( tmp_22 * tmp_60 - tmp_31 * tmp_62 ) * std::abs( tmp_22 * tmp_60 - tmp_31 * tmp_62 ) ) +
                      ( std::abs( tmp_22 * tmp_61 - tmp_44 * tmp_62 ) * std::abs( tmp_22 * tmp_61 - tmp_44 * tmp_62 ) ) +
                      ( std::abs( tmp_31 * tmp_61 - tmp_44 * tmp_60 ) * std::abs( tmp_31 * tmp_61 - tmp_44 * tmp_60 ) );
      real_t tmp_64 = 1.0 * std::pow( tmp_63, -0.25 );
      real_t tmp_65 =
          tmp_64 * ( tmp_1 * ( tmp_57 - 1.0 / 4.0 ) + tmp_13 * ( tmp_59 - 1.0 / 4.0 ) + tmp_4 * ( tmp_58 - 1.0 / 4.0 ) );
      real_t tmp_66 = 1.0 * std::pow( tmp_63, 1.0 / 2.0 );
      real_t tmp_67 = 0.0068572537431980923 * tmp_66;
      real_t tmp_68 = tmp_19 * ( 0.19601935860219369 * tmp_21 + 0.60796128279561268 * tmp_22 + tmp_23 );
      real_t tmp_69 = tmp_6 * tmp_68;
      real_t tmp_70 = tmp_26 * tmp_68;
      real_t tmp_71 = tmp_19 * ( 0.19601935860219369 * tmp_30 + 0.60796128279561268 * tmp_31 + tmp_32 );
      real_t tmp_72 = tmp_28 * tmp_71;
      real_t tmp_73 = tmp_35 * tmp_71;
      real_t tmp_74 = tmp_37 * tmp_68;
      real_t tmp_75 = tmp_39 * tmp_71;
      real_t tmp_76 = tmp_19 * ( 0.19601935860219369 * tmp_43 + 0.60796128279561268 * tmp_44 + tmp_45 );
      real_t tmp_77 = tmp_41 * tmp_76;
      real_t tmp_78 = tmp_48 * tmp_76;
      real_t tmp_79 = tmp_50 * tmp_76;
      real_t tmp_80 = -tmp_69 - tmp_70 - tmp_72 - tmp_73 - tmp_74 - tmp_75 - tmp_77 - tmp_78 - tmp_79 + 1;
      real_t tmp_81 = tmp_74 + tmp_75 + tmp_79;
      real_t tmp_82 = tmp_70 + tmp_73 + tmp_78;
      real_t tmp_83 = tmp_69 + tmp_72 + tmp_77;
      real_t tmp_84 =
          tmp_64 * ( tmp_1 * ( tmp_81 - 1.0 / 4.0 ) + tmp_13 * ( tmp_83 - 1.0 / 4.0 ) + tmp_4 * ( tmp_82 - 1.0 / 4.0 ) );
      real_t tmp_85  = 0.037198804536718075 * tmp_66;
      real_t tmp_86  = tmp_19 * ( 0.37605877282253791 * tmp_21 + 0.039308471900058539 * tmp_22 + tmp_23 );
      real_t tmp_87  = tmp_6 * tmp_86;
      real_t tmp_88  = tmp_26 * tmp_86;
      real_t tmp_89  = tmp_19 * ( 0.37605877282253791 * tmp_30 + 0.039308471900058539 * tmp_31 + tmp_32 );
      real_t tmp_90  = tmp_28 * tmp_89;
      real_t tmp_91  = tmp_35 * tmp_89;
      real_t tmp_92  = tmp_37 * tmp_86;
      real_t tmp_93  = tmp_39 * tmp_89;
      real_t tmp_94  = tmp_19 * ( 0.37605877282253791 * tmp_43 + 0.039308471900058539 * tmp_44 + tmp_45 );
      real_t tmp_95  = tmp_41 * tmp_94;
      real_t tmp_96  = tmp_48 * tmp_94;
      real_t tmp_97  = tmp_50 * tmp_94;
      real_t tmp_98  = -tmp_87 - tmp_88 - tmp_90 - tmp_91 - tmp_92 - tmp_93 - tmp_95 - tmp_96 - tmp_97 + 1;
      real_t tmp_99  = tmp_92 + tmp_93 + tmp_97;
      real_t tmp_100 = tmp_88 + tmp_91 + tmp_96;
      real_t tmp_101 = tmp_87 + tmp_90 + tmp_95;
      real_t tmp_102 =
          tmp_64 * ( tmp_1 * ( tmp_99 - 1.0 / 4.0 ) + tmp_13 * ( tmp_101 - 1.0 / 4.0 ) + tmp_4 * ( tmp_100 - 1.0 / 4.0 ) );
      real_t tmp_103 = 0.020848748529055869 * tmp_66;
      real_t tmp_104 = tmp_19 * ( 0.78764240869137092 * tmp_21 + 0.1711304259088916 * tmp_22 + tmp_23 );
      real_t tmp_105 = tmp_104 * tmp_6;
      real_t tmp_106 = tmp_104 * tmp_26;
      real_t tmp_107 = tmp_19 * ( 0.78764240869137092 * tmp_30 + 0.1711304259088916 * tmp_31 + tmp_32 );
      real_t tmp_108 = tmp_107 * tmp_28;
      real_t tmp_109 = tmp_107 * tmp_35;
      real_t tmp_110 = tmp_104 * tmp_37;
      real_t tmp_111 = tmp_107 * tmp_39;
      real_t tmp_112 = tmp_19 * ( 0.78764240869137092 * tmp_43 + 0.1711304259088916 * tmp_44 + tmp_45 );
      real_t tmp_113 = tmp_112 * tmp_41;
      real_t tmp_114 = tmp_112 * tmp_48;
      real_t tmp_115 = tmp_112 * tmp_50;
      real_t tmp_116 = -tmp_105 - tmp_106 - tmp_108 - tmp_109 - tmp_110 - tmp_111 - tmp_113 - tmp_114 - tmp_115 + 1;
      real_t tmp_117 = tmp_110 + tmp_111 + tmp_115;
      real_t tmp_118 = tmp_106 + tmp_109 + tmp_114;
      real_t tmp_119 = tmp_105 + tmp_108 + tmp_113;
      real_t tmp_120 =
          tmp_64 * ( tmp_1 * ( tmp_117 - 1.0 / 4.0 ) + tmp_13 * ( tmp_119 - 1.0 / 4.0 ) + tmp_4 * ( tmp_118 - 1.0 / 4.0 ) );
      real_t tmp_121 = 0.019202922745021479 * tmp_66;
      real_t tmp_122 = tmp_19 * ( 0.58463275527740355 * tmp_21 + 0.37605877282253791 * tmp_22 + tmp_23 );
      real_t tmp_123 = tmp_122 * tmp_6;
      real_t tmp_124 = tmp_122 * tmp_26;
      real_t tmp_125 = tmp_19 * ( 0.58463275527740355 * tmp_30 + 0.37605877282253791 * tmp_31 + tmp_32 );
      real_t tmp_126 = tmp_125 * tmp_28;
      real_t tmp_127 = tmp_125 * tmp_35;
      real_t tmp_128 = tmp_122 * tmp_37;
      real_t tmp_129 = tmp_125 * tmp_39;
      real_t tmp_130 = tmp_19 * ( 0.58463275527740355 * tmp_43 + 0.37605877282253791 * tmp_44 + tmp_45 );
      real_t tmp_131 = tmp_130 * tmp_41;
      real_t tmp_132 = tmp_130 * tmp_48;
      real_t tmp_133 = tmp_130 * tmp_50;
      real_t tmp_134 = -tmp_123 - tmp_124 - tmp_126 - tmp_127 - tmp_128 - tmp_129 - tmp_131 - tmp_132 - tmp_133 + 1;
      real_t tmp_135 = tmp_128 + tmp_129 + tmp_133;
      real_t tmp_136 = tmp_124 + tmp_127 + tmp_132;
      real_t tmp_137 = tmp_123 + tmp_126 + tmp_131;
      real_t tmp_138 =
          tmp_64 * ( tmp_1 * ( tmp_135 - 1.0 / 4.0 ) + tmp_13 * ( tmp_137 - 1.0 / 4.0 ) + tmp_4 * ( tmp_136 - 1.0 / 4.0 ) );
      real_t tmp_139 = 0.020848748529055869 * tmp_66;
      real_t tmp_140 = tmp_19 * ( 0.041227165399737475 * tmp_21 + 0.78764240869137092 * tmp_22 + tmp_23 );
      real_t tmp_141 = tmp_140 * tmp_6;
      real_t tmp_142 = tmp_140 * tmp_26;
      real_t tmp_143 = tmp_19 * ( 0.041227165399737475 * tmp_30 + 0.78764240869137092 * tmp_31 + tmp_32 );
      real_t tmp_144 = tmp_143 * tmp_28;
      real_t tmp_145 = tmp_143 * tmp_35;
      real_t tmp_146 = tmp_140 * tmp_37;
      real_t tmp_147 = tmp_143 * tmp_39;
      real_t tmp_148 = tmp_19 * ( 0.041227165399737475 * tmp_43 + 0.78764240869137092 * tmp_44 + tmp_45 );
      real_t tmp_149 = tmp_148 * tmp_41;
      real_t tmp_150 = tmp_148 * tmp_48;
      real_t tmp_151 = tmp_148 * tmp_50;
      real_t tmp_152 = -tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1;
      real_t tmp_153 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_154 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_155 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_156 =
          tmp_64 * ( tmp_1 * ( tmp_153 - 1.0 / 4.0 ) + tmp_13 * ( tmp_155 - 1.0 / 4.0 ) + tmp_4 * ( tmp_154 - 1.0 / 4.0 ) );
      real_t tmp_157 = 0.019202922745021479 * tmp_66;
      real_t tmp_158 = tmp_19 * ( 0.039308471900058539 * tmp_21 + 0.58463275527740355 * tmp_22 + tmp_23 );
      real_t tmp_159 = tmp_158 * tmp_6;
      real_t tmp_160 = tmp_158 * tmp_26;
      real_t tmp_161 = tmp_19 * ( 0.039308471900058539 * tmp_30 + 0.58463275527740355 * tmp_31 + tmp_32 );
      real_t tmp_162 = tmp_161 * tmp_28;
      real_t tmp_163 = tmp_161 * tmp_35;
      real_t tmp_164 = tmp_158 * tmp_37;
      real_t tmp_165 = tmp_161 * tmp_39;
      real_t tmp_166 = tmp_19 * ( 0.039308471900058539 * tmp_43 + 0.58463275527740355 * tmp_44 + tmp_45 );
      real_t tmp_167 = tmp_166 * tmp_41;
      real_t tmp_168 = tmp_166 * tmp_48;
      real_t tmp_169 = tmp_166 * tmp_50;
      real_t tmp_170 = -tmp_159 - tmp_160 - tmp_162 - tmp_163 - tmp_164 - tmp_165 - tmp_167 - tmp_168 - tmp_169 + 1;
      real_t tmp_171 = tmp_164 + tmp_165 + tmp_169;
      real_t tmp_172 = tmp_160 + tmp_163 + tmp_168;
      real_t tmp_173 = tmp_159 + tmp_162 + tmp_167;
      real_t tmp_174 =
          tmp_64 * ( tmp_1 * ( tmp_171 - 1.0 / 4.0 ) + tmp_13 * ( tmp_173 - 1.0 / 4.0 ) + tmp_4 * ( tmp_172 - 1.0 / 4.0 ) );
      real_t tmp_175 = 0.020848748529055869 * tmp_66;
      real_t tmp_176 = tmp_19 * ( 0.78764240869137092 * tmp_21 + 0.041227165399737475 * tmp_22 + tmp_23 );
      real_t tmp_177 = tmp_176 * tmp_6;
      real_t tmp_178 = tmp_176 * tmp_26;
      real_t tmp_179 = tmp_19 * ( 0.78764240869137092 * tmp_30 + 0.041227165399737475 * tmp_31 + tmp_32 );
      real_t tmp_180 = tmp_179 * tmp_28;
      real_t tmp_181 = tmp_179 * tmp_35;
      real_t tmp_182 = tmp_176 * tmp_37;
      real_t tmp_183 = tmp_179 * tmp_39;
      real_t tmp_184 = tmp_19 * ( 0.78764240869137092 * tmp_43 + 0.041227165399737475 * tmp_44 + tmp_45 );
      real_t tmp_185 = tmp_184 * tmp_41;
      real_t tmp_186 = tmp_184 * tmp_48;
      real_t tmp_187 = tmp_184 * tmp_50;
      real_t tmp_188 = -tmp_177 - tmp_178 - tmp_180 - tmp_181 - tmp_182 - tmp_183 - tmp_185 - tmp_186 - tmp_187 + 1;
      real_t tmp_189 = tmp_182 + tmp_183 + tmp_187;
      real_t tmp_190 = tmp_178 + tmp_181 + tmp_186;
      real_t tmp_191 = tmp_177 + tmp_180 + tmp_185;
      real_t tmp_192 =
          tmp_64 * ( tmp_1 * ( tmp_189 - 1.0 / 4.0 ) + tmp_13 * ( tmp_191 - 1.0 / 4.0 ) + tmp_4 * ( tmp_190 - 1.0 / 4.0 ) );
      real_t tmp_193 = 0.019202922745021479 * tmp_66;
      real_t tmp_194 = tmp_19 * ( 0.58463275527740355 * tmp_21 + 0.039308471900058539 * tmp_22 + tmp_23 );
      real_t tmp_195 = tmp_194 * tmp_6;
      real_t tmp_196 = tmp_194 * tmp_26;
      real_t tmp_197 = tmp_19 * ( 0.58463275527740355 * tmp_30 + 0.039308471900058539 * tmp_31 + tmp_32 );
      real_t tmp_198 = tmp_197 * tmp_28;
      real_t tmp_199 = tmp_197 * tmp_35;
      real_t tmp_200 = tmp_194 * tmp_37;
      real_t tmp_201 = tmp_197 * tmp_39;
      real_t tmp_202 = tmp_19 * ( 0.58463275527740355 * tmp_43 + 0.039308471900058539 * tmp_44 + tmp_45 );
      real_t tmp_203 = tmp_202 * tmp_41;
      real_t tmp_204 = tmp_202 * tmp_48;
      real_t tmp_205 = tmp_202 * tmp_50;
      real_t tmp_206 = -tmp_195 - tmp_196 - tmp_198 - tmp_199 - tmp_200 - tmp_201 - tmp_203 - tmp_204 - tmp_205 + 1;
      real_t tmp_207 = tmp_200 + tmp_201 + tmp_205;
      real_t tmp_208 = tmp_196 + tmp_199 + tmp_204;
      real_t tmp_209 = tmp_195 + tmp_198 + tmp_203;
      real_t tmp_210 =
          tmp_64 * ( tmp_1 * ( tmp_207 - 1.0 / 4.0 ) + tmp_13 * ( tmp_209 - 1.0 / 4.0 ) + tmp_4 * ( tmp_208 - 1.0 / 4.0 ) );
      real_t tmp_211 = 0.020848748529055869 * tmp_66;
      real_t tmp_212 = tmp_19 * ( 0.1711304259088916 * tmp_21 + 0.78764240869137092 * tmp_22 + tmp_23 );
      real_t tmp_213 = tmp_212 * tmp_6;
      real_t tmp_214 = tmp_212 * tmp_26;
      real_t tmp_215 = tmp_19 * ( 0.1711304259088916 * tmp_30 + 0.78764240869137092 * tmp_31 + tmp_32 );
      real_t tmp_216 = tmp_215 * tmp_28;
      real_t tmp_217 = tmp_215 * tmp_35;
      real_t tmp_218 = tmp_212 * tmp_37;
      real_t tmp_219 = tmp_215 * tmp_39;
      real_t tmp_220 = tmp_19 * ( 0.1711304259088916 * tmp_43 + 0.78764240869137092 * tmp_44 + tmp_45 );
      real_t tmp_221 = tmp_220 * tmp_41;
      real_t tmp_222 = tmp_220 * tmp_48;
      real_t tmp_223 = tmp_220 * tmp_50;
      real_t tmp_224 = -tmp_213 - tmp_214 - tmp_216 - tmp_217 - tmp_218 - tmp_219 - tmp_221 - tmp_222 - tmp_223 + 1;
      real_t tmp_225 = tmp_218 + tmp_219 + tmp_223;
      real_t tmp_226 = tmp_214 + tmp_217 + tmp_222;
      real_t tmp_227 = tmp_213 + tmp_216 + tmp_221;
      real_t tmp_228 =
          tmp_64 * ( tmp_1 * ( tmp_225 - 1.0 / 4.0 ) + tmp_13 * ( tmp_227 - 1.0 / 4.0 ) + tmp_4 * ( tmp_226 - 1.0 / 4.0 ) );
      real_t tmp_229 = 0.019202922745021479 * tmp_66;
      real_t tmp_230 = tmp_19 * ( 0.37605877282253791 * tmp_21 + 0.58463275527740355 * tmp_22 + tmp_23 );
      real_t tmp_231 = tmp_230 * tmp_6;
      real_t tmp_232 = tmp_230 * tmp_26;
      real_t tmp_233 = tmp_19 * ( 0.37605877282253791 * tmp_30 + 0.58463275527740355 * tmp_31 + tmp_32 );
      real_t tmp_234 = tmp_233 * tmp_28;
      real_t tmp_235 = tmp_233 * tmp_35;
      real_t tmp_236 = tmp_230 * tmp_37;
      real_t tmp_237 = tmp_233 * tmp_39;
      real_t tmp_238 = tmp_19 * ( 0.37605877282253791 * tmp_43 + 0.58463275527740355 * tmp_44 + tmp_45 );
      real_t tmp_239 = tmp_238 * tmp_41;
      real_t tmp_240 = tmp_238 * tmp_48;
      real_t tmp_241 = tmp_238 * tmp_50;
      real_t tmp_242 = -tmp_231 - tmp_232 - tmp_234 - tmp_235 - tmp_236 - tmp_237 - tmp_239 - tmp_240 - tmp_241 + 1;
      real_t tmp_243 = tmp_236 + tmp_237 + tmp_241;
      real_t tmp_244 = tmp_232 + tmp_235 + tmp_240;
      real_t tmp_245 = tmp_231 + tmp_234 + tmp_239;
      real_t tmp_246 =
          tmp_64 * ( tmp_1 * ( tmp_243 - 1.0 / 4.0 ) + tmp_13 * ( tmp_245 - 1.0 / 4.0 ) + tmp_4 * ( tmp_244 - 1.0 / 4.0 ) );
      real_t tmp_247 = 0.020848748529055869 * tmp_66;
      real_t tmp_248 = tmp_19 * ( 0.041227165399737475 * tmp_21 + 0.1711304259088916 * tmp_22 + tmp_23 );
      real_t tmp_249 = tmp_248 * tmp_6;
      real_t tmp_250 = tmp_248 * tmp_26;
      real_t tmp_251 = tmp_19 * ( 0.041227165399737475 * tmp_30 + 0.1711304259088916 * tmp_31 + tmp_32 );
      real_t tmp_252 = tmp_251 * tmp_28;
      real_t tmp_253 = tmp_251 * tmp_35;
      real_t tmp_254 = tmp_248 * tmp_37;
      real_t tmp_255 = tmp_251 * tmp_39;
      real_t tmp_256 = tmp_19 * ( 0.041227165399737475 * tmp_43 + 0.1711304259088916 * tmp_44 + tmp_45 );
      real_t tmp_257 = tmp_256 * tmp_41;
      real_t tmp_258 = tmp_256 * tmp_48;
      real_t tmp_259 = tmp_256 * tmp_50;
      real_t tmp_260 = -tmp_249 - tmp_250 - tmp_252 - tmp_253 - tmp_254 - tmp_255 - tmp_257 - tmp_258 - tmp_259 + 1;
      real_t tmp_261 = tmp_254 + tmp_255 + tmp_259;
      real_t tmp_262 = tmp_250 + tmp_253 + tmp_258;
      real_t tmp_263 = tmp_249 + tmp_252 + tmp_257;
      real_t tmp_264 =
          tmp_64 * ( tmp_1 * ( tmp_261 - 1.0 / 4.0 ) + tmp_13 * ( tmp_263 - 1.0 / 4.0 ) + tmp_4 * ( tmp_262 - 1.0 / 4.0 ) );
      real_t tmp_265 = 0.019202922745021479 * tmp_66;
      real_t tmp_266 = tmp_19 * ( 0.40446199974765351 * tmp_21 + 0.19107600050469298 * tmp_22 + tmp_23 );
      real_t tmp_267 = tmp_266 * tmp_6;
      real_t tmp_268 = tmp_26 * tmp_266;
      real_t tmp_269 = tmp_19 * ( 0.40446199974765351 * tmp_30 + 0.19107600050469298 * tmp_31 + tmp_32 );
      real_t tmp_270 = tmp_269 * tmp_28;
      real_t tmp_271 = tmp_269 * tmp_35;
      real_t tmp_272 = tmp_266 * tmp_37;
      real_t tmp_273 = tmp_269 * tmp_39;
      real_t tmp_274 = tmp_19 * ( 0.40446199974765351 * tmp_43 + 0.19107600050469298 * tmp_44 + tmp_45 );
      real_t tmp_275 = tmp_274 * tmp_41;
      real_t tmp_276 = tmp_274 * tmp_48;
      real_t tmp_277 = tmp_274 * tmp_50;
      real_t tmp_278 = -tmp_267 - tmp_268 - tmp_270 - tmp_271 - tmp_272 - tmp_273 - tmp_275 - tmp_276 - tmp_277 + 1;
      real_t tmp_279 = tmp_272 + tmp_273 + tmp_277;
      real_t tmp_280 = tmp_268 + tmp_271 + tmp_276;
      real_t tmp_281 = tmp_267 + tmp_270 + tmp_275;
      real_t tmp_282 =
          tmp_64 * ( tmp_1 * ( tmp_279 - 1.0 / 4.0 ) + tmp_13 * ( tmp_281 - 1.0 / 4.0 ) + tmp_4 * ( tmp_280 - 1.0 / 4.0 ) );
      real_t tmp_283 = 0.042507265838595799 * tmp_66;
      real_t tmp_284 = tmp_19 * ( 0.039308471900058539 * tmp_21 + 0.37605877282253791 * tmp_22 + tmp_23 );
      real_t tmp_285 = tmp_284 * tmp_6;
      real_t tmp_286 = tmp_26 * tmp_284;
      real_t tmp_287 = tmp_19 * ( 0.039308471900058539 * tmp_30 + 0.37605877282253791 * tmp_31 + tmp_32 );
      real_t tmp_288 = tmp_28 * tmp_287;
      real_t tmp_289 = tmp_287 * tmp_35;
      real_t tmp_290 = tmp_284 * tmp_37;
      real_t tmp_291 = tmp_287 * tmp_39;
      real_t tmp_292 = tmp_19 * ( 0.039308471900058539 * tmp_43 + 0.37605877282253791 * tmp_44 + tmp_45 );
      real_t tmp_293 = tmp_292 * tmp_41;
      real_t tmp_294 = tmp_292 * tmp_48;
      real_t tmp_295 = tmp_292 * tmp_50;
      real_t tmp_296 = -tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1;
      real_t tmp_297 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_298 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_299 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_300 =
          tmp_64 * ( tmp_1 * ( tmp_297 - 1.0 / 4.0 ) + tmp_13 * ( tmp_299 - 1.0 / 4.0 ) + tmp_4 * ( tmp_298 - 1.0 / 4.0 ) );
      real_t tmp_301 = 0.020848748529055869 * tmp_66;
      real_t tmp_302 = tmp_19 * ( 0.93718850182767688 * tmp_21 + 0.031405749086161582 * tmp_22 + tmp_23 );
      real_t tmp_303 = tmp_302 * tmp_6;
      real_t tmp_304 = tmp_26 * tmp_302;
      real_t tmp_305 = tmp_19 * ( 0.93718850182767688 * tmp_30 + 0.031405749086161582 * tmp_31 + tmp_32 );
      real_t tmp_306 = tmp_28 * tmp_305;
      real_t tmp_307 = tmp_305 * tmp_35;
      real_t tmp_308 = tmp_302 * tmp_37;
      real_t tmp_309 = tmp_305 * tmp_39;
      real_t tmp_310 = tmp_19 * ( 0.93718850182767688 * tmp_43 + 0.031405749086161582 * tmp_44 + tmp_45 );
      real_t tmp_311 = tmp_310 * tmp_41;
      real_t tmp_312 = tmp_310 * tmp_48;
      real_t tmp_313 = tmp_310 * tmp_50;
      real_t tmp_314 = -tmp_303 - tmp_304 - tmp_306 - tmp_307 - tmp_308 - tmp_309 - tmp_311 - tmp_312 - tmp_313 + 1;
      real_t tmp_315 = tmp_308 + tmp_309 + tmp_313;
      real_t tmp_316 = tmp_304 + tmp_307 + tmp_312;
      real_t tmp_317 = tmp_303 + tmp_306 + tmp_311;
      real_t tmp_318 =
          tmp_64 * ( tmp_1 * ( tmp_315 - 1.0 / 4.0 ) + tmp_13 * ( tmp_317 - 1.0 / 4.0 ) + tmp_4 * ( tmp_316 - 1.0 / 4.0 ) );
      real_t tmp_319 = 0.0068572537431980923 * tmp_66;
      real_t tmp_320 = tmp_19 * ( 0.60796128279561268 * tmp_21 + 0.19601935860219369 * tmp_22 + tmp_23 );
      real_t tmp_321 = tmp_320 * tmp_6;
      real_t tmp_322 = tmp_26 * tmp_320;
      real_t tmp_323 = tmp_19 * ( 0.60796128279561268 * tmp_30 + 0.19601935860219369 * tmp_31 + tmp_32 );
      real_t tmp_324 = tmp_28 * tmp_323;
      real_t tmp_325 = tmp_323 * tmp_35;
      real_t tmp_326 = tmp_320 * tmp_37;
      real_t tmp_327 = tmp_323 * tmp_39;
      real_t tmp_328 = tmp_19 * ( 0.60796128279561268 * tmp_43 + 0.19601935860219369 * tmp_44 + tmp_45 );
      real_t tmp_329 = tmp_328 * tmp_41;
      real_t tmp_330 = tmp_328 * tmp_48;
      real_t tmp_331 = tmp_328 * tmp_50;
      real_t tmp_332 = -tmp_321 - tmp_322 - tmp_324 - tmp_325 - tmp_326 - tmp_327 - tmp_329 - tmp_330 - tmp_331 + 1;
      real_t tmp_333 = tmp_326 + tmp_327 + tmp_331;
      real_t tmp_334 = tmp_322 + tmp_325 + tmp_330;
      real_t tmp_335 = tmp_321 + tmp_324 + tmp_329;
      real_t tmp_336 =
          tmp_64 * ( tmp_1 * ( tmp_333 - 1.0 / 4.0 ) + tmp_13 * ( tmp_335 - 1.0 / 4.0 ) + tmp_4 * ( tmp_334 - 1.0 / 4.0 ) );
      real_t tmp_337 = 0.037198804536718075 * tmp_66;
      real_t tmp_338 = tmp_19 * ( 0.19107600050469298 * tmp_21 + 0.40446199974765351 * tmp_22 + tmp_23 );
      real_t tmp_339 = tmp_338 * tmp_6;
      real_t tmp_340 = tmp_26 * tmp_338;
      real_t tmp_341 = tmp_19 * ( 0.19107600050469298 * tmp_30 + 0.40446199974765351 * tmp_31 + tmp_32 );
      real_t tmp_342 = tmp_28 * tmp_341;
      real_t tmp_343 = tmp_341 * tmp_35;
      real_t tmp_344 = tmp_338 * tmp_37;
      real_t tmp_345 = tmp_341 * tmp_39;
      real_t tmp_346 = tmp_19 * ( 0.19107600050469298 * tmp_43 + 0.40446199974765351 * tmp_44 + tmp_45 );
      real_t tmp_347 = tmp_346 * tmp_41;
      real_t tmp_348 = tmp_346 * tmp_48;
      real_t tmp_349 = tmp_346 * tmp_50;
      real_t tmp_350 = -tmp_339 - tmp_340 - tmp_342 - tmp_343 - tmp_344 - tmp_345 - tmp_347 - tmp_348 - tmp_349 + 1;
      real_t tmp_351 = tmp_344 + tmp_345 + tmp_349;
      real_t tmp_352 = tmp_340 + tmp_343 + tmp_348;
      real_t tmp_353 = tmp_339 + tmp_342 + tmp_347;
      real_t tmp_354 =
          tmp_64 * ( tmp_1 * ( tmp_351 - 1.0 / 4.0 ) + tmp_13 * ( tmp_353 - 1.0 / 4.0 ) + tmp_4 * ( tmp_352 - 1.0 / 4.0 ) );
      real_t tmp_355 = 0.042507265838595799 * tmp_66;
      real_t tmp_356 = tmp_19 * ( 0.031405749086161582 * tmp_21 + 0.031405749086161582 * tmp_22 + tmp_23 );
      real_t tmp_357 = tmp_356 * tmp_6;
      real_t tmp_358 = tmp_26 * tmp_356;
      real_t tmp_359 = tmp_19 * ( 0.031405749086161582 * tmp_30 + 0.031405749086161582 * tmp_31 + tmp_32 );
      real_t tmp_360 = tmp_28 * tmp_359;
      real_t tmp_361 = tmp_35 * tmp_359;
      real_t tmp_362 = tmp_356 * tmp_37;
      real_t tmp_363 = tmp_359 * tmp_39;
      real_t tmp_364 = tmp_19 * ( 0.031405749086161582 * tmp_43 + 0.031405749086161582 * tmp_44 + tmp_45 );
      real_t tmp_365 = tmp_364 * tmp_41;
      real_t tmp_366 = tmp_364 * tmp_48;
      real_t tmp_367 = tmp_364 * tmp_50;
      real_t tmp_368 = -tmp_357 - tmp_358 - tmp_360 - tmp_361 - tmp_362 - tmp_363 - tmp_365 - tmp_366 - tmp_367 + 1;
      real_t tmp_369 = tmp_362 + tmp_363 + tmp_367;
      real_t tmp_370 = tmp_358 + tmp_361 + tmp_366;
      real_t tmp_371 = tmp_357 + tmp_360 + tmp_365;
      real_t tmp_372 =
          tmp_64 * ( tmp_1 * ( tmp_369 - 1.0 / 4.0 ) + tmp_13 * ( tmp_371 - 1.0 / 4.0 ) + tmp_4 * ( tmp_370 - 1.0 / 4.0 ) );
      real_t tmp_373 = 0.0068572537431980923 * tmp_66;
      real_t tmp_374 = tmp_19 * ( 0.19601935860219369 * tmp_21 + 0.19601935860219369 * tmp_22 + tmp_23 );
      real_t tmp_375 = tmp_374 * tmp_6;
      real_t tmp_376 = tmp_26 * tmp_374;
      real_t tmp_377 = tmp_19 * ( 0.19601935860219369 * tmp_30 + 0.19601935860219369 * tmp_31 + tmp_32 );
      real_t tmp_378 = tmp_28 * tmp_377;
      real_t tmp_379 = tmp_35 * tmp_377;
      real_t tmp_380 = tmp_37 * tmp_374;
      real_t tmp_381 = tmp_377 * tmp_39;
      real_t tmp_382 = tmp_19 * ( 0.19601935860219369 * tmp_43 + 0.19601935860219369 * tmp_44 + tmp_45 );
      real_t tmp_383 = tmp_382 * tmp_41;
      real_t tmp_384 = tmp_382 * tmp_48;
      real_t tmp_385 = tmp_382 * tmp_50;
      real_t tmp_386 = -tmp_375 - tmp_376 - tmp_378 - tmp_379 - tmp_380 - tmp_381 - tmp_383 - tmp_384 - tmp_385 + 1;
      real_t tmp_387 = tmp_380 + tmp_381 + tmp_385;
      real_t tmp_388 = tmp_376 + tmp_379 + tmp_384;
      real_t tmp_389 = tmp_375 + tmp_378 + tmp_383;
      real_t tmp_390 =
          tmp_64 * ( tmp_1 * ( tmp_387 - 1.0 / 4.0 ) + tmp_13 * ( tmp_389 - 1.0 / 4.0 ) + tmp_4 * ( tmp_388 - 1.0 / 4.0 ) );
      real_t tmp_391 = 0.037198804536718075 * tmp_66;
      real_t tmp_392 = tmp_19 * ( 0.40446199974765351 * tmp_21 + 0.40446199974765351 * tmp_22 + tmp_23 );
      real_t tmp_393 = tmp_392 * tmp_6;
      real_t tmp_394 = tmp_26 * tmp_392;
      real_t tmp_395 = tmp_19 * ( 0.40446199974765351 * tmp_30 + 0.40446199974765351 * tmp_31 + tmp_32 );
      real_t tmp_396 = tmp_28 * tmp_395;
      real_t tmp_397 = tmp_35 * tmp_395;
      real_t tmp_398 = tmp_37 * tmp_392;
      real_t tmp_399 = tmp_39 * tmp_395;
      real_t tmp_400 = tmp_19 * ( 0.40446199974765351 * tmp_43 + 0.40446199974765351 * tmp_44 + tmp_45 );
      real_t tmp_401 = tmp_400 * tmp_41;
      real_t tmp_402 = tmp_400 * tmp_48;
      real_t tmp_403 = tmp_400 * tmp_50;
      real_t tmp_404 = -tmp_393 - tmp_394 - tmp_396 - tmp_397 - tmp_398 - tmp_399 - tmp_401 - tmp_402 - tmp_403 + 1;
      real_t tmp_405 = tmp_398 + tmp_399 + tmp_403;
      real_t tmp_406 = tmp_394 + tmp_397 + tmp_402;
      real_t tmp_407 = tmp_393 + tmp_396 + tmp_401;
      real_t tmp_408 =
          tmp_64 * ( tmp_1 * ( tmp_405 - 1.0 / 4.0 ) + tmp_13 * ( tmp_407 - 1.0 / 4.0 ) + tmp_4 * ( tmp_406 - 1.0 / 4.0 ) );
      real_t tmp_409 = 0.042507265838595799 * tmp_66;
      real_t tmp_410 = tmp_19 * ( 0.1711304259088916 * tmp_21 + 0.041227165399737475 * tmp_22 + tmp_23 );
      real_t tmp_411 = tmp_410 * tmp_6;
      real_t tmp_412 = tmp_26 * tmp_410;
      real_t tmp_413 = tmp_19 * ( 0.1711304259088916 * tmp_30 + 0.041227165399737475 * tmp_31 + tmp_32 );
      real_t tmp_414 = tmp_28 * tmp_413;
      real_t tmp_415 = tmp_35 * tmp_413;
      real_t tmp_416 = tmp_37 * tmp_410;
      real_t tmp_417 = tmp_39 * tmp_413;
      real_t tmp_418 = tmp_19 * ( 0.1711304259088916 * tmp_43 + 0.041227165399737475 * tmp_44 + tmp_45 );
      real_t tmp_419 = tmp_41 * tmp_418;
      real_t tmp_420 = tmp_418 * tmp_48;
      real_t tmp_421 = tmp_418 * tmp_50;
      real_t tmp_422 = -tmp_411 - tmp_412 - tmp_414 - tmp_415 - tmp_416 - tmp_417 - tmp_419 - tmp_420 - tmp_421 + 1;
      real_t tmp_423 = tmp_416 + tmp_417 + tmp_421;
      real_t tmp_424 = tmp_412 + tmp_415 + tmp_420;
      real_t tmp_425 = tmp_411 + tmp_414 + tmp_419;
      real_t tmp_426 =
          tmp_64 * ( tmp_1 * ( tmp_423 - 1.0 / 4.0 ) + tmp_13 * ( tmp_425 - 1.0 / 4.0 ) + tmp_4 * ( tmp_424 - 1.0 / 4.0 ) );
      real_t tmp_427 = 0.019202922745021479 * tmp_66;
      real_t a_0_0   = tmp_103 * ( tmp_102 * tmp_98 - tmp_56 * tmp_98 ) + tmp_121 * ( tmp_116 * tmp_120 - tmp_116 * tmp_56 ) +
                     tmp_139 * ( tmp_134 * tmp_138 - tmp_134 * tmp_56 ) + tmp_157 * ( tmp_152 * tmp_156 - tmp_152 * tmp_56 ) +
                     tmp_175 * ( tmp_170 * tmp_174 - tmp_170 * tmp_56 ) + tmp_193 * ( tmp_188 * tmp_192 - tmp_188 * tmp_56 ) +
                     tmp_211 * ( tmp_206 * tmp_210 - tmp_206 * tmp_56 ) + tmp_229 * ( tmp_224 * tmp_228 - tmp_224 * tmp_56 ) +
                     tmp_247 * ( tmp_242 * tmp_246 - tmp_242 * tmp_56 ) + tmp_265 * ( tmp_260 * tmp_264 - tmp_260 * tmp_56 ) +
                     tmp_283 * ( tmp_278 * tmp_282 - tmp_278 * tmp_56 ) + tmp_301 * ( tmp_296 * tmp_300 - tmp_296 * tmp_56 ) +
                     tmp_319 * ( tmp_314 * tmp_318 - tmp_314 * tmp_56 ) + tmp_337 * ( tmp_332 * tmp_336 - tmp_332 * tmp_56 ) +
                     tmp_355 * ( tmp_350 * tmp_354 - tmp_350 * tmp_56 ) + tmp_373 * ( tmp_368 * tmp_372 - tmp_368 * tmp_56 ) +
                     tmp_391 * ( tmp_386 * tmp_390 - tmp_386 * tmp_56 ) + tmp_409 * ( tmp_404 * tmp_408 - tmp_404 * tmp_56 ) +
                     tmp_427 * ( tmp_422 * tmp_426 - tmp_422 * tmp_56 ) + tmp_67 * ( -tmp_52 * tmp_56 + tmp_52 * tmp_65 ) +
                     tmp_85 * ( -tmp_56 * tmp_80 + tmp_80 * tmp_84 );
      real_t a_1_0 = tmp_103 * ( tmp_102 * tmp_99 - tmp_56 * tmp_99 ) + tmp_121 * ( tmp_117 * tmp_120 - tmp_117 * tmp_56 ) +
                     tmp_139 * ( tmp_135 * tmp_138 - tmp_135 * tmp_56 ) + tmp_157 * ( tmp_153 * tmp_156 - tmp_153 * tmp_56 ) +
                     tmp_175 * ( tmp_171 * tmp_174 - tmp_171 * tmp_56 ) + tmp_193 * ( tmp_189 * tmp_192 - tmp_189 * tmp_56 ) +
                     tmp_211 * ( tmp_207 * tmp_210 - tmp_207 * tmp_56 ) + tmp_229 * ( tmp_225 * tmp_228 - tmp_225 * tmp_56 ) +
                     tmp_247 * ( tmp_243 * tmp_246 - tmp_243 * tmp_56 ) + tmp_265 * ( tmp_261 * tmp_264 - tmp_261 * tmp_56 ) +
                     tmp_283 * ( tmp_279 * tmp_282 - tmp_279 * tmp_56 ) + tmp_301 * ( tmp_297 * tmp_300 - tmp_297 * tmp_56 ) +
                     tmp_319 * ( tmp_315 * tmp_318 - tmp_315 * tmp_56 ) + tmp_337 * ( tmp_333 * tmp_336 - tmp_333 * tmp_56 ) +
                     tmp_355 * ( tmp_351 * tmp_354 - tmp_351 * tmp_56 ) + tmp_373 * ( tmp_369 * tmp_372 - tmp_369 * tmp_56 ) +
                     tmp_391 * ( tmp_387 * tmp_390 - tmp_387 * tmp_56 ) + tmp_409 * ( tmp_405 * tmp_408 - tmp_405 * tmp_56 ) +
                     tmp_427 * ( tmp_423 * tmp_426 - tmp_423 * tmp_56 ) + tmp_67 * ( -tmp_56 * tmp_57 + tmp_57 * tmp_65 ) +
                     tmp_85 * ( -tmp_56 * tmp_81 + tmp_81 * tmp_84 );
      real_t a_2_0 = tmp_103 * ( tmp_100 * tmp_102 - tmp_100 * tmp_56 ) + tmp_121 * ( tmp_118 * tmp_120 - tmp_118 * tmp_56 ) +
                     tmp_139 * ( tmp_136 * tmp_138 - tmp_136 * tmp_56 ) + tmp_157 * ( tmp_154 * tmp_156 - tmp_154 * tmp_56 ) +
                     tmp_175 * ( tmp_172 * tmp_174 - tmp_172 * tmp_56 ) + tmp_193 * ( tmp_190 * tmp_192 - tmp_190 * tmp_56 ) +
                     tmp_211 * ( tmp_208 * tmp_210 - tmp_208 * tmp_56 ) + tmp_229 * ( tmp_226 * tmp_228 - tmp_226 * tmp_56 ) +
                     tmp_247 * ( tmp_244 * tmp_246 - tmp_244 * tmp_56 ) + tmp_265 * ( tmp_262 * tmp_264 - tmp_262 * tmp_56 ) +
                     tmp_283 * ( tmp_280 * tmp_282 - tmp_280 * tmp_56 ) + tmp_301 * ( tmp_298 * tmp_300 - tmp_298 * tmp_56 ) +
                     tmp_319 * ( tmp_316 * tmp_318 - tmp_316 * tmp_56 ) + tmp_337 * ( tmp_334 * tmp_336 - tmp_334 * tmp_56 ) +
                     tmp_355 * ( tmp_352 * tmp_354 - tmp_352 * tmp_56 ) + tmp_373 * ( tmp_370 * tmp_372 - tmp_370 * tmp_56 ) +
                     tmp_391 * ( tmp_388 * tmp_390 - tmp_388 * tmp_56 ) + tmp_409 * ( tmp_406 * tmp_408 - tmp_406 * tmp_56 ) +
                     tmp_427 * ( tmp_424 * tmp_426 - tmp_424 * tmp_56 ) + tmp_67 * ( -tmp_56 * tmp_58 + tmp_58 * tmp_65 ) +
                     tmp_85 * ( -tmp_56 * tmp_82 + tmp_82 * tmp_84 );
      real_t a_3_0 = tmp_103 * ( tmp_101 * tmp_102 - tmp_101 * tmp_56 ) + tmp_121 * ( tmp_119 * tmp_120 - tmp_119 * tmp_56 ) +
                     tmp_139 * ( tmp_137 * tmp_138 - tmp_137 * tmp_56 ) + tmp_157 * ( tmp_155 * tmp_156 - tmp_155 * tmp_56 ) +
                     tmp_175 * ( tmp_173 * tmp_174 - tmp_173 * tmp_56 ) + tmp_193 * ( tmp_191 * tmp_192 - tmp_191 * tmp_56 ) +
                     tmp_211 * ( tmp_209 * tmp_210 - tmp_209 * tmp_56 ) + tmp_229 * ( tmp_227 * tmp_228 - tmp_227 * tmp_56 ) +
                     tmp_247 * ( tmp_245 * tmp_246 - tmp_245 * tmp_56 ) + tmp_265 * ( tmp_263 * tmp_264 - tmp_263 * tmp_56 ) +
                     tmp_283 * ( tmp_281 * tmp_282 - tmp_281 * tmp_56 ) + tmp_301 * ( tmp_299 * tmp_300 - tmp_299 * tmp_56 ) +
                     tmp_319 * ( tmp_317 * tmp_318 - tmp_317 * tmp_56 ) + tmp_337 * ( tmp_335 * tmp_336 - tmp_335 * tmp_56 ) +
                     tmp_355 * ( tmp_353 * tmp_354 - tmp_353 * tmp_56 ) + tmp_373 * ( tmp_371 * tmp_372 - tmp_371 * tmp_56 ) +
                     tmp_391 * ( tmp_389 * tmp_390 - tmp_389 * tmp_56 ) + tmp_409 * ( tmp_407 * tmp_408 - tmp_407 * tmp_56 ) +
                     tmp_427 * ( tmp_425 * tmp_426 - tmp_425 * tmp_56 ) + tmp_67 * ( -tmp_56 * tmp_59 + tmp_59 * tmp_65 ) +
                     tmp_85 * ( -tmp_56 * tmp_83 + tmp_83 * tmp_84 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
      elMat( 3, 0 ) = a_3_0;
   }

   void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                  const std::vector< Point3D >& coordsElementOuter,
                                  const std::vector< Point3D >& coordsFacet,
                                  const Point3D&,
                                  const Point3D&,
                                  const Point3D&     outwardNormal,
                                  const DGBasisInfo& trialBasis,
                                  const DGBasisInfo& testBasis,
                                  int                trialDegree,
                                  int                testDegree,
                                  MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_0;
      real_t tmp_6  = p_affine_1_1 + tmp_2;
      real_t tmp_7  = tmp_5 * tmp_6;
      real_t tmp_8  = tmp_4 - tmp_7;
      real_t tmp_9  = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = p_affine_1_2 + tmp_9;
      real_t tmp_13 = tmp_12 * tmp_5;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = p_affine_2_2 + tmp_9;
      real_t tmp_16 = tmp_15 * tmp_6;
      real_t tmp_17 = tmp_1 * tmp_11;
      real_t tmp_18 = tmp_12 * tmp_14;
      real_t tmp_19 =
          1.0 / ( tmp_10 * tmp_4 - tmp_10 * tmp_7 + tmp_11 * tmp_13 + tmp_14 * tmp_16 - tmp_15 * tmp_17 - tmp_18 * tmp_3 );
      real_t tmp_20 = p_affine_8_2 + tmp_9;
      real_t tmp_21 = -p_affine_8_2;
      real_t tmp_22 = p_affine_9_2 + tmp_21;
      real_t tmp_23 = p_affine_10_2 + tmp_21;
      real_t tmp_24 = 0.031405749086161582 * tmp_22 + 0.93718850182767688 * tmp_23;
      real_t tmp_25 = tmp_19 * ( tmp_20 + tmp_24 );
      real_t tmp_26 = tmp_25 * tmp_8;
      real_t tmp_27 = tmp_14 * tmp_6 - tmp_17;
      real_t tmp_28 = tmp_25 * tmp_27;
      real_t tmp_29 = -tmp_1 * tmp_15 + tmp_13;
      real_t tmp_30 = p_affine_8_1 + tmp_2;
      real_t tmp_31 = -p_affine_8_1;
      real_t tmp_32 = p_affine_9_1 + tmp_31;
      real_t tmp_33 = p_affine_10_1 + tmp_31;
      real_t tmp_34 = 0.031405749086161582 * tmp_32 + 0.93718850182767688 * tmp_33;
      real_t tmp_35 = tmp_19 * ( tmp_30 + tmp_34 );
      real_t tmp_36 = tmp_29 * tmp_35;
      real_t tmp_37 = tmp_1 * tmp_10 - tmp_18;
      real_t tmp_38 = tmp_35 * tmp_37;
      real_t tmp_39 = tmp_11 * tmp_5 - tmp_14 * tmp_3;
      real_t tmp_40 = tmp_25 * tmp_39;
      real_t tmp_41 = -tmp_10 * tmp_5 + tmp_14 * tmp_15;
      real_t tmp_42 = tmp_35 * tmp_41;
      real_t tmp_43 = -tmp_12 * tmp_3 + tmp_16;
      real_t tmp_44 = p_affine_8_0 + tmp_0;
      real_t tmp_45 = -p_affine_8_0;
      real_t tmp_46 = p_affine_9_0 + tmp_45;
      real_t tmp_47 = p_affine_10_0 + tmp_45;
      real_t tmp_48 = 0.031405749086161582 * tmp_46 + 0.93718850182767688 * tmp_47;
      real_t tmp_49 = tmp_19 * ( tmp_44 + tmp_48 );
      real_t tmp_50 = tmp_43 * tmp_49;
      real_t tmp_51 = -tmp_10 * tmp_6 + tmp_11 * tmp_12;
      real_t tmp_52 = tmp_49 * tmp_51;
      real_t tmp_53 = tmp_10 * tmp_3 - tmp_11 * tmp_15;
      real_t tmp_54 = tmp_49 * tmp_53;
      real_t tmp_55 = -tmp_26 - tmp_28 - tmp_36 - tmp_38 - tmp_40 - tmp_42 - tmp_50 - tmp_52 - tmp_54 + 1;
      real_t tmp_56 = -p_affine_4_1;
      real_t tmp_57 = p_affine_6_1 + tmp_56;
      real_t tmp_58 = -p_affine_4_2;
      real_t tmp_59 = p_affine_7_2 + tmp_58;
      real_t tmp_60 = tmp_57 * tmp_59;
      real_t tmp_61 = p_affine_7_1 + tmp_56;
      real_t tmp_62 = p_affine_6_2 + tmp_58;
      real_t tmp_63 = tmp_61 * tmp_62;
      real_t tmp_64 = tmp_60 - tmp_63;
      real_t tmp_65 = -p_affine_4_0;
      real_t tmp_66 = p_affine_5_0 + tmp_65;
      real_t tmp_67 = p_affine_6_0 + tmp_65;
      real_t tmp_68 = p_affine_5_2 + tmp_58;
      real_t tmp_69 = tmp_61 * tmp_68;
      real_t tmp_70 = p_affine_7_0 + tmp_65;
      real_t tmp_71 = p_affine_5_1 + tmp_56;
      real_t tmp_72 = tmp_62 * tmp_71;
      real_t tmp_73 = tmp_59 * tmp_71;
      real_t tmp_74 = tmp_57 * tmp_68;
      real_t tmp_75 =
          1.0 / ( tmp_60 * tmp_66 - tmp_63 * tmp_66 + tmp_67 * tmp_69 - tmp_67 * tmp_73 + tmp_70 * tmp_72 - tmp_70 * tmp_74 );
      real_t tmp_76 = tmp_66 * tmp_75;
      real_t tmp_77 = tmp_69 - tmp_73;
      real_t tmp_78 = tmp_67 * tmp_75;
      real_t tmp_79 = tmp_72 - tmp_74;
      real_t tmp_80 = tmp_70 * tmp_75;
      real_t tmp_81 = -tmp_59 * tmp_67 + tmp_62 * tmp_70;
      real_t tmp_82 = tmp_59 * tmp_66 - tmp_68 * tmp_70;
      real_t tmp_83 = -tmp_62 * tmp_66 + tmp_67 * tmp_68;
      real_t tmp_84 = -tmp_57 * tmp_70 + tmp_61 * tmp_67;
      real_t tmp_85 = -tmp_61 * tmp_66 + tmp_70 * tmp_71;
      real_t tmp_86 = tmp_57 * tmp_66 - tmp_67 * tmp_71;
      real_t tmp_87 = 0.5 * p_affine_13_0 * ( tmp_64 * tmp_76 + tmp_77 * tmp_78 + tmp_79 * tmp_80 ) +
                      0.5 * p_affine_13_1 * ( tmp_76 * tmp_81 + tmp_78 * tmp_82 + tmp_80 * tmp_83 ) +
                      0.5 * p_affine_13_2 * ( tmp_76 * tmp_84 + tmp_78 * tmp_85 + tmp_80 * tmp_86 );
      real_t tmp_88 = p_affine_8_2 + tmp_58;
      real_t tmp_89 = tmp_75 * ( tmp_24 + tmp_88 );
      real_t tmp_90 = p_affine_8_1 + tmp_56;
      real_t tmp_91 = tmp_75 * ( tmp_34 + tmp_90 );
      real_t tmp_92 = p_affine_8_0 + tmp_65;
      real_t tmp_93 = tmp_75 * ( tmp_48 + tmp_92 );
      real_t tmp_94 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_95 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_96 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_97 = ( std::abs( tmp_23 * tmp_94 - tmp_33 * tmp_96 ) * std::abs( tmp_23 * tmp_94 - tmp_33 * tmp_96 ) ) +
                      ( std::abs( tmp_23 * tmp_95 - tmp_47 * tmp_96 ) * std::abs( tmp_23 * tmp_95 - tmp_47 * tmp_96 ) ) +
                      ( std::abs( tmp_33 * tmp_95 - tmp_47 * tmp_94 ) * std::abs( tmp_33 * tmp_95 - tmp_47 * tmp_94 ) );
      real_t tmp_98  = 1.0 * std::pow( tmp_97, -0.25 );
      real_t tmp_99  = tmp_98 * ( tmp_66 * ( tmp_64 * tmp_93 + tmp_81 * tmp_91 + tmp_84 * tmp_89 - 1.0 / 4.0 ) +
                                 tmp_67 * ( tmp_77 * tmp_93 + tmp_82 * tmp_91 + tmp_85 * tmp_89 - 1.0 / 4.0 ) +
                                 tmp_70 * ( tmp_79 * tmp_93 + tmp_83 * tmp_91 + tmp_86 * tmp_89 - 1.0 / 4.0 ) );
      real_t tmp_100 = 1.0 * std::pow( tmp_97, 1.0 / 2.0 );
      real_t tmp_101 = 0.0068572537431980923 * tmp_100;
      real_t tmp_102 = 0.19601935860219369 * tmp_22 + 0.60796128279561268 * tmp_23;
      real_t tmp_103 = tmp_19 * ( tmp_102 + tmp_20 );
      real_t tmp_104 = tmp_103 * tmp_8;
      real_t tmp_105 = tmp_103 * tmp_27;
      real_t tmp_106 = 0.19601935860219369 * tmp_32 + 0.60796128279561268 * tmp_33;
      real_t tmp_107 = tmp_19 * ( tmp_106 + tmp_30 );
      real_t tmp_108 = tmp_107 * tmp_29;
      real_t tmp_109 = tmp_107 * tmp_37;
      real_t tmp_110 = tmp_103 * tmp_39;
      real_t tmp_111 = tmp_107 * tmp_41;
      real_t tmp_112 = 0.19601935860219369 * tmp_46 + 0.60796128279561268 * tmp_47;
      real_t tmp_113 = tmp_19 * ( tmp_112 + tmp_44 );
      real_t tmp_114 = tmp_113 * tmp_43;
      real_t tmp_115 = tmp_113 * tmp_51;
      real_t tmp_116 = tmp_113 * tmp_53;
      real_t tmp_117 = -tmp_104 - tmp_105 - tmp_108 - tmp_109 - tmp_110 - tmp_111 - tmp_114 - tmp_115 - tmp_116 + 1;
      real_t tmp_118 = tmp_75 * ( tmp_102 + tmp_88 );
      real_t tmp_119 = tmp_75 * ( tmp_106 + tmp_90 );
      real_t tmp_120 = tmp_75 * ( tmp_112 + tmp_92 );
      real_t tmp_121 = tmp_98 * ( tmp_66 * ( tmp_118 * tmp_84 + tmp_119 * tmp_81 + tmp_120 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_118 * tmp_85 + tmp_119 * tmp_82 + tmp_120 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_118 * tmp_86 + tmp_119 * tmp_83 + tmp_120 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_122 = 0.037198804536718075 * tmp_100;
      real_t tmp_123 = 0.37605877282253791 * tmp_22 + 0.039308471900058539 * tmp_23;
      real_t tmp_124 = tmp_19 * ( tmp_123 + tmp_20 );
      real_t tmp_125 = tmp_124 * tmp_8;
      real_t tmp_126 = tmp_124 * tmp_27;
      real_t tmp_127 = 0.37605877282253791 * tmp_32 + 0.039308471900058539 * tmp_33;
      real_t tmp_128 = tmp_19 * ( tmp_127 + tmp_30 );
      real_t tmp_129 = tmp_128 * tmp_29;
      real_t tmp_130 = tmp_128 * tmp_37;
      real_t tmp_131 = tmp_124 * tmp_39;
      real_t tmp_132 = tmp_128 * tmp_41;
      real_t tmp_133 = 0.37605877282253791 * tmp_46 + 0.039308471900058539 * tmp_47;
      real_t tmp_134 = tmp_19 * ( tmp_133 + tmp_44 );
      real_t tmp_135 = tmp_134 * tmp_43;
      real_t tmp_136 = tmp_134 * tmp_51;
      real_t tmp_137 = tmp_134 * tmp_53;
      real_t tmp_138 = -tmp_125 - tmp_126 - tmp_129 - tmp_130 - tmp_131 - tmp_132 - tmp_135 - tmp_136 - tmp_137 + 1;
      real_t tmp_139 = tmp_75 * ( tmp_123 + tmp_88 );
      real_t tmp_140 = tmp_75 * ( tmp_127 + tmp_90 );
      real_t tmp_141 = tmp_75 * ( tmp_133 + tmp_92 );
      real_t tmp_142 = tmp_98 * ( tmp_66 * ( tmp_139 * tmp_84 + tmp_140 * tmp_81 + tmp_141 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_139 * tmp_85 + tmp_140 * tmp_82 + tmp_141 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_139 * tmp_86 + tmp_140 * tmp_83 + tmp_141 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_143 = 0.020848748529055869 * tmp_100;
      real_t tmp_144 = 0.78764240869137092 * tmp_22 + 0.1711304259088916 * tmp_23;
      real_t tmp_145 = tmp_19 * ( tmp_144 + tmp_20 );
      real_t tmp_146 = tmp_145 * tmp_8;
      real_t tmp_147 = tmp_145 * tmp_27;
      real_t tmp_148 = 0.78764240869137092 * tmp_32 + 0.1711304259088916 * tmp_33;
      real_t tmp_149 = tmp_19 * ( tmp_148 + tmp_30 );
      real_t tmp_150 = tmp_149 * tmp_29;
      real_t tmp_151 = tmp_149 * tmp_37;
      real_t tmp_152 = tmp_145 * tmp_39;
      real_t tmp_153 = tmp_149 * tmp_41;
      real_t tmp_154 = 0.78764240869137092 * tmp_46 + 0.1711304259088916 * tmp_47;
      real_t tmp_155 = tmp_19 * ( tmp_154 + tmp_44 );
      real_t tmp_156 = tmp_155 * tmp_43;
      real_t tmp_157 = tmp_155 * tmp_51;
      real_t tmp_158 = tmp_155 * tmp_53;
      real_t tmp_159 = -tmp_146 - tmp_147 - tmp_150 - tmp_151 - tmp_152 - tmp_153 - tmp_156 - tmp_157 - tmp_158 + 1;
      real_t tmp_160 = tmp_75 * ( tmp_144 + tmp_88 );
      real_t tmp_161 = tmp_75 * ( tmp_148 + tmp_90 );
      real_t tmp_162 = tmp_75 * ( tmp_154 + tmp_92 );
      real_t tmp_163 = tmp_98 * ( tmp_66 * ( tmp_160 * tmp_84 + tmp_161 * tmp_81 + tmp_162 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_160 * tmp_85 + tmp_161 * tmp_82 + tmp_162 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_160 * tmp_86 + tmp_161 * tmp_83 + tmp_162 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_164 = 0.019202922745021479 * tmp_100;
      real_t tmp_165 = 0.58463275527740355 * tmp_22 + 0.37605877282253791 * tmp_23;
      real_t tmp_166 = tmp_19 * ( tmp_165 + tmp_20 );
      real_t tmp_167 = tmp_166 * tmp_8;
      real_t tmp_168 = tmp_166 * tmp_27;
      real_t tmp_169 = 0.58463275527740355 * tmp_32 + 0.37605877282253791 * tmp_33;
      real_t tmp_170 = tmp_19 * ( tmp_169 + tmp_30 );
      real_t tmp_171 = tmp_170 * tmp_29;
      real_t tmp_172 = tmp_170 * tmp_37;
      real_t tmp_173 = tmp_166 * tmp_39;
      real_t tmp_174 = tmp_170 * tmp_41;
      real_t tmp_175 = 0.58463275527740355 * tmp_46 + 0.37605877282253791 * tmp_47;
      real_t tmp_176 = tmp_19 * ( tmp_175 + tmp_44 );
      real_t tmp_177 = tmp_176 * tmp_43;
      real_t tmp_178 = tmp_176 * tmp_51;
      real_t tmp_179 = tmp_176 * tmp_53;
      real_t tmp_180 = -tmp_167 - tmp_168 - tmp_171 - tmp_172 - tmp_173 - tmp_174 - tmp_177 - tmp_178 - tmp_179 + 1;
      real_t tmp_181 = tmp_75 * ( tmp_165 + tmp_88 );
      real_t tmp_182 = tmp_75 * ( tmp_169 + tmp_90 );
      real_t tmp_183 = tmp_75 * ( tmp_175 + tmp_92 );
      real_t tmp_184 = tmp_98 * ( tmp_66 * ( tmp_181 * tmp_84 + tmp_182 * tmp_81 + tmp_183 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_181 * tmp_85 + tmp_182 * tmp_82 + tmp_183 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_181 * tmp_86 + tmp_182 * tmp_83 + tmp_183 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_185 = 0.020848748529055869 * tmp_100;
      real_t tmp_186 = 0.041227165399737475 * tmp_22 + 0.78764240869137092 * tmp_23;
      real_t tmp_187 = tmp_19 * ( tmp_186 + tmp_20 );
      real_t tmp_188 = tmp_187 * tmp_8;
      real_t tmp_189 = tmp_187 * tmp_27;
      real_t tmp_190 = 0.041227165399737475 * tmp_32 + 0.78764240869137092 * tmp_33;
      real_t tmp_191 = tmp_19 * ( tmp_190 + tmp_30 );
      real_t tmp_192 = tmp_191 * tmp_29;
      real_t tmp_193 = tmp_191 * tmp_37;
      real_t tmp_194 = tmp_187 * tmp_39;
      real_t tmp_195 = tmp_191 * tmp_41;
      real_t tmp_196 = 0.041227165399737475 * tmp_46 + 0.78764240869137092 * tmp_47;
      real_t tmp_197 = tmp_19 * ( tmp_196 + tmp_44 );
      real_t tmp_198 = tmp_197 * tmp_43;
      real_t tmp_199 = tmp_197 * tmp_51;
      real_t tmp_200 = tmp_197 * tmp_53;
      real_t tmp_201 = -tmp_188 - tmp_189 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_198 - tmp_199 - tmp_200 + 1;
      real_t tmp_202 = tmp_75 * ( tmp_186 + tmp_88 );
      real_t tmp_203 = tmp_75 * ( tmp_190 + tmp_90 );
      real_t tmp_204 = tmp_75 * ( tmp_196 + tmp_92 );
      real_t tmp_205 = tmp_98 * ( tmp_66 * ( tmp_202 * tmp_84 + tmp_203 * tmp_81 + tmp_204 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_202 * tmp_85 + tmp_203 * tmp_82 + tmp_204 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_202 * tmp_86 + tmp_203 * tmp_83 + tmp_204 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_206 = 0.019202922745021479 * tmp_100;
      real_t tmp_207 = 0.039308471900058539 * tmp_22 + 0.58463275527740355 * tmp_23;
      real_t tmp_208 = tmp_19 * ( tmp_20 + tmp_207 );
      real_t tmp_209 = tmp_208 * tmp_8;
      real_t tmp_210 = tmp_208 * tmp_27;
      real_t tmp_211 = 0.039308471900058539 * tmp_32 + 0.58463275527740355 * tmp_33;
      real_t tmp_212 = tmp_19 * ( tmp_211 + tmp_30 );
      real_t tmp_213 = tmp_212 * tmp_29;
      real_t tmp_214 = tmp_212 * tmp_37;
      real_t tmp_215 = tmp_208 * tmp_39;
      real_t tmp_216 = tmp_212 * tmp_41;
      real_t tmp_217 = 0.039308471900058539 * tmp_46 + 0.58463275527740355 * tmp_47;
      real_t tmp_218 = tmp_19 * ( tmp_217 + tmp_44 );
      real_t tmp_219 = tmp_218 * tmp_43;
      real_t tmp_220 = tmp_218 * tmp_51;
      real_t tmp_221 = tmp_218 * tmp_53;
      real_t tmp_222 = -tmp_209 - tmp_210 - tmp_213 - tmp_214 - tmp_215 - tmp_216 - tmp_219 - tmp_220 - tmp_221 + 1;
      real_t tmp_223 = tmp_75 * ( tmp_207 + tmp_88 );
      real_t tmp_224 = tmp_75 * ( tmp_211 + tmp_90 );
      real_t tmp_225 = tmp_75 * ( tmp_217 + tmp_92 );
      real_t tmp_226 = tmp_98 * ( tmp_66 * ( tmp_223 * tmp_84 + tmp_224 * tmp_81 + tmp_225 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_223 * tmp_85 + tmp_224 * tmp_82 + tmp_225 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_223 * tmp_86 + tmp_224 * tmp_83 + tmp_225 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_227 = 0.020848748529055869 * tmp_100;
      real_t tmp_228 = 0.78764240869137092 * tmp_22 + 0.041227165399737475 * tmp_23;
      real_t tmp_229 = tmp_19 * ( tmp_20 + tmp_228 );
      real_t tmp_230 = tmp_229 * tmp_8;
      real_t tmp_231 = tmp_229 * tmp_27;
      real_t tmp_232 = 0.78764240869137092 * tmp_32 + 0.041227165399737475 * tmp_33;
      real_t tmp_233 = tmp_19 * ( tmp_232 + tmp_30 );
      real_t tmp_234 = tmp_233 * tmp_29;
      real_t tmp_235 = tmp_233 * tmp_37;
      real_t tmp_236 = tmp_229 * tmp_39;
      real_t tmp_237 = tmp_233 * tmp_41;
      real_t tmp_238 = 0.78764240869137092 * tmp_46 + 0.041227165399737475 * tmp_47;
      real_t tmp_239 = tmp_19 * ( tmp_238 + tmp_44 );
      real_t tmp_240 = tmp_239 * tmp_43;
      real_t tmp_241 = tmp_239 * tmp_51;
      real_t tmp_242 = tmp_239 * tmp_53;
      real_t tmp_243 = -tmp_230 - tmp_231 - tmp_234 - tmp_235 - tmp_236 - tmp_237 - tmp_240 - tmp_241 - tmp_242 + 1;
      real_t tmp_244 = tmp_75 * ( tmp_228 + tmp_88 );
      real_t tmp_245 = tmp_75 * ( tmp_232 + tmp_90 );
      real_t tmp_246 = tmp_75 * ( tmp_238 + tmp_92 );
      real_t tmp_247 = tmp_98 * ( tmp_66 * ( tmp_244 * tmp_84 + tmp_245 * tmp_81 + tmp_246 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_244 * tmp_85 + tmp_245 * tmp_82 + tmp_246 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_244 * tmp_86 + tmp_245 * tmp_83 + tmp_246 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_248 = 0.019202922745021479 * tmp_100;
      real_t tmp_249 = 0.58463275527740355 * tmp_22 + 0.039308471900058539 * tmp_23;
      real_t tmp_250 = tmp_19 * ( tmp_20 + tmp_249 );
      real_t tmp_251 = tmp_250 * tmp_8;
      real_t tmp_252 = tmp_250 * tmp_27;
      real_t tmp_253 = 0.58463275527740355 * tmp_32 + 0.039308471900058539 * tmp_33;
      real_t tmp_254 = tmp_19 * ( tmp_253 + tmp_30 );
      real_t tmp_255 = tmp_254 * tmp_29;
      real_t tmp_256 = tmp_254 * tmp_37;
      real_t tmp_257 = tmp_250 * tmp_39;
      real_t tmp_258 = tmp_254 * tmp_41;
      real_t tmp_259 = 0.58463275527740355 * tmp_46 + 0.039308471900058539 * tmp_47;
      real_t tmp_260 = tmp_19 * ( tmp_259 + tmp_44 );
      real_t tmp_261 = tmp_260 * tmp_43;
      real_t tmp_262 = tmp_260 * tmp_51;
      real_t tmp_263 = tmp_260 * tmp_53;
      real_t tmp_264 = -tmp_251 - tmp_252 - tmp_255 - tmp_256 - tmp_257 - tmp_258 - tmp_261 - tmp_262 - tmp_263 + 1;
      real_t tmp_265 = tmp_75 * ( tmp_249 + tmp_88 );
      real_t tmp_266 = tmp_75 * ( tmp_253 + tmp_90 );
      real_t tmp_267 = tmp_75 * ( tmp_259 + tmp_92 );
      real_t tmp_268 = tmp_98 * ( tmp_66 * ( tmp_265 * tmp_84 + tmp_266 * tmp_81 + tmp_267 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_265 * tmp_85 + tmp_266 * tmp_82 + tmp_267 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_265 * tmp_86 + tmp_266 * tmp_83 + tmp_267 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_269 = 0.020848748529055869 * tmp_100;
      real_t tmp_270 = 0.1711304259088916 * tmp_22 + 0.78764240869137092 * tmp_23;
      real_t tmp_271 = tmp_19 * ( tmp_20 + tmp_270 );
      real_t tmp_272 = tmp_271 * tmp_8;
      real_t tmp_273 = tmp_27 * tmp_271;
      real_t tmp_274 = 0.1711304259088916 * tmp_32 + 0.78764240869137092 * tmp_33;
      real_t tmp_275 = tmp_19 * ( tmp_274 + tmp_30 );
      real_t tmp_276 = tmp_275 * tmp_29;
      real_t tmp_277 = tmp_275 * tmp_37;
      real_t tmp_278 = tmp_271 * tmp_39;
      real_t tmp_279 = tmp_275 * tmp_41;
      real_t tmp_280 = 0.1711304259088916 * tmp_46 + 0.78764240869137092 * tmp_47;
      real_t tmp_281 = tmp_19 * ( tmp_280 + tmp_44 );
      real_t tmp_282 = tmp_281 * tmp_43;
      real_t tmp_283 = tmp_281 * tmp_51;
      real_t tmp_284 = tmp_281 * tmp_53;
      real_t tmp_285 = -tmp_272 - tmp_273 - tmp_276 - tmp_277 - tmp_278 - tmp_279 - tmp_282 - tmp_283 - tmp_284 + 1;
      real_t tmp_286 = tmp_75 * ( tmp_270 + tmp_88 );
      real_t tmp_287 = tmp_75 * ( tmp_274 + tmp_90 );
      real_t tmp_288 = tmp_75 * ( tmp_280 + tmp_92 );
      real_t tmp_289 = tmp_98 * ( tmp_66 * ( tmp_286 * tmp_84 + tmp_287 * tmp_81 + tmp_288 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_286 * tmp_85 + tmp_287 * tmp_82 + tmp_288 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_286 * tmp_86 + tmp_287 * tmp_83 + tmp_288 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_290 = 0.019202922745021479 * tmp_100;
      real_t tmp_291 = 0.37605877282253791 * tmp_22 + 0.58463275527740355 * tmp_23;
      real_t tmp_292 = tmp_19 * ( tmp_20 + tmp_291 );
      real_t tmp_293 = tmp_292 * tmp_8;
      real_t tmp_294 = tmp_27 * tmp_292;
      real_t tmp_295 = 0.37605877282253791 * tmp_32 + 0.58463275527740355 * tmp_33;
      real_t tmp_296 = tmp_19 * ( tmp_295 + tmp_30 );
      real_t tmp_297 = tmp_29 * tmp_296;
      real_t tmp_298 = tmp_296 * tmp_37;
      real_t tmp_299 = tmp_292 * tmp_39;
      real_t tmp_300 = tmp_296 * tmp_41;
      real_t tmp_301 = 0.37605877282253791 * tmp_46 + 0.58463275527740355 * tmp_47;
      real_t tmp_302 = tmp_19 * ( tmp_301 + tmp_44 );
      real_t tmp_303 = tmp_302 * tmp_43;
      real_t tmp_304 = tmp_302 * tmp_51;
      real_t tmp_305 = tmp_302 * tmp_53;
      real_t tmp_306 = -tmp_293 - tmp_294 - tmp_297 - tmp_298 - tmp_299 - tmp_300 - tmp_303 - tmp_304 - tmp_305 + 1;
      real_t tmp_307 = tmp_75 * ( tmp_291 + tmp_88 );
      real_t tmp_308 = tmp_75 * ( tmp_295 + tmp_90 );
      real_t tmp_309 = tmp_75 * ( tmp_301 + tmp_92 );
      real_t tmp_310 = tmp_98 * ( tmp_66 * ( tmp_307 * tmp_84 + tmp_308 * tmp_81 + tmp_309 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_307 * tmp_85 + tmp_308 * tmp_82 + tmp_309 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_307 * tmp_86 + tmp_308 * tmp_83 + tmp_309 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_311 = 0.020848748529055869 * tmp_100;
      real_t tmp_312 = 0.041227165399737475 * tmp_22 + 0.1711304259088916 * tmp_23;
      real_t tmp_313 = tmp_19 * ( tmp_20 + tmp_312 );
      real_t tmp_314 = tmp_313 * tmp_8;
      real_t tmp_315 = tmp_27 * tmp_313;
      real_t tmp_316 = 0.041227165399737475 * tmp_32 + 0.1711304259088916 * tmp_33;
      real_t tmp_317 = tmp_19 * ( tmp_30 + tmp_316 );
      real_t tmp_318 = tmp_29 * tmp_317;
      real_t tmp_319 = tmp_317 * tmp_37;
      real_t tmp_320 = tmp_313 * tmp_39;
      real_t tmp_321 = tmp_317 * tmp_41;
      real_t tmp_322 = 0.041227165399737475 * tmp_46 + 0.1711304259088916 * tmp_47;
      real_t tmp_323 = tmp_19 * ( tmp_322 + tmp_44 );
      real_t tmp_324 = tmp_323 * tmp_43;
      real_t tmp_325 = tmp_323 * tmp_51;
      real_t tmp_326 = tmp_323 * tmp_53;
      real_t tmp_327 = -tmp_314 - tmp_315 - tmp_318 - tmp_319 - tmp_320 - tmp_321 - tmp_324 - tmp_325 - tmp_326 + 1;
      real_t tmp_328 = tmp_75 * ( tmp_312 + tmp_88 );
      real_t tmp_329 = tmp_75 * ( tmp_316 + tmp_90 );
      real_t tmp_330 = tmp_75 * ( tmp_322 + tmp_92 );
      real_t tmp_331 = tmp_98 * ( tmp_66 * ( tmp_328 * tmp_84 + tmp_329 * tmp_81 + tmp_330 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_328 * tmp_85 + tmp_329 * tmp_82 + tmp_330 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_328 * tmp_86 + tmp_329 * tmp_83 + tmp_330 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_332 = 0.019202922745021479 * tmp_100;
      real_t tmp_333 = 0.40446199974765351 * tmp_22 + 0.19107600050469298 * tmp_23;
      real_t tmp_334 = tmp_19 * ( tmp_20 + tmp_333 );
      real_t tmp_335 = tmp_334 * tmp_8;
      real_t tmp_336 = tmp_27 * tmp_334;
      real_t tmp_337 = 0.40446199974765351 * tmp_32 + 0.19107600050469298 * tmp_33;
      real_t tmp_338 = tmp_19 * ( tmp_30 + tmp_337 );
      real_t tmp_339 = tmp_29 * tmp_338;
      real_t tmp_340 = tmp_338 * tmp_37;
      real_t tmp_341 = tmp_334 * tmp_39;
      real_t tmp_342 = tmp_338 * tmp_41;
      real_t tmp_343 = 0.40446199974765351 * tmp_46 + 0.19107600050469298 * tmp_47;
      real_t tmp_344 = tmp_19 * ( tmp_343 + tmp_44 );
      real_t tmp_345 = tmp_344 * tmp_43;
      real_t tmp_346 = tmp_344 * tmp_51;
      real_t tmp_347 = tmp_344 * tmp_53;
      real_t tmp_348 = -tmp_335 - tmp_336 - tmp_339 - tmp_340 - tmp_341 - tmp_342 - tmp_345 - tmp_346 - tmp_347 + 1;
      real_t tmp_349 = tmp_75 * ( tmp_333 + tmp_88 );
      real_t tmp_350 = tmp_75 * ( tmp_337 + tmp_90 );
      real_t tmp_351 = tmp_75 * ( tmp_343 + tmp_92 );
      real_t tmp_352 = tmp_98 * ( tmp_66 * ( tmp_349 * tmp_84 + tmp_350 * tmp_81 + tmp_351 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_349 * tmp_85 + tmp_350 * tmp_82 + tmp_351 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_349 * tmp_86 + tmp_350 * tmp_83 + tmp_351 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_353 = 0.042507265838595799 * tmp_100;
      real_t tmp_354 = 0.039308471900058539 * tmp_22 + 0.37605877282253791 * tmp_23;
      real_t tmp_355 = tmp_19 * ( tmp_20 + tmp_354 );
      real_t tmp_356 = tmp_355 * tmp_8;
      real_t tmp_357 = tmp_27 * tmp_355;
      real_t tmp_358 = 0.039308471900058539 * tmp_32 + 0.37605877282253791 * tmp_33;
      real_t tmp_359 = tmp_19 * ( tmp_30 + tmp_358 );
      real_t tmp_360 = tmp_29 * tmp_359;
      real_t tmp_361 = tmp_359 * tmp_37;
      real_t tmp_362 = tmp_355 * tmp_39;
      real_t tmp_363 = tmp_359 * tmp_41;
      real_t tmp_364 = 0.039308471900058539 * tmp_46 + 0.37605877282253791 * tmp_47;
      real_t tmp_365 = tmp_19 * ( tmp_364 + tmp_44 );
      real_t tmp_366 = tmp_365 * tmp_43;
      real_t tmp_367 = tmp_365 * tmp_51;
      real_t tmp_368 = tmp_365 * tmp_53;
      real_t tmp_369 = -tmp_356 - tmp_357 - tmp_360 - tmp_361 - tmp_362 - tmp_363 - tmp_366 - tmp_367 - tmp_368 + 1;
      real_t tmp_370 = tmp_75 * ( tmp_354 + tmp_88 );
      real_t tmp_371 = tmp_75 * ( tmp_358 + tmp_90 );
      real_t tmp_372 = tmp_75 * ( tmp_364 + tmp_92 );
      real_t tmp_373 = tmp_98 * ( tmp_66 * ( tmp_370 * tmp_84 + tmp_371 * tmp_81 + tmp_372 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_370 * tmp_85 + tmp_371 * tmp_82 + tmp_372 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_370 * tmp_86 + tmp_371 * tmp_83 + tmp_372 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_374 = 0.020848748529055869 * tmp_100;
      real_t tmp_375 = 0.93718850182767688 * tmp_22 + 0.031405749086161582 * tmp_23;
      real_t tmp_376 = tmp_19 * ( tmp_20 + tmp_375 );
      real_t tmp_377 = tmp_376 * tmp_8;
      real_t tmp_378 = tmp_27 * tmp_376;
      real_t tmp_379 = 0.93718850182767688 * tmp_32 + 0.031405749086161582 * tmp_33;
      real_t tmp_380 = tmp_19 * ( tmp_30 + tmp_379 );
      real_t tmp_381 = tmp_29 * tmp_380;
      real_t tmp_382 = tmp_37 * tmp_380;
      real_t tmp_383 = tmp_376 * tmp_39;
      real_t tmp_384 = tmp_380 * tmp_41;
      real_t tmp_385 = 0.93718850182767688 * tmp_46 + 0.031405749086161582 * tmp_47;
      real_t tmp_386 = tmp_19 * ( tmp_385 + tmp_44 );
      real_t tmp_387 = tmp_386 * tmp_43;
      real_t tmp_388 = tmp_386 * tmp_51;
      real_t tmp_389 = tmp_386 * tmp_53;
      real_t tmp_390 = -tmp_377 - tmp_378 - tmp_381 - tmp_382 - tmp_383 - tmp_384 - tmp_387 - tmp_388 - tmp_389 + 1;
      real_t tmp_391 = tmp_75 * ( tmp_375 + tmp_88 );
      real_t tmp_392 = tmp_75 * ( tmp_379 + tmp_90 );
      real_t tmp_393 = tmp_75 * ( tmp_385 + tmp_92 );
      real_t tmp_394 = tmp_98 * ( tmp_66 * ( tmp_391 * tmp_84 + tmp_392 * tmp_81 + tmp_393 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_391 * tmp_85 + tmp_392 * tmp_82 + tmp_393 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_391 * tmp_86 + tmp_392 * tmp_83 + tmp_393 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_395 = 0.0068572537431980923 * tmp_100;
      real_t tmp_396 = 0.60796128279561268 * tmp_22 + 0.19601935860219369 * tmp_23;
      real_t tmp_397 = tmp_19 * ( tmp_20 + tmp_396 );
      real_t tmp_398 = tmp_397 * tmp_8;
      real_t tmp_399 = tmp_27 * tmp_397;
      real_t tmp_400 = 0.60796128279561268 * tmp_32 + 0.19601935860219369 * tmp_33;
      real_t tmp_401 = tmp_19 * ( tmp_30 + tmp_400 );
      real_t tmp_402 = tmp_29 * tmp_401;
      real_t tmp_403 = tmp_37 * tmp_401;
      real_t tmp_404 = tmp_39 * tmp_397;
      real_t tmp_405 = tmp_401 * tmp_41;
      real_t tmp_406 = 0.60796128279561268 * tmp_46 + 0.19601935860219369 * tmp_47;
      real_t tmp_407 = tmp_19 * ( tmp_406 + tmp_44 );
      real_t tmp_408 = tmp_407 * tmp_43;
      real_t tmp_409 = tmp_407 * tmp_51;
      real_t tmp_410 = tmp_407 * tmp_53;
      real_t tmp_411 = -tmp_398 - tmp_399 - tmp_402 - tmp_403 - tmp_404 - tmp_405 - tmp_408 - tmp_409 - tmp_410 + 1;
      real_t tmp_412 = tmp_75 * ( tmp_396 + tmp_88 );
      real_t tmp_413 = tmp_75 * ( tmp_400 + tmp_90 );
      real_t tmp_414 = tmp_75 * ( tmp_406 + tmp_92 );
      real_t tmp_415 = tmp_98 * ( tmp_66 * ( tmp_412 * tmp_84 + tmp_413 * tmp_81 + tmp_414 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_412 * tmp_85 + tmp_413 * tmp_82 + tmp_414 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_412 * tmp_86 + tmp_413 * tmp_83 + tmp_414 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_416 = 0.037198804536718075 * tmp_100;
      real_t tmp_417 = 0.19107600050469298 * tmp_22 + 0.40446199974765351 * tmp_23;
      real_t tmp_418 = tmp_19 * ( tmp_20 + tmp_417 );
      real_t tmp_419 = tmp_418 * tmp_8;
      real_t tmp_420 = tmp_27 * tmp_418;
      real_t tmp_421 = 0.19107600050469298 * tmp_32 + 0.40446199974765351 * tmp_33;
      real_t tmp_422 = tmp_19 * ( tmp_30 + tmp_421 );
      real_t tmp_423 = tmp_29 * tmp_422;
      real_t tmp_424 = tmp_37 * tmp_422;
      real_t tmp_425 = tmp_39 * tmp_418;
      real_t tmp_426 = tmp_41 * tmp_422;
      real_t tmp_427 = 0.19107600050469298 * tmp_46 + 0.40446199974765351 * tmp_47;
      real_t tmp_428 = tmp_19 * ( tmp_427 + tmp_44 );
      real_t tmp_429 = tmp_428 * tmp_43;
      real_t tmp_430 = tmp_428 * tmp_51;
      real_t tmp_431 = tmp_428 * tmp_53;
      real_t tmp_432 = -tmp_419 - tmp_420 - tmp_423 - tmp_424 - tmp_425 - tmp_426 - tmp_429 - tmp_430 - tmp_431 + 1;
      real_t tmp_433 = tmp_75 * ( tmp_417 + tmp_88 );
      real_t tmp_434 = tmp_75 * ( tmp_421 + tmp_90 );
      real_t tmp_435 = tmp_75 * ( tmp_427 + tmp_92 );
      real_t tmp_436 = tmp_98 * ( tmp_66 * ( tmp_433 * tmp_84 + tmp_434 * tmp_81 + tmp_435 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_433 * tmp_85 + tmp_434 * tmp_82 + tmp_435 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_433 * tmp_86 + tmp_434 * tmp_83 + tmp_435 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_437 = 0.042507265838595799 * tmp_100;
      real_t tmp_438 = 0.031405749086161582 * tmp_22 + 0.031405749086161582 * tmp_23;
      real_t tmp_439 = tmp_19 * ( tmp_20 + tmp_438 );
      real_t tmp_440 = tmp_439 * tmp_8;
      real_t tmp_441 = tmp_27 * tmp_439;
      real_t tmp_442 = 0.031405749086161582 * tmp_32 + 0.031405749086161582 * tmp_33;
      real_t tmp_443 = tmp_19 * ( tmp_30 + tmp_442 );
      real_t tmp_444 = tmp_29 * tmp_443;
      real_t tmp_445 = tmp_37 * tmp_443;
      real_t tmp_446 = tmp_39 * tmp_439;
      real_t tmp_447 = tmp_41 * tmp_443;
      real_t tmp_448 = 0.031405749086161582 * tmp_46 + 0.031405749086161582 * tmp_47;
      real_t tmp_449 = tmp_19 * ( tmp_44 + tmp_448 );
      real_t tmp_450 = tmp_43 * tmp_449;
      real_t tmp_451 = tmp_449 * tmp_51;
      real_t tmp_452 = tmp_449 * tmp_53;
      real_t tmp_453 = -tmp_440 - tmp_441 - tmp_444 - tmp_445 - tmp_446 - tmp_447 - tmp_450 - tmp_451 - tmp_452 + 1;
      real_t tmp_454 = tmp_75 * ( tmp_438 + tmp_88 );
      real_t tmp_455 = tmp_75 * ( tmp_442 + tmp_90 );
      real_t tmp_456 = tmp_75 * ( tmp_448 + tmp_92 );
      real_t tmp_457 = tmp_98 * ( tmp_66 * ( tmp_454 * tmp_84 + tmp_455 * tmp_81 + tmp_456 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_454 * tmp_85 + tmp_455 * tmp_82 + tmp_456 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_454 * tmp_86 + tmp_455 * tmp_83 + tmp_456 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_458 = 0.0068572537431980923 * tmp_100;
      real_t tmp_459 = 0.19601935860219369 * tmp_22 + 0.19601935860219369 * tmp_23;
      real_t tmp_460 = tmp_19 * ( tmp_20 + tmp_459 );
      real_t tmp_461 = tmp_460 * tmp_8;
      real_t tmp_462 = tmp_27 * tmp_460;
      real_t tmp_463 = 0.19601935860219369 * tmp_32 + 0.19601935860219369 * tmp_33;
      real_t tmp_464 = tmp_19 * ( tmp_30 + tmp_463 );
      real_t tmp_465 = tmp_29 * tmp_464;
      real_t tmp_466 = tmp_37 * tmp_464;
      real_t tmp_467 = tmp_39 * tmp_460;
      real_t tmp_468 = tmp_41 * tmp_464;
      real_t tmp_469 = 0.19601935860219369 * tmp_46 + 0.19601935860219369 * tmp_47;
      real_t tmp_470 = tmp_19 * ( tmp_44 + tmp_469 );
      real_t tmp_471 = tmp_43 * tmp_470;
      real_t tmp_472 = tmp_470 * tmp_51;
      real_t tmp_473 = tmp_470 * tmp_53;
      real_t tmp_474 = -tmp_461 - tmp_462 - tmp_465 - tmp_466 - tmp_467 - tmp_468 - tmp_471 - tmp_472 - tmp_473 + 1;
      real_t tmp_475 = tmp_75 * ( tmp_459 + tmp_88 );
      real_t tmp_476 = tmp_75 * ( tmp_463 + tmp_90 );
      real_t tmp_477 = tmp_75 * ( tmp_469 + tmp_92 );
      real_t tmp_478 = tmp_98 * ( tmp_66 * ( tmp_475 * tmp_84 + tmp_476 * tmp_81 + tmp_477 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_475 * tmp_85 + tmp_476 * tmp_82 + tmp_477 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_475 * tmp_86 + tmp_476 * tmp_83 + tmp_477 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_479 = 0.037198804536718075 * tmp_100;
      real_t tmp_480 = 0.40446199974765351 * tmp_22 + 0.40446199974765351 * tmp_23;
      real_t tmp_481 = tmp_19 * ( tmp_20 + tmp_480 );
      real_t tmp_482 = tmp_481 * tmp_8;
      real_t tmp_483 = tmp_27 * tmp_481;
      real_t tmp_484 = 0.40446199974765351 * tmp_32 + 0.40446199974765351 * tmp_33;
      real_t tmp_485 = tmp_19 * ( tmp_30 + tmp_484 );
      real_t tmp_486 = tmp_29 * tmp_485;
      real_t tmp_487 = tmp_37 * tmp_485;
      real_t tmp_488 = tmp_39 * tmp_481;
      real_t tmp_489 = tmp_41 * tmp_485;
      real_t tmp_490 = 0.40446199974765351 * tmp_46 + 0.40446199974765351 * tmp_47;
      real_t tmp_491 = tmp_19 * ( tmp_44 + tmp_490 );
      real_t tmp_492 = tmp_43 * tmp_491;
      real_t tmp_493 = tmp_491 * tmp_51;
      real_t tmp_494 = tmp_491 * tmp_53;
      real_t tmp_495 = -tmp_482 - tmp_483 - tmp_486 - tmp_487 - tmp_488 - tmp_489 - tmp_492 - tmp_493 - tmp_494 + 1;
      real_t tmp_496 = tmp_75 * ( tmp_480 + tmp_88 );
      real_t tmp_497 = tmp_75 * ( tmp_484 + tmp_90 );
      real_t tmp_498 = tmp_75 * ( tmp_490 + tmp_92 );
      real_t tmp_499 = tmp_98 * ( tmp_66 * ( tmp_496 * tmp_84 + tmp_497 * tmp_81 + tmp_498 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_496 * tmp_85 + tmp_497 * tmp_82 + tmp_498 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_496 * tmp_86 + tmp_497 * tmp_83 + tmp_498 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_500 = 0.042507265838595799 * tmp_100;
      real_t tmp_501 = 0.1711304259088916 * tmp_22 + 0.041227165399737475 * tmp_23;
      real_t tmp_502 = tmp_19 * ( tmp_20 + tmp_501 );
      real_t tmp_503 = tmp_502 * tmp_8;
      real_t tmp_504 = tmp_27 * tmp_502;
      real_t tmp_505 = 0.1711304259088916 * tmp_32 + 0.041227165399737475 * tmp_33;
      real_t tmp_506 = tmp_19 * ( tmp_30 + tmp_505 );
      real_t tmp_507 = tmp_29 * tmp_506;
      real_t tmp_508 = tmp_37 * tmp_506;
      real_t tmp_509 = tmp_39 * tmp_502;
      real_t tmp_510 = tmp_41 * tmp_506;
      real_t tmp_511 = 0.1711304259088916 * tmp_46 + 0.041227165399737475 * tmp_47;
      real_t tmp_512 = tmp_19 * ( tmp_44 + tmp_511 );
      real_t tmp_513 = tmp_43 * tmp_512;
      real_t tmp_514 = tmp_51 * tmp_512;
      real_t tmp_515 = tmp_512 * tmp_53;
      real_t tmp_516 = -tmp_503 - tmp_504 - tmp_507 - tmp_508 - tmp_509 - tmp_510 - tmp_513 - tmp_514 - tmp_515 + 1;
      real_t tmp_517 = tmp_75 * ( tmp_501 + tmp_88 );
      real_t tmp_518 = tmp_75 * ( tmp_505 + tmp_90 );
      real_t tmp_519 = tmp_75 * ( tmp_511 + tmp_92 );
      real_t tmp_520 = tmp_98 * ( tmp_66 * ( tmp_517 * tmp_84 + tmp_518 * tmp_81 + tmp_519 * tmp_64 - 1.0 / 4.0 ) +
                                  tmp_67 * ( tmp_517 * tmp_85 + tmp_518 * tmp_82 + tmp_519 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_70 * ( tmp_517 * tmp_86 + tmp_518 * tmp_83 + tmp_519 * tmp_79 - 1.0 / 4.0 ) );
      real_t tmp_521 = 0.019202922745021479 * tmp_100;
      real_t tmp_522 = tmp_40 + tmp_42 + tmp_54;
      real_t tmp_523 = tmp_110 + tmp_111 + tmp_116;
      real_t tmp_524 = tmp_131 + tmp_132 + tmp_137;
      real_t tmp_525 = tmp_152 + tmp_153 + tmp_158;
      real_t tmp_526 = tmp_173 + tmp_174 + tmp_179;
      real_t tmp_527 = tmp_194 + tmp_195 + tmp_200;
      real_t tmp_528 = tmp_215 + tmp_216 + tmp_221;
      real_t tmp_529 = tmp_236 + tmp_237 + tmp_242;
      real_t tmp_530 = tmp_257 + tmp_258 + tmp_263;
      real_t tmp_531 = tmp_278 + tmp_279 + tmp_284;
      real_t tmp_532 = tmp_299 + tmp_300 + tmp_305;
      real_t tmp_533 = tmp_320 + tmp_321 + tmp_326;
      real_t tmp_534 = tmp_341 + tmp_342 + tmp_347;
      real_t tmp_535 = tmp_362 + tmp_363 + tmp_368;
      real_t tmp_536 = tmp_383 + tmp_384 + tmp_389;
      real_t tmp_537 = tmp_404 + tmp_405 + tmp_410;
      real_t tmp_538 = tmp_425 + tmp_426 + tmp_431;
      real_t tmp_539 = tmp_446 + tmp_447 + tmp_452;
      real_t tmp_540 = tmp_467 + tmp_468 + tmp_473;
      real_t tmp_541 = tmp_488 + tmp_489 + tmp_494;
      real_t tmp_542 = tmp_509 + tmp_510 + tmp_515;
      real_t tmp_543 = tmp_28 + tmp_38 + tmp_52;
      real_t tmp_544 = tmp_105 + tmp_109 + tmp_115;
      real_t tmp_545 = tmp_126 + tmp_130 + tmp_136;
      real_t tmp_546 = tmp_147 + tmp_151 + tmp_157;
      real_t tmp_547 = tmp_168 + tmp_172 + tmp_178;
      real_t tmp_548 = tmp_189 + tmp_193 + tmp_199;
      real_t tmp_549 = tmp_210 + tmp_214 + tmp_220;
      real_t tmp_550 = tmp_231 + tmp_235 + tmp_241;
      real_t tmp_551 = tmp_252 + tmp_256 + tmp_262;
      real_t tmp_552 = tmp_273 + tmp_277 + tmp_283;
      real_t tmp_553 = tmp_294 + tmp_298 + tmp_304;
      real_t tmp_554 = tmp_315 + tmp_319 + tmp_325;
      real_t tmp_555 = tmp_336 + tmp_340 + tmp_346;
      real_t tmp_556 = tmp_357 + tmp_361 + tmp_367;
      real_t tmp_557 = tmp_378 + tmp_382 + tmp_388;
      real_t tmp_558 = tmp_399 + tmp_403 + tmp_409;
      real_t tmp_559 = tmp_420 + tmp_424 + tmp_430;
      real_t tmp_560 = tmp_441 + tmp_445 + tmp_451;
      real_t tmp_561 = tmp_462 + tmp_466 + tmp_472;
      real_t tmp_562 = tmp_483 + tmp_487 + tmp_493;
      real_t tmp_563 = tmp_504 + tmp_508 + tmp_514;
      real_t tmp_564 = tmp_26 + tmp_36 + tmp_50;
      real_t tmp_565 = tmp_104 + tmp_108 + tmp_114;
      real_t tmp_566 = tmp_125 + tmp_129 + tmp_135;
      real_t tmp_567 = tmp_146 + tmp_150 + tmp_156;
      real_t tmp_568 = tmp_167 + tmp_171 + tmp_177;
      real_t tmp_569 = tmp_188 + tmp_192 + tmp_198;
      real_t tmp_570 = tmp_209 + tmp_213 + tmp_219;
      real_t tmp_571 = tmp_230 + tmp_234 + tmp_240;
      real_t tmp_572 = tmp_251 + tmp_255 + tmp_261;
      real_t tmp_573 = tmp_272 + tmp_276 + tmp_282;
      real_t tmp_574 = tmp_293 + tmp_297 + tmp_303;
      real_t tmp_575 = tmp_314 + tmp_318 + tmp_324;
      real_t tmp_576 = tmp_335 + tmp_339 + tmp_345;
      real_t tmp_577 = tmp_356 + tmp_360 + tmp_366;
      real_t tmp_578 = tmp_377 + tmp_381 + tmp_387;
      real_t tmp_579 = tmp_398 + tmp_402 + tmp_408;
      real_t tmp_580 = tmp_419 + tmp_423 + tmp_429;
      real_t tmp_581 = tmp_440 + tmp_444 + tmp_450;
      real_t tmp_582 = tmp_461 + tmp_465 + tmp_471;
      real_t tmp_583 = tmp_482 + tmp_486 + tmp_492;
      real_t tmp_584 = tmp_503 + tmp_507 + tmp_513;
      real_t a_0_0   = tmp_101 * ( -tmp_55 * tmp_87 - tmp_55 * tmp_99 ) + tmp_122 * ( -tmp_117 * tmp_121 - tmp_117 * tmp_87 ) +
                     tmp_143 * ( -tmp_138 * tmp_142 - tmp_138 * tmp_87 ) + tmp_164 * ( -tmp_159 * tmp_163 - tmp_159 * tmp_87 ) +
                     tmp_185 * ( -tmp_180 * tmp_184 - tmp_180 * tmp_87 ) + tmp_206 * ( -tmp_201 * tmp_205 - tmp_201 * tmp_87 ) +
                     tmp_227 * ( -tmp_222 * tmp_226 - tmp_222 * tmp_87 ) + tmp_248 * ( -tmp_243 * tmp_247 - tmp_243 * tmp_87 ) +
                     tmp_269 * ( -tmp_264 * tmp_268 - tmp_264 * tmp_87 ) + tmp_290 * ( -tmp_285 * tmp_289 - tmp_285 * tmp_87 ) +
                     tmp_311 * ( -tmp_306 * tmp_310 - tmp_306 * tmp_87 ) + tmp_332 * ( -tmp_327 * tmp_331 - tmp_327 * tmp_87 ) +
                     tmp_353 * ( -tmp_348 * tmp_352 - tmp_348 * tmp_87 ) + tmp_374 * ( -tmp_369 * tmp_373 - tmp_369 * tmp_87 ) +
                     tmp_395 * ( -tmp_390 * tmp_394 - tmp_390 * tmp_87 ) + tmp_416 * ( -tmp_411 * tmp_415 - tmp_411 * tmp_87 ) +
                     tmp_437 * ( -tmp_432 * tmp_436 - tmp_432 * tmp_87 ) + tmp_458 * ( -tmp_453 * tmp_457 - tmp_453 * tmp_87 ) +
                     tmp_479 * ( -tmp_474 * tmp_478 - tmp_474 * tmp_87 ) + tmp_500 * ( -tmp_495 * tmp_499 - tmp_495 * tmp_87 ) +
                     tmp_521 * ( -tmp_516 * tmp_520 - tmp_516 * tmp_87 );
      real_t a_1_0 = tmp_101 * ( -tmp_522 * tmp_87 - tmp_522 * tmp_99 ) + tmp_122 * ( -tmp_121 * tmp_523 - tmp_523 * tmp_87 ) +
                     tmp_143 * ( -tmp_142 * tmp_524 - tmp_524 * tmp_87 ) + tmp_164 * ( -tmp_163 * tmp_525 - tmp_525 * tmp_87 ) +
                     tmp_185 * ( -tmp_184 * tmp_526 - tmp_526 * tmp_87 ) + tmp_206 * ( -tmp_205 * tmp_527 - tmp_527 * tmp_87 ) +
                     tmp_227 * ( -tmp_226 * tmp_528 - tmp_528 * tmp_87 ) + tmp_248 * ( -tmp_247 * tmp_529 - tmp_529 * tmp_87 ) +
                     tmp_269 * ( -tmp_268 * tmp_530 - tmp_530 * tmp_87 ) + tmp_290 * ( -tmp_289 * tmp_531 - tmp_531 * tmp_87 ) +
                     tmp_311 * ( -tmp_310 * tmp_532 - tmp_532 * tmp_87 ) + tmp_332 * ( -tmp_331 * tmp_533 - tmp_533 * tmp_87 ) +
                     tmp_353 * ( -tmp_352 * tmp_534 - tmp_534 * tmp_87 ) + tmp_374 * ( -tmp_373 * tmp_535 - tmp_535 * tmp_87 ) +
                     tmp_395 * ( -tmp_394 * tmp_536 - tmp_536 * tmp_87 ) + tmp_416 * ( -tmp_415 * tmp_537 - tmp_537 * tmp_87 ) +
                     tmp_437 * ( -tmp_436 * tmp_538 - tmp_538 * tmp_87 ) + tmp_458 * ( -tmp_457 * tmp_539 - tmp_539 * tmp_87 ) +
                     tmp_479 * ( -tmp_478 * tmp_540 - tmp_540 * tmp_87 ) + tmp_500 * ( -tmp_499 * tmp_541 - tmp_541 * tmp_87 ) +
                     tmp_521 * ( -tmp_520 * tmp_542 - tmp_542 * tmp_87 );
      real_t a_2_0 = tmp_101 * ( -tmp_543 * tmp_87 - tmp_543 * tmp_99 ) + tmp_122 * ( -tmp_121 * tmp_544 - tmp_544 * tmp_87 ) +
                     tmp_143 * ( -tmp_142 * tmp_545 - tmp_545 * tmp_87 ) + tmp_164 * ( -tmp_163 * tmp_546 - tmp_546 * tmp_87 ) +
                     tmp_185 * ( -tmp_184 * tmp_547 - tmp_547 * tmp_87 ) + tmp_206 * ( -tmp_205 * tmp_548 - tmp_548 * tmp_87 ) +
                     tmp_227 * ( -tmp_226 * tmp_549 - tmp_549 * tmp_87 ) + tmp_248 * ( -tmp_247 * tmp_550 - tmp_550 * tmp_87 ) +
                     tmp_269 * ( -tmp_268 * tmp_551 - tmp_551 * tmp_87 ) + tmp_290 * ( -tmp_289 * tmp_552 - tmp_552 * tmp_87 ) +
                     tmp_311 * ( -tmp_310 * tmp_553 - tmp_553 * tmp_87 ) + tmp_332 * ( -tmp_331 * tmp_554 - tmp_554 * tmp_87 ) +
                     tmp_353 * ( -tmp_352 * tmp_555 - tmp_555 * tmp_87 ) + tmp_374 * ( -tmp_373 * tmp_556 - tmp_556 * tmp_87 ) +
                     tmp_395 * ( -tmp_394 * tmp_557 - tmp_557 * tmp_87 ) + tmp_416 * ( -tmp_415 * tmp_558 - tmp_558 * tmp_87 ) +
                     tmp_437 * ( -tmp_436 * tmp_559 - tmp_559 * tmp_87 ) + tmp_458 * ( -tmp_457 * tmp_560 - tmp_560 * tmp_87 ) +
                     tmp_479 * ( -tmp_478 * tmp_561 - tmp_561 * tmp_87 ) + tmp_500 * ( -tmp_499 * tmp_562 - tmp_562 * tmp_87 ) +
                     tmp_521 * ( -tmp_520 * tmp_563 - tmp_563 * tmp_87 );
      real_t a_3_0 = tmp_101 * ( -tmp_564 * tmp_87 - tmp_564 * tmp_99 ) + tmp_122 * ( -tmp_121 * tmp_565 - tmp_565 * tmp_87 ) +
                     tmp_143 * ( -tmp_142 * tmp_566 - tmp_566 * tmp_87 ) + tmp_164 * ( -tmp_163 * tmp_567 - tmp_567 * tmp_87 ) +
                     tmp_185 * ( -tmp_184 * tmp_568 - tmp_568 * tmp_87 ) + tmp_206 * ( -tmp_205 * tmp_569 - tmp_569 * tmp_87 ) +
                     tmp_227 * ( -tmp_226 * tmp_570 - tmp_570 * tmp_87 ) + tmp_248 * ( -tmp_247 * tmp_571 - tmp_571 * tmp_87 ) +
                     tmp_269 * ( -tmp_268 * tmp_572 - tmp_572 * tmp_87 ) + tmp_290 * ( -tmp_289 * tmp_573 - tmp_573 * tmp_87 ) +
                     tmp_311 * ( -tmp_310 * tmp_574 - tmp_574 * tmp_87 ) + tmp_332 * ( -tmp_331 * tmp_575 - tmp_575 * tmp_87 ) +
                     tmp_353 * ( -tmp_352 * tmp_576 - tmp_576 * tmp_87 ) + tmp_374 * ( -tmp_373 * tmp_577 - tmp_577 * tmp_87 ) +
                     tmp_395 * ( -tmp_394 * tmp_578 - tmp_578 * tmp_87 ) + tmp_416 * ( -tmp_415 * tmp_579 - tmp_579 * tmp_87 ) +
                     tmp_437 * ( -tmp_436 * tmp_580 - tmp_580 * tmp_87 ) + tmp_458 * ( -tmp_457 * tmp_581 - tmp_581 * tmp_87 ) +
                     tmp_479 * ( -tmp_478 * tmp_582 - tmp_582 * tmp_87 ) + tmp_500 * ( -tmp_499 * tmp_583 - tmp_583 * tmp_87 ) +
                     tmp_521 * ( -tmp_520 * tmp_584 - tmp_584 * tmp_87 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
      elMat( 3, 0 ) = a_3_0;
   }

   void integrateFacetDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                           const std::vector< Point3D >& coordsFacet,
                                           const Point3D&,
                                           const Point3D&     outwardNormal,
                                           const DGBasisInfo& trialBasis,
                                           const DGBasisInfo& testBasis,
                                           int                trialDegree,
                                           int                testDegree,
                                           MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = p_affine_2_0 + tmp_0;
      real_t tmp_5  = p_affine_1_1 + tmp_2;
      real_t tmp_6  = tmp_1 * tmp_3 - tmp_4 * tmp_5;
      real_t tmp_7  = -p_affine_0_2;
      real_t tmp_8  = p_affine_3_2 + tmp_7;
      real_t tmp_9  = tmp_3 * tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_7;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11 * tmp_4;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_13 * tmp_14;
      real_t tmp_16 = tmp_11 * tmp_14;
      real_t tmp_17 = tmp_4 * tmp_8;
      real_t tmp_18 = tmp_13 * tmp_3;
      real_t tmp_19 =
          1.0 / ( -tmp_1 * tmp_16 + tmp_1 * tmp_9 + tmp_10 * tmp_12 - tmp_10 * tmp_18 + tmp_15 * tmp_5 - tmp_17 * tmp_5 );
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19 * ( 0.031405749086161582 * tmp_21 + 0.93718850182767688 * tmp_22 + tmp_23 );
      real_t tmp_25 = tmp_24 * tmp_6;
      real_t tmp_26 = -tmp_1 * tmp_11 + tmp_13 * tmp_5;
      real_t tmp_27 = tmp_24 * tmp_26;
      real_t tmp_28 = -tmp_1 * tmp_14 + tmp_10 * tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19 * ( 0.031405749086161582 * tmp_30 + 0.93718850182767688 * tmp_31 + tmp_32 );
      real_t tmp_34 = tmp_28 * tmp_33;
      real_t tmp_35 = tmp_1 * tmp_8 - tmp_10 * tmp_13;
      real_t tmp_36 = tmp_33 * tmp_35;
      real_t tmp_37 = tmp_12 - tmp_18;
      real_t tmp_38 = tmp_24 * tmp_37;
      real_t tmp_39 = tmp_15 - tmp_17;
      real_t tmp_40 = tmp_33 * tmp_39;
      real_t tmp_41 = -tmp_10 * tmp_3 + tmp_14 * tmp_5;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19 * ( 0.031405749086161582 * tmp_43 + 0.93718850182767688 * tmp_44 + tmp_45 );
      real_t tmp_47 = tmp_41 * tmp_46;
      real_t tmp_48 = tmp_10 * tmp_11 - tmp_5 * tmp_8;
      real_t tmp_49 = tmp_46 * tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46 * tmp_50;
      real_t tmp_52 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_53 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_54 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_55 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_58 =
          1.0 * std::pow( ( std::abs( tmp_22 * tmp_55 - tmp_31 * tmp_57 ) * std::abs( tmp_22 * tmp_55 - tmp_31 * tmp_57 ) ) +
                              ( std::abs( tmp_22 * tmp_56 - tmp_44 * tmp_57 ) * std::abs( tmp_22 * tmp_56 - tmp_44 * tmp_57 ) ) +
                              ( std::abs( tmp_31 * tmp_56 - tmp_44 * tmp_55 ) * std::abs( tmp_31 * tmp_56 - tmp_44 * tmp_55 ) ),
                          0.25 );
      real_t tmp_59 = 0.0068572537431980923 * tmp_58 *
                      ( tmp_1 * ( tmp_52 - 1.0 / 4.0 ) + tmp_13 * ( tmp_54 - 1.0 / 4.0 ) + tmp_4 * ( tmp_53 - 1.0 / 4.0 ) );
      real_t tmp_60 = tmp_19 * ( 0.19601935860219369 * tmp_21 + 0.60796128279561268 * tmp_22 + tmp_23 );
      real_t tmp_61 = tmp_6 * tmp_60;
      real_t tmp_62 = tmp_26 * tmp_60;
      real_t tmp_63 = tmp_19 * ( 0.19601935860219369 * tmp_30 + 0.60796128279561268 * tmp_31 + tmp_32 );
      real_t tmp_64 = tmp_28 * tmp_63;
      real_t tmp_65 = tmp_35 * tmp_63;
      real_t tmp_66 = tmp_37 * tmp_60;
      real_t tmp_67 = tmp_39 * tmp_63;
      real_t tmp_68 = tmp_19 * ( 0.19601935860219369 * tmp_43 + 0.60796128279561268 * tmp_44 + tmp_45 );
      real_t tmp_69 = tmp_41 * tmp_68;
      real_t tmp_70 = tmp_48 * tmp_68;
      real_t tmp_71 = tmp_50 * tmp_68;
      real_t tmp_72 = tmp_66 + tmp_67 + tmp_71;
      real_t tmp_73 = tmp_62 + tmp_65 + tmp_70;
      real_t tmp_74 = tmp_61 + tmp_64 + tmp_69;
      real_t tmp_75 = 0.037198804536718075 * tmp_58 *
                      ( tmp_1 * ( tmp_72 - 1.0 / 4.0 ) + tmp_13 * ( tmp_74 - 1.0 / 4.0 ) + tmp_4 * ( tmp_73 - 1.0 / 4.0 ) );
      real_t tmp_76 = tmp_19 * ( 0.37605877282253791 * tmp_21 + 0.039308471900058539 * tmp_22 + tmp_23 );
      real_t tmp_77 = tmp_6 * tmp_76;
      real_t tmp_78 = tmp_26 * tmp_76;
      real_t tmp_79 = tmp_19 * ( 0.37605877282253791 * tmp_30 + 0.039308471900058539 * tmp_31 + tmp_32 );
      real_t tmp_80 = tmp_28 * tmp_79;
      real_t tmp_81 = tmp_35 * tmp_79;
      real_t tmp_82 = tmp_37 * tmp_76;
      real_t tmp_83 = tmp_39 * tmp_79;
      real_t tmp_84 = tmp_19 * ( 0.37605877282253791 * tmp_43 + 0.039308471900058539 * tmp_44 + tmp_45 );
      real_t tmp_85 = tmp_41 * tmp_84;
      real_t tmp_86 = tmp_48 * tmp_84;
      real_t tmp_87 = tmp_50 * tmp_84;
      real_t tmp_88 = tmp_82 + tmp_83 + tmp_87;
      real_t tmp_89 = tmp_78 + tmp_81 + tmp_86;
      real_t tmp_90 = tmp_77 + tmp_80 + tmp_85;
      real_t tmp_91 = 0.020848748529055869 * tmp_58 *
                      ( tmp_1 * ( tmp_88 - 1.0 / 4.0 ) + tmp_13 * ( tmp_90 - 1.0 / 4.0 ) + tmp_4 * ( tmp_89 - 1.0 / 4.0 ) );
      real_t tmp_92  = tmp_19 * ( 0.78764240869137092 * tmp_21 + 0.1711304259088916 * tmp_22 + tmp_23 );
      real_t tmp_93  = tmp_6 * tmp_92;
      real_t tmp_94  = tmp_26 * tmp_92;
      real_t tmp_95  = tmp_19 * ( 0.78764240869137092 * tmp_30 + 0.1711304259088916 * tmp_31 + tmp_32 );
      real_t tmp_96  = tmp_28 * tmp_95;
      real_t tmp_97  = tmp_35 * tmp_95;
      real_t tmp_98  = tmp_37 * tmp_92;
      real_t tmp_99  = tmp_39 * tmp_95;
      real_t tmp_100 = tmp_19 * ( 0.78764240869137092 * tmp_43 + 0.1711304259088916 * tmp_44 + tmp_45 );
      real_t tmp_101 = tmp_100 * tmp_41;
      real_t tmp_102 = tmp_100 * tmp_48;
      real_t tmp_103 = tmp_100 * tmp_50;
      real_t tmp_104 = tmp_103 + tmp_98 + tmp_99;
      real_t tmp_105 = tmp_102 + tmp_94 + tmp_97;
      real_t tmp_106 = tmp_101 + tmp_93 + tmp_96;
      real_t tmp_107 = 0.019202922745021479 * tmp_58 *
                       ( tmp_1 * ( tmp_104 - 1.0 / 4.0 ) + tmp_13 * ( tmp_106 - 1.0 / 4.0 ) + tmp_4 * ( tmp_105 - 1.0 / 4.0 ) );
      real_t tmp_108 = tmp_19 * ( 0.58463275527740355 * tmp_21 + 0.37605877282253791 * tmp_22 + tmp_23 );
      real_t tmp_109 = tmp_108 * tmp_6;
      real_t tmp_110 = tmp_108 * tmp_26;
      real_t tmp_111 = tmp_19 * ( 0.58463275527740355 * tmp_30 + 0.37605877282253791 * tmp_31 + tmp_32 );
      real_t tmp_112 = tmp_111 * tmp_28;
      real_t tmp_113 = tmp_111 * tmp_35;
      real_t tmp_114 = tmp_108 * tmp_37;
      real_t tmp_115 = tmp_111 * tmp_39;
      real_t tmp_116 = tmp_19 * ( 0.58463275527740355 * tmp_43 + 0.37605877282253791 * tmp_44 + tmp_45 );
      real_t tmp_117 = tmp_116 * tmp_41;
      real_t tmp_118 = tmp_116 * tmp_48;
      real_t tmp_119 = tmp_116 * tmp_50;
      real_t tmp_120 = tmp_114 + tmp_115 + tmp_119;
      real_t tmp_121 = tmp_110 + tmp_113 + tmp_118;
      real_t tmp_122 = tmp_109 + tmp_112 + tmp_117;
      real_t tmp_123 = 0.020848748529055869 * tmp_58 *
                       ( tmp_1 * ( tmp_120 - 1.0 / 4.0 ) + tmp_13 * ( tmp_122 - 1.0 / 4.0 ) + tmp_4 * ( tmp_121 - 1.0 / 4.0 ) );
      real_t tmp_124 = tmp_19 * ( 0.041227165399737475 * tmp_21 + 0.78764240869137092 * tmp_22 + tmp_23 );
      real_t tmp_125 = tmp_124 * tmp_6;
      real_t tmp_126 = tmp_124 * tmp_26;
      real_t tmp_127 = tmp_19 * ( 0.041227165399737475 * tmp_30 + 0.78764240869137092 * tmp_31 + tmp_32 );
      real_t tmp_128 = tmp_127 * tmp_28;
      real_t tmp_129 = tmp_127 * tmp_35;
      real_t tmp_130 = tmp_124 * tmp_37;
      real_t tmp_131 = tmp_127 * tmp_39;
      real_t tmp_132 = tmp_19 * ( 0.041227165399737475 * tmp_43 + 0.78764240869137092 * tmp_44 + tmp_45 );
      real_t tmp_133 = tmp_132 * tmp_41;
      real_t tmp_134 = tmp_132 * tmp_48;
      real_t tmp_135 = tmp_132 * tmp_50;
      real_t tmp_136 = tmp_130 + tmp_131 + tmp_135;
      real_t tmp_137 = tmp_126 + tmp_129 + tmp_134;
      real_t tmp_138 = tmp_125 + tmp_128 + tmp_133;
      real_t tmp_139 = 0.019202922745021479 * tmp_58 *
                       ( tmp_1 * ( tmp_136 - 1.0 / 4.0 ) + tmp_13 * ( tmp_138 - 1.0 / 4.0 ) + tmp_4 * ( tmp_137 - 1.0 / 4.0 ) );
      real_t tmp_140 = tmp_19 * ( 0.039308471900058539 * tmp_21 + 0.58463275527740355 * tmp_22 + tmp_23 );
      real_t tmp_141 = tmp_140 * tmp_6;
      real_t tmp_142 = tmp_140 * tmp_26;
      real_t tmp_143 = tmp_19 * ( 0.039308471900058539 * tmp_30 + 0.58463275527740355 * tmp_31 + tmp_32 );
      real_t tmp_144 = tmp_143 * tmp_28;
      real_t tmp_145 = tmp_143 * tmp_35;
      real_t tmp_146 = tmp_140 * tmp_37;
      real_t tmp_147 = tmp_143 * tmp_39;
      real_t tmp_148 = tmp_19 * ( 0.039308471900058539 * tmp_43 + 0.58463275527740355 * tmp_44 + tmp_45 );
      real_t tmp_149 = tmp_148 * tmp_41;
      real_t tmp_150 = tmp_148 * tmp_48;
      real_t tmp_151 = tmp_148 * tmp_50;
      real_t tmp_152 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_153 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_154 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_155 = 0.020848748529055869 * tmp_58 *
                       ( tmp_1 * ( tmp_152 - 1.0 / 4.0 ) + tmp_13 * ( tmp_154 - 1.0 / 4.0 ) + tmp_4 * ( tmp_153 - 1.0 / 4.0 ) );
      real_t tmp_156 = tmp_19 * ( 0.78764240869137092 * tmp_21 + 0.041227165399737475 * tmp_22 + tmp_23 );
      real_t tmp_157 = tmp_156 * tmp_6;
      real_t tmp_158 = tmp_156 * tmp_26;
      real_t tmp_159 = tmp_19 * ( 0.78764240869137092 * tmp_30 + 0.041227165399737475 * tmp_31 + tmp_32 );
      real_t tmp_160 = tmp_159 * tmp_28;
      real_t tmp_161 = tmp_159 * tmp_35;
      real_t tmp_162 = tmp_156 * tmp_37;
      real_t tmp_163 = tmp_159 * tmp_39;
      real_t tmp_164 = tmp_19 * ( 0.78764240869137092 * tmp_43 + 0.041227165399737475 * tmp_44 + tmp_45 );
      real_t tmp_165 = tmp_164 * tmp_41;
      real_t tmp_166 = tmp_164 * tmp_48;
      real_t tmp_167 = tmp_164 * tmp_50;
      real_t tmp_168 = tmp_162 + tmp_163 + tmp_167;
      real_t tmp_169 = tmp_158 + tmp_161 + tmp_166;
      real_t tmp_170 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_171 = 0.019202922745021479 * tmp_58 *
                       ( tmp_1 * ( tmp_168 - 1.0 / 4.0 ) + tmp_13 * ( tmp_170 - 1.0 / 4.0 ) + tmp_4 * ( tmp_169 - 1.0 / 4.0 ) );
      real_t tmp_172 = tmp_19 * ( 0.58463275527740355 * tmp_21 + 0.039308471900058539 * tmp_22 + tmp_23 );
      real_t tmp_173 = tmp_172 * tmp_6;
      real_t tmp_174 = tmp_172 * tmp_26;
      real_t tmp_175 = tmp_19 * ( 0.58463275527740355 * tmp_30 + 0.039308471900058539 * tmp_31 + tmp_32 );
      real_t tmp_176 = tmp_175 * tmp_28;
      real_t tmp_177 = tmp_175 * tmp_35;
      real_t tmp_178 = tmp_172 * tmp_37;
      real_t tmp_179 = tmp_175 * tmp_39;
      real_t tmp_180 = tmp_19 * ( 0.58463275527740355 * tmp_43 + 0.039308471900058539 * tmp_44 + tmp_45 );
      real_t tmp_181 = tmp_180 * tmp_41;
      real_t tmp_182 = tmp_180 * tmp_48;
      real_t tmp_183 = tmp_180 * tmp_50;
      real_t tmp_184 = tmp_178 + tmp_179 + tmp_183;
      real_t tmp_185 = tmp_174 + tmp_177 + tmp_182;
      real_t tmp_186 = tmp_173 + tmp_176 + tmp_181;
      real_t tmp_187 = 0.020848748529055869 * tmp_58 *
                       ( tmp_1 * ( tmp_184 - 1.0 / 4.0 ) + tmp_13 * ( tmp_186 - 1.0 / 4.0 ) + tmp_4 * ( tmp_185 - 1.0 / 4.0 ) );
      real_t tmp_188 = tmp_19 * ( 0.1711304259088916 * tmp_21 + 0.78764240869137092 * tmp_22 + tmp_23 );
      real_t tmp_189 = tmp_188 * tmp_6;
      real_t tmp_190 = tmp_188 * tmp_26;
      real_t tmp_191 = tmp_19 * ( 0.1711304259088916 * tmp_30 + 0.78764240869137092 * tmp_31 + tmp_32 );
      real_t tmp_192 = tmp_191 * tmp_28;
      real_t tmp_193 = tmp_191 * tmp_35;
      real_t tmp_194 = tmp_188 * tmp_37;
      real_t tmp_195 = tmp_191 * tmp_39;
      real_t tmp_196 = tmp_19 * ( 0.1711304259088916 * tmp_43 + 0.78764240869137092 * tmp_44 + tmp_45 );
      real_t tmp_197 = tmp_196 * tmp_41;
      real_t tmp_198 = tmp_196 * tmp_48;
      real_t tmp_199 = tmp_196 * tmp_50;
      real_t tmp_200 = tmp_194 + tmp_195 + tmp_199;
      real_t tmp_201 = tmp_190 + tmp_193 + tmp_198;
      real_t tmp_202 = tmp_189 + tmp_192 + tmp_197;
      real_t tmp_203 = 0.019202922745021479 * tmp_58 *
                       ( tmp_1 * ( tmp_200 - 1.0 / 4.0 ) + tmp_13 * ( tmp_202 - 1.0 / 4.0 ) + tmp_4 * ( tmp_201 - 1.0 / 4.0 ) );
      real_t tmp_204 = tmp_19 * ( 0.37605877282253791 * tmp_21 + 0.58463275527740355 * tmp_22 + tmp_23 );
      real_t tmp_205 = tmp_204 * tmp_6;
      real_t tmp_206 = tmp_204 * tmp_26;
      real_t tmp_207 = tmp_19 * ( 0.37605877282253791 * tmp_30 + 0.58463275527740355 * tmp_31 + tmp_32 );
      real_t tmp_208 = tmp_207 * tmp_28;
      real_t tmp_209 = tmp_207 * tmp_35;
      real_t tmp_210 = tmp_204 * tmp_37;
      real_t tmp_211 = tmp_207 * tmp_39;
      real_t tmp_212 = tmp_19 * ( 0.37605877282253791 * tmp_43 + 0.58463275527740355 * tmp_44 + tmp_45 );
      real_t tmp_213 = tmp_212 * tmp_41;
      real_t tmp_214 = tmp_212 * tmp_48;
      real_t tmp_215 = tmp_212 * tmp_50;
      real_t tmp_216 = tmp_210 + tmp_211 + tmp_215;
      real_t tmp_217 = tmp_206 + tmp_209 + tmp_214;
      real_t tmp_218 = tmp_205 + tmp_208 + tmp_213;
      real_t tmp_219 = 0.020848748529055869 * tmp_58 *
                       ( tmp_1 * ( tmp_216 - 1.0 / 4.0 ) + tmp_13 * ( tmp_218 - 1.0 / 4.0 ) + tmp_4 * ( tmp_217 - 1.0 / 4.0 ) );
      real_t tmp_220 = tmp_19 * ( 0.041227165399737475 * tmp_21 + 0.1711304259088916 * tmp_22 + tmp_23 );
      real_t tmp_221 = tmp_220 * tmp_6;
      real_t tmp_222 = tmp_220 * tmp_26;
      real_t tmp_223 = tmp_19 * ( 0.041227165399737475 * tmp_30 + 0.1711304259088916 * tmp_31 + tmp_32 );
      real_t tmp_224 = tmp_223 * tmp_28;
      real_t tmp_225 = tmp_223 * tmp_35;
      real_t tmp_226 = tmp_220 * tmp_37;
      real_t tmp_227 = tmp_223 * tmp_39;
      real_t tmp_228 = tmp_19 * ( 0.041227165399737475 * tmp_43 + 0.1711304259088916 * tmp_44 + tmp_45 );
      real_t tmp_229 = tmp_228 * tmp_41;
      real_t tmp_230 = tmp_228 * tmp_48;
      real_t tmp_231 = tmp_228 * tmp_50;
      real_t tmp_232 = tmp_226 + tmp_227 + tmp_231;
      real_t tmp_233 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_234 = tmp_221 + tmp_224 + tmp_229;
      real_t tmp_235 = 0.019202922745021479 * tmp_58 *
                       ( tmp_1 * ( tmp_232 - 1.0 / 4.0 ) + tmp_13 * ( tmp_234 - 1.0 / 4.0 ) + tmp_4 * ( tmp_233 - 1.0 / 4.0 ) );
      real_t tmp_236 = tmp_19 * ( 0.40446199974765351 * tmp_21 + 0.19107600050469298 * tmp_22 + tmp_23 );
      real_t tmp_237 = tmp_236 * tmp_6;
      real_t tmp_238 = tmp_236 * tmp_26;
      real_t tmp_239 = tmp_19 * ( 0.40446199974765351 * tmp_30 + 0.19107600050469298 * tmp_31 + tmp_32 );
      real_t tmp_240 = tmp_239 * tmp_28;
      real_t tmp_241 = tmp_239 * tmp_35;
      real_t tmp_242 = tmp_236 * tmp_37;
      real_t tmp_243 = tmp_239 * tmp_39;
      real_t tmp_244 = tmp_19 * ( 0.40446199974765351 * tmp_43 + 0.19107600050469298 * tmp_44 + tmp_45 );
      real_t tmp_245 = tmp_244 * tmp_41;
      real_t tmp_246 = tmp_244 * tmp_48;
      real_t tmp_247 = tmp_244 * tmp_50;
      real_t tmp_248 = tmp_242 + tmp_243 + tmp_247;
      real_t tmp_249 = tmp_238 + tmp_241 + tmp_246;
      real_t tmp_250 = tmp_237 + tmp_240 + tmp_245;
      real_t tmp_251 = 0.042507265838595799 * tmp_58 *
                       ( tmp_1 * ( tmp_248 - 1.0 / 4.0 ) + tmp_13 * ( tmp_250 - 1.0 / 4.0 ) + tmp_4 * ( tmp_249 - 1.0 / 4.0 ) );
      real_t tmp_252 = tmp_19 * ( 0.039308471900058539 * tmp_21 + 0.37605877282253791 * tmp_22 + tmp_23 );
      real_t tmp_253 = tmp_252 * tmp_6;
      real_t tmp_254 = tmp_252 * tmp_26;
      real_t tmp_255 = tmp_19 * ( 0.039308471900058539 * tmp_30 + 0.37605877282253791 * tmp_31 + tmp_32 );
      real_t tmp_256 = tmp_255 * tmp_28;
      real_t tmp_257 = tmp_255 * tmp_35;
      real_t tmp_258 = tmp_252 * tmp_37;
      real_t tmp_259 = tmp_255 * tmp_39;
      real_t tmp_260 = tmp_19 * ( 0.039308471900058539 * tmp_43 + 0.37605877282253791 * tmp_44 + tmp_45 );
      real_t tmp_261 = tmp_260 * tmp_41;
      real_t tmp_262 = tmp_260 * tmp_48;
      real_t tmp_263 = tmp_260 * tmp_50;
      real_t tmp_264 = tmp_258 + tmp_259 + tmp_263;
      real_t tmp_265 = tmp_254 + tmp_257 + tmp_262;
      real_t tmp_266 = tmp_253 + tmp_256 + tmp_261;
      real_t tmp_267 = 0.020848748529055869 * tmp_58 *
                       ( tmp_1 * ( tmp_264 - 1.0 / 4.0 ) + tmp_13 * ( tmp_266 - 1.0 / 4.0 ) + tmp_4 * ( tmp_265 - 1.0 / 4.0 ) );
      real_t tmp_268 = tmp_19 * ( 0.93718850182767688 * tmp_21 + 0.031405749086161582 * tmp_22 + tmp_23 );
      real_t tmp_269 = tmp_268 * tmp_6;
      real_t tmp_270 = tmp_26 * tmp_268;
      real_t tmp_271 = tmp_19 * ( 0.93718850182767688 * tmp_30 + 0.031405749086161582 * tmp_31 + tmp_32 );
      real_t tmp_272 = tmp_271 * tmp_28;
      real_t tmp_273 = tmp_271 * tmp_35;
      real_t tmp_274 = tmp_268 * tmp_37;
      real_t tmp_275 = tmp_271 * tmp_39;
      real_t tmp_276 = tmp_19 * ( 0.93718850182767688 * tmp_43 + 0.031405749086161582 * tmp_44 + tmp_45 );
      real_t tmp_277 = tmp_276 * tmp_41;
      real_t tmp_278 = tmp_276 * tmp_48;
      real_t tmp_279 = tmp_276 * tmp_50;
      real_t tmp_280 = tmp_274 + tmp_275 + tmp_279;
      real_t tmp_281 = tmp_270 + tmp_273 + tmp_278;
      real_t tmp_282 = tmp_269 + tmp_272 + tmp_277;
      real_t tmp_283 = 0.0068572537431980923 * tmp_58 *
                       ( tmp_1 * ( tmp_280 - 1.0 / 4.0 ) + tmp_13 * ( tmp_282 - 1.0 / 4.0 ) + tmp_4 * ( tmp_281 - 1.0 / 4.0 ) );
      real_t tmp_284 = tmp_19 * ( 0.60796128279561268 * tmp_21 + 0.19601935860219369 * tmp_22 + tmp_23 );
      real_t tmp_285 = tmp_284 * tmp_6;
      real_t tmp_286 = tmp_26 * tmp_284;
      real_t tmp_287 = tmp_19 * ( 0.60796128279561268 * tmp_30 + 0.19601935860219369 * tmp_31 + tmp_32 );
      real_t tmp_288 = tmp_28 * tmp_287;
      real_t tmp_289 = tmp_287 * tmp_35;
      real_t tmp_290 = tmp_284 * tmp_37;
      real_t tmp_291 = tmp_287 * tmp_39;
      real_t tmp_292 = tmp_19 * ( 0.60796128279561268 * tmp_43 + 0.19601935860219369 * tmp_44 + tmp_45 );
      real_t tmp_293 = tmp_292 * tmp_41;
      real_t tmp_294 = tmp_292 * tmp_48;
      real_t tmp_295 = tmp_292 * tmp_50;
      real_t tmp_296 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_297 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_298 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_299 = 0.037198804536718075 * tmp_58 *
                       ( tmp_1 * ( tmp_296 - 1.0 / 4.0 ) + tmp_13 * ( tmp_298 - 1.0 / 4.0 ) + tmp_4 * ( tmp_297 - 1.0 / 4.0 ) );
      real_t tmp_300 = tmp_19 * ( 0.19107600050469298 * tmp_21 + 0.40446199974765351 * tmp_22 + tmp_23 );
      real_t tmp_301 = tmp_300 * tmp_6;
      real_t tmp_302 = tmp_26 * tmp_300;
      real_t tmp_303 = tmp_19 * ( 0.19107600050469298 * tmp_30 + 0.40446199974765351 * tmp_31 + tmp_32 );
      real_t tmp_304 = tmp_28 * tmp_303;
      real_t tmp_305 = tmp_303 * tmp_35;
      real_t tmp_306 = tmp_300 * tmp_37;
      real_t tmp_307 = tmp_303 * tmp_39;
      real_t tmp_308 = tmp_19 * ( 0.19107600050469298 * tmp_43 + 0.40446199974765351 * tmp_44 + tmp_45 );
      real_t tmp_309 = tmp_308 * tmp_41;
      real_t tmp_310 = tmp_308 * tmp_48;
      real_t tmp_311 = tmp_308 * tmp_50;
      real_t tmp_312 = tmp_306 + tmp_307 + tmp_311;
      real_t tmp_313 = tmp_302 + tmp_305 + tmp_310;
      real_t tmp_314 = tmp_301 + tmp_304 + tmp_309;
      real_t tmp_315 = 0.042507265838595799 * tmp_58 *
                       ( tmp_1 * ( tmp_312 - 1.0 / 4.0 ) + tmp_13 * ( tmp_314 - 1.0 / 4.0 ) + tmp_4 * ( tmp_313 - 1.0 / 4.0 ) );
      real_t tmp_316 = tmp_19 * ( 0.031405749086161582 * tmp_21 + 0.031405749086161582 * tmp_22 + tmp_23 );
      real_t tmp_317 = tmp_316 * tmp_6;
      real_t tmp_318 = tmp_26 * tmp_316;
      real_t tmp_319 = tmp_19 * ( 0.031405749086161582 * tmp_30 + 0.031405749086161582 * tmp_31 + tmp_32 );
      real_t tmp_320 = tmp_28 * tmp_319;
      real_t tmp_321 = tmp_319 * tmp_35;
      real_t tmp_322 = tmp_316 * tmp_37;
      real_t tmp_323 = tmp_319 * tmp_39;
      real_t tmp_324 = tmp_19 * ( 0.031405749086161582 * tmp_43 + 0.031405749086161582 * tmp_44 + tmp_45 );
      real_t tmp_325 = tmp_324 * tmp_41;
      real_t tmp_326 = tmp_324 * tmp_48;
      real_t tmp_327 = tmp_324 * tmp_50;
      real_t tmp_328 = tmp_322 + tmp_323 + tmp_327;
      real_t tmp_329 = tmp_318 + tmp_321 + tmp_326;
      real_t tmp_330 = tmp_317 + tmp_320 + tmp_325;
      real_t tmp_331 = 0.0068572537431980923 * tmp_58 *
                       ( tmp_1 * ( tmp_328 - 1.0 / 4.0 ) + tmp_13 * ( tmp_330 - 1.0 / 4.0 ) + tmp_4 * ( tmp_329 - 1.0 / 4.0 ) );
      real_t tmp_332 = tmp_19 * ( 0.19601935860219369 * tmp_21 + 0.19601935860219369 * tmp_22 + tmp_23 );
      real_t tmp_333 = tmp_332 * tmp_6;
      real_t tmp_334 = tmp_26 * tmp_332;
      real_t tmp_335 = tmp_19 * ( 0.19601935860219369 * tmp_30 + 0.19601935860219369 * tmp_31 + tmp_32 );
      real_t tmp_336 = tmp_28 * tmp_335;
      real_t tmp_337 = tmp_335 * tmp_35;
      real_t tmp_338 = tmp_332 * tmp_37;
      real_t tmp_339 = tmp_335 * tmp_39;
      real_t tmp_340 = tmp_19 * ( 0.19601935860219369 * tmp_43 + 0.19601935860219369 * tmp_44 + tmp_45 );
      real_t tmp_341 = tmp_340 * tmp_41;
      real_t tmp_342 = tmp_340 * tmp_48;
      real_t tmp_343 = tmp_340 * tmp_50;
      real_t tmp_344 = tmp_338 + tmp_339 + tmp_343;
      real_t tmp_345 = tmp_334 + tmp_337 + tmp_342;
      real_t tmp_346 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_347 = 0.037198804536718075 * tmp_58 *
                       ( tmp_1 * ( tmp_344 - 1.0 / 4.0 ) + tmp_13 * ( tmp_346 - 1.0 / 4.0 ) + tmp_4 * ( tmp_345 - 1.0 / 4.0 ) );
      real_t tmp_348 = tmp_19 * ( 0.40446199974765351 * tmp_21 + 0.40446199974765351 * tmp_22 + tmp_23 );
      real_t tmp_349 = tmp_348 * tmp_6;
      real_t tmp_350 = tmp_26 * tmp_348;
      real_t tmp_351 = tmp_19 * ( 0.40446199974765351 * tmp_30 + 0.40446199974765351 * tmp_31 + tmp_32 );
      real_t tmp_352 = tmp_28 * tmp_351;
      real_t tmp_353 = tmp_35 * tmp_351;
      real_t tmp_354 = tmp_348 * tmp_37;
      real_t tmp_355 = tmp_351 * tmp_39;
      real_t tmp_356 = tmp_19 * ( 0.40446199974765351 * tmp_43 + 0.40446199974765351 * tmp_44 + tmp_45 );
      real_t tmp_357 = tmp_356 * tmp_41;
      real_t tmp_358 = tmp_356 * tmp_48;
      real_t tmp_359 = tmp_356 * tmp_50;
      real_t tmp_360 = tmp_354 + tmp_355 + tmp_359;
      real_t tmp_361 = tmp_350 + tmp_353 + tmp_358;
      real_t tmp_362 = tmp_349 + tmp_352 + tmp_357;
      real_t tmp_363 = 0.042507265838595799 * tmp_58 *
                       ( tmp_1 * ( tmp_360 - 1.0 / 4.0 ) + tmp_13 * ( tmp_362 - 1.0 / 4.0 ) + tmp_4 * ( tmp_361 - 1.0 / 4.0 ) );
      real_t tmp_364 = tmp_19 * ( 0.1711304259088916 * tmp_21 + 0.041227165399737475 * tmp_22 + tmp_23 );
      real_t tmp_365 = tmp_364 * tmp_6;
      real_t tmp_366 = tmp_26 * tmp_364;
      real_t tmp_367 = tmp_19 * ( 0.1711304259088916 * tmp_30 + 0.041227165399737475 * tmp_31 + tmp_32 );
      real_t tmp_368 = tmp_28 * tmp_367;
      real_t tmp_369 = tmp_35 * tmp_367;
      real_t tmp_370 = tmp_364 * tmp_37;
      real_t tmp_371 = tmp_367 * tmp_39;
      real_t tmp_372 = tmp_19 * ( 0.1711304259088916 * tmp_43 + 0.041227165399737475 * tmp_44 + tmp_45 );
      real_t tmp_373 = tmp_372 * tmp_41;
      real_t tmp_374 = tmp_372 * tmp_48;
      real_t tmp_375 = tmp_372 * tmp_50;
      real_t tmp_376 = tmp_370 + tmp_371 + tmp_375;
      real_t tmp_377 = tmp_366 + tmp_369 + tmp_374;
      real_t tmp_378 = tmp_365 + tmp_368 + tmp_373;
      real_t tmp_379 = 0.019202922745021479 * tmp_58 *
                       ( tmp_1 * ( tmp_376 - 1.0 / 4.0 ) + tmp_13 * ( tmp_378 - 1.0 / 4.0 ) + tmp_4 * ( tmp_377 - 1.0 / 4.0 ) );
      real_t a_0_0 = tmp_107 * ( -tmp_101 - tmp_102 - tmp_103 - tmp_93 - tmp_94 - tmp_96 - tmp_97 - tmp_98 - tmp_99 + 1 ) +
                     tmp_123 * ( -tmp_109 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_115 - tmp_117 - tmp_118 - tmp_119 + 1 ) +
                     tmp_139 * ( -tmp_125 - tmp_126 - tmp_128 - tmp_129 - tmp_130 - tmp_131 - tmp_133 - tmp_134 - tmp_135 + 1 ) +
                     tmp_155 * ( -tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1 ) +
                     tmp_171 * ( -tmp_157 - tmp_158 - tmp_160 - tmp_161 - tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 + 1 ) +
                     tmp_187 * ( -tmp_173 - tmp_174 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 + 1 ) +
                     tmp_203 * ( -tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_197 - tmp_198 - tmp_199 + 1 ) +
                     tmp_219 * ( -tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 - tmp_211 - tmp_213 - tmp_214 - tmp_215 + 1 ) +
                     tmp_235 * ( -tmp_221 - tmp_222 - tmp_224 - tmp_225 - tmp_226 - tmp_227 - tmp_229 - tmp_230 - tmp_231 + 1 ) +
                     tmp_251 * ( -tmp_237 - tmp_238 - tmp_240 - tmp_241 - tmp_242 - tmp_243 - tmp_245 - tmp_246 - tmp_247 + 1 ) +
                     tmp_267 * ( -tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1 ) +
                     tmp_283 * ( -tmp_269 - tmp_270 - tmp_272 - tmp_273 - tmp_274 - tmp_275 - tmp_277 - tmp_278 - tmp_279 + 1 ) +
                     tmp_299 * ( -tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1 ) +
                     tmp_315 * ( -tmp_301 - tmp_302 - tmp_304 - tmp_305 - tmp_306 - tmp_307 - tmp_309 - tmp_310 - tmp_311 + 1 ) +
                     tmp_331 * ( -tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1 ) +
                     tmp_347 * ( -tmp_333 - tmp_334 - tmp_336 - tmp_337 - tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 + 1 ) +
                     tmp_363 * ( -tmp_349 - tmp_350 - tmp_352 - tmp_353 - tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 + 1 ) +
                     tmp_379 * ( -tmp_365 - tmp_366 - tmp_368 - tmp_369 - tmp_370 - tmp_371 - tmp_373 - tmp_374 - tmp_375 + 1 ) +
                     tmp_59 * ( -tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1 ) +
                     tmp_75 * ( -tmp_61 - tmp_62 - tmp_64 - tmp_65 - tmp_66 - tmp_67 - tmp_69 - tmp_70 - tmp_71 + 1 ) +
                     tmp_91 * ( -tmp_77 - tmp_78 - tmp_80 - tmp_81 - tmp_82 - tmp_83 - tmp_85 - tmp_86 - tmp_87 + 1 );
      real_t a_1_0 = tmp_104 * tmp_107 + tmp_120 * tmp_123 + tmp_136 * tmp_139 + tmp_152 * tmp_155 + tmp_168 * tmp_171 +
                     tmp_184 * tmp_187 + tmp_200 * tmp_203 + tmp_216 * tmp_219 + tmp_232 * tmp_235 + tmp_248 * tmp_251 +
                     tmp_264 * tmp_267 + tmp_280 * tmp_283 + tmp_296 * tmp_299 + tmp_312 * tmp_315 + tmp_328 * tmp_331 +
                     tmp_344 * tmp_347 + tmp_360 * tmp_363 + tmp_376 * tmp_379 + tmp_52 * tmp_59 + tmp_72 * tmp_75 +
                     tmp_88 * tmp_91;
      real_t a_2_0 = tmp_105 * tmp_107 + tmp_121 * tmp_123 + tmp_137 * tmp_139 + tmp_153 * tmp_155 + tmp_169 * tmp_171 +
                     tmp_185 * tmp_187 + tmp_201 * tmp_203 + tmp_217 * tmp_219 + tmp_233 * tmp_235 + tmp_249 * tmp_251 +
                     tmp_265 * tmp_267 + tmp_281 * tmp_283 + tmp_297 * tmp_299 + tmp_313 * tmp_315 + tmp_329 * tmp_331 +
                     tmp_345 * tmp_347 + tmp_361 * tmp_363 + tmp_377 * tmp_379 + tmp_53 * tmp_59 + tmp_73 * tmp_75 +
                     tmp_89 * tmp_91;
      real_t a_3_0 = tmp_106 * tmp_107 + tmp_122 * tmp_123 + tmp_138 * tmp_139 + tmp_154 * tmp_155 + tmp_170 * tmp_171 +
                     tmp_186 * tmp_187 + tmp_202 * tmp_203 + tmp_218 * tmp_219 + tmp_234 * tmp_235 + tmp_250 * tmp_251 +
                     tmp_266 * tmp_267 + tmp_282 * tmp_283 + tmp_298 * tmp_299 + tmp_314 * tmp_315 + tmp_330 * tmp_331 +
                     tmp_346 * tmp_347 + tmp_362 * tmp_363 + tmp_378 * tmp_379 + tmp_54 * tmp_59 + tmp_74 * tmp_75 +
                     tmp_90 * tmp_91;
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
      elMat( 3, 0 ) = a_3_0;
   }
};

class EGIIPGVectorLaplaceFormP1E_1 : public hyteg::dg::DGForm
{
 protected:
   void integrateVolume2D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override
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

   virtual void integrateFacetInner2D( const std::vector< Point3D >& coordsElement,
                                       const std::vector< Point3D >& coordsFacet,
                                       const Point3D&                oppositeVertex,
                                       const Point3D&                outwardNormal,
                                       const DGBasisInfo&            trialBasis,
                                       const DGBasisInfo&            testBasis,
                                       int                           trialDegree,
                                       int                           testDegree,
                                       MatrixXr&                     elMat ) const override
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
      real_t tmp_7  = tmp_5 * tmp_6;
      real_t tmp_8  = p_affine_1_1 + tmp_1;
      real_t tmp_9  = 1.0 / ( tmp_7 - tmp_8 * ( p_affine_2_0 + tmp_4 ) );
      real_t tmp_10 = tmp_5 * tmp_9;
      real_t tmp_11 = tmp_10 * tmp_3;
      real_t tmp_12 = tmp_9 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_13 = tmp_12 * tmp_3;
      real_t tmp_14 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_15 = p_affine_6_0 + tmp_4;
      real_t tmp_16 = 0.046910077030668018 * tmp_14 + tmp_15;
      real_t tmp_17 = tmp_6 * tmp_9;
      real_t tmp_18 = tmp_16 * tmp_17;
      real_t tmp_19 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_20 = tmp_19 * tmp_9;
      real_t tmp_21 = tmp_16 * tmp_20;
      real_t tmp_22 = -tmp_11 - tmp_13 - tmp_18 - tmp_21 + 1;
      real_t tmp_23 =
          0.5 * p_affine_10_0 * ( tmp_17 * tmp_19 + tmp_17 * tmp_8 ) + 0.5 * p_affine_10_1 * ( tmp_12 * tmp_8 + tmp_7 * tmp_9 );
      real_t tmp_24 = std::abs( std::pow( ( tmp_0 * tmp_0 ) + ( tmp_14 * tmp_14 ), 1.0 / 2.0 ) );
      real_t tmp_25 = 1.0 / ( tmp_24 );
      real_t tmp_26 = tmp_13 + tmp_18;
      real_t tmp_27 = tmp_11 + tmp_21;
      real_t tmp_28 = tmp_25 * ( tmp_6 * ( tmp_27 - 1.0 / 3.0 ) + tmp_8 * ( tmp_26 - 1.0 / 3.0 ) );
      real_t tmp_29 = 0.11846344252809471 * tmp_24;
      real_t tmp_30 = 0.23076534494715845 * tmp_0 + tmp_2;
      real_t tmp_31 = tmp_10 * tmp_30;
      real_t tmp_32 = tmp_12 * tmp_30;
      real_t tmp_33 = 0.23076534494715845 * tmp_14 + tmp_15;
      real_t tmp_34 = tmp_17 * tmp_33;
      real_t tmp_35 = tmp_20 * tmp_33;
      real_t tmp_36 = -tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1;
      real_t tmp_37 = tmp_32 + tmp_34;
      real_t tmp_38 = tmp_31 + tmp_35;
      real_t tmp_39 = tmp_25 * ( tmp_6 * ( tmp_38 - 1.0 / 3.0 ) + tmp_8 * ( tmp_37 - 1.0 / 3.0 ) );
      real_t tmp_40 = 0.2393143352496831 * tmp_24;
      real_t tmp_41 = 0.5 * tmp_0 + tmp_2;
      real_t tmp_42 = tmp_10 * tmp_41;
      real_t tmp_43 = tmp_12 * tmp_41;
      real_t tmp_44 = 0.5 * tmp_14 + tmp_15;
      real_t tmp_45 = tmp_17 * tmp_44;
      real_t tmp_46 = tmp_20 * tmp_44;
      real_t tmp_47 = -tmp_42 - tmp_43 - tmp_45 - tmp_46 + 1;
      real_t tmp_48 = tmp_43 + tmp_45;
      real_t tmp_49 = tmp_42 + tmp_46;
      real_t tmp_50 = tmp_25 * ( tmp_6 * ( tmp_49 - 1.0 / 3.0 ) + tmp_8 * ( tmp_48 - 1.0 / 3.0 ) );
      real_t tmp_51 = 0.2844444444444445 * tmp_24;
      real_t tmp_52 = 0.7692346550528415 * tmp_0 + tmp_2;
      real_t tmp_53 = tmp_10 * tmp_52;
      real_t tmp_54 = tmp_12 * tmp_52;
      real_t tmp_55 = 0.7692346550528415 * tmp_14 + tmp_15;
      real_t tmp_56 = tmp_17 * tmp_55;
      real_t tmp_57 = tmp_20 * tmp_55;
      real_t tmp_58 = -tmp_53 - tmp_54 - tmp_56 - tmp_57 + 1;
      real_t tmp_59 = tmp_54 + tmp_56;
      real_t tmp_60 = tmp_53 + tmp_57;
      real_t tmp_61 = tmp_25 * ( tmp_6 * ( tmp_60 - 1.0 / 3.0 ) + tmp_8 * ( tmp_59 - 1.0 / 3.0 ) );
      real_t tmp_62 = 0.2393143352496831 * tmp_24;
      real_t tmp_63 = 0.95308992296933193 * tmp_0 + tmp_2;
      real_t tmp_64 = tmp_10 * tmp_63;
      real_t tmp_65 = tmp_12 * tmp_63;
      real_t tmp_66 = 0.95308992296933193 * tmp_14 + tmp_15;
      real_t tmp_67 = tmp_17 * tmp_66;
      real_t tmp_68 = tmp_20 * tmp_66;
      real_t tmp_69 = -tmp_64 - tmp_65 - tmp_67 - tmp_68 + 1;
      real_t tmp_70 = tmp_65 + tmp_67;
      real_t tmp_71 = tmp_64 + tmp_68;
      real_t tmp_72 = tmp_25 * ( tmp_6 * ( tmp_71 - 1.0 / 3.0 ) + tmp_8 * ( tmp_70 - 1.0 / 3.0 ) );
      real_t tmp_73 = 0.11846344252809471 * tmp_24;
      real_t a_0_0  = tmp_29 * ( -tmp_22 * tmp_23 + tmp_22 * tmp_28 ) + tmp_40 * ( -tmp_23 * tmp_36 + tmp_36 * tmp_39 ) +
                     tmp_51 * ( -tmp_23 * tmp_47 + tmp_47 * tmp_50 ) + tmp_62 * ( -tmp_23 * tmp_58 + tmp_58 * tmp_61 ) +
                     tmp_73 * ( -tmp_23 * tmp_69 + tmp_69 * tmp_72 );
      real_t a_1_0 = tmp_29 * ( -tmp_23 * tmp_26 + tmp_26 * tmp_28 ) + tmp_40 * ( -tmp_23 * tmp_37 + tmp_37 * tmp_39 ) +
                     tmp_51 * ( -tmp_23 * tmp_48 + tmp_48 * tmp_50 ) + tmp_62 * ( -tmp_23 * tmp_59 + tmp_59 * tmp_61 ) +
                     tmp_73 * ( -tmp_23 * tmp_70 + tmp_70 * tmp_72 );
      real_t a_2_0 = tmp_29 * ( -tmp_23 * tmp_27 + tmp_27 * tmp_28 ) + tmp_40 * ( -tmp_23 * tmp_38 + tmp_38 * tmp_39 ) +
                     tmp_51 * ( -tmp_23 * tmp_49 + tmp_49 * tmp_50 ) + tmp_62 * ( -tmp_23 * tmp_60 + tmp_60 * tmp_61 ) +
                     tmp_73 * ( -tmp_23 * tmp_71 + tmp_71 * tmp_72 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >& coordsElementInner,
                                          const std::vector< Point3D >& coordsElementOuter,
                                          const std::vector< Point3D >& coordsFacet,
                                          const Point3D&                oppositeVertexInnerElement,
                                          const Point3D&                oppositeVertexOuterElement,
                                          const Point3D&                outwardNormal,
                                          const DGBasisInfo&            trialBasis,
                                          const DGBasisInfo&            testBasis,
                                          int                           trialDegree,
                                          int                           testDegree,
                                          MatrixXr&                     elMat ) const override
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
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + 0.046910077030668018 * tmp_5;
      real_t tmp_7  = tmp_4 * ( tmp_2 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + 0.046910077030668018 * tmp_11;
      real_t tmp_13 = tmp_4 * ( tmp_0 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = -p_affine_3_1;
      real_t tmp_19 = p_affine_4_1 + tmp_18;
      real_t tmp_20 = p_affine_5_1 + tmp_18;
      real_t tmp_21 = -p_affine_3_0;
      real_t tmp_22 = p_affine_4_0 + tmp_21;
      real_t tmp_23 = tmp_20 * tmp_22;
      real_t tmp_24 = 1.0 / ( -tmp_19 * ( p_affine_5_0 + tmp_21 ) + tmp_23 );
      real_t tmp_25 = tmp_20 * tmp_24;
      real_t tmp_26 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_27 = tmp_24 * ( p_affine_3_0 - p_affine_5_0 );
      real_t tmp_28 = 0.5 * p_affine_10_0 * ( tmp_19 * tmp_25 + tmp_25 * tmp_26 ) +
                      0.5 * p_affine_10_1 * ( tmp_19 * tmp_27 + tmp_23 * tmp_24 );
      real_t tmp_29 = std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_30 = 1.0 / ( tmp_29 );
      real_t tmp_31 = tmp_18 + tmp_6;
      real_t tmp_32 = tmp_12 + tmp_21;
      real_t tmp_33 = tmp_22 * tmp_24;
      real_t tmp_34 = tmp_24 * tmp_26;
      real_t tmp_35 = tmp_30 * ( tmp_19 * ( tmp_25 * tmp_32 + tmp_27 * tmp_31 - 1.0 / 3.0 ) +
                                 tmp_20 * ( tmp_31 * tmp_33 + tmp_32 * tmp_34 - 1.0 / 3.0 ) );
      real_t tmp_36 = 0.11846344252809471 * tmp_29;
      real_t tmp_37 = p_affine_6_1 + 0.23076534494715845 * tmp_5;
      real_t tmp_38 = tmp_4 * ( tmp_2 + tmp_37 );
      real_t tmp_39 = tmp_1 * tmp_38;
      real_t tmp_40 = tmp_38 * tmp_9;
      real_t tmp_41 = p_affine_6_0 + 0.23076534494715845 * tmp_11;
      real_t tmp_42 = tmp_4 * ( tmp_0 + tmp_41 );
      real_t tmp_43 = tmp_3 * tmp_42;
      real_t tmp_44 = tmp_15 * tmp_42;
      real_t tmp_45 = -tmp_39 - tmp_40 - tmp_43 - tmp_44 + 1;
      real_t tmp_46 = tmp_18 + tmp_37;
      real_t tmp_47 = tmp_21 + tmp_41;
      real_t tmp_48 = tmp_30 * ( tmp_19 * ( tmp_25 * tmp_47 + tmp_27 * tmp_46 - 1.0 / 3.0 ) +
                                 tmp_20 * ( tmp_33 * tmp_46 + tmp_34 * tmp_47 - 1.0 / 3.0 ) );
      real_t tmp_49 = 0.2393143352496831 * tmp_29;
      real_t tmp_50 = p_affine_6_1 + 0.5 * tmp_5;
      real_t tmp_51 = tmp_4 * ( tmp_2 + tmp_50 );
      real_t tmp_52 = tmp_1 * tmp_51;
      real_t tmp_53 = tmp_51 * tmp_9;
      real_t tmp_54 = p_affine_6_0 + 0.5 * tmp_11;
      real_t tmp_55 = tmp_4 * ( tmp_0 + tmp_54 );
      real_t tmp_56 = tmp_3 * tmp_55;
      real_t tmp_57 = tmp_15 * tmp_55;
      real_t tmp_58 = -tmp_52 - tmp_53 - tmp_56 - tmp_57 + 1;
      real_t tmp_59 = tmp_18 + tmp_50;
      real_t tmp_60 = tmp_21 + tmp_54;
      real_t tmp_61 = tmp_30 * ( tmp_19 * ( tmp_25 * tmp_60 + tmp_27 * tmp_59 - 1.0 / 3.0 ) +
                                 tmp_20 * ( tmp_33 * tmp_59 + tmp_34 * tmp_60 - 1.0 / 3.0 ) );
      real_t tmp_62 = 0.2844444444444445 * tmp_29;
      real_t tmp_63 = p_affine_6_1 + 0.7692346550528415 * tmp_5;
      real_t tmp_64 = tmp_4 * ( tmp_2 + tmp_63 );
      real_t tmp_65 = tmp_1 * tmp_64;
      real_t tmp_66 = tmp_64 * tmp_9;
      real_t tmp_67 = p_affine_6_0 + 0.7692346550528415 * tmp_11;
      real_t tmp_68 = tmp_4 * ( tmp_0 + tmp_67 );
      real_t tmp_69 = tmp_3 * tmp_68;
      real_t tmp_70 = tmp_15 * tmp_68;
      real_t tmp_71 = -tmp_65 - tmp_66 - tmp_69 - tmp_70 + 1;
      real_t tmp_72 = tmp_18 + tmp_63;
      real_t tmp_73 = tmp_21 + tmp_67;
      real_t tmp_74 = tmp_30 * ( tmp_19 * ( tmp_25 * tmp_73 + tmp_27 * tmp_72 - 1.0 / 3.0 ) +
                                 tmp_20 * ( tmp_33 * tmp_72 + tmp_34 * tmp_73 - 1.0 / 3.0 ) );
      real_t tmp_75 = 0.2393143352496831 * tmp_29;
      real_t tmp_76 = p_affine_6_1 + 0.95308992296933193 * tmp_5;
      real_t tmp_77 = tmp_4 * ( tmp_2 + tmp_76 );
      real_t tmp_78 = tmp_1 * tmp_77;
      real_t tmp_79 = tmp_77 * tmp_9;
      real_t tmp_80 = p_affine_6_0 + 0.95308992296933193 * tmp_11;
      real_t tmp_81 = tmp_4 * ( tmp_0 + tmp_80 );
      real_t tmp_82 = tmp_3 * tmp_81;
      real_t tmp_83 = tmp_15 * tmp_81;
      real_t tmp_84 = -tmp_78 - tmp_79 - tmp_82 - tmp_83 + 1;
      real_t tmp_85 = tmp_18 + tmp_76;
      real_t tmp_86 = tmp_21 + tmp_80;
      real_t tmp_87 = tmp_30 * ( tmp_19 * ( tmp_25 * tmp_86 + tmp_27 * tmp_85 - 1.0 / 3.0 ) +
                                 tmp_20 * ( tmp_33 * tmp_85 + tmp_34 * tmp_86 - 1.0 / 3.0 ) );
      real_t tmp_88 = 0.11846344252809471 * tmp_29;
      real_t tmp_89 = tmp_10 + tmp_14;
      real_t tmp_90 = tmp_40 + tmp_43;
      real_t tmp_91 = tmp_53 + tmp_56;
      real_t tmp_92 = tmp_66 + tmp_69;
      real_t tmp_93 = tmp_79 + tmp_82;
      real_t tmp_94 = tmp_16 + tmp_8;
      real_t tmp_95 = tmp_39 + tmp_44;
      real_t tmp_96 = tmp_52 + tmp_57;
      real_t tmp_97 = tmp_65 + tmp_70;
      real_t tmp_98 = tmp_78 + tmp_83;
      real_t a_0_0  = tmp_36 * ( -tmp_17 * tmp_28 - tmp_17 * tmp_35 ) + tmp_49 * ( -tmp_28 * tmp_45 - tmp_45 * tmp_48 ) +
                     tmp_62 * ( -tmp_28 * tmp_58 - tmp_58 * tmp_61 ) + tmp_75 * ( -tmp_28 * tmp_71 - tmp_71 * tmp_74 ) +
                     tmp_88 * ( -tmp_28 * tmp_84 - tmp_84 * tmp_87 );
      real_t a_1_0 = tmp_36 * ( -tmp_28 * tmp_89 - tmp_35 * tmp_89 ) + tmp_49 * ( -tmp_28 * tmp_90 - tmp_48 * tmp_90 ) +
                     tmp_62 * ( -tmp_28 * tmp_91 - tmp_61 * tmp_91 ) + tmp_75 * ( -tmp_28 * tmp_92 - tmp_74 * tmp_92 ) +
                     tmp_88 * ( -tmp_28 * tmp_93 - tmp_87 * tmp_93 );
      real_t a_2_0 = tmp_36 * ( -tmp_28 * tmp_94 - tmp_35 * tmp_94 ) + tmp_49 * ( -tmp_28 * tmp_95 - tmp_48 * tmp_95 ) +
                     tmp_62 * ( -tmp_28 * tmp_96 - tmp_61 * tmp_96 ) + tmp_75 * ( -tmp_28 * tmp_97 - tmp_74 * tmp_97 ) +
                     tmp_88 * ( -tmp_28 * tmp_98 - tmp_87 * tmp_98 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                                   const std::vector< Point3D >& coordsFacet,
                                                   const Point3D&                oppositeVertex,
                                                   const Point3D&                outwardNormal,
                                                   const DGBasisInfo&            trialBasis,
                                                   const DGBasisInfo&            testBasis,
                                                   int                           trialDegree,
                                                   int                           testDegree,
                                                   MatrixXr&                     elMat ) const override
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
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = p_affine_1_1 + tmp_2;
      real_t tmp_5  = 1.0 / ( tmp_1 * tmp_3 - tmp_4 * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_6  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_7  = p_affine_6_1 + tmp_2;
      real_t tmp_8  = tmp_5 * ( 0.046910077030668018 * tmp_6 + tmp_7 );
      real_t tmp_9  = tmp_1 * tmp_8;
      real_t tmp_10 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_11 = tmp_10 * tmp_8;
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_0;
      real_t tmp_14 = tmp_5 * ( 0.046910077030668018 * tmp_12 + tmp_13 );
      real_t tmp_15 = tmp_14 * tmp_3;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_14 * tmp_16;
      real_t tmp_18 = tmp_11 + tmp_15;
      real_t tmp_19 = tmp_17 + tmp_9;
      real_t tmp_20 = 0.11846344252809471 * tmp_3 * ( tmp_19 - 1.0 / 3.0 ) + 0.11846344252809471 * tmp_4 * ( tmp_18 - 1.0 / 3.0 );
      real_t tmp_21 = tmp_5 * ( 0.23076534494715845 * tmp_6 + tmp_7 );
      real_t tmp_22 = tmp_1 * tmp_21;
      real_t tmp_23 = tmp_10 * tmp_21;
      real_t tmp_24 = tmp_5 * ( 0.23076534494715845 * tmp_12 + tmp_13 );
      real_t tmp_25 = tmp_24 * tmp_3;
      real_t tmp_26 = tmp_16 * tmp_24;
      real_t tmp_27 = tmp_23 + tmp_25;
      real_t tmp_28 = tmp_22 + tmp_26;
      real_t tmp_29 = 0.2393143352496831 * tmp_3 * ( tmp_28 - 1.0 / 3.0 ) + 0.2393143352496831 * tmp_4 * ( tmp_27 - 1.0 / 3.0 );
      real_t tmp_30 = tmp_5 * ( 0.5 * tmp_6 + tmp_7 );
      real_t tmp_31 = tmp_1 * tmp_30;
      real_t tmp_32 = tmp_10 * tmp_30;
      real_t tmp_33 = tmp_5 * ( 0.5 * tmp_12 + tmp_13 );
      real_t tmp_34 = tmp_3 * tmp_33;
      real_t tmp_35 = tmp_16 * tmp_33;
      real_t tmp_36 = tmp_32 + tmp_34;
      real_t tmp_37 = tmp_31 + tmp_35;
      real_t tmp_38 = 0.2844444444444445 * tmp_3 * ( tmp_37 - 1.0 / 3.0 ) + 0.2844444444444445 * tmp_4 * ( tmp_36 - 1.0 / 3.0 );
      real_t tmp_39 = tmp_5 * ( 0.7692346550528415 * tmp_6 + tmp_7 );
      real_t tmp_40 = tmp_1 * tmp_39;
      real_t tmp_41 = tmp_10 * tmp_39;
      real_t tmp_42 = tmp_5 * ( 0.7692346550528415 * tmp_12 + tmp_13 );
      real_t tmp_43 = tmp_3 * tmp_42;
      real_t tmp_44 = tmp_16 * tmp_42;
      real_t tmp_45 = tmp_41 + tmp_43;
      real_t tmp_46 = tmp_40 + tmp_44;
      real_t tmp_47 = 0.2393143352496831 * tmp_3 * ( tmp_46 - 1.0 / 3.0 ) + 0.2393143352496831 * tmp_4 * ( tmp_45 - 1.0 / 3.0 );
      real_t tmp_48 = tmp_5 * ( 0.95308992296933193 * tmp_6 + tmp_7 );
      real_t tmp_49 = tmp_1 * tmp_48;
      real_t tmp_50 = tmp_10 * tmp_48;
      real_t tmp_51 = tmp_5 * ( 0.95308992296933193 * tmp_12 + tmp_13 );
      real_t tmp_52 = tmp_3 * tmp_51;
      real_t tmp_53 = tmp_16 * tmp_51;
      real_t tmp_54 = tmp_50 + tmp_52;
      real_t tmp_55 = tmp_49 + tmp_53;
      real_t tmp_56 = 0.11846344252809471 * tmp_3 * ( tmp_55 - 1.0 / 3.0 ) + 0.11846344252809471 * tmp_4 * ( tmp_54 - 1.0 / 3.0 );
      real_t a_0_0  = tmp_20 * ( -tmp_11 - tmp_15 - tmp_17 - tmp_9 + 1 ) + tmp_29 * ( -tmp_22 - tmp_23 - tmp_25 - tmp_26 + 1 ) +
                     tmp_38 * ( -tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1 ) + tmp_47 * ( -tmp_40 - tmp_41 - tmp_43 - tmp_44 + 1 ) +
                     tmp_56 * ( -tmp_49 - tmp_50 - tmp_52 - tmp_53 + 1 );
      real_t a_1_0  = tmp_18 * tmp_20 + tmp_27 * tmp_29 + tmp_36 * tmp_38 + tmp_45 * tmp_47 + tmp_54 * tmp_56;
      real_t a_2_0  = tmp_19 * tmp_20 + tmp_28 * tmp_29 + tmp_37 * tmp_38 + tmp_46 * tmp_47 + tmp_55 * tmp_56;
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
   }

   void integrateRHSDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                         const std::vector< Point3D >& coordsFacet,
                                         const Point3D&                oppositeVertex,
                                         const Point3D&                outwardNormal,
                                         const DGBasisInfo&            basis,
                                         int                           degree,
                                         MatrixXr&                     elMat ) const override
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
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                         const std::vector< Point3D >& coordsFacet,
                                         const Point3D&                oppositeVertex,
                                         const Point3D&                outwardNormal,
                                         const DGBasisInfo&            basis,
                                         int                           degree,
                                         MatrixXr&                     elMat ) const override
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
   void integrateVolume3D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_0;
      real_t tmp_6  = p_affine_1_1 + tmp_2;
      real_t tmp_7  = tmp_5 * tmp_6;
      real_t tmp_8  = -p_affine_0_2;
      real_t tmp_9  = p_affine_3_2 + tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_8;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11 * tmp_5;
      real_t tmp_13 = p_affine_2_2 + tmp_8;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = tmp_14 * tmp_6;
      real_t tmp_16 = tmp_1 * tmp_11;
      real_t tmp_17 = tmp_14 * tmp_3;
      real_t tmp_18 =
          1.0 / ( tmp_10 * tmp_12 - tmp_10 * tmp_17 + tmp_13 * tmp_15 - tmp_13 * tmp_16 + tmp_4 * tmp_9 - tmp_7 * tmp_9 );
      real_t tmp_19 = tmp_18 * ( tmp_4 - tmp_7 );
      real_t tmp_20 = tmp_18 * ( tmp_15 - tmp_16 );
      real_t tmp_21 = tmp_18 * ( tmp_12 - tmp_17 );
      real_t tmp_22 = tmp_11 * tmp_19 + tmp_20 * tmp_3 + tmp_21 * tmp_6;
      real_t tmp_23 = tmp_18 * ( -tmp_1 * tmp_13 + tmp_10 * tmp_5 );
      real_t tmp_24 = tmp_18 * ( tmp_1 * tmp_9 - tmp_10 * tmp_14 );
      real_t tmp_25 = tmp_18 * ( tmp_13 * tmp_14 - tmp_5 * tmp_9 );
      real_t tmp_26 = tmp_11 * tmp_23 + tmp_24 * tmp_3 + tmp_25 * tmp_6;
      real_t tmp_27 = tmp_18 * ( -tmp_10 * tmp_3 + tmp_13 * tmp_6 );
      real_t tmp_28 = tmp_18 * ( tmp_10 * tmp_11 - tmp_6 * tmp_9 );
      real_t tmp_29 = tmp_18 * ( -tmp_11 * tmp_13 + tmp_3 * tmp_9 );
      real_t tmp_30 = tmp_11 * tmp_27 + tmp_28 * tmp_3 + tmp_29 * tmp_6;
      real_t tmp_31 = p_affine_0_0 * p_affine_1_1;
      real_t tmp_32 = p_affine_0_0 * p_affine_1_2;
      real_t tmp_33 = p_affine_2_1 * p_affine_3_2;
      real_t tmp_34 = p_affine_0_1 * p_affine_1_0;
      real_t tmp_35 = p_affine_0_1 * p_affine_1_2;
      real_t tmp_36 = p_affine_2_2 * p_affine_3_0;
      real_t tmp_37 = p_affine_0_2 * p_affine_1_0;
      real_t tmp_38 = p_affine_0_2 * p_affine_1_1;
      real_t tmp_39 = p_affine_2_0 * p_affine_3_1;
      real_t tmp_40 = p_affine_2_2 * p_affine_3_1;
      real_t tmp_41 = p_affine_2_0 * p_affine_3_2;
      real_t tmp_42 = p_affine_2_1 * p_affine_3_0;
      real_t tmp_43 = std::abs( p_affine_0_0 * tmp_33 - p_affine_0_0 * tmp_40 + p_affine_0_1 * tmp_36 - p_affine_0_1 * tmp_41 +
                                p_affine_0_2 * tmp_39 - p_affine_0_2 * tmp_42 - p_affine_1_0 * tmp_33 + p_affine_1_0 * tmp_40 -
                                p_affine_1_1 * tmp_36 + p_affine_1_1 * tmp_41 - p_affine_1_2 * tmp_39 + p_affine_1_2 * tmp_42 +
                                p_affine_2_0 * tmp_35 - p_affine_2_0 * tmp_38 - p_affine_2_1 * tmp_32 + p_affine_2_1 * tmp_37 +
                                p_affine_2_2 * tmp_31 - p_affine_2_2 * tmp_34 - p_affine_3_0 * tmp_35 + p_affine_3_0 * tmp_38 +
                                p_affine_3_1 * tmp_32 - p_affine_3_1 * tmp_37 - p_affine_3_2 * tmp_31 + p_affine_3_2 * tmp_34 );
      real_t tmp_44 = tmp_43 * ( tmp_22 * ( -tmp_19 - tmp_20 - tmp_21 ) + tmp_26 * ( -tmp_23 - tmp_24 - tmp_25 ) +
                                 tmp_30 * ( -tmp_27 - tmp_28 - tmp_29 ) );
      real_t tmp_45 = tmp_43 * ( tmp_21 * tmp_22 + tmp_25 * tmp_26 + tmp_29 * tmp_30 );
      real_t tmp_46 = tmp_43 * ( tmp_20 * tmp_22 + tmp_24 * tmp_26 + tmp_28 * tmp_30 );
      real_t tmp_47 = tmp_43 * ( tmp_19 * tmp_22 + tmp_23 * tmp_26 + tmp_27 * tmp_30 );
      real_t a_0_0  = 0.1666666666666668 * tmp_44;
      real_t a_1_0  = 0.1666666666666668 * tmp_45;
      real_t a_2_0  = 0.1666666666666668 * tmp_46;
      real_t a_3_0  = 0.1666666666666668 * tmp_47;
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
      elMat( 3, 0 ) = a_3_0;
   }

   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                               const std::vector< Point3D >& coordsFacet,
                               const Point3D&,
                               const Point3D&     outwardNormal,
                               const DGBasisInfo& trialBasis,
                               const DGBasisInfo& testBasis,
                               int                trialDegree,
                               int                testDegree,
                               MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = p_affine_2_0 + tmp_0;
      real_t tmp_5  = p_affine_1_1 + tmp_2;
      real_t tmp_6  = tmp_1 * tmp_3 - tmp_4 * tmp_5;
      real_t tmp_7  = -p_affine_0_2;
      real_t tmp_8  = p_affine_3_2 + tmp_7;
      real_t tmp_9  = tmp_3 * tmp_8;
      real_t tmp_10 = p_affine_3_1 + tmp_2;
      real_t tmp_11 = p_affine_1_2 + tmp_7;
      real_t tmp_12 = tmp_10 * tmp_11;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_14 * tmp_5;
      real_t tmp_16 = tmp_10 * tmp_14;
      real_t tmp_17 = tmp_5 * tmp_8;
      real_t tmp_18 = tmp_11 * tmp_3;
      real_t tmp_19 =
          1.0 / ( -tmp_1 * tmp_16 + tmp_1 * tmp_9 + tmp_12 * tmp_4 + tmp_13 * tmp_15 - tmp_13 * tmp_18 - tmp_17 * tmp_4 );
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19 * ( 0.031405749086161582 * tmp_21 + 0.93718850182767688 * tmp_22 + tmp_23 );
      real_t tmp_25 = tmp_24 * tmp_6;
      real_t tmp_26 = -tmp_1 * tmp_10 + tmp_13 * tmp_5;
      real_t tmp_27 = tmp_24 * tmp_26;
      real_t tmp_28 = -tmp_1 * tmp_14 + tmp_11 * tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19 * ( 0.031405749086161582 * tmp_30 + 0.93718850182767688 * tmp_31 + tmp_32 );
      real_t tmp_34 = tmp_28 * tmp_33;
      real_t tmp_35 = tmp_1 * tmp_8 - tmp_11 * tmp_13;
      real_t tmp_36 = tmp_33 * tmp_35;
      real_t tmp_37 = tmp_10 * tmp_4 - tmp_13 * tmp_3;
      real_t tmp_38 = tmp_24 * tmp_37;
      real_t tmp_39 = tmp_13 * tmp_14 - tmp_4 * tmp_8;
      real_t tmp_40 = tmp_33 * tmp_39;
      real_t tmp_41 = tmp_15 - tmp_18;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19 * ( 0.031405749086161582 * tmp_43 + 0.93718850182767688 * tmp_44 + tmp_45 );
      real_t tmp_47 = tmp_41 * tmp_46;
      real_t tmp_48 = tmp_12 - tmp_17;
      real_t tmp_49 = tmp_46 * tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46 * tmp_50;
      real_t tmp_52 = -tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1;
      real_t tmp_53 = tmp_19 * tmp_5;
      real_t tmp_54 = tmp_19 * tmp_3;
      real_t tmp_55 = tmp_10 * tmp_19;
      real_t tmp_56 = 0.5 * p_affine_13_0 * ( tmp_41 * tmp_55 + tmp_48 * tmp_54 + tmp_50 * tmp_53 ) +
                      0.5 * p_affine_13_1 * ( tmp_28 * tmp_55 + tmp_35 * tmp_54 + tmp_39 * tmp_53 ) +
                      0.5 * p_affine_13_2 * ( tmp_26 * tmp_54 + tmp_37 * tmp_53 + tmp_55 * tmp_6 );
      real_t tmp_57 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_58 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_59 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_60 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_61 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_62 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_63 = ( std::abs( tmp_22 * tmp_60 - tmp_31 * tmp_62 ) * std::abs( tmp_22 * tmp_60 - tmp_31 * tmp_62 ) ) +
                      ( std::abs( tmp_22 * tmp_61 - tmp_44 * tmp_62 ) * std::abs( tmp_22 * tmp_61 - tmp_44 * tmp_62 ) ) +
                      ( std::abs( tmp_31 * tmp_61 - tmp_44 * tmp_60 ) * std::abs( tmp_31 * tmp_61 - tmp_44 * tmp_60 ) );
      real_t tmp_64 = 1.0 * std::pow( tmp_63, -0.25 );
      real_t tmp_65 =
          tmp_64 * ( tmp_10 * ( tmp_59 - 1.0 / 4.0 ) + tmp_3 * ( tmp_58 - 1.0 / 4.0 ) + tmp_5 * ( tmp_57 - 1.0 / 4.0 ) );
      real_t tmp_66 = 1.0 * std::pow( tmp_63, 1.0 / 2.0 );
      real_t tmp_67 = 0.0068572537431980923 * tmp_66;
      real_t tmp_68 = tmp_19 * ( 0.19601935860219369 * tmp_21 + 0.60796128279561268 * tmp_22 + tmp_23 );
      real_t tmp_69 = tmp_6 * tmp_68;
      real_t tmp_70 = tmp_26 * tmp_68;
      real_t tmp_71 = tmp_19 * ( 0.19601935860219369 * tmp_30 + 0.60796128279561268 * tmp_31 + tmp_32 );
      real_t tmp_72 = tmp_28 * tmp_71;
      real_t tmp_73 = tmp_35 * tmp_71;
      real_t tmp_74 = tmp_37 * tmp_68;
      real_t tmp_75 = tmp_39 * tmp_71;
      real_t tmp_76 = tmp_19 * ( 0.19601935860219369 * tmp_43 + 0.60796128279561268 * tmp_44 + tmp_45 );
      real_t tmp_77 = tmp_41 * tmp_76;
      real_t tmp_78 = tmp_48 * tmp_76;
      real_t tmp_79 = tmp_50 * tmp_76;
      real_t tmp_80 = -tmp_69 - tmp_70 - tmp_72 - tmp_73 - tmp_74 - tmp_75 - tmp_77 - tmp_78 - tmp_79 + 1;
      real_t tmp_81 = tmp_74 + tmp_75 + tmp_79;
      real_t tmp_82 = tmp_70 + tmp_73 + tmp_78;
      real_t tmp_83 = tmp_69 + tmp_72 + tmp_77;
      real_t tmp_84 =
          tmp_64 * ( tmp_10 * ( tmp_83 - 1.0 / 4.0 ) + tmp_3 * ( tmp_82 - 1.0 / 4.0 ) + tmp_5 * ( tmp_81 - 1.0 / 4.0 ) );
      real_t tmp_85  = 0.037198804536718075 * tmp_66;
      real_t tmp_86  = tmp_19 * ( 0.37605877282253791 * tmp_21 + 0.039308471900058539 * tmp_22 + tmp_23 );
      real_t tmp_87  = tmp_6 * tmp_86;
      real_t tmp_88  = tmp_26 * tmp_86;
      real_t tmp_89  = tmp_19 * ( 0.37605877282253791 * tmp_30 + 0.039308471900058539 * tmp_31 + tmp_32 );
      real_t tmp_90  = tmp_28 * tmp_89;
      real_t tmp_91  = tmp_35 * tmp_89;
      real_t tmp_92  = tmp_37 * tmp_86;
      real_t tmp_93  = tmp_39 * tmp_89;
      real_t tmp_94  = tmp_19 * ( 0.37605877282253791 * tmp_43 + 0.039308471900058539 * tmp_44 + tmp_45 );
      real_t tmp_95  = tmp_41 * tmp_94;
      real_t tmp_96  = tmp_48 * tmp_94;
      real_t tmp_97  = tmp_50 * tmp_94;
      real_t tmp_98  = -tmp_87 - tmp_88 - tmp_90 - tmp_91 - tmp_92 - tmp_93 - tmp_95 - tmp_96 - tmp_97 + 1;
      real_t tmp_99  = tmp_92 + tmp_93 + tmp_97;
      real_t tmp_100 = tmp_88 + tmp_91 + tmp_96;
      real_t tmp_101 = tmp_87 + tmp_90 + tmp_95;
      real_t tmp_102 =
          tmp_64 * ( tmp_10 * ( tmp_101 - 1.0 / 4.0 ) + tmp_3 * ( tmp_100 - 1.0 / 4.0 ) + tmp_5 * ( tmp_99 - 1.0 / 4.0 ) );
      real_t tmp_103 = 0.020848748529055869 * tmp_66;
      real_t tmp_104 = tmp_19 * ( 0.78764240869137092 * tmp_21 + 0.1711304259088916 * tmp_22 + tmp_23 );
      real_t tmp_105 = tmp_104 * tmp_6;
      real_t tmp_106 = tmp_104 * tmp_26;
      real_t tmp_107 = tmp_19 * ( 0.78764240869137092 * tmp_30 + 0.1711304259088916 * tmp_31 + tmp_32 );
      real_t tmp_108 = tmp_107 * tmp_28;
      real_t tmp_109 = tmp_107 * tmp_35;
      real_t tmp_110 = tmp_104 * tmp_37;
      real_t tmp_111 = tmp_107 * tmp_39;
      real_t tmp_112 = tmp_19 * ( 0.78764240869137092 * tmp_43 + 0.1711304259088916 * tmp_44 + tmp_45 );
      real_t tmp_113 = tmp_112 * tmp_41;
      real_t tmp_114 = tmp_112 * tmp_48;
      real_t tmp_115 = tmp_112 * tmp_50;
      real_t tmp_116 = -tmp_105 - tmp_106 - tmp_108 - tmp_109 - tmp_110 - tmp_111 - tmp_113 - tmp_114 - tmp_115 + 1;
      real_t tmp_117 = tmp_110 + tmp_111 + tmp_115;
      real_t tmp_118 = tmp_106 + tmp_109 + tmp_114;
      real_t tmp_119 = tmp_105 + tmp_108 + tmp_113;
      real_t tmp_120 =
          tmp_64 * ( tmp_10 * ( tmp_119 - 1.0 / 4.0 ) + tmp_3 * ( tmp_118 - 1.0 / 4.0 ) + tmp_5 * ( tmp_117 - 1.0 / 4.0 ) );
      real_t tmp_121 = 0.019202922745021479 * tmp_66;
      real_t tmp_122 = tmp_19 * ( 0.58463275527740355 * tmp_21 + 0.37605877282253791 * tmp_22 + tmp_23 );
      real_t tmp_123 = tmp_122 * tmp_6;
      real_t tmp_124 = tmp_122 * tmp_26;
      real_t tmp_125 = tmp_19 * ( 0.58463275527740355 * tmp_30 + 0.37605877282253791 * tmp_31 + tmp_32 );
      real_t tmp_126 = tmp_125 * tmp_28;
      real_t tmp_127 = tmp_125 * tmp_35;
      real_t tmp_128 = tmp_122 * tmp_37;
      real_t tmp_129 = tmp_125 * tmp_39;
      real_t tmp_130 = tmp_19 * ( 0.58463275527740355 * tmp_43 + 0.37605877282253791 * tmp_44 + tmp_45 );
      real_t tmp_131 = tmp_130 * tmp_41;
      real_t tmp_132 = tmp_130 * tmp_48;
      real_t tmp_133 = tmp_130 * tmp_50;
      real_t tmp_134 = -tmp_123 - tmp_124 - tmp_126 - tmp_127 - tmp_128 - tmp_129 - tmp_131 - tmp_132 - tmp_133 + 1;
      real_t tmp_135 = tmp_128 + tmp_129 + tmp_133;
      real_t tmp_136 = tmp_124 + tmp_127 + tmp_132;
      real_t tmp_137 = tmp_123 + tmp_126 + tmp_131;
      real_t tmp_138 =
          tmp_64 * ( tmp_10 * ( tmp_137 - 1.0 / 4.0 ) + tmp_3 * ( tmp_136 - 1.0 / 4.0 ) + tmp_5 * ( tmp_135 - 1.0 / 4.0 ) );
      real_t tmp_139 = 0.020848748529055869 * tmp_66;
      real_t tmp_140 = tmp_19 * ( 0.041227165399737475 * tmp_21 + 0.78764240869137092 * tmp_22 + tmp_23 );
      real_t tmp_141 = tmp_140 * tmp_6;
      real_t tmp_142 = tmp_140 * tmp_26;
      real_t tmp_143 = tmp_19 * ( 0.041227165399737475 * tmp_30 + 0.78764240869137092 * tmp_31 + tmp_32 );
      real_t tmp_144 = tmp_143 * tmp_28;
      real_t tmp_145 = tmp_143 * tmp_35;
      real_t tmp_146 = tmp_140 * tmp_37;
      real_t tmp_147 = tmp_143 * tmp_39;
      real_t tmp_148 = tmp_19 * ( 0.041227165399737475 * tmp_43 + 0.78764240869137092 * tmp_44 + tmp_45 );
      real_t tmp_149 = tmp_148 * tmp_41;
      real_t tmp_150 = tmp_148 * tmp_48;
      real_t tmp_151 = tmp_148 * tmp_50;
      real_t tmp_152 = -tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1;
      real_t tmp_153 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_154 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_155 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_156 =
          tmp_64 * ( tmp_10 * ( tmp_155 - 1.0 / 4.0 ) + tmp_3 * ( tmp_154 - 1.0 / 4.0 ) + tmp_5 * ( tmp_153 - 1.0 / 4.0 ) );
      real_t tmp_157 = 0.019202922745021479 * tmp_66;
      real_t tmp_158 = tmp_19 * ( 0.039308471900058539 * tmp_21 + 0.58463275527740355 * tmp_22 + tmp_23 );
      real_t tmp_159 = tmp_158 * tmp_6;
      real_t tmp_160 = tmp_158 * tmp_26;
      real_t tmp_161 = tmp_19 * ( 0.039308471900058539 * tmp_30 + 0.58463275527740355 * tmp_31 + tmp_32 );
      real_t tmp_162 = tmp_161 * tmp_28;
      real_t tmp_163 = tmp_161 * tmp_35;
      real_t tmp_164 = tmp_158 * tmp_37;
      real_t tmp_165 = tmp_161 * tmp_39;
      real_t tmp_166 = tmp_19 * ( 0.039308471900058539 * tmp_43 + 0.58463275527740355 * tmp_44 + tmp_45 );
      real_t tmp_167 = tmp_166 * tmp_41;
      real_t tmp_168 = tmp_166 * tmp_48;
      real_t tmp_169 = tmp_166 * tmp_50;
      real_t tmp_170 = -tmp_159 - tmp_160 - tmp_162 - tmp_163 - tmp_164 - tmp_165 - tmp_167 - tmp_168 - tmp_169 + 1;
      real_t tmp_171 = tmp_164 + tmp_165 + tmp_169;
      real_t tmp_172 = tmp_160 + tmp_163 + tmp_168;
      real_t tmp_173 = tmp_159 + tmp_162 + tmp_167;
      real_t tmp_174 =
          tmp_64 * ( tmp_10 * ( tmp_173 - 1.0 / 4.0 ) + tmp_3 * ( tmp_172 - 1.0 / 4.0 ) + tmp_5 * ( tmp_171 - 1.0 / 4.0 ) );
      real_t tmp_175 = 0.020848748529055869 * tmp_66;
      real_t tmp_176 = tmp_19 * ( 0.78764240869137092 * tmp_21 + 0.041227165399737475 * tmp_22 + tmp_23 );
      real_t tmp_177 = tmp_176 * tmp_6;
      real_t tmp_178 = tmp_176 * tmp_26;
      real_t tmp_179 = tmp_19 * ( 0.78764240869137092 * tmp_30 + 0.041227165399737475 * tmp_31 + tmp_32 );
      real_t tmp_180 = tmp_179 * tmp_28;
      real_t tmp_181 = tmp_179 * tmp_35;
      real_t tmp_182 = tmp_176 * tmp_37;
      real_t tmp_183 = tmp_179 * tmp_39;
      real_t tmp_184 = tmp_19 * ( 0.78764240869137092 * tmp_43 + 0.041227165399737475 * tmp_44 + tmp_45 );
      real_t tmp_185 = tmp_184 * tmp_41;
      real_t tmp_186 = tmp_184 * tmp_48;
      real_t tmp_187 = tmp_184 * tmp_50;
      real_t tmp_188 = -tmp_177 - tmp_178 - tmp_180 - tmp_181 - tmp_182 - tmp_183 - tmp_185 - tmp_186 - tmp_187 + 1;
      real_t tmp_189 = tmp_182 + tmp_183 + tmp_187;
      real_t tmp_190 = tmp_178 + tmp_181 + tmp_186;
      real_t tmp_191 = tmp_177 + tmp_180 + tmp_185;
      real_t tmp_192 =
          tmp_64 * ( tmp_10 * ( tmp_191 - 1.0 / 4.0 ) + tmp_3 * ( tmp_190 - 1.0 / 4.0 ) + tmp_5 * ( tmp_189 - 1.0 / 4.0 ) );
      real_t tmp_193 = 0.019202922745021479 * tmp_66;
      real_t tmp_194 = tmp_19 * ( 0.58463275527740355 * tmp_21 + 0.039308471900058539 * tmp_22 + tmp_23 );
      real_t tmp_195 = tmp_194 * tmp_6;
      real_t tmp_196 = tmp_194 * tmp_26;
      real_t tmp_197 = tmp_19 * ( 0.58463275527740355 * tmp_30 + 0.039308471900058539 * tmp_31 + tmp_32 );
      real_t tmp_198 = tmp_197 * tmp_28;
      real_t tmp_199 = tmp_197 * tmp_35;
      real_t tmp_200 = tmp_194 * tmp_37;
      real_t tmp_201 = tmp_197 * tmp_39;
      real_t tmp_202 = tmp_19 * ( 0.58463275527740355 * tmp_43 + 0.039308471900058539 * tmp_44 + tmp_45 );
      real_t tmp_203 = tmp_202 * tmp_41;
      real_t tmp_204 = tmp_202 * tmp_48;
      real_t tmp_205 = tmp_202 * tmp_50;
      real_t tmp_206 = -tmp_195 - tmp_196 - tmp_198 - tmp_199 - tmp_200 - tmp_201 - tmp_203 - tmp_204 - tmp_205 + 1;
      real_t tmp_207 = tmp_200 + tmp_201 + tmp_205;
      real_t tmp_208 = tmp_196 + tmp_199 + tmp_204;
      real_t tmp_209 = tmp_195 + tmp_198 + tmp_203;
      real_t tmp_210 =
          tmp_64 * ( tmp_10 * ( tmp_209 - 1.0 / 4.0 ) + tmp_3 * ( tmp_208 - 1.0 / 4.0 ) + tmp_5 * ( tmp_207 - 1.0 / 4.0 ) );
      real_t tmp_211 = 0.020848748529055869 * tmp_66;
      real_t tmp_212 = tmp_19 * ( 0.1711304259088916 * tmp_21 + 0.78764240869137092 * tmp_22 + tmp_23 );
      real_t tmp_213 = tmp_212 * tmp_6;
      real_t tmp_214 = tmp_212 * tmp_26;
      real_t tmp_215 = tmp_19 * ( 0.1711304259088916 * tmp_30 + 0.78764240869137092 * tmp_31 + tmp_32 );
      real_t tmp_216 = tmp_215 * tmp_28;
      real_t tmp_217 = tmp_215 * tmp_35;
      real_t tmp_218 = tmp_212 * tmp_37;
      real_t tmp_219 = tmp_215 * tmp_39;
      real_t tmp_220 = tmp_19 * ( 0.1711304259088916 * tmp_43 + 0.78764240869137092 * tmp_44 + tmp_45 );
      real_t tmp_221 = tmp_220 * tmp_41;
      real_t tmp_222 = tmp_220 * tmp_48;
      real_t tmp_223 = tmp_220 * tmp_50;
      real_t tmp_224 = -tmp_213 - tmp_214 - tmp_216 - tmp_217 - tmp_218 - tmp_219 - tmp_221 - tmp_222 - tmp_223 + 1;
      real_t tmp_225 = tmp_218 + tmp_219 + tmp_223;
      real_t tmp_226 = tmp_214 + tmp_217 + tmp_222;
      real_t tmp_227 = tmp_213 + tmp_216 + tmp_221;
      real_t tmp_228 =
          tmp_64 * ( tmp_10 * ( tmp_227 - 1.0 / 4.0 ) + tmp_3 * ( tmp_226 - 1.0 / 4.0 ) + tmp_5 * ( tmp_225 - 1.0 / 4.0 ) );
      real_t tmp_229 = 0.019202922745021479 * tmp_66;
      real_t tmp_230 = tmp_19 * ( 0.37605877282253791 * tmp_21 + 0.58463275527740355 * tmp_22 + tmp_23 );
      real_t tmp_231 = tmp_230 * tmp_6;
      real_t tmp_232 = tmp_230 * tmp_26;
      real_t tmp_233 = tmp_19 * ( 0.37605877282253791 * tmp_30 + 0.58463275527740355 * tmp_31 + tmp_32 );
      real_t tmp_234 = tmp_233 * tmp_28;
      real_t tmp_235 = tmp_233 * tmp_35;
      real_t tmp_236 = tmp_230 * tmp_37;
      real_t tmp_237 = tmp_233 * tmp_39;
      real_t tmp_238 = tmp_19 * ( 0.37605877282253791 * tmp_43 + 0.58463275527740355 * tmp_44 + tmp_45 );
      real_t tmp_239 = tmp_238 * tmp_41;
      real_t tmp_240 = tmp_238 * tmp_48;
      real_t tmp_241 = tmp_238 * tmp_50;
      real_t tmp_242 = -tmp_231 - tmp_232 - tmp_234 - tmp_235 - tmp_236 - tmp_237 - tmp_239 - tmp_240 - tmp_241 + 1;
      real_t tmp_243 = tmp_236 + tmp_237 + tmp_241;
      real_t tmp_244 = tmp_232 + tmp_235 + tmp_240;
      real_t tmp_245 = tmp_231 + tmp_234 + tmp_239;
      real_t tmp_246 =
          tmp_64 * ( tmp_10 * ( tmp_245 - 1.0 / 4.0 ) + tmp_3 * ( tmp_244 - 1.0 / 4.0 ) + tmp_5 * ( tmp_243 - 1.0 / 4.0 ) );
      real_t tmp_247 = 0.020848748529055869 * tmp_66;
      real_t tmp_248 = tmp_19 * ( 0.041227165399737475 * tmp_21 + 0.1711304259088916 * tmp_22 + tmp_23 );
      real_t tmp_249 = tmp_248 * tmp_6;
      real_t tmp_250 = tmp_248 * tmp_26;
      real_t tmp_251 = tmp_19 * ( 0.041227165399737475 * tmp_30 + 0.1711304259088916 * tmp_31 + tmp_32 );
      real_t tmp_252 = tmp_251 * tmp_28;
      real_t tmp_253 = tmp_251 * tmp_35;
      real_t tmp_254 = tmp_248 * tmp_37;
      real_t tmp_255 = tmp_251 * tmp_39;
      real_t tmp_256 = tmp_19 * ( 0.041227165399737475 * tmp_43 + 0.1711304259088916 * tmp_44 + tmp_45 );
      real_t tmp_257 = tmp_256 * tmp_41;
      real_t tmp_258 = tmp_256 * tmp_48;
      real_t tmp_259 = tmp_256 * tmp_50;
      real_t tmp_260 = -tmp_249 - tmp_250 - tmp_252 - tmp_253 - tmp_254 - tmp_255 - tmp_257 - tmp_258 - tmp_259 + 1;
      real_t tmp_261 = tmp_254 + tmp_255 + tmp_259;
      real_t tmp_262 = tmp_250 + tmp_253 + tmp_258;
      real_t tmp_263 = tmp_249 + tmp_252 + tmp_257;
      real_t tmp_264 =
          tmp_64 * ( tmp_10 * ( tmp_263 - 1.0 / 4.0 ) + tmp_3 * ( tmp_262 - 1.0 / 4.0 ) + tmp_5 * ( tmp_261 - 1.0 / 4.0 ) );
      real_t tmp_265 = 0.019202922745021479 * tmp_66;
      real_t tmp_266 = tmp_19 * ( 0.40446199974765351 * tmp_21 + 0.19107600050469298 * tmp_22 + tmp_23 );
      real_t tmp_267 = tmp_266 * tmp_6;
      real_t tmp_268 = tmp_26 * tmp_266;
      real_t tmp_269 = tmp_19 * ( 0.40446199974765351 * tmp_30 + 0.19107600050469298 * tmp_31 + tmp_32 );
      real_t tmp_270 = tmp_269 * tmp_28;
      real_t tmp_271 = tmp_269 * tmp_35;
      real_t tmp_272 = tmp_266 * tmp_37;
      real_t tmp_273 = tmp_269 * tmp_39;
      real_t tmp_274 = tmp_19 * ( 0.40446199974765351 * tmp_43 + 0.19107600050469298 * tmp_44 + tmp_45 );
      real_t tmp_275 = tmp_274 * tmp_41;
      real_t tmp_276 = tmp_274 * tmp_48;
      real_t tmp_277 = tmp_274 * tmp_50;
      real_t tmp_278 = -tmp_267 - tmp_268 - tmp_270 - tmp_271 - tmp_272 - tmp_273 - tmp_275 - tmp_276 - tmp_277 + 1;
      real_t tmp_279 = tmp_272 + tmp_273 + tmp_277;
      real_t tmp_280 = tmp_268 + tmp_271 + tmp_276;
      real_t tmp_281 = tmp_267 + tmp_270 + tmp_275;
      real_t tmp_282 =
          tmp_64 * ( tmp_10 * ( tmp_281 - 1.0 / 4.0 ) + tmp_3 * ( tmp_280 - 1.0 / 4.0 ) + tmp_5 * ( tmp_279 - 1.0 / 4.0 ) );
      real_t tmp_283 = 0.042507265838595799 * tmp_66;
      real_t tmp_284 = tmp_19 * ( 0.039308471900058539 * tmp_21 + 0.37605877282253791 * tmp_22 + tmp_23 );
      real_t tmp_285 = tmp_284 * tmp_6;
      real_t tmp_286 = tmp_26 * tmp_284;
      real_t tmp_287 = tmp_19 * ( 0.039308471900058539 * tmp_30 + 0.37605877282253791 * tmp_31 + tmp_32 );
      real_t tmp_288 = tmp_28 * tmp_287;
      real_t tmp_289 = tmp_287 * tmp_35;
      real_t tmp_290 = tmp_284 * tmp_37;
      real_t tmp_291 = tmp_287 * tmp_39;
      real_t tmp_292 = tmp_19 * ( 0.039308471900058539 * tmp_43 + 0.37605877282253791 * tmp_44 + tmp_45 );
      real_t tmp_293 = tmp_292 * tmp_41;
      real_t tmp_294 = tmp_292 * tmp_48;
      real_t tmp_295 = tmp_292 * tmp_50;
      real_t tmp_296 = -tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1;
      real_t tmp_297 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_298 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_299 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_300 =
          tmp_64 * ( tmp_10 * ( tmp_299 - 1.0 / 4.0 ) + tmp_3 * ( tmp_298 - 1.0 / 4.0 ) + tmp_5 * ( tmp_297 - 1.0 / 4.0 ) );
      real_t tmp_301 = 0.020848748529055869 * tmp_66;
      real_t tmp_302 = tmp_19 * ( 0.93718850182767688 * tmp_21 + 0.031405749086161582 * tmp_22 + tmp_23 );
      real_t tmp_303 = tmp_302 * tmp_6;
      real_t tmp_304 = tmp_26 * tmp_302;
      real_t tmp_305 = tmp_19 * ( 0.93718850182767688 * tmp_30 + 0.031405749086161582 * tmp_31 + tmp_32 );
      real_t tmp_306 = tmp_28 * tmp_305;
      real_t tmp_307 = tmp_305 * tmp_35;
      real_t tmp_308 = tmp_302 * tmp_37;
      real_t tmp_309 = tmp_305 * tmp_39;
      real_t tmp_310 = tmp_19 * ( 0.93718850182767688 * tmp_43 + 0.031405749086161582 * tmp_44 + tmp_45 );
      real_t tmp_311 = tmp_310 * tmp_41;
      real_t tmp_312 = tmp_310 * tmp_48;
      real_t tmp_313 = tmp_310 * tmp_50;
      real_t tmp_314 = -tmp_303 - tmp_304 - tmp_306 - tmp_307 - tmp_308 - tmp_309 - tmp_311 - tmp_312 - tmp_313 + 1;
      real_t tmp_315 = tmp_308 + tmp_309 + tmp_313;
      real_t tmp_316 = tmp_304 + tmp_307 + tmp_312;
      real_t tmp_317 = tmp_303 + tmp_306 + tmp_311;
      real_t tmp_318 =
          tmp_64 * ( tmp_10 * ( tmp_317 - 1.0 / 4.0 ) + tmp_3 * ( tmp_316 - 1.0 / 4.0 ) + tmp_5 * ( tmp_315 - 1.0 / 4.0 ) );
      real_t tmp_319 = 0.0068572537431980923 * tmp_66;
      real_t tmp_320 = tmp_19 * ( 0.60796128279561268 * tmp_21 + 0.19601935860219369 * tmp_22 + tmp_23 );
      real_t tmp_321 = tmp_320 * tmp_6;
      real_t tmp_322 = tmp_26 * tmp_320;
      real_t tmp_323 = tmp_19 * ( 0.60796128279561268 * tmp_30 + 0.19601935860219369 * tmp_31 + tmp_32 );
      real_t tmp_324 = tmp_28 * tmp_323;
      real_t tmp_325 = tmp_323 * tmp_35;
      real_t tmp_326 = tmp_320 * tmp_37;
      real_t tmp_327 = tmp_323 * tmp_39;
      real_t tmp_328 = tmp_19 * ( 0.60796128279561268 * tmp_43 + 0.19601935860219369 * tmp_44 + tmp_45 );
      real_t tmp_329 = tmp_328 * tmp_41;
      real_t tmp_330 = tmp_328 * tmp_48;
      real_t tmp_331 = tmp_328 * tmp_50;
      real_t tmp_332 = -tmp_321 - tmp_322 - tmp_324 - tmp_325 - tmp_326 - tmp_327 - tmp_329 - tmp_330 - tmp_331 + 1;
      real_t tmp_333 = tmp_326 + tmp_327 + tmp_331;
      real_t tmp_334 = tmp_322 + tmp_325 + tmp_330;
      real_t tmp_335 = tmp_321 + tmp_324 + tmp_329;
      real_t tmp_336 =
          tmp_64 * ( tmp_10 * ( tmp_335 - 1.0 / 4.0 ) + tmp_3 * ( tmp_334 - 1.0 / 4.0 ) + tmp_5 * ( tmp_333 - 1.0 / 4.0 ) );
      real_t tmp_337 = 0.037198804536718075 * tmp_66;
      real_t tmp_338 = tmp_19 * ( 0.19107600050469298 * tmp_21 + 0.40446199974765351 * tmp_22 + tmp_23 );
      real_t tmp_339 = tmp_338 * tmp_6;
      real_t tmp_340 = tmp_26 * tmp_338;
      real_t tmp_341 = tmp_19 * ( 0.19107600050469298 * tmp_30 + 0.40446199974765351 * tmp_31 + tmp_32 );
      real_t tmp_342 = tmp_28 * tmp_341;
      real_t tmp_343 = tmp_341 * tmp_35;
      real_t tmp_344 = tmp_338 * tmp_37;
      real_t tmp_345 = tmp_341 * tmp_39;
      real_t tmp_346 = tmp_19 * ( 0.19107600050469298 * tmp_43 + 0.40446199974765351 * tmp_44 + tmp_45 );
      real_t tmp_347 = tmp_346 * tmp_41;
      real_t tmp_348 = tmp_346 * tmp_48;
      real_t tmp_349 = tmp_346 * tmp_50;
      real_t tmp_350 = -tmp_339 - tmp_340 - tmp_342 - tmp_343 - tmp_344 - tmp_345 - tmp_347 - tmp_348 - tmp_349 + 1;
      real_t tmp_351 = tmp_344 + tmp_345 + tmp_349;
      real_t tmp_352 = tmp_340 + tmp_343 + tmp_348;
      real_t tmp_353 = tmp_339 + tmp_342 + tmp_347;
      real_t tmp_354 =
          tmp_64 * ( tmp_10 * ( tmp_353 - 1.0 / 4.0 ) + tmp_3 * ( tmp_352 - 1.0 / 4.0 ) + tmp_5 * ( tmp_351 - 1.0 / 4.0 ) );
      real_t tmp_355 = 0.042507265838595799 * tmp_66;
      real_t tmp_356 = tmp_19 * ( 0.031405749086161582 * tmp_21 + 0.031405749086161582 * tmp_22 + tmp_23 );
      real_t tmp_357 = tmp_356 * tmp_6;
      real_t tmp_358 = tmp_26 * tmp_356;
      real_t tmp_359 = tmp_19 * ( 0.031405749086161582 * tmp_30 + 0.031405749086161582 * tmp_31 + tmp_32 );
      real_t tmp_360 = tmp_28 * tmp_359;
      real_t tmp_361 = tmp_35 * tmp_359;
      real_t tmp_362 = tmp_356 * tmp_37;
      real_t tmp_363 = tmp_359 * tmp_39;
      real_t tmp_364 = tmp_19 * ( 0.031405749086161582 * tmp_43 + 0.031405749086161582 * tmp_44 + tmp_45 );
      real_t tmp_365 = tmp_364 * tmp_41;
      real_t tmp_366 = tmp_364 * tmp_48;
      real_t tmp_367 = tmp_364 * tmp_50;
      real_t tmp_368 = -tmp_357 - tmp_358 - tmp_360 - tmp_361 - tmp_362 - tmp_363 - tmp_365 - tmp_366 - tmp_367 + 1;
      real_t tmp_369 = tmp_362 + tmp_363 + tmp_367;
      real_t tmp_370 = tmp_358 + tmp_361 + tmp_366;
      real_t tmp_371 = tmp_357 + tmp_360 + tmp_365;
      real_t tmp_372 =
          tmp_64 * ( tmp_10 * ( tmp_371 - 1.0 / 4.0 ) + tmp_3 * ( tmp_370 - 1.0 / 4.0 ) + tmp_5 * ( tmp_369 - 1.0 / 4.0 ) );
      real_t tmp_373 = 0.0068572537431980923 * tmp_66;
      real_t tmp_374 = tmp_19 * ( 0.19601935860219369 * tmp_21 + 0.19601935860219369 * tmp_22 + tmp_23 );
      real_t tmp_375 = tmp_374 * tmp_6;
      real_t tmp_376 = tmp_26 * tmp_374;
      real_t tmp_377 = tmp_19 * ( 0.19601935860219369 * tmp_30 + 0.19601935860219369 * tmp_31 + tmp_32 );
      real_t tmp_378 = tmp_28 * tmp_377;
      real_t tmp_379 = tmp_35 * tmp_377;
      real_t tmp_380 = tmp_37 * tmp_374;
      real_t tmp_381 = tmp_377 * tmp_39;
      real_t tmp_382 = tmp_19 * ( 0.19601935860219369 * tmp_43 + 0.19601935860219369 * tmp_44 + tmp_45 );
      real_t tmp_383 = tmp_382 * tmp_41;
      real_t tmp_384 = tmp_382 * tmp_48;
      real_t tmp_385 = tmp_382 * tmp_50;
      real_t tmp_386 = -tmp_375 - tmp_376 - tmp_378 - tmp_379 - tmp_380 - tmp_381 - tmp_383 - tmp_384 - tmp_385 + 1;
      real_t tmp_387 = tmp_380 + tmp_381 + tmp_385;
      real_t tmp_388 = tmp_376 + tmp_379 + tmp_384;
      real_t tmp_389 = tmp_375 + tmp_378 + tmp_383;
      real_t tmp_390 =
          tmp_64 * ( tmp_10 * ( tmp_389 - 1.0 / 4.0 ) + tmp_3 * ( tmp_388 - 1.0 / 4.0 ) + tmp_5 * ( tmp_387 - 1.0 / 4.0 ) );
      real_t tmp_391 = 0.037198804536718075 * tmp_66;
      real_t tmp_392 = tmp_19 * ( 0.40446199974765351 * tmp_21 + 0.40446199974765351 * tmp_22 + tmp_23 );
      real_t tmp_393 = tmp_392 * tmp_6;
      real_t tmp_394 = tmp_26 * tmp_392;
      real_t tmp_395 = tmp_19 * ( 0.40446199974765351 * tmp_30 + 0.40446199974765351 * tmp_31 + tmp_32 );
      real_t tmp_396 = tmp_28 * tmp_395;
      real_t tmp_397 = tmp_35 * tmp_395;
      real_t tmp_398 = tmp_37 * tmp_392;
      real_t tmp_399 = tmp_39 * tmp_395;
      real_t tmp_400 = tmp_19 * ( 0.40446199974765351 * tmp_43 + 0.40446199974765351 * tmp_44 + tmp_45 );
      real_t tmp_401 = tmp_400 * tmp_41;
      real_t tmp_402 = tmp_400 * tmp_48;
      real_t tmp_403 = tmp_400 * tmp_50;
      real_t tmp_404 = -tmp_393 - tmp_394 - tmp_396 - tmp_397 - tmp_398 - tmp_399 - tmp_401 - tmp_402 - tmp_403 + 1;
      real_t tmp_405 = tmp_398 + tmp_399 + tmp_403;
      real_t tmp_406 = tmp_394 + tmp_397 + tmp_402;
      real_t tmp_407 = tmp_393 + tmp_396 + tmp_401;
      real_t tmp_408 =
          tmp_64 * ( tmp_10 * ( tmp_407 - 1.0 / 4.0 ) + tmp_3 * ( tmp_406 - 1.0 / 4.0 ) + tmp_5 * ( tmp_405 - 1.0 / 4.0 ) );
      real_t tmp_409 = 0.042507265838595799 * tmp_66;
      real_t tmp_410 = tmp_19 * ( 0.1711304259088916 * tmp_21 + 0.041227165399737475 * tmp_22 + tmp_23 );
      real_t tmp_411 = tmp_410 * tmp_6;
      real_t tmp_412 = tmp_26 * tmp_410;
      real_t tmp_413 = tmp_19 * ( 0.1711304259088916 * tmp_30 + 0.041227165399737475 * tmp_31 + tmp_32 );
      real_t tmp_414 = tmp_28 * tmp_413;
      real_t tmp_415 = tmp_35 * tmp_413;
      real_t tmp_416 = tmp_37 * tmp_410;
      real_t tmp_417 = tmp_39 * tmp_413;
      real_t tmp_418 = tmp_19 * ( 0.1711304259088916 * tmp_43 + 0.041227165399737475 * tmp_44 + tmp_45 );
      real_t tmp_419 = tmp_41 * tmp_418;
      real_t tmp_420 = tmp_418 * tmp_48;
      real_t tmp_421 = tmp_418 * tmp_50;
      real_t tmp_422 = -tmp_411 - tmp_412 - tmp_414 - tmp_415 - tmp_416 - tmp_417 - tmp_419 - tmp_420 - tmp_421 + 1;
      real_t tmp_423 = tmp_416 + tmp_417 + tmp_421;
      real_t tmp_424 = tmp_412 + tmp_415 + tmp_420;
      real_t tmp_425 = tmp_411 + tmp_414 + tmp_419;
      real_t tmp_426 =
          tmp_64 * ( tmp_10 * ( tmp_425 - 1.0 / 4.0 ) + tmp_3 * ( tmp_424 - 1.0 / 4.0 ) + tmp_5 * ( tmp_423 - 1.0 / 4.0 ) );
      real_t tmp_427 = 0.019202922745021479 * tmp_66;
      real_t a_0_0   = tmp_103 * ( tmp_102 * tmp_98 - tmp_56 * tmp_98 ) + tmp_121 * ( tmp_116 * tmp_120 - tmp_116 * tmp_56 ) +
                     tmp_139 * ( tmp_134 * tmp_138 - tmp_134 * tmp_56 ) + tmp_157 * ( tmp_152 * tmp_156 - tmp_152 * tmp_56 ) +
                     tmp_175 * ( tmp_170 * tmp_174 - tmp_170 * tmp_56 ) + tmp_193 * ( tmp_188 * tmp_192 - tmp_188 * tmp_56 ) +
                     tmp_211 * ( tmp_206 * tmp_210 - tmp_206 * tmp_56 ) + tmp_229 * ( tmp_224 * tmp_228 - tmp_224 * tmp_56 ) +
                     tmp_247 * ( tmp_242 * tmp_246 - tmp_242 * tmp_56 ) + tmp_265 * ( tmp_260 * tmp_264 - tmp_260 * tmp_56 ) +
                     tmp_283 * ( tmp_278 * tmp_282 - tmp_278 * tmp_56 ) + tmp_301 * ( tmp_296 * tmp_300 - tmp_296 * tmp_56 ) +
                     tmp_319 * ( tmp_314 * tmp_318 - tmp_314 * tmp_56 ) + tmp_337 * ( tmp_332 * tmp_336 - tmp_332 * tmp_56 ) +
                     tmp_355 * ( tmp_350 * tmp_354 - tmp_350 * tmp_56 ) + tmp_373 * ( tmp_368 * tmp_372 - tmp_368 * tmp_56 ) +
                     tmp_391 * ( tmp_386 * tmp_390 - tmp_386 * tmp_56 ) + tmp_409 * ( tmp_404 * tmp_408 - tmp_404 * tmp_56 ) +
                     tmp_427 * ( tmp_422 * tmp_426 - tmp_422 * tmp_56 ) + tmp_67 * ( -tmp_52 * tmp_56 + tmp_52 * tmp_65 ) +
                     tmp_85 * ( -tmp_56 * tmp_80 + tmp_80 * tmp_84 );
      real_t a_1_0 = tmp_103 * ( tmp_102 * tmp_99 - tmp_56 * tmp_99 ) + tmp_121 * ( tmp_117 * tmp_120 - tmp_117 * tmp_56 ) +
                     tmp_139 * ( tmp_135 * tmp_138 - tmp_135 * tmp_56 ) + tmp_157 * ( tmp_153 * tmp_156 - tmp_153 * tmp_56 ) +
                     tmp_175 * ( tmp_171 * tmp_174 - tmp_171 * tmp_56 ) + tmp_193 * ( tmp_189 * tmp_192 - tmp_189 * tmp_56 ) +
                     tmp_211 * ( tmp_207 * tmp_210 - tmp_207 * tmp_56 ) + tmp_229 * ( tmp_225 * tmp_228 - tmp_225 * tmp_56 ) +
                     tmp_247 * ( tmp_243 * tmp_246 - tmp_243 * tmp_56 ) + tmp_265 * ( tmp_261 * tmp_264 - tmp_261 * tmp_56 ) +
                     tmp_283 * ( tmp_279 * tmp_282 - tmp_279 * tmp_56 ) + tmp_301 * ( tmp_297 * tmp_300 - tmp_297 * tmp_56 ) +
                     tmp_319 * ( tmp_315 * tmp_318 - tmp_315 * tmp_56 ) + tmp_337 * ( tmp_333 * tmp_336 - tmp_333 * tmp_56 ) +
                     tmp_355 * ( tmp_351 * tmp_354 - tmp_351 * tmp_56 ) + tmp_373 * ( tmp_369 * tmp_372 - tmp_369 * tmp_56 ) +
                     tmp_391 * ( tmp_387 * tmp_390 - tmp_387 * tmp_56 ) + tmp_409 * ( tmp_405 * tmp_408 - tmp_405 * tmp_56 ) +
                     tmp_427 * ( tmp_423 * tmp_426 - tmp_423 * tmp_56 ) + tmp_67 * ( -tmp_56 * tmp_57 + tmp_57 * tmp_65 ) +
                     tmp_85 * ( -tmp_56 * tmp_81 + tmp_81 * tmp_84 );
      real_t a_2_0 = tmp_103 * ( tmp_100 * tmp_102 - tmp_100 * tmp_56 ) + tmp_121 * ( tmp_118 * tmp_120 - tmp_118 * tmp_56 ) +
                     tmp_139 * ( tmp_136 * tmp_138 - tmp_136 * tmp_56 ) + tmp_157 * ( tmp_154 * tmp_156 - tmp_154 * tmp_56 ) +
                     tmp_175 * ( tmp_172 * tmp_174 - tmp_172 * tmp_56 ) + tmp_193 * ( tmp_190 * tmp_192 - tmp_190 * tmp_56 ) +
                     tmp_211 * ( tmp_208 * tmp_210 - tmp_208 * tmp_56 ) + tmp_229 * ( tmp_226 * tmp_228 - tmp_226 * tmp_56 ) +
                     tmp_247 * ( tmp_244 * tmp_246 - tmp_244 * tmp_56 ) + tmp_265 * ( tmp_262 * tmp_264 - tmp_262 * tmp_56 ) +
                     tmp_283 * ( tmp_280 * tmp_282 - tmp_280 * tmp_56 ) + tmp_301 * ( tmp_298 * tmp_300 - tmp_298 * tmp_56 ) +
                     tmp_319 * ( tmp_316 * tmp_318 - tmp_316 * tmp_56 ) + tmp_337 * ( tmp_334 * tmp_336 - tmp_334 * tmp_56 ) +
                     tmp_355 * ( tmp_352 * tmp_354 - tmp_352 * tmp_56 ) + tmp_373 * ( tmp_370 * tmp_372 - tmp_370 * tmp_56 ) +
                     tmp_391 * ( tmp_388 * tmp_390 - tmp_388 * tmp_56 ) + tmp_409 * ( tmp_406 * tmp_408 - tmp_406 * tmp_56 ) +
                     tmp_427 * ( tmp_424 * tmp_426 - tmp_424 * tmp_56 ) + tmp_67 * ( -tmp_56 * tmp_58 + tmp_58 * tmp_65 ) +
                     tmp_85 * ( -tmp_56 * tmp_82 + tmp_82 * tmp_84 );
      real_t a_3_0 = tmp_103 * ( tmp_101 * tmp_102 - tmp_101 * tmp_56 ) + tmp_121 * ( tmp_119 * tmp_120 - tmp_119 * tmp_56 ) +
                     tmp_139 * ( tmp_137 * tmp_138 - tmp_137 * tmp_56 ) + tmp_157 * ( tmp_155 * tmp_156 - tmp_155 * tmp_56 ) +
                     tmp_175 * ( tmp_173 * tmp_174 - tmp_173 * tmp_56 ) + tmp_193 * ( tmp_191 * tmp_192 - tmp_191 * tmp_56 ) +
                     tmp_211 * ( tmp_209 * tmp_210 - tmp_209 * tmp_56 ) + tmp_229 * ( tmp_227 * tmp_228 - tmp_227 * tmp_56 ) +
                     tmp_247 * ( tmp_245 * tmp_246 - tmp_245 * tmp_56 ) + tmp_265 * ( tmp_263 * tmp_264 - tmp_263 * tmp_56 ) +
                     tmp_283 * ( tmp_281 * tmp_282 - tmp_281 * tmp_56 ) + tmp_301 * ( tmp_299 * tmp_300 - tmp_299 * tmp_56 ) +
                     tmp_319 * ( tmp_317 * tmp_318 - tmp_317 * tmp_56 ) + tmp_337 * ( tmp_335 * tmp_336 - tmp_335 * tmp_56 ) +
                     tmp_355 * ( tmp_353 * tmp_354 - tmp_353 * tmp_56 ) + tmp_373 * ( tmp_371 * tmp_372 - tmp_371 * tmp_56 ) +
                     tmp_391 * ( tmp_389 * tmp_390 - tmp_389 * tmp_56 ) + tmp_409 * ( tmp_407 * tmp_408 - tmp_407 * tmp_56 ) +
                     tmp_427 * ( tmp_425 * tmp_426 - tmp_425 * tmp_56 ) + tmp_67 * ( -tmp_56 * tmp_59 + tmp_59 * tmp_65 ) +
                     tmp_85 * ( -tmp_56 * tmp_83 + tmp_83 * tmp_84 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
      elMat( 3, 0 ) = a_3_0;
   }

   void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                  const std::vector< Point3D >& coordsElementOuter,
                                  const std::vector< Point3D >& coordsFacet,
                                  const Point3D&,
                                  const Point3D&,
                                  const Point3D&     outwardNormal,
                                  const DGBasisInfo& trialBasis,
                                  const DGBasisInfo& testBasis,
                                  int                trialDegree,
                                  int                testDegree,
                                  MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_0;
      real_t tmp_6  = p_affine_1_1 + tmp_2;
      real_t tmp_7  = tmp_5 * tmp_6;
      real_t tmp_8  = tmp_4 - tmp_7;
      real_t tmp_9  = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = p_affine_1_2 + tmp_9;
      real_t tmp_13 = tmp_12 * tmp_5;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = p_affine_2_2 + tmp_9;
      real_t tmp_16 = tmp_15 * tmp_6;
      real_t tmp_17 = tmp_1 * tmp_11;
      real_t tmp_18 = tmp_12 * tmp_14;
      real_t tmp_19 =
          1.0 / ( tmp_10 * tmp_4 - tmp_10 * tmp_7 + tmp_11 * tmp_13 + tmp_14 * tmp_16 - tmp_15 * tmp_17 - tmp_18 * tmp_3 );
      real_t tmp_20 = p_affine_8_2 + tmp_9;
      real_t tmp_21 = -p_affine_8_2;
      real_t tmp_22 = p_affine_9_2 + tmp_21;
      real_t tmp_23 = p_affine_10_2 + tmp_21;
      real_t tmp_24 = 0.031405749086161582 * tmp_22 + 0.93718850182767688 * tmp_23;
      real_t tmp_25 = tmp_19 * ( tmp_20 + tmp_24 );
      real_t tmp_26 = tmp_25 * tmp_8;
      real_t tmp_27 = tmp_14 * tmp_6 - tmp_17;
      real_t tmp_28 = tmp_25 * tmp_27;
      real_t tmp_29 = -tmp_1 * tmp_15 + tmp_13;
      real_t tmp_30 = p_affine_8_1 + tmp_2;
      real_t tmp_31 = -p_affine_8_1;
      real_t tmp_32 = p_affine_9_1 + tmp_31;
      real_t tmp_33 = p_affine_10_1 + tmp_31;
      real_t tmp_34 = 0.031405749086161582 * tmp_32 + 0.93718850182767688 * tmp_33;
      real_t tmp_35 = tmp_19 * ( tmp_30 + tmp_34 );
      real_t tmp_36 = tmp_29 * tmp_35;
      real_t tmp_37 = tmp_1 * tmp_10 - tmp_18;
      real_t tmp_38 = tmp_35 * tmp_37;
      real_t tmp_39 = tmp_11 * tmp_5 - tmp_14 * tmp_3;
      real_t tmp_40 = tmp_25 * tmp_39;
      real_t tmp_41 = -tmp_10 * tmp_5 + tmp_14 * tmp_15;
      real_t tmp_42 = tmp_35 * tmp_41;
      real_t tmp_43 = -tmp_12 * tmp_3 + tmp_16;
      real_t tmp_44 = p_affine_8_0 + tmp_0;
      real_t tmp_45 = -p_affine_8_0;
      real_t tmp_46 = p_affine_9_0 + tmp_45;
      real_t tmp_47 = p_affine_10_0 + tmp_45;
      real_t tmp_48 = 0.031405749086161582 * tmp_46 + 0.93718850182767688 * tmp_47;
      real_t tmp_49 = tmp_19 * ( tmp_44 + tmp_48 );
      real_t tmp_50 = tmp_43 * tmp_49;
      real_t tmp_51 = -tmp_10 * tmp_6 + tmp_11 * tmp_12;
      real_t tmp_52 = tmp_49 * tmp_51;
      real_t tmp_53 = tmp_10 * tmp_3 - tmp_11 * tmp_15;
      real_t tmp_54 = tmp_49 * tmp_53;
      real_t tmp_55 = -tmp_26 - tmp_28 - tmp_36 - tmp_38 - tmp_40 - tmp_42 - tmp_50 - tmp_52 - tmp_54 + 1;
      real_t tmp_56 = -p_affine_4_1;
      real_t tmp_57 = p_affine_6_1 + tmp_56;
      real_t tmp_58 = -p_affine_4_2;
      real_t tmp_59 = p_affine_7_2 + tmp_58;
      real_t tmp_60 = tmp_57 * tmp_59;
      real_t tmp_61 = p_affine_7_1 + tmp_56;
      real_t tmp_62 = p_affine_6_2 + tmp_58;
      real_t tmp_63 = tmp_61 * tmp_62;
      real_t tmp_64 = tmp_60 - tmp_63;
      real_t tmp_65 = p_affine_5_1 + tmp_56;
      real_t tmp_66 = -p_affine_4_0;
      real_t tmp_67 = p_affine_5_0 + tmp_66;
      real_t tmp_68 = p_affine_6_0 + tmp_66;
      real_t tmp_69 = p_affine_5_2 + tmp_58;
      real_t tmp_70 = tmp_61 * tmp_69;
      real_t tmp_71 = p_affine_7_0 + tmp_66;
      real_t tmp_72 = tmp_62 * tmp_65;
      real_t tmp_73 = tmp_59 * tmp_65;
      real_t tmp_74 = tmp_57 * tmp_69;
      real_t tmp_75 =
          1.0 / ( tmp_60 * tmp_67 - tmp_63 * tmp_67 + tmp_68 * tmp_70 - tmp_68 * tmp_73 + tmp_71 * tmp_72 - tmp_71 * tmp_74 );
      real_t tmp_76 = tmp_65 * tmp_75;
      real_t tmp_77 = tmp_70 - tmp_73;
      real_t tmp_78 = tmp_57 * tmp_75;
      real_t tmp_79 = tmp_72 - tmp_74;
      real_t tmp_80 = tmp_61 * tmp_75;
      real_t tmp_81 = -tmp_59 * tmp_68 + tmp_62 * tmp_71;
      real_t tmp_82 = tmp_59 * tmp_67 - tmp_69 * tmp_71;
      real_t tmp_83 = -tmp_62 * tmp_67 + tmp_68 * tmp_69;
      real_t tmp_84 = -tmp_57 * tmp_71 + tmp_61 * tmp_68;
      real_t tmp_85 = -tmp_61 * tmp_67 + tmp_65 * tmp_71;
      real_t tmp_86 = tmp_57 * tmp_67 - tmp_65 * tmp_68;
      real_t tmp_87 = 0.5 * p_affine_13_0 * ( tmp_64 * tmp_76 + tmp_77 * tmp_78 + tmp_79 * tmp_80 ) +
                      0.5 * p_affine_13_1 * ( tmp_76 * tmp_81 + tmp_78 * tmp_82 + tmp_80 * tmp_83 ) +
                      0.5 * p_affine_13_2 * ( tmp_76 * tmp_84 + tmp_78 * tmp_85 + tmp_80 * tmp_86 );
      real_t tmp_88 = p_affine_8_2 + tmp_58;
      real_t tmp_89 = tmp_75 * ( tmp_24 + tmp_88 );
      real_t tmp_90 = p_affine_8_1 + tmp_56;
      real_t tmp_91 = tmp_75 * ( tmp_34 + tmp_90 );
      real_t tmp_92 = p_affine_8_0 + tmp_66;
      real_t tmp_93 = tmp_75 * ( tmp_48 + tmp_92 );
      real_t tmp_94 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_95 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_96 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_97 = ( std::abs( tmp_23 * tmp_94 - tmp_33 * tmp_96 ) * std::abs( tmp_23 * tmp_94 - tmp_33 * tmp_96 ) ) +
                      ( std::abs( tmp_23 * tmp_95 - tmp_47 * tmp_96 ) * std::abs( tmp_23 * tmp_95 - tmp_47 * tmp_96 ) ) +
                      ( std::abs( tmp_33 * tmp_95 - tmp_47 * tmp_94 ) * std::abs( tmp_33 * tmp_95 - tmp_47 * tmp_94 ) );
      real_t tmp_98  = 1.0 * std::pow( tmp_97, -0.25 );
      real_t tmp_99  = tmp_98 * ( tmp_57 * ( tmp_77 * tmp_93 + tmp_82 * tmp_91 + tmp_85 * tmp_89 - 1.0 / 4.0 ) +
                                 tmp_61 * ( tmp_79 * tmp_93 + tmp_83 * tmp_91 + tmp_86 * tmp_89 - 1.0 / 4.0 ) +
                                 tmp_65 * ( tmp_64 * tmp_93 + tmp_81 * tmp_91 + tmp_84 * tmp_89 - 1.0 / 4.0 ) );
      real_t tmp_100 = 1.0 * std::pow( tmp_97, 1.0 / 2.0 );
      real_t tmp_101 = 0.0068572537431980923 * tmp_100;
      real_t tmp_102 = 0.19601935860219369 * tmp_22 + 0.60796128279561268 * tmp_23;
      real_t tmp_103 = tmp_19 * ( tmp_102 + tmp_20 );
      real_t tmp_104 = tmp_103 * tmp_8;
      real_t tmp_105 = tmp_103 * tmp_27;
      real_t tmp_106 = 0.19601935860219369 * tmp_32 + 0.60796128279561268 * tmp_33;
      real_t tmp_107 = tmp_19 * ( tmp_106 + tmp_30 );
      real_t tmp_108 = tmp_107 * tmp_29;
      real_t tmp_109 = tmp_107 * tmp_37;
      real_t tmp_110 = tmp_103 * tmp_39;
      real_t tmp_111 = tmp_107 * tmp_41;
      real_t tmp_112 = 0.19601935860219369 * tmp_46 + 0.60796128279561268 * tmp_47;
      real_t tmp_113 = tmp_19 * ( tmp_112 + tmp_44 );
      real_t tmp_114 = tmp_113 * tmp_43;
      real_t tmp_115 = tmp_113 * tmp_51;
      real_t tmp_116 = tmp_113 * tmp_53;
      real_t tmp_117 = -tmp_104 - tmp_105 - tmp_108 - tmp_109 - tmp_110 - tmp_111 - tmp_114 - tmp_115 - tmp_116 + 1;
      real_t tmp_118 = tmp_75 * ( tmp_102 + tmp_88 );
      real_t tmp_119 = tmp_75 * ( tmp_106 + tmp_90 );
      real_t tmp_120 = tmp_75 * ( tmp_112 + tmp_92 );
      real_t tmp_121 = tmp_98 * ( tmp_57 * ( tmp_118 * tmp_85 + tmp_119 * tmp_82 + tmp_120 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_118 * tmp_86 + tmp_119 * tmp_83 + tmp_120 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_118 * tmp_84 + tmp_119 * tmp_81 + tmp_120 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_122 = 0.037198804536718075 * tmp_100;
      real_t tmp_123 = 0.37605877282253791 * tmp_22 + 0.039308471900058539 * tmp_23;
      real_t tmp_124 = tmp_19 * ( tmp_123 + tmp_20 );
      real_t tmp_125 = tmp_124 * tmp_8;
      real_t tmp_126 = tmp_124 * tmp_27;
      real_t tmp_127 = 0.37605877282253791 * tmp_32 + 0.039308471900058539 * tmp_33;
      real_t tmp_128 = tmp_19 * ( tmp_127 + tmp_30 );
      real_t tmp_129 = tmp_128 * tmp_29;
      real_t tmp_130 = tmp_128 * tmp_37;
      real_t tmp_131 = tmp_124 * tmp_39;
      real_t tmp_132 = tmp_128 * tmp_41;
      real_t tmp_133 = 0.37605877282253791 * tmp_46 + 0.039308471900058539 * tmp_47;
      real_t tmp_134 = tmp_19 * ( tmp_133 + tmp_44 );
      real_t tmp_135 = tmp_134 * tmp_43;
      real_t tmp_136 = tmp_134 * tmp_51;
      real_t tmp_137 = tmp_134 * tmp_53;
      real_t tmp_138 = -tmp_125 - tmp_126 - tmp_129 - tmp_130 - tmp_131 - tmp_132 - tmp_135 - tmp_136 - tmp_137 + 1;
      real_t tmp_139 = tmp_75 * ( tmp_123 + tmp_88 );
      real_t tmp_140 = tmp_75 * ( tmp_127 + tmp_90 );
      real_t tmp_141 = tmp_75 * ( tmp_133 + tmp_92 );
      real_t tmp_142 = tmp_98 * ( tmp_57 * ( tmp_139 * tmp_85 + tmp_140 * tmp_82 + tmp_141 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_139 * tmp_86 + tmp_140 * tmp_83 + tmp_141 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_139 * tmp_84 + tmp_140 * tmp_81 + tmp_141 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_143 = 0.020848748529055869 * tmp_100;
      real_t tmp_144 = 0.78764240869137092 * tmp_22 + 0.1711304259088916 * tmp_23;
      real_t tmp_145 = tmp_19 * ( tmp_144 + tmp_20 );
      real_t tmp_146 = tmp_145 * tmp_8;
      real_t tmp_147 = tmp_145 * tmp_27;
      real_t tmp_148 = 0.78764240869137092 * tmp_32 + 0.1711304259088916 * tmp_33;
      real_t tmp_149 = tmp_19 * ( tmp_148 + tmp_30 );
      real_t tmp_150 = tmp_149 * tmp_29;
      real_t tmp_151 = tmp_149 * tmp_37;
      real_t tmp_152 = tmp_145 * tmp_39;
      real_t tmp_153 = tmp_149 * tmp_41;
      real_t tmp_154 = 0.78764240869137092 * tmp_46 + 0.1711304259088916 * tmp_47;
      real_t tmp_155 = tmp_19 * ( tmp_154 + tmp_44 );
      real_t tmp_156 = tmp_155 * tmp_43;
      real_t tmp_157 = tmp_155 * tmp_51;
      real_t tmp_158 = tmp_155 * tmp_53;
      real_t tmp_159 = -tmp_146 - tmp_147 - tmp_150 - tmp_151 - tmp_152 - tmp_153 - tmp_156 - tmp_157 - tmp_158 + 1;
      real_t tmp_160 = tmp_75 * ( tmp_144 + tmp_88 );
      real_t tmp_161 = tmp_75 * ( tmp_148 + tmp_90 );
      real_t tmp_162 = tmp_75 * ( tmp_154 + tmp_92 );
      real_t tmp_163 = tmp_98 * ( tmp_57 * ( tmp_160 * tmp_85 + tmp_161 * tmp_82 + tmp_162 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_160 * tmp_86 + tmp_161 * tmp_83 + tmp_162 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_160 * tmp_84 + tmp_161 * tmp_81 + tmp_162 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_164 = 0.019202922745021479 * tmp_100;
      real_t tmp_165 = 0.58463275527740355 * tmp_22 + 0.37605877282253791 * tmp_23;
      real_t tmp_166 = tmp_19 * ( tmp_165 + tmp_20 );
      real_t tmp_167 = tmp_166 * tmp_8;
      real_t tmp_168 = tmp_166 * tmp_27;
      real_t tmp_169 = 0.58463275527740355 * tmp_32 + 0.37605877282253791 * tmp_33;
      real_t tmp_170 = tmp_19 * ( tmp_169 + tmp_30 );
      real_t tmp_171 = tmp_170 * tmp_29;
      real_t tmp_172 = tmp_170 * tmp_37;
      real_t tmp_173 = tmp_166 * tmp_39;
      real_t tmp_174 = tmp_170 * tmp_41;
      real_t tmp_175 = 0.58463275527740355 * tmp_46 + 0.37605877282253791 * tmp_47;
      real_t tmp_176 = tmp_19 * ( tmp_175 + tmp_44 );
      real_t tmp_177 = tmp_176 * tmp_43;
      real_t tmp_178 = tmp_176 * tmp_51;
      real_t tmp_179 = tmp_176 * tmp_53;
      real_t tmp_180 = -tmp_167 - tmp_168 - tmp_171 - tmp_172 - tmp_173 - tmp_174 - tmp_177 - tmp_178 - tmp_179 + 1;
      real_t tmp_181 = tmp_75 * ( tmp_165 + tmp_88 );
      real_t tmp_182 = tmp_75 * ( tmp_169 + tmp_90 );
      real_t tmp_183 = tmp_75 * ( tmp_175 + tmp_92 );
      real_t tmp_184 = tmp_98 * ( tmp_57 * ( tmp_181 * tmp_85 + tmp_182 * tmp_82 + tmp_183 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_181 * tmp_86 + tmp_182 * tmp_83 + tmp_183 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_181 * tmp_84 + tmp_182 * tmp_81 + tmp_183 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_185 = 0.020848748529055869 * tmp_100;
      real_t tmp_186 = 0.041227165399737475 * tmp_22 + 0.78764240869137092 * tmp_23;
      real_t tmp_187 = tmp_19 * ( tmp_186 + tmp_20 );
      real_t tmp_188 = tmp_187 * tmp_8;
      real_t tmp_189 = tmp_187 * tmp_27;
      real_t tmp_190 = 0.041227165399737475 * tmp_32 + 0.78764240869137092 * tmp_33;
      real_t tmp_191 = tmp_19 * ( tmp_190 + tmp_30 );
      real_t tmp_192 = tmp_191 * tmp_29;
      real_t tmp_193 = tmp_191 * tmp_37;
      real_t tmp_194 = tmp_187 * tmp_39;
      real_t tmp_195 = tmp_191 * tmp_41;
      real_t tmp_196 = 0.041227165399737475 * tmp_46 + 0.78764240869137092 * tmp_47;
      real_t tmp_197 = tmp_19 * ( tmp_196 + tmp_44 );
      real_t tmp_198 = tmp_197 * tmp_43;
      real_t tmp_199 = tmp_197 * tmp_51;
      real_t tmp_200 = tmp_197 * tmp_53;
      real_t tmp_201 = -tmp_188 - tmp_189 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_198 - tmp_199 - tmp_200 + 1;
      real_t tmp_202 = tmp_75 * ( tmp_186 + tmp_88 );
      real_t tmp_203 = tmp_75 * ( tmp_190 + tmp_90 );
      real_t tmp_204 = tmp_75 * ( tmp_196 + tmp_92 );
      real_t tmp_205 = tmp_98 * ( tmp_57 * ( tmp_202 * tmp_85 + tmp_203 * tmp_82 + tmp_204 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_202 * tmp_86 + tmp_203 * tmp_83 + tmp_204 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_202 * tmp_84 + tmp_203 * tmp_81 + tmp_204 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_206 = 0.019202922745021479 * tmp_100;
      real_t tmp_207 = 0.039308471900058539 * tmp_22 + 0.58463275527740355 * tmp_23;
      real_t tmp_208 = tmp_19 * ( tmp_20 + tmp_207 );
      real_t tmp_209 = tmp_208 * tmp_8;
      real_t tmp_210 = tmp_208 * tmp_27;
      real_t tmp_211 = 0.039308471900058539 * tmp_32 + 0.58463275527740355 * tmp_33;
      real_t tmp_212 = tmp_19 * ( tmp_211 + tmp_30 );
      real_t tmp_213 = tmp_212 * tmp_29;
      real_t tmp_214 = tmp_212 * tmp_37;
      real_t tmp_215 = tmp_208 * tmp_39;
      real_t tmp_216 = tmp_212 * tmp_41;
      real_t tmp_217 = 0.039308471900058539 * tmp_46 + 0.58463275527740355 * tmp_47;
      real_t tmp_218 = tmp_19 * ( tmp_217 + tmp_44 );
      real_t tmp_219 = tmp_218 * tmp_43;
      real_t tmp_220 = tmp_218 * tmp_51;
      real_t tmp_221 = tmp_218 * tmp_53;
      real_t tmp_222 = -tmp_209 - tmp_210 - tmp_213 - tmp_214 - tmp_215 - tmp_216 - tmp_219 - tmp_220 - tmp_221 + 1;
      real_t tmp_223 = tmp_75 * ( tmp_207 + tmp_88 );
      real_t tmp_224 = tmp_75 * ( tmp_211 + tmp_90 );
      real_t tmp_225 = tmp_75 * ( tmp_217 + tmp_92 );
      real_t tmp_226 = tmp_98 * ( tmp_57 * ( tmp_223 * tmp_85 + tmp_224 * tmp_82 + tmp_225 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_223 * tmp_86 + tmp_224 * tmp_83 + tmp_225 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_223 * tmp_84 + tmp_224 * tmp_81 + tmp_225 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_227 = 0.020848748529055869 * tmp_100;
      real_t tmp_228 = 0.78764240869137092 * tmp_22 + 0.041227165399737475 * tmp_23;
      real_t tmp_229 = tmp_19 * ( tmp_20 + tmp_228 );
      real_t tmp_230 = tmp_229 * tmp_8;
      real_t tmp_231 = tmp_229 * tmp_27;
      real_t tmp_232 = 0.78764240869137092 * tmp_32 + 0.041227165399737475 * tmp_33;
      real_t tmp_233 = tmp_19 * ( tmp_232 + tmp_30 );
      real_t tmp_234 = tmp_233 * tmp_29;
      real_t tmp_235 = tmp_233 * tmp_37;
      real_t tmp_236 = tmp_229 * tmp_39;
      real_t tmp_237 = tmp_233 * tmp_41;
      real_t tmp_238 = 0.78764240869137092 * tmp_46 + 0.041227165399737475 * tmp_47;
      real_t tmp_239 = tmp_19 * ( tmp_238 + tmp_44 );
      real_t tmp_240 = tmp_239 * tmp_43;
      real_t tmp_241 = tmp_239 * tmp_51;
      real_t tmp_242 = tmp_239 * tmp_53;
      real_t tmp_243 = -tmp_230 - tmp_231 - tmp_234 - tmp_235 - tmp_236 - tmp_237 - tmp_240 - tmp_241 - tmp_242 + 1;
      real_t tmp_244 = tmp_75 * ( tmp_228 + tmp_88 );
      real_t tmp_245 = tmp_75 * ( tmp_232 + tmp_90 );
      real_t tmp_246 = tmp_75 * ( tmp_238 + tmp_92 );
      real_t tmp_247 = tmp_98 * ( tmp_57 * ( tmp_244 * tmp_85 + tmp_245 * tmp_82 + tmp_246 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_244 * tmp_86 + tmp_245 * tmp_83 + tmp_246 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_244 * tmp_84 + tmp_245 * tmp_81 + tmp_246 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_248 = 0.019202922745021479 * tmp_100;
      real_t tmp_249 = 0.58463275527740355 * tmp_22 + 0.039308471900058539 * tmp_23;
      real_t tmp_250 = tmp_19 * ( tmp_20 + tmp_249 );
      real_t tmp_251 = tmp_250 * tmp_8;
      real_t tmp_252 = tmp_250 * tmp_27;
      real_t tmp_253 = 0.58463275527740355 * tmp_32 + 0.039308471900058539 * tmp_33;
      real_t tmp_254 = tmp_19 * ( tmp_253 + tmp_30 );
      real_t tmp_255 = tmp_254 * tmp_29;
      real_t tmp_256 = tmp_254 * tmp_37;
      real_t tmp_257 = tmp_250 * tmp_39;
      real_t tmp_258 = tmp_254 * tmp_41;
      real_t tmp_259 = 0.58463275527740355 * tmp_46 + 0.039308471900058539 * tmp_47;
      real_t tmp_260 = tmp_19 * ( tmp_259 + tmp_44 );
      real_t tmp_261 = tmp_260 * tmp_43;
      real_t tmp_262 = tmp_260 * tmp_51;
      real_t tmp_263 = tmp_260 * tmp_53;
      real_t tmp_264 = -tmp_251 - tmp_252 - tmp_255 - tmp_256 - tmp_257 - tmp_258 - tmp_261 - tmp_262 - tmp_263 + 1;
      real_t tmp_265 = tmp_75 * ( tmp_249 + tmp_88 );
      real_t tmp_266 = tmp_75 * ( tmp_253 + tmp_90 );
      real_t tmp_267 = tmp_75 * ( tmp_259 + tmp_92 );
      real_t tmp_268 = tmp_98 * ( tmp_57 * ( tmp_265 * tmp_85 + tmp_266 * tmp_82 + tmp_267 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_265 * tmp_86 + tmp_266 * tmp_83 + tmp_267 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_265 * tmp_84 + tmp_266 * tmp_81 + tmp_267 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_269 = 0.020848748529055869 * tmp_100;
      real_t tmp_270 = 0.1711304259088916 * tmp_22 + 0.78764240869137092 * tmp_23;
      real_t tmp_271 = tmp_19 * ( tmp_20 + tmp_270 );
      real_t tmp_272 = tmp_271 * tmp_8;
      real_t tmp_273 = tmp_27 * tmp_271;
      real_t tmp_274 = 0.1711304259088916 * tmp_32 + 0.78764240869137092 * tmp_33;
      real_t tmp_275 = tmp_19 * ( tmp_274 + tmp_30 );
      real_t tmp_276 = tmp_275 * tmp_29;
      real_t tmp_277 = tmp_275 * tmp_37;
      real_t tmp_278 = tmp_271 * tmp_39;
      real_t tmp_279 = tmp_275 * tmp_41;
      real_t tmp_280 = 0.1711304259088916 * tmp_46 + 0.78764240869137092 * tmp_47;
      real_t tmp_281 = tmp_19 * ( tmp_280 + tmp_44 );
      real_t tmp_282 = tmp_281 * tmp_43;
      real_t tmp_283 = tmp_281 * tmp_51;
      real_t tmp_284 = tmp_281 * tmp_53;
      real_t tmp_285 = -tmp_272 - tmp_273 - tmp_276 - tmp_277 - tmp_278 - tmp_279 - tmp_282 - tmp_283 - tmp_284 + 1;
      real_t tmp_286 = tmp_75 * ( tmp_270 + tmp_88 );
      real_t tmp_287 = tmp_75 * ( tmp_274 + tmp_90 );
      real_t tmp_288 = tmp_75 * ( tmp_280 + tmp_92 );
      real_t tmp_289 = tmp_98 * ( tmp_57 * ( tmp_286 * tmp_85 + tmp_287 * tmp_82 + tmp_288 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_286 * tmp_86 + tmp_287 * tmp_83 + tmp_288 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_286 * tmp_84 + tmp_287 * tmp_81 + tmp_288 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_290 = 0.019202922745021479 * tmp_100;
      real_t tmp_291 = 0.37605877282253791 * tmp_22 + 0.58463275527740355 * tmp_23;
      real_t tmp_292 = tmp_19 * ( tmp_20 + tmp_291 );
      real_t tmp_293 = tmp_292 * tmp_8;
      real_t tmp_294 = tmp_27 * tmp_292;
      real_t tmp_295 = 0.37605877282253791 * tmp_32 + 0.58463275527740355 * tmp_33;
      real_t tmp_296 = tmp_19 * ( tmp_295 + tmp_30 );
      real_t tmp_297 = tmp_29 * tmp_296;
      real_t tmp_298 = tmp_296 * tmp_37;
      real_t tmp_299 = tmp_292 * tmp_39;
      real_t tmp_300 = tmp_296 * tmp_41;
      real_t tmp_301 = 0.37605877282253791 * tmp_46 + 0.58463275527740355 * tmp_47;
      real_t tmp_302 = tmp_19 * ( tmp_301 + tmp_44 );
      real_t tmp_303 = tmp_302 * tmp_43;
      real_t tmp_304 = tmp_302 * tmp_51;
      real_t tmp_305 = tmp_302 * tmp_53;
      real_t tmp_306 = -tmp_293 - tmp_294 - tmp_297 - tmp_298 - tmp_299 - tmp_300 - tmp_303 - tmp_304 - tmp_305 + 1;
      real_t tmp_307 = tmp_75 * ( tmp_291 + tmp_88 );
      real_t tmp_308 = tmp_75 * ( tmp_295 + tmp_90 );
      real_t tmp_309 = tmp_75 * ( tmp_301 + tmp_92 );
      real_t tmp_310 = tmp_98 * ( tmp_57 * ( tmp_307 * tmp_85 + tmp_308 * tmp_82 + tmp_309 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_307 * tmp_86 + tmp_308 * tmp_83 + tmp_309 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_307 * tmp_84 + tmp_308 * tmp_81 + tmp_309 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_311 = 0.020848748529055869 * tmp_100;
      real_t tmp_312 = 0.041227165399737475 * tmp_22 + 0.1711304259088916 * tmp_23;
      real_t tmp_313 = tmp_19 * ( tmp_20 + tmp_312 );
      real_t tmp_314 = tmp_313 * tmp_8;
      real_t tmp_315 = tmp_27 * tmp_313;
      real_t tmp_316 = 0.041227165399737475 * tmp_32 + 0.1711304259088916 * tmp_33;
      real_t tmp_317 = tmp_19 * ( tmp_30 + tmp_316 );
      real_t tmp_318 = tmp_29 * tmp_317;
      real_t tmp_319 = tmp_317 * tmp_37;
      real_t tmp_320 = tmp_313 * tmp_39;
      real_t tmp_321 = tmp_317 * tmp_41;
      real_t tmp_322 = 0.041227165399737475 * tmp_46 + 0.1711304259088916 * tmp_47;
      real_t tmp_323 = tmp_19 * ( tmp_322 + tmp_44 );
      real_t tmp_324 = tmp_323 * tmp_43;
      real_t tmp_325 = tmp_323 * tmp_51;
      real_t tmp_326 = tmp_323 * tmp_53;
      real_t tmp_327 = -tmp_314 - tmp_315 - tmp_318 - tmp_319 - tmp_320 - tmp_321 - tmp_324 - tmp_325 - tmp_326 + 1;
      real_t tmp_328 = tmp_75 * ( tmp_312 + tmp_88 );
      real_t tmp_329 = tmp_75 * ( tmp_316 + tmp_90 );
      real_t tmp_330 = tmp_75 * ( tmp_322 + tmp_92 );
      real_t tmp_331 = tmp_98 * ( tmp_57 * ( tmp_328 * tmp_85 + tmp_329 * tmp_82 + tmp_330 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_328 * tmp_86 + tmp_329 * tmp_83 + tmp_330 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_328 * tmp_84 + tmp_329 * tmp_81 + tmp_330 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_332 = 0.019202922745021479 * tmp_100;
      real_t tmp_333 = 0.40446199974765351 * tmp_22 + 0.19107600050469298 * tmp_23;
      real_t tmp_334 = tmp_19 * ( tmp_20 + tmp_333 );
      real_t tmp_335 = tmp_334 * tmp_8;
      real_t tmp_336 = tmp_27 * tmp_334;
      real_t tmp_337 = 0.40446199974765351 * tmp_32 + 0.19107600050469298 * tmp_33;
      real_t tmp_338 = tmp_19 * ( tmp_30 + tmp_337 );
      real_t tmp_339 = tmp_29 * tmp_338;
      real_t tmp_340 = tmp_338 * tmp_37;
      real_t tmp_341 = tmp_334 * tmp_39;
      real_t tmp_342 = tmp_338 * tmp_41;
      real_t tmp_343 = 0.40446199974765351 * tmp_46 + 0.19107600050469298 * tmp_47;
      real_t tmp_344 = tmp_19 * ( tmp_343 + tmp_44 );
      real_t tmp_345 = tmp_344 * tmp_43;
      real_t tmp_346 = tmp_344 * tmp_51;
      real_t tmp_347 = tmp_344 * tmp_53;
      real_t tmp_348 = -tmp_335 - tmp_336 - tmp_339 - tmp_340 - tmp_341 - tmp_342 - tmp_345 - tmp_346 - tmp_347 + 1;
      real_t tmp_349 = tmp_75 * ( tmp_333 + tmp_88 );
      real_t tmp_350 = tmp_75 * ( tmp_337 + tmp_90 );
      real_t tmp_351 = tmp_75 * ( tmp_343 + tmp_92 );
      real_t tmp_352 = tmp_98 * ( tmp_57 * ( tmp_349 * tmp_85 + tmp_350 * tmp_82 + tmp_351 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_349 * tmp_86 + tmp_350 * tmp_83 + tmp_351 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_349 * tmp_84 + tmp_350 * tmp_81 + tmp_351 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_353 = 0.042507265838595799 * tmp_100;
      real_t tmp_354 = 0.039308471900058539 * tmp_22 + 0.37605877282253791 * tmp_23;
      real_t tmp_355 = tmp_19 * ( tmp_20 + tmp_354 );
      real_t tmp_356 = tmp_355 * tmp_8;
      real_t tmp_357 = tmp_27 * tmp_355;
      real_t tmp_358 = 0.039308471900058539 * tmp_32 + 0.37605877282253791 * tmp_33;
      real_t tmp_359 = tmp_19 * ( tmp_30 + tmp_358 );
      real_t tmp_360 = tmp_29 * tmp_359;
      real_t tmp_361 = tmp_359 * tmp_37;
      real_t tmp_362 = tmp_355 * tmp_39;
      real_t tmp_363 = tmp_359 * tmp_41;
      real_t tmp_364 = 0.039308471900058539 * tmp_46 + 0.37605877282253791 * tmp_47;
      real_t tmp_365 = tmp_19 * ( tmp_364 + tmp_44 );
      real_t tmp_366 = tmp_365 * tmp_43;
      real_t tmp_367 = tmp_365 * tmp_51;
      real_t tmp_368 = tmp_365 * tmp_53;
      real_t tmp_369 = -tmp_356 - tmp_357 - tmp_360 - tmp_361 - tmp_362 - tmp_363 - tmp_366 - tmp_367 - tmp_368 + 1;
      real_t tmp_370 = tmp_75 * ( tmp_354 + tmp_88 );
      real_t tmp_371 = tmp_75 * ( tmp_358 + tmp_90 );
      real_t tmp_372 = tmp_75 * ( tmp_364 + tmp_92 );
      real_t tmp_373 = tmp_98 * ( tmp_57 * ( tmp_370 * tmp_85 + tmp_371 * tmp_82 + tmp_372 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_370 * tmp_86 + tmp_371 * tmp_83 + tmp_372 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_370 * tmp_84 + tmp_371 * tmp_81 + tmp_372 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_374 = 0.020848748529055869 * tmp_100;
      real_t tmp_375 = 0.93718850182767688 * tmp_22 + 0.031405749086161582 * tmp_23;
      real_t tmp_376 = tmp_19 * ( tmp_20 + tmp_375 );
      real_t tmp_377 = tmp_376 * tmp_8;
      real_t tmp_378 = tmp_27 * tmp_376;
      real_t tmp_379 = 0.93718850182767688 * tmp_32 + 0.031405749086161582 * tmp_33;
      real_t tmp_380 = tmp_19 * ( tmp_30 + tmp_379 );
      real_t tmp_381 = tmp_29 * tmp_380;
      real_t tmp_382 = tmp_37 * tmp_380;
      real_t tmp_383 = tmp_376 * tmp_39;
      real_t tmp_384 = tmp_380 * tmp_41;
      real_t tmp_385 = 0.93718850182767688 * tmp_46 + 0.031405749086161582 * tmp_47;
      real_t tmp_386 = tmp_19 * ( tmp_385 + tmp_44 );
      real_t tmp_387 = tmp_386 * tmp_43;
      real_t tmp_388 = tmp_386 * tmp_51;
      real_t tmp_389 = tmp_386 * tmp_53;
      real_t tmp_390 = -tmp_377 - tmp_378 - tmp_381 - tmp_382 - tmp_383 - tmp_384 - tmp_387 - tmp_388 - tmp_389 + 1;
      real_t tmp_391 = tmp_75 * ( tmp_375 + tmp_88 );
      real_t tmp_392 = tmp_75 * ( tmp_379 + tmp_90 );
      real_t tmp_393 = tmp_75 * ( tmp_385 + tmp_92 );
      real_t tmp_394 = tmp_98 * ( tmp_57 * ( tmp_391 * tmp_85 + tmp_392 * tmp_82 + tmp_393 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_391 * tmp_86 + tmp_392 * tmp_83 + tmp_393 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_391 * tmp_84 + tmp_392 * tmp_81 + tmp_393 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_395 = 0.0068572537431980923 * tmp_100;
      real_t tmp_396 = 0.60796128279561268 * tmp_22 + 0.19601935860219369 * tmp_23;
      real_t tmp_397 = tmp_19 * ( tmp_20 + tmp_396 );
      real_t tmp_398 = tmp_397 * tmp_8;
      real_t tmp_399 = tmp_27 * tmp_397;
      real_t tmp_400 = 0.60796128279561268 * tmp_32 + 0.19601935860219369 * tmp_33;
      real_t tmp_401 = tmp_19 * ( tmp_30 + tmp_400 );
      real_t tmp_402 = tmp_29 * tmp_401;
      real_t tmp_403 = tmp_37 * tmp_401;
      real_t tmp_404 = tmp_39 * tmp_397;
      real_t tmp_405 = tmp_401 * tmp_41;
      real_t tmp_406 = 0.60796128279561268 * tmp_46 + 0.19601935860219369 * tmp_47;
      real_t tmp_407 = tmp_19 * ( tmp_406 + tmp_44 );
      real_t tmp_408 = tmp_407 * tmp_43;
      real_t tmp_409 = tmp_407 * tmp_51;
      real_t tmp_410 = tmp_407 * tmp_53;
      real_t tmp_411 = -tmp_398 - tmp_399 - tmp_402 - tmp_403 - tmp_404 - tmp_405 - tmp_408 - tmp_409 - tmp_410 + 1;
      real_t tmp_412 = tmp_75 * ( tmp_396 + tmp_88 );
      real_t tmp_413 = tmp_75 * ( tmp_400 + tmp_90 );
      real_t tmp_414 = tmp_75 * ( tmp_406 + tmp_92 );
      real_t tmp_415 = tmp_98 * ( tmp_57 * ( tmp_412 * tmp_85 + tmp_413 * tmp_82 + tmp_414 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_412 * tmp_86 + tmp_413 * tmp_83 + tmp_414 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_412 * tmp_84 + tmp_413 * tmp_81 + tmp_414 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_416 = 0.037198804536718075 * tmp_100;
      real_t tmp_417 = 0.19107600050469298 * tmp_22 + 0.40446199974765351 * tmp_23;
      real_t tmp_418 = tmp_19 * ( tmp_20 + tmp_417 );
      real_t tmp_419 = tmp_418 * tmp_8;
      real_t tmp_420 = tmp_27 * tmp_418;
      real_t tmp_421 = 0.19107600050469298 * tmp_32 + 0.40446199974765351 * tmp_33;
      real_t tmp_422 = tmp_19 * ( tmp_30 + tmp_421 );
      real_t tmp_423 = tmp_29 * tmp_422;
      real_t tmp_424 = tmp_37 * tmp_422;
      real_t tmp_425 = tmp_39 * tmp_418;
      real_t tmp_426 = tmp_41 * tmp_422;
      real_t tmp_427 = 0.19107600050469298 * tmp_46 + 0.40446199974765351 * tmp_47;
      real_t tmp_428 = tmp_19 * ( tmp_427 + tmp_44 );
      real_t tmp_429 = tmp_428 * tmp_43;
      real_t tmp_430 = tmp_428 * tmp_51;
      real_t tmp_431 = tmp_428 * tmp_53;
      real_t tmp_432 = -tmp_419 - tmp_420 - tmp_423 - tmp_424 - tmp_425 - tmp_426 - tmp_429 - tmp_430 - tmp_431 + 1;
      real_t tmp_433 = tmp_75 * ( tmp_417 + tmp_88 );
      real_t tmp_434 = tmp_75 * ( tmp_421 + tmp_90 );
      real_t tmp_435 = tmp_75 * ( tmp_427 + tmp_92 );
      real_t tmp_436 = tmp_98 * ( tmp_57 * ( tmp_433 * tmp_85 + tmp_434 * tmp_82 + tmp_435 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_433 * tmp_86 + tmp_434 * tmp_83 + tmp_435 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_433 * tmp_84 + tmp_434 * tmp_81 + tmp_435 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_437 = 0.042507265838595799 * tmp_100;
      real_t tmp_438 = 0.031405749086161582 * tmp_22 + 0.031405749086161582 * tmp_23;
      real_t tmp_439 = tmp_19 * ( tmp_20 + tmp_438 );
      real_t tmp_440 = tmp_439 * tmp_8;
      real_t tmp_441 = tmp_27 * tmp_439;
      real_t tmp_442 = 0.031405749086161582 * tmp_32 + 0.031405749086161582 * tmp_33;
      real_t tmp_443 = tmp_19 * ( tmp_30 + tmp_442 );
      real_t tmp_444 = tmp_29 * tmp_443;
      real_t tmp_445 = tmp_37 * tmp_443;
      real_t tmp_446 = tmp_39 * tmp_439;
      real_t tmp_447 = tmp_41 * tmp_443;
      real_t tmp_448 = 0.031405749086161582 * tmp_46 + 0.031405749086161582 * tmp_47;
      real_t tmp_449 = tmp_19 * ( tmp_44 + tmp_448 );
      real_t tmp_450 = tmp_43 * tmp_449;
      real_t tmp_451 = tmp_449 * tmp_51;
      real_t tmp_452 = tmp_449 * tmp_53;
      real_t tmp_453 = -tmp_440 - tmp_441 - tmp_444 - tmp_445 - tmp_446 - tmp_447 - tmp_450 - tmp_451 - tmp_452 + 1;
      real_t tmp_454 = tmp_75 * ( tmp_438 + tmp_88 );
      real_t tmp_455 = tmp_75 * ( tmp_442 + tmp_90 );
      real_t tmp_456 = tmp_75 * ( tmp_448 + tmp_92 );
      real_t tmp_457 = tmp_98 * ( tmp_57 * ( tmp_454 * tmp_85 + tmp_455 * tmp_82 + tmp_456 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_454 * tmp_86 + tmp_455 * tmp_83 + tmp_456 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_454 * tmp_84 + tmp_455 * tmp_81 + tmp_456 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_458 = 0.0068572537431980923 * tmp_100;
      real_t tmp_459 = 0.19601935860219369 * tmp_22 + 0.19601935860219369 * tmp_23;
      real_t tmp_460 = tmp_19 * ( tmp_20 + tmp_459 );
      real_t tmp_461 = tmp_460 * tmp_8;
      real_t tmp_462 = tmp_27 * tmp_460;
      real_t tmp_463 = 0.19601935860219369 * tmp_32 + 0.19601935860219369 * tmp_33;
      real_t tmp_464 = tmp_19 * ( tmp_30 + tmp_463 );
      real_t tmp_465 = tmp_29 * tmp_464;
      real_t tmp_466 = tmp_37 * tmp_464;
      real_t tmp_467 = tmp_39 * tmp_460;
      real_t tmp_468 = tmp_41 * tmp_464;
      real_t tmp_469 = 0.19601935860219369 * tmp_46 + 0.19601935860219369 * tmp_47;
      real_t tmp_470 = tmp_19 * ( tmp_44 + tmp_469 );
      real_t tmp_471 = tmp_43 * tmp_470;
      real_t tmp_472 = tmp_470 * tmp_51;
      real_t tmp_473 = tmp_470 * tmp_53;
      real_t tmp_474 = -tmp_461 - tmp_462 - tmp_465 - tmp_466 - tmp_467 - tmp_468 - tmp_471 - tmp_472 - tmp_473 + 1;
      real_t tmp_475 = tmp_75 * ( tmp_459 + tmp_88 );
      real_t tmp_476 = tmp_75 * ( tmp_463 + tmp_90 );
      real_t tmp_477 = tmp_75 * ( tmp_469 + tmp_92 );
      real_t tmp_478 = tmp_98 * ( tmp_57 * ( tmp_475 * tmp_85 + tmp_476 * tmp_82 + tmp_477 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_475 * tmp_86 + tmp_476 * tmp_83 + tmp_477 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_475 * tmp_84 + tmp_476 * tmp_81 + tmp_477 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_479 = 0.037198804536718075 * tmp_100;
      real_t tmp_480 = 0.40446199974765351 * tmp_22 + 0.40446199974765351 * tmp_23;
      real_t tmp_481 = tmp_19 * ( tmp_20 + tmp_480 );
      real_t tmp_482 = tmp_481 * tmp_8;
      real_t tmp_483 = tmp_27 * tmp_481;
      real_t tmp_484 = 0.40446199974765351 * tmp_32 + 0.40446199974765351 * tmp_33;
      real_t tmp_485 = tmp_19 * ( tmp_30 + tmp_484 );
      real_t tmp_486 = tmp_29 * tmp_485;
      real_t tmp_487 = tmp_37 * tmp_485;
      real_t tmp_488 = tmp_39 * tmp_481;
      real_t tmp_489 = tmp_41 * tmp_485;
      real_t tmp_490 = 0.40446199974765351 * tmp_46 + 0.40446199974765351 * tmp_47;
      real_t tmp_491 = tmp_19 * ( tmp_44 + tmp_490 );
      real_t tmp_492 = tmp_43 * tmp_491;
      real_t tmp_493 = tmp_491 * tmp_51;
      real_t tmp_494 = tmp_491 * tmp_53;
      real_t tmp_495 = -tmp_482 - tmp_483 - tmp_486 - tmp_487 - tmp_488 - tmp_489 - tmp_492 - tmp_493 - tmp_494 + 1;
      real_t tmp_496 = tmp_75 * ( tmp_480 + tmp_88 );
      real_t tmp_497 = tmp_75 * ( tmp_484 + tmp_90 );
      real_t tmp_498 = tmp_75 * ( tmp_490 + tmp_92 );
      real_t tmp_499 = tmp_98 * ( tmp_57 * ( tmp_496 * tmp_85 + tmp_497 * tmp_82 + tmp_498 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_496 * tmp_86 + tmp_497 * tmp_83 + tmp_498 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_496 * tmp_84 + tmp_497 * tmp_81 + tmp_498 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_500 = 0.042507265838595799 * tmp_100;
      real_t tmp_501 = 0.1711304259088916 * tmp_22 + 0.041227165399737475 * tmp_23;
      real_t tmp_502 = tmp_19 * ( tmp_20 + tmp_501 );
      real_t tmp_503 = tmp_502 * tmp_8;
      real_t tmp_504 = tmp_27 * tmp_502;
      real_t tmp_505 = 0.1711304259088916 * tmp_32 + 0.041227165399737475 * tmp_33;
      real_t tmp_506 = tmp_19 * ( tmp_30 + tmp_505 );
      real_t tmp_507 = tmp_29 * tmp_506;
      real_t tmp_508 = tmp_37 * tmp_506;
      real_t tmp_509 = tmp_39 * tmp_502;
      real_t tmp_510 = tmp_41 * tmp_506;
      real_t tmp_511 = 0.1711304259088916 * tmp_46 + 0.041227165399737475 * tmp_47;
      real_t tmp_512 = tmp_19 * ( tmp_44 + tmp_511 );
      real_t tmp_513 = tmp_43 * tmp_512;
      real_t tmp_514 = tmp_51 * tmp_512;
      real_t tmp_515 = tmp_512 * tmp_53;
      real_t tmp_516 = -tmp_503 - tmp_504 - tmp_507 - tmp_508 - tmp_509 - tmp_510 - tmp_513 - tmp_514 - tmp_515 + 1;
      real_t tmp_517 = tmp_75 * ( tmp_501 + tmp_88 );
      real_t tmp_518 = tmp_75 * ( tmp_505 + tmp_90 );
      real_t tmp_519 = tmp_75 * ( tmp_511 + tmp_92 );
      real_t tmp_520 = tmp_98 * ( tmp_57 * ( tmp_517 * tmp_85 + tmp_518 * tmp_82 + tmp_519 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_61 * ( tmp_517 * tmp_86 + tmp_518 * tmp_83 + tmp_519 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_517 * tmp_84 + tmp_518 * tmp_81 + tmp_519 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_521 = 0.019202922745021479 * tmp_100;
      real_t tmp_522 = tmp_40 + tmp_42 + tmp_54;
      real_t tmp_523 = tmp_110 + tmp_111 + tmp_116;
      real_t tmp_524 = tmp_131 + tmp_132 + tmp_137;
      real_t tmp_525 = tmp_152 + tmp_153 + tmp_158;
      real_t tmp_526 = tmp_173 + tmp_174 + tmp_179;
      real_t tmp_527 = tmp_194 + tmp_195 + tmp_200;
      real_t tmp_528 = tmp_215 + tmp_216 + tmp_221;
      real_t tmp_529 = tmp_236 + tmp_237 + tmp_242;
      real_t tmp_530 = tmp_257 + tmp_258 + tmp_263;
      real_t tmp_531 = tmp_278 + tmp_279 + tmp_284;
      real_t tmp_532 = tmp_299 + tmp_300 + tmp_305;
      real_t tmp_533 = tmp_320 + tmp_321 + tmp_326;
      real_t tmp_534 = tmp_341 + tmp_342 + tmp_347;
      real_t tmp_535 = tmp_362 + tmp_363 + tmp_368;
      real_t tmp_536 = tmp_383 + tmp_384 + tmp_389;
      real_t tmp_537 = tmp_404 + tmp_405 + tmp_410;
      real_t tmp_538 = tmp_425 + tmp_426 + tmp_431;
      real_t tmp_539 = tmp_446 + tmp_447 + tmp_452;
      real_t tmp_540 = tmp_467 + tmp_468 + tmp_473;
      real_t tmp_541 = tmp_488 + tmp_489 + tmp_494;
      real_t tmp_542 = tmp_509 + tmp_510 + tmp_515;
      real_t tmp_543 = tmp_28 + tmp_38 + tmp_52;
      real_t tmp_544 = tmp_105 + tmp_109 + tmp_115;
      real_t tmp_545 = tmp_126 + tmp_130 + tmp_136;
      real_t tmp_546 = tmp_147 + tmp_151 + tmp_157;
      real_t tmp_547 = tmp_168 + tmp_172 + tmp_178;
      real_t tmp_548 = tmp_189 + tmp_193 + tmp_199;
      real_t tmp_549 = tmp_210 + tmp_214 + tmp_220;
      real_t tmp_550 = tmp_231 + tmp_235 + tmp_241;
      real_t tmp_551 = tmp_252 + tmp_256 + tmp_262;
      real_t tmp_552 = tmp_273 + tmp_277 + tmp_283;
      real_t tmp_553 = tmp_294 + tmp_298 + tmp_304;
      real_t tmp_554 = tmp_315 + tmp_319 + tmp_325;
      real_t tmp_555 = tmp_336 + tmp_340 + tmp_346;
      real_t tmp_556 = tmp_357 + tmp_361 + tmp_367;
      real_t tmp_557 = tmp_378 + tmp_382 + tmp_388;
      real_t tmp_558 = tmp_399 + tmp_403 + tmp_409;
      real_t tmp_559 = tmp_420 + tmp_424 + tmp_430;
      real_t tmp_560 = tmp_441 + tmp_445 + tmp_451;
      real_t tmp_561 = tmp_462 + tmp_466 + tmp_472;
      real_t tmp_562 = tmp_483 + tmp_487 + tmp_493;
      real_t tmp_563 = tmp_504 + tmp_508 + tmp_514;
      real_t tmp_564 = tmp_26 + tmp_36 + tmp_50;
      real_t tmp_565 = tmp_104 + tmp_108 + tmp_114;
      real_t tmp_566 = tmp_125 + tmp_129 + tmp_135;
      real_t tmp_567 = tmp_146 + tmp_150 + tmp_156;
      real_t tmp_568 = tmp_167 + tmp_171 + tmp_177;
      real_t tmp_569 = tmp_188 + tmp_192 + tmp_198;
      real_t tmp_570 = tmp_209 + tmp_213 + tmp_219;
      real_t tmp_571 = tmp_230 + tmp_234 + tmp_240;
      real_t tmp_572 = tmp_251 + tmp_255 + tmp_261;
      real_t tmp_573 = tmp_272 + tmp_276 + tmp_282;
      real_t tmp_574 = tmp_293 + tmp_297 + tmp_303;
      real_t tmp_575 = tmp_314 + tmp_318 + tmp_324;
      real_t tmp_576 = tmp_335 + tmp_339 + tmp_345;
      real_t tmp_577 = tmp_356 + tmp_360 + tmp_366;
      real_t tmp_578 = tmp_377 + tmp_381 + tmp_387;
      real_t tmp_579 = tmp_398 + tmp_402 + tmp_408;
      real_t tmp_580 = tmp_419 + tmp_423 + tmp_429;
      real_t tmp_581 = tmp_440 + tmp_444 + tmp_450;
      real_t tmp_582 = tmp_461 + tmp_465 + tmp_471;
      real_t tmp_583 = tmp_482 + tmp_486 + tmp_492;
      real_t tmp_584 = tmp_503 + tmp_507 + tmp_513;
      real_t a_0_0   = tmp_101 * ( -tmp_55 * tmp_87 - tmp_55 * tmp_99 ) + tmp_122 * ( -tmp_117 * tmp_121 - tmp_117 * tmp_87 ) +
                     tmp_143 * ( -tmp_138 * tmp_142 - tmp_138 * tmp_87 ) + tmp_164 * ( -tmp_159 * tmp_163 - tmp_159 * tmp_87 ) +
                     tmp_185 * ( -tmp_180 * tmp_184 - tmp_180 * tmp_87 ) + tmp_206 * ( -tmp_201 * tmp_205 - tmp_201 * tmp_87 ) +
                     tmp_227 * ( -tmp_222 * tmp_226 - tmp_222 * tmp_87 ) + tmp_248 * ( -tmp_243 * tmp_247 - tmp_243 * tmp_87 ) +
                     tmp_269 * ( -tmp_264 * tmp_268 - tmp_264 * tmp_87 ) + tmp_290 * ( -tmp_285 * tmp_289 - tmp_285 * tmp_87 ) +
                     tmp_311 * ( -tmp_306 * tmp_310 - tmp_306 * tmp_87 ) + tmp_332 * ( -tmp_327 * tmp_331 - tmp_327 * tmp_87 ) +
                     tmp_353 * ( -tmp_348 * tmp_352 - tmp_348 * tmp_87 ) + tmp_374 * ( -tmp_369 * tmp_373 - tmp_369 * tmp_87 ) +
                     tmp_395 * ( -tmp_390 * tmp_394 - tmp_390 * tmp_87 ) + tmp_416 * ( -tmp_411 * tmp_415 - tmp_411 * tmp_87 ) +
                     tmp_437 * ( -tmp_432 * tmp_436 - tmp_432 * tmp_87 ) + tmp_458 * ( -tmp_453 * tmp_457 - tmp_453 * tmp_87 ) +
                     tmp_479 * ( -tmp_474 * tmp_478 - tmp_474 * tmp_87 ) + tmp_500 * ( -tmp_495 * tmp_499 - tmp_495 * tmp_87 ) +
                     tmp_521 * ( -tmp_516 * tmp_520 - tmp_516 * tmp_87 );
      real_t a_1_0 = tmp_101 * ( -tmp_522 * tmp_87 - tmp_522 * tmp_99 ) + tmp_122 * ( -tmp_121 * tmp_523 - tmp_523 * tmp_87 ) +
                     tmp_143 * ( -tmp_142 * tmp_524 - tmp_524 * tmp_87 ) + tmp_164 * ( -tmp_163 * tmp_525 - tmp_525 * tmp_87 ) +
                     tmp_185 * ( -tmp_184 * tmp_526 - tmp_526 * tmp_87 ) + tmp_206 * ( -tmp_205 * tmp_527 - tmp_527 * tmp_87 ) +
                     tmp_227 * ( -tmp_226 * tmp_528 - tmp_528 * tmp_87 ) + tmp_248 * ( -tmp_247 * tmp_529 - tmp_529 * tmp_87 ) +
                     tmp_269 * ( -tmp_268 * tmp_530 - tmp_530 * tmp_87 ) + tmp_290 * ( -tmp_289 * tmp_531 - tmp_531 * tmp_87 ) +
                     tmp_311 * ( -tmp_310 * tmp_532 - tmp_532 * tmp_87 ) + tmp_332 * ( -tmp_331 * tmp_533 - tmp_533 * tmp_87 ) +
                     tmp_353 * ( -tmp_352 * tmp_534 - tmp_534 * tmp_87 ) + tmp_374 * ( -tmp_373 * tmp_535 - tmp_535 * tmp_87 ) +
                     tmp_395 * ( -tmp_394 * tmp_536 - tmp_536 * tmp_87 ) + tmp_416 * ( -tmp_415 * tmp_537 - tmp_537 * tmp_87 ) +
                     tmp_437 * ( -tmp_436 * tmp_538 - tmp_538 * tmp_87 ) + tmp_458 * ( -tmp_457 * tmp_539 - tmp_539 * tmp_87 ) +
                     tmp_479 * ( -tmp_478 * tmp_540 - tmp_540 * tmp_87 ) + tmp_500 * ( -tmp_499 * tmp_541 - tmp_541 * tmp_87 ) +
                     tmp_521 * ( -tmp_520 * tmp_542 - tmp_542 * tmp_87 );
      real_t a_2_0 = tmp_101 * ( -tmp_543 * tmp_87 - tmp_543 * tmp_99 ) + tmp_122 * ( -tmp_121 * tmp_544 - tmp_544 * tmp_87 ) +
                     tmp_143 * ( -tmp_142 * tmp_545 - tmp_545 * tmp_87 ) + tmp_164 * ( -tmp_163 * tmp_546 - tmp_546 * tmp_87 ) +
                     tmp_185 * ( -tmp_184 * tmp_547 - tmp_547 * tmp_87 ) + tmp_206 * ( -tmp_205 * tmp_548 - tmp_548 * tmp_87 ) +
                     tmp_227 * ( -tmp_226 * tmp_549 - tmp_549 * tmp_87 ) + tmp_248 * ( -tmp_247 * tmp_550 - tmp_550 * tmp_87 ) +
                     tmp_269 * ( -tmp_268 * tmp_551 - tmp_551 * tmp_87 ) + tmp_290 * ( -tmp_289 * tmp_552 - tmp_552 * tmp_87 ) +
                     tmp_311 * ( -tmp_310 * tmp_553 - tmp_553 * tmp_87 ) + tmp_332 * ( -tmp_331 * tmp_554 - tmp_554 * tmp_87 ) +
                     tmp_353 * ( -tmp_352 * tmp_555 - tmp_555 * tmp_87 ) + tmp_374 * ( -tmp_373 * tmp_556 - tmp_556 * tmp_87 ) +
                     tmp_395 * ( -tmp_394 * tmp_557 - tmp_557 * tmp_87 ) + tmp_416 * ( -tmp_415 * tmp_558 - tmp_558 * tmp_87 ) +
                     tmp_437 * ( -tmp_436 * tmp_559 - tmp_559 * tmp_87 ) + tmp_458 * ( -tmp_457 * tmp_560 - tmp_560 * tmp_87 ) +
                     tmp_479 * ( -tmp_478 * tmp_561 - tmp_561 * tmp_87 ) + tmp_500 * ( -tmp_499 * tmp_562 - tmp_562 * tmp_87 ) +
                     tmp_521 * ( -tmp_520 * tmp_563 - tmp_563 * tmp_87 );
      real_t a_3_0 = tmp_101 * ( -tmp_564 * tmp_87 - tmp_564 * tmp_99 ) + tmp_122 * ( -tmp_121 * tmp_565 - tmp_565 * tmp_87 ) +
                     tmp_143 * ( -tmp_142 * tmp_566 - tmp_566 * tmp_87 ) + tmp_164 * ( -tmp_163 * tmp_567 - tmp_567 * tmp_87 ) +
                     tmp_185 * ( -tmp_184 * tmp_568 - tmp_568 * tmp_87 ) + tmp_206 * ( -tmp_205 * tmp_569 - tmp_569 * tmp_87 ) +
                     tmp_227 * ( -tmp_226 * tmp_570 - tmp_570 * tmp_87 ) + tmp_248 * ( -tmp_247 * tmp_571 - tmp_571 * tmp_87 ) +
                     tmp_269 * ( -tmp_268 * tmp_572 - tmp_572 * tmp_87 ) + tmp_290 * ( -tmp_289 * tmp_573 - tmp_573 * tmp_87 ) +
                     tmp_311 * ( -tmp_310 * tmp_574 - tmp_574 * tmp_87 ) + tmp_332 * ( -tmp_331 * tmp_575 - tmp_575 * tmp_87 ) +
                     tmp_353 * ( -tmp_352 * tmp_576 - tmp_576 * tmp_87 ) + tmp_374 * ( -tmp_373 * tmp_577 - tmp_577 * tmp_87 ) +
                     tmp_395 * ( -tmp_394 * tmp_578 - tmp_578 * tmp_87 ) + tmp_416 * ( -tmp_415 * tmp_579 - tmp_579 * tmp_87 ) +
                     tmp_437 * ( -tmp_436 * tmp_580 - tmp_580 * tmp_87 ) + tmp_458 * ( -tmp_457 * tmp_581 - tmp_581 * tmp_87 ) +
                     tmp_479 * ( -tmp_478 * tmp_582 - tmp_582 * tmp_87 ) + tmp_500 * ( -tmp_499 * tmp_583 - tmp_583 * tmp_87 ) +
                     tmp_521 * ( -tmp_520 * tmp_584 - tmp_584 * tmp_87 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
      elMat( 3, 0 ) = a_3_0;
   }

   void integrateFacetDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                           const std::vector< Point3D >& coordsFacet,
                                           const Point3D&,
                                           const Point3D&     outwardNormal,
                                           const DGBasisInfo& trialBasis,
                                           const DGBasisInfo& testBasis,
                                           int                trialDegree,
                                           int                testDegree,
                                           MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = p_affine_2_0 + tmp_0;
      real_t tmp_5  = p_affine_1_1 + tmp_2;
      real_t tmp_6  = tmp_1 * tmp_3 - tmp_4 * tmp_5;
      real_t tmp_7  = -p_affine_0_2;
      real_t tmp_8  = p_affine_3_2 + tmp_7;
      real_t tmp_9  = tmp_3 * tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_7;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11 * tmp_4;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_13 * tmp_14;
      real_t tmp_16 = tmp_11 * tmp_14;
      real_t tmp_17 = tmp_4 * tmp_8;
      real_t tmp_18 = tmp_13 * tmp_3;
      real_t tmp_19 =
          1.0 / ( -tmp_1 * tmp_16 + tmp_1 * tmp_9 + tmp_10 * tmp_12 - tmp_10 * tmp_18 + tmp_15 * tmp_5 - tmp_17 * tmp_5 );
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19 * ( 0.031405749086161582 * tmp_21 + 0.93718850182767688 * tmp_22 + tmp_23 );
      real_t tmp_25 = tmp_24 * tmp_6;
      real_t tmp_26 = -tmp_1 * tmp_11 + tmp_13 * tmp_5;
      real_t tmp_27 = tmp_24 * tmp_26;
      real_t tmp_28 = -tmp_1 * tmp_14 + tmp_10 * tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19 * ( 0.031405749086161582 * tmp_30 + 0.93718850182767688 * tmp_31 + tmp_32 );
      real_t tmp_34 = tmp_28 * tmp_33;
      real_t tmp_35 = tmp_1 * tmp_8 - tmp_10 * tmp_13;
      real_t tmp_36 = tmp_33 * tmp_35;
      real_t tmp_37 = tmp_12 - tmp_18;
      real_t tmp_38 = tmp_24 * tmp_37;
      real_t tmp_39 = tmp_15 - tmp_17;
      real_t tmp_40 = tmp_33 * tmp_39;
      real_t tmp_41 = -tmp_10 * tmp_3 + tmp_14 * tmp_5;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19 * ( 0.031405749086161582 * tmp_43 + 0.93718850182767688 * tmp_44 + tmp_45 );
      real_t tmp_47 = tmp_41 * tmp_46;
      real_t tmp_48 = tmp_10 * tmp_11 - tmp_5 * tmp_8;
      real_t tmp_49 = tmp_46 * tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46 * tmp_50;
      real_t tmp_52 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_53 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_54 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_55 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_58 =
          1.0 * std::pow( ( std::abs( tmp_22 * tmp_55 - tmp_31 * tmp_57 ) * std::abs( tmp_22 * tmp_55 - tmp_31 * tmp_57 ) ) +
                              ( std::abs( tmp_22 * tmp_56 - tmp_44 * tmp_57 ) * std::abs( tmp_22 * tmp_56 - tmp_44 * tmp_57 ) ) +
                              ( std::abs( tmp_31 * tmp_56 - tmp_44 * tmp_55 ) * std::abs( tmp_31 * tmp_56 - tmp_44 * tmp_55 ) ),
                          0.25 );
      real_t tmp_59 = 0.0068572537431980923 * tmp_58 *
                      ( tmp_11 * ( tmp_54 - 1.0 / 4.0 ) + tmp_3 * ( tmp_53 - 1.0 / 4.0 ) + tmp_5 * ( tmp_52 - 1.0 / 4.0 ) );
      real_t tmp_60 = tmp_19 * ( 0.19601935860219369 * tmp_21 + 0.60796128279561268 * tmp_22 + tmp_23 );
      real_t tmp_61 = tmp_6 * tmp_60;
      real_t tmp_62 = tmp_26 * tmp_60;
      real_t tmp_63 = tmp_19 * ( 0.19601935860219369 * tmp_30 + 0.60796128279561268 * tmp_31 + tmp_32 );
      real_t tmp_64 = tmp_28 * tmp_63;
      real_t tmp_65 = tmp_35 * tmp_63;
      real_t tmp_66 = tmp_37 * tmp_60;
      real_t tmp_67 = tmp_39 * tmp_63;
      real_t tmp_68 = tmp_19 * ( 0.19601935860219369 * tmp_43 + 0.60796128279561268 * tmp_44 + tmp_45 );
      real_t tmp_69 = tmp_41 * tmp_68;
      real_t tmp_70 = tmp_48 * tmp_68;
      real_t tmp_71 = tmp_50 * tmp_68;
      real_t tmp_72 = tmp_66 + tmp_67 + tmp_71;
      real_t tmp_73 = tmp_62 + tmp_65 + tmp_70;
      real_t tmp_74 = tmp_61 + tmp_64 + tmp_69;
      real_t tmp_75 = 0.037198804536718075 * tmp_58 *
                      ( tmp_11 * ( tmp_74 - 1.0 / 4.0 ) + tmp_3 * ( tmp_73 - 1.0 / 4.0 ) + tmp_5 * ( tmp_72 - 1.0 / 4.0 ) );
      real_t tmp_76 = tmp_19 * ( 0.37605877282253791 * tmp_21 + 0.039308471900058539 * tmp_22 + tmp_23 );
      real_t tmp_77 = tmp_6 * tmp_76;
      real_t tmp_78 = tmp_26 * tmp_76;
      real_t tmp_79 = tmp_19 * ( 0.37605877282253791 * tmp_30 + 0.039308471900058539 * tmp_31 + tmp_32 );
      real_t tmp_80 = tmp_28 * tmp_79;
      real_t tmp_81 = tmp_35 * tmp_79;
      real_t tmp_82 = tmp_37 * tmp_76;
      real_t tmp_83 = tmp_39 * tmp_79;
      real_t tmp_84 = tmp_19 * ( 0.37605877282253791 * tmp_43 + 0.039308471900058539 * tmp_44 + tmp_45 );
      real_t tmp_85 = tmp_41 * tmp_84;
      real_t tmp_86 = tmp_48 * tmp_84;
      real_t tmp_87 = tmp_50 * tmp_84;
      real_t tmp_88 = tmp_82 + tmp_83 + tmp_87;
      real_t tmp_89 = tmp_78 + tmp_81 + tmp_86;
      real_t tmp_90 = tmp_77 + tmp_80 + tmp_85;
      real_t tmp_91 = 0.020848748529055869 * tmp_58 *
                      ( tmp_11 * ( tmp_90 - 1.0 / 4.0 ) + tmp_3 * ( tmp_89 - 1.0 / 4.0 ) + tmp_5 * ( tmp_88 - 1.0 / 4.0 ) );
      real_t tmp_92  = tmp_19 * ( 0.78764240869137092 * tmp_21 + 0.1711304259088916 * tmp_22 + tmp_23 );
      real_t tmp_93  = tmp_6 * tmp_92;
      real_t tmp_94  = tmp_26 * tmp_92;
      real_t tmp_95  = tmp_19 * ( 0.78764240869137092 * tmp_30 + 0.1711304259088916 * tmp_31 + tmp_32 );
      real_t tmp_96  = tmp_28 * tmp_95;
      real_t tmp_97  = tmp_35 * tmp_95;
      real_t tmp_98  = tmp_37 * tmp_92;
      real_t tmp_99  = tmp_39 * tmp_95;
      real_t tmp_100 = tmp_19 * ( 0.78764240869137092 * tmp_43 + 0.1711304259088916 * tmp_44 + tmp_45 );
      real_t tmp_101 = tmp_100 * tmp_41;
      real_t tmp_102 = tmp_100 * tmp_48;
      real_t tmp_103 = tmp_100 * tmp_50;
      real_t tmp_104 = tmp_103 + tmp_98 + tmp_99;
      real_t tmp_105 = tmp_102 + tmp_94 + tmp_97;
      real_t tmp_106 = tmp_101 + tmp_93 + tmp_96;
      real_t tmp_107 = 0.019202922745021479 * tmp_58 *
                       ( tmp_11 * ( tmp_106 - 1.0 / 4.0 ) + tmp_3 * ( tmp_105 - 1.0 / 4.0 ) + tmp_5 * ( tmp_104 - 1.0 / 4.0 ) );
      real_t tmp_108 = tmp_19 * ( 0.58463275527740355 * tmp_21 + 0.37605877282253791 * tmp_22 + tmp_23 );
      real_t tmp_109 = tmp_108 * tmp_6;
      real_t tmp_110 = tmp_108 * tmp_26;
      real_t tmp_111 = tmp_19 * ( 0.58463275527740355 * tmp_30 + 0.37605877282253791 * tmp_31 + tmp_32 );
      real_t tmp_112 = tmp_111 * tmp_28;
      real_t tmp_113 = tmp_111 * tmp_35;
      real_t tmp_114 = tmp_108 * tmp_37;
      real_t tmp_115 = tmp_111 * tmp_39;
      real_t tmp_116 = tmp_19 * ( 0.58463275527740355 * tmp_43 + 0.37605877282253791 * tmp_44 + tmp_45 );
      real_t tmp_117 = tmp_116 * tmp_41;
      real_t tmp_118 = tmp_116 * tmp_48;
      real_t tmp_119 = tmp_116 * tmp_50;
      real_t tmp_120 = tmp_114 + tmp_115 + tmp_119;
      real_t tmp_121 = tmp_110 + tmp_113 + tmp_118;
      real_t tmp_122 = tmp_109 + tmp_112 + tmp_117;
      real_t tmp_123 = 0.020848748529055869 * tmp_58 *
                       ( tmp_11 * ( tmp_122 - 1.0 / 4.0 ) + tmp_3 * ( tmp_121 - 1.0 / 4.0 ) + tmp_5 * ( tmp_120 - 1.0 / 4.0 ) );
      real_t tmp_124 = tmp_19 * ( 0.041227165399737475 * tmp_21 + 0.78764240869137092 * tmp_22 + tmp_23 );
      real_t tmp_125 = tmp_124 * tmp_6;
      real_t tmp_126 = tmp_124 * tmp_26;
      real_t tmp_127 = tmp_19 * ( 0.041227165399737475 * tmp_30 + 0.78764240869137092 * tmp_31 + tmp_32 );
      real_t tmp_128 = tmp_127 * tmp_28;
      real_t tmp_129 = tmp_127 * tmp_35;
      real_t tmp_130 = tmp_124 * tmp_37;
      real_t tmp_131 = tmp_127 * tmp_39;
      real_t tmp_132 = tmp_19 * ( 0.041227165399737475 * tmp_43 + 0.78764240869137092 * tmp_44 + tmp_45 );
      real_t tmp_133 = tmp_132 * tmp_41;
      real_t tmp_134 = tmp_132 * tmp_48;
      real_t tmp_135 = tmp_132 * tmp_50;
      real_t tmp_136 = tmp_130 + tmp_131 + tmp_135;
      real_t tmp_137 = tmp_126 + tmp_129 + tmp_134;
      real_t tmp_138 = tmp_125 + tmp_128 + tmp_133;
      real_t tmp_139 = 0.019202922745021479 * tmp_58 *
                       ( tmp_11 * ( tmp_138 - 1.0 / 4.0 ) + tmp_3 * ( tmp_137 - 1.0 / 4.0 ) + tmp_5 * ( tmp_136 - 1.0 / 4.0 ) );
      real_t tmp_140 = tmp_19 * ( 0.039308471900058539 * tmp_21 + 0.58463275527740355 * tmp_22 + tmp_23 );
      real_t tmp_141 = tmp_140 * tmp_6;
      real_t tmp_142 = tmp_140 * tmp_26;
      real_t tmp_143 = tmp_19 * ( 0.039308471900058539 * tmp_30 + 0.58463275527740355 * tmp_31 + tmp_32 );
      real_t tmp_144 = tmp_143 * tmp_28;
      real_t tmp_145 = tmp_143 * tmp_35;
      real_t tmp_146 = tmp_140 * tmp_37;
      real_t tmp_147 = tmp_143 * tmp_39;
      real_t tmp_148 = tmp_19 * ( 0.039308471900058539 * tmp_43 + 0.58463275527740355 * tmp_44 + tmp_45 );
      real_t tmp_149 = tmp_148 * tmp_41;
      real_t tmp_150 = tmp_148 * tmp_48;
      real_t tmp_151 = tmp_148 * tmp_50;
      real_t tmp_152 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_153 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_154 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_155 = 0.020848748529055869 * tmp_58 *
                       ( tmp_11 * ( tmp_154 - 1.0 / 4.0 ) + tmp_3 * ( tmp_153 - 1.0 / 4.0 ) + tmp_5 * ( tmp_152 - 1.0 / 4.0 ) );
      real_t tmp_156 = tmp_19 * ( 0.78764240869137092 * tmp_21 + 0.041227165399737475 * tmp_22 + tmp_23 );
      real_t tmp_157 = tmp_156 * tmp_6;
      real_t tmp_158 = tmp_156 * tmp_26;
      real_t tmp_159 = tmp_19 * ( 0.78764240869137092 * tmp_30 + 0.041227165399737475 * tmp_31 + tmp_32 );
      real_t tmp_160 = tmp_159 * tmp_28;
      real_t tmp_161 = tmp_159 * tmp_35;
      real_t tmp_162 = tmp_156 * tmp_37;
      real_t tmp_163 = tmp_159 * tmp_39;
      real_t tmp_164 = tmp_19 * ( 0.78764240869137092 * tmp_43 + 0.041227165399737475 * tmp_44 + tmp_45 );
      real_t tmp_165 = tmp_164 * tmp_41;
      real_t tmp_166 = tmp_164 * tmp_48;
      real_t tmp_167 = tmp_164 * tmp_50;
      real_t tmp_168 = tmp_162 + tmp_163 + tmp_167;
      real_t tmp_169 = tmp_158 + tmp_161 + tmp_166;
      real_t tmp_170 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_171 = 0.019202922745021479 * tmp_58 *
                       ( tmp_11 * ( tmp_170 - 1.0 / 4.0 ) + tmp_3 * ( tmp_169 - 1.0 / 4.0 ) + tmp_5 * ( tmp_168 - 1.0 / 4.0 ) );
      real_t tmp_172 = tmp_19 * ( 0.58463275527740355 * tmp_21 + 0.039308471900058539 * tmp_22 + tmp_23 );
      real_t tmp_173 = tmp_172 * tmp_6;
      real_t tmp_174 = tmp_172 * tmp_26;
      real_t tmp_175 = tmp_19 * ( 0.58463275527740355 * tmp_30 + 0.039308471900058539 * tmp_31 + tmp_32 );
      real_t tmp_176 = tmp_175 * tmp_28;
      real_t tmp_177 = tmp_175 * tmp_35;
      real_t tmp_178 = tmp_172 * tmp_37;
      real_t tmp_179 = tmp_175 * tmp_39;
      real_t tmp_180 = tmp_19 * ( 0.58463275527740355 * tmp_43 + 0.039308471900058539 * tmp_44 + tmp_45 );
      real_t tmp_181 = tmp_180 * tmp_41;
      real_t tmp_182 = tmp_180 * tmp_48;
      real_t tmp_183 = tmp_180 * tmp_50;
      real_t tmp_184 = tmp_178 + tmp_179 + tmp_183;
      real_t tmp_185 = tmp_174 + tmp_177 + tmp_182;
      real_t tmp_186 = tmp_173 + tmp_176 + tmp_181;
      real_t tmp_187 = 0.020848748529055869 * tmp_58 *
                       ( tmp_11 * ( tmp_186 - 1.0 / 4.0 ) + tmp_3 * ( tmp_185 - 1.0 / 4.0 ) + tmp_5 * ( tmp_184 - 1.0 / 4.0 ) );
      real_t tmp_188 = tmp_19 * ( 0.1711304259088916 * tmp_21 + 0.78764240869137092 * tmp_22 + tmp_23 );
      real_t tmp_189 = tmp_188 * tmp_6;
      real_t tmp_190 = tmp_188 * tmp_26;
      real_t tmp_191 = tmp_19 * ( 0.1711304259088916 * tmp_30 + 0.78764240869137092 * tmp_31 + tmp_32 );
      real_t tmp_192 = tmp_191 * tmp_28;
      real_t tmp_193 = tmp_191 * tmp_35;
      real_t tmp_194 = tmp_188 * tmp_37;
      real_t tmp_195 = tmp_191 * tmp_39;
      real_t tmp_196 = tmp_19 * ( 0.1711304259088916 * tmp_43 + 0.78764240869137092 * tmp_44 + tmp_45 );
      real_t tmp_197 = tmp_196 * tmp_41;
      real_t tmp_198 = tmp_196 * tmp_48;
      real_t tmp_199 = tmp_196 * tmp_50;
      real_t tmp_200 = tmp_194 + tmp_195 + tmp_199;
      real_t tmp_201 = tmp_190 + tmp_193 + tmp_198;
      real_t tmp_202 = tmp_189 + tmp_192 + tmp_197;
      real_t tmp_203 = 0.019202922745021479 * tmp_58 *
                       ( tmp_11 * ( tmp_202 - 1.0 / 4.0 ) + tmp_3 * ( tmp_201 - 1.0 / 4.0 ) + tmp_5 * ( tmp_200 - 1.0 / 4.0 ) );
      real_t tmp_204 = tmp_19 * ( 0.37605877282253791 * tmp_21 + 0.58463275527740355 * tmp_22 + tmp_23 );
      real_t tmp_205 = tmp_204 * tmp_6;
      real_t tmp_206 = tmp_204 * tmp_26;
      real_t tmp_207 = tmp_19 * ( 0.37605877282253791 * tmp_30 + 0.58463275527740355 * tmp_31 + tmp_32 );
      real_t tmp_208 = tmp_207 * tmp_28;
      real_t tmp_209 = tmp_207 * tmp_35;
      real_t tmp_210 = tmp_204 * tmp_37;
      real_t tmp_211 = tmp_207 * tmp_39;
      real_t tmp_212 = tmp_19 * ( 0.37605877282253791 * tmp_43 + 0.58463275527740355 * tmp_44 + tmp_45 );
      real_t tmp_213 = tmp_212 * tmp_41;
      real_t tmp_214 = tmp_212 * tmp_48;
      real_t tmp_215 = tmp_212 * tmp_50;
      real_t tmp_216 = tmp_210 + tmp_211 + tmp_215;
      real_t tmp_217 = tmp_206 + tmp_209 + tmp_214;
      real_t tmp_218 = tmp_205 + tmp_208 + tmp_213;
      real_t tmp_219 = 0.020848748529055869 * tmp_58 *
                       ( tmp_11 * ( tmp_218 - 1.0 / 4.0 ) + tmp_3 * ( tmp_217 - 1.0 / 4.0 ) + tmp_5 * ( tmp_216 - 1.0 / 4.0 ) );
      real_t tmp_220 = tmp_19 * ( 0.041227165399737475 * tmp_21 + 0.1711304259088916 * tmp_22 + tmp_23 );
      real_t tmp_221 = tmp_220 * tmp_6;
      real_t tmp_222 = tmp_220 * tmp_26;
      real_t tmp_223 = tmp_19 * ( 0.041227165399737475 * tmp_30 + 0.1711304259088916 * tmp_31 + tmp_32 );
      real_t tmp_224 = tmp_223 * tmp_28;
      real_t tmp_225 = tmp_223 * tmp_35;
      real_t tmp_226 = tmp_220 * tmp_37;
      real_t tmp_227 = tmp_223 * tmp_39;
      real_t tmp_228 = tmp_19 * ( 0.041227165399737475 * tmp_43 + 0.1711304259088916 * tmp_44 + tmp_45 );
      real_t tmp_229 = tmp_228 * tmp_41;
      real_t tmp_230 = tmp_228 * tmp_48;
      real_t tmp_231 = tmp_228 * tmp_50;
      real_t tmp_232 = tmp_226 + tmp_227 + tmp_231;
      real_t tmp_233 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_234 = tmp_221 + tmp_224 + tmp_229;
      real_t tmp_235 = 0.019202922745021479 * tmp_58 *
                       ( tmp_11 * ( tmp_234 - 1.0 / 4.0 ) + tmp_3 * ( tmp_233 - 1.0 / 4.0 ) + tmp_5 * ( tmp_232 - 1.0 / 4.0 ) );
      real_t tmp_236 = tmp_19 * ( 0.40446199974765351 * tmp_21 + 0.19107600050469298 * tmp_22 + tmp_23 );
      real_t tmp_237 = tmp_236 * tmp_6;
      real_t tmp_238 = tmp_236 * tmp_26;
      real_t tmp_239 = tmp_19 * ( 0.40446199974765351 * tmp_30 + 0.19107600050469298 * tmp_31 + tmp_32 );
      real_t tmp_240 = tmp_239 * tmp_28;
      real_t tmp_241 = tmp_239 * tmp_35;
      real_t tmp_242 = tmp_236 * tmp_37;
      real_t tmp_243 = tmp_239 * tmp_39;
      real_t tmp_244 = tmp_19 * ( 0.40446199974765351 * tmp_43 + 0.19107600050469298 * tmp_44 + tmp_45 );
      real_t tmp_245 = tmp_244 * tmp_41;
      real_t tmp_246 = tmp_244 * tmp_48;
      real_t tmp_247 = tmp_244 * tmp_50;
      real_t tmp_248 = tmp_242 + tmp_243 + tmp_247;
      real_t tmp_249 = tmp_238 + tmp_241 + tmp_246;
      real_t tmp_250 = tmp_237 + tmp_240 + tmp_245;
      real_t tmp_251 = 0.042507265838595799 * tmp_58 *
                       ( tmp_11 * ( tmp_250 - 1.0 / 4.0 ) + tmp_3 * ( tmp_249 - 1.0 / 4.0 ) + tmp_5 * ( tmp_248 - 1.0 / 4.0 ) );
      real_t tmp_252 = tmp_19 * ( 0.039308471900058539 * tmp_21 + 0.37605877282253791 * tmp_22 + tmp_23 );
      real_t tmp_253 = tmp_252 * tmp_6;
      real_t tmp_254 = tmp_252 * tmp_26;
      real_t tmp_255 = tmp_19 * ( 0.039308471900058539 * tmp_30 + 0.37605877282253791 * tmp_31 + tmp_32 );
      real_t tmp_256 = tmp_255 * tmp_28;
      real_t tmp_257 = tmp_255 * tmp_35;
      real_t tmp_258 = tmp_252 * tmp_37;
      real_t tmp_259 = tmp_255 * tmp_39;
      real_t tmp_260 = tmp_19 * ( 0.039308471900058539 * tmp_43 + 0.37605877282253791 * tmp_44 + tmp_45 );
      real_t tmp_261 = tmp_260 * tmp_41;
      real_t tmp_262 = tmp_260 * tmp_48;
      real_t tmp_263 = tmp_260 * tmp_50;
      real_t tmp_264 = tmp_258 + tmp_259 + tmp_263;
      real_t tmp_265 = tmp_254 + tmp_257 + tmp_262;
      real_t tmp_266 = tmp_253 + tmp_256 + tmp_261;
      real_t tmp_267 = 0.020848748529055869 * tmp_58 *
                       ( tmp_11 * ( tmp_266 - 1.0 / 4.0 ) + tmp_3 * ( tmp_265 - 1.0 / 4.0 ) + tmp_5 * ( tmp_264 - 1.0 / 4.0 ) );
      real_t tmp_268 = tmp_19 * ( 0.93718850182767688 * tmp_21 + 0.031405749086161582 * tmp_22 + tmp_23 );
      real_t tmp_269 = tmp_268 * tmp_6;
      real_t tmp_270 = tmp_26 * tmp_268;
      real_t tmp_271 = tmp_19 * ( 0.93718850182767688 * tmp_30 + 0.031405749086161582 * tmp_31 + tmp_32 );
      real_t tmp_272 = tmp_271 * tmp_28;
      real_t tmp_273 = tmp_271 * tmp_35;
      real_t tmp_274 = tmp_268 * tmp_37;
      real_t tmp_275 = tmp_271 * tmp_39;
      real_t tmp_276 = tmp_19 * ( 0.93718850182767688 * tmp_43 + 0.031405749086161582 * tmp_44 + tmp_45 );
      real_t tmp_277 = tmp_276 * tmp_41;
      real_t tmp_278 = tmp_276 * tmp_48;
      real_t tmp_279 = tmp_276 * tmp_50;
      real_t tmp_280 = tmp_274 + tmp_275 + tmp_279;
      real_t tmp_281 = tmp_270 + tmp_273 + tmp_278;
      real_t tmp_282 = tmp_269 + tmp_272 + tmp_277;
      real_t tmp_283 = 0.0068572537431980923 * tmp_58 *
                       ( tmp_11 * ( tmp_282 - 1.0 / 4.0 ) + tmp_3 * ( tmp_281 - 1.0 / 4.0 ) + tmp_5 * ( tmp_280 - 1.0 / 4.0 ) );
      real_t tmp_284 = tmp_19 * ( 0.60796128279561268 * tmp_21 + 0.19601935860219369 * tmp_22 + tmp_23 );
      real_t tmp_285 = tmp_284 * tmp_6;
      real_t tmp_286 = tmp_26 * tmp_284;
      real_t tmp_287 = tmp_19 * ( 0.60796128279561268 * tmp_30 + 0.19601935860219369 * tmp_31 + tmp_32 );
      real_t tmp_288 = tmp_28 * tmp_287;
      real_t tmp_289 = tmp_287 * tmp_35;
      real_t tmp_290 = tmp_284 * tmp_37;
      real_t tmp_291 = tmp_287 * tmp_39;
      real_t tmp_292 = tmp_19 * ( 0.60796128279561268 * tmp_43 + 0.19601935860219369 * tmp_44 + tmp_45 );
      real_t tmp_293 = tmp_292 * tmp_41;
      real_t tmp_294 = tmp_292 * tmp_48;
      real_t tmp_295 = tmp_292 * tmp_50;
      real_t tmp_296 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_297 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_298 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_299 = 0.037198804536718075 * tmp_58 *
                       ( tmp_11 * ( tmp_298 - 1.0 / 4.0 ) + tmp_3 * ( tmp_297 - 1.0 / 4.0 ) + tmp_5 * ( tmp_296 - 1.0 / 4.0 ) );
      real_t tmp_300 = tmp_19 * ( 0.19107600050469298 * tmp_21 + 0.40446199974765351 * tmp_22 + tmp_23 );
      real_t tmp_301 = tmp_300 * tmp_6;
      real_t tmp_302 = tmp_26 * tmp_300;
      real_t tmp_303 = tmp_19 * ( 0.19107600050469298 * tmp_30 + 0.40446199974765351 * tmp_31 + tmp_32 );
      real_t tmp_304 = tmp_28 * tmp_303;
      real_t tmp_305 = tmp_303 * tmp_35;
      real_t tmp_306 = tmp_300 * tmp_37;
      real_t tmp_307 = tmp_303 * tmp_39;
      real_t tmp_308 = tmp_19 * ( 0.19107600050469298 * tmp_43 + 0.40446199974765351 * tmp_44 + tmp_45 );
      real_t tmp_309 = tmp_308 * tmp_41;
      real_t tmp_310 = tmp_308 * tmp_48;
      real_t tmp_311 = tmp_308 * tmp_50;
      real_t tmp_312 = tmp_306 + tmp_307 + tmp_311;
      real_t tmp_313 = tmp_302 + tmp_305 + tmp_310;
      real_t tmp_314 = tmp_301 + tmp_304 + tmp_309;
      real_t tmp_315 = 0.042507265838595799 * tmp_58 *
                       ( tmp_11 * ( tmp_314 - 1.0 / 4.0 ) + tmp_3 * ( tmp_313 - 1.0 / 4.0 ) + tmp_5 * ( tmp_312 - 1.0 / 4.0 ) );
      real_t tmp_316 = tmp_19 * ( 0.031405749086161582 * tmp_21 + 0.031405749086161582 * tmp_22 + tmp_23 );
      real_t tmp_317 = tmp_316 * tmp_6;
      real_t tmp_318 = tmp_26 * tmp_316;
      real_t tmp_319 = tmp_19 * ( 0.031405749086161582 * tmp_30 + 0.031405749086161582 * tmp_31 + tmp_32 );
      real_t tmp_320 = tmp_28 * tmp_319;
      real_t tmp_321 = tmp_319 * tmp_35;
      real_t tmp_322 = tmp_316 * tmp_37;
      real_t tmp_323 = tmp_319 * tmp_39;
      real_t tmp_324 = tmp_19 * ( 0.031405749086161582 * tmp_43 + 0.031405749086161582 * tmp_44 + tmp_45 );
      real_t tmp_325 = tmp_324 * tmp_41;
      real_t tmp_326 = tmp_324 * tmp_48;
      real_t tmp_327 = tmp_324 * tmp_50;
      real_t tmp_328 = tmp_322 + tmp_323 + tmp_327;
      real_t tmp_329 = tmp_318 + tmp_321 + tmp_326;
      real_t tmp_330 = tmp_317 + tmp_320 + tmp_325;
      real_t tmp_331 = 0.0068572537431980923 * tmp_58 *
                       ( tmp_11 * ( tmp_330 - 1.0 / 4.0 ) + tmp_3 * ( tmp_329 - 1.0 / 4.0 ) + tmp_5 * ( tmp_328 - 1.0 / 4.0 ) );
      real_t tmp_332 = tmp_19 * ( 0.19601935860219369 * tmp_21 + 0.19601935860219369 * tmp_22 + tmp_23 );
      real_t tmp_333 = tmp_332 * tmp_6;
      real_t tmp_334 = tmp_26 * tmp_332;
      real_t tmp_335 = tmp_19 * ( 0.19601935860219369 * tmp_30 + 0.19601935860219369 * tmp_31 + tmp_32 );
      real_t tmp_336 = tmp_28 * tmp_335;
      real_t tmp_337 = tmp_335 * tmp_35;
      real_t tmp_338 = tmp_332 * tmp_37;
      real_t tmp_339 = tmp_335 * tmp_39;
      real_t tmp_340 = tmp_19 * ( 0.19601935860219369 * tmp_43 + 0.19601935860219369 * tmp_44 + tmp_45 );
      real_t tmp_341 = tmp_340 * tmp_41;
      real_t tmp_342 = tmp_340 * tmp_48;
      real_t tmp_343 = tmp_340 * tmp_50;
      real_t tmp_344 = tmp_338 + tmp_339 + tmp_343;
      real_t tmp_345 = tmp_334 + tmp_337 + tmp_342;
      real_t tmp_346 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_347 = 0.037198804536718075 * tmp_58 *
                       ( tmp_11 * ( tmp_346 - 1.0 / 4.0 ) + tmp_3 * ( tmp_345 - 1.0 / 4.0 ) + tmp_5 * ( tmp_344 - 1.0 / 4.0 ) );
      real_t tmp_348 = tmp_19 * ( 0.40446199974765351 * tmp_21 + 0.40446199974765351 * tmp_22 + tmp_23 );
      real_t tmp_349 = tmp_348 * tmp_6;
      real_t tmp_350 = tmp_26 * tmp_348;
      real_t tmp_351 = tmp_19 * ( 0.40446199974765351 * tmp_30 + 0.40446199974765351 * tmp_31 + tmp_32 );
      real_t tmp_352 = tmp_28 * tmp_351;
      real_t tmp_353 = tmp_35 * tmp_351;
      real_t tmp_354 = tmp_348 * tmp_37;
      real_t tmp_355 = tmp_351 * tmp_39;
      real_t tmp_356 = tmp_19 * ( 0.40446199974765351 * tmp_43 + 0.40446199974765351 * tmp_44 + tmp_45 );
      real_t tmp_357 = tmp_356 * tmp_41;
      real_t tmp_358 = tmp_356 * tmp_48;
      real_t tmp_359 = tmp_356 * tmp_50;
      real_t tmp_360 = tmp_354 + tmp_355 + tmp_359;
      real_t tmp_361 = tmp_350 + tmp_353 + tmp_358;
      real_t tmp_362 = tmp_349 + tmp_352 + tmp_357;
      real_t tmp_363 = 0.042507265838595799 * tmp_58 *
                       ( tmp_11 * ( tmp_362 - 1.0 / 4.0 ) + tmp_3 * ( tmp_361 - 1.0 / 4.0 ) + tmp_5 * ( tmp_360 - 1.0 / 4.0 ) );
      real_t tmp_364 = tmp_19 * ( 0.1711304259088916 * tmp_21 + 0.041227165399737475 * tmp_22 + tmp_23 );
      real_t tmp_365 = tmp_364 * tmp_6;
      real_t tmp_366 = tmp_26 * tmp_364;
      real_t tmp_367 = tmp_19 * ( 0.1711304259088916 * tmp_30 + 0.041227165399737475 * tmp_31 + tmp_32 );
      real_t tmp_368 = tmp_28 * tmp_367;
      real_t tmp_369 = tmp_35 * tmp_367;
      real_t tmp_370 = tmp_364 * tmp_37;
      real_t tmp_371 = tmp_367 * tmp_39;
      real_t tmp_372 = tmp_19 * ( 0.1711304259088916 * tmp_43 + 0.041227165399737475 * tmp_44 + tmp_45 );
      real_t tmp_373 = tmp_372 * tmp_41;
      real_t tmp_374 = tmp_372 * tmp_48;
      real_t tmp_375 = tmp_372 * tmp_50;
      real_t tmp_376 = tmp_370 + tmp_371 + tmp_375;
      real_t tmp_377 = tmp_366 + tmp_369 + tmp_374;
      real_t tmp_378 = tmp_365 + tmp_368 + tmp_373;
      real_t tmp_379 = 0.019202922745021479 * tmp_58 *
                       ( tmp_11 * ( tmp_378 - 1.0 / 4.0 ) + tmp_3 * ( tmp_377 - 1.0 / 4.0 ) + tmp_5 * ( tmp_376 - 1.0 / 4.0 ) );
      real_t a_0_0 = tmp_107 * ( -tmp_101 - tmp_102 - tmp_103 - tmp_93 - tmp_94 - tmp_96 - tmp_97 - tmp_98 - tmp_99 + 1 ) +
                     tmp_123 * ( -tmp_109 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_115 - tmp_117 - tmp_118 - tmp_119 + 1 ) +
                     tmp_139 * ( -tmp_125 - tmp_126 - tmp_128 - tmp_129 - tmp_130 - tmp_131 - tmp_133 - tmp_134 - tmp_135 + 1 ) +
                     tmp_155 * ( -tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1 ) +
                     tmp_171 * ( -tmp_157 - tmp_158 - tmp_160 - tmp_161 - tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 + 1 ) +
                     tmp_187 * ( -tmp_173 - tmp_174 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 + 1 ) +
                     tmp_203 * ( -tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_197 - tmp_198 - tmp_199 + 1 ) +
                     tmp_219 * ( -tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 - tmp_211 - tmp_213 - tmp_214 - tmp_215 + 1 ) +
                     tmp_235 * ( -tmp_221 - tmp_222 - tmp_224 - tmp_225 - tmp_226 - tmp_227 - tmp_229 - tmp_230 - tmp_231 + 1 ) +
                     tmp_251 * ( -tmp_237 - tmp_238 - tmp_240 - tmp_241 - tmp_242 - tmp_243 - tmp_245 - tmp_246 - tmp_247 + 1 ) +
                     tmp_267 * ( -tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1 ) +
                     tmp_283 * ( -tmp_269 - tmp_270 - tmp_272 - tmp_273 - tmp_274 - tmp_275 - tmp_277 - tmp_278 - tmp_279 + 1 ) +
                     tmp_299 * ( -tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1 ) +
                     tmp_315 * ( -tmp_301 - tmp_302 - tmp_304 - tmp_305 - tmp_306 - tmp_307 - tmp_309 - tmp_310 - tmp_311 + 1 ) +
                     tmp_331 * ( -tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1 ) +
                     tmp_347 * ( -tmp_333 - tmp_334 - tmp_336 - tmp_337 - tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 + 1 ) +
                     tmp_363 * ( -tmp_349 - tmp_350 - tmp_352 - tmp_353 - tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 + 1 ) +
                     tmp_379 * ( -tmp_365 - tmp_366 - tmp_368 - tmp_369 - tmp_370 - tmp_371 - tmp_373 - tmp_374 - tmp_375 + 1 ) +
                     tmp_59 * ( -tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1 ) +
                     tmp_75 * ( -tmp_61 - tmp_62 - tmp_64 - tmp_65 - tmp_66 - tmp_67 - tmp_69 - tmp_70 - tmp_71 + 1 ) +
                     tmp_91 * ( -tmp_77 - tmp_78 - tmp_80 - tmp_81 - tmp_82 - tmp_83 - tmp_85 - tmp_86 - tmp_87 + 1 );
      real_t a_1_0 = tmp_104 * tmp_107 + tmp_120 * tmp_123 + tmp_136 * tmp_139 + tmp_152 * tmp_155 + tmp_168 * tmp_171 +
                     tmp_184 * tmp_187 + tmp_200 * tmp_203 + tmp_216 * tmp_219 + tmp_232 * tmp_235 + tmp_248 * tmp_251 +
                     tmp_264 * tmp_267 + tmp_280 * tmp_283 + tmp_296 * tmp_299 + tmp_312 * tmp_315 + tmp_328 * tmp_331 +
                     tmp_344 * tmp_347 + tmp_360 * tmp_363 + tmp_376 * tmp_379 + tmp_52 * tmp_59 + tmp_72 * tmp_75 +
                     tmp_88 * tmp_91;
      real_t a_2_0 = tmp_105 * tmp_107 + tmp_121 * tmp_123 + tmp_137 * tmp_139 + tmp_153 * tmp_155 + tmp_169 * tmp_171 +
                     tmp_185 * tmp_187 + tmp_201 * tmp_203 + tmp_217 * tmp_219 + tmp_233 * tmp_235 + tmp_249 * tmp_251 +
                     tmp_265 * tmp_267 + tmp_281 * tmp_283 + tmp_297 * tmp_299 + tmp_313 * tmp_315 + tmp_329 * tmp_331 +
                     tmp_345 * tmp_347 + tmp_361 * tmp_363 + tmp_377 * tmp_379 + tmp_53 * tmp_59 + tmp_73 * tmp_75 +
                     tmp_89 * tmp_91;
      real_t a_3_0 = tmp_106 * tmp_107 + tmp_122 * tmp_123 + tmp_138 * tmp_139 + tmp_154 * tmp_155 + tmp_170 * tmp_171 +
                     tmp_186 * tmp_187 + tmp_202 * tmp_203 + tmp_218 * tmp_219 + tmp_234 * tmp_235 + tmp_250 * tmp_251 +
                     tmp_266 * tmp_267 + tmp_282 * tmp_283 + tmp_298 * tmp_299 + tmp_314 * tmp_315 + tmp_330 * tmp_331 +
                     tmp_346 * tmp_347 + tmp_362 * tmp_363 + tmp_378 * tmp_379 + tmp_54 * tmp_59 + tmp_74 * tmp_75 +
                     tmp_90 * tmp_91;
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
      elMat( 3, 0 ) = a_3_0;
   }
};

class EGIIPGVectorLaplaceFormP1E_2 : public hyteg::dg::DGForm
{
 protected:
   void integrateVolume2D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( 2, testDegree ), trialBasis.numDoFsPerElement( 2, trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t a_0_0  = 0;
      real_t a_1_0  = 0;
      real_t a_2_0  = 0;
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >& coordsElement,
                                       const std::vector< Point3D >& coordsFacet,
                                       const Point3D&                oppositeVertex,
                                       const Point3D&                outwardNormal,
                                       const DGBasisInfo&            trialBasis,
                                       const DGBasisInfo&            testBasis,
                                       int                           trialDegree,
                                       int                           testDegree,
                                       MatrixXr&                     elMat ) const override
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
      real_t a_1_0  = 0;
      real_t a_2_0  = 0;
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >& coordsElementInner,
                                          const std::vector< Point3D >& coordsElementOuter,
                                          const std::vector< Point3D >& coordsFacet,
                                          const Point3D&                oppositeVertexInnerElement,
                                          const Point3D&                oppositeVertexOuterElement,
                                          const Point3D&                outwardNormal,
                                          const DGBasisInfo&            trialBasis,
                                          const DGBasisInfo&            testBasis,
                                          int                           trialDegree,
                                          int                           testDegree,
                                          MatrixXr&                     elMat ) const override
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
      real_t a_1_0  = 0;
      real_t a_2_0  = 0;
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                                   const std::vector< Point3D >& coordsFacet,
                                                   const Point3D&                oppositeVertex,
                                                   const Point3D&                outwardNormal,
                                                   const DGBasisInfo&            trialBasis,
                                                   const DGBasisInfo&            testBasis,
                                                   int                           trialDegree,
                                                   int                           testDegree,
                                                   MatrixXr&                     elMat ) const override
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
      real_t a_1_0  = 0;
      real_t a_2_0  = 0;
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
   }

   void integrateRHSDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                         const std::vector< Point3D >& coordsFacet,
                                         const Point3D&                oppositeVertex,
                                         const Point3D&                outwardNormal,
                                         const DGBasisInfo&            basis,
                                         int                           degree,
                                         MatrixXr&                     elMat ) const override
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
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                         const std::vector< Point3D >& coordsFacet,
                                         const Point3D&                oppositeVertex,
                                         const Point3D&                outwardNormal,
                                         const DGBasisInfo&            basis,
                                         int                           degree,
                                         MatrixXr&                     elMat ) const override
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
   void integrateVolume3D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_0;
      real_t tmp_6  = p_affine_1_1 + tmp_2;
      real_t tmp_7  = tmp_5 * tmp_6;
      real_t tmp_8  = -p_affine_0_2;
      real_t tmp_9  = p_affine_3_2 + tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_8;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11 * tmp_5;
      real_t tmp_13 = p_affine_2_2 + tmp_8;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = tmp_14 * tmp_6;
      real_t tmp_16 = tmp_1 * tmp_11;
      real_t tmp_17 = tmp_14 * tmp_3;
      real_t tmp_18 =
          1.0 / ( tmp_10 * tmp_12 - tmp_10 * tmp_17 + tmp_13 * tmp_15 - tmp_13 * tmp_16 + tmp_4 * tmp_9 - tmp_7 * tmp_9 );
      real_t tmp_19 = tmp_18 * ( tmp_4 - tmp_7 );
      real_t tmp_20 = tmp_18 * ( tmp_15 - tmp_16 );
      real_t tmp_21 = tmp_18 * ( tmp_12 - tmp_17 );
      real_t tmp_22 = tmp_10 * tmp_21 + tmp_13 * tmp_20 + tmp_19 * tmp_9;
      real_t tmp_23 = tmp_18 * ( -tmp_1 * tmp_13 + tmp_10 * tmp_5 );
      real_t tmp_24 = tmp_18 * ( tmp_1 * tmp_9 - tmp_10 * tmp_14 );
      real_t tmp_25 = tmp_18 * ( tmp_13 * tmp_14 - tmp_5 * tmp_9 );
      real_t tmp_26 = tmp_10 * tmp_25 + tmp_13 * tmp_24 + tmp_23 * tmp_9;
      real_t tmp_27 = tmp_18 * ( -tmp_10 * tmp_3 + tmp_13 * tmp_6 );
      real_t tmp_28 = tmp_18 * ( tmp_10 * tmp_11 - tmp_6 * tmp_9 );
      real_t tmp_29 = tmp_18 * ( -tmp_11 * tmp_13 + tmp_3 * tmp_9 );
      real_t tmp_30 = tmp_10 * tmp_29 + tmp_13 * tmp_28 + tmp_27 * tmp_9;
      real_t tmp_31 = p_affine_0_0 * p_affine_1_1;
      real_t tmp_32 = p_affine_0_0 * p_affine_1_2;
      real_t tmp_33 = p_affine_2_1 * p_affine_3_2;
      real_t tmp_34 = p_affine_0_1 * p_affine_1_0;
      real_t tmp_35 = p_affine_0_1 * p_affine_1_2;
      real_t tmp_36 = p_affine_2_2 * p_affine_3_0;
      real_t tmp_37 = p_affine_0_2 * p_affine_1_0;
      real_t tmp_38 = p_affine_0_2 * p_affine_1_1;
      real_t tmp_39 = p_affine_2_0 * p_affine_3_1;
      real_t tmp_40 = p_affine_2_2 * p_affine_3_1;
      real_t tmp_41 = p_affine_2_0 * p_affine_3_2;
      real_t tmp_42 = p_affine_2_1 * p_affine_3_0;
      real_t tmp_43 = std::abs( p_affine_0_0 * tmp_33 - p_affine_0_0 * tmp_40 + p_affine_0_1 * tmp_36 - p_affine_0_1 * tmp_41 +
                                p_affine_0_2 * tmp_39 - p_affine_0_2 * tmp_42 - p_affine_1_0 * tmp_33 + p_affine_1_0 * tmp_40 -
                                p_affine_1_1 * tmp_36 + p_affine_1_1 * tmp_41 - p_affine_1_2 * tmp_39 + p_affine_1_2 * tmp_42 +
                                p_affine_2_0 * tmp_35 - p_affine_2_0 * tmp_38 - p_affine_2_1 * tmp_32 + p_affine_2_1 * tmp_37 +
                                p_affine_2_2 * tmp_31 - p_affine_2_2 * tmp_34 - p_affine_3_0 * tmp_35 + p_affine_3_0 * tmp_38 +
                                p_affine_3_1 * tmp_32 - p_affine_3_1 * tmp_37 - p_affine_3_2 * tmp_31 + p_affine_3_2 * tmp_34 );
      real_t tmp_44 = tmp_43 * ( tmp_22 * ( -tmp_19 - tmp_20 - tmp_21 ) + tmp_26 * ( -tmp_23 - tmp_24 - tmp_25 ) +
                                 tmp_30 * ( -tmp_27 - tmp_28 - tmp_29 ) );
      real_t tmp_45 = tmp_43 * ( tmp_21 * tmp_22 + tmp_25 * tmp_26 + tmp_29 * tmp_30 );
      real_t tmp_46 = tmp_43 * ( tmp_20 * tmp_22 + tmp_24 * tmp_26 + tmp_28 * tmp_30 );
      real_t tmp_47 = tmp_43 * ( tmp_19 * tmp_22 + tmp_23 * tmp_26 + tmp_27 * tmp_30 );
      real_t a_0_0  = 0.1666666666666668 * tmp_44;
      real_t a_1_0  = 0.1666666666666668 * tmp_45;
      real_t a_2_0  = 0.1666666666666668 * tmp_46;
      real_t a_3_0  = 0.1666666666666668 * tmp_47;
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
      elMat( 3, 0 ) = a_3_0;
   }

   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                               const std::vector< Point3D >& coordsFacet,
                               const Point3D&,
                               const Point3D&     outwardNormal,
                               const DGBasisInfo& trialBasis,
                               const DGBasisInfo& testBasis,
                               int                trialDegree,
                               int                testDegree,
                               MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = p_affine_2_0 + tmp_0;
      real_t tmp_5  = p_affine_1_1 + tmp_2;
      real_t tmp_6  = tmp_1 * tmp_3 - tmp_4 * tmp_5;
      real_t tmp_7  = -p_affine_0_2;
      real_t tmp_8  = p_affine_3_2 + tmp_7;
      real_t tmp_9  = tmp_3 * tmp_8;
      real_t tmp_10 = p_affine_3_1 + tmp_2;
      real_t tmp_11 = p_affine_1_2 + tmp_7;
      real_t tmp_12 = tmp_10 * tmp_11;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_14 * tmp_5;
      real_t tmp_16 = tmp_10 * tmp_14;
      real_t tmp_17 = tmp_5 * tmp_8;
      real_t tmp_18 = tmp_11 * tmp_3;
      real_t tmp_19 =
          1.0 / ( -tmp_1 * tmp_16 + tmp_1 * tmp_9 + tmp_12 * tmp_4 + tmp_13 * tmp_15 - tmp_13 * tmp_18 - tmp_17 * tmp_4 );
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19 * ( 0.031405749086161582 * tmp_21 + 0.93718850182767688 * tmp_22 + tmp_23 );
      real_t tmp_25 = tmp_24 * tmp_6;
      real_t tmp_26 = -tmp_1 * tmp_10 + tmp_13 * tmp_5;
      real_t tmp_27 = tmp_24 * tmp_26;
      real_t tmp_28 = -tmp_1 * tmp_14 + tmp_11 * tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19 * ( 0.031405749086161582 * tmp_30 + 0.93718850182767688 * tmp_31 + tmp_32 );
      real_t tmp_34 = tmp_28 * tmp_33;
      real_t tmp_35 = tmp_1 * tmp_8 - tmp_11 * tmp_13;
      real_t tmp_36 = tmp_33 * tmp_35;
      real_t tmp_37 = tmp_10 * tmp_4 - tmp_13 * tmp_3;
      real_t tmp_38 = tmp_24 * tmp_37;
      real_t tmp_39 = tmp_13 * tmp_14 - tmp_4 * tmp_8;
      real_t tmp_40 = tmp_33 * tmp_39;
      real_t tmp_41 = tmp_15 - tmp_18;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19 * ( 0.031405749086161582 * tmp_43 + 0.93718850182767688 * tmp_44 + tmp_45 );
      real_t tmp_47 = tmp_41 * tmp_46;
      real_t tmp_48 = tmp_12 - tmp_17;
      real_t tmp_49 = tmp_46 * tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46 * tmp_50;
      real_t tmp_52 = -tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1;
      real_t tmp_53 = tmp_11 * tmp_19;
      real_t tmp_54 = tmp_14 * tmp_19;
      real_t tmp_55 = tmp_19 * tmp_8;
      real_t tmp_56 = 0.5 * p_affine_13_0 * ( tmp_41 * tmp_55 + tmp_48 * tmp_54 + tmp_50 * tmp_53 ) +
                      0.5 * p_affine_13_1 * ( tmp_28 * tmp_55 + tmp_35 * tmp_54 + tmp_39 * tmp_53 ) +
                      0.5 * p_affine_13_2 * ( tmp_26 * tmp_54 + tmp_37 * tmp_53 + tmp_55 * tmp_6 );
      real_t tmp_57 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_58 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_59 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_60 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_61 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_62 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_63 = ( std::abs( tmp_22 * tmp_60 - tmp_31 * tmp_62 ) * std::abs( tmp_22 * tmp_60 - tmp_31 * tmp_62 ) ) +
                      ( std::abs( tmp_22 * tmp_61 - tmp_44 * tmp_62 ) * std::abs( tmp_22 * tmp_61 - tmp_44 * tmp_62 ) ) +
                      ( std::abs( tmp_31 * tmp_61 - tmp_44 * tmp_60 ) * std::abs( tmp_31 * tmp_61 - tmp_44 * tmp_60 ) );
      real_t tmp_64 = 1.0 * std::pow( tmp_63, -0.25 );
      real_t tmp_65 =
          tmp_64 * ( tmp_11 * ( tmp_57 - 1.0 / 4.0 ) + tmp_14 * ( tmp_58 - 1.0 / 4.0 ) + tmp_8 * ( tmp_59 - 1.0 / 4.0 ) );
      real_t tmp_66 = 1.0 * std::pow( tmp_63, 1.0 / 2.0 );
      real_t tmp_67 = 0.0068572537431980923 * tmp_66;
      real_t tmp_68 = tmp_19 * ( 0.19601935860219369 * tmp_21 + 0.60796128279561268 * tmp_22 + tmp_23 );
      real_t tmp_69 = tmp_6 * tmp_68;
      real_t tmp_70 = tmp_26 * tmp_68;
      real_t tmp_71 = tmp_19 * ( 0.19601935860219369 * tmp_30 + 0.60796128279561268 * tmp_31 + tmp_32 );
      real_t tmp_72 = tmp_28 * tmp_71;
      real_t tmp_73 = tmp_35 * tmp_71;
      real_t tmp_74 = tmp_37 * tmp_68;
      real_t tmp_75 = tmp_39 * tmp_71;
      real_t tmp_76 = tmp_19 * ( 0.19601935860219369 * tmp_43 + 0.60796128279561268 * tmp_44 + tmp_45 );
      real_t tmp_77 = tmp_41 * tmp_76;
      real_t tmp_78 = tmp_48 * tmp_76;
      real_t tmp_79 = tmp_50 * tmp_76;
      real_t tmp_80 = -tmp_69 - tmp_70 - tmp_72 - tmp_73 - tmp_74 - tmp_75 - tmp_77 - tmp_78 - tmp_79 + 1;
      real_t tmp_81 = tmp_74 + tmp_75 + tmp_79;
      real_t tmp_82 = tmp_70 + tmp_73 + tmp_78;
      real_t tmp_83 = tmp_69 + tmp_72 + tmp_77;
      real_t tmp_84 =
          tmp_64 * ( tmp_11 * ( tmp_81 - 1.0 / 4.0 ) + tmp_14 * ( tmp_82 - 1.0 / 4.0 ) + tmp_8 * ( tmp_83 - 1.0 / 4.0 ) );
      real_t tmp_85  = 0.037198804536718075 * tmp_66;
      real_t tmp_86  = tmp_19 * ( 0.37605877282253791 * tmp_21 + 0.039308471900058539 * tmp_22 + tmp_23 );
      real_t tmp_87  = tmp_6 * tmp_86;
      real_t tmp_88  = tmp_26 * tmp_86;
      real_t tmp_89  = tmp_19 * ( 0.37605877282253791 * tmp_30 + 0.039308471900058539 * tmp_31 + tmp_32 );
      real_t tmp_90  = tmp_28 * tmp_89;
      real_t tmp_91  = tmp_35 * tmp_89;
      real_t tmp_92  = tmp_37 * tmp_86;
      real_t tmp_93  = tmp_39 * tmp_89;
      real_t tmp_94  = tmp_19 * ( 0.37605877282253791 * tmp_43 + 0.039308471900058539 * tmp_44 + tmp_45 );
      real_t tmp_95  = tmp_41 * tmp_94;
      real_t tmp_96  = tmp_48 * tmp_94;
      real_t tmp_97  = tmp_50 * tmp_94;
      real_t tmp_98  = -tmp_87 - tmp_88 - tmp_90 - tmp_91 - tmp_92 - tmp_93 - tmp_95 - tmp_96 - tmp_97 + 1;
      real_t tmp_99  = tmp_92 + tmp_93 + tmp_97;
      real_t tmp_100 = tmp_88 + tmp_91 + tmp_96;
      real_t tmp_101 = tmp_87 + tmp_90 + tmp_95;
      real_t tmp_102 =
          tmp_64 * ( tmp_11 * ( tmp_99 - 1.0 / 4.0 ) + tmp_14 * ( tmp_100 - 1.0 / 4.0 ) + tmp_8 * ( tmp_101 - 1.0 / 4.0 ) );
      real_t tmp_103 = 0.020848748529055869 * tmp_66;
      real_t tmp_104 = tmp_19 * ( 0.78764240869137092 * tmp_21 + 0.1711304259088916 * tmp_22 + tmp_23 );
      real_t tmp_105 = tmp_104 * tmp_6;
      real_t tmp_106 = tmp_104 * tmp_26;
      real_t tmp_107 = tmp_19 * ( 0.78764240869137092 * tmp_30 + 0.1711304259088916 * tmp_31 + tmp_32 );
      real_t tmp_108 = tmp_107 * tmp_28;
      real_t tmp_109 = tmp_107 * tmp_35;
      real_t tmp_110 = tmp_104 * tmp_37;
      real_t tmp_111 = tmp_107 * tmp_39;
      real_t tmp_112 = tmp_19 * ( 0.78764240869137092 * tmp_43 + 0.1711304259088916 * tmp_44 + tmp_45 );
      real_t tmp_113 = tmp_112 * tmp_41;
      real_t tmp_114 = tmp_112 * tmp_48;
      real_t tmp_115 = tmp_112 * tmp_50;
      real_t tmp_116 = -tmp_105 - tmp_106 - tmp_108 - tmp_109 - tmp_110 - tmp_111 - tmp_113 - tmp_114 - tmp_115 + 1;
      real_t tmp_117 = tmp_110 + tmp_111 + tmp_115;
      real_t tmp_118 = tmp_106 + tmp_109 + tmp_114;
      real_t tmp_119 = tmp_105 + tmp_108 + tmp_113;
      real_t tmp_120 =
          tmp_64 * ( tmp_11 * ( tmp_117 - 1.0 / 4.0 ) + tmp_14 * ( tmp_118 - 1.0 / 4.0 ) + tmp_8 * ( tmp_119 - 1.0 / 4.0 ) );
      real_t tmp_121 = 0.019202922745021479 * tmp_66;
      real_t tmp_122 = tmp_19 * ( 0.58463275527740355 * tmp_21 + 0.37605877282253791 * tmp_22 + tmp_23 );
      real_t tmp_123 = tmp_122 * tmp_6;
      real_t tmp_124 = tmp_122 * tmp_26;
      real_t tmp_125 = tmp_19 * ( 0.58463275527740355 * tmp_30 + 0.37605877282253791 * tmp_31 + tmp_32 );
      real_t tmp_126 = tmp_125 * tmp_28;
      real_t tmp_127 = tmp_125 * tmp_35;
      real_t tmp_128 = tmp_122 * tmp_37;
      real_t tmp_129 = tmp_125 * tmp_39;
      real_t tmp_130 = tmp_19 * ( 0.58463275527740355 * tmp_43 + 0.37605877282253791 * tmp_44 + tmp_45 );
      real_t tmp_131 = tmp_130 * tmp_41;
      real_t tmp_132 = tmp_130 * tmp_48;
      real_t tmp_133 = tmp_130 * tmp_50;
      real_t tmp_134 = -tmp_123 - tmp_124 - tmp_126 - tmp_127 - tmp_128 - tmp_129 - tmp_131 - tmp_132 - tmp_133 + 1;
      real_t tmp_135 = tmp_128 + tmp_129 + tmp_133;
      real_t tmp_136 = tmp_124 + tmp_127 + tmp_132;
      real_t tmp_137 = tmp_123 + tmp_126 + tmp_131;
      real_t tmp_138 =
          tmp_64 * ( tmp_11 * ( tmp_135 - 1.0 / 4.0 ) + tmp_14 * ( tmp_136 - 1.0 / 4.0 ) + tmp_8 * ( tmp_137 - 1.0 / 4.0 ) );
      real_t tmp_139 = 0.020848748529055869 * tmp_66;
      real_t tmp_140 = tmp_19 * ( 0.041227165399737475 * tmp_21 + 0.78764240869137092 * tmp_22 + tmp_23 );
      real_t tmp_141 = tmp_140 * tmp_6;
      real_t tmp_142 = tmp_140 * tmp_26;
      real_t tmp_143 = tmp_19 * ( 0.041227165399737475 * tmp_30 + 0.78764240869137092 * tmp_31 + tmp_32 );
      real_t tmp_144 = tmp_143 * tmp_28;
      real_t tmp_145 = tmp_143 * tmp_35;
      real_t tmp_146 = tmp_140 * tmp_37;
      real_t tmp_147 = tmp_143 * tmp_39;
      real_t tmp_148 = tmp_19 * ( 0.041227165399737475 * tmp_43 + 0.78764240869137092 * tmp_44 + tmp_45 );
      real_t tmp_149 = tmp_148 * tmp_41;
      real_t tmp_150 = tmp_148 * tmp_48;
      real_t tmp_151 = tmp_148 * tmp_50;
      real_t tmp_152 = -tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1;
      real_t tmp_153 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_154 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_155 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_156 =
          tmp_64 * ( tmp_11 * ( tmp_153 - 1.0 / 4.0 ) + tmp_14 * ( tmp_154 - 1.0 / 4.0 ) + tmp_8 * ( tmp_155 - 1.0 / 4.0 ) );
      real_t tmp_157 = 0.019202922745021479 * tmp_66;
      real_t tmp_158 = tmp_19 * ( 0.039308471900058539 * tmp_21 + 0.58463275527740355 * tmp_22 + tmp_23 );
      real_t tmp_159 = tmp_158 * tmp_6;
      real_t tmp_160 = tmp_158 * tmp_26;
      real_t tmp_161 = tmp_19 * ( 0.039308471900058539 * tmp_30 + 0.58463275527740355 * tmp_31 + tmp_32 );
      real_t tmp_162 = tmp_161 * tmp_28;
      real_t tmp_163 = tmp_161 * tmp_35;
      real_t tmp_164 = tmp_158 * tmp_37;
      real_t tmp_165 = tmp_161 * tmp_39;
      real_t tmp_166 = tmp_19 * ( 0.039308471900058539 * tmp_43 + 0.58463275527740355 * tmp_44 + tmp_45 );
      real_t tmp_167 = tmp_166 * tmp_41;
      real_t tmp_168 = tmp_166 * tmp_48;
      real_t tmp_169 = tmp_166 * tmp_50;
      real_t tmp_170 = -tmp_159 - tmp_160 - tmp_162 - tmp_163 - tmp_164 - tmp_165 - tmp_167 - tmp_168 - tmp_169 + 1;
      real_t tmp_171 = tmp_164 + tmp_165 + tmp_169;
      real_t tmp_172 = tmp_160 + tmp_163 + tmp_168;
      real_t tmp_173 = tmp_159 + tmp_162 + tmp_167;
      real_t tmp_174 =
          tmp_64 * ( tmp_11 * ( tmp_171 - 1.0 / 4.0 ) + tmp_14 * ( tmp_172 - 1.0 / 4.0 ) + tmp_8 * ( tmp_173 - 1.0 / 4.0 ) );
      real_t tmp_175 = 0.020848748529055869 * tmp_66;
      real_t tmp_176 = tmp_19 * ( 0.78764240869137092 * tmp_21 + 0.041227165399737475 * tmp_22 + tmp_23 );
      real_t tmp_177 = tmp_176 * tmp_6;
      real_t tmp_178 = tmp_176 * tmp_26;
      real_t tmp_179 = tmp_19 * ( 0.78764240869137092 * tmp_30 + 0.041227165399737475 * tmp_31 + tmp_32 );
      real_t tmp_180 = tmp_179 * tmp_28;
      real_t tmp_181 = tmp_179 * tmp_35;
      real_t tmp_182 = tmp_176 * tmp_37;
      real_t tmp_183 = tmp_179 * tmp_39;
      real_t tmp_184 = tmp_19 * ( 0.78764240869137092 * tmp_43 + 0.041227165399737475 * tmp_44 + tmp_45 );
      real_t tmp_185 = tmp_184 * tmp_41;
      real_t tmp_186 = tmp_184 * tmp_48;
      real_t tmp_187 = tmp_184 * tmp_50;
      real_t tmp_188 = -tmp_177 - tmp_178 - tmp_180 - tmp_181 - tmp_182 - tmp_183 - tmp_185 - tmp_186 - tmp_187 + 1;
      real_t tmp_189 = tmp_182 + tmp_183 + tmp_187;
      real_t tmp_190 = tmp_178 + tmp_181 + tmp_186;
      real_t tmp_191 = tmp_177 + tmp_180 + tmp_185;
      real_t tmp_192 =
          tmp_64 * ( tmp_11 * ( tmp_189 - 1.0 / 4.0 ) + tmp_14 * ( tmp_190 - 1.0 / 4.0 ) + tmp_8 * ( tmp_191 - 1.0 / 4.0 ) );
      real_t tmp_193 = 0.019202922745021479 * tmp_66;
      real_t tmp_194 = tmp_19 * ( 0.58463275527740355 * tmp_21 + 0.039308471900058539 * tmp_22 + tmp_23 );
      real_t tmp_195 = tmp_194 * tmp_6;
      real_t tmp_196 = tmp_194 * tmp_26;
      real_t tmp_197 = tmp_19 * ( 0.58463275527740355 * tmp_30 + 0.039308471900058539 * tmp_31 + tmp_32 );
      real_t tmp_198 = tmp_197 * tmp_28;
      real_t tmp_199 = tmp_197 * tmp_35;
      real_t tmp_200 = tmp_194 * tmp_37;
      real_t tmp_201 = tmp_197 * tmp_39;
      real_t tmp_202 = tmp_19 * ( 0.58463275527740355 * tmp_43 + 0.039308471900058539 * tmp_44 + tmp_45 );
      real_t tmp_203 = tmp_202 * tmp_41;
      real_t tmp_204 = tmp_202 * tmp_48;
      real_t tmp_205 = tmp_202 * tmp_50;
      real_t tmp_206 = -tmp_195 - tmp_196 - tmp_198 - tmp_199 - tmp_200 - tmp_201 - tmp_203 - tmp_204 - tmp_205 + 1;
      real_t tmp_207 = tmp_200 + tmp_201 + tmp_205;
      real_t tmp_208 = tmp_196 + tmp_199 + tmp_204;
      real_t tmp_209 = tmp_195 + tmp_198 + tmp_203;
      real_t tmp_210 =
          tmp_64 * ( tmp_11 * ( tmp_207 - 1.0 / 4.0 ) + tmp_14 * ( tmp_208 - 1.0 / 4.0 ) + tmp_8 * ( tmp_209 - 1.0 / 4.0 ) );
      real_t tmp_211 = 0.020848748529055869 * tmp_66;
      real_t tmp_212 = tmp_19 * ( 0.1711304259088916 * tmp_21 + 0.78764240869137092 * tmp_22 + tmp_23 );
      real_t tmp_213 = tmp_212 * tmp_6;
      real_t tmp_214 = tmp_212 * tmp_26;
      real_t tmp_215 = tmp_19 * ( 0.1711304259088916 * tmp_30 + 0.78764240869137092 * tmp_31 + tmp_32 );
      real_t tmp_216 = tmp_215 * tmp_28;
      real_t tmp_217 = tmp_215 * tmp_35;
      real_t tmp_218 = tmp_212 * tmp_37;
      real_t tmp_219 = tmp_215 * tmp_39;
      real_t tmp_220 = tmp_19 * ( 0.1711304259088916 * tmp_43 + 0.78764240869137092 * tmp_44 + tmp_45 );
      real_t tmp_221 = tmp_220 * tmp_41;
      real_t tmp_222 = tmp_220 * tmp_48;
      real_t tmp_223 = tmp_220 * tmp_50;
      real_t tmp_224 = -tmp_213 - tmp_214 - tmp_216 - tmp_217 - tmp_218 - tmp_219 - tmp_221 - tmp_222 - tmp_223 + 1;
      real_t tmp_225 = tmp_218 + tmp_219 + tmp_223;
      real_t tmp_226 = tmp_214 + tmp_217 + tmp_222;
      real_t tmp_227 = tmp_213 + tmp_216 + tmp_221;
      real_t tmp_228 =
          tmp_64 * ( tmp_11 * ( tmp_225 - 1.0 / 4.0 ) + tmp_14 * ( tmp_226 - 1.0 / 4.0 ) + tmp_8 * ( tmp_227 - 1.0 / 4.0 ) );
      real_t tmp_229 = 0.019202922745021479 * tmp_66;
      real_t tmp_230 = tmp_19 * ( 0.37605877282253791 * tmp_21 + 0.58463275527740355 * tmp_22 + tmp_23 );
      real_t tmp_231 = tmp_230 * tmp_6;
      real_t tmp_232 = tmp_230 * tmp_26;
      real_t tmp_233 = tmp_19 * ( 0.37605877282253791 * tmp_30 + 0.58463275527740355 * tmp_31 + tmp_32 );
      real_t tmp_234 = tmp_233 * tmp_28;
      real_t tmp_235 = tmp_233 * tmp_35;
      real_t tmp_236 = tmp_230 * tmp_37;
      real_t tmp_237 = tmp_233 * tmp_39;
      real_t tmp_238 = tmp_19 * ( 0.37605877282253791 * tmp_43 + 0.58463275527740355 * tmp_44 + tmp_45 );
      real_t tmp_239 = tmp_238 * tmp_41;
      real_t tmp_240 = tmp_238 * tmp_48;
      real_t tmp_241 = tmp_238 * tmp_50;
      real_t tmp_242 = -tmp_231 - tmp_232 - tmp_234 - tmp_235 - tmp_236 - tmp_237 - tmp_239 - tmp_240 - tmp_241 + 1;
      real_t tmp_243 = tmp_236 + tmp_237 + tmp_241;
      real_t tmp_244 = tmp_232 + tmp_235 + tmp_240;
      real_t tmp_245 = tmp_231 + tmp_234 + tmp_239;
      real_t tmp_246 =
          tmp_64 * ( tmp_11 * ( tmp_243 - 1.0 / 4.0 ) + tmp_14 * ( tmp_244 - 1.0 / 4.0 ) + tmp_8 * ( tmp_245 - 1.0 / 4.0 ) );
      real_t tmp_247 = 0.020848748529055869 * tmp_66;
      real_t tmp_248 = tmp_19 * ( 0.041227165399737475 * tmp_21 + 0.1711304259088916 * tmp_22 + tmp_23 );
      real_t tmp_249 = tmp_248 * tmp_6;
      real_t tmp_250 = tmp_248 * tmp_26;
      real_t tmp_251 = tmp_19 * ( 0.041227165399737475 * tmp_30 + 0.1711304259088916 * tmp_31 + tmp_32 );
      real_t tmp_252 = tmp_251 * tmp_28;
      real_t tmp_253 = tmp_251 * tmp_35;
      real_t tmp_254 = tmp_248 * tmp_37;
      real_t tmp_255 = tmp_251 * tmp_39;
      real_t tmp_256 = tmp_19 * ( 0.041227165399737475 * tmp_43 + 0.1711304259088916 * tmp_44 + tmp_45 );
      real_t tmp_257 = tmp_256 * tmp_41;
      real_t tmp_258 = tmp_256 * tmp_48;
      real_t tmp_259 = tmp_256 * tmp_50;
      real_t tmp_260 = -tmp_249 - tmp_250 - tmp_252 - tmp_253 - tmp_254 - tmp_255 - tmp_257 - tmp_258 - tmp_259 + 1;
      real_t tmp_261 = tmp_254 + tmp_255 + tmp_259;
      real_t tmp_262 = tmp_250 + tmp_253 + tmp_258;
      real_t tmp_263 = tmp_249 + tmp_252 + tmp_257;
      real_t tmp_264 =
          tmp_64 * ( tmp_11 * ( tmp_261 - 1.0 / 4.0 ) + tmp_14 * ( tmp_262 - 1.0 / 4.0 ) + tmp_8 * ( tmp_263 - 1.0 / 4.0 ) );
      real_t tmp_265 = 0.019202922745021479 * tmp_66;
      real_t tmp_266 = tmp_19 * ( 0.40446199974765351 * tmp_21 + 0.19107600050469298 * tmp_22 + tmp_23 );
      real_t tmp_267 = tmp_266 * tmp_6;
      real_t tmp_268 = tmp_26 * tmp_266;
      real_t tmp_269 = tmp_19 * ( 0.40446199974765351 * tmp_30 + 0.19107600050469298 * tmp_31 + tmp_32 );
      real_t tmp_270 = tmp_269 * tmp_28;
      real_t tmp_271 = tmp_269 * tmp_35;
      real_t tmp_272 = tmp_266 * tmp_37;
      real_t tmp_273 = tmp_269 * tmp_39;
      real_t tmp_274 = tmp_19 * ( 0.40446199974765351 * tmp_43 + 0.19107600050469298 * tmp_44 + tmp_45 );
      real_t tmp_275 = tmp_274 * tmp_41;
      real_t tmp_276 = tmp_274 * tmp_48;
      real_t tmp_277 = tmp_274 * tmp_50;
      real_t tmp_278 = -tmp_267 - tmp_268 - tmp_270 - tmp_271 - tmp_272 - tmp_273 - tmp_275 - tmp_276 - tmp_277 + 1;
      real_t tmp_279 = tmp_272 + tmp_273 + tmp_277;
      real_t tmp_280 = tmp_268 + tmp_271 + tmp_276;
      real_t tmp_281 = tmp_267 + tmp_270 + tmp_275;
      real_t tmp_282 =
          tmp_64 * ( tmp_11 * ( tmp_279 - 1.0 / 4.0 ) + tmp_14 * ( tmp_280 - 1.0 / 4.0 ) + tmp_8 * ( tmp_281 - 1.0 / 4.0 ) );
      real_t tmp_283 = 0.042507265838595799 * tmp_66;
      real_t tmp_284 = tmp_19 * ( 0.039308471900058539 * tmp_21 + 0.37605877282253791 * tmp_22 + tmp_23 );
      real_t tmp_285 = tmp_284 * tmp_6;
      real_t tmp_286 = tmp_26 * tmp_284;
      real_t tmp_287 = tmp_19 * ( 0.039308471900058539 * tmp_30 + 0.37605877282253791 * tmp_31 + tmp_32 );
      real_t tmp_288 = tmp_28 * tmp_287;
      real_t tmp_289 = tmp_287 * tmp_35;
      real_t tmp_290 = tmp_284 * tmp_37;
      real_t tmp_291 = tmp_287 * tmp_39;
      real_t tmp_292 = tmp_19 * ( 0.039308471900058539 * tmp_43 + 0.37605877282253791 * tmp_44 + tmp_45 );
      real_t tmp_293 = tmp_292 * tmp_41;
      real_t tmp_294 = tmp_292 * tmp_48;
      real_t tmp_295 = tmp_292 * tmp_50;
      real_t tmp_296 = -tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1;
      real_t tmp_297 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_298 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_299 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_300 =
          tmp_64 * ( tmp_11 * ( tmp_297 - 1.0 / 4.0 ) + tmp_14 * ( tmp_298 - 1.0 / 4.0 ) + tmp_8 * ( tmp_299 - 1.0 / 4.0 ) );
      real_t tmp_301 = 0.020848748529055869 * tmp_66;
      real_t tmp_302 = tmp_19 * ( 0.93718850182767688 * tmp_21 + 0.031405749086161582 * tmp_22 + tmp_23 );
      real_t tmp_303 = tmp_302 * tmp_6;
      real_t tmp_304 = tmp_26 * tmp_302;
      real_t tmp_305 = tmp_19 * ( 0.93718850182767688 * tmp_30 + 0.031405749086161582 * tmp_31 + tmp_32 );
      real_t tmp_306 = tmp_28 * tmp_305;
      real_t tmp_307 = tmp_305 * tmp_35;
      real_t tmp_308 = tmp_302 * tmp_37;
      real_t tmp_309 = tmp_305 * tmp_39;
      real_t tmp_310 = tmp_19 * ( 0.93718850182767688 * tmp_43 + 0.031405749086161582 * tmp_44 + tmp_45 );
      real_t tmp_311 = tmp_310 * tmp_41;
      real_t tmp_312 = tmp_310 * tmp_48;
      real_t tmp_313 = tmp_310 * tmp_50;
      real_t tmp_314 = -tmp_303 - tmp_304 - tmp_306 - tmp_307 - tmp_308 - tmp_309 - tmp_311 - tmp_312 - tmp_313 + 1;
      real_t tmp_315 = tmp_308 + tmp_309 + tmp_313;
      real_t tmp_316 = tmp_304 + tmp_307 + tmp_312;
      real_t tmp_317 = tmp_303 + tmp_306 + tmp_311;
      real_t tmp_318 =
          tmp_64 * ( tmp_11 * ( tmp_315 - 1.0 / 4.0 ) + tmp_14 * ( tmp_316 - 1.0 / 4.0 ) + tmp_8 * ( tmp_317 - 1.0 / 4.0 ) );
      real_t tmp_319 = 0.0068572537431980923 * tmp_66;
      real_t tmp_320 = tmp_19 * ( 0.60796128279561268 * tmp_21 + 0.19601935860219369 * tmp_22 + tmp_23 );
      real_t tmp_321 = tmp_320 * tmp_6;
      real_t tmp_322 = tmp_26 * tmp_320;
      real_t tmp_323 = tmp_19 * ( 0.60796128279561268 * tmp_30 + 0.19601935860219369 * tmp_31 + tmp_32 );
      real_t tmp_324 = tmp_28 * tmp_323;
      real_t tmp_325 = tmp_323 * tmp_35;
      real_t tmp_326 = tmp_320 * tmp_37;
      real_t tmp_327 = tmp_323 * tmp_39;
      real_t tmp_328 = tmp_19 * ( 0.60796128279561268 * tmp_43 + 0.19601935860219369 * tmp_44 + tmp_45 );
      real_t tmp_329 = tmp_328 * tmp_41;
      real_t tmp_330 = tmp_328 * tmp_48;
      real_t tmp_331 = tmp_328 * tmp_50;
      real_t tmp_332 = -tmp_321 - tmp_322 - tmp_324 - tmp_325 - tmp_326 - tmp_327 - tmp_329 - tmp_330 - tmp_331 + 1;
      real_t tmp_333 = tmp_326 + tmp_327 + tmp_331;
      real_t tmp_334 = tmp_322 + tmp_325 + tmp_330;
      real_t tmp_335 = tmp_321 + tmp_324 + tmp_329;
      real_t tmp_336 =
          tmp_64 * ( tmp_11 * ( tmp_333 - 1.0 / 4.0 ) + tmp_14 * ( tmp_334 - 1.0 / 4.0 ) + tmp_8 * ( tmp_335 - 1.0 / 4.0 ) );
      real_t tmp_337 = 0.037198804536718075 * tmp_66;
      real_t tmp_338 = tmp_19 * ( 0.19107600050469298 * tmp_21 + 0.40446199974765351 * tmp_22 + tmp_23 );
      real_t tmp_339 = tmp_338 * tmp_6;
      real_t tmp_340 = tmp_26 * tmp_338;
      real_t tmp_341 = tmp_19 * ( 0.19107600050469298 * tmp_30 + 0.40446199974765351 * tmp_31 + tmp_32 );
      real_t tmp_342 = tmp_28 * tmp_341;
      real_t tmp_343 = tmp_341 * tmp_35;
      real_t tmp_344 = tmp_338 * tmp_37;
      real_t tmp_345 = tmp_341 * tmp_39;
      real_t tmp_346 = tmp_19 * ( 0.19107600050469298 * tmp_43 + 0.40446199974765351 * tmp_44 + tmp_45 );
      real_t tmp_347 = tmp_346 * tmp_41;
      real_t tmp_348 = tmp_346 * tmp_48;
      real_t tmp_349 = tmp_346 * tmp_50;
      real_t tmp_350 = -tmp_339 - tmp_340 - tmp_342 - tmp_343 - tmp_344 - tmp_345 - tmp_347 - tmp_348 - tmp_349 + 1;
      real_t tmp_351 = tmp_344 + tmp_345 + tmp_349;
      real_t tmp_352 = tmp_340 + tmp_343 + tmp_348;
      real_t tmp_353 = tmp_339 + tmp_342 + tmp_347;
      real_t tmp_354 =
          tmp_64 * ( tmp_11 * ( tmp_351 - 1.0 / 4.0 ) + tmp_14 * ( tmp_352 - 1.0 / 4.0 ) + tmp_8 * ( tmp_353 - 1.0 / 4.0 ) );
      real_t tmp_355 = 0.042507265838595799 * tmp_66;
      real_t tmp_356 = tmp_19 * ( 0.031405749086161582 * tmp_21 + 0.031405749086161582 * tmp_22 + tmp_23 );
      real_t tmp_357 = tmp_356 * tmp_6;
      real_t tmp_358 = tmp_26 * tmp_356;
      real_t tmp_359 = tmp_19 * ( 0.031405749086161582 * tmp_30 + 0.031405749086161582 * tmp_31 + tmp_32 );
      real_t tmp_360 = tmp_28 * tmp_359;
      real_t tmp_361 = tmp_35 * tmp_359;
      real_t tmp_362 = tmp_356 * tmp_37;
      real_t tmp_363 = tmp_359 * tmp_39;
      real_t tmp_364 = tmp_19 * ( 0.031405749086161582 * tmp_43 + 0.031405749086161582 * tmp_44 + tmp_45 );
      real_t tmp_365 = tmp_364 * tmp_41;
      real_t tmp_366 = tmp_364 * tmp_48;
      real_t tmp_367 = tmp_364 * tmp_50;
      real_t tmp_368 = -tmp_357 - tmp_358 - tmp_360 - tmp_361 - tmp_362 - tmp_363 - tmp_365 - tmp_366 - tmp_367 + 1;
      real_t tmp_369 = tmp_362 + tmp_363 + tmp_367;
      real_t tmp_370 = tmp_358 + tmp_361 + tmp_366;
      real_t tmp_371 = tmp_357 + tmp_360 + tmp_365;
      real_t tmp_372 =
          tmp_64 * ( tmp_11 * ( tmp_369 - 1.0 / 4.0 ) + tmp_14 * ( tmp_370 - 1.0 / 4.0 ) + tmp_8 * ( tmp_371 - 1.0 / 4.0 ) );
      real_t tmp_373 = 0.0068572537431980923 * tmp_66;
      real_t tmp_374 = tmp_19 * ( 0.19601935860219369 * tmp_21 + 0.19601935860219369 * tmp_22 + tmp_23 );
      real_t tmp_375 = tmp_374 * tmp_6;
      real_t tmp_376 = tmp_26 * tmp_374;
      real_t tmp_377 = tmp_19 * ( 0.19601935860219369 * tmp_30 + 0.19601935860219369 * tmp_31 + tmp_32 );
      real_t tmp_378 = tmp_28 * tmp_377;
      real_t tmp_379 = tmp_35 * tmp_377;
      real_t tmp_380 = tmp_37 * tmp_374;
      real_t tmp_381 = tmp_377 * tmp_39;
      real_t tmp_382 = tmp_19 * ( 0.19601935860219369 * tmp_43 + 0.19601935860219369 * tmp_44 + tmp_45 );
      real_t tmp_383 = tmp_382 * tmp_41;
      real_t tmp_384 = tmp_382 * tmp_48;
      real_t tmp_385 = tmp_382 * tmp_50;
      real_t tmp_386 = -tmp_375 - tmp_376 - tmp_378 - tmp_379 - tmp_380 - tmp_381 - tmp_383 - tmp_384 - tmp_385 + 1;
      real_t tmp_387 = tmp_380 + tmp_381 + tmp_385;
      real_t tmp_388 = tmp_376 + tmp_379 + tmp_384;
      real_t tmp_389 = tmp_375 + tmp_378 + tmp_383;
      real_t tmp_390 =
          tmp_64 * ( tmp_11 * ( tmp_387 - 1.0 / 4.0 ) + tmp_14 * ( tmp_388 - 1.0 / 4.0 ) + tmp_8 * ( tmp_389 - 1.0 / 4.0 ) );
      real_t tmp_391 = 0.037198804536718075 * tmp_66;
      real_t tmp_392 = tmp_19 * ( 0.40446199974765351 * tmp_21 + 0.40446199974765351 * tmp_22 + tmp_23 );
      real_t tmp_393 = tmp_392 * tmp_6;
      real_t tmp_394 = tmp_26 * tmp_392;
      real_t tmp_395 = tmp_19 * ( 0.40446199974765351 * tmp_30 + 0.40446199974765351 * tmp_31 + tmp_32 );
      real_t tmp_396 = tmp_28 * tmp_395;
      real_t tmp_397 = tmp_35 * tmp_395;
      real_t tmp_398 = tmp_37 * tmp_392;
      real_t tmp_399 = tmp_39 * tmp_395;
      real_t tmp_400 = tmp_19 * ( 0.40446199974765351 * tmp_43 + 0.40446199974765351 * tmp_44 + tmp_45 );
      real_t tmp_401 = tmp_400 * tmp_41;
      real_t tmp_402 = tmp_400 * tmp_48;
      real_t tmp_403 = tmp_400 * tmp_50;
      real_t tmp_404 = -tmp_393 - tmp_394 - tmp_396 - tmp_397 - tmp_398 - tmp_399 - tmp_401 - tmp_402 - tmp_403 + 1;
      real_t tmp_405 = tmp_398 + tmp_399 + tmp_403;
      real_t tmp_406 = tmp_394 + tmp_397 + tmp_402;
      real_t tmp_407 = tmp_393 + tmp_396 + tmp_401;
      real_t tmp_408 =
          tmp_64 * ( tmp_11 * ( tmp_405 - 1.0 / 4.0 ) + tmp_14 * ( tmp_406 - 1.0 / 4.0 ) + tmp_8 * ( tmp_407 - 1.0 / 4.0 ) );
      real_t tmp_409 = 0.042507265838595799 * tmp_66;
      real_t tmp_410 = tmp_19 * ( 0.1711304259088916 * tmp_21 + 0.041227165399737475 * tmp_22 + tmp_23 );
      real_t tmp_411 = tmp_410 * tmp_6;
      real_t tmp_412 = tmp_26 * tmp_410;
      real_t tmp_413 = tmp_19 * ( 0.1711304259088916 * tmp_30 + 0.041227165399737475 * tmp_31 + tmp_32 );
      real_t tmp_414 = tmp_28 * tmp_413;
      real_t tmp_415 = tmp_35 * tmp_413;
      real_t tmp_416 = tmp_37 * tmp_410;
      real_t tmp_417 = tmp_39 * tmp_413;
      real_t tmp_418 = tmp_19 * ( 0.1711304259088916 * tmp_43 + 0.041227165399737475 * tmp_44 + tmp_45 );
      real_t tmp_419 = tmp_41 * tmp_418;
      real_t tmp_420 = tmp_418 * tmp_48;
      real_t tmp_421 = tmp_418 * tmp_50;
      real_t tmp_422 = -tmp_411 - tmp_412 - tmp_414 - tmp_415 - tmp_416 - tmp_417 - tmp_419 - tmp_420 - tmp_421 + 1;
      real_t tmp_423 = tmp_416 + tmp_417 + tmp_421;
      real_t tmp_424 = tmp_412 + tmp_415 + tmp_420;
      real_t tmp_425 = tmp_411 + tmp_414 + tmp_419;
      real_t tmp_426 =
          tmp_64 * ( tmp_11 * ( tmp_423 - 1.0 / 4.0 ) + tmp_14 * ( tmp_424 - 1.0 / 4.0 ) + tmp_8 * ( tmp_425 - 1.0 / 4.0 ) );
      real_t tmp_427 = 0.019202922745021479 * tmp_66;
      real_t a_0_0   = tmp_103 * ( tmp_102 * tmp_98 - tmp_56 * tmp_98 ) + tmp_121 * ( tmp_116 * tmp_120 - tmp_116 * tmp_56 ) +
                     tmp_139 * ( tmp_134 * tmp_138 - tmp_134 * tmp_56 ) + tmp_157 * ( tmp_152 * tmp_156 - tmp_152 * tmp_56 ) +
                     tmp_175 * ( tmp_170 * tmp_174 - tmp_170 * tmp_56 ) + tmp_193 * ( tmp_188 * tmp_192 - tmp_188 * tmp_56 ) +
                     tmp_211 * ( tmp_206 * tmp_210 - tmp_206 * tmp_56 ) + tmp_229 * ( tmp_224 * tmp_228 - tmp_224 * tmp_56 ) +
                     tmp_247 * ( tmp_242 * tmp_246 - tmp_242 * tmp_56 ) + tmp_265 * ( tmp_260 * tmp_264 - tmp_260 * tmp_56 ) +
                     tmp_283 * ( tmp_278 * tmp_282 - tmp_278 * tmp_56 ) + tmp_301 * ( tmp_296 * tmp_300 - tmp_296 * tmp_56 ) +
                     tmp_319 * ( tmp_314 * tmp_318 - tmp_314 * tmp_56 ) + tmp_337 * ( tmp_332 * tmp_336 - tmp_332 * tmp_56 ) +
                     tmp_355 * ( tmp_350 * tmp_354 - tmp_350 * tmp_56 ) + tmp_373 * ( tmp_368 * tmp_372 - tmp_368 * tmp_56 ) +
                     tmp_391 * ( tmp_386 * tmp_390 - tmp_386 * tmp_56 ) + tmp_409 * ( tmp_404 * tmp_408 - tmp_404 * tmp_56 ) +
                     tmp_427 * ( tmp_422 * tmp_426 - tmp_422 * tmp_56 ) + tmp_67 * ( -tmp_52 * tmp_56 + tmp_52 * tmp_65 ) +
                     tmp_85 * ( -tmp_56 * tmp_80 + tmp_80 * tmp_84 );
      real_t a_1_0 = tmp_103 * ( tmp_102 * tmp_99 - tmp_56 * tmp_99 ) + tmp_121 * ( tmp_117 * tmp_120 - tmp_117 * tmp_56 ) +
                     tmp_139 * ( tmp_135 * tmp_138 - tmp_135 * tmp_56 ) + tmp_157 * ( tmp_153 * tmp_156 - tmp_153 * tmp_56 ) +
                     tmp_175 * ( tmp_171 * tmp_174 - tmp_171 * tmp_56 ) + tmp_193 * ( tmp_189 * tmp_192 - tmp_189 * tmp_56 ) +
                     tmp_211 * ( tmp_207 * tmp_210 - tmp_207 * tmp_56 ) + tmp_229 * ( tmp_225 * tmp_228 - tmp_225 * tmp_56 ) +
                     tmp_247 * ( tmp_243 * tmp_246 - tmp_243 * tmp_56 ) + tmp_265 * ( tmp_261 * tmp_264 - tmp_261 * tmp_56 ) +
                     tmp_283 * ( tmp_279 * tmp_282 - tmp_279 * tmp_56 ) + tmp_301 * ( tmp_297 * tmp_300 - tmp_297 * tmp_56 ) +
                     tmp_319 * ( tmp_315 * tmp_318 - tmp_315 * tmp_56 ) + tmp_337 * ( tmp_333 * tmp_336 - tmp_333 * tmp_56 ) +
                     tmp_355 * ( tmp_351 * tmp_354 - tmp_351 * tmp_56 ) + tmp_373 * ( tmp_369 * tmp_372 - tmp_369 * tmp_56 ) +
                     tmp_391 * ( tmp_387 * tmp_390 - tmp_387 * tmp_56 ) + tmp_409 * ( tmp_405 * tmp_408 - tmp_405 * tmp_56 ) +
                     tmp_427 * ( tmp_423 * tmp_426 - tmp_423 * tmp_56 ) + tmp_67 * ( -tmp_56 * tmp_57 + tmp_57 * tmp_65 ) +
                     tmp_85 * ( -tmp_56 * tmp_81 + tmp_81 * tmp_84 );
      real_t a_2_0 = tmp_103 * ( tmp_100 * tmp_102 - tmp_100 * tmp_56 ) + tmp_121 * ( tmp_118 * tmp_120 - tmp_118 * tmp_56 ) +
                     tmp_139 * ( tmp_136 * tmp_138 - tmp_136 * tmp_56 ) + tmp_157 * ( tmp_154 * tmp_156 - tmp_154 * tmp_56 ) +
                     tmp_175 * ( tmp_172 * tmp_174 - tmp_172 * tmp_56 ) + tmp_193 * ( tmp_190 * tmp_192 - tmp_190 * tmp_56 ) +
                     tmp_211 * ( tmp_208 * tmp_210 - tmp_208 * tmp_56 ) + tmp_229 * ( tmp_226 * tmp_228 - tmp_226 * tmp_56 ) +
                     tmp_247 * ( tmp_244 * tmp_246 - tmp_244 * tmp_56 ) + tmp_265 * ( tmp_262 * tmp_264 - tmp_262 * tmp_56 ) +
                     tmp_283 * ( tmp_280 * tmp_282 - tmp_280 * tmp_56 ) + tmp_301 * ( tmp_298 * tmp_300 - tmp_298 * tmp_56 ) +
                     tmp_319 * ( tmp_316 * tmp_318 - tmp_316 * tmp_56 ) + tmp_337 * ( tmp_334 * tmp_336 - tmp_334 * tmp_56 ) +
                     tmp_355 * ( tmp_352 * tmp_354 - tmp_352 * tmp_56 ) + tmp_373 * ( tmp_370 * tmp_372 - tmp_370 * tmp_56 ) +
                     tmp_391 * ( tmp_388 * tmp_390 - tmp_388 * tmp_56 ) + tmp_409 * ( tmp_406 * tmp_408 - tmp_406 * tmp_56 ) +
                     tmp_427 * ( tmp_424 * tmp_426 - tmp_424 * tmp_56 ) + tmp_67 * ( -tmp_56 * tmp_58 + tmp_58 * tmp_65 ) +
                     tmp_85 * ( -tmp_56 * tmp_82 + tmp_82 * tmp_84 );
      real_t a_3_0 = tmp_103 * ( tmp_101 * tmp_102 - tmp_101 * tmp_56 ) + tmp_121 * ( tmp_119 * tmp_120 - tmp_119 * tmp_56 ) +
                     tmp_139 * ( tmp_137 * tmp_138 - tmp_137 * tmp_56 ) + tmp_157 * ( tmp_155 * tmp_156 - tmp_155 * tmp_56 ) +
                     tmp_175 * ( tmp_173 * tmp_174 - tmp_173 * tmp_56 ) + tmp_193 * ( tmp_191 * tmp_192 - tmp_191 * tmp_56 ) +
                     tmp_211 * ( tmp_209 * tmp_210 - tmp_209 * tmp_56 ) + tmp_229 * ( tmp_227 * tmp_228 - tmp_227 * tmp_56 ) +
                     tmp_247 * ( tmp_245 * tmp_246 - tmp_245 * tmp_56 ) + tmp_265 * ( tmp_263 * tmp_264 - tmp_263 * tmp_56 ) +
                     tmp_283 * ( tmp_281 * tmp_282 - tmp_281 * tmp_56 ) + tmp_301 * ( tmp_299 * tmp_300 - tmp_299 * tmp_56 ) +
                     tmp_319 * ( tmp_317 * tmp_318 - tmp_317 * tmp_56 ) + tmp_337 * ( tmp_335 * tmp_336 - tmp_335 * tmp_56 ) +
                     tmp_355 * ( tmp_353 * tmp_354 - tmp_353 * tmp_56 ) + tmp_373 * ( tmp_371 * tmp_372 - tmp_371 * tmp_56 ) +
                     tmp_391 * ( tmp_389 * tmp_390 - tmp_389 * tmp_56 ) + tmp_409 * ( tmp_407 * tmp_408 - tmp_407 * tmp_56 ) +
                     tmp_427 * ( tmp_425 * tmp_426 - tmp_425 * tmp_56 ) + tmp_67 * ( -tmp_56 * tmp_59 + tmp_59 * tmp_65 ) +
                     tmp_85 * ( -tmp_56 * tmp_83 + tmp_83 * tmp_84 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
      elMat( 3, 0 ) = a_3_0;
   }

   void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                  const std::vector< Point3D >& coordsElementOuter,
                                  const std::vector< Point3D >& coordsFacet,
                                  const Point3D&,
                                  const Point3D&,
                                  const Point3D&     outwardNormal,
                                  const DGBasisInfo& trialBasis,
                                  const DGBasisInfo& testBasis,
                                  int                trialDegree,
                                  int                testDegree,
                                  MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_0;
      real_t tmp_6  = p_affine_1_1 + tmp_2;
      real_t tmp_7  = tmp_5 * tmp_6;
      real_t tmp_8  = tmp_4 - tmp_7;
      real_t tmp_9  = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = p_affine_1_2 + tmp_9;
      real_t tmp_13 = tmp_12 * tmp_5;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = p_affine_2_2 + tmp_9;
      real_t tmp_16 = tmp_15 * tmp_6;
      real_t tmp_17 = tmp_1 * tmp_11;
      real_t tmp_18 = tmp_12 * tmp_14;
      real_t tmp_19 =
          1.0 / ( tmp_10 * tmp_4 - tmp_10 * tmp_7 + tmp_11 * tmp_13 + tmp_14 * tmp_16 - tmp_15 * tmp_17 - tmp_18 * tmp_3 );
      real_t tmp_20 = p_affine_8_2 + tmp_9;
      real_t tmp_21 = -p_affine_8_2;
      real_t tmp_22 = p_affine_9_2 + tmp_21;
      real_t tmp_23 = p_affine_10_2 + tmp_21;
      real_t tmp_24 = 0.031405749086161582 * tmp_22 + 0.93718850182767688 * tmp_23;
      real_t tmp_25 = tmp_19 * ( tmp_20 + tmp_24 );
      real_t tmp_26 = tmp_25 * tmp_8;
      real_t tmp_27 = tmp_14 * tmp_6 - tmp_17;
      real_t tmp_28 = tmp_25 * tmp_27;
      real_t tmp_29 = -tmp_1 * tmp_15 + tmp_13;
      real_t tmp_30 = p_affine_8_1 + tmp_2;
      real_t tmp_31 = -p_affine_8_1;
      real_t tmp_32 = p_affine_9_1 + tmp_31;
      real_t tmp_33 = p_affine_10_1 + tmp_31;
      real_t tmp_34 = 0.031405749086161582 * tmp_32 + 0.93718850182767688 * tmp_33;
      real_t tmp_35 = tmp_19 * ( tmp_30 + tmp_34 );
      real_t tmp_36 = tmp_29 * tmp_35;
      real_t tmp_37 = tmp_1 * tmp_10 - tmp_18;
      real_t tmp_38 = tmp_35 * tmp_37;
      real_t tmp_39 = tmp_11 * tmp_5 - tmp_14 * tmp_3;
      real_t tmp_40 = tmp_25 * tmp_39;
      real_t tmp_41 = -tmp_10 * tmp_5 + tmp_14 * tmp_15;
      real_t tmp_42 = tmp_35 * tmp_41;
      real_t tmp_43 = -tmp_12 * tmp_3 + tmp_16;
      real_t tmp_44 = p_affine_8_0 + tmp_0;
      real_t tmp_45 = -p_affine_8_0;
      real_t tmp_46 = p_affine_9_0 + tmp_45;
      real_t tmp_47 = p_affine_10_0 + tmp_45;
      real_t tmp_48 = 0.031405749086161582 * tmp_46 + 0.93718850182767688 * tmp_47;
      real_t tmp_49 = tmp_19 * ( tmp_44 + tmp_48 );
      real_t tmp_50 = tmp_43 * tmp_49;
      real_t tmp_51 = -tmp_10 * tmp_6 + tmp_11 * tmp_12;
      real_t tmp_52 = tmp_49 * tmp_51;
      real_t tmp_53 = tmp_10 * tmp_3 - tmp_11 * tmp_15;
      real_t tmp_54 = tmp_49 * tmp_53;
      real_t tmp_55 = -tmp_26 - tmp_28 - tmp_36 - tmp_38 - tmp_40 - tmp_42 - tmp_50 - tmp_52 - tmp_54 + 1;
      real_t tmp_56 = -p_affine_4_1;
      real_t tmp_57 = p_affine_6_1 + tmp_56;
      real_t tmp_58 = -p_affine_4_2;
      real_t tmp_59 = p_affine_7_2 + tmp_58;
      real_t tmp_60 = tmp_57 * tmp_59;
      real_t tmp_61 = p_affine_7_1 + tmp_56;
      real_t tmp_62 = p_affine_6_2 + tmp_58;
      real_t tmp_63 = tmp_61 * tmp_62;
      real_t tmp_64 = tmp_60 - tmp_63;
      real_t tmp_65 = p_affine_5_2 + tmp_58;
      real_t tmp_66 = -p_affine_4_0;
      real_t tmp_67 = p_affine_5_0 + tmp_66;
      real_t tmp_68 = p_affine_6_0 + tmp_66;
      real_t tmp_69 = tmp_61 * tmp_65;
      real_t tmp_70 = p_affine_7_0 + tmp_66;
      real_t tmp_71 = p_affine_5_1 + tmp_56;
      real_t tmp_72 = tmp_62 * tmp_71;
      real_t tmp_73 = tmp_59 * tmp_71;
      real_t tmp_74 = tmp_57 * tmp_65;
      real_t tmp_75 =
          1.0 / ( tmp_60 * tmp_67 - tmp_63 * tmp_67 + tmp_68 * tmp_69 - tmp_68 * tmp_73 + tmp_70 * tmp_72 - tmp_70 * tmp_74 );
      real_t tmp_76 = tmp_65 * tmp_75;
      real_t tmp_77 = tmp_69 - tmp_73;
      real_t tmp_78 = tmp_62 * tmp_75;
      real_t tmp_79 = tmp_72 - tmp_74;
      real_t tmp_80 = tmp_59 * tmp_75;
      real_t tmp_81 = -tmp_59 * tmp_68 + tmp_62 * tmp_70;
      real_t tmp_82 = tmp_59 * tmp_67 - tmp_65 * tmp_70;
      real_t tmp_83 = -tmp_62 * tmp_67 + tmp_65 * tmp_68;
      real_t tmp_84 = -tmp_57 * tmp_70 + tmp_61 * tmp_68;
      real_t tmp_85 = -tmp_61 * tmp_67 + tmp_70 * tmp_71;
      real_t tmp_86 = tmp_57 * tmp_67 - tmp_68 * tmp_71;
      real_t tmp_87 = 0.5 * p_affine_13_0 * ( tmp_64 * tmp_76 + tmp_77 * tmp_78 + tmp_79 * tmp_80 ) +
                      0.5 * p_affine_13_1 * ( tmp_76 * tmp_81 + tmp_78 * tmp_82 + tmp_80 * tmp_83 ) +
                      0.5 * p_affine_13_2 * ( tmp_76 * tmp_84 + tmp_78 * tmp_85 + tmp_80 * tmp_86 );
      real_t tmp_88 = p_affine_8_2 + tmp_58;
      real_t tmp_89 = tmp_75 * ( tmp_24 + tmp_88 );
      real_t tmp_90 = p_affine_8_1 + tmp_56;
      real_t tmp_91 = tmp_75 * ( tmp_34 + tmp_90 );
      real_t tmp_92 = p_affine_8_0 + tmp_66;
      real_t tmp_93 = tmp_75 * ( tmp_48 + tmp_92 );
      real_t tmp_94 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_95 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_96 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_97 = ( std::abs( tmp_23 * tmp_94 - tmp_33 * tmp_96 ) * std::abs( tmp_23 * tmp_94 - tmp_33 * tmp_96 ) ) +
                      ( std::abs( tmp_23 * tmp_95 - tmp_47 * tmp_96 ) * std::abs( tmp_23 * tmp_95 - tmp_47 * tmp_96 ) ) +
                      ( std::abs( tmp_33 * tmp_95 - tmp_47 * tmp_94 ) * std::abs( tmp_33 * tmp_95 - tmp_47 * tmp_94 ) );
      real_t tmp_98  = 1.0 * std::pow( tmp_97, -0.25 );
      real_t tmp_99  = tmp_98 * ( tmp_59 * ( tmp_79 * tmp_93 + tmp_83 * tmp_91 + tmp_86 * tmp_89 - 1.0 / 4.0 ) +
                                 tmp_62 * ( tmp_77 * tmp_93 + tmp_82 * tmp_91 + tmp_85 * tmp_89 - 1.0 / 4.0 ) +
                                 tmp_65 * ( tmp_64 * tmp_93 + tmp_81 * tmp_91 + tmp_84 * tmp_89 - 1.0 / 4.0 ) );
      real_t tmp_100 = 1.0 * std::pow( tmp_97, 1.0 / 2.0 );
      real_t tmp_101 = 0.0068572537431980923 * tmp_100;
      real_t tmp_102 = 0.19601935860219369 * tmp_22 + 0.60796128279561268 * tmp_23;
      real_t tmp_103 = tmp_19 * ( tmp_102 + tmp_20 );
      real_t tmp_104 = tmp_103 * tmp_8;
      real_t tmp_105 = tmp_103 * tmp_27;
      real_t tmp_106 = 0.19601935860219369 * tmp_32 + 0.60796128279561268 * tmp_33;
      real_t tmp_107 = tmp_19 * ( tmp_106 + tmp_30 );
      real_t tmp_108 = tmp_107 * tmp_29;
      real_t tmp_109 = tmp_107 * tmp_37;
      real_t tmp_110 = tmp_103 * tmp_39;
      real_t tmp_111 = tmp_107 * tmp_41;
      real_t tmp_112 = 0.19601935860219369 * tmp_46 + 0.60796128279561268 * tmp_47;
      real_t tmp_113 = tmp_19 * ( tmp_112 + tmp_44 );
      real_t tmp_114 = tmp_113 * tmp_43;
      real_t tmp_115 = tmp_113 * tmp_51;
      real_t tmp_116 = tmp_113 * tmp_53;
      real_t tmp_117 = -tmp_104 - tmp_105 - tmp_108 - tmp_109 - tmp_110 - tmp_111 - tmp_114 - tmp_115 - tmp_116 + 1;
      real_t tmp_118 = tmp_75 * ( tmp_102 + tmp_88 );
      real_t tmp_119 = tmp_75 * ( tmp_106 + tmp_90 );
      real_t tmp_120 = tmp_75 * ( tmp_112 + tmp_92 );
      real_t tmp_121 = tmp_98 * ( tmp_59 * ( tmp_118 * tmp_86 + tmp_119 * tmp_83 + tmp_120 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_118 * tmp_85 + tmp_119 * tmp_82 + tmp_120 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_118 * tmp_84 + tmp_119 * tmp_81 + tmp_120 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_122 = 0.037198804536718075 * tmp_100;
      real_t tmp_123 = 0.37605877282253791 * tmp_22 + 0.039308471900058539 * tmp_23;
      real_t tmp_124 = tmp_19 * ( tmp_123 + tmp_20 );
      real_t tmp_125 = tmp_124 * tmp_8;
      real_t tmp_126 = tmp_124 * tmp_27;
      real_t tmp_127 = 0.37605877282253791 * tmp_32 + 0.039308471900058539 * tmp_33;
      real_t tmp_128 = tmp_19 * ( tmp_127 + tmp_30 );
      real_t tmp_129 = tmp_128 * tmp_29;
      real_t tmp_130 = tmp_128 * tmp_37;
      real_t tmp_131 = tmp_124 * tmp_39;
      real_t tmp_132 = tmp_128 * tmp_41;
      real_t tmp_133 = 0.37605877282253791 * tmp_46 + 0.039308471900058539 * tmp_47;
      real_t tmp_134 = tmp_19 * ( tmp_133 + tmp_44 );
      real_t tmp_135 = tmp_134 * tmp_43;
      real_t tmp_136 = tmp_134 * tmp_51;
      real_t tmp_137 = tmp_134 * tmp_53;
      real_t tmp_138 = -tmp_125 - tmp_126 - tmp_129 - tmp_130 - tmp_131 - tmp_132 - tmp_135 - tmp_136 - tmp_137 + 1;
      real_t tmp_139 = tmp_75 * ( tmp_123 + tmp_88 );
      real_t tmp_140 = tmp_75 * ( tmp_127 + tmp_90 );
      real_t tmp_141 = tmp_75 * ( tmp_133 + tmp_92 );
      real_t tmp_142 = tmp_98 * ( tmp_59 * ( tmp_139 * tmp_86 + tmp_140 * tmp_83 + tmp_141 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_139 * tmp_85 + tmp_140 * tmp_82 + tmp_141 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_139 * tmp_84 + tmp_140 * tmp_81 + tmp_141 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_143 = 0.020848748529055869 * tmp_100;
      real_t tmp_144 = 0.78764240869137092 * tmp_22 + 0.1711304259088916 * tmp_23;
      real_t tmp_145 = tmp_19 * ( tmp_144 + tmp_20 );
      real_t tmp_146 = tmp_145 * tmp_8;
      real_t tmp_147 = tmp_145 * tmp_27;
      real_t tmp_148 = 0.78764240869137092 * tmp_32 + 0.1711304259088916 * tmp_33;
      real_t tmp_149 = tmp_19 * ( tmp_148 + tmp_30 );
      real_t tmp_150 = tmp_149 * tmp_29;
      real_t tmp_151 = tmp_149 * tmp_37;
      real_t tmp_152 = tmp_145 * tmp_39;
      real_t tmp_153 = tmp_149 * tmp_41;
      real_t tmp_154 = 0.78764240869137092 * tmp_46 + 0.1711304259088916 * tmp_47;
      real_t tmp_155 = tmp_19 * ( tmp_154 + tmp_44 );
      real_t tmp_156 = tmp_155 * tmp_43;
      real_t tmp_157 = tmp_155 * tmp_51;
      real_t tmp_158 = tmp_155 * tmp_53;
      real_t tmp_159 = -tmp_146 - tmp_147 - tmp_150 - tmp_151 - tmp_152 - tmp_153 - tmp_156 - tmp_157 - tmp_158 + 1;
      real_t tmp_160 = tmp_75 * ( tmp_144 + tmp_88 );
      real_t tmp_161 = tmp_75 * ( tmp_148 + tmp_90 );
      real_t tmp_162 = tmp_75 * ( tmp_154 + tmp_92 );
      real_t tmp_163 = tmp_98 * ( tmp_59 * ( tmp_160 * tmp_86 + tmp_161 * tmp_83 + tmp_162 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_160 * tmp_85 + tmp_161 * tmp_82 + tmp_162 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_160 * tmp_84 + tmp_161 * tmp_81 + tmp_162 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_164 = 0.019202922745021479 * tmp_100;
      real_t tmp_165 = 0.58463275527740355 * tmp_22 + 0.37605877282253791 * tmp_23;
      real_t tmp_166 = tmp_19 * ( tmp_165 + tmp_20 );
      real_t tmp_167 = tmp_166 * tmp_8;
      real_t tmp_168 = tmp_166 * tmp_27;
      real_t tmp_169 = 0.58463275527740355 * tmp_32 + 0.37605877282253791 * tmp_33;
      real_t tmp_170 = tmp_19 * ( tmp_169 + tmp_30 );
      real_t tmp_171 = tmp_170 * tmp_29;
      real_t tmp_172 = tmp_170 * tmp_37;
      real_t tmp_173 = tmp_166 * tmp_39;
      real_t tmp_174 = tmp_170 * tmp_41;
      real_t tmp_175 = 0.58463275527740355 * tmp_46 + 0.37605877282253791 * tmp_47;
      real_t tmp_176 = tmp_19 * ( tmp_175 + tmp_44 );
      real_t tmp_177 = tmp_176 * tmp_43;
      real_t tmp_178 = tmp_176 * tmp_51;
      real_t tmp_179 = tmp_176 * tmp_53;
      real_t tmp_180 = -tmp_167 - tmp_168 - tmp_171 - tmp_172 - tmp_173 - tmp_174 - tmp_177 - tmp_178 - tmp_179 + 1;
      real_t tmp_181 = tmp_75 * ( tmp_165 + tmp_88 );
      real_t tmp_182 = tmp_75 * ( tmp_169 + tmp_90 );
      real_t tmp_183 = tmp_75 * ( tmp_175 + tmp_92 );
      real_t tmp_184 = tmp_98 * ( tmp_59 * ( tmp_181 * tmp_86 + tmp_182 * tmp_83 + tmp_183 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_181 * tmp_85 + tmp_182 * tmp_82 + tmp_183 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_181 * tmp_84 + tmp_182 * tmp_81 + tmp_183 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_185 = 0.020848748529055869 * tmp_100;
      real_t tmp_186 = 0.041227165399737475 * tmp_22 + 0.78764240869137092 * tmp_23;
      real_t tmp_187 = tmp_19 * ( tmp_186 + tmp_20 );
      real_t tmp_188 = tmp_187 * tmp_8;
      real_t tmp_189 = tmp_187 * tmp_27;
      real_t tmp_190 = 0.041227165399737475 * tmp_32 + 0.78764240869137092 * tmp_33;
      real_t tmp_191 = tmp_19 * ( tmp_190 + tmp_30 );
      real_t tmp_192 = tmp_191 * tmp_29;
      real_t tmp_193 = tmp_191 * tmp_37;
      real_t tmp_194 = tmp_187 * tmp_39;
      real_t tmp_195 = tmp_191 * tmp_41;
      real_t tmp_196 = 0.041227165399737475 * tmp_46 + 0.78764240869137092 * tmp_47;
      real_t tmp_197 = tmp_19 * ( tmp_196 + tmp_44 );
      real_t tmp_198 = tmp_197 * tmp_43;
      real_t tmp_199 = tmp_197 * tmp_51;
      real_t tmp_200 = tmp_197 * tmp_53;
      real_t tmp_201 = -tmp_188 - tmp_189 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_198 - tmp_199 - tmp_200 + 1;
      real_t tmp_202 = tmp_75 * ( tmp_186 + tmp_88 );
      real_t tmp_203 = tmp_75 * ( tmp_190 + tmp_90 );
      real_t tmp_204 = tmp_75 * ( tmp_196 + tmp_92 );
      real_t tmp_205 = tmp_98 * ( tmp_59 * ( tmp_202 * tmp_86 + tmp_203 * tmp_83 + tmp_204 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_202 * tmp_85 + tmp_203 * tmp_82 + tmp_204 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_202 * tmp_84 + tmp_203 * tmp_81 + tmp_204 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_206 = 0.019202922745021479 * tmp_100;
      real_t tmp_207 = 0.039308471900058539 * tmp_22 + 0.58463275527740355 * tmp_23;
      real_t tmp_208 = tmp_19 * ( tmp_20 + tmp_207 );
      real_t tmp_209 = tmp_208 * tmp_8;
      real_t tmp_210 = tmp_208 * tmp_27;
      real_t tmp_211 = 0.039308471900058539 * tmp_32 + 0.58463275527740355 * tmp_33;
      real_t tmp_212 = tmp_19 * ( tmp_211 + tmp_30 );
      real_t tmp_213 = tmp_212 * tmp_29;
      real_t tmp_214 = tmp_212 * tmp_37;
      real_t tmp_215 = tmp_208 * tmp_39;
      real_t tmp_216 = tmp_212 * tmp_41;
      real_t tmp_217 = 0.039308471900058539 * tmp_46 + 0.58463275527740355 * tmp_47;
      real_t tmp_218 = tmp_19 * ( tmp_217 + tmp_44 );
      real_t tmp_219 = tmp_218 * tmp_43;
      real_t tmp_220 = tmp_218 * tmp_51;
      real_t tmp_221 = tmp_218 * tmp_53;
      real_t tmp_222 = -tmp_209 - tmp_210 - tmp_213 - tmp_214 - tmp_215 - tmp_216 - tmp_219 - tmp_220 - tmp_221 + 1;
      real_t tmp_223 = tmp_75 * ( tmp_207 + tmp_88 );
      real_t tmp_224 = tmp_75 * ( tmp_211 + tmp_90 );
      real_t tmp_225 = tmp_75 * ( tmp_217 + tmp_92 );
      real_t tmp_226 = tmp_98 * ( tmp_59 * ( tmp_223 * tmp_86 + tmp_224 * tmp_83 + tmp_225 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_223 * tmp_85 + tmp_224 * tmp_82 + tmp_225 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_223 * tmp_84 + tmp_224 * tmp_81 + tmp_225 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_227 = 0.020848748529055869 * tmp_100;
      real_t tmp_228 = 0.78764240869137092 * tmp_22 + 0.041227165399737475 * tmp_23;
      real_t tmp_229 = tmp_19 * ( tmp_20 + tmp_228 );
      real_t tmp_230 = tmp_229 * tmp_8;
      real_t tmp_231 = tmp_229 * tmp_27;
      real_t tmp_232 = 0.78764240869137092 * tmp_32 + 0.041227165399737475 * tmp_33;
      real_t tmp_233 = tmp_19 * ( tmp_232 + tmp_30 );
      real_t tmp_234 = tmp_233 * tmp_29;
      real_t tmp_235 = tmp_233 * tmp_37;
      real_t tmp_236 = tmp_229 * tmp_39;
      real_t tmp_237 = tmp_233 * tmp_41;
      real_t tmp_238 = 0.78764240869137092 * tmp_46 + 0.041227165399737475 * tmp_47;
      real_t tmp_239 = tmp_19 * ( tmp_238 + tmp_44 );
      real_t tmp_240 = tmp_239 * tmp_43;
      real_t tmp_241 = tmp_239 * tmp_51;
      real_t tmp_242 = tmp_239 * tmp_53;
      real_t tmp_243 = -tmp_230 - tmp_231 - tmp_234 - tmp_235 - tmp_236 - tmp_237 - tmp_240 - tmp_241 - tmp_242 + 1;
      real_t tmp_244 = tmp_75 * ( tmp_228 + tmp_88 );
      real_t tmp_245 = tmp_75 * ( tmp_232 + tmp_90 );
      real_t tmp_246 = tmp_75 * ( tmp_238 + tmp_92 );
      real_t tmp_247 = tmp_98 * ( tmp_59 * ( tmp_244 * tmp_86 + tmp_245 * tmp_83 + tmp_246 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_244 * tmp_85 + tmp_245 * tmp_82 + tmp_246 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_244 * tmp_84 + tmp_245 * tmp_81 + tmp_246 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_248 = 0.019202922745021479 * tmp_100;
      real_t tmp_249 = 0.58463275527740355 * tmp_22 + 0.039308471900058539 * tmp_23;
      real_t tmp_250 = tmp_19 * ( tmp_20 + tmp_249 );
      real_t tmp_251 = tmp_250 * tmp_8;
      real_t tmp_252 = tmp_250 * tmp_27;
      real_t tmp_253 = 0.58463275527740355 * tmp_32 + 0.039308471900058539 * tmp_33;
      real_t tmp_254 = tmp_19 * ( tmp_253 + tmp_30 );
      real_t tmp_255 = tmp_254 * tmp_29;
      real_t tmp_256 = tmp_254 * tmp_37;
      real_t tmp_257 = tmp_250 * tmp_39;
      real_t tmp_258 = tmp_254 * tmp_41;
      real_t tmp_259 = 0.58463275527740355 * tmp_46 + 0.039308471900058539 * tmp_47;
      real_t tmp_260 = tmp_19 * ( tmp_259 + tmp_44 );
      real_t tmp_261 = tmp_260 * tmp_43;
      real_t tmp_262 = tmp_260 * tmp_51;
      real_t tmp_263 = tmp_260 * tmp_53;
      real_t tmp_264 = -tmp_251 - tmp_252 - tmp_255 - tmp_256 - tmp_257 - tmp_258 - tmp_261 - tmp_262 - tmp_263 + 1;
      real_t tmp_265 = tmp_75 * ( tmp_249 + tmp_88 );
      real_t tmp_266 = tmp_75 * ( tmp_253 + tmp_90 );
      real_t tmp_267 = tmp_75 * ( tmp_259 + tmp_92 );
      real_t tmp_268 = tmp_98 * ( tmp_59 * ( tmp_265 * tmp_86 + tmp_266 * tmp_83 + tmp_267 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_265 * tmp_85 + tmp_266 * tmp_82 + tmp_267 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_265 * tmp_84 + tmp_266 * tmp_81 + tmp_267 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_269 = 0.020848748529055869 * tmp_100;
      real_t tmp_270 = 0.1711304259088916 * tmp_22 + 0.78764240869137092 * tmp_23;
      real_t tmp_271 = tmp_19 * ( tmp_20 + tmp_270 );
      real_t tmp_272 = tmp_271 * tmp_8;
      real_t tmp_273 = tmp_27 * tmp_271;
      real_t tmp_274 = 0.1711304259088916 * tmp_32 + 0.78764240869137092 * tmp_33;
      real_t tmp_275 = tmp_19 * ( tmp_274 + tmp_30 );
      real_t tmp_276 = tmp_275 * tmp_29;
      real_t tmp_277 = tmp_275 * tmp_37;
      real_t tmp_278 = tmp_271 * tmp_39;
      real_t tmp_279 = tmp_275 * tmp_41;
      real_t tmp_280 = 0.1711304259088916 * tmp_46 + 0.78764240869137092 * tmp_47;
      real_t tmp_281 = tmp_19 * ( tmp_280 + tmp_44 );
      real_t tmp_282 = tmp_281 * tmp_43;
      real_t tmp_283 = tmp_281 * tmp_51;
      real_t tmp_284 = tmp_281 * tmp_53;
      real_t tmp_285 = -tmp_272 - tmp_273 - tmp_276 - tmp_277 - tmp_278 - tmp_279 - tmp_282 - tmp_283 - tmp_284 + 1;
      real_t tmp_286 = tmp_75 * ( tmp_270 + tmp_88 );
      real_t tmp_287 = tmp_75 * ( tmp_274 + tmp_90 );
      real_t tmp_288 = tmp_75 * ( tmp_280 + tmp_92 );
      real_t tmp_289 = tmp_98 * ( tmp_59 * ( tmp_286 * tmp_86 + tmp_287 * tmp_83 + tmp_288 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_286 * tmp_85 + tmp_287 * tmp_82 + tmp_288 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_286 * tmp_84 + tmp_287 * tmp_81 + tmp_288 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_290 = 0.019202922745021479 * tmp_100;
      real_t tmp_291 = 0.37605877282253791 * tmp_22 + 0.58463275527740355 * tmp_23;
      real_t tmp_292 = tmp_19 * ( tmp_20 + tmp_291 );
      real_t tmp_293 = tmp_292 * tmp_8;
      real_t tmp_294 = tmp_27 * tmp_292;
      real_t tmp_295 = 0.37605877282253791 * tmp_32 + 0.58463275527740355 * tmp_33;
      real_t tmp_296 = tmp_19 * ( tmp_295 + tmp_30 );
      real_t tmp_297 = tmp_29 * tmp_296;
      real_t tmp_298 = tmp_296 * tmp_37;
      real_t tmp_299 = tmp_292 * tmp_39;
      real_t tmp_300 = tmp_296 * tmp_41;
      real_t tmp_301 = 0.37605877282253791 * tmp_46 + 0.58463275527740355 * tmp_47;
      real_t tmp_302 = tmp_19 * ( tmp_301 + tmp_44 );
      real_t tmp_303 = tmp_302 * tmp_43;
      real_t tmp_304 = tmp_302 * tmp_51;
      real_t tmp_305 = tmp_302 * tmp_53;
      real_t tmp_306 = -tmp_293 - tmp_294 - tmp_297 - tmp_298 - tmp_299 - tmp_300 - tmp_303 - tmp_304 - tmp_305 + 1;
      real_t tmp_307 = tmp_75 * ( tmp_291 + tmp_88 );
      real_t tmp_308 = tmp_75 * ( tmp_295 + tmp_90 );
      real_t tmp_309 = tmp_75 * ( tmp_301 + tmp_92 );
      real_t tmp_310 = tmp_98 * ( tmp_59 * ( tmp_307 * tmp_86 + tmp_308 * tmp_83 + tmp_309 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_307 * tmp_85 + tmp_308 * tmp_82 + tmp_309 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_307 * tmp_84 + tmp_308 * tmp_81 + tmp_309 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_311 = 0.020848748529055869 * tmp_100;
      real_t tmp_312 = 0.041227165399737475 * tmp_22 + 0.1711304259088916 * tmp_23;
      real_t tmp_313 = tmp_19 * ( tmp_20 + tmp_312 );
      real_t tmp_314 = tmp_313 * tmp_8;
      real_t tmp_315 = tmp_27 * tmp_313;
      real_t tmp_316 = 0.041227165399737475 * tmp_32 + 0.1711304259088916 * tmp_33;
      real_t tmp_317 = tmp_19 * ( tmp_30 + tmp_316 );
      real_t tmp_318 = tmp_29 * tmp_317;
      real_t tmp_319 = tmp_317 * tmp_37;
      real_t tmp_320 = tmp_313 * tmp_39;
      real_t tmp_321 = tmp_317 * tmp_41;
      real_t tmp_322 = 0.041227165399737475 * tmp_46 + 0.1711304259088916 * tmp_47;
      real_t tmp_323 = tmp_19 * ( tmp_322 + tmp_44 );
      real_t tmp_324 = tmp_323 * tmp_43;
      real_t tmp_325 = tmp_323 * tmp_51;
      real_t tmp_326 = tmp_323 * tmp_53;
      real_t tmp_327 = -tmp_314 - tmp_315 - tmp_318 - tmp_319 - tmp_320 - tmp_321 - tmp_324 - tmp_325 - tmp_326 + 1;
      real_t tmp_328 = tmp_75 * ( tmp_312 + tmp_88 );
      real_t tmp_329 = tmp_75 * ( tmp_316 + tmp_90 );
      real_t tmp_330 = tmp_75 * ( tmp_322 + tmp_92 );
      real_t tmp_331 = tmp_98 * ( tmp_59 * ( tmp_328 * tmp_86 + tmp_329 * tmp_83 + tmp_330 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_328 * tmp_85 + tmp_329 * tmp_82 + tmp_330 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_328 * tmp_84 + tmp_329 * tmp_81 + tmp_330 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_332 = 0.019202922745021479 * tmp_100;
      real_t tmp_333 = 0.40446199974765351 * tmp_22 + 0.19107600050469298 * tmp_23;
      real_t tmp_334 = tmp_19 * ( tmp_20 + tmp_333 );
      real_t tmp_335 = tmp_334 * tmp_8;
      real_t tmp_336 = tmp_27 * tmp_334;
      real_t tmp_337 = 0.40446199974765351 * tmp_32 + 0.19107600050469298 * tmp_33;
      real_t tmp_338 = tmp_19 * ( tmp_30 + tmp_337 );
      real_t tmp_339 = tmp_29 * tmp_338;
      real_t tmp_340 = tmp_338 * tmp_37;
      real_t tmp_341 = tmp_334 * tmp_39;
      real_t tmp_342 = tmp_338 * tmp_41;
      real_t tmp_343 = 0.40446199974765351 * tmp_46 + 0.19107600050469298 * tmp_47;
      real_t tmp_344 = tmp_19 * ( tmp_343 + tmp_44 );
      real_t tmp_345 = tmp_344 * tmp_43;
      real_t tmp_346 = tmp_344 * tmp_51;
      real_t tmp_347 = tmp_344 * tmp_53;
      real_t tmp_348 = -tmp_335 - tmp_336 - tmp_339 - tmp_340 - tmp_341 - tmp_342 - tmp_345 - tmp_346 - tmp_347 + 1;
      real_t tmp_349 = tmp_75 * ( tmp_333 + tmp_88 );
      real_t tmp_350 = tmp_75 * ( tmp_337 + tmp_90 );
      real_t tmp_351 = tmp_75 * ( tmp_343 + tmp_92 );
      real_t tmp_352 = tmp_98 * ( tmp_59 * ( tmp_349 * tmp_86 + tmp_350 * tmp_83 + tmp_351 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_349 * tmp_85 + tmp_350 * tmp_82 + tmp_351 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_349 * tmp_84 + tmp_350 * tmp_81 + tmp_351 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_353 = 0.042507265838595799 * tmp_100;
      real_t tmp_354 = 0.039308471900058539 * tmp_22 + 0.37605877282253791 * tmp_23;
      real_t tmp_355 = tmp_19 * ( tmp_20 + tmp_354 );
      real_t tmp_356 = tmp_355 * tmp_8;
      real_t tmp_357 = tmp_27 * tmp_355;
      real_t tmp_358 = 0.039308471900058539 * tmp_32 + 0.37605877282253791 * tmp_33;
      real_t tmp_359 = tmp_19 * ( tmp_30 + tmp_358 );
      real_t tmp_360 = tmp_29 * tmp_359;
      real_t tmp_361 = tmp_359 * tmp_37;
      real_t tmp_362 = tmp_355 * tmp_39;
      real_t tmp_363 = tmp_359 * tmp_41;
      real_t tmp_364 = 0.039308471900058539 * tmp_46 + 0.37605877282253791 * tmp_47;
      real_t tmp_365 = tmp_19 * ( tmp_364 + tmp_44 );
      real_t tmp_366 = tmp_365 * tmp_43;
      real_t tmp_367 = tmp_365 * tmp_51;
      real_t tmp_368 = tmp_365 * tmp_53;
      real_t tmp_369 = -tmp_356 - tmp_357 - tmp_360 - tmp_361 - tmp_362 - tmp_363 - tmp_366 - tmp_367 - tmp_368 + 1;
      real_t tmp_370 = tmp_75 * ( tmp_354 + tmp_88 );
      real_t tmp_371 = tmp_75 * ( tmp_358 + tmp_90 );
      real_t tmp_372 = tmp_75 * ( tmp_364 + tmp_92 );
      real_t tmp_373 = tmp_98 * ( tmp_59 * ( tmp_370 * tmp_86 + tmp_371 * tmp_83 + tmp_372 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_370 * tmp_85 + tmp_371 * tmp_82 + tmp_372 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_370 * tmp_84 + tmp_371 * tmp_81 + tmp_372 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_374 = 0.020848748529055869 * tmp_100;
      real_t tmp_375 = 0.93718850182767688 * tmp_22 + 0.031405749086161582 * tmp_23;
      real_t tmp_376 = tmp_19 * ( tmp_20 + tmp_375 );
      real_t tmp_377 = tmp_376 * tmp_8;
      real_t tmp_378 = tmp_27 * tmp_376;
      real_t tmp_379 = 0.93718850182767688 * tmp_32 + 0.031405749086161582 * tmp_33;
      real_t tmp_380 = tmp_19 * ( tmp_30 + tmp_379 );
      real_t tmp_381 = tmp_29 * tmp_380;
      real_t tmp_382 = tmp_37 * tmp_380;
      real_t tmp_383 = tmp_376 * tmp_39;
      real_t tmp_384 = tmp_380 * tmp_41;
      real_t tmp_385 = 0.93718850182767688 * tmp_46 + 0.031405749086161582 * tmp_47;
      real_t tmp_386 = tmp_19 * ( tmp_385 + tmp_44 );
      real_t tmp_387 = tmp_386 * tmp_43;
      real_t tmp_388 = tmp_386 * tmp_51;
      real_t tmp_389 = tmp_386 * tmp_53;
      real_t tmp_390 = -tmp_377 - tmp_378 - tmp_381 - tmp_382 - tmp_383 - tmp_384 - tmp_387 - tmp_388 - tmp_389 + 1;
      real_t tmp_391 = tmp_75 * ( tmp_375 + tmp_88 );
      real_t tmp_392 = tmp_75 * ( tmp_379 + tmp_90 );
      real_t tmp_393 = tmp_75 * ( tmp_385 + tmp_92 );
      real_t tmp_394 = tmp_98 * ( tmp_59 * ( tmp_391 * tmp_86 + tmp_392 * tmp_83 + tmp_393 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_391 * tmp_85 + tmp_392 * tmp_82 + tmp_393 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_391 * tmp_84 + tmp_392 * tmp_81 + tmp_393 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_395 = 0.0068572537431980923 * tmp_100;
      real_t tmp_396 = 0.60796128279561268 * tmp_22 + 0.19601935860219369 * tmp_23;
      real_t tmp_397 = tmp_19 * ( tmp_20 + tmp_396 );
      real_t tmp_398 = tmp_397 * tmp_8;
      real_t tmp_399 = tmp_27 * tmp_397;
      real_t tmp_400 = 0.60796128279561268 * tmp_32 + 0.19601935860219369 * tmp_33;
      real_t tmp_401 = tmp_19 * ( tmp_30 + tmp_400 );
      real_t tmp_402 = tmp_29 * tmp_401;
      real_t tmp_403 = tmp_37 * tmp_401;
      real_t tmp_404 = tmp_39 * tmp_397;
      real_t tmp_405 = tmp_401 * tmp_41;
      real_t tmp_406 = 0.60796128279561268 * tmp_46 + 0.19601935860219369 * tmp_47;
      real_t tmp_407 = tmp_19 * ( tmp_406 + tmp_44 );
      real_t tmp_408 = tmp_407 * tmp_43;
      real_t tmp_409 = tmp_407 * tmp_51;
      real_t tmp_410 = tmp_407 * tmp_53;
      real_t tmp_411 = -tmp_398 - tmp_399 - tmp_402 - tmp_403 - tmp_404 - tmp_405 - tmp_408 - tmp_409 - tmp_410 + 1;
      real_t tmp_412 = tmp_75 * ( tmp_396 + tmp_88 );
      real_t tmp_413 = tmp_75 * ( tmp_400 + tmp_90 );
      real_t tmp_414 = tmp_75 * ( tmp_406 + tmp_92 );
      real_t tmp_415 = tmp_98 * ( tmp_59 * ( tmp_412 * tmp_86 + tmp_413 * tmp_83 + tmp_414 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_412 * tmp_85 + tmp_413 * tmp_82 + tmp_414 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_412 * tmp_84 + tmp_413 * tmp_81 + tmp_414 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_416 = 0.037198804536718075 * tmp_100;
      real_t tmp_417 = 0.19107600050469298 * tmp_22 + 0.40446199974765351 * tmp_23;
      real_t tmp_418 = tmp_19 * ( tmp_20 + tmp_417 );
      real_t tmp_419 = tmp_418 * tmp_8;
      real_t tmp_420 = tmp_27 * tmp_418;
      real_t tmp_421 = 0.19107600050469298 * tmp_32 + 0.40446199974765351 * tmp_33;
      real_t tmp_422 = tmp_19 * ( tmp_30 + tmp_421 );
      real_t tmp_423 = tmp_29 * tmp_422;
      real_t tmp_424 = tmp_37 * tmp_422;
      real_t tmp_425 = tmp_39 * tmp_418;
      real_t tmp_426 = tmp_41 * tmp_422;
      real_t tmp_427 = 0.19107600050469298 * tmp_46 + 0.40446199974765351 * tmp_47;
      real_t tmp_428 = tmp_19 * ( tmp_427 + tmp_44 );
      real_t tmp_429 = tmp_428 * tmp_43;
      real_t tmp_430 = tmp_428 * tmp_51;
      real_t tmp_431 = tmp_428 * tmp_53;
      real_t tmp_432 = -tmp_419 - tmp_420 - tmp_423 - tmp_424 - tmp_425 - tmp_426 - tmp_429 - tmp_430 - tmp_431 + 1;
      real_t tmp_433 = tmp_75 * ( tmp_417 + tmp_88 );
      real_t tmp_434 = tmp_75 * ( tmp_421 + tmp_90 );
      real_t tmp_435 = tmp_75 * ( tmp_427 + tmp_92 );
      real_t tmp_436 = tmp_98 * ( tmp_59 * ( tmp_433 * tmp_86 + tmp_434 * tmp_83 + tmp_435 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_433 * tmp_85 + tmp_434 * tmp_82 + tmp_435 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_433 * tmp_84 + tmp_434 * tmp_81 + tmp_435 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_437 = 0.042507265838595799 * tmp_100;
      real_t tmp_438 = 0.031405749086161582 * tmp_22 + 0.031405749086161582 * tmp_23;
      real_t tmp_439 = tmp_19 * ( tmp_20 + tmp_438 );
      real_t tmp_440 = tmp_439 * tmp_8;
      real_t tmp_441 = tmp_27 * tmp_439;
      real_t tmp_442 = 0.031405749086161582 * tmp_32 + 0.031405749086161582 * tmp_33;
      real_t tmp_443 = tmp_19 * ( tmp_30 + tmp_442 );
      real_t tmp_444 = tmp_29 * tmp_443;
      real_t tmp_445 = tmp_37 * tmp_443;
      real_t tmp_446 = tmp_39 * tmp_439;
      real_t tmp_447 = tmp_41 * tmp_443;
      real_t tmp_448 = 0.031405749086161582 * tmp_46 + 0.031405749086161582 * tmp_47;
      real_t tmp_449 = tmp_19 * ( tmp_44 + tmp_448 );
      real_t tmp_450 = tmp_43 * tmp_449;
      real_t tmp_451 = tmp_449 * tmp_51;
      real_t tmp_452 = tmp_449 * tmp_53;
      real_t tmp_453 = -tmp_440 - tmp_441 - tmp_444 - tmp_445 - tmp_446 - tmp_447 - tmp_450 - tmp_451 - tmp_452 + 1;
      real_t tmp_454 = tmp_75 * ( tmp_438 + tmp_88 );
      real_t tmp_455 = tmp_75 * ( tmp_442 + tmp_90 );
      real_t tmp_456 = tmp_75 * ( tmp_448 + tmp_92 );
      real_t tmp_457 = tmp_98 * ( tmp_59 * ( tmp_454 * tmp_86 + tmp_455 * tmp_83 + tmp_456 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_454 * tmp_85 + tmp_455 * tmp_82 + tmp_456 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_454 * tmp_84 + tmp_455 * tmp_81 + tmp_456 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_458 = 0.0068572537431980923 * tmp_100;
      real_t tmp_459 = 0.19601935860219369 * tmp_22 + 0.19601935860219369 * tmp_23;
      real_t tmp_460 = tmp_19 * ( tmp_20 + tmp_459 );
      real_t tmp_461 = tmp_460 * tmp_8;
      real_t tmp_462 = tmp_27 * tmp_460;
      real_t tmp_463 = 0.19601935860219369 * tmp_32 + 0.19601935860219369 * tmp_33;
      real_t tmp_464 = tmp_19 * ( tmp_30 + tmp_463 );
      real_t tmp_465 = tmp_29 * tmp_464;
      real_t tmp_466 = tmp_37 * tmp_464;
      real_t tmp_467 = tmp_39 * tmp_460;
      real_t tmp_468 = tmp_41 * tmp_464;
      real_t tmp_469 = 0.19601935860219369 * tmp_46 + 0.19601935860219369 * tmp_47;
      real_t tmp_470 = tmp_19 * ( tmp_44 + tmp_469 );
      real_t tmp_471 = tmp_43 * tmp_470;
      real_t tmp_472 = tmp_470 * tmp_51;
      real_t tmp_473 = tmp_470 * tmp_53;
      real_t tmp_474 = -tmp_461 - tmp_462 - tmp_465 - tmp_466 - tmp_467 - tmp_468 - tmp_471 - tmp_472 - tmp_473 + 1;
      real_t tmp_475 = tmp_75 * ( tmp_459 + tmp_88 );
      real_t tmp_476 = tmp_75 * ( tmp_463 + tmp_90 );
      real_t tmp_477 = tmp_75 * ( tmp_469 + tmp_92 );
      real_t tmp_478 = tmp_98 * ( tmp_59 * ( tmp_475 * tmp_86 + tmp_476 * tmp_83 + tmp_477 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_475 * tmp_85 + tmp_476 * tmp_82 + tmp_477 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_475 * tmp_84 + tmp_476 * tmp_81 + tmp_477 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_479 = 0.037198804536718075 * tmp_100;
      real_t tmp_480 = 0.40446199974765351 * tmp_22 + 0.40446199974765351 * tmp_23;
      real_t tmp_481 = tmp_19 * ( tmp_20 + tmp_480 );
      real_t tmp_482 = tmp_481 * tmp_8;
      real_t tmp_483 = tmp_27 * tmp_481;
      real_t tmp_484 = 0.40446199974765351 * tmp_32 + 0.40446199974765351 * tmp_33;
      real_t tmp_485 = tmp_19 * ( tmp_30 + tmp_484 );
      real_t tmp_486 = tmp_29 * tmp_485;
      real_t tmp_487 = tmp_37 * tmp_485;
      real_t tmp_488 = tmp_39 * tmp_481;
      real_t tmp_489 = tmp_41 * tmp_485;
      real_t tmp_490 = 0.40446199974765351 * tmp_46 + 0.40446199974765351 * tmp_47;
      real_t tmp_491 = tmp_19 * ( tmp_44 + tmp_490 );
      real_t tmp_492 = tmp_43 * tmp_491;
      real_t tmp_493 = tmp_491 * tmp_51;
      real_t tmp_494 = tmp_491 * tmp_53;
      real_t tmp_495 = -tmp_482 - tmp_483 - tmp_486 - tmp_487 - tmp_488 - tmp_489 - tmp_492 - tmp_493 - tmp_494 + 1;
      real_t tmp_496 = tmp_75 * ( tmp_480 + tmp_88 );
      real_t tmp_497 = tmp_75 * ( tmp_484 + tmp_90 );
      real_t tmp_498 = tmp_75 * ( tmp_490 + tmp_92 );
      real_t tmp_499 = tmp_98 * ( tmp_59 * ( tmp_496 * tmp_86 + tmp_497 * tmp_83 + tmp_498 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_496 * tmp_85 + tmp_497 * tmp_82 + tmp_498 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_496 * tmp_84 + tmp_497 * tmp_81 + tmp_498 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_500 = 0.042507265838595799 * tmp_100;
      real_t tmp_501 = 0.1711304259088916 * tmp_22 + 0.041227165399737475 * tmp_23;
      real_t tmp_502 = tmp_19 * ( tmp_20 + tmp_501 );
      real_t tmp_503 = tmp_502 * tmp_8;
      real_t tmp_504 = tmp_27 * tmp_502;
      real_t tmp_505 = 0.1711304259088916 * tmp_32 + 0.041227165399737475 * tmp_33;
      real_t tmp_506 = tmp_19 * ( tmp_30 + tmp_505 );
      real_t tmp_507 = tmp_29 * tmp_506;
      real_t tmp_508 = tmp_37 * tmp_506;
      real_t tmp_509 = tmp_39 * tmp_502;
      real_t tmp_510 = tmp_41 * tmp_506;
      real_t tmp_511 = 0.1711304259088916 * tmp_46 + 0.041227165399737475 * tmp_47;
      real_t tmp_512 = tmp_19 * ( tmp_44 + tmp_511 );
      real_t tmp_513 = tmp_43 * tmp_512;
      real_t tmp_514 = tmp_51 * tmp_512;
      real_t tmp_515 = tmp_512 * tmp_53;
      real_t tmp_516 = -tmp_503 - tmp_504 - tmp_507 - tmp_508 - tmp_509 - tmp_510 - tmp_513 - tmp_514 - tmp_515 + 1;
      real_t tmp_517 = tmp_75 * ( tmp_501 + tmp_88 );
      real_t tmp_518 = tmp_75 * ( tmp_505 + tmp_90 );
      real_t tmp_519 = tmp_75 * ( tmp_511 + tmp_92 );
      real_t tmp_520 = tmp_98 * ( tmp_59 * ( tmp_517 * tmp_86 + tmp_518 * tmp_83 + tmp_519 * tmp_79 - 1.0 / 4.0 ) +
                                  tmp_62 * ( tmp_517 * tmp_85 + tmp_518 * tmp_82 + tmp_519 * tmp_77 - 1.0 / 4.0 ) +
                                  tmp_65 * ( tmp_517 * tmp_84 + tmp_518 * tmp_81 + tmp_519 * tmp_64 - 1.0 / 4.0 ) );
      real_t tmp_521 = 0.019202922745021479 * tmp_100;
      real_t tmp_522 = tmp_40 + tmp_42 + tmp_54;
      real_t tmp_523 = tmp_110 + tmp_111 + tmp_116;
      real_t tmp_524 = tmp_131 + tmp_132 + tmp_137;
      real_t tmp_525 = tmp_152 + tmp_153 + tmp_158;
      real_t tmp_526 = tmp_173 + tmp_174 + tmp_179;
      real_t tmp_527 = tmp_194 + tmp_195 + tmp_200;
      real_t tmp_528 = tmp_215 + tmp_216 + tmp_221;
      real_t tmp_529 = tmp_236 + tmp_237 + tmp_242;
      real_t tmp_530 = tmp_257 + tmp_258 + tmp_263;
      real_t tmp_531 = tmp_278 + tmp_279 + tmp_284;
      real_t tmp_532 = tmp_299 + tmp_300 + tmp_305;
      real_t tmp_533 = tmp_320 + tmp_321 + tmp_326;
      real_t tmp_534 = tmp_341 + tmp_342 + tmp_347;
      real_t tmp_535 = tmp_362 + tmp_363 + tmp_368;
      real_t tmp_536 = tmp_383 + tmp_384 + tmp_389;
      real_t tmp_537 = tmp_404 + tmp_405 + tmp_410;
      real_t tmp_538 = tmp_425 + tmp_426 + tmp_431;
      real_t tmp_539 = tmp_446 + tmp_447 + tmp_452;
      real_t tmp_540 = tmp_467 + tmp_468 + tmp_473;
      real_t tmp_541 = tmp_488 + tmp_489 + tmp_494;
      real_t tmp_542 = tmp_509 + tmp_510 + tmp_515;
      real_t tmp_543 = tmp_28 + tmp_38 + tmp_52;
      real_t tmp_544 = tmp_105 + tmp_109 + tmp_115;
      real_t tmp_545 = tmp_126 + tmp_130 + tmp_136;
      real_t tmp_546 = tmp_147 + tmp_151 + tmp_157;
      real_t tmp_547 = tmp_168 + tmp_172 + tmp_178;
      real_t tmp_548 = tmp_189 + tmp_193 + tmp_199;
      real_t tmp_549 = tmp_210 + tmp_214 + tmp_220;
      real_t tmp_550 = tmp_231 + tmp_235 + tmp_241;
      real_t tmp_551 = tmp_252 + tmp_256 + tmp_262;
      real_t tmp_552 = tmp_273 + tmp_277 + tmp_283;
      real_t tmp_553 = tmp_294 + tmp_298 + tmp_304;
      real_t tmp_554 = tmp_315 + tmp_319 + tmp_325;
      real_t tmp_555 = tmp_336 + tmp_340 + tmp_346;
      real_t tmp_556 = tmp_357 + tmp_361 + tmp_367;
      real_t tmp_557 = tmp_378 + tmp_382 + tmp_388;
      real_t tmp_558 = tmp_399 + tmp_403 + tmp_409;
      real_t tmp_559 = tmp_420 + tmp_424 + tmp_430;
      real_t tmp_560 = tmp_441 + tmp_445 + tmp_451;
      real_t tmp_561 = tmp_462 + tmp_466 + tmp_472;
      real_t tmp_562 = tmp_483 + tmp_487 + tmp_493;
      real_t tmp_563 = tmp_504 + tmp_508 + tmp_514;
      real_t tmp_564 = tmp_26 + tmp_36 + tmp_50;
      real_t tmp_565 = tmp_104 + tmp_108 + tmp_114;
      real_t tmp_566 = tmp_125 + tmp_129 + tmp_135;
      real_t tmp_567 = tmp_146 + tmp_150 + tmp_156;
      real_t tmp_568 = tmp_167 + tmp_171 + tmp_177;
      real_t tmp_569 = tmp_188 + tmp_192 + tmp_198;
      real_t tmp_570 = tmp_209 + tmp_213 + tmp_219;
      real_t tmp_571 = tmp_230 + tmp_234 + tmp_240;
      real_t tmp_572 = tmp_251 + tmp_255 + tmp_261;
      real_t tmp_573 = tmp_272 + tmp_276 + tmp_282;
      real_t tmp_574 = tmp_293 + tmp_297 + tmp_303;
      real_t tmp_575 = tmp_314 + tmp_318 + tmp_324;
      real_t tmp_576 = tmp_335 + tmp_339 + tmp_345;
      real_t tmp_577 = tmp_356 + tmp_360 + tmp_366;
      real_t tmp_578 = tmp_377 + tmp_381 + tmp_387;
      real_t tmp_579 = tmp_398 + tmp_402 + tmp_408;
      real_t tmp_580 = tmp_419 + tmp_423 + tmp_429;
      real_t tmp_581 = tmp_440 + tmp_444 + tmp_450;
      real_t tmp_582 = tmp_461 + tmp_465 + tmp_471;
      real_t tmp_583 = tmp_482 + tmp_486 + tmp_492;
      real_t tmp_584 = tmp_503 + tmp_507 + tmp_513;
      real_t a_0_0   = tmp_101 * ( -tmp_55 * tmp_87 - tmp_55 * tmp_99 ) + tmp_122 * ( -tmp_117 * tmp_121 - tmp_117 * tmp_87 ) +
                     tmp_143 * ( -tmp_138 * tmp_142 - tmp_138 * tmp_87 ) + tmp_164 * ( -tmp_159 * tmp_163 - tmp_159 * tmp_87 ) +
                     tmp_185 * ( -tmp_180 * tmp_184 - tmp_180 * tmp_87 ) + tmp_206 * ( -tmp_201 * tmp_205 - tmp_201 * tmp_87 ) +
                     tmp_227 * ( -tmp_222 * tmp_226 - tmp_222 * tmp_87 ) + tmp_248 * ( -tmp_243 * tmp_247 - tmp_243 * tmp_87 ) +
                     tmp_269 * ( -tmp_264 * tmp_268 - tmp_264 * tmp_87 ) + tmp_290 * ( -tmp_285 * tmp_289 - tmp_285 * tmp_87 ) +
                     tmp_311 * ( -tmp_306 * tmp_310 - tmp_306 * tmp_87 ) + tmp_332 * ( -tmp_327 * tmp_331 - tmp_327 * tmp_87 ) +
                     tmp_353 * ( -tmp_348 * tmp_352 - tmp_348 * tmp_87 ) + tmp_374 * ( -tmp_369 * tmp_373 - tmp_369 * tmp_87 ) +
                     tmp_395 * ( -tmp_390 * tmp_394 - tmp_390 * tmp_87 ) + tmp_416 * ( -tmp_411 * tmp_415 - tmp_411 * tmp_87 ) +
                     tmp_437 * ( -tmp_432 * tmp_436 - tmp_432 * tmp_87 ) + tmp_458 * ( -tmp_453 * tmp_457 - tmp_453 * tmp_87 ) +
                     tmp_479 * ( -tmp_474 * tmp_478 - tmp_474 * tmp_87 ) + tmp_500 * ( -tmp_495 * tmp_499 - tmp_495 * tmp_87 ) +
                     tmp_521 * ( -tmp_516 * tmp_520 - tmp_516 * tmp_87 );
      real_t a_1_0 = tmp_101 * ( -tmp_522 * tmp_87 - tmp_522 * tmp_99 ) + tmp_122 * ( -tmp_121 * tmp_523 - tmp_523 * tmp_87 ) +
                     tmp_143 * ( -tmp_142 * tmp_524 - tmp_524 * tmp_87 ) + tmp_164 * ( -tmp_163 * tmp_525 - tmp_525 * tmp_87 ) +
                     tmp_185 * ( -tmp_184 * tmp_526 - tmp_526 * tmp_87 ) + tmp_206 * ( -tmp_205 * tmp_527 - tmp_527 * tmp_87 ) +
                     tmp_227 * ( -tmp_226 * tmp_528 - tmp_528 * tmp_87 ) + tmp_248 * ( -tmp_247 * tmp_529 - tmp_529 * tmp_87 ) +
                     tmp_269 * ( -tmp_268 * tmp_530 - tmp_530 * tmp_87 ) + tmp_290 * ( -tmp_289 * tmp_531 - tmp_531 * tmp_87 ) +
                     tmp_311 * ( -tmp_310 * tmp_532 - tmp_532 * tmp_87 ) + tmp_332 * ( -tmp_331 * tmp_533 - tmp_533 * tmp_87 ) +
                     tmp_353 * ( -tmp_352 * tmp_534 - tmp_534 * tmp_87 ) + tmp_374 * ( -tmp_373 * tmp_535 - tmp_535 * tmp_87 ) +
                     tmp_395 * ( -tmp_394 * tmp_536 - tmp_536 * tmp_87 ) + tmp_416 * ( -tmp_415 * tmp_537 - tmp_537 * tmp_87 ) +
                     tmp_437 * ( -tmp_436 * tmp_538 - tmp_538 * tmp_87 ) + tmp_458 * ( -tmp_457 * tmp_539 - tmp_539 * tmp_87 ) +
                     tmp_479 * ( -tmp_478 * tmp_540 - tmp_540 * tmp_87 ) + tmp_500 * ( -tmp_499 * tmp_541 - tmp_541 * tmp_87 ) +
                     tmp_521 * ( -tmp_520 * tmp_542 - tmp_542 * tmp_87 );
      real_t a_2_0 = tmp_101 * ( -tmp_543 * tmp_87 - tmp_543 * tmp_99 ) + tmp_122 * ( -tmp_121 * tmp_544 - tmp_544 * tmp_87 ) +
                     tmp_143 * ( -tmp_142 * tmp_545 - tmp_545 * tmp_87 ) + tmp_164 * ( -tmp_163 * tmp_546 - tmp_546 * tmp_87 ) +
                     tmp_185 * ( -tmp_184 * tmp_547 - tmp_547 * tmp_87 ) + tmp_206 * ( -tmp_205 * tmp_548 - tmp_548 * tmp_87 ) +
                     tmp_227 * ( -tmp_226 * tmp_549 - tmp_549 * tmp_87 ) + tmp_248 * ( -tmp_247 * tmp_550 - tmp_550 * tmp_87 ) +
                     tmp_269 * ( -tmp_268 * tmp_551 - tmp_551 * tmp_87 ) + tmp_290 * ( -tmp_289 * tmp_552 - tmp_552 * tmp_87 ) +
                     tmp_311 * ( -tmp_310 * tmp_553 - tmp_553 * tmp_87 ) + tmp_332 * ( -tmp_331 * tmp_554 - tmp_554 * tmp_87 ) +
                     tmp_353 * ( -tmp_352 * tmp_555 - tmp_555 * tmp_87 ) + tmp_374 * ( -tmp_373 * tmp_556 - tmp_556 * tmp_87 ) +
                     tmp_395 * ( -tmp_394 * tmp_557 - tmp_557 * tmp_87 ) + tmp_416 * ( -tmp_415 * tmp_558 - tmp_558 * tmp_87 ) +
                     tmp_437 * ( -tmp_436 * tmp_559 - tmp_559 * tmp_87 ) + tmp_458 * ( -tmp_457 * tmp_560 - tmp_560 * tmp_87 ) +
                     tmp_479 * ( -tmp_478 * tmp_561 - tmp_561 * tmp_87 ) + tmp_500 * ( -tmp_499 * tmp_562 - tmp_562 * tmp_87 ) +
                     tmp_521 * ( -tmp_520 * tmp_563 - tmp_563 * tmp_87 );
      real_t a_3_0 = tmp_101 * ( -tmp_564 * tmp_87 - tmp_564 * tmp_99 ) + tmp_122 * ( -tmp_121 * tmp_565 - tmp_565 * tmp_87 ) +
                     tmp_143 * ( -tmp_142 * tmp_566 - tmp_566 * tmp_87 ) + tmp_164 * ( -tmp_163 * tmp_567 - tmp_567 * tmp_87 ) +
                     tmp_185 * ( -tmp_184 * tmp_568 - tmp_568 * tmp_87 ) + tmp_206 * ( -tmp_205 * tmp_569 - tmp_569 * tmp_87 ) +
                     tmp_227 * ( -tmp_226 * tmp_570 - tmp_570 * tmp_87 ) + tmp_248 * ( -tmp_247 * tmp_571 - tmp_571 * tmp_87 ) +
                     tmp_269 * ( -tmp_268 * tmp_572 - tmp_572 * tmp_87 ) + tmp_290 * ( -tmp_289 * tmp_573 - tmp_573 * tmp_87 ) +
                     tmp_311 * ( -tmp_310 * tmp_574 - tmp_574 * tmp_87 ) + tmp_332 * ( -tmp_331 * tmp_575 - tmp_575 * tmp_87 ) +
                     tmp_353 * ( -tmp_352 * tmp_576 - tmp_576 * tmp_87 ) + tmp_374 * ( -tmp_373 * tmp_577 - tmp_577 * tmp_87 ) +
                     tmp_395 * ( -tmp_394 * tmp_578 - tmp_578 * tmp_87 ) + tmp_416 * ( -tmp_415 * tmp_579 - tmp_579 * tmp_87 ) +
                     tmp_437 * ( -tmp_436 * tmp_580 - tmp_580 * tmp_87 ) + tmp_458 * ( -tmp_457 * tmp_581 - tmp_581 * tmp_87 ) +
                     tmp_479 * ( -tmp_478 * tmp_582 - tmp_582 * tmp_87 ) + tmp_500 * ( -tmp_499 * tmp_583 - tmp_583 * tmp_87 ) +
                     tmp_521 * ( -tmp_520 * tmp_584 - tmp_584 * tmp_87 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
      elMat( 3, 0 ) = a_3_0;
   }

   void integrateFacetDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                           const std::vector< Point3D >& coordsFacet,
                                           const Point3D&,
                                           const Point3D&     outwardNormal,
                                           const DGBasisInfo& trialBasis,
                                           const DGBasisInfo& testBasis,
                                           int                trialDegree,
                                           int                testDegree,
                                           MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = p_affine_2_0 + tmp_0;
      real_t tmp_5  = p_affine_1_1 + tmp_2;
      real_t tmp_6  = tmp_1 * tmp_3 - tmp_4 * tmp_5;
      real_t tmp_7  = -p_affine_0_2;
      real_t tmp_8  = p_affine_3_2 + tmp_7;
      real_t tmp_9  = tmp_3 * tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_7;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11 * tmp_4;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_13 * tmp_14;
      real_t tmp_16 = tmp_11 * tmp_14;
      real_t tmp_17 = tmp_4 * tmp_8;
      real_t tmp_18 = tmp_13 * tmp_3;
      real_t tmp_19 =
          1.0 / ( -tmp_1 * tmp_16 + tmp_1 * tmp_9 + tmp_10 * tmp_12 - tmp_10 * tmp_18 + tmp_15 * tmp_5 - tmp_17 * tmp_5 );
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_7;
      real_t tmp_24 = tmp_19 * ( 0.031405749086161582 * tmp_21 + 0.93718850182767688 * tmp_22 + tmp_23 );
      real_t tmp_25 = tmp_24 * tmp_6;
      real_t tmp_26 = -tmp_1 * tmp_11 + tmp_13 * tmp_5;
      real_t tmp_27 = tmp_24 * tmp_26;
      real_t tmp_28 = -tmp_1 * tmp_14 + tmp_10 * tmp_4;
      real_t tmp_29 = -p_affine_8_1;
      real_t tmp_30 = p_affine_9_1 + tmp_29;
      real_t tmp_31 = p_affine_10_1 + tmp_29;
      real_t tmp_32 = p_affine_8_1 + tmp_2;
      real_t tmp_33 = tmp_19 * ( 0.031405749086161582 * tmp_30 + 0.93718850182767688 * tmp_31 + tmp_32 );
      real_t tmp_34 = tmp_28 * tmp_33;
      real_t tmp_35 = tmp_1 * tmp_8 - tmp_10 * tmp_13;
      real_t tmp_36 = tmp_33 * tmp_35;
      real_t tmp_37 = tmp_12 - tmp_18;
      real_t tmp_38 = tmp_24 * tmp_37;
      real_t tmp_39 = tmp_15 - tmp_17;
      real_t tmp_40 = tmp_33 * tmp_39;
      real_t tmp_41 = -tmp_10 * tmp_3 + tmp_14 * tmp_5;
      real_t tmp_42 = -p_affine_8_0;
      real_t tmp_43 = p_affine_9_0 + tmp_42;
      real_t tmp_44 = p_affine_10_0 + tmp_42;
      real_t tmp_45 = p_affine_8_0 + tmp_0;
      real_t tmp_46 = tmp_19 * ( 0.031405749086161582 * tmp_43 + 0.93718850182767688 * tmp_44 + tmp_45 );
      real_t tmp_47 = tmp_41 * tmp_46;
      real_t tmp_48 = tmp_10 * tmp_11 - tmp_5 * tmp_8;
      real_t tmp_49 = tmp_46 * tmp_48;
      real_t tmp_50 = -tmp_16 + tmp_9;
      real_t tmp_51 = tmp_46 * tmp_50;
      real_t tmp_52 = tmp_38 + tmp_40 + tmp_51;
      real_t tmp_53 = tmp_27 + tmp_36 + tmp_49;
      real_t tmp_54 = tmp_25 + tmp_34 + tmp_47;
      real_t tmp_55 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_56 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_57 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_58 =
          1.0 * std::pow( ( std::abs( tmp_22 * tmp_55 - tmp_31 * tmp_57 ) * std::abs( tmp_22 * tmp_55 - tmp_31 * tmp_57 ) ) +
                              ( std::abs( tmp_22 * tmp_56 - tmp_44 * tmp_57 ) * std::abs( tmp_22 * tmp_56 - tmp_44 * tmp_57 ) ) +
                              ( std::abs( tmp_31 * tmp_56 - tmp_44 * tmp_55 ) * std::abs( tmp_31 * tmp_56 - tmp_44 * tmp_55 ) ),
                          0.25 );
      real_t tmp_59 = 0.0068572537431980923 * tmp_58 *
                      ( tmp_10 * ( tmp_52 - 1.0 / 4.0 ) + tmp_14 * ( tmp_53 - 1.0 / 4.0 ) + tmp_8 * ( tmp_54 - 1.0 / 4.0 ) );
      real_t tmp_60 = tmp_19 * ( 0.19601935860219369 * tmp_21 + 0.60796128279561268 * tmp_22 + tmp_23 );
      real_t tmp_61 = tmp_6 * tmp_60;
      real_t tmp_62 = tmp_26 * tmp_60;
      real_t tmp_63 = tmp_19 * ( 0.19601935860219369 * tmp_30 + 0.60796128279561268 * tmp_31 + tmp_32 );
      real_t tmp_64 = tmp_28 * tmp_63;
      real_t tmp_65 = tmp_35 * tmp_63;
      real_t tmp_66 = tmp_37 * tmp_60;
      real_t tmp_67 = tmp_39 * tmp_63;
      real_t tmp_68 = tmp_19 * ( 0.19601935860219369 * tmp_43 + 0.60796128279561268 * tmp_44 + tmp_45 );
      real_t tmp_69 = tmp_41 * tmp_68;
      real_t tmp_70 = tmp_48 * tmp_68;
      real_t tmp_71 = tmp_50 * tmp_68;
      real_t tmp_72 = tmp_66 + tmp_67 + tmp_71;
      real_t tmp_73 = tmp_62 + tmp_65 + tmp_70;
      real_t tmp_74 = tmp_61 + tmp_64 + tmp_69;
      real_t tmp_75 = 0.037198804536718075 * tmp_58 *
                      ( tmp_10 * ( tmp_72 - 1.0 / 4.0 ) + tmp_14 * ( tmp_73 - 1.0 / 4.0 ) + tmp_8 * ( tmp_74 - 1.0 / 4.0 ) );
      real_t tmp_76 = tmp_19 * ( 0.37605877282253791 * tmp_21 + 0.039308471900058539 * tmp_22 + tmp_23 );
      real_t tmp_77 = tmp_6 * tmp_76;
      real_t tmp_78 = tmp_26 * tmp_76;
      real_t tmp_79 = tmp_19 * ( 0.37605877282253791 * tmp_30 + 0.039308471900058539 * tmp_31 + tmp_32 );
      real_t tmp_80 = tmp_28 * tmp_79;
      real_t tmp_81 = tmp_35 * tmp_79;
      real_t tmp_82 = tmp_37 * tmp_76;
      real_t tmp_83 = tmp_39 * tmp_79;
      real_t tmp_84 = tmp_19 * ( 0.37605877282253791 * tmp_43 + 0.039308471900058539 * tmp_44 + tmp_45 );
      real_t tmp_85 = tmp_41 * tmp_84;
      real_t tmp_86 = tmp_48 * tmp_84;
      real_t tmp_87 = tmp_50 * tmp_84;
      real_t tmp_88 = tmp_82 + tmp_83 + tmp_87;
      real_t tmp_89 = tmp_78 + tmp_81 + tmp_86;
      real_t tmp_90 = tmp_77 + tmp_80 + tmp_85;
      real_t tmp_91 = 0.020848748529055869 * tmp_58 *
                      ( tmp_10 * ( tmp_88 - 1.0 / 4.0 ) + tmp_14 * ( tmp_89 - 1.0 / 4.0 ) + tmp_8 * ( tmp_90 - 1.0 / 4.0 ) );
      real_t tmp_92  = tmp_19 * ( 0.78764240869137092 * tmp_21 + 0.1711304259088916 * tmp_22 + tmp_23 );
      real_t tmp_93  = tmp_6 * tmp_92;
      real_t tmp_94  = tmp_26 * tmp_92;
      real_t tmp_95  = tmp_19 * ( 0.78764240869137092 * tmp_30 + 0.1711304259088916 * tmp_31 + tmp_32 );
      real_t tmp_96  = tmp_28 * tmp_95;
      real_t tmp_97  = tmp_35 * tmp_95;
      real_t tmp_98  = tmp_37 * tmp_92;
      real_t tmp_99  = tmp_39 * tmp_95;
      real_t tmp_100 = tmp_19 * ( 0.78764240869137092 * tmp_43 + 0.1711304259088916 * tmp_44 + tmp_45 );
      real_t tmp_101 = tmp_100 * tmp_41;
      real_t tmp_102 = tmp_100 * tmp_48;
      real_t tmp_103 = tmp_100 * tmp_50;
      real_t tmp_104 = tmp_103 + tmp_98 + tmp_99;
      real_t tmp_105 = tmp_102 + tmp_94 + tmp_97;
      real_t tmp_106 = tmp_101 + tmp_93 + tmp_96;
      real_t tmp_107 = 0.019202922745021479 * tmp_58 *
                       ( tmp_10 * ( tmp_104 - 1.0 / 4.0 ) + tmp_14 * ( tmp_105 - 1.0 / 4.0 ) + tmp_8 * ( tmp_106 - 1.0 / 4.0 ) );
      real_t tmp_108 = tmp_19 * ( 0.58463275527740355 * tmp_21 + 0.37605877282253791 * tmp_22 + tmp_23 );
      real_t tmp_109 = tmp_108 * tmp_6;
      real_t tmp_110 = tmp_108 * tmp_26;
      real_t tmp_111 = tmp_19 * ( 0.58463275527740355 * tmp_30 + 0.37605877282253791 * tmp_31 + tmp_32 );
      real_t tmp_112 = tmp_111 * tmp_28;
      real_t tmp_113 = tmp_111 * tmp_35;
      real_t tmp_114 = tmp_108 * tmp_37;
      real_t tmp_115 = tmp_111 * tmp_39;
      real_t tmp_116 = tmp_19 * ( 0.58463275527740355 * tmp_43 + 0.37605877282253791 * tmp_44 + tmp_45 );
      real_t tmp_117 = tmp_116 * tmp_41;
      real_t tmp_118 = tmp_116 * tmp_48;
      real_t tmp_119 = tmp_116 * tmp_50;
      real_t tmp_120 = tmp_114 + tmp_115 + tmp_119;
      real_t tmp_121 = tmp_110 + tmp_113 + tmp_118;
      real_t tmp_122 = tmp_109 + tmp_112 + tmp_117;
      real_t tmp_123 = 0.020848748529055869 * tmp_58 *
                       ( tmp_10 * ( tmp_120 - 1.0 / 4.0 ) + tmp_14 * ( tmp_121 - 1.0 / 4.0 ) + tmp_8 * ( tmp_122 - 1.0 / 4.0 ) );
      real_t tmp_124 = tmp_19 * ( 0.041227165399737475 * tmp_21 + 0.78764240869137092 * tmp_22 + tmp_23 );
      real_t tmp_125 = tmp_124 * tmp_6;
      real_t tmp_126 = tmp_124 * tmp_26;
      real_t tmp_127 = tmp_19 * ( 0.041227165399737475 * tmp_30 + 0.78764240869137092 * tmp_31 + tmp_32 );
      real_t tmp_128 = tmp_127 * tmp_28;
      real_t tmp_129 = tmp_127 * tmp_35;
      real_t tmp_130 = tmp_124 * tmp_37;
      real_t tmp_131 = tmp_127 * tmp_39;
      real_t tmp_132 = tmp_19 * ( 0.041227165399737475 * tmp_43 + 0.78764240869137092 * tmp_44 + tmp_45 );
      real_t tmp_133 = tmp_132 * tmp_41;
      real_t tmp_134 = tmp_132 * tmp_48;
      real_t tmp_135 = tmp_132 * tmp_50;
      real_t tmp_136 = tmp_130 + tmp_131 + tmp_135;
      real_t tmp_137 = tmp_126 + tmp_129 + tmp_134;
      real_t tmp_138 = tmp_125 + tmp_128 + tmp_133;
      real_t tmp_139 = 0.019202922745021479 * tmp_58 *
                       ( tmp_10 * ( tmp_136 - 1.0 / 4.0 ) + tmp_14 * ( tmp_137 - 1.0 / 4.0 ) + tmp_8 * ( tmp_138 - 1.0 / 4.0 ) );
      real_t tmp_140 = tmp_19 * ( 0.039308471900058539 * tmp_21 + 0.58463275527740355 * tmp_22 + tmp_23 );
      real_t tmp_141 = tmp_140 * tmp_6;
      real_t tmp_142 = tmp_140 * tmp_26;
      real_t tmp_143 = tmp_19 * ( 0.039308471900058539 * tmp_30 + 0.58463275527740355 * tmp_31 + tmp_32 );
      real_t tmp_144 = tmp_143 * tmp_28;
      real_t tmp_145 = tmp_143 * tmp_35;
      real_t tmp_146 = tmp_140 * tmp_37;
      real_t tmp_147 = tmp_143 * tmp_39;
      real_t tmp_148 = tmp_19 * ( 0.039308471900058539 * tmp_43 + 0.58463275527740355 * tmp_44 + tmp_45 );
      real_t tmp_149 = tmp_148 * tmp_41;
      real_t tmp_150 = tmp_148 * tmp_48;
      real_t tmp_151 = tmp_148 * tmp_50;
      real_t tmp_152 = tmp_146 + tmp_147 + tmp_151;
      real_t tmp_153 = tmp_142 + tmp_145 + tmp_150;
      real_t tmp_154 = tmp_141 + tmp_144 + tmp_149;
      real_t tmp_155 = 0.020848748529055869 * tmp_58 *
                       ( tmp_10 * ( tmp_152 - 1.0 / 4.0 ) + tmp_14 * ( tmp_153 - 1.0 / 4.0 ) + tmp_8 * ( tmp_154 - 1.0 / 4.0 ) );
      real_t tmp_156 = tmp_19 * ( 0.78764240869137092 * tmp_21 + 0.041227165399737475 * tmp_22 + tmp_23 );
      real_t tmp_157 = tmp_156 * tmp_6;
      real_t tmp_158 = tmp_156 * tmp_26;
      real_t tmp_159 = tmp_19 * ( 0.78764240869137092 * tmp_30 + 0.041227165399737475 * tmp_31 + tmp_32 );
      real_t tmp_160 = tmp_159 * tmp_28;
      real_t tmp_161 = tmp_159 * tmp_35;
      real_t tmp_162 = tmp_156 * tmp_37;
      real_t tmp_163 = tmp_159 * tmp_39;
      real_t tmp_164 = tmp_19 * ( 0.78764240869137092 * tmp_43 + 0.041227165399737475 * tmp_44 + tmp_45 );
      real_t tmp_165 = tmp_164 * tmp_41;
      real_t tmp_166 = tmp_164 * tmp_48;
      real_t tmp_167 = tmp_164 * tmp_50;
      real_t tmp_168 = tmp_162 + tmp_163 + tmp_167;
      real_t tmp_169 = tmp_158 + tmp_161 + tmp_166;
      real_t tmp_170 = tmp_157 + tmp_160 + tmp_165;
      real_t tmp_171 = 0.019202922745021479 * tmp_58 *
                       ( tmp_10 * ( tmp_168 - 1.0 / 4.0 ) + tmp_14 * ( tmp_169 - 1.0 / 4.0 ) + tmp_8 * ( tmp_170 - 1.0 / 4.0 ) );
      real_t tmp_172 = tmp_19 * ( 0.58463275527740355 * tmp_21 + 0.039308471900058539 * tmp_22 + tmp_23 );
      real_t tmp_173 = tmp_172 * tmp_6;
      real_t tmp_174 = tmp_172 * tmp_26;
      real_t tmp_175 = tmp_19 * ( 0.58463275527740355 * tmp_30 + 0.039308471900058539 * tmp_31 + tmp_32 );
      real_t tmp_176 = tmp_175 * tmp_28;
      real_t tmp_177 = tmp_175 * tmp_35;
      real_t tmp_178 = tmp_172 * tmp_37;
      real_t tmp_179 = tmp_175 * tmp_39;
      real_t tmp_180 = tmp_19 * ( 0.58463275527740355 * tmp_43 + 0.039308471900058539 * tmp_44 + tmp_45 );
      real_t tmp_181 = tmp_180 * tmp_41;
      real_t tmp_182 = tmp_180 * tmp_48;
      real_t tmp_183 = tmp_180 * tmp_50;
      real_t tmp_184 = tmp_178 + tmp_179 + tmp_183;
      real_t tmp_185 = tmp_174 + tmp_177 + tmp_182;
      real_t tmp_186 = tmp_173 + tmp_176 + tmp_181;
      real_t tmp_187 = 0.020848748529055869 * tmp_58 *
                       ( tmp_10 * ( tmp_184 - 1.0 / 4.0 ) + tmp_14 * ( tmp_185 - 1.0 / 4.0 ) + tmp_8 * ( tmp_186 - 1.0 / 4.0 ) );
      real_t tmp_188 = tmp_19 * ( 0.1711304259088916 * tmp_21 + 0.78764240869137092 * tmp_22 + tmp_23 );
      real_t tmp_189 = tmp_188 * tmp_6;
      real_t tmp_190 = tmp_188 * tmp_26;
      real_t tmp_191 = tmp_19 * ( 0.1711304259088916 * tmp_30 + 0.78764240869137092 * tmp_31 + tmp_32 );
      real_t tmp_192 = tmp_191 * tmp_28;
      real_t tmp_193 = tmp_191 * tmp_35;
      real_t tmp_194 = tmp_188 * tmp_37;
      real_t tmp_195 = tmp_191 * tmp_39;
      real_t tmp_196 = tmp_19 * ( 0.1711304259088916 * tmp_43 + 0.78764240869137092 * tmp_44 + tmp_45 );
      real_t tmp_197 = tmp_196 * tmp_41;
      real_t tmp_198 = tmp_196 * tmp_48;
      real_t tmp_199 = tmp_196 * tmp_50;
      real_t tmp_200 = tmp_194 + tmp_195 + tmp_199;
      real_t tmp_201 = tmp_190 + tmp_193 + tmp_198;
      real_t tmp_202 = tmp_189 + tmp_192 + tmp_197;
      real_t tmp_203 = 0.019202922745021479 * tmp_58 *
                       ( tmp_10 * ( tmp_200 - 1.0 / 4.0 ) + tmp_14 * ( tmp_201 - 1.0 / 4.0 ) + tmp_8 * ( tmp_202 - 1.0 / 4.0 ) );
      real_t tmp_204 = tmp_19 * ( 0.37605877282253791 * tmp_21 + 0.58463275527740355 * tmp_22 + tmp_23 );
      real_t tmp_205 = tmp_204 * tmp_6;
      real_t tmp_206 = tmp_204 * tmp_26;
      real_t tmp_207 = tmp_19 * ( 0.37605877282253791 * tmp_30 + 0.58463275527740355 * tmp_31 + tmp_32 );
      real_t tmp_208 = tmp_207 * tmp_28;
      real_t tmp_209 = tmp_207 * tmp_35;
      real_t tmp_210 = tmp_204 * tmp_37;
      real_t tmp_211 = tmp_207 * tmp_39;
      real_t tmp_212 = tmp_19 * ( 0.37605877282253791 * tmp_43 + 0.58463275527740355 * tmp_44 + tmp_45 );
      real_t tmp_213 = tmp_212 * tmp_41;
      real_t tmp_214 = tmp_212 * tmp_48;
      real_t tmp_215 = tmp_212 * tmp_50;
      real_t tmp_216 = tmp_210 + tmp_211 + tmp_215;
      real_t tmp_217 = tmp_206 + tmp_209 + tmp_214;
      real_t tmp_218 = tmp_205 + tmp_208 + tmp_213;
      real_t tmp_219 = 0.020848748529055869 * tmp_58 *
                       ( tmp_10 * ( tmp_216 - 1.0 / 4.0 ) + tmp_14 * ( tmp_217 - 1.0 / 4.0 ) + tmp_8 * ( tmp_218 - 1.0 / 4.0 ) );
      real_t tmp_220 = tmp_19 * ( 0.041227165399737475 * tmp_21 + 0.1711304259088916 * tmp_22 + tmp_23 );
      real_t tmp_221 = tmp_220 * tmp_6;
      real_t tmp_222 = tmp_220 * tmp_26;
      real_t tmp_223 = tmp_19 * ( 0.041227165399737475 * tmp_30 + 0.1711304259088916 * tmp_31 + tmp_32 );
      real_t tmp_224 = tmp_223 * tmp_28;
      real_t tmp_225 = tmp_223 * tmp_35;
      real_t tmp_226 = tmp_220 * tmp_37;
      real_t tmp_227 = tmp_223 * tmp_39;
      real_t tmp_228 = tmp_19 * ( 0.041227165399737475 * tmp_43 + 0.1711304259088916 * tmp_44 + tmp_45 );
      real_t tmp_229 = tmp_228 * tmp_41;
      real_t tmp_230 = tmp_228 * tmp_48;
      real_t tmp_231 = tmp_228 * tmp_50;
      real_t tmp_232 = tmp_226 + tmp_227 + tmp_231;
      real_t tmp_233 = tmp_222 + tmp_225 + tmp_230;
      real_t tmp_234 = tmp_221 + tmp_224 + tmp_229;
      real_t tmp_235 = 0.019202922745021479 * tmp_58 *
                       ( tmp_10 * ( tmp_232 - 1.0 / 4.0 ) + tmp_14 * ( tmp_233 - 1.0 / 4.0 ) + tmp_8 * ( tmp_234 - 1.0 / 4.0 ) );
      real_t tmp_236 = tmp_19 * ( 0.40446199974765351 * tmp_21 + 0.19107600050469298 * tmp_22 + tmp_23 );
      real_t tmp_237 = tmp_236 * tmp_6;
      real_t tmp_238 = tmp_236 * tmp_26;
      real_t tmp_239 = tmp_19 * ( 0.40446199974765351 * tmp_30 + 0.19107600050469298 * tmp_31 + tmp_32 );
      real_t tmp_240 = tmp_239 * tmp_28;
      real_t tmp_241 = tmp_239 * tmp_35;
      real_t tmp_242 = tmp_236 * tmp_37;
      real_t tmp_243 = tmp_239 * tmp_39;
      real_t tmp_244 = tmp_19 * ( 0.40446199974765351 * tmp_43 + 0.19107600050469298 * tmp_44 + tmp_45 );
      real_t tmp_245 = tmp_244 * tmp_41;
      real_t tmp_246 = tmp_244 * tmp_48;
      real_t tmp_247 = tmp_244 * tmp_50;
      real_t tmp_248 = tmp_242 + tmp_243 + tmp_247;
      real_t tmp_249 = tmp_238 + tmp_241 + tmp_246;
      real_t tmp_250 = tmp_237 + tmp_240 + tmp_245;
      real_t tmp_251 = 0.042507265838595799 * tmp_58 *
                       ( tmp_10 * ( tmp_248 - 1.0 / 4.0 ) + tmp_14 * ( tmp_249 - 1.0 / 4.0 ) + tmp_8 * ( tmp_250 - 1.0 / 4.0 ) );
      real_t tmp_252 = tmp_19 * ( 0.039308471900058539 * tmp_21 + 0.37605877282253791 * tmp_22 + tmp_23 );
      real_t tmp_253 = tmp_252 * tmp_6;
      real_t tmp_254 = tmp_252 * tmp_26;
      real_t tmp_255 = tmp_19 * ( 0.039308471900058539 * tmp_30 + 0.37605877282253791 * tmp_31 + tmp_32 );
      real_t tmp_256 = tmp_255 * tmp_28;
      real_t tmp_257 = tmp_255 * tmp_35;
      real_t tmp_258 = tmp_252 * tmp_37;
      real_t tmp_259 = tmp_255 * tmp_39;
      real_t tmp_260 = tmp_19 * ( 0.039308471900058539 * tmp_43 + 0.37605877282253791 * tmp_44 + tmp_45 );
      real_t tmp_261 = tmp_260 * tmp_41;
      real_t tmp_262 = tmp_260 * tmp_48;
      real_t tmp_263 = tmp_260 * tmp_50;
      real_t tmp_264 = tmp_258 + tmp_259 + tmp_263;
      real_t tmp_265 = tmp_254 + tmp_257 + tmp_262;
      real_t tmp_266 = tmp_253 + tmp_256 + tmp_261;
      real_t tmp_267 = 0.020848748529055869 * tmp_58 *
                       ( tmp_10 * ( tmp_264 - 1.0 / 4.0 ) + tmp_14 * ( tmp_265 - 1.0 / 4.0 ) + tmp_8 * ( tmp_266 - 1.0 / 4.0 ) );
      real_t tmp_268 = tmp_19 * ( 0.93718850182767688 * tmp_21 + 0.031405749086161582 * tmp_22 + tmp_23 );
      real_t tmp_269 = tmp_268 * tmp_6;
      real_t tmp_270 = tmp_26 * tmp_268;
      real_t tmp_271 = tmp_19 * ( 0.93718850182767688 * tmp_30 + 0.031405749086161582 * tmp_31 + tmp_32 );
      real_t tmp_272 = tmp_271 * tmp_28;
      real_t tmp_273 = tmp_271 * tmp_35;
      real_t tmp_274 = tmp_268 * tmp_37;
      real_t tmp_275 = tmp_271 * tmp_39;
      real_t tmp_276 = tmp_19 * ( 0.93718850182767688 * tmp_43 + 0.031405749086161582 * tmp_44 + tmp_45 );
      real_t tmp_277 = tmp_276 * tmp_41;
      real_t tmp_278 = tmp_276 * tmp_48;
      real_t tmp_279 = tmp_276 * tmp_50;
      real_t tmp_280 = tmp_274 + tmp_275 + tmp_279;
      real_t tmp_281 = tmp_270 + tmp_273 + tmp_278;
      real_t tmp_282 = tmp_269 + tmp_272 + tmp_277;
      real_t tmp_283 = 0.0068572537431980923 * tmp_58 *
                       ( tmp_10 * ( tmp_280 - 1.0 / 4.0 ) + tmp_14 * ( tmp_281 - 1.0 / 4.0 ) + tmp_8 * ( tmp_282 - 1.0 / 4.0 ) );
      real_t tmp_284 = tmp_19 * ( 0.60796128279561268 * tmp_21 + 0.19601935860219369 * tmp_22 + tmp_23 );
      real_t tmp_285 = tmp_284 * tmp_6;
      real_t tmp_286 = tmp_26 * tmp_284;
      real_t tmp_287 = tmp_19 * ( 0.60796128279561268 * tmp_30 + 0.19601935860219369 * tmp_31 + tmp_32 );
      real_t tmp_288 = tmp_28 * tmp_287;
      real_t tmp_289 = tmp_287 * tmp_35;
      real_t tmp_290 = tmp_284 * tmp_37;
      real_t tmp_291 = tmp_287 * tmp_39;
      real_t tmp_292 = tmp_19 * ( 0.60796128279561268 * tmp_43 + 0.19601935860219369 * tmp_44 + tmp_45 );
      real_t tmp_293 = tmp_292 * tmp_41;
      real_t tmp_294 = tmp_292 * tmp_48;
      real_t tmp_295 = tmp_292 * tmp_50;
      real_t tmp_296 = tmp_290 + tmp_291 + tmp_295;
      real_t tmp_297 = tmp_286 + tmp_289 + tmp_294;
      real_t tmp_298 = tmp_285 + tmp_288 + tmp_293;
      real_t tmp_299 = 0.037198804536718075 * tmp_58 *
                       ( tmp_10 * ( tmp_296 - 1.0 / 4.0 ) + tmp_14 * ( tmp_297 - 1.0 / 4.0 ) + tmp_8 * ( tmp_298 - 1.0 / 4.0 ) );
      real_t tmp_300 = tmp_19 * ( 0.19107600050469298 * tmp_21 + 0.40446199974765351 * tmp_22 + tmp_23 );
      real_t tmp_301 = tmp_300 * tmp_6;
      real_t tmp_302 = tmp_26 * tmp_300;
      real_t tmp_303 = tmp_19 * ( 0.19107600050469298 * tmp_30 + 0.40446199974765351 * tmp_31 + tmp_32 );
      real_t tmp_304 = tmp_28 * tmp_303;
      real_t tmp_305 = tmp_303 * tmp_35;
      real_t tmp_306 = tmp_300 * tmp_37;
      real_t tmp_307 = tmp_303 * tmp_39;
      real_t tmp_308 = tmp_19 * ( 0.19107600050469298 * tmp_43 + 0.40446199974765351 * tmp_44 + tmp_45 );
      real_t tmp_309 = tmp_308 * tmp_41;
      real_t tmp_310 = tmp_308 * tmp_48;
      real_t tmp_311 = tmp_308 * tmp_50;
      real_t tmp_312 = tmp_306 + tmp_307 + tmp_311;
      real_t tmp_313 = tmp_302 + tmp_305 + tmp_310;
      real_t tmp_314 = tmp_301 + tmp_304 + tmp_309;
      real_t tmp_315 = 0.042507265838595799 * tmp_58 *
                       ( tmp_10 * ( tmp_312 - 1.0 / 4.0 ) + tmp_14 * ( tmp_313 - 1.0 / 4.0 ) + tmp_8 * ( tmp_314 - 1.0 / 4.0 ) );
      real_t tmp_316 = tmp_19 * ( 0.031405749086161582 * tmp_21 + 0.031405749086161582 * tmp_22 + tmp_23 );
      real_t tmp_317 = tmp_316 * tmp_6;
      real_t tmp_318 = tmp_26 * tmp_316;
      real_t tmp_319 = tmp_19 * ( 0.031405749086161582 * tmp_30 + 0.031405749086161582 * tmp_31 + tmp_32 );
      real_t tmp_320 = tmp_28 * tmp_319;
      real_t tmp_321 = tmp_319 * tmp_35;
      real_t tmp_322 = tmp_316 * tmp_37;
      real_t tmp_323 = tmp_319 * tmp_39;
      real_t tmp_324 = tmp_19 * ( 0.031405749086161582 * tmp_43 + 0.031405749086161582 * tmp_44 + tmp_45 );
      real_t tmp_325 = tmp_324 * tmp_41;
      real_t tmp_326 = tmp_324 * tmp_48;
      real_t tmp_327 = tmp_324 * tmp_50;
      real_t tmp_328 = tmp_322 + tmp_323 + tmp_327;
      real_t tmp_329 = tmp_318 + tmp_321 + tmp_326;
      real_t tmp_330 = tmp_317 + tmp_320 + tmp_325;
      real_t tmp_331 = 0.0068572537431980923 * tmp_58 *
                       ( tmp_10 * ( tmp_328 - 1.0 / 4.0 ) + tmp_14 * ( tmp_329 - 1.0 / 4.0 ) + tmp_8 * ( tmp_330 - 1.0 / 4.0 ) );
      real_t tmp_332 = tmp_19 * ( 0.19601935860219369 * tmp_21 + 0.19601935860219369 * tmp_22 + tmp_23 );
      real_t tmp_333 = tmp_332 * tmp_6;
      real_t tmp_334 = tmp_26 * tmp_332;
      real_t tmp_335 = tmp_19 * ( 0.19601935860219369 * tmp_30 + 0.19601935860219369 * tmp_31 + tmp_32 );
      real_t tmp_336 = tmp_28 * tmp_335;
      real_t tmp_337 = tmp_335 * tmp_35;
      real_t tmp_338 = tmp_332 * tmp_37;
      real_t tmp_339 = tmp_335 * tmp_39;
      real_t tmp_340 = tmp_19 * ( 0.19601935860219369 * tmp_43 + 0.19601935860219369 * tmp_44 + tmp_45 );
      real_t tmp_341 = tmp_340 * tmp_41;
      real_t tmp_342 = tmp_340 * tmp_48;
      real_t tmp_343 = tmp_340 * tmp_50;
      real_t tmp_344 = tmp_338 + tmp_339 + tmp_343;
      real_t tmp_345 = tmp_334 + tmp_337 + tmp_342;
      real_t tmp_346 = tmp_333 + tmp_336 + tmp_341;
      real_t tmp_347 = 0.037198804536718075 * tmp_58 *
                       ( tmp_10 * ( tmp_344 - 1.0 / 4.0 ) + tmp_14 * ( tmp_345 - 1.0 / 4.0 ) + tmp_8 * ( tmp_346 - 1.0 / 4.0 ) );
      real_t tmp_348 = tmp_19 * ( 0.40446199974765351 * tmp_21 + 0.40446199974765351 * tmp_22 + tmp_23 );
      real_t tmp_349 = tmp_348 * tmp_6;
      real_t tmp_350 = tmp_26 * tmp_348;
      real_t tmp_351 = tmp_19 * ( 0.40446199974765351 * tmp_30 + 0.40446199974765351 * tmp_31 + tmp_32 );
      real_t tmp_352 = tmp_28 * tmp_351;
      real_t tmp_353 = tmp_35 * tmp_351;
      real_t tmp_354 = tmp_348 * tmp_37;
      real_t tmp_355 = tmp_351 * tmp_39;
      real_t tmp_356 = tmp_19 * ( 0.40446199974765351 * tmp_43 + 0.40446199974765351 * tmp_44 + tmp_45 );
      real_t tmp_357 = tmp_356 * tmp_41;
      real_t tmp_358 = tmp_356 * tmp_48;
      real_t tmp_359 = tmp_356 * tmp_50;
      real_t tmp_360 = tmp_354 + tmp_355 + tmp_359;
      real_t tmp_361 = tmp_350 + tmp_353 + tmp_358;
      real_t tmp_362 = tmp_349 + tmp_352 + tmp_357;
      real_t tmp_363 = 0.042507265838595799 * tmp_58 *
                       ( tmp_10 * ( tmp_360 - 1.0 / 4.0 ) + tmp_14 * ( tmp_361 - 1.0 / 4.0 ) + tmp_8 * ( tmp_362 - 1.0 / 4.0 ) );
      real_t tmp_364 = tmp_19 * ( 0.1711304259088916 * tmp_21 + 0.041227165399737475 * tmp_22 + tmp_23 );
      real_t tmp_365 = tmp_364 * tmp_6;
      real_t tmp_366 = tmp_26 * tmp_364;
      real_t tmp_367 = tmp_19 * ( 0.1711304259088916 * tmp_30 + 0.041227165399737475 * tmp_31 + tmp_32 );
      real_t tmp_368 = tmp_28 * tmp_367;
      real_t tmp_369 = tmp_35 * tmp_367;
      real_t tmp_370 = tmp_364 * tmp_37;
      real_t tmp_371 = tmp_367 * tmp_39;
      real_t tmp_372 = tmp_19 * ( 0.1711304259088916 * tmp_43 + 0.041227165399737475 * tmp_44 + tmp_45 );
      real_t tmp_373 = tmp_372 * tmp_41;
      real_t tmp_374 = tmp_372 * tmp_48;
      real_t tmp_375 = tmp_372 * tmp_50;
      real_t tmp_376 = tmp_370 + tmp_371 + tmp_375;
      real_t tmp_377 = tmp_366 + tmp_369 + tmp_374;
      real_t tmp_378 = tmp_365 + tmp_368 + tmp_373;
      real_t tmp_379 = 0.019202922745021479 * tmp_58 *
                       ( tmp_10 * ( tmp_376 - 1.0 / 4.0 ) + tmp_14 * ( tmp_377 - 1.0 / 4.0 ) + tmp_8 * ( tmp_378 - 1.0 / 4.0 ) );
      real_t a_0_0 = tmp_107 * ( -tmp_101 - tmp_102 - tmp_103 - tmp_93 - tmp_94 - tmp_96 - tmp_97 - tmp_98 - tmp_99 + 1 ) +
                     tmp_123 * ( -tmp_109 - tmp_110 - tmp_112 - tmp_113 - tmp_114 - tmp_115 - tmp_117 - tmp_118 - tmp_119 + 1 ) +
                     tmp_139 * ( -tmp_125 - tmp_126 - tmp_128 - tmp_129 - tmp_130 - tmp_131 - tmp_133 - tmp_134 - tmp_135 + 1 ) +
                     tmp_155 * ( -tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_147 - tmp_149 - tmp_150 - tmp_151 + 1 ) +
                     tmp_171 * ( -tmp_157 - tmp_158 - tmp_160 - tmp_161 - tmp_162 - tmp_163 - tmp_165 - tmp_166 - tmp_167 + 1 ) +
                     tmp_187 * ( -tmp_173 - tmp_174 - tmp_176 - tmp_177 - tmp_178 - tmp_179 - tmp_181 - tmp_182 - tmp_183 + 1 ) +
                     tmp_203 * ( -tmp_189 - tmp_190 - tmp_192 - tmp_193 - tmp_194 - tmp_195 - tmp_197 - tmp_198 - tmp_199 + 1 ) +
                     tmp_219 * ( -tmp_205 - tmp_206 - tmp_208 - tmp_209 - tmp_210 - tmp_211 - tmp_213 - tmp_214 - tmp_215 + 1 ) +
                     tmp_235 * ( -tmp_221 - tmp_222 - tmp_224 - tmp_225 - tmp_226 - tmp_227 - tmp_229 - tmp_230 - tmp_231 + 1 ) +
                     tmp_251 * ( -tmp_237 - tmp_238 - tmp_240 - tmp_241 - tmp_242 - tmp_243 - tmp_245 - tmp_246 - tmp_247 + 1 ) +
                     tmp_267 * ( -tmp_253 - tmp_254 - tmp_256 - tmp_257 - tmp_258 - tmp_259 - tmp_261 - tmp_262 - tmp_263 + 1 ) +
                     tmp_283 * ( -tmp_269 - tmp_270 - tmp_272 - tmp_273 - tmp_274 - tmp_275 - tmp_277 - tmp_278 - tmp_279 + 1 ) +
                     tmp_299 * ( -tmp_285 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_291 - tmp_293 - tmp_294 - tmp_295 + 1 ) +
                     tmp_315 * ( -tmp_301 - tmp_302 - tmp_304 - tmp_305 - tmp_306 - tmp_307 - tmp_309 - tmp_310 - tmp_311 + 1 ) +
                     tmp_331 * ( -tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 - tmp_325 - tmp_326 - tmp_327 + 1 ) +
                     tmp_347 * ( -tmp_333 - tmp_334 - tmp_336 - tmp_337 - tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 + 1 ) +
                     tmp_363 * ( -tmp_349 - tmp_350 - tmp_352 - tmp_353 - tmp_354 - tmp_355 - tmp_357 - tmp_358 - tmp_359 + 1 ) +
                     tmp_379 * ( -tmp_365 - tmp_366 - tmp_368 - tmp_369 - tmp_370 - tmp_371 - tmp_373 - tmp_374 - tmp_375 + 1 ) +
                     tmp_59 * ( -tmp_25 - tmp_27 - tmp_34 - tmp_36 - tmp_38 - tmp_40 - tmp_47 - tmp_49 - tmp_51 + 1 ) +
                     tmp_75 * ( -tmp_61 - tmp_62 - tmp_64 - tmp_65 - tmp_66 - tmp_67 - tmp_69 - tmp_70 - tmp_71 + 1 ) +
                     tmp_91 * ( -tmp_77 - tmp_78 - tmp_80 - tmp_81 - tmp_82 - tmp_83 - tmp_85 - tmp_86 - tmp_87 + 1 );
      real_t a_1_0 = tmp_104 * tmp_107 + tmp_120 * tmp_123 + tmp_136 * tmp_139 + tmp_152 * tmp_155 + tmp_168 * tmp_171 +
                     tmp_184 * tmp_187 + tmp_200 * tmp_203 + tmp_216 * tmp_219 + tmp_232 * tmp_235 + tmp_248 * tmp_251 +
                     tmp_264 * tmp_267 + tmp_280 * tmp_283 + tmp_296 * tmp_299 + tmp_312 * tmp_315 + tmp_328 * tmp_331 +
                     tmp_344 * tmp_347 + tmp_360 * tmp_363 + tmp_376 * tmp_379 + tmp_52 * tmp_59 + tmp_72 * tmp_75 +
                     tmp_88 * tmp_91;
      real_t a_2_0 = tmp_105 * tmp_107 + tmp_121 * tmp_123 + tmp_137 * tmp_139 + tmp_153 * tmp_155 + tmp_169 * tmp_171 +
                     tmp_185 * tmp_187 + tmp_201 * tmp_203 + tmp_217 * tmp_219 + tmp_233 * tmp_235 + tmp_249 * tmp_251 +
                     tmp_265 * tmp_267 + tmp_281 * tmp_283 + tmp_297 * tmp_299 + tmp_313 * tmp_315 + tmp_329 * tmp_331 +
                     tmp_345 * tmp_347 + tmp_361 * tmp_363 + tmp_377 * tmp_379 + tmp_53 * tmp_59 + tmp_73 * tmp_75 +
                     tmp_89 * tmp_91;
      real_t a_3_0 = tmp_106 * tmp_107 + tmp_122 * tmp_123 + tmp_138 * tmp_139 + tmp_154 * tmp_155 + tmp_170 * tmp_171 +
                     tmp_186 * tmp_187 + tmp_202 * tmp_203 + tmp_218 * tmp_219 + tmp_234 * tmp_235 + tmp_250 * tmp_251 +
                     tmp_266 * tmp_267 + tmp_282 * tmp_283 + tmp_298 * tmp_299 + tmp_314 * tmp_315 + tmp_330 * tmp_331 +
                     tmp_346 * tmp_347 + tmp_362 * tmp_363 + tmp_378 * tmp_379 + tmp_54 * tmp_59 + tmp_74 * tmp_75 +
                     tmp_90 * tmp_91;
      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
      elMat( 3, 0 ) = a_3_0;
   }
};

class EGIIPGVectorLaplaceFormEP1_0 : public hyteg::dg::DGForm
{
 protected:
   void integrateVolume2D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override
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

   virtual void integrateFacetInner2D( const std::vector< Point3D >& coordsElement,
                                       const std::vector< Point3D >& coordsFacet,
                                       const Point3D&                oppositeVertex,
                                       const Point3D&                outwardNormal,
                                       const DGBasisInfo&            trialBasis,
                                       const DGBasisInfo&            testBasis,
                                       int                           trialDegree,
                                       int                           testDegree,
                                       MatrixXr&                     elMat ) const override
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
      real_t tmp_5  = 0.046910077030668018 * tmp_2 + tmp_4;
      real_t tmp_6  = p_affine_2_1 + tmp_3;
      real_t tmp_7  = p_affine_2_0 + tmp_0;
      real_t tmp_8  = 1.0 / ( tmp_1 * tmp_6 - tmp_7 * ( p_affine_1_1 + tmp_3 ) );
      real_t tmp_9  = tmp_8 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_10 = tmp_5 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = 0.046910077030668018 * tmp_11 + tmp_12;
      real_t tmp_14 = tmp_6 * tmp_8;
      real_t tmp_15 = tmp_13 * tmp_14;
      real_t tmp_16 = tmp_10 + tmp_15;
      real_t tmp_17 = tmp_1 * tmp_8;
      real_t tmp_18 = tmp_17 * tmp_5;
      real_t tmp_19 = tmp_8 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_20 = tmp_13 * tmp_19;
      real_t tmp_21 = tmp_18 + tmp_20;
      real_t tmp_22 = tmp_1 * ( tmp_16 - 1.0 / 3.0 ) + tmp_7 * ( tmp_21 - 1.0 / 3.0 );
      real_t tmp_23 = 0.5 * p_affine_10_0 * ( -tmp_14 - tmp_19 ) + 0.5 * p_affine_10_1 * ( -tmp_17 - tmp_9 );
      real_t tmp_24 = std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_2 * tmp_2 ), 1.0 / 2.0 ) );
      real_t tmp_25 = 1.0 / ( tmp_24 );
      real_t tmp_26 = tmp_22 * tmp_25;
      real_t tmp_27 = 0.11846344252809471 * tmp_24;
      real_t tmp_28 = 0.23076534494715845 * tmp_2 + tmp_4;
      real_t tmp_29 = tmp_28 * tmp_9;
      real_t tmp_30 = 0.23076534494715845 * tmp_11 + tmp_12;
      real_t tmp_31 = tmp_14 * tmp_30;
      real_t tmp_32 = tmp_29 + tmp_31;
      real_t tmp_33 = tmp_17 * tmp_28;
      real_t tmp_34 = tmp_19 * tmp_30;
      real_t tmp_35 = tmp_33 + tmp_34;
      real_t tmp_36 = tmp_1 * ( tmp_32 - 1.0 / 3.0 ) + tmp_7 * ( tmp_35 - 1.0 / 3.0 );
      real_t tmp_37 = tmp_25 * tmp_36;
      real_t tmp_38 = 0.2393143352496831 * tmp_24;
      real_t tmp_39 = 0.5 * tmp_2 + tmp_4;
      real_t tmp_40 = tmp_39 * tmp_9;
      real_t tmp_41 = 0.5 * tmp_11 + tmp_12;
      real_t tmp_42 = tmp_14 * tmp_41;
      real_t tmp_43 = tmp_40 + tmp_42;
      real_t tmp_44 = tmp_17 * tmp_39;
      real_t tmp_45 = tmp_19 * tmp_41;
      real_t tmp_46 = tmp_44 + tmp_45;
      real_t tmp_47 = tmp_1 * ( tmp_43 - 1.0 / 3.0 ) + tmp_7 * ( tmp_46 - 1.0 / 3.0 );
      real_t tmp_48 = tmp_25 * tmp_47;
      real_t tmp_49 = 0.2844444444444445 * tmp_24;
      real_t tmp_50 = 0.7692346550528415 * tmp_2 + tmp_4;
      real_t tmp_51 = tmp_50 * tmp_9;
      real_t tmp_52 = 0.7692346550528415 * tmp_11 + tmp_12;
      real_t tmp_53 = tmp_14 * tmp_52;
      real_t tmp_54 = tmp_51 + tmp_53;
      real_t tmp_55 = tmp_17 * tmp_50;
      real_t tmp_56 = tmp_19 * tmp_52;
      real_t tmp_57 = tmp_55 + tmp_56;
      real_t tmp_58 = tmp_1 * ( tmp_54 - 1.0 / 3.0 ) + tmp_7 * ( tmp_57 - 1.0 / 3.0 );
      real_t tmp_59 = tmp_25 * tmp_58;
      real_t tmp_60 = 0.2393143352496831 * tmp_24;
      real_t tmp_61 = 0.95308992296933193 * tmp_2 + tmp_4;
      real_t tmp_62 = tmp_61 * tmp_9;
      real_t tmp_63 = 0.95308992296933193 * tmp_11 + tmp_12;
      real_t tmp_64 = tmp_14 * tmp_63;
      real_t tmp_65 = tmp_62 + tmp_64;
      real_t tmp_66 = tmp_17 * tmp_61;
      real_t tmp_67 = tmp_19 * tmp_63;
      real_t tmp_68 = tmp_66 + tmp_67;
      real_t tmp_69 = tmp_1 * ( tmp_65 - 1.0 / 3.0 ) + tmp_7 * ( tmp_68 - 1.0 / 3.0 );
      real_t tmp_70 = tmp_25 * tmp_69;
      real_t tmp_71 = 0.11846344252809471 * tmp_24;
      real_t tmp_72 = 0.5 * p_affine_10_0 * tmp_14 + 0.5 * p_affine_10_1 * tmp_9;
      real_t tmp_73 = 0.5 * p_affine_10_0 * tmp_19 + 0.5 * p_affine_10_1 * tmp_17;
      real_t a_0_0  = tmp_27 * ( -tmp_22 * tmp_23 + tmp_26 * ( -tmp_10 - tmp_15 - tmp_18 - tmp_20 + 1 ) ) +
                     tmp_38 * ( -tmp_23 * tmp_36 + tmp_37 * ( -tmp_29 - tmp_31 - tmp_33 - tmp_34 + 1 ) ) +
                     tmp_49 * ( -tmp_23 * tmp_47 + tmp_48 * ( -tmp_40 - tmp_42 - tmp_44 - tmp_45 + 1 ) ) +
                     tmp_60 * ( -tmp_23 * tmp_58 + tmp_59 * ( -tmp_51 - tmp_53 - tmp_55 - tmp_56 + 1 ) ) +
                     tmp_71 * ( -tmp_23 * tmp_69 + tmp_70 * ( -tmp_62 - tmp_64 - tmp_66 - tmp_67 + 1 ) );
      real_t a_0_1 = tmp_27 * ( tmp_16 * tmp_26 - tmp_22 * tmp_72 ) + tmp_38 * ( tmp_32 * tmp_37 - tmp_36 * tmp_72 ) +
                     tmp_49 * ( tmp_43 * tmp_48 - tmp_47 * tmp_72 ) + tmp_60 * ( tmp_54 * tmp_59 - tmp_58 * tmp_72 ) +
                     tmp_71 * ( tmp_65 * tmp_70 - tmp_69 * tmp_72 );
      real_t a_0_2 = tmp_27 * ( tmp_21 * tmp_26 - tmp_22 * tmp_73 ) + tmp_38 * ( tmp_35 * tmp_37 - tmp_36 * tmp_73 ) +
                     tmp_49 * ( tmp_46 * tmp_48 - tmp_47 * tmp_73 ) + tmp_60 * ( tmp_57 * tmp_59 - tmp_58 * tmp_73 ) +
                     tmp_71 * ( tmp_68 * tmp_70 - tmp_69 * tmp_73 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >& coordsElementInner,
                                          const std::vector< Point3D >& coordsElementOuter,
                                          const std::vector< Point3D >& coordsFacet,
                                          const Point3D&                oppositeVertexInnerElement,
                                          const Point3D&                oppositeVertexOuterElement,
                                          const Point3D&                outwardNormal,
                                          const DGBasisInfo&            trialBasis,
                                          const DGBasisInfo&            testBasis,
                                          int                           trialDegree,
                                          int                           testDegree,
                                          MatrixXr&                     elMat ) const override
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
      real_t tmp_2  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3  = -p_affine_0_1;
      real_t tmp_4  = p_affine_2_1 + tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_0;
      real_t tmp_6  = 1.0 / ( tmp_1 * tmp_4 - tmp_5 * ( p_affine_1_1 + tmp_3 ) );
      real_t tmp_7  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_8  = p_affine_6_1 + 0.046910077030668018 * tmp_7;
      real_t tmp_9  = tmp_6 * ( tmp_3 + tmp_8 );
      real_t tmp_10 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_11 = p_affine_6_0 + 0.046910077030668018 * tmp_10;
      real_t tmp_12 = tmp_6 * ( tmp_0 + tmp_11 );
      real_t tmp_13 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_14 =
          tmp_1 * ( tmp_12 * tmp_4 + tmp_2 * tmp_9 - 1.0 / 3.0 ) + tmp_5 * ( tmp_1 * tmp_9 + tmp_12 * tmp_13 - 1.0 / 3.0 );
      real_t tmp_15 = -p_affine_3_1;
      real_t tmp_16 = p_affine_5_1 + tmp_15;
      real_t tmp_17 = -p_affine_3_0;
      real_t tmp_18 = p_affine_4_0 + tmp_17;
      real_t tmp_19 = 1.0 / ( tmp_16 * tmp_18 - ( p_affine_4_1 + tmp_15 ) * ( p_affine_5_0 + tmp_17 ) );
      real_t tmp_20 = tmp_16 * tmp_19;
      real_t tmp_21 = tmp_19 * ( p_affine_3_1 - p_affine_4_1 );
      real_t tmp_22 = tmp_18 * tmp_19;
      real_t tmp_23 = tmp_19 * ( p_affine_3_0 - p_affine_5_0 );
      real_t tmp_24 = 0.5 * p_affine_10_0 * ( -tmp_20 - tmp_21 ) + 0.5 * p_affine_10_1 * ( -tmp_22 - tmp_23 );
      real_t tmp_25 = tmp_15 + tmp_8;
      real_t tmp_26 = tmp_22 * tmp_25;
      real_t tmp_27 = tmp_23 * tmp_25;
      real_t tmp_28 = tmp_11 + tmp_17;
      real_t tmp_29 = tmp_20 * tmp_28;
      real_t tmp_30 = tmp_21 * tmp_28;
      real_t tmp_31 = std::abs( std::pow( ( tmp_10 * tmp_10 ) + ( tmp_7 * tmp_7 ), 1.0 / 2.0 ) );
      real_t tmp_32 = 1.0 / ( tmp_31 );
      real_t tmp_33 = tmp_14 * tmp_32;
      real_t tmp_34 = 0.11846344252809471 * tmp_31;
      real_t tmp_35 = p_affine_6_1 + 0.23076534494715845 * tmp_7;
      real_t tmp_36 = tmp_6 * ( tmp_3 + tmp_35 );
      real_t tmp_37 = p_affine_6_0 + 0.23076534494715845 * tmp_10;
      real_t tmp_38 = tmp_6 * ( tmp_0 + tmp_37 );
      real_t tmp_39 =
          tmp_1 * ( tmp_2 * tmp_36 + tmp_38 * tmp_4 - 1.0 / 3.0 ) + tmp_5 * ( tmp_1 * tmp_36 + tmp_13 * tmp_38 - 1.0 / 3.0 );
      real_t tmp_40 = tmp_15 + tmp_35;
      real_t tmp_41 = tmp_22 * tmp_40;
      real_t tmp_42 = tmp_23 * tmp_40;
      real_t tmp_43 = tmp_17 + tmp_37;
      real_t tmp_44 = tmp_20 * tmp_43;
      real_t tmp_45 = tmp_21 * tmp_43;
      real_t tmp_46 = tmp_32 * tmp_39;
      real_t tmp_47 = 0.2393143352496831 * tmp_31;
      real_t tmp_48 = p_affine_6_1 + 0.5 * tmp_7;
      real_t tmp_49 = tmp_6 * ( tmp_3 + tmp_48 );
      real_t tmp_50 = p_affine_6_0 + 0.5 * tmp_10;
      real_t tmp_51 = tmp_6 * ( tmp_0 + tmp_50 );
      real_t tmp_52 =
          tmp_1 * ( tmp_2 * tmp_49 + tmp_4 * tmp_51 - 1.0 / 3.0 ) + tmp_5 * ( tmp_1 * tmp_49 + tmp_13 * tmp_51 - 1.0 / 3.0 );
      real_t tmp_53 = tmp_15 + tmp_48;
      real_t tmp_54 = tmp_22 * tmp_53;
      real_t tmp_55 = tmp_23 * tmp_53;
      real_t tmp_56 = tmp_17 + tmp_50;
      real_t tmp_57 = tmp_20 * tmp_56;
      real_t tmp_58 = tmp_21 * tmp_56;
      real_t tmp_59 = tmp_32 * tmp_52;
      real_t tmp_60 = 0.2844444444444445 * tmp_31;
      real_t tmp_61 = p_affine_6_1 + 0.7692346550528415 * tmp_7;
      real_t tmp_62 = tmp_6 * ( tmp_3 + tmp_61 );
      real_t tmp_63 = p_affine_6_0 + 0.7692346550528415 * tmp_10;
      real_t tmp_64 = tmp_6 * ( tmp_0 + tmp_63 );
      real_t tmp_65 =
          tmp_1 * ( tmp_2 * tmp_62 + tmp_4 * tmp_64 - 1.0 / 3.0 ) + tmp_5 * ( tmp_1 * tmp_62 + tmp_13 * tmp_64 - 1.0 / 3.0 );
      real_t tmp_66 = tmp_15 + tmp_61;
      real_t tmp_67 = tmp_22 * tmp_66;
      real_t tmp_68 = tmp_23 * tmp_66;
      real_t tmp_69 = tmp_17 + tmp_63;
      real_t tmp_70 = tmp_20 * tmp_69;
      real_t tmp_71 = tmp_21 * tmp_69;
      real_t tmp_72 = tmp_32 * tmp_65;
      real_t tmp_73 = 0.2393143352496831 * tmp_31;
      real_t tmp_74 = p_affine_6_1 + 0.95308992296933193 * tmp_7;
      real_t tmp_75 = tmp_6 * ( tmp_3 + tmp_74 );
      real_t tmp_76 = p_affine_6_0 + 0.95308992296933193 * tmp_10;
      real_t tmp_77 = tmp_6 * ( tmp_0 + tmp_76 );
      real_t tmp_78 =
          tmp_1 * ( tmp_2 * tmp_75 + tmp_4 * tmp_77 - 1.0 / 3.0 ) + tmp_5 * ( tmp_1 * tmp_75 + tmp_13 * tmp_77 - 1.0 / 3.0 );
      real_t tmp_79 = tmp_15 + tmp_74;
      real_t tmp_80 = tmp_22 * tmp_79;
      real_t tmp_81 = tmp_23 * tmp_79;
      real_t tmp_82 = tmp_17 + tmp_76;
      real_t tmp_83 = tmp_20 * tmp_82;
      real_t tmp_84 = tmp_21 * tmp_82;
      real_t tmp_85 = tmp_32 * tmp_78;
      real_t tmp_86 = 0.11846344252809471 * tmp_31;
      real_t tmp_87 = 0.5 * p_affine_10_0 * tmp_20 + 0.5 * p_affine_10_1 * tmp_23;
      real_t tmp_88 = 0.5 * p_affine_10_0 * tmp_21 + 0.5 * p_affine_10_1 * tmp_22;
      real_t a_0_0  = tmp_34 * ( -tmp_14 * tmp_24 - tmp_33 * ( -tmp_26 - tmp_27 - tmp_29 - tmp_30 + 1 ) ) +
                     tmp_47 * ( -tmp_24 * tmp_39 - tmp_46 * ( -tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1 ) ) +
                     tmp_60 * ( -tmp_24 * tmp_52 - tmp_59 * ( -tmp_54 - tmp_55 - tmp_57 - tmp_58 + 1 ) ) +
                     tmp_73 * ( -tmp_24 * tmp_65 - tmp_72 * ( -tmp_67 - tmp_68 - tmp_70 - tmp_71 + 1 ) ) +
                     tmp_86 * ( -tmp_24 * tmp_78 - tmp_85 * ( -tmp_80 - tmp_81 - tmp_83 - tmp_84 + 1 ) );
      real_t a_0_1 = tmp_34 * ( -tmp_14 * tmp_87 - tmp_33 * ( tmp_27 + tmp_29 ) ) +
                     tmp_47 * ( -tmp_39 * tmp_87 - tmp_46 * ( tmp_42 + tmp_44 ) ) +
                     tmp_60 * ( -tmp_52 * tmp_87 - tmp_59 * ( tmp_55 + tmp_57 ) ) +
                     tmp_73 * ( -tmp_65 * tmp_87 - tmp_72 * ( tmp_68 + tmp_70 ) ) +
                     tmp_86 * ( -tmp_78 * tmp_87 - tmp_85 * ( tmp_81 + tmp_83 ) );
      real_t a_0_2 = tmp_34 * ( -tmp_14 * tmp_88 - tmp_33 * ( tmp_26 + tmp_30 ) ) +
                     tmp_47 * ( -tmp_39 * tmp_88 - tmp_46 * ( tmp_41 + tmp_45 ) ) +
                     tmp_60 * ( -tmp_52 * tmp_88 - tmp_59 * ( tmp_54 + tmp_58 ) ) +
                     tmp_73 * ( -tmp_65 * tmp_88 - tmp_72 * ( tmp_67 + tmp_71 ) ) +
                     tmp_86 * ( -tmp_78 * tmp_88 - tmp_85 * ( tmp_80 + tmp_84 ) );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                                   const std::vector< Point3D >& coordsFacet,
                                                   const Point3D&                oppositeVertex,
                                                   const Point3D&                outwardNormal,
                                                   const DGBasisInfo&            trialBasis,
                                                   const DGBasisInfo&            testBasis,
                                                   int                           trialDegree,
                                                   int                           testDegree,
                                                   MatrixXr&                     elMat ) const override
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
      real_t tmp_4  = p_affine_2_0 + tmp_2;
      real_t tmp_5  = 1.0 / ( tmp_1 * tmp_3 - tmp_4 * ( p_affine_1_1 + tmp_0 ) );
      real_t tmp_6  = tmp_1 * tmp_5;
      real_t tmp_7  = tmp_5 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_8  = tmp_3 * tmp_5;
      real_t tmp_9  = tmp_5 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_10 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_11 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_12 = std::abs( std::pow( ( tmp_10 * tmp_10 ) + ( tmp_11 * tmp_11 ), 1.0 / 2.0 ) );
      real_t tmp_13 = tmp_12 * ( p_affine_10_0 * ( -tmp_6 - tmp_7 ) + p_affine_10_1 * ( -tmp_8 - tmp_9 ) );
      real_t tmp_14 = p_affine_6_1 + tmp_0;
      real_t tmp_15 = 0.046910077030668018 * tmp_11 + tmp_14;
      real_t tmp_16 = p_affine_6_0 + tmp_2;
      real_t tmp_17 = 0.046910077030668018 * tmp_10 + tmp_16;
      real_t tmp_18 = 0.11846344252809471 * tmp_3 * ( tmp_15 * tmp_9 + tmp_17 * tmp_6 - 1.0 / 3.0 ) +
                      0.11846344252809471 * tmp_4 * ( tmp_15 * tmp_8 + tmp_17 * tmp_7 - 1.0 / 3.0 );
      real_t tmp_19 = 0.23076534494715845 * tmp_11 + tmp_14;
      real_t tmp_20 = 0.23076534494715845 * tmp_10 + tmp_16;
      real_t tmp_21 = 0.2393143352496831 * tmp_3 * ( tmp_19 * tmp_9 + tmp_20 * tmp_6 - 1.0 / 3.0 ) +
                      0.2393143352496831 * tmp_4 * ( tmp_19 * tmp_8 + tmp_20 * tmp_7 - 1.0 / 3.0 );
      real_t tmp_22 = 0.5 * tmp_11 + tmp_14;
      real_t tmp_23 = 0.5 * tmp_10 + tmp_16;
      real_t tmp_24 = 0.2844444444444445 * tmp_3 * ( tmp_22 * tmp_9 + tmp_23 * tmp_6 - 1.0 / 3.0 ) +
                      0.2844444444444445 * tmp_4 * ( tmp_22 * tmp_8 + tmp_23 * tmp_7 - 1.0 / 3.0 );
      real_t tmp_25 = 0.7692346550528415 * tmp_11 + tmp_14;
      real_t tmp_26 = 0.7692346550528415 * tmp_10 + tmp_16;
      real_t tmp_27 = 0.2393143352496831 * tmp_3 * ( tmp_25 * tmp_9 + tmp_26 * tmp_6 - 1.0 / 3.0 ) +
                      0.2393143352496831 * tmp_4 * ( tmp_25 * tmp_8 + tmp_26 * tmp_7 - 1.0 / 3.0 );
      real_t tmp_28 = 0.95308992296933193 * tmp_11 + tmp_14;
      real_t tmp_29 = 0.95308992296933193 * tmp_10 + tmp_16;
      real_t tmp_30 = 0.11846344252809471 * tmp_3 * ( tmp_28 * tmp_9 + tmp_29 * tmp_6 - 1.0 / 3.0 ) +
                      0.11846344252809471 * tmp_4 * ( tmp_28 * tmp_8 + tmp_29 * tmp_7 - 1.0 / 3.0 );
      real_t tmp_31 = tmp_12 * ( p_affine_10_0 * tmp_6 + p_affine_10_1 * tmp_9 );
      real_t tmp_32 = tmp_12 * ( p_affine_10_0 * tmp_7 + p_affine_10_1 * tmp_8 );
      real_t a_0_0  = -tmp_13 * tmp_18 - tmp_13 * tmp_21 - tmp_13 * tmp_24 - tmp_13 * tmp_27 - tmp_13 * tmp_30;
      real_t a_0_1  = -tmp_18 * tmp_31 - tmp_21 * tmp_31 - tmp_24 * tmp_31 - tmp_27 * tmp_31 - tmp_30 * tmp_31;
      real_t a_0_2  = -tmp_18 * tmp_32 - tmp_21 * tmp_32 - tmp_24 * tmp_32 - tmp_27 * tmp_32 - tmp_30 * tmp_32;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
   }

   void integrateRHSDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                         const std::vector< Point3D >& coordsFacet,
                                         const Point3D&                oppositeVertex,
                                         const Point3D&                outwardNormal,
                                         const DGBasisInfo&            basis,
                                         int                           degree,
                                         MatrixXr&                     elMat ) const override
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
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                         const std::vector< Point3D >& coordsFacet,
                                         const Point3D&                oppositeVertex,
                                         const Point3D&                outwardNormal,
                                         const DGBasisInfo&            basis,
                                         int                           degree,
                                         MatrixXr&                     elMat ) const override
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
   void integrateVolume3D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_0;
      real_t tmp_6  = p_affine_1_1 + tmp_2;
      real_t tmp_7  = tmp_5 * tmp_6;
      real_t tmp_8  = -p_affine_0_2;
      real_t tmp_9  = p_affine_3_2 + tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_8;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11 * tmp_5;
      real_t tmp_13 = p_affine_2_2 + tmp_8;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = tmp_14 * tmp_6;
      real_t tmp_16 = tmp_1 * tmp_11;
      real_t tmp_17 = tmp_14 * tmp_3;
      real_t tmp_18 =
          1.0 / ( tmp_10 * tmp_12 - tmp_10 * tmp_17 + tmp_13 * tmp_15 - tmp_13 * tmp_16 + tmp_4 * tmp_9 - tmp_7 * tmp_9 );
      real_t tmp_19 = tmp_18 * ( tmp_4 - tmp_7 );
      real_t tmp_20 = tmp_18 * ( tmp_15 - tmp_16 );
      real_t tmp_21 = tmp_18 * ( tmp_12 - tmp_17 );
      real_t tmp_22 = tmp_1 * tmp_21 + tmp_14 * tmp_19 + tmp_20 * tmp_5;
      real_t tmp_23 = tmp_18 * ( -tmp_1 * tmp_13 + tmp_10 * tmp_5 );
      real_t tmp_24 = tmp_18 * ( tmp_1 * tmp_9 - tmp_10 * tmp_14 );
      real_t tmp_25 = tmp_18 * ( tmp_13 * tmp_14 - tmp_5 * tmp_9 );
      real_t tmp_26 = tmp_1 * tmp_25 + tmp_14 * tmp_23 + tmp_24 * tmp_5;
      real_t tmp_27 = tmp_18 * ( -tmp_10 * tmp_3 + tmp_13 * tmp_6 );
      real_t tmp_28 = tmp_18 * ( tmp_10 * tmp_11 - tmp_6 * tmp_9 );
      real_t tmp_29 = tmp_18 * ( -tmp_11 * tmp_13 + tmp_3 * tmp_9 );
      real_t tmp_30 = tmp_1 * tmp_29 + tmp_14 * tmp_27 + tmp_28 * tmp_5;
      real_t tmp_31 = p_affine_0_0 * p_affine_1_1;
      real_t tmp_32 = p_affine_0_0 * p_affine_1_2;
      real_t tmp_33 = p_affine_2_1 * p_affine_3_2;
      real_t tmp_34 = p_affine_0_1 * p_affine_1_0;
      real_t tmp_35 = p_affine_0_1 * p_affine_1_2;
      real_t tmp_36 = p_affine_2_2 * p_affine_3_0;
      real_t tmp_37 = p_affine_0_2 * p_affine_1_0;
      real_t tmp_38 = p_affine_0_2 * p_affine_1_1;
      real_t tmp_39 = p_affine_2_0 * p_affine_3_1;
      real_t tmp_40 = p_affine_2_2 * p_affine_3_1;
      real_t tmp_41 = p_affine_2_0 * p_affine_3_2;
      real_t tmp_42 = p_affine_2_1 * p_affine_3_0;
      real_t tmp_43 = std::abs( p_affine_0_0 * tmp_33 - p_affine_0_0 * tmp_40 + p_affine_0_1 * tmp_36 - p_affine_0_1 * tmp_41 +
                                p_affine_0_2 * tmp_39 - p_affine_0_2 * tmp_42 - p_affine_1_0 * tmp_33 + p_affine_1_0 * tmp_40 -
                                p_affine_1_1 * tmp_36 + p_affine_1_1 * tmp_41 - p_affine_1_2 * tmp_39 + p_affine_1_2 * tmp_42 +
                                p_affine_2_0 * tmp_35 - p_affine_2_0 * tmp_38 - p_affine_2_1 * tmp_32 + p_affine_2_1 * tmp_37 +
                                p_affine_2_2 * tmp_31 - p_affine_2_2 * tmp_34 - p_affine_3_0 * tmp_35 + p_affine_3_0 * tmp_38 +
                                p_affine_3_1 * tmp_32 - p_affine_3_1 * tmp_37 - p_affine_3_2 * tmp_31 + p_affine_3_2 * tmp_34 );
      real_t tmp_44 = tmp_43 * ( tmp_22 * ( -tmp_19 - tmp_20 - tmp_21 ) + tmp_26 * ( -tmp_23 - tmp_24 - tmp_25 ) +
                                 tmp_30 * ( -tmp_27 - tmp_28 - tmp_29 ) );
      real_t tmp_45 = tmp_43 * ( tmp_21 * tmp_22 + tmp_25 * tmp_26 + tmp_29 * tmp_30 );
      real_t tmp_46 = tmp_43 * ( tmp_20 * tmp_22 + tmp_24 * tmp_26 + tmp_28 * tmp_30 );
      real_t tmp_47 = tmp_43 * ( tmp_19 * tmp_22 + tmp_23 * tmp_26 + tmp_27 * tmp_30 );
      real_t a_0_0  = 0.1666666666666668 * tmp_44;
      real_t a_0_1  = 0.1666666666666668 * tmp_45;
      real_t a_0_2  = 0.1666666666666668 * tmp_46;
      real_t a_0_3  = 0.1666666666666668 * tmp_47;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 0, 3 ) = a_0_3;
   }

   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                               const std::vector< Point3D >& coordsFacet,
                               const Point3D&,
                               const Point3D&     outwardNormal,
                               const DGBasisInfo& trialBasis,
                               const DGBasisInfo& testBasis,
                               int                trialDegree,
                               int                testDegree,
                               MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_8_2;
      real_t tmp_3  = p_affine_9_2 + tmp_2;
      real_t tmp_4  = p_affine_10_2 + tmp_2;
      real_t tmp_5  = -p_affine_0_2;
      real_t tmp_6  = p_affine_8_2 + tmp_5;
      real_t tmp_7  = 0.031405749086161582 * tmp_3 + 0.93718850182767688 * tmp_4 + tmp_6;
      real_t tmp_8  = p_affine_2_0 + tmp_0;
      real_t tmp_9  = -p_affine_0_1;
      real_t tmp_10 = p_affine_3_1 + tmp_9;
      real_t tmp_11 = p_affine_3_0 + tmp_0;
      real_t tmp_12 = p_affine_2_1 + tmp_9;
      real_t tmp_13 = p_affine_3_2 + tmp_5;
      real_t tmp_14 = tmp_12 * tmp_13;
      real_t tmp_15 = p_affine_1_2 + tmp_5;
      real_t tmp_16 = tmp_10 * tmp_15;
      real_t tmp_17 = p_affine_1_1 + tmp_9;
      real_t tmp_18 = p_affine_2_2 + tmp_5;
      real_t tmp_19 = tmp_17 * tmp_18;
      real_t tmp_20 = tmp_10 * tmp_18;
      real_t tmp_21 = tmp_13 * tmp_17;
      real_t tmp_22 = tmp_12 * tmp_15;
      real_t tmp_23 =
          1.0 / ( tmp_1 * tmp_14 - tmp_1 * tmp_20 + tmp_11 * tmp_19 - tmp_11 * tmp_22 + tmp_16 * tmp_8 - tmp_21 * tmp_8 );
      real_t tmp_24 = tmp_23 * ( tmp_10 * tmp_8 - tmp_11 * tmp_12 );
      real_t tmp_25 = tmp_24 * tmp_7;
      real_t tmp_26 = -p_affine_8_1;
      real_t tmp_27 = p_affine_9_1 + tmp_26;
      real_t tmp_28 = p_affine_10_1 + tmp_26;
      real_t tmp_29 = p_affine_8_1 + tmp_9;
      real_t tmp_30 = 0.031405749086161582 * tmp_27 + 0.93718850182767688 * tmp_28 + tmp_29;
      real_t tmp_31 = tmp_23 * ( tmp_11 * tmp_18 - tmp_13 * tmp_8 );
      real_t tmp_32 = tmp_30 * tmp_31;
      real_t tmp_33 = -p_affine_8_0;
      real_t tmp_34 = p_affine_9_0 + tmp_33;
      real_t tmp_35 = p_affine_10_0 + tmp_33;
      real_t tmp_36 = p_affine_8_0 + tmp_0;
      real_t tmp_37 = 0.031405749086161582 * tmp_34 + 0.93718850182767688 * tmp_35 + tmp_36;
      real_t tmp_38 = tmp_23 * ( tmp_14 - tmp_20 );
      real_t tmp_39 = tmp_37 * tmp_38;
      real_t tmp_40 = tmp_25 + tmp_32 + tmp_39;
      real_t tmp_41 = tmp_23 * ( -tmp_1 * tmp_10 + tmp_11 * tmp_17 );
      real_t tmp_42 = tmp_41 * tmp_7;
      real_t tmp_43 = tmp_23 * ( tmp_1 * tmp_13 - tmp_11 * tmp_15 );
      real_t tmp_44 = tmp_30 * tmp_43;
      real_t tmp_45 = tmp_23 * ( tmp_16 - tmp_21 );
      real_t tmp_46 = tmp_37 * tmp_45;
      real_t tmp_47 = tmp_42 + tmp_44 + tmp_46;
      real_t tmp_48 = tmp_23 * ( tmp_1 * tmp_12 - tmp_17 * tmp_8 );
      real_t tmp_49 = tmp_48 * tmp_7;
      real_t tmp_50 = tmp_23 * ( -tmp_1 * tmp_18 + tmp_15 * tmp_8 );
      real_t tmp_51 = tmp_30 * tmp_50;
      real_t tmp_52 = tmp_23 * ( tmp_19 - tmp_22 );
      real_t tmp_53 = tmp_37 * tmp_52;
      real_t tmp_54 = tmp_49 + tmp_51 + tmp_53;
      real_t tmp_55 = tmp_1 * ( tmp_40 - 1.0 / 4.0 ) + tmp_11 * ( tmp_54 - 1.0 / 4.0 ) + tmp_8 * ( tmp_47 - 1.0 / 4.0 );
      real_t tmp_56 = 0.5 * p_affine_13_0 * ( -tmp_38 - tmp_45 - tmp_52 ) + 0.5 * p_affine_13_1 * ( -tmp_31 - tmp_43 - tmp_50 ) +
                      0.5 * p_affine_13_2 * ( -tmp_24 - tmp_41 - tmp_48 );
      real_t tmp_57 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_58 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_59 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_60 = ( std::abs( tmp_28 * tmp_58 - tmp_35 * tmp_57 ) * std::abs( tmp_28 * tmp_58 - tmp_35 * tmp_57 ) ) +
                      ( std::abs( tmp_28 * tmp_59 - tmp_4 * tmp_57 ) * std::abs( tmp_28 * tmp_59 - tmp_4 * tmp_57 ) ) +
                      ( std::abs( tmp_35 * tmp_59 - tmp_4 * tmp_58 ) * std::abs( tmp_35 * tmp_59 - tmp_4 * tmp_58 ) );
      real_t tmp_61  = 1.0 * std::pow( tmp_60, -0.25 );
      real_t tmp_62  = tmp_55 * tmp_61;
      real_t tmp_63  = 1.0 * std::pow( tmp_60, 1.0 / 2.0 );
      real_t tmp_64  = 0.0068572537431980923 * tmp_63;
      real_t tmp_65  = 0.19601935860219369 * tmp_3 + 0.60796128279561268 * tmp_4 + tmp_6;
      real_t tmp_66  = tmp_24 * tmp_65;
      real_t tmp_67  = 0.19601935860219369 * tmp_27 + 0.60796128279561268 * tmp_28 + tmp_29;
      real_t tmp_68  = tmp_31 * tmp_67;
      real_t tmp_69  = 0.19601935860219369 * tmp_34 + 0.60796128279561268 * tmp_35 + tmp_36;
      real_t tmp_70  = tmp_38 * tmp_69;
      real_t tmp_71  = tmp_66 + tmp_68 + tmp_70;
      real_t tmp_72  = tmp_41 * tmp_65;
      real_t tmp_73  = tmp_43 * tmp_67;
      real_t tmp_74  = tmp_45 * tmp_69;
      real_t tmp_75  = tmp_72 + tmp_73 + tmp_74;
      real_t tmp_76  = tmp_48 * tmp_65;
      real_t tmp_77  = tmp_50 * tmp_67;
      real_t tmp_78  = tmp_52 * tmp_69;
      real_t tmp_79  = tmp_76 + tmp_77 + tmp_78;
      real_t tmp_80  = tmp_1 * ( tmp_71 - 1.0 / 4.0 ) + tmp_11 * ( tmp_79 - 1.0 / 4.0 ) + tmp_8 * ( tmp_75 - 1.0 / 4.0 );
      real_t tmp_81  = tmp_61 * tmp_80;
      real_t tmp_82  = 0.037198804536718075 * tmp_63;
      real_t tmp_83  = 0.37605877282253791 * tmp_3 + 0.039308471900058539 * tmp_4 + tmp_6;
      real_t tmp_84  = tmp_24 * tmp_83;
      real_t tmp_85  = 0.37605877282253791 * tmp_27 + 0.039308471900058539 * tmp_28 + tmp_29;
      real_t tmp_86  = tmp_31 * tmp_85;
      real_t tmp_87  = 0.37605877282253791 * tmp_34 + 0.039308471900058539 * tmp_35 + tmp_36;
      real_t tmp_88  = tmp_38 * tmp_87;
      real_t tmp_89  = tmp_84 + tmp_86 + tmp_88;
      real_t tmp_90  = tmp_41 * tmp_83;
      real_t tmp_91  = tmp_43 * tmp_85;
      real_t tmp_92  = tmp_45 * tmp_87;
      real_t tmp_93  = tmp_90 + tmp_91 + tmp_92;
      real_t tmp_94  = tmp_48 * tmp_83;
      real_t tmp_95  = tmp_50 * tmp_85;
      real_t tmp_96  = tmp_52 * tmp_87;
      real_t tmp_97  = tmp_94 + tmp_95 + tmp_96;
      real_t tmp_98  = tmp_1 * ( tmp_89 - 1.0 / 4.0 ) + tmp_11 * ( tmp_97 - 1.0 / 4.0 ) + tmp_8 * ( tmp_93 - 1.0 / 4.0 );
      real_t tmp_99  = tmp_61 * tmp_98;
      real_t tmp_100 = 0.020848748529055869 * tmp_63;
      real_t tmp_101 = 0.78764240869137092 * tmp_3 + 0.1711304259088916 * tmp_4 + tmp_6;
      real_t tmp_102 = tmp_101 * tmp_24;
      real_t tmp_103 = 0.78764240869137092 * tmp_27 + 0.1711304259088916 * tmp_28 + tmp_29;
      real_t tmp_104 = tmp_103 * tmp_31;
      real_t tmp_105 = 0.78764240869137092 * tmp_34 + 0.1711304259088916 * tmp_35 + tmp_36;
      real_t tmp_106 = tmp_105 * tmp_38;
      real_t tmp_107 = tmp_102 + tmp_104 + tmp_106;
      real_t tmp_108 = tmp_101 * tmp_41;
      real_t tmp_109 = tmp_103 * tmp_43;
      real_t tmp_110 = tmp_105 * tmp_45;
      real_t tmp_111 = tmp_108 + tmp_109 + tmp_110;
      real_t tmp_112 = tmp_101 * tmp_48;
      real_t tmp_113 = tmp_103 * tmp_50;
      real_t tmp_114 = tmp_105 * tmp_52;
      real_t tmp_115 = tmp_112 + tmp_113 + tmp_114;
      real_t tmp_116 = tmp_1 * ( tmp_107 - 1.0 / 4.0 ) + tmp_11 * ( tmp_115 - 1.0 / 4.0 ) + tmp_8 * ( tmp_111 - 1.0 / 4.0 );
      real_t tmp_117 = tmp_116 * tmp_61;
      real_t tmp_118 = 0.019202922745021479 * tmp_63;
      real_t tmp_119 = 0.58463275527740355 * tmp_3 + 0.37605877282253791 * tmp_4 + tmp_6;
      real_t tmp_120 = tmp_119 * tmp_24;
      real_t tmp_121 = 0.58463275527740355 * tmp_27 + 0.37605877282253791 * tmp_28 + tmp_29;
      real_t tmp_122 = tmp_121 * tmp_31;
      real_t tmp_123 = 0.58463275527740355 * tmp_34 + 0.37605877282253791 * tmp_35 + tmp_36;
      real_t tmp_124 = tmp_123 * tmp_38;
      real_t tmp_125 = tmp_120 + tmp_122 + tmp_124;
      real_t tmp_126 = tmp_119 * tmp_41;
      real_t tmp_127 = tmp_121 * tmp_43;
      real_t tmp_128 = tmp_123 * tmp_45;
      real_t tmp_129 = tmp_126 + tmp_127 + tmp_128;
      real_t tmp_130 = tmp_119 * tmp_48;
      real_t tmp_131 = tmp_121 * tmp_50;
      real_t tmp_132 = tmp_123 * tmp_52;
      real_t tmp_133 = tmp_130 + tmp_131 + tmp_132;
      real_t tmp_134 = tmp_1 * ( tmp_125 - 1.0 / 4.0 ) + tmp_11 * ( tmp_133 - 1.0 / 4.0 ) + tmp_8 * ( tmp_129 - 1.0 / 4.0 );
      real_t tmp_135 = tmp_134 * tmp_61;
      real_t tmp_136 = 0.020848748529055869 * tmp_63;
      real_t tmp_137 = 0.041227165399737475 * tmp_3 + 0.78764240869137092 * tmp_4 + tmp_6;
      real_t tmp_138 = tmp_137 * tmp_24;
      real_t tmp_139 = 0.041227165399737475 * tmp_27 + 0.78764240869137092 * tmp_28 + tmp_29;
      real_t tmp_140 = tmp_139 * tmp_31;
      real_t tmp_141 = 0.041227165399737475 * tmp_34 + 0.78764240869137092 * tmp_35 + tmp_36;
      real_t tmp_142 = tmp_141 * tmp_38;
      real_t tmp_143 = tmp_138 + tmp_140 + tmp_142;
      real_t tmp_144 = tmp_137 * tmp_41;
      real_t tmp_145 = tmp_139 * tmp_43;
      real_t tmp_146 = tmp_141 * tmp_45;
      real_t tmp_147 = tmp_144 + tmp_145 + tmp_146;
      real_t tmp_148 = tmp_137 * tmp_48;
      real_t tmp_149 = tmp_139 * tmp_50;
      real_t tmp_150 = tmp_141 * tmp_52;
      real_t tmp_151 = tmp_148 + tmp_149 + tmp_150;
      real_t tmp_152 = tmp_1 * ( tmp_143 - 1.0 / 4.0 ) + tmp_11 * ( tmp_151 - 1.0 / 4.0 ) + tmp_8 * ( tmp_147 - 1.0 / 4.0 );
      real_t tmp_153 = tmp_152 * tmp_61;
      real_t tmp_154 = 0.019202922745021479 * tmp_63;
      real_t tmp_155 = 0.039308471900058539 * tmp_3 + 0.58463275527740355 * tmp_4 + tmp_6;
      real_t tmp_156 = tmp_155 * tmp_24;
      real_t tmp_157 = 0.039308471900058539 * tmp_27 + 0.58463275527740355 * tmp_28 + tmp_29;
      real_t tmp_158 = tmp_157 * tmp_31;
      real_t tmp_159 = 0.039308471900058539 * tmp_34 + 0.58463275527740355 * tmp_35 + tmp_36;
      real_t tmp_160 = tmp_159 * tmp_38;
      real_t tmp_161 = tmp_156 + tmp_158 + tmp_160;
      real_t tmp_162 = tmp_155 * tmp_41;
      real_t tmp_163 = tmp_157 * tmp_43;
      real_t tmp_164 = tmp_159 * tmp_45;
      real_t tmp_165 = tmp_162 + tmp_163 + tmp_164;
      real_t tmp_166 = tmp_155 * tmp_48;
      real_t tmp_167 = tmp_157 * tmp_50;
      real_t tmp_168 = tmp_159 * tmp_52;
      real_t tmp_169 = tmp_166 + tmp_167 + tmp_168;
      real_t tmp_170 = tmp_1 * ( tmp_161 - 1.0 / 4.0 ) + tmp_11 * ( tmp_169 - 1.0 / 4.0 ) + tmp_8 * ( tmp_165 - 1.0 / 4.0 );
      real_t tmp_171 = tmp_170 * tmp_61;
      real_t tmp_172 = 0.020848748529055869 * tmp_63;
      real_t tmp_173 = 0.78764240869137092 * tmp_3 + 0.041227165399737475 * tmp_4 + tmp_6;
      real_t tmp_174 = tmp_173 * tmp_24;
      real_t tmp_175 = 0.78764240869137092 * tmp_27 + 0.041227165399737475 * tmp_28 + tmp_29;
      real_t tmp_176 = tmp_175 * tmp_31;
      real_t tmp_177 = 0.78764240869137092 * tmp_34 + 0.041227165399737475 * tmp_35 + tmp_36;
      real_t tmp_178 = tmp_177 * tmp_38;
      real_t tmp_179 = tmp_174 + tmp_176 + tmp_178;
      real_t tmp_180 = tmp_173 * tmp_41;
      real_t tmp_181 = tmp_175 * tmp_43;
      real_t tmp_182 = tmp_177 * tmp_45;
      real_t tmp_183 = tmp_180 + tmp_181 + tmp_182;
      real_t tmp_184 = tmp_173 * tmp_48;
      real_t tmp_185 = tmp_175 * tmp_50;
      real_t tmp_186 = tmp_177 * tmp_52;
      real_t tmp_187 = tmp_184 + tmp_185 + tmp_186;
      real_t tmp_188 = tmp_1 * ( tmp_179 - 1.0 / 4.0 ) + tmp_11 * ( tmp_187 - 1.0 / 4.0 ) + tmp_8 * ( tmp_183 - 1.0 / 4.0 );
      real_t tmp_189 = tmp_188 * tmp_61;
      real_t tmp_190 = 0.019202922745021479 * tmp_63;
      real_t tmp_191 = 0.58463275527740355 * tmp_3 + 0.039308471900058539 * tmp_4 + tmp_6;
      real_t tmp_192 = tmp_191 * tmp_24;
      real_t tmp_193 = 0.58463275527740355 * tmp_27 + 0.039308471900058539 * tmp_28 + tmp_29;
      real_t tmp_194 = tmp_193 * tmp_31;
      real_t tmp_195 = 0.58463275527740355 * tmp_34 + 0.039308471900058539 * tmp_35 + tmp_36;
      real_t tmp_196 = tmp_195 * tmp_38;
      real_t tmp_197 = tmp_192 + tmp_194 + tmp_196;
      real_t tmp_198 = tmp_191 * tmp_41;
      real_t tmp_199 = tmp_193 * tmp_43;
      real_t tmp_200 = tmp_195 * tmp_45;
      real_t tmp_201 = tmp_198 + tmp_199 + tmp_200;
      real_t tmp_202 = tmp_191 * tmp_48;
      real_t tmp_203 = tmp_193 * tmp_50;
      real_t tmp_204 = tmp_195 * tmp_52;
      real_t tmp_205 = tmp_202 + tmp_203 + tmp_204;
      real_t tmp_206 = tmp_1 * ( tmp_197 - 1.0 / 4.0 ) + tmp_11 * ( tmp_205 - 1.0 / 4.0 ) + tmp_8 * ( tmp_201 - 1.0 / 4.0 );
      real_t tmp_207 = tmp_206 * tmp_61;
      real_t tmp_208 = 0.020848748529055869 * tmp_63;
      real_t tmp_209 = 0.1711304259088916 * tmp_3 + 0.78764240869137092 * tmp_4 + tmp_6;
      real_t tmp_210 = tmp_209 * tmp_24;
      real_t tmp_211 = 0.1711304259088916 * tmp_27 + 0.78764240869137092 * tmp_28 + tmp_29;
      real_t tmp_212 = tmp_211 * tmp_31;
      real_t tmp_213 = 0.1711304259088916 * tmp_34 + 0.78764240869137092 * tmp_35 + tmp_36;
      real_t tmp_214 = tmp_213 * tmp_38;
      real_t tmp_215 = tmp_210 + tmp_212 + tmp_214;
      real_t tmp_216 = tmp_209 * tmp_41;
      real_t tmp_217 = tmp_211 * tmp_43;
      real_t tmp_218 = tmp_213 * tmp_45;
      real_t tmp_219 = tmp_216 + tmp_217 + tmp_218;
      real_t tmp_220 = tmp_209 * tmp_48;
      real_t tmp_221 = tmp_211 * tmp_50;
      real_t tmp_222 = tmp_213 * tmp_52;
      real_t tmp_223 = tmp_220 + tmp_221 + tmp_222;
      real_t tmp_224 = tmp_1 * ( tmp_215 - 1.0 / 4.0 ) + tmp_11 * ( tmp_223 - 1.0 / 4.0 ) + tmp_8 * ( tmp_219 - 1.0 / 4.0 );
      real_t tmp_225 = tmp_224 * tmp_61;
      real_t tmp_226 = 0.019202922745021479 * tmp_63;
      real_t tmp_227 = 0.37605877282253791 * tmp_3 + 0.58463275527740355 * tmp_4 + tmp_6;
      real_t tmp_228 = tmp_227 * tmp_24;
      real_t tmp_229 = 0.37605877282253791 * tmp_27 + 0.58463275527740355 * tmp_28 + tmp_29;
      real_t tmp_230 = tmp_229 * tmp_31;
      real_t tmp_231 = 0.37605877282253791 * tmp_34 + 0.58463275527740355 * tmp_35 + tmp_36;
      real_t tmp_232 = tmp_231 * tmp_38;
      real_t tmp_233 = tmp_228 + tmp_230 + tmp_232;
      real_t tmp_234 = tmp_227 * tmp_41;
      real_t tmp_235 = tmp_229 * tmp_43;
      real_t tmp_236 = tmp_231 * tmp_45;
      real_t tmp_237 = tmp_234 + tmp_235 + tmp_236;
      real_t tmp_238 = tmp_227 * tmp_48;
      real_t tmp_239 = tmp_229 * tmp_50;
      real_t tmp_240 = tmp_231 * tmp_52;
      real_t tmp_241 = tmp_238 + tmp_239 + tmp_240;
      real_t tmp_242 = tmp_1 * ( tmp_233 - 1.0 / 4.0 ) + tmp_11 * ( tmp_241 - 1.0 / 4.0 ) + tmp_8 * ( tmp_237 - 1.0 / 4.0 );
      real_t tmp_243 = tmp_242 * tmp_61;
      real_t tmp_244 = 0.020848748529055869 * tmp_63;
      real_t tmp_245 = 0.041227165399737475 * tmp_3 + 0.1711304259088916 * tmp_4 + tmp_6;
      real_t tmp_246 = tmp_24 * tmp_245;
      real_t tmp_247 = 0.041227165399737475 * tmp_27 + 0.1711304259088916 * tmp_28 + tmp_29;
      real_t tmp_248 = tmp_247 * tmp_31;
      real_t tmp_249 = 0.041227165399737475 * tmp_34 + 0.1711304259088916 * tmp_35 + tmp_36;
      real_t tmp_250 = tmp_249 * tmp_38;
      real_t tmp_251 = tmp_246 + tmp_248 + tmp_250;
      real_t tmp_252 = tmp_245 * tmp_41;
      real_t tmp_253 = tmp_247 * tmp_43;
      real_t tmp_254 = tmp_249 * tmp_45;
      real_t tmp_255 = tmp_252 + tmp_253 + tmp_254;
      real_t tmp_256 = tmp_245 * tmp_48;
      real_t tmp_257 = tmp_247 * tmp_50;
      real_t tmp_258 = tmp_249 * tmp_52;
      real_t tmp_259 = tmp_256 + tmp_257 + tmp_258;
      real_t tmp_260 = tmp_1 * ( tmp_251 - 1.0 / 4.0 ) + tmp_11 * ( tmp_259 - 1.0 / 4.0 ) + tmp_8 * ( tmp_255 - 1.0 / 4.0 );
      real_t tmp_261 = tmp_260 * tmp_61;
      real_t tmp_262 = 0.019202922745021479 * tmp_63;
      real_t tmp_263 = 0.40446199974765351 * tmp_3 + 0.19107600050469298 * tmp_4 + tmp_6;
      real_t tmp_264 = tmp_24 * tmp_263;
      real_t tmp_265 = 0.40446199974765351 * tmp_27 + 0.19107600050469298 * tmp_28 + tmp_29;
      real_t tmp_266 = tmp_265 * tmp_31;
      real_t tmp_267 = 0.40446199974765351 * tmp_34 + 0.19107600050469298 * tmp_35 + tmp_36;
      real_t tmp_268 = tmp_267 * tmp_38;
      real_t tmp_269 = tmp_264 + tmp_266 + tmp_268;
      real_t tmp_270 = tmp_263 * tmp_41;
      real_t tmp_271 = tmp_265 * tmp_43;
      real_t tmp_272 = tmp_267 * tmp_45;
      real_t tmp_273 = tmp_270 + tmp_271 + tmp_272;
      real_t tmp_274 = tmp_263 * tmp_48;
      real_t tmp_275 = tmp_265 * tmp_50;
      real_t tmp_276 = tmp_267 * tmp_52;
      real_t tmp_277 = tmp_274 + tmp_275 + tmp_276;
      real_t tmp_278 = tmp_1 * ( tmp_269 - 1.0 / 4.0 ) + tmp_11 * ( tmp_277 - 1.0 / 4.0 ) + tmp_8 * ( tmp_273 - 1.0 / 4.0 );
      real_t tmp_279 = tmp_278 * tmp_61;
      real_t tmp_280 = 0.042507265838595799 * tmp_63;
      real_t tmp_281 = 0.039308471900058539 * tmp_3 + 0.37605877282253791 * tmp_4 + tmp_6;
      real_t tmp_282 = tmp_24 * tmp_281;
      real_t tmp_283 = 0.039308471900058539 * tmp_27 + 0.37605877282253791 * tmp_28 + tmp_29;
      real_t tmp_284 = tmp_283 * tmp_31;
      real_t tmp_285 = 0.039308471900058539 * tmp_34 + 0.37605877282253791 * tmp_35 + tmp_36;
      real_t tmp_286 = tmp_285 * tmp_38;
      real_t tmp_287 = tmp_282 + tmp_284 + tmp_286;
      real_t tmp_288 = tmp_281 * tmp_41;
      real_t tmp_289 = tmp_283 * tmp_43;
      real_t tmp_290 = tmp_285 * tmp_45;
      real_t tmp_291 = tmp_288 + tmp_289 + tmp_290;
      real_t tmp_292 = tmp_281 * tmp_48;
      real_t tmp_293 = tmp_283 * tmp_50;
      real_t tmp_294 = tmp_285 * tmp_52;
      real_t tmp_295 = tmp_292 + tmp_293 + tmp_294;
      real_t tmp_296 = tmp_1 * ( tmp_287 - 1.0 / 4.0 ) + tmp_11 * ( tmp_295 - 1.0 / 4.0 ) + tmp_8 * ( tmp_291 - 1.0 / 4.0 );
      real_t tmp_297 = tmp_296 * tmp_61;
      real_t tmp_298 = 0.020848748529055869 * tmp_63;
      real_t tmp_299 = 0.93718850182767688 * tmp_3 + 0.031405749086161582 * tmp_4 + tmp_6;
      real_t tmp_300 = tmp_24 * tmp_299;
      real_t tmp_301 = 0.93718850182767688 * tmp_27 + 0.031405749086161582 * tmp_28 + tmp_29;
      real_t tmp_302 = tmp_301 * tmp_31;
      real_t tmp_303 = 0.93718850182767688 * tmp_34 + 0.031405749086161582 * tmp_35 + tmp_36;
      real_t tmp_304 = tmp_303 * tmp_38;
      real_t tmp_305 = tmp_300 + tmp_302 + tmp_304;
      real_t tmp_306 = tmp_299 * tmp_41;
      real_t tmp_307 = tmp_301 * tmp_43;
      real_t tmp_308 = tmp_303 * tmp_45;
      real_t tmp_309 = tmp_306 + tmp_307 + tmp_308;
      real_t tmp_310 = tmp_299 * tmp_48;
      real_t tmp_311 = tmp_301 * tmp_50;
      real_t tmp_312 = tmp_303 * tmp_52;
      real_t tmp_313 = tmp_310 + tmp_311 + tmp_312;
      real_t tmp_314 = tmp_1 * ( tmp_305 - 1.0 / 4.0 ) + tmp_11 * ( tmp_313 - 1.0 / 4.0 ) + tmp_8 * ( tmp_309 - 1.0 / 4.0 );
      real_t tmp_315 = tmp_314 * tmp_61;
      real_t tmp_316 = 0.0068572537431980923 * tmp_63;
      real_t tmp_317 = 0.60796128279561268 * tmp_3 + 0.19601935860219369 * tmp_4 + tmp_6;
      real_t tmp_318 = tmp_24 * tmp_317;
      real_t tmp_319 = 0.60796128279561268 * tmp_27 + 0.19601935860219369 * tmp_28 + tmp_29;
      real_t tmp_320 = tmp_31 * tmp_319;
      real_t tmp_321 = 0.60796128279561268 * tmp_34 + 0.19601935860219369 * tmp_35 + tmp_36;
      real_t tmp_322 = tmp_321 * tmp_38;
      real_t tmp_323 = tmp_318 + tmp_320 + tmp_322;
      real_t tmp_324 = tmp_317 * tmp_41;
      real_t tmp_325 = tmp_319 * tmp_43;
      real_t tmp_326 = tmp_321 * tmp_45;
      real_t tmp_327 = tmp_324 + tmp_325 + tmp_326;
      real_t tmp_328 = tmp_317 * tmp_48;
      real_t tmp_329 = tmp_319 * tmp_50;
      real_t tmp_330 = tmp_321 * tmp_52;
      real_t tmp_331 = tmp_328 + tmp_329 + tmp_330;
      real_t tmp_332 = tmp_1 * ( tmp_323 - 1.0 / 4.0 ) + tmp_11 * ( tmp_331 - 1.0 / 4.0 ) + tmp_8 * ( tmp_327 - 1.0 / 4.0 );
      real_t tmp_333 = tmp_332 * tmp_61;
      real_t tmp_334 = 0.037198804536718075 * tmp_63;
      real_t tmp_335 = 0.19107600050469298 * tmp_3 + 0.40446199974765351 * tmp_4 + tmp_6;
      real_t tmp_336 = tmp_24 * tmp_335;
      real_t tmp_337 = 0.19107600050469298 * tmp_27 + 0.40446199974765351 * tmp_28 + tmp_29;
      real_t tmp_338 = tmp_31 * tmp_337;
      real_t tmp_339 = 0.19107600050469298 * tmp_34 + 0.40446199974765351 * tmp_35 + tmp_36;
      real_t tmp_340 = tmp_339 * tmp_38;
      real_t tmp_341 = tmp_336 + tmp_338 + tmp_340;
      real_t tmp_342 = tmp_335 * tmp_41;
      real_t tmp_343 = tmp_337 * tmp_43;
      real_t tmp_344 = tmp_339 * tmp_45;
      real_t tmp_345 = tmp_342 + tmp_343 + tmp_344;
      real_t tmp_346 = tmp_335 * tmp_48;
      real_t tmp_347 = tmp_337 * tmp_50;
      real_t tmp_348 = tmp_339 * tmp_52;
      real_t tmp_349 = tmp_346 + tmp_347 + tmp_348;
      real_t tmp_350 = tmp_1 * ( tmp_341 - 1.0 / 4.0 ) + tmp_11 * ( tmp_349 - 1.0 / 4.0 ) + tmp_8 * ( tmp_345 - 1.0 / 4.0 );
      real_t tmp_351 = tmp_350 * tmp_61;
      real_t tmp_352 = 0.042507265838595799 * tmp_63;
      real_t tmp_353 = 0.031405749086161582 * tmp_3 + 0.031405749086161582 * tmp_4 + tmp_6;
      real_t tmp_354 = tmp_24 * tmp_353;
      real_t tmp_355 = 0.031405749086161582 * tmp_27 + 0.031405749086161582 * tmp_28 + tmp_29;
      real_t tmp_356 = tmp_31 * tmp_355;
      real_t tmp_357 = 0.031405749086161582 * tmp_34 + 0.031405749086161582 * tmp_35 + tmp_36;
      real_t tmp_358 = tmp_357 * tmp_38;
      real_t tmp_359 = tmp_354 + tmp_356 + tmp_358;
      real_t tmp_360 = tmp_353 * tmp_41;
      real_t tmp_361 = tmp_355 * tmp_43;
      real_t tmp_362 = tmp_357 * tmp_45;
      real_t tmp_363 = tmp_360 + tmp_361 + tmp_362;
      real_t tmp_364 = tmp_353 * tmp_48;
      real_t tmp_365 = tmp_355 * tmp_50;
      real_t tmp_366 = tmp_357 * tmp_52;
      real_t tmp_367 = tmp_364 + tmp_365 + tmp_366;
      real_t tmp_368 = tmp_1 * ( tmp_359 - 1.0 / 4.0 ) + tmp_11 * ( tmp_367 - 1.0 / 4.0 ) + tmp_8 * ( tmp_363 - 1.0 / 4.0 );
      real_t tmp_369 = tmp_368 * tmp_61;
      real_t tmp_370 = 0.0068572537431980923 * tmp_63;
      real_t tmp_371 = 0.19601935860219369 * tmp_3 + 0.19601935860219369 * tmp_4 + tmp_6;
      real_t tmp_372 = tmp_24 * tmp_371;
      real_t tmp_373 = 0.19601935860219369 * tmp_27 + 0.19601935860219369 * tmp_28 + tmp_29;
      real_t tmp_374 = tmp_31 * tmp_373;
      real_t tmp_375 = 0.19601935860219369 * tmp_34 + 0.19601935860219369 * tmp_35 + tmp_36;
      real_t tmp_376 = tmp_375 * tmp_38;
      real_t tmp_377 = tmp_372 + tmp_374 + tmp_376;
      real_t tmp_378 = tmp_371 * tmp_41;
      real_t tmp_379 = tmp_373 * tmp_43;
      real_t tmp_380 = tmp_375 * tmp_45;
      real_t tmp_381 = tmp_378 + tmp_379 + tmp_380;
      real_t tmp_382 = tmp_371 * tmp_48;
      real_t tmp_383 = tmp_373 * tmp_50;
      real_t tmp_384 = tmp_375 * tmp_52;
      real_t tmp_385 = tmp_382 + tmp_383 + tmp_384;
      real_t tmp_386 = tmp_1 * ( tmp_377 - 1.0 / 4.0 ) + tmp_11 * ( tmp_385 - 1.0 / 4.0 ) + tmp_8 * ( tmp_381 - 1.0 / 4.0 );
      real_t tmp_387 = tmp_386 * tmp_61;
      real_t tmp_388 = 0.037198804536718075 * tmp_63;
      real_t tmp_389 = 0.40446199974765351 * tmp_3 + 0.40446199974765351 * tmp_4 + tmp_6;
      real_t tmp_390 = tmp_24 * tmp_389;
      real_t tmp_391 = 0.40446199974765351 * tmp_27 + 0.40446199974765351 * tmp_28 + tmp_29;
      real_t tmp_392 = tmp_31 * tmp_391;
      real_t tmp_393 = 0.40446199974765351 * tmp_34 + 0.40446199974765351 * tmp_35 + tmp_36;
      real_t tmp_394 = tmp_38 * tmp_393;
      real_t tmp_395 = tmp_390 + tmp_392 + tmp_394;
      real_t tmp_396 = tmp_389 * tmp_41;
      real_t tmp_397 = tmp_391 * tmp_43;
      real_t tmp_398 = tmp_393 * tmp_45;
      real_t tmp_399 = tmp_396 + tmp_397 + tmp_398;
      real_t tmp_400 = tmp_389 * tmp_48;
      real_t tmp_401 = tmp_391 * tmp_50;
      real_t tmp_402 = tmp_393 * tmp_52;
      real_t tmp_403 = tmp_400 + tmp_401 + tmp_402;
      real_t tmp_404 = tmp_1 * ( tmp_395 - 1.0 / 4.0 ) + tmp_11 * ( tmp_403 - 1.0 / 4.0 ) + tmp_8 * ( tmp_399 - 1.0 / 4.0 );
      real_t tmp_405 = tmp_404 * tmp_61;
      real_t tmp_406 = 0.042507265838595799 * tmp_63;
      real_t tmp_407 = 0.1711304259088916 * tmp_3 + 0.041227165399737475 * tmp_4 + tmp_6;
      real_t tmp_408 = tmp_24 * tmp_407;
      real_t tmp_409 = 0.1711304259088916 * tmp_27 + 0.041227165399737475 * tmp_28 + tmp_29;
      real_t tmp_410 = tmp_31 * tmp_409;
      real_t tmp_411 = 0.1711304259088916 * tmp_34 + 0.041227165399737475 * tmp_35 + tmp_36;
      real_t tmp_412 = tmp_38 * tmp_411;
      real_t tmp_413 = tmp_408 + tmp_410 + tmp_412;
      real_t tmp_414 = tmp_407 * tmp_41;
      real_t tmp_415 = tmp_409 * tmp_43;
      real_t tmp_416 = tmp_411 * tmp_45;
      real_t tmp_417 = tmp_414 + tmp_415 + tmp_416;
      real_t tmp_418 = tmp_407 * tmp_48;
      real_t tmp_419 = tmp_409 * tmp_50;
      real_t tmp_420 = tmp_411 * tmp_52;
      real_t tmp_421 = tmp_418 + tmp_419 + tmp_420;
      real_t tmp_422 = tmp_1 * ( tmp_413 - 1.0 / 4.0 ) + tmp_11 * ( tmp_421 - 1.0 / 4.0 ) + tmp_8 * ( tmp_417 - 1.0 / 4.0 );
      real_t tmp_423 = tmp_422 * tmp_61;
      real_t tmp_424 = 0.019202922745021479 * tmp_63;
      real_t tmp_425 = 0.5 * p_affine_13_0 * tmp_38 + 0.5 * p_affine_13_1 * tmp_31 + 0.5 * p_affine_13_2 * tmp_24;
      real_t tmp_426 = 0.5 * p_affine_13_0 * tmp_45 + 0.5 * p_affine_13_1 * tmp_43 + 0.5 * p_affine_13_2 * tmp_41;
      real_t tmp_427 = 0.5 * p_affine_13_0 * tmp_52 + 0.5 * p_affine_13_1 * tmp_50 + 0.5 * p_affine_13_2 * tmp_48;
      real_t a_0_0 =
          tmp_100 * ( -tmp_56 * tmp_98 +
                      tmp_99 * ( -tmp_84 - tmp_86 - tmp_88 - tmp_90 - tmp_91 - tmp_92 - tmp_94 - tmp_95 - tmp_96 + 1 ) ) +
          tmp_118 * ( -tmp_116 * tmp_56 + tmp_117 * ( -tmp_102 - tmp_104 - tmp_106 - tmp_108 - tmp_109 - tmp_110 - tmp_112 -
                                                      tmp_113 - tmp_114 + 1 ) ) +
          tmp_136 * ( -tmp_134 * tmp_56 + tmp_135 * ( -tmp_120 - tmp_122 - tmp_124 - tmp_126 - tmp_127 - tmp_128 - tmp_130 -
                                                      tmp_131 - tmp_132 + 1 ) ) +
          tmp_154 * ( -tmp_152 * tmp_56 + tmp_153 * ( -tmp_138 - tmp_140 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_148 -
                                                      tmp_149 - tmp_150 + 1 ) ) +
          tmp_172 * ( -tmp_170 * tmp_56 + tmp_171 * ( -tmp_156 - tmp_158 - tmp_160 - tmp_162 - tmp_163 - tmp_164 - tmp_166 -
                                                      tmp_167 - tmp_168 + 1 ) ) +
          tmp_190 * ( -tmp_188 * tmp_56 + tmp_189 * ( -tmp_174 - tmp_176 - tmp_178 - tmp_180 - tmp_181 - tmp_182 - tmp_184 -
                                                      tmp_185 - tmp_186 + 1 ) ) +
          tmp_208 * ( -tmp_206 * tmp_56 + tmp_207 * ( -tmp_192 - tmp_194 - tmp_196 - tmp_198 - tmp_199 - tmp_200 - tmp_202 -
                                                      tmp_203 - tmp_204 + 1 ) ) +
          tmp_226 * ( -tmp_224 * tmp_56 + tmp_225 * ( -tmp_210 - tmp_212 - tmp_214 - tmp_216 - tmp_217 - tmp_218 - tmp_220 -
                                                      tmp_221 - tmp_222 + 1 ) ) +
          tmp_244 * ( -tmp_242 * tmp_56 + tmp_243 * ( -tmp_228 - tmp_230 - tmp_232 - tmp_234 - tmp_235 - tmp_236 - tmp_238 -
                                                      tmp_239 - tmp_240 + 1 ) ) +
          tmp_262 * ( -tmp_260 * tmp_56 + tmp_261 * ( -tmp_246 - tmp_248 - tmp_250 - tmp_252 - tmp_253 - tmp_254 - tmp_256 -
                                                      tmp_257 - tmp_258 + 1 ) ) +
          tmp_280 * ( -tmp_278 * tmp_56 + tmp_279 * ( -tmp_264 - tmp_266 - tmp_268 - tmp_270 - tmp_271 - tmp_272 - tmp_274 -
                                                      tmp_275 - tmp_276 + 1 ) ) +
          tmp_298 * ( -tmp_296 * tmp_56 + tmp_297 * ( -tmp_282 - tmp_284 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_292 -
                                                      tmp_293 - tmp_294 + 1 ) ) +
          tmp_316 * ( -tmp_314 * tmp_56 + tmp_315 * ( -tmp_300 - tmp_302 - tmp_304 - tmp_306 - tmp_307 - tmp_308 - tmp_310 -
                                                      tmp_311 - tmp_312 + 1 ) ) +
          tmp_334 * ( -tmp_332 * tmp_56 + tmp_333 * ( -tmp_318 - tmp_320 - tmp_322 - tmp_324 - tmp_325 - tmp_326 - tmp_328 -
                                                      tmp_329 - tmp_330 + 1 ) ) +
          tmp_352 * ( -tmp_350 * tmp_56 + tmp_351 * ( -tmp_336 - tmp_338 - tmp_340 - tmp_342 - tmp_343 - tmp_344 - tmp_346 -
                                                      tmp_347 - tmp_348 + 1 ) ) +
          tmp_370 * ( -tmp_368 * tmp_56 + tmp_369 * ( -tmp_354 - tmp_356 - tmp_358 - tmp_360 - tmp_361 - tmp_362 - tmp_364 -
                                                      tmp_365 - tmp_366 + 1 ) ) +
          tmp_388 * ( -tmp_386 * tmp_56 + tmp_387 * ( -tmp_372 - tmp_374 - tmp_376 - tmp_378 - tmp_379 - tmp_380 - tmp_382 -
                                                      tmp_383 - tmp_384 + 1 ) ) +
          tmp_406 * ( -tmp_404 * tmp_56 + tmp_405 * ( -tmp_390 - tmp_392 - tmp_394 - tmp_396 - tmp_397 - tmp_398 - tmp_400 -
                                                      tmp_401 - tmp_402 + 1 ) ) +
          tmp_424 * ( -tmp_422 * tmp_56 + tmp_423 * ( -tmp_408 - tmp_410 - tmp_412 - tmp_414 - tmp_415 - tmp_416 - tmp_418 -
                                                      tmp_419 - tmp_420 + 1 ) ) +
          tmp_64 * ( -tmp_55 * tmp_56 +
                     tmp_62 * ( -tmp_25 - tmp_32 - tmp_39 - tmp_42 - tmp_44 - tmp_46 - tmp_49 - tmp_51 - tmp_53 + 1 ) ) +
          tmp_82 * ( -tmp_56 * tmp_80 +
                     tmp_81 * ( -tmp_66 - tmp_68 - tmp_70 - tmp_72 - tmp_73 - tmp_74 - tmp_76 - tmp_77 - tmp_78 + 1 ) );
      real_t a_0_1 = tmp_100 * ( -tmp_425 * tmp_98 + tmp_89 * tmp_99 ) + tmp_118 * ( tmp_107 * tmp_117 - tmp_116 * tmp_425 ) +
                     tmp_136 * ( tmp_125 * tmp_135 - tmp_134 * tmp_425 ) + tmp_154 * ( tmp_143 * tmp_153 - tmp_152 * tmp_425 ) +
                     tmp_172 * ( tmp_161 * tmp_171 - tmp_170 * tmp_425 ) + tmp_190 * ( tmp_179 * tmp_189 - tmp_188 * tmp_425 ) +
                     tmp_208 * ( tmp_197 * tmp_207 - tmp_206 * tmp_425 ) + tmp_226 * ( tmp_215 * tmp_225 - tmp_224 * tmp_425 ) +
                     tmp_244 * ( tmp_233 * tmp_243 - tmp_242 * tmp_425 ) + tmp_262 * ( tmp_251 * tmp_261 - tmp_260 * tmp_425 ) +
                     tmp_280 * ( tmp_269 * tmp_279 - tmp_278 * tmp_425 ) + tmp_298 * ( tmp_287 * tmp_297 - tmp_296 * tmp_425 ) +
                     tmp_316 * ( tmp_305 * tmp_315 - tmp_314 * tmp_425 ) + tmp_334 * ( tmp_323 * tmp_333 - tmp_332 * tmp_425 ) +
                     tmp_352 * ( tmp_341 * tmp_351 - tmp_350 * tmp_425 ) + tmp_370 * ( tmp_359 * tmp_369 - tmp_368 * tmp_425 ) +
                     tmp_388 * ( tmp_377 * tmp_387 - tmp_386 * tmp_425 ) + tmp_406 * ( tmp_395 * tmp_405 - tmp_404 * tmp_425 ) +
                     tmp_424 * ( tmp_413 * tmp_423 - tmp_422 * tmp_425 ) + tmp_64 * ( tmp_40 * tmp_62 - tmp_425 * tmp_55 ) +
                     tmp_82 * ( -tmp_425 * tmp_80 + tmp_71 * tmp_81 );
      real_t a_0_2 = tmp_100 * ( -tmp_426 * tmp_98 + tmp_93 * tmp_99 ) + tmp_118 * ( tmp_111 * tmp_117 - tmp_116 * tmp_426 ) +
                     tmp_136 * ( tmp_129 * tmp_135 - tmp_134 * tmp_426 ) + tmp_154 * ( tmp_147 * tmp_153 - tmp_152 * tmp_426 ) +
                     tmp_172 * ( tmp_165 * tmp_171 - tmp_170 * tmp_426 ) + tmp_190 * ( tmp_183 * tmp_189 - tmp_188 * tmp_426 ) +
                     tmp_208 * ( tmp_201 * tmp_207 - tmp_206 * tmp_426 ) + tmp_226 * ( tmp_219 * tmp_225 - tmp_224 * tmp_426 ) +
                     tmp_244 * ( tmp_237 * tmp_243 - tmp_242 * tmp_426 ) + tmp_262 * ( tmp_255 * tmp_261 - tmp_260 * tmp_426 ) +
                     tmp_280 * ( tmp_273 * tmp_279 - tmp_278 * tmp_426 ) + tmp_298 * ( tmp_291 * tmp_297 - tmp_296 * tmp_426 ) +
                     tmp_316 * ( tmp_309 * tmp_315 - tmp_314 * tmp_426 ) + tmp_334 * ( tmp_327 * tmp_333 - tmp_332 * tmp_426 ) +
                     tmp_352 * ( tmp_345 * tmp_351 - tmp_350 * tmp_426 ) + tmp_370 * ( tmp_363 * tmp_369 - tmp_368 * tmp_426 ) +
                     tmp_388 * ( tmp_381 * tmp_387 - tmp_386 * tmp_426 ) + tmp_406 * ( tmp_399 * tmp_405 - tmp_404 * tmp_426 ) +
                     tmp_424 * ( tmp_417 * tmp_423 - tmp_422 * tmp_426 ) + tmp_64 * ( -tmp_426 * tmp_55 + tmp_47 * tmp_62 ) +
                     tmp_82 * ( -tmp_426 * tmp_80 + tmp_75 * tmp_81 );
      real_t a_0_3 = tmp_100 * ( -tmp_427 * tmp_98 + tmp_97 * tmp_99 ) + tmp_118 * ( tmp_115 * tmp_117 - tmp_116 * tmp_427 ) +
                     tmp_136 * ( tmp_133 * tmp_135 - tmp_134 * tmp_427 ) + tmp_154 * ( tmp_151 * tmp_153 - tmp_152 * tmp_427 ) +
                     tmp_172 * ( tmp_169 * tmp_171 - tmp_170 * tmp_427 ) + tmp_190 * ( tmp_187 * tmp_189 - tmp_188 * tmp_427 ) +
                     tmp_208 * ( tmp_205 * tmp_207 - tmp_206 * tmp_427 ) + tmp_226 * ( tmp_223 * tmp_225 - tmp_224 * tmp_427 ) +
                     tmp_244 * ( tmp_241 * tmp_243 - tmp_242 * tmp_427 ) + tmp_262 * ( tmp_259 * tmp_261 - tmp_260 * tmp_427 ) +
                     tmp_280 * ( tmp_277 * tmp_279 - tmp_278 * tmp_427 ) + tmp_298 * ( tmp_295 * tmp_297 - tmp_296 * tmp_427 ) +
                     tmp_316 * ( tmp_313 * tmp_315 - tmp_314 * tmp_427 ) + tmp_334 * ( tmp_331 * tmp_333 - tmp_332 * tmp_427 ) +
                     tmp_352 * ( tmp_349 * tmp_351 - tmp_350 * tmp_427 ) + tmp_370 * ( tmp_367 * tmp_369 - tmp_368 * tmp_427 ) +
                     tmp_388 * ( tmp_385 * tmp_387 - tmp_386 * tmp_427 ) + tmp_406 * ( tmp_403 * tmp_405 - tmp_404 * tmp_427 ) +
                     tmp_424 * ( tmp_421 * tmp_423 - tmp_422 * tmp_427 ) + tmp_64 * ( -tmp_427 * tmp_55 + tmp_54 * tmp_62 ) +
                     tmp_82 * ( -tmp_427 * tmp_80 + tmp_79 * tmp_81 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 0, 3 ) = a_0_3;
   }

   void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                  const std::vector< Point3D >& coordsElementOuter,
                                  const std::vector< Point3D >& coordsFacet,
                                  const Point3D&,
                                  const Point3D&,
                                  const Point3D&     outwardNormal,
                                  const DGBasisInfo& trialBasis,
                                  const DGBasisInfo& testBasis,
                                  int                trialDegree,
                                  int                testDegree,
                                  MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = p_affine_2_0 + tmp_0;
      real_t tmp_3  = -p_affine_0_1;
      real_t tmp_4  = p_affine_3_1 + tmp_3;
      real_t tmp_5  = tmp_2 * tmp_4;
      real_t tmp_6  = p_affine_3_0 + tmp_0;
      real_t tmp_7  = p_affine_2_1 + tmp_3;
      real_t tmp_8  = tmp_6 * tmp_7;
      real_t tmp_9  = tmp_5 - tmp_8;
      real_t tmp_10 = -p_affine_0_2;
      real_t tmp_11 = p_affine_3_2 + tmp_10;
      real_t tmp_12 = tmp_11 * tmp_7;
      real_t tmp_13 = p_affine_1_2 + tmp_10;
      real_t tmp_14 = p_affine_1_1 + tmp_3;
      real_t tmp_15 = p_affine_2_2 + tmp_10;
      real_t tmp_16 = tmp_15 * tmp_6;
      real_t tmp_17 = tmp_15 * tmp_4;
      real_t tmp_18 = tmp_11 * tmp_2;
      real_t tmp_19 =
          1.0 / ( tmp_1 * tmp_12 - tmp_1 * tmp_17 + tmp_13 * tmp_5 - tmp_13 * tmp_8 + tmp_14 * tmp_16 - tmp_14 * tmp_18 );
      real_t tmp_20 = p_affine_8_2 + tmp_10;
      real_t tmp_21 = -p_affine_8_2;
      real_t tmp_22 = p_affine_9_2 + tmp_21;
      real_t tmp_23 = p_affine_10_2 + tmp_21;
      real_t tmp_24 = 0.031405749086161582 * tmp_22 + 0.93718850182767688 * tmp_23;
      real_t tmp_25 = tmp_19 * ( tmp_20 + tmp_24 );
      real_t tmp_26 = tmp_16 - tmp_18;
      real_t tmp_27 = p_affine_8_1 + tmp_3;
      real_t tmp_28 = -p_affine_8_1;
      real_t tmp_29 = p_affine_9_1 + tmp_28;
      real_t tmp_30 = p_affine_10_1 + tmp_28;
      real_t tmp_31 = 0.031405749086161582 * tmp_29 + 0.93718850182767688 * tmp_30;
      real_t tmp_32 = tmp_19 * ( tmp_27 + tmp_31 );
      real_t tmp_33 = tmp_12 - tmp_17;
      real_t tmp_34 = p_affine_8_0 + tmp_0;
      real_t tmp_35 = -p_affine_8_0;
      real_t tmp_36 = p_affine_9_0 + tmp_35;
      real_t tmp_37 = p_affine_10_0 + tmp_35;
      real_t tmp_38 = 0.031405749086161582 * tmp_36 + 0.93718850182767688 * tmp_37;
      real_t tmp_39 = tmp_19 * ( tmp_34 + tmp_38 );
      real_t tmp_40 = -tmp_1 * tmp_4 + tmp_14 * tmp_6;
      real_t tmp_41 = tmp_1 * tmp_11 - tmp_13 * tmp_6;
      real_t tmp_42 = -tmp_11 * tmp_14 + tmp_13 * tmp_4;
      real_t tmp_43 = tmp_1 * tmp_7 - tmp_14 * tmp_2;
      real_t tmp_44 = -tmp_1 * tmp_15 + tmp_13 * tmp_2;
      real_t tmp_45 = -tmp_13 * tmp_7 + tmp_14 * tmp_15;
      real_t tmp_46 = tmp_1 * ( tmp_25 * tmp_9 + tmp_26 * tmp_32 + tmp_33 * tmp_39 - 1.0 / 4.0 ) +
                      tmp_2 * ( tmp_25 * tmp_40 + tmp_32 * tmp_41 + tmp_39 * tmp_42 - 1.0 / 4.0 ) +
                      tmp_6 * ( tmp_25 * tmp_43 + tmp_32 * tmp_44 + tmp_39 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_47 = -p_affine_4_1;
      real_t tmp_48 = p_affine_5_1 + tmp_47;
      real_t tmp_49 = -p_affine_4_2;
      real_t tmp_50 = p_affine_6_2 + tmp_49;
      real_t tmp_51 = tmp_48 * tmp_50;
      real_t tmp_52 = p_affine_6_1 + tmp_47;
      real_t tmp_53 = p_affine_5_2 + tmp_49;
      real_t tmp_54 = tmp_52 * tmp_53;
      real_t tmp_55 = -p_affine_4_0;
      real_t tmp_56 = p_affine_5_0 + tmp_55;
      real_t tmp_57 = p_affine_7_2 + tmp_49;
      real_t tmp_58 = tmp_52 * tmp_57;
      real_t tmp_59 = p_affine_6_0 + tmp_55;
      real_t tmp_60 = p_affine_7_1 + tmp_47;
      real_t tmp_61 = tmp_53 * tmp_60;
      real_t tmp_62 = p_affine_7_0 + tmp_55;
      real_t tmp_63 = tmp_50 * tmp_60;
      real_t tmp_64 = tmp_48 * tmp_57;
      real_t tmp_65 =
          1.0 / ( tmp_51 * tmp_62 - tmp_54 * tmp_62 + tmp_56 * tmp_58 - tmp_56 * tmp_63 + tmp_59 * tmp_61 - tmp_59 * tmp_64 );
      real_t tmp_66 = tmp_65 * ( tmp_51 - tmp_54 );
      real_t tmp_67 = tmp_65 * ( tmp_61 - tmp_64 );
      real_t tmp_68 = tmp_65 * ( tmp_58 - tmp_63 );
      real_t tmp_69 = tmp_65 * ( -tmp_50 * tmp_56 + tmp_53 * tmp_59 );
      real_t tmp_70 = tmp_65 * ( -tmp_53 * tmp_62 + tmp_56 * tmp_57 );
      real_t tmp_71 = tmp_65 * ( tmp_50 * tmp_62 - tmp_57 * tmp_59 );
      real_t tmp_72 = tmp_65 * ( -tmp_48 * tmp_59 + tmp_52 * tmp_56 );
      real_t tmp_73 = tmp_65 * ( tmp_48 * tmp_62 - tmp_56 * tmp_60 );
      real_t tmp_74 = tmp_65 * ( -tmp_52 * tmp_62 + tmp_59 * tmp_60 );
      real_t tmp_75 = 0.5 * p_affine_13_0 * ( -tmp_66 - tmp_67 - tmp_68 ) + 0.5 * p_affine_13_1 * ( -tmp_69 - tmp_70 - tmp_71 ) +
                      0.5 * p_affine_13_2 * ( -tmp_72 - tmp_73 - tmp_74 );
      real_t tmp_76 = p_affine_8_2 + tmp_49;
      real_t tmp_77 = tmp_24 + tmp_76;
      real_t tmp_78 = tmp_72 * tmp_77;
      real_t tmp_79 = tmp_73 * tmp_77;
      real_t tmp_80 = p_affine_8_1 + tmp_47;
      real_t tmp_81 = tmp_31 + tmp_80;
      real_t tmp_82 = tmp_69 * tmp_81;
      real_t tmp_83 = tmp_70 * tmp_81;
      real_t tmp_84 = tmp_74 * tmp_77;
      real_t tmp_85 = tmp_71 * tmp_81;
      real_t tmp_86 = p_affine_8_0 + tmp_55;
      real_t tmp_87 = tmp_38 + tmp_86;
      real_t tmp_88 = tmp_66 * tmp_87;
      real_t tmp_89 = tmp_67 * tmp_87;
      real_t tmp_90 = tmp_68 * tmp_87;
      real_t tmp_91 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_92 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_93 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_94 = ( std::abs( tmp_23 * tmp_91 - tmp_30 * tmp_93 ) * std::abs( tmp_23 * tmp_91 - tmp_30 * tmp_93 ) ) +
                      ( std::abs( tmp_23 * tmp_92 - tmp_37 * tmp_93 ) * std::abs( tmp_23 * tmp_92 - tmp_37 * tmp_93 ) ) +
                      ( std::abs( tmp_30 * tmp_92 - tmp_37 * tmp_91 ) * std::abs( tmp_30 * tmp_92 - tmp_37 * tmp_91 ) );
      real_t tmp_95  = 1.0 * std::pow( tmp_94, -0.25 );
      real_t tmp_96  = tmp_46 * tmp_95;
      real_t tmp_97  = 1.0 * std::pow( tmp_94, 1.0 / 2.0 );
      real_t tmp_98  = 0.0068572537431980923 * tmp_97;
      real_t tmp_99  = 0.19601935860219369 * tmp_22 + 0.60796128279561268 * tmp_23;
      real_t tmp_100 = tmp_19 * ( tmp_20 + tmp_99 );
      real_t tmp_101 = 0.19601935860219369 * tmp_29 + 0.60796128279561268 * tmp_30;
      real_t tmp_102 = tmp_19 * ( tmp_101 + tmp_27 );
      real_t tmp_103 = 0.19601935860219369 * tmp_36 + 0.60796128279561268 * tmp_37;
      real_t tmp_104 = tmp_19 * ( tmp_103 + tmp_34 );
      real_t tmp_105 = tmp_1 * ( tmp_100 * tmp_9 + tmp_102 * tmp_26 + tmp_104 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_100 * tmp_40 + tmp_102 * tmp_41 + tmp_104 * tmp_42 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_100 * tmp_43 + tmp_102 * tmp_44 + tmp_104 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_106 = tmp_76 + tmp_99;
      real_t tmp_107 = tmp_106 * tmp_72;
      real_t tmp_108 = tmp_106 * tmp_73;
      real_t tmp_109 = tmp_101 + tmp_80;
      real_t tmp_110 = tmp_109 * tmp_69;
      real_t tmp_111 = tmp_109 * tmp_70;
      real_t tmp_112 = tmp_106 * tmp_74;
      real_t tmp_113 = tmp_109 * tmp_71;
      real_t tmp_114 = tmp_103 + tmp_86;
      real_t tmp_115 = tmp_114 * tmp_66;
      real_t tmp_116 = tmp_114 * tmp_67;
      real_t tmp_117 = tmp_114 * tmp_68;
      real_t tmp_118 = tmp_105 * tmp_95;
      real_t tmp_119 = 0.037198804536718075 * tmp_97;
      real_t tmp_120 = 0.37605877282253791 * tmp_22 + 0.039308471900058539 * tmp_23;
      real_t tmp_121 = tmp_19 * ( tmp_120 + tmp_20 );
      real_t tmp_122 = 0.37605877282253791 * tmp_29 + 0.039308471900058539 * tmp_30;
      real_t tmp_123 = tmp_19 * ( tmp_122 + tmp_27 );
      real_t tmp_124 = 0.37605877282253791 * tmp_36 + 0.039308471900058539 * tmp_37;
      real_t tmp_125 = tmp_19 * ( tmp_124 + tmp_34 );
      real_t tmp_126 = tmp_1 * ( tmp_121 * tmp_9 + tmp_123 * tmp_26 + tmp_125 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_121 * tmp_40 + tmp_123 * tmp_41 + tmp_125 * tmp_42 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_121 * tmp_43 + tmp_123 * tmp_44 + tmp_125 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_127 = tmp_120 + tmp_76;
      real_t tmp_128 = tmp_127 * tmp_72;
      real_t tmp_129 = tmp_127 * tmp_73;
      real_t tmp_130 = tmp_122 + tmp_80;
      real_t tmp_131 = tmp_130 * tmp_69;
      real_t tmp_132 = tmp_130 * tmp_70;
      real_t tmp_133 = tmp_127 * tmp_74;
      real_t tmp_134 = tmp_130 * tmp_71;
      real_t tmp_135 = tmp_124 + tmp_86;
      real_t tmp_136 = tmp_135 * tmp_66;
      real_t tmp_137 = tmp_135 * tmp_67;
      real_t tmp_138 = tmp_135 * tmp_68;
      real_t tmp_139 = tmp_126 * tmp_95;
      real_t tmp_140 = 0.020848748529055869 * tmp_97;
      real_t tmp_141 = 0.78764240869137092 * tmp_22 + 0.1711304259088916 * tmp_23;
      real_t tmp_142 = tmp_19 * ( tmp_141 + tmp_20 );
      real_t tmp_143 = 0.78764240869137092 * tmp_29 + 0.1711304259088916 * tmp_30;
      real_t tmp_144 = tmp_19 * ( tmp_143 + tmp_27 );
      real_t tmp_145 = 0.78764240869137092 * tmp_36 + 0.1711304259088916 * tmp_37;
      real_t tmp_146 = tmp_19 * ( tmp_145 + tmp_34 );
      real_t tmp_147 = tmp_1 * ( tmp_142 * tmp_9 + tmp_144 * tmp_26 + tmp_146 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_142 * tmp_40 + tmp_144 * tmp_41 + tmp_146 * tmp_42 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_142 * tmp_43 + tmp_144 * tmp_44 + tmp_146 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_148 = tmp_141 + tmp_76;
      real_t tmp_149 = tmp_148 * tmp_72;
      real_t tmp_150 = tmp_148 * tmp_73;
      real_t tmp_151 = tmp_143 + tmp_80;
      real_t tmp_152 = tmp_151 * tmp_69;
      real_t tmp_153 = tmp_151 * tmp_70;
      real_t tmp_154 = tmp_148 * tmp_74;
      real_t tmp_155 = tmp_151 * tmp_71;
      real_t tmp_156 = tmp_145 + tmp_86;
      real_t tmp_157 = tmp_156 * tmp_66;
      real_t tmp_158 = tmp_156 * tmp_67;
      real_t tmp_159 = tmp_156 * tmp_68;
      real_t tmp_160 = tmp_147 * tmp_95;
      real_t tmp_161 = 0.019202922745021479 * tmp_97;
      real_t tmp_162 = 0.58463275527740355 * tmp_22 + 0.37605877282253791 * tmp_23;
      real_t tmp_163 = tmp_19 * ( tmp_162 + tmp_20 );
      real_t tmp_164 = 0.58463275527740355 * tmp_29 + 0.37605877282253791 * tmp_30;
      real_t tmp_165 = tmp_19 * ( tmp_164 + tmp_27 );
      real_t tmp_166 = 0.58463275527740355 * tmp_36 + 0.37605877282253791 * tmp_37;
      real_t tmp_167 = tmp_19 * ( tmp_166 + tmp_34 );
      real_t tmp_168 = tmp_1 * ( tmp_163 * tmp_9 + tmp_165 * tmp_26 + tmp_167 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_163 * tmp_40 + tmp_165 * tmp_41 + tmp_167 * tmp_42 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_163 * tmp_43 + tmp_165 * tmp_44 + tmp_167 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_169 = tmp_162 + tmp_76;
      real_t tmp_170 = tmp_169 * tmp_72;
      real_t tmp_171 = tmp_169 * tmp_73;
      real_t tmp_172 = tmp_164 + tmp_80;
      real_t tmp_173 = tmp_172 * tmp_69;
      real_t tmp_174 = tmp_172 * tmp_70;
      real_t tmp_175 = tmp_169 * tmp_74;
      real_t tmp_176 = tmp_172 * tmp_71;
      real_t tmp_177 = tmp_166 + tmp_86;
      real_t tmp_178 = tmp_177 * tmp_66;
      real_t tmp_179 = tmp_177 * tmp_67;
      real_t tmp_180 = tmp_177 * tmp_68;
      real_t tmp_181 = tmp_168 * tmp_95;
      real_t tmp_182 = 0.020848748529055869 * tmp_97;
      real_t tmp_183 = 0.041227165399737475 * tmp_22 + 0.78764240869137092 * tmp_23;
      real_t tmp_184 = tmp_19 * ( tmp_183 + tmp_20 );
      real_t tmp_185 = 0.041227165399737475 * tmp_29 + 0.78764240869137092 * tmp_30;
      real_t tmp_186 = tmp_19 * ( tmp_185 + tmp_27 );
      real_t tmp_187 = 0.041227165399737475 * tmp_36 + 0.78764240869137092 * tmp_37;
      real_t tmp_188 = tmp_19 * ( tmp_187 + tmp_34 );
      real_t tmp_189 = tmp_1 * ( tmp_184 * tmp_9 + tmp_186 * tmp_26 + tmp_188 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_184 * tmp_40 + tmp_186 * tmp_41 + tmp_188 * tmp_42 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_184 * tmp_43 + tmp_186 * tmp_44 + tmp_188 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_190 = tmp_183 + tmp_76;
      real_t tmp_191 = tmp_190 * tmp_72;
      real_t tmp_192 = tmp_190 * tmp_73;
      real_t tmp_193 = tmp_185 + tmp_80;
      real_t tmp_194 = tmp_193 * tmp_69;
      real_t tmp_195 = tmp_193 * tmp_70;
      real_t tmp_196 = tmp_190 * tmp_74;
      real_t tmp_197 = tmp_193 * tmp_71;
      real_t tmp_198 = tmp_187 + tmp_86;
      real_t tmp_199 = tmp_198 * tmp_66;
      real_t tmp_200 = tmp_198 * tmp_67;
      real_t tmp_201 = tmp_198 * tmp_68;
      real_t tmp_202 = tmp_189 * tmp_95;
      real_t tmp_203 = 0.019202922745021479 * tmp_97;
      real_t tmp_204 = 0.039308471900058539 * tmp_22 + 0.58463275527740355 * tmp_23;
      real_t tmp_205 = tmp_19 * ( tmp_20 + tmp_204 );
      real_t tmp_206 = 0.039308471900058539 * tmp_29 + 0.58463275527740355 * tmp_30;
      real_t tmp_207 = tmp_19 * ( tmp_206 + tmp_27 );
      real_t tmp_208 = 0.039308471900058539 * tmp_36 + 0.58463275527740355 * tmp_37;
      real_t tmp_209 = tmp_19 * ( tmp_208 + tmp_34 );
      real_t tmp_210 = tmp_1 * ( tmp_205 * tmp_9 + tmp_207 * tmp_26 + tmp_209 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_205 * tmp_40 + tmp_207 * tmp_41 + tmp_209 * tmp_42 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_205 * tmp_43 + tmp_207 * tmp_44 + tmp_209 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_211 = tmp_204 + tmp_76;
      real_t tmp_212 = tmp_211 * tmp_72;
      real_t tmp_213 = tmp_211 * tmp_73;
      real_t tmp_214 = tmp_206 + tmp_80;
      real_t tmp_215 = tmp_214 * tmp_69;
      real_t tmp_216 = tmp_214 * tmp_70;
      real_t tmp_217 = tmp_211 * tmp_74;
      real_t tmp_218 = tmp_214 * tmp_71;
      real_t tmp_219 = tmp_208 + tmp_86;
      real_t tmp_220 = tmp_219 * tmp_66;
      real_t tmp_221 = tmp_219 * tmp_67;
      real_t tmp_222 = tmp_219 * tmp_68;
      real_t tmp_223 = tmp_210 * tmp_95;
      real_t tmp_224 = 0.020848748529055869 * tmp_97;
      real_t tmp_225 = 0.78764240869137092 * tmp_22 + 0.041227165399737475 * tmp_23;
      real_t tmp_226 = tmp_19 * ( tmp_20 + tmp_225 );
      real_t tmp_227 = 0.78764240869137092 * tmp_29 + 0.041227165399737475 * tmp_30;
      real_t tmp_228 = tmp_19 * ( tmp_227 + tmp_27 );
      real_t tmp_229 = 0.78764240869137092 * tmp_36 + 0.041227165399737475 * tmp_37;
      real_t tmp_230 = tmp_19 * ( tmp_229 + tmp_34 );
      real_t tmp_231 = tmp_1 * ( tmp_226 * tmp_9 + tmp_228 * tmp_26 + tmp_230 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_226 * tmp_40 + tmp_228 * tmp_41 + tmp_230 * tmp_42 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_226 * tmp_43 + tmp_228 * tmp_44 + tmp_230 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_232 = tmp_225 + tmp_76;
      real_t tmp_233 = tmp_232 * tmp_72;
      real_t tmp_234 = tmp_232 * tmp_73;
      real_t tmp_235 = tmp_227 + tmp_80;
      real_t tmp_236 = tmp_235 * tmp_69;
      real_t tmp_237 = tmp_235 * tmp_70;
      real_t tmp_238 = tmp_232 * tmp_74;
      real_t tmp_239 = tmp_235 * tmp_71;
      real_t tmp_240 = tmp_229 + tmp_86;
      real_t tmp_241 = tmp_240 * tmp_66;
      real_t tmp_242 = tmp_240 * tmp_67;
      real_t tmp_243 = tmp_240 * tmp_68;
      real_t tmp_244 = tmp_231 * tmp_95;
      real_t tmp_245 = 0.019202922745021479 * tmp_97;
      real_t tmp_246 = 0.58463275527740355 * tmp_22 + 0.039308471900058539 * tmp_23;
      real_t tmp_247 = tmp_19 * ( tmp_20 + tmp_246 );
      real_t tmp_248 = 0.58463275527740355 * tmp_29 + 0.039308471900058539 * tmp_30;
      real_t tmp_249 = tmp_19 * ( tmp_248 + tmp_27 );
      real_t tmp_250 = 0.58463275527740355 * tmp_36 + 0.039308471900058539 * tmp_37;
      real_t tmp_251 = tmp_19 * ( tmp_250 + tmp_34 );
      real_t tmp_252 = tmp_1 * ( tmp_247 * tmp_9 + tmp_249 * tmp_26 + tmp_251 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_247 * tmp_40 + tmp_249 * tmp_41 + tmp_251 * tmp_42 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_247 * tmp_43 + tmp_249 * tmp_44 + tmp_251 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_253 = tmp_246 + tmp_76;
      real_t tmp_254 = tmp_253 * tmp_72;
      real_t tmp_255 = tmp_253 * tmp_73;
      real_t tmp_256 = tmp_248 + tmp_80;
      real_t tmp_257 = tmp_256 * tmp_69;
      real_t tmp_258 = tmp_256 * tmp_70;
      real_t tmp_259 = tmp_253 * tmp_74;
      real_t tmp_260 = tmp_256 * tmp_71;
      real_t tmp_261 = tmp_250 + tmp_86;
      real_t tmp_262 = tmp_261 * tmp_66;
      real_t tmp_263 = tmp_261 * tmp_67;
      real_t tmp_264 = tmp_261 * tmp_68;
      real_t tmp_265 = tmp_252 * tmp_95;
      real_t tmp_266 = 0.020848748529055869 * tmp_97;
      real_t tmp_267 = 0.1711304259088916 * tmp_22 + 0.78764240869137092 * tmp_23;
      real_t tmp_268 = tmp_19 * ( tmp_20 + tmp_267 );
      real_t tmp_269 = 0.1711304259088916 * tmp_29 + 0.78764240869137092 * tmp_30;
      real_t tmp_270 = tmp_19 * ( tmp_269 + tmp_27 );
      real_t tmp_271 = 0.1711304259088916 * tmp_36 + 0.78764240869137092 * tmp_37;
      real_t tmp_272 = tmp_19 * ( tmp_271 + tmp_34 );
      real_t tmp_273 = tmp_1 * ( tmp_26 * tmp_270 + tmp_268 * tmp_9 + tmp_272 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_268 * tmp_40 + tmp_270 * tmp_41 + tmp_272 * tmp_42 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_268 * tmp_43 + tmp_270 * tmp_44 + tmp_272 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_274 = tmp_267 + tmp_76;
      real_t tmp_275 = tmp_274 * tmp_72;
      real_t tmp_276 = tmp_274 * tmp_73;
      real_t tmp_277 = tmp_269 + tmp_80;
      real_t tmp_278 = tmp_277 * tmp_69;
      real_t tmp_279 = tmp_277 * tmp_70;
      real_t tmp_280 = tmp_274 * tmp_74;
      real_t tmp_281 = tmp_277 * tmp_71;
      real_t tmp_282 = tmp_271 + tmp_86;
      real_t tmp_283 = tmp_282 * tmp_66;
      real_t tmp_284 = tmp_282 * tmp_67;
      real_t tmp_285 = tmp_282 * tmp_68;
      real_t tmp_286 = tmp_273 * tmp_95;
      real_t tmp_287 = 0.019202922745021479 * tmp_97;
      real_t tmp_288 = 0.37605877282253791 * tmp_22 + 0.58463275527740355 * tmp_23;
      real_t tmp_289 = tmp_19 * ( tmp_20 + tmp_288 );
      real_t tmp_290 = 0.37605877282253791 * tmp_29 + 0.58463275527740355 * tmp_30;
      real_t tmp_291 = tmp_19 * ( tmp_27 + tmp_290 );
      real_t tmp_292 = 0.37605877282253791 * tmp_36 + 0.58463275527740355 * tmp_37;
      real_t tmp_293 = tmp_19 * ( tmp_292 + tmp_34 );
      real_t tmp_294 = tmp_1 * ( tmp_26 * tmp_291 + tmp_289 * tmp_9 + tmp_293 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_289 * tmp_40 + tmp_291 * tmp_41 + tmp_293 * tmp_42 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_289 * tmp_43 + tmp_291 * tmp_44 + tmp_293 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_295 = tmp_288 + tmp_76;
      real_t tmp_296 = tmp_295 * tmp_72;
      real_t tmp_297 = tmp_295 * tmp_73;
      real_t tmp_298 = tmp_290 + tmp_80;
      real_t tmp_299 = tmp_298 * tmp_69;
      real_t tmp_300 = tmp_298 * tmp_70;
      real_t tmp_301 = tmp_295 * tmp_74;
      real_t tmp_302 = tmp_298 * tmp_71;
      real_t tmp_303 = tmp_292 + tmp_86;
      real_t tmp_304 = tmp_303 * tmp_66;
      real_t tmp_305 = tmp_303 * tmp_67;
      real_t tmp_306 = tmp_303 * tmp_68;
      real_t tmp_307 = tmp_294 * tmp_95;
      real_t tmp_308 = 0.020848748529055869 * tmp_97;
      real_t tmp_309 = 0.041227165399737475 * tmp_22 + 0.1711304259088916 * tmp_23;
      real_t tmp_310 = tmp_19 * ( tmp_20 + tmp_309 );
      real_t tmp_311 = 0.041227165399737475 * tmp_29 + 0.1711304259088916 * tmp_30;
      real_t tmp_312 = tmp_19 * ( tmp_27 + tmp_311 );
      real_t tmp_313 = 0.041227165399737475 * tmp_36 + 0.1711304259088916 * tmp_37;
      real_t tmp_314 = tmp_19 * ( tmp_313 + tmp_34 );
      real_t tmp_315 = tmp_1 * ( tmp_26 * tmp_312 + tmp_310 * tmp_9 + tmp_314 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_310 * tmp_40 + tmp_312 * tmp_41 + tmp_314 * tmp_42 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_310 * tmp_43 + tmp_312 * tmp_44 + tmp_314 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_316 = tmp_309 + tmp_76;
      real_t tmp_317 = tmp_316 * tmp_72;
      real_t tmp_318 = tmp_316 * tmp_73;
      real_t tmp_319 = tmp_311 + tmp_80;
      real_t tmp_320 = tmp_319 * tmp_69;
      real_t tmp_321 = tmp_319 * tmp_70;
      real_t tmp_322 = tmp_316 * tmp_74;
      real_t tmp_323 = tmp_319 * tmp_71;
      real_t tmp_324 = tmp_313 + tmp_86;
      real_t tmp_325 = tmp_324 * tmp_66;
      real_t tmp_326 = tmp_324 * tmp_67;
      real_t tmp_327 = tmp_324 * tmp_68;
      real_t tmp_328 = tmp_315 * tmp_95;
      real_t tmp_329 = 0.019202922745021479 * tmp_97;
      real_t tmp_330 = 0.40446199974765351 * tmp_22 + 0.19107600050469298 * tmp_23;
      real_t tmp_331 = tmp_19 * ( tmp_20 + tmp_330 );
      real_t tmp_332 = 0.40446199974765351 * tmp_29 + 0.19107600050469298 * tmp_30;
      real_t tmp_333 = tmp_19 * ( tmp_27 + tmp_332 );
      real_t tmp_334 = 0.40446199974765351 * tmp_36 + 0.19107600050469298 * tmp_37;
      real_t tmp_335 = tmp_19 * ( tmp_334 + tmp_34 );
      real_t tmp_336 = tmp_1 * ( tmp_26 * tmp_333 + tmp_33 * tmp_335 + tmp_331 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_331 * tmp_40 + tmp_333 * tmp_41 + tmp_335 * tmp_42 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_331 * tmp_43 + tmp_333 * tmp_44 + tmp_335 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_337 = tmp_330 + tmp_76;
      real_t tmp_338 = tmp_337 * tmp_72;
      real_t tmp_339 = tmp_337 * tmp_73;
      real_t tmp_340 = tmp_332 + tmp_80;
      real_t tmp_341 = tmp_340 * tmp_69;
      real_t tmp_342 = tmp_340 * tmp_70;
      real_t tmp_343 = tmp_337 * tmp_74;
      real_t tmp_344 = tmp_340 * tmp_71;
      real_t tmp_345 = tmp_334 + tmp_86;
      real_t tmp_346 = tmp_345 * tmp_66;
      real_t tmp_347 = tmp_345 * tmp_67;
      real_t tmp_348 = tmp_345 * tmp_68;
      real_t tmp_349 = tmp_336 * tmp_95;
      real_t tmp_350 = 0.042507265838595799 * tmp_97;
      real_t tmp_351 = 0.039308471900058539 * tmp_22 + 0.37605877282253791 * tmp_23;
      real_t tmp_352 = tmp_19 * ( tmp_20 + tmp_351 );
      real_t tmp_353 = 0.039308471900058539 * tmp_29 + 0.37605877282253791 * tmp_30;
      real_t tmp_354 = tmp_19 * ( tmp_27 + tmp_353 );
      real_t tmp_355 = 0.039308471900058539 * tmp_36 + 0.37605877282253791 * tmp_37;
      real_t tmp_356 = tmp_19 * ( tmp_34 + tmp_355 );
      real_t tmp_357 = tmp_1 * ( tmp_26 * tmp_354 + tmp_33 * tmp_356 + tmp_352 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_352 * tmp_40 + tmp_354 * tmp_41 + tmp_356 * tmp_42 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_352 * tmp_43 + tmp_354 * tmp_44 + tmp_356 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_358 = tmp_351 + tmp_76;
      real_t tmp_359 = tmp_358 * tmp_72;
      real_t tmp_360 = tmp_358 * tmp_73;
      real_t tmp_361 = tmp_353 + tmp_80;
      real_t tmp_362 = tmp_361 * tmp_69;
      real_t tmp_363 = tmp_361 * tmp_70;
      real_t tmp_364 = tmp_358 * tmp_74;
      real_t tmp_365 = tmp_361 * tmp_71;
      real_t tmp_366 = tmp_355 + tmp_86;
      real_t tmp_367 = tmp_366 * tmp_66;
      real_t tmp_368 = tmp_366 * tmp_67;
      real_t tmp_369 = tmp_366 * tmp_68;
      real_t tmp_370 = tmp_357 * tmp_95;
      real_t tmp_371 = 0.020848748529055869 * tmp_97;
      real_t tmp_372 = 0.93718850182767688 * tmp_22 + 0.031405749086161582 * tmp_23;
      real_t tmp_373 = tmp_19 * ( tmp_20 + tmp_372 );
      real_t tmp_374 = 0.93718850182767688 * tmp_29 + 0.031405749086161582 * tmp_30;
      real_t tmp_375 = tmp_19 * ( tmp_27 + tmp_374 );
      real_t tmp_376 = 0.93718850182767688 * tmp_36 + 0.031405749086161582 * tmp_37;
      real_t tmp_377 = tmp_19 * ( tmp_34 + tmp_376 );
      real_t tmp_378 = tmp_1 * ( tmp_26 * tmp_375 + tmp_33 * tmp_377 + tmp_373 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_373 * tmp_40 + tmp_375 * tmp_41 + tmp_377 * tmp_42 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_373 * tmp_43 + tmp_375 * tmp_44 + tmp_377 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_379 = tmp_372 + tmp_76;
      real_t tmp_380 = tmp_379 * tmp_72;
      real_t tmp_381 = tmp_379 * tmp_73;
      real_t tmp_382 = tmp_374 + tmp_80;
      real_t tmp_383 = tmp_382 * tmp_69;
      real_t tmp_384 = tmp_382 * tmp_70;
      real_t tmp_385 = tmp_379 * tmp_74;
      real_t tmp_386 = tmp_382 * tmp_71;
      real_t tmp_387 = tmp_376 + tmp_86;
      real_t tmp_388 = tmp_387 * tmp_66;
      real_t tmp_389 = tmp_387 * tmp_67;
      real_t tmp_390 = tmp_387 * tmp_68;
      real_t tmp_391 = tmp_378 * tmp_95;
      real_t tmp_392 = 0.0068572537431980923 * tmp_97;
      real_t tmp_393 = 0.60796128279561268 * tmp_22 + 0.19601935860219369 * tmp_23;
      real_t tmp_394 = tmp_19 * ( tmp_20 + tmp_393 );
      real_t tmp_395 = 0.60796128279561268 * tmp_29 + 0.19601935860219369 * tmp_30;
      real_t tmp_396 = tmp_19 * ( tmp_27 + tmp_395 );
      real_t tmp_397 = 0.60796128279561268 * tmp_36 + 0.19601935860219369 * tmp_37;
      real_t tmp_398 = tmp_19 * ( tmp_34 + tmp_397 );
      real_t tmp_399 = tmp_1 * ( tmp_26 * tmp_396 + tmp_33 * tmp_398 + tmp_394 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_394 * tmp_40 + tmp_396 * tmp_41 + tmp_398 * tmp_42 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_394 * tmp_43 + tmp_396 * tmp_44 + tmp_398 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_400 = tmp_393 + tmp_76;
      real_t tmp_401 = tmp_400 * tmp_72;
      real_t tmp_402 = tmp_400 * tmp_73;
      real_t tmp_403 = tmp_395 + tmp_80;
      real_t tmp_404 = tmp_403 * tmp_69;
      real_t tmp_405 = tmp_403 * tmp_70;
      real_t tmp_406 = tmp_400 * tmp_74;
      real_t tmp_407 = tmp_403 * tmp_71;
      real_t tmp_408 = tmp_397 + tmp_86;
      real_t tmp_409 = tmp_408 * tmp_66;
      real_t tmp_410 = tmp_408 * tmp_67;
      real_t tmp_411 = tmp_408 * tmp_68;
      real_t tmp_412 = tmp_399 * tmp_95;
      real_t tmp_413 = 0.037198804536718075 * tmp_97;
      real_t tmp_414 = 0.19107600050469298 * tmp_22 + 0.40446199974765351 * tmp_23;
      real_t tmp_415 = tmp_19 * ( tmp_20 + tmp_414 );
      real_t tmp_416 = 0.19107600050469298 * tmp_29 + 0.40446199974765351 * tmp_30;
      real_t tmp_417 = tmp_19 * ( tmp_27 + tmp_416 );
      real_t tmp_418 = 0.19107600050469298 * tmp_36 + 0.40446199974765351 * tmp_37;
      real_t tmp_419 = tmp_19 * ( tmp_34 + tmp_418 );
      real_t tmp_420 = tmp_1 * ( tmp_26 * tmp_417 + tmp_33 * tmp_419 + tmp_415 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_40 * tmp_415 + tmp_41 * tmp_417 + tmp_419 * tmp_42 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_415 * tmp_43 + tmp_417 * tmp_44 + tmp_419 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_421 = tmp_414 + tmp_76;
      real_t tmp_422 = tmp_421 * tmp_72;
      real_t tmp_423 = tmp_421 * tmp_73;
      real_t tmp_424 = tmp_416 + tmp_80;
      real_t tmp_425 = tmp_424 * tmp_69;
      real_t tmp_426 = tmp_424 * tmp_70;
      real_t tmp_427 = tmp_421 * tmp_74;
      real_t tmp_428 = tmp_424 * tmp_71;
      real_t tmp_429 = tmp_418 + tmp_86;
      real_t tmp_430 = tmp_429 * tmp_66;
      real_t tmp_431 = tmp_429 * tmp_67;
      real_t tmp_432 = tmp_429 * tmp_68;
      real_t tmp_433 = tmp_420 * tmp_95;
      real_t tmp_434 = 0.042507265838595799 * tmp_97;
      real_t tmp_435 = 0.031405749086161582 * tmp_22 + 0.031405749086161582 * tmp_23;
      real_t tmp_436 = tmp_19 * ( tmp_20 + tmp_435 );
      real_t tmp_437 = 0.031405749086161582 * tmp_29 + 0.031405749086161582 * tmp_30;
      real_t tmp_438 = tmp_19 * ( tmp_27 + tmp_437 );
      real_t tmp_439 = 0.031405749086161582 * tmp_36 + 0.031405749086161582 * tmp_37;
      real_t tmp_440 = tmp_19 * ( tmp_34 + tmp_439 );
      real_t tmp_441 = tmp_1 * ( tmp_26 * tmp_438 + tmp_33 * tmp_440 + tmp_436 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_40 * tmp_436 + tmp_41 * tmp_438 + tmp_42 * tmp_440 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_43 * tmp_436 + tmp_438 * tmp_44 + tmp_440 * tmp_45 - 1.0 / 4.0 );
      real_t tmp_442 = tmp_435 + tmp_76;
      real_t tmp_443 = tmp_442 * tmp_72;
      real_t tmp_444 = tmp_442 * tmp_73;
      real_t tmp_445 = tmp_437 + tmp_80;
      real_t tmp_446 = tmp_445 * tmp_69;
      real_t tmp_447 = tmp_445 * tmp_70;
      real_t tmp_448 = tmp_442 * tmp_74;
      real_t tmp_449 = tmp_445 * tmp_71;
      real_t tmp_450 = tmp_439 + tmp_86;
      real_t tmp_451 = tmp_450 * tmp_66;
      real_t tmp_452 = tmp_450 * tmp_67;
      real_t tmp_453 = tmp_450 * tmp_68;
      real_t tmp_454 = tmp_441 * tmp_95;
      real_t tmp_455 = 0.0068572537431980923 * tmp_97;
      real_t tmp_456 = 0.19601935860219369 * tmp_22 + 0.19601935860219369 * tmp_23;
      real_t tmp_457 = tmp_19 * ( tmp_20 + tmp_456 );
      real_t tmp_458 = 0.19601935860219369 * tmp_29 + 0.19601935860219369 * tmp_30;
      real_t tmp_459 = tmp_19 * ( tmp_27 + tmp_458 );
      real_t tmp_460 = 0.19601935860219369 * tmp_36 + 0.19601935860219369 * tmp_37;
      real_t tmp_461 = tmp_19 * ( tmp_34 + tmp_460 );
      real_t tmp_462 = tmp_1 * ( tmp_26 * tmp_459 + tmp_33 * tmp_461 + tmp_457 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_40 * tmp_457 + tmp_41 * tmp_459 + tmp_42 * tmp_461 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_43 * tmp_457 + tmp_44 * tmp_459 + tmp_45 * tmp_461 - 1.0 / 4.0 );
      real_t tmp_463 = tmp_456 + tmp_76;
      real_t tmp_464 = tmp_463 * tmp_72;
      real_t tmp_465 = tmp_463 * tmp_73;
      real_t tmp_466 = tmp_458 + tmp_80;
      real_t tmp_467 = tmp_466 * tmp_69;
      real_t tmp_468 = tmp_466 * tmp_70;
      real_t tmp_469 = tmp_463 * tmp_74;
      real_t tmp_470 = tmp_466 * tmp_71;
      real_t tmp_471 = tmp_460 + tmp_86;
      real_t tmp_472 = tmp_471 * tmp_66;
      real_t tmp_473 = tmp_471 * tmp_67;
      real_t tmp_474 = tmp_471 * tmp_68;
      real_t tmp_475 = tmp_462 * tmp_95;
      real_t tmp_476 = 0.037198804536718075 * tmp_97;
      real_t tmp_477 = 0.40446199974765351 * tmp_22 + 0.40446199974765351 * tmp_23;
      real_t tmp_478 = tmp_19 * ( tmp_20 + tmp_477 );
      real_t tmp_479 = 0.40446199974765351 * tmp_29 + 0.40446199974765351 * tmp_30;
      real_t tmp_480 = tmp_19 * ( tmp_27 + tmp_479 );
      real_t tmp_481 = 0.40446199974765351 * tmp_36 + 0.40446199974765351 * tmp_37;
      real_t tmp_482 = tmp_19 * ( tmp_34 + tmp_481 );
      real_t tmp_483 = tmp_1 * ( tmp_26 * tmp_480 + tmp_33 * tmp_482 + tmp_478 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_40 * tmp_478 + tmp_41 * tmp_480 + tmp_42 * tmp_482 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_43 * tmp_478 + tmp_44 * tmp_480 + tmp_45 * tmp_482 - 1.0 / 4.0 );
      real_t tmp_484 = tmp_477 + tmp_76;
      real_t tmp_485 = tmp_484 * tmp_72;
      real_t tmp_486 = tmp_484 * tmp_73;
      real_t tmp_487 = tmp_479 + tmp_80;
      real_t tmp_488 = tmp_487 * tmp_69;
      real_t tmp_489 = tmp_487 * tmp_70;
      real_t tmp_490 = tmp_484 * tmp_74;
      real_t tmp_491 = tmp_487 * tmp_71;
      real_t tmp_492 = tmp_481 + tmp_86;
      real_t tmp_493 = tmp_492 * tmp_66;
      real_t tmp_494 = tmp_492 * tmp_67;
      real_t tmp_495 = tmp_492 * tmp_68;
      real_t tmp_496 = tmp_483 * tmp_95;
      real_t tmp_497 = 0.042507265838595799 * tmp_97;
      real_t tmp_498 = 0.1711304259088916 * tmp_22 + 0.041227165399737475 * tmp_23;
      real_t tmp_499 = tmp_19 * ( tmp_20 + tmp_498 );
      real_t tmp_500 = 0.1711304259088916 * tmp_29 + 0.041227165399737475 * tmp_30;
      real_t tmp_501 = tmp_19 * ( tmp_27 + tmp_500 );
      real_t tmp_502 = 0.1711304259088916 * tmp_36 + 0.041227165399737475 * tmp_37;
      real_t tmp_503 = tmp_19 * ( tmp_34 + tmp_502 );
      real_t tmp_504 = tmp_1 * ( tmp_26 * tmp_501 + tmp_33 * tmp_503 + tmp_499 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_2 * ( tmp_40 * tmp_499 + tmp_41 * tmp_501 + tmp_42 * tmp_503 - 1.0 / 4.0 ) +
                       tmp_6 * ( tmp_43 * tmp_499 + tmp_44 * tmp_501 + tmp_45 * tmp_503 - 1.0 / 4.0 );
      real_t tmp_505 = tmp_498 + tmp_76;
      real_t tmp_506 = tmp_505 * tmp_72;
      real_t tmp_507 = tmp_505 * tmp_73;
      real_t tmp_508 = tmp_500 + tmp_80;
      real_t tmp_509 = tmp_508 * tmp_69;
      real_t tmp_510 = tmp_508 * tmp_70;
      real_t tmp_511 = tmp_505 * tmp_74;
      real_t tmp_512 = tmp_508 * tmp_71;
      real_t tmp_513 = tmp_502 + tmp_86;
      real_t tmp_514 = tmp_513 * tmp_66;
      real_t tmp_515 = tmp_513 * tmp_67;
      real_t tmp_516 = tmp_513 * tmp_68;
      real_t tmp_517 = tmp_504 * tmp_95;
      real_t tmp_518 = 0.019202922745021479 * tmp_97;
      real_t tmp_519 = 0.5 * p_affine_13_0 * tmp_68 + 0.5 * p_affine_13_1 * tmp_71 + 0.5 * p_affine_13_2 * tmp_74;
      real_t tmp_520 = 0.5 * p_affine_13_0 * tmp_67 + 0.5 * p_affine_13_1 * tmp_70 + 0.5 * p_affine_13_2 * tmp_73;
      real_t tmp_521 = 0.5 * p_affine_13_0 * tmp_66 + 0.5 * p_affine_13_1 * tmp_69 + 0.5 * p_affine_13_2 * tmp_72;
      real_t a_0_0   = tmp_119 * ( -tmp_105 * tmp_75 - tmp_118 * ( -tmp_107 - tmp_108 - tmp_110 - tmp_111 - tmp_112 - tmp_113 -
                                                                 tmp_115 - tmp_116 - tmp_117 + 1 ) ) +
                     tmp_140 * ( -tmp_126 * tmp_75 - tmp_139 * ( -tmp_128 - tmp_129 - tmp_131 - tmp_132 - tmp_133 - tmp_134 -
                                                                 tmp_136 - tmp_137 - tmp_138 + 1 ) ) +
                     tmp_161 * ( -tmp_147 * tmp_75 - tmp_160 * ( -tmp_149 - tmp_150 - tmp_152 - tmp_153 - tmp_154 - tmp_155 -
                                                                 tmp_157 - tmp_158 - tmp_159 + 1 ) ) +
                     tmp_182 * ( -tmp_168 * tmp_75 - tmp_181 * ( -tmp_170 - tmp_171 - tmp_173 - tmp_174 - tmp_175 - tmp_176 -
                                                                 tmp_178 - tmp_179 - tmp_180 + 1 ) ) +
                     tmp_203 * ( -tmp_189 * tmp_75 - tmp_202 * ( -tmp_191 - tmp_192 - tmp_194 - tmp_195 - tmp_196 - tmp_197 -
                                                                 tmp_199 - tmp_200 - tmp_201 + 1 ) ) +
                     tmp_224 * ( -tmp_210 * tmp_75 - tmp_223 * ( -tmp_212 - tmp_213 - tmp_215 - tmp_216 - tmp_217 - tmp_218 -
                                                                 tmp_220 - tmp_221 - tmp_222 + 1 ) ) +
                     tmp_245 * ( -tmp_231 * tmp_75 - tmp_244 * ( -tmp_233 - tmp_234 - tmp_236 - tmp_237 - tmp_238 - tmp_239 -
                                                                 tmp_241 - tmp_242 - tmp_243 + 1 ) ) +
                     tmp_266 * ( -tmp_252 * tmp_75 - tmp_265 * ( -tmp_254 - tmp_255 - tmp_257 - tmp_258 - tmp_259 - tmp_260 -
                                                                 tmp_262 - tmp_263 - tmp_264 + 1 ) ) +
                     tmp_287 * ( -tmp_273 * tmp_75 - tmp_286 * ( -tmp_275 - tmp_276 - tmp_278 - tmp_279 - tmp_280 - tmp_281 -
                                                                 tmp_283 - tmp_284 - tmp_285 + 1 ) ) +
                     tmp_308 * ( -tmp_294 * tmp_75 - tmp_307 * ( -tmp_296 - tmp_297 - tmp_299 - tmp_300 - tmp_301 - tmp_302 -
                                                                 tmp_304 - tmp_305 - tmp_306 + 1 ) ) +
                     tmp_329 * ( -tmp_315 * tmp_75 - tmp_328 * ( -tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 -
                                                                 tmp_325 - tmp_326 - tmp_327 + 1 ) ) +
                     tmp_350 * ( -tmp_336 * tmp_75 - tmp_349 * ( -tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 - tmp_344 -
                                                                 tmp_346 - tmp_347 - tmp_348 + 1 ) ) +
                     tmp_371 * ( -tmp_357 * tmp_75 - tmp_370 * ( -tmp_359 - tmp_360 - tmp_362 - tmp_363 - tmp_364 - tmp_365 -
                                                                 tmp_367 - tmp_368 - tmp_369 + 1 ) ) +
                     tmp_392 * ( -tmp_378 * tmp_75 - tmp_391 * ( -tmp_380 - tmp_381 - tmp_383 - tmp_384 - tmp_385 - tmp_386 -
                                                                 tmp_388 - tmp_389 - tmp_390 + 1 ) ) +
                     tmp_413 * ( -tmp_399 * tmp_75 - tmp_412 * ( -tmp_401 - tmp_402 - tmp_404 - tmp_405 - tmp_406 - tmp_407 -
                                                                 tmp_409 - tmp_410 - tmp_411 + 1 ) ) +
                     tmp_434 * ( -tmp_420 * tmp_75 - tmp_433 * ( -tmp_422 - tmp_423 - tmp_425 - tmp_426 - tmp_427 - tmp_428 -
                                                                 tmp_430 - tmp_431 - tmp_432 + 1 ) ) +
                     tmp_455 * ( -tmp_441 * tmp_75 - tmp_454 * ( -tmp_443 - tmp_444 - tmp_446 - tmp_447 - tmp_448 - tmp_449 -
                                                                 tmp_451 - tmp_452 - tmp_453 + 1 ) ) +
                     tmp_476 * ( -tmp_462 * tmp_75 - tmp_475 * ( -tmp_464 - tmp_465 - tmp_467 - tmp_468 - tmp_469 - tmp_470 -
                                                                 tmp_472 - tmp_473 - tmp_474 + 1 ) ) +
                     tmp_497 * ( -tmp_483 * tmp_75 - tmp_496 * ( -tmp_485 - tmp_486 - tmp_488 - tmp_489 - tmp_490 - tmp_491 -
                                                                 tmp_493 - tmp_494 - tmp_495 + 1 ) ) +
                     tmp_518 * ( -tmp_504 * tmp_75 - tmp_517 * ( -tmp_506 - tmp_507 - tmp_509 - tmp_510 - tmp_511 - tmp_512 -
                                                                 tmp_514 - tmp_515 - tmp_516 + 1 ) ) +
                     tmp_98 * ( -tmp_46 * tmp_75 - tmp_96 * ( -tmp_78 - tmp_79 - tmp_82 - tmp_83 - tmp_84 - tmp_85 - tmp_88 -
                                                              tmp_89 - tmp_90 + 1 ) );
      real_t a_0_1 = tmp_119 * ( -tmp_105 * tmp_519 - tmp_118 * ( tmp_112 + tmp_113 + tmp_117 ) ) +
                     tmp_140 * ( -tmp_126 * tmp_519 - tmp_139 * ( tmp_133 + tmp_134 + tmp_138 ) ) +
                     tmp_161 * ( -tmp_147 * tmp_519 - tmp_160 * ( tmp_154 + tmp_155 + tmp_159 ) ) +
                     tmp_182 * ( -tmp_168 * tmp_519 - tmp_181 * ( tmp_175 + tmp_176 + tmp_180 ) ) +
                     tmp_203 * ( -tmp_189 * tmp_519 - tmp_202 * ( tmp_196 + tmp_197 + tmp_201 ) ) +
                     tmp_224 * ( -tmp_210 * tmp_519 - tmp_223 * ( tmp_217 + tmp_218 + tmp_222 ) ) +
                     tmp_245 * ( -tmp_231 * tmp_519 - tmp_244 * ( tmp_238 + tmp_239 + tmp_243 ) ) +
                     tmp_266 * ( -tmp_252 * tmp_519 - tmp_265 * ( tmp_259 + tmp_260 + tmp_264 ) ) +
                     tmp_287 * ( -tmp_273 * tmp_519 - tmp_286 * ( tmp_280 + tmp_281 + tmp_285 ) ) +
                     tmp_308 * ( -tmp_294 * tmp_519 - tmp_307 * ( tmp_301 + tmp_302 + tmp_306 ) ) +
                     tmp_329 * ( -tmp_315 * tmp_519 - tmp_328 * ( tmp_322 + tmp_323 + tmp_327 ) ) +
                     tmp_350 * ( -tmp_336 * tmp_519 - tmp_349 * ( tmp_343 + tmp_344 + tmp_348 ) ) +
                     tmp_371 * ( -tmp_357 * tmp_519 - tmp_370 * ( tmp_364 + tmp_365 + tmp_369 ) ) +
                     tmp_392 * ( -tmp_378 * tmp_519 - tmp_391 * ( tmp_385 + tmp_386 + tmp_390 ) ) +
                     tmp_413 * ( -tmp_399 * tmp_519 - tmp_412 * ( tmp_406 + tmp_407 + tmp_411 ) ) +
                     tmp_434 * ( -tmp_420 * tmp_519 - tmp_433 * ( tmp_427 + tmp_428 + tmp_432 ) ) +
                     tmp_455 * ( -tmp_441 * tmp_519 - tmp_454 * ( tmp_448 + tmp_449 + tmp_453 ) ) +
                     tmp_476 * ( -tmp_462 * tmp_519 - tmp_475 * ( tmp_469 + tmp_470 + tmp_474 ) ) +
                     tmp_497 * ( -tmp_483 * tmp_519 - tmp_496 * ( tmp_490 + tmp_491 + tmp_495 ) ) +
                     tmp_518 * ( -tmp_504 * tmp_519 - tmp_517 * ( tmp_511 + tmp_512 + tmp_516 ) ) +
                     tmp_98 * ( -tmp_46 * tmp_519 - tmp_96 * ( tmp_84 + tmp_85 + tmp_90 ) );
      real_t a_0_2 = tmp_119 * ( -tmp_105 * tmp_520 - tmp_118 * ( tmp_108 + tmp_111 + tmp_116 ) ) +
                     tmp_140 * ( -tmp_126 * tmp_520 - tmp_139 * ( tmp_129 + tmp_132 + tmp_137 ) ) +
                     tmp_161 * ( -tmp_147 * tmp_520 - tmp_160 * ( tmp_150 + tmp_153 + tmp_158 ) ) +
                     tmp_182 * ( -tmp_168 * tmp_520 - tmp_181 * ( tmp_171 + tmp_174 + tmp_179 ) ) +
                     tmp_203 * ( -tmp_189 * tmp_520 - tmp_202 * ( tmp_192 + tmp_195 + tmp_200 ) ) +
                     tmp_224 * ( -tmp_210 * tmp_520 - tmp_223 * ( tmp_213 + tmp_216 + tmp_221 ) ) +
                     tmp_245 * ( -tmp_231 * tmp_520 - tmp_244 * ( tmp_234 + tmp_237 + tmp_242 ) ) +
                     tmp_266 * ( -tmp_252 * tmp_520 - tmp_265 * ( tmp_255 + tmp_258 + tmp_263 ) ) +
                     tmp_287 * ( -tmp_273 * tmp_520 - tmp_286 * ( tmp_276 + tmp_279 + tmp_284 ) ) +
                     tmp_308 * ( -tmp_294 * tmp_520 - tmp_307 * ( tmp_297 + tmp_300 + tmp_305 ) ) +
                     tmp_329 * ( -tmp_315 * tmp_520 - tmp_328 * ( tmp_318 + tmp_321 + tmp_326 ) ) +
                     tmp_350 * ( -tmp_336 * tmp_520 - tmp_349 * ( tmp_339 + tmp_342 + tmp_347 ) ) +
                     tmp_371 * ( -tmp_357 * tmp_520 - tmp_370 * ( tmp_360 + tmp_363 + tmp_368 ) ) +
                     tmp_392 * ( -tmp_378 * tmp_520 - tmp_391 * ( tmp_381 + tmp_384 + tmp_389 ) ) +
                     tmp_413 * ( -tmp_399 * tmp_520 - tmp_412 * ( tmp_402 + tmp_405 + tmp_410 ) ) +
                     tmp_434 * ( -tmp_420 * tmp_520 - tmp_433 * ( tmp_423 + tmp_426 + tmp_431 ) ) +
                     tmp_455 * ( -tmp_441 * tmp_520 - tmp_454 * ( tmp_444 + tmp_447 + tmp_452 ) ) +
                     tmp_476 * ( -tmp_462 * tmp_520 - tmp_475 * ( tmp_465 + tmp_468 + tmp_473 ) ) +
                     tmp_497 * ( -tmp_483 * tmp_520 - tmp_496 * ( tmp_486 + tmp_489 + tmp_494 ) ) +
                     tmp_518 * ( -tmp_504 * tmp_520 - tmp_517 * ( tmp_507 + tmp_510 + tmp_515 ) ) +
                     tmp_98 * ( -tmp_46 * tmp_520 - tmp_96 * ( tmp_79 + tmp_83 + tmp_89 ) );
      real_t a_0_3 = tmp_119 * ( -tmp_105 * tmp_521 - tmp_118 * ( tmp_107 + tmp_110 + tmp_115 ) ) +
                     tmp_140 * ( -tmp_126 * tmp_521 - tmp_139 * ( tmp_128 + tmp_131 + tmp_136 ) ) +
                     tmp_161 * ( -tmp_147 * tmp_521 - tmp_160 * ( tmp_149 + tmp_152 + tmp_157 ) ) +
                     tmp_182 * ( -tmp_168 * tmp_521 - tmp_181 * ( tmp_170 + tmp_173 + tmp_178 ) ) +
                     tmp_203 * ( -tmp_189 * tmp_521 - tmp_202 * ( tmp_191 + tmp_194 + tmp_199 ) ) +
                     tmp_224 * ( -tmp_210 * tmp_521 - tmp_223 * ( tmp_212 + tmp_215 + tmp_220 ) ) +
                     tmp_245 * ( -tmp_231 * tmp_521 - tmp_244 * ( tmp_233 + tmp_236 + tmp_241 ) ) +
                     tmp_266 * ( -tmp_252 * tmp_521 - tmp_265 * ( tmp_254 + tmp_257 + tmp_262 ) ) +
                     tmp_287 * ( -tmp_273 * tmp_521 - tmp_286 * ( tmp_275 + tmp_278 + tmp_283 ) ) +
                     tmp_308 * ( -tmp_294 * tmp_521 - tmp_307 * ( tmp_296 + tmp_299 + tmp_304 ) ) +
                     tmp_329 * ( -tmp_315 * tmp_521 - tmp_328 * ( tmp_317 + tmp_320 + tmp_325 ) ) +
                     tmp_350 * ( -tmp_336 * tmp_521 - tmp_349 * ( tmp_338 + tmp_341 + tmp_346 ) ) +
                     tmp_371 * ( -tmp_357 * tmp_521 - tmp_370 * ( tmp_359 + tmp_362 + tmp_367 ) ) +
                     tmp_392 * ( -tmp_378 * tmp_521 - tmp_391 * ( tmp_380 + tmp_383 + tmp_388 ) ) +
                     tmp_413 * ( -tmp_399 * tmp_521 - tmp_412 * ( tmp_401 + tmp_404 + tmp_409 ) ) +
                     tmp_434 * ( -tmp_420 * tmp_521 - tmp_433 * ( tmp_422 + tmp_425 + tmp_430 ) ) +
                     tmp_455 * ( -tmp_441 * tmp_521 - tmp_454 * ( tmp_443 + tmp_446 + tmp_451 ) ) +
                     tmp_476 * ( -tmp_462 * tmp_521 - tmp_475 * ( tmp_464 + tmp_467 + tmp_472 ) ) +
                     tmp_497 * ( -tmp_483 * tmp_521 - tmp_496 * ( tmp_485 + tmp_488 + tmp_493 ) ) +
                     tmp_518 * ( -tmp_504 * tmp_521 - tmp_517 * ( tmp_506 + tmp_509 + tmp_514 ) ) +
                     tmp_98 * ( -tmp_46 * tmp_521 - tmp_96 * ( tmp_78 + tmp_82 + tmp_88 ) );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 0, 3 ) = a_0_3;
   }

   void integrateFacetDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                           const std::vector< Point3D >& coordsFacet,
                                           const Point3D&,
                                           const Point3D&     outwardNormal,
                                           const DGBasisInfo& trialBasis,
                                           const DGBasisInfo& testBasis,
                                           int                trialDegree,
                                           int                testDegree,
                                           MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_1_1 + tmp_0;
      real_t tmp_2  = -p_affine_0_2;
      real_t tmp_3  = p_affine_2_2 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_1 + tmp_0;
      real_t tmp_6  = p_affine_1_2 + tmp_2;
      real_t tmp_7  = tmp_5 * tmp_6;
      real_t tmp_8  = -p_affine_0_0;
      real_t tmp_9  = p_affine_1_0 + tmp_8;
      real_t tmp_10 = p_affine_3_2 + tmp_2;
      real_t tmp_11 = tmp_10 * tmp_5;
      real_t tmp_12 = p_affine_2_0 + tmp_8;
      real_t tmp_13 = p_affine_3_1 + tmp_0;
      real_t tmp_14 = tmp_13 * tmp_6;
      real_t tmp_15 = p_affine_3_0 + tmp_8;
      real_t tmp_16 = tmp_13 * tmp_3;
      real_t tmp_17 = tmp_1 * tmp_10;
      real_t tmp_18 =
          1.0 / ( tmp_11 * tmp_9 + tmp_12 * tmp_14 - tmp_12 * tmp_17 + tmp_15 * tmp_4 - tmp_15 * tmp_7 - tmp_16 * tmp_9 );
      real_t tmp_19 = tmp_18 * ( tmp_4 - tmp_7 );
      real_t tmp_20 = tmp_18 * ( tmp_14 - tmp_17 );
      real_t tmp_21 = tmp_18 * ( tmp_11 - tmp_16 );
      real_t tmp_22 = tmp_18 * ( tmp_12 * tmp_6 - tmp_3 * tmp_9 );
      real_t tmp_23 = tmp_18 * ( tmp_10 * tmp_9 - tmp_15 * tmp_6 );
      real_t tmp_24 = tmp_18 * ( -tmp_10 * tmp_12 + tmp_15 * tmp_3 );
      real_t tmp_25 = tmp_18 * ( -tmp_1 * tmp_12 + tmp_5 * tmp_9 );
      real_t tmp_26 = tmp_18 * ( tmp_1 * tmp_15 - tmp_13 * tmp_9 );
      real_t tmp_27 = tmp_18 * ( tmp_12 * tmp_13 - tmp_15 * tmp_5 );
      real_t tmp_28 = -p_affine_8_0;
      real_t tmp_29 = p_affine_10_0 + tmp_28;
      real_t tmp_30 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_31 = -p_affine_8_1;
      real_t tmp_32 = p_affine_10_1 + tmp_31;
      real_t tmp_33 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_34 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_35 = -p_affine_8_2;
      real_t tmp_36 = p_affine_10_2 + tmp_35;
      real_t tmp_37 =
          1.0 * std::pow( ( std::abs( tmp_29 * tmp_30 - tmp_32 * tmp_33 ) * std::abs( tmp_29 * tmp_30 - tmp_32 * tmp_33 ) ) +
                              ( std::abs( tmp_29 * tmp_34 - tmp_33 * tmp_36 ) * std::abs( tmp_29 * tmp_34 - tmp_33 * tmp_36 ) ) +
                              ( std::abs( tmp_30 * tmp_36 - tmp_32 * tmp_34 ) * std::abs( tmp_30 * tmp_36 - tmp_32 * tmp_34 ) ),
                          1.0 / 2.0 );
      real_t tmp_38 = tmp_37 * ( p_affine_13_0 * ( -tmp_19 - tmp_20 - tmp_21 ) + p_affine_13_1 * ( -tmp_22 - tmp_23 - tmp_24 ) +
                                 p_affine_13_2 * ( -tmp_25 - tmp_26 - tmp_27 ) );
      real_t tmp_39 = p_affine_9_2 + tmp_35;
      real_t tmp_40 = p_affine_8_2 + tmp_2;
      real_t tmp_41 = 0.93718850182767688 * tmp_36 + 0.031405749086161582 * tmp_39 + tmp_40;
      real_t tmp_42 = p_affine_9_1 + tmp_31;
      real_t tmp_43 = p_affine_8_1 + tmp_0;
      real_t tmp_44 = 0.93718850182767688 * tmp_32 + 0.031405749086161582 * tmp_42 + tmp_43;
      real_t tmp_45 = p_affine_9_0 + tmp_28;
      real_t tmp_46 = p_affine_8_0 + tmp_8;
      real_t tmp_47 = 0.93718850182767688 * tmp_29 + 0.031405749086161582 * tmp_45 + tmp_46;
      real_t tmp_48 = 0.0068572537431980923 * tmp_12 * ( tmp_20 * tmp_47 + tmp_23 * tmp_44 + tmp_26 * tmp_41 - 1.0 / 4.0 ) +
                      0.0068572537431980923 * tmp_15 * ( tmp_19 * tmp_47 + tmp_22 * tmp_44 + tmp_25 * tmp_41 - 1.0 / 4.0 ) +
                      0.0068572537431980923 * tmp_9 * ( tmp_21 * tmp_47 + tmp_24 * tmp_44 + tmp_27 * tmp_41 - 1.0 / 4.0 );
      real_t tmp_49 = 0.60796128279561268 * tmp_36 + 0.19601935860219369 * tmp_39 + tmp_40;
      real_t tmp_50 = 0.60796128279561268 * tmp_32 + 0.19601935860219369 * tmp_42 + tmp_43;
      real_t tmp_51 = 0.60796128279561268 * tmp_29 + 0.19601935860219369 * tmp_45 + tmp_46;
      real_t tmp_52 = 0.037198804536718075 * tmp_12 * ( tmp_20 * tmp_51 + tmp_23 * tmp_50 + tmp_26 * tmp_49 - 1.0 / 4.0 ) +
                      0.037198804536718075 * tmp_15 * ( tmp_19 * tmp_51 + tmp_22 * tmp_50 + tmp_25 * tmp_49 - 1.0 / 4.0 ) +
                      0.037198804536718075 * tmp_9 * ( tmp_21 * tmp_51 + tmp_24 * tmp_50 + tmp_27 * tmp_49 - 1.0 / 4.0 );
      real_t tmp_53 = 0.039308471900058539 * tmp_36 + 0.37605877282253791 * tmp_39 + tmp_40;
      real_t tmp_54 = 0.039308471900058539 * tmp_32 + 0.37605877282253791 * tmp_42 + tmp_43;
      real_t tmp_55 = 0.039308471900058539 * tmp_29 + 0.37605877282253791 * tmp_45 + tmp_46;
      real_t tmp_56 = 0.020848748529055869 * tmp_12 * ( tmp_20 * tmp_55 + tmp_23 * tmp_54 + tmp_26 * tmp_53 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_15 * ( tmp_19 * tmp_55 + tmp_22 * tmp_54 + tmp_25 * tmp_53 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_9 * ( tmp_21 * tmp_55 + tmp_24 * tmp_54 + tmp_27 * tmp_53 - 1.0 / 4.0 );
      real_t tmp_57 = 0.1711304259088916 * tmp_36 + 0.78764240869137092 * tmp_39 + tmp_40;
      real_t tmp_58 = 0.1711304259088916 * tmp_32 + 0.78764240869137092 * tmp_42 + tmp_43;
      real_t tmp_59 = 0.1711304259088916 * tmp_29 + 0.78764240869137092 * tmp_45 + tmp_46;
      real_t tmp_60 = 0.019202922745021479 * tmp_12 * ( tmp_20 * tmp_59 + tmp_23 * tmp_58 + tmp_26 * tmp_57 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_15 * ( tmp_19 * tmp_59 + tmp_22 * tmp_58 + tmp_25 * tmp_57 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_9 * ( tmp_21 * tmp_59 + tmp_24 * tmp_58 + tmp_27 * tmp_57 - 1.0 / 4.0 );
      real_t tmp_61 = 0.37605877282253791 * tmp_36 + 0.58463275527740355 * tmp_39 + tmp_40;
      real_t tmp_62 = 0.37605877282253791 * tmp_32 + 0.58463275527740355 * tmp_42 + tmp_43;
      real_t tmp_63 = 0.37605877282253791 * tmp_29 + 0.58463275527740355 * tmp_45 + tmp_46;
      real_t tmp_64 = 0.020848748529055869 * tmp_12 * ( tmp_20 * tmp_63 + tmp_23 * tmp_62 + tmp_26 * tmp_61 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_15 * ( tmp_19 * tmp_63 + tmp_22 * tmp_62 + tmp_25 * tmp_61 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_9 * ( tmp_21 * tmp_63 + tmp_24 * tmp_62 + tmp_27 * tmp_61 - 1.0 / 4.0 );
      real_t tmp_65 = 0.78764240869137092 * tmp_36 + 0.041227165399737475 * tmp_39 + tmp_40;
      real_t tmp_66 = 0.78764240869137092 * tmp_32 + 0.041227165399737475 * tmp_42 + tmp_43;
      real_t tmp_67 = 0.78764240869137092 * tmp_29 + 0.041227165399737475 * tmp_45 + tmp_46;
      real_t tmp_68 = 0.019202922745021479 * tmp_12 * ( tmp_20 * tmp_67 + tmp_23 * tmp_66 + tmp_26 * tmp_65 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_15 * ( tmp_19 * tmp_67 + tmp_22 * tmp_66 + tmp_25 * tmp_65 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_9 * ( tmp_21 * tmp_67 + tmp_24 * tmp_66 + tmp_27 * tmp_65 - 1.0 / 4.0 );
      real_t tmp_69 = 0.58463275527740355 * tmp_36 + 0.039308471900058539 * tmp_39 + tmp_40;
      real_t tmp_70 = 0.58463275527740355 * tmp_32 + 0.039308471900058539 * tmp_42 + tmp_43;
      real_t tmp_71 = 0.58463275527740355 * tmp_29 + 0.039308471900058539 * tmp_45 + tmp_46;
      real_t tmp_72 = 0.020848748529055869 * tmp_12 * ( tmp_20 * tmp_71 + tmp_23 * tmp_70 + tmp_26 * tmp_69 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_15 * ( tmp_19 * tmp_71 + tmp_22 * tmp_70 + tmp_25 * tmp_69 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_9 * ( tmp_21 * tmp_71 + tmp_24 * tmp_70 + tmp_27 * tmp_69 - 1.0 / 4.0 );
      real_t tmp_73 = 0.041227165399737475 * tmp_36 + 0.78764240869137092 * tmp_39 + tmp_40;
      real_t tmp_74 = 0.041227165399737475 * tmp_32 + 0.78764240869137092 * tmp_42 + tmp_43;
      real_t tmp_75 = 0.041227165399737475 * tmp_29 + 0.78764240869137092 * tmp_45 + tmp_46;
      real_t tmp_76 = 0.019202922745021479 * tmp_12 * ( tmp_20 * tmp_75 + tmp_23 * tmp_74 + tmp_26 * tmp_73 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_15 * ( tmp_19 * tmp_75 + tmp_22 * tmp_74 + tmp_25 * tmp_73 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_9 * ( tmp_21 * tmp_75 + tmp_24 * tmp_74 + tmp_27 * tmp_73 - 1.0 / 4.0 );
      real_t tmp_77 = 0.039308471900058539 * tmp_36 + 0.58463275527740355 * tmp_39 + tmp_40;
      real_t tmp_78 = 0.039308471900058539 * tmp_32 + 0.58463275527740355 * tmp_42 + tmp_43;
      real_t tmp_79 = 0.039308471900058539 * tmp_29 + 0.58463275527740355 * tmp_45 + tmp_46;
      real_t tmp_80 = 0.020848748529055869 * tmp_12 * ( tmp_20 * tmp_79 + tmp_23 * tmp_78 + tmp_26 * tmp_77 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_15 * ( tmp_19 * tmp_79 + tmp_22 * tmp_78 + tmp_25 * tmp_77 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_9 * ( tmp_21 * tmp_79 + tmp_24 * tmp_78 + tmp_27 * tmp_77 - 1.0 / 4.0 );
      real_t tmp_81 = 0.78764240869137092 * tmp_36 + 0.1711304259088916 * tmp_39 + tmp_40;
      real_t tmp_82 = 0.78764240869137092 * tmp_32 + 0.1711304259088916 * tmp_42 + tmp_43;
      real_t tmp_83 = 0.78764240869137092 * tmp_29 + 0.1711304259088916 * tmp_45 + tmp_46;
      real_t tmp_84 = 0.019202922745021479 * tmp_12 * ( tmp_20 * tmp_83 + tmp_23 * tmp_82 + tmp_26 * tmp_81 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_15 * ( tmp_19 * tmp_83 + tmp_22 * tmp_82 + tmp_25 * tmp_81 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_9 * ( tmp_21 * tmp_83 + tmp_24 * tmp_82 + tmp_27 * tmp_81 - 1.0 / 4.0 );
      real_t tmp_85 = 0.58463275527740355 * tmp_36 + 0.37605877282253791 * tmp_39 + tmp_40;
      real_t tmp_86 = 0.58463275527740355 * tmp_32 + 0.37605877282253791 * tmp_42 + tmp_43;
      real_t tmp_87 = 0.58463275527740355 * tmp_29 + 0.37605877282253791 * tmp_45 + tmp_46;
      real_t tmp_88 = 0.020848748529055869 * tmp_12 * ( tmp_20 * tmp_87 + tmp_23 * tmp_86 + tmp_26 * tmp_85 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_15 * ( tmp_19 * tmp_87 + tmp_22 * tmp_86 + tmp_25 * tmp_85 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_9 * ( tmp_21 * tmp_87 + tmp_24 * tmp_86 + tmp_27 * tmp_85 - 1.0 / 4.0 );
      real_t tmp_89 = 0.1711304259088916 * tmp_36 + 0.041227165399737475 * tmp_39 + tmp_40;
      real_t tmp_90 = 0.1711304259088916 * tmp_32 + 0.041227165399737475 * tmp_42 + tmp_43;
      real_t tmp_91 = 0.1711304259088916 * tmp_29 + 0.041227165399737475 * tmp_45 + tmp_46;
      real_t tmp_92 = 0.019202922745021479 * tmp_12 * ( tmp_20 * tmp_91 + tmp_23 * tmp_90 + tmp_26 * tmp_89 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_15 * ( tmp_19 * tmp_91 + tmp_22 * tmp_90 + tmp_25 * tmp_89 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_9 * ( tmp_21 * tmp_91 + tmp_24 * tmp_90 + tmp_27 * tmp_89 - 1.0 / 4.0 );
      real_t tmp_93 = 0.19107600050469298 * tmp_36 + 0.40446199974765351 * tmp_39 + tmp_40;
      real_t tmp_94 = 0.19107600050469298 * tmp_32 + 0.40446199974765351 * tmp_42 + tmp_43;
      real_t tmp_95 = 0.19107600050469298 * tmp_29 + 0.40446199974765351 * tmp_45 + tmp_46;
      real_t tmp_96 = 0.042507265838595799 * tmp_12 * ( tmp_20 * tmp_95 + tmp_23 * tmp_94 + tmp_26 * tmp_93 - 1.0 / 4.0 ) +
                      0.042507265838595799 * tmp_15 * ( tmp_19 * tmp_95 + tmp_22 * tmp_94 + tmp_25 * tmp_93 - 1.0 / 4.0 ) +
                      0.042507265838595799 * tmp_9 * ( tmp_21 * tmp_95 + tmp_24 * tmp_94 + tmp_27 * tmp_93 - 1.0 / 4.0 );
      real_t tmp_97  = 0.37605877282253791 * tmp_36 + 0.039308471900058539 * tmp_39 + tmp_40;
      real_t tmp_98  = 0.37605877282253791 * tmp_32 + 0.039308471900058539 * tmp_42 + tmp_43;
      real_t tmp_99  = 0.37605877282253791 * tmp_29 + 0.039308471900058539 * tmp_45 + tmp_46;
      real_t tmp_100 = 0.020848748529055869 * tmp_12 * ( tmp_20 * tmp_99 + tmp_23 * tmp_98 + tmp_26 * tmp_97 - 1.0 / 4.0 ) +
                       0.020848748529055869 * tmp_15 * ( tmp_19 * tmp_99 + tmp_22 * tmp_98 + tmp_25 * tmp_97 - 1.0 / 4.0 ) +
                       0.020848748529055869 * tmp_9 * ( tmp_21 * tmp_99 + tmp_24 * tmp_98 + tmp_27 * tmp_97 - 1.0 / 4.0 );
      real_t tmp_101 = 0.031405749086161582 * tmp_36 + 0.93718850182767688 * tmp_39 + tmp_40;
      real_t tmp_102 = 0.031405749086161582 * tmp_32 + 0.93718850182767688 * tmp_42 + tmp_43;
      real_t tmp_103 = 0.031405749086161582 * tmp_29 + 0.93718850182767688 * tmp_45 + tmp_46;
      real_t tmp_104 = 0.0068572537431980923 * tmp_12 * ( tmp_101 * tmp_26 + tmp_102 * tmp_23 + tmp_103 * tmp_20 - 1.0 / 4.0 ) +
                       0.0068572537431980923 * tmp_15 * ( tmp_101 * tmp_25 + tmp_102 * tmp_22 + tmp_103 * tmp_19 - 1.0 / 4.0 ) +
                       0.0068572537431980923 * tmp_9 * ( tmp_101 * tmp_27 + tmp_102 * tmp_24 + tmp_103 * tmp_21 - 1.0 / 4.0 );
      real_t tmp_105 = 0.19601935860219369 * tmp_36 + 0.60796128279561268 * tmp_39 + tmp_40;
      real_t tmp_106 = 0.19601935860219369 * tmp_32 + 0.60796128279561268 * tmp_42 + tmp_43;
      real_t tmp_107 = 0.19601935860219369 * tmp_29 + 0.60796128279561268 * tmp_45 + tmp_46;
      real_t tmp_108 = 0.037198804536718075 * tmp_12 * ( tmp_105 * tmp_26 + tmp_106 * tmp_23 + tmp_107 * tmp_20 - 1.0 / 4.0 ) +
                       0.037198804536718075 * tmp_15 * ( tmp_105 * tmp_25 + tmp_106 * tmp_22 + tmp_107 * tmp_19 - 1.0 / 4.0 ) +
                       0.037198804536718075 * tmp_9 * ( tmp_105 * tmp_27 + tmp_106 * tmp_24 + tmp_107 * tmp_21 - 1.0 / 4.0 );
      real_t tmp_109 = 0.40446199974765351 * tmp_36 + 0.19107600050469298 * tmp_39 + tmp_40;
      real_t tmp_110 = 0.40446199974765351 * tmp_32 + 0.19107600050469298 * tmp_42 + tmp_43;
      real_t tmp_111 = 0.40446199974765351 * tmp_29 + 0.19107600050469298 * tmp_45 + tmp_46;
      real_t tmp_112 = 0.042507265838595799 * tmp_12 * ( tmp_109 * tmp_26 + tmp_110 * tmp_23 + tmp_111 * tmp_20 - 1.0 / 4.0 ) +
                       0.042507265838595799 * tmp_15 * ( tmp_109 * tmp_25 + tmp_110 * tmp_22 + tmp_111 * tmp_19 - 1.0 / 4.0 ) +
                       0.042507265838595799 * tmp_9 * ( tmp_109 * tmp_27 + tmp_110 * tmp_24 + tmp_111 * tmp_21 - 1.0 / 4.0 );
      real_t tmp_113 = 0.031405749086161582 * tmp_36 + 0.031405749086161582 * tmp_39 + tmp_40;
      real_t tmp_114 = 0.031405749086161582 * tmp_32 + 0.031405749086161582 * tmp_42 + tmp_43;
      real_t tmp_115 = 0.031405749086161582 * tmp_29 + 0.031405749086161582 * tmp_45 + tmp_46;
      real_t tmp_116 = 0.0068572537431980923 * tmp_12 * ( tmp_113 * tmp_26 + tmp_114 * tmp_23 + tmp_115 * tmp_20 - 1.0 / 4.0 ) +
                       0.0068572537431980923 * tmp_15 * ( tmp_113 * tmp_25 + tmp_114 * tmp_22 + tmp_115 * tmp_19 - 1.0 / 4.0 ) +
                       0.0068572537431980923 * tmp_9 * ( tmp_113 * tmp_27 + tmp_114 * tmp_24 + tmp_115 * tmp_21 - 1.0 / 4.0 );
      real_t tmp_117 = 0.19601935860219369 * tmp_36 + 0.19601935860219369 * tmp_39 + tmp_40;
      real_t tmp_118 = 0.19601935860219369 * tmp_32 + 0.19601935860219369 * tmp_42 + tmp_43;
      real_t tmp_119 = 0.19601935860219369 * tmp_29 + 0.19601935860219369 * tmp_45 + tmp_46;
      real_t tmp_120 = 0.037198804536718075 * tmp_12 * ( tmp_117 * tmp_26 + tmp_118 * tmp_23 + tmp_119 * tmp_20 - 1.0 / 4.0 ) +
                       0.037198804536718075 * tmp_15 * ( tmp_117 * tmp_25 + tmp_118 * tmp_22 + tmp_119 * tmp_19 - 1.0 / 4.0 ) +
                       0.037198804536718075 * tmp_9 * ( tmp_117 * tmp_27 + tmp_118 * tmp_24 + tmp_119 * tmp_21 - 1.0 / 4.0 );
      real_t tmp_121 = 0.40446199974765351 * tmp_36 + 0.40446199974765351 * tmp_39 + tmp_40;
      real_t tmp_122 = 0.40446199974765351 * tmp_32 + 0.40446199974765351 * tmp_42 + tmp_43;
      real_t tmp_123 = 0.40446199974765351 * tmp_29 + 0.40446199974765351 * tmp_45 + tmp_46;
      real_t tmp_124 = 0.042507265838595799 * tmp_12 * ( tmp_121 * tmp_26 + tmp_122 * tmp_23 + tmp_123 * tmp_20 - 1.0 / 4.0 ) +
                       0.042507265838595799 * tmp_15 * ( tmp_121 * tmp_25 + tmp_122 * tmp_22 + tmp_123 * tmp_19 - 1.0 / 4.0 ) +
                       0.042507265838595799 * tmp_9 * ( tmp_121 * tmp_27 + tmp_122 * tmp_24 + tmp_123 * tmp_21 - 1.0 / 4.0 );
      real_t tmp_125 = 0.041227165399737475 * tmp_36 + 0.1711304259088916 * tmp_39 + tmp_40;
      real_t tmp_126 = 0.041227165399737475 * tmp_32 + 0.1711304259088916 * tmp_42 + tmp_43;
      real_t tmp_127 = 0.041227165399737475 * tmp_29 + 0.1711304259088916 * tmp_45 + tmp_46;
      real_t tmp_128 = 0.019202922745021479 * tmp_12 * ( tmp_125 * tmp_26 + tmp_126 * tmp_23 + tmp_127 * tmp_20 - 1.0 / 4.0 ) +
                       0.019202922745021479 * tmp_15 * ( tmp_125 * tmp_25 + tmp_126 * tmp_22 + tmp_127 * tmp_19 - 1.0 / 4.0 ) +
                       0.019202922745021479 * tmp_9 * ( tmp_125 * tmp_27 + tmp_126 * tmp_24 + tmp_127 * tmp_21 - 1.0 / 4.0 );
      real_t tmp_129 = tmp_37 * ( p_affine_13_0 * tmp_21 + p_affine_13_1 * tmp_24 + p_affine_13_2 * tmp_27 );
      real_t tmp_130 = tmp_37 * ( p_affine_13_0 * tmp_20 + p_affine_13_1 * tmp_23 + p_affine_13_2 * tmp_26 );
      real_t tmp_131 = tmp_37 * ( p_affine_13_0 * tmp_19 + p_affine_13_1 * tmp_22 + p_affine_13_2 * tmp_25 );
      real_t a_0_0   = -tmp_100 * tmp_38 - tmp_104 * tmp_38 - tmp_108 * tmp_38 - tmp_112 * tmp_38 - tmp_116 * tmp_38 -
                     tmp_120 * tmp_38 - tmp_124 * tmp_38 - tmp_128 * tmp_38 - tmp_38 * tmp_48 - tmp_38 * tmp_52 -
                     tmp_38 * tmp_56 - tmp_38 * tmp_60 - tmp_38 * tmp_64 - tmp_38 * tmp_68 - tmp_38 * tmp_72 - tmp_38 * tmp_76 -
                     tmp_38 * tmp_80 - tmp_38 * tmp_84 - tmp_38 * tmp_88 - tmp_38 * tmp_92 - tmp_38 * tmp_96;
      real_t a_0_1 = -tmp_100 * tmp_129 - tmp_104 * tmp_129 - tmp_108 * tmp_129 - tmp_112 * tmp_129 - tmp_116 * tmp_129 -
                     tmp_120 * tmp_129 - tmp_124 * tmp_129 - tmp_128 * tmp_129 - tmp_129 * tmp_48 - tmp_129 * tmp_52 -
                     tmp_129 * tmp_56 - tmp_129 * tmp_60 - tmp_129 * tmp_64 - tmp_129 * tmp_68 - tmp_129 * tmp_72 -
                     tmp_129 * tmp_76 - tmp_129 * tmp_80 - tmp_129 * tmp_84 - tmp_129 * tmp_88 - tmp_129 * tmp_92 -
                     tmp_129 * tmp_96;
      real_t a_0_2 = -tmp_100 * tmp_130 - tmp_104 * tmp_130 - tmp_108 * tmp_130 - tmp_112 * tmp_130 - tmp_116 * tmp_130 -
                     tmp_120 * tmp_130 - tmp_124 * tmp_130 - tmp_128 * tmp_130 - tmp_130 * tmp_48 - tmp_130 * tmp_52 -
                     tmp_130 * tmp_56 - tmp_130 * tmp_60 - tmp_130 * tmp_64 - tmp_130 * tmp_68 - tmp_130 * tmp_72 -
                     tmp_130 * tmp_76 - tmp_130 * tmp_80 - tmp_130 * tmp_84 - tmp_130 * tmp_88 - tmp_130 * tmp_92 -
                     tmp_130 * tmp_96;
      real_t a_0_3 = -tmp_100 * tmp_131 - tmp_104 * tmp_131 - tmp_108 * tmp_131 - tmp_112 * tmp_131 - tmp_116 * tmp_131 -
                     tmp_120 * tmp_131 - tmp_124 * tmp_131 - tmp_128 * tmp_131 - tmp_131 * tmp_48 - tmp_131 * tmp_52 -
                     tmp_131 * tmp_56 - tmp_131 * tmp_60 - tmp_131 * tmp_64 - tmp_131 * tmp_68 - tmp_131 * tmp_72 -
                     tmp_131 * tmp_76 - tmp_131 * tmp_80 - tmp_131 * tmp_84 - tmp_131 * tmp_88 - tmp_131 * tmp_92 -
                     tmp_131 * tmp_96;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 0, 3 ) = a_0_3;
   }
};

class EGIIPGVectorLaplaceFormEP1_1 : public hyteg::dg::DGForm
{
 protected:
   void integrateVolume2D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override
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

   virtual void integrateFacetInner2D( const std::vector< Point3D >& coordsElement,
                                       const std::vector< Point3D >& coordsFacet,
                                       const Point3D&                oppositeVertex,
                                       const Point3D&                outwardNormal,
                                       const DGBasisInfo&            trialBasis,
                                       const DGBasisInfo&            testBasis,
                                       int                           trialDegree,
                                       int                           testDegree,
                                       MatrixXr&                     elMat ) const override
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
      real_t tmp_4  = 0.046910077030668018 * tmp_2 + tmp_3;
      real_t tmp_5  = -p_affine_0_0;
      real_t tmp_6  = p_affine_1_0 + tmp_5;
      real_t tmp_7  = p_affine_2_1 + tmp_0;
      real_t tmp_8  = 1.0 / ( -tmp_1 * ( p_affine_2_0 + tmp_5 ) + tmp_6 * tmp_7 );
      real_t tmp_9  = tmp_8 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_10 = tmp_4 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_5;
      real_t tmp_13 = 0.046910077030668018 * tmp_11 + tmp_12;
      real_t tmp_14 = tmp_7 * tmp_8;
      real_t tmp_15 = tmp_13 * tmp_14;
      real_t tmp_16 = tmp_10 + tmp_15;
      real_t tmp_17 = tmp_6 * tmp_8;
      real_t tmp_18 = tmp_17 * tmp_4;
      real_t tmp_19 = tmp_8 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_20 = tmp_13 * tmp_19;
      real_t tmp_21 = tmp_18 + tmp_20;
      real_t tmp_22 = tmp_1 * ( tmp_16 - 1.0 / 3.0 ) + tmp_7 * ( tmp_21 - 1.0 / 3.0 );
      real_t tmp_23 = 0.5 * p_affine_10_0 * ( -tmp_14 - tmp_19 ) + 0.5 * p_affine_10_1 * ( -tmp_17 - tmp_9 );
      real_t tmp_24 = std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_2 * tmp_2 ), 1.0 / 2.0 ) );
      real_t tmp_25 = 1.0 / ( tmp_24 );
      real_t tmp_26 = tmp_22 * tmp_25;
      real_t tmp_27 = 0.11846344252809471 * tmp_24;
      real_t tmp_28 = 0.23076534494715845 * tmp_2 + tmp_3;
      real_t tmp_29 = tmp_28 * tmp_9;
      real_t tmp_30 = 0.23076534494715845 * tmp_11 + tmp_12;
      real_t tmp_31 = tmp_14 * tmp_30;
      real_t tmp_32 = tmp_29 + tmp_31;
      real_t tmp_33 = tmp_17 * tmp_28;
      real_t tmp_34 = tmp_19 * tmp_30;
      real_t tmp_35 = tmp_33 + tmp_34;
      real_t tmp_36 = tmp_1 * ( tmp_32 - 1.0 / 3.0 ) + tmp_7 * ( tmp_35 - 1.0 / 3.0 );
      real_t tmp_37 = tmp_25 * tmp_36;
      real_t tmp_38 = 0.2393143352496831 * tmp_24;
      real_t tmp_39 = 0.5 * tmp_2 + tmp_3;
      real_t tmp_40 = tmp_39 * tmp_9;
      real_t tmp_41 = 0.5 * tmp_11 + tmp_12;
      real_t tmp_42 = tmp_14 * tmp_41;
      real_t tmp_43 = tmp_40 + tmp_42;
      real_t tmp_44 = tmp_17 * tmp_39;
      real_t tmp_45 = tmp_19 * tmp_41;
      real_t tmp_46 = tmp_44 + tmp_45;
      real_t tmp_47 = tmp_1 * ( tmp_43 - 1.0 / 3.0 ) + tmp_7 * ( tmp_46 - 1.0 / 3.0 );
      real_t tmp_48 = tmp_25 * tmp_47;
      real_t tmp_49 = 0.2844444444444445 * tmp_24;
      real_t tmp_50 = 0.7692346550528415 * tmp_2 + tmp_3;
      real_t tmp_51 = tmp_50 * tmp_9;
      real_t tmp_52 = 0.7692346550528415 * tmp_11 + tmp_12;
      real_t tmp_53 = tmp_14 * tmp_52;
      real_t tmp_54 = tmp_51 + tmp_53;
      real_t tmp_55 = tmp_17 * tmp_50;
      real_t tmp_56 = tmp_19 * tmp_52;
      real_t tmp_57 = tmp_55 + tmp_56;
      real_t tmp_58 = tmp_1 * ( tmp_54 - 1.0 / 3.0 ) + tmp_7 * ( tmp_57 - 1.0 / 3.0 );
      real_t tmp_59 = tmp_25 * tmp_58;
      real_t tmp_60 = 0.2393143352496831 * tmp_24;
      real_t tmp_61 = 0.95308992296933193 * tmp_2 + tmp_3;
      real_t tmp_62 = tmp_61 * tmp_9;
      real_t tmp_63 = 0.95308992296933193 * tmp_11 + tmp_12;
      real_t tmp_64 = tmp_14 * tmp_63;
      real_t tmp_65 = tmp_62 + tmp_64;
      real_t tmp_66 = tmp_17 * tmp_61;
      real_t tmp_67 = tmp_19 * tmp_63;
      real_t tmp_68 = tmp_66 + tmp_67;
      real_t tmp_69 = tmp_1 * ( tmp_65 - 1.0 / 3.0 ) + tmp_7 * ( tmp_68 - 1.0 / 3.0 );
      real_t tmp_70 = tmp_25 * tmp_69;
      real_t tmp_71 = 0.11846344252809471 * tmp_24;
      real_t tmp_72 = 0.5 * p_affine_10_0 * tmp_14 + 0.5 * p_affine_10_1 * tmp_9;
      real_t tmp_73 = 0.5 * p_affine_10_0 * tmp_19 + 0.5 * p_affine_10_1 * tmp_17;
      real_t a_0_0  = tmp_27 * ( -tmp_22 * tmp_23 + tmp_26 * ( -tmp_10 - tmp_15 - tmp_18 - tmp_20 + 1 ) ) +
                     tmp_38 * ( -tmp_23 * tmp_36 + tmp_37 * ( -tmp_29 - tmp_31 - tmp_33 - tmp_34 + 1 ) ) +
                     tmp_49 * ( -tmp_23 * tmp_47 + tmp_48 * ( -tmp_40 - tmp_42 - tmp_44 - tmp_45 + 1 ) ) +
                     tmp_60 * ( -tmp_23 * tmp_58 + tmp_59 * ( -tmp_51 - tmp_53 - tmp_55 - tmp_56 + 1 ) ) +
                     tmp_71 * ( -tmp_23 * tmp_69 + tmp_70 * ( -tmp_62 - tmp_64 - tmp_66 - tmp_67 + 1 ) );
      real_t a_0_1 = tmp_27 * ( tmp_16 * tmp_26 - tmp_22 * tmp_72 ) + tmp_38 * ( tmp_32 * tmp_37 - tmp_36 * tmp_72 ) +
                     tmp_49 * ( tmp_43 * tmp_48 - tmp_47 * tmp_72 ) + tmp_60 * ( tmp_54 * tmp_59 - tmp_58 * tmp_72 ) +
                     tmp_71 * ( tmp_65 * tmp_70 - tmp_69 * tmp_72 );
      real_t a_0_2 = tmp_27 * ( tmp_21 * tmp_26 - tmp_22 * tmp_73 ) + tmp_38 * ( tmp_35 * tmp_37 - tmp_36 * tmp_73 ) +
                     tmp_49 * ( tmp_46 * tmp_48 - tmp_47 * tmp_73 ) + tmp_60 * ( tmp_57 * tmp_59 - tmp_58 * tmp_73 ) +
                     tmp_71 * ( tmp_68 * tmp_70 - tmp_69 * tmp_73 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >& coordsElementInner,
                                          const std::vector< Point3D >& coordsElementOuter,
                                          const std::vector< Point3D >& coordsFacet,
                                          const Point3D&                oppositeVertexInnerElement,
                                          const Point3D&                oppositeVertexOuterElement,
                                          const Point3D&                outwardNormal,
                                          const DGBasisInfo&            trialBasis,
                                          const DGBasisInfo&            testBasis,
                                          int                           trialDegree,
                                          int                           testDegree,
                                          MatrixXr&                     elMat ) const override
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
      real_t tmp_2  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3  = -p_affine_0_0;
      real_t tmp_4  = p_affine_1_0 + tmp_3;
      real_t tmp_5  = p_affine_2_1 + tmp_0;
      real_t tmp_6  = 1.0 / ( -tmp_1 * ( p_affine_2_0 + tmp_3 ) + tmp_4 * tmp_5 );
      real_t tmp_7  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_8  = p_affine_6_1 + 0.046910077030668018 * tmp_7;
      real_t tmp_9  = tmp_6 * ( tmp_0 + tmp_8 );
      real_t tmp_10 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_11 = p_affine_6_0 + 0.046910077030668018 * tmp_10;
      real_t tmp_12 = tmp_6 * ( tmp_11 + tmp_3 );
      real_t tmp_13 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_14 =
          tmp_1 * ( tmp_12 * tmp_5 + tmp_2 * tmp_9 - 1.0 / 3.0 ) + tmp_5 * ( tmp_12 * tmp_13 + tmp_4 * tmp_9 - 1.0 / 3.0 );
      real_t tmp_15 = -p_affine_3_1;
      real_t tmp_16 = p_affine_5_1 + tmp_15;
      real_t tmp_17 = -p_affine_3_0;
      real_t tmp_18 = p_affine_4_0 + tmp_17;
      real_t tmp_19 = 1.0 / ( tmp_16 * tmp_18 - ( p_affine_4_1 + tmp_15 ) * ( p_affine_5_0 + tmp_17 ) );
      real_t tmp_20 = tmp_16 * tmp_19;
      real_t tmp_21 = tmp_19 * ( p_affine_3_1 - p_affine_4_1 );
      real_t tmp_22 = tmp_18 * tmp_19;
      real_t tmp_23 = tmp_19 * ( p_affine_3_0 - p_affine_5_0 );
      real_t tmp_24 = 0.5 * p_affine_10_0 * ( -tmp_20 - tmp_21 ) + 0.5 * p_affine_10_1 * ( -tmp_22 - tmp_23 );
      real_t tmp_25 = tmp_15 + tmp_8;
      real_t tmp_26 = tmp_22 * tmp_25;
      real_t tmp_27 = tmp_23 * tmp_25;
      real_t tmp_28 = tmp_11 + tmp_17;
      real_t tmp_29 = tmp_20 * tmp_28;
      real_t tmp_30 = tmp_21 * tmp_28;
      real_t tmp_31 = std::abs( std::pow( ( tmp_10 * tmp_10 ) + ( tmp_7 * tmp_7 ), 1.0 / 2.0 ) );
      real_t tmp_32 = 1.0 / ( tmp_31 );
      real_t tmp_33 = tmp_14 * tmp_32;
      real_t tmp_34 = 0.11846344252809471 * tmp_31;
      real_t tmp_35 = p_affine_6_1 + 0.23076534494715845 * tmp_7;
      real_t tmp_36 = tmp_6 * ( tmp_0 + tmp_35 );
      real_t tmp_37 = p_affine_6_0 + 0.23076534494715845 * tmp_10;
      real_t tmp_38 = tmp_6 * ( tmp_3 + tmp_37 );
      real_t tmp_39 =
          tmp_1 * ( tmp_2 * tmp_36 + tmp_38 * tmp_5 - 1.0 / 3.0 ) + tmp_5 * ( tmp_13 * tmp_38 + tmp_36 * tmp_4 - 1.0 / 3.0 );
      real_t tmp_40 = tmp_15 + tmp_35;
      real_t tmp_41 = tmp_22 * tmp_40;
      real_t tmp_42 = tmp_23 * tmp_40;
      real_t tmp_43 = tmp_17 + tmp_37;
      real_t tmp_44 = tmp_20 * tmp_43;
      real_t tmp_45 = tmp_21 * tmp_43;
      real_t tmp_46 = tmp_32 * tmp_39;
      real_t tmp_47 = 0.2393143352496831 * tmp_31;
      real_t tmp_48 = p_affine_6_1 + 0.5 * tmp_7;
      real_t tmp_49 = tmp_6 * ( tmp_0 + tmp_48 );
      real_t tmp_50 = p_affine_6_0 + 0.5 * tmp_10;
      real_t tmp_51 = tmp_6 * ( tmp_3 + tmp_50 );
      real_t tmp_52 =
          tmp_1 * ( tmp_2 * tmp_49 + tmp_5 * tmp_51 - 1.0 / 3.0 ) + tmp_5 * ( tmp_13 * tmp_51 + tmp_4 * tmp_49 - 1.0 / 3.0 );
      real_t tmp_53 = tmp_15 + tmp_48;
      real_t tmp_54 = tmp_22 * tmp_53;
      real_t tmp_55 = tmp_23 * tmp_53;
      real_t tmp_56 = tmp_17 + tmp_50;
      real_t tmp_57 = tmp_20 * tmp_56;
      real_t tmp_58 = tmp_21 * tmp_56;
      real_t tmp_59 = tmp_32 * tmp_52;
      real_t tmp_60 = 0.2844444444444445 * tmp_31;
      real_t tmp_61 = p_affine_6_1 + 0.7692346550528415 * tmp_7;
      real_t tmp_62 = tmp_6 * ( tmp_0 + tmp_61 );
      real_t tmp_63 = p_affine_6_0 + 0.7692346550528415 * tmp_10;
      real_t tmp_64 = tmp_6 * ( tmp_3 + tmp_63 );
      real_t tmp_65 =
          tmp_1 * ( tmp_2 * tmp_62 + tmp_5 * tmp_64 - 1.0 / 3.0 ) + tmp_5 * ( tmp_13 * tmp_64 + tmp_4 * tmp_62 - 1.0 / 3.0 );
      real_t tmp_66 = tmp_15 + tmp_61;
      real_t tmp_67 = tmp_22 * tmp_66;
      real_t tmp_68 = tmp_23 * tmp_66;
      real_t tmp_69 = tmp_17 + tmp_63;
      real_t tmp_70 = tmp_20 * tmp_69;
      real_t tmp_71 = tmp_21 * tmp_69;
      real_t tmp_72 = tmp_32 * tmp_65;
      real_t tmp_73 = 0.2393143352496831 * tmp_31;
      real_t tmp_74 = p_affine_6_1 + 0.95308992296933193 * tmp_7;
      real_t tmp_75 = tmp_6 * ( tmp_0 + tmp_74 );
      real_t tmp_76 = p_affine_6_0 + 0.95308992296933193 * tmp_10;
      real_t tmp_77 = tmp_6 * ( tmp_3 + tmp_76 );
      real_t tmp_78 =
          tmp_1 * ( tmp_2 * tmp_75 + tmp_5 * tmp_77 - 1.0 / 3.0 ) + tmp_5 * ( tmp_13 * tmp_77 + tmp_4 * tmp_75 - 1.0 / 3.0 );
      real_t tmp_79 = tmp_15 + tmp_74;
      real_t tmp_80 = tmp_22 * tmp_79;
      real_t tmp_81 = tmp_23 * tmp_79;
      real_t tmp_82 = tmp_17 + tmp_76;
      real_t tmp_83 = tmp_20 * tmp_82;
      real_t tmp_84 = tmp_21 * tmp_82;
      real_t tmp_85 = tmp_32 * tmp_78;
      real_t tmp_86 = 0.11846344252809471 * tmp_31;
      real_t tmp_87 = 0.5 * p_affine_10_0 * tmp_20 + 0.5 * p_affine_10_1 * tmp_23;
      real_t tmp_88 = 0.5 * p_affine_10_0 * tmp_21 + 0.5 * p_affine_10_1 * tmp_22;
      real_t a_0_0  = tmp_34 * ( -tmp_14 * tmp_24 - tmp_33 * ( -tmp_26 - tmp_27 - tmp_29 - tmp_30 + 1 ) ) +
                     tmp_47 * ( -tmp_24 * tmp_39 - tmp_46 * ( -tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1 ) ) +
                     tmp_60 * ( -tmp_24 * tmp_52 - tmp_59 * ( -tmp_54 - tmp_55 - tmp_57 - tmp_58 + 1 ) ) +
                     tmp_73 * ( -tmp_24 * tmp_65 - tmp_72 * ( -tmp_67 - tmp_68 - tmp_70 - tmp_71 + 1 ) ) +
                     tmp_86 * ( -tmp_24 * tmp_78 - tmp_85 * ( -tmp_80 - tmp_81 - tmp_83 - tmp_84 + 1 ) );
      real_t a_0_1 = tmp_34 * ( -tmp_14 * tmp_87 - tmp_33 * ( tmp_27 + tmp_29 ) ) +
                     tmp_47 * ( -tmp_39 * tmp_87 - tmp_46 * ( tmp_42 + tmp_44 ) ) +
                     tmp_60 * ( -tmp_52 * tmp_87 - tmp_59 * ( tmp_55 + tmp_57 ) ) +
                     tmp_73 * ( -tmp_65 * tmp_87 - tmp_72 * ( tmp_68 + tmp_70 ) ) +
                     tmp_86 * ( -tmp_78 * tmp_87 - tmp_85 * ( tmp_81 + tmp_83 ) );
      real_t a_0_2 = tmp_34 * ( -tmp_14 * tmp_88 - tmp_33 * ( tmp_26 + tmp_30 ) ) +
                     tmp_47 * ( -tmp_39 * tmp_88 - tmp_46 * ( tmp_41 + tmp_45 ) ) +
                     tmp_60 * ( -tmp_52 * tmp_88 - tmp_59 * ( tmp_54 + tmp_58 ) ) +
                     tmp_73 * ( -tmp_65 * tmp_88 - tmp_72 * ( tmp_67 + tmp_71 ) ) +
                     tmp_86 * ( -tmp_78 * tmp_88 - tmp_85 * ( tmp_80 + tmp_84 ) );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                                   const std::vector< Point3D >& coordsFacet,
                                                   const Point3D&                oppositeVertex,
                                                   const Point3D&                outwardNormal,
                                                   const DGBasisInfo&            trialBasis,
                                                   const DGBasisInfo&            testBasis,
                                                   int                           trialDegree,
                                                   int                           testDegree,
                                                   MatrixXr&                     elMat ) const override
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
      real_t tmp_4  = p_affine_1_1 + tmp_0;
      real_t tmp_5  = 1.0 / ( tmp_1 * tmp_3 - tmp_4 * ( p_affine_2_0 + tmp_2 ) );
      real_t tmp_6  = tmp_1 * tmp_5;
      real_t tmp_7  = tmp_5 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_8  = tmp_3 * tmp_5;
      real_t tmp_9  = tmp_5 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_10 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_11 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_12 = std::abs( std::pow( ( tmp_10 * tmp_10 ) + ( tmp_11 * tmp_11 ), 1.0 / 2.0 ) );
      real_t tmp_13 = tmp_12 * ( p_affine_10_0 * ( -tmp_6 - tmp_7 ) + p_affine_10_1 * ( -tmp_8 - tmp_9 ) );
      real_t tmp_14 = p_affine_6_1 + tmp_0;
      real_t tmp_15 = 0.046910077030668018 * tmp_11 + tmp_14;
      real_t tmp_16 = p_affine_6_0 + tmp_2;
      real_t tmp_17 = 0.046910077030668018 * tmp_10 + tmp_16;
      real_t tmp_18 = 0.11846344252809471 * tmp_1 * ( tmp_15 * tmp_8 + tmp_17 * tmp_7 - 1.0 / 3.0 ) +
                      0.11846344252809471 * tmp_4 * ( tmp_15 * tmp_9 + tmp_17 * tmp_6 - 1.0 / 3.0 );
      real_t tmp_19 = 0.23076534494715845 * tmp_11 + tmp_14;
      real_t tmp_20 = 0.23076534494715845 * tmp_10 + tmp_16;
      real_t tmp_21 = 0.2393143352496831 * tmp_1 * ( tmp_19 * tmp_8 + tmp_20 * tmp_7 - 1.0 / 3.0 ) +
                      0.2393143352496831 * tmp_4 * ( tmp_19 * tmp_9 + tmp_20 * tmp_6 - 1.0 / 3.0 );
      real_t tmp_22 = 0.5 * tmp_11 + tmp_14;
      real_t tmp_23 = 0.5 * tmp_10 + tmp_16;
      real_t tmp_24 = 0.2844444444444445 * tmp_1 * ( tmp_22 * tmp_8 + tmp_23 * tmp_7 - 1.0 / 3.0 ) +
                      0.2844444444444445 * tmp_4 * ( tmp_22 * tmp_9 + tmp_23 * tmp_6 - 1.0 / 3.0 );
      real_t tmp_25 = 0.7692346550528415 * tmp_11 + tmp_14;
      real_t tmp_26 = 0.7692346550528415 * tmp_10 + tmp_16;
      real_t tmp_27 = 0.2393143352496831 * tmp_1 * ( tmp_25 * tmp_8 + tmp_26 * tmp_7 - 1.0 / 3.0 ) +
                      0.2393143352496831 * tmp_4 * ( tmp_25 * tmp_9 + tmp_26 * tmp_6 - 1.0 / 3.0 );
      real_t tmp_28 = 0.95308992296933193 * tmp_11 + tmp_14;
      real_t tmp_29 = 0.95308992296933193 * tmp_10 + tmp_16;
      real_t tmp_30 = 0.11846344252809471 * tmp_1 * ( tmp_28 * tmp_8 + tmp_29 * tmp_7 - 1.0 / 3.0 ) +
                      0.11846344252809471 * tmp_4 * ( tmp_28 * tmp_9 + tmp_29 * tmp_6 - 1.0 / 3.0 );
      real_t tmp_31 = tmp_12 * ( p_affine_10_0 * tmp_6 + p_affine_10_1 * tmp_9 );
      real_t tmp_32 = tmp_12 * ( p_affine_10_0 * tmp_7 + p_affine_10_1 * tmp_8 );
      real_t a_0_0  = -tmp_13 * tmp_18 - tmp_13 * tmp_21 - tmp_13 * tmp_24 - tmp_13 * tmp_27 - tmp_13 * tmp_30;
      real_t a_0_1  = -tmp_18 * tmp_31 - tmp_21 * tmp_31 - tmp_24 * tmp_31 - tmp_27 * tmp_31 - tmp_30 * tmp_31;
      real_t a_0_2  = -tmp_18 * tmp_32 - tmp_21 * tmp_32 - tmp_24 * tmp_32 - tmp_27 * tmp_32 - tmp_30 * tmp_32;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
   }

   void integrateRHSDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                         const std::vector< Point3D >& coordsFacet,
                                         const Point3D&                oppositeVertex,
                                         const Point3D&                outwardNormal,
                                         const DGBasisInfo&            basis,
                                         int                           degree,
                                         MatrixXr&                     elMat ) const override
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
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                         const std::vector< Point3D >& coordsFacet,
                                         const Point3D&                oppositeVertex,
                                         const Point3D&                outwardNormal,
                                         const DGBasisInfo&            basis,
                                         int                           degree,
                                         MatrixXr&                     elMat ) const override
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
   void integrateVolume3D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_0;
      real_t tmp_6  = p_affine_1_1 + tmp_2;
      real_t tmp_7  = tmp_5 * tmp_6;
      real_t tmp_8  = -p_affine_0_2;
      real_t tmp_9  = p_affine_3_2 + tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_8;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11 * tmp_5;
      real_t tmp_13 = p_affine_2_2 + tmp_8;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = tmp_14 * tmp_6;
      real_t tmp_16 = tmp_1 * tmp_11;
      real_t tmp_17 = tmp_14 * tmp_3;
      real_t tmp_18 =
          1.0 / ( tmp_10 * tmp_12 - tmp_10 * tmp_17 + tmp_13 * tmp_15 - tmp_13 * tmp_16 + tmp_4 * tmp_9 - tmp_7 * tmp_9 );
      real_t tmp_19 = tmp_18 * ( tmp_4 - tmp_7 );
      real_t tmp_20 = tmp_18 * ( tmp_15 - tmp_16 );
      real_t tmp_21 = tmp_18 * ( tmp_12 - tmp_17 );
      real_t tmp_22 = tmp_11 * tmp_19 + tmp_20 * tmp_3 + tmp_21 * tmp_6;
      real_t tmp_23 = tmp_18 * ( -tmp_1 * tmp_13 + tmp_10 * tmp_5 );
      real_t tmp_24 = tmp_18 * ( tmp_1 * tmp_9 - tmp_10 * tmp_14 );
      real_t tmp_25 = tmp_18 * ( tmp_13 * tmp_14 - tmp_5 * tmp_9 );
      real_t tmp_26 = tmp_11 * tmp_23 + tmp_24 * tmp_3 + tmp_25 * tmp_6;
      real_t tmp_27 = tmp_18 * ( -tmp_10 * tmp_3 + tmp_13 * tmp_6 );
      real_t tmp_28 = tmp_18 * ( tmp_10 * tmp_11 - tmp_6 * tmp_9 );
      real_t tmp_29 = tmp_18 * ( -tmp_11 * tmp_13 + tmp_3 * tmp_9 );
      real_t tmp_30 = tmp_11 * tmp_27 + tmp_28 * tmp_3 + tmp_29 * tmp_6;
      real_t tmp_31 = p_affine_0_0 * p_affine_1_1;
      real_t tmp_32 = p_affine_0_0 * p_affine_1_2;
      real_t tmp_33 = p_affine_2_1 * p_affine_3_2;
      real_t tmp_34 = p_affine_0_1 * p_affine_1_0;
      real_t tmp_35 = p_affine_0_1 * p_affine_1_2;
      real_t tmp_36 = p_affine_2_2 * p_affine_3_0;
      real_t tmp_37 = p_affine_0_2 * p_affine_1_0;
      real_t tmp_38 = p_affine_0_2 * p_affine_1_1;
      real_t tmp_39 = p_affine_2_0 * p_affine_3_1;
      real_t tmp_40 = p_affine_2_2 * p_affine_3_1;
      real_t tmp_41 = p_affine_2_0 * p_affine_3_2;
      real_t tmp_42 = p_affine_2_1 * p_affine_3_0;
      real_t tmp_43 = std::abs( p_affine_0_0 * tmp_33 - p_affine_0_0 * tmp_40 + p_affine_0_1 * tmp_36 - p_affine_0_1 * tmp_41 +
                                p_affine_0_2 * tmp_39 - p_affine_0_2 * tmp_42 - p_affine_1_0 * tmp_33 + p_affine_1_0 * tmp_40 -
                                p_affine_1_1 * tmp_36 + p_affine_1_1 * tmp_41 - p_affine_1_2 * tmp_39 + p_affine_1_2 * tmp_42 +
                                p_affine_2_0 * tmp_35 - p_affine_2_0 * tmp_38 - p_affine_2_1 * tmp_32 + p_affine_2_1 * tmp_37 +
                                p_affine_2_2 * tmp_31 - p_affine_2_2 * tmp_34 - p_affine_3_0 * tmp_35 + p_affine_3_0 * tmp_38 +
                                p_affine_3_1 * tmp_32 - p_affine_3_1 * tmp_37 - p_affine_3_2 * tmp_31 + p_affine_3_2 * tmp_34 );
      real_t tmp_44 = tmp_43 * ( tmp_22 * ( -tmp_19 - tmp_20 - tmp_21 ) + tmp_26 * ( -tmp_23 - tmp_24 - tmp_25 ) +
                                 tmp_30 * ( -tmp_27 - tmp_28 - tmp_29 ) );
      real_t tmp_45 = tmp_43 * ( tmp_21 * tmp_22 + tmp_25 * tmp_26 + tmp_29 * tmp_30 );
      real_t tmp_46 = tmp_43 * ( tmp_20 * tmp_22 + tmp_24 * tmp_26 + tmp_28 * tmp_30 );
      real_t tmp_47 = tmp_43 * ( tmp_19 * tmp_22 + tmp_23 * tmp_26 + tmp_27 * tmp_30 );
      real_t a_0_0  = 0.1666666666666668 * tmp_44;
      real_t a_0_1  = 0.1666666666666668 * tmp_45;
      real_t a_0_2  = 0.1666666666666668 * tmp_46;
      real_t a_0_3  = 0.1666666666666668 * tmp_47;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 0, 3 ) = a_0_3;
   }

   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                               const std::vector< Point3D >& coordsFacet,
                               const Point3D&,
                               const Point3D&     outwardNormal,
                               const DGBasisInfo& trialBasis,
                               const DGBasisInfo& testBasis,
                               int                trialDegree,
                               int                testDegree,
                               MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_1_1 + tmp_0;
      real_t tmp_2  = -p_affine_8_2;
      real_t tmp_3  = p_affine_9_2 + tmp_2;
      real_t tmp_4  = p_affine_10_2 + tmp_2;
      real_t tmp_5  = -p_affine_0_2;
      real_t tmp_6  = p_affine_8_2 + tmp_5;
      real_t tmp_7  = 0.031405749086161582 * tmp_3 + 0.93718850182767688 * tmp_4 + tmp_6;
      real_t tmp_8  = -p_affine_0_0;
      real_t tmp_9  = p_affine_2_0 + tmp_8;
      real_t tmp_10 = p_affine_3_1 + tmp_0;
      real_t tmp_11 = p_affine_3_0 + tmp_8;
      real_t tmp_12 = p_affine_2_1 + tmp_0;
      real_t tmp_13 = p_affine_1_0 + tmp_8;
      real_t tmp_14 = p_affine_3_2 + tmp_5;
      real_t tmp_15 = tmp_12 * tmp_14;
      real_t tmp_16 = p_affine_1_2 + tmp_5;
      real_t tmp_17 = tmp_10 * tmp_16;
      real_t tmp_18 = p_affine_2_2 + tmp_5;
      real_t tmp_19 = tmp_1 * tmp_18;
      real_t tmp_20 = tmp_10 * tmp_18;
      real_t tmp_21 = tmp_1 * tmp_14;
      real_t tmp_22 = tmp_12 * tmp_16;
      real_t tmp_23 =
          1.0 / ( tmp_11 * tmp_19 - tmp_11 * tmp_22 + tmp_13 * tmp_15 - tmp_13 * tmp_20 + tmp_17 * tmp_9 - tmp_21 * tmp_9 );
      real_t tmp_24 = tmp_23 * ( tmp_10 * tmp_9 - tmp_11 * tmp_12 );
      real_t tmp_25 = tmp_24 * tmp_7;
      real_t tmp_26 = -p_affine_8_1;
      real_t tmp_27 = p_affine_9_1 + tmp_26;
      real_t tmp_28 = p_affine_10_1 + tmp_26;
      real_t tmp_29 = p_affine_8_1 + tmp_0;
      real_t tmp_30 = 0.031405749086161582 * tmp_27 + 0.93718850182767688 * tmp_28 + tmp_29;
      real_t tmp_31 = tmp_23 * ( tmp_11 * tmp_18 - tmp_14 * tmp_9 );
      real_t tmp_32 = tmp_30 * tmp_31;
      real_t tmp_33 = -p_affine_8_0;
      real_t tmp_34 = p_affine_9_0 + tmp_33;
      real_t tmp_35 = p_affine_10_0 + tmp_33;
      real_t tmp_36 = p_affine_8_0 + tmp_8;
      real_t tmp_37 = 0.031405749086161582 * tmp_34 + 0.93718850182767688 * tmp_35 + tmp_36;
      real_t tmp_38 = tmp_23 * ( tmp_15 - tmp_20 );
      real_t tmp_39 = tmp_37 * tmp_38;
      real_t tmp_40 = tmp_25 + tmp_32 + tmp_39;
      real_t tmp_41 = tmp_23 * ( tmp_1 * tmp_11 - tmp_10 * tmp_13 );
      real_t tmp_42 = tmp_41 * tmp_7;
      real_t tmp_43 = tmp_23 * ( -tmp_11 * tmp_16 + tmp_13 * tmp_14 );
      real_t tmp_44 = tmp_30 * tmp_43;
      real_t tmp_45 = tmp_23 * ( tmp_17 - tmp_21 );
      real_t tmp_46 = tmp_37 * tmp_45;
      real_t tmp_47 = tmp_42 + tmp_44 + tmp_46;
      real_t tmp_48 = tmp_23 * ( -tmp_1 * tmp_9 + tmp_12 * tmp_13 );
      real_t tmp_49 = tmp_48 * tmp_7;
      real_t tmp_50 = tmp_23 * ( -tmp_13 * tmp_18 + tmp_16 * tmp_9 );
      real_t tmp_51 = tmp_30 * tmp_50;
      real_t tmp_52 = tmp_23 * ( tmp_19 - tmp_22 );
      real_t tmp_53 = tmp_37 * tmp_52;
      real_t tmp_54 = tmp_49 + tmp_51 + tmp_53;
      real_t tmp_55 = tmp_1 * ( tmp_40 - 1.0 / 4.0 ) + tmp_10 * ( tmp_54 - 1.0 / 4.0 ) + tmp_12 * ( tmp_47 - 1.0 / 4.0 );
      real_t tmp_56 = 0.5 * p_affine_13_0 * ( -tmp_38 - tmp_45 - tmp_52 ) + 0.5 * p_affine_13_1 * ( -tmp_31 - tmp_43 - tmp_50 ) +
                      0.5 * p_affine_13_2 * ( -tmp_24 - tmp_41 - tmp_48 );
      real_t tmp_57 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_58 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_59 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_60 = ( std::abs( tmp_28 * tmp_58 - tmp_35 * tmp_57 ) * std::abs( tmp_28 * tmp_58 - tmp_35 * tmp_57 ) ) +
                      ( std::abs( tmp_28 * tmp_59 - tmp_4 * tmp_57 ) * std::abs( tmp_28 * tmp_59 - tmp_4 * tmp_57 ) ) +
                      ( std::abs( tmp_35 * tmp_59 - tmp_4 * tmp_58 ) * std::abs( tmp_35 * tmp_59 - tmp_4 * tmp_58 ) );
      real_t tmp_61  = 1.0 * std::pow( tmp_60, -0.25 );
      real_t tmp_62  = tmp_55 * tmp_61;
      real_t tmp_63  = 1.0 * std::pow( tmp_60, 1.0 / 2.0 );
      real_t tmp_64  = 0.0068572537431980923 * tmp_63;
      real_t tmp_65  = 0.19601935860219369 * tmp_3 + 0.60796128279561268 * tmp_4 + tmp_6;
      real_t tmp_66  = tmp_24 * tmp_65;
      real_t tmp_67  = 0.19601935860219369 * tmp_27 + 0.60796128279561268 * tmp_28 + tmp_29;
      real_t tmp_68  = tmp_31 * tmp_67;
      real_t tmp_69  = 0.19601935860219369 * tmp_34 + 0.60796128279561268 * tmp_35 + tmp_36;
      real_t tmp_70  = tmp_38 * tmp_69;
      real_t tmp_71  = tmp_66 + tmp_68 + tmp_70;
      real_t tmp_72  = tmp_41 * tmp_65;
      real_t tmp_73  = tmp_43 * tmp_67;
      real_t tmp_74  = tmp_45 * tmp_69;
      real_t tmp_75  = tmp_72 + tmp_73 + tmp_74;
      real_t tmp_76  = tmp_48 * tmp_65;
      real_t tmp_77  = tmp_50 * tmp_67;
      real_t tmp_78  = tmp_52 * tmp_69;
      real_t tmp_79  = tmp_76 + tmp_77 + tmp_78;
      real_t tmp_80  = tmp_1 * ( tmp_71 - 1.0 / 4.0 ) + tmp_10 * ( tmp_79 - 1.0 / 4.0 ) + tmp_12 * ( tmp_75 - 1.0 / 4.0 );
      real_t tmp_81  = tmp_61 * tmp_80;
      real_t tmp_82  = 0.037198804536718075 * tmp_63;
      real_t tmp_83  = 0.37605877282253791 * tmp_3 + 0.039308471900058539 * tmp_4 + tmp_6;
      real_t tmp_84  = tmp_24 * tmp_83;
      real_t tmp_85  = 0.37605877282253791 * tmp_27 + 0.039308471900058539 * tmp_28 + tmp_29;
      real_t tmp_86  = tmp_31 * tmp_85;
      real_t tmp_87  = 0.37605877282253791 * tmp_34 + 0.039308471900058539 * tmp_35 + tmp_36;
      real_t tmp_88  = tmp_38 * tmp_87;
      real_t tmp_89  = tmp_84 + tmp_86 + tmp_88;
      real_t tmp_90  = tmp_41 * tmp_83;
      real_t tmp_91  = tmp_43 * tmp_85;
      real_t tmp_92  = tmp_45 * tmp_87;
      real_t tmp_93  = tmp_90 + tmp_91 + tmp_92;
      real_t tmp_94  = tmp_48 * tmp_83;
      real_t tmp_95  = tmp_50 * tmp_85;
      real_t tmp_96  = tmp_52 * tmp_87;
      real_t tmp_97  = tmp_94 + tmp_95 + tmp_96;
      real_t tmp_98  = tmp_1 * ( tmp_89 - 1.0 / 4.0 ) + tmp_10 * ( tmp_97 - 1.0 / 4.0 ) + tmp_12 * ( tmp_93 - 1.0 / 4.0 );
      real_t tmp_99  = tmp_61 * tmp_98;
      real_t tmp_100 = 0.020848748529055869 * tmp_63;
      real_t tmp_101 = 0.78764240869137092 * tmp_3 + 0.1711304259088916 * tmp_4 + tmp_6;
      real_t tmp_102 = tmp_101 * tmp_24;
      real_t tmp_103 = 0.78764240869137092 * tmp_27 + 0.1711304259088916 * tmp_28 + tmp_29;
      real_t tmp_104 = tmp_103 * tmp_31;
      real_t tmp_105 = 0.78764240869137092 * tmp_34 + 0.1711304259088916 * tmp_35 + tmp_36;
      real_t tmp_106 = tmp_105 * tmp_38;
      real_t tmp_107 = tmp_102 + tmp_104 + tmp_106;
      real_t tmp_108 = tmp_101 * tmp_41;
      real_t tmp_109 = tmp_103 * tmp_43;
      real_t tmp_110 = tmp_105 * tmp_45;
      real_t tmp_111 = tmp_108 + tmp_109 + tmp_110;
      real_t tmp_112 = tmp_101 * tmp_48;
      real_t tmp_113 = tmp_103 * tmp_50;
      real_t tmp_114 = tmp_105 * tmp_52;
      real_t tmp_115 = tmp_112 + tmp_113 + tmp_114;
      real_t tmp_116 = tmp_1 * ( tmp_107 - 1.0 / 4.0 ) + tmp_10 * ( tmp_115 - 1.0 / 4.0 ) + tmp_12 * ( tmp_111 - 1.0 / 4.0 );
      real_t tmp_117 = tmp_116 * tmp_61;
      real_t tmp_118 = 0.019202922745021479 * tmp_63;
      real_t tmp_119 = 0.58463275527740355 * tmp_3 + 0.37605877282253791 * tmp_4 + tmp_6;
      real_t tmp_120 = tmp_119 * tmp_24;
      real_t tmp_121 = 0.58463275527740355 * tmp_27 + 0.37605877282253791 * tmp_28 + tmp_29;
      real_t tmp_122 = tmp_121 * tmp_31;
      real_t tmp_123 = 0.58463275527740355 * tmp_34 + 0.37605877282253791 * tmp_35 + tmp_36;
      real_t tmp_124 = tmp_123 * tmp_38;
      real_t tmp_125 = tmp_120 + tmp_122 + tmp_124;
      real_t tmp_126 = tmp_119 * tmp_41;
      real_t tmp_127 = tmp_121 * tmp_43;
      real_t tmp_128 = tmp_123 * tmp_45;
      real_t tmp_129 = tmp_126 + tmp_127 + tmp_128;
      real_t tmp_130 = tmp_119 * tmp_48;
      real_t tmp_131 = tmp_121 * tmp_50;
      real_t tmp_132 = tmp_123 * tmp_52;
      real_t tmp_133 = tmp_130 + tmp_131 + tmp_132;
      real_t tmp_134 = tmp_1 * ( tmp_125 - 1.0 / 4.0 ) + tmp_10 * ( tmp_133 - 1.0 / 4.0 ) + tmp_12 * ( tmp_129 - 1.0 / 4.0 );
      real_t tmp_135 = tmp_134 * tmp_61;
      real_t tmp_136 = 0.020848748529055869 * tmp_63;
      real_t tmp_137 = 0.041227165399737475 * tmp_3 + 0.78764240869137092 * tmp_4 + tmp_6;
      real_t tmp_138 = tmp_137 * tmp_24;
      real_t tmp_139 = 0.041227165399737475 * tmp_27 + 0.78764240869137092 * tmp_28 + tmp_29;
      real_t tmp_140 = tmp_139 * tmp_31;
      real_t tmp_141 = 0.041227165399737475 * tmp_34 + 0.78764240869137092 * tmp_35 + tmp_36;
      real_t tmp_142 = tmp_141 * tmp_38;
      real_t tmp_143 = tmp_138 + tmp_140 + tmp_142;
      real_t tmp_144 = tmp_137 * tmp_41;
      real_t tmp_145 = tmp_139 * tmp_43;
      real_t tmp_146 = tmp_141 * tmp_45;
      real_t tmp_147 = tmp_144 + tmp_145 + tmp_146;
      real_t tmp_148 = tmp_137 * tmp_48;
      real_t tmp_149 = tmp_139 * tmp_50;
      real_t tmp_150 = tmp_141 * tmp_52;
      real_t tmp_151 = tmp_148 + tmp_149 + tmp_150;
      real_t tmp_152 = tmp_1 * ( tmp_143 - 1.0 / 4.0 ) + tmp_10 * ( tmp_151 - 1.0 / 4.0 ) + tmp_12 * ( tmp_147 - 1.0 / 4.0 );
      real_t tmp_153 = tmp_152 * tmp_61;
      real_t tmp_154 = 0.019202922745021479 * tmp_63;
      real_t tmp_155 = 0.039308471900058539 * tmp_3 + 0.58463275527740355 * tmp_4 + tmp_6;
      real_t tmp_156 = tmp_155 * tmp_24;
      real_t tmp_157 = 0.039308471900058539 * tmp_27 + 0.58463275527740355 * tmp_28 + tmp_29;
      real_t tmp_158 = tmp_157 * tmp_31;
      real_t tmp_159 = 0.039308471900058539 * tmp_34 + 0.58463275527740355 * tmp_35 + tmp_36;
      real_t tmp_160 = tmp_159 * tmp_38;
      real_t tmp_161 = tmp_156 + tmp_158 + tmp_160;
      real_t tmp_162 = tmp_155 * tmp_41;
      real_t tmp_163 = tmp_157 * tmp_43;
      real_t tmp_164 = tmp_159 * tmp_45;
      real_t tmp_165 = tmp_162 + tmp_163 + tmp_164;
      real_t tmp_166 = tmp_155 * tmp_48;
      real_t tmp_167 = tmp_157 * tmp_50;
      real_t tmp_168 = tmp_159 * tmp_52;
      real_t tmp_169 = tmp_166 + tmp_167 + tmp_168;
      real_t tmp_170 = tmp_1 * ( tmp_161 - 1.0 / 4.0 ) + tmp_10 * ( tmp_169 - 1.0 / 4.0 ) + tmp_12 * ( tmp_165 - 1.0 / 4.0 );
      real_t tmp_171 = tmp_170 * tmp_61;
      real_t tmp_172 = 0.020848748529055869 * tmp_63;
      real_t tmp_173 = 0.78764240869137092 * tmp_3 + 0.041227165399737475 * tmp_4 + tmp_6;
      real_t tmp_174 = tmp_173 * tmp_24;
      real_t tmp_175 = 0.78764240869137092 * tmp_27 + 0.041227165399737475 * tmp_28 + tmp_29;
      real_t tmp_176 = tmp_175 * tmp_31;
      real_t tmp_177 = 0.78764240869137092 * tmp_34 + 0.041227165399737475 * tmp_35 + tmp_36;
      real_t tmp_178 = tmp_177 * tmp_38;
      real_t tmp_179 = tmp_174 + tmp_176 + tmp_178;
      real_t tmp_180 = tmp_173 * tmp_41;
      real_t tmp_181 = tmp_175 * tmp_43;
      real_t tmp_182 = tmp_177 * tmp_45;
      real_t tmp_183 = tmp_180 + tmp_181 + tmp_182;
      real_t tmp_184 = tmp_173 * tmp_48;
      real_t tmp_185 = tmp_175 * tmp_50;
      real_t tmp_186 = tmp_177 * tmp_52;
      real_t tmp_187 = tmp_184 + tmp_185 + tmp_186;
      real_t tmp_188 = tmp_1 * ( tmp_179 - 1.0 / 4.0 ) + tmp_10 * ( tmp_187 - 1.0 / 4.0 ) + tmp_12 * ( tmp_183 - 1.0 / 4.0 );
      real_t tmp_189 = tmp_188 * tmp_61;
      real_t tmp_190 = 0.019202922745021479 * tmp_63;
      real_t tmp_191 = 0.58463275527740355 * tmp_3 + 0.039308471900058539 * tmp_4 + tmp_6;
      real_t tmp_192 = tmp_191 * tmp_24;
      real_t tmp_193 = 0.58463275527740355 * tmp_27 + 0.039308471900058539 * tmp_28 + tmp_29;
      real_t tmp_194 = tmp_193 * tmp_31;
      real_t tmp_195 = 0.58463275527740355 * tmp_34 + 0.039308471900058539 * tmp_35 + tmp_36;
      real_t tmp_196 = tmp_195 * tmp_38;
      real_t tmp_197 = tmp_192 + tmp_194 + tmp_196;
      real_t tmp_198 = tmp_191 * tmp_41;
      real_t tmp_199 = tmp_193 * tmp_43;
      real_t tmp_200 = tmp_195 * tmp_45;
      real_t tmp_201 = tmp_198 + tmp_199 + tmp_200;
      real_t tmp_202 = tmp_191 * tmp_48;
      real_t tmp_203 = tmp_193 * tmp_50;
      real_t tmp_204 = tmp_195 * tmp_52;
      real_t tmp_205 = tmp_202 + tmp_203 + tmp_204;
      real_t tmp_206 = tmp_1 * ( tmp_197 - 1.0 / 4.0 ) + tmp_10 * ( tmp_205 - 1.0 / 4.0 ) + tmp_12 * ( tmp_201 - 1.0 / 4.0 );
      real_t tmp_207 = tmp_206 * tmp_61;
      real_t tmp_208 = 0.020848748529055869 * tmp_63;
      real_t tmp_209 = 0.1711304259088916 * tmp_3 + 0.78764240869137092 * tmp_4 + tmp_6;
      real_t tmp_210 = tmp_209 * tmp_24;
      real_t tmp_211 = 0.1711304259088916 * tmp_27 + 0.78764240869137092 * tmp_28 + tmp_29;
      real_t tmp_212 = tmp_211 * tmp_31;
      real_t tmp_213 = 0.1711304259088916 * tmp_34 + 0.78764240869137092 * tmp_35 + tmp_36;
      real_t tmp_214 = tmp_213 * tmp_38;
      real_t tmp_215 = tmp_210 + tmp_212 + tmp_214;
      real_t tmp_216 = tmp_209 * tmp_41;
      real_t tmp_217 = tmp_211 * tmp_43;
      real_t tmp_218 = tmp_213 * tmp_45;
      real_t tmp_219 = tmp_216 + tmp_217 + tmp_218;
      real_t tmp_220 = tmp_209 * tmp_48;
      real_t tmp_221 = tmp_211 * tmp_50;
      real_t tmp_222 = tmp_213 * tmp_52;
      real_t tmp_223 = tmp_220 + tmp_221 + tmp_222;
      real_t tmp_224 = tmp_1 * ( tmp_215 - 1.0 / 4.0 ) + tmp_10 * ( tmp_223 - 1.0 / 4.0 ) + tmp_12 * ( tmp_219 - 1.0 / 4.0 );
      real_t tmp_225 = tmp_224 * tmp_61;
      real_t tmp_226 = 0.019202922745021479 * tmp_63;
      real_t tmp_227 = 0.37605877282253791 * tmp_3 + 0.58463275527740355 * tmp_4 + tmp_6;
      real_t tmp_228 = tmp_227 * tmp_24;
      real_t tmp_229 = 0.37605877282253791 * tmp_27 + 0.58463275527740355 * tmp_28 + tmp_29;
      real_t tmp_230 = tmp_229 * tmp_31;
      real_t tmp_231 = 0.37605877282253791 * tmp_34 + 0.58463275527740355 * tmp_35 + tmp_36;
      real_t tmp_232 = tmp_231 * tmp_38;
      real_t tmp_233 = tmp_228 + tmp_230 + tmp_232;
      real_t tmp_234 = tmp_227 * tmp_41;
      real_t tmp_235 = tmp_229 * tmp_43;
      real_t tmp_236 = tmp_231 * tmp_45;
      real_t tmp_237 = tmp_234 + tmp_235 + tmp_236;
      real_t tmp_238 = tmp_227 * tmp_48;
      real_t tmp_239 = tmp_229 * tmp_50;
      real_t tmp_240 = tmp_231 * tmp_52;
      real_t tmp_241 = tmp_238 + tmp_239 + tmp_240;
      real_t tmp_242 = tmp_1 * ( tmp_233 - 1.0 / 4.0 ) + tmp_10 * ( tmp_241 - 1.0 / 4.0 ) + tmp_12 * ( tmp_237 - 1.0 / 4.0 );
      real_t tmp_243 = tmp_242 * tmp_61;
      real_t tmp_244 = 0.020848748529055869 * tmp_63;
      real_t tmp_245 = 0.041227165399737475 * tmp_3 + 0.1711304259088916 * tmp_4 + tmp_6;
      real_t tmp_246 = tmp_24 * tmp_245;
      real_t tmp_247 = 0.041227165399737475 * tmp_27 + 0.1711304259088916 * tmp_28 + tmp_29;
      real_t tmp_248 = tmp_247 * tmp_31;
      real_t tmp_249 = 0.041227165399737475 * tmp_34 + 0.1711304259088916 * tmp_35 + tmp_36;
      real_t tmp_250 = tmp_249 * tmp_38;
      real_t tmp_251 = tmp_246 + tmp_248 + tmp_250;
      real_t tmp_252 = tmp_245 * tmp_41;
      real_t tmp_253 = tmp_247 * tmp_43;
      real_t tmp_254 = tmp_249 * tmp_45;
      real_t tmp_255 = tmp_252 + tmp_253 + tmp_254;
      real_t tmp_256 = tmp_245 * tmp_48;
      real_t tmp_257 = tmp_247 * tmp_50;
      real_t tmp_258 = tmp_249 * tmp_52;
      real_t tmp_259 = tmp_256 + tmp_257 + tmp_258;
      real_t tmp_260 = tmp_1 * ( tmp_251 - 1.0 / 4.0 ) + tmp_10 * ( tmp_259 - 1.0 / 4.0 ) + tmp_12 * ( tmp_255 - 1.0 / 4.0 );
      real_t tmp_261 = tmp_260 * tmp_61;
      real_t tmp_262 = 0.019202922745021479 * tmp_63;
      real_t tmp_263 = 0.40446199974765351 * tmp_3 + 0.19107600050469298 * tmp_4 + tmp_6;
      real_t tmp_264 = tmp_24 * tmp_263;
      real_t tmp_265 = 0.40446199974765351 * tmp_27 + 0.19107600050469298 * tmp_28 + tmp_29;
      real_t tmp_266 = tmp_265 * tmp_31;
      real_t tmp_267 = 0.40446199974765351 * tmp_34 + 0.19107600050469298 * tmp_35 + tmp_36;
      real_t tmp_268 = tmp_267 * tmp_38;
      real_t tmp_269 = tmp_264 + tmp_266 + tmp_268;
      real_t tmp_270 = tmp_263 * tmp_41;
      real_t tmp_271 = tmp_265 * tmp_43;
      real_t tmp_272 = tmp_267 * tmp_45;
      real_t tmp_273 = tmp_270 + tmp_271 + tmp_272;
      real_t tmp_274 = tmp_263 * tmp_48;
      real_t tmp_275 = tmp_265 * tmp_50;
      real_t tmp_276 = tmp_267 * tmp_52;
      real_t tmp_277 = tmp_274 + tmp_275 + tmp_276;
      real_t tmp_278 = tmp_1 * ( tmp_269 - 1.0 / 4.0 ) + tmp_10 * ( tmp_277 - 1.0 / 4.0 ) + tmp_12 * ( tmp_273 - 1.0 / 4.0 );
      real_t tmp_279 = tmp_278 * tmp_61;
      real_t tmp_280 = 0.042507265838595799 * tmp_63;
      real_t tmp_281 = 0.039308471900058539 * tmp_3 + 0.37605877282253791 * tmp_4 + tmp_6;
      real_t tmp_282 = tmp_24 * tmp_281;
      real_t tmp_283 = 0.039308471900058539 * tmp_27 + 0.37605877282253791 * tmp_28 + tmp_29;
      real_t tmp_284 = tmp_283 * tmp_31;
      real_t tmp_285 = 0.039308471900058539 * tmp_34 + 0.37605877282253791 * tmp_35 + tmp_36;
      real_t tmp_286 = tmp_285 * tmp_38;
      real_t tmp_287 = tmp_282 + tmp_284 + tmp_286;
      real_t tmp_288 = tmp_281 * tmp_41;
      real_t tmp_289 = tmp_283 * tmp_43;
      real_t tmp_290 = tmp_285 * tmp_45;
      real_t tmp_291 = tmp_288 + tmp_289 + tmp_290;
      real_t tmp_292 = tmp_281 * tmp_48;
      real_t tmp_293 = tmp_283 * tmp_50;
      real_t tmp_294 = tmp_285 * tmp_52;
      real_t tmp_295 = tmp_292 + tmp_293 + tmp_294;
      real_t tmp_296 = tmp_1 * ( tmp_287 - 1.0 / 4.0 ) + tmp_10 * ( tmp_295 - 1.0 / 4.0 ) + tmp_12 * ( tmp_291 - 1.0 / 4.0 );
      real_t tmp_297 = tmp_296 * tmp_61;
      real_t tmp_298 = 0.020848748529055869 * tmp_63;
      real_t tmp_299 = 0.93718850182767688 * tmp_3 + 0.031405749086161582 * tmp_4 + tmp_6;
      real_t tmp_300 = tmp_24 * tmp_299;
      real_t tmp_301 = 0.93718850182767688 * tmp_27 + 0.031405749086161582 * tmp_28 + tmp_29;
      real_t tmp_302 = tmp_301 * tmp_31;
      real_t tmp_303 = 0.93718850182767688 * tmp_34 + 0.031405749086161582 * tmp_35 + tmp_36;
      real_t tmp_304 = tmp_303 * tmp_38;
      real_t tmp_305 = tmp_300 + tmp_302 + tmp_304;
      real_t tmp_306 = tmp_299 * tmp_41;
      real_t tmp_307 = tmp_301 * tmp_43;
      real_t tmp_308 = tmp_303 * tmp_45;
      real_t tmp_309 = tmp_306 + tmp_307 + tmp_308;
      real_t tmp_310 = tmp_299 * tmp_48;
      real_t tmp_311 = tmp_301 * tmp_50;
      real_t tmp_312 = tmp_303 * tmp_52;
      real_t tmp_313 = tmp_310 + tmp_311 + tmp_312;
      real_t tmp_314 = tmp_1 * ( tmp_305 - 1.0 / 4.0 ) + tmp_10 * ( tmp_313 - 1.0 / 4.0 ) + tmp_12 * ( tmp_309 - 1.0 / 4.0 );
      real_t tmp_315 = tmp_314 * tmp_61;
      real_t tmp_316 = 0.0068572537431980923 * tmp_63;
      real_t tmp_317 = 0.60796128279561268 * tmp_3 + 0.19601935860219369 * tmp_4 + tmp_6;
      real_t tmp_318 = tmp_24 * tmp_317;
      real_t tmp_319 = 0.60796128279561268 * tmp_27 + 0.19601935860219369 * tmp_28 + tmp_29;
      real_t tmp_320 = tmp_31 * tmp_319;
      real_t tmp_321 = 0.60796128279561268 * tmp_34 + 0.19601935860219369 * tmp_35 + tmp_36;
      real_t tmp_322 = tmp_321 * tmp_38;
      real_t tmp_323 = tmp_318 + tmp_320 + tmp_322;
      real_t tmp_324 = tmp_317 * tmp_41;
      real_t tmp_325 = tmp_319 * tmp_43;
      real_t tmp_326 = tmp_321 * tmp_45;
      real_t tmp_327 = tmp_324 + tmp_325 + tmp_326;
      real_t tmp_328 = tmp_317 * tmp_48;
      real_t tmp_329 = tmp_319 * tmp_50;
      real_t tmp_330 = tmp_321 * tmp_52;
      real_t tmp_331 = tmp_328 + tmp_329 + tmp_330;
      real_t tmp_332 = tmp_1 * ( tmp_323 - 1.0 / 4.0 ) + tmp_10 * ( tmp_331 - 1.0 / 4.0 ) + tmp_12 * ( tmp_327 - 1.0 / 4.0 );
      real_t tmp_333 = tmp_332 * tmp_61;
      real_t tmp_334 = 0.037198804536718075 * tmp_63;
      real_t tmp_335 = 0.19107600050469298 * tmp_3 + 0.40446199974765351 * tmp_4 + tmp_6;
      real_t tmp_336 = tmp_24 * tmp_335;
      real_t tmp_337 = 0.19107600050469298 * tmp_27 + 0.40446199974765351 * tmp_28 + tmp_29;
      real_t tmp_338 = tmp_31 * tmp_337;
      real_t tmp_339 = 0.19107600050469298 * tmp_34 + 0.40446199974765351 * tmp_35 + tmp_36;
      real_t tmp_340 = tmp_339 * tmp_38;
      real_t tmp_341 = tmp_336 + tmp_338 + tmp_340;
      real_t tmp_342 = tmp_335 * tmp_41;
      real_t tmp_343 = tmp_337 * tmp_43;
      real_t tmp_344 = tmp_339 * tmp_45;
      real_t tmp_345 = tmp_342 + tmp_343 + tmp_344;
      real_t tmp_346 = tmp_335 * tmp_48;
      real_t tmp_347 = tmp_337 * tmp_50;
      real_t tmp_348 = tmp_339 * tmp_52;
      real_t tmp_349 = tmp_346 + tmp_347 + tmp_348;
      real_t tmp_350 = tmp_1 * ( tmp_341 - 1.0 / 4.0 ) + tmp_10 * ( tmp_349 - 1.0 / 4.0 ) + tmp_12 * ( tmp_345 - 1.0 / 4.0 );
      real_t tmp_351 = tmp_350 * tmp_61;
      real_t tmp_352 = 0.042507265838595799 * tmp_63;
      real_t tmp_353 = 0.031405749086161582 * tmp_3 + 0.031405749086161582 * tmp_4 + tmp_6;
      real_t tmp_354 = tmp_24 * tmp_353;
      real_t tmp_355 = 0.031405749086161582 * tmp_27 + 0.031405749086161582 * tmp_28 + tmp_29;
      real_t tmp_356 = tmp_31 * tmp_355;
      real_t tmp_357 = 0.031405749086161582 * tmp_34 + 0.031405749086161582 * tmp_35 + tmp_36;
      real_t tmp_358 = tmp_357 * tmp_38;
      real_t tmp_359 = tmp_354 + tmp_356 + tmp_358;
      real_t tmp_360 = tmp_353 * tmp_41;
      real_t tmp_361 = tmp_355 * tmp_43;
      real_t tmp_362 = tmp_357 * tmp_45;
      real_t tmp_363 = tmp_360 + tmp_361 + tmp_362;
      real_t tmp_364 = tmp_353 * tmp_48;
      real_t tmp_365 = tmp_355 * tmp_50;
      real_t tmp_366 = tmp_357 * tmp_52;
      real_t tmp_367 = tmp_364 + tmp_365 + tmp_366;
      real_t tmp_368 = tmp_1 * ( tmp_359 - 1.0 / 4.0 ) + tmp_10 * ( tmp_367 - 1.0 / 4.0 ) + tmp_12 * ( tmp_363 - 1.0 / 4.0 );
      real_t tmp_369 = tmp_368 * tmp_61;
      real_t tmp_370 = 0.0068572537431980923 * tmp_63;
      real_t tmp_371 = 0.19601935860219369 * tmp_3 + 0.19601935860219369 * tmp_4 + tmp_6;
      real_t tmp_372 = tmp_24 * tmp_371;
      real_t tmp_373 = 0.19601935860219369 * tmp_27 + 0.19601935860219369 * tmp_28 + tmp_29;
      real_t tmp_374 = tmp_31 * tmp_373;
      real_t tmp_375 = 0.19601935860219369 * tmp_34 + 0.19601935860219369 * tmp_35 + tmp_36;
      real_t tmp_376 = tmp_375 * tmp_38;
      real_t tmp_377 = tmp_372 + tmp_374 + tmp_376;
      real_t tmp_378 = tmp_371 * tmp_41;
      real_t tmp_379 = tmp_373 * tmp_43;
      real_t tmp_380 = tmp_375 * tmp_45;
      real_t tmp_381 = tmp_378 + tmp_379 + tmp_380;
      real_t tmp_382 = tmp_371 * tmp_48;
      real_t tmp_383 = tmp_373 * tmp_50;
      real_t tmp_384 = tmp_375 * tmp_52;
      real_t tmp_385 = tmp_382 + tmp_383 + tmp_384;
      real_t tmp_386 = tmp_1 * ( tmp_377 - 1.0 / 4.0 ) + tmp_10 * ( tmp_385 - 1.0 / 4.0 ) + tmp_12 * ( tmp_381 - 1.0 / 4.0 );
      real_t tmp_387 = tmp_386 * tmp_61;
      real_t tmp_388 = 0.037198804536718075 * tmp_63;
      real_t tmp_389 = 0.40446199974765351 * tmp_3 + 0.40446199974765351 * tmp_4 + tmp_6;
      real_t tmp_390 = tmp_24 * tmp_389;
      real_t tmp_391 = 0.40446199974765351 * tmp_27 + 0.40446199974765351 * tmp_28 + tmp_29;
      real_t tmp_392 = tmp_31 * tmp_391;
      real_t tmp_393 = 0.40446199974765351 * tmp_34 + 0.40446199974765351 * tmp_35 + tmp_36;
      real_t tmp_394 = tmp_38 * tmp_393;
      real_t tmp_395 = tmp_390 + tmp_392 + tmp_394;
      real_t tmp_396 = tmp_389 * tmp_41;
      real_t tmp_397 = tmp_391 * tmp_43;
      real_t tmp_398 = tmp_393 * tmp_45;
      real_t tmp_399 = tmp_396 + tmp_397 + tmp_398;
      real_t tmp_400 = tmp_389 * tmp_48;
      real_t tmp_401 = tmp_391 * tmp_50;
      real_t tmp_402 = tmp_393 * tmp_52;
      real_t tmp_403 = tmp_400 + tmp_401 + tmp_402;
      real_t tmp_404 = tmp_1 * ( tmp_395 - 1.0 / 4.0 ) + tmp_10 * ( tmp_403 - 1.0 / 4.0 ) + tmp_12 * ( tmp_399 - 1.0 / 4.0 );
      real_t tmp_405 = tmp_404 * tmp_61;
      real_t tmp_406 = 0.042507265838595799 * tmp_63;
      real_t tmp_407 = 0.1711304259088916 * tmp_3 + 0.041227165399737475 * tmp_4 + tmp_6;
      real_t tmp_408 = tmp_24 * tmp_407;
      real_t tmp_409 = 0.1711304259088916 * tmp_27 + 0.041227165399737475 * tmp_28 + tmp_29;
      real_t tmp_410 = tmp_31 * tmp_409;
      real_t tmp_411 = 0.1711304259088916 * tmp_34 + 0.041227165399737475 * tmp_35 + tmp_36;
      real_t tmp_412 = tmp_38 * tmp_411;
      real_t tmp_413 = tmp_408 + tmp_410 + tmp_412;
      real_t tmp_414 = tmp_407 * tmp_41;
      real_t tmp_415 = tmp_409 * tmp_43;
      real_t tmp_416 = tmp_411 * tmp_45;
      real_t tmp_417 = tmp_414 + tmp_415 + tmp_416;
      real_t tmp_418 = tmp_407 * tmp_48;
      real_t tmp_419 = tmp_409 * tmp_50;
      real_t tmp_420 = tmp_411 * tmp_52;
      real_t tmp_421 = tmp_418 + tmp_419 + tmp_420;
      real_t tmp_422 = tmp_1 * ( tmp_413 - 1.0 / 4.0 ) + tmp_10 * ( tmp_421 - 1.0 / 4.0 ) + tmp_12 * ( tmp_417 - 1.0 / 4.0 );
      real_t tmp_423 = tmp_422 * tmp_61;
      real_t tmp_424 = 0.019202922745021479 * tmp_63;
      real_t tmp_425 = 0.5 * p_affine_13_0 * tmp_38 + 0.5 * p_affine_13_1 * tmp_31 + 0.5 * p_affine_13_2 * tmp_24;
      real_t tmp_426 = 0.5 * p_affine_13_0 * tmp_45 + 0.5 * p_affine_13_1 * tmp_43 + 0.5 * p_affine_13_2 * tmp_41;
      real_t tmp_427 = 0.5 * p_affine_13_0 * tmp_52 + 0.5 * p_affine_13_1 * tmp_50 + 0.5 * p_affine_13_2 * tmp_48;
      real_t a_0_0 =
          tmp_100 * ( -tmp_56 * tmp_98 +
                      tmp_99 * ( -tmp_84 - tmp_86 - tmp_88 - tmp_90 - tmp_91 - tmp_92 - tmp_94 - tmp_95 - tmp_96 + 1 ) ) +
          tmp_118 * ( -tmp_116 * tmp_56 + tmp_117 * ( -tmp_102 - tmp_104 - tmp_106 - tmp_108 - tmp_109 - tmp_110 - tmp_112 -
                                                      tmp_113 - tmp_114 + 1 ) ) +
          tmp_136 * ( -tmp_134 * tmp_56 + tmp_135 * ( -tmp_120 - tmp_122 - tmp_124 - tmp_126 - tmp_127 - tmp_128 - tmp_130 -
                                                      tmp_131 - tmp_132 + 1 ) ) +
          tmp_154 * ( -tmp_152 * tmp_56 + tmp_153 * ( -tmp_138 - tmp_140 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_148 -
                                                      tmp_149 - tmp_150 + 1 ) ) +
          tmp_172 * ( -tmp_170 * tmp_56 + tmp_171 * ( -tmp_156 - tmp_158 - tmp_160 - tmp_162 - tmp_163 - tmp_164 - tmp_166 -
                                                      tmp_167 - tmp_168 + 1 ) ) +
          tmp_190 * ( -tmp_188 * tmp_56 + tmp_189 * ( -tmp_174 - tmp_176 - tmp_178 - tmp_180 - tmp_181 - tmp_182 - tmp_184 -
                                                      tmp_185 - tmp_186 + 1 ) ) +
          tmp_208 * ( -tmp_206 * tmp_56 + tmp_207 * ( -tmp_192 - tmp_194 - tmp_196 - tmp_198 - tmp_199 - tmp_200 - tmp_202 -
                                                      tmp_203 - tmp_204 + 1 ) ) +
          tmp_226 * ( -tmp_224 * tmp_56 + tmp_225 * ( -tmp_210 - tmp_212 - tmp_214 - tmp_216 - tmp_217 - tmp_218 - tmp_220 -
                                                      tmp_221 - tmp_222 + 1 ) ) +
          tmp_244 * ( -tmp_242 * tmp_56 + tmp_243 * ( -tmp_228 - tmp_230 - tmp_232 - tmp_234 - tmp_235 - tmp_236 - tmp_238 -
                                                      tmp_239 - tmp_240 + 1 ) ) +
          tmp_262 * ( -tmp_260 * tmp_56 + tmp_261 * ( -tmp_246 - tmp_248 - tmp_250 - tmp_252 - tmp_253 - tmp_254 - tmp_256 -
                                                      tmp_257 - tmp_258 + 1 ) ) +
          tmp_280 * ( -tmp_278 * tmp_56 + tmp_279 * ( -tmp_264 - tmp_266 - tmp_268 - tmp_270 - tmp_271 - tmp_272 - tmp_274 -
                                                      tmp_275 - tmp_276 + 1 ) ) +
          tmp_298 * ( -tmp_296 * tmp_56 + tmp_297 * ( -tmp_282 - tmp_284 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_292 -
                                                      tmp_293 - tmp_294 + 1 ) ) +
          tmp_316 * ( -tmp_314 * tmp_56 + tmp_315 * ( -tmp_300 - tmp_302 - tmp_304 - tmp_306 - tmp_307 - tmp_308 - tmp_310 -
                                                      tmp_311 - tmp_312 + 1 ) ) +
          tmp_334 * ( -tmp_332 * tmp_56 + tmp_333 * ( -tmp_318 - tmp_320 - tmp_322 - tmp_324 - tmp_325 - tmp_326 - tmp_328 -
                                                      tmp_329 - tmp_330 + 1 ) ) +
          tmp_352 * ( -tmp_350 * tmp_56 + tmp_351 * ( -tmp_336 - tmp_338 - tmp_340 - tmp_342 - tmp_343 - tmp_344 - tmp_346 -
                                                      tmp_347 - tmp_348 + 1 ) ) +
          tmp_370 * ( -tmp_368 * tmp_56 + tmp_369 * ( -tmp_354 - tmp_356 - tmp_358 - tmp_360 - tmp_361 - tmp_362 - tmp_364 -
                                                      tmp_365 - tmp_366 + 1 ) ) +
          tmp_388 * ( -tmp_386 * tmp_56 + tmp_387 * ( -tmp_372 - tmp_374 - tmp_376 - tmp_378 - tmp_379 - tmp_380 - tmp_382 -
                                                      tmp_383 - tmp_384 + 1 ) ) +
          tmp_406 * ( -tmp_404 * tmp_56 + tmp_405 * ( -tmp_390 - tmp_392 - tmp_394 - tmp_396 - tmp_397 - tmp_398 - tmp_400 -
                                                      tmp_401 - tmp_402 + 1 ) ) +
          tmp_424 * ( -tmp_422 * tmp_56 + tmp_423 * ( -tmp_408 - tmp_410 - tmp_412 - tmp_414 - tmp_415 - tmp_416 - tmp_418 -
                                                      tmp_419 - tmp_420 + 1 ) ) +
          tmp_64 * ( -tmp_55 * tmp_56 +
                     tmp_62 * ( -tmp_25 - tmp_32 - tmp_39 - tmp_42 - tmp_44 - tmp_46 - tmp_49 - tmp_51 - tmp_53 + 1 ) ) +
          tmp_82 * ( -tmp_56 * tmp_80 +
                     tmp_81 * ( -tmp_66 - tmp_68 - tmp_70 - tmp_72 - tmp_73 - tmp_74 - tmp_76 - tmp_77 - tmp_78 + 1 ) );
      real_t a_0_1 = tmp_100 * ( -tmp_425 * tmp_98 + tmp_89 * tmp_99 ) + tmp_118 * ( tmp_107 * tmp_117 - tmp_116 * tmp_425 ) +
                     tmp_136 * ( tmp_125 * tmp_135 - tmp_134 * tmp_425 ) + tmp_154 * ( tmp_143 * tmp_153 - tmp_152 * tmp_425 ) +
                     tmp_172 * ( tmp_161 * tmp_171 - tmp_170 * tmp_425 ) + tmp_190 * ( tmp_179 * tmp_189 - tmp_188 * tmp_425 ) +
                     tmp_208 * ( tmp_197 * tmp_207 - tmp_206 * tmp_425 ) + tmp_226 * ( tmp_215 * tmp_225 - tmp_224 * tmp_425 ) +
                     tmp_244 * ( tmp_233 * tmp_243 - tmp_242 * tmp_425 ) + tmp_262 * ( tmp_251 * tmp_261 - tmp_260 * tmp_425 ) +
                     tmp_280 * ( tmp_269 * tmp_279 - tmp_278 * tmp_425 ) + tmp_298 * ( tmp_287 * tmp_297 - tmp_296 * tmp_425 ) +
                     tmp_316 * ( tmp_305 * tmp_315 - tmp_314 * tmp_425 ) + tmp_334 * ( tmp_323 * tmp_333 - tmp_332 * tmp_425 ) +
                     tmp_352 * ( tmp_341 * tmp_351 - tmp_350 * tmp_425 ) + tmp_370 * ( tmp_359 * tmp_369 - tmp_368 * tmp_425 ) +
                     tmp_388 * ( tmp_377 * tmp_387 - tmp_386 * tmp_425 ) + tmp_406 * ( tmp_395 * tmp_405 - tmp_404 * tmp_425 ) +
                     tmp_424 * ( tmp_413 * tmp_423 - tmp_422 * tmp_425 ) + tmp_64 * ( tmp_40 * tmp_62 - tmp_425 * tmp_55 ) +
                     tmp_82 * ( -tmp_425 * tmp_80 + tmp_71 * tmp_81 );
      real_t a_0_2 = tmp_100 * ( -tmp_426 * tmp_98 + tmp_93 * tmp_99 ) + tmp_118 * ( tmp_111 * tmp_117 - tmp_116 * tmp_426 ) +
                     tmp_136 * ( tmp_129 * tmp_135 - tmp_134 * tmp_426 ) + tmp_154 * ( tmp_147 * tmp_153 - tmp_152 * tmp_426 ) +
                     tmp_172 * ( tmp_165 * tmp_171 - tmp_170 * tmp_426 ) + tmp_190 * ( tmp_183 * tmp_189 - tmp_188 * tmp_426 ) +
                     tmp_208 * ( tmp_201 * tmp_207 - tmp_206 * tmp_426 ) + tmp_226 * ( tmp_219 * tmp_225 - tmp_224 * tmp_426 ) +
                     tmp_244 * ( tmp_237 * tmp_243 - tmp_242 * tmp_426 ) + tmp_262 * ( tmp_255 * tmp_261 - tmp_260 * tmp_426 ) +
                     tmp_280 * ( tmp_273 * tmp_279 - tmp_278 * tmp_426 ) + tmp_298 * ( tmp_291 * tmp_297 - tmp_296 * tmp_426 ) +
                     tmp_316 * ( tmp_309 * tmp_315 - tmp_314 * tmp_426 ) + tmp_334 * ( tmp_327 * tmp_333 - tmp_332 * tmp_426 ) +
                     tmp_352 * ( tmp_345 * tmp_351 - tmp_350 * tmp_426 ) + tmp_370 * ( tmp_363 * tmp_369 - tmp_368 * tmp_426 ) +
                     tmp_388 * ( tmp_381 * tmp_387 - tmp_386 * tmp_426 ) + tmp_406 * ( tmp_399 * tmp_405 - tmp_404 * tmp_426 ) +
                     tmp_424 * ( tmp_417 * tmp_423 - tmp_422 * tmp_426 ) + tmp_64 * ( -tmp_426 * tmp_55 + tmp_47 * tmp_62 ) +
                     tmp_82 * ( -tmp_426 * tmp_80 + tmp_75 * tmp_81 );
      real_t a_0_3 = tmp_100 * ( -tmp_427 * tmp_98 + tmp_97 * tmp_99 ) + tmp_118 * ( tmp_115 * tmp_117 - tmp_116 * tmp_427 ) +
                     tmp_136 * ( tmp_133 * tmp_135 - tmp_134 * tmp_427 ) + tmp_154 * ( tmp_151 * tmp_153 - tmp_152 * tmp_427 ) +
                     tmp_172 * ( tmp_169 * tmp_171 - tmp_170 * tmp_427 ) + tmp_190 * ( tmp_187 * tmp_189 - tmp_188 * tmp_427 ) +
                     tmp_208 * ( tmp_205 * tmp_207 - tmp_206 * tmp_427 ) + tmp_226 * ( tmp_223 * tmp_225 - tmp_224 * tmp_427 ) +
                     tmp_244 * ( tmp_241 * tmp_243 - tmp_242 * tmp_427 ) + tmp_262 * ( tmp_259 * tmp_261 - tmp_260 * tmp_427 ) +
                     tmp_280 * ( tmp_277 * tmp_279 - tmp_278 * tmp_427 ) + tmp_298 * ( tmp_295 * tmp_297 - tmp_296 * tmp_427 ) +
                     tmp_316 * ( tmp_313 * tmp_315 - tmp_314 * tmp_427 ) + tmp_334 * ( tmp_331 * tmp_333 - tmp_332 * tmp_427 ) +
                     tmp_352 * ( tmp_349 * tmp_351 - tmp_350 * tmp_427 ) + tmp_370 * ( tmp_367 * tmp_369 - tmp_368 * tmp_427 ) +
                     tmp_388 * ( tmp_385 * tmp_387 - tmp_386 * tmp_427 ) + tmp_406 * ( tmp_403 * tmp_405 - tmp_404 * tmp_427 ) +
                     tmp_424 * ( tmp_421 * tmp_423 - tmp_422 * tmp_427 ) + tmp_64 * ( -tmp_427 * tmp_55 + tmp_54 * tmp_62 ) +
                     tmp_82 * ( -tmp_427 * tmp_80 + tmp_79 * tmp_81 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 0, 3 ) = a_0_3;
   }

   void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                  const std::vector< Point3D >& coordsElementOuter,
                                  const std::vector< Point3D >& coordsFacet,
                                  const Point3D&,
                                  const Point3D&,
                                  const Point3D&     outwardNormal,
                                  const DGBasisInfo& trialBasis,
                                  const DGBasisInfo& testBasis,
                                  int                trialDegree,
                                  int                testDegree,
                                  MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_1_1 + tmp_0;
      real_t tmp_2  = -p_affine_0_0;
      real_t tmp_3  = p_affine_2_0 + tmp_2;
      real_t tmp_4  = p_affine_3_1 + tmp_0;
      real_t tmp_5  = tmp_3 * tmp_4;
      real_t tmp_6  = p_affine_3_0 + tmp_2;
      real_t tmp_7  = p_affine_2_1 + tmp_0;
      real_t tmp_8  = tmp_6 * tmp_7;
      real_t tmp_9  = tmp_5 - tmp_8;
      real_t tmp_10 = p_affine_1_0 + tmp_2;
      real_t tmp_11 = -p_affine_0_2;
      real_t tmp_12 = p_affine_3_2 + tmp_11;
      real_t tmp_13 = tmp_12 * tmp_7;
      real_t tmp_14 = p_affine_1_2 + tmp_11;
      real_t tmp_15 = p_affine_2_2 + tmp_11;
      real_t tmp_16 = tmp_15 * tmp_6;
      real_t tmp_17 = tmp_15 * tmp_4;
      real_t tmp_18 = tmp_12 * tmp_3;
      real_t tmp_19 =
          1.0 / ( tmp_1 * tmp_16 - tmp_1 * tmp_18 + tmp_10 * tmp_13 - tmp_10 * tmp_17 + tmp_14 * tmp_5 - tmp_14 * tmp_8 );
      real_t tmp_20 = p_affine_8_2 + tmp_11;
      real_t tmp_21 = -p_affine_8_2;
      real_t tmp_22 = p_affine_9_2 + tmp_21;
      real_t tmp_23 = p_affine_10_2 + tmp_21;
      real_t tmp_24 = 0.031405749086161582 * tmp_22 + 0.93718850182767688 * tmp_23;
      real_t tmp_25 = tmp_19 * ( tmp_20 + tmp_24 );
      real_t tmp_26 = tmp_16 - tmp_18;
      real_t tmp_27 = p_affine_8_1 + tmp_0;
      real_t tmp_28 = -p_affine_8_1;
      real_t tmp_29 = p_affine_9_1 + tmp_28;
      real_t tmp_30 = p_affine_10_1 + tmp_28;
      real_t tmp_31 = 0.031405749086161582 * tmp_29 + 0.93718850182767688 * tmp_30;
      real_t tmp_32 = tmp_19 * ( tmp_27 + tmp_31 );
      real_t tmp_33 = tmp_13 - tmp_17;
      real_t tmp_34 = p_affine_8_0 + tmp_2;
      real_t tmp_35 = -p_affine_8_0;
      real_t tmp_36 = p_affine_9_0 + tmp_35;
      real_t tmp_37 = p_affine_10_0 + tmp_35;
      real_t tmp_38 = 0.031405749086161582 * tmp_36 + 0.93718850182767688 * tmp_37;
      real_t tmp_39 = tmp_19 * ( tmp_34 + tmp_38 );
      real_t tmp_40 = tmp_1 * tmp_6 - tmp_10 * tmp_4;
      real_t tmp_41 = tmp_10 * tmp_12 - tmp_14 * tmp_6;
      real_t tmp_42 = -tmp_1 * tmp_12 + tmp_14 * tmp_4;
      real_t tmp_43 = -tmp_1 * tmp_3 + tmp_10 * tmp_7;
      real_t tmp_44 = -tmp_10 * tmp_15 + tmp_14 * tmp_3;
      real_t tmp_45 = tmp_1 * tmp_15 - tmp_14 * tmp_7;
      real_t tmp_46 = tmp_1 * ( tmp_25 * tmp_9 + tmp_26 * tmp_32 + tmp_33 * tmp_39 - 1.0 / 4.0 ) +
                      tmp_4 * ( tmp_25 * tmp_43 + tmp_32 * tmp_44 + tmp_39 * tmp_45 - 1.0 / 4.0 ) +
                      tmp_7 * ( tmp_25 * tmp_40 + tmp_32 * tmp_41 + tmp_39 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_47 = -p_affine_4_1;
      real_t tmp_48 = p_affine_5_1 + tmp_47;
      real_t tmp_49 = -p_affine_4_2;
      real_t tmp_50 = p_affine_6_2 + tmp_49;
      real_t tmp_51 = tmp_48 * tmp_50;
      real_t tmp_52 = p_affine_6_1 + tmp_47;
      real_t tmp_53 = p_affine_5_2 + tmp_49;
      real_t tmp_54 = tmp_52 * tmp_53;
      real_t tmp_55 = -p_affine_4_0;
      real_t tmp_56 = p_affine_5_0 + tmp_55;
      real_t tmp_57 = p_affine_7_2 + tmp_49;
      real_t tmp_58 = tmp_52 * tmp_57;
      real_t tmp_59 = p_affine_6_0 + tmp_55;
      real_t tmp_60 = p_affine_7_1 + tmp_47;
      real_t tmp_61 = tmp_53 * tmp_60;
      real_t tmp_62 = p_affine_7_0 + tmp_55;
      real_t tmp_63 = tmp_50 * tmp_60;
      real_t tmp_64 = tmp_48 * tmp_57;
      real_t tmp_65 =
          1.0 / ( tmp_51 * tmp_62 - tmp_54 * tmp_62 + tmp_56 * tmp_58 - tmp_56 * tmp_63 + tmp_59 * tmp_61 - tmp_59 * tmp_64 );
      real_t tmp_66 = tmp_65 * ( tmp_51 - tmp_54 );
      real_t tmp_67 = tmp_65 * ( tmp_61 - tmp_64 );
      real_t tmp_68 = tmp_65 * ( tmp_58 - tmp_63 );
      real_t tmp_69 = tmp_65 * ( -tmp_50 * tmp_56 + tmp_53 * tmp_59 );
      real_t tmp_70 = tmp_65 * ( -tmp_53 * tmp_62 + tmp_56 * tmp_57 );
      real_t tmp_71 = tmp_65 * ( tmp_50 * tmp_62 - tmp_57 * tmp_59 );
      real_t tmp_72 = tmp_65 * ( -tmp_48 * tmp_59 + tmp_52 * tmp_56 );
      real_t tmp_73 = tmp_65 * ( tmp_48 * tmp_62 - tmp_56 * tmp_60 );
      real_t tmp_74 = tmp_65 * ( -tmp_52 * tmp_62 + tmp_59 * tmp_60 );
      real_t tmp_75 = 0.5 * p_affine_13_0 * ( -tmp_66 - tmp_67 - tmp_68 ) + 0.5 * p_affine_13_1 * ( -tmp_69 - tmp_70 - tmp_71 ) +
                      0.5 * p_affine_13_2 * ( -tmp_72 - tmp_73 - tmp_74 );
      real_t tmp_76 = p_affine_8_2 + tmp_49;
      real_t tmp_77 = tmp_24 + tmp_76;
      real_t tmp_78 = tmp_72 * tmp_77;
      real_t tmp_79 = tmp_73 * tmp_77;
      real_t tmp_80 = p_affine_8_1 + tmp_47;
      real_t tmp_81 = tmp_31 + tmp_80;
      real_t tmp_82 = tmp_69 * tmp_81;
      real_t tmp_83 = tmp_70 * tmp_81;
      real_t tmp_84 = tmp_74 * tmp_77;
      real_t tmp_85 = tmp_71 * tmp_81;
      real_t tmp_86 = p_affine_8_0 + tmp_55;
      real_t tmp_87 = tmp_38 + tmp_86;
      real_t tmp_88 = tmp_66 * tmp_87;
      real_t tmp_89 = tmp_67 * tmp_87;
      real_t tmp_90 = tmp_68 * tmp_87;
      real_t tmp_91 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_92 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_93 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_94 = ( std::abs( tmp_23 * tmp_91 - tmp_30 * tmp_93 ) * std::abs( tmp_23 * tmp_91 - tmp_30 * tmp_93 ) ) +
                      ( std::abs( tmp_23 * tmp_92 - tmp_37 * tmp_93 ) * std::abs( tmp_23 * tmp_92 - tmp_37 * tmp_93 ) ) +
                      ( std::abs( tmp_30 * tmp_92 - tmp_37 * tmp_91 ) * std::abs( tmp_30 * tmp_92 - tmp_37 * tmp_91 ) );
      real_t tmp_95  = 1.0 * std::pow( tmp_94, -0.25 );
      real_t tmp_96  = tmp_46 * tmp_95;
      real_t tmp_97  = 1.0 * std::pow( tmp_94, 1.0 / 2.0 );
      real_t tmp_98  = 0.0068572537431980923 * tmp_97;
      real_t tmp_99  = 0.19601935860219369 * tmp_22 + 0.60796128279561268 * tmp_23;
      real_t tmp_100 = tmp_19 * ( tmp_20 + tmp_99 );
      real_t tmp_101 = 0.19601935860219369 * tmp_29 + 0.60796128279561268 * tmp_30;
      real_t tmp_102 = tmp_19 * ( tmp_101 + tmp_27 );
      real_t tmp_103 = 0.19601935860219369 * tmp_36 + 0.60796128279561268 * tmp_37;
      real_t tmp_104 = tmp_19 * ( tmp_103 + tmp_34 );
      real_t tmp_105 = tmp_1 * ( tmp_100 * tmp_9 + tmp_102 * tmp_26 + tmp_104 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_100 * tmp_43 + tmp_102 * tmp_44 + tmp_104 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_100 * tmp_40 + tmp_102 * tmp_41 + tmp_104 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_106 = tmp_76 + tmp_99;
      real_t tmp_107 = tmp_106 * tmp_72;
      real_t tmp_108 = tmp_106 * tmp_73;
      real_t tmp_109 = tmp_101 + tmp_80;
      real_t tmp_110 = tmp_109 * tmp_69;
      real_t tmp_111 = tmp_109 * tmp_70;
      real_t tmp_112 = tmp_106 * tmp_74;
      real_t tmp_113 = tmp_109 * tmp_71;
      real_t tmp_114 = tmp_103 + tmp_86;
      real_t tmp_115 = tmp_114 * tmp_66;
      real_t tmp_116 = tmp_114 * tmp_67;
      real_t tmp_117 = tmp_114 * tmp_68;
      real_t tmp_118 = tmp_105 * tmp_95;
      real_t tmp_119 = 0.037198804536718075 * tmp_97;
      real_t tmp_120 = 0.37605877282253791 * tmp_22 + 0.039308471900058539 * tmp_23;
      real_t tmp_121 = tmp_19 * ( tmp_120 + tmp_20 );
      real_t tmp_122 = 0.37605877282253791 * tmp_29 + 0.039308471900058539 * tmp_30;
      real_t tmp_123 = tmp_19 * ( tmp_122 + tmp_27 );
      real_t tmp_124 = 0.37605877282253791 * tmp_36 + 0.039308471900058539 * tmp_37;
      real_t tmp_125 = tmp_19 * ( tmp_124 + tmp_34 );
      real_t tmp_126 = tmp_1 * ( tmp_121 * tmp_9 + tmp_123 * tmp_26 + tmp_125 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_121 * tmp_43 + tmp_123 * tmp_44 + tmp_125 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_121 * tmp_40 + tmp_123 * tmp_41 + tmp_125 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_127 = tmp_120 + tmp_76;
      real_t tmp_128 = tmp_127 * tmp_72;
      real_t tmp_129 = tmp_127 * tmp_73;
      real_t tmp_130 = tmp_122 + tmp_80;
      real_t tmp_131 = tmp_130 * tmp_69;
      real_t tmp_132 = tmp_130 * tmp_70;
      real_t tmp_133 = tmp_127 * tmp_74;
      real_t tmp_134 = tmp_130 * tmp_71;
      real_t tmp_135 = tmp_124 + tmp_86;
      real_t tmp_136 = tmp_135 * tmp_66;
      real_t tmp_137 = tmp_135 * tmp_67;
      real_t tmp_138 = tmp_135 * tmp_68;
      real_t tmp_139 = tmp_126 * tmp_95;
      real_t tmp_140 = 0.020848748529055869 * tmp_97;
      real_t tmp_141 = 0.78764240869137092 * tmp_22 + 0.1711304259088916 * tmp_23;
      real_t tmp_142 = tmp_19 * ( tmp_141 + tmp_20 );
      real_t tmp_143 = 0.78764240869137092 * tmp_29 + 0.1711304259088916 * tmp_30;
      real_t tmp_144 = tmp_19 * ( tmp_143 + tmp_27 );
      real_t tmp_145 = 0.78764240869137092 * tmp_36 + 0.1711304259088916 * tmp_37;
      real_t tmp_146 = tmp_19 * ( tmp_145 + tmp_34 );
      real_t tmp_147 = tmp_1 * ( tmp_142 * tmp_9 + tmp_144 * tmp_26 + tmp_146 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_142 * tmp_43 + tmp_144 * tmp_44 + tmp_146 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_142 * tmp_40 + tmp_144 * tmp_41 + tmp_146 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_148 = tmp_141 + tmp_76;
      real_t tmp_149 = tmp_148 * tmp_72;
      real_t tmp_150 = tmp_148 * tmp_73;
      real_t tmp_151 = tmp_143 + tmp_80;
      real_t tmp_152 = tmp_151 * tmp_69;
      real_t tmp_153 = tmp_151 * tmp_70;
      real_t tmp_154 = tmp_148 * tmp_74;
      real_t tmp_155 = tmp_151 * tmp_71;
      real_t tmp_156 = tmp_145 + tmp_86;
      real_t tmp_157 = tmp_156 * tmp_66;
      real_t tmp_158 = tmp_156 * tmp_67;
      real_t tmp_159 = tmp_156 * tmp_68;
      real_t tmp_160 = tmp_147 * tmp_95;
      real_t tmp_161 = 0.019202922745021479 * tmp_97;
      real_t tmp_162 = 0.58463275527740355 * tmp_22 + 0.37605877282253791 * tmp_23;
      real_t tmp_163 = tmp_19 * ( tmp_162 + tmp_20 );
      real_t tmp_164 = 0.58463275527740355 * tmp_29 + 0.37605877282253791 * tmp_30;
      real_t tmp_165 = tmp_19 * ( tmp_164 + tmp_27 );
      real_t tmp_166 = 0.58463275527740355 * tmp_36 + 0.37605877282253791 * tmp_37;
      real_t tmp_167 = tmp_19 * ( tmp_166 + tmp_34 );
      real_t tmp_168 = tmp_1 * ( tmp_163 * tmp_9 + tmp_165 * tmp_26 + tmp_167 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_163 * tmp_43 + tmp_165 * tmp_44 + tmp_167 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_163 * tmp_40 + tmp_165 * tmp_41 + tmp_167 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_169 = tmp_162 + tmp_76;
      real_t tmp_170 = tmp_169 * tmp_72;
      real_t tmp_171 = tmp_169 * tmp_73;
      real_t tmp_172 = tmp_164 + tmp_80;
      real_t tmp_173 = tmp_172 * tmp_69;
      real_t tmp_174 = tmp_172 * tmp_70;
      real_t tmp_175 = tmp_169 * tmp_74;
      real_t tmp_176 = tmp_172 * tmp_71;
      real_t tmp_177 = tmp_166 + tmp_86;
      real_t tmp_178 = tmp_177 * tmp_66;
      real_t tmp_179 = tmp_177 * tmp_67;
      real_t tmp_180 = tmp_177 * tmp_68;
      real_t tmp_181 = tmp_168 * tmp_95;
      real_t tmp_182 = 0.020848748529055869 * tmp_97;
      real_t tmp_183 = 0.041227165399737475 * tmp_22 + 0.78764240869137092 * tmp_23;
      real_t tmp_184 = tmp_19 * ( tmp_183 + tmp_20 );
      real_t tmp_185 = 0.041227165399737475 * tmp_29 + 0.78764240869137092 * tmp_30;
      real_t tmp_186 = tmp_19 * ( tmp_185 + tmp_27 );
      real_t tmp_187 = 0.041227165399737475 * tmp_36 + 0.78764240869137092 * tmp_37;
      real_t tmp_188 = tmp_19 * ( tmp_187 + tmp_34 );
      real_t tmp_189 = tmp_1 * ( tmp_184 * tmp_9 + tmp_186 * tmp_26 + tmp_188 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_184 * tmp_43 + tmp_186 * tmp_44 + tmp_188 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_184 * tmp_40 + tmp_186 * tmp_41 + tmp_188 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_190 = tmp_183 + tmp_76;
      real_t tmp_191 = tmp_190 * tmp_72;
      real_t tmp_192 = tmp_190 * tmp_73;
      real_t tmp_193 = tmp_185 + tmp_80;
      real_t tmp_194 = tmp_193 * tmp_69;
      real_t tmp_195 = tmp_193 * tmp_70;
      real_t tmp_196 = tmp_190 * tmp_74;
      real_t tmp_197 = tmp_193 * tmp_71;
      real_t tmp_198 = tmp_187 + tmp_86;
      real_t tmp_199 = tmp_198 * tmp_66;
      real_t tmp_200 = tmp_198 * tmp_67;
      real_t tmp_201 = tmp_198 * tmp_68;
      real_t tmp_202 = tmp_189 * tmp_95;
      real_t tmp_203 = 0.019202922745021479 * tmp_97;
      real_t tmp_204 = 0.039308471900058539 * tmp_22 + 0.58463275527740355 * tmp_23;
      real_t tmp_205 = tmp_19 * ( tmp_20 + tmp_204 );
      real_t tmp_206 = 0.039308471900058539 * tmp_29 + 0.58463275527740355 * tmp_30;
      real_t tmp_207 = tmp_19 * ( tmp_206 + tmp_27 );
      real_t tmp_208 = 0.039308471900058539 * tmp_36 + 0.58463275527740355 * tmp_37;
      real_t tmp_209 = tmp_19 * ( tmp_208 + tmp_34 );
      real_t tmp_210 = tmp_1 * ( tmp_205 * tmp_9 + tmp_207 * tmp_26 + tmp_209 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_205 * tmp_43 + tmp_207 * tmp_44 + tmp_209 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_205 * tmp_40 + tmp_207 * tmp_41 + tmp_209 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_211 = tmp_204 + tmp_76;
      real_t tmp_212 = tmp_211 * tmp_72;
      real_t tmp_213 = tmp_211 * tmp_73;
      real_t tmp_214 = tmp_206 + tmp_80;
      real_t tmp_215 = tmp_214 * tmp_69;
      real_t tmp_216 = tmp_214 * tmp_70;
      real_t tmp_217 = tmp_211 * tmp_74;
      real_t tmp_218 = tmp_214 * tmp_71;
      real_t tmp_219 = tmp_208 + tmp_86;
      real_t tmp_220 = tmp_219 * tmp_66;
      real_t tmp_221 = tmp_219 * tmp_67;
      real_t tmp_222 = tmp_219 * tmp_68;
      real_t tmp_223 = tmp_210 * tmp_95;
      real_t tmp_224 = 0.020848748529055869 * tmp_97;
      real_t tmp_225 = 0.78764240869137092 * tmp_22 + 0.041227165399737475 * tmp_23;
      real_t tmp_226 = tmp_19 * ( tmp_20 + tmp_225 );
      real_t tmp_227 = 0.78764240869137092 * tmp_29 + 0.041227165399737475 * tmp_30;
      real_t tmp_228 = tmp_19 * ( tmp_227 + tmp_27 );
      real_t tmp_229 = 0.78764240869137092 * tmp_36 + 0.041227165399737475 * tmp_37;
      real_t tmp_230 = tmp_19 * ( tmp_229 + tmp_34 );
      real_t tmp_231 = tmp_1 * ( tmp_226 * tmp_9 + tmp_228 * tmp_26 + tmp_230 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_226 * tmp_43 + tmp_228 * tmp_44 + tmp_230 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_226 * tmp_40 + tmp_228 * tmp_41 + tmp_230 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_232 = tmp_225 + tmp_76;
      real_t tmp_233 = tmp_232 * tmp_72;
      real_t tmp_234 = tmp_232 * tmp_73;
      real_t tmp_235 = tmp_227 + tmp_80;
      real_t tmp_236 = tmp_235 * tmp_69;
      real_t tmp_237 = tmp_235 * tmp_70;
      real_t tmp_238 = tmp_232 * tmp_74;
      real_t tmp_239 = tmp_235 * tmp_71;
      real_t tmp_240 = tmp_229 + tmp_86;
      real_t tmp_241 = tmp_240 * tmp_66;
      real_t tmp_242 = tmp_240 * tmp_67;
      real_t tmp_243 = tmp_240 * tmp_68;
      real_t tmp_244 = tmp_231 * tmp_95;
      real_t tmp_245 = 0.019202922745021479 * tmp_97;
      real_t tmp_246 = 0.58463275527740355 * tmp_22 + 0.039308471900058539 * tmp_23;
      real_t tmp_247 = tmp_19 * ( tmp_20 + tmp_246 );
      real_t tmp_248 = 0.58463275527740355 * tmp_29 + 0.039308471900058539 * tmp_30;
      real_t tmp_249 = tmp_19 * ( tmp_248 + tmp_27 );
      real_t tmp_250 = 0.58463275527740355 * tmp_36 + 0.039308471900058539 * tmp_37;
      real_t tmp_251 = tmp_19 * ( tmp_250 + tmp_34 );
      real_t tmp_252 = tmp_1 * ( tmp_247 * tmp_9 + tmp_249 * tmp_26 + tmp_251 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_247 * tmp_43 + tmp_249 * tmp_44 + tmp_251 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_247 * tmp_40 + tmp_249 * tmp_41 + tmp_251 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_253 = tmp_246 + tmp_76;
      real_t tmp_254 = tmp_253 * tmp_72;
      real_t tmp_255 = tmp_253 * tmp_73;
      real_t tmp_256 = tmp_248 + tmp_80;
      real_t tmp_257 = tmp_256 * tmp_69;
      real_t tmp_258 = tmp_256 * tmp_70;
      real_t tmp_259 = tmp_253 * tmp_74;
      real_t tmp_260 = tmp_256 * tmp_71;
      real_t tmp_261 = tmp_250 + tmp_86;
      real_t tmp_262 = tmp_261 * tmp_66;
      real_t tmp_263 = tmp_261 * tmp_67;
      real_t tmp_264 = tmp_261 * tmp_68;
      real_t tmp_265 = tmp_252 * tmp_95;
      real_t tmp_266 = 0.020848748529055869 * tmp_97;
      real_t tmp_267 = 0.1711304259088916 * tmp_22 + 0.78764240869137092 * tmp_23;
      real_t tmp_268 = tmp_19 * ( tmp_20 + tmp_267 );
      real_t tmp_269 = 0.1711304259088916 * tmp_29 + 0.78764240869137092 * tmp_30;
      real_t tmp_270 = tmp_19 * ( tmp_269 + tmp_27 );
      real_t tmp_271 = 0.1711304259088916 * tmp_36 + 0.78764240869137092 * tmp_37;
      real_t tmp_272 = tmp_19 * ( tmp_271 + tmp_34 );
      real_t tmp_273 = tmp_1 * ( tmp_26 * tmp_270 + tmp_268 * tmp_9 + tmp_272 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_268 * tmp_43 + tmp_270 * tmp_44 + tmp_272 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_268 * tmp_40 + tmp_270 * tmp_41 + tmp_272 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_274 = tmp_267 + tmp_76;
      real_t tmp_275 = tmp_274 * tmp_72;
      real_t tmp_276 = tmp_274 * tmp_73;
      real_t tmp_277 = tmp_269 + tmp_80;
      real_t tmp_278 = tmp_277 * tmp_69;
      real_t tmp_279 = tmp_277 * tmp_70;
      real_t tmp_280 = tmp_274 * tmp_74;
      real_t tmp_281 = tmp_277 * tmp_71;
      real_t tmp_282 = tmp_271 + tmp_86;
      real_t tmp_283 = tmp_282 * tmp_66;
      real_t tmp_284 = tmp_282 * tmp_67;
      real_t tmp_285 = tmp_282 * tmp_68;
      real_t tmp_286 = tmp_273 * tmp_95;
      real_t tmp_287 = 0.019202922745021479 * tmp_97;
      real_t tmp_288 = 0.37605877282253791 * tmp_22 + 0.58463275527740355 * tmp_23;
      real_t tmp_289 = tmp_19 * ( tmp_20 + tmp_288 );
      real_t tmp_290 = 0.37605877282253791 * tmp_29 + 0.58463275527740355 * tmp_30;
      real_t tmp_291 = tmp_19 * ( tmp_27 + tmp_290 );
      real_t tmp_292 = 0.37605877282253791 * tmp_36 + 0.58463275527740355 * tmp_37;
      real_t tmp_293 = tmp_19 * ( tmp_292 + tmp_34 );
      real_t tmp_294 = tmp_1 * ( tmp_26 * tmp_291 + tmp_289 * tmp_9 + tmp_293 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_289 * tmp_43 + tmp_291 * tmp_44 + tmp_293 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_289 * tmp_40 + tmp_291 * tmp_41 + tmp_293 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_295 = tmp_288 + tmp_76;
      real_t tmp_296 = tmp_295 * tmp_72;
      real_t tmp_297 = tmp_295 * tmp_73;
      real_t tmp_298 = tmp_290 + tmp_80;
      real_t tmp_299 = tmp_298 * tmp_69;
      real_t tmp_300 = tmp_298 * tmp_70;
      real_t tmp_301 = tmp_295 * tmp_74;
      real_t tmp_302 = tmp_298 * tmp_71;
      real_t tmp_303 = tmp_292 + tmp_86;
      real_t tmp_304 = tmp_303 * tmp_66;
      real_t tmp_305 = tmp_303 * tmp_67;
      real_t tmp_306 = tmp_303 * tmp_68;
      real_t tmp_307 = tmp_294 * tmp_95;
      real_t tmp_308 = 0.020848748529055869 * tmp_97;
      real_t tmp_309 = 0.041227165399737475 * tmp_22 + 0.1711304259088916 * tmp_23;
      real_t tmp_310 = tmp_19 * ( tmp_20 + tmp_309 );
      real_t tmp_311 = 0.041227165399737475 * tmp_29 + 0.1711304259088916 * tmp_30;
      real_t tmp_312 = tmp_19 * ( tmp_27 + tmp_311 );
      real_t tmp_313 = 0.041227165399737475 * tmp_36 + 0.1711304259088916 * tmp_37;
      real_t tmp_314 = tmp_19 * ( tmp_313 + tmp_34 );
      real_t tmp_315 = tmp_1 * ( tmp_26 * tmp_312 + tmp_310 * tmp_9 + tmp_314 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_310 * tmp_43 + tmp_312 * tmp_44 + tmp_314 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_310 * tmp_40 + tmp_312 * tmp_41 + tmp_314 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_316 = tmp_309 + tmp_76;
      real_t tmp_317 = tmp_316 * tmp_72;
      real_t tmp_318 = tmp_316 * tmp_73;
      real_t tmp_319 = tmp_311 + tmp_80;
      real_t tmp_320 = tmp_319 * tmp_69;
      real_t tmp_321 = tmp_319 * tmp_70;
      real_t tmp_322 = tmp_316 * tmp_74;
      real_t tmp_323 = tmp_319 * tmp_71;
      real_t tmp_324 = tmp_313 + tmp_86;
      real_t tmp_325 = tmp_324 * tmp_66;
      real_t tmp_326 = tmp_324 * tmp_67;
      real_t tmp_327 = tmp_324 * tmp_68;
      real_t tmp_328 = tmp_315 * tmp_95;
      real_t tmp_329 = 0.019202922745021479 * tmp_97;
      real_t tmp_330 = 0.40446199974765351 * tmp_22 + 0.19107600050469298 * tmp_23;
      real_t tmp_331 = tmp_19 * ( tmp_20 + tmp_330 );
      real_t tmp_332 = 0.40446199974765351 * tmp_29 + 0.19107600050469298 * tmp_30;
      real_t tmp_333 = tmp_19 * ( tmp_27 + tmp_332 );
      real_t tmp_334 = 0.40446199974765351 * tmp_36 + 0.19107600050469298 * tmp_37;
      real_t tmp_335 = tmp_19 * ( tmp_334 + tmp_34 );
      real_t tmp_336 = tmp_1 * ( tmp_26 * tmp_333 + tmp_33 * tmp_335 + tmp_331 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_331 * tmp_43 + tmp_333 * tmp_44 + tmp_335 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_331 * tmp_40 + tmp_333 * tmp_41 + tmp_335 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_337 = tmp_330 + tmp_76;
      real_t tmp_338 = tmp_337 * tmp_72;
      real_t tmp_339 = tmp_337 * tmp_73;
      real_t tmp_340 = tmp_332 + tmp_80;
      real_t tmp_341 = tmp_340 * tmp_69;
      real_t tmp_342 = tmp_340 * tmp_70;
      real_t tmp_343 = tmp_337 * tmp_74;
      real_t tmp_344 = tmp_340 * tmp_71;
      real_t tmp_345 = tmp_334 + tmp_86;
      real_t tmp_346 = tmp_345 * tmp_66;
      real_t tmp_347 = tmp_345 * tmp_67;
      real_t tmp_348 = tmp_345 * tmp_68;
      real_t tmp_349 = tmp_336 * tmp_95;
      real_t tmp_350 = 0.042507265838595799 * tmp_97;
      real_t tmp_351 = 0.039308471900058539 * tmp_22 + 0.37605877282253791 * tmp_23;
      real_t tmp_352 = tmp_19 * ( tmp_20 + tmp_351 );
      real_t tmp_353 = 0.039308471900058539 * tmp_29 + 0.37605877282253791 * tmp_30;
      real_t tmp_354 = tmp_19 * ( tmp_27 + tmp_353 );
      real_t tmp_355 = 0.039308471900058539 * tmp_36 + 0.37605877282253791 * tmp_37;
      real_t tmp_356 = tmp_19 * ( tmp_34 + tmp_355 );
      real_t tmp_357 = tmp_1 * ( tmp_26 * tmp_354 + tmp_33 * tmp_356 + tmp_352 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_352 * tmp_43 + tmp_354 * tmp_44 + tmp_356 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_352 * tmp_40 + tmp_354 * tmp_41 + tmp_356 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_358 = tmp_351 + tmp_76;
      real_t tmp_359 = tmp_358 * tmp_72;
      real_t tmp_360 = tmp_358 * tmp_73;
      real_t tmp_361 = tmp_353 + tmp_80;
      real_t tmp_362 = tmp_361 * tmp_69;
      real_t tmp_363 = tmp_361 * tmp_70;
      real_t tmp_364 = tmp_358 * tmp_74;
      real_t tmp_365 = tmp_361 * tmp_71;
      real_t tmp_366 = tmp_355 + tmp_86;
      real_t tmp_367 = tmp_366 * tmp_66;
      real_t tmp_368 = tmp_366 * tmp_67;
      real_t tmp_369 = tmp_366 * tmp_68;
      real_t tmp_370 = tmp_357 * tmp_95;
      real_t tmp_371 = 0.020848748529055869 * tmp_97;
      real_t tmp_372 = 0.93718850182767688 * tmp_22 + 0.031405749086161582 * tmp_23;
      real_t tmp_373 = tmp_19 * ( tmp_20 + tmp_372 );
      real_t tmp_374 = 0.93718850182767688 * tmp_29 + 0.031405749086161582 * tmp_30;
      real_t tmp_375 = tmp_19 * ( tmp_27 + tmp_374 );
      real_t tmp_376 = 0.93718850182767688 * tmp_36 + 0.031405749086161582 * tmp_37;
      real_t tmp_377 = tmp_19 * ( tmp_34 + tmp_376 );
      real_t tmp_378 = tmp_1 * ( tmp_26 * tmp_375 + tmp_33 * tmp_377 + tmp_373 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_373 * tmp_43 + tmp_375 * tmp_44 + tmp_377 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_373 * tmp_40 + tmp_375 * tmp_41 + tmp_377 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_379 = tmp_372 + tmp_76;
      real_t tmp_380 = tmp_379 * tmp_72;
      real_t tmp_381 = tmp_379 * tmp_73;
      real_t tmp_382 = tmp_374 + tmp_80;
      real_t tmp_383 = tmp_382 * tmp_69;
      real_t tmp_384 = tmp_382 * tmp_70;
      real_t tmp_385 = tmp_379 * tmp_74;
      real_t tmp_386 = tmp_382 * tmp_71;
      real_t tmp_387 = tmp_376 + tmp_86;
      real_t tmp_388 = tmp_387 * tmp_66;
      real_t tmp_389 = tmp_387 * tmp_67;
      real_t tmp_390 = tmp_387 * tmp_68;
      real_t tmp_391 = tmp_378 * tmp_95;
      real_t tmp_392 = 0.0068572537431980923 * tmp_97;
      real_t tmp_393 = 0.60796128279561268 * tmp_22 + 0.19601935860219369 * tmp_23;
      real_t tmp_394 = tmp_19 * ( tmp_20 + tmp_393 );
      real_t tmp_395 = 0.60796128279561268 * tmp_29 + 0.19601935860219369 * tmp_30;
      real_t tmp_396 = tmp_19 * ( tmp_27 + tmp_395 );
      real_t tmp_397 = 0.60796128279561268 * tmp_36 + 0.19601935860219369 * tmp_37;
      real_t tmp_398 = tmp_19 * ( tmp_34 + tmp_397 );
      real_t tmp_399 = tmp_1 * ( tmp_26 * tmp_396 + tmp_33 * tmp_398 + tmp_394 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_394 * tmp_43 + tmp_396 * tmp_44 + tmp_398 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_394 * tmp_40 + tmp_396 * tmp_41 + tmp_398 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_400 = tmp_393 + tmp_76;
      real_t tmp_401 = tmp_400 * tmp_72;
      real_t tmp_402 = tmp_400 * tmp_73;
      real_t tmp_403 = tmp_395 + tmp_80;
      real_t tmp_404 = tmp_403 * tmp_69;
      real_t tmp_405 = tmp_403 * tmp_70;
      real_t tmp_406 = tmp_400 * tmp_74;
      real_t tmp_407 = tmp_403 * tmp_71;
      real_t tmp_408 = tmp_397 + tmp_86;
      real_t tmp_409 = tmp_408 * tmp_66;
      real_t tmp_410 = tmp_408 * tmp_67;
      real_t tmp_411 = tmp_408 * tmp_68;
      real_t tmp_412 = tmp_399 * tmp_95;
      real_t tmp_413 = 0.037198804536718075 * tmp_97;
      real_t tmp_414 = 0.19107600050469298 * tmp_22 + 0.40446199974765351 * tmp_23;
      real_t tmp_415 = tmp_19 * ( tmp_20 + tmp_414 );
      real_t tmp_416 = 0.19107600050469298 * tmp_29 + 0.40446199974765351 * tmp_30;
      real_t tmp_417 = tmp_19 * ( tmp_27 + tmp_416 );
      real_t tmp_418 = 0.19107600050469298 * tmp_36 + 0.40446199974765351 * tmp_37;
      real_t tmp_419 = tmp_19 * ( tmp_34 + tmp_418 );
      real_t tmp_420 = tmp_1 * ( tmp_26 * tmp_417 + tmp_33 * tmp_419 + tmp_415 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_415 * tmp_43 + tmp_417 * tmp_44 + tmp_419 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_40 * tmp_415 + tmp_41 * tmp_417 + tmp_419 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_421 = tmp_414 + tmp_76;
      real_t tmp_422 = tmp_421 * tmp_72;
      real_t tmp_423 = tmp_421 * tmp_73;
      real_t tmp_424 = tmp_416 + tmp_80;
      real_t tmp_425 = tmp_424 * tmp_69;
      real_t tmp_426 = tmp_424 * tmp_70;
      real_t tmp_427 = tmp_421 * tmp_74;
      real_t tmp_428 = tmp_424 * tmp_71;
      real_t tmp_429 = tmp_418 + tmp_86;
      real_t tmp_430 = tmp_429 * tmp_66;
      real_t tmp_431 = tmp_429 * tmp_67;
      real_t tmp_432 = tmp_429 * tmp_68;
      real_t tmp_433 = tmp_420 * tmp_95;
      real_t tmp_434 = 0.042507265838595799 * tmp_97;
      real_t tmp_435 = 0.031405749086161582 * tmp_22 + 0.031405749086161582 * tmp_23;
      real_t tmp_436 = tmp_19 * ( tmp_20 + tmp_435 );
      real_t tmp_437 = 0.031405749086161582 * tmp_29 + 0.031405749086161582 * tmp_30;
      real_t tmp_438 = tmp_19 * ( tmp_27 + tmp_437 );
      real_t tmp_439 = 0.031405749086161582 * tmp_36 + 0.031405749086161582 * tmp_37;
      real_t tmp_440 = tmp_19 * ( tmp_34 + tmp_439 );
      real_t tmp_441 = tmp_1 * ( tmp_26 * tmp_438 + tmp_33 * tmp_440 + tmp_436 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_43 * tmp_436 + tmp_438 * tmp_44 + tmp_440 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_40 * tmp_436 + tmp_41 * tmp_438 + tmp_42 * tmp_440 - 1.0 / 4.0 );
      real_t tmp_442 = tmp_435 + tmp_76;
      real_t tmp_443 = tmp_442 * tmp_72;
      real_t tmp_444 = tmp_442 * tmp_73;
      real_t tmp_445 = tmp_437 + tmp_80;
      real_t tmp_446 = tmp_445 * tmp_69;
      real_t tmp_447 = tmp_445 * tmp_70;
      real_t tmp_448 = tmp_442 * tmp_74;
      real_t tmp_449 = tmp_445 * tmp_71;
      real_t tmp_450 = tmp_439 + tmp_86;
      real_t tmp_451 = tmp_450 * tmp_66;
      real_t tmp_452 = tmp_450 * tmp_67;
      real_t tmp_453 = tmp_450 * tmp_68;
      real_t tmp_454 = tmp_441 * tmp_95;
      real_t tmp_455 = 0.0068572537431980923 * tmp_97;
      real_t tmp_456 = 0.19601935860219369 * tmp_22 + 0.19601935860219369 * tmp_23;
      real_t tmp_457 = tmp_19 * ( tmp_20 + tmp_456 );
      real_t tmp_458 = 0.19601935860219369 * tmp_29 + 0.19601935860219369 * tmp_30;
      real_t tmp_459 = tmp_19 * ( tmp_27 + tmp_458 );
      real_t tmp_460 = 0.19601935860219369 * tmp_36 + 0.19601935860219369 * tmp_37;
      real_t tmp_461 = tmp_19 * ( tmp_34 + tmp_460 );
      real_t tmp_462 = tmp_1 * ( tmp_26 * tmp_459 + tmp_33 * tmp_461 + tmp_457 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_43 * tmp_457 + tmp_44 * tmp_459 + tmp_45 * tmp_461 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_40 * tmp_457 + tmp_41 * tmp_459 + tmp_42 * tmp_461 - 1.0 / 4.0 );
      real_t tmp_463 = tmp_456 + tmp_76;
      real_t tmp_464 = tmp_463 * tmp_72;
      real_t tmp_465 = tmp_463 * tmp_73;
      real_t tmp_466 = tmp_458 + tmp_80;
      real_t tmp_467 = tmp_466 * tmp_69;
      real_t tmp_468 = tmp_466 * tmp_70;
      real_t tmp_469 = tmp_463 * tmp_74;
      real_t tmp_470 = tmp_466 * tmp_71;
      real_t tmp_471 = tmp_460 + tmp_86;
      real_t tmp_472 = tmp_471 * tmp_66;
      real_t tmp_473 = tmp_471 * tmp_67;
      real_t tmp_474 = tmp_471 * tmp_68;
      real_t tmp_475 = tmp_462 * tmp_95;
      real_t tmp_476 = 0.037198804536718075 * tmp_97;
      real_t tmp_477 = 0.40446199974765351 * tmp_22 + 0.40446199974765351 * tmp_23;
      real_t tmp_478 = tmp_19 * ( tmp_20 + tmp_477 );
      real_t tmp_479 = 0.40446199974765351 * tmp_29 + 0.40446199974765351 * tmp_30;
      real_t tmp_480 = tmp_19 * ( tmp_27 + tmp_479 );
      real_t tmp_481 = 0.40446199974765351 * tmp_36 + 0.40446199974765351 * tmp_37;
      real_t tmp_482 = tmp_19 * ( tmp_34 + tmp_481 );
      real_t tmp_483 = tmp_1 * ( tmp_26 * tmp_480 + tmp_33 * tmp_482 + tmp_478 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_43 * tmp_478 + tmp_44 * tmp_480 + tmp_45 * tmp_482 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_40 * tmp_478 + tmp_41 * tmp_480 + tmp_42 * tmp_482 - 1.0 / 4.0 );
      real_t tmp_484 = tmp_477 + tmp_76;
      real_t tmp_485 = tmp_484 * tmp_72;
      real_t tmp_486 = tmp_484 * tmp_73;
      real_t tmp_487 = tmp_479 + tmp_80;
      real_t tmp_488 = tmp_487 * tmp_69;
      real_t tmp_489 = tmp_487 * tmp_70;
      real_t tmp_490 = tmp_484 * tmp_74;
      real_t tmp_491 = tmp_487 * tmp_71;
      real_t tmp_492 = tmp_481 + tmp_86;
      real_t tmp_493 = tmp_492 * tmp_66;
      real_t tmp_494 = tmp_492 * tmp_67;
      real_t tmp_495 = tmp_492 * tmp_68;
      real_t tmp_496 = tmp_483 * tmp_95;
      real_t tmp_497 = 0.042507265838595799 * tmp_97;
      real_t tmp_498 = 0.1711304259088916 * tmp_22 + 0.041227165399737475 * tmp_23;
      real_t tmp_499 = tmp_19 * ( tmp_20 + tmp_498 );
      real_t tmp_500 = 0.1711304259088916 * tmp_29 + 0.041227165399737475 * tmp_30;
      real_t tmp_501 = tmp_19 * ( tmp_27 + tmp_500 );
      real_t tmp_502 = 0.1711304259088916 * tmp_36 + 0.041227165399737475 * tmp_37;
      real_t tmp_503 = tmp_19 * ( tmp_34 + tmp_502 );
      real_t tmp_504 = tmp_1 * ( tmp_26 * tmp_501 + tmp_33 * tmp_503 + tmp_499 * tmp_9 - 1.0 / 4.0 ) +
                       tmp_4 * ( tmp_43 * tmp_499 + tmp_44 * tmp_501 + tmp_45 * tmp_503 - 1.0 / 4.0 ) +
                       tmp_7 * ( tmp_40 * tmp_499 + tmp_41 * tmp_501 + tmp_42 * tmp_503 - 1.0 / 4.0 );
      real_t tmp_505 = tmp_498 + tmp_76;
      real_t tmp_506 = tmp_505 * tmp_72;
      real_t tmp_507 = tmp_505 * tmp_73;
      real_t tmp_508 = tmp_500 + tmp_80;
      real_t tmp_509 = tmp_508 * tmp_69;
      real_t tmp_510 = tmp_508 * tmp_70;
      real_t tmp_511 = tmp_505 * tmp_74;
      real_t tmp_512 = tmp_508 * tmp_71;
      real_t tmp_513 = tmp_502 + tmp_86;
      real_t tmp_514 = tmp_513 * tmp_66;
      real_t tmp_515 = tmp_513 * tmp_67;
      real_t tmp_516 = tmp_513 * tmp_68;
      real_t tmp_517 = tmp_504 * tmp_95;
      real_t tmp_518 = 0.019202922745021479 * tmp_97;
      real_t tmp_519 = 0.5 * p_affine_13_0 * tmp_68 + 0.5 * p_affine_13_1 * tmp_71 + 0.5 * p_affine_13_2 * tmp_74;
      real_t tmp_520 = 0.5 * p_affine_13_0 * tmp_67 + 0.5 * p_affine_13_1 * tmp_70 + 0.5 * p_affine_13_2 * tmp_73;
      real_t tmp_521 = 0.5 * p_affine_13_0 * tmp_66 + 0.5 * p_affine_13_1 * tmp_69 + 0.5 * p_affine_13_2 * tmp_72;
      real_t a_0_0   = tmp_119 * ( -tmp_105 * tmp_75 - tmp_118 * ( -tmp_107 - tmp_108 - tmp_110 - tmp_111 - tmp_112 - tmp_113 -
                                                                 tmp_115 - tmp_116 - tmp_117 + 1 ) ) +
                     tmp_140 * ( -tmp_126 * tmp_75 - tmp_139 * ( -tmp_128 - tmp_129 - tmp_131 - tmp_132 - tmp_133 - tmp_134 -
                                                                 tmp_136 - tmp_137 - tmp_138 + 1 ) ) +
                     tmp_161 * ( -tmp_147 * tmp_75 - tmp_160 * ( -tmp_149 - tmp_150 - tmp_152 - tmp_153 - tmp_154 - tmp_155 -
                                                                 tmp_157 - tmp_158 - tmp_159 + 1 ) ) +
                     tmp_182 * ( -tmp_168 * tmp_75 - tmp_181 * ( -tmp_170 - tmp_171 - tmp_173 - tmp_174 - tmp_175 - tmp_176 -
                                                                 tmp_178 - tmp_179 - tmp_180 + 1 ) ) +
                     tmp_203 * ( -tmp_189 * tmp_75 - tmp_202 * ( -tmp_191 - tmp_192 - tmp_194 - tmp_195 - tmp_196 - tmp_197 -
                                                                 tmp_199 - tmp_200 - tmp_201 + 1 ) ) +
                     tmp_224 * ( -tmp_210 * tmp_75 - tmp_223 * ( -tmp_212 - tmp_213 - tmp_215 - tmp_216 - tmp_217 - tmp_218 -
                                                                 tmp_220 - tmp_221 - tmp_222 + 1 ) ) +
                     tmp_245 * ( -tmp_231 * tmp_75 - tmp_244 * ( -tmp_233 - tmp_234 - tmp_236 - tmp_237 - tmp_238 - tmp_239 -
                                                                 tmp_241 - tmp_242 - tmp_243 + 1 ) ) +
                     tmp_266 * ( -tmp_252 * tmp_75 - tmp_265 * ( -tmp_254 - tmp_255 - tmp_257 - tmp_258 - tmp_259 - tmp_260 -
                                                                 tmp_262 - tmp_263 - tmp_264 + 1 ) ) +
                     tmp_287 * ( -tmp_273 * tmp_75 - tmp_286 * ( -tmp_275 - tmp_276 - tmp_278 - tmp_279 - tmp_280 - tmp_281 -
                                                                 tmp_283 - tmp_284 - tmp_285 + 1 ) ) +
                     tmp_308 * ( -tmp_294 * tmp_75 - tmp_307 * ( -tmp_296 - tmp_297 - tmp_299 - tmp_300 - tmp_301 - tmp_302 -
                                                                 tmp_304 - tmp_305 - tmp_306 + 1 ) ) +
                     tmp_329 * ( -tmp_315 * tmp_75 - tmp_328 * ( -tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 -
                                                                 tmp_325 - tmp_326 - tmp_327 + 1 ) ) +
                     tmp_350 * ( -tmp_336 * tmp_75 - tmp_349 * ( -tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 - tmp_344 -
                                                                 tmp_346 - tmp_347 - tmp_348 + 1 ) ) +
                     tmp_371 * ( -tmp_357 * tmp_75 - tmp_370 * ( -tmp_359 - tmp_360 - tmp_362 - tmp_363 - tmp_364 - tmp_365 -
                                                                 tmp_367 - tmp_368 - tmp_369 + 1 ) ) +
                     tmp_392 * ( -tmp_378 * tmp_75 - tmp_391 * ( -tmp_380 - tmp_381 - tmp_383 - tmp_384 - tmp_385 - tmp_386 -
                                                                 tmp_388 - tmp_389 - tmp_390 + 1 ) ) +
                     tmp_413 * ( -tmp_399 * tmp_75 - tmp_412 * ( -tmp_401 - tmp_402 - tmp_404 - tmp_405 - tmp_406 - tmp_407 -
                                                                 tmp_409 - tmp_410 - tmp_411 + 1 ) ) +
                     tmp_434 * ( -tmp_420 * tmp_75 - tmp_433 * ( -tmp_422 - tmp_423 - tmp_425 - tmp_426 - tmp_427 - tmp_428 -
                                                                 tmp_430 - tmp_431 - tmp_432 + 1 ) ) +
                     tmp_455 * ( -tmp_441 * tmp_75 - tmp_454 * ( -tmp_443 - tmp_444 - tmp_446 - tmp_447 - tmp_448 - tmp_449 -
                                                                 tmp_451 - tmp_452 - tmp_453 + 1 ) ) +
                     tmp_476 * ( -tmp_462 * tmp_75 - tmp_475 * ( -tmp_464 - tmp_465 - tmp_467 - tmp_468 - tmp_469 - tmp_470 -
                                                                 tmp_472 - tmp_473 - tmp_474 + 1 ) ) +
                     tmp_497 * ( -tmp_483 * tmp_75 - tmp_496 * ( -tmp_485 - tmp_486 - tmp_488 - tmp_489 - tmp_490 - tmp_491 -
                                                                 tmp_493 - tmp_494 - tmp_495 + 1 ) ) +
                     tmp_518 * ( -tmp_504 * tmp_75 - tmp_517 * ( -tmp_506 - tmp_507 - tmp_509 - tmp_510 - tmp_511 - tmp_512 -
                                                                 tmp_514 - tmp_515 - tmp_516 + 1 ) ) +
                     tmp_98 * ( -tmp_46 * tmp_75 - tmp_96 * ( -tmp_78 - tmp_79 - tmp_82 - tmp_83 - tmp_84 - tmp_85 - tmp_88 -
                                                              tmp_89 - tmp_90 + 1 ) );
      real_t a_0_1 = tmp_119 * ( -tmp_105 * tmp_519 - tmp_118 * ( tmp_112 + tmp_113 + tmp_117 ) ) +
                     tmp_140 * ( -tmp_126 * tmp_519 - tmp_139 * ( tmp_133 + tmp_134 + tmp_138 ) ) +
                     tmp_161 * ( -tmp_147 * tmp_519 - tmp_160 * ( tmp_154 + tmp_155 + tmp_159 ) ) +
                     tmp_182 * ( -tmp_168 * tmp_519 - tmp_181 * ( tmp_175 + tmp_176 + tmp_180 ) ) +
                     tmp_203 * ( -tmp_189 * tmp_519 - tmp_202 * ( tmp_196 + tmp_197 + tmp_201 ) ) +
                     tmp_224 * ( -tmp_210 * tmp_519 - tmp_223 * ( tmp_217 + tmp_218 + tmp_222 ) ) +
                     tmp_245 * ( -tmp_231 * tmp_519 - tmp_244 * ( tmp_238 + tmp_239 + tmp_243 ) ) +
                     tmp_266 * ( -tmp_252 * tmp_519 - tmp_265 * ( tmp_259 + tmp_260 + tmp_264 ) ) +
                     tmp_287 * ( -tmp_273 * tmp_519 - tmp_286 * ( tmp_280 + tmp_281 + tmp_285 ) ) +
                     tmp_308 * ( -tmp_294 * tmp_519 - tmp_307 * ( tmp_301 + tmp_302 + tmp_306 ) ) +
                     tmp_329 * ( -tmp_315 * tmp_519 - tmp_328 * ( tmp_322 + tmp_323 + tmp_327 ) ) +
                     tmp_350 * ( -tmp_336 * tmp_519 - tmp_349 * ( tmp_343 + tmp_344 + tmp_348 ) ) +
                     tmp_371 * ( -tmp_357 * tmp_519 - tmp_370 * ( tmp_364 + tmp_365 + tmp_369 ) ) +
                     tmp_392 * ( -tmp_378 * tmp_519 - tmp_391 * ( tmp_385 + tmp_386 + tmp_390 ) ) +
                     tmp_413 * ( -tmp_399 * tmp_519 - tmp_412 * ( tmp_406 + tmp_407 + tmp_411 ) ) +
                     tmp_434 * ( -tmp_420 * tmp_519 - tmp_433 * ( tmp_427 + tmp_428 + tmp_432 ) ) +
                     tmp_455 * ( -tmp_441 * tmp_519 - tmp_454 * ( tmp_448 + tmp_449 + tmp_453 ) ) +
                     tmp_476 * ( -tmp_462 * tmp_519 - tmp_475 * ( tmp_469 + tmp_470 + tmp_474 ) ) +
                     tmp_497 * ( -tmp_483 * tmp_519 - tmp_496 * ( tmp_490 + tmp_491 + tmp_495 ) ) +
                     tmp_518 * ( -tmp_504 * tmp_519 - tmp_517 * ( tmp_511 + tmp_512 + tmp_516 ) ) +
                     tmp_98 * ( -tmp_46 * tmp_519 - tmp_96 * ( tmp_84 + tmp_85 + tmp_90 ) );
      real_t a_0_2 = tmp_119 * ( -tmp_105 * tmp_520 - tmp_118 * ( tmp_108 + tmp_111 + tmp_116 ) ) +
                     tmp_140 * ( -tmp_126 * tmp_520 - tmp_139 * ( tmp_129 + tmp_132 + tmp_137 ) ) +
                     tmp_161 * ( -tmp_147 * tmp_520 - tmp_160 * ( tmp_150 + tmp_153 + tmp_158 ) ) +
                     tmp_182 * ( -tmp_168 * tmp_520 - tmp_181 * ( tmp_171 + tmp_174 + tmp_179 ) ) +
                     tmp_203 * ( -tmp_189 * tmp_520 - tmp_202 * ( tmp_192 + tmp_195 + tmp_200 ) ) +
                     tmp_224 * ( -tmp_210 * tmp_520 - tmp_223 * ( tmp_213 + tmp_216 + tmp_221 ) ) +
                     tmp_245 * ( -tmp_231 * tmp_520 - tmp_244 * ( tmp_234 + tmp_237 + tmp_242 ) ) +
                     tmp_266 * ( -tmp_252 * tmp_520 - tmp_265 * ( tmp_255 + tmp_258 + tmp_263 ) ) +
                     tmp_287 * ( -tmp_273 * tmp_520 - tmp_286 * ( tmp_276 + tmp_279 + tmp_284 ) ) +
                     tmp_308 * ( -tmp_294 * tmp_520 - tmp_307 * ( tmp_297 + tmp_300 + tmp_305 ) ) +
                     tmp_329 * ( -tmp_315 * tmp_520 - tmp_328 * ( tmp_318 + tmp_321 + tmp_326 ) ) +
                     tmp_350 * ( -tmp_336 * tmp_520 - tmp_349 * ( tmp_339 + tmp_342 + tmp_347 ) ) +
                     tmp_371 * ( -tmp_357 * tmp_520 - tmp_370 * ( tmp_360 + tmp_363 + tmp_368 ) ) +
                     tmp_392 * ( -tmp_378 * tmp_520 - tmp_391 * ( tmp_381 + tmp_384 + tmp_389 ) ) +
                     tmp_413 * ( -tmp_399 * tmp_520 - tmp_412 * ( tmp_402 + tmp_405 + tmp_410 ) ) +
                     tmp_434 * ( -tmp_420 * tmp_520 - tmp_433 * ( tmp_423 + tmp_426 + tmp_431 ) ) +
                     tmp_455 * ( -tmp_441 * tmp_520 - tmp_454 * ( tmp_444 + tmp_447 + tmp_452 ) ) +
                     tmp_476 * ( -tmp_462 * tmp_520 - tmp_475 * ( tmp_465 + tmp_468 + tmp_473 ) ) +
                     tmp_497 * ( -tmp_483 * tmp_520 - tmp_496 * ( tmp_486 + tmp_489 + tmp_494 ) ) +
                     tmp_518 * ( -tmp_504 * tmp_520 - tmp_517 * ( tmp_507 + tmp_510 + tmp_515 ) ) +
                     tmp_98 * ( -tmp_46 * tmp_520 - tmp_96 * ( tmp_79 + tmp_83 + tmp_89 ) );
      real_t a_0_3 = tmp_119 * ( -tmp_105 * tmp_521 - tmp_118 * ( tmp_107 + tmp_110 + tmp_115 ) ) +
                     tmp_140 * ( -tmp_126 * tmp_521 - tmp_139 * ( tmp_128 + tmp_131 + tmp_136 ) ) +
                     tmp_161 * ( -tmp_147 * tmp_521 - tmp_160 * ( tmp_149 + tmp_152 + tmp_157 ) ) +
                     tmp_182 * ( -tmp_168 * tmp_521 - tmp_181 * ( tmp_170 + tmp_173 + tmp_178 ) ) +
                     tmp_203 * ( -tmp_189 * tmp_521 - tmp_202 * ( tmp_191 + tmp_194 + tmp_199 ) ) +
                     tmp_224 * ( -tmp_210 * tmp_521 - tmp_223 * ( tmp_212 + tmp_215 + tmp_220 ) ) +
                     tmp_245 * ( -tmp_231 * tmp_521 - tmp_244 * ( tmp_233 + tmp_236 + tmp_241 ) ) +
                     tmp_266 * ( -tmp_252 * tmp_521 - tmp_265 * ( tmp_254 + tmp_257 + tmp_262 ) ) +
                     tmp_287 * ( -tmp_273 * tmp_521 - tmp_286 * ( tmp_275 + tmp_278 + tmp_283 ) ) +
                     tmp_308 * ( -tmp_294 * tmp_521 - tmp_307 * ( tmp_296 + tmp_299 + tmp_304 ) ) +
                     tmp_329 * ( -tmp_315 * tmp_521 - tmp_328 * ( tmp_317 + tmp_320 + tmp_325 ) ) +
                     tmp_350 * ( -tmp_336 * tmp_521 - tmp_349 * ( tmp_338 + tmp_341 + tmp_346 ) ) +
                     tmp_371 * ( -tmp_357 * tmp_521 - tmp_370 * ( tmp_359 + tmp_362 + tmp_367 ) ) +
                     tmp_392 * ( -tmp_378 * tmp_521 - tmp_391 * ( tmp_380 + tmp_383 + tmp_388 ) ) +
                     tmp_413 * ( -tmp_399 * tmp_521 - tmp_412 * ( tmp_401 + tmp_404 + tmp_409 ) ) +
                     tmp_434 * ( -tmp_420 * tmp_521 - tmp_433 * ( tmp_422 + tmp_425 + tmp_430 ) ) +
                     tmp_455 * ( -tmp_441 * tmp_521 - tmp_454 * ( tmp_443 + tmp_446 + tmp_451 ) ) +
                     tmp_476 * ( -tmp_462 * tmp_521 - tmp_475 * ( tmp_464 + tmp_467 + tmp_472 ) ) +
                     tmp_497 * ( -tmp_483 * tmp_521 - tmp_496 * ( tmp_485 + tmp_488 + tmp_493 ) ) +
                     tmp_518 * ( -tmp_504 * tmp_521 - tmp_517 * ( tmp_506 + tmp_509 + tmp_514 ) ) +
                     tmp_98 * ( -tmp_46 * tmp_521 - tmp_96 * ( tmp_78 + tmp_82 + tmp_88 ) );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 0, 3 ) = a_0_3;
   }

   void integrateFacetDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                           const std::vector< Point3D >& coordsFacet,
                                           const Point3D&,
                                           const Point3D&     outwardNormal,
                                           const DGBasisInfo& trialBasis,
                                           const DGBasisInfo& testBasis,
                                           int                trialDegree,
                                           int                testDegree,
                                           MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_1_1 + tmp_0;
      real_t tmp_2  = -p_affine_0_2;
      real_t tmp_3  = p_affine_2_2 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_1 + tmp_0;
      real_t tmp_6  = p_affine_1_2 + tmp_2;
      real_t tmp_7  = tmp_5 * tmp_6;
      real_t tmp_8  = -p_affine_0_0;
      real_t tmp_9  = p_affine_1_0 + tmp_8;
      real_t tmp_10 = p_affine_3_2 + tmp_2;
      real_t tmp_11 = tmp_10 * tmp_5;
      real_t tmp_12 = p_affine_2_0 + tmp_8;
      real_t tmp_13 = p_affine_3_1 + tmp_0;
      real_t tmp_14 = tmp_13 * tmp_6;
      real_t tmp_15 = p_affine_3_0 + tmp_8;
      real_t tmp_16 = tmp_13 * tmp_3;
      real_t tmp_17 = tmp_1 * tmp_10;
      real_t tmp_18 =
          1.0 / ( tmp_11 * tmp_9 + tmp_12 * tmp_14 - tmp_12 * tmp_17 + tmp_15 * tmp_4 - tmp_15 * tmp_7 - tmp_16 * tmp_9 );
      real_t tmp_19 = tmp_18 * ( tmp_4 - tmp_7 );
      real_t tmp_20 = tmp_18 * ( tmp_14 - tmp_17 );
      real_t tmp_21 = tmp_18 * ( tmp_11 - tmp_16 );
      real_t tmp_22 = tmp_18 * ( tmp_12 * tmp_6 - tmp_3 * tmp_9 );
      real_t tmp_23 = tmp_18 * ( tmp_10 * tmp_9 - tmp_15 * tmp_6 );
      real_t tmp_24 = tmp_18 * ( -tmp_10 * tmp_12 + tmp_15 * tmp_3 );
      real_t tmp_25 = tmp_18 * ( -tmp_1 * tmp_12 + tmp_5 * tmp_9 );
      real_t tmp_26 = tmp_18 * ( tmp_1 * tmp_15 - tmp_13 * tmp_9 );
      real_t tmp_27 = tmp_18 * ( tmp_12 * tmp_13 - tmp_15 * tmp_5 );
      real_t tmp_28 = -p_affine_8_0;
      real_t tmp_29 = p_affine_10_0 + tmp_28;
      real_t tmp_30 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_31 = -p_affine_8_1;
      real_t tmp_32 = p_affine_10_1 + tmp_31;
      real_t tmp_33 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_34 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_35 = -p_affine_8_2;
      real_t tmp_36 = p_affine_10_2 + tmp_35;
      real_t tmp_37 =
          1.0 * std::pow( ( std::abs( tmp_29 * tmp_30 - tmp_32 * tmp_33 ) * std::abs( tmp_29 * tmp_30 - tmp_32 * tmp_33 ) ) +
                              ( std::abs( tmp_29 * tmp_34 - tmp_33 * tmp_36 ) * std::abs( tmp_29 * tmp_34 - tmp_33 * tmp_36 ) ) +
                              ( std::abs( tmp_30 * tmp_36 - tmp_32 * tmp_34 ) * std::abs( tmp_30 * tmp_36 - tmp_32 * tmp_34 ) ),
                          1.0 / 2.0 );
      real_t tmp_38 = tmp_37 * ( p_affine_13_0 * ( -tmp_19 - tmp_20 - tmp_21 ) + p_affine_13_1 * ( -tmp_22 - tmp_23 - tmp_24 ) +
                                 p_affine_13_2 * ( -tmp_25 - tmp_26 - tmp_27 ) );
      real_t tmp_39 = p_affine_9_2 + tmp_35;
      real_t tmp_40 = p_affine_8_2 + tmp_2;
      real_t tmp_41 = 0.93718850182767688 * tmp_36 + 0.031405749086161582 * tmp_39 + tmp_40;
      real_t tmp_42 = p_affine_9_1 + tmp_31;
      real_t tmp_43 = p_affine_8_1 + tmp_0;
      real_t tmp_44 = 0.93718850182767688 * tmp_32 + 0.031405749086161582 * tmp_42 + tmp_43;
      real_t tmp_45 = p_affine_9_0 + tmp_28;
      real_t tmp_46 = p_affine_8_0 + tmp_8;
      real_t tmp_47 = 0.93718850182767688 * tmp_29 + 0.031405749086161582 * tmp_45 + tmp_46;
      real_t tmp_48 = 0.0068572537431980923 * tmp_1 * ( tmp_21 * tmp_47 + tmp_24 * tmp_44 + tmp_27 * tmp_41 - 1.0 / 4.0 ) +
                      0.0068572537431980923 * tmp_13 * ( tmp_19 * tmp_47 + tmp_22 * tmp_44 + tmp_25 * tmp_41 - 1.0 / 4.0 ) +
                      0.0068572537431980923 * tmp_5 * ( tmp_20 * tmp_47 + tmp_23 * tmp_44 + tmp_26 * tmp_41 - 1.0 / 4.0 );
      real_t tmp_49 = 0.60796128279561268 * tmp_36 + 0.19601935860219369 * tmp_39 + tmp_40;
      real_t tmp_50 = 0.60796128279561268 * tmp_32 + 0.19601935860219369 * tmp_42 + tmp_43;
      real_t tmp_51 = 0.60796128279561268 * tmp_29 + 0.19601935860219369 * tmp_45 + tmp_46;
      real_t tmp_52 = 0.037198804536718075 * tmp_1 * ( tmp_21 * tmp_51 + tmp_24 * tmp_50 + tmp_27 * tmp_49 - 1.0 / 4.0 ) +
                      0.037198804536718075 * tmp_13 * ( tmp_19 * tmp_51 + tmp_22 * tmp_50 + tmp_25 * tmp_49 - 1.0 / 4.0 ) +
                      0.037198804536718075 * tmp_5 * ( tmp_20 * tmp_51 + tmp_23 * tmp_50 + tmp_26 * tmp_49 - 1.0 / 4.0 );
      real_t tmp_53 = 0.039308471900058539 * tmp_36 + 0.37605877282253791 * tmp_39 + tmp_40;
      real_t tmp_54 = 0.039308471900058539 * tmp_32 + 0.37605877282253791 * tmp_42 + tmp_43;
      real_t tmp_55 = 0.039308471900058539 * tmp_29 + 0.37605877282253791 * tmp_45 + tmp_46;
      real_t tmp_56 = 0.020848748529055869 * tmp_1 * ( tmp_21 * tmp_55 + tmp_24 * tmp_54 + tmp_27 * tmp_53 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_13 * ( tmp_19 * tmp_55 + tmp_22 * tmp_54 + tmp_25 * tmp_53 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_5 * ( tmp_20 * tmp_55 + tmp_23 * tmp_54 + tmp_26 * tmp_53 - 1.0 / 4.0 );
      real_t tmp_57 = 0.1711304259088916 * tmp_36 + 0.78764240869137092 * tmp_39 + tmp_40;
      real_t tmp_58 = 0.1711304259088916 * tmp_32 + 0.78764240869137092 * tmp_42 + tmp_43;
      real_t tmp_59 = 0.1711304259088916 * tmp_29 + 0.78764240869137092 * tmp_45 + tmp_46;
      real_t tmp_60 = 0.019202922745021479 * tmp_1 * ( tmp_21 * tmp_59 + tmp_24 * tmp_58 + tmp_27 * tmp_57 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_13 * ( tmp_19 * tmp_59 + tmp_22 * tmp_58 + tmp_25 * tmp_57 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_5 * ( tmp_20 * tmp_59 + tmp_23 * tmp_58 + tmp_26 * tmp_57 - 1.0 / 4.0 );
      real_t tmp_61 = 0.37605877282253791 * tmp_36 + 0.58463275527740355 * tmp_39 + tmp_40;
      real_t tmp_62 = 0.37605877282253791 * tmp_32 + 0.58463275527740355 * tmp_42 + tmp_43;
      real_t tmp_63 = 0.37605877282253791 * tmp_29 + 0.58463275527740355 * tmp_45 + tmp_46;
      real_t tmp_64 = 0.020848748529055869 * tmp_1 * ( tmp_21 * tmp_63 + tmp_24 * tmp_62 + tmp_27 * tmp_61 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_13 * ( tmp_19 * tmp_63 + tmp_22 * tmp_62 + tmp_25 * tmp_61 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_5 * ( tmp_20 * tmp_63 + tmp_23 * tmp_62 + tmp_26 * tmp_61 - 1.0 / 4.0 );
      real_t tmp_65 = 0.78764240869137092 * tmp_36 + 0.041227165399737475 * tmp_39 + tmp_40;
      real_t tmp_66 = 0.78764240869137092 * tmp_32 + 0.041227165399737475 * tmp_42 + tmp_43;
      real_t tmp_67 = 0.78764240869137092 * tmp_29 + 0.041227165399737475 * tmp_45 + tmp_46;
      real_t tmp_68 = 0.019202922745021479 * tmp_1 * ( tmp_21 * tmp_67 + tmp_24 * tmp_66 + tmp_27 * tmp_65 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_13 * ( tmp_19 * tmp_67 + tmp_22 * tmp_66 + tmp_25 * tmp_65 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_5 * ( tmp_20 * tmp_67 + tmp_23 * tmp_66 + tmp_26 * tmp_65 - 1.0 / 4.0 );
      real_t tmp_69 = 0.58463275527740355 * tmp_36 + 0.039308471900058539 * tmp_39 + tmp_40;
      real_t tmp_70 = 0.58463275527740355 * tmp_32 + 0.039308471900058539 * tmp_42 + tmp_43;
      real_t tmp_71 = 0.58463275527740355 * tmp_29 + 0.039308471900058539 * tmp_45 + tmp_46;
      real_t tmp_72 = 0.020848748529055869 * tmp_1 * ( tmp_21 * tmp_71 + tmp_24 * tmp_70 + tmp_27 * tmp_69 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_13 * ( tmp_19 * tmp_71 + tmp_22 * tmp_70 + tmp_25 * tmp_69 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_5 * ( tmp_20 * tmp_71 + tmp_23 * tmp_70 + tmp_26 * tmp_69 - 1.0 / 4.0 );
      real_t tmp_73 = 0.041227165399737475 * tmp_36 + 0.78764240869137092 * tmp_39 + tmp_40;
      real_t tmp_74 = 0.041227165399737475 * tmp_32 + 0.78764240869137092 * tmp_42 + tmp_43;
      real_t tmp_75 = 0.041227165399737475 * tmp_29 + 0.78764240869137092 * tmp_45 + tmp_46;
      real_t tmp_76 = 0.019202922745021479 * tmp_1 * ( tmp_21 * tmp_75 + tmp_24 * tmp_74 + tmp_27 * tmp_73 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_13 * ( tmp_19 * tmp_75 + tmp_22 * tmp_74 + tmp_25 * tmp_73 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_5 * ( tmp_20 * tmp_75 + tmp_23 * tmp_74 + tmp_26 * tmp_73 - 1.0 / 4.0 );
      real_t tmp_77 = 0.039308471900058539 * tmp_36 + 0.58463275527740355 * tmp_39 + tmp_40;
      real_t tmp_78 = 0.039308471900058539 * tmp_32 + 0.58463275527740355 * tmp_42 + tmp_43;
      real_t tmp_79 = 0.039308471900058539 * tmp_29 + 0.58463275527740355 * tmp_45 + tmp_46;
      real_t tmp_80 = 0.020848748529055869 * tmp_1 * ( tmp_21 * tmp_79 + tmp_24 * tmp_78 + tmp_27 * tmp_77 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_13 * ( tmp_19 * tmp_79 + tmp_22 * tmp_78 + tmp_25 * tmp_77 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_5 * ( tmp_20 * tmp_79 + tmp_23 * tmp_78 + tmp_26 * tmp_77 - 1.0 / 4.0 );
      real_t tmp_81 = 0.78764240869137092 * tmp_36 + 0.1711304259088916 * tmp_39 + tmp_40;
      real_t tmp_82 = 0.78764240869137092 * tmp_32 + 0.1711304259088916 * tmp_42 + tmp_43;
      real_t tmp_83 = 0.78764240869137092 * tmp_29 + 0.1711304259088916 * tmp_45 + tmp_46;
      real_t tmp_84 = 0.019202922745021479 * tmp_1 * ( tmp_21 * tmp_83 + tmp_24 * tmp_82 + tmp_27 * tmp_81 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_13 * ( tmp_19 * tmp_83 + tmp_22 * tmp_82 + tmp_25 * tmp_81 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_5 * ( tmp_20 * tmp_83 + tmp_23 * tmp_82 + tmp_26 * tmp_81 - 1.0 / 4.0 );
      real_t tmp_85 = 0.58463275527740355 * tmp_36 + 0.37605877282253791 * tmp_39 + tmp_40;
      real_t tmp_86 = 0.58463275527740355 * tmp_32 + 0.37605877282253791 * tmp_42 + tmp_43;
      real_t tmp_87 = 0.58463275527740355 * tmp_29 + 0.37605877282253791 * tmp_45 + tmp_46;
      real_t tmp_88 = 0.020848748529055869 * tmp_1 * ( tmp_21 * tmp_87 + tmp_24 * tmp_86 + tmp_27 * tmp_85 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_13 * ( tmp_19 * tmp_87 + tmp_22 * tmp_86 + tmp_25 * tmp_85 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_5 * ( tmp_20 * tmp_87 + tmp_23 * tmp_86 + tmp_26 * tmp_85 - 1.0 / 4.0 );
      real_t tmp_89 = 0.1711304259088916 * tmp_36 + 0.041227165399737475 * tmp_39 + tmp_40;
      real_t tmp_90 = 0.1711304259088916 * tmp_32 + 0.041227165399737475 * tmp_42 + tmp_43;
      real_t tmp_91 = 0.1711304259088916 * tmp_29 + 0.041227165399737475 * tmp_45 + tmp_46;
      real_t tmp_92 = 0.019202922745021479 * tmp_1 * ( tmp_21 * tmp_91 + tmp_24 * tmp_90 + tmp_27 * tmp_89 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_13 * ( tmp_19 * tmp_91 + tmp_22 * tmp_90 + tmp_25 * tmp_89 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_5 * ( tmp_20 * tmp_91 + tmp_23 * tmp_90 + tmp_26 * tmp_89 - 1.0 / 4.0 );
      real_t tmp_93 = 0.19107600050469298 * tmp_36 + 0.40446199974765351 * tmp_39 + tmp_40;
      real_t tmp_94 = 0.19107600050469298 * tmp_32 + 0.40446199974765351 * tmp_42 + tmp_43;
      real_t tmp_95 = 0.19107600050469298 * tmp_29 + 0.40446199974765351 * tmp_45 + tmp_46;
      real_t tmp_96 = 0.042507265838595799 * tmp_1 * ( tmp_21 * tmp_95 + tmp_24 * tmp_94 + tmp_27 * tmp_93 - 1.0 / 4.0 ) +
                      0.042507265838595799 * tmp_13 * ( tmp_19 * tmp_95 + tmp_22 * tmp_94 + tmp_25 * tmp_93 - 1.0 / 4.0 ) +
                      0.042507265838595799 * tmp_5 * ( tmp_20 * tmp_95 + tmp_23 * tmp_94 + tmp_26 * tmp_93 - 1.0 / 4.0 );
      real_t tmp_97  = 0.37605877282253791 * tmp_36 + 0.039308471900058539 * tmp_39 + tmp_40;
      real_t tmp_98  = 0.37605877282253791 * tmp_32 + 0.039308471900058539 * tmp_42 + tmp_43;
      real_t tmp_99  = 0.37605877282253791 * tmp_29 + 0.039308471900058539 * tmp_45 + tmp_46;
      real_t tmp_100 = 0.020848748529055869 * tmp_1 * ( tmp_21 * tmp_99 + tmp_24 * tmp_98 + tmp_27 * tmp_97 - 1.0 / 4.0 ) +
                       0.020848748529055869 * tmp_13 * ( tmp_19 * tmp_99 + tmp_22 * tmp_98 + tmp_25 * tmp_97 - 1.0 / 4.0 ) +
                       0.020848748529055869 * tmp_5 * ( tmp_20 * tmp_99 + tmp_23 * tmp_98 + tmp_26 * tmp_97 - 1.0 / 4.0 );
      real_t tmp_101 = 0.031405749086161582 * tmp_36 + 0.93718850182767688 * tmp_39 + tmp_40;
      real_t tmp_102 = 0.031405749086161582 * tmp_32 + 0.93718850182767688 * tmp_42 + tmp_43;
      real_t tmp_103 = 0.031405749086161582 * tmp_29 + 0.93718850182767688 * tmp_45 + tmp_46;
      real_t tmp_104 = 0.0068572537431980923 * tmp_1 * ( tmp_101 * tmp_27 + tmp_102 * tmp_24 + tmp_103 * tmp_21 - 1.0 / 4.0 ) +
                       0.0068572537431980923 * tmp_13 * ( tmp_101 * tmp_25 + tmp_102 * tmp_22 + tmp_103 * tmp_19 - 1.0 / 4.0 ) +
                       0.0068572537431980923 * tmp_5 * ( tmp_101 * tmp_26 + tmp_102 * tmp_23 + tmp_103 * tmp_20 - 1.0 / 4.0 );
      real_t tmp_105 = 0.19601935860219369 * tmp_36 + 0.60796128279561268 * tmp_39 + tmp_40;
      real_t tmp_106 = 0.19601935860219369 * tmp_32 + 0.60796128279561268 * tmp_42 + tmp_43;
      real_t tmp_107 = 0.19601935860219369 * tmp_29 + 0.60796128279561268 * tmp_45 + tmp_46;
      real_t tmp_108 = 0.037198804536718075 * tmp_1 * ( tmp_105 * tmp_27 + tmp_106 * tmp_24 + tmp_107 * tmp_21 - 1.0 / 4.0 ) +
                       0.037198804536718075 * tmp_13 * ( tmp_105 * tmp_25 + tmp_106 * tmp_22 + tmp_107 * tmp_19 - 1.0 / 4.0 ) +
                       0.037198804536718075 * tmp_5 * ( tmp_105 * tmp_26 + tmp_106 * tmp_23 + tmp_107 * tmp_20 - 1.0 / 4.0 );
      real_t tmp_109 = 0.40446199974765351 * tmp_36 + 0.19107600050469298 * tmp_39 + tmp_40;
      real_t tmp_110 = 0.40446199974765351 * tmp_32 + 0.19107600050469298 * tmp_42 + tmp_43;
      real_t tmp_111 = 0.40446199974765351 * tmp_29 + 0.19107600050469298 * tmp_45 + tmp_46;
      real_t tmp_112 = 0.042507265838595799 * tmp_1 * ( tmp_109 * tmp_27 + tmp_110 * tmp_24 + tmp_111 * tmp_21 - 1.0 / 4.0 ) +
                       0.042507265838595799 * tmp_13 * ( tmp_109 * tmp_25 + tmp_110 * tmp_22 + tmp_111 * tmp_19 - 1.0 / 4.0 ) +
                       0.042507265838595799 * tmp_5 * ( tmp_109 * tmp_26 + tmp_110 * tmp_23 + tmp_111 * tmp_20 - 1.0 / 4.0 );
      real_t tmp_113 = 0.031405749086161582 * tmp_36 + 0.031405749086161582 * tmp_39 + tmp_40;
      real_t tmp_114 = 0.031405749086161582 * tmp_32 + 0.031405749086161582 * tmp_42 + tmp_43;
      real_t tmp_115 = 0.031405749086161582 * tmp_29 + 0.031405749086161582 * tmp_45 + tmp_46;
      real_t tmp_116 = 0.0068572537431980923 * tmp_1 * ( tmp_113 * tmp_27 + tmp_114 * tmp_24 + tmp_115 * tmp_21 - 1.0 / 4.0 ) +
                       0.0068572537431980923 * tmp_13 * ( tmp_113 * tmp_25 + tmp_114 * tmp_22 + tmp_115 * tmp_19 - 1.0 / 4.0 ) +
                       0.0068572537431980923 * tmp_5 * ( tmp_113 * tmp_26 + tmp_114 * tmp_23 + tmp_115 * tmp_20 - 1.0 / 4.0 );
      real_t tmp_117 = 0.19601935860219369 * tmp_36 + 0.19601935860219369 * tmp_39 + tmp_40;
      real_t tmp_118 = 0.19601935860219369 * tmp_32 + 0.19601935860219369 * tmp_42 + tmp_43;
      real_t tmp_119 = 0.19601935860219369 * tmp_29 + 0.19601935860219369 * tmp_45 + tmp_46;
      real_t tmp_120 = 0.037198804536718075 * tmp_1 * ( tmp_117 * tmp_27 + tmp_118 * tmp_24 + tmp_119 * tmp_21 - 1.0 / 4.0 ) +
                       0.037198804536718075 * tmp_13 * ( tmp_117 * tmp_25 + tmp_118 * tmp_22 + tmp_119 * tmp_19 - 1.0 / 4.0 ) +
                       0.037198804536718075 * tmp_5 * ( tmp_117 * tmp_26 + tmp_118 * tmp_23 + tmp_119 * tmp_20 - 1.0 / 4.0 );
      real_t tmp_121 = 0.40446199974765351 * tmp_36 + 0.40446199974765351 * tmp_39 + tmp_40;
      real_t tmp_122 = 0.40446199974765351 * tmp_32 + 0.40446199974765351 * tmp_42 + tmp_43;
      real_t tmp_123 = 0.40446199974765351 * tmp_29 + 0.40446199974765351 * tmp_45 + tmp_46;
      real_t tmp_124 = 0.042507265838595799 * tmp_1 * ( tmp_121 * tmp_27 + tmp_122 * tmp_24 + tmp_123 * tmp_21 - 1.0 / 4.0 ) +
                       0.042507265838595799 * tmp_13 * ( tmp_121 * tmp_25 + tmp_122 * tmp_22 + tmp_123 * tmp_19 - 1.0 / 4.0 ) +
                       0.042507265838595799 * tmp_5 * ( tmp_121 * tmp_26 + tmp_122 * tmp_23 + tmp_123 * tmp_20 - 1.0 / 4.0 );
      real_t tmp_125 = 0.041227165399737475 * tmp_36 + 0.1711304259088916 * tmp_39 + tmp_40;
      real_t tmp_126 = 0.041227165399737475 * tmp_32 + 0.1711304259088916 * tmp_42 + tmp_43;
      real_t tmp_127 = 0.041227165399737475 * tmp_29 + 0.1711304259088916 * tmp_45 + tmp_46;
      real_t tmp_128 = 0.019202922745021479 * tmp_1 * ( tmp_125 * tmp_27 + tmp_126 * tmp_24 + tmp_127 * tmp_21 - 1.0 / 4.0 ) +
                       0.019202922745021479 * tmp_13 * ( tmp_125 * tmp_25 + tmp_126 * tmp_22 + tmp_127 * tmp_19 - 1.0 / 4.0 ) +
                       0.019202922745021479 * tmp_5 * ( tmp_125 * tmp_26 + tmp_126 * tmp_23 + tmp_127 * tmp_20 - 1.0 / 4.0 );
      real_t tmp_129 = tmp_37 * ( p_affine_13_0 * tmp_21 + p_affine_13_1 * tmp_24 + p_affine_13_2 * tmp_27 );
      real_t tmp_130 = tmp_37 * ( p_affine_13_0 * tmp_20 + p_affine_13_1 * tmp_23 + p_affine_13_2 * tmp_26 );
      real_t tmp_131 = tmp_37 * ( p_affine_13_0 * tmp_19 + p_affine_13_1 * tmp_22 + p_affine_13_2 * tmp_25 );
      real_t a_0_0   = -tmp_100 * tmp_38 - tmp_104 * tmp_38 - tmp_108 * tmp_38 - tmp_112 * tmp_38 - tmp_116 * tmp_38 -
                     tmp_120 * tmp_38 - tmp_124 * tmp_38 - tmp_128 * tmp_38 - tmp_38 * tmp_48 - tmp_38 * tmp_52 -
                     tmp_38 * tmp_56 - tmp_38 * tmp_60 - tmp_38 * tmp_64 - tmp_38 * tmp_68 - tmp_38 * tmp_72 - tmp_38 * tmp_76 -
                     tmp_38 * tmp_80 - tmp_38 * tmp_84 - tmp_38 * tmp_88 - tmp_38 * tmp_92 - tmp_38 * tmp_96;
      real_t a_0_1 = -tmp_100 * tmp_129 - tmp_104 * tmp_129 - tmp_108 * tmp_129 - tmp_112 * tmp_129 - tmp_116 * tmp_129 -
                     tmp_120 * tmp_129 - tmp_124 * tmp_129 - tmp_128 * tmp_129 - tmp_129 * tmp_48 - tmp_129 * tmp_52 -
                     tmp_129 * tmp_56 - tmp_129 * tmp_60 - tmp_129 * tmp_64 - tmp_129 * tmp_68 - tmp_129 * tmp_72 -
                     tmp_129 * tmp_76 - tmp_129 * tmp_80 - tmp_129 * tmp_84 - tmp_129 * tmp_88 - tmp_129 * tmp_92 -
                     tmp_129 * tmp_96;
      real_t a_0_2 = -tmp_100 * tmp_130 - tmp_104 * tmp_130 - tmp_108 * tmp_130 - tmp_112 * tmp_130 - tmp_116 * tmp_130 -
                     tmp_120 * tmp_130 - tmp_124 * tmp_130 - tmp_128 * tmp_130 - tmp_130 * tmp_48 - tmp_130 * tmp_52 -
                     tmp_130 * tmp_56 - tmp_130 * tmp_60 - tmp_130 * tmp_64 - tmp_130 * tmp_68 - tmp_130 * tmp_72 -
                     tmp_130 * tmp_76 - tmp_130 * tmp_80 - tmp_130 * tmp_84 - tmp_130 * tmp_88 - tmp_130 * tmp_92 -
                     tmp_130 * tmp_96;
      real_t a_0_3 = -tmp_100 * tmp_131 - tmp_104 * tmp_131 - tmp_108 * tmp_131 - tmp_112 * tmp_131 - tmp_116 * tmp_131 -
                     tmp_120 * tmp_131 - tmp_124 * tmp_131 - tmp_128 * tmp_131 - tmp_131 * tmp_48 - tmp_131 * tmp_52 -
                     tmp_131 * tmp_56 - tmp_131 * tmp_60 - tmp_131 * tmp_64 - tmp_131 * tmp_68 - tmp_131 * tmp_72 -
                     tmp_131 * tmp_76 - tmp_131 * tmp_80 - tmp_131 * tmp_84 - tmp_131 * tmp_88 - tmp_131 * tmp_92 -
                     tmp_131 * tmp_96;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 0, 3 ) = a_0_3;
   }
};

class EGIIPGVectorLaplaceFormEP1_2 : public hyteg::dg::DGForm
{
 protected:
   void integrateVolume2D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override
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
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >& coordsElement,
                                       const std::vector< Point3D >& coordsFacet,
                                       const Point3D&                oppositeVertex,
                                       const Point3D&                outwardNormal,
                                       const DGBasisInfo&            trialBasis,
                                       const DGBasisInfo&            testBasis,
                                       int                           trialDegree,
                                       int                           testDegree,
                                       MatrixXr&                     elMat ) const override
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
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >& coordsElementInner,
                                          const std::vector< Point3D >& coordsElementOuter,
                                          const std::vector< Point3D >& coordsFacet,
                                          const Point3D&                oppositeVertexInnerElement,
                                          const Point3D&                oppositeVertexOuterElement,
                                          const Point3D&                outwardNormal,
                                          const DGBasisInfo&            trialBasis,
                                          const DGBasisInfo&            testBasis,
                                          int                           trialDegree,
                                          int                           testDegree,
                                          MatrixXr&                     elMat ) const override
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
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                                   const std::vector< Point3D >& coordsFacet,
                                                   const Point3D&                oppositeVertex,
                                                   const Point3D&                outwardNormal,
                                                   const DGBasisInfo&            trialBasis,
                                                   const DGBasisInfo&            testBasis,
                                                   int                           trialDegree,
                                                   int                           testDegree,
                                                   MatrixXr&                     elMat ) const override
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
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
   }

   void integrateRHSDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                         const std::vector< Point3D >& coordsFacet,
                                         const Point3D&                oppositeVertex,
                                         const Point3D&                outwardNormal,
                                         const DGBasisInfo&            basis,
                                         int                           degree,
                                         MatrixXr&                     elMat ) const override
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
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                         const std::vector< Point3D >& coordsFacet,
                                         const Point3D&                oppositeVertex,
                                         const Point3D&                outwardNormal,
                                         const DGBasisInfo&            basis,
                                         int                           degree,
                                         MatrixXr&                     elMat ) const override
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
   void integrateVolume3D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_0;
      real_t tmp_6  = p_affine_1_1 + tmp_2;
      real_t tmp_7  = tmp_5 * tmp_6;
      real_t tmp_8  = -p_affine_0_2;
      real_t tmp_9  = p_affine_3_2 + tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_8;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11 * tmp_5;
      real_t tmp_13 = p_affine_2_2 + tmp_8;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = tmp_14 * tmp_6;
      real_t tmp_16 = tmp_1 * tmp_11;
      real_t tmp_17 = tmp_14 * tmp_3;
      real_t tmp_18 =
          1.0 / ( tmp_10 * tmp_12 - tmp_10 * tmp_17 + tmp_13 * tmp_15 - tmp_13 * tmp_16 + tmp_4 * tmp_9 - tmp_7 * tmp_9 );
      real_t tmp_19 = tmp_18 * ( tmp_4 - tmp_7 );
      real_t tmp_20 = tmp_18 * ( tmp_15 - tmp_16 );
      real_t tmp_21 = tmp_18 * ( tmp_12 - tmp_17 );
      real_t tmp_22 = tmp_10 * tmp_21 + tmp_13 * tmp_20 + tmp_19 * tmp_9;
      real_t tmp_23 = tmp_18 * ( -tmp_1 * tmp_13 + tmp_10 * tmp_5 );
      real_t tmp_24 = tmp_18 * ( tmp_1 * tmp_9 - tmp_10 * tmp_14 );
      real_t tmp_25 = tmp_18 * ( tmp_13 * tmp_14 - tmp_5 * tmp_9 );
      real_t tmp_26 = tmp_10 * tmp_25 + tmp_13 * tmp_24 + tmp_23 * tmp_9;
      real_t tmp_27 = tmp_18 * ( -tmp_10 * tmp_3 + tmp_13 * tmp_6 );
      real_t tmp_28 = tmp_18 * ( tmp_10 * tmp_11 - tmp_6 * tmp_9 );
      real_t tmp_29 = tmp_18 * ( -tmp_11 * tmp_13 + tmp_3 * tmp_9 );
      real_t tmp_30 = tmp_10 * tmp_29 + tmp_13 * tmp_28 + tmp_27 * tmp_9;
      real_t tmp_31 = p_affine_0_0 * p_affine_1_1;
      real_t tmp_32 = p_affine_0_0 * p_affine_1_2;
      real_t tmp_33 = p_affine_2_1 * p_affine_3_2;
      real_t tmp_34 = p_affine_0_1 * p_affine_1_0;
      real_t tmp_35 = p_affine_0_1 * p_affine_1_2;
      real_t tmp_36 = p_affine_2_2 * p_affine_3_0;
      real_t tmp_37 = p_affine_0_2 * p_affine_1_0;
      real_t tmp_38 = p_affine_0_2 * p_affine_1_1;
      real_t tmp_39 = p_affine_2_0 * p_affine_3_1;
      real_t tmp_40 = p_affine_2_2 * p_affine_3_1;
      real_t tmp_41 = p_affine_2_0 * p_affine_3_2;
      real_t tmp_42 = p_affine_2_1 * p_affine_3_0;
      real_t tmp_43 = std::abs( p_affine_0_0 * tmp_33 - p_affine_0_0 * tmp_40 + p_affine_0_1 * tmp_36 - p_affine_0_1 * tmp_41 +
                                p_affine_0_2 * tmp_39 - p_affine_0_2 * tmp_42 - p_affine_1_0 * tmp_33 + p_affine_1_0 * tmp_40 -
                                p_affine_1_1 * tmp_36 + p_affine_1_1 * tmp_41 - p_affine_1_2 * tmp_39 + p_affine_1_2 * tmp_42 +
                                p_affine_2_0 * tmp_35 - p_affine_2_0 * tmp_38 - p_affine_2_1 * tmp_32 + p_affine_2_1 * tmp_37 +
                                p_affine_2_2 * tmp_31 - p_affine_2_2 * tmp_34 - p_affine_3_0 * tmp_35 + p_affine_3_0 * tmp_38 +
                                p_affine_3_1 * tmp_32 - p_affine_3_1 * tmp_37 - p_affine_3_2 * tmp_31 + p_affine_3_2 * tmp_34 );
      real_t tmp_44 = tmp_43 * ( tmp_22 * ( -tmp_19 - tmp_20 - tmp_21 ) + tmp_26 * ( -tmp_23 - tmp_24 - tmp_25 ) +
                                 tmp_30 * ( -tmp_27 - tmp_28 - tmp_29 ) );
      real_t tmp_45 = tmp_43 * ( tmp_21 * tmp_22 + tmp_25 * tmp_26 + tmp_29 * tmp_30 );
      real_t tmp_46 = tmp_43 * ( tmp_20 * tmp_22 + tmp_24 * tmp_26 + tmp_28 * tmp_30 );
      real_t tmp_47 = tmp_43 * ( tmp_19 * tmp_22 + tmp_23 * tmp_26 + tmp_27 * tmp_30 );
      real_t a_0_0  = 0.1666666666666668 * tmp_44;
      real_t a_0_1  = 0.1666666666666668 * tmp_45;
      real_t a_0_2  = 0.1666666666666668 * tmp_46;
      real_t a_0_3  = 0.1666666666666668 * tmp_47;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 0, 3 ) = a_0_3;
   }

   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                               const std::vector< Point3D >& coordsFacet,
                               const Point3D&,
                               const Point3D&     outwardNormal,
                               const DGBasisInfo& trialBasis,
                               const DGBasisInfo& testBasis,
                               int                trialDegree,
                               int                testDegree,
                               MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_2;
      real_t tmp_1  = p_affine_1_2 + tmp_0;
      real_t tmp_2  = -p_affine_8_2;
      real_t tmp_3  = p_affine_9_2 + tmp_2;
      real_t tmp_4  = p_affine_10_2 + tmp_2;
      real_t tmp_5  = p_affine_8_2 + tmp_0;
      real_t tmp_6  = 0.031405749086161582 * tmp_3 + 0.93718850182767688 * tmp_4 + tmp_5;
      real_t tmp_7  = -p_affine_0_0;
      real_t tmp_8  = p_affine_2_0 + tmp_7;
      real_t tmp_9  = -p_affine_0_1;
      real_t tmp_10 = p_affine_3_1 + tmp_9;
      real_t tmp_11 = p_affine_3_0 + tmp_7;
      real_t tmp_12 = p_affine_2_1 + tmp_9;
      real_t tmp_13 = p_affine_1_0 + tmp_7;
      real_t tmp_14 = p_affine_3_2 + tmp_0;
      real_t tmp_15 = tmp_12 * tmp_14;
      real_t tmp_16 = tmp_1 * tmp_10;
      real_t tmp_17 = p_affine_1_1 + tmp_9;
      real_t tmp_18 = p_affine_2_2 + tmp_0;
      real_t tmp_19 = tmp_17 * tmp_18;
      real_t tmp_20 = tmp_10 * tmp_18;
      real_t tmp_21 = tmp_14 * tmp_17;
      real_t tmp_22 = tmp_1 * tmp_12;
      real_t tmp_23 =
          1.0 / ( tmp_11 * tmp_19 - tmp_11 * tmp_22 + tmp_13 * tmp_15 - tmp_13 * tmp_20 + tmp_16 * tmp_8 - tmp_21 * tmp_8 );
      real_t tmp_24 = tmp_23 * ( tmp_10 * tmp_8 - tmp_11 * tmp_12 );
      real_t tmp_25 = tmp_24 * tmp_6;
      real_t tmp_26 = -p_affine_8_1;
      real_t tmp_27 = p_affine_9_1 + tmp_26;
      real_t tmp_28 = p_affine_10_1 + tmp_26;
      real_t tmp_29 = p_affine_8_1 + tmp_9;
      real_t tmp_30 = 0.031405749086161582 * tmp_27 + 0.93718850182767688 * tmp_28 + tmp_29;
      real_t tmp_31 = tmp_23 * ( tmp_11 * tmp_18 - tmp_14 * tmp_8 );
      real_t tmp_32 = tmp_30 * tmp_31;
      real_t tmp_33 = -p_affine_8_0;
      real_t tmp_34 = p_affine_9_0 + tmp_33;
      real_t tmp_35 = p_affine_10_0 + tmp_33;
      real_t tmp_36 = p_affine_8_0 + tmp_7;
      real_t tmp_37 = 0.031405749086161582 * tmp_34 + 0.93718850182767688 * tmp_35 + tmp_36;
      real_t tmp_38 = tmp_23 * ( tmp_15 - tmp_20 );
      real_t tmp_39 = tmp_37 * tmp_38;
      real_t tmp_40 = tmp_25 + tmp_32 + tmp_39;
      real_t tmp_41 = tmp_23 * ( -tmp_10 * tmp_13 + tmp_11 * tmp_17 );
      real_t tmp_42 = tmp_41 * tmp_6;
      real_t tmp_43 = tmp_23 * ( -tmp_1 * tmp_11 + tmp_13 * tmp_14 );
      real_t tmp_44 = tmp_30 * tmp_43;
      real_t tmp_45 = tmp_23 * ( tmp_16 - tmp_21 );
      real_t tmp_46 = tmp_37 * tmp_45;
      real_t tmp_47 = tmp_42 + tmp_44 + tmp_46;
      real_t tmp_48 = tmp_23 * ( tmp_12 * tmp_13 - tmp_17 * tmp_8 );
      real_t tmp_49 = tmp_48 * tmp_6;
      real_t tmp_50 = tmp_23 * ( tmp_1 * tmp_8 - tmp_13 * tmp_18 );
      real_t tmp_51 = tmp_30 * tmp_50;
      real_t tmp_52 = tmp_23 * ( tmp_19 - tmp_22 );
      real_t tmp_53 = tmp_37 * tmp_52;
      real_t tmp_54 = tmp_49 + tmp_51 + tmp_53;
      real_t tmp_55 = tmp_1 * ( tmp_40 - 1.0 / 4.0 ) + tmp_14 * ( tmp_54 - 1.0 / 4.0 ) + tmp_18 * ( tmp_47 - 1.0 / 4.0 );
      real_t tmp_56 = 0.5 * p_affine_13_0 * ( -tmp_38 - tmp_45 - tmp_52 ) + 0.5 * p_affine_13_1 * ( -tmp_31 - tmp_43 - tmp_50 ) +
                      0.5 * p_affine_13_2 * ( -tmp_24 - tmp_41 - tmp_48 );
      real_t tmp_57 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_58 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_59 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_60 = ( std::abs( tmp_28 * tmp_58 - tmp_35 * tmp_57 ) * std::abs( tmp_28 * tmp_58 - tmp_35 * tmp_57 ) ) +
                      ( std::abs( tmp_28 * tmp_59 - tmp_4 * tmp_57 ) * std::abs( tmp_28 * tmp_59 - tmp_4 * tmp_57 ) ) +
                      ( std::abs( tmp_35 * tmp_59 - tmp_4 * tmp_58 ) * std::abs( tmp_35 * tmp_59 - tmp_4 * tmp_58 ) );
      real_t tmp_61  = 1.0 * std::pow( tmp_60, -0.25 );
      real_t tmp_62  = tmp_55 * tmp_61;
      real_t tmp_63  = 1.0 * std::pow( tmp_60, 1.0 / 2.0 );
      real_t tmp_64  = 0.0068572537431980923 * tmp_63;
      real_t tmp_65  = 0.19601935860219369 * tmp_3 + 0.60796128279561268 * tmp_4 + tmp_5;
      real_t tmp_66  = tmp_24 * tmp_65;
      real_t tmp_67  = 0.19601935860219369 * tmp_27 + 0.60796128279561268 * tmp_28 + tmp_29;
      real_t tmp_68  = tmp_31 * tmp_67;
      real_t tmp_69  = 0.19601935860219369 * tmp_34 + 0.60796128279561268 * tmp_35 + tmp_36;
      real_t tmp_70  = tmp_38 * tmp_69;
      real_t tmp_71  = tmp_66 + tmp_68 + tmp_70;
      real_t tmp_72  = tmp_41 * tmp_65;
      real_t tmp_73  = tmp_43 * tmp_67;
      real_t tmp_74  = tmp_45 * tmp_69;
      real_t tmp_75  = tmp_72 + tmp_73 + tmp_74;
      real_t tmp_76  = tmp_48 * tmp_65;
      real_t tmp_77  = tmp_50 * tmp_67;
      real_t tmp_78  = tmp_52 * tmp_69;
      real_t tmp_79  = tmp_76 + tmp_77 + tmp_78;
      real_t tmp_80  = tmp_1 * ( tmp_71 - 1.0 / 4.0 ) + tmp_14 * ( tmp_79 - 1.0 / 4.0 ) + tmp_18 * ( tmp_75 - 1.0 / 4.0 );
      real_t tmp_81  = tmp_61 * tmp_80;
      real_t tmp_82  = 0.037198804536718075 * tmp_63;
      real_t tmp_83  = 0.37605877282253791 * tmp_3 + 0.039308471900058539 * tmp_4 + tmp_5;
      real_t tmp_84  = tmp_24 * tmp_83;
      real_t tmp_85  = 0.37605877282253791 * tmp_27 + 0.039308471900058539 * tmp_28 + tmp_29;
      real_t tmp_86  = tmp_31 * tmp_85;
      real_t tmp_87  = 0.37605877282253791 * tmp_34 + 0.039308471900058539 * tmp_35 + tmp_36;
      real_t tmp_88  = tmp_38 * tmp_87;
      real_t tmp_89  = tmp_84 + tmp_86 + tmp_88;
      real_t tmp_90  = tmp_41 * tmp_83;
      real_t tmp_91  = tmp_43 * tmp_85;
      real_t tmp_92  = tmp_45 * tmp_87;
      real_t tmp_93  = tmp_90 + tmp_91 + tmp_92;
      real_t tmp_94  = tmp_48 * tmp_83;
      real_t tmp_95  = tmp_50 * tmp_85;
      real_t tmp_96  = tmp_52 * tmp_87;
      real_t tmp_97  = tmp_94 + tmp_95 + tmp_96;
      real_t tmp_98  = tmp_1 * ( tmp_89 - 1.0 / 4.0 ) + tmp_14 * ( tmp_97 - 1.0 / 4.0 ) + tmp_18 * ( tmp_93 - 1.0 / 4.0 );
      real_t tmp_99  = tmp_61 * tmp_98;
      real_t tmp_100 = 0.020848748529055869 * tmp_63;
      real_t tmp_101 = 0.78764240869137092 * tmp_3 + 0.1711304259088916 * tmp_4 + tmp_5;
      real_t tmp_102 = tmp_101 * tmp_24;
      real_t tmp_103 = 0.78764240869137092 * tmp_27 + 0.1711304259088916 * tmp_28 + tmp_29;
      real_t tmp_104 = tmp_103 * tmp_31;
      real_t tmp_105 = 0.78764240869137092 * tmp_34 + 0.1711304259088916 * tmp_35 + tmp_36;
      real_t tmp_106 = tmp_105 * tmp_38;
      real_t tmp_107 = tmp_102 + tmp_104 + tmp_106;
      real_t tmp_108 = tmp_101 * tmp_41;
      real_t tmp_109 = tmp_103 * tmp_43;
      real_t tmp_110 = tmp_105 * tmp_45;
      real_t tmp_111 = tmp_108 + tmp_109 + tmp_110;
      real_t tmp_112 = tmp_101 * tmp_48;
      real_t tmp_113 = tmp_103 * tmp_50;
      real_t tmp_114 = tmp_105 * tmp_52;
      real_t tmp_115 = tmp_112 + tmp_113 + tmp_114;
      real_t tmp_116 = tmp_1 * ( tmp_107 - 1.0 / 4.0 ) + tmp_14 * ( tmp_115 - 1.0 / 4.0 ) + tmp_18 * ( tmp_111 - 1.0 / 4.0 );
      real_t tmp_117 = tmp_116 * tmp_61;
      real_t tmp_118 = 0.019202922745021479 * tmp_63;
      real_t tmp_119 = 0.58463275527740355 * tmp_3 + 0.37605877282253791 * tmp_4 + tmp_5;
      real_t tmp_120 = tmp_119 * tmp_24;
      real_t tmp_121 = 0.58463275527740355 * tmp_27 + 0.37605877282253791 * tmp_28 + tmp_29;
      real_t tmp_122 = tmp_121 * tmp_31;
      real_t tmp_123 = 0.58463275527740355 * tmp_34 + 0.37605877282253791 * tmp_35 + tmp_36;
      real_t tmp_124 = tmp_123 * tmp_38;
      real_t tmp_125 = tmp_120 + tmp_122 + tmp_124;
      real_t tmp_126 = tmp_119 * tmp_41;
      real_t tmp_127 = tmp_121 * tmp_43;
      real_t tmp_128 = tmp_123 * tmp_45;
      real_t tmp_129 = tmp_126 + tmp_127 + tmp_128;
      real_t tmp_130 = tmp_119 * tmp_48;
      real_t tmp_131 = tmp_121 * tmp_50;
      real_t tmp_132 = tmp_123 * tmp_52;
      real_t tmp_133 = tmp_130 + tmp_131 + tmp_132;
      real_t tmp_134 = tmp_1 * ( tmp_125 - 1.0 / 4.0 ) + tmp_14 * ( tmp_133 - 1.0 / 4.0 ) + tmp_18 * ( tmp_129 - 1.0 / 4.0 );
      real_t tmp_135 = tmp_134 * tmp_61;
      real_t tmp_136 = 0.020848748529055869 * tmp_63;
      real_t tmp_137 = 0.041227165399737475 * tmp_3 + 0.78764240869137092 * tmp_4 + tmp_5;
      real_t tmp_138 = tmp_137 * tmp_24;
      real_t tmp_139 = 0.041227165399737475 * tmp_27 + 0.78764240869137092 * tmp_28 + tmp_29;
      real_t tmp_140 = tmp_139 * tmp_31;
      real_t tmp_141 = 0.041227165399737475 * tmp_34 + 0.78764240869137092 * tmp_35 + tmp_36;
      real_t tmp_142 = tmp_141 * tmp_38;
      real_t tmp_143 = tmp_138 + tmp_140 + tmp_142;
      real_t tmp_144 = tmp_137 * tmp_41;
      real_t tmp_145 = tmp_139 * tmp_43;
      real_t tmp_146 = tmp_141 * tmp_45;
      real_t tmp_147 = tmp_144 + tmp_145 + tmp_146;
      real_t tmp_148 = tmp_137 * tmp_48;
      real_t tmp_149 = tmp_139 * tmp_50;
      real_t tmp_150 = tmp_141 * tmp_52;
      real_t tmp_151 = tmp_148 + tmp_149 + tmp_150;
      real_t tmp_152 = tmp_1 * ( tmp_143 - 1.0 / 4.0 ) + tmp_14 * ( tmp_151 - 1.0 / 4.0 ) + tmp_18 * ( tmp_147 - 1.0 / 4.0 );
      real_t tmp_153 = tmp_152 * tmp_61;
      real_t tmp_154 = 0.019202922745021479 * tmp_63;
      real_t tmp_155 = 0.039308471900058539 * tmp_3 + 0.58463275527740355 * tmp_4 + tmp_5;
      real_t tmp_156 = tmp_155 * tmp_24;
      real_t tmp_157 = 0.039308471900058539 * tmp_27 + 0.58463275527740355 * tmp_28 + tmp_29;
      real_t tmp_158 = tmp_157 * tmp_31;
      real_t tmp_159 = 0.039308471900058539 * tmp_34 + 0.58463275527740355 * tmp_35 + tmp_36;
      real_t tmp_160 = tmp_159 * tmp_38;
      real_t tmp_161 = tmp_156 + tmp_158 + tmp_160;
      real_t tmp_162 = tmp_155 * tmp_41;
      real_t tmp_163 = tmp_157 * tmp_43;
      real_t tmp_164 = tmp_159 * tmp_45;
      real_t tmp_165 = tmp_162 + tmp_163 + tmp_164;
      real_t tmp_166 = tmp_155 * tmp_48;
      real_t tmp_167 = tmp_157 * tmp_50;
      real_t tmp_168 = tmp_159 * tmp_52;
      real_t tmp_169 = tmp_166 + tmp_167 + tmp_168;
      real_t tmp_170 = tmp_1 * ( tmp_161 - 1.0 / 4.0 ) + tmp_14 * ( tmp_169 - 1.0 / 4.0 ) + tmp_18 * ( tmp_165 - 1.0 / 4.0 );
      real_t tmp_171 = tmp_170 * tmp_61;
      real_t tmp_172 = 0.020848748529055869 * tmp_63;
      real_t tmp_173 = 0.78764240869137092 * tmp_3 + 0.041227165399737475 * tmp_4 + tmp_5;
      real_t tmp_174 = tmp_173 * tmp_24;
      real_t tmp_175 = 0.78764240869137092 * tmp_27 + 0.041227165399737475 * tmp_28 + tmp_29;
      real_t tmp_176 = tmp_175 * tmp_31;
      real_t tmp_177 = 0.78764240869137092 * tmp_34 + 0.041227165399737475 * tmp_35 + tmp_36;
      real_t tmp_178 = tmp_177 * tmp_38;
      real_t tmp_179 = tmp_174 + tmp_176 + tmp_178;
      real_t tmp_180 = tmp_173 * tmp_41;
      real_t tmp_181 = tmp_175 * tmp_43;
      real_t tmp_182 = tmp_177 * tmp_45;
      real_t tmp_183 = tmp_180 + tmp_181 + tmp_182;
      real_t tmp_184 = tmp_173 * tmp_48;
      real_t tmp_185 = tmp_175 * tmp_50;
      real_t tmp_186 = tmp_177 * tmp_52;
      real_t tmp_187 = tmp_184 + tmp_185 + tmp_186;
      real_t tmp_188 = tmp_1 * ( tmp_179 - 1.0 / 4.0 ) + tmp_14 * ( tmp_187 - 1.0 / 4.0 ) + tmp_18 * ( tmp_183 - 1.0 / 4.0 );
      real_t tmp_189 = tmp_188 * tmp_61;
      real_t tmp_190 = 0.019202922745021479 * tmp_63;
      real_t tmp_191 = 0.58463275527740355 * tmp_3 + 0.039308471900058539 * tmp_4 + tmp_5;
      real_t tmp_192 = tmp_191 * tmp_24;
      real_t tmp_193 = 0.58463275527740355 * tmp_27 + 0.039308471900058539 * tmp_28 + tmp_29;
      real_t tmp_194 = tmp_193 * tmp_31;
      real_t tmp_195 = 0.58463275527740355 * tmp_34 + 0.039308471900058539 * tmp_35 + tmp_36;
      real_t tmp_196 = tmp_195 * tmp_38;
      real_t tmp_197 = tmp_192 + tmp_194 + tmp_196;
      real_t tmp_198 = tmp_191 * tmp_41;
      real_t tmp_199 = tmp_193 * tmp_43;
      real_t tmp_200 = tmp_195 * tmp_45;
      real_t tmp_201 = tmp_198 + tmp_199 + tmp_200;
      real_t tmp_202 = tmp_191 * tmp_48;
      real_t tmp_203 = tmp_193 * tmp_50;
      real_t tmp_204 = tmp_195 * tmp_52;
      real_t tmp_205 = tmp_202 + tmp_203 + tmp_204;
      real_t tmp_206 = tmp_1 * ( tmp_197 - 1.0 / 4.0 ) + tmp_14 * ( tmp_205 - 1.0 / 4.0 ) + tmp_18 * ( tmp_201 - 1.0 / 4.0 );
      real_t tmp_207 = tmp_206 * tmp_61;
      real_t tmp_208 = 0.020848748529055869 * tmp_63;
      real_t tmp_209 = 0.1711304259088916 * tmp_3 + 0.78764240869137092 * tmp_4 + tmp_5;
      real_t tmp_210 = tmp_209 * tmp_24;
      real_t tmp_211 = 0.1711304259088916 * tmp_27 + 0.78764240869137092 * tmp_28 + tmp_29;
      real_t tmp_212 = tmp_211 * tmp_31;
      real_t tmp_213 = 0.1711304259088916 * tmp_34 + 0.78764240869137092 * tmp_35 + tmp_36;
      real_t tmp_214 = tmp_213 * tmp_38;
      real_t tmp_215 = tmp_210 + tmp_212 + tmp_214;
      real_t tmp_216 = tmp_209 * tmp_41;
      real_t tmp_217 = tmp_211 * tmp_43;
      real_t tmp_218 = tmp_213 * tmp_45;
      real_t tmp_219 = tmp_216 + tmp_217 + tmp_218;
      real_t tmp_220 = tmp_209 * tmp_48;
      real_t tmp_221 = tmp_211 * tmp_50;
      real_t tmp_222 = tmp_213 * tmp_52;
      real_t tmp_223 = tmp_220 + tmp_221 + tmp_222;
      real_t tmp_224 = tmp_1 * ( tmp_215 - 1.0 / 4.0 ) + tmp_14 * ( tmp_223 - 1.0 / 4.0 ) + tmp_18 * ( tmp_219 - 1.0 / 4.0 );
      real_t tmp_225 = tmp_224 * tmp_61;
      real_t tmp_226 = 0.019202922745021479 * tmp_63;
      real_t tmp_227 = 0.37605877282253791 * tmp_3 + 0.58463275527740355 * tmp_4 + tmp_5;
      real_t tmp_228 = tmp_227 * tmp_24;
      real_t tmp_229 = 0.37605877282253791 * tmp_27 + 0.58463275527740355 * tmp_28 + tmp_29;
      real_t tmp_230 = tmp_229 * tmp_31;
      real_t tmp_231 = 0.37605877282253791 * tmp_34 + 0.58463275527740355 * tmp_35 + tmp_36;
      real_t tmp_232 = tmp_231 * tmp_38;
      real_t tmp_233 = tmp_228 + tmp_230 + tmp_232;
      real_t tmp_234 = tmp_227 * tmp_41;
      real_t tmp_235 = tmp_229 * tmp_43;
      real_t tmp_236 = tmp_231 * tmp_45;
      real_t tmp_237 = tmp_234 + tmp_235 + tmp_236;
      real_t tmp_238 = tmp_227 * tmp_48;
      real_t tmp_239 = tmp_229 * tmp_50;
      real_t tmp_240 = tmp_231 * tmp_52;
      real_t tmp_241 = tmp_238 + tmp_239 + tmp_240;
      real_t tmp_242 = tmp_1 * ( tmp_233 - 1.0 / 4.0 ) + tmp_14 * ( tmp_241 - 1.0 / 4.0 ) + tmp_18 * ( tmp_237 - 1.0 / 4.0 );
      real_t tmp_243 = tmp_242 * tmp_61;
      real_t tmp_244 = 0.020848748529055869 * tmp_63;
      real_t tmp_245 = 0.041227165399737475 * tmp_3 + 0.1711304259088916 * tmp_4 + tmp_5;
      real_t tmp_246 = tmp_24 * tmp_245;
      real_t tmp_247 = 0.041227165399737475 * tmp_27 + 0.1711304259088916 * tmp_28 + tmp_29;
      real_t tmp_248 = tmp_247 * tmp_31;
      real_t tmp_249 = 0.041227165399737475 * tmp_34 + 0.1711304259088916 * tmp_35 + tmp_36;
      real_t tmp_250 = tmp_249 * tmp_38;
      real_t tmp_251 = tmp_246 + tmp_248 + tmp_250;
      real_t tmp_252 = tmp_245 * tmp_41;
      real_t tmp_253 = tmp_247 * tmp_43;
      real_t tmp_254 = tmp_249 * tmp_45;
      real_t tmp_255 = tmp_252 + tmp_253 + tmp_254;
      real_t tmp_256 = tmp_245 * tmp_48;
      real_t tmp_257 = tmp_247 * tmp_50;
      real_t tmp_258 = tmp_249 * tmp_52;
      real_t tmp_259 = tmp_256 + tmp_257 + tmp_258;
      real_t tmp_260 = tmp_1 * ( tmp_251 - 1.0 / 4.0 ) + tmp_14 * ( tmp_259 - 1.0 / 4.0 ) + tmp_18 * ( tmp_255 - 1.0 / 4.0 );
      real_t tmp_261 = tmp_260 * tmp_61;
      real_t tmp_262 = 0.019202922745021479 * tmp_63;
      real_t tmp_263 = 0.40446199974765351 * tmp_3 + 0.19107600050469298 * tmp_4 + tmp_5;
      real_t tmp_264 = tmp_24 * tmp_263;
      real_t tmp_265 = 0.40446199974765351 * tmp_27 + 0.19107600050469298 * tmp_28 + tmp_29;
      real_t tmp_266 = tmp_265 * tmp_31;
      real_t tmp_267 = 0.40446199974765351 * tmp_34 + 0.19107600050469298 * tmp_35 + tmp_36;
      real_t tmp_268 = tmp_267 * tmp_38;
      real_t tmp_269 = tmp_264 + tmp_266 + tmp_268;
      real_t tmp_270 = tmp_263 * tmp_41;
      real_t tmp_271 = tmp_265 * tmp_43;
      real_t tmp_272 = tmp_267 * tmp_45;
      real_t tmp_273 = tmp_270 + tmp_271 + tmp_272;
      real_t tmp_274 = tmp_263 * tmp_48;
      real_t tmp_275 = tmp_265 * tmp_50;
      real_t tmp_276 = tmp_267 * tmp_52;
      real_t tmp_277 = tmp_274 + tmp_275 + tmp_276;
      real_t tmp_278 = tmp_1 * ( tmp_269 - 1.0 / 4.0 ) + tmp_14 * ( tmp_277 - 1.0 / 4.0 ) + tmp_18 * ( tmp_273 - 1.0 / 4.0 );
      real_t tmp_279 = tmp_278 * tmp_61;
      real_t tmp_280 = 0.042507265838595799 * tmp_63;
      real_t tmp_281 = 0.039308471900058539 * tmp_3 + 0.37605877282253791 * tmp_4 + tmp_5;
      real_t tmp_282 = tmp_24 * tmp_281;
      real_t tmp_283 = 0.039308471900058539 * tmp_27 + 0.37605877282253791 * tmp_28 + tmp_29;
      real_t tmp_284 = tmp_283 * tmp_31;
      real_t tmp_285 = 0.039308471900058539 * tmp_34 + 0.37605877282253791 * tmp_35 + tmp_36;
      real_t tmp_286 = tmp_285 * tmp_38;
      real_t tmp_287 = tmp_282 + tmp_284 + tmp_286;
      real_t tmp_288 = tmp_281 * tmp_41;
      real_t tmp_289 = tmp_283 * tmp_43;
      real_t tmp_290 = tmp_285 * tmp_45;
      real_t tmp_291 = tmp_288 + tmp_289 + tmp_290;
      real_t tmp_292 = tmp_281 * tmp_48;
      real_t tmp_293 = tmp_283 * tmp_50;
      real_t tmp_294 = tmp_285 * tmp_52;
      real_t tmp_295 = tmp_292 + tmp_293 + tmp_294;
      real_t tmp_296 = tmp_1 * ( tmp_287 - 1.0 / 4.0 ) + tmp_14 * ( tmp_295 - 1.0 / 4.0 ) + tmp_18 * ( tmp_291 - 1.0 / 4.0 );
      real_t tmp_297 = tmp_296 * tmp_61;
      real_t tmp_298 = 0.020848748529055869 * tmp_63;
      real_t tmp_299 = 0.93718850182767688 * tmp_3 + 0.031405749086161582 * tmp_4 + tmp_5;
      real_t tmp_300 = tmp_24 * tmp_299;
      real_t tmp_301 = 0.93718850182767688 * tmp_27 + 0.031405749086161582 * tmp_28 + tmp_29;
      real_t tmp_302 = tmp_301 * tmp_31;
      real_t tmp_303 = 0.93718850182767688 * tmp_34 + 0.031405749086161582 * tmp_35 + tmp_36;
      real_t tmp_304 = tmp_303 * tmp_38;
      real_t tmp_305 = tmp_300 + tmp_302 + tmp_304;
      real_t tmp_306 = tmp_299 * tmp_41;
      real_t tmp_307 = tmp_301 * tmp_43;
      real_t tmp_308 = tmp_303 * tmp_45;
      real_t tmp_309 = tmp_306 + tmp_307 + tmp_308;
      real_t tmp_310 = tmp_299 * tmp_48;
      real_t tmp_311 = tmp_301 * tmp_50;
      real_t tmp_312 = tmp_303 * tmp_52;
      real_t tmp_313 = tmp_310 + tmp_311 + tmp_312;
      real_t tmp_314 = tmp_1 * ( tmp_305 - 1.0 / 4.0 ) + tmp_14 * ( tmp_313 - 1.0 / 4.0 ) + tmp_18 * ( tmp_309 - 1.0 / 4.0 );
      real_t tmp_315 = tmp_314 * tmp_61;
      real_t tmp_316 = 0.0068572537431980923 * tmp_63;
      real_t tmp_317 = 0.60796128279561268 * tmp_3 + 0.19601935860219369 * tmp_4 + tmp_5;
      real_t tmp_318 = tmp_24 * tmp_317;
      real_t tmp_319 = 0.60796128279561268 * tmp_27 + 0.19601935860219369 * tmp_28 + tmp_29;
      real_t tmp_320 = tmp_31 * tmp_319;
      real_t tmp_321 = 0.60796128279561268 * tmp_34 + 0.19601935860219369 * tmp_35 + tmp_36;
      real_t tmp_322 = tmp_321 * tmp_38;
      real_t tmp_323 = tmp_318 + tmp_320 + tmp_322;
      real_t tmp_324 = tmp_317 * tmp_41;
      real_t tmp_325 = tmp_319 * tmp_43;
      real_t tmp_326 = tmp_321 * tmp_45;
      real_t tmp_327 = tmp_324 + tmp_325 + tmp_326;
      real_t tmp_328 = tmp_317 * tmp_48;
      real_t tmp_329 = tmp_319 * tmp_50;
      real_t tmp_330 = tmp_321 * tmp_52;
      real_t tmp_331 = tmp_328 + tmp_329 + tmp_330;
      real_t tmp_332 = tmp_1 * ( tmp_323 - 1.0 / 4.0 ) + tmp_14 * ( tmp_331 - 1.0 / 4.0 ) + tmp_18 * ( tmp_327 - 1.0 / 4.0 );
      real_t tmp_333 = tmp_332 * tmp_61;
      real_t tmp_334 = 0.037198804536718075 * tmp_63;
      real_t tmp_335 = 0.19107600050469298 * tmp_3 + 0.40446199974765351 * tmp_4 + tmp_5;
      real_t tmp_336 = tmp_24 * tmp_335;
      real_t tmp_337 = 0.19107600050469298 * tmp_27 + 0.40446199974765351 * tmp_28 + tmp_29;
      real_t tmp_338 = tmp_31 * tmp_337;
      real_t tmp_339 = 0.19107600050469298 * tmp_34 + 0.40446199974765351 * tmp_35 + tmp_36;
      real_t tmp_340 = tmp_339 * tmp_38;
      real_t tmp_341 = tmp_336 + tmp_338 + tmp_340;
      real_t tmp_342 = tmp_335 * tmp_41;
      real_t tmp_343 = tmp_337 * tmp_43;
      real_t tmp_344 = tmp_339 * tmp_45;
      real_t tmp_345 = tmp_342 + tmp_343 + tmp_344;
      real_t tmp_346 = tmp_335 * tmp_48;
      real_t tmp_347 = tmp_337 * tmp_50;
      real_t tmp_348 = tmp_339 * tmp_52;
      real_t tmp_349 = tmp_346 + tmp_347 + tmp_348;
      real_t tmp_350 = tmp_1 * ( tmp_341 - 1.0 / 4.0 ) + tmp_14 * ( tmp_349 - 1.0 / 4.0 ) + tmp_18 * ( tmp_345 - 1.0 / 4.0 );
      real_t tmp_351 = tmp_350 * tmp_61;
      real_t tmp_352 = 0.042507265838595799 * tmp_63;
      real_t tmp_353 = 0.031405749086161582 * tmp_3 + 0.031405749086161582 * tmp_4 + tmp_5;
      real_t tmp_354 = tmp_24 * tmp_353;
      real_t tmp_355 = 0.031405749086161582 * tmp_27 + 0.031405749086161582 * tmp_28 + tmp_29;
      real_t tmp_356 = tmp_31 * tmp_355;
      real_t tmp_357 = 0.031405749086161582 * tmp_34 + 0.031405749086161582 * tmp_35 + tmp_36;
      real_t tmp_358 = tmp_357 * tmp_38;
      real_t tmp_359 = tmp_354 + tmp_356 + tmp_358;
      real_t tmp_360 = tmp_353 * tmp_41;
      real_t tmp_361 = tmp_355 * tmp_43;
      real_t tmp_362 = tmp_357 * tmp_45;
      real_t tmp_363 = tmp_360 + tmp_361 + tmp_362;
      real_t tmp_364 = tmp_353 * tmp_48;
      real_t tmp_365 = tmp_355 * tmp_50;
      real_t tmp_366 = tmp_357 * tmp_52;
      real_t tmp_367 = tmp_364 + tmp_365 + tmp_366;
      real_t tmp_368 = tmp_1 * ( tmp_359 - 1.0 / 4.0 ) + tmp_14 * ( tmp_367 - 1.0 / 4.0 ) + tmp_18 * ( tmp_363 - 1.0 / 4.0 );
      real_t tmp_369 = tmp_368 * tmp_61;
      real_t tmp_370 = 0.0068572537431980923 * tmp_63;
      real_t tmp_371 = 0.19601935860219369 * tmp_3 + 0.19601935860219369 * tmp_4 + tmp_5;
      real_t tmp_372 = tmp_24 * tmp_371;
      real_t tmp_373 = 0.19601935860219369 * tmp_27 + 0.19601935860219369 * tmp_28 + tmp_29;
      real_t tmp_374 = tmp_31 * tmp_373;
      real_t tmp_375 = 0.19601935860219369 * tmp_34 + 0.19601935860219369 * tmp_35 + tmp_36;
      real_t tmp_376 = tmp_375 * tmp_38;
      real_t tmp_377 = tmp_372 + tmp_374 + tmp_376;
      real_t tmp_378 = tmp_371 * tmp_41;
      real_t tmp_379 = tmp_373 * tmp_43;
      real_t tmp_380 = tmp_375 * tmp_45;
      real_t tmp_381 = tmp_378 + tmp_379 + tmp_380;
      real_t tmp_382 = tmp_371 * tmp_48;
      real_t tmp_383 = tmp_373 * tmp_50;
      real_t tmp_384 = tmp_375 * tmp_52;
      real_t tmp_385 = tmp_382 + tmp_383 + tmp_384;
      real_t tmp_386 = tmp_1 * ( tmp_377 - 1.0 / 4.0 ) + tmp_14 * ( tmp_385 - 1.0 / 4.0 ) + tmp_18 * ( tmp_381 - 1.0 / 4.0 );
      real_t tmp_387 = tmp_386 * tmp_61;
      real_t tmp_388 = 0.037198804536718075 * tmp_63;
      real_t tmp_389 = 0.40446199974765351 * tmp_3 + 0.40446199974765351 * tmp_4 + tmp_5;
      real_t tmp_390 = tmp_24 * tmp_389;
      real_t tmp_391 = 0.40446199974765351 * tmp_27 + 0.40446199974765351 * tmp_28 + tmp_29;
      real_t tmp_392 = tmp_31 * tmp_391;
      real_t tmp_393 = 0.40446199974765351 * tmp_34 + 0.40446199974765351 * tmp_35 + tmp_36;
      real_t tmp_394 = tmp_38 * tmp_393;
      real_t tmp_395 = tmp_390 + tmp_392 + tmp_394;
      real_t tmp_396 = tmp_389 * tmp_41;
      real_t tmp_397 = tmp_391 * tmp_43;
      real_t tmp_398 = tmp_393 * tmp_45;
      real_t tmp_399 = tmp_396 + tmp_397 + tmp_398;
      real_t tmp_400 = tmp_389 * tmp_48;
      real_t tmp_401 = tmp_391 * tmp_50;
      real_t tmp_402 = tmp_393 * tmp_52;
      real_t tmp_403 = tmp_400 + tmp_401 + tmp_402;
      real_t tmp_404 = tmp_1 * ( tmp_395 - 1.0 / 4.0 ) + tmp_14 * ( tmp_403 - 1.0 / 4.0 ) + tmp_18 * ( tmp_399 - 1.0 / 4.0 );
      real_t tmp_405 = tmp_404 * tmp_61;
      real_t tmp_406 = 0.042507265838595799 * tmp_63;
      real_t tmp_407 = 0.1711304259088916 * tmp_3 + 0.041227165399737475 * tmp_4 + tmp_5;
      real_t tmp_408 = tmp_24 * tmp_407;
      real_t tmp_409 = 0.1711304259088916 * tmp_27 + 0.041227165399737475 * tmp_28 + tmp_29;
      real_t tmp_410 = tmp_31 * tmp_409;
      real_t tmp_411 = 0.1711304259088916 * tmp_34 + 0.041227165399737475 * tmp_35 + tmp_36;
      real_t tmp_412 = tmp_38 * tmp_411;
      real_t tmp_413 = tmp_408 + tmp_410 + tmp_412;
      real_t tmp_414 = tmp_407 * tmp_41;
      real_t tmp_415 = tmp_409 * tmp_43;
      real_t tmp_416 = tmp_411 * tmp_45;
      real_t tmp_417 = tmp_414 + tmp_415 + tmp_416;
      real_t tmp_418 = tmp_407 * tmp_48;
      real_t tmp_419 = tmp_409 * tmp_50;
      real_t tmp_420 = tmp_411 * tmp_52;
      real_t tmp_421 = tmp_418 + tmp_419 + tmp_420;
      real_t tmp_422 = tmp_1 * ( tmp_413 - 1.0 / 4.0 ) + tmp_14 * ( tmp_421 - 1.0 / 4.0 ) + tmp_18 * ( tmp_417 - 1.0 / 4.0 );
      real_t tmp_423 = tmp_422 * tmp_61;
      real_t tmp_424 = 0.019202922745021479 * tmp_63;
      real_t tmp_425 = 0.5 * p_affine_13_0 * tmp_38 + 0.5 * p_affine_13_1 * tmp_31 + 0.5 * p_affine_13_2 * tmp_24;
      real_t tmp_426 = 0.5 * p_affine_13_0 * tmp_45 + 0.5 * p_affine_13_1 * tmp_43 + 0.5 * p_affine_13_2 * tmp_41;
      real_t tmp_427 = 0.5 * p_affine_13_0 * tmp_52 + 0.5 * p_affine_13_1 * tmp_50 + 0.5 * p_affine_13_2 * tmp_48;
      real_t a_0_0 =
          tmp_100 * ( -tmp_56 * tmp_98 +
                      tmp_99 * ( -tmp_84 - tmp_86 - tmp_88 - tmp_90 - tmp_91 - tmp_92 - tmp_94 - tmp_95 - tmp_96 + 1 ) ) +
          tmp_118 * ( -tmp_116 * tmp_56 + tmp_117 * ( -tmp_102 - tmp_104 - tmp_106 - tmp_108 - tmp_109 - tmp_110 - tmp_112 -
                                                      tmp_113 - tmp_114 + 1 ) ) +
          tmp_136 * ( -tmp_134 * tmp_56 + tmp_135 * ( -tmp_120 - tmp_122 - tmp_124 - tmp_126 - tmp_127 - tmp_128 - tmp_130 -
                                                      tmp_131 - tmp_132 + 1 ) ) +
          tmp_154 * ( -tmp_152 * tmp_56 + tmp_153 * ( -tmp_138 - tmp_140 - tmp_142 - tmp_144 - tmp_145 - tmp_146 - tmp_148 -
                                                      tmp_149 - tmp_150 + 1 ) ) +
          tmp_172 * ( -tmp_170 * tmp_56 + tmp_171 * ( -tmp_156 - tmp_158 - tmp_160 - tmp_162 - tmp_163 - tmp_164 - tmp_166 -
                                                      tmp_167 - tmp_168 + 1 ) ) +
          tmp_190 * ( -tmp_188 * tmp_56 + tmp_189 * ( -tmp_174 - tmp_176 - tmp_178 - tmp_180 - tmp_181 - tmp_182 - tmp_184 -
                                                      tmp_185 - tmp_186 + 1 ) ) +
          tmp_208 * ( -tmp_206 * tmp_56 + tmp_207 * ( -tmp_192 - tmp_194 - tmp_196 - tmp_198 - tmp_199 - tmp_200 - tmp_202 -
                                                      tmp_203 - tmp_204 + 1 ) ) +
          tmp_226 * ( -tmp_224 * tmp_56 + tmp_225 * ( -tmp_210 - tmp_212 - tmp_214 - tmp_216 - tmp_217 - tmp_218 - tmp_220 -
                                                      tmp_221 - tmp_222 + 1 ) ) +
          tmp_244 * ( -tmp_242 * tmp_56 + tmp_243 * ( -tmp_228 - tmp_230 - tmp_232 - tmp_234 - tmp_235 - tmp_236 - tmp_238 -
                                                      tmp_239 - tmp_240 + 1 ) ) +
          tmp_262 * ( -tmp_260 * tmp_56 + tmp_261 * ( -tmp_246 - tmp_248 - tmp_250 - tmp_252 - tmp_253 - tmp_254 - tmp_256 -
                                                      tmp_257 - tmp_258 + 1 ) ) +
          tmp_280 * ( -tmp_278 * tmp_56 + tmp_279 * ( -tmp_264 - tmp_266 - tmp_268 - tmp_270 - tmp_271 - tmp_272 - tmp_274 -
                                                      tmp_275 - tmp_276 + 1 ) ) +
          tmp_298 * ( -tmp_296 * tmp_56 + tmp_297 * ( -tmp_282 - tmp_284 - tmp_286 - tmp_288 - tmp_289 - tmp_290 - tmp_292 -
                                                      tmp_293 - tmp_294 + 1 ) ) +
          tmp_316 * ( -tmp_314 * tmp_56 + tmp_315 * ( -tmp_300 - tmp_302 - tmp_304 - tmp_306 - tmp_307 - tmp_308 - tmp_310 -
                                                      tmp_311 - tmp_312 + 1 ) ) +
          tmp_334 * ( -tmp_332 * tmp_56 + tmp_333 * ( -tmp_318 - tmp_320 - tmp_322 - tmp_324 - tmp_325 - tmp_326 - tmp_328 -
                                                      tmp_329 - tmp_330 + 1 ) ) +
          tmp_352 * ( -tmp_350 * tmp_56 + tmp_351 * ( -tmp_336 - tmp_338 - tmp_340 - tmp_342 - tmp_343 - tmp_344 - tmp_346 -
                                                      tmp_347 - tmp_348 + 1 ) ) +
          tmp_370 * ( -tmp_368 * tmp_56 + tmp_369 * ( -tmp_354 - tmp_356 - tmp_358 - tmp_360 - tmp_361 - tmp_362 - tmp_364 -
                                                      tmp_365 - tmp_366 + 1 ) ) +
          tmp_388 * ( -tmp_386 * tmp_56 + tmp_387 * ( -tmp_372 - tmp_374 - tmp_376 - tmp_378 - tmp_379 - tmp_380 - tmp_382 -
                                                      tmp_383 - tmp_384 + 1 ) ) +
          tmp_406 * ( -tmp_404 * tmp_56 + tmp_405 * ( -tmp_390 - tmp_392 - tmp_394 - tmp_396 - tmp_397 - tmp_398 - tmp_400 -
                                                      tmp_401 - tmp_402 + 1 ) ) +
          tmp_424 * ( -tmp_422 * tmp_56 + tmp_423 * ( -tmp_408 - tmp_410 - tmp_412 - tmp_414 - tmp_415 - tmp_416 - tmp_418 -
                                                      tmp_419 - tmp_420 + 1 ) ) +
          tmp_64 * ( -tmp_55 * tmp_56 +
                     tmp_62 * ( -tmp_25 - tmp_32 - tmp_39 - tmp_42 - tmp_44 - tmp_46 - tmp_49 - tmp_51 - tmp_53 + 1 ) ) +
          tmp_82 * ( -tmp_56 * tmp_80 +
                     tmp_81 * ( -tmp_66 - tmp_68 - tmp_70 - tmp_72 - tmp_73 - tmp_74 - tmp_76 - tmp_77 - tmp_78 + 1 ) );
      real_t a_0_1 = tmp_100 * ( -tmp_425 * tmp_98 + tmp_89 * tmp_99 ) + tmp_118 * ( tmp_107 * tmp_117 - tmp_116 * tmp_425 ) +
                     tmp_136 * ( tmp_125 * tmp_135 - tmp_134 * tmp_425 ) + tmp_154 * ( tmp_143 * tmp_153 - tmp_152 * tmp_425 ) +
                     tmp_172 * ( tmp_161 * tmp_171 - tmp_170 * tmp_425 ) + tmp_190 * ( tmp_179 * tmp_189 - tmp_188 * tmp_425 ) +
                     tmp_208 * ( tmp_197 * tmp_207 - tmp_206 * tmp_425 ) + tmp_226 * ( tmp_215 * tmp_225 - tmp_224 * tmp_425 ) +
                     tmp_244 * ( tmp_233 * tmp_243 - tmp_242 * tmp_425 ) + tmp_262 * ( tmp_251 * tmp_261 - tmp_260 * tmp_425 ) +
                     tmp_280 * ( tmp_269 * tmp_279 - tmp_278 * tmp_425 ) + tmp_298 * ( tmp_287 * tmp_297 - tmp_296 * tmp_425 ) +
                     tmp_316 * ( tmp_305 * tmp_315 - tmp_314 * tmp_425 ) + tmp_334 * ( tmp_323 * tmp_333 - tmp_332 * tmp_425 ) +
                     tmp_352 * ( tmp_341 * tmp_351 - tmp_350 * tmp_425 ) + tmp_370 * ( tmp_359 * tmp_369 - tmp_368 * tmp_425 ) +
                     tmp_388 * ( tmp_377 * tmp_387 - tmp_386 * tmp_425 ) + tmp_406 * ( tmp_395 * tmp_405 - tmp_404 * tmp_425 ) +
                     tmp_424 * ( tmp_413 * tmp_423 - tmp_422 * tmp_425 ) + tmp_64 * ( tmp_40 * tmp_62 - tmp_425 * tmp_55 ) +
                     tmp_82 * ( -tmp_425 * tmp_80 + tmp_71 * tmp_81 );
      real_t a_0_2 = tmp_100 * ( -tmp_426 * tmp_98 + tmp_93 * tmp_99 ) + tmp_118 * ( tmp_111 * tmp_117 - tmp_116 * tmp_426 ) +
                     tmp_136 * ( tmp_129 * tmp_135 - tmp_134 * tmp_426 ) + tmp_154 * ( tmp_147 * tmp_153 - tmp_152 * tmp_426 ) +
                     tmp_172 * ( tmp_165 * tmp_171 - tmp_170 * tmp_426 ) + tmp_190 * ( tmp_183 * tmp_189 - tmp_188 * tmp_426 ) +
                     tmp_208 * ( tmp_201 * tmp_207 - tmp_206 * tmp_426 ) + tmp_226 * ( tmp_219 * tmp_225 - tmp_224 * tmp_426 ) +
                     tmp_244 * ( tmp_237 * tmp_243 - tmp_242 * tmp_426 ) + tmp_262 * ( tmp_255 * tmp_261 - tmp_260 * tmp_426 ) +
                     tmp_280 * ( tmp_273 * tmp_279 - tmp_278 * tmp_426 ) + tmp_298 * ( tmp_291 * tmp_297 - tmp_296 * tmp_426 ) +
                     tmp_316 * ( tmp_309 * tmp_315 - tmp_314 * tmp_426 ) + tmp_334 * ( tmp_327 * tmp_333 - tmp_332 * tmp_426 ) +
                     tmp_352 * ( tmp_345 * tmp_351 - tmp_350 * tmp_426 ) + tmp_370 * ( tmp_363 * tmp_369 - tmp_368 * tmp_426 ) +
                     tmp_388 * ( tmp_381 * tmp_387 - tmp_386 * tmp_426 ) + tmp_406 * ( tmp_399 * tmp_405 - tmp_404 * tmp_426 ) +
                     tmp_424 * ( tmp_417 * tmp_423 - tmp_422 * tmp_426 ) + tmp_64 * ( -tmp_426 * tmp_55 + tmp_47 * tmp_62 ) +
                     tmp_82 * ( -tmp_426 * tmp_80 + tmp_75 * tmp_81 );
      real_t a_0_3 = tmp_100 * ( -tmp_427 * tmp_98 + tmp_97 * tmp_99 ) + tmp_118 * ( tmp_115 * tmp_117 - tmp_116 * tmp_427 ) +
                     tmp_136 * ( tmp_133 * tmp_135 - tmp_134 * tmp_427 ) + tmp_154 * ( tmp_151 * tmp_153 - tmp_152 * tmp_427 ) +
                     tmp_172 * ( tmp_169 * tmp_171 - tmp_170 * tmp_427 ) + tmp_190 * ( tmp_187 * tmp_189 - tmp_188 * tmp_427 ) +
                     tmp_208 * ( tmp_205 * tmp_207 - tmp_206 * tmp_427 ) + tmp_226 * ( tmp_223 * tmp_225 - tmp_224 * tmp_427 ) +
                     tmp_244 * ( tmp_241 * tmp_243 - tmp_242 * tmp_427 ) + tmp_262 * ( tmp_259 * tmp_261 - tmp_260 * tmp_427 ) +
                     tmp_280 * ( tmp_277 * tmp_279 - tmp_278 * tmp_427 ) + tmp_298 * ( tmp_295 * tmp_297 - tmp_296 * tmp_427 ) +
                     tmp_316 * ( tmp_313 * tmp_315 - tmp_314 * tmp_427 ) + tmp_334 * ( tmp_331 * tmp_333 - tmp_332 * tmp_427 ) +
                     tmp_352 * ( tmp_349 * tmp_351 - tmp_350 * tmp_427 ) + tmp_370 * ( tmp_367 * tmp_369 - tmp_368 * tmp_427 ) +
                     tmp_388 * ( tmp_385 * tmp_387 - tmp_386 * tmp_427 ) + tmp_406 * ( tmp_403 * tmp_405 - tmp_404 * tmp_427 ) +
                     tmp_424 * ( tmp_421 * tmp_423 - tmp_422 * tmp_427 ) + tmp_64 * ( -tmp_427 * tmp_55 + tmp_54 * tmp_62 ) +
                     tmp_82 * ( -tmp_427 * tmp_80 + tmp_79 * tmp_81 );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 0, 3 ) = a_0_3;
   }

   void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                  const std::vector< Point3D >& coordsElementOuter,
                                  const std::vector< Point3D >& coordsFacet,
                                  const Point3D&,
                                  const Point3D&,
                                  const Point3D&     outwardNormal,
                                  const DGBasisInfo& trialBasis,
                                  const DGBasisInfo& testBasis,
                                  int                trialDegree,
                                  int                testDegree,
                                  MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_2;
      real_t tmp_1  = p_affine_1_2 + tmp_0;
      real_t tmp_2  = -p_affine_0_0;
      real_t tmp_3  = p_affine_2_0 + tmp_2;
      real_t tmp_4  = -p_affine_0_1;
      real_t tmp_5  = p_affine_3_1 + tmp_4;
      real_t tmp_6  = tmp_3 * tmp_5;
      real_t tmp_7  = p_affine_3_0 + tmp_2;
      real_t tmp_8  = p_affine_2_1 + tmp_4;
      real_t tmp_9  = tmp_7 * tmp_8;
      real_t tmp_10 = tmp_6 - tmp_9;
      real_t tmp_11 = p_affine_1_0 + tmp_2;
      real_t tmp_12 = p_affine_3_2 + tmp_0;
      real_t tmp_13 = tmp_12 * tmp_8;
      real_t tmp_14 = p_affine_1_1 + tmp_4;
      real_t tmp_15 = p_affine_2_2 + tmp_0;
      real_t tmp_16 = tmp_15 * tmp_7;
      real_t tmp_17 = tmp_15 * tmp_5;
      real_t tmp_18 = tmp_12 * tmp_3;
      real_t tmp_19 =
          1.0 / ( tmp_1 * tmp_6 - tmp_1 * tmp_9 + tmp_11 * tmp_13 - tmp_11 * tmp_17 + tmp_14 * tmp_16 - tmp_14 * tmp_18 );
      real_t tmp_20 = p_affine_8_2 + tmp_0;
      real_t tmp_21 = -p_affine_8_2;
      real_t tmp_22 = p_affine_9_2 + tmp_21;
      real_t tmp_23 = p_affine_10_2 + tmp_21;
      real_t tmp_24 = 0.031405749086161582 * tmp_22 + 0.93718850182767688 * tmp_23;
      real_t tmp_25 = tmp_19 * ( tmp_20 + tmp_24 );
      real_t tmp_26 = tmp_16 - tmp_18;
      real_t tmp_27 = p_affine_8_1 + tmp_4;
      real_t tmp_28 = -p_affine_8_1;
      real_t tmp_29 = p_affine_9_1 + tmp_28;
      real_t tmp_30 = p_affine_10_1 + tmp_28;
      real_t tmp_31 = 0.031405749086161582 * tmp_29 + 0.93718850182767688 * tmp_30;
      real_t tmp_32 = tmp_19 * ( tmp_27 + tmp_31 );
      real_t tmp_33 = tmp_13 - tmp_17;
      real_t tmp_34 = p_affine_8_0 + tmp_2;
      real_t tmp_35 = -p_affine_8_0;
      real_t tmp_36 = p_affine_9_0 + tmp_35;
      real_t tmp_37 = p_affine_10_0 + tmp_35;
      real_t tmp_38 = 0.031405749086161582 * tmp_36 + 0.93718850182767688 * tmp_37;
      real_t tmp_39 = tmp_19 * ( tmp_34 + tmp_38 );
      real_t tmp_40 = -tmp_11 * tmp_5 + tmp_14 * tmp_7;
      real_t tmp_41 = -tmp_1 * tmp_7 + tmp_11 * tmp_12;
      real_t tmp_42 = tmp_1 * tmp_5 - tmp_12 * tmp_14;
      real_t tmp_43 = tmp_11 * tmp_8 - tmp_14 * tmp_3;
      real_t tmp_44 = tmp_1 * tmp_3 - tmp_11 * tmp_15;
      real_t tmp_45 = -tmp_1 * tmp_8 + tmp_14 * tmp_15;
      real_t tmp_46 = tmp_1 * ( tmp_10 * tmp_25 + tmp_26 * tmp_32 + tmp_33 * tmp_39 - 1.0 / 4.0 ) +
                      tmp_12 * ( tmp_25 * tmp_43 + tmp_32 * tmp_44 + tmp_39 * tmp_45 - 1.0 / 4.0 ) +
                      tmp_15 * ( tmp_25 * tmp_40 + tmp_32 * tmp_41 + tmp_39 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_47 = -p_affine_4_1;
      real_t tmp_48 = p_affine_5_1 + tmp_47;
      real_t tmp_49 = -p_affine_4_2;
      real_t tmp_50 = p_affine_6_2 + tmp_49;
      real_t tmp_51 = tmp_48 * tmp_50;
      real_t tmp_52 = p_affine_6_1 + tmp_47;
      real_t tmp_53 = p_affine_5_2 + tmp_49;
      real_t tmp_54 = tmp_52 * tmp_53;
      real_t tmp_55 = -p_affine_4_0;
      real_t tmp_56 = p_affine_5_0 + tmp_55;
      real_t tmp_57 = p_affine_7_2 + tmp_49;
      real_t tmp_58 = tmp_52 * tmp_57;
      real_t tmp_59 = p_affine_6_0 + tmp_55;
      real_t tmp_60 = p_affine_7_1 + tmp_47;
      real_t tmp_61 = tmp_53 * tmp_60;
      real_t tmp_62 = p_affine_7_0 + tmp_55;
      real_t tmp_63 = tmp_50 * tmp_60;
      real_t tmp_64 = tmp_48 * tmp_57;
      real_t tmp_65 =
          1.0 / ( tmp_51 * tmp_62 - tmp_54 * tmp_62 + tmp_56 * tmp_58 - tmp_56 * tmp_63 + tmp_59 * tmp_61 - tmp_59 * tmp_64 );
      real_t tmp_66 = tmp_65 * ( tmp_51 - tmp_54 );
      real_t tmp_67 = tmp_65 * ( tmp_61 - tmp_64 );
      real_t tmp_68 = tmp_65 * ( tmp_58 - tmp_63 );
      real_t tmp_69 = tmp_65 * ( -tmp_50 * tmp_56 + tmp_53 * tmp_59 );
      real_t tmp_70 = tmp_65 * ( -tmp_53 * tmp_62 + tmp_56 * tmp_57 );
      real_t tmp_71 = tmp_65 * ( tmp_50 * tmp_62 - tmp_57 * tmp_59 );
      real_t tmp_72 = tmp_65 * ( -tmp_48 * tmp_59 + tmp_52 * tmp_56 );
      real_t tmp_73 = tmp_65 * ( tmp_48 * tmp_62 - tmp_56 * tmp_60 );
      real_t tmp_74 = tmp_65 * ( -tmp_52 * tmp_62 + tmp_59 * tmp_60 );
      real_t tmp_75 = 0.5 * p_affine_13_0 * ( -tmp_66 - tmp_67 - tmp_68 ) + 0.5 * p_affine_13_1 * ( -tmp_69 - tmp_70 - tmp_71 ) +
                      0.5 * p_affine_13_2 * ( -tmp_72 - tmp_73 - tmp_74 );
      real_t tmp_76 = p_affine_8_2 + tmp_49;
      real_t tmp_77 = tmp_24 + tmp_76;
      real_t tmp_78 = tmp_72 * tmp_77;
      real_t tmp_79 = tmp_73 * tmp_77;
      real_t tmp_80 = p_affine_8_1 + tmp_47;
      real_t tmp_81 = tmp_31 + tmp_80;
      real_t tmp_82 = tmp_69 * tmp_81;
      real_t tmp_83 = tmp_70 * tmp_81;
      real_t tmp_84 = tmp_74 * tmp_77;
      real_t tmp_85 = tmp_71 * tmp_81;
      real_t tmp_86 = p_affine_8_0 + tmp_55;
      real_t tmp_87 = tmp_38 + tmp_86;
      real_t tmp_88 = tmp_66 * tmp_87;
      real_t tmp_89 = tmp_67 * tmp_87;
      real_t tmp_90 = tmp_68 * tmp_87;
      real_t tmp_91 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_92 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_93 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_94 = ( std::abs( tmp_23 * tmp_91 - tmp_30 * tmp_93 ) * std::abs( tmp_23 * tmp_91 - tmp_30 * tmp_93 ) ) +
                      ( std::abs( tmp_23 * tmp_92 - tmp_37 * tmp_93 ) * std::abs( tmp_23 * tmp_92 - tmp_37 * tmp_93 ) ) +
                      ( std::abs( tmp_30 * tmp_92 - tmp_37 * tmp_91 ) * std::abs( tmp_30 * tmp_92 - tmp_37 * tmp_91 ) );
      real_t tmp_95  = 1.0 * std::pow( tmp_94, -0.25 );
      real_t tmp_96  = tmp_46 * tmp_95;
      real_t tmp_97  = 1.0 * std::pow( tmp_94, 1.0 / 2.0 );
      real_t tmp_98  = 0.0068572537431980923 * tmp_97;
      real_t tmp_99  = 0.19601935860219369 * tmp_22 + 0.60796128279561268 * tmp_23;
      real_t tmp_100 = tmp_19 * ( tmp_20 + tmp_99 );
      real_t tmp_101 = 0.19601935860219369 * tmp_29 + 0.60796128279561268 * tmp_30;
      real_t tmp_102 = tmp_19 * ( tmp_101 + tmp_27 );
      real_t tmp_103 = 0.19601935860219369 * tmp_36 + 0.60796128279561268 * tmp_37;
      real_t tmp_104 = tmp_19 * ( tmp_103 + tmp_34 );
      real_t tmp_105 = tmp_1 * ( tmp_10 * tmp_100 + tmp_102 * tmp_26 + tmp_104 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_100 * tmp_43 + tmp_102 * tmp_44 + tmp_104 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_100 * tmp_40 + tmp_102 * tmp_41 + tmp_104 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_106 = tmp_76 + tmp_99;
      real_t tmp_107 = tmp_106 * tmp_72;
      real_t tmp_108 = tmp_106 * tmp_73;
      real_t tmp_109 = tmp_101 + tmp_80;
      real_t tmp_110 = tmp_109 * tmp_69;
      real_t tmp_111 = tmp_109 * tmp_70;
      real_t tmp_112 = tmp_106 * tmp_74;
      real_t tmp_113 = tmp_109 * tmp_71;
      real_t tmp_114 = tmp_103 + tmp_86;
      real_t tmp_115 = tmp_114 * tmp_66;
      real_t tmp_116 = tmp_114 * tmp_67;
      real_t tmp_117 = tmp_114 * tmp_68;
      real_t tmp_118 = tmp_105 * tmp_95;
      real_t tmp_119 = 0.037198804536718075 * tmp_97;
      real_t tmp_120 = 0.37605877282253791 * tmp_22 + 0.039308471900058539 * tmp_23;
      real_t tmp_121 = tmp_19 * ( tmp_120 + tmp_20 );
      real_t tmp_122 = 0.37605877282253791 * tmp_29 + 0.039308471900058539 * tmp_30;
      real_t tmp_123 = tmp_19 * ( tmp_122 + tmp_27 );
      real_t tmp_124 = 0.37605877282253791 * tmp_36 + 0.039308471900058539 * tmp_37;
      real_t tmp_125 = tmp_19 * ( tmp_124 + tmp_34 );
      real_t tmp_126 = tmp_1 * ( tmp_10 * tmp_121 + tmp_123 * tmp_26 + tmp_125 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_121 * tmp_43 + tmp_123 * tmp_44 + tmp_125 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_121 * tmp_40 + tmp_123 * tmp_41 + tmp_125 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_127 = tmp_120 + tmp_76;
      real_t tmp_128 = tmp_127 * tmp_72;
      real_t tmp_129 = tmp_127 * tmp_73;
      real_t tmp_130 = tmp_122 + tmp_80;
      real_t tmp_131 = tmp_130 * tmp_69;
      real_t tmp_132 = tmp_130 * tmp_70;
      real_t tmp_133 = tmp_127 * tmp_74;
      real_t tmp_134 = tmp_130 * tmp_71;
      real_t tmp_135 = tmp_124 + tmp_86;
      real_t tmp_136 = tmp_135 * tmp_66;
      real_t tmp_137 = tmp_135 * tmp_67;
      real_t tmp_138 = tmp_135 * tmp_68;
      real_t tmp_139 = tmp_126 * tmp_95;
      real_t tmp_140 = 0.020848748529055869 * tmp_97;
      real_t tmp_141 = 0.78764240869137092 * tmp_22 + 0.1711304259088916 * tmp_23;
      real_t tmp_142 = tmp_19 * ( tmp_141 + tmp_20 );
      real_t tmp_143 = 0.78764240869137092 * tmp_29 + 0.1711304259088916 * tmp_30;
      real_t tmp_144 = tmp_19 * ( tmp_143 + tmp_27 );
      real_t tmp_145 = 0.78764240869137092 * tmp_36 + 0.1711304259088916 * tmp_37;
      real_t tmp_146 = tmp_19 * ( tmp_145 + tmp_34 );
      real_t tmp_147 = tmp_1 * ( tmp_10 * tmp_142 + tmp_144 * tmp_26 + tmp_146 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_142 * tmp_43 + tmp_144 * tmp_44 + tmp_146 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_142 * tmp_40 + tmp_144 * tmp_41 + tmp_146 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_148 = tmp_141 + tmp_76;
      real_t tmp_149 = tmp_148 * tmp_72;
      real_t tmp_150 = tmp_148 * tmp_73;
      real_t tmp_151 = tmp_143 + tmp_80;
      real_t tmp_152 = tmp_151 * tmp_69;
      real_t tmp_153 = tmp_151 * tmp_70;
      real_t tmp_154 = tmp_148 * tmp_74;
      real_t tmp_155 = tmp_151 * tmp_71;
      real_t tmp_156 = tmp_145 + tmp_86;
      real_t tmp_157 = tmp_156 * tmp_66;
      real_t tmp_158 = tmp_156 * tmp_67;
      real_t tmp_159 = tmp_156 * tmp_68;
      real_t tmp_160 = tmp_147 * tmp_95;
      real_t tmp_161 = 0.019202922745021479 * tmp_97;
      real_t tmp_162 = 0.58463275527740355 * tmp_22 + 0.37605877282253791 * tmp_23;
      real_t tmp_163 = tmp_19 * ( tmp_162 + tmp_20 );
      real_t tmp_164 = 0.58463275527740355 * tmp_29 + 0.37605877282253791 * tmp_30;
      real_t tmp_165 = tmp_19 * ( tmp_164 + tmp_27 );
      real_t tmp_166 = 0.58463275527740355 * tmp_36 + 0.37605877282253791 * tmp_37;
      real_t tmp_167 = tmp_19 * ( tmp_166 + tmp_34 );
      real_t tmp_168 = tmp_1 * ( tmp_10 * tmp_163 + tmp_165 * tmp_26 + tmp_167 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_163 * tmp_43 + tmp_165 * tmp_44 + tmp_167 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_163 * tmp_40 + tmp_165 * tmp_41 + tmp_167 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_169 = tmp_162 + tmp_76;
      real_t tmp_170 = tmp_169 * tmp_72;
      real_t tmp_171 = tmp_169 * tmp_73;
      real_t tmp_172 = tmp_164 + tmp_80;
      real_t tmp_173 = tmp_172 * tmp_69;
      real_t tmp_174 = tmp_172 * tmp_70;
      real_t tmp_175 = tmp_169 * tmp_74;
      real_t tmp_176 = tmp_172 * tmp_71;
      real_t tmp_177 = tmp_166 + tmp_86;
      real_t tmp_178 = tmp_177 * tmp_66;
      real_t tmp_179 = tmp_177 * tmp_67;
      real_t tmp_180 = tmp_177 * tmp_68;
      real_t tmp_181 = tmp_168 * tmp_95;
      real_t tmp_182 = 0.020848748529055869 * tmp_97;
      real_t tmp_183 = 0.041227165399737475 * tmp_22 + 0.78764240869137092 * tmp_23;
      real_t tmp_184 = tmp_19 * ( tmp_183 + tmp_20 );
      real_t tmp_185 = 0.041227165399737475 * tmp_29 + 0.78764240869137092 * tmp_30;
      real_t tmp_186 = tmp_19 * ( tmp_185 + tmp_27 );
      real_t tmp_187 = 0.041227165399737475 * tmp_36 + 0.78764240869137092 * tmp_37;
      real_t tmp_188 = tmp_19 * ( tmp_187 + tmp_34 );
      real_t tmp_189 = tmp_1 * ( tmp_10 * tmp_184 + tmp_186 * tmp_26 + tmp_188 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_184 * tmp_43 + tmp_186 * tmp_44 + tmp_188 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_184 * tmp_40 + tmp_186 * tmp_41 + tmp_188 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_190 = tmp_183 + tmp_76;
      real_t tmp_191 = tmp_190 * tmp_72;
      real_t tmp_192 = tmp_190 * tmp_73;
      real_t tmp_193 = tmp_185 + tmp_80;
      real_t tmp_194 = tmp_193 * tmp_69;
      real_t tmp_195 = tmp_193 * tmp_70;
      real_t tmp_196 = tmp_190 * tmp_74;
      real_t tmp_197 = tmp_193 * tmp_71;
      real_t tmp_198 = tmp_187 + tmp_86;
      real_t tmp_199 = tmp_198 * tmp_66;
      real_t tmp_200 = tmp_198 * tmp_67;
      real_t tmp_201 = tmp_198 * tmp_68;
      real_t tmp_202 = tmp_189 * tmp_95;
      real_t tmp_203 = 0.019202922745021479 * tmp_97;
      real_t tmp_204 = 0.039308471900058539 * tmp_22 + 0.58463275527740355 * tmp_23;
      real_t tmp_205 = tmp_19 * ( tmp_20 + tmp_204 );
      real_t tmp_206 = 0.039308471900058539 * tmp_29 + 0.58463275527740355 * tmp_30;
      real_t tmp_207 = tmp_19 * ( tmp_206 + tmp_27 );
      real_t tmp_208 = 0.039308471900058539 * tmp_36 + 0.58463275527740355 * tmp_37;
      real_t tmp_209 = tmp_19 * ( tmp_208 + tmp_34 );
      real_t tmp_210 = tmp_1 * ( tmp_10 * tmp_205 + tmp_207 * tmp_26 + tmp_209 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_205 * tmp_43 + tmp_207 * tmp_44 + tmp_209 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_205 * tmp_40 + tmp_207 * tmp_41 + tmp_209 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_211 = tmp_204 + tmp_76;
      real_t tmp_212 = tmp_211 * tmp_72;
      real_t tmp_213 = tmp_211 * tmp_73;
      real_t tmp_214 = tmp_206 + tmp_80;
      real_t tmp_215 = tmp_214 * tmp_69;
      real_t tmp_216 = tmp_214 * tmp_70;
      real_t tmp_217 = tmp_211 * tmp_74;
      real_t tmp_218 = tmp_214 * tmp_71;
      real_t tmp_219 = tmp_208 + tmp_86;
      real_t tmp_220 = tmp_219 * tmp_66;
      real_t tmp_221 = tmp_219 * tmp_67;
      real_t tmp_222 = tmp_219 * tmp_68;
      real_t tmp_223 = tmp_210 * tmp_95;
      real_t tmp_224 = 0.020848748529055869 * tmp_97;
      real_t tmp_225 = 0.78764240869137092 * tmp_22 + 0.041227165399737475 * tmp_23;
      real_t tmp_226 = tmp_19 * ( tmp_20 + tmp_225 );
      real_t tmp_227 = 0.78764240869137092 * tmp_29 + 0.041227165399737475 * tmp_30;
      real_t tmp_228 = tmp_19 * ( tmp_227 + tmp_27 );
      real_t tmp_229 = 0.78764240869137092 * tmp_36 + 0.041227165399737475 * tmp_37;
      real_t tmp_230 = tmp_19 * ( tmp_229 + tmp_34 );
      real_t tmp_231 = tmp_1 * ( tmp_10 * tmp_226 + tmp_228 * tmp_26 + tmp_230 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_226 * tmp_43 + tmp_228 * tmp_44 + tmp_230 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_226 * tmp_40 + tmp_228 * tmp_41 + tmp_230 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_232 = tmp_225 + tmp_76;
      real_t tmp_233 = tmp_232 * tmp_72;
      real_t tmp_234 = tmp_232 * tmp_73;
      real_t tmp_235 = tmp_227 + tmp_80;
      real_t tmp_236 = tmp_235 * tmp_69;
      real_t tmp_237 = tmp_235 * tmp_70;
      real_t tmp_238 = tmp_232 * tmp_74;
      real_t tmp_239 = tmp_235 * tmp_71;
      real_t tmp_240 = tmp_229 + tmp_86;
      real_t tmp_241 = tmp_240 * tmp_66;
      real_t tmp_242 = tmp_240 * tmp_67;
      real_t tmp_243 = tmp_240 * tmp_68;
      real_t tmp_244 = tmp_231 * tmp_95;
      real_t tmp_245 = 0.019202922745021479 * tmp_97;
      real_t tmp_246 = 0.58463275527740355 * tmp_22 + 0.039308471900058539 * tmp_23;
      real_t tmp_247 = tmp_19 * ( tmp_20 + tmp_246 );
      real_t tmp_248 = 0.58463275527740355 * tmp_29 + 0.039308471900058539 * tmp_30;
      real_t tmp_249 = tmp_19 * ( tmp_248 + tmp_27 );
      real_t tmp_250 = 0.58463275527740355 * tmp_36 + 0.039308471900058539 * tmp_37;
      real_t tmp_251 = tmp_19 * ( tmp_250 + tmp_34 );
      real_t tmp_252 = tmp_1 * ( tmp_10 * tmp_247 + tmp_249 * tmp_26 + tmp_251 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_247 * tmp_43 + tmp_249 * tmp_44 + tmp_251 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_247 * tmp_40 + tmp_249 * tmp_41 + tmp_251 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_253 = tmp_246 + tmp_76;
      real_t tmp_254 = tmp_253 * tmp_72;
      real_t tmp_255 = tmp_253 * tmp_73;
      real_t tmp_256 = tmp_248 + tmp_80;
      real_t tmp_257 = tmp_256 * tmp_69;
      real_t tmp_258 = tmp_256 * tmp_70;
      real_t tmp_259 = tmp_253 * tmp_74;
      real_t tmp_260 = tmp_256 * tmp_71;
      real_t tmp_261 = tmp_250 + tmp_86;
      real_t tmp_262 = tmp_261 * tmp_66;
      real_t tmp_263 = tmp_261 * tmp_67;
      real_t tmp_264 = tmp_261 * tmp_68;
      real_t tmp_265 = tmp_252 * tmp_95;
      real_t tmp_266 = 0.020848748529055869 * tmp_97;
      real_t tmp_267 = 0.1711304259088916 * tmp_22 + 0.78764240869137092 * tmp_23;
      real_t tmp_268 = tmp_19 * ( tmp_20 + tmp_267 );
      real_t tmp_269 = 0.1711304259088916 * tmp_29 + 0.78764240869137092 * tmp_30;
      real_t tmp_270 = tmp_19 * ( tmp_269 + tmp_27 );
      real_t tmp_271 = 0.1711304259088916 * tmp_36 + 0.78764240869137092 * tmp_37;
      real_t tmp_272 = tmp_19 * ( tmp_271 + tmp_34 );
      real_t tmp_273 = tmp_1 * ( tmp_10 * tmp_268 + tmp_26 * tmp_270 + tmp_272 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_268 * tmp_43 + tmp_270 * tmp_44 + tmp_272 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_268 * tmp_40 + tmp_270 * tmp_41 + tmp_272 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_274 = tmp_267 + tmp_76;
      real_t tmp_275 = tmp_274 * tmp_72;
      real_t tmp_276 = tmp_274 * tmp_73;
      real_t tmp_277 = tmp_269 + tmp_80;
      real_t tmp_278 = tmp_277 * tmp_69;
      real_t tmp_279 = tmp_277 * tmp_70;
      real_t tmp_280 = tmp_274 * tmp_74;
      real_t tmp_281 = tmp_277 * tmp_71;
      real_t tmp_282 = tmp_271 + tmp_86;
      real_t tmp_283 = tmp_282 * tmp_66;
      real_t tmp_284 = tmp_282 * tmp_67;
      real_t tmp_285 = tmp_282 * tmp_68;
      real_t tmp_286 = tmp_273 * tmp_95;
      real_t tmp_287 = 0.019202922745021479 * tmp_97;
      real_t tmp_288 = 0.37605877282253791 * tmp_22 + 0.58463275527740355 * tmp_23;
      real_t tmp_289 = tmp_19 * ( tmp_20 + tmp_288 );
      real_t tmp_290 = 0.37605877282253791 * tmp_29 + 0.58463275527740355 * tmp_30;
      real_t tmp_291 = tmp_19 * ( tmp_27 + tmp_290 );
      real_t tmp_292 = 0.37605877282253791 * tmp_36 + 0.58463275527740355 * tmp_37;
      real_t tmp_293 = tmp_19 * ( tmp_292 + tmp_34 );
      real_t tmp_294 = tmp_1 * ( tmp_10 * tmp_289 + tmp_26 * tmp_291 + tmp_293 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_289 * tmp_43 + tmp_291 * tmp_44 + tmp_293 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_289 * tmp_40 + tmp_291 * tmp_41 + tmp_293 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_295 = tmp_288 + tmp_76;
      real_t tmp_296 = tmp_295 * tmp_72;
      real_t tmp_297 = tmp_295 * tmp_73;
      real_t tmp_298 = tmp_290 + tmp_80;
      real_t tmp_299 = tmp_298 * tmp_69;
      real_t tmp_300 = tmp_298 * tmp_70;
      real_t tmp_301 = tmp_295 * tmp_74;
      real_t tmp_302 = tmp_298 * tmp_71;
      real_t tmp_303 = tmp_292 + tmp_86;
      real_t tmp_304 = tmp_303 * tmp_66;
      real_t tmp_305 = tmp_303 * tmp_67;
      real_t tmp_306 = tmp_303 * tmp_68;
      real_t tmp_307 = tmp_294 * tmp_95;
      real_t tmp_308 = 0.020848748529055869 * tmp_97;
      real_t tmp_309 = 0.041227165399737475 * tmp_22 + 0.1711304259088916 * tmp_23;
      real_t tmp_310 = tmp_19 * ( tmp_20 + tmp_309 );
      real_t tmp_311 = 0.041227165399737475 * tmp_29 + 0.1711304259088916 * tmp_30;
      real_t tmp_312 = tmp_19 * ( tmp_27 + tmp_311 );
      real_t tmp_313 = 0.041227165399737475 * tmp_36 + 0.1711304259088916 * tmp_37;
      real_t tmp_314 = tmp_19 * ( tmp_313 + tmp_34 );
      real_t tmp_315 = tmp_1 * ( tmp_10 * tmp_310 + tmp_26 * tmp_312 + tmp_314 * tmp_33 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_310 * tmp_43 + tmp_312 * tmp_44 + tmp_314 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_310 * tmp_40 + tmp_312 * tmp_41 + tmp_314 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_316 = tmp_309 + tmp_76;
      real_t tmp_317 = tmp_316 * tmp_72;
      real_t tmp_318 = tmp_316 * tmp_73;
      real_t tmp_319 = tmp_311 + tmp_80;
      real_t tmp_320 = tmp_319 * tmp_69;
      real_t tmp_321 = tmp_319 * tmp_70;
      real_t tmp_322 = tmp_316 * tmp_74;
      real_t tmp_323 = tmp_319 * tmp_71;
      real_t tmp_324 = tmp_313 + tmp_86;
      real_t tmp_325 = tmp_324 * tmp_66;
      real_t tmp_326 = tmp_324 * tmp_67;
      real_t tmp_327 = tmp_324 * tmp_68;
      real_t tmp_328 = tmp_315 * tmp_95;
      real_t tmp_329 = 0.019202922745021479 * tmp_97;
      real_t tmp_330 = 0.40446199974765351 * tmp_22 + 0.19107600050469298 * tmp_23;
      real_t tmp_331 = tmp_19 * ( tmp_20 + tmp_330 );
      real_t tmp_332 = 0.40446199974765351 * tmp_29 + 0.19107600050469298 * tmp_30;
      real_t tmp_333 = tmp_19 * ( tmp_27 + tmp_332 );
      real_t tmp_334 = 0.40446199974765351 * tmp_36 + 0.19107600050469298 * tmp_37;
      real_t tmp_335 = tmp_19 * ( tmp_334 + tmp_34 );
      real_t tmp_336 = tmp_1 * ( tmp_10 * tmp_331 + tmp_26 * tmp_333 + tmp_33 * tmp_335 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_331 * tmp_43 + tmp_333 * tmp_44 + tmp_335 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_331 * tmp_40 + tmp_333 * tmp_41 + tmp_335 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_337 = tmp_330 + tmp_76;
      real_t tmp_338 = tmp_337 * tmp_72;
      real_t tmp_339 = tmp_337 * tmp_73;
      real_t tmp_340 = tmp_332 + tmp_80;
      real_t tmp_341 = tmp_340 * tmp_69;
      real_t tmp_342 = tmp_340 * tmp_70;
      real_t tmp_343 = tmp_337 * tmp_74;
      real_t tmp_344 = tmp_340 * tmp_71;
      real_t tmp_345 = tmp_334 + tmp_86;
      real_t tmp_346 = tmp_345 * tmp_66;
      real_t tmp_347 = tmp_345 * tmp_67;
      real_t tmp_348 = tmp_345 * tmp_68;
      real_t tmp_349 = tmp_336 * tmp_95;
      real_t tmp_350 = 0.042507265838595799 * tmp_97;
      real_t tmp_351 = 0.039308471900058539 * tmp_22 + 0.37605877282253791 * tmp_23;
      real_t tmp_352 = tmp_19 * ( tmp_20 + tmp_351 );
      real_t tmp_353 = 0.039308471900058539 * tmp_29 + 0.37605877282253791 * tmp_30;
      real_t tmp_354 = tmp_19 * ( tmp_27 + tmp_353 );
      real_t tmp_355 = 0.039308471900058539 * tmp_36 + 0.37605877282253791 * tmp_37;
      real_t tmp_356 = tmp_19 * ( tmp_34 + tmp_355 );
      real_t tmp_357 = tmp_1 * ( tmp_10 * tmp_352 + tmp_26 * tmp_354 + tmp_33 * tmp_356 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_352 * tmp_43 + tmp_354 * tmp_44 + tmp_356 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_352 * tmp_40 + tmp_354 * tmp_41 + tmp_356 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_358 = tmp_351 + tmp_76;
      real_t tmp_359 = tmp_358 * tmp_72;
      real_t tmp_360 = tmp_358 * tmp_73;
      real_t tmp_361 = tmp_353 + tmp_80;
      real_t tmp_362 = tmp_361 * tmp_69;
      real_t tmp_363 = tmp_361 * tmp_70;
      real_t tmp_364 = tmp_358 * tmp_74;
      real_t tmp_365 = tmp_361 * tmp_71;
      real_t tmp_366 = tmp_355 + tmp_86;
      real_t tmp_367 = tmp_366 * tmp_66;
      real_t tmp_368 = tmp_366 * tmp_67;
      real_t tmp_369 = tmp_366 * tmp_68;
      real_t tmp_370 = tmp_357 * tmp_95;
      real_t tmp_371 = 0.020848748529055869 * tmp_97;
      real_t tmp_372 = 0.93718850182767688 * tmp_22 + 0.031405749086161582 * tmp_23;
      real_t tmp_373 = tmp_19 * ( tmp_20 + tmp_372 );
      real_t tmp_374 = 0.93718850182767688 * tmp_29 + 0.031405749086161582 * tmp_30;
      real_t tmp_375 = tmp_19 * ( tmp_27 + tmp_374 );
      real_t tmp_376 = 0.93718850182767688 * tmp_36 + 0.031405749086161582 * tmp_37;
      real_t tmp_377 = tmp_19 * ( tmp_34 + tmp_376 );
      real_t tmp_378 = tmp_1 * ( tmp_10 * tmp_373 + tmp_26 * tmp_375 + tmp_33 * tmp_377 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_373 * tmp_43 + tmp_375 * tmp_44 + tmp_377 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_373 * tmp_40 + tmp_375 * tmp_41 + tmp_377 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_379 = tmp_372 + tmp_76;
      real_t tmp_380 = tmp_379 * tmp_72;
      real_t tmp_381 = tmp_379 * tmp_73;
      real_t tmp_382 = tmp_374 + tmp_80;
      real_t tmp_383 = tmp_382 * tmp_69;
      real_t tmp_384 = tmp_382 * tmp_70;
      real_t tmp_385 = tmp_379 * tmp_74;
      real_t tmp_386 = tmp_382 * tmp_71;
      real_t tmp_387 = tmp_376 + tmp_86;
      real_t tmp_388 = tmp_387 * tmp_66;
      real_t tmp_389 = tmp_387 * tmp_67;
      real_t tmp_390 = tmp_387 * tmp_68;
      real_t tmp_391 = tmp_378 * tmp_95;
      real_t tmp_392 = 0.0068572537431980923 * tmp_97;
      real_t tmp_393 = 0.60796128279561268 * tmp_22 + 0.19601935860219369 * tmp_23;
      real_t tmp_394 = tmp_19 * ( tmp_20 + tmp_393 );
      real_t tmp_395 = 0.60796128279561268 * tmp_29 + 0.19601935860219369 * tmp_30;
      real_t tmp_396 = tmp_19 * ( tmp_27 + tmp_395 );
      real_t tmp_397 = 0.60796128279561268 * tmp_36 + 0.19601935860219369 * tmp_37;
      real_t tmp_398 = tmp_19 * ( tmp_34 + tmp_397 );
      real_t tmp_399 = tmp_1 * ( tmp_10 * tmp_394 + tmp_26 * tmp_396 + tmp_33 * tmp_398 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_394 * tmp_43 + tmp_396 * tmp_44 + tmp_398 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_394 * tmp_40 + tmp_396 * tmp_41 + tmp_398 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_400 = tmp_393 + tmp_76;
      real_t tmp_401 = tmp_400 * tmp_72;
      real_t tmp_402 = tmp_400 * tmp_73;
      real_t tmp_403 = tmp_395 + tmp_80;
      real_t tmp_404 = tmp_403 * tmp_69;
      real_t tmp_405 = tmp_403 * tmp_70;
      real_t tmp_406 = tmp_400 * tmp_74;
      real_t tmp_407 = tmp_403 * tmp_71;
      real_t tmp_408 = tmp_397 + tmp_86;
      real_t tmp_409 = tmp_408 * tmp_66;
      real_t tmp_410 = tmp_408 * tmp_67;
      real_t tmp_411 = tmp_408 * tmp_68;
      real_t tmp_412 = tmp_399 * tmp_95;
      real_t tmp_413 = 0.037198804536718075 * tmp_97;
      real_t tmp_414 = 0.19107600050469298 * tmp_22 + 0.40446199974765351 * tmp_23;
      real_t tmp_415 = tmp_19 * ( tmp_20 + tmp_414 );
      real_t tmp_416 = 0.19107600050469298 * tmp_29 + 0.40446199974765351 * tmp_30;
      real_t tmp_417 = tmp_19 * ( tmp_27 + tmp_416 );
      real_t tmp_418 = 0.19107600050469298 * tmp_36 + 0.40446199974765351 * tmp_37;
      real_t tmp_419 = tmp_19 * ( tmp_34 + tmp_418 );
      real_t tmp_420 = tmp_1 * ( tmp_10 * tmp_415 + tmp_26 * tmp_417 + tmp_33 * tmp_419 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_415 * tmp_43 + tmp_417 * tmp_44 + tmp_419 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_40 * tmp_415 + tmp_41 * tmp_417 + tmp_419 * tmp_42 - 1.0 / 4.0 );
      real_t tmp_421 = tmp_414 + tmp_76;
      real_t tmp_422 = tmp_421 * tmp_72;
      real_t tmp_423 = tmp_421 * tmp_73;
      real_t tmp_424 = tmp_416 + tmp_80;
      real_t tmp_425 = tmp_424 * tmp_69;
      real_t tmp_426 = tmp_424 * tmp_70;
      real_t tmp_427 = tmp_421 * tmp_74;
      real_t tmp_428 = tmp_424 * tmp_71;
      real_t tmp_429 = tmp_418 + tmp_86;
      real_t tmp_430 = tmp_429 * tmp_66;
      real_t tmp_431 = tmp_429 * tmp_67;
      real_t tmp_432 = tmp_429 * tmp_68;
      real_t tmp_433 = tmp_420 * tmp_95;
      real_t tmp_434 = 0.042507265838595799 * tmp_97;
      real_t tmp_435 = 0.031405749086161582 * tmp_22 + 0.031405749086161582 * tmp_23;
      real_t tmp_436 = tmp_19 * ( tmp_20 + tmp_435 );
      real_t tmp_437 = 0.031405749086161582 * tmp_29 + 0.031405749086161582 * tmp_30;
      real_t tmp_438 = tmp_19 * ( tmp_27 + tmp_437 );
      real_t tmp_439 = 0.031405749086161582 * tmp_36 + 0.031405749086161582 * tmp_37;
      real_t tmp_440 = tmp_19 * ( tmp_34 + tmp_439 );
      real_t tmp_441 = tmp_1 * ( tmp_10 * tmp_436 + tmp_26 * tmp_438 + tmp_33 * tmp_440 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_43 * tmp_436 + tmp_438 * tmp_44 + tmp_440 * tmp_45 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_40 * tmp_436 + tmp_41 * tmp_438 + tmp_42 * tmp_440 - 1.0 / 4.0 );
      real_t tmp_442 = tmp_435 + tmp_76;
      real_t tmp_443 = tmp_442 * tmp_72;
      real_t tmp_444 = tmp_442 * tmp_73;
      real_t tmp_445 = tmp_437 + tmp_80;
      real_t tmp_446 = tmp_445 * tmp_69;
      real_t tmp_447 = tmp_445 * tmp_70;
      real_t tmp_448 = tmp_442 * tmp_74;
      real_t tmp_449 = tmp_445 * tmp_71;
      real_t tmp_450 = tmp_439 + tmp_86;
      real_t tmp_451 = tmp_450 * tmp_66;
      real_t tmp_452 = tmp_450 * tmp_67;
      real_t tmp_453 = tmp_450 * tmp_68;
      real_t tmp_454 = tmp_441 * tmp_95;
      real_t tmp_455 = 0.0068572537431980923 * tmp_97;
      real_t tmp_456 = 0.19601935860219369 * tmp_22 + 0.19601935860219369 * tmp_23;
      real_t tmp_457 = tmp_19 * ( tmp_20 + tmp_456 );
      real_t tmp_458 = 0.19601935860219369 * tmp_29 + 0.19601935860219369 * tmp_30;
      real_t tmp_459 = tmp_19 * ( tmp_27 + tmp_458 );
      real_t tmp_460 = 0.19601935860219369 * tmp_36 + 0.19601935860219369 * tmp_37;
      real_t tmp_461 = tmp_19 * ( tmp_34 + tmp_460 );
      real_t tmp_462 = tmp_1 * ( tmp_10 * tmp_457 + tmp_26 * tmp_459 + tmp_33 * tmp_461 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_43 * tmp_457 + tmp_44 * tmp_459 + tmp_45 * tmp_461 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_40 * tmp_457 + tmp_41 * tmp_459 + tmp_42 * tmp_461 - 1.0 / 4.0 );
      real_t tmp_463 = tmp_456 + tmp_76;
      real_t tmp_464 = tmp_463 * tmp_72;
      real_t tmp_465 = tmp_463 * tmp_73;
      real_t tmp_466 = tmp_458 + tmp_80;
      real_t tmp_467 = tmp_466 * tmp_69;
      real_t tmp_468 = tmp_466 * tmp_70;
      real_t tmp_469 = tmp_463 * tmp_74;
      real_t tmp_470 = tmp_466 * tmp_71;
      real_t tmp_471 = tmp_460 + tmp_86;
      real_t tmp_472 = tmp_471 * tmp_66;
      real_t tmp_473 = tmp_471 * tmp_67;
      real_t tmp_474 = tmp_471 * tmp_68;
      real_t tmp_475 = tmp_462 * tmp_95;
      real_t tmp_476 = 0.037198804536718075 * tmp_97;
      real_t tmp_477 = 0.40446199974765351 * tmp_22 + 0.40446199974765351 * tmp_23;
      real_t tmp_478 = tmp_19 * ( tmp_20 + tmp_477 );
      real_t tmp_479 = 0.40446199974765351 * tmp_29 + 0.40446199974765351 * tmp_30;
      real_t tmp_480 = tmp_19 * ( tmp_27 + tmp_479 );
      real_t tmp_481 = 0.40446199974765351 * tmp_36 + 0.40446199974765351 * tmp_37;
      real_t tmp_482 = tmp_19 * ( tmp_34 + tmp_481 );
      real_t tmp_483 = tmp_1 * ( tmp_10 * tmp_478 + tmp_26 * tmp_480 + tmp_33 * tmp_482 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_43 * tmp_478 + tmp_44 * tmp_480 + tmp_45 * tmp_482 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_40 * tmp_478 + tmp_41 * tmp_480 + tmp_42 * tmp_482 - 1.0 / 4.0 );
      real_t tmp_484 = tmp_477 + tmp_76;
      real_t tmp_485 = tmp_484 * tmp_72;
      real_t tmp_486 = tmp_484 * tmp_73;
      real_t tmp_487 = tmp_479 + tmp_80;
      real_t tmp_488 = tmp_487 * tmp_69;
      real_t tmp_489 = tmp_487 * tmp_70;
      real_t tmp_490 = tmp_484 * tmp_74;
      real_t tmp_491 = tmp_487 * tmp_71;
      real_t tmp_492 = tmp_481 + tmp_86;
      real_t tmp_493 = tmp_492 * tmp_66;
      real_t tmp_494 = tmp_492 * tmp_67;
      real_t tmp_495 = tmp_492 * tmp_68;
      real_t tmp_496 = tmp_483 * tmp_95;
      real_t tmp_497 = 0.042507265838595799 * tmp_97;
      real_t tmp_498 = 0.1711304259088916 * tmp_22 + 0.041227165399737475 * tmp_23;
      real_t tmp_499 = tmp_19 * ( tmp_20 + tmp_498 );
      real_t tmp_500 = 0.1711304259088916 * tmp_29 + 0.041227165399737475 * tmp_30;
      real_t tmp_501 = tmp_19 * ( tmp_27 + tmp_500 );
      real_t tmp_502 = 0.1711304259088916 * tmp_36 + 0.041227165399737475 * tmp_37;
      real_t tmp_503 = tmp_19 * ( tmp_34 + tmp_502 );
      real_t tmp_504 = tmp_1 * ( tmp_10 * tmp_499 + tmp_26 * tmp_501 + tmp_33 * tmp_503 - 1.0 / 4.0 ) +
                       tmp_12 * ( tmp_43 * tmp_499 + tmp_44 * tmp_501 + tmp_45 * tmp_503 - 1.0 / 4.0 ) +
                       tmp_15 * ( tmp_40 * tmp_499 + tmp_41 * tmp_501 + tmp_42 * tmp_503 - 1.0 / 4.0 );
      real_t tmp_505 = tmp_498 + tmp_76;
      real_t tmp_506 = tmp_505 * tmp_72;
      real_t tmp_507 = tmp_505 * tmp_73;
      real_t tmp_508 = tmp_500 + tmp_80;
      real_t tmp_509 = tmp_508 * tmp_69;
      real_t tmp_510 = tmp_508 * tmp_70;
      real_t tmp_511 = tmp_505 * tmp_74;
      real_t tmp_512 = tmp_508 * tmp_71;
      real_t tmp_513 = tmp_502 + tmp_86;
      real_t tmp_514 = tmp_513 * tmp_66;
      real_t tmp_515 = tmp_513 * tmp_67;
      real_t tmp_516 = tmp_513 * tmp_68;
      real_t tmp_517 = tmp_504 * tmp_95;
      real_t tmp_518 = 0.019202922745021479 * tmp_97;
      real_t tmp_519 = 0.5 * p_affine_13_0 * tmp_68 + 0.5 * p_affine_13_1 * tmp_71 + 0.5 * p_affine_13_2 * tmp_74;
      real_t tmp_520 = 0.5 * p_affine_13_0 * tmp_67 + 0.5 * p_affine_13_1 * tmp_70 + 0.5 * p_affine_13_2 * tmp_73;
      real_t tmp_521 = 0.5 * p_affine_13_0 * tmp_66 + 0.5 * p_affine_13_1 * tmp_69 + 0.5 * p_affine_13_2 * tmp_72;
      real_t a_0_0   = tmp_119 * ( -tmp_105 * tmp_75 - tmp_118 * ( -tmp_107 - tmp_108 - tmp_110 - tmp_111 - tmp_112 - tmp_113 -
                                                                 tmp_115 - tmp_116 - tmp_117 + 1 ) ) +
                     tmp_140 * ( -tmp_126 * tmp_75 - tmp_139 * ( -tmp_128 - tmp_129 - tmp_131 - tmp_132 - tmp_133 - tmp_134 -
                                                                 tmp_136 - tmp_137 - tmp_138 + 1 ) ) +
                     tmp_161 * ( -tmp_147 * tmp_75 - tmp_160 * ( -tmp_149 - tmp_150 - tmp_152 - tmp_153 - tmp_154 - tmp_155 -
                                                                 tmp_157 - tmp_158 - tmp_159 + 1 ) ) +
                     tmp_182 * ( -tmp_168 * tmp_75 - tmp_181 * ( -tmp_170 - tmp_171 - tmp_173 - tmp_174 - tmp_175 - tmp_176 -
                                                                 tmp_178 - tmp_179 - tmp_180 + 1 ) ) +
                     tmp_203 * ( -tmp_189 * tmp_75 - tmp_202 * ( -tmp_191 - tmp_192 - tmp_194 - tmp_195 - tmp_196 - tmp_197 -
                                                                 tmp_199 - tmp_200 - tmp_201 + 1 ) ) +
                     tmp_224 * ( -tmp_210 * tmp_75 - tmp_223 * ( -tmp_212 - tmp_213 - tmp_215 - tmp_216 - tmp_217 - tmp_218 -
                                                                 tmp_220 - tmp_221 - tmp_222 + 1 ) ) +
                     tmp_245 * ( -tmp_231 * tmp_75 - tmp_244 * ( -tmp_233 - tmp_234 - tmp_236 - tmp_237 - tmp_238 - tmp_239 -
                                                                 tmp_241 - tmp_242 - tmp_243 + 1 ) ) +
                     tmp_266 * ( -tmp_252 * tmp_75 - tmp_265 * ( -tmp_254 - tmp_255 - tmp_257 - tmp_258 - tmp_259 - tmp_260 -
                                                                 tmp_262 - tmp_263 - tmp_264 + 1 ) ) +
                     tmp_287 * ( -tmp_273 * tmp_75 - tmp_286 * ( -tmp_275 - tmp_276 - tmp_278 - tmp_279 - tmp_280 - tmp_281 -
                                                                 tmp_283 - tmp_284 - tmp_285 + 1 ) ) +
                     tmp_308 * ( -tmp_294 * tmp_75 - tmp_307 * ( -tmp_296 - tmp_297 - tmp_299 - tmp_300 - tmp_301 - tmp_302 -
                                                                 tmp_304 - tmp_305 - tmp_306 + 1 ) ) +
                     tmp_329 * ( -tmp_315 * tmp_75 - tmp_328 * ( -tmp_317 - tmp_318 - tmp_320 - tmp_321 - tmp_322 - tmp_323 -
                                                                 tmp_325 - tmp_326 - tmp_327 + 1 ) ) +
                     tmp_350 * ( -tmp_336 * tmp_75 - tmp_349 * ( -tmp_338 - tmp_339 - tmp_341 - tmp_342 - tmp_343 - tmp_344 -
                                                                 tmp_346 - tmp_347 - tmp_348 + 1 ) ) +
                     tmp_371 * ( -tmp_357 * tmp_75 - tmp_370 * ( -tmp_359 - tmp_360 - tmp_362 - tmp_363 - tmp_364 - tmp_365 -
                                                                 tmp_367 - tmp_368 - tmp_369 + 1 ) ) +
                     tmp_392 * ( -tmp_378 * tmp_75 - tmp_391 * ( -tmp_380 - tmp_381 - tmp_383 - tmp_384 - tmp_385 - tmp_386 -
                                                                 tmp_388 - tmp_389 - tmp_390 + 1 ) ) +
                     tmp_413 * ( -tmp_399 * tmp_75 - tmp_412 * ( -tmp_401 - tmp_402 - tmp_404 - tmp_405 - tmp_406 - tmp_407 -
                                                                 tmp_409 - tmp_410 - tmp_411 + 1 ) ) +
                     tmp_434 * ( -tmp_420 * tmp_75 - tmp_433 * ( -tmp_422 - tmp_423 - tmp_425 - tmp_426 - tmp_427 - tmp_428 -
                                                                 tmp_430 - tmp_431 - tmp_432 + 1 ) ) +
                     tmp_455 * ( -tmp_441 * tmp_75 - tmp_454 * ( -tmp_443 - tmp_444 - tmp_446 - tmp_447 - tmp_448 - tmp_449 -
                                                                 tmp_451 - tmp_452 - tmp_453 + 1 ) ) +
                     tmp_476 * ( -tmp_462 * tmp_75 - tmp_475 * ( -tmp_464 - tmp_465 - tmp_467 - tmp_468 - tmp_469 - tmp_470 -
                                                                 tmp_472 - tmp_473 - tmp_474 + 1 ) ) +
                     tmp_497 * ( -tmp_483 * tmp_75 - tmp_496 * ( -tmp_485 - tmp_486 - tmp_488 - tmp_489 - tmp_490 - tmp_491 -
                                                                 tmp_493 - tmp_494 - tmp_495 + 1 ) ) +
                     tmp_518 * ( -tmp_504 * tmp_75 - tmp_517 * ( -tmp_506 - tmp_507 - tmp_509 - tmp_510 - tmp_511 - tmp_512 -
                                                                 tmp_514 - tmp_515 - tmp_516 + 1 ) ) +
                     tmp_98 * ( -tmp_46 * tmp_75 - tmp_96 * ( -tmp_78 - tmp_79 - tmp_82 - tmp_83 - tmp_84 - tmp_85 - tmp_88 -
                                                              tmp_89 - tmp_90 + 1 ) );
      real_t a_0_1 = tmp_119 * ( -tmp_105 * tmp_519 - tmp_118 * ( tmp_112 + tmp_113 + tmp_117 ) ) +
                     tmp_140 * ( -tmp_126 * tmp_519 - tmp_139 * ( tmp_133 + tmp_134 + tmp_138 ) ) +
                     tmp_161 * ( -tmp_147 * tmp_519 - tmp_160 * ( tmp_154 + tmp_155 + tmp_159 ) ) +
                     tmp_182 * ( -tmp_168 * tmp_519 - tmp_181 * ( tmp_175 + tmp_176 + tmp_180 ) ) +
                     tmp_203 * ( -tmp_189 * tmp_519 - tmp_202 * ( tmp_196 + tmp_197 + tmp_201 ) ) +
                     tmp_224 * ( -tmp_210 * tmp_519 - tmp_223 * ( tmp_217 + tmp_218 + tmp_222 ) ) +
                     tmp_245 * ( -tmp_231 * tmp_519 - tmp_244 * ( tmp_238 + tmp_239 + tmp_243 ) ) +
                     tmp_266 * ( -tmp_252 * tmp_519 - tmp_265 * ( tmp_259 + tmp_260 + tmp_264 ) ) +
                     tmp_287 * ( -tmp_273 * tmp_519 - tmp_286 * ( tmp_280 + tmp_281 + tmp_285 ) ) +
                     tmp_308 * ( -tmp_294 * tmp_519 - tmp_307 * ( tmp_301 + tmp_302 + tmp_306 ) ) +
                     tmp_329 * ( -tmp_315 * tmp_519 - tmp_328 * ( tmp_322 + tmp_323 + tmp_327 ) ) +
                     tmp_350 * ( -tmp_336 * tmp_519 - tmp_349 * ( tmp_343 + tmp_344 + tmp_348 ) ) +
                     tmp_371 * ( -tmp_357 * tmp_519 - tmp_370 * ( tmp_364 + tmp_365 + tmp_369 ) ) +
                     tmp_392 * ( -tmp_378 * tmp_519 - tmp_391 * ( tmp_385 + tmp_386 + tmp_390 ) ) +
                     tmp_413 * ( -tmp_399 * tmp_519 - tmp_412 * ( tmp_406 + tmp_407 + tmp_411 ) ) +
                     tmp_434 * ( -tmp_420 * tmp_519 - tmp_433 * ( tmp_427 + tmp_428 + tmp_432 ) ) +
                     tmp_455 * ( -tmp_441 * tmp_519 - tmp_454 * ( tmp_448 + tmp_449 + tmp_453 ) ) +
                     tmp_476 * ( -tmp_462 * tmp_519 - tmp_475 * ( tmp_469 + tmp_470 + tmp_474 ) ) +
                     tmp_497 * ( -tmp_483 * tmp_519 - tmp_496 * ( tmp_490 + tmp_491 + tmp_495 ) ) +
                     tmp_518 * ( -tmp_504 * tmp_519 - tmp_517 * ( tmp_511 + tmp_512 + tmp_516 ) ) +
                     tmp_98 * ( -tmp_46 * tmp_519 - tmp_96 * ( tmp_84 + tmp_85 + tmp_90 ) );
      real_t a_0_2 = tmp_119 * ( -tmp_105 * tmp_520 - tmp_118 * ( tmp_108 + tmp_111 + tmp_116 ) ) +
                     tmp_140 * ( -tmp_126 * tmp_520 - tmp_139 * ( tmp_129 + tmp_132 + tmp_137 ) ) +
                     tmp_161 * ( -tmp_147 * tmp_520 - tmp_160 * ( tmp_150 + tmp_153 + tmp_158 ) ) +
                     tmp_182 * ( -tmp_168 * tmp_520 - tmp_181 * ( tmp_171 + tmp_174 + tmp_179 ) ) +
                     tmp_203 * ( -tmp_189 * tmp_520 - tmp_202 * ( tmp_192 + tmp_195 + tmp_200 ) ) +
                     tmp_224 * ( -tmp_210 * tmp_520 - tmp_223 * ( tmp_213 + tmp_216 + tmp_221 ) ) +
                     tmp_245 * ( -tmp_231 * tmp_520 - tmp_244 * ( tmp_234 + tmp_237 + tmp_242 ) ) +
                     tmp_266 * ( -tmp_252 * tmp_520 - tmp_265 * ( tmp_255 + tmp_258 + tmp_263 ) ) +
                     tmp_287 * ( -tmp_273 * tmp_520 - tmp_286 * ( tmp_276 + tmp_279 + tmp_284 ) ) +
                     tmp_308 * ( -tmp_294 * tmp_520 - tmp_307 * ( tmp_297 + tmp_300 + tmp_305 ) ) +
                     tmp_329 * ( -tmp_315 * tmp_520 - tmp_328 * ( tmp_318 + tmp_321 + tmp_326 ) ) +
                     tmp_350 * ( -tmp_336 * tmp_520 - tmp_349 * ( tmp_339 + tmp_342 + tmp_347 ) ) +
                     tmp_371 * ( -tmp_357 * tmp_520 - tmp_370 * ( tmp_360 + tmp_363 + tmp_368 ) ) +
                     tmp_392 * ( -tmp_378 * tmp_520 - tmp_391 * ( tmp_381 + tmp_384 + tmp_389 ) ) +
                     tmp_413 * ( -tmp_399 * tmp_520 - tmp_412 * ( tmp_402 + tmp_405 + tmp_410 ) ) +
                     tmp_434 * ( -tmp_420 * tmp_520 - tmp_433 * ( tmp_423 + tmp_426 + tmp_431 ) ) +
                     tmp_455 * ( -tmp_441 * tmp_520 - tmp_454 * ( tmp_444 + tmp_447 + tmp_452 ) ) +
                     tmp_476 * ( -tmp_462 * tmp_520 - tmp_475 * ( tmp_465 + tmp_468 + tmp_473 ) ) +
                     tmp_497 * ( -tmp_483 * tmp_520 - tmp_496 * ( tmp_486 + tmp_489 + tmp_494 ) ) +
                     tmp_518 * ( -tmp_504 * tmp_520 - tmp_517 * ( tmp_507 + tmp_510 + tmp_515 ) ) +
                     tmp_98 * ( -tmp_46 * tmp_520 - tmp_96 * ( tmp_79 + tmp_83 + tmp_89 ) );
      real_t a_0_3 = tmp_119 * ( -tmp_105 * tmp_521 - tmp_118 * ( tmp_107 + tmp_110 + tmp_115 ) ) +
                     tmp_140 * ( -tmp_126 * tmp_521 - tmp_139 * ( tmp_128 + tmp_131 + tmp_136 ) ) +
                     tmp_161 * ( -tmp_147 * tmp_521 - tmp_160 * ( tmp_149 + tmp_152 + tmp_157 ) ) +
                     tmp_182 * ( -tmp_168 * tmp_521 - tmp_181 * ( tmp_170 + tmp_173 + tmp_178 ) ) +
                     tmp_203 * ( -tmp_189 * tmp_521 - tmp_202 * ( tmp_191 + tmp_194 + tmp_199 ) ) +
                     tmp_224 * ( -tmp_210 * tmp_521 - tmp_223 * ( tmp_212 + tmp_215 + tmp_220 ) ) +
                     tmp_245 * ( -tmp_231 * tmp_521 - tmp_244 * ( tmp_233 + tmp_236 + tmp_241 ) ) +
                     tmp_266 * ( -tmp_252 * tmp_521 - tmp_265 * ( tmp_254 + tmp_257 + tmp_262 ) ) +
                     tmp_287 * ( -tmp_273 * tmp_521 - tmp_286 * ( tmp_275 + tmp_278 + tmp_283 ) ) +
                     tmp_308 * ( -tmp_294 * tmp_521 - tmp_307 * ( tmp_296 + tmp_299 + tmp_304 ) ) +
                     tmp_329 * ( -tmp_315 * tmp_521 - tmp_328 * ( tmp_317 + tmp_320 + tmp_325 ) ) +
                     tmp_350 * ( -tmp_336 * tmp_521 - tmp_349 * ( tmp_338 + tmp_341 + tmp_346 ) ) +
                     tmp_371 * ( -tmp_357 * tmp_521 - tmp_370 * ( tmp_359 + tmp_362 + tmp_367 ) ) +
                     tmp_392 * ( -tmp_378 * tmp_521 - tmp_391 * ( tmp_380 + tmp_383 + tmp_388 ) ) +
                     tmp_413 * ( -tmp_399 * tmp_521 - tmp_412 * ( tmp_401 + tmp_404 + tmp_409 ) ) +
                     tmp_434 * ( -tmp_420 * tmp_521 - tmp_433 * ( tmp_422 + tmp_425 + tmp_430 ) ) +
                     tmp_455 * ( -tmp_441 * tmp_521 - tmp_454 * ( tmp_443 + tmp_446 + tmp_451 ) ) +
                     tmp_476 * ( -tmp_462 * tmp_521 - tmp_475 * ( tmp_464 + tmp_467 + tmp_472 ) ) +
                     tmp_497 * ( -tmp_483 * tmp_521 - tmp_496 * ( tmp_485 + tmp_488 + tmp_493 ) ) +
                     tmp_518 * ( -tmp_504 * tmp_521 - tmp_517 * ( tmp_506 + tmp_509 + tmp_514 ) ) +
                     tmp_98 * ( -tmp_46 * tmp_521 - tmp_96 * ( tmp_78 + tmp_82 + tmp_88 ) );
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 0, 3 ) = a_0_3;
   }

   void integrateFacetDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                           const std::vector< Point3D >& coordsFacet,
                                           const Point3D&,
                                           const Point3D&     outwardNormal,
                                           const DGBasisInfo& trialBasis,
                                           const DGBasisInfo& testBasis,
                                           int                trialDegree,
                                           int                testDegree,
                                           MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_1_1 + tmp_0;
      real_t tmp_2  = -p_affine_0_2;
      real_t tmp_3  = p_affine_2_2 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_1 + tmp_0;
      real_t tmp_6  = p_affine_1_2 + tmp_2;
      real_t tmp_7  = tmp_5 * tmp_6;
      real_t tmp_8  = -p_affine_0_0;
      real_t tmp_9  = p_affine_1_0 + tmp_8;
      real_t tmp_10 = p_affine_3_2 + tmp_2;
      real_t tmp_11 = tmp_10 * tmp_5;
      real_t tmp_12 = p_affine_2_0 + tmp_8;
      real_t tmp_13 = p_affine_3_1 + tmp_0;
      real_t tmp_14 = tmp_13 * tmp_6;
      real_t tmp_15 = p_affine_3_0 + tmp_8;
      real_t tmp_16 = tmp_13 * tmp_3;
      real_t tmp_17 = tmp_1 * tmp_10;
      real_t tmp_18 =
          1.0 / ( tmp_11 * tmp_9 + tmp_12 * tmp_14 - tmp_12 * tmp_17 + tmp_15 * tmp_4 - tmp_15 * tmp_7 - tmp_16 * tmp_9 );
      real_t tmp_19 = tmp_18 * ( tmp_4 - tmp_7 );
      real_t tmp_20 = tmp_18 * ( tmp_14 - tmp_17 );
      real_t tmp_21 = tmp_18 * ( tmp_11 - tmp_16 );
      real_t tmp_22 = tmp_18 * ( tmp_12 * tmp_6 - tmp_3 * tmp_9 );
      real_t tmp_23 = tmp_18 * ( tmp_10 * tmp_9 - tmp_15 * tmp_6 );
      real_t tmp_24 = tmp_18 * ( -tmp_10 * tmp_12 + tmp_15 * tmp_3 );
      real_t tmp_25 = tmp_18 * ( -tmp_1 * tmp_12 + tmp_5 * tmp_9 );
      real_t tmp_26 = tmp_18 * ( tmp_1 * tmp_15 - tmp_13 * tmp_9 );
      real_t tmp_27 = tmp_18 * ( tmp_12 * tmp_13 - tmp_15 * tmp_5 );
      real_t tmp_28 = -p_affine_8_0;
      real_t tmp_29 = p_affine_10_0 + tmp_28;
      real_t tmp_30 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_31 = -p_affine_8_1;
      real_t tmp_32 = p_affine_10_1 + tmp_31;
      real_t tmp_33 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_34 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_35 = -p_affine_8_2;
      real_t tmp_36 = p_affine_10_2 + tmp_35;
      real_t tmp_37 =
          1.0 * std::pow( ( std::abs( tmp_29 * tmp_30 - tmp_32 * tmp_33 ) * std::abs( tmp_29 * tmp_30 - tmp_32 * tmp_33 ) ) +
                              ( std::abs( tmp_29 * tmp_34 - tmp_33 * tmp_36 ) * std::abs( tmp_29 * tmp_34 - tmp_33 * tmp_36 ) ) +
                              ( std::abs( tmp_30 * tmp_36 - tmp_32 * tmp_34 ) * std::abs( tmp_30 * tmp_36 - tmp_32 * tmp_34 ) ),
                          1.0 / 2.0 );
      real_t tmp_38 = tmp_37 * ( p_affine_13_0 * ( -tmp_19 - tmp_20 - tmp_21 ) + p_affine_13_1 * ( -tmp_22 - tmp_23 - tmp_24 ) +
                                 p_affine_13_2 * ( -tmp_25 - tmp_26 - tmp_27 ) );
      real_t tmp_39 = p_affine_9_2 + tmp_35;
      real_t tmp_40 = p_affine_8_2 + tmp_2;
      real_t tmp_41 = 0.93718850182767688 * tmp_36 + 0.031405749086161582 * tmp_39 + tmp_40;
      real_t tmp_42 = p_affine_9_1 + tmp_31;
      real_t tmp_43 = p_affine_8_1 + tmp_0;
      real_t tmp_44 = 0.93718850182767688 * tmp_32 + 0.031405749086161582 * tmp_42 + tmp_43;
      real_t tmp_45 = p_affine_9_0 + tmp_28;
      real_t tmp_46 = p_affine_8_0 + tmp_8;
      real_t tmp_47 = 0.93718850182767688 * tmp_29 + 0.031405749086161582 * tmp_45 + tmp_46;
      real_t tmp_48 = 0.0068572537431980923 * tmp_10 * ( tmp_19 * tmp_47 + tmp_22 * tmp_44 + tmp_25 * tmp_41 - 1.0 / 4.0 ) +
                      0.0068572537431980923 * tmp_3 * ( tmp_20 * tmp_47 + tmp_23 * tmp_44 + tmp_26 * tmp_41 - 1.0 / 4.0 ) +
                      0.0068572537431980923 * tmp_6 * ( tmp_21 * tmp_47 + tmp_24 * tmp_44 + tmp_27 * tmp_41 - 1.0 / 4.0 );
      real_t tmp_49 = 0.60796128279561268 * tmp_36 + 0.19601935860219369 * tmp_39 + tmp_40;
      real_t tmp_50 = 0.60796128279561268 * tmp_32 + 0.19601935860219369 * tmp_42 + tmp_43;
      real_t tmp_51 = 0.60796128279561268 * tmp_29 + 0.19601935860219369 * tmp_45 + tmp_46;
      real_t tmp_52 = 0.037198804536718075 * tmp_10 * ( tmp_19 * tmp_51 + tmp_22 * tmp_50 + tmp_25 * tmp_49 - 1.0 / 4.0 ) +
                      0.037198804536718075 * tmp_3 * ( tmp_20 * tmp_51 + tmp_23 * tmp_50 + tmp_26 * tmp_49 - 1.0 / 4.0 ) +
                      0.037198804536718075 * tmp_6 * ( tmp_21 * tmp_51 + tmp_24 * tmp_50 + tmp_27 * tmp_49 - 1.0 / 4.0 );
      real_t tmp_53 = 0.039308471900058539 * tmp_36 + 0.37605877282253791 * tmp_39 + tmp_40;
      real_t tmp_54 = 0.039308471900058539 * tmp_32 + 0.37605877282253791 * tmp_42 + tmp_43;
      real_t tmp_55 = 0.039308471900058539 * tmp_29 + 0.37605877282253791 * tmp_45 + tmp_46;
      real_t tmp_56 = 0.020848748529055869 * tmp_10 * ( tmp_19 * tmp_55 + tmp_22 * tmp_54 + tmp_25 * tmp_53 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_3 * ( tmp_20 * tmp_55 + tmp_23 * tmp_54 + tmp_26 * tmp_53 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_6 * ( tmp_21 * tmp_55 + tmp_24 * tmp_54 + tmp_27 * tmp_53 - 1.0 / 4.0 );
      real_t tmp_57 = 0.1711304259088916 * tmp_36 + 0.78764240869137092 * tmp_39 + tmp_40;
      real_t tmp_58 = 0.1711304259088916 * tmp_32 + 0.78764240869137092 * tmp_42 + tmp_43;
      real_t tmp_59 = 0.1711304259088916 * tmp_29 + 0.78764240869137092 * tmp_45 + tmp_46;
      real_t tmp_60 = 0.019202922745021479 * tmp_10 * ( tmp_19 * tmp_59 + tmp_22 * tmp_58 + tmp_25 * tmp_57 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_3 * ( tmp_20 * tmp_59 + tmp_23 * tmp_58 + tmp_26 * tmp_57 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_6 * ( tmp_21 * tmp_59 + tmp_24 * tmp_58 + tmp_27 * tmp_57 - 1.0 / 4.0 );
      real_t tmp_61 = 0.37605877282253791 * tmp_36 + 0.58463275527740355 * tmp_39 + tmp_40;
      real_t tmp_62 = 0.37605877282253791 * tmp_32 + 0.58463275527740355 * tmp_42 + tmp_43;
      real_t tmp_63 = 0.37605877282253791 * tmp_29 + 0.58463275527740355 * tmp_45 + tmp_46;
      real_t tmp_64 = 0.020848748529055869 * tmp_10 * ( tmp_19 * tmp_63 + tmp_22 * tmp_62 + tmp_25 * tmp_61 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_3 * ( tmp_20 * tmp_63 + tmp_23 * tmp_62 + tmp_26 * tmp_61 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_6 * ( tmp_21 * tmp_63 + tmp_24 * tmp_62 + tmp_27 * tmp_61 - 1.0 / 4.0 );
      real_t tmp_65 = 0.78764240869137092 * tmp_36 + 0.041227165399737475 * tmp_39 + tmp_40;
      real_t tmp_66 = 0.78764240869137092 * tmp_32 + 0.041227165399737475 * tmp_42 + tmp_43;
      real_t tmp_67 = 0.78764240869137092 * tmp_29 + 0.041227165399737475 * tmp_45 + tmp_46;
      real_t tmp_68 = 0.019202922745021479 * tmp_10 * ( tmp_19 * tmp_67 + tmp_22 * tmp_66 + tmp_25 * tmp_65 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_3 * ( tmp_20 * tmp_67 + tmp_23 * tmp_66 + tmp_26 * tmp_65 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_6 * ( tmp_21 * tmp_67 + tmp_24 * tmp_66 + tmp_27 * tmp_65 - 1.0 / 4.0 );
      real_t tmp_69 = 0.58463275527740355 * tmp_36 + 0.039308471900058539 * tmp_39 + tmp_40;
      real_t tmp_70 = 0.58463275527740355 * tmp_32 + 0.039308471900058539 * tmp_42 + tmp_43;
      real_t tmp_71 = 0.58463275527740355 * tmp_29 + 0.039308471900058539 * tmp_45 + tmp_46;
      real_t tmp_72 = 0.020848748529055869 * tmp_10 * ( tmp_19 * tmp_71 + tmp_22 * tmp_70 + tmp_25 * tmp_69 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_3 * ( tmp_20 * tmp_71 + tmp_23 * tmp_70 + tmp_26 * tmp_69 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_6 * ( tmp_21 * tmp_71 + tmp_24 * tmp_70 + tmp_27 * tmp_69 - 1.0 / 4.0 );
      real_t tmp_73 = 0.041227165399737475 * tmp_36 + 0.78764240869137092 * tmp_39 + tmp_40;
      real_t tmp_74 = 0.041227165399737475 * tmp_32 + 0.78764240869137092 * tmp_42 + tmp_43;
      real_t tmp_75 = 0.041227165399737475 * tmp_29 + 0.78764240869137092 * tmp_45 + tmp_46;
      real_t tmp_76 = 0.019202922745021479 * tmp_10 * ( tmp_19 * tmp_75 + tmp_22 * tmp_74 + tmp_25 * tmp_73 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_3 * ( tmp_20 * tmp_75 + tmp_23 * tmp_74 + tmp_26 * tmp_73 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_6 * ( tmp_21 * tmp_75 + tmp_24 * tmp_74 + tmp_27 * tmp_73 - 1.0 / 4.0 );
      real_t tmp_77 = 0.039308471900058539 * tmp_36 + 0.58463275527740355 * tmp_39 + tmp_40;
      real_t tmp_78 = 0.039308471900058539 * tmp_32 + 0.58463275527740355 * tmp_42 + tmp_43;
      real_t tmp_79 = 0.039308471900058539 * tmp_29 + 0.58463275527740355 * tmp_45 + tmp_46;
      real_t tmp_80 = 0.020848748529055869 * tmp_10 * ( tmp_19 * tmp_79 + tmp_22 * tmp_78 + tmp_25 * tmp_77 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_3 * ( tmp_20 * tmp_79 + tmp_23 * tmp_78 + tmp_26 * tmp_77 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_6 * ( tmp_21 * tmp_79 + tmp_24 * tmp_78 + tmp_27 * tmp_77 - 1.0 / 4.0 );
      real_t tmp_81 = 0.78764240869137092 * tmp_36 + 0.1711304259088916 * tmp_39 + tmp_40;
      real_t tmp_82 = 0.78764240869137092 * tmp_32 + 0.1711304259088916 * tmp_42 + tmp_43;
      real_t tmp_83 = 0.78764240869137092 * tmp_29 + 0.1711304259088916 * tmp_45 + tmp_46;
      real_t tmp_84 = 0.019202922745021479 * tmp_10 * ( tmp_19 * tmp_83 + tmp_22 * tmp_82 + tmp_25 * tmp_81 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_3 * ( tmp_20 * tmp_83 + tmp_23 * tmp_82 + tmp_26 * tmp_81 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_6 * ( tmp_21 * tmp_83 + tmp_24 * tmp_82 + tmp_27 * tmp_81 - 1.0 / 4.0 );
      real_t tmp_85 = 0.58463275527740355 * tmp_36 + 0.37605877282253791 * tmp_39 + tmp_40;
      real_t tmp_86 = 0.58463275527740355 * tmp_32 + 0.37605877282253791 * tmp_42 + tmp_43;
      real_t tmp_87 = 0.58463275527740355 * tmp_29 + 0.37605877282253791 * tmp_45 + tmp_46;
      real_t tmp_88 = 0.020848748529055869 * tmp_10 * ( tmp_19 * tmp_87 + tmp_22 * tmp_86 + tmp_25 * tmp_85 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_3 * ( tmp_20 * tmp_87 + tmp_23 * tmp_86 + tmp_26 * tmp_85 - 1.0 / 4.0 ) +
                      0.020848748529055869 * tmp_6 * ( tmp_21 * tmp_87 + tmp_24 * tmp_86 + tmp_27 * tmp_85 - 1.0 / 4.0 );
      real_t tmp_89 = 0.1711304259088916 * tmp_36 + 0.041227165399737475 * tmp_39 + tmp_40;
      real_t tmp_90 = 0.1711304259088916 * tmp_32 + 0.041227165399737475 * tmp_42 + tmp_43;
      real_t tmp_91 = 0.1711304259088916 * tmp_29 + 0.041227165399737475 * tmp_45 + tmp_46;
      real_t tmp_92 = 0.019202922745021479 * tmp_10 * ( tmp_19 * tmp_91 + tmp_22 * tmp_90 + tmp_25 * tmp_89 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_3 * ( tmp_20 * tmp_91 + tmp_23 * tmp_90 + tmp_26 * tmp_89 - 1.0 / 4.0 ) +
                      0.019202922745021479 * tmp_6 * ( tmp_21 * tmp_91 + tmp_24 * tmp_90 + tmp_27 * tmp_89 - 1.0 / 4.0 );
      real_t tmp_93 = 0.19107600050469298 * tmp_36 + 0.40446199974765351 * tmp_39 + tmp_40;
      real_t tmp_94 = 0.19107600050469298 * tmp_32 + 0.40446199974765351 * tmp_42 + tmp_43;
      real_t tmp_95 = 0.19107600050469298 * tmp_29 + 0.40446199974765351 * tmp_45 + tmp_46;
      real_t tmp_96 = 0.042507265838595799 * tmp_10 * ( tmp_19 * tmp_95 + tmp_22 * tmp_94 + tmp_25 * tmp_93 - 1.0 / 4.0 ) +
                      0.042507265838595799 * tmp_3 * ( tmp_20 * tmp_95 + tmp_23 * tmp_94 + tmp_26 * tmp_93 - 1.0 / 4.0 ) +
                      0.042507265838595799 * tmp_6 * ( tmp_21 * tmp_95 + tmp_24 * tmp_94 + tmp_27 * tmp_93 - 1.0 / 4.0 );
      real_t tmp_97  = 0.37605877282253791 * tmp_36 + 0.039308471900058539 * tmp_39 + tmp_40;
      real_t tmp_98  = 0.37605877282253791 * tmp_32 + 0.039308471900058539 * tmp_42 + tmp_43;
      real_t tmp_99  = 0.37605877282253791 * tmp_29 + 0.039308471900058539 * tmp_45 + tmp_46;
      real_t tmp_100 = 0.020848748529055869 * tmp_10 * ( tmp_19 * tmp_99 + tmp_22 * tmp_98 + tmp_25 * tmp_97 - 1.0 / 4.0 ) +
                       0.020848748529055869 * tmp_3 * ( tmp_20 * tmp_99 + tmp_23 * tmp_98 + tmp_26 * tmp_97 - 1.0 / 4.0 ) +
                       0.020848748529055869 * tmp_6 * ( tmp_21 * tmp_99 + tmp_24 * tmp_98 + tmp_27 * tmp_97 - 1.0 / 4.0 );
      real_t tmp_101 = 0.031405749086161582 * tmp_36 + 0.93718850182767688 * tmp_39 + tmp_40;
      real_t tmp_102 = 0.031405749086161582 * tmp_32 + 0.93718850182767688 * tmp_42 + tmp_43;
      real_t tmp_103 = 0.031405749086161582 * tmp_29 + 0.93718850182767688 * tmp_45 + tmp_46;
      real_t tmp_104 = 0.0068572537431980923 * tmp_10 * ( tmp_101 * tmp_25 + tmp_102 * tmp_22 + tmp_103 * tmp_19 - 1.0 / 4.0 ) +
                       0.0068572537431980923 * tmp_3 * ( tmp_101 * tmp_26 + tmp_102 * tmp_23 + tmp_103 * tmp_20 - 1.0 / 4.0 ) +
                       0.0068572537431980923 * tmp_6 * ( tmp_101 * tmp_27 + tmp_102 * tmp_24 + tmp_103 * tmp_21 - 1.0 / 4.0 );
      real_t tmp_105 = 0.19601935860219369 * tmp_36 + 0.60796128279561268 * tmp_39 + tmp_40;
      real_t tmp_106 = 0.19601935860219369 * tmp_32 + 0.60796128279561268 * tmp_42 + tmp_43;
      real_t tmp_107 = 0.19601935860219369 * tmp_29 + 0.60796128279561268 * tmp_45 + tmp_46;
      real_t tmp_108 = 0.037198804536718075 * tmp_10 * ( tmp_105 * tmp_25 + tmp_106 * tmp_22 + tmp_107 * tmp_19 - 1.0 / 4.0 ) +
                       0.037198804536718075 * tmp_3 * ( tmp_105 * tmp_26 + tmp_106 * tmp_23 + tmp_107 * tmp_20 - 1.0 / 4.0 ) +
                       0.037198804536718075 * tmp_6 * ( tmp_105 * tmp_27 + tmp_106 * tmp_24 + tmp_107 * tmp_21 - 1.0 / 4.0 );
      real_t tmp_109 = 0.40446199974765351 * tmp_36 + 0.19107600050469298 * tmp_39 + tmp_40;
      real_t tmp_110 = 0.40446199974765351 * tmp_32 + 0.19107600050469298 * tmp_42 + tmp_43;
      real_t tmp_111 = 0.40446199974765351 * tmp_29 + 0.19107600050469298 * tmp_45 + tmp_46;
      real_t tmp_112 = 0.042507265838595799 * tmp_10 * ( tmp_109 * tmp_25 + tmp_110 * tmp_22 + tmp_111 * tmp_19 - 1.0 / 4.0 ) +
                       0.042507265838595799 * tmp_3 * ( tmp_109 * tmp_26 + tmp_110 * tmp_23 + tmp_111 * tmp_20 - 1.0 / 4.0 ) +
                       0.042507265838595799 * tmp_6 * ( tmp_109 * tmp_27 + tmp_110 * tmp_24 + tmp_111 * tmp_21 - 1.0 / 4.0 );
      real_t tmp_113 = 0.031405749086161582 * tmp_36 + 0.031405749086161582 * tmp_39 + tmp_40;
      real_t tmp_114 = 0.031405749086161582 * tmp_32 + 0.031405749086161582 * tmp_42 + tmp_43;
      real_t tmp_115 = 0.031405749086161582 * tmp_29 + 0.031405749086161582 * tmp_45 + tmp_46;
      real_t tmp_116 = 0.0068572537431980923 * tmp_10 * ( tmp_113 * tmp_25 + tmp_114 * tmp_22 + tmp_115 * tmp_19 - 1.0 / 4.0 ) +
                       0.0068572537431980923 * tmp_3 * ( tmp_113 * tmp_26 + tmp_114 * tmp_23 + tmp_115 * tmp_20 - 1.0 / 4.0 ) +
                       0.0068572537431980923 * tmp_6 * ( tmp_113 * tmp_27 + tmp_114 * tmp_24 + tmp_115 * tmp_21 - 1.0 / 4.0 );
      real_t tmp_117 = 0.19601935860219369 * tmp_36 + 0.19601935860219369 * tmp_39 + tmp_40;
      real_t tmp_118 = 0.19601935860219369 * tmp_32 + 0.19601935860219369 * tmp_42 + tmp_43;
      real_t tmp_119 = 0.19601935860219369 * tmp_29 + 0.19601935860219369 * tmp_45 + tmp_46;
      real_t tmp_120 = 0.037198804536718075 * tmp_10 * ( tmp_117 * tmp_25 + tmp_118 * tmp_22 + tmp_119 * tmp_19 - 1.0 / 4.0 ) +
                       0.037198804536718075 * tmp_3 * ( tmp_117 * tmp_26 + tmp_118 * tmp_23 + tmp_119 * tmp_20 - 1.0 / 4.0 ) +
                       0.037198804536718075 * tmp_6 * ( tmp_117 * tmp_27 + tmp_118 * tmp_24 + tmp_119 * tmp_21 - 1.0 / 4.0 );
      real_t tmp_121 = 0.40446199974765351 * tmp_36 + 0.40446199974765351 * tmp_39 + tmp_40;
      real_t tmp_122 = 0.40446199974765351 * tmp_32 + 0.40446199974765351 * tmp_42 + tmp_43;
      real_t tmp_123 = 0.40446199974765351 * tmp_29 + 0.40446199974765351 * tmp_45 + tmp_46;
      real_t tmp_124 = 0.042507265838595799 * tmp_10 * ( tmp_121 * tmp_25 + tmp_122 * tmp_22 + tmp_123 * tmp_19 - 1.0 / 4.0 ) +
                       0.042507265838595799 * tmp_3 * ( tmp_121 * tmp_26 + tmp_122 * tmp_23 + tmp_123 * tmp_20 - 1.0 / 4.0 ) +
                       0.042507265838595799 * tmp_6 * ( tmp_121 * tmp_27 + tmp_122 * tmp_24 + tmp_123 * tmp_21 - 1.0 / 4.0 );
      real_t tmp_125 = 0.041227165399737475 * tmp_36 + 0.1711304259088916 * tmp_39 + tmp_40;
      real_t tmp_126 = 0.041227165399737475 * tmp_32 + 0.1711304259088916 * tmp_42 + tmp_43;
      real_t tmp_127 = 0.041227165399737475 * tmp_29 + 0.1711304259088916 * tmp_45 + tmp_46;
      real_t tmp_128 = 0.019202922745021479 * tmp_10 * ( tmp_125 * tmp_25 + tmp_126 * tmp_22 + tmp_127 * tmp_19 - 1.0 / 4.0 ) +
                       0.019202922745021479 * tmp_3 * ( tmp_125 * tmp_26 + tmp_126 * tmp_23 + tmp_127 * tmp_20 - 1.0 / 4.0 ) +
                       0.019202922745021479 * tmp_6 * ( tmp_125 * tmp_27 + tmp_126 * tmp_24 + tmp_127 * tmp_21 - 1.0 / 4.0 );
      real_t tmp_129 = tmp_37 * ( p_affine_13_0 * tmp_21 + p_affine_13_1 * tmp_24 + p_affine_13_2 * tmp_27 );
      real_t tmp_130 = tmp_37 * ( p_affine_13_0 * tmp_20 + p_affine_13_1 * tmp_23 + p_affine_13_2 * tmp_26 );
      real_t tmp_131 = tmp_37 * ( p_affine_13_0 * tmp_19 + p_affine_13_1 * tmp_22 + p_affine_13_2 * tmp_25 );
      real_t a_0_0   = -tmp_100 * tmp_38 - tmp_104 * tmp_38 - tmp_108 * tmp_38 - tmp_112 * tmp_38 - tmp_116 * tmp_38 -
                     tmp_120 * tmp_38 - tmp_124 * tmp_38 - tmp_128 * tmp_38 - tmp_38 * tmp_48 - tmp_38 * tmp_52 -
                     tmp_38 * tmp_56 - tmp_38 * tmp_60 - tmp_38 * tmp_64 - tmp_38 * tmp_68 - tmp_38 * tmp_72 - tmp_38 * tmp_76 -
                     tmp_38 * tmp_80 - tmp_38 * tmp_84 - tmp_38 * tmp_88 - tmp_38 * tmp_92 - tmp_38 * tmp_96;
      real_t a_0_1 = -tmp_100 * tmp_129 - tmp_104 * tmp_129 - tmp_108 * tmp_129 - tmp_112 * tmp_129 - tmp_116 * tmp_129 -
                     tmp_120 * tmp_129 - tmp_124 * tmp_129 - tmp_128 * tmp_129 - tmp_129 * tmp_48 - tmp_129 * tmp_52 -
                     tmp_129 * tmp_56 - tmp_129 * tmp_60 - tmp_129 * tmp_64 - tmp_129 * tmp_68 - tmp_129 * tmp_72 -
                     tmp_129 * tmp_76 - tmp_129 * tmp_80 - tmp_129 * tmp_84 - tmp_129 * tmp_88 - tmp_129 * tmp_92 -
                     tmp_129 * tmp_96;
      real_t a_0_2 = -tmp_100 * tmp_130 - tmp_104 * tmp_130 - tmp_108 * tmp_130 - tmp_112 * tmp_130 - tmp_116 * tmp_130 -
                     tmp_120 * tmp_130 - tmp_124 * tmp_130 - tmp_128 * tmp_130 - tmp_130 * tmp_48 - tmp_130 * tmp_52 -
                     tmp_130 * tmp_56 - tmp_130 * tmp_60 - tmp_130 * tmp_64 - tmp_130 * tmp_68 - tmp_130 * tmp_72 -
                     tmp_130 * tmp_76 - tmp_130 * tmp_80 - tmp_130 * tmp_84 - tmp_130 * tmp_88 - tmp_130 * tmp_92 -
                     tmp_130 * tmp_96;
      real_t a_0_3 = -tmp_100 * tmp_131 - tmp_104 * tmp_131 - tmp_108 * tmp_131 - tmp_112 * tmp_131 - tmp_116 * tmp_131 -
                     tmp_120 * tmp_131 - tmp_124 * tmp_131 - tmp_128 * tmp_131 - tmp_131 * tmp_48 - tmp_131 * tmp_52 -
                     tmp_131 * tmp_56 - tmp_131 * tmp_60 - tmp_131 * tmp_64 - tmp_131 * tmp_68 - tmp_131 * tmp_72 -
                     tmp_131 * tmp_76 - tmp_131 * tmp_80 - tmp_131 * tmp_84 - tmp_131 * tmp_88 - tmp_131 * tmp_92 -
                     tmp_131 * tmp_96;
      elMat( 0, 0 ) = a_0_0;
      elMat( 0, 1 ) = a_0_1;
      elMat( 0, 2 ) = a_0_2;
      elMat( 0, 3 ) = a_0_3;
   }
};

class EGIIPGVectorLaplaceFormEE : public hyteg::dg::DGForm
{
 protected:
   void integrateVolume2D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override
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

   virtual void integrateFacetInner2D( const std::vector< Point3D >& coordsElement,
                                       const std::vector< Point3D >& coordsFacet,
                                       const Point3D&                oppositeVertex,
                                       const Point3D&                outwardNormal,
                                       const DGBasisInfo&            trialBasis,
                                       const DGBasisInfo&            testBasis,
                                       int                           trialDegree,
                                       int                           testDegree,
                                       MatrixXr&                     elMat ) const override
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
      real_t tmp_13 = tmp_11 * ( 0.046910077030668018 * tmp_1 + tmp_12 );
      real_t tmp_14 = p_affine_6_0 + tmp_3;
      real_t tmp_15 = tmp_11 * ( 0.046910077030668018 * tmp_0 + tmp_14 );
      real_t tmp_16 = tmp_13 * tmp_5 + tmp_15 * tmp_7 - 1.0 / 3.0;
      real_t tmp_17 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_18 = tmp_13 * tmp_4 + tmp_15 * tmp_17 - 1.0 / 3.0;
      real_t tmp_19 = tmp_16 * tmp_4 + tmp_18 * tmp_9;
      real_t tmp_20 = tmp_11 * tmp_8;
      real_t tmp_21 = tmp_11 * tmp_9;
      real_t tmp_22 = tmp_11 * tmp_5;
      real_t tmp_23 =
          0.5 * p_affine_10_0 * ( tmp_17 * tmp_21 + tmp_20 ) + 0.5 * p_affine_10_1 * ( tmp_21 * tmp_4 + tmp_22 * tmp_4 );
      real_t tmp_24 = tmp_10 * tmp_16 + tmp_18 * tmp_7;
      real_t tmp_25 = tmp_11 * tmp_7;
      real_t tmp_26 =
          0.5 * p_affine_10_0 * ( tmp_10 * tmp_25 + tmp_17 * tmp_25 ) + 0.5 * p_affine_10_1 * ( tmp_10 * tmp_22 + tmp_20 );
      real_t tmp_27 = 1.0 / ( tmp_2 );
      real_t tmp_28 = 0.23076534494715845 * tmp_1 + tmp_12;
      real_t tmp_29 = 0.23076534494715845 * tmp_0 + tmp_14;
      real_t tmp_30 = tmp_22 * tmp_28 + tmp_25 * tmp_29 - 1.0 / 3.0;
      real_t tmp_31 = tmp_11 * tmp_4;
      real_t tmp_32 = tmp_11 * tmp_17;
      real_t tmp_33 = tmp_28 * tmp_31 + tmp_29 * tmp_32 - 1.0 / 3.0;
      real_t tmp_34 = tmp_30 * tmp_4 + tmp_33 * tmp_9;
      real_t tmp_35 = tmp_10 * tmp_30 + tmp_33 * tmp_7;
      real_t tmp_36 = 0.5 * tmp_1 + tmp_12;
      real_t tmp_37 = 0.5 * tmp_0 + tmp_14;
      real_t tmp_38 = tmp_22 * tmp_36 + tmp_25 * tmp_37 - 1.0 / 3.0;
      real_t tmp_39 = tmp_31 * tmp_36 + tmp_32 * tmp_37 - 1.0 / 3.0;
      real_t tmp_40 = tmp_38 * tmp_4 + tmp_39 * tmp_9;
      real_t tmp_41 = tmp_10 * tmp_38 + tmp_39 * tmp_7;
      real_t tmp_42 = 0.7692346550528415 * tmp_1 + tmp_12;
      real_t tmp_43 = 0.7692346550528415 * tmp_0 + tmp_14;
      real_t tmp_44 = tmp_22 * tmp_42 + tmp_25 * tmp_43 - 1.0 / 3.0;
      real_t tmp_45 = tmp_31 * tmp_42 + tmp_32 * tmp_43 - 1.0 / 3.0;
      real_t tmp_46 = tmp_4 * tmp_44 + tmp_45 * tmp_9;
      real_t tmp_47 = tmp_10 * tmp_44 + tmp_45 * tmp_7;
      real_t tmp_48 = 0.95308992296933193 * tmp_1 + tmp_12;
      real_t tmp_49 = 0.95308992296933193 * tmp_0 + tmp_14;
      real_t tmp_50 = tmp_22 * tmp_48 + tmp_25 * tmp_49 - 1.0 / 3.0;
      real_t tmp_51 = tmp_31 * tmp_48 + tmp_32 * tmp_49 - 1.0 / 3.0;
      real_t tmp_52 = tmp_4 * tmp_50 + tmp_51 * tmp_9;
      real_t tmp_53 = tmp_10 * tmp_50 + tmp_51 * tmp_7;
      real_t a_0_0  = 0.11846344252809471 * tmp_2 *
                         ( -tmp_19 * tmp_23 - tmp_24 * tmp_26 + tmp_27 * ( ( tmp_19 * tmp_19 ) + ( tmp_24 * tmp_24 ) ) ) +
                     0.2393143352496831 * tmp_2 *
                         ( -tmp_23 * tmp_34 - tmp_26 * tmp_35 + tmp_27 * ( ( tmp_34 * tmp_34 ) + ( tmp_35 * tmp_35 ) ) ) +
                     0.2844444444444445 * tmp_2 *
                         ( -tmp_23 * tmp_40 - tmp_26 * tmp_41 + tmp_27 * ( ( tmp_40 * tmp_40 ) + ( tmp_41 * tmp_41 ) ) ) +
                     0.2393143352496831 * tmp_2 *
                         ( -tmp_23 * tmp_46 - tmp_26 * tmp_47 + tmp_27 * ( ( tmp_46 * tmp_46 ) + ( tmp_47 * tmp_47 ) ) ) +
                     0.11846344252809471 * tmp_2 *
                         ( -tmp_23 * tmp_52 - tmp_26 * tmp_53 + tmp_27 * ( ( tmp_52 * tmp_52 ) + ( tmp_53 * tmp_53 ) ) );
      elMat( 0, 0 ) = a_0_0;
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >& coordsElementInner,
                                          const std::vector< Point3D >& coordsElementOuter,
                                          const std::vector< Point3D >& coordsFacet,
                                          const Point3D&                oppositeVertexInnerElement,
                                          const Point3D&                oppositeVertexOuterElement,
                                          const Point3D&                outwardNormal,
                                          const DGBasisInfo&            trialBasis,
                                          const DGBasisInfo&            testBasis,
                                          int                           trialDegree,
                                          int                           testDegree,
                                          MatrixXr&                     elMat ) const override
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
      real_t tmp_3  = -p_affine_0_0;
      real_t tmp_4  = p_affine_1_0 + tmp_3;
      real_t tmp_5  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_6  = -p_affine_0_1;
      real_t tmp_7  = p_affine_2_1 + tmp_6;
      real_t tmp_8  = p_affine_2_0 + tmp_3;
      real_t tmp_9  = p_affine_1_1 + tmp_6;
      real_t tmp_10 = 1.0 / ( tmp_4 * tmp_7 - tmp_8 * tmp_9 );
      real_t tmp_11 = p_affine_6_1 + 0.046910077030668018 * tmp_1;
      real_t tmp_12 = tmp_10 * ( tmp_11 + tmp_6 );
      real_t tmp_13 = p_affine_6_0 + 0.046910077030668018 * tmp_0;
      real_t tmp_14 = tmp_10 * ( tmp_13 + tmp_3 );
      real_t tmp_15 = tmp_12 * tmp_5 + tmp_14 * tmp_7 - 1.0 / 3.0;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_12 * tmp_4 + tmp_14 * tmp_16 - 1.0 / 3.0;
      real_t tmp_18 = tmp_15 * tmp_4 + tmp_17 * tmp_8;
      real_t tmp_19 = -p_affine_3_0;
      real_t tmp_20 = p_affine_4_0 + tmp_19;
      real_t tmp_21 = -p_affine_3_1;
      real_t tmp_22 = p_affine_5_1 + tmp_21;
      real_t tmp_23 = tmp_20 * tmp_22;
      real_t tmp_24 = p_affine_5_0 + tmp_19;
      real_t tmp_25 = p_affine_4_1 + tmp_21;
      real_t tmp_26 = 1.0 / ( tmp_23 - tmp_24 * tmp_25 );
      real_t tmp_27 = tmp_23 * tmp_26;
      real_t tmp_28 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_29 = tmp_24 * tmp_26;
      real_t tmp_30 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_31 = tmp_26 * tmp_30;
      real_t tmp_32 =
          0.5 * p_affine_10_0 * ( tmp_27 + tmp_28 * tmp_29 ) + 0.5 * p_affine_10_1 * ( tmp_20 * tmp_29 + tmp_20 * tmp_31 );
      real_t tmp_33 = tmp_15 * tmp_9 + tmp_17 * tmp_7;
      real_t tmp_34 = tmp_22 * tmp_26;
      real_t tmp_35 =
          0.5 * p_affine_10_0 * ( tmp_25 * tmp_34 + tmp_28 * tmp_34 ) + 0.5 * p_affine_10_1 * ( tmp_25 * tmp_31 + tmp_27 );
      real_t tmp_36 = 1.0 / ( tmp_2 );
      real_t tmp_37 = tmp_26 * ( tmp_11 + tmp_21 );
      real_t tmp_38 = tmp_26 * ( tmp_13 + tmp_19 );
      real_t tmp_39 = tmp_22 * tmp_38 + tmp_30 * tmp_37 - 1.0 / 3.0;
      real_t tmp_40 = tmp_20 * tmp_37 + tmp_28 * tmp_38 - 1.0 / 3.0;
      real_t tmp_41 = p_affine_6_1 + 0.23076534494715845 * tmp_1;
      real_t tmp_42 = tmp_10 * ( tmp_41 + tmp_6 );
      real_t tmp_43 = p_affine_6_0 + 0.23076534494715845 * tmp_0;
      real_t tmp_44 = tmp_10 * ( tmp_3 + tmp_43 );
      real_t tmp_45 = tmp_42 * tmp_5 + tmp_44 * tmp_7 - 1.0 / 3.0;
      real_t tmp_46 = tmp_16 * tmp_44 + tmp_4 * tmp_42 - 1.0 / 3.0;
      real_t tmp_47 = tmp_4 * tmp_45 + tmp_46 * tmp_8;
      real_t tmp_48 = tmp_45 * tmp_9 + tmp_46 * tmp_7;
      real_t tmp_49 = tmp_21 + tmp_41;
      real_t tmp_50 = tmp_19 + tmp_43;
      real_t tmp_51 = tmp_31 * tmp_49 + tmp_34 * tmp_50 - 1.0 / 3.0;
      real_t tmp_52 = tmp_20 * tmp_26;
      real_t tmp_53 = tmp_26 * tmp_28;
      real_t tmp_54 = tmp_49 * tmp_52 + tmp_50 * tmp_53 - 1.0 / 3.0;
      real_t tmp_55 = p_affine_6_1 + 0.5 * tmp_1;
      real_t tmp_56 = tmp_10 * ( tmp_55 + tmp_6 );
      real_t tmp_57 = p_affine_6_0 + 0.5 * tmp_0;
      real_t tmp_58 = tmp_10 * ( tmp_3 + tmp_57 );
      real_t tmp_59 = tmp_5 * tmp_56 + tmp_58 * tmp_7 - 1.0 / 3.0;
      real_t tmp_60 = tmp_16 * tmp_58 + tmp_4 * tmp_56 - 1.0 / 3.0;
      real_t tmp_61 = tmp_4 * tmp_59 + tmp_60 * tmp_8;
      real_t tmp_62 = tmp_59 * tmp_9 + tmp_60 * tmp_7;
      real_t tmp_63 = tmp_21 + tmp_55;
      real_t tmp_64 = tmp_19 + tmp_57;
      real_t tmp_65 = tmp_31 * tmp_63 + tmp_34 * tmp_64 - 1.0 / 3.0;
      real_t tmp_66 = tmp_52 * tmp_63 + tmp_53 * tmp_64 - 1.0 / 3.0;
      real_t tmp_67 = p_affine_6_1 + 0.7692346550528415 * tmp_1;
      real_t tmp_68 = tmp_10 * ( tmp_6 + tmp_67 );
      real_t tmp_69 = p_affine_6_0 + 0.7692346550528415 * tmp_0;
      real_t tmp_70 = tmp_10 * ( tmp_3 + tmp_69 );
      real_t tmp_71 = tmp_5 * tmp_68 + tmp_7 * tmp_70 - 1.0 / 3.0;
      real_t tmp_72 = tmp_16 * tmp_70 + tmp_4 * tmp_68 - 1.0 / 3.0;
      real_t tmp_73 = tmp_4 * tmp_71 + tmp_72 * tmp_8;
      real_t tmp_74 = tmp_7 * tmp_72 + tmp_71 * tmp_9;
      real_t tmp_75 = tmp_21 + tmp_67;
      real_t tmp_76 = tmp_19 + tmp_69;
      real_t tmp_77 = tmp_31 * tmp_75 + tmp_34 * tmp_76 - 1.0 / 3.0;
      real_t tmp_78 = tmp_52 * tmp_75 + tmp_53 * tmp_76 - 1.0 / 3.0;
      real_t tmp_79 = p_affine_6_1 + 0.95308992296933193 * tmp_1;
      real_t tmp_80 = tmp_10 * ( tmp_6 + tmp_79 );
      real_t tmp_81 = p_affine_6_0 + 0.95308992296933193 * tmp_0;
      real_t tmp_82 = tmp_10 * ( tmp_3 + tmp_81 );
      real_t tmp_83 = tmp_5 * tmp_80 + tmp_7 * tmp_82 - 1.0 / 3.0;
      real_t tmp_84 = tmp_16 * tmp_82 + tmp_4 * tmp_80 - 1.0 / 3.0;
      real_t tmp_85 = tmp_4 * tmp_83 + tmp_8 * tmp_84;
      real_t tmp_86 = tmp_7 * tmp_84 + tmp_83 * tmp_9;
      real_t tmp_87 = tmp_21 + tmp_79;
      real_t tmp_88 = tmp_19 + tmp_81;
      real_t tmp_89 = tmp_31 * tmp_87 + tmp_34 * tmp_88 - 1.0 / 3.0;
      real_t tmp_90 = tmp_52 * tmp_87 + tmp_53 * tmp_88 - 1.0 / 3.0;
      real_t a_0_0 =
          0.11846344252809471 * tmp_2 *
              ( -tmp_18 * tmp_32 - tmp_33 * tmp_35 -
                tmp_36 * ( tmp_18 * ( tmp_20 * tmp_39 + tmp_24 * tmp_40 ) + tmp_33 * ( tmp_22 * tmp_40 + tmp_25 * tmp_39 ) ) ) +
          0.2393143352496831 * tmp_2 *
              ( -tmp_32 * tmp_47 - tmp_35 * tmp_48 -
                tmp_36 * ( tmp_47 * ( tmp_20 * tmp_51 + tmp_24 * tmp_54 ) + tmp_48 * ( tmp_22 * tmp_54 + tmp_25 * tmp_51 ) ) ) +
          0.2844444444444445 * tmp_2 *
              ( -tmp_32 * tmp_61 - tmp_35 * tmp_62 -
                tmp_36 * ( tmp_61 * ( tmp_20 * tmp_65 + tmp_24 * tmp_66 ) + tmp_62 * ( tmp_22 * tmp_66 + tmp_25 * tmp_65 ) ) ) +
          0.2393143352496831 * tmp_2 *
              ( -tmp_32 * tmp_73 - tmp_35 * tmp_74 -
                tmp_36 * ( tmp_73 * ( tmp_20 * tmp_77 + tmp_24 * tmp_78 ) + tmp_74 * ( tmp_22 * tmp_78 + tmp_25 * tmp_77 ) ) ) +
          0.11846344252809471 * tmp_2 *
              ( -tmp_32 * tmp_85 - tmp_35 * tmp_86 -
                tmp_36 * ( tmp_85 * ( tmp_20 * tmp_89 + tmp_24 * tmp_90 ) + tmp_86 * ( tmp_22 * tmp_90 + tmp_25 * tmp_89 ) ) );
      elMat( 0, 0 ) = a_0_0;
   };

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                                   const std::vector< Point3D >& coordsFacet,
                                                   const Point3D&                oppositeVertex,
                                                   const Point3D&                outwardNormal,
                                                   const DGBasisInfo&            trialBasis,
                                                   const DGBasisInfo&            testBasis,
                                                   int                           trialDegree,
                                                   int                           testDegree,
                                                   MatrixXr&                     elMat ) const override
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
      real_t tmp_2  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3  = -p_affine_0_1;
      real_t tmp_4  = p_affine_2_1 + tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_0;
      real_t tmp_6  = p_affine_1_1 + tmp_3;
      real_t tmp_7  = 1.0 / ( tmp_1 * tmp_4 - tmp_5 * tmp_6 );
      real_t tmp_8  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9  = p_affine_6_1 + tmp_3;
      real_t tmp_10 = tmp_7 * ( 0.046910077030668018 * tmp_8 + tmp_9 );
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_7 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_10 * tmp_2 + tmp_13 * tmp_4 - 1.0 / 3.0;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_1 * tmp_10 + tmp_13 * tmp_15 - 1.0 / 3.0;
      real_t tmp_17 = tmp_7 * ( 0.23076534494715845 * tmp_8 + tmp_9 );
      real_t tmp_18 = tmp_7 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_19 = tmp_17 * tmp_2 + tmp_18 * tmp_4 - 1.0 / 3.0;
      real_t tmp_20 = tmp_1 * tmp_17 + tmp_15 * tmp_18 - 1.0 / 3.0;
      real_t tmp_21 = tmp_7 * ( 0.5 * tmp_8 + tmp_9 );
      real_t tmp_22 = tmp_7 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_23 = tmp_2 * tmp_21 + tmp_22 * tmp_4 - 1.0 / 3.0;
      real_t tmp_24 = tmp_1 * tmp_21 + tmp_15 * tmp_22 - 1.0 / 3.0;
      real_t tmp_25 = tmp_7 * ( 0.7692346550528415 * tmp_8 + tmp_9 );
      real_t tmp_26 = tmp_7 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_27 = tmp_2 * tmp_25 + tmp_26 * tmp_4 - 1.0 / 3.0;
      real_t tmp_28 = tmp_1 * tmp_25 + tmp_15 * tmp_26 - 1.0 / 3.0;
      real_t tmp_29 = tmp_7 * ( 0.95308992296933193 * tmp_8 + tmp_9 );
      real_t tmp_30 = tmp_7 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_31 = tmp_2 * tmp_29 + tmp_30 * tmp_4 - 1.0 / 3.0;
      real_t tmp_32 = tmp_1 * tmp_29 + tmp_15 * tmp_30 - 1.0 / 3.0;
      real_t a_0_0  = 0.11846344252809471 * ( ( tmp_1 * tmp_14 + tmp_16 * tmp_5 ) * ( tmp_1 * tmp_14 + tmp_16 * tmp_5 ) ) +
                     0.2393143352496831 * ( ( tmp_1 * tmp_19 + tmp_20 * tmp_5 ) * ( tmp_1 * tmp_19 + tmp_20 * tmp_5 ) ) +
                     0.2844444444444445 * ( ( tmp_1 * tmp_23 + tmp_24 * tmp_5 ) * ( tmp_1 * tmp_23 + tmp_24 * tmp_5 ) ) +
                     0.2393143352496831 * ( ( tmp_1 * tmp_27 + tmp_28 * tmp_5 ) * ( tmp_1 * tmp_27 + tmp_28 * tmp_5 ) ) +
                     0.11846344252809471 * ( ( tmp_1 * tmp_31 + tmp_32 * tmp_5 ) * ( tmp_1 * tmp_31 + tmp_32 * tmp_5 ) ) +
                     0.11846344252809471 * ( ( tmp_14 * tmp_6 + tmp_16 * tmp_4 ) * ( tmp_14 * tmp_6 + tmp_16 * tmp_4 ) ) +
                     0.2393143352496831 * ( ( tmp_19 * tmp_6 + tmp_20 * tmp_4 ) * ( tmp_19 * tmp_6 + tmp_20 * tmp_4 ) ) +
                     0.2844444444444445 * ( ( tmp_23 * tmp_6 + tmp_24 * tmp_4 ) * ( tmp_23 * tmp_6 + tmp_24 * tmp_4 ) ) +
                     0.2393143352496831 * ( ( tmp_27 * tmp_6 + tmp_28 * tmp_4 ) * ( tmp_27 * tmp_6 + tmp_28 * tmp_4 ) ) +
                     0.11846344252809471 * ( ( tmp_31 * tmp_6 + tmp_32 * tmp_4 ) * ( tmp_31 * tmp_6 + tmp_32 * tmp_4 ) );
      elMat( 0, 0 ) = a_0_0;
   }

   void integrateRHSDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                         const std::vector< Point3D >& coordsFacet,
                                         const Point3D&                oppositeVertex,
                                         const Point3D&                outwardNormal,
                                         const DGBasisInfo&            basis,
                                         int                           degree,
                                         MatrixXr&                     elMat ) const override
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
   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                         const std::vector< Point3D >& coordsFacet,
                                         const Point3D&                oppositeVertex,
                                         const Point3D&                outwardNormal,
                                         const DGBasisInfo&            basis,
                                         int                           degree,
                                         MatrixXr&                     elMat ) const override
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
   void integrateVolume3D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_2_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_3_1 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_3_0 + tmp_0;
      real_t tmp_6  = p_affine_2_1 + tmp_2;
      real_t tmp_7  = tmp_5 * tmp_6;
      real_t tmp_8  = tmp_4 - tmp_7;
      real_t tmp_9  = p_affine_1_0 + tmp_0;
      real_t tmp_10 = -p_affine_0_2;
      real_t tmp_11 = p_affine_3_2 + tmp_10;
      real_t tmp_12 = tmp_6 * tmp_9;
      real_t tmp_13 = p_affine_1_2 + tmp_10;
      real_t tmp_14 = p_affine_2_2 + tmp_10;
      real_t tmp_15 = p_affine_1_1 + tmp_2;
      real_t tmp_16 = tmp_15 * tmp_5;
      real_t tmp_17 = tmp_3 * tmp_9;
      real_t tmp_18 = tmp_1 * tmp_15;
      real_t tmp_19 =
          1.0 / ( tmp_11 * tmp_12 - tmp_11 * tmp_18 + tmp_13 * tmp_4 - tmp_13 * tmp_7 + tmp_14 * tmp_16 - tmp_14 * tmp_17 );
      real_t tmp_20 = tmp_19 * tmp_9;
      real_t tmp_21 = tmp_16 - tmp_17;
      real_t tmp_22 = tmp_1 * tmp_19;
      real_t tmp_23 = tmp_12 - tmp_18;
      real_t tmp_24 = tmp_19 * tmp_5;
      real_t tmp_25 = -tmp_1 * tmp_11 + tmp_14 * tmp_5;
      real_t tmp_26 = tmp_11 * tmp_9 - tmp_13 * tmp_5;
      real_t tmp_27 = tmp_1 * tmp_13 - tmp_14 * tmp_9;
      real_t tmp_28 = tmp_11 * tmp_6 - tmp_14 * tmp_3;
      real_t tmp_29 = -tmp_11 * tmp_15 + tmp_13 * tmp_3;
      real_t tmp_30 = -tmp_13 * tmp_6 + tmp_14 * tmp_15;
      real_t tmp_31 = tmp_15 * tmp_19;
      real_t tmp_32 = tmp_19 * tmp_6;
      real_t tmp_33 = tmp_19 * tmp_3;
      real_t tmp_34 = tmp_13 * tmp_19;
      real_t tmp_35 = tmp_14 * tmp_19;
      real_t tmp_36 = tmp_11 * tmp_19;
      real_t tmp_37 = p_affine_0_0 * p_affine_1_1;
      real_t tmp_38 = p_affine_0_0 * p_affine_1_2;
      real_t tmp_39 = p_affine_2_1 * p_affine_3_2;
      real_t tmp_40 = p_affine_0_1 * p_affine_1_0;
      real_t tmp_41 = p_affine_0_1 * p_affine_1_2;
      real_t tmp_42 = p_affine_2_2 * p_affine_3_0;
      real_t tmp_43 = p_affine_0_2 * p_affine_1_0;
      real_t tmp_44 = p_affine_0_2 * p_affine_1_1;
      real_t tmp_45 = p_affine_2_0 * p_affine_3_1;
      real_t tmp_46 = p_affine_2_2 * p_affine_3_1;
      real_t tmp_47 = p_affine_2_0 * p_affine_3_2;
      real_t tmp_48 = p_affine_2_1 * p_affine_3_0;
      real_t tmp_49 =
          ( ( ( tmp_20 * tmp_25 + tmp_22 * tmp_26 + tmp_24 * tmp_27 ) *
              ( tmp_20 * tmp_25 + tmp_22 * tmp_26 + tmp_24 * tmp_27 ) ) +
            ( ( tmp_20 * tmp_28 + tmp_22 * tmp_29 + tmp_24 * tmp_30 ) *
              ( tmp_20 * tmp_28 + tmp_22 * tmp_29 + tmp_24 * tmp_30 ) ) +
            ( ( tmp_20 * tmp_8 + tmp_21 * tmp_22 + tmp_23 * tmp_24 ) * ( tmp_20 * tmp_8 + tmp_21 * tmp_22 + tmp_23 * tmp_24 ) ) +
            ( ( tmp_21 * tmp_32 + tmp_23 * tmp_33 + tmp_31 * tmp_8 ) * ( tmp_21 * tmp_32 + tmp_23 * tmp_33 + tmp_31 * tmp_8 ) ) +
            ( ( tmp_21 * tmp_35 + tmp_23 * tmp_36 + tmp_34 * tmp_8 ) * ( tmp_21 * tmp_35 + tmp_23 * tmp_36 + tmp_34 * tmp_8 ) ) +
            ( ( tmp_25 * tmp_31 + tmp_26 * tmp_32 + tmp_27 * tmp_33 ) *
              ( tmp_25 * tmp_31 + tmp_26 * tmp_32 + tmp_27 * tmp_33 ) ) +
            ( ( tmp_25 * tmp_34 + tmp_26 * tmp_35 + tmp_27 * tmp_36 ) *
              ( tmp_25 * tmp_34 + tmp_26 * tmp_35 + tmp_27 * tmp_36 ) ) +
            ( ( tmp_28 * tmp_31 + tmp_29 * tmp_32 + tmp_30 * tmp_33 ) *
              ( tmp_28 * tmp_31 + tmp_29 * tmp_32 + tmp_30 * tmp_33 ) ) +
            ( ( tmp_28 * tmp_34 + tmp_29 * tmp_35 + tmp_30 * tmp_36 ) *
              ( tmp_28 * tmp_34 + tmp_29 * tmp_35 + tmp_30 * tmp_36 ) ) ) *
          std::abs( p_affine_0_0 * tmp_39 - p_affine_0_0 * tmp_46 + p_affine_0_1 * tmp_42 - p_affine_0_1 * tmp_47 +
                    p_affine_0_2 * tmp_45 - p_affine_0_2 * tmp_48 - p_affine_1_0 * tmp_39 + p_affine_1_0 * tmp_46 -
                    p_affine_1_1 * tmp_42 + p_affine_1_1 * tmp_47 - p_affine_1_2 * tmp_45 + p_affine_1_2 * tmp_48 +
                    p_affine_2_0 * tmp_41 - p_affine_2_0 * tmp_44 - p_affine_2_1 * tmp_38 + p_affine_2_1 * tmp_43 +
                    p_affine_2_2 * tmp_37 - p_affine_2_2 * tmp_40 - p_affine_3_0 * tmp_41 + p_affine_3_0 * tmp_44 +
                    p_affine_3_1 * tmp_38 - p_affine_3_1 * tmp_43 - p_affine_3_2 * tmp_37 + p_affine_3_2 * tmp_40 );
      real_t a_0_0  = 0.1666666666666668 * tmp_49;
      elMat( 0, 0 ) = a_0_0;
   }

   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                               const std::vector< Point3D >& coordsFacet,
                               const Point3D&,
                               const Point3D&     outwardNormal,
                               const DGBasisInfo& trialBasis,
                               const DGBasisInfo& testBasis,
                               int                trialDegree,
                               int                testDegree,
                               MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = p_affine_2_0 + tmp_0;
      real_t tmp_3  = -p_affine_0_1;
      real_t tmp_4  = p_affine_3_1 + tmp_3;
      real_t tmp_5  = tmp_2 * tmp_4;
      real_t tmp_6  = p_affine_3_0 + tmp_0;
      real_t tmp_7  = p_affine_2_1 + tmp_3;
      real_t tmp_8  = tmp_6 * tmp_7;
      real_t tmp_9  = tmp_5 - tmp_8;
      real_t tmp_10 = -p_affine_0_2;
      real_t tmp_11 = p_affine_3_2 + tmp_10;
      real_t tmp_12 = tmp_11 * tmp_7;
      real_t tmp_13 = p_affine_1_2 + tmp_10;
      real_t tmp_14 = p_affine_1_1 + tmp_3;
      real_t tmp_15 = p_affine_2_2 + tmp_10;
      real_t tmp_16 = tmp_15 * tmp_6;
      real_t tmp_17 = tmp_15 * tmp_4;
      real_t tmp_18 = tmp_11 * tmp_2;
      real_t tmp_19 =
          1.0 / ( tmp_1 * tmp_12 - tmp_1 * tmp_17 + tmp_13 * tmp_5 - tmp_13 * tmp_8 + tmp_14 * tmp_16 - tmp_14 * tmp_18 );
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_10;
      real_t tmp_24 = tmp_19 * ( 0.031405749086161582 * tmp_21 + 0.93718850182767688 * tmp_22 + tmp_23 );
      real_t tmp_25 = tmp_16 - tmp_18;
      real_t tmp_26 = -p_affine_8_1;
      real_t tmp_27 = p_affine_9_1 + tmp_26;
      real_t tmp_28 = p_affine_10_1 + tmp_26;
      real_t tmp_29 = p_affine_8_1 + tmp_3;
      real_t tmp_30 = tmp_19 * ( 0.031405749086161582 * tmp_27 + 0.93718850182767688 * tmp_28 + tmp_29 );
      real_t tmp_31 = tmp_12 - tmp_17;
      real_t tmp_32 = -p_affine_8_0;
      real_t tmp_33 = p_affine_9_0 + tmp_32;
      real_t tmp_34 = p_affine_10_0 + tmp_32;
      real_t tmp_35 = p_affine_8_0 + tmp_0;
      real_t tmp_36 = tmp_19 * ( 0.031405749086161582 * tmp_33 + 0.93718850182767688 * tmp_34 + tmp_35 );
      real_t tmp_37 = tmp_24 * tmp_9 + tmp_25 * tmp_30 + tmp_31 * tmp_36 - 1.0 / 4.0;
      real_t tmp_38 = -tmp_1 * tmp_4 + tmp_14 * tmp_6;
      real_t tmp_39 = tmp_1 * tmp_11 - tmp_13 * tmp_6;
      real_t tmp_40 = -tmp_11 * tmp_14 + tmp_13 * tmp_4;
      real_t tmp_41 = tmp_24 * tmp_38 + tmp_30 * tmp_39 + tmp_36 * tmp_40 - 1.0 / 4.0;
      real_t tmp_42 = tmp_1 * tmp_7 - tmp_14 * tmp_2;
      real_t tmp_43 = -tmp_1 * tmp_15 + tmp_13 * tmp_2;
      real_t tmp_44 = -tmp_13 * tmp_7 + tmp_14 * tmp_15;
      real_t tmp_45 = tmp_24 * tmp_42 + tmp_30 * tmp_43 + tmp_36 * tmp_44 - 1.0 / 4.0;
      real_t tmp_46 = tmp_1 * tmp_37 + tmp_2 * tmp_41 + tmp_45 * tmp_6;
      real_t tmp_47 = tmp_1 * tmp_19;
      real_t tmp_48 = tmp_19 * tmp_2;
      real_t tmp_49 = tmp_19 * tmp_6;
      real_t tmp_50 = 0.5 * p_affine_13_0 * ( tmp_31 * tmp_47 + tmp_40 * tmp_48 + tmp_44 * tmp_49 ) +
                      0.5 * p_affine_13_1 * ( tmp_25 * tmp_47 + tmp_39 * tmp_48 + tmp_43 * tmp_49 ) +
                      0.5 * p_affine_13_2 * ( tmp_38 * tmp_48 + tmp_42 * tmp_49 + tmp_47 * tmp_9 );
      real_t tmp_51 = tmp_14 * tmp_37 + tmp_4 * tmp_45 + tmp_41 * tmp_7;
      real_t tmp_52 = tmp_14 * tmp_19;
      real_t tmp_53 = tmp_19 * tmp_7;
      real_t tmp_54 = tmp_19 * tmp_4;
      real_t tmp_55 = 0.5 * p_affine_13_0 * ( tmp_31 * tmp_52 + tmp_40 * tmp_53 + tmp_44 * tmp_54 ) +
                      0.5 * p_affine_13_1 * ( tmp_25 * tmp_52 + tmp_39 * tmp_53 + tmp_43 * tmp_54 ) +
                      0.5 * p_affine_13_2 * ( tmp_38 * tmp_53 + tmp_42 * tmp_54 + tmp_52 * tmp_9 );
      real_t tmp_56 = tmp_11 * tmp_45 + tmp_13 * tmp_37 + tmp_15 * tmp_41;
      real_t tmp_57 = tmp_13 * tmp_19;
      real_t tmp_58 = tmp_15 * tmp_19;
      real_t tmp_59 = tmp_11 * tmp_19;
      real_t tmp_60 = 0.5 * p_affine_13_0 * ( tmp_31 * tmp_57 + tmp_40 * tmp_58 + tmp_44 * tmp_59 ) +
                      0.5 * p_affine_13_1 * ( tmp_25 * tmp_57 + tmp_39 * tmp_58 + tmp_43 * tmp_59 ) +
                      0.5 * p_affine_13_2 * ( tmp_38 * tmp_58 + tmp_42 * tmp_59 + tmp_57 * tmp_9 );
      real_t tmp_61 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_62 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_63 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_64 = ( std::abs( tmp_22 * tmp_61 - tmp_28 * tmp_63 ) * std::abs( tmp_22 * tmp_61 - tmp_28 * tmp_63 ) ) +
                      ( std::abs( tmp_22 * tmp_62 - tmp_34 * tmp_63 ) * std::abs( tmp_22 * tmp_62 - tmp_34 * tmp_63 ) ) +
                      ( std::abs( tmp_28 * tmp_62 - tmp_34 * tmp_61 ) * std::abs( tmp_28 * tmp_62 - tmp_34 * tmp_61 ) );
      real_t tmp_65  = 1.0 * std::pow( tmp_64, -0.25 );
      real_t tmp_66  = 1.0 * std::pow( tmp_64, 1.0 / 2.0 );
      real_t tmp_67  = tmp_19 * ( 0.19601935860219369 * tmp_21 + 0.60796128279561268 * tmp_22 + tmp_23 );
      real_t tmp_68  = tmp_19 * ( 0.19601935860219369 * tmp_27 + 0.60796128279561268 * tmp_28 + tmp_29 );
      real_t tmp_69  = tmp_19 * ( 0.19601935860219369 * tmp_33 + 0.60796128279561268 * tmp_34 + tmp_35 );
      real_t tmp_70  = tmp_25 * tmp_68 + tmp_31 * tmp_69 + tmp_67 * tmp_9 - 1.0 / 4.0;
      real_t tmp_71  = tmp_38 * tmp_67 + tmp_39 * tmp_68 + tmp_40 * tmp_69 - 1.0 / 4.0;
      real_t tmp_72  = tmp_42 * tmp_67 + tmp_43 * tmp_68 + tmp_44 * tmp_69 - 1.0 / 4.0;
      real_t tmp_73  = tmp_1 * tmp_70 + tmp_2 * tmp_71 + tmp_6 * tmp_72;
      real_t tmp_74  = tmp_14 * tmp_70 + tmp_4 * tmp_72 + tmp_7 * tmp_71;
      real_t tmp_75  = tmp_11 * tmp_72 + tmp_13 * tmp_70 + tmp_15 * tmp_71;
      real_t tmp_76  = tmp_19 * ( 0.37605877282253791 * tmp_21 + 0.039308471900058539 * tmp_22 + tmp_23 );
      real_t tmp_77  = tmp_19 * ( 0.37605877282253791 * tmp_27 + 0.039308471900058539 * tmp_28 + tmp_29 );
      real_t tmp_78  = tmp_19 * ( 0.37605877282253791 * tmp_33 + 0.039308471900058539 * tmp_34 + tmp_35 );
      real_t tmp_79  = tmp_25 * tmp_77 + tmp_31 * tmp_78 + tmp_76 * tmp_9 - 1.0 / 4.0;
      real_t tmp_80  = tmp_38 * tmp_76 + tmp_39 * tmp_77 + tmp_40 * tmp_78 - 1.0 / 4.0;
      real_t tmp_81  = tmp_42 * tmp_76 + tmp_43 * tmp_77 + tmp_44 * tmp_78 - 1.0 / 4.0;
      real_t tmp_82  = tmp_1 * tmp_79 + tmp_2 * tmp_80 + tmp_6 * tmp_81;
      real_t tmp_83  = tmp_14 * tmp_79 + tmp_4 * tmp_81 + tmp_7 * tmp_80;
      real_t tmp_84  = tmp_11 * tmp_81 + tmp_13 * tmp_79 + tmp_15 * tmp_80;
      real_t tmp_85  = tmp_19 * ( 0.78764240869137092 * tmp_21 + 0.1711304259088916 * tmp_22 + tmp_23 );
      real_t tmp_86  = tmp_19 * ( 0.78764240869137092 * tmp_27 + 0.1711304259088916 * tmp_28 + tmp_29 );
      real_t tmp_87  = tmp_19 * ( 0.78764240869137092 * tmp_33 + 0.1711304259088916 * tmp_34 + tmp_35 );
      real_t tmp_88  = tmp_25 * tmp_86 + tmp_31 * tmp_87 + tmp_85 * tmp_9 - 1.0 / 4.0;
      real_t tmp_89  = tmp_38 * tmp_85 + tmp_39 * tmp_86 + tmp_40 * tmp_87 - 1.0 / 4.0;
      real_t tmp_90  = tmp_42 * tmp_85 + tmp_43 * tmp_86 + tmp_44 * tmp_87 - 1.0 / 4.0;
      real_t tmp_91  = tmp_1 * tmp_88 + tmp_2 * tmp_89 + tmp_6 * tmp_90;
      real_t tmp_92  = tmp_14 * tmp_88 + tmp_4 * tmp_90 + tmp_7 * tmp_89;
      real_t tmp_93  = tmp_11 * tmp_90 + tmp_13 * tmp_88 + tmp_15 * tmp_89;
      real_t tmp_94  = tmp_19 * ( 0.58463275527740355 * tmp_21 + 0.37605877282253791 * tmp_22 + tmp_23 );
      real_t tmp_95  = tmp_19 * ( 0.58463275527740355 * tmp_27 + 0.37605877282253791 * tmp_28 + tmp_29 );
      real_t tmp_96  = tmp_19 * ( 0.58463275527740355 * tmp_33 + 0.37605877282253791 * tmp_34 + tmp_35 );
      real_t tmp_97  = tmp_25 * tmp_95 + tmp_31 * tmp_96 + tmp_9 * tmp_94 - 1.0 / 4.0;
      real_t tmp_98  = tmp_38 * tmp_94 + tmp_39 * tmp_95 + tmp_40 * tmp_96 - 1.0 / 4.0;
      real_t tmp_99  = tmp_42 * tmp_94 + tmp_43 * tmp_95 + tmp_44 * tmp_96 - 1.0 / 4.0;
      real_t tmp_100 = tmp_1 * tmp_97 + tmp_2 * tmp_98 + tmp_6 * tmp_99;
      real_t tmp_101 = tmp_14 * tmp_97 + tmp_4 * tmp_99 + tmp_7 * tmp_98;
      real_t tmp_102 = tmp_11 * tmp_99 + tmp_13 * tmp_97 + tmp_15 * tmp_98;
      real_t tmp_103 = tmp_19 * ( 0.041227165399737475 * tmp_21 + 0.78764240869137092 * tmp_22 + tmp_23 );
      real_t tmp_104 = tmp_19 * ( 0.041227165399737475 * tmp_27 + 0.78764240869137092 * tmp_28 + tmp_29 );
      real_t tmp_105 = tmp_19 * ( 0.041227165399737475 * tmp_33 + 0.78764240869137092 * tmp_34 + tmp_35 );
      real_t tmp_106 = tmp_103 * tmp_9 + tmp_104 * tmp_25 + tmp_105 * tmp_31 - 1.0 / 4.0;
      real_t tmp_107 = tmp_103 * tmp_38 + tmp_104 * tmp_39 + tmp_105 * tmp_40 - 1.0 / 4.0;
      real_t tmp_108 = tmp_103 * tmp_42 + tmp_104 * tmp_43 + tmp_105 * tmp_44 - 1.0 / 4.0;
      real_t tmp_109 = tmp_1 * tmp_106 + tmp_107 * tmp_2 + tmp_108 * tmp_6;
      real_t tmp_110 = tmp_106 * tmp_14 + tmp_107 * tmp_7 + tmp_108 * tmp_4;
      real_t tmp_111 = tmp_106 * tmp_13 + tmp_107 * tmp_15 + tmp_108 * tmp_11;
      real_t tmp_112 = tmp_19 * ( 0.039308471900058539 * tmp_21 + 0.58463275527740355 * tmp_22 + tmp_23 );
      real_t tmp_113 = tmp_19 * ( 0.039308471900058539 * tmp_27 + 0.58463275527740355 * tmp_28 + tmp_29 );
      real_t tmp_114 = tmp_19 * ( 0.039308471900058539 * tmp_33 + 0.58463275527740355 * tmp_34 + tmp_35 );
      real_t tmp_115 = tmp_112 * tmp_9 + tmp_113 * tmp_25 + tmp_114 * tmp_31 - 1.0 / 4.0;
      real_t tmp_116 = tmp_112 * tmp_38 + tmp_113 * tmp_39 + tmp_114 * tmp_40 - 1.0 / 4.0;
      real_t tmp_117 = tmp_112 * tmp_42 + tmp_113 * tmp_43 + tmp_114 * tmp_44 - 1.0 / 4.0;
      real_t tmp_118 = tmp_1 * tmp_115 + tmp_116 * tmp_2 + tmp_117 * tmp_6;
      real_t tmp_119 = tmp_115 * tmp_14 + tmp_116 * tmp_7 + tmp_117 * tmp_4;
      real_t tmp_120 = tmp_11 * tmp_117 + tmp_115 * tmp_13 + tmp_116 * tmp_15;
      real_t tmp_121 = tmp_19 * ( 0.78764240869137092 * tmp_21 + 0.041227165399737475 * tmp_22 + tmp_23 );
      real_t tmp_122 = tmp_19 * ( 0.78764240869137092 * tmp_27 + 0.041227165399737475 * tmp_28 + tmp_29 );
      real_t tmp_123 = tmp_19 * ( 0.78764240869137092 * tmp_33 + 0.041227165399737475 * tmp_34 + tmp_35 );
      real_t tmp_124 = tmp_121 * tmp_9 + tmp_122 * tmp_25 + tmp_123 * tmp_31 - 1.0 / 4.0;
      real_t tmp_125 = tmp_121 * tmp_38 + tmp_122 * tmp_39 + tmp_123 * tmp_40 - 1.0 / 4.0;
      real_t tmp_126 = tmp_121 * tmp_42 + tmp_122 * tmp_43 + tmp_123 * tmp_44 - 1.0 / 4.0;
      real_t tmp_127 = tmp_1 * tmp_124 + tmp_125 * tmp_2 + tmp_126 * tmp_6;
      real_t tmp_128 = tmp_124 * tmp_14 + tmp_125 * tmp_7 + tmp_126 * tmp_4;
      real_t tmp_129 = tmp_11 * tmp_126 + tmp_124 * tmp_13 + tmp_125 * tmp_15;
      real_t tmp_130 = tmp_19 * ( 0.58463275527740355 * tmp_21 + 0.039308471900058539 * tmp_22 + tmp_23 );
      real_t tmp_131 = tmp_19 * ( 0.58463275527740355 * tmp_27 + 0.039308471900058539 * tmp_28 + tmp_29 );
      real_t tmp_132 = tmp_19 * ( 0.58463275527740355 * tmp_33 + 0.039308471900058539 * tmp_34 + tmp_35 );
      real_t tmp_133 = tmp_130 * tmp_9 + tmp_131 * tmp_25 + tmp_132 * tmp_31 - 1.0 / 4.0;
      real_t tmp_134 = tmp_130 * tmp_38 + tmp_131 * tmp_39 + tmp_132 * tmp_40 - 1.0 / 4.0;
      real_t tmp_135 = tmp_130 * tmp_42 + tmp_131 * tmp_43 + tmp_132 * tmp_44 - 1.0 / 4.0;
      real_t tmp_136 = tmp_1 * tmp_133 + tmp_134 * tmp_2 + tmp_135 * tmp_6;
      real_t tmp_137 = tmp_133 * tmp_14 + tmp_134 * tmp_7 + tmp_135 * tmp_4;
      real_t tmp_138 = tmp_11 * tmp_135 + tmp_13 * tmp_133 + tmp_134 * tmp_15;
      real_t tmp_139 = tmp_19 * ( 0.1711304259088916 * tmp_21 + 0.78764240869137092 * tmp_22 + tmp_23 );
      real_t tmp_140 = tmp_19 * ( 0.1711304259088916 * tmp_27 + 0.78764240869137092 * tmp_28 + tmp_29 );
      real_t tmp_141 = tmp_19 * ( 0.1711304259088916 * tmp_33 + 0.78764240869137092 * tmp_34 + tmp_35 );
      real_t tmp_142 = tmp_139 * tmp_9 + tmp_140 * tmp_25 + tmp_141 * tmp_31 - 1.0 / 4.0;
      real_t tmp_143 = tmp_139 * tmp_38 + tmp_140 * tmp_39 + tmp_141 * tmp_40 - 1.0 / 4.0;
      real_t tmp_144 = tmp_139 * tmp_42 + tmp_140 * tmp_43 + tmp_141 * tmp_44 - 1.0 / 4.0;
      real_t tmp_145 = tmp_1 * tmp_142 + tmp_143 * tmp_2 + tmp_144 * tmp_6;
      real_t tmp_146 = tmp_14 * tmp_142 + tmp_143 * tmp_7 + tmp_144 * tmp_4;
      real_t tmp_147 = tmp_11 * tmp_144 + tmp_13 * tmp_142 + tmp_143 * tmp_15;
      real_t tmp_148 = tmp_19 * ( 0.37605877282253791 * tmp_21 + 0.58463275527740355 * tmp_22 + tmp_23 );
      real_t tmp_149 = tmp_19 * ( 0.37605877282253791 * tmp_27 + 0.58463275527740355 * tmp_28 + tmp_29 );
      real_t tmp_150 = tmp_19 * ( 0.37605877282253791 * tmp_33 + 0.58463275527740355 * tmp_34 + tmp_35 );
      real_t tmp_151 = tmp_148 * tmp_9 + tmp_149 * tmp_25 + tmp_150 * tmp_31 - 1.0 / 4.0;
      real_t tmp_152 = tmp_148 * tmp_38 + tmp_149 * tmp_39 + tmp_150 * tmp_40 - 1.0 / 4.0;
      real_t tmp_153 = tmp_148 * tmp_42 + tmp_149 * tmp_43 + tmp_150 * tmp_44 - 1.0 / 4.0;
      real_t tmp_154 = tmp_1 * tmp_151 + tmp_152 * tmp_2 + tmp_153 * tmp_6;
      real_t tmp_155 = tmp_14 * tmp_151 + tmp_152 * tmp_7 + tmp_153 * tmp_4;
      real_t tmp_156 = tmp_11 * tmp_153 + tmp_13 * tmp_151 + tmp_15 * tmp_152;
      real_t tmp_157 = tmp_19 * ( 0.041227165399737475 * tmp_21 + 0.1711304259088916 * tmp_22 + tmp_23 );
      real_t tmp_158 = tmp_19 * ( 0.041227165399737475 * tmp_27 + 0.1711304259088916 * tmp_28 + tmp_29 );
      real_t tmp_159 = tmp_19 * ( 0.041227165399737475 * tmp_33 + 0.1711304259088916 * tmp_34 + tmp_35 );
      real_t tmp_160 = tmp_157 * tmp_9 + tmp_158 * tmp_25 + tmp_159 * tmp_31 - 1.0 / 4.0;
      real_t tmp_161 = tmp_157 * tmp_38 + tmp_158 * tmp_39 + tmp_159 * tmp_40 - 1.0 / 4.0;
      real_t tmp_162 = tmp_157 * tmp_42 + tmp_158 * tmp_43 + tmp_159 * tmp_44 - 1.0 / 4.0;
      real_t tmp_163 = tmp_1 * tmp_160 + tmp_161 * tmp_2 + tmp_162 * tmp_6;
      real_t tmp_164 = tmp_14 * tmp_160 + tmp_161 * tmp_7 + tmp_162 * tmp_4;
      real_t tmp_165 = tmp_11 * tmp_162 + tmp_13 * tmp_160 + tmp_15 * tmp_161;
      real_t tmp_166 = tmp_19 * ( 0.40446199974765351 * tmp_21 + 0.19107600050469298 * tmp_22 + tmp_23 );
      real_t tmp_167 = tmp_19 * ( 0.40446199974765351 * tmp_27 + 0.19107600050469298 * tmp_28 + tmp_29 );
      real_t tmp_168 = tmp_19 * ( 0.40446199974765351 * tmp_33 + 0.19107600050469298 * tmp_34 + tmp_35 );
      real_t tmp_169 = tmp_166 * tmp_9 + tmp_167 * tmp_25 + tmp_168 * tmp_31 - 1.0 / 4.0;
      real_t tmp_170 = tmp_166 * tmp_38 + tmp_167 * tmp_39 + tmp_168 * tmp_40 - 1.0 / 4.0;
      real_t tmp_171 = tmp_166 * tmp_42 + tmp_167 * tmp_43 + tmp_168 * tmp_44 - 1.0 / 4.0;
      real_t tmp_172 = tmp_1 * tmp_169 + tmp_170 * tmp_2 + tmp_171 * tmp_6;
      real_t tmp_173 = tmp_14 * tmp_169 + tmp_170 * tmp_7 + tmp_171 * tmp_4;
      real_t tmp_174 = tmp_11 * tmp_171 + tmp_13 * tmp_169 + tmp_15 * tmp_170;
      real_t tmp_175 = tmp_19 * ( 0.039308471900058539 * tmp_21 + 0.37605877282253791 * tmp_22 + tmp_23 );
      real_t tmp_176 = tmp_19 * ( 0.039308471900058539 * tmp_27 + 0.37605877282253791 * tmp_28 + tmp_29 );
      real_t tmp_177 = tmp_19 * ( 0.039308471900058539 * tmp_33 + 0.37605877282253791 * tmp_34 + tmp_35 );
      real_t tmp_178 = tmp_175 * tmp_9 + tmp_176 * tmp_25 + tmp_177 * tmp_31 - 1.0 / 4.0;
      real_t tmp_179 = tmp_175 * tmp_38 + tmp_176 * tmp_39 + tmp_177 * tmp_40 - 1.0 / 4.0;
      real_t tmp_180 = tmp_175 * tmp_42 + tmp_176 * tmp_43 + tmp_177 * tmp_44 - 1.0 / 4.0;
      real_t tmp_181 = tmp_1 * tmp_178 + tmp_179 * tmp_2 + tmp_180 * tmp_6;
      real_t tmp_182 = tmp_14 * tmp_178 + tmp_179 * tmp_7 + tmp_180 * tmp_4;
      real_t tmp_183 = tmp_11 * tmp_180 + tmp_13 * tmp_178 + tmp_15 * tmp_179;
      real_t tmp_184 = tmp_19 * ( 0.93718850182767688 * tmp_21 + 0.031405749086161582 * tmp_22 + tmp_23 );
      real_t tmp_185 = tmp_19 * ( 0.93718850182767688 * tmp_27 + 0.031405749086161582 * tmp_28 + tmp_29 );
      real_t tmp_186 = tmp_19 * ( 0.93718850182767688 * tmp_33 + 0.031405749086161582 * tmp_34 + tmp_35 );
      real_t tmp_187 = tmp_184 * tmp_9 + tmp_185 * tmp_25 + tmp_186 * tmp_31 - 1.0 / 4.0;
      real_t tmp_188 = tmp_184 * tmp_38 + tmp_185 * tmp_39 + tmp_186 * tmp_40 - 1.0 / 4.0;
      real_t tmp_189 = tmp_184 * tmp_42 + tmp_185 * tmp_43 + tmp_186 * tmp_44 - 1.0 / 4.0;
      real_t tmp_190 = tmp_1 * tmp_187 + tmp_188 * tmp_2 + tmp_189 * tmp_6;
      real_t tmp_191 = tmp_14 * tmp_187 + tmp_188 * tmp_7 + tmp_189 * tmp_4;
      real_t tmp_192 = tmp_11 * tmp_189 + tmp_13 * tmp_187 + tmp_15 * tmp_188;
      real_t tmp_193 = tmp_19 * ( 0.60796128279561268 * tmp_21 + 0.19601935860219369 * tmp_22 + tmp_23 );
      real_t tmp_194 = tmp_19 * ( 0.60796128279561268 * tmp_27 + 0.19601935860219369 * tmp_28 + tmp_29 );
      real_t tmp_195 = tmp_19 * ( 0.60796128279561268 * tmp_33 + 0.19601935860219369 * tmp_34 + tmp_35 );
      real_t tmp_196 = tmp_193 * tmp_9 + tmp_194 * tmp_25 + tmp_195 * tmp_31 - 1.0 / 4.0;
      real_t tmp_197 = tmp_193 * tmp_38 + tmp_194 * tmp_39 + tmp_195 * tmp_40 - 1.0 / 4.0;
      real_t tmp_198 = tmp_193 * tmp_42 + tmp_194 * tmp_43 + tmp_195 * tmp_44 - 1.0 / 4.0;
      real_t tmp_199 = tmp_1 * tmp_196 + tmp_197 * tmp_2 + tmp_198 * tmp_6;
      real_t tmp_200 = tmp_14 * tmp_196 + tmp_197 * tmp_7 + tmp_198 * tmp_4;
      real_t tmp_201 = tmp_11 * tmp_198 + tmp_13 * tmp_196 + tmp_15 * tmp_197;
      real_t tmp_202 = tmp_19 * ( 0.19107600050469298 * tmp_21 + 0.40446199974765351 * tmp_22 + tmp_23 );
      real_t tmp_203 = tmp_19 * ( 0.19107600050469298 * tmp_27 + 0.40446199974765351 * tmp_28 + tmp_29 );
      real_t tmp_204 = tmp_19 * ( 0.19107600050469298 * tmp_33 + 0.40446199974765351 * tmp_34 + tmp_35 );
      real_t tmp_205 = tmp_202 * tmp_9 + tmp_203 * tmp_25 + tmp_204 * tmp_31 - 1.0 / 4.0;
      real_t tmp_206 = tmp_202 * tmp_38 + tmp_203 * tmp_39 + tmp_204 * tmp_40 - 1.0 / 4.0;
      real_t tmp_207 = tmp_202 * tmp_42 + tmp_203 * tmp_43 + tmp_204 * tmp_44 - 1.0 / 4.0;
      real_t tmp_208 = tmp_1 * tmp_205 + tmp_2 * tmp_206 + tmp_207 * tmp_6;
      real_t tmp_209 = tmp_14 * tmp_205 + tmp_206 * tmp_7 + tmp_207 * tmp_4;
      real_t tmp_210 = tmp_11 * tmp_207 + tmp_13 * tmp_205 + tmp_15 * tmp_206;
      real_t tmp_211 = tmp_19 * ( 0.031405749086161582 * tmp_21 + 0.031405749086161582 * tmp_22 + tmp_23 );
      real_t tmp_212 = tmp_19 * ( 0.031405749086161582 * tmp_27 + 0.031405749086161582 * tmp_28 + tmp_29 );
      real_t tmp_213 = tmp_19 * ( 0.031405749086161582 * tmp_33 + 0.031405749086161582 * tmp_34 + tmp_35 );
      real_t tmp_214 = tmp_211 * tmp_9 + tmp_212 * tmp_25 + tmp_213 * tmp_31 - 1.0 / 4.0;
      real_t tmp_215 = tmp_211 * tmp_38 + tmp_212 * tmp_39 + tmp_213 * tmp_40 - 1.0 / 4.0;
      real_t tmp_216 = tmp_211 * tmp_42 + tmp_212 * tmp_43 + tmp_213 * tmp_44 - 1.0 / 4.0;
      real_t tmp_217 = tmp_1 * tmp_214 + tmp_2 * tmp_215 + tmp_216 * tmp_6;
      real_t tmp_218 = tmp_14 * tmp_214 + tmp_215 * tmp_7 + tmp_216 * tmp_4;
      real_t tmp_219 = tmp_11 * tmp_216 + tmp_13 * tmp_214 + tmp_15 * tmp_215;
      real_t tmp_220 = tmp_19 * ( 0.19601935860219369 * tmp_21 + 0.19601935860219369 * tmp_22 + tmp_23 );
      real_t tmp_221 = tmp_19 * ( 0.19601935860219369 * tmp_27 + 0.19601935860219369 * tmp_28 + tmp_29 );
      real_t tmp_222 = tmp_19 * ( 0.19601935860219369 * tmp_33 + 0.19601935860219369 * tmp_34 + tmp_35 );
      real_t tmp_223 = tmp_220 * tmp_9 + tmp_221 * tmp_25 + tmp_222 * tmp_31 - 1.0 / 4.0;
      real_t tmp_224 = tmp_220 * tmp_38 + tmp_221 * tmp_39 + tmp_222 * tmp_40 - 1.0 / 4.0;
      real_t tmp_225 = tmp_220 * tmp_42 + tmp_221 * tmp_43 + tmp_222 * tmp_44 - 1.0 / 4.0;
      real_t tmp_226 = tmp_1 * tmp_223 + tmp_2 * tmp_224 + tmp_225 * tmp_6;
      real_t tmp_227 = tmp_14 * tmp_223 + tmp_224 * tmp_7 + tmp_225 * tmp_4;
      real_t tmp_228 = tmp_11 * tmp_225 + tmp_13 * tmp_223 + tmp_15 * tmp_224;
      real_t tmp_229 = tmp_19 * ( 0.40446199974765351 * tmp_21 + 0.40446199974765351 * tmp_22 + tmp_23 );
      real_t tmp_230 = tmp_19 * ( 0.40446199974765351 * tmp_27 + 0.40446199974765351 * tmp_28 + tmp_29 );
      real_t tmp_231 = tmp_19 * ( 0.40446199974765351 * tmp_33 + 0.40446199974765351 * tmp_34 + tmp_35 );
      real_t tmp_232 = tmp_229 * tmp_9 + tmp_230 * tmp_25 + tmp_231 * tmp_31 - 1.0 / 4.0;
      real_t tmp_233 = tmp_229 * tmp_38 + tmp_230 * tmp_39 + tmp_231 * tmp_40 - 1.0 / 4.0;
      real_t tmp_234 = tmp_229 * tmp_42 + tmp_230 * tmp_43 + tmp_231 * tmp_44 - 1.0 / 4.0;
      real_t tmp_235 = tmp_1 * tmp_232 + tmp_2 * tmp_233 + tmp_234 * tmp_6;
      real_t tmp_236 = tmp_14 * tmp_232 + tmp_233 * tmp_7 + tmp_234 * tmp_4;
      real_t tmp_237 = tmp_11 * tmp_234 + tmp_13 * tmp_232 + tmp_15 * tmp_233;
      real_t tmp_238 = tmp_19 * ( 0.1711304259088916 * tmp_21 + 0.041227165399737475 * tmp_22 + tmp_23 );
      real_t tmp_239 = tmp_19 * ( 0.1711304259088916 * tmp_27 + 0.041227165399737475 * tmp_28 + tmp_29 );
      real_t tmp_240 = tmp_19 * ( 0.1711304259088916 * tmp_33 + 0.041227165399737475 * tmp_34 + tmp_35 );
      real_t tmp_241 = tmp_238 * tmp_9 + tmp_239 * tmp_25 + tmp_240 * tmp_31 - 1.0 / 4.0;
      real_t tmp_242 = tmp_238 * tmp_38 + tmp_239 * tmp_39 + tmp_240 * tmp_40 - 1.0 / 4.0;
      real_t tmp_243 = tmp_238 * tmp_42 + tmp_239 * tmp_43 + tmp_240 * tmp_44 - 1.0 / 4.0;
      real_t tmp_244 = tmp_1 * tmp_241 + tmp_2 * tmp_242 + tmp_243 * tmp_6;
      real_t tmp_245 = tmp_14 * tmp_241 + tmp_242 * tmp_7 + tmp_243 * tmp_4;
      real_t tmp_246 = tmp_11 * tmp_243 + tmp_13 * tmp_241 + tmp_15 * tmp_242;
      real_t a_0_0   = 0.020848748529055869 * tmp_66 *
                         ( -tmp_100 * tmp_50 - tmp_101 * tmp_55 - tmp_102 * tmp_60 +
                           tmp_65 * ( ( tmp_100 * tmp_100 ) + ( tmp_101 * tmp_101 ) + ( tmp_102 * tmp_102 ) ) ) +
                     0.019202922745021479 * tmp_66 *
                         ( -tmp_109 * tmp_50 - tmp_110 * tmp_55 - tmp_111 * tmp_60 +
                           tmp_65 * ( ( tmp_109 * tmp_109 ) + ( tmp_110 * tmp_110 ) + ( tmp_111 * tmp_111 ) ) ) +
                     0.020848748529055869 * tmp_66 *
                         ( -tmp_118 * tmp_50 - tmp_119 * tmp_55 - tmp_120 * tmp_60 +
                           tmp_65 * ( ( tmp_118 * tmp_118 ) + ( tmp_119 * tmp_119 ) + ( tmp_120 * tmp_120 ) ) ) +
                     0.019202922745021479 * tmp_66 *
                         ( -tmp_127 * tmp_50 - tmp_128 * tmp_55 - tmp_129 * tmp_60 +
                           tmp_65 * ( ( tmp_127 * tmp_127 ) + ( tmp_128 * tmp_128 ) + ( tmp_129 * tmp_129 ) ) ) +
                     0.020848748529055869 * tmp_66 *
                         ( -tmp_136 * tmp_50 - tmp_137 * tmp_55 - tmp_138 * tmp_60 +
                           tmp_65 * ( ( tmp_136 * tmp_136 ) + ( tmp_137 * tmp_137 ) + ( tmp_138 * tmp_138 ) ) ) +
                     0.019202922745021479 * tmp_66 *
                         ( -tmp_145 * tmp_50 - tmp_146 * tmp_55 - tmp_147 * tmp_60 +
                           tmp_65 * ( ( tmp_145 * tmp_145 ) + ( tmp_146 * tmp_146 ) + ( tmp_147 * tmp_147 ) ) ) +
                     0.020848748529055869 * tmp_66 *
                         ( -tmp_154 * tmp_50 - tmp_155 * tmp_55 - tmp_156 * tmp_60 +
                           tmp_65 * ( ( tmp_154 * tmp_154 ) + ( tmp_155 * tmp_155 ) + ( tmp_156 * tmp_156 ) ) ) +
                     0.019202922745021479 * tmp_66 *
                         ( -tmp_163 * tmp_50 - tmp_164 * tmp_55 - tmp_165 * tmp_60 +
                           tmp_65 * ( ( tmp_163 * tmp_163 ) + ( tmp_164 * tmp_164 ) + ( tmp_165 * tmp_165 ) ) ) +
                     0.042507265838595799 * tmp_66 *
                         ( -tmp_172 * tmp_50 - tmp_173 * tmp_55 - tmp_174 * tmp_60 +
                           tmp_65 * ( ( tmp_172 * tmp_172 ) + ( tmp_173 * tmp_173 ) + ( tmp_174 * tmp_174 ) ) ) +
                     0.020848748529055869 * tmp_66 *
                         ( -tmp_181 * tmp_50 - tmp_182 * tmp_55 - tmp_183 * tmp_60 +
                           tmp_65 * ( ( tmp_181 * tmp_181 ) + ( tmp_182 * tmp_182 ) + ( tmp_183 * tmp_183 ) ) ) +
                     0.0068572537431980923 * tmp_66 *
                         ( -tmp_190 * tmp_50 - tmp_191 * tmp_55 - tmp_192 * tmp_60 +
                           tmp_65 * ( ( tmp_190 * tmp_190 ) + ( tmp_191 * tmp_191 ) + ( tmp_192 * tmp_192 ) ) ) +
                     0.037198804536718075 * tmp_66 *
                         ( -tmp_199 * tmp_50 - tmp_200 * tmp_55 - tmp_201 * tmp_60 +
                           tmp_65 * ( ( tmp_199 * tmp_199 ) + ( tmp_200 * tmp_200 ) + ( tmp_201 * tmp_201 ) ) ) +
                     0.042507265838595799 * tmp_66 *
                         ( -tmp_208 * tmp_50 - tmp_209 * tmp_55 - tmp_210 * tmp_60 +
                           tmp_65 * ( ( tmp_208 * tmp_208 ) + ( tmp_209 * tmp_209 ) + ( tmp_210 * tmp_210 ) ) ) +
                     0.0068572537431980923 * tmp_66 *
                         ( -tmp_217 * tmp_50 - tmp_218 * tmp_55 - tmp_219 * tmp_60 +
                           tmp_65 * ( ( tmp_217 * tmp_217 ) + ( tmp_218 * tmp_218 ) + ( tmp_219 * tmp_219 ) ) ) +
                     0.037198804536718075 * tmp_66 *
                         ( -tmp_226 * tmp_50 - tmp_227 * tmp_55 - tmp_228 * tmp_60 +
                           tmp_65 * ( ( tmp_226 * tmp_226 ) + ( tmp_227 * tmp_227 ) + ( tmp_228 * tmp_228 ) ) ) +
                     0.042507265838595799 * tmp_66 *
                         ( -tmp_235 * tmp_50 - tmp_236 * tmp_55 - tmp_237 * tmp_60 +
                           tmp_65 * ( ( tmp_235 * tmp_235 ) + ( tmp_236 * tmp_236 ) + ( tmp_237 * tmp_237 ) ) ) +
                     0.019202922745021479 * tmp_66 *
                         ( -tmp_244 * tmp_50 - tmp_245 * tmp_55 - tmp_246 * tmp_60 +
                           tmp_65 * ( ( tmp_244 * tmp_244 ) + ( tmp_245 * tmp_245 ) + ( tmp_246 * tmp_246 ) ) ) +
                     0.0068572537431980923 * tmp_66 *
                         ( -tmp_46 * tmp_50 - tmp_51 * tmp_55 - tmp_56 * tmp_60 +
                           tmp_65 * ( ( tmp_46 * tmp_46 ) + ( tmp_51 * tmp_51 ) + ( tmp_56 * tmp_56 ) ) ) +
                     0.037198804536718075 * tmp_66 *
                         ( -tmp_50 * tmp_73 - tmp_55 * tmp_74 - tmp_60 * tmp_75 +
                           tmp_65 * ( ( tmp_73 * tmp_73 ) + ( tmp_74 * tmp_74 ) + ( tmp_75 * tmp_75 ) ) ) +
                     0.020848748529055869 * tmp_66 *
                         ( -tmp_50 * tmp_82 - tmp_55 * tmp_83 - tmp_60 * tmp_84 +
                           tmp_65 * ( ( tmp_82 * tmp_82 ) + ( tmp_83 * tmp_83 ) + ( tmp_84 * tmp_84 ) ) ) +
                     0.019202922745021479 * tmp_66 *
                         ( -tmp_50 * tmp_91 - tmp_55 * tmp_92 - tmp_60 * tmp_93 +
                           tmp_65 * ( ( tmp_91 * tmp_91 ) + ( tmp_92 * tmp_92 ) + ( tmp_93 * tmp_93 ) ) );
      elMat( 0, 0 ) = a_0_0;
   }

   void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                  const std::vector< Point3D >& coordsElementOuter,
                                  const std::vector< Point3D >& coordsFacet,
                                  const Point3D&,
                                  const Point3D&,
                                  const Point3D&     outwardNormal,
                                  const DGBasisInfo& trialBasis,
                                  const DGBasisInfo& testBasis,
                                  int                trialDegree,
                                  int                testDegree,
                                  MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = p_affine_2_0 + tmp_0;
      real_t tmp_3  = -p_affine_0_1;
      real_t tmp_4  = p_affine_3_1 + tmp_3;
      real_t tmp_5  = tmp_2 * tmp_4;
      real_t tmp_6  = p_affine_3_0 + tmp_0;
      real_t tmp_7  = p_affine_2_1 + tmp_3;
      real_t tmp_8  = tmp_6 * tmp_7;
      real_t tmp_9  = tmp_5 - tmp_8;
      real_t tmp_10 = -p_affine_0_2;
      real_t tmp_11 = p_affine_3_2 + tmp_10;
      real_t tmp_12 = tmp_11 * tmp_7;
      real_t tmp_13 = p_affine_1_2 + tmp_10;
      real_t tmp_14 = p_affine_1_1 + tmp_3;
      real_t tmp_15 = p_affine_2_2 + tmp_10;
      real_t tmp_16 = tmp_15 * tmp_6;
      real_t tmp_17 = tmp_15 * tmp_4;
      real_t tmp_18 = tmp_11 * tmp_2;
      real_t tmp_19 =
          1.0 / ( tmp_1 * tmp_12 - tmp_1 * tmp_17 + tmp_13 * tmp_5 - tmp_13 * tmp_8 + tmp_14 * tmp_16 - tmp_14 * tmp_18 );
      real_t tmp_20 = p_affine_8_2 + tmp_10;
      real_t tmp_21 = -p_affine_8_2;
      real_t tmp_22 = p_affine_9_2 + tmp_21;
      real_t tmp_23 = p_affine_10_2 + tmp_21;
      real_t tmp_24 = 0.031405749086161582 * tmp_22 + 0.93718850182767688 * tmp_23;
      real_t tmp_25 = tmp_19 * ( tmp_20 + tmp_24 );
      real_t tmp_26 = tmp_16 - tmp_18;
      real_t tmp_27 = p_affine_8_1 + tmp_3;
      real_t tmp_28 = -p_affine_8_1;
      real_t tmp_29 = p_affine_9_1 + tmp_28;
      real_t tmp_30 = p_affine_10_1 + tmp_28;
      real_t tmp_31 = 0.031405749086161582 * tmp_29 + 0.93718850182767688 * tmp_30;
      real_t tmp_32 = tmp_19 * ( tmp_27 + tmp_31 );
      real_t tmp_33 = tmp_12 - tmp_17;
      real_t tmp_34 = p_affine_8_0 + tmp_0;
      real_t tmp_35 = -p_affine_8_0;
      real_t tmp_36 = p_affine_9_0 + tmp_35;
      real_t tmp_37 = p_affine_10_0 + tmp_35;
      real_t tmp_38 = 0.031405749086161582 * tmp_36 + 0.93718850182767688 * tmp_37;
      real_t tmp_39 = tmp_19 * ( tmp_34 + tmp_38 );
      real_t tmp_40 = tmp_25 * tmp_9 + tmp_26 * tmp_32 + tmp_33 * tmp_39 - 1.0 / 4.0;
      real_t tmp_41 = -tmp_1 * tmp_4 + tmp_14 * tmp_6;
      real_t tmp_42 = tmp_1 * tmp_11 - tmp_13 * tmp_6;
      real_t tmp_43 = -tmp_11 * tmp_14 + tmp_13 * tmp_4;
      real_t tmp_44 = tmp_25 * tmp_41 + tmp_32 * tmp_42 + tmp_39 * tmp_43 - 1.0 / 4.0;
      real_t tmp_45 = tmp_1 * tmp_7 - tmp_14 * tmp_2;
      real_t tmp_46 = -tmp_1 * tmp_15 + tmp_13 * tmp_2;
      real_t tmp_47 = -tmp_13 * tmp_7 + tmp_14 * tmp_15;
      real_t tmp_48 = tmp_25 * tmp_45 + tmp_32 * tmp_46 + tmp_39 * tmp_47 - 1.0 / 4.0;
      real_t tmp_49 = tmp_1 * tmp_40 + tmp_2 * tmp_44 + tmp_48 * tmp_6;
      real_t tmp_50 = -p_affine_4_1;
      real_t tmp_51 = p_affine_6_1 + tmp_50;
      real_t tmp_52 = -p_affine_4_2;
      real_t tmp_53 = p_affine_7_2 + tmp_52;
      real_t tmp_54 = tmp_51 * tmp_53;
      real_t tmp_55 = p_affine_7_1 + tmp_50;
      real_t tmp_56 = p_affine_6_2 + tmp_52;
      real_t tmp_57 = tmp_55 * tmp_56;
      real_t tmp_58 = tmp_54 - tmp_57;
      real_t tmp_59 = -p_affine_4_0;
      real_t tmp_60 = p_affine_5_0 + tmp_59;
      real_t tmp_61 = p_affine_6_0 + tmp_59;
      real_t tmp_62 = p_affine_5_2 + tmp_52;
      real_t tmp_63 = tmp_55 * tmp_62;
      real_t tmp_64 = p_affine_7_0 + tmp_59;
      real_t tmp_65 = p_affine_5_1 + tmp_50;
      real_t tmp_66 = tmp_56 * tmp_65;
      real_t tmp_67 = tmp_53 * tmp_65;
      real_t tmp_68 = tmp_51 * tmp_62;
      real_t tmp_69 =
          1.0 / ( tmp_54 * tmp_60 - tmp_57 * tmp_60 + tmp_61 * tmp_63 - tmp_61 * tmp_67 + tmp_64 * tmp_66 - tmp_64 * tmp_68 );
      real_t tmp_70 = tmp_60 * tmp_69;
      real_t tmp_71 = tmp_63 - tmp_67;
      real_t tmp_72 = tmp_61 * tmp_69;
      real_t tmp_73 = tmp_66 - tmp_68;
      real_t tmp_74 = tmp_64 * tmp_69;
      real_t tmp_75 = -tmp_53 * tmp_61 + tmp_56 * tmp_64;
      real_t tmp_76 = tmp_53 * tmp_60 - tmp_62 * tmp_64;
      real_t tmp_77 = -tmp_56 * tmp_60 + tmp_61 * tmp_62;
      real_t tmp_78 = -tmp_51 * tmp_64 + tmp_55 * tmp_61;
      real_t tmp_79 = -tmp_55 * tmp_60 + tmp_64 * tmp_65;
      real_t tmp_80 = tmp_51 * tmp_60 - tmp_61 * tmp_65;
      real_t tmp_81 = 0.5 * p_affine_13_0 * ( tmp_58 * tmp_70 + tmp_71 * tmp_72 + tmp_73 * tmp_74 ) +
                      0.5 * p_affine_13_1 * ( tmp_70 * tmp_75 + tmp_72 * tmp_76 + tmp_74 * tmp_77 ) +
                      0.5 * p_affine_13_2 * ( tmp_70 * tmp_78 + tmp_72 * tmp_79 + tmp_74 * tmp_80 );
      real_t tmp_82 = tmp_14 * tmp_40 + tmp_4 * tmp_48 + tmp_44 * tmp_7;
      real_t tmp_83 = tmp_65 * tmp_69;
      real_t tmp_84 = tmp_51 * tmp_69;
      real_t tmp_85 = tmp_55 * tmp_69;
      real_t tmp_86 = 0.5 * p_affine_13_0 * ( tmp_58 * tmp_83 + tmp_71 * tmp_84 + tmp_73 * tmp_85 ) +
                      0.5 * p_affine_13_1 * ( tmp_75 * tmp_83 + tmp_76 * tmp_84 + tmp_77 * tmp_85 ) +
                      0.5 * p_affine_13_2 * ( tmp_78 * tmp_83 + tmp_79 * tmp_84 + tmp_80 * tmp_85 );
      real_t tmp_87 = tmp_11 * tmp_48 + tmp_13 * tmp_40 + tmp_15 * tmp_44;
      real_t tmp_88 = tmp_62 * tmp_69;
      real_t tmp_89 = tmp_56 * tmp_69;
      real_t tmp_90 = tmp_53 * tmp_69;
      real_t tmp_91 = 0.5 * p_affine_13_0 * ( tmp_58 * tmp_88 + tmp_71 * tmp_89 + tmp_73 * tmp_90 ) +
                      0.5 * p_affine_13_1 * ( tmp_75 * tmp_88 + tmp_76 * tmp_89 + tmp_77 * tmp_90 ) +
                      0.5 * p_affine_13_2 * ( tmp_78 * tmp_88 + tmp_79 * tmp_89 + tmp_80 * tmp_90 );
      real_t tmp_92  = p_affine_8_2 + tmp_52;
      real_t tmp_93  = tmp_69 * ( tmp_24 + tmp_92 );
      real_t tmp_94  = p_affine_8_1 + tmp_50;
      real_t tmp_95  = tmp_69 * ( tmp_31 + tmp_94 );
      real_t tmp_96  = p_affine_8_0 + tmp_59;
      real_t tmp_97  = tmp_69 * ( tmp_38 + tmp_96 );
      real_t tmp_98  = tmp_58 * tmp_97 + tmp_75 * tmp_95 + tmp_78 * tmp_93 - 1.0 / 4.0;
      real_t tmp_99  = tmp_71 * tmp_97 + tmp_76 * tmp_95 + tmp_79 * tmp_93 - 1.0 / 4.0;
      real_t tmp_100 = tmp_73 * tmp_97 + tmp_77 * tmp_95 + tmp_80 * tmp_93 - 1.0 / 4.0;
      real_t tmp_101 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_102 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_103 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_104 = ( std::abs( tmp_101 * tmp_23 - tmp_103 * tmp_30 ) * std::abs( tmp_101 * tmp_23 - tmp_103 * tmp_30 ) ) +
                       ( std::abs( tmp_101 * tmp_37 - tmp_102 * tmp_30 ) * std::abs( tmp_101 * tmp_37 - tmp_102 * tmp_30 ) ) +
                       ( std::abs( tmp_102 * tmp_23 - tmp_103 * tmp_37 ) * std::abs( tmp_102 * tmp_23 - tmp_103 * tmp_37 ) );
      real_t tmp_105 = 1.0 * std::pow( tmp_104, -0.25 );
      real_t tmp_106 = 1.0 * std::pow( tmp_104, 1.0 / 2.0 );
      real_t tmp_107 = 0.19601935860219369 * tmp_22 + 0.60796128279561268 * tmp_23;
      real_t tmp_108 = tmp_19 * ( tmp_107 + tmp_20 );
      real_t tmp_109 = 0.19601935860219369 * tmp_29 + 0.60796128279561268 * tmp_30;
      real_t tmp_110 = tmp_19 * ( tmp_109 + tmp_27 );
      real_t tmp_111 = 0.19601935860219369 * tmp_36 + 0.60796128279561268 * tmp_37;
      real_t tmp_112 = tmp_19 * ( tmp_111 + tmp_34 );
      real_t tmp_113 = tmp_108 * tmp_9 + tmp_110 * tmp_26 + tmp_112 * tmp_33 - 1.0 / 4.0;
      real_t tmp_114 = tmp_108 * tmp_41 + tmp_110 * tmp_42 + tmp_112 * tmp_43 - 1.0 / 4.0;
      real_t tmp_115 = tmp_108 * tmp_45 + tmp_110 * tmp_46 + tmp_112 * tmp_47 - 1.0 / 4.0;
      real_t tmp_116 = tmp_1 * tmp_113 + tmp_114 * tmp_2 + tmp_115 * tmp_6;
      real_t tmp_117 = tmp_113 * tmp_14 + tmp_114 * tmp_7 + tmp_115 * tmp_4;
      real_t tmp_118 = tmp_11 * tmp_115 + tmp_113 * tmp_13 + tmp_114 * tmp_15;
      real_t tmp_119 = tmp_69 * ( tmp_107 + tmp_92 );
      real_t tmp_120 = tmp_69 * ( tmp_109 + tmp_94 );
      real_t tmp_121 = tmp_69 * ( tmp_111 + tmp_96 );
      real_t tmp_122 = tmp_119 * tmp_78 + tmp_120 * tmp_75 + tmp_121 * tmp_58 - 1.0 / 4.0;
      real_t tmp_123 = tmp_119 * tmp_79 + tmp_120 * tmp_76 + tmp_121 * tmp_71 - 1.0 / 4.0;
      real_t tmp_124 = tmp_119 * tmp_80 + tmp_120 * tmp_77 + tmp_121 * tmp_73 - 1.0 / 4.0;
      real_t tmp_125 = 0.37605877282253791 * tmp_22 + 0.039308471900058539 * tmp_23;
      real_t tmp_126 = tmp_19 * ( tmp_125 + tmp_20 );
      real_t tmp_127 = 0.37605877282253791 * tmp_29 + 0.039308471900058539 * tmp_30;
      real_t tmp_128 = tmp_19 * ( tmp_127 + tmp_27 );
      real_t tmp_129 = 0.37605877282253791 * tmp_36 + 0.039308471900058539 * tmp_37;
      real_t tmp_130 = tmp_19 * ( tmp_129 + tmp_34 );
      real_t tmp_131 = tmp_126 * tmp_9 + tmp_128 * tmp_26 + tmp_130 * tmp_33 - 1.0 / 4.0;
      real_t tmp_132 = tmp_126 * tmp_41 + tmp_128 * tmp_42 + tmp_130 * tmp_43 - 1.0 / 4.0;
      real_t tmp_133 = tmp_126 * tmp_45 + tmp_128 * tmp_46 + tmp_130 * tmp_47 - 1.0 / 4.0;
      real_t tmp_134 = tmp_1 * tmp_131 + tmp_132 * tmp_2 + tmp_133 * tmp_6;
      real_t tmp_135 = tmp_131 * tmp_14 + tmp_132 * tmp_7 + tmp_133 * tmp_4;
      real_t tmp_136 = tmp_11 * tmp_133 + tmp_13 * tmp_131 + tmp_132 * tmp_15;
      real_t tmp_137 = tmp_69 * ( tmp_125 + tmp_92 );
      real_t tmp_138 = tmp_69 * ( tmp_127 + tmp_94 );
      real_t tmp_139 = tmp_69 * ( tmp_129 + tmp_96 );
      real_t tmp_140 = tmp_137 * tmp_78 + tmp_138 * tmp_75 + tmp_139 * tmp_58 - 1.0 / 4.0;
      real_t tmp_141 = tmp_137 * tmp_79 + tmp_138 * tmp_76 + tmp_139 * tmp_71 - 1.0 / 4.0;
      real_t tmp_142 = tmp_137 * tmp_80 + tmp_138 * tmp_77 + tmp_139 * tmp_73 - 1.0 / 4.0;
      real_t tmp_143 = 0.78764240869137092 * tmp_22 + 0.1711304259088916 * tmp_23;
      real_t tmp_144 = tmp_19 * ( tmp_143 + tmp_20 );
      real_t tmp_145 = 0.78764240869137092 * tmp_29 + 0.1711304259088916 * tmp_30;
      real_t tmp_146 = tmp_19 * ( tmp_145 + tmp_27 );
      real_t tmp_147 = 0.78764240869137092 * tmp_36 + 0.1711304259088916 * tmp_37;
      real_t tmp_148 = tmp_19 * ( tmp_147 + tmp_34 );
      real_t tmp_149 = tmp_144 * tmp_9 + tmp_146 * tmp_26 + tmp_148 * tmp_33 - 1.0 / 4.0;
      real_t tmp_150 = tmp_144 * tmp_41 + tmp_146 * tmp_42 + tmp_148 * tmp_43 - 1.0 / 4.0;
      real_t tmp_151 = tmp_144 * tmp_45 + tmp_146 * tmp_46 + tmp_148 * tmp_47 - 1.0 / 4.0;
      real_t tmp_152 = tmp_1 * tmp_149 + tmp_150 * tmp_2 + tmp_151 * tmp_6;
      real_t tmp_153 = tmp_14 * tmp_149 + tmp_150 * tmp_7 + tmp_151 * tmp_4;
      real_t tmp_154 = tmp_11 * tmp_151 + tmp_13 * tmp_149 + tmp_15 * tmp_150;
      real_t tmp_155 = tmp_69 * ( tmp_143 + tmp_92 );
      real_t tmp_156 = tmp_69 * ( tmp_145 + tmp_94 );
      real_t tmp_157 = tmp_69 * ( tmp_147 + tmp_96 );
      real_t tmp_158 = tmp_155 * tmp_78 + tmp_156 * tmp_75 + tmp_157 * tmp_58 - 1.0 / 4.0;
      real_t tmp_159 = tmp_155 * tmp_79 + tmp_156 * tmp_76 + tmp_157 * tmp_71 - 1.0 / 4.0;
      real_t tmp_160 = tmp_155 * tmp_80 + tmp_156 * tmp_77 + tmp_157 * tmp_73 - 1.0 / 4.0;
      real_t tmp_161 = 0.58463275527740355 * tmp_22 + 0.37605877282253791 * tmp_23;
      real_t tmp_162 = tmp_19 * ( tmp_161 + tmp_20 );
      real_t tmp_163 = 0.58463275527740355 * tmp_29 + 0.37605877282253791 * tmp_30;
      real_t tmp_164 = tmp_19 * ( tmp_163 + tmp_27 );
      real_t tmp_165 = 0.58463275527740355 * tmp_36 + 0.37605877282253791 * tmp_37;
      real_t tmp_166 = tmp_19 * ( tmp_165 + tmp_34 );
      real_t tmp_167 = tmp_162 * tmp_9 + tmp_164 * tmp_26 + tmp_166 * tmp_33 - 1.0 / 4.0;
      real_t tmp_168 = tmp_162 * tmp_41 + tmp_164 * tmp_42 + tmp_166 * tmp_43 - 1.0 / 4.0;
      real_t tmp_169 = tmp_162 * tmp_45 + tmp_164 * tmp_46 + tmp_166 * tmp_47 - 1.0 / 4.0;
      real_t tmp_170 = tmp_1 * tmp_167 + tmp_168 * tmp_2 + tmp_169 * tmp_6;
      real_t tmp_171 = tmp_14 * tmp_167 + tmp_168 * tmp_7 + tmp_169 * tmp_4;
      real_t tmp_172 = tmp_11 * tmp_169 + tmp_13 * tmp_167 + tmp_15 * tmp_168;
      real_t tmp_173 = tmp_69 * ( tmp_161 + tmp_92 );
      real_t tmp_174 = tmp_69 * ( tmp_163 + tmp_94 );
      real_t tmp_175 = tmp_69 * ( tmp_165 + tmp_96 );
      real_t tmp_176 = tmp_173 * tmp_78 + tmp_174 * tmp_75 + tmp_175 * tmp_58 - 1.0 / 4.0;
      real_t tmp_177 = tmp_173 * tmp_79 + tmp_174 * tmp_76 + tmp_175 * tmp_71 - 1.0 / 4.0;
      real_t tmp_178 = tmp_173 * tmp_80 + tmp_174 * tmp_77 + tmp_175 * tmp_73 - 1.0 / 4.0;
      real_t tmp_179 = 0.041227165399737475 * tmp_22 + 0.78764240869137092 * tmp_23;
      real_t tmp_180 = tmp_19 * ( tmp_179 + tmp_20 );
      real_t tmp_181 = 0.041227165399737475 * tmp_29 + 0.78764240869137092 * tmp_30;
      real_t tmp_182 = tmp_19 * ( tmp_181 + tmp_27 );
      real_t tmp_183 = 0.041227165399737475 * tmp_36 + 0.78764240869137092 * tmp_37;
      real_t tmp_184 = tmp_19 * ( tmp_183 + tmp_34 );
      real_t tmp_185 = tmp_180 * tmp_9 + tmp_182 * tmp_26 + tmp_184 * tmp_33 - 1.0 / 4.0;
      real_t tmp_186 = tmp_180 * tmp_41 + tmp_182 * tmp_42 + tmp_184 * tmp_43 - 1.0 / 4.0;
      real_t tmp_187 = tmp_180 * tmp_45 + tmp_182 * tmp_46 + tmp_184 * tmp_47 - 1.0 / 4.0;
      real_t tmp_188 = tmp_1 * tmp_185 + tmp_186 * tmp_2 + tmp_187 * tmp_6;
      real_t tmp_189 = tmp_14 * tmp_185 + tmp_186 * tmp_7 + tmp_187 * tmp_4;
      real_t tmp_190 = tmp_11 * tmp_187 + tmp_13 * tmp_185 + tmp_15 * tmp_186;
      real_t tmp_191 = tmp_69 * ( tmp_179 + tmp_92 );
      real_t tmp_192 = tmp_69 * ( tmp_181 + tmp_94 );
      real_t tmp_193 = tmp_69 * ( tmp_183 + tmp_96 );
      real_t tmp_194 = tmp_191 * tmp_78 + tmp_192 * tmp_75 + tmp_193 * tmp_58 - 1.0 / 4.0;
      real_t tmp_195 = tmp_191 * tmp_79 + tmp_192 * tmp_76 + tmp_193 * tmp_71 - 1.0 / 4.0;
      real_t tmp_196 = tmp_191 * tmp_80 + tmp_192 * tmp_77 + tmp_193 * tmp_73 - 1.0 / 4.0;
      real_t tmp_197 = 0.039308471900058539 * tmp_22 + 0.58463275527740355 * tmp_23;
      real_t tmp_198 = tmp_19 * ( tmp_197 + tmp_20 );
      real_t tmp_199 = 0.039308471900058539 * tmp_29 + 0.58463275527740355 * tmp_30;
      real_t tmp_200 = tmp_19 * ( tmp_199 + tmp_27 );
      real_t tmp_201 = 0.039308471900058539 * tmp_36 + 0.58463275527740355 * tmp_37;
      real_t tmp_202 = tmp_19 * ( tmp_201 + tmp_34 );
      real_t tmp_203 = tmp_198 * tmp_9 + tmp_200 * tmp_26 + tmp_202 * tmp_33 - 1.0 / 4.0;
      real_t tmp_204 = tmp_198 * tmp_41 + tmp_200 * tmp_42 + tmp_202 * tmp_43 - 1.0 / 4.0;
      real_t tmp_205 = tmp_198 * tmp_45 + tmp_200 * tmp_46 + tmp_202 * tmp_47 - 1.0 / 4.0;
      real_t tmp_206 = tmp_1 * tmp_203 + tmp_2 * tmp_204 + tmp_205 * tmp_6;
      real_t tmp_207 = tmp_14 * tmp_203 + tmp_204 * tmp_7 + tmp_205 * tmp_4;
      real_t tmp_208 = tmp_11 * tmp_205 + tmp_13 * tmp_203 + tmp_15 * tmp_204;
      real_t tmp_209 = tmp_69 * ( tmp_197 + tmp_92 );
      real_t tmp_210 = tmp_69 * ( tmp_199 + tmp_94 );
      real_t tmp_211 = tmp_69 * ( tmp_201 + tmp_96 );
      real_t tmp_212 = tmp_209 * tmp_78 + tmp_210 * tmp_75 + tmp_211 * tmp_58 - 1.0 / 4.0;
      real_t tmp_213 = tmp_209 * tmp_79 + tmp_210 * tmp_76 + tmp_211 * tmp_71 - 1.0 / 4.0;
      real_t tmp_214 = tmp_209 * tmp_80 + tmp_210 * tmp_77 + tmp_211 * tmp_73 - 1.0 / 4.0;
      real_t tmp_215 = 0.78764240869137092 * tmp_22 + 0.041227165399737475 * tmp_23;
      real_t tmp_216 = tmp_19 * ( tmp_20 + tmp_215 );
      real_t tmp_217 = 0.78764240869137092 * tmp_29 + 0.041227165399737475 * tmp_30;
      real_t tmp_218 = tmp_19 * ( tmp_217 + tmp_27 );
      real_t tmp_219 = 0.78764240869137092 * tmp_36 + 0.041227165399737475 * tmp_37;
      real_t tmp_220 = tmp_19 * ( tmp_219 + tmp_34 );
      real_t tmp_221 = tmp_216 * tmp_9 + tmp_218 * tmp_26 + tmp_220 * tmp_33 - 1.0 / 4.0;
      real_t tmp_222 = tmp_216 * tmp_41 + tmp_218 * tmp_42 + tmp_220 * tmp_43 - 1.0 / 4.0;
      real_t tmp_223 = tmp_216 * tmp_45 + tmp_218 * tmp_46 + tmp_220 * tmp_47 - 1.0 / 4.0;
      real_t tmp_224 = tmp_1 * tmp_221 + tmp_2 * tmp_222 + tmp_223 * tmp_6;
      real_t tmp_225 = tmp_14 * tmp_221 + tmp_222 * tmp_7 + tmp_223 * tmp_4;
      real_t tmp_226 = tmp_11 * tmp_223 + tmp_13 * tmp_221 + tmp_15 * tmp_222;
      real_t tmp_227 = tmp_69 * ( tmp_215 + tmp_92 );
      real_t tmp_228 = tmp_69 * ( tmp_217 + tmp_94 );
      real_t tmp_229 = tmp_69 * ( tmp_219 + tmp_96 );
      real_t tmp_230 = tmp_227 * tmp_78 + tmp_228 * tmp_75 + tmp_229 * tmp_58 - 1.0 / 4.0;
      real_t tmp_231 = tmp_227 * tmp_79 + tmp_228 * tmp_76 + tmp_229 * tmp_71 - 1.0 / 4.0;
      real_t tmp_232 = tmp_227 * tmp_80 + tmp_228 * tmp_77 + tmp_229 * tmp_73 - 1.0 / 4.0;
      real_t tmp_233 = 0.58463275527740355 * tmp_22 + 0.039308471900058539 * tmp_23;
      real_t tmp_234 = tmp_19 * ( tmp_20 + tmp_233 );
      real_t tmp_235 = 0.58463275527740355 * tmp_29 + 0.039308471900058539 * tmp_30;
      real_t tmp_236 = tmp_19 * ( tmp_235 + tmp_27 );
      real_t tmp_237 = 0.58463275527740355 * tmp_36 + 0.039308471900058539 * tmp_37;
      real_t tmp_238 = tmp_19 * ( tmp_237 + tmp_34 );
      real_t tmp_239 = tmp_234 * tmp_9 + tmp_236 * tmp_26 + tmp_238 * tmp_33 - 1.0 / 4.0;
      real_t tmp_240 = tmp_234 * tmp_41 + tmp_236 * tmp_42 + tmp_238 * tmp_43 - 1.0 / 4.0;
      real_t tmp_241 = tmp_234 * tmp_45 + tmp_236 * tmp_46 + tmp_238 * tmp_47 - 1.0 / 4.0;
      real_t tmp_242 = tmp_1 * tmp_239 + tmp_2 * tmp_240 + tmp_241 * tmp_6;
      real_t tmp_243 = tmp_14 * tmp_239 + tmp_240 * tmp_7 + tmp_241 * tmp_4;
      real_t tmp_244 = tmp_11 * tmp_241 + tmp_13 * tmp_239 + tmp_15 * tmp_240;
      real_t tmp_245 = tmp_69 * ( tmp_233 + tmp_92 );
      real_t tmp_246 = tmp_69 * ( tmp_235 + tmp_94 );
      real_t tmp_247 = tmp_69 * ( tmp_237 + tmp_96 );
      real_t tmp_248 = tmp_245 * tmp_78 + tmp_246 * tmp_75 + tmp_247 * tmp_58 - 1.0 / 4.0;
      real_t tmp_249 = tmp_245 * tmp_79 + tmp_246 * tmp_76 + tmp_247 * tmp_71 - 1.0 / 4.0;
      real_t tmp_250 = tmp_245 * tmp_80 + tmp_246 * tmp_77 + tmp_247 * tmp_73 - 1.0 / 4.0;
      real_t tmp_251 = 0.1711304259088916 * tmp_22 + 0.78764240869137092 * tmp_23;
      real_t tmp_252 = tmp_19 * ( tmp_20 + tmp_251 );
      real_t tmp_253 = 0.1711304259088916 * tmp_29 + 0.78764240869137092 * tmp_30;
      real_t tmp_254 = tmp_19 * ( tmp_253 + tmp_27 );
      real_t tmp_255 = 0.1711304259088916 * tmp_36 + 0.78764240869137092 * tmp_37;
      real_t tmp_256 = tmp_19 * ( tmp_255 + tmp_34 );
      real_t tmp_257 = tmp_252 * tmp_9 + tmp_254 * tmp_26 + tmp_256 * tmp_33 - 1.0 / 4.0;
      real_t tmp_258 = tmp_252 * tmp_41 + tmp_254 * tmp_42 + tmp_256 * tmp_43 - 1.0 / 4.0;
      real_t tmp_259 = tmp_252 * tmp_45 + tmp_254 * tmp_46 + tmp_256 * tmp_47 - 1.0 / 4.0;
      real_t tmp_260 = tmp_1 * tmp_257 + tmp_2 * tmp_258 + tmp_259 * tmp_6;
      real_t tmp_261 = tmp_14 * tmp_257 + tmp_258 * tmp_7 + tmp_259 * tmp_4;
      real_t tmp_262 = tmp_11 * tmp_259 + tmp_13 * tmp_257 + tmp_15 * tmp_258;
      real_t tmp_263 = tmp_69 * ( tmp_251 + tmp_92 );
      real_t tmp_264 = tmp_69 * ( tmp_253 + tmp_94 );
      real_t tmp_265 = tmp_69 * ( tmp_255 + tmp_96 );
      real_t tmp_266 = tmp_263 * tmp_78 + tmp_264 * tmp_75 + tmp_265 * tmp_58 - 1.0 / 4.0;
      real_t tmp_267 = tmp_263 * tmp_79 + tmp_264 * tmp_76 + tmp_265 * tmp_71 - 1.0 / 4.0;
      real_t tmp_268 = tmp_263 * tmp_80 + tmp_264 * tmp_77 + tmp_265 * tmp_73 - 1.0 / 4.0;
      real_t tmp_269 = 0.37605877282253791 * tmp_22 + 0.58463275527740355 * tmp_23;
      real_t tmp_270 = tmp_19 * ( tmp_20 + tmp_269 );
      real_t tmp_271 = 0.37605877282253791 * tmp_29 + 0.58463275527740355 * tmp_30;
      real_t tmp_272 = tmp_19 * ( tmp_27 + tmp_271 );
      real_t tmp_273 = 0.37605877282253791 * tmp_36 + 0.58463275527740355 * tmp_37;
      real_t tmp_274 = tmp_19 * ( tmp_273 + tmp_34 );
      real_t tmp_275 = tmp_26 * tmp_272 + tmp_270 * tmp_9 + tmp_274 * tmp_33 - 1.0 / 4.0;
      real_t tmp_276 = tmp_270 * tmp_41 + tmp_272 * tmp_42 + tmp_274 * tmp_43 - 1.0 / 4.0;
      real_t tmp_277 = tmp_270 * tmp_45 + tmp_272 * tmp_46 + tmp_274 * tmp_47 - 1.0 / 4.0;
      real_t tmp_278 = tmp_1 * tmp_275 + tmp_2 * tmp_276 + tmp_277 * tmp_6;
      real_t tmp_279 = tmp_14 * tmp_275 + tmp_276 * tmp_7 + tmp_277 * tmp_4;
      real_t tmp_280 = tmp_11 * tmp_277 + tmp_13 * tmp_275 + tmp_15 * tmp_276;
      real_t tmp_281 = tmp_69 * ( tmp_269 + tmp_92 );
      real_t tmp_282 = tmp_69 * ( tmp_271 + tmp_94 );
      real_t tmp_283 = tmp_69 * ( tmp_273 + tmp_96 );
      real_t tmp_284 = tmp_281 * tmp_78 + tmp_282 * tmp_75 + tmp_283 * tmp_58 - 1.0 / 4.0;
      real_t tmp_285 = tmp_281 * tmp_79 + tmp_282 * tmp_76 + tmp_283 * tmp_71 - 1.0 / 4.0;
      real_t tmp_286 = tmp_281 * tmp_80 + tmp_282 * tmp_77 + tmp_283 * tmp_73 - 1.0 / 4.0;
      real_t tmp_287 = 0.041227165399737475 * tmp_22 + 0.1711304259088916 * tmp_23;
      real_t tmp_288 = tmp_19 * ( tmp_20 + tmp_287 );
      real_t tmp_289 = 0.041227165399737475 * tmp_29 + 0.1711304259088916 * tmp_30;
      real_t tmp_290 = tmp_19 * ( tmp_27 + tmp_289 );
      real_t tmp_291 = 0.041227165399737475 * tmp_36 + 0.1711304259088916 * tmp_37;
      real_t tmp_292 = tmp_19 * ( tmp_291 + tmp_34 );
      real_t tmp_293 = tmp_26 * tmp_290 + tmp_288 * tmp_9 + tmp_292 * tmp_33 - 1.0 / 4.0;
      real_t tmp_294 = tmp_288 * tmp_41 + tmp_290 * tmp_42 + tmp_292 * tmp_43 - 1.0 / 4.0;
      real_t tmp_295 = tmp_288 * tmp_45 + tmp_290 * tmp_46 + tmp_292 * tmp_47 - 1.0 / 4.0;
      real_t tmp_296 = tmp_1 * tmp_293 + tmp_2 * tmp_294 + tmp_295 * tmp_6;
      real_t tmp_297 = tmp_14 * tmp_293 + tmp_294 * tmp_7 + tmp_295 * tmp_4;
      real_t tmp_298 = tmp_11 * tmp_295 + tmp_13 * tmp_293 + tmp_15 * tmp_294;
      real_t tmp_299 = tmp_69 * ( tmp_287 + tmp_92 );
      real_t tmp_300 = tmp_69 * ( tmp_289 + tmp_94 );
      real_t tmp_301 = tmp_69 * ( tmp_291 + tmp_96 );
      real_t tmp_302 = tmp_299 * tmp_78 + tmp_300 * tmp_75 + tmp_301 * tmp_58 - 1.0 / 4.0;
      real_t tmp_303 = tmp_299 * tmp_79 + tmp_300 * tmp_76 + tmp_301 * tmp_71 - 1.0 / 4.0;
      real_t tmp_304 = tmp_299 * tmp_80 + tmp_300 * tmp_77 + tmp_301 * tmp_73 - 1.0 / 4.0;
      real_t tmp_305 = 0.40446199974765351 * tmp_22 + 0.19107600050469298 * tmp_23;
      real_t tmp_306 = tmp_19 * ( tmp_20 + tmp_305 );
      real_t tmp_307 = 0.40446199974765351 * tmp_29 + 0.19107600050469298 * tmp_30;
      real_t tmp_308 = tmp_19 * ( tmp_27 + tmp_307 );
      real_t tmp_309 = 0.40446199974765351 * tmp_36 + 0.19107600050469298 * tmp_37;
      real_t tmp_310 = tmp_19 * ( tmp_309 + tmp_34 );
      real_t tmp_311 = tmp_26 * tmp_308 + tmp_306 * tmp_9 + tmp_310 * tmp_33 - 1.0 / 4.0;
      real_t tmp_312 = tmp_306 * tmp_41 + tmp_308 * tmp_42 + tmp_310 * tmp_43 - 1.0 / 4.0;
      real_t tmp_313 = tmp_306 * tmp_45 + tmp_308 * tmp_46 + tmp_310 * tmp_47 - 1.0 / 4.0;
      real_t tmp_314 = tmp_1 * tmp_311 + tmp_2 * tmp_312 + tmp_313 * tmp_6;
      real_t tmp_315 = tmp_14 * tmp_311 + tmp_312 * tmp_7 + tmp_313 * tmp_4;
      real_t tmp_316 = tmp_11 * tmp_313 + tmp_13 * tmp_311 + tmp_15 * tmp_312;
      real_t tmp_317 = tmp_69 * ( tmp_305 + tmp_92 );
      real_t tmp_318 = tmp_69 * ( tmp_307 + tmp_94 );
      real_t tmp_319 = tmp_69 * ( tmp_309 + tmp_96 );
      real_t tmp_320 = tmp_317 * tmp_78 + tmp_318 * tmp_75 + tmp_319 * tmp_58 - 1.0 / 4.0;
      real_t tmp_321 = tmp_317 * tmp_79 + tmp_318 * tmp_76 + tmp_319 * tmp_71 - 1.0 / 4.0;
      real_t tmp_322 = tmp_317 * tmp_80 + tmp_318 * tmp_77 + tmp_319 * tmp_73 - 1.0 / 4.0;
      real_t tmp_323 = 0.039308471900058539 * tmp_22 + 0.37605877282253791 * tmp_23;
      real_t tmp_324 = tmp_19 * ( tmp_20 + tmp_323 );
      real_t tmp_325 = 0.039308471900058539 * tmp_29 + 0.37605877282253791 * tmp_30;
      real_t tmp_326 = tmp_19 * ( tmp_27 + tmp_325 );
      real_t tmp_327 = 0.039308471900058539 * tmp_36 + 0.37605877282253791 * tmp_37;
      real_t tmp_328 = tmp_19 * ( tmp_327 + tmp_34 );
      real_t tmp_329 = tmp_26 * tmp_326 + tmp_324 * tmp_9 + tmp_328 * tmp_33 - 1.0 / 4.0;
      real_t tmp_330 = tmp_324 * tmp_41 + tmp_326 * tmp_42 + tmp_328 * tmp_43 - 1.0 / 4.0;
      real_t tmp_331 = tmp_324 * tmp_45 + tmp_326 * tmp_46 + tmp_328 * tmp_47 - 1.0 / 4.0;
      real_t tmp_332 = tmp_1 * tmp_329 + tmp_2 * tmp_330 + tmp_331 * tmp_6;
      real_t tmp_333 = tmp_14 * tmp_329 + tmp_330 * tmp_7 + tmp_331 * tmp_4;
      real_t tmp_334 = tmp_11 * tmp_331 + tmp_13 * tmp_329 + tmp_15 * tmp_330;
      real_t tmp_335 = tmp_69 * ( tmp_323 + tmp_92 );
      real_t tmp_336 = tmp_69 * ( tmp_325 + tmp_94 );
      real_t tmp_337 = tmp_69 * ( tmp_327 + tmp_96 );
      real_t tmp_338 = tmp_335 * tmp_78 + tmp_336 * tmp_75 + tmp_337 * tmp_58 - 1.0 / 4.0;
      real_t tmp_339 = tmp_335 * tmp_79 + tmp_336 * tmp_76 + tmp_337 * tmp_71 - 1.0 / 4.0;
      real_t tmp_340 = tmp_335 * tmp_80 + tmp_336 * tmp_77 + tmp_337 * tmp_73 - 1.0 / 4.0;
      real_t tmp_341 = 0.93718850182767688 * tmp_22 + 0.031405749086161582 * tmp_23;
      real_t tmp_342 = tmp_19 * ( tmp_20 + tmp_341 );
      real_t tmp_343 = 0.93718850182767688 * tmp_29 + 0.031405749086161582 * tmp_30;
      real_t tmp_344 = tmp_19 * ( tmp_27 + tmp_343 );
      real_t tmp_345 = 0.93718850182767688 * tmp_36 + 0.031405749086161582 * tmp_37;
      real_t tmp_346 = tmp_19 * ( tmp_34 + tmp_345 );
      real_t tmp_347 = tmp_26 * tmp_344 + tmp_33 * tmp_346 + tmp_342 * tmp_9 - 1.0 / 4.0;
      real_t tmp_348 = tmp_342 * tmp_41 + tmp_344 * tmp_42 + tmp_346 * tmp_43 - 1.0 / 4.0;
      real_t tmp_349 = tmp_342 * tmp_45 + tmp_344 * tmp_46 + tmp_346 * tmp_47 - 1.0 / 4.0;
      real_t tmp_350 = tmp_1 * tmp_347 + tmp_2 * tmp_348 + tmp_349 * tmp_6;
      real_t tmp_351 = tmp_14 * tmp_347 + tmp_348 * tmp_7 + tmp_349 * tmp_4;
      real_t tmp_352 = tmp_11 * tmp_349 + tmp_13 * tmp_347 + tmp_15 * tmp_348;
      real_t tmp_353 = tmp_69 * ( tmp_341 + tmp_92 );
      real_t tmp_354 = tmp_69 * ( tmp_343 + tmp_94 );
      real_t tmp_355 = tmp_69 * ( tmp_345 + tmp_96 );
      real_t tmp_356 = tmp_353 * tmp_78 + tmp_354 * tmp_75 + tmp_355 * tmp_58 - 1.0 / 4.0;
      real_t tmp_357 = tmp_353 * tmp_79 + tmp_354 * tmp_76 + tmp_355 * tmp_71 - 1.0 / 4.0;
      real_t tmp_358 = tmp_353 * tmp_80 + tmp_354 * tmp_77 + tmp_355 * tmp_73 - 1.0 / 4.0;
      real_t tmp_359 = 0.60796128279561268 * tmp_22 + 0.19601935860219369 * tmp_23;
      real_t tmp_360 = tmp_19 * ( tmp_20 + tmp_359 );
      real_t tmp_361 = 0.60796128279561268 * tmp_29 + 0.19601935860219369 * tmp_30;
      real_t tmp_362 = tmp_19 * ( tmp_27 + tmp_361 );
      real_t tmp_363 = 0.60796128279561268 * tmp_36 + 0.19601935860219369 * tmp_37;
      real_t tmp_364 = tmp_19 * ( tmp_34 + tmp_363 );
      real_t tmp_365 = tmp_26 * tmp_362 + tmp_33 * tmp_364 + tmp_360 * tmp_9 - 1.0 / 4.0;
      real_t tmp_366 = tmp_360 * tmp_41 + tmp_362 * tmp_42 + tmp_364 * tmp_43 - 1.0 / 4.0;
      real_t tmp_367 = tmp_360 * tmp_45 + tmp_362 * tmp_46 + tmp_364 * tmp_47 - 1.0 / 4.0;
      real_t tmp_368 = tmp_1 * tmp_365 + tmp_2 * tmp_366 + tmp_367 * tmp_6;
      real_t tmp_369 = tmp_14 * tmp_365 + tmp_366 * tmp_7 + tmp_367 * tmp_4;
      real_t tmp_370 = tmp_11 * tmp_367 + tmp_13 * tmp_365 + tmp_15 * tmp_366;
      real_t tmp_371 = tmp_69 * ( tmp_359 + tmp_92 );
      real_t tmp_372 = tmp_69 * ( tmp_361 + tmp_94 );
      real_t tmp_373 = tmp_69 * ( tmp_363 + tmp_96 );
      real_t tmp_374 = tmp_371 * tmp_78 + tmp_372 * tmp_75 + tmp_373 * tmp_58 - 1.0 / 4.0;
      real_t tmp_375 = tmp_371 * tmp_79 + tmp_372 * tmp_76 + tmp_373 * tmp_71 - 1.0 / 4.0;
      real_t tmp_376 = tmp_371 * tmp_80 + tmp_372 * tmp_77 + tmp_373 * tmp_73 - 1.0 / 4.0;
      real_t tmp_377 = 0.19107600050469298 * tmp_22 + 0.40446199974765351 * tmp_23;
      real_t tmp_378 = tmp_19 * ( tmp_20 + tmp_377 );
      real_t tmp_379 = 0.19107600050469298 * tmp_29 + 0.40446199974765351 * tmp_30;
      real_t tmp_380 = tmp_19 * ( tmp_27 + tmp_379 );
      real_t tmp_381 = 0.19107600050469298 * tmp_36 + 0.40446199974765351 * tmp_37;
      real_t tmp_382 = tmp_19 * ( tmp_34 + tmp_381 );
      real_t tmp_383 = tmp_26 * tmp_380 + tmp_33 * tmp_382 + tmp_378 * tmp_9 - 1.0 / 4.0;
      real_t tmp_384 = tmp_378 * tmp_41 + tmp_380 * tmp_42 + tmp_382 * tmp_43 - 1.0 / 4.0;
      real_t tmp_385 = tmp_378 * tmp_45 + tmp_380 * tmp_46 + tmp_382 * tmp_47 - 1.0 / 4.0;
      real_t tmp_386 = tmp_1 * tmp_383 + tmp_2 * tmp_384 + tmp_385 * tmp_6;
      real_t tmp_387 = tmp_14 * tmp_383 + tmp_384 * tmp_7 + tmp_385 * tmp_4;
      real_t tmp_388 = tmp_11 * tmp_385 + tmp_13 * tmp_383 + tmp_15 * tmp_384;
      real_t tmp_389 = tmp_69 * ( tmp_377 + tmp_92 );
      real_t tmp_390 = tmp_69 * ( tmp_379 + tmp_94 );
      real_t tmp_391 = tmp_69 * ( tmp_381 + tmp_96 );
      real_t tmp_392 = tmp_389 * tmp_78 + tmp_390 * tmp_75 + tmp_391 * tmp_58 - 1.0 / 4.0;
      real_t tmp_393 = tmp_389 * tmp_79 + tmp_390 * tmp_76 + tmp_391 * tmp_71 - 1.0 / 4.0;
      real_t tmp_394 = tmp_389 * tmp_80 + tmp_390 * tmp_77 + tmp_391 * tmp_73 - 1.0 / 4.0;
      real_t tmp_395 = 0.031405749086161582 * tmp_22 + 0.031405749086161582 * tmp_23;
      real_t tmp_396 = tmp_19 * ( tmp_20 + tmp_395 );
      real_t tmp_397 = 0.031405749086161582 * tmp_29 + 0.031405749086161582 * tmp_30;
      real_t tmp_398 = tmp_19 * ( tmp_27 + tmp_397 );
      real_t tmp_399 = 0.031405749086161582 * tmp_36 + 0.031405749086161582 * tmp_37;
      real_t tmp_400 = tmp_19 * ( tmp_34 + tmp_399 );
      real_t tmp_401 = tmp_26 * tmp_398 + tmp_33 * tmp_400 + tmp_396 * tmp_9 - 1.0 / 4.0;
      real_t tmp_402 = tmp_396 * tmp_41 + tmp_398 * tmp_42 + tmp_400 * tmp_43 - 1.0 / 4.0;
      real_t tmp_403 = tmp_396 * tmp_45 + tmp_398 * tmp_46 + tmp_400 * tmp_47 - 1.0 / 4.0;
      real_t tmp_404 = tmp_1 * tmp_401 + tmp_2 * tmp_402 + tmp_403 * tmp_6;
      real_t tmp_405 = tmp_14 * tmp_401 + tmp_4 * tmp_403 + tmp_402 * tmp_7;
      real_t tmp_406 = tmp_11 * tmp_403 + tmp_13 * tmp_401 + tmp_15 * tmp_402;
      real_t tmp_407 = tmp_69 * ( tmp_395 + tmp_92 );
      real_t tmp_408 = tmp_69 * ( tmp_397 + tmp_94 );
      real_t tmp_409 = tmp_69 * ( tmp_399 + tmp_96 );
      real_t tmp_410 = tmp_407 * tmp_78 + tmp_408 * tmp_75 + tmp_409 * tmp_58 - 1.0 / 4.0;
      real_t tmp_411 = tmp_407 * tmp_79 + tmp_408 * tmp_76 + tmp_409 * tmp_71 - 1.0 / 4.0;
      real_t tmp_412 = tmp_407 * tmp_80 + tmp_408 * tmp_77 + tmp_409 * tmp_73 - 1.0 / 4.0;
      real_t tmp_413 = 0.19601935860219369 * tmp_22 + 0.19601935860219369 * tmp_23;
      real_t tmp_414 = tmp_19 * ( tmp_20 + tmp_413 );
      real_t tmp_415 = 0.19601935860219369 * tmp_29 + 0.19601935860219369 * tmp_30;
      real_t tmp_416 = tmp_19 * ( tmp_27 + tmp_415 );
      real_t tmp_417 = 0.19601935860219369 * tmp_36 + 0.19601935860219369 * tmp_37;
      real_t tmp_418 = tmp_19 * ( tmp_34 + tmp_417 );
      real_t tmp_419 = tmp_26 * tmp_416 + tmp_33 * tmp_418 + tmp_414 * tmp_9 - 1.0 / 4.0;
      real_t tmp_420 = tmp_41 * tmp_414 + tmp_416 * tmp_42 + tmp_418 * tmp_43 - 1.0 / 4.0;
      real_t tmp_421 = tmp_414 * tmp_45 + tmp_416 * tmp_46 + tmp_418 * tmp_47 - 1.0 / 4.0;
      real_t tmp_422 = tmp_1 * tmp_419 + tmp_2 * tmp_420 + tmp_421 * tmp_6;
      real_t tmp_423 = tmp_14 * tmp_419 + tmp_4 * tmp_421 + tmp_420 * tmp_7;
      real_t tmp_424 = tmp_11 * tmp_421 + tmp_13 * tmp_419 + tmp_15 * tmp_420;
      real_t tmp_425 = tmp_69 * ( tmp_413 + tmp_92 );
      real_t tmp_426 = tmp_69 * ( tmp_415 + tmp_94 );
      real_t tmp_427 = tmp_69 * ( tmp_417 + tmp_96 );
      real_t tmp_428 = tmp_425 * tmp_78 + tmp_426 * tmp_75 + tmp_427 * tmp_58 - 1.0 / 4.0;
      real_t tmp_429 = tmp_425 * tmp_79 + tmp_426 * tmp_76 + tmp_427 * tmp_71 - 1.0 / 4.0;
      real_t tmp_430 = tmp_425 * tmp_80 + tmp_426 * tmp_77 + tmp_427 * tmp_73 - 1.0 / 4.0;
      real_t tmp_431 = 0.40446199974765351 * tmp_22 + 0.40446199974765351 * tmp_23;
      real_t tmp_432 = tmp_19 * ( tmp_20 + tmp_431 );
      real_t tmp_433 = 0.40446199974765351 * tmp_29 + 0.40446199974765351 * tmp_30;
      real_t tmp_434 = tmp_19 * ( tmp_27 + tmp_433 );
      real_t tmp_435 = 0.40446199974765351 * tmp_36 + 0.40446199974765351 * tmp_37;
      real_t tmp_436 = tmp_19 * ( tmp_34 + tmp_435 );
      real_t tmp_437 = tmp_26 * tmp_434 + tmp_33 * tmp_436 + tmp_432 * tmp_9 - 1.0 / 4.0;
      real_t tmp_438 = tmp_41 * tmp_432 + tmp_42 * tmp_434 + tmp_43 * tmp_436 - 1.0 / 4.0;
      real_t tmp_439 = tmp_432 * tmp_45 + tmp_434 * tmp_46 + tmp_436 * tmp_47 - 1.0 / 4.0;
      real_t tmp_440 = tmp_1 * tmp_437 + tmp_2 * tmp_438 + tmp_439 * tmp_6;
      real_t tmp_441 = tmp_14 * tmp_437 + tmp_4 * tmp_439 + tmp_438 * tmp_7;
      real_t tmp_442 = tmp_11 * tmp_439 + tmp_13 * tmp_437 + tmp_15 * tmp_438;
      real_t tmp_443 = tmp_69 * ( tmp_431 + tmp_92 );
      real_t tmp_444 = tmp_69 * ( tmp_433 + tmp_94 );
      real_t tmp_445 = tmp_69 * ( tmp_435 + tmp_96 );
      real_t tmp_446 = tmp_443 * tmp_78 + tmp_444 * tmp_75 + tmp_445 * tmp_58 - 1.0 / 4.0;
      real_t tmp_447 = tmp_443 * tmp_79 + tmp_444 * tmp_76 + tmp_445 * tmp_71 - 1.0 / 4.0;
      real_t tmp_448 = tmp_443 * tmp_80 + tmp_444 * tmp_77 + tmp_445 * tmp_73 - 1.0 / 4.0;
      real_t tmp_449 = 0.1711304259088916 * tmp_22 + 0.041227165399737475 * tmp_23;
      real_t tmp_450 = tmp_19 * ( tmp_20 + tmp_449 );
      real_t tmp_451 = 0.1711304259088916 * tmp_29 + 0.041227165399737475 * tmp_30;
      real_t tmp_452 = tmp_19 * ( tmp_27 + tmp_451 );
      real_t tmp_453 = 0.1711304259088916 * tmp_36 + 0.041227165399737475 * tmp_37;
      real_t tmp_454 = tmp_19 * ( tmp_34 + tmp_453 );
      real_t tmp_455 = tmp_26 * tmp_452 + tmp_33 * tmp_454 + tmp_450 * tmp_9 - 1.0 / 4.0;
      real_t tmp_456 = tmp_41 * tmp_450 + tmp_42 * tmp_452 + tmp_43 * tmp_454 - 1.0 / 4.0;
      real_t tmp_457 = tmp_45 * tmp_450 + tmp_452 * tmp_46 + tmp_454 * tmp_47 - 1.0 / 4.0;
      real_t tmp_458 = tmp_1 * tmp_455 + tmp_2 * tmp_456 + tmp_457 * tmp_6;
      real_t tmp_459 = tmp_14 * tmp_455 + tmp_4 * tmp_457 + tmp_456 * tmp_7;
      real_t tmp_460 = tmp_11 * tmp_457 + tmp_13 * tmp_455 + tmp_15 * tmp_456;
      real_t tmp_461 = tmp_69 * ( tmp_449 + tmp_92 );
      real_t tmp_462 = tmp_69 * ( tmp_451 + tmp_94 );
      real_t tmp_463 = tmp_69 * ( tmp_453 + tmp_96 );
      real_t tmp_464 = tmp_461 * tmp_78 + tmp_462 * tmp_75 + tmp_463 * tmp_58 - 1.0 / 4.0;
      real_t tmp_465 = tmp_461 * tmp_79 + tmp_462 * tmp_76 + tmp_463 * tmp_71 - 1.0 / 4.0;
      real_t tmp_466 = tmp_461 * tmp_80 + tmp_462 * tmp_77 + tmp_463 * tmp_73 - 1.0 / 4.0;
      real_t a_0_0   = 0.037198804536718075 * tmp_106 *
                         ( -tmp_105 * ( tmp_116 * ( tmp_122 * tmp_60 + tmp_123 * tmp_61 + tmp_124 * tmp_64 ) +
                                        tmp_117 * ( tmp_122 * tmp_65 + tmp_123 * tmp_51 + tmp_124 * tmp_55 ) +
                                        tmp_118 * ( tmp_122 * tmp_62 + tmp_123 * tmp_56 + tmp_124 * tmp_53 ) ) -
                           tmp_116 * tmp_81 - tmp_117 * tmp_86 - tmp_118 * tmp_91 ) +
                     0.020848748529055869 * tmp_106 *
                         ( -tmp_105 * ( tmp_134 * ( tmp_140 * tmp_60 + tmp_141 * tmp_61 + tmp_142 * tmp_64 ) +
                                        tmp_135 * ( tmp_140 * tmp_65 + tmp_141 * tmp_51 + tmp_142 * tmp_55 ) +
                                        tmp_136 * ( tmp_140 * tmp_62 + tmp_141 * tmp_56 + tmp_142 * tmp_53 ) ) -
                           tmp_134 * tmp_81 - tmp_135 * tmp_86 - tmp_136 * tmp_91 ) +
                     0.019202922745021479 * tmp_106 *
                         ( -tmp_105 * ( tmp_152 * ( tmp_158 * tmp_60 + tmp_159 * tmp_61 + tmp_160 * tmp_64 ) +
                                        tmp_153 * ( tmp_158 * tmp_65 + tmp_159 * tmp_51 + tmp_160 * tmp_55 ) +
                                        tmp_154 * ( tmp_158 * tmp_62 + tmp_159 * tmp_56 + tmp_160 * tmp_53 ) ) -
                           tmp_152 * tmp_81 - tmp_153 * tmp_86 - tmp_154 * tmp_91 ) +
                     0.020848748529055869 * tmp_106 *
                         ( -tmp_105 * ( tmp_170 * ( tmp_176 * tmp_60 + tmp_177 * tmp_61 + tmp_178 * tmp_64 ) +
                                        tmp_171 * ( tmp_176 * tmp_65 + tmp_177 * tmp_51 + tmp_178 * tmp_55 ) +
                                        tmp_172 * ( tmp_176 * tmp_62 + tmp_177 * tmp_56 + tmp_178 * tmp_53 ) ) -
                           tmp_170 * tmp_81 - tmp_171 * tmp_86 - tmp_172 * tmp_91 ) +
                     0.019202922745021479 * tmp_106 *
                         ( -tmp_105 * ( tmp_188 * ( tmp_194 * tmp_60 + tmp_195 * tmp_61 + tmp_196 * tmp_64 ) +
                                        tmp_189 * ( tmp_194 * tmp_65 + tmp_195 * tmp_51 + tmp_196 * tmp_55 ) +
                                        tmp_190 * ( tmp_194 * tmp_62 + tmp_195 * tmp_56 + tmp_196 * tmp_53 ) ) -
                           tmp_188 * tmp_81 - tmp_189 * tmp_86 - tmp_190 * tmp_91 ) +
                     0.020848748529055869 * tmp_106 *
                         ( -tmp_105 * ( tmp_206 * ( tmp_212 * tmp_60 + tmp_213 * tmp_61 + tmp_214 * tmp_64 ) +
                                        tmp_207 * ( tmp_212 * tmp_65 + tmp_213 * tmp_51 + tmp_214 * tmp_55 ) +
                                        tmp_208 * ( tmp_212 * tmp_62 + tmp_213 * tmp_56 + tmp_214 * tmp_53 ) ) -
                           tmp_206 * tmp_81 - tmp_207 * tmp_86 - tmp_208 * tmp_91 ) +
                     0.019202922745021479 * tmp_106 *
                         ( -tmp_105 * ( tmp_224 * ( tmp_230 * tmp_60 + tmp_231 * tmp_61 + tmp_232 * tmp_64 ) +
                                        tmp_225 * ( tmp_230 * tmp_65 + tmp_231 * tmp_51 + tmp_232 * tmp_55 ) +
                                        tmp_226 * ( tmp_230 * tmp_62 + tmp_231 * tmp_56 + tmp_232 * tmp_53 ) ) -
                           tmp_224 * tmp_81 - tmp_225 * tmp_86 - tmp_226 * tmp_91 ) +
                     0.020848748529055869 * tmp_106 *
                         ( -tmp_105 * ( tmp_242 * ( tmp_248 * tmp_60 + tmp_249 * tmp_61 + tmp_250 * tmp_64 ) +
                                        tmp_243 * ( tmp_248 * tmp_65 + tmp_249 * tmp_51 + tmp_250 * tmp_55 ) +
                                        tmp_244 * ( tmp_248 * tmp_62 + tmp_249 * tmp_56 + tmp_250 * tmp_53 ) ) -
                           tmp_242 * tmp_81 - tmp_243 * tmp_86 - tmp_244 * tmp_91 ) +
                     0.019202922745021479 * tmp_106 *
                         ( -tmp_105 * ( tmp_260 * ( tmp_266 * tmp_60 + tmp_267 * tmp_61 + tmp_268 * tmp_64 ) +
                                        tmp_261 * ( tmp_266 * tmp_65 + tmp_267 * tmp_51 + tmp_268 * tmp_55 ) +
                                        tmp_262 * ( tmp_266 * tmp_62 + tmp_267 * tmp_56 + tmp_268 * tmp_53 ) ) -
                           tmp_260 * tmp_81 - tmp_261 * tmp_86 - tmp_262 * tmp_91 ) +
                     0.020848748529055869 * tmp_106 *
                         ( -tmp_105 * ( tmp_278 * ( tmp_284 * tmp_60 + tmp_285 * tmp_61 + tmp_286 * tmp_64 ) +
                                        tmp_279 * ( tmp_284 * tmp_65 + tmp_285 * tmp_51 + tmp_286 * tmp_55 ) +
                                        tmp_280 * ( tmp_284 * tmp_62 + tmp_285 * tmp_56 + tmp_286 * tmp_53 ) ) -
                           tmp_278 * tmp_81 - tmp_279 * tmp_86 - tmp_280 * tmp_91 ) +
                     0.019202922745021479 * tmp_106 *
                         ( -tmp_105 * ( tmp_296 * ( tmp_302 * tmp_60 + tmp_303 * tmp_61 + tmp_304 * tmp_64 ) +
                                        tmp_297 * ( tmp_302 * tmp_65 + tmp_303 * tmp_51 + tmp_304 * tmp_55 ) +
                                        tmp_298 * ( tmp_302 * tmp_62 + tmp_303 * tmp_56 + tmp_304 * tmp_53 ) ) -
                           tmp_296 * tmp_81 - tmp_297 * tmp_86 - tmp_298 * tmp_91 ) +
                     0.042507265838595799 * tmp_106 *
                         ( -tmp_105 * ( tmp_314 * ( tmp_320 * tmp_60 + tmp_321 * tmp_61 + tmp_322 * tmp_64 ) +
                                        tmp_315 * ( tmp_320 * tmp_65 + tmp_321 * tmp_51 + tmp_322 * tmp_55 ) +
                                        tmp_316 * ( tmp_320 * tmp_62 + tmp_321 * tmp_56 + tmp_322 * tmp_53 ) ) -
                           tmp_314 * tmp_81 - tmp_315 * tmp_86 - tmp_316 * tmp_91 ) +
                     0.020848748529055869 * tmp_106 *
                         ( -tmp_105 * ( tmp_332 * ( tmp_338 * tmp_60 + tmp_339 * tmp_61 + tmp_340 * tmp_64 ) +
                                        tmp_333 * ( tmp_338 * tmp_65 + tmp_339 * tmp_51 + tmp_340 * tmp_55 ) +
                                        tmp_334 * ( tmp_338 * tmp_62 + tmp_339 * tmp_56 + tmp_340 * tmp_53 ) ) -
                           tmp_332 * tmp_81 - tmp_333 * tmp_86 - tmp_334 * tmp_91 ) +
                     0.0068572537431980923 * tmp_106 *
                         ( -tmp_105 * ( tmp_350 * ( tmp_356 * tmp_60 + tmp_357 * tmp_61 + tmp_358 * tmp_64 ) +
                                        tmp_351 * ( tmp_356 * tmp_65 + tmp_357 * tmp_51 + tmp_358 * tmp_55 ) +
                                        tmp_352 * ( tmp_356 * tmp_62 + tmp_357 * tmp_56 + tmp_358 * tmp_53 ) ) -
                           tmp_350 * tmp_81 - tmp_351 * tmp_86 - tmp_352 * tmp_91 ) +
                     0.037198804536718075 * tmp_106 *
                         ( -tmp_105 * ( tmp_368 * ( tmp_374 * tmp_60 + tmp_375 * tmp_61 + tmp_376 * tmp_64 ) +
                                        tmp_369 * ( tmp_374 * tmp_65 + tmp_375 * tmp_51 + tmp_376 * tmp_55 ) +
                                        tmp_370 * ( tmp_374 * tmp_62 + tmp_375 * tmp_56 + tmp_376 * tmp_53 ) ) -
                           tmp_368 * tmp_81 - tmp_369 * tmp_86 - tmp_370 * tmp_91 ) +
                     0.042507265838595799 * tmp_106 *
                         ( -tmp_105 * ( tmp_386 * ( tmp_392 * tmp_60 + tmp_393 * tmp_61 + tmp_394 * tmp_64 ) +
                                        tmp_387 * ( tmp_392 * tmp_65 + tmp_393 * tmp_51 + tmp_394 * tmp_55 ) +
                                        tmp_388 * ( tmp_392 * tmp_62 + tmp_393 * tmp_56 + tmp_394 * tmp_53 ) ) -
                           tmp_386 * tmp_81 - tmp_387 * tmp_86 - tmp_388 * tmp_91 ) +
                     0.0068572537431980923 * tmp_106 *
                         ( -tmp_105 * ( tmp_404 * ( tmp_410 * tmp_60 + tmp_411 * tmp_61 + tmp_412 * tmp_64 ) +
                                        tmp_405 * ( tmp_410 * tmp_65 + tmp_411 * tmp_51 + tmp_412 * tmp_55 ) +
                                        tmp_406 * ( tmp_410 * tmp_62 + tmp_411 * tmp_56 + tmp_412 * tmp_53 ) ) -
                           tmp_404 * tmp_81 - tmp_405 * tmp_86 - tmp_406 * tmp_91 ) +
                     0.037198804536718075 * tmp_106 *
                         ( -tmp_105 * ( tmp_422 * ( tmp_428 * tmp_60 + tmp_429 * tmp_61 + tmp_430 * tmp_64 ) +
                                        tmp_423 * ( tmp_428 * tmp_65 + tmp_429 * tmp_51 + tmp_430 * tmp_55 ) +
                                        tmp_424 * ( tmp_428 * tmp_62 + tmp_429 * tmp_56 + tmp_430 * tmp_53 ) ) -
                           tmp_422 * tmp_81 - tmp_423 * tmp_86 - tmp_424 * tmp_91 ) +
                     0.042507265838595799 * tmp_106 *
                         ( -tmp_105 * ( tmp_440 * ( tmp_446 * tmp_60 + tmp_447 * tmp_61 + tmp_448 * tmp_64 ) +
                                        tmp_441 * ( tmp_446 * tmp_65 + tmp_447 * tmp_51 + tmp_448 * tmp_55 ) +
                                        tmp_442 * ( tmp_446 * tmp_62 + tmp_447 * tmp_56 + tmp_448 * tmp_53 ) ) -
                           tmp_440 * tmp_81 - tmp_441 * tmp_86 - tmp_442 * tmp_91 ) +
                     0.019202922745021479 * tmp_106 *
                         ( -tmp_105 * ( tmp_458 * ( tmp_464 * tmp_60 + tmp_465 * tmp_61 + tmp_466 * tmp_64 ) +
                                        tmp_459 * ( tmp_464 * tmp_65 + tmp_465 * tmp_51 + tmp_466 * tmp_55 ) +
                                        tmp_460 * ( tmp_464 * tmp_62 + tmp_465 * tmp_56 + tmp_466 * tmp_53 ) ) -
                           tmp_458 * tmp_81 - tmp_459 * tmp_86 - tmp_460 * tmp_91 ) +
                     0.0068572537431980923 * tmp_106 *
                         ( -tmp_105 * ( tmp_49 * ( tmp_100 * tmp_64 + tmp_60 * tmp_98 + tmp_61 * tmp_99 ) +
                                        tmp_82 * ( tmp_100 * tmp_55 + tmp_51 * tmp_99 + tmp_65 * tmp_98 ) +
                                        tmp_87 * ( tmp_100 * tmp_53 + tmp_56 * tmp_99 + tmp_62 * tmp_98 ) ) -
                           tmp_49 * tmp_81 - tmp_82 * tmp_86 - tmp_87 * tmp_91 );
      elMat( 0, 0 ) = a_0_0;
   }

   void integrateFacetDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                           const std::vector< Point3D >& coordsFacet,
                                           const Point3D&,
                                           const Point3D&     outwardNormal,
                                           const DGBasisInfo& trialBasis,
                                           const DGBasisInfo& testBasis,
                                           int                trialDegree,
                                           int                testDegree,
                                           MatrixXr&          elMat ) const override
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = p_affine_2_0 + tmp_0;
      real_t tmp_3  = -p_affine_0_1;
      real_t tmp_4  = p_affine_3_1 + tmp_3;
      real_t tmp_5  = tmp_2 * tmp_4;
      real_t tmp_6  = p_affine_3_0 + tmp_0;
      real_t tmp_7  = p_affine_2_1 + tmp_3;
      real_t tmp_8  = tmp_6 * tmp_7;
      real_t tmp_9  = tmp_5 - tmp_8;
      real_t tmp_10 = -p_affine_0_2;
      real_t tmp_11 = p_affine_3_2 + tmp_10;
      real_t tmp_12 = tmp_11 * tmp_7;
      real_t tmp_13 = p_affine_1_2 + tmp_10;
      real_t tmp_14 = p_affine_1_1 + tmp_3;
      real_t tmp_15 = p_affine_2_2 + tmp_10;
      real_t tmp_16 = tmp_15 * tmp_6;
      real_t tmp_17 = tmp_15 * tmp_4;
      real_t tmp_18 = tmp_11 * tmp_2;
      real_t tmp_19 =
          1.0 / ( tmp_1 * tmp_12 - tmp_1 * tmp_17 + tmp_13 * tmp_5 - tmp_13 * tmp_8 + tmp_14 * tmp_16 - tmp_14 * tmp_18 );
      real_t tmp_20 = -p_affine_8_2;
      real_t tmp_21 = p_affine_9_2 + tmp_20;
      real_t tmp_22 = p_affine_10_2 + tmp_20;
      real_t tmp_23 = p_affine_8_2 + tmp_10;
      real_t tmp_24 = tmp_19 * ( 0.031405749086161582 * tmp_21 + 0.93718850182767688 * tmp_22 + tmp_23 );
      real_t tmp_25 = tmp_16 - tmp_18;
      real_t tmp_26 = -p_affine_8_1;
      real_t tmp_27 = p_affine_9_1 + tmp_26;
      real_t tmp_28 = p_affine_10_1 + tmp_26;
      real_t tmp_29 = p_affine_8_1 + tmp_3;
      real_t tmp_30 = tmp_19 * ( 0.031405749086161582 * tmp_27 + 0.93718850182767688 * tmp_28 + tmp_29 );
      real_t tmp_31 = tmp_12 - tmp_17;
      real_t tmp_32 = -p_affine_8_0;
      real_t tmp_33 = p_affine_9_0 + tmp_32;
      real_t tmp_34 = p_affine_10_0 + tmp_32;
      real_t tmp_35 = p_affine_8_0 + tmp_0;
      real_t tmp_36 = tmp_19 * ( 0.031405749086161582 * tmp_33 + 0.93718850182767688 * tmp_34 + tmp_35 );
      real_t tmp_37 = tmp_24 * tmp_9 + tmp_25 * tmp_30 + tmp_31 * tmp_36 - 1.0 / 4.0;
      real_t tmp_38 = -tmp_1 * tmp_4 + tmp_14 * tmp_6;
      real_t tmp_39 = tmp_1 * tmp_11 - tmp_13 * tmp_6;
      real_t tmp_40 = -tmp_11 * tmp_14 + tmp_13 * tmp_4;
      real_t tmp_41 = tmp_24 * tmp_38 + tmp_30 * tmp_39 + tmp_36 * tmp_40 - 1.0 / 4.0;
      real_t tmp_42 = tmp_1 * tmp_7 - tmp_14 * tmp_2;
      real_t tmp_43 = -tmp_1 * tmp_15 + tmp_13 * tmp_2;
      real_t tmp_44 = -tmp_13 * tmp_7 + tmp_14 * tmp_15;
      real_t tmp_45 = tmp_24 * tmp_42 + tmp_30 * tmp_43 + tmp_36 * tmp_44 - 1.0 / 4.0;
      real_t tmp_46 = p_affine_8_1 - p_affine_9_1;
      real_t tmp_47 = p_affine_8_0 - p_affine_9_0;
      real_t tmp_48 = p_affine_8_2 - p_affine_9_2;
      real_t tmp_49 =
          1.0 * std::pow( ( std::abs( tmp_22 * tmp_46 - tmp_28 * tmp_48 ) * std::abs( tmp_22 * tmp_46 - tmp_28 * tmp_48 ) ) +
                              ( std::abs( tmp_22 * tmp_47 - tmp_34 * tmp_48 ) * std::abs( tmp_22 * tmp_47 - tmp_34 * tmp_48 ) ) +
                              ( std::abs( tmp_28 * tmp_47 - tmp_34 * tmp_46 ) * std::abs( tmp_28 * tmp_47 - tmp_34 * tmp_46 ) ),
                          0.25 );
      real_t tmp_50  = tmp_19 * ( 0.19601935860219369 * tmp_21 + 0.60796128279561268 * tmp_22 + tmp_23 );
      real_t tmp_51  = tmp_19 * ( 0.19601935860219369 * tmp_27 + 0.60796128279561268 * tmp_28 + tmp_29 );
      real_t tmp_52  = tmp_19 * ( 0.19601935860219369 * tmp_33 + 0.60796128279561268 * tmp_34 + tmp_35 );
      real_t tmp_53  = tmp_25 * tmp_51 + tmp_31 * tmp_52 + tmp_50 * tmp_9 - 1.0 / 4.0;
      real_t tmp_54  = tmp_38 * tmp_50 + tmp_39 * tmp_51 + tmp_40 * tmp_52 - 1.0 / 4.0;
      real_t tmp_55  = tmp_42 * tmp_50 + tmp_43 * tmp_51 + tmp_44 * tmp_52 - 1.0 / 4.0;
      real_t tmp_56  = tmp_19 * ( 0.37605877282253791 * tmp_21 + 0.039308471900058539 * tmp_22 + tmp_23 );
      real_t tmp_57  = tmp_19 * ( 0.37605877282253791 * tmp_27 + 0.039308471900058539 * tmp_28 + tmp_29 );
      real_t tmp_58  = tmp_19 * ( 0.37605877282253791 * tmp_33 + 0.039308471900058539 * tmp_34 + tmp_35 );
      real_t tmp_59  = tmp_25 * tmp_57 + tmp_31 * tmp_58 + tmp_56 * tmp_9 - 1.0 / 4.0;
      real_t tmp_60  = tmp_38 * tmp_56 + tmp_39 * tmp_57 + tmp_40 * tmp_58 - 1.0 / 4.0;
      real_t tmp_61  = tmp_42 * tmp_56 + tmp_43 * tmp_57 + tmp_44 * tmp_58 - 1.0 / 4.0;
      real_t tmp_62  = tmp_19 * ( 0.78764240869137092 * tmp_21 + 0.1711304259088916 * tmp_22 + tmp_23 );
      real_t tmp_63  = tmp_19 * ( 0.78764240869137092 * tmp_27 + 0.1711304259088916 * tmp_28 + tmp_29 );
      real_t tmp_64  = tmp_19 * ( 0.78764240869137092 * tmp_33 + 0.1711304259088916 * tmp_34 + tmp_35 );
      real_t tmp_65  = tmp_25 * tmp_63 + tmp_31 * tmp_64 + tmp_62 * tmp_9 - 1.0 / 4.0;
      real_t tmp_66  = tmp_38 * tmp_62 + tmp_39 * tmp_63 + tmp_40 * tmp_64 - 1.0 / 4.0;
      real_t tmp_67  = tmp_42 * tmp_62 + tmp_43 * tmp_63 + tmp_44 * tmp_64 - 1.0 / 4.0;
      real_t tmp_68  = tmp_19 * ( 0.58463275527740355 * tmp_21 + 0.37605877282253791 * tmp_22 + tmp_23 );
      real_t tmp_69  = tmp_19 * ( 0.58463275527740355 * tmp_27 + 0.37605877282253791 * tmp_28 + tmp_29 );
      real_t tmp_70  = tmp_19 * ( 0.58463275527740355 * tmp_33 + 0.37605877282253791 * tmp_34 + tmp_35 );
      real_t tmp_71  = tmp_25 * tmp_69 + tmp_31 * tmp_70 + tmp_68 * tmp_9 - 1.0 / 4.0;
      real_t tmp_72  = tmp_38 * tmp_68 + tmp_39 * tmp_69 + tmp_40 * tmp_70 - 1.0 / 4.0;
      real_t tmp_73  = tmp_42 * tmp_68 + tmp_43 * tmp_69 + tmp_44 * tmp_70 - 1.0 / 4.0;
      real_t tmp_74  = tmp_19 * ( 0.041227165399737475 * tmp_21 + 0.78764240869137092 * tmp_22 + tmp_23 );
      real_t tmp_75  = tmp_19 * ( 0.041227165399737475 * tmp_27 + 0.78764240869137092 * tmp_28 + tmp_29 );
      real_t tmp_76  = tmp_19 * ( 0.041227165399737475 * tmp_33 + 0.78764240869137092 * tmp_34 + tmp_35 );
      real_t tmp_77  = tmp_25 * tmp_75 + tmp_31 * tmp_76 + tmp_74 * tmp_9 - 1.0 / 4.0;
      real_t tmp_78  = tmp_38 * tmp_74 + tmp_39 * tmp_75 + tmp_40 * tmp_76 - 1.0 / 4.0;
      real_t tmp_79  = tmp_42 * tmp_74 + tmp_43 * tmp_75 + tmp_44 * tmp_76 - 1.0 / 4.0;
      real_t tmp_80  = tmp_19 * ( 0.039308471900058539 * tmp_21 + 0.58463275527740355 * tmp_22 + tmp_23 );
      real_t tmp_81  = tmp_19 * ( 0.039308471900058539 * tmp_27 + 0.58463275527740355 * tmp_28 + tmp_29 );
      real_t tmp_82  = tmp_19 * ( 0.039308471900058539 * tmp_33 + 0.58463275527740355 * tmp_34 + tmp_35 );
      real_t tmp_83  = tmp_25 * tmp_81 + tmp_31 * tmp_82 + tmp_80 * tmp_9 - 1.0 / 4.0;
      real_t tmp_84  = tmp_38 * tmp_80 + tmp_39 * tmp_81 + tmp_40 * tmp_82 - 1.0 / 4.0;
      real_t tmp_85  = tmp_42 * tmp_80 + tmp_43 * tmp_81 + tmp_44 * tmp_82 - 1.0 / 4.0;
      real_t tmp_86  = tmp_19 * ( 0.78764240869137092 * tmp_21 + 0.041227165399737475 * tmp_22 + tmp_23 );
      real_t tmp_87  = tmp_19 * ( 0.78764240869137092 * tmp_27 + 0.041227165399737475 * tmp_28 + tmp_29 );
      real_t tmp_88  = tmp_19 * ( 0.78764240869137092 * tmp_33 + 0.041227165399737475 * tmp_34 + tmp_35 );
      real_t tmp_89  = tmp_25 * tmp_87 + tmp_31 * tmp_88 + tmp_86 * tmp_9 - 1.0 / 4.0;
      real_t tmp_90  = tmp_38 * tmp_86 + tmp_39 * tmp_87 + tmp_40 * tmp_88 - 1.0 / 4.0;
      real_t tmp_91  = tmp_42 * tmp_86 + tmp_43 * tmp_87 + tmp_44 * tmp_88 - 1.0 / 4.0;
      real_t tmp_92  = tmp_19 * ( 0.58463275527740355 * tmp_21 + 0.039308471900058539 * tmp_22 + tmp_23 );
      real_t tmp_93  = tmp_19 * ( 0.58463275527740355 * tmp_27 + 0.039308471900058539 * tmp_28 + tmp_29 );
      real_t tmp_94  = tmp_19 * ( 0.58463275527740355 * tmp_33 + 0.039308471900058539 * tmp_34 + tmp_35 );
      real_t tmp_95  = tmp_25 * tmp_93 + tmp_31 * tmp_94 + tmp_9 * tmp_92 - 1.0 / 4.0;
      real_t tmp_96  = tmp_38 * tmp_92 + tmp_39 * tmp_93 + tmp_40 * tmp_94 - 1.0 / 4.0;
      real_t tmp_97  = tmp_42 * tmp_92 + tmp_43 * tmp_93 + tmp_44 * tmp_94 - 1.0 / 4.0;
      real_t tmp_98  = tmp_19 * ( 0.1711304259088916 * tmp_21 + 0.78764240869137092 * tmp_22 + tmp_23 );
      real_t tmp_99  = tmp_19 * ( 0.1711304259088916 * tmp_27 + 0.78764240869137092 * tmp_28 + tmp_29 );
      real_t tmp_100 = tmp_19 * ( 0.1711304259088916 * tmp_33 + 0.78764240869137092 * tmp_34 + tmp_35 );
      real_t tmp_101 = tmp_100 * tmp_31 + tmp_25 * tmp_99 + tmp_9 * tmp_98 - 1.0 / 4.0;
      real_t tmp_102 = tmp_100 * tmp_40 + tmp_38 * tmp_98 + tmp_39 * tmp_99 - 1.0 / 4.0;
      real_t tmp_103 = tmp_100 * tmp_44 + tmp_42 * tmp_98 + tmp_43 * tmp_99 - 1.0 / 4.0;
      real_t tmp_104 = tmp_19 * ( 0.37605877282253791 * tmp_21 + 0.58463275527740355 * tmp_22 + tmp_23 );
      real_t tmp_105 = tmp_19 * ( 0.37605877282253791 * tmp_27 + 0.58463275527740355 * tmp_28 + tmp_29 );
      real_t tmp_106 = tmp_19 * ( 0.37605877282253791 * tmp_33 + 0.58463275527740355 * tmp_34 + tmp_35 );
      real_t tmp_107 = tmp_104 * tmp_9 + tmp_105 * tmp_25 + tmp_106 * tmp_31 - 1.0 / 4.0;
      real_t tmp_108 = tmp_104 * tmp_38 + tmp_105 * tmp_39 + tmp_106 * tmp_40 - 1.0 / 4.0;
      real_t tmp_109 = tmp_104 * tmp_42 + tmp_105 * tmp_43 + tmp_106 * tmp_44 - 1.0 / 4.0;
      real_t tmp_110 = tmp_19 * ( 0.041227165399737475 * tmp_21 + 0.1711304259088916 * tmp_22 + tmp_23 );
      real_t tmp_111 = tmp_19 * ( 0.041227165399737475 * tmp_27 + 0.1711304259088916 * tmp_28 + tmp_29 );
      real_t tmp_112 = tmp_19 * ( 0.041227165399737475 * tmp_33 + 0.1711304259088916 * tmp_34 + tmp_35 );
      real_t tmp_113 = tmp_110 * tmp_9 + tmp_111 * tmp_25 + tmp_112 * tmp_31 - 1.0 / 4.0;
      real_t tmp_114 = tmp_110 * tmp_38 + tmp_111 * tmp_39 + tmp_112 * tmp_40 - 1.0 / 4.0;
      real_t tmp_115 = tmp_110 * tmp_42 + tmp_111 * tmp_43 + tmp_112 * tmp_44 - 1.0 / 4.0;
      real_t tmp_116 = tmp_19 * ( 0.40446199974765351 * tmp_21 + 0.19107600050469298 * tmp_22 + tmp_23 );
      real_t tmp_117 = tmp_19 * ( 0.40446199974765351 * tmp_27 + 0.19107600050469298 * tmp_28 + tmp_29 );
      real_t tmp_118 = tmp_19 * ( 0.40446199974765351 * tmp_33 + 0.19107600050469298 * tmp_34 + tmp_35 );
      real_t tmp_119 = tmp_116 * tmp_9 + tmp_117 * tmp_25 + tmp_118 * tmp_31 - 1.0 / 4.0;
      real_t tmp_120 = tmp_116 * tmp_38 + tmp_117 * tmp_39 + tmp_118 * tmp_40 - 1.0 / 4.0;
      real_t tmp_121 = tmp_116 * tmp_42 + tmp_117 * tmp_43 + tmp_118 * tmp_44 - 1.0 / 4.0;
      real_t tmp_122 = tmp_19 * ( 0.039308471900058539 * tmp_21 + 0.37605877282253791 * tmp_22 + tmp_23 );
      real_t tmp_123 = tmp_19 * ( 0.039308471900058539 * tmp_27 + 0.37605877282253791 * tmp_28 + tmp_29 );
      real_t tmp_124 = tmp_19 * ( 0.039308471900058539 * tmp_33 + 0.37605877282253791 * tmp_34 + tmp_35 );
      real_t tmp_125 = tmp_122 * tmp_9 + tmp_123 * tmp_25 + tmp_124 * tmp_31 - 1.0 / 4.0;
      real_t tmp_126 = tmp_122 * tmp_38 + tmp_123 * tmp_39 + tmp_124 * tmp_40 - 1.0 / 4.0;
      real_t tmp_127 = tmp_122 * tmp_42 + tmp_123 * tmp_43 + tmp_124 * tmp_44 - 1.0 / 4.0;
      real_t tmp_128 = tmp_19 * ( 0.93718850182767688 * tmp_21 + 0.031405749086161582 * tmp_22 + tmp_23 );
      real_t tmp_129 = tmp_19 * ( 0.93718850182767688 * tmp_27 + 0.031405749086161582 * tmp_28 + tmp_29 );
      real_t tmp_130 = tmp_19 * ( 0.93718850182767688 * tmp_33 + 0.031405749086161582 * tmp_34 + tmp_35 );
      real_t tmp_131 = tmp_128 * tmp_9 + tmp_129 * tmp_25 + tmp_130 * tmp_31 - 1.0 / 4.0;
      real_t tmp_132 = tmp_128 * tmp_38 + tmp_129 * tmp_39 + tmp_130 * tmp_40 - 1.0 / 4.0;
      real_t tmp_133 = tmp_128 * tmp_42 + tmp_129 * tmp_43 + tmp_130 * tmp_44 - 1.0 / 4.0;
      real_t tmp_134 = tmp_19 * ( 0.60796128279561268 * tmp_21 + 0.19601935860219369 * tmp_22 + tmp_23 );
      real_t tmp_135 = tmp_19 * ( 0.60796128279561268 * tmp_27 + 0.19601935860219369 * tmp_28 + tmp_29 );
      real_t tmp_136 = tmp_19 * ( 0.60796128279561268 * tmp_33 + 0.19601935860219369 * tmp_34 + tmp_35 );
      real_t tmp_137 = tmp_134 * tmp_9 + tmp_135 * tmp_25 + tmp_136 * tmp_31 - 1.0 / 4.0;
      real_t tmp_138 = tmp_134 * tmp_38 + tmp_135 * tmp_39 + tmp_136 * tmp_40 - 1.0 / 4.0;
      real_t tmp_139 = tmp_134 * tmp_42 + tmp_135 * tmp_43 + tmp_136 * tmp_44 - 1.0 / 4.0;
      real_t tmp_140 = tmp_19 * ( 0.19107600050469298 * tmp_21 + 0.40446199974765351 * tmp_22 + tmp_23 );
      real_t tmp_141 = tmp_19 * ( 0.19107600050469298 * tmp_27 + 0.40446199974765351 * tmp_28 + tmp_29 );
      real_t tmp_142 = tmp_19 * ( 0.19107600050469298 * tmp_33 + 0.40446199974765351 * tmp_34 + tmp_35 );
      real_t tmp_143 = tmp_140 * tmp_9 + tmp_141 * tmp_25 + tmp_142 * tmp_31 - 1.0 / 4.0;
      real_t tmp_144 = tmp_140 * tmp_38 + tmp_141 * tmp_39 + tmp_142 * tmp_40 - 1.0 / 4.0;
      real_t tmp_145 = tmp_140 * tmp_42 + tmp_141 * tmp_43 + tmp_142 * tmp_44 - 1.0 / 4.0;
      real_t tmp_146 = tmp_19 * ( 0.031405749086161582 * tmp_21 + 0.031405749086161582 * tmp_22 + tmp_23 );
      real_t tmp_147 = tmp_19 * ( 0.031405749086161582 * tmp_27 + 0.031405749086161582 * tmp_28 + tmp_29 );
      real_t tmp_148 = tmp_19 * ( 0.031405749086161582 * tmp_33 + 0.031405749086161582 * tmp_34 + tmp_35 );
      real_t tmp_149 = tmp_146 * tmp_9 + tmp_147 * tmp_25 + tmp_148 * tmp_31 - 1.0 / 4.0;
      real_t tmp_150 = tmp_146 * tmp_38 + tmp_147 * tmp_39 + tmp_148 * tmp_40 - 1.0 / 4.0;
      real_t tmp_151 = tmp_146 * tmp_42 + tmp_147 * tmp_43 + tmp_148 * tmp_44 - 1.0 / 4.0;
      real_t tmp_152 = tmp_19 * ( 0.19601935860219369 * tmp_21 + 0.19601935860219369 * tmp_22 + tmp_23 );
      real_t tmp_153 = tmp_19 * ( 0.19601935860219369 * tmp_27 + 0.19601935860219369 * tmp_28 + tmp_29 );
      real_t tmp_154 = tmp_19 * ( 0.19601935860219369 * tmp_33 + 0.19601935860219369 * tmp_34 + tmp_35 );
      real_t tmp_155 = tmp_152 * tmp_9 + tmp_153 * tmp_25 + tmp_154 * tmp_31 - 1.0 / 4.0;
      real_t tmp_156 = tmp_152 * tmp_38 + tmp_153 * tmp_39 + tmp_154 * tmp_40 - 1.0 / 4.0;
      real_t tmp_157 = tmp_152 * tmp_42 + tmp_153 * tmp_43 + tmp_154 * tmp_44 - 1.0 / 4.0;
      real_t tmp_158 = tmp_19 * ( 0.40446199974765351 * tmp_21 + 0.40446199974765351 * tmp_22 + tmp_23 );
      real_t tmp_159 = tmp_19 * ( 0.40446199974765351 * tmp_27 + 0.40446199974765351 * tmp_28 + tmp_29 );
      real_t tmp_160 = tmp_19 * ( 0.40446199974765351 * tmp_33 + 0.40446199974765351 * tmp_34 + tmp_35 );
      real_t tmp_161 = tmp_158 * tmp_9 + tmp_159 * tmp_25 + tmp_160 * tmp_31 - 1.0 / 4.0;
      real_t tmp_162 = tmp_158 * tmp_38 + tmp_159 * tmp_39 + tmp_160 * tmp_40 - 1.0 / 4.0;
      real_t tmp_163 = tmp_158 * tmp_42 + tmp_159 * tmp_43 + tmp_160 * tmp_44 - 1.0 / 4.0;
      real_t tmp_164 = tmp_19 * ( 0.1711304259088916 * tmp_21 + 0.041227165399737475 * tmp_22 + tmp_23 );
      real_t tmp_165 = tmp_19 * ( 0.1711304259088916 * tmp_27 + 0.041227165399737475 * tmp_28 + tmp_29 );
      real_t tmp_166 = tmp_19 * ( 0.1711304259088916 * tmp_33 + 0.041227165399737475 * tmp_34 + tmp_35 );
      real_t tmp_167 = tmp_164 * tmp_9 + tmp_165 * tmp_25 + tmp_166 * tmp_31 - 1.0 / 4.0;
      real_t tmp_168 = tmp_164 * tmp_38 + tmp_165 * tmp_39 + tmp_166 * tmp_40 - 1.0 / 4.0;
      real_t tmp_169 = tmp_164 * tmp_42 + tmp_165 * tmp_43 + tmp_166 * tmp_44 - 1.0 / 4.0;
      real_t a_0_0 =
          0.019202922745021479 * tmp_49 *
              ( ( ( tmp_1 * tmp_101 + tmp_102 * tmp_2 + tmp_103 * tmp_6 ) *
                  ( tmp_1 * tmp_101 + tmp_102 * tmp_2 + tmp_103 * tmp_6 ) ) +
                ( ( tmp_101 * tmp_13 + tmp_102 * tmp_15 + tmp_103 * tmp_11 ) *
                  ( tmp_101 * tmp_13 + tmp_102 * tmp_15 + tmp_103 * tmp_11 ) ) +
                ( ( tmp_101 * tmp_14 + tmp_102 * tmp_7 + tmp_103 * tmp_4 ) *
                  ( tmp_101 * tmp_14 + tmp_102 * tmp_7 + tmp_103 * tmp_4 ) ) ) +
          0.020848748529055869 * tmp_49 *
              ( ( ( tmp_1 * tmp_107 + tmp_108 * tmp_2 + tmp_109 * tmp_6 ) *
                  ( tmp_1 * tmp_107 + tmp_108 * tmp_2 + tmp_109 * tmp_6 ) ) +
                ( ( tmp_107 * tmp_13 + tmp_108 * tmp_15 + tmp_109 * tmp_11 ) *
                  ( tmp_107 * tmp_13 + tmp_108 * tmp_15 + tmp_109 * tmp_11 ) ) +
                ( ( tmp_107 * tmp_14 + tmp_108 * tmp_7 + tmp_109 * tmp_4 ) *
                  ( tmp_107 * tmp_14 + tmp_108 * tmp_7 + tmp_109 * tmp_4 ) ) ) +
          0.019202922745021479 * tmp_49 *
              ( ( ( tmp_1 * tmp_113 + tmp_114 * tmp_2 + tmp_115 * tmp_6 ) *
                  ( tmp_1 * tmp_113 + tmp_114 * tmp_2 + tmp_115 * tmp_6 ) ) +
                ( ( tmp_11 * tmp_115 + tmp_113 * tmp_13 + tmp_114 * tmp_15 ) *
                  ( tmp_11 * tmp_115 + tmp_113 * tmp_13 + tmp_114 * tmp_15 ) ) +
                ( ( tmp_113 * tmp_14 + tmp_114 * tmp_7 + tmp_115 * tmp_4 ) *
                  ( tmp_113 * tmp_14 + tmp_114 * tmp_7 + tmp_115 * tmp_4 ) ) ) +
          0.042507265838595799 * tmp_49 *
              ( ( ( tmp_1 * tmp_119 + tmp_120 * tmp_2 + tmp_121 * tmp_6 ) *
                  ( tmp_1 * tmp_119 + tmp_120 * tmp_2 + tmp_121 * tmp_6 ) ) +
                ( ( tmp_11 * tmp_121 + tmp_119 * tmp_13 + tmp_120 * tmp_15 ) *
                  ( tmp_11 * tmp_121 + tmp_119 * tmp_13 + tmp_120 * tmp_15 ) ) +
                ( ( tmp_119 * tmp_14 + tmp_120 * tmp_7 + tmp_121 * tmp_4 ) *
                  ( tmp_119 * tmp_14 + tmp_120 * tmp_7 + tmp_121 * tmp_4 ) ) ) +
          0.020848748529055869 * tmp_49 *
              ( ( ( tmp_1 * tmp_125 + tmp_126 * tmp_2 + tmp_127 * tmp_6 ) *
                  ( tmp_1 * tmp_125 + tmp_126 * tmp_2 + tmp_127 * tmp_6 ) ) +
                ( ( tmp_11 * tmp_127 + tmp_125 * tmp_13 + tmp_126 * tmp_15 ) *
                  ( tmp_11 * tmp_127 + tmp_125 * tmp_13 + tmp_126 * tmp_15 ) ) +
                ( ( tmp_125 * tmp_14 + tmp_126 * tmp_7 + tmp_127 * tmp_4 ) *
                  ( tmp_125 * tmp_14 + tmp_126 * tmp_7 + tmp_127 * tmp_4 ) ) ) +
          0.0068572537431980923 * tmp_49 *
              ( ( ( tmp_1 * tmp_131 + tmp_132 * tmp_2 + tmp_133 * tmp_6 ) *
                  ( tmp_1 * tmp_131 + tmp_132 * tmp_2 + tmp_133 * tmp_6 ) ) +
                ( ( tmp_11 * tmp_133 + tmp_13 * tmp_131 + tmp_132 * tmp_15 ) *
                  ( tmp_11 * tmp_133 + tmp_13 * tmp_131 + tmp_132 * tmp_15 ) ) +
                ( ( tmp_131 * tmp_14 + tmp_132 * tmp_7 + tmp_133 * tmp_4 ) *
                  ( tmp_131 * tmp_14 + tmp_132 * tmp_7 + tmp_133 * tmp_4 ) ) ) +
          0.037198804536718075 * tmp_49 *
              ( ( ( tmp_1 * tmp_137 + tmp_138 * tmp_2 + tmp_139 * tmp_6 ) *
                  ( tmp_1 * tmp_137 + tmp_138 * tmp_2 + tmp_139 * tmp_6 ) ) +
                ( ( tmp_11 * tmp_139 + tmp_13 * tmp_137 + tmp_138 * tmp_15 ) *
                  ( tmp_11 * tmp_139 + tmp_13 * tmp_137 + tmp_138 * tmp_15 ) ) +
                ( ( tmp_137 * tmp_14 + tmp_138 * tmp_7 + tmp_139 * tmp_4 ) *
                  ( tmp_137 * tmp_14 + tmp_138 * tmp_7 + tmp_139 * tmp_4 ) ) ) +
          0.042507265838595799 * tmp_49 *
              ( ( ( tmp_1 * tmp_143 + tmp_144 * tmp_2 + tmp_145 * tmp_6 ) *
                  ( tmp_1 * tmp_143 + tmp_144 * tmp_2 + tmp_145 * tmp_6 ) ) +
                ( ( tmp_11 * tmp_145 + tmp_13 * tmp_143 + tmp_144 * tmp_15 ) *
                  ( tmp_11 * tmp_145 + tmp_13 * tmp_143 + tmp_144 * tmp_15 ) ) +
                ( ( tmp_14 * tmp_143 + tmp_144 * tmp_7 + tmp_145 * tmp_4 ) *
                  ( tmp_14 * tmp_143 + tmp_144 * tmp_7 + tmp_145 * tmp_4 ) ) ) +
          0.0068572537431980923 * tmp_49 *
              ( ( ( tmp_1 * tmp_149 + tmp_150 * tmp_2 + tmp_151 * tmp_6 ) *
                  ( tmp_1 * tmp_149 + tmp_150 * tmp_2 + tmp_151 * tmp_6 ) ) +
                ( ( tmp_11 * tmp_151 + tmp_13 * tmp_149 + tmp_15 * tmp_150 ) *
                  ( tmp_11 * tmp_151 + tmp_13 * tmp_149 + tmp_15 * tmp_150 ) ) +
                ( ( tmp_14 * tmp_149 + tmp_150 * tmp_7 + tmp_151 * tmp_4 ) *
                  ( tmp_14 * tmp_149 + tmp_150 * tmp_7 + tmp_151 * tmp_4 ) ) ) +
          0.037198804536718075 * tmp_49 *
              ( ( ( tmp_1 * tmp_155 + tmp_156 * tmp_2 + tmp_157 * tmp_6 ) *
                  ( tmp_1 * tmp_155 + tmp_156 * tmp_2 + tmp_157 * tmp_6 ) ) +
                ( ( tmp_11 * tmp_157 + tmp_13 * tmp_155 + tmp_15 * tmp_156 ) *
                  ( tmp_11 * tmp_157 + tmp_13 * tmp_155 + tmp_15 * tmp_156 ) ) +
                ( ( tmp_14 * tmp_155 + tmp_156 * tmp_7 + tmp_157 * tmp_4 ) *
                  ( tmp_14 * tmp_155 + tmp_156 * tmp_7 + tmp_157 * tmp_4 ) ) ) +
          0.042507265838595799 * tmp_49 *
              ( ( ( tmp_1 * tmp_161 + tmp_162 * tmp_2 + tmp_163 * tmp_6 ) *
                  ( tmp_1 * tmp_161 + tmp_162 * tmp_2 + tmp_163 * tmp_6 ) ) +
                ( ( tmp_11 * tmp_163 + tmp_13 * tmp_161 + tmp_15 * tmp_162 ) *
                  ( tmp_11 * tmp_163 + tmp_13 * tmp_161 + tmp_15 * tmp_162 ) ) +
                ( ( tmp_14 * tmp_161 + tmp_162 * tmp_7 + tmp_163 * tmp_4 ) *
                  ( tmp_14 * tmp_161 + tmp_162 * tmp_7 + tmp_163 * tmp_4 ) ) ) +
          0.019202922745021479 * tmp_49 *
              ( ( ( tmp_1 * tmp_167 + tmp_168 * tmp_2 + tmp_169 * tmp_6 ) *
                  ( tmp_1 * tmp_167 + tmp_168 * tmp_2 + tmp_169 * tmp_6 ) ) +
                ( ( tmp_11 * tmp_169 + tmp_13 * tmp_167 + tmp_15 * tmp_168 ) *
                  ( tmp_11 * tmp_169 + tmp_13 * tmp_167 + tmp_15 * tmp_168 ) ) +
                ( ( tmp_14 * tmp_167 + tmp_168 * tmp_7 + tmp_169 * tmp_4 ) *
                  ( tmp_14 * tmp_167 + tmp_168 * tmp_7 + tmp_169 * tmp_4 ) ) ) +
          0.0068572537431980923 * tmp_49 *
              ( ( ( tmp_1 * tmp_37 + tmp_2 * tmp_41 + tmp_45 * tmp_6 ) * ( tmp_1 * tmp_37 + tmp_2 * tmp_41 + tmp_45 * tmp_6 ) ) +
                ( ( tmp_11 * tmp_45 + tmp_13 * tmp_37 + tmp_15 * tmp_41 ) *
                  ( tmp_11 * tmp_45 + tmp_13 * tmp_37 + tmp_15 * tmp_41 ) ) +
                ( ( tmp_14 * tmp_37 + tmp_4 * tmp_45 + tmp_41 * tmp_7 ) *
                  ( tmp_14 * tmp_37 + tmp_4 * tmp_45 + tmp_41 * tmp_7 ) ) ) +
          0.037198804536718075 * tmp_49 *
              ( ( ( tmp_1 * tmp_53 + tmp_2 * tmp_54 + tmp_55 * tmp_6 ) * ( tmp_1 * tmp_53 + tmp_2 * tmp_54 + tmp_55 * tmp_6 ) ) +
                ( ( tmp_11 * tmp_55 + tmp_13 * tmp_53 + tmp_15 * tmp_54 ) *
                  ( tmp_11 * tmp_55 + tmp_13 * tmp_53 + tmp_15 * tmp_54 ) ) +
                ( ( tmp_14 * tmp_53 + tmp_4 * tmp_55 + tmp_54 * tmp_7 ) *
                  ( tmp_14 * tmp_53 + tmp_4 * tmp_55 + tmp_54 * tmp_7 ) ) ) +
          0.020848748529055869 * tmp_49 *
              ( ( ( tmp_1 * tmp_59 + tmp_2 * tmp_60 + tmp_6 * tmp_61 ) * ( tmp_1 * tmp_59 + tmp_2 * tmp_60 + tmp_6 * tmp_61 ) ) +
                ( ( tmp_11 * tmp_61 + tmp_13 * tmp_59 + tmp_15 * tmp_60 ) *
                  ( tmp_11 * tmp_61 + tmp_13 * tmp_59 + tmp_15 * tmp_60 ) ) +
                ( ( tmp_14 * tmp_59 + tmp_4 * tmp_61 + tmp_60 * tmp_7 ) *
                  ( tmp_14 * tmp_59 + tmp_4 * tmp_61 + tmp_60 * tmp_7 ) ) ) +
          0.019202922745021479 * tmp_49 *
              ( ( ( tmp_1 * tmp_65 + tmp_2 * tmp_66 + tmp_6 * tmp_67 ) * ( tmp_1 * tmp_65 + tmp_2 * tmp_66 + tmp_6 * tmp_67 ) ) +
                ( ( tmp_11 * tmp_67 + tmp_13 * tmp_65 + tmp_15 * tmp_66 ) *
                  ( tmp_11 * tmp_67 + tmp_13 * tmp_65 + tmp_15 * tmp_66 ) ) +
                ( ( tmp_14 * tmp_65 + tmp_4 * tmp_67 + tmp_66 * tmp_7 ) *
                  ( tmp_14 * tmp_65 + tmp_4 * tmp_67 + tmp_66 * tmp_7 ) ) ) +
          0.020848748529055869 * tmp_49 *
              ( ( ( tmp_1 * tmp_71 + tmp_2 * tmp_72 + tmp_6 * tmp_73 ) * ( tmp_1 * tmp_71 + tmp_2 * tmp_72 + tmp_6 * tmp_73 ) ) +
                ( ( tmp_11 * tmp_73 + tmp_13 * tmp_71 + tmp_15 * tmp_72 ) *
                  ( tmp_11 * tmp_73 + tmp_13 * tmp_71 + tmp_15 * tmp_72 ) ) +
                ( ( tmp_14 * tmp_71 + tmp_4 * tmp_73 + tmp_7 * tmp_72 ) *
                  ( tmp_14 * tmp_71 + tmp_4 * tmp_73 + tmp_7 * tmp_72 ) ) ) +
          0.019202922745021479 * tmp_49 *
              ( ( ( tmp_1 * tmp_77 + tmp_2 * tmp_78 + tmp_6 * tmp_79 ) * ( tmp_1 * tmp_77 + tmp_2 * tmp_78 + tmp_6 * tmp_79 ) ) +
                ( ( tmp_11 * tmp_79 + tmp_13 * tmp_77 + tmp_15 * tmp_78 ) *
                  ( tmp_11 * tmp_79 + tmp_13 * tmp_77 + tmp_15 * tmp_78 ) ) +
                ( ( tmp_14 * tmp_77 + tmp_4 * tmp_79 + tmp_7 * tmp_78 ) *
                  ( tmp_14 * tmp_77 + tmp_4 * tmp_79 + tmp_7 * tmp_78 ) ) ) +
          0.020848748529055869 * tmp_49 *
              ( ( ( tmp_1 * tmp_83 + tmp_2 * tmp_84 + tmp_6 * tmp_85 ) * ( tmp_1 * tmp_83 + tmp_2 * tmp_84 + tmp_6 * tmp_85 ) ) +
                ( ( tmp_11 * tmp_85 + tmp_13 * tmp_83 + tmp_15 * tmp_84 ) *
                  ( tmp_11 * tmp_85 + tmp_13 * tmp_83 + tmp_15 * tmp_84 ) ) +
                ( ( tmp_14 * tmp_83 + tmp_4 * tmp_85 + tmp_7 * tmp_84 ) *
                  ( tmp_14 * tmp_83 + tmp_4 * tmp_85 + tmp_7 * tmp_84 ) ) ) +
          0.019202922745021479 * tmp_49 *
              ( ( ( tmp_1 * tmp_89 + tmp_2 * tmp_90 + tmp_6 * tmp_91 ) * ( tmp_1 * tmp_89 + tmp_2 * tmp_90 + tmp_6 * tmp_91 ) ) +
                ( ( tmp_11 * tmp_91 + tmp_13 * tmp_89 + tmp_15 * tmp_90 ) *
                  ( tmp_11 * tmp_91 + tmp_13 * tmp_89 + tmp_15 * tmp_90 ) ) +
                ( ( tmp_14 * tmp_89 + tmp_4 * tmp_91 + tmp_7 * tmp_90 ) *
                  ( tmp_14 * tmp_89 + tmp_4 * tmp_91 + tmp_7 * tmp_90 ) ) ) +
          0.020848748529055869 * tmp_49 *
              ( ( ( tmp_1 * tmp_95 + tmp_2 * tmp_96 + tmp_6 * tmp_97 ) * ( tmp_1 * tmp_95 + tmp_2 * tmp_96 + tmp_6 * tmp_97 ) ) +
                ( ( tmp_11 * tmp_97 + tmp_13 * tmp_95 + tmp_15 * tmp_96 ) *
                  ( tmp_11 * tmp_97 + tmp_13 * tmp_95 + tmp_15 * tmp_96 ) ) +
                ( ( tmp_14 * tmp_95 + tmp_4 * tmp_97 + tmp_7 * tmp_96 ) *
                  ( tmp_14 * tmp_95 + tmp_4 * tmp_97 + tmp_7 * tmp_96 ) ) );
      elMat( 0, 0 ) = a_0_0;
   }
};

} // namespace eg
} // namespace dg
} // namespace hyteg

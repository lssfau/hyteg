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
#include "hyteg/types/matrix.hpp"
#include "hyteg/types/pointnd.hpp"

#include "Eigen/Eigen"

namespace hyteg {
namespace dg {

class DGDiffusionForm_Example : public DGForm
{
 public:
   void integrateVolume( const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >&    coords,
                         const DGBasisInfo&                                       trialBasis,
                         const DGBasisInfo&                                       testBasis,
                         int                                                      trialDegree,
                         int                                                      testDegree,
                         Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const override
   {
      elMat.resize( testBasis.numDoFsPerElement( testDegree ), trialBasis.numDoFsPerElement( trialDegree ) );

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      real_t tmp_0  = p_affine_0_1 * p_affine_1_0;
      real_t tmp_1  = p_affine_0_0 * p_affine_1_1;
      real_t tmp_2  = 2 * tmp_1;
      real_t tmp_3  = 2 * p_affine_0_0;
      real_t tmp_4  = p_affine_1_0 * p_affine_2_1;
      real_t tmp_5  = p_affine_0_1 * tmp_4;
      real_t tmp_6  = p_affine_0_1 * p_affine_2_0;
      real_t tmp_7  = p_affine_0_0 * p_affine_2_1;
      real_t tmp_8  = 2 * tmp_6;
      real_t tmp_9  = p_affine_2_0 * p_affine_2_1;
      real_t tmp_10 = p_affine_1_0 * p_affine_1_1;
      real_t tmp_11 = p_affine_1_1 * p_affine_2_0;
      real_t tmp_12 = tmp_11 * tmp_4;
      real_t tmp_13 = ( p_affine_0_0 * p_affine_0_0 );
      real_t tmp_14 = ( p_affine_1_1 * p_affine_1_1 );
      real_t tmp_15 = tmp_13 * tmp_14;
      real_t tmp_16 = ( p_affine_2_1 * p_affine_2_1 );
      real_t tmp_17 = tmp_13 * tmp_16;
      real_t tmp_18 = ( p_affine_0_1 * p_affine_0_1 );
      real_t tmp_19 = ( p_affine_1_0 * p_affine_1_0 );
      real_t tmp_20 = tmp_18 * tmp_19;
      real_t tmp_21 = ( p_affine_2_0 * p_affine_2_0 );
      real_t tmp_22 = tmp_18 * tmp_21;
      real_t tmp_23 = tmp_16 * tmp_19;
      real_t tmp_24 = tmp_14 * tmp_21;
      real_t tmp_25 = p_affine_1_0 * tmp_16;
      real_t tmp_26 = p_affine_2_0 * tmp_14;
      real_t tmp_27 = p_affine_1_1 * p_affine_2_1;
      real_t tmp_28 = tmp_13 * tmp_27;
      real_t tmp_29 = 2 * p_affine_0_1;
      real_t tmp_30 = p_affine_2_1 * tmp_19;
      real_t tmp_31 = p_affine_1_1 * tmp_21;
      real_t tmp_32 = p_affine_1_0 * p_affine_2_0;
      real_t tmp_33 = tmp_18 * tmp_32;
      real_t tmp_34 = std::abs( tmp_0 - tmp_1 + tmp_11 - tmp_4 - tmp_6 + tmp_7 );
      real_t tmp_35 =
          tmp_34 / ( -tmp_0 * tmp_2 + tmp_10 * tmp_8 - 2 * tmp_12 + tmp_15 + tmp_17 + tmp_2 * tmp_4 + tmp_2 * tmp_6 +
                     tmp_2 * tmp_9 + tmp_20 + tmp_22 + tmp_23 + tmp_24 - tmp_25 * tmp_3 - tmp_26 * tmp_3 - 2 * tmp_28 -
                     tmp_29 * tmp_30 - tmp_29 * tmp_31 + tmp_3 * tmp_5 - 2 * tmp_33 + tmp_4 * tmp_8 - tmp_7 * tmp_8 );
      real_t tmp_36 = tmp_32 * tmp_35;
      real_t tmp_37 = tmp_27 * tmp_35;
      real_t tmp_38 = 4 * tmp_1;
      real_t tmp_39 = 4 * p_affine_0_0;
      real_t tmp_40 = 4 * tmp_6;
      real_t tmp_41 = 4 * p_affine_0_1;
      real_t tmp_42 = tmp_34 / ( -tmp_0 * tmp_38 + tmp_10 * tmp_40 - 4 * tmp_12 + 2 * tmp_15 + 2 * tmp_17 + 2 * tmp_20 +
                                 2 * tmp_22 + 2 * tmp_23 + 2 * tmp_24 - tmp_25 * tmp_39 - tmp_26 * tmp_39 - 4 * tmp_28 -
                                 tmp_30 * tmp_41 - tmp_31 * tmp_41 - 4 * tmp_33 + tmp_38 * tmp_4 + tmp_38 * tmp_6 +
                                 tmp_38 * tmp_9 + tmp_39 * tmp_5 + tmp_4 * tmp_40 - tmp_40 * tmp_7 );
      real_t tmp_43 = tmp_32 * tmp_42;
      real_t tmp_44 = tmp_27 * tmp_42;
      real_t tmp_45 = tmp_21 * tmp_35;
      real_t tmp_46 = tmp_16 * tmp_35;
      real_t tmp_47 = tmp_21 * tmp_42;
      real_t tmp_48 = tmp_16 * tmp_42;
      real_t tmp_49 = tmp_45 + tmp_46 - tmp_47 - tmp_48;
      real_t tmp_50 = tmp_19 * tmp_35;
      real_t tmp_51 = tmp_14 * tmp_35;
      real_t tmp_52 = tmp_19 * tmp_42;
      real_t tmp_53 = tmp_14 * tmp_42;
      real_t tmp_54 = tmp_50 + tmp_51 - tmp_52 - tmp_53;
      real_t tmp_55 = p_affine_0_0 * p_affine_1_0;
      real_t tmp_56 = tmp_42 * tmp_55;
      real_t tmp_57 = p_affine_0_1 * p_affine_1_1;
      real_t tmp_58 = tmp_42 * tmp_57;
      real_t tmp_59 = tmp_35 * tmp_55;
      real_t tmp_60 = tmp_35 * tmp_57;
      real_t tmp_61 = tmp_36 + tmp_37 - tmp_43 - tmp_44;
      real_t tmp_62 = p_affine_0_0 * p_affine_2_0;
      real_t tmp_63 = tmp_35 * tmp_62;
      real_t tmp_64 = p_affine_0_1 * p_affine_2_1;
      real_t tmp_65 = tmp_35 * tmp_64;
      real_t tmp_66 = tmp_42 * tmp_62;
      real_t tmp_67 = tmp_42 * tmp_64;
      real_t tmp_68 = tmp_63 + tmp_65 - tmp_66 - tmp_67;
      real_t tmp_69 = -tmp_45 - tmp_46 + tmp_47 + tmp_48 + tmp_56 + tmp_58 - tmp_59 - tmp_60 + tmp_61 + tmp_68;
      real_t tmp_70 = -tmp_56 - tmp_58 + tmp_59 + tmp_60;
      real_t tmp_71 = -tmp_50 - tmp_51 + tmp_52 + tmp_53 + tmp_61 - tmp_63 - tmp_65 + tmp_66 + tmp_67 + tmp_70;
      real_t tmp_72 = p_affine_2_0 * tmp_3;
      real_t tmp_73 = p_affine_2_1 * tmp_29;
      real_t tmp_74 = tmp_13 * tmp_35;
      real_t tmp_75 = tmp_18 * tmp_35;
      real_t tmp_76 = tmp_13 * tmp_42;
      real_t tmp_77 = tmp_18 * tmp_42;
      real_t tmp_78 = tmp_74 + tmp_75 - tmp_76 - tmp_77;
      real_t tmp_79 = -tmp_36 - tmp_37 + tmp_43 + tmp_44 + tmp_68 + tmp_70 - tmp_74 - tmp_75 + tmp_76 + tmp_77;
      real_t tmp_80 = p_affine_1_0 * tmp_3;
      real_t tmp_81 = p_affine_1_1 * tmp_29;
      real_t a_0_0  = -2 * tmp_36 - 2 * tmp_37 + 2 * tmp_43 + 2 * tmp_44 + tmp_49 + tmp_54;
      real_t a_0_1  = tmp_69;
      real_t a_0_2  = tmp_71;
      real_t a_1_0  = tmp_69;
      real_t a_1_1  = -tmp_35 * tmp_72 - tmp_35 * tmp_73 + tmp_42 * tmp_72 + tmp_42 * tmp_73 + tmp_49 + tmp_78;
      real_t a_1_2  = tmp_79;
      real_t a_2_0  = tmp_71;
      real_t a_2_1  = tmp_79;
      real_t a_2_2  = -tmp_35 * tmp_80 - tmp_35 * tmp_81 + tmp_42 * tmp_80 + tmp_42 * tmp_81 + tmp_54 + tmp_78;

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

   virtual void integrateFacetInner( const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >&    coordsElement,
                                     const std::array< Eigen::Matrix< real_t, 2, 1 >, 2 >&    coordsFacet,
                                     const Eigen::Matrix< real_t, 2, 1 >&                     outwardNormal,
                                     const DGBasisInfo&                                       trialBasis,
                                     const DGBasisInfo&                                       testBasis,
                                     int                                                      trialDegree,
                                     int                                                      testDegree,
                                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
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

      const auto p_affine_8_0 = outwardNormal( 0 );
      const auto p_affine_8_1 = outwardNormal( 1 );

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
      real_t tmp_10 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_11 = p_affine_6_0 + tmp_4;
      real_t tmp_12 = 0.21132486540518713 * tmp_10 + tmp_11;
      real_t tmp_13 = tmp_7 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_14 = tmp_12 * tmp_13;
      real_t tmp_15 = std::pow( ( tmp_0 * tmp_0 ) + ( tmp_10 * tmp_10 ), 1.0 / 2.0 );
      real_t tmp_16 = tmp_14 + tmp_9;
      real_t tmp_17 = tmp_15 * tmp_16;
      real_t tmp_18 = tmp_7 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_19 = tmp_18 * ( -p_affine_8_1 * tmp_17 + tmp_3 );
      real_t tmp_20 = tmp_6 * tmp_7;
      real_t tmp_21 = tmp_20 * ( -p_affine_8_0 * tmp_17 + tmp_12 );
      real_t tmp_22 = -tmp_14 - tmp_19 - tmp_21 - tmp_9 + 1;
      real_t tmp_23 = p_affine_8_0 * ( -tmp_18 - tmp_20 ) + p_affine_8_1 * ( -tmp_13 - tmp_8 );
      real_t tmp_24 = 1.0 * tmp_23;
      real_t tmp_25 = 0.78867513459481287 * tmp_0 + tmp_2;
      real_t tmp_26 = tmp_25 * tmp_8;
      real_t tmp_27 = 0.78867513459481287 * tmp_10 + tmp_11;
      real_t tmp_28 = tmp_13 * tmp_27;
      real_t tmp_29 = tmp_26 + tmp_28;
      real_t tmp_30 = tmp_15 * tmp_29;
      real_t tmp_31 = tmp_18 * ( -p_affine_8_1 * tmp_30 + tmp_25 );
      real_t tmp_32 = tmp_20 * ( -p_affine_8_0 * tmp_30 + tmp_27 );
      real_t tmp_33 = -tmp_26 - tmp_28 - tmp_31 - tmp_32 + 1;
      real_t tmp_34 = tmp_19 + tmp_21;
      real_t tmp_35 = 0.5 * tmp_23;
      real_t tmp_36 = p_affine_8_0 * tmp_20 + p_affine_8_1 * tmp_13;
      real_t tmp_37 = 0.5 * tmp_36;
      real_t tmp_38 = 6 * tmp_22;
      real_t tmp_39 = tmp_31 + tmp_32;
      real_t tmp_40 = 6 * tmp_33;
      real_t tmp_41 = -0.5 * tmp_22 * tmp_37 - 0.5 * tmp_33 * tmp_37 - 0.5 * tmp_34 * tmp_35 + 0.5 * tmp_34 * tmp_38 -
                      0.5 * tmp_35 * tmp_39 + 0.5 * tmp_39 * tmp_40;
      real_t tmp_42 = p_affine_8_0 * tmp_18 + p_affine_8_1 * tmp_8;
      real_t tmp_43 = 0.5 * tmp_42;
      real_t tmp_44 = -0.5 * tmp_16 * tmp_35 + 0.5 * tmp_16 * tmp_38 - 0.5 * tmp_22 * tmp_43 - 0.5 * tmp_29 * tmp_35 +
                      0.5 * tmp_29 * tmp_40 - 0.5 * tmp_33 * tmp_43;
      real_t tmp_45 = 1.0 * tmp_36;
      real_t tmp_46 = 3.0 * tmp_16 * tmp_34 - 0.5 * tmp_16 * tmp_37 - 0.5 * tmp_29 * tmp_37 + 3.0 * tmp_29 * tmp_39 -
                      0.5 * tmp_34 * tmp_43 - 0.5 * tmp_39 * tmp_43;
      real_t tmp_47 = 1.0 * tmp_42;
      real_t a_0_0  = 3.0 * ( tmp_22 * tmp_22 ) - 0.5 * tmp_22 * tmp_24 - 0.5 * tmp_24 * tmp_33 + 3.0 * ( tmp_33 * tmp_33 );
      real_t a_0_1  = tmp_41;
      real_t a_0_2  = tmp_44;
      real_t a_1_0  = tmp_41;
      real_t a_1_1  = 3.0 * ( tmp_34 * tmp_34 ) - 0.5 * tmp_34 * tmp_45 + 3.0 * ( tmp_39 * tmp_39 ) - 0.5 * tmp_39 * tmp_45;
      real_t a_1_2  = tmp_46;
      real_t a_2_0  = tmp_44;
      real_t a_2_1  = tmp_46;
      real_t a_2_2  = 3.0 * ( tmp_16 * tmp_16 ) - 0.5 * tmp_16 * tmp_47 + 3.0 * ( tmp_29 * tmp_29 ) - 0.5 * tmp_29 * tmp_47;

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

   virtual void integrateFacetCoupling( const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >&    coordsElementInner,
                                        const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >&    coordsElementOuter,
                                        const std::array< Eigen::Matrix< real_t, 2, 1 >, 2 >&    coordsFacet,
                                        const Eigen::Matrix< real_t, 2, 1 >&                     outwardNormal,
                                        const DGBasisInfo&                                       trialBasis,
                                        const DGBasisInfo&                                       testBasis,
                                        int                                                      trialDegree,
                                        int                                                      testDegree,
                                        Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      elMat.resize( testBasis.numDoFsPerElement( testDegree ), trialBasis.numDoFsPerElement( trialDegree ) );

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

      const auto p_affine_8_0 = outwardNormal( 0 );
      const auto p_affine_8_1 = outwardNormal( 1 );

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
      real_t tmp_10 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_11 = p_affine_6_0 + 0.21132486540518713 * tmp_10;
      real_t tmp_12 = tmp_11 + tmp_4;
      real_t tmp_13 = tmp_7 * ( p_affine_3_1 - p_affine_4_1 );
      real_t tmp_14 = tmp_12 * tmp_13;
      real_t tmp_15 = std::pow( ( tmp_1 * tmp_1 ) + ( tmp_10 * tmp_10 ), 1.0 / 2.0 );
      real_t tmp_16 = tmp_14 + tmp_9;
      real_t tmp_17 = tmp_15 * tmp_16;
      real_t tmp_18 = tmp_7 * ( p_affine_3_0 - p_affine_5_0 );
      real_t tmp_19 = tmp_18 * ( p_affine_8_1 * tmp_17 + tmp_3 );
      real_t tmp_20 = tmp_6 * tmp_7;
      real_t tmp_21 = tmp_20 * ( p_affine_8_0 * tmp_17 + tmp_12 );
      real_t tmp_22 = -tmp_14 - tmp_19 - tmp_21 - tmp_9 + 1;
      real_t tmp_23 = -p_affine_0_0;
      real_t tmp_24 = p_affine_1_0 + tmp_23;
      real_t tmp_25 = -p_affine_0_1;
      real_t tmp_26 = p_affine_2_1 + tmp_25;
      real_t tmp_27 = 1.0 / ( tmp_24 * tmp_26 - ( p_affine_1_1 + tmp_25 ) * ( p_affine_2_0 + tmp_23 ) );
      real_t tmp_28 = tmp_27 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_29 = tmp_26 * tmp_27;
      real_t tmp_30 = tmp_24 * tmp_27;
      real_t tmp_31 = tmp_27 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_32 = 0.5 * p_affine_8_0 * ( -tmp_28 - tmp_29 ) + 0.5 * p_affine_8_1 * ( -tmp_30 - tmp_31 );
      real_t tmp_33 = tmp_2 + tmp_25;
      real_t tmp_34 = tmp_30 * tmp_33;
      real_t tmp_35 = tmp_11 + tmp_23;
      real_t tmp_36 = tmp_31 * tmp_35;
      real_t tmp_37 = tmp_34 + tmp_36;
      real_t tmp_38 = tmp_15 * tmp_37;
      real_t tmp_39 = tmp_28 * ( -p_affine_8_1 * tmp_38 + tmp_33 );
      real_t tmp_40 = tmp_29 * ( -p_affine_8_0 * tmp_38 + tmp_35 );
      real_t tmp_41 = -tmp_34 - tmp_36 - tmp_39 - tmp_40 + 1;
      real_t tmp_42 = 0.5 * p_affine_8_0 * ( -tmp_18 - tmp_20 ) + 0.5 * p_affine_8_1 * ( -tmp_13 - tmp_8 );
      real_t tmp_43 = 6 * tmp_41;
      real_t tmp_44 = p_affine_6_1 + 0.78867513459481287 * tmp_1;
      real_t tmp_45 = tmp_0 + tmp_44;
      real_t tmp_46 = tmp_45 * tmp_8;
      real_t tmp_47 = p_affine_6_0 + 0.78867513459481287 * tmp_10;
      real_t tmp_48 = tmp_4 + tmp_47;
      real_t tmp_49 = tmp_13 * tmp_48;
      real_t tmp_50 = tmp_46 + tmp_49;
      real_t tmp_51 = tmp_15 * tmp_50;
      real_t tmp_52 = tmp_18 * ( p_affine_8_1 * tmp_51 + tmp_45 );
      real_t tmp_53 = tmp_20 * ( p_affine_8_0 * tmp_51 + tmp_48 );
      real_t tmp_54 = -tmp_46 - tmp_49 - tmp_52 - tmp_53 + 1;
      real_t tmp_55 = tmp_25 + tmp_44;
      real_t tmp_56 = tmp_30 * tmp_55;
      real_t tmp_57 = tmp_23 + tmp_47;
      real_t tmp_58 = tmp_31 * tmp_57;
      real_t tmp_59 = tmp_56 + tmp_58;
      real_t tmp_60 = tmp_15 * tmp_59;
      real_t tmp_61 = tmp_28 * ( -p_affine_8_1 * tmp_60 + tmp_55 );
      real_t tmp_62 = tmp_29 * ( -p_affine_8_0 * tmp_60 + tmp_57 );
      real_t tmp_63 = -tmp_56 - tmp_58 - tmp_61 - tmp_62 + 1;
      real_t tmp_64 = 6 * tmp_63;
      real_t tmp_65 = tmp_19 + tmp_21;
      real_t tmp_66 = 0.5 * p_affine_8_0 * tmp_20 + 0.5 * p_affine_8_1 * tmp_13;
      real_t tmp_67 = tmp_52 + tmp_53;
      real_t tmp_68 = 0.5 * p_affine_8_0 * tmp_18 + 0.5 * p_affine_8_1 * tmp_8;
      real_t tmp_69 = tmp_39 + tmp_40;
      real_t tmp_70 = 0.5 * p_affine_8_0 * tmp_29 + 0.5 * p_affine_8_1 * tmp_31;
      real_t tmp_71 = 6 * tmp_69;
      real_t tmp_72 = tmp_61 + tmp_62;
      real_t tmp_73 = 6 * tmp_72;
      real_t tmp_74 = 0.5 * p_affine_8_0 * tmp_28 + 0.5 * p_affine_8_1 * tmp_30;
      real_t tmp_75 = 6 * tmp_37;
      real_t tmp_76 = 6 * tmp_59;
      real_t a_0_0  = -0.5 * tmp_22 * tmp_32 - 0.5 * tmp_22 * tmp_43 - 0.5 * tmp_32 * tmp_54 + 0.5 * tmp_41 * tmp_42 +
                     0.5 * tmp_42 * tmp_63 - 0.5 * tmp_54 * tmp_64;
      real_t a_0_1 = -0.5 * tmp_32 * tmp_65 - 0.5 * tmp_32 * tmp_67 + 0.5 * tmp_41 * tmp_66 - 0.5 * tmp_43 * tmp_65 +
                     0.5 * tmp_63 * tmp_66 - 0.5 * tmp_64 * tmp_67;
      real_t a_0_2 = -0.5 * tmp_16 * tmp_32 - 0.5 * tmp_16 * tmp_43 - 0.5 * tmp_32 * tmp_50 + 0.5 * tmp_41 * tmp_68 -
                     0.5 * tmp_50 * tmp_64 + 0.5 * tmp_63 * tmp_68;
      real_t a_1_0 = -0.5 * tmp_22 * tmp_70 - 0.5 * tmp_22 * tmp_71 + 0.5 * tmp_42 * tmp_69 + 0.5 * tmp_42 * tmp_72 -
                     0.5 * tmp_54 * tmp_70 - 0.5 * tmp_54 * tmp_73;
      real_t a_1_1 = -0.5 * tmp_65 * tmp_70 - 0.5 * tmp_65 * tmp_71 + 0.5 * tmp_66 * tmp_69 + 0.5 * tmp_66 * tmp_72 -
                     0.5 * tmp_67 * tmp_70 - 0.5 * tmp_67 * tmp_73;
      real_t a_1_2 = -0.5 * tmp_16 * tmp_70 - 0.5 * tmp_16 * tmp_71 - 0.5 * tmp_50 * tmp_70 - 0.5 * tmp_50 * tmp_73 +
                     0.5 * tmp_68 * tmp_69 + 0.5 * tmp_68 * tmp_72;
      real_t a_2_0 = -0.5 * tmp_22 * tmp_74 - 0.5 * tmp_22 * tmp_75 + 0.5 * tmp_37 * tmp_42 + 0.5 * tmp_42 * tmp_59 -
                     0.5 * tmp_54 * tmp_74 - 0.5 * tmp_54 * tmp_76;
      real_t a_2_1 = 0.5 * tmp_37 * tmp_66 + 0.5 * tmp_59 * tmp_66 - 0.5 * tmp_65 * tmp_74 - 0.5 * tmp_65 * tmp_75 -
                     0.5 * tmp_67 * tmp_74 - 0.5 * tmp_67 * tmp_76;
      real_t a_2_2 = -0.5 * tmp_16 * tmp_74 - 0.5 * tmp_16 * tmp_75 + 0.5 * tmp_37 * tmp_68 - 0.5 * tmp_50 * tmp_74 -
                     0.5 * tmp_50 * tmp_76 + 0.5 * tmp_59 * tmp_68;

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

   virtual void integrateFacetDirichletBoundary( const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >&    coordsElement,
                                                 const std::array< Eigen::Matrix< real_t, 2, 1 >, 2 >&    coordsFacet,
                                                 const Eigen::Matrix< real_t, 2, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
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

      const auto p_affine_8_0 = outwardNormal( 0 );
      const auto p_affine_8_1 = outwardNormal( 1 );

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
      real_t tmp_10 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_11 = p_affine_6_0 + tmp_4;
      real_t tmp_12 = 0.21132486540518713 * tmp_10 + tmp_11;
      real_t tmp_13 = tmp_7 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_14 = tmp_12 * tmp_13;
      real_t tmp_15 = tmp_14 + tmp_9;
      real_t tmp_16 = std::pow( ( tmp_0 * tmp_0 ) + ( tmp_10 * tmp_10 ), 1.0 / 2.0 );
      real_t tmp_17 = tmp_15 * tmp_16;
      real_t tmp_18 = tmp_7 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_19 = tmp_18 * ( -p_affine_8_1 * tmp_17 + tmp_3 );
      real_t tmp_20 = tmp_6 * tmp_7;
      real_t tmp_21 = tmp_20 * ( -p_affine_8_0 * tmp_17 + tmp_12 );
      real_t tmp_22 = -tmp_14 - tmp_19 - tmp_21 - tmp_9 + 1;
      real_t tmp_23 = p_affine_8_0 * ( -tmp_18 - tmp_20 ) + p_affine_8_1 * ( -tmp_13 - tmp_8 );
      real_t tmp_24 = 2 * tmp_23;
      real_t tmp_25 = 0.78867513459481287 * tmp_0 + tmp_2;
      real_t tmp_26 = tmp_25 * tmp_8;
      real_t tmp_27 = 0.78867513459481287 * tmp_10 + tmp_11;
      real_t tmp_28 = tmp_13 * tmp_27;
      real_t tmp_29 = tmp_26 + tmp_28;
      real_t tmp_30 = tmp_16 * tmp_29;
      real_t tmp_31 = tmp_18 * ( -p_affine_8_1 * tmp_30 + tmp_25 );
      real_t tmp_32 = tmp_20 * ( -p_affine_8_0 * tmp_30 + tmp_27 );
      real_t tmp_33 = -tmp_26 - tmp_28 - tmp_31 - tmp_32 + 1;
      real_t tmp_34 = tmp_19 + tmp_21;
      real_t tmp_35 = p_affine_8_0 * tmp_20 + p_affine_8_1 * tmp_13;
      real_t tmp_36 = 6 * tmp_22;
      real_t tmp_37 = tmp_31 + tmp_32;
      real_t tmp_38 = 6 * tmp_33;
      real_t tmp_39 = -0.5 * tmp_22 * tmp_35 - 0.5 * tmp_23 * tmp_34 - 0.5 * tmp_23 * tmp_37 - 0.5 * tmp_33 * tmp_35 +
                      0.5 * tmp_34 * tmp_36 + 0.5 * tmp_37 * tmp_38;
      real_t tmp_40 = p_affine_8_0 * tmp_18 + p_affine_8_1 * tmp_8;
      real_t tmp_41 = -0.5 * tmp_15 * tmp_23 + 0.5 * tmp_15 * tmp_36 - 0.5 * tmp_22 * tmp_40 - 0.5 * tmp_23 * tmp_29 +
                      0.5 * tmp_29 * tmp_38 - 0.5 * tmp_33 * tmp_40;
      real_t tmp_42 = 2 * tmp_35;
      real_t tmp_43 = 3.0 * tmp_15 * tmp_34 - 0.5 * tmp_15 * tmp_35 - 0.5 * tmp_29 * tmp_35 + 3.0 * tmp_29 * tmp_37 -
                      0.5 * tmp_34 * tmp_40 - 0.5 * tmp_37 * tmp_40;
      real_t tmp_44 = 2 * tmp_40;
      real_t a_0_0  = 3.0 * ( tmp_22 * tmp_22 ) - 0.5 * tmp_22 * tmp_24 - 0.5 * tmp_24 * tmp_33 + 3.0 * ( tmp_33 * tmp_33 );
      real_t a_0_1  = tmp_39;
      real_t a_0_2  = tmp_41;
      real_t a_1_0  = tmp_39;
      real_t a_1_1  = 3.0 * ( tmp_34 * tmp_34 ) - 0.5 * tmp_34 * tmp_42 + 3.0 * ( tmp_37 * tmp_37 ) - 0.5 * tmp_37 * tmp_42;
      real_t a_1_2  = tmp_43;
      real_t a_2_0  = tmp_41;
      real_t a_2_1  = tmp_43;
      real_t a_2_2  = 3.0 * ( tmp_15 * tmp_15 ) - 0.5 * tmp_15 * tmp_44 + 3.0 * ( tmp_29 * tmp_29 ) - 0.5 * tmp_29 * tmp_44;

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
};

} // namespace dg
} // namespace hyteg
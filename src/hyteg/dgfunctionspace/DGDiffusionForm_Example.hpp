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

using walberla::real_c;

class DGDiffusionForm_Example : public DGForm2D
{
 public:
   DGDiffusionForm_Example()
   : callback_Scalar_Variable_Coefficient_2D_g( []( const Point3D& ) { return real_c( 0 ); } )
   {}

   DGDiffusionForm_Example( const std::function< real_t( const Point3D& ) >& _callback_Scalar_Variable_Coefficient_2D_g )
   : callback_Scalar_Variable_Coefficient_2D_g( _callback_Scalar_Variable_Coefficient_2D_g )
   {}

 protected:
   void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
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
      real_t tmp_3  = 0.1127016653792583 * tmp_0 + tmp_2;
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
      real_t tmp_14 = 0.1127016653792583 * tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6 * tmp_7;
      real_t tmp_16 = tmp_14 * tmp_15;
      real_t tmp_17 = tmp_7 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_18 = tmp_14 * tmp_17;
      real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20 = std::abs( std::pow( ( tmp_0 * tmp_0 ) + ( tmp_12 * tmp_12 ), 1.0 / 2.0 ) );
      real_t tmp_21 = 6 / tmp_20;
      real_t tmp_22 = 0.5 * p_affine_10_0;
      real_t tmp_23 = tmp_22 * ( -tmp_15 - tmp_17 );
      real_t tmp_24 = 0.5 * p_affine_10_1;
      real_t tmp_25 = tmp_24 * ( -tmp_10 - tmp_8 );
      real_t tmp_26 = -tmp_23 - tmp_25;
      real_t tmp_27 = tmp_23 + tmp_25;
      real_t tmp_28 = 0.27777777777777785 * tmp_20;
      real_t tmp_29 = 0.5 * tmp_0 + tmp_2;
      real_t tmp_30 = tmp_29 * tmp_8;
      real_t tmp_31 = tmp_10 * tmp_29;
      real_t tmp_32 = 0.5 * tmp_12 + tmp_13;
      real_t tmp_33 = tmp_15 * tmp_32;
      real_t tmp_34 = tmp_17 * tmp_32;
      real_t tmp_35 = -tmp_30 - tmp_31 - tmp_33 - tmp_34 + 1;
      real_t tmp_36 = 0.44444444444444442 * tmp_20;
      real_t tmp_37 = 0.8872983346207417 * tmp_0 + tmp_2;
      real_t tmp_38 = tmp_37 * tmp_8;
      real_t tmp_39 = tmp_10 * tmp_37;
      real_t tmp_40 = 0.8872983346207417 * tmp_12 + tmp_13;
      real_t tmp_41 = tmp_15 * tmp_40;
      real_t tmp_42 = tmp_17 * tmp_40;
      real_t tmp_43 = -tmp_38 - tmp_39 - tmp_41 - tmp_42 + 1;
      real_t tmp_44 = 0.27777777777777785 * tmp_20;
      real_t tmp_45 = tmp_11 + tmp_16;
      real_t tmp_46 = tmp_15 * tmp_22;
      real_t tmp_47 = tmp_10 * tmp_24;
      real_t tmp_48 = tmp_46 + tmp_47;
      real_t tmp_49 = tmp_19 * tmp_21;
      real_t tmp_50 = tmp_45 * tmp_49;
      real_t tmp_51 = tmp_31 + tmp_33;
      real_t tmp_52 = tmp_21 * tmp_35;
      real_t tmp_53 = tmp_51 * tmp_52;
      real_t tmp_54 = tmp_39 + tmp_41;
      real_t tmp_55 = tmp_21 * tmp_43;
      real_t tmp_56 = tmp_54 * tmp_55;
      real_t tmp_57 = tmp_18 + tmp_9;
      real_t tmp_58 = tmp_17 * tmp_22;
      real_t tmp_59 = tmp_24 * tmp_8;
      real_t tmp_60 = tmp_58 + tmp_59;
      real_t tmp_61 = tmp_49 * tmp_57;
      real_t tmp_62 = tmp_30 + tmp_34;
      real_t tmp_63 = tmp_52 * tmp_62;
      real_t tmp_64 = tmp_38 + tmp_42;
      real_t tmp_65 = tmp_55 * tmp_64;
      real_t tmp_66 = -tmp_46 - tmp_47;
      real_t tmp_67 = tmp_21 * tmp_45 * tmp_57;
      real_t tmp_68 = tmp_21 * tmp_51 * tmp_62;
      real_t tmp_69 = tmp_21 * tmp_54 * tmp_64;
      real_t tmp_70 = -tmp_58 - tmp_59;
      real_t a_0_0  = tmp_28 * ( ( tmp_19 * tmp_19 ) * tmp_21 + tmp_19 * tmp_26 - tmp_19 * tmp_27 ) +
                     tmp_36 * ( tmp_21 * ( tmp_35 * tmp_35 ) + tmp_26 * tmp_35 - tmp_27 * tmp_35 ) +
                     tmp_44 * ( tmp_21 * ( tmp_43 * tmp_43 ) + tmp_26 * tmp_43 - tmp_27 * tmp_43 );
      real_t a_0_1 = tmp_28 * ( -tmp_19 * tmp_48 + tmp_26 * tmp_45 + tmp_50 ) +
                     tmp_36 * ( tmp_26 * tmp_51 - tmp_35 * tmp_48 + tmp_53 ) +
                     tmp_44 * ( tmp_26 * tmp_54 - tmp_43 * tmp_48 + tmp_56 );
      real_t a_0_2 = tmp_28 * ( -tmp_19 * tmp_60 + tmp_26 * tmp_57 + tmp_61 ) +
                     tmp_36 * ( tmp_26 * tmp_62 - tmp_35 * tmp_60 + tmp_63 ) +
                     tmp_44 * ( tmp_26 * tmp_64 - tmp_43 * tmp_60 + tmp_65 );
      real_t a_1_0 = tmp_28 * ( tmp_19 * tmp_66 - tmp_27 * tmp_45 + tmp_50 ) +
                     tmp_36 * ( -tmp_27 * tmp_51 + tmp_35 * tmp_66 + tmp_53 ) +
                     tmp_44 * ( -tmp_27 * tmp_54 + tmp_43 * tmp_66 + tmp_56 );
      real_t a_1_1 = tmp_28 * ( tmp_21 * ( tmp_45 * tmp_45 ) - tmp_45 * tmp_48 + tmp_45 * tmp_66 ) +
                     tmp_36 * ( tmp_21 * ( tmp_51 * tmp_51 ) - tmp_48 * tmp_51 + tmp_51 * tmp_66 ) +
                     tmp_44 * ( tmp_21 * ( tmp_54 * tmp_54 ) - tmp_48 * tmp_54 + tmp_54 * tmp_66 );
      real_t a_1_2 = tmp_28 * ( -tmp_45 * tmp_60 + tmp_57 * tmp_66 + tmp_67 ) +
                     tmp_36 * ( -tmp_51 * tmp_60 + tmp_62 * tmp_66 + tmp_68 ) +
                     tmp_44 * ( -tmp_54 * tmp_60 + tmp_64 * tmp_66 + tmp_69 );
      real_t a_2_0 = tmp_28 * ( tmp_19 * tmp_70 - tmp_27 * tmp_57 + tmp_61 ) +
                     tmp_36 * ( -tmp_27 * tmp_62 + tmp_35 * tmp_70 + tmp_63 ) +
                     tmp_44 * ( -tmp_27 * tmp_64 + tmp_43 * tmp_70 + tmp_65 );
      real_t a_2_1 = tmp_28 * ( tmp_45 * tmp_70 - tmp_48 * tmp_57 + tmp_67 ) +
                     tmp_36 * ( -tmp_48 * tmp_62 + tmp_51 * tmp_70 + tmp_68 ) +
                     tmp_44 * ( -tmp_48 * tmp_64 + tmp_54 * tmp_70 + tmp_69 );
      real_t a_2_2 = tmp_28 * ( tmp_21 * ( tmp_57 * tmp_57 ) - tmp_57 * tmp_60 + tmp_57 * tmp_70 ) +
                     tmp_36 * ( tmp_21 * ( tmp_62 * tmp_62 ) - tmp_60 * tmp_62 + tmp_62 * tmp_70 ) +
                     tmp_44 * ( tmp_21 * ( tmp_64 * tmp_64 ) - tmp_60 * tmp_64 + tmp_64 * tmp_70 );

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

      const auto p_affine_8_0 = oppositeVertexInnerElement( 0 );
      const auto p_affine_8_1 = oppositeVertexInnerElement( 1 );

      const auto p_affine_9_0 = oppositeVertexOuterElement( 0 );
      const auto p_affine_9_1 = oppositeVertexOuterElement( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t tmp_0   = -p_affine_0_1;
      real_t tmp_1   = p_affine_2_1 + tmp_0;
      real_t tmp_2   = -p_affine_0_0;
      real_t tmp_3   = p_affine_1_0 + tmp_2;
      real_t tmp_4   = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_0 ) * ( p_affine_2_0 + tmp_2 ) );
      real_t tmp_5   = tmp_1 * tmp_4;
      real_t tmp_6   = tmp_4 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_7   = 0.5 * p_affine_10_0;
      real_t tmp_8   = tmp_3 * tmp_4;
      real_t tmp_9   = tmp_4 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_10  = 0.5 * p_affine_10_1;
      real_t tmp_11  = tmp_10 * ( -tmp_8 - tmp_9 ) + tmp_7 * ( -tmp_5 - tmp_6 );
      real_t tmp_12  = -p_affine_3_1;
      real_t tmp_13  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_14  = p_affine_6_1 + 0.1127016653792583 * tmp_13;
      real_t tmp_15  = tmp_12 + tmp_14;
      real_t tmp_16  = -p_affine_3_0;
      real_t tmp_17  = p_affine_4_0 + tmp_16;
      real_t tmp_18  = p_affine_5_1 + tmp_12;
      real_t tmp_19  = 1.0 / ( tmp_17 * tmp_18 - ( p_affine_4_1 + tmp_12 ) * ( p_affine_5_0 + tmp_16 ) );
      real_t tmp_20  = tmp_17 * tmp_19;
      real_t tmp_21  = tmp_15 * tmp_20;
      real_t tmp_22  = tmp_19 * ( p_affine_3_0 - p_affine_5_0 );
      real_t tmp_23  = tmp_15 * tmp_22;
      real_t tmp_24  = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_25  = p_affine_6_0 + 0.1127016653792583 * tmp_24;
      real_t tmp_26  = tmp_16 + tmp_25;
      real_t tmp_27  = tmp_18 * tmp_19;
      real_t tmp_28  = tmp_26 * tmp_27;
      real_t tmp_29  = tmp_19 * ( p_affine_3_1 - p_affine_4_1 );
      real_t tmp_30  = tmp_26 * tmp_29;
      real_t tmp_31  = -tmp_21 - tmp_23 - tmp_28 - tmp_30 + 1;
      real_t tmp_32  = tmp_10 * ( -tmp_20 - tmp_22 ) + tmp_7 * ( -tmp_27 - tmp_29 );
      real_t tmp_33  = tmp_0 + tmp_14;
      real_t tmp_34  = tmp_33 * tmp_8;
      real_t tmp_35  = tmp_33 * tmp_9;
      real_t tmp_36  = tmp_2 + tmp_25;
      real_t tmp_37  = tmp_36 * tmp_5;
      real_t tmp_38  = tmp_36 * tmp_6;
      real_t tmp_39  = -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1;
      real_t tmp_40  = std::abs( std::pow( ( tmp_13 * tmp_13 ) + ( tmp_24 * tmp_24 ), 1.0 / 2.0 ) );
      real_t tmp_41  = 6 / tmp_40;
      real_t tmp_42  = tmp_39 * tmp_41;
      real_t tmp_43  = 0.27777777777777785 * tmp_40;
      real_t tmp_44  = p_affine_6_1 + 0.5 * tmp_13;
      real_t tmp_45  = tmp_12 + tmp_44;
      real_t tmp_46  = tmp_20 * tmp_45;
      real_t tmp_47  = tmp_22 * tmp_45;
      real_t tmp_48  = p_affine_6_0 + 0.5 * tmp_24;
      real_t tmp_49  = tmp_16 + tmp_48;
      real_t tmp_50  = tmp_27 * tmp_49;
      real_t tmp_51  = tmp_29 * tmp_49;
      real_t tmp_52  = -tmp_46 - tmp_47 - tmp_50 - tmp_51 + 1;
      real_t tmp_53  = tmp_0 + tmp_44;
      real_t tmp_54  = tmp_53 * tmp_8;
      real_t tmp_55  = tmp_53 * tmp_9;
      real_t tmp_56  = tmp_2 + tmp_48;
      real_t tmp_57  = tmp_5 * tmp_56;
      real_t tmp_58  = tmp_56 * tmp_6;
      real_t tmp_59  = -tmp_54 - tmp_55 - tmp_57 - tmp_58 + 1;
      real_t tmp_60  = tmp_41 * tmp_59;
      real_t tmp_61  = 0.44444444444444442 * tmp_40;
      real_t tmp_62  = p_affine_6_1 + 0.8872983346207417 * tmp_13;
      real_t tmp_63  = tmp_12 + tmp_62;
      real_t tmp_64  = tmp_20 * tmp_63;
      real_t tmp_65  = tmp_22 * tmp_63;
      real_t tmp_66  = p_affine_6_0 + 0.8872983346207417 * tmp_24;
      real_t tmp_67  = tmp_16 + tmp_66;
      real_t tmp_68  = tmp_27 * tmp_67;
      real_t tmp_69  = tmp_29 * tmp_67;
      real_t tmp_70  = -tmp_64 - tmp_65 - tmp_68 - tmp_69 + 1;
      real_t tmp_71  = tmp_0 + tmp_62;
      real_t tmp_72  = tmp_71 * tmp_8;
      real_t tmp_73  = tmp_71 * tmp_9;
      real_t tmp_74  = tmp_2 + tmp_66;
      real_t tmp_75  = tmp_5 * tmp_74;
      real_t tmp_76  = tmp_6 * tmp_74;
      real_t tmp_77  = -tmp_72 - tmp_73 - tmp_75 - tmp_76 + 1;
      real_t tmp_78  = tmp_41 * tmp_77;
      real_t tmp_79  = 0.27777777777777785 * tmp_40;
      real_t tmp_80  = tmp_23 + tmp_28;
      real_t tmp_81  = tmp_10 * tmp_22 + tmp_27 * tmp_7;
      real_t tmp_82  = tmp_47 + tmp_50;
      real_t tmp_83  = tmp_65 + tmp_68;
      real_t tmp_84  = tmp_21 + tmp_30;
      real_t tmp_85  = tmp_10 * tmp_20 + tmp_29 * tmp_7;
      real_t tmp_86  = tmp_46 + tmp_51;
      real_t tmp_87  = tmp_64 + tmp_69;
      real_t tmp_88  = tmp_35 + tmp_37;
      real_t tmp_89  = tmp_10 * tmp_9 + tmp_5 * tmp_7;
      real_t tmp_90  = tmp_41 * tmp_88;
      real_t tmp_91  = tmp_55 + tmp_57;
      real_t tmp_92  = tmp_41 * tmp_91;
      real_t tmp_93  = tmp_73 + tmp_75;
      real_t tmp_94  = tmp_41 * tmp_93;
      real_t tmp_95  = tmp_34 + tmp_38;
      real_t tmp_96  = tmp_10 * tmp_8 + tmp_6 * tmp_7;
      real_t tmp_97  = tmp_41 * tmp_95;
      real_t tmp_98  = tmp_54 + tmp_58;
      real_t tmp_99  = tmp_41 * tmp_98;
      real_t tmp_100 = tmp_72 + tmp_76;
      real_t tmp_101 = tmp_100 * tmp_41;
      real_t a_0_0   = tmp_43 * ( tmp_11 * tmp_31 - tmp_31 * tmp_42 - tmp_32 * tmp_39 ) +
                     tmp_61 * ( tmp_11 * tmp_52 - tmp_32 * tmp_59 - tmp_52 * tmp_60 ) +
                     tmp_79 * ( tmp_11 * tmp_70 - tmp_32 * tmp_77 - tmp_70 * tmp_78 );
      real_t a_0_1 = tmp_43 * ( tmp_11 * tmp_80 - tmp_39 * tmp_81 - tmp_42 * tmp_80 ) +
                     tmp_61 * ( tmp_11 * tmp_82 - tmp_59 * tmp_81 - tmp_60 * tmp_82 ) +
                     tmp_79 * ( tmp_11 * tmp_83 - tmp_77 * tmp_81 - tmp_78 * tmp_83 );
      real_t a_0_2 = tmp_43 * ( tmp_11 * tmp_84 - tmp_39 * tmp_85 - tmp_42 * tmp_84 ) +
                     tmp_61 * ( tmp_11 * tmp_86 - tmp_59 * tmp_85 - tmp_60 * tmp_86 ) +
                     tmp_79 * ( tmp_11 * tmp_87 - tmp_77 * tmp_85 - tmp_78 * tmp_87 );
      real_t a_1_0 = tmp_43 * ( tmp_31 * tmp_89 - tmp_31 * tmp_90 - tmp_32 * tmp_88 ) +
                     tmp_61 * ( -tmp_32 * tmp_91 + tmp_52 * tmp_89 - tmp_52 * tmp_92 ) +
                     tmp_79 * ( -tmp_32 * tmp_93 + tmp_70 * tmp_89 - tmp_70 * tmp_94 );
      real_t a_1_1 = tmp_43 * ( tmp_80 * tmp_89 - tmp_80 * tmp_90 - tmp_81 * tmp_88 ) +
                     tmp_61 * ( -tmp_81 * tmp_91 + tmp_82 * tmp_89 - tmp_82 * tmp_92 ) +
                     tmp_79 * ( -tmp_81 * tmp_93 + tmp_83 * tmp_89 - tmp_83 * tmp_94 );
      real_t a_1_2 = tmp_43 * ( tmp_84 * tmp_89 - tmp_84 * tmp_90 - tmp_85 * tmp_88 ) +
                     tmp_61 * ( -tmp_85 * tmp_91 + tmp_86 * tmp_89 - tmp_86 * tmp_92 ) +
                     tmp_79 * ( -tmp_85 * tmp_93 + tmp_87 * tmp_89 - tmp_87 * tmp_94 );
      real_t a_2_0 = tmp_43 * ( tmp_31 * tmp_96 - tmp_31 * tmp_97 - tmp_32 * tmp_95 ) +
                     tmp_61 * ( -tmp_32 * tmp_98 + tmp_52 * tmp_96 - tmp_52 * tmp_99 ) +
                     tmp_79 * ( -tmp_100 * tmp_32 - tmp_101 * tmp_70 + tmp_70 * tmp_96 );
      real_t a_2_1 = tmp_43 * ( tmp_80 * tmp_96 - tmp_80 * tmp_97 - tmp_81 * tmp_95 ) +
                     tmp_61 * ( -tmp_81 * tmp_98 + tmp_82 * tmp_96 - tmp_82 * tmp_99 ) +
                     tmp_79 * ( -tmp_100 * tmp_81 - tmp_101 * tmp_83 + tmp_83 * tmp_96 );
      real_t a_2_2 = tmp_43 * ( tmp_84 * tmp_96 - tmp_84 * tmp_97 - tmp_85 * tmp_95 ) +
                     tmp_61 * ( -tmp_85 * tmp_98 + tmp_86 * tmp_96 - tmp_86 * tmp_99 ) +
                     tmp_79 * ( -tmp_100 * tmp_85 - tmp_101 * tmp_87 + tmp_87 * tmp_96 );

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
      real_t tmp_3 = 0.1127016653792583*tmp_0 + tmp_2;
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
      real_t tmp_14 = 0.1127016653792583*tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6*tmp_7;
      real_t tmp_16 = tmp_14*tmp_15;
      real_t tmp_17 = tmp_7*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_18 = tmp_14*tmp_17;
      real_t tmp_19 = -tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1;
      real_t tmp_20 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_12*tmp_12), 1.0/2.0));
      real_t tmp_21 = 24/tmp_20;
      real_t tmp_22 = p_affine_10_0*(-tmp_15 - tmp_17);
      real_t tmp_23 = p_affine_10_1*(-tmp_10 - tmp_8);
      real_t tmp_24 = tmp_22 + tmp_23;
      real_t tmp_25 = -tmp_22 - tmp_23;
      real_t tmp_26 = 0.27777777777777785*tmp_20;
      real_t tmp_27 = 0.5*tmp_0 + tmp_2;
      real_t tmp_28 = tmp_27*tmp_8;
      real_t tmp_29 = tmp_10*tmp_27;
      real_t tmp_30 = 0.5*tmp_12 + tmp_13;
      real_t tmp_31 = tmp_15*tmp_30;
      real_t tmp_32 = tmp_17*tmp_30;
      real_t tmp_33 = -tmp_28 - tmp_29 - tmp_31 - tmp_32 + 1;
      real_t tmp_34 = 0.44444444444444442*tmp_20;
      real_t tmp_35 = 0.8872983346207417*tmp_0 + tmp_2;
      real_t tmp_36 = tmp_35*tmp_8;
      real_t tmp_37 = tmp_10*tmp_35;
      real_t tmp_38 = 0.8872983346207417*tmp_12 + tmp_13;
      real_t tmp_39 = tmp_15*tmp_38;
      real_t tmp_40 = tmp_17*tmp_38;
      real_t tmp_41 = -tmp_36 - tmp_37 - tmp_39 - tmp_40 + 1;
      real_t tmp_42 = 0.27777777777777785*tmp_20;
      real_t tmp_43 = tmp_11 + tmp_16;
      real_t tmp_44 = p_affine_10_0*tmp_15;
      real_t tmp_45 = p_affine_10_1*tmp_10;
      real_t tmp_46 = tmp_44 + tmp_45;
      real_t tmp_47 = tmp_19*tmp_21;
      real_t tmp_48 = tmp_43*tmp_47;
      real_t tmp_49 = tmp_29 + tmp_31;
      real_t tmp_50 = tmp_21*tmp_33;
      real_t tmp_51 = tmp_49*tmp_50;
      real_t tmp_52 = tmp_37 + tmp_39;
      real_t tmp_53 = tmp_21*tmp_41;
      real_t tmp_54 = tmp_52*tmp_53;
      real_t tmp_55 = tmp_18 + tmp_9;
      real_t tmp_56 = p_affine_10_0*tmp_17;
      real_t tmp_57 = p_affine_10_1*tmp_8;
      real_t tmp_58 = tmp_56 + tmp_57;
      real_t tmp_59 = tmp_47*tmp_55;
      real_t tmp_60 = tmp_28 + tmp_32;
      real_t tmp_61 = tmp_50*tmp_60;
      real_t tmp_62 = tmp_36 + tmp_40;
      real_t tmp_63 = tmp_53*tmp_62;
      real_t tmp_64 = -tmp_44 - tmp_45;
      real_t tmp_65 = tmp_21*tmp_43*tmp_55;
      real_t tmp_66 = tmp_21*tmp_49*tmp_60;
      real_t tmp_67 = tmp_21*tmp_52*tmp_62;
      real_t tmp_68 = -tmp_56 - tmp_57;
      real_t a_0_0 = tmp_26*((tmp_19*tmp_19)*tmp_21 - tmp_19*tmp_24 + tmp_19*tmp_25) + tmp_34*(tmp_21*(tmp_33*tmp_33) - tmp_24*tmp_33 + tmp_25*tmp_33) + tmp_42*(tmp_21*(tmp_41*tmp_41) - tmp_24*tmp_41 + tmp_25*tmp_41);
      real_t a_0_1 = tmp_26*(-tmp_19*tmp_46 + tmp_25*tmp_43 + tmp_48) + tmp_34*(tmp_25*tmp_49 - tmp_33*tmp_46 + tmp_51) + tmp_42*(tmp_25*tmp_52 - tmp_41*tmp_46 + tmp_54);
      real_t a_0_2 = tmp_26*(-tmp_19*tmp_58 + tmp_25*tmp_55 + tmp_59) + tmp_34*(tmp_25*tmp_60 - tmp_33*tmp_58 + tmp_61) + tmp_42*(tmp_25*tmp_62 - tmp_41*tmp_58 + tmp_63);
      real_t a_1_0 = tmp_26*(tmp_19*tmp_64 - tmp_24*tmp_43 + tmp_48) + tmp_34*(-tmp_24*tmp_49 + tmp_33*tmp_64 + tmp_51) + tmp_42*(-tmp_24*tmp_52 + tmp_41*tmp_64 + tmp_54);
      real_t a_1_1 = tmp_26*(tmp_21*(tmp_43*tmp_43) - tmp_43*tmp_46 + tmp_43*tmp_64) + tmp_34*(tmp_21*(tmp_49*tmp_49) - tmp_46*tmp_49 + tmp_49*tmp_64) + tmp_42*(tmp_21*(tmp_52*tmp_52) - tmp_46*tmp_52 + tmp_52*tmp_64);
      real_t a_1_2 = tmp_26*(-tmp_43*tmp_58 + tmp_55*tmp_64 + tmp_65) + tmp_34*(-tmp_49*tmp_58 + tmp_60*tmp_64 + tmp_66) + tmp_42*(-tmp_52*tmp_58 + tmp_62*tmp_64 + tmp_67);
      real_t a_2_0 = tmp_26*(tmp_19*tmp_68 - tmp_24*tmp_55 + tmp_59) + tmp_34*(-tmp_24*tmp_60 + tmp_33*tmp_68 + tmp_61) + tmp_42*(-tmp_24*tmp_62 + tmp_41*tmp_68 + tmp_63);
      real_t a_2_1 = tmp_26*(tmp_43*tmp_68 - tmp_46*tmp_55 + tmp_65) + tmp_34*(-tmp_46*tmp_60 + tmp_49*tmp_68 + tmp_66) + tmp_42*(-tmp_46*tmp_62 + tmp_52*tmp_68 + tmp_67);
      real_t a_2_2 = tmp_26*(tmp_21*(tmp_55*tmp_55) - tmp_55*tmp_58 + tmp_55*tmp_68) + tmp_34*(tmp_21*(tmp_60*tmp_60) - tmp_58*tmp_60 + tmp_60*tmp_68) + tmp_42*(tmp_21*(tmp_62*tmp_62) - tmp_58*tmp_62 + tmp_62*tmp_68);


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

   virtual void integrateRHSDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
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

      const auto p_affine_8_0 = oppositeVertex( 0 );
      const auto p_affine_8_1 = oppositeVertex( 1 );

      const auto p_affine_10_0 = outwardNormal( 0 );
      const auto p_affine_10_1 = outwardNormal( 1 );

      real_t Scalar_Variable_Coefficient_2D_g_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_g_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_g_out0_id2 = 0;
      Scalar_Variable_Coefficient_2D_g( 0.8872983346207417*p_affine_6_0 + 0.1127016653792583*p_affine_7_0, 0.8872983346207417*p_affine_6_1 + 0.1127016653792583*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g_out0_id0 );
      Scalar_Variable_Coefficient_2D_g( 0.5*p_affine_6_0 + 0.5*p_affine_7_0, 0.5*p_affine_6_1 + 0.5*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g_out0_id1 );
      Scalar_Variable_Coefficient_2D_g( 0.1127016653792583*p_affine_6_0 + 0.8872983346207417*p_affine_7_0, 0.1127016653792583*p_affine_6_1 + 0.8872983346207417*p_affine_7_1, &Scalar_Variable_Coefficient_2D_g_out0_id2 );
      real_t tmp_0 = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_1 = -p_affine_0_1;
      real_t tmp_2 = p_affine_6_1 + tmp_1;
      real_t tmp_3 = 0.1127016653792583*tmp_0 + tmp_2;
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
      real_t tmp_14 = 0.1127016653792583*tmp_12 + tmp_13;
      real_t tmp_15 = tmp_6*tmp_7;
      real_t tmp_16 = tmp_14*tmp_15;
      real_t tmp_17 = tmp_7*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_18 = tmp_14*tmp_17;
      real_t tmp_19 = std::abs(std::pow((tmp_0*tmp_0) + (tmp_12*tmp_12), 1.0/2.0));
      real_t tmp_20 = 24/tmp_19;
      real_t tmp_21 = -p_affine_10_0*(-tmp_15 - tmp_17) - p_affine_10_1*(-tmp_10 - tmp_8);
      real_t tmp_22 = 0.27777777777777785*Scalar_Variable_Coefficient_2D_g_out0_id0*tmp_19;
      real_t tmp_23 = 0.5*tmp_0 + tmp_2;
      real_t tmp_24 = tmp_23*tmp_8;
      real_t tmp_25 = tmp_10*tmp_23;
      real_t tmp_26 = 0.5*tmp_12 + tmp_13;
      real_t tmp_27 = tmp_15*tmp_26;
      real_t tmp_28 = tmp_17*tmp_26;
      real_t tmp_29 = 0.44444444444444442*Scalar_Variable_Coefficient_2D_g_out0_id1*tmp_19;
      real_t tmp_30 = 0.8872983346207417*tmp_0 + tmp_2;
      real_t tmp_31 = tmp_30*tmp_8;
      real_t tmp_32 = tmp_10*tmp_30;
      real_t tmp_33 = 0.8872983346207417*tmp_12 + tmp_13;
      real_t tmp_34 = tmp_15*tmp_33;
      real_t tmp_35 = tmp_17*tmp_33;
      real_t tmp_36 = 0.27777777777777785*Scalar_Variable_Coefficient_2D_g_out0_id2*tmp_19;
      real_t tmp_37 = -p_affine_10_0*tmp_15 - p_affine_10_1*tmp_10;
      real_t tmp_38 = -p_affine_10_0*tmp_17 - p_affine_10_1*tmp_8;
      real_t a_0_0 = tmp_22*(tmp_20*(-tmp_11 - tmp_16 - tmp_18 - tmp_9 + 1) + tmp_21) + tmp_29*(tmp_20*(-tmp_24 - tmp_25 - tmp_27 - tmp_28 + 1) + tmp_21) + tmp_36*(tmp_20*(-tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1) + tmp_21);
      real_t a_1_0 = tmp_22*(tmp_20*(tmp_11 + tmp_16) + tmp_37) + tmp_29*(tmp_20*(tmp_25 + tmp_27) + tmp_37) + tmp_36*(tmp_20*(tmp_32 + tmp_34) + tmp_37);
      real_t a_2_0 = tmp_22*(tmp_20*(tmp_18 + tmp_9) + tmp_38) + tmp_29*(tmp_20*(tmp_24 + tmp_28) + tmp_38) + tmp_36*(tmp_20*(tmp_31 + tmp_35) + tmp_38);

      elMat( 0, 0 ) = a_0_0;
      elMat( 1, 0 ) = a_1_0;
      elMat( 2, 0 ) = a_2_0;
   }

 private:
   void Scalar_Variable_Coefficient_2D_g( real_t in_0, real_t in_1, real_t* out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_g( Point3D( { in_0, in_1, 0 } ) );
   }

   std::function< real_t( const Point3D& ) > callback_Scalar_Variable_Coefficient_2D_g;
};

} // namespace dg
} // namespace hyteg
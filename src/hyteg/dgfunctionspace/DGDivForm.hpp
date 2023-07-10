
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
#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"


#include "Eigen/Eigen"

namespace hyteg {
namespace dg {

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

      real_t tmp_0  = 0.063089014491502282;
      real_t tmp_1  = -p_affine_0_1;
      real_t tmp_2  = p_affine_2_1 + tmp_1;
      real_t tmp_3  = -p_affine_0_0;
      real_t tmp_4  = 1.0 / ( tmp_2 * ( p_affine_1_0 + tmp_3 ) - ( p_affine_1_1 + tmp_1 ) * ( p_affine_2_0 + tmp_3 ) );
      real_t tmp_5  = tmp_2 * tmp_4;
      real_t tmp_6  = tmp_4 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_7  = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                               p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_8  = tmp_7 * ( tmp_5 + tmp_6 );
      real_t tmp_9  = 0.025422453185103409 * tmp_8;
      real_t tmp_10 = 0.24928674517091043;
      real_t tmp_11 = 0.058393137863189684 * tmp_8;
      real_t tmp_12 = 0.63650249912139867;
      real_t tmp_13 = 0.041425537809186785 * tmp_8;
      real_t tmp_14 = 0.053145049844816938;
      real_t tmp_15 = 0.041425537809186785 * tmp_8;
      real_t tmp_16 = 0.063089014491502227;
      real_t tmp_17 = 0.025422453185103409 * tmp_8;
      real_t tmp_18 = 0.24928674517091043;
      real_t tmp_19 = 0.058393137863189684 * tmp_8;
      real_t tmp_20 = 0.87382197101699566;
      real_t tmp_21 = 0.025422453185103409 * tmp_8;
      real_t tmp_22 = 0.50142650965817914;
      real_t tmp_23 = 0.058393137863189684 * tmp_8;
      real_t tmp_24 = 0.053145049844816938;
      real_t tmp_25 = 0.041425537809186785 * tmp_8;
      real_t tmp_26 = 0.63650249912139867;
      real_t tmp_27 = 0.041425537809186785 * tmp_8;
      real_t tmp_28 = 0.31035245103378439;
      real_t tmp_29 = 0.041425537809186785 * tmp_8;
      real_t tmp_30 = 0.31035245103378439;
      real_t tmp_31 = 0.041425537809186785 * tmp_8;
      real_t tmp_32 = tmp_5 * tmp_7;
      real_t tmp_33 = 0.025422453185103409 * tmp_0;
      real_t tmp_34 = 0.058393137863189684 * tmp_10;
      real_t tmp_35 = 0.041425537809186785 * tmp_12;
      real_t tmp_36 = 0.041425537809186785 * tmp_14;
      real_t tmp_37 = 0.025422453185103409 * tmp_16;
      real_t tmp_38 = 0.058393137863189684 * tmp_18;
      real_t tmp_39 = 0.025422453185103409 * tmp_20;
      real_t tmp_40 = 0.058393137863189684 * tmp_22;
      real_t tmp_41 = 0.041425537809186785 * tmp_24;
      real_t tmp_42 = 0.041425537809186785 * tmp_26;
      real_t tmp_43 = 0.041425537809186785 * tmp_28;
      real_t tmp_44 = 0.041425537809186785 * tmp_30;
      real_t tmp_45 = tmp_6 * tmp_7;
      real_t tmp_46 = 0.001603877517404526;
      real_t tmp_47 = 0.012856517194473826;
      real_t tmp_48 = 0.026367458342995378;
      real_t tmp_49 = 0.014556635278230808;
      real_t tmp_50 = 0.022214698150294358;
      real_t tmp_51 = 0.029279867306728068;
      real_t tmp_52 = 0.001603877517404526;
      real_t tmp_53 = 0.014556635278230808;
      real_t tmp_54 = 0.012856517194473826;
      real_t tmp_55 = 0.0022015622717175805;
      real_t tmp_56 = 0.026367458342995378;
      real_t tmp_57 = 0.0022015622717175805;
      real_t tmp_58 = 0.022214698150294358;
      real_t tmp_59 = 0.0022015622717175805;
      real_t tmp_60 = 0.012856517194473826;
      real_t tmp_61 = 0.029279867306728068;
      real_t tmp_62 = 0.001603877517404526;
      real_t tmp_63 = 0.014556635278230808;
      real_t tmp_64 = 0.001603877517404526;
      real_t tmp_65 = 0.014556635278230808;
      real_t tmp_66 = 0.026367458342995378;
      real_t tmp_67 = 0.012856517194473826;
      real_t tmp_68 = 0.0022015622717175805;
      real_t tmp_69 = 0.026367458342995378;
      real_t a_0_0  = tmp_0 * tmp_9 + tmp_10 * tmp_11 + tmp_12 * tmp_13 + tmp_14 * tmp_15 + tmp_16 * tmp_17 + tmp_18 * tmp_19 +
                     tmp_20 * tmp_21 + tmp_22 * tmp_23 + tmp_24 * tmp_25 + tmp_26 * tmp_27 + tmp_28 * tmp_29 + tmp_30 * tmp_31;
      real_t a_0_1 = -tmp_32 * tmp_33 - tmp_32 * tmp_34 - tmp_32 * tmp_35 - tmp_32 * tmp_36 - tmp_32 * tmp_37 - tmp_32 * tmp_38 -
                     tmp_32 * tmp_39 - tmp_32 * tmp_40 - tmp_32 * tmp_41 - tmp_32 * tmp_42 - tmp_32 * tmp_43 - tmp_32 * tmp_44;
      real_t a_0_2 = -tmp_33 * tmp_45 - tmp_34 * tmp_45 - tmp_35 * tmp_45 - tmp_36 * tmp_45 - tmp_37 * tmp_45 - tmp_38 * tmp_45 -
                     tmp_39 * tmp_45 - tmp_40 * tmp_45 - tmp_41 * tmp_45 - tmp_42 * tmp_45 - tmp_43 * tmp_45 - tmp_44 * tmp_45;
      real_t a_1_0 = 0.24928674517091043 * tmp_11 + 0.31035245103378439 * tmp_13 + 0.63650249912139867 * tmp_15 +
                     0.87382197101699555 * tmp_17 + 0.50142650965817914 * tmp_19 + 0.063089014491502227 * tmp_21 +
                     0.24928674517091043 * tmp_23 + 0.31035245103378439 * tmp_25 + 0.053145049844816938 * tmp_27 +
                     0.63650249912139867 * tmp_29 + 0.053145049844816938 * tmp_31 + 0.063089014491502227 * tmp_9;
      real_t a_1_1 = -tmp_32 * tmp_46 - tmp_32 * tmp_47 - tmp_32 * tmp_48 - tmp_32 * tmp_49 - tmp_32 * tmp_50 - tmp_32 * tmp_51 -
                     tmp_32 * tmp_52 - tmp_32 * tmp_53 - tmp_32 * tmp_54 - tmp_32 * tmp_55 - tmp_32 * tmp_56 - tmp_32 * tmp_57;
      real_t a_1_2 = -tmp_45 * tmp_46 - tmp_45 * tmp_47 - tmp_45 * tmp_48 - tmp_45 * tmp_49 - tmp_45 * tmp_50 - tmp_45 * tmp_51 -
                     tmp_45 * tmp_52 - tmp_45 * tmp_53 - tmp_45 * tmp_54 - tmp_45 * tmp_55 - tmp_45 * tmp_56 - tmp_45 * tmp_57;
      real_t a_2_0 = 0.50142650965817914 * tmp_11 + 0.053145049844816938 * tmp_13 + 0.31035245103378439 * tmp_15 +
                     0.063089014491502227 * tmp_17 + 0.24928674517091043 * tmp_19 + 0.063089014491502227 * tmp_21 +
                     0.24928674517091043 * tmp_23 + 0.63650249912139867 * tmp_25 + 0.31035245103378439 * tmp_27 +
                     0.053145049844816938 * tmp_29 + 0.63650249912139867 * tmp_31 + 0.87382197101699555 * tmp_9;
      real_t a_2_1 = -tmp_32 * tmp_58 - tmp_32 * tmp_59 - tmp_32 * tmp_60 - tmp_32 * tmp_61 - tmp_32 * tmp_62 - tmp_32 * tmp_63 -
                     tmp_32 * tmp_64 - tmp_32 * tmp_65 - tmp_32 * tmp_66 - tmp_32 * tmp_67 - tmp_32 * tmp_68 - tmp_32 * tmp_69;
      real_t a_2_2 = -tmp_45 * tmp_58 - tmp_45 * tmp_59 - tmp_45 * tmp_60 - tmp_45 * tmp_61 - tmp_45 * tmp_62 - tmp_45 * tmp_63 -
                     tmp_45 * tmp_64 - tmp_45 * tmp_65 - tmp_45 * tmp_66 - tmp_45 * tmp_67 - tmp_45 * tmp_68 - tmp_45 * tmp_69;
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + tmp_2;
      real_t tmp_7  = tmp_4 * ( 0.046910077030668018 * tmp_5 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = 0.5 * p_affine_10_0 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_19 = 0.11846344252809471 * tmp_18;
      real_t tmp_20 = tmp_4 * ( 0.23076534494715845 * tmp_5 + tmp_6 );
      real_t tmp_21 = tmp_1 * tmp_20;
      real_t tmp_22 = tmp_20 * tmp_9;
      real_t tmp_23 = tmp_4 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_24 = tmp_23 * tmp_3;
      real_t tmp_25 = tmp_15 * tmp_23;
      real_t tmp_26 = -tmp_21 - tmp_22 - tmp_24 - tmp_25 + 1;
      real_t tmp_27 = 0.2393143352496831 * tmp_18;
      real_t tmp_28 = tmp_4 * ( 0.5 * tmp_5 + tmp_6 );
      real_t tmp_29 = tmp_1 * tmp_28;
      real_t tmp_30 = tmp_28 * tmp_9;
      real_t tmp_31 = tmp_4 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_32 = tmp_3 * tmp_31;
      real_t tmp_33 = tmp_15 * tmp_31;
      real_t tmp_34 = -tmp_29 - tmp_30 - tmp_32 - tmp_33 + 1;
      real_t tmp_35 = 0.2844444444444445 * tmp_18;
      real_t tmp_36 = tmp_4 * ( 0.7692346550528415 * tmp_5 + tmp_6 );
      real_t tmp_37 = tmp_1 * tmp_36;
      real_t tmp_38 = tmp_36 * tmp_9;
      real_t tmp_39 = tmp_4 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_40 = tmp_3 * tmp_39;
      real_t tmp_41 = tmp_15 * tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = 0.2393143352496831 * tmp_18;
      real_t tmp_44 = tmp_4 * ( 0.95308992296933193 * tmp_5 + tmp_6 );
      real_t tmp_45 = tmp_1 * tmp_44;
      real_t tmp_46 = tmp_44 * tmp_9;
      real_t tmp_47 = tmp_4 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_48 = tmp_3 * tmp_47;
      real_t tmp_49 = tmp_15 * tmp_47;
      real_t tmp_50 = -tmp_45 - tmp_46 - tmp_48 - tmp_49 + 1;
      real_t tmp_51 = 0.11846344252809471 * tmp_18;
      real_t tmp_52 = tmp_10 + tmp_14;
      real_t tmp_53 = tmp_17 * tmp_19;
      real_t tmp_54 = tmp_22 + tmp_24;
      real_t tmp_55 = tmp_26 * tmp_27;
      real_t tmp_56 = tmp_30 + tmp_32;
      real_t tmp_57 = tmp_34 * tmp_35;
      real_t tmp_58 = tmp_38 + tmp_40;
      real_t tmp_59 = tmp_42 * tmp_43;
      real_t tmp_60 = tmp_46 + tmp_48;
      real_t tmp_61 = tmp_50 * tmp_51;
      real_t tmp_62 = tmp_52 * tmp_53 + tmp_54 * tmp_55 + tmp_56 * tmp_57 + tmp_58 * tmp_59 + tmp_60 * tmp_61;
      real_t tmp_63 = tmp_16 + tmp_8;
      real_t tmp_64 = tmp_21 + tmp_25;
      real_t tmp_65 = tmp_29 + tmp_33;
      real_t tmp_66 = tmp_37 + tmp_41;
      real_t tmp_67 = tmp_45 + tmp_49;
      real_t tmp_68 = tmp_53 * tmp_63 + tmp_55 * tmp_64 + tmp_57 * tmp_65 + tmp_59 * tmp_66 + tmp_61 * tmp_67;
      real_t tmp_69 = tmp_19 * tmp_52 * tmp_63 + tmp_27 * tmp_54 * tmp_64 + tmp_35 * tmp_56 * tmp_65 + tmp_43 * tmp_58 * tmp_66 +
                      tmp_51 * tmp_60 * tmp_67;
      real_t a_0_0 = ( tmp_17 * tmp_17 ) * tmp_19 + ( tmp_26 * tmp_26 ) * tmp_27 + ( tmp_34 * tmp_34 ) * tmp_35 +
                     ( tmp_42 * tmp_42 ) * tmp_43 + ( tmp_50 * tmp_50 ) * tmp_51;
      real_t a_0_1 = tmp_62;
      real_t a_0_2 = tmp_68;
      real_t a_1_0 = tmp_62;
      real_t a_1_1 = tmp_19 * ( tmp_52 * tmp_52 ) + tmp_27 * ( tmp_54 * tmp_54 ) + tmp_35 * ( tmp_56 * tmp_56 ) +
                     tmp_43 * ( tmp_58 * tmp_58 ) + tmp_51 * ( tmp_60 * tmp_60 );
      real_t a_1_2 = tmp_69;
      real_t a_2_0 = tmp_68;
      real_t a_2_1 = tmp_69;
      real_t a_2_2 = tmp_19 * ( tmp_63 * tmp_63 ) + tmp_27 * ( tmp_64 * tmp_64 ) + tmp_35 * ( tmp_65 * tmp_65 ) +
                     tmp_43 * ( tmp_66 * tmp_66 ) + tmp_51 * ( tmp_67 * tmp_67 );
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

      real_t tmp_0   = -p_affine_3_0;
      real_t tmp_1   = p_affine_4_0 + tmp_0;
      real_t tmp_2   = -p_affine_3_1;
      real_t tmp_3   = p_affine_5_1 + tmp_2;
      real_t tmp_4   = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_4_1 + tmp_2 ) * ( p_affine_5_0 + tmp_0 ) );
      real_t tmp_5   = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6   = p_affine_6_1 + 0.046910077030668018 * tmp_5;
      real_t tmp_7   = tmp_4 * ( tmp_2 + tmp_6 );
      real_t tmp_8   = tmp_1 * tmp_7;
      real_t tmp_9   = p_affine_3_0 - p_affine_5_0;
      real_t tmp_10  = tmp_7 * tmp_9;
      real_t tmp_11  = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12  = p_affine_6_0 + 0.046910077030668018 * tmp_11;
      real_t tmp_13  = tmp_4 * ( tmp_0 + tmp_12 );
      real_t tmp_14  = tmp_13 * tmp_3;
      real_t tmp_15  = p_affine_3_1 - p_affine_4_1;
      real_t tmp_16  = tmp_13 * tmp_15;
      real_t tmp_17  = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18  = -p_affine_0_0;
      real_t tmp_19  = p_affine_1_0 + tmp_18;
      real_t tmp_20  = -p_affine_0_1;
      real_t tmp_21  = p_affine_2_1 + tmp_20;
      real_t tmp_22  = 1.0 / ( tmp_19 * tmp_21 - ( p_affine_1_1 + tmp_20 ) * ( p_affine_2_0 + tmp_18 ) );
      real_t tmp_23  = tmp_22 * ( tmp_20 + tmp_6 );
      real_t tmp_24  = tmp_19 * tmp_23;
      real_t tmp_25  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_26  = tmp_23 * tmp_25;
      real_t tmp_27  = tmp_22 * ( tmp_12 + tmp_18 );
      real_t tmp_28  = tmp_21 * tmp_27;
      real_t tmp_29  = p_affine_0_1 - p_affine_1_1;
      real_t tmp_30  = tmp_27 * tmp_29;
      real_t tmp_31  = 0.5 * p_affine_10_0 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_32  = 0.11846344252809471 * tmp_31;
      real_t tmp_33  = tmp_32 * ( -tmp_24 - tmp_26 - tmp_28 - tmp_30 + 1 );
      real_t tmp_34  = p_affine_6_1 + 0.23076534494715845 * tmp_5;
      real_t tmp_35  = tmp_4 * ( tmp_2 + tmp_34 );
      real_t tmp_36  = tmp_1 * tmp_35;
      real_t tmp_37  = tmp_35 * tmp_9;
      real_t tmp_38  = p_affine_6_0 + 0.23076534494715845 * tmp_11;
      real_t tmp_39  = tmp_4 * ( tmp_0 + tmp_38 );
      real_t tmp_40  = tmp_3 * tmp_39;
      real_t tmp_41  = tmp_15 * tmp_39;
      real_t tmp_42  = -tmp_36 - tmp_37 - tmp_40 - tmp_41 + 1;
      real_t tmp_43  = tmp_22 * ( tmp_20 + tmp_34 );
      real_t tmp_44  = tmp_19 * tmp_43;
      real_t tmp_45  = tmp_25 * tmp_43;
      real_t tmp_46  = tmp_22 * ( tmp_18 + tmp_38 );
      real_t tmp_47  = tmp_21 * tmp_46;
      real_t tmp_48  = tmp_29 * tmp_46;
      real_t tmp_49  = 0.2393143352496831 * tmp_31;
      real_t tmp_50  = tmp_49 * ( -tmp_44 - tmp_45 - tmp_47 - tmp_48 + 1 );
      real_t tmp_51  = p_affine_6_1 + 0.5 * tmp_5;
      real_t tmp_52  = tmp_4 * ( tmp_2 + tmp_51 );
      real_t tmp_53  = tmp_1 * tmp_52;
      real_t tmp_54  = tmp_52 * tmp_9;
      real_t tmp_55  = p_affine_6_0 + 0.5 * tmp_11;
      real_t tmp_56  = tmp_4 * ( tmp_0 + tmp_55 );
      real_t tmp_57  = tmp_3 * tmp_56;
      real_t tmp_58  = tmp_15 * tmp_56;
      real_t tmp_59  = -tmp_53 - tmp_54 - tmp_57 - tmp_58 + 1;
      real_t tmp_60  = tmp_22 * ( tmp_20 + tmp_51 );
      real_t tmp_61  = tmp_19 * tmp_60;
      real_t tmp_62  = tmp_25 * tmp_60;
      real_t tmp_63  = tmp_22 * ( tmp_18 + tmp_55 );
      real_t tmp_64  = tmp_21 * tmp_63;
      real_t tmp_65  = tmp_29 * tmp_63;
      real_t tmp_66  = 0.2844444444444445 * tmp_31;
      real_t tmp_67  = tmp_66 * ( -tmp_61 - tmp_62 - tmp_64 - tmp_65 + 1 );
      real_t tmp_68  = p_affine_6_1 + 0.7692346550528415 * tmp_5;
      real_t tmp_69  = tmp_4 * ( tmp_2 + tmp_68 );
      real_t tmp_70  = tmp_1 * tmp_69;
      real_t tmp_71  = tmp_69 * tmp_9;
      real_t tmp_72  = p_affine_6_0 + 0.7692346550528415 * tmp_11;
      real_t tmp_73  = tmp_4 * ( tmp_0 + tmp_72 );
      real_t tmp_74  = tmp_3 * tmp_73;
      real_t tmp_75  = tmp_15 * tmp_73;
      real_t tmp_76  = -tmp_70 - tmp_71 - tmp_74 - tmp_75 + 1;
      real_t tmp_77  = tmp_22 * ( tmp_20 + tmp_68 );
      real_t tmp_78  = tmp_19 * tmp_77;
      real_t tmp_79  = tmp_25 * tmp_77;
      real_t tmp_80  = tmp_22 * ( tmp_18 + tmp_72 );
      real_t tmp_81  = tmp_21 * tmp_80;
      real_t tmp_82  = tmp_29 * tmp_80;
      real_t tmp_83  = 0.2393143352496831 * tmp_31;
      real_t tmp_84  = tmp_83 * ( -tmp_78 - tmp_79 - tmp_81 - tmp_82 + 1 );
      real_t tmp_85  = p_affine_6_1 + 0.95308992296933193 * tmp_5;
      real_t tmp_86  = tmp_4 * ( tmp_2 + tmp_85 );
      real_t tmp_87  = tmp_1 * tmp_86;
      real_t tmp_88  = tmp_86 * tmp_9;
      real_t tmp_89  = p_affine_6_0 + 0.95308992296933193 * tmp_11;
      real_t tmp_90  = tmp_4 * ( tmp_0 + tmp_89 );
      real_t tmp_91  = tmp_3 * tmp_90;
      real_t tmp_92  = tmp_15 * tmp_90;
      real_t tmp_93  = -tmp_87 - tmp_88 - tmp_91 - tmp_92 + 1;
      real_t tmp_94  = tmp_22 * ( tmp_20 + tmp_85 );
      real_t tmp_95  = tmp_19 * tmp_94;
      real_t tmp_96  = tmp_25 * tmp_94;
      real_t tmp_97  = tmp_22 * ( tmp_18 + tmp_89 );
      real_t tmp_98  = tmp_21 * tmp_97;
      real_t tmp_99  = tmp_29 * tmp_97;
      real_t tmp_100 = 0.11846344252809471 * tmp_31;
      real_t tmp_101 = tmp_100 * ( -tmp_95 - tmp_96 - tmp_98 - tmp_99 + 1 );
      real_t tmp_102 = tmp_10 + tmp_14;
      real_t tmp_103 = tmp_37 + tmp_40;
      real_t tmp_104 = tmp_54 + tmp_57;
      real_t tmp_105 = tmp_71 + tmp_74;
      real_t tmp_106 = tmp_88 + tmp_91;
      real_t tmp_107 = tmp_16 + tmp_8;
      real_t tmp_108 = tmp_36 + tmp_41;
      real_t tmp_109 = tmp_53 + tmp_58;
      real_t tmp_110 = tmp_70 + tmp_75;
      real_t tmp_111 = tmp_87 + tmp_92;
      real_t tmp_112 = tmp_32 * ( tmp_26 + tmp_28 );
      real_t tmp_113 = tmp_49 * ( tmp_45 + tmp_47 );
      real_t tmp_114 = tmp_66 * ( tmp_62 + tmp_64 );
      real_t tmp_115 = tmp_83 * ( tmp_79 + tmp_81 );
      real_t tmp_116 = tmp_100 * ( tmp_96 + tmp_98 );
      real_t tmp_117 = tmp_32 * ( tmp_24 + tmp_30 );
      real_t tmp_118 = tmp_49 * ( tmp_44 + tmp_48 );
      real_t tmp_119 = tmp_66 * ( tmp_61 + tmp_65 );
      real_t tmp_120 = tmp_83 * ( tmp_78 + tmp_82 );
      real_t tmp_121 = tmp_100 * ( tmp_95 + tmp_99 );
      real_t a_0_0   = -tmp_101 * tmp_93 - tmp_17 * tmp_33 - tmp_42 * tmp_50 - tmp_59 * tmp_67 - tmp_76 * tmp_84;
      real_t a_0_1   = -tmp_101 * tmp_106 - tmp_102 * tmp_33 - tmp_103 * tmp_50 - tmp_104 * tmp_67 - tmp_105 * tmp_84;
      real_t a_0_2   = -tmp_101 * tmp_111 - tmp_107 * tmp_33 - tmp_108 * tmp_50 - tmp_109 * tmp_67 - tmp_110 * tmp_84;
      real_t a_1_0   = -tmp_112 * tmp_17 - tmp_113 * tmp_42 - tmp_114 * tmp_59 - tmp_115 * tmp_76 - tmp_116 * tmp_93;
      real_t a_1_1   = -tmp_102 * tmp_112 - tmp_103 * tmp_113 - tmp_104 * tmp_114 - tmp_105 * tmp_115 - tmp_106 * tmp_116;
      real_t a_1_2   = -tmp_107 * tmp_112 - tmp_108 * tmp_113 - tmp_109 * tmp_114 - tmp_110 * tmp_115 - tmp_111 * tmp_116;
      real_t a_2_0   = -tmp_117 * tmp_17 - tmp_118 * tmp_42 - tmp_119 * tmp_59 - tmp_120 * tmp_76 - tmp_121 * tmp_93;
      real_t a_2_1   = -tmp_102 * tmp_117 - tmp_103 * tmp_118 - tmp_104 * tmp_119 - tmp_105 * tmp_120 - tmp_106 * tmp_121;
      real_t a_2_2   = -tmp_107 * tmp_117 - tmp_108 * tmp_118 - tmp_109 * tmp_119 - tmp_110 * tmp_120 - tmp_111 * tmp_121;
      elMat( 0, 0 )  = a_0_0;
      elMat( 0, 1 )  = a_0_1;
      elMat( 0, 2 )  = a_0_2;
      elMat( 1, 0 )  = a_1_0;
      elMat( 1, 1 )  = a_1_1;
      elMat( 1, 2 )  = a_1_2;
      elMat( 2, 0 )  = a_2_0;
      elMat( 2, 1 )  = a_2_1;
      elMat( 2, 2 )  = a_2_2;
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

      real_t tmp_0  = 0.063089014491502282;
      real_t tmp_1  = -p_affine_0_0;
      real_t tmp_2  = p_affine_1_0 + tmp_1;
      real_t tmp_3  = -p_affine_0_1;
      real_t tmp_4  = 1.0 / ( tmp_2 * ( p_affine_2_1 + tmp_3 ) - ( p_affine_1_1 + tmp_3 ) * ( p_affine_2_0 + tmp_1 ) );
      real_t tmp_5  = tmp_2 * tmp_4;
      real_t tmp_6  = tmp_4 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_7  = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                               p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_8  = tmp_7 * ( tmp_5 + tmp_6 );
      real_t tmp_9  = 0.025422453185103409 * tmp_8;
      real_t tmp_10 = 0.24928674517091043;
      real_t tmp_11 = 0.058393137863189684 * tmp_8;
      real_t tmp_12 = 0.63650249912139867;
      real_t tmp_13 = 0.041425537809186785 * tmp_8;
      real_t tmp_14 = 0.053145049844816938;
      real_t tmp_15 = 0.041425537809186785 * tmp_8;
      real_t tmp_16 = 0.063089014491502227;
      real_t tmp_17 = 0.025422453185103409 * tmp_8;
      real_t tmp_18 = 0.24928674517091043;
      real_t tmp_19 = 0.058393137863189684 * tmp_8;
      real_t tmp_20 = 0.87382197101699566;
      real_t tmp_21 = 0.025422453185103409 * tmp_8;
      real_t tmp_22 = 0.50142650965817914;
      real_t tmp_23 = 0.058393137863189684 * tmp_8;
      real_t tmp_24 = 0.053145049844816938;
      real_t tmp_25 = 0.041425537809186785 * tmp_8;
      real_t tmp_26 = 0.63650249912139867;
      real_t tmp_27 = 0.041425537809186785 * tmp_8;
      real_t tmp_28 = 0.31035245103378439;
      real_t tmp_29 = 0.041425537809186785 * tmp_8;
      real_t tmp_30 = 0.31035245103378439;
      real_t tmp_31 = 0.041425537809186785 * tmp_8;
      real_t tmp_32 = tmp_6 * tmp_7;
      real_t tmp_33 = 0.025422453185103409 * tmp_0;
      real_t tmp_34 = 0.058393137863189684 * tmp_10;
      real_t tmp_35 = 0.041425537809186785 * tmp_12;
      real_t tmp_36 = 0.041425537809186785 * tmp_14;
      real_t tmp_37 = 0.025422453185103409 * tmp_16;
      real_t tmp_38 = 0.058393137863189684 * tmp_18;
      real_t tmp_39 = 0.025422453185103409 * tmp_20;
      real_t tmp_40 = 0.058393137863189684 * tmp_22;
      real_t tmp_41 = 0.041425537809186785 * tmp_24;
      real_t tmp_42 = 0.041425537809186785 * tmp_26;
      real_t tmp_43 = 0.041425537809186785 * tmp_28;
      real_t tmp_44 = 0.041425537809186785 * tmp_30;
      real_t tmp_45 = tmp_5 * tmp_7;
      real_t tmp_46 = 0.001603877517404526;
      real_t tmp_47 = 0.012856517194473826;
      real_t tmp_48 = 0.026367458342995378;
      real_t tmp_49 = 0.014556635278230808;
      real_t tmp_50 = 0.022214698150294358;
      real_t tmp_51 = 0.029279867306728068;
      real_t tmp_52 = 0.001603877517404526;
      real_t tmp_53 = 0.014556635278230808;
      real_t tmp_54 = 0.012856517194473826;
      real_t tmp_55 = 0.0022015622717175805;
      real_t tmp_56 = 0.026367458342995378;
      real_t tmp_57 = 0.0022015622717175805;
      real_t tmp_58 = 0.022214698150294358;
      real_t tmp_59 = 0.0022015622717175805;
      real_t tmp_60 = 0.012856517194473826;
      real_t tmp_61 = 0.029279867306728068;
      real_t tmp_62 = 0.001603877517404526;
      real_t tmp_63 = 0.014556635278230808;
      real_t tmp_64 = 0.001603877517404526;
      real_t tmp_65 = 0.014556635278230808;
      real_t tmp_66 = 0.026367458342995378;
      real_t tmp_67 = 0.012856517194473826;
      real_t tmp_68 = 0.0022015622717175805;
      real_t tmp_69 = 0.026367458342995378;
      real_t a_0_0  = tmp_0 * tmp_9 + tmp_10 * tmp_11 + tmp_12 * tmp_13 + tmp_14 * tmp_15 + tmp_16 * tmp_17 + tmp_18 * tmp_19 +
                     tmp_20 * tmp_21 + tmp_22 * tmp_23 + tmp_24 * tmp_25 + tmp_26 * tmp_27 + tmp_28 * tmp_29 + tmp_30 * tmp_31;
      real_t a_0_1 = -tmp_32 * tmp_33 - tmp_32 * tmp_34 - tmp_32 * tmp_35 - tmp_32 * tmp_36 - tmp_32 * tmp_37 - tmp_32 * tmp_38 -
                     tmp_32 * tmp_39 - tmp_32 * tmp_40 - tmp_32 * tmp_41 - tmp_32 * tmp_42 - tmp_32 * tmp_43 - tmp_32 * tmp_44;
      real_t a_0_2 = -tmp_33 * tmp_45 - tmp_34 * tmp_45 - tmp_35 * tmp_45 - tmp_36 * tmp_45 - tmp_37 * tmp_45 - tmp_38 * tmp_45 -
                     tmp_39 * tmp_45 - tmp_40 * tmp_45 - tmp_41 * tmp_45 - tmp_42 * tmp_45 - tmp_43 * tmp_45 - tmp_44 * tmp_45;
      real_t a_1_0 = 0.24928674517091043 * tmp_11 + 0.31035245103378439 * tmp_13 + 0.63650249912139867 * tmp_15 +
                     0.87382197101699555 * tmp_17 + 0.50142650965817914 * tmp_19 + 0.063089014491502227 * tmp_21 +
                     0.24928674517091043 * tmp_23 + 0.31035245103378439 * tmp_25 + 0.053145049844816938 * tmp_27 +
                     0.63650249912139867 * tmp_29 + 0.053145049844816938 * tmp_31 + 0.063089014491502227 * tmp_9;
      real_t a_1_1 = -tmp_32 * tmp_46 - tmp_32 * tmp_47 - tmp_32 * tmp_48 - tmp_32 * tmp_49 - tmp_32 * tmp_50 - tmp_32 * tmp_51 -
                     tmp_32 * tmp_52 - tmp_32 * tmp_53 - tmp_32 * tmp_54 - tmp_32 * tmp_55 - tmp_32 * tmp_56 - tmp_32 * tmp_57;
      real_t a_1_2 = -tmp_45 * tmp_46 - tmp_45 * tmp_47 - tmp_45 * tmp_48 - tmp_45 * tmp_49 - tmp_45 * tmp_50 - tmp_45 * tmp_51 -
                     tmp_45 * tmp_52 - tmp_45 * tmp_53 - tmp_45 * tmp_54 - tmp_45 * tmp_55 - tmp_45 * tmp_56 - tmp_45 * tmp_57;
      real_t a_2_0 = 0.50142650965817914 * tmp_11 + 0.053145049844816938 * tmp_13 + 0.31035245103378439 * tmp_15 +
                     0.063089014491502227 * tmp_17 + 0.24928674517091043 * tmp_19 + 0.063089014491502227 * tmp_21 +
                     0.24928674517091043 * tmp_23 + 0.63650249912139867 * tmp_25 + 0.31035245103378439 * tmp_27 +
                     0.053145049844816938 * tmp_29 + 0.63650249912139867 * tmp_31 + 0.87382197101699555 * tmp_9;
      real_t a_2_1 = -tmp_32 * tmp_58 - tmp_32 * tmp_59 - tmp_32 * tmp_60 - tmp_32 * tmp_61 - tmp_32 * tmp_62 - tmp_32 * tmp_63 -
                     tmp_32 * tmp_64 - tmp_32 * tmp_65 - tmp_32 * tmp_66 - tmp_32 * tmp_67 - tmp_32 * tmp_68 - tmp_32 * tmp_69;
      real_t a_2_2 = -tmp_45 * tmp_58 - tmp_45 * tmp_59 - tmp_45 * tmp_60 - tmp_45 * tmp_61 - tmp_45 * tmp_62 - tmp_45 * tmp_63 -
                     tmp_45 * tmp_64 - tmp_45 * tmp_65 - tmp_45 * tmp_66 - tmp_45 * tmp_67 - tmp_45 * tmp_68 - tmp_45 * tmp_69;
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + tmp_2;
      real_t tmp_7  = tmp_4 * ( 0.046910077030668018 * tmp_5 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = 0.5 * p_affine_10_1 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_19 = 0.11846344252809471 * tmp_18;
      real_t tmp_20 = tmp_4 * ( 0.23076534494715845 * tmp_5 + tmp_6 );
      real_t tmp_21 = tmp_1 * tmp_20;
      real_t tmp_22 = tmp_20 * tmp_9;
      real_t tmp_23 = tmp_4 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_24 = tmp_23 * tmp_3;
      real_t tmp_25 = tmp_15 * tmp_23;
      real_t tmp_26 = -tmp_21 - tmp_22 - tmp_24 - tmp_25 + 1;
      real_t tmp_27 = 0.2393143352496831 * tmp_18;
      real_t tmp_28 = tmp_4 * ( 0.5 * tmp_5 + tmp_6 );
      real_t tmp_29 = tmp_1 * tmp_28;
      real_t tmp_30 = tmp_28 * tmp_9;
      real_t tmp_31 = tmp_4 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_32 = tmp_3 * tmp_31;
      real_t tmp_33 = tmp_15 * tmp_31;
      real_t tmp_34 = -tmp_29 - tmp_30 - tmp_32 - tmp_33 + 1;
      real_t tmp_35 = 0.2844444444444445 * tmp_18;
      real_t tmp_36 = tmp_4 * ( 0.7692346550528415 * tmp_5 + tmp_6 );
      real_t tmp_37 = tmp_1 * tmp_36;
      real_t tmp_38 = tmp_36 * tmp_9;
      real_t tmp_39 = tmp_4 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_40 = tmp_3 * tmp_39;
      real_t tmp_41 = tmp_15 * tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = 0.2393143352496831 * tmp_18;
      real_t tmp_44 = tmp_4 * ( 0.95308992296933193 * tmp_5 + tmp_6 );
      real_t tmp_45 = tmp_1 * tmp_44;
      real_t tmp_46 = tmp_44 * tmp_9;
      real_t tmp_47 = tmp_4 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_48 = tmp_3 * tmp_47;
      real_t tmp_49 = tmp_15 * tmp_47;
      real_t tmp_50 = -tmp_45 - tmp_46 - tmp_48 - tmp_49 + 1;
      real_t tmp_51 = 0.11846344252809471 * tmp_18;
      real_t tmp_52 = tmp_10 + tmp_14;
      real_t tmp_53 = tmp_17 * tmp_19;
      real_t tmp_54 = tmp_22 + tmp_24;
      real_t tmp_55 = tmp_26 * tmp_27;
      real_t tmp_56 = tmp_30 + tmp_32;
      real_t tmp_57 = tmp_34 * tmp_35;
      real_t tmp_58 = tmp_38 + tmp_40;
      real_t tmp_59 = tmp_42 * tmp_43;
      real_t tmp_60 = tmp_46 + tmp_48;
      real_t tmp_61 = tmp_50 * tmp_51;
      real_t tmp_62 = tmp_52 * tmp_53 + tmp_54 * tmp_55 + tmp_56 * tmp_57 + tmp_58 * tmp_59 + tmp_60 * tmp_61;
      real_t tmp_63 = tmp_16 + tmp_8;
      real_t tmp_64 = tmp_21 + tmp_25;
      real_t tmp_65 = tmp_29 + tmp_33;
      real_t tmp_66 = tmp_37 + tmp_41;
      real_t tmp_67 = tmp_45 + tmp_49;
      real_t tmp_68 = tmp_53 * tmp_63 + tmp_55 * tmp_64 + tmp_57 * tmp_65 + tmp_59 * tmp_66 + tmp_61 * tmp_67;
      real_t tmp_69 = tmp_19 * tmp_52 * tmp_63 + tmp_27 * tmp_54 * tmp_64 + tmp_35 * tmp_56 * tmp_65 + tmp_43 * tmp_58 * tmp_66 +
                      tmp_51 * tmp_60 * tmp_67;
      real_t a_0_0 = ( tmp_17 * tmp_17 ) * tmp_19 + ( tmp_26 * tmp_26 ) * tmp_27 + ( tmp_34 * tmp_34 ) * tmp_35 +
                     ( tmp_42 * tmp_42 ) * tmp_43 + ( tmp_50 * tmp_50 ) * tmp_51;
      real_t a_0_1 = tmp_62;
      real_t a_0_2 = tmp_68;
      real_t a_1_0 = tmp_62;
      real_t a_1_1 = tmp_19 * ( tmp_52 * tmp_52 ) + tmp_27 * ( tmp_54 * tmp_54 ) + tmp_35 * ( tmp_56 * tmp_56 ) +
                     tmp_43 * ( tmp_58 * tmp_58 ) + tmp_51 * ( tmp_60 * tmp_60 );
      real_t a_1_2 = tmp_69;
      real_t a_2_0 = tmp_68;
      real_t a_2_1 = tmp_69;
      real_t a_2_2 = tmp_19 * ( tmp_63 * tmp_63 ) + tmp_27 * ( tmp_64 * tmp_64 ) + tmp_35 * ( tmp_65 * tmp_65 ) +
                     tmp_43 * ( tmp_66 * tmp_66 ) + tmp_51 * ( tmp_67 * tmp_67 );
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

      real_t tmp_0   = -p_affine_3_0;
      real_t tmp_1   = p_affine_4_0 + tmp_0;
      real_t tmp_2   = -p_affine_3_1;
      real_t tmp_3   = p_affine_5_1 + tmp_2;
      real_t tmp_4   = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_4_1 + tmp_2 ) * ( p_affine_5_0 + tmp_0 ) );
      real_t tmp_5   = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6   = p_affine_6_1 + 0.046910077030668018 * tmp_5;
      real_t tmp_7   = tmp_4 * ( tmp_2 + tmp_6 );
      real_t tmp_8   = tmp_1 * tmp_7;
      real_t tmp_9   = p_affine_3_0 - p_affine_5_0;
      real_t tmp_10  = tmp_7 * tmp_9;
      real_t tmp_11  = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12  = p_affine_6_0 + 0.046910077030668018 * tmp_11;
      real_t tmp_13  = tmp_4 * ( tmp_0 + tmp_12 );
      real_t tmp_14  = tmp_13 * tmp_3;
      real_t tmp_15  = p_affine_3_1 - p_affine_4_1;
      real_t tmp_16  = tmp_13 * tmp_15;
      real_t tmp_17  = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18  = -p_affine_0_0;
      real_t tmp_19  = p_affine_1_0 + tmp_18;
      real_t tmp_20  = -p_affine_0_1;
      real_t tmp_21  = p_affine_2_1 + tmp_20;
      real_t tmp_22  = 1.0 / ( tmp_19 * tmp_21 - ( p_affine_1_1 + tmp_20 ) * ( p_affine_2_0 + tmp_18 ) );
      real_t tmp_23  = tmp_22 * ( tmp_20 + tmp_6 );
      real_t tmp_24  = tmp_19 * tmp_23;
      real_t tmp_25  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_26  = tmp_23 * tmp_25;
      real_t tmp_27  = tmp_22 * ( tmp_12 + tmp_18 );
      real_t tmp_28  = tmp_21 * tmp_27;
      real_t tmp_29  = p_affine_0_1 - p_affine_1_1;
      real_t tmp_30  = tmp_27 * tmp_29;
      real_t tmp_31  = 0.5 * p_affine_10_1 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_32  = 0.11846344252809471 * tmp_31;
      real_t tmp_33  = tmp_32 * ( -tmp_24 - tmp_26 - tmp_28 - tmp_30 + 1 );
      real_t tmp_34  = p_affine_6_1 + 0.23076534494715845 * tmp_5;
      real_t tmp_35  = tmp_4 * ( tmp_2 + tmp_34 );
      real_t tmp_36  = tmp_1 * tmp_35;
      real_t tmp_37  = tmp_35 * tmp_9;
      real_t tmp_38  = p_affine_6_0 + 0.23076534494715845 * tmp_11;
      real_t tmp_39  = tmp_4 * ( tmp_0 + tmp_38 );
      real_t tmp_40  = tmp_3 * tmp_39;
      real_t tmp_41  = tmp_15 * tmp_39;
      real_t tmp_42  = -tmp_36 - tmp_37 - tmp_40 - tmp_41 + 1;
      real_t tmp_43  = tmp_22 * ( tmp_20 + tmp_34 );
      real_t tmp_44  = tmp_19 * tmp_43;
      real_t tmp_45  = tmp_25 * tmp_43;
      real_t tmp_46  = tmp_22 * ( tmp_18 + tmp_38 );
      real_t tmp_47  = tmp_21 * tmp_46;
      real_t tmp_48  = tmp_29 * tmp_46;
      real_t tmp_49  = 0.2393143352496831 * tmp_31;
      real_t tmp_50  = tmp_49 * ( -tmp_44 - tmp_45 - tmp_47 - tmp_48 + 1 );
      real_t tmp_51  = p_affine_6_1 + 0.5 * tmp_5;
      real_t tmp_52  = tmp_4 * ( tmp_2 + tmp_51 );
      real_t tmp_53  = tmp_1 * tmp_52;
      real_t tmp_54  = tmp_52 * tmp_9;
      real_t tmp_55  = p_affine_6_0 + 0.5 * tmp_11;
      real_t tmp_56  = tmp_4 * ( tmp_0 + tmp_55 );
      real_t tmp_57  = tmp_3 * tmp_56;
      real_t tmp_58  = tmp_15 * tmp_56;
      real_t tmp_59  = -tmp_53 - tmp_54 - tmp_57 - tmp_58 + 1;
      real_t tmp_60  = tmp_22 * ( tmp_20 + tmp_51 );
      real_t tmp_61  = tmp_19 * tmp_60;
      real_t tmp_62  = tmp_25 * tmp_60;
      real_t tmp_63  = tmp_22 * ( tmp_18 + tmp_55 );
      real_t tmp_64  = tmp_21 * tmp_63;
      real_t tmp_65  = tmp_29 * tmp_63;
      real_t tmp_66  = 0.2844444444444445 * tmp_31;
      real_t tmp_67  = tmp_66 * ( -tmp_61 - tmp_62 - tmp_64 - tmp_65 + 1 );
      real_t tmp_68  = p_affine_6_1 + 0.7692346550528415 * tmp_5;
      real_t tmp_69  = tmp_4 * ( tmp_2 + tmp_68 );
      real_t tmp_70  = tmp_1 * tmp_69;
      real_t tmp_71  = tmp_69 * tmp_9;
      real_t tmp_72  = p_affine_6_0 + 0.7692346550528415 * tmp_11;
      real_t tmp_73  = tmp_4 * ( tmp_0 + tmp_72 );
      real_t tmp_74  = tmp_3 * tmp_73;
      real_t tmp_75  = tmp_15 * tmp_73;
      real_t tmp_76  = -tmp_70 - tmp_71 - tmp_74 - tmp_75 + 1;
      real_t tmp_77  = tmp_22 * ( tmp_20 + tmp_68 );
      real_t tmp_78  = tmp_19 * tmp_77;
      real_t tmp_79  = tmp_25 * tmp_77;
      real_t tmp_80  = tmp_22 * ( tmp_18 + tmp_72 );
      real_t tmp_81  = tmp_21 * tmp_80;
      real_t tmp_82  = tmp_29 * tmp_80;
      real_t tmp_83  = 0.2393143352496831 * tmp_31;
      real_t tmp_84  = tmp_83 * ( -tmp_78 - tmp_79 - tmp_81 - tmp_82 + 1 );
      real_t tmp_85  = p_affine_6_1 + 0.95308992296933193 * tmp_5;
      real_t tmp_86  = tmp_4 * ( tmp_2 + tmp_85 );
      real_t tmp_87  = tmp_1 * tmp_86;
      real_t tmp_88  = tmp_86 * tmp_9;
      real_t tmp_89  = p_affine_6_0 + 0.95308992296933193 * tmp_11;
      real_t tmp_90  = tmp_4 * ( tmp_0 + tmp_89 );
      real_t tmp_91  = tmp_3 * tmp_90;
      real_t tmp_92  = tmp_15 * tmp_90;
      real_t tmp_93  = -tmp_87 - tmp_88 - tmp_91 - tmp_92 + 1;
      real_t tmp_94  = tmp_22 * ( tmp_20 + tmp_85 );
      real_t tmp_95  = tmp_19 * tmp_94;
      real_t tmp_96  = tmp_25 * tmp_94;
      real_t tmp_97  = tmp_22 * ( tmp_18 + tmp_89 );
      real_t tmp_98  = tmp_21 * tmp_97;
      real_t tmp_99  = tmp_29 * tmp_97;
      real_t tmp_100 = 0.11846344252809471 * tmp_31;
      real_t tmp_101 = tmp_100 * ( -tmp_95 - tmp_96 - tmp_98 - tmp_99 + 1 );
      real_t tmp_102 = tmp_10 + tmp_14;
      real_t tmp_103 = tmp_37 + tmp_40;
      real_t tmp_104 = tmp_54 + tmp_57;
      real_t tmp_105 = tmp_71 + tmp_74;
      real_t tmp_106 = tmp_88 + tmp_91;
      real_t tmp_107 = tmp_16 + tmp_8;
      real_t tmp_108 = tmp_36 + tmp_41;
      real_t tmp_109 = tmp_53 + tmp_58;
      real_t tmp_110 = tmp_70 + tmp_75;
      real_t tmp_111 = tmp_87 + tmp_92;
      real_t tmp_112 = tmp_32 * ( tmp_26 + tmp_28 );
      real_t tmp_113 = tmp_49 * ( tmp_45 + tmp_47 );
      real_t tmp_114 = tmp_66 * ( tmp_62 + tmp_64 );
      real_t tmp_115 = tmp_83 * ( tmp_79 + tmp_81 );
      real_t tmp_116 = tmp_100 * ( tmp_96 + tmp_98 );
      real_t tmp_117 = tmp_32 * ( tmp_24 + tmp_30 );
      real_t tmp_118 = tmp_49 * ( tmp_44 + tmp_48 );
      real_t tmp_119 = tmp_66 * ( tmp_61 + tmp_65 );
      real_t tmp_120 = tmp_83 * ( tmp_78 + tmp_82 );
      real_t tmp_121 = tmp_100 * ( tmp_95 + tmp_99 );
      real_t a_0_0   = -tmp_101 * tmp_93 - tmp_17 * tmp_33 - tmp_42 * tmp_50 - tmp_59 * tmp_67 - tmp_76 * tmp_84;
      real_t a_0_1   = -tmp_101 * tmp_106 - tmp_102 * tmp_33 - tmp_103 * tmp_50 - tmp_104 * tmp_67 - tmp_105 * tmp_84;
      real_t a_0_2   = -tmp_101 * tmp_111 - tmp_107 * tmp_33 - tmp_108 * tmp_50 - tmp_109 * tmp_67 - tmp_110 * tmp_84;
      real_t a_1_0   = -tmp_112 * tmp_17 - tmp_113 * tmp_42 - tmp_114 * tmp_59 - tmp_115 * tmp_76 - tmp_116 * tmp_93;
      real_t a_1_1   = -tmp_102 * tmp_112 - tmp_103 * tmp_113 - tmp_104 * tmp_114 - tmp_105 * tmp_115 - tmp_106 * tmp_116;
      real_t a_1_2   = -tmp_107 * tmp_112 - tmp_108 * tmp_113 - tmp_109 * tmp_114 - tmp_110 * tmp_115 - tmp_111 * tmp_116;
      real_t a_2_0   = -tmp_117 * tmp_17 - tmp_118 * tmp_42 - tmp_119 * tmp_59 - tmp_120 * tmp_76 - tmp_121 * tmp_93;
      real_t a_2_1   = -tmp_102 * tmp_117 - tmp_103 * tmp_118 - tmp_104 * tmp_119 - tmp_105 * tmp_120 - tmp_106 * tmp_121;
      real_t a_2_2   = -tmp_107 * tmp_117 - tmp_108 * tmp_118 - tmp_109 * tmp_119 - tmp_110 * tmp_120 - tmp_111 * tmp_121;
      elMat( 0, 0 )  = a_0_0;
      elMat( 0, 1 )  = a_0_1;
      elMat( 0, 2 )  = a_0_2;
      elMat( 1, 0 )  = a_1_0;
      elMat( 1, 1 )  = a_1_1;
      elMat( 1, 2 )  = a_1_2;
      elMat( 2, 0 )  = a_2_0;
      elMat( 2, 1 )  = a_2_1;
      elMat( 2, 2 )  = a_2_2;
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

      real_t tmp_0  = 0.063089014491502282;
      real_t tmp_1  = -p_affine_0_1;
      real_t tmp_2  = p_affine_2_1 + tmp_1;
      real_t tmp_3  = -p_affine_0_0;
      real_t tmp_4  = 1.0 / ( tmp_2 * ( p_affine_1_0 + tmp_3 ) - ( p_affine_1_1 + tmp_1 ) * ( p_affine_2_0 + tmp_3 ) );
      real_t tmp_5  = tmp_2 * tmp_4;
      real_t tmp_6  = tmp_4 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_7  = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                               p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_8  = tmp_7 * ( tmp_5 + tmp_6 );
      real_t tmp_9  = 0.025422453185103409 * tmp_8;
      real_t tmp_10 = 0.24928674517091043;
      real_t tmp_11 = 0.058393137863189684 * tmp_8;
      real_t tmp_12 = 0.63650249912139867;
      real_t tmp_13 = 0.041425537809186785 * tmp_8;
      real_t tmp_14 = 0.053145049844816938;
      real_t tmp_15 = 0.041425537809186785 * tmp_8;
      real_t tmp_16 = 0.063089014491502227;
      real_t tmp_17 = 0.025422453185103409 * tmp_8;
      real_t tmp_18 = 0.24928674517091043;
      real_t tmp_19 = 0.058393137863189684 * tmp_8;
      real_t tmp_20 = 0.87382197101699566;
      real_t tmp_21 = 0.025422453185103409 * tmp_8;
      real_t tmp_22 = 0.50142650965817914;
      real_t tmp_23 = 0.058393137863189684 * tmp_8;
      real_t tmp_24 = 0.053145049844816938;
      real_t tmp_25 = 0.041425537809186785 * tmp_8;
      real_t tmp_26 = 0.63650249912139867;
      real_t tmp_27 = 0.041425537809186785 * tmp_8;
      real_t tmp_28 = 0.31035245103378439;
      real_t tmp_29 = 0.041425537809186785 * tmp_8;
      real_t tmp_30 = 0.31035245103378439;
      real_t tmp_31 = 0.041425537809186785 * tmp_8;
      real_t tmp_32 = tmp_5 * tmp_7;
      real_t tmp_33 = 0.025422453185103409 * tmp_32;
      real_t tmp_34 = 0.058393137863189684 * tmp_32;
      real_t tmp_35 = 0.041425537809186785 * tmp_32;
      real_t tmp_36 = 0.041425537809186785 * tmp_32;
      real_t tmp_37 = 0.025422453185103409 * tmp_32;
      real_t tmp_38 = 0.058393137863189684 * tmp_32;
      real_t tmp_39 = 0.025422453185103409 * tmp_32;
      real_t tmp_40 = 0.058393137863189684 * tmp_32;
      real_t tmp_41 = 0.041425537809186785 * tmp_32;
      real_t tmp_42 = 0.041425537809186785 * tmp_32;
      real_t tmp_43 = 0.041425537809186785 * tmp_32;
      real_t tmp_44 = 0.041425537809186785 * tmp_32;
      real_t tmp_45 = tmp_6 * tmp_7;
      real_t tmp_46 = 0.025422453185103409 * tmp_45;
      real_t tmp_47 = 0.058393137863189684 * tmp_45;
      real_t tmp_48 = 0.041425537809186785 * tmp_45;
      real_t tmp_49 = 0.041425537809186785 * tmp_45;
      real_t tmp_50 = 0.025422453185103409 * tmp_45;
      real_t tmp_51 = 0.058393137863189684 * tmp_45;
      real_t tmp_52 = 0.025422453185103409 * tmp_45;
      real_t tmp_53 = 0.058393137863189684 * tmp_45;
      real_t tmp_54 = 0.041425537809186785 * tmp_45;
      real_t tmp_55 = 0.041425537809186785 * tmp_45;
      real_t tmp_56 = 0.041425537809186785 * tmp_45;
      real_t tmp_57 = 0.041425537809186785 * tmp_45;
      real_t a_0_0  = tmp_0 * tmp_9 + tmp_10 * tmp_11 + tmp_12 * tmp_13 + tmp_14 * tmp_15 + tmp_16 * tmp_17 + tmp_18 * tmp_19 +
                     tmp_20 * tmp_21 + tmp_22 * tmp_23 + tmp_24 * tmp_25 + tmp_26 * tmp_27 + tmp_28 * tmp_29 + tmp_30 * tmp_31;
      real_t a_0_1 = 0.24928674517091043 * tmp_11 + 0.31035245103378439 * tmp_13 + 0.63650249912139867 * tmp_15 +
                     0.87382197101699555 * tmp_17 + 0.50142650965817914 * tmp_19 + 0.063089014491502227 * tmp_21 +
                     0.24928674517091043 * tmp_23 + 0.31035245103378439 * tmp_25 + 0.053145049844816938 * tmp_27 +
                     0.63650249912139867 * tmp_29 + 0.053145049844816938 * tmp_31 + 0.063089014491502227 * tmp_9;
      real_t a_0_2 = 0.50142650965817914 * tmp_11 + 0.053145049844816938 * tmp_13 + 0.31035245103378439 * tmp_15 +
                     0.063089014491502227 * tmp_17 + 0.24928674517091043 * tmp_19 + 0.063089014491502227 * tmp_21 +
                     0.24928674517091043 * tmp_23 + 0.63650249912139867 * tmp_25 + 0.31035245103378439 * tmp_27 +
                     0.053145049844816938 * tmp_29 + 0.63650249912139867 * tmp_31 + 0.87382197101699555 * tmp_9;
      real_t a_1_0 = -tmp_0 * tmp_33 - tmp_10 * tmp_34 - tmp_12 * tmp_35 - tmp_14 * tmp_36 - tmp_16 * tmp_37 - tmp_18 * tmp_38 -
                     tmp_20 * tmp_39 - tmp_22 * tmp_40 - tmp_24 * tmp_41 - tmp_26 * tmp_42 - tmp_28 * tmp_43 - tmp_30 * tmp_44;
      real_t a_1_1 = -0.063089014491502227 * tmp_33 - 0.24928674517091043 * tmp_34 - 0.31035245103378439 * tmp_35 -
                     0.63650249912139867 * tmp_36 - 0.87382197101699555 * tmp_37 - 0.50142650965817914 * tmp_38 -
                     0.063089014491502227 * tmp_39 - 0.24928674517091043 * tmp_40 - 0.31035245103378439 * tmp_41 -
                     0.053145049844816938 * tmp_42 - 0.63650249912139867 * tmp_43 - 0.053145049844816938 * tmp_44;
      real_t a_1_2 = -0.87382197101699555 * tmp_33 - 0.50142650965817914 * tmp_34 - 0.053145049844816938 * tmp_35 -
                     0.31035245103378439 * tmp_36 - 0.063089014491502227 * tmp_37 - 0.24928674517091043 * tmp_38 -
                     0.063089014491502227 * tmp_39 - 0.24928674517091043 * tmp_40 - 0.63650249912139867 * tmp_41 -
                     0.31035245103378439 * tmp_42 - 0.053145049844816938 * tmp_43 - 0.63650249912139867 * tmp_44;
      real_t a_2_0 = -tmp_0 * tmp_46 - tmp_10 * tmp_47 - tmp_12 * tmp_48 - tmp_14 * tmp_49 - tmp_16 * tmp_50 - tmp_18 * tmp_51 -
                     tmp_20 * tmp_52 - tmp_22 * tmp_53 - tmp_24 * tmp_54 - tmp_26 * tmp_55 - tmp_28 * tmp_56 - tmp_30 * tmp_57;
      real_t a_2_1 = -0.063089014491502227 * tmp_46 - 0.24928674517091043 * tmp_47 - 0.31035245103378439 * tmp_48 -
                     0.63650249912139867 * tmp_49 - 0.87382197101699555 * tmp_50 - 0.50142650965817914 * tmp_51 -
                     0.063089014491502227 * tmp_52 - 0.24928674517091043 * tmp_53 - 0.31035245103378439 * tmp_54 -
                     0.053145049844816938 * tmp_55 - 0.63650249912139867 * tmp_56 - 0.053145049844816938 * tmp_57;
      real_t a_2_2 = -0.87382197101699555 * tmp_46 - 0.50142650965817914 * tmp_47 - 0.053145049844816938 * tmp_48 -
                     0.31035245103378439 * tmp_49 - 0.063089014491502227 * tmp_50 - 0.24928674517091043 * tmp_51 -
                     0.063089014491502227 * tmp_52 - 0.24928674517091043 * tmp_53 - 0.63650249912139867 * tmp_54 -
                     0.31035245103378439 * tmp_55 - 0.053145049844816938 * tmp_56 - 0.63650249912139867 * tmp_57;
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + tmp_2;
      real_t tmp_7  = tmp_4 * ( 0.046910077030668018 * tmp_5 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = 0.5 * p_affine_10_0 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_19 = 0.11846344252809471 * tmp_18;
      real_t tmp_20 = tmp_4 * ( 0.23076534494715845 * tmp_5 + tmp_6 );
      real_t tmp_21 = tmp_1 * tmp_20;
      real_t tmp_22 = tmp_20 * tmp_9;
      real_t tmp_23 = tmp_4 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_24 = tmp_23 * tmp_3;
      real_t tmp_25 = tmp_15 * tmp_23;
      real_t tmp_26 = -tmp_21 - tmp_22 - tmp_24 - tmp_25 + 1;
      real_t tmp_27 = 0.2393143352496831 * tmp_18;
      real_t tmp_28 = tmp_4 * ( 0.5 * tmp_5 + tmp_6 );
      real_t tmp_29 = tmp_1 * tmp_28;
      real_t tmp_30 = tmp_28 * tmp_9;
      real_t tmp_31 = tmp_4 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_32 = tmp_3 * tmp_31;
      real_t tmp_33 = tmp_15 * tmp_31;
      real_t tmp_34 = -tmp_29 - tmp_30 - tmp_32 - tmp_33 + 1;
      real_t tmp_35 = 0.2844444444444445 * tmp_18;
      real_t tmp_36 = tmp_4 * ( 0.7692346550528415 * tmp_5 + tmp_6 );
      real_t tmp_37 = tmp_1 * tmp_36;
      real_t tmp_38 = tmp_36 * tmp_9;
      real_t tmp_39 = tmp_4 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_40 = tmp_3 * tmp_39;
      real_t tmp_41 = tmp_15 * tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = 0.2393143352496831 * tmp_18;
      real_t tmp_44 = tmp_4 * ( 0.95308992296933193 * tmp_5 + tmp_6 );
      real_t tmp_45 = tmp_1 * tmp_44;
      real_t tmp_46 = tmp_44 * tmp_9;
      real_t tmp_47 = tmp_4 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_48 = tmp_3 * tmp_47;
      real_t tmp_49 = tmp_15 * tmp_47;
      real_t tmp_50 = -tmp_45 - tmp_46 - tmp_48 - tmp_49 + 1;
      real_t tmp_51 = 0.11846344252809471 * tmp_18;
      real_t tmp_52 = tmp_10 + tmp_14;
      real_t tmp_53 = tmp_17 * tmp_19;
      real_t tmp_54 = tmp_22 + tmp_24;
      real_t tmp_55 = tmp_26 * tmp_27;
      real_t tmp_56 = tmp_30 + tmp_32;
      real_t tmp_57 = tmp_34 * tmp_35;
      real_t tmp_58 = tmp_38 + tmp_40;
      real_t tmp_59 = tmp_42 * tmp_43;
      real_t tmp_60 = tmp_46 + tmp_48;
      real_t tmp_61 = tmp_50 * tmp_51;
      real_t tmp_62 = tmp_52 * tmp_53 + tmp_54 * tmp_55 + tmp_56 * tmp_57 + tmp_58 * tmp_59 + tmp_60 * tmp_61;
      real_t tmp_63 = tmp_16 + tmp_8;
      real_t tmp_64 = tmp_21 + tmp_25;
      real_t tmp_65 = tmp_29 + tmp_33;
      real_t tmp_66 = tmp_37 + tmp_41;
      real_t tmp_67 = tmp_45 + tmp_49;
      real_t tmp_68 = tmp_53 * tmp_63 + tmp_55 * tmp_64 + tmp_57 * tmp_65 + tmp_59 * tmp_66 + tmp_61 * tmp_67;
      real_t tmp_69 = tmp_19 * tmp_52 * tmp_63 + tmp_27 * tmp_54 * tmp_64 + tmp_35 * tmp_56 * tmp_65 + tmp_43 * tmp_58 * tmp_66 +
                      tmp_51 * tmp_60 * tmp_67;
      real_t a_0_0 = ( tmp_17 * tmp_17 ) * tmp_19 + ( tmp_26 * tmp_26 ) * tmp_27 + ( tmp_34 * tmp_34 ) * tmp_35 +
                     ( tmp_42 * tmp_42 ) * tmp_43 + ( tmp_50 * tmp_50 ) * tmp_51;
      real_t a_0_1 = tmp_62;
      real_t a_0_2 = tmp_68;
      real_t a_1_0 = tmp_62;
      real_t a_1_1 = tmp_19 * ( tmp_52 * tmp_52 ) + tmp_27 * ( tmp_54 * tmp_54 ) + tmp_35 * ( tmp_56 * tmp_56 ) +
                     tmp_43 * ( tmp_58 * tmp_58 ) + tmp_51 * ( tmp_60 * tmp_60 );
      real_t a_1_2 = tmp_69;
      real_t a_2_0 = tmp_68;
      real_t a_2_1 = tmp_69;
      real_t a_2_2 = tmp_19 * ( tmp_63 * tmp_63 ) + tmp_27 * ( tmp_64 * tmp_64 ) + tmp_35 * ( tmp_65 * tmp_65 ) +
                     tmp_43 * ( tmp_66 * tmp_66 ) + tmp_51 * ( tmp_67 * tmp_67 );
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

      real_t tmp_0   = -p_affine_3_0;
      real_t tmp_1   = p_affine_4_0 + tmp_0;
      real_t tmp_2   = -p_affine_3_1;
      real_t tmp_3   = p_affine_5_1 + tmp_2;
      real_t tmp_4   = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_4_1 + tmp_2 ) * ( p_affine_5_0 + tmp_0 ) );
      real_t tmp_5   = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6   = p_affine_6_1 + 0.046910077030668018 * tmp_5;
      real_t tmp_7   = tmp_4 * ( tmp_2 + tmp_6 );
      real_t tmp_8   = tmp_1 * tmp_7;
      real_t tmp_9   = p_affine_3_0 - p_affine_5_0;
      real_t tmp_10  = tmp_7 * tmp_9;
      real_t tmp_11  = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12  = p_affine_6_0 + 0.046910077030668018 * tmp_11;
      real_t tmp_13  = tmp_4 * ( tmp_0 + tmp_12 );
      real_t tmp_14  = tmp_13 * tmp_3;
      real_t tmp_15  = p_affine_3_1 - p_affine_4_1;
      real_t tmp_16  = tmp_13 * tmp_15;
      real_t tmp_17  = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18  = -p_affine_0_0;
      real_t tmp_19  = p_affine_1_0 + tmp_18;
      real_t tmp_20  = -p_affine_0_1;
      real_t tmp_21  = p_affine_2_1 + tmp_20;
      real_t tmp_22  = 1.0 / ( tmp_19 * tmp_21 - ( p_affine_1_1 + tmp_20 ) * ( p_affine_2_0 + tmp_18 ) );
      real_t tmp_23  = tmp_22 * ( tmp_20 + tmp_6 );
      real_t tmp_24  = tmp_19 * tmp_23;
      real_t tmp_25  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_26  = tmp_23 * tmp_25;
      real_t tmp_27  = tmp_22 * ( tmp_12 + tmp_18 );
      real_t tmp_28  = tmp_21 * tmp_27;
      real_t tmp_29  = p_affine_0_1 - p_affine_1_1;
      real_t tmp_30  = tmp_27 * tmp_29;
      real_t tmp_31  = 0.5 * p_affine_10_0 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_32  = 0.11846344252809471 * tmp_31;
      real_t tmp_33  = tmp_32 * ( -tmp_24 - tmp_26 - tmp_28 - tmp_30 + 1 );
      real_t tmp_34  = p_affine_6_1 + 0.23076534494715845 * tmp_5;
      real_t tmp_35  = tmp_4 * ( tmp_2 + tmp_34 );
      real_t tmp_36  = tmp_1 * tmp_35;
      real_t tmp_37  = tmp_35 * tmp_9;
      real_t tmp_38  = p_affine_6_0 + 0.23076534494715845 * tmp_11;
      real_t tmp_39  = tmp_4 * ( tmp_0 + tmp_38 );
      real_t tmp_40  = tmp_3 * tmp_39;
      real_t tmp_41  = tmp_15 * tmp_39;
      real_t tmp_42  = -tmp_36 - tmp_37 - tmp_40 - tmp_41 + 1;
      real_t tmp_43  = tmp_22 * ( tmp_20 + tmp_34 );
      real_t tmp_44  = tmp_19 * tmp_43;
      real_t tmp_45  = tmp_25 * tmp_43;
      real_t tmp_46  = tmp_22 * ( tmp_18 + tmp_38 );
      real_t tmp_47  = tmp_21 * tmp_46;
      real_t tmp_48  = tmp_29 * tmp_46;
      real_t tmp_49  = 0.2393143352496831 * tmp_31;
      real_t tmp_50  = tmp_49 * ( -tmp_44 - tmp_45 - tmp_47 - tmp_48 + 1 );
      real_t tmp_51  = p_affine_6_1 + 0.5 * tmp_5;
      real_t tmp_52  = tmp_4 * ( tmp_2 + tmp_51 );
      real_t tmp_53  = tmp_1 * tmp_52;
      real_t tmp_54  = tmp_52 * tmp_9;
      real_t tmp_55  = p_affine_6_0 + 0.5 * tmp_11;
      real_t tmp_56  = tmp_4 * ( tmp_0 + tmp_55 );
      real_t tmp_57  = tmp_3 * tmp_56;
      real_t tmp_58  = tmp_15 * tmp_56;
      real_t tmp_59  = -tmp_53 - tmp_54 - tmp_57 - tmp_58 + 1;
      real_t tmp_60  = tmp_22 * ( tmp_20 + tmp_51 );
      real_t tmp_61  = tmp_19 * tmp_60;
      real_t tmp_62  = tmp_25 * tmp_60;
      real_t tmp_63  = tmp_22 * ( tmp_18 + tmp_55 );
      real_t tmp_64  = tmp_21 * tmp_63;
      real_t tmp_65  = tmp_29 * tmp_63;
      real_t tmp_66  = 0.2844444444444445 * tmp_31;
      real_t tmp_67  = tmp_66 * ( -tmp_61 - tmp_62 - tmp_64 - tmp_65 + 1 );
      real_t tmp_68  = p_affine_6_1 + 0.7692346550528415 * tmp_5;
      real_t tmp_69  = tmp_4 * ( tmp_2 + tmp_68 );
      real_t tmp_70  = tmp_1 * tmp_69;
      real_t tmp_71  = tmp_69 * tmp_9;
      real_t tmp_72  = p_affine_6_0 + 0.7692346550528415 * tmp_11;
      real_t tmp_73  = tmp_4 * ( tmp_0 + tmp_72 );
      real_t tmp_74  = tmp_3 * tmp_73;
      real_t tmp_75  = tmp_15 * tmp_73;
      real_t tmp_76  = -tmp_70 - tmp_71 - tmp_74 - tmp_75 + 1;
      real_t tmp_77  = tmp_22 * ( tmp_20 + tmp_68 );
      real_t tmp_78  = tmp_19 * tmp_77;
      real_t tmp_79  = tmp_25 * tmp_77;
      real_t tmp_80  = tmp_22 * ( tmp_18 + tmp_72 );
      real_t tmp_81  = tmp_21 * tmp_80;
      real_t tmp_82  = tmp_29 * tmp_80;
      real_t tmp_83  = 0.2393143352496831 * tmp_31;
      real_t tmp_84  = tmp_83 * ( -tmp_78 - tmp_79 - tmp_81 - tmp_82 + 1 );
      real_t tmp_85  = p_affine_6_1 + 0.95308992296933193 * tmp_5;
      real_t tmp_86  = tmp_4 * ( tmp_2 + tmp_85 );
      real_t tmp_87  = tmp_1 * tmp_86;
      real_t tmp_88  = tmp_86 * tmp_9;
      real_t tmp_89  = p_affine_6_0 + 0.95308992296933193 * tmp_11;
      real_t tmp_90  = tmp_4 * ( tmp_0 + tmp_89 );
      real_t tmp_91  = tmp_3 * tmp_90;
      real_t tmp_92  = tmp_15 * tmp_90;
      real_t tmp_93  = -tmp_87 - tmp_88 - tmp_91 - tmp_92 + 1;
      real_t tmp_94  = tmp_22 * ( tmp_20 + tmp_85 );
      real_t tmp_95  = tmp_19 * tmp_94;
      real_t tmp_96  = tmp_25 * tmp_94;
      real_t tmp_97  = tmp_22 * ( tmp_18 + tmp_89 );
      real_t tmp_98  = tmp_21 * tmp_97;
      real_t tmp_99  = tmp_29 * tmp_97;
      real_t tmp_100 = 0.11846344252809471 * tmp_31;
      real_t tmp_101 = tmp_100 * ( -tmp_95 - tmp_96 - tmp_98 - tmp_99 + 1 );
      real_t tmp_102 = tmp_10 + tmp_14;
      real_t tmp_103 = tmp_37 + tmp_40;
      real_t tmp_104 = tmp_54 + tmp_57;
      real_t tmp_105 = tmp_71 + tmp_74;
      real_t tmp_106 = tmp_88 + tmp_91;
      real_t tmp_107 = tmp_16 + tmp_8;
      real_t tmp_108 = tmp_36 + tmp_41;
      real_t tmp_109 = tmp_53 + tmp_58;
      real_t tmp_110 = tmp_70 + tmp_75;
      real_t tmp_111 = tmp_87 + tmp_92;
      real_t tmp_112 = tmp_32 * ( tmp_26 + tmp_28 );
      real_t tmp_113 = tmp_49 * ( tmp_45 + tmp_47 );
      real_t tmp_114 = tmp_66 * ( tmp_62 + tmp_64 );
      real_t tmp_115 = tmp_83 * ( tmp_79 + tmp_81 );
      real_t tmp_116 = tmp_100 * ( tmp_96 + tmp_98 );
      real_t tmp_117 = tmp_32 * ( tmp_24 + tmp_30 );
      real_t tmp_118 = tmp_49 * ( tmp_44 + tmp_48 );
      real_t tmp_119 = tmp_66 * ( tmp_61 + tmp_65 );
      real_t tmp_120 = tmp_83 * ( tmp_78 + tmp_82 );
      real_t tmp_121 = tmp_100 * ( tmp_95 + tmp_99 );
      real_t a_0_0   = tmp_101 * tmp_93 + tmp_17 * tmp_33 + tmp_42 * tmp_50 + tmp_59 * tmp_67 + tmp_76 * tmp_84;
      real_t a_0_1   = tmp_101 * tmp_106 + tmp_102 * tmp_33 + tmp_103 * tmp_50 + tmp_104 * tmp_67 + tmp_105 * tmp_84;
      real_t a_0_2   = tmp_101 * tmp_111 + tmp_107 * tmp_33 + tmp_108 * tmp_50 + tmp_109 * tmp_67 + tmp_110 * tmp_84;
      real_t a_1_0   = tmp_112 * tmp_17 + tmp_113 * tmp_42 + tmp_114 * tmp_59 + tmp_115 * tmp_76 + tmp_116 * tmp_93;
      real_t a_1_1   = tmp_102 * tmp_112 + tmp_103 * tmp_113 + tmp_104 * tmp_114 + tmp_105 * tmp_115 + tmp_106 * tmp_116;
      real_t a_1_2   = tmp_107 * tmp_112 + tmp_108 * tmp_113 + tmp_109 * tmp_114 + tmp_110 * tmp_115 + tmp_111 * tmp_116;
      real_t a_2_0   = tmp_117 * tmp_17 + tmp_118 * tmp_42 + tmp_119 * tmp_59 + tmp_120 * tmp_76 + tmp_121 * tmp_93;
      real_t a_2_1   = tmp_102 * tmp_117 + tmp_103 * tmp_118 + tmp_104 * tmp_119 + tmp_105 * tmp_120 + tmp_106 * tmp_121;
      real_t a_2_2   = tmp_107 * tmp_117 + tmp_108 * tmp_118 + tmp_109 * tmp_119 + tmp_110 * tmp_120 + tmp_111 * tmp_121;
      elMat( 0, 0 )  = a_0_0;
      elMat( 0, 1 )  = a_0_1;
      elMat( 0, 2 )  = a_0_2;
      elMat( 1, 0 )  = a_1_0;
      elMat( 1, 1 )  = a_1_1;
      elMat( 1, 2 )  = a_1_2;
      elMat( 2, 0 )  = a_2_0;
      elMat( 2, 1 )  = a_2_1;
      elMat( 2, 2 )  = a_2_2;
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + tmp_2;
      real_t tmp_7  = tmp_4 * ( 0.046910077030668018 * tmp_5 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = p_affine_10_0 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_19 = 0.11846344252809471 * tmp_18;
      real_t tmp_20 = tmp_4 * ( 0.23076534494715845 * tmp_5 + tmp_6 );
      real_t tmp_21 = tmp_1 * tmp_20;
      real_t tmp_22 = tmp_20 * tmp_9;
      real_t tmp_23 = tmp_4 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_24 = tmp_23 * tmp_3;
      real_t tmp_25 = tmp_15 * tmp_23;
      real_t tmp_26 = -tmp_21 - tmp_22 - tmp_24 - tmp_25 + 1;
      real_t tmp_27 = 0.2393143352496831 * tmp_18;
      real_t tmp_28 = tmp_4 * ( 0.5 * tmp_5 + tmp_6 );
      real_t tmp_29 = tmp_1 * tmp_28;
      real_t tmp_30 = tmp_28 * tmp_9;
      real_t tmp_31 = tmp_4 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_32 = tmp_3 * tmp_31;
      real_t tmp_33 = tmp_15 * tmp_31;
      real_t tmp_34 = -tmp_29 - tmp_30 - tmp_32 - tmp_33 + 1;
      real_t tmp_35 = 0.2844444444444445 * tmp_18;
      real_t tmp_36 = tmp_4 * ( 0.7692346550528415 * tmp_5 + tmp_6 );
      real_t tmp_37 = tmp_1 * tmp_36;
      real_t tmp_38 = tmp_36 * tmp_9;
      real_t tmp_39 = tmp_4 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_40 = tmp_3 * tmp_39;
      real_t tmp_41 = tmp_15 * tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = 0.2393143352496831 * tmp_18;
      real_t tmp_44 = tmp_4 * ( 0.95308992296933193 * tmp_5 + tmp_6 );
      real_t tmp_45 = tmp_1 * tmp_44;
      real_t tmp_46 = tmp_44 * tmp_9;
      real_t tmp_47 = tmp_4 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_48 = tmp_3 * tmp_47;
      real_t tmp_49 = tmp_15 * tmp_47;
      real_t tmp_50 = -tmp_45 - tmp_46 - tmp_48 - tmp_49 + 1;
      real_t tmp_51 = 0.11846344252809471 * tmp_18;
      real_t tmp_52 = tmp_10 + tmp_14;
      real_t tmp_53 = tmp_17 * tmp_19;
      real_t tmp_54 = tmp_22 + tmp_24;
      real_t tmp_55 = tmp_26 * tmp_27;
      real_t tmp_56 = tmp_30 + tmp_32;
      real_t tmp_57 = tmp_34 * tmp_35;
      real_t tmp_58 = tmp_38 + tmp_40;
      real_t tmp_59 = tmp_42 * tmp_43;
      real_t tmp_60 = tmp_46 + tmp_48;
      real_t tmp_61 = tmp_50 * tmp_51;
      real_t tmp_62 = tmp_52 * tmp_53 + tmp_54 * tmp_55 + tmp_56 * tmp_57 + tmp_58 * tmp_59 + tmp_60 * tmp_61;
      real_t tmp_63 = tmp_16 + tmp_8;
      real_t tmp_64 = tmp_21 + tmp_25;
      real_t tmp_65 = tmp_29 + tmp_33;
      real_t tmp_66 = tmp_37 + tmp_41;
      real_t tmp_67 = tmp_45 + tmp_49;
      real_t tmp_68 = tmp_53 * tmp_63 + tmp_55 * tmp_64 + tmp_57 * tmp_65 + tmp_59 * tmp_66 + tmp_61 * tmp_67;
      real_t tmp_69 = tmp_19 * tmp_52 * tmp_63 + tmp_27 * tmp_54 * tmp_64 + tmp_35 * tmp_56 * tmp_65 + tmp_43 * tmp_58 * tmp_66 +
                      tmp_51 * tmp_60 * tmp_67;
      real_t a_0_0 = ( tmp_17 * tmp_17 ) * tmp_19 + ( tmp_26 * tmp_26 ) * tmp_27 + ( tmp_34 * tmp_34 ) * tmp_35 +
                     ( tmp_42 * tmp_42 ) * tmp_43 + ( tmp_50 * tmp_50 ) * tmp_51;
      real_t a_0_1 = tmp_62;
      real_t a_0_2 = tmp_68;
      real_t a_1_0 = tmp_62;
      real_t a_1_1 = tmp_19 * ( tmp_52 * tmp_52 ) + tmp_27 * ( tmp_54 * tmp_54 ) + tmp_35 * ( tmp_56 * tmp_56 ) +
                     tmp_43 * ( tmp_58 * tmp_58 ) + tmp_51 * ( tmp_60 * tmp_60 );
      real_t a_1_2 = tmp_69;
      real_t a_2_0 = tmp_68;
      real_t a_2_1 = tmp_69;
      real_t a_2_2 = tmp_19 * ( tmp_63 * tmp_63 ) + tmp_27 * ( tmp_64 * tmp_64 ) + tmp_35 * ( tmp_65 * tmp_65 ) +
                     tmp_43 * ( tmp_66 * tmp_66 ) + tmp_51 * ( tmp_67 * tmp_67 );
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

      real_t tmp_0  = 0.063089014491502282;
      real_t tmp_1  = -p_affine_0_0;
      real_t tmp_2  = p_affine_1_0 + tmp_1;
      real_t tmp_3  = -p_affine_0_1;
      real_t tmp_4  = 1.0 / ( tmp_2 * ( p_affine_2_1 + tmp_3 ) - ( p_affine_1_1 + tmp_3 ) * ( p_affine_2_0 + tmp_1 ) );
      real_t tmp_5  = tmp_2 * tmp_4;
      real_t tmp_6  = tmp_4 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_7  = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                               p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_8  = tmp_7 * ( tmp_5 + tmp_6 );
      real_t tmp_9  = 0.025422453185103409 * tmp_8;
      real_t tmp_10 = 0.24928674517091043;
      real_t tmp_11 = 0.058393137863189684 * tmp_8;
      real_t tmp_12 = 0.63650249912139867;
      real_t tmp_13 = 0.041425537809186785 * tmp_8;
      real_t tmp_14 = 0.053145049844816938;
      real_t tmp_15 = 0.041425537809186785 * tmp_8;
      real_t tmp_16 = 0.063089014491502227;
      real_t tmp_17 = 0.025422453185103409 * tmp_8;
      real_t tmp_18 = 0.24928674517091043;
      real_t tmp_19 = 0.058393137863189684 * tmp_8;
      real_t tmp_20 = 0.87382197101699566;
      real_t tmp_21 = 0.025422453185103409 * tmp_8;
      real_t tmp_22 = 0.50142650965817914;
      real_t tmp_23 = 0.058393137863189684 * tmp_8;
      real_t tmp_24 = 0.053145049844816938;
      real_t tmp_25 = 0.041425537809186785 * tmp_8;
      real_t tmp_26 = 0.63650249912139867;
      real_t tmp_27 = 0.041425537809186785 * tmp_8;
      real_t tmp_28 = 0.31035245103378439;
      real_t tmp_29 = 0.041425537809186785 * tmp_8;
      real_t tmp_30 = 0.31035245103378439;
      real_t tmp_31 = 0.041425537809186785 * tmp_8;
      real_t tmp_32 = tmp_6 * tmp_7;
      real_t tmp_33 = 0.025422453185103409 * tmp_32;
      real_t tmp_34 = 0.058393137863189684 * tmp_32;
      real_t tmp_35 = 0.041425537809186785 * tmp_32;
      real_t tmp_36 = 0.041425537809186785 * tmp_32;
      real_t tmp_37 = 0.025422453185103409 * tmp_32;
      real_t tmp_38 = 0.058393137863189684 * tmp_32;
      real_t tmp_39 = 0.025422453185103409 * tmp_32;
      real_t tmp_40 = 0.058393137863189684 * tmp_32;
      real_t tmp_41 = 0.041425537809186785 * tmp_32;
      real_t tmp_42 = 0.041425537809186785 * tmp_32;
      real_t tmp_43 = 0.041425537809186785 * tmp_32;
      real_t tmp_44 = 0.041425537809186785 * tmp_32;
      real_t tmp_45 = tmp_5 * tmp_7;
      real_t tmp_46 = 0.025422453185103409 * tmp_45;
      real_t tmp_47 = 0.058393137863189684 * tmp_45;
      real_t tmp_48 = 0.041425537809186785 * tmp_45;
      real_t tmp_49 = 0.041425537809186785 * tmp_45;
      real_t tmp_50 = 0.025422453185103409 * tmp_45;
      real_t tmp_51 = 0.058393137863189684 * tmp_45;
      real_t tmp_52 = 0.025422453185103409 * tmp_45;
      real_t tmp_53 = 0.058393137863189684 * tmp_45;
      real_t tmp_54 = 0.041425537809186785 * tmp_45;
      real_t tmp_55 = 0.041425537809186785 * tmp_45;
      real_t tmp_56 = 0.041425537809186785 * tmp_45;
      real_t tmp_57 = 0.041425537809186785 * tmp_45;
      real_t a_0_0  = tmp_0 * tmp_9 + tmp_10 * tmp_11 + tmp_12 * tmp_13 + tmp_14 * tmp_15 + tmp_16 * tmp_17 + tmp_18 * tmp_19 +
                     tmp_20 * tmp_21 + tmp_22 * tmp_23 + tmp_24 * tmp_25 + tmp_26 * tmp_27 + tmp_28 * tmp_29 + tmp_30 * tmp_31;
      real_t a_0_1 = 0.24928674517091043 * tmp_11 + 0.31035245103378439 * tmp_13 + 0.63650249912139867 * tmp_15 +
                     0.87382197101699555 * tmp_17 + 0.50142650965817914 * tmp_19 + 0.063089014491502227 * tmp_21 +
                     0.24928674517091043 * tmp_23 + 0.31035245103378439 * tmp_25 + 0.053145049844816938 * tmp_27 +
                     0.63650249912139867 * tmp_29 + 0.053145049844816938 * tmp_31 + 0.063089014491502227 * tmp_9;
      real_t a_0_2 = 0.50142650965817914 * tmp_11 + 0.053145049844816938 * tmp_13 + 0.31035245103378439 * tmp_15 +
                     0.063089014491502227 * tmp_17 + 0.24928674517091043 * tmp_19 + 0.063089014491502227 * tmp_21 +
                     0.24928674517091043 * tmp_23 + 0.63650249912139867 * tmp_25 + 0.31035245103378439 * tmp_27 +
                     0.053145049844816938 * tmp_29 + 0.63650249912139867 * tmp_31 + 0.87382197101699555 * tmp_9;
      real_t a_1_0 = -tmp_0 * tmp_33 - tmp_10 * tmp_34 - tmp_12 * tmp_35 - tmp_14 * tmp_36 - tmp_16 * tmp_37 - tmp_18 * tmp_38 -
                     tmp_20 * tmp_39 - tmp_22 * tmp_40 - tmp_24 * tmp_41 - tmp_26 * tmp_42 - tmp_28 * tmp_43 - tmp_30 * tmp_44;
      real_t a_1_1 = -0.063089014491502227 * tmp_33 - 0.24928674517091043 * tmp_34 - 0.31035245103378439 * tmp_35 -
                     0.63650249912139867 * tmp_36 - 0.87382197101699555 * tmp_37 - 0.50142650965817914 * tmp_38 -
                     0.063089014491502227 * tmp_39 - 0.24928674517091043 * tmp_40 - 0.31035245103378439 * tmp_41 -
                     0.053145049844816938 * tmp_42 - 0.63650249912139867 * tmp_43 - 0.053145049844816938 * tmp_44;
      real_t a_1_2 = -0.87382197101699555 * tmp_33 - 0.50142650965817914 * tmp_34 - 0.053145049844816938 * tmp_35 -
                     0.31035245103378439 * tmp_36 - 0.063089014491502227 * tmp_37 - 0.24928674517091043 * tmp_38 -
                     0.063089014491502227 * tmp_39 - 0.24928674517091043 * tmp_40 - 0.63650249912139867 * tmp_41 -
                     0.31035245103378439 * tmp_42 - 0.053145049844816938 * tmp_43 - 0.63650249912139867 * tmp_44;
      real_t a_2_0 = -tmp_0 * tmp_46 - tmp_10 * tmp_47 - tmp_12 * tmp_48 - tmp_14 * tmp_49 - tmp_16 * tmp_50 - tmp_18 * tmp_51 -
                     tmp_20 * tmp_52 - tmp_22 * tmp_53 - tmp_24 * tmp_54 - tmp_26 * tmp_55 - tmp_28 * tmp_56 - tmp_30 * tmp_57;
      real_t a_2_1 = -0.063089014491502227 * tmp_46 - 0.24928674517091043 * tmp_47 - 0.31035245103378439 * tmp_48 -
                     0.63650249912139867 * tmp_49 - 0.87382197101699555 * tmp_50 - 0.50142650965817914 * tmp_51 -
                     0.063089014491502227 * tmp_52 - 0.24928674517091043 * tmp_53 - 0.31035245103378439 * tmp_54 -
                     0.053145049844816938 * tmp_55 - 0.63650249912139867 * tmp_56 - 0.053145049844816938 * tmp_57;
      real_t a_2_2 = -0.87382197101699555 * tmp_46 - 0.50142650965817914 * tmp_47 - 0.053145049844816938 * tmp_48 -
                     0.31035245103378439 * tmp_49 - 0.063089014491502227 * tmp_50 - 0.24928674517091043 * tmp_51 -
                     0.063089014491502227 * tmp_52 - 0.24928674517091043 * tmp_53 - 0.63650249912139867 * tmp_54 -
                     0.31035245103378439 * tmp_55 - 0.053145049844816938 * tmp_56 - 0.63650249912139867 * tmp_57;
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + tmp_2;
      real_t tmp_7  = tmp_4 * ( 0.046910077030668018 * tmp_5 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = 0.5 * p_affine_10_1 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_19 = 0.11846344252809471 * tmp_18;
      real_t tmp_20 = tmp_4 * ( 0.23076534494715845 * tmp_5 + tmp_6 );
      real_t tmp_21 = tmp_1 * tmp_20;
      real_t tmp_22 = tmp_20 * tmp_9;
      real_t tmp_23 = tmp_4 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_24 = tmp_23 * tmp_3;
      real_t tmp_25 = tmp_15 * tmp_23;
      real_t tmp_26 = -tmp_21 - tmp_22 - tmp_24 - tmp_25 + 1;
      real_t tmp_27 = 0.2393143352496831 * tmp_18;
      real_t tmp_28 = tmp_4 * ( 0.5 * tmp_5 + tmp_6 );
      real_t tmp_29 = tmp_1 * tmp_28;
      real_t tmp_30 = tmp_28 * tmp_9;
      real_t tmp_31 = tmp_4 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_32 = tmp_3 * tmp_31;
      real_t tmp_33 = tmp_15 * tmp_31;
      real_t tmp_34 = -tmp_29 - tmp_30 - tmp_32 - tmp_33 + 1;
      real_t tmp_35 = 0.2844444444444445 * tmp_18;
      real_t tmp_36 = tmp_4 * ( 0.7692346550528415 * tmp_5 + tmp_6 );
      real_t tmp_37 = tmp_1 * tmp_36;
      real_t tmp_38 = tmp_36 * tmp_9;
      real_t tmp_39 = tmp_4 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_40 = tmp_3 * tmp_39;
      real_t tmp_41 = tmp_15 * tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = 0.2393143352496831 * tmp_18;
      real_t tmp_44 = tmp_4 * ( 0.95308992296933193 * tmp_5 + tmp_6 );
      real_t tmp_45 = tmp_1 * tmp_44;
      real_t tmp_46 = tmp_44 * tmp_9;
      real_t tmp_47 = tmp_4 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_48 = tmp_3 * tmp_47;
      real_t tmp_49 = tmp_15 * tmp_47;
      real_t tmp_50 = -tmp_45 - tmp_46 - tmp_48 - tmp_49 + 1;
      real_t tmp_51 = 0.11846344252809471 * tmp_18;
      real_t tmp_52 = tmp_10 + tmp_14;
      real_t tmp_53 = tmp_17 * tmp_19;
      real_t tmp_54 = tmp_22 + tmp_24;
      real_t tmp_55 = tmp_26 * tmp_27;
      real_t tmp_56 = tmp_30 + tmp_32;
      real_t tmp_57 = tmp_34 * tmp_35;
      real_t tmp_58 = tmp_38 + tmp_40;
      real_t tmp_59 = tmp_42 * tmp_43;
      real_t tmp_60 = tmp_46 + tmp_48;
      real_t tmp_61 = tmp_50 * tmp_51;
      real_t tmp_62 = tmp_52 * tmp_53 + tmp_54 * tmp_55 + tmp_56 * tmp_57 + tmp_58 * tmp_59 + tmp_60 * tmp_61;
      real_t tmp_63 = tmp_16 + tmp_8;
      real_t tmp_64 = tmp_21 + tmp_25;
      real_t tmp_65 = tmp_29 + tmp_33;
      real_t tmp_66 = tmp_37 + tmp_41;
      real_t tmp_67 = tmp_45 + tmp_49;
      real_t tmp_68 = tmp_53 * tmp_63 + tmp_55 * tmp_64 + tmp_57 * tmp_65 + tmp_59 * tmp_66 + tmp_61 * tmp_67;
      real_t tmp_69 = tmp_19 * tmp_52 * tmp_63 + tmp_27 * tmp_54 * tmp_64 + tmp_35 * tmp_56 * tmp_65 + tmp_43 * tmp_58 * tmp_66 +
                      tmp_51 * tmp_60 * tmp_67;
      real_t a_0_0 = ( tmp_17 * tmp_17 ) * tmp_19 + ( tmp_26 * tmp_26 ) * tmp_27 + ( tmp_34 * tmp_34 ) * tmp_35 +
                     ( tmp_42 * tmp_42 ) * tmp_43 + ( tmp_50 * tmp_50 ) * tmp_51;
      real_t a_0_1 = tmp_62;
      real_t a_0_2 = tmp_68;
      real_t a_1_0 = tmp_62;
      real_t a_1_1 = tmp_19 * ( tmp_52 * tmp_52 ) + tmp_27 * ( tmp_54 * tmp_54 ) + tmp_35 * ( tmp_56 * tmp_56 ) +
                     tmp_43 * ( tmp_58 * tmp_58 ) + tmp_51 * ( tmp_60 * tmp_60 );
      real_t a_1_2 = tmp_69;
      real_t a_2_0 = tmp_68;
      real_t a_2_1 = tmp_69;
      real_t a_2_2 = tmp_19 * ( tmp_63 * tmp_63 ) + tmp_27 * ( tmp_64 * tmp_64 ) + tmp_35 * ( tmp_65 * tmp_65 ) +
                     tmp_43 * ( tmp_66 * tmp_66 ) + tmp_51 * ( tmp_67 * tmp_67 );
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

      real_t tmp_0   = -p_affine_3_0;
      real_t tmp_1   = p_affine_4_0 + tmp_0;
      real_t tmp_2   = -p_affine_3_1;
      real_t tmp_3   = p_affine_5_1 + tmp_2;
      real_t tmp_4   = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_4_1 + tmp_2 ) * ( p_affine_5_0 + tmp_0 ) );
      real_t tmp_5   = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6   = p_affine_6_1 + 0.046910077030668018 * tmp_5;
      real_t tmp_7   = tmp_4 * ( tmp_2 + tmp_6 );
      real_t tmp_8   = tmp_1 * tmp_7;
      real_t tmp_9   = p_affine_3_0 - p_affine_5_0;
      real_t tmp_10  = tmp_7 * tmp_9;
      real_t tmp_11  = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12  = p_affine_6_0 + 0.046910077030668018 * tmp_11;
      real_t tmp_13  = tmp_4 * ( tmp_0 + tmp_12 );
      real_t tmp_14  = tmp_13 * tmp_3;
      real_t tmp_15  = p_affine_3_1 - p_affine_4_1;
      real_t tmp_16  = tmp_13 * tmp_15;
      real_t tmp_17  = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18  = -p_affine_0_0;
      real_t tmp_19  = p_affine_1_0 + tmp_18;
      real_t tmp_20  = -p_affine_0_1;
      real_t tmp_21  = p_affine_2_1 + tmp_20;
      real_t tmp_22  = 1.0 / ( tmp_19 * tmp_21 - ( p_affine_1_1 + tmp_20 ) * ( p_affine_2_0 + tmp_18 ) );
      real_t tmp_23  = tmp_22 * ( tmp_20 + tmp_6 );
      real_t tmp_24  = tmp_19 * tmp_23;
      real_t tmp_25  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_26  = tmp_23 * tmp_25;
      real_t tmp_27  = tmp_22 * ( tmp_12 + tmp_18 );
      real_t tmp_28  = tmp_21 * tmp_27;
      real_t tmp_29  = p_affine_0_1 - p_affine_1_1;
      real_t tmp_30  = tmp_27 * tmp_29;
      real_t tmp_31  = 0.5 * p_affine_10_1 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_32  = 0.11846344252809471 * tmp_31;
      real_t tmp_33  = tmp_32 * ( -tmp_24 - tmp_26 - tmp_28 - tmp_30 + 1 );
      real_t tmp_34  = p_affine_6_1 + 0.23076534494715845 * tmp_5;
      real_t tmp_35  = tmp_4 * ( tmp_2 + tmp_34 );
      real_t tmp_36  = tmp_1 * tmp_35;
      real_t tmp_37  = tmp_35 * tmp_9;
      real_t tmp_38  = p_affine_6_0 + 0.23076534494715845 * tmp_11;
      real_t tmp_39  = tmp_4 * ( tmp_0 + tmp_38 );
      real_t tmp_40  = tmp_3 * tmp_39;
      real_t tmp_41  = tmp_15 * tmp_39;
      real_t tmp_42  = -tmp_36 - tmp_37 - tmp_40 - tmp_41 + 1;
      real_t tmp_43  = tmp_22 * ( tmp_20 + tmp_34 );
      real_t tmp_44  = tmp_19 * tmp_43;
      real_t tmp_45  = tmp_25 * tmp_43;
      real_t tmp_46  = tmp_22 * ( tmp_18 + tmp_38 );
      real_t tmp_47  = tmp_21 * tmp_46;
      real_t tmp_48  = tmp_29 * tmp_46;
      real_t tmp_49  = 0.2393143352496831 * tmp_31;
      real_t tmp_50  = tmp_49 * ( -tmp_44 - tmp_45 - tmp_47 - tmp_48 + 1 );
      real_t tmp_51  = p_affine_6_1 + 0.5 * tmp_5;
      real_t tmp_52  = tmp_4 * ( tmp_2 + tmp_51 );
      real_t tmp_53  = tmp_1 * tmp_52;
      real_t tmp_54  = tmp_52 * tmp_9;
      real_t tmp_55  = p_affine_6_0 + 0.5 * tmp_11;
      real_t tmp_56  = tmp_4 * ( tmp_0 + tmp_55 );
      real_t tmp_57  = tmp_3 * tmp_56;
      real_t tmp_58  = tmp_15 * tmp_56;
      real_t tmp_59  = -tmp_53 - tmp_54 - tmp_57 - tmp_58 + 1;
      real_t tmp_60  = tmp_22 * ( tmp_20 + tmp_51 );
      real_t tmp_61  = tmp_19 * tmp_60;
      real_t tmp_62  = tmp_25 * tmp_60;
      real_t tmp_63  = tmp_22 * ( tmp_18 + tmp_55 );
      real_t tmp_64  = tmp_21 * tmp_63;
      real_t tmp_65  = tmp_29 * tmp_63;
      real_t tmp_66  = 0.2844444444444445 * tmp_31;
      real_t tmp_67  = tmp_66 * ( -tmp_61 - tmp_62 - tmp_64 - tmp_65 + 1 );
      real_t tmp_68  = p_affine_6_1 + 0.7692346550528415 * tmp_5;
      real_t tmp_69  = tmp_4 * ( tmp_2 + tmp_68 );
      real_t tmp_70  = tmp_1 * tmp_69;
      real_t tmp_71  = tmp_69 * tmp_9;
      real_t tmp_72  = p_affine_6_0 + 0.7692346550528415 * tmp_11;
      real_t tmp_73  = tmp_4 * ( tmp_0 + tmp_72 );
      real_t tmp_74  = tmp_3 * tmp_73;
      real_t tmp_75  = tmp_15 * tmp_73;
      real_t tmp_76  = -tmp_70 - tmp_71 - tmp_74 - tmp_75 + 1;
      real_t tmp_77  = tmp_22 * ( tmp_20 + tmp_68 );
      real_t tmp_78  = tmp_19 * tmp_77;
      real_t tmp_79  = tmp_25 * tmp_77;
      real_t tmp_80  = tmp_22 * ( tmp_18 + tmp_72 );
      real_t tmp_81  = tmp_21 * tmp_80;
      real_t tmp_82  = tmp_29 * tmp_80;
      real_t tmp_83  = 0.2393143352496831 * tmp_31;
      real_t tmp_84  = tmp_83 * ( -tmp_78 - tmp_79 - tmp_81 - tmp_82 + 1 );
      real_t tmp_85  = p_affine_6_1 + 0.95308992296933193 * tmp_5;
      real_t tmp_86  = tmp_4 * ( tmp_2 + tmp_85 );
      real_t tmp_87  = tmp_1 * tmp_86;
      real_t tmp_88  = tmp_86 * tmp_9;
      real_t tmp_89  = p_affine_6_0 + 0.95308992296933193 * tmp_11;
      real_t tmp_90  = tmp_4 * ( tmp_0 + tmp_89 );
      real_t tmp_91  = tmp_3 * tmp_90;
      real_t tmp_92  = tmp_15 * tmp_90;
      real_t tmp_93  = -tmp_87 - tmp_88 - tmp_91 - tmp_92 + 1;
      real_t tmp_94  = tmp_22 * ( tmp_20 + tmp_85 );
      real_t tmp_95  = tmp_19 * tmp_94;
      real_t tmp_96  = tmp_25 * tmp_94;
      real_t tmp_97  = tmp_22 * ( tmp_18 + tmp_89 );
      real_t tmp_98  = tmp_21 * tmp_97;
      real_t tmp_99  = tmp_29 * tmp_97;
      real_t tmp_100 = 0.11846344252809471 * tmp_31;
      real_t tmp_101 = tmp_100 * ( -tmp_95 - tmp_96 - tmp_98 - tmp_99 + 1 );
      real_t tmp_102 = tmp_10 + tmp_14;
      real_t tmp_103 = tmp_37 + tmp_40;
      real_t tmp_104 = tmp_54 + tmp_57;
      real_t tmp_105 = tmp_71 + tmp_74;
      real_t tmp_106 = tmp_88 + tmp_91;
      real_t tmp_107 = tmp_16 + tmp_8;
      real_t tmp_108 = tmp_36 + tmp_41;
      real_t tmp_109 = tmp_53 + tmp_58;
      real_t tmp_110 = tmp_70 + tmp_75;
      real_t tmp_111 = tmp_87 + tmp_92;
      real_t tmp_112 = tmp_32 * ( tmp_26 + tmp_28 );
      real_t tmp_113 = tmp_49 * ( tmp_45 + tmp_47 );
      real_t tmp_114 = tmp_66 * ( tmp_62 + tmp_64 );
      real_t tmp_115 = tmp_83 * ( tmp_79 + tmp_81 );
      real_t tmp_116 = tmp_100 * ( tmp_96 + tmp_98 );
      real_t tmp_117 = tmp_32 * ( tmp_24 + tmp_30 );
      real_t tmp_118 = tmp_49 * ( tmp_44 + tmp_48 );
      real_t tmp_119 = tmp_66 * ( tmp_61 + tmp_65 );
      real_t tmp_120 = tmp_83 * ( tmp_78 + tmp_82 );
      real_t tmp_121 = tmp_100 * ( tmp_95 + tmp_99 );
      real_t a_0_0   = tmp_101 * tmp_93 + tmp_17 * tmp_33 + tmp_42 * tmp_50 + tmp_59 * tmp_67 + tmp_76 * tmp_84;
      real_t a_0_1   = tmp_101 * tmp_106 + tmp_102 * tmp_33 + tmp_103 * tmp_50 + tmp_104 * tmp_67 + tmp_105 * tmp_84;
      real_t a_0_2   = tmp_101 * tmp_111 + tmp_107 * tmp_33 + tmp_108 * tmp_50 + tmp_109 * tmp_67 + tmp_110 * tmp_84;
      real_t a_1_0   = tmp_112 * tmp_17 + tmp_113 * tmp_42 + tmp_114 * tmp_59 + tmp_115 * tmp_76 + tmp_116 * tmp_93;
      real_t a_1_1   = tmp_102 * tmp_112 + tmp_103 * tmp_113 + tmp_104 * tmp_114 + tmp_105 * tmp_115 + tmp_106 * tmp_116;
      real_t a_1_2   = tmp_107 * tmp_112 + tmp_108 * tmp_113 + tmp_109 * tmp_114 + tmp_110 * tmp_115 + tmp_111 * tmp_116;
      real_t a_2_0   = tmp_117 * tmp_17 + tmp_118 * tmp_42 + tmp_119 * tmp_59 + tmp_120 * tmp_76 + tmp_121 * tmp_93;
      real_t a_2_1   = tmp_102 * tmp_117 + tmp_103 * tmp_118 + tmp_104 * tmp_119 + tmp_105 * tmp_120 + tmp_106 * tmp_121;
      real_t a_2_2   = tmp_107 * tmp_117 + tmp_108 * tmp_118 + tmp_109 * tmp_119 + tmp_110 * tmp_120 + tmp_111 * tmp_121;
      elMat( 0, 0 )  = a_0_0;
      elMat( 0, 1 )  = a_0_1;
      elMat( 0, 2 )  = a_0_2;
      elMat( 1, 0 )  = a_1_0;
      elMat( 1, 1 )  = a_1_1;
      elMat( 1, 2 )  = a_1_2;
      elMat( 2, 0 )  = a_2_0;
      elMat( 2, 1 )  = a_2_1;
      elMat( 2, 2 )  = a_2_2;
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + tmp_2;
      real_t tmp_7  = tmp_4 * ( 0.046910077030668018 * tmp_5 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1;
      real_t tmp_18 = p_affine_10_1 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_19 = 0.11846344252809471 * tmp_18;
      real_t tmp_20 = tmp_4 * ( 0.23076534494715845 * tmp_5 + tmp_6 );
      real_t tmp_21 = tmp_1 * tmp_20;
      real_t tmp_22 = tmp_20 * tmp_9;
      real_t tmp_23 = tmp_4 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_24 = tmp_23 * tmp_3;
      real_t tmp_25 = tmp_15 * tmp_23;
      real_t tmp_26 = -tmp_21 - tmp_22 - tmp_24 - tmp_25 + 1;
      real_t tmp_27 = 0.2393143352496831 * tmp_18;
      real_t tmp_28 = tmp_4 * ( 0.5 * tmp_5 + tmp_6 );
      real_t tmp_29 = tmp_1 * tmp_28;
      real_t tmp_30 = tmp_28 * tmp_9;
      real_t tmp_31 = tmp_4 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_32 = tmp_3 * tmp_31;
      real_t tmp_33 = tmp_15 * tmp_31;
      real_t tmp_34 = -tmp_29 - tmp_30 - tmp_32 - tmp_33 + 1;
      real_t tmp_35 = 0.2844444444444445 * tmp_18;
      real_t tmp_36 = tmp_4 * ( 0.7692346550528415 * tmp_5 + tmp_6 );
      real_t tmp_37 = tmp_1 * tmp_36;
      real_t tmp_38 = tmp_36 * tmp_9;
      real_t tmp_39 = tmp_4 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_40 = tmp_3 * tmp_39;
      real_t tmp_41 = tmp_15 * tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = 0.2393143352496831 * tmp_18;
      real_t tmp_44 = tmp_4 * ( 0.95308992296933193 * tmp_5 + tmp_6 );
      real_t tmp_45 = tmp_1 * tmp_44;
      real_t tmp_46 = tmp_44 * tmp_9;
      real_t tmp_47 = tmp_4 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_48 = tmp_3 * tmp_47;
      real_t tmp_49 = tmp_15 * tmp_47;
      real_t tmp_50 = -tmp_45 - tmp_46 - tmp_48 - tmp_49 + 1;
      real_t tmp_51 = 0.11846344252809471 * tmp_18;
      real_t tmp_52 = tmp_10 + tmp_14;
      real_t tmp_53 = tmp_17 * tmp_19;
      real_t tmp_54 = tmp_22 + tmp_24;
      real_t tmp_55 = tmp_26 * tmp_27;
      real_t tmp_56 = tmp_30 + tmp_32;
      real_t tmp_57 = tmp_34 * tmp_35;
      real_t tmp_58 = tmp_38 + tmp_40;
      real_t tmp_59 = tmp_42 * tmp_43;
      real_t tmp_60 = tmp_46 + tmp_48;
      real_t tmp_61 = tmp_50 * tmp_51;
      real_t tmp_62 = tmp_52 * tmp_53 + tmp_54 * tmp_55 + tmp_56 * tmp_57 + tmp_58 * tmp_59 + tmp_60 * tmp_61;
      real_t tmp_63 = tmp_16 + tmp_8;
      real_t tmp_64 = tmp_21 + tmp_25;
      real_t tmp_65 = tmp_29 + tmp_33;
      real_t tmp_66 = tmp_37 + tmp_41;
      real_t tmp_67 = tmp_45 + tmp_49;
      real_t tmp_68 = tmp_53 * tmp_63 + tmp_55 * tmp_64 + tmp_57 * tmp_65 + tmp_59 * tmp_66 + tmp_61 * tmp_67;
      real_t tmp_69 = tmp_19 * tmp_52 * tmp_63 + tmp_27 * tmp_54 * tmp_64 + tmp_35 * tmp_56 * tmp_65 + tmp_43 * tmp_58 * tmp_66 +
                      tmp_51 * tmp_60 * tmp_67;
      real_t a_0_0 = ( tmp_17 * tmp_17 ) * tmp_19 + ( tmp_26 * tmp_26 ) * tmp_27 + ( tmp_34 * tmp_34 ) * tmp_35 +
                     ( tmp_42 * tmp_42 ) * tmp_43 + ( tmp_50 * tmp_50 ) * tmp_51;
      real_t a_0_1 = tmp_62;
      real_t a_0_2 = tmp_68;
      real_t a_1_0 = tmp_62;
      real_t a_1_1 = tmp_19 * ( tmp_52 * tmp_52 ) + tmp_27 * ( tmp_54 * tmp_54 ) + tmp_35 * ( tmp_56 * tmp_56 ) +
                     tmp_43 * ( tmp_58 * tmp_58 ) + tmp_51 * ( tmp_60 * tmp_60 );
      real_t a_1_2 = tmp_69;
      real_t a_2_0 = tmp_68;
      real_t a_2_1 = tmp_69;
      real_t a_2_2 = tmp_19 * ( tmp_63 * tmp_63 ) + tmp_27 * ( tmp_64 * tmp_64 ) + tmp_35 * ( tmp_65 * tmp_65 ) +
                     tmp_43 * ( tmp_66 * tmp_66 ) + tmp_51 * ( tmp_67 * tmp_67 );
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


class DGDivFormP0P1_0 : public hyteg::dg::DGForm2D
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

      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_2_1 + tmp_0;
      real_t tmp_2  = -p_affine_0_0;
      real_t tmp_3  = 1.0 / ( tmp_1 * ( p_affine_1_0 + tmp_2 ) - ( p_affine_1_1 + tmp_0 ) * ( p_affine_2_0 + tmp_2 ) );
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = tmp_3 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_6  = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                               p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_7  = tmp_6 * ( tmp_4 + tmp_5 );
      real_t tmp_8  = tmp_4 * tmp_6;
      real_t tmp_9  = tmp_5 * tmp_6;
      real_t a_0_0  = 0.5 * tmp_7;
      real_t a_0_1  = -0.5 * tmp_8;
      real_t a_0_2  = -0.5 * tmp_9;
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
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + tmp_2;
      real_t tmp_7  = tmp_4 * ( 0.046910077030668018 * tmp_5 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = 0.5 * p_affine_10_0 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_18 = 0.11846344252809471 * tmp_17;
      real_t tmp_19 = tmp_4 * ( 0.23076534494715845 * tmp_5 + tmp_6 );
      real_t tmp_20 = tmp_1 * tmp_19;
      real_t tmp_21 = tmp_19 * tmp_9;
      real_t tmp_22 = tmp_4 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_23 = tmp_22 * tmp_3;
      real_t tmp_24 = tmp_15 * tmp_22;
      real_t tmp_25 = 0.2393143352496831 * tmp_17;
      real_t tmp_26 = tmp_4 * ( 0.5 * tmp_5 + tmp_6 );
      real_t tmp_27 = tmp_1 * tmp_26;
      real_t tmp_28 = tmp_26 * tmp_9;
      real_t tmp_29 = tmp_4 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_30 = tmp_29 * tmp_3;
      real_t tmp_31 = tmp_15 * tmp_29;
      real_t tmp_32 = 0.2844444444444445 * tmp_17;
      real_t tmp_33 = tmp_4 * ( 0.7692346550528415 * tmp_5 + tmp_6 );
      real_t tmp_34 = tmp_1 * tmp_33;
      real_t tmp_35 = tmp_33 * tmp_9;
      real_t tmp_36 = tmp_4 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_37 = tmp_3 * tmp_36;
      real_t tmp_38 = tmp_15 * tmp_36;
      real_t tmp_39 = 0.2393143352496831 * tmp_17;
      real_t tmp_40 = tmp_4 * ( 0.95308992296933193 * tmp_5 + tmp_6 );
      real_t tmp_41 = tmp_1 * tmp_40;
      real_t tmp_42 = tmp_40 * tmp_9;
      real_t tmp_43 = tmp_4 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_44 = tmp_3 * tmp_43;
      real_t tmp_45 = tmp_15 * tmp_43;
      real_t tmp_46 = 0.11846344252809471 * tmp_17;
      real_t a_0_0  = tmp_18 * ( -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1 ) + tmp_25 * ( -tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1 ) +
                     tmp_32 * ( -tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1 ) + tmp_39 * ( -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1 ) +
                     tmp_46 * ( -tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1 );
      real_t a_0_1 = tmp_18 * ( tmp_10 + tmp_14 ) + tmp_25 * ( tmp_21 + tmp_23 ) + tmp_32 * ( tmp_28 + tmp_30 ) +
                     tmp_39 * ( tmp_35 + tmp_37 ) + tmp_46 * ( tmp_42 + tmp_44 );
      real_t a_0_2 = tmp_18 * ( tmp_16 + tmp_8 ) + tmp_25 * ( tmp_20 + tmp_24 ) + tmp_32 * ( tmp_27 + tmp_31 ) +
                     tmp_39 * ( tmp_34 + tmp_38 ) + tmp_46 * ( tmp_41 + tmp_45 );
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

      real_t tmp_0  = -p_affine_3_0;
      real_t tmp_1  = p_affine_4_0 + tmp_0;
      real_t tmp_2  = -p_affine_3_1;
      real_t tmp_3  = p_affine_5_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_4_1 + tmp_2 ) * ( p_affine_5_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + tmp_2;
      real_t tmp_7  = tmp_4 * ( 0.046910077030668018 * tmp_5 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_3_0 - p_affine_5_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = 0.5 * p_affine_10_0 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_18 = 0.11846344252809471 * tmp_17;
      real_t tmp_19 = tmp_4 * ( 0.23076534494715845 * tmp_5 + tmp_6 );
      real_t tmp_20 = tmp_1 * tmp_19;
      real_t tmp_21 = tmp_19 * tmp_9;
      real_t tmp_22 = tmp_4 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_23 = tmp_22 * tmp_3;
      real_t tmp_24 = tmp_15 * tmp_22;
      real_t tmp_25 = 0.2393143352496831 * tmp_17;
      real_t tmp_26 = tmp_4 * ( 0.5 * tmp_5 + tmp_6 );
      real_t tmp_27 = tmp_1 * tmp_26;
      real_t tmp_28 = tmp_26 * tmp_9;
      real_t tmp_29 = tmp_4 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_30 = tmp_29 * tmp_3;
      real_t tmp_31 = tmp_15 * tmp_29;
      real_t tmp_32 = 0.2844444444444445 * tmp_17;
      real_t tmp_33 = tmp_4 * ( 0.7692346550528415 * tmp_5 + tmp_6 );
      real_t tmp_34 = tmp_1 * tmp_33;
      real_t tmp_35 = tmp_33 * tmp_9;
      real_t tmp_36 = tmp_4 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_37 = tmp_3 * tmp_36;
      real_t tmp_38 = tmp_15 * tmp_36;
      real_t tmp_39 = 0.2393143352496831 * tmp_17;
      real_t tmp_40 = tmp_4 * ( 0.95308992296933193 * tmp_5 + tmp_6 );
      real_t tmp_41 = tmp_1 * tmp_40;
      real_t tmp_42 = tmp_40 * tmp_9;
      real_t tmp_43 = tmp_4 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_44 = tmp_3 * tmp_43;
      real_t tmp_45 = tmp_15 * tmp_43;
      real_t tmp_46 = 0.11846344252809471 * tmp_17;
      real_t a_0_0  = -tmp_18 * ( -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1 ) - tmp_25 * ( -tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1 ) -
                     tmp_32 * ( -tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1 ) - tmp_39 * ( -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1 ) -
                     tmp_46 * ( -tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1 );
      real_t a_0_1 = -tmp_18 * ( tmp_10 + tmp_14 ) - tmp_25 * ( tmp_21 + tmp_23 ) - tmp_32 * ( tmp_28 + tmp_30 ) -
                     tmp_39 * ( tmp_35 + tmp_37 ) - tmp_46 * ( tmp_42 + tmp_44 );
      real_t a_0_2 = -tmp_18 * ( tmp_16 + tmp_8 ) - tmp_25 * ( tmp_20 + tmp_24 ) - tmp_32 * ( tmp_27 + tmp_31 ) -
                     tmp_39 * ( tmp_34 + tmp_38 ) - tmp_46 * ( tmp_41 + tmp_45 );
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

      real_t a_0_0  = 0;
      real_t a_0_1  = 0;
      real_t a_0_2  = 0;
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

class DGDivFormP0P1_1 : public hyteg::dg::DGForm2D
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
      real_t tmp_3  = 1.0 / ( tmp_1 * ( p_affine_2_1 + tmp_2 ) - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = tmp_3 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_6  = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                               p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_7  = tmp_6 * ( tmp_4 + tmp_5 );
      real_t tmp_8  = tmp_5 * tmp_6;
      real_t tmp_9  = tmp_4 * tmp_6;
      real_t a_0_0  = 0.5 * tmp_7;
      real_t a_0_1  = -0.5 * tmp_8;
      real_t a_0_2  = -0.5 * tmp_9;
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
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + tmp_2;
      real_t tmp_7  = tmp_4 * ( 0.046910077030668018 * tmp_5 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = 0.5 * p_affine_10_1 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_18 = 0.11846344252809471 * tmp_17;
      real_t tmp_19 = tmp_4 * ( 0.23076534494715845 * tmp_5 + tmp_6 );
      real_t tmp_20 = tmp_1 * tmp_19;
      real_t tmp_21 = tmp_19 * tmp_9;
      real_t tmp_22 = tmp_4 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_23 = tmp_22 * tmp_3;
      real_t tmp_24 = tmp_15 * tmp_22;
      real_t tmp_25 = 0.2393143352496831 * tmp_17;
      real_t tmp_26 = tmp_4 * ( 0.5 * tmp_5 + tmp_6 );
      real_t tmp_27 = tmp_1 * tmp_26;
      real_t tmp_28 = tmp_26 * tmp_9;
      real_t tmp_29 = tmp_4 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_30 = tmp_29 * tmp_3;
      real_t tmp_31 = tmp_15 * tmp_29;
      real_t tmp_32 = 0.2844444444444445 * tmp_17;
      real_t tmp_33 = tmp_4 * ( 0.7692346550528415 * tmp_5 + tmp_6 );
      real_t tmp_34 = tmp_1 * tmp_33;
      real_t tmp_35 = tmp_33 * tmp_9;
      real_t tmp_36 = tmp_4 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_37 = tmp_3 * tmp_36;
      real_t tmp_38 = tmp_15 * tmp_36;
      real_t tmp_39 = 0.2393143352496831 * tmp_17;
      real_t tmp_40 = tmp_4 * ( 0.95308992296933193 * tmp_5 + tmp_6 );
      real_t tmp_41 = tmp_1 * tmp_40;
      real_t tmp_42 = tmp_40 * tmp_9;
      real_t tmp_43 = tmp_4 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_44 = tmp_3 * tmp_43;
      real_t tmp_45 = tmp_15 * tmp_43;
      real_t tmp_46 = 0.11846344252809471 * tmp_17;
      real_t a_0_0  = tmp_18 * ( -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1 ) + tmp_25 * ( -tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1 ) +
                     tmp_32 * ( -tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1 ) + tmp_39 * ( -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1 ) +
                     tmp_46 * ( -tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1 );
      real_t a_0_1 = tmp_18 * ( tmp_10 + tmp_14 ) + tmp_25 * ( tmp_21 + tmp_23 ) + tmp_32 * ( tmp_28 + tmp_30 ) +
                     tmp_39 * ( tmp_35 + tmp_37 ) + tmp_46 * ( tmp_42 + tmp_44 );
      real_t a_0_2 = tmp_18 * ( tmp_16 + tmp_8 ) + tmp_25 * ( tmp_20 + tmp_24 ) + tmp_32 * ( tmp_27 + tmp_31 ) +
                     tmp_39 * ( tmp_34 + tmp_38 ) + tmp_46 * ( tmp_41 + tmp_45 );
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

      real_t tmp_0  = -p_affine_3_0;
      real_t tmp_1  = p_affine_4_0 + tmp_0;
      real_t tmp_2  = -p_affine_3_1;
      real_t tmp_3  = p_affine_5_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_4_1 + tmp_2 ) * ( p_affine_5_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + tmp_2;
      real_t tmp_7  = tmp_4 * ( 0.046910077030668018 * tmp_5 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_3_0 - p_affine_5_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = 0.5 * p_affine_10_1 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_18 = 0.11846344252809471 * tmp_17;
      real_t tmp_19 = tmp_4 * ( 0.23076534494715845 * tmp_5 + tmp_6 );
      real_t tmp_20 = tmp_1 * tmp_19;
      real_t tmp_21 = tmp_19 * tmp_9;
      real_t tmp_22 = tmp_4 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_23 = tmp_22 * tmp_3;
      real_t tmp_24 = tmp_15 * tmp_22;
      real_t tmp_25 = 0.2393143352496831 * tmp_17;
      real_t tmp_26 = tmp_4 * ( 0.5 * tmp_5 + tmp_6 );
      real_t tmp_27 = tmp_1 * tmp_26;
      real_t tmp_28 = tmp_26 * tmp_9;
      real_t tmp_29 = tmp_4 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_30 = tmp_29 * tmp_3;
      real_t tmp_31 = tmp_15 * tmp_29;
      real_t tmp_32 = 0.2844444444444445 * tmp_17;
      real_t tmp_33 = tmp_4 * ( 0.7692346550528415 * tmp_5 + tmp_6 );
      real_t tmp_34 = tmp_1 * tmp_33;
      real_t tmp_35 = tmp_33 * tmp_9;
      real_t tmp_36 = tmp_4 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_37 = tmp_3 * tmp_36;
      real_t tmp_38 = tmp_15 * tmp_36;
      real_t tmp_39 = 0.2393143352496831 * tmp_17;
      real_t tmp_40 = tmp_4 * ( 0.95308992296933193 * tmp_5 + tmp_6 );
      real_t tmp_41 = tmp_1 * tmp_40;
      real_t tmp_42 = tmp_40 * tmp_9;
      real_t tmp_43 = tmp_4 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_44 = tmp_3 * tmp_43;
      real_t tmp_45 = tmp_15 * tmp_43;
      real_t tmp_46 = 0.11846344252809471 * tmp_17;
      real_t a_0_0  = -tmp_18 * ( -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1 ) - tmp_25 * ( -tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1 ) -
                     tmp_32 * ( -tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1 ) - tmp_39 * ( -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1 ) -
                     tmp_46 * ( -tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1 );
      real_t a_0_1 = -tmp_18 * ( tmp_10 + tmp_14 ) - tmp_25 * ( tmp_21 + tmp_23 ) - tmp_32 * ( tmp_28 + tmp_30 ) -
                     tmp_39 * ( tmp_35 + tmp_37 ) - tmp_46 * ( tmp_42 + tmp_44 );
      real_t a_0_2 = -tmp_18 * ( tmp_16 + tmp_8 ) - tmp_25 * ( tmp_20 + tmp_24 ) - tmp_32 * ( tmp_27 + tmp_31 ) -
                     tmp_39 * ( tmp_34 + tmp_38 ) - tmp_46 * ( tmp_41 + tmp_45 );
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

      real_t a_0_0  = 0;
      real_t a_0_1  = 0;
      real_t a_0_2  = 0;
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


class DGDivtFormP1P0_0 : public hyteg::dg::DGForm2D
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

      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_2_1 + tmp_0;
      real_t tmp_2  = -p_affine_0_0;
      real_t tmp_3  = 1.0 / ( tmp_1 * ( p_affine_1_0 + tmp_2 ) - ( p_affine_1_1 + tmp_0 ) * ( p_affine_2_0 + tmp_2 ) );
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = tmp_3 * ( p_affine_0_1 - p_affine_1_1 );
      real_t tmp_6  = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                               p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_7  = tmp_6 * ( tmp_4 + tmp_5 );
      real_t tmp_8  = tmp_4 * tmp_6;
      real_t tmp_9  = tmp_5 * tmp_6;
      real_t a_0_0  = 0.5 * tmp_7;
      real_t a_1_0  = -0.5 * tmp_8;
      real_t a_2_0  = -0.5 * tmp_9;
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
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + tmp_2;
      real_t tmp_7  = tmp_4 * ( 0.046910077030668018 * tmp_5 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = 0.5 * p_affine_10_0 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_18 = 0.11846344252809471 * tmp_17;
      real_t tmp_19 = tmp_4 * ( 0.23076534494715845 * tmp_5 + tmp_6 );
      real_t tmp_20 = tmp_1 * tmp_19;
      real_t tmp_21 = tmp_19 * tmp_9;
      real_t tmp_22 = tmp_4 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_23 = tmp_22 * tmp_3;
      real_t tmp_24 = tmp_15 * tmp_22;
      real_t tmp_25 = 0.2393143352496831 * tmp_17;
      real_t tmp_26 = tmp_4 * ( 0.5 * tmp_5 + tmp_6 );
      real_t tmp_27 = tmp_1 * tmp_26;
      real_t tmp_28 = tmp_26 * tmp_9;
      real_t tmp_29 = tmp_4 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_30 = tmp_29 * tmp_3;
      real_t tmp_31 = tmp_15 * tmp_29;
      real_t tmp_32 = 0.2844444444444445 * tmp_17;
      real_t tmp_33 = tmp_4 * ( 0.7692346550528415 * tmp_5 + tmp_6 );
      real_t tmp_34 = tmp_1 * tmp_33;
      real_t tmp_35 = tmp_33 * tmp_9;
      real_t tmp_36 = tmp_4 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_37 = tmp_3 * tmp_36;
      real_t tmp_38 = tmp_15 * tmp_36;
      real_t tmp_39 = 0.2393143352496831 * tmp_17;
      real_t tmp_40 = tmp_4 * ( 0.95308992296933193 * tmp_5 + tmp_6 );
      real_t tmp_41 = tmp_1 * tmp_40;
      real_t tmp_42 = tmp_40 * tmp_9;
      real_t tmp_43 = tmp_4 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_44 = tmp_3 * tmp_43;
      real_t tmp_45 = tmp_15 * tmp_43;
      real_t tmp_46 = 0.11846344252809471 * tmp_17;
      real_t a_0_0  = tmp_18 * ( -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1 ) + tmp_25 * ( -tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1 ) +
                     tmp_32 * ( -tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1 ) + tmp_39 * ( -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1 ) +
                     tmp_46 * ( -tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1 );
      real_t a_1_0 = tmp_18 * ( tmp_10 + tmp_14 ) + tmp_25 * ( tmp_21 + tmp_23 ) + tmp_32 * ( tmp_28 + tmp_30 ) +
                     tmp_39 * ( tmp_35 + tmp_37 ) + tmp_46 * ( tmp_42 + tmp_44 );
      real_t a_2_0 = tmp_18 * ( tmp_16 + tmp_8 ) + tmp_25 * ( tmp_20 + tmp_24 ) + tmp_32 * ( tmp_27 + tmp_31 ) +
                     tmp_39 * ( tmp_34 + tmp_38 ) + tmp_46 * ( tmp_41 + tmp_45 );
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + tmp_2;
      real_t tmp_7  = tmp_4 * ( 0.046910077030668018 * tmp_5 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = 0.5 * p_affine_10_0 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_18 = 0.11846344252809471 * tmp_17;
      real_t tmp_19 = tmp_4 * ( 0.23076534494715845 * tmp_5 + tmp_6 );
      real_t tmp_20 = tmp_1 * tmp_19;
      real_t tmp_21 = tmp_19 * tmp_9;
      real_t tmp_22 = tmp_4 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_23 = tmp_22 * tmp_3;
      real_t tmp_24 = tmp_15 * tmp_22;
      real_t tmp_25 = 0.2393143352496831 * tmp_17;
      real_t tmp_26 = tmp_4 * ( 0.5 * tmp_5 + tmp_6 );
      real_t tmp_27 = tmp_1 * tmp_26;
      real_t tmp_28 = tmp_26 * tmp_9;
      real_t tmp_29 = tmp_4 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_30 = tmp_29 * tmp_3;
      real_t tmp_31 = tmp_15 * tmp_29;
      real_t tmp_32 = 0.2844444444444445 * tmp_17;
      real_t tmp_33 = tmp_4 * ( 0.7692346550528415 * tmp_5 + tmp_6 );
      real_t tmp_34 = tmp_1 * tmp_33;
      real_t tmp_35 = tmp_33 * tmp_9;
      real_t tmp_36 = tmp_4 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_37 = tmp_3 * tmp_36;
      real_t tmp_38 = tmp_15 * tmp_36;
      real_t tmp_39 = 0.2393143352496831 * tmp_17;
      real_t tmp_40 = tmp_4 * ( 0.95308992296933193 * tmp_5 + tmp_6 );
      real_t tmp_41 = tmp_1 * tmp_40;
      real_t tmp_42 = tmp_40 * tmp_9;
      real_t tmp_43 = tmp_4 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_44 = tmp_3 * tmp_43;
      real_t tmp_45 = tmp_15 * tmp_43;
      real_t tmp_46 = 0.11846344252809471 * tmp_17;
      real_t a_0_0  = tmp_18 * ( -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1 ) + tmp_25 * ( -tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1 ) +
                     tmp_32 * ( -tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1 ) + tmp_39 * ( -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1 ) +
                     tmp_46 * ( -tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1 );
      real_t a_1_0 = tmp_18 * ( tmp_10 + tmp_14 ) + tmp_25 * ( tmp_21 + tmp_23 ) + tmp_32 * ( tmp_28 + tmp_30 ) +
                     tmp_39 * ( tmp_35 + tmp_37 ) + tmp_46 * ( tmp_42 + tmp_44 );
      real_t a_2_0 = tmp_18 * ( tmp_16 + tmp_8 ) + tmp_25 * ( tmp_20 + tmp_24 ) + tmp_32 * ( tmp_27 + tmp_31 ) +
                     tmp_39 * ( tmp_34 + tmp_38 ) + tmp_46 * ( tmp_41 + tmp_45 );
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + tmp_2;
      real_t tmp_7  = tmp_4 * ( 0.046910077030668018 * tmp_5 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = p_affine_10_0 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_18 = 0.11846344252809471 * tmp_17;
      real_t tmp_19 = tmp_4 * ( 0.23076534494715845 * tmp_5 + tmp_6 );
      real_t tmp_20 = tmp_1 * tmp_19;
      real_t tmp_21 = tmp_19 * tmp_9;
      real_t tmp_22 = tmp_4 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_23 = tmp_22 * tmp_3;
      real_t tmp_24 = tmp_15 * tmp_22;
      real_t tmp_25 = 0.2393143352496831 * tmp_17;
      real_t tmp_26 = tmp_4 * ( 0.5 * tmp_5 + tmp_6 );
      real_t tmp_27 = tmp_1 * tmp_26;
      real_t tmp_28 = tmp_26 * tmp_9;
      real_t tmp_29 = tmp_4 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_30 = tmp_29 * tmp_3;
      real_t tmp_31 = tmp_15 * tmp_29;
      real_t tmp_32 = 0.2844444444444445 * tmp_17;
      real_t tmp_33 = tmp_4 * ( 0.7692346550528415 * tmp_5 + tmp_6 );
      real_t tmp_34 = tmp_1 * tmp_33;
      real_t tmp_35 = tmp_33 * tmp_9;
      real_t tmp_36 = tmp_4 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_37 = tmp_3 * tmp_36;
      real_t tmp_38 = tmp_15 * tmp_36;
      real_t tmp_39 = 0.2393143352496831 * tmp_17;
      real_t tmp_40 = tmp_4 * ( 0.95308992296933193 * tmp_5 + tmp_6 );
      real_t tmp_41 = tmp_1 * tmp_40;
      real_t tmp_42 = tmp_40 * tmp_9;
      real_t tmp_43 = tmp_4 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_44 = tmp_3 * tmp_43;
      real_t tmp_45 = tmp_15 * tmp_43;
      real_t tmp_46 = 0.11846344252809471 * tmp_17;
      real_t a_0_0  = tmp_18 * ( -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1 ) + tmp_25 * ( -tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1 ) +
                     tmp_32 * ( -tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1 ) + tmp_39 * ( -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1 ) +
                     tmp_46 * ( -tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1 );
      real_t a_1_0 = tmp_18 * ( tmp_10 + tmp_14 ) + tmp_25 * ( tmp_21 + tmp_23 ) + tmp_32 * ( tmp_28 + tmp_30 ) +
                     tmp_39 * ( tmp_35 + tmp_37 ) + tmp_46 * ( tmp_42 + tmp_44 );
      real_t a_2_0 = tmp_18 * ( tmp_16 + tmp_8 ) + tmp_25 * ( tmp_20 + tmp_24 ) + tmp_32 * ( tmp_27 + tmp_31 ) +
                     tmp_39 * ( tmp_34 + tmp_38 ) + tmp_46 * ( tmp_41 + tmp_45 );
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

class DGDivtFormP1P0_1 : public hyteg::dg::DGForm2D
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
      real_t tmp_3  = 1.0 / ( tmp_1 * ( p_affine_2_1 + tmp_2 ) - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = tmp_3 * ( p_affine_0_0 - p_affine_2_0 );
      real_t tmp_6  = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                               p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_7  = tmp_6 * ( tmp_4 + tmp_5 );
      real_t tmp_8  = tmp_5 * tmp_6;
      real_t tmp_9  = tmp_4 * tmp_6;
      real_t a_0_0  = 0.5 * tmp_7;
      real_t a_1_0  = -0.5 * tmp_8;
      real_t a_2_0  = -0.5 * tmp_9;
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
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + tmp_2;
      real_t tmp_7  = tmp_4 * ( 0.046910077030668018 * tmp_5 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = 0.5 * p_affine_10_1 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_18 = 0.11846344252809471 * tmp_17;
      real_t tmp_19 = tmp_4 * ( 0.23076534494715845 * tmp_5 + tmp_6 );
      real_t tmp_20 = tmp_1 * tmp_19;
      real_t tmp_21 = tmp_19 * tmp_9;
      real_t tmp_22 = tmp_4 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_23 = tmp_22 * tmp_3;
      real_t tmp_24 = tmp_15 * tmp_22;
      real_t tmp_25 = 0.2393143352496831 * tmp_17;
      real_t tmp_26 = tmp_4 * ( 0.5 * tmp_5 + tmp_6 );
      real_t tmp_27 = tmp_1 * tmp_26;
      real_t tmp_28 = tmp_26 * tmp_9;
      real_t tmp_29 = tmp_4 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_30 = tmp_29 * tmp_3;
      real_t tmp_31 = tmp_15 * tmp_29;
      real_t tmp_32 = 0.2844444444444445 * tmp_17;
      real_t tmp_33 = tmp_4 * ( 0.7692346550528415 * tmp_5 + tmp_6 );
      real_t tmp_34 = tmp_1 * tmp_33;
      real_t tmp_35 = tmp_33 * tmp_9;
      real_t tmp_36 = tmp_4 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_37 = tmp_3 * tmp_36;
      real_t tmp_38 = tmp_15 * tmp_36;
      real_t tmp_39 = 0.2393143352496831 * tmp_17;
      real_t tmp_40 = tmp_4 * ( 0.95308992296933193 * tmp_5 + tmp_6 );
      real_t tmp_41 = tmp_1 * tmp_40;
      real_t tmp_42 = tmp_40 * tmp_9;
      real_t tmp_43 = tmp_4 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_44 = tmp_3 * tmp_43;
      real_t tmp_45 = tmp_15 * tmp_43;
      real_t tmp_46 = 0.11846344252809471 * tmp_17;
      real_t a_0_0  = tmp_18 * ( -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1 ) + tmp_25 * ( -tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1 ) +
                     tmp_32 * ( -tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1 ) + tmp_39 * ( -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1 ) +
                     tmp_46 * ( -tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1 );
      real_t a_1_0 = tmp_18 * ( tmp_10 + tmp_14 ) + tmp_25 * ( tmp_21 + tmp_23 ) + tmp_32 * ( tmp_28 + tmp_30 ) +
                     tmp_39 * ( tmp_35 + tmp_37 ) + tmp_46 * ( tmp_42 + tmp_44 );
      real_t a_2_0 = tmp_18 * ( tmp_16 + tmp_8 ) + tmp_25 * ( tmp_20 + tmp_24 ) + tmp_32 * ( tmp_27 + tmp_31 ) +
                     tmp_39 * ( tmp_34 + tmp_38 ) + tmp_46 * ( tmp_41 + tmp_45 );
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + tmp_2;
      real_t tmp_7  = tmp_4 * ( 0.046910077030668018 * tmp_5 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = 0.5 * p_affine_10_1 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_18 = 0.11846344252809471 * tmp_17;
      real_t tmp_19 = tmp_4 * ( 0.23076534494715845 * tmp_5 + tmp_6 );
      real_t tmp_20 = tmp_1 * tmp_19;
      real_t tmp_21 = tmp_19 * tmp_9;
      real_t tmp_22 = tmp_4 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_23 = tmp_22 * tmp_3;
      real_t tmp_24 = tmp_15 * tmp_22;
      real_t tmp_25 = 0.2393143352496831 * tmp_17;
      real_t tmp_26 = tmp_4 * ( 0.5 * tmp_5 + tmp_6 );
      real_t tmp_27 = tmp_1 * tmp_26;
      real_t tmp_28 = tmp_26 * tmp_9;
      real_t tmp_29 = tmp_4 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_30 = tmp_29 * tmp_3;
      real_t tmp_31 = tmp_15 * tmp_29;
      real_t tmp_32 = 0.2844444444444445 * tmp_17;
      real_t tmp_33 = tmp_4 * ( 0.7692346550528415 * tmp_5 + tmp_6 );
      real_t tmp_34 = tmp_1 * tmp_33;
      real_t tmp_35 = tmp_33 * tmp_9;
      real_t tmp_36 = tmp_4 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_37 = tmp_3 * tmp_36;
      real_t tmp_38 = tmp_15 * tmp_36;
      real_t tmp_39 = 0.2393143352496831 * tmp_17;
      real_t tmp_40 = tmp_4 * ( 0.95308992296933193 * tmp_5 + tmp_6 );
      real_t tmp_41 = tmp_1 * tmp_40;
      real_t tmp_42 = tmp_40 * tmp_9;
      real_t tmp_43 = tmp_4 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_44 = tmp_3 * tmp_43;
      real_t tmp_45 = tmp_15 * tmp_43;
      real_t tmp_46 = 0.11846344252809471 * tmp_17;
      real_t a_0_0  = tmp_18 * ( -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1 ) + tmp_25 * ( -tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1 ) +
                     tmp_32 * ( -tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1 ) + tmp_39 * ( -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1 ) +
                     tmp_46 * ( -tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1 );
      real_t a_1_0 = tmp_18 * ( tmp_10 + tmp_14 ) + tmp_25 * ( tmp_21 + tmp_23 ) + tmp_32 * ( tmp_28 + tmp_30 ) +
                     tmp_39 * ( tmp_35 + tmp_37 ) + tmp_46 * ( tmp_42 + tmp_44 );
      real_t a_2_0 = tmp_18 * ( tmp_16 + tmp_8 ) + tmp_25 * ( tmp_20 + tmp_24 ) + tmp_32 * ( tmp_27 + tmp_31 ) +
                     tmp_39 * ( tmp_34 + tmp_38 ) + tmp_46 * ( tmp_41 + tmp_45 );
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

      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = 1.0 / ( tmp_1 * tmp_3 - ( p_affine_1_1 + tmp_2 ) * ( p_affine_2_0 + tmp_0 ) );
      real_t tmp_5  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_6  = p_affine_6_1 + tmp_2;
      real_t tmp_7  = tmp_4 * ( 0.046910077030668018 * tmp_5 + tmp_6 );
      real_t tmp_8  = tmp_1 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_7 * tmp_9;
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_0;
      real_t tmp_13 = tmp_4 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_13 * tmp_3;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_13 * tmp_15;
      real_t tmp_17 = p_affine_10_1 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_5 * tmp_5 ), 1.0 / 2.0 ) );
      real_t tmp_18 = 0.11846344252809471 * tmp_17;
      real_t tmp_19 = tmp_4 * ( 0.23076534494715845 * tmp_5 + tmp_6 );
      real_t tmp_20 = tmp_1 * tmp_19;
      real_t tmp_21 = tmp_19 * tmp_9;
      real_t tmp_22 = tmp_4 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_23 = tmp_22 * tmp_3;
      real_t tmp_24 = tmp_15 * tmp_22;
      real_t tmp_25 = 0.2393143352496831 * tmp_17;
      real_t tmp_26 = tmp_4 * ( 0.5 * tmp_5 + tmp_6 );
      real_t tmp_27 = tmp_1 * tmp_26;
      real_t tmp_28 = tmp_26 * tmp_9;
      real_t tmp_29 = tmp_4 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_30 = tmp_29 * tmp_3;
      real_t tmp_31 = tmp_15 * tmp_29;
      real_t tmp_32 = 0.2844444444444445 * tmp_17;
      real_t tmp_33 = tmp_4 * ( 0.7692346550528415 * tmp_5 + tmp_6 );
      real_t tmp_34 = tmp_1 * tmp_33;
      real_t tmp_35 = tmp_33 * tmp_9;
      real_t tmp_36 = tmp_4 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_37 = tmp_3 * tmp_36;
      real_t tmp_38 = tmp_15 * tmp_36;
      real_t tmp_39 = 0.2393143352496831 * tmp_17;
      real_t tmp_40 = tmp_4 * ( 0.95308992296933193 * tmp_5 + tmp_6 );
      real_t tmp_41 = tmp_1 * tmp_40;
      real_t tmp_42 = tmp_40 * tmp_9;
      real_t tmp_43 = tmp_4 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_44 = tmp_3 * tmp_43;
      real_t tmp_45 = tmp_15 * tmp_43;
      real_t tmp_46 = 0.11846344252809471 * tmp_17;
      real_t a_0_0  = tmp_18 * ( -tmp_10 - tmp_14 - tmp_16 - tmp_8 + 1 ) + tmp_25 * ( -tmp_20 - tmp_21 - tmp_23 - tmp_24 + 1 ) +
                     tmp_32 * ( -tmp_27 - tmp_28 - tmp_30 - tmp_31 + 1 ) + tmp_39 * ( -tmp_34 - tmp_35 - tmp_37 - tmp_38 + 1 ) +
                     tmp_46 * ( -tmp_41 - tmp_42 - tmp_44 - tmp_45 + 1 );
      real_t a_1_0 = tmp_18 * ( tmp_10 + tmp_14 ) + tmp_25 * ( tmp_21 + tmp_23 ) + tmp_32 * ( tmp_28 + tmp_30 ) +
                     tmp_39 * ( tmp_35 + tmp_37 ) + tmp_46 * ( tmp_42 + tmp_44 );
      real_t a_2_0 = tmp_18 * ( tmp_16 + tmp_8 ) + tmp_25 * ( tmp_20 + tmp_24 ) + tmp_32 * ( tmp_27 + tmp_31 ) +
                     tmp_39 * ( tmp_34 + tmp_38 ) + tmp_46 * ( tmp_41 + tmp_45 );
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
} // namespace dg
} // namespace hyteg

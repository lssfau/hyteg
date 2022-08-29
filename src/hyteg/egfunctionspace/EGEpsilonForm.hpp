
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
namespace eg {

class EGEpsilonFormEP1_0 : public hyteg::dg::DGForm2D
{
 public:
   EGEpsilonFormEP1_0( std::function< real_t( const Point3D& ) > _callback2D )
   : callback2D( _callback2D )
   {}

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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id5  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id6  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id7  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id8  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id9  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id11 = 0;
      Scalar_Variable_Coefficient_2D_mu(
          0.063089014491502282 * p_affine_0_0 + 0.063089014491502227 * p_affine_1_0 + 0.87382197101699555 * p_affine_2_0,
          0.063089014491502282 * p_affine_0_1 + 0.063089014491502227 * p_affine_1_1 + 0.87382197101699555 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu(
          0.24928674517091043 * p_affine_0_0 + 0.24928674517091043 * p_affine_1_0 + 0.50142650965817914 * p_affine_2_0,
          0.24928674517091043 * p_affine_0_1 + 0.24928674517091043 * p_affine_1_1 + 0.50142650965817914 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu(
          0.63650249912139867 * p_affine_0_0 + 0.31035245103378439 * p_affine_1_0 + 0.053145049844816938 * p_affine_2_0,
          0.63650249912139867 * p_affine_0_1 + 0.31035245103378439 * p_affine_1_1 + 0.053145049844816938 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu(
          0.053145049844816938 * p_affine_0_0 + 0.63650249912139867 * p_affine_1_0 + 0.31035245103378439 * p_affine_2_0,
          0.053145049844816938 * p_affine_0_1 + 0.63650249912139867 * p_affine_1_1 + 0.31035245103378439 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu(
          0.063089014491502227 * p_affine_0_0 + 0.87382197101699555 * p_affine_1_0 + 0.063089014491502227 * p_affine_2_0,
          0.063089014491502227 * p_affine_0_1 + 0.87382197101699555 * p_affine_1_1 + 0.063089014491502227 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      Scalar_Variable_Coefficient_2D_mu(
          0.24928674517091043 * p_affine_0_0 + 0.50142650965817914 * p_affine_1_0 + 0.24928674517091043 * p_affine_2_0,
          0.24928674517091043 * p_affine_0_1 + 0.50142650965817914 * p_affine_1_1 + 0.24928674517091043 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id5 );
      Scalar_Variable_Coefficient_2D_mu(
          0.87382197101699566 * p_affine_0_0 + 0.063089014491502227 * p_affine_1_0 + 0.063089014491502227 * p_affine_2_0,
          0.87382197101699566 * p_affine_0_1 + 0.063089014491502227 * p_affine_1_1 + 0.063089014491502227 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id6 );
      Scalar_Variable_Coefficient_2D_mu(
          0.50142650965817914 * p_affine_0_0 + 0.24928674517091043 * p_affine_1_0 + 0.24928674517091043 * p_affine_2_0,
          0.50142650965817914 * p_affine_0_1 + 0.24928674517091043 * p_affine_1_1 + 0.24928674517091043 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id7 );
      Scalar_Variable_Coefficient_2D_mu(
          0.053145049844816938 * p_affine_0_0 + 0.31035245103378439 * p_affine_1_0 + 0.63650249912139867 * p_affine_2_0,
          0.053145049844816938 * p_affine_0_1 + 0.31035245103378439 * p_affine_1_1 + 0.63650249912139867 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id8 );
      Scalar_Variable_Coefficient_2D_mu(
          0.63650249912139867 * p_affine_0_0 + 0.053145049844816938 * p_affine_1_0 + 0.31035245103378439 * p_affine_2_0,
          0.63650249912139867 * p_affine_0_1 + 0.053145049844816938 * p_affine_1_1 + 0.31035245103378439 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id9 );
      Scalar_Variable_Coefficient_2D_mu(
          0.31035245103378439 * p_affine_0_0 + 0.63650249912139867 * p_affine_1_0 + 0.053145049844816938 * p_affine_2_0,
          0.31035245103378439 * p_affine_0_1 + 0.63650249912139867 * p_affine_1_1 + 0.053145049844816938 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id10 );
      Scalar_Variable_Coefficient_2D_mu(
          0.31035245103378439 * p_affine_0_0 + 0.053145049844816938 * p_affine_1_0 + 0.63650249912139867 * p_affine_2_0,
          0.31035245103378439 * p_affine_0_1 + 0.053145049844816938 * p_affine_1_1 + 0.63650249912139867 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id11 );
      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_0;
      real_t tmp_6  = p_affine_1_1 + tmp_2;
      real_t tmp_7  = 1.0 / ( tmp_4 - tmp_5 * tmp_6 );
      real_t tmp_8  = 1.0 * tmp_7;
      real_t tmp_9  = p_affine_0_1 - p_affine_1_1;
      real_t tmp_10 = tmp_8 * tmp_9;
      real_t tmp_11 = tmp_10 * tmp_5 + tmp_4 * tmp_8;
      real_t tmp_12 = Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_11;
      real_t tmp_13 = -2 * tmp_10 - 2 * tmp_3 * tmp_8;
      real_t tmp_14 = 0.5 * tmp_7;
      real_t tmp_15 = tmp_1 * tmp_14;
      real_t tmp_16 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_17 = tmp_14 * tmp_16;
      real_t tmp_18 = tmp_14 * tmp_3;
      real_t tmp_19 = tmp_1 * tmp_17 + tmp_15 * tmp_5 + tmp_18 * tmp_6 + tmp_18 * tmp_9;
      real_t tmp_20 = Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_19;
      real_t tmp_21 = -4 * tmp_15 - 4 * tmp_17;
      real_t tmp_22 = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                                p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_23 = 0.025422453185103409 * tmp_22;
      real_t tmp_24 = tmp_11 * tmp_13;
      real_t tmp_25 = tmp_19 * tmp_21;
      real_t tmp_26 = 0.058393137863189684 * tmp_22;
      real_t tmp_27 = 0.041425537809186785 * tmp_22;
      real_t tmp_28 = 0.041425537809186785 * tmp_22;
      real_t tmp_29 = 0.025422453185103409 * tmp_22;
      real_t tmp_30 = 0.058393137863189684 * tmp_22;
      real_t tmp_31 = 0.025422453185103409 * tmp_22;
      real_t tmp_32 = 0.058393137863189684 * tmp_22;
      real_t tmp_33 = 0.041425537809186785 * tmp_22;
      real_t tmp_34 = 0.041425537809186785 * tmp_22;
      real_t tmp_35 = 0.041425537809186785 * tmp_22;
      real_t tmp_36 = 0.041425537809186785 * tmp_22;
      real_t tmp_37 = 2.0 * tmp_7;
      real_t tmp_38 = tmp_3 * tmp_37;
      real_t tmp_39 = tmp_16 * tmp_37;
      real_t tmp_40 = tmp_11 * tmp_38;
      real_t tmp_41 = tmp_19 * tmp_39;
      real_t tmp_42 = tmp_37 * tmp_9;
      real_t tmp_43 = tmp_1 * tmp_37;
      real_t tmp_44 = tmp_11 * tmp_42;
      real_t tmp_45 = tmp_19 * tmp_43;
      real_t a_0_0 =
          tmp_23 * ( tmp_12 * tmp_13 + tmp_20 * tmp_21 ) +
          tmp_26 * ( Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_25 ) +
          tmp_27 * ( Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_25 ) +
          tmp_28 * ( Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_25 ) +
          tmp_29 * ( Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_25 ) +
          tmp_30 * ( Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_25 ) +
          tmp_31 * ( Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_25 ) +
          tmp_32 * ( Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_25 ) +
          tmp_33 * ( Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_25 ) +
          tmp_34 * ( Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_25 ) +
          tmp_35 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_25 ) +
          tmp_36 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_25 );
      real_t a_0_1 =
          tmp_23 * ( tmp_12 * tmp_38 + tmp_20 * tmp_39 ) +
          tmp_26 * ( Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_41 ) +
          tmp_27 * ( Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_41 ) +
          tmp_28 * ( Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_41 ) +
          tmp_29 * ( Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_41 ) +
          tmp_30 * ( Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_41 ) +
          tmp_31 * ( Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_41 ) +
          tmp_32 * ( Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_41 ) +
          tmp_33 * ( Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_41 ) +
          tmp_34 * ( Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_41 ) +
          tmp_35 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_41 ) +
          tmp_36 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_41 );
      real_t a_0_2 =
          tmp_23 * ( tmp_12 * tmp_42 + tmp_20 * tmp_43 ) +
          tmp_26 * ( Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_45 ) +
          tmp_27 * ( Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_45 ) +
          tmp_28 * ( Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_45 ) +
          tmp_29 * ( Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_45 ) +
          tmp_30 * ( Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_45 ) +
          tmp_31 * ( Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_45 ) +
          tmp_32 * ( Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_45 ) +
          tmp_33 * ( Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_45 ) +
          tmp_34 * ( Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_45 ) +
          tmp_35 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_45 ) +
          tmp_36 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_45 );
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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.95308992296933193 * p_affine_6_0 + 0.046910077030668018 * p_affine_7_0,
                                         0.95308992296933193 * p_affine_6_1 + 0.046910077030668018 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.7692346550528415 * p_affine_6_0 + 0.23076534494715845 * p_affine_7_0,
                                         0.7692346550528415 * p_affine_6_1 + 0.23076534494715845 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.5 * p_affine_6_0 + 0.5 * p_affine_7_0,
                                         0.5 * p_affine_6_1 + 0.5 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.2307653449471585 * p_affine_6_0 + 0.7692346550528415 * p_affine_7_0,
                                         0.2307653449471585 * p_affine_6_1 + 0.7692346550528415 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.046910077030668074 * p_affine_6_0 + 0.95308992296933193 * p_affine_7_0,
                                         0.046910077030668074 * p_affine_6_1 + 0.95308992296933193 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_1_1 + tmp_0;
      real_t tmp_2  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3  = -p_affine_0_0;
      real_t tmp_4  = p_affine_1_0 + tmp_3;
      real_t tmp_5  = p_affine_2_1 + tmp_0;
      real_t tmp_6  = tmp_4 * tmp_5;
      real_t tmp_7  = p_affine_2_0 + tmp_3;
      real_t tmp_8  = 1.0 / ( -tmp_1 * tmp_7 + tmp_6 );
      real_t tmp_9  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + tmp_0;
      real_t tmp_11 = tmp_8 * ( tmp_10 + 0.046910077030668018 * tmp_9 );
      real_t tmp_12 = tmp_11 * tmp_2;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_3;
      real_t tmp_15 = tmp_8 * ( 0.046910077030668018 * tmp_13 + tmp_14 );
      real_t tmp_16 = tmp_15 * tmp_5;
      real_t tmp_17 = tmp_12 + tmp_16;
      real_t tmp_18 = tmp_17 - 1.0 / 3.0;
      real_t tmp_19 = tmp_11 * tmp_4;
      real_t tmp_20 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_21 = tmp_15 * tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_22 - 1.0 / 3.0;
      real_t tmp_24 = p_affine_10_0 * ( tmp_1 * tmp_18 + tmp_23 * tmp_5 );
      real_t tmp_25 = 0.5 * tmp_8;
      real_t tmp_26 = tmp_25 * tmp_4;
      real_t tmp_27 = tmp_2 * tmp_25;
      real_t tmp_28 = -tmp_26 - tmp_27;
      real_t tmp_29 = 0.5 * tmp_28;
      real_t tmp_30 = tmp_18 * tmp_4 + tmp_23 * tmp_7;
      real_t tmp_31 = 1.0 * tmp_8;
      real_t tmp_32 = tmp_31 * tmp_5;
      real_t tmp_33 = tmp_20 * tmp_31;
      real_t tmp_34 = 0.5 * p_affine_10_0 * ( -tmp_32 - tmp_33 ) + 0.5 * p_affine_10_1 * tmp_28;
      real_t tmp_35 = -tmp_12 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_36 = std::abs( std::pow( ( tmp_13 * tmp_13 ) + ( tmp_9 * tmp_9 ), 1.0 / 2.0 ) );
      real_t tmp_37 = 6 / tmp_36;
      real_t tmp_38 = tmp_30 * tmp_37;
      real_t tmp_39 = tmp_25 * tmp_5;
      real_t tmp_40 = 0.5 * p_affine_10_0 * ( tmp_31 * tmp_6 + tmp_33 * tmp_7 ) +
                      0.5 * p_affine_10_1 * ( tmp_1 * tmp_39 + tmp_20 * tmp_39 + tmp_26 * tmp_7 + tmp_27 * tmp_4 );
      real_t tmp_41  = 2 * tmp_36;
      real_t tmp_42  = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_41;
      real_t tmp_43  = tmp_8 * ( tmp_10 + 0.23076534494715845 * tmp_9 );
      real_t tmp_44  = tmp_2 * tmp_43;
      real_t tmp_45  = tmp_8 * ( 0.23076534494715845 * tmp_13 + tmp_14 );
      real_t tmp_46  = tmp_45 * tmp_5;
      real_t tmp_47  = tmp_44 + tmp_46;
      real_t tmp_48  = tmp_47 - 1.0 / 3.0;
      real_t tmp_49  = tmp_4 * tmp_43;
      real_t tmp_50  = tmp_20 * tmp_45;
      real_t tmp_51  = tmp_49 + tmp_50;
      real_t tmp_52  = tmp_51 - 1.0 / 3.0;
      real_t tmp_53  = tmp_1 * tmp_48 + tmp_5 * tmp_52;
      real_t tmp_54  = p_affine_10_0 * tmp_29;
      real_t tmp_55  = tmp_4 * tmp_48 + tmp_52 * tmp_7;
      real_t tmp_56  = -tmp_44 - tmp_46 - tmp_49 - tmp_50 + 1;
      real_t tmp_57  = tmp_37 * tmp_55;
      real_t tmp_58  = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_41;
      real_t tmp_59  = tmp_8 * ( tmp_10 + 0.5 * tmp_9 );
      real_t tmp_60  = tmp_2 * tmp_59;
      real_t tmp_61  = tmp_8 * ( 0.5 * tmp_13 + tmp_14 );
      real_t tmp_62  = tmp_5 * tmp_61;
      real_t tmp_63  = tmp_60 + tmp_62;
      real_t tmp_64  = tmp_63 - 1.0 / 3.0;
      real_t tmp_65  = tmp_4 * tmp_59;
      real_t tmp_66  = tmp_20 * tmp_61;
      real_t tmp_67  = tmp_65 + tmp_66;
      real_t tmp_68  = tmp_67 - 1.0 / 3.0;
      real_t tmp_69  = tmp_1 * tmp_64 + tmp_5 * tmp_68;
      real_t tmp_70  = tmp_4 * tmp_64 + tmp_68 * tmp_7;
      real_t tmp_71  = -tmp_60 - tmp_62 - tmp_65 - tmp_66 + 1;
      real_t tmp_72  = tmp_37 * tmp_70;
      real_t tmp_73  = 0.2844444444444445 * Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_41;
      real_t tmp_74  = tmp_8 * ( tmp_10 + 0.7692346550528415 * tmp_9 );
      real_t tmp_75  = tmp_2 * tmp_74;
      real_t tmp_76  = tmp_8 * ( 0.7692346550528415 * tmp_13 + tmp_14 );
      real_t tmp_77  = tmp_5 * tmp_76;
      real_t tmp_78  = tmp_75 + tmp_77;
      real_t tmp_79  = tmp_78 - 1.0 / 3.0;
      real_t tmp_80  = tmp_4 * tmp_74;
      real_t tmp_81  = tmp_20 * tmp_76;
      real_t tmp_82  = tmp_80 + tmp_81;
      real_t tmp_83  = tmp_82 - 1.0 / 3.0;
      real_t tmp_84  = tmp_1 * tmp_79 + tmp_5 * tmp_83;
      real_t tmp_85  = tmp_4 * tmp_79 + tmp_7 * tmp_83;
      real_t tmp_86  = -tmp_75 - tmp_77 - tmp_80 - tmp_81 + 1;
      real_t tmp_87  = tmp_37 * tmp_85;
      real_t tmp_88  = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_41;
      real_t tmp_89  = tmp_8 * ( tmp_10 + 0.95308992296933193 * tmp_9 );
      real_t tmp_90  = tmp_2 * tmp_89;
      real_t tmp_91  = tmp_8 * ( 0.95308992296933193 * tmp_13 + tmp_14 );
      real_t tmp_92  = tmp_5 * tmp_91;
      real_t tmp_93  = tmp_90 + tmp_92;
      real_t tmp_94  = tmp_93 - 1.0 / 3.0;
      real_t tmp_95  = tmp_4 * tmp_89;
      real_t tmp_96  = tmp_20 * tmp_91;
      real_t tmp_97  = tmp_95 + tmp_96;
      real_t tmp_98  = tmp_97 - 1.0 / 3.0;
      real_t tmp_99  = tmp_1 * tmp_94 + tmp_5 * tmp_98;
      real_t tmp_100 = tmp_4 * tmp_94 + tmp_7 * tmp_98;
      real_t tmp_101 = -tmp_90 - tmp_92 - tmp_95 - tmp_96 + 1;
      real_t tmp_102 = tmp_100 * tmp_37;
      real_t tmp_103 = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_41;
      real_t tmp_104 = 0.25 * tmp_8;
      real_t tmp_105 = tmp_104 * tmp_2;
      real_t tmp_106 = 0.5 * p_affine_10_0 * tmp_32 + 0.5 * p_affine_10_1 * tmp_27;
      real_t tmp_107 = p_affine_10_0 * tmp_105;
      real_t tmp_108 = tmp_104 * tmp_4;
      real_t tmp_109 = 0.5 * p_affine_10_0 * tmp_33 + 0.5 * p_affine_10_1 * tmp_26;
      real_t tmp_110 = p_affine_10_0 * tmp_108;
      real_t a_0_0   = tmp_103 * ( -tmp_100 * tmp_34 + tmp_101 * tmp_102 - tmp_101 * tmp_40 - tmp_54 * tmp_99 ) +
                     tmp_42 * ( -tmp_24 * tmp_29 - tmp_30 * tmp_34 + tmp_35 * tmp_38 - tmp_35 * tmp_40 ) +
                     tmp_58 * ( -tmp_34 * tmp_55 - tmp_40 * tmp_56 - tmp_53 * tmp_54 + tmp_56 * tmp_57 ) +
                     tmp_73 * ( -tmp_34 * tmp_70 - tmp_40 * tmp_71 - tmp_54 * tmp_69 + tmp_71 * tmp_72 ) +
                     tmp_88 * ( -tmp_34 * tmp_85 - tmp_40 * tmp_86 - tmp_54 * tmp_84 + tmp_86 * tmp_87 );
      real_t a_0_1 = tmp_103 * ( -tmp_100 * tmp_106 + tmp_102 * tmp_93 - tmp_107 * tmp_99 - tmp_40 * tmp_93 ) +
                     tmp_42 * ( -tmp_105 * tmp_24 - tmp_106 * tmp_30 + tmp_17 * tmp_38 - tmp_17 * tmp_40 ) +
                     tmp_58 * ( -tmp_106 * tmp_55 - tmp_107 * tmp_53 - tmp_40 * tmp_47 + tmp_47 * tmp_57 ) +
                     tmp_73 * ( -tmp_106 * tmp_70 - tmp_107 * tmp_69 - tmp_40 * tmp_63 + tmp_63 * tmp_72 ) +
                     tmp_88 * ( -tmp_106 * tmp_85 - tmp_107 * tmp_84 - tmp_40 * tmp_78 + tmp_78 * tmp_87 );
      real_t a_0_2 = tmp_103 * ( -tmp_100 * tmp_109 + tmp_102 * tmp_97 - tmp_110 * tmp_99 - tmp_40 * tmp_97 ) +
                     tmp_42 * ( -tmp_108 * tmp_24 - tmp_109 * tmp_30 + tmp_22 * tmp_38 - tmp_22 * tmp_40 ) +
                     tmp_58 * ( -tmp_109 * tmp_55 - tmp_110 * tmp_53 - tmp_40 * tmp_51 + tmp_51 * tmp_57 ) +
                     tmp_73 * ( -tmp_109 * tmp_70 - tmp_110 * tmp_69 - tmp_40 * tmp_67 + tmp_67 * tmp_72 ) +
                     tmp_88 * ( -tmp_109 * tmp_85 - tmp_110 * tmp_84 - tmp_40 * tmp_82 + tmp_82 * tmp_87 );
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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.95308992296933193 * p_affine_6_0 + 0.046910077030668018 * p_affine_7_0,
                                         0.95308992296933193 * p_affine_6_1 + 0.046910077030668018 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.7692346550528415 * p_affine_6_0 + 0.23076534494715845 * p_affine_7_0,
                                         0.7692346550528415 * p_affine_6_1 + 0.23076534494715845 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.5 * p_affine_6_0 + 0.5 * p_affine_7_0,
                                         0.5 * p_affine_6_1 + 0.5 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.2307653449471585 * p_affine_6_0 + 0.7692346550528415 * p_affine_7_0,
                                         0.2307653449471585 * p_affine_6_1 + 0.7692346550528415 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.046910077030668074 * p_affine_6_0 + 0.95308992296933193 * p_affine_7_0,
                                         0.046910077030668074 * p_affine_6_1 + 0.95308992296933193 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_1_1 + tmp_0;
      real_t tmp_2  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3  = -p_affine_0_0;
      real_t tmp_4  = p_affine_1_0 + tmp_3;
      real_t tmp_5  = p_affine_2_1 + tmp_0;
      real_t tmp_6  = tmp_4 * tmp_5;
      real_t tmp_7  = p_affine_2_0 + tmp_3;
      real_t tmp_8  = 1.0 / ( -tmp_1 * tmp_7 + tmp_6 );
      real_t tmp_9  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + 0.046910077030668018 * tmp_9;
      real_t tmp_11 = tmp_8 * ( tmp_0 + tmp_10 );
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + 0.046910077030668018 * tmp_12;
      real_t tmp_14 = tmp_8 * ( tmp_13 + tmp_3 );
      real_t tmp_15 = tmp_11 * tmp_2 + tmp_14 * tmp_5 - 1.0 / 3.0;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_11 * tmp_4 + tmp_14 * tmp_16 - 1.0 / 3.0;
      real_t tmp_18 = p_affine_10_0 * ( tmp_1 * tmp_15 + tmp_17 * tmp_5 );
      real_t tmp_19 = -p_affine_3_0;
      real_t tmp_20 = p_affine_4_0 + tmp_19;
      real_t tmp_21 = -p_affine_3_1;
      real_t tmp_22 = p_affine_5_1 + tmp_21;
      real_t tmp_23 = 1.0 / ( tmp_20 * tmp_22 - ( p_affine_4_1 + tmp_21 ) * ( p_affine_5_0 + tmp_19 ) );
      real_t tmp_24 = 0.5 * tmp_23;
      real_t tmp_25 = tmp_20 * tmp_24;
      real_t tmp_26 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_27 = tmp_24 * tmp_26;
      real_t tmp_28 = -tmp_25 - tmp_27;
      real_t tmp_29 = 0.5 * tmp_28;
      real_t tmp_30 = tmp_15 * tmp_4 + tmp_17 * tmp_7;
      real_t tmp_31 = 1.0 * tmp_23;
      real_t tmp_32 = tmp_22 * tmp_31;
      real_t tmp_33 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_34 = tmp_31 * tmp_33;
      real_t tmp_35 = 0.5 * p_affine_10_0 * ( -tmp_32 - tmp_34 ) + 0.5 * p_affine_10_1 * tmp_28;
      real_t tmp_36 = tmp_23 * ( tmp_10 + tmp_21 );
      real_t tmp_37 = tmp_20 * tmp_36;
      real_t tmp_38 = tmp_26 * tmp_36;
      real_t tmp_39 = tmp_23 * ( tmp_13 + tmp_19 );
      real_t tmp_40 = tmp_22 * tmp_39;
      real_t tmp_41 = tmp_33 * tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = std::abs( std::pow( ( tmp_12 * tmp_12 ) + ( tmp_9 * tmp_9 ), 1.0 / 2.0 ) );
      real_t tmp_44 = 6 / tmp_43;
      real_t tmp_45 = tmp_30 * tmp_44;
      real_t tmp_46 = 1.0 * tmp_8;
      real_t tmp_47 = 0.5 * tmp_8;
      real_t tmp_48 = tmp_4 * tmp_47;
      real_t tmp_49 = tmp_47 * tmp_5;
      real_t tmp_50 = 0.5 * p_affine_10_0 * ( tmp_16 * tmp_46 * tmp_7 + tmp_46 * tmp_6 ) +
                      0.5 * p_affine_10_1 * ( tmp_1 * tmp_49 + tmp_16 * tmp_49 + tmp_2 * tmp_48 + tmp_48 * tmp_7 );
      real_t tmp_51  = 2 * tmp_43;
      real_t tmp_52  = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_51;
      real_t tmp_53  = p_affine_6_1 + 0.23076534494715845 * tmp_9;
      real_t tmp_54  = tmp_8 * ( tmp_0 + tmp_53 );
      real_t tmp_55  = p_affine_6_0 + 0.23076534494715845 * tmp_12;
      real_t tmp_56  = tmp_8 * ( tmp_3 + tmp_55 );
      real_t tmp_57  = tmp_2 * tmp_54 + tmp_5 * tmp_56 - 1.0 / 3.0;
      real_t tmp_58  = tmp_16 * tmp_56 + tmp_4 * tmp_54 - 1.0 / 3.0;
      real_t tmp_59  = tmp_1 * tmp_57 + tmp_5 * tmp_58;
      real_t tmp_60  = p_affine_10_0 * tmp_29;
      real_t tmp_61  = tmp_4 * tmp_57 + tmp_58 * tmp_7;
      real_t tmp_62  = tmp_23 * ( tmp_21 + tmp_53 );
      real_t tmp_63  = tmp_20 * tmp_62;
      real_t tmp_64  = tmp_26 * tmp_62;
      real_t tmp_65  = tmp_23 * ( tmp_19 + tmp_55 );
      real_t tmp_66  = tmp_22 * tmp_65;
      real_t tmp_67  = tmp_33 * tmp_65;
      real_t tmp_68  = -tmp_63 - tmp_64 - tmp_66 - tmp_67 + 1;
      real_t tmp_69  = tmp_44 * tmp_61;
      real_t tmp_70  = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_51;
      real_t tmp_71  = p_affine_6_1 + 0.5 * tmp_9;
      real_t tmp_72  = tmp_8 * ( tmp_0 + tmp_71 );
      real_t tmp_73  = p_affine_6_0 + 0.5 * tmp_12;
      real_t tmp_74  = tmp_8 * ( tmp_3 + tmp_73 );
      real_t tmp_75  = tmp_2 * tmp_72 + tmp_5 * tmp_74 - 1.0 / 3.0;
      real_t tmp_76  = tmp_16 * tmp_74 + tmp_4 * tmp_72 - 1.0 / 3.0;
      real_t tmp_77  = tmp_1 * tmp_75 + tmp_5 * tmp_76;
      real_t tmp_78  = tmp_4 * tmp_75 + tmp_7 * tmp_76;
      real_t tmp_79  = tmp_23 * ( tmp_21 + tmp_71 );
      real_t tmp_80  = tmp_20 * tmp_79;
      real_t tmp_81  = tmp_26 * tmp_79;
      real_t tmp_82  = tmp_23 * ( tmp_19 + tmp_73 );
      real_t tmp_83  = tmp_22 * tmp_82;
      real_t tmp_84  = tmp_33 * tmp_82;
      real_t tmp_85  = -tmp_80 - tmp_81 - tmp_83 - tmp_84 + 1;
      real_t tmp_86  = tmp_44 * tmp_78;
      real_t tmp_87  = 0.2844444444444445 * Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_51;
      real_t tmp_88  = p_affine_6_1 + 0.7692346550528415 * tmp_9;
      real_t tmp_89  = tmp_8 * ( tmp_0 + tmp_88 );
      real_t tmp_90  = p_affine_6_0 + 0.7692346550528415 * tmp_12;
      real_t tmp_91  = tmp_8 * ( tmp_3 + tmp_90 );
      real_t tmp_92  = tmp_2 * tmp_89 + tmp_5 * tmp_91 - 1.0 / 3.0;
      real_t tmp_93  = tmp_16 * tmp_91 + tmp_4 * tmp_89 - 1.0 / 3.0;
      real_t tmp_94  = tmp_1 * tmp_92 + tmp_5 * tmp_93;
      real_t tmp_95  = tmp_4 * tmp_92 + tmp_7 * tmp_93;
      real_t tmp_96  = tmp_23 * ( tmp_21 + tmp_88 );
      real_t tmp_97  = tmp_20 * tmp_96;
      real_t tmp_98  = tmp_26 * tmp_96;
      real_t tmp_99  = tmp_23 * ( tmp_19 + tmp_90 );
      real_t tmp_100 = tmp_22 * tmp_99;
      real_t tmp_101 = tmp_33 * tmp_99;
      real_t tmp_102 = -tmp_100 - tmp_101 - tmp_97 - tmp_98 + 1;
      real_t tmp_103 = tmp_44 * tmp_95;
      real_t tmp_104 = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_51;
      real_t tmp_105 = p_affine_6_1 + 0.95308992296933193 * tmp_9;
      real_t tmp_106 = tmp_8 * ( tmp_0 + tmp_105 );
      real_t tmp_107 = p_affine_6_0 + 0.95308992296933193 * tmp_12;
      real_t tmp_108 = tmp_8 * ( tmp_107 + tmp_3 );
      real_t tmp_109 = tmp_106 * tmp_2 + tmp_108 * tmp_5 - 1.0 / 3.0;
      real_t tmp_110 = tmp_106 * tmp_4 + tmp_108 * tmp_16 - 1.0 / 3.0;
      real_t tmp_111 = tmp_1 * tmp_109 + tmp_110 * tmp_5;
      real_t tmp_112 = tmp_109 * tmp_4 + tmp_110 * tmp_7;
      real_t tmp_113 = tmp_23 * ( tmp_105 + tmp_21 );
      real_t tmp_114 = tmp_113 * tmp_20;
      real_t tmp_115 = tmp_113 * tmp_26;
      real_t tmp_116 = tmp_23 * ( tmp_107 + tmp_19 );
      real_t tmp_117 = tmp_116 * tmp_22;
      real_t tmp_118 = tmp_116 * tmp_33;
      real_t tmp_119 = -tmp_114 - tmp_115 - tmp_117 - tmp_118 + 1;
      real_t tmp_120 = tmp_112 * tmp_44;
      real_t tmp_121 = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_51;
      real_t tmp_122 = 0.25 * tmp_23;
      real_t tmp_123 = tmp_122 * tmp_26;
      real_t tmp_124 = 0.5 * p_affine_10_0 * tmp_32 + 0.5 * p_affine_10_1 * tmp_27;
      real_t tmp_125 = tmp_38 + tmp_40;
      real_t tmp_126 = p_affine_10_0 * tmp_123;
      real_t tmp_127 = tmp_64 + tmp_66;
      real_t tmp_128 = tmp_81 + tmp_83;
      real_t tmp_129 = tmp_100 + tmp_98;
      real_t tmp_130 = tmp_115 + tmp_117;
      real_t tmp_131 = tmp_122 * tmp_20;
      real_t tmp_132 = 0.5 * p_affine_10_0 * tmp_34 + 0.5 * p_affine_10_1 * tmp_25;
      real_t tmp_133 = tmp_37 + tmp_41;
      real_t tmp_134 = p_affine_10_0 * tmp_131;
      real_t tmp_135 = tmp_63 + tmp_67;
      real_t tmp_136 = tmp_80 + tmp_84;
      real_t tmp_137 = tmp_101 + tmp_97;
      real_t tmp_138 = tmp_114 + tmp_118;
      real_t a_0_0   = tmp_104 * ( -tmp_102 * tmp_103 + tmp_102 * tmp_50 - tmp_35 * tmp_95 - tmp_60 * tmp_94 ) +
                     tmp_121 * ( -tmp_111 * tmp_60 - tmp_112 * tmp_35 - tmp_119 * tmp_120 + tmp_119 * tmp_50 ) +
                     tmp_52 * ( -tmp_18 * tmp_29 - tmp_30 * tmp_35 - tmp_42 * tmp_45 + tmp_42 * tmp_50 ) +
                     tmp_70 * ( -tmp_35 * tmp_61 + tmp_50 * tmp_68 - tmp_59 * tmp_60 - tmp_68 * tmp_69 ) +
                     tmp_87 * ( -tmp_35 * tmp_78 + tmp_50 * tmp_85 - tmp_60 * tmp_77 - tmp_85 * tmp_86 );
      real_t a_0_1 = tmp_104 * ( -tmp_103 * tmp_129 - tmp_124 * tmp_95 - tmp_126 * tmp_94 + tmp_129 * tmp_50 ) +
                     tmp_121 * ( -tmp_111 * tmp_126 - tmp_112 * tmp_124 - tmp_120 * tmp_130 + tmp_130 * tmp_50 ) +
                     tmp_52 * ( -tmp_123 * tmp_18 - tmp_124 * tmp_30 - tmp_125 * tmp_45 + tmp_125 * tmp_50 ) +
                     tmp_70 * ( -tmp_124 * tmp_61 - tmp_126 * tmp_59 + tmp_127 * tmp_50 - tmp_127 * tmp_69 ) +
                     tmp_87 * ( -tmp_124 * tmp_78 - tmp_126 * tmp_77 + tmp_128 * tmp_50 - tmp_128 * tmp_86 );
      real_t a_0_2 = tmp_104 * ( -tmp_103 * tmp_137 - tmp_132 * tmp_95 - tmp_134 * tmp_94 + tmp_137 * tmp_50 ) +
                     tmp_121 * ( -tmp_111 * tmp_134 - tmp_112 * tmp_132 - tmp_120 * tmp_138 + tmp_138 * tmp_50 ) +
                     tmp_52 * ( -tmp_131 * tmp_18 - tmp_132 * tmp_30 - tmp_133 * tmp_45 + tmp_133 * tmp_50 ) +
                     tmp_70 * ( -tmp_132 * tmp_61 - tmp_134 * tmp_59 + tmp_135 * tmp_50 - tmp_135 * tmp_69 ) +
                     tmp_87 * ( -tmp_132 * tmp_78 - tmp_134 * tmp_77 + tmp_136 * tmp_50 - tmp_136 * tmp_86 );
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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.95308992296933193 * p_affine_6_0 + 0.046910077030668018 * p_affine_7_0,
                                         0.95308992296933193 * p_affine_6_1 + 0.046910077030668018 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.7692346550528415 * p_affine_6_0 + 0.23076534494715845 * p_affine_7_0,
                                         0.7692346550528415 * p_affine_6_1 + 0.23076534494715845 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.5 * p_affine_6_0 + 0.5 * p_affine_7_0,
                                         0.5 * p_affine_6_1 + 0.5 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.2307653449471585 * p_affine_6_0 + 0.7692346550528415 * p_affine_7_0,
                                         0.2307653449471585 * p_affine_6_1 + 0.7692346550528415 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.046910077030668074 * p_affine_6_0 + 0.95308992296933193 * p_affine_7_0,
                                         0.046910077030668074 * p_affine_6_1 + 0.95308992296933193 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_1_1 + tmp_0;
      real_t tmp_2  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3  = -p_affine_0_0;
      real_t tmp_4  = p_affine_1_0 + tmp_3;
      real_t tmp_5  = p_affine_2_1 + tmp_0;
      real_t tmp_6  = p_affine_2_0 + tmp_3;
      real_t tmp_7  = 1.0 / ( -tmp_1 * tmp_6 + tmp_4 * tmp_5 );
      real_t tmp_8  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_9  = p_affine_6_1 + tmp_0;
      real_t tmp_10 = tmp_7 * ( 0.046910077030668018 * tmp_8 + tmp_9 );
      real_t tmp_11 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_12 = p_affine_6_0 + tmp_3;
      real_t tmp_13 = tmp_7 * ( 0.046910077030668018 * tmp_11 + tmp_12 );
      real_t tmp_14 = tmp_10 * tmp_2 + tmp_13 * tmp_5 - 1.0 / 3.0;
      real_t tmp_15 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_16 = tmp_10 * tmp_4 + tmp_13 * tmp_15 - 1.0 / 3.0;
      real_t tmp_17 = tmp_1 * tmp_14 + tmp_16 * tmp_5;
      real_t tmp_18 = 0.5 * tmp_7;
      real_t tmp_19 = tmp_18 * tmp_4;
      real_t tmp_20 = tmp_18 * tmp_2;
      real_t tmp_21 = -tmp_19 - tmp_20;
      real_t tmp_22 = p_affine_10_0 * tmp_21;
      real_t tmp_23 = 1.0 * tmp_7;
      real_t tmp_24 = tmp_23 * tmp_5;
      real_t tmp_25 = tmp_15 * tmp_23;
      real_t tmp_26 = p_affine_10_0 * ( -tmp_24 - tmp_25 ) + p_affine_10_1 * tmp_21;
      real_t tmp_27 = tmp_14 * tmp_4 + tmp_16 * tmp_6;
      real_t tmp_28 = 2 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_8 * tmp_8 ), 1.0 / 2.0 ) );
      real_t tmp_29 = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_28;
      real_t tmp_30 = tmp_7 * ( 0.23076534494715845 * tmp_8 + tmp_9 );
      real_t tmp_31 = tmp_7 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_32 = tmp_2 * tmp_30 + tmp_31 * tmp_5 - 1.0 / 3.0;
      real_t tmp_33 = tmp_15 * tmp_31 + tmp_30 * tmp_4 - 1.0 / 3.0;
      real_t tmp_34 = tmp_1 * tmp_32 + tmp_33 * tmp_5;
      real_t tmp_35 = tmp_32 * tmp_4 + tmp_33 * tmp_6;
      real_t tmp_36 = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_28;
      real_t tmp_37 = tmp_7 * ( 0.5 * tmp_8 + tmp_9 );
      real_t tmp_38 = tmp_7 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_39 = tmp_2 * tmp_37 + tmp_38 * tmp_5 - 1.0 / 3.0;
      real_t tmp_40 = tmp_15 * tmp_38 + tmp_37 * tmp_4 - 1.0 / 3.0;
      real_t tmp_41 = tmp_1 * tmp_39 + tmp_40 * tmp_5;
      real_t tmp_42 = tmp_39 * tmp_4 + tmp_40 * tmp_6;
      real_t tmp_43 = 0.2844444444444445 * Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_28;
      real_t tmp_44 = tmp_7 * ( 0.7692346550528415 * tmp_8 + tmp_9 );
      real_t tmp_45 = tmp_7 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_46 = tmp_2 * tmp_44 + tmp_45 * tmp_5 - 1.0 / 3.0;
      real_t tmp_47 = tmp_15 * tmp_45 + tmp_4 * tmp_44 - 1.0 / 3.0;
      real_t tmp_48 = tmp_1 * tmp_46 + tmp_47 * tmp_5;
      real_t tmp_49 = tmp_4 * tmp_46 + tmp_47 * tmp_6;
      real_t tmp_50 = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_28;
      real_t tmp_51 = tmp_7 * ( 0.95308992296933193 * tmp_8 + tmp_9 );
      real_t tmp_52 = tmp_7 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_53 = tmp_2 * tmp_51 + tmp_5 * tmp_52 - 1.0 / 3.0;
      real_t tmp_54 = tmp_15 * tmp_52 + tmp_4 * tmp_51 - 1.0 / 3.0;
      real_t tmp_55 = tmp_1 * tmp_53 + tmp_5 * tmp_54;
      real_t tmp_56 = tmp_4 * tmp_53 + tmp_54 * tmp_6;
      real_t tmp_57 = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_28;
      real_t tmp_58 = p_affine_10_0 * tmp_20;
      real_t tmp_59 = p_affine_10_0 * tmp_24 + p_affine_10_1 * tmp_20;
      real_t tmp_60 = p_affine_10_0 * tmp_19;
      real_t tmp_61 = p_affine_10_0 * tmp_25 + p_affine_10_1 * tmp_19;
      real_t a_0_0  = tmp_29 * ( -tmp_17 * tmp_22 - tmp_26 * tmp_27 ) + tmp_36 * ( -tmp_22 * tmp_34 - tmp_26 * tmp_35 ) +
                     tmp_43 * ( -tmp_22 * tmp_41 - tmp_26 * tmp_42 ) + tmp_50 * ( -tmp_22 * tmp_48 - tmp_26 * tmp_49 ) +
                     tmp_57 * ( -tmp_22 * tmp_55 - tmp_26 * tmp_56 );
      real_t a_0_1 = tmp_29 * ( -tmp_17 * tmp_58 - tmp_27 * tmp_59 ) + tmp_36 * ( -tmp_34 * tmp_58 - tmp_35 * tmp_59 ) +
                     tmp_43 * ( -tmp_41 * tmp_58 - tmp_42 * tmp_59 ) + tmp_50 * ( -tmp_48 * tmp_58 - tmp_49 * tmp_59 ) +
                     tmp_57 * ( -tmp_55 * tmp_58 - tmp_56 * tmp_59 );
      real_t a_0_2 = tmp_29 * ( -tmp_17 * tmp_60 - tmp_27 * tmp_61 ) + tmp_36 * ( -tmp_34 * tmp_60 - tmp_35 * tmp_61 ) +
                     tmp_43 * ( -tmp_41 * tmp_60 - tmp_42 * tmp_61 ) + tmp_50 * ( -tmp_48 * tmp_60 - tmp_49 * tmp_61 ) +
                     tmp_57 * ( -tmp_55 * tmp_60 - tmp_56 * tmp_61 );
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

 private:
   void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t* out_0 ) const
   {
      *out_0 = callback2D( Point3D( { in_0, in_1 } ) );
   }
   std::function< real_t( const Point3D& ) > callback2D;
};

class EGEpsilonFormP1E_0 : public hyteg::dg::DGForm2D
{
 public:
   EGEpsilonFormP1E_0( std::function< real_t( const Point3D& ) > _callback2D )
   : callback2D( _callback2D )
   {}

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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id5  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id6  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id7  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id8  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id9  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id11 = 0;
      Scalar_Variable_Coefficient_2D_mu(
          0.063089014491502282 * p_affine_0_0 + 0.063089014491502227 * p_affine_1_0 + 0.87382197101699555 * p_affine_2_0,
          0.063089014491502282 * p_affine_0_1 + 0.063089014491502227 * p_affine_1_1 + 0.87382197101699555 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu(
          0.24928674517091043 * p_affine_0_0 + 0.24928674517091043 * p_affine_1_0 + 0.50142650965817914 * p_affine_2_0,
          0.24928674517091043 * p_affine_0_1 + 0.24928674517091043 * p_affine_1_1 + 0.50142650965817914 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu(
          0.63650249912139867 * p_affine_0_0 + 0.31035245103378439 * p_affine_1_0 + 0.053145049844816938 * p_affine_2_0,
          0.63650249912139867 * p_affine_0_1 + 0.31035245103378439 * p_affine_1_1 + 0.053145049844816938 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu(
          0.053145049844816938 * p_affine_0_0 + 0.63650249912139867 * p_affine_1_0 + 0.31035245103378439 * p_affine_2_0,
          0.053145049844816938 * p_affine_0_1 + 0.63650249912139867 * p_affine_1_1 + 0.31035245103378439 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu(
          0.063089014491502227 * p_affine_0_0 + 0.87382197101699555 * p_affine_1_0 + 0.063089014491502227 * p_affine_2_0,
          0.063089014491502227 * p_affine_0_1 + 0.87382197101699555 * p_affine_1_1 + 0.063089014491502227 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      Scalar_Variable_Coefficient_2D_mu(
          0.24928674517091043 * p_affine_0_0 + 0.50142650965817914 * p_affine_1_0 + 0.24928674517091043 * p_affine_2_0,
          0.24928674517091043 * p_affine_0_1 + 0.50142650965817914 * p_affine_1_1 + 0.24928674517091043 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id5 );
      Scalar_Variable_Coefficient_2D_mu(
          0.87382197101699566 * p_affine_0_0 + 0.063089014491502227 * p_affine_1_0 + 0.063089014491502227 * p_affine_2_0,
          0.87382197101699566 * p_affine_0_1 + 0.063089014491502227 * p_affine_1_1 + 0.063089014491502227 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id6 );
      Scalar_Variable_Coefficient_2D_mu(
          0.50142650965817914 * p_affine_0_0 + 0.24928674517091043 * p_affine_1_0 + 0.24928674517091043 * p_affine_2_0,
          0.50142650965817914 * p_affine_0_1 + 0.24928674517091043 * p_affine_1_1 + 0.24928674517091043 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id7 );
      Scalar_Variable_Coefficient_2D_mu(
          0.053145049844816938 * p_affine_0_0 + 0.31035245103378439 * p_affine_1_0 + 0.63650249912139867 * p_affine_2_0,
          0.053145049844816938 * p_affine_0_1 + 0.31035245103378439 * p_affine_1_1 + 0.63650249912139867 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id8 );
      Scalar_Variable_Coefficient_2D_mu(
          0.63650249912139867 * p_affine_0_0 + 0.053145049844816938 * p_affine_1_0 + 0.31035245103378439 * p_affine_2_0,
          0.63650249912139867 * p_affine_0_1 + 0.053145049844816938 * p_affine_1_1 + 0.31035245103378439 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id9 );
      Scalar_Variable_Coefficient_2D_mu(
          0.31035245103378439 * p_affine_0_0 + 0.63650249912139867 * p_affine_1_0 + 0.053145049844816938 * p_affine_2_0,
          0.31035245103378439 * p_affine_0_1 + 0.63650249912139867 * p_affine_1_1 + 0.053145049844816938 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id10 );
      Scalar_Variable_Coefficient_2D_mu(
          0.31035245103378439 * p_affine_0_0 + 0.053145049844816938 * p_affine_1_0 + 0.63650249912139867 * p_affine_2_0,
          0.31035245103378439 * p_affine_0_1 + 0.053145049844816938 * p_affine_1_1 + 0.63650249912139867 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id11 );
      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_0;
      real_t tmp_6  = p_affine_1_1 + tmp_2;
      real_t tmp_7  = 1.0 / ( tmp_4 - tmp_5 * tmp_6 );
      real_t tmp_8  = 1.0 * tmp_7;
      real_t tmp_9  = p_affine_0_1 - p_affine_1_1;
      real_t tmp_10 = tmp_8 * tmp_9;
      real_t tmp_11 = tmp_10 * tmp_5 + tmp_4 * tmp_8;
      real_t tmp_12 = Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_11;
      real_t tmp_13 = -2 * tmp_10 - 2 * tmp_3 * tmp_8;
      real_t tmp_14 = 0.5 * tmp_7;
      real_t tmp_15 = tmp_1 * tmp_14;
      real_t tmp_16 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_17 = tmp_14 * tmp_16;
      real_t tmp_18 = tmp_14 * tmp_3;
      real_t tmp_19 = tmp_1 * tmp_17 + tmp_15 * tmp_5 + tmp_18 * tmp_6 + tmp_18 * tmp_9;
      real_t tmp_20 = Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_19;
      real_t tmp_21 = -4 * tmp_15 - 4 * tmp_17;
      real_t tmp_22 = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                                p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_23 = 0.025422453185103409 * tmp_22;
      real_t tmp_24 = tmp_11 * tmp_13;
      real_t tmp_25 = tmp_19 * tmp_21;
      real_t tmp_26 = 0.058393137863189684 * tmp_22;
      real_t tmp_27 = 0.041425537809186785 * tmp_22;
      real_t tmp_28 = 0.041425537809186785 * tmp_22;
      real_t tmp_29 = 0.025422453185103409 * tmp_22;
      real_t tmp_30 = 0.058393137863189684 * tmp_22;
      real_t tmp_31 = 0.025422453185103409 * tmp_22;
      real_t tmp_32 = 0.058393137863189684 * tmp_22;
      real_t tmp_33 = 0.041425537809186785 * tmp_22;
      real_t tmp_34 = 0.041425537809186785 * tmp_22;
      real_t tmp_35 = 0.041425537809186785 * tmp_22;
      real_t tmp_36 = 0.041425537809186785 * tmp_22;
      real_t tmp_37 = 2.0 * tmp_7;
      real_t tmp_38 = tmp_3 * tmp_37;
      real_t tmp_39 = tmp_16 * tmp_37;
      real_t tmp_40 = tmp_11 * tmp_38;
      real_t tmp_41 = tmp_19 * tmp_39;
      real_t tmp_42 = tmp_37 * tmp_9;
      real_t tmp_43 = tmp_1 * tmp_37;
      real_t tmp_44 = tmp_11 * tmp_42;
      real_t tmp_45 = tmp_19 * tmp_43;
      real_t a_0_0 =
          tmp_23 * ( tmp_12 * tmp_13 + tmp_20 * tmp_21 ) +
          tmp_26 * ( Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_25 ) +
          tmp_27 * ( Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_25 ) +
          tmp_28 * ( Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_25 ) +
          tmp_29 * ( Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_25 ) +
          tmp_30 * ( Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_25 ) +
          tmp_31 * ( Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_25 ) +
          tmp_32 * ( Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_25 ) +
          tmp_33 * ( Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_25 ) +
          tmp_34 * ( Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_25 ) +
          tmp_35 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_25 ) +
          tmp_36 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_25 );
      real_t a_1_0 =
          tmp_23 * ( tmp_12 * tmp_38 + tmp_20 * tmp_39 ) +
          tmp_26 * ( Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_41 ) +
          tmp_27 * ( Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_41 ) +
          tmp_28 * ( Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_41 ) +
          tmp_29 * ( Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_41 ) +
          tmp_30 * ( Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_41 ) +
          tmp_31 * ( Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_41 ) +
          tmp_32 * ( Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_41 ) +
          tmp_33 * ( Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_41 ) +
          tmp_34 * ( Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_41 ) +
          tmp_35 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_41 ) +
          tmp_36 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_41 );
      real_t a_2_0 =
          tmp_23 * ( tmp_12 * tmp_42 + tmp_20 * tmp_43 ) +
          tmp_26 * ( Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_45 ) +
          tmp_27 * ( Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_45 ) +
          tmp_28 * ( Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_45 ) +
          tmp_29 * ( Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_45 ) +
          tmp_30 * ( Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_45 ) +
          tmp_31 * ( Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_45 ) +
          tmp_32 * ( Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_45 ) +
          tmp_33 * ( Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_45 ) +
          tmp_34 * ( Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_45 ) +
          tmp_35 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_45 ) +
          tmp_36 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_45 );
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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.95308992296933193 * p_affine_6_0 + 0.046910077030668018 * p_affine_7_0,
                                         0.95308992296933193 * p_affine_6_1 + 0.046910077030668018 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.7692346550528415 * p_affine_6_0 + 0.23076534494715845 * p_affine_7_0,
                                         0.7692346550528415 * p_affine_6_1 + 0.23076534494715845 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.5 * p_affine_6_0 + 0.5 * p_affine_7_0,
                                         0.5 * p_affine_6_1 + 0.5 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.2307653449471585 * p_affine_6_0 + 0.7692346550528415 * p_affine_7_0,
                                         0.2307653449471585 * p_affine_6_1 + 0.7692346550528415 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.046910077030668074 * p_affine_6_0 + 0.95308992296933193 * p_affine_7_0,
                                         0.046910077030668074 * p_affine_6_1 + 0.95308992296933193 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      real_t tmp_0  = -p_affine_0_1;
      real_t tmp_1  = p_affine_1_1 + tmp_0;
      real_t tmp_2  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3  = -p_affine_0_0;
      real_t tmp_4  = p_affine_1_0 + tmp_3;
      real_t tmp_5  = p_affine_2_1 + tmp_0;
      real_t tmp_6  = tmp_4 * tmp_5;
      real_t tmp_7  = p_affine_2_0 + tmp_3;
      real_t tmp_8  = 1.0 / ( -tmp_1 * tmp_7 + tmp_6 );
      real_t tmp_9  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + tmp_0;
      real_t tmp_11 = tmp_8 * ( tmp_10 + 0.046910077030668018 * tmp_9 );
      real_t tmp_12 = tmp_11 * tmp_2;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_3;
      real_t tmp_15 = tmp_8 * ( 0.046910077030668018 * tmp_13 + tmp_14 );
      real_t tmp_16 = tmp_15 * tmp_5;
      real_t tmp_17 = tmp_12 + tmp_16;
      real_t tmp_18 = tmp_17 - 1.0 / 3.0;
      real_t tmp_19 = tmp_11 * tmp_4;
      real_t tmp_20 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_21 = tmp_15 * tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_22 - 1.0 / 3.0;
      real_t tmp_24 = p_affine_10_0 * ( tmp_1 * tmp_18 + tmp_23 * tmp_5 );
      real_t tmp_25 = 0.5 * tmp_8;
      real_t tmp_26 = tmp_25 * tmp_4;
      real_t tmp_27 = tmp_2 * tmp_25;
      real_t tmp_28 = -tmp_26 - tmp_27;
      real_t tmp_29 = 0.5 * tmp_28;
      real_t tmp_30 = tmp_18 * tmp_4 + tmp_23 * tmp_7;
      real_t tmp_31 = 1.0 * tmp_8;
      real_t tmp_32 = tmp_31 * tmp_5;
      real_t tmp_33 = tmp_20 * tmp_31;
      real_t tmp_34 = 0.5 * p_affine_10_0 * ( -tmp_32 - tmp_33 ) + 0.5 * p_affine_10_1 * tmp_28;
      real_t tmp_35 = -tmp_12 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_36 = std::abs( std::pow( ( tmp_13 * tmp_13 ) + ( tmp_9 * tmp_9 ), 1.0 / 2.0 ) );
      real_t tmp_37 = 6 / tmp_36;
      real_t tmp_38 = tmp_30 * tmp_37;
      real_t tmp_39 = tmp_25 * tmp_5;
      real_t tmp_40 = 0.5 * p_affine_10_0 * ( tmp_31 * tmp_6 + tmp_33 * tmp_7 ) +
                      0.5 * p_affine_10_1 * ( tmp_1 * tmp_39 + tmp_20 * tmp_39 + tmp_26 * tmp_7 + tmp_27 * tmp_4 );
      real_t tmp_41  = 2 * tmp_36;
      real_t tmp_42  = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_41;
      real_t tmp_43  = tmp_8 * ( tmp_10 + 0.23076534494715845 * tmp_9 );
      real_t tmp_44  = tmp_2 * tmp_43;
      real_t tmp_45  = tmp_8 * ( 0.23076534494715845 * tmp_13 + tmp_14 );
      real_t tmp_46  = tmp_45 * tmp_5;
      real_t tmp_47  = tmp_44 + tmp_46;
      real_t tmp_48  = tmp_47 - 1.0 / 3.0;
      real_t tmp_49  = tmp_4 * tmp_43;
      real_t tmp_50  = tmp_20 * tmp_45;
      real_t tmp_51  = tmp_49 + tmp_50;
      real_t tmp_52  = tmp_51 - 1.0 / 3.0;
      real_t tmp_53  = tmp_1 * tmp_48 + tmp_5 * tmp_52;
      real_t tmp_54  = p_affine_10_0 * tmp_29;
      real_t tmp_55  = tmp_4 * tmp_48 + tmp_52 * tmp_7;
      real_t tmp_56  = -tmp_44 - tmp_46 - tmp_49 - tmp_50 + 1;
      real_t tmp_57  = tmp_37 * tmp_55;
      real_t tmp_58  = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_41;
      real_t tmp_59  = tmp_8 * ( tmp_10 + 0.5 * tmp_9 );
      real_t tmp_60  = tmp_2 * tmp_59;
      real_t tmp_61  = tmp_8 * ( 0.5 * tmp_13 + tmp_14 );
      real_t tmp_62  = tmp_5 * tmp_61;
      real_t tmp_63  = tmp_60 + tmp_62;
      real_t tmp_64  = tmp_63 - 1.0 / 3.0;
      real_t tmp_65  = tmp_4 * tmp_59;
      real_t tmp_66  = tmp_20 * tmp_61;
      real_t tmp_67  = tmp_65 + tmp_66;
      real_t tmp_68  = tmp_67 - 1.0 / 3.0;
      real_t tmp_69  = tmp_1 * tmp_64 + tmp_5 * tmp_68;
      real_t tmp_70  = tmp_4 * tmp_64 + tmp_68 * tmp_7;
      real_t tmp_71  = -tmp_60 - tmp_62 - tmp_65 - tmp_66 + 1;
      real_t tmp_72  = tmp_37 * tmp_70;
      real_t tmp_73  = 0.2844444444444445 * Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_41;
      real_t tmp_74  = tmp_8 * ( tmp_10 + 0.7692346550528415 * tmp_9 );
      real_t tmp_75  = tmp_2 * tmp_74;
      real_t tmp_76  = tmp_8 * ( 0.7692346550528415 * tmp_13 + tmp_14 );
      real_t tmp_77  = tmp_5 * tmp_76;
      real_t tmp_78  = tmp_75 + tmp_77;
      real_t tmp_79  = tmp_78 - 1.0 / 3.0;
      real_t tmp_80  = tmp_4 * tmp_74;
      real_t tmp_81  = tmp_20 * tmp_76;
      real_t tmp_82  = tmp_80 + tmp_81;
      real_t tmp_83  = tmp_82 - 1.0 / 3.0;
      real_t tmp_84  = tmp_1 * tmp_79 + tmp_5 * tmp_83;
      real_t tmp_85  = tmp_4 * tmp_79 + tmp_7 * tmp_83;
      real_t tmp_86  = -tmp_75 - tmp_77 - tmp_80 - tmp_81 + 1;
      real_t tmp_87  = tmp_37 * tmp_85;
      real_t tmp_88  = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_41;
      real_t tmp_89  = tmp_8 * ( tmp_10 + 0.95308992296933193 * tmp_9 );
      real_t tmp_90  = tmp_2 * tmp_89;
      real_t tmp_91  = tmp_8 * ( 0.95308992296933193 * tmp_13 + tmp_14 );
      real_t tmp_92  = tmp_5 * tmp_91;
      real_t tmp_93  = tmp_90 + tmp_92;
      real_t tmp_94  = tmp_93 - 1.0 / 3.0;
      real_t tmp_95  = tmp_4 * tmp_89;
      real_t tmp_96  = tmp_20 * tmp_91;
      real_t tmp_97  = tmp_95 + tmp_96;
      real_t tmp_98  = tmp_97 - 1.0 / 3.0;
      real_t tmp_99  = tmp_1 * tmp_94 + tmp_5 * tmp_98;
      real_t tmp_100 = tmp_4 * tmp_94 + tmp_7 * tmp_98;
      real_t tmp_101 = -tmp_90 - tmp_92 - tmp_95 - tmp_96 + 1;
      real_t tmp_102 = tmp_100 * tmp_37;
      real_t tmp_103 = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_41;
      real_t tmp_104 = 0.25 * tmp_8;
      real_t tmp_105 = tmp_104 * tmp_2;
      real_t tmp_106 = 0.5 * p_affine_10_0 * tmp_32 + 0.5 * p_affine_10_1 * tmp_27;
      real_t tmp_107 = p_affine_10_0 * tmp_105;
      real_t tmp_108 = tmp_104 * tmp_4;
      real_t tmp_109 = 0.5 * p_affine_10_0 * tmp_33 + 0.5 * p_affine_10_1 * tmp_26;
      real_t tmp_110 = p_affine_10_0 * tmp_108;
      real_t a_0_0   = tmp_103 * ( -tmp_100 * tmp_34 + tmp_101 * tmp_102 - tmp_101 * tmp_40 - tmp_54 * tmp_99 ) +
                     tmp_42 * ( -tmp_24 * tmp_29 - tmp_30 * tmp_34 + tmp_35 * tmp_38 - tmp_35 * tmp_40 ) +
                     tmp_58 * ( -tmp_34 * tmp_55 - tmp_40 * tmp_56 - tmp_53 * tmp_54 + tmp_56 * tmp_57 ) +
                     tmp_73 * ( -tmp_34 * tmp_70 - tmp_40 * tmp_71 - tmp_54 * tmp_69 + tmp_71 * tmp_72 ) +
                     tmp_88 * ( -tmp_34 * tmp_85 - tmp_40 * tmp_86 - tmp_54 * tmp_84 + tmp_86 * tmp_87 );
      real_t a_1_0 = tmp_103 * ( -tmp_100 * tmp_106 + tmp_102 * tmp_93 - tmp_107 * tmp_99 - tmp_40 * tmp_93 ) +
                     tmp_42 * ( -tmp_105 * tmp_24 - tmp_106 * tmp_30 + tmp_17 * tmp_38 - tmp_17 * tmp_40 ) +
                     tmp_58 * ( -tmp_106 * tmp_55 - tmp_107 * tmp_53 - tmp_40 * tmp_47 + tmp_47 * tmp_57 ) +
                     tmp_73 * ( -tmp_106 * tmp_70 - tmp_107 * tmp_69 - tmp_40 * tmp_63 + tmp_63 * tmp_72 ) +
                     tmp_88 * ( -tmp_106 * tmp_85 - tmp_107 * tmp_84 - tmp_40 * tmp_78 + tmp_78 * tmp_87 );
      real_t a_2_0 = tmp_103 * ( -tmp_100 * tmp_109 + tmp_102 * tmp_97 - tmp_110 * tmp_99 - tmp_40 * tmp_97 ) +
                     tmp_42 * ( -tmp_108 * tmp_24 - tmp_109 * tmp_30 + tmp_22 * tmp_38 - tmp_22 * tmp_40 ) +
                     tmp_58 * ( -tmp_109 * tmp_55 - tmp_110 * tmp_53 - tmp_40 * tmp_51 + tmp_51 * tmp_57 ) +
                     tmp_73 * ( -tmp_109 * tmp_70 - tmp_110 * tmp_69 - tmp_40 * tmp_67 + tmp_67 * tmp_72 ) +
                     tmp_88 * ( -tmp_109 * tmp_85 - tmp_110 * tmp_84 - tmp_40 * tmp_82 + tmp_82 * tmp_87 );
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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.95308992296933193 * p_affine_6_0 + 0.046910077030668018 * p_affine_7_0,
                                         0.95308992296933193 * p_affine_6_1 + 0.046910077030668018 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.7692346550528415 * p_affine_6_0 + 0.23076534494715845 * p_affine_7_0,
                                         0.7692346550528415 * p_affine_6_1 + 0.23076534494715845 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.5 * p_affine_6_0 + 0.5 * p_affine_7_0,
                                         0.5 * p_affine_6_1 + 0.5 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.2307653449471585 * p_affine_6_0 + 0.7692346550528415 * p_affine_7_0,
                                         0.2307653449471585 * p_affine_6_1 + 0.7692346550528415 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.046910077030668074 * p_affine_6_0 + 0.95308992296933193 * p_affine_7_0,
                                         0.046910077030668074 * p_affine_6_1 + 0.95308992296933193 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      real_t tmp_0  = -p_affine_3_1;
      real_t tmp_1  = p_affine_4_1 + tmp_0;
      real_t tmp_2  = p_affine_3_0 - p_affine_5_0;
      real_t tmp_3  = -p_affine_3_0;
      real_t tmp_4  = p_affine_4_0 + tmp_3;
      real_t tmp_5  = p_affine_5_1 + tmp_0;
      real_t tmp_6  = tmp_4 * tmp_5;
      real_t tmp_7  = p_affine_5_0 + tmp_3;
      real_t tmp_8  = 1.0 / ( -tmp_1 * tmp_7 + tmp_6 );
      real_t tmp_9  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + 0.046910077030668018 * tmp_9;
      real_t tmp_11 = tmp_8 * ( tmp_0 + tmp_10 );
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + 0.046910077030668018 * tmp_12;
      real_t tmp_14 = tmp_8 * ( tmp_13 + tmp_3 );
      real_t tmp_15 = tmp_11 * tmp_2 + tmp_14 * tmp_5 - 1.0 / 3.0;
      real_t tmp_16 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_17 = tmp_11 * tmp_4 + tmp_14 * tmp_16 - 1.0 / 3.0;
      real_t tmp_18 = p_affine_10_0 * ( tmp_1 * tmp_15 + tmp_17 * tmp_5 );
      real_t tmp_19 = -p_affine_0_0;
      real_t tmp_20 = p_affine_1_0 + tmp_19;
      real_t tmp_21 = -p_affine_0_1;
      real_t tmp_22 = p_affine_2_1 + tmp_21;
      real_t tmp_23 = 1.0 / ( tmp_20 * tmp_22 - ( p_affine_1_1 + tmp_21 ) * ( p_affine_2_0 + tmp_19 ) );
      real_t tmp_24 = 0.5 * tmp_23;
      real_t tmp_25 = tmp_20 * tmp_24;
      real_t tmp_26 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_27 = tmp_24 * tmp_26;
      real_t tmp_28 = -tmp_25 - tmp_27;
      real_t tmp_29 = 0.5 * tmp_28;
      real_t tmp_30 = tmp_15 * tmp_4 + tmp_17 * tmp_7;
      real_t tmp_31 = 1.0 * tmp_23;
      real_t tmp_32 = tmp_22 * tmp_31;
      real_t tmp_33 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_34 = tmp_31 * tmp_33;
      real_t tmp_35 = 0.5 * p_affine_10_0 * ( -tmp_32 - tmp_34 ) + 0.5 * p_affine_10_1 * tmp_28;
      real_t tmp_36 = tmp_23 * ( tmp_10 + tmp_21 );
      real_t tmp_37 = tmp_20 * tmp_36;
      real_t tmp_38 = tmp_26 * tmp_36;
      real_t tmp_39 = tmp_23 * ( tmp_13 + tmp_19 );
      real_t tmp_40 = tmp_22 * tmp_39;
      real_t tmp_41 = tmp_33 * tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = std::abs( std::pow( ( tmp_12 * tmp_12 ) + ( tmp_9 * tmp_9 ), 1.0 / 2.0 ) );
      real_t tmp_44 = 6 / tmp_43;
      real_t tmp_45 = tmp_30 * tmp_44;
      real_t tmp_46 = 1.0 * tmp_8;
      real_t tmp_47 = 0.5 * tmp_8;
      real_t tmp_48 = tmp_4 * tmp_47;
      real_t tmp_49 = tmp_47 * tmp_5;
      real_t tmp_50 = 0.5 * p_affine_10_0 * ( tmp_16 * tmp_46 * tmp_7 + tmp_46 * tmp_6 ) +
                      0.5 * p_affine_10_1 * ( tmp_1 * tmp_49 + tmp_16 * tmp_49 + tmp_2 * tmp_48 + tmp_48 * tmp_7 );
      real_t tmp_51  = 2 * tmp_43;
      real_t tmp_52  = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_51;
      real_t tmp_53  = p_affine_6_1 + 0.23076534494715845 * tmp_9;
      real_t tmp_54  = tmp_8 * ( tmp_0 + tmp_53 );
      real_t tmp_55  = p_affine_6_0 + 0.23076534494715845 * tmp_12;
      real_t tmp_56  = tmp_8 * ( tmp_3 + tmp_55 );
      real_t tmp_57  = tmp_2 * tmp_54 + tmp_5 * tmp_56 - 1.0 / 3.0;
      real_t tmp_58  = tmp_16 * tmp_56 + tmp_4 * tmp_54 - 1.0 / 3.0;
      real_t tmp_59  = tmp_1 * tmp_57 + tmp_5 * tmp_58;
      real_t tmp_60  = p_affine_10_0 * tmp_29;
      real_t tmp_61  = tmp_4 * tmp_57 + tmp_58 * tmp_7;
      real_t tmp_62  = tmp_23 * ( tmp_21 + tmp_53 );
      real_t tmp_63  = tmp_20 * tmp_62;
      real_t tmp_64  = tmp_26 * tmp_62;
      real_t tmp_65  = tmp_23 * ( tmp_19 + tmp_55 );
      real_t tmp_66  = tmp_22 * tmp_65;
      real_t tmp_67  = tmp_33 * tmp_65;
      real_t tmp_68  = -tmp_63 - tmp_64 - tmp_66 - tmp_67 + 1;
      real_t tmp_69  = tmp_44 * tmp_61;
      real_t tmp_70  = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_51;
      real_t tmp_71  = p_affine_6_1 + 0.5 * tmp_9;
      real_t tmp_72  = tmp_8 * ( tmp_0 + tmp_71 );
      real_t tmp_73  = p_affine_6_0 + 0.5 * tmp_12;
      real_t tmp_74  = tmp_8 * ( tmp_3 + tmp_73 );
      real_t tmp_75  = tmp_2 * tmp_72 + tmp_5 * tmp_74 - 1.0 / 3.0;
      real_t tmp_76  = tmp_16 * tmp_74 + tmp_4 * tmp_72 - 1.0 / 3.0;
      real_t tmp_77  = tmp_1 * tmp_75 + tmp_5 * tmp_76;
      real_t tmp_78  = tmp_4 * tmp_75 + tmp_7 * tmp_76;
      real_t tmp_79  = tmp_23 * ( tmp_21 + tmp_71 );
      real_t tmp_80  = tmp_20 * tmp_79;
      real_t tmp_81  = tmp_26 * tmp_79;
      real_t tmp_82  = tmp_23 * ( tmp_19 + tmp_73 );
      real_t tmp_83  = tmp_22 * tmp_82;
      real_t tmp_84  = tmp_33 * tmp_82;
      real_t tmp_85  = -tmp_80 - tmp_81 - tmp_83 - tmp_84 + 1;
      real_t tmp_86  = tmp_44 * tmp_78;
      real_t tmp_87  = 0.2844444444444445 * Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_51;
      real_t tmp_88  = p_affine_6_1 + 0.7692346550528415 * tmp_9;
      real_t tmp_89  = tmp_8 * ( tmp_0 + tmp_88 );
      real_t tmp_90  = p_affine_6_0 + 0.7692346550528415 * tmp_12;
      real_t tmp_91  = tmp_8 * ( tmp_3 + tmp_90 );
      real_t tmp_92  = tmp_2 * tmp_89 + tmp_5 * tmp_91 - 1.0 / 3.0;
      real_t tmp_93  = tmp_16 * tmp_91 + tmp_4 * tmp_89 - 1.0 / 3.0;
      real_t tmp_94  = tmp_1 * tmp_92 + tmp_5 * tmp_93;
      real_t tmp_95  = tmp_4 * tmp_92 + tmp_7 * tmp_93;
      real_t tmp_96  = tmp_23 * ( tmp_21 + tmp_88 );
      real_t tmp_97  = tmp_20 * tmp_96;
      real_t tmp_98  = tmp_26 * tmp_96;
      real_t tmp_99  = tmp_23 * ( tmp_19 + tmp_90 );
      real_t tmp_100 = tmp_22 * tmp_99;
      real_t tmp_101 = tmp_33 * tmp_99;
      real_t tmp_102 = -tmp_100 - tmp_101 - tmp_97 - tmp_98 + 1;
      real_t tmp_103 = tmp_44 * tmp_95;
      real_t tmp_104 = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_51;
      real_t tmp_105 = p_affine_6_1 + 0.95308992296933193 * tmp_9;
      real_t tmp_106 = tmp_8 * ( tmp_0 + tmp_105 );
      real_t tmp_107 = p_affine_6_0 + 0.95308992296933193 * tmp_12;
      real_t tmp_108 = tmp_8 * ( tmp_107 + tmp_3 );
      real_t tmp_109 = tmp_106 * tmp_2 + tmp_108 * tmp_5 - 1.0 / 3.0;
      real_t tmp_110 = tmp_106 * tmp_4 + tmp_108 * tmp_16 - 1.0 / 3.0;
      real_t tmp_111 = tmp_1 * tmp_109 + tmp_110 * tmp_5;
      real_t tmp_112 = tmp_109 * tmp_4 + tmp_110 * tmp_7;
      real_t tmp_113 = tmp_23 * ( tmp_105 + tmp_21 );
      real_t tmp_114 = tmp_113 * tmp_20;
      real_t tmp_115 = tmp_113 * tmp_26;
      real_t tmp_116 = tmp_23 * ( tmp_107 + tmp_19 );
      real_t tmp_117 = tmp_116 * tmp_22;
      real_t tmp_118 = tmp_116 * tmp_33;
      real_t tmp_119 = -tmp_114 - tmp_115 - tmp_117 - tmp_118 + 1;
      real_t tmp_120 = tmp_112 * tmp_44;
      real_t tmp_121 = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_51;
      real_t tmp_122 = 0.25 * tmp_23;
      real_t tmp_123 = tmp_122 * tmp_26;
      real_t tmp_124 = 0.5 * p_affine_10_0 * tmp_32 + 0.5 * p_affine_10_1 * tmp_27;
      real_t tmp_125 = tmp_38 + tmp_40;
      real_t tmp_126 = p_affine_10_0 * tmp_123;
      real_t tmp_127 = tmp_64 + tmp_66;
      real_t tmp_128 = tmp_81 + tmp_83;
      real_t tmp_129 = tmp_100 + tmp_98;
      real_t tmp_130 = tmp_115 + tmp_117;
      real_t tmp_131 = tmp_122 * tmp_20;
      real_t tmp_132 = 0.5 * p_affine_10_0 * tmp_34 + 0.5 * p_affine_10_1 * tmp_25;
      real_t tmp_133 = tmp_37 + tmp_41;
      real_t tmp_134 = p_affine_10_0 * tmp_131;
      real_t tmp_135 = tmp_63 + tmp_67;
      real_t tmp_136 = tmp_80 + tmp_84;
      real_t tmp_137 = tmp_101 + tmp_97;
      real_t tmp_138 = tmp_114 + tmp_118;
      real_t a_0_0   = tmp_104 * ( -tmp_102 * tmp_103 - tmp_102 * tmp_50 + tmp_35 * tmp_95 + tmp_60 * tmp_94 ) +
                     tmp_121 * ( tmp_111 * tmp_60 + tmp_112 * tmp_35 - tmp_119 * tmp_120 - tmp_119 * tmp_50 ) +
                     tmp_52 * ( tmp_18 * tmp_29 + tmp_30 * tmp_35 - tmp_42 * tmp_45 - tmp_42 * tmp_50 ) +
                     tmp_70 * ( tmp_35 * tmp_61 - tmp_50 * tmp_68 + tmp_59 * tmp_60 - tmp_68 * tmp_69 ) +
                     tmp_87 * ( tmp_35 * tmp_78 - tmp_50 * tmp_85 + tmp_60 * tmp_77 - tmp_85 * tmp_86 );
      real_t a_1_0 = tmp_104 * ( -tmp_103 * tmp_129 + tmp_124 * tmp_95 + tmp_126 * tmp_94 - tmp_129 * tmp_50 ) +
                     tmp_121 * ( tmp_111 * tmp_126 + tmp_112 * tmp_124 - tmp_120 * tmp_130 - tmp_130 * tmp_50 ) +
                     tmp_52 * ( tmp_123 * tmp_18 + tmp_124 * tmp_30 - tmp_125 * tmp_45 - tmp_125 * tmp_50 ) +
                     tmp_70 * ( tmp_124 * tmp_61 + tmp_126 * tmp_59 - tmp_127 * tmp_50 - tmp_127 * tmp_69 ) +
                     tmp_87 * ( tmp_124 * tmp_78 + tmp_126 * tmp_77 - tmp_128 * tmp_50 - tmp_128 * tmp_86 );
      real_t a_2_0 = tmp_104 * ( -tmp_103 * tmp_137 + tmp_132 * tmp_95 + tmp_134 * tmp_94 - tmp_137 * tmp_50 ) +
                     tmp_121 * ( tmp_111 * tmp_134 + tmp_112 * tmp_132 - tmp_120 * tmp_138 - tmp_138 * tmp_50 ) +
                     tmp_52 * ( tmp_131 * tmp_18 + tmp_132 * tmp_30 - tmp_133 * tmp_45 - tmp_133 * tmp_50 ) +
                     tmp_70 * ( tmp_132 * tmp_61 + tmp_134 * tmp_59 - tmp_135 * tmp_50 - tmp_135 * tmp_69 ) +
                     tmp_87 * ( tmp_132 * tmp_78 + tmp_134 * tmp_77 - tmp_136 * tmp_50 - tmp_136 * tmp_86 );
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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.95308992296933193 * p_affine_6_0 + 0.046910077030668018 * p_affine_7_0,
                                         0.95308992296933193 * p_affine_6_1 + 0.046910077030668018 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.7692346550528415 * p_affine_6_0 + 0.23076534494715845 * p_affine_7_0,
                                         0.7692346550528415 * p_affine_6_1 + 0.23076534494715845 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.5 * p_affine_6_0 + 0.5 * p_affine_7_0,
                                         0.5 * p_affine_6_1 + 0.5 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.2307653449471585 * p_affine_6_0 + 0.7692346550528415 * p_affine_7_0,
                                         0.2307653449471585 * p_affine_6_1 + 0.7692346550528415 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.046910077030668074 * p_affine_6_0 + 0.95308992296933193 * p_affine_7_0,
                                         0.046910077030668074 * p_affine_6_1 + 0.95308992296933193 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
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
      real_t tmp_20 = 1.4215613103371365 * Scalar_Variable_Coefficient_2D_mu_out0_id0 *
                      ( tmp_1 * ( tmp_18 - 1.0 / 3.0 ) + tmp_4 * ( tmp_19 - 1.0 / 3.0 ) );
      real_t tmp_21 = tmp_5 * ( 0.23076534494715845 * tmp_6 + tmp_7 );
      real_t tmp_22 = tmp_1 * tmp_21;
      real_t tmp_23 = tmp_10 * tmp_21;
      real_t tmp_24 = tmp_5 * ( 0.23076534494715845 * tmp_12 + tmp_13 );
      real_t tmp_25 = tmp_24 * tmp_3;
      real_t tmp_26 = tmp_16 * tmp_24;
      real_t tmp_27 = tmp_23 + tmp_25;
      real_t tmp_28 = tmp_22 + tmp_26;
      real_t tmp_29 = 2.8717720229961969 * Scalar_Variable_Coefficient_2D_mu_out0_id1 *
                      ( tmp_1 * ( tmp_27 - 1.0 / 3.0 ) + tmp_4 * ( tmp_28 - 1.0 / 3.0 ) );
      real_t tmp_30 = tmp_5 * ( 0.5 * tmp_6 + tmp_7 );
      real_t tmp_31 = tmp_1 * tmp_30;
      real_t tmp_32 = tmp_10 * tmp_30;
      real_t tmp_33 = tmp_5 * ( 0.5 * tmp_12 + tmp_13 );
      real_t tmp_34 = tmp_3 * tmp_33;
      real_t tmp_35 = tmp_16 * tmp_33;
      real_t tmp_36 = tmp_32 + tmp_34;
      real_t tmp_37 = tmp_31 + tmp_35;
      real_t tmp_38 = 3.413333333333334 * Scalar_Variable_Coefficient_2D_mu_out0_id2 *
                      ( tmp_1 * ( tmp_36 - 1.0 / 3.0 ) + tmp_4 * ( tmp_37 - 1.0 / 3.0 ) );
      real_t tmp_39 = tmp_5 * ( 0.7692346550528415 * tmp_6 + tmp_7 );
      real_t tmp_40 = tmp_1 * tmp_39;
      real_t tmp_41 = tmp_10 * tmp_39;
      real_t tmp_42 = tmp_5 * ( 0.7692346550528415 * tmp_12 + tmp_13 );
      real_t tmp_43 = tmp_3 * tmp_42;
      real_t tmp_44 = tmp_16 * tmp_42;
      real_t tmp_45 = tmp_41 + tmp_43;
      real_t tmp_46 = tmp_40 + tmp_44;
      real_t tmp_47 = 2.8717720229961969 * Scalar_Variable_Coefficient_2D_mu_out0_id3 *
                      ( tmp_1 * ( tmp_45 - 1.0 / 3.0 ) + tmp_4 * ( tmp_46 - 1.0 / 3.0 ) );
      real_t tmp_48 = tmp_5 * ( 0.95308992296933193 * tmp_6 + tmp_7 );
      real_t tmp_49 = tmp_1 * tmp_48;
      real_t tmp_50 = tmp_10 * tmp_48;
      real_t tmp_51 = tmp_5 * ( 0.95308992296933193 * tmp_12 + tmp_13 );
      real_t tmp_52 = tmp_3 * tmp_51;
      real_t tmp_53 = tmp_16 * tmp_51;
      real_t tmp_54 = tmp_50 + tmp_52;
      real_t tmp_55 = tmp_49 + tmp_53;
      real_t tmp_56 = 1.4215613103371365 * Scalar_Variable_Coefficient_2D_mu_out0_id4 *
                      ( tmp_1 * ( tmp_54 - 1.0 / 3.0 ) + tmp_4 * ( tmp_55 - 1.0 / 3.0 ) );
      real_t a_0_0 = tmp_20 * ( -tmp_11 - tmp_15 - tmp_17 - tmp_9 + 1 ) + tmp_29 * ( -tmp_22 - tmp_23 - tmp_25 - tmp_26 + 1 ) +
                     tmp_38 * ( -tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1 ) + tmp_47 * ( -tmp_40 - tmp_41 - tmp_43 - tmp_44 + 1 ) +
                     tmp_56 * ( -tmp_49 - tmp_50 - tmp_52 - tmp_53 + 1 );
      real_t a_1_0  = tmp_18 * tmp_20 + tmp_27 * tmp_29 + tmp_36 * tmp_38 + tmp_45 * tmp_47 + tmp_54 * tmp_56;
      real_t a_2_0  = tmp_19 * tmp_20 + tmp_28 * tmp_29 + tmp_37 * tmp_38 + tmp_46 * tmp_47 + tmp_55 * tmp_56;
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

 private:
   void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t* out_0 ) const
   {
      *out_0 = callback2D( Point3D( { in_0, in_1 } ) );
   }
   std::function< real_t( const Point3D& ) > callback2D;
};

class EGEpsilonFormEP1_1 : public hyteg::dg::DGForm2D
{
 public:
   EGEpsilonFormEP1_1( std::function< real_t( const Point3D& ) > _callback2D )
   : callback2D( _callback2D )
   {}
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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id5  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id6  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id7  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id8  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id9  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id11 = 0;
      Scalar_Variable_Coefficient_2D_mu(
          0.063089014491502282 * p_affine_0_0 + 0.063089014491502227 * p_affine_1_0 + 0.87382197101699555 * p_affine_2_0,
          0.063089014491502282 * p_affine_0_1 + 0.063089014491502227 * p_affine_1_1 + 0.87382197101699555 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu(
          0.24928674517091043 * p_affine_0_0 + 0.24928674517091043 * p_affine_1_0 + 0.50142650965817914 * p_affine_2_0,
          0.24928674517091043 * p_affine_0_1 + 0.24928674517091043 * p_affine_1_1 + 0.50142650965817914 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu(
          0.63650249912139867 * p_affine_0_0 + 0.31035245103378439 * p_affine_1_0 + 0.053145049844816938 * p_affine_2_0,
          0.63650249912139867 * p_affine_0_1 + 0.31035245103378439 * p_affine_1_1 + 0.053145049844816938 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu(
          0.053145049844816938 * p_affine_0_0 + 0.63650249912139867 * p_affine_1_0 + 0.31035245103378439 * p_affine_2_0,
          0.053145049844816938 * p_affine_0_1 + 0.63650249912139867 * p_affine_1_1 + 0.31035245103378439 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu(
          0.063089014491502227 * p_affine_0_0 + 0.87382197101699555 * p_affine_1_0 + 0.063089014491502227 * p_affine_2_0,
          0.063089014491502227 * p_affine_0_1 + 0.87382197101699555 * p_affine_1_1 + 0.063089014491502227 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      Scalar_Variable_Coefficient_2D_mu(
          0.24928674517091043 * p_affine_0_0 + 0.50142650965817914 * p_affine_1_0 + 0.24928674517091043 * p_affine_2_0,
          0.24928674517091043 * p_affine_0_1 + 0.50142650965817914 * p_affine_1_1 + 0.24928674517091043 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id5 );
      Scalar_Variable_Coefficient_2D_mu(
          0.87382197101699566 * p_affine_0_0 + 0.063089014491502227 * p_affine_1_0 + 0.063089014491502227 * p_affine_2_0,
          0.87382197101699566 * p_affine_0_1 + 0.063089014491502227 * p_affine_1_1 + 0.063089014491502227 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id6 );
      Scalar_Variable_Coefficient_2D_mu(
          0.50142650965817914 * p_affine_0_0 + 0.24928674517091043 * p_affine_1_0 + 0.24928674517091043 * p_affine_2_0,
          0.50142650965817914 * p_affine_0_1 + 0.24928674517091043 * p_affine_1_1 + 0.24928674517091043 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id7 );
      Scalar_Variable_Coefficient_2D_mu(
          0.053145049844816938 * p_affine_0_0 + 0.31035245103378439 * p_affine_1_0 + 0.63650249912139867 * p_affine_2_0,
          0.053145049844816938 * p_affine_0_1 + 0.31035245103378439 * p_affine_1_1 + 0.63650249912139867 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id8 );
      Scalar_Variable_Coefficient_2D_mu(
          0.63650249912139867 * p_affine_0_0 + 0.053145049844816938 * p_affine_1_0 + 0.31035245103378439 * p_affine_2_0,
          0.63650249912139867 * p_affine_0_1 + 0.053145049844816938 * p_affine_1_1 + 0.31035245103378439 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id9 );
      Scalar_Variable_Coefficient_2D_mu(
          0.31035245103378439 * p_affine_0_0 + 0.63650249912139867 * p_affine_1_0 + 0.053145049844816938 * p_affine_2_0,
          0.31035245103378439 * p_affine_0_1 + 0.63650249912139867 * p_affine_1_1 + 0.053145049844816938 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id10 );
      Scalar_Variable_Coefficient_2D_mu(
          0.31035245103378439 * p_affine_0_0 + 0.053145049844816938 * p_affine_1_0 + 0.63650249912139867 * p_affine_2_0,
          0.31035245103378439 * p_affine_0_1 + 0.053145049844816938 * p_affine_1_1 + 0.63650249912139867 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id11 );
      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_0;
      real_t tmp_6  = p_affine_1_1 + tmp_2;
      real_t tmp_7  = 1.0 / ( tmp_4 - tmp_5 * tmp_6 );
      real_t tmp_8  = 1.0 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_8 * tmp_9;
      real_t tmp_11 = tmp_10 * tmp_6 + tmp_4 * tmp_8;
      real_t tmp_12 = Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_11;
      real_t tmp_13 = -2 * tmp_1 * tmp_8 - 2 * tmp_10;
      real_t tmp_14 = 0.5 * tmp_7;
      real_t tmp_15 = tmp_1 * tmp_14;
      real_t tmp_16 = tmp_14 * tmp_3;
      real_t tmp_17 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_18 = tmp_14 * tmp_17;
      real_t tmp_19 = tmp_15 * tmp_5 + tmp_15 * tmp_9 + tmp_16 * tmp_6 + tmp_18 * tmp_3;
      real_t tmp_20 = Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_19;
      real_t tmp_21 = -4 * tmp_16 - 4 * tmp_18;
      real_t tmp_22 = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                                p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_23 = 0.025422453185103409 * tmp_22;
      real_t tmp_24 = tmp_11 * tmp_13;
      real_t tmp_25 = tmp_19 * tmp_21;
      real_t tmp_26 = 0.058393137863189684 * tmp_22;
      real_t tmp_27 = 0.041425537809186785 * tmp_22;
      real_t tmp_28 = 0.041425537809186785 * tmp_22;
      real_t tmp_29 = 0.025422453185103409 * tmp_22;
      real_t tmp_30 = 0.058393137863189684 * tmp_22;
      real_t tmp_31 = 0.025422453185103409 * tmp_22;
      real_t tmp_32 = 0.058393137863189684 * tmp_22;
      real_t tmp_33 = 0.041425537809186785 * tmp_22;
      real_t tmp_34 = 0.041425537809186785 * tmp_22;
      real_t tmp_35 = 0.041425537809186785 * tmp_22;
      real_t tmp_36 = 0.041425537809186785 * tmp_22;
      real_t tmp_37 = 2.0 * tmp_7;
      real_t tmp_38 = tmp_37 * tmp_9;
      real_t tmp_39 = tmp_3 * tmp_37;
      real_t tmp_40 = tmp_11 * tmp_38;
      real_t tmp_41 = tmp_19 * tmp_39;
      real_t tmp_42 = tmp_1 * tmp_37;
      real_t tmp_43 = tmp_17 * tmp_37;
      real_t tmp_44 = tmp_11 * tmp_42;
      real_t tmp_45 = tmp_19 * tmp_43;
      real_t a_0_0 =
          tmp_23 * ( tmp_12 * tmp_13 + tmp_20 * tmp_21 ) +
          tmp_26 * ( Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_25 ) +
          tmp_27 * ( Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_25 ) +
          tmp_28 * ( Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_25 ) +
          tmp_29 * ( Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_25 ) +
          tmp_30 * ( Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_25 ) +
          tmp_31 * ( Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_25 ) +
          tmp_32 * ( Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_25 ) +
          tmp_33 * ( Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_25 ) +
          tmp_34 * ( Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_25 ) +
          tmp_35 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_25 ) +
          tmp_36 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_25 );
      real_t a_0_1 =
          tmp_23 * ( tmp_12 * tmp_38 + tmp_20 * tmp_39 ) +
          tmp_26 * ( Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_41 ) +
          tmp_27 * ( Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_41 ) +
          tmp_28 * ( Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_41 ) +
          tmp_29 * ( Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_41 ) +
          tmp_30 * ( Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_41 ) +
          tmp_31 * ( Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_41 ) +
          tmp_32 * ( Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_41 ) +
          tmp_33 * ( Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_41 ) +
          tmp_34 * ( Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_41 ) +
          tmp_35 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_41 ) +
          tmp_36 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_41 );
      real_t a_0_2 =
          tmp_23 * ( tmp_12 * tmp_42 + tmp_20 * tmp_43 ) +
          tmp_26 * ( Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_45 ) +
          tmp_27 * ( Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_45 ) +
          tmp_28 * ( Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_45 ) +
          tmp_29 * ( Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_45 ) +
          tmp_30 * ( Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_45 ) +
          tmp_31 * ( Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_45 ) +
          tmp_32 * ( Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_45 ) +
          tmp_33 * ( Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_45 ) +
          tmp_34 * ( Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_45 ) +
          tmp_35 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_45 ) +
          tmp_36 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_45 );
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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.95308992296933193 * p_affine_6_0 + 0.046910077030668018 * p_affine_7_0,
                                         0.95308992296933193 * p_affine_6_1 + 0.046910077030668018 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.7692346550528415 * p_affine_6_0 + 0.23076534494715845 * p_affine_7_0,
                                         0.7692346550528415 * p_affine_6_1 + 0.23076534494715845 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.5 * p_affine_6_0 + 0.5 * p_affine_7_0,
                                         0.5 * p_affine_6_1 + 0.5 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.2307653449471585 * p_affine_6_0 + 0.7692346550528415 * p_affine_7_0,
                                         0.2307653449471585 * p_affine_6_1 + 0.7692346550528415 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.046910077030668074 * p_affine_6_0 + 0.95308992296933193 * p_affine_7_0,
                                         0.046910077030668074 * p_affine_6_1 + 0.95308992296933193 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3  = -p_affine_0_1;
      real_t tmp_4  = p_affine_2_1 + tmp_3;
      real_t tmp_5  = tmp_1 * tmp_4;
      real_t tmp_6  = p_affine_2_0 + tmp_0;
      real_t tmp_7  = p_affine_1_1 + tmp_3;
      real_t tmp_8  = 1.0 / ( tmp_5 - tmp_6 * tmp_7 );
      real_t tmp_9  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + tmp_3;
      real_t tmp_11 = tmp_8 * ( tmp_10 + 0.046910077030668018 * tmp_9 );
      real_t tmp_12 = tmp_11 * tmp_2;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_0;
      real_t tmp_15 = tmp_8 * ( 0.046910077030668018 * tmp_13 + tmp_14 );
      real_t tmp_16 = tmp_15 * tmp_4;
      real_t tmp_17 = tmp_12 + tmp_16;
      real_t tmp_18 = tmp_17 - 1.0 / 3.0;
      real_t tmp_19 = tmp_1 * tmp_11;
      real_t tmp_20 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_21 = tmp_15 * tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_22 - 1.0 / 3.0;
      real_t tmp_24 = p_affine_10_1 * ( tmp_1 * tmp_18 + tmp_23 * tmp_6 );
      real_t tmp_25 = 0.5 * tmp_8;
      real_t tmp_26 = tmp_25 * tmp_4;
      real_t tmp_27 = tmp_20 * tmp_25;
      real_t tmp_28 = -tmp_26 - tmp_27;
      real_t tmp_29 = 0.5 * tmp_28;
      real_t tmp_30 = tmp_18 * tmp_7 + tmp_23 * tmp_4;
      real_t tmp_31 = 1.0 * tmp_8;
      real_t tmp_32 = tmp_1 * tmp_31;
      real_t tmp_33 = tmp_2 * tmp_31;
      real_t tmp_34 = 0.5 * p_affine_10_0 * tmp_28 + 0.5 * p_affine_10_1 * ( -tmp_32 - tmp_33 );
      real_t tmp_35 = -tmp_12 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_36 = std::abs( std::pow( ( tmp_13 * tmp_13 ) + ( tmp_9 * tmp_9 ), 1.0 / 2.0 ) );
      real_t tmp_37 = 6 / tmp_36;
      real_t tmp_38 = tmp_30 * tmp_37;
      real_t tmp_39 = tmp_1 * tmp_25;
      real_t tmp_40 = 0.5 * p_affine_10_0 * ( tmp_2 * tmp_39 + tmp_26 * tmp_7 + tmp_27 * tmp_4 + tmp_39 * tmp_6 ) +
                      0.5 * p_affine_10_1 * ( tmp_31 * tmp_5 + tmp_33 * tmp_7 );
      real_t tmp_41  = 2 * tmp_36;
      real_t tmp_42  = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_41;
      real_t tmp_43  = tmp_8 * ( tmp_10 + 0.23076534494715845 * tmp_9 );
      real_t tmp_44  = tmp_2 * tmp_43;
      real_t tmp_45  = tmp_8 * ( 0.23076534494715845 * tmp_13 + tmp_14 );
      real_t tmp_46  = tmp_4 * tmp_45;
      real_t tmp_47  = tmp_44 + tmp_46;
      real_t tmp_48  = tmp_47 - 1.0 / 3.0;
      real_t tmp_49  = tmp_1 * tmp_43;
      real_t tmp_50  = tmp_20 * tmp_45;
      real_t tmp_51  = tmp_49 + tmp_50;
      real_t tmp_52  = tmp_51 - 1.0 / 3.0;
      real_t tmp_53  = tmp_1 * tmp_48 + tmp_52 * tmp_6;
      real_t tmp_54  = p_affine_10_1 * tmp_29;
      real_t tmp_55  = tmp_4 * tmp_52 + tmp_48 * tmp_7;
      real_t tmp_56  = -tmp_44 - tmp_46 - tmp_49 - tmp_50 + 1;
      real_t tmp_57  = tmp_37 * tmp_55;
      real_t tmp_58  = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_41;
      real_t tmp_59  = tmp_8 * ( tmp_10 + 0.5 * tmp_9 );
      real_t tmp_60  = tmp_2 * tmp_59;
      real_t tmp_61  = tmp_8 * ( 0.5 * tmp_13 + tmp_14 );
      real_t tmp_62  = tmp_4 * tmp_61;
      real_t tmp_63  = tmp_60 + tmp_62;
      real_t tmp_64  = tmp_63 - 1.0 / 3.0;
      real_t tmp_65  = tmp_1 * tmp_59;
      real_t tmp_66  = tmp_20 * tmp_61;
      real_t tmp_67  = tmp_65 + tmp_66;
      real_t tmp_68  = tmp_67 - 1.0 / 3.0;
      real_t tmp_69  = tmp_1 * tmp_64 + tmp_6 * tmp_68;
      real_t tmp_70  = tmp_4 * tmp_68 + tmp_64 * tmp_7;
      real_t tmp_71  = -tmp_60 - tmp_62 - tmp_65 - tmp_66 + 1;
      real_t tmp_72  = tmp_37 * tmp_70;
      real_t tmp_73  = 0.2844444444444445 * Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_41;
      real_t tmp_74  = tmp_8 * ( tmp_10 + 0.7692346550528415 * tmp_9 );
      real_t tmp_75  = tmp_2 * tmp_74;
      real_t tmp_76  = tmp_8 * ( 0.7692346550528415 * tmp_13 + tmp_14 );
      real_t tmp_77  = tmp_4 * tmp_76;
      real_t tmp_78  = tmp_75 + tmp_77;
      real_t tmp_79  = tmp_78 - 1.0 / 3.0;
      real_t tmp_80  = tmp_1 * tmp_74;
      real_t tmp_81  = tmp_20 * tmp_76;
      real_t tmp_82  = tmp_80 + tmp_81;
      real_t tmp_83  = tmp_82 - 1.0 / 3.0;
      real_t tmp_84  = tmp_1 * tmp_79 + tmp_6 * tmp_83;
      real_t tmp_85  = tmp_4 * tmp_83 + tmp_7 * tmp_79;
      real_t tmp_86  = -tmp_75 - tmp_77 - tmp_80 - tmp_81 + 1;
      real_t tmp_87  = tmp_37 * tmp_85;
      real_t tmp_88  = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_41;
      real_t tmp_89  = tmp_8 * ( tmp_10 + 0.95308992296933193 * tmp_9 );
      real_t tmp_90  = tmp_2 * tmp_89;
      real_t tmp_91  = tmp_8 * ( 0.95308992296933193 * tmp_13 + tmp_14 );
      real_t tmp_92  = tmp_4 * tmp_91;
      real_t tmp_93  = tmp_90 + tmp_92;
      real_t tmp_94  = tmp_93 - 1.0 / 3.0;
      real_t tmp_95  = tmp_1 * tmp_89;
      real_t tmp_96  = tmp_20 * tmp_91;
      real_t tmp_97  = tmp_95 + tmp_96;
      real_t tmp_98  = tmp_97 - 1.0 / 3.0;
      real_t tmp_99  = tmp_1 * tmp_94 + tmp_6 * tmp_98;
      real_t tmp_100 = tmp_4 * tmp_98 + tmp_7 * tmp_94;
      real_t tmp_101 = -tmp_90 - tmp_92 - tmp_95 - tmp_96 + 1;
      real_t tmp_102 = tmp_100 * tmp_37;
      real_t tmp_103 = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_41;
      real_t tmp_104 = 0.25 * tmp_8;
      real_t tmp_105 = tmp_104 * tmp_4;
      real_t tmp_106 = 0.5 * p_affine_10_0 * tmp_26 + 0.5 * p_affine_10_1 * tmp_33;
      real_t tmp_107 = p_affine_10_1 * tmp_105;
      real_t tmp_108 = tmp_104 * tmp_20;
      real_t tmp_109 = 0.5 * p_affine_10_0 * tmp_27 + 0.5 * p_affine_10_1 * tmp_32;
      real_t tmp_110 = p_affine_10_1 * tmp_108;
      real_t a_0_0   = tmp_103 * ( -tmp_100 * tmp_34 + tmp_101 * tmp_102 - tmp_101 * tmp_40 - tmp_54 * tmp_99 ) +
                     tmp_42 * ( -tmp_24 * tmp_29 - tmp_30 * tmp_34 + tmp_35 * tmp_38 - tmp_35 * tmp_40 ) +
                     tmp_58 * ( -tmp_34 * tmp_55 - tmp_40 * tmp_56 - tmp_53 * tmp_54 + tmp_56 * tmp_57 ) +
                     tmp_73 * ( -tmp_34 * tmp_70 - tmp_40 * tmp_71 - tmp_54 * tmp_69 + tmp_71 * tmp_72 ) +
                     tmp_88 * ( -tmp_34 * tmp_85 - tmp_40 * tmp_86 - tmp_54 * tmp_84 + tmp_86 * tmp_87 );
      real_t a_0_1 = tmp_103 * ( -tmp_100 * tmp_106 + tmp_102 * tmp_93 - tmp_107 * tmp_99 - tmp_40 * tmp_93 ) +
                     tmp_42 * ( -tmp_105 * tmp_24 - tmp_106 * tmp_30 + tmp_17 * tmp_38 - tmp_17 * tmp_40 ) +
                     tmp_58 * ( -tmp_106 * tmp_55 - tmp_107 * tmp_53 - tmp_40 * tmp_47 + tmp_47 * tmp_57 ) +
                     tmp_73 * ( -tmp_106 * tmp_70 - tmp_107 * tmp_69 - tmp_40 * tmp_63 + tmp_63 * tmp_72 ) +
                     tmp_88 * ( -tmp_106 * tmp_85 - tmp_107 * tmp_84 - tmp_40 * tmp_78 + tmp_78 * tmp_87 );
      real_t a_0_2 = tmp_103 * ( -tmp_100 * tmp_109 + tmp_102 * tmp_97 - tmp_110 * tmp_99 - tmp_40 * tmp_97 ) +
                     tmp_42 * ( -tmp_108 * tmp_24 - tmp_109 * tmp_30 + tmp_22 * tmp_38 - tmp_22 * tmp_40 ) +
                     tmp_58 * ( -tmp_109 * tmp_55 - tmp_110 * tmp_53 - tmp_40 * tmp_51 + tmp_51 * tmp_57 ) +
                     tmp_73 * ( -tmp_109 * tmp_70 - tmp_110 * tmp_69 - tmp_40 * tmp_67 + tmp_67 * tmp_72 ) +
                     tmp_88 * ( -tmp_109 * tmp_85 - tmp_110 * tmp_84 - tmp_40 * tmp_82 + tmp_82 * tmp_87 );
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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.95308992296933193 * p_affine_6_0 + 0.046910077030668018 * p_affine_7_0,
                                         0.95308992296933193 * p_affine_6_1 + 0.046910077030668018 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.7692346550528415 * p_affine_6_0 + 0.23076534494715845 * p_affine_7_0,
                                         0.7692346550528415 * p_affine_6_1 + 0.23076534494715845 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.5 * p_affine_6_0 + 0.5 * p_affine_7_0,
                                         0.5 * p_affine_6_1 + 0.5 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.2307653449471585 * p_affine_6_0 + 0.7692346550528415 * p_affine_7_0,
                                         0.2307653449471585 * p_affine_6_1 + 0.7692346550528415 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.046910077030668074 * p_affine_6_0 + 0.95308992296933193 * p_affine_7_0,
                                         0.046910077030668074 * p_affine_6_1 + 0.95308992296933193 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3  = -p_affine_0_1;
      real_t tmp_4  = p_affine_2_1 + tmp_3;
      real_t tmp_5  = tmp_1 * tmp_4;
      real_t tmp_6  = p_affine_2_0 + tmp_0;
      real_t tmp_7  = p_affine_1_1 + tmp_3;
      real_t tmp_8  = 1.0 / ( tmp_5 - tmp_6 * tmp_7 );
      real_t tmp_9  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + 0.046910077030668018 * tmp_9;
      real_t tmp_11 = tmp_8 * ( tmp_10 + tmp_3 );
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + 0.046910077030668018 * tmp_12;
      real_t tmp_14 = tmp_8 * ( tmp_0 + tmp_13 );
      real_t tmp_15 = tmp_11 * tmp_2 + tmp_14 * tmp_4 - 1.0 / 3.0;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_1 * tmp_11 + tmp_14 * tmp_16 - 1.0 / 3.0;
      real_t tmp_18 = p_affine_10_1 * ( tmp_1 * tmp_15 + tmp_17 * tmp_6 );
      real_t tmp_19 = -p_affine_3_1;
      real_t tmp_20 = p_affine_5_1 + tmp_19;
      real_t tmp_21 = -p_affine_3_0;
      real_t tmp_22 = p_affine_4_0 + tmp_21;
      real_t tmp_23 = 1.0 / ( tmp_20 * tmp_22 - ( p_affine_4_1 + tmp_19 ) * ( p_affine_5_0 + tmp_21 ) );
      real_t tmp_24 = 0.5 * tmp_23;
      real_t tmp_25 = tmp_20 * tmp_24;
      real_t tmp_26 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_27 = tmp_24 * tmp_26;
      real_t tmp_28 = -tmp_25 - tmp_27;
      real_t tmp_29 = 0.5 * tmp_28;
      real_t tmp_30 = tmp_15 * tmp_7 + tmp_17 * tmp_4;
      real_t tmp_31 = 1.0 * tmp_23;
      real_t tmp_32 = tmp_22 * tmp_31;
      real_t tmp_33 = p_affine_3_0 - p_affine_5_0;
      real_t tmp_34 = tmp_31 * tmp_33;
      real_t tmp_35 = 0.5 * p_affine_10_0 * tmp_28 + 0.5 * p_affine_10_1 * ( -tmp_32 - tmp_34 );
      real_t tmp_36 = tmp_23 * ( tmp_10 + tmp_19 );
      real_t tmp_37 = tmp_22 * tmp_36;
      real_t tmp_38 = tmp_33 * tmp_36;
      real_t tmp_39 = tmp_23 * ( tmp_13 + tmp_21 );
      real_t tmp_40 = tmp_20 * tmp_39;
      real_t tmp_41 = tmp_26 * tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = std::abs( std::pow( ( tmp_12 * tmp_12 ) + ( tmp_9 * tmp_9 ), 1.0 / 2.0 ) );
      real_t tmp_44 = 6 / tmp_43;
      real_t tmp_45 = tmp_30 * tmp_44;
      real_t tmp_46 = 1.0 * tmp_8;
      real_t tmp_47 = 0.5 * tmp_8;
      real_t tmp_48 = tmp_1 * tmp_47;
      real_t tmp_49 = tmp_4 * tmp_47;
      real_t tmp_50 = 0.5 * p_affine_10_0 * ( tmp_16 * tmp_49 + tmp_2 * tmp_48 + tmp_48 * tmp_6 + tmp_49 * tmp_7 ) +
                      0.5 * p_affine_10_1 * ( tmp_2 * tmp_46 * tmp_7 + tmp_46 * tmp_5 );
      real_t tmp_51  = 2 * tmp_43;
      real_t tmp_52  = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_51;
      real_t tmp_53  = p_affine_6_1 + 0.23076534494715845 * tmp_9;
      real_t tmp_54  = tmp_8 * ( tmp_3 + tmp_53 );
      real_t tmp_55  = p_affine_6_0 + 0.23076534494715845 * tmp_12;
      real_t tmp_56  = tmp_8 * ( tmp_0 + tmp_55 );
      real_t tmp_57  = tmp_2 * tmp_54 + tmp_4 * tmp_56 - 1.0 / 3.0;
      real_t tmp_58  = tmp_1 * tmp_54 + tmp_16 * tmp_56 - 1.0 / 3.0;
      real_t tmp_59  = tmp_1 * tmp_57 + tmp_58 * tmp_6;
      real_t tmp_60  = p_affine_10_1 * tmp_29;
      real_t tmp_61  = tmp_4 * tmp_58 + tmp_57 * tmp_7;
      real_t tmp_62  = tmp_23 * ( tmp_19 + tmp_53 );
      real_t tmp_63  = tmp_22 * tmp_62;
      real_t tmp_64  = tmp_33 * tmp_62;
      real_t tmp_65  = tmp_23 * ( tmp_21 + tmp_55 );
      real_t tmp_66  = tmp_20 * tmp_65;
      real_t tmp_67  = tmp_26 * tmp_65;
      real_t tmp_68  = -tmp_63 - tmp_64 - tmp_66 - tmp_67 + 1;
      real_t tmp_69  = tmp_44 * tmp_61;
      real_t tmp_70  = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_51;
      real_t tmp_71  = p_affine_6_1 + 0.5 * tmp_9;
      real_t tmp_72  = tmp_8 * ( tmp_3 + tmp_71 );
      real_t tmp_73  = p_affine_6_0 + 0.5 * tmp_12;
      real_t tmp_74  = tmp_8 * ( tmp_0 + tmp_73 );
      real_t tmp_75  = tmp_2 * tmp_72 + tmp_4 * tmp_74 - 1.0 / 3.0;
      real_t tmp_76  = tmp_1 * tmp_72 + tmp_16 * tmp_74 - 1.0 / 3.0;
      real_t tmp_77  = tmp_1 * tmp_75 + tmp_6 * tmp_76;
      real_t tmp_78  = tmp_4 * tmp_76 + tmp_7 * tmp_75;
      real_t tmp_79  = tmp_23 * ( tmp_19 + tmp_71 );
      real_t tmp_80  = tmp_22 * tmp_79;
      real_t tmp_81  = tmp_33 * tmp_79;
      real_t tmp_82  = tmp_23 * ( tmp_21 + tmp_73 );
      real_t tmp_83  = tmp_20 * tmp_82;
      real_t tmp_84  = tmp_26 * tmp_82;
      real_t tmp_85  = -tmp_80 - tmp_81 - tmp_83 - tmp_84 + 1;
      real_t tmp_86  = tmp_44 * tmp_78;
      real_t tmp_87  = 0.2844444444444445 * Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_51;
      real_t tmp_88  = p_affine_6_1 + 0.7692346550528415 * tmp_9;
      real_t tmp_89  = tmp_8 * ( tmp_3 + tmp_88 );
      real_t tmp_90  = p_affine_6_0 + 0.7692346550528415 * tmp_12;
      real_t tmp_91  = tmp_8 * ( tmp_0 + tmp_90 );
      real_t tmp_92  = tmp_2 * tmp_89 + tmp_4 * tmp_91 - 1.0 / 3.0;
      real_t tmp_93  = tmp_1 * tmp_89 + tmp_16 * tmp_91 - 1.0 / 3.0;
      real_t tmp_94  = tmp_1 * tmp_92 + tmp_6 * tmp_93;
      real_t tmp_95  = tmp_4 * tmp_93 + tmp_7 * tmp_92;
      real_t tmp_96  = tmp_23 * ( tmp_19 + tmp_88 );
      real_t tmp_97  = tmp_22 * tmp_96;
      real_t tmp_98  = tmp_33 * tmp_96;
      real_t tmp_99  = tmp_23 * ( tmp_21 + tmp_90 );
      real_t tmp_100 = tmp_20 * tmp_99;
      real_t tmp_101 = tmp_26 * tmp_99;
      real_t tmp_102 = -tmp_100 - tmp_101 - tmp_97 - tmp_98 + 1;
      real_t tmp_103 = tmp_44 * tmp_95;
      real_t tmp_104 = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_51;
      real_t tmp_105 = p_affine_6_1 + 0.95308992296933193 * tmp_9;
      real_t tmp_106 = tmp_8 * ( tmp_105 + tmp_3 );
      real_t tmp_107 = p_affine_6_0 + 0.95308992296933193 * tmp_12;
      real_t tmp_108 = tmp_8 * ( tmp_0 + tmp_107 );
      real_t tmp_109 = tmp_106 * tmp_2 + tmp_108 * tmp_4 - 1.0 / 3.0;
      real_t tmp_110 = tmp_1 * tmp_106 + tmp_108 * tmp_16 - 1.0 / 3.0;
      real_t tmp_111 = tmp_1 * tmp_109 + tmp_110 * tmp_6;
      real_t tmp_112 = tmp_109 * tmp_7 + tmp_110 * tmp_4;
      real_t tmp_113 = tmp_23 * ( tmp_105 + tmp_19 );
      real_t tmp_114 = tmp_113 * tmp_22;
      real_t tmp_115 = tmp_113 * tmp_33;
      real_t tmp_116 = tmp_23 * ( tmp_107 + tmp_21 );
      real_t tmp_117 = tmp_116 * tmp_20;
      real_t tmp_118 = tmp_116 * tmp_26;
      real_t tmp_119 = -tmp_114 - tmp_115 - tmp_117 - tmp_118 + 1;
      real_t tmp_120 = tmp_112 * tmp_44;
      real_t tmp_121 = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_51;
      real_t tmp_122 = 0.25 * tmp_23;
      real_t tmp_123 = tmp_122 * tmp_20;
      real_t tmp_124 = 0.5 * p_affine_10_0 * tmp_25 + 0.5 * p_affine_10_1 * tmp_34;
      real_t tmp_125 = tmp_38 + tmp_40;
      real_t tmp_126 = p_affine_10_1 * tmp_123;
      real_t tmp_127 = tmp_64 + tmp_66;
      real_t tmp_128 = tmp_81 + tmp_83;
      real_t tmp_129 = tmp_100 + tmp_98;
      real_t tmp_130 = tmp_115 + tmp_117;
      real_t tmp_131 = tmp_122 * tmp_26;
      real_t tmp_132 = 0.5 * p_affine_10_0 * tmp_27 + 0.5 * p_affine_10_1 * tmp_32;
      real_t tmp_133 = tmp_37 + tmp_41;
      real_t tmp_134 = p_affine_10_1 * tmp_131;
      real_t tmp_135 = tmp_63 + tmp_67;
      real_t tmp_136 = tmp_80 + tmp_84;
      real_t tmp_137 = tmp_101 + tmp_97;
      real_t tmp_138 = tmp_114 + tmp_118;
      real_t a_0_0   = tmp_104 * ( -tmp_102 * tmp_103 + tmp_102 * tmp_50 - tmp_35 * tmp_95 - tmp_60 * tmp_94 ) +
                     tmp_121 * ( -tmp_111 * tmp_60 - tmp_112 * tmp_35 - tmp_119 * tmp_120 + tmp_119 * tmp_50 ) +
                     tmp_52 * ( -tmp_18 * tmp_29 - tmp_30 * tmp_35 - tmp_42 * tmp_45 + tmp_42 * tmp_50 ) +
                     tmp_70 * ( -tmp_35 * tmp_61 + tmp_50 * tmp_68 - tmp_59 * tmp_60 - tmp_68 * tmp_69 ) +
                     tmp_87 * ( -tmp_35 * tmp_78 + tmp_50 * tmp_85 - tmp_60 * tmp_77 - tmp_85 * tmp_86 );
      real_t a_0_1 = tmp_104 * ( -tmp_103 * tmp_129 - tmp_124 * tmp_95 - tmp_126 * tmp_94 + tmp_129 * tmp_50 ) +
                     tmp_121 * ( -tmp_111 * tmp_126 - tmp_112 * tmp_124 - tmp_120 * tmp_130 + tmp_130 * tmp_50 ) +
                     tmp_52 * ( -tmp_123 * tmp_18 - tmp_124 * tmp_30 - tmp_125 * tmp_45 + tmp_125 * tmp_50 ) +
                     tmp_70 * ( -tmp_124 * tmp_61 - tmp_126 * tmp_59 + tmp_127 * tmp_50 - tmp_127 * tmp_69 ) +
                     tmp_87 * ( -tmp_124 * tmp_78 - tmp_126 * tmp_77 + tmp_128 * tmp_50 - tmp_128 * tmp_86 );
      real_t a_0_2 = tmp_104 * ( -tmp_103 * tmp_137 - tmp_132 * tmp_95 - tmp_134 * tmp_94 + tmp_137 * tmp_50 ) +
                     tmp_121 * ( -tmp_111 * tmp_134 - tmp_112 * tmp_132 - tmp_120 * tmp_138 + tmp_138 * tmp_50 ) +
                     tmp_52 * ( -tmp_131 * tmp_18 - tmp_132 * tmp_30 - tmp_133 * tmp_45 + tmp_133 * tmp_50 ) +
                     tmp_70 * ( -tmp_132 * tmp_61 - tmp_134 * tmp_59 + tmp_135 * tmp_50 - tmp_135 * tmp_69 ) +
                     tmp_87 * ( -tmp_132 * tmp_78 - tmp_134 * tmp_77 + tmp_136 * tmp_50 - tmp_136 * tmp_86 );
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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.95308992296933193 * p_affine_6_0 + 0.046910077030668018 * p_affine_7_0,
                                         0.95308992296933193 * p_affine_6_1 + 0.046910077030668018 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.7692346550528415 * p_affine_6_0 + 0.23076534494715845 * p_affine_7_0,
                                         0.7692346550528415 * p_affine_6_1 + 0.23076534494715845 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.5 * p_affine_6_0 + 0.5 * p_affine_7_0,
                                         0.5 * p_affine_6_1 + 0.5 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.2307653449471585 * p_affine_6_0 + 0.7692346550528415 * p_affine_7_0,
                                         0.2307653449471585 * p_affine_6_1 + 0.7692346550528415 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.046910077030668074 * p_affine_6_0 + 0.95308992296933193 * p_affine_7_0,
                                         0.046910077030668074 * p_affine_6_1 + 0.95308992296933193 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
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
      real_t tmp_17 = tmp_1 * tmp_14 + tmp_16 * tmp_5;
      real_t tmp_18 = 0.5 * tmp_7;
      real_t tmp_19 = tmp_18 * tmp_4;
      real_t tmp_20 = tmp_15 * tmp_18;
      real_t tmp_21 = -tmp_19 - tmp_20;
      real_t tmp_22 = p_affine_10_1 * tmp_21;
      real_t tmp_23 = 1.0 * tmp_7;
      real_t tmp_24 = tmp_1 * tmp_23;
      real_t tmp_25 = tmp_2 * tmp_23;
      real_t tmp_26 = p_affine_10_0 * tmp_21 + p_affine_10_1 * ( -tmp_24 - tmp_25 );
      real_t tmp_27 = tmp_14 * tmp_6 + tmp_16 * tmp_4;
      real_t tmp_28 = 2 * std::abs( std::pow( ( tmp_11 * tmp_11 ) + ( tmp_8 * tmp_8 ), 1.0 / 2.0 ) );
      real_t tmp_29 = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_28;
      real_t tmp_30 = tmp_7 * ( 0.23076534494715845 * tmp_8 + tmp_9 );
      real_t tmp_31 = tmp_7 * ( 0.23076534494715845 * tmp_11 + tmp_12 );
      real_t tmp_32 = tmp_2 * tmp_30 + tmp_31 * tmp_4 - 1.0 / 3.0;
      real_t tmp_33 = tmp_1 * tmp_30 + tmp_15 * tmp_31 - 1.0 / 3.0;
      real_t tmp_34 = tmp_1 * tmp_32 + tmp_33 * tmp_5;
      real_t tmp_35 = tmp_32 * tmp_6 + tmp_33 * tmp_4;
      real_t tmp_36 = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_28;
      real_t tmp_37 = tmp_7 * ( 0.5 * tmp_8 + tmp_9 );
      real_t tmp_38 = tmp_7 * ( 0.5 * tmp_11 + tmp_12 );
      real_t tmp_39 = tmp_2 * tmp_37 + tmp_38 * tmp_4 - 1.0 / 3.0;
      real_t tmp_40 = tmp_1 * tmp_37 + tmp_15 * tmp_38 - 1.0 / 3.0;
      real_t tmp_41 = tmp_1 * tmp_39 + tmp_40 * tmp_5;
      real_t tmp_42 = tmp_39 * tmp_6 + tmp_4 * tmp_40;
      real_t tmp_43 = 0.2844444444444445 * Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_28;
      real_t tmp_44 = tmp_7 * ( 0.7692346550528415 * tmp_8 + tmp_9 );
      real_t tmp_45 = tmp_7 * ( 0.7692346550528415 * tmp_11 + tmp_12 );
      real_t tmp_46 = tmp_2 * tmp_44 + tmp_4 * tmp_45 - 1.0 / 3.0;
      real_t tmp_47 = tmp_1 * tmp_44 + tmp_15 * tmp_45 - 1.0 / 3.0;
      real_t tmp_48 = tmp_1 * tmp_46 + tmp_47 * tmp_5;
      real_t tmp_49 = tmp_4 * tmp_47 + tmp_46 * tmp_6;
      real_t tmp_50 = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_28;
      real_t tmp_51 = tmp_7 * ( 0.95308992296933193 * tmp_8 + tmp_9 );
      real_t tmp_52 = tmp_7 * ( 0.95308992296933193 * tmp_11 + tmp_12 );
      real_t tmp_53 = tmp_2 * tmp_51 + tmp_4 * tmp_52 - 1.0 / 3.0;
      real_t tmp_54 = tmp_1 * tmp_51 + tmp_15 * tmp_52 - 1.0 / 3.0;
      real_t tmp_55 = tmp_1 * tmp_53 + tmp_5 * tmp_54;
      real_t tmp_56 = tmp_4 * tmp_54 + tmp_53 * tmp_6;
      real_t tmp_57 = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_28;
      real_t tmp_58 = p_affine_10_1 * tmp_19;
      real_t tmp_59 = p_affine_10_0 * tmp_19 + p_affine_10_1 * tmp_25;
      real_t tmp_60 = p_affine_10_1 * tmp_20;
      real_t tmp_61 = p_affine_10_0 * tmp_20 + p_affine_10_1 * tmp_24;
      real_t a_0_0  = tmp_29 * ( -tmp_17 * tmp_22 - tmp_26 * tmp_27 ) + tmp_36 * ( -tmp_22 * tmp_34 - tmp_26 * tmp_35 ) +
                     tmp_43 * ( -tmp_22 * tmp_41 - tmp_26 * tmp_42 ) + tmp_50 * ( -tmp_22 * tmp_48 - tmp_26 * tmp_49 ) +
                     tmp_57 * ( -tmp_22 * tmp_55 - tmp_26 * tmp_56 );
      real_t a_0_1 = tmp_29 * ( -tmp_17 * tmp_58 - tmp_27 * tmp_59 ) + tmp_36 * ( -tmp_34 * tmp_58 - tmp_35 * tmp_59 ) +
                     tmp_43 * ( -tmp_41 * tmp_58 - tmp_42 * tmp_59 ) + tmp_50 * ( -tmp_48 * tmp_58 - tmp_49 * tmp_59 ) +
                     tmp_57 * ( -tmp_55 * tmp_58 - tmp_56 * tmp_59 );
      real_t a_0_2 = tmp_29 * ( -tmp_17 * tmp_60 - tmp_27 * tmp_61 ) + tmp_36 * ( -tmp_34 * tmp_60 - tmp_35 * tmp_61 ) +
                     tmp_43 * ( -tmp_41 * tmp_60 - tmp_42 * tmp_61 ) + tmp_50 * ( -tmp_48 * tmp_60 - tmp_49 * tmp_61 ) +
                     tmp_57 * ( -tmp_55 * tmp_60 - tmp_56 * tmp_61 );
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

 private:
   void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t* out_0 ) const
   {
      *out_0 = callback2D( Point3D( { in_0, in_1 } ) );
   }
   std::function< real_t( const Point3D& ) > callback2D;
};

class EGEpsilonFormP1E_1 : public hyteg::dg::DGForm2D
{
 public:
   EGEpsilonFormP1E_1( std::function< real_t( const Point3D& ) > _callback2D )
   : callback2D( _callback2D )
   {}

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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id5  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id6  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id7  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id8  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id9  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id11 = 0;
      Scalar_Variable_Coefficient_2D_mu(
          0.063089014491502282 * p_affine_0_0 + 0.063089014491502227 * p_affine_1_0 + 0.87382197101699555 * p_affine_2_0,
          0.063089014491502282 * p_affine_0_1 + 0.063089014491502227 * p_affine_1_1 + 0.87382197101699555 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu(
          0.24928674517091043 * p_affine_0_0 + 0.24928674517091043 * p_affine_1_0 + 0.50142650965817914 * p_affine_2_0,
          0.24928674517091043 * p_affine_0_1 + 0.24928674517091043 * p_affine_1_1 + 0.50142650965817914 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu(
          0.63650249912139867 * p_affine_0_0 + 0.31035245103378439 * p_affine_1_0 + 0.053145049844816938 * p_affine_2_0,
          0.63650249912139867 * p_affine_0_1 + 0.31035245103378439 * p_affine_1_1 + 0.053145049844816938 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu(
          0.053145049844816938 * p_affine_0_0 + 0.63650249912139867 * p_affine_1_0 + 0.31035245103378439 * p_affine_2_0,
          0.053145049844816938 * p_affine_0_1 + 0.63650249912139867 * p_affine_1_1 + 0.31035245103378439 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu(
          0.063089014491502227 * p_affine_0_0 + 0.87382197101699555 * p_affine_1_0 + 0.063089014491502227 * p_affine_2_0,
          0.063089014491502227 * p_affine_0_1 + 0.87382197101699555 * p_affine_1_1 + 0.063089014491502227 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      Scalar_Variable_Coefficient_2D_mu(
          0.24928674517091043 * p_affine_0_0 + 0.50142650965817914 * p_affine_1_0 + 0.24928674517091043 * p_affine_2_0,
          0.24928674517091043 * p_affine_0_1 + 0.50142650965817914 * p_affine_1_1 + 0.24928674517091043 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id5 );
      Scalar_Variable_Coefficient_2D_mu(
          0.87382197101699566 * p_affine_0_0 + 0.063089014491502227 * p_affine_1_0 + 0.063089014491502227 * p_affine_2_0,
          0.87382197101699566 * p_affine_0_1 + 0.063089014491502227 * p_affine_1_1 + 0.063089014491502227 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id6 );
      Scalar_Variable_Coefficient_2D_mu(
          0.50142650965817914 * p_affine_0_0 + 0.24928674517091043 * p_affine_1_0 + 0.24928674517091043 * p_affine_2_0,
          0.50142650965817914 * p_affine_0_1 + 0.24928674517091043 * p_affine_1_1 + 0.24928674517091043 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id7 );
      Scalar_Variable_Coefficient_2D_mu(
          0.053145049844816938 * p_affine_0_0 + 0.31035245103378439 * p_affine_1_0 + 0.63650249912139867 * p_affine_2_0,
          0.053145049844816938 * p_affine_0_1 + 0.31035245103378439 * p_affine_1_1 + 0.63650249912139867 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id8 );
      Scalar_Variable_Coefficient_2D_mu(
          0.63650249912139867 * p_affine_0_0 + 0.053145049844816938 * p_affine_1_0 + 0.31035245103378439 * p_affine_2_0,
          0.63650249912139867 * p_affine_0_1 + 0.053145049844816938 * p_affine_1_1 + 0.31035245103378439 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id9 );
      Scalar_Variable_Coefficient_2D_mu(
          0.31035245103378439 * p_affine_0_0 + 0.63650249912139867 * p_affine_1_0 + 0.053145049844816938 * p_affine_2_0,
          0.31035245103378439 * p_affine_0_1 + 0.63650249912139867 * p_affine_1_1 + 0.053145049844816938 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id10 );
      Scalar_Variable_Coefficient_2D_mu(
          0.31035245103378439 * p_affine_0_0 + 0.053145049844816938 * p_affine_1_0 + 0.63650249912139867 * p_affine_2_0,
          0.31035245103378439 * p_affine_0_1 + 0.053145049844816938 * p_affine_1_1 + 0.63650249912139867 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id11 );
      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = -p_affine_0_1;
      real_t tmp_3  = p_affine_2_1 + tmp_2;
      real_t tmp_4  = tmp_1 * tmp_3;
      real_t tmp_5  = p_affine_2_0 + tmp_0;
      real_t tmp_6  = p_affine_1_1 + tmp_2;
      real_t tmp_7  = 1.0 / ( tmp_4 - tmp_5 * tmp_6 );
      real_t tmp_8  = 1.0 * tmp_7;
      real_t tmp_9  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = tmp_8 * tmp_9;
      real_t tmp_11 = tmp_10 * tmp_6 + tmp_4 * tmp_8;
      real_t tmp_12 = Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_11;
      real_t tmp_13 = -2 * tmp_1 * tmp_8 - 2 * tmp_10;
      real_t tmp_14 = 0.5 * tmp_7;
      real_t tmp_15 = tmp_1 * tmp_14;
      real_t tmp_16 = tmp_14 * tmp_3;
      real_t tmp_17 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_18 = tmp_14 * tmp_17;
      real_t tmp_19 = tmp_15 * tmp_5 + tmp_15 * tmp_9 + tmp_16 * tmp_6 + tmp_18 * tmp_3;
      real_t tmp_20 = Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_19;
      real_t tmp_21 = -4 * tmp_16 - 4 * tmp_18;
      real_t tmp_22 = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                                p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_23 = 0.025422453185103409 * tmp_22;
      real_t tmp_24 = tmp_11 * tmp_13;
      real_t tmp_25 = tmp_19 * tmp_21;
      real_t tmp_26 = 0.058393137863189684 * tmp_22;
      real_t tmp_27 = 0.041425537809186785 * tmp_22;
      real_t tmp_28 = 0.041425537809186785 * tmp_22;
      real_t tmp_29 = 0.025422453185103409 * tmp_22;
      real_t tmp_30 = 0.058393137863189684 * tmp_22;
      real_t tmp_31 = 0.025422453185103409 * tmp_22;
      real_t tmp_32 = 0.058393137863189684 * tmp_22;
      real_t tmp_33 = 0.041425537809186785 * tmp_22;
      real_t tmp_34 = 0.041425537809186785 * tmp_22;
      real_t tmp_35 = 0.041425537809186785 * tmp_22;
      real_t tmp_36 = 0.041425537809186785 * tmp_22;
      real_t tmp_37 = 2.0 * tmp_7;
      real_t tmp_38 = tmp_37 * tmp_9;
      real_t tmp_39 = tmp_3 * tmp_37;
      real_t tmp_40 = tmp_11 * tmp_38;
      real_t tmp_41 = tmp_19 * tmp_39;
      real_t tmp_42 = tmp_1 * tmp_37;
      real_t tmp_43 = tmp_17 * tmp_37;
      real_t tmp_44 = tmp_11 * tmp_42;
      real_t tmp_45 = tmp_19 * tmp_43;
      real_t a_0_0 =
          tmp_23 * ( tmp_12 * tmp_13 + tmp_20 * tmp_21 ) +
          tmp_26 * ( Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_25 ) +
          tmp_27 * ( Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_25 ) +
          tmp_28 * ( Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_25 ) +
          tmp_29 * ( Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_25 ) +
          tmp_30 * ( Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_25 ) +
          tmp_31 * ( Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_25 ) +
          tmp_32 * ( Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_25 ) +
          tmp_33 * ( Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_25 ) +
          tmp_34 * ( Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_25 ) +
          tmp_35 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_25 ) +
          tmp_36 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_24 + Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_25 );
      real_t a_1_0 =
          tmp_23 * ( tmp_12 * tmp_38 + tmp_20 * tmp_39 ) +
          tmp_26 * ( Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_41 ) +
          tmp_27 * ( Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_41 ) +
          tmp_28 * ( Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_41 ) +
          tmp_29 * ( Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_41 ) +
          tmp_30 * ( Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_41 ) +
          tmp_31 * ( Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_41 ) +
          tmp_32 * ( Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_41 ) +
          tmp_33 * ( Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_41 ) +
          tmp_34 * ( Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_41 ) +
          tmp_35 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_41 ) +
          tmp_36 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_40 + Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_41 );
      real_t a_2_0 =
          tmp_23 * ( tmp_12 * tmp_42 + tmp_20 * tmp_43 ) +
          tmp_26 * ( Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_45 ) +
          tmp_27 * ( Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_45 ) +
          tmp_28 * ( Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_45 ) +
          tmp_29 * ( Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_45 ) +
          tmp_30 * ( Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_45 ) +
          tmp_31 * ( Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_45 ) +
          tmp_32 * ( Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_45 ) +
          tmp_33 * ( Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_45 ) +
          tmp_34 * ( Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_45 ) +
          tmp_35 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_45 ) +
          tmp_36 *
              ( Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_44 + Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_45 );
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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.95308992296933193 * p_affine_6_0 + 0.046910077030668018 * p_affine_7_0,
                                         0.95308992296933193 * p_affine_6_1 + 0.046910077030668018 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.7692346550528415 * p_affine_6_0 + 0.23076534494715845 * p_affine_7_0,
                                         0.7692346550528415 * p_affine_6_1 + 0.23076534494715845 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.5 * p_affine_6_0 + 0.5 * p_affine_7_0,
                                         0.5 * p_affine_6_1 + 0.5 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.2307653449471585 * p_affine_6_0 + 0.7692346550528415 * p_affine_7_0,
                                         0.2307653449471585 * p_affine_6_1 + 0.7692346550528415 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.046910077030668074 * p_affine_6_0 + 0.95308992296933193 * p_affine_7_0,
                                         0.046910077030668074 * p_affine_6_1 + 0.95308992296933193 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3  = -p_affine_0_1;
      real_t tmp_4  = p_affine_2_1 + tmp_3;
      real_t tmp_5  = tmp_1 * tmp_4;
      real_t tmp_6  = p_affine_2_0 + tmp_0;
      real_t tmp_7  = p_affine_1_1 + tmp_3;
      real_t tmp_8  = 1.0 / ( tmp_5 - tmp_6 * tmp_7 );
      real_t tmp_9  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + tmp_3;
      real_t tmp_11 = tmp_8 * ( tmp_10 + 0.046910077030668018 * tmp_9 );
      real_t tmp_12 = tmp_11 * tmp_2;
      real_t tmp_13 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_14 = p_affine_6_0 + tmp_0;
      real_t tmp_15 = tmp_8 * ( 0.046910077030668018 * tmp_13 + tmp_14 );
      real_t tmp_16 = tmp_15 * tmp_4;
      real_t tmp_17 = tmp_12 + tmp_16;
      real_t tmp_18 = tmp_17 - 1.0 / 3.0;
      real_t tmp_19 = tmp_1 * tmp_11;
      real_t tmp_20 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_21 = tmp_15 * tmp_20;
      real_t tmp_22 = tmp_19 + tmp_21;
      real_t tmp_23 = tmp_22 - 1.0 / 3.0;
      real_t tmp_24 = p_affine_10_1 * ( tmp_1 * tmp_18 + tmp_23 * tmp_6 );
      real_t tmp_25 = 0.5 * tmp_8;
      real_t tmp_26 = tmp_25 * tmp_4;
      real_t tmp_27 = tmp_20 * tmp_25;
      real_t tmp_28 = -tmp_26 - tmp_27;
      real_t tmp_29 = 0.5 * tmp_28;
      real_t tmp_30 = tmp_18 * tmp_7 + tmp_23 * tmp_4;
      real_t tmp_31 = 1.0 * tmp_8;
      real_t tmp_32 = tmp_1 * tmp_31;
      real_t tmp_33 = tmp_2 * tmp_31;
      real_t tmp_34 = 0.5 * p_affine_10_0 * tmp_28 + 0.5 * p_affine_10_1 * ( -tmp_32 - tmp_33 );
      real_t tmp_35 = -tmp_12 - tmp_16 - tmp_19 - tmp_21 + 1;
      real_t tmp_36 = std::abs( std::pow( ( tmp_13 * tmp_13 ) + ( tmp_9 * tmp_9 ), 1.0 / 2.0 ) );
      real_t tmp_37 = 6 / tmp_36;
      real_t tmp_38 = tmp_30 * tmp_37;
      real_t tmp_39 = tmp_1 * tmp_25;
      real_t tmp_40 = 0.5 * p_affine_10_0 * ( tmp_2 * tmp_39 + tmp_26 * tmp_7 + tmp_27 * tmp_4 + tmp_39 * tmp_6 ) +
                      0.5 * p_affine_10_1 * ( tmp_31 * tmp_5 + tmp_33 * tmp_7 );
      real_t tmp_41  = 2 * tmp_36;
      real_t tmp_42  = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_41;
      real_t tmp_43  = tmp_8 * ( tmp_10 + 0.23076534494715845 * tmp_9 );
      real_t tmp_44  = tmp_2 * tmp_43;
      real_t tmp_45  = tmp_8 * ( 0.23076534494715845 * tmp_13 + tmp_14 );
      real_t tmp_46  = tmp_4 * tmp_45;
      real_t tmp_47  = tmp_44 + tmp_46;
      real_t tmp_48  = tmp_47 - 1.0 / 3.0;
      real_t tmp_49  = tmp_1 * tmp_43;
      real_t tmp_50  = tmp_20 * tmp_45;
      real_t tmp_51  = tmp_49 + tmp_50;
      real_t tmp_52  = tmp_51 - 1.0 / 3.0;
      real_t tmp_53  = tmp_1 * tmp_48 + tmp_52 * tmp_6;
      real_t tmp_54  = p_affine_10_1 * tmp_29;
      real_t tmp_55  = tmp_4 * tmp_52 + tmp_48 * tmp_7;
      real_t tmp_56  = -tmp_44 - tmp_46 - tmp_49 - tmp_50 + 1;
      real_t tmp_57  = tmp_37 * tmp_55;
      real_t tmp_58  = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_41;
      real_t tmp_59  = tmp_8 * ( tmp_10 + 0.5 * tmp_9 );
      real_t tmp_60  = tmp_2 * tmp_59;
      real_t tmp_61  = tmp_8 * ( 0.5 * tmp_13 + tmp_14 );
      real_t tmp_62  = tmp_4 * tmp_61;
      real_t tmp_63  = tmp_60 + tmp_62;
      real_t tmp_64  = tmp_63 - 1.0 / 3.0;
      real_t tmp_65  = tmp_1 * tmp_59;
      real_t tmp_66  = tmp_20 * tmp_61;
      real_t tmp_67  = tmp_65 + tmp_66;
      real_t tmp_68  = tmp_67 - 1.0 / 3.0;
      real_t tmp_69  = tmp_1 * tmp_64 + tmp_6 * tmp_68;
      real_t tmp_70  = tmp_4 * tmp_68 + tmp_64 * tmp_7;
      real_t tmp_71  = -tmp_60 - tmp_62 - tmp_65 - tmp_66 + 1;
      real_t tmp_72  = tmp_37 * tmp_70;
      real_t tmp_73  = 0.2844444444444445 * Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_41;
      real_t tmp_74  = tmp_8 * ( tmp_10 + 0.7692346550528415 * tmp_9 );
      real_t tmp_75  = tmp_2 * tmp_74;
      real_t tmp_76  = tmp_8 * ( 0.7692346550528415 * tmp_13 + tmp_14 );
      real_t tmp_77  = tmp_4 * tmp_76;
      real_t tmp_78  = tmp_75 + tmp_77;
      real_t tmp_79  = tmp_78 - 1.0 / 3.0;
      real_t tmp_80  = tmp_1 * tmp_74;
      real_t tmp_81  = tmp_20 * tmp_76;
      real_t tmp_82  = tmp_80 + tmp_81;
      real_t tmp_83  = tmp_82 - 1.0 / 3.0;
      real_t tmp_84  = tmp_1 * tmp_79 + tmp_6 * tmp_83;
      real_t tmp_85  = tmp_4 * tmp_83 + tmp_7 * tmp_79;
      real_t tmp_86  = -tmp_75 - tmp_77 - tmp_80 - tmp_81 + 1;
      real_t tmp_87  = tmp_37 * tmp_85;
      real_t tmp_88  = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_41;
      real_t tmp_89  = tmp_8 * ( tmp_10 + 0.95308992296933193 * tmp_9 );
      real_t tmp_90  = tmp_2 * tmp_89;
      real_t tmp_91  = tmp_8 * ( 0.95308992296933193 * tmp_13 + tmp_14 );
      real_t tmp_92  = tmp_4 * tmp_91;
      real_t tmp_93  = tmp_90 + tmp_92;
      real_t tmp_94  = tmp_93 - 1.0 / 3.0;
      real_t tmp_95  = tmp_1 * tmp_89;
      real_t tmp_96  = tmp_20 * tmp_91;
      real_t tmp_97  = tmp_95 + tmp_96;
      real_t tmp_98  = tmp_97 - 1.0 / 3.0;
      real_t tmp_99  = tmp_1 * tmp_94 + tmp_6 * tmp_98;
      real_t tmp_100 = tmp_4 * tmp_98 + tmp_7 * tmp_94;
      real_t tmp_101 = -tmp_90 - tmp_92 - tmp_95 - tmp_96 + 1;
      real_t tmp_102 = tmp_100 * tmp_37;
      real_t tmp_103 = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_41;
      real_t tmp_104 = 0.25 * tmp_8;
      real_t tmp_105 = tmp_104 * tmp_4;
      real_t tmp_106 = 0.5 * p_affine_10_0 * tmp_26 + 0.5 * p_affine_10_1 * tmp_33;
      real_t tmp_107 = p_affine_10_1 * tmp_105;
      real_t tmp_108 = tmp_104 * tmp_20;
      real_t tmp_109 = 0.5 * p_affine_10_0 * tmp_27 + 0.5 * p_affine_10_1 * tmp_32;
      real_t tmp_110 = p_affine_10_1 * tmp_108;
      real_t a_0_0   = tmp_103 * ( -tmp_100 * tmp_34 + tmp_101 * tmp_102 - tmp_101 * tmp_40 - tmp_54 * tmp_99 ) +
                     tmp_42 * ( -tmp_24 * tmp_29 - tmp_30 * tmp_34 + tmp_35 * tmp_38 - tmp_35 * tmp_40 ) +
                     tmp_58 * ( -tmp_34 * tmp_55 - tmp_40 * tmp_56 - tmp_53 * tmp_54 + tmp_56 * tmp_57 ) +
                     tmp_73 * ( -tmp_34 * tmp_70 - tmp_40 * tmp_71 - tmp_54 * tmp_69 + tmp_71 * tmp_72 ) +
                     tmp_88 * ( -tmp_34 * tmp_85 - tmp_40 * tmp_86 - tmp_54 * tmp_84 + tmp_86 * tmp_87 );
      real_t a_1_0 = tmp_103 * ( -tmp_100 * tmp_106 + tmp_102 * tmp_93 - tmp_107 * tmp_99 - tmp_40 * tmp_93 ) +
                     tmp_42 * ( -tmp_105 * tmp_24 - tmp_106 * tmp_30 + tmp_17 * tmp_38 - tmp_17 * tmp_40 ) +
                     tmp_58 * ( -tmp_106 * tmp_55 - tmp_107 * tmp_53 - tmp_40 * tmp_47 + tmp_47 * tmp_57 ) +
                     tmp_73 * ( -tmp_106 * tmp_70 - tmp_107 * tmp_69 - tmp_40 * tmp_63 + tmp_63 * tmp_72 ) +
                     tmp_88 * ( -tmp_106 * tmp_85 - tmp_107 * tmp_84 - tmp_40 * tmp_78 + tmp_78 * tmp_87 );
      real_t a_2_0 = tmp_103 * ( -tmp_100 * tmp_109 + tmp_102 * tmp_97 - tmp_110 * tmp_99 - tmp_40 * tmp_97 ) +
                     tmp_42 * ( -tmp_108 * tmp_24 - tmp_109 * tmp_30 + tmp_22 * tmp_38 - tmp_22 * tmp_40 ) +
                     tmp_58 * ( -tmp_109 * tmp_55 - tmp_110 * tmp_53 - tmp_40 * tmp_51 + tmp_51 * tmp_57 ) +
                     tmp_73 * ( -tmp_109 * tmp_70 - tmp_110 * tmp_69 - tmp_40 * tmp_67 + tmp_67 * tmp_72 ) +
                     tmp_88 * ( -tmp_109 * tmp_85 - tmp_110 * tmp_84 - tmp_40 * tmp_82 + tmp_82 * tmp_87 );
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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.95308992296933193 * p_affine_6_0 + 0.046910077030668018 * p_affine_7_0,
                                         0.95308992296933193 * p_affine_6_1 + 0.046910077030668018 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.7692346550528415 * p_affine_6_0 + 0.23076534494715845 * p_affine_7_0,
                                         0.7692346550528415 * p_affine_6_1 + 0.23076534494715845 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.5 * p_affine_6_0 + 0.5 * p_affine_7_0,
                                         0.5 * p_affine_6_1 + 0.5 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.2307653449471585 * p_affine_6_0 + 0.7692346550528415 * p_affine_7_0,
                                         0.2307653449471585 * p_affine_6_1 + 0.7692346550528415 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.046910077030668074 * p_affine_6_0 + 0.95308992296933193 * p_affine_7_0,
                                         0.046910077030668074 * p_affine_6_1 + 0.95308992296933193 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      real_t tmp_0  = -p_affine_3_0;
      real_t tmp_1  = p_affine_4_0 + tmp_0;
      real_t tmp_2  = p_affine_3_0 - p_affine_5_0;
      real_t tmp_3  = -p_affine_3_1;
      real_t tmp_4  = p_affine_5_1 + tmp_3;
      real_t tmp_5  = tmp_1 * tmp_4;
      real_t tmp_6  = p_affine_5_0 + tmp_0;
      real_t tmp_7  = p_affine_4_1 + tmp_3;
      real_t tmp_8  = 1.0 / ( tmp_5 - tmp_6 * tmp_7 );
      real_t tmp_9  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + 0.046910077030668018 * tmp_9;
      real_t tmp_11 = tmp_8 * ( tmp_10 + tmp_3 );
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + 0.046910077030668018 * tmp_12;
      real_t tmp_14 = tmp_8 * ( tmp_0 + tmp_13 );
      real_t tmp_15 = tmp_11 * tmp_2 + tmp_14 * tmp_4 - 1.0 / 3.0;
      real_t tmp_16 = p_affine_3_1 - p_affine_4_1;
      real_t tmp_17 = tmp_1 * tmp_11 + tmp_14 * tmp_16 - 1.0 / 3.0;
      real_t tmp_18 = p_affine_10_1 * ( tmp_1 * tmp_15 + tmp_17 * tmp_6 );
      real_t tmp_19 = -p_affine_0_1;
      real_t tmp_20 = p_affine_2_1 + tmp_19;
      real_t tmp_21 = -p_affine_0_0;
      real_t tmp_22 = p_affine_1_0 + tmp_21;
      real_t tmp_23 = 1.0 / ( tmp_20 * tmp_22 - ( p_affine_1_1 + tmp_19 ) * ( p_affine_2_0 + tmp_21 ) );
      real_t tmp_24 = 0.5 * tmp_23;
      real_t tmp_25 = tmp_20 * tmp_24;
      real_t tmp_26 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_27 = tmp_24 * tmp_26;
      real_t tmp_28 = -tmp_25 - tmp_27;
      real_t tmp_29 = 0.5 * tmp_28;
      real_t tmp_30 = tmp_15 * tmp_7 + tmp_17 * tmp_4;
      real_t tmp_31 = 1.0 * tmp_23;
      real_t tmp_32 = tmp_22 * tmp_31;
      real_t tmp_33 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_34 = tmp_31 * tmp_33;
      real_t tmp_35 = 0.5 * p_affine_10_0 * tmp_28 + 0.5 * p_affine_10_1 * ( -tmp_32 - tmp_34 );
      real_t tmp_36 = tmp_23 * ( tmp_10 + tmp_19 );
      real_t tmp_37 = tmp_22 * tmp_36;
      real_t tmp_38 = tmp_33 * tmp_36;
      real_t tmp_39 = tmp_23 * ( tmp_13 + tmp_21 );
      real_t tmp_40 = tmp_20 * tmp_39;
      real_t tmp_41 = tmp_26 * tmp_39;
      real_t tmp_42 = -tmp_37 - tmp_38 - tmp_40 - tmp_41 + 1;
      real_t tmp_43 = std::abs( std::pow( ( tmp_12 * tmp_12 ) + ( tmp_9 * tmp_9 ), 1.0 / 2.0 ) );
      real_t tmp_44 = 6 / tmp_43;
      real_t tmp_45 = tmp_30 * tmp_44;
      real_t tmp_46 = 1.0 * tmp_8;
      real_t tmp_47 = 0.5 * tmp_8;
      real_t tmp_48 = tmp_1 * tmp_47;
      real_t tmp_49 = tmp_4 * tmp_47;
      real_t tmp_50 = 0.5 * p_affine_10_0 * ( tmp_16 * tmp_49 + tmp_2 * tmp_48 + tmp_48 * tmp_6 + tmp_49 * tmp_7 ) +
                      0.5 * p_affine_10_1 * ( tmp_2 * tmp_46 * tmp_7 + tmp_46 * tmp_5 );
      real_t tmp_51  = 2 * tmp_43;
      real_t tmp_52  = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_51;
      real_t tmp_53  = p_affine_6_1 + 0.23076534494715845 * tmp_9;
      real_t tmp_54  = tmp_8 * ( tmp_3 + tmp_53 );
      real_t tmp_55  = p_affine_6_0 + 0.23076534494715845 * tmp_12;
      real_t tmp_56  = tmp_8 * ( tmp_0 + tmp_55 );
      real_t tmp_57  = tmp_2 * tmp_54 + tmp_4 * tmp_56 - 1.0 / 3.0;
      real_t tmp_58  = tmp_1 * tmp_54 + tmp_16 * tmp_56 - 1.0 / 3.0;
      real_t tmp_59  = tmp_1 * tmp_57 + tmp_58 * tmp_6;
      real_t tmp_60  = p_affine_10_1 * tmp_29;
      real_t tmp_61  = tmp_4 * tmp_58 + tmp_57 * tmp_7;
      real_t tmp_62  = tmp_23 * ( tmp_19 + tmp_53 );
      real_t tmp_63  = tmp_22 * tmp_62;
      real_t tmp_64  = tmp_33 * tmp_62;
      real_t tmp_65  = tmp_23 * ( tmp_21 + tmp_55 );
      real_t tmp_66  = tmp_20 * tmp_65;
      real_t tmp_67  = tmp_26 * tmp_65;
      real_t tmp_68  = -tmp_63 - tmp_64 - tmp_66 - tmp_67 + 1;
      real_t tmp_69  = tmp_44 * tmp_61;
      real_t tmp_70  = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_51;
      real_t tmp_71  = p_affine_6_1 + 0.5 * tmp_9;
      real_t tmp_72  = tmp_8 * ( tmp_3 + tmp_71 );
      real_t tmp_73  = p_affine_6_0 + 0.5 * tmp_12;
      real_t tmp_74  = tmp_8 * ( tmp_0 + tmp_73 );
      real_t tmp_75  = tmp_2 * tmp_72 + tmp_4 * tmp_74 - 1.0 / 3.0;
      real_t tmp_76  = tmp_1 * tmp_72 + tmp_16 * tmp_74 - 1.0 / 3.0;
      real_t tmp_77  = tmp_1 * tmp_75 + tmp_6 * tmp_76;
      real_t tmp_78  = tmp_4 * tmp_76 + tmp_7 * tmp_75;
      real_t tmp_79  = tmp_23 * ( tmp_19 + tmp_71 );
      real_t tmp_80  = tmp_22 * tmp_79;
      real_t tmp_81  = tmp_33 * tmp_79;
      real_t tmp_82  = tmp_23 * ( tmp_21 + tmp_73 );
      real_t tmp_83  = tmp_20 * tmp_82;
      real_t tmp_84  = tmp_26 * tmp_82;
      real_t tmp_85  = -tmp_80 - tmp_81 - tmp_83 - tmp_84 + 1;
      real_t tmp_86  = tmp_44 * tmp_78;
      real_t tmp_87  = 0.2844444444444445 * Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_51;
      real_t tmp_88  = p_affine_6_1 + 0.7692346550528415 * tmp_9;
      real_t tmp_89  = tmp_8 * ( tmp_3 + tmp_88 );
      real_t tmp_90  = p_affine_6_0 + 0.7692346550528415 * tmp_12;
      real_t tmp_91  = tmp_8 * ( tmp_0 + tmp_90 );
      real_t tmp_92  = tmp_2 * tmp_89 + tmp_4 * tmp_91 - 1.0 / 3.0;
      real_t tmp_93  = tmp_1 * tmp_89 + tmp_16 * tmp_91 - 1.0 / 3.0;
      real_t tmp_94  = tmp_1 * tmp_92 + tmp_6 * tmp_93;
      real_t tmp_95  = tmp_4 * tmp_93 + tmp_7 * tmp_92;
      real_t tmp_96  = tmp_23 * ( tmp_19 + tmp_88 );
      real_t tmp_97  = tmp_22 * tmp_96;
      real_t tmp_98  = tmp_33 * tmp_96;
      real_t tmp_99  = tmp_23 * ( tmp_21 + tmp_90 );
      real_t tmp_100 = tmp_20 * tmp_99;
      real_t tmp_101 = tmp_26 * tmp_99;
      real_t tmp_102 = -tmp_100 - tmp_101 - tmp_97 - tmp_98 + 1;
      real_t tmp_103 = tmp_44 * tmp_95;
      real_t tmp_104 = 0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_51;
      real_t tmp_105 = p_affine_6_1 + 0.95308992296933193 * tmp_9;
      real_t tmp_106 = tmp_8 * ( tmp_105 + tmp_3 );
      real_t tmp_107 = p_affine_6_0 + 0.95308992296933193 * tmp_12;
      real_t tmp_108 = tmp_8 * ( tmp_0 + tmp_107 );
      real_t tmp_109 = tmp_106 * tmp_2 + tmp_108 * tmp_4 - 1.0 / 3.0;
      real_t tmp_110 = tmp_1 * tmp_106 + tmp_108 * tmp_16 - 1.0 / 3.0;
      real_t tmp_111 = tmp_1 * tmp_109 + tmp_110 * tmp_6;
      real_t tmp_112 = tmp_109 * tmp_7 + tmp_110 * tmp_4;
      real_t tmp_113 = tmp_23 * ( tmp_105 + tmp_19 );
      real_t tmp_114 = tmp_113 * tmp_22;
      real_t tmp_115 = tmp_113 * tmp_33;
      real_t tmp_116 = tmp_23 * ( tmp_107 + tmp_21 );
      real_t tmp_117 = tmp_116 * tmp_20;
      real_t tmp_118 = tmp_116 * tmp_26;
      real_t tmp_119 = -tmp_114 - tmp_115 - tmp_117 - tmp_118 + 1;
      real_t tmp_120 = tmp_112 * tmp_44;
      real_t tmp_121 = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_51;
      real_t tmp_122 = 0.25 * tmp_23;
      real_t tmp_123 = tmp_122 * tmp_20;
      real_t tmp_124 = 0.5 * p_affine_10_0 * tmp_25 + 0.5 * p_affine_10_1 * tmp_34;
      real_t tmp_125 = tmp_38 + tmp_40;
      real_t tmp_126 = p_affine_10_1 * tmp_123;
      real_t tmp_127 = tmp_64 + tmp_66;
      real_t tmp_128 = tmp_81 + tmp_83;
      real_t tmp_129 = tmp_100 + tmp_98;
      real_t tmp_130 = tmp_115 + tmp_117;
      real_t tmp_131 = tmp_122 * tmp_26;
      real_t tmp_132 = 0.5 * p_affine_10_0 * tmp_27 + 0.5 * p_affine_10_1 * tmp_32;
      real_t tmp_133 = tmp_37 + tmp_41;
      real_t tmp_134 = p_affine_10_1 * tmp_131;
      real_t tmp_135 = tmp_63 + tmp_67;
      real_t tmp_136 = tmp_80 + tmp_84;
      real_t tmp_137 = tmp_101 + tmp_97;
      real_t tmp_138 = tmp_114 + tmp_118;
      real_t a_0_0   = tmp_104 * ( -tmp_102 * tmp_103 - tmp_102 * tmp_50 + tmp_35 * tmp_95 + tmp_60 * tmp_94 ) +
                     tmp_121 * ( tmp_111 * tmp_60 + tmp_112 * tmp_35 - tmp_119 * tmp_120 - tmp_119 * tmp_50 ) +
                     tmp_52 * ( tmp_18 * tmp_29 + tmp_30 * tmp_35 - tmp_42 * tmp_45 - tmp_42 * tmp_50 ) +
                     tmp_70 * ( tmp_35 * tmp_61 - tmp_50 * tmp_68 + tmp_59 * tmp_60 - tmp_68 * tmp_69 ) +
                     tmp_87 * ( tmp_35 * tmp_78 - tmp_50 * tmp_85 + tmp_60 * tmp_77 - tmp_85 * tmp_86 );
      real_t a_1_0 = tmp_104 * ( -tmp_103 * tmp_129 + tmp_124 * tmp_95 + tmp_126 * tmp_94 - tmp_129 * tmp_50 ) +
                     tmp_121 * ( tmp_111 * tmp_126 + tmp_112 * tmp_124 - tmp_120 * tmp_130 - tmp_130 * tmp_50 ) +
                     tmp_52 * ( tmp_123 * tmp_18 + tmp_124 * tmp_30 - tmp_125 * tmp_45 - tmp_125 * tmp_50 ) +
                     tmp_70 * ( tmp_124 * tmp_61 + tmp_126 * tmp_59 - tmp_127 * tmp_50 - tmp_127 * tmp_69 ) +
                     tmp_87 * ( tmp_124 * tmp_78 + tmp_126 * tmp_77 - tmp_128 * tmp_50 - tmp_128 * tmp_86 );
      real_t a_2_0 = tmp_104 * ( -tmp_103 * tmp_137 + tmp_132 * tmp_95 + tmp_134 * tmp_94 - tmp_137 * tmp_50 ) +
                     tmp_121 * ( tmp_111 * tmp_134 + tmp_112 * tmp_132 - tmp_120 * tmp_138 - tmp_138 * tmp_50 ) +
                     tmp_52 * ( tmp_131 * tmp_18 + tmp_132 * tmp_30 - tmp_133 * tmp_45 - tmp_133 * tmp_50 ) +
                     tmp_70 * ( tmp_132 * tmp_61 + tmp_134 * tmp_59 - tmp_135 * tmp_50 - tmp_135 * tmp_69 ) +
                     tmp_87 * ( tmp_132 * tmp_78 + tmp_134 * tmp_77 - tmp_136 * tmp_50 - tmp_136 * tmp_86 );
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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.95308992296933193 * p_affine_6_0 + 0.046910077030668018 * p_affine_7_0,
                                         0.95308992296933193 * p_affine_6_1 + 0.046910077030668018 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.7692346550528415 * p_affine_6_0 + 0.23076534494715845 * p_affine_7_0,
                                         0.7692346550528415 * p_affine_6_1 + 0.23076534494715845 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.5 * p_affine_6_0 + 0.5 * p_affine_7_0,
                                         0.5 * p_affine_6_1 + 0.5 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.2307653449471585 * p_affine_6_0 + 0.7692346550528415 * p_affine_7_0,
                                         0.2307653449471585 * p_affine_6_1 + 0.7692346550528415 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.046910077030668074 * p_affine_6_0 + 0.95308992296933193 * p_affine_7_0,
                                         0.046910077030668074 * p_affine_6_1 + 0.95308992296933193 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
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
      real_t tmp_20 = 1.4215613103371365 * Scalar_Variable_Coefficient_2D_mu_out0_id0 *
                      ( tmp_3 * ( tmp_19 - 1.0 / 3.0 ) + tmp_4 * ( tmp_18 - 1.0 / 3.0 ) );
      real_t tmp_21 = tmp_5 * ( 0.23076534494715845 * tmp_6 + tmp_7 );
      real_t tmp_22 = tmp_1 * tmp_21;
      real_t tmp_23 = tmp_10 * tmp_21;
      real_t tmp_24 = tmp_5 * ( 0.23076534494715845 * tmp_12 + tmp_13 );
      real_t tmp_25 = tmp_24 * tmp_3;
      real_t tmp_26 = tmp_16 * tmp_24;
      real_t tmp_27 = tmp_23 + tmp_25;
      real_t tmp_28 = tmp_22 + tmp_26;
      real_t tmp_29 = 2.8717720229961969 * Scalar_Variable_Coefficient_2D_mu_out0_id1 *
                      ( tmp_3 * ( tmp_28 - 1.0 / 3.0 ) + tmp_4 * ( tmp_27 - 1.0 / 3.0 ) );
      real_t tmp_30 = tmp_5 * ( 0.5 * tmp_6 + tmp_7 );
      real_t tmp_31 = tmp_1 * tmp_30;
      real_t tmp_32 = tmp_10 * tmp_30;
      real_t tmp_33 = tmp_5 * ( 0.5 * tmp_12 + tmp_13 );
      real_t tmp_34 = tmp_3 * tmp_33;
      real_t tmp_35 = tmp_16 * tmp_33;
      real_t tmp_36 = tmp_32 + tmp_34;
      real_t tmp_37 = tmp_31 + tmp_35;
      real_t tmp_38 = 3.413333333333334 * Scalar_Variable_Coefficient_2D_mu_out0_id2 *
                      ( tmp_3 * ( tmp_37 - 1.0 / 3.0 ) + tmp_4 * ( tmp_36 - 1.0 / 3.0 ) );
      real_t tmp_39 = tmp_5 * ( 0.7692346550528415 * tmp_6 + tmp_7 );
      real_t tmp_40 = tmp_1 * tmp_39;
      real_t tmp_41 = tmp_10 * tmp_39;
      real_t tmp_42 = tmp_5 * ( 0.7692346550528415 * tmp_12 + tmp_13 );
      real_t tmp_43 = tmp_3 * tmp_42;
      real_t tmp_44 = tmp_16 * tmp_42;
      real_t tmp_45 = tmp_41 + tmp_43;
      real_t tmp_46 = tmp_40 + tmp_44;
      real_t tmp_47 = 2.8717720229961969 * Scalar_Variable_Coefficient_2D_mu_out0_id3 *
                      ( tmp_3 * ( tmp_46 - 1.0 / 3.0 ) + tmp_4 * ( tmp_45 - 1.0 / 3.0 ) );
      real_t tmp_48 = tmp_5 * ( 0.95308992296933193 * tmp_6 + tmp_7 );
      real_t tmp_49 = tmp_1 * tmp_48;
      real_t tmp_50 = tmp_10 * tmp_48;
      real_t tmp_51 = tmp_5 * ( 0.95308992296933193 * tmp_12 + tmp_13 );
      real_t tmp_52 = tmp_3 * tmp_51;
      real_t tmp_53 = tmp_16 * tmp_51;
      real_t tmp_54 = tmp_50 + tmp_52;
      real_t tmp_55 = tmp_49 + tmp_53;
      real_t tmp_56 = 1.4215613103371365 * Scalar_Variable_Coefficient_2D_mu_out0_id4 *
                      ( tmp_3 * ( tmp_55 - 1.0 / 3.0 ) + tmp_4 * ( tmp_54 - 1.0 / 3.0 ) );
      real_t a_0_0 = tmp_20 * ( -tmp_11 - tmp_15 - tmp_17 - tmp_9 + 1 ) + tmp_29 * ( -tmp_22 - tmp_23 - tmp_25 - tmp_26 + 1 ) +
                     tmp_38 * ( -tmp_31 - tmp_32 - tmp_34 - tmp_35 + 1 ) + tmp_47 * ( -tmp_40 - tmp_41 - tmp_43 - tmp_44 + 1 ) +
                     tmp_56 * ( -tmp_49 - tmp_50 - tmp_52 - tmp_53 + 1 );
      real_t a_1_0  = tmp_18 * tmp_20 + tmp_27 * tmp_29 + tmp_36 * tmp_38 + tmp_45 * tmp_47 + tmp_54 * tmp_56;
      real_t a_2_0  = tmp_19 * tmp_20 + tmp_28 * tmp_29 + tmp_37 * tmp_38 + tmp_46 * tmp_47 + tmp_55 * tmp_56;
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

 private:
   void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t* out_0 ) const
   {
      *out_0 = callback2D( Point3D( { in_0, in_1 } ) );
   }
   std::function< real_t( const Point3D& ) > callback2D;
};

class EGEpsilonFormEE : public hyteg::dg::DGForm2D
{
 public:
   EGEpsilonFormEE( std::function< real_t( const Point3D& ) > _callback2D )
   : callback2D( _callback2D )
   {}

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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id5  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id6  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id7  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id8  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id9  = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id11 = 0;
      Scalar_Variable_Coefficient_2D_mu(
          0.063089014491502282 * p_affine_0_0 + 0.063089014491502227 * p_affine_1_0 + 0.87382197101699555 * p_affine_2_0,
          0.063089014491502282 * p_affine_0_1 + 0.063089014491502227 * p_affine_1_1 + 0.87382197101699555 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu(
          0.24928674517091043 * p_affine_0_0 + 0.24928674517091043 * p_affine_1_0 + 0.50142650965817914 * p_affine_2_0,
          0.24928674517091043 * p_affine_0_1 + 0.24928674517091043 * p_affine_1_1 + 0.50142650965817914 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu(
          0.63650249912139867 * p_affine_0_0 + 0.31035245103378439 * p_affine_1_0 + 0.053145049844816938 * p_affine_2_0,
          0.63650249912139867 * p_affine_0_1 + 0.31035245103378439 * p_affine_1_1 + 0.053145049844816938 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu(
          0.053145049844816938 * p_affine_0_0 + 0.63650249912139867 * p_affine_1_0 + 0.31035245103378439 * p_affine_2_0,
          0.053145049844816938 * p_affine_0_1 + 0.63650249912139867 * p_affine_1_1 + 0.31035245103378439 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu(
          0.063089014491502227 * p_affine_0_0 + 0.87382197101699555 * p_affine_1_0 + 0.063089014491502227 * p_affine_2_0,
          0.063089014491502227 * p_affine_0_1 + 0.87382197101699555 * p_affine_1_1 + 0.063089014491502227 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      Scalar_Variable_Coefficient_2D_mu(
          0.24928674517091043 * p_affine_0_0 + 0.50142650965817914 * p_affine_1_0 + 0.24928674517091043 * p_affine_2_0,
          0.24928674517091043 * p_affine_0_1 + 0.50142650965817914 * p_affine_1_1 + 0.24928674517091043 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id5 );
      Scalar_Variable_Coefficient_2D_mu(
          0.87382197101699566 * p_affine_0_0 + 0.063089014491502227 * p_affine_1_0 + 0.063089014491502227 * p_affine_2_0,
          0.87382197101699566 * p_affine_0_1 + 0.063089014491502227 * p_affine_1_1 + 0.063089014491502227 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id6 );
      Scalar_Variable_Coefficient_2D_mu(
          0.50142650965817914 * p_affine_0_0 + 0.24928674517091043 * p_affine_1_0 + 0.24928674517091043 * p_affine_2_0,
          0.50142650965817914 * p_affine_0_1 + 0.24928674517091043 * p_affine_1_1 + 0.24928674517091043 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id7 );
      Scalar_Variable_Coefficient_2D_mu(
          0.053145049844816938 * p_affine_0_0 + 0.31035245103378439 * p_affine_1_0 + 0.63650249912139867 * p_affine_2_0,
          0.053145049844816938 * p_affine_0_1 + 0.31035245103378439 * p_affine_1_1 + 0.63650249912139867 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id8 );
      Scalar_Variable_Coefficient_2D_mu(
          0.63650249912139867 * p_affine_0_0 + 0.053145049844816938 * p_affine_1_0 + 0.31035245103378439 * p_affine_2_0,
          0.63650249912139867 * p_affine_0_1 + 0.053145049844816938 * p_affine_1_1 + 0.31035245103378439 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id9 );
      Scalar_Variable_Coefficient_2D_mu(
          0.31035245103378439 * p_affine_0_0 + 0.63650249912139867 * p_affine_1_0 + 0.053145049844816938 * p_affine_2_0,
          0.31035245103378439 * p_affine_0_1 + 0.63650249912139867 * p_affine_1_1 + 0.053145049844816938 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id10 );
      Scalar_Variable_Coefficient_2D_mu(
          0.31035245103378439 * p_affine_0_0 + 0.053145049844816938 * p_affine_1_0 + 0.63650249912139867 * p_affine_2_0,
          0.31035245103378439 * p_affine_0_1 + 0.053145049844816938 * p_affine_1_1 + 0.63650249912139867 * p_affine_2_1,
          &Scalar_Variable_Coefficient_2D_mu_out0_id11 );
      real_t tmp_0  = std::abs( p_affine_0_0 * p_affine_1_1 - p_affine_0_0 * p_affine_2_1 - p_affine_0_1 * p_affine_1_0 +
                               p_affine_0_1 * p_affine_2_0 + p_affine_1_0 * p_affine_2_1 - p_affine_1_1 * p_affine_2_0 );
      real_t tmp_1  = -p_affine_0_0;
      real_t tmp_2  = p_affine_1_0 + tmp_1;
      real_t tmp_3  = -p_affine_0_1;
      real_t tmp_4  = p_affine_2_1 + tmp_3;
      real_t tmp_5  = tmp_2 * tmp_4;
      real_t tmp_6  = p_affine_2_0 + tmp_1;
      real_t tmp_7  = p_affine_1_1 + tmp_3;
      real_t tmp_8  = 1.0 / ( tmp_5 - tmp_6 * tmp_7 );
      real_t tmp_9  = 1.0 * tmp_8;
      real_t tmp_10 = tmp_5 * tmp_9;
      real_t tmp_11 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_12 = ( ( tmp_10 + tmp_11 * tmp_6 * tmp_9 ) * ( tmp_10 + tmp_11 * tmp_6 * tmp_9 ) );
      real_t tmp_13 = 2 * Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_14 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_15 = ( ( tmp_10 + tmp_14 * tmp_7 * tmp_9 ) * ( tmp_10 + tmp_14 * tmp_7 * tmp_9 ) );
      real_t tmp_16 = tmp_2 * tmp_8;
      real_t tmp_17 = tmp_4 * tmp_8;
      real_t tmp_18 = 1.0 * ( ( tmp_11 * tmp_17 + tmp_14 * tmp_16 + tmp_16 * tmp_6 + tmp_17 * tmp_7 ) *
                              ( tmp_11 * tmp_17 + tmp_14 * tmp_16 + tmp_16 * tmp_6 + tmp_17 * tmp_7 ) );
      real_t tmp_19 = 2 * Scalar_Variable_Coefficient_2D_mu_out0_id1;
      real_t tmp_20 = 2 * Scalar_Variable_Coefficient_2D_mu_out0_id2;
      real_t tmp_21 = 2 * Scalar_Variable_Coefficient_2D_mu_out0_id3;
      real_t tmp_22 = 2 * Scalar_Variable_Coefficient_2D_mu_out0_id4;
      real_t tmp_23 = 2 * Scalar_Variable_Coefficient_2D_mu_out0_id5;
      real_t tmp_24 = 2 * Scalar_Variable_Coefficient_2D_mu_out0_id6;
      real_t tmp_25 = 2 * Scalar_Variable_Coefficient_2D_mu_out0_id7;
      real_t tmp_26 = 2 * Scalar_Variable_Coefficient_2D_mu_out0_id8;
      real_t tmp_27 = 2 * Scalar_Variable_Coefficient_2D_mu_out0_id9;
      real_t tmp_28 = 2 * Scalar_Variable_Coefficient_2D_mu_out0_id10;
      real_t tmp_29 = 2 * Scalar_Variable_Coefficient_2D_mu_out0_id11;
      real_t a_0_0  = 0.025422453185103409 * tmp_0 *
                         ( Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_18 + tmp_12 * tmp_13 + tmp_13 * tmp_15 ) +
                     0.058393137863189684 * tmp_0 *
                         ( Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_18 + tmp_12 * tmp_19 + tmp_15 * tmp_19 ) +
                     0.041425537809186785 * tmp_0 *
                         ( Scalar_Variable_Coefficient_2D_mu_out0_id10 * tmp_18 + tmp_12 * tmp_28 + tmp_15 * tmp_28 ) +
                     0.041425537809186785 * tmp_0 *
                         ( Scalar_Variable_Coefficient_2D_mu_out0_id11 * tmp_18 + tmp_12 * tmp_29 + tmp_15 * tmp_29 ) +
                     0.041425537809186785 * tmp_0 *
                         ( Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_18 + tmp_12 * tmp_20 + tmp_15 * tmp_20 ) +
                     0.041425537809186785 * tmp_0 *
                         ( Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_18 + tmp_12 * tmp_21 + tmp_15 * tmp_21 ) +
                     0.025422453185103409 * tmp_0 *
                         ( Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_18 + tmp_12 * tmp_22 + tmp_15 * tmp_22 ) +
                     0.058393137863189684 * tmp_0 *
                         ( Scalar_Variable_Coefficient_2D_mu_out0_id5 * tmp_18 + tmp_12 * tmp_23 + tmp_15 * tmp_23 ) +
                     0.025422453185103409 * tmp_0 *
                         ( Scalar_Variable_Coefficient_2D_mu_out0_id6 * tmp_18 + tmp_12 * tmp_24 + tmp_15 * tmp_24 ) +
                     0.058393137863189684 * tmp_0 *
                         ( Scalar_Variable_Coefficient_2D_mu_out0_id7 * tmp_18 + tmp_12 * tmp_25 + tmp_15 * tmp_25 ) +
                     0.041425537809186785 * tmp_0 *
                         ( Scalar_Variable_Coefficient_2D_mu_out0_id8 * tmp_18 + tmp_12 * tmp_26 + tmp_15 * tmp_26 ) +
                     0.041425537809186785 * tmp_0 *
                         ( Scalar_Variable_Coefficient_2D_mu_out0_id9 * tmp_18 + tmp_12 * tmp_27 + tmp_15 * tmp_27 );
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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.95308992296933193 * p_affine_6_0 + 0.046910077030668018 * p_affine_7_0,
                                         0.95308992296933193 * p_affine_6_1 + 0.046910077030668018 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.7692346550528415 * p_affine_6_0 + 0.23076534494715845 * p_affine_7_0,
                                         0.7692346550528415 * p_affine_6_1 + 0.23076534494715845 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.5 * p_affine_6_0 + 0.5 * p_affine_7_0,
                                         0.5 * p_affine_6_1 + 0.5 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.2307653449471585 * p_affine_6_0 + 0.7692346550528415 * p_affine_7_0,
                                         0.2307653449471585 * p_affine_6_1 + 0.7692346550528415 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.046910077030668074 * p_affine_6_0 + 0.95308992296933193 * p_affine_7_0,
                                         0.046910077030668074 * p_affine_6_1 + 0.95308992296933193 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      real_t tmp_0  = -p_affine_0_0;
      real_t tmp_1  = p_affine_1_0 + tmp_0;
      real_t tmp_2  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_3  = -p_affine_0_1;
      real_t tmp_4  = p_affine_2_1 + tmp_3;
      real_t tmp_5  = tmp_1 * tmp_4;
      real_t tmp_6  = p_affine_2_0 + tmp_0;
      real_t tmp_7  = p_affine_1_1 + tmp_3;
      real_t tmp_8  = 1.0 / ( tmp_5 - tmp_6 * tmp_7 );
      real_t tmp_9  = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10 = p_affine_6_1 + tmp_3;
      real_t tmp_11 = tmp_8 * ( tmp_10 + 0.046910077030668018 * tmp_9 );
      real_t tmp_12 = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13 = p_affine_6_0 + tmp_0;
      real_t tmp_14 = tmp_8 * ( 0.046910077030668018 * tmp_12 + tmp_13 );
      real_t tmp_15 = tmp_11 * tmp_2 + tmp_14 * tmp_4 - 1.0 / 3.0;
      real_t tmp_16 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_17 = tmp_1 * tmp_11 + tmp_14 * tmp_16 - 1.0 / 3.0;
      real_t tmp_18 = tmp_1 * tmp_15 + tmp_17 * tmp_6;
      real_t tmp_19 = tmp_15 * tmp_7 + tmp_17 * tmp_4;
      real_t tmp_20 = std::abs( std::pow( ( tmp_12 * tmp_12 ) + ( tmp_9 * tmp_9 ), 1.0 / 2.0 ) );
      real_t tmp_21 = 6 / tmp_20;
      real_t tmp_22 = 1.0 * tmp_8;
      real_t tmp_23 = tmp_22 * tmp_5;
      real_t tmp_24 = 0.5 * tmp_8;
      real_t tmp_25 = tmp_1 * tmp_24;
      real_t tmp_26 = tmp_24 * tmp_4;
      real_t tmp_27 = tmp_16 * tmp_26 + tmp_2 * tmp_25 + tmp_25 * tmp_6 + tmp_26 * tmp_7;
      real_t tmp_28 = 1.0 * p_affine_10_0 * ( tmp_16 * tmp_22 * tmp_6 + tmp_23 ) + 1.0 * p_affine_10_1 * tmp_27;
      real_t tmp_29 = 1.0 * p_affine_10_0 * tmp_27 + 1.0 * p_affine_10_1 * ( tmp_2 * tmp_22 * tmp_7 + tmp_23 );
      real_t tmp_30 = 2 * tmp_20;
      real_t tmp_31 = tmp_8 * ( tmp_10 + 0.23076534494715845 * tmp_9 );
      real_t tmp_32 = tmp_8 * ( 0.23076534494715845 * tmp_12 + tmp_13 );
      real_t tmp_33 = tmp_2 * tmp_31 + tmp_32 * tmp_4 - 1.0 / 3.0;
      real_t tmp_34 = tmp_1 * tmp_31 + tmp_16 * tmp_32 - 1.0 / 3.0;
      real_t tmp_35 = tmp_1 * tmp_33 + tmp_34 * tmp_6;
      real_t tmp_36 = tmp_33 * tmp_7 + tmp_34 * tmp_4;
      real_t tmp_37 = tmp_8 * ( tmp_10 + 0.5 * tmp_9 );
      real_t tmp_38 = tmp_8 * ( 0.5 * tmp_12 + tmp_13 );
      real_t tmp_39 = tmp_2 * tmp_37 + tmp_38 * tmp_4 - 1.0 / 3.0;
      real_t tmp_40 = tmp_1 * tmp_37 + tmp_16 * tmp_38 - 1.0 / 3.0;
      real_t tmp_41 = tmp_1 * tmp_39 + tmp_40 * tmp_6;
      real_t tmp_42 = tmp_39 * tmp_7 + tmp_4 * tmp_40;
      real_t tmp_43 = tmp_8 * ( tmp_10 + 0.7692346550528415 * tmp_9 );
      real_t tmp_44 = tmp_8 * ( 0.7692346550528415 * tmp_12 + tmp_13 );
      real_t tmp_45 = tmp_2 * tmp_43 + tmp_4 * tmp_44 - 1.0 / 3.0;
      real_t tmp_46 = tmp_1 * tmp_43 + tmp_16 * tmp_44 - 1.0 / 3.0;
      real_t tmp_47 = tmp_1 * tmp_45 + tmp_46 * tmp_6;
      real_t tmp_48 = tmp_4 * tmp_46 + tmp_45 * tmp_7;
      real_t tmp_49 = tmp_8 * ( tmp_10 + 0.95308992296933193 * tmp_9 );
      real_t tmp_50 = tmp_8 * ( 0.95308992296933193 * tmp_12 + tmp_13 );
      real_t tmp_51 = tmp_2 * tmp_49 + tmp_4 * tmp_50 - 1.0 / 3.0;
      real_t tmp_52 = tmp_1 * tmp_49 + tmp_16 * tmp_50 - 1.0 / 3.0;
      real_t tmp_53 = tmp_1 * tmp_51 + tmp_52 * tmp_6;
      real_t tmp_54 = tmp_4 * tmp_52 + tmp_51 * tmp_7;
      real_t a_0_0  = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_30 *
                         ( -tmp_18 * tmp_28 - tmp_19 * tmp_29 + tmp_21 * ( ( tmp_18 * tmp_18 ) + ( tmp_19 * tmp_19 ) ) ) +
                     0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_30 *
                         ( tmp_21 * ( ( tmp_35 * tmp_35 ) + ( tmp_36 * tmp_36 ) ) - tmp_28 * tmp_35 - tmp_29 * tmp_36 ) +
                     0.2844444444444445 * Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_30 *
                         ( tmp_21 * ( ( tmp_41 * tmp_41 ) + ( tmp_42 * tmp_42 ) ) - tmp_28 * tmp_41 - tmp_29 * tmp_42 ) +
                     0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_30 *
                         ( tmp_21 * ( ( tmp_47 * tmp_47 ) + ( tmp_48 * tmp_48 ) ) - tmp_28 * tmp_47 - tmp_29 * tmp_48 ) +
                     0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_30 *
                         ( tmp_21 * ( ( tmp_53 * tmp_53 ) + ( tmp_54 * tmp_54 ) ) - tmp_28 * tmp_53 - tmp_29 * tmp_54 );
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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.95308992296933193 * p_affine_6_0 + 0.046910077030668018 * p_affine_7_0,
                                         0.95308992296933193 * p_affine_6_1 + 0.046910077030668018 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.7692346550528415 * p_affine_6_0 + 0.23076534494715845 * p_affine_7_0,
                                         0.7692346550528415 * p_affine_6_1 + 0.23076534494715845 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.5 * p_affine_6_0 + 0.5 * p_affine_7_0,
                                         0.5 * p_affine_6_1 + 0.5 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.2307653449471585 * p_affine_6_0 + 0.7692346550528415 * p_affine_7_0,
                                         0.2307653449471585 * p_affine_6_1 + 0.7692346550528415 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.046910077030668074 * p_affine_6_0 + 0.95308992296933193 * p_affine_7_0,
                                         0.046910077030668074 * p_affine_6_1 + 0.95308992296933193 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
      real_t tmp_0   = -p_affine_3_0;
      real_t tmp_1   = p_affine_4_0 + tmp_0;
      real_t tmp_2   = p_affine_3_0 - p_affine_5_0;
      real_t tmp_3   = -p_affine_3_1;
      real_t tmp_4   = p_affine_5_1 + tmp_3;
      real_t tmp_5   = tmp_1 * tmp_4;
      real_t tmp_6   = p_affine_5_0 + tmp_0;
      real_t tmp_7   = p_affine_4_1 + tmp_3;
      real_t tmp_8   = 1.0 / ( tmp_5 - tmp_6 * tmp_7 );
      real_t tmp_9   = -p_affine_6_1 + p_affine_7_1;
      real_t tmp_10  = p_affine_6_1 + 0.046910077030668018 * tmp_9;
      real_t tmp_11  = tmp_8 * ( tmp_10 + tmp_3 );
      real_t tmp_12  = -p_affine_6_0 + p_affine_7_0;
      real_t tmp_13  = p_affine_6_0 + 0.046910077030668018 * tmp_12;
      real_t tmp_14  = tmp_8 * ( tmp_0 + tmp_13 );
      real_t tmp_15  = tmp_11 * tmp_2 + tmp_14 * tmp_4 - 1.0 / 3.0;
      real_t tmp_16  = p_affine_3_1 - p_affine_4_1;
      real_t tmp_17  = tmp_1 * tmp_11 + tmp_14 * tmp_16 - 1.0 / 3.0;
      real_t tmp_18  = tmp_1 * tmp_15 + tmp_17 * tmp_6;
      real_t tmp_19  = -p_affine_0_0;
      real_t tmp_20  = p_affine_1_0 + tmp_19;
      real_t tmp_21  = -p_affine_0_1;
      real_t tmp_22  = p_affine_2_1 + tmp_21;
      real_t tmp_23  = tmp_20 * tmp_22;
      real_t tmp_24  = p_affine_2_0 + tmp_19;
      real_t tmp_25  = p_affine_1_1 + tmp_21;
      real_t tmp_26  = 1.0 / ( tmp_23 - tmp_24 * tmp_25 );
      real_t tmp_27  = 1.0 * tmp_26;
      real_t tmp_28  = tmp_23 * tmp_27;
      real_t tmp_29  = p_affine_0_1 - p_affine_1_1;
      real_t tmp_30  = 0.5 * tmp_26;
      real_t tmp_31  = tmp_20 * tmp_30;
      real_t tmp_32  = p_affine_0_0 - p_affine_2_0;
      real_t tmp_33  = tmp_22 * tmp_30;
      real_t tmp_34  = tmp_24 * tmp_31 + tmp_25 * tmp_33 + tmp_29 * tmp_33 + tmp_31 * tmp_32;
      real_t tmp_35  = 0.5 * p_affine_10_0 * ( tmp_24 * tmp_27 * tmp_29 + tmp_28 ) + 0.5 * p_affine_10_1 * tmp_34;
      real_t tmp_36  = tmp_26 * ( tmp_10 + tmp_21 );
      real_t tmp_37  = tmp_26 * ( tmp_13 + tmp_19 );
      real_t tmp_38  = tmp_22 * tmp_37 + tmp_32 * tmp_36 - 1.0 / 3.0;
      real_t tmp_39  = tmp_20 * tmp_36 + tmp_29 * tmp_37 - 1.0 / 3.0;
      real_t tmp_40  = tmp_20 * tmp_38 + tmp_24 * tmp_39;
      real_t tmp_41  = 1.0 * tmp_8;
      real_t tmp_42  = tmp_41 * tmp_5;
      real_t tmp_43  = 0.5 * tmp_8;
      real_t tmp_44  = tmp_1 * tmp_43;
      real_t tmp_45  = tmp_4 * tmp_43;
      real_t tmp_46  = tmp_16 * tmp_45 + tmp_2 * tmp_44 + tmp_44 * tmp_6 + tmp_45 * tmp_7;
      real_t tmp_47  = 0.5 * p_affine_10_0 * ( tmp_16 * tmp_41 * tmp_6 + tmp_42 ) + 0.5 * p_affine_10_1 * tmp_46;
      real_t tmp_48  = tmp_15 * tmp_7 + tmp_17 * tmp_4;
      real_t tmp_49  = 0.5 * p_affine_10_0 * tmp_34 + 0.5 * p_affine_10_1 * ( tmp_25 * tmp_27 * tmp_32 + tmp_28 );
      real_t tmp_50  = tmp_22 * tmp_39 + tmp_25 * tmp_38;
      real_t tmp_51  = 0.5 * p_affine_10_0 * tmp_46 + 0.5 * p_affine_10_1 * ( tmp_2 * tmp_41 * tmp_7 + tmp_42 );
      real_t tmp_52  = std::abs( std::pow( ( tmp_12 * tmp_12 ) + ( tmp_9 * tmp_9 ), 1.0 / 2.0 ) );
      real_t tmp_53  = 6 / tmp_52;
      real_t tmp_54  = 2 * tmp_52;
      real_t tmp_55  = p_affine_6_1 + 0.23076534494715845 * tmp_9;
      real_t tmp_56  = tmp_8 * ( tmp_3 + tmp_55 );
      real_t tmp_57  = p_affine_6_0 + 0.23076534494715845 * tmp_12;
      real_t tmp_58  = tmp_8 * ( tmp_0 + tmp_57 );
      real_t tmp_59  = tmp_2 * tmp_56 + tmp_4 * tmp_58 - 1.0 / 3.0;
      real_t tmp_60  = tmp_1 * tmp_56 + tmp_16 * tmp_58 - 1.0 / 3.0;
      real_t tmp_61  = tmp_1 * tmp_59 + tmp_6 * tmp_60;
      real_t tmp_62  = tmp_26 * ( tmp_21 + tmp_55 );
      real_t tmp_63  = tmp_26 * ( tmp_19 + tmp_57 );
      real_t tmp_64  = tmp_22 * tmp_63 + tmp_32 * tmp_62 - 1.0 / 3.0;
      real_t tmp_65  = tmp_20 * tmp_62 + tmp_29 * tmp_63 - 1.0 / 3.0;
      real_t tmp_66  = tmp_20 * tmp_64 + tmp_24 * tmp_65;
      real_t tmp_67  = tmp_4 * tmp_60 + tmp_59 * tmp_7;
      real_t tmp_68  = tmp_22 * tmp_65 + tmp_25 * tmp_64;
      real_t tmp_69  = p_affine_6_1 + 0.5 * tmp_9;
      real_t tmp_70  = tmp_8 * ( tmp_3 + tmp_69 );
      real_t tmp_71  = p_affine_6_0 + 0.5 * tmp_12;
      real_t tmp_72  = tmp_8 * ( tmp_0 + tmp_71 );
      real_t tmp_73  = tmp_2 * tmp_70 + tmp_4 * tmp_72 - 1.0 / 3.0;
      real_t tmp_74  = tmp_1 * tmp_70 + tmp_16 * tmp_72 - 1.0 / 3.0;
      real_t tmp_75  = tmp_1 * tmp_73 + tmp_6 * tmp_74;
      real_t tmp_76  = tmp_26 * ( tmp_21 + tmp_69 );
      real_t tmp_77  = tmp_26 * ( tmp_19 + tmp_71 );
      real_t tmp_78  = tmp_22 * tmp_77 + tmp_32 * tmp_76 - 1.0 / 3.0;
      real_t tmp_79  = tmp_20 * tmp_76 + tmp_29 * tmp_77 - 1.0 / 3.0;
      real_t tmp_80  = tmp_20 * tmp_78 + tmp_24 * tmp_79;
      real_t tmp_81  = tmp_4 * tmp_74 + tmp_7 * tmp_73;
      real_t tmp_82  = tmp_22 * tmp_79 + tmp_25 * tmp_78;
      real_t tmp_83  = p_affine_6_1 + 0.7692346550528415 * tmp_9;
      real_t tmp_84  = tmp_8 * ( tmp_3 + tmp_83 );
      real_t tmp_85  = p_affine_6_0 + 0.7692346550528415 * tmp_12;
      real_t tmp_86  = tmp_8 * ( tmp_0 + tmp_85 );
      real_t tmp_87  = tmp_2 * tmp_84 + tmp_4 * tmp_86 - 1.0 / 3.0;
      real_t tmp_88  = tmp_1 * tmp_84 + tmp_16 * tmp_86 - 1.0 / 3.0;
      real_t tmp_89  = tmp_1 * tmp_87 + tmp_6 * tmp_88;
      real_t tmp_90  = tmp_26 * ( tmp_21 + tmp_83 );
      real_t tmp_91  = tmp_26 * ( tmp_19 + tmp_85 );
      real_t tmp_92  = tmp_22 * tmp_91 + tmp_32 * tmp_90 - 1.0 / 3.0;
      real_t tmp_93  = tmp_20 * tmp_90 + tmp_29 * tmp_91 - 1.0 / 3.0;
      real_t tmp_94  = tmp_20 * tmp_92 + tmp_24 * tmp_93;
      real_t tmp_95  = tmp_4 * tmp_88 + tmp_7 * tmp_87;
      real_t tmp_96  = tmp_22 * tmp_93 + tmp_25 * tmp_92;
      real_t tmp_97  = p_affine_6_1 + 0.95308992296933193 * tmp_9;
      real_t tmp_98  = tmp_8 * ( tmp_3 + tmp_97 );
      real_t tmp_99  = p_affine_6_0 + 0.95308992296933193 * tmp_12;
      real_t tmp_100 = tmp_8 * ( tmp_0 + tmp_99 );
      real_t tmp_101 = tmp_100 * tmp_4 + tmp_2 * tmp_98 - 1.0 / 3.0;
      real_t tmp_102 = tmp_1 * tmp_98 + tmp_100 * tmp_16 - 1.0 / 3.0;
      real_t tmp_103 = tmp_1 * tmp_101 + tmp_102 * tmp_6;
      real_t tmp_104 = tmp_26 * ( tmp_21 + tmp_97 );
      real_t tmp_105 = tmp_26 * ( tmp_19 + tmp_99 );
      real_t tmp_106 = tmp_104 * tmp_32 + tmp_105 * tmp_22 - 1.0 / 3.0;
      real_t tmp_107 = tmp_104 * tmp_20 + tmp_105 * tmp_29 - 1.0 / 3.0;
      real_t tmp_108 = tmp_106 * tmp_20 + tmp_107 * tmp_24;
      real_t tmp_109 = tmp_101 * tmp_7 + tmp_102 * tmp_4;
      real_t tmp_110 = tmp_106 * tmp_25 + tmp_107 * tmp_22;
      real_t a_0_0   = 0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id0 * tmp_54 *
                         ( tmp_18 * tmp_35 - tmp_40 * tmp_47 + tmp_48 * tmp_49 - tmp_50 * tmp_51 -
                           tmp_53 * ( tmp_18 * tmp_40 + tmp_48 * tmp_50 ) ) +
                     0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id1 * tmp_54 *
                         ( tmp_35 * tmp_61 - tmp_47 * tmp_66 + tmp_49 * tmp_67 - tmp_51 * tmp_68 -
                           tmp_53 * ( tmp_61 * tmp_66 + tmp_67 * tmp_68 ) ) +
                     0.2844444444444445 * Scalar_Variable_Coefficient_2D_mu_out0_id2 * tmp_54 *
                         ( tmp_35 * tmp_75 - tmp_47 * tmp_80 + tmp_49 * tmp_81 - tmp_51 * tmp_82 -
                           tmp_53 * ( tmp_75 * tmp_80 + tmp_81 * tmp_82 ) ) +
                     0.2393143352496831 * Scalar_Variable_Coefficient_2D_mu_out0_id3 * tmp_54 *
                         ( tmp_35 * tmp_89 - tmp_47 * tmp_94 + tmp_49 * tmp_95 - tmp_51 * tmp_96 -
                           tmp_53 * ( tmp_89 * tmp_94 + tmp_95 * tmp_96 ) ) +
                     0.11846344252809471 * Scalar_Variable_Coefficient_2D_mu_out0_id4 * tmp_54 *
                         ( tmp_103 * tmp_35 - tmp_108 * tmp_47 + tmp_109 * tmp_49 - tmp_110 * tmp_51 -
                           tmp_53 * ( tmp_103 * tmp_108 + tmp_109 * tmp_110 ) );
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

      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id4 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.95308992296933193 * p_affine_6_0 + 0.046910077030668018 * p_affine_7_0,
                                         0.95308992296933193 * p_affine_6_1 + 0.046910077030668018 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      Scalar_Variable_Coefficient_2D_mu( 0.7692346550528415 * p_affine_6_0 + 0.23076534494715845 * p_affine_7_0,
                                         0.7692346550528415 * p_affine_6_1 + 0.23076534494715845 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id1 );
      Scalar_Variable_Coefficient_2D_mu( 0.5 * p_affine_6_0 + 0.5 * p_affine_7_0,
                                         0.5 * p_affine_6_1 + 0.5 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id2 );
      Scalar_Variable_Coefficient_2D_mu( 0.2307653449471585 * p_affine_6_0 + 0.7692346550528415 * p_affine_7_0,
                                         0.2307653449471585 * p_affine_6_1 + 0.7692346550528415 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id3 );
      Scalar_Variable_Coefficient_2D_mu( 0.046910077030668074 * p_affine_6_0 + 0.95308992296933193 * p_affine_7_0,
                                         0.046910077030668074 * p_affine_6_1 + 0.95308992296933193 * p_affine_7_1,
                                         &Scalar_Variable_Coefficient_2D_mu_out0_id4 );
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
      real_t a_0_0  = 1.4215613103371365 * Scalar_Variable_Coefficient_2D_mu_out0_id0 *
                         ( ( ( tmp_1 * tmp_14 + tmp_16 * tmp_5 ) * ( tmp_1 * tmp_14 + tmp_16 * tmp_5 ) ) +
                           ( ( tmp_14 * tmp_6 + tmp_16 * tmp_4 ) * ( tmp_14 * tmp_6 + tmp_16 * tmp_4 ) ) ) +
                     2.8717720229961969 * Scalar_Variable_Coefficient_2D_mu_out0_id1 *
                         ( ( ( tmp_1 * tmp_19 + tmp_20 * tmp_5 ) * ( tmp_1 * tmp_19 + tmp_20 * tmp_5 ) ) +
                           ( ( tmp_19 * tmp_6 + tmp_20 * tmp_4 ) * ( tmp_19 * tmp_6 + tmp_20 * tmp_4 ) ) ) +
                     3.413333333333334 * Scalar_Variable_Coefficient_2D_mu_out0_id2 *
                         ( ( ( tmp_1 * tmp_23 + tmp_24 * tmp_5 ) * ( tmp_1 * tmp_23 + tmp_24 * tmp_5 ) ) +
                           ( ( tmp_23 * tmp_6 + tmp_24 * tmp_4 ) * ( tmp_23 * tmp_6 + tmp_24 * tmp_4 ) ) ) +
                     2.8717720229961969 * Scalar_Variable_Coefficient_2D_mu_out0_id3 *
                         ( ( ( tmp_1 * tmp_27 + tmp_28 * tmp_5 ) * ( tmp_1 * tmp_27 + tmp_28 * tmp_5 ) ) +
                           ( ( tmp_27 * tmp_6 + tmp_28 * tmp_4 ) * ( tmp_27 * tmp_6 + tmp_28 * tmp_4 ) ) ) +
                     1.4215613103371365 * Scalar_Variable_Coefficient_2D_mu_out0_id4 *
                         ( ( ( tmp_1 * tmp_31 + tmp_32 * tmp_5 ) * ( tmp_1 * tmp_31 + tmp_32 * tmp_5 ) ) +
                           ( ( tmp_31 * tmp_6 + tmp_32 * tmp_4 ) * ( tmp_31 * tmp_6 + tmp_32 * tmp_4 ) ) );
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

 private:
   void Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t* out_0 ) const
   {
      *out_0 = callback2D( Point3D( { in_0, in_1 } ) );
   }
   std::function< real_t( const Point3D& ) > callback2D;
};

} // namespace eg
} // namespace dg
} // namespace hyteg

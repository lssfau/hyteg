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

/*
 * The entire file was generated with the HyTeG form generator.
 *
 * Software:
 *
 * - quadpy version: 0.16.5
 *
 * Avoid modifying this file. If buggy, consider fixing the generator itself.
 */

#include "p1_invk_mass_affine_q4.hpp"

namespace hyteg {
namespace forms {

   void p1_invk_mass_affine_q4::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Scalar_Variable_Coefficient_2D_k_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id5 = 0;
      Scalar_Variable_Coefficient_2D_k( 0.091576213509770743*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.81684757298045851*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.81684757298045851*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id0 );
      Scalar_Variable_Coefficient_2D_k( 0.44594849091596489*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.10810301816807022*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.10810301816807022*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id1 );
      Scalar_Variable_Coefficient_2D_k( 0.091576213509770743*p_affine_0_0 + 0.81684757298045851*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.81684757298045851*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id2 );
      Scalar_Variable_Coefficient_2D_k( 0.44594849091596489*p_affine_0_0 + 0.10810301816807022*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.10810301816807022*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id3 );
      Scalar_Variable_Coefficient_2D_k( 0.81684757298045851*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.81684757298045851*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id4 );
      Scalar_Variable_Coefficient_2D_k( 0.10810301816807022*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.10810301816807022*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id5 );
      real_t tmp_0 = 0.091576213509770743;
      real_t tmp_1 = 1.0*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_2 = 0.054975871827660928*tmp_1/Scalar_Variable_Coefficient_2D_k_out0_id0;
      real_t tmp_3 = 0.44594849091596489;
      real_t tmp_4 = 0.11169079483900572*tmp_1/Scalar_Variable_Coefficient_2D_k_out0_id1;
      real_t tmp_5 = 0.091576213509770743;
      real_t tmp_6 = 0.054975871827660928*tmp_1/Scalar_Variable_Coefficient_2D_k_out0_id2;
      real_t tmp_7 = 0.44594849091596489;
      real_t tmp_8 = 0.11169079483900572*tmp_1/Scalar_Variable_Coefficient_2D_k_out0_id3;
      real_t tmp_9 = 0.81684757298045851;
      real_t tmp_10 = 0.054975871827660928*tmp_1/Scalar_Variable_Coefficient_2D_k_out0_id4;
      real_t tmp_11 = 0.10810301816807022;
      real_t tmp_12 = 0.11169079483900572*tmp_1/Scalar_Variable_Coefficient_2D_k_out0_id5;
      real_t tmp_13 = tmp_0*tmp_2;
      real_t tmp_14 = tmp_3*tmp_4;
      real_t tmp_15 = tmp_5*tmp_6;
      real_t tmp_16 = tmp_7*tmp_8;
      real_t tmp_17 = tmp_10*tmp_9;
      real_t tmp_18 = tmp_11*tmp_12;
      real_t tmp_19 = 0.091576213509770743*tmp_13 + 0.44594849091596489*tmp_14 + 0.81684757298045851*tmp_15 + 0.10810301816807022*tmp_16 + 0.091576213509770743*tmp_17 + 0.44594849091596489*tmp_18;
      real_t tmp_20 = 0.81684757298045851*tmp_13 + 0.10810301816807022*tmp_14 + 0.091576213509770743*tmp_15 + 0.44594849091596489*tmp_16 + 0.091576213509770743*tmp_17 + 0.44594849091596489*tmp_18;
      real_t tmp_21 = 0.0083862028807871174*tmp_10 + 0.19887005655022641*tmp_12 + 0.074803807748196505*tmp_2 + 0.04820837781551205*tmp_4 + 0.074803807748196505*tmp_6 + 0.04820837781551205*tmp_8;
      real_t a_0_0 = (tmp_0*tmp_0)*tmp_2 + tmp_10*(tmp_9*tmp_9) + (tmp_11*tmp_11)*tmp_12 + (tmp_3*tmp_3)*tmp_4 + (tmp_5*tmp_5)*tmp_6 + (tmp_7*tmp_7)*tmp_8;
      real_t a_0_1 = tmp_19;
      real_t a_0_2 = tmp_20;
      real_t a_1_0 = tmp_19;
      real_t a_1_1 = 0.0083862028807871174*tmp_10 + 0.19887005655022641*tmp_12 + 0.0083862028807871174*tmp_2 + 0.19887005655022641*tmp_4 + 0.6672399574840655*tmp_6 + 0.01168626253704612*tmp_8;
      real_t a_1_2 = tmp_21;
      real_t a_2_0 = tmp_20;
      real_t a_2_1 = tmp_21;
      real_t a_2_2 = 0.0083862028807871174*tmp_10 + 0.19887005655022641*tmp_12 + 0.6672399574840655*tmp_2 + 0.01168626253704612*tmp_4 + 0.0083862028807871174*tmp_6 + 0.19887005655022641*tmp_8;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(1, 0)) = a_1_0;
      (elMat(1, 1)) = a_1_1;
      (elMat(1, 2)) = a_1_2;
      (elMat(2, 0)) = a_2_0;
      (elMat(2, 1)) = a_2_1;
      (elMat(2, 2)) = a_2_2;
   }

   void p1_invk_mass_affine_q4::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Scalar_Variable_Coefficient_2D_k_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id5 = 0;
      Scalar_Variable_Coefficient_2D_k( 0.091576213509770743*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.81684757298045851*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.81684757298045851*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id0 );
      Scalar_Variable_Coefficient_2D_k( 0.44594849091596489*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.10810301816807022*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.10810301816807022*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id1 );
      Scalar_Variable_Coefficient_2D_k( 0.091576213509770743*p_affine_0_0 + 0.81684757298045851*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.81684757298045851*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id2 );
      Scalar_Variable_Coefficient_2D_k( 0.44594849091596489*p_affine_0_0 + 0.10810301816807022*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.10810301816807022*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id3 );
      Scalar_Variable_Coefficient_2D_k( 0.81684757298045851*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.81684757298045851*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id4 );
      Scalar_Variable_Coefficient_2D_k( 0.10810301816807022*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.10810301816807022*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id5 );
      real_t tmp_0 = 0.091576213509770743;
      real_t tmp_1 = 1.0*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_2 = 0.054975871827660928*tmp_1/Scalar_Variable_Coefficient_2D_k_out0_id0;
      real_t tmp_3 = 0.44594849091596489;
      real_t tmp_4 = 0.11169079483900572*tmp_1/Scalar_Variable_Coefficient_2D_k_out0_id1;
      real_t tmp_5 = 0.091576213509770743;
      real_t tmp_6 = 0.054975871827660928*tmp_1/Scalar_Variable_Coefficient_2D_k_out0_id2;
      real_t tmp_7 = 0.44594849091596489;
      real_t tmp_8 = 0.11169079483900572*tmp_1/Scalar_Variable_Coefficient_2D_k_out0_id3;
      real_t tmp_9 = 0.81684757298045851;
      real_t tmp_10 = 0.054975871827660928*tmp_1/Scalar_Variable_Coefficient_2D_k_out0_id4;
      real_t tmp_11 = 0.10810301816807022;
      real_t tmp_12 = 0.11169079483900572*tmp_1/Scalar_Variable_Coefficient_2D_k_out0_id5;
      real_t tmp_13 = tmp_0*tmp_2;
      real_t tmp_14 = tmp_3*tmp_4;
      real_t tmp_15 = tmp_5*tmp_6;
      real_t tmp_16 = tmp_7*tmp_8;
      real_t tmp_17 = tmp_10*tmp_9;
      real_t tmp_18 = tmp_11*tmp_12;
      real_t a_0_0 = (tmp_0*tmp_0)*tmp_2 + tmp_10*(tmp_9*tmp_9) + (tmp_11*tmp_11)*tmp_12 + (tmp_3*tmp_3)*tmp_4 + (tmp_5*tmp_5)*tmp_6 + (tmp_7*tmp_7)*tmp_8;
      real_t a_0_1 = 0.091576213509770743*tmp_13 + 0.44594849091596489*tmp_14 + 0.81684757298045851*tmp_15 + 0.10810301816807022*tmp_16 + 0.091576213509770743*tmp_17 + 0.44594849091596489*tmp_18;
      real_t a_0_2 = 0.81684757298045851*tmp_13 + 0.10810301816807022*tmp_14 + 0.091576213509770743*tmp_15 + 0.44594849091596489*tmp_16 + 0.091576213509770743*tmp_17 + 0.44594849091596489*tmp_18;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_invk_mass_affine_q4::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_0_2 = coords[0][2];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_1_2 = coords[1][2];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t p_affine_2_2 = coords[2][2];
      real_t p_affine_3_0 = coords[3][0];
      real_t p_affine_3_1 = coords[3][1];
      real_t p_affine_3_2 = coords[3][2];
      real_t Scalar_Variable_Coefficient_3D_k_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id13 = 0;
      Scalar_Variable_Coefficient_3D_k( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.067342242210098213*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.067342242210098213*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.067342242210098213*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id0 );
      Scalar_Variable_Coefficient_3D_k( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.72179424906732625*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.72179424906732625*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.72179424906732625*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id1 );
      Scalar_Variable_Coefficient_3D_k( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id2 );
      Scalar_Variable_Coefficient_3D_k( 0.045503704125649629*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.045503704125649629*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.045503704125649629*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id3 );
      Scalar_Variable_Coefficient_3D_k( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id4 );
      Scalar_Variable_Coefficient_3D_k( 0.45449629587435036*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id5 );
      Scalar_Variable_Coefficient_3D_k( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.067342242210098213*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.067342242210098213*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.067342242210098213*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id6 );
      Scalar_Variable_Coefficient_3D_k( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.72179424906732625*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.72179424906732625*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.72179424906732625*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id7 );
      Scalar_Variable_Coefficient_3D_k( 0.3108859192633005*p_affine_0_0 + 0.067342242210098213*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.067342242210098213*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.067342242210098213*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id8 );
      Scalar_Variable_Coefficient_3D_k( 0.092735250310891248*p_affine_0_0 + 0.72179424906732625*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.72179424906732625*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.72179424906732625*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id9 );
      Scalar_Variable_Coefficient_3D_k( 0.067342242210098102*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.067342242210098102*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.067342242210098102*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id10 );
      Scalar_Variable_Coefficient_3D_k( 0.72179424906732625*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.72179424906732625*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.72179424906732625*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id11 );
      Scalar_Variable_Coefficient_3D_k( 0.045503704125649636*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id12 );
      Scalar_Variable_Coefficient_3D_k( 0.045503704125649636*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id13 );
      real_t tmp_0 = 0.3108859192633005;
      real_t tmp_1 = p_affine_0_0*p_affine_1_1;
      real_t tmp_2 = p_affine_0_0*p_affine_1_2;
      real_t tmp_3 = p_affine_2_1*p_affine_3_2;
      real_t tmp_4 = p_affine_0_1*p_affine_1_0;
      real_t tmp_5 = p_affine_0_1*p_affine_1_2;
      real_t tmp_6 = p_affine_2_2*p_affine_3_0;
      real_t tmp_7 = p_affine_0_2*p_affine_1_0;
      real_t tmp_8 = p_affine_0_2*p_affine_1_1;
      real_t tmp_9 = p_affine_2_0*p_affine_3_1;
      real_t tmp_10 = p_affine_2_2*p_affine_3_1;
      real_t tmp_11 = p_affine_2_0*p_affine_3_2;
      real_t tmp_12 = p_affine_2_1*p_affine_3_0;
      real_t tmp_13 = 1.0*std::abs(p_affine_0_0*tmp_10 - p_affine_0_0*tmp_3 + p_affine_0_1*tmp_11 - p_affine_0_1*tmp_6 + p_affine_0_2*tmp_12 - p_affine_0_2*tmp_9 - p_affine_1_0*tmp_10 + p_affine_1_0*tmp_3 - p_affine_1_1*tmp_11 + p_affine_1_1*tmp_6 - p_affine_1_2*tmp_12 + p_affine_1_2*tmp_9 - p_affine_2_0*tmp_5 + p_affine_2_0*tmp_8 + p_affine_2_1*tmp_2 - p_affine_2_1*tmp_7 - p_affine_2_2*tmp_1 + p_affine_2_2*tmp_4 + p_affine_3_0*tmp_5 - p_affine_3_0*tmp_8 - p_affine_3_1*tmp_2 + p_affine_3_1*tmp_7 + p_affine_3_2*tmp_1 - p_affine_3_2*tmp_4);
      real_t tmp_14 = 0.018781320953002646*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id0;
      real_t tmp_15 = 0.092735250310891248;
      real_t tmp_16 = 0.012248840519393657*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id1;
      real_t tmp_17 = 0.067342242210098102;
      real_t tmp_18 = 0.018781320953002646*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id10;
      real_t tmp_19 = 0.72179424906732625;
      real_t tmp_20 = 0.012248840519393657*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id11;
      real_t tmp_21 = 0.045503704125649636;
      real_t tmp_22 = 0.0070910034628469103*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id12;
      real_t tmp_23 = 0.045503704125649636;
      real_t tmp_24 = 0.0070910034628469103*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id13;
      real_t tmp_25 = 0.45449629587435036;
      real_t tmp_26 = 0.0070910034628469103*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id2;
      real_t tmp_27 = 0.045503704125649629;
      real_t tmp_28 = 0.0070910034628469103*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id3;
      real_t tmp_29 = 0.45449629587435036;
      real_t tmp_30 = 0.0070910034628469103*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id4;
      real_t tmp_31 = 0.45449629587435036;
      real_t tmp_32 = 0.0070910034628469103*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id5;
      real_t tmp_33 = 0.3108859192633005;
      real_t tmp_34 = 0.018781320953002646*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id6;
      real_t tmp_35 = 0.092735250310891248;
      real_t tmp_36 = 0.012248840519393657*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id7;
      real_t tmp_37 = 0.3108859192633005;
      real_t tmp_38 = 0.018781320953002646*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id8;
      real_t tmp_39 = 0.092735250310891248;
      real_t tmp_40 = 0.012248840519393657*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id9;
      real_t tmp_41 = tmp_0*tmp_14;
      real_t tmp_42 = tmp_15*tmp_16;
      real_t tmp_43 = tmp_17*tmp_18;
      real_t tmp_44 = tmp_19*tmp_20;
      real_t tmp_45 = tmp_21*tmp_22;
      real_t tmp_46 = tmp_23*tmp_24;
      real_t tmp_47 = tmp_25*tmp_26;
      real_t tmp_48 = tmp_27*tmp_28;
      real_t tmp_49 = tmp_29*tmp_30;
      real_t tmp_50 = tmp_31*tmp_32;
      real_t tmp_51 = tmp_33*tmp_34;
      real_t tmp_52 = tmp_35*tmp_36;
      real_t tmp_53 = tmp_37*tmp_38;
      real_t tmp_54 = tmp_39*tmp_40;
      real_t tmp_55 = 0.31088591926330061*tmp_41 + 0.092735250310891248*tmp_42 + 0.31088591926330061*tmp_43 + 0.092735250310891248*tmp_44 + 0.045503704125649642*tmp_45 + 0.45449629587435036*tmp_46 + 0.045503704125649642*tmp_47 + 0.45449629587435036*tmp_48 + 0.045503704125649642*tmp_49 + 0.45449629587435036*tmp_50 + 0.31088591926330061*tmp_51 + 0.092735250310891248*tmp_52 + 0.067342242210098213*tmp_53 + 0.72179424906732625*tmp_54;
      real_t tmp_56 = 0.31088591926330061*tmp_41 + 0.092735250310891248*tmp_42 + 0.31088591926330061*tmp_43 + 0.092735250310891248*tmp_44 + 0.45449629587435036*tmp_45 + 0.045503704125649642*tmp_46 + 0.045503704125649642*tmp_47 + 0.45449629587435036*tmp_48 + 0.45449629587435036*tmp_49 + 0.045503704125649642*tmp_50 + 0.067342242210098213*tmp_51 + 0.72179424906732625*tmp_52 + 0.31088591926330061*tmp_53 + 0.092735250310891248*tmp_54;
      real_t tmp_57 = 0.067342242210098213*tmp_41 + 0.72179424906732625*tmp_42 + 0.31088591926330061*tmp_43 + 0.092735250310891248*tmp_44 + 0.45449629587435036*tmp_45 + 0.45449629587435036*tmp_46 + 0.45449629587435036*tmp_47 + 0.045503704125649642*tmp_48 + 0.045503704125649642*tmp_49 + 0.045503704125649642*tmp_50 + 0.31088591926330061*tmp_51 + 0.092735250310891248*tmp_52 + 0.31088591926330061*tmp_53 + 0.092735250310891248*tmp_54;
      real_t tmp_58 = 0.31088591926330061*tmp_14;
      real_t tmp_59 = 0.092735250310891248*tmp_16;
      real_t tmp_60 = 0.31088591926330061*tmp_18;
      real_t tmp_61 = 0.092735250310891248*tmp_20;
      real_t tmp_62 = 0.045503704125649642*tmp_22;
      real_t tmp_63 = 0.45449629587435036*tmp_24;
      real_t tmp_64 = 0.045503704125649642*tmp_26;
      real_t tmp_65 = 0.45449629587435036*tmp_28;
      real_t tmp_66 = 0.045503704125649642*tmp_30;
      real_t tmp_67 = 0.45449629587435036*tmp_32;
      real_t tmp_68 = 0.31088591926330061*tmp_34;
      real_t tmp_69 = 0.092735250310891248*tmp_36;
      real_t tmp_70 = 0.067342242210098213*tmp_38;
      real_t tmp_71 = 0.72179424906732625*tmp_40;
      real_t tmp_72 = 0.31088591926330061*tmp_58 + 0.092735250310891248*tmp_59 + 0.31088591926330061*tmp_60 + 0.092735250310891248*tmp_61 + 0.45449629587435036*tmp_62 + 0.045503704125649642*tmp_63 + 0.045503704125649642*tmp_64 + 0.45449629587435036*tmp_65 + 0.45449629587435036*tmp_66 + 0.045503704125649642*tmp_67 + 0.067342242210098213*tmp_68 + 0.72179424906732625*tmp_69 + 0.31088591926330061*tmp_70 + 0.092735250310891248*tmp_71;
      real_t tmp_73 = 0.067342242210098213*tmp_58 + 0.72179424906732625*tmp_59 + 0.31088591926330061*tmp_60 + 0.092735250310891248*tmp_61 + 0.45449629587435036*tmp_62 + 0.45449629587435036*tmp_63 + 0.45449629587435036*tmp_64 + 0.045503704125649642*tmp_65 + 0.045503704125649642*tmp_66 + 0.045503704125649642*tmp_67 + 0.31088591926330061*tmp_68 + 0.092735250310891248*tmp_69 + 0.31088591926330061*tmp_70 + 0.092735250310891248*tmp_71;
      real_t tmp_74 = 0.020935754874738227*tmp_14 + 0.066935770360220276*tmp_16 + 0.096650054796187462*tmp_18 + 0.0085998266502236558*tmp_20 + 0.20656688296350503*tmp_22 + 0.020681264973670156*tmp_24 + 0.020681264973670156*tmp_26 + 0.020681264973670156*tmp_28 + 0.020681264973670156*tmp_30 + 0.0020705870891546642*tmp_32 + 0.020935754874738227*tmp_34 + 0.066935770360220276*tmp_36 + 0.096650054796187462*tmp_38 + 0.0085998266502236558*tmp_40;
      real_t a_0_0 = (tmp_0*tmp_0)*tmp_14 + (tmp_15*tmp_15)*tmp_16 + (tmp_17*tmp_17)*tmp_18 + (tmp_19*tmp_19)*tmp_20 + (tmp_21*tmp_21)*tmp_22 + (tmp_23*tmp_23)*tmp_24 + (tmp_25*tmp_25)*tmp_26 + (tmp_27*tmp_27)*tmp_28 + (tmp_29*tmp_29)*tmp_30 + (tmp_31*tmp_31)*tmp_32 + (tmp_33*tmp_33)*tmp_34 + (tmp_35*tmp_35)*tmp_36 + (tmp_37*tmp_37)*tmp_38 + (tmp_39*tmp_39)*tmp_40;
      real_t a_0_1 = tmp_55;
      real_t a_0_2 = tmp_56;
      real_t a_0_3 = tmp_57;
      real_t a_1_0 = tmp_55;
      real_t a_1_1 = 0.096650054796187462*tmp_14 + 0.0085998266502236558*tmp_16 + 0.096650054796187462*tmp_18 + 0.0085998266502236558*tmp_20 + 0.0020705870891546642*tmp_22 + 0.20656688296350503*tmp_24 + 0.0020705870891546642*tmp_26 + 0.20656688296350503*tmp_28 + 0.0020705870891546642*tmp_30 + 0.20656688296350503*tmp_32 + 0.096650054796187462*tmp_34 + 0.0085998266502236558*tmp_36 + 0.0045349775858835335*tmp_38 + 0.5209869379866654*tmp_40;
      real_t a_1_2 = tmp_72;
      real_t a_1_3 = tmp_73;
      real_t a_2_0 = tmp_56;
      real_t a_2_1 = tmp_72;
      real_t a_2_2 = 0.096650054796187462*tmp_14 + 0.0085998266502236558*tmp_16 + 0.096650054796187462*tmp_18 + 0.0085998266502236558*tmp_20 + 0.20656688296350503*tmp_22 + 0.0020705870891546642*tmp_24 + 0.0020705870891546642*tmp_26 + 0.20656688296350503*tmp_28 + 0.20656688296350503*tmp_30 + 0.0020705870891546642*tmp_32 + 0.0045349775858835335*tmp_34 + 0.5209869379866654*tmp_36 + 0.096650054796187462*tmp_38 + 0.0085998266502236558*tmp_40;
      real_t a_2_3 = tmp_74;
      real_t a_3_0 = tmp_57;
      real_t a_3_1 = tmp_73;
      real_t a_3_2 = tmp_74;
      real_t a_3_3 = 0.0045349775858835335*tmp_14 + 0.5209869379866654*tmp_16 + 0.096650054796187462*tmp_18 + 0.0085998266502236558*tmp_20 + 0.20656688296350503*tmp_22 + 0.20656688296350503*tmp_24 + 0.20656688296350503*tmp_26 + 0.0020705870891546642*tmp_28 + 0.0020705870891546642*tmp_30 + 0.0020705870891546642*tmp_32 + 0.096650054796187462*tmp_34 + 0.0085998266502236558*tmp_36 + 0.096650054796187462*tmp_38 + 0.0085998266502236558*tmp_40;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
      (elMat(1, 0)) = a_1_0;
      (elMat(1, 1)) = a_1_1;
      (elMat(1, 2)) = a_1_2;
      (elMat(1, 3)) = a_1_3;
      (elMat(2, 0)) = a_2_0;
      (elMat(2, 1)) = a_2_1;
      (elMat(2, 2)) = a_2_2;
      (elMat(2, 3)) = a_2_3;
      (elMat(3, 0)) = a_3_0;
      (elMat(3, 1)) = a_3_1;
      (elMat(3, 2)) = a_3_2;
      (elMat(3, 3)) = a_3_3;
   }

   void p1_invk_mass_affine_q4::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_0_2 = coords[0][2];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_1_2 = coords[1][2];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t p_affine_2_2 = coords[2][2];
      real_t p_affine_3_0 = coords[3][0];
      real_t p_affine_3_1 = coords[3][1];
      real_t p_affine_3_2 = coords[3][2];
      real_t Scalar_Variable_Coefficient_3D_k_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id11 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id12 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id13 = 0;
      Scalar_Variable_Coefficient_3D_k( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.067342242210098213*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.067342242210098213*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.067342242210098213*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id0 );
      Scalar_Variable_Coefficient_3D_k( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.72179424906732625*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.72179424906732625*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.72179424906732625*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id1 );
      Scalar_Variable_Coefficient_3D_k( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id2 );
      Scalar_Variable_Coefficient_3D_k( 0.045503704125649629*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.045503704125649629*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.045503704125649629*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id3 );
      Scalar_Variable_Coefficient_3D_k( 0.45449629587435036*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id4 );
      Scalar_Variable_Coefficient_3D_k( 0.45449629587435036*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.045503704125649642*p_affine_3_0, 0.45449629587435036*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.045503704125649642*p_affine_3_1, 0.45449629587435036*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.045503704125649642*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id5 );
      Scalar_Variable_Coefficient_3D_k( 0.3108859192633005*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.067342242210098213*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.067342242210098213*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.067342242210098213*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id6 );
      Scalar_Variable_Coefficient_3D_k( 0.092735250310891248*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.72179424906732625*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.72179424906732625*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.72179424906732625*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id7 );
      Scalar_Variable_Coefficient_3D_k( 0.3108859192633005*p_affine_0_0 + 0.067342242210098213*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.3108859192633005*p_affine_0_1 + 0.067342242210098213*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.3108859192633005*p_affine_0_2 + 0.067342242210098213*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id8 );
      Scalar_Variable_Coefficient_3D_k( 0.092735250310891248*p_affine_0_0 + 0.72179424906732625*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.092735250310891248*p_affine_0_1 + 0.72179424906732625*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.092735250310891248*p_affine_0_2 + 0.72179424906732625*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id9 );
      Scalar_Variable_Coefficient_3D_k( 0.067342242210098102*p_affine_0_0 + 0.31088591926330061*p_affine_1_0 + 0.31088591926330061*p_affine_2_0 + 0.31088591926330061*p_affine_3_0, 0.067342242210098102*p_affine_0_1 + 0.31088591926330061*p_affine_1_1 + 0.31088591926330061*p_affine_2_1 + 0.31088591926330061*p_affine_3_1, 0.067342242210098102*p_affine_0_2 + 0.31088591926330061*p_affine_1_2 + 0.31088591926330061*p_affine_2_2 + 0.31088591926330061*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id10 );
      Scalar_Variable_Coefficient_3D_k( 0.72179424906732625*p_affine_0_0 + 0.092735250310891248*p_affine_1_0 + 0.092735250310891248*p_affine_2_0 + 0.092735250310891248*p_affine_3_0, 0.72179424906732625*p_affine_0_1 + 0.092735250310891248*p_affine_1_1 + 0.092735250310891248*p_affine_2_1 + 0.092735250310891248*p_affine_3_1, 0.72179424906732625*p_affine_0_2 + 0.092735250310891248*p_affine_1_2 + 0.092735250310891248*p_affine_2_2 + 0.092735250310891248*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id11 );
      Scalar_Variable_Coefficient_3D_k( 0.045503704125649636*p_affine_0_0 + 0.045503704125649642*p_affine_1_0 + 0.45449629587435036*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.045503704125649642*p_affine_1_1 + 0.45449629587435036*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.045503704125649642*p_affine_1_2 + 0.45449629587435036*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id12 );
      Scalar_Variable_Coefficient_3D_k( 0.045503704125649636*p_affine_0_0 + 0.45449629587435036*p_affine_1_0 + 0.045503704125649642*p_affine_2_0 + 0.45449629587435036*p_affine_3_0, 0.045503704125649636*p_affine_0_1 + 0.45449629587435036*p_affine_1_1 + 0.045503704125649642*p_affine_2_1 + 0.45449629587435036*p_affine_3_1, 0.045503704125649636*p_affine_0_2 + 0.45449629587435036*p_affine_1_2 + 0.045503704125649642*p_affine_2_2 + 0.45449629587435036*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id13 );
      real_t tmp_0 = 0.3108859192633005;
      real_t tmp_1 = p_affine_0_0*p_affine_1_1;
      real_t tmp_2 = p_affine_0_0*p_affine_1_2;
      real_t tmp_3 = p_affine_2_1*p_affine_3_2;
      real_t tmp_4 = p_affine_0_1*p_affine_1_0;
      real_t tmp_5 = p_affine_0_1*p_affine_1_2;
      real_t tmp_6 = p_affine_2_2*p_affine_3_0;
      real_t tmp_7 = p_affine_0_2*p_affine_1_0;
      real_t tmp_8 = p_affine_0_2*p_affine_1_1;
      real_t tmp_9 = p_affine_2_0*p_affine_3_1;
      real_t tmp_10 = p_affine_2_2*p_affine_3_1;
      real_t tmp_11 = p_affine_2_0*p_affine_3_2;
      real_t tmp_12 = p_affine_2_1*p_affine_3_0;
      real_t tmp_13 = 1.0*std::abs(p_affine_0_0*tmp_10 - p_affine_0_0*tmp_3 + p_affine_0_1*tmp_11 - p_affine_0_1*tmp_6 + p_affine_0_2*tmp_12 - p_affine_0_2*tmp_9 - p_affine_1_0*tmp_10 + p_affine_1_0*tmp_3 - p_affine_1_1*tmp_11 + p_affine_1_1*tmp_6 - p_affine_1_2*tmp_12 + p_affine_1_2*tmp_9 - p_affine_2_0*tmp_5 + p_affine_2_0*tmp_8 + p_affine_2_1*tmp_2 - p_affine_2_1*tmp_7 - p_affine_2_2*tmp_1 + p_affine_2_2*tmp_4 + p_affine_3_0*tmp_5 - p_affine_3_0*tmp_8 - p_affine_3_1*tmp_2 + p_affine_3_1*tmp_7 + p_affine_3_2*tmp_1 - p_affine_3_2*tmp_4);
      real_t tmp_14 = 0.018781320953002646*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id0;
      real_t tmp_15 = 0.092735250310891248;
      real_t tmp_16 = 0.012248840519393657*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id1;
      real_t tmp_17 = 0.067342242210098102;
      real_t tmp_18 = 0.018781320953002646*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id10;
      real_t tmp_19 = 0.72179424906732625;
      real_t tmp_20 = 0.012248840519393657*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id11;
      real_t tmp_21 = 0.045503704125649636;
      real_t tmp_22 = 0.0070910034628469103*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id12;
      real_t tmp_23 = 0.045503704125649636;
      real_t tmp_24 = 0.0070910034628469103*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id13;
      real_t tmp_25 = 0.45449629587435036;
      real_t tmp_26 = 0.0070910034628469103*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id2;
      real_t tmp_27 = 0.045503704125649629;
      real_t tmp_28 = 0.0070910034628469103*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id3;
      real_t tmp_29 = 0.45449629587435036;
      real_t tmp_30 = 0.0070910034628469103*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id4;
      real_t tmp_31 = 0.45449629587435036;
      real_t tmp_32 = 0.0070910034628469103*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id5;
      real_t tmp_33 = 0.3108859192633005;
      real_t tmp_34 = 0.018781320953002646*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id6;
      real_t tmp_35 = 0.092735250310891248;
      real_t tmp_36 = 0.012248840519393657*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id7;
      real_t tmp_37 = 0.3108859192633005;
      real_t tmp_38 = 0.018781320953002646*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id8;
      real_t tmp_39 = 0.092735250310891248;
      real_t tmp_40 = 0.012248840519393657*tmp_13/Scalar_Variable_Coefficient_3D_k_out0_id9;
      real_t tmp_41 = tmp_0*tmp_14;
      real_t tmp_42 = tmp_15*tmp_16;
      real_t tmp_43 = tmp_17*tmp_18;
      real_t tmp_44 = tmp_19*tmp_20;
      real_t tmp_45 = tmp_21*tmp_22;
      real_t tmp_46 = tmp_23*tmp_24;
      real_t tmp_47 = tmp_25*tmp_26;
      real_t tmp_48 = tmp_27*tmp_28;
      real_t tmp_49 = tmp_29*tmp_30;
      real_t tmp_50 = tmp_31*tmp_32;
      real_t tmp_51 = tmp_33*tmp_34;
      real_t tmp_52 = tmp_35*tmp_36;
      real_t tmp_53 = tmp_37*tmp_38;
      real_t tmp_54 = tmp_39*tmp_40;
      real_t a_0_0 = (tmp_0*tmp_0)*tmp_14 + (tmp_15*tmp_15)*tmp_16 + (tmp_17*tmp_17)*tmp_18 + (tmp_19*tmp_19)*tmp_20 + (tmp_21*tmp_21)*tmp_22 + (tmp_23*tmp_23)*tmp_24 + (tmp_25*tmp_25)*tmp_26 + (tmp_27*tmp_27)*tmp_28 + (tmp_29*tmp_29)*tmp_30 + (tmp_31*tmp_31)*tmp_32 + (tmp_33*tmp_33)*tmp_34 + (tmp_35*tmp_35)*tmp_36 + (tmp_37*tmp_37)*tmp_38 + (tmp_39*tmp_39)*tmp_40;
      real_t a_0_1 = 0.31088591926330061*tmp_41 + 0.092735250310891248*tmp_42 + 0.31088591926330061*tmp_43 + 0.092735250310891248*tmp_44 + 0.045503704125649642*tmp_45 + 0.45449629587435036*tmp_46 + 0.045503704125649642*tmp_47 + 0.45449629587435036*tmp_48 + 0.045503704125649642*tmp_49 + 0.45449629587435036*tmp_50 + 0.31088591926330061*tmp_51 + 0.092735250310891248*tmp_52 + 0.067342242210098213*tmp_53 + 0.72179424906732625*tmp_54;
      real_t a_0_2 = 0.31088591926330061*tmp_41 + 0.092735250310891248*tmp_42 + 0.31088591926330061*tmp_43 + 0.092735250310891248*tmp_44 + 0.45449629587435036*tmp_45 + 0.045503704125649642*tmp_46 + 0.045503704125649642*tmp_47 + 0.45449629587435036*tmp_48 + 0.45449629587435036*tmp_49 + 0.045503704125649642*tmp_50 + 0.067342242210098213*tmp_51 + 0.72179424906732625*tmp_52 + 0.31088591926330061*tmp_53 + 0.092735250310891248*tmp_54;
      real_t a_0_3 = 0.067342242210098213*tmp_41 + 0.72179424906732625*tmp_42 + 0.31088591926330061*tmp_43 + 0.092735250310891248*tmp_44 + 0.45449629587435036*tmp_45 + 0.45449629587435036*tmp_46 + 0.45449629587435036*tmp_47 + 0.045503704125649642*tmp_48 + 0.045503704125649642*tmp_49 + 0.045503704125649642*tmp_50 + 0.31088591926330061*tmp_51 + 0.092735250310891248*tmp_52 + 0.31088591926330061*tmp_53 + 0.092735250310891248*tmp_54;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_invk_mass_affine_q4::Scalar_Variable_Coefficient_2D_k( real_t in_0, real_t in_1, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_k( Point3D( {in_0, in_1, 0} ) );
   }

   void p1_invk_mass_affine_q4::Scalar_Variable_Coefficient_3D_k( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_k( Point3D( {in_0, in_1, in_2} ) );
   }

} // namespace forms
} // namespace hyteg

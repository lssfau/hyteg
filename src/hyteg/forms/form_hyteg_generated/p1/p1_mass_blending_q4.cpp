/*
 * Copyright (c) 2017-2021 Nils Kohl.
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

#include "p1_mass_blending_q4.hpp"

namespace hyteg {
namespace forms {

   void p1_mass_blending_q4::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Blending_DF_Triangle_0_0 = 0;
      real_t Blending_DF_Triangle_0_1 = 0;
      real_t Blending_DF_Triangle_0_2 = 0;
      real_t Blending_DF_Triangle_0_3 = 0;
      real_t Blending_DF_Triangle_1_0 = 0;
      real_t Blending_DF_Triangle_1_1 = 0;
      real_t Blending_DF_Triangle_1_2 = 0;
      real_t Blending_DF_Triangle_1_3 = 0;
      real_t Blending_DF_Triangle_2_0 = 0;
      real_t Blending_DF_Triangle_2_1 = 0;
      real_t Blending_DF_Triangle_2_2 = 0;
      real_t Blending_DF_Triangle_2_3 = 0;
      real_t Blending_DF_Triangle_3_0 = 0;
      real_t Blending_DF_Triangle_3_1 = 0;
      real_t Blending_DF_Triangle_3_2 = 0;
      real_t Blending_DF_Triangle_3_3 = 0;
      real_t Blending_DF_Triangle_4_0 = 0;
      real_t Blending_DF_Triangle_4_1 = 0;
      real_t Blending_DF_Triangle_4_2 = 0;
      real_t Blending_DF_Triangle_4_3 = 0;
      real_t Blending_DF_Triangle_5_0 = 0;
      real_t Blending_DF_Triangle_5_1 = 0;
      real_t Blending_DF_Triangle_5_2 = 0;
      real_t Blending_DF_Triangle_5_3 = 0;
      real_t q_p_0_0 = 0.091576213509770743;
      real_t q_p_0_1 = 0.81684757298045851;
      real_t q_p_1_0 = 0.44594849091596489;
      real_t q_p_1_1 = 0.10810301816807022;
      real_t q_p_2_0 = 0.81684757298045851;
      real_t q_p_2_1 = 0.091576213509770743;
      real_t q_p_3_0 = 0.10810301816807022;
      real_t q_p_3_1 = 0.44594849091596489;
      real_t q_p_4_0 = 0.091576213509770743;
      real_t q_p_4_1 = 0.091576213509770743;
      real_t q_p_5_0 = 0.44594849091596489;
      real_t q_p_5_1 = 0.44594849091596489;
      real_t w_p_0 = 0.054975871827660928;
      real_t w_p_1 = 0.11169079483900572;
      real_t w_p_2 = 0.054975871827660928;
      real_t w_p_3 = 0.11169079483900572;
      real_t w_p_4 = 0.054975871827660928;
      real_t w_p_5 = 0.11169079483900572;
      Blending_DF_Triangle( -p_affine_0_0*q_p_0_0 - p_affine_0_0*q_p_0_1 + p_affine_0_0 + p_affine_1_0*q_p_0_0 + p_affine_2_0*q_p_0_1, -p_affine_0_1*q_p_0_0 - p_affine_0_1*q_p_0_1 + p_affine_0_1 + p_affine_1_1*q_p_0_0 + p_affine_2_1*q_p_0_1, &Blending_DF_Triangle_0_0, &Blending_DF_Triangle_0_1, &Blending_DF_Triangle_0_2, &Blending_DF_Triangle_0_3 );
      Blending_DF_Triangle( -p_affine_0_0*q_p_1_0 - p_affine_0_0*q_p_1_1 + p_affine_0_0 + p_affine_1_0*q_p_1_0 + p_affine_2_0*q_p_1_1, -p_affine_0_1*q_p_1_0 - p_affine_0_1*q_p_1_1 + p_affine_0_1 + p_affine_1_1*q_p_1_0 + p_affine_2_1*q_p_1_1, &Blending_DF_Triangle_1_0, &Blending_DF_Triangle_1_1, &Blending_DF_Triangle_1_2, &Blending_DF_Triangle_1_3 );
      Blending_DF_Triangle( -p_affine_0_0*q_p_2_0 - p_affine_0_0*q_p_2_1 + p_affine_0_0 + p_affine_1_0*q_p_2_0 + p_affine_2_0*q_p_2_1, -p_affine_0_1*q_p_2_0 - p_affine_0_1*q_p_2_1 + p_affine_0_1 + p_affine_1_1*q_p_2_0 + p_affine_2_1*q_p_2_1, &Blending_DF_Triangle_2_0, &Blending_DF_Triangle_2_1, &Blending_DF_Triangle_2_2, &Blending_DF_Triangle_2_3 );
      Blending_DF_Triangle( -p_affine_0_0*q_p_3_0 - p_affine_0_0*q_p_3_1 + p_affine_0_0 + p_affine_1_0*q_p_3_0 + p_affine_2_0*q_p_3_1, -p_affine_0_1*q_p_3_0 - p_affine_0_1*q_p_3_1 + p_affine_0_1 + p_affine_1_1*q_p_3_0 + p_affine_2_1*q_p_3_1, &Blending_DF_Triangle_3_0, &Blending_DF_Triangle_3_1, &Blending_DF_Triangle_3_2, &Blending_DF_Triangle_3_3 );
      Blending_DF_Triangle( -p_affine_0_0*q_p_4_0 - p_affine_0_0*q_p_4_1 + p_affine_0_0 + p_affine_1_0*q_p_4_0 + p_affine_2_0*q_p_4_1, -p_affine_0_1*q_p_4_0 - p_affine_0_1*q_p_4_1 + p_affine_0_1 + p_affine_1_1*q_p_4_0 + p_affine_2_1*q_p_4_1, &Blending_DF_Triangle_4_0, &Blending_DF_Triangle_4_1, &Blending_DF_Triangle_4_2, &Blending_DF_Triangle_4_3 );
      Blending_DF_Triangle( -p_affine_0_0*q_p_5_0 - p_affine_0_0*q_p_5_1 + p_affine_0_0 + p_affine_1_0*q_p_5_0 + p_affine_2_0*q_p_5_1, -p_affine_0_1*q_p_5_0 - p_affine_0_1*q_p_5_1 + p_affine_0_1 + p_affine_1_1*q_p_5_0 + p_affine_2_1*q_p_5_1, &Blending_DF_Triangle_5_0, &Blending_DF_Triangle_5_1, &Blending_DF_Triangle_5_2, &Blending_DF_Triangle_5_3 );
      real_t tmp_0 = -q_p_0_0 - q_p_0_1 + 1;
      real_t tmp_1 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_2 = tmp_1*w_p_0*std::abs(Blending_DF_Triangle_0_0*Blending_DF_Triangle_0_3 - Blending_DF_Triangle_0_1*Blending_DF_Triangle_0_2);
      real_t tmp_3 = -q_p_1_0 - q_p_1_1 + 1;
      real_t tmp_4 = tmp_1*w_p_1*std::abs(Blending_DF_Triangle_1_0*Blending_DF_Triangle_1_3 - Blending_DF_Triangle_1_1*Blending_DF_Triangle_1_2);
      real_t tmp_5 = -q_p_2_0 - q_p_2_1 + 1;
      real_t tmp_6 = tmp_1*w_p_2*std::abs(Blending_DF_Triangle_2_0*Blending_DF_Triangle_2_3 - Blending_DF_Triangle_2_1*Blending_DF_Triangle_2_2);
      real_t tmp_7 = -q_p_3_0 - q_p_3_1 + 1;
      real_t tmp_8 = tmp_1*w_p_3*std::abs(Blending_DF_Triangle_3_0*Blending_DF_Triangle_3_3 - Blending_DF_Triangle_3_1*Blending_DF_Triangle_3_2);
      real_t tmp_9 = -q_p_4_0 - q_p_4_1 + 1;
      real_t tmp_10 = tmp_1*w_p_4*std::abs(Blending_DF_Triangle_4_0*Blending_DF_Triangle_4_3 - Blending_DF_Triangle_4_1*Blending_DF_Triangle_4_2);
      real_t tmp_11 = -q_p_5_0 - q_p_5_1 + 1;
      real_t tmp_12 = tmp_1*w_p_5*std::abs(Blending_DF_Triangle_5_0*Blending_DF_Triangle_5_3 - Blending_DF_Triangle_5_1*Blending_DF_Triangle_5_2);
      real_t tmp_13 = tmp_0*tmp_2;
      real_t tmp_14 = tmp_3*tmp_4;
      real_t tmp_15 = tmp_5*tmp_6;
      real_t tmp_16 = tmp_7*tmp_8;
      real_t tmp_17 = tmp_10*tmp_9;
      real_t tmp_18 = tmp_11*tmp_12;
      real_t tmp_19 = q_p_0_0*tmp_13 + q_p_1_0*tmp_14 + q_p_2_0*tmp_15 + q_p_3_0*tmp_16 + q_p_4_0*tmp_17 + q_p_5_0*tmp_18;
      real_t tmp_20 = q_p_0_1*tmp_13 + q_p_1_1*tmp_14 + q_p_2_1*tmp_15 + q_p_3_1*tmp_16 + q_p_4_1*tmp_17 + q_p_5_1*tmp_18;
      real_t tmp_21 = q_p_0_0*q_p_0_1*tmp_2 + q_p_1_0*q_p_1_1*tmp_4 + q_p_2_0*q_p_2_1*tmp_6 + q_p_3_0*q_p_3_1*tmp_8 + q_p_4_0*q_p_4_1*tmp_10 + q_p_5_0*q_p_5_1*tmp_12;
      real_t a_0_0 = (tmp_0*tmp_0)*tmp_2 + tmp_10*(tmp_9*tmp_9) + (tmp_11*tmp_11)*tmp_12 + (tmp_3*tmp_3)*tmp_4 + (tmp_5*tmp_5)*tmp_6 + (tmp_7*tmp_7)*tmp_8;
      real_t a_0_1 = tmp_19;
      real_t a_0_2 = tmp_20;
      real_t a_1_0 = tmp_19;
      real_t a_1_1 = (q_p_0_0*q_p_0_0)*tmp_2 + (q_p_1_0*q_p_1_0)*tmp_4 + (q_p_2_0*q_p_2_0)*tmp_6 + (q_p_3_0*q_p_3_0)*tmp_8 + (q_p_4_0*q_p_4_0)*tmp_10 + (q_p_5_0*q_p_5_0)*tmp_12;
      real_t a_1_2 = tmp_21;
      real_t a_2_0 = tmp_20;
      real_t a_2_1 = tmp_21;
      real_t a_2_2 = (q_p_0_1*q_p_0_1)*tmp_2 + (q_p_1_1*q_p_1_1)*tmp_4 + (q_p_2_1*q_p_2_1)*tmp_6 + (q_p_3_1*q_p_3_1)*tmp_8 + (q_p_4_1*q_p_4_1)*tmp_10 + (q_p_5_1*q_p_5_1)*tmp_12;
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

   void p1_mass_blending_q4::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_0_0 = 0;
      real_t Blending_DF_Tetrahedron_0_1 = 0;
      real_t Blending_DF_Tetrahedron_0_2 = 0;
      real_t Blending_DF_Tetrahedron_0_3 = 0;
      real_t Blending_DF_Tetrahedron_0_4 = 0;
      real_t Blending_DF_Tetrahedron_0_5 = 0;
      real_t Blending_DF_Tetrahedron_0_6 = 0;
      real_t Blending_DF_Tetrahedron_0_7 = 0;
      real_t Blending_DF_Tetrahedron_0_8 = 0;
      real_t Blending_DF_Tetrahedron_1_0 = 0;
      real_t Blending_DF_Tetrahedron_1_1 = 0;
      real_t Blending_DF_Tetrahedron_1_2 = 0;
      real_t Blending_DF_Tetrahedron_1_3 = 0;
      real_t Blending_DF_Tetrahedron_1_4 = 0;
      real_t Blending_DF_Tetrahedron_1_5 = 0;
      real_t Blending_DF_Tetrahedron_1_6 = 0;
      real_t Blending_DF_Tetrahedron_1_7 = 0;
      real_t Blending_DF_Tetrahedron_1_8 = 0;
      real_t Blending_DF_Tetrahedron_2_0 = 0;
      real_t Blending_DF_Tetrahedron_2_1 = 0;
      real_t Blending_DF_Tetrahedron_2_2 = 0;
      real_t Blending_DF_Tetrahedron_2_3 = 0;
      real_t Blending_DF_Tetrahedron_2_4 = 0;
      real_t Blending_DF_Tetrahedron_2_5 = 0;
      real_t Blending_DF_Tetrahedron_2_6 = 0;
      real_t Blending_DF_Tetrahedron_2_7 = 0;
      real_t Blending_DF_Tetrahedron_2_8 = 0;
      real_t Blending_DF_Tetrahedron_3_0 = 0;
      real_t Blending_DF_Tetrahedron_3_1 = 0;
      real_t Blending_DF_Tetrahedron_3_2 = 0;
      real_t Blending_DF_Tetrahedron_3_3 = 0;
      real_t Blending_DF_Tetrahedron_3_4 = 0;
      real_t Blending_DF_Tetrahedron_3_5 = 0;
      real_t Blending_DF_Tetrahedron_3_6 = 0;
      real_t Blending_DF_Tetrahedron_3_7 = 0;
      real_t Blending_DF_Tetrahedron_3_8 = 0;
      real_t Blending_DF_Tetrahedron_4_0 = 0;
      real_t Blending_DF_Tetrahedron_4_1 = 0;
      real_t Blending_DF_Tetrahedron_4_2 = 0;
      real_t Blending_DF_Tetrahedron_4_3 = 0;
      real_t Blending_DF_Tetrahedron_4_4 = 0;
      real_t Blending_DF_Tetrahedron_4_5 = 0;
      real_t Blending_DF_Tetrahedron_4_6 = 0;
      real_t Blending_DF_Tetrahedron_4_7 = 0;
      real_t Blending_DF_Tetrahedron_4_8 = 0;
      real_t Blending_DF_Tetrahedron_5_0 = 0;
      real_t Blending_DF_Tetrahedron_5_1 = 0;
      real_t Blending_DF_Tetrahedron_5_2 = 0;
      real_t Blending_DF_Tetrahedron_5_3 = 0;
      real_t Blending_DF_Tetrahedron_5_4 = 0;
      real_t Blending_DF_Tetrahedron_5_5 = 0;
      real_t Blending_DF_Tetrahedron_5_6 = 0;
      real_t Blending_DF_Tetrahedron_5_7 = 0;
      real_t Blending_DF_Tetrahedron_5_8 = 0;
      real_t Blending_DF_Tetrahedron_6_0 = 0;
      real_t Blending_DF_Tetrahedron_6_1 = 0;
      real_t Blending_DF_Tetrahedron_6_2 = 0;
      real_t Blending_DF_Tetrahedron_6_3 = 0;
      real_t Blending_DF_Tetrahedron_6_4 = 0;
      real_t Blending_DF_Tetrahedron_6_5 = 0;
      real_t Blending_DF_Tetrahedron_6_6 = 0;
      real_t Blending_DF_Tetrahedron_6_7 = 0;
      real_t Blending_DF_Tetrahedron_6_8 = 0;
      real_t Blending_DF_Tetrahedron_7_0 = 0;
      real_t Blending_DF_Tetrahedron_7_1 = 0;
      real_t Blending_DF_Tetrahedron_7_2 = 0;
      real_t Blending_DF_Tetrahedron_7_3 = 0;
      real_t Blending_DF_Tetrahedron_7_4 = 0;
      real_t Blending_DF_Tetrahedron_7_5 = 0;
      real_t Blending_DF_Tetrahedron_7_6 = 0;
      real_t Blending_DF_Tetrahedron_7_7 = 0;
      real_t Blending_DF_Tetrahedron_7_8 = 0;
      real_t Blending_DF_Tetrahedron_8_0 = 0;
      real_t Blending_DF_Tetrahedron_8_1 = 0;
      real_t Blending_DF_Tetrahedron_8_2 = 0;
      real_t Blending_DF_Tetrahedron_8_3 = 0;
      real_t Blending_DF_Tetrahedron_8_4 = 0;
      real_t Blending_DF_Tetrahedron_8_5 = 0;
      real_t Blending_DF_Tetrahedron_8_6 = 0;
      real_t Blending_DF_Tetrahedron_8_7 = 0;
      real_t Blending_DF_Tetrahedron_8_8 = 0;
      real_t Blending_DF_Tetrahedron_9_0 = 0;
      real_t Blending_DF_Tetrahedron_9_1 = 0;
      real_t Blending_DF_Tetrahedron_9_2 = 0;
      real_t Blending_DF_Tetrahedron_9_3 = 0;
      real_t Blending_DF_Tetrahedron_9_4 = 0;
      real_t Blending_DF_Tetrahedron_9_5 = 0;
      real_t Blending_DF_Tetrahedron_9_6 = 0;
      real_t Blending_DF_Tetrahedron_9_7 = 0;
      real_t Blending_DF_Tetrahedron_9_8 = 0;
      real_t Blending_DF_Tetrahedron_10_0 = 0;
      real_t Blending_DF_Tetrahedron_10_1 = 0;
      real_t Blending_DF_Tetrahedron_10_2 = 0;
      real_t Blending_DF_Tetrahedron_10_3 = 0;
      real_t Blending_DF_Tetrahedron_10_4 = 0;
      real_t Blending_DF_Tetrahedron_10_5 = 0;
      real_t Blending_DF_Tetrahedron_10_6 = 0;
      real_t Blending_DF_Tetrahedron_10_7 = 0;
      real_t Blending_DF_Tetrahedron_10_8 = 0;
      real_t q_p_0_0 = 0.040490506727590428;
      real_t q_p_0_1 = 0.01356070187980288;
      real_t q_p_0_2 = 0.77125473269537614;
      real_t q_p_1_0 = 0.75250850700965499;
      real_t q_p_1_1 = 0.068099370938206658;
      real_t q_p_1_2 = 0.097987203649279112;
      real_t q_p_2_0 = 0.067223294893383398;
      real_t q_p_2_1 = 0.035183929773598722;
      real_t q_p_2_2 = 0.15636389323939531;
      real_t q_p_3_0 = 0.41926631387951302;
      real_t q_p_3_1 = 0.047781435559086663;
      real_t q_p_3_2 = 0.47961101102565512;
      real_t q_p_4_0 = 0.45076587609127677;
      real_t q_p_4_1 = 0.059456616299433829;
      real_t q_p_4_2 = 0.056824017127933668;
      real_t q_p_5_0 = 0.12941137378891041;
      real_t q_p_5_1 = 0.33019041483746447;
      real_t q_p_5_2 = 0.0023910074574393651;
      real_t q_p_6_0 = 0.1215419913339278;
      real_t q_p_6_1 = 0.30649398842969028;
      real_t q_p_6_2 = 0.56297276014304609;
      real_t q_p_7_0 = 0.097204644587583267;
      real_t q_p_7_1 = 0.68439041545304002;
      real_t q_p_7_2 = 0.11180076739738309;
      real_t q_p_8_0 = 0.029569495206479609;
      real_t q_p_8_1 = 0.31790356021339461;
      real_t q_p_8_2 = 0.32329398483747901;
      real_t q_p_9_0 = 0.43271023904776862;
      real_t q_p_9_1 = 0.35382323920929709;
      real_t q_p_9_2 = 0.1096224053319412;
      real_t q_p_10_0 = 0.24027666492807259;
      real_t q_p_10_1 = 0.12680172591539199;
      real_t q_p_10_2 = 0.32847320672203839;
      real_t w_p_0 = 0.006541848487473325;
      real_t w_p_1 = 0.0092122281926561474;
      real_t w_p_2 = 0.0092322998119293929;
      real_t w_p_3 = 0.0099888641910932524;
      real_t w_p_4 = 0.011578327656272558;
      real_t w_p_5 = 0.012693785874259723;
      real_t w_p_6 = 0.013237780011337548;
      real_t w_p_7 = 0.017744672359248346;
      real_t w_p_8 = 0.018372372071416277;
      real_t w_p_9 = 0.025829352669374347;
      real_t w_p_10 = 0.032235135341605742;
      Blending_DF_Tetrahedron( -p_affine_0_0*q_p_0_0 - p_affine_0_0*q_p_0_1 - p_affine_0_0*q_p_0_2 + p_affine_0_0 + p_affine_1_0*q_p_0_0 + p_affine_2_0*q_p_0_1 + p_affine_3_0*q_p_0_2, -p_affine_0_1*q_p_0_0 - p_affine_0_1*q_p_0_1 - p_affine_0_1*q_p_0_2 + p_affine_0_1 + p_affine_1_1*q_p_0_0 + p_affine_2_1*q_p_0_1 + p_affine_3_1*q_p_0_2, -p_affine_0_2*q_p_0_0 - p_affine_0_2*q_p_0_1 - p_affine_0_2*q_p_0_2 + p_affine_0_2 + p_affine_1_2*q_p_0_0 + p_affine_2_2*q_p_0_1 + p_affine_3_2*q_p_0_2, &Blending_DF_Tetrahedron_0_0, &Blending_DF_Tetrahedron_0_1, &Blending_DF_Tetrahedron_0_2, &Blending_DF_Tetrahedron_0_3, &Blending_DF_Tetrahedron_0_4, &Blending_DF_Tetrahedron_0_5, &Blending_DF_Tetrahedron_0_6, &Blending_DF_Tetrahedron_0_7, &Blending_DF_Tetrahedron_0_8 );
      Blending_DF_Tetrahedron( -p_affine_0_0*q_p_1_0 - p_affine_0_0*q_p_1_1 - p_affine_0_0*q_p_1_2 + p_affine_0_0 + p_affine_1_0*q_p_1_0 + p_affine_2_0*q_p_1_1 + p_affine_3_0*q_p_1_2, -p_affine_0_1*q_p_1_0 - p_affine_0_1*q_p_1_1 - p_affine_0_1*q_p_1_2 + p_affine_0_1 + p_affine_1_1*q_p_1_0 + p_affine_2_1*q_p_1_1 + p_affine_3_1*q_p_1_2, -p_affine_0_2*q_p_1_0 - p_affine_0_2*q_p_1_1 - p_affine_0_2*q_p_1_2 + p_affine_0_2 + p_affine_1_2*q_p_1_0 + p_affine_2_2*q_p_1_1 + p_affine_3_2*q_p_1_2, &Blending_DF_Tetrahedron_1_0, &Blending_DF_Tetrahedron_1_1, &Blending_DF_Tetrahedron_1_2, &Blending_DF_Tetrahedron_1_3, &Blending_DF_Tetrahedron_1_4, &Blending_DF_Tetrahedron_1_5, &Blending_DF_Tetrahedron_1_6, &Blending_DF_Tetrahedron_1_7, &Blending_DF_Tetrahedron_1_8 );
      Blending_DF_Tetrahedron( -p_affine_0_0*q_p_10_0 - p_affine_0_0*q_p_10_1 - p_affine_0_0*q_p_10_2 + p_affine_0_0 + p_affine_1_0*q_p_10_0 + p_affine_2_0*q_p_10_1 + p_affine_3_0*q_p_10_2, -p_affine_0_1*q_p_10_0 - p_affine_0_1*q_p_10_1 - p_affine_0_1*q_p_10_2 + p_affine_0_1 + p_affine_1_1*q_p_10_0 + p_affine_2_1*q_p_10_1 + p_affine_3_1*q_p_10_2, -p_affine_0_2*q_p_10_0 - p_affine_0_2*q_p_10_1 - p_affine_0_2*q_p_10_2 + p_affine_0_2 + p_affine_1_2*q_p_10_0 + p_affine_2_2*q_p_10_1 + p_affine_3_2*q_p_10_2, &Blending_DF_Tetrahedron_2_0, &Blending_DF_Tetrahedron_2_1, &Blending_DF_Tetrahedron_2_2, &Blending_DF_Tetrahedron_2_3, &Blending_DF_Tetrahedron_2_4, &Blending_DF_Tetrahedron_2_5, &Blending_DF_Tetrahedron_2_6, &Blending_DF_Tetrahedron_2_7, &Blending_DF_Tetrahedron_2_8 );
      Blending_DF_Tetrahedron( -p_affine_0_0*q_p_2_0 - p_affine_0_0*q_p_2_1 - p_affine_0_0*q_p_2_2 + p_affine_0_0 + p_affine_1_0*q_p_2_0 + p_affine_2_0*q_p_2_1 + p_affine_3_0*q_p_2_2, -p_affine_0_1*q_p_2_0 - p_affine_0_1*q_p_2_1 - p_affine_0_1*q_p_2_2 + p_affine_0_1 + p_affine_1_1*q_p_2_0 + p_affine_2_1*q_p_2_1 + p_affine_3_1*q_p_2_2, -p_affine_0_2*q_p_2_0 - p_affine_0_2*q_p_2_1 - p_affine_0_2*q_p_2_2 + p_affine_0_2 + p_affine_1_2*q_p_2_0 + p_affine_2_2*q_p_2_1 + p_affine_3_2*q_p_2_2, &Blending_DF_Tetrahedron_3_0, &Blending_DF_Tetrahedron_3_1, &Blending_DF_Tetrahedron_3_2, &Blending_DF_Tetrahedron_3_3, &Blending_DF_Tetrahedron_3_4, &Blending_DF_Tetrahedron_3_5, &Blending_DF_Tetrahedron_3_6, &Blending_DF_Tetrahedron_3_7, &Blending_DF_Tetrahedron_3_8 );
      Blending_DF_Tetrahedron( -p_affine_0_0*q_p_3_0 - p_affine_0_0*q_p_3_1 - p_affine_0_0*q_p_3_2 + p_affine_0_0 + p_affine_1_0*q_p_3_0 + p_affine_2_0*q_p_3_1 + p_affine_3_0*q_p_3_2, -p_affine_0_1*q_p_3_0 - p_affine_0_1*q_p_3_1 - p_affine_0_1*q_p_3_2 + p_affine_0_1 + p_affine_1_1*q_p_3_0 + p_affine_2_1*q_p_3_1 + p_affine_3_1*q_p_3_2, -p_affine_0_2*q_p_3_0 - p_affine_0_2*q_p_3_1 - p_affine_0_2*q_p_3_2 + p_affine_0_2 + p_affine_1_2*q_p_3_0 + p_affine_2_2*q_p_3_1 + p_affine_3_2*q_p_3_2, &Blending_DF_Tetrahedron_4_0, &Blending_DF_Tetrahedron_4_1, &Blending_DF_Tetrahedron_4_2, &Blending_DF_Tetrahedron_4_3, &Blending_DF_Tetrahedron_4_4, &Blending_DF_Tetrahedron_4_5, &Blending_DF_Tetrahedron_4_6, &Blending_DF_Tetrahedron_4_7, &Blending_DF_Tetrahedron_4_8 );
      Blending_DF_Tetrahedron( -p_affine_0_0*q_p_4_0 - p_affine_0_0*q_p_4_1 - p_affine_0_0*q_p_4_2 + p_affine_0_0 + p_affine_1_0*q_p_4_0 + p_affine_2_0*q_p_4_1 + p_affine_3_0*q_p_4_2, -p_affine_0_1*q_p_4_0 - p_affine_0_1*q_p_4_1 - p_affine_0_1*q_p_4_2 + p_affine_0_1 + p_affine_1_1*q_p_4_0 + p_affine_2_1*q_p_4_1 + p_affine_3_1*q_p_4_2, -p_affine_0_2*q_p_4_0 - p_affine_0_2*q_p_4_1 - p_affine_0_2*q_p_4_2 + p_affine_0_2 + p_affine_1_2*q_p_4_0 + p_affine_2_2*q_p_4_1 + p_affine_3_2*q_p_4_2, &Blending_DF_Tetrahedron_5_0, &Blending_DF_Tetrahedron_5_1, &Blending_DF_Tetrahedron_5_2, &Blending_DF_Tetrahedron_5_3, &Blending_DF_Tetrahedron_5_4, &Blending_DF_Tetrahedron_5_5, &Blending_DF_Tetrahedron_5_6, &Blending_DF_Tetrahedron_5_7, &Blending_DF_Tetrahedron_5_8 );
      Blending_DF_Tetrahedron( -p_affine_0_0*q_p_5_0 - p_affine_0_0*q_p_5_1 - p_affine_0_0*q_p_5_2 + p_affine_0_0 + p_affine_1_0*q_p_5_0 + p_affine_2_0*q_p_5_1 + p_affine_3_0*q_p_5_2, -p_affine_0_1*q_p_5_0 - p_affine_0_1*q_p_5_1 - p_affine_0_1*q_p_5_2 + p_affine_0_1 + p_affine_1_1*q_p_5_0 + p_affine_2_1*q_p_5_1 + p_affine_3_1*q_p_5_2, -p_affine_0_2*q_p_5_0 - p_affine_0_2*q_p_5_1 - p_affine_0_2*q_p_5_2 + p_affine_0_2 + p_affine_1_2*q_p_5_0 + p_affine_2_2*q_p_5_1 + p_affine_3_2*q_p_5_2, &Blending_DF_Tetrahedron_6_0, &Blending_DF_Tetrahedron_6_1, &Blending_DF_Tetrahedron_6_2, &Blending_DF_Tetrahedron_6_3, &Blending_DF_Tetrahedron_6_4, &Blending_DF_Tetrahedron_6_5, &Blending_DF_Tetrahedron_6_6, &Blending_DF_Tetrahedron_6_7, &Blending_DF_Tetrahedron_6_8 );
      Blending_DF_Tetrahedron( -p_affine_0_0*q_p_6_0 - p_affine_0_0*q_p_6_1 - p_affine_0_0*q_p_6_2 + p_affine_0_0 + p_affine_1_0*q_p_6_0 + p_affine_2_0*q_p_6_1 + p_affine_3_0*q_p_6_2, -p_affine_0_1*q_p_6_0 - p_affine_0_1*q_p_6_1 - p_affine_0_1*q_p_6_2 + p_affine_0_1 + p_affine_1_1*q_p_6_0 + p_affine_2_1*q_p_6_1 + p_affine_3_1*q_p_6_2, -p_affine_0_2*q_p_6_0 - p_affine_0_2*q_p_6_1 - p_affine_0_2*q_p_6_2 + p_affine_0_2 + p_affine_1_2*q_p_6_0 + p_affine_2_2*q_p_6_1 + p_affine_3_2*q_p_6_2, &Blending_DF_Tetrahedron_7_0, &Blending_DF_Tetrahedron_7_1, &Blending_DF_Tetrahedron_7_2, &Blending_DF_Tetrahedron_7_3, &Blending_DF_Tetrahedron_7_4, &Blending_DF_Tetrahedron_7_5, &Blending_DF_Tetrahedron_7_6, &Blending_DF_Tetrahedron_7_7, &Blending_DF_Tetrahedron_7_8 );
      Blending_DF_Tetrahedron( -p_affine_0_0*q_p_7_0 - p_affine_0_0*q_p_7_1 - p_affine_0_0*q_p_7_2 + p_affine_0_0 + p_affine_1_0*q_p_7_0 + p_affine_2_0*q_p_7_1 + p_affine_3_0*q_p_7_2, -p_affine_0_1*q_p_7_0 - p_affine_0_1*q_p_7_1 - p_affine_0_1*q_p_7_2 + p_affine_0_1 + p_affine_1_1*q_p_7_0 + p_affine_2_1*q_p_7_1 + p_affine_3_1*q_p_7_2, -p_affine_0_2*q_p_7_0 - p_affine_0_2*q_p_7_1 - p_affine_0_2*q_p_7_2 + p_affine_0_2 + p_affine_1_2*q_p_7_0 + p_affine_2_2*q_p_7_1 + p_affine_3_2*q_p_7_2, &Blending_DF_Tetrahedron_8_0, &Blending_DF_Tetrahedron_8_1, &Blending_DF_Tetrahedron_8_2, &Blending_DF_Tetrahedron_8_3, &Blending_DF_Tetrahedron_8_4, &Blending_DF_Tetrahedron_8_5, &Blending_DF_Tetrahedron_8_6, &Blending_DF_Tetrahedron_8_7, &Blending_DF_Tetrahedron_8_8 );
      Blending_DF_Tetrahedron( -p_affine_0_0*q_p_8_0 - p_affine_0_0*q_p_8_1 - p_affine_0_0*q_p_8_2 + p_affine_0_0 + p_affine_1_0*q_p_8_0 + p_affine_2_0*q_p_8_1 + p_affine_3_0*q_p_8_2, -p_affine_0_1*q_p_8_0 - p_affine_0_1*q_p_8_1 - p_affine_0_1*q_p_8_2 + p_affine_0_1 + p_affine_1_1*q_p_8_0 + p_affine_2_1*q_p_8_1 + p_affine_3_1*q_p_8_2, -p_affine_0_2*q_p_8_0 - p_affine_0_2*q_p_8_1 - p_affine_0_2*q_p_8_2 + p_affine_0_2 + p_affine_1_2*q_p_8_0 + p_affine_2_2*q_p_8_1 + p_affine_3_2*q_p_8_2, &Blending_DF_Tetrahedron_9_0, &Blending_DF_Tetrahedron_9_1, &Blending_DF_Tetrahedron_9_2, &Blending_DF_Tetrahedron_9_3, &Blending_DF_Tetrahedron_9_4, &Blending_DF_Tetrahedron_9_5, &Blending_DF_Tetrahedron_9_6, &Blending_DF_Tetrahedron_9_7, &Blending_DF_Tetrahedron_9_8 );
      Blending_DF_Tetrahedron( -p_affine_0_0*q_p_9_0 - p_affine_0_0*q_p_9_1 - p_affine_0_0*q_p_9_2 + p_affine_0_0 + p_affine_1_0*q_p_9_0 + p_affine_2_0*q_p_9_1 + p_affine_3_0*q_p_9_2, -p_affine_0_1*q_p_9_0 - p_affine_0_1*q_p_9_1 - p_affine_0_1*q_p_9_2 + p_affine_0_1 + p_affine_1_1*q_p_9_0 + p_affine_2_1*q_p_9_1 + p_affine_3_1*q_p_9_2, -p_affine_0_2*q_p_9_0 - p_affine_0_2*q_p_9_1 - p_affine_0_2*q_p_9_2 + p_affine_0_2 + p_affine_1_2*q_p_9_0 + p_affine_2_2*q_p_9_1 + p_affine_3_2*q_p_9_2, &Blending_DF_Tetrahedron_10_0, &Blending_DF_Tetrahedron_10_1, &Blending_DF_Tetrahedron_10_2, &Blending_DF_Tetrahedron_10_3, &Blending_DF_Tetrahedron_10_4, &Blending_DF_Tetrahedron_10_5, &Blending_DF_Tetrahedron_10_6, &Blending_DF_Tetrahedron_10_7, &Blending_DF_Tetrahedron_10_8 );
      real_t tmp_0 = -q_p_0_0 - q_p_0_1 - q_p_0_2 + 1;
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
      real_t tmp_13 = std::abs(p_affine_0_0*tmp_10 - p_affine_0_0*tmp_3 + p_affine_0_1*tmp_11 - p_affine_0_1*tmp_6 + p_affine_0_2*tmp_12 - p_affine_0_2*tmp_9 - p_affine_1_0*tmp_10 + p_affine_1_0*tmp_3 - p_affine_1_1*tmp_11 + p_affine_1_1*tmp_6 - p_affine_1_2*tmp_12 + p_affine_1_2*tmp_9 - p_affine_2_0*tmp_5 + p_affine_2_0*tmp_8 + p_affine_2_1*tmp_2 - p_affine_2_1*tmp_7 - p_affine_2_2*tmp_1 + p_affine_2_2*tmp_4 + p_affine_3_0*tmp_5 - p_affine_3_0*tmp_8 - p_affine_3_1*tmp_2 + p_affine_3_1*tmp_7 + p_affine_3_2*tmp_1 - p_affine_3_2*tmp_4);
      real_t tmp_14 = tmp_13*w_p_0*std::abs(Blending_DF_Tetrahedron_0_0*Blending_DF_Tetrahedron_0_4*Blending_DF_Tetrahedron_0_8 - Blending_DF_Tetrahedron_0_0*Blending_DF_Tetrahedron_0_5*Blending_DF_Tetrahedron_0_7 - Blending_DF_Tetrahedron_0_1*Blending_DF_Tetrahedron_0_3*Blending_DF_Tetrahedron_0_8 + Blending_DF_Tetrahedron_0_1*Blending_DF_Tetrahedron_0_5*Blending_DF_Tetrahedron_0_6 + Blending_DF_Tetrahedron_0_2*Blending_DF_Tetrahedron_0_3*Blending_DF_Tetrahedron_0_7 - Blending_DF_Tetrahedron_0_2*Blending_DF_Tetrahedron_0_4*Blending_DF_Tetrahedron_0_6);
      real_t tmp_15 = -q_p_1_0 - q_p_1_1 - q_p_1_2 + 1;
      real_t tmp_16 = tmp_13*w_p_1*std::abs(Blending_DF_Tetrahedron_1_0*Blending_DF_Tetrahedron_1_4*Blending_DF_Tetrahedron_1_8 - Blending_DF_Tetrahedron_1_0*Blending_DF_Tetrahedron_1_5*Blending_DF_Tetrahedron_1_7 - Blending_DF_Tetrahedron_1_1*Blending_DF_Tetrahedron_1_3*Blending_DF_Tetrahedron_1_8 + Blending_DF_Tetrahedron_1_1*Blending_DF_Tetrahedron_1_5*Blending_DF_Tetrahedron_1_6 + Blending_DF_Tetrahedron_1_2*Blending_DF_Tetrahedron_1_3*Blending_DF_Tetrahedron_1_7 - Blending_DF_Tetrahedron_1_2*Blending_DF_Tetrahedron_1_4*Blending_DF_Tetrahedron_1_6);
      real_t tmp_17 = -q_p_10_0 - q_p_10_1 - q_p_10_2 + 1;
      real_t tmp_18 = tmp_13*w_p_10*std::abs(Blending_DF_Tetrahedron_2_0*Blending_DF_Tetrahedron_2_4*Blending_DF_Tetrahedron_2_8 - Blending_DF_Tetrahedron_2_0*Blending_DF_Tetrahedron_2_5*Blending_DF_Tetrahedron_2_7 - Blending_DF_Tetrahedron_2_1*Blending_DF_Tetrahedron_2_3*Blending_DF_Tetrahedron_2_8 + Blending_DF_Tetrahedron_2_1*Blending_DF_Tetrahedron_2_5*Blending_DF_Tetrahedron_2_6 + Blending_DF_Tetrahedron_2_2*Blending_DF_Tetrahedron_2_3*Blending_DF_Tetrahedron_2_7 - Blending_DF_Tetrahedron_2_2*Blending_DF_Tetrahedron_2_4*Blending_DF_Tetrahedron_2_6);
      real_t tmp_19 = -q_p_2_0 - q_p_2_1 - q_p_2_2 + 1;
      real_t tmp_20 = tmp_13*w_p_2*std::abs(Blending_DF_Tetrahedron_3_0*Blending_DF_Tetrahedron_3_4*Blending_DF_Tetrahedron_3_8 - Blending_DF_Tetrahedron_3_0*Blending_DF_Tetrahedron_3_5*Blending_DF_Tetrahedron_3_7 - Blending_DF_Tetrahedron_3_1*Blending_DF_Tetrahedron_3_3*Blending_DF_Tetrahedron_3_8 + Blending_DF_Tetrahedron_3_1*Blending_DF_Tetrahedron_3_5*Blending_DF_Tetrahedron_3_6 + Blending_DF_Tetrahedron_3_2*Blending_DF_Tetrahedron_3_3*Blending_DF_Tetrahedron_3_7 - Blending_DF_Tetrahedron_3_2*Blending_DF_Tetrahedron_3_4*Blending_DF_Tetrahedron_3_6);
      real_t tmp_21 = -q_p_3_0 - q_p_3_1 - q_p_3_2 + 1;
      real_t tmp_22 = tmp_13*w_p_3*std::abs(Blending_DF_Tetrahedron_4_0*Blending_DF_Tetrahedron_4_4*Blending_DF_Tetrahedron_4_8 - Blending_DF_Tetrahedron_4_0*Blending_DF_Tetrahedron_4_5*Blending_DF_Tetrahedron_4_7 - Blending_DF_Tetrahedron_4_1*Blending_DF_Tetrahedron_4_3*Blending_DF_Tetrahedron_4_8 + Blending_DF_Tetrahedron_4_1*Blending_DF_Tetrahedron_4_5*Blending_DF_Tetrahedron_4_6 + Blending_DF_Tetrahedron_4_2*Blending_DF_Tetrahedron_4_3*Blending_DF_Tetrahedron_4_7 - Blending_DF_Tetrahedron_4_2*Blending_DF_Tetrahedron_4_4*Blending_DF_Tetrahedron_4_6);
      real_t tmp_23 = -q_p_4_0 - q_p_4_1 - q_p_4_2 + 1;
      real_t tmp_24 = tmp_13*w_p_4*std::abs(Blending_DF_Tetrahedron_5_0*Blending_DF_Tetrahedron_5_4*Blending_DF_Tetrahedron_5_8 - Blending_DF_Tetrahedron_5_0*Blending_DF_Tetrahedron_5_5*Blending_DF_Tetrahedron_5_7 - Blending_DF_Tetrahedron_5_1*Blending_DF_Tetrahedron_5_3*Blending_DF_Tetrahedron_5_8 + Blending_DF_Tetrahedron_5_1*Blending_DF_Tetrahedron_5_5*Blending_DF_Tetrahedron_5_6 + Blending_DF_Tetrahedron_5_2*Blending_DF_Tetrahedron_5_3*Blending_DF_Tetrahedron_5_7 - Blending_DF_Tetrahedron_5_2*Blending_DF_Tetrahedron_5_4*Blending_DF_Tetrahedron_5_6);
      real_t tmp_25 = -q_p_5_0 - q_p_5_1 - q_p_5_2 + 1;
      real_t tmp_26 = tmp_13*w_p_5*std::abs(Blending_DF_Tetrahedron_6_0*Blending_DF_Tetrahedron_6_4*Blending_DF_Tetrahedron_6_8 - Blending_DF_Tetrahedron_6_0*Blending_DF_Tetrahedron_6_5*Blending_DF_Tetrahedron_6_7 - Blending_DF_Tetrahedron_6_1*Blending_DF_Tetrahedron_6_3*Blending_DF_Tetrahedron_6_8 + Blending_DF_Tetrahedron_6_1*Blending_DF_Tetrahedron_6_5*Blending_DF_Tetrahedron_6_6 + Blending_DF_Tetrahedron_6_2*Blending_DF_Tetrahedron_6_3*Blending_DF_Tetrahedron_6_7 - Blending_DF_Tetrahedron_6_2*Blending_DF_Tetrahedron_6_4*Blending_DF_Tetrahedron_6_6);
      real_t tmp_27 = -q_p_6_0 - q_p_6_1 - q_p_6_2 + 1;
      real_t tmp_28 = tmp_13*w_p_6*std::abs(Blending_DF_Tetrahedron_7_0*Blending_DF_Tetrahedron_7_4*Blending_DF_Tetrahedron_7_8 - Blending_DF_Tetrahedron_7_0*Blending_DF_Tetrahedron_7_5*Blending_DF_Tetrahedron_7_7 - Blending_DF_Tetrahedron_7_1*Blending_DF_Tetrahedron_7_3*Blending_DF_Tetrahedron_7_8 + Blending_DF_Tetrahedron_7_1*Blending_DF_Tetrahedron_7_5*Blending_DF_Tetrahedron_7_6 + Blending_DF_Tetrahedron_7_2*Blending_DF_Tetrahedron_7_3*Blending_DF_Tetrahedron_7_7 - Blending_DF_Tetrahedron_7_2*Blending_DF_Tetrahedron_7_4*Blending_DF_Tetrahedron_7_6);
      real_t tmp_29 = -q_p_7_0 - q_p_7_1 - q_p_7_2 + 1;
      real_t tmp_30 = tmp_13*w_p_7*std::abs(Blending_DF_Tetrahedron_8_0*Blending_DF_Tetrahedron_8_4*Blending_DF_Tetrahedron_8_8 - Blending_DF_Tetrahedron_8_0*Blending_DF_Tetrahedron_8_5*Blending_DF_Tetrahedron_8_7 - Blending_DF_Tetrahedron_8_1*Blending_DF_Tetrahedron_8_3*Blending_DF_Tetrahedron_8_8 + Blending_DF_Tetrahedron_8_1*Blending_DF_Tetrahedron_8_5*Blending_DF_Tetrahedron_8_6 + Blending_DF_Tetrahedron_8_2*Blending_DF_Tetrahedron_8_3*Blending_DF_Tetrahedron_8_7 - Blending_DF_Tetrahedron_8_2*Blending_DF_Tetrahedron_8_4*Blending_DF_Tetrahedron_8_6);
      real_t tmp_31 = -q_p_8_0 - q_p_8_1 - q_p_8_2 + 1;
      real_t tmp_32 = tmp_13*w_p_8*std::abs(Blending_DF_Tetrahedron_9_0*Blending_DF_Tetrahedron_9_4*Blending_DF_Tetrahedron_9_8 - Blending_DF_Tetrahedron_9_0*Blending_DF_Tetrahedron_9_5*Blending_DF_Tetrahedron_9_7 - Blending_DF_Tetrahedron_9_1*Blending_DF_Tetrahedron_9_3*Blending_DF_Tetrahedron_9_8 + Blending_DF_Tetrahedron_9_1*Blending_DF_Tetrahedron_9_5*Blending_DF_Tetrahedron_9_6 + Blending_DF_Tetrahedron_9_2*Blending_DF_Tetrahedron_9_3*Blending_DF_Tetrahedron_9_7 - Blending_DF_Tetrahedron_9_2*Blending_DF_Tetrahedron_9_4*Blending_DF_Tetrahedron_9_6);
      real_t tmp_33 = -q_p_9_0 - q_p_9_1 - q_p_9_2 + 1;
      real_t tmp_34 = tmp_13*w_p_9*std::abs(Blending_DF_Tetrahedron_10_0*Blending_DF_Tetrahedron_10_4*Blending_DF_Tetrahedron_10_8 - Blending_DF_Tetrahedron_10_0*Blending_DF_Tetrahedron_10_5*Blending_DF_Tetrahedron_10_7 - Blending_DF_Tetrahedron_10_1*Blending_DF_Tetrahedron_10_3*Blending_DF_Tetrahedron_10_8 + Blending_DF_Tetrahedron_10_1*Blending_DF_Tetrahedron_10_5*Blending_DF_Tetrahedron_10_6 + Blending_DF_Tetrahedron_10_2*Blending_DF_Tetrahedron_10_3*Blending_DF_Tetrahedron_10_7 - Blending_DF_Tetrahedron_10_2*Blending_DF_Tetrahedron_10_4*Blending_DF_Tetrahedron_10_6);
      real_t tmp_35 = tmp_0*tmp_14;
      real_t tmp_36 = tmp_17*tmp_18;
      real_t tmp_37 = tmp_15*tmp_16;
      real_t tmp_38 = tmp_19*tmp_20;
      real_t tmp_39 = tmp_21*tmp_22;
      real_t tmp_40 = tmp_23*tmp_24;
      real_t tmp_41 = tmp_25*tmp_26;
      real_t tmp_42 = tmp_27*tmp_28;
      real_t tmp_43 = tmp_29*tmp_30;
      real_t tmp_44 = tmp_31*tmp_32;
      real_t tmp_45 = tmp_33*tmp_34;
      real_t tmp_46 = q_p_0_0*tmp_35 + q_p_10_0*tmp_36 + q_p_1_0*tmp_37 + q_p_2_0*tmp_38 + q_p_3_0*tmp_39 + q_p_4_0*tmp_40 + q_p_5_0*tmp_41 + q_p_6_0*tmp_42 + q_p_7_0*tmp_43 + q_p_8_0*tmp_44 + q_p_9_0*tmp_45;
      real_t tmp_47 = q_p_0_1*tmp_35 + q_p_10_1*tmp_36 + q_p_1_1*tmp_37 + q_p_2_1*tmp_38 + q_p_3_1*tmp_39 + q_p_4_1*tmp_40 + q_p_5_1*tmp_41 + q_p_6_1*tmp_42 + q_p_7_1*tmp_43 + q_p_8_1*tmp_44 + q_p_9_1*tmp_45;
      real_t tmp_48 = q_p_0_2*tmp_35 + q_p_10_2*tmp_36 + q_p_1_2*tmp_37 + q_p_2_2*tmp_38 + q_p_3_2*tmp_39 + q_p_4_2*tmp_40 + q_p_5_2*tmp_41 + q_p_6_2*tmp_42 + q_p_7_2*tmp_43 + q_p_8_2*tmp_44 + q_p_9_2*tmp_45;
      real_t tmp_49 = q_p_0_0*tmp_14;
      real_t tmp_50 = q_p_10_0*tmp_18;
      real_t tmp_51 = q_p_1_0*tmp_16;
      real_t tmp_52 = q_p_2_0*tmp_20;
      real_t tmp_53 = q_p_3_0*tmp_22;
      real_t tmp_54 = q_p_4_0*tmp_24;
      real_t tmp_55 = q_p_5_0*tmp_26;
      real_t tmp_56 = q_p_6_0*tmp_28;
      real_t tmp_57 = q_p_7_0*tmp_30;
      real_t tmp_58 = q_p_8_0*tmp_32;
      real_t tmp_59 = q_p_9_0*tmp_34;
      real_t tmp_60 = q_p_0_1*tmp_49 + q_p_10_1*tmp_50 + q_p_1_1*tmp_51 + q_p_2_1*tmp_52 + q_p_3_1*tmp_53 + q_p_4_1*tmp_54 + q_p_5_1*tmp_55 + q_p_6_1*tmp_56 + q_p_7_1*tmp_57 + q_p_8_1*tmp_58 + q_p_9_1*tmp_59;
      real_t tmp_61 = q_p_0_2*tmp_49 + q_p_10_2*tmp_50 + q_p_1_2*tmp_51 + q_p_2_2*tmp_52 + q_p_3_2*tmp_53 + q_p_4_2*tmp_54 + q_p_5_2*tmp_55 + q_p_6_2*tmp_56 + q_p_7_2*tmp_57 + q_p_8_2*tmp_58 + q_p_9_2*tmp_59;
      real_t tmp_62 = q_p_0_1*q_p_0_2*tmp_14 + q_p_10_1*q_p_10_2*tmp_18 + q_p_1_1*q_p_1_2*tmp_16 + q_p_2_1*q_p_2_2*tmp_20 + q_p_3_1*q_p_3_2*tmp_22 + q_p_4_1*q_p_4_2*tmp_24 + q_p_5_1*q_p_5_2*tmp_26 + q_p_6_1*q_p_6_2*tmp_28 + q_p_7_1*q_p_7_2*tmp_30 + q_p_8_1*q_p_8_2*tmp_32 + q_p_9_1*q_p_9_2*tmp_34;
      real_t a_0_0 = (tmp_0*tmp_0)*tmp_14 + (tmp_15*tmp_15)*tmp_16 + (tmp_17*tmp_17)*tmp_18 + (tmp_19*tmp_19)*tmp_20 + (tmp_21*tmp_21)*tmp_22 + (tmp_23*tmp_23)*tmp_24 + (tmp_25*tmp_25)*tmp_26 + (tmp_27*tmp_27)*tmp_28 + (tmp_29*tmp_29)*tmp_30 + (tmp_31*tmp_31)*tmp_32 + (tmp_33*tmp_33)*tmp_34;
      real_t a_0_1 = tmp_46;
      real_t a_0_2 = tmp_47;
      real_t a_0_3 = tmp_48;
      real_t a_1_0 = tmp_46;
      real_t a_1_1 = (q_p_0_0*q_p_0_0)*tmp_14 + (q_p_10_0*q_p_10_0)*tmp_18 + (q_p_1_0*q_p_1_0)*tmp_16 + (q_p_2_0*q_p_2_0)*tmp_20 + (q_p_3_0*q_p_3_0)*tmp_22 + (q_p_4_0*q_p_4_0)*tmp_24 + (q_p_5_0*q_p_5_0)*tmp_26 + (q_p_6_0*q_p_6_0)*tmp_28 + (q_p_7_0*q_p_7_0)*tmp_30 + (q_p_8_0*q_p_8_0)*tmp_32 + (q_p_9_0*q_p_9_0)*tmp_34;
      real_t a_1_2 = tmp_60;
      real_t a_1_3 = tmp_61;
      real_t a_2_0 = tmp_47;
      real_t a_2_1 = tmp_60;
      real_t a_2_2 = (q_p_0_1*q_p_0_1)*tmp_14 + (q_p_10_1*q_p_10_1)*tmp_18 + (q_p_1_1*q_p_1_1)*tmp_16 + (q_p_2_1*q_p_2_1)*tmp_20 + (q_p_3_1*q_p_3_1)*tmp_22 + (q_p_4_1*q_p_4_1)*tmp_24 + (q_p_5_1*q_p_5_1)*tmp_26 + (q_p_6_1*q_p_6_1)*tmp_28 + (q_p_7_1*q_p_7_1)*tmp_30 + (q_p_8_1*q_p_8_1)*tmp_32 + (q_p_9_1*q_p_9_1)*tmp_34;
      real_t a_2_3 = tmp_62;
      real_t a_3_0 = tmp_48;
      real_t a_3_1 = tmp_61;
      real_t a_3_2 = tmp_62;
      real_t a_3_3 = (q_p_0_2*q_p_0_2)*tmp_14 + (q_p_10_2*q_p_10_2)*tmp_18 + (q_p_1_2*q_p_1_2)*tmp_16 + (q_p_2_2*q_p_2_2)*tmp_20 + (q_p_3_2*q_p_3_2)*tmp_22 + (q_p_4_2*q_p_4_2)*tmp_24 + (q_p_5_2*q_p_5_2)*tmp_26 + (q_p_6_2*q_p_6_2)*tmp_28 + (q_p_7_2*q_p_7_2)*tmp_30 + (q_p_8_2*q_p_8_2)*tmp_32 + (q_p_9_2*q_p_9_2)*tmp_34;
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

   void p1_mass_blending_q4::Blending_DF_Triangle( real_t in_0, real_t in_1, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3 ) const
   {
      Point3D  mappedPt( {in_0, in_1, 0} );
      Matrix2r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 1, 0 );
      *out_3 = DPsi( 1, 1 );
   }

   void p1_mass_blending_q4::Blending_DF_Tetrahedron( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
   {
      Point3D  mappedPt( {in_0, in_1, in_2} );
      Matrix3r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 0, 2 );
      *out_3 = DPsi( 1, 0 );
      *out_4 = DPsi( 1, 1 );
      *out_5 = DPsi( 1, 2 );
      *out_6 = DPsi( 2, 0 );
      *out_7 = DPsi( 2, 1 );
      *out_8 = DPsi( 2, 2 );
   }

} // namespace forms
} // namespace hyteg

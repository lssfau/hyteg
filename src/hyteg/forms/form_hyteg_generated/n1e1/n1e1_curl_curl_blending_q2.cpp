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
 * Avoid modifying this file. If buggy, consider fixing the generator itself.
 */

#include "n1e1_curl_curl_blending_q2.hpp"

namespace hyteg {
namespace forms {

   void n1e1_curl_curl_blending_q2::integrateAll( const std::array< Point3D, 4 >& coords, const std::array< int, 6 >& edgeDirections, Matrix< real_t, 6, 6 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id1 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id1 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id1 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id1 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id1 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id1 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id1 = 0;
      real_t Blending_DF_Tetrahedron_blend_out0_id2 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id2 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id2 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id2 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id2 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id2 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id2 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id2 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id2 = 0;
      real_t Blending_DF_Tetrahedron_blend_out0_id3 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id3 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id3 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id3 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id3 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id3 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id3 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id3 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id3 = 0;
      Blending_DF_Tetrahedron_blend( 0.13819660112501042*p_affine_0_0 + 0.13819660112501059*p_affine_1_0 + 0.13819660112501059*p_affine_2_0 + 0.58541019662496829*p_affine_3_0, 0.13819660112501042*p_affine_0_1 + 0.13819660112501059*p_affine_1_1 + 0.13819660112501059*p_affine_2_1 + 0.58541019662496829*p_affine_3_1, 0.13819660112501042*p_affine_0_2 + 0.13819660112501059*p_affine_1_2 + 0.13819660112501059*p_affine_2_2 + 0.58541019662496829*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_DF_Tetrahedron_blend( 0.13819660112501048*p_affine_0_0 + 0.13819660112501059*p_affine_1_0 + 0.58541019662496829*p_affine_2_0 + 0.13819660112501059*p_affine_3_0, 0.13819660112501048*p_affine_0_1 + 0.13819660112501059*p_affine_1_1 + 0.58541019662496829*p_affine_2_1 + 0.13819660112501059*p_affine_3_1, 0.13819660112501048*p_affine_0_2 + 0.13819660112501059*p_affine_1_2 + 0.58541019662496829*p_affine_2_2 + 0.13819660112501059*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id1, &Blending_DF_Tetrahedron_blend_out1_id1, &Blending_DF_Tetrahedron_blend_out2_id1, &Blending_DF_Tetrahedron_blend_out3_id1, &Blending_DF_Tetrahedron_blend_out4_id1, &Blending_DF_Tetrahedron_blend_out5_id1, &Blending_DF_Tetrahedron_blend_out6_id1, &Blending_DF_Tetrahedron_blend_out7_id1, &Blending_DF_Tetrahedron_blend_out8_id1 );
      Blending_DF_Tetrahedron_blend( 0.13819660112501053*p_affine_0_0 + 0.58541019662496829*p_affine_1_0 + 0.13819660112501059*p_affine_2_0 + 0.13819660112501059*p_affine_3_0, 0.13819660112501053*p_affine_0_1 + 0.58541019662496829*p_affine_1_1 + 0.13819660112501059*p_affine_2_1 + 0.13819660112501059*p_affine_3_1, 0.13819660112501053*p_affine_0_2 + 0.58541019662496829*p_affine_1_2 + 0.13819660112501059*p_affine_2_2 + 0.13819660112501059*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id2, &Blending_DF_Tetrahedron_blend_out1_id2, &Blending_DF_Tetrahedron_blend_out2_id2, &Blending_DF_Tetrahedron_blend_out3_id2, &Blending_DF_Tetrahedron_blend_out4_id2, &Blending_DF_Tetrahedron_blend_out5_id2, &Blending_DF_Tetrahedron_blend_out6_id2, &Blending_DF_Tetrahedron_blend_out7_id2, &Blending_DF_Tetrahedron_blend_out8_id2 );
      Blending_DF_Tetrahedron_blend( 0.58541019662496807*p_affine_0_0 + 0.13819660112501059*p_affine_1_0 + 0.13819660112501059*p_affine_2_0 + 0.13819660112501059*p_affine_3_0, 0.58541019662496807*p_affine_0_1 + 0.13819660112501059*p_affine_1_1 + 0.13819660112501059*p_affine_2_1 + 0.13819660112501059*p_affine_3_1, 0.58541019662496807*p_affine_0_2 + 0.13819660112501059*p_affine_1_2 + 0.13819660112501059*p_affine_2_2 + 0.13819660112501059*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id3, &Blending_DF_Tetrahedron_blend_out1_id3, &Blending_DF_Tetrahedron_blend_out2_id3, &Blending_DF_Tetrahedron_blend_out3_id3, &Blending_DF_Tetrahedron_blend_out4_id3, &Blending_DF_Tetrahedron_blend_out5_id3, &Blending_DF_Tetrahedron_blend_out6_id3, &Blending_DF_Tetrahedron_blend_out7_id3, &Blending_DF_Tetrahedron_blend_out8_id3 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_1_1 + tmp_2;
      real_t tmp_4 = -p_affine_0_2;
      real_t tmp_5 = p_affine_1_2 + tmp_4;
      real_t tmp_6 = Blending_DF_Tetrahedron_blend_out0_id0*tmp_1 + Blending_DF_Tetrahedron_blend_out1_id0*tmp_3 + Blending_DF_Tetrahedron_blend_out2_id0*tmp_5;
      real_t tmp_7 = 4*(walberla::real_c(edgeDirections[0])*walberla::real_c(edgeDirections[0]));
      real_t tmp_8 = Blending_DF_Tetrahedron_blend_out3_id0*tmp_1 + Blending_DF_Tetrahedron_blend_out4_id0*tmp_3 + Blending_DF_Tetrahedron_blend_out5_id0*tmp_5;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_1 + Blending_DF_Tetrahedron_blend_out7_id0*tmp_3 + Blending_DF_Tetrahedron_blend_out8_id0*tmp_5;
      real_t tmp_10 = p_affine_0_0*p_affine_1_1;
      real_t tmp_11 = p_affine_0_0*p_affine_1_2;
      real_t tmp_12 = p_affine_2_1*p_affine_3_2;
      real_t tmp_13 = p_affine_0_1*p_affine_1_0;
      real_t tmp_14 = p_affine_0_1*p_affine_1_2;
      real_t tmp_15 = p_affine_2_2*p_affine_3_0;
      real_t tmp_16 = p_affine_0_2*p_affine_1_0;
      real_t tmp_17 = p_affine_0_2*p_affine_1_1;
      real_t tmp_18 = p_affine_2_0*p_affine_3_1;
      real_t tmp_19 = p_affine_2_2*p_affine_3_1;
      real_t tmp_20 = p_affine_2_0*p_affine_3_2;
      real_t tmp_21 = p_affine_2_1*p_affine_3_0;
      real_t tmp_22 = 1.0 / (std::abs(p_affine_0_0*tmp_12 - p_affine_0_0*tmp_19 + p_affine_0_1*tmp_15 - p_affine_0_1*tmp_20 + p_affine_0_2*tmp_18 - p_affine_0_2*tmp_21 - p_affine_1_0*tmp_12 + p_affine_1_0*tmp_19 - p_affine_1_1*tmp_15 + p_affine_1_1*tmp_20 - p_affine_1_2*tmp_18 + p_affine_1_2*tmp_21 + p_affine_2_0*tmp_14 - p_affine_2_0*tmp_17 - p_affine_2_1*tmp_11 + p_affine_2_1*tmp_16 + p_affine_2_2*tmp_10 - p_affine_2_2*tmp_13 - p_affine_3_0*tmp_14 + p_affine_3_0*tmp_17 + p_affine_3_1*tmp_11 - p_affine_3_1*tmp_16 - p_affine_3_2*tmp_10 + p_affine_3_2*tmp_13));
      real_t tmp_23 = 0.041666666666666657*tmp_22/std::abs(Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0 + Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_24 = Blending_DF_Tetrahedron_blend_out0_id1*tmp_1 + Blending_DF_Tetrahedron_blend_out1_id1*tmp_3 + Blending_DF_Tetrahedron_blend_out2_id1*tmp_5;
      real_t tmp_25 = Blending_DF_Tetrahedron_blend_out3_id1*tmp_1 + Blending_DF_Tetrahedron_blend_out4_id1*tmp_3 + Blending_DF_Tetrahedron_blend_out5_id1*tmp_5;
      real_t tmp_26 = Blending_DF_Tetrahedron_blend_out6_id1*tmp_1 + Blending_DF_Tetrahedron_blend_out7_id1*tmp_3 + Blending_DF_Tetrahedron_blend_out8_id1*tmp_5;
      real_t tmp_27 = 0.041666666666666657*tmp_22/std::abs(Blending_DF_Tetrahedron_blend_out0_id1*Blending_DF_Tetrahedron_blend_out4_id1*Blending_DF_Tetrahedron_blend_out8_id1 - Blending_DF_Tetrahedron_blend_out0_id1*Blending_DF_Tetrahedron_blend_out5_id1*Blending_DF_Tetrahedron_blend_out7_id1 - Blending_DF_Tetrahedron_blend_out1_id1*Blending_DF_Tetrahedron_blend_out3_id1*Blending_DF_Tetrahedron_blend_out8_id1 + Blending_DF_Tetrahedron_blend_out1_id1*Blending_DF_Tetrahedron_blend_out5_id1*Blending_DF_Tetrahedron_blend_out6_id1 + Blending_DF_Tetrahedron_blend_out2_id1*Blending_DF_Tetrahedron_blend_out3_id1*Blending_DF_Tetrahedron_blend_out7_id1 - Blending_DF_Tetrahedron_blend_out2_id1*Blending_DF_Tetrahedron_blend_out4_id1*Blending_DF_Tetrahedron_blend_out6_id1);
      real_t tmp_28 = Blending_DF_Tetrahedron_blend_out0_id2*tmp_1 + Blending_DF_Tetrahedron_blend_out1_id2*tmp_3 + Blending_DF_Tetrahedron_blend_out2_id2*tmp_5;
      real_t tmp_29 = Blending_DF_Tetrahedron_blend_out3_id2*tmp_1 + Blending_DF_Tetrahedron_blend_out4_id2*tmp_3 + Blending_DF_Tetrahedron_blend_out5_id2*tmp_5;
      real_t tmp_30 = Blending_DF_Tetrahedron_blend_out6_id2*tmp_1 + Blending_DF_Tetrahedron_blend_out7_id2*tmp_3 + Blending_DF_Tetrahedron_blend_out8_id2*tmp_5;
      real_t tmp_31 = 0.041666666666666657*tmp_22/std::abs(Blending_DF_Tetrahedron_blend_out0_id2*Blending_DF_Tetrahedron_blend_out4_id2*Blending_DF_Tetrahedron_blend_out8_id2 - Blending_DF_Tetrahedron_blend_out0_id2*Blending_DF_Tetrahedron_blend_out5_id2*Blending_DF_Tetrahedron_blend_out7_id2 - Blending_DF_Tetrahedron_blend_out1_id2*Blending_DF_Tetrahedron_blend_out3_id2*Blending_DF_Tetrahedron_blend_out8_id2 + Blending_DF_Tetrahedron_blend_out1_id2*Blending_DF_Tetrahedron_blend_out5_id2*Blending_DF_Tetrahedron_blend_out6_id2 + Blending_DF_Tetrahedron_blend_out2_id2*Blending_DF_Tetrahedron_blend_out3_id2*Blending_DF_Tetrahedron_blend_out7_id2 - Blending_DF_Tetrahedron_blend_out2_id2*Blending_DF_Tetrahedron_blend_out4_id2*Blending_DF_Tetrahedron_blend_out6_id2);
      real_t tmp_32 = Blending_DF_Tetrahedron_blend_out0_id3*tmp_1 + Blending_DF_Tetrahedron_blend_out1_id3*tmp_3 + Blending_DF_Tetrahedron_blend_out2_id3*tmp_5;
      real_t tmp_33 = Blending_DF_Tetrahedron_blend_out3_id3*tmp_1 + Blending_DF_Tetrahedron_blend_out4_id3*tmp_3 + Blending_DF_Tetrahedron_blend_out5_id3*tmp_5;
      real_t tmp_34 = Blending_DF_Tetrahedron_blend_out6_id3*tmp_1 + Blending_DF_Tetrahedron_blend_out7_id3*tmp_3 + Blending_DF_Tetrahedron_blend_out8_id3*tmp_5;
      real_t tmp_35 = 0.041666666666666657*tmp_22/std::abs(Blending_DF_Tetrahedron_blend_out0_id3*Blending_DF_Tetrahedron_blend_out4_id3*Blending_DF_Tetrahedron_blend_out8_id3 - Blending_DF_Tetrahedron_blend_out0_id3*Blending_DF_Tetrahedron_blend_out5_id3*Blending_DF_Tetrahedron_blend_out7_id3 - Blending_DF_Tetrahedron_blend_out1_id3*Blending_DF_Tetrahedron_blend_out3_id3*Blending_DF_Tetrahedron_blend_out8_id3 + Blending_DF_Tetrahedron_blend_out1_id3*Blending_DF_Tetrahedron_blend_out5_id3*Blending_DF_Tetrahedron_blend_out6_id3 + Blending_DF_Tetrahedron_blend_out2_id3*Blending_DF_Tetrahedron_blend_out3_id3*Blending_DF_Tetrahedron_blend_out7_id3 - Blending_DF_Tetrahedron_blend_out2_id3*Blending_DF_Tetrahedron_blend_out4_id3*Blending_DF_Tetrahedron_blend_out6_id3);
      real_t tmp_36 = p_affine_2_0 + tmp_0;
      real_t tmp_37 = p_affine_2_1 + tmp_2;
      real_t tmp_38 = p_affine_2_2 + tmp_4;
      real_t tmp_39 = Blending_DF_Tetrahedron_blend_out0_id0*tmp_36 + Blending_DF_Tetrahedron_blend_out1_id0*tmp_37 + Blending_DF_Tetrahedron_blend_out2_id0*tmp_38;
      real_t tmp_40 = 4*walberla::real_c(edgeDirections[0]);
      real_t tmp_41 = tmp_40*walberla::real_c(edgeDirections[1]);
      real_t tmp_42 = Blending_DF_Tetrahedron_blend_out3_id0*tmp_36 + Blending_DF_Tetrahedron_blend_out4_id0*tmp_37 + Blending_DF_Tetrahedron_blend_out5_id0*tmp_38;
      real_t tmp_43 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_36 + Blending_DF_Tetrahedron_blend_out7_id0*tmp_37 + Blending_DF_Tetrahedron_blend_out8_id0*tmp_38;
      real_t tmp_44 = Blending_DF_Tetrahedron_blend_out0_id1*tmp_36 + Blending_DF_Tetrahedron_blend_out1_id1*tmp_37 + Blending_DF_Tetrahedron_blend_out2_id1*tmp_38;
      real_t tmp_45 = Blending_DF_Tetrahedron_blend_out3_id1*tmp_36 + Blending_DF_Tetrahedron_blend_out4_id1*tmp_37 + Blending_DF_Tetrahedron_blend_out5_id1*tmp_38;
      real_t tmp_46 = Blending_DF_Tetrahedron_blend_out6_id1*tmp_36 + Blending_DF_Tetrahedron_blend_out7_id1*tmp_37 + Blending_DF_Tetrahedron_blend_out8_id1*tmp_38;
      real_t tmp_47 = Blending_DF_Tetrahedron_blend_out0_id2*tmp_36 + Blending_DF_Tetrahedron_blend_out1_id2*tmp_37 + Blending_DF_Tetrahedron_blend_out2_id2*tmp_38;
      real_t tmp_48 = Blending_DF_Tetrahedron_blend_out3_id2*tmp_36 + Blending_DF_Tetrahedron_blend_out4_id2*tmp_37 + Blending_DF_Tetrahedron_blend_out5_id2*tmp_38;
      real_t tmp_49 = Blending_DF_Tetrahedron_blend_out6_id2*tmp_36 + Blending_DF_Tetrahedron_blend_out7_id2*tmp_37 + Blending_DF_Tetrahedron_blend_out8_id2*tmp_38;
      real_t tmp_50 = Blending_DF_Tetrahedron_blend_out0_id3*tmp_36 + Blending_DF_Tetrahedron_blend_out1_id3*tmp_37 + Blending_DF_Tetrahedron_blend_out2_id3*tmp_38;
      real_t tmp_51 = Blending_DF_Tetrahedron_blend_out3_id3*tmp_36 + Blending_DF_Tetrahedron_blend_out4_id3*tmp_37 + Blending_DF_Tetrahedron_blend_out5_id3*tmp_38;
      real_t tmp_52 = Blending_DF_Tetrahedron_blend_out6_id3*tmp_36 + Blending_DF_Tetrahedron_blend_out7_id3*tmp_37 + Blending_DF_Tetrahedron_blend_out8_id3*tmp_38;
      real_t tmp_53 = tmp_23*(-tmp_39*tmp_41*tmp_6 - tmp_41*tmp_42*tmp_8 - tmp_41*tmp_43*tmp_9) + tmp_27*(-tmp_24*tmp_41*tmp_44 - tmp_25*tmp_41*tmp_45 - tmp_26*tmp_41*tmp_46) + tmp_31*(-tmp_28*tmp_41*tmp_47 - tmp_29*tmp_41*tmp_48 - tmp_30*tmp_41*tmp_49) + tmp_35*(-tmp_32*tmp_41*tmp_50 - tmp_33*tmp_41*tmp_51 - tmp_34*tmp_41*tmp_52);
      real_t tmp_54 = p_affine_3_0 + tmp_0;
      real_t tmp_55 = p_affine_3_1 + tmp_2;
      real_t tmp_56 = p_affine_3_2 + tmp_4;
      real_t tmp_57 = Blending_DF_Tetrahedron_blend_out0_id0*tmp_54 + Blending_DF_Tetrahedron_blend_out1_id0*tmp_55 + Blending_DF_Tetrahedron_blend_out2_id0*tmp_56;
      real_t tmp_58 = tmp_40*walberla::real_c(edgeDirections[2]);
      real_t tmp_59 = Blending_DF_Tetrahedron_blend_out3_id0*tmp_54 + Blending_DF_Tetrahedron_blend_out4_id0*tmp_55 + Blending_DF_Tetrahedron_blend_out5_id0*tmp_56;
      real_t tmp_60 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_54 + Blending_DF_Tetrahedron_blend_out7_id0*tmp_55 + Blending_DF_Tetrahedron_blend_out8_id0*tmp_56;
      real_t tmp_61 = Blending_DF_Tetrahedron_blend_out0_id1*tmp_54 + Blending_DF_Tetrahedron_blend_out1_id1*tmp_55 + Blending_DF_Tetrahedron_blend_out2_id1*tmp_56;
      real_t tmp_62 = Blending_DF_Tetrahedron_blend_out3_id1*tmp_54 + Blending_DF_Tetrahedron_blend_out4_id1*tmp_55 + Blending_DF_Tetrahedron_blend_out5_id1*tmp_56;
      real_t tmp_63 = Blending_DF_Tetrahedron_blend_out6_id1*tmp_54 + Blending_DF_Tetrahedron_blend_out7_id1*tmp_55 + Blending_DF_Tetrahedron_blend_out8_id1*tmp_56;
      real_t tmp_64 = Blending_DF_Tetrahedron_blend_out0_id2*tmp_54 + Blending_DF_Tetrahedron_blend_out1_id2*tmp_55 + Blending_DF_Tetrahedron_blend_out2_id2*tmp_56;
      real_t tmp_65 = Blending_DF_Tetrahedron_blend_out3_id2*tmp_54 + Blending_DF_Tetrahedron_blend_out4_id2*tmp_55 + Blending_DF_Tetrahedron_blend_out5_id2*tmp_56;
      real_t tmp_66 = Blending_DF_Tetrahedron_blend_out6_id2*tmp_54 + Blending_DF_Tetrahedron_blend_out7_id2*tmp_55 + Blending_DF_Tetrahedron_blend_out8_id2*tmp_56;
      real_t tmp_67 = Blending_DF_Tetrahedron_blend_out0_id3*tmp_54 + Blending_DF_Tetrahedron_blend_out1_id3*tmp_55 + Blending_DF_Tetrahedron_blend_out2_id3*tmp_56;
      real_t tmp_68 = Blending_DF_Tetrahedron_blend_out3_id3*tmp_54 + Blending_DF_Tetrahedron_blend_out4_id3*tmp_55 + Blending_DF_Tetrahedron_blend_out5_id3*tmp_56;
      real_t tmp_69 = Blending_DF_Tetrahedron_blend_out6_id3*tmp_54 + Blending_DF_Tetrahedron_blend_out7_id3*tmp_55 + Blending_DF_Tetrahedron_blend_out8_id3*tmp_56;
      real_t tmp_70 = tmp_23*(tmp_57*tmp_58*tmp_6 + tmp_58*tmp_59*tmp_8 + tmp_58*tmp_60*tmp_9) + tmp_27*(tmp_24*tmp_58*tmp_61 + tmp_25*tmp_58*tmp_62 + tmp_26*tmp_58*tmp_63) + tmp_31*(tmp_28*tmp_58*tmp_64 + tmp_29*tmp_58*tmp_65 + tmp_30*tmp_58*tmp_66) + tmp_35*(tmp_32*tmp_58*tmp_67 + tmp_33*tmp_58*tmp_68 + tmp_34*tmp_58*tmp_69);
      real_t tmp_71 = 2*walberla::real_c(edgeDirections[3]);
      real_t tmp_72 = tmp_39*tmp_71 - tmp_6*tmp_71;
      real_t tmp_73 = 2*walberla::real_c(edgeDirections[0]);
      real_t tmp_74 = tmp_6*tmp_73;
      real_t tmp_75 = tmp_42*tmp_71 - tmp_71*tmp_8;
      real_t tmp_76 = tmp_73*tmp_8;
      real_t tmp_77 = tmp_43*tmp_71 - tmp_71*tmp_9;
      real_t tmp_78 = tmp_73*tmp_9;
      real_t tmp_79 = -tmp_24*tmp_71 + tmp_44*tmp_71;
      real_t tmp_80 = tmp_24*tmp_73;
      real_t tmp_81 = -tmp_25*tmp_71 + tmp_45*tmp_71;
      real_t tmp_82 = tmp_25*tmp_73;
      real_t tmp_83 = -tmp_26*tmp_71 + tmp_46*tmp_71;
      real_t tmp_84 = tmp_26*tmp_73;
      real_t tmp_85 = -tmp_28*tmp_71 + tmp_47*tmp_71;
      real_t tmp_86 = tmp_28*tmp_73;
      real_t tmp_87 = -tmp_29*tmp_71 + tmp_48*tmp_71;
      real_t tmp_88 = tmp_29*tmp_73;
      real_t tmp_89 = -tmp_30*tmp_71 + tmp_49*tmp_71;
      real_t tmp_90 = tmp_30*tmp_73;
      real_t tmp_91 = -tmp_32*tmp_71 + tmp_50*tmp_71;
      real_t tmp_92 = tmp_32*tmp_73;
      real_t tmp_93 = -tmp_33*tmp_71 + tmp_51*tmp_71;
      real_t tmp_94 = tmp_33*tmp_73;
      real_t tmp_95 = -tmp_34*tmp_71 + tmp_52*tmp_71;
      real_t tmp_96 = tmp_34*tmp_73;
      real_t tmp_97 = tmp_23*(tmp_72*tmp_74 + tmp_75*tmp_76 + tmp_77*tmp_78) + tmp_27*(tmp_79*tmp_80 + tmp_81*tmp_82 + tmp_83*tmp_84) + tmp_31*(tmp_85*tmp_86 + tmp_87*tmp_88 + tmp_89*tmp_90) + tmp_35*(tmp_91*tmp_92 + tmp_93*tmp_94 + tmp_95*tmp_96);
      real_t tmp_98 = 2*walberla::real_c(edgeDirections[4]);
      real_t tmp_99 = -tmp_57*tmp_98 + tmp_6*tmp_98;
      real_t tmp_100 = -tmp_59*tmp_98 + tmp_8*tmp_98;
      real_t tmp_101 = -tmp_60*tmp_98 + tmp_9*tmp_98;
      real_t tmp_102 = tmp_24*tmp_98 - tmp_61*tmp_98;
      real_t tmp_103 = tmp_25*tmp_98 - tmp_62*tmp_98;
      real_t tmp_104 = tmp_26*tmp_98 - tmp_63*tmp_98;
      real_t tmp_105 = tmp_28*tmp_98 - tmp_64*tmp_98;
      real_t tmp_106 = tmp_29*tmp_98 - tmp_65*tmp_98;
      real_t tmp_107 = tmp_30*tmp_98 - tmp_66*tmp_98;
      real_t tmp_108 = tmp_32*tmp_98 - tmp_67*tmp_98;
      real_t tmp_109 = tmp_33*tmp_98 - tmp_68*tmp_98;
      real_t tmp_110 = tmp_34*tmp_98 - tmp_69*tmp_98;
      real_t tmp_111 = tmp_23*(tmp_100*tmp_76 + tmp_101*tmp_78 + tmp_74*tmp_99) + tmp_27*(tmp_102*tmp_80 + tmp_103*tmp_82 + tmp_104*tmp_84) + tmp_31*(tmp_105*tmp_86 + tmp_106*tmp_88 + tmp_107*tmp_90) + tmp_35*(tmp_108*tmp_92 + tmp_109*tmp_94 + tmp_110*tmp_96);
      real_t tmp_112 = 2*walberla::real_c(edgeDirections[5]);
      real_t tmp_113 = -tmp_112*tmp_39 + tmp_112*tmp_57;
      real_t tmp_114 = -tmp_112*tmp_42 + tmp_112*tmp_59;
      real_t tmp_115 = -tmp_112*tmp_43 + tmp_112*tmp_60;
      real_t tmp_116 = -tmp_112*tmp_44 + tmp_112*tmp_61;
      real_t tmp_117 = -tmp_112*tmp_45 + tmp_112*tmp_62;
      real_t tmp_118 = -tmp_112*tmp_46 + tmp_112*tmp_63;
      real_t tmp_119 = -tmp_112*tmp_47 + tmp_112*tmp_64;
      real_t tmp_120 = -tmp_112*tmp_48 + tmp_112*tmp_65;
      real_t tmp_121 = -tmp_112*tmp_49 + tmp_112*tmp_66;
      real_t tmp_122 = -tmp_112*tmp_50 + tmp_112*tmp_67;
      real_t tmp_123 = -tmp_112*tmp_51 + tmp_112*tmp_68;
      real_t tmp_124 = -tmp_112*tmp_52 + tmp_112*tmp_69;
      real_t tmp_125 = tmp_23*(tmp_113*tmp_74 + tmp_114*tmp_76 + tmp_115*tmp_78) + tmp_27*(tmp_116*tmp_80 + tmp_117*tmp_82 + tmp_118*tmp_84) + tmp_31*(tmp_119*tmp_86 + tmp_120*tmp_88 + tmp_121*tmp_90) + tmp_35*(tmp_122*tmp_92 + tmp_123*tmp_94 + tmp_124*tmp_96);
      real_t tmp_126 = 4*(walberla::real_c(edgeDirections[1])*walberla::real_c(edgeDirections[1]));
      real_t tmp_127 = 4*walberla::real_c(edgeDirections[1])*walberla::real_c(edgeDirections[2]);
      real_t tmp_128 = tmp_23*(-tmp_127*tmp_39*tmp_57 - tmp_127*tmp_42*tmp_59 - tmp_127*tmp_43*tmp_60) + tmp_27*(-tmp_127*tmp_44*tmp_61 - tmp_127*tmp_45*tmp_62 - tmp_127*tmp_46*tmp_63) + tmp_31*(-tmp_127*tmp_47*tmp_64 - tmp_127*tmp_48*tmp_65 - tmp_127*tmp_49*tmp_66) + tmp_35*(-tmp_127*tmp_50*tmp_67 - tmp_127*tmp_51*tmp_68 - tmp_127*tmp_52*tmp_69);
      real_t tmp_129 = 2*walberla::real_c(edgeDirections[1]);
      real_t tmp_130 = tmp_129*tmp_39;
      real_t tmp_131 = tmp_129*tmp_42;
      real_t tmp_132 = tmp_129*tmp_43;
      real_t tmp_133 = tmp_129*tmp_44;
      real_t tmp_134 = tmp_129*tmp_45;
      real_t tmp_135 = tmp_129*tmp_46;
      real_t tmp_136 = tmp_129*tmp_47;
      real_t tmp_137 = tmp_129*tmp_48;
      real_t tmp_138 = tmp_129*tmp_49;
      real_t tmp_139 = tmp_129*tmp_50;
      real_t tmp_140 = tmp_129*tmp_51;
      real_t tmp_141 = tmp_129*tmp_52;
      real_t tmp_142 = tmp_23*(-tmp_130*tmp_72 - tmp_131*tmp_75 - tmp_132*tmp_77) + tmp_27*(-tmp_133*tmp_79 - tmp_134*tmp_81 - tmp_135*tmp_83) + tmp_31*(-tmp_136*tmp_85 - tmp_137*tmp_87 - tmp_138*tmp_89) + tmp_35*(-tmp_139*tmp_91 - tmp_140*tmp_93 - tmp_141*tmp_95);
      real_t tmp_143 = tmp_23*(-tmp_100*tmp_131 - tmp_101*tmp_132 - tmp_130*tmp_99) + tmp_27*(-tmp_102*tmp_133 - tmp_103*tmp_134 - tmp_104*tmp_135) + tmp_31*(-tmp_105*tmp_136 - tmp_106*tmp_137 - tmp_107*tmp_138) + tmp_35*(-tmp_108*tmp_139 - tmp_109*tmp_140 - tmp_110*tmp_141);
      real_t tmp_144 = tmp_23*(-tmp_113*tmp_130 - tmp_114*tmp_131 - tmp_115*tmp_132) + tmp_27*(-tmp_116*tmp_133 - tmp_117*tmp_134 - tmp_118*tmp_135) + tmp_31*(-tmp_119*tmp_136 - tmp_120*tmp_137 - tmp_121*tmp_138) + tmp_35*(-tmp_122*tmp_139 - tmp_123*tmp_140 - tmp_124*tmp_141);
      real_t tmp_145 = 4*(walberla::real_c(edgeDirections[2])*walberla::real_c(edgeDirections[2]));
      real_t tmp_146 = 2*walberla::real_c(edgeDirections[2]);
      real_t tmp_147 = tmp_146*tmp_57;
      real_t tmp_148 = tmp_146*tmp_59;
      real_t tmp_149 = tmp_146*tmp_60;
      real_t tmp_150 = tmp_146*tmp_61;
      real_t tmp_151 = tmp_146*tmp_62;
      real_t tmp_152 = tmp_146*tmp_63;
      real_t tmp_153 = tmp_146*tmp_64;
      real_t tmp_154 = tmp_146*tmp_65;
      real_t tmp_155 = tmp_146*tmp_66;
      real_t tmp_156 = tmp_146*tmp_67;
      real_t tmp_157 = tmp_146*tmp_68;
      real_t tmp_158 = tmp_146*tmp_69;
      real_t tmp_159 = tmp_23*(tmp_147*tmp_72 + tmp_148*tmp_75 + tmp_149*tmp_77) + tmp_27*(tmp_150*tmp_79 + tmp_151*tmp_81 + tmp_152*tmp_83) + tmp_31*(tmp_153*tmp_85 + tmp_154*tmp_87 + tmp_155*tmp_89) + tmp_35*(tmp_156*tmp_91 + tmp_157*tmp_93 + tmp_158*tmp_95);
      real_t tmp_160 = tmp_23*(tmp_100*tmp_148 + tmp_101*tmp_149 + tmp_147*tmp_99) + tmp_27*(tmp_102*tmp_150 + tmp_103*tmp_151 + tmp_104*tmp_152) + tmp_31*(tmp_105*tmp_153 + tmp_106*tmp_154 + tmp_107*tmp_155) + tmp_35*(tmp_108*tmp_156 + tmp_109*tmp_157 + tmp_110*tmp_158);
      real_t tmp_161 = tmp_23*(tmp_113*tmp_147 + tmp_114*tmp_148 + tmp_115*tmp_149) + tmp_27*(tmp_116*tmp_150 + tmp_117*tmp_151 + tmp_118*tmp_152) + tmp_31*(tmp_119*tmp_153 + tmp_120*tmp_154 + tmp_121*tmp_155) + tmp_35*(tmp_122*tmp_156 + tmp_123*tmp_157 + tmp_124*tmp_158);
      real_t tmp_162 = tmp_23*(tmp_100*tmp_75 + tmp_101*tmp_77 + tmp_72*tmp_99) + tmp_27*(tmp_102*tmp_79 + tmp_103*tmp_81 + tmp_104*tmp_83) + tmp_31*(tmp_105*tmp_85 + tmp_106*tmp_87 + tmp_107*tmp_89) + tmp_35*(tmp_108*tmp_91 + tmp_109*tmp_93 + tmp_110*tmp_95);
      real_t tmp_163 = tmp_23*(tmp_113*tmp_72 + tmp_114*tmp_75 + tmp_115*tmp_77) + tmp_27*(tmp_116*tmp_79 + tmp_117*tmp_81 + tmp_118*tmp_83) + tmp_31*(tmp_119*tmp_85 + tmp_120*tmp_87 + tmp_121*tmp_89) + tmp_35*(tmp_122*tmp_91 + tmp_123*tmp_93 + tmp_124*tmp_95);
      real_t tmp_164 = tmp_23*(tmp_100*tmp_114 + tmp_101*tmp_115 + tmp_113*tmp_99) + tmp_27*(tmp_102*tmp_116 + tmp_103*tmp_117 + tmp_104*tmp_118) + tmp_31*(tmp_105*tmp_119 + tmp_106*tmp_120 + tmp_107*tmp_121) + tmp_35*(tmp_108*tmp_122 + tmp_109*tmp_123 + tmp_110*tmp_124);
      real_t a_0_0 = tmp_23*((tmp_6*tmp_6)*tmp_7 + tmp_7*(tmp_8*tmp_8) + tmp_7*(tmp_9*tmp_9)) + tmp_27*((tmp_24*tmp_24)*tmp_7 + (tmp_25*tmp_25)*tmp_7 + (tmp_26*tmp_26)*tmp_7) + tmp_31*((tmp_28*tmp_28)*tmp_7 + (tmp_29*tmp_29)*tmp_7 + (tmp_30*tmp_30)*tmp_7) + tmp_35*((tmp_32*tmp_32)*tmp_7 + (tmp_33*tmp_33)*tmp_7 + (tmp_34*tmp_34)*tmp_7);
      real_t a_0_1 = tmp_53;
      real_t a_0_2 = tmp_70;
      real_t a_0_3 = tmp_97;
      real_t a_0_4 = tmp_111;
      real_t a_0_5 = tmp_125;
      real_t a_1_0 = tmp_53;
      real_t a_1_1 = tmp_23*(tmp_126*(tmp_39*tmp_39) + tmp_126*(tmp_42*tmp_42) + tmp_126*(tmp_43*tmp_43)) + tmp_27*(tmp_126*(tmp_44*tmp_44) + tmp_126*(tmp_45*tmp_45) + tmp_126*(tmp_46*tmp_46)) + tmp_31*(tmp_126*(tmp_47*tmp_47) + tmp_126*(tmp_48*tmp_48) + tmp_126*(tmp_49*tmp_49)) + tmp_35*(tmp_126*(tmp_50*tmp_50) + tmp_126*(tmp_51*tmp_51) + tmp_126*(tmp_52*tmp_52));
      real_t a_1_2 = tmp_128;
      real_t a_1_3 = tmp_142;
      real_t a_1_4 = tmp_143;
      real_t a_1_5 = tmp_144;
      real_t a_2_0 = tmp_70;
      real_t a_2_1 = tmp_128;
      real_t a_2_2 = tmp_23*(tmp_145*(tmp_57*tmp_57) + tmp_145*(tmp_59*tmp_59) + tmp_145*(tmp_60*tmp_60)) + tmp_27*(tmp_145*(tmp_61*tmp_61) + tmp_145*(tmp_62*tmp_62) + tmp_145*(tmp_63*tmp_63)) + tmp_31*(tmp_145*(tmp_64*tmp_64) + tmp_145*(tmp_65*tmp_65) + tmp_145*(tmp_66*tmp_66)) + tmp_35*(tmp_145*(tmp_67*tmp_67) + tmp_145*(tmp_68*tmp_68) + tmp_145*(tmp_69*tmp_69));
      real_t a_2_3 = tmp_159;
      real_t a_2_4 = tmp_160;
      real_t a_2_5 = tmp_161;
      real_t a_3_0 = tmp_97;
      real_t a_3_1 = tmp_142;
      real_t a_3_2 = tmp_159;
      real_t a_3_3 = tmp_23*((tmp_72*tmp_72) + (tmp_75*tmp_75) + (tmp_77*tmp_77)) + tmp_27*((tmp_79*tmp_79) + (tmp_81*tmp_81) + (tmp_83*tmp_83)) + tmp_31*((tmp_85*tmp_85) + (tmp_87*tmp_87) + (tmp_89*tmp_89)) + tmp_35*((tmp_91*tmp_91) + (tmp_93*tmp_93) + (tmp_95*tmp_95));
      real_t a_3_4 = tmp_162;
      real_t a_3_5 = tmp_163;
      real_t a_4_0 = tmp_111;
      real_t a_4_1 = tmp_143;
      real_t a_4_2 = tmp_160;
      real_t a_4_3 = tmp_162;
      real_t a_4_4 = tmp_23*((tmp_100*tmp_100) + (tmp_101*tmp_101) + (tmp_99*tmp_99)) + tmp_27*((tmp_102*tmp_102) + (tmp_103*tmp_103) + (tmp_104*tmp_104)) + tmp_31*((tmp_105*tmp_105) + (tmp_106*tmp_106) + (tmp_107*tmp_107)) + tmp_35*((tmp_108*tmp_108) + (tmp_109*tmp_109) + (tmp_110*tmp_110));
      real_t a_4_5 = tmp_164;
      real_t a_5_0 = tmp_125;
      real_t a_5_1 = tmp_144;
      real_t a_5_2 = tmp_161;
      real_t a_5_3 = tmp_163;
      real_t a_5_4 = tmp_164;
      real_t a_5_5 = tmp_23*((tmp_113*tmp_113) + (tmp_114*tmp_114) + (tmp_115*tmp_115)) + tmp_27*((tmp_116*tmp_116) + (tmp_117*tmp_117) + (tmp_118*tmp_118)) + tmp_31*((tmp_119*tmp_119) + (tmp_120*tmp_120) + (tmp_121*tmp_121)) + tmp_35*((tmp_122*tmp_122) + (tmp_123*tmp_123) + (tmp_124*tmp_124));
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
      (elMat(0, 4)) = a_0_4;
      (elMat(0, 5)) = a_0_5;
      (elMat(1, 0)) = a_1_0;
      (elMat(1, 1)) = a_1_1;
      (elMat(1, 2)) = a_1_2;
      (elMat(1, 3)) = a_1_3;
      (elMat(1, 4)) = a_1_4;
      (elMat(1, 5)) = a_1_5;
      (elMat(2, 0)) = a_2_0;
      (elMat(2, 1)) = a_2_1;
      (elMat(2, 2)) = a_2_2;
      (elMat(2, 3)) = a_2_3;
      (elMat(2, 4)) = a_2_4;
      (elMat(2, 5)) = a_2_5;
      (elMat(3, 0)) = a_3_0;
      (elMat(3, 1)) = a_3_1;
      (elMat(3, 2)) = a_3_2;
      (elMat(3, 3)) = a_3_3;
      (elMat(3, 4)) = a_3_4;
      (elMat(3, 5)) = a_3_5;
      (elMat(4, 0)) = a_4_0;
      (elMat(4, 1)) = a_4_1;
      (elMat(4, 2)) = a_4_2;
      (elMat(4, 3)) = a_4_3;
      (elMat(4, 4)) = a_4_4;
      (elMat(4, 5)) = a_4_5;
      (elMat(5, 0)) = a_5_0;
      (elMat(5, 1)) = a_5_1;
      (elMat(5, 2)) = a_5_2;
      (elMat(5, 3)) = a_5_3;
      (elMat(5, 4)) = a_5_4;
      (elMat(5, 5)) = a_5_5;
   }

   void n1e1_curl_curl_blending_q2::Blending_DF_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
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

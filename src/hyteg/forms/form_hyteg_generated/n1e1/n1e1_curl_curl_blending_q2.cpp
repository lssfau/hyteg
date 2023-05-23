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

   void n1e1_curl_curl_blending_q2::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 6, 6 >& elMat ) const
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
      real_t tmp_1 = 2*p_affine_1_0 + 2*tmp_0;
      real_t tmp_2 = Blending_DF_Tetrahedron_blend_out0_id0*tmp_1;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = 2*p_affine_1_1 + 2*tmp_3;
      real_t tmp_5 = Blending_DF_Tetrahedron_blend_out1_id0*tmp_4;
      real_t tmp_6 = -p_affine_0_2;
      real_t tmp_7 = 2*p_affine_1_2 + 2*tmp_6;
      real_t tmp_8 = Blending_DF_Tetrahedron_blend_out2_id0*tmp_7;
      real_t tmp_9 = tmp_2 + tmp_5 + tmp_8;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out3_id0*tmp_1;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out4_id0*tmp_4;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out5_id0*tmp_7;
      real_t tmp_13 = tmp_10 + tmp_11 + tmp_12;
      real_t tmp_14 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_1;
      real_t tmp_15 = Blending_DF_Tetrahedron_blend_out7_id0*tmp_4;
      real_t tmp_16 = Blending_DF_Tetrahedron_blend_out8_id0*tmp_7;
      real_t tmp_17 = tmp_14 + tmp_15 + tmp_16;
      real_t tmp_18 = p_affine_0_0*p_affine_1_1;
      real_t tmp_19 = p_affine_0_0*p_affine_1_2;
      real_t tmp_20 = p_affine_2_1*p_affine_3_2;
      real_t tmp_21 = p_affine_0_1*p_affine_1_0;
      real_t tmp_22 = p_affine_0_1*p_affine_1_2;
      real_t tmp_23 = p_affine_2_2*p_affine_3_0;
      real_t tmp_24 = p_affine_0_2*p_affine_1_0;
      real_t tmp_25 = p_affine_0_2*p_affine_1_1;
      real_t tmp_26 = p_affine_2_0*p_affine_3_1;
      real_t tmp_27 = p_affine_2_2*p_affine_3_1;
      real_t tmp_28 = p_affine_2_0*p_affine_3_2;
      real_t tmp_29 = p_affine_2_1*p_affine_3_0;
      real_t tmp_30 = 1.0 / (std::abs(p_affine_0_0*tmp_20 - p_affine_0_0*tmp_27 + p_affine_0_1*tmp_23 - p_affine_0_1*tmp_28 + p_affine_0_2*tmp_26 - p_affine_0_2*tmp_29 - p_affine_1_0*tmp_20 + p_affine_1_0*tmp_27 - p_affine_1_1*tmp_23 + p_affine_1_1*tmp_28 - p_affine_1_2*tmp_26 + p_affine_1_2*tmp_29 + p_affine_2_0*tmp_22 - p_affine_2_0*tmp_25 - p_affine_2_1*tmp_19 + p_affine_2_1*tmp_24 + p_affine_2_2*tmp_18 - p_affine_2_2*tmp_21 - p_affine_3_0*tmp_22 + p_affine_3_0*tmp_25 + p_affine_3_1*tmp_19 - p_affine_3_1*tmp_24 - p_affine_3_2*tmp_18 + p_affine_3_2*tmp_21));
      real_t tmp_31 = 0.041666666666666657*tmp_30/std::abs(Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0 + Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_32 = Blending_DF_Tetrahedron_blend_out0_id1*tmp_1;
      real_t tmp_33 = Blending_DF_Tetrahedron_blend_out1_id1*tmp_4;
      real_t tmp_34 = Blending_DF_Tetrahedron_blend_out2_id1*tmp_7;
      real_t tmp_35 = tmp_32 + tmp_33 + tmp_34;
      real_t tmp_36 = Blending_DF_Tetrahedron_blend_out3_id1*tmp_1;
      real_t tmp_37 = Blending_DF_Tetrahedron_blend_out4_id1*tmp_4;
      real_t tmp_38 = Blending_DF_Tetrahedron_blend_out5_id1*tmp_7;
      real_t tmp_39 = tmp_36 + tmp_37 + tmp_38;
      real_t tmp_40 = Blending_DF_Tetrahedron_blend_out6_id1*tmp_1;
      real_t tmp_41 = Blending_DF_Tetrahedron_blend_out7_id1*tmp_4;
      real_t tmp_42 = Blending_DF_Tetrahedron_blend_out8_id1*tmp_7;
      real_t tmp_43 = tmp_40 + tmp_41 + tmp_42;
      real_t tmp_44 = 0.041666666666666657*tmp_30/std::abs(Blending_DF_Tetrahedron_blend_out0_id1*Blending_DF_Tetrahedron_blend_out4_id1*Blending_DF_Tetrahedron_blend_out8_id1 - Blending_DF_Tetrahedron_blend_out0_id1*Blending_DF_Tetrahedron_blend_out5_id1*Blending_DF_Tetrahedron_blend_out7_id1 - Blending_DF_Tetrahedron_blend_out1_id1*Blending_DF_Tetrahedron_blend_out3_id1*Blending_DF_Tetrahedron_blend_out8_id1 + Blending_DF_Tetrahedron_blend_out1_id1*Blending_DF_Tetrahedron_blend_out5_id1*Blending_DF_Tetrahedron_blend_out6_id1 + Blending_DF_Tetrahedron_blend_out2_id1*Blending_DF_Tetrahedron_blend_out3_id1*Blending_DF_Tetrahedron_blend_out7_id1 - Blending_DF_Tetrahedron_blend_out2_id1*Blending_DF_Tetrahedron_blend_out4_id1*Blending_DF_Tetrahedron_blend_out6_id1);
      real_t tmp_45 = Blending_DF_Tetrahedron_blend_out0_id2*tmp_1;
      real_t tmp_46 = Blending_DF_Tetrahedron_blend_out1_id2*tmp_4;
      real_t tmp_47 = Blending_DF_Tetrahedron_blend_out2_id2*tmp_7;
      real_t tmp_48 = tmp_45 + tmp_46 + tmp_47;
      real_t tmp_49 = Blending_DF_Tetrahedron_blend_out3_id2*tmp_1;
      real_t tmp_50 = Blending_DF_Tetrahedron_blend_out4_id2*tmp_4;
      real_t tmp_51 = Blending_DF_Tetrahedron_blend_out5_id2*tmp_7;
      real_t tmp_52 = tmp_49 + tmp_50 + tmp_51;
      real_t tmp_53 = Blending_DF_Tetrahedron_blend_out6_id2*tmp_1;
      real_t tmp_54 = Blending_DF_Tetrahedron_blend_out7_id2*tmp_4;
      real_t tmp_55 = Blending_DF_Tetrahedron_blend_out8_id2*tmp_7;
      real_t tmp_56 = tmp_53 + tmp_54 + tmp_55;
      real_t tmp_57 = 0.041666666666666657*tmp_30/std::abs(Blending_DF_Tetrahedron_blend_out0_id2*Blending_DF_Tetrahedron_blend_out4_id2*Blending_DF_Tetrahedron_blend_out8_id2 - Blending_DF_Tetrahedron_blend_out0_id2*Blending_DF_Tetrahedron_blend_out5_id2*Blending_DF_Tetrahedron_blend_out7_id2 - Blending_DF_Tetrahedron_blend_out1_id2*Blending_DF_Tetrahedron_blend_out3_id2*Blending_DF_Tetrahedron_blend_out8_id2 + Blending_DF_Tetrahedron_blend_out1_id2*Blending_DF_Tetrahedron_blend_out5_id2*Blending_DF_Tetrahedron_blend_out6_id2 + Blending_DF_Tetrahedron_blend_out2_id2*Blending_DF_Tetrahedron_blend_out3_id2*Blending_DF_Tetrahedron_blend_out7_id2 - Blending_DF_Tetrahedron_blend_out2_id2*Blending_DF_Tetrahedron_blend_out4_id2*Blending_DF_Tetrahedron_blend_out6_id2);
      real_t tmp_58 = Blending_DF_Tetrahedron_blend_out0_id3*tmp_1;
      real_t tmp_59 = Blending_DF_Tetrahedron_blend_out1_id3*tmp_4;
      real_t tmp_60 = Blending_DF_Tetrahedron_blend_out2_id3*tmp_7;
      real_t tmp_61 = tmp_58 + tmp_59 + tmp_60;
      real_t tmp_62 = Blending_DF_Tetrahedron_blend_out3_id3*tmp_1;
      real_t tmp_63 = Blending_DF_Tetrahedron_blend_out4_id3*tmp_4;
      real_t tmp_64 = Blending_DF_Tetrahedron_blend_out5_id3*tmp_7;
      real_t tmp_65 = tmp_62 + tmp_63 + tmp_64;
      real_t tmp_66 = Blending_DF_Tetrahedron_blend_out6_id3*tmp_1;
      real_t tmp_67 = Blending_DF_Tetrahedron_blend_out7_id3*tmp_4;
      real_t tmp_68 = Blending_DF_Tetrahedron_blend_out8_id3*tmp_7;
      real_t tmp_69 = tmp_66 + tmp_67 + tmp_68;
      real_t tmp_70 = 0.041666666666666657*tmp_30/std::abs(Blending_DF_Tetrahedron_blend_out0_id3*Blending_DF_Tetrahedron_blend_out4_id3*Blending_DF_Tetrahedron_blend_out8_id3 - Blending_DF_Tetrahedron_blend_out0_id3*Blending_DF_Tetrahedron_blend_out5_id3*Blending_DF_Tetrahedron_blend_out7_id3 - Blending_DF_Tetrahedron_blend_out1_id3*Blending_DF_Tetrahedron_blend_out3_id3*Blending_DF_Tetrahedron_blend_out8_id3 + Blending_DF_Tetrahedron_blend_out1_id3*Blending_DF_Tetrahedron_blend_out5_id3*Blending_DF_Tetrahedron_blend_out6_id3 + Blending_DF_Tetrahedron_blend_out2_id3*Blending_DF_Tetrahedron_blend_out3_id3*Blending_DF_Tetrahedron_blend_out7_id3 - Blending_DF_Tetrahedron_blend_out2_id3*Blending_DF_Tetrahedron_blend_out4_id3*Blending_DF_Tetrahedron_blend_out6_id3);
      real_t tmp_71 = 2*p_affine_2_0 + 2*tmp_0;
      real_t tmp_72 = Blending_DF_Tetrahedron_blend_out0_id0*tmp_71;
      real_t tmp_73 = 2*p_affine_2_1 + 2*tmp_3;
      real_t tmp_74 = Blending_DF_Tetrahedron_blend_out1_id0*tmp_73;
      real_t tmp_75 = 2*p_affine_2_2 + 2*tmp_6;
      real_t tmp_76 = Blending_DF_Tetrahedron_blend_out2_id0*tmp_75;
      real_t tmp_77 = -tmp_72 - tmp_74 - tmp_76;
      real_t tmp_78 = Blending_DF_Tetrahedron_blend_out3_id0*tmp_71;
      real_t tmp_79 = Blending_DF_Tetrahedron_blend_out4_id0*tmp_73;
      real_t tmp_80 = Blending_DF_Tetrahedron_blend_out5_id0*tmp_75;
      real_t tmp_81 = -tmp_78 - tmp_79 - tmp_80;
      real_t tmp_82 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_71;
      real_t tmp_83 = Blending_DF_Tetrahedron_blend_out7_id0*tmp_73;
      real_t tmp_84 = Blending_DF_Tetrahedron_blend_out8_id0*tmp_75;
      real_t tmp_85 = -tmp_82 - tmp_83 - tmp_84;
      real_t tmp_86 = Blending_DF_Tetrahedron_blend_out0_id1*tmp_71;
      real_t tmp_87 = Blending_DF_Tetrahedron_blend_out1_id1*tmp_73;
      real_t tmp_88 = Blending_DF_Tetrahedron_blend_out2_id1*tmp_75;
      real_t tmp_89 = -tmp_86 - tmp_87 - tmp_88;
      real_t tmp_90 = Blending_DF_Tetrahedron_blend_out3_id1*tmp_71;
      real_t tmp_91 = Blending_DF_Tetrahedron_blend_out4_id1*tmp_73;
      real_t tmp_92 = Blending_DF_Tetrahedron_blend_out5_id1*tmp_75;
      real_t tmp_93 = -tmp_90 - tmp_91 - tmp_92;
      real_t tmp_94 = Blending_DF_Tetrahedron_blend_out6_id1*tmp_71;
      real_t tmp_95 = Blending_DF_Tetrahedron_blend_out7_id1*tmp_73;
      real_t tmp_96 = Blending_DF_Tetrahedron_blend_out8_id1*tmp_75;
      real_t tmp_97 = -tmp_94 - tmp_95 - tmp_96;
      real_t tmp_98 = Blending_DF_Tetrahedron_blend_out0_id2*tmp_71;
      real_t tmp_99 = Blending_DF_Tetrahedron_blend_out1_id2*tmp_73;
      real_t tmp_100 = Blending_DF_Tetrahedron_blend_out2_id2*tmp_75;
      real_t tmp_101 = -tmp_100 - tmp_98 - tmp_99;
      real_t tmp_102 = Blending_DF_Tetrahedron_blend_out3_id2*tmp_71;
      real_t tmp_103 = Blending_DF_Tetrahedron_blend_out4_id2*tmp_73;
      real_t tmp_104 = Blending_DF_Tetrahedron_blend_out5_id2*tmp_75;
      real_t tmp_105 = -tmp_102 - tmp_103 - tmp_104;
      real_t tmp_106 = Blending_DF_Tetrahedron_blend_out6_id2*tmp_71;
      real_t tmp_107 = Blending_DF_Tetrahedron_blend_out7_id2*tmp_73;
      real_t tmp_108 = Blending_DF_Tetrahedron_blend_out8_id2*tmp_75;
      real_t tmp_109 = -tmp_106 - tmp_107 - tmp_108;
      real_t tmp_110 = Blending_DF_Tetrahedron_blend_out0_id3*tmp_71;
      real_t tmp_111 = Blending_DF_Tetrahedron_blend_out1_id3*tmp_73;
      real_t tmp_112 = Blending_DF_Tetrahedron_blend_out2_id3*tmp_75;
      real_t tmp_113 = -tmp_110 - tmp_111 - tmp_112;
      real_t tmp_114 = Blending_DF_Tetrahedron_blend_out3_id3*tmp_71;
      real_t tmp_115 = Blending_DF_Tetrahedron_blend_out4_id3*tmp_73;
      real_t tmp_116 = Blending_DF_Tetrahedron_blend_out5_id3*tmp_75;
      real_t tmp_117 = -tmp_114 - tmp_115 - tmp_116;
      real_t tmp_118 = Blending_DF_Tetrahedron_blend_out6_id3*tmp_71;
      real_t tmp_119 = Blending_DF_Tetrahedron_blend_out7_id3*tmp_73;
      real_t tmp_120 = Blending_DF_Tetrahedron_blend_out8_id3*tmp_75;
      real_t tmp_121 = -tmp_118 - tmp_119 - tmp_120;
      real_t tmp_122 = tmp_31*(tmp_13*tmp_81 + tmp_17*tmp_85 + tmp_77*tmp_9) + tmp_44*(tmp_35*tmp_89 + tmp_39*tmp_93 + tmp_43*tmp_97) + tmp_57*(tmp_101*tmp_48 + tmp_105*tmp_52 + tmp_109*tmp_56) + tmp_70*(tmp_113*tmp_61 + tmp_117*tmp_65 + tmp_121*tmp_69);
      real_t tmp_123 = 2*p_affine_3_0 + 2*tmp_0;
      real_t tmp_124 = Blending_DF_Tetrahedron_blend_out0_id0*tmp_123;
      real_t tmp_125 = 2*p_affine_3_1 + 2*tmp_3;
      real_t tmp_126 = Blending_DF_Tetrahedron_blend_out1_id0*tmp_125;
      real_t tmp_127 = 2*p_affine_3_2 + 2*tmp_6;
      real_t tmp_128 = Blending_DF_Tetrahedron_blend_out2_id0*tmp_127;
      real_t tmp_129 = tmp_124 + tmp_126 + tmp_128;
      real_t tmp_130 = Blending_DF_Tetrahedron_blend_out3_id0*tmp_123;
      real_t tmp_131 = Blending_DF_Tetrahedron_blend_out4_id0*tmp_125;
      real_t tmp_132 = Blending_DF_Tetrahedron_blend_out5_id0*tmp_127;
      real_t tmp_133 = tmp_130 + tmp_131 + tmp_132;
      real_t tmp_134 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_123;
      real_t tmp_135 = Blending_DF_Tetrahedron_blend_out7_id0*tmp_125;
      real_t tmp_136 = Blending_DF_Tetrahedron_blend_out8_id0*tmp_127;
      real_t tmp_137 = tmp_134 + tmp_135 + tmp_136;
      real_t tmp_138 = Blending_DF_Tetrahedron_blend_out0_id1*tmp_123;
      real_t tmp_139 = Blending_DF_Tetrahedron_blend_out1_id1*tmp_125;
      real_t tmp_140 = Blending_DF_Tetrahedron_blend_out2_id1*tmp_127;
      real_t tmp_141 = tmp_138 + tmp_139 + tmp_140;
      real_t tmp_142 = Blending_DF_Tetrahedron_blend_out3_id1*tmp_123;
      real_t tmp_143 = Blending_DF_Tetrahedron_blend_out4_id1*tmp_125;
      real_t tmp_144 = Blending_DF_Tetrahedron_blend_out5_id1*tmp_127;
      real_t tmp_145 = tmp_142 + tmp_143 + tmp_144;
      real_t tmp_146 = Blending_DF_Tetrahedron_blend_out6_id1*tmp_123;
      real_t tmp_147 = Blending_DF_Tetrahedron_blend_out7_id1*tmp_125;
      real_t tmp_148 = Blending_DF_Tetrahedron_blend_out8_id1*tmp_127;
      real_t tmp_149 = tmp_146 + tmp_147 + tmp_148;
      real_t tmp_150 = Blending_DF_Tetrahedron_blend_out0_id2*tmp_123;
      real_t tmp_151 = Blending_DF_Tetrahedron_blend_out1_id2*tmp_125;
      real_t tmp_152 = Blending_DF_Tetrahedron_blend_out2_id2*tmp_127;
      real_t tmp_153 = tmp_150 + tmp_151 + tmp_152;
      real_t tmp_154 = Blending_DF_Tetrahedron_blend_out3_id2*tmp_123;
      real_t tmp_155 = Blending_DF_Tetrahedron_blend_out4_id2*tmp_125;
      real_t tmp_156 = Blending_DF_Tetrahedron_blend_out5_id2*tmp_127;
      real_t tmp_157 = tmp_154 + tmp_155 + tmp_156;
      real_t tmp_158 = Blending_DF_Tetrahedron_blend_out6_id2*tmp_123;
      real_t tmp_159 = Blending_DF_Tetrahedron_blend_out7_id2*tmp_125;
      real_t tmp_160 = Blending_DF_Tetrahedron_blend_out8_id2*tmp_127;
      real_t tmp_161 = tmp_158 + tmp_159 + tmp_160;
      real_t tmp_162 = Blending_DF_Tetrahedron_blend_out0_id3*tmp_123;
      real_t tmp_163 = Blending_DF_Tetrahedron_blend_out1_id3*tmp_125;
      real_t tmp_164 = Blending_DF_Tetrahedron_blend_out2_id3*tmp_127;
      real_t tmp_165 = tmp_162 + tmp_163 + tmp_164;
      real_t tmp_166 = Blending_DF_Tetrahedron_blend_out3_id3*tmp_123;
      real_t tmp_167 = Blending_DF_Tetrahedron_blend_out4_id3*tmp_125;
      real_t tmp_168 = Blending_DF_Tetrahedron_blend_out5_id3*tmp_127;
      real_t tmp_169 = tmp_166 + tmp_167 + tmp_168;
      real_t tmp_170 = Blending_DF_Tetrahedron_blend_out6_id3*tmp_123;
      real_t tmp_171 = Blending_DF_Tetrahedron_blend_out7_id3*tmp_125;
      real_t tmp_172 = Blending_DF_Tetrahedron_blend_out8_id3*tmp_127;
      real_t tmp_173 = tmp_170 + tmp_171 + tmp_172;
      real_t tmp_174 = tmp_31*(tmp_129*tmp_9 + tmp_13*tmp_133 + tmp_137*tmp_17) + tmp_44*(tmp_141*tmp_35 + tmp_145*tmp_39 + tmp_149*tmp_43) + tmp_57*(tmp_153*tmp_48 + tmp_157*tmp_52 + tmp_161*tmp_56) + tmp_70*(tmp_165*tmp_61 + tmp_169*tmp_65 + tmp_173*tmp_69);
      real_t tmp_175 = -tmp_2 - tmp_5 + tmp_72 + tmp_74 + tmp_76 - tmp_8;
      real_t tmp_176 = -tmp_10 - tmp_11 - tmp_12 + tmp_78 + tmp_79 + tmp_80;
      real_t tmp_177 = -tmp_14 - tmp_15 - tmp_16 + tmp_82 + tmp_83 + tmp_84;
      real_t tmp_178 = -tmp_32 - tmp_33 - tmp_34 + tmp_86 + tmp_87 + tmp_88;
      real_t tmp_179 = -tmp_36 - tmp_37 - tmp_38 + tmp_90 + tmp_91 + tmp_92;
      real_t tmp_180 = -tmp_40 - tmp_41 - tmp_42 + tmp_94 + tmp_95 + tmp_96;
      real_t tmp_181 = tmp_100 - tmp_45 - tmp_46 - tmp_47 + tmp_98 + tmp_99;
      real_t tmp_182 = tmp_102 + tmp_103 + tmp_104 - tmp_49 - tmp_50 - tmp_51;
      real_t tmp_183 = tmp_106 + tmp_107 + tmp_108 - tmp_53 - tmp_54 - tmp_55;
      real_t tmp_184 = tmp_110 + tmp_111 + tmp_112 - tmp_58 - tmp_59 - tmp_60;
      real_t tmp_185 = tmp_114 + tmp_115 + tmp_116 - tmp_62 - tmp_63 - tmp_64;
      real_t tmp_186 = tmp_118 + tmp_119 + tmp_120 - tmp_66 - tmp_67 - tmp_68;
      real_t tmp_187 = tmp_31*(tmp_13*tmp_176 + tmp_17*tmp_177 + tmp_175*tmp_9) + tmp_44*(tmp_178*tmp_35 + tmp_179*tmp_39 + tmp_180*tmp_43) + tmp_57*(tmp_181*tmp_48 + tmp_182*tmp_52 + tmp_183*tmp_56) + tmp_70*(tmp_184*tmp_61 + tmp_185*tmp_65 + tmp_186*tmp_69);
      real_t tmp_188 = -tmp_124 - tmp_126 - tmp_128 + tmp_9;
      real_t tmp_189 = tmp_13 - tmp_130 - tmp_131 - tmp_132;
      real_t tmp_190 = -tmp_134 - tmp_135 - tmp_136 + tmp_17;
      real_t tmp_191 = -tmp_138 - tmp_139 - tmp_140 + tmp_35;
      real_t tmp_192 = -tmp_142 - tmp_143 - tmp_144 + tmp_39;
      real_t tmp_193 = -tmp_146 - tmp_147 - tmp_148 + tmp_43;
      real_t tmp_194 = -tmp_150 - tmp_151 - tmp_152 + tmp_48;
      real_t tmp_195 = -tmp_154 - tmp_155 - tmp_156 + tmp_52;
      real_t tmp_196 = -tmp_158 - tmp_159 - tmp_160 + tmp_56;
      real_t tmp_197 = -tmp_162 - tmp_163 - tmp_164 + tmp_61;
      real_t tmp_198 = -tmp_166 - tmp_167 - tmp_168 + tmp_65;
      real_t tmp_199 = -tmp_170 - tmp_171 - tmp_172 + tmp_69;
      real_t tmp_200 = tmp_31*(tmp_13*tmp_189 + tmp_17*tmp_190 + tmp_188*tmp_9) + tmp_44*(tmp_191*tmp_35 + tmp_192*tmp_39 + tmp_193*tmp_43) + tmp_57*(tmp_194*tmp_48 + tmp_195*tmp_52 + tmp_196*tmp_56) + tmp_70*(tmp_197*tmp_61 + tmp_198*tmp_65 + tmp_199*tmp_69);
      real_t tmp_201 = tmp_129 + tmp_77;
      real_t tmp_202 = tmp_133 + tmp_81;
      real_t tmp_203 = tmp_137 + tmp_85;
      real_t tmp_204 = tmp_141 + tmp_89;
      real_t tmp_205 = tmp_145 + tmp_93;
      real_t tmp_206 = tmp_149 + tmp_97;
      real_t tmp_207 = tmp_101 + tmp_153;
      real_t tmp_208 = tmp_105 + tmp_157;
      real_t tmp_209 = tmp_109 + tmp_161;
      real_t tmp_210 = tmp_113 + tmp_165;
      real_t tmp_211 = tmp_117 + tmp_169;
      real_t tmp_212 = tmp_121 + tmp_173;
      real_t tmp_213 = tmp_31*(tmp_13*tmp_202 + tmp_17*tmp_203 + tmp_201*tmp_9) + tmp_44*(tmp_204*tmp_35 + tmp_205*tmp_39 + tmp_206*tmp_43) + tmp_57*(tmp_207*tmp_48 + tmp_208*tmp_52 + tmp_209*tmp_56) + tmp_70*(tmp_210*tmp_61 + tmp_211*tmp_65 + tmp_212*tmp_69);
      real_t tmp_214 = tmp_31*(tmp_129*tmp_77 + tmp_133*tmp_81 + tmp_137*tmp_85) + tmp_44*(tmp_141*tmp_89 + tmp_145*tmp_93 + tmp_149*tmp_97) + tmp_57*(tmp_101*tmp_153 + tmp_105*tmp_157 + tmp_109*tmp_161) + tmp_70*(tmp_113*tmp_165 + tmp_117*tmp_169 + tmp_121*tmp_173);
      real_t tmp_215 = tmp_31*(tmp_175*tmp_77 + tmp_176*tmp_81 + tmp_177*tmp_85) + tmp_44*(tmp_178*tmp_89 + tmp_179*tmp_93 + tmp_180*tmp_97) + tmp_57*(tmp_101*tmp_181 + tmp_105*tmp_182 + tmp_109*tmp_183) + tmp_70*(tmp_113*tmp_184 + tmp_117*tmp_185 + tmp_121*tmp_186);
      real_t tmp_216 = tmp_31*(tmp_188*tmp_77 + tmp_189*tmp_81 + tmp_190*tmp_85) + tmp_44*(tmp_191*tmp_89 + tmp_192*tmp_93 + tmp_193*tmp_97) + tmp_57*(tmp_101*tmp_194 + tmp_105*tmp_195 + tmp_109*tmp_196) + tmp_70*(tmp_113*tmp_197 + tmp_117*tmp_198 + tmp_121*tmp_199);
      real_t tmp_217 = tmp_31*(tmp_201*tmp_77 + tmp_202*tmp_81 + tmp_203*tmp_85) + tmp_44*(tmp_204*tmp_89 + tmp_205*tmp_93 + tmp_206*tmp_97) + tmp_57*(tmp_101*tmp_207 + tmp_105*tmp_208 + tmp_109*tmp_209) + tmp_70*(tmp_113*tmp_210 + tmp_117*tmp_211 + tmp_121*tmp_212);
      real_t tmp_218 = tmp_31*(tmp_129*tmp_175 + tmp_133*tmp_176 + tmp_137*tmp_177) + tmp_44*(tmp_141*tmp_178 + tmp_145*tmp_179 + tmp_149*tmp_180) + tmp_57*(tmp_153*tmp_181 + tmp_157*tmp_182 + tmp_161*tmp_183) + tmp_70*(tmp_165*tmp_184 + tmp_169*tmp_185 + tmp_173*tmp_186);
      real_t tmp_219 = tmp_31*(tmp_129*tmp_188 + tmp_133*tmp_189 + tmp_137*tmp_190) + tmp_44*(tmp_141*tmp_191 + tmp_145*tmp_192 + tmp_149*tmp_193) + tmp_57*(tmp_153*tmp_194 + tmp_157*tmp_195 + tmp_161*tmp_196) + tmp_70*(tmp_165*tmp_197 + tmp_169*tmp_198 + tmp_173*tmp_199);
      real_t tmp_220 = tmp_31*(tmp_129*tmp_201 + tmp_133*tmp_202 + tmp_137*tmp_203) + tmp_44*(tmp_141*tmp_204 + tmp_145*tmp_205 + tmp_149*tmp_206) + tmp_57*(tmp_153*tmp_207 + tmp_157*tmp_208 + tmp_161*tmp_209) + tmp_70*(tmp_165*tmp_210 + tmp_169*tmp_211 + tmp_173*tmp_212);
      real_t tmp_221 = tmp_31*(tmp_175*tmp_188 + tmp_176*tmp_189 + tmp_177*tmp_190) + tmp_44*(tmp_178*tmp_191 + tmp_179*tmp_192 + tmp_180*tmp_193) + tmp_57*(tmp_181*tmp_194 + tmp_182*tmp_195 + tmp_183*tmp_196) + tmp_70*(tmp_184*tmp_197 + tmp_185*tmp_198 + tmp_186*tmp_199);
      real_t tmp_222 = tmp_31*(tmp_175*tmp_201 + tmp_176*tmp_202 + tmp_177*tmp_203) + tmp_44*(tmp_178*tmp_204 + tmp_179*tmp_205 + tmp_180*tmp_206) + tmp_57*(tmp_181*tmp_207 + tmp_182*tmp_208 + tmp_183*tmp_209) + tmp_70*(tmp_184*tmp_210 + tmp_185*tmp_211 + tmp_186*tmp_212);
      real_t tmp_223 = tmp_31*(tmp_188*tmp_201 + tmp_189*tmp_202 + tmp_190*tmp_203) + tmp_44*(tmp_191*tmp_204 + tmp_192*tmp_205 + tmp_193*tmp_206) + tmp_57*(tmp_194*tmp_207 + tmp_195*tmp_208 + tmp_196*tmp_209) + tmp_70*(tmp_197*tmp_210 + tmp_198*tmp_211 + tmp_199*tmp_212);
      real_t a_0_0 = tmp_31*((tmp_13*tmp_13) + (tmp_17*tmp_17) + (tmp_9*tmp_9)) + tmp_44*((tmp_35*tmp_35) + (tmp_39*tmp_39) + (tmp_43*tmp_43)) + tmp_57*((tmp_48*tmp_48) + (tmp_52*tmp_52) + (tmp_56*tmp_56)) + tmp_70*((tmp_61*tmp_61) + (tmp_65*tmp_65) + (tmp_69*tmp_69));
      real_t a_0_1 = tmp_122;
      real_t a_0_2 = tmp_174;
      real_t a_0_3 = tmp_187;
      real_t a_0_4 = tmp_200;
      real_t a_0_5 = tmp_213;
      real_t a_1_0 = tmp_122;
      real_t a_1_1 = tmp_31*((tmp_77*tmp_77) + (tmp_81*tmp_81) + (tmp_85*tmp_85)) + tmp_44*((tmp_89*tmp_89) + (tmp_93*tmp_93) + (tmp_97*tmp_97)) + tmp_57*((tmp_101*tmp_101) + (tmp_105*tmp_105) + (tmp_109*tmp_109)) + tmp_70*((tmp_113*tmp_113) + (tmp_117*tmp_117) + (tmp_121*tmp_121));
      real_t a_1_2 = tmp_214;
      real_t a_1_3 = tmp_215;
      real_t a_1_4 = tmp_216;
      real_t a_1_5 = tmp_217;
      real_t a_2_0 = tmp_174;
      real_t a_2_1 = tmp_214;
      real_t a_2_2 = tmp_31*((tmp_129*tmp_129) + (tmp_133*tmp_133) + (tmp_137*tmp_137)) + tmp_44*((tmp_141*tmp_141) + (tmp_145*tmp_145) + (tmp_149*tmp_149)) + tmp_57*((tmp_153*tmp_153) + (tmp_157*tmp_157) + (tmp_161*tmp_161)) + tmp_70*((tmp_165*tmp_165) + (tmp_169*tmp_169) + (tmp_173*tmp_173));
      real_t a_2_3 = tmp_218;
      real_t a_2_4 = tmp_219;
      real_t a_2_5 = tmp_220;
      real_t a_3_0 = tmp_187;
      real_t a_3_1 = tmp_215;
      real_t a_3_2 = tmp_218;
      real_t a_3_3 = tmp_31*((tmp_175*tmp_175) + (tmp_176*tmp_176) + (tmp_177*tmp_177)) + tmp_44*((tmp_178*tmp_178) + (tmp_179*tmp_179) + (tmp_180*tmp_180)) + tmp_57*((tmp_181*tmp_181) + (tmp_182*tmp_182) + (tmp_183*tmp_183)) + tmp_70*((tmp_184*tmp_184) + (tmp_185*tmp_185) + (tmp_186*tmp_186));
      real_t a_3_4 = tmp_221;
      real_t a_3_5 = tmp_222;
      real_t a_4_0 = tmp_200;
      real_t a_4_1 = tmp_216;
      real_t a_4_2 = tmp_219;
      real_t a_4_3 = tmp_221;
      real_t a_4_4 = tmp_31*((tmp_188*tmp_188) + (tmp_189*tmp_189) + (tmp_190*tmp_190)) + tmp_44*((tmp_191*tmp_191) + (tmp_192*tmp_192) + (tmp_193*tmp_193)) + tmp_57*((tmp_194*tmp_194) + (tmp_195*tmp_195) + (tmp_196*tmp_196)) + tmp_70*((tmp_197*tmp_197) + (tmp_198*tmp_198) + (tmp_199*tmp_199));
      real_t a_4_5 = tmp_223;
      real_t a_5_0 = tmp_213;
      real_t a_5_1 = tmp_217;
      real_t a_5_2 = tmp_220;
      real_t a_5_3 = tmp_222;
      real_t a_5_4 = tmp_223;
      real_t a_5_5 = tmp_31*((tmp_201*tmp_201) + (tmp_202*tmp_202) + (tmp_203*tmp_203)) + tmp_44*((tmp_204*tmp_204) + (tmp_205*tmp_205) + (tmp_206*tmp_206)) + tmp_57*((tmp_207*tmp_207) + (tmp_208*tmp_208) + (tmp_209*tmp_209)) + tmp_70*((tmp_210*tmp_210) + (tmp_211*tmp_211) + (tmp_212*tmp_212));
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

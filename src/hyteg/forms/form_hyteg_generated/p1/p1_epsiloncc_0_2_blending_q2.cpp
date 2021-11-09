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

#include "p1_epsiloncc_0_2_blending_q2.hpp"

namespace hyteg {
namespace forms {

   void p1_epsiloncc_0_2_blending_q2::integrateAll( const std::array< Point3D, 3 >& , Matrix< real_t, 3, 3 >&  ) const
   {
      
   }

   void p1_epsiloncc_0_2_blending_q2::integrateRow0( const std::array< Point3D, 3 >& , Matrix< real_t, 1, 3 >&  ) const
   {
      
   }

   void p1_epsiloncc_0_2_blending_q2::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      Blending_DF_Tetrahedron( 0.13819660112501042*p_affine_0_0 + 0.13819660112501059*p_affine_1_0 + 0.13819660112501059*p_affine_2_0 + 0.58541019662496829*p_affine_3_0, 0.13819660112501042*p_affine_0_1 + 0.13819660112501059*p_affine_1_1 + 0.13819660112501059*p_affine_2_1 + 0.58541019662496829*p_affine_3_1, 0.13819660112501042*p_affine_0_2 + 0.13819660112501059*p_affine_1_2 + 0.13819660112501059*p_affine_2_2 + 0.58541019662496829*p_affine_3_2, &Blending_DF_Tetrahedron_0_0, &Blending_DF_Tetrahedron_0_1, &Blending_DF_Tetrahedron_0_2, &Blending_DF_Tetrahedron_0_3, &Blending_DF_Tetrahedron_0_4, &Blending_DF_Tetrahedron_0_5, &Blending_DF_Tetrahedron_0_6, &Blending_DF_Tetrahedron_0_7, &Blending_DF_Tetrahedron_0_8 );
      Blending_DF_Tetrahedron( 0.13819660112501048*p_affine_0_0 + 0.13819660112501059*p_affine_1_0 + 0.58541019662496829*p_affine_2_0 + 0.13819660112501059*p_affine_3_0, 0.13819660112501048*p_affine_0_1 + 0.13819660112501059*p_affine_1_1 + 0.58541019662496829*p_affine_2_1 + 0.13819660112501059*p_affine_3_1, 0.13819660112501048*p_affine_0_2 + 0.13819660112501059*p_affine_1_2 + 0.58541019662496829*p_affine_2_2 + 0.13819660112501059*p_affine_3_2, &Blending_DF_Tetrahedron_1_0, &Blending_DF_Tetrahedron_1_1, &Blending_DF_Tetrahedron_1_2, &Blending_DF_Tetrahedron_1_3, &Blending_DF_Tetrahedron_1_4, &Blending_DF_Tetrahedron_1_5, &Blending_DF_Tetrahedron_1_6, &Blending_DF_Tetrahedron_1_7, &Blending_DF_Tetrahedron_1_8 );
      Blending_DF_Tetrahedron( 0.13819660112501053*p_affine_0_0 + 0.58541019662496829*p_affine_1_0 + 0.13819660112501059*p_affine_2_0 + 0.13819660112501059*p_affine_3_0, 0.13819660112501053*p_affine_0_1 + 0.58541019662496829*p_affine_1_1 + 0.13819660112501059*p_affine_2_1 + 0.13819660112501059*p_affine_3_1, 0.13819660112501053*p_affine_0_2 + 0.58541019662496829*p_affine_1_2 + 0.13819660112501059*p_affine_2_2 + 0.13819660112501059*p_affine_3_2, &Blending_DF_Tetrahedron_2_0, &Blending_DF_Tetrahedron_2_1, &Blending_DF_Tetrahedron_2_2, &Blending_DF_Tetrahedron_2_3, &Blending_DF_Tetrahedron_2_4, &Blending_DF_Tetrahedron_2_5, &Blending_DF_Tetrahedron_2_6, &Blending_DF_Tetrahedron_2_7, &Blending_DF_Tetrahedron_2_8 );
      Blending_DF_Tetrahedron( 0.58541019662496807*p_affine_0_0 + 0.13819660112501059*p_affine_1_0 + 0.13819660112501059*p_affine_2_0 + 0.13819660112501059*p_affine_3_0, 0.58541019662496807*p_affine_0_1 + 0.13819660112501059*p_affine_1_1 + 0.13819660112501059*p_affine_2_1 + 0.13819660112501059*p_affine_3_1, 0.58541019662496807*p_affine_0_2 + 0.13819660112501059*p_affine_1_2 + 0.13819660112501059*p_affine_2_2 + 0.13819660112501059*p_affine_3_2, &Blending_DF_Tetrahedron_3_0, &Blending_DF_Tetrahedron_3_1, &Blending_DF_Tetrahedron_3_2, &Blending_DF_Tetrahedron_3_3, &Blending_DF_Tetrahedron_3_4, &Blending_DF_Tetrahedron_3_5, &Blending_DF_Tetrahedron_3_6, &Blending_DF_Tetrahedron_3_7, &Blending_DF_Tetrahedron_3_8 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = Blending_DF_Tetrahedron_0_0*Blending_DF_Tetrahedron_0_4;
      real_t tmp_10 = Blending_DF_Tetrahedron_0_1*Blending_DF_Tetrahedron_0_5;
      real_t tmp_11 = Blending_DF_Tetrahedron_0_2*Blending_DF_Tetrahedron_0_3;
      real_t tmp_12 = Blending_DF_Tetrahedron_0_0*Blending_DF_Tetrahedron_0_5;
      real_t tmp_13 = Blending_DF_Tetrahedron_0_1*Blending_DF_Tetrahedron_0_3;
      real_t tmp_14 = Blending_DF_Tetrahedron_0_2*Blending_DF_Tetrahedron_0_4;
      real_t tmp_15 = Blending_DF_Tetrahedron_0_6*tmp_10 - Blending_DF_Tetrahedron_0_6*tmp_14 + Blending_DF_Tetrahedron_0_7*tmp_11 - Blending_DF_Tetrahedron_0_7*tmp_12 - Blending_DF_Tetrahedron_0_8*tmp_13 + Blending_DF_Tetrahedron_0_8*tmp_9;
      real_t tmp_16 = -p_affine_0_2;
      real_t tmp_17 = p_affine_3_2 + tmp_16;
      real_t tmp_18 = p_affine_1_2 + tmp_16;
      real_t tmp_19 = p_affine_3_1 + tmp_2;
      real_t tmp_20 = tmp_19*tmp_5;
      real_t tmp_21 = p_affine_2_2 + tmp_16;
      real_t tmp_22 = p_affine_3_0 + tmp_0;
      real_t tmp_23 = tmp_22*tmp_6;
      real_t tmp_24 = tmp_1*tmp_19;
      real_t tmp_25 = tmp_22*tmp_3;
      real_t tmp_26 = 1.0 / (tmp_17*tmp_4 - tmp_17*tmp_7 + tmp_18*tmp_20 - tmp_18*tmp_25 + tmp_21*tmp_23 - tmp_21*tmp_24);
      real_t tmp_27 = tmp_26/tmp_15;
      real_t tmp_28 = tmp_27*(Blending_DF_Tetrahedron_0_3*Blending_DF_Tetrahedron_0_7 - Blending_DF_Tetrahedron_0_4*Blending_DF_Tetrahedron_0_6);
      real_t tmp_29 = tmp_28*tmp_8;
      real_t tmp_30 = tmp_23 - tmp_24;
      real_t tmp_31 = tmp_28*tmp_30;
      real_t tmp_32 = tmp_20 - tmp_25;
      real_t tmp_33 = tmp_28*tmp_32;
      real_t tmp_34 = -tmp_1*tmp_21 + tmp_18*tmp_5;
      real_t tmp_35 = tmp_27*(-Blending_DF_Tetrahedron_0_3*Blending_DF_Tetrahedron_0_8 + Blending_DF_Tetrahedron_0_5*Blending_DF_Tetrahedron_0_6);
      real_t tmp_36 = tmp_34*tmp_35;
      real_t tmp_37 = tmp_1*tmp_17 - tmp_18*tmp_22;
      real_t tmp_38 = tmp_35*tmp_37;
      real_t tmp_39 = -tmp_17*tmp_5 + tmp_21*tmp_22;
      real_t tmp_40 = tmp_35*tmp_39;
      real_t tmp_41 = -tmp_18*tmp_3 + tmp_21*tmp_6;
      real_t tmp_42 = tmp_27*(Blending_DF_Tetrahedron_0_4*Blending_DF_Tetrahedron_0_8 - Blending_DF_Tetrahedron_0_5*Blending_DF_Tetrahedron_0_7);
      real_t tmp_43 = tmp_41*tmp_42;
      real_t tmp_44 = -tmp_17*tmp_6 + tmp_18*tmp_19;
      real_t tmp_45 = tmp_42*tmp_44;
      real_t tmp_46 = tmp_17*tmp_3 - tmp_19*tmp_21;
      real_t tmp_47 = tmp_42*tmp_46;
      real_t tmp_48 = -1.0*tmp_29 - 1.0*tmp_31 - 1.0*tmp_33 - 1.0*tmp_36 - 1.0*tmp_38 - 1.0*tmp_40 - 1.0*tmp_43 - 1.0*tmp_45 - 1.0*tmp_47;
      real_t tmp_49 = 0.5*tmp_27;
      real_t tmp_50 = tmp_49*(-tmp_13 + tmp_9);
      real_t tmp_51 = tmp_50*tmp_8;
      real_t tmp_52 = tmp_30*tmp_50;
      real_t tmp_53 = tmp_32*tmp_50;
      real_t tmp_54 = tmp_49*(tmp_11 - tmp_12);
      real_t tmp_55 = tmp_34*tmp_54;
      real_t tmp_56 = tmp_37*tmp_54;
      real_t tmp_57 = tmp_39*tmp_54;
      real_t tmp_58 = tmp_49*(tmp_10 - tmp_14);
      real_t tmp_59 = tmp_41*tmp_58;
      real_t tmp_60 = tmp_44*tmp_58;
      real_t tmp_61 = tmp_46*tmp_58;
      real_t tmp_62 = p_affine_0_0*p_affine_1_1;
      real_t tmp_63 = p_affine_0_0*p_affine_1_2;
      real_t tmp_64 = p_affine_2_1*p_affine_3_2;
      real_t tmp_65 = p_affine_0_1*p_affine_1_0;
      real_t tmp_66 = p_affine_0_1*p_affine_1_2;
      real_t tmp_67 = p_affine_2_2*p_affine_3_0;
      real_t tmp_68 = p_affine_0_2*p_affine_1_0;
      real_t tmp_69 = p_affine_0_2*p_affine_1_1;
      real_t tmp_70 = p_affine_2_0*p_affine_3_1;
      real_t tmp_71 = p_affine_2_2*p_affine_3_1;
      real_t tmp_72 = p_affine_2_0*p_affine_3_2;
      real_t tmp_73 = p_affine_2_1*p_affine_3_0;
      real_t tmp_74 = 2*std::abs(p_affine_0_0*tmp_64 - p_affine_0_0*tmp_71 + p_affine_0_1*tmp_67 - p_affine_0_1*tmp_72 + p_affine_0_2*tmp_70 - p_affine_0_2*tmp_73 - p_affine_1_0*tmp_64 + p_affine_1_0*tmp_71 - p_affine_1_1*tmp_67 + p_affine_1_1*tmp_72 - p_affine_1_2*tmp_70 + p_affine_1_2*tmp_73 + p_affine_2_0*tmp_66 - p_affine_2_0*tmp_69 - p_affine_2_1*tmp_63 + p_affine_2_1*tmp_68 + p_affine_2_2*tmp_62 - p_affine_2_2*tmp_65 - p_affine_3_0*tmp_66 + p_affine_3_0*tmp_69 + p_affine_3_1*tmp_63 - p_affine_3_1*tmp_68 - p_affine_3_2*tmp_62 + p_affine_3_2*tmp_65);
      real_t tmp_75 = 0.041666666666666657*tmp_74*std::abs(tmp_15);
      real_t tmp_76 = tmp_75*(-tmp_51 - tmp_52 - tmp_53 - tmp_55 - tmp_56 - tmp_57 - tmp_59 - tmp_60 - tmp_61);
      real_t tmp_77 = Blending_DF_Tetrahedron_1_0*Blending_DF_Tetrahedron_1_4;
      real_t tmp_78 = Blending_DF_Tetrahedron_1_1*Blending_DF_Tetrahedron_1_5;
      real_t tmp_79 = Blending_DF_Tetrahedron_1_2*Blending_DF_Tetrahedron_1_3;
      real_t tmp_80 = Blending_DF_Tetrahedron_1_0*Blending_DF_Tetrahedron_1_5;
      real_t tmp_81 = Blending_DF_Tetrahedron_1_1*Blending_DF_Tetrahedron_1_3;
      real_t tmp_82 = Blending_DF_Tetrahedron_1_2*Blending_DF_Tetrahedron_1_4;
      real_t tmp_83 = Blending_DF_Tetrahedron_1_6*tmp_78 - Blending_DF_Tetrahedron_1_6*tmp_82 + Blending_DF_Tetrahedron_1_7*tmp_79 - Blending_DF_Tetrahedron_1_7*tmp_80 + Blending_DF_Tetrahedron_1_8*tmp_77 - Blending_DF_Tetrahedron_1_8*tmp_81;
      real_t tmp_84 = tmp_26/tmp_83;
      real_t tmp_85 = tmp_84*(Blending_DF_Tetrahedron_1_3*Blending_DF_Tetrahedron_1_7 - Blending_DF_Tetrahedron_1_4*Blending_DF_Tetrahedron_1_6);
      real_t tmp_86 = tmp_8*tmp_85;
      real_t tmp_87 = tmp_30*tmp_85;
      real_t tmp_88 = tmp_32*tmp_85;
      real_t tmp_89 = tmp_84*(-Blending_DF_Tetrahedron_1_3*Blending_DF_Tetrahedron_1_8 + Blending_DF_Tetrahedron_1_5*Blending_DF_Tetrahedron_1_6);
      real_t tmp_90 = tmp_34*tmp_89;
      real_t tmp_91 = tmp_37*tmp_89;
      real_t tmp_92 = tmp_39*tmp_89;
      real_t tmp_93 = tmp_84*(Blending_DF_Tetrahedron_1_4*Blending_DF_Tetrahedron_1_8 - Blending_DF_Tetrahedron_1_5*Blending_DF_Tetrahedron_1_7);
      real_t tmp_94 = tmp_41*tmp_93;
      real_t tmp_95 = tmp_44*tmp_93;
      real_t tmp_96 = tmp_46*tmp_93;
      real_t tmp_97 = -1.0*tmp_86 - 1.0*tmp_87 - 1.0*tmp_88 - 1.0*tmp_90 - 1.0*tmp_91 - 1.0*tmp_92 - 1.0*tmp_94 - 1.0*tmp_95 - 1.0*tmp_96;
      real_t tmp_98 = 0.5*tmp_84;
      real_t tmp_99 = tmp_98*(tmp_77 - tmp_81);
      real_t tmp_100 = tmp_8*tmp_99;
      real_t tmp_101 = tmp_30*tmp_99;
      real_t tmp_102 = tmp_32*tmp_99;
      real_t tmp_103 = tmp_98*(tmp_79 - tmp_80);
      real_t tmp_104 = tmp_103*tmp_34;
      real_t tmp_105 = tmp_103*tmp_37;
      real_t tmp_106 = tmp_103*tmp_39;
      real_t tmp_107 = tmp_98*(tmp_78 - tmp_82);
      real_t tmp_108 = tmp_107*tmp_41;
      real_t tmp_109 = tmp_107*tmp_44;
      real_t tmp_110 = tmp_107*tmp_46;
      real_t tmp_111 = 0.041666666666666657*tmp_74*std::abs(tmp_83);
      real_t tmp_112 = tmp_111*(-tmp_100 - tmp_101 - tmp_102 - tmp_104 - tmp_105 - tmp_106 - tmp_108 - tmp_109 - tmp_110);
      real_t tmp_113 = Blending_DF_Tetrahedron_2_0*Blending_DF_Tetrahedron_2_4;
      real_t tmp_114 = Blending_DF_Tetrahedron_2_1*Blending_DF_Tetrahedron_2_5;
      real_t tmp_115 = Blending_DF_Tetrahedron_2_2*Blending_DF_Tetrahedron_2_3;
      real_t tmp_116 = Blending_DF_Tetrahedron_2_0*Blending_DF_Tetrahedron_2_5;
      real_t tmp_117 = Blending_DF_Tetrahedron_2_1*Blending_DF_Tetrahedron_2_3;
      real_t tmp_118 = Blending_DF_Tetrahedron_2_2*Blending_DF_Tetrahedron_2_4;
      real_t tmp_119 = Blending_DF_Tetrahedron_2_6*tmp_114 - Blending_DF_Tetrahedron_2_6*tmp_118 + Blending_DF_Tetrahedron_2_7*tmp_115 - Blending_DF_Tetrahedron_2_7*tmp_116 + Blending_DF_Tetrahedron_2_8*tmp_113 - Blending_DF_Tetrahedron_2_8*tmp_117;
      real_t tmp_120 = tmp_26/tmp_119;
      real_t tmp_121 = tmp_120*(Blending_DF_Tetrahedron_2_3*Blending_DF_Tetrahedron_2_7 - Blending_DF_Tetrahedron_2_4*Blending_DF_Tetrahedron_2_6);
      real_t tmp_122 = tmp_121*tmp_8;
      real_t tmp_123 = tmp_121*tmp_30;
      real_t tmp_124 = tmp_121*tmp_32;
      real_t tmp_125 = tmp_120*(-Blending_DF_Tetrahedron_2_3*Blending_DF_Tetrahedron_2_8 + Blending_DF_Tetrahedron_2_5*Blending_DF_Tetrahedron_2_6);
      real_t tmp_126 = tmp_125*tmp_34;
      real_t tmp_127 = tmp_125*tmp_37;
      real_t tmp_128 = tmp_125*tmp_39;
      real_t tmp_129 = tmp_120*(Blending_DF_Tetrahedron_2_4*Blending_DF_Tetrahedron_2_8 - Blending_DF_Tetrahedron_2_5*Blending_DF_Tetrahedron_2_7);
      real_t tmp_130 = tmp_129*tmp_41;
      real_t tmp_131 = tmp_129*tmp_44;
      real_t tmp_132 = tmp_129*tmp_46;
      real_t tmp_133 = -1.0*tmp_122 - 1.0*tmp_123 - 1.0*tmp_124 - 1.0*tmp_126 - 1.0*tmp_127 - 1.0*tmp_128 - 1.0*tmp_130 - 1.0*tmp_131 - 1.0*tmp_132;
      real_t tmp_134 = 0.5*tmp_120;
      real_t tmp_135 = tmp_134*(tmp_113 - tmp_117);
      real_t tmp_136 = tmp_135*tmp_8;
      real_t tmp_137 = tmp_135*tmp_30;
      real_t tmp_138 = tmp_135*tmp_32;
      real_t tmp_139 = tmp_134*(tmp_115 - tmp_116);
      real_t tmp_140 = tmp_139*tmp_34;
      real_t tmp_141 = tmp_139*tmp_37;
      real_t tmp_142 = tmp_139*tmp_39;
      real_t tmp_143 = tmp_134*(tmp_114 - tmp_118);
      real_t tmp_144 = tmp_143*tmp_41;
      real_t tmp_145 = tmp_143*tmp_44;
      real_t tmp_146 = tmp_143*tmp_46;
      real_t tmp_147 = 0.041666666666666657*tmp_74*std::abs(tmp_119);
      real_t tmp_148 = tmp_147*(-tmp_136 - tmp_137 - tmp_138 - tmp_140 - tmp_141 - tmp_142 - tmp_144 - tmp_145 - tmp_146);
      real_t tmp_149 = Blending_DF_Tetrahedron_3_0*Blending_DF_Tetrahedron_3_4;
      real_t tmp_150 = Blending_DF_Tetrahedron_3_1*Blending_DF_Tetrahedron_3_5;
      real_t tmp_151 = Blending_DF_Tetrahedron_3_2*Blending_DF_Tetrahedron_3_3;
      real_t tmp_152 = Blending_DF_Tetrahedron_3_0*Blending_DF_Tetrahedron_3_5;
      real_t tmp_153 = Blending_DF_Tetrahedron_3_1*Blending_DF_Tetrahedron_3_3;
      real_t tmp_154 = Blending_DF_Tetrahedron_3_2*Blending_DF_Tetrahedron_3_4;
      real_t tmp_155 = Blending_DF_Tetrahedron_3_6*tmp_150 - Blending_DF_Tetrahedron_3_6*tmp_154 + Blending_DF_Tetrahedron_3_7*tmp_151 - Blending_DF_Tetrahedron_3_7*tmp_152 + Blending_DF_Tetrahedron_3_8*tmp_149 - Blending_DF_Tetrahedron_3_8*tmp_153;
      real_t tmp_156 = tmp_26/tmp_155;
      real_t tmp_157 = tmp_156*(Blending_DF_Tetrahedron_3_3*Blending_DF_Tetrahedron_3_7 - Blending_DF_Tetrahedron_3_4*Blending_DF_Tetrahedron_3_6);
      real_t tmp_158 = tmp_157*tmp_8;
      real_t tmp_159 = tmp_157*tmp_30;
      real_t tmp_160 = tmp_157*tmp_32;
      real_t tmp_161 = tmp_156*(-Blending_DF_Tetrahedron_3_3*Blending_DF_Tetrahedron_3_8 + Blending_DF_Tetrahedron_3_5*Blending_DF_Tetrahedron_3_6);
      real_t tmp_162 = tmp_161*tmp_34;
      real_t tmp_163 = tmp_161*tmp_37;
      real_t tmp_164 = tmp_161*tmp_39;
      real_t tmp_165 = tmp_156*(Blending_DF_Tetrahedron_3_4*Blending_DF_Tetrahedron_3_8 - Blending_DF_Tetrahedron_3_5*Blending_DF_Tetrahedron_3_7);
      real_t tmp_166 = tmp_165*tmp_41;
      real_t tmp_167 = tmp_165*tmp_44;
      real_t tmp_168 = tmp_165*tmp_46;
      real_t tmp_169 = -1.0*tmp_158 - 1.0*tmp_159 - 1.0*tmp_160 - 1.0*tmp_162 - 1.0*tmp_163 - 1.0*tmp_164 - 1.0*tmp_166 - 1.0*tmp_167 - 1.0*tmp_168;
      real_t tmp_170 = 0.5*tmp_156;
      real_t tmp_171 = tmp_170*(tmp_149 - tmp_153);
      real_t tmp_172 = tmp_171*tmp_8;
      real_t tmp_173 = tmp_171*tmp_30;
      real_t tmp_174 = tmp_171*tmp_32;
      real_t tmp_175 = tmp_170*(tmp_151 - tmp_152);
      real_t tmp_176 = tmp_175*tmp_34;
      real_t tmp_177 = tmp_175*tmp_37;
      real_t tmp_178 = tmp_175*tmp_39;
      real_t tmp_179 = tmp_170*(tmp_150 - tmp_154);
      real_t tmp_180 = tmp_179*tmp_41;
      real_t tmp_181 = tmp_179*tmp_44;
      real_t tmp_182 = tmp_179*tmp_46;
      real_t tmp_183 = 0.041666666666666657*tmp_74*std::abs(tmp_155);
      real_t tmp_184 = tmp_183*(-tmp_172 - tmp_173 - tmp_174 - tmp_176 - tmp_177 - tmp_178 - tmp_180 - tmp_181 - tmp_182);
      real_t tmp_185 = tmp_33 + tmp_40 + tmp_47;
      real_t tmp_186 = tmp_88 + tmp_92 + tmp_96;
      real_t tmp_187 = tmp_124 + tmp_128 + tmp_132;
      real_t tmp_188 = tmp_160 + tmp_164 + tmp_168;
      real_t tmp_189 = tmp_31 + tmp_38 + tmp_45;
      real_t tmp_190 = tmp_87 + tmp_91 + tmp_95;
      real_t tmp_191 = tmp_123 + tmp_127 + tmp_131;
      real_t tmp_192 = tmp_159 + tmp_163 + tmp_167;
      real_t tmp_193 = tmp_29 + tmp_36 + tmp_43;
      real_t tmp_194 = tmp_86 + tmp_90 + tmp_94;
      real_t tmp_195 = tmp_122 + tmp_126 + tmp_130;
      real_t tmp_196 = tmp_158 + tmp_162 + tmp_166;
      real_t tmp_197 = tmp_75*(tmp_53 + tmp_57 + tmp_61);
      real_t tmp_198 = tmp_111*(tmp_102 + tmp_106 + tmp_110);
      real_t tmp_199 = tmp_147*(tmp_138 + tmp_142 + tmp_146);
      real_t tmp_200 = tmp_183*(tmp_174 + tmp_178 + tmp_182);
      real_t tmp_201 = tmp_75*(tmp_52 + tmp_56 + tmp_60);
      real_t tmp_202 = tmp_111*(tmp_101 + tmp_105 + tmp_109);
      real_t tmp_203 = tmp_147*(tmp_137 + tmp_141 + tmp_145);
      real_t tmp_204 = tmp_183*(tmp_173 + tmp_177 + tmp_181);
      real_t tmp_205 = tmp_75*(tmp_51 + tmp_55 + tmp_59);
      real_t tmp_206 = tmp_111*(tmp_100 + tmp_104 + tmp_108);
      real_t tmp_207 = tmp_147*(tmp_136 + tmp_140 + tmp_144);
      real_t tmp_208 = tmp_183*(tmp_172 + tmp_176 + tmp_180);
      real_t a_0_0 = tmp_112*tmp_97 + tmp_133*tmp_148 + tmp_169*tmp_184 + tmp_48*tmp_76;
      real_t a_0_1 = tmp_112*tmp_186 + tmp_148*tmp_187 + tmp_184*tmp_188 + tmp_185*tmp_76;
      real_t a_0_2 = tmp_112*tmp_190 + tmp_148*tmp_191 + tmp_184*tmp_192 + tmp_189*tmp_76;
      real_t a_0_3 = tmp_112*tmp_194 + tmp_148*tmp_195 + tmp_184*tmp_196 + tmp_193*tmp_76;
      real_t a_1_0 = tmp_133*tmp_199 + tmp_169*tmp_200 + tmp_197*tmp_48 + tmp_198*tmp_97;
      real_t a_1_1 = tmp_185*tmp_197 + tmp_186*tmp_198 + tmp_187*tmp_199 + tmp_188*tmp_200;
      real_t a_1_2 = tmp_189*tmp_197 + tmp_190*tmp_198 + tmp_191*tmp_199 + tmp_192*tmp_200;
      real_t a_1_3 = tmp_193*tmp_197 + tmp_194*tmp_198 + tmp_195*tmp_199 + tmp_196*tmp_200;
      real_t a_2_0 = tmp_133*tmp_203 + tmp_169*tmp_204 + tmp_201*tmp_48 + tmp_202*tmp_97;
      real_t a_2_1 = tmp_185*tmp_201 + tmp_186*tmp_202 + tmp_187*tmp_203 + tmp_188*tmp_204;
      real_t a_2_2 = tmp_189*tmp_201 + tmp_190*tmp_202 + tmp_191*tmp_203 + tmp_192*tmp_204;
      real_t a_2_3 = tmp_193*tmp_201 + tmp_194*tmp_202 + tmp_195*tmp_203 + tmp_196*tmp_204;
      real_t a_3_0 = tmp_133*tmp_207 + tmp_169*tmp_208 + tmp_205*tmp_48 + tmp_206*tmp_97;
      real_t a_3_1 = tmp_185*tmp_205 + tmp_186*tmp_206 + tmp_187*tmp_207 + tmp_188*tmp_208;
      real_t a_3_2 = tmp_189*tmp_205 + tmp_190*tmp_206 + tmp_191*tmp_207 + tmp_192*tmp_208;
      real_t a_3_3 = tmp_193*tmp_205 + tmp_194*tmp_206 + tmp_195*tmp_207 + tmp_196*tmp_208;
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

   void p1_epsiloncc_0_2_blending_q2::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      Blending_DF_Tetrahedron( 0.13819660112501042*p_affine_0_0 + 0.13819660112501059*p_affine_1_0 + 0.13819660112501059*p_affine_2_0 + 0.58541019662496829*p_affine_3_0, 0.13819660112501042*p_affine_0_1 + 0.13819660112501059*p_affine_1_1 + 0.13819660112501059*p_affine_2_1 + 0.58541019662496829*p_affine_3_1, 0.13819660112501042*p_affine_0_2 + 0.13819660112501059*p_affine_1_2 + 0.13819660112501059*p_affine_2_2 + 0.58541019662496829*p_affine_3_2, &Blending_DF_Tetrahedron_0_0, &Blending_DF_Tetrahedron_0_1, &Blending_DF_Tetrahedron_0_2, &Blending_DF_Tetrahedron_0_3, &Blending_DF_Tetrahedron_0_4, &Blending_DF_Tetrahedron_0_5, &Blending_DF_Tetrahedron_0_6, &Blending_DF_Tetrahedron_0_7, &Blending_DF_Tetrahedron_0_8 );
      Blending_DF_Tetrahedron( 0.13819660112501048*p_affine_0_0 + 0.13819660112501059*p_affine_1_0 + 0.58541019662496829*p_affine_2_0 + 0.13819660112501059*p_affine_3_0, 0.13819660112501048*p_affine_0_1 + 0.13819660112501059*p_affine_1_1 + 0.58541019662496829*p_affine_2_1 + 0.13819660112501059*p_affine_3_1, 0.13819660112501048*p_affine_0_2 + 0.13819660112501059*p_affine_1_2 + 0.58541019662496829*p_affine_2_2 + 0.13819660112501059*p_affine_3_2, &Blending_DF_Tetrahedron_1_0, &Blending_DF_Tetrahedron_1_1, &Blending_DF_Tetrahedron_1_2, &Blending_DF_Tetrahedron_1_3, &Blending_DF_Tetrahedron_1_4, &Blending_DF_Tetrahedron_1_5, &Blending_DF_Tetrahedron_1_6, &Blending_DF_Tetrahedron_1_7, &Blending_DF_Tetrahedron_1_8 );
      Blending_DF_Tetrahedron( 0.13819660112501053*p_affine_0_0 + 0.58541019662496829*p_affine_1_0 + 0.13819660112501059*p_affine_2_0 + 0.13819660112501059*p_affine_3_0, 0.13819660112501053*p_affine_0_1 + 0.58541019662496829*p_affine_1_1 + 0.13819660112501059*p_affine_2_1 + 0.13819660112501059*p_affine_3_1, 0.13819660112501053*p_affine_0_2 + 0.58541019662496829*p_affine_1_2 + 0.13819660112501059*p_affine_2_2 + 0.13819660112501059*p_affine_3_2, &Blending_DF_Tetrahedron_2_0, &Blending_DF_Tetrahedron_2_1, &Blending_DF_Tetrahedron_2_2, &Blending_DF_Tetrahedron_2_3, &Blending_DF_Tetrahedron_2_4, &Blending_DF_Tetrahedron_2_5, &Blending_DF_Tetrahedron_2_6, &Blending_DF_Tetrahedron_2_7, &Blending_DF_Tetrahedron_2_8 );
      Blending_DF_Tetrahedron( 0.58541019662496807*p_affine_0_0 + 0.13819660112501059*p_affine_1_0 + 0.13819660112501059*p_affine_2_0 + 0.13819660112501059*p_affine_3_0, 0.58541019662496807*p_affine_0_1 + 0.13819660112501059*p_affine_1_1 + 0.13819660112501059*p_affine_2_1 + 0.13819660112501059*p_affine_3_1, 0.58541019662496807*p_affine_0_2 + 0.13819660112501059*p_affine_1_2 + 0.13819660112501059*p_affine_2_2 + 0.13819660112501059*p_affine_3_2, &Blending_DF_Tetrahedron_3_0, &Blending_DF_Tetrahedron_3_1, &Blending_DF_Tetrahedron_3_2, &Blending_DF_Tetrahedron_3_3, &Blending_DF_Tetrahedron_3_4, &Blending_DF_Tetrahedron_3_5, &Blending_DF_Tetrahedron_3_6, &Blending_DF_Tetrahedron_3_7, &Blending_DF_Tetrahedron_3_8 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = Blending_DF_Tetrahedron_0_0*Blending_DF_Tetrahedron_0_4;
      real_t tmp_10 = Blending_DF_Tetrahedron_0_1*Blending_DF_Tetrahedron_0_5;
      real_t tmp_11 = Blending_DF_Tetrahedron_0_2*Blending_DF_Tetrahedron_0_3;
      real_t tmp_12 = Blending_DF_Tetrahedron_0_0*Blending_DF_Tetrahedron_0_5;
      real_t tmp_13 = Blending_DF_Tetrahedron_0_1*Blending_DF_Tetrahedron_0_3;
      real_t tmp_14 = Blending_DF_Tetrahedron_0_2*Blending_DF_Tetrahedron_0_4;
      real_t tmp_15 = Blending_DF_Tetrahedron_0_6*tmp_10 - Blending_DF_Tetrahedron_0_6*tmp_14 + Blending_DF_Tetrahedron_0_7*tmp_11 - Blending_DF_Tetrahedron_0_7*tmp_12 - Blending_DF_Tetrahedron_0_8*tmp_13 + Blending_DF_Tetrahedron_0_8*tmp_9;
      real_t tmp_16 = -p_affine_0_2;
      real_t tmp_17 = p_affine_3_2 + tmp_16;
      real_t tmp_18 = p_affine_1_2 + tmp_16;
      real_t tmp_19 = p_affine_3_1 + tmp_2;
      real_t tmp_20 = tmp_19*tmp_5;
      real_t tmp_21 = p_affine_2_2 + tmp_16;
      real_t tmp_22 = p_affine_3_0 + tmp_0;
      real_t tmp_23 = tmp_22*tmp_6;
      real_t tmp_24 = tmp_1*tmp_19;
      real_t tmp_25 = tmp_22*tmp_3;
      real_t tmp_26 = 1.0 / (tmp_17*tmp_4 - tmp_17*tmp_7 + tmp_18*tmp_20 - tmp_18*tmp_25 + tmp_21*tmp_23 - tmp_21*tmp_24);
      real_t tmp_27 = tmp_26/tmp_15;
      real_t tmp_28 = tmp_27*(Blending_DF_Tetrahedron_0_3*Blending_DF_Tetrahedron_0_7 - Blending_DF_Tetrahedron_0_4*Blending_DF_Tetrahedron_0_6);
      real_t tmp_29 = tmp_28*tmp_8;
      real_t tmp_30 = tmp_23 - tmp_24;
      real_t tmp_31 = tmp_28*tmp_30;
      real_t tmp_32 = tmp_20 - tmp_25;
      real_t tmp_33 = tmp_28*tmp_32;
      real_t tmp_34 = -tmp_1*tmp_21 + tmp_18*tmp_5;
      real_t tmp_35 = tmp_27*(-Blending_DF_Tetrahedron_0_3*Blending_DF_Tetrahedron_0_8 + Blending_DF_Tetrahedron_0_5*Blending_DF_Tetrahedron_0_6);
      real_t tmp_36 = tmp_34*tmp_35;
      real_t tmp_37 = tmp_1*tmp_17 - tmp_18*tmp_22;
      real_t tmp_38 = tmp_35*tmp_37;
      real_t tmp_39 = -tmp_17*tmp_5 + tmp_21*tmp_22;
      real_t tmp_40 = tmp_35*tmp_39;
      real_t tmp_41 = -tmp_18*tmp_3 + tmp_21*tmp_6;
      real_t tmp_42 = tmp_27*(Blending_DF_Tetrahedron_0_4*Blending_DF_Tetrahedron_0_8 - Blending_DF_Tetrahedron_0_5*Blending_DF_Tetrahedron_0_7);
      real_t tmp_43 = tmp_41*tmp_42;
      real_t tmp_44 = -tmp_17*tmp_6 + tmp_18*tmp_19;
      real_t tmp_45 = tmp_42*tmp_44;
      real_t tmp_46 = tmp_17*tmp_3 - tmp_19*tmp_21;
      real_t tmp_47 = tmp_42*tmp_46;
      real_t tmp_48 = 0.5*tmp_27;
      real_t tmp_49 = tmp_48*(-tmp_13 + tmp_9);
      real_t tmp_50 = tmp_48*(tmp_11 - tmp_12);
      real_t tmp_51 = tmp_48*(tmp_10 - tmp_14);
      real_t tmp_52 = p_affine_0_0*p_affine_1_1;
      real_t tmp_53 = p_affine_0_0*p_affine_1_2;
      real_t tmp_54 = p_affine_2_1*p_affine_3_2;
      real_t tmp_55 = p_affine_0_1*p_affine_1_0;
      real_t tmp_56 = p_affine_0_1*p_affine_1_2;
      real_t tmp_57 = p_affine_2_2*p_affine_3_0;
      real_t tmp_58 = p_affine_0_2*p_affine_1_0;
      real_t tmp_59 = p_affine_0_2*p_affine_1_1;
      real_t tmp_60 = p_affine_2_0*p_affine_3_1;
      real_t tmp_61 = p_affine_2_2*p_affine_3_1;
      real_t tmp_62 = p_affine_2_0*p_affine_3_2;
      real_t tmp_63 = p_affine_2_1*p_affine_3_0;
      real_t tmp_64 = 2*std::abs(p_affine_0_0*tmp_54 - p_affine_0_0*tmp_61 + p_affine_0_1*tmp_57 - p_affine_0_1*tmp_62 + p_affine_0_2*tmp_60 - p_affine_0_2*tmp_63 - p_affine_1_0*tmp_54 + p_affine_1_0*tmp_61 - p_affine_1_1*tmp_57 + p_affine_1_1*tmp_62 - p_affine_1_2*tmp_60 + p_affine_1_2*tmp_63 + p_affine_2_0*tmp_56 - p_affine_2_0*tmp_59 - p_affine_2_1*tmp_53 + p_affine_2_1*tmp_58 + p_affine_2_2*tmp_52 - p_affine_2_2*tmp_55 - p_affine_3_0*tmp_56 + p_affine_3_0*tmp_59 + p_affine_3_1*tmp_53 - p_affine_3_1*tmp_58 - p_affine_3_2*tmp_52 + p_affine_3_2*tmp_55);
      real_t tmp_65 = 0.041666666666666657*tmp_64*(-tmp_30*tmp_49 - tmp_32*tmp_49 - tmp_34*tmp_50 - tmp_37*tmp_50 - tmp_39*tmp_50 - tmp_41*tmp_51 - tmp_44*tmp_51 - tmp_46*tmp_51 - tmp_49*tmp_8)*std::abs(tmp_15);
      real_t tmp_66 = Blending_DF_Tetrahedron_1_0*Blending_DF_Tetrahedron_1_4;
      real_t tmp_67 = Blending_DF_Tetrahedron_1_1*Blending_DF_Tetrahedron_1_5;
      real_t tmp_68 = Blending_DF_Tetrahedron_1_2*Blending_DF_Tetrahedron_1_3;
      real_t tmp_69 = Blending_DF_Tetrahedron_1_0*Blending_DF_Tetrahedron_1_5;
      real_t tmp_70 = Blending_DF_Tetrahedron_1_1*Blending_DF_Tetrahedron_1_3;
      real_t tmp_71 = Blending_DF_Tetrahedron_1_2*Blending_DF_Tetrahedron_1_4;
      real_t tmp_72 = Blending_DF_Tetrahedron_1_6*tmp_67 - Blending_DF_Tetrahedron_1_6*tmp_71 + Blending_DF_Tetrahedron_1_7*tmp_68 - Blending_DF_Tetrahedron_1_7*tmp_69 + Blending_DF_Tetrahedron_1_8*tmp_66 - Blending_DF_Tetrahedron_1_8*tmp_70;
      real_t tmp_73 = tmp_26/tmp_72;
      real_t tmp_74 = tmp_73*(Blending_DF_Tetrahedron_1_3*Blending_DF_Tetrahedron_1_7 - Blending_DF_Tetrahedron_1_4*Blending_DF_Tetrahedron_1_6);
      real_t tmp_75 = tmp_74*tmp_8;
      real_t tmp_76 = tmp_30*tmp_74;
      real_t tmp_77 = tmp_32*tmp_74;
      real_t tmp_78 = tmp_73*(-Blending_DF_Tetrahedron_1_3*Blending_DF_Tetrahedron_1_8 + Blending_DF_Tetrahedron_1_5*Blending_DF_Tetrahedron_1_6);
      real_t tmp_79 = tmp_34*tmp_78;
      real_t tmp_80 = tmp_37*tmp_78;
      real_t tmp_81 = tmp_39*tmp_78;
      real_t tmp_82 = tmp_73*(Blending_DF_Tetrahedron_1_4*Blending_DF_Tetrahedron_1_8 - Blending_DF_Tetrahedron_1_5*Blending_DF_Tetrahedron_1_7);
      real_t tmp_83 = tmp_41*tmp_82;
      real_t tmp_84 = tmp_44*tmp_82;
      real_t tmp_85 = tmp_46*tmp_82;
      real_t tmp_86 = 0.5*tmp_73;
      real_t tmp_87 = tmp_86*(tmp_66 - tmp_70);
      real_t tmp_88 = tmp_86*(tmp_68 - tmp_69);
      real_t tmp_89 = tmp_86*(tmp_67 - tmp_71);
      real_t tmp_90 = 0.041666666666666657*tmp_64*(-tmp_30*tmp_87 - tmp_32*tmp_87 - tmp_34*tmp_88 - tmp_37*tmp_88 - tmp_39*tmp_88 - tmp_41*tmp_89 - tmp_44*tmp_89 - tmp_46*tmp_89 - tmp_8*tmp_87)*std::abs(tmp_72);
      real_t tmp_91 = Blending_DF_Tetrahedron_2_0*Blending_DF_Tetrahedron_2_4;
      real_t tmp_92 = Blending_DF_Tetrahedron_2_1*Blending_DF_Tetrahedron_2_5;
      real_t tmp_93 = Blending_DF_Tetrahedron_2_2*Blending_DF_Tetrahedron_2_3;
      real_t tmp_94 = Blending_DF_Tetrahedron_2_0*Blending_DF_Tetrahedron_2_5;
      real_t tmp_95 = Blending_DF_Tetrahedron_2_1*Blending_DF_Tetrahedron_2_3;
      real_t tmp_96 = Blending_DF_Tetrahedron_2_2*Blending_DF_Tetrahedron_2_4;
      real_t tmp_97 = Blending_DF_Tetrahedron_2_6*tmp_92 - Blending_DF_Tetrahedron_2_6*tmp_96 + Blending_DF_Tetrahedron_2_7*tmp_93 - Blending_DF_Tetrahedron_2_7*tmp_94 + Blending_DF_Tetrahedron_2_8*tmp_91 - Blending_DF_Tetrahedron_2_8*tmp_95;
      real_t tmp_98 = tmp_26/tmp_97;
      real_t tmp_99 = tmp_98*(Blending_DF_Tetrahedron_2_3*Blending_DF_Tetrahedron_2_7 - Blending_DF_Tetrahedron_2_4*Blending_DF_Tetrahedron_2_6);
      real_t tmp_100 = tmp_8*tmp_99;
      real_t tmp_101 = tmp_30*tmp_99;
      real_t tmp_102 = tmp_32*tmp_99;
      real_t tmp_103 = tmp_98*(-Blending_DF_Tetrahedron_2_3*Blending_DF_Tetrahedron_2_8 + Blending_DF_Tetrahedron_2_5*Blending_DF_Tetrahedron_2_6);
      real_t tmp_104 = tmp_103*tmp_34;
      real_t tmp_105 = tmp_103*tmp_37;
      real_t tmp_106 = tmp_103*tmp_39;
      real_t tmp_107 = tmp_98*(Blending_DF_Tetrahedron_2_4*Blending_DF_Tetrahedron_2_8 - Blending_DF_Tetrahedron_2_5*Blending_DF_Tetrahedron_2_7);
      real_t tmp_108 = tmp_107*tmp_41;
      real_t tmp_109 = tmp_107*tmp_44;
      real_t tmp_110 = tmp_107*tmp_46;
      real_t tmp_111 = 0.5*tmp_98;
      real_t tmp_112 = tmp_111*(tmp_91 - tmp_95);
      real_t tmp_113 = tmp_111*(tmp_93 - tmp_94);
      real_t tmp_114 = tmp_111*(tmp_92 - tmp_96);
      real_t tmp_115 = 0.041666666666666657*tmp_64*(-tmp_112*tmp_30 - tmp_112*tmp_32 - tmp_112*tmp_8 - tmp_113*tmp_34 - tmp_113*tmp_37 - tmp_113*tmp_39 - tmp_114*tmp_41 - tmp_114*tmp_44 - tmp_114*tmp_46)*std::abs(tmp_97);
      real_t tmp_116 = Blending_DF_Tetrahedron_3_0*Blending_DF_Tetrahedron_3_4;
      real_t tmp_117 = Blending_DF_Tetrahedron_3_1*Blending_DF_Tetrahedron_3_5;
      real_t tmp_118 = Blending_DF_Tetrahedron_3_2*Blending_DF_Tetrahedron_3_3;
      real_t tmp_119 = Blending_DF_Tetrahedron_3_0*Blending_DF_Tetrahedron_3_5;
      real_t tmp_120 = Blending_DF_Tetrahedron_3_1*Blending_DF_Tetrahedron_3_3;
      real_t tmp_121 = Blending_DF_Tetrahedron_3_2*Blending_DF_Tetrahedron_3_4;
      real_t tmp_122 = Blending_DF_Tetrahedron_3_6*tmp_117 - Blending_DF_Tetrahedron_3_6*tmp_121 + Blending_DF_Tetrahedron_3_7*tmp_118 - Blending_DF_Tetrahedron_3_7*tmp_119 + Blending_DF_Tetrahedron_3_8*tmp_116 - Blending_DF_Tetrahedron_3_8*tmp_120;
      real_t tmp_123 = tmp_26/tmp_122;
      real_t tmp_124 = tmp_123*(Blending_DF_Tetrahedron_3_3*Blending_DF_Tetrahedron_3_7 - Blending_DF_Tetrahedron_3_4*Blending_DF_Tetrahedron_3_6);
      real_t tmp_125 = tmp_124*tmp_8;
      real_t tmp_126 = tmp_124*tmp_30;
      real_t tmp_127 = tmp_124*tmp_32;
      real_t tmp_128 = tmp_123*(-Blending_DF_Tetrahedron_3_3*Blending_DF_Tetrahedron_3_8 + Blending_DF_Tetrahedron_3_5*Blending_DF_Tetrahedron_3_6);
      real_t tmp_129 = tmp_128*tmp_34;
      real_t tmp_130 = tmp_128*tmp_37;
      real_t tmp_131 = tmp_128*tmp_39;
      real_t tmp_132 = tmp_123*(Blending_DF_Tetrahedron_3_4*Blending_DF_Tetrahedron_3_8 - Blending_DF_Tetrahedron_3_5*Blending_DF_Tetrahedron_3_7);
      real_t tmp_133 = tmp_132*tmp_41;
      real_t tmp_134 = tmp_132*tmp_44;
      real_t tmp_135 = tmp_132*tmp_46;
      real_t tmp_136 = 0.5*tmp_123;
      real_t tmp_137 = tmp_136*(tmp_116 - tmp_120);
      real_t tmp_138 = tmp_136*(tmp_118 - tmp_119);
      real_t tmp_139 = tmp_136*(tmp_117 - tmp_121);
      real_t tmp_140 = 0.041666666666666657*tmp_64*(-tmp_137*tmp_30 - tmp_137*tmp_32 - tmp_137*tmp_8 - tmp_138*tmp_34 - tmp_138*tmp_37 - tmp_138*tmp_39 - tmp_139*tmp_41 - tmp_139*tmp_44 - tmp_139*tmp_46)*std::abs(tmp_122);
      real_t a_0_0 = tmp_115*(-1.0*tmp_100 - 1.0*tmp_101 - 1.0*tmp_102 - 1.0*tmp_104 - 1.0*tmp_105 - 1.0*tmp_106 - 1.0*tmp_108 - 1.0*tmp_109 - 1.0*tmp_110) + tmp_140*(-1.0*tmp_125 - 1.0*tmp_126 - 1.0*tmp_127 - 1.0*tmp_129 - 1.0*tmp_130 - 1.0*tmp_131 - 1.0*tmp_133 - 1.0*tmp_134 - 1.0*tmp_135) + tmp_65*(-1.0*tmp_29 - 1.0*tmp_31 - 1.0*tmp_33 - 1.0*tmp_36 - 1.0*tmp_38 - 1.0*tmp_40 - 1.0*tmp_43 - 1.0*tmp_45 - 1.0*tmp_47) + tmp_90*(-1.0*tmp_75 - 1.0*tmp_76 - 1.0*tmp_77 - 1.0*tmp_79 - 1.0*tmp_80 - 1.0*tmp_81 - 1.0*tmp_83 - 1.0*tmp_84 - 1.0*tmp_85);
      real_t a_0_1 = tmp_115*(tmp_102 + tmp_106 + tmp_110) + tmp_140*(tmp_127 + tmp_131 + tmp_135) + tmp_65*(tmp_33 + tmp_40 + tmp_47) + tmp_90*(tmp_77 + tmp_81 + tmp_85);
      real_t a_0_2 = tmp_115*(tmp_101 + tmp_105 + tmp_109) + tmp_140*(tmp_126 + tmp_130 + tmp_134) + tmp_65*(tmp_31 + tmp_38 + tmp_45) + tmp_90*(tmp_76 + tmp_80 + tmp_84);
      real_t a_0_3 = tmp_115*(tmp_100 + tmp_104 + tmp_108) + tmp_140*(tmp_125 + tmp_129 + tmp_133) + tmp_65*(tmp_29 + tmp_36 + tmp_43) + tmp_90*(tmp_75 + tmp_79 + tmp_83);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_epsiloncc_0_2_blending_q2::Blending_DF_Tetrahedron( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
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

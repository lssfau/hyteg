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

#include "p1_to_p2_divt_1_affine_q2.hpp"

namespace hyteg {
namespace forms {

   void p1_to_p2_divt_1_affine_q2::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_0 = 0.16666666666666674;
      real_t tmp_1 = -p_affine_0_0;
      real_t tmp_2 = p_affine_1_0 + tmp_1;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = 1.0 / (tmp_2*(p_affine_2_1 + tmp_3) - (p_affine_1_1 + tmp_3)*(p_affine_2_0 + tmp_1));
      real_t tmp_5 = 0.66666666666666663;
      real_t tmp_6 = 2.6666666666666665;
      real_t tmp_7 = tmp_4*(tmp_5 + tmp_6 - 3);
      real_t tmp_8 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_9 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_10 = 0.16666666666666666*tmp_9;
      real_t tmp_11 = tmp_10*(-tmp_2*tmp_7 - tmp_7*tmp_8);
      real_t tmp_12 = 0.16666666666666671;
      real_t tmp_13 = 2.6666666666666665;
      real_t tmp_14 = 0.66666666666666663;
      real_t tmp_15 = tmp_4*(tmp_13 + tmp_14 - 3);
      real_t tmp_16 = 0.16666666666666666*tmp_9;
      real_t tmp_17 = tmp_16*(-tmp_15*tmp_2 - tmp_15*tmp_8);
      real_t tmp_18 = 0.66666666666666674;
      real_t tmp_19 = 0.66666666666666663;
      real_t tmp_20 = 0.66666666666666663;
      real_t tmp_21 = tmp_4*(tmp_19 + tmp_20 - 3);
      real_t tmp_22 = 0.16666666666666666*tmp_9;
      real_t tmp_23 = tmp_22*(-tmp_2*tmp_21 - tmp_21*tmp_8);
      real_t tmp_24 = tmp_0*tmp_10;
      real_t tmp_25 = tmp_4*tmp_8;
      real_t tmp_26 = tmp_25*(tmp_5 - 1);
      real_t tmp_27 = tmp_12*tmp_16;
      real_t tmp_28 = tmp_25*(tmp_13 - 1);
      real_t tmp_29 = tmp_18*tmp_22;
      real_t tmp_30 = tmp_25*(tmp_19 - 1);
      real_t tmp_31 = tmp_10*tmp_26;
      real_t tmp_32 = tmp_16*tmp_28;
      real_t tmp_33 = tmp_22*tmp_30;
      real_t tmp_34 = tmp_2*tmp_4;
      real_t tmp_35 = tmp_34*(tmp_6 - 1);
      real_t tmp_36 = tmp_34*(tmp_14 - 1);
      real_t tmp_37 = tmp_34*(tmp_20 - 1);
      real_t tmp_38 = tmp_10*tmp_35;
      real_t tmp_39 = tmp_16*tmp_36;
      real_t tmp_40 = tmp_22*tmp_37;
      real_t tmp_41 = tmp_34*tmp_5;
      real_t tmp_42 = tmp_25*tmp_6;
      real_t tmp_43 = tmp_10*(-tmp_41 - tmp_42);
      real_t tmp_44 = tmp_13*tmp_34;
      real_t tmp_45 = tmp_14*tmp_25;
      real_t tmp_46 = tmp_16*(-tmp_44 - tmp_45);
      real_t tmp_47 = tmp_19*tmp_34;
      real_t tmp_48 = tmp_20*tmp_25;
      real_t tmp_49 = tmp_22*(-tmp_47 - tmp_48);
      real_t tmp_50 = tmp_10*(-tmp_34*(-tmp_5 - 1.333333333333333) + tmp_42);
      real_t tmp_51 = tmp_16*(-tmp_34*(2.666666666666667 - tmp_13) + tmp_45);
      real_t tmp_52 = tmp_22*(-tmp_34*(2.666666666666667 - tmp_19) + tmp_48);
      real_t tmp_53 = tmp_10*(-tmp_25*(2.666666666666667 - tmp_6) + tmp_41);
      real_t tmp_54 = tmp_16*(-tmp_25*(-tmp_14 - 1.333333333333333) + tmp_44);
      real_t tmp_55 = tmp_22*(-tmp_25*(2.666666666666667 - tmp_20) + tmp_47);
      real_t a_0_0 = tmp_0*tmp_11 + tmp_12*tmp_17 + tmp_18*tmp_23;
      real_t a_0_1 = 0.16666666666666666*tmp_11 + 0.66666666666666663*tmp_17 + 0.16666666666666666*tmp_23;
      real_t a_0_2 = 0.66666666666666663*tmp_11 + 0.16666666666666666*tmp_17 + 0.16666666666666666*tmp_23;
      real_t a_1_0 = -tmp_24*tmp_26 - tmp_27*tmp_28 - tmp_29*tmp_30;
      real_t a_1_1 = -0.16666666666666666*tmp_31 - 0.66666666666666663*tmp_32 - 0.16666666666666666*tmp_33;
      real_t a_1_2 = -0.66666666666666663*tmp_31 - 0.16666666666666666*tmp_32 - 0.16666666666666666*tmp_33;
      real_t a_2_0 = -tmp_24*tmp_35 - tmp_27*tmp_36 - tmp_29*tmp_37;
      real_t a_2_1 = -0.16666666666666666*tmp_38 - 0.66666666666666663*tmp_39 - 0.16666666666666666*tmp_40;
      real_t a_2_2 = -0.66666666666666663*tmp_38 - 0.16666666666666666*tmp_39 - 0.16666666666666666*tmp_40;
      real_t a_3_0 = tmp_0*tmp_43 + tmp_12*tmp_46 + tmp_18*tmp_49;
      real_t a_3_1 = 0.16666666666666666*tmp_43 + 0.66666666666666663*tmp_46 + 0.16666666666666666*tmp_49;
      real_t a_3_2 = 0.66666666666666663*tmp_43 + 0.16666666666666666*tmp_46 + 0.16666666666666666*tmp_49;
      real_t a_4_0 = tmp_0*tmp_50 + tmp_12*tmp_51 + tmp_18*tmp_52;
      real_t a_4_1 = 0.16666666666666666*tmp_50 + 0.66666666666666663*tmp_51 + 0.16666666666666666*tmp_52;
      real_t a_4_2 = 0.66666666666666663*tmp_50 + 0.16666666666666666*tmp_51 + 0.16666666666666666*tmp_52;
      real_t a_5_0 = tmp_0*tmp_53 + tmp_12*tmp_54 + tmp_18*tmp_55;
      real_t a_5_1 = 0.16666666666666666*tmp_53 + 0.66666666666666663*tmp_54 + 0.16666666666666666*tmp_55;
      real_t a_5_2 = 0.66666666666666663*tmp_53 + 0.16666666666666666*tmp_54 + 0.16666666666666666*tmp_55;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(1, 0)) = a_1_0;
      (elMat(1, 1)) = a_1_1;
      (elMat(1, 2)) = a_1_2;
      (elMat(2, 0)) = a_2_0;
      (elMat(2, 1)) = a_2_1;
      (elMat(2, 2)) = a_2_2;
      (elMat(3, 0)) = a_3_0;
      (elMat(3, 1)) = a_3_1;
      (elMat(3, 2)) = a_3_2;
      (elMat(4, 0)) = a_4_0;
      (elMat(4, 1)) = a_4_1;
      (elMat(4, 2)) = a_4_2;
      (elMat(5, 0)) = a_5_0;
      (elMat(5, 1)) = a_5_1;
      (elMat(5, 2)) = a_5_2;
   }

   void p1_to_p2_divt_1_affine_q2::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = 1.0 / (tmp_1*(p_affine_2_1 + tmp_2) - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_4 = 0.33333333333333304*tmp_3;
      real_t tmp_5 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_6 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_7 = 0.16666666666666666*tmp_6*(-tmp_1*tmp_4 - tmp_4*tmp_5);
      real_t tmp_8 = 0.33333333333333315*tmp_3;
      real_t tmp_9 = 0.16666666666666666*tmp_6*(-tmp_1*tmp_8 - tmp_5*tmp_8);
      real_t tmp_10 = -1.666666666666667*tmp_3;
      real_t tmp_11 = 0.16666666666666666*tmp_6*(-tmp_1*tmp_10 - tmp_10*tmp_5);
      real_t a_0_0 = 0.66666666666666674*tmp_11 + 0.16666666666666674*tmp_7 + 0.16666666666666671*tmp_9;
      real_t a_0_1 = 0.16666666666666666*tmp_11 + 0.16666666666666666*tmp_7 + 0.66666666666666663*tmp_9;
      real_t a_0_2 = 0.16666666666666666*tmp_11 + 0.66666666666666663*tmp_7 + 0.16666666666666666*tmp_9;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_to_p2_divt_1_affine_q2::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 4 >& elMat ) const
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
      real_t tmp_0 = 0.13819660112501042;
      real_t tmp_1 = -p_affine_0_0;
      real_t tmp_2 = p_affine_2_0 + tmp_1;
      real_t tmp_3 = -p_affine_0_2;
      real_t tmp_4 = p_affine_1_2 + tmp_3;
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = p_affine_1_0 + tmp_1;
      real_t tmp_7 = p_affine_2_2 + tmp_3;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = tmp_5 - tmp_8;
      real_t tmp_10 = -p_affine_0_1;
      real_t tmp_11 = p_affine_2_1 + tmp_10;
      real_t tmp_12 = p_affine_3_2 + tmp_3;
      real_t tmp_13 = tmp_12*tmp_6;
      real_t tmp_14 = p_affine_3_1 + tmp_10;
      real_t tmp_15 = p_affine_1_1 + tmp_10;
      real_t tmp_16 = p_affine_3_0 + tmp_1;
      real_t tmp_17 = tmp_16*tmp_7;
      real_t tmp_18 = tmp_12*tmp_2;
      real_t tmp_19 = tmp_16*tmp_4;
      real_t tmp_20 = 1.0 / (tmp_11*tmp_13 - tmp_11*tmp_19 + tmp_14*tmp_5 - tmp_14*tmp_8 + tmp_15*tmp_17 - tmp_15*tmp_18);
      real_t tmp_21 = 0.55278640450004235;
      real_t tmp_22 = 0.55278640450004235;
      real_t tmp_23 = 2.3416407864998732;
      real_t tmp_24 = tmp_20*(tmp_21 + tmp_22 + tmp_23 - 3.0);
      real_t tmp_25 = tmp_13 - tmp_19;
      real_t tmp_26 = tmp_17 - tmp_18;
      real_t tmp_27 = p_affine_0_0*p_affine_1_1;
      real_t tmp_28 = p_affine_0_0*p_affine_1_2;
      real_t tmp_29 = p_affine_2_1*p_affine_3_2;
      real_t tmp_30 = p_affine_0_1*p_affine_1_0;
      real_t tmp_31 = p_affine_0_1*p_affine_1_2;
      real_t tmp_32 = p_affine_2_2*p_affine_3_0;
      real_t tmp_33 = p_affine_0_2*p_affine_1_0;
      real_t tmp_34 = p_affine_0_2*p_affine_1_1;
      real_t tmp_35 = p_affine_2_0*p_affine_3_1;
      real_t tmp_36 = p_affine_2_2*p_affine_3_1;
      real_t tmp_37 = p_affine_2_0*p_affine_3_2;
      real_t tmp_38 = p_affine_2_1*p_affine_3_0;
      real_t tmp_39 = std::abs(p_affine_0_0*tmp_29 - p_affine_0_0*tmp_36 + p_affine_0_1*tmp_32 - p_affine_0_1*tmp_37 + p_affine_0_2*tmp_35 - p_affine_0_2*tmp_38 - p_affine_1_0*tmp_29 + p_affine_1_0*tmp_36 - p_affine_1_1*tmp_32 + p_affine_1_1*tmp_37 - p_affine_1_2*tmp_35 + p_affine_1_2*tmp_38 + p_affine_2_0*tmp_31 - p_affine_2_0*tmp_34 - p_affine_2_1*tmp_28 + p_affine_2_1*tmp_33 + p_affine_2_2*tmp_27 - p_affine_2_2*tmp_30 - p_affine_3_0*tmp_31 + p_affine_3_0*tmp_34 + p_affine_3_1*tmp_28 - p_affine_3_1*tmp_33 - p_affine_3_2*tmp_27 + p_affine_3_2*tmp_30);
      real_t tmp_40 = 0.041666666666666657*tmp_39;
      real_t tmp_41 = tmp_40*(-tmp_24*tmp_25 - tmp_24*tmp_26 - tmp_24*tmp_9);
      real_t tmp_42 = 0.13819660112501048;
      real_t tmp_43 = 0.55278640450004235;
      real_t tmp_44 = 2.3416407864998732;
      real_t tmp_45 = 0.55278640450004235;
      real_t tmp_46 = tmp_20*(tmp_43 + tmp_44 + tmp_45 - 3.0);
      real_t tmp_47 = 0.041666666666666657*tmp_39;
      real_t tmp_48 = tmp_47*(-tmp_25*tmp_46 - tmp_26*tmp_46 - tmp_46*tmp_9);
      real_t tmp_49 = 0.13819660112501053;
      real_t tmp_50 = 2.3416407864998732;
      real_t tmp_51 = 0.55278640450004235;
      real_t tmp_52 = 0.55278640450004235;
      real_t tmp_53 = tmp_20*(tmp_50 + tmp_51 + tmp_52 - 3.0);
      real_t tmp_54 = 0.041666666666666657*tmp_39;
      real_t tmp_55 = tmp_54*(-tmp_25*tmp_53 - tmp_26*tmp_53 - tmp_53*tmp_9);
      real_t tmp_56 = 0.58541019662496807;
      real_t tmp_57 = 0.55278640450004235;
      real_t tmp_58 = 0.55278640450004235;
      real_t tmp_59 = 0.55278640450004235;
      real_t tmp_60 = tmp_20*(tmp_57 + tmp_58 + tmp_59 - 3.0);
      real_t tmp_61 = 0.041666666666666657*tmp_39;
      real_t tmp_62 = tmp_61*(-tmp_25*tmp_60 - tmp_26*tmp_60 - tmp_60*tmp_9);
      real_t tmp_63 = tmp_0*tmp_40;
      real_t tmp_64 = tmp_20*tmp_26;
      real_t tmp_65 = tmp_64*(tmp_21 - 1.0);
      real_t tmp_66 = tmp_42*tmp_47;
      real_t tmp_67 = tmp_64*(tmp_43 - 1.0);
      real_t tmp_68 = tmp_49*tmp_54;
      real_t tmp_69 = tmp_64*(tmp_50 - 1.0);
      real_t tmp_70 = tmp_56*tmp_61;
      real_t tmp_71 = tmp_64*(tmp_57 - 1.0);
      real_t tmp_72 = tmp_40*tmp_65;
      real_t tmp_73 = tmp_47*tmp_67;
      real_t tmp_74 = tmp_54*tmp_69;
      real_t tmp_75 = tmp_61*tmp_71;
      real_t tmp_76 = tmp_20*tmp_25;
      real_t tmp_77 = tmp_76*(tmp_22 - 1.0);
      real_t tmp_78 = tmp_76*(tmp_44 - 1.0);
      real_t tmp_79 = tmp_76*(tmp_51 - 1.0);
      real_t tmp_80 = tmp_76*(tmp_58 - 1.0);
      real_t tmp_81 = tmp_40*tmp_77;
      real_t tmp_82 = tmp_47*tmp_78;
      real_t tmp_83 = tmp_54*tmp_79;
      real_t tmp_84 = tmp_61*tmp_80;
      real_t tmp_85 = tmp_20*tmp_9;
      real_t tmp_86 = tmp_85*(tmp_23 - 1.0);
      real_t tmp_87 = tmp_85*(tmp_45 - 1.0);
      real_t tmp_88 = tmp_85*(tmp_52 - 1.0);
      real_t tmp_89 = tmp_85*(tmp_59 - 1.0);
      real_t tmp_90 = tmp_40*tmp_86;
      real_t tmp_91 = tmp_47*tmp_87;
      real_t tmp_92 = tmp_54*tmp_88;
      real_t tmp_93 = tmp_61*tmp_89;
      real_t tmp_94 = tmp_22*tmp_85;
      real_t tmp_95 = tmp_23*tmp_76;
      real_t tmp_96 = tmp_40*(-tmp_94 - tmp_95);
      real_t tmp_97 = tmp_44*tmp_85;
      real_t tmp_98 = tmp_45*tmp_76;
      real_t tmp_99 = tmp_47*(-tmp_97 - tmp_98);
      real_t tmp_100 = tmp_51*tmp_85;
      real_t tmp_101 = tmp_52*tmp_76;
      real_t tmp_102 = tmp_54*(-tmp_100 - tmp_101);
      real_t tmp_103 = tmp_58*tmp_85;
      real_t tmp_104 = tmp_59*tmp_76;
      real_t tmp_105 = tmp_61*(-tmp_103 - tmp_104);
      real_t tmp_106 = tmp_21*tmp_85;
      real_t tmp_107 = tmp_23*tmp_64;
      real_t tmp_108 = tmp_40*(-tmp_106 - tmp_107);
      real_t tmp_109 = tmp_43*tmp_85;
      real_t tmp_110 = tmp_45*tmp_64;
      real_t tmp_111 = tmp_47*(-tmp_109 - tmp_110);
      real_t tmp_112 = tmp_50*tmp_85;
      real_t tmp_113 = tmp_52*tmp_64;
      real_t tmp_114 = tmp_54*(-tmp_112 - tmp_113);
      real_t tmp_115 = tmp_57*tmp_85;
      real_t tmp_116 = tmp_59*tmp_64;
      real_t tmp_117 = tmp_61*(-tmp_115 - tmp_116);
      real_t tmp_118 = tmp_21*tmp_76;
      real_t tmp_119 = tmp_22*tmp_64;
      real_t tmp_120 = tmp_40*(-tmp_118 - tmp_119);
      real_t tmp_121 = tmp_43*tmp_76;
      real_t tmp_122 = tmp_44*tmp_64;
      real_t tmp_123 = tmp_47*(-tmp_121 - tmp_122);
      real_t tmp_124 = tmp_50*tmp_76;
      real_t tmp_125 = tmp_51*tmp_64;
      real_t tmp_126 = tmp_54*(-tmp_124 - tmp_125);
      real_t tmp_127 = tmp_57*tmp_76;
      real_t tmp_128 = tmp_58*tmp_64;
      real_t tmp_129 = tmp_61*(-tmp_127 - tmp_128);
      real_t tmp_130 = -tmp_22;
      real_t tmp_131 = 4.0 - tmp_21;
      real_t tmp_132 = tmp_40*(tmp_107 - tmp_85*(tmp_130 + tmp_131 - 4.6832815729997463) + tmp_95);
      real_t tmp_133 = -tmp_44;
      real_t tmp_134 = 4.0 - tmp_43;
      real_t tmp_135 = tmp_47*(tmp_110 - tmp_85*(tmp_133 + tmp_134 - 1.1055728090000847) + tmp_98);
      real_t tmp_136 = -tmp_51;
      real_t tmp_137 = 4.0 - tmp_50;
      real_t tmp_138 = tmp_54*(tmp_101 + tmp_113 - tmp_85*(tmp_136 + tmp_137 - 1.1055728090000847));
      real_t tmp_139 = -tmp_58;
      real_t tmp_140 = 4.0 - tmp_57;
      real_t tmp_141 = tmp_61*(tmp_104 + tmp_116 - tmp_85*(tmp_139 + tmp_140 - 1.1055728090000847));
      real_t tmp_142 = -tmp_23;
      real_t tmp_143 = tmp_40*(tmp_119 - tmp_76*(tmp_131 + tmp_142 - 1.1055728090000847) + tmp_94);
      real_t tmp_144 = -tmp_45;
      real_t tmp_145 = tmp_47*(tmp_122 - tmp_76*(tmp_134 + tmp_144 - 4.6832815729997463) + tmp_97);
      real_t tmp_146 = -tmp_52;
      real_t tmp_147 = tmp_54*(tmp_100 + tmp_125 - tmp_76*(tmp_137 + tmp_146 - 1.1055728090000847));
      real_t tmp_148 = -tmp_59;
      real_t tmp_149 = tmp_61*(tmp_103 + tmp_128 - tmp_76*(tmp_140 + tmp_148 - 1.1055728090000847));
      real_t tmp_150 = tmp_40*(tmp_106 + tmp_118 - tmp_64*(tmp_130 + tmp_142 + 2.8944271909999153));
      real_t tmp_151 = tmp_47*(tmp_109 + tmp_121 - tmp_64*(tmp_133 + tmp_144 + 2.8944271909999153));
      real_t tmp_152 = tmp_54*(tmp_112 + tmp_124 - tmp_64*(tmp_136 + tmp_146 - 0.68328157299974634));
      real_t tmp_153 = tmp_61*(tmp_115 + tmp_127 - tmp_64*(tmp_139 + tmp_148 + 2.8944271909999153));
      real_t a_0_0 = tmp_0*tmp_41 + tmp_42*tmp_48 + tmp_49*tmp_55 + tmp_56*tmp_62;
      real_t a_0_1 = 0.13819660112501059*tmp_41 + 0.13819660112501059*tmp_48 + 0.58541019662496829*tmp_55 + 0.13819660112501059*tmp_62;
      real_t a_0_2 = 0.13819660112501059*tmp_41 + 0.58541019662496829*tmp_48 + 0.13819660112501059*tmp_55 + 0.13819660112501059*tmp_62;
      real_t a_0_3 = 0.58541019662496829*tmp_41 + 0.13819660112501059*tmp_48 + 0.13819660112501059*tmp_55 + 0.13819660112501059*tmp_62;
      real_t a_1_0 = -tmp_63*tmp_65 - tmp_66*tmp_67 - tmp_68*tmp_69 - tmp_70*tmp_71;
      real_t a_1_1 = -0.13819660112501059*tmp_72 - 0.13819660112501059*tmp_73 - 0.58541019662496829*tmp_74 - 0.13819660112501059*tmp_75;
      real_t a_1_2 = -0.13819660112501059*tmp_72 - 0.58541019662496829*tmp_73 - 0.13819660112501059*tmp_74 - 0.13819660112501059*tmp_75;
      real_t a_1_3 = -0.58541019662496829*tmp_72 - 0.13819660112501059*tmp_73 - 0.13819660112501059*tmp_74 - 0.13819660112501059*tmp_75;
      real_t a_2_0 = -tmp_63*tmp_77 - tmp_66*tmp_78 - tmp_68*tmp_79 - tmp_70*tmp_80;
      real_t a_2_1 = -0.13819660112501059*tmp_81 - 0.13819660112501059*tmp_82 - 0.58541019662496829*tmp_83 - 0.13819660112501059*tmp_84;
      real_t a_2_2 = -0.13819660112501059*tmp_81 - 0.58541019662496829*tmp_82 - 0.13819660112501059*tmp_83 - 0.13819660112501059*tmp_84;
      real_t a_2_3 = -0.58541019662496829*tmp_81 - 0.13819660112501059*tmp_82 - 0.13819660112501059*tmp_83 - 0.13819660112501059*tmp_84;
      real_t a_3_0 = -tmp_63*tmp_86 - tmp_66*tmp_87 - tmp_68*tmp_88 - tmp_70*tmp_89;
      real_t a_3_1 = -0.13819660112501059*tmp_90 - 0.13819660112501059*tmp_91 - 0.58541019662496829*tmp_92 - 0.13819660112501059*tmp_93;
      real_t a_3_2 = -0.13819660112501059*tmp_90 - 0.58541019662496829*tmp_91 - 0.13819660112501059*tmp_92 - 0.13819660112501059*tmp_93;
      real_t a_3_3 = -0.58541019662496829*tmp_90 - 0.13819660112501059*tmp_91 - 0.13819660112501059*tmp_92 - 0.13819660112501059*tmp_93;
      real_t a_4_0 = tmp_0*tmp_96 + tmp_102*tmp_49 + tmp_105*tmp_56 + tmp_42*tmp_99;
      real_t a_4_1 = 0.58541019662496829*tmp_102 + 0.13819660112501059*tmp_105 + 0.13819660112501059*tmp_96 + 0.13819660112501059*tmp_99;
      real_t a_4_2 = 0.13819660112501059*tmp_102 + 0.13819660112501059*tmp_105 + 0.13819660112501059*tmp_96 + 0.58541019662496829*tmp_99;
      real_t a_4_3 = 0.13819660112501059*tmp_102 + 0.13819660112501059*tmp_105 + 0.58541019662496829*tmp_96 + 0.13819660112501059*tmp_99;
      real_t a_5_0 = tmp_0*tmp_108 + tmp_111*tmp_42 + tmp_114*tmp_49 + tmp_117*tmp_56;
      real_t a_5_1 = 0.13819660112501059*tmp_108 + 0.13819660112501059*tmp_111 + 0.58541019662496829*tmp_114 + 0.13819660112501059*tmp_117;
      real_t a_5_2 = 0.13819660112501059*tmp_108 + 0.58541019662496829*tmp_111 + 0.13819660112501059*tmp_114 + 0.13819660112501059*tmp_117;
      real_t a_5_3 = 0.58541019662496829*tmp_108 + 0.13819660112501059*tmp_111 + 0.13819660112501059*tmp_114 + 0.13819660112501059*tmp_117;
      real_t a_6_0 = tmp_0*tmp_120 + tmp_123*tmp_42 + tmp_126*tmp_49 + tmp_129*tmp_56;
      real_t a_6_1 = 0.13819660112501059*tmp_120 + 0.13819660112501059*tmp_123 + 0.58541019662496829*tmp_126 + 0.13819660112501059*tmp_129;
      real_t a_6_2 = 0.13819660112501059*tmp_120 + 0.58541019662496829*tmp_123 + 0.13819660112501059*tmp_126 + 0.13819660112501059*tmp_129;
      real_t a_6_3 = 0.58541019662496829*tmp_120 + 0.13819660112501059*tmp_123 + 0.13819660112501059*tmp_126 + 0.13819660112501059*tmp_129;
      real_t a_7_0 = tmp_0*tmp_132 + tmp_135*tmp_42 + tmp_138*tmp_49 + tmp_141*tmp_56;
      real_t a_7_1 = 0.13819660112501059*tmp_132 + 0.13819660112501059*tmp_135 + 0.58541019662496829*tmp_138 + 0.13819660112501059*tmp_141;
      real_t a_7_2 = 0.13819660112501059*tmp_132 + 0.58541019662496829*tmp_135 + 0.13819660112501059*tmp_138 + 0.13819660112501059*tmp_141;
      real_t a_7_3 = 0.58541019662496829*tmp_132 + 0.13819660112501059*tmp_135 + 0.13819660112501059*tmp_138 + 0.13819660112501059*tmp_141;
      real_t a_8_0 = tmp_0*tmp_143 + tmp_145*tmp_42 + tmp_147*tmp_49 + tmp_149*tmp_56;
      real_t a_8_1 = 0.13819660112501059*tmp_143 + 0.13819660112501059*tmp_145 + 0.58541019662496829*tmp_147 + 0.13819660112501059*tmp_149;
      real_t a_8_2 = 0.13819660112501059*tmp_143 + 0.58541019662496829*tmp_145 + 0.13819660112501059*tmp_147 + 0.13819660112501059*tmp_149;
      real_t a_8_3 = 0.58541019662496829*tmp_143 + 0.13819660112501059*tmp_145 + 0.13819660112501059*tmp_147 + 0.13819660112501059*tmp_149;
      real_t a_9_0 = tmp_0*tmp_150 + tmp_151*tmp_42 + tmp_152*tmp_49 + tmp_153*tmp_56;
      real_t a_9_1 = 0.13819660112501059*tmp_150 + 0.13819660112501059*tmp_151 + 0.58541019662496829*tmp_152 + 0.13819660112501059*tmp_153;
      real_t a_9_2 = 0.13819660112501059*tmp_150 + 0.58541019662496829*tmp_151 + 0.13819660112501059*tmp_152 + 0.13819660112501059*tmp_153;
      real_t a_9_3 = 0.58541019662496829*tmp_150 + 0.13819660112501059*tmp_151 + 0.13819660112501059*tmp_152 + 0.13819660112501059*tmp_153;
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
      (elMat(4, 0)) = a_4_0;
      (elMat(4, 1)) = a_4_1;
      (elMat(4, 2)) = a_4_2;
      (elMat(4, 3)) = a_4_3;
      (elMat(5, 0)) = a_5_0;
      (elMat(5, 1)) = a_5_1;
      (elMat(5, 2)) = a_5_2;
      (elMat(5, 3)) = a_5_3;
      (elMat(6, 0)) = a_6_0;
      (elMat(6, 1)) = a_6_1;
      (elMat(6, 2)) = a_6_2;
      (elMat(6, 3)) = a_6_3;
      (elMat(7, 0)) = a_7_0;
      (elMat(7, 1)) = a_7_1;
      (elMat(7, 2)) = a_7_2;
      (elMat(7, 3)) = a_7_3;
      (elMat(8, 0)) = a_8_0;
      (elMat(8, 1)) = a_8_1;
      (elMat(8, 2)) = a_8_2;
      (elMat(8, 3)) = a_8_3;
      (elMat(9, 0)) = a_9_0;
      (elMat(9, 1)) = a_9_1;
      (elMat(9, 2)) = a_9_2;
      (elMat(9, 3)) = a_9_3;
   }

   void p1_to_p2_divt_1_affine_q2::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_2_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_2;
      real_t tmp_3 = p_affine_1_2 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_1_0 + tmp_0;
      real_t tmp_6 = p_affine_2_2 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_0_1;
      real_t tmp_10 = p_affine_2_1 + tmp_9;
      real_t tmp_11 = p_affine_3_2 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = p_affine_3_1 + tmp_9;
      real_t tmp_14 = p_affine_1_1 + tmp_9;
      real_t tmp_15 = p_affine_3_0 + tmp_0;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_1*tmp_11;
      real_t tmp_18 = tmp_15*tmp_3;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_12 - tmp_10*tmp_18 + tmp_13*tmp_4 - tmp_13*tmp_7 + tmp_14*tmp_16 - tmp_14*tmp_17);
      real_t tmp_20 = 0.44721359549995809*tmp_19;
      real_t tmp_21 = tmp_12 - tmp_18;
      real_t tmp_22 = tmp_16 - tmp_17;
      real_t tmp_23 = p_affine_0_0*p_affine_1_1;
      real_t tmp_24 = p_affine_0_0*p_affine_1_2;
      real_t tmp_25 = p_affine_2_1*p_affine_3_2;
      real_t tmp_26 = p_affine_0_1*p_affine_1_0;
      real_t tmp_27 = p_affine_0_1*p_affine_1_2;
      real_t tmp_28 = p_affine_2_2*p_affine_3_0;
      real_t tmp_29 = p_affine_0_2*p_affine_1_0;
      real_t tmp_30 = p_affine_0_2*p_affine_1_1;
      real_t tmp_31 = p_affine_2_0*p_affine_3_1;
      real_t tmp_32 = p_affine_2_2*p_affine_3_1;
      real_t tmp_33 = p_affine_2_0*p_affine_3_2;
      real_t tmp_34 = p_affine_2_1*p_affine_3_0;
      real_t tmp_35 = std::abs(p_affine_0_0*tmp_25 - p_affine_0_0*tmp_32 + p_affine_0_1*tmp_28 - p_affine_0_1*tmp_33 + p_affine_0_2*tmp_31 - p_affine_0_2*tmp_34 - p_affine_1_0*tmp_25 + p_affine_1_0*tmp_32 - p_affine_1_1*tmp_28 + p_affine_1_1*tmp_33 - p_affine_1_2*tmp_31 + p_affine_1_2*tmp_34 + p_affine_2_0*tmp_27 - p_affine_2_0*tmp_30 - p_affine_2_1*tmp_24 + p_affine_2_1*tmp_29 + p_affine_2_2*tmp_23 - p_affine_2_2*tmp_26 - p_affine_3_0*tmp_27 + p_affine_3_0*tmp_30 + p_affine_3_1*tmp_24 - p_affine_3_1*tmp_29 - p_affine_3_2*tmp_23 + p_affine_3_2*tmp_26);
      real_t tmp_36 = 0.041666666666666657*tmp_35*(-tmp_20*tmp_21 - tmp_20*tmp_22 - tmp_20*tmp_8);
      real_t tmp_37 = 0.44721359549995809*tmp_19;
      real_t tmp_38 = 0.041666666666666657*tmp_35*(-tmp_21*tmp_37 - tmp_22*tmp_37 - tmp_37*tmp_8);
      real_t tmp_39 = 0.44721359549995787*tmp_19;
      real_t tmp_40 = 0.041666666666666657*tmp_35*(-tmp_21*tmp_39 - tmp_22*tmp_39 - tmp_39*tmp_8);
      real_t tmp_41 = -1.3416407864998727*tmp_19;
      real_t tmp_42 = 0.041666666666666657*tmp_35*(-tmp_21*tmp_41 - tmp_22*tmp_41 - tmp_41*tmp_8);
      real_t a_0_0 = 0.13819660112501042*tmp_36 + 0.13819660112501048*tmp_38 + 0.13819660112501053*tmp_40 + 0.58541019662496807*tmp_42;
      real_t a_0_1 = 0.13819660112501059*tmp_36 + 0.13819660112501059*tmp_38 + 0.58541019662496829*tmp_40 + 0.13819660112501059*tmp_42;
      real_t a_0_2 = 0.13819660112501059*tmp_36 + 0.58541019662496829*tmp_38 + 0.13819660112501059*tmp_40 + 0.13819660112501059*tmp_42;
      real_t a_0_3 = 0.58541019662496829*tmp_36 + 0.13819660112501059*tmp_38 + 0.13819660112501059*tmp_40 + 0.13819660112501059*tmp_42;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

} // namespace forms
} // namespace hyteg

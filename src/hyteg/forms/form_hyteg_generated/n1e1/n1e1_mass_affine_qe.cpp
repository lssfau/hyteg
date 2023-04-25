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

#include "n1e1_mass_affine_qe.hpp"

namespace hyteg {
namespace forms {

   void n1e1_mass_affine_qe::integrateAll( const std::array< Point3D, 4 >& coords, const std::array< int, 6 >& edgeDirections, Matrix< real_t, 6, 6 >& elMat ) const
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
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = (tmp_8*tmp_8);
      real_t tmp_10 = -p_affine_0_2;
      real_t tmp_11 = p_affine_3_2 + tmp_10;
      real_t tmp_12 = p_affine_3_1 + tmp_2;
      real_t tmp_13 = p_affine_1_2 + tmp_10;
      real_t tmp_14 = tmp_13*tmp_5;
      real_t tmp_15 = p_affine_3_0 + tmp_0;
      real_t tmp_16 = p_affine_2_2 + tmp_10;
      real_t tmp_17 = tmp_16*tmp_6;
      real_t tmp_18 = tmp_1*tmp_12;
      real_t tmp_19 = tmp_13*tmp_15;
      real_t tmp_20 = p_affine_0_0*p_affine_1_1;
      real_t tmp_21 = p_affine_0_0*p_affine_1_2;
      real_t tmp_22 = p_affine_2_1*p_affine_3_2;
      real_t tmp_23 = p_affine_0_1*p_affine_1_0;
      real_t tmp_24 = p_affine_0_1*p_affine_1_2;
      real_t tmp_25 = p_affine_2_2*p_affine_3_0;
      real_t tmp_26 = p_affine_0_2*p_affine_1_0;
      real_t tmp_27 = p_affine_0_2*p_affine_1_1;
      real_t tmp_28 = p_affine_2_0*p_affine_3_1;
      real_t tmp_29 = p_affine_2_2*p_affine_3_1;
      real_t tmp_30 = p_affine_2_0*p_affine_3_2;
      real_t tmp_31 = p_affine_2_1*p_affine_3_0;
      real_t tmp_32 = std::abs(p_affine_0_0*tmp_22 - p_affine_0_0*tmp_29 + p_affine_0_1*tmp_25 - p_affine_0_1*tmp_30 + p_affine_0_2*tmp_28 - p_affine_0_2*tmp_31 - p_affine_1_0*tmp_22 + p_affine_1_0*tmp_29 - p_affine_1_1*tmp_25 + p_affine_1_1*tmp_30 - p_affine_1_2*tmp_28 + p_affine_1_2*tmp_31 + p_affine_2_0*tmp_24 - p_affine_2_0*tmp_27 - p_affine_2_1*tmp_21 + p_affine_2_1*tmp_26 + p_affine_2_2*tmp_20 - p_affine_2_2*tmp_23 - p_affine_3_0*tmp_24 + p_affine_3_0*tmp_27 + p_affine_3_1*tmp_21 - p_affine_3_1*tmp_26 - p_affine_3_2*tmp_20 + p_affine_3_2*tmp_23)/((tmp_11*tmp_4 - tmp_11*tmp_7 + tmp_12*tmp_14 + tmp_15*tmp_17 - tmp_16*tmp_18 - tmp_19*tmp_3)*(tmp_11*tmp_4 - tmp_11*tmp_7 + tmp_12*tmp_14 + tmp_15*tmp_17 - tmp_16*tmp_18 - tmp_19*tmp_3));
      real_t tmp_33 = (1.0/60.0)*tmp_32;
      real_t tmp_34 = tmp_33*(walberla::real_c(edgeDirections[0])*walberla::real_c(edgeDirections[0]));
      real_t tmp_35 = tmp_15*tmp_6 - tmp_18;
      real_t tmp_36 = (tmp_35*tmp_35);
      real_t tmp_37 = -tmp_1*tmp_16 + tmp_14;
      real_t tmp_38 = (tmp_37*tmp_37);
      real_t tmp_39 = tmp_1*tmp_11 - tmp_19;
      real_t tmp_40 = (tmp_39*tmp_39);
      real_t tmp_41 = -tmp_13*tmp_3 + tmp_17;
      real_t tmp_42 = (tmp_41*tmp_41);
      real_t tmp_43 = -tmp_11*tmp_6 + tmp_12*tmp_13;
      real_t tmp_44 = (tmp_43*tmp_43);
      real_t tmp_45 = tmp_35*tmp_8;
      real_t tmp_46 = tmp_37*tmp_39;
      real_t tmp_47 = tmp_41*tmp_43;
      real_t tmp_48 = tmp_9*walberla::real_c(edgeDirections[0]);
      real_t tmp_49 = (1.0/120.0)*tmp_32;
      real_t tmp_50 = tmp_49*walberla::real_c(edgeDirections[1]);
      real_t tmp_51 = tmp_50*walberla::real_c(edgeDirections[0]);
      real_t tmp_52 = tmp_12*tmp_5 - tmp_15*tmp_3;
      real_t tmp_53 = tmp_52*tmp_8;
      real_t tmp_54 = tmp_35*tmp_52;
      real_t tmp_55 = tmp_33*walberla::real_c(edgeDirections[0]);
      real_t tmp_56 = tmp_55*walberla::real_c(edgeDirections[1]);
      real_t tmp_57 = -tmp_11*tmp_5 + tmp_15*tmp_16;
      real_t tmp_58 = tmp_37*tmp_57;
      real_t tmp_59 = tmp_39*tmp_57;
      real_t tmp_60 = tmp_59*walberla::real_c(edgeDirections[1]);
      real_t tmp_61 = tmp_11*tmp_3 - tmp_12*tmp_16;
      real_t tmp_62 = tmp_41*tmp_61;
      real_t tmp_63 = tmp_43*tmp_61;
      real_t tmp_64 = tmp_38*tmp_51 + tmp_42*tmp_51 - tmp_45*tmp_51 - tmp_46*tmp_51 - tmp_47*tmp_51 + tmp_48*tmp_50 - tmp_51*tmp_53 - tmp_51*tmp_58 - tmp_51*tmp_62 + tmp_54*tmp_56 + tmp_55*tmp_60 + tmp_56*tmp_63;
      real_t tmp_65 = tmp_49*walberla::real_c(edgeDirections[2]);
      real_t tmp_66 = tmp_65*walberla::real_c(edgeDirections[0]);
      real_t tmp_67 = tmp_40*walberla::real_c(edgeDirections[0]);
      real_t tmp_68 = tmp_55*walberla::real_c(edgeDirections[2]);
      real_t tmp_69 = -tmp_36*tmp_66 - tmp_44*tmp_66 + tmp_45*tmp_66 + tmp_46*tmp_66 + tmp_47*tmp_66 - tmp_53*tmp_68 + tmp_54*tmp_66 - tmp_58*tmp_68 + tmp_59*tmp_66 - tmp_62*tmp_68 + tmp_63*tmp_66 - tmp_65*tmp_67;
      real_t tmp_70 = tmp_33*walberla::real_c(edgeDirections[3]);
      real_t tmp_71 = tmp_70*walberla::real_c(edgeDirections[0]);
      real_t tmp_72 = tmp_49*walberla::real_c(edgeDirections[3]);
      real_t tmp_73 = tmp_72*walberla::real_c(edgeDirections[0]);
      real_t tmp_74 = -tmp_36*tmp_71 + tmp_38*tmp_71 + tmp_42*tmp_71 - tmp_44*tmp_71 - tmp_45*tmp_71 - tmp_46*tmp_71 - tmp_47*tmp_71 + tmp_48*tmp_70 + tmp_53*tmp_73 - tmp_54*tmp_71 + tmp_58*tmp_73 - tmp_59*tmp_71 + tmp_62*tmp_73 - tmp_63*tmp_71 - tmp_67*tmp_70;
      real_t tmp_75 = tmp_33*walberla::real_c(edgeDirections[4]);
      real_t tmp_76 = tmp_75*walberla::real_c(edgeDirections[0]);
      real_t tmp_77 = tmp_49*walberla::real_c(edgeDirections[4]);
      real_t tmp_78 = tmp_77*walberla::real_c(edgeDirections[0]);
      real_t tmp_79 = tmp_59*walberla::real_c(edgeDirections[0]);
      real_t tmp_80 = -tmp_36*tmp_76 + tmp_38*tmp_76 + tmp_42*tmp_76 - tmp_44*tmp_76 + tmp_45*tmp_76 + tmp_46*tmp_76 + tmp_47*tmp_76 + tmp_48*tmp_75 + tmp_53*tmp_76 - tmp_54*tmp_78 + tmp_58*tmp_76 + tmp_62*tmp_76 - tmp_63*tmp_78 - tmp_67*tmp_75 - tmp_77*tmp_79;
      real_t tmp_81 = tmp_49*walberla::real_c(edgeDirections[5]);
      real_t tmp_82 = tmp_81*walberla::real_c(edgeDirections[0]);
      real_t tmp_83 = tmp_33*walberla::real_c(edgeDirections[5]);
      real_t tmp_84 = tmp_83*walberla::real_c(edgeDirections[0]);
      real_t tmp_85 = -tmp_36*tmp_82 + tmp_38*tmp_82 + tmp_42*tmp_82 - tmp_44*tmp_82 + tmp_48*tmp_81 + tmp_53*tmp_84 - tmp_54*tmp_84 + tmp_58*tmp_84 + tmp_62*tmp_84 - tmp_63*tmp_84 - tmp_67*tmp_81 - tmp_79*tmp_83;
      real_t tmp_86 = tmp_33*(walberla::real_c(edgeDirections[1])*walberla::real_c(edgeDirections[1]));
      real_t tmp_87 = (tmp_52*tmp_52);
      real_t tmp_88 = (tmp_57*tmp_57);
      real_t tmp_89 = (tmp_61*tmp_61);
      real_t tmp_90 = tmp_50*walberla::real_c(edgeDirections[2]);
      real_t tmp_91 = tmp_45*walberla::real_c(edgeDirections[1]);
      real_t tmp_92 = tmp_33*walberla::real_c(edgeDirections[2]);
      real_t tmp_93 = tmp_92*walberla::real_c(edgeDirections[1]);
      real_t tmp_94 = tmp_46*tmp_93 + tmp_47*tmp_93 - tmp_53*tmp_90 - tmp_54*tmp_90 - tmp_58*tmp_90 - tmp_59*tmp_90 - tmp_62*tmp_90 - tmp_63*tmp_90 + tmp_87*tmp_90 + tmp_88*tmp_90 + tmp_89*tmp_90 + tmp_91*tmp_92;
      real_t tmp_95 = tmp_70*walberla::real_c(edgeDirections[1]);
      real_t tmp_96 = tmp_50*walberla::real_c(edgeDirections[3]);
      real_t tmp_97 = tmp_38*tmp_95 + tmp_42*tmp_95 + tmp_45*tmp_96 + tmp_46*tmp_96 + tmp_47*tmp_96 - tmp_53*tmp_95 - tmp_54*tmp_95 - tmp_58*tmp_95 - tmp_59*tmp_95 - tmp_62*tmp_95 - tmp_63*tmp_95 - tmp_87*tmp_95 - tmp_88*tmp_95 - tmp_89*tmp_95 + tmp_9*tmp_95;
      real_t tmp_98 = tmp_50*walberla::real_c(edgeDirections[4]);
      real_t tmp_99 = tmp_75*walberla::real_c(edgeDirections[1]);
      real_t tmp_100 = tmp_38*tmp_98 + tmp_42*tmp_98 + tmp_46*tmp_99 + tmp_47*tmp_99 - tmp_54*tmp_99 - tmp_60*tmp_75 - tmp_63*tmp_99 + tmp_75*tmp_91 - tmp_87*tmp_98 - tmp_88*tmp_98 - tmp_89*tmp_98 + tmp_9*tmp_98;
      real_t tmp_101 = tmp_83*walberla::real_c(edgeDirections[1]);
      real_t tmp_102 = tmp_50*walberla::real_c(edgeDirections[5]);
      real_t tmp_103 = tmp_101*tmp_38 + tmp_101*tmp_42 + tmp_101*tmp_45 + tmp_101*tmp_46 + tmp_101*tmp_47 + tmp_101*tmp_53 + tmp_101*tmp_58 + tmp_101*tmp_62 - tmp_101*tmp_87 - tmp_101*tmp_88 - tmp_101*tmp_89 + tmp_101*tmp_9 - tmp_102*tmp_54 - tmp_102*tmp_59 - tmp_102*tmp_63;
      real_t tmp_104 = tmp_33*(walberla::real_c(edgeDirections[2])*walberla::real_c(edgeDirections[2]));
      real_t tmp_105 = tmp_65*walberla::real_c(edgeDirections[3]);
      real_t tmp_106 = tmp_70*walberla::real_c(edgeDirections[2]);
      real_t tmp_107 = tmp_105*tmp_36 + tmp_105*tmp_40 + tmp_105*tmp_44 - tmp_105*tmp_87 - tmp_105*tmp_88 - tmp_105*tmp_89 + tmp_106*tmp_45 + tmp_106*tmp_46 + tmp_106*tmp_47 - tmp_106*tmp_53 - tmp_106*tmp_58 - tmp_106*tmp_62;
      real_t tmp_108 = tmp_75*walberla::real_c(edgeDirections[2]);
      real_t tmp_109 = tmp_65*walberla::real_c(edgeDirections[4]);
      real_t tmp_110 = tmp_108*tmp_36 + tmp_108*tmp_40 + tmp_108*tmp_44 - tmp_108*tmp_53 - tmp_108*tmp_54 - tmp_108*tmp_58 - tmp_108*tmp_59 - tmp_108*tmp_62 - tmp_108*tmp_63 - tmp_108*tmp_87 - tmp_108*tmp_88 - tmp_108*tmp_89 + tmp_109*tmp_45 + tmp_109*tmp_46 + tmp_109*tmp_47;
      real_t tmp_111 = tmp_83*walberla::real_c(edgeDirections[2]);
      real_t tmp_112 = tmp_65*walberla::real_c(edgeDirections[5]);
      real_t tmp_113 = tmp_111*tmp_36 + tmp_111*tmp_40 + tmp_111*tmp_44 + tmp_111*tmp_45 + tmp_111*tmp_46 + tmp_111*tmp_47 + tmp_111*tmp_54 + tmp_111*tmp_59 + tmp_111*tmp_63 - tmp_111*tmp_87 - tmp_111*tmp_88 - tmp_111*tmp_89 - tmp_112*tmp_53 - tmp_112*tmp_58 - tmp_112*tmp_62;
      real_t tmp_114 = (walberla::real_c(edgeDirections[3])*walberla::real_c(edgeDirections[3]));
      real_t tmp_115 = (1.0/20.0)*tmp_32;
      real_t tmp_116 = tmp_114*tmp_115;
      real_t tmp_117 = tmp_114*tmp_33;
      real_t tmp_118 = (1.0/30.0)*tmp_32;
      real_t tmp_119 = tmp_114*tmp_118;
      real_t tmp_120 = tmp_70*walberla::real_c(edgeDirections[4]);
      real_t tmp_121 = tmp_72*walberla::real_c(edgeDirections[4]);
      real_t tmp_122 = walberla::real_c(edgeDirections[3])*walberla::real_c(edgeDirections[4]);
      real_t tmp_123 = tmp_115*tmp_122;
      real_t tmp_124 = (1.0/40.0)*tmp_32;
      real_t tmp_125 = tmp_122*tmp_124;
      real_t tmp_126 = tmp_120*tmp_36 + tmp_120*tmp_38 + tmp_120*tmp_40 + tmp_120*tmp_42 + tmp_120*tmp_44 + tmp_120*tmp_9 + tmp_121*tmp_87 + tmp_121*tmp_88 + tmp_121*tmp_89 + tmp_123*tmp_45 + tmp_123*tmp_46 + tmp_123*tmp_47 + tmp_125*tmp_53 + tmp_125*tmp_54 + tmp_125*tmp_58 + tmp_125*tmp_59 + tmp_125*tmp_62 + tmp_125*tmp_63;
      real_t tmp_127 = tmp_70*walberla::real_c(edgeDirections[5]);
      real_t tmp_128 = tmp_81*walberla::real_c(edgeDirections[3]);
      real_t tmp_129 = walberla::real_c(edgeDirections[3])*walberla::real_c(edgeDirections[5]);
      real_t tmp_130 = tmp_124*tmp_129;
      real_t tmp_131 = tmp_115*tmp_129;
      real_t tmp_132 = tmp_127*tmp_38 + tmp_127*tmp_42 + tmp_127*tmp_87 + tmp_127*tmp_88 + tmp_127*tmp_89 + tmp_127*tmp_9 + tmp_128*tmp_36 + tmp_128*tmp_40 + tmp_128*tmp_44 + tmp_130*tmp_45 + tmp_130*tmp_46 + tmp_130*tmp_47 + tmp_130*tmp_54 + tmp_130*tmp_59 + tmp_130*tmp_63 + tmp_131*tmp_53 + tmp_131*tmp_58 + tmp_131*tmp_62;
      real_t tmp_133 = (walberla::real_c(edgeDirections[4])*walberla::real_c(edgeDirections[4]));
      real_t tmp_134 = tmp_133*tmp_33;
      real_t tmp_135 = tmp_115*tmp_133;
      real_t tmp_136 = tmp_118*tmp_133;
      real_t tmp_137 = tmp_81*walberla::real_c(edgeDirections[4]);
      real_t tmp_138 = tmp_75*walberla::real_c(edgeDirections[5]);
      real_t tmp_139 = walberla::real_c(edgeDirections[4])*walberla::real_c(edgeDirections[5]);
      real_t tmp_140 = tmp_124*tmp_139;
      real_t tmp_141 = tmp_115*tmp_139;
      real_t tmp_142 = tmp_137*tmp_38 + tmp_137*tmp_42 + tmp_137*tmp_9 + tmp_138*tmp_36 + tmp_138*tmp_40 + tmp_138*tmp_44 + tmp_138*tmp_87 + tmp_138*tmp_88 + tmp_138*tmp_89 + tmp_140*tmp_45 + tmp_140*tmp_46 + tmp_140*tmp_47 + tmp_140*tmp_53 + tmp_140*tmp_58 + tmp_140*tmp_62 + tmp_141*tmp_54 + tmp_141*tmp_59 + tmp_141*tmp_63;
      real_t tmp_143 = (walberla::real_c(edgeDirections[5])*walberla::real_c(edgeDirections[5]));
      real_t tmp_144 = tmp_143*tmp_33;
      real_t tmp_145 = tmp_115*tmp_143;
      real_t tmp_146 = tmp_118*tmp_143;
      real_t a_0_0 = tmp_34*tmp_36 + tmp_34*tmp_38 + tmp_34*tmp_40 + tmp_34*tmp_42 + tmp_34*tmp_44 - tmp_34*tmp_45 - tmp_34*tmp_46 - tmp_34*tmp_47 + tmp_34*tmp_9;
      real_t a_0_1 = tmp_64;
      real_t a_0_2 = tmp_69;
      real_t a_0_3 = tmp_74;
      real_t a_0_4 = tmp_80;
      real_t a_0_5 = tmp_85;
      real_t a_1_0 = tmp_64;
      real_t a_1_1 = tmp_38*tmp_86 + tmp_42*tmp_86 - tmp_53*tmp_86 - tmp_58*tmp_86 - tmp_62*tmp_86 + tmp_86*tmp_87 + tmp_86*tmp_88 + tmp_86*tmp_89 + tmp_86*tmp_9;
      real_t a_1_2 = tmp_94;
      real_t a_1_3 = tmp_97;
      real_t a_1_4 = tmp_100;
      real_t a_1_5 = tmp_103;
      real_t a_2_0 = tmp_69;
      real_t a_2_1 = tmp_94;
      real_t a_2_2 = tmp_104*tmp_36 + tmp_104*tmp_40 + tmp_104*tmp_44 - tmp_104*tmp_54 - tmp_104*tmp_59 - tmp_104*tmp_63 + tmp_104*tmp_87 + tmp_104*tmp_88 + tmp_104*tmp_89;
      real_t a_2_3 = tmp_107;
      real_t a_2_4 = tmp_110;
      real_t a_2_5 = tmp_113;
      real_t a_3_0 = tmp_74;
      real_t a_3_1 = tmp_97;
      real_t a_3_2 = tmp_107;
      real_t a_3_3 = tmp_116*tmp_38 + tmp_116*tmp_42 + tmp_116*tmp_45 + tmp_116*tmp_46 + tmp_116*tmp_47 + tmp_116*tmp_53 + tmp_116*tmp_58 + tmp_116*tmp_62 + tmp_116*tmp_9 + tmp_117*tmp_36 + tmp_117*tmp_40 + tmp_117*tmp_44 + tmp_117*tmp_87 + tmp_117*tmp_88 + tmp_117*tmp_89 + tmp_119*tmp_54 + tmp_119*tmp_59 + tmp_119*tmp_63;
      real_t a_3_4 = tmp_126;
      real_t a_3_5 = tmp_132;
      real_t a_4_0 = tmp_80;
      real_t a_4_1 = tmp_100;
      real_t a_4_2 = tmp_110;
      real_t a_4_3 = tmp_126;
      real_t a_4_4 = tmp_134*tmp_38 + tmp_134*tmp_42 + tmp_134*tmp_87 + tmp_134*tmp_88 + tmp_134*tmp_89 + tmp_134*tmp_9 + tmp_135*tmp_36 + tmp_135*tmp_40 + tmp_135*tmp_44 + tmp_135*tmp_45 + tmp_135*tmp_46 + tmp_135*tmp_47 + tmp_135*tmp_54 + tmp_135*tmp_59 + tmp_135*tmp_63 + tmp_136*tmp_53 + tmp_136*tmp_58 + tmp_136*tmp_62;
      real_t a_4_5 = tmp_142;
      real_t a_5_0 = tmp_85;
      real_t a_5_1 = tmp_103;
      real_t a_5_2 = tmp_113;
      real_t a_5_3 = tmp_132;
      real_t a_5_4 = tmp_142;
      real_t a_5_5 = tmp_144*tmp_36 + tmp_144*tmp_38 + tmp_144*tmp_40 + tmp_144*tmp_42 + tmp_144*tmp_44 + tmp_144*tmp_9 + tmp_145*tmp_53 + tmp_145*tmp_54 + tmp_145*tmp_58 + tmp_145*tmp_59 + tmp_145*tmp_62 + tmp_145*tmp_63 + tmp_145*tmp_87 + tmp_145*tmp_88 + tmp_145*tmp_89 + tmp_146*tmp_45 + tmp_146*tmp_46 + tmp_146*tmp_47;
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

} // namespace forms
} // namespace hyteg

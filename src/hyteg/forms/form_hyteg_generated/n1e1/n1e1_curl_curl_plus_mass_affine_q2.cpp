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

#include "n1e1_curl_curl_plus_mass_affine_q2.hpp"

namespace hyteg {
namespace forms {

   void n1e1_curl_curl_plus_mass_affine_q2::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 6, 6 >& elMat ) const
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
      real_t jac_affine_0_0 = -1.0*p_affine_0_0 + 1.0*p_affine_1_0;
      real_t jac_affine_0_1 = -1.0*p_affine_0_0 + 1.0*p_affine_2_0;
      real_t jac_affine_0_2 = -1.0*p_affine_0_0 + 1.0*p_affine_3_0;
      real_t jac_affine_1_0 = -1.0*p_affine_0_1 + 1.0*p_affine_1_1;
      real_t jac_affine_1_1 = -1.0*p_affine_0_1 + 1.0*p_affine_2_1;
      real_t jac_affine_1_2 = -1.0*p_affine_0_1 + 1.0*p_affine_3_1;
      real_t jac_affine_2_0 = -1.0*p_affine_0_2 + 1.0*p_affine_1_2;
      real_t jac_affine_2_1 = -1.0*p_affine_0_2 + 1.0*p_affine_2_2;
      real_t jac_affine_2_2 = -1.0*p_affine_0_2 + 1.0*p_affine_3_2;
      real_t tmp_73 = 1.0*jac_affine_1_1*jac_affine_2_2 - 1.0*jac_affine_1_2*jac_affine_2_1;
      real_t tmp_56 = -1.0*jac_affine_0_0*jac_affine_2_1 + 1.0*jac_affine_0_1*jac_affine_2_0;
      real_t tmp_57 = 1.0*jac_affine_1_0*jac_affine_2_1 - 1.0*jac_affine_1_1*jac_affine_2_0;
      real_t tmp_55 = 1.0*jac_affine_0_0*jac_affine_1_1 - 1.0*jac_affine_0_1*jac_affine_1_0;
      real_t tmp_37 = jac_affine_0_2*tmp_57 + jac_affine_1_2*tmp_56 + jac_affine_2_2*tmp_55;
      real_t np_0_arg_0 = tmp_37;
      real_t np_0 = 1.0 / (np_0_arg_0);
      real_t jac_affine_inv_0_0 = np_0*tmp_73;
      real_t tmp_74 = -1.0*jac_affine_0_1*jac_affine_2_2 + 1.0*jac_affine_0_2*jac_affine_2_1;
      real_t jac_affine_inv_0_1 = np_0*tmp_74;
      real_t tmp_75 = 1.0*jac_affine_0_1*jac_affine_1_2 - 1.0*jac_affine_0_2*jac_affine_1_1;
      real_t jac_affine_inv_0_2 = np_0*tmp_75;
      real_t tmp_76 = -1.0*jac_affine_1_0*jac_affine_2_2 + 1.0*jac_affine_1_2*jac_affine_2_0;
      real_t jac_affine_inv_1_0 = np_0*tmp_76;
      real_t tmp_77 = 1.0*jac_affine_0_0*jac_affine_2_2 - 1.0*jac_affine_0_2*jac_affine_2_0;
      real_t jac_affine_inv_1_1 = np_0*tmp_77;
      real_t tmp_78 = -1.0*jac_affine_0_0*jac_affine_1_2 + 1.0*jac_affine_0_2*jac_affine_1_0;
      real_t jac_affine_inv_1_2 = np_0*tmp_78;
      real_t jac_affine_inv_2_0 = np_0*tmp_57;
      real_t jac_affine_inv_2_1 = np_0*tmp_56;
      real_t jac_affine_inv_2_2 = np_0*tmp_55;
      real_t np_1_arg_0 = tmp_37;
      real_t np_1 = std::abs(np_1_arg_0);
      real_t abs_det_jac_affine = 1.0*np_1;
      real_t np_2 = 1.0 / (abs_det_jac_affine);
      real_t tmp_36 = (jac_affine_0_0*jac_affine_0_0) + (jac_affine_1_0*jac_affine_1_0) + (jac_affine_2_0*jac_affine_2_0);
      real_t tmp_3 = 0.16666666666666657;
      real_t tmp_12 = 0.016666666666666642;
      real_t tmp_6 = (jac_affine_inv_1_0*jac_affine_inv_1_0) + (jac_affine_inv_1_1*jac_affine_inv_1_1) + (jac_affine_inv_1_2*jac_affine_inv_1_2);
      real_t tmp_4 = (jac_affine_inv_2_0*jac_affine_inv_2_0) + (jac_affine_inv_2_1*jac_affine_inv_2_1) + (jac_affine_inv_2_2*jac_affine_inv_2_2);
      real_t tmp_11 = 0.016666666666666653;
      real_t tmp_53 = tmp_11*tmp_4 + tmp_12*tmp_6;
      real_t tmp_2 = jac_affine_inv_1_0*jac_affine_inv_2_0 + jac_affine_inv_1_1*jac_affine_inv_2_1 + jac_affine_inv_1_2*jac_affine_inv_2_2;
      real_t tmp_8 = 0.0083333333333333193;
      real_t tmp_81 = -2.0*tmp_2*tmp_8 + 1.0*tmp_53;
      real_t tmp_22 = abs_det_jac_affine*tmp_81 + 4.0*np_2*tmp_3*tmp_36;
      real_t a_0_0 = tmp_22;
      real_t tmp_30 = jac_affine_0_0*jac_affine_0_1 + jac_affine_1_0*jac_affine_1_1 + jac_affine_2_0*jac_affine_2_1;
      real_t tmp_9 = 0.0083333333333333228;
      real_t tmp_1 = jac_affine_inv_0_0*jac_affine_inv_2_0 + jac_affine_inv_0_1*jac_affine_inv_2_1 + jac_affine_inv_0_2*jac_affine_inv_2_2;
      real_t tmp_38 = tmp_1*tmp_8 + tmp_2*tmp_9;
      real_t tmp_0 = jac_affine_inv_0_0*jac_affine_inv_1_0 + jac_affine_inv_0_1*jac_affine_inv_1_1 + jac_affine_inv_0_2*jac_affine_inv_1_2;
      real_t tmp_7 = 0.0083333333333333263;
      real_t tmp_39 = tmp_0*tmp_12 + tmp_4*tmp_7;
      real_t tmp_87 = -1.0*tmp_38 + 1.0*tmp_39;
      real_t tmp_20 = abs_det_jac_affine*tmp_87 - 4.0*np_2*tmp_3*tmp_30;
      real_t a_0_1 = tmp_20;
      real_t tmp_40 = tmp_1*tmp_11 + tmp_6*tmp_9;
      real_t tmp_41 = tmp_0*tmp_8 + tmp_2*tmp_7;
      real_t tmp_93 = -1.0*tmp_40 + 1.0*tmp_41;
      real_t tmp_29 = jac_affine_0_0*jac_affine_0_2 + jac_affine_1_0*jac_affine_1_2 + jac_affine_2_0*jac_affine_2_2;
      real_t tmp_16 = abs_det_jac_affine*tmp_93 + 4.0*np_2*tmp_29*tmp_3;
      real_t a_0_2 = tmp_16;
      real_t tmp_85 = 1.0*tmp_38 - 1.0*tmp_39;
      real_t tmp_18 = abs_det_jac_affine*tmp_85 + 4.0*np_2*tmp_3*tmp_30;
      real_t tmp_13 = 0.041666666666666602;
      real_t tmp_89 = -1.0*tmp_13 + 2.0*tmp_8;
      real_t tmp_14 = 0.041666666666666644;
      real_t tmp_66 = 1.0*tmp_14*tmp_4 + tmp_2*tmp_89 - 1.0*tmp_53;
      real_t tmp_26 = abs_det_jac_affine*tmp_66 - 4.0*np_2*tmp_3*tmp_36 + tmp_18;
      real_t a_0_3 = tmp_26;
      real_t tmp_86 = -1.0*tmp_13*tmp_6 + 1.0*tmp_14*tmp_2;
      real_t tmp_92 = 1.0*tmp_40 - 1.0*tmp_41;
      real_t tmp_21 = abs_det_jac_affine*tmp_92 - 4.0*np_2*tmp_29*tmp_3;
      real_t tmp_47 = abs_det_jac_affine*tmp_86 + tmp_21 + tmp_22;
      real_t a_0_4 = tmp_47;
      real_t tmp_45 = -1.0*abs_det_jac_affine*tmp_0*tmp_13 + tmp_20;
      real_t tmp_60 = 1.0*abs_det_jac_affine*tmp_1*tmp_14 + tmp_16 + tmp_45;
      real_t a_0_5 = tmp_60;
      real_t a_1_0 = tmp_20;
      real_t tmp_34 = (jac_affine_0_1*jac_affine_0_1) + (jac_affine_1_1*jac_affine_1_1) + (jac_affine_2_1*jac_affine_2_1);
      real_t tmp_10 = 0.016666666666666656;
      real_t tmp_5 = (jac_affine_inv_0_0*jac_affine_inv_0_0) + (jac_affine_inv_0_1*jac_affine_inv_0_1) + (jac_affine_inv_0_2*jac_affine_inv_0_2);
      real_t tmp_54 = tmp_10*tmp_4 + tmp_12*tmp_5;
      real_t tmp_90 = -2.0*tmp_1*tmp_9 + 1.0*tmp_54;
      real_t tmp_24 = abs_det_jac_affine*tmp_90 + 4.0*np_2*tmp_3*tmp_34;
      real_t a_1_1 = tmp_24;
      real_t tmp_43 = tmp_0*tmp_9 + tmp_1*tmp_7;
      real_t tmp_42 = tmp_10*tmp_2 + tmp_5*tmp_8;
      real_t tmp_83 = 1.0*tmp_42 - 1.0*tmp_43;
      real_t tmp_28 = jac_affine_0_1*jac_affine_0_2 + jac_affine_1_1*jac_affine_1_2 + jac_affine_2_1*jac_affine_2_2;
      real_t tmp_19 = abs_det_jac_affine*tmp_83 - 4.0*np_2*tmp_28*tmp_3;
      real_t a_1_2 = tmp_19;
      real_t tmp_15 = 0.041666666666666644;
      real_t tmp_103 = -1.0*tmp_13 + 2.0*tmp_9;
      real_t tmp_64 = tmp_1*tmp_103 + 1.0*tmp_15*tmp_4 - 1.0*tmp_54;
      real_t tmp_25 = abs_det_jac_affine*tmp_64 - 4.0*np_2*tmp_3*tmp_34 + tmp_18;
      real_t a_1_3 = tmp_25;
      real_t tmp_68 = 1.0*tmp_15*tmp_2 - 1.0*tmp_42 + 1.0*tmp_43;
      real_t tmp_17 = abs_det_jac_affine*tmp_68 + 4.0*np_2*tmp_28*tmp_3;
      real_t tmp_100 = tmp_17 + tmp_45;
      real_t a_1_4 = tmp_100;
      real_t tmp_91 = 1.0*tmp_1*tmp_15 - 1.0*tmp_13*tmp_5;
      real_t tmp_46 = abs_det_jac_affine*tmp_91 + tmp_19 + tmp_24;
      real_t a_1_5 = tmp_46;
      real_t a_2_0 = tmp_16;
      real_t a_2_1 = tmp_19;
      real_t tmp_52 = tmp_10*tmp_6 + tmp_11*tmp_5;
      real_t tmp_82 = -2.0*tmp_0*tmp_7 + 1.0*tmp_52;
      real_t tmp_35 = (jac_affine_0_2*jac_affine_0_2) + (jac_affine_1_2*jac_affine_1_2) + (jac_affine_2_2*jac_affine_2_2);
      real_t tmp_23 = abs_det_jac_affine*tmp_82 + 4.0*np_2*tmp_3*tmp_35;
      real_t a_2_2 = tmp_23;
      real_t tmp_98 = tmp_17 + tmp_21;
      real_t tmp_59 = -1.0*abs_det_jac_affine*tmp_1*tmp_14 + tmp_98;
      real_t a_2_3 = tmp_59;
      real_t tmp_84 = -1.0*tmp_14 + 2.0*tmp_7;
      real_t tmp_65 = tmp_0*tmp_84 + 1.0*tmp_15*tmp_6 - 1.0*tmp_52;
      real_t tmp_27 = abs_det_jac_affine*tmp_65 - 4.0*np_2*tmp_3*tmp_35 + tmp_16;
      real_t a_2_4 = tmp_27;
      real_t tmp_88 = 1.0*tmp_0*tmp_15 - 1.0*tmp_14*tmp_5;
      real_t tmp_48 = abs_det_jac_affine*tmp_88 + tmp_19 + tmp_23;
      real_t a_2_5 = tmp_48;
      real_t a_3_0 = tmp_26;
      real_t a_3_1 = tmp_25;
      real_t a_3_2 = tmp_59;
      real_t tmp_97 = tmp_1 + tmp_2;
      real_t tmp_70 = tmp_13*tmp_97 + tmp_39;
      real_t tmp_61 = -2.0*tmp_15 + 1.0*tmp_3;
      real_t tmp_72 = -2.0*tmp_14;
      real_t tmp_80 = tmp_61 + tmp_72;
      real_t tmp_49 = -2.0*tmp_38 + tmp_4*tmp_80 + 2.0*tmp_70;
      real_t a_3_3 = abs_det_jac_affine*tmp_49 - 8.0*np_2*tmp_3*tmp_30 + tmp_22 + tmp_24;
      real_t tmp_104 = tmp_0 + tmp_6;
      real_t tmp_79 = 1.0*tmp_1 - 1.0*tmp_2;
      real_t tmp_31 = 1.0*tmp_104*tmp_13 + tmp_14*tmp_79 + tmp_2*tmp_61;
      real_t tmp_63 = abs_det_jac_affine*tmp_31 + tmp_16 + tmp_19 + tmp_26;
      real_t a_3_4 = tmp_63;
      real_t tmp_101 = tmp_0 + tmp_5;
      real_t tmp_71 = -1.0*tmp_15 + 1.0*tmp_3 + tmp_72;
      real_t tmp_32 = tmp_1*tmp_71 + 1.0*tmp_101*tmp_13;
      real_t tmp_58 = abs_det_jac_affine*tmp_32 + tmp_25 + tmp_98;
      real_t a_3_5 = tmp_58;
      real_t a_4_0 = tmp_47;
      real_t a_4_1 = tmp_100;
      real_t a_4_2 = tmp_27;
      real_t a_4_3 = tmp_63;
      real_t tmp_96 = tmp_0 + tmp_2;
      real_t tmp_69 = tmp_14*tmp_96 + tmp_40;
      real_t tmp_94 = tmp_15*tmp_6 + tmp_41;
      real_t tmp_44 = -2.0*tmp_13 + 1.0*tmp_3;
      real_t tmp_50 = tmp_44*tmp_6 + 2.0*tmp_69 - 2.0*tmp_94;
      real_t a_4_4 = abs_det_jac_affine*tmp_50 - 8.0*np_2*tmp_29*tmp_3 + tmp_22 + tmp_23;
      real_t tmp_102 = tmp_1 + tmp_5;
      real_t tmp_95 = -1.0*tmp_15 + tmp_44;
      real_t tmp_33 = tmp_0*tmp_95 + 1.0*tmp_102*tmp_14;
      real_t tmp_62 = abs_det_jac_affine*tmp_33 + tmp_17 + tmp_20 + tmp_27;
      real_t a_4_5 = tmp_62;
      real_t a_5_0 = tmp_60;
      real_t a_5_1 = tmp_46;
      real_t a_5_2 = tmp_48;
      real_t a_5_3 = tmp_58;
      real_t a_5_4 = tmp_62;
      real_t tmp_99 = tmp_0 + tmp_1;
      real_t tmp_67 = tmp_15*tmp_99 + tmp_42;
      real_t tmp_105 = tmp_44 + tmp_72;
      real_t tmp_51 = tmp_105*tmp_5 - 2.0*tmp_43 + 2.0*tmp_67;
      real_t a_5_5 = abs_det_jac_affine*tmp_51 - 8.0*np_2*tmp_28*tmp_3 + tmp_23 + tmp_24;
      elMat(0,0) = a_0_0;
      elMat(0,1) = a_0_1;
      elMat(0,2) = a_0_2;
      elMat(0,3) = a_0_3;
      elMat(0,4) = a_0_4;
      elMat(0,5) = a_0_5;
      elMat(1,0) = a_1_0;
      elMat(1,1) = a_1_1;
      elMat(1,2) = a_1_2;
      elMat(1,3) = a_1_3;
      elMat(1,4) = a_1_4;
      elMat(1,5) = a_1_5;
      elMat(2,0) = a_2_0;
      elMat(2,1) = a_2_1;
      elMat(2,2) = a_2_2;
      elMat(2,3) = a_2_3;
      elMat(2,4) = a_2_4;
      elMat(2,5) = a_2_5;
      elMat(3,0) = a_3_0;
      elMat(3,1) = a_3_1;
      elMat(3,2) = a_3_2;
      elMat(3,3) = a_3_3;
      elMat(3,4) = a_3_4;
      elMat(3,5) = a_3_5;
      elMat(4,0) = a_4_0;
      elMat(4,1) = a_4_1;
      elMat(4,2) = a_4_2;
      elMat(4,3) = a_4_3;
      elMat(4,4) = a_4_4;
      elMat(4,5) = a_4_5;
      elMat(5,0) = a_5_0;
      elMat(5,1) = a_5_1;
      elMat(5,2) = a_5_2;
      elMat(5,3) = a_5_3;
      elMat(5,4) = a_5_4;
      elMat(5,5) = a_5_5;
   }

} // namespace forms
} // namespace hyteg

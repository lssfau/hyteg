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

#include "n1e1_curl_curl_affine_qe.hpp"

namespace hyteg {
namespace forms {

   void n1e1_curl_curl_affine_qe::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 6, 6 >& elMat ) const
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
      real_t tmp_0 = p_affine_0_0*p_affine_1_1;
      real_t tmp_1 = p_affine_0_0*p_affine_1_2;
      real_t tmp_2 = p_affine_2_1*p_affine_3_2;
      real_t tmp_3 = p_affine_0_1*p_affine_1_0;
      real_t tmp_4 = p_affine_0_1*p_affine_1_2;
      real_t tmp_5 = p_affine_2_2*p_affine_3_0;
      real_t tmp_6 = p_affine_0_2*p_affine_1_0;
      real_t tmp_7 = p_affine_0_2*p_affine_1_1;
      real_t tmp_8 = p_affine_2_0*p_affine_3_1;
      real_t tmp_9 = p_affine_2_2*p_affine_3_1;
      real_t tmp_10 = p_affine_2_0*p_affine_3_2;
      real_t tmp_11 = p_affine_2_1*p_affine_3_0;
      real_t tmp_12 = 1.0 / (std::abs(p_affine_0_0*tmp_2 - p_affine_0_0*tmp_9 - p_affine_0_1*tmp_10 + p_affine_0_1*tmp_5 - p_affine_0_2*tmp_11 + p_affine_0_2*tmp_8 - p_affine_1_0*tmp_2 + p_affine_1_0*tmp_9 + p_affine_1_1*tmp_10 - p_affine_1_1*tmp_5 + p_affine_1_2*tmp_11 - p_affine_1_2*tmp_8 + p_affine_2_0*tmp_4 - p_affine_2_0*tmp_7 - p_affine_2_1*tmp_1 + p_affine_2_1*tmp_6 + p_affine_2_2*tmp_0 - p_affine_2_2*tmp_3 - p_affine_3_0*tmp_4 + p_affine_3_0*tmp_7 + p_affine_3_1*tmp_1 - p_affine_3_1*tmp_6 - p_affine_3_2*tmp_0 + p_affine_3_2*tmp_3));
      real_t tmp_13 = (4.0/3.0)*tmp_12;
      real_t tmp_14 = p_affine_0_0*p_affine_1_0;
      real_t tmp_15 = p_affine_0_1*p_affine_1_1;
      real_t tmp_16 = p_affine_0_2*p_affine_1_2;
      real_t tmp_17 = (2.0/3.0)*tmp_12;
      real_t tmp_18 = (p_affine_0_0*p_affine_0_0)*tmp_17;
      real_t tmp_19 = (p_affine_0_1*p_affine_0_1)*tmp_17;
      real_t tmp_20 = (p_affine_0_2*p_affine_0_2)*tmp_17;
      real_t tmp_21 = tmp_18 + tmp_19 + tmp_20;
      real_t tmp_22 = (p_affine_1_0*p_affine_1_0)*tmp_17;
      real_t tmp_23 = (p_affine_1_1*p_affine_1_1)*tmp_17;
      real_t tmp_24 = (p_affine_1_2*p_affine_1_2)*tmp_17;
      real_t tmp_25 = tmp_22 + tmp_23 + tmp_24;
      real_t tmp_26 = tmp_14*tmp_17;
      real_t tmp_27 = tmp_15*tmp_17;
      real_t tmp_28 = tmp_16*tmp_17;
      real_t tmp_29 = tmp_26 + tmp_27 + tmp_28;
      real_t tmp_30 = p_affine_2_0*tmp_17;
      real_t tmp_31 = p_affine_0_0*tmp_30;
      real_t tmp_32 = p_affine_2_1*tmp_17;
      real_t tmp_33 = p_affine_0_1*tmp_32;
      real_t tmp_34 = p_affine_2_2*tmp_17;
      real_t tmp_35 = p_affine_0_2*tmp_34;
      real_t tmp_36 = tmp_31 + tmp_33 + tmp_35;
      real_t tmp_37 = p_affine_1_0*tmp_30;
      real_t tmp_38 = p_affine_1_1*tmp_32;
      real_t tmp_39 = p_affine_1_2*tmp_34;
      real_t tmp_40 = -tmp_37 - tmp_38 - tmp_39;
      real_t tmp_41 = tmp_36 + tmp_40;
      real_t tmp_42 = -tmp_18 - tmp_19 - tmp_20;
      real_t tmp_43 = tmp_29 + tmp_41 + tmp_42;
      real_t tmp_44 = -tmp_26 - tmp_27 - tmp_28;
      real_t tmp_45 = p_affine_3_0*tmp_17;
      real_t tmp_46 = p_affine_0_0*tmp_45;
      real_t tmp_47 = p_affine_3_1*tmp_17;
      real_t tmp_48 = p_affine_0_1*tmp_47;
      real_t tmp_49 = p_affine_3_2*tmp_17;
      real_t tmp_50 = p_affine_0_2*tmp_49;
      real_t tmp_51 = -tmp_46 - tmp_48 - tmp_50;
      real_t tmp_52 = p_affine_1_0*tmp_45;
      real_t tmp_53 = p_affine_1_1*tmp_47;
      real_t tmp_54 = p_affine_1_2*tmp_49;
      real_t tmp_55 = tmp_52 + tmp_53 + tmp_54;
      real_t tmp_56 = tmp_51 + tmp_55;
      real_t tmp_57 = tmp_21 + tmp_44 + tmp_56;
      real_t tmp_58 = -tmp_31 - tmp_33 - tmp_35;
      real_t tmp_59 = tmp_29 + tmp_58;
      real_t tmp_60 = tmp_37 + tmp_38 + tmp_39;
      real_t tmp_61 = -tmp_22 - tmp_23 - tmp_24 + tmp_60;
      real_t tmp_62 = tmp_59 + tmp_61;
      real_t tmp_63 = -tmp_52 - tmp_53 - tmp_54;
      real_t tmp_64 = tmp_46 + tmp_48 + tmp_50;
      real_t tmp_65 = tmp_44 + tmp_64;
      real_t tmp_66 = tmp_25 + tmp_63 + tmp_65;
      real_t tmp_67 = tmp_41 + tmp_56;
      real_t tmp_68 = p_affine_0_0*tmp_13;
      real_t tmp_69 = p_affine_0_1*tmp_13;
      real_t tmp_70 = p_affine_0_2*tmp_13;
      real_t tmp_71 = (p_affine_2_0*p_affine_2_0)*tmp_17;
      real_t tmp_72 = (p_affine_2_1*p_affine_2_1)*tmp_17;
      real_t tmp_73 = (p_affine_2_2*p_affine_2_2)*tmp_17;
      real_t tmp_74 = tmp_71 + tmp_72 + tmp_73;
      real_t tmp_75 = p_affine_3_0*tmp_30;
      real_t tmp_76 = p_affine_3_1*tmp_32;
      real_t tmp_77 = p_affine_3_2*tmp_34;
      real_t tmp_78 = -tmp_75 - tmp_76 - tmp_77;
      real_t tmp_79 = tmp_64 + tmp_78;
      real_t tmp_80 = tmp_36 + tmp_42 + tmp_79;
      real_t tmp_81 = tmp_60 - tmp_71 - tmp_72 - tmp_73;
      real_t tmp_82 = tmp_36 + tmp_44 + tmp_81;
      real_t tmp_83 = tmp_75 + tmp_76 + tmp_77;
      real_t tmp_84 = tmp_40 + tmp_83;
      real_t tmp_85 = tmp_29 + tmp_51 + tmp_84;
      real_t tmp_86 = tmp_58 + tmp_74 + tmp_79;
      real_t tmp_87 = (p_affine_3_0*p_affine_3_0)*tmp_17;
      real_t tmp_88 = (p_affine_3_1*p_affine_3_1)*tmp_17;
      real_t tmp_89 = (p_affine_3_2*p_affine_3_2)*tmp_17;
      real_t tmp_90 = tmp_87 + tmp_88 + tmp_89;
      real_t tmp_91 = tmp_63 + tmp_83;
      real_t tmp_92 = tmp_59 + tmp_91;
      real_t tmp_93 = tmp_55 - tmp_87 - tmp_88 - tmp_89;
      real_t tmp_94 = tmp_65 + tmp_93;
      real_t tmp_95 = tmp_36 + tmp_51 + tmp_78 + tmp_90;
      real_t tmp_96 = p_affine_1_0*tmp_13;
      real_t tmp_97 = p_affine_1_1*tmp_13;
      real_t tmp_98 = p_affine_1_2*tmp_13;
      real_t tmp_99 = tmp_55 + tmp_61 + tmp_78;
      real_t tmp_100 = tmp_81 + tmp_91;
      real_t tmp_101 = tmp_84 + tmp_93;
      real_t a_0_0 = -tmp_13*tmp_14 - tmp_13*tmp_15 - tmp_13*tmp_16 + tmp_21 + tmp_25;
      real_t a_0_1 = tmp_43;
      real_t a_0_2 = tmp_57;
      real_t a_0_3 = tmp_62;
      real_t a_0_4 = tmp_66;
      real_t a_0_5 = tmp_67;
      real_t a_1_0 = tmp_43;
      real_t a_1_1 = -p_affine_2_0*tmp_68 - p_affine_2_1*tmp_69 - p_affine_2_2*tmp_70 + tmp_21 + tmp_74;
      real_t a_1_2 = tmp_80;
      real_t a_1_3 = tmp_82;
      real_t a_1_4 = tmp_85;
      real_t a_1_5 = tmp_86;
      real_t a_2_0 = tmp_57;
      real_t a_2_1 = tmp_80;
      real_t a_2_2 = -p_affine_3_0*tmp_68 - p_affine_3_1*tmp_69 - p_affine_3_2*tmp_70 + tmp_21 + tmp_90;
      real_t a_2_3 = tmp_92;
      real_t a_2_4 = tmp_94;
      real_t a_2_5 = tmp_95;
      real_t a_3_0 = tmp_62;
      real_t a_3_1 = tmp_82;
      real_t a_3_2 = tmp_92;
      real_t a_3_3 = -p_affine_2_0*tmp_96 - p_affine_2_1*tmp_97 - p_affine_2_2*tmp_98 + tmp_25 + tmp_74;
      real_t a_3_4 = tmp_99;
      real_t a_3_5 = tmp_100;
      real_t a_4_0 = tmp_66;
      real_t a_4_1 = tmp_85;
      real_t a_4_2 = tmp_94;
      real_t a_4_3 = tmp_99;
      real_t a_4_4 = -p_affine_3_0*tmp_96 - p_affine_3_1*tmp_97 - p_affine_3_2*tmp_98 + tmp_25 + tmp_90;
      real_t a_4_5 = tmp_101;
      real_t a_5_0 = tmp_67;
      real_t a_5_1 = tmp_86;
      real_t a_5_2 = tmp_95;
      real_t a_5_3 = tmp_100;
      real_t a_5_4 = tmp_101;
      real_t a_5_5 = -p_affine_2_0*p_affine_3_0*tmp_13 - p_affine_2_1*p_affine_3_1*tmp_13 - p_affine_2_2*p_affine_3_2*tmp_13 + tmp_74 + tmp_90;
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

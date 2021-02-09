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

#include "p1_epsiloncc_1_2_affine_q2.hpp"

namespace hyteg {
namespace forms {

   void p1_epsiloncc_1_2_affine_q2::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t q_p_0_0 = 0.16666666666666666;
      real_t q_p_0_1 = 0.66666666666666663;
      real_t q_p_1_0 = 0.66666666666666663;
      real_t q_p_1_1 = 0.16666666666666666;
      real_t q_p_2_0 = 0.16666666666666666;
      real_t q_p_2_1 = 0.16666666666666666;
      real_t w_p_0 = 0.16666666666666666;
      real_t w_p_1 = 0.16666666666666666;
      real_t w_p_2 = 0.16666666666666666;
      real_t a_0_0 = 0;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      real_t a_1_0 = 0;
      real_t a_1_1 = 0;
      real_t a_1_2 = 0;
      real_t a_2_0 = 0;
      real_t a_2_1 = 0;
      real_t a_2_2 = 0;
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

   void p1_epsiloncc_1_2_affine_q2::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t q_p_0_0 = 0.13819660112501059;
      real_t q_p_0_1 = 0.13819660112501059;
      real_t q_p_0_2 = 0.58541019662496829;
      real_t q_p_1_0 = 0.13819660112501059;
      real_t q_p_1_1 = 0.58541019662496829;
      real_t q_p_1_2 = 0.13819660112501059;
      real_t q_p_2_0 = 0.58541019662496829;
      real_t q_p_2_1 = 0.13819660112501059;
      real_t q_p_2_2 = 0.13819660112501059;
      real_t q_p_3_0 = 0.13819660112501059;
      real_t q_p_3_1 = 0.13819660112501059;
      real_t q_p_3_2 = 0.13819660112501059;
      real_t w_p_0 = 0.041666666666666657;
      real_t w_p_1 = 0.041666666666666657;
      real_t w_p_2 = 0.041666666666666657;
      real_t w_p_3 = 0.041666666666666657;
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = p_affine_1_2 + tmp_9;
      real_t tmp_12 = p_affine_3_1 + tmp_2;
      real_t tmp_13 = tmp_12*tmp_5;
      real_t tmp_14 = p_affine_2_2 + tmp_9;
      real_t tmp_15 = p_affine_3_0 + tmp_0;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t tmp_17 = tmp_1*tmp_12;
      real_t tmp_18 = tmp_15*tmp_3;
      real_t tmp_19 = tmp_10*tmp_4 - tmp_10*tmp_7 + tmp_11*tmp_13 - tmp_11*tmp_18 + tmp_14*tmp_16 - tmp_14*tmp_17;
      real_t tmp_20 = 1.0 / (tmp_19);
      real_t tmp_21 = 0.5*tmp_20;
      real_t tmp_22 = tmp_16 - tmp_17;
      real_t tmp_23 = tmp_13 - tmp_18;
      real_t tmp_24 = p_affine_0_0*p_affine_1_1;
      real_t tmp_25 = p_affine_0_0*p_affine_1_2;
      real_t tmp_26 = p_affine_2_1*p_affine_3_2;
      real_t tmp_27 = p_affine_0_1*p_affine_1_0;
      real_t tmp_28 = p_affine_0_1*p_affine_1_2;
      real_t tmp_29 = p_affine_2_2*p_affine_3_0;
      real_t tmp_30 = p_affine_0_2*p_affine_1_0;
      real_t tmp_31 = p_affine_0_2*p_affine_1_1;
      real_t tmp_32 = p_affine_2_0*p_affine_3_1;
      real_t tmp_33 = p_affine_2_2*p_affine_3_1;
      real_t tmp_34 = p_affine_2_0*p_affine_3_2;
      real_t tmp_35 = p_affine_2_1*p_affine_3_0;
      real_t tmp_36 = std::abs(p_affine_0_0*tmp_26 - p_affine_0_0*tmp_33 + p_affine_0_1*tmp_29 - p_affine_0_1*tmp_34 + p_affine_0_2*tmp_32 - p_affine_0_2*tmp_35 - p_affine_1_0*tmp_26 + p_affine_1_0*tmp_33 - p_affine_1_1*tmp_29 + p_affine_1_1*tmp_34 - p_affine_1_2*tmp_32 + p_affine_1_2*tmp_35 + p_affine_2_0*tmp_28 - p_affine_2_0*tmp_31 - p_affine_2_1*tmp_25 + p_affine_2_1*tmp_30 + p_affine_2_2*tmp_24 - p_affine_2_2*tmp_27 - p_affine_3_0*tmp_28 + p_affine_3_0*tmp_31 + p_affine_3_1*tmp_25 - p_affine_3_1*tmp_30 - p_affine_3_2*tmp_24 + p_affine_3_2*tmp_27);
      real_t tmp_37 = tmp_36*(-tmp_21*tmp_22 - tmp_21*tmp_23 - tmp_21*tmp_8);
      real_t tmp_38 = tmp_37*w_p_0;
      real_t tmp_39 = -tmp_1*tmp_14 + tmp_11*tmp_5;
      real_t tmp_40 = 1.0*tmp_20;
      real_t tmp_41 = tmp_1*tmp_10 - tmp_11*tmp_15;
      real_t tmp_42 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_43 = -tmp_39*tmp_40 - tmp_40*tmp_41 - tmp_40*tmp_42;
      real_t tmp_44 = 2*tmp_43;
      real_t tmp_45 = tmp_37*tmp_44;
      real_t tmp_46 = 2.0*tmp_20;
      real_t tmp_47 = tmp_42*tmp_46;
      real_t tmp_48 = tmp_37*tmp_47;
      real_t tmp_49 = tmp_41*tmp_46;
      real_t tmp_50 = tmp_37*tmp_49;
      real_t tmp_51 = tmp_39*tmp_46;
      real_t tmp_52 = tmp_37*tmp_51;
      real_t tmp_53 = tmp_23*tmp_36;
      real_t tmp_54 = tmp_53*w_p_0;
      real_t tmp_55 = tmp_40*tmp_43;
      real_t tmp_56 = tmp_53*tmp_55;
      real_t tmp_57 = 1.0/(tmp_19*tmp_19);
      real_t tmp_58 = tmp_42*tmp_57;
      real_t tmp_59 = tmp_53*tmp_58;
      real_t tmp_60 = tmp_41*tmp_57;
      real_t tmp_61 = tmp_53*tmp_60;
      real_t tmp_62 = tmp_39*tmp_57;
      real_t tmp_63 = tmp_53*tmp_62;
      real_t tmp_64 = tmp_22*tmp_36;
      real_t tmp_65 = tmp_55*tmp_64;
      real_t tmp_66 = tmp_58*tmp_64;
      real_t tmp_67 = tmp_60*tmp_64;
      real_t tmp_68 = tmp_62*tmp_64;
      real_t tmp_69 = tmp_36*tmp_8;
      real_t tmp_70 = tmp_55*tmp_69;
      real_t tmp_71 = tmp_58*tmp_69;
      real_t tmp_72 = tmp_60*tmp_69;
      real_t tmp_73 = tmp_62*tmp_69;
      real_t a_0_0 = tmp_38*tmp_44 + tmp_45*w_p_1 + tmp_45*w_p_2 + tmp_45*w_p_3;
      real_t a_0_1 = tmp_38*tmp_47 + tmp_48*w_p_1 + tmp_48*w_p_2 + tmp_48*w_p_3;
      real_t a_0_2 = tmp_38*tmp_49 + tmp_50*w_p_1 + tmp_50*w_p_2 + tmp_50*w_p_3;
      real_t a_0_3 = tmp_38*tmp_51 + tmp_52*w_p_1 + tmp_52*w_p_2 + tmp_52*w_p_3;
      real_t a_1_0 = tmp_54*tmp_55 + tmp_56*w_p_1 + tmp_56*w_p_2 + tmp_56*w_p_3;
      real_t a_1_1 = tmp_54*tmp_58 + tmp_59*w_p_1 + tmp_59*w_p_2 + tmp_59*w_p_3;
      real_t a_1_2 = tmp_54*tmp_60 + tmp_61*w_p_1 + tmp_61*w_p_2 + tmp_61*w_p_3;
      real_t a_1_3 = tmp_54*tmp_62 + tmp_63*w_p_1 + tmp_63*w_p_2 + tmp_63*w_p_3;
      real_t a_2_0 = tmp_65*w_p_0 + tmp_65*w_p_1 + tmp_65*w_p_2 + tmp_65*w_p_3;
      real_t a_2_1 = tmp_66*w_p_0 + tmp_66*w_p_1 + tmp_66*w_p_2 + tmp_66*w_p_3;
      real_t a_2_2 = tmp_67*w_p_0 + tmp_67*w_p_1 + tmp_67*w_p_2 + tmp_67*w_p_3;
      real_t a_2_3 = tmp_68*w_p_0 + tmp_68*w_p_1 + tmp_68*w_p_2 + tmp_68*w_p_3;
      real_t a_3_0 = tmp_70*w_p_0 + tmp_70*w_p_1 + tmp_70*w_p_2 + tmp_70*w_p_3;
      real_t a_3_1 = tmp_71*w_p_0 + tmp_71*w_p_1 + tmp_71*w_p_2 + tmp_71*w_p_3;
      real_t a_3_2 = tmp_72*w_p_0 + tmp_72*w_p_1 + tmp_72*w_p_2 + tmp_72*w_p_3;
      real_t a_3_3 = tmp_73*w_p_0 + tmp_73*w_p_1 + tmp_73*w_p_2 + tmp_73*w_p_3;
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

} // namespace forms
} // namespace hyteg

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

#include "p1_epsiloncc_1_0_affine_q2.hpp"

namespace hyteg {
namespace forms {

   void p1_epsiloncc_1_0_affine_q2::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_0;
      real_t tmp_3 = p_affine_1_0 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = tmp_4 - (p_affine_1_1 + tmp_0)*(p_affine_2_0 + tmp_2);
      real_t tmp_6 = 1.0 / (tmp_5);
      real_t tmp_7 = 0.5*tmp_6;
      real_t tmp_8 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_9 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_10 = tmp_9*(-tmp_1*tmp_7 - tmp_7*tmp_8);
      real_t tmp_11 = 0.16666666666666666*tmp_10;
      real_t tmp_12 = 1.0*tmp_6;
      real_t tmp_13 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_14 = -tmp_12*tmp_13 - tmp_12*tmp_3;
      real_t tmp_15 = 2*tmp_14;
      real_t tmp_16 = tmp_10*tmp_15;
      real_t tmp_17 = 2.0*tmp_6;
      real_t tmp_18 = tmp_13*tmp_17;
      real_t tmp_19 = tmp_10*tmp_18;
      real_t tmp_20 = tmp_17*tmp_3;
      real_t tmp_21 = tmp_10*tmp_20;
      real_t tmp_22 = tmp_1*tmp_9;
      real_t tmp_23 = 0.16666666666666666*tmp_22;
      real_t tmp_24 = tmp_12*tmp_14;
      real_t tmp_25 = tmp_22*tmp_24;
      real_t tmp_26 = 1.0/(tmp_5*tmp_5);
      real_t tmp_27 = tmp_13*tmp_26;
      real_t tmp_28 = tmp_22*tmp_27;
      real_t tmp_29 = 0.16666666666666666*tmp_9;
      real_t tmp_30 = tmp_26*tmp_4;
      real_t tmp_31 = tmp_30*tmp_9;
      real_t tmp_32 = tmp_24*tmp_8;
      real_t tmp_33 = tmp_32*tmp_9;
      real_t tmp_34 = tmp_27*tmp_8;
      real_t tmp_35 = tmp_34*tmp_9;
      real_t tmp_36 = tmp_26*tmp_3*tmp_8;
      real_t tmp_37 = tmp_36*tmp_9;
      real_t a_0_0 = tmp_11*tmp_15 + 0.33333333333333331*tmp_16;
      real_t a_0_1 = tmp_11*tmp_18 + 0.33333333333333331*tmp_19;
      real_t a_0_2 = tmp_11*tmp_20 + 0.33333333333333331*tmp_21;
      real_t a_1_0 = tmp_23*tmp_24 + 0.33333333333333331*tmp_25;
      real_t a_1_1 = tmp_23*tmp_27 + 0.33333333333333331*tmp_28;
      real_t a_1_2 = tmp_29*tmp_30 + 0.33333333333333331*tmp_31;
      real_t a_2_0 = tmp_29*tmp_32 + 0.33333333333333331*tmp_33;
      real_t a_2_1 = tmp_29*tmp_34 + 0.33333333333333331*tmp_35;
      real_t a_2_2 = tmp_29*tmp_36 + 0.33333333333333331*tmp_37;
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

   void p1_epsiloncc_1_0_affine_q2::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_0;
      real_t tmp_3 = p_affine_1_0 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_0)*(p_affine_2_0 + tmp_2));
      real_t tmp_5 = 0.5*tmp_4;
      real_t tmp_6 = (-tmp_1*tmp_5 - tmp_5*(p_affine_0_1 - p_affine_1_1))*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_7 = 0.16666666666666666*tmp_6;
      real_t tmp_8 = 1.0*tmp_4;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = -2*tmp_3*tmp_8 - 2*tmp_8*tmp_9;
      real_t tmp_11 = tmp_10*tmp_6;
      real_t tmp_12 = 2.0*tmp_4;
      real_t tmp_13 = tmp_12*tmp_9;
      real_t tmp_14 = tmp_13*tmp_6;
      real_t tmp_15 = tmp_12*tmp_3;
      real_t tmp_16 = tmp_15*tmp_6;
      real_t a_0_0 = tmp_10*tmp_7 + 0.33333333333333331*tmp_11;
      real_t a_0_1 = tmp_13*tmp_7 + 0.33333333333333331*tmp_14;
      real_t a_0_2 = tmp_15*tmp_7 + 0.33333333333333331*tmp_16;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_epsiloncc_1_0_affine_q2::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_1_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_2;
      real_t tmp_3 = p_affine_2_2 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_1 + tmp_0;
      real_t tmp_6 = p_affine_1_2 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = -p_affine_0_0;
      real_t tmp_10 = p_affine_1_0 + tmp_9;
      real_t tmp_11 = p_affine_3_2 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = p_affine_2_0 + tmp_9;
      real_t tmp_14 = p_affine_3_1 + tmp_0;
      real_t tmp_15 = tmp_14*tmp_6;
      real_t tmp_16 = p_affine_3_0 + tmp_9;
      real_t tmp_17 = tmp_14*tmp_3;
      real_t tmp_18 = tmp_1*tmp_11;
      real_t tmp_19 = tmp_10*tmp_12 - tmp_10*tmp_17 + tmp_13*tmp_15 - tmp_13*tmp_18 + tmp_16*tmp_4 - tmp_16*tmp_7;
      real_t tmp_20 = 1.0 / (tmp_19);
      real_t tmp_21 = 0.5*tmp_20;
      real_t tmp_22 = tmp_15 - tmp_18;
      real_t tmp_23 = tmp_12 - tmp_17;
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
      real_t tmp_38 = 0.041666666666666657*tmp_37;
      real_t tmp_39 = -tmp_10*tmp_3 + tmp_13*tmp_6;
      real_t tmp_40 = 1.0*tmp_20;
      real_t tmp_41 = tmp_10*tmp_11 - tmp_16*tmp_6;
      real_t tmp_42 = -tmp_11*tmp_13 + tmp_16*tmp_3;
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
      real_t tmp_54 = 0.041666666666666657*tmp_53;
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
      real_t a_0_0 = tmp_38*tmp_44 + 0.12499999999999997*tmp_45;
      real_t a_0_1 = tmp_38*tmp_47 + 0.12499999999999997*tmp_48;
      real_t a_0_2 = tmp_38*tmp_49 + 0.12499999999999997*tmp_50;
      real_t a_0_3 = tmp_38*tmp_51 + 0.12499999999999997*tmp_52;
      real_t a_1_0 = tmp_54*tmp_55 + 0.12499999999999997*tmp_56;
      real_t a_1_1 = tmp_54*tmp_58 + 0.12499999999999997*tmp_59;
      real_t a_1_2 = tmp_54*tmp_60 + 0.12499999999999997*tmp_61;
      real_t a_1_3 = tmp_54*tmp_62 + 0.12499999999999997*tmp_63;
      real_t a_2_0 = 0.16666666666666663*tmp_65;
      real_t a_2_1 = 0.16666666666666663*tmp_66;
      real_t a_2_2 = 0.16666666666666663*tmp_67;
      real_t a_2_3 = 0.16666666666666663*tmp_68;
      real_t a_3_0 = 0.16666666666666663*tmp_70;
      real_t a_3_1 = 0.16666666666666663*tmp_71;
      real_t a_3_2 = 0.16666666666666663*tmp_72;
      real_t a_3_3 = 0.16666666666666663*tmp_73;
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

   void p1_epsiloncc_1_0_affine_q2::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_1_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_2;
      real_t tmp_3 = p_affine_2_2 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_1 + tmp_0;
      real_t tmp_6 = p_affine_1_2 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = -p_affine_0_0;
      real_t tmp_9 = p_affine_1_0 + tmp_8;
      real_t tmp_10 = p_affine_3_2 + tmp_2;
      real_t tmp_11 = tmp_10*tmp_5;
      real_t tmp_12 = p_affine_2_0 + tmp_8;
      real_t tmp_13 = p_affine_3_1 + tmp_0;
      real_t tmp_14 = tmp_13*tmp_6;
      real_t tmp_15 = p_affine_3_0 + tmp_8;
      real_t tmp_16 = tmp_13*tmp_3;
      real_t tmp_17 = tmp_1*tmp_10;
      real_t tmp_18 = 1.0 / (tmp_11*tmp_9 + tmp_12*tmp_14 - tmp_12*tmp_17 + tmp_15*tmp_4 - tmp_15*tmp_7 - tmp_16*tmp_9);
      real_t tmp_19 = 0.5*tmp_18;
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
      real_t tmp_32 = (-tmp_19*(tmp_11 - tmp_16) - tmp_19*(tmp_14 - tmp_17) - tmp_19*(tmp_4 - tmp_7))*std::abs(p_affine_0_0*tmp_22 - p_affine_0_0*tmp_29 + p_affine_0_1*tmp_25 - p_affine_0_1*tmp_30 + p_affine_0_2*tmp_28 - p_affine_0_2*tmp_31 - p_affine_1_0*tmp_22 + p_affine_1_0*tmp_29 - p_affine_1_1*tmp_25 + p_affine_1_1*tmp_30 - p_affine_1_2*tmp_28 + p_affine_1_2*tmp_31 + p_affine_2_0*tmp_24 - p_affine_2_0*tmp_27 - p_affine_2_1*tmp_21 + p_affine_2_1*tmp_26 + p_affine_2_2*tmp_20 - p_affine_2_2*tmp_23 - p_affine_3_0*tmp_24 + p_affine_3_0*tmp_27 + p_affine_3_1*tmp_21 - p_affine_3_1*tmp_26 - p_affine_3_2*tmp_20 + p_affine_3_2*tmp_23);
      real_t tmp_33 = 0.041666666666666657*tmp_32;
      real_t tmp_34 = tmp_12*tmp_6 - tmp_3*tmp_9;
      real_t tmp_35 = 1.0*tmp_18;
      real_t tmp_36 = tmp_10*tmp_9 - tmp_15*tmp_6;
      real_t tmp_37 = -tmp_10*tmp_12 + tmp_15*tmp_3;
      real_t tmp_38 = -2*tmp_34*tmp_35 - 2*tmp_35*tmp_36 - 2*tmp_35*tmp_37;
      real_t tmp_39 = tmp_32*tmp_38;
      real_t tmp_40 = 2.0*tmp_18;
      real_t tmp_41 = tmp_37*tmp_40;
      real_t tmp_42 = tmp_32*tmp_41;
      real_t tmp_43 = tmp_36*tmp_40;
      real_t tmp_44 = tmp_32*tmp_43;
      real_t tmp_45 = tmp_34*tmp_40;
      real_t tmp_46 = tmp_32*tmp_45;
      real_t a_0_0 = tmp_33*tmp_38 + 0.12499999999999997*tmp_39;
      real_t a_0_1 = tmp_33*tmp_41 + 0.12499999999999997*tmp_42;
      real_t a_0_2 = tmp_33*tmp_43 + 0.12499999999999997*tmp_44;
      real_t a_0_3 = tmp_33*tmp_45 + 0.12499999999999997*tmp_46;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

} // namespace forms
} // namespace hyteg

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

#include "p1_epsiloncc_0_0_affine_q2.hpp"

namespace hyteg {
namespace forms {

   void p1_epsiloncc_0_0_affine_q2::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
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
      real_t tmp_4 = tmp_1*tmp_3 - (p_affine_1_1 + tmp_0)*(p_affine_2_0 + tmp_2);
      real_t tmp_5 = 1.0 / (tmp_4);
      real_t tmp_6 = 2.0*tmp_5;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_9 = tmp_6*tmp_8;
      real_t tmp_10 = -tmp_7 - tmp_9;
      real_t tmp_11 = 1.0*tmp_5;
      real_t tmp_12 = tmp_1*tmp_11;
      real_t tmp_13 = tmp_11*tmp_8;
      real_t tmp_14 = -tmp_12 - tmp_13;
      real_t tmp_15 = tmp_11*tmp_3;
      real_t tmp_16 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_17 = tmp_11*tmp_16;
      real_t tmp_18 = -tmp_15 - tmp_17;
      real_t tmp_19 = 0.5*tmp_5;
      real_t tmp_20 = -tmp_16*tmp_19 - tmp_19*tmp_3;
      real_t tmp_21 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_22 = tmp_21*(tmp_10*tmp_14 + 2*tmp_18*tmp_20);
      real_t tmp_23 = tmp_20*tmp_6;
      real_t tmp_24 = tmp_21*(tmp_14*tmp_7 + tmp_16*tmp_23);
      real_t tmp_25 = tmp_21*(tmp_14*tmp_9 + tmp_23*tmp_3);
      real_t tmp_26 = tmp_21*(tmp_10*tmp_12 + tmp_17*tmp_18);
      real_t tmp_27 = 1.0 / (tmp_4*tmp_4);
      real_t tmp_28 = 1.0*tmp_27;
      real_t tmp_29 = 2.0*tmp_27;
      real_t tmp_30 = tmp_21*((tmp_1*tmp_1)*tmp_29 + (tmp_16*tmp_16)*tmp_28);
      real_t tmp_31 = tmp_21*(tmp_1*tmp_29*tmp_8 + tmp_16*tmp_28*tmp_3);
      real_t tmp_32 = 0.5*tmp_31;
      real_t tmp_33 = tmp_21*(tmp_10*tmp_13 + tmp_15*tmp_18);
      real_t tmp_34 = tmp_21*(tmp_28*(tmp_3*tmp_3) + tmp_29*(tmp_8*tmp_8));
      real_t a_0_0 = 0.5*tmp_22;
      real_t a_0_1 = 0.5*tmp_24;
      real_t a_0_2 = 0.5*tmp_25;
      real_t a_1_0 = 0.5*tmp_26;
      real_t a_1_1 = 0.5*tmp_30;
      real_t a_1_2 = tmp_32;
      real_t a_2_0 = 0.5*tmp_33;
      real_t a_2_1 = tmp_32;
      real_t a_2_2 = 0.5*tmp_34;
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

   void p1_epsiloncc_0_0_affine_q2::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
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
      real_t tmp_5 = 2.0*tmp_4;
      real_t tmp_6 = tmp_1*tmp_5;
      real_t tmp_7 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_8 = tmp_5*tmp_7;
      real_t tmp_9 = 1.0*tmp_4;
      real_t tmp_10 = -tmp_1*tmp_9 - tmp_7*tmp_9;
      real_t tmp_11 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_12 = 0.5*tmp_4;
      real_t tmp_13 = -tmp_11*tmp_12 - tmp_12*tmp_3;
      real_t tmp_14 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_15 = tmp_14*(tmp_10*(-tmp_6 - tmp_8) + 2*tmp_13*(-tmp_11*tmp_9 - tmp_3*tmp_9));
      real_t tmp_16 = tmp_13*tmp_5;
      real_t tmp_17 = tmp_14*(tmp_10*tmp_6 + tmp_11*tmp_16);
      real_t tmp_18 = tmp_14*(tmp_10*tmp_8 + tmp_16*tmp_3);
      real_t a_0_0 = 0.5*tmp_15;
      real_t a_0_1 = 0.5*tmp_17;
      real_t a_0_2 = 0.5*tmp_18;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_epsiloncc_0_0_affine_q2::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t tmp_21 = 2.0*tmp_20;
      real_t tmp_22 = tmp_21*tmp_8;
      real_t tmp_23 = tmp_15 - tmp_18;
      real_t tmp_24 = tmp_21*tmp_23;
      real_t tmp_25 = tmp_12 - tmp_17;
      real_t tmp_26 = tmp_21*tmp_25;
      real_t tmp_27 = -tmp_22 - tmp_24 - tmp_26;
      real_t tmp_28 = 1.0*tmp_20;
      real_t tmp_29 = tmp_28*tmp_8;
      real_t tmp_30 = tmp_23*tmp_28;
      real_t tmp_31 = tmp_25*tmp_28;
      real_t tmp_32 = -tmp_29 - tmp_30 - tmp_31;
      real_t tmp_33 = -tmp_1*tmp_13 + tmp_10*tmp_5;
      real_t tmp_34 = tmp_28*tmp_33;
      real_t tmp_35 = tmp_1*tmp_16 - tmp_10*tmp_14;
      real_t tmp_36 = tmp_28*tmp_35;
      real_t tmp_37 = tmp_13*tmp_14 - tmp_16*tmp_5;
      real_t tmp_38 = tmp_28*tmp_37;
      real_t tmp_39 = -tmp_34 - tmp_36 - tmp_38;
      real_t tmp_40 = 0.5*tmp_20;
      real_t tmp_41 = -tmp_33*tmp_40 - tmp_35*tmp_40 - tmp_37*tmp_40;
      real_t tmp_42 = -tmp_10*tmp_3 + tmp_13*tmp_6;
      real_t tmp_43 = tmp_28*tmp_42;
      real_t tmp_44 = tmp_10*tmp_11 - tmp_16*tmp_6;
      real_t tmp_45 = tmp_28*tmp_44;
      real_t tmp_46 = -tmp_11*tmp_13 + tmp_16*tmp_3;
      real_t tmp_47 = tmp_28*tmp_46;
      real_t tmp_48 = -tmp_43 - tmp_45 - tmp_47;
      real_t tmp_49 = -tmp_40*tmp_42 - tmp_40*tmp_44 - tmp_40*tmp_46;
      real_t tmp_50 = p_affine_0_0*p_affine_1_1;
      real_t tmp_51 = p_affine_0_0*p_affine_1_2;
      real_t tmp_52 = p_affine_2_1*p_affine_3_2;
      real_t tmp_53 = p_affine_0_1*p_affine_1_0;
      real_t tmp_54 = p_affine_0_1*p_affine_1_2;
      real_t tmp_55 = p_affine_2_2*p_affine_3_0;
      real_t tmp_56 = p_affine_0_2*p_affine_1_0;
      real_t tmp_57 = p_affine_0_2*p_affine_1_1;
      real_t tmp_58 = p_affine_2_0*p_affine_3_1;
      real_t tmp_59 = p_affine_2_2*p_affine_3_1;
      real_t tmp_60 = p_affine_2_0*p_affine_3_2;
      real_t tmp_61 = p_affine_2_1*p_affine_3_0;
      real_t tmp_62 = std::abs(p_affine_0_0*tmp_52 - p_affine_0_0*tmp_59 + p_affine_0_1*tmp_55 - p_affine_0_1*tmp_60 + p_affine_0_2*tmp_58 - p_affine_0_2*tmp_61 - p_affine_1_0*tmp_52 + p_affine_1_0*tmp_59 - p_affine_1_1*tmp_55 + p_affine_1_1*tmp_60 - p_affine_1_2*tmp_58 + p_affine_1_2*tmp_61 + p_affine_2_0*tmp_54 - p_affine_2_0*tmp_57 - p_affine_2_1*tmp_51 + p_affine_2_1*tmp_56 + p_affine_2_2*tmp_50 - p_affine_2_2*tmp_53 - p_affine_3_0*tmp_54 + p_affine_3_0*tmp_57 + p_affine_3_1*tmp_51 - p_affine_3_1*tmp_56 - p_affine_3_2*tmp_50 + p_affine_3_2*tmp_53);
      real_t tmp_63 = tmp_62*(tmp_27*tmp_32 + 2*tmp_39*tmp_41 + 2*tmp_48*tmp_49);
      real_t tmp_64 = tmp_21*tmp_41;
      real_t tmp_65 = tmp_21*tmp_49;
      real_t tmp_66 = tmp_62*(tmp_26*tmp_32 + tmp_37*tmp_64 + tmp_46*tmp_65);
      real_t tmp_67 = tmp_62*(tmp_24*tmp_32 + tmp_35*tmp_64 + tmp_44*tmp_65);
      real_t tmp_68 = tmp_62*(tmp_22*tmp_32 + tmp_33*tmp_64 + tmp_42*tmp_65);
      real_t tmp_69 = tmp_62*(tmp_27*tmp_31 + tmp_38*tmp_39 + tmp_47*tmp_48);
      real_t tmp_70 = 1.0 / (tmp_19*tmp_19);
      real_t tmp_71 = 1.0*tmp_70;
      real_t tmp_72 = 2.0*tmp_70;
      real_t tmp_73 = tmp_62*((tmp_25*tmp_25)*tmp_72 + (tmp_37*tmp_37)*tmp_71 + (tmp_46*tmp_46)*tmp_71);
      real_t tmp_74 = tmp_37*tmp_71;
      real_t tmp_75 = tmp_46*tmp_71;
      real_t tmp_76 = tmp_25*tmp_72;
      real_t tmp_77 = tmp_62*(tmp_23*tmp_76 + tmp_35*tmp_74 + tmp_44*tmp_75);
      real_t tmp_78 = 0.16666666666666663*tmp_77;
      real_t tmp_79 = tmp_62*(tmp_33*tmp_74 + tmp_42*tmp_75 + tmp_76*tmp_8);
      real_t tmp_80 = 0.16666666666666663*tmp_79;
      real_t tmp_81 = tmp_62*(tmp_27*tmp_30 + tmp_36*tmp_39 + tmp_45*tmp_48);
      real_t tmp_82 = tmp_62*((tmp_23*tmp_23)*tmp_72 + (tmp_35*tmp_35)*tmp_71 + (tmp_44*tmp_44)*tmp_71);
      real_t tmp_83 = tmp_62*(tmp_23*tmp_72*tmp_8 + tmp_33*tmp_35*tmp_71 + tmp_42*tmp_44*tmp_71);
      real_t tmp_84 = 0.16666666666666663*tmp_83;
      real_t tmp_85 = tmp_62*(tmp_27*tmp_29 + tmp_34*tmp_39 + tmp_43*tmp_48);
      real_t tmp_86 = tmp_62*((tmp_33*tmp_33)*tmp_71 + (tmp_42*tmp_42)*tmp_71 + tmp_72*(tmp_8*tmp_8));
      real_t a_0_0 = 0.16666666666666663*tmp_63;
      real_t a_0_1 = 0.16666666666666663*tmp_66;
      real_t a_0_2 = 0.16666666666666663*tmp_67;
      real_t a_0_3 = 0.16666666666666663*tmp_68;
      real_t a_1_0 = 0.16666666666666663*tmp_69;
      real_t a_1_1 = 0.16666666666666663*tmp_73;
      real_t a_1_2 = tmp_78;
      real_t a_1_3 = tmp_80;
      real_t a_2_0 = 0.16666666666666663*tmp_81;
      real_t a_2_1 = tmp_78;
      real_t a_2_2 = 0.16666666666666663*tmp_82;
      real_t a_2_3 = tmp_84;
      real_t a_3_0 = 0.16666666666666663*tmp_85;
      real_t a_3_1 = tmp_80;
      real_t a_3_2 = tmp_84;
      real_t a_3_3 = 0.16666666666666663*tmp_86;
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

   void p1_epsiloncc_0_0_affine_q2::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t tmp_19 = 1.0 / (tmp_10*tmp_12 - tmp_10*tmp_17 + tmp_13*tmp_15 - tmp_13*tmp_18 + tmp_16*tmp_4 - tmp_16*tmp_7);
      real_t tmp_20 = 2.0*tmp_19;
      real_t tmp_21 = tmp_20*tmp_8;
      real_t tmp_22 = tmp_15 - tmp_18;
      real_t tmp_23 = tmp_20*tmp_22;
      real_t tmp_24 = tmp_12 - tmp_17;
      real_t tmp_25 = tmp_20*tmp_24;
      real_t tmp_26 = 1.0*tmp_19;
      real_t tmp_27 = -tmp_22*tmp_26 - tmp_24*tmp_26 - tmp_26*tmp_8;
      real_t tmp_28 = -tmp_1*tmp_13 + tmp_10*tmp_5;
      real_t tmp_29 = tmp_1*tmp_16 - tmp_10*tmp_14;
      real_t tmp_30 = tmp_13*tmp_14 - tmp_16*tmp_5;
      real_t tmp_31 = 0.5*tmp_19;
      real_t tmp_32 = -tmp_28*tmp_31 - tmp_29*tmp_31 - tmp_30*tmp_31;
      real_t tmp_33 = -tmp_10*tmp_3 + tmp_13*tmp_6;
      real_t tmp_34 = tmp_10*tmp_11 - tmp_16*tmp_6;
      real_t tmp_35 = -tmp_11*tmp_13 + tmp_16*tmp_3;
      real_t tmp_36 = -tmp_31*tmp_33 - tmp_31*tmp_34 - tmp_31*tmp_35;
      real_t tmp_37 = p_affine_0_0*p_affine_1_1;
      real_t tmp_38 = p_affine_0_0*p_affine_1_2;
      real_t tmp_39 = p_affine_2_1*p_affine_3_2;
      real_t tmp_40 = p_affine_0_1*p_affine_1_0;
      real_t tmp_41 = p_affine_0_1*p_affine_1_2;
      real_t tmp_42 = p_affine_2_2*p_affine_3_0;
      real_t tmp_43 = p_affine_0_2*p_affine_1_0;
      real_t tmp_44 = p_affine_0_2*p_affine_1_1;
      real_t tmp_45 = p_affine_2_0*p_affine_3_1;
      real_t tmp_46 = p_affine_2_2*p_affine_3_1;
      real_t tmp_47 = p_affine_2_0*p_affine_3_2;
      real_t tmp_48 = p_affine_2_1*p_affine_3_0;
      real_t tmp_49 = std::abs(p_affine_0_0*tmp_39 - p_affine_0_0*tmp_46 + p_affine_0_1*tmp_42 - p_affine_0_1*tmp_47 + p_affine_0_2*tmp_45 - p_affine_0_2*tmp_48 - p_affine_1_0*tmp_39 + p_affine_1_0*tmp_46 - p_affine_1_1*tmp_42 + p_affine_1_1*tmp_47 - p_affine_1_2*tmp_45 + p_affine_1_2*tmp_48 + p_affine_2_0*tmp_41 - p_affine_2_0*tmp_44 - p_affine_2_1*tmp_38 + p_affine_2_1*tmp_43 + p_affine_2_2*tmp_37 - p_affine_2_2*tmp_40 - p_affine_3_0*tmp_41 + p_affine_3_0*tmp_44 + p_affine_3_1*tmp_38 - p_affine_3_1*tmp_43 - p_affine_3_2*tmp_37 + p_affine_3_2*tmp_40);
      real_t tmp_50 = tmp_49*(tmp_27*(-tmp_21 - tmp_23 - tmp_25) + 2*tmp_32*(-tmp_26*tmp_28 - tmp_26*tmp_29 - tmp_26*tmp_30) + 2*tmp_36*(-tmp_26*tmp_33 - tmp_26*tmp_34 - tmp_26*tmp_35));
      real_t tmp_51 = tmp_20*tmp_32;
      real_t tmp_52 = tmp_20*tmp_36;
      real_t tmp_53 = tmp_49*(tmp_25*tmp_27 + tmp_30*tmp_51 + tmp_35*tmp_52);
      real_t tmp_54 = tmp_49*(tmp_23*tmp_27 + tmp_29*tmp_51 + tmp_34*tmp_52);
      real_t tmp_55 = tmp_49*(tmp_21*tmp_27 + tmp_28*tmp_51 + tmp_33*tmp_52);
      real_t a_0_0 = 0.16666666666666663*tmp_50;
      real_t a_0_1 = 0.16666666666666663*tmp_53;
      real_t a_0_2 = 0.16666666666666663*tmp_54;
      real_t a_0_3 = 0.16666666666666663*tmp_55;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

} // namespace forms
} // namespace hyteg

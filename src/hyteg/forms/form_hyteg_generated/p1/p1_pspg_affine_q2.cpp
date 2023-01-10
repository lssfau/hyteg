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

#include "p1_pspg_affine_q2.hpp"

namespace hyteg {
namespace forms {

   void p1_pspg_affine_q2::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_0 = p_affine_0_1*p_affine_1_0;
      real_t tmp_1 = p_affine_0_0*p_affine_1_1;
      real_t tmp_2 = 2.0*tmp_1;
      real_t tmp_3 = p_affine_1_0*p_affine_2_1;
      real_t tmp_4 = 2.0*p_affine_0_0;
      real_t tmp_5 = p_affine_0_1*p_affine_2_0;
      real_t tmp_6 = p_affine_0_0*p_affine_2_1;
      real_t tmp_7 = 2.0*tmp_5;
      real_t tmp_8 = 2.0*p_affine_2_0;
      real_t tmp_9 = 2.0*p_affine_1_1;
      real_t tmp_10 = p_affine_1_1*p_affine_2_0;
      real_t tmp_11 = (p_affine_2_1*p_affine_2_1);
      real_t tmp_12 = (p_affine_1_1*p_affine_1_1);
      real_t tmp_13 = (p_affine_0_0*p_affine_0_0);
      real_t tmp_14 = (p_affine_1_0*p_affine_1_0);
      real_t tmp_15 = 2.0*p_affine_0_1;
      real_t tmp_16 = (p_affine_2_0*p_affine_2_0);
      real_t tmp_17 = (p_affine_0_1*p_affine_0_1);
      real_t tmp_18 = 1.0*tmp_13;
      real_t tmp_19 = 1.0*tmp_17;
      real_t tmp_20 = (std::abs(tmp_0 - tmp_1 + tmp_10 - tmp_3 - tmp_5 + tmp_6)*std::abs(tmp_0 - tmp_1 + tmp_10 - tmp_3 - tmp_5 + tmp_6))/(p_affine_0_1*tmp_3*tmp_4 - p_affine_1_0*tmp_11*tmp_4 - p_affine_1_0*tmp_17*tmp_8 + p_affine_1_0*tmp_5*tmp_9 - p_affine_1_1*tmp_15*tmp_16 - p_affine_2_0*tmp_12*tmp_4 + p_affine_2_1*tmp_1*tmp_8 - p_affine_2_1*tmp_13*tmp_9 - p_affine_2_1*tmp_14*tmp_15 - tmp_0*tmp_2 - 2.0*tmp_10*tmp_3 + 1.0*tmp_11*tmp_14 + tmp_11*tmp_18 + 1.0*tmp_12*tmp_16 + tmp_12*tmp_18 + tmp_14*tmp_19 + tmp_16*tmp_19 + tmp_2*tmp_3 + tmp_2*tmp_5 + tmp_3*tmp_7 - tmp_6*tmp_7);
      real_t tmp_21 = 0.10000000000000001*tmp_20;
      real_t tmp_22 = p_affine_1_0*p_affine_2_0;
      real_t tmp_23 = p_affine_1_1*p_affine_2_1;
      real_t tmp_24 = 0.050000000000000003*tmp_20;
      real_t tmp_25 = tmp_16*tmp_24;
      real_t tmp_26 = tmp_11*tmp_24;
      real_t tmp_27 = -tmp_25 - tmp_26;
      real_t tmp_28 = tmp_14*tmp_24;
      real_t tmp_29 = tmp_12*tmp_24;
      real_t tmp_30 = -tmp_28 - tmp_29;
      real_t tmp_31 = p_affine_0_0*tmp_24;
      real_t tmp_32 = p_affine_1_0*tmp_31;
      real_t tmp_33 = p_affine_0_1*tmp_24;
      real_t tmp_34 = p_affine_1_1*tmp_33;
      real_t tmp_35 = tmp_22*tmp_24;
      real_t tmp_36 = tmp_23*tmp_24;
      real_t tmp_37 = -tmp_35 - tmp_36;
      real_t tmp_38 = p_affine_2_0*tmp_31;
      real_t tmp_39 = p_affine_2_1*tmp_33;
      real_t tmp_40 = -tmp_38 - tmp_39;
      real_t tmp_41 = tmp_25 + tmp_26 + tmp_32 + tmp_34 + tmp_37 + tmp_40;
      real_t tmp_42 = -tmp_32 - tmp_34;
      real_t tmp_43 = tmp_28 + tmp_29 + tmp_37 + tmp_38 + tmp_39 + tmp_42;
      real_t tmp_44 = p_affine_0_0*tmp_21;
      real_t tmp_45 = p_affine_0_1*tmp_21;
      real_t tmp_46 = tmp_13*tmp_24;
      real_t tmp_47 = tmp_17*tmp_24;
      real_t tmp_48 = -tmp_46 - tmp_47;
      real_t tmp_49 = tmp_35 + tmp_36 + tmp_40 + tmp_42 + tmp_46 + tmp_47;
      real_t a_0_0 = tmp_21*tmp_22 + tmp_21*tmp_23 + tmp_27 + tmp_30;
      real_t a_0_1 = tmp_41;
      real_t a_0_2 = tmp_43;
      real_t a_1_0 = tmp_41;
      real_t a_1_1 = p_affine_2_0*tmp_44 + p_affine_2_1*tmp_45 + tmp_27 + tmp_48;
      real_t a_1_2 = tmp_49;
      real_t a_2_0 = tmp_43;
      real_t a_2_1 = tmp_49;
      real_t a_2_2 = p_affine_1_0*tmp_44 + p_affine_1_1*tmp_45 + tmp_30 + tmp_48;
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

   void p1_pspg_affine_q2::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_0 = p_affine_0_1*p_affine_1_0;
      real_t tmp_1 = p_affine_0_0*p_affine_1_1;
      real_t tmp_2 = 2.0*tmp_1;
      real_t tmp_3 = p_affine_1_0*p_affine_2_1;
      real_t tmp_4 = 2.0*p_affine_0_0;
      real_t tmp_5 = p_affine_0_1*p_affine_2_0;
      real_t tmp_6 = p_affine_0_0*p_affine_2_1;
      real_t tmp_7 = 2.0*tmp_5;
      real_t tmp_8 = 2.0*p_affine_2_0;
      real_t tmp_9 = 2.0*p_affine_1_1;
      real_t tmp_10 = p_affine_1_1*p_affine_2_0;
      real_t tmp_11 = (p_affine_2_1*p_affine_2_1);
      real_t tmp_12 = (p_affine_1_1*p_affine_1_1);
      real_t tmp_13 = (p_affine_0_0*p_affine_0_0);
      real_t tmp_14 = (p_affine_1_0*p_affine_1_0);
      real_t tmp_15 = 2.0*p_affine_0_1;
      real_t tmp_16 = (p_affine_2_0*p_affine_2_0);
      real_t tmp_17 = (p_affine_0_1*p_affine_0_1);
      real_t tmp_18 = 1.0*tmp_13;
      real_t tmp_19 = 1.0*tmp_17;
      real_t tmp_20 = (std::abs(tmp_0 - tmp_1 + tmp_10 - tmp_3 - tmp_5 + tmp_6)*std::abs(tmp_0 - tmp_1 + tmp_10 - tmp_3 - tmp_5 + tmp_6))/(p_affine_0_1*tmp_3*tmp_4 - p_affine_1_0*tmp_11*tmp_4 - p_affine_1_0*tmp_17*tmp_8 + p_affine_1_0*tmp_5*tmp_9 - p_affine_1_1*tmp_15*tmp_16 - p_affine_2_0*tmp_12*tmp_4 + p_affine_2_1*tmp_1*tmp_8 - p_affine_2_1*tmp_13*tmp_9 - p_affine_2_1*tmp_14*tmp_15 - tmp_0*tmp_2 - 2.0*tmp_10*tmp_3 + 1.0*tmp_11*tmp_14 + tmp_11*tmp_18 + 1.0*tmp_12*tmp_16 + tmp_12*tmp_18 + tmp_14*tmp_19 + tmp_16*tmp_19 + tmp_2*tmp_3 + tmp_2*tmp_5 + tmp_3*tmp_7 - tmp_6*tmp_7);
      real_t tmp_21 = 0.10000000000000001*tmp_20;
      real_t tmp_22 = p_affine_1_0*p_affine_2_0;
      real_t tmp_23 = p_affine_1_1*p_affine_2_1;
      real_t tmp_24 = 0.050000000000000003*tmp_20;
      real_t tmp_25 = tmp_14*tmp_24;
      real_t tmp_26 = tmp_12*tmp_24;
      real_t tmp_27 = tmp_16*tmp_24;
      real_t tmp_28 = tmp_11*tmp_24;
      real_t tmp_29 = p_affine_0_0*tmp_24;
      real_t tmp_30 = p_affine_1_0*tmp_29;
      real_t tmp_31 = p_affine_0_1*tmp_24;
      real_t tmp_32 = p_affine_1_1*tmp_31;
      real_t tmp_33 = p_affine_2_0*tmp_29;
      real_t tmp_34 = p_affine_2_1*tmp_31;
      real_t tmp_35 = -tmp_22*tmp_24 - tmp_23*tmp_24;
      real_t a_0_0 = tmp_21*tmp_22 + tmp_21*tmp_23 - tmp_25 - tmp_26 - tmp_27 - tmp_28;
      real_t a_0_1 = tmp_27 + tmp_28 + tmp_30 + tmp_32 - tmp_33 - tmp_34 + tmp_35;
      real_t a_0_2 = tmp_25 + tmp_26 - tmp_30 - tmp_32 + tmp_33 + tmp_34 + tmp_35;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_pspg_affine_q2::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t tmp_21 = tmp_20*tmp_8;
      real_t tmp_22 = tmp_16 - tmp_17;
      real_t tmp_23 = tmp_20*tmp_22;
      real_t tmp_24 = tmp_13 - tmp_18;
      real_t tmp_25 = tmp_20*tmp_24;
      real_t tmp_26 = -tmp_21 - tmp_23 - tmp_25;
      real_t tmp_27 = -tmp_1*tmp_14 + tmp_11*tmp_5;
      real_t tmp_28 = tmp_20*tmp_27;
      real_t tmp_29 = tmp_1*tmp_10 - tmp_11*tmp_15;
      real_t tmp_30 = tmp_20*tmp_29;
      real_t tmp_31 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_32 = tmp_20*tmp_31;
      real_t tmp_33 = -tmp_28 - tmp_30 - tmp_32;
      real_t tmp_34 = -tmp_11*tmp_3 + tmp_14*tmp_6;
      real_t tmp_35 = tmp_20*tmp_34;
      real_t tmp_36 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_37 = tmp_20*tmp_36;
      real_t tmp_38 = tmp_10*tmp_3 - tmp_12*tmp_14;
      real_t tmp_39 = tmp_20*tmp_38;
      real_t tmp_40 = -tmp_35 - tmp_37 - tmp_39;
      real_t tmp_41 = p_affine_0_0*p_affine_1_1;
      real_t tmp_42 = p_affine_0_0*p_affine_1_2;
      real_t tmp_43 = p_affine_2_1*p_affine_3_2;
      real_t tmp_44 = p_affine_0_1*p_affine_1_0;
      real_t tmp_45 = p_affine_0_1*p_affine_1_2;
      real_t tmp_46 = p_affine_2_2*p_affine_3_0;
      real_t tmp_47 = p_affine_0_2*p_affine_1_0;
      real_t tmp_48 = p_affine_0_2*p_affine_1_1;
      real_t tmp_49 = p_affine_2_0*p_affine_3_1;
      real_t tmp_50 = p_affine_2_2*p_affine_3_1;
      real_t tmp_51 = p_affine_2_0*p_affine_3_2;
      real_t tmp_52 = p_affine_2_1*p_affine_3_0;
      real_t tmp_53 = 0.025237786011557493*std::pow(std::abs(p_affine_0_0*tmp_43 - p_affine_0_0*tmp_50 + p_affine_0_1*tmp_46 - p_affine_0_1*tmp_51 + p_affine_0_2*tmp_49 - p_affine_0_2*tmp_52 - p_affine_1_0*tmp_43 + p_affine_1_0*tmp_50 - p_affine_1_1*tmp_46 + p_affine_1_1*tmp_51 - p_affine_1_2*tmp_49 + p_affine_1_2*tmp_52 + p_affine_2_0*tmp_45 - p_affine_2_0*tmp_48 - p_affine_2_1*tmp_42 + p_affine_2_1*tmp_47 + p_affine_2_2*tmp_41 - p_affine_2_2*tmp_44 - p_affine_3_0*tmp_45 + p_affine_3_0*tmp_48 + p_affine_3_1*tmp_42 - p_affine_3_1*tmp_47 - p_affine_3_2*tmp_41 + p_affine_3_2*tmp_44), 1.6666666666666665);
      real_t tmp_54 = tmp_53*((tmp_26*tmp_26) + (tmp_33*tmp_33) + (tmp_40*tmp_40));
      real_t tmp_55 = tmp_53*(tmp_25*tmp_26 + tmp_32*tmp_33 + tmp_39*tmp_40);
      real_t tmp_56 = -0.16666666666666663*tmp_55;
      real_t tmp_57 = tmp_53*(tmp_23*tmp_26 + tmp_30*tmp_33 + tmp_37*tmp_40);
      real_t tmp_58 = -0.16666666666666663*tmp_57;
      real_t tmp_59 = tmp_53*(tmp_21*tmp_26 + tmp_28*tmp_33 + tmp_35*tmp_40);
      real_t tmp_60 = -0.16666666666666663*tmp_59;
      real_t tmp_61 = 1.0 / (tmp_19*tmp_19);
      real_t tmp_62 = tmp_53*((tmp_24*tmp_24)*tmp_61 + (tmp_31*tmp_31)*tmp_61 + (tmp_38*tmp_38)*tmp_61);
      real_t tmp_63 = tmp_24*tmp_61;
      real_t tmp_64 = tmp_31*tmp_61;
      real_t tmp_65 = tmp_38*tmp_61;
      real_t tmp_66 = tmp_53*(tmp_22*tmp_63 + tmp_29*tmp_64 + tmp_36*tmp_65);
      real_t tmp_67 = -0.16666666666666663*tmp_66;
      real_t tmp_68 = tmp_53*(tmp_27*tmp_64 + tmp_34*tmp_65 + tmp_63*tmp_8);
      real_t tmp_69 = -0.16666666666666663*tmp_68;
      real_t tmp_70 = tmp_53*((tmp_22*tmp_22)*tmp_61 + (tmp_29*tmp_29)*tmp_61 + (tmp_36*tmp_36)*tmp_61);
      real_t tmp_71 = tmp_53*(tmp_22*tmp_61*tmp_8 + tmp_27*tmp_29*tmp_61 + tmp_34*tmp_36*tmp_61);
      real_t tmp_72 = -0.16666666666666663*tmp_71;
      real_t tmp_73 = tmp_53*((tmp_27*tmp_27)*tmp_61 + (tmp_34*tmp_34)*tmp_61 + tmp_61*(tmp_8*tmp_8));
      real_t a_0_0 = -0.16666666666666663*tmp_54;
      real_t a_0_1 = tmp_56;
      real_t a_0_2 = tmp_58;
      real_t a_0_3 = tmp_60;
      real_t a_1_0 = tmp_56;
      real_t a_1_1 = -0.16666666666666663*tmp_62;
      real_t a_1_2 = tmp_67;
      real_t a_1_3 = tmp_69;
      real_t a_2_0 = tmp_58;
      real_t a_2_1 = tmp_67;
      real_t a_2_2 = -0.16666666666666663*tmp_70;
      real_t a_2_3 = tmp_72;
      real_t a_3_0 = tmp_60;
      real_t a_3_1 = tmp_69;
      real_t a_3_2 = tmp_72;
      real_t a_3_3 = -0.16666666666666663*tmp_73;
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

   void p1_pspg_affine_q2::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t tmp_8 = -p_affine_0_2;
      real_t tmp_9 = p_affine_3_2 + tmp_8;
      real_t tmp_10 = p_affine_1_2 + tmp_8;
      real_t tmp_11 = p_affine_3_1 + tmp_2;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = p_affine_2_2 + tmp_8;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = tmp_14*tmp_6;
      real_t tmp_16 = tmp_1*tmp_11;
      real_t tmp_17 = tmp_14*tmp_3;
      real_t tmp_18 = 1.0 / (tmp_10*tmp_12 - tmp_10*tmp_17 + tmp_13*tmp_15 - tmp_13*tmp_16 + tmp_4*tmp_9 - tmp_7*tmp_9);
      real_t tmp_19 = tmp_18*(tmp_4 - tmp_7);
      real_t tmp_20 = tmp_18*(tmp_15 - tmp_16);
      real_t tmp_21 = tmp_18*(tmp_12 - tmp_17);
      real_t tmp_22 = -tmp_19 - tmp_20 - tmp_21;
      real_t tmp_23 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_5);
      real_t tmp_24 = tmp_18*(tmp_1*tmp_9 - tmp_10*tmp_14);
      real_t tmp_25 = tmp_18*(tmp_13*tmp_14 - tmp_5*tmp_9);
      real_t tmp_26 = -tmp_23 - tmp_24 - tmp_25;
      real_t tmp_27 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_6);
      real_t tmp_28 = tmp_18*(tmp_10*tmp_11 - tmp_6*tmp_9);
      real_t tmp_29 = tmp_18*(-tmp_11*tmp_13 + tmp_3*tmp_9);
      real_t tmp_30 = -tmp_27 - tmp_28 - tmp_29;
      real_t tmp_31 = p_affine_0_0*p_affine_1_1;
      real_t tmp_32 = p_affine_0_0*p_affine_1_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_2;
      real_t tmp_34 = p_affine_0_1*p_affine_1_0;
      real_t tmp_35 = p_affine_0_1*p_affine_1_2;
      real_t tmp_36 = p_affine_2_2*p_affine_3_0;
      real_t tmp_37 = p_affine_0_2*p_affine_1_0;
      real_t tmp_38 = p_affine_0_2*p_affine_1_1;
      real_t tmp_39 = p_affine_2_0*p_affine_3_1;
      real_t tmp_40 = p_affine_2_2*p_affine_3_1;
      real_t tmp_41 = p_affine_2_0*p_affine_3_2;
      real_t tmp_42 = p_affine_2_1*p_affine_3_0;
      real_t tmp_43 = 0.025237786011557493*std::pow(std::abs(p_affine_0_0*tmp_33 - p_affine_0_0*tmp_40 + p_affine_0_1*tmp_36 - p_affine_0_1*tmp_41 + p_affine_0_2*tmp_39 - p_affine_0_2*tmp_42 - p_affine_1_0*tmp_33 + p_affine_1_0*tmp_40 - p_affine_1_1*tmp_36 + p_affine_1_1*tmp_41 - p_affine_1_2*tmp_39 + p_affine_1_2*tmp_42 + p_affine_2_0*tmp_35 - p_affine_2_0*tmp_38 - p_affine_2_1*tmp_32 + p_affine_2_1*tmp_37 + p_affine_2_2*tmp_31 - p_affine_2_2*tmp_34 - p_affine_3_0*tmp_35 + p_affine_3_0*tmp_38 + p_affine_3_1*tmp_32 - p_affine_3_1*tmp_37 - p_affine_3_2*tmp_31 + p_affine_3_2*tmp_34), 1.6666666666666665);
      real_t tmp_44 = tmp_43*((tmp_22*tmp_22) + (tmp_26*tmp_26) + (tmp_30*tmp_30));
      real_t tmp_45 = tmp_43*(tmp_21*tmp_22 + tmp_25*tmp_26 + tmp_29*tmp_30);
      real_t tmp_46 = tmp_43*(tmp_20*tmp_22 + tmp_24*tmp_26 + tmp_28*tmp_30);
      real_t tmp_47 = tmp_43*(tmp_19*tmp_22 + tmp_23*tmp_26 + tmp_27*tmp_30);
      real_t a_0_0 = -0.16666666666666663*tmp_44;
      real_t a_0_1 = -0.16666666666666663*tmp_45;
      real_t a_0_2 = -0.16666666666666663*tmp_46;
      real_t a_0_3 = -0.16666666666666663*tmp_47;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

} // namespace forms
} // namespace hyteg

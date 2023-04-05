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

#include "p1_div_k_grad_affine_q1.hpp"

namespace hyteg {
namespace forms {

   void p1_div_k_grad_affine_q1::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Scalar_Variable_Coefficient_2D_0_0 = 0;
      Scalar_Variable_Coefficient_2D( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Scalar_Variable_Coefficient_2D_0_0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0);
      real_t tmp_5 = 1.0 / (tmp_4);
      real_t tmp_6 = tmp_1*tmp_5;
      real_t tmp_7 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_8 = tmp_5*tmp_7;
      real_t tmp_9 = -tmp_6 - tmp_8;
      real_t tmp_10 = Scalar_Variable_Coefficient_2D_0_0*tmp_6;
      real_t tmp_11 = Scalar_Variable_Coefficient_2D_0_0*tmp_8;
      real_t tmp_12 = -tmp_10 - tmp_11;
      real_t tmp_13 = tmp_3*tmp_5;
      real_t tmp_14 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_15 = tmp_14*tmp_5;
      real_t tmp_16 = -tmp_13 - tmp_15;
      real_t tmp_17 = Scalar_Variable_Coefficient_2D_0_0*tmp_13;
      real_t tmp_18 = Scalar_Variable_Coefficient_2D_0_0*tmp_15;
      real_t tmp_19 = -tmp_17 - tmp_18;
      real_t tmp_20 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_21 = Scalar_Variable_Coefficient_2D_0_0/(tmp_4*tmp_4);
      real_t tmp_22 = tmp_20*(tmp_1*tmp_21*tmp_7 + tmp_14*tmp_21*tmp_3);
      real_t a_0_0 = tmp_20*(tmp_12*tmp_9 + tmp_16*tmp_19);
      real_t a_0_1 = tmp_20*(tmp_11*tmp_9 + tmp_16*tmp_17);
      real_t a_0_2 = tmp_20*(tmp_10*tmp_9 + tmp_16*tmp_18);
      real_t a_1_0 = tmp_20*(tmp_12*tmp_8 + tmp_13*tmp_19);
      real_t a_1_1 = tmp_20*(tmp_21*(tmp_3*tmp_3) + tmp_21*(tmp_7*tmp_7));
      real_t a_1_2 = tmp_22;
      real_t a_2_0 = tmp_20*(tmp_12*tmp_6 + tmp_15*tmp_19);
      real_t a_2_1 = tmp_22;
      real_t a_2_2 = tmp_20*((tmp_1*tmp_1)*tmp_21 + (tmp_14*tmp_14)*tmp_21);
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

   void p1_div_k_grad_affine_q1::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Scalar_Variable_Coefficient_2D_0_0 = 0;
      Scalar_Variable_Coefficient_2D( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Scalar_Variable_Coefficient_2D_0_0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = tmp_1*tmp_4;
      real_t tmp_6 = tmp_4*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_7 = -tmp_5 - tmp_6;
      real_t tmp_8 = Scalar_Variable_Coefficient_2D_0_0*tmp_5;
      real_t tmp_9 = Scalar_Variable_Coefficient_2D_0_0*tmp_6;
      real_t tmp_10 = tmp_3*tmp_4;
      real_t tmp_11 = tmp_4*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_12 = -tmp_10 - tmp_11;
      real_t tmp_13 = Scalar_Variable_Coefficient_2D_0_0*tmp_10;
      real_t tmp_14 = Scalar_Variable_Coefficient_2D_0_0*tmp_11;
      real_t tmp_15 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = tmp_15*(tmp_12*(-tmp_13 - tmp_14) + tmp_7*(-tmp_8 - tmp_9));
      real_t a_0_1 = tmp_15*(tmp_12*tmp_13 + tmp_7*tmp_9);
      real_t a_0_2 = tmp_15*(tmp_12*tmp_14 + tmp_7*tmp_8);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_div_k_grad_affine_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_0_0 = 0;
      Scalar_Variable_Coefficient_3D( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_0_0 );
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
      real_t tmp_27 = Scalar_Variable_Coefficient_3D_0_0*tmp_21;
      real_t tmp_28 = Scalar_Variable_Coefficient_3D_0_0*tmp_23;
      real_t tmp_29 = Scalar_Variable_Coefficient_3D_0_0*tmp_25;
      real_t tmp_30 = -tmp_27 - tmp_28 - tmp_29;
      real_t tmp_31 = -tmp_1*tmp_14 + tmp_11*tmp_5;
      real_t tmp_32 = tmp_20*tmp_31;
      real_t tmp_33 = tmp_1*tmp_10 - tmp_11*tmp_15;
      real_t tmp_34 = tmp_20*tmp_33;
      real_t tmp_35 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_36 = tmp_20*tmp_35;
      real_t tmp_37 = -tmp_32 - tmp_34 - tmp_36;
      real_t tmp_38 = Scalar_Variable_Coefficient_3D_0_0*tmp_32;
      real_t tmp_39 = Scalar_Variable_Coefficient_3D_0_0*tmp_34;
      real_t tmp_40 = Scalar_Variable_Coefficient_3D_0_0*tmp_36;
      real_t tmp_41 = -tmp_38 - tmp_39 - tmp_40;
      real_t tmp_42 = -tmp_11*tmp_3 + tmp_14*tmp_6;
      real_t tmp_43 = tmp_20*tmp_42;
      real_t tmp_44 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_45 = tmp_20*tmp_44;
      real_t tmp_46 = tmp_10*tmp_3 - tmp_12*tmp_14;
      real_t tmp_47 = tmp_20*tmp_46;
      real_t tmp_48 = -tmp_43 - tmp_45 - tmp_47;
      real_t tmp_49 = Scalar_Variable_Coefficient_3D_0_0*tmp_43;
      real_t tmp_50 = Scalar_Variable_Coefficient_3D_0_0*tmp_45;
      real_t tmp_51 = Scalar_Variable_Coefficient_3D_0_0*tmp_47;
      real_t tmp_52 = -tmp_49 - tmp_50 - tmp_51;
      real_t tmp_53 = p_affine_0_0*p_affine_1_1;
      real_t tmp_54 = p_affine_0_0*p_affine_1_2;
      real_t tmp_55 = p_affine_2_1*p_affine_3_2;
      real_t tmp_56 = p_affine_0_1*p_affine_1_0;
      real_t tmp_57 = p_affine_0_1*p_affine_1_2;
      real_t tmp_58 = p_affine_2_2*p_affine_3_0;
      real_t tmp_59 = p_affine_0_2*p_affine_1_0;
      real_t tmp_60 = p_affine_0_2*p_affine_1_1;
      real_t tmp_61 = p_affine_2_0*p_affine_3_1;
      real_t tmp_62 = p_affine_2_2*p_affine_3_1;
      real_t tmp_63 = p_affine_2_0*p_affine_3_2;
      real_t tmp_64 = p_affine_2_1*p_affine_3_0;
      real_t tmp_65 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_55 - p_affine_0_0*tmp_62 + p_affine_0_1*tmp_58 - p_affine_0_1*tmp_63 + p_affine_0_2*tmp_61 - p_affine_0_2*tmp_64 - p_affine_1_0*tmp_55 + p_affine_1_0*tmp_62 - p_affine_1_1*tmp_58 + p_affine_1_1*tmp_63 - p_affine_1_2*tmp_61 + p_affine_1_2*tmp_64 + p_affine_2_0*tmp_57 - p_affine_2_0*tmp_60 - p_affine_2_1*tmp_54 + p_affine_2_1*tmp_59 + p_affine_2_2*tmp_53 - p_affine_2_2*tmp_56 - p_affine_3_0*tmp_57 + p_affine_3_0*tmp_60 + p_affine_3_1*tmp_54 - p_affine_3_1*tmp_59 - p_affine_3_2*tmp_53 + p_affine_3_2*tmp_56);
      real_t tmp_66 = Scalar_Variable_Coefficient_3D_0_0/(tmp_19*tmp_19);
      real_t tmp_67 = tmp_24*tmp_66;
      real_t tmp_68 = tmp_35*tmp_66;
      real_t tmp_69 = tmp_46*tmp_66;
      real_t tmp_70 = tmp_65*(tmp_22*tmp_67 + tmp_33*tmp_68 + tmp_44*tmp_69);
      real_t tmp_71 = tmp_65*(tmp_31*tmp_68 + tmp_42*tmp_69 + tmp_67*tmp_8);
      real_t tmp_72 = tmp_65*(tmp_22*tmp_66*tmp_8 + tmp_31*tmp_33*tmp_66 + tmp_42*tmp_44*tmp_66);
      real_t a_0_0 = tmp_65*(tmp_26*tmp_30 + tmp_37*tmp_41 + tmp_48*tmp_52);
      real_t a_0_1 = tmp_65*(tmp_26*tmp_29 + tmp_37*tmp_40 + tmp_48*tmp_51);
      real_t a_0_2 = tmp_65*(tmp_26*tmp_28 + tmp_37*tmp_39 + tmp_48*tmp_50);
      real_t a_0_3 = tmp_65*(tmp_26*tmp_27 + tmp_37*tmp_38 + tmp_48*tmp_49);
      real_t a_1_0 = tmp_65*(tmp_25*tmp_30 + tmp_36*tmp_41 + tmp_47*tmp_52);
      real_t a_1_1 = tmp_65*((tmp_24*tmp_24)*tmp_66 + (tmp_35*tmp_35)*tmp_66 + (tmp_46*tmp_46)*tmp_66);
      real_t a_1_2 = tmp_70;
      real_t a_1_3 = tmp_71;
      real_t a_2_0 = tmp_65*(tmp_23*tmp_30 + tmp_34*tmp_41 + tmp_45*tmp_52);
      real_t a_2_1 = tmp_70;
      real_t a_2_2 = tmp_65*((tmp_22*tmp_22)*tmp_66 + (tmp_33*tmp_33)*tmp_66 + (tmp_44*tmp_44)*tmp_66);
      real_t a_2_3 = tmp_72;
      real_t a_3_0 = tmp_65*(tmp_21*tmp_30 + tmp_32*tmp_41 + tmp_43*tmp_52);
      real_t a_3_1 = tmp_71;
      real_t a_3_2 = tmp_72;
      real_t a_3_3 = tmp_65*((tmp_31*tmp_31)*tmp_66 + (tmp_42*tmp_42)*tmp_66 + tmp_66*(tmp_8*tmp_8));
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

   void p1_div_k_grad_affine_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_0_0 = 0;
      Scalar_Variable_Coefficient_3D( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_0_0 );
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
      real_t tmp_23 = Scalar_Variable_Coefficient_3D_0_0*tmp_19;
      real_t tmp_24 = Scalar_Variable_Coefficient_3D_0_0*tmp_20;
      real_t tmp_25 = Scalar_Variable_Coefficient_3D_0_0*tmp_21;
      real_t tmp_26 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_5);
      real_t tmp_27 = tmp_18*(tmp_1*tmp_9 - tmp_10*tmp_14);
      real_t tmp_28 = tmp_18*(tmp_13*tmp_14 - tmp_5*tmp_9);
      real_t tmp_29 = -tmp_26 - tmp_27 - tmp_28;
      real_t tmp_30 = Scalar_Variable_Coefficient_3D_0_0*tmp_26;
      real_t tmp_31 = Scalar_Variable_Coefficient_3D_0_0*tmp_27;
      real_t tmp_32 = Scalar_Variable_Coefficient_3D_0_0*tmp_28;
      real_t tmp_33 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_6);
      real_t tmp_34 = tmp_18*(tmp_10*tmp_11 - tmp_6*tmp_9);
      real_t tmp_35 = tmp_18*(-tmp_11*tmp_13 + tmp_3*tmp_9);
      real_t tmp_36 = -tmp_33 - tmp_34 - tmp_35;
      real_t tmp_37 = Scalar_Variable_Coefficient_3D_0_0*tmp_33;
      real_t tmp_38 = Scalar_Variable_Coefficient_3D_0_0*tmp_34;
      real_t tmp_39 = Scalar_Variable_Coefficient_3D_0_0*tmp_35;
      real_t tmp_40 = p_affine_0_0*p_affine_1_1;
      real_t tmp_41 = p_affine_0_0*p_affine_1_2;
      real_t tmp_42 = p_affine_2_1*p_affine_3_2;
      real_t tmp_43 = p_affine_0_1*p_affine_1_0;
      real_t tmp_44 = p_affine_0_1*p_affine_1_2;
      real_t tmp_45 = p_affine_2_2*p_affine_3_0;
      real_t tmp_46 = p_affine_0_2*p_affine_1_0;
      real_t tmp_47 = p_affine_0_2*p_affine_1_1;
      real_t tmp_48 = p_affine_2_0*p_affine_3_1;
      real_t tmp_49 = p_affine_2_2*p_affine_3_1;
      real_t tmp_50 = p_affine_2_0*p_affine_3_2;
      real_t tmp_51 = p_affine_2_1*p_affine_3_0;
      real_t tmp_52 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_42 - p_affine_0_0*tmp_49 + p_affine_0_1*tmp_45 - p_affine_0_1*tmp_50 + p_affine_0_2*tmp_48 - p_affine_0_2*tmp_51 - p_affine_1_0*tmp_42 + p_affine_1_0*tmp_49 - p_affine_1_1*tmp_45 + p_affine_1_1*tmp_50 - p_affine_1_2*tmp_48 + p_affine_1_2*tmp_51 + p_affine_2_0*tmp_44 - p_affine_2_0*tmp_47 - p_affine_2_1*tmp_41 + p_affine_2_1*tmp_46 + p_affine_2_2*tmp_40 - p_affine_2_2*tmp_43 - p_affine_3_0*tmp_44 + p_affine_3_0*tmp_47 + p_affine_3_1*tmp_41 - p_affine_3_1*tmp_46 - p_affine_3_2*tmp_40 + p_affine_3_2*tmp_43);
      real_t a_0_0 = tmp_52*(tmp_22*(-tmp_23 - tmp_24 - tmp_25) + tmp_29*(-tmp_30 - tmp_31 - tmp_32) + tmp_36*(-tmp_37 - tmp_38 - tmp_39));
      real_t a_0_1 = tmp_52*(tmp_22*tmp_25 + tmp_29*tmp_32 + tmp_36*tmp_39);
      real_t a_0_2 = tmp_52*(tmp_22*tmp_24 + tmp_29*tmp_31 + tmp_36*tmp_38);
      real_t a_0_3 = tmp_52*(tmp_22*tmp_23 + tmp_29*tmp_30 + tmp_36*tmp_37);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_div_k_grad_affine_q1::Scalar_Variable_Coefficient_2D( real_t in_0, real_t in_1, real_t * out_0 ) const
   {
      *out_0 = callback2D( Point3D( in_0, in_1, 0 ) );
   }

   void p1_div_k_grad_affine_q1::Scalar_Variable_Coefficient_3D( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback3D( Point3D( in_0, in_1, in_2 ) );
   }

} // namespace forms
} // namespace hyteg

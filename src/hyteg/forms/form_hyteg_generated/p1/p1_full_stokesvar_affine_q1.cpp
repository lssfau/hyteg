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

#include "p1_full_stokesvar_affine_q1.hpp"

namespace hyteg {
namespace forms {

   void p1_full_stokesvar_0_0_affine_q1::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0);
      real_t tmp_5 = 1.0 / (tmp_4);
      real_t tmp_6 = tmp_1*tmp_5;
      real_t tmp_7 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_8 = tmp_5*tmp_7;
      real_t tmp_9 = 1.0*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_10 = tmp_3*tmp_5;
      real_t tmp_11 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = -tmp_10 - tmp_12;
      real_t tmp_14 = (2.0/3.0)*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_15 = -1.0*tmp_10 - 1.0*tmp_12;
      real_t tmp_16 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_17 = 2.0*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_18 = tmp_17*(-0.5*tmp_6 - 0.5*tmp_8);
      real_t tmp_19 = tmp_13*tmp_14;
      real_t tmp_20 = tmp_15*tmp_17;
      real_t tmp_21 = tmp_16*(-tmp_10*tmp_19 + tmp_10*tmp_20 + tmp_18*tmp_8);
      real_t tmp_22 = tmp_16*(-tmp_12*tmp_19 + tmp_12*tmp_20 + tmp_18*tmp_6);
      real_t tmp_23 = 1.0 / (tmp_4*tmp_4);
      real_t tmp_24 = tmp_23*tmp_9;
      real_t tmp_25 = 1.3333333333333335*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_23;
      real_t tmp_26 = tmp_16*(tmp_1*tmp_24*tmp_7 + tmp_11*tmp_25*tmp_3);
      real_t a_0_0 = tmp_16*(2*Scalar_Variable_Coefficient_2D_mu_out0_id0*(tmp_15*tmp_15) - (tmp_13*tmp_13)*tmp_14 + tmp_9*((-tmp_6 - tmp_8)*(-tmp_6 - tmp_8)));
      real_t a_0_1 = tmp_21;
      real_t a_0_2 = tmp_22;
      real_t a_1_0 = tmp_21;
      real_t a_1_1 = tmp_16*(tmp_24*(tmp_7*tmp_7) + tmp_25*(tmp_3*tmp_3));
      real_t a_1_2 = tmp_26;
      real_t a_2_0 = tmp_22;
      real_t a_2_1 = tmp_26;
      real_t a_2_2 = tmp_16*((tmp_1*tmp_1)*tmp_24 + (tmp_11*tmp_11)*tmp_25);
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

   void p1_full_stokesvar_0_0_affine_q1::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = tmp_1*tmp_4;
      real_t tmp_6 = tmp_4*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_7 = tmp_3*tmp_4;
      real_t tmp_8 = tmp_4*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_9 = -tmp_7 - tmp_8;
      real_t tmp_10 = (2.0/3.0)*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_11 = -1.0*tmp_7 - 1.0*tmp_8;
      real_t tmp_12 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_13 = 2.0*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_14 = tmp_13*(-0.5*tmp_5 - 0.5*tmp_6);
      real_t tmp_15 = tmp_10*tmp_9;
      real_t tmp_16 = tmp_11*tmp_13;
      real_t a_0_0 = tmp_12*(2*Scalar_Variable_Coefficient_2D_mu_out0_id0*(tmp_11*tmp_11) + 1.0*Scalar_Variable_Coefficient_2D_mu_out0_id0*((-tmp_5 - tmp_6)*(-tmp_5 - tmp_6)) - tmp_10*(tmp_9*tmp_9));
      real_t a_0_1 = tmp_12*(tmp_14*tmp_6 - tmp_15*tmp_7 + tmp_16*tmp_7);
      real_t a_0_2 = tmp_12*(tmp_14*tmp_5 - tmp_15*tmp_8 + tmp_16*tmp_8);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_full_stokesvar_0_0_affine_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_2;
      real_t tmp_8 = p_affine_3_2 + tmp_7;
      real_t tmp_9 = tmp_3*tmp_8;
      real_t tmp_10 = p_affine_3_1 + tmp_2;
      real_t tmp_11 = p_affine_1_2 + tmp_7;
      real_t tmp_12 = tmp_10*tmp_11;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_14*tmp_5;
      real_t tmp_16 = tmp_10*tmp_14;
      real_t tmp_17 = tmp_5*tmp_8;
      real_t tmp_18 = tmp_11*tmp_3;
      real_t tmp_19 = -tmp_1*tmp_16 + tmp_1*tmp_9 + tmp_12*tmp_4 + tmp_13*tmp_15 - tmp_13*tmp_18 - tmp_17*tmp_4;
      real_t tmp_20 = 1.0 / (tmp_19);
      real_t tmp_21 = tmp_20*tmp_6;
      real_t tmp_22 = -tmp_1*tmp_10 + tmp_13*tmp_5;
      real_t tmp_23 = tmp_20*tmp_22;
      real_t tmp_24 = tmp_10*tmp_4 - tmp_13*tmp_3;
      real_t tmp_25 = tmp_20*tmp_24;
      real_t tmp_26 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_27 = -tmp_1*tmp_14 + tmp_11*tmp_4;
      real_t tmp_28 = tmp_20*tmp_27;
      real_t tmp_29 = tmp_1*tmp_8 - tmp_11*tmp_13;
      real_t tmp_30 = tmp_20*tmp_29;
      real_t tmp_31 = tmp_13*tmp_14 - tmp_4*tmp_8;
      real_t tmp_32 = tmp_20*tmp_31;
      real_t tmp_33 = tmp_15 - tmp_18;
      real_t tmp_34 = tmp_20*tmp_33;
      real_t tmp_35 = tmp_12 - tmp_17;
      real_t tmp_36 = tmp_20*tmp_35;
      real_t tmp_37 = -tmp_16 + tmp_9;
      real_t tmp_38 = tmp_20*tmp_37;
      real_t tmp_39 = -tmp_34 - tmp_36 - tmp_38;
      real_t tmp_40 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_41 = -1.0*tmp_34 - 1.0*tmp_36 - 1.0*tmp_38;
      real_t tmp_42 = p_affine_0_0*p_affine_1_1;
      real_t tmp_43 = p_affine_0_0*p_affine_1_2;
      real_t tmp_44 = p_affine_2_1*p_affine_3_2;
      real_t tmp_45 = p_affine_0_1*p_affine_1_0;
      real_t tmp_46 = p_affine_0_1*p_affine_1_2;
      real_t tmp_47 = p_affine_2_2*p_affine_3_0;
      real_t tmp_48 = p_affine_0_2*p_affine_1_0;
      real_t tmp_49 = p_affine_0_2*p_affine_1_1;
      real_t tmp_50 = p_affine_2_0*p_affine_3_1;
      real_t tmp_51 = p_affine_2_2*p_affine_3_1;
      real_t tmp_52 = p_affine_2_0*p_affine_3_2;
      real_t tmp_53 = p_affine_2_1*p_affine_3_0;
      real_t tmp_54 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_44 - p_affine_0_0*tmp_51 + p_affine_0_1*tmp_47 - p_affine_0_1*tmp_52 + p_affine_0_2*tmp_50 - p_affine_0_2*tmp_53 - p_affine_1_0*tmp_44 + p_affine_1_0*tmp_51 - p_affine_1_1*tmp_47 + p_affine_1_1*tmp_52 - p_affine_1_2*tmp_50 + p_affine_1_2*tmp_53 + p_affine_2_0*tmp_46 - p_affine_2_0*tmp_49 - p_affine_2_1*tmp_43 + p_affine_2_1*tmp_48 + p_affine_2_2*tmp_42 - p_affine_2_2*tmp_45 - p_affine_3_0*tmp_46 + p_affine_3_0*tmp_49 + p_affine_3_1*tmp_43 - p_affine_3_1*tmp_48 - p_affine_3_2*tmp_42 + p_affine_3_2*tmp_45);
      real_t tmp_55 = 2.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_56 = tmp_55*(-0.5*tmp_21 - 0.5*tmp_23 - 0.5*tmp_25);
      real_t tmp_57 = tmp_55*(-0.5*tmp_28 - 0.5*tmp_30 - 0.5*tmp_32);
      real_t tmp_58 = tmp_39*tmp_40;
      real_t tmp_59 = tmp_41*tmp_55;
      real_t tmp_60 = tmp_54*(tmp_25*tmp_56 + tmp_32*tmp_57 - tmp_38*tmp_58 + tmp_38*tmp_59);
      real_t tmp_61 = tmp_54*(tmp_23*tmp_56 + tmp_30*tmp_57 - tmp_36*tmp_58 + tmp_36*tmp_59);
      real_t tmp_62 = tmp_54*(tmp_21*tmp_56 + tmp_28*tmp_57 - tmp_34*tmp_58 + tmp_34*tmp_59);
      real_t tmp_63 = 1.0 / (tmp_19*tmp_19);
      real_t tmp_64 = tmp_26*tmp_63;
      real_t tmp_65 = 1.3333333333333335*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_63;
      real_t tmp_66 = tmp_24*tmp_64;
      real_t tmp_67 = tmp_31*tmp_64;
      real_t tmp_68 = tmp_37*tmp_65;
      real_t tmp_69 = tmp_54*(tmp_22*tmp_66 + tmp_29*tmp_67 + tmp_35*tmp_68);
      real_t tmp_70 = tmp_54*(tmp_27*tmp_67 + tmp_33*tmp_68 + tmp_6*tmp_66);
      real_t tmp_71 = tmp_54*(tmp_22*tmp_6*tmp_64 + tmp_27*tmp_29*tmp_64 + tmp_33*tmp_35*tmp_65);
      real_t a_0_0 = tmp_54*(2*Scalar_Variable_Coefficient_3D_mu_out0_id0*(tmp_41*tmp_41) + tmp_26*((-tmp_21 - tmp_23 - tmp_25)*(-tmp_21 - tmp_23 - tmp_25)) + tmp_26*((-tmp_28 - tmp_30 - tmp_32)*(-tmp_28 - tmp_30 - tmp_32)) - (tmp_39*tmp_39)*tmp_40);
      real_t a_0_1 = tmp_60;
      real_t a_0_2 = tmp_61;
      real_t a_0_3 = tmp_62;
      real_t a_1_0 = tmp_60;
      real_t a_1_1 = tmp_54*((tmp_24*tmp_24)*tmp_64 + (tmp_31*tmp_31)*tmp_64 + (tmp_37*tmp_37)*tmp_65);
      real_t a_1_2 = tmp_69;
      real_t a_1_3 = tmp_70;
      real_t a_2_0 = tmp_61;
      real_t a_2_1 = tmp_69;
      real_t a_2_2 = tmp_54*((tmp_22*tmp_22)*tmp_64 + (tmp_29*tmp_29)*tmp_64 + (tmp_35*tmp_35)*tmp_65);
      real_t a_2_3 = tmp_71;
      real_t a_3_0 = tmp_62;
      real_t a_3_1 = tmp_70;
      real_t a_3_2 = tmp_71;
      real_t a_3_3 = tmp_54*((tmp_27*tmp_27)*tmp_64 + (tmp_33*tmp_33)*tmp_65 + (tmp_6*tmp_6)*tmp_64);
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

   void p1_full_stokesvar_0_0_affine_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = -p_affine_0_2;
      real_t tmp_7 = p_affine_3_2 + tmp_6;
      real_t tmp_8 = tmp_3*tmp_7;
      real_t tmp_9 = p_affine_3_1 + tmp_2;
      real_t tmp_10 = p_affine_1_2 + tmp_6;
      real_t tmp_11 = tmp_10*tmp_9;
      real_t tmp_12 = p_affine_3_0 + tmp_0;
      real_t tmp_13 = p_affine_2_2 + tmp_6;
      real_t tmp_14 = tmp_13*tmp_5;
      real_t tmp_15 = tmp_13*tmp_9;
      real_t tmp_16 = tmp_5*tmp_7;
      real_t tmp_17 = tmp_10*tmp_3;
      real_t tmp_18 = 1.0 / (-tmp_1*tmp_15 + tmp_1*tmp_8 + tmp_11*tmp_4 + tmp_12*tmp_14 - tmp_12*tmp_17 - tmp_16*tmp_4);
      real_t tmp_19 = tmp_18*(tmp_1*tmp_3 - tmp_4*tmp_5);
      real_t tmp_20 = tmp_18*(-tmp_1*tmp_9 + tmp_12*tmp_5);
      real_t tmp_21 = tmp_18*(-tmp_12*tmp_3 + tmp_4*tmp_9);
      real_t tmp_22 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_23 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_4);
      real_t tmp_24 = tmp_18*(tmp_1*tmp_7 - tmp_10*tmp_12);
      real_t tmp_25 = tmp_18*(tmp_12*tmp_13 - tmp_4*tmp_7);
      real_t tmp_26 = tmp_18*(tmp_14 - tmp_17);
      real_t tmp_27 = tmp_18*(tmp_11 - tmp_16);
      real_t tmp_28 = tmp_18*(-tmp_15 + tmp_8);
      real_t tmp_29 = -tmp_26 - tmp_27 - tmp_28;
      real_t tmp_30 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_31 = -1.0*tmp_26 - 1.0*tmp_27 - 1.0*tmp_28;
      real_t tmp_32 = p_affine_0_0*p_affine_1_1;
      real_t tmp_33 = p_affine_0_0*p_affine_1_2;
      real_t tmp_34 = p_affine_2_1*p_affine_3_2;
      real_t tmp_35 = p_affine_0_1*p_affine_1_0;
      real_t tmp_36 = p_affine_0_1*p_affine_1_2;
      real_t tmp_37 = p_affine_2_2*p_affine_3_0;
      real_t tmp_38 = p_affine_0_2*p_affine_1_0;
      real_t tmp_39 = p_affine_0_2*p_affine_1_1;
      real_t tmp_40 = p_affine_2_0*p_affine_3_1;
      real_t tmp_41 = p_affine_2_2*p_affine_3_1;
      real_t tmp_42 = p_affine_2_0*p_affine_3_2;
      real_t tmp_43 = p_affine_2_1*p_affine_3_0;
      real_t tmp_44 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_34 - p_affine_0_0*tmp_41 + p_affine_0_1*tmp_37 - p_affine_0_1*tmp_42 + p_affine_0_2*tmp_40 - p_affine_0_2*tmp_43 - p_affine_1_0*tmp_34 + p_affine_1_0*tmp_41 - p_affine_1_1*tmp_37 + p_affine_1_1*tmp_42 - p_affine_1_2*tmp_40 + p_affine_1_2*tmp_43 + p_affine_2_0*tmp_36 - p_affine_2_0*tmp_39 - p_affine_2_1*tmp_33 + p_affine_2_1*tmp_38 + p_affine_2_2*tmp_32 - p_affine_2_2*tmp_35 - p_affine_3_0*tmp_36 + p_affine_3_0*tmp_39 + p_affine_3_1*tmp_33 - p_affine_3_1*tmp_38 - p_affine_3_2*tmp_32 + p_affine_3_2*tmp_35);
      real_t tmp_45 = 2.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_46 = tmp_45*(-0.5*tmp_19 - 0.5*tmp_20 - 0.5*tmp_21);
      real_t tmp_47 = tmp_45*(-0.5*tmp_23 - 0.5*tmp_24 - 0.5*tmp_25);
      real_t tmp_48 = tmp_29*tmp_30;
      real_t tmp_49 = tmp_31*tmp_45;
      real_t a_0_0 = tmp_44*(2*Scalar_Variable_Coefficient_3D_mu_out0_id0*(tmp_31*tmp_31) + tmp_22*((-tmp_19 - tmp_20 - tmp_21)*(-tmp_19 - tmp_20 - tmp_21)) + tmp_22*((-tmp_23 - tmp_24 - tmp_25)*(-tmp_23 - tmp_24 - tmp_25)) - (tmp_29*tmp_29)*tmp_30);
      real_t a_0_1 = tmp_44*(tmp_21*tmp_46 + tmp_25*tmp_47 - tmp_28*tmp_48 + tmp_28*tmp_49);
      real_t a_0_2 = tmp_44*(tmp_20*tmp_46 + tmp_24*tmp_47 - tmp_27*tmp_48 + tmp_27*tmp_49);
      real_t a_0_3 = tmp_44*(tmp_19*tmp_46 + tmp_23*tmp_47 - tmp_26*tmp_48 + tmp_26*tmp_49);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_full_stokesvar_0_0_affine_q1::Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_mu( Point3D( {in_0, in_1, 0} ) );
   }

   void p1_full_stokesvar_0_0_affine_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

   void p1_full_stokesvar_0_1_affine_q1::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = tmp_4 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0);
      real_t tmp_6 = 1.0 / (tmp_5);
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_9 = tmp_6*tmp_8;
      real_t tmp_10 = -tmp_7 - tmp_9;
      real_t tmp_11 = tmp_3*tmp_6;
      real_t tmp_12 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_13 = tmp_12*tmp_6;
      real_t tmp_14 = (2.0/3.0)*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_15 = tmp_14*(-tmp_11 - tmp_13);
      real_t tmp_16 = -0.5*tmp_11 - 0.5*tmp_13;
      real_t tmp_17 = Scalar_Variable_Coefficient_2D_mu_out0_id0*(-0.5*tmp_7 - 0.5*tmp_9);
      real_t tmp_18 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_19 = 2.0*tmp_17;
      real_t tmp_20 = 2.0*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_16;
      real_t tmp_21 = tmp_10*tmp_14;
      real_t tmp_22 = 1.0 / (tmp_5*tmp_5);
      real_t tmp_23 = 0.33333333333333337*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_18*tmp_22;
      real_t tmp_24 = tmp_14*tmp_22;
      real_t tmp_25 = tmp_12*tmp_8;
      real_t tmp_26 = 1.0*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_22;
      real_t a_0_0 = tmp_18*(-tmp_10*tmp_15 + 4*tmp_16*tmp_17);
      real_t a_0_1 = tmp_18*(tmp_11*tmp_19 - tmp_15*tmp_9);
      real_t a_0_2 = tmp_18*(tmp_13*tmp_19 - tmp_15*tmp_7);
      real_t a_1_0 = tmp_18*(-tmp_11*tmp_21 + tmp_20*tmp_9);
      real_t a_1_1 = tmp_23*tmp_3*tmp_8;
      real_t a_1_2 = tmp_18*(-tmp_24*tmp_4 + tmp_25*tmp_26);
      real_t a_2_0 = tmp_18*(-tmp_13*tmp_21 + tmp_20*tmp_7);
      real_t a_2_1 = tmp_18*(-tmp_24*tmp_25 + tmp_26*tmp_4);
      real_t a_2_2 = tmp_1*tmp_12*tmp_23;
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

   void p1_full_stokesvar_0_1_affine_q1::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = tmp_1*tmp_4;
      real_t tmp_6 = tmp_4*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_7 = tmp_3*tmp_4;
      real_t tmp_8 = tmp_4*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_9 = (2.0/3.0)*Scalar_Variable_Coefficient_2D_mu_out0_id0*(-tmp_7 - tmp_8);
      real_t tmp_10 = Scalar_Variable_Coefficient_2D_mu_out0_id0*(-0.5*tmp_5 - 0.5*tmp_6);
      real_t tmp_11 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_12 = 2.0*tmp_10;
      real_t a_0_0 = tmp_11*(4*tmp_10*(-0.5*tmp_7 - 0.5*tmp_8) - tmp_9*(-tmp_5 - tmp_6));
      real_t a_0_1 = tmp_11*(tmp_12*tmp_7 - tmp_6*tmp_9);
      real_t a_0_2 = tmp_11*(tmp_12*tmp_8 - tmp_5*tmp_9);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_full_stokesvar_0_1_affine_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
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
      real_t tmp_19 = tmp_10*tmp_12 - tmp_10*tmp_18 + tmp_13*tmp_4 - tmp_13*tmp_7 + tmp_14*tmp_16 - tmp_14*tmp_17;
      real_t tmp_20 = 1.0 / (tmp_19);
      real_t tmp_21 = tmp_20*tmp_8;
      real_t tmp_22 = tmp_12 - tmp_18;
      real_t tmp_23 = tmp_20*tmp_22;
      real_t tmp_24 = tmp_16 - tmp_17;
      real_t tmp_25 = tmp_20*tmp_24;
      real_t tmp_26 = -tmp_21 - tmp_23 - tmp_25;
      real_t tmp_27 = -tmp_10*tmp_3 + tmp_14*tmp_6;
      real_t tmp_28 = tmp_20*tmp_27;
      real_t tmp_29 = -tmp_11*tmp_14 + tmp_13*tmp_3;
      real_t tmp_30 = tmp_20*tmp_29;
      real_t tmp_31 = tmp_10*tmp_11 - tmp_13*tmp_6;
      real_t tmp_32 = tmp_20*tmp_31;
      real_t tmp_33 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_34 = tmp_33*(-tmp_28 - tmp_30 - tmp_32);
      real_t tmp_35 = -0.5*tmp_28 - 0.5*tmp_30 - 0.5*tmp_32;
      real_t tmp_36 = Scalar_Variable_Coefficient_3D_mu_out0_id0*(-0.5*tmp_21 - 0.5*tmp_23 - 0.5*tmp_25);
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
      real_t tmp_49 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_39 - p_affine_0_0*tmp_46 + p_affine_0_1*tmp_42 - p_affine_0_1*tmp_47 + p_affine_0_2*tmp_45 - p_affine_0_2*tmp_48 - p_affine_1_0*tmp_39 + p_affine_1_0*tmp_46 - p_affine_1_1*tmp_42 + p_affine_1_1*tmp_47 - p_affine_1_2*tmp_45 + p_affine_1_2*tmp_48 + p_affine_2_0*tmp_41 - p_affine_2_0*tmp_44 - p_affine_2_1*tmp_38 + p_affine_2_1*tmp_43 + p_affine_2_2*tmp_37 - p_affine_2_2*tmp_40 - p_affine_3_0*tmp_41 + p_affine_3_0*tmp_44 + p_affine_3_1*tmp_38 - p_affine_3_1*tmp_43 - p_affine_3_2*tmp_37 + p_affine_3_2*tmp_40);
      real_t tmp_50 = 2.0*tmp_36;
      real_t tmp_51 = 2.0*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_35;
      real_t tmp_52 = tmp_26*tmp_33;
      real_t tmp_53 = 1.0 / (tmp_19*tmp_19);
      real_t tmp_54 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_53;
      real_t tmp_55 = 0.33333333333333337*tmp_49*tmp_54;
      real_t tmp_56 = tmp_33*tmp_53;
      real_t tmp_57 = tmp_31*tmp_56;
      real_t tmp_58 = tmp_24*tmp_29;
      real_t tmp_59 = 1.0*tmp_54;
      real_t tmp_60 = tmp_27*tmp_59;
      real_t tmp_61 = tmp_31*tmp_59;
      real_t tmp_62 = tmp_29*tmp_8;
      real_t tmp_63 = tmp_27*tmp_56;
      real_t a_0_0 = tmp_49*(-tmp_26*tmp_34 + 4*tmp_35*tmp_36);
      real_t a_0_1 = tmp_49*(-tmp_25*tmp_34 + tmp_32*tmp_50);
      real_t a_0_2 = tmp_49*(-tmp_23*tmp_34 + tmp_30*tmp_50);
      real_t a_0_3 = tmp_49*(-tmp_21*tmp_34 + tmp_28*tmp_50);
      real_t a_1_0 = tmp_49*(tmp_25*tmp_51 - tmp_32*tmp_52);
      real_t a_1_1 = tmp_24*tmp_31*tmp_55;
      real_t a_1_2 = tmp_49*(-tmp_22*tmp_57 + tmp_58*tmp_59);
      real_t a_1_3 = tmp_49*(tmp_24*tmp_60 - tmp_57*tmp_8);
      real_t a_2_0 = tmp_49*(tmp_23*tmp_51 - tmp_30*tmp_52);
      real_t a_2_1 = tmp_49*(tmp_22*tmp_61 - tmp_56*tmp_58);
      real_t a_2_2 = tmp_22*tmp_29*tmp_55;
      real_t a_2_3 = tmp_49*(tmp_22*tmp_60 - tmp_56*tmp_62);
      real_t a_3_0 = tmp_49*(tmp_21*tmp_51 - tmp_28*tmp_52);
      real_t a_3_1 = tmp_49*(-tmp_24*tmp_63 + tmp_61*tmp_8);
      real_t a_3_2 = tmp_49*(-tmp_22*tmp_63 + tmp_59*tmp_62);
      real_t a_3_3 = tmp_27*tmp_55*tmp_8;
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

   void p1_full_stokesvar_0_1_affine_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_2_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_2;
      real_t tmp_3 = p_affine_1_2 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_1_0 + tmp_0;
      real_t tmp_6 = p_affine_2_2 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = -p_affine_0_1;
      real_t tmp_9 = p_affine_2_1 + tmp_8;
      real_t tmp_10 = p_affine_3_2 + tmp_2;
      real_t tmp_11 = tmp_10*tmp_5;
      real_t tmp_12 = p_affine_3_1 + tmp_8;
      real_t tmp_13 = p_affine_1_1 + tmp_8;
      real_t tmp_14 = p_affine_3_0 + tmp_0;
      real_t tmp_15 = tmp_14*tmp_6;
      real_t tmp_16 = tmp_1*tmp_10;
      real_t tmp_17 = tmp_14*tmp_3;
      real_t tmp_18 = 1.0 / (tmp_11*tmp_9 + tmp_12*tmp_4 - tmp_12*tmp_7 + tmp_13*tmp_15 - tmp_13*tmp_16 - tmp_17*tmp_9);
      real_t tmp_19 = tmp_18*(tmp_4 - tmp_7);
      real_t tmp_20 = tmp_18*(tmp_11 - tmp_17);
      real_t tmp_21 = tmp_18*(tmp_15 - tmp_16);
      real_t tmp_22 = tmp_18*(tmp_13*tmp_6 - tmp_3*tmp_9);
      real_t tmp_23 = tmp_18*(-tmp_10*tmp_13 + tmp_12*tmp_3);
      real_t tmp_24 = tmp_18*(tmp_10*tmp_9 - tmp_12*tmp_6);
      real_t tmp_25 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0*(-tmp_22 - tmp_23 - tmp_24);
      real_t tmp_26 = Scalar_Variable_Coefficient_3D_mu_out0_id0*(-0.5*tmp_19 - 0.5*tmp_20 - 0.5*tmp_21);
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
      real_t tmp_39 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_29 - p_affine_0_0*tmp_36 + p_affine_0_1*tmp_32 - p_affine_0_1*tmp_37 + p_affine_0_2*tmp_35 - p_affine_0_2*tmp_38 - p_affine_1_0*tmp_29 + p_affine_1_0*tmp_36 - p_affine_1_1*tmp_32 + p_affine_1_1*tmp_37 - p_affine_1_2*tmp_35 + p_affine_1_2*tmp_38 + p_affine_2_0*tmp_31 - p_affine_2_0*tmp_34 - p_affine_2_1*tmp_28 + p_affine_2_1*tmp_33 + p_affine_2_2*tmp_27 - p_affine_2_2*tmp_30 - p_affine_3_0*tmp_31 + p_affine_3_0*tmp_34 + p_affine_3_1*tmp_28 - p_affine_3_1*tmp_33 - p_affine_3_2*tmp_27 + p_affine_3_2*tmp_30);
      real_t tmp_40 = 2.0*tmp_26;
      real_t a_0_0 = tmp_39*(-tmp_25*(-tmp_19 - tmp_20 - tmp_21) + 4*tmp_26*(-0.5*tmp_22 - 0.5*tmp_23 - 0.5*tmp_24));
      real_t a_0_1 = tmp_39*(-tmp_21*tmp_25 + tmp_24*tmp_40);
      real_t a_0_2 = tmp_39*(-tmp_20*tmp_25 + tmp_23*tmp_40);
      real_t a_0_3 = tmp_39*(-tmp_19*tmp_25 + tmp_22*tmp_40);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_full_stokesvar_0_1_affine_q1::Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_mu( Point3D( {in_0, in_1, 0} ) );
   }

   void p1_full_stokesvar_0_1_affine_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

   void p1_full_stokesvar_0_2_affine_q1::integrateAll( const std::array< Point3D, 3 >& , Matrix< real_t, 3, 3 >&  ) const
   {
      
   }

   void p1_full_stokesvar_0_2_affine_q1::integrateRow0( const std::array< Point3D, 3 >& , Matrix< real_t, 1, 3 >&  ) const
   {
      
   }

   void p1_full_stokesvar_0_2_affine_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
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
      real_t tmp_27 = -tmp_11*tmp_3 + tmp_14*tmp_6;
      real_t tmp_28 = tmp_20*tmp_27;
      real_t tmp_29 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_30 = tmp_20*tmp_29;
      real_t tmp_31 = tmp_10*tmp_3 - tmp_12*tmp_14;
      real_t tmp_32 = tmp_20*tmp_31;
      real_t tmp_33 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_34 = tmp_33*(-tmp_28 - tmp_30 - tmp_32);
      real_t tmp_35 = -0.5*tmp_28 - 0.5*tmp_30 - 0.5*tmp_32;
      real_t tmp_36 = Scalar_Variable_Coefficient_3D_mu_out0_id0*(-0.5*tmp_21 - 0.5*tmp_23 - 0.5*tmp_25);
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
      real_t tmp_49 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_39 - p_affine_0_0*tmp_46 + p_affine_0_1*tmp_42 - p_affine_0_1*tmp_47 + p_affine_0_2*tmp_45 - p_affine_0_2*tmp_48 - p_affine_1_0*tmp_39 + p_affine_1_0*tmp_46 - p_affine_1_1*tmp_42 + p_affine_1_1*tmp_47 - p_affine_1_2*tmp_45 + p_affine_1_2*tmp_48 + p_affine_2_0*tmp_41 - p_affine_2_0*tmp_44 - p_affine_2_1*tmp_38 + p_affine_2_1*tmp_43 + p_affine_2_2*tmp_37 - p_affine_2_2*tmp_40 - p_affine_3_0*tmp_41 + p_affine_3_0*tmp_44 + p_affine_3_1*tmp_38 - p_affine_3_1*tmp_43 - p_affine_3_2*tmp_37 + p_affine_3_2*tmp_40);
      real_t tmp_50 = 2.0*tmp_36;
      real_t tmp_51 = 2.0*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_35;
      real_t tmp_52 = tmp_26*tmp_33;
      real_t tmp_53 = 1.0 / (tmp_19*tmp_19);
      real_t tmp_54 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_53;
      real_t tmp_55 = 0.33333333333333337*tmp_49*tmp_54;
      real_t tmp_56 = tmp_33*tmp_53;
      real_t tmp_57 = tmp_31*tmp_56;
      real_t tmp_58 = tmp_24*tmp_29;
      real_t tmp_59 = 1.0*tmp_54;
      real_t tmp_60 = tmp_27*tmp_59;
      real_t tmp_61 = tmp_31*tmp_59;
      real_t tmp_62 = tmp_29*tmp_8;
      real_t tmp_63 = tmp_27*tmp_56;
      real_t a_0_0 = tmp_49*(-tmp_26*tmp_34 + 4*tmp_35*tmp_36);
      real_t a_0_1 = tmp_49*(-tmp_25*tmp_34 + tmp_32*tmp_50);
      real_t a_0_2 = tmp_49*(-tmp_23*tmp_34 + tmp_30*tmp_50);
      real_t a_0_3 = tmp_49*(-tmp_21*tmp_34 + tmp_28*tmp_50);
      real_t a_1_0 = tmp_49*(tmp_25*tmp_51 - tmp_32*tmp_52);
      real_t a_1_1 = tmp_24*tmp_31*tmp_55;
      real_t a_1_2 = tmp_49*(-tmp_22*tmp_57 + tmp_58*tmp_59);
      real_t a_1_3 = tmp_49*(tmp_24*tmp_60 - tmp_57*tmp_8);
      real_t a_2_0 = tmp_49*(tmp_23*tmp_51 - tmp_30*tmp_52);
      real_t a_2_1 = tmp_49*(tmp_22*tmp_61 - tmp_56*tmp_58);
      real_t a_2_2 = tmp_22*tmp_29*tmp_55;
      real_t a_2_3 = tmp_49*(tmp_22*tmp_60 - tmp_56*tmp_62);
      real_t a_3_0 = tmp_49*(tmp_21*tmp_51 - tmp_28*tmp_52);
      real_t a_3_1 = tmp_49*(-tmp_24*tmp_63 + tmp_61*tmp_8);
      real_t a_3_2 = tmp_49*(-tmp_22*tmp_63 + tmp_59*tmp_62);
      real_t a_3_3 = tmp_27*tmp_55*tmp_8;
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

   void p1_full_stokesvar_0_2_affine_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
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
      real_t tmp_22 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_6);
      real_t tmp_23 = tmp_18*(tmp_10*tmp_11 - tmp_6*tmp_9);
      real_t tmp_24 = tmp_18*(-tmp_11*tmp_13 + tmp_3*tmp_9);
      real_t tmp_25 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0*(-tmp_22 - tmp_23 - tmp_24);
      real_t tmp_26 = Scalar_Variable_Coefficient_3D_mu_out0_id0*(-0.5*tmp_19 - 0.5*tmp_20 - 0.5*tmp_21);
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
      real_t tmp_39 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_29 - p_affine_0_0*tmp_36 + p_affine_0_1*tmp_32 - p_affine_0_1*tmp_37 + p_affine_0_2*tmp_35 - p_affine_0_2*tmp_38 - p_affine_1_0*tmp_29 + p_affine_1_0*tmp_36 - p_affine_1_1*tmp_32 + p_affine_1_1*tmp_37 - p_affine_1_2*tmp_35 + p_affine_1_2*tmp_38 + p_affine_2_0*tmp_31 - p_affine_2_0*tmp_34 - p_affine_2_1*tmp_28 + p_affine_2_1*tmp_33 + p_affine_2_2*tmp_27 - p_affine_2_2*tmp_30 - p_affine_3_0*tmp_31 + p_affine_3_0*tmp_34 + p_affine_3_1*tmp_28 - p_affine_3_1*tmp_33 - p_affine_3_2*tmp_27 + p_affine_3_2*tmp_30);
      real_t tmp_40 = 2.0*tmp_26;
      real_t a_0_0 = tmp_39*(-tmp_25*(-tmp_19 - tmp_20 - tmp_21) + 4*tmp_26*(-0.5*tmp_22 - 0.5*tmp_23 - 0.5*tmp_24));
      real_t a_0_1 = tmp_39*(-tmp_21*tmp_25 + tmp_24*tmp_40);
      real_t a_0_2 = tmp_39*(-tmp_20*tmp_25 + tmp_23*tmp_40);
      real_t a_0_3 = tmp_39*(-tmp_19*tmp_25 + tmp_22*tmp_40);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_full_stokesvar_0_2_affine_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

   void p1_full_stokesvar_1_0_affine_q1::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_0;
      real_t tmp_3 = p_affine_1_0 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = tmp_4 - (p_affine_1_1 + tmp_0)*(p_affine_2_0 + tmp_2);
      real_t tmp_6 = 1.0 / (tmp_5);
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_9 = tmp_6*tmp_8;
      real_t tmp_10 = -tmp_7 - tmp_9;
      real_t tmp_11 = tmp_3*tmp_6;
      real_t tmp_12 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_13 = tmp_12*tmp_6;
      real_t tmp_14 = (2.0/3.0)*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_15 = tmp_14*(-tmp_11 - tmp_13);
      real_t tmp_16 = -0.5*tmp_11 - 0.5*tmp_13;
      real_t tmp_17 = Scalar_Variable_Coefficient_2D_mu_out0_id0*(-0.5*tmp_7 - 0.5*tmp_9);
      real_t tmp_18 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_19 = 2.0*tmp_17;
      real_t tmp_20 = tmp_10*tmp_14;
      real_t tmp_21 = 2.0*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_16;
      real_t tmp_22 = 1.0 / (tmp_5*tmp_5);
      real_t tmp_23 = Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_22;
      real_t tmp_24 = 0.33333333333333337*tmp_18*tmp_23;
      real_t tmp_25 = 1.0*tmp_23;
      real_t tmp_26 = tmp_14*tmp_22;
      real_t tmp_27 = tmp_12*tmp_8;
      real_t a_0_0 = tmp_18*(-tmp_10*tmp_15 + 4*tmp_16*tmp_17);
      real_t a_0_1 = tmp_18*(tmp_13*tmp_19 - tmp_15*tmp_7);
      real_t a_0_2 = tmp_18*(tmp_11*tmp_19 - tmp_15*tmp_9);
      real_t a_1_0 = tmp_18*(-tmp_13*tmp_20 + tmp_21*tmp_7);
      real_t a_1_1 = tmp_1*tmp_12*tmp_24;
      real_t a_1_2 = tmp_18*(tmp_25*tmp_4 - tmp_26*tmp_27);
      real_t a_2_0 = tmp_18*(-tmp_11*tmp_20 + tmp_21*tmp_9);
      real_t a_2_1 = tmp_18*(tmp_25*tmp_27 - tmp_26*tmp_4);
      real_t a_2_2 = tmp_24*tmp_3*tmp_8;
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

   void p1_full_stokesvar_1_0_affine_q1::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_0;
      real_t tmp_3 = p_affine_1_0 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_0)*(p_affine_2_0 + tmp_2));
      real_t tmp_5 = tmp_1*tmp_4;
      real_t tmp_6 = tmp_4*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_7 = tmp_3*tmp_4;
      real_t tmp_8 = tmp_4*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_9 = (2.0/3.0)*Scalar_Variable_Coefficient_2D_mu_out0_id0*(-tmp_7 - tmp_8);
      real_t tmp_10 = Scalar_Variable_Coefficient_2D_mu_out0_id0*(-0.5*tmp_5 - 0.5*tmp_6);
      real_t tmp_11 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_12 = 2.0*tmp_10;
      real_t a_0_0 = tmp_11*(4*tmp_10*(-0.5*tmp_7 - 0.5*tmp_8) - tmp_9*(-tmp_5 - tmp_6));
      real_t a_0_1 = tmp_11*(tmp_12*tmp_8 - tmp_5*tmp_9);
      real_t a_0_2 = tmp_11*(tmp_12*tmp_7 - tmp_6*tmp_9);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_full_stokesvar_1_0_affine_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_1_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_2;
      real_t tmp_3 = p_affine_2_2 + tmp_2;
      real_t tmp_4 = p_affine_2_1 + tmp_0;
      real_t tmp_5 = p_affine_1_2 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_0;
      real_t tmp_8 = p_affine_1_0 + tmp_7;
      real_t tmp_9 = p_affine_3_2 + tmp_2;
      real_t tmp_10 = tmp_8*tmp_9;
      real_t tmp_11 = p_affine_3_1 + tmp_0;
      real_t tmp_12 = p_affine_2_0 + tmp_7;
      real_t tmp_13 = tmp_12*tmp_5;
      real_t tmp_14 = p_affine_3_0 + tmp_7;
      real_t tmp_15 = tmp_14*tmp_3;
      real_t tmp_16 = tmp_3*tmp_8;
      real_t tmp_17 = tmp_12*tmp_9;
      real_t tmp_18 = tmp_14*tmp_5;
      real_t tmp_19 = tmp_1*tmp_15 - tmp_1*tmp_17 + tmp_10*tmp_4 + tmp_11*tmp_13 - tmp_11*tmp_16 - tmp_18*tmp_4;
      real_t tmp_20 = 1.0 / (tmp_19);
      real_t tmp_21 = tmp_20*tmp_6;
      real_t tmp_22 = -tmp_1*tmp_9 + tmp_11*tmp_5;
      real_t tmp_23 = tmp_20*tmp_22;
      real_t tmp_24 = -tmp_11*tmp_3 + tmp_4*tmp_9;
      real_t tmp_25 = tmp_20*tmp_24;
      real_t tmp_26 = -tmp_21 - tmp_23 - tmp_25;
      real_t tmp_27 = tmp_13 - tmp_16;
      real_t tmp_28 = tmp_20*tmp_27;
      real_t tmp_29 = tmp_10 - tmp_18;
      real_t tmp_30 = tmp_20*tmp_29;
      real_t tmp_31 = tmp_15 - tmp_17;
      real_t tmp_32 = tmp_20*tmp_31;
      real_t tmp_33 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_34 = tmp_33*(-tmp_28 - tmp_30 - tmp_32);
      real_t tmp_35 = -0.5*tmp_28 - 0.5*tmp_30 - 0.5*tmp_32;
      real_t tmp_36 = Scalar_Variable_Coefficient_3D_mu_out0_id0*(-0.5*tmp_21 - 0.5*tmp_23 - 0.5*tmp_25);
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
      real_t tmp_49 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_39 - p_affine_0_0*tmp_46 + p_affine_0_1*tmp_42 - p_affine_0_1*tmp_47 + p_affine_0_2*tmp_45 - p_affine_0_2*tmp_48 - p_affine_1_0*tmp_39 + p_affine_1_0*tmp_46 - p_affine_1_1*tmp_42 + p_affine_1_1*tmp_47 - p_affine_1_2*tmp_45 + p_affine_1_2*tmp_48 + p_affine_2_0*tmp_41 - p_affine_2_0*tmp_44 - p_affine_2_1*tmp_38 + p_affine_2_1*tmp_43 + p_affine_2_2*tmp_37 - p_affine_2_2*tmp_40 - p_affine_3_0*tmp_41 + p_affine_3_0*tmp_44 + p_affine_3_1*tmp_38 - p_affine_3_1*tmp_43 - p_affine_3_2*tmp_37 + p_affine_3_2*tmp_40);
      real_t tmp_50 = 2.0*tmp_36;
      real_t tmp_51 = tmp_26*tmp_33;
      real_t tmp_52 = 2.0*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_35;
      real_t tmp_53 = 1.0 / (tmp_19*tmp_19);
      real_t tmp_54 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_53;
      real_t tmp_55 = 0.33333333333333337*tmp_49*tmp_54;
      real_t tmp_56 = tmp_24*tmp_29;
      real_t tmp_57 = 1.0*tmp_54;
      real_t tmp_58 = tmp_33*tmp_53;
      real_t tmp_59 = tmp_31*tmp_58;
      real_t tmp_60 = tmp_27*tmp_57;
      real_t tmp_61 = tmp_31*tmp_57;
      real_t tmp_62 = tmp_29*tmp_6;
      real_t tmp_63 = tmp_27*tmp_58;
      real_t a_0_0 = tmp_49*(-tmp_26*tmp_34 + 4*tmp_35*tmp_36);
      real_t a_0_1 = tmp_49*(-tmp_25*tmp_34 + tmp_32*tmp_50);
      real_t a_0_2 = tmp_49*(-tmp_23*tmp_34 + tmp_30*tmp_50);
      real_t a_0_3 = tmp_49*(-tmp_21*tmp_34 + tmp_28*tmp_50);
      real_t a_1_0 = tmp_49*(tmp_25*tmp_52 - tmp_32*tmp_51);
      real_t a_1_1 = tmp_24*tmp_31*tmp_55;
      real_t a_1_2 = tmp_49*(-tmp_22*tmp_59 + tmp_56*tmp_57);
      real_t a_1_3 = tmp_49*(tmp_24*tmp_60 - tmp_59*tmp_6);
      real_t a_2_0 = tmp_49*(tmp_23*tmp_52 - tmp_30*tmp_51);
      real_t a_2_1 = tmp_49*(tmp_22*tmp_61 - tmp_56*tmp_58);
      real_t a_2_2 = tmp_22*tmp_29*tmp_55;
      real_t a_2_3 = tmp_49*(tmp_22*tmp_60 - tmp_58*tmp_62);
      real_t a_3_0 = tmp_49*(tmp_21*tmp_52 - tmp_28*tmp_51);
      real_t a_3_1 = tmp_49*(-tmp_24*tmp_63 + tmp_6*tmp_61);
      real_t a_3_2 = tmp_49*(-tmp_22*tmp_63 + tmp_57*tmp_62);
      real_t a_3_3 = tmp_27*tmp_55*tmp_6;
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

   void p1_full_stokesvar_1_0_affine_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_1_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_2;
      real_t tmp_3 = p_affine_2_2 + tmp_2;
      real_t tmp_4 = p_affine_2_1 + tmp_0;
      real_t tmp_5 = p_affine_1_2 + tmp_2;
      real_t tmp_6 = -p_affine_0_0;
      real_t tmp_7 = p_affine_1_0 + tmp_6;
      real_t tmp_8 = p_affine_3_2 + tmp_2;
      real_t tmp_9 = tmp_7*tmp_8;
      real_t tmp_10 = p_affine_3_1 + tmp_0;
      real_t tmp_11 = p_affine_2_0 + tmp_6;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = p_affine_3_0 + tmp_6;
      real_t tmp_14 = tmp_13*tmp_3;
      real_t tmp_15 = tmp_3*tmp_7;
      real_t tmp_16 = tmp_11*tmp_8;
      real_t tmp_17 = tmp_13*tmp_5;
      real_t tmp_18 = 1.0 / (tmp_1*tmp_14 - tmp_1*tmp_16 + tmp_10*tmp_12 - tmp_10*tmp_15 - tmp_17*tmp_4 + tmp_4*tmp_9);
      real_t tmp_19 = tmp_18*(tmp_1*tmp_3 - tmp_4*tmp_5);
      real_t tmp_20 = tmp_18*(-tmp_1*tmp_8 + tmp_10*tmp_5);
      real_t tmp_21 = tmp_18*(-tmp_10*tmp_3 + tmp_4*tmp_8);
      real_t tmp_22 = tmp_18*(tmp_12 - tmp_15);
      real_t tmp_23 = tmp_18*(-tmp_17 + tmp_9);
      real_t tmp_24 = tmp_18*(tmp_14 - tmp_16);
      real_t tmp_25 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0*(-tmp_22 - tmp_23 - tmp_24);
      real_t tmp_26 = Scalar_Variable_Coefficient_3D_mu_out0_id0*(-0.5*tmp_19 - 0.5*tmp_20 - 0.5*tmp_21);
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
      real_t tmp_39 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_29 - p_affine_0_0*tmp_36 + p_affine_0_1*tmp_32 - p_affine_0_1*tmp_37 + p_affine_0_2*tmp_35 - p_affine_0_2*tmp_38 - p_affine_1_0*tmp_29 + p_affine_1_0*tmp_36 - p_affine_1_1*tmp_32 + p_affine_1_1*tmp_37 - p_affine_1_2*tmp_35 + p_affine_1_2*tmp_38 + p_affine_2_0*tmp_31 - p_affine_2_0*tmp_34 - p_affine_2_1*tmp_28 + p_affine_2_1*tmp_33 + p_affine_2_2*tmp_27 - p_affine_2_2*tmp_30 - p_affine_3_0*tmp_31 + p_affine_3_0*tmp_34 + p_affine_3_1*tmp_28 - p_affine_3_1*tmp_33 - p_affine_3_2*tmp_27 + p_affine_3_2*tmp_30);
      real_t tmp_40 = 2.0*tmp_26;
      real_t a_0_0 = tmp_39*(-tmp_25*(-tmp_19 - tmp_20 - tmp_21) + 4*tmp_26*(-0.5*tmp_22 - 0.5*tmp_23 - 0.5*tmp_24));
      real_t a_0_1 = tmp_39*(-tmp_21*tmp_25 + tmp_24*tmp_40);
      real_t a_0_2 = tmp_39*(-tmp_20*tmp_25 + tmp_23*tmp_40);
      real_t a_0_3 = tmp_39*(-tmp_19*tmp_25 + tmp_22*tmp_40);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_full_stokesvar_1_0_affine_q1::Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_mu( Point3D( {in_0, in_1, 0} ) );
   }

   void p1_full_stokesvar_1_0_affine_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

   void p1_full_stokesvar_1_1_affine_q1::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
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
      real_t tmp_10 = (2.0/3.0)*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_11 = -1.0*tmp_6 - 1.0*tmp_8;
      real_t tmp_12 = tmp_3*tmp_5;
      real_t tmp_13 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_14 = tmp_13*tmp_5;
      real_t tmp_15 = 1.0*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_16 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_17 = tmp_10*tmp_9;
      real_t tmp_18 = 2.0*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_19 = tmp_11*tmp_18;
      real_t tmp_20 = tmp_18*(-0.5*tmp_12 - 0.5*tmp_14);
      real_t tmp_21 = tmp_16*(tmp_12*tmp_20 - tmp_17*tmp_8 + tmp_19*tmp_8);
      real_t tmp_22 = tmp_16*(tmp_14*tmp_20 - tmp_17*tmp_6 + tmp_19*tmp_6);
      real_t tmp_23 = 1.0 / (tmp_4*tmp_4);
      real_t tmp_24 = 1.3333333333333335*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_23;
      real_t tmp_25 = tmp_15*tmp_23;
      real_t tmp_26 = tmp_16*(tmp_1*tmp_24*tmp_7 + tmp_13*tmp_25*tmp_3);
      real_t a_0_0 = tmp_16*(2*Scalar_Variable_Coefficient_2D_mu_out0_id0*(tmp_11*tmp_11) - tmp_10*(tmp_9*tmp_9) + tmp_15*((-tmp_12 - tmp_14)*(-tmp_12 - tmp_14)));
      real_t a_0_1 = tmp_21;
      real_t a_0_2 = tmp_22;
      real_t a_1_0 = tmp_21;
      real_t a_1_1 = tmp_16*(tmp_24*(tmp_7*tmp_7) + tmp_25*(tmp_3*tmp_3));
      real_t a_1_2 = tmp_26;
      real_t a_2_0 = tmp_22;
      real_t a_2_1 = tmp_26;
      real_t a_2_2 = tmp_16*((tmp_1*tmp_1)*tmp_24 + (tmp_13*tmp_13)*tmp_25);
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

   void p1_full_stokesvar_1_1_affine_q1::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_2D_mu( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = tmp_1*tmp_4;
      real_t tmp_6 = tmp_4*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_7 = -tmp_5 - tmp_6;
      real_t tmp_8 = (2.0/3.0)*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_9 = -1.0*tmp_5 - 1.0*tmp_6;
      real_t tmp_10 = tmp_3*tmp_4;
      real_t tmp_11 = tmp_4*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_12 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_13 = tmp_7*tmp_8;
      real_t tmp_14 = 2.0*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_15 = tmp_14*tmp_9;
      real_t tmp_16 = tmp_14*(-0.5*tmp_10 - 0.5*tmp_11);
      real_t a_0_0 = tmp_12*(2*Scalar_Variable_Coefficient_2D_mu_out0_id0*(tmp_9*tmp_9) + 1.0*Scalar_Variable_Coefficient_2D_mu_out0_id0*((-tmp_10 - tmp_11)*(-tmp_10 - tmp_11)) - (tmp_7*tmp_7)*tmp_8);
      real_t a_0_1 = tmp_12*(tmp_10*tmp_16 - tmp_13*tmp_6 + tmp_15*tmp_6);
      real_t a_0_2 = tmp_12*(tmp_11*tmp_16 - tmp_13*tmp_5 + tmp_15*tmp_5);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_full_stokesvar_1_1_affine_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = -p_affine_0_2;
      real_t tmp_8 = p_affine_3_2 + tmp_7;
      real_t tmp_9 = tmp_1*tmp_8;
      real_t tmp_10 = p_affine_3_1 + tmp_2;
      real_t tmp_11 = p_affine_1_2 + tmp_7;
      real_t tmp_12 = tmp_11*tmp_4;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_2_2 + tmp_7;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_1*tmp_14;
      real_t tmp_17 = tmp_4*tmp_8;
      real_t tmp_18 = tmp_11*tmp_13;
      real_t tmp_19 = tmp_10*tmp_12 - tmp_10*tmp_16 + tmp_15*tmp_5 - tmp_17*tmp_5 - tmp_18*tmp_3 + tmp_3*tmp_9;
      real_t tmp_20 = 1.0 / (tmp_19);
      real_t tmp_21 = tmp_20*tmp_6;
      real_t tmp_22 = -tmp_1*tmp_10 + tmp_13*tmp_5;
      real_t tmp_23 = tmp_20*tmp_22;
      real_t tmp_24 = tmp_10*tmp_4 - tmp_13*tmp_3;
      real_t tmp_25 = tmp_20*tmp_24;
      real_t tmp_26 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_27 = tmp_12 - tmp_16;
      real_t tmp_28 = tmp_20*tmp_27;
      real_t tmp_29 = -tmp_18 + tmp_9;
      real_t tmp_30 = tmp_20*tmp_29;
      real_t tmp_31 = tmp_15 - tmp_17;
      real_t tmp_32 = tmp_20*tmp_31;
      real_t tmp_33 = -tmp_28 - tmp_30 - tmp_32;
      real_t tmp_34 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_35 = -1.0*tmp_28 - 1.0*tmp_30 - 1.0*tmp_32;
      real_t tmp_36 = -tmp_11*tmp_3 + tmp_14*tmp_5;
      real_t tmp_37 = tmp_20*tmp_36;
      real_t tmp_38 = tmp_10*tmp_11 - tmp_5*tmp_8;
      real_t tmp_39 = tmp_20*tmp_38;
      real_t tmp_40 = -tmp_10*tmp_14 + tmp_3*tmp_8;
      real_t tmp_41 = tmp_20*tmp_40;
      real_t tmp_42 = p_affine_0_0*p_affine_1_1;
      real_t tmp_43 = p_affine_0_0*p_affine_1_2;
      real_t tmp_44 = p_affine_2_1*p_affine_3_2;
      real_t tmp_45 = p_affine_0_1*p_affine_1_0;
      real_t tmp_46 = p_affine_0_1*p_affine_1_2;
      real_t tmp_47 = p_affine_2_2*p_affine_3_0;
      real_t tmp_48 = p_affine_0_2*p_affine_1_0;
      real_t tmp_49 = p_affine_0_2*p_affine_1_1;
      real_t tmp_50 = p_affine_2_0*p_affine_3_1;
      real_t tmp_51 = p_affine_2_2*p_affine_3_1;
      real_t tmp_52 = p_affine_2_0*p_affine_3_2;
      real_t tmp_53 = p_affine_2_1*p_affine_3_0;
      real_t tmp_54 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_44 - p_affine_0_0*tmp_51 + p_affine_0_1*tmp_47 - p_affine_0_1*tmp_52 + p_affine_0_2*tmp_50 - p_affine_0_2*tmp_53 - p_affine_1_0*tmp_44 + p_affine_1_0*tmp_51 - p_affine_1_1*tmp_47 + p_affine_1_1*tmp_52 - p_affine_1_2*tmp_50 + p_affine_1_2*tmp_53 + p_affine_2_0*tmp_46 - p_affine_2_0*tmp_49 - p_affine_2_1*tmp_43 + p_affine_2_1*tmp_48 + p_affine_2_2*tmp_42 - p_affine_2_2*tmp_45 - p_affine_3_0*tmp_46 + p_affine_3_0*tmp_49 + p_affine_3_1*tmp_43 - p_affine_3_1*tmp_48 - p_affine_3_2*tmp_42 + p_affine_3_2*tmp_45);
      real_t tmp_55 = 2.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_56 = tmp_55*(-0.5*tmp_21 - 0.5*tmp_23 - 0.5*tmp_25);
      real_t tmp_57 = tmp_33*tmp_34;
      real_t tmp_58 = tmp_35*tmp_55;
      real_t tmp_59 = tmp_55*(-0.5*tmp_37 - 0.5*tmp_39 - 0.5*tmp_41);
      real_t tmp_60 = tmp_54*(tmp_25*tmp_56 - tmp_32*tmp_57 + tmp_32*tmp_58 + tmp_41*tmp_59);
      real_t tmp_61 = tmp_54*(tmp_23*tmp_56 - tmp_30*tmp_57 + tmp_30*tmp_58 + tmp_39*tmp_59);
      real_t tmp_62 = tmp_54*(tmp_21*tmp_56 - tmp_28*tmp_57 + tmp_28*tmp_58 + tmp_37*tmp_59);
      real_t tmp_63 = 1.0 / (tmp_19*tmp_19);
      real_t tmp_64 = tmp_26*tmp_63;
      real_t tmp_65 = 1.3333333333333335*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_63;
      real_t tmp_66 = tmp_24*tmp_64;
      real_t tmp_67 = tmp_31*tmp_65;
      real_t tmp_68 = tmp_40*tmp_64;
      real_t tmp_69 = tmp_54*(tmp_22*tmp_66 + tmp_29*tmp_67 + tmp_38*tmp_68);
      real_t tmp_70 = tmp_54*(tmp_27*tmp_67 + tmp_36*tmp_68 + tmp_6*tmp_66);
      real_t tmp_71 = tmp_54*(tmp_22*tmp_6*tmp_64 + tmp_27*tmp_29*tmp_65 + tmp_36*tmp_38*tmp_64);
      real_t a_0_0 = tmp_54*(2*Scalar_Variable_Coefficient_3D_mu_out0_id0*(tmp_35*tmp_35) + tmp_26*((-tmp_21 - tmp_23 - tmp_25)*(-tmp_21 - tmp_23 - tmp_25)) + tmp_26*((-tmp_37 - tmp_39 - tmp_41)*(-tmp_37 - tmp_39 - tmp_41)) - (tmp_33*tmp_33)*tmp_34);
      real_t a_0_1 = tmp_60;
      real_t a_0_2 = tmp_61;
      real_t a_0_3 = tmp_62;
      real_t a_1_0 = tmp_60;
      real_t a_1_1 = tmp_54*((tmp_24*tmp_24)*tmp_64 + (tmp_31*tmp_31)*tmp_65 + (tmp_40*tmp_40)*tmp_64);
      real_t a_1_2 = tmp_69;
      real_t a_1_3 = tmp_70;
      real_t a_2_0 = tmp_61;
      real_t a_2_1 = tmp_69;
      real_t a_2_2 = tmp_54*((tmp_22*tmp_22)*tmp_64 + (tmp_29*tmp_29)*tmp_65 + (tmp_38*tmp_38)*tmp_64);
      real_t a_2_3 = tmp_71;
      real_t a_3_0 = tmp_62;
      real_t a_3_1 = tmp_70;
      real_t a_3_2 = tmp_71;
      real_t a_3_3 = tmp_54*((tmp_27*tmp_27)*tmp_65 + (tmp_36*tmp_36)*tmp_64 + (tmp_6*tmp_6)*tmp_64);
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

   void p1_full_stokesvar_1_1_affine_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = -p_affine_0_2;
      real_t tmp_7 = p_affine_3_2 + tmp_6;
      real_t tmp_8 = tmp_1*tmp_7;
      real_t tmp_9 = p_affine_3_1 + tmp_2;
      real_t tmp_10 = p_affine_1_2 + tmp_6;
      real_t tmp_11 = tmp_10*tmp_4;
      real_t tmp_12 = p_affine_3_0 + tmp_0;
      real_t tmp_13 = p_affine_2_2 + tmp_6;
      real_t tmp_14 = tmp_12*tmp_13;
      real_t tmp_15 = tmp_1*tmp_13;
      real_t tmp_16 = tmp_4*tmp_7;
      real_t tmp_17 = tmp_10*tmp_12;
      real_t tmp_18 = 1.0 / (tmp_11*tmp_9 + tmp_14*tmp_5 - tmp_15*tmp_9 - tmp_16*tmp_5 - tmp_17*tmp_3 + tmp_3*tmp_8);
      real_t tmp_19 = tmp_18*(tmp_1*tmp_3 - tmp_4*tmp_5);
      real_t tmp_20 = tmp_18*(-tmp_1*tmp_9 + tmp_12*tmp_5);
      real_t tmp_21 = tmp_18*(-tmp_12*tmp_3 + tmp_4*tmp_9);
      real_t tmp_22 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_23 = tmp_18*(tmp_11 - tmp_15);
      real_t tmp_24 = tmp_18*(-tmp_17 + tmp_8);
      real_t tmp_25 = tmp_18*(tmp_14 - tmp_16);
      real_t tmp_26 = -tmp_23 - tmp_24 - tmp_25;
      real_t tmp_27 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_28 = -1.0*tmp_23 - 1.0*tmp_24 - 1.0*tmp_25;
      real_t tmp_29 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_5);
      real_t tmp_30 = tmp_18*(tmp_10*tmp_9 - tmp_5*tmp_7);
      real_t tmp_31 = tmp_18*(-tmp_13*tmp_9 + tmp_3*tmp_7);
      real_t tmp_32 = p_affine_0_0*p_affine_1_1;
      real_t tmp_33 = p_affine_0_0*p_affine_1_2;
      real_t tmp_34 = p_affine_2_1*p_affine_3_2;
      real_t tmp_35 = p_affine_0_1*p_affine_1_0;
      real_t tmp_36 = p_affine_0_1*p_affine_1_2;
      real_t tmp_37 = p_affine_2_2*p_affine_3_0;
      real_t tmp_38 = p_affine_0_2*p_affine_1_0;
      real_t tmp_39 = p_affine_0_2*p_affine_1_1;
      real_t tmp_40 = p_affine_2_0*p_affine_3_1;
      real_t tmp_41 = p_affine_2_2*p_affine_3_1;
      real_t tmp_42 = p_affine_2_0*p_affine_3_2;
      real_t tmp_43 = p_affine_2_1*p_affine_3_0;
      real_t tmp_44 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_34 - p_affine_0_0*tmp_41 + p_affine_0_1*tmp_37 - p_affine_0_1*tmp_42 + p_affine_0_2*tmp_40 - p_affine_0_2*tmp_43 - p_affine_1_0*tmp_34 + p_affine_1_0*tmp_41 - p_affine_1_1*tmp_37 + p_affine_1_1*tmp_42 - p_affine_1_2*tmp_40 + p_affine_1_2*tmp_43 + p_affine_2_0*tmp_36 - p_affine_2_0*tmp_39 - p_affine_2_1*tmp_33 + p_affine_2_1*tmp_38 + p_affine_2_2*tmp_32 - p_affine_2_2*tmp_35 - p_affine_3_0*tmp_36 + p_affine_3_0*tmp_39 + p_affine_3_1*tmp_33 - p_affine_3_1*tmp_38 - p_affine_3_2*tmp_32 + p_affine_3_2*tmp_35);
      real_t tmp_45 = 2.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_46 = tmp_45*(-0.5*tmp_19 - 0.5*tmp_20 - 0.5*tmp_21);
      real_t tmp_47 = tmp_26*tmp_27;
      real_t tmp_48 = tmp_28*tmp_45;
      real_t tmp_49 = tmp_45*(-0.5*tmp_29 - 0.5*tmp_30 - 0.5*tmp_31);
      real_t a_0_0 = tmp_44*(2*Scalar_Variable_Coefficient_3D_mu_out0_id0*(tmp_28*tmp_28) + tmp_22*((-tmp_19 - tmp_20 - tmp_21)*(-tmp_19 - tmp_20 - tmp_21)) + tmp_22*((-tmp_29 - tmp_30 - tmp_31)*(-tmp_29 - tmp_30 - tmp_31)) - (tmp_26*tmp_26)*tmp_27);
      real_t a_0_1 = tmp_44*(tmp_21*tmp_46 - tmp_25*tmp_47 + tmp_25*tmp_48 + tmp_31*tmp_49);
      real_t a_0_2 = tmp_44*(tmp_20*tmp_46 - tmp_24*tmp_47 + tmp_24*tmp_48 + tmp_30*tmp_49);
      real_t a_0_3 = tmp_44*(tmp_19*tmp_46 - tmp_23*tmp_47 + tmp_23*tmp_48 + tmp_29*tmp_49);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_full_stokesvar_1_1_affine_q1::Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_mu( Point3D( {in_0, in_1, 0} ) );
   }

   void p1_full_stokesvar_1_1_affine_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

   void p1_full_stokesvar_1_2_affine_q1::integrateAll( const std::array< Point3D, 3 >& , Matrix< real_t, 3, 3 >&  ) const
   {
      
   }

   void p1_full_stokesvar_1_2_affine_q1::integrateRow0( const std::array< Point3D, 3 >& , Matrix< real_t, 1, 3 >&  ) const
   {
      
   }

   void p1_full_stokesvar_1_2_affine_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
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
      real_t tmp_33 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_34 = tmp_33*(-tmp_28 - tmp_30 - tmp_32);
      real_t tmp_35 = -0.5*tmp_28 - 0.5*tmp_30 - 0.5*tmp_32;
      real_t tmp_36 = Scalar_Variable_Coefficient_3D_mu_out0_id0*(-0.5*tmp_21 - 0.5*tmp_23 - 0.5*tmp_25);
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
      real_t tmp_49 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_39 - p_affine_0_0*tmp_46 + p_affine_0_1*tmp_42 - p_affine_0_1*tmp_47 + p_affine_0_2*tmp_45 - p_affine_0_2*tmp_48 - p_affine_1_0*tmp_39 + p_affine_1_0*tmp_46 - p_affine_1_1*tmp_42 + p_affine_1_1*tmp_47 - p_affine_1_2*tmp_45 + p_affine_1_2*tmp_48 + p_affine_2_0*tmp_41 - p_affine_2_0*tmp_44 - p_affine_2_1*tmp_38 + p_affine_2_1*tmp_43 + p_affine_2_2*tmp_37 - p_affine_2_2*tmp_40 - p_affine_3_0*tmp_41 + p_affine_3_0*tmp_44 + p_affine_3_1*tmp_38 - p_affine_3_1*tmp_43 - p_affine_3_2*tmp_37 + p_affine_3_2*tmp_40);
      real_t tmp_50 = 2.0*tmp_36;
      real_t tmp_51 = 2.0*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_35;
      real_t tmp_52 = tmp_26*tmp_33;
      real_t tmp_53 = 1.0 / (tmp_19*tmp_19);
      real_t tmp_54 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_53;
      real_t tmp_55 = 0.33333333333333337*tmp_49*tmp_54;
      real_t tmp_56 = tmp_33*tmp_53;
      real_t tmp_57 = tmp_31*tmp_56;
      real_t tmp_58 = tmp_24*tmp_29;
      real_t tmp_59 = 1.0*tmp_54;
      real_t tmp_60 = tmp_27*tmp_59;
      real_t tmp_61 = tmp_31*tmp_59;
      real_t tmp_62 = tmp_29*tmp_8;
      real_t tmp_63 = tmp_27*tmp_56;
      real_t a_0_0 = tmp_49*(-tmp_26*tmp_34 + 4*tmp_35*tmp_36);
      real_t a_0_1 = tmp_49*(-tmp_25*tmp_34 + tmp_32*tmp_50);
      real_t a_0_2 = tmp_49*(-tmp_23*tmp_34 + tmp_30*tmp_50);
      real_t a_0_3 = tmp_49*(-tmp_21*tmp_34 + tmp_28*tmp_50);
      real_t a_1_0 = tmp_49*(tmp_25*tmp_51 - tmp_32*tmp_52);
      real_t a_1_1 = tmp_24*tmp_31*tmp_55;
      real_t a_1_2 = tmp_49*(-tmp_22*tmp_57 + tmp_58*tmp_59);
      real_t a_1_3 = tmp_49*(tmp_24*tmp_60 - tmp_57*tmp_8);
      real_t a_2_0 = tmp_49*(tmp_23*tmp_51 - tmp_30*tmp_52);
      real_t a_2_1 = tmp_49*(tmp_22*tmp_61 - tmp_56*tmp_58);
      real_t a_2_2 = tmp_22*tmp_29*tmp_55;
      real_t a_2_3 = tmp_49*(tmp_22*tmp_60 - tmp_56*tmp_62);
      real_t a_3_0 = tmp_49*(tmp_21*tmp_51 - tmp_28*tmp_52);
      real_t a_3_1 = tmp_49*(-tmp_24*tmp_63 + tmp_61*tmp_8);
      real_t a_3_2 = tmp_49*(-tmp_22*tmp_63 + tmp_59*tmp_62);
      real_t a_3_3 = tmp_27*tmp_55*tmp_8;
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

   void p1_full_stokesvar_1_2_affine_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
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
      real_t tmp_22 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_5);
      real_t tmp_23 = tmp_18*(tmp_1*tmp_9 - tmp_10*tmp_14);
      real_t tmp_24 = tmp_18*(tmp_13*tmp_14 - tmp_5*tmp_9);
      real_t tmp_25 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0*(-tmp_22 - tmp_23 - tmp_24);
      real_t tmp_26 = Scalar_Variable_Coefficient_3D_mu_out0_id0*(-0.5*tmp_19 - 0.5*tmp_20 - 0.5*tmp_21);
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
      real_t tmp_39 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_29 - p_affine_0_0*tmp_36 + p_affine_0_1*tmp_32 - p_affine_0_1*tmp_37 + p_affine_0_2*tmp_35 - p_affine_0_2*tmp_38 - p_affine_1_0*tmp_29 + p_affine_1_0*tmp_36 - p_affine_1_1*tmp_32 + p_affine_1_1*tmp_37 - p_affine_1_2*tmp_35 + p_affine_1_2*tmp_38 + p_affine_2_0*tmp_31 - p_affine_2_0*tmp_34 - p_affine_2_1*tmp_28 + p_affine_2_1*tmp_33 + p_affine_2_2*tmp_27 - p_affine_2_2*tmp_30 - p_affine_3_0*tmp_31 + p_affine_3_0*tmp_34 + p_affine_3_1*tmp_28 - p_affine_3_1*tmp_33 - p_affine_3_2*tmp_27 + p_affine_3_2*tmp_30);
      real_t tmp_40 = 2.0*tmp_26;
      real_t a_0_0 = tmp_39*(-tmp_25*(-tmp_19 - tmp_20 - tmp_21) + 4*tmp_26*(-0.5*tmp_22 - 0.5*tmp_23 - 0.5*tmp_24));
      real_t a_0_1 = tmp_39*(-tmp_21*tmp_25 + tmp_24*tmp_40);
      real_t a_0_2 = tmp_39*(-tmp_20*tmp_25 + tmp_23*tmp_40);
      real_t a_0_3 = tmp_39*(-tmp_19*tmp_25 + tmp_22*tmp_40);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_full_stokesvar_1_2_affine_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

   void p1_full_stokesvar_2_0_affine_q1::integrateAll( const std::array< Point3D, 3 >& , Matrix< real_t, 3, 3 >&  ) const
   {
      
   }

   void p1_full_stokesvar_2_0_affine_q1::integrateRow0( const std::array< Point3D, 3 >& , Matrix< real_t, 1, 3 >&  ) const
   {
      
   }

   void p1_full_stokesvar_2_0_affine_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_1_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_2;
      real_t tmp_3 = p_affine_2_2 + tmp_2;
      real_t tmp_4 = p_affine_2_1 + tmp_0;
      real_t tmp_5 = p_affine_1_2 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = p_affine_3_2 + tmp_2;
      real_t tmp_8 = -p_affine_0_0;
      real_t tmp_9 = p_affine_1_0 + tmp_8;
      real_t tmp_10 = tmp_4*tmp_9;
      real_t tmp_11 = p_affine_2_0 + tmp_8;
      real_t tmp_12 = p_affine_3_1 + tmp_0;
      real_t tmp_13 = tmp_11*tmp_12;
      real_t tmp_14 = p_affine_3_0 + tmp_8;
      real_t tmp_15 = tmp_1*tmp_14;
      real_t tmp_16 = tmp_12*tmp_9;
      real_t tmp_17 = tmp_1*tmp_11;
      real_t tmp_18 = tmp_14*tmp_4;
      real_t tmp_19 = tmp_10*tmp_7 + tmp_13*tmp_5 + tmp_15*tmp_3 - tmp_16*tmp_3 - tmp_17*tmp_7 - tmp_18*tmp_5;
      real_t tmp_20 = 1.0 / (tmp_19);
      real_t tmp_21 = tmp_20*tmp_6;
      real_t tmp_22 = -tmp_1*tmp_7 + tmp_12*tmp_5;
      real_t tmp_23 = tmp_20*tmp_22;
      real_t tmp_24 = -tmp_12*tmp_3 + tmp_4*tmp_7;
      real_t tmp_25 = tmp_20*tmp_24;
      real_t tmp_26 = -tmp_21 - tmp_23 - tmp_25;
      real_t tmp_27 = tmp_10 - tmp_17;
      real_t tmp_28 = tmp_20*tmp_27;
      real_t tmp_29 = tmp_15 - tmp_16;
      real_t tmp_30 = tmp_20*tmp_29;
      real_t tmp_31 = tmp_13 - tmp_18;
      real_t tmp_32 = tmp_20*tmp_31;
      real_t tmp_33 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_34 = tmp_33*(-tmp_28 - tmp_30 - tmp_32);
      real_t tmp_35 = -0.5*tmp_28 - 0.5*tmp_30 - 0.5*tmp_32;
      real_t tmp_36 = Scalar_Variable_Coefficient_3D_mu_out0_id0*(-0.5*tmp_21 - 0.5*tmp_23 - 0.5*tmp_25);
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
      real_t tmp_49 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_39 - p_affine_0_0*tmp_46 + p_affine_0_1*tmp_42 - p_affine_0_1*tmp_47 + p_affine_0_2*tmp_45 - p_affine_0_2*tmp_48 - p_affine_1_0*tmp_39 + p_affine_1_0*tmp_46 - p_affine_1_1*tmp_42 + p_affine_1_1*tmp_47 - p_affine_1_2*tmp_45 + p_affine_1_2*tmp_48 + p_affine_2_0*tmp_41 - p_affine_2_0*tmp_44 - p_affine_2_1*tmp_38 + p_affine_2_1*tmp_43 + p_affine_2_2*tmp_37 - p_affine_2_2*tmp_40 - p_affine_3_0*tmp_41 + p_affine_3_0*tmp_44 + p_affine_3_1*tmp_38 - p_affine_3_1*tmp_43 - p_affine_3_2*tmp_37 + p_affine_3_2*tmp_40);
      real_t tmp_50 = 2.0*tmp_36;
      real_t tmp_51 = tmp_26*tmp_33;
      real_t tmp_52 = 2.0*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_35;
      real_t tmp_53 = 1.0 / (tmp_19*tmp_19);
      real_t tmp_54 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_53;
      real_t tmp_55 = 0.33333333333333337*tmp_49*tmp_54;
      real_t tmp_56 = tmp_24*tmp_29;
      real_t tmp_57 = 1.0*tmp_54;
      real_t tmp_58 = tmp_33*tmp_53;
      real_t tmp_59 = tmp_31*tmp_58;
      real_t tmp_60 = tmp_27*tmp_57;
      real_t tmp_61 = tmp_31*tmp_57;
      real_t tmp_62 = tmp_29*tmp_6;
      real_t tmp_63 = tmp_27*tmp_58;
      real_t a_0_0 = tmp_49*(-tmp_26*tmp_34 + 4*tmp_35*tmp_36);
      real_t a_0_1 = tmp_49*(-tmp_25*tmp_34 + tmp_32*tmp_50);
      real_t a_0_2 = tmp_49*(-tmp_23*tmp_34 + tmp_30*tmp_50);
      real_t a_0_3 = tmp_49*(-tmp_21*tmp_34 + tmp_28*tmp_50);
      real_t a_1_0 = tmp_49*(tmp_25*tmp_52 - tmp_32*tmp_51);
      real_t a_1_1 = tmp_24*tmp_31*tmp_55;
      real_t a_1_2 = tmp_49*(-tmp_22*tmp_59 + tmp_56*tmp_57);
      real_t a_1_3 = tmp_49*(tmp_24*tmp_60 - tmp_59*tmp_6);
      real_t a_2_0 = tmp_49*(tmp_23*tmp_52 - tmp_30*tmp_51);
      real_t a_2_1 = tmp_49*(tmp_22*tmp_61 - tmp_56*tmp_58);
      real_t a_2_2 = tmp_22*tmp_29*tmp_55;
      real_t a_2_3 = tmp_49*(tmp_22*tmp_60 - tmp_58*tmp_62);
      real_t a_3_0 = tmp_49*(tmp_21*tmp_52 - tmp_28*tmp_51);
      real_t a_3_1 = tmp_49*(-tmp_24*tmp_63 + tmp_6*tmp_61);
      real_t a_3_2 = tmp_49*(-tmp_22*tmp_63 + tmp_57*tmp_62);
      real_t a_3_3 = tmp_27*tmp_55*tmp_6;
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

   void p1_full_stokesvar_2_0_affine_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_1_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_2;
      real_t tmp_3 = p_affine_2_2 + tmp_2;
      real_t tmp_4 = p_affine_2_1 + tmp_0;
      real_t tmp_5 = p_affine_1_2 + tmp_2;
      real_t tmp_6 = p_affine_3_2 + tmp_2;
      real_t tmp_7 = -p_affine_0_0;
      real_t tmp_8 = p_affine_1_0 + tmp_7;
      real_t tmp_9 = tmp_4*tmp_8;
      real_t tmp_10 = p_affine_2_0 + tmp_7;
      real_t tmp_11 = p_affine_3_1 + tmp_0;
      real_t tmp_12 = tmp_10*tmp_11;
      real_t tmp_13 = p_affine_3_0 + tmp_7;
      real_t tmp_14 = tmp_1*tmp_13;
      real_t tmp_15 = tmp_11*tmp_8;
      real_t tmp_16 = tmp_1*tmp_10;
      real_t tmp_17 = tmp_13*tmp_4;
      real_t tmp_18 = 1.0 / (tmp_12*tmp_5 + tmp_14*tmp_3 - tmp_15*tmp_3 - tmp_16*tmp_6 - tmp_17*tmp_5 + tmp_6*tmp_9);
      real_t tmp_19 = tmp_18*(tmp_1*tmp_3 - tmp_4*tmp_5);
      real_t tmp_20 = tmp_18*(-tmp_1*tmp_6 + tmp_11*tmp_5);
      real_t tmp_21 = tmp_18*(-tmp_11*tmp_3 + tmp_4*tmp_6);
      real_t tmp_22 = tmp_18*(-tmp_16 + tmp_9);
      real_t tmp_23 = tmp_18*(tmp_14 - tmp_15);
      real_t tmp_24 = tmp_18*(tmp_12 - tmp_17);
      real_t tmp_25 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0*(-tmp_22 - tmp_23 - tmp_24);
      real_t tmp_26 = Scalar_Variable_Coefficient_3D_mu_out0_id0*(-0.5*tmp_19 - 0.5*tmp_20 - 0.5*tmp_21);
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
      real_t tmp_39 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_29 - p_affine_0_0*tmp_36 + p_affine_0_1*tmp_32 - p_affine_0_1*tmp_37 + p_affine_0_2*tmp_35 - p_affine_0_2*tmp_38 - p_affine_1_0*tmp_29 + p_affine_1_0*tmp_36 - p_affine_1_1*tmp_32 + p_affine_1_1*tmp_37 - p_affine_1_2*tmp_35 + p_affine_1_2*tmp_38 + p_affine_2_0*tmp_31 - p_affine_2_0*tmp_34 - p_affine_2_1*tmp_28 + p_affine_2_1*tmp_33 + p_affine_2_2*tmp_27 - p_affine_2_2*tmp_30 - p_affine_3_0*tmp_31 + p_affine_3_0*tmp_34 + p_affine_3_1*tmp_28 - p_affine_3_1*tmp_33 - p_affine_3_2*tmp_27 + p_affine_3_2*tmp_30);
      real_t tmp_40 = 2.0*tmp_26;
      real_t a_0_0 = tmp_39*(-tmp_25*(-tmp_19 - tmp_20 - tmp_21) + 4*tmp_26*(-0.5*tmp_22 - 0.5*tmp_23 - 0.5*tmp_24));
      real_t a_0_1 = tmp_39*(-tmp_21*tmp_25 + tmp_24*tmp_40);
      real_t a_0_2 = tmp_39*(-tmp_20*tmp_25 + tmp_23*tmp_40);
      real_t a_0_3 = tmp_39*(-tmp_19*tmp_25 + tmp_22*tmp_40);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_full_stokesvar_2_0_affine_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

   void p1_full_stokesvar_2_1_affine_q1::integrateAll( const std::array< Point3D, 3 >& , Matrix< real_t, 3, 3 >&  ) const
   {
      
   }

   void p1_full_stokesvar_2_1_affine_q1::integrateRow0( const std::array< Point3D, 3 >& , Matrix< real_t, 1, 3 >&  ) const
   {
      
   }

   void p1_full_stokesvar_2_1_affine_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_2_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_2;
      real_t tmp_3 = p_affine_1_2 + tmp_2;
      real_t tmp_4 = p_affine_1_0 + tmp_0;
      real_t tmp_5 = p_affine_2_2 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = p_affine_3_2 + tmp_2;
      real_t tmp_8 = -p_affine_0_1;
      real_t tmp_9 = p_affine_2_1 + tmp_8;
      real_t tmp_10 = tmp_4*tmp_9;
      real_t tmp_11 = p_affine_3_1 + tmp_8;
      real_t tmp_12 = tmp_1*tmp_11;
      real_t tmp_13 = p_affine_3_0 + tmp_0;
      real_t tmp_14 = p_affine_1_1 + tmp_8;
      real_t tmp_15 = tmp_13*tmp_14;
      real_t tmp_16 = tmp_11*tmp_4;
      real_t tmp_17 = tmp_1*tmp_14;
      real_t tmp_18 = tmp_13*tmp_9;
      real_t tmp_19 = tmp_10*tmp_7 + tmp_12*tmp_3 + tmp_15*tmp_5 - tmp_16*tmp_5 - tmp_17*tmp_7 - tmp_18*tmp_3;
      real_t tmp_20 = 1.0 / (tmp_19);
      real_t tmp_21 = tmp_20*tmp_6;
      real_t tmp_22 = -tmp_13*tmp_3 + tmp_4*tmp_7;
      real_t tmp_23 = tmp_20*tmp_22;
      real_t tmp_24 = -tmp_1*tmp_7 + tmp_13*tmp_5;
      real_t tmp_25 = tmp_20*tmp_24;
      real_t tmp_26 = -tmp_21 - tmp_23 - tmp_25;
      real_t tmp_27 = tmp_10 - tmp_17;
      real_t tmp_28 = tmp_20*tmp_27;
      real_t tmp_29 = tmp_15 - tmp_16;
      real_t tmp_30 = tmp_20*tmp_29;
      real_t tmp_31 = tmp_12 - tmp_18;
      real_t tmp_32 = tmp_20*tmp_31;
      real_t tmp_33 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_34 = tmp_33*(-tmp_28 - tmp_30 - tmp_32);
      real_t tmp_35 = -0.5*tmp_28 - 0.5*tmp_30 - 0.5*tmp_32;
      real_t tmp_36 = Scalar_Variable_Coefficient_3D_mu_out0_id0*(-0.5*tmp_21 - 0.5*tmp_23 - 0.5*tmp_25);
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
      real_t tmp_49 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_39 - p_affine_0_0*tmp_46 + p_affine_0_1*tmp_42 - p_affine_0_1*tmp_47 + p_affine_0_2*tmp_45 - p_affine_0_2*tmp_48 - p_affine_1_0*tmp_39 + p_affine_1_0*tmp_46 - p_affine_1_1*tmp_42 + p_affine_1_1*tmp_47 - p_affine_1_2*tmp_45 + p_affine_1_2*tmp_48 + p_affine_2_0*tmp_41 - p_affine_2_0*tmp_44 - p_affine_2_1*tmp_38 + p_affine_2_1*tmp_43 + p_affine_2_2*tmp_37 - p_affine_2_2*tmp_40 - p_affine_3_0*tmp_41 + p_affine_3_0*tmp_44 + p_affine_3_1*tmp_38 - p_affine_3_1*tmp_43 - p_affine_3_2*tmp_37 + p_affine_3_2*tmp_40);
      real_t tmp_50 = 2.0*tmp_36;
      real_t tmp_51 = tmp_26*tmp_33;
      real_t tmp_52 = 2.0*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_35;
      real_t tmp_53 = 1.0 / (tmp_19*tmp_19);
      real_t tmp_54 = Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_53;
      real_t tmp_55 = 0.33333333333333337*tmp_49*tmp_54;
      real_t tmp_56 = tmp_24*tmp_29;
      real_t tmp_57 = 1.0*tmp_54;
      real_t tmp_58 = tmp_33*tmp_53;
      real_t tmp_59 = tmp_31*tmp_58;
      real_t tmp_60 = tmp_27*tmp_57;
      real_t tmp_61 = tmp_31*tmp_57;
      real_t tmp_62 = tmp_29*tmp_6;
      real_t tmp_63 = tmp_27*tmp_58;
      real_t a_0_0 = tmp_49*(-tmp_26*tmp_34 + 4*tmp_35*tmp_36);
      real_t a_0_1 = tmp_49*(-tmp_25*tmp_34 + tmp_32*tmp_50);
      real_t a_0_2 = tmp_49*(-tmp_23*tmp_34 + tmp_30*tmp_50);
      real_t a_0_3 = tmp_49*(-tmp_21*tmp_34 + tmp_28*tmp_50);
      real_t a_1_0 = tmp_49*(tmp_25*tmp_52 - tmp_32*tmp_51);
      real_t a_1_1 = tmp_24*tmp_31*tmp_55;
      real_t a_1_2 = tmp_49*(-tmp_22*tmp_59 + tmp_56*tmp_57);
      real_t a_1_3 = tmp_49*(tmp_24*tmp_60 - tmp_59*tmp_6);
      real_t a_2_0 = tmp_49*(tmp_23*tmp_52 - tmp_30*tmp_51);
      real_t a_2_1 = tmp_49*(tmp_22*tmp_61 - tmp_56*tmp_58);
      real_t a_2_2 = tmp_22*tmp_29*tmp_55;
      real_t a_2_3 = tmp_49*(tmp_22*tmp_60 - tmp_58*tmp_62);
      real_t a_3_0 = tmp_49*(tmp_21*tmp_52 - tmp_28*tmp_51);
      real_t a_3_1 = tmp_49*(-tmp_24*tmp_63 + tmp_6*tmp_61);
      real_t a_3_2 = tmp_49*(-tmp_22*tmp_63 + tmp_57*tmp_62);
      real_t a_3_3 = tmp_27*tmp_55*tmp_6;
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

   void p1_full_stokesvar_2_1_affine_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_2_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_2;
      real_t tmp_3 = p_affine_1_2 + tmp_2;
      real_t tmp_4 = p_affine_1_0 + tmp_0;
      real_t tmp_5 = p_affine_2_2 + tmp_2;
      real_t tmp_6 = p_affine_3_2 + tmp_2;
      real_t tmp_7 = -p_affine_0_1;
      real_t tmp_8 = p_affine_2_1 + tmp_7;
      real_t tmp_9 = tmp_4*tmp_8;
      real_t tmp_10 = p_affine_3_1 + tmp_7;
      real_t tmp_11 = tmp_1*tmp_10;
      real_t tmp_12 = p_affine_3_0 + tmp_0;
      real_t tmp_13 = p_affine_1_1 + tmp_7;
      real_t tmp_14 = tmp_12*tmp_13;
      real_t tmp_15 = tmp_10*tmp_4;
      real_t tmp_16 = tmp_1*tmp_13;
      real_t tmp_17 = tmp_12*tmp_8;
      real_t tmp_18 = 1.0 / (tmp_11*tmp_3 + tmp_14*tmp_5 - tmp_15*tmp_5 - tmp_16*tmp_6 - tmp_17*tmp_3 + tmp_6*tmp_9);
      real_t tmp_19 = tmp_18*(tmp_1*tmp_3 - tmp_4*tmp_5);
      real_t tmp_20 = tmp_18*(-tmp_12*tmp_3 + tmp_4*tmp_6);
      real_t tmp_21 = tmp_18*(-tmp_1*tmp_6 + tmp_12*tmp_5);
      real_t tmp_22 = tmp_18*(-tmp_16 + tmp_9);
      real_t tmp_23 = tmp_18*(tmp_14 - tmp_15);
      real_t tmp_24 = tmp_18*(tmp_11 - tmp_17);
      real_t tmp_25 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0*(-tmp_22 - tmp_23 - tmp_24);
      real_t tmp_26 = Scalar_Variable_Coefficient_3D_mu_out0_id0*(-0.5*tmp_19 - 0.5*tmp_20 - 0.5*tmp_21);
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
      real_t tmp_39 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_29 - p_affine_0_0*tmp_36 + p_affine_0_1*tmp_32 - p_affine_0_1*tmp_37 + p_affine_0_2*tmp_35 - p_affine_0_2*tmp_38 - p_affine_1_0*tmp_29 + p_affine_1_0*tmp_36 - p_affine_1_1*tmp_32 + p_affine_1_1*tmp_37 - p_affine_1_2*tmp_35 + p_affine_1_2*tmp_38 + p_affine_2_0*tmp_31 - p_affine_2_0*tmp_34 - p_affine_2_1*tmp_28 + p_affine_2_1*tmp_33 + p_affine_2_2*tmp_27 - p_affine_2_2*tmp_30 - p_affine_3_0*tmp_31 + p_affine_3_0*tmp_34 + p_affine_3_1*tmp_28 - p_affine_3_1*tmp_33 - p_affine_3_2*tmp_27 + p_affine_3_2*tmp_30);
      real_t tmp_40 = 2.0*tmp_26;
      real_t a_0_0 = tmp_39*(-tmp_25*(-tmp_19 - tmp_20 - tmp_21) + 4*tmp_26*(-0.5*tmp_22 - 0.5*tmp_23 - 0.5*tmp_24));
      real_t a_0_1 = tmp_39*(-tmp_21*tmp_25 + tmp_24*tmp_40);
      real_t a_0_2 = tmp_39*(-tmp_20*tmp_25 + tmp_23*tmp_40);
      real_t a_0_3 = tmp_39*(-tmp_19*tmp_25 + tmp_22*tmp_40);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_full_stokesvar_2_1_affine_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

   void p1_full_stokesvar_2_2_affine_q1::integrateAll( const std::array< Point3D, 3 >& , Matrix< real_t, 3, 3 >&  ) const
   {
      
   }

   void p1_full_stokesvar_2_2_affine_q1::integrateRow0( const std::array< Point3D, 3 >& , Matrix< real_t, 1, 3 >&  ) const
   {
      
   }

   void p1_full_stokesvar_2_2_affine_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
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
      real_t tmp_27 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_28 = -1.0*tmp_21 - 1.0*tmp_23 - 1.0*tmp_25;
      real_t tmp_29 = -tmp_1*tmp_14 + tmp_11*tmp_5;
      real_t tmp_30 = tmp_20*tmp_29;
      real_t tmp_31 = tmp_1*tmp_10 - tmp_11*tmp_15;
      real_t tmp_32 = tmp_20*tmp_31;
      real_t tmp_33 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_34 = tmp_20*tmp_33;
      real_t tmp_35 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_36 = -tmp_11*tmp_3 + tmp_14*tmp_6;
      real_t tmp_37 = tmp_20*tmp_36;
      real_t tmp_38 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_39 = tmp_20*tmp_38;
      real_t tmp_40 = tmp_10*tmp_3 - tmp_12*tmp_14;
      real_t tmp_41 = tmp_20*tmp_40;
      real_t tmp_42 = p_affine_0_0*p_affine_1_1;
      real_t tmp_43 = p_affine_0_0*p_affine_1_2;
      real_t tmp_44 = p_affine_2_1*p_affine_3_2;
      real_t tmp_45 = p_affine_0_1*p_affine_1_0;
      real_t tmp_46 = p_affine_0_1*p_affine_1_2;
      real_t tmp_47 = p_affine_2_2*p_affine_3_0;
      real_t tmp_48 = p_affine_0_2*p_affine_1_0;
      real_t tmp_49 = p_affine_0_2*p_affine_1_1;
      real_t tmp_50 = p_affine_2_0*p_affine_3_1;
      real_t tmp_51 = p_affine_2_2*p_affine_3_1;
      real_t tmp_52 = p_affine_2_0*p_affine_3_2;
      real_t tmp_53 = p_affine_2_1*p_affine_3_0;
      real_t tmp_54 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_44 - p_affine_0_0*tmp_51 + p_affine_0_1*tmp_47 - p_affine_0_1*tmp_52 + p_affine_0_2*tmp_50 - p_affine_0_2*tmp_53 - p_affine_1_0*tmp_44 + p_affine_1_0*tmp_51 - p_affine_1_1*tmp_47 + p_affine_1_1*tmp_52 - p_affine_1_2*tmp_50 + p_affine_1_2*tmp_53 + p_affine_2_0*tmp_46 - p_affine_2_0*tmp_49 - p_affine_2_1*tmp_43 + p_affine_2_1*tmp_48 + p_affine_2_2*tmp_42 - p_affine_2_2*tmp_45 - p_affine_3_0*tmp_46 + p_affine_3_0*tmp_49 + p_affine_3_1*tmp_43 - p_affine_3_1*tmp_48 - p_affine_3_2*tmp_42 + p_affine_3_2*tmp_45);
      real_t tmp_55 = tmp_26*tmp_27;
      real_t tmp_56 = 2.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_57 = tmp_28*tmp_56;
      real_t tmp_58 = tmp_56*(-0.5*tmp_30 - 0.5*tmp_32 - 0.5*tmp_34);
      real_t tmp_59 = tmp_56*(-0.5*tmp_37 - 0.5*tmp_39 - 0.5*tmp_41);
      real_t tmp_60 = tmp_54*(-tmp_25*tmp_55 + tmp_25*tmp_57 + tmp_34*tmp_58 + tmp_41*tmp_59);
      real_t tmp_61 = tmp_54*(-tmp_23*tmp_55 + tmp_23*tmp_57 + tmp_32*tmp_58 + tmp_39*tmp_59);
      real_t tmp_62 = tmp_54*(-tmp_21*tmp_55 + tmp_21*tmp_57 + tmp_30*tmp_58 + tmp_37*tmp_59);
      real_t tmp_63 = 1.0 / (tmp_19*tmp_19);
      real_t tmp_64 = 1.3333333333333335*Scalar_Variable_Coefficient_3D_mu_out0_id0*tmp_63;
      real_t tmp_65 = tmp_35*tmp_63;
      real_t tmp_66 = tmp_24*tmp_64;
      real_t tmp_67 = tmp_33*tmp_65;
      real_t tmp_68 = tmp_40*tmp_65;
      real_t tmp_69 = tmp_54*(tmp_22*tmp_66 + tmp_31*tmp_67 + tmp_38*tmp_68);
      real_t tmp_70 = tmp_54*(tmp_29*tmp_67 + tmp_36*tmp_68 + tmp_66*tmp_8);
      real_t tmp_71 = tmp_54*(tmp_22*tmp_64*tmp_8 + tmp_29*tmp_31*tmp_65 + tmp_36*tmp_38*tmp_65);
      real_t a_0_0 = tmp_54*(2*Scalar_Variable_Coefficient_3D_mu_out0_id0*(tmp_28*tmp_28) - (tmp_26*tmp_26)*tmp_27 + tmp_35*((-tmp_30 - tmp_32 - tmp_34)*(-tmp_30 - tmp_32 - tmp_34)) + tmp_35*((-tmp_37 - tmp_39 - tmp_41)*(-tmp_37 - tmp_39 - tmp_41)));
      real_t a_0_1 = tmp_60;
      real_t a_0_2 = tmp_61;
      real_t a_0_3 = tmp_62;
      real_t a_1_0 = tmp_60;
      real_t a_1_1 = tmp_54*((tmp_24*tmp_24)*tmp_64 + (tmp_33*tmp_33)*tmp_65 + (tmp_40*tmp_40)*tmp_65);
      real_t a_1_2 = tmp_69;
      real_t a_1_3 = tmp_70;
      real_t a_2_0 = tmp_61;
      real_t a_2_1 = tmp_69;
      real_t a_2_2 = tmp_54*((tmp_22*tmp_22)*tmp_64 + (tmp_31*tmp_31)*tmp_65 + (tmp_38*tmp_38)*tmp_65);
      real_t a_2_3 = tmp_71;
      real_t a_3_0 = tmp_62;
      real_t a_3_1 = tmp_70;
      real_t a_3_2 = tmp_71;
      real_t a_3_3 = tmp_54*((tmp_29*tmp_29)*tmp_65 + (tmp_36*tmp_36)*tmp_65 + tmp_64*(tmp_8*tmp_8));
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

   void p1_full_stokesvar_2_2_affine_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Scalar_Variable_Coefficient_3D_mu( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
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
      real_t tmp_23 = (2.0/3.0)*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_24 = -1.0*tmp_19 - 1.0*tmp_20 - 1.0*tmp_21;
      real_t tmp_25 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_5);
      real_t tmp_26 = tmp_18*(tmp_1*tmp_9 - tmp_10*tmp_14);
      real_t tmp_27 = tmp_18*(tmp_13*tmp_14 - tmp_5*tmp_9);
      real_t tmp_28 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_29 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_6);
      real_t tmp_30 = tmp_18*(tmp_10*tmp_11 - tmp_6*tmp_9);
      real_t tmp_31 = tmp_18*(-tmp_11*tmp_13 + tmp_3*tmp_9);
      real_t tmp_32 = p_affine_0_0*p_affine_1_1;
      real_t tmp_33 = p_affine_0_0*p_affine_1_2;
      real_t tmp_34 = p_affine_2_1*p_affine_3_2;
      real_t tmp_35 = p_affine_0_1*p_affine_1_0;
      real_t tmp_36 = p_affine_0_1*p_affine_1_2;
      real_t tmp_37 = p_affine_2_2*p_affine_3_0;
      real_t tmp_38 = p_affine_0_2*p_affine_1_0;
      real_t tmp_39 = p_affine_0_2*p_affine_1_1;
      real_t tmp_40 = p_affine_2_0*p_affine_3_1;
      real_t tmp_41 = p_affine_2_2*p_affine_3_1;
      real_t tmp_42 = p_affine_2_0*p_affine_3_2;
      real_t tmp_43 = p_affine_2_1*p_affine_3_0;
      real_t tmp_44 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_34 - p_affine_0_0*tmp_41 + p_affine_0_1*tmp_37 - p_affine_0_1*tmp_42 + p_affine_0_2*tmp_40 - p_affine_0_2*tmp_43 - p_affine_1_0*tmp_34 + p_affine_1_0*tmp_41 - p_affine_1_1*tmp_37 + p_affine_1_1*tmp_42 - p_affine_1_2*tmp_40 + p_affine_1_2*tmp_43 + p_affine_2_0*tmp_36 - p_affine_2_0*tmp_39 - p_affine_2_1*tmp_33 + p_affine_2_1*tmp_38 + p_affine_2_2*tmp_32 - p_affine_2_2*tmp_35 - p_affine_3_0*tmp_36 + p_affine_3_0*tmp_39 + p_affine_3_1*tmp_33 - p_affine_3_1*tmp_38 - p_affine_3_2*tmp_32 + p_affine_3_2*tmp_35);
      real_t tmp_45 = tmp_22*tmp_23;
      real_t tmp_46 = 2.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_47 = tmp_24*tmp_46;
      real_t tmp_48 = tmp_46*(-0.5*tmp_25 - 0.5*tmp_26 - 0.5*tmp_27);
      real_t tmp_49 = tmp_46*(-0.5*tmp_29 - 0.5*tmp_30 - 0.5*tmp_31);
      real_t a_0_0 = tmp_44*(2*Scalar_Variable_Coefficient_3D_mu_out0_id0*(tmp_24*tmp_24) - (tmp_22*tmp_22)*tmp_23 + tmp_28*((-tmp_25 - tmp_26 - tmp_27)*(-tmp_25 - tmp_26 - tmp_27)) + tmp_28*((-tmp_29 - tmp_30 - tmp_31)*(-tmp_29 - tmp_30 - tmp_31)));
      real_t a_0_1 = tmp_44*(-tmp_21*tmp_45 + tmp_21*tmp_47 + tmp_27*tmp_48 + tmp_31*tmp_49);
      real_t a_0_2 = tmp_44*(-tmp_20*tmp_45 + tmp_20*tmp_47 + tmp_26*tmp_48 + tmp_30*tmp_49);
      real_t a_0_3 = tmp_44*(-tmp_19*tmp_45 + tmp_19*tmp_47 + tmp_25*tmp_48 + tmp_29*tmp_49);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_full_stokesvar_2_2_affine_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

} // namespace forms
} // namespace hyteg

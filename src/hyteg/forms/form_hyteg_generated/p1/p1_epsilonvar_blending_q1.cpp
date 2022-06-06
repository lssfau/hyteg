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
 * Software:
 *
 * - quadpy version: 0.16.5
 *
 * Avoid modifying this file. If buggy, consider fixing the generator itself.
 */

#include "p1_epsilonvar_blending_q1.hpp"

namespace hyteg {
namespace forms {

   void p1_epsilonvar_0_0_blending_q1::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Blending_DF_Triangle_blend_out0_id0 = 0;
      real_t Blending_DF_Triangle_blend_out1_id0 = 0;
      real_t Blending_DF_Triangle_blend_out2_id0 = 0;
      real_t Blending_DF_Triangle_blend_out3_id0 = 0;
      real_t Blending_F_Triangle_blend_out0_id1 = 0;
      real_t Blending_F_Triangle_blend_out1_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      Blending_DF_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_DF_Triangle_blend_out0_id0, &Blending_DF_Triangle_blend_out1_id0, &Blending_DF_Triangle_blend_out2_id0, &Blending_DF_Triangle_blend_out3_id0 );
      Blending_F_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_F_Triangle_blend_out0_id1, &Blending_F_Triangle_blend_out1_id1 );
      Scalar_Variable_Coefficient_2D_mu( Blending_F_Triangle_blend_out0_id1, Blending_F_Triangle_blend_out1_id1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = Blending_DF_Triangle_blend_out0_id0*Blending_DF_Triangle_blend_out3_id0 - Blending_DF_Triangle_blend_out1_id0*Blending_DF_Triangle_blend_out2_id0;
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = p_affine_1_0 + tmp_3;
      real_t tmp_5 = 1/(tmp_2*(tmp_1*tmp_4 - (p_affine_1_1 + tmp_0)*(p_affine_2_0 + tmp_3)));
      real_t tmp_6 = Blending_DF_Triangle_blend_out1_id0*tmp_5;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_9 = tmp_6*tmp_8;
      real_t tmp_10 = Blending_DF_Triangle_blend_out0_id0*tmp_5;
      real_t tmp_11 = tmp_10*tmp_4;
      real_t tmp_12 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_13 = tmp_10*tmp_12;
      real_t tmp_14 = 1.0*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_15 = 1.0*tmp_5;
      real_t tmp_16 = Blending_DF_Triangle_blend_out2_id0*tmp_15;
      real_t tmp_17 = tmp_16*tmp_4;
      real_t tmp_18 = tmp_12*tmp_16;
      real_t tmp_19 = Blending_DF_Triangle_blend_out3_id0*tmp_15;
      real_t tmp_20 = tmp_1*tmp_19;
      real_t tmp_21 = tmp_19*tmp_8;
      real_t tmp_22 = tmp_17 + tmp_18 - tmp_20 - tmp_21;
      real_t tmp_23 = 2*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_24 = 0.5*std::abs(tmp_2)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_25 = 0.5*tmp_13;
      real_t tmp_26 = 0.5*tmp_7;
      real_t tmp_27 = tmp_25 - tmp_26;
      real_t tmp_28 = 0.5*tmp_11;
      real_t tmp_29 = 0.5*tmp_9;
      real_t tmp_30 = 4*Scalar_Variable_Coefficient_2D_mu_out0_id0*(-tmp_25 + tmp_26 - tmp_28 + tmp_29);
      real_t tmp_31 = -tmp_18 + tmp_20;
      real_t tmp_32 = tmp_22*tmp_23;
      real_t tmp_33 = tmp_24*(tmp_27*tmp_30 + tmp_31*tmp_32);
      real_t tmp_34 = tmp_28 - tmp_29;
      real_t tmp_35 = -tmp_17 + tmp_21;
      real_t tmp_36 = tmp_24*(tmp_30*tmp_34 + tmp_32*tmp_35);
      real_t tmp_37 = tmp_24*(4*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_27*tmp_34 + tmp_23*tmp_31*tmp_35);
      real_t a_0_0 = tmp_24*(tmp_14*((-tmp_11 - tmp_13 + tmp_7 + tmp_9)*(-tmp_11 - tmp_13 + tmp_7 + tmp_9)) + (tmp_22*tmp_22)*tmp_23);
      real_t a_0_1 = tmp_33;
      real_t a_0_2 = tmp_36;
      real_t a_1_0 = tmp_33;
      real_t a_1_1 = tmp_24*(tmp_14*((tmp_13 - tmp_7)*(tmp_13 - tmp_7)) + tmp_23*(tmp_31*tmp_31));
      real_t a_1_2 = tmp_37;
      real_t a_2_0 = tmp_36;
      real_t a_2_1 = tmp_37;
      real_t a_2_2 = tmp_24*(tmp_14*((tmp_11 - tmp_9)*(tmp_11 - tmp_9)) + tmp_23*(tmp_35*tmp_35));
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

   void p1_epsilonvar_0_0_blending_q1::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Blending_DF_Triangle_blend_out0_id0 = 0;
      real_t Blending_DF_Triangle_blend_out1_id0 = 0;
      real_t Blending_DF_Triangle_blend_out2_id0 = 0;
      real_t Blending_DF_Triangle_blend_out3_id0 = 0;
      real_t Blending_F_Triangle_blend_out0_id1 = 0;
      real_t Blending_F_Triangle_blend_out1_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      Blending_DF_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_DF_Triangle_blend_out0_id0, &Blending_DF_Triangle_blend_out1_id0, &Blending_DF_Triangle_blend_out2_id0, &Blending_DF_Triangle_blend_out3_id0 );
      Blending_F_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_F_Triangle_blend_out0_id1, &Blending_F_Triangle_blend_out1_id1 );
      Scalar_Variable_Coefficient_2D_mu( Blending_F_Triangle_blend_out0_id1, Blending_F_Triangle_blend_out1_id1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = Blending_DF_Triangle_blend_out0_id0*Blending_DF_Triangle_blend_out3_id0 - Blending_DF_Triangle_blend_out1_id0*Blending_DF_Triangle_blend_out2_id0;
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = p_affine_1_0 + tmp_3;
      real_t tmp_5 = 1/(tmp_2*(tmp_1*tmp_4 - (p_affine_1_1 + tmp_0)*(p_affine_2_0 + tmp_3)));
      real_t tmp_6 = Blending_DF_Triangle_blend_out1_id0*tmp_5;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_9 = tmp_6*tmp_8;
      real_t tmp_10 = Blending_DF_Triangle_blend_out0_id0*tmp_5;
      real_t tmp_11 = tmp_10*tmp_4;
      real_t tmp_12 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_13 = tmp_10*tmp_12;
      real_t tmp_14 = 1.0*tmp_5;
      real_t tmp_15 = Blending_DF_Triangle_blend_out2_id0*tmp_14;
      real_t tmp_16 = tmp_15*tmp_4;
      real_t tmp_17 = tmp_12*tmp_15;
      real_t tmp_18 = Blending_DF_Triangle_blend_out3_id0*tmp_14;
      real_t tmp_19 = tmp_1*tmp_18;
      real_t tmp_20 = tmp_18*tmp_8;
      real_t tmp_21 = tmp_16 + tmp_17 - tmp_19 - tmp_20;
      real_t tmp_22 = 2*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_23 = 0.5*std::abs(tmp_2)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_24 = 0.5*tmp_13;
      real_t tmp_25 = 0.5*tmp_7;
      real_t tmp_26 = 0.5*tmp_11;
      real_t tmp_27 = 0.5*tmp_9;
      real_t tmp_28 = 4*Scalar_Variable_Coefficient_2D_mu_out0_id0*(-tmp_24 + tmp_25 - tmp_26 + tmp_27);
      real_t tmp_29 = tmp_21*tmp_22;
      real_t a_0_0 = tmp_23*(1.0*Scalar_Variable_Coefficient_2D_mu_out0_id0*((-tmp_11 - tmp_13 + tmp_7 + tmp_9)*(-tmp_11 - tmp_13 + tmp_7 + tmp_9)) + (tmp_21*tmp_21)*tmp_22);
      real_t a_0_1 = tmp_23*(tmp_28*(tmp_24 - tmp_25) + tmp_29*(-tmp_17 + tmp_19));
      real_t a_0_2 = tmp_23*(tmp_28*(tmp_26 - tmp_27) + tmp_29*(-tmp_16 + tmp_20));
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_epsilonvar_0_0_blending_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_14 = Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_15 = -Blending_DF_Tetrahedron_blend_out0_id0*tmp_12 + Blending_DF_Tetrahedron_blend_out0_id0*tmp_9 + Blending_DF_Tetrahedron_blend_out1_id0*tmp_10 - Blending_DF_Tetrahedron_blend_out1_id0*tmp_13 + Blending_DF_Tetrahedron_blend_out2_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out2_id0*tmp_14;
      real_t tmp_16 = -p_affine_0_2;
      real_t tmp_17 = p_affine_3_2 + tmp_16;
      real_t tmp_18 = p_affine_1_2 + tmp_16;
      real_t tmp_19 = p_affine_3_1 + tmp_2;
      real_t tmp_20 = tmp_19*tmp_5;
      real_t tmp_21 = p_affine_2_2 + tmp_16;
      real_t tmp_22 = p_affine_3_0 + tmp_0;
      real_t tmp_23 = tmp_22*tmp_6;
      real_t tmp_24 = tmp_1*tmp_19;
      real_t tmp_25 = tmp_22*tmp_3;
      real_t tmp_26 = 1/(tmp_15*(tmp_17*tmp_4 - tmp_17*tmp_7 + tmp_18*tmp_20 - tmp_18*tmp_25 + tmp_21*tmp_23 - tmp_21*tmp_24));
      real_t tmp_27 = tmp_26*(Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0 - Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0);
      real_t tmp_28 = tmp_27*tmp_8;
      real_t tmp_29 = tmp_23 - tmp_24;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = tmp_20 - tmp_25;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = -tmp_1*tmp_21 + tmp_18*tmp_5;
      real_t tmp_34 = tmp_26*(-Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0 + Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0);
      real_t tmp_35 = tmp_33*tmp_34;
      real_t tmp_36 = tmp_1*tmp_17 - tmp_18*tmp_22;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = -tmp_17*tmp_5 + tmp_21*tmp_22;
      real_t tmp_39 = tmp_34*tmp_38;
      real_t tmp_40 = -tmp_18*tmp_3 + tmp_21*tmp_6;
      real_t tmp_41 = tmp_26*(Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0 - Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0);
      real_t tmp_42 = tmp_40*tmp_41;
      real_t tmp_43 = -tmp_17*tmp_6 + tmp_18*tmp_19;
      real_t tmp_44 = tmp_41*tmp_43;
      real_t tmp_45 = tmp_17*tmp_3 - tmp_19*tmp_21;
      real_t tmp_46 = tmp_41*tmp_45;
      real_t tmp_47 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_48 = tmp_26*(-Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out7_id0 + Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_49 = tmp_48*tmp_8;
      real_t tmp_50 = tmp_29*tmp_48;
      real_t tmp_51 = tmp_31*tmp_48;
      real_t tmp_52 = tmp_26*(Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_53 = tmp_33*tmp_52;
      real_t tmp_54 = tmp_36*tmp_52;
      real_t tmp_55 = tmp_38*tmp_52;
      real_t tmp_56 = tmp_26*(-Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out7_id0);
      real_t tmp_57 = tmp_40*tmp_56;
      real_t tmp_58 = tmp_43*tmp_56;
      real_t tmp_59 = tmp_45*tmp_56;
      real_t tmp_60 = 1.0*tmp_26;
      real_t tmp_61 = tmp_60*(tmp_11 - tmp_14);
      real_t tmp_62 = tmp_61*tmp_8;
      real_t tmp_63 = tmp_29*tmp_61;
      real_t tmp_64 = tmp_31*tmp_61;
      real_t tmp_65 = tmp_60*(tmp_10 - tmp_13);
      real_t tmp_66 = tmp_33*tmp_65;
      real_t tmp_67 = tmp_36*tmp_65;
      real_t tmp_68 = tmp_38*tmp_65;
      real_t tmp_69 = tmp_60*(-tmp_12 + tmp_9);
      real_t tmp_70 = tmp_40*tmp_69;
      real_t tmp_71 = tmp_43*tmp_69;
      real_t tmp_72 = tmp_45*tmp_69;
      real_t tmp_73 = -tmp_62 - tmp_63 - tmp_64 - tmp_66 - tmp_67 - tmp_68 - tmp_70 - tmp_71 - tmp_72;
      real_t tmp_74 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_75 = p_affine_0_0*p_affine_1_1;
      real_t tmp_76 = p_affine_0_0*p_affine_1_2;
      real_t tmp_77 = p_affine_2_1*p_affine_3_2;
      real_t tmp_78 = p_affine_0_1*p_affine_1_0;
      real_t tmp_79 = p_affine_0_1*p_affine_1_2;
      real_t tmp_80 = p_affine_2_2*p_affine_3_0;
      real_t tmp_81 = p_affine_0_2*p_affine_1_0;
      real_t tmp_82 = p_affine_0_2*p_affine_1_1;
      real_t tmp_83 = p_affine_2_0*p_affine_3_1;
      real_t tmp_84 = p_affine_2_2*p_affine_3_1;
      real_t tmp_85 = p_affine_2_0*p_affine_3_2;
      real_t tmp_86 = p_affine_2_1*p_affine_3_0;
      real_t tmp_87 = 0.16666666666666663*std::abs(tmp_15)*std::abs(p_affine_0_0*tmp_77 - p_affine_0_0*tmp_84 + p_affine_0_1*tmp_80 - p_affine_0_1*tmp_85 + p_affine_0_2*tmp_83 - p_affine_0_2*tmp_86 - p_affine_1_0*tmp_77 + p_affine_1_0*tmp_84 - p_affine_1_1*tmp_80 + p_affine_1_1*tmp_85 - p_affine_1_2*tmp_83 + p_affine_1_2*tmp_86 + p_affine_2_0*tmp_79 - p_affine_2_0*tmp_82 - p_affine_2_1*tmp_76 + p_affine_2_1*tmp_81 + p_affine_2_2*tmp_75 - p_affine_2_2*tmp_78 - p_affine_3_0*tmp_79 + p_affine_3_0*tmp_82 + p_affine_3_1*tmp_76 - p_affine_3_1*tmp_81 - p_affine_3_2*tmp_75 + p_affine_3_2*tmp_78);
      real_t tmp_88 = 0.5*tmp_32;
      real_t tmp_89 = 0.5*tmp_39;
      real_t tmp_90 = 0.5*tmp_46;
      real_t tmp_91 = tmp_88 + tmp_89 + tmp_90;
      real_t tmp_92 = 0.5*tmp_28;
      real_t tmp_93 = 0.5*tmp_30;
      real_t tmp_94 = 0.5*tmp_35;
      real_t tmp_95 = 0.5*tmp_37;
      real_t tmp_96 = 0.5*tmp_42;
      real_t tmp_97 = 0.5*tmp_44;
      real_t tmp_98 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_99 = tmp_98*(-tmp_88 - tmp_89 - tmp_90 - tmp_92 - tmp_93 - tmp_94 - tmp_95 - tmp_96 - tmp_97);
      real_t tmp_100 = 0.5*tmp_51;
      real_t tmp_101 = 0.5*tmp_55;
      real_t tmp_102 = 0.5*tmp_59;
      real_t tmp_103 = tmp_100 + tmp_101 + tmp_102;
      real_t tmp_104 = 0.5*tmp_49;
      real_t tmp_105 = 0.5*tmp_50;
      real_t tmp_106 = 0.5*tmp_53;
      real_t tmp_107 = 0.5*tmp_54;
      real_t tmp_108 = 0.5*tmp_57;
      real_t tmp_109 = 0.5*tmp_58;
      real_t tmp_110 = tmp_98*(-tmp_100 - tmp_101 - tmp_102 - tmp_104 - tmp_105 - tmp_106 - tmp_107 - tmp_108 - tmp_109);
      real_t tmp_111 = tmp_64 + tmp_68 + tmp_72;
      real_t tmp_112 = tmp_73*tmp_74;
      real_t tmp_113 = tmp_87*(tmp_103*tmp_110 + tmp_111*tmp_112 + tmp_91*tmp_99);
      real_t tmp_114 = tmp_93 + tmp_95 + tmp_97;
      real_t tmp_115 = tmp_105 + tmp_107 + tmp_109;
      real_t tmp_116 = tmp_63 + tmp_67 + tmp_71;
      real_t tmp_117 = tmp_87*(tmp_110*tmp_115 + tmp_112*tmp_116 + tmp_114*tmp_99);
      real_t tmp_118 = tmp_92 + tmp_94 + tmp_96;
      real_t tmp_119 = tmp_104 + tmp_106 + tmp_108;
      real_t tmp_120 = tmp_62 + tmp_66 + tmp_70;
      real_t tmp_121 = tmp_87*(tmp_110*tmp_119 + tmp_112*tmp_120 + tmp_118*tmp_99);
      real_t tmp_122 = tmp_91*tmp_98;
      real_t tmp_123 = tmp_103*tmp_98;
      real_t tmp_124 = tmp_111*tmp_74;
      real_t tmp_125 = tmp_87*(tmp_114*tmp_122 + tmp_115*tmp_123 + tmp_116*tmp_124);
      real_t tmp_126 = tmp_87*(tmp_118*tmp_122 + tmp_119*tmp_123 + tmp_120*tmp_124);
      real_t tmp_127 = tmp_87*(tmp_114*tmp_118*tmp_98 + tmp_115*tmp_119*tmp_98 + tmp_116*tmp_120*tmp_74);
      real_t a_0_0 = tmp_87*(tmp_47*((-tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46)*(-tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46)) + tmp_47*((-tmp_49 - tmp_50 - tmp_51 - tmp_53 - tmp_54 - tmp_55 - tmp_57 - tmp_58 - tmp_59)*(-tmp_49 - tmp_50 - tmp_51 - tmp_53 - tmp_54 - tmp_55 - tmp_57 - tmp_58 - tmp_59)) + (tmp_73*tmp_73)*tmp_74);
      real_t a_0_1 = tmp_113;
      real_t a_0_2 = tmp_117;
      real_t a_0_3 = tmp_121;
      real_t a_1_0 = tmp_113;
      real_t a_1_1 = tmp_87*((tmp_111*tmp_111)*tmp_74 + tmp_47*((tmp_32 + tmp_39 + tmp_46)*(tmp_32 + tmp_39 + tmp_46)) + tmp_47*((tmp_51 + tmp_55 + tmp_59)*(tmp_51 + tmp_55 + tmp_59)));
      real_t a_1_2 = tmp_125;
      real_t a_1_3 = tmp_126;
      real_t a_2_0 = tmp_117;
      real_t a_2_1 = tmp_125;
      real_t a_2_2 = tmp_87*((tmp_116*tmp_116)*tmp_74 + tmp_47*((tmp_30 + tmp_37 + tmp_44)*(tmp_30 + tmp_37 + tmp_44)) + tmp_47*((tmp_50 + tmp_54 + tmp_58)*(tmp_50 + tmp_54 + tmp_58)));
      real_t a_2_3 = tmp_127;
      real_t a_3_0 = tmp_121;
      real_t a_3_1 = tmp_126;
      real_t a_3_2 = tmp_127;
      real_t a_3_3 = tmp_87*((tmp_120*tmp_120)*tmp_74 + tmp_47*((tmp_28 + tmp_35 + tmp_42)*(tmp_28 + tmp_35 + tmp_42)) + tmp_47*((tmp_49 + tmp_53 + tmp_57)*(tmp_49 + tmp_53 + tmp_57)));
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

   void p1_epsilonvar_0_0_blending_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_14 = Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_15 = -Blending_DF_Tetrahedron_blend_out0_id0*tmp_12 + Blending_DF_Tetrahedron_blend_out0_id0*tmp_9 + Blending_DF_Tetrahedron_blend_out1_id0*tmp_10 - Blending_DF_Tetrahedron_blend_out1_id0*tmp_13 + Blending_DF_Tetrahedron_blend_out2_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out2_id0*tmp_14;
      real_t tmp_16 = -p_affine_0_2;
      real_t tmp_17 = p_affine_3_2 + tmp_16;
      real_t tmp_18 = p_affine_1_2 + tmp_16;
      real_t tmp_19 = p_affine_3_1 + tmp_2;
      real_t tmp_20 = tmp_19*tmp_5;
      real_t tmp_21 = p_affine_2_2 + tmp_16;
      real_t tmp_22 = p_affine_3_0 + tmp_0;
      real_t tmp_23 = tmp_22*tmp_6;
      real_t tmp_24 = tmp_1*tmp_19;
      real_t tmp_25 = tmp_22*tmp_3;
      real_t tmp_26 = 1/(tmp_15*(tmp_17*tmp_4 - tmp_17*tmp_7 + tmp_18*tmp_20 - tmp_18*tmp_25 + tmp_21*tmp_23 - tmp_21*tmp_24));
      real_t tmp_27 = tmp_26*(Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0 - Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0);
      real_t tmp_28 = tmp_27*tmp_8;
      real_t tmp_29 = tmp_23 - tmp_24;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = tmp_20 - tmp_25;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = -tmp_1*tmp_21 + tmp_18*tmp_5;
      real_t tmp_34 = tmp_26*(-Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0 + Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0);
      real_t tmp_35 = tmp_33*tmp_34;
      real_t tmp_36 = tmp_1*tmp_17 - tmp_18*tmp_22;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = -tmp_17*tmp_5 + tmp_21*tmp_22;
      real_t tmp_39 = tmp_34*tmp_38;
      real_t tmp_40 = -tmp_18*tmp_3 + tmp_21*tmp_6;
      real_t tmp_41 = tmp_26*(Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0 - Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0);
      real_t tmp_42 = tmp_40*tmp_41;
      real_t tmp_43 = -tmp_17*tmp_6 + tmp_18*tmp_19;
      real_t tmp_44 = tmp_41*tmp_43;
      real_t tmp_45 = tmp_17*tmp_3 - tmp_19*tmp_21;
      real_t tmp_46 = tmp_41*tmp_45;
      real_t tmp_47 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_48 = tmp_26*(-Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out7_id0 + Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_49 = tmp_48*tmp_8;
      real_t tmp_50 = tmp_29*tmp_48;
      real_t tmp_51 = tmp_31*tmp_48;
      real_t tmp_52 = tmp_26*(Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_53 = tmp_33*tmp_52;
      real_t tmp_54 = tmp_36*tmp_52;
      real_t tmp_55 = tmp_38*tmp_52;
      real_t tmp_56 = tmp_26*(-Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out7_id0);
      real_t tmp_57 = tmp_40*tmp_56;
      real_t tmp_58 = tmp_43*tmp_56;
      real_t tmp_59 = tmp_45*tmp_56;
      real_t tmp_60 = 1.0*tmp_26;
      real_t tmp_61 = tmp_60*(tmp_11 - tmp_14);
      real_t tmp_62 = tmp_61*tmp_8;
      real_t tmp_63 = tmp_29*tmp_61;
      real_t tmp_64 = tmp_31*tmp_61;
      real_t tmp_65 = tmp_60*(tmp_10 - tmp_13);
      real_t tmp_66 = tmp_33*tmp_65;
      real_t tmp_67 = tmp_36*tmp_65;
      real_t tmp_68 = tmp_38*tmp_65;
      real_t tmp_69 = tmp_60*(-tmp_12 + tmp_9);
      real_t tmp_70 = tmp_40*tmp_69;
      real_t tmp_71 = tmp_43*tmp_69;
      real_t tmp_72 = tmp_45*tmp_69;
      real_t tmp_73 = -tmp_62 - tmp_63 - tmp_64 - tmp_66 - tmp_67 - tmp_68 - tmp_70 - tmp_71 - tmp_72;
      real_t tmp_74 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_75 = p_affine_0_0*p_affine_1_1;
      real_t tmp_76 = p_affine_0_0*p_affine_1_2;
      real_t tmp_77 = p_affine_2_1*p_affine_3_2;
      real_t tmp_78 = p_affine_0_1*p_affine_1_0;
      real_t tmp_79 = p_affine_0_1*p_affine_1_2;
      real_t tmp_80 = p_affine_2_2*p_affine_3_0;
      real_t tmp_81 = p_affine_0_2*p_affine_1_0;
      real_t tmp_82 = p_affine_0_2*p_affine_1_1;
      real_t tmp_83 = p_affine_2_0*p_affine_3_1;
      real_t tmp_84 = p_affine_2_2*p_affine_3_1;
      real_t tmp_85 = p_affine_2_0*p_affine_3_2;
      real_t tmp_86 = p_affine_2_1*p_affine_3_0;
      real_t tmp_87 = 0.16666666666666663*std::abs(tmp_15)*std::abs(p_affine_0_0*tmp_77 - p_affine_0_0*tmp_84 + p_affine_0_1*tmp_80 - p_affine_0_1*tmp_85 + p_affine_0_2*tmp_83 - p_affine_0_2*tmp_86 - p_affine_1_0*tmp_77 + p_affine_1_0*tmp_84 - p_affine_1_1*tmp_80 + p_affine_1_1*tmp_85 - p_affine_1_2*tmp_83 + p_affine_1_2*tmp_86 + p_affine_2_0*tmp_79 - p_affine_2_0*tmp_82 - p_affine_2_1*tmp_76 + p_affine_2_1*tmp_81 + p_affine_2_2*tmp_75 - p_affine_2_2*tmp_78 - p_affine_3_0*tmp_79 + p_affine_3_0*tmp_82 + p_affine_3_1*tmp_76 - p_affine_3_1*tmp_81 - p_affine_3_2*tmp_75 + p_affine_3_2*tmp_78);
      real_t tmp_88 = 0.5*tmp_32;
      real_t tmp_89 = 0.5*tmp_39;
      real_t tmp_90 = 0.5*tmp_46;
      real_t tmp_91 = 0.5*tmp_28;
      real_t tmp_92 = 0.5*tmp_30;
      real_t tmp_93 = 0.5*tmp_35;
      real_t tmp_94 = 0.5*tmp_37;
      real_t tmp_95 = 0.5*tmp_42;
      real_t tmp_96 = 0.5*tmp_44;
      real_t tmp_97 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_98 = tmp_97*(-tmp_88 - tmp_89 - tmp_90 - tmp_91 - tmp_92 - tmp_93 - tmp_94 - tmp_95 - tmp_96);
      real_t tmp_99 = 0.5*tmp_51;
      real_t tmp_100 = 0.5*tmp_55;
      real_t tmp_101 = 0.5*tmp_59;
      real_t tmp_102 = 0.5*tmp_49;
      real_t tmp_103 = 0.5*tmp_50;
      real_t tmp_104 = 0.5*tmp_53;
      real_t tmp_105 = 0.5*tmp_54;
      real_t tmp_106 = 0.5*tmp_57;
      real_t tmp_107 = 0.5*tmp_58;
      real_t tmp_108 = tmp_97*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_104 - tmp_105 - tmp_106 - tmp_107 - tmp_99);
      real_t tmp_109 = tmp_73*tmp_74;
      real_t a_0_0 = tmp_87*(tmp_47*((-tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46)*(-tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46)) + tmp_47*((-tmp_49 - tmp_50 - tmp_51 - tmp_53 - tmp_54 - tmp_55 - tmp_57 - tmp_58 - tmp_59)*(-tmp_49 - tmp_50 - tmp_51 - tmp_53 - tmp_54 - tmp_55 - tmp_57 - tmp_58 - tmp_59)) + (tmp_73*tmp_73)*tmp_74);
      real_t a_0_1 = tmp_87*(tmp_108*(tmp_100 + tmp_101 + tmp_99) + tmp_109*(tmp_64 + tmp_68 + tmp_72) + tmp_98*(tmp_88 + tmp_89 + tmp_90));
      real_t a_0_2 = tmp_87*(tmp_108*(tmp_103 + tmp_105 + tmp_107) + tmp_109*(tmp_63 + tmp_67 + tmp_71) + tmp_98*(tmp_92 + tmp_94 + tmp_96));
      real_t a_0_3 = tmp_87*(tmp_108*(tmp_102 + tmp_104 + tmp_106) + tmp_109*(tmp_62 + tmp_66 + tmp_70) + tmp_98*(tmp_91 + tmp_93 + tmp_95));
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_epsilonvar_0_0_blending_q1::Blending_DF_Triangle_blend( real_t in_0, real_t in_1, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3 ) const
   {
      Point3D  mappedPt( {in_0, in_1, 0} );
      Matrix2r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 1, 0 );
      *out_3 = DPsi( 1, 1 );
   }

   void p1_epsilonvar_0_0_blending_q1::Blending_F_Triangle_blend( real_t in_0, real_t in_1, real_t * out_0, real_t * out_1 ) const
   {
      Point3D  in( {in_0, in_1, 0} );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
   }

   void p1_epsilonvar_0_0_blending_q1::Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_mu( Point3D( {in_0, in_1, 0} ) );
   }

   void p1_epsilonvar_0_0_blending_q1::Blending_DF_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
   {
      Point3D  mappedPt( {in_0, in_1, in_2} );
      Matrix3r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 0, 2 );
      *out_3 = DPsi( 1, 0 );
      *out_4 = DPsi( 1, 1 );
      *out_5 = DPsi( 1, 2 );
      *out_6 = DPsi( 2, 0 );
      *out_7 = DPsi( 2, 1 );
      *out_8 = DPsi( 2, 2 );
   }

   void p1_epsilonvar_0_0_blending_q1::Blending_F_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2 ) const
   {
      Point3D  in( {in_0, in_1, in_2} );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
      *out_2 = out[2];
   }

   void p1_epsilonvar_0_0_blending_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

   void p1_epsilonvar_0_1_blending_q1::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Blending_DF_Triangle_blend_out0_id0 = 0;
      real_t Blending_DF_Triangle_blend_out1_id0 = 0;
      real_t Blending_DF_Triangle_blend_out2_id0 = 0;
      real_t Blending_DF_Triangle_blend_out3_id0 = 0;
      real_t Blending_F_Triangle_blend_out0_id1 = 0;
      real_t Blending_F_Triangle_blend_out1_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      Blending_DF_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_DF_Triangle_blend_out0_id0, &Blending_DF_Triangle_blend_out1_id0, &Blending_DF_Triangle_blend_out2_id0, &Blending_DF_Triangle_blend_out3_id0 );
      Blending_F_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_F_Triangle_blend_out0_id1, &Blending_F_Triangle_blend_out1_id1 );
      Scalar_Variable_Coefficient_2D_mu( Blending_F_Triangle_blend_out0_id1, Blending_F_Triangle_blend_out1_id1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = Blending_DF_Triangle_blend_out0_id0*Blending_DF_Triangle_blend_out3_id0 - Blending_DF_Triangle_blend_out1_id0*Blending_DF_Triangle_blend_out2_id0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = 0.5/(tmp_2*(tmp_1*tmp_4 - (p_affine_1_1 + tmp_3)*(p_affine_2_0 + tmp_0)));
      real_t tmp_6 = Blending_DF_Triangle_blend_out2_id0*tmp_5;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_9 = tmp_6*tmp_8;
      real_t tmp_10 = Blending_DF_Triangle_blend_out3_id0*tmp_5;
      real_t tmp_11 = tmp_10*tmp_4;
      real_t tmp_12 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_13 = tmp_10*tmp_12;
      real_t tmp_14 = -tmp_11 - tmp_13 + tmp_7 + tmp_9;
      real_t tmp_15 = Blending_DF_Triangle_blend_out0_id0*tmp_5;
      real_t tmp_16 = tmp_1*tmp_15;
      real_t tmp_17 = tmp_15*tmp_8;
      real_t tmp_18 = Blending_DF_Triangle_blend_out1_id0*tmp_5;
      real_t tmp_19 = tmp_18*tmp_4;
      real_t tmp_20 = tmp_12*tmp_18;
      real_t tmp_21 = 2.0*Scalar_Variable_Coefficient_2D_mu_out0_id0*std::abs(tmp_2)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_22 = tmp_21*(-tmp_16 - tmp_17 + tmp_19 + tmp_20);
      real_t tmp_23 = tmp_11 - tmp_9;
      real_t tmp_24 = tmp_13 - tmp_7;
      real_t tmp_25 = tmp_21*(tmp_17 - tmp_19);
      real_t tmp_26 = tmp_21*(tmp_16 - tmp_20);
      real_t a_0_0 = tmp_14*tmp_22;
      real_t a_0_1 = tmp_22*tmp_23;
      real_t a_0_2 = tmp_22*tmp_24;
      real_t a_1_0 = tmp_14*tmp_25;
      real_t a_1_1 = tmp_23*tmp_25;
      real_t a_1_2 = tmp_24*tmp_25;
      real_t a_2_0 = tmp_14*tmp_26;
      real_t a_2_1 = tmp_23*tmp_26;
      real_t a_2_2 = tmp_24*tmp_26;
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

   void p1_epsilonvar_0_1_blending_q1::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Blending_DF_Triangle_blend_out0_id0 = 0;
      real_t Blending_DF_Triangle_blend_out1_id0 = 0;
      real_t Blending_DF_Triangle_blend_out2_id0 = 0;
      real_t Blending_DF_Triangle_blend_out3_id0 = 0;
      real_t Blending_F_Triangle_blend_out0_id1 = 0;
      real_t Blending_F_Triangle_blend_out1_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      Blending_DF_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_DF_Triangle_blend_out0_id0, &Blending_DF_Triangle_blend_out1_id0, &Blending_DF_Triangle_blend_out2_id0, &Blending_DF_Triangle_blend_out3_id0 );
      Blending_F_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_F_Triangle_blend_out0_id1, &Blending_F_Triangle_blend_out1_id1 );
      Scalar_Variable_Coefficient_2D_mu( Blending_F_Triangle_blend_out0_id1, Blending_F_Triangle_blend_out1_id1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = Blending_DF_Triangle_blend_out0_id0*Blending_DF_Triangle_blend_out3_id0 - Blending_DF_Triangle_blend_out1_id0*Blending_DF_Triangle_blend_out2_id0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = 0.5/(tmp_2*(tmp_1*tmp_4 - (p_affine_1_1 + tmp_3)*(p_affine_2_0 + tmp_0)));
      real_t tmp_6 = Blending_DF_Triangle_blend_out2_id0*tmp_5;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_9 = tmp_6*tmp_8;
      real_t tmp_10 = Blending_DF_Triangle_blend_out3_id0*tmp_5;
      real_t tmp_11 = tmp_10*tmp_4;
      real_t tmp_12 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_13 = tmp_10*tmp_12;
      real_t tmp_14 = Blending_DF_Triangle_blend_out0_id0*tmp_5;
      real_t tmp_15 = Blending_DF_Triangle_blend_out1_id0*tmp_5;
      real_t tmp_16 = 2.0*Scalar_Variable_Coefficient_2D_mu_out0_id0*(-tmp_1*tmp_14 + tmp_12*tmp_15 - tmp_14*tmp_8 + tmp_15*tmp_4)*std::abs(tmp_2)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = tmp_16*(-tmp_11 - tmp_13 + tmp_7 + tmp_9);
      real_t a_0_1 = tmp_16*(tmp_11 - tmp_9);
      real_t a_0_2 = tmp_16*(tmp_13 - tmp_7);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_epsilonvar_0_1_blending_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_8 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_13 = -Blending_DF_Tetrahedron_blend_out3_id0*tmp_11 + Blending_DF_Tetrahedron_blend_out3_id0*tmp_9 - Blending_DF_Tetrahedron_blend_out4_id0*tmp_12 + Blending_DF_Tetrahedron_blend_out4_id0*tmp_7 - Blending_DF_Tetrahedron_blend_out5_id0*tmp_10 + Blending_DF_Tetrahedron_blend_out5_id0*tmp_8;
      real_t tmp_14 = -p_affine_0_2;
      real_t tmp_15 = p_affine_3_2 + tmp_14;
      real_t tmp_16 = tmp_1*tmp_15;
      real_t tmp_17 = p_affine_3_1 + tmp_2;
      real_t tmp_18 = p_affine_1_2 + tmp_14;
      real_t tmp_19 = tmp_18*tmp_4;
      real_t tmp_20 = p_affine_3_0 + tmp_0;
      real_t tmp_21 = p_affine_2_2 + tmp_14;
      real_t tmp_22 = tmp_20*tmp_21;
      real_t tmp_23 = tmp_1*tmp_21;
      real_t tmp_24 = tmp_15*tmp_4;
      real_t tmp_25 = tmp_18*tmp_20;
      real_t tmp_26 = 0.5/(tmp_13*(tmp_16*tmp_3 + tmp_17*tmp_19 - tmp_17*tmp_23 + tmp_22*tmp_5 - tmp_24*tmp_5 - tmp_25*tmp_3));
      real_t tmp_27 = tmp_26*(Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_28 = tmp_27*tmp_6;
      real_t tmp_29 = -tmp_1*tmp_17 + tmp_20*tmp_5;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = tmp_17*tmp_4 - tmp_20*tmp_3;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = tmp_19 - tmp_23;
      real_t tmp_34 = tmp_26*(-Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_35 = tmp_33*tmp_34;
      real_t tmp_36 = tmp_16 - tmp_25;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = tmp_22 - tmp_24;
      real_t tmp_39 = tmp_34*tmp_38;
      real_t tmp_40 = -tmp_18*tmp_3 + tmp_21*tmp_5;
      real_t tmp_41 = tmp_26*(Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0);
      real_t tmp_42 = tmp_40*tmp_41;
      real_t tmp_43 = -tmp_15*tmp_5 + tmp_17*tmp_18;
      real_t tmp_44 = tmp_41*tmp_43;
      real_t tmp_45 = tmp_15*tmp_3 - tmp_17*tmp_21;
      real_t tmp_46 = tmp_41*tmp_45;
      real_t tmp_47 = -tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46;
      real_t tmp_48 = tmp_26*(-tmp_10 + tmp_8);
      real_t tmp_49 = tmp_48*tmp_6;
      real_t tmp_50 = tmp_29*tmp_48;
      real_t tmp_51 = tmp_31*tmp_48;
      real_t tmp_52 = tmp_26*(-tmp_12 + tmp_7);
      real_t tmp_53 = tmp_33*tmp_52;
      real_t tmp_54 = tmp_36*tmp_52;
      real_t tmp_55 = tmp_38*tmp_52;
      real_t tmp_56 = tmp_26*(-tmp_11 + tmp_9);
      real_t tmp_57 = tmp_40*tmp_56;
      real_t tmp_58 = tmp_43*tmp_56;
      real_t tmp_59 = tmp_45*tmp_56;
      real_t tmp_60 = p_affine_0_0*p_affine_1_1;
      real_t tmp_61 = p_affine_0_0*p_affine_1_2;
      real_t tmp_62 = p_affine_2_1*p_affine_3_2;
      real_t tmp_63 = p_affine_0_1*p_affine_1_0;
      real_t tmp_64 = p_affine_0_1*p_affine_1_2;
      real_t tmp_65 = p_affine_2_2*p_affine_3_0;
      real_t tmp_66 = p_affine_0_2*p_affine_1_0;
      real_t tmp_67 = p_affine_0_2*p_affine_1_1;
      real_t tmp_68 = p_affine_2_0*p_affine_3_1;
      real_t tmp_69 = p_affine_2_2*p_affine_3_1;
      real_t tmp_70 = p_affine_2_0*p_affine_3_2;
      real_t tmp_71 = p_affine_2_1*p_affine_3_0;
      real_t tmp_72 = 0.66666666666666652*Scalar_Variable_Coefficient_3D_mu_out0_id0*std::abs(tmp_13)*std::abs(p_affine_0_0*tmp_62 - p_affine_0_0*tmp_69 + p_affine_0_1*tmp_65 - p_affine_0_1*tmp_70 + p_affine_0_2*tmp_68 - p_affine_0_2*tmp_71 - p_affine_1_0*tmp_62 + p_affine_1_0*tmp_69 - p_affine_1_1*tmp_65 + p_affine_1_1*tmp_70 - p_affine_1_2*tmp_68 + p_affine_1_2*tmp_71 + p_affine_2_0*tmp_64 - p_affine_2_0*tmp_67 - p_affine_2_1*tmp_61 + p_affine_2_1*tmp_66 + p_affine_2_2*tmp_60 - p_affine_2_2*tmp_63 - p_affine_3_0*tmp_64 + p_affine_3_0*tmp_67 + p_affine_3_1*tmp_61 - p_affine_3_1*tmp_66 - p_affine_3_2*tmp_60 + p_affine_3_2*tmp_63);
      real_t tmp_73 = tmp_72*(-tmp_49 - tmp_50 - tmp_51 - tmp_53 - tmp_54 - tmp_55 - tmp_57 - tmp_58 - tmp_59);
      real_t tmp_74 = tmp_32 + tmp_39 + tmp_46;
      real_t tmp_75 = tmp_30 + tmp_37 + tmp_44;
      real_t tmp_76 = tmp_28 + tmp_35 + tmp_42;
      real_t tmp_77 = tmp_72*(tmp_51 + tmp_55 + tmp_59);
      real_t tmp_78 = tmp_72*(tmp_50 + tmp_54 + tmp_58);
      real_t tmp_79 = tmp_72*(tmp_49 + tmp_53 + tmp_57);
      real_t a_0_0 = tmp_47*tmp_73;
      real_t a_0_1 = tmp_73*tmp_74;
      real_t a_0_2 = tmp_73*tmp_75;
      real_t a_0_3 = tmp_73*tmp_76;
      real_t a_1_0 = tmp_47*tmp_77;
      real_t a_1_1 = tmp_74*tmp_77;
      real_t a_1_2 = tmp_75*tmp_77;
      real_t a_1_3 = tmp_76*tmp_77;
      real_t a_2_0 = tmp_47*tmp_78;
      real_t a_2_1 = tmp_74*tmp_78;
      real_t a_2_2 = tmp_75*tmp_78;
      real_t a_2_3 = tmp_76*tmp_78;
      real_t a_3_0 = tmp_47*tmp_79;
      real_t a_3_1 = tmp_74*tmp_79;
      real_t a_3_2 = tmp_75*tmp_79;
      real_t a_3_3 = tmp_76*tmp_79;
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

   void p1_epsilonvar_0_1_blending_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_8 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_13 = -Blending_DF_Tetrahedron_blend_out3_id0*tmp_11 + Blending_DF_Tetrahedron_blend_out3_id0*tmp_9 - Blending_DF_Tetrahedron_blend_out4_id0*tmp_12 + Blending_DF_Tetrahedron_blend_out4_id0*tmp_7 - Blending_DF_Tetrahedron_blend_out5_id0*tmp_10 + Blending_DF_Tetrahedron_blend_out5_id0*tmp_8;
      real_t tmp_14 = -p_affine_0_2;
      real_t tmp_15 = p_affine_3_2 + tmp_14;
      real_t tmp_16 = tmp_1*tmp_15;
      real_t tmp_17 = p_affine_3_1 + tmp_2;
      real_t tmp_18 = p_affine_1_2 + tmp_14;
      real_t tmp_19 = tmp_18*tmp_4;
      real_t tmp_20 = p_affine_3_0 + tmp_0;
      real_t tmp_21 = p_affine_2_2 + tmp_14;
      real_t tmp_22 = tmp_20*tmp_21;
      real_t tmp_23 = tmp_1*tmp_21;
      real_t tmp_24 = tmp_15*tmp_4;
      real_t tmp_25 = tmp_18*tmp_20;
      real_t tmp_26 = 0.5/(tmp_13*(tmp_16*tmp_3 + tmp_17*tmp_19 - tmp_17*tmp_23 + tmp_22*tmp_5 - tmp_24*tmp_5 - tmp_25*tmp_3));
      real_t tmp_27 = tmp_26*(Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_28 = tmp_27*tmp_6;
      real_t tmp_29 = -tmp_1*tmp_17 + tmp_20*tmp_5;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = tmp_17*tmp_4 - tmp_20*tmp_3;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = tmp_19 - tmp_23;
      real_t tmp_34 = tmp_26*(-Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_35 = tmp_33*tmp_34;
      real_t tmp_36 = tmp_16 - tmp_25;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = tmp_22 - tmp_24;
      real_t tmp_39 = tmp_34*tmp_38;
      real_t tmp_40 = -tmp_18*tmp_3 + tmp_21*tmp_5;
      real_t tmp_41 = tmp_26*(Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0);
      real_t tmp_42 = tmp_40*tmp_41;
      real_t tmp_43 = -tmp_15*tmp_5 + tmp_17*tmp_18;
      real_t tmp_44 = tmp_41*tmp_43;
      real_t tmp_45 = tmp_15*tmp_3 - tmp_17*tmp_21;
      real_t tmp_46 = tmp_41*tmp_45;
      real_t tmp_47 = tmp_26*(-tmp_10 + tmp_8);
      real_t tmp_48 = tmp_26*(-tmp_12 + tmp_7);
      real_t tmp_49 = tmp_26*(-tmp_11 + tmp_9);
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
      real_t tmp_62 = 0.66666666666666652*Scalar_Variable_Coefficient_3D_mu_out0_id0*(-tmp_29*tmp_47 - tmp_31*tmp_47 - tmp_33*tmp_48 - tmp_36*tmp_48 - tmp_38*tmp_48 - tmp_40*tmp_49 - tmp_43*tmp_49 - tmp_45*tmp_49 - tmp_47*tmp_6)*std::abs(tmp_13)*std::abs(p_affine_0_0*tmp_52 - p_affine_0_0*tmp_59 + p_affine_0_1*tmp_55 - p_affine_0_1*tmp_60 + p_affine_0_2*tmp_58 - p_affine_0_2*tmp_61 - p_affine_1_0*tmp_52 + p_affine_1_0*tmp_59 - p_affine_1_1*tmp_55 + p_affine_1_1*tmp_60 - p_affine_1_2*tmp_58 + p_affine_1_2*tmp_61 + p_affine_2_0*tmp_54 - p_affine_2_0*tmp_57 - p_affine_2_1*tmp_51 + p_affine_2_1*tmp_56 + p_affine_2_2*tmp_50 - p_affine_2_2*tmp_53 - p_affine_3_0*tmp_54 + p_affine_3_0*tmp_57 + p_affine_3_1*tmp_51 - p_affine_3_1*tmp_56 - p_affine_3_2*tmp_50 + p_affine_3_2*tmp_53);
      real_t a_0_0 = tmp_62*(-tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46);
      real_t a_0_1 = tmp_62*(tmp_32 + tmp_39 + tmp_46);
      real_t a_0_2 = tmp_62*(tmp_30 + tmp_37 + tmp_44);
      real_t a_0_3 = tmp_62*(tmp_28 + tmp_35 + tmp_42);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_epsilonvar_0_1_blending_q1::Blending_DF_Triangle_blend( real_t in_0, real_t in_1, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3 ) const
   {
      Point3D  mappedPt( {in_0, in_1, 0} );
      Matrix2r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 1, 0 );
      *out_3 = DPsi( 1, 1 );
   }

   void p1_epsilonvar_0_1_blending_q1::Blending_F_Triangle_blend( real_t in_0, real_t in_1, real_t * out_0, real_t * out_1 ) const
   {
      Point3D  in( {in_0, in_1, 0} );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
   }

   void p1_epsilonvar_0_1_blending_q1::Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_mu( Point3D( {in_0, in_1, 0} ) );
   }

   void p1_epsilonvar_0_1_blending_q1::Blending_DF_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
   {
      Point3D  mappedPt( {in_0, in_1, in_2} );
      Matrix3r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 0, 2 );
      *out_3 = DPsi( 1, 0 );
      *out_4 = DPsi( 1, 1 );
      *out_5 = DPsi( 1, 2 );
      *out_6 = DPsi( 2, 0 );
      *out_7 = DPsi( 2, 1 );
      *out_8 = DPsi( 2, 2 );
   }

   void p1_epsilonvar_0_1_blending_q1::Blending_F_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2 ) const
   {
      Point3D  in( {in_0, in_1, in_2} );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
      *out_2 = out[2];
   }

   void p1_epsilonvar_0_1_blending_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

   void p1_epsilonvar_0_2_blending_q1::integrateAll( const std::array< Point3D, 3 >& , Matrix< real_t, 3, 3 >&  ) const
   {
      
   }

   void p1_epsilonvar_0_2_blending_q1::integrateRow0( const std::array< Point3D, 3 >& , Matrix< real_t, 1, 3 >&  ) const
   {
      
   }

   void p1_epsilonvar_0_2_blending_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_14 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_15 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_10 - Blending_DF_Tetrahedron_blend_out6_id0*tmp_14 + Blending_DF_Tetrahedron_blend_out7_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out7_id0*tmp_12 - Blending_DF_Tetrahedron_blend_out8_id0*tmp_13 + Blending_DF_Tetrahedron_blend_out8_id0*tmp_9;
      real_t tmp_16 = -p_affine_0_2;
      real_t tmp_17 = p_affine_3_2 + tmp_16;
      real_t tmp_18 = p_affine_1_2 + tmp_16;
      real_t tmp_19 = p_affine_3_1 + tmp_2;
      real_t tmp_20 = tmp_19*tmp_5;
      real_t tmp_21 = p_affine_2_2 + tmp_16;
      real_t tmp_22 = p_affine_3_0 + tmp_0;
      real_t tmp_23 = tmp_22*tmp_6;
      real_t tmp_24 = tmp_1*tmp_19;
      real_t tmp_25 = tmp_22*tmp_3;
      real_t tmp_26 = 0.5/(tmp_15*(tmp_17*tmp_4 - tmp_17*tmp_7 + tmp_18*tmp_20 - tmp_18*tmp_25 + tmp_21*tmp_23 - tmp_21*tmp_24));
      real_t tmp_27 = tmp_26*(Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_28 = tmp_27*tmp_8;
      real_t tmp_29 = tmp_23 - tmp_24;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = tmp_20 - tmp_25;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = -tmp_1*tmp_21 + tmp_18*tmp_5;
      real_t tmp_34 = tmp_26*(-Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_35 = tmp_33*tmp_34;
      real_t tmp_36 = tmp_1*tmp_17 - tmp_18*tmp_22;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = -tmp_17*tmp_5 + tmp_21*tmp_22;
      real_t tmp_39 = tmp_34*tmp_38;
      real_t tmp_40 = -tmp_18*tmp_3 + tmp_21*tmp_6;
      real_t tmp_41 = tmp_26*(Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0);
      real_t tmp_42 = tmp_40*tmp_41;
      real_t tmp_43 = -tmp_17*tmp_6 + tmp_18*tmp_19;
      real_t tmp_44 = tmp_41*tmp_43;
      real_t tmp_45 = tmp_17*tmp_3 - tmp_19*tmp_21;
      real_t tmp_46 = tmp_41*tmp_45;
      real_t tmp_47 = -tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46;
      real_t tmp_48 = tmp_26*(-tmp_13 + tmp_9);
      real_t tmp_49 = tmp_48*tmp_8;
      real_t tmp_50 = tmp_29*tmp_48;
      real_t tmp_51 = tmp_31*tmp_48;
      real_t tmp_52 = tmp_26*(tmp_11 - tmp_12);
      real_t tmp_53 = tmp_33*tmp_52;
      real_t tmp_54 = tmp_36*tmp_52;
      real_t tmp_55 = tmp_38*tmp_52;
      real_t tmp_56 = tmp_26*(tmp_10 - tmp_14);
      real_t tmp_57 = tmp_40*tmp_56;
      real_t tmp_58 = tmp_43*tmp_56;
      real_t tmp_59 = tmp_45*tmp_56;
      real_t tmp_60 = p_affine_0_0*p_affine_1_1;
      real_t tmp_61 = p_affine_0_0*p_affine_1_2;
      real_t tmp_62 = p_affine_2_1*p_affine_3_2;
      real_t tmp_63 = p_affine_0_1*p_affine_1_0;
      real_t tmp_64 = p_affine_0_1*p_affine_1_2;
      real_t tmp_65 = p_affine_2_2*p_affine_3_0;
      real_t tmp_66 = p_affine_0_2*p_affine_1_0;
      real_t tmp_67 = p_affine_0_2*p_affine_1_1;
      real_t tmp_68 = p_affine_2_0*p_affine_3_1;
      real_t tmp_69 = p_affine_2_2*p_affine_3_1;
      real_t tmp_70 = p_affine_2_0*p_affine_3_2;
      real_t tmp_71 = p_affine_2_1*p_affine_3_0;
      real_t tmp_72 = 0.66666666666666652*Scalar_Variable_Coefficient_3D_mu_out0_id0*std::abs(tmp_15)*std::abs(p_affine_0_0*tmp_62 - p_affine_0_0*tmp_69 + p_affine_0_1*tmp_65 - p_affine_0_1*tmp_70 + p_affine_0_2*tmp_68 - p_affine_0_2*tmp_71 - p_affine_1_0*tmp_62 + p_affine_1_0*tmp_69 - p_affine_1_1*tmp_65 + p_affine_1_1*tmp_70 - p_affine_1_2*tmp_68 + p_affine_1_2*tmp_71 + p_affine_2_0*tmp_64 - p_affine_2_0*tmp_67 - p_affine_2_1*tmp_61 + p_affine_2_1*tmp_66 + p_affine_2_2*tmp_60 - p_affine_2_2*tmp_63 - p_affine_3_0*tmp_64 + p_affine_3_0*tmp_67 + p_affine_3_1*tmp_61 - p_affine_3_1*tmp_66 - p_affine_3_2*tmp_60 + p_affine_3_2*tmp_63);
      real_t tmp_73 = tmp_72*(-tmp_49 - tmp_50 - tmp_51 - tmp_53 - tmp_54 - tmp_55 - tmp_57 - tmp_58 - tmp_59);
      real_t tmp_74 = tmp_32 + tmp_39 + tmp_46;
      real_t tmp_75 = tmp_30 + tmp_37 + tmp_44;
      real_t tmp_76 = tmp_28 + tmp_35 + tmp_42;
      real_t tmp_77 = tmp_72*(tmp_51 + tmp_55 + tmp_59);
      real_t tmp_78 = tmp_72*(tmp_50 + tmp_54 + tmp_58);
      real_t tmp_79 = tmp_72*(tmp_49 + tmp_53 + tmp_57);
      real_t a_0_0 = tmp_47*tmp_73;
      real_t a_0_1 = tmp_73*tmp_74;
      real_t a_0_2 = tmp_73*tmp_75;
      real_t a_0_3 = tmp_73*tmp_76;
      real_t a_1_0 = tmp_47*tmp_77;
      real_t a_1_1 = tmp_74*tmp_77;
      real_t a_1_2 = tmp_75*tmp_77;
      real_t a_1_3 = tmp_76*tmp_77;
      real_t a_2_0 = tmp_47*tmp_78;
      real_t a_2_1 = tmp_74*tmp_78;
      real_t a_2_2 = tmp_75*tmp_78;
      real_t a_2_3 = tmp_76*tmp_78;
      real_t a_3_0 = tmp_47*tmp_79;
      real_t a_3_1 = tmp_74*tmp_79;
      real_t a_3_2 = tmp_75*tmp_79;
      real_t a_3_3 = tmp_76*tmp_79;
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

   void p1_epsilonvar_0_2_blending_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_14 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_15 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_10 - Blending_DF_Tetrahedron_blend_out6_id0*tmp_14 + Blending_DF_Tetrahedron_blend_out7_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out7_id0*tmp_12 - Blending_DF_Tetrahedron_blend_out8_id0*tmp_13 + Blending_DF_Tetrahedron_blend_out8_id0*tmp_9;
      real_t tmp_16 = -p_affine_0_2;
      real_t tmp_17 = p_affine_3_2 + tmp_16;
      real_t tmp_18 = p_affine_1_2 + tmp_16;
      real_t tmp_19 = p_affine_3_1 + tmp_2;
      real_t tmp_20 = tmp_19*tmp_5;
      real_t tmp_21 = p_affine_2_2 + tmp_16;
      real_t tmp_22 = p_affine_3_0 + tmp_0;
      real_t tmp_23 = tmp_22*tmp_6;
      real_t tmp_24 = tmp_1*tmp_19;
      real_t tmp_25 = tmp_22*tmp_3;
      real_t tmp_26 = 0.5/(tmp_15*(tmp_17*tmp_4 - tmp_17*tmp_7 + tmp_18*tmp_20 - tmp_18*tmp_25 + tmp_21*tmp_23 - tmp_21*tmp_24));
      real_t tmp_27 = tmp_26*(Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_28 = tmp_27*tmp_8;
      real_t tmp_29 = tmp_23 - tmp_24;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = tmp_20 - tmp_25;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = -tmp_1*tmp_21 + tmp_18*tmp_5;
      real_t tmp_34 = tmp_26*(-Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_35 = tmp_33*tmp_34;
      real_t tmp_36 = tmp_1*tmp_17 - tmp_18*tmp_22;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = -tmp_17*tmp_5 + tmp_21*tmp_22;
      real_t tmp_39 = tmp_34*tmp_38;
      real_t tmp_40 = -tmp_18*tmp_3 + tmp_21*tmp_6;
      real_t tmp_41 = tmp_26*(Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0);
      real_t tmp_42 = tmp_40*tmp_41;
      real_t tmp_43 = -tmp_17*tmp_6 + tmp_18*tmp_19;
      real_t tmp_44 = tmp_41*tmp_43;
      real_t tmp_45 = tmp_17*tmp_3 - tmp_19*tmp_21;
      real_t tmp_46 = tmp_41*tmp_45;
      real_t tmp_47 = tmp_26*(-tmp_13 + tmp_9);
      real_t tmp_48 = tmp_26*(tmp_11 - tmp_12);
      real_t tmp_49 = tmp_26*(tmp_10 - tmp_14);
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
      real_t tmp_62 = 0.66666666666666652*Scalar_Variable_Coefficient_3D_mu_out0_id0*(-tmp_29*tmp_47 - tmp_31*tmp_47 - tmp_33*tmp_48 - tmp_36*tmp_48 - tmp_38*tmp_48 - tmp_40*tmp_49 - tmp_43*tmp_49 - tmp_45*tmp_49 - tmp_47*tmp_8)*std::abs(tmp_15)*std::abs(p_affine_0_0*tmp_52 - p_affine_0_0*tmp_59 + p_affine_0_1*tmp_55 - p_affine_0_1*tmp_60 + p_affine_0_2*tmp_58 - p_affine_0_2*tmp_61 - p_affine_1_0*tmp_52 + p_affine_1_0*tmp_59 - p_affine_1_1*tmp_55 + p_affine_1_1*tmp_60 - p_affine_1_2*tmp_58 + p_affine_1_2*tmp_61 + p_affine_2_0*tmp_54 - p_affine_2_0*tmp_57 - p_affine_2_1*tmp_51 + p_affine_2_1*tmp_56 + p_affine_2_2*tmp_50 - p_affine_2_2*tmp_53 - p_affine_3_0*tmp_54 + p_affine_3_0*tmp_57 + p_affine_3_1*tmp_51 - p_affine_3_1*tmp_56 - p_affine_3_2*tmp_50 + p_affine_3_2*tmp_53);
      real_t a_0_0 = tmp_62*(-tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46);
      real_t a_0_1 = tmp_62*(tmp_32 + tmp_39 + tmp_46);
      real_t a_0_2 = tmp_62*(tmp_30 + tmp_37 + tmp_44);
      real_t a_0_3 = tmp_62*(tmp_28 + tmp_35 + tmp_42);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_epsilonvar_0_2_blending_q1::Blending_DF_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
   {
      Point3D  mappedPt( {in_0, in_1, in_2} );
      Matrix3r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 0, 2 );
      *out_3 = DPsi( 1, 0 );
      *out_4 = DPsi( 1, 1 );
      *out_5 = DPsi( 1, 2 );
      *out_6 = DPsi( 2, 0 );
      *out_7 = DPsi( 2, 1 );
      *out_8 = DPsi( 2, 2 );
   }

   void p1_epsilonvar_0_2_blending_q1::Blending_F_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2 ) const
   {
      Point3D  in( {in_0, in_1, in_2} );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
      *out_2 = out[2];
   }

   void p1_epsilonvar_0_2_blending_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

   void p1_epsilonvar_1_0_blending_q1::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Blending_DF_Triangle_blend_out0_id0 = 0;
      real_t Blending_DF_Triangle_blend_out1_id0 = 0;
      real_t Blending_DF_Triangle_blend_out2_id0 = 0;
      real_t Blending_DF_Triangle_blend_out3_id0 = 0;
      real_t Blending_F_Triangle_blend_out0_id1 = 0;
      real_t Blending_F_Triangle_blend_out1_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      Blending_DF_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_DF_Triangle_blend_out0_id0, &Blending_DF_Triangle_blend_out1_id0, &Blending_DF_Triangle_blend_out2_id0, &Blending_DF_Triangle_blend_out3_id0 );
      Blending_F_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_F_Triangle_blend_out0_id1, &Blending_F_Triangle_blend_out1_id1 );
      Scalar_Variable_Coefficient_2D_mu( Blending_F_Triangle_blend_out0_id1, Blending_F_Triangle_blend_out1_id1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = Blending_DF_Triangle_blend_out0_id0*Blending_DF_Triangle_blend_out3_id0 - Blending_DF_Triangle_blend_out1_id0*Blending_DF_Triangle_blend_out2_id0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = 0.5/(tmp_2*(tmp_1*tmp_4 - (p_affine_1_1 + tmp_3)*(p_affine_2_0 + tmp_0)));
      real_t tmp_6 = Blending_DF_Triangle_blend_out0_id0*tmp_5;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_9 = tmp_6*tmp_8;
      real_t tmp_10 = Blending_DF_Triangle_blend_out1_id0*tmp_5;
      real_t tmp_11 = tmp_10*tmp_4;
      real_t tmp_12 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_13 = tmp_10*tmp_12;
      real_t tmp_14 = tmp_11 + tmp_13 - tmp_7 - tmp_9;
      real_t tmp_15 = Blending_DF_Triangle_blend_out2_id0*tmp_5;
      real_t tmp_16 = tmp_1*tmp_15;
      real_t tmp_17 = tmp_15*tmp_8;
      real_t tmp_18 = Blending_DF_Triangle_blend_out3_id0*tmp_5;
      real_t tmp_19 = tmp_18*tmp_4;
      real_t tmp_20 = tmp_12*tmp_18;
      real_t tmp_21 = 2.0*Scalar_Variable_Coefficient_2D_mu_out0_id0*std::abs(tmp_2)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_22 = tmp_21*(tmp_16 + tmp_17 - tmp_19 - tmp_20);
      real_t tmp_23 = -tmp_11 + tmp_9;
      real_t tmp_24 = -tmp_13 + tmp_7;
      real_t tmp_25 = tmp_21*(-tmp_17 + tmp_19);
      real_t tmp_26 = tmp_21*(-tmp_16 + tmp_20);
      real_t a_0_0 = tmp_14*tmp_22;
      real_t a_0_1 = tmp_22*tmp_23;
      real_t a_0_2 = tmp_22*tmp_24;
      real_t a_1_0 = tmp_14*tmp_25;
      real_t a_1_1 = tmp_23*tmp_25;
      real_t a_1_2 = tmp_24*tmp_25;
      real_t a_2_0 = tmp_14*tmp_26;
      real_t a_2_1 = tmp_23*tmp_26;
      real_t a_2_2 = tmp_24*tmp_26;
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

   void p1_epsilonvar_1_0_blending_q1::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Blending_DF_Triangle_blend_out0_id0 = 0;
      real_t Blending_DF_Triangle_blend_out1_id0 = 0;
      real_t Blending_DF_Triangle_blend_out2_id0 = 0;
      real_t Blending_DF_Triangle_blend_out3_id0 = 0;
      real_t Blending_F_Triangle_blend_out0_id1 = 0;
      real_t Blending_F_Triangle_blend_out1_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      Blending_DF_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_DF_Triangle_blend_out0_id0, &Blending_DF_Triangle_blend_out1_id0, &Blending_DF_Triangle_blend_out2_id0, &Blending_DF_Triangle_blend_out3_id0 );
      Blending_F_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_F_Triangle_blend_out0_id1, &Blending_F_Triangle_blend_out1_id1 );
      Scalar_Variable_Coefficient_2D_mu( Blending_F_Triangle_blend_out0_id1, Blending_F_Triangle_blend_out1_id1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = Blending_DF_Triangle_blend_out0_id0*Blending_DF_Triangle_blend_out3_id0 - Blending_DF_Triangle_blend_out1_id0*Blending_DF_Triangle_blend_out2_id0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = 0.5/(tmp_2*(tmp_1*tmp_4 - (p_affine_1_1 + tmp_3)*(p_affine_2_0 + tmp_0)));
      real_t tmp_6 = Blending_DF_Triangle_blend_out0_id0*tmp_5;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_9 = tmp_6*tmp_8;
      real_t tmp_10 = Blending_DF_Triangle_blend_out1_id0*tmp_5;
      real_t tmp_11 = tmp_10*tmp_4;
      real_t tmp_12 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_13 = tmp_10*tmp_12;
      real_t tmp_14 = Blending_DF_Triangle_blend_out2_id0*tmp_5;
      real_t tmp_15 = Blending_DF_Triangle_blend_out3_id0*tmp_5;
      real_t tmp_16 = 2.0*Scalar_Variable_Coefficient_2D_mu_out0_id0*(tmp_1*tmp_14 - tmp_12*tmp_15 + tmp_14*tmp_8 - tmp_15*tmp_4)*std::abs(tmp_2)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = tmp_16*(tmp_11 + tmp_13 - tmp_7 - tmp_9);
      real_t a_0_1 = tmp_16*(-tmp_11 + tmp_9);
      real_t a_0_2 = tmp_16*(-tmp_13 + tmp_7);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_epsilonvar_1_0_blending_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_8 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out3_id0*tmp_10 - Blending_DF_Tetrahedron_blend_out3_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out4_id0*tmp_12 + Blending_DF_Tetrahedron_blend_out4_id0*tmp_9 + Blending_DF_Tetrahedron_blend_out5_id0*tmp_7 - Blending_DF_Tetrahedron_blend_out5_id0*tmp_8;
      real_t tmp_14 = -p_affine_0_2;
      real_t tmp_15 = p_affine_3_2 + tmp_14;
      real_t tmp_16 = tmp_1*tmp_15;
      real_t tmp_17 = p_affine_3_1 + tmp_2;
      real_t tmp_18 = p_affine_1_2 + tmp_14;
      real_t tmp_19 = tmp_18*tmp_4;
      real_t tmp_20 = p_affine_3_0 + tmp_0;
      real_t tmp_21 = p_affine_2_2 + tmp_14;
      real_t tmp_22 = tmp_20*tmp_21;
      real_t tmp_23 = tmp_1*tmp_21;
      real_t tmp_24 = tmp_15*tmp_4;
      real_t tmp_25 = tmp_18*tmp_20;
      real_t tmp_26 = 0.5/(tmp_13*(tmp_16*tmp_3 + tmp_17*tmp_19 - tmp_17*tmp_23 + tmp_22*tmp_5 - tmp_24*tmp_5 - tmp_25*tmp_3));
      real_t tmp_27 = tmp_26*(tmp_7 - tmp_8);
      real_t tmp_28 = tmp_27*tmp_6;
      real_t tmp_29 = -tmp_1*tmp_17 + tmp_20*tmp_5;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = tmp_17*tmp_4 - tmp_20*tmp_3;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = tmp_19 - tmp_23;
      real_t tmp_34 = tmp_26*(-tmp_12 + tmp_9);
      real_t tmp_35 = tmp_33*tmp_34;
      real_t tmp_36 = tmp_16 - tmp_25;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = tmp_22 - tmp_24;
      real_t tmp_39 = tmp_34*tmp_38;
      real_t tmp_40 = -tmp_18*tmp_3 + tmp_21*tmp_5;
      real_t tmp_41 = tmp_26*(tmp_10 - tmp_11);
      real_t tmp_42 = tmp_40*tmp_41;
      real_t tmp_43 = -tmp_15*tmp_5 + tmp_17*tmp_18;
      real_t tmp_44 = tmp_41*tmp_43;
      real_t tmp_45 = tmp_15*tmp_3 - tmp_17*tmp_21;
      real_t tmp_46 = tmp_41*tmp_45;
      real_t tmp_47 = -tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46;
      real_t tmp_48 = tmp_26*(Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_49 = tmp_48*tmp_6;
      real_t tmp_50 = tmp_29*tmp_48;
      real_t tmp_51 = tmp_31*tmp_48;
      real_t tmp_52 = tmp_26*(-Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_53 = tmp_33*tmp_52;
      real_t tmp_54 = tmp_36*tmp_52;
      real_t tmp_55 = tmp_38*tmp_52;
      real_t tmp_56 = tmp_26*(Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0);
      real_t tmp_57 = tmp_40*tmp_56;
      real_t tmp_58 = tmp_43*tmp_56;
      real_t tmp_59 = tmp_45*tmp_56;
      real_t tmp_60 = p_affine_0_0*p_affine_1_1;
      real_t tmp_61 = p_affine_0_0*p_affine_1_2;
      real_t tmp_62 = p_affine_2_1*p_affine_3_2;
      real_t tmp_63 = p_affine_0_1*p_affine_1_0;
      real_t tmp_64 = p_affine_0_1*p_affine_1_2;
      real_t tmp_65 = p_affine_2_2*p_affine_3_0;
      real_t tmp_66 = p_affine_0_2*p_affine_1_0;
      real_t tmp_67 = p_affine_0_2*p_affine_1_1;
      real_t tmp_68 = p_affine_2_0*p_affine_3_1;
      real_t tmp_69 = p_affine_2_2*p_affine_3_1;
      real_t tmp_70 = p_affine_2_0*p_affine_3_2;
      real_t tmp_71 = p_affine_2_1*p_affine_3_0;
      real_t tmp_72 = 0.66666666666666652*Scalar_Variable_Coefficient_3D_mu_out0_id0*std::abs(tmp_13)*std::abs(p_affine_0_0*tmp_62 - p_affine_0_0*tmp_69 + p_affine_0_1*tmp_65 - p_affine_0_1*tmp_70 + p_affine_0_2*tmp_68 - p_affine_0_2*tmp_71 - p_affine_1_0*tmp_62 + p_affine_1_0*tmp_69 - p_affine_1_1*tmp_65 + p_affine_1_1*tmp_70 - p_affine_1_2*tmp_68 + p_affine_1_2*tmp_71 + p_affine_2_0*tmp_64 - p_affine_2_0*tmp_67 - p_affine_2_1*tmp_61 + p_affine_2_1*tmp_66 + p_affine_2_2*tmp_60 - p_affine_2_2*tmp_63 - p_affine_3_0*tmp_64 + p_affine_3_0*tmp_67 + p_affine_3_1*tmp_61 - p_affine_3_1*tmp_66 - p_affine_3_2*tmp_60 + p_affine_3_2*tmp_63);
      real_t tmp_73 = tmp_72*(-tmp_49 - tmp_50 - tmp_51 - tmp_53 - tmp_54 - tmp_55 - tmp_57 - tmp_58 - tmp_59);
      real_t tmp_74 = tmp_32 + tmp_39 + tmp_46;
      real_t tmp_75 = tmp_30 + tmp_37 + tmp_44;
      real_t tmp_76 = tmp_28 + tmp_35 + tmp_42;
      real_t tmp_77 = tmp_72*(tmp_51 + tmp_55 + tmp_59);
      real_t tmp_78 = tmp_72*(tmp_50 + tmp_54 + tmp_58);
      real_t tmp_79 = tmp_72*(tmp_49 + tmp_53 + tmp_57);
      real_t a_0_0 = tmp_47*tmp_73;
      real_t a_0_1 = tmp_73*tmp_74;
      real_t a_0_2 = tmp_73*tmp_75;
      real_t a_0_3 = tmp_73*tmp_76;
      real_t a_1_0 = tmp_47*tmp_77;
      real_t a_1_1 = tmp_74*tmp_77;
      real_t a_1_2 = tmp_75*tmp_77;
      real_t a_1_3 = tmp_76*tmp_77;
      real_t a_2_0 = tmp_47*tmp_78;
      real_t a_2_1 = tmp_74*tmp_78;
      real_t a_2_2 = tmp_75*tmp_78;
      real_t a_2_3 = tmp_76*tmp_78;
      real_t a_3_0 = tmp_47*tmp_79;
      real_t a_3_1 = tmp_74*tmp_79;
      real_t a_3_2 = tmp_75*tmp_79;
      real_t a_3_3 = tmp_76*tmp_79;
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

   void p1_epsilonvar_1_0_blending_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_8 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out3_id0*tmp_10 - Blending_DF_Tetrahedron_blend_out3_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out4_id0*tmp_12 + Blending_DF_Tetrahedron_blend_out4_id0*tmp_9 + Blending_DF_Tetrahedron_blend_out5_id0*tmp_7 - Blending_DF_Tetrahedron_blend_out5_id0*tmp_8;
      real_t tmp_14 = -p_affine_0_2;
      real_t tmp_15 = p_affine_3_2 + tmp_14;
      real_t tmp_16 = tmp_1*tmp_15;
      real_t tmp_17 = p_affine_3_1 + tmp_2;
      real_t tmp_18 = p_affine_1_2 + tmp_14;
      real_t tmp_19 = tmp_18*tmp_4;
      real_t tmp_20 = p_affine_3_0 + tmp_0;
      real_t tmp_21 = p_affine_2_2 + tmp_14;
      real_t tmp_22 = tmp_20*tmp_21;
      real_t tmp_23 = tmp_1*tmp_21;
      real_t tmp_24 = tmp_15*tmp_4;
      real_t tmp_25 = tmp_18*tmp_20;
      real_t tmp_26 = 0.5/(tmp_13*(tmp_16*tmp_3 + tmp_17*tmp_19 - tmp_17*tmp_23 + tmp_22*tmp_5 - tmp_24*tmp_5 - tmp_25*tmp_3));
      real_t tmp_27 = tmp_26*(tmp_7 - tmp_8);
      real_t tmp_28 = tmp_27*tmp_6;
      real_t tmp_29 = -tmp_1*tmp_17 + tmp_20*tmp_5;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = tmp_17*tmp_4 - tmp_20*tmp_3;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = tmp_19 - tmp_23;
      real_t tmp_34 = tmp_26*(-tmp_12 + tmp_9);
      real_t tmp_35 = tmp_33*tmp_34;
      real_t tmp_36 = tmp_16 - tmp_25;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = tmp_22 - tmp_24;
      real_t tmp_39 = tmp_34*tmp_38;
      real_t tmp_40 = -tmp_18*tmp_3 + tmp_21*tmp_5;
      real_t tmp_41 = tmp_26*(tmp_10 - tmp_11);
      real_t tmp_42 = tmp_40*tmp_41;
      real_t tmp_43 = -tmp_15*tmp_5 + tmp_17*tmp_18;
      real_t tmp_44 = tmp_41*tmp_43;
      real_t tmp_45 = tmp_15*tmp_3 - tmp_17*tmp_21;
      real_t tmp_46 = tmp_41*tmp_45;
      real_t tmp_47 = tmp_26*(Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_48 = tmp_26*(-Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_49 = tmp_26*(Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0);
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
      real_t tmp_62 = 0.66666666666666652*Scalar_Variable_Coefficient_3D_mu_out0_id0*(-tmp_29*tmp_47 - tmp_31*tmp_47 - tmp_33*tmp_48 - tmp_36*tmp_48 - tmp_38*tmp_48 - tmp_40*tmp_49 - tmp_43*tmp_49 - tmp_45*tmp_49 - tmp_47*tmp_6)*std::abs(tmp_13)*std::abs(p_affine_0_0*tmp_52 - p_affine_0_0*tmp_59 + p_affine_0_1*tmp_55 - p_affine_0_1*tmp_60 + p_affine_0_2*tmp_58 - p_affine_0_2*tmp_61 - p_affine_1_0*tmp_52 + p_affine_1_0*tmp_59 - p_affine_1_1*tmp_55 + p_affine_1_1*tmp_60 - p_affine_1_2*tmp_58 + p_affine_1_2*tmp_61 + p_affine_2_0*tmp_54 - p_affine_2_0*tmp_57 - p_affine_2_1*tmp_51 + p_affine_2_1*tmp_56 + p_affine_2_2*tmp_50 - p_affine_2_2*tmp_53 - p_affine_3_0*tmp_54 + p_affine_3_0*tmp_57 + p_affine_3_1*tmp_51 - p_affine_3_1*tmp_56 - p_affine_3_2*tmp_50 + p_affine_3_2*tmp_53);
      real_t a_0_0 = tmp_62*(-tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46);
      real_t a_0_1 = tmp_62*(tmp_32 + tmp_39 + tmp_46);
      real_t a_0_2 = tmp_62*(tmp_30 + tmp_37 + tmp_44);
      real_t a_0_3 = tmp_62*(tmp_28 + tmp_35 + tmp_42);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_epsilonvar_1_0_blending_q1::Blending_DF_Triangle_blend( real_t in_0, real_t in_1, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3 ) const
   {
      Point3D  mappedPt( {in_0, in_1, 0} );
      Matrix2r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 1, 0 );
      *out_3 = DPsi( 1, 1 );
   }

   void p1_epsilonvar_1_0_blending_q1::Blending_F_Triangle_blend( real_t in_0, real_t in_1, real_t * out_0, real_t * out_1 ) const
   {
      Point3D  in( {in_0, in_1, 0} );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
   }

   void p1_epsilonvar_1_0_blending_q1::Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_mu( Point3D( {in_0, in_1, 0} ) );
   }

   void p1_epsilonvar_1_0_blending_q1::Blending_DF_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
   {
      Point3D  mappedPt( {in_0, in_1, in_2} );
      Matrix3r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 0, 2 );
      *out_3 = DPsi( 1, 0 );
      *out_4 = DPsi( 1, 1 );
      *out_5 = DPsi( 1, 2 );
      *out_6 = DPsi( 2, 0 );
      *out_7 = DPsi( 2, 1 );
      *out_8 = DPsi( 2, 2 );
   }

   void p1_epsilonvar_1_0_blending_q1::Blending_F_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2 ) const
   {
      Point3D  in( {in_0, in_1, in_2} );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
      *out_2 = out[2];
   }

   void p1_epsilonvar_1_0_blending_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

   void p1_epsilonvar_1_1_blending_q1::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Blending_DF_Triangle_blend_out0_id0 = 0;
      real_t Blending_DF_Triangle_blend_out1_id0 = 0;
      real_t Blending_DF_Triangle_blend_out2_id0 = 0;
      real_t Blending_DF_Triangle_blend_out3_id0 = 0;
      real_t Blending_F_Triangle_blend_out0_id1 = 0;
      real_t Blending_F_Triangle_blend_out1_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      Blending_DF_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_DF_Triangle_blend_out0_id0, &Blending_DF_Triangle_blend_out1_id0, &Blending_DF_Triangle_blend_out2_id0, &Blending_DF_Triangle_blend_out3_id0 );
      Blending_F_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_F_Triangle_blend_out0_id1, &Blending_F_Triangle_blend_out1_id1 );
      Scalar_Variable_Coefficient_2D_mu( Blending_F_Triangle_blend_out0_id1, Blending_F_Triangle_blend_out1_id1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = Blending_DF_Triangle_blend_out0_id0*Blending_DF_Triangle_blend_out3_id0 - Blending_DF_Triangle_blend_out1_id0*Blending_DF_Triangle_blend_out2_id0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = 1/(tmp_2*(tmp_1*tmp_4 - (p_affine_1_1 + tmp_3)*(p_affine_2_0 + tmp_0)));
      real_t tmp_6 = Blending_DF_Triangle_blend_out2_id0*tmp_5;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_9 = tmp_6*tmp_8;
      real_t tmp_10 = Blending_DF_Triangle_blend_out3_id0*tmp_5;
      real_t tmp_11 = tmp_10*tmp_4;
      real_t tmp_12 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_13 = tmp_10*tmp_12;
      real_t tmp_14 = 1.0*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_15 = 1.0*tmp_5;
      real_t tmp_16 = Blending_DF_Triangle_blend_out0_id0*tmp_15;
      real_t tmp_17 = tmp_1*tmp_16;
      real_t tmp_18 = tmp_16*tmp_8;
      real_t tmp_19 = Blending_DF_Triangle_blend_out1_id0*tmp_15;
      real_t tmp_20 = tmp_19*tmp_4;
      real_t tmp_21 = tmp_12*tmp_19;
      real_t tmp_22 = -tmp_17 - tmp_18 + tmp_20 + tmp_21;
      real_t tmp_23 = 2*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_24 = 0.5*std::abs(tmp_2)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_25 = tmp_18 - tmp_20;
      real_t tmp_26 = tmp_22*tmp_23;
      real_t tmp_27 = 0.5*tmp_9;
      real_t tmp_28 = 0.5*tmp_11;
      real_t tmp_29 = -tmp_27 + tmp_28;
      real_t tmp_30 = 0.5*tmp_7;
      real_t tmp_31 = 0.5*tmp_13;
      real_t tmp_32 = 4*Scalar_Variable_Coefficient_2D_mu_out0_id0*(tmp_27 - tmp_28 + tmp_30 - tmp_31);
      real_t tmp_33 = tmp_24*(tmp_25*tmp_26 + tmp_29*tmp_32);
      real_t tmp_34 = tmp_17 - tmp_21;
      real_t tmp_35 = -tmp_30 + tmp_31;
      real_t tmp_36 = tmp_24*(tmp_26*tmp_34 + tmp_32*tmp_35);
      real_t tmp_37 = tmp_24*(4*Scalar_Variable_Coefficient_2D_mu_out0_id0*tmp_29*tmp_35 + tmp_23*tmp_25*tmp_34);
      real_t a_0_0 = tmp_24*(tmp_14*((-tmp_11 - tmp_13 + tmp_7 + tmp_9)*(-tmp_11 - tmp_13 + tmp_7 + tmp_9)) + (tmp_22*tmp_22)*tmp_23);
      real_t a_0_1 = tmp_33;
      real_t a_0_2 = tmp_36;
      real_t a_1_0 = tmp_33;
      real_t a_1_1 = tmp_24*(tmp_14*((tmp_11 - tmp_9)*(tmp_11 - tmp_9)) + tmp_23*(tmp_25*tmp_25));
      real_t a_1_2 = tmp_37;
      real_t a_2_0 = tmp_36;
      real_t a_2_1 = tmp_37;
      real_t a_2_2 = tmp_24*(tmp_14*((tmp_13 - tmp_7)*(tmp_13 - tmp_7)) + tmp_23*(tmp_34*tmp_34));
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

   void p1_epsilonvar_1_1_blending_q1::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Blending_DF_Triangle_blend_out0_id0 = 0;
      real_t Blending_DF_Triangle_blend_out1_id0 = 0;
      real_t Blending_DF_Triangle_blend_out2_id0 = 0;
      real_t Blending_DF_Triangle_blend_out3_id0 = 0;
      real_t Blending_F_Triangle_blend_out0_id1 = 0;
      real_t Blending_F_Triangle_blend_out1_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_mu_out0_id0 = 0;
      Blending_DF_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_DF_Triangle_blend_out0_id0, &Blending_DF_Triangle_blend_out1_id0, &Blending_DF_Triangle_blend_out2_id0, &Blending_DF_Triangle_blend_out3_id0 );
      Blending_F_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_F_Triangle_blend_out0_id1, &Blending_F_Triangle_blend_out1_id1 );
      Scalar_Variable_Coefficient_2D_mu( Blending_F_Triangle_blend_out0_id1, Blending_F_Triangle_blend_out1_id1, &Scalar_Variable_Coefficient_2D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = Blending_DF_Triangle_blend_out0_id0*Blending_DF_Triangle_blend_out3_id0 - Blending_DF_Triangle_blend_out1_id0*Blending_DF_Triangle_blend_out2_id0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = 1/(tmp_2*(tmp_1*tmp_4 - (p_affine_1_1 + tmp_3)*(p_affine_2_0 + tmp_0)));
      real_t tmp_6 = Blending_DF_Triangle_blend_out2_id0*tmp_5;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_9 = tmp_6*tmp_8;
      real_t tmp_10 = Blending_DF_Triangle_blend_out3_id0*tmp_5;
      real_t tmp_11 = tmp_10*tmp_4;
      real_t tmp_12 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_13 = tmp_10*tmp_12;
      real_t tmp_14 = 1.0*tmp_5;
      real_t tmp_15 = Blending_DF_Triangle_blend_out0_id0*tmp_14;
      real_t tmp_16 = tmp_1*tmp_15;
      real_t tmp_17 = tmp_15*tmp_8;
      real_t tmp_18 = Blending_DF_Triangle_blend_out1_id0*tmp_14;
      real_t tmp_19 = tmp_18*tmp_4;
      real_t tmp_20 = tmp_12*tmp_18;
      real_t tmp_21 = -tmp_16 - tmp_17 + tmp_19 + tmp_20;
      real_t tmp_22 = 2*Scalar_Variable_Coefficient_2D_mu_out0_id0;
      real_t tmp_23 = 0.5*std::abs(tmp_2)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_24 = tmp_21*tmp_22;
      real_t tmp_25 = 0.5*tmp_9;
      real_t tmp_26 = 0.5*tmp_11;
      real_t tmp_27 = 0.5*tmp_7;
      real_t tmp_28 = 0.5*tmp_13;
      real_t tmp_29 = 4*Scalar_Variable_Coefficient_2D_mu_out0_id0*(tmp_25 - tmp_26 + tmp_27 - tmp_28);
      real_t a_0_0 = tmp_23*(1.0*Scalar_Variable_Coefficient_2D_mu_out0_id0*((-tmp_11 - tmp_13 + tmp_7 + tmp_9)*(-tmp_11 - tmp_13 + tmp_7 + tmp_9)) + (tmp_21*tmp_21)*tmp_22);
      real_t a_0_1 = tmp_23*(tmp_24*(tmp_17 - tmp_19) + tmp_29*(-tmp_25 + tmp_26));
      real_t a_0_2 = tmp_23*(tmp_24*(tmp_16 - tmp_20) + tmp_29*(-tmp_27 + tmp_28));
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_epsilonvar_1_1_blending_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_8 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_13 = -Blending_DF_Tetrahedron_blend_out3_id0*tmp_11 + Blending_DF_Tetrahedron_blend_out3_id0*tmp_9 - Blending_DF_Tetrahedron_blend_out4_id0*tmp_12 + Blending_DF_Tetrahedron_blend_out4_id0*tmp_7 - Blending_DF_Tetrahedron_blend_out5_id0*tmp_10 + Blending_DF_Tetrahedron_blend_out5_id0*tmp_8;
      real_t tmp_14 = -p_affine_0_2;
      real_t tmp_15 = p_affine_3_2 + tmp_14;
      real_t tmp_16 = tmp_1*tmp_15;
      real_t tmp_17 = p_affine_3_1 + tmp_2;
      real_t tmp_18 = p_affine_1_2 + tmp_14;
      real_t tmp_19 = tmp_18*tmp_4;
      real_t tmp_20 = p_affine_3_0 + tmp_0;
      real_t tmp_21 = p_affine_2_2 + tmp_14;
      real_t tmp_22 = tmp_20*tmp_21;
      real_t tmp_23 = tmp_1*tmp_21;
      real_t tmp_24 = tmp_15*tmp_4;
      real_t tmp_25 = tmp_18*tmp_20;
      real_t tmp_26 = 1/(tmp_13*(tmp_16*tmp_3 + tmp_17*tmp_19 - tmp_17*tmp_23 + tmp_22*tmp_5 - tmp_24*tmp_5 - tmp_25*tmp_3));
      real_t tmp_27 = tmp_26*(Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0 - Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0);
      real_t tmp_28 = tmp_27*tmp_6;
      real_t tmp_29 = -tmp_1*tmp_17 + tmp_20*tmp_5;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = tmp_17*tmp_4 - tmp_20*tmp_3;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = tmp_19 - tmp_23;
      real_t tmp_34 = tmp_26*(-Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0 + Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0);
      real_t tmp_35 = tmp_33*tmp_34;
      real_t tmp_36 = tmp_16 - tmp_25;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = tmp_22 - tmp_24;
      real_t tmp_39 = tmp_34*tmp_38;
      real_t tmp_40 = -tmp_18*tmp_3 + tmp_21*tmp_5;
      real_t tmp_41 = tmp_26*(Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0 - Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0);
      real_t tmp_42 = tmp_40*tmp_41;
      real_t tmp_43 = -tmp_15*tmp_5 + tmp_17*tmp_18;
      real_t tmp_44 = tmp_41*tmp_43;
      real_t tmp_45 = tmp_15*tmp_3 - tmp_17*tmp_21;
      real_t tmp_46 = tmp_41*tmp_45;
      real_t tmp_47 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_48 = 1.0*tmp_26;
      real_t tmp_49 = tmp_48*(-tmp_10 + tmp_8);
      real_t tmp_50 = tmp_49*tmp_6;
      real_t tmp_51 = tmp_29*tmp_49;
      real_t tmp_52 = tmp_31*tmp_49;
      real_t tmp_53 = tmp_48*(-tmp_12 + tmp_7);
      real_t tmp_54 = tmp_33*tmp_53;
      real_t tmp_55 = tmp_36*tmp_53;
      real_t tmp_56 = tmp_38*tmp_53;
      real_t tmp_57 = tmp_48*(-tmp_11 + tmp_9);
      real_t tmp_58 = tmp_40*tmp_57;
      real_t tmp_59 = tmp_43*tmp_57;
      real_t tmp_60 = tmp_45*tmp_57;
      real_t tmp_61 = -tmp_50 - tmp_51 - tmp_52 - tmp_54 - tmp_55 - tmp_56 - tmp_58 - tmp_59 - tmp_60;
      real_t tmp_62 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_63 = tmp_26*(Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_64 = tmp_6*tmp_63;
      real_t tmp_65 = tmp_29*tmp_63;
      real_t tmp_66 = tmp_31*tmp_63;
      real_t tmp_67 = tmp_26*(-Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_68 = tmp_33*tmp_67;
      real_t tmp_69 = tmp_36*tmp_67;
      real_t tmp_70 = tmp_38*tmp_67;
      real_t tmp_71 = tmp_26*(Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0);
      real_t tmp_72 = tmp_40*tmp_71;
      real_t tmp_73 = tmp_43*tmp_71;
      real_t tmp_74 = tmp_45*tmp_71;
      real_t tmp_75 = p_affine_0_0*p_affine_1_1;
      real_t tmp_76 = p_affine_0_0*p_affine_1_2;
      real_t tmp_77 = p_affine_2_1*p_affine_3_2;
      real_t tmp_78 = p_affine_0_1*p_affine_1_0;
      real_t tmp_79 = p_affine_0_1*p_affine_1_2;
      real_t tmp_80 = p_affine_2_2*p_affine_3_0;
      real_t tmp_81 = p_affine_0_2*p_affine_1_0;
      real_t tmp_82 = p_affine_0_2*p_affine_1_1;
      real_t tmp_83 = p_affine_2_0*p_affine_3_1;
      real_t tmp_84 = p_affine_2_2*p_affine_3_1;
      real_t tmp_85 = p_affine_2_0*p_affine_3_2;
      real_t tmp_86 = p_affine_2_1*p_affine_3_0;
      real_t tmp_87 = 0.16666666666666663*std::abs(tmp_13)*std::abs(p_affine_0_0*tmp_77 - p_affine_0_0*tmp_84 + p_affine_0_1*tmp_80 - p_affine_0_1*tmp_85 + p_affine_0_2*tmp_83 - p_affine_0_2*tmp_86 - p_affine_1_0*tmp_77 + p_affine_1_0*tmp_84 - p_affine_1_1*tmp_80 + p_affine_1_1*tmp_85 - p_affine_1_2*tmp_83 + p_affine_1_2*tmp_86 + p_affine_2_0*tmp_79 - p_affine_2_0*tmp_82 - p_affine_2_1*tmp_76 + p_affine_2_1*tmp_81 + p_affine_2_2*tmp_75 - p_affine_2_2*tmp_78 - p_affine_3_0*tmp_79 + p_affine_3_0*tmp_82 + p_affine_3_1*tmp_76 - p_affine_3_1*tmp_81 - p_affine_3_2*tmp_75 + p_affine_3_2*tmp_78);
      real_t tmp_88 = 0.5*tmp_32;
      real_t tmp_89 = 0.5*tmp_39;
      real_t tmp_90 = 0.5*tmp_46;
      real_t tmp_91 = tmp_88 + tmp_89 + tmp_90;
      real_t tmp_92 = 0.5*tmp_28;
      real_t tmp_93 = 0.5*tmp_30;
      real_t tmp_94 = 0.5*tmp_35;
      real_t tmp_95 = 0.5*tmp_37;
      real_t tmp_96 = 0.5*tmp_42;
      real_t tmp_97 = 0.5*tmp_44;
      real_t tmp_98 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_99 = tmp_98*(-tmp_88 - tmp_89 - tmp_90 - tmp_92 - tmp_93 - tmp_94 - tmp_95 - tmp_96 - tmp_97);
      real_t tmp_100 = tmp_52 + tmp_56 + tmp_60;
      real_t tmp_101 = tmp_61*tmp_62;
      real_t tmp_102 = 0.5*tmp_66;
      real_t tmp_103 = 0.5*tmp_70;
      real_t tmp_104 = 0.5*tmp_74;
      real_t tmp_105 = tmp_102 + tmp_103 + tmp_104;
      real_t tmp_106 = 0.5*tmp_64;
      real_t tmp_107 = 0.5*tmp_65;
      real_t tmp_108 = 0.5*tmp_68;
      real_t tmp_109 = 0.5*tmp_69;
      real_t tmp_110 = 0.5*tmp_72;
      real_t tmp_111 = 0.5*tmp_73;
      real_t tmp_112 = tmp_98*(-tmp_102 - tmp_103 - tmp_104 - tmp_106 - tmp_107 - tmp_108 - tmp_109 - tmp_110 - tmp_111);
      real_t tmp_113 = tmp_87*(tmp_100*tmp_101 + tmp_105*tmp_112 + tmp_91*tmp_99);
      real_t tmp_114 = tmp_93 + tmp_95 + tmp_97;
      real_t tmp_115 = tmp_51 + tmp_55 + tmp_59;
      real_t tmp_116 = tmp_107 + tmp_109 + tmp_111;
      real_t tmp_117 = tmp_87*(tmp_101*tmp_115 + tmp_112*tmp_116 + tmp_114*tmp_99);
      real_t tmp_118 = tmp_92 + tmp_94 + tmp_96;
      real_t tmp_119 = tmp_50 + tmp_54 + tmp_58;
      real_t tmp_120 = tmp_106 + tmp_108 + tmp_110;
      real_t tmp_121 = tmp_87*(tmp_101*tmp_119 + tmp_112*tmp_120 + tmp_118*tmp_99);
      real_t tmp_122 = tmp_91*tmp_98;
      real_t tmp_123 = tmp_100*tmp_62;
      real_t tmp_124 = tmp_105*tmp_98;
      real_t tmp_125 = tmp_87*(tmp_114*tmp_122 + tmp_115*tmp_123 + tmp_116*tmp_124);
      real_t tmp_126 = tmp_87*(tmp_118*tmp_122 + tmp_119*tmp_123 + tmp_120*tmp_124);
      real_t tmp_127 = tmp_87*(tmp_114*tmp_118*tmp_98 + tmp_115*tmp_119*tmp_62 + tmp_116*tmp_120*tmp_98);
      real_t a_0_0 = tmp_87*(tmp_47*((-tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46)*(-tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46)) + tmp_47*((-tmp_64 - tmp_65 - tmp_66 - tmp_68 - tmp_69 - tmp_70 - tmp_72 - tmp_73 - tmp_74)*(-tmp_64 - tmp_65 - tmp_66 - tmp_68 - tmp_69 - tmp_70 - tmp_72 - tmp_73 - tmp_74)) + (tmp_61*tmp_61)*tmp_62);
      real_t a_0_1 = tmp_113;
      real_t a_0_2 = tmp_117;
      real_t a_0_3 = tmp_121;
      real_t a_1_0 = tmp_113;
      real_t a_1_1 = tmp_87*((tmp_100*tmp_100)*tmp_62 + tmp_47*((tmp_32 + tmp_39 + tmp_46)*(tmp_32 + tmp_39 + tmp_46)) + tmp_47*((tmp_66 + tmp_70 + tmp_74)*(tmp_66 + tmp_70 + tmp_74)));
      real_t a_1_2 = tmp_125;
      real_t a_1_3 = tmp_126;
      real_t a_2_0 = tmp_117;
      real_t a_2_1 = tmp_125;
      real_t a_2_2 = tmp_87*((tmp_115*tmp_115)*tmp_62 + tmp_47*((tmp_30 + tmp_37 + tmp_44)*(tmp_30 + tmp_37 + tmp_44)) + tmp_47*((tmp_65 + tmp_69 + tmp_73)*(tmp_65 + tmp_69 + tmp_73)));
      real_t a_2_3 = tmp_127;
      real_t a_3_0 = tmp_121;
      real_t a_3_1 = tmp_126;
      real_t a_3_2 = tmp_127;
      real_t a_3_3 = tmp_87*((tmp_119*tmp_119)*tmp_62 + tmp_47*((tmp_28 + tmp_35 + tmp_42)*(tmp_28 + tmp_35 + tmp_42)) + tmp_47*((tmp_64 + tmp_68 + tmp_72)*(tmp_64 + tmp_68 + tmp_72)));
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

   void p1_epsilonvar_1_1_blending_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = tmp_1*tmp_3 - tmp_4*tmp_5;
      real_t tmp_7 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_8 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_13 = -Blending_DF_Tetrahedron_blend_out3_id0*tmp_11 + Blending_DF_Tetrahedron_blend_out3_id0*tmp_9 - Blending_DF_Tetrahedron_blend_out4_id0*tmp_12 + Blending_DF_Tetrahedron_blend_out4_id0*tmp_7 - Blending_DF_Tetrahedron_blend_out5_id0*tmp_10 + Blending_DF_Tetrahedron_blend_out5_id0*tmp_8;
      real_t tmp_14 = -p_affine_0_2;
      real_t tmp_15 = p_affine_3_2 + tmp_14;
      real_t tmp_16 = tmp_1*tmp_15;
      real_t tmp_17 = p_affine_3_1 + tmp_2;
      real_t tmp_18 = p_affine_1_2 + tmp_14;
      real_t tmp_19 = tmp_18*tmp_4;
      real_t tmp_20 = p_affine_3_0 + tmp_0;
      real_t tmp_21 = p_affine_2_2 + tmp_14;
      real_t tmp_22 = tmp_20*tmp_21;
      real_t tmp_23 = tmp_1*tmp_21;
      real_t tmp_24 = tmp_15*tmp_4;
      real_t tmp_25 = tmp_18*tmp_20;
      real_t tmp_26 = 1/(tmp_13*(tmp_16*tmp_3 + tmp_17*tmp_19 - tmp_17*tmp_23 + tmp_22*tmp_5 - tmp_24*tmp_5 - tmp_25*tmp_3));
      real_t tmp_27 = tmp_26*(Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0 - Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0);
      real_t tmp_28 = tmp_27*tmp_6;
      real_t tmp_29 = -tmp_1*tmp_17 + tmp_20*tmp_5;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = tmp_17*tmp_4 - tmp_20*tmp_3;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = tmp_19 - tmp_23;
      real_t tmp_34 = tmp_26*(-Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0 + Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0);
      real_t tmp_35 = tmp_33*tmp_34;
      real_t tmp_36 = tmp_16 - tmp_25;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = tmp_22 - tmp_24;
      real_t tmp_39 = tmp_34*tmp_38;
      real_t tmp_40 = -tmp_18*tmp_3 + tmp_21*tmp_5;
      real_t tmp_41 = tmp_26*(Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0 - Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0);
      real_t tmp_42 = tmp_40*tmp_41;
      real_t tmp_43 = -tmp_15*tmp_5 + tmp_17*tmp_18;
      real_t tmp_44 = tmp_41*tmp_43;
      real_t tmp_45 = tmp_15*tmp_3 - tmp_17*tmp_21;
      real_t tmp_46 = tmp_41*tmp_45;
      real_t tmp_47 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_48 = 1.0*tmp_26;
      real_t tmp_49 = tmp_48*(-tmp_10 + tmp_8);
      real_t tmp_50 = tmp_49*tmp_6;
      real_t tmp_51 = tmp_29*tmp_49;
      real_t tmp_52 = tmp_31*tmp_49;
      real_t tmp_53 = tmp_48*(-tmp_12 + tmp_7);
      real_t tmp_54 = tmp_33*tmp_53;
      real_t tmp_55 = tmp_36*tmp_53;
      real_t tmp_56 = tmp_38*tmp_53;
      real_t tmp_57 = tmp_48*(-tmp_11 + tmp_9);
      real_t tmp_58 = tmp_40*tmp_57;
      real_t tmp_59 = tmp_43*tmp_57;
      real_t tmp_60 = tmp_45*tmp_57;
      real_t tmp_61 = -tmp_50 - tmp_51 - tmp_52 - tmp_54 - tmp_55 - tmp_56 - tmp_58 - tmp_59 - tmp_60;
      real_t tmp_62 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_63 = tmp_26*(Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_64 = tmp_6*tmp_63;
      real_t tmp_65 = tmp_29*tmp_63;
      real_t tmp_66 = tmp_31*tmp_63;
      real_t tmp_67 = tmp_26*(-Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_68 = tmp_33*tmp_67;
      real_t tmp_69 = tmp_36*tmp_67;
      real_t tmp_70 = tmp_38*tmp_67;
      real_t tmp_71 = tmp_26*(Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0);
      real_t tmp_72 = tmp_40*tmp_71;
      real_t tmp_73 = tmp_43*tmp_71;
      real_t tmp_74 = tmp_45*tmp_71;
      real_t tmp_75 = p_affine_0_0*p_affine_1_1;
      real_t tmp_76 = p_affine_0_0*p_affine_1_2;
      real_t tmp_77 = p_affine_2_1*p_affine_3_2;
      real_t tmp_78 = p_affine_0_1*p_affine_1_0;
      real_t tmp_79 = p_affine_0_1*p_affine_1_2;
      real_t tmp_80 = p_affine_2_2*p_affine_3_0;
      real_t tmp_81 = p_affine_0_2*p_affine_1_0;
      real_t tmp_82 = p_affine_0_2*p_affine_1_1;
      real_t tmp_83 = p_affine_2_0*p_affine_3_1;
      real_t tmp_84 = p_affine_2_2*p_affine_3_1;
      real_t tmp_85 = p_affine_2_0*p_affine_3_2;
      real_t tmp_86 = p_affine_2_1*p_affine_3_0;
      real_t tmp_87 = 0.16666666666666663*std::abs(tmp_13)*std::abs(p_affine_0_0*tmp_77 - p_affine_0_0*tmp_84 + p_affine_0_1*tmp_80 - p_affine_0_1*tmp_85 + p_affine_0_2*tmp_83 - p_affine_0_2*tmp_86 - p_affine_1_0*tmp_77 + p_affine_1_0*tmp_84 - p_affine_1_1*tmp_80 + p_affine_1_1*tmp_85 - p_affine_1_2*tmp_83 + p_affine_1_2*tmp_86 + p_affine_2_0*tmp_79 - p_affine_2_0*tmp_82 - p_affine_2_1*tmp_76 + p_affine_2_1*tmp_81 + p_affine_2_2*tmp_75 - p_affine_2_2*tmp_78 - p_affine_3_0*tmp_79 + p_affine_3_0*tmp_82 + p_affine_3_1*tmp_76 - p_affine_3_1*tmp_81 - p_affine_3_2*tmp_75 + p_affine_3_2*tmp_78);
      real_t tmp_88 = 0.5*tmp_32;
      real_t tmp_89 = 0.5*tmp_39;
      real_t tmp_90 = 0.5*tmp_46;
      real_t tmp_91 = 0.5*tmp_28;
      real_t tmp_92 = 0.5*tmp_30;
      real_t tmp_93 = 0.5*tmp_35;
      real_t tmp_94 = 0.5*tmp_37;
      real_t tmp_95 = 0.5*tmp_42;
      real_t tmp_96 = 0.5*tmp_44;
      real_t tmp_97 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_98 = tmp_97*(-tmp_88 - tmp_89 - tmp_90 - tmp_91 - tmp_92 - tmp_93 - tmp_94 - tmp_95 - tmp_96);
      real_t tmp_99 = tmp_61*tmp_62;
      real_t tmp_100 = 0.5*tmp_66;
      real_t tmp_101 = 0.5*tmp_70;
      real_t tmp_102 = 0.5*tmp_74;
      real_t tmp_103 = 0.5*tmp_64;
      real_t tmp_104 = 0.5*tmp_65;
      real_t tmp_105 = 0.5*tmp_68;
      real_t tmp_106 = 0.5*tmp_69;
      real_t tmp_107 = 0.5*tmp_72;
      real_t tmp_108 = 0.5*tmp_73;
      real_t tmp_109 = tmp_97*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_104 - tmp_105 - tmp_106 - tmp_107 - tmp_108);
      real_t a_0_0 = tmp_87*(tmp_47*((-tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46)*(-tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46)) + tmp_47*((-tmp_64 - tmp_65 - tmp_66 - tmp_68 - tmp_69 - tmp_70 - tmp_72 - tmp_73 - tmp_74)*(-tmp_64 - tmp_65 - tmp_66 - tmp_68 - tmp_69 - tmp_70 - tmp_72 - tmp_73 - tmp_74)) + (tmp_61*tmp_61)*tmp_62);
      real_t a_0_1 = tmp_87*(tmp_109*(tmp_100 + tmp_101 + tmp_102) + tmp_98*(tmp_88 + tmp_89 + tmp_90) + tmp_99*(tmp_52 + tmp_56 + tmp_60));
      real_t a_0_2 = tmp_87*(tmp_109*(tmp_104 + tmp_106 + tmp_108) + tmp_98*(tmp_92 + tmp_94 + tmp_96) + tmp_99*(tmp_51 + tmp_55 + tmp_59));
      real_t a_0_3 = tmp_87*(tmp_109*(tmp_103 + tmp_105 + tmp_107) + tmp_98*(tmp_91 + tmp_93 + tmp_95) + tmp_99*(tmp_50 + tmp_54 + tmp_58));
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_epsilonvar_1_1_blending_q1::Blending_DF_Triangle_blend( real_t in_0, real_t in_1, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3 ) const
   {
      Point3D  mappedPt( {in_0, in_1, 0} );
      Matrix2r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 1, 0 );
      *out_3 = DPsi( 1, 1 );
   }

   void p1_epsilonvar_1_1_blending_q1::Blending_F_Triangle_blend( real_t in_0, real_t in_1, real_t * out_0, real_t * out_1 ) const
   {
      Point3D  in( {in_0, in_1, 0} );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
   }

   void p1_epsilonvar_1_1_blending_q1::Scalar_Variable_Coefficient_2D_mu( real_t in_0, real_t in_1, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_mu( Point3D( {in_0, in_1, 0} ) );
   }

   void p1_epsilonvar_1_1_blending_q1::Blending_DF_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
   {
      Point3D  mappedPt( {in_0, in_1, in_2} );
      Matrix3r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 0, 2 );
      *out_3 = DPsi( 1, 0 );
      *out_4 = DPsi( 1, 1 );
      *out_5 = DPsi( 1, 2 );
      *out_6 = DPsi( 2, 0 );
      *out_7 = DPsi( 2, 1 );
      *out_8 = DPsi( 2, 2 );
   }

   void p1_epsilonvar_1_1_blending_q1::Blending_F_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2 ) const
   {
      Point3D  in( {in_0, in_1, in_2} );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
      *out_2 = out[2];
   }

   void p1_epsilonvar_1_1_blending_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

   void p1_epsilonvar_1_2_blending_q1::integrateAll( const std::array< Point3D, 3 >& , Matrix< real_t, 3, 3 >&  ) const
   {
      
   }

   void p1_epsilonvar_1_2_blending_q1::integrateRow0( const std::array< Point3D, 3 >& , Matrix< real_t, 1, 3 >&  ) const
   {
      
   }

   void p1_epsilonvar_1_2_blending_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_14 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_15 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_10 - Blending_DF_Tetrahedron_blend_out6_id0*tmp_14 + Blending_DF_Tetrahedron_blend_out7_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out7_id0*tmp_12 - Blending_DF_Tetrahedron_blend_out8_id0*tmp_13 + Blending_DF_Tetrahedron_blend_out8_id0*tmp_9;
      real_t tmp_16 = -p_affine_0_2;
      real_t tmp_17 = p_affine_3_2 + tmp_16;
      real_t tmp_18 = p_affine_1_2 + tmp_16;
      real_t tmp_19 = p_affine_3_1 + tmp_2;
      real_t tmp_20 = tmp_19*tmp_5;
      real_t tmp_21 = p_affine_2_2 + tmp_16;
      real_t tmp_22 = p_affine_3_0 + tmp_0;
      real_t tmp_23 = tmp_22*tmp_6;
      real_t tmp_24 = tmp_1*tmp_19;
      real_t tmp_25 = tmp_22*tmp_3;
      real_t tmp_26 = 0.5/(tmp_15*(tmp_17*tmp_4 - tmp_17*tmp_7 + tmp_18*tmp_20 - tmp_18*tmp_25 + tmp_21*tmp_23 - tmp_21*tmp_24));
      real_t tmp_27 = tmp_26*(-Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out7_id0 + Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_28 = tmp_27*tmp_8;
      real_t tmp_29 = tmp_23 - tmp_24;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = tmp_20 - tmp_25;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = -tmp_1*tmp_21 + tmp_18*tmp_5;
      real_t tmp_34 = tmp_26*(Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_35 = tmp_33*tmp_34;
      real_t tmp_36 = tmp_1*tmp_17 - tmp_18*tmp_22;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = -tmp_17*tmp_5 + tmp_21*tmp_22;
      real_t tmp_39 = tmp_34*tmp_38;
      real_t tmp_40 = -tmp_18*tmp_3 + tmp_21*tmp_6;
      real_t tmp_41 = tmp_26*(-Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out7_id0);
      real_t tmp_42 = tmp_40*tmp_41;
      real_t tmp_43 = -tmp_17*tmp_6 + tmp_18*tmp_19;
      real_t tmp_44 = tmp_41*tmp_43;
      real_t tmp_45 = tmp_17*tmp_3 - tmp_19*tmp_21;
      real_t tmp_46 = tmp_41*tmp_45;
      real_t tmp_47 = -tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46;
      real_t tmp_48 = tmp_26*(-tmp_13 + tmp_9);
      real_t tmp_49 = tmp_48*tmp_8;
      real_t tmp_50 = tmp_29*tmp_48;
      real_t tmp_51 = tmp_31*tmp_48;
      real_t tmp_52 = tmp_26*(tmp_11 - tmp_12);
      real_t tmp_53 = tmp_33*tmp_52;
      real_t tmp_54 = tmp_36*tmp_52;
      real_t tmp_55 = tmp_38*tmp_52;
      real_t tmp_56 = tmp_26*(tmp_10 - tmp_14);
      real_t tmp_57 = tmp_40*tmp_56;
      real_t tmp_58 = tmp_43*tmp_56;
      real_t tmp_59 = tmp_45*tmp_56;
      real_t tmp_60 = p_affine_0_0*p_affine_1_1;
      real_t tmp_61 = p_affine_0_0*p_affine_1_2;
      real_t tmp_62 = p_affine_2_1*p_affine_3_2;
      real_t tmp_63 = p_affine_0_1*p_affine_1_0;
      real_t tmp_64 = p_affine_0_1*p_affine_1_2;
      real_t tmp_65 = p_affine_2_2*p_affine_3_0;
      real_t tmp_66 = p_affine_0_2*p_affine_1_0;
      real_t tmp_67 = p_affine_0_2*p_affine_1_1;
      real_t tmp_68 = p_affine_2_0*p_affine_3_1;
      real_t tmp_69 = p_affine_2_2*p_affine_3_1;
      real_t tmp_70 = p_affine_2_0*p_affine_3_2;
      real_t tmp_71 = p_affine_2_1*p_affine_3_0;
      real_t tmp_72 = 0.66666666666666652*Scalar_Variable_Coefficient_3D_mu_out0_id0*std::abs(tmp_15)*std::abs(p_affine_0_0*tmp_62 - p_affine_0_0*tmp_69 + p_affine_0_1*tmp_65 - p_affine_0_1*tmp_70 + p_affine_0_2*tmp_68 - p_affine_0_2*tmp_71 - p_affine_1_0*tmp_62 + p_affine_1_0*tmp_69 - p_affine_1_1*tmp_65 + p_affine_1_1*tmp_70 - p_affine_1_2*tmp_68 + p_affine_1_2*tmp_71 + p_affine_2_0*tmp_64 - p_affine_2_0*tmp_67 - p_affine_2_1*tmp_61 + p_affine_2_1*tmp_66 + p_affine_2_2*tmp_60 - p_affine_2_2*tmp_63 - p_affine_3_0*tmp_64 + p_affine_3_0*tmp_67 + p_affine_3_1*tmp_61 - p_affine_3_1*tmp_66 - p_affine_3_2*tmp_60 + p_affine_3_2*tmp_63);
      real_t tmp_73 = tmp_72*(-tmp_49 - tmp_50 - tmp_51 - tmp_53 - tmp_54 - tmp_55 - tmp_57 - tmp_58 - tmp_59);
      real_t tmp_74 = tmp_32 + tmp_39 + tmp_46;
      real_t tmp_75 = tmp_30 + tmp_37 + tmp_44;
      real_t tmp_76 = tmp_28 + tmp_35 + tmp_42;
      real_t tmp_77 = tmp_72*(tmp_51 + tmp_55 + tmp_59);
      real_t tmp_78 = tmp_72*(tmp_50 + tmp_54 + tmp_58);
      real_t tmp_79 = tmp_72*(tmp_49 + tmp_53 + tmp_57);
      real_t a_0_0 = tmp_47*tmp_73;
      real_t a_0_1 = tmp_73*tmp_74;
      real_t a_0_2 = tmp_73*tmp_75;
      real_t a_0_3 = tmp_73*tmp_76;
      real_t a_1_0 = tmp_47*tmp_77;
      real_t a_1_1 = tmp_74*tmp_77;
      real_t a_1_2 = tmp_75*tmp_77;
      real_t a_1_3 = tmp_76*tmp_77;
      real_t a_2_0 = tmp_47*tmp_78;
      real_t a_2_1 = tmp_74*tmp_78;
      real_t a_2_2 = tmp_75*tmp_78;
      real_t a_2_3 = tmp_76*tmp_78;
      real_t a_3_0 = tmp_47*tmp_79;
      real_t a_3_1 = tmp_74*tmp_79;
      real_t a_3_2 = tmp_75*tmp_79;
      real_t a_3_3 = tmp_76*tmp_79;
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

   void p1_epsilonvar_1_2_blending_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_14 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_15 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_10 - Blending_DF_Tetrahedron_blend_out6_id0*tmp_14 + Blending_DF_Tetrahedron_blend_out7_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out7_id0*tmp_12 - Blending_DF_Tetrahedron_blend_out8_id0*tmp_13 + Blending_DF_Tetrahedron_blend_out8_id0*tmp_9;
      real_t tmp_16 = -p_affine_0_2;
      real_t tmp_17 = p_affine_3_2 + tmp_16;
      real_t tmp_18 = p_affine_1_2 + tmp_16;
      real_t tmp_19 = p_affine_3_1 + tmp_2;
      real_t tmp_20 = tmp_19*tmp_5;
      real_t tmp_21 = p_affine_2_2 + tmp_16;
      real_t tmp_22 = p_affine_3_0 + tmp_0;
      real_t tmp_23 = tmp_22*tmp_6;
      real_t tmp_24 = tmp_1*tmp_19;
      real_t tmp_25 = tmp_22*tmp_3;
      real_t tmp_26 = 0.5/(tmp_15*(tmp_17*tmp_4 - tmp_17*tmp_7 + tmp_18*tmp_20 - tmp_18*tmp_25 + tmp_21*tmp_23 - tmp_21*tmp_24));
      real_t tmp_27 = tmp_26*(-Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out7_id0 + Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_28 = tmp_27*tmp_8;
      real_t tmp_29 = tmp_23 - tmp_24;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = tmp_20 - tmp_25;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = -tmp_1*tmp_21 + tmp_18*tmp_5;
      real_t tmp_34 = tmp_26*(Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_35 = tmp_33*tmp_34;
      real_t tmp_36 = tmp_1*tmp_17 - tmp_18*tmp_22;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = -tmp_17*tmp_5 + tmp_21*tmp_22;
      real_t tmp_39 = tmp_34*tmp_38;
      real_t tmp_40 = -tmp_18*tmp_3 + tmp_21*tmp_6;
      real_t tmp_41 = tmp_26*(-Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out7_id0);
      real_t tmp_42 = tmp_40*tmp_41;
      real_t tmp_43 = -tmp_17*tmp_6 + tmp_18*tmp_19;
      real_t tmp_44 = tmp_41*tmp_43;
      real_t tmp_45 = tmp_17*tmp_3 - tmp_19*tmp_21;
      real_t tmp_46 = tmp_41*tmp_45;
      real_t tmp_47 = tmp_26*(-tmp_13 + tmp_9);
      real_t tmp_48 = tmp_26*(tmp_11 - tmp_12);
      real_t tmp_49 = tmp_26*(tmp_10 - tmp_14);
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
      real_t tmp_62 = 0.66666666666666652*Scalar_Variable_Coefficient_3D_mu_out0_id0*(-tmp_29*tmp_47 - tmp_31*tmp_47 - tmp_33*tmp_48 - tmp_36*tmp_48 - tmp_38*tmp_48 - tmp_40*tmp_49 - tmp_43*tmp_49 - tmp_45*tmp_49 - tmp_47*tmp_8)*std::abs(tmp_15)*std::abs(p_affine_0_0*tmp_52 - p_affine_0_0*tmp_59 + p_affine_0_1*tmp_55 - p_affine_0_1*tmp_60 + p_affine_0_2*tmp_58 - p_affine_0_2*tmp_61 - p_affine_1_0*tmp_52 + p_affine_1_0*tmp_59 - p_affine_1_1*tmp_55 + p_affine_1_1*tmp_60 - p_affine_1_2*tmp_58 + p_affine_1_2*tmp_61 + p_affine_2_0*tmp_54 - p_affine_2_0*tmp_57 - p_affine_2_1*tmp_51 + p_affine_2_1*tmp_56 + p_affine_2_2*tmp_50 - p_affine_2_2*tmp_53 - p_affine_3_0*tmp_54 + p_affine_3_0*tmp_57 + p_affine_3_1*tmp_51 - p_affine_3_1*tmp_56 - p_affine_3_2*tmp_50 + p_affine_3_2*tmp_53);
      real_t a_0_0 = tmp_62*(-tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46);
      real_t a_0_1 = tmp_62*(tmp_32 + tmp_39 + tmp_46);
      real_t a_0_2 = tmp_62*(tmp_30 + tmp_37 + tmp_44);
      real_t a_0_3 = tmp_62*(tmp_28 + tmp_35 + tmp_42);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_epsilonvar_1_2_blending_q1::Blending_DF_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
   {
      Point3D  mappedPt( {in_0, in_1, in_2} );
      Matrix3r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 0, 2 );
      *out_3 = DPsi( 1, 0 );
      *out_4 = DPsi( 1, 1 );
      *out_5 = DPsi( 1, 2 );
      *out_6 = DPsi( 2, 0 );
      *out_7 = DPsi( 2, 1 );
      *out_8 = DPsi( 2, 2 );
   }

   void p1_epsilonvar_1_2_blending_q1::Blending_F_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2 ) const
   {
      Point3D  in( {in_0, in_1, in_2} );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
      *out_2 = out[2];
   }

   void p1_epsilonvar_1_2_blending_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

   void p1_epsilonvar_2_0_blending_q1::integrateAll( const std::array< Point3D, 3 >& , Matrix< real_t, 3, 3 >&  ) const
   {
      
   }

   void p1_epsilonvar_2_0_blending_q1::integrateRow0( const std::array< Point3D, 3 >& , Matrix< real_t, 1, 3 >&  ) const
   {
      
   }

   void p1_epsilonvar_2_0_blending_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_14 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_15 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out6_id0*tmp_14 + Blending_DF_Tetrahedron_blend_out7_id0*tmp_12 - Blending_DF_Tetrahedron_blend_out7_id0*tmp_13 - Blending_DF_Tetrahedron_blend_out8_id0*tmp_10 + Blending_DF_Tetrahedron_blend_out8_id0*tmp_9;
      real_t tmp_16 = -p_affine_0_2;
      real_t tmp_17 = p_affine_3_2 + tmp_16;
      real_t tmp_18 = p_affine_1_2 + tmp_16;
      real_t tmp_19 = p_affine_3_1 + tmp_2;
      real_t tmp_20 = tmp_19*tmp_5;
      real_t tmp_21 = p_affine_2_2 + tmp_16;
      real_t tmp_22 = p_affine_3_0 + tmp_0;
      real_t tmp_23 = tmp_22*tmp_6;
      real_t tmp_24 = tmp_1*tmp_19;
      real_t tmp_25 = tmp_22*tmp_3;
      real_t tmp_26 = 0.5/(tmp_15*(tmp_17*tmp_4 - tmp_17*tmp_7 + tmp_18*tmp_20 - tmp_18*tmp_25 + tmp_21*tmp_23 - tmp_21*tmp_24));
      real_t tmp_27 = tmp_26*(-tmp_10 + tmp_9);
      real_t tmp_28 = tmp_27*tmp_8;
      real_t tmp_29 = tmp_23 - tmp_24;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = tmp_20 - tmp_25;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = -tmp_1*tmp_21 + tmp_18*tmp_5;
      real_t tmp_34 = tmp_26*(tmp_12 - tmp_13);
      real_t tmp_35 = tmp_33*tmp_34;
      real_t tmp_36 = tmp_1*tmp_17 - tmp_18*tmp_22;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = -tmp_17*tmp_5 + tmp_21*tmp_22;
      real_t tmp_39 = tmp_34*tmp_38;
      real_t tmp_40 = -tmp_18*tmp_3 + tmp_21*tmp_6;
      real_t tmp_41 = tmp_26*(tmp_11 - tmp_14);
      real_t tmp_42 = tmp_40*tmp_41;
      real_t tmp_43 = -tmp_17*tmp_6 + tmp_18*tmp_19;
      real_t tmp_44 = tmp_41*tmp_43;
      real_t tmp_45 = tmp_17*tmp_3 - tmp_19*tmp_21;
      real_t tmp_46 = tmp_41*tmp_45;
      real_t tmp_47 = -tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46;
      real_t tmp_48 = tmp_26*(Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_49 = tmp_48*tmp_8;
      real_t tmp_50 = tmp_29*tmp_48;
      real_t tmp_51 = tmp_31*tmp_48;
      real_t tmp_52 = tmp_26*(-Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_53 = tmp_33*tmp_52;
      real_t tmp_54 = tmp_36*tmp_52;
      real_t tmp_55 = tmp_38*tmp_52;
      real_t tmp_56 = tmp_26*(Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0);
      real_t tmp_57 = tmp_40*tmp_56;
      real_t tmp_58 = tmp_43*tmp_56;
      real_t tmp_59 = tmp_45*tmp_56;
      real_t tmp_60 = p_affine_0_0*p_affine_1_1;
      real_t tmp_61 = p_affine_0_0*p_affine_1_2;
      real_t tmp_62 = p_affine_2_1*p_affine_3_2;
      real_t tmp_63 = p_affine_0_1*p_affine_1_0;
      real_t tmp_64 = p_affine_0_1*p_affine_1_2;
      real_t tmp_65 = p_affine_2_2*p_affine_3_0;
      real_t tmp_66 = p_affine_0_2*p_affine_1_0;
      real_t tmp_67 = p_affine_0_2*p_affine_1_1;
      real_t tmp_68 = p_affine_2_0*p_affine_3_1;
      real_t tmp_69 = p_affine_2_2*p_affine_3_1;
      real_t tmp_70 = p_affine_2_0*p_affine_3_2;
      real_t tmp_71 = p_affine_2_1*p_affine_3_0;
      real_t tmp_72 = 0.66666666666666652*Scalar_Variable_Coefficient_3D_mu_out0_id0*std::abs(tmp_15)*std::abs(p_affine_0_0*tmp_62 - p_affine_0_0*tmp_69 + p_affine_0_1*tmp_65 - p_affine_0_1*tmp_70 + p_affine_0_2*tmp_68 - p_affine_0_2*tmp_71 - p_affine_1_0*tmp_62 + p_affine_1_0*tmp_69 - p_affine_1_1*tmp_65 + p_affine_1_1*tmp_70 - p_affine_1_2*tmp_68 + p_affine_1_2*tmp_71 + p_affine_2_0*tmp_64 - p_affine_2_0*tmp_67 - p_affine_2_1*tmp_61 + p_affine_2_1*tmp_66 + p_affine_2_2*tmp_60 - p_affine_2_2*tmp_63 - p_affine_3_0*tmp_64 + p_affine_3_0*tmp_67 + p_affine_3_1*tmp_61 - p_affine_3_1*tmp_66 - p_affine_3_2*tmp_60 + p_affine_3_2*tmp_63);
      real_t tmp_73 = tmp_72*(-tmp_49 - tmp_50 - tmp_51 - tmp_53 - tmp_54 - tmp_55 - tmp_57 - tmp_58 - tmp_59);
      real_t tmp_74 = tmp_32 + tmp_39 + tmp_46;
      real_t tmp_75 = tmp_30 + tmp_37 + tmp_44;
      real_t tmp_76 = tmp_28 + tmp_35 + tmp_42;
      real_t tmp_77 = tmp_72*(tmp_51 + tmp_55 + tmp_59);
      real_t tmp_78 = tmp_72*(tmp_50 + tmp_54 + tmp_58);
      real_t tmp_79 = tmp_72*(tmp_49 + tmp_53 + tmp_57);
      real_t a_0_0 = tmp_47*tmp_73;
      real_t a_0_1 = tmp_73*tmp_74;
      real_t a_0_2 = tmp_73*tmp_75;
      real_t a_0_3 = tmp_73*tmp_76;
      real_t a_1_0 = tmp_47*tmp_77;
      real_t a_1_1 = tmp_74*tmp_77;
      real_t a_1_2 = tmp_75*tmp_77;
      real_t a_1_3 = tmp_76*tmp_77;
      real_t a_2_0 = tmp_47*tmp_78;
      real_t a_2_1 = tmp_74*tmp_78;
      real_t a_2_2 = tmp_75*tmp_78;
      real_t a_2_3 = tmp_76*tmp_78;
      real_t a_3_0 = tmp_47*tmp_79;
      real_t a_3_1 = tmp_74*tmp_79;
      real_t a_3_2 = tmp_75*tmp_79;
      real_t a_3_3 = tmp_76*tmp_79;
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

   void p1_epsilonvar_2_0_blending_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_14 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_15 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out6_id0*tmp_14 + Blending_DF_Tetrahedron_blend_out7_id0*tmp_12 - Blending_DF_Tetrahedron_blend_out7_id0*tmp_13 - Blending_DF_Tetrahedron_blend_out8_id0*tmp_10 + Blending_DF_Tetrahedron_blend_out8_id0*tmp_9;
      real_t tmp_16 = -p_affine_0_2;
      real_t tmp_17 = p_affine_3_2 + tmp_16;
      real_t tmp_18 = p_affine_1_2 + tmp_16;
      real_t tmp_19 = p_affine_3_1 + tmp_2;
      real_t tmp_20 = tmp_19*tmp_5;
      real_t tmp_21 = p_affine_2_2 + tmp_16;
      real_t tmp_22 = p_affine_3_0 + tmp_0;
      real_t tmp_23 = tmp_22*tmp_6;
      real_t tmp_24 = tmp_1*tmp_19;
      real_t tmp_25 = tmp_22*tmp_3;
      real_t tmp_26 = 0.5/(tmp_15*(tmp_17*tmp_4 - tmp_17*tmp_7 + tmp_18*tmp_20 - tmp_18*tmp_25 + tmp_21*tmp_23 - tmp_21*tmp_24));
      real_t tmp_27 = tmp_26*(-tmp_10 + tmp_9);
      real_t tmp_28 = tmp_27*tmp_8;
      real_t tmp_29 = tmp_23 - tmp_24;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = tmp_20 - tmp_25;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = -tmp_1*tmp_21 + tmp_18*tmp_5;
      real_t tmp_34 = tmp_26*(tmp_12 - tmp_13);
      real_t tmp_35 = tmp_33*tmp_34;
      real_t tmp_36 = tmp_1*tmp_17 - tmp_18*tmp_22;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = -tmp_17*tmp_5 + tmp_21*tmp_22;
      real_t tmp_39 = tmp_34*tmp_38;
      real_t tmp_40 = -tmp_18*tmp_3 + tmp_21*tmp_6;
      real_t tmp_41 = tmp_26*(tmp_11 - tmp_14);
      real_t tmp_42 = tmp_40*tmp_41;
      real_t tmp_43 = -tmp_17*tmp_6 + tmp_18*tmp_19;
      real_t tmp_44 = tmp_41*tmp_43;
      real_t tmp_45 = tmp_17*tmp_3 - tmp_19*tmp_21;
      real_t tmp_46 = tmp_41*tmp_45;
      real_t tmp_47 = tmp_26*(Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_48 = tmp_26*(-Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_49 = tmp_26*(Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0);
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
      real_t tmp_62 = 0.66666666666666652*Scalar_Variable_Coefficient_3D_mu_out0_id0*(-tmp_29*tmp_47 - tmp_31*tmp_47 - tmp_33*tmp_48 - tmp_36*tmp_48 - tmp_38*tmp_48 - tmp_40*tmp_49 - tmp_43*tmp_49 - tmp_45*tmp_49 - tmp_47*tmp_8)*std::abs(tmp_15)*std::abs(p_affine_0_0*tmp_52 - p_affine_0_0*tmp_59 + p_affine_0_1*tmp_55 - p_affine_0_1*tmp_60 + p_affine_0_2*tmp_58 - p_affine_0_2*tmp_61 - p_affine_1_0*tmp_52 + p_affine_1_0*tmp_59 - p_affine_1_1*tmp_55 + p_affine_1_1*tmp_60 - p_affine_1_2*tmp_58 + p_affine_1_2*tmp_61 + p_affine_2_0*tmp_54 - p_affine_2_0*tmp_57 - p_affine_2_1*tmp_51 + p_affine_2_1*tmp_56 + p_affine_2_2*tmp_50 - p_affine_2_2*tmp_53 - p_affine_3_0*tmp_54 + p_affine_3_0*tmp_57 + p_affine_3_1*tmp_51 - p_affine_3_1*tmp_56 - p_affine_3_2*tmp_50 + p_affine_3_2*tmp_53);
      real_t a_0_0 = tmp_62*(-tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46);
      real_t a_0_1 = tmp_62*(tmp_32 + tmp_39 + tmp_46);
      real_t a_0_2 = tmp_62*(tmp_30 + tmp_37 + tmp_44);
      real_t a_0_3 = tmp_62*(tmp_28 + tmp_35 + tmp_42);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_epsilonvar_2_0_blending_q1::Blending_DF_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
   {
      Point3D  mappedPt( {in_0, in_1, in_2} );
      Matrix3r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 0, 2 );
      *out_3 = DPsi( 1, 0 );
      *out_4 = DPsi( 1, 1 );
      *out_5 = DPsi( 1, 2 );
      *out_6 = DPsi( 2, 0 );
      *out_7 = DPsi( 2, 1 );
      *out_8 = DPsi( 2, 2 );
   }

   void p1_epsilonvar_2_0_blending_q1::Blending_F_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2 ) const
   {
      Point3D  in( {in_0, in_1, in_2} );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
      *out_2 = out[2];
   }

   void p1_epsilonvar_2_0_blending_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

   void p1_epsilonvar_2_1_blending_q1::integrateAll( const std::array< Point3D, 3 >& , Matrix< real_t, 3, 3 >&  ) const
   {
      
   }

   void p1_epsilonvar_2_1_blending_q1::integrateRow0( const std::array< Point3D, 3 >& , Matrix< real_t, 1, 3 >&  ) const
   {
      
   }

   void p1_epsilonvar_2_1_blending_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_14 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_15 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out6_id0*tmp_14 + Blending_DF_Tetrahedron_blend_out7_id0*tmp_12 - Blending_DF_Tetrahedron_blend_out7_id0*tmp_13 - Blending_DF_Tetrahedron_blend_out8_id0*tmp_10 + Blending_DF_Tetrahedron_blend_out8_id0*tmp_9;
      real_t tmp_16 = -p_affine_0_2;
      real_t tmp_17 = p_affine_3_2 + tmp_16;
      real_t tmp_18 = p_affine_1_2 + tmp_16;
      real_t tmp_19 = p_affine_3_1 + tmp_2;
      real_t tmp_20 = tmp_19*tmp_5;
      real_t tmp_21 = p_affine_2_2 + tmp_16;
      real_t tmp_22 = p_affine_3_0 + tmp_0;
      real_t tmp_23 = tmp_22*tmp_6;
      real_t tmp_24 = tmp_1*tmp_19;
      real_t tmp_25 = tmp_22*tmp_3;
      real_t tmp_26 = 0.5/(tmp_15*(tmp_17*tmp_4 - tmp_17*tmp_7 + tmp_18*tmp_20 - tmp_18*tmp_25 + tmp_21*tmp_23 - tmp_21*tmp_24));
      real_t tmp_27 = tmp_26*(-tmp_10 + tmp_9);
      real_t tmp_28 = tmp_27*tmp_8;
      real_t tmp_29 = tmp_23 - tmp_24;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = tmp_20 - tmp_25;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = -tmp_1*tmp_21 + tmp_18*tmp_5;
      real_t tmp_34 = tmp_26*(tmp_12 - tmp_13);
      real_t tmp_35 = tmp_33*tmp_34;
      real_t tmp_36 = tmp_1*tmp_17 - tmp_18*tmp_22;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = -tmp_17*tmp_5 + tmp_21*tmp_22;
      real_t tmp_39 = tmp_34*tmp_38;
      real_t tmp_40 = -tmp_18*tmp_3 + tmp_21*tmp_6;
      real_t tmp_41 = tmp_26*(tmp_11 - tmp_14);
      real_t tmp_42 = tmp_40*tmp_41;
      real_t tmp_43 = -tmp_17*tmp_6 + tmp_18*tmp_19;
      real_t tmp_44 = tmp_41*tmp_43;
      real_t tmp_45 = tmp_17*tmp_3 - tmp_19*tmp_21;
      real_t tmp_46 = tmp_41*tmp_45;
      real_t tmp_47 = -tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46;
      real_t tmp_48 = tmp_26*(-Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out7_id0 + Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_49 = tmp_48*tmp_8;
      real_t tmp_50 = tmp_29*tmp_48;
      real_t tmp_51 = tmp_31*tmp_48;
      real_t tmp_52 = tmp_26*(Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_53 = tmp_33*tmp_52;
      real_t tmp_54 = tmp_36*tmp_52;
      real_t tmp_55 = tmp_38*tmp_52;
      real_t tmp_56 = tmp_26*(-Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out7_id0);
      real_t tmp_57 = tmp_40*tmp_56;
      real_t tmp_58 = tmp_43*tmp_56;
      real_t tmp_59 = tmp_45*tmp_56;
      real_t tmp_60 = p_affine_0_0*p_affine_1_1;
      real_t tmp_61 = p_affine_0_0*p_affine_1_2;
      real_t tmp_62 = p_affine_2_1*p_affine_3_2;
      real_t tmp_63 = p_affine_0_1*p_affine_1_0;
      real_t tmp_64 = p_affine_0_1*p_affine_1_2;
      real_t tmp_65 = p_affine_2_2*p_affine_3_0;
      real_t tmp_66 = p_affine_0_2*p_affine_1_0;
      real_t tmp_67 = p_affine_0_2*p_affine_1_1;
      real_t tmp_68 = p_affine_2_0*p_affine_3_1;
      real_t tmp_69 = p_affine_2_2*p_affine_3_1;
      real_t tmp_70 = p_affine_2_0*p_affine_3_2;
      real_t tmp_71 = p_affine_2_1*p_affine_3_0;
      real_t tmp_72 = 0.66666666666666652*Scalar_Variable_Coefficient_3D_mu_out0_id0*std::abs(tmp_15)*std::abs(p_affine_0_0*tmp_62 - p_affine_0_0*tmp_69 + p_affine_0_1*tmp_65 - p_affine_0_1*tmp_70 + p_affine_0_2*tmp_68 - p_affine_0_2*tmp_71 - p_affine_1_0*tmp_62 + p_affine_1_0*tmp_69 - p_affine_1_1*tmp_65 + p_affine_1_1*tmp_70 - p_affine_1_2*tmp_68 + p_affine_1_2*tmp_71 + p_affine_2_0*tmp_64 - p_affine_2_0*tmp_67 - p_affine_2_1*tmp_61 + p_affine_2_1*tmp_66 + p_affine_2_2*tmp_60 - p_affine_2_2*tmp_63 - p_affine_3_0*tmp_64 + p_affine_3_0*tmp_67 + p_affine_3_1*tmp_61 - p_affine_3_1*tmp_66 - p_affine_3_2*tmp_60 + p_affine_3_2*tmp_63);
      real_t tmp_73 = tmp_72*(-tmp_49 - tmp_50 - tmp_51 - tmp_53 - tmp_54 - tmp_55 - tmp_57 - tmp_58 - tmp_59);
      real_t tmp_74 = tmp_32 + tmp_39 + tmp_46;
      real_t tmp_75 = tmp_30 + tmp_37 + tmp_44;
      real_t tmp_76 = tmp_28 + tmp_35 + tmp_42;
      real_t tmp_77 = tmp_72*(tmp_51 + tmp_55 + tmp_59);
      real_t tmp_78 = tmp_72*(tmp_50 + tmp_54 + tmp_58);
      real_t tmp_79 = tmp_72*(tmp_49 + tmp_53 + tmp_57);
      real_t a_0_0 = tmp_47*tmp_73;
      real_t a_0_1 = tmp_73*tmp_74;
      real_t a_0_2 = tmp_73*tmp_75;
      real_t a_0_3 = tmp_73*tmp_76;
      real_t a_1_0 = tmp_47*tmp_77;
      real_t a_1_1 = tmp_74*tmp_77;
      real_t a_1_2 = tmp_75*tmp_77;
      real_t a_1_3 = tmp_76*tmp_77;
      real_t a_2_0 = tmp_47*tmp_78;
      real_t a_2_1 = tmp_74*tmp_78;
      real_t a_2_2 = tmp_75*tmp_78;
      real_t a_2_3 = tmp_76*tmp_78;
      real_t a_3_0 = tmp_47*tmp_79;
      real_t a_3_1 = tmp_74*tmp_79;
      real_t a_3_2 = tmp_75*tmp_79;
      real_t a_3_3 = tmp_76*tmp_79;
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

   void p1_epsilonvar_2_1_blending_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_14 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_15 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out6_id0*tmp_14 + Blending_DF_Tetrahedron_blend_out7_id0*tmp_12 - Blending_DF_Tetrahedron_blend_out7_id0*tmp_13 - Blending_DF_Tetrahedron_blend_out8_id0*tmp_10 + Blending_DF_Tetrahedron_blend_out8_id0*tmp_9;
      real_t tmp_16 = -p_affine_0_2;
      real_t tmp_17 = p_affine_3_2 + tmp_16;
      real_t tmp_18 = p_affine_1_2 + tmp_16;
      real_t tmp_19 = p_affine_3_1 + tmp_2;
      real_t tmp_20 = tmp_19*tmp_5;
      real_t tmp_21 = p_affine_2_2 + tmp_16;
      real_t tmp_22 = p_affine_3_0 + tmp_0;
      real_t tmp_23 = tmp_22*tmp_6;
      real_t tmp_24 = tmp_1*tmp_19;
      real_t tmp_25 = tmp_22*tmp_3;
      real_t tmp_26 = 0.5/(tmp_15*(tmp_17*tmp_4 - tmp_17*tmp_7 + tmp_18*tmp_20 - tmp_18*tmp_25 + tmp_21*tmp_23 - tmp_21*tmp_24));
      real_t tmp_27 = tmp_26*(-tmp_10 + tmp_9);
      real_t tmp_28 = tmp_27*tmp_8;
      real_t tmp_29 = tmp_23 - tmp_24;
      real_t tmp_30 = tmp_27*tmp_29;
      real_t tmp_31 = tmp_20 - tmp_25;
      real_t tmp_32 = tmp_27*tmp_31;
      real_t tmp_33 = -tmp_1*tmp_21 + tmp_18*tmp_5;
      real_t tmp_34 = tmp_26*(tmp_12 - tmp_13);
      real_t tmp_35 = tmp_33*tmp_34;
      real_t tmp_36 = tmp_1*tmp_17 - tmp_18*tmp_22;
      real_t tmp_37 = tmp_34*tmp_36;
      real_t tmp_38 = -tmp_17*tmp_5 + tmp_21*tmp_22;
      real_t tmp_39 = tmp_34*tmp_38;
      real_t tmp_40 = -tmp_18*tmp_3 + tmp_21*tmp_6;
      real_t tmp_41 = tmp_26*(tmp_11 - tmp_14);
      real_t tmp_42 = tmp_40*tmp_41;
      real_t tmp_43 = -tmp_17*tmp_6 + tmp_18*tmp_19;
      real_t tmp_44 = tmp_41*tmp_43;
      real_t tmp_45 = tmp_17*tmp_3 - tmp_19*tmp_21;
      real_t tmp_46 = tmp_41*tmp_45;
      real_t tmp_47 = tmp_26*(-Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out7_id0 + Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_48 = tmp_26*(Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_49 = tmp_26*(-Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out7_id0);
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
      real_t tmp_62 = 0.66666666666666652*Scalar_Variable_Coefficient_3D_mu_out0_id0*(-tmp_29*tmp_47 - tmp_31*tmp_47 - tmp_33*tmp_48 - tmp_36*tmp_48 - tmp_38*tmp_48 - tmp_40*tmp_49 - tmp_43*tmp_49 - tmp_45*tmp_49 - tmp_47*tmp_8)*std::abs(tmp_15)*std::abs(p_affine_0_0*tmp_52 - p_affine_0_0*tmp_59 + p_affine_0_1*tmp_55 - p_affine_0_1*tmp_60 + p_affine_0_2*tmp_58 - p_affine_0_2*tmp_61 - p_affine_1_0*tmp_52 + p_affine_1_0*tmp_59 - p_affine_1_1*tmp_55 + p_affine_1_1*tmp_60 - p_affine_1_2*tmp_58 + p_affine_1_2*tmp_61 + p_affine_2_0*tmp_54 - p_affine_2_0*tmp_57 - p_affine_2_1*tmp_51 + p_affine_2_1*tmp_56 + p_affine_2_2*tmp_50 - p_affine_2_2*tmp_53 - p_affine_3_0*tmp_54 + p_affine_3_0*tmp_57 + p_affine_3_1*tmp_51 - p_affine_3_1*tmp_56 - p_affine_3_2*tmp_50 + p_affine_3_2*tmp_53);
      real_t a_0_0 = tmp_62*(-tmp_28 - tmp_30 - tmp_32 - tmp_35 - tmp_37 - tmp_39 - tmp_42 - tmp_44 - tmp_46);
      real_t a_0_1 = tmp_62*(tmp_32 + tmp_39 + tmp_46);
      real_t a_0_2 = tmp_62*(tmp_30 + tmp_37 + tmp_44);
      real_t a_0_3 = tmp_62*(tmp_28 + tmp_35 + tmp_42);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_epsilonvar_2_1_blending_q1::Blending_DF_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
   {
      Point3D  mappedPt( {in_0, in_1, in_2} );
      Matrix3r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 0, 2 );
      *out_3 = DPsi( 1, 0 );
      *out_4 = DPsi( 1, 1 );
      *out_5 = DPsi( 1, 2 );
      *out_6 = DPsi( 2, 0 );
      *out_7 = DPsi( 2, 1 );
      *out_8 = DPsi( 2, 2 );
   }

   void p1_epsilonvar_2_1_blending_q1::Blending_F_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2 ) const
   {
      Point3D  in( {in_0, in_1, in_2} );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
      *out_2 = out[2];
   }

   void p1_epsilonvar_2_1_blending_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

   void p1_epsilonvar_2_2_blending_q1::integrateAll( const std::array< Point3D, 3 >& , Matrix< real_t, 3, 3 >&  ) const
   {
      
   }

   void p1_epsilonvar_2_2_blending_q1::integrateRow0( const std::array< Point3D, 3 >& , Matrix< real_t, 1, 3 >&  ) const
   {
      
   }

   void p1_epsilonvar_2_2_blending_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_14 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_15 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out6_id0*tmp_14 + Blending_DF_Tetrahedron_blend_out7_id0*tmp_12 - Blending_DF_Tetrahedron_blend_out7_id0*tmp_13 - Blending_DF_Tetrahedron_blend_out8_id0*tmp_10 + Blending_DF_Tetrahedron_blend_out8_id0*tmp_9;
      real_t tmp_16 = -p_affine_0_2;
      real_t tmp_17 = p_affine_3_2 + tmp_16;
      real_t tmp_18 = p_affine_1_2 + tmp_16;
      real_t tmp_19 = p_affine_3_1 + tmp_2;
      real_t tmp_20 = tmp_19*tmp_5;
      real_t tmp_21 = p_affine_2_2 + tmp_16;
      real_t tmp_22 = p_affine_3_0 + tmp_0;
      real_t tmp_23 = tmp_22*tmp_6;
      real_t tmp_24 = tmp_1*tmp_19;
      real_t tmp_25 = tmp_22*tmp_3;
      real_t tmp_26 = 1/(tmp_15*(tmp_17*tmp_4 - tmp_17*tmp_7 + tmp_18*tmp_20 - tmp_18*tmp_25 + tmp_21*tmp_23 - tmp_21*tmp_24));
      real_t tmp_27 = 1.0*tmp_26;
      real_t tmp_28 = tmp_27*(-tmp_10 + tmp_9);
      real_t tmp_29 = tmp_28*tmp_8;
      real_t tmp_30 = tmp_23 - tmp_24;
      real_t tmp_31 = tmp_28*tmp_30;
      real_t tmp_32 = tmp_20 - tmp_25;
      real_t tmp_33 = tmp_28*tmp_32;
      real_t tmp_34 = -tmp_1*tmp_21 + tmp_18*tmp_5;
      real_t tmp_35 = tmp_27*(tmp_12 - tmp_13);
      real_t tmp_36 = tmp_34*tmp_35;
      real_t tmp_37 = tmp_1*tmp_17 - tmp_18*tmp_22;
      real_t tmp_38 = tmp_35*tmp_37;
      real_t tmp_39 = -tmp_17*tmp_5 + tmp_21*tmp_22;
      real_t tmp_40 = tmp_35*tmp_39;
      real_t tmp_41 = -tmp_18*tmp_3 + tmp_21*tmp_6;
      real_t tmp_42 = tmp_27*(tmp_11 - tmp_14);
      real_t tmp_43 = tmp_41*tmp_42;
      real_t tmp_44 = -tmp_17*tmp_6 + tmp_18*tmp_19;
      real_t tmp_45 = tmp_42*tmp_44;
      real_t tmp_46 = tmp_17*tmp_3 - tmp_19*tmp_21;
      real_t tmp_47 = tmp_42*tmp_46;
      real_t tmp_48 = -tmp_29 - tmp_31 - tmp_33 - tmp_36 - tmp_38 - tmp_40 - tmp_43 - tmp_45 - tmp_47;
      real_t tmp_49 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_50 = tmp_26*(-Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out7_id0 + Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_51 = tmp_50*tmp_8;
      real_t tmp_52 = tmp_30*tmp_50;
      real_t tmp_53 = tmp_32*tmp_50;
      real_t tmp_54 = tmp_26*(Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_55 = tmp_34*tmp_54;
      real_t tmp_56 = tmp_37*tmp_54;
      real_t tmp_57 = tmp_39*tmp_54;
      real_t tmp_58 = tmp_26*(-Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out7_id0);
      real_t tmp_59 = tmp_41*tmp_58;
      real_t tmp_60 = tmp_44*tmp_58;
      real_t tmp_61 = tmp_46*tmp_58;
      real_t tmp_62 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_63 = tmp_26*(Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_64 = tmp_63*tmp_8;
      real_t tmp_65 = tmp_30*tmp_63;
      real_t tmp_66 = tmp_32*tmp_63;
      real_t tmp_67 = tmp_26*(-Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_68 = tmp_34*tmp_67;
      real_t tmp_69 = tmp_37*tmp_67;
      real_t tmp_70 = tmp_39*tmp_67;
      real_t tmp_71 = tmp_26*(Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0);
      real_t tmp_72 = tmp_41*tmp_71;
      real_t tmp_73 = tmp_44*tmp_71;
      real_t tmp_74 = tmp_46*tmp_71;
      real_t tmp_75 = p_affine_0_0*p_affine_1_1;
      real_t tmp_76 = p_affine_0_0*p_affine_1_2;
      real_t tmp_77 = p_affine_2_1*p_affine_3_2;
      real_t tmp_78 = p_affine_0_1*p_affine_1_0;
      real_t tmp_79 = p_affine_0_1*p_affine_1_2;
      real_t tmp_80 = p_affine_2_2*p_affine_3_0;
      real_t tmp_81 = p_affine_0_2*p_affine_1_0;
      real_t tmp_82 = p_affine_0_2*p_affine_1_1;
      real_t tmp_83 = p_affine_2_0*p_affine_3_1;
      real_t tmp_84 = p_affine_2_2*p_affine_3_1;
      real_t tmp_85 = p_affine_2_0*p_affine_3_2;
      real_t tmp_86 = p_affine_2_1*p_affine_3_0;
      real_t tmp_87 = 0.16666666666666663*std::abs(tmp_15)*std::abs(p_affine_0_0*tmp_77 - p_affine_0_0*tmp_84 + p_affine_0_1*tmp_80 - p_affine_0_1*tmp_85 + p_affine_0_2*tmp_83 - p_affine_0_2*tmp_86 - p_affine_1_0*tmp_77 + p_affine_1_0*tmp_84 - p_affine_1_1*tmp_80 + p_affine_1_1*tmp_85 - p_affine_1_2*tmp_83 + p_affine_1_2*tmp_86 + p_affine_2_0*tmp_79 - p_affine_2_0*tmp_82 - p_affine_2_1*tmp_76 + p_affine_2_1*tmp_81 + p_affine_2_2*tmp_75 - p_affine_2_2*tmp_78 - p_affine_3_0*tmp_79 + p_affine_3_0*tmp_82 + p_affine_3_1*tmp_76 - p_affine_3_1*tmp_81 - p_affine_3_2*tmp_75 + p_affine_3_2*tmp_78);
      real_t tmp_88 = tmp_33 + tmp_40 + tmp_47;
      real_t tmp_89 = tmp_48*tmp_49;
      real_t tmp_90 = 0.5*tmp_53;
      real_t tmp_91 = 0.5*tmp_57;
      real_t tmp_92 = 0.5*tmp_61;
      real_t tmp_93 = tmp_90 + tmp_91 + tmp_92;
      real_t tmp_94 = 0.5*tmp_51;
      real_t tmp_95 = 0.5*tmp_52;
      real_t tmp_96 = 0.5*tmp_55;
      real_t tmp_97 = 0.5*tmp_56;
      real_t tmp_98 = 0.5*tmp_59;
      real_t tmp_99 = 0.5*tmp_60;
      real_t tmp_100 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_101 = tmp_100*(-tmp_90 - tmp_91 - tmp_92 - tmp_94 - tmp_95 - tmp_96 - tmp_97 - tmp_98 - tmp_99);
      real_t tmp_102 = 0.5*tmp_66;
      real_t tmp_103 = 0.5*tmp_70;
      real_t tmp_104 = 0.5*tmp_74;
      real_t tmp_105 = tmp_102 + tmp_103 + tmp_104;
      real_t tmp_106 = 0.5*tmp_64;
      real_t tmp_107 = 0.5*tmp_65;
      real_t tmp_108 = 0.5*tmp_68;
      real_t tmp_109 = 0.5*tmp_69;
      real_t tmp_110 = 0.5*tmp_72;
      real_t tmp_111 = 0.5*tmp_73;
      real_t tmp_112 = tmp_100*(-tmp_102 - tmp_103 - tmp_104 - tmp_106 - tmp_107 - tmp_108 - tmp_109 - tmp_110 - tmp_111);
      real_t tmp_113 = tmp_87*(tmp_101*tmp_93 + tmp_105*tmp_112 + tmp_88*tmp_89);
      real_t tmp_114 = tmp_31 + tmp_38 + tmp_45;
      real_t tmp_115 = tmp_95 + tmp_97 + tmp_99;
      real_t tmp_116 = tmp_107 + tmp_109 + tmp_111;
      real_t tmp_117 = tmp_87*(tmp_101*tmp_115 + tmp_112*tmp_116 + tmp_114*tmp_89);
      real_t tmp_118 = tmp_29 + tmp_36 + tmp_43;
      real_t tmp_119 = tmp_94 + tmp_96 + tmp_98;
      real_t tmp_120 = tmp_106 + tmp_108 + tmp_110;
      real_t tmp_121 = tmp_87*(tmp_101*tmp_119 + tmp_112*tmp_120 + tmp_118*tmp_89);
      real_t tmp_122 = tmp_49*tmp_88;
      real_t tmp_123 = tmp_100*tmp_93;
      real_t tmp_124 = tmp_100*tmp_105;
      real_t tmp_125 = tmp_87*(tmp_114*tmp_122 + tmp_115*tmp_123 + tmp_116*tmp_124);
      real_t tmp_126 = tmp_87*(tmp_118*tmp_122 + tmp_119*tmp_123 + tmp_120*tmp_124);
      real_t tmp_127 = tmp_87*(tmp_100*tmp_115*tmp_119 + tmp_100*tmp_116*tmp_120 + tmp_114*tmp_118*tmp_49);
      real_t a_0_0 = tmp_87*((tmp_48*tmp_48)*tmp_49 + tmp_62*((-tmp_51 - tmp_52 - tmp_53 - tmp_55 - tmp_56 - tmp_57 - tmp_59 - tmp_60 - tmp_61)*(-tmp_51 - tmp_52 - tmp_53 - tmp_55 - tmp_56 - tmp_57 - tmp_59 - tmp_60 - tmp_61)) + tmp_62*((-tmp_64 - tmp_65 - tmp_66 - tmp_68 - tmp_69 - tmp_70 - tmp_72 - tmp_73 - tmp_74)*(-tmp_64 - tmp_65 - tmp_66 - tmp_68 - tmp_69 - tmp_70 - tmp_72 - tmp_73 - tmp_74)));
      real_t a_0_1 = tmp_113;
      real_t a_0_2 = tmp_117;
      real_t a_0_3 = tmp_121;
      real_t a_1_0 = tmp_113;
      real_t a_1_1 = tmp_87*(tmp_49*(tmp_88*tmp_88) + tmp_62*((tmp_53 + tmp_57 + tmp_61)*(tmp_53 + tmp_57 + tmp_61)) + tmp_62*((tmp_66 + tmp_70 + tmp_74)*(tmp_66 + tmp_70 + tmp_74)));
      real_t a_1_2 = tmp_125;
      real_t a_1_3 = tmp_126;
      real_t a_2_0 = tmp_117;
      real_t a_2_1 = tmp_125;
      real_t a_2_2 = tmp_87*((tmp_114*tmp_114)*tmp_49 + tmp_62*((tmp_52 + tmp_56 + tmp_60)*(tmp_52 + tmp_56 + tmp_60)) + tmp_62*((tmp_65 + tmp_69 + tmp_73)*(tmp_65 + tmp_69 + tmp_73)));
      real_t a_2_3 = tmp_127;
      real_t a_3_0 = tmp_121;
      real_t a_3_1 = tmp_126;
      real_t a_3_2 = tmp_127;
      real_t a_3_3 = tmp_87*((tmp_118*tmp_118)*tmp_49 + tmp_62*((tmp_51 + tmp_55 + tmp_59)*(tmp_51 + tmp_55 + tmp_59)) + tmp_62*((tmp_64 + tmp_68 + tmp_72)*(tmp_64 + tmp_68 + tmp_72)));
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

   void p1_epsilonvar_2_2_blending_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_mu_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_mu( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_mu_out0_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_14 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_15 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out6_id0*tmp_14 + Blending_DF_Tetrahedron_blend_out7_id0*tmp_12 - Blending_DF_Tetrahedron_blend_out7_id0*tmp_13 - Blending_DF_Tetrahedron_blend_out8_id0*tmp_10 + Blending_DF_Tetrahedron_blend_out8_id0*tmp_9;
      real_t tmp_16 = -p_affine_0_2;
      real_t tmp_17 = p_affine_3_2 + tmp_16;
      real_t tmp_18 = p_affine_1_2 + tmp_16;
      real_t tmp_19 = p_affine_3_1 + tmp_2;
      real_t tmp_20 = tmp_19*tmp_5;
      real_t tmp_21 = p_affine_2_2 + tmp_16;
      real_t tmp_22 = p_affine_3_0 + tmp_0;
      real_t tmp_23 = tmp_22*tmp_6;
      real_t tmp_24 = tmp_1*tmp_19;
      real_t tmp_25 = tmp_22*tmp_3;
      real_t tmp_26 = 1/(tmp_15*(tmp_17*tmp_4 - tmp_17*tmp_7 + tmp_18*tmp_20 - tmp_18*tmp_25 + tmp_21*tmp_23 - tmp_21*tmp_24));
      real_t tmp_27 = 1.0*tmp_26;
      real_t tmp_28 = tmp_27*(-tmp_10 + tmp_9);
      real_t tmp_29 = tmp_28*tmp_8;
      real_t tmp_30 = tmp_23 - tmp_24;
      real_t tmp_31 = tmp_28*tmp_30;
      real_t tmp_32 = tmp_20 - tmp_25;
      real_t tmp_33 = tmp_28*tmp_32;
      real_t tmp_34 = -tmp_1*tmp_21 + tmp_18*tmp_5;
      real_t tmp_35 = tmp_27*(tmp_12 - tmp_13);
      real_t tmp_36 = tmp_34*tmp_35;
      real_t tmp_37 = tmp_1*tmp_17 - tmp_18*tmp_22;
      real_t tmp_38 = tmp_35*tmp_37;
      real_t tmp_39 = -tmp_17*tmp_5 + tmp_21*tmp_22;
      real_t tmp_40 = tmp_35*tmp_39;
      real_t tmp_41 = -tmp_18*tmp_3 + tmp_21*tmp_6;
      real_t tmp_42 = tmp_27*(tmp_11 - tmp_14);
      real_t tmp_43 = tmp_41*tmp_42;
      real_t tmp_44 = -tmp_17*tmp_6 + tmp_18*tmp_19;
      real_t tmp_45 = tmp_42*tmp_44;
      real_t tmp_46 = tmp_17*tmp_3 - tmp_19*tmp_21;
      real_t tmp_47 = tmp_42*tmp_46;
      real_t tmp_48 = -tmp_29 - tmp_31 - tmp_33 - tmp_36 - tmp_38 - tmp_40 - tmp_43 - tmp_45 - tmp_47;
      real_t tmp_49 = 2*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_50 = tmp_26*(-Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out7_id0 + Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_51 = tmp_50*tmp_8;
      real_t tmp_52 = tmp_30*tmp_50;
      real_t tmp_53 = tmp_32*tmp_50;
      real_t tmp_54 = tmp_26*(Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_55 = tmp_34*tmp_54;
      real_t tmp_56 = tmp_37*tmp_54;
      real_t tmp_57 = tmp_39*tmp_54;
      real_t tmp_58 = tmp_26*(-Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out7_id0);
      real_t tmp_59 = tmp_41*tmp_58;
      real_t tmp_60 = tmp_44*tmp_58;
      real_t tmp_61 = tmp_46*tmp_58;
      real_t tmp_62 = 1.0*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_63 = tmp_26*(Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_64 = tmp_63*tmp_8;
      real_t tmp_65 = tmp_30*tmp_63;
      real_t tmp_66 = tmp_32*tmp_63;
      real_t tmp_67 = tmp_26*(-Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0);
      real_t tmp_68 = tmp_34*tmp_67;
      real_t tmp_69 = tmp_37*tmp_67;
      real_t tmp_70 = tmp_39*tmp_67;
      real_t tmp_71 = tmp_26*(Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0);
      real_t tmp_72 = tmp_41*tmp_71;
      real_t tmp_73 = tmp_44*tmp_71;
      real_t tmp_74 = tmp_46*tmp_71;
      real_t tmp_75 = p_affine_0_0*p_affine_1_1;
      real_t tmp_76 = p_affine_0_0*p_affine_1_2;
      real_t tmp_77 = p_affine_2_1*p_affine_3_2;
      real_t tmp_78 = p_affine_0_1*p_affine_1_0;
      real_t tmp_79 = p_affine_0_1*p_affine_1_2;
      real_t tmp_80 = p_affine_2_2*p_affine_3_0;
      real_t tmp_81 = p_affine_0_2*p_affine_1_0;
      real_t tmp_82 = p_affine_0_2*p_affine_1_1;
      real_t tmp_83 = p_affine_2_0*p_affine_3_1;
      real_t tmp_84 = p_affine_2_2*p_affine_3_1;
      real_t tmp_85 = p_affine_2_0*p_affine_3_2;
      real_t tmp_86 = p_affine_2_1*p_affine_3_0;
      real_t tmp_87 = 0.16666666666666663*std::abs(tmp_15)*std::abs(p_affine_0_0*tmp_77 - p_affine_0_0*tmp_84 + p_affine_0_1*tmp_80 - p_affine_0_1*tmp_85 + p_affine_0_2*tmp_83 - p_affine_0_2*tmp_86 - p_affine_1_0*tmp_77 + p_affine_1_0*tmp_84 - p_affine_1_1*tmp_80 + p_affine_1_1*tmp_85 - p_affine_1_2*tmp_83 + p_affine_1_2*tmp_86 + p_affine_2_0*tmp_79 - p_affine_2_0*tmp_82 - p_affine_2_1*tmp_76 + p_affine_2_1*tmp_81 + p_affine_2_2*tmp_75 - p_affine_2_2*tmp_78 - p_affine_3_0*tmp_79 + p_affine_3_0*tmp_82 + p_affine_3_1*tmp_76 - p_affine_3_1*tmp_81 - p_affine_3_2*tmp_75 + p_affine_3_2*tmp_78);
      real_t tmp_88 = tmp_48*tmp_49;
      real_t tmp_89 = 0.5*tmp_53;
      real_t tmp_90 = 0.5*tmp_57;
      real_t tmp_91 = 0.5*tmp_61;
      real_t tmp_92 = 0.5*tmp_51;
      real_t tmp_93 = 0.5*tmp_52;
      real_t tmp_94 = 0.5*tmp_55;
      real_t tmp_95 = 0.5*tmp_56;
      real_t tmp_96 = 0.5*tmp_59;
      real_t tmp_97 = 0.5*tmp_60;
      real_t tmp_98 = 4*Scalar_Variable_Coefficient_3D_mu_out0_id0;
      real_t tmp_99 = tmp_98*(-tmp_89 - tmp_90 - tmp_91 - tmp_92 - tmp_93 - tmp_94 - tmp_95 - tmp_96 - tmp_97);
      real_t tmp_100 = 0.5*tmp_66;
      real_t tmp_101 = 0.5*tmp_70;
      real_t tmp_102 = 0.5*tmp_74;
      real_t tmp_103 = 0.5*tmp_64;
      real_t tmp_104 = 0.5*tmp_65;
      real_t tmp_105 = 0.5*tmp_68;
      real_t tmp_106 = 0.5*tmp_69;
      real_t tmp_107 = 0.5*tmp_72;
      real_t tmp_108 = 0.5*tmp_73;
      real_t tmp_109 = tmp_98*(-tmp_100 - tmp_101 - tmp_102 - tmp_103 - tmp_104 - tmp_105 - tmp_106 - tmp_107 - tmp_108);
      real_t a_0_0 = tmp_87*((tmp_48*tmp_48)*tmp_49 + tmp_62*((-tmp_51 - tmp_52 - tmp_53 - tmp_55 - tmp_56 - tmp_57 - tmp_59 - tmp_60 - tmp_61)*(-tmp_51 - tmp_52 - tmp_53 - tmp_55 - tmp_56 - tmp_57 - tmp_59 - tmp_60 - tmp_61)) + tmp_62*((-tmp_64 - tmp_65 - tmp_66 - tmp_68 - tmp_69 - tmp_70 - tmp_72 - tmp_73 - tmp_74)*(-tmp_64 - tmp_65 - tmp_66 - tmp_68 - tmp_69 - tmp_70 - tmp_72 - tmp_73 - tmp_74)));
      real_t a_0_1 = tmp_87*(tmp_109*(tmp_100 + tmp_101 + tmp_102) + tmp_88*(tmp_33 + tmp_40 + tmp_47) + tmp_99*(tmp_89 + tmp_90 + tmp_91));
      real_t a_0_2 = tmp_87*(tmp_109*(tmp_104 + tmp_106 + tmp_108) + tmp_88*(tmp_31 + tmp_38 + tmp_45) + tmp_99*(tmp_93 + tmp_95 + tmp_97));
      real_t a_0_3 = tmp_87*(tmp_109*(tmp_103 + tmp_105 + tmp_107) + tmp_88*(tmp_29 + tmp_36 + tmp_43) + tmp_99*(tmp_92 + tmp_94 + tmp_96));
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_epsilonvar_2_2_blending_q1::Blending_DF_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
   {
      Point3D  mappedPt( {in_0, in_1, in_2} );
      Matrix3r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 0, 2 );
      *out_3 = DPsi( 1, 0 );
      *out_4 = DPsi( 1, 1 );
      *out_5 = DPsi( 1, 2 );
      *out_6 = DPsi( 2, 0 );
      *out_7 = DPsi( 2, 1 );
      *out_8 = DPsi( 2, 2 );
   }

   void p1_epsilonvar_2_2_blending_q1::Blending_F_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2 ) const
   {
      Point3D  in( {in_0, in_1, in_2} );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
      *out_2 = out[2];
   }

   void p1_epsilonvar_2_2_blending_q1::Scalar_Variable_Coefficient_3D_mu( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_mu( Point3D( {in_0, in_1, in_2} ) );
   }

} // namespace forms
} // namespace hyteg

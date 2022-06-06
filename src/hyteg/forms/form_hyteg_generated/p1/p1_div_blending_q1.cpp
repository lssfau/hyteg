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

#include "p1_div_blending_q1.hpp"

namespace hyteg {
namespace forms {

   void p1_div_0_blending_q1::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
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
      Blending_DF_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_DF_Triangle_blend_out0_id0, &Blending_DF_Triangle_blend_out1_id0, &Blending_DF_Triangle_blend_out2_id0, &Blending_DF_Triangle_blend_out3_id0 );
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = Blending_DF_Triangle_blend_out0_id0*Blending_DF_Triangle_blend_out3_id0 - Blending_DF_Triangle_blend_out1_id0*Blending_DF_Triangle_blend_out2_id0;
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = p_affine_1_0 + tmp_3;
      real_t tmp_5 = 1/(tmp_2*(tmp_1*tmp_4 - (p_affine_1_1 + tmp_0)*(p_affine_2_0 + tmp_3)));
      real_t tmp_6 = Blending_DF_Triangle_blend_out3_id0*tmp_5;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = tmp_6*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_9 = Blending_DF_Triangle_blend_out2_id0*tmp_5;
      real_t tmp_10 = tmp_4*tmp_9;
      real_t tmp_11 = tmp_9*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_12 = -tmp_10 - tmp_11 + tmp_7 + tmp_8;
      real_t tmp_13 = 0.5*std::abs(tmp_2)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_14 = 0.33333333333333343*tmp_13;
      real_t tmp_15 = tmp_11 - tmp_7;
      real_t tmp_16 = tmp_10 - tmp_8;
      real_t tmp_17 = 0.33333333333333331*tmp_13;
      real_t tmp_18 = 0.33333333333333331*tmp_13;
      real_t a_0_0 = tmp_12*tmp_14;
      real_t a_0_1 = tmp_14*tmp_15;
      real_t a_0_2 = tmp_14*tmp_16;
      real_t a_1_0 = tmp_12*tmp_17;
      real_t a_1_1 = tmp_15*tmp_17;
      real_t a_1_2 = tmp_16*tmp_17;
      real_t a_2_0 = tmp_12*tmp_18;
      real_t a_2_1 = tmp_15*tmp_18;
      real_t a_2_2 = tmp_16*tmp_18;
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

   void p1_div_0_blending_q1::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
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
      Blending_DF_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_DF_Triangle_blend_out0_id0, &Blending_DF_Triangle_blend_out1_id0, &Blending_DF_Triangle_blend_out2_id0, &Blending_DF_Triangle_blend_out3_id0 );
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = Blending_DF_Triangle_blend_out0_id0*Blending_DF_Triangle_blend_out3_id0 - Blending_DF_Triangle_blend_out1_id0*Blending_DF_Triangle_blend_out2_id0;
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = p_affine_1_0 + tmp_3;
      real_t tmp_5 = 1/(tmp_2*(tmp_1*tmp_4 - (p_affine_1_1 + tmp_0)*(p_affine_2_0 + tmp_3)));
      real_t tmp_6 = Blending_DF_Triangle_blend_out3_id0*tmp_5;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = tmp_6*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_9 = Blending_DF_Triangle_blend_out2_id0*tmp_5;
      real_t tmp_10 = tmp_4*tmp_9;
      real_t tmp_11 = tmp_9*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_12 = 0.16666666666666671*std::abs(tmp_2)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = tmp_12*(-tmp_10 - tmp_11 + tmp_7 + tmp_8);
      real_t a_0_1 = tmp_12*(tmp_11 - tmp_7);
      real_t a_0_2 = tmp_12*(tmp_10 - tmp_8);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_div_0_blending_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_14 = Blending_DF_Tetrahedron_blend_out0_id0*tmp_10 - Blending_DF_Tetrahedron_blend_out0_id0*tmp_12 + Blending_DF_Tetrahedron_blend_out1_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out1_id0*tmp_13 + Blending_DF_Tetrahedron_blend_out2_id0*tmp_8 - Blending_DF_Tetrahedron_blend_out2_id0*tmp_9;
      real_t tmp_15 = -p_affine_0_2;
      real_t tmp_16 = p_affine_3_2 + tmp_15;
      real_t tmp_17 = p_affine_1_2 + tmp_15;
      real_t tmp_18 = p_affine_3_1 + tmp_2;
      real_t tmp_19 = tmp_18*tmp_5;
      real_t tmp_20 = p_affine_2_2 + tmp_15;
      real_t tmp_21 = p_affine_3_0 + tmp_0;
      real_t tmp_22 = tmp_21*tmp_6;
      real_t tmp_23 = tmp_1*tmp_18;
      real_t tmp_24 = tmp_21*tmp_3;
      real_t tmp_25 = 1/(tmp_14*(tmp_16*tmp_4 - tmp_16*tmp_7 + tmp_17*tmp_19 - tmp_17*tmp_24 + tmp_20*tmp_22 - tmp_20*tmp_23));
      real_t tmp_26 = tmp_25*(tmp_8 - tmp_9);
      real_t tmp_27 = tmp_26*(tmp_4 - tmp_7);
      real_t tmp_28 = tmp_26*(tmp_22 - tmp_23);
      real_t tmp_29 = tmp_26*(tmp_19 - tmp_24);
      real_t tmp_30 = tmp_25*(tmp_11 - tmp_13);
      real_t tmp_31 = tmp_30*(-tmp_1*tmp_20 + tmp_17*tmp_5);
      real_t tmp_32 = tmp_30*(tmp_1*tmp_16 - tmp_17*tmp_21);
      real_t tmp_33 = tmp_30*(-tmp_16*tmp_5 + tmp_20*tmp_21);
      real_t tmp_34 = tmp_25*(tmp_10 - tmp_12);
      real_t tmp_35 = tmp_34*(-tmp_17*tmp_3 + tmp_20*tmp_6);
      real_t tmp_36 = tmp_34*(-tmp_16*tmp_6 + tmp_17*tmp_18);
      real_t tmp_37 = tmp_34*(tmp_16*tmp_3 - tmp_18*tmp_20);
      real_t tmp_38 = tmp_27 + tmp_28 + tmp_29 + tmp_31 + tmp_32 + tmp_33 + tmp_35 + tmp_36 + tmp_37;
      real_t tmp_39 = p_affine_0_0*p_affine_1_1;
      real_t tmp_40 = p_affine_0_0*p_affine_1_2;
      real_t tmp_41 = p_affine_2_1*p_affine_3_2;
      real_t tmp_42 = p_affine_0_1*p_affine_1_0;
      real_t tmp_43 = p_affine_0_1*p_affine_1_2;
      real_t tmp_44 = p_affine_2_2*p_affine_3_0;
      real_t tmp_45 = p_affine_0_2*p_affine_1_0;
      real_t tmp_46 = p_affine_0_2*p_affine_1_1;
      real_t tmp_47 = p_affine_2_0*p_affine_3_1;
      real_t tmp_48 = p_affine_2_2*p_affine_3_1;
      real_t tmp_49 = p_affine_2_0*p_affine_3_2;
      real_t tmp_50 = p_affine_2_1*p_affine_3_0;
      real_t tmp_51 = 0.16666666666666663*std::abs(tmp_14)*std::abs(p_affine_0_0*tmp_41 - p_affine_0_0*tmp_48 + p_affine_0_1*tmp_44 - p_affine_0_1*tmp_49 + p_affine_0_2*tmp_47 - p_affine_0_2*tmp_50 - p_affine_1_0*tmp_41 + p_affine_1_0*tmp_48 - p_affine_1_1*tmp_44 + p_affine_1_1*tmp_49 - p_affine_1_2*tmp_47 + p_affine_1_2*tmp_50 + p_affine_2_0*tmp_43 - p_affine_2_0*tmp_46 - p_affine_2_1*tmp_40 + p_affine_2_1*tmp_45 + p_affine_2_2*tmp_39 - p_affine_2_2*tmp_42 - p_affine_3_0*tmp_43 + p_affine_3_0*tmp_46 + p_affine_3_1*tmp_40 - p_affine_3_1*tmp_45 - p_affine_3_2*tmp_39 + p_affine_3_2*tmp_42);
      real_t tmp_52 = 0.25*tmp_51;
      real_t tmp_53 = -tmp_29 - tmp_33 - tmp_37;
      real_t tmp_54 = -tmp_28 - tmp_32 - tmp_36;
      real_t tmp_55 = -tmp_27 - tmp_31 - tmp_35;
      real_t tmp_56 = 0.25*tmp_51;
      real_t tmp_57 = 0.25*tmp_51;
      real_t tmp_58 = 0.25*tmp_51;
      real_t a_0_0 = tmp_38*tmp_52;
      real_t a_0_1 = tmp_52*tmp_53;
      real_t a_0_2 = tmp_52*tmp_54;
      real_t a_0_3 = tmp_52*tmp_55;
      real_t a_1_0 = tmp_38*tmp_56;
      real_t a_1_1 = tmp_53*tmp_56;
      real_t a_1_2 = tmp_54*tmp_56;
      real_t a_1_3 = tmp_55*tmp_56;
      real_t a_2_0 = tmp_38*tmp_57;
      real_t a_2_1 = tmp_53*tmp_57;
      real_t a_2_2 = tmp_54*tmp_57;
      real_t a_2_3 = tmp_55*tmp_57;
      real_t a_3_0 = tmp_38*tmp_58;
      real_t a_3_1 = tmp_53*tmp_58;
      real_t a_3_2 = tmp_54*tmp_58;
      real_t a_3_3 = tmp_55*tmp_58;
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

   void p1_div_0_blending_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_14 = Blending_DF_Tetrahedron_blend_out0_id0*tmp_10 - Blending_DF_Tetrahedron_blend_out0_id0*tmp_12 + Blending_DF_Tetrahedron_blend_out1_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out1_id0*tmp_13 + Blending_DF_Tetrahedron_blend_out2_id0*tmp_8 - Blending_DF_Tetrahedron_blend_out2_id0*tmp_9;
      real_t tmp_15 = -p_affine_0_2;
      real_t tmp_16 = p_affine_3_2 + tmp_15;
      real_t tmp_17 = p_affine_1_2 + tmp_15;
      real_t tmp_18 = p_affine_3_1 + tmp_2;
      real_t tmp_19 = tmp_18*tmp_5;
      real_t tmp_20 = p_affine_2_2 + tmp_15;
      real_t tmp_21 = p_affine_3_0 + tmp_0;
      real_t tmp_22 = tmp_21*tmp_6;
      real_t tmp_23 = tmp_1*tmp_18;
      real_t tmp_24 = tmp_21*tmp_3;
      real_t tmp_25 = 1/(tmp_14*(tmp_16*tmp_4 - tmp_16*tmp_7 + tmp_17*tmp_19 - tmp_17*tmp_24 + tmp_20*tmp_22 - tmp_20*tmp_23));
      real_t tmp_26 = tmp_25*(tmp_8 - tmp_9);
      real_t tmp_27 = tmp_26*(tmp_4 - tmp_7);
      real_t tmp_28 = tmp_26*(tmp_22 - tmp_23);
      real_t tmp_29 = tmp_26*(tmp_19 - tmp_24);
      real_t tmp_30 = tmp_25*(tmp_11 - tmp_13);
      real_t tmp_31 = tmp_30*(-tmp_1*tmp_20 + tmp_17*tmp_5);
      real_t tmp_32 = tmp_30*(tmp_1*tmp_16 - tmp_17*tmp_21);
      real_t tmp_33 = tmp_30*(-tmp_16*tmp_5 + tmp_20*tmp_21);
      real_t tmp_34 = tmp_25*(tmp_10 - tmp_12);
      real_t tmp_35 = tmp_34*(-tmp_17*tmp_3 + tmp_20*tmp_6);
      real_t tmp_36 = tmp_34*(-tmp_16*tmp_6 + tmp_17*tmp_18);
      real_t tmp_37 = tmp_34*(tmp_16*tmp_3 - tmp_18*tmp_20);
      real_t tmp_38 = p_affine_0_0*p_affine_1_1;
      real_t tmp_39 = p_affine_0_0*p_affine_1_2;
      real_t tmp_40 = p_affine_2_1*p_affine_3_2;
      real_t tmp_41 = p_affine_0_1*p_affine_1_0;
      real_t tmp_42 = p_affine_0_1*p_affine_1_2;
      real_t tmp_43 = p_affine_2_2*p_affine_3_0;
      real_t tmp_44 = p_affine_0_2*p_affine_1_0;
      real_t tmp_45 = p_affine_0_2*p_affine_1_1;
      real_t tmp_46 = p_affine_2_0*p_affine_3_1;
      real_t tmp_47 = p_affine_2_2*p_affine_3_1;
      real_t tmp_48 = p_affine_2_0*p_affine_3_2;
      real_t tmp_49 = p_affine_2_1*p_affine_3_0;
      real_t tmp_50 = 0.041666666666666657*std::abs(tmp_14)*std::abs(p_affine_0_0*tmp_40 - p_affine_0_0*tmp_47 + p_affine_0_1*tmp_43 - p_affine_0_1*tmp_48 + p_affine_0_2*tmp_46 - p_affine_0_2*tmp_49 - p_affine_1_0*tmp_40 + p_affine_1_0*tmp_47 - p_affine_1_1*tmp_43 + p_affine_1_1*tmp_48 - p_affine_1_2*tmp_46 + p_affine_1_2*tmp_49 + p_affine_2_0*tmp_42 - p_affine_2_0*tmp_45 - p_affine_2_1*tmp_39 + p_affine_2_1*tmp_44 + p_affine_2_2*tmp_38 - p_affine_2_2*tmp_41 - p_affine_3_0*tmp_42 + p_affine_3_0*tmp_45 + p_affine_3_1*tmp_39 - p_affine_3_1*tmp_44 - p_affine_3_2*tmp_38 + p_affine_3_2*tmp_41);
      real_t a_0_0 = tmp_50*(tmp_27 + tmp_28 + tmp_29 + tmp_31 + tmp_32 + tmp_33 + tmp_35 + tmp_36 + tmp_37);
      real_t a_0_1 = tmp_50*(-tmp_29 - tmp_33 - tmp_37);
      real_t a_0_2 = tmp_50*(-tmp_28 - tmp_32 - tmp_36);
      real_t a_0_3 = tmp_50*(-tmp_27 - tmp_31 - tmp_35);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_div_0_blending_q1::Blending_DF_Triangle_blend( real_t in_0, real_t in_1, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3 ) const
   {
      Point3D  mappedPt( {in_0, in_1, 0} );
      Matrix2r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 1, 0 );
      *out_3 = DPsi( 1, 1 );
   }

   void p1_div_0_blending_q1::Blending_DF_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
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

   void p1_div_1_blending_q1::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
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
      Blending_DF_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_DF_Triangle_blend_out0_id0, &Blending_DF_Triangle_blend_out1_id0, &Blending_DF_Triangle_blend_out2_id0, &Blending_DF_Triangle_blend_out3_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = Blending_DF_Triangle_blend_out0_id0*Blending_DF_Triangle_blend_out3_id0 - Blending_DF_Triangle_blend_out1_id0*Blending_DF_Triangle_blend_out2_id0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = 1/(tmp_2*(tmp_1*tmp_4 - (p_affine_1_1 + tmp_3)*(p_affine_2_0 + tmp_0)));
      real_t tmp_6 = Blending_DF_Triangle_blend_out0_id0*tmp_5;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = tmp_6*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_9 = Blending_DF_Triangle_blend_out1_id0*tmp_5;
      real_t tmp_10 = tmp_4*tmp_9;
      real_t tmp_11 = tmp_9*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_12 = -tmp_10 - tmp_11 + tmp_7 + tmp_8;
      real_t tmp_13 = 0.5*std::abs(tmp_2)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_14 = 0.33333333333333343*tmp_13;
      real_t tmp_15 = tmp_10 - tmp_8;
      real_t tmp_16 = tmp_11 - tmp_7;
      real_t tmp_17 = 0.33333333333333331*tmp_13;
      real_t tmp_18 = 0.33333333333333331*tmp_13;
      real_t a_0_0 = tmp_12*tmp_14;
      real_t a_0_1 = tmp_14*tmp_15;
      real_t a_0_2 = tmp_14*tmp_16;
      real_t a_1_0 = tmp_12*tmp_17;
      real_t a_1_1 = tmp_15*tmp_17;
      real_t a_1_2 = tmp_16*tmp_17;
      real_t a_2_0 = tmp_12*tmp_18;
      real_t a_2_1 = tmp_15*tmp_18;
      real_t a_2_2 = tmp_16*tmp_18;
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

   void p1_div_1_blending_q1::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
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
      Blending_DF_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_DF_Triangle_blend_out0_id0, &Blending_DF_Triangle_blend_out1_id0, &Blending_DF_Triangle_blend_out2_id0, &Blending_DF_Triangle_blend_out3_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = Blending_DF_Triangle_blend_out0_id0*Blending_DF_Triangle_blend_out3_id0 - Blending_DF_Triangle_blend_out1_id0*Blending_DF_Triangle_blend_out2_id0;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = 1/(tmp_2*(tmp_1*tmp_4 - (p_affine_1_1 + tmp_3)*(p_affine_2_0 + tmp_0)));
      real_t tmp_6 = Blending_DF_Triangle_blend_out0_id0*tmp_5;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = tmp_6*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_9 = Blending_DF_Triangle_blend_out1_id0*tmp_5;
      real_t tmp_10 = tmp_4*tmp_9;
      real_t tmp_11 = tmp_9*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_12 = 0.16666666666666671*std::abs(tmp_2)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = tmp_12*(-tmp_10 - tmp_11 + tmp_7 + tmp_8);
      real_t a_0_1 = tmp_12*(tmp_10 - tmp_8);
      real_t a_0_2 = tmp_12*(tmp_11 - tmp_7);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_div_1_blending_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_7 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_8 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_12 = -Blending_DF_Tetrahedron_blend_out3_id0*tmp_10 + Blending_DF_Tetrahedron_blend_out3_id0*tmp_9 - Blending_DF_Tetrahedron_blend_out4_id0*tmp_11 + Blending_DF_Tetrahedron_blend_out4_id0*tmp_8 + Blending_DF_Tetrahedron_blend_out5_id0*tmp_6 - Blending_DF_Tetrahedron_blend_out5_id0*tmp_7;
      real_t tmp_13 = -p_affine_0_2;
      real_t tmp_14 = p_affine_3_2 + tmp_13;
      real_t tmp_15 = tmp_1*tmp_14;
      real_t tmp_16 = p_affine_3_1 + tmp_2;
      real_t tmp_17 = p_affine_1_2 + tmp_13;
      real_t tmp_18 = tmp_17*tmp_4;
      real_t tmp_19 = p_affine_3_0 + tmp_0;
      real_t tmp_20 = p_affine_2_2 + tmp_13;
      real_t tmp_21 = tmp_19*tmp_20;
      real_t tmp_22 = tmp_1*tmp_20;
      real_t tmp_23 = tmp_14*tmp_4;
      real_t tmp_24 = tmp_17*tmp_19;
      real_t tmp_25 = 1/(tmp_12*(tmp_15*tmp_3 + tmp_16*tmp_18 - tmp_16*tmp_22 + tmp_21*tmp_5 - tmp_23*tmp_5 - tmp_24*tmp_3));
      real_t tmp_26 = tmp_25*(tmp_6 - tmp_7);
      real_t tmp_27 = tmp_26*(tmp_1*tmp_3 - tmp_4*tmp_5);
      real_t tmp_28 = tmp_26*(-tmp_1*tmp_16 + tmp_19*tmp_5);
      real_t tmp_29 = tmp_26*(tmp_16*tmp_4 - tmp_19*tmp_3);
      real_t tmp_30 = tmp_25*(-tmp_11 + tmp_8);
      real_t tmp_31 = tmp_30*(tmp_18 - tmp_22);
      real_t tmp_32 = tmp_30*(tmp_15 - tmp_24);
      real_t tmp_33 = tmp_30*(tmp_21 - tmp_23);
      real_t tmp_34 = tmp_25*(-tmp_10 + tmp_9);
      real_t tmp_35 = tmp_34*(-tmp_17*tmp_3 + tmp_20*tmp_5);
      real_t tmp_36 = tmp_34*(-tmp_14*tmp_5 + tmp_16*tmp_17);
      real_t tmp_37 = tmp_34*(tmp_14*tmp_3 - tmp_16*tmp_20);
      real_t tmp_38 = tmp_27 + tmp_28 + tmp_29 + tmp_31 + tmp_32 + tmp_33 + tmp_35 + tmp_36 + tmp_37;
      real_t tmp_39 = p_affine_0_0*p_affine_1_1;
      real_t tmp_40 = p_affine_0_0*p_affine_1_2;
      real_t tmp_41 = p_affine_2_1*p_affine_3_2;
      real_t tmp_42 = p_affine_0_1*p_affine_1_0;
      real_t tmp_43 = p_affine_0_1*p_affine_1_2;
      real_t tmp_44 = p_affine_2_2*p_affine_3_0;
      real_t tmp_45 = p_affine_0_2*p_affine_1_0;
      real_t tmp_46 = p_affine_0_2*p_affine_1_1;
      real_t tmp_47 = p_affine_2_0*p_affine_3_1;
      real_t tmp_48 = p_affine_2_2*p_affine_3_1;
      real_t tmp_49 = p_affine_2_0*p_affine_3_2;
      real_t tmp_50 = p_affine_2_1*p_affine_3_0;
      real_t tmp_51 = 0.16666666666666663*std::abs(tmp_12)*std::abs(p_affine_0_0*tmp_41 - p_affine_0_0*tmp_48 + p_affine_0_1*tmp_44 - p_affine_0_1*tmp_49 + p_affine_0_2*tmp_47 - p_affine_0_2*tmp_50 - p_affine_1_0*tmp_41 + p_affine_1_0*tmp_48 - p_affine_1_1*tmp_44 + p_affine_1_1*tmp_49 - p_affine_1_2*tmp_47 + p_affine_1_2*tmp_50 + p_affine_2_0*tmp_43 - p_affine_2_0*tmp_46 - p_affine_2_1*tmp_40 + p_affine_2_1*tmp_45 + p_affine_2_2*tmp_39 - p_affine_2_2*tmp_42 - p_affine_3_0*tmp_43 + p_affine_3_0*tmp_46 + p_affine_3_1*tmp_40 - p_affine_3_1*tmp_45 - p_affine_3_2*tmp_39 + p_affine_3_2*tmp_42);
      real_t tmp_52 = 0.25*tmp_51;
      real_t tmp_53 = -tmp_29 - tmp_33 - tmp_37;
      real_t tmp_54 = -tmp_28 - tmp_32 - tmp_36;
      real_t tmp_55 = -tmp_27 - tmp_31 - tmp_35;
      real_t tmp_56 = 0.25*tmp_51;
      real_t tmp_57 = 0.25*tmp_51;
      real_t tmp_58 = 0.25*tmp_51;
      real_t a_0_0 = tmp_38*tmp_52;
      real_t a_0_1 = tmp_52*tmp_53;
      real_t a_0_2 = tmp_52*tmp_54;
      real_t a_0_3 = tmp_52*tmp_55;
      real_t a_1_0 = tmp_38*tmp_56;
      real_t a_1_1 = tmp_53*tmp_56;
      real_t a_1_2 = tmp_54*tmp_56;
      real_t a_1_3 = tmp_55*tmp_56;
      real_t a_2_0 = tmp_38*tmp_57;
      real_t a_2_1 = tmp_53*tmp_57;
      real_t a_2_2 = tmp_54*tmp_57;
      real_t a_2_3 = tmp_55*tmp_57;
      real_t a_3_0 = tmp_38*tmp_58;
      real_t a_3_1 = tmp_53*tmp_58;
      real_t a_3_2 = tmp_54*tmp_58;
      real_t a_3_3 = tmp_55*tmp_58;
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

   void p1_div_1_blending_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = p_affine_2_0 + tmp_0;
      real_t tmp_5 = p_affine_1_1 + tmp_2;
      real_t tmp_6 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_7 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_8 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out7_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out8_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out6_id0;
      real_t tmp_12 = -Blending_DF_Tetrahedron_blend_out3_id0*tmp_10 + Blending_DF_Tetrahedron_blend_out3_id0*tmp_9 - Blending_DF_Tetrahedron_blend_out4_id0*tmp_11 + Blending_DF_Tetrahedron_blend_out4_id0*tmp_8 + Blending_DF_Tetrahedron_blend_out5_id0*tmp_6 - Blending_DF_Tetrahedron_blend_out5_id0*tmp_7;
      real_t tmp_13 = -p_affine_0_2;
      real_t tmp_14 = p_affine_3_2 + tmp_13;
      real_t tmp_15 = tmp_1*tmp_14;
      real_t tmp_16 = p_affine_3_1 + tmp_2;
      real_t tmp_17 = p_affine_1_2 + tmp_13;
      real_t tmp_18 = tmp_17*tmp_4;
      real_t tmp_19 = p_affine_3_0 + tmp_0;
      real_t tmp_20 = p_affine_2_2 + tmp_13;
      real_t tmp_21 = tmp_19*tmp_20;
      real_t tmp_22 = tmp_1*tmp_20;
      real_t tmp_23 = tmp_14*tmp_4;
      real_t tmp_24 = tmp_17*tmp_19;
      real_t tmp_25 = 1/(tmp_12*(tmp_15*tmp_3 + tmp_16*tmp_18 - tmp_16*tmp_22 + tmp_21*tmp_5 - tmp_23*tmp_5 - tmp_24*tmp_3));
      real_t tmp_26 = tmp_25*(tmp_6 - tmp_7);
      real_t tmp_27 = tmp_26*(tmp_1*tmp_3 - tmp_4*tmp_5);
      real_t tmp_28 = tmp_26*(-tmp_1*tmp_16 + tmp_19*tmp_5);
      real_t tmp_29 = tmp_26*(tmp_16*tmp_4 - tmp_19*tmp_3);
      real_t tmp_30 = tmp_25*(-tmp_11 + tmp_8);
      real_t tmp_31 = tmp_30*(tmp_18 - tmp_22);
      real_t tmp_32 = tmp_30*(tmp_15 - tmp_24);
      real_t tmp_33 = tmp_30*(tmp_21 - tmp_23);
      real_t tmp_34 = tmp_25*(-tmp_10 + tmp_9);
      real_t tmp_35 = tmp_34*(-tmp_17*tmp_3 + tmp_20*tmp_5);
      real_t tmp_36 = tmp_34*(-tmp_14*tmp_5 + tmp_16*tmp_17);
      real_t tmp_37 = tmp_34*(tmp_14*tmp_3 - tmp_16*tmp_20);
      real_t tmp_38 = p_affine_0_0*p_affine_1_1;
      real_t tmp_39 = p_affine_0_0*p_affine_1_2;
      real_t tmp_40 = p_affine_2_1*p_affine_3_2;
      real_t tmp_41 = p_affine_0_1*p_affine_1_0;
      real_t tmp_42 = p_affine_0_1*p_affine_1_2;
      real_t tmp_43 = p_affine_2_2*p_affine_3_0;
      real_t tmp_44 = p_affine_0_2*p_affine_1_0;
      real_t tmp_45 = p_affine_0_2*p_affine_1_1;
      real_t tmp_46 = p_affine_2_0*p_affine_3_1;
      real_t tmp_47 = p_affine_2_2*p_affine_3_1;
      real_t tmp_48 = p_affine_2_0*p_affine_3_2;
      real_t tmp_49 = p_affine_2_1*p_affine_3_0;
      real_t tmp_50 = 0.041666666666666657*std::abs(tmp_12)*std::abs(p_affine_0_0*tmp_40 - p_affine_0_0*tmp_47 + p_affine_0_1*tmp_43 - p_affine_0_1*tmp_48 + p_affine_0_2*tmp_46 - p_affine_0_2*tmp_49 - p_affine_1_0*tmp_40 + p_affine_1_0*tmp_47 - p_affine_1_1*tmp_43 + p_affine_1_1*tmp_48 - p_affine_1_2*tmp_46 + p_affine_1_2*tmp_49 + p_affine_2_0*tmp_42 - p_affine_2_0*tmp_45 - p_affine_2_1*tmp_39 + p_affine_2_1*tmp_44 + p_affine_2_2*tmp_38 - p_affine_2_2*tmp_41 - p_affine_3_0*tmp_42 + p_affine_3_0*tmp_45 + p_affine_3_1*tmp_39 - p_affine_3_1*tmp_44 - p_affine_3_2*tmp_38 + p_affine_3_2*tmp_41);
      real_t a_0_0 = tmp_50*(tmp_27 + tmp_28 + tmp_29 + tmp_31 + tmp_32 + tmp_33 + tmp_35 + tmp_36 + tmp_37);
      real_t a_0_1 = tmp_50*(-tmp_29 - tmp_33 - tmp_37);
      real_t a_0_2 = tmp_50*(-tmp_28 - tmp_32 - tmp_36);
      real_t a_0_3 = tmp_50*(-tmp_27 - tmp_31 - tmp_35);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_div_1_blending_q1::Blending_DF_Triangle_blend( real_t in_0, real_t in_1, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3 ) const
   {
      Point3D  mappedPt( {in_0, in_1, 0} );
      Matrix2r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 1, 0 );
      *out_3 = DPsi( 1, 1 );
   }

   void p1_div_1_blending_q1::Blending_DF_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
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

   void p1_div_2_blending_q1::integrateAll( const std::array< Point3D, 3 >& , Matrix< real_t, 3, 3 >&  ) const
   {
      
   }

   void p1_div_2_blending_q1::integrateRow0( const std::array< Point3D, 3 >& , Matrix< real_t, 1, 3 >&  ) const
   {
      
   }

   void p1_div_2_blending_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_14 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_10 - Blending_DF_Tetrahedron_blend_out6_id0*tmp_13 + Blending_DF_Tetrahedron_blend_out7_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out7_id0*tmp_12 + Blending_DF_Tetrahedron_blend_out8_id0*tmp_8 - Blending_DF_Tetrahedron_blend_out8_id0*tmp_9;
      real_t tmp_15 = -p_affine_0_2;
      real_t tmp_16 = p_affine_3_2 + tmp_15;
      real_t tmp_17 = p_affine_1_2 + tmp_15;
      real_t tmp_18 = p_affine_3_1 + tmp_2;
      real_t tmp_19 = tmp_18*tmp_5;
      real_t tmp_20 = p_affine_2_2 + tmp_15;
      real_t tmp_21 = p_affine_3_0 + tmp_0;
      real_t tmp_22 = tmp_21*tmp_6;
      real_t tmp_23 = tmp_1*tmp_18;
      real_t tmp_24 = tmp_21*tmp_3;
      real_t tmp_25 = 1/(tmp_14*(tmp_16*tmp_4 - tmp_16*tmp_7 + tmp_17*tmp_19 - tmp_17*tmp_24 + tmp_20*tmp_22 - tmp_20*tmp_23));
      real_t tmp_26 = tmp_25*(tmp_8 - tmp_9);
      real_t tmp_27 = tmp_26*(tmp_4 - tmp_7);
      real_t tmp_28 = tmp_26*(tmp_22 - tmp_23);
      real_t tmp_29 = tmp_26*(tmp_19 - tmp_24);
      real_t tmp_30 = tmp_25*(tmp_11 - tmp_12);
      real_t tmp_31 = tmp_30*(-tmp_1*tmp_20 + tmp_17*tmp_5);
      real_t tmp_32 = tmp_30*(tmp_1*tmp_16 - tmp_17*tmp_21);
      real_t tmp_33 = tmp_30*(-tmp_16*tmp_5 + tmp_20*tmp_21);
      real_t tmp_34 = tmp_25*(tmp_10 - tmp_13);
      real_t tmp_35 = tmp_34*(-tmp_17*tmp_3 + tmp_20*tmp_6);
      real_t tmp_36 = tmp_34*(-tmp_16*tmp_6 + tmp_17*tmp_18);
      real_t tmp_37 = tmp_34*(tmp_16*tmp_3 - tmp_18*tmp_20);
      real_t tmp_38 = tmp_27 + tmp_28 + tmp_29 + tmp_31 + tmp_32 + tmp_33 + tmp_35 + tmp_36 + tmp_37;
      real_t tmp_39 = p_affine_0_0*p_affine_1_1;
      real_t tmp_40 = p_affine_0_0*p_affine_1_2;
      real_t tmp_41 = p_affine_2_1*p_affine_3_2;
      real_t tmp_42 = p_affine_0_1*p_affine_1_0;
      real_t tmp_43 = p_affine_0_1*p_affine_1_2;
      real_t tmp_44 = p_affine_2_2*p_affine_3_0;
      real_t tmp_45 = p_affine_0_2*p_affine_1_0;
      real_t tmp_46 = p_affine_0_2*p_affine_1_1;
      real_t tmp_47 = p_affine_2_0*p_affine_3_1;
      real_t tmp_48 = p_affine_2_2*p_affine_3_1;
      real_t tmp_49 = p_affine_2_0*p_affine_3_2;
      real_t tmp_50 = p_affine_2_1*p_affine_3_0;
      real_t tmp_51 = 0.16666666666666663*std::abs(tmp_14)*std::abs(p_affine_0_0*tmp_41 - p_affine_0_0*tmp_48 + p_affine_0_1*tmp_44 - p_affine_0_1*tmp_49 + p_affine_0_2*tmp_47 - p_affine_0_2*tmp_50 - p_affine_1_0*tmp_41 + p_affine_1_0*tmp_48 - p_affine_1_1*tmp_44 + p_affine_1_1*tmp_49 - p_affine_1_2*tmp_47 + p_affine_1_2*tmp_50 + p_affine_2_0*tmp_43 - p_affine_2_0*tmp_46 - p_affine_2_1*tmp_40 + p_affine_2_1*tmp_45 + p_affine_2_2*tmp_39 - p_affine_2_2*tmp_42 - p_affine_3_0*tmp_43 + p_affine_3_0*tmp_46 + p_affine_3_1*tmp_40 - p_affine_3_1*tmp_45 - p_affine_3_2*tmp_39 + p_affine_3_2*tmp_42);
      real_t tmp_52 = 0.25*tmp_51;
      real_t tmp_53 = -tmp_29 - tmp_33 - tmp_37;
      real_t tmp_54 = -tmp_28 - tmp_32 - tmp_36;
      real_t tmp_55 = -tmp_27 - tmp_31 - tmp_35;
      real_t tmp_56 = 0.25*tmp_51;
      real_t tmp_57 = 0.25*tmp_51;
      real_t tmp_58 = 0.25*tmp_51;
      real_t a_0_0 = tmp_38*tmp_52;
      real_t a_0_1 = tmp_52*tmp_53;
      real_t a_0_2 = tmp_52*tmp_54;
      real_t a_0_3 = tmp_52*tmp_55;
      real_t a_1_0 = tmp_38*tmp_56;
      real_t a_1_1 = tmp_53*tmp_56;
      real_t a_1_2 = tmp_54*tmp_56;
      real_t a_1_3 = tmp_55*tmp_56;
      real_t a_2_0 = tmp_38*tmp_57;
      real_t a_2_1 = tmp_53*tmp_57;
      real_t a_2_2 = tmp_54*tmp_57;
      real_t a_2_3 = tmp_55*tmp_57;
      real_t a_3_0 = tmp_38*tmp_58;
      real_t a_3_1 = tmp_53*tmp_58;
      real_t a_3_2 = tmp_54*tmp_58;
      real_t a_3_3 = tmp_55*tmp_58;
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

   void p1_div_2_blending_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_9 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_10 = Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_11 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0;
      real_t tmp_12 = Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0;
      real_t tmp_13 = Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0;
      real_t tmp_14 = Blending_DF_Tetrahedron_blend_out6_id0*tmp_10 - Blending_DF_Tetrahedron_blend_out6_id0*tmp_13 + Blending_DF_Tetrahedron_blend_out7_id0*tmp_11 - Blending_DF_Tetrahedron_blend_out7_id0*tmp_12 + Blending_DF_Tetrahedron_blend_out8_id0*tmp_8 - Blending_DF_Tetrahedron_blend_out8_id0*tmp_9;
      real_t tmp_15 = -p_affine_0_2;
      real_t tmp_16 = p_affine_3_2 + tmp_15;
      real_t tmp_17 = p_affine_1_2 + tmp_15;
      real_t tmp_18 = p_affine_3_1 + tmp_2;
      real_t tmp_19 = tmp_18*tmp_5;
      real_t tmp_20 = p_affine_2_2 + tmp_15;
      real_t tmp_21 = p_affine_3_0 + tmp_0;
      real_t tmp_22 = tmp_21*tmp_6;
      real_t tmp_23 = tmp_1*tmp_18;
      real_t tmp_24 = tmp_21*tmp_3;
      real_t tmp_25 = 1/(tmp_14*(tmp_16*tmp_4 - tmp_16*tmp_7 + tmp_17*tmp_19 - tmp_17*tmp_24 + tmp_20*tmp_22 - tmp_20*tmp_23));
      real_t tmp_26 = tmp_25*(tmp_8 - tmp_9);
      real_t tmp_27 = tmp_26*(tmp_4 - tmp_7);
      real_t tmp_28 = tmp_26*(tmp_22 - tmp_23);
      real_t tmp_29 = tmp_26*(tmp_19 - tmp_24);
      real_t tmp_30 = tmp_25*(tmp_11 - tmp_12);
      real_t tmp_31 = tmp_30*(-tmp_1*tmp_20 + tmp_17*tmp_5);
      real_t tmp_32 = tmp_30*(tmp_1*tmp_16 - tmp_17*tmp_21);
      real_t tmp_33 = tmp_30*(-tmp_16*tmp_5 + tmp_20*tmp_21);
      real_t tmp_34 = tmp_25*(tmp_10 - tmp_13);
      real_t tmp_35 = tmp_34*(-tmp_17*tmp_3 + tmp_20*tmp_6);
      real_t tmp_36 = tmp_34*(-tmp_16*tmp_6 + tmp_17*tmp_18);
      real_t tmp_37 = tmp_34*(tmp_16*tmp_3 - tmp_18*tmp_20);
      real_t tmp_38 = p_affine_0_0*p_affine_1_1;
      real_t tmp_39 = p_affine_0_0*p_affine_1_2;
      real_t tmp_40 = p_affine_2_1*p_affine_3_2;
      real_t tmp_41 = p_affine_0_1*p_affine_1_0;
      real_t tmp_42 = p_affine_0_1*p_affine_1_2;
      real_t tmp_43 = p_affine_2_2*p_affine_3_0;
      real_t tmp_44 = p_affine_0_2*p_affine_1_0;
      real_t tmp_45 = p_affine_0_2*p_affine_1_1;
      real_t tmp_46 = p_affine_2_0*p_affine_3_1;
      real_t tmp_47 = p_affine_2_2*p_affine_3_1;
      real_t tmp_48 = p_affine_2_0*p_affine_3_2;
      real_t tmp_49 = p_affine_2_1*p_affine_3_0;
      real_t tmp_50 = 0.041666666666666657*std::abs(tmp_14)*std::abs(p_affine_0_0*tmp_40 - p_affine_0_0*tmp_47 + p_affine_0_1*tmp_43 - p_affine_0_1*tmp_48 + p_affine_0_2*tmp_46 - p_affine_0_2*tmp_49 - p_affine_1_0*tmp_40 + p_affine_1_0*tmp_47 - p_affine_1_1*tmp_43 + p_affine_1_1*tmp_48 - p_affine_1_2*tmp_46 + p_affine_1_2*tmp_49 + p_affine_2_0*tmp_42 - p_affine_2_0*tmp_45 - p_affine_2_1*tmp_39 + p_affine_2_1*tmp_44 + p_affine_2_2*tmp_38 - p_affine_2_2*tmp_41 - p_affine_3_0*tmp_42 + p_affine_3_0*tmp_45 + p_affine_3_1*tmp_39 - p_affine_3_1*tmp_44 - p_affine_3_2*tmp_38 + p_affine_3_2*tmp_41);
      real_t a_0_0 = tmp_50*(tmp_27 + tmp_28 + tmp_29 + tmp_31 + tmp_32 + tmp_33 + tmp_35 + tmp_36 + tmp_37);
      real_t a_0_1 = tmp_50*(-tmp_29 - tmp_33 - tmp_37);
      real_t a_0_2 = tmp_50*(-tmp_28 - tmp_32 - tmp_36);
      real_t a_0_3 = tmp_50*(-tmp_27 - tmp_31 - tmp_35);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_div_2_blending_q1::Blending_DF_Triangle_blend( real_t in_0, real_t in_1, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3 ) const
   {
      Point3D  mappedPt( {in_0, in_1, 0} );
      Matrix2r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 1, 0 );
      *out_3 = DPsi( 1, 1 );
   }

   void p1_div_2_blending_q1::Blending_DF_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
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

} // namespace forms
} // namespace hyteg

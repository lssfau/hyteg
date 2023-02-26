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

#include "p1_div_k_grad_blending_q1.hpp"

namespace hyteg {
namespace forms {

   void p1_div_k_grad_blending_q1::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Blending_F_Triangle_0_0 = 0;
      real_t Blending_F_Triangle_0_1 = 0;
      real_t Blending_DF_Triangle_0_0 = 0;
      real_t Blending_DF_Triangle_0_1 = 0;
      real_t Blending_DF_Triangle_0_2 = 0;
      real_t Blending_DF_Triangle_0_3 = 0;
      real_t Scalar_Variable_Coefficient_2D_0_0 = 0;
      Blending_F_Triangle( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_F_Triangle_0_0, &Blending_F_Triangle_0_1 );
      Blending_DF_Triangle( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_DF_Triangle_0_0, &Blending_DF_Triangle_0_1, &Blending_DF_Triangle_0_2, &Blending_DF_Triangle_0_3 );
      Scalar_Variable_Coefficient_2D( Blending_F_Triangle_0_0, Blending_F_Triangle_0_1, &Scalar_Variable_Coefficient_2D_0_0 );
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = Blending_DF_Triangle_0_0*Blending_DF_Triangle_0_3 - Blending_DF_Triangle_0_1*Blending_DF_Triangle_0_2;
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = p_affine_1_0 + tmp_3;
      real_t tmp_5 = 1/(tmp_2*(tmp_1*tmp_4 - (p_affine_1_1 + tmp_0)*(p_affine_2_0 + tmp_3)));
      real_t tmp_6 = Blending_DF_Triangle_0_1*tmp_5;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_9 = tmp_6*tmp_8;
      real_t tmp_10 = Blending_DF_Triangle_0_0*tmp_5;
      real_t tmp_11 = tmp_10*tmp_4;
      real_t tmp_12 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_13 = tmp_10*tmp_12;
      real_t tmp_14 = -tmp_11 - tmp_13 + tmp_7 + tmp_9;
      real_t tmp_15 = Scalar_Variable_Coefficient_2D_0_0*tmp_7;
      real_t tmp_16 = Scalar_Variable_Coefficient_2D_0_0*tmp_9;
      real_t tmp_17 = Scalar_Variable_Coefficient_2D_0_0*tmp_11;
      real_t tmp_18 = Scalar_Variable_Coefficient_2D_0_0*tmp_13;
      real_t tmp_19 = tmp_15 + tmp_16 - tmp_17 - tmp_18;
      real_t tmp_20 = Blending_DF_Triangle_0_2*tmp_5;
      real_t tmp_21 = tmp_20*tmp_4;
      real_t tmp_22 = tmp_12*tmp_20;
      real_t tmp_23 = Blending_DF_Triangle_0_3*tmp_5;
      real_t tmp_24 = tmp_1*tmp_23;
      real_t tmp_25 = tmp_23*tmp_8;
      real_t tmp_26 = tmp_21 + tmp_22 - tmp_24 - tmp_25;
      real_t tmp_27 = Scalar_Variable_Coefficient_2D_0_0*tmp_21;
      real_t tmp_28 = Scalar_Variable_Coefficient_2D_0_0*tmp_22;
      real_t tmp_29 = Scalar_Variable_Coefficient_2D_0_0*tmp_24;
      real_t tmp_30 = Scalar_Variable_Coefficient_2D_0_0*tmp_25;
      real_t tmp_31 = tmp_27 + tmp_28 - tmp_29 - tmp_30;
      real_t tmp_32 = 0.5*std::abs(tmp_2)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_33 = -tmp_15 + tmp_18;
      real_t tmp_34 = -tmp_28 + tmp_29;
      real_t tmp_35 = -tmp_16 + tmp_17;
      real_t tmp_36 = -tmp_27 + tmp_30;
      real_t tmp_37 = tmp_13 - tmp_7;
      real_t tmp_38 = -tmp_22 + tmp_24;
      real_t tmp_39 = tmp_11 - tmp_9;
      real_t tmp_40 = -tmp_21 + tmp_25;
      real_t a_0_0 = tmp_32*(tmp_14*tmp_19 + tmp_26*tmp_31);
      real_t a_0_1 = tmp_32*(tmp_14*tmp_33 + tmp_26*tmp_34);
      real_t a_0_2 = tmp_32*(tmp_14*tmp_35 + tmp_26*tmp_36);
      real_t a_1_0 = tmp_32*(tmp_19*tmp_37 + tmp_31*tmp_38);
      real_t a_1_1 = tmp_32*(tmp_33*tmp_37 + tmp_34*tmp_38);
      real_t a_1_2 = tmp_32*(tmp_35*tmp_37 + tmp_36*tmp_38);
      real_t a_2_0 = tmp_32*(tmp_19*tmp_39 + tmp_31*tmp_40);
      real_t a_2_1 = tmp_32*(tmp_33*tmp_39 + tmp_34*tmp_40);
      real_t a_2_2 = tmp_32*(tmp_35*tmp_39 + tmp_36*tmp_40);
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

   void p1_div_k_grad_blending_q1::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Blending_F_Triangle_0_0 = 0;
      real_t Blending_F_Triangle_0_1 = 0;
      real_t Blending_DF_Triangle_0_0 = 0;
      real_t Blending_DF_Triangle_0_1 = 0;
      real_t Blending_DF_Triangle_0_2 = 0;
      real_t Blending_DF_Triangle_0_3 = 0;
      real_t Scalar_Variable_Coefficient_2D_0_0 = 0;
      Blending_F_Triangle( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_F_Triangle_0_0, &Blending_F_Triangle_0_1 );
      Blending_DF_Triangle( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_DF_Triangle_0_0, &Blending_DF_Triangle_0_1, &Blending_DF_Triangle_0_2, &Blending_DF_Triangle_0_3 );
      Scalar_Variable_Coefficient_2D( Blending_F_Triangle_0_0, Blending_F_Triangle_0_1, &Scalar_Variable_Coefficient_2D_0_0 );
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = Blending_DF_Triangle_0_0*Blending_DF_Triangle_0_3 - Blending_DF_Triangle_0_1*Blending_DF_Triangle_0_2;
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = p_affine_1_0 + tmp_3;
      real_t tmp_5 = 1/(tmp_2*(tmp_1*tmp_4 - (p_affine_1_1 + tmp_0)*(p_affine_2_0 + tmp_3)));
      real_t tmp_6 = Blending_DF_Triangle_0_1*tmp_5;
      real_t tmp_7 = tmp_1*tmp_6;
      real_t tmp_8 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_9 = tmp_6*tmp_8;
      real_t tmp_10 = Blending_DF_Triangle_0_0*tmp_5;
      real_t tmp_11 = tmp_10*tmp_4;
      real_t tmp_12 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_13 = tmp_10*tmp_12;
      real_t tmp_14 = -tmp_11 - tmp_13 + tmp_7 + tmp_9;
      real_t tmp_15 = Scalar_Variable_Coefficient_2D_0_0*tmp_7;
      real_t tmp_16 = Scalar_Variable_Coefficient_2D_0_0*tmp_9;
      real_t tmp_17 = Scalar_Variable_Coefficient_2D_0_0*tmp_11;
      real_t tmp_18 = Scalar_Variable_Coefficient_2D_0_0*tmp_13;
      real_t tmp_19 = Blending_DF_Triangle_0_2*tmp_5;
      real_t tmp_20 = tmp_19*tmp_4;
      real_t tmp_21 = tmp_12*tmp_19;
      real_t tmp_22 = Blending_DF_Triangle_0_3*tmp_5;
      real_t tmp_23 = tmp_1*tmp_22;
      real_t tmp_24 = tmp_22*tmp_8;
      real_t tmp_25 = tmp_20 + tmp_21 - tmp_23 - tmp_24;
      real_t tmp_26 = Scalar_Variable_Coefficient_2D_0_0*tmp_20;
      real_t tmp_27 = Scalar_Variable_Coefficient_2D_0_0*tmp_21;
      real_t tmp_28 = Scalar_Variable_Coefficient_2D_0_0*tmp_23;
      real_t tmp_29 = Scalar_Variable_Coefficient_2D_0_0*tmp_24;
      real_t tmp_30 = 0.5*std::abs(tmp_2)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = tmp_30*(tmp_14*(tmp_15 + tmp_16 - tmp_17 - tmp_18) + tmp_25*(tmp_26 + tmp_27 - tmp_28 - tmp_29));
      real_t a_0_1 = tmp_30*(tmp_14*(-tmp_15 + tmp_18) + tmp_25*(-tmp_27 + tmp_28));
      real_t a_0_2 = tmp_30*(tmp_14*(-tmp_16 + tmp_17) + tmp_25*(-tmp_26 + tmp_29));
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_div_k_grad_blending_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Blending_F_Tetrahedron_0_0 = 0;
      real_t Blending_F_Tetrahedron_0_1 = 0;
      real_t Blending_F_Tetrahedron_0_2 = 0;
      real_t Blending_DF_Tetrahedron_0_0 = 0;
      real_t Blending_DF_Tetrahedron_0_1 = 0;
      real_t Blending_DF_Tetrahedron_0_2 = 0;
      real_t Blending_DF_Tetrahedron_0_3 = 0;
      real_t Blending_DF_Tetrahedron_0_4 = 0;
      real_t Blending_DF_Tetrahedron_0_5 = 0;
      real_t Blending_DF_Tetrahedron_0_6 = 0;
      real_t Blending_DF_Tetrahedron_0_7 = 0;
      real_t Blending_DF_Tetrahedron_0_8 = 0;
      real_t Scalar_Variable_Coefficient_3D_0_0 = 0;
      Blending_F_Tetrahedron( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_0_0, &Blending_F_Tetrahedron_0_1, &Blending_F_Tetrahedron_0_2 );
      Blending_DF_Tetrahedron( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_0_0, &Blending_DF_Tetrahedron_0_1, &Blending_DF_Tetrahedron_0_2, &Blending_DF_Tetrahedron_0_3, &Blending_DF_Tetrahedron_0_4, &Blending_DF_Tetrahedron_0_5, &Blending_DF_Tetrahedron_0_6, &Blending_DF_Tetrahedron_0_7, &Blending_DF_Tetrahedron_0_8 );
      Scalar_Variable_Coefficient_3D( Blending_F_Tetrahedron_0_0, Blending_F_Tetrahedron_0_1, Blending_F_Tetrahedron_0_2, &Scalar_Variable_Coefficient_3D_0_0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = Blending_DF_Tetrahedron_0_0*Blending_DF_Tetrahedron_0_4;
      real_t tmp_10 = Blending_DF_Tetrahedron_0_1*Blending_DF_Tetrahedron_0_3;
      real_t tmp_11 = Blending_DF_Tetrahedron_0_1*Blending_DF_Tetrahedron_0_5;
      real_t tmp_12 = Blending_DF_Tetrahedron_0_2*Blending_DF_Tetrahedron_0_3;
      real_t tmp_13 = Blending_DF_Tetrahedron_0_0*Blending_DF_Tetrahedron_0_5;
      real_t tmp_14 = Blending_DF_Tetrahedron_0_2*Blending_DF_Tetrahedron_0_4;
      real_t tmp_15 = Blending_DF_Tetrahedron_0_6*tmp_11 - Blending_DF_Tetrahedron_0_6*tmp_14 + Blending_DF_Tetrahedron_0_7*tmp_12 - Blending_DF_Tetrahedron_0_7*tmp_13 - Blending_DF_Tetrahedron_0_8*tmp_10 + Blending_DF_Tetrahedron_0_8*tmp_9;
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
      real_t tmp_48 = Scalar_Variable_Coefficient_3D_0_0*tmp_28;
      real_t tmp_49 = Scalar_Variable_Coefficient_3D_0_0*tmp_30;
      real_t tmp_50 = Scalar_Variable_Coefficient_3D_0_0*tmp_32;
      real_t tmp_51 = Scalar_Variable_Coefficient_3D_0_0*tmp_35;
      real_t tmp_52 = Scalar_Variable_Coefficient_3D_0_0*tmp_37;
      real_t tmp_53 = Scalar_Variable_Coefficient_3D_0_0*tmp_39;
      real_t tmp_54 = Scalar_Variable_Coefficient_3D_0_0*tmp_42;
      real_t tmp_55 = Scalar_Variable_Coefficient_3D_0_0*tmp_44;
      real_t tmp_56 = Scalar_Variable_Coefficient_3D_0_0*tmp_46;
      real_t tmp_57 = -tmp_48 - tmp_49 - tmp_50 - tmp_51 - tmp_52 - tmp_53 - tmp_54 - tmp_55 - tmp_56;
      real_t tmp_58 = tmp_26*(-Blending_DF_Tetrahedron_0_0*Blending_DF_Tetrahedron_0_7 + Blending_DF_Tetrahedron_0_1*Blending_DF_Tetrahedron_0_6);
      real_t tmp_59 = tmp_58*tmp_8;
      real_t tmp_60 = tmp_29*tmp_58;
      real_t tmp_61 = tmp_31*tmp_58;
      real_t tmp_62 = tmp_26*(Blending_DF_Tetrahedron_0_0*Blending_DF_Tetrahedron_0_8 - Blending_DF_Tetrahedron_0_2*Blending_DF_Tetrahedron_0_6);
      real_t tmp_63 = tmp_33*tmp_62;
      real_t tmp_64 = tmp_36*tmp_62;
      real_t tmp_65 = tmp_38*tmp_62;
      real_t tmp_66 = tmp_26*(-Blending_DF_Tetrahedron_0_1*Blending_DF_Tetrahedron_0_8 + Blending_DF_Tetrahedron_0_2*Blending_DF_Tetrahedron_0_7);
      real_t tmp_67 = tmp_40*tmp_66;
      real_t tmp_68 = tmp_43*tmp_66;
      real_t tmp_69 = tmp_45*tmp_66;
      real_t tmp_70 = -tmp_59 - tmp_60 - tmp_61 - tmp_63 - tmp_64 - tmp_65 - tmp_67 - tmp_68 - tmp_69;
      real_t tmp_71 = Scalar_Variable_Coefficient_3D_0_0*tmp_59;
      real_t tmp_72 = Scalar_Variable_Coefficient_3D_0_0*tmp_60;
      real_t tmp_73 = Scalar_Variable_Coefficient_3D_0_0*tmp_61;
      real_t tmp_74 = Scalar_Variable_Coefficient_3D_0_0*tmp_63;
      real_t tmp_75 = Scalar_Variable_Coefficient_3D_0_0*tmp_64;
      real_t tmp_76 = Scalar_Variable_Coefficient_3D_0_0*tmp_65;
      real_t tmp_77 = Scalar_Variable_Coefficient_3D_0_0*tmp_67;
      real_t tmp_78 = Scalar_Variable_Coefficient_3D_0_0*tmp_68;
      real_t tmp_79 = Scalar_Variable_Coefficient_3D_0_0*tmp_69;
      real_t tmp_80 = -tmp_71 - tmp_72 - tmp_73 - tmp_74 - tmp_75 - tmp_76 - tmp_77 - tmp_78 - tmp_79;
      real_t tmp_81 = tmp_26*(Blending_DF_Tetrahedron_0_3*Blending_DF_Tetrahedron_0_7 - Blending_DF_Tetrahedron_0_4*Blending_DF_Tetrahedron_0_6);
      real_t tmp_82 = tmp_8*tmp_81;
      real_t tmp_83 = tmp_29*tmp_81;
      real_t tmp_84 = tmp_31*tmp_81;
      real_t tmp_85 = tmp_26*(-Blending_DF_Tetrahedron_0_3*Blending_DF_Tetrahedron_0_8 + Blending_DF_Tetrahedron_0_5*Blending_DF_Tetrahedron_0_6);
      real_t tmp_86 = tmp_33*tmp_85;
      real_t tmp_87 = tmp_36*tmp_85;
      real_t tmp_88 = tmp_38*tmp_85;
      real_t tmp_89 = tmp_26*(Blending_DF_Tetrahedron_0_4*Blending_DF_Tetrahedron_0_8 - Blending_DF_Tetrahedron_0_5*Blending_DF_Tetrahedron_0_7);
      real_t tmp_90 = tmp_40*tmp_89;
      real_t tmp_91 = tmp_43*tmp_89;
      real_t tmp_92 = tmp_45*tmp_89;
      real_t tmp_93 = -tmp_82 - tmp_83 - tmp_84 - tmp_86 - tmp_87 - tmp_88 - tmp_90 - tmp_91 - tmp_92;
      real_t tmp_94 = Scalar_Variable_Coefficient_3D_0_0*tmp_82;
      real_t tmp_95 = Scalar_Variable_Coefficient_3D_0_0*tmp_83;
      real_t tmp_96 = Scalar_Variable_Coefficient_3D_0_0*tmp_84;
      real_t tmp_97 = Scalar_Variable_Coefficient_3D_0_0*tmp_86;
      real_t tmp_98 = Scalar_Variable_Coefficient_3D_0_0*tmp_87;
      real_t tmp_99 = Scalar_Variable_Coefficient_3D_0_0*tmp_88;
      real_t tmp_100 = Scalar_Variable_Coefficient_3D_0_0*tmp_90;
      real_t tmp_101 = Scalar_Variable_Coefficient_3D_0_0*tmp_91;
      real_t tmp_102 = Scalar_Variable_Coefficient_3D_0_0*tmp_92;
      real_t tmp_103 = -tmp_100 - tmp_101 - tmp_102 - tmp_94 - tmp_95 - tmp_96 - tmp_97 - tmp_98 - tmp_99;
      real_t tmp_104 = p_affine_0_0*p_affine_1_1;
      real_t tmp_105 = p_affine_0_0*p_affine_1_2;
      real_t tmp_106 = p_affine_2_1*p_affine_3_2;
      real_t tmp_107 = p_affine_0_1*p_affine_1_0;
      real_t tmp_108 = p_affine_0_1*p_affine_1_2;
      real_t tmp_109 = p_affine_2_2*p_affine_3_0;
      real_t tmp_110 = p_affine_0_2*p_affine_1_0;
      real_t tmp_111 = p_affine_0_2*p_affine_1_1;
      real_t tmp_112 = p_affine_2_0*p_affine_3_1;
      real_t tmp_113 = p_affine_2_2*p_affine_3_1;
      real_t tmp_114 = p_affine_2_0*p_affine_3_2;
      real_t tmp_115 = p_affine_2_1*p_affine_3_0;
      real_t tmp_116 = 0.16666666666666663*std::abs(tmp_15)*std::abs(p_affine_0_0*tmp_106 - p_affine_0_0*tmp_113 + p_affine_0_1*tmp_109 - p_affine_0_1*tmp_114 + p_affine_0_2*tmp_112 - p_affine_0_2*tmp_115 - p_affine_1_0*tmp_106 + p_affine_1_0*tmp_113 - p_affine_1_1*tmp_109 + p_affine_1_1*tmp_114 - p_affine_1_2*tmp_112 + p_affine_1_2*tmp_115 + p_affine_2_0*tmp_108 - p_affine_2_0*tmp_111 - p_affine_2_1*tmp_105 + p_affine_2_1*tmp_110 + p_affine_2_2*tmp_104 - p_affine_2_2*tmp_107 - p_affine_3_0*tmp_108 + p_affine_3_0*tmp_111 + p_affine_3_1*tmp_105 - p_affine_3_1*tmp_110 - p_affine_3_2*tmp_104 + p_affine_3_2*tmp_107);
      real_t tmp_117 = tmp_50 + tmp_53 + tmp_56;
      real_t tmp_118 = tmp_73 + tmp_76 + tmp_79;
      real_t tmp_119 = tmp_102 + tmp_96 + tmp_99;
      real_t tmp_120 = tmp_49 + tmp_52 + tmp_55;
      real_t tmp_121 = tmp_72 + tmp_75 + tmp_78;
      real_t tmp_122 = tmp_101 + tmp_95 + tmp_98;
      real_t tmp_123 = tmp_48 + tmp_51 + tmp_54;
      real_t tmp_124 = tmp_71 + tmp_74 + tmp_77;
      real_t tmp_125 = tmp_100 + tmp_94 + tmp_97;
      real_t tmp_126 = tmp_32 + tmp_39 + tmp_46;
      real_t tmp_127 = tmp_61 + tmp_65 + tmp_69;
      real_t tmp_128 = tmp_84 + tmp_88 + tmp_92;
      real_t tmp_129 = tmp_30 + tmp_37 + tmp_44;
      real_t tmp_130 = tmp_60 + tmp_64 + tmp_68;
      real_t tmp_131 = tmp_83 + tmp_87 + tmp_91;
      real_t tmp_132 = tmp_28 + tmp_35 + tmp_42;
      real_t tmp_133 = tmp_59 + tmp_63 + tmp_67;
      real_t tmp_134 = tmp_82 + tmp_86 + tmp_90;
      real_t a_0_0 = tmp_116*(tmp_103*tmp_93 + tmp_47*tmp_57 + tmp_70*tmp_80);
      real_t a_0_1 = tmp_116*(tmp_117*tmp_47 + tmp_118*tmp_70 + tmp_119*tmp_93);
      real_t a_0_2 = tmp_116*(tmp_120*tmp_47 + tmp_121*tmp_70 + tmp_122*tmp_93);
      real_t a_0_3 = tmp_116*(tmp_123*tmp_47 + tmp_124*tmp_70 + tmp_125*tmp_93);
      real_t a_1_0 = tmp_116*(tmp_103*tmp_128 + tmp_126*tmp_57 + tmp_127*tmp_80);
      real_t a_1_1 = tmp_116*(tmp_117*tmp_126 + tmp_118*tmp_127 + tmp_119*tmp_128);
      real_t a_1_2 = tmp_116*(tmp_120*tmp_126 + tmp_121*tmp_127 + tmp_122*tmp_128);
      real_t a_1_3 = tmp_116*(tmp_123*tmp_126 + tmp_124*tmp_127 + tmp_125*tmp_128);
      real_t a_2_0 = tmp_116*(tmp_103*tmp_131 + tmp_129*tmp_57 + tmp_130*tmp_80);
      real_t a_2_1 = tmp_116*(tmp_117*tmp_129 + tmp_118*tmp_130 + tmp_119*tmp_131);
      real_t a_2_2 = tmp_116*(tmp_120*tmp_129 + tmp_121*tmp_130 + tmp_122*tmp_131);
      real_t a_2_3 = tmp_116*(tmp_123*tmp_129 + tmp_124*tmp_130 + tmp_125*tmp_131);
      real_t a_3_0 = tmp_116*(tmp_103*tmp_134 + tmp_132*tmp_57 + tmp_133*tmp_80);
      real_t a_3_1 = tmp_116*(tmp_117*tmp_132 + tmp_118*tmp_133 + tmp_119*tmp_134);
      real_t a_3_2 = tmp_116*(tmp_120*tmp_132 + tmp_121*tmp_133 + tmp_122*tmp_134);
      real_t a_3_3 = tmp_116*(tmp_123*tmp_132 + tmp_124*tmp_133 + tmp_125*tmp_134);
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

   void p1_div_k_grad_blending_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Blending_F_Tetrahedron_0_0 = 0;
      real_t Blending_F_Tetrahedron_0_1 = 0;
      real_t Blending_F_Tetrahedron_0_2 = 0;
      real_t Blending_DF_Tetrahedron_0_0 = 0;
      real_t Blending_DF_Tetrahedron_0_1 = 0;
      real_t Blending_DF_Tetrahedron_0_2 = 0;
      real_t Blending_DF_Tetrahedron_0_3 = 0;
      real_t Blending_DF_Tetrahedron_0_4 = 0;
      real_t Blending_DF_Tetrahedron_0_5 = 0;
      real_t Blending_DF_Tetrahedron_0_6 = 0;
      real_t Blending_DF_Tetrahedron_0_7 = 0;
      real_t Blending_DF_Tetrahedron_0_8 = 0;
      real_t Scalar_Variable_Coefficient_3D_0_0 = 0;
      Blending_F_Tetrahedron( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_0_0, &Blending_F_Tetrahedron_0_1, &Blending_F_Tetrahedron_0_2 );
      Blending_DF_Tetrahedron( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_0_0, &Blending_DF_Tetrahedron_0_1, &Blending_DF_Tetrahedron_0_2, &Blending_DF_Tetrahedron_0_3, &Blending_DF_Tetrahedron_0_4, &Blending_DF_Tetrahedron_0_5, &Blending_DF_Tetrahedron_0_6, &Blending_DF_Tetrahedron_0_7, &Blending_DF_Tetrahedron_0_8 );
      Scalar_Variable_Coefficient_3D( Blending_F_Tetrahedron_0_0, Blending_F_Tetrahedron_0_1, Blending_F_Tetrahedron_0_2, &Scalar_Variable_Coefficient_3D_0_0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_0 + tmp_0;
      real_t tmp_6 = p_affine_1_1 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = tmp_4 - tmp_7;
      real_t tmp_9 = Blending_DF_Tetrahedron_0_0*Blending_DF_Tetrahedron_0_4;
      real_t tmp_10 = Blending_DF_Tetrahedron_0_1*Blending_DF_Tetrahedron_0_3;
      real_t tmp_11 = Blending_DF_Tetrahedron_0_1*Blending_DF_Tetrahedron_0_5;
      real_t tmp_12 = Blending_DF_Tetrahedron_0_2*Blending_DF_Tetrahedron_0_3;
      real_t tmp_13 = Blending_DF_Tetrahedron_0_0*Blending_DF_Tetrahedron_0_5;
      real_t tmp_14 = Blending_DF_Tetrahedron_0_2*Blending_DF_Tetrahedron_0_4;
      real_t tmp_15 = Blending_DF_Tetrahedron_0_6*tmp_11 - Blending_DF_Tetrahedron_0_6*tmp_14 + Blending_DF_Tetrahedron_0_7*tmp_12 - Blending_DF_Tetrahedron_0_7*tmp_13 - Blending_DF_Tetrahedron_0_8*tmp_10 + Blending_DF_Tetrahedron_0_8*tmp_9;
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
      real_t tmp_48 = Scalar_Variable_Coefficient_3D_0_0*tmp_28;
      real_t tmp_49 = Scalar_Variable_Coefficient_3D_0_0*tmp_30;
      real_t tmp_50 = Scalar_Variable_Coefficient_3D_0_0*tmp_32;
      real_t tmp_51 = Scalar_Variable_Coefficient_3D_0_0*tmp_35;
      real_t tmp_52 = Scalar_Variable_Coefficient_3D_0_0*tmp_37;
      real_t tmp_53 = Scalar_Variable_Coefficient_3D_0_0*tmp_39;
      real_t tmp_54 = Scalar_Variable_Coefficient_3D_0_0*tmp_42;
      real_t tmp_55 = Scalar_Variable_Coefficient_3D_0_0*tmp_44;
      real_t tmp_56 = Scalar_Variable_Coefficient_3D_0_0*tmp_46;
      real_t tmp_57 = tmp_26*(-Blending_DF_Tetrahedron_0_0*Blending_DF_Tetrahedron_0_7 + Blending_DF_Tetrahedron_0_1*Blending_DF_Tetrahedron_0_6);
      real_t tmp_58 = tmp_57*tmp_8;
      real_t tmp_59 = tmp_29*tmp_57;
      real_t tmp_60 = tmp_31*tmp_57;
      real_t tmp_61 = tmp_26*(Blending_DF_Tetrahedron_0_0*Blending_DF_Tetrahedron_0_8 - Blending_DF_Tetrahedron_0_2*Blending_DF_Tetrahedron_0_6);
      real_t tmp_62 = tmp_33*tmp_61;
      real_t tmp_63 = tmp_36*tmp_61;
      real_t tmp_64 = tmp_38*tmp_61;
      real_t tmp_65 = tmp_26*(-Blending_DF_Tetrahedron_0_1*Blending_DF_Tetrahedron_0_8 + Blending_DF_Tetrahedron_0_2*Blending_DF_Tetrahedron_0_7);
      real_t tmp_66 = tmp_40*tmp_65;
      real_t tmp_67 = tmp_43*tmp_65;
      real_t tmp_68 = tmp_45*tmp_65;
      real_t tmp_69 = -tmp_58 - tmp_59 - tmp_60 - tmp_62 - tmp_63 - tmp_64 - tmp_66 - tmp_67 - tmp_68;
      real_t tmp_70 = Scalar_Variable_Coefficient_3D_0_0*tmp_58;
      real_t tmp_71 = Scalar_Variable_Coefficient_3D_0_0*tmp_59;
      real_t tmp_72 = Scalar_Variable_Coefficient_3D_0_0*tmp_60;
      real_t tmp_73 = Scalar_Variable_Coefficient_3D_0_0*tmp_62;
      real_t tmp_74 = Scalar_Variable_Coefficient_3D_0_0*tmp_63;
      real_t tmp_75 = Scalar_Variable_Coefficient_3D_0_0*tmp_64;
      real_t tmp_76 = Scalar_Variable_Coefficient_3D_0_0*tmp_66;
      real_t tmp_77 = Scalar_Variable_Coefficient_3D_0_0*tmp_67;
      real_t tmp_78 = Scalar_Variable_Coefficient_3D_0_0*tmp_68;
      real_t tmp_79 = tmp_26*(Blending_DF_Tetrahedron_0_3*Blending_DF_Tetrahedron_0_7 - Blending_DF_Tetrahedron_0_4*Blending_DF_Tetrahedron_0_6);
      real_t tmp_80 = tmp_79*tmp_8;
      real_t tmp_81 = tmp_29*tmp_79;
      real_t tmp_82 = tmp_31*tmp_79;
      real_t tmp_83 = tmp_26*(-Blending_DF_Tetrahedron_0_3*Blending_DF_Tetrahedron_0_8 + Blending_DF_Tetrahedron_0_5*Blending_DF_Tetrahedron_0_6);
      real_t tmp_84 = tmp_33*tmp_83;
      real_t tmp_85 = tmp_36*tmp_83;
      real_t tmp_86 = tmp_38*tmp_83;
      real_t tmp_87 = tmp_26*(Blending_DF_Tetrahedron_0_4*Blending_DF_Tetrahedron_0_8 - Blending_DF_Tetrahedron_0_5*Blending_DF_Tetrahedron_0_7);
      real_t tmp_88 = tmp_40*tmp_87;
      real_t tmp_89 = tmp_43*tmp_87;
      real_t tmp_90 = tmp_45*tmp_87;
      real_t tmp_91 = -tmp_80 - tmp_81 - tmp_82 - tmp_84 - tmp_85 - tmp_86 - tmp_88 - tmp_89 - tmp_90;
      real_t tmp_92 = Scalar_Variable_Coefficient_3D_0_0*tmp_80;
      real_t tmp_93 = Scalar_Variable_Coefficient_3D_0_0*tmp_81;
      real_t tmp_94 = Scalar_Variable_Coefficient_3D_0_0*tmp_82;
      real_t tmp_95 = Scalar_Variable_Coefficient_3D_0_0*tmp_84;
      real_t tmp_96 = Scalar_Variable_Coefficient_3D_0_0*tmp_85;
      real_t tmp_97 = Scalar_Variable_Coefficient_3D_0_0*tmp_86;
      real_t tmp_98 = Scalar_Variable_Coefficient_3D_0_0*tmp_88;
      real_t tmp_99 = Scalar_Variable_Coefficient_3D_0_0*tmp_89;
      real_t tmp_100 = Scalar_Variable_Coefficient_3D_0_0*tmp_90;
      real_t tmp_101 = p_affine_0_0*p_affine_1_1;
      real_t tmp_102 = p_affine_0_0*p_affine_1_2;
      real_t tmp_103 = p_affine_2_1*p_affine_3_2;
      real_t tmp_104 = p_affine_0_1*p_affine_1_0;
      real_t tmp_105 = p_affine_0_1*p_affine_1_2;
      real_t tmp_106 = p_affine_2_2*p_affine_3_0;
      real_t tmp_107 = p_affine_0_2*p_affine_1_0;
      real_t tmp_108 = p_affine_0_2*p_affine_1_1;
      real_t tmp_109 = p_affine_2_0*p_affine_3_1;
      real_t tmp_110 = p_affine_2_2*p_affine_3_1;
      real_t tmp_111 = p_affine_2_0*p_affine_3_2;
      real_t tmp_112 = p_affine_2_1*p_affine_3_0;
      real_t tmp_113 = 0.16666666666666663*std::abs(tmp_15)*std::abs(p_affine_0_0*tmp_103 - p_affine_0_0*tmp_110 + p_affine_0_1*tmp_106 - p_affine_0_1*tmp_111 + p_affine_0_2*tmp_109 - p_affine_0_2*tmp_112 - p_affine_1_0*tmp_103 + p_affine_1_0*tmp_110 - p_affine_1_1*tmp_106 + p_affine_1_1*tmp_111 - p_affine_1_2*tmp_109 + p_affine_1_2*tmp_112 + p_affine_2_0*tmp_105 - p_affine_2_0*tmp_108 - p_affine_2_1*tmp_102 + p_affine_2_1*tmp_107 + p_affine_2_2*tmp_101 - p_affine_2_2*tmp_104 - p_affine_3_0*tmp_105 + p_affine_3_0*tmp_108 + p_affine_3_1*tmp_102 - p_affine_3_1*tmp_107 - p_affine_3_2*tmp_101 + p_affine_3_2*tmp_104);
      real_t a_0_0 = tmp_113*(tmp_47*(-tmp_48 - tmp_49 - tmp_50 - tmp_51 - tmp_52 - tmp_53 - tmp_54 - tmp_55 - tmp_56) + tmp_69*(-tmp_70 - tmp_71 - tmp_72 - tmp_73 - tmp_74 - tmp_75 - tmp_76 - tmp_77 - tmp_78) + tmp_91*(-tmp_100 - tmp_92 - tmp_93 - tmp_94 - tmp_95 - tmp_96 - tmp_97 - tmp_98 - tmp_99));
      real_t a_0_1 = tmp_113*(tmp_47*(tmp_50 + tmp_53 + tmp_56) + tmp_69*(tmp_72 + tmp_75 + tmp_78) + tmp_91*(tmp_100 + tmp_94 + tmp_97));
      real_t a_0_2 = tmp_113*(tmp_47*(tmp_49 + tmp_52 + tmp_55) + tmp_69*(tmp_71 + tmp_74 + tmp_77) + tmp_91*(tmp_93 + tmp_96 + tmp_99));
      real_t a_0_3 = tmp_113*(tmp_47*(tmp_48 + tmp_51 + tmp_54) + tmp_69*(tmp_70 + tmp_73 + tmp_76) + tmp_91*(tmp_92 + tmp_95 + tmp_98));
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_div_k_grad_blending_q1::Blending_F_Triangle( real_t in_0, real_t in_1, real_t * out_0, real_t * out_1 ) const
   {
      Point3D  in( in_0, in_1, 0 );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
   }

   void p1_div_k_grad_blending_q1::Blending_DF_Triangle( real_t in_0, real_t in_1, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3 ) const
   {
      Point3D  mappedPt( in_0, in_1, 0 );
      Matrix2r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 1, 0 );
      *out_3 = DPsi( 1, 1 );
   }

   void p1_div_k_grad_blending_q1::Scalar_Variable_Coefficient_2D( real_t in_0, real_t in_1, real_t * out_0 ) const
   {
      *out_0 = callback2D( Point3D( in_0, in_1, 0 ) );
   }

   void p1_div_k_grad_blending_q1::Blending_F_Tetrahedron( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2 ) const
   {
      Point3D  in( in_0, in_1, in_2 );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
      *out_2 = out[2];
   }

   void p1_div_k_grad_blending_q1::Blending_DF_Tetrahedron( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
   {
      Point3D  mappedPt( in_0, in_1, in_2 );
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

   void p1_div_k_grad_blending_q1::Scalar_Variable_Coefficient_3D( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback3D( Point3D( in_0, in_1, in_2 ) );
   }

} // namespace forms
} // namespace hyteg

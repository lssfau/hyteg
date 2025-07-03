/*
* Copyright (c) 2025 Marcus Mohr.
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
* The computational kernels in thie file were generated with the HyTeG Operator Generator.
* Interfaces needed slight modifications, as HOG support for DG1/P2+ pair is not complete, yet.
*/

#include "dg1_to_p2_plus_bubble_divt_affine_q3.hpp"

namespace hyteg {
namespace forms {

   void dg1_to_p2_plus_bubble_divt_0_affine_q3::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 7, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_2 = -0.66666666666666674;
      real_t tmp_3 = -tmp_2 - 0.33333333333333331;
      real_t tmp_4 = 1.3333333333333333;
      real_t tmp_5 = 1.3333333333333333;
      real_t tmp_6 = tmp_4 + tmp_5 - 3;
      real_t tmp_9 = -0.40000000000000002;
      real_t tmp_10 = -tmp_9 - 0.20000000000000001;
      real_t tmp_11 = 0.80000000000000004;
      real_t tmp_12 = 2.3999999999999999;
      real_t tmp_13 = tmp_11 + tmp_12 - 3;
      real_t tmp_16 = -0.80000000000000004;
      real_t tmp_17 = -tmp_16 - 0.59999999999999998;
      real_t tmp_18 = 2.3999999999999999;
      real_t tmp_19 = 0.80000000000000004;
      real_t tmp_20 = tmp_18 + tmp_19 - 3;
      real_t tmp_23 = -0.80000000000000004;
      real_t tmp_24 = -tmp_23 - 0.20000000000000001;
      real_t tmp_25 = 0.80000000000000004;
      real_t tmp_26 = 0.80000000000000004;
      real_t tmp_27 = tmp_25 + tmp_26 - 3;
      real_t jac_affine_0_0 = -p_affine_0_0 + p_affine_1_0;
      real_t jac_affine_0_1 = -p_affine_0_0 + p_affine_2_0;
      real_t jac_affine_1_0 = -p_affine_0_1 + p_affine_1_1;
      real_t jac_affine_1_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_0 = jac_affine_0_0*jac_affine_1_1 - jac_affine_0_1*jac_affine_1_0;
      real_t tmp_1 = 1.0 / (tmp_0);
      real_t jac_affine_inv_0_0 = jac_affine_1_1*tmp_1;
      real_t tmp_31 = jac_affine_inv_0_0*(tmp_4 - 1);
      real_t tmp_33 = jac_affine_inv_0_0*(tmp_11 - 1);
      real_t tmp_35 = jac_affine_inv_0_0*(tmp_18 - 1);
      real_t tmp_37 = jac_affine_inv_0_0*(tmp_25 - 1);
      real_t tmp_50 = jac_affine_inv_0_0*tmp_5;
      real_t tmp_53 = jac_affine_inv_0_0*tmp_12;
      real_t tmp_56 = jac_affine_inv_0_0*tmp_19;
      real_t tmp_59 = jac_affine_inv_0_0*tmp_26;
      real_t tmp_70 = 27*jac_affine_inv_0_0;
      real_t jac_affine_inv_1_0 = -jac_affine_1_0*tmp_1;
      real_t tmp_42 = jac_affine_inv_1_0*(tmp_5 - 1);
      real_t tmp_43 = jac_affine_inv_1_0*(tmp_12 - 1);
      real_t tmp_44 = jac_affine_inv_1_0*(tmp_19 - 1);
      real_t tmp_45 = jac_affine_inv_1_0*(tmp_26 - 1);
      real_t tmp_51 = jac_affine_inv_1_0*tmp_4;
      real_t tmp_54 = jac_affine_inv_1_0*tmp_11;
      real_t tmp_57 = jac_affine_inv_1_0*tmp_18;
      real_t tmp_60 = jac_affine_inv_1_0*tmp_25;
      real_t tmp_71 = 27*jac_affine_inv_1_0;
      real_t abs_det_jac_affine = std::abs(tmp_0);
      real_t tmp_7 = -0.28125*abs_det_jac_affine;
      real_t tmp_8 = tmp_7*(jac_affine_inv_0_0*tmp_6 + jac_affine_inv_1_0*tmp_6);
      real_t tmp_14 = 0.26041666666666669*abs_det_jac_affine;
      real_t tmp_15 = tmp_14*(jac_affine_inv_0_0*tmp_13 + jac_affine_inv_1_0*tmp_13);
      real_t tmp_21 = 0.26041666666666669*abs_det_jac_affine;
      real_t tmp_22 = tmp_21*(jac_affine_inv_0_0*tmp_20 + jac_affine_inv_1_0*tmp_20);
      real_t tmp_28 = 0.26041666666666669*abs_det_jac_affine;
      real_t tmp_29 = tmp_28*(jac_affine_inv_0_0*tmp_27 + jac_affine_inv_1_0*tmp_27);
      real_t tmp_30 = tmp_3*tmp_7;
      real_t tmp_32 = tmp_10*tmp_14;
      real_t tmp_34 = tmp_17*tmp_21;
      real_t tmp_36 = tmp_24*tmp_28;
      real_t tmp_38 = tmp_31*tmp_7;
      real_t tmp_39 = tmp_14*tmp_33;
      real_t tmp_40 = tmp_21*tmp_35;
      real_t tmp_41 = tmp_28*tmp_37;
      real_t tmp_46 = tmp_42*tmp_7;
      real_t tmp_47 = tmp_14*tmp_43;
      real_t tmp_48 = tmp_21*tmp_44;
      real_t tmp_49 = tmp_28*tmp_45;
      real_t tmp_52 = tmp_7*(tmp_50 + tmp_51);
      real_t tmp_55 = tmp_14*(tmp_53 + tmp_54);
      real_t tmp_58 = tmp_21*(tmp_56 + tmp_57);
      real_t tmp_61 = tmp_28*(tmp_59 + tmp_60);
      real_t tmp_62 = tmp_7*(jac_affine_inv_1_0*(1.3333333333333335 - tmp_4) - tmp_50);
      real_t tmp_63 = tmp_14*(jac_affine_inv_1_0*(-tmp_11 - 0.79999999999999982) - tmp_53);
      real_t tmp_64 = tmp_21*(jac_affine_inv_1_0*(2.3999999999999999 - tmp_18) - tmp_56);
      real_t tmp_65 = tmp_28*(jac_affine_inv_1_0*(2.3999999999999999 - tmp_25) - tmp_59);
      real_t tmp_66 = tmp_7*(jac_affine_inv_0_0*(1.3333333333333335 - tmp_5) - tmp_51);
      real_t tmp_67 = tmp_14*(jac_affine_inv_0_0*(2.3999999999999999 - tmp_12) - tmp_54);
      real_t tmp_68 = tmp_21*(jac_affine_inv_0_0*(-tmp_19 - 0.79999999999999982) - tmp_57);
      real_t tmp_69 = tmp_28*(jac_affine_inv_0_0*(2.3999999999999999 - tmp_26) - tmp_60);
      real_t tmp_72 = tmp_7*(0.33333333333333331*tmp_70*(-tmp_2 - 0.66666666666666663) + 3.7007434154171883e-17*tmp_71);
      real_t tmp_73 = tmp_14*(0.59999999999999998*tmp_70*(-tmp_9 - 0.40000000000000002) - 0.079999999999999988*tmp_71);
      real_t tmp_74 = 0.20000000000000001*tmp_21*tmp_70*(-tmp_16 - 1.2);
      real_t tmp_75 = tmp_28*(0.20000000000000001*tmp_70*(-tmp_23 - 0.40000000000000002) + 0.080000000000000016*tmp_71);
      real_t a_0_0 = -tmp_10*tmp_15 - tmp_17*tmp_22 - tmp_24*tmp_29 - tmp_3*tmp_8;
      real_t a_0_1 = -0.20000000000000001*tmp_15 - 0.59999999999999998*tmp_22 - 0.20000000000000001*tmp_29 - 0.33333333333333331*tmp_8;
      real_t a_0_2 = -0.59999999999999998*tmp_15 - 0.20000000000000001*tmp_22 - 0.20000000000000001*tmp_29 - 0.33333333333333331*tmp_8;
      real_t a_1_0 = -tmp_30*tmp_31 - tmp_32*tmp_33 - tmp_34*tmp_35 - tmp_36*tmp_37;
      real_t a_1_1 = -0.33333333333333331*tmp_38 - 0.20000000000000001*tmp_39 - 0.59999999999999998*tmp_40 - 0.20000000000000001*tmp_41;
      real_t a_1_2 = -0.33333333333333331*tmp_38 - 0.59999999999999998*tmp_39 - 0.20000000000000001*tmp_40 - 0.20000000000000001*tmp_41;
      real_t a_2_0 = -tmp_30*tmp_42 - tmp_32*tmp_43 - tmp_34*tmp_44 - tmp_36*tmp_45;
      real_t a_2_1 = -0.33333333333333331*tmp_46 - 0.20000000000000001*tmp_47 - 0.59999999999999998*tmp_48 - 0.20000000000000001*tmp_49;
      real_t a_2_2 = -0.33333333333333331*tmp_46 - 0.59999999999999998*tmp_47 - 0.20000000000000001*tmp_48 - 0.20000000000000001*tmp_49;
      real_t a_3_0 = -tmp_10*tmp_55 - tmp_17*tmp_58 - tmp_24*tmp_61 - tmp_3*tmp_52;
      real_t a_3_1 = -0.33333333333333331*tmp_52 - 0.20000000000000001*tmp_55 - 0.59999999999999998*tmp_58 - 0.20000000000000001*tmp_61;
      real_t a_3_2 = -0.33333333333333331*tmp_52 - 0.59999999999999998*tmp_55 - 0.20000000000000001*tmp_58 - 0.20000000000000001*tmp_61;
      real_t a_4_0 = -tmp_10*tmp_63 - tmp_17*tmp_64 - tmp_24*tmp_65 - tmp_3*tmp_62;
      real_t a_4_1 = -0.33333333333333331*tmp_62 - 0.20000000000000001*tmp_63 - 0.59999999999999998*tmp_64 - 0.20000000000000001*tmp_65;
      real_t a_4_2 = -0.33333333333333331*tmp_62 - 0.59999999999999998*tmp_63 - 0.20000000000000001*tmp_64 - 0.20000000000000001*tmp_65;
      real_t a_5_0 = -tmp_10*tmp_67 - tmp_17*tmp_68 - tmp_24*tmp_69 - tmp_3*tmp_66;
      real_t a_5_1 = -0.33333333333333331*tmp_66 - 0.20000000000000001*tmp_67 - 0.59999999999999998*tmp_68 - 0.20000000000000001*tmp_69;
      real_t a_5_2 = -0.33333333333333331*tmp_66 - 0.59999999999999998*tmp_67 - 0.20000000000000001*tmp_68 - 0.20000000000000001*tmp_69;
      real_t a_6_0 = -tmp_10*tmp_73 - tmp_17*tmp_74 - tmp_24*tmp_75 - tmp_3*tmp_72;
      real_t a_6_1 = -0.33333333333333331*tmp_72 - 0.20000000000000001*tmp_73 - 0.59999999999999998*tmp_74 - 0.20000000000000001*tmp_75;
      real_t a_6_2 = -0.33333333333333331*tmp_72 - 0.59999999999999998*tmp_73 - 0.20000000000000001*tmp_74 - 0.20000000000000001*tmp_75;
      elMat(0,0) = a_0_0;
      elMat(0,1) = a_0_1;
      elMat(0,2) = a_0_2;
      elMat(1,0) = a_1_0;
      elMat(1,1) = a_1_1;
      elMat(1,2) = a_1_2;
      elMat(2,0) = a_2_0;
      elMat(2,1) = a_2_1;
      elMat(2,2) = a_2_2;
      elMat(3,0) = a_3_0;
      elMat(3,1) = a_3_1;
      elMat(3,2) = a_3_2;
      elMat(4,0) = a_4_0;
      elMat(4,1) = a_4_1;
      elMat(4,2) = a_4_2;
      elMat(5,0) = a_5_0;
      elMat(5,1) = a_5_1;
      elMat(5,2) = a_5_2;
      elMat(6,0) = a_6_0;
      elMat(6,1) = a_6_1;
      elMat(6,2) = a_6_2;
   }

   void dg1_to_p2_plus_bubble_divt_0_affine_q3::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_2 = -0.33333333333333348;
      real_t tmp_4 = 0.19999999999999973;
      real_t tmp_6 = 0.19999999999999996;
      real_t tmp_8 = -1.4000000000000001;
      real_t jac_affine_0_0 = -p_affine_0_0 + p_affine_1_0;
      real_t jac_affine_0_1 = -p_affine_0_0 + p_affine_2_0;
      real_t jac_affine_1_0 = -p_affine_0_1 + p_affine_1_1;
      real_t jac_affine_1_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_0 = jac_affine_0_0*jac_affine_1_1 - jac_affine_0_1*jac_affine_1_0;
      real_t tmp_1 = 1.0 / (tmp_0);
      real_t jac_affine_inv_0_0 = jac_affine_1_1*tmp_1;
      real_t jac_affine_inv_1_0 = -jac_affine_1_0*tmp_1;
      real_t abs_det_jac_affine = std::abs(tmp_0);
      real_t tmp_3 = -0.28125*abs_det_jac_affine*(jac_affine_inv_0_0*tmp_2 + jac_affine_inv_1_0*tmp_2);
      real_t tmp_5 = 0.26041666666666669*abs_det_jac_affine*(jac_affine_inv_0_0*tmp_4 + jac_affine_inv_1_0*tmp_4);
      real_t tmp_7 = 0.26041666666666669*abs_det_jac_affine*(jac_affine_inv_0_0*tmp_6 + jac_affine_inv_1_0*tmp_6);
      real_t tmp_9 = 0.26041666666666669*abs_det_jac_affine*(jac_affine_inv_0_0*tmp_8 + jac_affine_inv_1_0*tmp_8);
      real_t a_0_0 = -0.33333333333333343*tmp_3 - 0.20000000000000007*tmp_5 - 0.20000000000000001*tmp_7 - 0.60000000000000009*tmp_9;
      real_t a_0_1 = -0.33333333333333331*tmp_3 - 0.20000000000000001*tmp_5 - 0.59999999999999998*tmp_7 - 0.20000000000000001*tmp_9;
      real_t a_0_2 = -0.33333333333333331*tmp_3 - 0.59999999999999998*tmp_5 - 0.20000000000000001*tmp_7 - 0.20000000000000001*tmp_9;
      elMat(0,0) = a_0_0;
      elMat(0,1) = a_0_1;
      elMat(0,2) = a_0_2;
   }

   void dg1_to_p2_plus_bubble_divt_1_affine_q3::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 7, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_2 = -0.66666666666666674;
      real_t tmp_3 = -tmp_2 - 0.33333333333333331;
      real_t tmp_4 = 1.3333333333333333;
      real_t tmp_5 = 1.3333333333333333;
      real_t tmp_6 = tmp_4 + tmp_5 - 3;
      real_t tmp_9 = -0.40000000000000002;
      real_t tmp_10 = -tmp_9 - 0.20000000000000001;
      real_t tmp_11 = 0.80000000000000004;
      real_t tmp_12 = 2.3999999999999999;
      real_t tmp_13 = tmp_11 + tmp_12 - 3;
      real_t tmp_16 = -0.80000000000000004;
      real_t tmp_17 = -tmp_16 - 0.59999999999999998;
      real_t tmp_18 = 2.3999999999999999;
      real_t tmp_19 = 0.80000000000000004;
      real_t tmp_20 = tmp_18 + tmp_19 - 3;
      real_t tmp_23 = -0.80000000000000004;
      real_t tmp_24 = -tmp_23 - 0.20000000000000001;
      real_t tmp_25 = 0.80000000000000004;
      real_t tmp_26 = 0.80000000000000004;
      real_t tmp_27 = tmp_25 + tmp_26 - 3;
      real_t jac_affine_0_0 = -p_affine_0_0 + p_affine_1_0;
      real_t jac_affine_0_1 = -p_affine_0_0 + p_affine_2_0;
      real_t jac_affine_1_0 = -p_affine_0_1 + p_affine_1_1;
      real_t jac_affine_1_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_0 = jac_affine_0_0*jac_affine_1_1 - jac_affine_0_1*jac_affine_1_0;
      real_t tmp_1 = 1.0 / (tmp_0);
      real_t jac_affine_inv_0_1 = -jac_affine_0_1*tmp_1;
      real_t tmp_31 = jac_affine_inv_0_1*(tmp_4 - 1);
      real_t tmp_33 = jac_affine_inv_0_1*(tmp_11 - 1);
      real_t tmp_35 = jac_affine_inv_0_1*(tmp_18 - 1);
      real_t tmp_37 = jac_affine_inv_0_1*(tmp_25 - 1);
      real_t tmp_50 = jac_affine_inv_0_1*tmp_5;
      real_t tmp_53 = jac_affine_inv_0_1*tmp_12;
      real_t tmp_56 = jac_affine_inv_0_1*tmp_19;
      real_t tmp_59 = jac_affine_inv_0_1*tmp_26;
      real_t tmp_70 = 27*jac_affine_inv_0_1;
      real_t jac_affine_inv_1_1 = jac_affine_0_0*tmp_1;
      real_t tmp_42 = jac_affine_inv_1_1*(tmp_5 - 1);
      real_t tmp_43 = jac_affine_inv_1_1*(tmp_12 - 1);
      real_t tmp_44 = jac_affine_inv_1_1*(tmp_19 - 1);
      real_t tmp_45 = jac_affine_inv_1_1*(tmp_26 - 1);
      real_t tmp_51 = jac_affine_inv_1_1*tmp_4;
      real_t tmp_54 = jac_affine_inv_1_1*tmp_11;
      real_t tmp_57 = jac_affine_inv_1_1*tmp_18;
      real_t tmp_60 = jac_affine_inv_1_1*tmp_25;
      real_t tmp_71 = 27*jac_affine_inv_1_1;
      real_t abs_det_jac_affine = std::abs(tmp_0);
      real_t tmp_7 = -0.28125*abs_det_jac_affine;
      real_t tmp_8 = tmp_7*(jac_affine_inv_0_1*tmp_6 + jac_affine_inv_1_1*tmp_6);
      real_t tmp_14 = 0.26041666666666669*abs_det_jac_affine;
      real_t tmp_15 = tmp_14*(jac_affine_inv_0_1*tmp_13 + jac_affine_inv_1_1*tmp_13);
      real_t tmp_21 = 0.26041666666666669*abs_det_jac_affine;
      real_t tmp_22 = tmp_21*(jac_affine_inv_0_1*tmp_20 + jac_affine_inv_1_1*tmp_20);
      real_t tmp_28 = 0.26041666666666669*abs_det_jac_affine;
      real_t tmp_29 = tmp_28*(jac_affine_inv_0_1*tmp_27 + jac_affine_inv_1_1*tmp_27);
      real_t tmp_30 = tmp_3*tmp_7;
      real_t tmp_32 = tmp_10*tmp_14;
      real_t tmp_34 = tmp_17*tmp_21;
      real_t tmp_36 = tmp_24*tmp_28;
      real_t tmp_38 = tmp_31*tmp_7;
      real_t tmp_39 = tmp_14*tmp_33;
      real_t tmp_40 = tmp_21*tmp_35;
      real_t tmp_41 = tmp_28*tmp_37;
      real_t tmp_46 = tmp_42*tmp_7;
      real_t tmp_47 = tmp_14*tmp_43;
      real_t tmp_48 = tmp_21*tmp_44;
      real_t tmp_49 = tmp_28*tmp_45;
      real_t tmp_52 = tmp_7*(tmp_50 + tmp_51);
      real_t tmp_55 = tmp_14*(tmp_53 + tmp_54);
      real_t tmp_58 = tmp_21*(tmp_56 + tmp_57);
      real_t tmp_61 = tmp_28*(tmp_59 + tmp_60);
      real_t tmp_62 = tmp_7*(jac_affine_inv_1_1*(1.3333333333333335 - tmp_4) - tmp_50);
      real_t tmp_63 = tmp_14*(jac_affine_inv_1_1*(-tmp_11 - 0.79999999999999982) - tmp_53);
      real_t tmp_64 = tmp_21*(jac_affine_inv_1_1*(2.3999999999999999 - tmp_18) - tmp_56);
      real_t tmp_65 = tmp_28*(jac_affine_inv_1_1*(2.3999999999999999 - tmp_25) - tmp_59);
      real_t tmp_66 = tmp_7*(jac_affine_inv_0_1*(1.3333333333333335 - tmp_5) - tmp_51);
      real_t tmp_67 = tmp_14*(jac_affine_inv_0_1*(2.3999999999999999 - tmp_12) - tmp_54);
      real_t tmp_68 = tmp_21*(jac_affine_inv_0_1*(-tmp_19 - 0.79999999999999982) - tmp_57);
      real_t tmp_69 = tmp_28*(jac_affine_inv_0_1*(2.3999999999999999 - tmp_26) - tmp_60);
      real_t tmp_72 = tmp_7*(0.33333333333333331*tmp_70*(-tmp_2 - 0.66666666666666663) + 3.7007434154171883e-17*tmp_71);
      real_t tmp_73 = tmp_14*(0.59999999999999998*tmp_70*(-tmp_9 - 0.40000000000000002) - 0.079999999999999988*tmp_71);
      real_t tmp_74 = 0.20000000000000001*tmp_21*tmp_70*(-tmp_16 - 1.2);
      real_t tmp_75 = tmp_28*(0.20000000000000001*tmp_70*(-tmp_23 - 0.40000000000000002) + 0.080000000000000016*tmp_71);
      real_t a_0_0 = -tmp_10*tmp_15 - tmp_17*tmp_22 - tmp_24*tmp_29 - tmp_3*tmp_8;
      real_t a_0_1 = -0.20000000000000001*tmp_15 - 0.59999999999999998*tmp_22 - 0.20000000000000001*tmp_29 - 0.33333333333333331*tmp_8;
      real_t a_0_2 = -0.59999999999999998*tmp_15 - 0.20000000000000001*tmp_22 - 0.20000000000000001*tmp_29 - 0.33333333333333331*tmp_8;
      real_t a_1_0 = -tmp_30*tmp_31 - tmp_32*tmp_33 - tmp_34*tmp_35 - tmp_36*tmp_37;
      real_t a_1_1 = -0.33333333333333331*tmp_38 - 0.20000000000000001*tmp_39 - 0.59999999999999998*tmp_40 - 0.20000000000000001*tmp_41;
      real_t a_1_2 = -0.33333333333333331*tmp_38 - 0.59999999999999998*tmp_39 - 0.20000000000000001*tmp_40 - 0.20000000000000001*tmp_41;
      real_t a_2_0 = -tmp_30*tmp_42 - tmp_32*tmp_43 - tmp_34*tmp_44 - tmp_36*tmp_45;
      real_t a_2_1 = -0.33333333333333331*tmp_46 - 0.20000000000000001*tmp_47 - 0.59999999999999998*tmp_48 - 0.20000000000000001*tmp_49;
      real_t a_2_2 = -0.33333333333333331*tmp_46 - 0.59999999999999998*tmp_47 - 0.20000000000000001*tmp_48 - 0.20000000000000001*tmp_49;
      real_t a_3_0 = -tmp_10*tmp_55 - tmp_17*tmp_58 - tmp_24*tmp_61 - tmp_3*tmp_52;
      real_t a_3_1 = -0.33333333333333331*tmp_52 - 0.20000000000000001*tmp_55 - 0.59999999999999998*tmp_58 - 0.20000000000000001*tmp_61;
      real_t a_3_2 = -0.33333333333333331*tmp_52 - 0.59999999999999998*tmp_55 - 0.20000000000000001*tmp_58 - 0.20000000000000001*tmp_61;
      real_t a_4_0 = -tmp_10*tmp_63 - tmp_17*tmp_64 - tmp_24*tmp_65 - tmp_3*tmp_62;
      real_t a_4_1 = -0.33333333333333331*tmp_62 - 0.20000000000000001*tmp_63 - 0.59999999999999998*tmp_64 - 0.20000000000000001*tmp_65;
      real_t a_4_2 = -0.33333333333333331*tmp_62 - 0.59999999999999998*tmp_63 - 0.20000000000000001*tmp_64 - 0.20000000000000001*tmp_65;
      real_t a_5_0 = -tmp_10*tmp_67 - tmp_17*tmp_68 - tmp_24*tmp_69 - tmp_3*tmp_66;
      real_t a_5_1 = -0.33333333333333331*tmp_66 - 0.20000000000000001*tmp_67 - 0.59999999999999998*tmp_68 - 0.20000000000000001*tmp_69;
      real_t a_5_2 = -0.33333333333333331*tmp_66 - 0.59999999999999998*tmp_67 - 0.20000000000000001*tmp_68 - 0.20000000000000001*tmp_69;
      real_t a_6_0 = -tmp_10*tmp_73 - tmp_17*tmp_74 - tmp_24*tmp_75 - tmp_3*tmp_72;
      real_t a_6_1 = -0.33333333333333331*tmp_72 - 0.20000000000000001*tmp_73 - 0.59999999999999998*tmp_74 - 0.20000000000000001*tmp_75;
      real_t a_6_2 = -0.33333333333333331*tmp_72 - 0.59999999999999998*tmp_73 - 0.20000000000000001*tmp_74 - 0.20000000000000001*tmp_75;
      elMat(0,0) = a_0_0;
      elMat(0,1) = a_0_1;
      elMat(0,2) = a_0_2;
      elMat(1,0) = a_1_0;
      elMat(1,1) = a_1_1;
      elMat(1,2) = a_1_2;
      elMat(2,0) = a_2_0;
      elMat(2,1) = a_2_1;
      elMat(2,2) = a_2_2;
      elMat(3,0) = a_3_0;
      elMat(3,1) = a_3_1;
      elMat(3,2) = a_3_2;
      elMat(4,0) = a_4_0;
      elMat(4,1) = a_4_1;
      elMat(4,2) = a_4_2;
      elMat(5,0) = a_5_0;
      elMat(5,1) = a_5_1;
      elMat(5,2) = a_5_2;
      elMat(6,0) = a_6_0;
      elMat(6,1) = a_6_1;
      elMat(6,2) = a_6_2;
   }

   void dg1_to_p2_plus_bubble_divt_1_affine_q3::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_2 = -0.33333333333333348;
      real_t tmp_4 = 0.19999999999999973;
      real_t tmp_6 = 0.19999999999999996;
      real_t tmp_8 = -1.4000000000000001;
      real_t jac_affine_0_0 = -p_affine_0_0 + p_affine_1_0;
      real_t jac_affine_0_1 = -p_affine_0_0 + p_affine_2_0;
      real_t jac_affine_1_0 = -p_affine_0_1 + p_affine_1_1;
      real_t jac_affine_1_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_0 = jac_affine_0_0*jac_affine_1_1 - jac_affine_0_1*jac_affine_1_0;
      real_t tmp_1 = 1.0 / (tmp_0);
      real_t jac_affine_inv_0_1 = -jac_affine_0_1*tmp_1;
      real_t jac_affine_inv_1_1 = jac_affine_0_0*tmp_1;
      real_t abs_det_jac_affine = std::abs(tmp_0);
      real_t tmp_3 = -0.28125*abs_det_jac_affine*(jac_affine_inv_0_1*tmp_2 + jac_affine_inv_1_1*tmp_2);
      real_t tmp_5 = 0.26041666666666669*abs_det_jac_affine*(jac_affine_inv_0_1*tmp_4 + jac_affine_inv_1_1*tmp_4);
      real_t tmp_7 = 0.26041666666666669*abs_det_jac_affine*(jac_affine_inv_0_1*tmp_6 + jac_affine_inv_1_1*tmp_6);
      real_t tmp_9 = 0.26041666666666669*abs_det_jac_affine*(jac_affine_inv_0_1*tmp_8 + jac_affine_inv_1_1*tmp_8);
      real_t a_0_0 = -0.33333333333333343*tmp_3 - 0.20000000000000007*tmp_5 - 0.20000000000000001*tmp_7 - 0.60000000000000009*tmp_9;
      real_t a_0_1 = -0.33333333333333331*tmp_3 - 0.20000000000000001*tmp_5 - 0.59999999999999998*tmp_7 - 0.20000000000000001*tmp_9;
      real_t a_0_2 = -0.33333333333333331*tmp_3 - 0.59999999999999998*tmp_5 - 0.20000000000000001*tmp_7 - 0.20000000000000001*tmp_9;
      elMat(0,0) = a_0_0;
      elMat(0,1) = a_0_1;
      elMat(0,2) = a_0_2;
   }



} // namespace forms
} // namespace hyteg

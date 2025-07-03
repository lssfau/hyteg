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

#include "p2_plus_bubble_to_dg1_div_affine_q3.hpp"

namespace hyteg {
namespace forms {

   void p2_plus_bubble_to_dg1_div_0_affine_q3::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 7 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_2 = 1.3333333333333333;
      real_t tmp_3 = 1.3333333333333333;
      real_t tmp_4 = tmp_2 + tmp_3 - 3;
      real_t tmp_6 = -0.66666666666666674;
      real_t tmp_9 = 0.80000000000000004;
      real_t tmp_10 = 2.3999999999999999;
      real_t tmp_11 = tmp_10 + tmp_9 - 3;
      real_t tmp_13 = -0.40000000000000002;
      real_t tmp_16 = 2.3999999999999999;
      real_t tmp_17 = 0.80000000000000004;
      real_t tmp_18 = tmp_16 + tmp_17 - 3;
      real_t tmp_20 = -0.80000000000000004;
      real_t tmp_23 = 0.80000000000000004;
      real_t tmp_24 = 0.80000000000000004;
      real_t tmp_25 = tmp_23 + tmp_24 - 3;
      real_t tmp_27 = -0.80000000000000004;
      real_t jac_affine_0_0 = -p_affine_0_0 + p_affine_1_0;
      real_t jac_affine_0_1 = -p_affine_0_0 + p_affine_2_0;
      real_t jac_affine_1_0 = -p_affine_0_1 + p_affine_1_1;
      real_t jac_affine_1_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_0 = jac_affine_0_0*jac_affine_1_1 - jac_affine_0_1*jac_affine_1_0;
      real_t tmp_1 = 1.0 / (tmp_0);
      real_t jac_affine_inv_0_0 = jac_affine_1_1*tmp_1;
      real_t tmp_30 = jac_affine_inv_0_0*(tmp_2 - 1);
      real_t tmp_31 = jac_affine_inv_0_0*(tmp_9 - 1);
      real_t tmp_32 = jac_affine_inv_0_0*(tmp_16 - 1);
      real_t tmp_33 = jac_affine_inv_0_0*(tmp_23 - 1);
      real_t tmp_38 = jac_affine_inv_0_0*tmp_3;
      real_t tmp_41 = jac_affine_inv_0_0*tmp_10;
      real_t tmp_44 = jac_affine_inv_0_0*tmp_17;
      real_t tmp_47 = jac_affine_inv_0_0*tmp_24;
      real_t tmp_58 = 27*jac_affine_inv_0_0;
      real_t jac_affine_inv_1_0 = -jac_affine_1_0*tmp_1;
      real_t tmp_5 = jac_affine_inv_0_0*tmp_4 + jac_affine_inv_1_0*tmp_4;
      real_t tmp_12 = jac_affine_inv_0_0*tmp_11 + jac_affine_inv_1_0*tmp_11;
      real_t tmp_19 = jac_affine_inv_0_0*tmp_18 + jac_affine_inv_1_0*tmp_18;
      real_t tmp_26 = jac_affine_inv_0_0*tmp_25 + jac_affine_inv_1_0*tmp_25;
      real_t tmp_34 = jac_affine_inv_1_0*(tmp_3 - 1);
      real_t tmp_35 = jac_affine_inv_1_0*(tmp_10 - 1);
      real_t tmp_36 = jac_affine_inv_1_0*(tmp_17 - 1);
      real_t tmp_37 = jac_affine_inv_1_0*(tmp_24 - 1);
      real_t tmp_39 = jac_affine_inv_1_0*tmp_2;
      real_t tmp_40 = tmp_38 + tmp_39;
      real_t tmp_42 = jac_affine_inv_1_0*tmp_9;
      real_t tmp_43 = tmp_41 + tmp_42;
      real_t tmp_45 = jac_affine_inv_1_0*tmp_16;
      real_t tmp_46 = tmp_44 + tmp_45;
      real_t tmp_48 = jac_affine_inv_1_0*tmp_23;
      real_t tmp_49 = tmp_47 + tmp_48;
      real_t tmp_50 = jac_affine_inv_1_0*(1.3333333333333335 - tmp_2) - tmp_38;
      real_t tmp_51 = jac_affine_inv_1_0*(-tmp_9 - 0.79999999999999982) - tmp_41;
      real_t tmp_52 = jac_affine_inv_1_0*(2.3999999999999999 - tmp_16) - tmp_44;
      real_t tmp_53 = jac_affine_inv_1_0*(2.3999999999999999 - tmp_23) - tmp_47;
      real_t tmp_54 = jac_affine_inv_0_0*(1.3333333333333335 - tmp_3) - tmp_39;
      real_t tmp_55 = jac_affine_inv_0_0*(2.3999999999999999 - tmp_10) - tmp_42;
      real_t tmp_56 = jac_affine_inv_0_0*(-tmp_17 - 0.79999999999999982) - tmp_45;
      real_t tmp_57 = jac_affine_inv_0_0*(2.3999999999999999 - tmp_24) - tmp_48;
      real_t tmp_59 = 27*jac_affine_inv_1_0;
      real_t tmp_60 = 0.33333333333333331*tmp_58*(-tmp_6 - 0.66666666666666663) + 3.7007434154171883e-17*tmp_59;
      real_t tmp_61 = 0.59999999999999998*tmp_58*(-tmp_13 - 0.40000000000000002) - 0.079999999999999988*tmp_59;
      real_t tmp_62 = 0.20000000000000001*tmp_58*(-tmp_20 - 1.2);
      real_t tmp_63 = 0.20000000000000001*tmp_58*(-tmp_27 - 0.40000000000000002) + 0.080000000000000016*tmp_59;
      real_t abs_det_jac_affine = std::abs(tmp_0);
      real_t tmp_7 = -0.28125*abs_det_jac_affine;
      real_t tmp_8 = tmp_7*(-tmp_6 - 0.33333333333333331);
      real_t tmp_14 = 0.26041666666666669*abs_det_jac_affine;
      real_t tmp_15 = tmp_14*(-tmp_13 - 0.20000000000000001);
      real_t tmp_21 = 0.26041666666666669*abs_det_jac_affine;
      real_t tmp_22 = tmp_21*(-tmp_20 - 0.59999999999999998);
      real_t tmp_28 = 0.26041666666666669*abs_det_jac_affine;
      real_t tmp_29 = tmp_28*(-tmp_27 - 0.20000000000000001);
      real_t tmp_64 = 0.33333333333333331*tmp_7;
      real_t tmp_65 = 0.20000000000000001*tmp_14;
      real_t tmp_66 = 0.59999999999999998*tmp_21;
      real_t tmp_67 = 0.20000000000000001*tmp_28;
      real_t tmp_68 = 0.33333333333333331*tmp_7;
      real_t tmp_69 = 0.59999999999999998*tmp_14;
      real_t tmp_70 = 0.20000000000000001*tmp_21;
      real_t tmp_71 = 0.20000000000000001*tmp_28;
      real_t a_0_0 = -tmp_12*tmp_15 - tmp_19*tmp_22 - tmp_26*tmp_29 - tmp_5*tmp_8;
      real_t a_0_1 = -tmp_15*tmp_31 - tmp_22*tmp_32 - tmp_29*tmp_33 - tmp_30*tmp_8;
      real_t a_0_2 = -tmp_15*tmp_35 - tmp_22*tmp_36 - tmp_29*tmp_37 - tmp_34*tmp_8;
      real_t a_0_3 = -tmp_15*tmp_43 - tmp_22*tmp_46 - tmp_29*tmp_49 - tmp_40*tmp_8;
      real_t a_0_4 = -tmp_15*tmp_51 - tmp_22*tmp_52 - tmp_29*tmp_53 - tmp_50*tmp_8;
      real_t a_0_5 = -tmp_15*tmp_55 - tmp_22*tmp_56 - tmp_29*tmp_57 - tmp_54*tmp_8;
      real_t a_0_6 = -tmp_15*tmp_61 - tmp_22*tmp_62 - tmp_29*tmp_63 - tmp_60*tmp_8;
      real_t a_1_0 = -tmp_12*tmp_65 - tmp_19*tmp_66 - tmp_26*tmp_67 - tmp_5*tmp_64;
      real_t a_1_1 = -tmp_30*tmp_64 - tmp_31*tmp_65 - tmp_32*tmp_66 - tmp_33*tmp_67;
      real_t a_1_2 = -tmp_34*tmp_64 - tmp_35*tmp_65 - tmp_36*tmp_66 - tmp_37*tmp_67;
      real_t a_1_3 = -tmp_40*tmp_64 - tmp_43*tmp_65 - tmp_46*tmp_66 - tmp_49*tmp_67;
      real_t a_1_4 = -tmp_50*tmp_64 - tmp_51*tmp_65 - tmp_52*tmp_66 - tmp_53*tmp_67;
      real_t a_1_5 = -tmp_54*tmp_64 - tmp_55*tmp_65 - tmp_56*tmp_66 - tmp_57*tmp_67;
      real_t a_1_6 = -tmp_60*tmp_64 - tmp_61*tmp_65 - tmp_62*tmp_66 - tmp_63*tmp_67;
      real_t a_2_0 = -tmp_12*tmp_69 - tmp_19*tmp_70 - tmp_26*tmp_71 - tmp_5*tmp_68;
      real_t a_2_1 = -tmp_30*tmp_68 - tmp_31*tmp_69 - tmp_32*tmp_70 - tmp_33*tmp_71;
      real_t a_2_2 = -tmp_34*tmp_68 - tmp_35*tmp_69 - tmp_36*tmp_70 - tmp_37*tmp_71;
      real_t a_2_3 = -tmp_40*tmp_68 - tmp_43*tmp_69 - tmp_46*tmp_70 - tmp_49*tmp_71;
      real_t a_2_4 = -tmp_50*tmp_68 - tmp_51*tmp_69 - tmp_52*tmp_70 - tmp_53*tmp_71;
      real_t a_2_5 = -tmp_54*tmp_68 - tmp_55*tmp_69 - tmp_56*tmp_70 - tmp_57*tmp_71;
      real_t a_2_6 = -tmp_60*tmp_68 - tmp_61*tmp_69 - tmp_62*tmp_70 - tmp_63*tmp_71;
      elMat(0,0) = a_0_0;
      elMat(0,1) = a_0_1;
      elMat(0,2) = a_0_2;
      elMat(0,3) = a_0_3;
      elMat(0,4) = a_0_4;
      elMat(0,5) = a_0_5;
      elMat(0,6) = a_0_6;
      elMat(1,0) = a_1_0;
      elMat(1,1) = a_1_1;
      elMat(1,2) = a_1_2;
      elMat(1,3) = a_1_3;
      elMat(1,4) = a_1_4;
      elMat(1,5) = a_1_5;
      elMat(1,6) = a_1_6;
      elMat(2,0) = a_2_0;
      elMat(2,1) = a_2_1;
      elMat(2,2) = a_2_2;
      elMat(2,3) = a_2_3;
      elMat(2,4) = a_2_4;
      elMat(2,5) = a_2_5;
      elMat(2,6) = a_2_6;
   }

   void p2_plus_bubble_to_dg1_div_0_affine_q3::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 7 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_2 = 1.3333333333333333;
      real_t tmp_3 = 1.3333333333333333;
      real_t tmp_4 = tmp_2 + tmp_3 - 3;
      real_t tmp_5 = -0.66666666666666674;
      real_t tmp_7 = 0.80000000000000004;
      real_t tmp_8 = 2.3999999999999999;
      real_t tmp_9 = tmp_7 + tmp_8 - 3;
      real_t tmp_10 = -0.40000000000000002;
      real_t tmp_12 = 2.3999999999999999;
      real_t tmp_13 = 0.80000000000000004;
      real_t tmp_14 = tmp_12 + tmp_13 - 3;
      real_t tmp_15 = -0.80000000000000004;
      real_t tmp_17 = 0.80000000000000004;
      real_t tmp_18 = 0.80000000000000004;
      real_t tmp_19 = tmp_17 + tmp_18 - 3;
      real_t tmp_20 = -0.80000000000000004;
      real_t jac_affine_0_0 = -p_affine_0_0 + p_affine_1_0;
      real_t jac_affine_0_1 = -p_affine_0_0 + p_affine_2_0;
      real_t jac_affine_1_0 = -p_affine_0_1 + p_affine_1_1;
      real_t jac_affine_1_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_0 = jac_affine_0_0*jac_affine_1_1 - jac_affine_0_1*jac_affine_1_0;
      real_t tmp_1 = 1.0 / (tmp_0);
      real_t jac_affine_inv_0_0 = jac_affine_1_1*tmp_1;
      real_t tmp_22 = jac_affine_inv_0_0*tmp_3;
      real_t tmp_24 = jac_affine_inv_0_0*tmp_8;
      real_t tmp_26 = jac_affine_inv_0_0*tmp_13;
      real_t tmp_28 = jac_affine_inv_0_0*tmp_18;
      real_t tmp_30 = 27*jac_affine_inv_0_0;
      real_t jac_affine_inv_1_0 = -jac_affine_1_0*tmp_1;
      real_t tmp_23 = jac_affine_inv_1_0*tmp_2;
      real_t tmp_25 = jac_affine_inv_1_0*tmp_7;
      real_t tmp_27 = jac_affine_inv_1_0*tmp_12;
      real_t tmp_29 = jac_affine_inv_1_0*tmp_17;
      real_t tmp_31 = 27*jac_affine_inv_1_0;
      real_t abs_det_jac_affine = std::abs(tmp_0);
      real_t tmp_6 = -0.28125*abs_det_jac_affine*(-tmp_5 - 0.33333333333333331);
      real_t tmp_11 = 0.26041666666666669*abs_det_jac_affine*(-tmp_10 - 0.20000000000000001);
      real_t tmp_16 = 0.26041666666666669*abs_det_jac_affine*(-tmp_15 - 0.59999999999999998);
      real_t tmp_21 = 0.26041666666666669*abs_det_jac_affine*(-tmp_20 - 0.20000000000000001);
      real_t a_0_0 = -tmp_11*(jac_affine_inv_0_0*tmp_9 + jac_affine_inv_1_0*tmp_9) - tmp_16*(jac_affine_inv_0_0*tmp_14 + jac_affine_inv_1_0*tmp_14) - tmp_21*(jac_affine_inv_0_0*tmp_19 + jac_affine_inv_1_0*tmp_19) - tmp_6*(jac_affine_inv_0_0*tmp_4 + jac_affine_inv_1_0*tmp_4);
      real_t a_0_1 = -jac_affine_inv_0_0*tmp_11*(tmp_7 - 1) - jac_affine_inv_0_0*tmp_16*(tmp_12 - 1) - jac_affine_inv_0_0*tmp_21*(tmp_17 - 1) - jac_affine_inv_0_0*tmp_6*(tmp_2 - 1);
      real_t a_0_2 = -jac_affine_inv_1_0*tmp_11*(tmp_8 - 1) - jac_affine_inv_1_0*tmp_16*(tmp_13 - 1) - jac_affine_inv_1_0*tmp_21*(tmp_18 - 1) - jac_affine_inv_1_0*tmp_6*(tmp_3 - 1);
      real_t a_0_3 = -tmp_11*(tmp_24 + tmp_25) - tmp_16*(tmp_26 + tmp_27) - tmp_21*(tmp_28 + tmp_29) - tmp_6*(tmp_22 + tmp_23);
      real_t a_0_4 = -tmp_11*(jac_affine_inv_1_0*(-tmp_7 - 0.79999999999999982) - tmp_24) - tmp_16*(jac_affine_inv_1_0*(2.3999999999999999 - tmp_12) - tmp_26) - tmp_21*(jac_affine_inv_1_0*(2.3999999999999999 - tmp_17) - tmp_28) - tmp_6*(jac_affine_inv_1_0*(1.3333333333333335 - tmp_2) - tmp_22);
      real_t a_0_5 = -tmp_11*(jac_affine_inv_0_0*(2.3999999999999999 - tmp_8) - tmp_25) - tmp_16*(jac_affine_inv_0_0*(-tmp_13 - 0.79999999999999982) - tmp_27) - tmp_21*(jac_affine_inv_0_0*(2.3999999999999999 - tmp_18) - tmp_29) - tmp_6*(jac_affine_inv_0_0*(1.3333333333333335 - tmp_3) - tmp_23);
      real_t a_0_6 = -tmp_11*(0.59999999999999998*tmp_30*(-tmp_10 - 0.40000000000000002) - 0.079999999999999988*tmp_31) - 0.20000000000000001*tmp_16*tmp_30*(-tmp_15 - 1.2) - tmp_21*(0.20000000000000001*tmp_30*(-tmp_20 - 0.40000000000000002) + 0.080000000000000016*tmp_31) - tmp_6*(0.33333333333333331*tmp_30*(-tmp_5 - 0.66666666666666663) + 3.7007434154171883e-17*tmp_31);
      elMat(0,0) = a_0_0;
      elMat(0,1) = a_0_1;
      elMat(0,2) = a_0_2;
      elMat(0,3) = a_0_3;
      elMat(0,4) = a_0_4;
      elMat(0,5) = a_0_5;
      elMat(0,6) = a_0_6;
   }

   void p2_plus_bubble_to_dg1_div_1_affine_q3::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 7 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_2 = 1.3333333333333333;
      real_t tmp_3 = 1.3333333333333333;
      real_t tmp_4 = tmp_2 + tmp_3 - 3;
      real_t tmp_6 = -0.66666666666666674;
      real_t tmp_9 = 0.80000000000000004;
      real_t tmp_10 = 2.3999999999999999;
      real_t tmp_11 = tmp_10 + tmp_9 - 3;
      real_t tmp_13 = -0.40000000000000002;
      real_t tmp_16 = 2.3999999999999999;
      real_t tmp_17 = 0.80000000000000004;
      real_t tmp_18 = tmp_16 + tmp_17 - 3;
      real_t tmp_20 = -0.80000000000000004;
      real_t tmp_23 = 0.80000000000000004;
      real_t tmp_24 = 0.80000000000000004;
      real_t tmp_25 = tmp_23 + tmp_24 - 3;
      real_t tmp_27 = -0.80000000000000004;
      real_t jac_affine_0_0 = -p_affine_0_0 + p_affine_1_0;
      real_t jac_affine_0_1 = -p_affine_0_0 + p_affine_2_0;
      real_t jac_affine_1_0 = -p_affine_0_1 + p_affine_1_1;
      real_t jac_affine_1_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_0 = jac_affine_0_0*jac_affine_1_1 - jac_affine_0_1*jac_affine_1_0;
      real_t tmp_1 = 1.0 / (tmp_0);
      real_t jac_affine_inv_0_1 = -jac_affine_0_1*tmp_1;
      real_t tmp_30 = jac_affine_inv_0_1*(tmp_2 - 1);
      real_t tmp_31 = jac_affine_inv_0_1*(tmp_9 - 1);
      real_t tmp_32 = jac_affine_inv_0_1*(tmp_16 - 1);
      real_t tmp_33 = jac_affine_inv_0_1*(tmp_23 - 1);
      real_t tmp_38 = jac_affine_inv_0_1*tmp_3;
      real_t tmp_41 = jac_affine_inv_0_1*tmp_10;
      real_t tmp_44 = jac_affine_inv_0_1*tmp_17;
      real_t tmp_47 = jac_affine_inv_0_1*tmp_24;
      real_t tmp_58 = 27*jac_affine_inv_0_1;
      real_t jac_affine_inv_1_1 = jac_affine_0_0*tmp_1;
      real_t tmp_5 = jac_affine_inv_0_1*tmp_4 + jac_affine_inv_1_1*tmp_4;
      real_t tmp_12 = jac_affine_inv_0_1*tmp_11 + jac_affine_inv_1_1*tmp_11;
      real_t tmp_19 = jac_affine_inv_0_1*tmp_18 + jac_affine_inv_1_1*tmp_18;
      real_t tmp_26 = jac_affine_inv_0_1*tmp_25 + jac_affine_inv_1_1*tmp_25;
      real_t tmp_34 = jac_affine_inv_1_1*(tmp_3 - 1);
      real_t tmp_35 = jac_affine_inv_1_1*(tmp_10 - 1);
      real_t tmp_36 = jac_affine_inv_1_1*(tmp_17 - 1);
      real_t tmp_37 = jac_affine_inv_1_1*(tmp_24 - 1);
      real_t tmp_39 = jac_affine_inv_1_1*tmp_2;
      real_t tmp_40 = tmp_38 + tmp_39;
      real_t tmp_42 = jac_affine_inv_1_1*tmp_9;
      real_t tmp_43 = tmp_41 + tmp_42;
      real_t tmp_45 = jac_affine_inv_1_1*tmp_16;
      real_t tmp_46 = tmp_44 + tmp_45;
      real_t tmp_48 = jac_affine_inv_1_1*tmp_23;
      real_t tmp_49 = tmp_47 + tmp_48;
      real_t tmp_50 = jac_affine_inv_1_1*(1.3333333333333335 - tmp_2) - tmp_38;
      real_t tmp_51 = jac_affine_inv_1_1*(-tmp_9 - 0.79999999999999982) - tmp_41;
      real_t tmp_52 = jac_affine_inv_1_1*(2.3999999999999999 - tmp_16) - tmp_44;
      real_t tmp_53 = jac_affine_inv_1_1*(2.3999999999999999 - tmp_23) - tmp_47;
      real_t tmp_54 = jac_affine_inv_0_1*(1.3333333333333335 - tmp_3) - tmp_39;
      real_t tmp_55 = jac_affine_inv_0_1*(2.3999999999999999 - tmp_10) - tmp_42;
      real_t tmp_56 = jac_affine_inv_0_1*(-tmp_17 - 0.79999999999999982) - tmp_45;
      real_t tmp_57 = jac_affine_inv_0_1*(2.3999999999999999 - tmp_24) - tmp_48;
      real_t tmp_59 = 27*jac_affine_inv_1_1;
      real_t tmp_60 = 0.33333333333333331*tmp_58*(-tmp_6 - 0.66666666666666663) + 3.7007434154171883e-17*tmp_59;
      real_t tmp_61 = 0.59999999999999998*tmp_58*(-tmp_13 - 0.40000000000000002) - 0.079999999999999988*tmp_59;
      real_t tmp_62 = 0.20000000000000001*tmp_58*(-tmp_20 - 1.2);
      real_t tmp_63 = 0.20000000000000001*tmp_58*(-tmp_27 - 0.40000000000000002) + 0.080000000000000016*tmp_59;
      real_t abs_det_jac_affine = std::abs(tmp_0);
      real_t tmp_7 = -0.28125*abs_det_jac_affine;
      real_t tmp_8 = tmp_7*(-tmp_6 - 0.33333333333333331);
      real_t tmp_14 = 0.26041666666666669*abs_det_jac_affine;
      real_t tmp_15 = tmp_14*(-tmp_13 - 0.20000000000000001);
      real_t tmp_21 = 0.26041666666666669*abs_det_jac_affine;
      real_t tmp_22 = tmp_21*(-tmp_20 - 0.59999999999999998);
      real_t tmp_28 = 0.26041666666666669*abs_det_jac_affine;
      real_t tmp_29 = tmp_28*(-tmp_27 - 0.20000000000000001);
      real_t tmp_64 = 0.33333333333333331*tmp_7;
      real_t tmp_65 = 0.20000000000000001*tmp_14;
      real_t tmp_66 = 0.59999999999999998*tmp_21;
      real_t tmp_67 = 0.20000000000000001*tmp_28;
      real_t tmp_68 = 0.33333333333333331*tmp_7;
      real_t tmp_69 = 0.59999999999999998*tmp_14;
      real_t tmp_70 = 0.20000000000000001*tmp_21;
      real_t tmp_71 = 0.20000000000000001*tmp_28;
      real_t a_0_0 = -tmp_12*tmp_15 - tmp_19*tmp_22 - tmp_26*tmp_29 - tmp_5*tmp_8;
      real_t a_0_1 = -tmp_15*tmp_31 - tmp_22*tmp_32 - tmp_29*tmp_33 - tmp_30*tmp_8;
      real_t a_0_2 = -tmp_15*tmp_35 - tmp_22*tmp_36 - tmp_29*tmp_37 - tmp_34*tmp_8;
      real_t a_0_3 = -tmp_15*tmp_43 - tmp_22*tmp_46 - tmp_29*tmp_49 - tmp_40*tmp_8;
      real_t a_0_4 = -tmp_15*tmp_51 - tmp_22*tmp_52 - tmp_29*tmp_53 - tmp_50*tmp_8;
      real_t a_0_5 = -tmp_15*tmp_55 - tmp_22*tmp_56 - tmp_29*tmp_57 - tmp_54*tmp_8;
      real_t a_0_6 = -tmp_15*tmp_61 - tmp_22*tmp_62 - tmp_29*tmp_63 - tmp_60*tmp_8;
      real_t a_1_0 = -tmp_12*tmp_65 - tmp_19*tmp_66 - tmp_26*tmp_67 - tmp_5*tmp_64;
      real_t a_1_1 = -tmp_30*tmp_64 - tmp_31*tmp_65 - tmp_32*tmp_66 - tmp_33*tmp_67;
      real_t a_1_2 = -tmp_34*tmp_64 - tmp_35*tmp_65 - tmp_36*tmp_66 - tmp_37*tmp_67;
      real_t a_1_3 = -tmp_40*tmp_64 - tmp_43*tmp_65 - tmp_46*tmp_66 - tmp_49*tmp_67;
      real_t a_1_4 = -tmp_50*tmp_64 - tmp_51*tmp_65 - tmp_52*tmp_66 - tmp_53*tmp_67;
      real_t a_1_5 = -tmp_54*tmp_64 - tmp_55*tmp_65 - tmp_56*tmp_66 - tmp_57*tmp_67;
      real_t a_1_6 = -tmp_60*tmp_64 - tmp_61*tmp_65 - tmp_62*tmp_66 - tmp_63*tmp_67;
      real_t a_2_0 = -tmp_12*tmp_69 - tmp_19*tmp_70 - tmp_26*tmp_71 - tmp_5*tmp_68;
      real_t a_2_1 = -tmp_30*tmp_68 - tmp_31*tmp_69 - tmp_32*tmp_70 - tmp_33*tmp_71;
      real_t a_2_2 = -tmp_34*tmp_68 - tmp_35*tmp_69 - tmp_36*tmp_70 - tmp_37*tmp_71;
      real_t a_2_3 = -tmp_40*tmp_68 - tmp_43*tmp_69 - tmp_46*tmp_70 - tmp_49*tmp_71;
      real_t a_2_4 = -tmp_50*tmp_68 - tmp_51*tmp_69 - tmp_52*tmp_70 - tmp_53*tmp_71;
      real_t a_2_5 = -tmp_54*tmp_68 - tmp_55*tmp_69 - tmp_56*tmp_70 - tmp_57*tmp_71;
      real_t a_2_6 = -tmp_60*tmp_68 - tmp_61*tmp_69 - tmp_62*tmp_70 - tmp_63*tmp_71;
      elMat(0,0) = a_0_0;
      elMat(0,1) = a_0_1;
      elMat(0,2) = a_0_2;
      elMat(0,3) = a_0_3;
      elMat(0,4) = a_0_4;
      elMat(0,5) = a_0_5;
      elMat(0,6) = a_0_6;
      elMat(1,0) = a_1_0;
      elMat(1,1) = a_1_1;
      elMat(1,2) = a_1_2;
      elMat(1,3) = a_1_3;
      elMat(1,4) = a_1_4;
      elMat(1,5) = a_1_5;
      elMat(1,6) = a_1_6;
      elMat(2,0) = a_2_0;
      elMat(2,1) = a_2_1;
      elMat(2,2) = a_2_2;
      elMat(2,3) = a_2_3;
      elMat(2,4) = a_2_4;
      elMat(2,5) = a_2_5;
      elMat(2,6) = a_2_6;
   }

   void p2_plus_bubble_to_dg1_div_1_affine_q3::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 7 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_2 = 1.3333333333333333;
      real_t tmp_3 = 1.3333333333333333;
      real_t tmp_4 = tmp_2 + tmp_3 - 3;
      real_t tmp_5 = -0.66666666666666674;
      real_t tmp_7 = 0.80000000000000004;
      real_t tmp_8 = 2.3999999999999999;
      real_t tmp_9 = tmp_7 + tmp_8 - 3;
      real_t tmp_10 = -0.40000000000000002;
      real_t tmp_12 = 2.3999999999999999;
      real_t tmp_13 = 0.80000000000000004;
      real_t tmp_14 = tmp_12 + tmp_13 - 3;
      real_t tmp_15 = -0.80000000000000004;
      real_t tmp_17 = 0.80000000000000004;
      real_t tmp_18 = 0.80000000000000004;
      real_t tmp_19 = tmp_17 + tmp_18 - 3;
      real_t tmp_20 = -0.80000000000000004;
      real_t jac_affine_0_0 = -p_affine_0_0 + p_affine_1_0;
      real_t jac_affine_0_1 = -p_affine_0_0 + p_affine_2_0;
      real_t jac_affine_1_0 = -p_affine_0_1 + p_affine_1_1;
      real_t jac_affine_1_1 = -p_affine_0_1 + p_affine_2_1;
      real_t tmp_0 = jac_affine_0_0*jac_affine_1_1 - jac_affine_0_1*jac_affine_1_0;
      real_t tmp_1 = 1.0 / (tmp_0);
      real_t jac_affine_inv_0_1 = -jac_affine_0_1*tmp_1;
      real_t tmp_22 = jac_affine_inv_0_1*tmp_3;
      real_t tmp_24 = jac_affine_inv_0_1*tmp_8;
      real_t tmp_26 = jac_affine_inv_0_1*tmp_13;
      real_t tmp_28 = jac_affine_inv_0_1*tmp_18;
      real_t tmp_30 = 27*jac_affine_inv_0_1;
      real_t jac_affine_inv_1_1 = jac_affine_0_0*tmp_1;
      real_t tmp_23 = jac_affine_inv_1_1*tmp_2;
      real_t tmp_25 = jac_affine_inv_1_1*tmp_7;
      real_t tmp_27 = jac_affine_inv_1_1*tmp_12;
      real_t tmp_29 = jac_affine_inv_1_1*tmp_17;
      real_t tmp_31 = 27*jac_affine_inv_1_1;
      real_t abs_det_jac_affine = std::abs(tmp_0);
      real_t tmp_6 = -0.28125*abs_det_jac_affine*(-tmp_5 - 0.33333333333333331);
      real_t tmp_11 = 0.26041666666666669*abs_det_jac_affine*(-tmp_10 - 0.20000000000000001);
      real_t tmp_16 = 0.26041666666666669*abs_det_jac_affine*(-tmp_15 - 0.59999999999999998);
      real_t tmp_21 = 0.26041666666666669*abs_det_jac_affine*(-tmp_20 - 0.20000000000000001);
      real_t a_0_0 = -tmp_11*(jac_affine_inv_0_1*tmp_9 + jac_affine_inv_1_1*tmp_9) - tmp_16*(jac_affine_inv_0_1*tmp_14 + jac_affine_inv_1_1*tmp_14) - tmp_21*(jac_affine_inv_0_1*tmp_19 + jac_affine_inv_1_1*tmp_19) - tmp_6*(jac_affine_inv_0_1*tmp_4 + jac_affine_inv_1_1*tmp_4);
      real_t a_0_1 = -jac_affine_inv_0_1*tmp_11*(tmp_7 - 1) - jac_affine_inv_0_1*tmp_16*(tmp_12 - 1) - jac_affine_inv_0_1*tmp_21*(tmp_17 - 1) - jac_affine_inv_0_1*tmp_6*(tmp_2 - 1);
      real_t a_0_2 = -jac_affine_inv_1_1*tmp_11*(tmp_8 - 1) - jac_affine_inv_1_1*tmp_16*(tmp_13 - 1) - jac_affine_inv_1_1*tmp_21*(tmp_18 - 1) - jac_affine_inv_1_1*tmp_6*(tmp_3 - 1);
      real_t a_0_3 = -tmp_11*(tmp_24 + tmp_25) - tmp_16*(tmp_26 + tmp_27) - tmp_21*(tmp_28 + tmp_29) - tmp_6*(tmp_22 + tmp_23);
      real_t a_0_4 = -tmp_11*(jac_affine_inv_1_1*(-tmp_7 - 0.79999999999999982) - tmp_24) - tmp_16*(jac_affine_inv_1_1*(2.3999999999999999 - tmp_12) - tmp_26) - tmp_21*(jac_affine_inv_1_1*(2.3999999999999999 - tmp_17) - tmp_28) - tmp_6*(jac_affine_inv_1_1*(1.3333333333333335 - tmp_2) - tmp_22);
      real_t a_0_5 = -tmp_11*(jac_affine_inv_0_1*(2.3999999999999999 - tmp_8) - tmp_25) - tmp_16*(jac_affine_inv_0_1*(-tmp_13 - 0.79999999999999982) - tmp_27) - tmp_21*(jac_affine_inv_0_1*(2.3999999999999999 - tmp_18) - tmp_29) - tmp_6*(jac_affine_inv_0_1*(1.3333333333333335 - tmp_3) - tmp_23);
      real_t a_0_6 = -tmp_11*(0.59999999999999998*tmp_30*(-tmp_10 - 0.40000000000000002) - 0.079999999999999988*tmp_31) - 0.20000000000000001*tmp_16*tmp_30*(-tmp_15 - 1.2) - tmp_21*(0.20000000000000001*tmp_30*(-tmp_20 - 0.40000000000000002) + 0.080000000000000016*tmp_31) - tmp_6*(0.33333333333333331*tmp_30*(-tmp_5 - 0.66666666666666663) + 3.7007434154171883e-17*tmp_31);
      elMat(0,0) = a_0_0;
      elMat(0,1) = a_0_1;
      elMat(0,2) = a_0_2;
      elMat(0,3) = a_0_3;
      elMat(0,4) = a_0_4;
      elMat(0,5) = a_0_5;
      elMat(0,6) = a_0_6;
   }



} // namespace forms
} // namespace hyteg

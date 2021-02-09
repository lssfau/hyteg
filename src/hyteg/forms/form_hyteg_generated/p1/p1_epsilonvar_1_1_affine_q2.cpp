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

#include "p1_epsilonvar_1_1_affine_q2.hpp"

namespace hyteg {
namespace forms {

   void p1_epsilonvar_1_1_affine_q2::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Scalar_Variable_Coefficient_2D_0_0 = 0;
      real_t Scalar_Variable_Coefficient_2D_1_0 = 0;
      real_t Scalar_Variable_Coefficient_2D_2_0 = 0;
      real_t q_p_0_0 = 0.16666666666666666;
      real_t q_p_0_1 = 0.66666666666666663;
      real_t q_p_1_0 = 0.66666666666666663;
      real_t q_p_1_1 = 0.16666666666666666;
      real_t q_p_2_0 = 0.16666666666666666;
      real_t q_p_2_1 = 0.16666666666666666;
      real_t w_p_0 = 0.16666666666666666;
      real_t w_p_1 = 0.16666666666666666;
      real_t w_p_2 = 0.16666666666666666;
      Scalar_Variable_Coefficient_2D( p_affine_0_0 + q_p_0_0*(-p_affine_0_0 + p_affine_1_0) + q_p_0_1*(-p_affine_0_0 + p_affine_2_0), p_affine_0_1 + q_p_0_0*(-p_affine_0_1 + p_affine_1_1) + q_p_0_1*(-p_affine_0_1 + p_affine_2_1), &Scalar_Variable_Coefficient_2D_0_0 );
      Scalar_Variable_Coefficient_2D( p_affine_0_0 + q_p_1_0*(-p_affine_0_0 + p_affine_1_0) + q_p_1_1*(-p_affine_0_0 + p_affine_2_0), p_affine_0_1 + q_p_1_0*(-p_affine_0_1 + p_affine_1_1) + q_p_1_1*(-p_affine_0_1 + p_affine_2_1), &Scalar_Variable_Coefficient_2D_1_0 );
      Scalar_Variable_Coefficient_2D( p_affine_0_0 + q_p_2_0*(-p_affine_0_0 + p_affine_1_0) + q_p_2_1*(-p_affine_0_0 + p_affine_2_0), p_affine_0_1 + q_p_2_0*(-p_affine_0_1 + p_affine_1_1) + q_p_2_1*(-p_affine_0_1 + p_affine_2_1), &Scalar_Variable_Coefficient_2D_2_0 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0);
      real_t tmp_5 = 1.0 / (tmp_4);
      real_t tmp_6 = 1.0*tmp_5;
      real_t tmp_7 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_8 = -tmp_1*tmp_6 - tmp_6*tmp_7;
      real_t tmp_9 = 2*(tmp_8*tmp_8);
      real_t tmp_10 = tmp_3*tmp_5;
      real_t tmp_11 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = 1.0*((-tmp_10 - tmp_12)*(-tmp_10 - tmp_12));
      real_t tmp_14 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_15 = tmp_14*w_p_0;
      real_t tmp_16 = tmp_14*w_p_1;
      real_t tmp_17 = tmp_14*w_p_2;
      real_t tmp_18 = 2.0*Scalar_Variable_Coefficient_2D_0_0;
      real_t tmp_19 = tmp_5*tmp_8;
      real_t tmp_20 = tmp_19*tmp_7;
      real_t tmp_21 = -0.5*tmp_10 - 0.5*tmp_12;
      real_t tmp_22 = tmp_10*tmp_21;
      real_t tmp_23 = 2.0*Scalar_Variable_Coefficient_2D_1_0;
      real_t tmp_24 = 2.0*Scalar_Variable_Coefficient_2D_2_0;
      real_t tmp_25 = tmp_15*(tmp_18*tmp_20 + tmp_18*tmp_22) + tmp_16*(tmp_20*tmp_23 + tmp_22*tmp_23) + tmp_17*(tmp_20*tmp_24 + tmp_22*tmp_24);
      real_t tmp_26 = tmp_1*tmp_19;
      real_t tmp_27 = tmp_12*tmp_21;
      real_t tmp_28 = tmp_15*(tmp_18*tmp_26 + tmp_18*tmp_27) + tmp_16*(tmp_23*tmp_26 + tmp_23*tmp_27) + tmp_17*(tmp_24*tmp_26 + tmp_24*tmp_27);
      real_t tmp_29 = 1.0 / (tmp_4*tmp_4);
      real_t tmp_30 = Scalar_Variable_Coefficient_2D_0_0*tmp_29;
      real_t tmp_31 = 2.0*(tmp_7*tmp_7);
      real_t tmp_32 = 1.0*(tmp_3*tmp_3);
      real_t tmp_33 = Scalar_Variable_Coefficient_2D_1_0*tmp_29;
      real_t tmp_34 = Scalar_Variable_Coefficient_2D_2_0*tmp_29;
      real_t tmp_35 = 2.0*tmp_1*tmp_7;
      real_t tmp_36 = 1.0*tmp_11*tmp_3;
      real_t tmp_37 = tmp_15*(tmp_30*tmp_35 + tmp_30*tmp_36) + tmp_16*(tmp_33*tmp_35 + tmp_33*tmp_36) + tmp_17*(tmp_34*tmp_35 + tmp_34*tmp_36);
      real_t tmp_38 = 2.0*(tmp_1*tmp_1);
      real_t tmp_39 = 1.0*(tmp_11*tmp_11);
      real_t a_0_0 = tmp_15*(Scalar_Variable_Coefficient_2D_0_0*tmp_13 + Scalar_Variable_Coefficient_2D_0_0*tmp_9) + tmp_16*(Scalar_Variable_Coefficient_2D_1_0*tmp_13 + Scalar_Variable_Coefficient_2D_1_0*tmp_9) + tmp_17*(Scalar_Variable_Coefficient_2D_2_0*tmp_13 + Scalar_Variable_Coefficient_2D_2_0*tmp_9);
      real_t a_0_1 = tmp_25;
      real_t a_0_2 = tmp_28;
      real_t a_1_0 = tmp_25;
      real_t a_1_1 = tmp_15*(tmp_30*tmp_31 + tmp_30*tmp_32) + tmp_16*(tmp_31*tmp_33 + tmp_32*tmp_33) + tmp_17*(tmp_31*tmp_34 + tmp_32*tmp_34);
      real_t a_1_2 = tmp_37;
      real_t a_2_0 = tmp_28;
      real_t a_2_1 = tmp_37;
      real_t a_2_2 = tmp_15*(tmp_30*tmp_38 + tmp_30*tmp_39) + tmp_16*(tmp_33*tmp_38 + tmp_33*tmp_39) + tmp_17*(tmp_34*tmp_38 + tmp_34*tmp_39);
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

   void p1_epsilonvar_1_1_affine_q2::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_1_0 = 0;
      real_t Scalar_Variable_Coefficient_3D_2_0 = 0;
      real_t Scalar_Variable_Coefficient_3D_3_0 = 0;
      real_t q_p_0_0 = 0.13819660112501059;
      real_t q_p_0_1 = 0.13819660112501059;
      real_t q_p_0_2 = 0.58541019662496829;
      real_t q_p_1_0 = 0.13819660112501059;
      real_t q_p_1_1 = 0.58541019662496829;
      real_t q_p_1_2 = 0.13819660112501059;
      real_t q_p_2_0 = 0.58541019662496829;
      real_t q_p_2_1 = 0.13819660112501059;
      real_t q_p_2_2 = 0.13819660112501059;
      real_t q_p_3_0 = 0.13819660112501059;
      real_t q_p_3_1 = 0.13819660112501059;
      real_t q_p_3_2 = 0.13819660112501059;
      real_t w_p_0 = 0.041666666666666657;
      real_t w_p_1 = 0.041666666666666657;
      real_t w_p_2 = 0.041666666666666657;
      real_t w_p_3 = 0.041666666666666657;
      Scalar_Variable_Coefficient_3D( p_affine_0_0 + q_p_0_0*(-p_affine_0_0 + p_affine_1_0) + q_p_0_1*(-p_affine_0_0 + p_affine_2_0) + q_p_0_2*(-p_affine_0_0 + p_affine_3_0), p_affine_0_1 + q_p_0_0*(-p_affine_0_1 + p_affine_1_1) + q_p_0_1*(-p_affine_0_1 + p_affine_2_1) + q_p_0_2*(-p_affine_0_1 + p_affine_3_1), p_affine_0_2 + q_p_0_0*(-p_affine_0_2 + p_affine_1_2) + q_p_0_1*(-p_affine_0_2 + p_affine_2_2) + q_p_0_2*(-p_affine_0_2 + p_affine_3_2), &Scalar_Variable_Coefficient_3D_0_0 );
      Scalar_Variable_Coefficient_3D( p_affine_0_0 + q_p_1_0*(-p_affine_0_0 + p_affine_1_0) + q_p_1_1*(-p_affine_0_0 + p_affine_2_0) + q_p_1_2*(-p_affine_0_0 + p_affine_3_0), p_affine_0_1 + q_p_1_0*(-p_affine_0_1 + p_affine_1_1) + q_p_1_1*(-p_affine_0_1 + p_affine_2_1) + q_p_1_2*(-p_affine_0_1 + p_affine_3_1), p_affine_0_2 + q_p_1_0*(-p_affine_0_2 + p_affine_1_2) + q_p_1_1*(-p_affine_0_2 + p_affine_2_2) + q_p_1_2*(-p_affine_0_2 + p_affine_3_2), &Scalar_Variable_Coefficient_3D_1_0 );
      Scalar_Variable_Coefficient_3D( p_affine_0_0 + q_p_2_0*(-p_affine_0_0 + p_affine_1_0) + q_p_2_1*(-p_affine_0_0 + p_affine_2_0) + q_p_2_2*(-p_affine_0_0 + p_affine_3_0), p_affine_0_1 + q_p_2_0*(-p_affine_0_1 + p_affine_1_1) + q_p_2_1*(-p_affine_0_1 + p_affine_2_1) + q_p_2_2*(-p_affine_0_1 + p_affine_3_1), p_affine_0_2 + q_p_2_0*(-p_affine_0_2 + p_affine_1_2) + q_p_2_1*(-p_affine_0_2 + p_affine_2_2) + q_p_2_2*(-p_affine_0_2 + p_affine_3_2), &Scalar_Variable_Coefficient_3D_2_0 );
      Scalar_Variable_Coefficient_3D( p_affine_0_0 + q_p_3_0*(-p_affine_0_0 + p_affine_1_0) + q_p_3_1*(-p_affine_0_0 + p_affine_2_0) + q_p_3_2*(-p_affine_0_0 + p_affine_3_0), p_affine_0_1 + q_p_3_0*(-p_affine_0_1 + p_affine_1_1) + q_p_3_1*(-p_affine_0_1 + p_affine_2_1) + q_p_3_2*(-p_affine_0_1 + p_affine_3_1), p_affine_0_2 + q_p_3_0*(-p_affine_0_2 + p_affine_1_2) + q_p_3_1*(-p_affine_0_2 + p_affine_2_2) + q_p_3_2*(-p_affine_0_2 + p_affine_3_2), &Scalar_Variable_Coefficient_3D_3_0 );
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
      real_t tmp_26 = ((-tmp_21 - tmp_23 - tmp_25)*(-tmp_21 - tmp_23 - tmp_25));
      real_t tmp_27 = 1.0*Scalar_Variable_Coefficient_3D_0_0;
      real_t tmp_28 = tmp_12 - tmp_16;
      real_t tmp_29 = 1.0*tmp_20;
      real_t tmp_30 = -tmp_18 + tmp_9;
      real_t tmp_31 = tmp_15 - tmp_17;
      real_t tmp_32 = -tmp_28*tmp_29 - tmp_29*tmp_30 - tmp_29*tmp_31;
      real_t tmp_33 = 2*(tmp_32*tmp_32);
      real_t tmp_34 = -tmp_11*tmp_3 + tmp_14*tmp_5;
      real_t tmp_35 = tmp_20*tmp_34;
      real_t tmp_36 = tmp_10*tmp_11 - tmp_5*tmp_8;
      real_t tmp_37 = tmp_20*tmp_36;
      real_t tmp_38 = -tmp_10*tmp_14 + tmp_3*tmp_8;
      real_t tmp_39 = tmp_20*tmp_38;
      real_t tmp_40 = ((-tmp_35 - tmp_37 - tmp_39)*(-tmp_35 - tmp_37 - tmp_39));
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
      real_t tmp_53 = std::abs(p_affine_0_0*tmp_43 - p_affine_0_0*tmp_50 + p_affine_0_1*tmp_46 - p_affine_0_1*tmp_51 + p_affine_0_2*tmp_49 - p_affine_0_2*tmp_52 - p_affine_1_0*tmp_43 + p_affine_1_0*tmp_50 - p_affine_1_1*tmp_46 + p_affine_1_1*tmp_51 - p_affine_1_2*tmp_49 + p_affine_1_2*tmp_52 + p_affine_2_0*tmp_45 - p_affine_2_0*tmp_48 - p_affine_2_1*tmp_42 + p_affine_2_1*tmp_47 + p_affine_2_2*tmp_41 - p_affine_2_2*tmp_44 - p_affine_3_0*tmp_45 + p_affine_3_0*tmp_48 + p_affine_3_1*tmp_42 - p_affine_3_1*tmp_47 - p_affine_3_2*tmp_41 + p_affine_3_2*tmp_44);
      real_t tmp_54 = tmp_53*w_p_0;
      real_t tmp_55 = 1.0*Scalar_Variable_Coefficient_3D_1_0;
      real_t tmp_56 = tmp_53*w_p_1;
      real_t tmp_57 = 1.0*Scalar_Variable_Coefficient_3D_2_0;
      real_t tmp_58 = tmp_53*w_p_2;
      real_t tmp_59 = 1.0*Scalar_Variable_Coefficient_3D_3_0;
      real_t tmp_60 = tmp_53*w_p_3;
      real_t tmp_61 = 2.0*Scalar_Variable_Coefficient_3D_0_0;
      real_t tmp_62 = -0.5*tmp_21 - 0.5*tmp_23 - 0.5*tmp_25;
      real_t tmp_63 = tmp_25*tmp_62;
      real_t tmp_64 = tmp_20*tmp_32;
      real_t tmp_65 = tmp_31*tmp_64;
      real_t tmp_66 = -0.5*tmp_35 - 0.5*tmp_37 - 0.5*tmp_39;
      real_t tmp_67 = tmp_39*tmp_66;
      real_t tmp_68 = 2.0*Scalar_Variable_Coefficient_3D_1_0;
      real_t tmp_69 = 2.0*Scalar_Variable_Coefficient_3D_2_0;
      real_t tmp_70 = 2.0*Scalar_Variable_Coefficient_3D_3_0;
      real_t tmp_71 = tmp_54*(tmp_61*tmp_63 + tmp_61*tmp_65 + tmp_61*tmp_67) + tmp_56*(tmp_63*tmp_68 + tmp_65*tmp_68 + tmp_67*tmp_68) + tmp_58*(tmp_63*tmp_69 + tmp_65*tmp_69 + tmp_67*tmp_69) + tmp_60*(tmp_63*tmp_70 + tmp_65*tmp_70 + tmp_67*tmp_70);
      real_t tmp_72 = tmp_23*tmp_62;
      real_t tmp_73 = tmp_30*tmp_64;
      real_t tmp_74 = tmp_37*tmp_66;
      real_t tmp_75 = tmp_54*(tmp_61*tmp_72 + tmp_61*tmp_73 + tmp_61*tmp_74) + tmp_56*(tmp_68*tmp_72 + tmp_68*tmp_73 + tmp_68*tmp_74) + tmp_58*(tmp_69*tmp_72 + tmp_69*tmp_73 + tmp_69*tmp_74) + tmp_60*(tmp_70*tmp_72 + tmp_70*tmp_73 + tmp_70*tmp_74);
      real_t tmp_76 = tmp_21*tmp_62;
      real_t tmp_77 = tmp_28*tmp_64;
      real_t tmp_78 = tmp_35*tmp_66;
      real_t tmp_79 = tmp_54*(tmp_61*tmp_76 + tmp_61*tmp_77 + tmp_61*tmp_78) + tmp_56*(tmp_68*tmp_76 + tmp_68*tmp_77 + tmp_68*tmp_78) + tmp_58*(tmp_69*tmp_76 + tmp_69*tmp_77 + tmp_69*tmp_78) + tmp_60*(tmp_70*tmp_76 + tmp_70*tmp_77 + tmp_70*tmp_78);
      real_t tmp_80 = (tmp_24*tmp_24);
      real_t tmp_81 = 1.0 / (tmp_19*tmp_19);
      real_t tmp_82 = tmp_27*tmp_81;
      real_t tmp_83 = 2.0*tmp_81;
      real_t tmp_84 = (tmp_31*tmp_31)*tmp_83;
      real_t tmp_85 = (tmp_38*tmp_38);
      real_t tmp_86 = tmp_55*tmp_81;
      real_t tmp_87 = tmp_57*tmp_81;
      real_t tmp_88 = tmp_59*tmp_81;
      real_t tmp_89 = tmp_22*tmp_24;
      real_t tmp_90 = tmp_31*tmp_83;
      real_t tmp_91 = tmp_30*tmp_90;
      real_t tmp_92 = tmp_36*tmp_38;
      real_t tmp_93 = tmp_54*(Scalar_Variable_Coefficient_3D_0_0*tmp_91 + tmp_82*tmp_89 + tmp_82*tmp_92) + tmp_56*(Scalar_Variable_Coefficient_3D_1_0*tmp_91 + tmp_86*tmp_89 + tmp_86*tmp_92) + tmp_58*(Scalar_Variable_Coefficient_3D_2_0*tmp_91 + tmp_87*tmp_89 + tmp_87*tmp_92) + tmp_60*(Scalar_Variable_Coefficient_3D_3_0*tmp_91 + tmp_88*tmp_89 + tmp_88*tmp_92);
      real_t tmp_94 = tmp_24*tmp_6;
      real_t tmp_95 = tmp_28*tmp_90;
      real_t tmp_96 = tmp_34*tmp_38;
      real_t tmp_97 = tmp_54*(Scalar_Variable_Coefficient_3D_0_0*tmp_95 + tmp_82*tmp_94 + tmp_82*tmp_96) + tmp_56*(Scalar_Variable_Coefficient_3D_1_0*tmp_95 + tmp_86*tmp_94 + tmp_86*tmp_96) + tmp_58*(Scalar_Variable_Coefficient_3D_2_0*tmp_95 + tmp_87*tmp_94 + tmp_87*tmp_96) + tmp_60*(Scalar_Variable_Coefficient_3D_3_0*tmp_95 + tmp_88*tmp_94 + tmp_88*tmp_96);
      real_t tmp_98 = (tmp_22*tmp_22);
      real_t tmp_99 = (tmp_30*tmp_30)*tmp_83;
      real_t tmp_100 = (tmp_36*tmp_36);
      real_t tmp_101 = tmp_22*tmp_6;
      real_t tmp_102 = tmp_28*tmp_30*tmp_83;
      real_t tmp_103 = tmp_34*tmp_36;
      real_t tmp_104 = tmp_54*(Scalar_Variable_Coefficient_3D_0_0*tmp_102 + tmp_101*tmp_82 + tmp_103*tmp_82) + tmp_56*(Scalar_Variable_Coefficient_3D_1_0*tmp_102 + tmp_101*tmp_86 + tmp_103*tmp_86) + tmp_58*(Scalar_Variable_Coefficient_3D_2_0*tmp_102 + tmp_101*tmp_87 + tmp_103*tmp_87) + tmp_60*(Scalar_Variable_Coefficient_3D_3_0*tmp_102 + tmp_101*tmp_88 + tmp_103*tmp_88);
      real_t tmp_105 = (tmp_6*tmp_6);
      real_t tmp_106 = (tmp_28*tmp_28)*tmp_83;
      real_t tmp_107 = (tmp_34*tmp_34);
      real_t a_0_0 = tmp_54*(Scalar_Variable_Coefficient_3D_0_0*tmp_33 + tmp_26*tmp_27 + tmp_27*tmp_40) + tmp_56*(Scalar_Variable_Coefficient_3D_1_0*tmp_33 + tmp_26*tmp_55 + tmp_40*tmp_55) + tmp_58*(Scalar_Variable_Coefficient_3D_2_0*tmp_33 + tmp_26*tmp_57 + tmp_40*tmp_57) + tmp_60*(Scalar_Variable_Coefficient_3D_3_0*tmp_33 + tmp_26*tmp_59 + tmp_40*tmp_59);
      real_t a_0_1 = tmp_71;
      real_t a_0_2 = tmp_75;
      real_t a_0_3 = tmp_79;
      real_t a_1_0 = tmp_71;
      real_t a_1_1 = tmp_54*(Scalar_Variable_Coefficient_3D_0_0*tmp_84 + tmp_80*tmp_82 + tmp_82*tmp_85) + tmp_56*(Scalar_Variable_Coefficient_3D_1_0*tmp_84 + tmp_80*tmp_86 + tmp_85*tmp_86) + tmp_58*(Scalar_Variable_Coefficient_3D_2_0*tmp_84 + tmp_80*tmp_87 + tmp_85*tmp_87) + tmp_60*(Scalar_Variable_Coefficient_3D_3_0*tmp_84 + tmp_80*tmp_88 + tmp_85*tmp_88);
      real_t a_1_2 = tmp_93;
      real_t a_1_3 = tmp_97;
      real_t a_2_0 = tmp_75;
      real_t a_2_1 = tmp_93;
      real_t a_2_2 = tmp_54*(Scalar_Variable_Coefficient_3D_0_0*tmp_99 + tmp_100*tmp_82 + tmp_82*tmp_98) + tmp_56*(Scalar_Variable_Coefficient_3D_1_0*tmp_99 + tmp_100*tmp_86 + tmp_86*tmp_98) + tmp_58*(Scalar_Variable_Coefficient_3D_2_0*tmp_99 + tmp_100*tmp_87 + tmp_87*tmp_98) + tmp_60*(Scalar_Variable_Coefficient_3D_3_0*tmp_99 + tmp_100*tmp_88 + tmp_88*tmp_98);
      real_t a_2_3 = tmp_104;
      real_t a_3_0 = tmp_79;
      real_t a_3_1 = tmp_97;
      real_t a_3_2 = tmp_104;
      real_t a_3_3 = tmp_54*(Scalar_Variable_Coefficient_3D_0_0*tmp_106 + tmp_105*tmp_82 + tmp_107*tmp_82) + tmp_56*(Scalar_Variable_Coefficient_3D_1_0*tmp_106 + tmp_105*tmp_86 + tmp_107*tmp_86) + tmp_58*(Scalar_Variable_Coefficient_3D_2_0*tmp_106 + tmp_105*tmp_87 + tmp_107*tmp_87) + tmp_60*(Scalar_Variable_Coefficient_3D_3_0*tmp_106 + tmp_105*tmp_88 + tmp_107*tmp_88);
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

   void p1_epsilonvar_1_1_affine_q2::Scalar_Variable_Coefficient_2D( real_t in_0, real_t in_1, real_t * out_0 ) const
   {
      *out_0 = callback2D( Point3D( {in_0, in_1, 0} ) );
   }

   void p1_epsilonvar_1_1_affine_q2::Scalar_Variable_Coefficient_3D( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback3D( Point3D( {in_0, in_1, in_2} ) );
   }

} // namespace forms
} // namespace hyteg

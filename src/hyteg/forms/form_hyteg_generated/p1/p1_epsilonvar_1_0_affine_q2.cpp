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

#include "p1_epsilonvar_1_0_affine_q2.hpp"

namespace hyteg {
namespace forms {

   void p1_epsilonvar_1_0_affine_q2::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
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
      real_t tmp_0 = Scalar_Variable_Coefficient_2D_0_0*w_p_0;
      real_t tmp_1 = -p_affine_0_0;
      real_t tmp_2 = p_affine_1_0 + tmp_1;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = tmp_5 - (p_affine_1_1 + tmp_3)*(p_affine_2_0 + tmp_1);
      real_t tmp_7 = 1.0 / (tmp_6);
      real_t tmp_8 = 0.5*tmp_7;
      real_t tmp_9 = p_affine_0_0 - p_affine_2_0;
      real_t tmp_10 = -tmp_2*tmp_8 - tmp_8*tmp_9;
      real_t tmp_11 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_12 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_13 = tmp_12*(-tmp_11*tmp_8 - tmp_4*tmp_8);
      real_t tmp_14 = 4*tmp_10*tmp_13;
      real_t tmp_15 = Scalar_Variable_Coefficient_2D_1_0*w_p_1;
      real_t tmp_16 = Scalar_Variable_Coefficient_2D_2_0*w_p_2;
      real_t tmp_17 = tmp_0*tmp_9;
      real_t tmp_18 = 2.0*tmp_7;
      real_t tmp_19 = tmp_13*tmp_18;
      real_t tmp_20 = tmp_19*tmp_9;
      real_t tmp_21 = tmp_19*tmp_2;
      real_t tmp_22 = tmp_12*tmp_4;
      real_t tmp_23 = tmp_10*tmp_18;
      real_t tmp_24 = tmp_22*tmp_23;
      real_t tmp_25 = 1.0/(tmp_6*tmp_6);
      real_t tmp_26 = tmp_22*tmp_25;
      real_t tmp_27 = tmp_26*tmp_9;
      real_t tmp_28 = tmp_0*tmp_12;
      real_t tmp_29 = tmp_25*tmp_5;
      real_t tmp_30 = tmp_12*tmp_29;
      real_t tmp_31 = tmp_11*tmp_23;
      real_t tmp_32 = tmp_12*tmp_31;
      real_t tmp_33 = tmp_11*tmp_25;
      real_t tmp_34 = tmp_12*tmp_33;
      real_t tmp_35 = tmp_34*tmp_9;
      real_t tmp_36 = tmp_2*tmp_34;
      real_t a_0_0 = tmp_0*tmp_14 + tmp_14*tmp_15 + tmp_14*tmp_16;
      real_t a_0_1 = tmp_15*tmp_20 + tmp_16*tmp_20 + tmp_17*tmp_19;
      real_t a_0_2 = tmp_0*tmp_21 + tmp_15*tmp_21 + tmp_16*tmp_21;
      real_t a_1_0 = tmp_0*tmp_24 + tmp_15*tmp_24 + tmp_16*tmp_24;
      real_t a_1_1 = tmp_15*tmp_27 + tmp_16*tmp_27 + tmp_17*tmp_26;
      real_t a_1_2 = tmp_15*tmp_30 + tmp_16*tmp_30 + tmp_28*tmp_29;
      real_t a_2_0 = tmp_15*tmp_32 + tmp_16*tmp_32 + tmp_28*tmp_31;
      real_t a_2_1 = tmp_15*tmp_35 + tmp_16*tmp_35 + tmp_17*tmp_34;
      real_t a_2_2 = tmp_15*tmp_36 + tmp_16*tmp_36 + tmp_2*tmp_28*tmp_33;
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

   void p1_epsilonvar_1_0_affine_q2::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t tmp_0 = Scalar_Variable_Coefficient_3D_0_0*w_p_0;
      real_t tmp_1 = -p_affine_0_0;
      real_t tmp_2 = p_affine_2_0 + tmp_1;
      real_t tmp_3 = -p_affine_0_2;
      real_t tmp_4 = p_affine_1_2 + tmp_3;
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = p_affine_1_0 + tmp_1;
      real_t tmp_7 = p_affine_2_2 + tmp_3;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = tmp_5 - tmp_8;
      real_t tmp_10 = -p_affine_0_1;
      real_t tmp_11 = p_affine_2_1 + tmp_10;
      real_t tmp_12 = p_affine_3_2 + tmp_3;
      real_t tmp_13 = tmp_12*tmp_6;
      real_t tmp_14 = p_affine_3_1 + tmp_10;
      real_t tmp_15 = p_affine_1_1 + tmp_10;
      real_t tmp_16 = p_affine_3_0 + tmp_1;
      real_t tmp_17 = tmp_16*tmp_7;
      real_t tmp_18 = tmp_12*tmp_2;
      real_t tmp_19 = tmp_16*tmp_4;
      real_t tmp_20 = tmp_11*tmp_13 - tmp_11*tmp_19 + tmp_14*tmp_5 - tmp_14*tmp_8 + tmp_15*tmp_17 - tmp_15*tmp_18;
      real_t tmp_21 = 1.0 / (tmp_20);
      real_t tmp_22 = 0.5*tmp_21;
      real_t tmp_23 = tmp_13 - tmp_19;
      real_t tmp_24 = tmp_17 - tmp_18;
      real_t tmp_25 = -tmp_22*tmp_23 - tmp_22*tmp_24 - tmp_22*tmp_9;
      real_t tmp_26 = -tmp_11*tmp_4 + tmp_15*tmp_7;
      real_t tmp_27 = -tmp_12*tmp_15 + tmp_14*tmp_4;
      real_t tmp_28 = tmp_11*tmp_12 - tmp_14*tmp_7;
      real_t tmp_29 = p_affine_0_0*p_affine_1_1;
      real_t tmp_30 = p_affine_0_0*p_affine_1_2;
      real_t tmp_31 = p_affine_2_1*p_affine_3_2;
      real_t tmp_32 = p_affine_0_1*p_affine_1_0;
      real_t tmp_33 = p_affine_0_1*p_affine_1_2;
      real_t tmp_34 = p_affine_2_2*p_affine_3_0;
      real_t tmp_35 = p_affine_0_2*p_affine_1_0;
      real_t tmp_36 = p_affine_0_2*p_affine_1_1;
      real_t tmp_37 = p_affine_2_0*p_affine_3_1;
      real_t tmp_38 = p_affine_2_2*p_affine_3_1;
      real_t tmp_39 = p_affine_2_0*p_affine_3_2;
      real_t tmp_40 = p_affine_2_1*p_affine_3_0;
      real_t tmp_41 = std::abs(p_affine_0_0*tmp_31 - p_affine_0_0*tmp_38 + p_affine_0_1*tmp_34 - p_affine_0_1*tmp_39 + p_affine_0_2*tmp_37 - p_affine_0_2*tmp_40 - p_affine_1_0*tmp_31 + p_affine_1_0*tmp_38 - p_affine_1_1*tmp_34 + p_affine_1_1*tmp_39 - p_affine_1_2*tmp_37 + p_affine_1_2*tmp_40 + p_affine_2_0*tmp_33 - p_affine_2_0*tmp_36 - p_affine_2_1*tmp_30 + p_affine_2_1*tmp_35 + p_affine_2_2*tmp_29 - p_affine_2_2*tmp_32 - p_affine_3_0*tmp_33 + p_affine_3_0*tmp_36 + p_affine_3_1*tmp_30 - p_affine_3_1*tmp_35 - p_affine_3_2*tmp_29 + p_affine_3_2*tmp_32);
      real_t tmp_42 = tmp_41*(-tmp_22*tmp_26 - tmp_22*tmp_27 - tmp_22*tmp_28);
      real_t tmp_43 = 4*tmp_25*tmp_42;
      real_t tmp_44 = Scalar_Variable_Coefficient_3D_1_0*w_p_1;
      real_t tmp_45 = Scalar_Variable_Coefficient_3D_2_0*w_p_2;
      real_t tmp_46 = Scalar_Variable_Coefficient_3D_3_0*w_p_3;
      real_t tmp_47 = tmp_0*tmp_24;
      real_t tmp_48 = 2.0*tmp_21;
      real_t tmp_49 = tmp_42*tmp_48;
      real_t tmp_50 = tmp_24*tmp_49;
      real_t tmp_51 = tmp_23*tmp_49;
      real_t tmp_52 = tmp_49*tmp_9;
      real_t tmp_53 = tmp_28*tmp_41;
      real_t tmp_54 = tmp_0*tmp_53;
      real_t tmp_55 = tmp_25*tmp_48;
      real_t tmp_56 = tmp_53*tmp_55;
      real_t tmp_57 = 1.0/(tmp_20*tmp_20);
      real_t tmp_58 = tmp_53*tmp_57;
      real_t tmp_59 = tmp_24*tmp_58;
      real_t tmp_60 = tmp_54*tmp_57;
      real_t tmp_61 = tmp_23*tmp_58;
      real_t tmp_62 = tmp_58*tmp_9;
      real_t tmp_63 = tmp_27*tmp_41;
      real_t tmp_64 = tmp_55*tmp_63;
      real_t tmp_65 = tmp_57*tmp_63;
      real_t tmp_66 = tmp_24*tmp_65;
      real_t tmp_67 = tmp_23*tmp_65;
      real_t tmp_68 = tmp_65*tmp_9;
      real_t tmp_69 = tmp_26*tmp_41;
      real_t tmp_70 = tmp_55*tmp_69;
      real_t tmp_71 = tmp_57*tmp_69;
      real_t tmp_72 = tmp_24*tmp_71;
      real_t tmp_73 = tmp_23*tmp_71;
      real_t tmp_74 = tmp_71*tmp_9;
      real_t a_0_0 = tmp_0*tmp_43 + tmp_43*tmp_44 + tmp_43*tmp_45 + tmp_43*tmp_46;
      real_t a_0_1 = tmp_44*tmp_50 + tmp_45*tmp_50 + tmp_46*tmp_50 + tmp_47*tmp_49;
      real_t a_0_2 = tmp_0*tmp_51 + tmp_44*tmp_51 + tmp_45*tmp_51 + tmp_46*tmp_51;
      real_t a_0_3 = tmp_0*tmp_52 + tmp_44*tmp_52 + tmp_45*tmp_52 + tmp_46*tmp_52;
      real_t a_1_0 = tmp_44*tmp_56 + tmp_45*tmp_56 + tmp_46*tmp_56 + tmp_54*tmp_55;
      real_t a_1_1 = tmp_44*tmp_59 + tmp_45*tmp_59 + tmp_46*tmp_59 + tmp_47*tmp_58;
      real_t a_1_2 = tmp_23*tmp_60 + tmp_44*tmp_61 + tmp_45*tmp_61 + tmp_46*tmp_61;
      real_t a_1_3 = tmp_44*tmp_62 + tmp_45*tmp_62 + tmp_46*tmp_62 + tmp_60*tmp_9;
      real_t a_2_0 = tmp_0*tmp_64 + tmp_44*tmp_64 + tmp_45*tmp_64 + tmp_46*tmp_64;
      real_t a_2_1 = tmp_44*tmp_66 + tmp_45*tmp_66 + tmp_46*tmp_66 + tmp_47*tmp_65;
      real_t a_2_2 = tmp_0*tmp_67 + tmp_44*tmp_67 + tmp_45*tmp_67 + tmp_46*tmp_67;
      real_t a_2_3 = tmp_0*tmp_68 + tmp_44*tmp_68 + tmp_45*tmp_68 + tmp_46*tmp_68;
      real_t a_3_0 = tmp_0*tmp_70 + tmp_44*tmp_70 + tmp_45*tmp_70 + tmp_46*tmp_70;
      real_t a_3_1 = tmp_44*tmp_72 + tmp_45*tmp_72 + tmp_46*tmp_72 + tmp_47*tmp_71;
      real_t a_3_2 = tmp_0*tmp_73 + tmp_44*tmp_73 + tmp_45*tmp_73 + tmp_46*tmp_73;
      real_t a_3_3 = tmp_0*tmp_74 + tmp_44*tmp_74 + tmp_45*tmp_74 + tmp_46*tmp_74;
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

   void p1_epsilonvar_1_0_affine_q2::Scalar_Variable_Coefficient_2D( real_t in_0, real_t in_1, real_t * out_0 ) const
   {
      *out_0 = callback2D( Point3D( {in_0, in_1, 0} ) );
   }

   void p1_epsilonvar_1_0_affine_q2::Scalar_Variable_Coefficient_3D( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback3D( Point3D( {in_0, in_1, in_2} ) );
   }

} // namespace forms
} // namespace hyteg

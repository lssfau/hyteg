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

#include "p1_diffusion_affine_q1.hpp"

namespace hyteg {
namespace forms {

   void p1_diffusion_affine_q1::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
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
      real_t tmp_10 = tmp_3*tmp_5;
      real_t tmp_11 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_12 = tmp_11*tmp_5;
      real_t tmp_13 = -tmp_10 - tmp_12;
      real_t tmp_14 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_15 = tmp_14*(tmp_10*tmp_13 + tmp_8*tmp_9);
      real_t tmp_16 = tmp_14*(tmp_12*tmp_13 + tmp_6*tmp_9);
      real_t tmp_17 = 1.0 / (tmp_4*tmp_4);
      real_t tmp_18 = tmp_14*(tmp_1*tmp_17*tmp_7 + tmp_11*tmp_17*tmp_3);
      real_t a_0_0 = tmp_14*((tmp_13*tmp_13) + (tmp_9*tmp_9));
      real_t a_0_1 = tmp_15;
      real_t a_0_2 = tmp_16;
      real_t a_1_0 = tmp_15;
      real_t a_1_1 = tmp_14*(tmp_17*(tmp_3*tmp_3) + tmp_17*(tmp_7*tmp_7));
      real_t a_1_2 = tmp_18;
      real_t a_2_0 = tmp_16;
      real_t a_2_1 = tmp_18;
      real_t a_2_2 = tmp_14*((tmp_1*tmp_1)*tmp_17 + (tmp_11*tmp_11)*tmp_17);
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

   void p1_diffusion_affine_q1::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = tmp_1*tmp_4;
      real_t tmp_6 = tmp_4*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_7 = -tmp_5 - tmp_6;
      real_t tmp_8 = tmp_3*tmp_4;
      real_t tmp_9 = tmp_4*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_10 = -tmp_8 - tmp_9;
      real_t tmp_11 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = tmp_11*((tmp_10*tmp_10) + (tmp_7*tmp_7));
      real_t a_0_1 = tmp_11*(tmp_10*tmp_8 + tmp_6*tmp_7);
      real_t a_0_2 = tmp_11*(tmp_10*tmp_9 + tmp_5*tmp_7);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_diffusion_affine_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t tmp_53 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_43 - p_affine_0_0*tmp_50 + p_affine_0_1*tmp_46 - p_affine_0_1*tmp_51 + p_affine_0_2*tmp_49 - p_affine_0_2*tmp_52 - p_affine_1_0*tmp_43 + p_affine_1_0*tmp_50 - p_affine_1_1*tmp_46 + p_affine_1_1*tmp_51 - p_affine_1_2*tmp_49 + p_affine_1_2*tmp_52 + p_affine_2_0*tmp_45 - p_affine_2_0*tmp_48 - p_affine_2_1*tmp_42 + p_affine_2_1*tmp_47 + p_affine_2_2*tmp_41 - p_affine_2_2*tmp_44 - p_affine_3_0*tmp_45 + p_affine_3_0*tmp_48 + p_affine_3_1*tmp_42 - p_affine_3_1*tmp_47 - p_affine_3_2*tmp_41 + p_affine_3_2*tmp_44);
      real_t tmp_54 = tmp_53*(tmp_25*tmp_26 + tmp_32*tmp_33 + tmp_39*tmp_40);
      real_t tmp_55 = tmp_53*(tmp_23*tmp_26 + tmp_30*tmp_33 + tmp_37*tmp_40);
      real_t tmp_56 = tmp_53*(tmp_21*tmp_26 + tmp_28*tmp_33 + tmp_35*tmp_40);
      real_t tmp_57 = 1.0 / (tmp_19*tmp_19);
      real_t tmp_58 = tmp_24*tmp_57;
      real_t tmp_59 = tmp_31*tmp_57;
      real_t tmp_60 = tmp_38*tmp_57;
      real_t tmp_61 = tmp_53*(tmp_22*tmp_58 + tmp_29*tmp_59 + tmp_36*tmp_60);
      real_t tmp_62 = tmp_53*(tmp_27*tmp_59 + tmp_34*tmp_60 + tmp_58*tmp_8);
      real_t tmp_63 = tmp_53*(tmp_22*tmp_57*tmp_8 + tmp_27*tmp_29*tmp_57 + tmp_34*tmp_36*tmp_57);
      real_t a_0_0 = tmp_53*((tmp_26*tmp_26) + (tmp_33*tmp_33) + (tmp_40*tmp_40));
      real_t a_0_1 = tmp_54;
      real_t a_0_2 = tmp_55;
      real_t a_0_3 = tmp_56;
      real_t a_1_0 = tmp_54;
      real_t a_1_1 = tmp_53*((tmp_24*tmp_24)*tmp_57 + (tmp_31*tmp_31)*tmp_57 + (tmp_38*tmp_38)*tmp_57);
      real_t a_1_2 = tmp_61;
      real_t a_1_3 = tmp_62;
      real_t a_2_0 = tmp_55;
      real_t a_2_1 = tmp_61;
      real_t a_2_2 = tmp_53*((tmp_22*tmp_22)*tmp_57 + (tmp_29*tmp_29)*tmp_57 + (tmp_36*tmp_36)*tmp_57);
      real_t a_2_3 = tmp_63;
      real_t a_3_0 = tmp_56;
      real_t a_3_1 = tmp_62;
      real_t a_3_2 = tmp_63;
      real_t a_3_3 = tmp_53*((tmp_27*tmp_27)*tmp_57 + (tmp_34*tmp_34)*tmp_57 + tmp_57*(tmp_8*tmp_8));
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

   void p1_diffusion_affine_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t tmp_43 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_33 - p_affine_0_0*tmp_40 + p_affine_0_1*tmp_36 - p_affine_0_1*tmp_41 + p_affine_0_2*tmp_39 - p_affine_0_2*tmp_42 - p_affine_1_0*tmp_33 + p_affine_1_0*tmp_40 - p_affine_1_1*tmp_36 + p_affine_1_1*tmp_41 - p_affine_1_2*tmp_39 + p_affine_1_2*tmp_42 + p_affine_2_0*tmp_35 - p_affine_2_0*tmp_38 - p_affine_2_1*tmp_32 + p_affine_2_1*tmp_37 + p_affine_2_2*tmp_31 - p_affine_2_2*tmp_34 - p_affine_3_0*tmp_35 + p_affine_3_0*tmp_38 + p_affine_3_1*tmp_32 - p_affine_3_1*tmp_37 - p_affine_3_2*tmp_31 + p_affine_3_2*tmp_34);
      real_t a_0_0 = tmp_43*((tmp_22*tmp_22) + (tmp_26*tmp_26) + (tmp_30*tmp_30));
      real_t a_0_1 = tmp_43*(tmp_21*tmp_22 + tmp_25*tmp_26 + tmp_29*tmp_30);
      real_t a_0_2 = tmp_43*(tmp_20*tmp_22 + tmp_24*tmp_26 + tmp_28*tmp_30);
      real_t a_0_3 = tmp_43*(tmp_19*tmp_22 + tmp_23*tmp_26 + tmp_27*tmp_30);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

} // namespace forms
} // namespace hyteg

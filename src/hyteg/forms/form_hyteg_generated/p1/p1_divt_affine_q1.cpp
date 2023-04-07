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

#include "p1_divt_affine_q1.hpp"

namespace hyteg {
namespace forms {

   void p1_divt_0_affine_q1::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_0 = 0.33333333333333343;
      real_t tmp_1 = -p_affine_0_1;
      real_t tmp_2 = p_affine_2_1 + tmp_1;
      real_t tmp_3 = -p_affine_0_0;
      real_t tmp_4 = 1.0 / (tmp_2*(p_affine_1_0 + tmp_3) - (p_affine_1_1 + tmp_1)*(p_affine_2_0 + tmp_3));
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = tmp_4*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_7 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_8 = tmp_7*(tmp_5 + tmp_6);
      real_t tmp_9 = tmp_5*tmp_7;
      real_t tmp_10 = tmp_6*tmp_7;
      real_t a_0_0 = tmp_0*tmp_8;
      real_t a_0_1 = 0.33333333333333331*tmp_8;
      real_t a_0_2 = 0.33333333333333331*tmp_8;
      real_t a_1_0 = -tmp_0*tmp_9;
      real_t a_1_1 = -0.33333333333333331*tmp_9;
      real_t a_1_2 = -0.33333333333333331*tmp_9;
      real_t a_2_0 = -tmp_0*tmp_10;
      real_t a_2_1 = -0.33333333333333331*tmp_10;
      real_t a_2_2 = -0.33333333333333331*tmp_10;
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

   void p1_divt_0_affine_q1::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_2_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_0;
      real_t tmp_3 = 1.0 / (tmp_1*(p_affine_1_0 + tmp_2) - (p_affine_1_1 + tmp_0)*(p_affine_2_0 + tmp_2));
      real_t tmp_4 = 0.5*(tmp_1*tmp_3 + tmp_3*(p_affine_0_1 - p_affine_1_1))*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = 0.33333333333333343*tmp_4;
      real_t a_0_1 = 0.33333333333333331*tmp_4;
      real_t a_0_2 = 0.33333333333333331*tmp_4;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_divt_0_affine_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t tmp_0 = 0.25;
      real_t tmp_1 = -p_affine_0_1;
      real_t tmp_2 = p_affine_1_1 + tmp_1;
      real_t tmp_3 = -p_affine_0_2;
      real_t tmp_4 = p_affine_2_2 + tmp_3;
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = p_affine_2_1 + tmp_1;
      real_t tmp_7 = p_affine_1_2 + tmp_3;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = -p_affine_0_0;
      real_t tmp_10 = p_affine_1_0 + tmp_9;
      real_t tmp_11 = p_affine_3_2 + tmp_3;
      real_t tmp_12 = tmp_11*tmp_6;
      real_t tmp_13 = p_affine_2_0 + tmp_9;
      real_t tmp_14 = p_affine_3_1 + tmp_1;
      real_t tmp_15 = tmp_14*tmp_7;
      real_t tmp_16 = p_affine_3_0 + tmp_9;
      real_t tmp_17 = tmp_14*tmp_4;
      real_t tmp_18 = tmp_11*tmp_2;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_12 - tmp_10*tmp_17 + tmp_13*tmp_15 - tmp_13*tmp_18 + tmp_16*tmp_5 - tmp_16*tmp_8);
      real_t tmp_20 = tmp_19*(tmp_5 - tmp_8);
      real_t tmp_21 = tmp_19*(tmp_15 - tmp_18);
      real_t tmp_22 = tmp_19*(tmp_12 - tmp_17);
      real_t tmp_23 = p_affine_0_0*p_affine_1_1;
      real_t tmp_24 = p_affine_0_0*p_affine_1_2;
      real_t tmp_25 = p_affine_2_1*p_affine_3_2;
      real_t tmp_26 = p_affine_0_1*p_affine_1_0;
      real_t tmp_27 = p_affine_0_1*p_affine_1_2;
      real_t tmp_28 = p_affine_2_2*p_affine_3_0;
      real_t tmp_29 = p_affine_0_2*p_affine_1_0;
      real_t tmp_30 = p_affine_0_2*p_affine_1_1;
      real_t tmp_31 = p_affine_2_0*p_affine_3_1;
      real_t tmp_32 = p_affine_2_2*p_affine_3_1;
      real_t tmp_33 = p_affine_2_0*p_affine_3_2;
      real_t tmp_34 = p_affine_2_1*p_affine_3_0;
      real_t tmp_35 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_25 - p_affine_0_0*tmp_32 + p_affine_0_1*tmp_28 - p_affine_0_1*tmp_33 + p_affine_0_2*tmp_31 - p_affine_0_2*tmp_34 - p_affine_1_0*tmp_25 + p_affine_1_0*tmp_32 - p_affine_1_1*tmp_28 + p_affine_1_1*tmp_33 - p_affine_1_2*tmp_31 + p_affine_1_2*tmp_34 + p_affine_2_0*tmp_27 - p_affine_2_0*tmp_30 - p_affine_2_1*tmp_24 + p_affine_2_1*tmp_29 + p_affine_2_2*tmp_23 - p_affine_2_2*tmp_26 - p_affine_3_0*tmp_27 + p_affine_3_0*tmp_30 + p_affine_3_1*tmp_24 - p_affine_3_1*tmp_29 - p_affine_3_2*tmp_23 + p_affine_3_2*tmp_26);
      real_t tmp_36 = tmp_35*(tmp_20 + tmp_21 + tmp_22);
      real_t tmp_37 = tmp_22*tmp_35;
      real_t tmp_38 = tmp_21*tmp_35;
      real_t tmp_39 = tmp_20*tmp_35;
      real_t a_0_0 = tmp_0*tmp_36;
      real_t a_0_1 = 0.25*tmp_36;
      real_t a_0_2 = 0.25*tmp_36;
      real_t a_0_3 = 0.25*tmp_36;
      real_t a_1_0 = -tmp_0*tmp_37;
      real_t a_1_1 = -0.25*tmp_37;
      real_t a_1_2 = -0.25*tmp_37;
      real_t a_1_3 = -0.25*tmp_37;
      real_t a_2_0 = -tmp_0*tmp_38;
      real_t a_2_1 = -0.25*tmp_38;
      real_t a_2_2 = -0.25*tmp_38;
      real_t a_2_3 = -0.25*tmp_38;
      real_t a_3_0 = -tmp_0*tmp_39;
      real_t a_3_1 = -0.25*tmp_39;
      real_t a_3_2 = -0.25*tmp_39;
      real_t a_3_3 = -0.25*tmp_39;
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

   void p1_divt_0_affine_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t tmp_0 = -p_affine_0_1;
      real_t tmp_1 = p_affine_1_1 + tmp_0;
      real_t tmp_2 = -p_affine_0_2;
      real_t tmp_3 = p_affine_2_2 + tmp_2;
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = p_affine_2_1 + tmp_0;
      real_t tmp_6 = p_affine_1_2 + tmp_2;
      real_t tmp_7 = tmp_5*tmp_6;
      real_t tmp_8 = -p_affine_0_0;
      real_t tmp_9 = p_affine_1_0 + tmp_8;
      real_t tmp_10 = p_affine_3_2 + tmp_2;
      real_t tmp_11 = tmp_10*tmp_5;
      real_t tmp_12 = p_affine_2_0 + tmp_8;
      real_t tmp_13 = p_affine_3_1 + tmp_0;
      real_t tmp_14 = tmp_13*tmp_6;
      real_t tmp_15 = p_affine_3_0 + tmp_8;
      real_t tmp_16 = tmp_13*tmp_3;
      real_t tmp_17 = tmp_1*tmp_10;
      real_t tmp_18 = 1.0 / (tmp_11*tmp_9 + tmp_12*tmp_14 - tmp_12*tmp_17 + tmp_15*tmp_4 - tmp_15*tmp_7 - tmp_16*tmp_9);
      real_t tmp_19 = p_affine_0_0*p_affine_1_1;
      real_t tmp_20 = p_affine_0_0*p_affine_1_2;
      real_t tmp_21 = p_affine_2_1*p_affine_3_2;
      real_t tmp_22 = p_affine_0_1*p_affine_1_0;
      real_t tmp_23 = p_affine_0_1*p_affine_1_2;
      real_t tmp_24 = p_affine_2_2*p_affine_3_0;
      real_t tmp_25 = p_affine_0_2*p_affine_1_0;
      real_t tmp_26 = p_affine_0_2*p_affine_1_1;
      real_t tmp_27 = p_affine_2_0*p_affine_3_1;
      real_t tmp_28 = p_affine_2_2*p_affine_3_1;
      real_t tmp_29 = p_affine_2_0*p_affine_3_2;
      real_t tmp_30 = p_affine_2_1*p_affine_3_0;
      real_t tmp_31 = 0.16666666666666663*(tmp_18*(tmp_11 - tmp_16) + tmp_18*(tmp_14 - tmp_17) + tmp_18*(tmp_4 - tmp_7))*std::abs(p_affine_0_0*tmp_21 - p_affine_0_0*tmp_28 + p_affine_0_1*tmp_24 - p_affine_0_1*tmp_29 + p_affine_0_2*tmp_27 - p_affine_0_2*tmp_30 - p_affine_1_0*tmp_21 + p_affine_1_0*tmp_28 - p_affine_1_1*tmp_24 + p_affine_1_1*tmp_29 - p_affine_1_2*tmp_27 + p_affine_1_2*tmp_30 + p_affine_2_0*tmp_23 - p_affine_2_0*tmp_26 - p_affine_2_1*tmp_20 + p_affine_2_1*tmp_25 + p_affine_2_2*tmp_19 - p_affine_2_2*tmp_22 - p_affine_3_0*tmp_23 + p_affine_3_0*tmp_26 + p_affine_3_1*tmp_20 - p_affine_3_1*tmp_25 - p_affine_3_2*tmp_19 + p_affine_3_2*tmp_22);
      real_t a_0_0 = 0.25*tmp_31;
      real_t a_0_1 = 0.25*tmp_31;
      real_t a_0_2 = 0.25*tmp_31;
      real_t a_0_3 = 0.25*tmp_31;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_divt_1_affine_q1::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_0 = 0.33333333333333343;
      real_t tmp_1 = -p_affine_0_0;
      real_t tmp_2 = p_affine_1_0 + tmp_1;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = 1.0 / (tmp_2*(p_affine_2_1 + tmp_3) - (p_affine_1_1 + tmp_3)*(p_affine_2_0 + tmp_1));
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = tmp_4*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_7 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_8 = tmp_7*(tmp_5 + tmp_6);
      real_t tmp_9 = tmp_6*tmp_7;
      real_t tmp_10 = tmp_5*tmp_7;
      real_t a_0_0 = tmp_0*tmp_8;
      real_t a_0_1 = 0.33333333333333331*tmp_8;
      real_t a_0_2 = 0.33333333333333331*tmp_8;
      real_t a_1_0 = -tmp_0*tmp_9;
      real_t a_1_1 = -0.33333333333333331*tmp_9;
      real_t a_1_2 = -0.33333333333333331*tmp_9;
      real_t a_2_0 = -tmp_0*tmp_10;
      real_t a_2_1 = -0.33333333333333331*tmp_10;
      real_t a_2_2 = -0.33333333333333331*tmp_10;
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

   void p1_divt_1_affine_q1::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
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
      real_t tmp_3 = 1.0 / (tmp_1*(p_affine_2_1 + tmp_2) - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_4 = 0.5*(tmp_1*tmp_3 + tmp_3*(p_affine_0_0 - p_affine_2_0))*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = 0.33333333333333343*tmp_4;
      real_t a_0_1 = 0.33333333333333331*tmp_4;
      real_t a_0_2 = 0.33333333333333331*tmp_4;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_divt_1_affine_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t tmp_0 = 0.25;
      real_t tmp_1 = -p_affine_0_0;
      real_t tmp_2 = p_affine_2_0 + tmp_1;
      real_t tmp_3 = -p_affine_0_2;
      real_t tmp_4 = p_affine_1_2 + tmp_3;
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = p_affine_1_0 + tmp_1;
      real_t tmp_7 = p_affine_2_2 + tmp_3;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = -p_affine_0_1;
      real_t tmp_10 = p_affine_2_1 + tmp_9;
      real_t tmp_11 = p_affine_3_2 + tmp_3;
      real_t tmp_12 = tmp_11*tmp_6;
      real_t tmp_13 = p_affine_3_1 + tmp_9;
      real_t tmp_14 = p_affine_1_1 + tmp_9;
      real_t tmp_15 = p_affine_3_0 + tmp_1;
      real_t tmp_16 = tmp_15*tmp_7;
      real_t tmp_17 = tmp_11*tmp_2;
      real_t tmp_18 = tmp_15*tmp_4;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_12 - tmp_10*tmp_18 + tmp_13*tmp_5 - tmp_13*tmp_8 + tmp_14*tmp_16 - tmp_14*tmp_17);
      real_t tmp_20 = tmp_19*(tmp_5 - tmp_8);
      real_t tmp_21 = tmp_19*(tmp_12 - tmp_18);
      real_t tmp_22 = tmp_19*(tmp_16 - tmp_17);
      real_t tmp_23 = p_affine_0_0*p_affine_1_1;
      real_t tmp_24 = p_affine_0_0*p_affine_1_2;
      real_t tmp_25 = p_affine_2_1*p_affine_3_2;
      real_t tmp_26 = p_affine_0_1*p_affine_1_0;
      real_t tmp_27 = p_affine_0_1*p_affine_1_2;
      real_t tmp_28 = p_affine_2_2*p_affine_3_0;
      real_t tmp_29 = p_affine_0_2*p_affine_1_0;
      real_t tmp_30 = p_affine_0_2*p_affine_1_1;
      real_t tmp_31 = p_affine_2_0*p_affine_3_1;
      real_t tmp_32 = p_affine_2_2*p_affine_3_1;
      real_t tmp_33 = p_affine_2_0*p_affine_3_2;
      real_t tmp_34 = p_affine_2_1*p_affine_3_0;
      real_t tmp_35 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_25 - p_affine_0_0*tmp_32 + p_affine_0_1*tmp_28 - p_affine_0_1*tmp_33 + p_affine_0_2*tmp_31 - p_affine_0_2*tmp_34 - p_affine_1_0*tmp_25 + p_affine_1_0*tmp_32 - p_affine_1_1*tmp_28 + p_affine_1_1*tmp_33 - p_affine_1_2*tmp_31 + p_affine_1_2*tmp_34 + p_affine_2_0*tmp_27 - p_affine_2_0*tmp_30 - p_affine_2_1*tmp_24 + p_affine_2_1*tmp_29 + p_affine_2_2*tmp_23 - p_affine_2_2*tmp_26 - p_affine_3_0*tmp_27 + p_affine_3_0*tmp_30 + p_affine_3_1*tmp_24 - p_affine_3_1*tmp_29 - p_affine_3_2*tmp_23 + p_affine_3_2*tmp_26);
      real_t tmp_36 = tmp_35*(tmp_20 + tmp_21 + tmp_22);
      real_t tmp_37 = tmp_22*tmp_35;
      real_t tmp_38 = tmp_21*tmp_35;
      real_t tmp_39 = tmp_20*tmp_35;
      real_t a_0_0 = tmp_0*tmp_36;
      real_t a_0_1 = 0.25*tmp_36;
      real_t a_0_2 = 0.25*tmp_36;
      real_t a_0_3 = 0.25*tmp_36;
      real_t a_1_0 = -tmp_0*tmp_37;
      real_t a_1_1 = -0.25*tmp_37;
      real_t a_1_2 = -0.25*tmp_37;
      real_t a_1_3 = -0.25*tmp_37;
      real_t a_2_0 = -tmp_0*tmp_38;
      real_t a_2_1 = -0.25*tmp_38;
      real_t a_2_2 = -0.25*tmp_38;
      real_t a_2_3 = -0.25*tmp_38;
      real_t a_3_0 = -tmp_0*tmp_39;
      real_t a_3_1 = -0.25*tmp_39;
      real_t a_3_2 = -0.25*tmp_39;
      real_t a_3_3 = -0.25*tmp_39;
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

   void p1_divt_1_affine_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t tmp_19 = p_affine_0_0*p_affine_1_1;
      real_t tmp_20 = p_affine_0_0*p_affine_1_2;
      real_t tmp_21 = p_affine_2_1*p_affine_3_2;
      real_t tmp_22 = p_affine_0_1*p_affine_1_0;
      real_t tmp_23 = p_affine_0_1*p_affine_1_2;
      real_t tmp_24 = p_affine_2_2*p_affine_3_0;
      real_t tmp_25 = p_affine_0_2*p_affine_1_0;
      real_t tmp_26 = p_affine_0_2*p_affine_1_1;
      real_t tmp_27 = p_affine_2_0*p_affine_3_1;
      real_t tmp_28 = p_affine_2_2*p_affine_3_1;
      real_t tmp_29 = p_affine_2_0*p_affine_3_2;
      real_t tmp_30 = p_affine_2_1*p_affine_3_0;
      real_t tmp_31 = 0.16666666666666663*(tmp_18*(tmp_11 - tmp_17) + tmp_18*(tmp_15 - tmp_16) + tmp_18*(tmp_4 - tmp_7))*std::abs(p_affine_0_0*tmp_21 - p_affine_0_0*tmp_28 + p_affine_0_1*tmp_24 - p_affine_0_1*tmp_29 + p_affine_0_2*tmp_27 - p_affine_0_2*tmp_30 - p_affine_1_0*tmp_21 + p_affine_1_0*tmp_28 - p_affine_1_1*tmp_24 + p_affine_1_1*tmp_29 - p_affine_1_2*tmp_27 + p_affine_1_2*tmp_30 + p_affine_2_0*tmp_23 - p_affine_2_0*tmp_26 - p_affine_2_1*tmp_20 + p_affine_2_1*tmp_25 + p_affine_2_2*tmp_19 - p_affine_2_2*tmp_22 - p_affine_3_0*tmp_23 + p_affine_3_0*tmp_26 + p_affine_3_1*tmp_20 - p_affine_3_1*tmp_25 - p_affine_3_2*tmp_19 + p_affine_3_2*tmp_22);
      real_t a_0_0 = 0.25*tmp_31;
      real_t a_0_1 = 0.25*tmp_31;
      real_t a_0_2 = 0.25*tmp_31;
      real_t a_0_3 = 0.25*tmp_31;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_divt_2_affine_q1::integrateAll( const std::array< Point3D, 3 >& , Matrix< real_t, 3, 3 >&  ) const
   {
      
   }

   void p1_divt_2_affine_q1::integrateRow0( const std::array< Point3D, 3 >& , Matrix< real_t, 1, 3 >&  ) const
   {
      
   }

   void p1_divt_2_affine_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t tmp_0 = 0.25;
      real_t tmp_1 = -p_affine_0_0;
      real_t tmp_2 = p_affine_1_0 + tmp_1;
      real_t tmp_3 = -p_affine_0_1;
      real_t tmp_4 = p_affine_2_1 + tmp_3;
      real_t tmp_5 = tmp_2*tmp_4;
      real_t tmp_6 = p_affine_2_0 + tmp_1;
      real_t tmp_7 = p_affine_1_1 + tmp_3;
      real_t tmp_8 = tmp_6*tmp_7;
      real_t tmp_9 = -p_affine_0_2;
      real_t tmp_10 = p_affine_3_2 + tmp_9;
      real_t tmp_11 = p_affine_1_2 + tmp_9;
      real_t tmp_12 = p_affine_3_1 + tmp_3;
      real_t tmp_13 = tmp_12*tmp_6;
      real_t tmp_14 = p_affine_2_2 + tmp_9;
      real_t tmp_15 = p_affine_3_0 + tmp_1;
      real_t tmp_16 = tmp_15*tmp_7;
      real_t tmp_17 = tmp_12*tmp_2;
      real_t tmp_18 = tmp_15*tmp_4;
      real_t tmp_19 = 1.0 / (tmp_10*tmp_5 - tmp_10*tmp_8 + tmp_11*tmp_13 - tmp_11*tmp_18 + tmp_14*tmp_16 - tmp_14*tmp_17);
      real_t tmp_20 = tmp_19*(tmp_5 - tmp_8);
      real_t tmp_21 = tmp_19*(tmp_16 - tmp_17);
      real_t tmp_22 = tmp_19*(tmp_13 - tmp_18);
      real_t tmp_23 = p_affine_0_0*p_affine_1_1;
      real_t tmp_24 = p_affine_0_0*p_affine_1_2;
      real_t tmp_25 = p_affine_2_1*p_affine_3_2;
      real_t tmp_26 = p_affine_0_1*p_affine_1_0;
      real_t tmp_27 = p_affine_0_1*p_affine_1_2;
      real_t tmp_28 = p_affine_2_2*p_affine_3_0;
      real_t tmp_29 = p_affine_0_2*p_affine_1_0;
      real_t tmp_30 = p_affine_0_2*p_affine_1_1;
      real_t tmp_31 = p_affine_2_0*p_affine_3_1;
      real_t tmp_32 = p_affine_2_2*p_affine_3_1;
      real_t tmp_33 = p_affine_2_0*p_affine_3_2;
      real_t tmp_34 = p_affine_2_1*p_affine_3_0;
      real_t tmp_35 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_25 - p_affine_0_0*tmp_32 + p_affine_0_1*tmp_28 - p_affine_0_1*tmp_33 + p_affine_0_2*tmp_31 - p_affine_0_2*tmp_34 - p_affine_1_0*tmp_25 + p_affine_1_0*tmp_32 - p_affine_1_1*tmp_28 + p_affine_1_1*tmp_33 - p_affine_1_2*tmp_31 + p_affine_1_2*tmp_34 + p_affine_2_0*tmp_27 - p_affine_2_0*tmp_30 - p_affine_2_1*tmp_24 + p_affine_2_1*tmp_29 + p_affine_2_2*tmp_23 - p_affine_2_2*tmp_26 - p_affine_3_0*tmp_27 + p_affine_3_0*tmp_30 + p_affine_3_1*tmp_24 - p_affine_3_1*tmp_29 - p_affine_3_2*tmp_23 + p_affine_3_2*tmp_26);
      real_t tmp_36 = tmp_35*(tmp_20 + tmp_21 + tmp_22);
      real_t tmp_37 = tmp_22*tmp_35;
      real_t tmp_38 = tmp_21*tmp_35;
      real_t tmp_39 = tmp_20*tmp_35;
      real_t a_0_0 = tmp_0*tmp_36;
      real_t a_0_1 = 0.25*tmp_36;
      real_t a_0_2 = 0.25*tmp_36;
      real_t a_0_3 = 0.25*tmp_36;
      real_t a_1_0 = -tmp_0*tmp_37;
      real_t a_1_1 = -0.25*tmp_37;
      real_t a_1_2 = -0.25*tmp_37;
      real_t a_1_3 = -0.25*tmp_37;
      real_t a_2_0 = -tmp_0*tmp_38;
      real_t a_2_1 = -0.25*tmp_38;
      real_t a_2_2 = -0.25*tmp_38;
      real_t a_2_3 = -0.25*tmp_38;
      real_t a_3_0 = -tmp_0*tmp_39;
      real_t a_3_1 = -0.25*tmp_39;
      real_t a_3_2 = -0.25*tmp_39;
      real_t a_3_3 = -0.25*tmp_39;
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

   void p1_divt_2_affine_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t tmp_19 = p_affine_0_0*p_affine_1_1;
      real_t tmp_20 = p_affine_0_0*p_affine_1_2;
      real_t tmp_21 = p_affine_2_1*p_affine_3_2;
      real_t tmp_22 = p_affine_0_1*p_affine_1_0;
      real_t tmp_23 = p_affine_0_1*p_affine_1_2;
      real_t tmp_24 = p_affine_2_2*p_affine_3_0;
      real_t tmp_25 = p_affine_0_2*p_affine_1_0;
      real_t tmp_26 = p_affine_0_2*p_affine_1_1;
      real_t tmp_27 = p_affine_2_0*p_affine_3_1;
      real_t tmp_28 = p_affine_2_2*p_affine_3_1;
      real_t tmp_29 = p_affine_2_0*p_affine_3_2;
      real_t tmp_30 = p_affine_2_1*p_affine_3_0;
      real_t tmp_31 = 0.16666666666666663*(tmp_18*(tmp_12 - tmp_17) + tmp_18*(tmp_15 - tmp_16) + tmp_18*(tmp_4 - tmp_7))*std::abs(p_affine_0_0*tmp_21 - p_affine_0_0*tmp_28 + p_affine_0_1*tmp_24 - p_affine_0_1*tmp_29 + p_affine_0_2*tmp_27 - p_affine_0_2*tmp_30 - p_affine_1_0*tmp_21 + p_affine_1_0*tmp_28 - p_affine_1_1*tmp_24 + p_affine_1_1*tmp_29 - p_affine_1_2*tmp_27 + p_affine_1_2*tmp_30 + p_affine_2_0*tmp_23 - p_affine_2_0*tmp_26 - p_affine_2_1*tmp_20 + p_affine_2_1*tmp_25 + p_affine_2_2*tmp_19 - p_affine_2_2*tmp_22 - p_affine_3_0*tmp_23 + p_affine_3_0*tmp_26 + p_affine_3_1*tmp_20 - p_affine_3_1*tmp_25 - p_affine_3_2*tmp_19 + p_affine_3_2*tmp_22);
      real_t a_0_0 = 0.25*tmp_31;
      real_t a_0_1 = 0.25*tmp_31;
      real_t a_0_2 = 0.25*tmp_31;
      real_t a_0_3 = 0.25*tmp_31;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

} // namespace forms
} // namespace hyteg

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

#include "p1_to_p0_div_affine_q0.hpp"

namespace hyteg {
namespace forms {

   void p1_to_p0_div_0_affine_q0::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
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
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = tmp_3*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_6 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = tmp_6*(tmp_4 + tmp_5);
      real_t a_0_1 = -tmp_4*tmp_6;
      real_t a_0_2 = -tmp_5*tmp_6;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_to_p0_div_0_affine_q0::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
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
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = tmp_3*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_6 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = tmp_6*(tmp_4 + tmp_5);
      real_t a_0_1 = -tmp_4*tmp_6;
      real_t a_0_2 = -tmp_5*tmp_6;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_to_p0_div_0_affine_q0::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t tmp_19 = tmp_18*(tmp_4 - tmp_7);
      real_t tmp_20 = tmp_18*(tmp_14 - tmp_17);
      real_t tmp_21 = tmp_18*(tmp_11 - tmp_16);
      real_t tmp_22 = p_affine_0_0*p_affine_1_1;
      real_t tmp_23 = p_affine_0_0*p_affine_1_2;
      real_t tmp_24 = p_affine_2_1*p_affine_3_2;
      real_t tmp_25 = p_affine_0_1*p_affine_1_0;
      real_t tmp_26 = p_affine_0_1*p_affine_1_2;
      real_t tmp_27 = p_affine_2_2*p_affine_3_0;
      real_t tmp_28 = p_affine_0_2*p_affine_1_0;
      real_t tmp_29 = p_affine_0_2*p_affine_1_1;
      real_t tmp_30 = p_affine_2_0*p_affine_3_1;
      real_t tmp_31 = p_affine_2_2*p_affine_3_1;
      real_t tmp_32 = p_affine_2_0*p_affine_3_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_0;
      real_t tmp_34 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_24 - p_affine_0_0*tmp_31 + p_affine_0_1*tmp_27 - p_affine_0_1*tmp_32 + p_affine_0_2*tmp_30 - p_affine_0_2*tmp_33 - p_affine_1_0*tmp_24 + p_affine_1_0*tmp_31 - p_affine_1_1*tmp_27 + p_affine_1_1*tmp_32 - p_affine_1_2*tmp_30 + p_affine_1_2*tmp_33 + p_affine_2_0*tmp_26 - p_affine_2_0*tmp_29 - p_affine_2_1*tmp_23 + p_affine_2_1*tmp_28 + p_affine_2_2*tmp_22 - p_affine_2_2*tmp_25 - p_affine_3_0*tmp_26 + p_affine_3_0*tmp_29 + p_affine_3_1*tmp_23 - p_affine_3_1*tmp_28 - p_affine_3_2*tmp_22 + p_affine_3_2*tmp_25);
      real_t a_0_0 = tmp_34*(tmp_19 + tmp_20 + tmp_21);
      real_t a_0_1 = -tmp_21*tmp_34;
      real_t a_0_2 = -tmp_20*tmp_34;
      real_t a_0_3 = -tmp_19*tmp_34;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_to_p0_div_0_affine_q0::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t tmp_19 = tmp_18*(tmp_4 - tmp_7);
      real_t tmp_20 = tmp_18*(tmp_14 - tmp_17);
      real_t tmp_21 = tmp_18*(tmp_11 - tmp_16);
      real_t tmp_22 = p_affine_0_0*p_affine_1_1;
      real_t tmp_23 = p_affine_0_0*p_affine_1_2;
      real_t tmp_24 = p_affine_2_1*p_affine_3_2;
      real_t tmp_25 = p_affine_0_1*p_affine_1_0;
      real_t tmp_26 = p_affine_0_1*p_affine_1_2;
      real_t tmp_27 = p_affine_2_2*p_affine_3_0;
      real_t tmp_28 = p_affine_0_2*p_affine_1_0;
      real_t tmp_29 = p_affine_0_2*p_affine_1_1;
      real_t tmp_30 = p_affine_2_0*p_affine_3_1;
      real_t tmp_31 = p_affine_2_2*p_affine_3_1;
      real_t tmp_32 = p_affine_2_0*p_affine_3_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_0;
      real_t tmp_34 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_24 - p_affine_0_0*tmp_31 + p_affine_0_1*tmp_27 - p_affine_0_1*tmp_32 + p_affine_0_2*tmp_30 - p_affine_0_2*tmp_33 - p_affine_1_0*tmp_24 + p_affine_1_0*tmp_31 - p_affine_1_1*tmp_27 + p_affine_1_1*tmp_32 - p_affine_1_2*tmp_30 + p_affine_1_2*tmp_33 + p_affine_2_0*tmp_26 - p_affine_2_0*tmp_29 - p_affine_2_1*tmp_23 + p_affine_2_1*tmp_28 + p_affine_2_2*tmp_22 - p_affine_2_2*tmp_25 - p_affine_3_0*tmp_26 + p_affine_3_0*tmp_29 + p_affine_3_1*tmp_23 - p_affine_3_1*tmp_28 - p_affine_3_2*tmp_22 + p_affine_3_2*tmp_25);
      real_t a_0_0 = tmp_34*(tmp_19 + tmp_20 + tmp_21);
      real_t a_0_1 = -tmp_21*tmp_34;
      real_t a_0_2 = -tmp_20*tmp_34;
      real_t a_0_3 = -tmp_19*tmp_34;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_to_p0_div_1_affine_q0::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
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
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = tmp_3*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_6 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = tmp_6*(tmp_4 + tmp_5);
      real_t a_0_1 = -tmp_5*tmp_6;
      real_t a_0_2 = -tmp_4*tmp_6;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_to_p0_div_1_affine_q0::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
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
      real_t tmp_4 = tmp_1*tmp_3;
      real_t tmp_5 = tmp_3*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_6 = 0.5*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = tmp_6*(tmp_4 + tmp_5);
      real_t a_0_1 = -tmp_5*tmp_6;
      real_t a_0_2 = -tmp_4*tmp_6;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_to_p0_div_1_affine_q0::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t tmp_19 = tmp_18*(tmp_4 - tmp_7);
      real_t tmp_20 = tmp_18*(tmp_11 - tmp_17);
      real_t tmp_21 = tmp_18*(tmp_15 - tmp_16);
      real_t tmp_22 = p_affine_0_0*p_affine_1_1;
      real_t tmp_23 = p_affine_0_0*p_affine_1_2;
      real_t tmp_24 = p_affine_2_1*p_affine_3_2;
      real_t tmp_25 = p_affine_0_1*p_affine_1_0;
      real_t tmp_26 = p_affine_0_1*p_affine_1_2;
      real_t tmp_27 = p_affine_2_2*p_affine_3_0;
      real_t tmp_28 = p_affine_0_2*p_affine_1_0;
      real_t tmp_29 = p_affine_0_2*p_affine_1_1;
      real_t tmp_30 = p_affine_2_0*p_affine_3_1;
      real_t tmp_31 = p_affine_2_2*p_affine_3_1;
      real_t tmp_32 = p_affine_2_0*p_affine_3_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_0;
      real_t tmp_34 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_24 - p_affine_0_0*tmp_31 + p_affine_0_1*tmp_27 - p_affine_0_1*tmp_32 + p_affine_0_2*tmp_30 - p_affine_0_2*tmp_33 - p_affine_1_0*tmp_24 + p_affine_1_0*tmp_31 - p_affine_1_1*tmp_27 + p_affine_1_1*tmp_32 - p_affine_1_2*tmp_30 + p_affine_1_2*tmp_33 + p_affine_2_0*tmp_26 - p_affine_2_0*tmp_29 - p_affine_2_1*tmp_23 + p_affine_2_1*tmp_28 + p_affine_2_2*tmp_22 - p_affine_2_2*tmp_25 - p_affine_3_0*tmp_26 + p_affine_3_0*tmp_29 + p_affine_3_1*tmp_23 - p_affine_3_1*tmp_28 - p_affine_3_2*tmp_22 + p_affine_3_2*tmp_25);
      real_t a_0_0 = tmp_34*(tmp_19 + tmp_20 + tmp_21);
      real_t a_0_1 = -tmp_21*tmp_34;
      real_t a_0_2 = -tmp_20*tmp_34;
      real_t a_0_3 = -tmp_19*tmp_34;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_to_p0_div_1_affine_q0::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t tmp_19 = tmp_18*(tmp_4 - tmp_7);
      real_t tmp_20 = tmp_18*(tmp_11 - tmp_17);
      real_t tmp_21 = tmp_18*(tmp_15 - tmp_16);
      real_t tmp_22 = p_affine_0_0*p_affine_1_1;
      real_t tmp_23 = p_affine_0_0*p_affine_1_2;
      real_t tmp_24 = p_affine_2_1*p_affine_3_2;
      real_t tmp_25 = p_affine_0_1*p_affine_1_0;
      real_t tmp_26 = p_affine_0_1*p_affine_1_2;
      real_t tmp_27 = p_affine_2_2*p_affine_3_0;
      real_t tmp_28 = p_affine_0_2*p_affine_1_0;
      real_t tmp_29 = p_affine_0_2*p_affine_1_1;
      real_t tmp_30 = p_affine_2_0*p_affine_3_1;
      real_t tmp_31 = p_affine_2_2*p_affine_3_1;
      real_t tmp_32 = p_affine_2_0*p_affine_3_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_0;
      real_t tmp_34 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_24 - p_affine_0_0*tmp_31 + p_affine_0_1*tmp_27 - p_affine_0_1*tmp_32 + p_affine_0_2*tmp_30 - p_affine_0_2*tmp_33 - p_affine_1_0*tmp_24 + p_affine_1_0*tmp_31 - p_affine_1_1*tmp_27 + p_affine_1_1*tmp_32 - p_affine_1_2*tmp_30 + p_affine_1_2*tmp_33 + p_affine_2_0*tmp_26 - p_affine_2_0*tmp_29 - p_affine_2_1*tmp_23 + p_affine_2_1*tmp_28 + p_affine_2_2*tmp_22 - p_affine_2_2*tmp_25 - p_affine_3_0*tmp_26 + p_affine_3_0*tmp_29 + p_affine_3_1*tmp_23 - p_affine_3_1*tmp_28 - p_affine_3_2*tmp_22 + p_affine_3_2*tmp_25);
      real_t a_0_0 = tmp_34*(tmp_19 + tmp_20 + tmp_21);
      real_t a_0_1 = -tmp_21*tmp_34;
      real_t a_0_2 = -tmp_20*tmp_34;
      real_t a_0_3 = -tmp_19*tmp_34;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_to_p0_div_2_affine_q0::integrateAll( const std::array< Point3D, 3 >& , Matrix< real_t, 1, 3 >&  ) const
   {
      
   }

   void p1_to_p0_div_2_affine_q0::integrateRow0( const std::array< Point3D, 3 >& , Matrix< real_t, 1, 3 >&  ) const
   {
      
   }

   void p1_to_p0_div_2_affine_q0::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t tmp_22 = p_affine_0_0*p_affine_1_1;
      real_t tmp_23 = p_affine_0_0*p_affine_1_2;
      real_t tmp_24 = p_affine_2_1*p_affine_3_2;
      real_t tmp_25 = p_affine_0_1*p_affine_1_0;
      real_t tmp_26 = p_affine_0_1*p_affine_1_2;
      real_t tmp_27 = p_affine_2_2*p_affine_3_0;
      real_t tmp_28 = p_affine_0_2*p_affine_1_0;
      real_t tmp_29 = p_affine_0_2*p_affine_1_1;
      real_t tmp_30 = p_affine_2_0*p_affine_3_1;
      real_t tmp_31 = p_affine_2_2*p_affine_3_1;
      real_t tmp_32 = p_affine_2_0*p_affine_3_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_0;
      real_t tmp_34 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_24 - p_affine_0_0*tmp_31 + p_affine_0_1*tmp_27 - p_affine_0_1*tmp_32 + p_affine_0_2*tmp_30 - p_affine_0_2*tmp_33 - p_affine_1_0*tmp_24 + p_affine_1_0*tmp_31 - p_affine_1_1*tmp_27 + p_affine_1_1*tmp_32 - p_affine_1_2*tmp_30 + p_affine_1_2*tmp_33 + p_affine_2_0*tmp_26 - p_affine_2_0*tmp_29 - p_affine_2_1*tmp_23 + p_affine_2_1*tmp_28 + p_affine_2_2*tmp_22 - p_affine_2_2*tmp_25 - p_affine_3_0*tmp_26 + p_affine_3_0*tmp_29 + p_affine_3_1*tmp_23 - p_affine_3_1*tmp_28 - p_affine_3_2*tmp_22 + p_affine_3_2*tmp_25);
      real_t a_0_0 = tmp_34*(tmp_19 + tmp_20 + tmp_21);
      real_t a_0_1 = -tmp_21*tmp_34;
      real_t a_0_2 = -tmp_20*tmp_34;
      real_t a_0_3 = -tmp_19*tmp_34;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_to_p0_div_2_affine_q0::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t tmp_22 = p_affine_0_0*p_affine_1_1;
      real_t tmp_23 = p_affine_0_0*p_affine_1_2;
      real_t tmp_24 = p_affine_2_1*p_affine_3_2;
      real_t tmp_25 = p_affine_0_1*p_affine_1_0;
      real_t tmp_26 = p_affine_0_1*p_affine_1_2;
      real_t tmp_27 = p_affine_2_2*p_affine_3_0;
      real_t tmp_28 = p_affine_0_2*p_affine_1_0;
      real_t tmp_29 = p_affine_0_2*p_affine_1_1;
      real_t tmp_30 = p_affine_2_0*p_affine_3_1;
      real_t tmp_31 = p_affine_2_2*p_affine_3_1;
      real_t tmp_32 = p_affine_2_0*p_affine_3_2;
      real_t tmp_33 = p_affine_2_1*p_affine_3_0;
      real_t tmp_34 = 0.16666666666666663*std::abs(p_affine_0_0*tmp_24 - p_affine_0_0*tmp_31 + p_affine_0_1*tmp_27 - p_affine_0_1*tmp_32 + p_affine_0_2*tmp_30 - p_affine_0_2*tmp_33 - p_affine_1_0*tmp_24 + p_affine_1_0*tmp_31 - p_affine_1_1*tmp_27 + p_affine_1_1*tmp_32 - p_affine_1_2*tmp_30 + p_affine_1_2*tmp_33 + p_affine_2_0*tmp_26 - p_affine_2_0*tmp_29 - p_affine_2_1*tmp_23 + p_affine_2_1*tmp_28 + p_affine_2_2*tmp_22 - p_affine_2_2*tmp_25 - p_affine_3_0*tmp_26 + p_affine_3_0*tmp_29 + p_affine_3_1*tmp_23 - p_affine_3_1*tmp_28 - p_affine_3_2*tmp_22 + p_affine_3_2*tmp_25);
      real_t a_0_0 = tmp_34*(tmp_19 + tmp_20 + tmp_21);
      real_t a_0_1 = -tmp_21*tmp_34;
      real_t a_0_2 = -tmp_20*tmp_34;
      real_t a_0_3 = -tmp_19*tmp_34;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

} // namespace forms
} // namespace hyteg

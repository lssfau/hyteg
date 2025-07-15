/*
* Copyright (c) 2023 Andreas Burkhart
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

// This file has been generated with the AHFC. If buggy try fixing the generator itself.

#include "p2_to_p1_divT_affine_q2.hpp"

namespace hyteg {
namespace forms {

   void p2_to_p1_divT_0_affine_q2::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 6 >& elMat ) const
   {
      const real_t tmp0 = -coords[0][1] + coords[2][1];
      const real_t tmp1 = coords[0][1] - coords[1][1];
      const real_t tmp2 = 1.0 / (tmp0*(-coords[0][0] + coords[1][0]) + tmp1*(-coords[0][0] + coords[2][0]));
      const real_t tmp3 = tmp0*tmp2;
      const real_t tmp4 = tmp1*tmp2;
      const real_t tmp5 = std::abs(coords[0][0]*coords[1][1] - coords[0][0]*coords[2][1] - coords[0][1]*coords[1][0] + coords[0][1]*coords[2][0] + coords[1][0]*coords[2][1] - coords[1][1]*coords[2][0]);
      const real_t tmp6 = tmp5*(tmp3 + tmp4);
      const real_t tmp7 = tmp3*tmp5;
      const real_t tmp8 = tmp4*tmp5;
      elMat(0,0) = 5.8980598183211441e-17*tmp6;
      elMat(0,1) = 2.0816681711721685e-17*tmp6;
      elMat(0,2) = 5.5511151231257827e-17*tmp6;
      elMat(0,3) = 0.16666666666666674*tmp6;
      elMat(0,4) = 0.16666666666666669*tmp6;
      elMat(0,5) = 0.16666666666666677*tmp6;
      elMat(1,0) = -5.8980598183211441e-17*tmp7;
      elMat(1,1) = -2.0816681711721685e-17*tmp7;
      elMat(1,2) = -5.5511151231257827e-17*tmp7;
      elMat(1,3) = -0.16666666666666674*tmp7;
      elMat(1,4) = -0.16666666666666669*tmp7;
      elMat(1,5) = -0.16666666666666677*tmp7;
      elMat(2,0) = -5.8980598183211441e-17*tmp8;
      elMat(2,1) = -2.0816681711721685e-17*tmp8;
      elMat(2,2) = -5.5511151231257827e-17*tmp8;
      elMat(2,3) = -0.16666666666666674*tmp8;
      elMat(2,4) = -0.16666666666666669*tmp8;
      elMat(2,5) = -0.16666666666666677*tmp8;
   }

   void p2_to_p1_divT_0_affine_q2::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const   {
      const real_t tmp0 = -coords[0][1] + coords[2][1];
      const real_t tmp1 = coords[0][1] - coords[1][1];
      const real_t tmp2 = 1.0 / (tmp0*(-coords[0][0] + coords[1][0]) + tmp1*(-coords[0][0] + coords[2][0]));
      const real_t tmp3 = (tmp0*tmp2 + tmp1*tmp2)*std::abs(coords[0][0]*coords[1][1] - coords[0][0]*coords[2][1] - coords[0][1]*coords[1][0] + coords[0][1]*coords[2][0] + coords[1][0]*coords[2][1] - coords[1][1]*coords[2][0]);
      elMat(0,0) = 5.8980598183211441e-17*tmp3;
      elMat(0,1) = 2.0816681711721685e-17*tmp3;
      elMat(0,2) = 5.5511151231257827e-17*tmp3;
      elMat(0,3) = 0.16666666666666674*tmp3;
      elMat(0,4) = 0.16666666666666669*tmp3;
      elMat(0,5) = 0.16666666666666677*tmp3;
   }

   void p2_to_p1_divT_0_affine_q2::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 10 >& elMat ) const
   {
      const real_t tmp0 = -coords[0][1] + coords[1][1];
      const real_t tmp1 = -coords[0][2] + coords[2][2];
      const real_t tmp2 = tmp0*tmp1;
      const real_t tmp3 = -coords[0][1] + coords[2][1];
      const real_t tmp4 = -coords[0][2] + coords[1][2];
      const real_t tmp5 = tmp3*tmp4;
      const real_t tmp6 = -coords[0][0] + coords[1][0];
      const real_t tmp7 = -coords[0][2] + coords[3][2];
      const real_t tmp8 = tmp3*tmp7;
      const real_t tmp9 = -coords[0][0] + coords[2][0];
      const real_t tmp10 = -coords[0][1] + coords[3][1];
      const real_t tmp11 = -coords[0][0] + coords[3][0];
      const real_t tmp12 = tmp1*tmp10;
      const real_t tmp13 = tmp0*tmp7;
      const real_t tmp14 = 1.0 / (tmp10*tmp4*tmp9 + tmp11*tmp2 - tmp11*tmp5 - tmp12*tmp6 - tmp13*tmp9 + tmp6*tmp8);
      const real_t tmp15 = tmp14*(tmp2 - tmp5);
      const real_t tmp16 = tmp14*(tmp10*tmp4 - tmp13);
      const real_t tmp17 = tmp14*(-tmp12 + tmp8);
      const real_t tmp18 = coords[0][0]*coords[1][1];
      const real_t tmp19 = coords[0][0]*coords[1][2];
      const real_t tmp20 = coords[2][1]*coords[3][2];
      const real_t tmp21 = coords[0][1]*coords[1][0];
      const real_t tmp22 = coords[0][1]*coords[1][2];
      const real_t tmp23 = coords[2][2]*coords[3][0];
      const real_t tmp24 = coords[0][2]*coords[1][0];
      const real_t tmp25 = coords[0][2]*coords[1][1];
      const real_t tmp26 = coords[2][0]*coords[3][1];
      const real_t tmp27 = coords[2][2]*coords[3][1];
      const real_t tmp28 = coords[2][0]*coords[3][2];
      const real_t tmp29 = coords[2][1]*coords[3][0];
      const real_t tmp30 = std::abs(coords[0][0]*tmp20 - coords[0][0]*tmp27 + coords[0][1]*tmp23 - coords[0][1]*tmp28 + coords[0][2]*tmp26 - coords[0][2]*tmp29 - coords[1][0]*tmp20 + coords[1][0]*tmp27 - coords[1][1]*tmp23 + coords[1][1]*tmp28 - coords[1][2]*tmp26 + coords[1][2]*tmp29 + coords[2][0]*tmp22 - coords[2][0]*tmp25 - coords[2][1]*tmp19 + coords[2][1]*tmp24 + coords[2][2]*tmp18 - coords[2][2]*tmp21 - coords[3][0]*tmp22 + coords[3][0]*tmp25 + coords[3][1]*tmp19 - coords[3][1]*tmp24 - coords[3][2]*tmp18 + coords[3][2]*tmp21);
      const real_t tmp31 = tmp30*(tmp15 + tmp16 + tmp17);
      const real_t tmp32 = tmp17*tmp30;
      const real_t tmp33 = tmp16*tmp30;
      const real_t tmp34 = tmp15*tmp30;
      elMat(0,0) = -0.0083333333333333297*tmp31;
      elMat(0,1) = -0.0083333333333333315*tmp31;
      elMat(0,2) = -0.0083333333333333384*tmp31;
      elMat(0,3) = -0.0083333333333333245*tmp31;
      elMat(0,4) = 0.033333333333333277*tmp31;
      elMat(0,5) = 0.033333333333333291*tmp31;
      elMat(0,6) = 0.033333333333333305*tmp31;
      elMat(0,7) = 0.033333333333333284*tmp31;
      elMat(0,8) = 0.033333333333333375*tmp31;
      elMat(0,9) = 0.033333333333333368*tmp31;
      elMat(1,0) = 0.0083333333333333297*tmp32;
      elMat(1,1) = 0.0083333333333333315*tmp32;
      elMat(1,2) = 0.0083333333333333384*tmp32;
      elMat(1,3) = 0.0083333333333333245*tmp32;
      elMat(1,4) = -0.033333333333333277*tmp32;
      elMat(1,5) = -0.033333333333333291*tmp32;
      elMat(1,6) = -0.033333333333333305*tmp32;
      elMat(1,7) = -0.033333333333333284*tmp32;
      elMat(1,8) = -0.033333333333333375*tmp32;
      elMat(1,9) = -0.033333333333333368*tmp32;
      elMat(2,0) = 0.0083333333333333297*tmp33;
      elMat(2,1) = 0.0083333333333333315*tmp33;
      elMat(2,2) = 0.0083333333333333384*tmp33;
      elMat(2,3) = 0.0083333333333333245*tmp33;
      elMat(2,4) = -0.033333333333333277*tmp33;
      elMat(2,5) = -0.033333333333333291*tmp33;
      elMat(2,6) = -0.033333333333333305*tmp33;
      elMat(2,7) = -0.033333333333333284*tmp33;
      elMat(2,8) = -0.033333333333333375*tmp33;
      elMat(2,9) = -0.033333333333333368*tmp33;
      elMat(3,0) = 0.0083333333333333297*tmp34;
      elMat(3,1) = 0.0083333333333333315*tmp34;
      elMat(3,2) = 0.0083333333333333384*tmp34;
      elMat(3,3) = 0.0083333333333333245*tmp34;
      elMat(3,4) = -0.033333333333333277*tmp34;
      elMat(3,5) = -0.033333333333333291*tmp34;
      elMat(3,6) = -0.033333333333333305*tmp34;
      elMat(3,7) = -0.033333333333333284*tmp34;
      elMat(3,8) = -0.033333333333333375*tmp34;
      elMat(3,9) = -0.033333333333333368*tmp34;
   }

   void p2_to_p1_divT_0_affine_q2::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const   {
      const real_t tmp0 = -coords[0][1] + coords[1][1];
      const real_t tmp1 = -coords[0][2] + coords[2][2];
      const real_t tmp2 = tmp0*tmp1;
      const real_t tmp3 = -coords[0][1] + coords[2][1];
      const real_t tmp4 = -coords[0][2] + coords[1][2];
      const real_t tmp5 = tmp3*tmp4;
      const real_t tmp6 = -coords[0][0] + coords[1][0];
      const real_t tmp7 = -coords[0][2] + coords[3][2];
      const real_t tmp8 = tmp3*tmp7;
      const real_t tmp9 = -coords[0][0] + coords[2][0];
      const real_t tmp10 = -coords[0][1] + coords[3][1];
      const real_t tmp11 = -coords[0][0] + coords[3][0];
      const real_t tmp12 = tmp1*tmp10;
      const real_t tmp13 = tmp0*tmp7;
      const real_t tmp14 = 1.0 / (tmp10*tmp4*tmp9 + tmp11*tmp2 - tmp11*tmp5 - tmp12*tmp6 - tmp13*tmp9 + tmp6*tmp8);
      const real_t tmp15 = coords[0][0]*coords[1][1];
      const real_t tmp16 = coords[0][0]*coords[1][2];
      const real_t tmp17 = coords[2][1]*coords[3][2];
      const real_t tmp18 = coords[0][1]*coords[1][0];
      const real_t tmp19 = coords[0][1]*coords[1][2];
      const real_t tmp20 = coords[2][2]*coords[3][0];
      const real_t tmp21 = coords[0][2]*coords[1][0];
      const real_t tmp22 = coords[0][2]*coords[1][1];
      const real_t tmp23 = coords[2][0]*coords[3][1];
      const real_t tmp24 = coords[2][2]*coords[3][1];
      const real_t tmp25 = coords[2][0]*coords[3][2];
      const real_t tmp26 = coords[2][1]*coords[3][0];
      const real_t tmp27 = (tmp14*(-tmp12 + tmp8) + tmp14*(tmp2 - tmp5) + tmp14*(tmp10*tmp4 - tmp13))*std::abs(coords[0][0]*tmp17 - coords[0][0]*tmp24 + coords[0][1]*tmp20 - coords[0][1]*tmp25 + coords[0][2]*tmp23 - coords[0][2]*tmp26 - coords[1][0]*tmp17 + coords[1][0]*tmp24 - coords[1][1]*tmp20 + coords[1][1]*tmp25 - coords[1][2]*tmp23 + coords[1][2]*tmp26 + coords[2][0]*tmp19 - coords[2][0]*tmp22 - coords[2][1]*tmp16 + coords[2][1]*tmp21 + coords[2][2]*tmp15 - coords[2][2]*tmp18 - coords[3][0]*tmp19 + coords[3][0]*tmp22 + coords[3][1]*tmp16 - coords[3][1]*tmp21 - coords[3][2]*tmp15 + coords[3][2]*tmp18);
      elMat(0,0) = -0.0083333333333333297*tmp27;
      elMat(0,1) = -0.0083333333333333315*tmp27;
      elMat(0,2) = -0.0083333333333333384*tmp27;
      elMat(0,3) = -0.0083333333333333245*tmp27;
      elMat(0,4) = 0.033333333333333277*tmp27;
      elMat(0,5) = 0.033333333333333291*tmp27;
      elMat(0,6) = 0.033333333333333305*tmp27;
      elMat(0,7) = 0.033333333333333284*tmp27;
      elMat(0,8) = 0.033333333333333375*tmp27;
      elMat(0,9) = 0.033333333333333368*tmp27;
   }

   void p2_to_p1_divT_1_affine_q2::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 6 >& elMat ) const
   {
      const real_t tmp0 = -coords[0][0] + coords[1][0];
      const real_t tmp1 = coords[0][0] - coords[2][0];
      const real_t tmp2 = 1.0 / (tmp0*(-coords[0][1] + coords[2][1]) + tmp1*(-coords[0][1] + coords[1][1]));
      const real_t tmp3 = tmp0*tmp2;
      const real_t tmp4 = tmp1*tmp2;
      const real_t tmp5 = std::abs(coords[0][0]*coords[1][1] - coords[0][0]*coords[2][1] - coords[0][1]*coords[1][0] + coords[0][1]*coords[2][0] + coords[1][0]*coords[2][1] - coords[1][1]*coords[2][0]);
      const real_t tmp6 = tmp5*(tmp3 + tmp4);
      const real_t tmp7 = tmp4*tmp5;
      const real_t tmp8 = tmp3*tmp5;
      elMat(0,0) = 5.8980598183211441e-17*tmp6;
      elMat(0,1) = 2.0816681711721685e-17*tmp6;
      elMat(0,2) = 5.5511151231257827e-17*tmp6;
      elMat(0,3) = 0.16666666666666674*tmp6;
      elMat(0,4) = 0.16666666666666669*tmp6;
      elMat(0,5) = 0.16666666666666677*tmp6;
      elMat(1,0) = -5.8980598183211441e-17*tmp7;
      elMat(1,1) = -2.0816681711721685e-17*tmp7;
      elMat(1,2) = -5.5511151231257827e-17*tmp7;
      elMat(1,3) = -0.16666666666666674*tmp7;
      elMat(1,4) = -0.16666666666666669*tmp7;
      elMat(1,5) = -0.16666666666666677*tmp7;
      elMat(2,0) = -5.8980598183211441e-17*tmp8;
      elMat(2,1) = -2.0816681711721685e-17*tmp8;
      elMat(2,2) = -5.5511151231257827e-17*tmp8;
      elMat(2,3) = -0.16666666666666674*tmp8;
      elMat(2,4) = -0.16666666666666669*tmp8;
      elMat(2,5) = -0.16666666666666677*tmp8;
   }

   void p2_to_p1_divT_1_affine_q2::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const   {
      const real_t tmp0 = -coords[0][0] + coords[1][0];
      const real_t tmp1 = coords[0][0] - coords[2][0];
      const real_t tmp2 = 1.0 / (tmp0*(-coords[0][1] + coords[2][1]) + tmp1*(-coords[0][1] + coords[1][1]));
      const real_t tmp3 = (tmp0*tmp2 + tmp1*tmp2)*std::abs(coords[0][0]*coords[1][1] - coords[0][0]*coords[2][1] - coords[0][1]*coords[1][0] + coords[0][1]*coords[2][0] + coords[1][0]*coords[2][1] - coords[1][1]*coords[2][0]);
      elMat(0,0) = 5.8980598183211441e-17*tmp3;
      elMat(0,1) = 2.0816681711721685e-17*tmp3;
      elMat(0,2) = 5.5511151231257827e-17*tmp3;
      elMat(0,3) = 0.16666666666666674*tmp3;
      elMat(0,4) = 0.16666666666666669*tmp3;
      elMat(0,5) = 0.16666666666666677*tmp3;
   }

   void p2_to_p1_divT_1_affine_q2::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 10 >& elMat ) const
   {
      const real_t tmp0 = -coords[0][0] + coords[1][0];
      const real_t tmp1 = -coords[0][2] + coords[2][2];
      const real_t tmp2 = tmp0*tmp1;
      const real_t tmp3 = -coords[0][0] + coords[2][0];
      const real_t tmp4 = -coords[0][2] + coords[1][2];
      const real_t tmp5 = -coords[0][1] + coords[2][1];
      const real_t tmp6 = -coords[0][2] + coords[3][2];
      const real_t tmp7 = tmp0*tmp6;
      const real_t tmp8 = -coords[0][1] + coords[3][1];
      const real_t tmp9 = -coords[0][1] + coords[1][1];
      const real_t tmp10 = -coords[0][0] + coords[3][0];
      const real_t tmp11 = tmp3*tmp6;
      const real_t tmp12 = tmp10*tmp4;
      const real_t tmp13 = 1.0 / (tmp1*tmp10*tmp9 - tmp11*tmp9 - tmp12*tmp5 - tmp2*tmp8 + tmp3*tmp4*tmp8 + tmp5*tmp7);
      const real_t tmp14 = tmp13*(-tmp2 + tmp3*tmp4);
      const real_t tmp15 = tmp13*(-tmp12 + tmp7);
      const real_t tmp16 = tmp13*(tmp1*tmp10 - tmp11);
      const real_t tmp17 = coords[0][0]*coords[1][1];
      const real_t tmp18 = coords[0][0]*coords[1][2];
      const real_t tmp19 = coords[2][1]*coords[3][2];
      const real_t tmp20 = coords[0][1]*coords[1][0];
      const real_t tmp21 = coords[0][1]*coords[1][2];
      const real_t tmp22 = coords[2][2]*coords[3][0];
      const real_t tmp23 = coords[0][2]*coords[1][0];
      const real_t tmp24 = coords[0][2]*coords[1][1];
      const real_t tmp25 = coords[2][0]*coords[3][1];
      const real_t tmp26 = coords[2][2]*coords[3][1];
      const real_t tmp27 = coords[2][0]*coords[3][2];
      const real_t tmp28 = coords[2][1]*coords[3][0];
      const real_t tmp29 = std::abs(coords[0][0]*tmp19 - coords[0][0]*tmp26 + coords[0][1]*tmp22 - coords[0][1]*tmp27 + coords[0][2]*tmp25 - coords[0][2]*tmp28 - coords[1][0]*tmp19 + coords[1][0]*tmp26 - coords[1][1]*tmp22 + coords[1][1]*tmp27 - coords[1][2]*tmp25 + coords[1][2]*tmp28 + coords[2][0]*tmp21 - coords[2][0]*tmp24 - coords[2][1]*tmp18 + coords[2][1]*tmp23 + coords[2][2]*tmp17 - coords[2][2]*tmp20 - coords[3][0]*tmp21 + coords[3][0]*tmp24 + coords[3][1]*tmp18 - coords[3][1]*tmp23 - coords[3][2]*tmp17 + coords[3][2]*tmp20);
      const real_t tmp30 = tmp29*(tmp14 + tmp15 + tmp16);
      const real_t tmp31 = tmp16*tmp29;
      const real_t tmp32 = tmp15*tmp29;
      const real_t tmp33 = tmp14*tmp29;
      elMat(0,0) = -0.0083333333333333297*tmp30;
      elMat(0,1) = -0.0083333333333333315*tmp30;
      elMat(0,2) = -0.0083333333333333384*tmp30;
      elMat(0,3) = -0.0083333333333333245*tmp30;
      elMat(0,4) = 0.033333333333333277*tmp30;
      elMat(0,5) = 0.033333333333333291*tmp30;
      elMat(0,6) = 0.033333333333333305*tmp30;
      elMat(0,7) = 0.033333333333333284*tmp30;
      elMat(0,8) = 0.033333333333333375*tmp30;
      elMat(0,9) = 0.033333333333333368*tmp30;
      elMat(1,0) = 0.0083333333333333297*tmp31;
      elMat(1,1) = 0.0083333333333333315*tmp31;
      elMat(1,2) = 0.0083333333333333384*tmp31;
      elMat(1,3) = 0.0083333333333333245*tmp31;
      elMat(1,4) = -0.033333333333333277*tmp31;
      elMat(1,5) = -0.033333333333333291*tmp31;
      elMat(1,6) = -0.033333333333333305*tmp31;
      elMat(1,7) = -0.033333333333333284*tmp31;
      elMat(1,8) = -0.033333333333333375*tmp31;
      elMat(1,9) = -0.033333333333333368*tmp31;
      elMat(2,0) = 0.0083333333333333297*tmp32;
      elMat(2,1) = 0.0083333333333333315*tmp32;
      elMat(2,2) = 0.0083333333333333384*tmp32;
      elMat(2,3) = 0.0083333333333333245*tmp32;
      elMat(2,4) = -0.033333333333333277*tmp32;
      elMat(2,5) = -0.033333333333333291*tmp32;
      elMat(2,6) = -0.033333333333333305*tmp32;
      elMat(2,7) = -0.033333333333333284*tmp32;
      elMat(2,8) = -0.033333333333333375*tmp32;
      elMat(2,9) = -0.033333333333333368*tmp32;
      elMat(3,0) = 0.0083333333333333297*tmp33;
      elMat(3,1) = 0.0083333333333333315*tmp33;
      elMat(3,2) = 0.0083333333333333384*tmp33;
      elMat(3,3) = 0.0083333333333333245*tmp33;
      elMat(3,4) = -0.033333333333333277*tmp33;
      elMat(3,5) = -0.033333333333333291*tmp33;
      elMat(3,6) = -0.033333333333333305*tmp33;
      elMat(3,7) = -0.033333333333333284*tmp33;
      elMat(3,8) = -0.033333333333333375*tmp33;
      elMat(3,9) = -0.033333333333333368*tmp33;
   }

   void p2_to_p1_divT_1_affine_q2::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const   {
      const real_t tmp0 = -coords[0][0] + coords[1][0];
      const real_t tmp1 = -coords[0][2] + coords[2][2];
      const real_t tmp2 = tmp0*tmp1;
      const real_t tmp3 = -coords[0][0] + coords[2][0];
      const real_t tmp4 = -coords[0][2] + coords[1][2];
      const real_t tmp5 = -coords[0][1] + coords[2][1];
      const real_t tmp6 = -coords[0][2] + coords[3][2];
      const real_t tmp7 = tmp0*tmp6;
      const real_t tmp8 = -coords[0][1] + coords[3][1];
      const real_t tmp9 = -coords[0][1] + coords[1][1];
      const real_t tmp10 = -coords[0][0] + coords[3][0];
      const real_t tmp11 = tmp3*tmp6;
      const real_t tmp12 = tmp10*tmp4;
      const real_t tmp13 = 1.0 / (tmp1*tmp10*tmp9 - tmp11*tmp9 - tmp12*tmp5 - tmp2*tmp8 + tmp3*tmp4*tmp8 + tmp5*tmp7);
      const real_t tmp14 = coords[0][0]*coords[1][1];
      const real_t tmp15 = coords[0][0]*coords[1][2];
      const real_t tmp16 = coords[2][1]*coords[3][2];
      const real_t tmp17 = coords[0][1]*coords[1][0];
      const real_t tmp18 = coords[0][1]*coords[1][2];
      const real_t tmp19 = coords[2][2]*coords[3][0];
      const real_t tmp20 = coords[0][2]*coords[1][0];
      const real_t tmp21 = coords[0][2]*coords[1][1];
      const real_t tmp22 = coords[2][0]*coords[3][1];
      const real_t tmp23 = coords[2][2]*coords[3][1];
      const real_t tmp24 = coords[2][0]*coords[3][2];
      const real_t tmp25 = coords[2][1]*coords[3][0];
      const real_t tmp26 = (tmp13*(-tmp12 + tmp7) + tmp13*(-tmp2 + tmp3*tmp4) + tmp13*(tmp1*tmp10 - tmp11))*std::abs(coords[0][0]*tmp16 - coords[0][0]*tmp23 + coords[0][1]*tmp19 - coords[0][1]*tmp24 + coords[0][2]*tmp22 - coords[0][2]*tmp25 - coords[1][0]*tmp16 + coords[1][0]*tmp23 - coords[1][1]*tmp19 + coords[1][1]*tmp24 - coords[1][2]*tmp22 + coords[1][2]*tmp25 + coords[2][0]*tmp18 - coords[2][0]*tmp21 - coords[2][1]*tmp15 + coords[2][1]*tmp20 + coords[2][2]*tmp14 - coords[2][2]*tmp17 - coords[3][0]*tmp18 + coords[3][0]*tmp21 + coords[3][1]*tmp15 - coords[3][1]*tmp20 - coords[3][2]*tmp14 + coords[3][2]*tmp17);
      elMat(0,0) = -0.0083333333333333297*tmp26;
      elMat(0,1) = -0.0083333333333333315*tmp26;
      elMat(0,2) = -0.0083333333333333384*tmp26;
      elMat(0,3) = -0.0083333333333333245*tmp26;
      elMat(0,4) = 0.033333333333333277*tmp26;
      elMat(0,5) = 0.033333333333333291*tmp26;
      elMat(0,6) = 0.033333333333333305*tmp26;
      elMat(0,7) = 0.033333333333333284*tmp26;
      elMat(0,8) = 0.033333333333333375*tmp26;
      elMat(0,9) = 0.033333333333333368*tmp26;
   }

   void p2_to_p1_divT_2_affine_q2::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 6 >& elMat ) const
   {
      (void)(coords);
      elMat(0,0) = 0;
      elMat(0,1) = 0;
      elMat(0,2) = 0;
      elMat(0,3) = 0;
      elMat(0,4) = 0;
      elMat(0,5) = 0;
      elMat(1,0) = 0;
      elMat(1,1) = 0;
      elMat(1,2) = 0;
      elMat(1,3) = 0;
      elMat(1,4) = 0;
      elMat(1,5) = 0;
      elMat(2,0) = 0;
      elMat(2,1) = 0;
      elMat(2,2) = 0;
      elMat(2,3) = 0;
      elMat(2,4) = 0;
      elMat(2,5) = 0;
   }

   void p2_to_p1_divT_2_affine_q2::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const   {
      (void)(coords);
      elMat(0,0) = 0;
      elMat(0,1) = 0;
      elMat(0,2) = 0;
      elMat(0,3) = 0;
      elMat(0,4) = 0;
      elMat(0,5) = 0;
   }

   void p2_to_p1_divT_2_affine_q2::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 10 >& elMat ) const
   {
      const real_t tmp0 = -coords[0][0] + coords[1][0];
      const real_t tmp1 = -coords[0][1] + coords[2][1];
      const real_t tmp2 = tmp0*tmp1;
      const real_t tmp3 = -coords[0][0] + coords[2][0];
      const real_t tmp4 = -coords[0][1] + coords[1][1];
      const real_t tmp5 = tmp3*tmp4;
      const real_t tmp6 = -coords[0][2] + coords[3][2];
      const real_t tmp7 = -coords[0][2] + coords[1][2];
      const real_t tmp8 = -coords[0][1] + coords[3][1];
      const real_t tmp9 = tmp3*tmp8;
      const real_t tmp10 = -coords[0][2] + coords[2][2];
      const real_t tmp11 = -coords[0][0] + coords[3][0];
      const real_t tmp12 = tmp0*tmp8;
      const real_t tmp13 = tmp1*tmp11;
      const real_t tmp14 = 1.0 / (tmp10*tmp11*tmp4 - tmp10*tmp12 - tmp13*tmp7 + tmp2*tmp6 - tmp5*tmp6 + tmp7*tmp9);
      const real_t tmp15 = tmp14*(tmp2 - tmp5);
      const real_t tmp16 = tmp14*(tmp11*tmp4 - tmp12);
      const real_t tmp17 = tmp14*(-tmp13 + tmp9);
      const real_t tmp18 = coords[0][0]*coords[1][1];
      const real_t tmp19 = coords[0][0]*coords[1][2];
      const real_t tmp20 = coords[2][1]*coords[3][2];
      const real_t tmp21 = coords[0][1]*coords[1][0];
      const real_t tmp22 = coords[0][1]*coords[1][2];
      const real_t tmp23 = coords[2][2]*coords[3][0];
      const real_t tmp24 = coords[0][2]*coords[1][0];
      const real_t tmp25 = coords[0][2]*coords[1][1];
      const real_t tmp26 = coords[2][0]*coords[3][1];
      const real_t tmp27 = coords[2][2]*coords[3][1];
      const real_t tmp28 = coords[2][0]*coords[3][2];
      const real_t tmp29 = coords[2][1]*coords[3][0];
      const real_t tmp30 = std::abs(coords[0][0]*tmp20 - coords[0][0]*tmp27 + coords[0][1]*tmp23 - coords[0][1]*tmp28 + coords[0][2]*tmp26 - coords[0][2]*tmp29 - coords[1][0]*tmp20 + coords[1][0]*tmp27 - coords[1][1]*tmp23 + coords[1][1]*tmp28 - coords[1][2]*tmp26 + coords[1][2]*tmp29 + coords[2][0]*tmp22 - coords[2][0]*tmp25 - coords[2][1]*tmp19 + coords[2][1]*tmp24 + coords[2][2]*tmp18 - coords[2][2]*tmp21 - coords[3][0]*tmp22 + coords[3][0]*tmp25 + coords[3][1]*tmp19 - coords[3][1]*tmp24 - coords[3][2]*tmp18 + coords[3][2]*tmp21);
      const real_t tmp31 = tmp30*(tmp15 + tmp16 + tmp17);
      const real_t tmp32 = tmp17*tmp30;
      const real_t tmp33 = tmp16*tmp30;
      const real_t tmp34 = tmp15*tmp30;
      elMat(0,0) = -0.0083333333333333297*tmp31;
      elMat(0,1) = -0.0083333333333333315*tmp31;
      elMat(0,2) = -0.0083333333333333384*tmp31;
      elMat(0,3) = -0.0083333333333333245*tmp31;
      elMat(0,4) = 0.033333333333333277*tmp31;
      elMat(0,5) = 0.033333333333333291*tmp31;
      elMat(0,6) = 0.033333333333333305*tmp31;
      elMat(0,7) = 0.033333333333333284*tmp31;
      elMat(0,8) = 0.033333333333333375*tmp31;
      elMat(0,9) = 0.033333333333333368*tmp31;
      elMat(1,0) = 0.0083333333333333297*tmp32;
      elMat(1,1) = 0.0083333333333333315*tmp32;
      elMat(1,2) = 0.0083333333333333384*tmp32;
      elMat(1,3) = 0.0083333333333333245*tmp32;
      elMat(1,4) = -0.033333333333333277*tmp32;
      elMat(1,5) = -0.033333333333333291*tmp32;
      elMat(1,6) = -0.033333333333333305*tmp32;
      elMat(1,7) = -0.033333333333333284*tmp32;
      elMat(1,8) = -0.033333333333333375*tmp32;
      elMat(1,9) = -0.033333333333333368*tmp32;
      elMat(2,0) = 0.0083333333333333297*tmp33;
      elMat(2,1) = 0.0083333333333333315*tmp33;
      elMat(2,2) = 0.0083333333333333384*tmp33;
      elMat(2,3) = 0.0083333333333333245*tmp33;
      elMat(2,4) = -0.033333333333333277*tmp33;
      elMat(2,5) = -0.033333333333333291*tmp33;
      elMat(2,6) = -0.033333333333333305*tmp33;
      elMat(2,7) = -0.033333333333333284*tmp33;
      elMat(2,8) = -0.033333333333333375*tmp33;
      elMat(2,9) = -0.033333333333333368*tmp33;
      elMat(3,0) = 0.0083333333333333297*tmp34;
      elMat(3,1) = 0.0083333333333333315*tmp34;
      elMat(3,2) = 0.0083333333333333384*tmp34;
      elMat(3,3) = 0.0083333333333333245*tmp34;
      elMat(3,4) = -0.033333333333333277*tmp34;
      elMat(3,5) = -0.033333333333333291*tmp34;
      elMat(3,6) = -0.033333333333333305*tmp34;
      elMat(3,7) = -0.033333333333333284*tmp34;
      elMat(3,8) = -0.033333333333333375*tmp34;
      elMat(3,9) = -0.033333333333333368*tmp34;
   }

   void p2_to_p1_divT_2_affine_q2::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const   {
      const real_t tmp0 = -coords[0][0] + coords[1][0];
      const real_t tmp1 = -coords[0][1] + coords[2][1];
      const real_t tmp2 = tmp0*tmp1;
      const real_t tmp3 = -coords[0][0] + coords[2][0];
      const real_t tmp4 = -coords[0][1] + coords[1][1];
      const real_t tmp5 = tmp3*tmp4;
      const real_t tmp6 = -coords[0][2] + coords[3][2];
      const real_t tmp7 = -coords[0][2] + coords[1][2];
      const real_t tmp8 = -coords[0][1] + coords[3][1];
      const real_t tmp9 = tmp3*tmp8;
      const real_t tmp10 = -coords[0][2] + coords[2][2];
      const real_t tmp11 = -coords[0][0] + coords[3][0];
      const real_t tmp12 = tmp0*tmp8;
      const real_t tmp13 = tmp1*tmp11;
      const real_t tmp14 = 1.0 / (tmp10*tmp11*tmp4 - tmp10*tmp12 - tmp13*tmp7 + tmp2*tmp6 - tmp5*tmp6 + tmp7*tmp9);
      const real_t tmp15 = coords[0][0]*coords[1][1];
      const real_t tmp16 = coords[0][0]*coords[1][2];
      const real_t tmp17 = coords[2][1]*coords[3][2];
      const real_t tmp18 = coords[0][1]*coords[1][0];
      const real_t tmp19 = coords[0][1]*coords[1][2];
      const real_t tmp20 = coords[2][2]*coords[3][0];
      const real_t tmp21 = coords[0][2]*coords[1][0];
      const real_t tmp22 = coords[0][2]*coords[1][1];
      const real_t tmp23 = coords[2][0]*coords[3][1];
      const real_t tmp24 = coords[2][2]*coords[3][1];
      const real_t tmp25 = coords[2][0]*coords[3][2];
      const real_t tmp26 = coords[2][1]*coords[3][0];
      const real_t tmp27 = (tmp14*(-tmp13 + tmp9) + tmp14*(tmp2 - tmp5) + tmp14*(tmp11*tmp4 - tmp12))*std::abs(coords[0][0]*tmp17 - coords[0][0]*tmp24 + coords[0][1]*tmp20 - coords[0][1]*tmp25 + coords[0][2]*tmp23 - coords[0][2]*tmp26 - coords[1][0]*tmp17 + coords[1][0]*tmp24 - coords[1][1]*tmp20 + coords[1][1]*tmp25 - coords[1][2]*tmp23 + coords[1][2]*tmp26 + coords[2][0]*tmp19 - coords[2][0]*tmp22 - coords[2][1]*tmp16 + coords[2][1]*tmp21 + coords[2][2]*tmp15 - coords[2][2]*tmp18 - coords[3][0]*tmp19 + coords[3][0]*tmp22 + coords[3][1]*tmp16 - coords[3][1]*tmp21 - coords[3][2]*tmp15 + coords[3][2]*tmp18);
      elMat(0,0) = -0.0083333333333333297*tmp27;
      elMat(0,1) = -0.0083333333333333315*tmp27;
      elMat(0,2) = -0.0083333333333333384*tmp27;
      elMat(0,3) = -0.0083333333333333245*tmp27;
      elMat(0,4) = 0.033333333333333277*tmp27;
      elMat(0,5) = 0.033333333333333291*tmp27;
      elMat(0,6) = 0.033333333333333305*tmp27;
      elMat(0,7) = 0.033333333333333284*tmp27;
      elMat(0,8) = 0.033333333333333375*tmp27;
      elMat(0,9) = 0.033333333333333368*tmp27;
   }

} // namespace forms
} // namespace hyteg
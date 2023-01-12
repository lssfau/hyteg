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

#include "p2_mass_affine_qe.hpp"

namespace hyteg {
namespace forms {

   void p2_mass_affine_qe::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 6, 6 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_0 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_1 = (1.0/60.0)*tmp_0;
      real_t tmp_2 = -1.0/360.0*tmp_0;
      real_t tmp_3 = -1.0/90.0*tmp_0;
      real_t tmp_4 = (4.0/45.0)*tmp_0;
      real_t tmp_5 = (2.0/45.0)*tmp_0;
      real_t a_0_0 = tmp_1;
      real_t a_0_1 = tmp_2;
      real_t a_0_2 = tmp_2;
      real_t a_0_3 = tmp_3;
      real_t a_0_4 = 0;
      real_t a_0_5 = 0;
      real_t a_1_0 = tmp_2;
      real_t a_1_1 = tmp_1;
      real_t a_1_2 = tmp_2;
      real_t a_1_3 = 0;
      real_t a_1_4 = tmp_3;
      real_t a_1_5 = 0;
      real_t a_2_0 = tmp_2;
      real_t a_2_1 = tmp_2;
      real_t a_2_2 = tmp_1;
      real_t a_2_3 = 0;
      real_t a_2_4 = 0;
      real_t a_2_5 = tmp_3;
      real_t a_3_0 = tmp_3;
      real_t a_3_1 = 0;
      real_t a_3_2 = 0;
      real_t a_3_3 = tmp_4;
      real_t a_3_4 = tmp_5;
      real_t a_3_5 = tmp_5;
      real_t a_4_0 = 0;
      real_t a_4_1 = tmp_3;
      real_t a_4_2 = 0;
      real_t a_4_3 = tmp_5;
      real_t a_4_4 = tmp_4;
      real_t a_4_5 = tmp_5;
      real_t a_5_0 = 0;
      real_t a_5_1 = 0;
      real_t a_5_2 = tmp_3;
      real_t a_5_3 = tmp_5;
      real_t a_5_4 = tmp_5;
      real_t a_5_5 = tmp_4;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
      (elMat(0, 4)) = a_0_4;
      (elMat(0, 5)) = a_0_5;
      (elMat(1, 0)) = a_1_0;
      (elMat(1, 1)) = a_1_1;
      (elMat(1, 2)) = a_1_2;
      (elMat(1, 3)) = a_1_3;
      (elMat(1, 4)) = a_1_4;
      (elMat(1, 5)) = a_1_5;
      (elMat(2, 0)) = a_2_0;
      (elMat(2, 1)) = a_2_1;
      (elMat(2, 2)) = a_2_2;
      (elMat(2, 3)) = a_2_3;
      (elMat(2, 4)) = a_2_4;
      (elMat(2, 5)) = a_2_5;
      (elMat(3, 0)) = a_3_0;
      (elMat(3, 1)) = a_3_1;
      (elMat(3, 2)) = a_3_2;
      (elMat(3, 3)) = a_3_3;
      (elMat(3, 4)) = a_3_4;
      (elMat(3, 5)) = a_3_5;
      (elMat(4, 0)) = a_4_0;
      (elMat(4, 1)) = a_4_1;
      (elMat(4, 2)) = a_4_2;
      (elMat(4, 3)) = a_4_3;
      (elMat(4, 4)) = a_4_4;
      (elMat(4, 5)) = a_4_5;
      (elMat(5, 0)) = a_5_0;
      (elMat(5, 1)) = a_5_1;
      (elMat(5, 2)) = a_5_2;
      (elMat(5, 3)) = a_5_3;
      (elMat(5, 4)) = a_5_4;
      (elMat(5, 5)) = a_5_5;
   }

   void p2_mass_affine_qe::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 6 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t tmp_0 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_1 = -1.0/360.0*tmp_0;
      real_t a_0_0 = (1.0/60.0)*tmp_0;
      real_t a_0_1 = tmp_1;
      real_t a_0_2 = tmp_1;
      real_t a_0_3 = -1.0/90.0*tmp_0;
      real_t a_0_4 = 0;
      real_t a_0_5 = 0;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
      (elMat(0, 4)) = a_0_4;
      (elMat(0, 5)) = a_0_5;
   }

   void p2_mass_affine_qe::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 10, 10 >& elMat ) const
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
      real_t tmp_0 = p_affine_0_0*p_affine_1_1;
      real_t tmp_1 = p_affine_0_0*p_affine_1_2;
      real_t tmp_2 = p_affine_2_1*p_affine_3_2;
      real_t tmp_3 = p_affine_0_1*p_affine_1_0;
      real_t tmp_4 = p_affine_0_1*p_affine_1_2;
      real_t tmp_5 = p_affine_2_2*p_affine_3_0;
      real_t tmp_6 = p_affine_0_2*p_affine_1_0;
      real_t tmp_7 = p_affine_0_2*p_affine_1_1;
      real_t tmp_8 = p_affine_2_0*p_affine_3_1;
      real_t tmp_9 = p_affine_2_2*p_affine_3_1;
      real_t tmp_10 = p_affine_2_0*p_affine_3_2;
      real_t tmp_11 = p_affine_2_1*p_affine_3_0;
      real_t tmp_12 = std::abs(p_affine_0_0*tmp_2 - p_affine_0_0*tmp_9 - p_affine_0_1*tmp_10 + p_affine_0_1*tmp_5 - p_affine_0_2*tmp_11 + p_affine_0_2*tmp_8 - p_affine_1_0*tmp_2 + p_affine_1_0*tmp_9 + p_affine_1_1*tmp_10 - p_affine_1_1*tmp_5 + p_affine_1_2*tmp_11 - p_affine_1_2*tmp_8 + p_affine_2_0*tmp_4 - p_affine_2_0*tmp_7 - p_affine_2_1*tmp_1 + p_affine_2_1*tmp_6 + p_affine_2_2*tmp_0 - p_affine_2_2*tmp_3 - p_affine_3_0*tmp_4 + p_affine_3_0*tmp_7 + p_affine_3_1*tmp_1 - p_affine_3_1*tmp_6 - p_affine_3_2*tmp_0 + p_affine_3_2*tmp_3);
      real_t tmp_13 = 0.00039682539682555154*tmp_12;
      real_t tmp_14 = 0.0003968253968249201*tmp_12;
      real_t tmp_15 = 0.00039682539682550644*tmp_12;
      real_t tmp_16 = -0.0023809523809521678*tmp_12;
      real_t tmp_17 = -0.0023809523809530386*tmp_12;
      real_t tmp_18 = -0.0023809523809526765*tmp_12;
      real_t tmp_19 = -0.001587301587302116*tmp_12;
      real_t tmp_20 = -0.0015873015872993196*tmp_12;
      real_t tmp_21 = -0.0015873015873007629*tmp_12;
      real_t tmp_22 = 0.00039682539682555848*tmp_12;
      real_t tmp_23 = 0.00039682539682542317*tmp_12;
      real_t tmp_24 = -0.0015873015873014776*tmp_12;
      real_t tmp_25 = -0.0015873015873012625*tmp_12;
      real_t tmp_26 = -0.0023809523809528639*tmp_12;
      real_t tmp_27 = -0.002380952380952546*tmp_12;
      real_t tmp_28 = -0.0015873015873028584*tmp_12;
      real_t tmp_29 = 0.00039682539682520113*tmp_12;
      real_t tmp_30 = -0.0015873015873023241*tmp_12;
      real_t tmp_31 = -0.0023809523809513802*tmp_12;
      real_t tmp_32 = -0.0015873015873015262*tmp_12;
      real_t tmp_33 = -0.0023809523809518798*tmp_12;
      real_t tmp_34 = -0.0015873015873017482*tmp_12;
      real_t tmp_35 = -0.0015873015873018037*tmp_12;
      real_t tmp_36 = -0.0023809523809533717*tmp_12;
      real_t tmp_37 = -0.0023809523809519839*tmp_12;
      real_t tmp_38 = 0.0063492063492061046*tmp_12;
      real_t tmp_39 = 0.0063492063492064932*tmp_12;
      real_t tmp_40 = 0.0031746031746029968*tmp_12;
      real_t tmp_41 = 0.012698412698413652*tmp_12;
      real_t tmp_42 = 0.0063492063492061254*tmp_12;
      real_t tmp_43 = 0.0063492063492057715*tmp_12;
      real_t tmp_44 = 0.0063492063492064724*tmp_12;
      real_t tmp_45 = 0.012698412698412737*tmp_12;
      real_t tmp_46 = 0.0063492063492060144*tmp_12;
      real_t tmp_47 = 0.0063492063492060907*tmp_12;
      real_t tmp_48 = 0.0063492063492058826*tmp_12;
      real_t tmp_49 = 0.0063492063492071801*tmp_12;
      real_t a_0_0 = 0.0023809523809515398*tmp_12;
      real_t a_0_1 = tmp_13;
      real_t a_0_2 = tmp_14;
      real_t a_0_3 = tmp_15;
      real_t a_0_4 = tmp_16;
      real_t a_0_5 = tmp_17;
      real_t a_0_6 = tmp_18;
      real_t a_0_7 = tmp_19;
      real_t a_0_8 = tmp_20;
      real_t a_0_9 = tmp_21;
      real_t a_1_0 = tmp_13;
      real_t a_1_1 = 0.0023809523809514843*tmp_12;
      real_t a_1_2 = tmp_22;
      real_t a_1_3 = tmp_23;
      real_t a_1_4 = tmp_16;
      real_t a_1_5 = tmp_24;
      real_t a_1_6 = tmp_25;
      real_t a_1_7 = tmp_26;
      real_t a_1_8 = tmp_27;
      real_t a_1_9 = tmp_28;
      real_t a_2_0 = tmp_14;
      real_t a_2_1 = tmp_22;
      real_t a_2_2 = 0.002380952380951637*tmp_12;
      real_t a_2_3 = tmp_29;
      real_t a_2_4 = tmp_30;
      real_t a_2_5 = tmp_31;
      real_t a_2_6 = tmp_32;
      real_t a_2_7 = tmp_31;
      real_t a_2_8 = tmp_32;
      real_t a_2_9 = tmp_33;
      real_t a_3_0 = tmp_15;
      real_t a_3_1 = tmp_23;
      real_t a_3_2 = tmp_29;
      real_t a_3_3 = 0.0023809523809525945*tmp_12;
      real_t a_3_4 = tmp_34;
      real_t a_3_5 = tmp_35;
      real_t a_3_6 = tmp_36;
      real_t a_3_7 = tmp_35;
      real_t a_3_8 = tmp_36;
      real_t a_3_9 = tmp_37;
      real_t a_4_0 = tmp_16;
      real_t a_4_1 = tmp_16;
      real_t a_4_2 = tmp_30;
      real_t a_4_3 = tmp_34;
      real_t a_4_4 = 0.012698412698413319*tmp_12;
      real_t a_4_5 = tmp_38;
      real_t a_4_6 = tmp_39;
      real_t a_4_7 = tmp_38;
      real_t a_4_8 = tmp_39;
      real_t a_4_9 = tmp_40;
      real_t a_5_0 = tmp_17;
      real_t a_5_1 = tmp_24;
      real_t a_5_2 = tmp_31;
      real_t a_5_3 = tmp_35;
      real_t a_5_4 = tmp_38;
      real_t a_5_5 = tmp_41;
      real_t a_5_6 = tmp_42;
      real_t a_5_7 = tmp_43;
      real_t a_5_8 = tmp_40;
      real_t a_5_9 = tmp_44;
      real_t a_6_0 = tmp_18;
      real_t a_6_1 = tmp_25;
      real_t a_6_2 = tmp_32;
      real_t a_6_3 = tmp_36;
      real_t a_6_4 = tmp_39;
      real_t a_6_5 = tmp_42;
      real_t a_6_6 = tmp_45;
      real_t a_6_7 = tmp_40;
      real_t a_6_8 = tmp_46;
      real_t a_6_9 = tmp_47;
      real_t a_7_0 = tmp_19;
      real_t a_7_1 = tmp_26;
      real_t a_7_2 = tmp_31;
      real_t a_7_3 = tmp_35;
      real_t a_7_4 = tmp_38;
      real_t a_7_5 = tmp_43;
      real_t a_7_6 = tmp_40;
      real_t a_7_7 = tmp_41;
      real_t a_7_8 = tmp_42;
      real_t a_7_9 = tmp_48;
      real_t a_8_0 = tmp_20;
      real_t a_8_1 = tmp_27;
      real_t a_8_2 = tmp_32;
      real_t a_8_3 = tmp_36;
      real_t a_8_4 = tmp_39;
      real_t a_8_5 = tmp_40;
      real_t a_8_6 = tmp_46;
      real_t a_8_7 = tmp_42;
      real_t a_8_8 = tmp_45;
      real_t a_8_9 = tmp_49;
      real_t a_9_0 = tmp_21;
      real_t a_9_1 = tmp_28;
      real_t a_9_2 = tmp_33;
      real_t a_9_3 = tmp_37;
      real_t a_9_4 = tmp_40;
      real_t a_9_5 = tmp_44;
      real_t a_9_6 = tmp_47;
      real_t a_9_7 = tmp_48;
      real_t a_9_8 = tmp_49;
      real_t a_9_9 = 0.012698412698411696*tmp_12;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
      (elMat(0, 4)) = a_0_4;
      (elMat(0, 5)) = a_0_5;
      (elMat(0, 6)) = a_0_6;
      (elMat(0, 7)) = a_0_7;
      (elMat(0, 8)) = a_0_8;
      (elMat(0, 9)) = a_0_9;
      (elMat(1, 0)) = a_1_0;
      (elMat(1, 1)) = a_1_1;
      (elMat(1, 2)) = a_1_2;
      (elMat(1, 3)) = a_1_3;
      (elMat(1, 4)) = a_1_4;
      (elMat(1, 5)) = a_1_5;
      (elMat(1, 6)) = a_1_6;
      (elMat(1, 7)) = a_1_7;
      (elMat(1, 8)) = a_1_8;
      (elMat(1, 9)) = a_1_9;
      (elMat(2, 0)) = a_2_0;
      (elMat(2, 1)) = a_2_1;
      (elMat(2, 2)) = a_2_2;
      (elMat(2, 3)) = a_2_3;
      (elMat(2, 4)) = a_2_4;
      (elMat(2, 5)) = a_2_5;
      (elMat(2, 6)) = a_2_6;
      (elMat(2, 7)) = a_2_7;
      (elMat(2, 8)) = a_2_8;
      (elMat(2, 9)) = a_2_9;
      (elMat(3, 0)) = a_3_0;
      (elMat(3, 1)) = a_3_1;
      (elMat(3, 2)) = a_3_2;
      (elMat(3, 3)) = a_3_3;
      (elMat(3, 4)) = a_3_4;
      (elMat(3, 5)) = a_3_5;
      (elMat(3, 6)) = a_3_6;
      (elMat(3, 7)) = a_3_7;
      (elMat(3, 8)) = a_3_8;
      (elMat(3, 9)) = a_3_9;
      (elMat(4, 0)) = a_4_0;
      (elMat(4, 1)) = a_4_1;
      (elMat(4, 2)) = a_4_2;
      (elMat(4, 3)) = a_4_3;
      (elMat(4, 4)) = a_4_4;
      (elMat(4, 5)) = a_4_5;
      (elMat(4, 6)) = a_4_6;
      (elMat(4, 7)) = a_4_7;
      (elMat(4, 8)) = a_4_8;
      (elMat(4, 9)) = a_4_9;
      (elMat(5, 0)) = a_5_0;
      (elMat(5, 1)) = a_5_1;
      (elMat(5, 2)) = a_5_2;
      (elMat(5, 3)) = a_5_3;
      (elMat(5, 4)) = a_5_4;
      (elMat(5, 5)) = a_5_5;
      (elMat(5, 6)) = a_5_6;
      (elMat(5, 7)) = a_5_7;
      (elMat(5, 8)) = a_5_8;
      (elMat(5, 9)) = a_5_9;
      (elMat(6, 0)) = a_6_0;
      (elMat(6, 1)) = a_6_1;
      (elMat(6, 2)) = a_6_2;
      (elMat(6, 3)) = a_6_3;
      (elMat(6, 4)) = a_6_4;
      (elMat(6, 5)) = a_6_5;
      (elMat(6, 6)) = a_6_6;
      (elMat(6, 7)) = a_6_7;
      (elMat(6, 8)) = a_6_8;
      (elMat(6, 9)) = a_6_9;
      (elMat(7, 0)) = a_7_0;
      (elMat(7, 1)) = a_7_1;
      (elMat(7, 2)) = a_7_2;
      (elMat(7, 3)) = a_7_3;
      (elMat(7, 4)) = a_7_4;
      (elMat(7, 5)) = a_7_5;
      (elMat(7, 6)) = a_7_6;
      (elMat(7, 7)) = a_7_7;
      (elMat(7, 8)) = a_7_8;
      (elMat(7, 9)) = a_7_9;
      (elMat(8, 0)) = a_8_0;
      (elMat(8, 1)) = a_8_1;
      (elMat(8, 2)) = a_8_2;
      (elMat(8, 3)) = a_8_3;
      (elMat(8, 4)) = a_8_4;
      (elMat(8, 5)) = a_8_5;
      (elMat(8, 6)) = a_8_6;
      (elMat(8, 7)) = a_8_7;
      (elMat(8, 8)) = a_8_8;
      (elMat(8, 9)) = a_8_9;
      (elMat(9, 0)) = a_9_0;
      (elMat(9, 1)) = a_9_1;
      (elMat(9, 2)) = a_9_2;
      (elMat(9, 3)) = a_9_3;
      (elMat(9, 4)) = a_9_4;
      (elMat(9, 5)) = a_9_5;
      (elMat(9, 6)) = a_9_6;
      (elMat(9, 7)) = a_9_7;
      (elMat(9, 8)) = a_9_8;
      (elMat(9, 9)) = a_9_9;
   }

   void p2_mass_affine_qe::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 10 >& elMat ) const
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
      real_t tmp_0 = p_affine_0_0*p_affine_1_1;
      real_t tmp_1 = p_affine_0_0*p_affine_1_2;
      real_t tmp_2 = p_affine_2_1*p_affine_3_2;
      real_t tmp_3 = p_affine_0_1*p_affine_1_0;
      real_t tmp_4 = p_affine_0_1*p_affine_1_2;
      real_t tmp_5 = p_affine_2_2*p_affine_3_0;
      real_t tmp_6 = p_affine_0_2*p_affine_1_0;
      real_t tmp_7 = p_affine_0_2*p_affine_1_1;
      real_t tmp_8 = p_affine_2_0*p_affine_3_1;
      real_t tmp_9 = p_affine_2_2*p_affine_3_1;
      real_t tmp_10 = p_affine_2_0*p_affine_3_2;
      real_t tmp_11 = p_affine_2_1*p_affine_3_0;
      real_t tmp_12 = std::abs(p_affine_0_0*tmp_2 - p_affine_0_0*tmp_9 - p_affine_0_1*tmp_10 + p_affine_0_1*tmp_5 - p_affine_0_2*tmp_11 + p_affine_0_2*tmp_8 - p_affine_1_0*tmp_2 + p_affine_1_0*tmp_9 + p_affine_1_1*tmp_10 - p_affine_1_1*tmp_5 + p_affine_1_2*tmp_11 - p_affine_1_2*tmp_8 + p_affine_2_0*tmp_4 - p_affine_2_0*tmp_7 - p_affine_2_1*tmp_1 + p_affine_2_1*tmp_6 + p_affine_2_2*tmp_0 - p_affine_2_2*tmp_3 - p_affine_3_0*tmp_4 + p_affine_3_0*tmp_7 + p_affine_3_1*tmp_1 - p_affine_3_1*tmp_6 - p_affine_3_2*tmp_0 + p_affine_3_2*tmp_3);
      real_t a_0_0 = 0.0023809523809515398*tmp_12;
      real_t a_0_1 = 0.00039682539682555154*tmp_12;
      real_t a_0_2 = 0.0003968253968249201*tmp_12;
      real_t a_0_3 = 0.00039682539682550644*tmp_12;
      real_t a_0_4 = -0.0023809523809521678*tmp_12;
      real_t a_0_5 = -0.0023809523809530386*tmp_12;
      real_t a_0_6 = -0.0023809523809526765*tmp_12;
      real_t a_0_7 = -0.001587301587302116*tmp_12;
      real_t a_0_8 = -0.0015873015872993196*tmp_12;
      real_t a_0_9 = -0.0015873015873007629*tmp_12;
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
      (elMat(0, 4)) = a_0_4;
      (elMat(0, 5)) = a_0_5;
      (elMat(0, 6)) = a_0_6;
      (elMat(0, 7)) = a_0_7;
      (elMat(0, 8)) = a_0_8;
      (elMat(0, 9)) = a_0_9;
   }

} // namespace forms
} // namespace hyteg

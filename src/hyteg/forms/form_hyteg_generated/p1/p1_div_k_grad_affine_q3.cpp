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

#include "p1_div_k_grad_affine_q3.hpp"

namespace hyteg {
namespace forms {

   void p1_div_k_grad_affine_q3::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 3, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Scalar_Variable_Coefficient_2D_k_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id3 = 0;
      Scalar_Variable_Coefficient_2D_k( 0.66639024601470143*p_affine_0_0 + 0.17855872826361643*p_affine_1_0 + 0.1550510257216822*p_affine_2_0, 0.66639024601470143*p_affine_0_1 + 0.17855872826361643*p_affine_1_1 + 0.1550510257216822*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id0 );
      Scalar_Variable_Coefficient_2D_k( 0.28001991549907412*p_affine_0_0 + 0.075031110222608124*p_affine_1_0 + 0.64494897427831777*p_affine_2_0, 0.28001991549907412*p_affine_0_1 + 0.075031110222608124*p_affine_1_1 + 0.64494897427831777*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id1 );
      Scalar_Variable_Coefficient_2D_k( 0.17855872826361638*p_affine_0_0 + 0.66639024601470143*p_affine_1_0 + 0.1550510257216822*p_affine_2_0, 0.17855872826361638*p_affine_0_1 + 0.66639024601470143*p_affine_1_1 + 0.1550510257216822*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id2 );
      Scalar_Variable_Coefficient_2D_k( 0.075031110222608111*p_affine_0_0 + 0.28001991549907407*p_affine_1_0 + 0.64494897427831777*p_affine_2_0, 0.075031110222608111*p_affine_0_1 + 0.28001991549907407*p_affine_1_1 + 0.64494897427831777*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id3 );
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
      real_t tmp_10 = Scalar_Variable_Coefficient_2D_k_out0_id0*tmp_6;
      real_t tmp_11 = Scalar_Variable_Coefficient_2D_k_out0_id0*tmp_8;
      real_t tmp_12 = -tmp_10 - tmp_11;
      real_t tmp_13 = tmp_3*tmp_5;
      real_t tmp_14 = p_affine_0_1 - p_affine_1_1;
      real_t tmp_15 = tmp_14*tmp_5;
      real_t tmp_16 = -tmp_13 - tmp_15;
      real_t tmp_17 = Scalar_Variable_Coefficient_2D_k_out0_id0*tmp_13;
      real_t tmp_18 = Scalar_Variable_Coefficient_2D_k_out0_id0*tmp_15;
      real_t tmp_19 = -tmp_17 - tmp_18;
      real_t tmp_20 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_21 = 0.15902069087198858*tmp_20;
      real_t tmp_22 = Scalar_Variable_Coefficient_2D_k_out0_id1*tmp_6;
      real_t tmp_23 = Scalar_Variable_Coefficient_2D_k_out0_id1*tmp_8;
      real_t tmp_24 = -tmp_22 - tmp_23;
      real_t tmp_25 = Scalar_Variable_Coefficient_2D_k_out0_id1*tmp_13;
      real_t tmp_26 = Scalar_Variable_Coefficient_2D_k_out0_id1*tmp_15;
      real_t tmp_27 = -tmp_25 - tmp_26;
      real_t tmp_28 = 0.090979309128011415*tmp_20;
      real_t tmp_29 = Scalar_Variable_Coefficient_2D_k_out0_id2*tmp_6;
      real_t tmp_30 = Scalar_Variable_Coefficient_2D_k_out0_id2*tmp_8;
      real_t tmp_31 = -tmp_29 - tmp_30;
      real_t tmp_32 = Scalar_Variable_Coefficient_2D_k_out0_id2*tmp_13;
      real_t tmp_33 = Scalar_Variable_Coefficient_2D_k_out0_id2*tmp_15;
      real_t tmp_34 = -tmp_32 - tmp_33;
      real_t tmp_35 = 0.15902069087198858*tmp_20;
      real_t tmp_36 = Scalar_Variable_Coefficient_2D_k_out0_id3*tmp_6;
      real_t tmp_37 = Scalar_Variable_Coefficient_2D_k_out0_id3*tmp_8;
      real_t tmp_38 = -tmp_36 - tmp_37;
      real_t tmp_39 = Scalar_Variable_Coefficient_2D_k_out0_id3*tmp_13;
      real_t tmp_40 = Scalar_Variable_Coefficient_2D_k_out0_id3*tmp_15;
      real_t tmp_41 = -tmp_39 - tmp_40;
      real_t tmp_42 = 0.090979309128011415*tmp_20;
      real_t tmp_43 = (tmp_7*tmp_7);
      real_t tmp_44 = 1.0 / (tmp_4*tmp_4);
      real_t tmp_45 = Scalar_Variable_Coefficient_2D_k_out0_id0*tmp_44;
      real_t tmp_46 = (tmp_3*tmp_3);
      real_t tmp_47 = Scalar_Variable_Coefficient_2D_k_out0_id1*tmp_44;
      real_t tmp_48 = Scalar_Variable_Coefficient_2D_k_out0_id2*tmp_44;
      real_t tmp_49 = Scalar_Variable_Coefficient_2D_k_out0_id3*tmp_44;
      real_t tmp_50 = tmp_1*tmp_7;
      real_t tmp_51 = tmp_14*tmp_3;
      real_t tmp_52 = tmp_21*(tmp_45*tmp_50 + tmp_45*tmp_51) + tmp_28*(tmp_47*tmp_50 + tmp_47*tmp_51) + tmp_35*(tmp_48*tmp_50 + tmp_48*tmp_51) + tmp_42*(tmp_49*tmp_50 + tmp_49*tmp_51);
      real_t tmp_53 = (tmp_1*tmp_1);
      real_t tmp_54 = (tmp_14*tmp_14);
      real_t a_0_0 = tmp_21*(tmp_12*tmp_9 + tmp_16*tmp_19) + tmp_28*(tmp_16*tmp_27 + tmp_24*tmp_9) + tmp_35*(tmp_16*tmp_34 + tmp_31*tmp_9) + tmp_42*(tmp_16*tmp_41 + tmp_38*tmp_9);
      real_t a_0_1 = tmp_21*(tmp_11*tmp_9 + tmp_16*tmp_17) + tmp_28*(tmp_16*tmp_25 + tmp_23*tmp_9) + tmp_35*(tmp_16*tmp_32 + tmp_30*tmp_9) + tmp_42*(tmp_16*tmp_39 + tmp_37*tmp_9);
      real_t a_0_2 = tmp_21*(tmp_10*tmp_9 + tmp_16*tmp_18) + tmp_28*(tmp_16*tmp_26 + tmp_22*tmp_9) + tmp_35*(tmp_16*tmp_33 + tmp_29*tmp_9) + tmp_42*(tmp_16*tmp_40 + tmp_36*tmp_9);
      real_t a_1_0 = tmp_21*(tmp_12*tmp_8 + tmp_13*tmp_19) + tmp_28*(tmp_13*tmp_27 + tmp_24*tmp_8) + tmp_35*(tmp_13*tmp_34 + tmp_31*tmp_8) + tmp_42*(tmp_13*tmp_41 + tmp_38*tmp_8);
      real_t a_1_1 = tmp_21*(tmp_43*tmp_45 + tmp_45*tmp_46) + tmp_28*(tmp_43*tmp_47 + tmp_46*tmp_47) + tmp_35*(tmp_43*tmp_48 + tmp_46*tmp_48) + tmp_42*(tmp_43*tmp_49 + tmp_46*tmp_49);
      real_t a_1_2 = tmp_52;
      real_t a_2_0 = tmp_21*(tmp_12*tmp_6 + tmp_15*tmp_19) + tmp_28*(tmp_15*tmp_27 + tmp_24*tmp_6) + tmp_35*(tmp_15*tmp_34 + tmp_31*tmp_6) + tmp_42*(tmp_15*tmp_41 + tmp_38*tmp_6);
      real_t a_2_1 = tmp_52;
      real_t a_2_2 = tmp_21*(tmp_45*tmp_53 + tmp_45*tmp_54) + tmp_28*(tmp_47*tmp_53 + tmp_47*tmp_54) + tmp_35*(tmp_48*tmp_53 + tmp_48*tmp_54) + tmp_42*(tmp_49*tmp_53 + tmp_49*tmp_54);
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

   void p1_div_k_grad_affine_q3::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Scalar_Variable_Coefficient_2D_k_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id3 = 0;
      Scalar_Variable_Coefficient_2D_k( 0.66639024601470143*p_affine_0_0 + 0.17855872826361643*p_affine_1_0 + 0.1550510257216822*p_affine_2_0, 0.66639024601470143*p_affine_0_1 + 0.17855872826361643*p_affine_1_1 + 0.1550510257216822*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id0 );
      Scalar_Variable_Coefficient_2D_k( 0.28001991549907412*p_affine_0_0 + 0.075031110222608124*p_affine_1_0 + 0.64494897427831777*p_affine_2_0, 0.28001991549907412*p_affine_0_1 + 0.075031110222608124*p_affine_1_1 + 0.64494897427831777*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id1 );
      Scalar_Variable_Coefficient_2D_k( 0.17855872826361638*p_affine_0_0 + 0.66639024601470143*p_affine_1_0 + 0.1550510257216822*p_affine_2_0, 0.17855872826361638*p_affine_0_1 + 0.66639024601470143*p_affine_1_1 + 0.1550510257216822*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id2 );
      Scalar_Variable_Coefficient_2D_k( 0.075031110222608111*p_affine_0_0 + 0.28001991549907407*p_affine_1_0 + 0.64494897427831777*p_affine_2_0, 0.075031110222608111*p_affine_0_1 + 0.28001991549907407*p_affine_1_1 + 0.64494897427831777*p_affine_2_1, &Scalar_Variable_Coefficient_2D_k_out0_id3 );
      real_t tmp_0 = -p_affine_0_0;
      real_t tmp_1 = p_affine_1_0 + tmp_0;
      real_t tmp_2 = -p_affine_0_1;
      real_t tmp_3 = p_affine_2_1 + tmp_2;
      real_t tmp_4 = 1.0 / (tmp_1*tmp_3 - (p_affine_1_1 + tmp_2)*(p_affine_2_0 + tmp_0));
      real_t tmp_5 = tmp_1*tmp_4;
      real_t tmp_6 = tmp_4*(p_affine_0_0 - p_affine_2_0);
      real_t tmp_7 = -tmp_5 - tmp_6;
      real_t tmp_8 = Scalar_Variable_Coefficient_2D_k_out0_id0*tmp_5;
      real_t tmp_9 = Scalar_Variable_Coefficient_2D_k_out0_id0*tmp_6;
      real_t tmp_10 = tmp_3*tmp_4;
      real_t tmp_11 = tmp_4*(p_affine_0_1 - p_affine_1_1);
      real_t tmp_12 = -tmp_10 - tmp_11;
      real_t tmp_13 = Scalar_Variable_Coefficient_2D_k_out0_id0*tmp_10;
      real_t tmp_14 = Scalar_Variable_Coefficient_2D_k_out0_id0*tmp_11;
      real_t tmp_15 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t tmp_16 = 0.15902069087198858*tmp_15;
      real_t tmp_17 = Scalar_Variable_Coefficient_2D_k_out0_id1*tmp_5;
      real_t tmp_18 = Scalar_Variable_Coefficient_2D_k_out0_id1*tmp_6;
      real_t tmp_19 = Scalar_Variable_Coefficient_2D_k_out0_id1*tmp_10;
      real_t tmp_20 = Scalar_Variable_Coefficient_2D_k_out0_id1*tmp_11;
      real_t tmp_21 = 0.090979309128011415*tmp_15;
      real_t tmp_22 = Scalar_Variable_Coefficient_2D_k_out0_id2*tmp_5;
      real_t tmp_23 = Scalar_Variable_Coefficient_2D_k_out0_id2*tmp_6;
      real_t tmp_24 = Scalar_Variable_Coefficient_2D_k_out0_id2*tmp_10;
      real_t tmp_25 = Scalar_Variable_Coefficient_2D_k_out0_id2*tmp_11;
      real_t tmp_26 = 0.15902069087198858*tmp_15;
      real_t tmp_27 = Scalar_Variable_Coefficient_2D_k_out0_id3*tmp_5;
      real_t tmp_28 = Scalar_Variable_Coefficient_2D_k_out0_id3*tmp_6;
      real_t tmp_29 = Scalar_Variable_Coefficient_2D_k_out0_id3*tmp_10;
      real_t tmp_30 = Scalar_Variable_Coefficient_2D_k_out0_id3*tmp_11;
      real_t tmp_31 = 0.090979309128011415*tmp_15;
      real_t a_0_0 = tmp_16*(tmp_12*(-tmp_13 - tmp_14) + tmp_7*(-tmp_8 - tmp_9)) + tmp_21*(tmp_12*(-tmp_19 - tmp_20) + tmp_7*(-tmp_17 - tmp_18)) + tmp_26*(tmp_12*(-tmp_24 - tmp_25) + tmp_7*(-tmp_22 - tmp_23)) + tmp_31*(tmp_12*(-tmp_29 - tmp_30) + tmp_7*(-tmp_27 - tmp_28));
      real_t a_0_1 = tmp_16*(tmp_12*tmp_13 + tmp_7*tmp_9) + tmp_21*(tmp_12*tmp_19 + tmp_18*tmp_7) + tmp_26*(tmp_12*tmp_24 + tmp_23*tmp_7) + tmp_31*(tmp_12*tmp_29 + tmp_28*tmp_7);
      real_t a_0_2 = tmp_16*(tmp_12*tmp_14 + tmp_7*tmp_8) + tmp_21*(tmp_12*tmp_20 + tmp_17*tmp_7) + tmp_26*(tmp_12*tmp_25 + tmp_22*tmp_7) + tmp_31*(tmp_12*tmp_30 + tmp_27*tmp_7);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
   }

   void p1_div_k_grad_affine_q3::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_k_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id4 = 0;
      Scalar_Variable_Coefficient_3D_k( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id0 );
      Scalar_Variable_Coefficient_3D_k( 0.16666666666666674*p_affine_0_0 + 0.16666666666666666*p_affine_1_0 + 0.16666666666666666*p_affine_2_0 + 0.5*p_affine_3_0, 0.16666666666666674*p_affine_0_1 + 0.16666666666666666*p_affine_1_1 + 0.16666666666666666*p_affine_2_1 + 0.5*p_affine_3_1, 0.16666666666666674*p_affine_0_2 + 0.16666666666666666*p_affine_1_2 + 0.16666666666666666*p_affine_2_2 + 0.5*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id1 );
      Scalar_Variable_Coefficient_3D_k( 0.16666666666666671*p_affine_0_0 + 0.16666666666666666*p_affine_1_0 + 0.5*p_affine_2_0 + 0.16666666666666666*p_affine_3_0, 0.16666666666666671*p_affine_0_1 + 0.16666666666666666*p_affine_1_1 + 0.5*p_affine_2_1 + 0.16666666666666666*p_affine_3_1, 0.16666666666666671*p_affine_0_2 + 0.16666666666666666*p_affine_1_2 + 0.5*p_affine_2_2 + 0.16666666666666666*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id2 );
      Scalar_Variable_Coefficient_3D_k( 0.16666666666666671*p_affine_0_0 + 0.5*p_affine_1_0 + 0.16666666666666666*p_affine_2_0 + 0.16666666666666666*p_affine_3_0, 0.16666666666666671*p_affine_0_1 + 0.5*p_affine_1_1 + 0.16666666666666666*p_affine_2_1 + 0.16666666666666666*p_affine_3_1, 0.16666666666666671*p_affine_0_2 + 0.5*p_affine_1_2 + 0.16666666666666666*p_affine_2_2 + 0.16666666666666666*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id3 );
      Scalar_Variable_Coefficient_3D_k( 0.50000000000000011*p_affine_0_0 + 0.16666666666666666*p_affine_1_0 + 0.16666666666666666*p_affine_2_0 + 0.16666666666666666*p_affine_3_0, 0.50000000000000011*p_affine_0_1 + 0.16666666666666666*p_affine_1_1 + 0.16666666666666666*p_affine_2_1 + 0.16666666666666666*p_affine_3_1, 0.50000000000000011*p_affine_0_2 + 0.16666666666666666*p_affine_1_2 + 0.16666666666666666*p_affine_2_2 + 0.16666666666666666*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id4 );
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
      real_t tmp_27 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_21;
      real_t tmp_28 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_23;
      real_t tmp_29 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_25;
      real_t tmp_30 = -tmp_27 - tmp_28 - tmp_29;
      real_t tmp_31 = -tmp_1*tmp_14 + tmp_11*tmp_5;
      real_t tmp_32 = tmp_20*tmp_31;
      real_t tmp_33 = tmp_1*tmp_10 - tmp_11*tmp_15;
      real_t tmp_34 = tmp_20*tmp_33;
      real_t tmp_35 = -tmp_10*tmp_5 + tmp_14*tmp_15;
      real_t tmp_36 = tmp_20*tmp_35;
      real_t tmp_37 = -tmp_32 - tmp_34 - tmp_36;
      real_t tmp_38 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_32;
      real_t tmp_39 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_34;
      real_t tmp_40 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_36;
      real_t tmp_41 = -tmp_38 - tmp_39 - tmp_40;
      real_t tmp_42 = -tmp_11*tmp_3 + tmp_14*tmp_6;
      real_t tmp_43 = tmp_20*tmp_42;
      real_t tmp_44 = -tmp_10*tmp_6 + tmp_11*tmp_12;
      real_t tmp_45 = tmp_20*tmp_44;
      real_t tmp_46 = tmp_10*tmp_3 - tmp_12*tmp_14;
      real_t tmp_47 = tmp_20*tmp_46;
      real_t tmp_48 = -tmp_43 - tmp_45 - tmp_47;
      real_t tmp_49 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_43;
      real_t tmp_50 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_45;
      real_t tmp_51 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_47;
      real_t tmp_52 = -tmp_49 - tmp_50 - tmp_51;
      real_t tmp_53 = p_affine_0_0*p_affine_1_1;
      real_t tmp_54 = p_affine_0_0*p_affine_1_2;
      real_t tmp_55 = p_affine_2_1*p_affine_3_2;
      real_t tmp_56 = p_affine_0_1*p_affine_1_0;
      real_t tmp_57 = p_affine_0_1*p_affine_1_2;
      real_t tmp_58 = p_affine_2_2*p_affine_3_0;
      real_t tmp_59 = p_affine_0_2*p_affine_1_0;
      real_t tmp_60 = p_affine_0_2*p_affine_1_1;
      real_t tmp_61 = p_affine_2_0*p_affine_3_1;
      real_t tmp_62 = p_affine_2_2*p_affine_3_1;
      real_t tmp_63 = p_affine_2_0*p_affine_3_2;
      real_t tmp_64 = p_affine_2_1*p_affine_3_0;
      real_t tmp_65 = std::abs(p_affine_0_0*tmp_55 - p_affine_0_0*tmp_62 + p_affine_0_1*tmp_58 - p_affine_0_1*tmp_63 + p_affine_0_2*tmp_61 - p_affine_0_2*tmp_64 - p_affine_1_0*tmp_55 + p_affine_1_0*tmp_62 - p_affine_1_1*tmp_58 + p_affine_1_1*tmp_63 - p_affine_1_2*tmp_61 + p_affine_1_2*tmp_64 + p_affine_2_0*tmp_57 - p_affine_2_0*tmp_60 - p_affine_2_1*tmp_54 + p_affine_2_1*tmp_59 + p_affine_2_2*tmp_53 - p_affine_2_2*tmp_56 - p_affine_3_0*tmp_57 + p_affine_3_0*tmp_60 + p_affine_3_1*tmp_54 - p_affine_3_1*tmp_59 - p_affine_3_2*tmp_53 + p_affine_3_2*tmp_56);
      real_t tmp_66 = -0.1333333333333333*tmp_65;
      real_t tmp_67 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_21;
      real_t tmp_68 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_23;
      real_t tmp_69 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_25;
      real_t tmp_70 = -tmp_67 - tmp_68 - tmp_69;
      real_t tmp_71 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_32;
      real_t tmp_72 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_34;
      real_t tmp_73 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_36;
      real_t tmp_74 = -tmp_71 - tmp_72 - tmp_73;
      real_t tmp_75 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_43;
      real_t tmp_76 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_45;
      real_t tmp_77 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_47;
      real_t tmp_78 = -tmp_75 - tmp_76 - tmp_77;
      real_t tmp_79 = 0.074999999999999983*tmp_65;
      real_t tmp_80 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_21;
      real_t tmp_81 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_23;
      real_t tmp_82 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_25;
      real_t tmp_83 = -tmp_80 - tmp_81 - tmp_82;
      real_t tmp_84 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_32;
      real_t tmp_85 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_34;
      real_t tmp_86 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_36;
      real_t tmp_87 = -tmp_84 - tmp_85 - tmp_86;
      real_t tmp_88 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_43;
      real_t tmp_89 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_45;
      real_t tmp_90 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_47;
      real_t tmp_91 = -tmp_88 - tmp_89 - tmp_90;
      real_t tmp_92 = 0.074999999999999983*tmp_65;
      real_t tmp_93 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_21;
      real_t tmp_94 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_23;
      real_t tmp_95 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_25;
      real_t tmp_96 = -tmp_93 - tmp_94 - tmp_95;
      real_t tmp_97 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_32;
      real_t tmp_98 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_34;
      real_t tmp_99 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_36;
      real_t tmp_100 = -tmp_97 - tmp_98 - tmp_99;
      real_t tmp_101 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_43;
      real_t tmp_102 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_45;
      real_t tmp_103 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_47;
      real_t tmp_104 = -tmp_101 - tmp_102 - tmp_103;
      real_t tmp_105 = 0.074999999999999983*tmp_65;
      real_t tmp_106 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_21;
      real_t tmp_107 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_23;
      real_t tmp_108 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_25;
      real_t tmp_109 = -tmp_106 - tmp_107 - tmp_108;
      real_t tmp_110 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_32;
      real_t tmp_111 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_34;
      real_t tmp_112 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_36;
      real_t tmp_113 = -tmp_110 - tmp_111 - tmp_112;
      real_t tmp_114 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_43;
      real_t tmp_115 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_45;
      real_t tmp_116 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_47;
      real_t tmp_117 = -tmp_114 - tmp_115 - tmp_116;
      real_t tmp_118 = 0.074999999999999983*tmp_65;
      real_t tmp_119 = (tmp_24*tmp_24);
      real_t tmp_120 = 1.0 / (tmp_19*tmp_19);
      real_t tmp_121 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_120;
      real_t tmp_122 = (tmp_35*tmp_35);
      real_t tmp_123 = (tmp_46*tmp_46);
      real_t tmp_124 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_120;
      real_t tmp_125 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_120;
      real_t tmp_126 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_120;
      real_t tmp_127 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_120;
      real_t tmp_128 = tmp_22*tmp_24;
      real_t tmp_129 = tmp_33*tmp_35;
      real_t tmp_130 = tmp_44*tmp_46;
      real_t tmp_131 = tmp_105*(tmp_126*tmp_128 + tmp_126*tmp_129 + tmp_126*tmp_130) + tmp_118*(tmp_127*tmp_128 + tmp_127*tmp_129 + tmp_127*tmp_130) + tmp_66*(tmp_121*tmp_128 + tmp_121*tmp_129 + tmp_121*tmp_130) + tmp_79*(tmp_124*tmp_128 + tmp_124*tmp_129 + tmp_124*tmp_130) + tmp_92*(tmp_125*tmp_128 + tmp_125*tmp_129 + tmp_125*tmp_130);
      real_t tmp_132 = tmp_24*tmp_8;
      real_t tmp_133 = tmp_31*tmp_35;
      real_t tmp_134 = tmp_42*tmp_46;
      real_t tmp_135 = tmp_105*(tmp_126*tmp_132 + tmp_126*tmp_133 + tmp_126*tmp_134) + tmp_118*(tmp_127*tmp_132 + tmp_127*tmp_133 + tmp_127*tmp_134) + tmp_66*(tmp_121*tmp_132 + tmp_121*tmp_133 + tmp_121*tmp_134) + tmp_79*(tmp_124*tmp_132 + tmp_124*tmp_133 + tmp_124*tmp_134) + tmp_92*(tmp_125*tmp_132 + tmp_125*tmp_133 + tmp_125*tmp_134);
      real_t tmp_136 = (tmp_22*tmp_22);
      real_t tmp_137 = (tmp_33*tmp_33);
      real_t tmp_138 = (tmp_44*tmp_44);
      real_t tmp_139 = tmp_22*tmp_8;
      real_t tmp_140 = tmp_31*tmp_33;
      real_t tmp_141 = tmp_42*tmp_44;
      real_t tmp_142 = tmp_105*(tmp_126*tmp_139 + tmp_126*tmp_140 + tmp_126*tmp_141) + tmp_118*(tmp_127*tmp_139 + tmp_127*tmp_140 + tmp_127*tmp_141) + tmp_66*(tmp_121*tmp_139 + tmp_121*tmp_140 + tmp_121*tmp_141) + tmp_79*(tmp_124*tmp_139 + tmp_124*tmp_140 + tmp_124*tmp_141) + tmp_92*(tmp_125*tmp_139 + tmp_125*tmp_140 + tmp_125*tmp_141);
      real_t tmp_143 = (tmp_8*tmp_8);
      real_t tmp_144 = (tmp_31*tmp_31);
      real_t tmp_145 = (tmp_42*tmp_42);
      real_t a_0_0 = tmp_105*(tmp_100*tmp_37 + tmp_104*tmp_48 + tmp_26*tmp_96) + tmp_118*(tmp_109*tmp_26 + tmp_113*tmp_37 + tmp_117*tmp_48) + tmp_66*(tmp_26*tmp_30 + tmp_37*tmp_41 + tmp_48*tmp_52) + tmp_79*(tmp_26*tmp_70 + tmp_37*tmp_74 + tmp_48*tmp_78) + tmp_92*(tmp_26*tmp_83 + tmp_37*tmp_87 + tmp_48*tmp_91);
      real_t a_0_1 = tmp_105*(tmp_103*tmp_48 + tmp_26*tmp_95 + tmp_37*tmp_99) + tmp_118*(tmp_108*tmp_26 + tmp_112*tmp_37 + tmp_116*tmp_48) + tmp_66*(tmp_26*tmp_29 + tmp_37*tmp_40 + tmp_48*tmp_51) + tmp_79*(tmp_26*tmp_69 + tmp_37*tmp_73 + tmp_48*tmp_77) + tmp_92*(tmp_26*tmp_82 + tmp_37*tmp_86 + tmp_48*tmp_90);
      real_t a_0_2 = tmp_105*(tmp_102*tmp_48 + tmp_26*tmp_94 + tmp_37*tmp_98) + tmp_118*(tmp_107*tmp_26 + tmp_111*tmp_37 + tmp_115*tmp_48) + tmp_66*(tmp_26*tmp_28 + tmp_37*tmp_39 + tmp_48*tmp_50) + tmp_79*(tmp_26*tmp_68 + tmp_37*tmp_72 + tmp_48*tmp_76) + tmp_92*(tmp_26*tmp_81 + tmp_37*tmp_85 + tmp_48*tmp_89);
      real_t a_0_3 = tmp_105*(tmp_101*tmp_48 + tmp_26*tmp_93 + tmp_37*tmp_97) + tmp_118*(tmp_106*tmp_26 + tmp_110*tmp_37 + tmp_114*tmp_48) + tmp_66*(tmp_26*tmp_27 + tmp_37*tmp_38 + tmp_48*tmp_49) + tmp_79*(tmp_26*tmp_67 + tmp_37*tmp_71 + tmp_48*tmp_75) + tmp_92*(tmp_26*tmp_80 + tmp_37*tmp_84 + tmp_48*tmp_88);
      real_t a_1_0 = tmp_105*(tmp_100*tmp_36 + tmp_104*tmp_47 + tmp_25*tmp_96) + tmp_118*(tmp_109*tmp_25 + tmp_113*tmp_36 + tmp_117*tmp_47) + tmp_66*(tmp_25*tmp_30 + tmp_36*tmp_41 + tmp_47*tmp_52) + tmp_79*(tmp_25*tmp_70 + tmp_36*tmp_74 + tmp_47*tmp_78) + tmp_92*(tmp_25*tmp_83 + tmp_36*tmp_87 + tmp_47*tmp_91);
      real_t a_1_1 = tmp_105*(tmp_119*tmp_126 + tmp_122*tmp_126 + tmp_123*tmp_126) + tmp_118*(tmp_119*tmp_127 + tmp_122*tmp_127 + tmp_123*tmp_127) + tmp_66*(tmp_119*tmp_121 + tmp_121*tmp_122 + tmp_121*tmp_123) + tmp_79*(tmp_119*tmp_124 + tmp_122*tmp_124 + tmp_123*tmp_124) + tmp_92*(tmp_119*tmp_125 + tmp_122*tmp_125 + tmp_123*tmp_125);
      real_t a_1_2 = tmp_131;
      real_t a_1_3 = tmp_135;
      real_t a_2_0 = tmp_105*(tmp_100*tmp_34 + tmp_104*tmp_45 + tmp_23*tmp_96) + tmp_118*(tmp_109*tmp_23 + tmp_113*tmp_34 + tmp_117*tmp_45) + tmp_66*(tmp_23*tmp_30 + tmp_34*tmp_41 + tmp_45*tmp_52) + tmp_79*(tmp_23*tmp_70 + tmp_34*tmp_74 + tmp_45*tmp_78) + tmp_92*(tmp_23*tmp_83 + tmp_34*tmp_87 + tmp_45*tmp_91);
      real_t a_2_1 = tmp_131;
      real_t a_2_2 = tmp_105*(tmp_126*tmp_136 + tmp_126*tmp_137 + tmp_126*tmp_138) + tmp_118*(tmp_127*tmp_136 + tmp_127*tmp_137 + tmp_127*tmp_138) + tmp_66*(tmp_121*tmp_136 + tmp_121*tmp_137 + tmp_121*tmp_138) + tmp_79*(tmp_124*tmp_136 + tmp_124*tmp_137 + tmp_124*tmp_138) + tmp_92*(tmp_125*tmp_136 + tmp_125*tmp_137 + tmp_125*tmp_138);
      real_t a_2_3 = tmp_142;
      real_t a_3_0 = tmp_105*(tmp_100*tmp_32 + tmp_104*tmp_43 + tmp_21*tmp_96) + tmp_118*(tmp_109*tmp_21 + tmp_113*tmp_32 + tmp_117*tmp_43) + tmp_66*(tmp_21*tmp_30 + tmp_32*tmp_41 + tmp_43*tmp_52) + tmp_79*(tmp_21*tmp_70 + tmp_32*tmp_74 + tmp_43*tmp_78) + tmp_92*(tmp_21*tmp_83 + tmp_32*tmp_87 + tmp_43*tmp_91);
      real_t a_3_1 = tmp_135;
      real_t a_3_2 = tmp_142;
      real_t a_3_3 = tmp_105*(tmp_126*tmp_143 + tmp_126*tmp_144 + tmp_126*tmp_145) + tmp_118*(tmp_127*tmp_143 + tmp_127*tmp_144 + tmp_127*tmp_145) + tmp_66*(tmp_121*tmp_143 + tmp_121*tmp_144 + tmp_121*tmp_145) + tmp_79*(tmp_124*tmp_143 + tmp_124*tmp_144 + tmp_124*tmp_145) + tmp_92*(tmp_125*tmp_143 + tmp_125*tmp_144 + tmp_125*tmp_145);
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

   void p1_div_k_grad_affine_q3::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t Scalar_Variable_Coefficient_3D_k_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id4 = 0;
      Scalar_Variable_Coefficient_3D_k( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id0 );
      Scalar_Variable_Coefficient_3D_k( 0.16666666666666674*p_affine_0_0 + 0.16666666666666666*p_affine_1_0 + 0.16666666666666666*p_affine_2_0 + 0.5*p_affine_3_0, 0.16666666666666674*p_affine_0_1 + 0.16666666666666666*p_affine_1_1 + 0.16666666666666666*p_affine_2_1 + 0.5*p_affine_3_1, 0.16666666666666674*p_affine_0_2 + 0.16666666666666666*p_affine_1_2 + 0.16666666666666666*p_affine_2_2 + 0.5*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id1 );
      Scalar_Variable_Coefficient_3D_k( 0.16666666666666671*p_affine_0_0 + 0.16666666666666666*p_affine_1_0 + 0.5*p_affine_2_0 + 0.16666666666666666*p_affine_3_0, 0.16666666666666671*p_affine_0_1 + 0.16666666666666666*p_affine_1_1 + 0.5*p_affine_2_1 + 0.16666666666666666*p_affine_3_1, 0.16666666666666671*p_affine_0_2 + 0.16666666666666666*p_affine_1_2 + 0.5*p_affine_2_2 + 0.16666666666666666*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id2 );
      Scalar_Variable_Coefficient_3D_k( 0.16666666666666671*p_affine_0_0 + 0.5*p_affine_1_0 + 0.16666666666666666*p_affine_2_0 + 0.16666666666666666*p_affine_3_0, 0.16666666666666671*p_affine_0_1 + 0.5*p_affine_1_1 + 0.16666666666666666*p_affine_2_1 + 0.16666666666666666*p_affine_3_1, 0.16666666666666671*p_affine_0_2 + 0.5*p_affine_1_2 + 0.16666666666666666*p_affine_2_2 + 0.16666666666666666*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id3 );
      Scalar_Variable_Coefficient_3D_k( 0.50000000000000011*p_affine_0_0 + 0.16666666666666666*p_affine_1_0 + 0.16666666666666666*p_affine_2_0 + 0.16666666666666666*p_affine_3_0, 0.50000000000000011*p_affine_0_1 + 0.16666666666666666*p_affine_1_1 + 0.16666666666666666*p_affine_2_1 + 0.16666666666666666*p_affine_3_1, 0.50000000000000011*p_affine_0_2 + 0.16666666666666666*p_affine_1_2 + 0.16666666666666666*p_affine_2_2 + 0.16666666666666666*p_affine_3_2, &Scalar_Variable_Coefficient_3D_k_out0_id4 );
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
      real_t tmp_23 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_19;
      real_t tmp_24 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_20;
      real_t tmp_25 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_21;
      real_t tmp_26 = tmp_18*(-tmp_1*tmp_13 + tmp_10*tmp_5);
      real_t tmp_27 = tmp_18*(tmp_1*tmp_9 - tmp_10*tmp_14);
      real_t tmp_28 = tmp_18*(tmp_13*tmp_14 - tmp_5*tmp_9);
      real_t tmp_29 = -tmp_26 - tmp_27 - tmp_28;
      real_t tmp_30 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_26;
      real_t tmp_31 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_27;
      real_t tmp_32 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_28;
      real_t tmp_33 = tmp_18*(-tmp_10*tmp_3 + tmp_13*tmp_6);
      real_t tmp_34 = tmp_18*(tmp_10*tmp_11 - tmp_6*tmp_9);
      real_t tmp_35 = tmp_18*(-tmp_11*tmp_13 + tmp_3*tmp_9);
      real_t tmp_36 = -tmp_33 - tmp_34 - tmp_35;
      real_t tmp_37 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_33;
      real_t tmp_38 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_34;
      real_t tmp_39 = Scalar_Variable_Coefficient_3D_k_out0_id0*tmp_35;
      real_t tmp_40 = p_affine_0_0*p_affine_1_1;
      real_t tmp_41 = p_affine_0_0*p_affine_1_2;
      real_t tmp_42 = p_affine_2_1*p_affine_3_2;
      real_t tmp_43 = p_affine_0_1*p_affine_1_0;
      real_t tmp_44 = p_affine_0_1*p_affine_1_2;
      real_t tmp_45 = p_affine_2_2*p_affine_3_0;
      real_t tmp_46 = p_affine_0_2*p_affine_1_0;
      real_t tmp_47 = p_affine_0_2*p_affine_1_1;
      real_t tmp_48 = p_affine_2_0*p_affine_3_1;
      real_t tmp_49 = p_affine_2_2*p_affine_3_1;
      real_t tmp_50 = p_affine_2_0*p_affine_3_2;
      real_t tmp_51 = p_affine_2_1*p_affine_3_0;
      real_t tmp_52 = std::abs(p_affine_0_0*tmp_42 - p_affine_0_0*tmp_49 + p_affine_0_1*tmp_45 - p_affine_0_1*tmp_50 + p_affine_0_2*tmp_48 - p_affine_0_2*tmp_51 - p_affine_1_0*tmp_42 + p_affine_1_0*tmp_49 - p_affine_1_1*tmp_45 + p_affine_1_1*tmp_50 - p_affine_1_2*tmp_48 + p_affine_1_2*tmp_51 + p_affine_2_0*tmp_44 - p_affine_2_0*tmp_47 - p_affine_2_1*tmp_41 + p_affine_2_1*tmp_46 + p_affine_2_2*tmp_40 - p_affine_2_2*tmp_43 - p_affine_3_0*tmp_44 + p_affine_3_0*tmp_47 + p_affine_3_1*tmp_41 - p_affine_3_1*tmp_46 - p_affine_3_2*tmp_40 + p_affine_3_2*tmp_43);
      real_t tmp_53 = -0.1333333333333333*tmp_52;
      real_t tmp_54 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_19;
      real_t tmp_55 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_20;
      real_t tmp_56 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_21;
      real_t tmp_57 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_26;
      real_t tmp_58 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_27;
      real_t tmp_59 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_28;
      real_t tmp_60 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_33;
      real_t tmp_61 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_34;
      real_t tmp_62 = Scalar_Variable_Coefficient_3D_k_out0_id1*tmp_35;
      real_t tmp_63 = 0.074999999999999983*tmp_52;
      real_t tmp_64 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_19;
      real_t tmp_65 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_20;
      real_t tmp_66 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_21;
      real_t tmp_67 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_26;
      real_t tmp_68 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_27;
      real_t tmp_69 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_28;
      real_t tmp_70 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_33;
      real_t tmp_71 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_34;
      real_t tmp_72 = Scalar_Variable_Coefficient_3D_k_out0_id2*tmp_35;
      real_t tmp_73 = 0.074999999999999983*tmp_52;
      real_t tmp_74 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_19;
      real_t tmp_75 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_20;
      real_t tmp_76 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_21;
      real_t tmp_77 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_26;
      real_t tmp_78 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_27;
      real_t tmp_79 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_28;
      real_t tmp_80 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_33;
      real_t tmp_81 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_34;
      real_t tmp_82 = Scalar_Variable_Coefficient_3D_k_out0_id3*tmp_35;
      real_t tmp_83 = 0.074999999999999983*tmp_52;
      real_t tmp_84 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_19;
      real_t tmp_85 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_20;
      real_t tmp_86 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_21;
      real_t tmp_87 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_26;
      real_t tmp_88 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_27;
      real_t tmp_89 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_28;
      real_t tmp_90 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_33;
      real_t tmp_91 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_34;
      real_t tmp_92 = Scalar_Variable_Coefficient_3D_k_out0_id4*tmp_35;
      real_t tmp_93 = 0.074999999999999983*tmp_52;
      real_t a_0_0 = tmp_53*(tmp_22*(-tmp_23 - tmp_24 - tmp_25) + tmp_29*(-tmp_30 - tmp_31 - tmp_32) + tmp_36*(-tmp_37 - tmp_38 - tmp_39)) + tmp_63*(tmp_22*(-tmp_54 - tmp_55 - tmp_56) + tmp_29*(-tmp_57 - tmp_58 - tmp_59) + tmp_36*(-tmp_60 - tmp_61 - tmp_62)) + tmp_73*(tmp_22*(-tmp_64 - tmp_65 - tmp_66) + tmp_29*(-tmp_67 - tmp_68 - tmp_69) + tmp_36*(-tmp_70 - tmp_71 - tmp_72)) + tmp_83*(tmp_22*(-tmp_74 - tmp_75 - tmp_76) + tmp_29*(-tmp_77 - tmp_78 - tmp_79) + tmp_36*(-tmp_80 - tmp_81 - tmp_82)) + tmp_93*(tmp_22*(-tmp_84 - tmp_85 - tmp_86) + tmp_29*(-tmp_87 - tmp_88 - tmp_89) + tmp_36*(-tmp_90 - tmp_91 - tmp_92));
      real_t a_0_1 = tmp_53*(tmp_22*tmp_25 + tmp_29*tmp_32 + tmp_36*tmp_39) + tmp_63*(tmp_22*tmp_56 + tmp_29*tmp_59 + tmp_36*tmp_62) + tmp_73*(tmp_22*tmp_66 + tmp_29*tmp_69 + tmp_36*tmp_72) + tmp_83*(tmp_22*tmp_76 + tmp_29*tmp_79 + tmp_36*tmp_82) + tmp_93*(tmp_22*tmp_86 + tmp_29*tmp_89 + tmp_36*tmp_92);
      real_t a_0_2 = tmp_53*(tmp_22*tmp_24 + tmp_29*tmp_31 + tmp_36*tmp_38) + tmp_63*(tmp_22*tmp_55 + tmp_29*tmp_58 + tmp_36*tmp_61) + tmp_73*(tmp_22*tmp_65 + tmp_29*tmp_68 + tmp_36*tmp_71) + tmp_83*(tmp_22*tmp_75 + tmp_29*tmp_78 + tmp_36*tmp_81) + tmp_93*(tmp_22*tmp_85 + tmp_29*tmp_88 + tmp_36*tmp_91);
      real_t a_0_3 = tmp_53*(tmp_22*tmp_23 + tmp_29*tmp_30 + tmp_36*tmp_37) + tmp_63*(tmp_22*tmp_54 + tmp_29*tmp_57 + tmp_36*tmp_60) + tmp_73*(tmp_22*tmp_64 + tmp_29*tmp_67 + tmp_36*tmp_70) + tmp_83*(tmp_22*tmp_74 + tmp_29*tmp_77 + tmp_36*tmp_80) + tmp_93*(tmp_22*tmp_84 + tmp_29*tmp_87 + tmp_36*tmp_90);
      (elMat(0, 0)) = a_0_0;
      (elMat(0, 1)) = a_0_1;
      (elMat(0, 2)) = a_0_2;
      (elMat(0, 3)) = a_0_3;
   }

   void p1_div_k_grad_affine_q3::Scalar_Variable_Coefficient_2D_k( real_t in_0, real_t in_1, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_k( Point3D( {in_0, in_1, 0} ) );
   }

   void p1_div_k_grad_affine_q3::Scalar_Variable_Coefficient_3D_k( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_k( Point3D( {in_0, in_1, in_2} ) );
   }

} // namespace forms
} // namespace hyteg

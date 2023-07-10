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

#include "n1e1_curl_curl_affine_q0.hpp"

namespace hyteg {
namespace forms {

   void n1e1_curl_curl_affine_q0::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 6, 6 >& elMat ) const
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
      real_t tmp_0 = 2*p_affine_0_0;
      real_t tmp_1 = -tmp_0;
      real_t tmp_2 = 2*p_affine_1_0;
      real_t tmp_3 = tmp_1 + tmp_2;
      real_t tmp_4 = 2*p_affine_0_1;
      real_t tmp_5 = -tmp_4;
      real_t tmp_6 = 2*p_affine_1_1;
      real_t tmp_7 = tmp_5 + tmp_6;
      real_t tmp_8 = 2*p_affine_0_2;
      real_t tmp_9 = -tmp_8;
      real_t tmp_10 = 2*p_affine_1_2;
      real_t tmp_11 = tmp_10 + tmp_9;
      real_t tmp_12 = p_affine_0_0*p_affine_1_1;
      real_t tmp_13 = p_affine_0_0*p_affine_1_2;
      real_t tmp_14 = p_affine_2_1*p_affine_3_2;
      real_t tmp_15 = p_affine_0_1*p_affine_1_0;
      real_t tmp_16 = p_affine_0_1*p_affine_1_2;
      real_t tmp_17 = p_affine_2_2*p_affine_3_0;
      real_t tmp_18 = p_affine_0_2*p_affine_1_0;
      real_t tmp_19 = p_affine_0_2*p_affine_1_1;
      real_t tmp_20 = p_affine_2_0*p_affine_3_1;
      real_t tmp_21 = p_affine_2_2*p_affine_3_1;
      real_t tmp_22 = p_affine_2_0*p_affine_3_2;
      real_t tmp_23 = p_affine_2_1*p_affine_3_0;
      real_t tmp_24 = 0.16666666666666663/std::abs(p_affine_0_0*tmp_14 - p_affine_0_0*tmp_21 + p_affine_0_1*tmp_17 - p_affine_0_1*tmp_22 + p_affine_0_2*tmp_20 - p_affine_0_2*tmp_23 - p_affine_1_0*tmp_14 + p_affine_1_0*tmp_21 - p_affine_1_1*tmp_17 + p_affine_1_1*tmp_22 - p_affine_1_2*tmp_20 + p_affine_1_2*tmp_23 + p_affine_2_0*tmp_16 - p_affine_2_0*tmp_19 - p_affine_2_1*tmp_13 + p_affine_2_1*tmp_18 + p_affine_2_2*tmp_12 - p_affine_2_2*tmp_15 - p_affine_3_0*tmp_16 + p_affine_3_0*tmp_19 + p_affine_3_1*tmp_13 - p_affine_3_1*tmp_18 - p_affine_3_2*tmp_12 + p_affine_3_2*tmp_15);
      real_t tmp_25 = 2*p_affine_2_0;
      real_t tmp_26 = -tmp_25;
      real_t tmp_27 = tmp_0 + tmp_26;
      real_t tmp_28 = 2*p_affine_2_1;
      real_t tmp_29 = -tmp_28;
      real_t tmp_30 = tmp_29 + tmp_4;
      real_t tmp_31 = 2*p_affine_2_2;
      real_t tmp_32 = -tmp_31;
      real_t tmp_33 = tmp_32 + tmp_8;
      real_t tmp_34 = tmp_24*(tmp_11*tmp_33 + tmp_27*tmp_3 + tmp_30*tmp_7);
      real_t tmp_35 = 2*p_affine_3_0;
      real_t tmp_36 = tmp_1 + tmp_35;
      real_t tmp_37 = 2*p_affine_3_1;
      real_t tmp_38 = tmp_37 + tmp_5;
      real_t tmp_39 = 2*p_affine_3_2;
      real_t tmp_40 = tmp_39 + tmp_9;
      real_t tmp_41 = tmp_24*(tmp_11*tmp_40 + tmp_3*tmp_36 + tmp_38*tmp_7);
      real_t tmp_42 = -tmp_2 + tmp_25;
      real_t tmp_43 = tmp_28 - tmp_6;
      real_t tmp_44 = -tmp_10 + tmp_31;
      real_t tmp_45 = tmp_24*(tmp_11*tmp_44 + tmp_3*tmp_42 + tmp_43*tmp_7);
      real_t tmp_46 = tmp_2 - tmp_35;
      real_t tmp_47 = -tmp_37 + tmp_6;
      real_t tmp_48 = tmp_10 - tmp_39;
      real_t tmp_49 = tmp_24*(tmp_11*tmp_48 + tmp_3*tmp_46 + tmp_47*tmp_7);
      real_t tmp_50 = tmp_26 + tmp_35;
      real_t tmp_51 = tmp_29 + tmp_37;
      real_t tmp_52 = tmp_32 + tmp_39;
      real_t tmp_53 = tmp_24*(tmp_11*tmp_52 + tmp_3*tmp_50 + tmp_51*tmp_7);
      real_t tmp_54 = tmp_24*(tmp_27*tmp_36 + tmp_30*tmp_38 + tmp_33*tmp_40);
      real_t tmp_55 = tmp_24*(tmp_27*tmp_42 + tmp_30*tmp_43 + tmp_33*tmp_44);
      real_t tmp_56 = tmp_24*(tmp_27*tmp_46 + tmp_30*tmp_47 + tmp_33*tmp_48);
      real_t tmp_57 = tmp_24*(tmp_27*tmp_50 + tmp_30*tmp_51 + tmp_33*tmp_52);
      real_t tmp_58 = tmp_24*(tmp_36*tmp_42 + tmp_38*tmp_43 + tmp_40*tmp_44);
      real_t tmp_59 = tmp_24*(tmp_36*tmp_46 + tmp_38*tmp_47 + tmp_40*tmp_48);
      real_t tmp_60 = tmp_24*(tmp_36*tmp_50 + tmp_38*tmp_51 + tmp_40*tmp_52);
      real_t tmp_61 = tmp_24*(tmp_42*tmp_46 + tmp_43*tmp_47 + tmp_44*tmp_48);
      real_t tmp_62 = tmp_24*(tmp_42*tmp_50 + tmp_43*tmp_51 + tmp_44*tmp_52);
      real_t tmp_63 = tmp_24*(tmp_46*tmp_50 + tmp_47*tmp_51 + tmp_48*tmp_52);
      real_t a_0_0 = tmp_24*((tmp_11*tmp_11) + (tmp_3*tmp_3) + (tmp_7*tmp_7));
      real_t a_0_1 = tmp_34;
      real_t a_0_2 = tmp_41;
      real_t a_0_3 = tmp_45;
      real_t a_0_4 = tmp_49;
      real_t a_0_5 = tmp_53;
      real_t a_1_0 = tmp_34;
      real_t a_1_1 = tmp_24*((tmp_27*tmp_27) + (tmp_30*tmp_30) + (tmp_33*tmp_33));
      real_t a_1_2 = tmp_54;
      real_t a_1_3 = tmp_55;
      real_t a_1_4 = tmp_56;
      real_t a_1_5 = tmp_57;
      real_t a_2_0 = tmp_41;
      real_t a_2_1 = tmp_54;
      real_t a_2_2 = tmp_24*((tmp_36*tmp_36) + (tmp_38*tmp_38) + (tmp_40*tmp_40));
      real_t a_2_3 = tmp_58;
      real_t a_2_4 = tmp_59;
      real_t a_2_5 = tmp_60;
      real_t a_3_0 = tmp_45;
      real_t a_3_1 = tmp_55;
      real_t a_3_2 = tmp_58;
      real_t a_3_3 = tmp_24*((tmp_42*tmp_42) + (tmp_43*tmp_43) + (tmp_44*tmp_44));
      real_t a_3_4 = tmp_61;
      real_t a_3_5 = tmp_62;
      real_t a_4_0 = tmp_49;
      real_t a_4_1 = tmp_56;
      real_t a_4_2 = tmp_59;
      real_t a_4_3 = tmp_61;
      real_t a_4_4 = tmp_24*((tmp_46*tmp_46) + (tmp_47*tmp_47) + (tmp_48*tmp_48));
      real_t a_4_5 = tmp_63;
      real_t a_5_0 = tmp_53;
      real_t a_5_1 = tmp_57;
      real_t a_5_2 = tmp_60;
      real_t a_5_3 = tmp_62;
      real_t a_5_4 = tmp_63;
      real_t a_5_5 = tmp_24*((tmp_50*tmp_50) + (tmp_51*tmp_51) + (tmp_52*tmp_52));
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

} // namespace forms
} // namespace hyteg

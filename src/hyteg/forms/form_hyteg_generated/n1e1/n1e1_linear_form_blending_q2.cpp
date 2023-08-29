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

#include "n1e1_linear_form_blending_q2.hpp"

namespace hyteg {
namespace forms {

   void n1e1_linear_form_blending_q2::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 6, 6 >& elMat ) const
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
      real_t Blending_F_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id2 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id2 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id2 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id3 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id3 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id3 = 0;
      real_t Vector_Variable_Coefficient_3D_k_out0_id0 = 0;
      real_t Vector_Variable_Coefficient_3D_k_out1_id0 = 0;
      real_t Vector_Variable_Coefficient_3D_k_out2_id0 = 0;
      real_t Vector_Variable_Coefficient_3D_k_out0_id1 = 0;
      real_t Vector_Variable_Coefficient_3D_k_out1_id1 = 0;
      real_t Vector_Variable_Coefficient_3D_k_out2_id1 = 0;
      real_t Vector_Variable_Coefficient_3D_k_out0_id2 = 0;
      real_t Vector_Variable_Coefficient_3D_k_out1_id2 = 0;
      real_t Vector_Variable_Coefficient_3D_k_out2_id2 = 0;
      real_t Vector_Variable_Coefficient_3D_k_out0_id3 = 0;
      real_t Vector_Variable_Coefficient_3D_k_out1_id3 = 0;
      real_t Vector_Variable_Coefficient_3D_k_out2_id3 = 0;
      Blending_F_Tetrahedron_blend( 0.44773255210137247*p_affine_0_0 + 0.18002969351036546*p_affine_1_0 + 0.36531451881463461*p_affine_2_0 + 0.0069232355736274509*p_affine_3_0, 0.44773255210137247*p_affine_0_1 + 0.18002969351036546*p_affine_1_1 + 0.36531451881463461*p_affine_2_1 + 0.0069232355736274509*p_affine_3_1, 0.44773255210137247*p_affine_0_2 + 0.18002969351036546*p_affine_1_2 + 0.36531451881463461*p_affine_2_2 + 0.0069232355736274509*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id0, &Blending_F_Tetrahedron_blend_out1_id0, &Blending_F_Tetrahedron_blend_out2_id0 );
      Blending_F_Tetrahedron_blend( 0.0048399363458718758*p_affine_0_0 + 0.15593312049918584*p_affine_1_0 + 0.45746158708559559*p_affine_2_0 + 0.3817653560693467*p_affine_3_0, 0.0048399363458718758*p_affine_0_1 + 0.15593312049918584*p_affine_1_1 + 0.45746158708559559*p_affine_2_1 + 0.3817653560693467*p_affine_3_1, 0.0048399363458718758*p_affine_0_2 + 0.15593312049918584*p_affine_1_2 + 0.45746158708559559*p_affine_2_2 + 0.3817653560693467*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Blending_F_Tetrahedron_blend( 0.35284634870858683*p_affine_0_0 + 0.21607642918484793*p_affine_1_0 + 0.00037551502872928966*p_affine_2_0 + 0.43070170707783589*p_affine_3_0, 0.35284634870858683*p_affine_0_1 + 0.21607642918484793*p_affine_1_1 + 0.00037551502872928966*p_affine_2_1 + 0.43070170707783589*p_affine_3_1, 0.35284634870858683*p_affine_0_2 + 0.21607642918484793*p_affine_1_2 + 0.00037551502872928966*p_affine_2_2 + 0.43070170707783589*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id2, &Blending_F_Tetrahedron_blend_out1_id2, &Blending_F_Tetrahedron_blend_out2_id2 );
      Blending_F_Tetrahedron_blend( 0.014827610062423724*p_affine_0_0 + 0.82157254096761967*p_affine_1_0 + 0.12366680032845823*p_affine_2_0 + 0.039933048641498381*p_affine_3_0, 0.014827610062423724*p_affine_0_1 + 0.82157254096761967*p_affine_1_1 + 0.12366680032845823*p_affine_2_1 + 0.039933048641498381*p_affine_3_1, 0.014827610062423724*p_affine_0_2 + 0.82157254096761967*p_affine_1_2 + 0.12366680032845823*p_affine_2_2 + 0.039933048641498381*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id3, &Blending_F_Tetrahedron_blend_out1_id3, &Blending_F_Tetrahedron_blend_out2_id3 );
      Vector_Variable_Coefficient_3D_k( Blending_F_Tetrahedron_blend_out0_id0, Blending_F_Tetrahedron_blend_out1_id0, Blending_F_Tetrahedron_blend_out2_id0, &Vector_Variable_Coefficient_3D_k_out0_id0, &Vector_Variable_Coefficient_3D_k_out1_id0, &Vector_Variable_Coefficient_3D_k_out2_id0 );
      Vector_Variable_Coefficient_3D_k( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Vector_Variable_Coefficient_3D_k_out0_id1, &Vector_Variable_Coefficient_3D_k_out1_id1, &Vector_Variable_Coefficient_3D_k_out2_id1 );
      Vector_Variable_Coefficient_3D_k( Blending_F_Tetrahedron_blend_out0_id2, Blending_F_Tetrahedron_blend_out1_id2, Blending_F_Tetrahedron_blend_out2_id2, &Vector_Variable_Coefficient_3D_k_out0_id2, &Vector_Variable_Coefficient_3D_k_out1_id2, &Vector_Variable_Coefficient_3D_k_out2_id2 );
      Vector_Variable_Coefficient_3D_k( Blending_F_Tetrahedron_blend_out0_id3, Blending_F_Tetrahedron_blend_out1_id3, Blending_F_Tetrahedron_blend_out2_id3, &Vector_Variable_Coefficient_3D_k_out0_id3, &Vector_Variable_Coefficient_3D_k_out1_id3, &Vector_Variable_Coefficient_3D_k_out2_id3 );
      real_t jac_blending_id0_0_0 = 0;
      real_t jac_blending_id0_0_1 = 0;
      real_t jac_blending_id0_0_2 = 0;
      real_t jac_blending_id0_1_0 = 0;
      real_t jac_blending_id0_1_1 = 0;
      real_t jac_blending_id0_1_2 = 0;
      real_t jac_blending_id0_2_0 = 0;
      real_t jac_blending_id0_2_1 = 0;
      real_t jac_blending_id0_2_2 = 0;
      Blending_DF_Tetrahedron_jac_blending( 0.44773255210137247*p_affine_0_0 + 0.18002969351036546*p_affine_1_0 + 0.36531451881463461*p_affine_2_0 + 0.0069232355736274509*p_affine_3_0, 0.44773255210137247*p_affine_0_1 + 0.18002969351036546*p_affine_1_1 + 0.36531451881463461*p_affine_2_1 + 0.0069232355736274509*p_affine_3_1, 0.44773255210137247*p_affine_0_2 + 0.18002969351036546*p_affine_1_2 + 0.36531451881463461*p_affine_2_2 + 0.0069232355736274509*p_affine_3_2, &jac_blending_id0_0_0, &jac_blending_id0_0_1, &jac_blending_id0_0_2, &jac_blending_id0_1_0, &jac_blending_id0_1_1, &jac_blending_id0_1_2, &jac_blending_id0_2_0, &jac_blending_id0_2_1, &jac_blending_id0_2_2 );
      real_t jac_blending_id1_0_0 = 0;
      real_t jac_blending_id1_0_1 = 0;
      real_t jac_blending_id1_0_2 = 0;
      real_t jac_blending_id1_1_0 = 0;
      real_t jac_blending_id1_1_1 = 0;
      real_t jac_blending_id1_1_2 = 0;
      real_t jac_blending_id1_2_0 = 0;
      real_t jac_blending_id1_2_1 = 0;
      real_t jac_blending_id1_2_2 = 0;
      Blending_DF_Tetrahedron_jac_blending( 0.0048399363458718758*p_affine_0_0 + 0.15593312049918584*p_affine_1_0 + 0.45746158708559559*p_affine_2_0 + 0.3817653560693467*p_affine_3_0, 0.0048399363458718758*p_affine_0_1 + 0.15593312049918584*p_affine_1_1 + 0.45746158708559559*p_affine_2_1 + 0.3817653560693467*p_affine_3_1, 0.0048399363458718758*p_affine_0_2 + 0.15593312049918584*p_affine_1_2 + 0.45746158708559559*p_affine_2_2 + 0.3817653560693467*p_affine_3_2, &jac_blending_id1_0_0, &jac_blending_id1_0_1, &jac_blending_id1_0_2, &jac_blending_id1_1_0, &jac_blending_id1_1_1, &jac_blending_id1_1_2, &jac_blending_id1_2_0, &jac_blending_id1_2_1, &jac_blending_id1_2_2 );
      real_t jac_blending_id2_0_0 = 0;
      real_t jac_blending_id2_0_1 = 0;
      real_t jac_blending_id2_0_2 = 0;
      real_t jac_blending_id2_1_0 = 0;
      real_t jac_blending_id2_1_1 = 0;
      real_t jac_blending_id2_1_2 = 0;
      real_t jac_blending_id2_2_0 = 0;
      real_t jac_blending_id2_2_1 = 0;
      real_t jac_blending_id2_2_2 = 0;
      Blending_DF_Tetrahedron_jac_blending( 0.35284634870858683*p_affine_0_0 + 0.21607642918484793*p_affine_1_0 + 0.00037551502872928966*p_affine_2_0 + 0.43070170707783589*p_affine_3_0, 0.35284634870858683*p_affine_0_1 + 0.21607642918484793*p_affine_1_1 + 0.00037551502872928966*p_affine_2_1 + 0.43070170707783589*p_affine_3_1, 0.35284634870858683*p_affine_0_2 + 0.21607642918484793*p_affine_1_2 + 0.00037551502872928966*p_affine_2_2 + 0.43070170707783589*p_affine_3_2, &jac_blending_id2_0_0, &jac_blending_id2_0_1, &jac_blending_id2_0_2, &jac_blending_id2_1_0, &jac_blending_id2_1_1, &jac_blending_id2_1_2, &jac_blending_id2_2_0, &jac_blending_id2_2_1, &jac_blending_id2_2_2 );
      real_t jac_blending_id3_0_0 = 0;
      real_t jac_blending_id3_0_1 = 0;
      real_t jac_blending_id3_0_2 = 0;
      real_t jac_blending_id3_1_0 = 0;
      real_t jac_blending_id3_1_1 = 0;
      real_t jac_blending_id3_1_2 = 0;
      real_t jac_blending_id3_2_0 = 0;
      real_t jac_blending_id3_2_1 = 0;
      real_t jac_blending_id3_2_2 = 0;
      Blending_DF_Tetrahedron_jac_blending( 0.014827610062423724*p_affine_0_0 + 0.82157254096761967*p_affine_1_0 + 0.12366680032845823*p_affine_2_0 + 0.039933048641498381*p_affine_3_0, 0.014827610062423724*p_affine_0_1 + 0.82157254096761967*p_affine_1_1 + 0.12366680032845823*p_affine_2_1 + 0.039933048641498381*p_affine_3_1, 0.014827610062423724*p_affine_0_2 + 0.82157254096761967*p_affine_1_2 + 0.12366680032845823*p_affine_2_2 + 0.039933048641498381*p_affine_3_2, &jac_blending_id3_0_0, &jac_blending_id3_0_1, &jac_blending_id3_0_2, &jac_blending_id3_1_0, &jac_blending_id3_1_1, &jac_blending_id3_1_2, &jac_blending_id3_2_0, &jac_blending_id3_2_1, &jac_blending_id3_2_2 );
      real_t jac_affine_0_0 = -1.0*p_affine_0_0 + 1.0*p_affine_1_0;
      real_t jac_affine_0_1 = -1.0*p_affine_0_0 + 1.0*p_affine_2_0;
      real_t jac_affine_0_2 = -1.0*p_affine_0_0 + 1.0*p_affine_3_0;
      real_t jac_affine_1_0 = -1.0*p_affine_0_1 + 1.0*p_affine_1_1;
      real_t jac_affine_1_1 = -1.0*p_affine_0_1 + 1.0*p_affine_2_1;
      real_t jac_affine_1_2 = -1.0*p_affine_0_1 + 1.0*p_affine_3_1;
      real_t jac_affine_2_0 = -1.0*p_affine_0_2 + 1.0*p_affine_1_2;
      real_t jac_affine_2_1 = -1.0*p_affine_0_2 + 1.0*p_affine_2_2;
      real_t jac_affine_2_2 = -1.0*p_affine_0_2 + 1.0*p_affine_3_2;
      real_t jac_id0_0_0 = 1.0*jac_affine_0_0*jac_blending_id0_0_0 + 1.0*jac_affine_1_0*jac_blending_id0_0_1 + 1.0*jac_affine_2_0*jac_blending_id0_0_2;
      real_t jac_id0_0_1 = 1.0*jac_affine_0_1*jac_blending_id0_0_0 + 1.0*jac_affine_1_1*jac_blending_id0_0_1 + 1.0*jac_affine_2_1*jac_blending_id0_0_2;
      real_t jac_id0_0_2 = 1.0*jac_affine_0_2*jac_blending_id0_0_0 + 1.0*jac_affine_1_2*jac_blending_id0_0_1 + 1.0*jac_affine_2_2*jac_blending_id0_0_2;
      real_t jac_id0_1_0 = 1.0*jac_affine_0_0*jac_blending_id0_1_0 + 1.0*jac_affine_1_0*jac_blending_id0_1_1 + 1.0*jac_affine_2_0*jac_blending_id0_1_2;
      real_t jac_id0_1_1 = 1.0*jac_affine_0_1*jac_blending_id0_1_0 + 1.0*jac_affine_1_1*jac_blending_id0_1_1 + 1.0*jac_affine_2_1*jac_blending_id0_1_2;
      real_t jac_id0_1_2 = 1.0*jac_affine_0_2*jac_blending_id0_1_0 + 1.0*jac_affine_1_2*jac_blending_id0_1_1 + 1.0*jac_affine_2_2*jac_blending_id0_1_2;
      real_t jac_id0_2_0 = 1.0*jac_affine_0_0*jac_blending_id0_2_0 + 1.0*jac_affine_1_0*jac_blending_id0_2_1 + 1.0*jac_affine_2_0*jac_blending_id0_2_2;
      real_t jac_id0_2_1 = 1.0*jac_affine_0_1*jac_blending_id0_2_0 + 1.0*jac_affine_1_1*jac_blending_id0_2_1 + 1.0*jac_affine_2_1*jac_blending_id0_2_2;
      real_t jac_id0_2_2 = 1.0*jac_affine_0_2*jac_blending_id0_2_0 + 1.0*jac_affine_1_2*jac_blending_id0_2_1 + 1.0*jac_affine_2_2*jac_blending_id0_2_2;
      real_t tmp_34 = 1.0*jac_id0_0_0*jac_id0_1_1 - 1.0*jac_id0_0_1*jac_id0_1_0;
      real_t tmp_36 = 1.0*jac_id0_1_0*jac_id0_2_1 - 1.0*jac_id0_1_1*jac_id0_2_0;
      real_t tmp_35 = -1.0*jac_id0_0_0*jac_id0_2_1 + 1.0*jac_id0_0_1*jac_id0_2_0;
      real_t tmp_23 = tmp_34*jac_id0_2_2 + tmp_35*jac_id0_1_2 + tmp_36*jac_id0_0_2;
      real_t tmp_np_0_arg_0 = tmp_23;
      real_t tmp_np_0 = std::abs(tmp_np_0_arg_0);
      real_t abs_det_jac_id0 = 1.0*tmp_np_0;
      real_t tmp_np_1_arg_0 = tmp_23;
      real_t tmp_np_1 = 1.0 / (tmp_np_1_arg_0);
      real_t tmp_37 = 1.0*jac_id0_1_1*jac_id0_2_2 - 1.0*jac_id0_1_2*jac_id0_2_1;
      real_t jac_inv_id0_0_0 = tmp_37*tmp_np_1;
      real_t tmp_38 = -1.0*jac_id0_0_1*jac_id0_2_2 + 1.0*jac_id0_0_2*jac_id0_2_1;
      real_t jac_inv_id0_0_1 = tmp_38*tmp_np_1;
      real_t tmp_39 = 1.0*jac_id0_0_1*jac_id0_1_2 - 1.0*jac_id0_0_2*jac_id0_1_1;
      real_t jac_inv_id0_0_2 = tmp_39*tmp_np_1;
      real_t tmp_40 = -1.0*jac_id0_1_0*jac_id0_2_2 + 1.0*jac_id0_1_2*jac_id0_2_0;
      real_t jac_inv_id0_1_0 = tmp_40*tmp_np_1;
      real_t tmp_41 = 1.0*jac_id0_0_0*jac_id0_2_2 - 1.0*jac_id0_0_2*jac_id0_2_0;
      real_t jac_inv_id0_1_1 = tmp_41*tmp_np_1;
      real_t tmp_42 = -1.0*jac_id0_0_0*jac_id0_1_2 + 1.0*jac_id0_0_2*jac_id0_1_0;
      real_t jac_inv_id0_1_2 = tmp_42*tmp_np_1;
      real_t jac_inv_id0_2_0 = tmp_36*tmp_np_1;
      real_t jac_inv_id0_2_1 = tmp_35*tmp_np_1;
      real_t jac_inv_id0_2_2 = tmp_34*tmp_np_1;
      real_t jac_id1_0_0 = 1.0*jac_affine_0_0*jac_blending_id1_0_0 + 1.0*jac_affine_1_0*jac_blending_id1_0_1 + 1.0*jac_affine_2_0*jac_blending_id1_0_2;
      real_t jac_id1_0_1 = 1.0*jac_affine_0_1*jac_blending_id1_0_0 + 1.0*jac_affine_1_1*jac_blending_id1_0_1 + 1.0*jac_affine_2_1*jac_blending_id1_0_2;
      real_t jac_id1_0_2 = 1.0*jac_affine_0_2*jac_blending_id1_0_0 + 1.0*jac_affine_1_2*jac_blending_id1_0_1 + 1.0*jac_affine_2_2*jac_blending_id1_0_2;
      real_t jac_id1_1_0 = 1.0*jac_affine_0_0*jac_blending_id1_1_0 + 1.0*jac_affine_1_0*jac_blending_id1_1_1 + 1.0*jac_affine_2_0*jac_blending_id1_1_2;
      real_t jac_id1_1_1 = 1.0*jac_affine_0_1*jac_blending_id1_1_0 + 1.0*jac_affine_1_1*jac_blending_id1_1_1 + 1.0*jac_affine_2_1*jac_blending_id1_1_2;
      real_t jac_id1_1_2 = 1.0*jac_affine_0_2*jac_blending_id1_1_0 + 1.0*jac_affine_1_2*jac_blending_id1_1_1 + 1.0*jac_affine_2_2*jac_blending_id1_1_2;
      real_t jac_id1_2_0 = 1.0*jac_affine_0_0*jac_blending_id1_2_0 + 1.0*jac_affine_1_0*jac_blending_id1_2_1 + 1.0*jac_affine_2_0*jac_blending_id1_2_2;
      real_t jac_id1_2_1 = 1.0*jac_affine_0_1*jac_blending_id1_2_0 + 1.0*jac_affine_1_1*jac_blending_id1_2_1 + 1.0*jac_affine_2_1*jac_blending_id1_2_2;
      real_t jac_id1_2_2 = 1.0*jac_affine_0_2*jac_blending_id1_2_0 + 1.0*jac_affine_1_2*jac_blending_id1_2_1 + 1.0*jac_affine_2_2*jac_blending_id1_2_2;
      real_t tmp_33 = 1.0*jac_id1_1_0*jac_id1_2_1 - 1.0*jac_id1_1_1*jac_id1_2_0;
      real_t tmp_32 = -1.0*jac_id1_0_0*jac_id1_2_1 + 1.0*jac_id1_0_1*jac_id1_2_0;
      real_t tmp_31 = 1.0*jac_id1_0_0*jac_id1_1_1 - 1.0*jac_id1_0_1*jac_id1_1_0;
      real_t tmp_22 = tmp_31*jac_id1_2_2 + tmp_32*jac_id1_1_2 + tmp_33*jac_id1_0_2;
      real_t tmp_np_2_arg_0 = tmp_22;
      real_t tmp_np_2 = std::abs(tmp_np_2_arg_0);
      real_t abs_det_jac_id1 = 1.0*tmp_np_2;
      real_t tmp_np_3_arg_0 = tmp_22;
      real_t tmp_np_3 = 1.0 / (tmp_np_3_arg_0);
      real_t tmp_43 = 1.0*jac_id1_1_1*jac_id1_2_2 - 1.0*jac_id1_1_2*jac_id1_2_1;
      real_t jac_inv_id1_0_0 = tmp_43*tmp_np_3;
      real_t tmp_44 = -1.0*jac_id1_0_1*jac_id1_2_2 + 1.0*jac_id1_0_2*jac_id1_2_1;
      real_t jac_inv_id1_0_1 = tmp_44*tmp_np_3;
      real_t tmp_45 = 1.0*jac_id1_0_1*jac_id1_1_2 - 1.0*jac_id1_0_2*jac_id1_1_1;
      real_t jac_inv_id1_0_2 = tmp_45*tmp_np_3;
      real_t tmp_46 = -1.0*jac_id1_1_0*jac_id1_2_2 + 1.0*jac_id1_1_2*jac_id1_2_0;
      real_t jac_inv_id1_1_0 = tmp_46*tmp_np_3;
      real_t tmp_47 = 1.0*jac_id1_0_0*jac_id1_2_2 - 1.0*jac_id1_0_2*jac_id1_2_0;
      real_t jac_inv_id1_1_1 = tmp_47*tmp_np_3;
      real_t tmp_48 = -1.0*jac_id1_0_0*jac_id1_1_2 + 1.0*jac_id1_0_2*jac_id1_1_0;
      real_t jac_inv_id1_1_2 = tmp_48*tmp_np_3;
      real_t jac_inv_id1_2_0 = tmp_33*tmp_np_3;
      real_t jac_inv_id1_2_1 = tmp_32*tmp_np_3;
      real_t jac_inv_id1_2_2 = tmp_31*tmp_np_3;
      real_t jac_id2_0_0 = 1.0*jac_affine_0_0*jac_blending_id2_0_0 + 1.0*jac_affine_1_0*jac_blending_id2_0_1 + 1.0*jac_affine_2_0*jac_blending_id2_0_2;
      real_t jac_id2_0_1 = 1.0*jac_affine_0_1*jac_blending_id2_0_0 + 1.0*jac_affine_1_1*jac_blending_id2_0_1 + 1.0*jac_affine_2_1*jac_blending_id2_0_2;
      real_t jac_id2_0_2 = 1.0*jac_affine_0_2*jac_blending_id2_0_0 + 1.0*jac_affine_1_2*jac_blending_id2_0_1 + 1.0*jac_affine_2_2*jac_blending_id2_0_2;
      real_t jac_id2_1_0 = 1.0*jac_affine_0_0*jac_blending_id2_1_0 + 1.0*jac_affine_1_0*jac_blending_id2_1_1 + 1.0*jac_affine_2_0*jac_blending_id2_1_2;
      real_t jac_id2_1_1 = 1.0*jac_affine_0_1*jac_blending_id2_1_0 + 1.0*jac_affine_1_1*jac_blending_id2_1_1 + 1.0*jac_affine_2_1*jac_blending_id2_1_2;
      real_t jac_id2_1_2 = 1.0*jac_affine_0_2*jac_blending_id2_1_0 + 1.0*jac_affine_1_2*jac_blending_id2_1_1 + 1.0*jac_affine_2_2*jac_blending_id2_1_2;
      real_t jac_id2_2_0 = 1.0*jac_affine_0_0*jac_blending_id2_2_0 + 1.0*jac_affine_1_0*jac_blending_id2_2_1 + 1.0*jac_affine_2_0*jac_blending_id2_2_2;
      real_t jac_id2_2_1 = 1.0*jac_affine_0_1*jac_blending_id2_2_0 + 1.0*jac_affine_1_1*jac_blending_id2_2_1 + 1.0*jac_affine_2_1*jac_blending_id2_2_2;
      real_t jac_id2_2_2 = 1.0*jac_affine_0_2*jac_blending_id2_2_0 + 1.0*jac_affine_1_2*jac_blending_id2_2_1 + 1.0*jac_affine_2_2*jac_blending_id2_2_2;
      real_t tmp_28 = 1.0*jac_id2_0_0*jac_id2_1_1 - 1.0*jac_id2_0_1*jac_id2_1_0;
      real_t tmp_30 = 1.0*jac_id2_1_0*jac_id2_2_1 - 1.0*jac_id2_1_1*jac_id2_2_0;
      real_t tmp_29 = -1.0*jac_id2_0_0*jac_id2_2_1 + 1.0*jac_id2_0_1*jac_id2_2_0;
      real_t tmp_24 = tmp_28*jac_id2_2_2 + tmp_29*jac_id2_1_2 + tmp_30*jac_id2_0_2;
      real_t tmp_np_4_arg_0 = tmp_24;
      real_t tmp_np_4 = std::abs(tmp_np_4_arg_0);
      real_t abs_det_jac_id2 = 1.0*tmp_np_4;
      real_t tmp_np_5_arg_0 = tmp_24;
      real_t tmp_np_5 = 1.0 / (tmp_np_5_arg_0);
      real_t tmp_49 = 1.0*jac_id2_1_1*jac_id2_2_2 - 1.0*jac_id2_1_2*jac_id2_2_1;
      real_t jac_inv_id2_0_0 = tmp_49*tmp_np_5;
      real_t tmp_50 = -1.0*jac_id2_0_1*jac_id2_2_2 + 1.0*jac_id2_0_2*jac_id2_2_1;
      real_t jac_inv_id2_0_1 = tmp_50*tmp_np_5;
      real_t tmp_51 = 1.0*jac_id2_0_1*jac_id2_1_2 - 1.0*jac_id2_0_2*jac_id2_1_1;
      real_t jac_inv_id2_0_2 = tmp_51*tmp_np_5;
      real_t tmp_52 = -1.0*jac_id2_1_0*jac_id2_2_2 + 1.0*jac_id2_1_2*jac_id2_2_0;
      real_t jac_inv_id2_1_0 = tmp_52*tmp_np_5;
      real_t tmp_53 = 1.0*jac_id2_0_0*jac_id2_2_2 - 1.0*jac_id2_0_2*jac_id2_2_0;
      real_t jac_inv_id2_1_1 = tmp_53*tmp_np_5;
      real_t tmp_54 = -1.0*jac_id2_0_0*jac_id2_1_2 + 1.0*jac_id2_0_2*jac_id2_1_0;
      real_t jac_inv_id2_1_2 = tmp_54*tmp_np_5;
      real_t jac_inv_id2_2_0 = tmp_30*tmp_np_5;
      real_t jac_inv_id2_2_1 = tmp_29*tmp_np_5;
      real_t jac_inv_id2_2_2 = tmp_28*tmp_np_5;
      real_t jac_id3_0_0 = 1.0*jac_affine_0_0*jac_blending_id3_0_0 + 1.0*jac_affine_1_0*jac_blending_id3_0_1 + 1.0*jac_affine_2_0*jac_blending_id3_0_2;
      real_t jac_id3_0_1 = 1.0*jac_affine_0_1*jac_blending_id3_0_0 + 1.0*jac_affine_1_1*jac_blending_id3_0_1 + 1.0*jac_affine_2_1*jac_blending_id3_0_2;
      real_t jac_id3_0_2 = 1.0*jac_affine_0_2*jac_blending_id3_0_0 + 1.0*jac_affine_1_2*jac_blending_id3_0_1 + 1.0*jac_affine_2_2*jac_blending_id3_0_2;
      real_t jac_id3_1_0 = 1.0*jac_affine_0_0*jac_blending_id3_1_0 + 1.0*jac_affine_1_0*jac_blending_id3_1_1 + 1.0*jac_affine_2_0*jac_blending_id3_1_2;
      real_t jac_id3_1_1 = 1.0*jac_affine_0_1*jac_blending_id3_1_0 + 1.0*jac_affine_1_1*jac_blending_id3_1_1 + 1.0*jac_affine_2_1*jac_blending_id3_1_2;
      real_t jac_id3_1_2 = 1.0*jac_affine_0_2*jac_blending_id3_1_0 + 1.0*jac_affine_1_2*jac_blending_id3_1_1 + 1.0*jac_affine_2_2*jac_blending_id3_1_2;
      real_t jac_id3_2_0 = 1.0*jac_affine_0_0*jac_blending_id3_2_0 + 1.0*jac_affine_1_0*jac_blending_id3_2_1 + 1.0*jac_affine_2_0*jac_blending_id3_2_2;
      real_t jac_id3_2_1 = 1.0*jac_affine_0_1*jac_blending_id3_2_0 + 1.0*jac_affine_1_1*jac_blending_id3_2_1 + 1.0*jac_affine_2_1*jac_blending_id3_2_2;
      real_t jac_id3_2_2 = 1.0*jac_affine_0_2*jac_blending_id3_2_0 + 1.0*jac_affine_1_2*jac_blending_id3_2_1 + 1.0*jac_affine_2_2*jac_blending_id3_2_2;
      real_t tmp_25 = 1.0*jac_id3_0_0*jac_id3_1_1 - 1.0*jac_id3_0_1*jac_id3_1_0;
      real_t tmp_27 = 1.0*jac_id3_1_0*jac_id3_2_1 - 1.0*jac_id3_1_1*jac_id3_2_0;
      real_t tmp_26 = -1.0*jac_id3_0_0*jac_id3_2_1 + 1.0*jac_id3_0_1*jac_id3_2_0;
      real_t tmp_21 = tmp_25*jac_id3_2_2 + tmp_26*jac_id3_1_2 + tmp_27*jac_id3_0_2;
      real_t tmp_np_6_arg_0 = tmp_21;
      real_t tmp_np_6 = std::abs(tmp_np_6_arg_0);
      real_t abs_det_jac_id3 = 1.0*tmp_np_6;
      real_t tmp_np_7_arg_0 = tmp_21;
      real_t tmp_np_7 = 1.0 / (tmp_np_7_arg_0);
      real_t tmp_55 = 1.0*jac_id3_1_1*jac_id3_2_2 - 1.0*jac_id3_1_2*jac_id3_2_1;
      real_t jac_inv_id3_0_0 = tmp_55*tmp_np_7;
      real_t tmp_56 = -1.0*jac_id3_0_1*jac_id3_2_2 + 1.0*jac_id3_0_2*jac_id3_2_1;
      real_t jac_inv_id3_0_1 = tmp_56*tmp_np_7;
      real_t tmp_57 = 1.0*jac_id3_0_1*jac_id3_1_2 - 1.0*jac_id3_0_2*jac_id3_1_1;
      real_t jac_inv_id3_0_2 = tmp_57*tmp_np_7;
      real_t tmp_58 = -1.0*jac_id3_1_0*jac_id3_2_2 + 1.0*jac_id3_1_2*jac_id3_2_0;
      real_t jac_inv_id3_1_0 = tmp_58*tmp_np_7;
      real_t tmp_59 = 1.0*jac_id3_0_0*jac_id3_2_2 - 1.0*jac_id3_0_2*jac_id3_2_0;
      real_t jac_inv_id3_1_1 = tmp_59*tmp_np_7;
      real_t tmp_60 = -1.0*jac_id3_0_0*jac_id3_1_2 + 1.0*jac_id3_0_2*jac_id3_1_0;
      real_t jac_inv_id3_1_2 = tmp_60*tmp_np_7;
      real_t jac_inv_id3_2_0 = tmp_27*tmp_np_7;
      real_t jac_inv_id3_2_1 = tmp_26*tmp_np_7;
      real_t jac_inv_id3_2_2 = tmp_25*tmp_np_7;
      real_t tmp_17 = Vector_Variable_Coefficient_3D_k_out0_id0*jac_inv_id0_1_0 + Vector_Variable_Coefficient_3D_k_out1_id0*jac_inv_id0_1_1 + Vector_Variable_Coefficient_3D_k_out2_id0*jac_inv_id0_1_2;
      real_t tmp_20 = Vector_Variable_Coefficient_3D_k_out0_id1*jac_inv_id1_1_0 + Vector_Variable_Coefficient_3D_k_out1_id1*jac_inv_id1_1_1 + Vector_Variable_Coefficient_3D_k_out2_id1*jac_inv_id1_1_2;
      real_t tmp_18 = Vector_Variable_Coefficient_3D_k_out0_id3*jac_inv_id3_1_0 + Vector_Variable_Coefficient_3D_k_out1_id3*jac_inv_id3_1_1 + Vector_Variable_Coefficient_3D_k_out2_id3*jac_inv_id3_1_2;
      real_t tmp_19 = Vector_Variable_Coefficient_3D_k_out0_id2*jac_inv_id2_1_0 + Vector_Variable_Coefficient_3D_k_out1_id2*jac_inv_id2_1_1 + Vector_Variable_Coefficient_3D_k_out2_id2*jac_inv_id2_1_2;
      real_t tmp_5 = 0.00034676287630628196*abs_det_jac_id0*tmp_17 + 0.01773793680464952*abs_det_jac_id1*tmp_20 + 0.022905717123111669*abs_det_jac_id2*tmp_19 + 0.00067624986259913489*abs_det_jac_id3*tmp_18;
      real_t tmp_16 = Vector_Variable_Coefficient_3D_k_out0_id0*jac_inv_id0_2_0 + Vector_Variable_Coefficient_3D_k_out1_id0*jac_inv_id0_2_1 + Vector_Variable_Coefficient_3D_k_out2_id0*jac_inv_id0_2_2;
      real_t tmp_14 = Vector_Variable_Coefficient_3D_k_out0_id2*jac_inv_id2_2_0 + Vector_Variable_Coefficient_3D_k_out1_id2*jac_inv_id2_2_1 + Vector_Variable_Coefficient_3D_k_out2_id2*jac_inv_id2_2_2;
      real_t tmp_13 = Vector_Variable_Coefficient_3D_k_out0_id3*jac_inv_id3_2_0 + Vector_Variable_Coefficient_3D_k_out1_id3*jac_inv_id3_2_1 + Vector_Variable_Coefficient_3D_k_out2_id3*jac_inv_id3_2_2;
      real_t tmp_15 = Vector_Variable_Coefficient_3D_k_out0_id1*jac_inv_id1_2_0 + Vector_Variable_Coefficient_3D_k_out1_id1*jac_inv_id1_2_1 + Vector_Variable_Coefficient_3D_k_out2_id1*jac_inv_id1_2_2;
      real_t tmp_6 = 0.018297443724601584*abs_det_jac_id0*tmp_16 + 0.021255005445818931*abs_det_jac_id1*tmp_15 + 1.9970761392863065e-5*abs_det_jac_id2*tmp_14 + 0.0020942467348532634*abs_det_jac_id3*tmp_13;
      real_t tmp_2 = -1.0*tmp_5 + 1.0*tmp_6;
      real_t a_0_0 = tmp_2;
      real_t a_0_1 = 0;
      real_t a_0_2 = 0;
      real_t a_0_3 = 0;
      real_t a_0_4 = 0;
      real_t a_0_5 = 0;
      real_t a_1_0 = 0;
      real_t tmp_10 = Vector_Variable_Coefficient_3D_k_out0_id2*jac_inv_id2_0_0 + Vector_Variable_Coefficient_3D_k_out1_id2*jac_inv_id2_0_1 + Vector_Variable_Coefficient_3D_k_out2_id2*jac_inv_id2_0_2;
      real_t tmp_12 = Vector_Variable_Coefficient_3D_k_out0_id0*jac_inv_id0_0_0 + Vector_Variable_Coefficient_3D_k_out1_id0*jac_inv_id0_0_1 + Vector_Variable_Coefficient_3D_k_out2_id0*jac_inv_id0_0_2;
      real_t tmp_9 = Vector_Variable_Coefficient_3D_k_out0_id3*jac_inv_id3_0_0 + Vector_Variable_Coefficient_3D_k_out1_id3*jac_inv_id3_0_1 + Vector_Variable_Coefficient_3D_k_out2_id3*jac_inv_id3_0_2;
      real_t tmp_11 = Vector_Variable_Coefficient_3D_k_out0_id1*jac_inv_id1_0_0 + Vector_Variable_Coefficient_3D_k_out1_id1*jac_inv_id1_0_1 + Vector_Variable_Coefficient_3D_k_out2_id1*jac_inv_id1_0_2;
      real_t tmp_3 = 0.00034676287630628196*abs_det_jac_id0*tmp_12 + 0.01773793680464952*abs_det_jac_id1*tmp_11 + 0.022905717123111669*abs_det_jac_id2*tmp_10 + 0.00067624986259913489*abs_det_jac_id3*tmp_9;
      real_t tmp_4 = 0.0090171154337138296*abs_det_jac_id0*tmp_16 + 0.0072451095763229299*abs_det_jac_id1*tmp_15 + 0.011491446359616473*abs_det_jac_id2*tmp_14 + 0.013912995297013415*abs_det_jac_id3*tmp_13;
      real_t tmp_1 = -1.0*tmp_3 + 1.0*tmp_4;
      real_t a_1_1 = tmp_1;
      real_t a_1_2 = 0;
      real_t a_1_3 = 0;
      real_t a_1_4 = 0;
      real_t a_1_5 = 0;
      real_t a_2_0 = 0;
      real_t a_2_1 = 0;
      real_t tmp_8 = 0.018297443724601584*abs_det_jac_id0*tmp_12 + 0.021255005445818931*abs_det_jac_id1*tmp_11 + 1.9970761392863065e-5*abs_det_jac_id2*tmp_10 + 0.0020942467348532634*abs_det_jac_id3*tmp_9;
      real_t tmp_7 = 0.0090171154337138296*abs_det_jac_id0*tmp_17 + 0.0072451095763229299*abs_det_jac_id1*tmp_20 + 0.011491446359616473*abs_det_jac_id2*tmp_19 + 0.013912995297013415*abs_det_jac_id3*tmp_18;
      real_t tmp_0 = 1.0*tmp_7 - 1.0*tmp_8;
      real_t a_2_2 = tmp_0;
      real_t a_2_3 = 0;
      real_t a_2_4 = 0;
      real_t a_2_5 = 0;
      real_t a_3_0 = 0;
      real_t a_3_1 = 0;
      real_t a_3_2 = 0;
      real_t a_3_3 = 0.050086823222829389*abs_det_jac_id0*tmp_16 + 0.046462929447761279*abs_det_jac_id1*tmp_15 + 0.05318232258357912*abs_det_jac_id2*tmp_14 + 0.016934591412496786*abs_det_jac_id3*tmp_13 + 1.0*tmp_3 - 1.0*tmp_4 + 1.0*tmp_5 - 1.0*tmp_6;
      real_t a_3_4 = 0;
      real_t a_3_5 = 0;
      real_t a_4_0 = 0;
      real_t a_4_1 = 0;
      real_t a_4_2 = 0;
      real_t a_4_3 = 0;
      real_t a_4_4 = 0.050086823222829389*abs_det_jac_id0*tmp_17 + 0.046462929447761279*abs_det_jac_id1*tmp_20 + 0.05318232258357912*abs_det_jac_id2*tmp_19 + 0.016934591412496786*abs_det_jac_id3*tmp_18 + tmp_2 - 1.0*tmp_7 + 1.0*tmp_8;
      real_t a_4_5 = 0;
      real_t a_5_0 = 0;
      real_t a_5_1 = 0;
      real_t a_5_2 = 0;
      real_t a_5_3 = 0;
      real_t a_5_4 = 0;
      real_t a_5_5 = 0.050086823222829389*abs_det_jac_id0*tmp_12 + 0.046462929447761279*abs_det_jac_id1*tmp_11 + 0.05318232258357912*abs_det_jac_id2*tmp_10 + 0.016934591412496786*abs_det_jac_id3*tmp_9 + tmp_0 + tmp_1;
      elMat(0,0) = a_0_0;
      elMat(0,1) = a_0_1;
      elMat(0,2) = a_0_2;
      elMat(0,3) = a_0_3;
      elMat(0,4) = a_0_4;
      elMat(0,5) = a_0_5;
      elMat(1,0) = a_1_0;
      elMat(1,1) = a_1_1;
      elMat(1,2) = a_1_2;
      elMat(1,3) = a_1_3;
      elMat(1,4) = a_1_4;
      elMat(1,5) = a_1_5;
      elMat(2,0) = a_2_0;
      elMat(2,1) = a_2_1;
      elMat(2,2) = a_2_2;
      elMat(2,3) = a_2_3;
      elMat(2,4) = a_2_4;
      elMat(2,5) = a_2_5;
      elMat(3,0) = a_3_0;
      elMat(3,1) = a_3_1;
      elMat(3,2) = a_3_2;
      elMat(3,3) = a_3_3;
      elMat(3,4) = a_3_4;
      elMat(3,5) = a_3_5;
      elMat(4,0) = a_4_0;
      elMat(4,1) = a_4_1;
      elMat(4,2) = a_4_2;
      elMat(4,3) = a_4_3;
      elMat(4,4) = a_4_4;
      elMat(4,5) = a_4_5;
      elMat(5,0) = a_5_0;
      elMat(5,1) = a_5_1;
      elMat(5,2) = a_5_2;
      elMat(5,3) = a_5_3;
      elMat(5,4) = a_5_4;
      elMat(5,5) = a_5_5;
   }

   void n1e1_linear_form_blending_q2::Blending_F_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2 ) const
   {
      Point3D  in( {in_0, in_1, in_2} );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
      *out_2 = out[2];
   }

   void n1e1_linear_form_blending_q2::Vector_Variable_Coefficient_3D_k( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2 ) const
   {
      Point3D result_Vector_Variable_Coefficient_3D_k = callback_Vector_Variable_Coefficient_3D_k( Point3D{ in_0, in_1, in_2 } );
      *out_0 = result_Vector_Variable_Coefficient_3D_k[0];
      *out_1 = result_Vector_Variable_Coefficient_3D_k[1];
      *out_2 = result_Vector_Variable_Coefficient_3D_k[2];
   }

   void n1e1_linear_form_blending_q2::Blending_DF_Tetrahedron_jac_blending( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
   {
      Point3D  mappedPt( {in_0, in_1, in_2} );
      Matrix3r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 0, 2 );
      *out_3 = DPsi( 1, 0 );
      *out_4 = DPsi( 1, 1 );
      *out_5 = DPsi( 1, 2 );
      *out_6 = DPsi( 2, 0 );
      *out_7 = DPsi( 2, 1 );
      *out_8 = DPsi( 2, 2 );
   }

} // namespace forms
} // namespace hyteg

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

#include "p1_diffusion_blending_q2.hpp"
#include "core/Abort.h"

namespace hyteg {
namespace forms {

   void p1_diffusion_blending_q2::integrateAll( const std::array< Point3D, 3 >&, Matrix< real_t, 3, 3 >& ) const
   {
      WALBERLA_ABORT("Not implemented!")
   }

   void p1_diffusion_blending_q2::integrateRow0( const std::array< Point3D, 3 >&, Matrix< real_t, 1, 3 >& ) const
   {
      WALBERLA_ABORT("Not implemented!")
   }

   void p1_diffusion_blending_q2::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 4, 4 >& elMat ) const
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
      real_t jac_id0_1_2 = 1.0*jac_affine_0_2*jac_blending_id0_1_0 + 1.0*jac_affine_1_2*jac_blending_id0_1_1 + 1.0*jac_affine_2_2*jac_blending_id0_1_2;
      real_t jac_id0_2_1 = 1.0*jac_affine_0_1*jac_blending_id0_2_0 + 1.0*jac_affine_1_1*jac_blending_id0_2_1 + 1.0*jac_affine_2_1*jac_blending_id0_2_2;
      real_t jac_id0_1_1 = 1.0*jac_affine_0_1*jac_blending_id0_1_0 + 1.0*jac_affine_1_1*jac_blending_id0_1_1 + 1.0*jac_affine_2_1*jac_blending_id0_1_2;
      real_t jac_id0_1_0 = 1.0*jac_affine_0_0*jac_blending_id0_1_0 + 1.0*jac_affine_1_0*jac_blending_id0_1_1 + 1.0*jac_affine_2_0*jac_blending_id0_1_2;
      real_t jac_id0_2_0 = 1.0*jac_affine_0_0*jac_blending_id0_2_0 + 1.0*jac_affine_1_0*jac_blending_id0_2_1 + 1.0*jac_affine_2_0*jac_blending_id0_2_2;
      real_t tmp_21 = 1.0*jac_id0_1_0*jac_id0_2_1 - 1.0*jac_id0_1_1*jac_id0_2_0;
      real_t jac_id0_0_2 = 1.0*jac_affine_0_2*jac_blending_id0_0_0 + 1.0*jac_affine_1_2*jac_blending_id0_0_1 + 1.0*jac_affine_2_2*jac_blending_id0_0_2;
      real_t jac_id0_2_2 = 1.0*jac_affine_0_2*jac_blending_id0_2_0 + 1.0*jac_affine_1_2*jac_blending_id0_2_1 + 1.0*jac_affine_2_2*jac_blending_id0_2_2;
      real_t jac_id0_0_1 = 1.0*jac_affine_0_1*jac_blending_id0_0_0 + 1.0*jac_affine_1_1*jac_blending_id0_0_1 + 1.0*jac_affine_2_1*jac_blending_id0_0_2;
      real_t jac_id0_0_0 = 1.0*jac_affine_0_0*jac_blending_id0_0_0 + 1.0*jac_affine_1_0*jac_blending_id0_0_1 + 1.0*jac_affine_2_0*jac_blending_id0_0_2;
      real_t tmp_19 = 1.0*jac_id0_0_0*jac_id0_1_1 - 1.0*jac_id0_0_1*jac_id0_1_0;
      real_t tmp_20 = -1.0*jac_id0_0_0*jac_id0_2_1 + 1.0*jac_id0_0_1*jac_id0_2_0;
      real_t tmp_7 = jac_id0_0_2*tmp_21 + jac_id0_1_2*tmp_20 + jac_id0_2_2*tmp_19;
      real_t np_0_arg_0 = tmp_7;
      real_t np_0 = std::abs(np_0_arg_0);
      real_t abs_det_jac_id0 = 1.0*np_0;
      real_t tmp_49 = 1.0*jac_id0_1_1*jac_id0_2_2 - 1.0*jac_id0_1_2*jac_id0_2_1;
      real_t np_1_arg_0 = tmp_7;
      real_t np_1 = 1.0 / (np_1_arg_0);
      real_t jac_inv_id0_0_0 = np_1*tmp_49;
      real_t tmp_50 = -1.0*jac_id0_0_1*jac_id0_2_2 + 1.0*jac_id0_0_2*jac_id0_2_1;
      real_t jac_inv_id0_0_1 = np_1*tmp_50;
      real_t tmp_51 = 1.0*jac_id0_0_1*jac_id0_1_2 - 1.0*jac_id0_0_2*jac_id0_1_1;
      real_t jac_inv_id0_0_2 = np_1*tmp_51;
      real_t tmp_52 = -1.0*jac_id0_1_0*jac_id0_2_2 + 1.0*jac_id0_1_2*jac_id0_2_0;
      real_t jac_inv_id0_1_0 = np_1*tmp_52;
      real_t tmp_53 = 1.0*jac_id0_0_0*jac_id0_2_2 - 1.0*jac_id0_0_2*jac_id0_2_0;
      real_t jac_inv_id0_1_1 = np_1*tmp_53;
      real_t tmp_54 = -1.0*jac_id0_0_0*jac_id0_1_2 + 1.0*jac_id0_0_2*jac_id0_1_0;
      real_t jac_inv_id0_1_2 = np_1*tmp_54;
      real_t jac_inv_id0_2_0 = np_1*tmp_21;
      real_t jac_inv_id0_2_1 = np_1*tmp_20;
      real_t jac_inv_id0_2_2 = np_1*tmp_19;
      real_t jac_id1_2_0 = 1.0*jac_affine_0_0*jac_blending_id1_2_0 + 1.0*jac_affine_1_0*jac_blending_id1_2_1 + 1.0*jac_affine_2_0*jac_blending_id1_2_2;
      real_t jac_id1_0_0 = 1.0*jac_affine_0_0*jac_blending_id1_0_0 + 1.0*jac_affine_1_0*jac_blending_id1_0_1 + 1.0*jac_affine_2_0*jac_blending_id1_0_2;
      real_t jac_id1_0_1 = 1.0*jac_affine_0_1*jac_blending_id1_0_0 + 1.0*jac_affine_1_1*jac_blending_id1_0_1 + 1.0*jac_affine_2_1*jac_blending_id1_0_2;
      real_t jac_id1_2_1 = 1.0*jac_affine_0_1*jac_blending_id1_2_0 + 1.0*jac_affine_1_1*jac_blending_id1_2_1 + 1.0*jac_affine_2_1*jac_blending_id1_2_2;
      real_t tmp_17 = -1.0*jac_id1_0_0*jac_id1_2_1 + 1.0*jac_id1_0_1*jac_id1_2_0;
      real_t jac_id1_0_2 = 1.0*jac_affine_0_2*jac_blending_id1_0_0 + 1.0*jac_affine_1_2*jac_blending_id1_0_1 + 1.0*jac_affine_2_2*jac_blending_id1_0_2;
      real_t jac_id1_2_2 = 1.0*jac_affine_0_2*jac_blending_id1_2_0 + 1.0*jac_affine_1_2*jac_blending_id1_2_1 + 1.0*jac_affine_2_2*jac_blending_id1_2_2;
      real_t jac_id1_1_1 = 1.0*jac_affine_0_1*jac_blending_id1_1_0 + 1.0*jac_affine_1_1*jac_blending_id1_1_1 + 1.0*jac_affine_2_1*jac_blending_id1_1_2;
      real_t jac_id1_1_0 = 1.0*jac_affine_0_0*jac_blending_id1_1_0 + 1.0*jac_affine_1_0*jac_blending_id1_1_1 + 1.0*jac_affine_2_0*jac_blending_id1_1_2;
      real_t tmp_18 = 1.0*jac_id1_1_0*jac_id1_2_1 - 1.0*jac_id1_1_1*jac_id1_2_0;
      real_t tmp_16 = 1.0*jac_id1_0_0*jac_id1_1_1 - 1.0*jac_id1_0_1*jac_id1_1_0;
      real_t jac_id1_1_2 = 1.0*jac_affine_0_2*jac_blending_id1_1_0 + 1.0*jac_affine_1_2*jac_blending_id1_1_1 + 1.0*jac_affine_2_2*jac_blending_id1_1_2;
      real_t tmp_6 = jac_id1_0_2*tmp_18 + jac_id1_1_2*tmp_17 + jac_id1_2_2*tmp_16;
      real_t np_2_arg_0 = tmp_6;
      real_t np_2 = std::abs(np_2_arg_0);
      real_t abs_det_jac_id1 = 1.0*np_2;
      real_t tmp_55 = 1.0*jac_id1_1_1*jac_id1_2_2 - 1.0*jac_id1_1_2*jac_id1_2_1;
      real_t np_3_arg_0 = tmp_6;
      real_t np_3 = 1.0 / (np_3_arg_0);
      real_t jac_inv_id1_0_0 = np_3*tmp_55;
      real_t tmp_56 = -1.0*jac_id1_0_1*jac_id1_2_2 + 1.0*jac_id1_0_2*jac_id1_2_1;
      real_t jac_inv_id1_0_1 = np_3*tmp_56;
      real_t tmp_57 = 1.0*jac_id1_0_1*jac_id1_1_2 - 1.0*jac_id1_0_2*jac_id1_1_1;
      real_t jac_inv_id1_0_2 = np_3*tmp_57;
      real_t tmp_58 = -1.0*jac_id1_1_0*jac_id1_2_2 + 1.0*jac_id1_1_2*jac_id1_2_0;
      real_t jac_inv_id1_1_0 = np_3*tmp_58;
      real_t tmp_59 = 1.0*jac_id1_0_0*jac_id1_2_2 - 1.0*jac_id1_0_2*jac_id1_2_0;
      real_t jac_inv_id1_1_1 = np_3*tmp_59;
      real_t tmp_60 = -1.0*jac_id1_0_0*jac_id1_1_2 + 1.0*jac_id1_0_2*jac_id1_1_0;
      real_t jac_inv_id1_1_2 = np_3*tmp_60;
      real_t jac_inv_id1_2_0 = np_3*tmp_18;
      real_t jac_inv_id1_2_1 = np_3*tmp_17;
      real_t jac_inv_id1_2_2 = np_3*tmp_16;
      real_t jac_id2_1_1 = 1.0*jac_affine_0_1*jac_blending_id2_1_0 + 1.0*jac_affine_1_1*jac_blending_id2_1_1 + 1.0*jac_affine_2_1*jac_blending_id2_1_2;
      real_t jac_id2_1_0 = 1.0*jac_affine_0_0*jac_blending_id2_1_0 + 1.0*jac_affine_1_0*jac_blending_id2_1_1 + 1.0*jac_affine_2_0*jac_blending_id2_1_2;
      real_t jac_id2_0_1 = 1.0*jac_affine_0_1*jac_blending_id2_0_0 + 1.0*jac_affine_1_1*jac_blending_id2_0_1 + 1.0*jac_affine_2_1*jac_blending_id2_0_2;
      real_t jac_id2_0_0 = 1.0*jac_affine_0_0*jac_blending_id2_0_0 + 1.0*jac_affine_1_0*jac_blending_id2_0_1 + 1.0*jac_affine_2_0*jac_blending_id2_0_2;
      real_t tmp_13 = 1.0*jac_id2_0_0*jac_id2_1_1 - 1.0*jac_id2_0_1*jac_id2_1_0;
      real_t jac_id2_2_2 = 1.0*jac_affine_0_2*jac_blending_id2_2_0 + 1.0*jac_affine_1_2*jac_blending_id2_2_1 + 1.0*jac_affine_2_2*jac_blending_id2_2_2;
      real_t jac_id2_0_2 = 1.0*jac_affine_0_2*jac_blending_id2_0_0 + 1.0*jac_affine_1_2*jac_blending_id2_0_1 + 1.0*jac_affine_2_2*jac_blending_id2_0_2;
      real_t jac_id2_1_2 = 1.0*jac_affine_0_2*jac_blending_id2_1_0 + 1.0*jac_affine_1_2*jac_blending_id2_1_1 + 1.0*jac_affine_2_2*jac_blending_id2_1_2;
      real_t jac_id2_2_1 = 1.0*jac_affine_0_1*jac_blending_id2_2_0 + 1.0*jac_affine_1_1*jac_blending_id2_2_1 + 1.0*jac_affine_2_1*jac_blending_id2_2_2;
      real_t jac_id2_2_0 = 1.0*jac_affine_0_0*jac_blending_id2_2_0 + 1.0*jac_affine_1_0*jac_blending_id2_2_1 + 1.0*jac_affine_2_0*jac_blending_id2_2_2;
      real_t tmp_15 = 1.0*jac_id2_1_0*jac_id2_2_1 - 1.0*jac_id2_1_1*jac_id2_2_0;
      real_t tmp_14 = -1.0*jac_id2_0_0*jac_id2_2_1 + 1.0*jac_id2_0_1*jac_id2_2_0;
      real_t tmp_9 = jac_id2_0_2*tmp_15 + jac_id2_1_2*tmp_14 + jac_id2_2_2*tmp_13;
      real_t np_4_arg_0 = tmp_9;
      real_t np_4 = std::abs(np_4_arg_0);
      real_t abs_det_jac_id2 = 1.0*np_4;
      real_t tmp_61 = 1.0*jac_id2_1_1*jac_id2_2_2 - 1.0*jac_id2_1_2*jac_id2_2_1;
      real_t np_5_arg_0 = tmp_9;
      real_t np_5 = 1.0 / (np_5_arg_0);
      real_t jac_inv_id2_0_0 = np_5*tmp_61;
      real_t tmp_62 = -1.0*jac_id2_0_1*jac_id2_2_2 + 1.0*jac_id2_0_2*jac_id2_2_1;
      real_t jac_inv_id2_0_1 = np_5*tmp_62;
      real_t tmp_63 = 1.0*jac_id2_0_1*jac_id2_1_2 - 1.0*jac_id2_0_2*jac_id2_1_1;
      real_t jac_inv_id2_0_2 = np_5*tmp_63;
      real_t tmp_64 = -1.0*jac_id2_1_0*jac_id2_2_2 + 1.0*jac_id2_1_2*jac_id2_2_0;
      real_t jac_inv_id2_1_0 = np_5*tmp_64;
      real_t tmp_65 = 1.0*jac_id2_0_0*jac_id2_2_2 - 1.0*jac_id2_0_2*jac_id2_2_0;
      real_t jac_inv_id2_1_1 = np_5*tmp_65;
      real_t tmp_66 = -1.0*jac_id2_0_0*jac_id2_1_2 + 1.0*jac_id2_0_2*jac_id2_1_0;
      real_t jac_inv_id2_1_2 = np_5*tmp_66;
      real_t jac_inv_id2_2_0 = np_5*tmp_15;
      real_t jac_inv_id2_2_1 = np_5*tmp_14;
      real_t jac_inv_id2_2_2 = np_5*tmp_13;
      real_t jac_id3_1_2 = 1.0*jac_affine_0_2*jac_blending_id3_1_0 + 1.0*jac_affine_1_2*jac_blending_id3_1_1 + 1.0*jac_affine_2_2*jac_blending_id3_1_2;
      real_t jac_id3_2_2 = 1.0*jac_affine_0_2*jac_blending_id3_2_0 + 1.0*jac_affine_1_2*jac_blending_id3_2_1 + 1.0*jac_affine_2_2*jac_blending_id3_2_2;
      real_t jac_id3_0_0 = 1.0*jac_affine_0_0*jac_blending_id3_0_0 + 1.0*jac_affine_1_0*jac_blending_id3_0_1 + 1.0*jac_affine_2_0*jac_blending_id3_0_2;
      real_t jac_id3_1_1 = 1.0*jac_affine_0_1*jac_blending_id3_1_0 + 1.0*jac_affine_1_1*jac_blending_id3_1_1 + 1.0*jac_affine_2_1*jac_blending_id3_1_2;
      real_t jac_id3_1_0 = 1.0*jac_affine_0_0*jac_blending_id3_1_0 + 1.0*jac_affine_1_0*jac_blending_id3_1_1 + 1.0*jac_affine_2_0*jac_blending_id3_1_2;
      real_t jac_id3_0_1 = 1.0*jac_affine_0_1*jac_blending_id3_0_0 + 1.0*jac_affine_1_1*jac_blending_id3_0_1 + 1.0*jac_affine_2_1*jac_blending_id3_0_2;
      real_t tmp_10 = 1.0*jac_id3_0_0*jac_id3_1_1 - 1.0*jac_id3_0_1*jac_id3_1_0;
      real_t jac_id3_2_1 = 1.0*jac_affine_0_1*jac_blending_id3_2_0 + 1.0*jac_affine_1_1*jac_blending_id3_2_1 + 1.0*jac_affine_2_1*jac_blending_id3_2_2;
      real_t jac_id3_2_0 = 1.0*jac_affine_0_0*jac_blending_id3_2_0 + 1.0*jac_affine_1_0*jac_blending_id3_2_1 + 1.0*jac_affine_2_0*jac_blending_id3_2_2;
      real_t tmp_12 = 1.0*jac_id3_1_0*jac_id3_2_1 - 1.0*jac_id3_1_1*jac_id3_2_0;
      real_t tmp_11 = -1.0*jac_id3_0_0*jac_id3_2_1 + 1.0*jac_id3_0_1*jac_id3_2_0;
      real_t jac_id3_0_2 = 1.0*jac_affine_0_2*jac_blending_id3_0_0 + 1.0*jac_affine_1_2*jac_blending_id3_0_1 + 1.0*jac_affine_2_2*jac_blending_id3_0_2;
      real_t tmp_8 = jac_id3_0_2*tmp_12 + jac_id3_1_2*tmp_11 + jac_id3_2_2*tmp_10;
      real_t np_6_arg_0 = tmp_8;
      real_t np_6 = std::abs(np_6_arg_0);
      real_t abs_det_jac_id3 = 1.0*np_6;
      real_t np_7_arg_0 = tmp_8;
      real_t np_7 = 1.0 / (np_7_arg_0);
      real_t tmp_67 = 1.0*jac_id3_1_1*jac_id3_2_2 - 1.0*jac_id3_1_2*jac_id3_2_1;
      real_t jac_inv_id3_0_0 = np_7*tmp_67;
      real_t tmp_68 = -1.0*jac_id3_0_1*jac_id3_2_2 + 1.0*jac_id3_0_2*jac_id3_2_1;
      real_t jac_inv_id3_0_1 = np_7*tmp_68;
      real_t tmp_69 = 1.0*jac_id3_0_1*jac_id3_1_2 - 1.0*jac_id3_0_2*jac_id3_1_1;
      real_t jac_inv_id3_0_2 = np_7*tmp_69;
      real_t tmp_70 = -1.0*jac_id3_1_0*jac_id3_2_2 + 1.0*jac_id3_1_2*jac_id3_2_0;
      real_t jac_inv_id3_1_0 = np_7*tmp_70;
      real_t tmp_71 = 1.0*jac_id3_0_0*jac_id3_2_2 - 1.0*jac_id3_0_2*jac_id3_2_0;
      real_t jac_inv_id3_1_1 = np_7*tmp_71;
      real_t tmp_72 = -1.0*jac_id3_0_0*jac_id3_1_2 + 1.0*jac_id3_0_2*jac_id3_1_0;
      real_t jac_inv_id3_1_2 = np_7*tmp_72;
      real_t jac_inv_id3_2_0 = np_7*tmp_12;
      real_t jac_inv_id3_2_1 = np_7*tmp_11;
      real_t jac_inv_id3_2_2 = np_7*tmp_10;
      real_t tmp_25 = (jac_inv_id1_0_0*jac_inv_id1_0_0) + (jac_inv_id1_0_1*jac_inv_id1_0_1) + (jac_inv_id1_0_2*jac_inv_id1_0_2);
      real_t tmp_23 = (jac_inv_id0_0_0*jac_inv_id0_0_0) + (jac_inv_id0_0_1*jac_inv_id0_0_1) + (jac_inv_id0_0_2*jac_inv_id0_0_2);
      real_t tmp_27 = (jac_inv_id2_0_0*jac_inv_id2_0_0) + (jac_inv_id2_0_1*jac_inv_id2_0_1) + (jac_inv_id2_0_2*jac_inv_id2_0_2);
      real_t tmp_29 = (jac_inv_id3_0_0*jac_inv_id3_0_0) + (jac_inv_id3_0_1*jac_inv_id3_0_1) + (jac_inv_id3_0_2*jac_inv_id3_0_2);
      real_t tmp_3 = 0.050086823222829389*abs_det_jac_id0*tmp_23 + 0.046462929447761279*abs_det_jac_id1*tmp_25 + 0.05318232258357912*abs_det_jac_id2*tmp_27 + 0.016934591412496786*abs_det_jac_id3*tmp_29;
      real_t tmp_35 = jac_inv_id2_1_0*jac_inv_id2_2_0 + jac_inv_id2_1_1*jac_inv_id2_2_1 + jac_inv_id2_1_2*jac_inv_id2_2_2;
      real_t tmp_37 = jac_inv_id3_1_0*jac_inv_id3_2_0 + jac_inv_id3_1_1*jac_inv_id3_2_1 + jac_inv_id3_1_2*jac_inv_id3_2_2;
      real_t tmp_31 = jac_inv_id1_1_0*jac_inv_id1_2_0 + jac_inv_id1_1_1*jac_inv_id1_2_1 + jac_inv_id1_1_2*jac_inv_id1_2_2;
      real_t tmp_33 = jac_inv_id0_1_0*jac_inv_id0_2_0 + jac_inv_id0_1_1*jac_inv_id0_2_1 + jac_inv_id0_1_2*jac_inv_id0_2_2;
      real_t tmp_0 = 0.050086823222829389*abs_det_jac_id0*tmp_33 + 0.046462929447761279*abs_det_jac_id1*tmp_31 + 0.05318232258357912*abs_det_jac_id2*tmp_35 + 0.016934591412496786*abs_det_jac_id3*tmp_37;
      real_t tmp_41 = jac_inv_id1_0_0*jac_inv_id1_2_0 + jac_inv_id1_0_1*jac_inv_id1_2_1 + jac_inv_id1_0_2*jac_inv_id1_2_2;
      real_t tmp_39 = jac_inv_id0_0_0*jac_inv_id0_2_0 + jac_inv_id0_0_1*jac_inv_id0_2_1 + jac_inv_id0_0_2*jac_inv_id0_2_2;
      real_t tmp_43 = jac_inv_id2_0_0*jac_inv_id2_2_0 + jac_inv_id2_0_1*jac_inv_id2_2_1 + jac_inv_id2_0_2*jac_inv_id2_2_2;
      real_t tmp_45 = jac_inv_id3_0_0*jac_inv_id3_2_0 + jac_inv_id3_0_1*jac_inv_id3_2_1 + jac_inv_id3_0_2*jac_inv_id3_2_2;
      real_t tmp_1 = 0.050086823222829389*abs_det_jac_id0*tmp_39 + 0.046462929447761279*abs_det_jac_id1*tmp_41 + 0.05318232258357912*abs_det_jac_id2*tmp_43 + 0.016934591412496786*abs_det_jac_id3*tmp_45;
      real_t tmp_36 = jac_inv_id3_0_0*jac_inv_id3_1_0 + jac_inv_id3_0_1*jac_inv_id3_1_1 + jac_inv_id3_0_2*jac_inv_id3_1_2;
      real_t tmp_32 = jac_inv_id1_0_0*jac_inv_id1_1_0 + jac_inv_id1_0_1*jac_inv_id1_1_1 + jac_inv_id1_0_2*jac_inv_id1_1_2;
      real_t tmp_30 = jac_inv_id0_0_0*jac_inv_id0_1_0 + jac_inv_id0_0_1*jac_inv_id0_1_1 + jac_inv_id0_0_2*jac_inv_id0_1_2;
      real_t tmp_34 = jac_inv_id2_0_0*jac_inv_id2_1_0 + jac_inv_id2_0_1*jac_inv_id2_1_1 + jac_inv_id2_0_2*jac_inv_id2_1_2;
      real_t tmp_2 = 0.050086823222829389*abs_det_jac_id0*tmp_30 + 0.046462929447761279*abs_det_jac_id1*tmp_32 + 0.05318232258357912*abs_det_jac_id2*tmp_34 + 0.016934591412496786*abs_det_jac_id3*tmp_36;
      real_t tmp_46 = tmp_1 + tmp_2;
      real_t tmp_74 = tmp_0 + tmp_46;
      real_t tmp_44 = (jac_inv_id3_1_0*jac_inv_id3_1_0) + (jac_inv_id3_1_1*jac_inv_id3_1_1) + (jac_inv_id3_1_2*jac_inv_id3_1_2);
      real_t tmp_38 = (jac_inv_id0_1_0*jac_inv_id0_1_0) + (jac_inv_id0_1_1*jac_inv_id0_1_1) + (jac_inv_id0_1_2*jac_inv_id0_1_2);
      real_t tmp_42 = (jac_inv_id2_1_0*jac_inv_id2_1_0) + (jac_inv_id2_1_1*jac_inv_id2_1_1) + (jac_inv_id2_1_2*jac_inv_id2_1_2);
      real_t tmp_40 = (jac_inv_id1_1_0*jac_inv_id1_1_0) + (jac_inv_id1_1_1*jac_inv_id1_1_1) + (jac_inv_id1_1_2*jac_inv_id1_1_2);
      real_t tmp_5 = 0.050086823222829389*abs_det_jac_id0*tmp_38 + 0.046462929447761279*abs_det_jac_id1*tmp_40 + 0.05318232258357912*abs_det_jac_id2*tmp_42 + 0.016934591412496786*abs_det_jac_id3*tmp_44;
      real_t tmp_28 = (jac_inv_id3_2_0*jac_inv_id3_2_0) + (jac_inv_id3_2_1*jac_inv_id3_2_1) + (jac_inv_id3_2_2*jac_inv_id3_2_2);
      real_t tmp_22 = (jac_inv_id0_2_0*jac_inv_id0_2_0) + (jac_inv_id0_2_1*jac_inv_id0_2_1) + (jac_inv_id0_2_2*jac_inv_id0_2_2);
      real_t tmp_26 = (jac_inv_id2_2_0*jac_inv_id2_2_0) + (jac_inv_id2_2_1*jac_inv_id2_2_1) + (jac_inv_id2_2_2*jac_inv_id2_2_2);
      real_t tmp_24 = (jac_inv_id1_2_0*jac_inv_id1_2_0) + (jac_inv_id1_2_1*jac_inv_id1_2_1) + (jac_inv_id1_2_2*jac_inv_id1_2_2);
      real_t tmp_4 = 0.050086823222829389*abs_det_jac_id0*tmp_22 + 0.046462929447761279*abs_det_jac_id1*tmp_24 + 0.05318232258357912*abs_det_jac_id2*tmp_26 + 0.016934591412496786*abs_det_jac_id3*tmp_28;
      real_t a_0_0 = 1.0*tmp_3 + 1.0*tmp_4 + 1.0*tmp_5 + 2.0*tmp_74;
      real_t tmp_73 = tmp_3 + tmp_46;
      real_t a_0_1 = -1.0*tmp_73;
      real_t tmp_47 = tmp_0 + tmp_2 + tmp_5;
      real_t a_0_2 = -1.0*tmp_47;
      real_t tmp_48 = tmp_0 + tmp_1 + tmp_4;
      real_t a_0_3 = -1.0*tmp_48;
      real_t a_1_0 = -1.0*tmp_73;
      real_t a_1_1 = 1.0*tmp_3;
      real_t a_1_2 = 1.0*tmp_2;
      real_t a_1_3 = 1.0*tmp_1;
      real_t a_2_0 = -1.0*tmp_47;
      real_t a_2_1 = 1.0*tmp_2;
      real_t a_2_2 = 1.0*tmp_5;
      real_t a_2_3 = 1.0*tmp_0;
      real_t a_3_0 = -1.0*tmp_48;
      real_t a_3_1 = 1.0*tmp_1;
      real_t a_3_2 = 1.0*tmp_0;
      real_t a_3_3 = 1.0*tmp_4;
      elMat(0,0) = a_0_0;
      elMat(0,1) = a_0_1;
      elMat(0,2) = a_0_2;
      elMat(0,3) = a_0_3;
      elMat(1,0) = a_1_0;
      elMat(1,1) = a_1_1;
      elMat(1,2) = a_1_2;
      elMat(1,3) = a_1_3;
      elMat(2,0) = a_2_0;
      elMat(2,1) = a_2_1;
      elMat(2,2) = a_2_2;
      elMat(2,3) = a_2_3;
      elMat(3,0) = a_3_0;
      elMat(3,1) = a_3_1;
      elMat(3,2) = a_3_2;
      elMat(3,3) = a_3_3;
   }

   void p1_diffusion_blending_q2::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
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
      real_t jac_id0_1_2 = 1.0*jac_affine_0_2*jac_blending_id0_1_0 + 1.0*jac_affine_1_2*jac_blending_id0_1_1 + 1.0*jac_affine_2_2*jac_blending_id0_1_2;
      real_t jac_id0_2_1 = 1.0*jac_affine_0_1*jac_blending_id0_2_0 + 1.0*jac_affine_1_1*jac_blending_id0_2_1 + 1.0*jac_affine_2_1*jac_blending_id0_2_2;
      real_t jac_id0_0_1 = 1.0*jac_affine_0_1*jac_blending_id0_0_0 + 1.0*jac_affine_1_1*jac_blending_id0_0_1 + 1.0*jac_affine_2_1*jac_blending_id0_0_2;
      real_t jac_id0_2_0 = 1.0*jac_affine_0_0*jac_blending_id0_2_0 + 1.0*jac_affine_1_0*jac_blending_id0_2_1 + 1.0*jac_affine_2_0*jac_blending_id0_2_2;
      real_t jac_id0_0_0 = 1.0*jac_affine_0_0*jac_blending_id0_0_0 + 1.0*jac_affine_1_0*jac_blending_id0_0_1 + 1.0*jac_affine_2_0*jac_blending_id0_0_2;
      real_t tmp_21 = -1.0*jac_id0_0_0*jac_id0_2_1 + 1.0*jac_id0_0_1*jac_id0_2_0;
      real_t jac_id0_0_2 = 1.0*jac_affine_0_2*jac_blending_id0_0_0 + 1.0*jac_affine_1_2*jac_blending_id0_0_1 + 1.0*jac_affine_2_2*jac_blending_id0_0_2;
      real_t jac_id0_1_1 = 1.0*jac_affine_0_1*jac_blending_id0_1_0 + 1.0*jac_affine_1_1*jac_blending_id0_1_1 + 1.0*jac_affine_2_1*jac_blending_id0_1_2;
      real_t jac_id0_1_0 = 1.0*jac_affine_0_0*jac_blending_id0_1_0 + 1.0*jac_affine_1_0*jac_blending_id0_1_1 + 1.0*jac_affine_2_0*jac_blending_id0_1_2;
      real_t tmp_22 = 1.0*jac_id0_1_0*jac_id0_2_1 - 1.0*jac_id0_1_1*jac_id0_2_0;
      real_t jac_id0_2_2 = 1.0*jac_affine_0_2*jac_blending_id0_2_0 + 1.0*jac_affine_1_2*jac_blending_id0_2_1 + 1.0*jac_affine_2_2*jac_blending_id0_2_2;
      real_t tmp_20 = 1.0*jac_id0_0_0*jac_id0_1_1 - 1.0*jac_id0_0_1*jac_id0_1_0;
      real_t tmp_10 = jac_id0_0_2*tmp_22 + jac_id0_1_2*tmp_21 + jac_id0_2_2*tmp_20;
      real_t np_0_arg_0 = tmp_10;
      real_t np_0 = std::abs(np_0_arg_0);
      real_t abs_det_jac_id0 = 1.0*np_0;
      real_t np_1_arg_0 = tmp_10;
      real_t np_1 = 1.0 / (np_1_arg_0);
      real_t tmp_47 = 1.0*jac_id0_1_1*jac_id0_2_2 - 1.0*jac_id0_1_2*jac_id0_2_1;
      real_t jac_inv_id0_0_0 = np_1*tmp_47;
      real_t tmp_48 = -1.0*jac_id0_0_1*jac_id0_2_2 + 1.0*jac_id0_0_2*jac_id0_2_1;
      real_t jac_inv_id0_0_1 = np_1*tmp_48;
      real_t tmp_49 = 1.0*jac_id0_0_1*jac_id0_1_2 - 1.0*jac_id0_0_2*jac_id0_1_1;
      real_t jac_inv_id0_0_2 = np_1*tmp_49;
      real_t tmp_50 = -1.0*jac_id0_1_0*jac_id0_2_2 + 1.0*jac_id0_1_2*jac_id0_2_0;
      real_t jac_inv_id0_1_0 = np_1*tmp_50;
      real_t tmp_51 = 1.0*jac_id0_0_0*jac_id0_2_2 - 1.0*jac_id0_0_2*jac_id0_2_0;
      real_t jac_inv_id0_1_1 = np_1*tmp_51;
      real_t tmp_52 = -1.0*jac_id0_0_0*jac_id0_1_2 + 1.0*jac_id0_0_2*jac_id0_1_0;
      real_t jac_inv_id0_1_2 = np_1*tmp_52;
      real_t jac_inv_id0_2_0 = np_1*tmp_22;
      real_t jac_inv_id0_2_1 = np_1*tmp_21;
      real_t jac_inv_id0_2_2 = np_1*tmp_20;
      real_t jac_id1_1_1 = 1.0*jac_affine_0_1*jac_blending_id1_1_0 + 1.0*jac_affine_1_1*jac_blending_id1_1_1 + 1.0*jac_affine_2_1*jac_blending_id1_1_2;
      real_t jac_id1_0_0 = 1.0*jac_affine_0_0*jac_blending_id1_0_0 + 1.0*jac_affine_1_0*jac_blending_id1_0_1 + 1.0*jac_affine_2_0*jac_blending_id1_0_2;
      real_t jac_id1_1_0 = 1.0*jac_affine_0_0*jac_blending_id1_1_0 + 1.0*jac_affine_1_0*jac_blending_id1_1_1 + 1.0*jac_affine_2_0*jac_blending_id1_1_2;
      real_t jac_id1_0_1 = 1.0*jac_affine_0_1*jac_blending_id1_0_0 + 1.0*jac_affine_1_1*jac_blending_id1_0_1 + 1.0*jac_affine_2_1*jac_blending_id1_0_2;
      real_t tmp_17 = 1.0*jac_id1_0_0*jac_id1_1_1 - 1.0*jac_id1_0_1*jac_id1_1_0;
      real_t jac_id1_0_2 = 1.0*jac_affine_0_2*jac_blending_id1_0_0 + 1.0*jac_affine_1_2*jac_blending_id1_0_1 + 1.0*jac_affine_2_2*jac_blending_id1_0_2;
      real_t jac_id1_2_2 = 1.0*jac_affine_0_2*jac_blending_id1_2_0 + 1.0*jac_affine_1_2*jac_blending_id1_2_1 + 1.0*jac_affine_2_2*jac_blending_id1_2_2;
      real_t jac_id1_2_0 = 1.0*jac_affine_0_0*jac_blending_id1_2_0 + 1.0*jac_affine_1_0*jac_blending_id1_2_1 + 1.0*jac_affine_2_0*jac_blending_id1_2_2;
      real_t jac_id1_2_1 = 1.0*jac_affine_0_1*jac_blending_id1_2_0 + 1.0*jac_affine_1_1*jac_blending_id1_2_1 + 1.0*jac_affine_2_1*jac_blending_id1_2_2;
      real_t tmp_18 = -1.0*jac_id1_0_0*jac_id1_2_1 + 1.0*jac_id1_0_1*jac_id1_2_0;
      real_t tmp_19 = 1.0*jac_id1_1_0*jac_id1_2_1 - 1.0*jac_id1_1_1*jac_id1_2_0;
      real_t jac_id1_1_2 = 1.0*jac_affine_0_2*jac_blending_id1_1_0 + 1.0*jac_affine_1_2*jac_blending_id1_1_1 + 1.0*jac_affine_2_2*jac_blending_id1_1_2;
      real_t tmp_9 = jac_id1_0_2*tmp_19 + jac_id1_1_2*tmp_18 + jac_id1_2_2*tmp_17;
      real_t np_2_arg_0 = tmp_9;
      real_t np_2 = std::abs(np_2_arg_0);
      real_t abs_det_jac_id1 = 1.0*np_2;
      real_t tmp_53 = 1.0*jac_id1_1_1*jac_id1_2_2 - 1.0*jac_id1_1_2*jac_id1_2_1;
      real_t np_3_arg_0 = tmp_9;
      real_t np_3 = 1.0 / (np_3_arg_0);
      real_t jac_inv_id1_0_0 = np_3*tmp_53;
      real_t tmp_54 = -1.0*jac_id1_0_1*jac_id1_2_2 + 1.0*jac_id1_0_2*jac_id1_2_1;
      real_t jac_inv_id1_0_1 = np_3*tmp_54;
      real_t tmp_55 = 1.0*jac_id1_0_1*jac_id1_1_2 - 1.0*jac_id1_0_2*jac_id1_1_1;
      real_t jac_inv_id1_0_2 = np_3*tmp_55;
      real_t tmp_56 = -1.0*jac_id1_1_0*jac_id1_2_2 + 1.0*jac_id1_1_2*jac_id1_2_0;
      real_t jac_inv_id1_1_0 = np_3*tmp_56;
      real_t tmp_57 = 1.0*jac_id1_0_0*jac_id1_2_2 - 1.0*jac_id1_0_2*jac_id1_2_0;
      real_t jac_inv_id1_1_1 = np_3*tmp_57;
      real_t tmp_58 = -1.0*jac_id1_0_0*jac_id1_1_2 + 1.0*jac_id1_0_2*jac_id1_1_0;
      real_t jac_inv_id1_1_2 = np_3*tmp_58;
      real_t jac_inv_id1_2_0 = np_3*tmp_19;
      real_t jac_inv_id1_2_1 = np_3*tmp_18;
      real_t jac_inv_id1_2_2 = np_3*tmp_17;
      real_t jac_id2_2_1 = 1.0*jac_affine_0_1*jac_blending_id2_2_0 + 1.0*jac_affine_1_1*jac_blending_id2_2_1 + 1.0*jac_affine_2_1*jac_blending_id2_2_2;
      real_t jac_id2_1_1 = 1.0*jac_affine_0_1*jac_blending_id2_1_0 + 1.0*jac_affine_1_1*jac_blending_id2_1_1 + 1.0*jac_affine_2_1*jac_blending_id2_1_2;
      real_t jac_id2_1_0 = 1.0*jac_affine_0_0*jac_blending_id2_1_0 + 1.0*jac_affine_1_0*jac_blending_id2_1_1 + 1.0*jac_affine_2_0*jac_blending_id2_1_2;
      real_t jac_id2_2_0 = 1.0*jac_affine_0_0*jac_blending_id2_2_0 + 1.0*jac_affine_1_0*jac_blending_id2_2_1 + 1.0*jac_affine_2_0*jac_blending_id2_2_2;
      real_t tmp_16 = 1.0*jac_id2_1_0*jac_id2_2_1 - 1.0*jac_id2_1_1*jac_id2_2_0;
      real_t jac_id2_2_2 = 1.0*jac_affine_0_2*jac_blending_id2_2_0 + 1.0*jac_affine_1_2*jac_blending_id2_2_1 + 1.0*jac_affine_2_2*jac_blending_id2_2_2;
      real_t jac_id2_0_2 = 1.0*jac_affine_0_2*jac_blending_id2_0_0 + 1.0*jac_affine_1_2*jac_blending_id2_0_1 + 1.0*jac_affine_2_2*jac_blending_id2_0_2;
      real_t jac_id2_1_2 = 1.0*jac_affine_0_2*jac_blending_id2_1_0 + 1.0*jac_affine_1_2*jac_blending_id2_1_1 + 1.0*jac_affine_2_2*jac_blending_id2_1_2;
      real_t jac_id2_0_1 = 1.0*jac_affine_0_1*jac_blending_id2_0_0 + 1.0*jac_affine_1_1*jac_blending_id2_0_1 + 1.0*jac_affine_2_1*jac_blending_id2_0_2;
      real_t jac_id2_0_0 = 1.0*jac_affine_0_0*jac_blending_id2_0_0 + 1.0*jac_affine_1_0*jac_blending_id2_0_1 + 1.0*jac_affine_2_0*jac_blending_id2_0_2;
      real_t tmp_15 = -1.0*jac_id2_0_0*jac_id2_2_1 + 1.0*jac_id2_0_1*jac_id2_2_0;
      real_t tmp_14 = 1.0*jac_id2_0_0*jac_id2_1_1 - 1.0*jac_id2_0_1*jac_id2_1_0;
      real_t tmp_8 = jac_id2_0_2*tmp_16 + jac_id2_1_2*tmp_15 + jac_id2_2_2*tmp_14;
      real_t np_4_arg_0 = tmp_8;
      real_t np_4 = std::abs(np_4_arg_0);
      real_t abs_det_jac_id2 = 1.0*np_4;
      real_t np_5_arg_0 = tmp_8;
      real_t np_5 = 1.0 / (np_5_arg_0);
      real_t tmp_59 = 1.0*jac_id2_1_1*jac_id2_2_2 - 1.0*jac_id2_1_2*jac_id2_2_1;
      real_t jac_inv_id2_0_0 = np_5*tmp_59;
      real_t tmp_60 = -1.0*jac_id2_0_1*jac_id2_2_2 + 1.0*jac_id2_0_2*jac_id2_2_1;
      real_t jac_inv_id2_0_1 = np_5*tmp_60;
      real_t tmp_61 = 1.0*jac_id2_0_1*jac_id2_1_2 - 1.0*jac_id2_0_2*jac_id2_1_1;
      real_t jac_inv_id2_0_2 = np_5*tmp_61;
      real_t tmp_62 = -1.0*jac_id2_1_0*jac_id2_2_2 + 1.0*jac_id2_1_2*jac_id2_2_0;
      real_t jac_inv_id2_1_0 = np_5*tmp_62;
      real_t tmp_63 = 1.0*jac_id2_0_0*jac_id2_2_2 - 1.0*jac_id2_0_2*jac_id2_2_0;
      real_t jac_inv_id2_1_1 = np_5*tmp_63;
      real_t tmp_64 = -1.0*jac_id2_0_0*jac_id2_1_2 + 1.0*jac_id2_0_2*jac_id2_1_0;
      real_t jac_inv_id2_1_2 = np_5*tmp_64;
      real_t jac_inv_id2_2_0 = np_5*tmp_16;
      real_t jac_inv_id2_2_1 = np_5*tmp_15;
      real_t jac_inv_id2_2_2 = np_5*tmp_14;
      real_t jac_id3_1_2 = 1.0*jac_affine_0_2*jac_blending_id3_1_0 + 1.0*jac_affine_1_2*jac_blending_id3_1_1 + 1.0*jac_affine_2_2*jac_blending_id3_1_2;
      real_t jac_id3_2_1 = 1.0*jac_affine_0_1*jac_blending_id3_2_0 + 1.0*jac_affine_1_1*jac_blending_id3_2_1 + 1.0*jac_affine_2_1*jac_blending_id3_2_2;
      real_t jac_id3_1_1 = 1.0*jac_affine_0_1*jac_blending_id3_1_0 + 1.0*jac_affine_1_1*jac_blending_id3_1_1 + 1.0*jac_affine_2_1*jac_blending_id3_1_2;
      real_t jac_id3_2_0 = 1.0*jac_affine_0_0*jac_blending_id3_2_0 + 1.0*jac_affine_1_0*jac_blending_id3_2_1 + 1.0*jac_affine_2_0*jac_blending_id3_2_2;
      real_t jac_id3_1_0 = 1.0*jac_affine_0_0*jac_blending_id3_1_0 + 1.0*jac_affine_1_0*jac_blending_id3_1_1 + 1.0*jac_affine_2_0*jac_blending_id3_1_2;
      real_t tmp_13 = 1.0*jac_id3_1_0*jac_id3_2_1 - 1.0*jac_id3_1_1*jac_id3_2_0;
      real_t jac_id3_2_2 = 1.0*jac_affine_0_2*jac_blending_id3_2_0 + 1.0*jac_affine_1_2*jac_blending_id3_2_1 + 1.0*jac_affine_2_2*jac_blending_id3_2_2;
      real_t jac_id3_0_0 = 1.0*jac_affine_0_0*jac_blending_id3_0_0 + 1.0*jac_affine_1_0*jac_blending_id3_0_1 + 1.0*jac_affine_2_0*jac_blending_id3_0_2;
      real_t jac_id3_0_1 = 1.0*jac_affine_0_1*jac_blending_id3_0_0 + 1.0*jac_affine_1_1*jac_blending_id3_0_1 + 1.0*jac_affine_2_1*jac_blending_id3_0_2;
      real_t tmp_12 = -1.0*jac_id3_0_0*jac_id3_2_1 + 1.0*jac_id3_0_1*jac_id3_2_0;
      real_t tmp_11 = 1.0*jac_id3_0_0*jac_id3_1_1 - 1.0*jac_id3_0_1*jac_id3_1_0;
      real_t jac_id3_0_2 = 1.0*jac_affine_0_2*jac_blending_id3_0_0 + 1.0*jac_affine_1_2*jac_blending_id3_0_1 + 1.0*jac_affine_2_2*jac_blending_id3_0_2;
      real_t tmp_7 = jac_id3_0_2*tmp_13 + jac_id3_1_2*tmp_12 + jac_id3_2_2*tmp_11;
      real_t np_6_arg_0 = tmp_7;
      real_t np_6 = std::abs(np_6_arg_0);
      real_t abs_det_jac_id3 = 1.0*np_6;
      real_t np_7_arg_0 = tmp_7;
      real_t np_7 = 1.0 / (np_7_arg_0);
      real_t tmp_65 = 1.0*jac_id3_1_1*jac_id3_2_2 - 1.0*jac_id3_1_2*jac_id3_2_1;
      real_t jac_inv_id3_0_0 = np_7*tmp_65;
      real_t tmp_66 = -1.0*jac_id3_0_1*jac_id3_2_2 + 1.0*jac_id3_0_2*jac_id3_2_1;
      real_t jac_inv_id3_0_1 = np_7*tmp_66;
      real_t tmp_67 = 1.0*jac_id3_0_1*jac_id3_1_2 - 1.0*jac_id3_0_2*jac_id3_1_1;
      real_t jac_inv_id3_0_2 = np_7*tmp_67;
      real_t tmp_68 = -1.0*jac_id3_1_0*jac_id3_2_2 + 1.0*jac_id3_1_2*jac_id3_2_0;
      real_t jac_inv_id3_1_0 = np_7*tmp_68;
      real_t tmp_69 = 1.0*jac_id3_0_0*jac_id3_2_2 - 1.0*jac_id3_0_2*jac_id3_2_0;
      real_t jac_inv_id3_1_1 = np_7*tmp_69;
      real_t tmp_70 = -1.0*jac_id3_0_0*jac_id3_1_2 + 1.0*jac_id3_0_2*jac_id3_1_0;
      real_t jac_inv_id3_1_2 = np_7*tmp_70;
      real_t jac_inv_id3_2_0 = np_7*tmp_13;
      real_t jac_inv_id3_2_1 = np_7*tmp_12;
      real_t jac_inv_id3_2_2 = np_7*tmp_11;
      real_t tmp_28 = (jac_inv_id2_1_0*jac_inv_id2_1_0) + (jac_inv_id2_1_1*jac_inv_id2_1_1) + (jac_inv_id2_1_2*jac_inv_id2_1_2);
      real_t tmp_26 = (jac_inv_id1_1_0*jac_inv_id1_1_0) + (jac_inv_id1_1_1*jac_inv_id1_1_1) + (jac_inv_id1_1_2*jac_inv_id1_1_2);
      real_t tmp_30 = (jac_inv_id3_1_0*jac_inv_id3_1_0) + (jac_inv_id3_1_1*jac_inv_id3_1_1) + (jac_inv_id3_1_2*jac_inv_id3_1_2);
      real_t tmp_24 = (jac_inv_id0_1_0*jac_inv_id0_1_0) + (jac_inv_id0_1_1*jac_inv_id0_1_1) + (jac_inv_id0_1_2*jac_inv_id0_1_2);
      real_t tmp_3 = 0.050086823222829389*abs_det_jac_id0*tmp_24 + 0.046462929447761279*abs_det_jac_id1*tmp_26 + 0.05318232258357912*abs_det_jac_id2*tmp_28 + 0.016934591412496786*abs_det_jac_id3*tmp_30;
      real_t tmp_35 = (jac_inv_id2_0_0*jac_inv_id2_0_0) + (jac_inv_id2_0_1*jac_inv_id2_0_1) + (jac_inv_id2_0_2*jac_inv_id2_0_2);
      real_t tmp_37 = (jac_inv_id3_0_0*jac_inv_id3_0_0) + (jac_inv_id3_0_1*jac_inv_id3_0_1) + (jac_inv_id3_0_2*jac_inv_id3_0_2);
      real_t tmp_31 = (jac_inv_id0_0_0*jac_inv_id0_0_0) + (jac_inv_id0_0_1*jac_inv_id0_0_1) + (jac_inv_id0_0_2*jac_inv_id0_0_2);
      real_t tmp_33 = (jac_inv_id1_0_0*jac_inv_id1_0_0) + (jac_inv_id1_0_1*jac_inv_id1_0_1) + (jac_inv_id1_0_2*jac_inv_id1_0_2);
      real_t tmp_2 = 0.050086823222829389*abs_det_jac_id0*tmp_31 + 0.046462929447761279*abs_det_jac_id1*tmp_33 + 0.05318232258357912*abs_det_jac_id2*tmp_35 + 0.016934591412496786*abs_det_jac_id3*tmp_37;
      real_t tmp_44 = (jac_inv_id2_2_0*jac_inv_id2_2_0) + (jac_inv_id2_2_1*jac_inv_id2_2_1) + (jac_inv_id2_2_2*jac_inv_id2_2_2);
      real_t tmp_42 = (jac_inv_id1_2_0*jac_inv_id1_2_0) + (jac_inv_id1_2_1*jac_inv_id1_2_1) + (jac_inv_id1_2_2*jac_inv_id1_2_2);
      real_t tmp_40 = (jac_inv_id0_2_0*jac_inv_id0_2_0) + (jac_inv_id0_2_1*jac_inv_id0_2_1) + (jac_inv_id0_2_2*jac_inv_id0_2_2);
      real_t tmp_46 = (jac_inv_id3_2_0*jac_inv_id3_2_0) + (jac_inv_id3_2_1*jac_inv_id3_2_1) + (jac_inv_id3_2_2*jac_inv_id3_2_2);
      real_t tmp_6 = 0.050086823222829389*abs_det_jac_id0*tmp_40 + 0.046462929447761279*abs_det_jac_id1*tmp_42 + 0.05318232258357912*abs_det_jac_id2*tmp_44 + 0.016934591412496786*abs_det_jac_id3*tmp_46;
      real_t tmp_36 = jac_inv_id2_0_0*jac_inv_id2_1_0 + jac_inv_id2_0_1*jac_inv_id2_1_1 + jac_inv_id2_0_2*jac_inv_id2_1_2;
      real_t tmp_32 = jac_inv_id0_0_0*jac_inv_id0_1_0 + jac_inv_id0_0_1*jac_inv_id0_1_1 + jac_inv_id0_0_2*jac_inv_id0_1_2;
      real_t tmp_38 = jac_inv_id3_0_0*jac_inv_id3_1_0 + jac_inv_id3_0_1*jac_inv_id3_1_1 + jac_inv_id3_0_2*jac_inv_id3_1_2;
      real_t tmp_34 = jac_inv_id1_0_0*jac_inv_id1_1_0 + jac_inv_id1_0_1*jac_inv_id1_1_1 + jac_inv_id1_0_2*jac_inv_id1_1_2;
      real_t tmp_1 = 0.050086823222829389*abs_det_jac_id0*tmp_32 + 0.046462929447761279*abs_det_jac_id1*tmp_34 + 0.05318232258357912*abs_det_jac_id2*tmp_36 + 0.016934591412496786*abs_det_jac_id3*tmp_38;
      real_t tmp_25 = jac_inv_id1_1_0*jac_inv_id1_2_0 + jac_inv_id1_1_1*jac_inv_id1_2_1 + jac_inv_id1_1_2*jac_inv_id1_2_2;
      real_t tmp_23 = jac_inv_id0_1_0*jac_inv_id0_2_0 + jac_inv_id0_1_1*jac_inv_id0_2_1 + jac_inv_id0_1_2*jac_inv_id0_2_2;
      real_t tmp_27 = jac_inv_id2_1_0*jac_inv_id2_2_0 + jac_inv_id2_1_1*jac_inv_id2_2_1 + jac_inv_id2_1_2*jac_inv_id2_2_2;
      real_t tmp_29 = jac_inv_id3_1_0*jac_inv_id3_2_0 + jac_inv_id3_1_1*jac_inv_id3_2_1 + jac_inv_id3_1_2*jac_inv_id3_2_2;
      real_t tmp_5 = 0.050086823222829389*abs_det_jac_id0*tmp_23 + 0.046462929447761279*abs_det_jac_id1*tmp_25 + 0.05318232258357912*abs_det_jac_id2*tmp_27 + 0.016934591412496786*abs_det_jac_id3*tmp_29;
      real_t tmp_41 = jac_inv_id1_0_0*jac_inv_id1_2_0 + jac_inv_id1_0_1*jac_inv_id1_2_1 + jac_inv_id1_0_2*jac_inv_id1_2_2;
      real_t tmp_39 = jac_inv_id0_0_0*jac_inv_id0_2_0 + jac_inv_id0_0_1*jac_inv_id0_2_1 + jac_inv_id0_0_2*jac_inv_id0_2_2;
      real_t tmp_43 = jac_inv_id2_0_0*jac_inv_id2_2_0 + jac_inv_id2_0_1*jac_inv_id2_2_1 + jac_inv_id2_0_2*jac_inv_id2_2_2;
      real_t tmp_45 = jac_inv_id3_0_0*jac_inv_id3_2_0 + jac_inv_id3_0_1*jac_inv_id3_2_1 + jac_inv_id3_0_2*jac_inv_id3_2_2;
      real_t tmp_4 = 0.050086823222829389*abs_det_jac_id0*tmp_39 + 0.046462929447761279*abs_det_jac_id1*tmp_41 + 0.05318232258357912*abs_det_jac_id2*tmp_43 + 0.016934591412496786*abs_det_jac_id3*tmp_45;
      real_t tmp_0 = tmp_4 + tmp_5;
      real_t tmp_71 = tmp_0 + tmp_1;
      real_t a_0_0 = 1.0*tmp_2 + 1.0*tmp_3 + 1.0*tmp_6 + 2.0*tmp_71;
      real_t a_0_1 = -1.0*tmp_1 - 1.0*tmp_2 - 1.0*tmp_4;
      real_t a_0_2 = -1.0*tmp_1 - 1.0*tmp_3 - 1.0*tmp_5;
      real_t a_0_3 = -1.0*tmp_0 - 1.0*tmp_6;
      elMat(0,0) = a_0_0;
      elMat(0,1) = a_0_1;
      elMat(0,2) = a_0_2;
      elMat(0,3) = a_0_3;
   }

   void p1_diffusion_blending_q2::Blending_DF_Tetrahedron_jac_blending( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
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

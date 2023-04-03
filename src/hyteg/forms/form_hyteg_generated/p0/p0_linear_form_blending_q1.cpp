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

#include "p0_linear_form_blending_q1.hpp"

namespace hyteg {
namespace forms {

   void p0_linear_form_blending_q1::integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 1 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Blending_DF_Triangle_blend_out0_id0 = 0;
      real_t Blending_DF_Triangle_blend_out1_id0 = 0;
      real_t Blending_DF_Triangle_blend_out2_id0 = 0;
      real_t Blending_DF_Triangle_blend_out3_id0 = 0;
      real_t Blending_F_Triangle_blend_out0_id1 = 0;
      real_t Blending_F_Triangle_blend_out1_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id0 = 0;
      Blending_DF_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_DF_Triangle_blend_out0_id0, &Blending_DF_Triangle_blend_out1_id0, &Blending_DF_Triangle_blend_out2_id0, &Blending_DF_Triangle_blend_out3_id0 );
      Blending_F_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_F_Triangle_blend_out0_id1, &Blending_F_Triangle_blend_out1_id1 );
      Scalar_Variable_Coefficient_2D_k( Blending_F_Triangle_blend_out0_id1, Blending_F_Triangle_blend_out1_id1, &Scalar_Variable_Coefficient_2D_k_out0_id0 );
      real_t a_0_0 = 0.5*Scalar_Variable_Coefficient_2D_k_out0_id0*std::abs(Blending_DF_Triangle_blend_out0_id0*Blending_DF_Triangle_blend_out3_id0 - Blending_DF_Triangle_blend_out1_id0*Blending_DF_Triangle_blend_out2_id0)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      (elMat(0, 0)) = a_0_0;
   }

   void p0_linear_form_blending_q1::integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 1 >& elMat ) const
   {
      real_t p_affine_0_0 = coords[0][0];
      real_t p_affine_0_1 = coords[0][1];
      real_t p_affine_1_0 = coords[1][0];
      real_t p_affine_1_1 = coords[1][1];
      real_t p_affine_2_0 = coords[2][0];
      real_t p_affine_2_1 = coords[2][1];
      real_t Blending_DF_Triangle_blend_out0_id0 = 0;
      real_t Blending_DF_Triangle_blend_out1_id0 = 0;
      real_t Blending_DF_Triangle_blend_out2_id0 = 0;
      real_t Blending_DF_Triangle_blend_out3_id0 = 0;
      real_t Blending_F_Triangle_blend_out0_id1 = 0;
      real_t Blending_F_Triangle_blend_out1_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_k_out0_id0 = 0;
      Blending_DF_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_DF_Triangle_blend_out0_id0, &Blending_DF_Triangle_blend_out1_id0, &Blending_DF_Triangle_blend_out2_id0, &Blending_DF_Triangle_blend_out3_id0 );
      Blending_F_Triangle_blend( 0.33333333333333343*p_affine_0_0 + 0.33333333333333331*p_affine_1_0 + 0.33333333333333331*p_affine_2_0, 0.33333333333333343*p_affine_0_1 + 0.33333333333333331*p_affine_1_1 + 0.33333333333333331*p_affine_2_1, &Blending_F_Triangle_blend_out0_id1, &Blending_F_Triangle_blend_out1_id1 );
      Scalar_Variable_Coefficient_2D_k( Blending_F_Triangle_blend_out0_id1, Blending_F_Triangle_blend_out1_id1, &Scalar_Variable_Coefficient_2D_k_out0_id0 );
      real_t a_0_0 = 0.5*Scalar_Variable_Coefficient_2D_k_out0_id0*std::abs(Blending_DF_Triangle_blend_out0_id0*Blending_DF_Triangle_blend_out3_id0 - Blending_DF_Triangle_blend_out1_id0*Blending_DF_Triangle_blend_out2_id0)*std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      (elMat(0, 0)) = a_0_0;
   }

   void p0_linear_form_blending_q1::integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 1 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_k( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_k_out0_id0 );
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
      real_t a_0_0 = 0.16666666666666663*Scalar_Variable_Coefficient_3D_k_out0_id0*std::abs(Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0 + Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0)*std::abs(p_affine_0_0*tmp_2 - p_affine_0_0*tmp_9 - p_affine_0_1*tmp_10 + p_affine_0_1*tmp_5 - p_affine_0_2*tmp_11 + p_affine_0_2*tmp_8 - p_affine_1_0*tmp_2 + p_affine_1_0*tmp_9 + p_affine_1_1*tmp_10 - p_affine_1_1*tmp_5 + p_affine_1_2*tmp_11 - p_affine_1_2*tmp_8 + p_affine_2_0*tmp_4 - p_affine_2_0*tmp_7 - p_affine_2_1*tmp_1 + p_affine_2_1*tmp_6 + p_affine_2_2*tmp_0 - p_affine_2_2*tmp_3 - p_affine_3_0*tmp_4 + p_affine_3_0*tmp_7 + p_affine_3_1*tmp_1 - p_affine_3_1*tmp_6 - p_affine_3_2*tmp_0 + p_affine_3_2*tmp_3);
      (elMat(0, 0)) = a_0_0;
   }

   void p0_linear_form_blending_q1::integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 1 >& elMat ) const
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
      real_t Blending_DF_Tetrahedron_blend_out0_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out1_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out2_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out3_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out4_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out5_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out6_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out7_id0 = 0;
      real_t Blending_DF_Tetrahedron_blend_out8_id0 = 0;
      real_t Blending_F_Tetrahedron_blend_out0_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out1_id1 = 0;
      real_t Blending_F_Tetrahedron_blend_out2_id1 = 0;
      real_t Scalar_Variable_Coefficient_3D_k_out0_id0 = 0;
      Blending_DF_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_DF_Tetrahedron_blend_out0_id0, &Blending_DF_Tetrahedron_blend_out1_id0, &Blending_DF_Tetrahedron_blend_out2_id0, &Blending_DF_Tetrahedron_blend_out3_id0, &Blending_DF_Tetrahedron_blend_out4_id0, &Blending_DF_Tetrahedron_blend_out5_id0, &Blending_DF_Tetrahedron_blend_out6_id0, &Blending_DF_Tetrahedron_blend_out7_id0, &Blending_DF_Tetrahedron_blend_out8_id0 );
      Blending_F_Tetrahedron_blend( 0.25*p_affine_0_0 + 0.25*p_affine_1_0 + 0.25*p_affine_2_0 + 0.25*p_affine_3_0, 0.25*p_affine_0_1 + 0.25*p_affine_1_1 + 0.25*p_affine_2_1 + 0.25*p_affine_3_1, 0.25*p_affine_0_2 + 0.25*p_affine_1_2 + 0.25*p_affine_2_2 + 0.25*p_affine_3_2, &Blending_F_Tetrahedron_blend_out0_id1, &Blending_F_Tetrahedron_blend_out1_id1, &Blending_F_Tetrahedron_blend_out2_id1 );
      Scalar_Variable_Coefficient_3D_k( Blending_F_Tetrahedron_blend_out0_id1, Blending_F_Tetrahedron_blend_out1_id1, Blending_F_Tetrahedron_blend_out2_id1, &Scalar_Variable_Coefficient_3D_k_out0_id0 );
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
      real_t a_0_0 = 0.16666666666666663*Scalar_Variable_Coefficient_3D_k_out0_id0*std::abs(Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out8_id0 - Blending_DF_Tetrahedron_blend_out0_id0*Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out8_id0 + Blending_DF_Tetrahedron_blend_out1_id0*Blending_DF_Tetrahedron_blend_out5_id0*Blending_DF_Tetrahedron_blend_out6_id0 + Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out3_id0*Blending_DF_Tetrahedron_blend_out7_id0 - Blending_DF_Tetrahedron_blend_out2_id0*Blending_DF_Tetrahedron_blend_out4_id0*Blending_DF_Tetrahedron_blend_out6_id0)*std::abs(p_affine_0_0*tmp_2 - p_affine_0_0*tmp_9 - p_affine_0_1*tmp_10 + p_affine_0_1*tmp_5 - p_affine_0_2*tmp_11 + p_affine_0_2*tmp_8 - p_affine_1_0*tmp_2 + p_affine_1_0*tmp_9 + p_affine_1_1*tmp_10 - p_affine_1_1*tmp_5 + p_affine_1_2*tmp_11 - p_affine_1_2*tmp_8 + p_affine_2_0*tmp_4 - p_affine_2_0*tmp_7 - p_affine_2_1*tmp_1 + p_affine_2_1*tmp_6 + p_affine_2_2*tmp_0 - p_affine_2_2*tmp_3 - p_affine_3_0*tmp_4 + p_affine_3_0*tmp_7 + p_affine_3_1*tmp_1 - p_affine_3_1*tmp_6 - p_affine_3_2*tmp_0 + p_affine_3_2*tmp_3);
      (elMat(0, 0)) = a_0_0;
   }

   void p0_linear_form_blending_q1::Blending_DF_Triangle_blend( real_t in_0, real_t in_1, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3 ) const
   {
      Point3D  mappedPt( {in_0, in_1, 0} );
      Matrix2r DPsi;
      geometryMap_->evalDF( mappedPt, DPsi );
      *out_0 = DPsi( 0, 0 );
      *out_1 = DPsi( 0, 1 );
      *out_2 = DPsi( 1, 0 );
      *out_3 = DPsi( 1, 1 );
   }

   void p0_linear_form_blending_q1::Blending_F_Triangle_blend( real_t in_0, real_t in_1, real_t * out_0, real_t * out_1 ) const
   {
      Point3D  in( {in_0, in_1, 0} );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
   }

   void p0_linear_form_blending_q1::Scalar_Variable_Coefficient_2D_k( real_t in_0, real_t in_1, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_k( Point3D( {in_0, in_1, 0} ) );
   }

   void p0_linear_form_blending_q1::Blending_DF_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2, real_t * out_3, real_t * out_4, real_t * out_5, real_t * out_6, real_t * out_7, real_t * out_8 ) const
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

   void p0_linear_form_blending_q1::Blending_F_Tetrahedron_blend( real_t in_0, real_t in_1, real_t in_2, real_t * out_0, real_t * out_1, real_t * out_2 ) const
   {
      Point3D  in( {in_0, in_1, in_2} );
      Point3D out;
      geometryMap_->evalF( in, out );
      *out_0 = out[0];
      *out_1 = out[1];
      *out_2 = out[2];
   }

   void p0_linear_form_blending_q1::Scalar_Variable_Coefficient_3D_k( real_t in_0, real_t in_1, real_t in_2, real_t * out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_k( Point3D( {in_0, in_1, in_2} ) );
   }

} // namespace forms
} // namespace hyteg

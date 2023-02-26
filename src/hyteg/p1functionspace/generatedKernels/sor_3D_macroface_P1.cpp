/*
 * Copyright (c) 2019 Nils Kohl, Dominik Thoennes.
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

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "sor_3D_macroface_P1.hpp"

namespace hyteg {
namespace vertexdof {
namespace macroface {
namespace generated {

void sor_3D_macroface_P1(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double * RESTRICT _data_vertexFaceDst_gl1, double const * RESTRICT const _data_vertexFaceRhs, int level, int64_t neighbor_cell_0_local_vertex_id_0, int64_t neighbor_cell_0_local_vertex_id_1, int64_t neighbor_cell_0_local_vertex_id_2, int64_t neighbor_cell_1_local_vertex_id_0, int64_t neighbor_cell_1_local_vertex_id_1, int64_t neighbor_cell_1_local_vertex_id_2, double relax, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_0, std::map< hyteg::indexing::Index, double > v2v_cell_stencil_fused_face_1)
{
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_012_012(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_012_013(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_012_021(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_012_023(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_012_031(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_012_032(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_012_102(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_012_103(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_012_120(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_012_123(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_012_130(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_012_132(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_012_201(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_012_203(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_012_210(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_012_213(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_012_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_012_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_012_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_012_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_012_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_012_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_012_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_012_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_013_013(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_013_021(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_013_023(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_013_031(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_013_032(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_013_102(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_013_103(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_013_120(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_013_123(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_013_130(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_013_132(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_013_201(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_013_203(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_013_210(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_013_213(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_013_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_013_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_013_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_013_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_013_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_013_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_013_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_013_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_021_021(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_021_023(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_021_031(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_021_032(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_021_102(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_021_103(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_021_120(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_021_123(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_021_130(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_021_132(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_021_201(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_021_203(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_021_210(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_021_213(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_021_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_021_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_021_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_021_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_021_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_021_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_021_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_021_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_023_023(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_023_031(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_023_032(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_023_102(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_023_103(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_023_120(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_023_123(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_023_130(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_023_132(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_023_201(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_023_203(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_023_210(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_023_213(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_023_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_023_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_023_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_023_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_023_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_023_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_023_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_023_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_031_031(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_031_032(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_031_102(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_031_103(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_031_120(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_031_123(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_031_130(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_031_132(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_031_201(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_031_203(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_031_210(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_031_213(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_031_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_031_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_031_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_031_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_031_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_031_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_031_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_031_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_032_032(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_032_102(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_032_103(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_032_120(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_032_123(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_032_130(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_032_132(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_032_201(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_032_203(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_032_210(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_032_213(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_032_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_032_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_032_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_032_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_032_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_032_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_032_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_032_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_102_102(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_102_103(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_102_120(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_102_123(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_102_130(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_102_132(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_102_201(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_102_203(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_102_210(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_102_213(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_102_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_102_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_102_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_102_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_102_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_102_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_102_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_102_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_103_103(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_103_120(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_103_123(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_103_130(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_103_132(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_103_201(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_103_203(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_103_210(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_103_213(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_103_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_103_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_103_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_103_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_103_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_103_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_103_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_103_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_120_120(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_120_123(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_120_130(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_120_132(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_120_201(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_120_203(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_120_210(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_120_213(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_120_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_120_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_120_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_120_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_120_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_120_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_120_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_120_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_123_123(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_123_130(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_123_132(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_123_201(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_123_203(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_123_210(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_123_213(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_123_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_123_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_123_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_123_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_123_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_123_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_123_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_123_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_130_130(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_130_132(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_130_201(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_130_203(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_130_210(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_130_213(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_130_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_130_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_130_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_130_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_130_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_130_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_130_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_130_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_132_132(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_132_201(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_132_203(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_132_210(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_132_213(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_132_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_132_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_132_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_132_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_132_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_132_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_132_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_132_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_201_201(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_201_203(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_201_210(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_201_213(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_201_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_201_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_201_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_201_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_201_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_201_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_201_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_201_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_203_203(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_203_210(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_203_213(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_203_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_203_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_203_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_203_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_203_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_203_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_203_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_203_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_210_210(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_210_213(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_210_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_210_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_210_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_210_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_210_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_210_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_210_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_210_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_impl_213_213(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_213_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_213_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_213_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_213_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_213_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_213_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_213_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_213_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_230_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_230_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_230_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_230_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_230_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_230_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_230_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_230_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_impl_231_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_231_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_231_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_231_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_231_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_231_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_231_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_301_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_301_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_301_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_301_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_301_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_301_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_302_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_302_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_302_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_302_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_302_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_310_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_310_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_310_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_310_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_312_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_312_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_312_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_320_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_320_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_impl_321_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
      
      return;
   } 
}


} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg
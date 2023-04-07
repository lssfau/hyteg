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

#include "sor_3D_macroface_P2_update_edgedofs_one_sided.hpp"

namespace hyteg {
namespace P2 {
namespace macroface {
namespace generated {

void sor_3D_macroface_P2_update_edgedofs_one_sided(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceDst_gl0_X, double * RESTRICT _data_edgeFaceDst_gl0_XY, double * RESTRICT _data_edgeFaceDst_gl0_XYZ, double * RESTRICT _data_edgeFaceDst_gl0_XZ, double * RESTRICT _data_edgeFaceDst_gl0_Y, double * RESTRICT _data_edgeFaceDst_gl0_YZ, double * RESTRICT _data_edgeFaceDst_gl0_Z, double const * RESTRICT const _data_edgeFaceRhs_X, double const * RESTRICT const _data_edgeFaceRhs_XY, double const * RESTRICT const _data_edgeFaceRhs_Y, double const * RESTRICT const _data_vertexFaceDst, double const * RESTRICT const _data_vertexFaceDst_gl0, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil_fused_face_0, int level, int64_t neighbor_cell_0_local_vertex_id_0, int64_t neighbor_cell_0_local_vertex_id_1, int64_t neighbor_cell_0_local_vertex_id_2, double relax, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > v2e_cell_stencil_fused_face_0)
{
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_012(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_013(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_021(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_023(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_031(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_032(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_102(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_103(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_120(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_123(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_130(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_132(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_201(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_203(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_210(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_one_sided_impl_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, e2e_cell_stencil_fused_face_0, level, relax, v2e_cell_stencil_fused_face_0);
      
      return;
   } 
}


} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hyteg
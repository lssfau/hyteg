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

#include "apply_3D_macroface_one_sided_edgedof_to_edgedof_add.hpp"

namespace hyteg {
namespace edgedof {
namespace macroface {
namespace generated {

void apply_3D_macroface_one_sided_edgedof_to_edgedof_add(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * RESTRICT const _data_edgeFaceSrc_gl0_X, double const * RESTRICT const _data_edgeFaceSrc_gl0_XY, double const * RESTRICT const _data_edgeFaceSrc_gl0_XYZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_XZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_Y, double const * RESTRICT const _data_edgeFaceSrc_gl0_YZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_Z, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level, int64_t neighbor_cell_local_vertex_id_0, int64_t neighbor_cell_local_vertex_id_1, int64_t neighbor_cell_local_vertex_id_2)
{
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_012(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_013(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_021(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_023(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_031(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_032(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_102(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_103(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_120(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_123(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_130(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_132(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_201(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_203(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_210(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
      
      return;
   } 
}


} // namespace generated
} // namespace macroface
} // namespace edgedof
} // namespace hyteg
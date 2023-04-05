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

#include "sor_3D_macroface_P2_update_edgedofs.hpp"

namespace hyteg {
namespace P2 {
namespace macroface {
namespace generated {

void sor_3D_macroface_P2_update_edgedofs(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceDst_gl0_X, double * RESTRICT _data_edgeFaceDst_gl0_XY, double * RESTRICT _data_edgeFaceDst_gl0_XYZ, double * RESTRICT _data_edgeFaceDst_gl0_XZ, double * RESTRICT _data_edgeFaceDst_gl0_Y, double * RESTRICT _data_edgeFaceDst_gl0_YZ, double * RESTRICT _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_edgeFaceDst_gl1_X, double * RESTRICT _data_edgeFaceDst_gl1_XY, double * RESTRICT _data_edgeFaceDst_gl1_XYZ, double * RESTRICT _data_edgeFaceDst_gl1_XZ, double * RESTRICT _data_edgeFaceDst_gl1_Y, double * RESTRICT _data_edgeFaceDst_gl1_YZ, double * RESTRICT _data_edgeFaceDst_gl1_Z, double const * RESTRICT const _data_edgeFaceRhs_X, double const * RESTRICT const _data_edgeFaceRhs_XY, double const * RESTRICT const _data_edgeFaceRhs_Y, double const * RESTRICT const _data_vertexFaceDst, double const * RESTRICT const _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceDst_gl1, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil_fused_face_0, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil_fused_face_1, int level, int64_t neighbor_cell_0_local_vertex_id_0, int64_t neighbor_cell_0_local_vertex_id_1, int64_t neighbor_cell_0_local_vertex_id_2, int64_t neighbor_cell_1_local_vertex_id_0, int64_t neighbor_cell_1_local_vertex_id_1, int64_t neighbor_cell_1_local_vertex_id_2, double relax, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > v2e_cell_stencil_fused_face_0, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > v2e_cell_stencil_fused_face_1)
{
  double e2e_cell_stencil_fused_face_0_ptr[280];
  double e2e_cell_stencil_fused_face_1_ptr[280];
  double v2e_cell_stencil_fused_face_0_ptr[280];
  double v2e_cell_stencil_fused_face_1_ptr[280];
  e2e_cell_stencil_fused_face_0_ptr[0] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
  e2e_cell_stencil_fused_face_1_ptr[0] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
  e2e_cell_stencil_fused_face_0_ptr[1] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[1] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[2] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[2] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[3] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[3] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[4] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 1, -1 }];
  e2e_cell_stencil_fused_face_1_ptr[4] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 1, -1 }];
  e2e_cell_stencil_fused_face_0_ptr[5] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[5] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[6] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
  e2e_cell_stencil_fused_face_1_ptr[6] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
  e2e_cell_stencil_fused_face_0_ptr[7] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, -1 }];
  e2e_cell_stencil_fused_face_1_ptr[7] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, -1 }];
  e2e_cell_stencil_fused_face_0_ptr[8] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[8] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[9] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
  e2e_cell_stencil_fused_face_1_ptr[9] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
  e2e_cell_stencil_fused_face_0_ptr[10] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[10] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[11] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, -1 }];
  e2e_cell_stencil_fused_face_1_ptr[11] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, -1 }];
  e2e_cell_stencil_fused_face_0_ptr[12] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[12] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[13] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[13] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[14] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[14] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[15] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[15] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[16] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
  e2e_cell_stencil_fused_face_1_ptr[16] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
  e2e_cell_stencil_fused_face_0_ptr[17] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[17] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[18] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[18] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[19] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[19] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[20] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
  e2e_cell_stencil_fused_face_1_ptr[20] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
  e2e_cell_stencil_fused_face_0_ptr[21] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[21] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[22] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[22] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[23] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[23] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[24] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
  e2e_cell_stencil_fused_face_1_ptr[24] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
  e2e_cell_stencil_fused_face_0_ptr[25] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[25] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[26] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[26] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[27] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[27] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[28] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[28] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[29] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[29] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[30] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[30] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[31] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
  e2e_cell_stencil_fused_face_1_ptr[31] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
  e2e_cell_stencil_fused_face_0_ptr[32] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[32] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[33] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
  e2e_cell_stencil_fused_face_1_ptr[33] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
  e2e_cell_stencil_fused_face_0_ptr[34] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
  e2e_cell_stencil_fused_face_1_ptr[34] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
  e2e_cell_stencil_fused_face_0_ptr[35] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[35] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[36] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
  e2e_cell_stencil_fused_face_1_ptr[36] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
  e2e_cell_stencil_fused_face_0_ptr[37] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[37] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[38] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[38] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[39] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
  e2e_cell_stencil_fused_face_1_ptr[39] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
  e2e_cell_stencil_fused_face_0_ptr[40] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[40] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[41] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[41] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[42] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[42] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[43] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[43] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[44] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[44] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[45] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
  e2e_cell_stencil_fused_face_1_ptr[45] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
  e2e_cell_stencil_fused_face_0_ptr[46] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[46] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[47] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
  e2e_cell_stencil_fused_face_1_ptr[47] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
  e2e_cell_stencil_fused_face_0_ptr[48] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[48] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[49] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[49] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[50] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[50] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[51] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 1 }];
  e2e_cell_stencil_fused_face_1_ptr[51] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 1 }];
  e2e_cell_stencil_fused_face_0_ptr[52] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[52] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[53] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 1 }];
  e2e_cell_stencil_fused_face_1_ptr[53] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 1 }];
  e2e_cell_stencil_fused_face_0_ptr[54] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[54] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[55] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[55] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[56] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[56] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[57] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[57] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[58] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
  e2e_cell_stencil_fused_face_1_ptr[58] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
  e2e_cell_stencil_fused_face_0_ptr[59] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[59] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[60] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[60] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[61] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 1 }];
  e2e_cell_stencil_fused_face_1_ptr[61] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 1 }];
  e2e_cell_stencil_fused_face_0_ptr[62] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[62] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[63] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[63] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[64] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[64] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[65] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[65] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[66] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[66] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[67] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
  e2e_cell_stencil_fused_face_1_ptr[67] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
  e2e_cell_stencil_fused_face_0_ptr[68] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[68] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[69] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[69] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[70] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[70] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[71] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[71] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[72] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[72] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[73] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
  e2e_cell_stencil_fused_face_1_ptr[73] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
  e2e_cell_stencil_fused_face_0_ptr[74] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[74] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[75] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[75] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[76] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[76] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[77] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
  e2e_cell_stencil_fused_face_1_ptr[77] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
  e2e_cell_stencil_fused_face_0_ptr[78] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[78] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[79] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[79] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[80] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[80] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[81] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
  e2e_cell_stencil_fused_face_1_ptr[81] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
  e2e_cell_stencil_fused_face_0_ptr[82] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[82] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[83] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[83] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[84] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[84] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[85] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
  e2e_cell_stencil_fused_face_1_ptr[85] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
  e2e_cell_stencil_fused_face_0_ptr[86] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[86] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[87] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
  e2e_cell_stencil_fused_face_1_ptr[87] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
  e2e_cell_stencil_fused_face_0_ptr[88] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[88] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[89] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[89] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[90] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[90] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[91] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[91] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[92] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[92] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[93] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
  e2e_cell_stencil_fused_face_1_ptr[93] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
  e2e_cell_stencil_fused_face_0_ptr[94] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[94] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_0_ptr[95] = e2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
  e2e_cell_stencil_fused_face_1_ptr[95] = e2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[0] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[0] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[1] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 1 }];
  v2e_cell_stencil_fused_face_1_ptr[1] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 1 }];
  v2e_cell_stencil_fused_face_0_ptr[2] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 1, -1 }];
  v2e_cell_stencil_fused_face_1_ptr[2] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 1, -1 }];
  v2e_cell_stencil_fused_face_0_ptr[3] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 1, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[3] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 1, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[4] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 1, 0, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[4] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][{ 1, 0, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[5] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 1, 0, 1 }];
  v2e_cell_stencil_fused_face_1_ptr[5] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][{ 1, 0, 1 }];
  v2e_cell_stencil_fused_face_0_ptr[6] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 1, 1, -1 }];
  v2e_cell_stencil_fused_face_1_ptr[6] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][{ 1, 1, -1 }];
  v2e_cell_stencil_fused_face_0_ptr[7] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XY][{ 1, 1, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[7] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XY][{ 1, 1, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[8] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[8] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[9] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 1 }];
  v2e_cell_stencil_fused_face_1_ptr[9] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 1 }];
  v2e_cell_stencil_fused_face_0_ptr[10] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 1, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[10] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 1, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[11] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 1, -1, 1 }];
  v2e_cell_stencil_fused_face_1_ptr[11] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 1, -1, 1 }];
  v2e_cell_stencil_fused_face_0_ptr[12] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 1, 0, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[12] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 1, 0, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[13] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 1, 0, 1 }];
  v2e_cell_stencil_fused_face_1_ptr[13] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 1, 0, 1 }];
  v2e_cell_stencil_fused_face_0_ptr[14] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[14] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[15] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
  v2e_cell_stencil_fused_face_1_ptr[15] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
  v2e_cell_stencil_fused_face_0_ptr[16] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, -1 }];
  v2e_cell_stencil_fused_face_1_ptr[16] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, -1 }];
  v2e_cell_stencil_fused_face_0_ptr[17] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[17] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[18] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 1, -1, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[18] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][{ 1, -1, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[19] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 1, -1, 1 }];
  v2e_cell_stencil_fused_face_1_ptr[19] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][{ 1, -1, 1 }];
  v2e_cell_stencil_fused_face_0_ptr[20] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 1, 0, -1 }];
  v2e_cell_stencil_fused_face_1_ptr[20] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][{ 1, 0, -1 }];
  v2e_cell_stencil_fused_face_0_ptr[21] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::X][{ 1, 0, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[21] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::X][{ 1, 0, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[22] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 1, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[22] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 1, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[23] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 1, 1 }];
  v2e_cell_stencil_fused_face_1_ptr[23] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 1, 1 }];
  v2e_cell_stencil_fused_face_0_ptr[24] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[24] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[25] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 1 }];
  v2e_cell_stencil_fused_face_1_ptr[25] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 1 }];
  v2e_cell_stencil_fused_face_0_ptr[26] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 1, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[26] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 1, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[27] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 1, 1 }];
  v2e_cell_stencil_fused_face_1_ptr[27] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 1, 1 }];
  v2e_cell_stencil_fused_face_0_ptr[28] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[28] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[29] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, 1 }];
  v2e_cell_stencil_fused_face_1_ptr[29] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, 1 }];
  v2e_cell_stencil_fused_face_0_ptr[30] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 1, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[30] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 1, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[31] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[31] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[32] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 1 }];
  v2e_cell_stencil_fused_face_1_ptr[32] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 1 }];
  v2e_cell_stencil_fused_face_0_ptr[33] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 1, -1 }];
  v2e_cell_stencil_fused_face_1_ptr[33] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 1, -1 }];
  v2e_cell_stencil_fused_face_0_ptr[34] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 1, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[34] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 1, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[35] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, 0, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[35] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, 0, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[36] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 1 }];
  v2e_cell_stencil_fused_face_1_ptr[36] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 1 }];
  v2e_cell_stencil_fused_face_0_ptr[37] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[37] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[38] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 1 }];
  v2e_cell_stencil_fused_face_1_ptr[38] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 1 }];
  v2e_cell_stencil_fused_face_0_ptr[39] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[39] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[40] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 1 }];
  v2e_cell_stencil_fused_face_1_ptr[40] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 1 }];
  v2e_cell_stencil_fused_face_0_ptr[41] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[41] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, 0 }];
  v2e_cell_stencil_fused_face_0_ptr[42] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 1 }];
  v2e_cell_stencil_fused_face_1_ptr[42] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 1 }];
  v2e_cell_stencil_fused_face_0_ptr[43] = v2e_cell_stencil_fused_face_0[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
  v2e_cell_stencil_fused_face_1_ptr[43] = v2e_cell_stencil_fused_face_1[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_012(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_013(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_021(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_023(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_031(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_032(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_102(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_103(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_120(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_123(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_130(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_132(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_201(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_203(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_210(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_012_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_013(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_021(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_023(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_031(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_032(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_102(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_103(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_120(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_123(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_130(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_132(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_201(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_203(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_210(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_013_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_021(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_023(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_031(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_032(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_102(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_103(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_120(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_123(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_130(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_132(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_201(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_203(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_210(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_021_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_023(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_031(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_032(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_102(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_103(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_120(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_123(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_130(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_132(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_201(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_203(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_210(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_023_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_031(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_032(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_102(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_103(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_120(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_123(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_130(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_132(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_201(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_203(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_210(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_031_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_032(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_102(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_103(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_120(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_123(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_130(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_132(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_201(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_203(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_210(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_032_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_102(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_103(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_120(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_123(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_130(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_132(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_201(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_203(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_210(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_102_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_103_103(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_103_120(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_103_123(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_103_130(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_103_132(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_103_201(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_103_203(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_103_210(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_103_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_103_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_103_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_103_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_103_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_103_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_103_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_103_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_103_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_120_120(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_120_123(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_120_130(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_120_132(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_120_201(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_120_203(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_120_210(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_120_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_120_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_120_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_120_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_120_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_120_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_120_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_120_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_120_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_123_123(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_123_130(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_123_132(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_123_201(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_123_203(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_123_210(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_123_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_123_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_123_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_123_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_123_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_123_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_123_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_123_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_123_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_130_130(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_130_132(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_130_201(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_130_203(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_130_210(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_130_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_130_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_130_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_130_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_130_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_130_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_130_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_130_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_130_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_132_132(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_132_201(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_132_203(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_132_210(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_132_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_132_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_132_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_132_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_132_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_132_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_132_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_132_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_132_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_201_201(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_201_203(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_201_210(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_201_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_201_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_201_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_201_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_201_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_201_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_201_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_201_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_201_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_203_203(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_203_210(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_203_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_203_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_203_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_203_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_203_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_203_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_203_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_203_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_203_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_210_210(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_210_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_210_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_210_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_210_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_210_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_210_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_210_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_210_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_210_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_213_213(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_213_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_213_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_213_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_213_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_213_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_213_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_213_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_213_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_230_230(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_230_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_230_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_230_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_230_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_230_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_230_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_230_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_231_231(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_231_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_231_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_231_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_231_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_231_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_231_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_301_301(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_301_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_301_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_301_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_301_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_301_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_302_302(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_302_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_302_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_302_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_302_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_310_310(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_310_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_310_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_310_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_312_312(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_1_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_312_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_312_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((0) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_320_320(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_320_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_1_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_1_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_1_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P2_update_edgedofs_impl_321_321(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_edgeFaceRhs_X, _data_edgeFaceRhs_XY, _data_edgeFaceRhs_Y, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, e2e_cell_stencil_fused_face_0_ptr, e2e_cell_stencil_fused_face_1_ptr, level, relax, v2e_cell_stencil_fused_face_0_ptr, v2e_cell_stencil_fused_face_1_ptr);
      
      return;
   } 
}


} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hyteg
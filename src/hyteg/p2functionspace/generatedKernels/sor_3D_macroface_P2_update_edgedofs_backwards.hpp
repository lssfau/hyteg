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

#pragma once
#include "hyteg/edgedofspace/EdgeDoFOrientation.hpp"
#include "hyteg/indexing/Common.hpp"
#include <map>
#include "core/Macros.h"
#define RESTRICT WALBERLA_RESTRICT
#include "sor_3D_macroface_P2_update_edgedofs_backwards_impl.hpp"

namespace hyteg {
namespace P2 {
namespace macroface {
namespace generated {

void sor_3D_macroface_P2_update_edgedofs_backwards(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceDst_gl0_X, double * RESTRICT _data_edgeFaceDst_gl0_XY, double * RESTRICT _data_edgeFaceDst_gl0_XYZ, double * RESTRICT _data_edgeFaceDst_gl0_XZ, double * RESTRICT _data_edgeFaceDst_gl0_Y, double * RESTRICT _data_edgeFaceDst_gl0_YZ, double * RESTRICT _data_edgeFaceDst_gl0_Z, double * RESTRICT _data_edgeFaceDst_gl1_X, double * RESTRICT _data_edgeFaceDst_gl1_XY, double * RESTRICT _data_edgeFaceDst_gl1_XYZ, double * RESTRICT _data_edgeFaceDst_gl1_XZ, double * RESTRICT _data_edgeFaceDst_gl1_Y, double * RESTRICT _data_edgeFaceDst_gl1_YZ, double * RESTRICT _data_edgeFaceDst_gl1_Z, double const * RESTRICT const _data_edgeFaceRhs_X, double const * RESTRICT const _data_edgeFaceRhs_XY, double const * RESTRICT const _data_edgeFaceRhs_Y, double const * RESTRICT const _data_vertexFaceDst, double const * RESTRICT const _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceDst_gl1, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil_fused_face_0, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil_fused_face_1, int level, int64_t neighbor_cell_0_local_vertex_id_0, int64_t neighbor_cell_0_local_vertex_id_1, int64_t neighbor_cell_0_local_vertex_id_2, int64_t neighbor_cell_1_local_vertex_id_0, int64_t neighbor_cell_1_local_vertex_id_1, int64_t neighbor_cell_1_local_vertex_id_2, double relax, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > v2e_cell_stencil_fused_face_0, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > v2e_cell_stencil_fused_face_1);

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hyteg
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

namespace hyteg {
namespace EdgeDoFToVertexDoF {
namespace generated {

void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_012(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_013(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_021(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_023(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_031(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_032(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_102(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_103(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_120(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_123(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_130(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_132(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_201(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_203(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_210(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_213(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_230(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_231(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_301(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_302(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_310(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_312(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_320(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);
void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_321(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int level, int64_t num_neighbor_faces);

} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hyteg
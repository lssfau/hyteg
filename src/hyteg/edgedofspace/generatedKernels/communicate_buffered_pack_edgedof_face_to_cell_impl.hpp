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
#include "core/Macros.h"
#include "hyteg/edgedofspace/EdgeDoFOrientation.hpp"
#include "hyteg/indexing/Common.hpp"
#include <map>
#define RESTRICT WALBERLA_RESTRICT

namespace hyteg {
namespace edgedof {
namespace comm {
namespace generated {

void communicate_buffered_pack_edgedof_face_to_cell_impl_012(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_013(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_021(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_023(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_031(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_032(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_102(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_103(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_120(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_123(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_130(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_132(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_201(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_203(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_210(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_213(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_230(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_231(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_301(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_302(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_310(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_312(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_320(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);
void communicate_buffered_pack_edgedof_face_to_cell_impl_321(double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, double * RESTRICT _data_send_buffer, int level, int64_t send_buffer_first_element_idx);

} // namespace generated
} // namespace comm
} // namespace edgedof
} // namespace hyteg
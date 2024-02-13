/*
 * Copyright (c) 2019-2023 Nils Kohl, Dominik Thoennes, Michael Zikeli.
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
#include <map>
#include "core/Macros.h"
#include "hyteg/edgedofspace/EdgeDoFOrientation.hpp"
#include "hyteg/indexing/Common.hpp"
#define RESTRICT WALBERLA_RESTRICT

namespace hyteg {
namespace vertexdof {
namespace macroface {
namespace generated {

template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, int64_t neighbor_cell_local_vertex_id_0, int64_t neighbor_cell_local_vertex_id_1, int64_t neighbor_cell_local_vertex_id_2, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_012(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_013(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_021(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_023(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_031(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_032(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_102(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_103(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_120(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_123(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_130(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_132(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_201(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_203(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_210(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_213(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_230(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_231(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_301(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_302(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_310(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_312(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_320(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);
template < typename ValueType >
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add_321(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, ValueType >& p1FaceStencil);

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg
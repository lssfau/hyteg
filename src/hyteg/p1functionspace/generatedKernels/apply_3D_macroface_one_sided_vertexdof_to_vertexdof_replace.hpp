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
namespace vertexdof {
namespace macroface {
namespace generated {

void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, int64_t neighbor_cell_local_vertex_id_0, int64_t neighbor_cell_local_vertex_id_1, int64_t neighbor_cell_local_vertex_id_2, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_012(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_013(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_021(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_023(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_031(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_032(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_102(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_103(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_120(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_123(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_130(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_132(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_201(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_203(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_210(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_213(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_230(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_231(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_301(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_302(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_310(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_312(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_320(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);
void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_321(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int level, std::map< hyteg::indexing::Index, double > p1FaceStencil);

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg
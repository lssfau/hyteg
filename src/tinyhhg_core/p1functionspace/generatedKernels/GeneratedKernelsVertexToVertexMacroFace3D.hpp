
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include "core/DataTypes.h"

#include "core/Macros.h"

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"

#include <map>

#define RESTRICT WALBERLA_RESTRICT

namespace hhg {
namespace vertexdof {
namespace macroface {
namespace generated {

void apply_3D_macroface_vertexdof_to_vertexdof_replace(double * RESTRICT _data_p1CellDst, double const * RESTRICT const _data_p1CellSrc, double const * RESTRICT const _data_p1CellSrc_gl_cell_0, double const * RESTRICT const _data_p1CellSrc_gl_cell_1, int64_t level, int64_t neighbor_cell_0_local_vertex_id_0, int64_t neighbor_cell_0_local_vertex_id_1, int64_t neighbor_cell_0_local_vertex_id_2, int64_t neighbor_cell_1_local_vertex_id_0, int64_t neighbor_cell_1_local_vertex_id_1, int64_t neighbor_cell_1_local_vertex_id_2, std::map< walberla::uint_t, std::map< hhg::indexing::IndexIncrement, double > > p1CellStencil);

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hhg
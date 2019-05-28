
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once
#include "all.hpp"
#include "core/DataTypes.h"
#include "core/Macros.h"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include <map>
#define RESTRICT WALBERLA_RESTRICT

namespace hhg {
namespace edgedof {
namespace comm {
namespace generated {

void communicate_directly_edgedof_face_to_cell(double * RESTRICT _data_edge_cell_dst_X, double * RESTRICT _data_edge_cell_dst_XY, double * RESTRICT _data_edge_cell_dst_XZ, double * RESTRICT _data_edge_cell_dst_Y, double * RESTRICT _data_edge_cell_dst_YZ, double * RESTRICT _data_edge_cell_dst_Z, double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, int32_t level, int64_t neighbor_cell_local_vertex_id_0, int64_t neighbor_cell_local_vertex_id_1, int64_t neighbor_cell_local_vertex_id_2);

} // namespace generated
} // namespace comm
} // namespace edgedof
} // namespace hhg
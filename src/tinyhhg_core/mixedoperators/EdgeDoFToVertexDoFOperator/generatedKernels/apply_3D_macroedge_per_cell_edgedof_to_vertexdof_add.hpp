
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
namespace EdgeDoFToVertexDoF {
namespace generated {

void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int32_t level, int64_t neighbor_cell_local_vertex_id_0, int64_t neighbor_cell_local_vertex_id_1, int64_t neighbor_cell_local_vertex_id_2, int64_t num_neighbor_faces);

} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hhg

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once
#include "all.hpp"
#include "core/Macros.h"
#include "hyteg/edgedofspace/EdgeDoFOrientation.hpp"
#include "hyteg/indexing/Common.hpp"
#include <map>
#define RESTRICT WALBERLA_RESTRICT

namespace hyteg {
namespace EdgeDoFToVertexDoF {
namespace generated {

void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_301(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::IndexIncrement, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int32_t level, int64_t num_neighbor_faces);

} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hyteg
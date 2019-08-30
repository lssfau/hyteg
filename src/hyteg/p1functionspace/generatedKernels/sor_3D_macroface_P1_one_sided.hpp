
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
namespace vertexdof {
namespace macroface {
namespace generated {

void sor_3D_macroface_P1_one_sided(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int32_t level, int64_t neighbor_cell_0_local_vertex_id_0, int64_t neighbor_cell_0_local_vertex_id_1, int64_t neighbor_cell_0_local_vertex_id_2, double relax, std::map< hyteg::indexing::IndexIncrement, double > v2v_cell_stencil_fused_face_0);

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg
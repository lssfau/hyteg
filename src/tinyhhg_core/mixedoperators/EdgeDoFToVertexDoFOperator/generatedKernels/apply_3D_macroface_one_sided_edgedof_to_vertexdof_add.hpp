
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once
#include "all.hpp"
#include "core/Macros.h"
#include "tinyhhg_core/edgedofspace/EdgeDoFOrientation.hpp"
#include "tinyhhg_core/indexing/Common.hpp"
#include <map>
#define RESTRICT WALBERLA_RESTRICT

namespace hyteg {
namespace EdgeDoFToVertexDoF {
namespace generated {

void apply_3D_macroface_one_sided_edgedof_to_vertexdof_add(double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * RESTRICT const _data_edgeFaceSrc_gl0_X, double const * RESTRICT const _data_edgeFaceSrc_gl0_XY, double const * RESTRICT const _data_edgeFaceSrc_gl0_XYZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_XZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_Y, double const * RESTRICT const _data_edgeFaceSrc_gl0_YZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_Z, double * RESTRICT _data_vertexFaceDst, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::IndexIncrement, double > > e2v_cell_stencil, int32_t level, int64_t neighbor_cell_local_vertex_id_0, int64_t neighbor_cell_local_vertex_id_1, int64_t neighbor_cell_local_vertex_id_2);

} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hyteg
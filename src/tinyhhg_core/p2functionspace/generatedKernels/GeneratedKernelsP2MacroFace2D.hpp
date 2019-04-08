
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include "core/DataTypes.h"

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"

#include <map>

namespace hhg {
namespace P2 {
namespace macroface {
namespace generated {

void sor_2D_macroface_P2_update_vertexdofs(double * _data_edgeFaceDst, double * const _data_edge_stencil_at_vertex, double * _data_vertexFaceDst, double * _data_vertexFaceRhs, double * const _data_vertex_stencil_at_vertex, int64_t level, double relax);

void sor_2D_macroface_P2_update_edgedofs(double * _data_edgeFaceDst, double * _data_edgeFaceRhs, double * const _data_edge_stencil_at_edge_x, double * const _data_edge_stencil_at_edge_xy, double * const _data_edge_stencil_at_edge_y, double * _data_vertexFaceDst, double * const _data_vertex_stencil_at_edge_x, double * const _data_vertex_stencil_at_edge_xy, double * const _data_vertex_stencil_at_edge_y, int64_t level, double relax);

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hhg
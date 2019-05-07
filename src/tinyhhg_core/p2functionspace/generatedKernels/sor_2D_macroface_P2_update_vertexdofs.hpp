
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
namespace P2 {
namespace macroface {
namespace generated {

void sor_2D_macroface_P2_update_vertexdofs(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edge_stencil_at_vertex, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceRhs, double const * RESTRICT const _data_vertex_stencil_at_vertex, int64_t level, double relax);

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hhg
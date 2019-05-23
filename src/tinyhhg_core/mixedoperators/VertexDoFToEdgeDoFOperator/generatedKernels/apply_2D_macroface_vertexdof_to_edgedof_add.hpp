
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
namespace VertexDoFToEdgeDoF {
namespace generated {

void apply_2D_macroface_vertexdof_to_edgedof_add(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_vertexFaceSrc, double const * RESTRICT const _data_vertexToDiagonalEdgeFaceStencil, double const * RESTRICT const _data_vertexToHorizontalEdgeFaceStencil, double const * RESTRICT const _data_vertexToVerticalEdgeFaceStencil, int32_t level);

} // namespace generated
} // namespace VertexDoFToEdgeDoF
} // namespace hhg
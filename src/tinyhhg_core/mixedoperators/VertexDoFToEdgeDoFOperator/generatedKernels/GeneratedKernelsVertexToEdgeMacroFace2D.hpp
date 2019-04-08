
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include "core/DataTypes.h"

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"

#include <map>

namespace hhg {
namespace VertexDoFToEdgeDoF {
namespace generated {

void apply_2D_macroface_vertexdof_to_edgedof_replace(double * _data_edgeFaceDst, double * _data_vertexFaceSrc, double * const _data_vertexToDiagonalEdgeFaceStencil, double * const _data_vertexToHorizontalEdgeFaceStencil, double * const _data_vertexToVerticalEdgeFaceStencil, int64_t level);

void apply_2D_macroface_vertexdof_to_edgedof_add(double * _data_edgeFaceDst, double * _data_vertexFaceSrc, double * const _data_vertexToDiagonalEdgeFaceStencil, double * const _data_vertexToHorizontalEdgeFaceStencil, double * const _data_vertexToVerticalEdgeFaceStencil, int64_t level);

} // namespace generated
} // namespace VertexDoFToEdgeDoF
} // namespace hhg

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include "core/DataTypes.h"

namespace hhg {
namespace VertexDoFToEdgeDoF {
namespace generated {

void apply_2D_macroface_vertexdof_to_edgedof_replace(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil, int64_t level);

void apply_2D_macroface_vertexdof_to_edgedof_add(double * fd_edgeFaceDst, double * fd_vertexFaceSrc, double * fd_vertexToDiagonalEdgeFaceStencil, double * fd_vertexToHorizontalEdgeFaceStencil, double * fd_vertexToVerticalEdgeFaceStencil, int64_t level);

} // namespace generated
} // namespace VertexDoFToEdgeDoF
} // namespace hhg
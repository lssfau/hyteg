
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include "core/DataTypes.h"

#include "core/Macros.h"

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"

#include <map>

#define RESTRICT WALBERLA_RESTRICT

namespace hhg {
namespace EdgeDoFToVertexDoF {
namespace generated {

void apply_2D_macroface_edgedof_to_vertexdof_replace(double const * RESTRICT const _data_edgeFaceSrc, double const * const _data_edgeToVertexFaceStencil, double * RESTRICT _data_p1FaceDst, int64_t level);

void apply_2D_macroface_edgedof_to_vertexdof_add(double const * RESTRICT const _data_edgeFaceSrc, double const * const _data_edgeToVertexFaceStencil, double * RESTRICT _data_p1FaceDst, int64_t level);

} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hhg
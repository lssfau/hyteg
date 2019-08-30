
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
namespace VertexDoFToEdgeDoF {
namespace generated {

void apply_3D_macroface_one_sided_vertexdof_to_edgedof_replace_031(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_vertexFaceSrc, double const * RESTRICT const _data_vertexFaceSrc_gl0, int32_t level, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::IndexIncrement, double > > v2e_cell_stencil);

} // namespace generated
} // namespace VertexDoFToEdgeDoF
} // namespace hyteg
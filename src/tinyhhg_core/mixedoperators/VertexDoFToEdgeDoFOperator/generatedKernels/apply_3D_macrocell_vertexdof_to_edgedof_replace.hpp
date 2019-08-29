
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

namespace hhg {
namespace VertexDoFToEdgeDoF {
namespace generated {

void apply_3D_macrocell_vertexdof_to_edgedof_replace(double * RESTRICT _data_edgeCellDst_X, double * RESTRICT _data_edgeCellDst_XY, double * RESTRICT _data_edgeCellDst_XYZ, double * RESTRICT _data_edgeCellDst_XZ, double * RESTRICT _data_edgeCellDst_Y, double * RESTRICT _data_edgeCellDst_YZ, double * RESTRICT _data_edgeCellDst_Z, double const * RESTRICT const _data_vertexCellSrc, int32_t level, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > v2eStencilMap);

} // namespace generated
} // namespace VertexDoFToEdgeDoF
} // namespace hhg
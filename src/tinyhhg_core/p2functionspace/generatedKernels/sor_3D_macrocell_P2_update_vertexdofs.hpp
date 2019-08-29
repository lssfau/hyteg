
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
namespace P2 {
namespace macrocell {
namespace generated {

void sor_3D_macrocell_P2_update_vertexdofs(double const * RESTRICT const _data_edgeCellDst_X, double const * RESTRICT const _data_edgeCellDst_XY, double const * RESTRICT const _data_edgeCellDst_XYZ, double const * RESTRICT const _data_edgeCellDst_XZ, double const * RESTRICT const _data_edgeCellDst_Y, double const * RESTRICT const _data_edgeCellDst_YZ, double const * RESTRICT const _data_edgeCellDst_Z, double * RESTRICT _data_vertexCellDst, double const * RESTRICT const _data_vertexCellRhs, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2vStencilMap, int32_t level, double relax, std::map< hhg::indexing::IndexIncrement, double > v2vStencilMap);

} // namespace generated
} // namespace macrocell
} // namespace P2
} // namespace hhg
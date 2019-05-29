
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
namespace macrocell {
namespace generated {

void sor_3D_macrocell_P2_update_edgedofs_by_type_XZ(double const * RESTRICT const _data_edgeCellDst_X, double const * RESTRICT const _data_edgeCellDst_XY, double const * RESTRICT const _data_edgeCellDst_XYZ, double * RESTRICT _data_edgeCellDst_XZ, double const * RESTRICT const _data_edgeCellDst_Y, double const * RESTRICT const _data_edgeCellDst_YZ, double const * RESTRICT const _data_edgeCellDst_Z, double const * RESTRICT const _data_edgeCellRhs_XZ, double const * RESTRICT const _data_vertexCellDst, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > > e2eStencilMap, int32_t level, double relax, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > v2eStencilMap);

} // namespace generated
} // namespace macrocell
} // namespace P2
} // namespace hhg
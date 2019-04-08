
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

void apply_3D_macrocell_vertexdof_to_edgedof_replace(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_vertexCellSrc, int64_t level, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > v2eStencilMap);

void apply_3D_macrocell_vertexdof_to_edgedof_add(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_vertexCellSrc, int64_t level, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > v2eStencilMap);

} // namespace generated
} // namespace VertexDoFToEdgeDoF
} // namespace hhg
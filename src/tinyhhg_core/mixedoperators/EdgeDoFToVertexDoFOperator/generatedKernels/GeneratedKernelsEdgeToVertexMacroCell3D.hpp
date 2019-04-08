
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include "core/DataTypes.h"

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"

#include <map>

namespace hhg {
namespace EdgeDoFToVertexDoF {
namespace generated {

void apply_3D_macrocell_edgedof_to_vertexdof_replace(double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double * _data_vertexCellDst, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2vStencilMap, int64_t level);

void apply_3D_macrocell_edgedof_to_vertexdof_add(double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double * _data_vertexCellDst, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2vStencilMap, int64_t level);

} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hhg
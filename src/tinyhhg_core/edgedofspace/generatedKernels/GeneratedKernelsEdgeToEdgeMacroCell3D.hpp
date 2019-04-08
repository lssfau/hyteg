
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include "core/DataTypes.h"

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"

#include <map>

namespace hhg {
namespace edgedof {
namespace macrocell {
namespace generated {

void assign_3D_macrocell_edgedof_1_rhs_function(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, double c, int64_t level);

void apply_3D_macrocell_edgedof_to_edgedof_replace(double * _data_edgeCellDst_X, double * _data_edgeCellDst_XY, double * _data_edgeCellDst_XYZ, double * _data_edgeCellDst_XZ, double * _data_edgeCellDst_Y, double * _data_edgeCellDst_YZ, double * _data_edgeCellDst_Z, double * _data_edgeCellSrc_X, double * _data_edgeCellSrc_XY, double * _data_edgeCellSrc_XYZ, double * _data_edgeCellSrc_XZ, double * _data_edgeCellSrc_Y, double * _data_edgeCellSrc_YZ, double * _data_edgeCellSrc_Z, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > > e2eStencilMap, int64_t level);

} // namespace generated
} // namespace macrocell
} // namespace edgedof
} // namespace hhg

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
namespace edgedof {
namespace macroface {
namespace generated {

void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_301(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * RESTRICT const _data_edgeFaceSrc_gl0_X, double const * RESTRICT const _data_edgeFaceSrc_gl0_XY, double const * RESTRICT const _data_edgeFaceSrc_gl0_XYZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_XZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_YZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_Z, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > > e2e_cell_stencil, int64_t level);

} // namespace generated
} // namespace macroface
} // namespace edgedof
} // namespace hhg
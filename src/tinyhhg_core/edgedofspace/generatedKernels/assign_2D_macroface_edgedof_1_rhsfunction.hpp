
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once
#include "all.hpp"
#include "core/Macros.h"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include <map>
#define RESTRICT WALBERLA_RESTRICT

namespace hhg {
namespace edgedof {
namespace macroface {
namespace generated {

void assign_2D_macroface_edgedof_1_rhs_function(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceSrc_X, double * RESTRICT _data_edgeFaceSrc_XY, double * RESTRICT _data_edgeFaceSrc_Y, double c, int32_t level);

} // namespace generated
} // namespace macroface
} // namespace edgedof
} // namespace hhg
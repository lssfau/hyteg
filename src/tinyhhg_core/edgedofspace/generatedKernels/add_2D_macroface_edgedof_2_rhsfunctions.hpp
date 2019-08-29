
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
namespace edgedof {
namespace macroface {
namespace generated {

void add_2D_macroface_edgedof_2_rhs_functions(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double * RESTRICT _data_edgeFaceSrc0_X, double * RESTRICT _data_edgeFaceSrc0_XY, double * RESTRICT _data_edgeFaceSrc0_Y, double * RESTRICT _data_edgeFaceSrc1_X, double * RESTRICT _data_edgeFaceSrc1_XY, double * RESTRICT _data_edgeFaceSrc1_Y, double c0, double c1, int32_t level);

} // namespace generated
} // namespace macroface
} // namespace edgedof
} // namespace hhg
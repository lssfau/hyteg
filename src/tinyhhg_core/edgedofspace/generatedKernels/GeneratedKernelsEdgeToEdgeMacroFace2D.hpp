
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include "core/DataTypes.h"

#include "core/Macros.h"

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"

#include <map>

#define RESTRICT WALBERLA_RESTRICT

namespace hhg {
namespace edgedof {
namespace macroface {
namespace generated {

void apply_2D_macroface_edgedof_to_edgedof_replace(double * RESTRICT _data_edgeFaceDst, double const * RESTRICT const _data_edgeFaceSrc, double const * const _data_edgeToDiagonalEdgeFaceStencil, double const * const _data_edgeToHorizontalEdgeFaceStencil, double const * const _data_edgeToVerticalEdgeFaceStencil, int64_t level);

void apply_2D_macroface_edgedof_to_edgedof_add(double * RESTRICT _data_edgeFaceDst, double const * RESTRICT const _data_edgeFaceSrc, double const * const _data_edgeToDiagonalEdgeFaceStencil, double const * const _data_edgeToHorizontalEdgeFaceStencil, double const * const _data_edgeToVerticalEdgeFaceStencil, int64_t level);

void assign_2D_macroface_edgedof_1_rhs_function(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceSrc, double c, int64_t level);

void assign_2D_macroface_edgedof_2_rhs_functions(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceSrc0, double * RESTRICT _data_edgeFaceSrc1, double c0, double c1, int64_t level);

void assign_2D_macroface_edgedof_3_rhs_functions(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceSrc0, double * RESTRICT _data_edgeFaceSrc1, double * RESTRICT _data_edgeFaceSrc2, double c0, double c1, double c2, int64_t level);

} // namespace generated
} // namespace macroface
} // namespace edgedof
} // namespace hhg
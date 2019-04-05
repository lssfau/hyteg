
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include "core/DataTypes.h"

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"

#include <map>

namespace hhg {
namespace edgedof {
namespace macroface {
namespace generated {

void apply_2D_macroface_edgedof_to_edgedof_replace(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double * const _data_edgeToDiagonalEdgeFaceStencil, double * const _data_edgeToHorizontalEdgeFaceStencil, double * const _data_edgeToVerticalEdgeFaceStencil, int64_t level);

void apply_2D_macroface_edgedof_to_edgedof_add(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double * const _data_edgeToDiagonalEdgeFaceStencil, double * const _data_edgeToHorizontalEdgeFaceStencil, double * const _data_edgeToVerticalEdgeFaceStencil, int64_t level);

void assign_2D_macroface_edgedof_1_rhs_function(double * _data_edgeFaceDst, double * _data_edgeFaceSrc, double c, int64_t level);

void assign_2D_macroface_edgedof_2_rhs_functions(double * _data_edgeFaceDst, double * _data_edgeFaceSrc0, double * _data_edgeFaceSrc1, double c0, double c1, int64_t level);

void assign_2D_macroface_edgedof_3_rhs_functions(double * _data_edgeFaceDst, double * _data_edgeFaceSrc0, double * _data_edgeFaceSrc1, double * _data_edgeFaceSrc2, double c0, double c1, double c2, int64_t level);

} // namespace generated
} // namespace macroface
} // namespace edgedof
} // namespace hhg
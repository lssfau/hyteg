
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
namespace vertexdof {
namespace macroface {
namespace generated {

void apply_2D_macroface_vertexdof_to_vertexdof_replace(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * const _data_p1FaceStencil, int64_t level);

void apply_2D_macroface_vertexdof_to_vertexdof_add(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * const _data_p1FaceStencil, int64_t level);

void gaussseidel_2D_macroface_vertexdof_to_vertexdof(double * RESTRICT _data_p1FaceDst, double * RESTRICT _data_p1FaceRhs, double const * const _data_p1FaceStencil, int64_t level);

void assign_2D_macroface_vertexdof_1_rhs_function(double * RESTRICT _data_p1FaceDst, double * RESTRICT _data_p1FaceSrc, double c, int64_t level);

void assign_2D_macroface_vertexdof_2_rhs_functions(double * RESTRICT _data_p1FaceDst, double * RESTRICT _data_p1FaceSrc0, double * RESTRICT _data_p1FaceSrc1, double c0, double c1, int64_t level);

void assign_2D_macroface_vertexdof_3_rhs_functions(double * RESTRICT _data_p1FaceDst, double * RESTRICT _data_p1FaceSrc0, double * RESTRICT _data_p1FaceSrc1, double * RESTRICT _data_p1FaceSrc3, double c0, double c1, double c2, int64_t level);

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hhg
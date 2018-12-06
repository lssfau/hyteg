
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include "core/DataTypes.h"

namespace hhg {
namespace vertexdof {
namespace macroface {
namespace generated {

void apply_2D_macroface_vertexdof_to_vertexdof_replace(double * _data_p1FaceDst, double * _data_p1FaceSrc, double * const _data_p1FaceStencil, int64_t level);

void apply_2D_macroface_vertexdof_to_vertexdof_add(double * _data_p1FaceDst, double * _data_p1FaceSrc, double * const _data_p1FaceStencil, int64_t level);

void gaussseidel_2D_macroface_vertexdof_to_vertexdof(double * _data_p1FaceDst, double * _data_p1FaceRhs, double * const _data_p1FaceStencil, int64_t level);

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hhg
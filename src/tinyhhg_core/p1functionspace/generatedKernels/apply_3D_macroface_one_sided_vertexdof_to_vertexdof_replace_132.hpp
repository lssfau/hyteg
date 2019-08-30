
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

namespace hyteg {
namespace vertexdof {
namespace macroface {
namespace generated {

void apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace_132(double * RESTRICT _data_p1FaceDst, double const * RESTRICT const _data_p1FaceSrc, double const * RESTRICT const _data_p1FaceSrc_gl0, int32_t level, std::map< hyteg::indexing::IndexIncrement, double > p1FaceStencil);

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg

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
namespace vertexdof {
namespace macroface {
namespace generated {

void sor_3D_macroface_P1_one_sided_impl_312(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int32_t level, double relax, std::map< hhg::indexing::IndexIncrement, double > v2v_cell_stencil_fused_face_0);

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hhg
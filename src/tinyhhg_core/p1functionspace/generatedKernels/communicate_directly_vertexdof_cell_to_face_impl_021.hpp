
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
namespace vertexdof {
namespace comm {
namespace generated {

void communicate_directly_vertexdof_cell_to_face_impl_021(double const * RESTRICT const _data_p1_cell_src, double * RESTRICT _data_p1_face_dst_gl0, int32_t level);

} // namespace generated
} // namespace comm
} // namespace vertexdof
} // namespace hhg
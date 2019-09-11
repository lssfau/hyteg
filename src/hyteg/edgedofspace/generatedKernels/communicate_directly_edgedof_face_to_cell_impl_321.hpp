
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once
#include "all.hpp"
#include "core/Macros.h"
#include "hyteg/edgedofspace/EdgeDoFOrientation.hpp"
#include "hyteg/indexing/Common.hpp"
#include <map>
#define RESTRICT WALBERLA_RESTRICT

namespace hyteg {
namespace edgedof {
namespace comm {
namespace generated {

void communicate_directly_edgedof_face_to_cell_impl_321(double * RESTRICT _data_edge_cell_dst_XY, double * RESTRICT _data_edge_cell_dst_XZ, double * RESTRICT _data_edge_cell_dst_YZ, double const * RESTRICT const _data_edge_face_src_X, double const * RESTRICT const _data_edge_face_src_XY, double const * RESTRICT const _data_edge_face_src_Y, int32_t level);

} // namespace generated
} // namespace comm
} // namespace edgedof
} // namespace hyteg
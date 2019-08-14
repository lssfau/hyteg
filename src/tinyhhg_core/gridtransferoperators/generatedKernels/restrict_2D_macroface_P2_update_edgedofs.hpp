
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
namespace P2 {
namespace macroface {
namespace generated {

void restrict_2D_macroface_P2_update_edgedofs(double * RESTRICT _data_edgeCoarseDst_X, double * RESTRICT _data_edgeCoarseDst_XY, double * RESTRICT _data_edgeCoarseDst_Y, double * RESTRICT _data_edgeFineSrc_X, double * RESTRICT _data_edgeFineSrc_XY, double * RESTRICT _data_edgeFineSrc_Y, double * RESTRICT _data_vertexFineSrc, int32_t coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2);

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hhg
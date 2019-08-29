
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
namespace P2 {
namespace macroface {
namespace generated {

void prolongate_2D_macroface_P2_push_from_edgedofs(double * RESTRICT _data_edgeCoarseSrc_X, double * RESTRICT _data_edgeCoarseSrc_XY, double * RESTRICT _data_edgeCoarseSrc_Y, double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexFineDst, int32_t coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2);

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hhg
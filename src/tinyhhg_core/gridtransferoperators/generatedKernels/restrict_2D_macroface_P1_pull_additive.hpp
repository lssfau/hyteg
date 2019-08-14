
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
namespace macroface {
namespace generated {

void restrict_2D_macroface_P1_pull_additive(double * RESTRICT _data_vertexCoarseDst, double const * RESTRICT const _data_vertexFineSrc, int32_t coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2);

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hhg
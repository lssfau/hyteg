
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
namespace P2 {
namespace macroface {
namespace generated {

void restrict_2D_macroface_P2_update_vertexdofs(double * RESTRICT _data_edgeFineSrc_X, double * RESTRICT _data_edgeFineSrc_XY, double * RESTRICT _data_edgeFineSrc_Y, double * RESTRICT _data_vertexCoarseDst, double * RESTRICT _data_vertexFineSrc, int32_t coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2);

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hyteg
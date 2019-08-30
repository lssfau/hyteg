
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
namespace vertexdof {
namespace macrocell {
namespace generated {

void restrict_3D_macrocell_P1_pull_additive(double * RESTRICT _data_vertexCoarseDst, double const * RESTRICT const _data_vertexFineSrc, int32_t coarse_level, double num_neighbor_cells_edge0, double num_neighbor_cells_edge1, double num_neighbor_cells_edge2, double num_neighbor_cells_edge3, double num_neighbor_cells_edge4, double num_neighbor_cells_edge5, double num_neighbor_cells_face0, double num_neighbor_cells_face1, double num_neighbor_cells_face2, double num_neighbor_cells_face3, double num_neighbor_cells_vertex0, double num_neighbor_cells_vertex1, double num_neighbor_cells_vertex2, double num_neighbor_cells_vertex3);

} // namespace generated
} // namespace macrocell
} // namespace vertexdof
} // namespace hyteg
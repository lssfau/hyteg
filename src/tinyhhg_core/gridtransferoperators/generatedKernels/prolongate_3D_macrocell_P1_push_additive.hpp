
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
namespace macrocell {
namespace generated {

void prolongate_3D_macrocell_P1_push_additive(double const * RESTRICT const _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, int32_t coarse_level, double num_neighbor_cells_edge0, double num_neighbor_cells_edge1, double num_neighbor_cells_edge2, double num_neighbor_cells_edge3, double num_neighbor_cells_edge4, double num_neighbor_cells_edge5, double num_neighbor_cells_face0, double num_neighbor_cells_face1, double num_neighbor_cells_face2, double num_neighbor_cells_face3, double num_neighbor_cells_vertex0, double num_neighbor_cells_vertex1, double num_neighbor_cells_vertex2, double num_neighbor_cells_vertex3);

} // namespace generated
} // namespace macrocell
} // namespace vertexdof
} // namespace hhg

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
namespace macrocell {
namespace generated {

void restrict_3D_macrocell_P2_update_vertexdofs(double const * RESTRICT const _data_edgeFineSrc_X, double const * RESTRICT const _data_edgeFineSrc_XY, double const * RESTRICT const _data_edgeFineSrc_XYZ, double const * RESTRICT const _data_edgeFineSrc_XZ, double const * RESTRICT const _data_edgeFineSrc_Y, double const * RESTRICT const _data_edgeFineSrc_YZ, double const * RESTRICT const _data_edgeFineSrc_Z, double * RESTRICT _data_vertexCoarseDst, double const * RESTRICT const _data_vertexFineSrc, int32_t coarse_level, double num_neighbor_cells_edge0, double num_neighbor_cells_edge1, double num_neighbor_cells_edge2, double num_neighbor_cells_edge3, double num_neighbor_cells_edge4, double num_neighbor_cells_edge5, double num_neighbor_cells_face0, double num_neighbor_cells_face1, double num_neighbor_cells_face2, double num_neighbor_cells_face3, double num_neighbor_cells_vertex0, double num_neighbor_cells_vertex1, double num_neighbor_cells_vertex2, double num_neighbor_cells_vertex3);

} // namespace generated
} // namespace macrocell
} // namespace P2
} // namespace hhg
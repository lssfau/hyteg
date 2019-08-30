
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

namespace hyteg {
namespace P2 {
namespace macrocell {
namespace generated {

void restrict_3D_macrocell_P2_update_edgedofs(double * RESTRICT _data_edgeCoarseDst_X, double * RESTRICT _data_edgeCoarseDst_XY, double * RESTRICT _data_edgeCoarseDst_XYZ, double * RESTRICT _data_edgeCoarseDst_XZ, double * RESTRICT _data_edgeCoarseDst_Y, double * RESTRICT _data_edgeCoarseDst_YZ, double * RESTRICT _data_edgeCoarseDst_Z, double const * RESTRICT const _data_edgeFineSrc_X, double const * RESTRICT const _data_edgeFineSrc_XY, double const * RESTRICT const _data_edgeFineSrc_XYZ, double const * RESTRICT const _data_edgeFineSrc_XZ, double const * RESTRICT const _data_edgeFineSrc_Y, double const * RESTRICT const _data_edgeFineSrc_YZ, double const * RESTRICT const _data_edgeFineSrc_Z, double const * RESTRICT const _data_vertexFineSrc, int32_t coarse_level, double num_neighbor_cells_edge0, double num_neighbor_cells_edge1, double num_neighbor_cells_edge2, double num_neighbor_cells_edge3, double num_neighbor_cells_edge4, double num_neighbor_cells_edge5, double num_neighbor_cells_face0, double num_neighbor_cells_face1, double num_neighbor_cells_face2, double num_neighbor_cells_face3);

} // namespace generated
} // namespace macrocell
} // namespace P2
} // namespace hyteg
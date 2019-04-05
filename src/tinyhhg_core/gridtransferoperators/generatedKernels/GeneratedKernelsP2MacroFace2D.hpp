
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include "core/DataTypes.h"

#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"

#include <map>

namespace hhg {
namespace P2 {
namespace macroface {
namespace generated {

void restrict_2D_macroface_P2_update_vertexdofs(double * _data_edgeFineSrc, double * _data_vertexCoarseDst, double * _data_vertexFineSrc, int64_t coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2);

void restrict_2D_macroface_P2_update_edgedofs(double * _data_edgeCoarseDst, double * _data_edgeFineSrc, double * _data_vertexFineSrc, int64_t coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2);

void prolongate_2D_macroface_P2_push_from_vertexdofs(double * _data_edgeFineDst, double * _data_vertexCoarseSrc, double * _data_vertexFineDst, int64_t coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2);

void prolongate_2D_macroface_P2_push_from_edgedofs(double * _data_edgeCoarseSrc, double * _data_edgeFineDst, double * _data_vertexFineDst, int64_t coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2);

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hhg
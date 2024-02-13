/*
 * Copyright (c) 2019 Nils Kohl, Dominik Thoennes.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "sor_3D_macroedge_P2_update_vertexdofs_impl.hpp"

namespace hyteg {
namespace P2 {
namespace macroedge {
namespace generated {

static void sor_3D_macroedge_P2_update_vertexdofs_impl_012_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_16 = v2v_cell_stencil[{ -1, 0, 0 }];
   const double xi_17 = v2v_cell_stencil[{ -1, 0, 1 }];
   const double xi_18 = v2v_cell_stencil[{ -1, 1, 0 }];
   const double xi_19 = v2v_cell_stencil[{ 0, 0, 1 }];
   const double xi_20 = v2v_cell_stencil[{ 0, 1, 0 }];
   const double xi_21 = v2v_cell_stencil[{ 1, 0, 0 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_11*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_12*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_13*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_14*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_15*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_16*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_17*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_18*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_19*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_2*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_20*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_21*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_3*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_4*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_5*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_6*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_7*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_8*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_9*_data_edgedof_macroedge_src[micro_edge_index_x];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_012(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_012_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_013_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_16 = v2v_cell_stencil[{ -1, 0, 0 }];
   const double xi_17 = v2v_cell_stencil[{ -1, 0, 1 }];
   const double xi_18 = v2v_cell_stencil[{ -1, 1, 0 }];
   const double xi_19 = v2v_cell_stencil[{ 0, 0, 1 }];
   const double xi_20 = v2v_cell_stencil[{ 0, 1, 0 }];
   const double xi_21 = v2v_cell_stencil[{ 1, 0, 0 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_11*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_12*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_13*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_14*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_15*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_16*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_17*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_18*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_19*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_2*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_20*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_21*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_3*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_4*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_5*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_6*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_7*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_8*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_9*_data_edgedof_macroedge_src[micro_edge_index_x];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_013(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_013_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_021_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_16 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_17 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_18 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_19 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_20 = v2v_cell_stencil[{ 0, -1, 0 }];
   const double xi_21 = v2v_cell_stencil[{ 0, -1, 1 }];
   const double xi_22 = v2v_cell_stencil[{ 0, 0, 1 }];
   const double xi_23 = v2v_cell_stencil[{ 0, 1, 0 }];
   const double xi_24 = v2v_cell_stencil[{ 1, -1, 0 }];
   const double xi_25 = v2v_cell_stencil[{ 1, -1, 1 }];
   const double xi_26 = v2v_cell_stencil[{ 1, 0, 0 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_11*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_12*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 5*(1 << (level)) - 6] + xi_13*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_14*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_15*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_16*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_17*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_18*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_19*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 3*(1 << (level)) - 4] + xi_2*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_20*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_21*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_22*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_23*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_24*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_25*_data_vertexdof_macroedge_src[id_cell*((1 << (level)) - 1) + micro_edge_index_x + num_neighbor_faces*(1 << (level)) + (1 << (level))] + xi_26*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_3*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 4*(1 << (level)) - 5] + xi_4*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_5*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_6*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_7*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_8*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 2*(1 << (level)) - 3] + xi_9*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_021(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_021_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_023_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_16 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_17 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_18 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_19 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_20 = v2v_cell_stencil[{ 0, -1, 0 }];
   const double xi_21 = v2v_cell_stencil[{ 0, -1, 1 }];
   const double xi_22 = v2v_cell_stencil[{ 0, 0, 1 }];
   const double xi_23 = v2v_cell_stencil[{ 0, 1, 0 }];
   const double xi_24 = v2v_cell_stencil[{ 1, -1, 0 }];
   const double xi_25 = v2v_cell_stencil[{ 1, -1, 1 }];
   const double xi_26 = v2v_cell_stencil[{ 1, 0, 0 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_11*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_12*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 4*(1 << (level)) - 5] + xi_13*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_14*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_15*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_16*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_17*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_18*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_19*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 2*(1 << (level)) - 3] + xi_2*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_20*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_21*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_22*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_23*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_24*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_25*_data_vertexdof_macroedge_src[id_cell*((1 << (level)) - 1) + micro_edge_index_x + num_neighbor_faces*(1 << (level)) + (1 << (level))] + xi_26*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_3*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 5*(1 << (level)) - 6] + xi_4*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_5*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_6*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_7*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_8*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 3*(1 << (level)) - 4] + xi_9*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_023(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_023_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_031_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_16 = v2v_cell_stencil[{ 0, 0, -1 }];
   const double xi_17 = v2v_cell_stencil[{ 0, 0, 1 }];
   const double xi_18 = v2v_cell_stencil[{ 0, 1, -1 }];
   const double xi_19 = v2v_cell_stencil[{ 0, 1, 0 }];
   const double xi_20 = v2v_cell_stencil[{ 1, 0, -1 }];
   const double xi_21 = v2v_cell_stencil[{ 1, 0, 0 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_11*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_12*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_13*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_14*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_15*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_16*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_17*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_18*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_19*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_2*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_20*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_21*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_3*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_4*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_5*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_6*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_7*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_8*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_9*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_031(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_031_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_032_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_16 = v2v_cell_stencil[{ 0, 0, -1 }];
   const double xi_17 = v2v_cell_stencil[{ 0, 0, 1 }];
   const double xi_18 = v2v_cell_stencil[{ 0, 1, -1 }];
   const double xi_19 = v2v_cell_stencil[{ 0, 1, 0 }];
   const double xi_20 = v2v_cell_stencil[{ 1, 0, -1 }];
   const double xi_21 = v2v_cell_stencil[{ 1, 0, 0 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_11*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_12*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_13*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_14*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_15*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_16*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_17*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_18*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_19*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_2*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_20*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_21*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_3*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_4*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_5*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_6*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_7*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_8*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_9*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_032(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_032_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_102_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_16 = v2v_cell_stencil[{ -1, 0, 0 }];
   const double xi_17 = v2v_cell_stencil[{ -1, 0, 1 }];
   const double xi_18 = v2v_cell_stencil[{ -1, 1, 0 }];
   const double xi_19 = v2v_cell_stencil[{ 0, 0, 1 }];
   const double xi_20 = v2v_cell_stencil[{ 0, 1, 0 }];
   const double xi_21 = v2v_cell_stencil[{ 1, 0, 0 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_11*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_12*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_13*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_14*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_15*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_16*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_17*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_18*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_19*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_2*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_20*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_21*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_3*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_4*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_5*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_6*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_7*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_8*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_9*_data_edgedof_macroedge_src[micro_edge_index_x - 1];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_102(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_102_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_103_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_16 = v2v_cell_stencil[{ -1, 0, 0 }];
   const double xi_17 = v2v_cell_stencil[{ -1, 0, 1 }];
   const double xi_18 = v2v_cell_stencil[{ -1, 1, 0 }];
   const double xi_19 = v2v_cell_stencil[{ 0, 0, 1 }];
   const double xi_20 = v2v_cell_stencil[{ 0, 1, 0 }];
   const double xi_21 = v2v_cell_stencil[{ 1, 0, 0 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_11*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_12*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_13*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_14*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_15*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_16*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_17*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_18*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_19*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_2*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_20*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_21*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_3*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_4*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_5*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_6*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_7*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_8*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_9*_data_edgedof_macroedge_src[micro_edge_index_x - 1];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_103(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_103_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_120_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_16 = v2v_cell_stencil[{ -1, 0, 0 }];
   const double xi_17 = v2v_cell_stencil[{ -1, 0, 1 }];
   const double xi_18 = v2v_cell_stencil[{ -1, 1, 0 }];
   const double xi_19 = v2v_cell_stencil[{ 0, -1, 0 }];
   const double xi_20 = v2v_cell_stencil[{ 0, -1, 1 }];
   const double xi_21 = v2v_cell_stencil[{ 1, -1, 0 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_11*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_12*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_13*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_14*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_15*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_16*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_17*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_18*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_19*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_2*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_20*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_21*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_3*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_4*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_5*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_6*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_7*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_8*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_9*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_120(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_120_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_123_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_16 = v2v_cell_stencil[{ -1, 0, 0 }];
   const double xi_17 = v2v_cell_stencil[{ -1, 0, 1 }];
   const double xi_18 = v2v_cell_stencil[{ -1, 1, 0 }];
   const double xi_19 = v2v_cell_stencil[{ 0, -1, 0 }];
   const double xi_20 = v2v_cell_stencil[{ 0, -1, 1 }];
   const double xi_21 = v2v_cell_stencil[{ 1, -1, 0 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_11*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_12*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_13*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_14*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_15*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_16*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_17*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_18*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_19*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_2*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_20*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_21*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_3*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_4*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_5*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_6*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_7*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_8*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_9*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_123(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_123_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_130_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_16 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_17 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_18 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_19 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_20 = v2v_cell_stencil[{ -1, 0, 0 }];
   const double xi_21 = v2v_cell_stencil[{ -1, 0, 1 }];
   const double xi_22 = v2v_cell_stencil[{ -1, 1, -1 }];
   const double xi_23 = v2v_cell_stencil[{ -1, 1, 0 }];
   const double xi_24 = v2v_cell_stencil[{ 0, 0, -1 }];
   const double xi_25 = v2v_cell_stencil[{ 0, 1, -1 }];
   const double xi_26 = v2v_cell_stencil[{ 1, 0, -1 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 2*(1 << (level)) - 3] + xi_11*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_12*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 5*(1 << (level)) - 6] + xi_13*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_14*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_15*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_16*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_17*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_18*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 4*(1 << (level)) - 5] + xi_19*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_2*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 3*(1 << (level)) - 4] + xi_20*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_21*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_22*_data_vertexdof_macroedge_src[id_cell*((1 << (level)) - 1) + micro_edge_index_x + num_neighbor_faces*(1 << (level)) + (1 << (level))] + xi_23*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_24*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_25*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_26*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_3*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_4*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_5*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_6*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_7*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_8*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_9*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_130(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_130_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_132_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_16 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_17 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_18 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_19 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_20 = v2v_cell_stencil[{ -1, 0, 0 }];
   const double xi_21 = v2v_cell_stencil[{ -1, 0, 1 }];
   const double xi_22 = v2v_cell_stencil[{ -1, 1, -1 }];
   const double xi_23 = v2v_cell_stencil[{ -1, 1, 0 }];
   const double xi_24 = v2v_cell_stencil[{ 0, 0, -1 }];
   const double xi_25 = v2v_cell_stencil[{ 0, 1, -1 }];
   const double xi_26 = v2v_cell_stencil[{ 1, 0, -1 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 3*(1 << (level)) - 4] + xi_11*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_12*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 4*(1 << (level)) - 5] + xi_13*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_14*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_15*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_16*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_17*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_18*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 5*(1 << (level)) - 6] + xi_19*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_2*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 2*(1 << (level)) - 3] + xi_20*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_21*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_22*_data_vertexdof_macroedge_src[id_cell*((1 << (level)) - 1) + micro_edge_index_x + num_neighbor_faces*(1 << (level)) + (1 << (level))] + xi_23*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_24*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_25*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_26*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_3*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_4*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_5*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_6*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_7*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_8*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_9*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_132(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_132_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_201_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_16 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_17 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_18 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_19 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_20 = v2v_cell_stencil[{ 0, -1, 0 }];
   const double xi_21 = v2v_cell_stencil[{ 0, -1, 1 }];
   const double xi_22 = v2v_cell_stencil[{ 0, 0, 1 }];
   const double xi_23 = v2v_cell_stencil[{ 0, 1, 0 }];
   const double xi_24 = v2v_cell_stencil[{ 1, -1, 0 }];
   const double xi_25 = v2v_cell_stencil[{ 1, -1, 1 }];
   const double xi_26 = v2v_cell_stencil[{ 1, 0, 0 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_11*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_12*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 3*(1 << (level)) - 4] + xi_13*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_14*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_15*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_16*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_17*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_18*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_19*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 5*(1 << (level)) - 6] + xi_2*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_20*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_21*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_22*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_23*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_24*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_25*_data_vertexdof_macroedge_src[id_cell*((1 << (level)) - 1) + micro_edge_index_x + num_neighbor_faces*(1 << (level)) + (1 << (level))] + xi_26*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_3*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 2*(1 << (level)) - 3] + xi_4*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_5*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_6*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_7*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_8*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 4*(1 << (level)) - 5] + xi_9*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_201(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_201_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_203_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_16 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_17 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_18 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_19 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_20 = v2v_cell_stencil[{ 0, -1, 0 }];
   const double xi_21 = v2v_cell_stencil[{ 0, -1, 1 }];
   const double xi_22 = v2v_cell_stencil[{ 0, 0, 1 }];
   const double xi_23 = v2v_cell_stencil[{ 0, 1, 0 }];
   const double xi_24 = v2v_cell_stencil[{ 1, -1, 0 }];
   const double xi_25 = v2v_cell_stencil[{ 1, -1, 1 }];
   const double xi_26 = v2v_cell_stencil[{ 1, 0, 0 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_11*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_12*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 2*(1 << (level)) - 3] + xi_13*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_14*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_15*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_16*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_17*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_18*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_19*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 4*(1 << (level)) - 5] + xi_2*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_20*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_21*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_22*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_23*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_24*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_25*_data_vertexdof_macroedge_src[id_cell*((1 << (level)) - 1) + micro_edge_index_x + num_neighbor_faces*(1 << (level)) + (1 << (level))] + xi_26*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_3*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 3*(1 << (level)) - 4] + xi_4*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_5*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_6*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_7*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_8*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 5*(1 << (level)) - 6] + xi_9*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_203(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_203_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_210_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_16 = v2v_cell_stencil[{ -1, 0, 0 }];
   const double xi_17 = v2v_cell_stencil[{ -1, 0, 1 }];
   const double xi_18 = v2v_cell_stencil[{ -1, 1, 0 }];
   const double xi_19 = v2v_cell_stencil[{ 0, -1, 0 }];
   const double xi_20 = v2v_cell_stencil[{ 0, -1, 1 }];
   const double xi_21 = v2v_cell_stencil[{ 1, -1, 0 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_11*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_12*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_13*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_14*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_15*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_16*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_17*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_18*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_19*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_2*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_20*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_21*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_3*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_4*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_5*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_6*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_7*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_8*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_9*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_210(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_210_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_213_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_16 = v2v_cell_stencil[{ -1, 0, 0 }];
   const double xi_17 = v2v_cell_stencil[{ -1, 0, 1 }];
   const double xi_18 = v2v_cell_stencil[{ -1, 1, 0 }];
   const double xi_19 = v2v_cell_stencil[{ 0, -1, 0 }];
   const double xi_20 = v2v_cell_stencil[{ 0, -1, 1 }];
   const double xi_21 = v2v_cell_stencil[{ 1, -1, 0 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_11*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_12*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_13*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_14*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_15*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_16*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_17*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_18*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_19*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_2*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_20*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_21*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_3*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_4*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_5*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_6*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_7*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_8*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_9*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_213(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_213_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_230_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_16 = v2v_cell_stencil[{ 0, -1, 0 }];
   const double xi_17 = v2v_cell_stencil[{ 0, -1, 1 }];
   const double xi_18 = v2v_cell_stencil[{ 0, 0, -1 }];
   const double xi_19 = v2v_cell_stencil[{ 0, 1, -1 }];
   const double xi_20 = v2v_cell_stencil[{ 1, -1, 0 }];
   const double xi_21 = v2v_cell_stencil[{ 1, 0, -1 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_11*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_12*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_13*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_14*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_15*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_16*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_17*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_18*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_19*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_2*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_20*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_21*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_3*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_4*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_5*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_6*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_7*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_8*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_9*_data_edgedof_macroedge_src[micro_edge_index_x];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_230(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_230_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_231_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_16 = v2v_cell_stencil[{ 0, -1, 0 }];
   const double xi_17 = v2v_cell_stencil[{ 0, -1, 1 }];
   const double xi_18 = v2v_cell_stencil[{ 0, 0, -1 }];
   const double xi_19 = v2v_cell_stencil[{ 0, 1, -1 }];
   const double xi_20 = v2v_cell_stencil[{ 1, -1, 0 }];
   const double xi_21 = v2v_cell_stencil[{ 1, 0, -1 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_11*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_12*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_13*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_14*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_15*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_16*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_17*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_18*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_19*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_2*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_20*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_21*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_3*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_4*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_5*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_6*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_7*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_8*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_9*_data_edgedof_macroedge_src[micro_edge_index_x];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_231(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_231_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_301_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_16 = v2v_cell_stencil[{ 0, 0, -1 }];
   const double xi_17 = v2v_cell_stencil[{ 0, 0, 1 }];
   const double xi_18 = v2v_cell_stencil[{ 0, 1, -1 }];
   const double xi_19 = v2v_cell_stencil[{ 0, 1, 0 }];
   const double xi_20 = v2v_cell_stencil[{ 1, 0, -1 }];
   const double xi_21 = v2v_cell_stencil[{ 1, 0, 0 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_11*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_12*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_13*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_14*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_15*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_16*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_17*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_18*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_19*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_2*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_20*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_21*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_3*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_4*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_5*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_6*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_7*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_8*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_9*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_301(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_301_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_302_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_16 = v2v_cell_stencil[{ 0, 0, -1 }];
   const double xi_17 = v2v_cell_stencil[{ 0, 0, 1 }];
   const double xi_18 = v2v_cell_stencil[{ 0, 1, -1 }];
   const double xi_19 = v2v_cell_stencil[{ 0, 1, 0 }];
   const double xi_20 = v2v_cell_stencil[{ 1, 0, -1 }];
   const double xi_21 = v2v_cell_stencil[{ 1, 0, 0 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_11*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_12*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_13*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_14*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_15*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_16*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_17*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_18*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_19*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_2*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_20*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_21*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_3*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_4*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_5*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_6*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_7*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_8*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_9*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_302(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_302_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_310_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_16 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_17 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_18 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_19 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_20 = v2v_cell_stencil[{ -1, 0, 0 }];
   const double xi_21 = v2v_cell_stencil[{ -1, 0, 1 }];
   const double xi_22 = v2v_cell_stencil[{ -1, 1, -1 }];
   const double xi_23 = v2v_cell_stencil[{ -1, 1, 0 }];
   const double xi_24 = v2v_cell_stencil[{ 0, 0, -1 }];
   const double xi_25 = v2v_cell_stencil[{ 0, 1, -1 }];
   const double xi_26 = v2v_cell_stencil[{ 1, 0, -1 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 4*(1 << (level)) - 5] + xi_11*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_12*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 3*(1 << (level)) - 4] + xi_13*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_14*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_15*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_16*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_17*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_18*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 2*(1 << (level)) - 3] + xi_19*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_2*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 5*(1 << (level)) - 6] + xi_20*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_21*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_22*_data_vertexdof_macroedge_src[id_cell*((1 << (level)) - 1) + micro_edge_index_x + num_neighbor_faces*(1 << (level)) + (1 << (level))] + xi_23*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_24*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_25*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_26*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_3*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_4*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_5*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_6*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_7*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_8*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_9*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_310(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_310_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_312_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_16 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_17 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_18 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_19 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_20 = v2v_cell_stencil[{ -1, 0, 0 }];
   const double xi_21 = v2v_cell_stencil[{ -1, 0, 1 }];
   const double xi_22 = v2v_cell_stencil[{ -1, 1, -1 }];
   const double xi_23 = v2v_cell_stencil[{ -1, 1, 0 }];
   const double xi_24 = v2v_cell_stencil[{ 0, 0, -1 }];
   const double xi_25 = v2v_cell_stencil[{ 0, 1, -1 }];
   const double xi_26 = v2v_cell_stencil[{ 1, 0, -1 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 5*(1 << (level)) - 6] + xi_11*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_12*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 2*(1 << (level)) - 3] + xi_13*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_14*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_15*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_16*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_17*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_18*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 3*(1 << (level)) - 4] + xi_19*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_2*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 4*(1 << (level)) - 5] + xi_20*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_21*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_22*_data_vertexdof_macroedge_src[id_cell*((1 << (level)) - 1) + micro_edge_index_x + num_neighbor_faces*(1 << (level)) + (1 << (level))] + xi_23*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_24*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_25*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_26*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_3*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_4*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_5*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_6*_data_edgedof_macroedge_src[micro_edge_index_x - 1] + xi_7*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_8*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_9*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_312(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_312_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_320_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_16 = v2v_cell_stencil[{ 0, -1, 0 }];
   const double xi_17 = v2v_cell_stencil[{ 0, -1, 1 }];
   const double xi_18 = v2v_cell_stencil[{ 0, 0, -1 }];
   const double xi_19 = v2v_cell_stencil[{ 0, 1, -1 }];
   const double xi_20 = v2v_cell_stencil[{ 1, -1, 0 }];
   const double xi_21 = v2v_cell_stencil[{ 1, 0, -1 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_11*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_12*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_13*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_14*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_15*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_16*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_17*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_18*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_19*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_2*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_20*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_21*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_3*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_4*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_5*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_6*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_7*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_8*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_9*_data_edgedof_macroedge_src[micro_edge_index_x - 1];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_320(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_320_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    
static void sor_3D_macroedge_P2_update_vertexdofs_impl_321_level_any(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
   const double xi_1 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_2 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_3 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_4 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_5 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_6 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_7 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_8 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_9 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_10 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_11 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_12 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_13 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_14 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_15 = e2v_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_16 = v2v_cell_stencil[{ 0, -1, 0 }];
   const double xi_17 = v2v_cell_stencil[{ 0, -1, 1 }];
   const double xi_18 = v2v_cell_stencil[{ 0, 0, -1 }];
   const double xi_19 = v2v_cell_stencil[{ 0, 1, -1 }];
   const double xi_20 = v2v_cell_stencil[{ 1, -1, 0 }];
   const double xi_21 = v2v_cell_stencil[{ 1, 0, -1 }];
   _data_vertexdof_macroedge_dst[0] = xi_1*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7] + xi_10*_data_edgedof_macroedge_src[micro_edge_index_x] + xi_11*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_12*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_13*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_14*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_15*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_16*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_17*_data_vertexdof_macroedge_src[micro_edge_index_x - 1] + xi_18*_data_vertexdof_macroedge_src[id_face_1*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_19*_data_vertexdof_macroedge_src[micro_edge_index_x + 1] + xi_2*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 2] + xi_20*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level))] + xi_21*_data_vertexdof_macroedge_src[id_face_0*(1 << (level)) + micro_edge_index_x + (1 << (level)) + 1] + xi_3*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 2*(1 << (level)) - 1] + xi_4*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 2] + xi_5*_data_edgedof_macroedge_src[id_face_0*(3*(1 << (level)) - 1) + micro_edge_index_x + 3*(1 << (level)) - 1] + xi_6*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8] + xi_7*_data_edgedof_macroedge_src[id_cell*(7*(1 << (level)) - 7) + micro_edge_index_x + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7] + xi_8*_data_edgedof_macroedge_src[id_face_1*(3*(1 << (level)) - 1) + micro_edge_index_x + (1 << (level)) - 1] + xi_9*_data_edgedof_macroedge_src[micro_edge_index_x - 1];
   {
      
   }
}


void sor_3D_macroedge_P2_update_vertexdofs_impl_321(double const * RESTRICT const _data_edgedof_macroedge_src, double * RESTRICT _data_vertexdof_macroedge_dst, double const * RESTRICT const _data_vertexdof_macroedge_src, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2v_cell_stencil, int64_t id_cell, int64_t id_face_0, int64_t id_face_1, int level, int64_t micro_edge_index_x, int64_t num_neighbor_faces, std::map< hyteg::indexing::Index, double > v2v_cell_stencil)
{
    switch( level )
    {

    default:
        sor_3D_macroedge_P2_update_vertexdofs_impl_321_level_any(_data_edgedof_macroedge_src, _data_vertexdof_macroedge_dst, _data_vertexdof_macroedge_src, e2v_cell_stencil, id_cell, id_face_0, id_face_1, level, micro_edge_index_x, num_neighbor_faces, v2v_cell_stencil);
        break;
    }
}
    

} // namespace generated
} // namespace macroedge
} // namespace P2
} // namespace hyteg
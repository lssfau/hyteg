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

#include "communicate_directly_edgedof_cell_to_face_part_2_impl.hpp"

namespace hyteg {
namespace edgedof {
namespace comm {
namespace generated {

static void communicate_directly_edgedof_cell_to_face_part_2_impl_012_level_any(double const * RESTRICT const _data_edge_cell_src_XZ, double const * RESTRICT const _data_edge_cell_src_YZ, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_012(double const * RESTRICT const _data_edge_cell_src_XZ, double const * RESTRICT const _data_edge_cell_src_YZ, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_012_level_any(_data_edge_cell_src_XZ, _data_edge_cell_src_YZ, _data_edge_cell_src_Z, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_013_level_any(double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_YZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_013(double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_YZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_013_level_any(_data_edge_cell_src_XY, _data_edge_cell_src_Y, _data_edge_cell_src_YZ, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_021_level_any(double const * RESTRICT const _data_edge_cell_src_XZ, double const * RESTRICT const _data_edge_cell_src_YZ, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_021(double const * RESTRICT const _data_edge_cell_src_XZ, double const * RESTRICT const _data_edge_cell_src_YZ, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_021_level_any(_data_edge_cell_src_XZ, _data_edge_cell_src_YZ, _data_edge_cell_src_Z, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_023_level_any(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_XZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_023(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_XZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_023_level_any(_data_edge_cell_src_X, _data_edge_cell_src_XY, _data_edge_cell_src_XZ, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_031_level_any(double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_YZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_031(double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_YZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_031_level_any(_data_edge_cell_src_XY, _data_edge_cell_src_Y, _data_edge_cell_src_YZ, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_032_level_any(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_XZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(-ctr_1 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_032(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_XZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_032_level_any(_data_edge_cell_src_X, _data_edge_cell_src_XY, _data_edge_cell_src_XZ, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_102_level_any(double const * RESTRICT const _data_edge_cell_src_XZ, double const * RESTRICT const _data_edge_cell_src_YZ, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + (1 << (level)) - 1];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_102(double const * RESTRICT const _data_edge_cell_src_XZ, double const * RESTRICT const _data_edge_cell_src_YZ, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_102_level_any(_data_edge_cell_src_XZ, _data_edge_cell_src_YZ, _data_edge_cell_src_Z, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_103_level_any(double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_YZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_103(double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_YZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_103_level_any(_data_edge_cell_src_XY, _data_edge_cell_src_Y, _data_edge_cell_src_YZ, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_120_level_any(double const * RESTRICT const _data_edge_cell_src_XZ, double const * RESTRICT const _data_edge_cell_src_YZ, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + (1 << (level)) - 1];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_120(double const * RESTRICT const _data_edge_cell_src_XZ, double const * RESTRICT const _data_edge_cell_src_YZ, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_120_level_any(_data_edge_cell_src_XZ, _data_edge_cell_src_YZ, _data_edge_cell_src_Z, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_123_level_any(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*(-ctr_2 + (1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_123(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_123_level_any(_data_edge_cell_src_X, _data_edge_cell_src_Y, _data_edge_cell_src_Z, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_130_level_any(double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_YZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 - ctr_2 - ((0) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_130(double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_YZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_130_level_any(_data_edge_cell_src_XY, _data_edge_cell_src_Y, _data_edge_cell_src_YZ, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_132_level_any(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_132(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_132_level_any(_data_edge_cell_src_X, _data_edge_cell_src_Y, _data_edge_cell_src_Z, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_201_level_any(double const * RESTRICT const _data_edge_cell_src_XZ, double const * RESTRICT const _data_edge_cell_src_YZ, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_201(double const * RESTRICT const _data_edge_cell_src_XZ, double const * RESTRICT const _data_edge_cell_src_YZ, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_201_level_any(_data_edge_cell_src_XZ, _data_edge_cell_src_YZ, _data_edge_cell_src_Z, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_203_level_any(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_XZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[(-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_203(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_XZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_203_level_any(_data_edge_cell_src_X, _data_edge_cell_src_XY, _data_edge_cell_src_XZ, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_210_level_any(double const * RESTRICT const _data_edge_cell_src_XZ, double const * RESTRICT const _data_edge_cell_src_YZ, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_210(double const * RESTRICT const _data_edge_cell_src_XZ, double const * RESTRICT const _data_edge_cell_src_YZ, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_210_level_any(_data_edge_cell_src_XZ, _data_edge_cell_src_YZ, _data_edge_cell_src_Z, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_213_level_any(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + (-ctr_2 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_2 + (1 << (level)))*(-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_213(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_213_level_any(_data_edge_cell_src_X, _data_edge_cell_src_Y, _data_edge_cell_src_Z, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_230_level_any(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_XZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_230(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_XZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_230_level_any(_data_edge_cell_src_X, _data_edge_cell_src_XY, _data_edge_cell_src_XZ, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_231_level_any(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2 + (-ctr_1 + (1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) - (((-ctr_1 + (1 << (level)))*(-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_231(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_231_level_any(_data_edge_cell_src_X, _data_edge_cell_src_Y, _data_edge_cell_src_Z, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_301_level_any(double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_YZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_301(double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_YZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_301_level_any(_data_edge_cell_src_XY, _data_edge_cell_src_Y, _data_edge_cell_src_YZ, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_302_level_any(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_XZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_302(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_XZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_302_level_any(_data_edge_cell_src_X, _data_edge_cell_src_XY, _data_edge_cell_src_XZ, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_310_level_any(double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_YZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_YZ[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 - ((0) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_310(double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_YZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_310_level_any(_data_edge_cell_src_XY, _data_edge_cell_src_Y, _data_edge_cell_src_YZ, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_312_level_any(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_312(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_312_level_any(_data_edge_cell_src_X, _data_edge_cell_src_Y, _data_edge_cell_src_Z, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_320_level_any(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_XZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XZ[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_XY[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(ctr_1 + ctr_2 + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_320(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_XY, double const * RESTRICT const _data_edge_cell_src_XZ, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_320_level_any(_data_edge_cell_src_X, _data_edge_cell_src_XY, _data_edge_cell_src_XZ, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    
static void communicate_directly_edgedof_cell_to_face_part_2_impl_321_level_any(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_edge_face_dst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Z[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_Y[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
            _data_edge_face_dst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_edge_cell_src_X[ctr_1*(ctr_1 + ctr_2 + 2) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((ctr_1 + ctr_2 + 1)*(ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)) / (6))];
         }
      }
   }
}


void communicate_directly_edgedof_cell_to_face_part_2_impl_321(double const * RESTRICT const _data_edge_cell_src_X, double const * RESTRICT const _data_edge_cell_src_Y, double const * RESTRICT const _data_edge_cell_src_Z, double * RESTRICT _data_edge_face_dst_gl0_XZ, double * RESTRICT _data_edge_face_dst_gl0_YZ, double * RESTRICT _data_edge_face_dst_gl0_Z, int level)
{
    switch( level )
    {

    default:
        communicate_directly_edgedof_cell_to_face_part_2_impl_321_level_any(_data_edge_cell_src_X, _data_edge_cell_src_Y, _data_edge_cell_src_Z, _data_edge_face_dst_gl0_XZ, _data_edge_face_dst_gl0_YZ, _data_edge_face_dst_gl0_Z, level);
        break;
    }
}
    

} // namespace generated
} // namespace comm
} // namespace edgedof
} // namespace hyteg
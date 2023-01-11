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

#include "prolongate_2D_macroface_P2_push_from_edgedofs_level_0_to_1.hpp"

namespace hyteg {
namespace P2 {
namespace macroface {
namespace generated {

static void prolongate_2D_macroface_P2_push_from_edgedofs_level_0_to_1_level_any(real_t * RESTRICT _data_edgeCoarseSrc_X, real_t * RESTRICT _data_edgeCoarseSrc_XY, real_t * RESTRICT _data_edgeCoarseSrc_Y, real_t * RESTRICT _data_edgeFineDst_X, real_t * RESTRICT _data_edgeFineDst_XY, real_t * RESTRICT _data_edgeFineDst_Y, real_t * RESTRICT _data_vertexFineDst, int coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
   const double xi_5 = 1 / (num_neighbor_faces_edge0);
   const double xi_6 = 1 / (num_neighbor_faces_edge1);
   const double xi_7 = 1 / (num_neighbor_faces_edge2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_5*0.75*_data_edgeCoarseSrc_X[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_5*0.75*_data_edgeCoarseSrc_X[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = 0.25*_data_edgeCoarseSrc_X[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = 0.5*_data_edgeCoarseSrc_X[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 0.5*_data_edgeCoarseSrc_X[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = 1.0*xi_5*_data_edgeCoarseSrc_X[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_7*0.75*_data_edgeCoarseSrc_XY[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_7*0.75*_data_edgeCoarseSrc_XY[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 0.25*_data_edgeCoarseSrc_XY[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = 0.5*_data_edgeCoarseSrc_XY[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = 0.5*_data_edgeCoarseSrc_XY[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1] = 1.0*xi_7*_data_edgeCoarseSrc_XY[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_6*0.75*_data_edgeCoarseSrc_Y[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_6*0.75*_data_edgeCoarseSrc_Y[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 0.5*_data_edgeCoarseSrc_Y[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = 0.25*_data_edgeCoarseSrc_Y[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = 0.5*_data_edgeCoarseSrc_Y[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = 1.0*xi_6*_data_edgeCoarseSrc_Y[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void prolongate_2D_macroface_P2_push_from_edgedofs_level_0_to_1(real_t * RESTRICT _data_edgeCoarseSrc_X, real_t * RESTRICT _data_edgeCoarseSrc_XY, real_t * RESTRICT _data_edgeCoarseSrc_Y, real_t * RESTRICT _data_edgeFineDst_X, real_t * RESTRICT _data_edgeFineDst_XY, real_t * RESTRICT _data_edgeFineDst_Y, real_t * RESTRICT _data_vertexFineDst, int coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
    switch( coarse_level )
    {

    default:
        prolongate_2D_macroface_P2_push_from_edgedofs_level_0_to_1_level_any(_data_edgeCoarseSrc_X, _data_edgeCoarseSrc_XY, _data_edgeCoarseSrc_Y, _data_edgeFineDst_X, _data_edgeFineDst_XY, _data_edgeFineDst_Y, _data_vertexFineDst, coarse_level, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hyteg
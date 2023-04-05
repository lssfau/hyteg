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

#include "restrict_2D_macroface_P2_update_edgedofs_level_0_to_1.hpp"

namespace hyteg {
namespace P2 {
namespace macroface {
namespace generated {

static void restrict_2D_macroface_P2_update_edgedofs_level_0_to_1_level_any(real_t * RESTRICT _data_edgeCoarseDst_X, real_t * RESTRICT _data_edgeCoarseDst_XY, real_t * RESTRICT _data_edgeCoarseDst_Y, real_t * RESTRICT _data_edgeFineSrc_X, real_t * RESTRICT _data_edgeFineSrc_XY, real_t * RESTRICT _data_edgeFineSrc_Y, real_t * RESTRICT _data_vertexFineSrc, int coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
   const double xi_6 = 1 / (num_neighbor_faces_edge0);
   const double xi_7 = 1 / (num_neighbor_faces_edge1);
   const double xi_8 = 1 / (num_neighbor_faces_edge2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_10 = 0.25*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_23 = 0.5*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_18 = 0.5*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_16 = 0.25*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_24 = 0.5*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_22 = 0.25*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_13 = 1.0*xi_6*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_14 = xi_6*0.75*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_15 = xi_6*0.75*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_25 = 1.0*xi_7*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_26 = xi_7*0.75*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_27 = xi_7*0.75*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_19 = 1.0*xi_8*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
         const double xi_20 = xi_8*0.75*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_21 = xi_8*0.75*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         _data_edgeCoarseDst_X[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_10 + xi_13 + xi_14 + xi_15 + xi_18 + xi_23 );
         _data_edgeCoarseDst_XY[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_16 + xi_18 + xi_19 + xi_20 + xi_21 + xi_24 );
         _data_edgeCoarseDst_Y[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 );
      }
   }
}


void restrict_2D_macroface_P2_update_edgedofs_level_0_to_1(real_t * RESTRICT _data_edgeCoarseDst_X, real_t * RESTRICT _data_edgeCoarseDst_XY, real_t * RESTRICT _data_edgeCoarseDst_Y, real_t * RESTRICT _data_edgeFineSrc_X, real_t * RESTRICT _data_edgeFineSrc_XY, real_t * RESTRICT _data_edgeFineSrc_Y, real_t * RESTRICT _data_vertexFineSrc, int coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
    switch( coarse_level )
    {

    default:
        restrict_2D_macroface_P2_update_edgedofs_level_0_to_1_level_any(_data_edgeCoarseDst_X, _data_edgeCoarseDst_XY, _data_edgeCoarseDst_Y, _data_edgeFineSrc_X, _data_edgeFineSrc_XY, _data_edgeFineSrc_Y, _data_vertexFineSrc, coarse_level, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hyteg
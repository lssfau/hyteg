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

#include "restrict_2D_macroface_P1_pull_additive.hpp"

namespace hyteg {
namespace vertexdof {
namespace macroface {
namespace generated {

static void restrict_2D_macroface_P1_pull_additive_level_any(real_t * RESTRICT _data_vertexCoarseDst, real_t const * RESTRICT const _data_vertexFineSrc, int coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_13 = 1 / (num_neighbor_faces_edge0);
   const double xi_14 = 1 / (num_neighbor_faces_edge1);
   const double xi_15 = 1 / (num_neighbor_faces_vertex0);
   const double xi_24 = 1 / (num_neighbor_faces_edge2);
   const double xi_25 = 1 / (num_neighbor_faces_vertex1);
   const double xi_35 = 1 / (num_neighbor_faces_vertex2);
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_17 = xi_13*0.5*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_18 = xi_14*0.5*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_19 = 1.0*xi_15*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            _data_vertexCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_17 + xi_18 + xi_19;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (coarse_level)); ctr_1 += 1)
         {
            const double xi_49 = 0.5*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_50 = 0.5*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_51 = xi_13*0.5*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_52 = xi_13*0.5*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_53 = 1.0*xi_13*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            _data_vertexCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_49 + xi_50 + xi_51 + xi_52 + xi_53;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (coarse_level)); ctr_1 < (1 << (coarse_level)) + 1; ctr_1 += 1)
         {
            const double xi_27 = xi_13*0.5*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_28 = xi_24*0.5*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_29 = 1.0*xi_25*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            _data_vertexCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_27 + xi_28 + xi_29;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (coarse_level)); ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_59 = 0.5*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_60 = 0.5*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_61 = xi_14*0.5*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_62 = xi_14*0.5*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_63 = 1.0*xi_14*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            _data_vertexCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_59 + xi_60 + xi_61 + xi_62 + xi_63;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (coarse_level)); ctr_1 += 1)
         {
            const double xi_3 = 0.5*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_4 = 0.5*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_5 = 0.5*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_6 = 0.5*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_7 = 0.5*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_8 = 0.5*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_9 = 1.0*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            _data_vertexCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4 + xi_5 + xi_6 + xi_7 + xi_8 + xi_9;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (coarse_level)); ctr_1 < -ctr_2 + (1 << (coarse_level)) + 1; ctr_1 += 1)
         {
            const double xi_69 = 0.5*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_70 = 0.5*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_71 = xi_24*0.5*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_72 = xi_24*0.5*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_73 = 1.0*xi_24*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            _data_vertexCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_69 + xi_70 + xi_71 + xi_72 + xi_73;
         }
      }
      for (int ctr_2 = (1 << (coarse_level)); ctr_2 < (1 << (coarse_level)) + 1; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_37 = xi_14*0.5*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_38 = xi_24*0.5*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_39 = 1.0*xi_35*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            _data_vertexCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_37 + xi_38 + xi_39;
         }
      }
   }
}


void restrict_2D_macroface_P1_pull_additive(real_t * RESTRICT _data_vertexCoarseDst, real_t const * RESTRICT const _data_vertexFineSrc, int coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
    switch( coarse_level )
    {

    default:
        restrict_2D_macroface_P1_pull_additive_level_any(_data_vertexCoarseDst, _data_vertexFineSrc, coarse_level, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg
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

#include "prolongate_2D_macroface_P1_push_additive.hpp"

namespace hyteg {
namespace vertexdof {
namespace macroface {
namespace generated {

static void prolongate_2D_macroface_P1_push_additive_level_any(double const * RESTRICT const _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, int coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_18 = 1 / (num_neighbor_faces_edge0);
   const double xi_19 = 1 / (num_neighbor_faces_edge1);
   const double xi_20 = 1 / (num_neighbor_faces_vertex0);
   const double xi_30 = 1 / (num_neighbor_faces_edge2);
   const double xi_31 = 1 / (num_neighbor_faces_vertex1);
   const double xi_42 = 1 / (num_neighbor_faces_vertex2);
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_23 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_25 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_22 = xi_18*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_24 = xi_19*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_22 + xi_23;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_24 + xi_25;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_20*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (coarse_level)); ctr_1 += 1)
         {
            const double xi_58 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_61 = 0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_62 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_64 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_63 = xi_18*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_58 + xi_63;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_60 + xi_61;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1] = xi_61 + xi_62;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_63 + xi_64;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_18*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (coarse_level)); ctr_1 < (1 << (coarse_level)) + 1; ctr_1 += 1)
         {
            const double xi_34 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_36 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_35 = xi_18*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_33 = xi_30*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1] = xi_33 + xi_34;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_35 + xi_36;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_31*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (coarse_level)); ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_74 = 0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_71 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_73 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_75 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_77 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_76 = xi_19*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_71 + xi_74;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_73 + xi_76;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_74 + xi_75;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_76 + xi_77;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_19*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (coarse_level)); ctr_1 += 1)
         {
            const double xi_13 = 0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_4 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_6 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_8 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_10 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_12 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_14 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_13 + xi_4;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_13 + xi_6;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1] = xi_13 + xi_8;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_10 + xi_13;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_12 + xi_13;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_13 + xi_14;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (coarse_level)); ctr_1 < -ctr_2 + (1 << (coarse_level)) + 1; ctr_1 += 1)
         {
            const double xi_84 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_87 = 0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_86 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_88 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_90 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_89 = xi_30*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1] = xi_84 + xi_89;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_86 + xi_87;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_87 + xi_88;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_89 + xi_90;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_30*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = (1 << (coarse_level)); ctr_2 < (1 << (coarse_level)) + 1; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_45 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_47 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_44 = xi_19*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_46 = xi_30*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_44 + xi_45;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_46 + xi_47;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_42*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}


void prolongate_2D_macroface_P1_push_additive(double const * RESTRICT const _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, int coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
    switch( coarse_level )
    {

    default:
        prolongate_2D_macroface_P1_push_additive_level_any(_data_vertexCoarseSrc, _data_vertexFineDst, coarse_level, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg
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
   const double xi_20 = 1 / (num_neighbor_faces_edge0);
   const double xi_21 = 1 / (num_neighbor_faces_edge1);
   const double xi_22 = 1 / (num_neighbor_faces_vertex0);
   const double xi_34 = 1 / (num_neighbor_faces_edge2);
   const double xi_35 = 1 / (num_neighbor_faces_vertex1);
   const double xi_48 = 1 / (num_neighbor_faces_vertex2);
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_25 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_27 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_29 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_24 = xi_20*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_26 = xi_21*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_28 = 1.0*xi_22*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_24 + xi_25;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_26 + xi_27;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_28 + xi_29;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (coarse_level)); ctr_1 += 1)
         {
            const double xi_66 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_69 = 0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_68 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_70 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_72 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_74 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_71 = xi_20*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_73 = 1.0*xi_20*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_66 + xi_71;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_68 + xi_69;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1] = xi_69 + xi_70;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_71 + xi_72;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_73 + xi_74;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (coarse_level)); ctr_1 < (1 << (coarse_level)) + 1; ctr_1 += 1)
         {
            const double xi_38 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_40 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_42 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_39 = xi_20*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_37 = xi_34*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_41 = 1.0*xi_35*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1] = xi_37 + xi_38;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_39 + xi_40;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_41 + xi_42;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (coarse_level)); ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_84 = 0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_81 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_83 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_85 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_87 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_89 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_86 = xi_21*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_88 = 1.0*xi_21*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_81 + xi_84;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_83 + xi_86;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_84 + xi_85;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_86 + xi_87;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_88 + xi_89;
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
            const double xi_15 = 1.0*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_16 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_13 + xi_4;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_13 + xi_6;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1] = xi_13 + xi_8;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_10 + xi_13;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_12 + xi_13;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_13 + xi_14;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_15 + xi_16;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (coarse_level)); ctr_1 < -ctr_2 + (1 << (coarse_level)) + 1; ctr_1 += 1)
         {
            const double xi_96 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_99 = 0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_98 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_100 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_102 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_104 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_101 = xi_34*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_103 = 1.0*xi_34*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1] = xi_101 + xi_96;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_98 + xi_99;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_100 + xi_99;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_101 + xi_102;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_103 + xi_104;
         }
      }
      for (int ctr_2 = (1 << (coarse_level)); ctr_2 < (1 << (coarse_level)) + 1; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_51 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_53 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_55 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_50 = xi_21*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_52 = xi_34*0.5*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_54 = 1.0*xi_48*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_50 + xi_51;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_52 + xi_53;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_54 + xi_55;
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


static void prolongate_2D_macroface_P1_push_additive_level_any(float const * RESTRICT const _data_vertexCoarseSrc, float * RESTRICT _data_vertexFineDst, int coarse_level, float num_neighbor_faces_edge0, float num_neighbor_faces_edge1, float num_neighbor_faces_edge2, float num_neighbor_faces_vertex0, float num_neighbor_faces_vertex1, float num_neighbor_faces_vertex2)
{
    const float xi_20 = 1 / (num_neighbor_faces_edge0);
    const float xi_21 = 1 / (num_neighbor_faces_edge1);
    const float xi_22 = 1 / (num_neighbor_faces_vertex0);
    const float xi_34 = 1 / (num_neighbor_faces_edge2);
    const float xi_35 = 1 / (num_neighbor_faces_vertex1);
    const float xi_48 = 1 / (num_neighbor_faces_vertex2);
    {
        for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
        {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const float xi_25 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const float xi_27 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const float xi_29 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const float xi_24 = xi_20*static_cast< float >(0.5)*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const float xi_26 = xi_21*static_cast< float >(0.5)*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const float xi_28 = static_cast< float >(1.0)*xi_22*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_24 + xi_25;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_26 + xi_27;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_28 + xi_29;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (coarse_level)); ctr_1 += 1)
         {
            const float xi_66 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const float xi_69 = static_cast< float >(0.5)*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const float xi_68 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const float xi_70 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const float xi_72 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const float xi_74 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const float xi_71 = xi_20*static_cast< float >(0.5)*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const float xi_73 = static_cast< float >(1.0)*xi_20*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_66 + xi_71;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_68 + xi_69;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1] = xi_69 + xi_70;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_71 + xi_72;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_73 + xi_74;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (coarse_level)); ctr_1 < (1 << (coarse_level)) + 1; ctr_1 += 1)
         {
            const float xi_38 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const float xi_40 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const float xi_42 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const float xi_39 = xi_20*static_cast< float >(0.5)*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const float xi_37 = xi_34*static_cast< float >(0.5)*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const float xi_41 = static_cast< float >(1.0)*xi_35*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1] = xi_37 + xi_38;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_39 + xi_40;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_41 + xi_42;
         }
        }
        for (int ctr_2 = 1; ctr_2 < (1 << (coarse_level)); ctr_2 += 1)
        {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const float xi_84 = static_cast< float >(0.5)*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const float xi_81 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const float xi_83 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const float xi_85 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const float xi_87 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const float xi_89 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const float xi_86 = xi_21*static_cast< float >(0.5)*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const float xi_88 = static_cast< float >(1.0)*xi_21*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_81 + xi_84;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_83 + xi_86;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_84 + xi_85;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_86 + xi_87;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_88 + xi_89;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (coarse_level)); ctr_1 += 1)
         {
            const float xi_13 = static_cast< float >(0.5)*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const float xi_4 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const float xi_6 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const float xi_8 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const float xi_10 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const float xi_12 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const float xi_14 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const float xi_15 = static_cast< float >(1.0)*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const float xi_16 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_13 + xi_4;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_13 + xi_6;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1] = xi_13 + xi_8;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_10 + xi_13;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_12 + xi_13;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_13 + xi_14;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_15 + xi_16;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (coarse_level)); ctr_1 < -ctr_2 + (1 << (coarse_level)) + 1; ctr_1 += 1)
         {
            const float xi_96 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const float xi_99 = static_cast< float >(0.5)*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const float xi_98 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const float xi_100 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const float xi_102 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const float xi_104 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const float xi_101 = xi_34*static_cast< float >(0.5)*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const float xi_103 = static_cast< float >(1.0)*xi_34*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1] = xi_101 + xi_96;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_98 + xi_99;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_100 + xi_99;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_101 + xi_102;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_103 + xi_104;
         }
        }
        for (int ctr_2 = (1 << (coarse_level)); ctr_2 < (1 << (coarse_level)) + 1; ctr_2 += 1)
        {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const float xi_51 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const float xi_53 = _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const float xi_55 = _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const float xi_50 = xi_21*static_cast< float >(0.5)*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const float xi_52 = xi_34*static_cast< float >(0.5)*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const float xi_54 = static_cast< float >(1.0)*xi_48*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_50 + xi_51;
            _data_vertexFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_52 + xi_53;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_54 + xi_55;
         }
        }
    }
}


void prolongate_2D_macroface_P1_push_additive(float const * RESTRICT const _data_vertexCoarseSrc, float * RESTRICT _data_vertexFineDst, int coarse_level, float num_neighbor_faces_edge0, float num_neighbor_faces_edge1, float num_neighbor_faces_edge2, float num_neighbor_faces_vertex0, float num_neighbor_faces_vertex1, float num_neighbor_faces_vertex2)
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
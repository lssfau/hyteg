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

#include "restrict_2D_macroface_P2_update_vertexdofs.hpp"

namespace hyteg {
namespace P2 {
namespace macroface {
namespace generated {

static void restrict_2D_macroface_P2_update_vertexdofs_level_any(real_t * RESTRICT _data_edgeFineSrc_X, real_t * RESTRICT _data_edgeFineSrc_XY, real_t * RESTRICT _data_edgeFineSrc_Y, real_t * RESTRICT _data_vertexCoarseDst, real_t * RESTRICT _data_vertexFineSrc, int coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_91 = 1 / (num_neighbor_faces_edge0);
   const double xi_92 = 1 / (num_neighbor_faces_edge1);
   const double xi_93 = 1 / (num_neighbor_faces_vertex0);
   const double xi_106 = 1 / (num_neighbor_faces_edge2);
   const double xi_107 = 1 / (num_neighbor_faces_vertex1);
   const double xi_121 = 1 / (num_neighbor_faces_vertex2);
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_95 = -0.125*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_96 = -0.125*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_98 = xi_91*0.375*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_100 = xi_91*-0.125*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_99 = xi_92*0.375*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_101 = xi_92*-0.125*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_97 = 1.0*xi_93*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            _data_vertexCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_101 + xi_95 + xi_96 + xi_97 + xi_98 + xi_99;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (coarse_level)); ctr_1 += 1)
         {
            const double xi_33 = 0.375*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_34 = 0.375*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_35 = -0.125*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_36 = -0.125*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_37 = -0.125*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_38 = -0.125*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_39 = -0.125*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_40 = -0.125*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_41 = -0.125*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_42 = -0.125*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_43 = 1.0*xi_91*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_44 = xi_91*0.375*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_45 = xi_91*0.375*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_46 = xi_91*-0.125*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_47 = xi_91*-0.125*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            _data_vertexCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (coarse_level)); ctr_1 < (1 << (coarse_level)) + 1; ctr_1 += 1)
         {
            const double xi_109 = -0.125*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_110 = -0.125*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_112 = xi_91*0.375*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_114 = xi_91*-0.125*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_113 = xi_106*0.375*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_115 = xi_106*-0.125*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_111 = 1.0*xi_107*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            _data_vertexCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (coarse_level)); ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_53 = 0.375*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_54 = 0.375*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_55 = -0.125*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_56 = -0.125*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_57 = -0.125*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_58 = -0.125*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_59 = -0.125*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_60 = -0.125*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_61 = -0.125*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_62 = -0.125*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_63 = 1.0*xi_92*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_64 = xi_92*0.375*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_65 = xi_92*0.375*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_66 = xi_92*-0.125*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_67 = xi_92*-0.125*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            _data_vertexCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (coarse_level)); ctr_1 += 1)
         {
            const double xi_3 = 0.375*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_4 = 0.375*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_5 = 0.375*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_6 = 0.375*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_7 = 0.375*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_8 = 0.375*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_9 = -0.125*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_10 = -0.125*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_11 = -0.125*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_12 = -0.125*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_13 = -0.125*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_14 = -0.125*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_15 = -0.125*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1];
            const double xi_16 = -0.125*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_17 = -0.125*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_18 = -0.125*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_19 = -0.125*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_20 = -0.125*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_21 = -0.125*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1];
            const double xi_22 = -0.125*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_23 = -0.125*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_24 = -0.125*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_25 = -0.125*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_26 = -0.125*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_27 = _data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            _data_vertexCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_3 + xi_4 + xi_5 + xi_6 + xi_7 + xi_8 + xi_9;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (coarse_level)); ctr_1 < -ctr_2 + (1 << (coarse_level)) + 1; ctr_1 += 1)
         {
            const double xi_73 = 0.375*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_74 = 0.375*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_75 = -0.125*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_76 = -0.125*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_77 = -0.125*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1];
            const double xi_78 = -0.125*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_79 = -0.125*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_80 = -0.125*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1];
            const double xi_81 = -0.125*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_82 = -0.125*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_83 = 1.0*xi_106*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_84 = xi_106*0.375*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_85 = xi_106*0.375*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_86 = xi_106*-0.125*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_87 = xi_106*-0.125*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            _data_vertexCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_73 + xi_74 + xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80 + xi_81 + xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87;
         }
      }
      for (int ctr_2 = (1 << (coarse_level)); ctr_2 < (1 << (coarse_level)) + 1; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_123 = -0.125*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_124 = -0.125*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_126 = xi_92*0.375*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_128 = xi_92*-0.125*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_127 = xi_106*0.375*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_129 = xi_106*-0.125*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_125 = 1.0*xi_121*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            _data_vertexCoarseDst[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_123 + xi_124 + xi_125 + xi_126 + xi_127 + xi_128 + xi_129;
         }
      }
   }
}


void restrict_2D_macroface_P2_update_vertexdofs(real_t * RESTRICT _data_edgeFineSrc_X, real_t * RESTRICT _data_edgeFineSrc_XY, real_t * RESTRICT _data_edgeFineSrc_Y, real_t * RESTRICT _data_vertexCoarseDst, real_t * RESTRICT _data_vertexFineSrc, int coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
    switch( coarse_level )
    {

    default:
        restrict_2D_macroface_P2_update_vertexdofs_level_any(_data_edgeFineSrc_X, _data_edgeFineSrc_XY, _data_edgeFineSrc_Y, _data_vertexCoarseDst, _data_vertexFineSrc, coarse_level, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hyteg
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

#include "restrict_2D_macroface_P2_update_edgedofs.hpp"

namespace hyteg {
namespace P2 {
namespace macroface {
namespace generated {

static void restrict_2D_macroface_P2_update_edgedofs_level_any(real_t * RESTRICT _data_edgeCoarseDst_X, real_t * RESTRICT _data_edgeCoarseDst_XY, real_t * RESTRICT _data_edgeCoarseDst_Y, real_t * RESTRICT _data_edgeFineSrc_X, real_t * RESTRICT _data_edgeFineSrc_XY, real_t * RESTRICT _data_edgeFineSrc_Y, real_t * RESTRICT _data_vertexFineSrc, int coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
   const double xi_33 = 1 / (num_neighbor_faces_edge0);
   const double xi_34 = 1 / (num_neighbor_faces_edge1);
   const double xi_61 = 1 / (num_neighbor_faces_edge2);
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_36 = 0.25*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_52 = 0.5*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_46 = 0.5*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_42 = 0.25*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_43 = 0.25*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_53 = 0.5*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_45 = 0.5*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_47 = 0.5*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_48 = 0.75*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_49 = 0.75*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_50 = _data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_51 = 0.25*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_39 = 1.0*xi_33*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_40 = xi_33*0.75*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_41 = xi_33*0.75*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_54 = 1.0*xi_34*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_55 = xi_34*0.75*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_56 = xi_34*0.75*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            _data_edgeCoarseDst_X[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_36 + xi_39 + xi_40 + xi_41 + xi_46 + xi_52 );
            _data_edgeCoarseDst_XY[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_42 + xi_43 + xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_53 );
            _data_edgeCoarseDst_Y[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56 );
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (coarse_level)) - 1; ctr_1 += 1)
         {
            const double xi_141 = 0.25*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_159 = 0.5*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_151 = 0.5*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_147 = 0.25*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_148 = 0.25*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_161 = 0.5*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_150 = 0.5*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_152 = 0.5*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_153 = 0.75*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_154 = 0.75*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_155 = _data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_156 = 0.25*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_157 = 0.25*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_158 = 0.5*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_160 = 0.5*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_162 = 0.75*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_163 = 0.75*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_164 = _data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_144 = 1.0*xi_33*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_145 = xi_33*0.75*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_146 = xi_33*0.75*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            _data_edgeCoarseDst_X[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_141 + xi_144 + xi_145 + xi_146 + xi_151 + xi_159 );
            _data_edgeCoarseDst_XY[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_147 + xi_148 + xi_150 + xi_151 + xi_152 + xi_153 + xi_154 + xi_155 + xi_161 );
            _data_edgeCoarseDst_Y[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_156 + xi_157 + xi_158 + xi_159 + xi_160 + xi_161 + xi_162 + xi_163 + xi_164 );
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (coarse_level)) - 1; ctr_1 < (1 << (coarse_level)); ctr_1 += 1)
         {
            const double xi_63 = 0.25*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_78 = 0.5*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_71 = 0.5*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_69 = 0.25*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_80 = 0.5*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_75 = 0.25*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_76 = 0.25*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_77 = 0.5*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_79 = 0.5*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_81 = 0.75*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_82 = 0.75*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_83 = _data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_66 = 1.0*xi_33*_data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_67 = xi_33*0.75*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_68 = xi_33*0.75*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_72 = 1.0*xi_61*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_73 = xi_61*0.75*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_74 = xi_61*0.75*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            _data_edgeCoarseDst_X[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_63 + xi_66 + xi_67 + xi_68 + xi_71 + xi_78 );
            _data_edgeCoarseDst_XY[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_69 + xi_71 + xi_72 + xi_73 + xi_74 + xi_80 );
            _data_edgeCoarseDst_Y[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80 + xi_81 + xi_82 + xi_83 );
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (coarse_level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_170 = 0.25*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_171 = 0.25*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_189 = 0.5*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_173 = 0.5*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_174 = 0.5*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_183 = 0.5*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_176 = 0.75*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_177 = 0.75*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_178 = _data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_179 = 0.25*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_180 = 0.25*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_190 = 0.5*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_182 = 0.5*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_184 = 0.5*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_185 = 0.75*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_186 = 0.75*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_187 = _data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_188 = 0.25*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_191 = 1.0*xi_34*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_192 = xi_34*0.75*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_193 = xi_34*0.75*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            _data_edgeCoarseDst_X[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_170 + xi_171 + xi_173 + xi_174 + xi_176 + xi_177 + xi_178 + xi_183 + xi_189 );
            _data_edgeCoarseDst_XY[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_179 + xi_180 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_190 );
            _data_edgeCoarseDst_Y[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 );
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (coarse_level)) - 1; ctr_1 += 1)
         {
            const double xi_3 = 0.25*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_4 = 0.25*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_24 = 0.5*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_6 = 0.5*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_7 = 0.5*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_16 = 0.5*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_9 = 0.75*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_10 = 0.75*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_11 = _data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_12 = 0.25*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_13 = 0.25*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_26 = 0.5*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_15 = 0.5*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_17 = 0.5*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_18 = 0.75*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_19 = 0.75*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_20 = _data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_21 = 0.25*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_22 = 0.25*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_23 = 0.5*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_25 = 0.5*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_27 = 0.75*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_28 = 0.75*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_29 = _data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            _data_edgeCoarseDst_X[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_10 + xi_11 + xi_16 + xi_24 + xi_3 + xi_4 + xi_6 + xi_7 + xi_9 );
            _data_edgeCoarseDst_XY[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_12 + xi_13 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19 + xi_20 + xi_26 );
            _data_edgeCoarseDst_Y[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 );
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (coarse_level)) - 1; ctr_1 < -ctr_2 + (1 << (coarse_level)); ctr_1 += 1)
         {
            const double xi_199 = 0.25*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_200 = 0.25*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_217 = 0.5*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_202 = 0.5*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_203 = 0.5*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_210 = 0.5*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_205 = 0.75*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_206 = 0.75*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_207 = _data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_208 = 0.25*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_219 = 0.5*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_214 = 0.25*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_215 = 0.25*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_216 = 0.5*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_218 = 0.5*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_220 = 0.75*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_221 = 0.75*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_222 = _data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_211 = 1.0*xi_61*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_212 = xi_61*0.75*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_213 = xi_61*0.75*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            _data_edgeCoarseDst_X[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_199 + xi_200 + xi_202 + xi_203 + xi_205 + xi_206 + xi_207 + xi_210 + xi_217 );
            _data_edgeCoarseDst_XY[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_208 + xi_210 + xi_211 + xi_212 + xi_213 + xi_219 );
            _data_edgeCoarseDst_Y[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_214 + xi_215 + xi_216 + xi_217 + xi_218 + xi_219 + xi_220 + xi_221 + xi_222 );
         }
      }
      for (int ctr_2 = (1 << (coarse_level)) - 1; ctr_2 < (1 << (coarse_level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_90 = 0.25*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_91 = 0.25*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_106 = 0.5*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_93 = 0.5*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_94 = 0.5*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_101 = 0.5*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_96 = 0.75*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_97 = 0.75*_data_edgeFineSrc_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_98 = _data_vertexFineSrc[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_99 = 0.25*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_107 = 0.5*_data_edgeFineSrc_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_105 = 0.25*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_108 = 1.0*xi_34*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_109 = xi_34*0.75*_data_edgeFineSrc_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_110 = xi_34*0.75*_data_edgeFineSrc_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_102 = 1.0*xi_61*_data_vertexFineSrc[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 2) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1];
            const double xi_103 = xi_61*0.75*_data_edgeFineSrc_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_104 = xi_61*0.75*_data_edgeFineSrc_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            _data_edgeCoarseDst_X[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_101 + xi_106 + xi_90 + xi_91 + xi_93 + xi_94 + xi_96 + xi_97 + xi_98 );
            _data_edgeCoarseDst_XY[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_101 + xi_102 + xi_103 + xi_104 + xi_107 + xi_99 );
            _data_edgeCoarseDst_Y[ctr_1 + ctr_2*((1 << (coarse_level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = real_c( xi_105 + xi_106 + xi_107 + xi_108 + xi_109 + xi_110 );
         }
      }
   }
}


void restrict_2D_macroface_P2_update_edgedofs(real_t * RESTRICT _data_edgeCoarseDst_X, real_t * RESTRICT _data_edgeCoarseDst_XY, real_t * RESTRICT _data_edgeCoarseDst_Y, real_t * RESTRICT _data_edgeFineSrc_X, real_t * RESTRICT _data_edgeFineSrc_XY, real_t * RESTRICT _data_edgeFineSrc_Y, real_t * RESTRICT _data_vertexFineSrc, int coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2)
{
    switch( coarse_level )
    {

    default:
        restrict_2D_macroface_P2_update_edgedofs_level_any(_data_edgeCoarseDst_X, _data_edgeCoarseDst_XY, _data_edgeCoarseDst_Y, _data_edgeFineSrc_X, _data_edgeFineSrc_XY, _data_edgeFineSrc_Y, _data_vertexFineSrc, coarse_level, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hyteg
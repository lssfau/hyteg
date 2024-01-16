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

#include "apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl.hpp"

namespace hyteg {
namespace edgedof {
namespace macroface {
namespace generated {

static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_012_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_80 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, 0, 0 }];
   const double xi_81 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_82 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_199 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_200 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_201 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_202 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_203 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_204 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_85 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_87 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_89 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_76*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_91 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_94 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_95 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_96 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_176 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_177 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_178 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_179 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_180 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_181 = xi_76*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_182 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_183 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_184 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_185 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_186 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_187 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_188 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_189 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_190 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_192 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_193 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_195 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_196 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_118 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_110 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_111 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_112 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_113 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_114 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_115 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_116 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_117 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_225 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_226 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_227 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_228 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_231 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_233 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_235 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_236 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_237 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_239 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_240 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_241 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_76*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_243 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_244 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_245 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_246 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_247 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_248 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_225 + xi_226 + xi_227 + xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246 + xi_247 + xi_248;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_41 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_45 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_47 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_38 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_48 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_50 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_52 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_53 = xi_76*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_54 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_55 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_56 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_57 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_58 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_59 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_60 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_61 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_62 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_63 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_64 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_65 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_66 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_67 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_68 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_274 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_276 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_278 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_280 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_282 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_283 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_284 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_285 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_286 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_287 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_288 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_289 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_290 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_291 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_292 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_293 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_294 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282 + xi_283 + xi_284 + xi_285;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_135 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_136 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_137 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_138 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_141 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_142 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_143 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_145 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_012(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_012_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_013_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   const double xi_162 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_163 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_164 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_199 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_200 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_201 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_202 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_203 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_204 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_82 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_83 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_84 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_85 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_86 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_87 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_89 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_176 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_177 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_178 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_179 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_180 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_181 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_182 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_183 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_184 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_185 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_186 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_187 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_188 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_189 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_190 = xi_167*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_191 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_192 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_193 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_194 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_195 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_196 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_112 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_113 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_114 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_115 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_116 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_117 = xi_167*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_118 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_107 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_108 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_109 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_110 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_111 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_222 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_223 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_224 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_225 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_226 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_227 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_228 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_231 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_234 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_236 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_237 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_238 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_239 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_240 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_241 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_222 + xi_223 + xi_224 + xi_225 + xi_226 + xi_227 + xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_234 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_41 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_45 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_48 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_50 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_52 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_53 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_54 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_55 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_56 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_57 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_58 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_59 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_60 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_61 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_62 = xi_167*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_63 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_64 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_65 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_66 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_67 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_68 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_271 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_272 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_273 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_274 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_276 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_277 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_278 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_283 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_284 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_285 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_286 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_287 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_288 = xi_167*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_289 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_290 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_291 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_292 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_293 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_294 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_271 + xi_272 + xi_273 + xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_283 + xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_135 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_136 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_137 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_138 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_141 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_142 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_013(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_013_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_021_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_80 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, 0, 0 }];
   const double xi_81 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_82 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_174 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_175 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_176 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_211 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_212 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_213 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_85 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_87 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_89 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_91 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_94 = xi_80*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_95 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_96 = xi_82*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_179 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_180 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_181 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_182 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_183 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_184 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_185 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_186 = xi_78*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_187 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_188 = xi_80*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_189 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_190 = xi_82*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_192 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_193 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_194 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_195 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_196 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_197 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_198 = xi_172*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_199 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_200 = xi_174*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_201 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_202 = xi_176*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196 + xi_197 + xi_198 + xi_199 + xi_200 + xi_201 + xi_202;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_124 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_113 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_114 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_115 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_116 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_117 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_118 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_119 = xi_172*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_120 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_121 = xi_174*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_122 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_123 = xi_176*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122 + xi_123 + xi_124;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_228 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_231 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_233 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_236 = xi_213*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_237 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_239 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_240 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_241 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_243 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_244 = xi_78*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_245 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_246 = xi_80*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_247 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_248 = xi_82*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246 + xi_247 + xi_248;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_41 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_44 = xi_213*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_45 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_48 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_50 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_51 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_52 = xi_78*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_53 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_54 = xi_80*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_55 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_56 = xi_82*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_57 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_58 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_59 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_60 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_61 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_62 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_63 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_64 = xi_172*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_65 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_66 = xi_174*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_67 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_68 = xi_176*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_274 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_276 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_278 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_279 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_282 = xi_213*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_283 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_284 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_285 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_286 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_287 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_288 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_289 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_290 = xi_172*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_291 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_292 = xi_174*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_293 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_294 = xi_176*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_283 + xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_138 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_141 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_142 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_143 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_146 = xi_213*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_021(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_021_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_023_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 1 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 1 }];
   const double xi_80 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_81 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, 0 }];
   const double xi_82 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_174 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_175 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_176 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_211 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_212 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_213 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_85 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_87 = xi_73*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_89 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_76*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_91 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_94 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_95 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_96 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_179 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_180 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_181 = xi_73*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_182 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_183 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_184 = xi_76*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_185 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_186 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_187 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_188 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_189 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_190 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_192 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_193 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_194 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_195 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_196 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_197 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_198 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_199 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_200 = xi_174*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_201 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_202 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196 + xi_197 + xi_198 + xi_199 + xi_200 + xi_201 + xi_202;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_124 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_113 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_114 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_115 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_116 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_117 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_118 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_119 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_120 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_121 = xi_174*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_122 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_123 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122 + xi_123 + xi_124;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_228 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_231 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_235 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_236 = xi_213*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_237 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_239 = xi_73*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_240 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_241 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_76*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_243 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_244 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_245 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_246 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_247 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_248 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246 + xi_247 + xi_248;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_41 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_43 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_213*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_45 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_73*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_48 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_50 = xi_76*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_52 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_53 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_54 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_55 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_56 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_57 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_58 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_59 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_60 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_61 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_62 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_63 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_64 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_65 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_66 = xi_174*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_67 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_68 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_274 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_276 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_278 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_281 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_213*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_283 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_284 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_285 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_286 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_287 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_288 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_289 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_290 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_291 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_292 = xi_174*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_293 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_294 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_283 + xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_138 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_141 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_142 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_145 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_213*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_023(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_023_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_031_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   const double xi_162 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_163 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_164 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_199 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_200 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_201 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_202 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_203 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_204 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_82 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_83 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_84 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_85 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_87 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_89 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_90 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_176 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_177 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_178 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_179 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_180 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_181 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_182 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_183 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_184 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_185 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_186 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_187 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_188 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_189 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_190 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_192 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_193 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_194 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_195 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_196 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_112 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_113 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_114 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_115 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_116 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_117 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_118 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_107 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_108 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_109 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_110 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_111 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_222 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_223 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_224 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_225 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_226 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_227 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_228 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_231 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_234 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_236 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_237 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_239 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_240 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_241 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_242 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_222 + xi_223 + xi_224 + xi_225 + xi_226 + xi_227 + xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_234 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_41 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_45 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_48 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_50 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_52 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_53 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_54 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_55 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_56 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_57 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_58 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_59 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_60 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_61 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_62 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_63 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_64 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_65 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_66 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_67 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_68 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_271 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_272 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_273 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_274 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_276 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_277 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_278 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_283 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_284 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_285 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_286 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_287 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_288 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_289 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_290 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_291 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_292 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_293 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_294 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_271 + xi_272 + xi_273 + xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_283 + xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_135 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_136 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_137 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_138 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_141 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_142 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_031(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_031_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_032_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 1 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 1 }];
   const double xi_80 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_81 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, 0 }];
   const double xi_82 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_199 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_200 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_201 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_202 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_203 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_204 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_85 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_87 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_88 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_89 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_91 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_94 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_95 = xi_81*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_96 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_176 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_177 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_178 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_179 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_180 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_181 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_182 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_183 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_184 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_185 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_186 = xi_81*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_187 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_188 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_189 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_190 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_192 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_193 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_194 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_195 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_196 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_118 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_110 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_111 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_112 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_113 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_114 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_115 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_116 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_117 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_225 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_226 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_227 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_228 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_230 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_231 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_235 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_236 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_237 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_239 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_240 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_241 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_243 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_244 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_245 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_246 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_247 = xi_81*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_248 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_225 + xi_226 + xi_227 + xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246 + xi_247 + xi_248;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_41 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_44 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_45 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_38 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_48 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_50 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_51 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_52 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_53 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_54 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_55 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_56 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_57 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_58 = xi_81*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_59 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_60 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_61 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_62 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_63 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_64 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_65 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_66 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_67 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_68 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_274 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_276 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_278 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_279 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_283 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_284 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_285 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_286 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_287 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_288 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_289 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_290 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_291 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_292 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_293 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_294 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282 + xi_283 + xi_284 + xi_285;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_135 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_136 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_137 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_138 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_140 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_141 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_142 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_145 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_032(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_032_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_102_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_162 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_163 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_164 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, 0, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   const double xi_199 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_200 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_201 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_202 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_203 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_204 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_82 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_83 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_84 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_85 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_87 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_88 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_89 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_176 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_177 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_178 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_179 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_180 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_181 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_182 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_183 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_184 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_185 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_186 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_187 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_188 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_189 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_190 = xi_167*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_191 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_192 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_193 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_194 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_195 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_196 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_112 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_113 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_114 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_115 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_116 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_117 = xi_167*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_118 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_107 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_108 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_109 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_110 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_111 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_222 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_223 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_224 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_225 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_226 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_227 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_228 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_231 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_236 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_237 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_239 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_240 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_241 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_222 + xi_223 + xi_224 + xi_225 + xi_226 + xi_227 + xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_234 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_41 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_45 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_48 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_50 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_52 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_53 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_54 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_55 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_56 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_57 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_58 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_59 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_60 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_61 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_62 = xi_167*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_63 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_64 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_65 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_66 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_67 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_68 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_271 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_272 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_273 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_274 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_276 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_277 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_278 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_283 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_284 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_285 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_286 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_287 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_288 = xi_167*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_289 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_290 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_291 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_292 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_293 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_294 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_271 + xi_272 + xi_273 + xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_283 + xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_135 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_136 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_137 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_138 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_141 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_142 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_102(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_102_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_103_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_80 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_81 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_82 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   const double xi_199 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_200 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_201 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_202 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_203 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_204 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_85 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_87 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_89 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_76*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_91 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_94 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_95 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_96 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_176 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_177 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_178 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_179 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_180 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_181 = xi_76*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_182 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_183 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_184 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_185 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_186 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_187 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_188 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_189 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_190 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_192 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_193 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_194 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_195 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_196 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_118 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_110 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_111 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_112 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_113 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_114 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_115 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_116 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_117 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_225 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_226 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_227 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_228 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_229 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_231 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_232 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_236 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_237 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_239 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_240 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_241 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_76*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_243 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_244 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_245 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_246 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_247 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_248 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_225 + xi_226 + xi_227 + xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246 + xi_247 + xi_248;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_41 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_43 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_45 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_46 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_48 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_50 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_52 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_53 = xi_76*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_54 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_55 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_56 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_57 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_58 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_59 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_60 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_61 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_62 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_63 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_64 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_65 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_66 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_67 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_68 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_274 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_276 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_278 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_280 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_281 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_283 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_284 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_285 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_286 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_287 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_288 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_289 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_290 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_291 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_292 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_293 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_294 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282 + xi_283 + xi_284 + xi_285;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_135 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_136 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_137 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_138 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_139 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_141 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_142 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_103(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_103_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_120_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_162 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_163 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_164 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_199 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_200 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_201 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_202 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_203 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_204 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, 0 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_82 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_83 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_84 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_85 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_87 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_89 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_176 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_177 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_178 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_179 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_180 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_181 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_182 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_183 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_184 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_185 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_186 = xi_163*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_187 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_188 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_189 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_190 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_192 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_193 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_194 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_195 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_196 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_112 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_113 = xi_163*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_114 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_115 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_116 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_117 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_118 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_107 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_108 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_109 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_110 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_111 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_222 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_223 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_224 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_225 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_226 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_227 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_228 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_231 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_236 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_237 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_239 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_240 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_241 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_222 + xi_223 + xi_224 + xi_225 + xi_226 + xi_227 + xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_234 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_41 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_45 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_48 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_50 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_51 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_52 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_53 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_54 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_55 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_56 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_57 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_58 = xi_163*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_59 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_60 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_61 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_62 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_63 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_64 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_65 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_66 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_67 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_68 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_271 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_272 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_273 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_274 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_276 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_277 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_278 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_283 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_284 = xi_163*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_285 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_286 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_287 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_288 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_289 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_290 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_291 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_292 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_293 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_294 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_271 + xi_272 + xi_273 + xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_283 + xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_135 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_136 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_137 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_138 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_141 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_142 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_120(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_120_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_123_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 1 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, 0 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_80 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_81 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, 0 }];
   const double xi_82 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_199 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_200 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_201 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_202 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 1, -1 }];
   const double xi_203 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_204 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, -1 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, -1 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_85 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_87 = xi_73*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_88 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_89 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_90 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_91 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_94 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_95 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_96 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_176 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_177 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_178 = xi_73*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_179 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_180 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_181 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_182 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_183 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_184 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_185 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_186 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_187 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_188 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_189 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_190 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_191 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_192 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_193 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_194 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_195 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_196 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_118 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_110 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_111 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_112 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_113 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_114 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_115 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_116 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_117 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_225 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_226 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_227 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_228 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_229 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_231 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_236 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_237 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_239 = xi_73*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_240 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_241 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_242 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_243 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_244 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_245 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_246 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_247 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_248 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_225 + xi_226 + xi_227 + xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246 + xi_247 + xi_248;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_41 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_43 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_45 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_48 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_50 = xi_73*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_51 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_52 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_53 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_54 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_55 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_56 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_57 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_58 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_59 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_60 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_61 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_62 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_63 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_64 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_65 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_66 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_67 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_68 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_274 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_276 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_278 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_280 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_283 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_284 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_285 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_286 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_287 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_288 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_289 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_290 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_291 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_292 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_293 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_294 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282 + xi_283 + xi_284 + xi_285;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_135 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_136 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_137 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_138 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_139 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_141 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_142 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_123(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_123_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_130_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_80 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_81 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_82 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_174 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_175 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_176 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_211 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_212 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_213 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_85 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_87 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_74*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_89 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_90 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_91 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_94 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_95 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_96 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_179 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_180 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_181 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_182 = xi_74*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_183 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_184 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_185 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_186 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_187 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_188 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_189 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_190 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_192 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_193 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_194 = xi_168*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_195 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_196 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_197 = xi_171*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_198 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_199 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_200 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_201 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_202 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196 + xi_197 + xi_198 + xi_199 + xi_200 + xi_201 + xi_202;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_124 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_113 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_114 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_115 = xi_168*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_116 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_117 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_118 = xi_171*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_119 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_120 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_121 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_122 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_123 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122 + xi_123 + xi_124;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_228 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_231 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_233 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_236 = xi_213*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_237 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_239 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_240 = xi_74*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_241 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_242 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_243 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_244 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_245 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_246 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_247 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_248 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246 + xi_247 + xi_248;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_39 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_41 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_213*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_45 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_48 = xi_74*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_49 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_50 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_51 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_52 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_53 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_54 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_55 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_56 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_57 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_58 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_59 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_60 = xi_168*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_61 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_63 = xi_171*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_64 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_65 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_66 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_67 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_68 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_274 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_276 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_277 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_278 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_279 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_213*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_283 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_284 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_285 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_286 = xi_168*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_287 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_288 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_289 = xi_171*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_290 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_291 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_292 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_293 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_294 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_283 + xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_138 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_141 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_142 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_143 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_213*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_130(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_130_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_132_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 1 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, 0 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_80 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_81 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, 0 }];
   const double xi_82 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 1, -1 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, -1 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, -1 }];
   const double xi_174 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_175 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_176 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_211 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_212 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_213 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_85 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_87 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_88 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_89 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_90 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_91 = xi_77*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_94 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_95 = xi_81*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_96 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_179 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_180 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_181 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_182 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_183 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_184 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_185 = xi_77*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_186 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_187 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_188 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_189 = xi_81*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_190 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_192 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_193 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = xi_168*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_195 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_196 = xi_170*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_197 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_198 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_199 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_200 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_201 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_202 = xi_176*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196 + xi_197 + xi_198 + xi_199 + xi_200 + xi_201 + xi_202;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_124 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_113 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_114 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_115 = xi_168*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_116 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_117 = xi_170*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_118 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_119 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_120 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_121 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_122 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_123 = xi_176*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122 + xi_123 + xi_124;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_228 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_230 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_231 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_232 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_211*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_236 = xi_213*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_237 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_239 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_240 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_241 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_242 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_243 = xi_77*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_244 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_245 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_246 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_247 = xi_81*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_248 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246 + xi_247 + xi_248;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_38 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_40 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_41 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_211*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_213*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_45 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_48 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_50 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_51 = xi_77*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_52 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_53 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_54 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_55 = xi_81*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_56 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_57 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_58 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_59 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = xi_168*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_61 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_62 = xi_170*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_63 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_64 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_65 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_66 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_67 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_68 = xi_176*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_274 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_276 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_278 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_211*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_213*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_283 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_284 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_285 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_286 = xi_168*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_287 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_288 = xi_170*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_289 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_290 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_291 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_292 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_293 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_294 = xi_176*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_283 + xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_138 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_140 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_141 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_142 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_211*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_213*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_132(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_132_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_201_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_80 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_81 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_82 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_174 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, 0, 0 }];
   const double xi_175 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_176 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_211 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_212 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_213 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_85 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_87 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_89 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_91 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_94 = xi_80*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_95 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_96 = xi_82*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_179 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_180 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_181 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_182 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_183 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_184 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_185 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_186 = xi_78*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_187 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_188 = xi_80*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_189 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_190 = xi_82*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_192 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_193 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_194 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_195 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_196 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_197 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_198 = xi_172*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_199 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_200 = xi_174*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_201 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_202 = xi_176*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196 + xi_197 + xi_198 + xi_199 + xi_200 + xi_201 + xi_202;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_124 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_113 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_114 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_115 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_116 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_117 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_118 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_119 = xi_172*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_120 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_121 = xi_174*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_122 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_123 = xi_176*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122 + xi_123 + xi_124;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_228 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_230 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_231 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_236 = xi_213*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_237 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_239 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_240 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_241 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_243 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_244 = xi_78*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_245 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_246 = xi_80*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_247 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_248 = xi_82*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246 + xi_247 + xi_248;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_38 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_41 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_44 = xi_213*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_45 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_47 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_48 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_50 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_52 = xi_78*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_53 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_54 = xi_80*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_55 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_56 = xi_82*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_57 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_58 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_59 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_60 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_61 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_62 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_63 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_64 = xi_172*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_65 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_66 = xi_174*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_67 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_68 = xi_176*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_274 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_276 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_278 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_282 = xi_213*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_283 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_284 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_285 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_286 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_287 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_288 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_289 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_290 = xi_172*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_291 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_292 = xi_174*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_293 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_294 = xi_176*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_283 + xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_138 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_140 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_141 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_142 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_146 = xi_213*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_201(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_201_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_203_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_80 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_81 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_82 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 1 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 1 }];
   const double xi_174 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_175 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, 0 }];
   const double xi_176 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_211 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_212 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_213 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_85 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_87 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_89 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_91 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_94 = xi_80*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_95 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_96 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_179 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_180 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_181 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_182 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_183 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_184 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_185 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_186 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_187 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_188 = xi_80*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_189 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_190 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_192 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_193 = xi_167*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_194 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_195 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_196 = xi_170*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_197 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_198 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_199 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_200 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_201 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_202 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196 + xi_197 + xi_198 + xi_199 + xi_200 + xi_201 + xi_202;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_124 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_113 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_114 = xi_167*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_115 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_116 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_117 = xi_170*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_118 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_119 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_120 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_121 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_122 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_123 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122 + xi_123 + xi_124;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_228 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_231 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_233 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_235 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_236 = xi_213*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_237 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_239 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_240 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_241 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_243 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_244 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_245 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_246 = xi_80*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_247 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_248 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246 + xi_247 + xi_248;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_41 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_43 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_213*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_45 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_48 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_50 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_52 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_53 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_54 = xi_80*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_55 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_56 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_57 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_58 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_59 = xi_167*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_60 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_61 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_62 = xi_170*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_63 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_64 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_65 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_66 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_67 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_68 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_274 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_276 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_278 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_279 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_281 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_213*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_283 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_284 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_285 = xi_167*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_286 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_287 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_288 = xi_170*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_289 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_290 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_291 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_292 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_293 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_294 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_283 + xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_138 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_141 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_142 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_143 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_145 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_213*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_203(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_203_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_210_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_80 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_81 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_82 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_199 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_200 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_201 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_202 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_203 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_204 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, 0 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_85 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_87 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_89 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_91 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_94 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_95 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_96 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_176 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_177 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_178 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_179 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_180 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_181 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_182 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_183 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_184 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_185 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_186 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_187 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_188 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_189 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_190 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_191 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_192 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_193 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_194 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_195 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_196 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_118 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_110 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_111 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_112 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_113 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_114 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_115 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_116 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_117 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_225 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_226 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_227 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_228 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_231 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_233 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_235 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_236 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_237 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_239 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_240 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_241 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_243 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_244 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_245 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_246 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_247 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_248 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_225 + xi_226 + xi_227 + xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246 + xi_247 + xi_248;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_41 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_45 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_47 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_38 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_48 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_50 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_52 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_53 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_54 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_55 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_56 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_57 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_58 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_59 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_60 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_61 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_62 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_63 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_64 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_65 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_66 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_67 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_68 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_274 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_276 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_278 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_282 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_283 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_284 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_285 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_286 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_287 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_288 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_289 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_290 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_291 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_292 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_293 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_294 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282 + xi_283 + xi_284 + xi_285;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_135 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_136 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_137 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_138 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_141 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_142 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_143 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_145 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_210(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_210_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_213_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_162 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_163 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_164 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 1 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_199 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_200 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_201 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_202 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 1, -1 }];
   const double xi_203 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_204 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, -1 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, -1 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_82 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_83 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_84 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_85 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_87 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_89 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_176 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_177 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_178 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_179 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_180 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_181 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_182 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_183 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_184 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_185 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_186 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_187 = xi_164*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_188 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_189 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_190 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_192 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_193 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_194 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_195 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_196 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_112 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_113 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_114 = xi_164*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_115 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_116 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_117 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_118 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_107 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_108 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_109 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_110 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_111 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_222 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_223 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_224 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_225 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_226 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_227 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_228 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_231 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_236 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_237 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_239 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_240 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_241 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_222 + xi_223 + xi_224 + xi_225 + xi_226 + xi_227 + xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_234 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_41 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_45 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_37 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_48 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_50 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_51 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_52 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_53 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_54 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_55 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_56 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_57 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_58 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_59 = xi_164*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_60 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_61 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_62 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_63 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_64 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_65 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_66 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_67 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_68 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_271 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_272 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_273 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_274 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_276 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_278 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_280 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_283 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_284 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_285 = xi_164*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_286 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_287 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_288 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_289 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_290 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_291 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_292 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_293 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_294 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_271 + xi_272 + xi_273 + xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_283 + xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_135 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_136 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_137 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_138 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_141 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_142 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_144 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_213(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_213_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_230_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_80 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_81 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_82 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_199 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_200 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_201 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 1 }];
   const double xi_202 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_203 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_204 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 1 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_85 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_87 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_89 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_90 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_91 = xi_77*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_94 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_95 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_96 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_176 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_177 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_178 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_179 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_180 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_181 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_182 = xi_77*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_183 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_184 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_185 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_186 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_187 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_188 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_189 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_190 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_192 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_193 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_195 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_196 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_118 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_110 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_111 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_112 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_113 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_114 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_115 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_116 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_117 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_225 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_226 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_227 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_228 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_231 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_232 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_234 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_236 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_237 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_239 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_240 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_241 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_242 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_243 = xi_77*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_244 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_245 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_246 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_247 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_248 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_225 + xi_226 + xi_227 + xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246 + xi_247 + xi_248;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_41 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_42 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_45 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_46 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_37 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_48 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_50 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_52 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_53 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_54 = xi_77*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_55 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_56 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_57 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_58 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_59 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_60 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_61 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_62 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_63 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_64 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_65 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_66 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_67 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_68 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_274 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_276 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_277 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_278 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_280 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_281 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_283 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_284 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_285 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_286 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_287 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_288 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_289 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_290 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_291 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_292 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_293 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_294 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282 + xi_283 + xi_284 + xi_285;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_135 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_136 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_137 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_138 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_141 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_142 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_144 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_230(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_230_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_231_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_162 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_163 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_164 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 1, -1 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, -1 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, -1 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_199 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_200 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_201 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 1 }];
   const double xi_202 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_203 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, 0 }];
   const double xi_204 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_82 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_83 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_84 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_85 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_87 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_88 = xi_77*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_89 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_176 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_177 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_178 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_179 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_180 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_181 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_182 = xi_77*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_183 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_184 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_185 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_186 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_187 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_188 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_189 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_190 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_191 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_192 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_193 = xi_170*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_194 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_195 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_196 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_112 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_113 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_114 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_115 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_116 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_117 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_118 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_107 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_108 = xi_170*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_109 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_110 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_111 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_222 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_223 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_224 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_225 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_226 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_227 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_228 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_231 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_236 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_237 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_239 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_240 = xi_77*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_241 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_222 + xi_223 + xi_224 + xi_225 + xi_226 + xi_227 + xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_234 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_41 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_42 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_45 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_37 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_48 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_50 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_52 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_53 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_54 = xi_77*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_55 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_56 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_57 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_58 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_59 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_60 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_61 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_62 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_63 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_64 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_65 = xi_170*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_66 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_67 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_68 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_271 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_272 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_273 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_274 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_276 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_278 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_280 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_283 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_284 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_285 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_286 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_287 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_288 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_289 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_290 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_291 = xi_170*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_292 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_293 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_294 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_271 + xi_272 + xi_273 + xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_283 + xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_135 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_136 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_137 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_138 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_141 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_142 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_144 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_231(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_231_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_301_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_80 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_81 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_82 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   const double xi_199 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_200 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_201 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_202 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_203 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_204 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_85 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_87 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_89 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_91 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_94 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_95 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_96 = xi_82*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_176 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_177 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_178 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_179 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_180 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_181 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_182 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_183 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_184 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_185 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_186 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_187 = xi_82*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_188 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_189 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_190 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_192 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_193 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_194 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_195 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_196 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_118 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_110 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_111 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_112 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_113 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_114 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_115 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_116 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_117 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_225 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_226 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_227 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_228 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_229 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_231 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_234 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_236 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_237 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_239 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_240 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_241 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_243 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_244 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_245 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_246 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_247 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_248 = xi_82*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_225 + xi_226 + xi_227 + xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246 + xi_247 + xi_248;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_41 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_43 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_45 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_37 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_48 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_50 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_52 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_53 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_54 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_55 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_56 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_57 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_58 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_59 = xi_82*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_60 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_61 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_62 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_63 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_64 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_65 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_66 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_67 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_68 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_274 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_276 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_278 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_283 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_284 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_285 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_286 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_287 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_288 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_289 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_290 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_291 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_292 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_293 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_294 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282 + xi_283 + xi_284 + xi_285;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_135 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_136 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_137 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_138 = xi_202*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_139 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_141 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_142 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_144 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_301(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_301_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_302_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_162 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_163 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_164 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 1 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 1 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   const double xi_199 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_200 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_201 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_202 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_203 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_204 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_82 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_83 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_84 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_85 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_87 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_89 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_90 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_176 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_177 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_178 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_179 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_180 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_181 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_182 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_183 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_184 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_185 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_186 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_187 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_188 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_189 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_190 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_191 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_192 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_193 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_194 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_195 = xi_172*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_196 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_112 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_113 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_114 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_115 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_116 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_117 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_118 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_107 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_108 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_109 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_110 = xi_172*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_111 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_222 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_223 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_224 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_225 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_226 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_227 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_228 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_229 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_231 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_234 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_236 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_237 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_239 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_240 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_241 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_242 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_222 + xi_223 + xi_224 + xi_225 + xi_226 + xi_227 + xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_234 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_41 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_45 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_46 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_48 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_50 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_52 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_53 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_54 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_55 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_56 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_57 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_58 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_59 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_60 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_61 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_62 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_63 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_64 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_65 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_66 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_67 = xi_172*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_68 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_271 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_272 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_273 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_274 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_276 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_278 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_283 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_284 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_285 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_286 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_287 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_288 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_289 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_290 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_291 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_292 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_293 = xi_172*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_294 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_271 + xi_272 + xi_273 + xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_283 + xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_135 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_136 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_137 = xi_201*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_138 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_141 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_142 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_208*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_302(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_302_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_310_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_80 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_81 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_82 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::X][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_174 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_175 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_176 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_211 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_212 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_213 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_85 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_87 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_74*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_89 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_90 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_91 = xi_77*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_94 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_95 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_96 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_179 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_180 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_181 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_182 = xi_74*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_183 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_184 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_185 = xi_77*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_186 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_187 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_188 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_189 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_190 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_191 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_192 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_193 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_194 = xi_168*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_195 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_196 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_197 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_198 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_199 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_200 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_201 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_202 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196 + xi_197 + xi_198 + xi_199 + xi_200 + xi_201 + xi_202;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_124 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_113 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_114 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_115 = xi_168*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_116 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_117 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_118 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_119 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_120 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_121 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_122 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_123 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122 + xi_123 + xi_124;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_228 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_231 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_236 = xi_213*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_237 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_239 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_240 = xi_74*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_241 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_242 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_243 = xi_77*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_244 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_245 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_246 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_247 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_248 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246 + xi_247 + xi_248;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_39 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_41 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_213*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_45 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_72*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_48 = xi_74*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_49 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_50 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_77*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_52 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_53 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_54 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_55 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_56 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_57 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_58 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_59 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_60 = xi_168*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_61 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_63 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_64 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_65 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_66 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_67 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_68 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_274 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_276 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_277 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_278 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_213*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_283 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_284 = xi_166*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_285 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_286 = xi_168*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_287 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_288 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_289 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_290 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_291 = xi_173*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_292 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_293 = xi_175*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_294 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_283 + xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_138 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_141 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_142 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_211*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_213*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_310(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_310_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_312_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 1, -1 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, -1 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, -1 }];
   const double xi_80 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_81 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_82 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 1 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_174 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_175 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, 0 }];
   const double xi_176 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_211 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_212 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_213 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_85 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_87 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_88 = xi_74*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_89 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_76*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_91 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_94 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_95 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_96 = xi_82*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_179 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_180 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_181 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_182 = xi_74*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_183 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_184 = xi_76*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_185 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_186 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_187 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_188 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_189 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_190 = xi_82*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_192 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_193 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_194 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_195 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_196 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_197 = xi_171*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_198 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_199 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_200 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_201 = xi_175*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_202 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196 + xi_197 + xi_198 + xi_199 + xi_200 + xi_201 + xi_202;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_124 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_113 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_114 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_115 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_116 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_117 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_118 = xi_171*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_119 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_120 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_121 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_122 = xi_175*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_123 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122 + xi_123 + xi_124;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_228 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_231 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_232 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_211*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_235 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_236 = xi_213*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_237 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_239 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_240 = xi_74*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_241 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_76*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_243 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_244 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_245 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_246 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_247 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_248 = xi_82*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246 + xi_247 + xi_248;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_39 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_40 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_41 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_211*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_43 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_213*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_45 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_46 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_48 = xi_74*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_49 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_50 = xi_76*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_52 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_53 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_54 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_55 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_56 = xi_82*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_57 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_58 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_59 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_60 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_61 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_62 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_63 = xi_171*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_64 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_65 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_66 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_67 = xi_175*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_68 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_274 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_276 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_278 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_211*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_281 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_213*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_283 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_284 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_285 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_286 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_287 = xi_169*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_288 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_289 = xi_171*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_290 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_291 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_292 = xi_174*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_293 = xi_175*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_294 = xi_176*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_283 + xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_138 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_206*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_141 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_142 = xi_209*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_211*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_145 = xi_212*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_146 = xi_213*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_312(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_312_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_320_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Y][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_162 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_163 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_164 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::Z][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_199 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, 0 }];
   const double xi_200 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_201 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 1 }];
   const double xi_202 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_203 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_204 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 1 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_82 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_83 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_84 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_85 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_87 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_88 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_89 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_176 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_177 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_178 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_179 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_180 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_181 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_182 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_183 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_184 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_185 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_186 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_187 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_188 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_189 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_190 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_168*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_192 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_193 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_194 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_195 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_196 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_112 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_113 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_114 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_115 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_116 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_117 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_118 = xi_168*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_107 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_108 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_109 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_110 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_111 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_222 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_223 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_224 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_225 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_226 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_227 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_228 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_229 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_231 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_233 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_236 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_237 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_239 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_240 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_241 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_222 + xi_223 + xi_224 + xi_225 + xi_226 + xi_227 + xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_234 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_41 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_44 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_45 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_46 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_37 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_39 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_48 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_50 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_52 = xi_75*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_53 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_54 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_55 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_56 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_57 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_58 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_59 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_60 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_61 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_62 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_63 = xi_168*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_64 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_65 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_66 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_67 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_68 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_271 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_272 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_273 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_274 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_276 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_278 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_282 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_283 = xi_162*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_284 = xi_163*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_285 = xi_164*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_286 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_287 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_288 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_289 = xi_168*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_290 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_291 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_292 = xi_171*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_293 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_294 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_271 + xi_272 + xi_273 + xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_283 + xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_135 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_136 = xi_200*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_137 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_138 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_203*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_140 = xi_204*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_141 = xi_205*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_142 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_207*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_144 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_146 = xi_210*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_320(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_320_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    
static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_321_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
   const double xi_71 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_72 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_73 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 1, -1 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, -1 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, -1 }];
   const double xi_80 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_81 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_82 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_165 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_166 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_167 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_170 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_199 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_200 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_201 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 1 }];
   const double xi_202 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_203 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, 0 }];
   const double xi_204 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_205 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_206 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_207 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_208 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_209 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, 0 }];
   const double xi_210 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_85 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_87 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_89 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_91 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_94 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_95 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_96 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_176 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_177 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_178 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_179 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_180 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_181 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_182 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_183 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_184 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_185 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_186 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_187 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_188 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_189 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_190 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_192 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_193 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = xi_171*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_195 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_196 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_118 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_110 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_111 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_112 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_113 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_114 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_115 = xi_171*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_116 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_117 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_225 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_226 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_227 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_228 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_229 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_230 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_231 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_232 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_233 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_234 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_235 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_236 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_237 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_239 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_240 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_241 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_243 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_244 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_245 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_246 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_247 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_248 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_225 + xi_226 + xi_227 + xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242 + xi_243 + xi_244 + xi_245 + xi_246 + xi_247 + xi_248;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_41 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_44 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_45 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_46 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_37 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_39 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_48 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_50 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_52 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_53 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_54 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_55 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_56 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_57 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_58 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_59 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_60 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_61 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_62 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_63 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_64 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_65 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_66 = xi_171*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_67 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_68 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_274 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_276 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_278 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_279 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_280 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_281 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_282 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_283 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_284 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_285 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_286 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_287 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_288 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_289 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_290 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_291 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_292 = xi_171*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_293 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_294 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282 + xi_283 + xi_284 + xi_285;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_135 = xi_199*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_136 = xi_200*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_137 = xi_201*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_138 = xi_202*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_203*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_140 = xi_204*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_141 = xi_205*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_142 = xi_206*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_207*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_144 = xi_208*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_209*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_146 = xi_210*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_321(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > > e2e_cell_stencil, int level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_321_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace edgedof
} // namespace hyteg
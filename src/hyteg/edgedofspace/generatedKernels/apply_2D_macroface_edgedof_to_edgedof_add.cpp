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

#include "apply_2D_macroface_edgedof_to_edgedof_add.hpp"

namespace hyteg {
namespace edgedof {
namespace macroface {
namespace generated {

static void apply_2D_macroface_edgedof_to_edgedof_add_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * RESTRICT const _data_edgeToDiagonalEdgeFaceStencil, double const * RESTRICT const _data_edgeToHorizontalEdgeFaceStencil, double const * RESTRICT const _data_edgeToVerticalEdgeFaceStencil, int level)
{
   const double xi_37 = _data_edgeToDiagonalEdgeFaceStencil[0];
   const double xi_38 = _data_edgeToDiagonalEdgeFaceStencil[1];
   const double xi_39 = _data_edgeToDiagonalEdgeFaceStencil[3];
   const double xi_40 = _data_edgeToDiagonalEdgeFaceStencil[4];
   const double xi_41 = _data_edgeToDiagonalEdgeFaceStencil[2];
   const double xi_91 = _data_edgeToVerticalEdgeFaceStencil[4];
   const double xi_92 = _data_edgeToVerticalEdgeFaceStencil[2];
   const double xi_93 = _data_edgeToVerticalEdgeFaceStencil[3];
   const double xi_94 = _data_edgeToVerticalEdgeFaceStencil[1];
   const double xi_95 = _data_edgeToVerticalEdgeFaceStencil[0];
   const double xi_112 = _data_edgeToHorizontalEdgeFaceStencil[1];
   const double xi_113 = _data_edgeToHorizontalEdgeFaceStencil[3];
   const double xi_114 = _data_edgeToHorizontalEdgeFaceStencil[0];
   const double xi_115 = _data_edgeToHorizontalEdgeFaceStencil[4];
   const double xi_116 = _data_edgeToHorizontalEdgeFaceStencil[2];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_50 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_45 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_46 = xi_38*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_47 = xi_39*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_48 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_49 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_104 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_110 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_99 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_100 = xi_38*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_101 = xi_39*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_102 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_103 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_105 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_106 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_107 = xi_93*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_108 = xi_94*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_109 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_101 + xi_102 + xi_103 + xi_104 + xi_99;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_105 + xi_106 + xi_107 + xi_108 + xi_109 + xi_110;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_65 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_61 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_93*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_63 = xi_94*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_64 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_130 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_136 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_125 = xi_112*_data_edgeFaceSrc_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
            const double xi_126 = xi_113*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_127 = xi_114*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_128 = xi_115*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_129 = xi_116*_data_edgeFaceSrc_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_131 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_132 = xi_38*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_133 = xi_39*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_134 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_135 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_125 + xi_126 + xi_127 + xi_128 + xi_129 + xi_130;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_136;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_23 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_29 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_35 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_18 = xi_112*_data_edgeFaceSrc_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
            const double xi_19 = xi_113*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = xi_114*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_21 = xi_115*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_22 = xi_116*_data_edgeFaceSrc_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_24 = xi_37*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_25 = xi_38*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_26 = xi_39*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_27 = xi_40*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_28 = xi_41*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_30 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_31 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_32 = xi_93*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_33 = xi_94*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_34 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_156 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_162 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_151 = xi_112*_data_edgeFaceSrc_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
            const double xi_152 = xi_113*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_153 = xi_114*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_154 = xi_115*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_155 = xi_116*_data_edgeFaceSrc_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_157 = xi_91*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_158 = xi_92*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_159 = xi_93*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_160 = xi_94*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_161 = xi_95*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_151 + xi_152 + xi_153 + xi_154 + xi_155 + xi_156;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_157 + xi_158 + xi_159 + xi_160 + xi_161 + xi_162;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_80 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_75 = xi_112*_data_edgeFaceSrc_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
            const double xi_76 = xi_113*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_77 = xi_114*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_78 = xi_115*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_79 = xi_116*_data_edgeFaceSrc_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80;
         }
      }
   }
}


void apply_2D_macroface_edgedof_to_edgedof_add(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * RESTRICT const _data_edgeToDiagonalEdgeFaceStencil, double const * RESTRICT const _data_edgeToHorizontalEdgeFaceStencil, double const * RESTRICT const _data_edgeToVerticalEdgeFaceStencil, int level)
{
    switch( level )
    {

    default:
        apply_2D_macroface_edgedof_to_edgedof_add_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeToDiagonalEdgeFaceStencil, _data_edgeToHorizontalEdgeFaceStencil, _data_edgeToVerticalEdgeFaceStencil, level);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace edgedof
} // namespace hyteg
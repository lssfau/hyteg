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

#include "sor_3D_macrocell_P2_update_vertexdofs.hpp"

namespace hyteg {
namespace P2 {
namespace macrocell {
namespace generated {

static void sor_3D_macrocell_P2_update_vertexdofs_level_any(double const * RESTRICT const _data_edgeCellDst_X, double const * RESTRICT const _data_edgeCellDst_XY, double const * RESTRICT const _data_edgeCellDst_XYZ, double const * RESTRICT const _data_edgeCellDst_XZ, double const * RESTRICT const _data_edgeCellDst_Y, double const * RESTRICT const _data_edgeCellDst_YZ, double const * RESTRICT const _data_edgeCellDst_Z, double * RESTRICT _data_vertexCellDst, double const * RESTRICT const _data_vertexCellRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2vStencilMap, int level, double relax, std::map< hyteg::indexing::Index, double > v2vStencilMap)
{
   const double xi_69 = 1.0;
   const double xi_70 = -relax;
   const double xi_1 = v2vStencilMap[{ 0, 0, 0 }];
   const double xi_67 = 1 / (xi_1);
   const double xi_2 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_3 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_4 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_5 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_6 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_7 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_8 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_9 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_10 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_11 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_12 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_13 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_14 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_15 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_16 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_17 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_18 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_19 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_20 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_21 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_22 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_23 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_24 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_25 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_26 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_27 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_28 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_29 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_30 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_31 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_32 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_33 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_34 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_35 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_36 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_37 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_38 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_39 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_40 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_41 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_42 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_43 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_44 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_45 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_46 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_47 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_48 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_49 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_50 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_51 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_52 = v2vStencilMap[{ -1, 0, 0 }];
   const double xi_53 = v2vStencilMap[{ -1, 0, 1 }];
   const double xi_54 = v2vStencilMap[{ -1, 1, -1 }];
   const double xi_55 = v2vStencilMap[{ -1, 1, 0 }];
   const double xi_56 = v2vStencilMap[{ 0, -1, 0 }];
   const double xi_57 = v2vStencilMap[{ 0, -1, 1 }];
   const double xi_58 = v2vStencilMap[{ 0, 0, -1 }];
   const double xi_59 = v2vStencilMap[{ 0, 0, 1 }];
   const double xi_60 = v2vStencilMap[{ 0, 1, -1 }];
   const double xi_61 = v2vStencilMap[{ 0, 1, 0 }];
   const double xi_62 = v2vStencilMap[{ 1, -1, 0 }];
   const double xi_63 = v2vStencilMap[{ 1, -1, 1 }];
   const double xi_64 = v2vStencilMap[{ 1, 0, -1 }];
   const double xi_65 = v2vStencilMap[{ 1, 0, 0 }];
   for (int ctr_3 = 1; ctr_3 < (1 << (level)); ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + (1 << (level)); ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_135 = _data_vertexCellRhs[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_81 = -xi_2*_data_edgeCellDst_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6)) - 1];
            const double xi_92 = -xi_3*_data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6)) - 1];
            const double xi_103 = -xi_4*_data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6)) - 1];
            const double xi_114 = -xi_5*_data_edgeCellDst_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
            const double xi_125 = -xi_6*_data_edgeCellDst_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
            const double xi_132 = -xi_7*_data_edgeCellDst_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
            const double xi_133 = -xi_8*_data_edgeCellDst_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_134 = -xi_9*_data_edgeCellDst_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_71 = -xi_10*_data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_72 = -xi_11*_data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_73 = -xi_12*_data_edgeCellDst_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_74 = -xi_13*_data_edgeCellDst_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_75 = -xi_14*_data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_76 = -xi_15*_data_edgeCellDst_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_77 = -xi_16*_data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_78 = -xi_17*_data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_79 = -xi_18*_data_edgeCellDst_XZ[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_80 = -xi_19*_data_edgeCellDst_XZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_82 = -xi_20*_data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_83 = -xi_21*_data_edgeCellDst_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_84 = -xi_22*_data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_85 = -xi_23*_data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_86 = -xi_24*_data_edgeCellDst_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_87 = -xi_25*_data_edgeCellDst_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_88 = -xi_26*_data_edgeCellDst_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_89 = -xi_27*_data_edgeCellDst_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_90 = -xi_28*_data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_91 = -xi_29*_data_edgeCellDst_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_93 = -xi_30*_data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_94 = -xi_31*_data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_95 = -xi_32*_data_edgeCellDst_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_96 = -xi_33*_data_edgeCellDst_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_97 = -xi_34*_data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_98 = -xi_35*_data_edgeCellDst_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_99 = -xi_36*_data_edgeCellDst_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) + 1];
            const double xi_100 = -xi_37*_data_edgeCellDst_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            const double xi_101 = -xi_38*_data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_102 = -xi_39*_data_edgeCellDst_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_104 = -xi_40*_data_edgeCellDst_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_105 = -xi_41*_data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_106 = -xi_42*_data_edgeCellDst_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_107 = -xi_43*_data_edgeCellDst_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            const double xi_108 = -xi_44*_data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_109 = -xi_45*_data_edgeCellDst_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_110 = -xi_46*_data_edgeCellDst_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_111 = -xi_47*_data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_112 = -xi_48*_data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_113 = -xi_49*_data_edgeCellDst_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_115 = -xi_50*_data_edgeCellDst_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            const double xi_116 = -xi_51*_data_edgeCellDst_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) + 1];
            const double xi_117 = -xi_52*_data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_118 = -xi_53*_data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - 1];
            const double xi_119 = -xi_54*_data_vertexCellDst[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 3) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)*(-ctr_3 + (1 << (level)) + 4)) / (6)) - 1];
            const double xi_120 = -xi_55*_data_vertexCellDst[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_121 = -xi_56*_data_vertexCellDst[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_122 = -xi_57*_data_vertexCellDst[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6))];
            const double xi_123 = -xi_58*_data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)*(-ctr_3 + (1 << (level)) + 4)) / (6))];
            const double xi_124 = -xi_59*_data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6))];
            const double xi_126 = -xi_60*_data_vertexCellDst[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 3) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)*(-ctr_3 + (1 << (level)) + 4)) / (6))];
            const double xi_127 = -xi_61*_data_vertexCellDst[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_128 = -xi_62*_data_vertexCellDst[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) + 1];
            const double xi_129 = -xi_63*_data_vertexCellDst[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) + 1];
            const double xi_130 = -xi_64*_data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)*(-ctr_3 + (1 << (level)) + 4)) / (6)) + 1];
            const double xi_131 = -xi_65*_data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) + 1];
            _data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))] = relax*xi_67*(xi_100 + xi_101 + xi_102 + xi_103 + xi_104 + xi_105 + xi_106 + xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122 + xi_123 + xi_124 + xi_125 + xi_126 + xi_127 + xi_128 + xi_129 + xi_130 + xi_131 + xi_132 + xi_133 + xi_134 + xi_135 + xi_71 + xi_72 + xi_73 + xi_74 + xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80 + xi_81 + xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96 + xi_97 + xi_98 + xi_99) + (xi_69 + xi_70)*_data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
         }
      }
   }
}


void sor_3D_macrocell_P2_update_vertexdofs(double const * RESTRICT const _data_edgeCellDst_X, double const * RESTRICT const _data_edgeCellDst_XY, double const * RESTRICT const _data_edgeCellDst_XYZ, double const * RESTRICT const _data_edgeCellDst_XZ, double const * RESTRICT const _data_edgeCellDst_Y, double const * RESTRICT const _data_edgeCellDst_YZ, double const * RESTRICT const _data_edgeCellDst_Z, double * RESTRICT _data_vertexCellDst, double const * RESTRICT const _data_vertexCellRhs, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2vStencilMap, int level, double relax, std::map< hyteg::indexing::Index, double > v2vStencilMap)
{
    switch( level )
    {

    default:
        sor_3D_macrocell_P2_update_vertexdofs_level_any(_data_edgeCellDst_X, _data_edgeCellDst_XY, _data_edgeCellDst_XYZ, _data_edgeCellDst_XZ, _data_edgeCellDst_Y, _data_edgeCellDst_YZ, _data_edgeCellDst_Z, _data_vertexCellDst, _data_vertexCellRhs, e2vStencilMap, level, relax, v2vStencilMap);
        break;
    }
}
    

} // namespace generated
} // namespace macrocell
} // namespace P2
} // namespace hyteg
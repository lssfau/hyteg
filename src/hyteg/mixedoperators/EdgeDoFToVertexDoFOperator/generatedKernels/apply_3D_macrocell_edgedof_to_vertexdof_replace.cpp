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

#include "apply_3D_macrocell_edgedof_to_vertexdof_replace.hpp"

namespace hyteg {
namespace EdgeDoFToVertexDoF {
namespace generated {

static void apply_3D_macrocell_edgedof_to_vertexdof_replace_level_any(double const * RESTRICT const _data_edgeCellSrc_X, double const * RESTRICT const _data_edgeCellSrc_XY, double const * RESTRICT const _data_edgeCellSrc_XYZ, double const * RESTRICT const _data_edgeCellSrc_XZ, double const * RESTRICT const _data_edgeCellSrc_Y, double const * RESTRICT const _data_edgeCellSrc_YZ, double const * RESTRICT const _data_edgeCellSrc_Z, double * RESTRICT _data_vertexCellDst, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2vStencilMap, int level)
{
   const double xi_1 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_2 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_3 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_5 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_6 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_7 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_8 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_9 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_10 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_11 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_12 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_13 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_14 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_15 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_16 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_17 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_18 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_19 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_20 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_21 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_22 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_23 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_24 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_25 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_26 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_27 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_28 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_29 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_30 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_31 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_32 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_33 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_34 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_35 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_36 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_37 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_38 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_39 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_40 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_41 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_42 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_43 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_44 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_45 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_46 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_47 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_48 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_49 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_50 = e2vStencilMap[hyteg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   for (int ctr_3 = 1; ctr_3 < (1 << (level)); ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + (1 << (level)); ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_53 = xi_1*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6)) - 1];
            const double xi_64 = xi_2*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6)) - 1];
            const double xi_75 = xi_3*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6)) - 1];
            const double xi_86 = xi_4*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
            const double xi_97 = xi_5*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
            const double xi_99 = xi_6*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
            const double xi_100 = xi_7*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_101 = xi_8*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_102 = xi_9*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_54 = xi_10*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_55 = xi_11*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_56 = xi_12*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_57 = xi_13*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_58 = xi_14*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_59 = xi_15*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_60 = xi_16*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_61 = xi_17*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_62 = xi_18*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_63 = xi_19*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_65 = xi_20*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_66 = xi_21*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_67 = xi_22*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_68 = xi_23*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_69 = xi_24*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_70 = xi_25*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_71 = xi_26*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_72 = xi_27*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_73 = xi_28*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_74 = xi_29*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_76 = xi_30*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_77 = xi_31*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_78 = xi_32*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_79 = xi_33*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_80 = xi_34*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_81 = xi_35*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) + 1];
            const double xi_82 = xi_36*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            const double xi_83 = xi_37*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_84 = xi_38*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_85 = xi_39*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_87 = xi_40*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_88 = xi_41*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_89 = xi_42*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            const double xi_90 = xi_43*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_91 = xi_44*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_92 = xi_45*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_93 = xi_46*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_94 = xi_47*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_95 = xi_48*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_96 = xi_49*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            const double xi_98 = xi_50*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) + 1];
            _data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))] = xi_100 + xi_101 + xi_102 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68 + xi_69 + xi_70 + xi_71 + xi_72 + xi_73 + xi_74 + xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80 + xi_81 + xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96 + xi_97 + xi_98 + xi_99;
         }
      }
   }
}


void apply_3D_macrocell_edgedof_to_vertexdof_replace(double const * RESTRICT const _data_edgeCellSrc_X, double const * RESTRICT const _data_edgeCellSrc_XY, double const * RESTRICT const _data_edgeCellSrc_XYZ, double const * RESTRICT const _data_edgeCellSrc_XZ, double const * RESTRICT const _data_edgeCellSrc_Y, double const * RESTRICT const _data_edgeCellSrc_YZ, double const * RESTRICT const _data_edgeCellSrc_Z, double * RESTRICT _data_vertexCellDst, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::Index, double > > e2vStencilMap, int level)
{
    switch( level )
    {

    default:
        apply_3D_macrocell_edgedof_to_vertexdof_replace_level_any(_data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, _data_vertexCellDst, e2vStencilMap, level);
        break;
    }
}
    

} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hyteg
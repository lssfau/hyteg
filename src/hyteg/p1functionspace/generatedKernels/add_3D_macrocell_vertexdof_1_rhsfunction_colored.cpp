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

#include "add_3D_macrocell_vertexdof_1_rhsfunction_colored.hpp"

namespace hyteg {
namespace vertexdof {
namespace macrocell {
namespace generated {

static void add_3D_macrocell_vertexdof_1_rhs_function_colored_level_any(double * RESTRICT _data_p1FaceDst_group_0, double * RESTRICT _data_p1FaceDst_group_1, double * RESTRICT _data_p1FaceDst_group_2, double * RESTRICT _data_p1FaceDst_group_3, double * RESTRICT _data_p1FaceDst_group_4, double * RESTRICT _data_p1FaceDst_group_5, double * RESTRICT _data_p1FaceDst_group_6, double * RESTRICT _data_p1FaceDst_group_7, double const * RESTRICT const _data_p1FaceSrc_group_0_const, double const * RESTRICT const _data_p1FaceSrc_group_1_const, double const * RESTRICT const _data_p1FaceSrc_group_2_const, double const * RESTRICT const _data_p1FaceSrc_group_3_const, double const * RESTRICT const _data_p1FaceSrc_group_4_const, double const * RESTRICT const _data_p1FaceSrc_group_5_const, double const * RESTRICT const _data_p1FaceSrc_group_6_const, double const * RESTRICT const _data_p1FaceSrc_group_7_const, double c, int32_t level)
{
   {
      for (int ctr_3 = 0; ctr_3 < 1; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // vertex 0
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               const double xi_23 = c*_data_p1FaceSrc_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
               const double xi_24 = _data_p1FaceDst_group_7[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
               _data_p1FaceDst_group_7[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))] = xi_23 + xi_24;
            }
            // edge 0
            for (int ctr_1 = 1; ctr_1 < (1 << (level - 1)) - 1; ctr_1 += 1)
            {
               const double xi_41 = c*_data_p1FaceSrc_group_6_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_42 = _data_p1FaceDst_group_6[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_43 = c*_data_p1FaceSrc_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
               const double xi_44 = _data_p1FaceDst_group_7[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
               _data_p1FaceDst_group_6[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_41 + xi_42;
               _data_p1FaceDst_group_7[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))] = xi_43 + xi_44;
            }
         }
         for (int ctr_2 = 1; ctr_2 < -ctr_3 + (1 << (level - 1)) - 1; ctr_2 += 1)
         {
            // edge 1
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               const double xi_49 = c*_data_p1FaceSrc_group_5_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_50 = _data_p1FaceDst_group_5[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_51 = c*_data_p1FaceSrc_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
               const double xi_52 = _data_p1FaceDst_group_7[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
               _data_p1FaceDst_group_5[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_49 + xi_50;
               _data_p1FaceDst_group_7[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))] = xi_51 + xi_52;
            }
            // face 0
            for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level - 1)) - 1; ctr_1 += 1)
            {
               const double xi_83 = c*_data_p1FaceSrc_group_4_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_84 = _data_p1FaceDst_group_4[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_85 = c*_data_p1FaceSrc_group_5_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_86 = _data_p1FaceDst_group_5[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_87 = c*_data_p1FaceSrc_group_6_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_88 = _data_p1FaceDst_group_6[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_89 = c*_data_p1FaceSrc_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
               const double xi_90 = _data_p1FaceDst_group_7[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
               _data_p1FaceDst_group_4[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_83 + xi_84;
               _data_p1FaceDst_group_5[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_85 + xi_86;
               _data_p1FaceDst_group_6[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_87 + xi_88;
               _data_p1FaceDst_group_7[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))] = xi_89 + xi_90;
            }
            // edge 2
            for (int ctr_1 = -ctr_2 + (1 << (level - 1)) - 1; ctr_1 < -ctr_2 + (1 << (level - 1)); ctr_1 += 1)
            {
               const double xi_57 = c*_data_p1FaceSrc_group_4_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_58 = _data_p1FaceDst_group_4[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               _data_p1FaceDst_group_4[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_57 + xi_58;
            }
         }
      }
      for (int ctr_3 = 1; ctr_3 < (1 << (level - 1)) - 1; ctr_3 += 1)
      {
         for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
         {
            // edge 3
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               const double xi_63 = c*_data_p1FaceSrc_group_3_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_64 = _data_p1FaceDst_group_3[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_65 = c*_data_p1FaceSrc_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
               const double xi_66 = _data_p1FaceDst_group_7[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
               _data_p1FaceDst_group_3[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_63 + xi_64;
               _data_p1FaceDst_group_7[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))] = xi_65 + xi_66;
            }
            // face 1
            for (int ctr_1 = 1; ctr_1 < -ctr_3 + (1 << (level - 1)) - 1; ctr_1 += 1)
            {
               const double xi_95 = c*_data_p1FaceSrc_group_2_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_96 = _data_p1FaceDst_group_2[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_97 = c*_data_p1FaceSrc_group_3_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_98 = _data_p1FaceDst_group_3[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_99 = c*_data_p1FaceSrc_group_6_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_100 = _data_p1FaceDst_group_6[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_101 = c*_data_p1FaceSrc_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
               const double xi_102 = _data_p1FaceDst_group_7[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
               _data_p1FaceDst_group_2[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_95 + xi_96;
               _data_p1FaceDst_group_3[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_97 + xi_98;
               _data_p1FaceDst_group_6[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_100 + xi_99;
               _data_p1FaceDst_group_7[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))] = xi_101 + xi_102;
            }
            // edge 4
            for (int ctr_1 = -ctr_3 + (1 << (level - 1)) - 1; ctr_1 < -ctr_3 + (1 << (level - 1)); ctr_1 += 1)
            {
               const double xi_71 = c*_data_p1FaceSrc_group_2_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_72 = _data_p1FaceDst_group_2[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               _data_p1FaceDst_group_2[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_71 + xi_72;
            }
         }
         for (int ctr_2 = 1; ctr_2 < -ctr_3 + (1 << (level - 1)) - 1; ctr_2 += 1)
         {
            // face 2
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               const double xi_107 = c*_data_p1FaceSrc_group_1_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_108 = _data_p1FaceDst_group_1[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_109 = c*_data_p1FaceSrc_group_3_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_110 = _data_p1FaceDst_group_3[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_111 = c*_data_p1FaceSrc_group_5_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_112 = _data_p1FaceDst_group_5[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_113 = c*_data_p1FaceSrc_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
               const double xi_114 = _data_p1FaceDst_group_7[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
               _data_p1FaceDst_group_1[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_107 + xi_108;
               _data_p1FaceDst_group_3[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_109 + xi_110;
               _data_p1FaceDst_group_5[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_111 + xi_112;
               _data_p1FaceDst_group_7[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))] = xi_113 + xi_114;
            }
            // cell (inner)
            for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + (1 << (level - 1)) - 1; ctr_1 += 1)
            {
               const double xi_3 = c*_data_p1FaceSrc_group_0_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
               const double xi_4 = _data_p1FaceDst_group_0[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
               const double xi_5 = c*_data_p1FaceSrc_group_1_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_6 = _data_p1FaceDst_group_1[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_7 = c*_data_p1FaceSrc_group_2_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_8 = _data_p1FaceDst_group_2[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_9 = c*_data_p1FaceSrc_group_3_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_10 = _data_p1FaceDst_group_3[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_11 = c*_data_p1FaceSrc_group_4_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_12 = _data_p1FaceDst_group_4[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_13 = c*_data_p1FaceSrc_group_5_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_14 = _data_p1FaceDst_group_5[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_15 = c*_data_p1FaceSrc_group_6_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_16 = _data_p1FaceDst_group_6[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_17 = c*_data_p1FaceSrc_group_7_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
               const double xi_18 = _data_p1FaceDst_group_7[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))];
               _data_p1FaceDst_group_0[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))] = xi_3 + xi_4;
               _data_p1FaceDst_group_1[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_5 + xi_6;
               _data_p1FaceDst_group_2[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_7 + xi_8;
               _data_p1FaceDst_group_3[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_10 + xi_9;
               _data_p1FaceDst_group_4[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_11 + xi_12;
               _data_p1FaceDst_group_5[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_13 + xi_14;
               _data_p1FaceDst_group_6[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_15 + xi_16;
               _data_p1FaceDst_group_7[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) - 1)*(-ctr_3 + (1 << (level - 1)) + 1)) / (6)) + ((((1 << (level - 1)) - 1)*((1 << (level - 1)) + 1)*(1 << (level - 1))) / (6))] = xi_17 + xi_18;
            }
            // face 3
            for (int ctr_1 = -ctr_2 - ctr_3 + (1 << (level - 1)) - 1; ctr_1 < -ctr_2 - ctr_3 + (1 << (level - 1)); ctr_1 += 1)
            {
               const double xi_119 = c*_data_p1FaceSrc_group_0_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
               const double xi_120 = _data_p1FaceDst_group_0[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))];
               const double xi_121 = c*_data_p1FaceSrc_group_1_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_122 = _data_p1FaceDst_group_1[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_123 = c*_data_p1FaceSrc_group_2_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_124 = _data_p1FaceDst_group_2[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_125 = c*_data_p1FaceSrc_group_4_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_126 = _data_p1FaceDst_group_4[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               _data_p1FaceDst_group_0[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*((1 << (level - 1)) + 3)) / (6)) - (((-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)*(-ctr_3 + (1 << (level - 1)) + 3)) / (6))] = xi_119 + xi_120;
               _data_p1FaceDst_group_1[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_121 + xi_122;
               _data_p1FaceDst_group_2[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_123 + xi_124;
               _data_p1FaceDst_group_4[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_125 + xi_126;
            }
         }
         for (int ctr_2 = -ctr_3 + (1 << (level - 1)) - 1; ctr_2 < -ctr_3 + (1 << (level - 1)); ctr_2 += 1)
         {
            // edge 5
            for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
            {
               const double xi_77 = c*_data_p1FaceSrc_group_1_const[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               const double xi_78 = _data_p1FaceDst_group_1[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))];
               _data_p1FaceDst_group_1[ctr_1 + ctr_2*(-ctr_3 + (1 << (level - 1)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level - 1)))*(-ctr_3 + (1 << (level - 1)) + 1)*(-ctr_3 + (1 << (level - 1)) + 2)) / (6)) + ((((1 << (level - 1)) + 1)*((1 << (level - 1)) + 2)*(1 << (level - 1))) / (6))] = xi_77 + xi_78;
            }
         }
      }
   }
}


void add_3D_macrocell_vertexdof_1_rhs_function_colored(double * RESTRICT _data_p1FaceDst_group_0, double * RESTRICT _data_p1FaceDst_group_1, double * RESTRICT _data_p1FaceDst_group_2, double * RESTRICT _data_p1FaceDst_group_3, double * RESTRICT _data_p1FaceDst_group_4, double * RESTRICT _data_p1FaceDst_group_5, double * RESTRICT _data_p1FaceDst_group_6, double * RESTRICT _data_p1FaceDst_group_7, double const * RESTRICT const _data_p1FaceSrc_group_0_const, double const * RESTRICT const _data_p1FaceSrc_group_1_const, double const * RESTRICT const _data_p1FaceSrc_group_2_const, double const * RESTRICT const _data_p1FaceSrc_group_3_const, double const * RESTRICT const _data_p1FaceSrc_group_4_const, double const * RESTRICT const _data_p1FaceSrc_group_5_const, double const * RESTRICT const _data_p1FaceSrc_group_6_const, double const * RESTRICT const _data_p1FaceSrc_group_7_const, double c, int32_t level)
{
    switch( level )
    {

    default:
        add_3D_macrocell_vertexdof_1_rhs_function_colored_level_any(_data_p1FaceDst_group_0, _data_p1FaceDst_group_1, _data_p1FaceDst_group_2, _data_p1FaceDst_group_3, _data_p1FaceDst_group_4, _data_p1FaceDst_group_5, _data_p1FaceDst_group_6, _data_p1FaceDst_group_7, _data_p1FaceSrc_group_0_const, _data_p1FaceSrc_group_1_const, _data_p1FaceSrc_group_2_const, _data_p1FaceSrc_group_3_const, _data_p1FaceSrc_group_4_const, _data_p1FaceSrc_group_5_const, _data_p1FaceSrc_group_6_const, _data_p1FaceSrc_group_7_const, c, level);
        break;
    }
}
    

} // namespace generated
} // namespace macrocell
} // namespace vertexdof
} // namespace hyteg
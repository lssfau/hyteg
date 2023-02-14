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

#include "add_2D_macroface_vertexdof_1_rhsfunction.hpp"

namespace hyteg {
namespace vertexdof {
namespace macroface {
namespace generated {

static void add_2D_macroface_vertexdof_1_rhs_function_level_any(double * RESTRICT _data_p1FaceDst, double * RESTRICT _data_p1FaceSrc, double c, int level)
{
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_3 = c*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4;
      }
   }
}


void add_2D_macroface_vertexdof_1_rhs_function(double * RESTRICT _data_p1FaceDst, double * RESTRICT _data_p1FaceSrc, double c, int level)
{
    switch( level )
    {

    default:
        add_2D_macroface_vertexdof_1_rhs_function_level_any(_data_p1FaceDst, _data_p1FaceSrc, c, level);
        break;
    }
}

static void add_2D_macroface_vertexdof_1_rhs_function_level_any(float * RESTRICT _data_p1FaceDst, float * RESTRICT _data_p1FaceSrc, float c, int level)
{
    for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
    {
        // inner triangle
        for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
        {
         const float xi_3 = c*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const float xi_4 = _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4;
        }
    }
}


void add_2D_macroface_vertexdof_1_rhs_function(float * RESTRICT _data_p1FaceDst, float * RESTRICT _data_p1FaceSrc, float c, int level)
{
    switch( level )
    {

    default:
        add_2D_macroface_vertexdof_1_rhs_function_level_any(_data_p1FaceDst, _data_p1FaceSrc, c, level);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg
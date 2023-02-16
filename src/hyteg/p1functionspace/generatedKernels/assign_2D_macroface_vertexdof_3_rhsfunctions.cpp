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

#include "assign_2D_macroface_vertexdof_3_rhsfunctions.hpp"

namespace hyteg {
namespace vertexdof {
namespace macroface {
namespace generated {

static void assign_2D_macroface_vertexdof_3_rhs_functions_level_any(double * RESTRICT _data_p1FaceDst, double * RESTRICT _data_p1FaceSrc0, double * RESTRICT _data_p1FaceSrc1, double * RESTRICT _data_p1FaceSrc3, double c0, double c1, double c2, int level)
{
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_3 = c0*_data_p1FaceSrc0[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = c1*_data_p1FaceSrc1[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_5 = c2*_data_p1FaceSrc3[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4 + xi_5;
      }
   }
}


void assign_2D_macroface_vertexdof_3_rhs_functions(double * RESTRICT _data_p1FaceDst, double * RESTRICT _data_p1FaceSrc0, double * RESTRICT _data_p1FaceSrc1, double * RESTRICT _data_p1FaceSrc3, double c0, double c1, double c2, int level)
{
    switch( level )
    {

    default:
        assign_2D_macroface_vertexdof_3_rhs_functions_level_any(_data_p1FaceDst, _data_p1FaceSrc0, _data_p1FaceSrc1, _data_p1FaceSrc3, c0, c1, c2, level);
        break;
    }
}

static void assign_2D_macroface_vertexdof_3_rhs_functions_level_any(float * RESTRICT _data_p1FaceDst, float * RESTRICT _data_p1FaceSrc0, float * RESTRICT _data_p1FaceSrc1, float * RESTRICT _data_p1FaceSrc3, float c0, float c1, float c2, int level)
{
    for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
    {
        // inner triangle
        for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
        {
         const float xi_3 = c0*_data_p1FaceSrc0[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const float xi_4 = c1*_data_p1FaceSrc1[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const float xi_5 = c2*_data_p1FaceSrc3[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_3 + xi_4 + xi_5;
        }
    }
}


void assign_2D_macroface_vertexdof_3_rhs_functions(float * RESTRICT _data_p1FaceDst, float * RESTRICT _data_p1FaceSrc0, float * RESTRICT _data_p1FaceSrc1, float * RESTRICT _data_p1FaceSrc3, float c0, float c1, float c2, int level)
{
    switch( level )
    {

    default:
        assign_2D_macroface_vertexdof_3_rhs_functions_level_any(_data_p1FaceDst, _data_p1FaceSrc0, _data_p1FaceSrc1, _data_p1FaceSrc3, c0, c1, c2, level);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg
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

#include "sor_2D_macroface_vertexdof_to_vertexdof_backwards.hpp"

namespace hyteg {
namespace vertexdof {
namespace macroface {
namespace generated {

static void sor_2D_macroface_vertexdof_to_vertexdof_backwards_level_any(double * RESTRICT _data_p1FaceDst, double * RESTRICT _data_p1FaceRhs, double const * RESTRICT const _data_p1FaceStencil, int level, double relax)
{
   const double xi_11 = 1.0;
   const double xi_12 = -relax;
   const double xi_0 = _data_p1FaceStencil[3];
   const double xi_9 = 1 / (xi_0);
   const double xi_1 = _data_p1FaceStencil[2];
   const double xi_2 = _data_p1FaceStencil[5];
   const double xi_3 = _data_p1FaceStencil[0];
   const double xi_4 = _data_p1FaceStencil[6];
   const double xi_5 = _data_p1FaceStencil[1];
   const double xi_6 = _data_p1FaceStencil[4];
   for (int ctr_2 = (1 << (level)) - 1; ctr_2 >= 1; ctr_2 += -1)
   {
      // inner triangle
      for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 >= 1; ctr_1 += -1)
      {
         const double xi_19 = _data_p1FaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_13 = -xi_1*_data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_14 = -xi_2*_data_p1FaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_15 = -xi_3*_data_p1FaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_16 = -xi_4*_data_p1FaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_17 = -xi_5*_data_p1FaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_18 = -xi_6*_data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_9*(xi_13 + xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19) + (xi_11 + xi_12)*_data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_2D_macroface_vertexdof_to_vertexdof_backwards(double * RESTRICT _data_p1FaceDst, double * RESTRICT _data_p1FaceRhs, double const * RESTRICT const _data_p1FaceStencil, int level, double relax)
{
    switch( level )
    {

    default:
        sor_2D_macroface_vertexdof_to_vertexdof_backwards_level_any(_data_p1FaceDst, _data_p1FaceRhs, _data_p1FaceStencil, level, relax);
        break;
    }
}

static void sor_2D_macroface_vertexdof_to_vertexdof_backwards_level_any(float * RESTRICT _data_p1FaceDst, float * RESTRICT _data_p1FaceRhs, float const * RESTRICT const _data_p1FaceStencil, int level, float relax)
{
    const float xi_11 = 1.0;
    const float xi_12 = -relax;
    const float xi_0 = _data_p1FaceStencil[3];
    const float xi_9 = 1 / (xi_0);
    const float xi_1 = _data_p1FaceStencil[2];
    const float xi_2 = _data_p1FaceStencil[5];
    const float xi_3 = _data_p1FaceStencil[0];
    const float xi_4 = _data_p1FaceStencil[6];
    const float xi_5 = _data_p1FaceStencil[1];
    const float xi_6 = _data_p1FaceStencil[4];
    for (int ctr_2 = (1 << (level)) - 1; ctr_2 >= 1; ctr_2 += -1)
    {
        // inner triangle
        for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 >= 1; ctr_1 += -1)
        {
         const float xi_19 = _data_p1FaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const float xi_13 = -xi_1*_data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const float xi_14 = -xi_2*_data_p1FaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const float xi_15 = -xi_3*_data_p1FaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const float xi_16 = -xi_4*_data_p1FaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const float xi_17 = -xi_5*_data_p1FaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const float xi_18 = -xi_6*_data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_9*(xi_13 + xi_14 + xi_15 + xi_16 + xi_17 + xi_18 + xi_19) + (xi_11 + xi_12)*_data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
        }
    }
}


void sor_2D_macroface_vertexdof_to_vertexdof_backwards(float * RESTRICT _data_p1FaceDst, float * RESTRICT _data_p1FaceRhs, float const * RESTRICT const _data_p1FaceStencil, int level, float relax)
{
    switch( level )
    {

    default:
        sor_2D_macroface_vertexdof_to_vertexdof_backwards_level_any(_data_p1FaceDst, _data_p1FaceRhs, _data_p1FaceStencil, level, relax);
        break;
    }
}

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg
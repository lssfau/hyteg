/*
 * Copyright (c) 2019-2023 Nils Kohl, Dominik Thoennes, Michael Zikeli.
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

#include "apply_2D_macroface_vertexdof_to_vertexdof_replace.hpp"
#include "core/DataTypes.h"

namespace hyteg {
namespace vertexdof {
namespace macroface {
namespace generated {

template < typename ValueType >
static void apply_2D_macroface_vertexdof_to_vertexdof_replace_level_any(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceStencil, int level)
{
   const auto xi_0 = _data_p1FaceStencil[2];
   const auto xi_1 = _data_p1FaceStencil[5];
   const auto xi_2 = _data_p1FaceStencil[0];
   const auto xi_3 = _data_p1FaceStencil[3];
   const auto xi_4 = _data_p1FaceStencil[6];
   const auto xi_5 = _data_p1FaceStencil[1];
   const auto xi_6 = _data_p1FaceStencil[4];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const auto xi_10 = xi_0*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const auto xi_11 = xi_1*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const auto xi_12 = xi_2*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const auto xi_13 = xi_3*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const auto xi_14 = xi_4*_data_p1FaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const auto xi_15 = xi_5*_data_p1FaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const auto xi_16 = xi_6*_data_p1FaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_p1FaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_10 + xi_11 + xi_12 + xi_13 + xi_14 + xi_15 + xi_16;
      }
   }
}


template < typename ValueType >
void apply_2D_macroface_vertexdof_to_vertexdof_replace(ValueType * RESTRICT _data_p1FaceDst, ValueType const * RESTRICT const _data_p1FaceSrc, ValueType const * RESTRICT const _data_p1FaceStencil, int level)
{
    switch( level )
    {

    default:
        apply_2D_macroface_vertexdof_to_vertexdof_replace_level_any(_data_p1FaceDst, _data_p1FaceSrc, _data_p1FaceStencil, level);
        break;
    }
}

// ========================
//  explicit instantiation
// ========================
template void apply_2D_macroface_vertexdof_to_vertexdof_replace<walberla::float64>(walberla::float64 * RESTRICT _data_p1FaceDst, walberla::float64 const * RESTRICT const _data_p1FaceSrc, walberla::float64 const * RESTRICT const _data_p1FaceStencil, int level);
template void apply_2D_macroface_vertexdof_to_vertexdof_replace<walberla::float32>(walberla::float32 * RESTRICT _data_p1FaceDst, walberla::float32 const * RESTRICT const _data_p1FaceSrc, walberla::float32 const * RESTRICT const _data_p1FaceStencil, int level);
#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
template void apply_2D_macroface_vertexdof_to_vertexdof_replace<walberla::float16>(walberla::float16 * RESTRICT _data_p1FaceDst, walberla::float16 const * RESTRICT const _data_p1FaceSrc, walberla::float16 const * RESTRICT const _data_p1FaceStencil, int level);
#endif

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg
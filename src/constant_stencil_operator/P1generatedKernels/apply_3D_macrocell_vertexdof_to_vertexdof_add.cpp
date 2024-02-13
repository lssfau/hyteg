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

#include "apply_3D_macrocell_vertexdof_to_vertexdof_add.hpp"
#include "core/DataTypes.h"

namespace hyteg {
namespace vertexdof {
namespace macrocell {
namespace generated {

template < typename ValueType >
static void apply_3D_macrocell_vertexdof_to_vertexdof_add_level_any(ValueType * RESTRICT _data_p1CellDstAdd, ValueType const * RESTRICT const _data_p1CellSrcAdd, int level, std::map< hyteg::indexing::Index, ValueType >& p1CellStencil)
{
   const auto xi_1 = p1CellStencil[{ -1, 0, 0 }];
   const auto xi_2 = p1CellStencil[{ -1, 0, 1 }];
   const auto xi_3 = p1CellStencil[{ -1, 1, -1 }];
   const auto xi_4 = p1CellStencil[{ -1, 1, 0 }];
   const auto xi_5 = p1CellStencil[{ 0, -1, 0 }];
   const auto xi_6 = p1CellStencil[{ 0, -1, 1 }];
   const auto xi_7 = p1CellStencil[{ 0, 0, -1 }];
   const auto xi_8 = p1CellStencil[{ 0, 0, 0 }];
   const auto xi_9 = p1CellStencil[{ 0, 0, 1 }];
   const auto xi_10 = p1CellStencil[{ 0, 1, -1 }];
   const auto xi_11 = p1CellStencil[{ 0, 1, 0 }];
   const auto xi_12 = p1CellStencil[{ 1, -1, 0 }];
   const auto xi_13 = p1CellStencil[{ 1, -1, 1 }];
   const auto xi_14 = p1CellStencil[{ 1, 0, -1 }];
   const auto xi_15 = p1CellStencil[{ 1, 0, 0 }];
   for (int ctr_3 = 1; ctr_3 < (1 << (level)); ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + (1 << (level)); ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + (1 << (level)); ctr_1 += 1)
         {
            const auto xi_33 = _data_p1CellDstAdd[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const auto xi_18 = xi_1*_data_p1CellSrcAdd[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const auto xi_25 = xi_2*_data_p1CellSrcAdd[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - 1];
            const auto xi_26 = xi_3*_data_p1CellSrcAdd[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 3) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)*(-ctr_3 + (1 << (level)) + 4)) / (6)) - 1];
            const auto xi_27 = xi_4*_data_p1CellSrcAdd[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const auto xi_28 = xi_5*_data_p1CellSrcAdd[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const auto xi_29 = xi_6*_data_p1CellSrcAdd[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6))];
            const auto xi_30 = xi_7*_data_p1CellSrcAdd[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)*(-ctr_3 + (1 << (level)) + 4)) / (6))];
            const auto xi_31 = xi_8*_data_p1CellSrcAdd[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const auto xi_32 = xi_9*_data_p1CellSrcAdd[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6))];
            const auto xi_19 = xi_10*_data_p1CellSrcAdd[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 3) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)*(-ctr_3 + (1 << (level)) + 4)) / (6))];
            const auto xi_20 = xi_11*_data_p1CellSrcAdd[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const auto xi_21 = xi_12*_data_p1CellSrcAdd[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) + 1];
            const auto xi_22 = xi_13*_data_p1CellSrcAdd[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) + 1];
            const auto xi_23 = xi_14*_data_p1CellSrcAdd[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)*(-ctr_3 + (1 << (level)) + 4)) / (6)) + 1];
            const auto xi_24 = xi_15*_data_p1CellSrcAdd[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) + 1];
            _data_p1CellDstAdd[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))] = xi_18 + xi_19 + xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33;
         }
      }
   }
}

template < typename ValueType >
void apply_3D_macrocell_vertexdof_to_vertexdof_add(ValueType * RESTRICT _data_p1CellDstAdd, ValueType const * RESTRICT const _data_p1CellSrcAdd, int level, std::map< hyteg::indexing::Index, ValueType >& p1CellStencil)
{
    switch( level )
    {

    default:
        apply_3D_macrocell_vertexdof_to_vertexdof_add_level_any(_data_p1CellDstAdd, _data_p1CellSrcAdd, level, p1CellStencil);
        break;
    }
}

// ========================
//  explicit instantiation
// ========================
template void apply_3D_macrocell_vertexdof_to_vertexdof_add<walberla::float64>(walberla::float64 * RESTRICT _data_p1CellDstAdd, walberla::float64 const * RESTRICT const _data_p1CellSrcAdd, int level, std::map< hyteg::indexing::Index, walberla::float64 >& p1CellStencil);
template void apply_3D_macrocell_vertexdof_to_vertexdof_add<walberla::float32>(walberla::float32 * RESTRICT _data_p1CellDstAdd, walberla::float32 const * RESTRICT const _data_p1CellSrcAdd, int level, std::map< hyteg::indexing::Index, walberla::float32 >& p1CellStencil);
#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
template void apply_3D_macrocell_vertexdof_to_vertexdof_add<walberla::float16>(walberla::float16 * RESTRICT _data_p1CellDstAdd, walberla::float16 const * RESTRICT const _data_p1CellSrcAdd, int level, std::map< hyteg::indexing::Index, walberla::float16 >& p1CellStencil);
#endif

} // namespace generated
} // namespace macrocell
} // namespace vertexdof
} // namespace hyteg
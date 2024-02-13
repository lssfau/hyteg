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

#include "communicate_directly_vertexdof_cell_to_face_impl.hpp"
#include "core/DataTypes.h"

namespace hyteg {
namespace vertexdof {
namespace comm {
namespace generated {

template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_012_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_012(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_012_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_013_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 - ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level)) + 2];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 - ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level)) + 2];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 - ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level)) + 2];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 - ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level)) + 2];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 - ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level)) + 2];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 - ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level)) + 2];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 - ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level)) + 2];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_013(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_013_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_021_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*((1 << (level)) + 1) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_021(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_021_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_023_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_023(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_023_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_031_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + (1 << (level)) + 2];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + (1 << (level)) + 2];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + (1 << (level)) + 2];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + (1 << (level)) + 2];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + (1 << (level)) + 2];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + (1 << (level)) + 2];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + (1 << (level)) + 2];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_031(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_031_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_032_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2*(-ctr_1 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2*(-ctr_1 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2*(-ctr_1 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2*(-ctr_1 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2*(-ctr_1 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2*(-ctr_1 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2*(-ctr_1 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 1];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_032(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_032_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_102_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2*((1 << (level)) + 1) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_102(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_102_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_103_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 - 2*ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 2*(1 << (level)) + 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 - 2*ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 2*(1 << (level)) + 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 - 2*ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 2*(1 << (level)) + 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 - 2*ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 2*(1 << (level)) + 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 - 2*ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 2*(1 << (level)) + 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 - 2*ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 2*(1 << (level)) + 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 - 2*ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 2*(1 << (level)) + 1];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_103(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_103_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_120_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*((1 << (level)) + 1) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + (1 << (level)) - 1];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_120(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_120_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_123_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level)) - 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level)) - 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level)) - 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level)) - 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(-ctr_2 + (1 << (level)) + 2) - ctr_1 - ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + (1 << (level)) - 1];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_123(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_123_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_130_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-2*ctr_1 - ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 2*(1 << (level)) + 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-2*ctr_1 - ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 2*(1 << (level)) + 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-2*ctr_1 - ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 2*(1 << (level)) + 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-2*ctr_1 - ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 2*(1 << (level)) + 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-2*ctr_1 - ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 2*(1 << (level)) + 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-2*ctr_1 - ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 2*(1 << (level)) + 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-2*ctr_1 - ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 2*(1 << (level)) + 1];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_130(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_130_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_132_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 2) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + (1 << (level)) - 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 2) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + (1 << (level)) - 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 2) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 2) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + (1 << (level)) - 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 2) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + (1 << (level)) - 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 2) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + (1 << (level)) - 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[-ctr_1 + ctr_2*(-ctr_1 + (1 << (level)) + 2) - ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + (1 << (level)) - 1];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_132(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_132_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_201_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_201(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_201_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_203_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[(-ctr_2 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[(-ctr_2 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[(-ctr_2 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[(-ctr_2 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[(-ctr_2 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[(-ctr_2 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[(-ctr_2 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6)) + 1];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_203(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_203_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_210_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ((1 << (level)) + 1)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_210(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_210_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_213_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + (-ctr_2 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + (-ctr_2 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + (-ctr_2 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + (-ctr_2 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + (-ctr_2 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + (-ctr_2 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + (-ctr_2 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_2 + (1 << (level)) + 1)*(-ctr_2 + (1 << (level)) + 2)*(-ctr_2 + (1 << (level)) + 3)) / (6))];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_213(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_213_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_230_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6)) + 1];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_230(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_230_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_231_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2 + (-ctr_1 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2 + (-ctr_1 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2 + (-ctr_1 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2 + (-ctr_1 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2 + (-ctr_1 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2 + (-ctr_1 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2 + (-ctr_1 + (1 << (level)) + 2)*(-ctr_1 - ctr_2 + (1 << (level)) - 1) - (((-ctr_1 - ctr_2 + (1 << (level)))*(-ctr_1 - ctr_2 + (1 << (level)) - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_1 + (1 << (level)) + 1)*(-ctr_1 + (1 << (level)) + 2)*(-ctr_1 + (1 << (level)) + 3)) / (6))];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_231(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_231_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_301_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + 2*ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 3];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + 2*ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 3];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + 2*ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 3];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + 2*ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 3];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + 2*ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 3];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + 2*ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 3];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + 2*ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 3];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_301(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_301_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_302_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2*(ctr_1 + ctr_2 + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2*(ctr_1 + ctr_2 + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2*(ctr_1 + ctr_2 + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2*(ctr_1 + ctr_2 + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2*(ctr_1 + ctr_2 + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2*(ctr_1 + ctr_2 + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_2*(ctr_1 + ctr_2 + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 1];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_302(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_302_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_310_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[2*ctr_1 + ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 3];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[2*ctr_1 + ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 3];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[2*ctr_1 + ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 3];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[2*ctr_1 + ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 3];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[2*ctr_1 + ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 3];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[2*ctr_1 + ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 3];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[2*ctr_1 + ctr_2 - ((2) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 3];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_310(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_310_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_312_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1 + ctr_2*(ctr_1 + ctr_2 + 3) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6))];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_312(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_312_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_320_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(ctr_1 + ctr_2 + 3) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 1];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(ctr_1 + ctr_2 + 3) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 1];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(ctr_1 + ctr_2 + 3) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 1];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(ctr_1 + ctr_2 + 3) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 1];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(ctr_1 + ctr_2 + 3) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 1];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(ctr_1 + ctr_2 + 3) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 1];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(ctr_1 + ctr_2 + 3) - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6)) + 1];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_320(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_320_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}


template < typename ValueType >
static void communicate_directly_vertexdof_cell_to_face_impl_321_level_any(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(ctr_1 + ctr_2 + 3) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(ctr_1 + ctr_2 + 3) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(ctr_1 + ctr_2 + 3) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(ctr_1 + ctr_2 + 3) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(ctr_1 + ctr_2 + 3) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(ctr_1 + ctr_2 + 3) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6))];
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            _data_p1_face_dst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = _data_p1_cell_src[ctr_1*(ctr_1 + ctr_2 + 3) + ctr_2 - ((ctr_1*(ctr_1 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((ctr_1 + ctr_2 + 2)*(ctr_1 + ctr_2 + 3)*(ctr_1 + ctr_2 + 4)) / (6))];
         }
      }
   }
}

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face_impl_321(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level)
{
    switch( level )
    {

    default:
        communicate_directly_vertexdof_cell_to_face_impl_321_level_any(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
        break;
    }
}

// ========================
//  explicit instantiation
// ========================
// double
template void communicate_directly_vertexdof_cell_to_face_impl_012<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_013<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_021<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_023<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_031<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_032<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_102<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_103<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_120<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_123<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_130<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_132<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_201<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_203<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_210<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_213<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_230<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_231<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_301<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_302<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_310<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_312<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_320<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_321<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level);

// single
template void communicate_directly_vertexdof_cell_to_face_impl_012<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_013<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_021<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_023<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_031<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_032<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_102<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_103<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_120<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_123<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_130<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_132<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_201<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_203<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_210<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_213<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_230<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_231<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_301<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_302<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_310<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_312<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_320<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_321<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level);

// half
#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
template void communicate_directly_vertexdof_cell_to_face_impl_012<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_013<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_021<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_023<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_031<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_032<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_102<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_103<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_120<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_123<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_130<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_132<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_201<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_203<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_210<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_213<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_230<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_231<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_301<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_302<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_310<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_312<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_320<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
template void communicate_directly_vertexdof_cell_to_face_impl_321<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level);
#endif

} // namespace generated
} // namespace comm
} // namespace vertexdof
} // namespace hyteg
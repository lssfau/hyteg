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

#include "communicate_directly_vertexdof_cell_to_face.hpp"
#include "core/DataTypes.h"

namespace hyteg {
namespace vertexdof {
namespace comm {
namespace generated {

template < typename ValueType >
void communicate_directly_vertexdof_cell_to_face(ValueType const * RESTRICT const _data_p1_cell_src, ValueType * RESTRICT _data_p1_face_dst_gl0, int level, int neighbor_cell_local_vertex_id_0, int neighbor_cell_local_vertex_id_1, int neighbor_cell_local_vertex_id_2)
{
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_2)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_012(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_013(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_1)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_021(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_023(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_031(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_032(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_2)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_102(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_103(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_1)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_120(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_123(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_130(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_132(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_0)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_201(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_203(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_0)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_210(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_213(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_230(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_231(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_301(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_302(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_310(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_312(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_320(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      communicate_directly_vertexdof_cell_to_face_impl_321(_data_p1_cell_src, _data_p1_face_dst_gl0, level);
      
      return;
   } 
}

// ========================
//  explicit instantiation
// ========================
template void communicate_directly_vertexdof_cell_to_face<walberla::float64>(walberla::float64 const * RESTRICT const _data_p1_cell_src, walberla::float64 * RESTRICT _data_p1_face_dst_gl0, int level, int neighbor_cell_local_vertex_id_0, int neighbor_cell_local_vertex_id_1, int neighbor_cell_local_vertex_id_2);
template void communicate_directly_vertexdof_cell_to_face<walberla::float32>(walberla::float32 const * RESTRICT const _data_p1_cell_src, walberla::float32 * RESTRICT _data_p1_face_dst_gl0, int level, int neighbor_cell_local_vertex_id_0, int neighbor_cell_local_vertex_id_1, int neighbor_cell_local_vertex_id_2);
#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
template void communicate_directly_vertexdof_cell_to_face<walberla::float16>(walberla::float16 const * RESTRICT const _data_p1_cell_src, walberla::float16 * RESTRICT _data_p1_face_dst_gl0, int level, int neighbor_cell_local_vertex_id_0, int neighbor_cell_local_vertex_id_1, int neighbor_cell_local_vertex_id_2);
#endif

} // namespace generated
} // namespace comm
} // namespace vertexdof
} // namespace hyteg
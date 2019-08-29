
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "communicate_directly_vertexdof_cell_to_face.hpp"

namespace hhg {
namespace vertexdof {
namespace comm {
namespace generated {

void communicate_directly_vertexdof_cell_to_face(double const * RESTRICT const _data_p1_cell_src, double * RESTRICT _data_p1_face_dst_gl0, int32_t level, int64_t neighbor_cell_local_vertex_id_0, int64_t neighbor_cell_local_vertex_id_1, int64_t neighbor_cell_local_vertex_id_2)
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


} // namespace generated
} // namespace comm
} // namespace vertexdof
} // namespace hhg
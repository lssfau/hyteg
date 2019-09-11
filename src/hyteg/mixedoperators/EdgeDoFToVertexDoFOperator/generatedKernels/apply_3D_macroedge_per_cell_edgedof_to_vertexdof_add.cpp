
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add.hpp"

namespace hyteg {
namespace EdgeDoFToVertexDoF {
namespace generated {

void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::IndexIncrement, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int32_t level, int64_t neighbor_cell_local_vertex_id_0, int64_t neighbor_cell_local_vertex_id_1, int64_t neighbor_cell_local_vertex_id_2, int64_t num_neighbor_faces)
{
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_012(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_013(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_021(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_023(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_031(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_032(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_102(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_103(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_120(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_123(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_130(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_132(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_201(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_203(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_210(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_213(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_230(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_231(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_301(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_302(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_310(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_312(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_320(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_321(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
      
      return;
   } 
}


} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hyteg
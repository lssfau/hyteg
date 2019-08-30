
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "sor_3D_macroface_P1_one_sided.hpp"

namespace hyteg {
namespace vertexdof {
namespace macroface {
namespace generated {

void sor_3D_macroface_P1_one_sided(double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double const * RESTRICT const _data_vertexFaceRhs, int32_t level, int64_t neighbor_cell_0_local_vertex_id_0, int64_t neighbor_cell_0_local_vertex_id_1, int64_t neighbor_cell_0_local_vertex_id_2, double relax, std::map< hyteg::indexing::IndexIncrement, double > v2v_cell_stencil_fused_face_0)
{
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_012(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_013(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_021(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_023(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_031(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_032(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_102(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_103(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_120(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_123(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_130(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_0)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_132(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_201(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_203(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_210(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_2)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_213(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_230(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_0)) && ((3) == (neighbor_cell_0_local_vertex_id_1)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_231(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((1) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_301(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_302(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((1) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_310(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_1)) && ((2) == (neighbor_cell_0_local_vertex_id_2)) && ((3) == (neighbor_cell_0_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_312(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((0) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_320(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
   if (((1) == (neighbor_cell_0_local_vertex_id_2)) && ((2) == (neighbor_cell_0_local_vertex_id_1)) && ((3) == (neighbor_cell_0_local_vertex_id_0)))
   {
      
      sor_3D_macroface_P1_one_sided_impl_321(_data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceRhs, level, relax, v2v_cell_stencil_fused_face_0);
      
      return;
   } 
}


} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg
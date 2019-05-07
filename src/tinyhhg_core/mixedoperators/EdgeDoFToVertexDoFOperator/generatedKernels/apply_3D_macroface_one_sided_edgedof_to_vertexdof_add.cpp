
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "apply_3D_macroface_one_sided_edgedof_to_vertexdof_add.hpp"

namespace hhg {
namespace EdgeDoFToVertexDoF {
namespace generated {

void apply_3D_macroface_one_sided_edgedof_to_vertexdof_add(double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * RESTRICT const _data_edgeFaceSrc_gl0_X, double const * RESTRICT const _data_edgeFaceSrc_gl0_XY, double const * RESTRICT const _data_edgeFaceSrc_gl0_XYZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_XZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_Y, double const * RESTRICT const _data_edgeFaceSrc_gl0_YZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_Z, double * RESTRICT _data_vertexFaceDst, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2v_cell_stencil, int64_t level, int64_t neighbor_cell_local_vertex_id_0, int64_t neighbor_cell_local_vertex_id_1, int64_t neighbor_cell_local_vertex_id_2)
{
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_012(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_013(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_021(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_023(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_031(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_032(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_102(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_103(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_120(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_123(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_130(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_0)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_132(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_201(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_203(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_210(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_2)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_213(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_230(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_0)) && ((3) == (neighbor_cell_local_vertex_id_1)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_231(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((1) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_301(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_302(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((1) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_310(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_1)) && ((2) == (neighbor_cell_local_vertex_id_2)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_312(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((0) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_320(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
   if (((1) == (neighbor_cell_local_vertex_id_2)) && ((2) == (neighbor_cell_local_vertex_id_1)) && ((3) == (neighbor_cell_local_vertex_id_0)))
   {
      
      apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_321(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
      
      return;
   } 
}


} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hhg
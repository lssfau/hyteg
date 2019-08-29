
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_130.hpp"

namespace hhg {
namespace EdgeDoFToVertexDoF {
namespace generated {

static void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_130_level_any(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int32_t level, int64_t num_neighbor_faces)
{
   const double xi_1 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_2 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_3 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_4 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_5 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_6 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_7 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_8 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_9 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_11 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_12 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_13 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_14 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_15 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_16 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_17 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_18 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_19 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   // inner edge
   for (int ctr_1 = 1; ctr_1 < (1 << (level)); ctr_1 += 1)
   {
      const double xi_41 = _data_vertexEdgeDst[ctr_1];
      const double xi_22 = xi_1*_data_edgeEdgeSrc[cell_id*(7*(1 << (level)) - 7) + ctr_1 + num_neighbor_faces*(3*(1 << (level)) - 1) + 6*(1 << (level)) - 7];
      const double xi_33 = xi_2*_data_edgeEdgeSrc[cell_id*(7*(1 << (level)) - 7) + ctr_1 + num_neighbor_faces*(3*(1 << (level)) - 1) + 3*(1 << (level)) - 4];
      const double xi_34 = xi_3*_data_edgeEdgeSrc[ctr_1 + face_id_1*(3*(1 << (level)) - 1) + 3*(1 << (level)) - 1];
      const double xi_35 = xi_4*_data_edgeEdgeSrc[ctr_1 + face_id_1*(3*(1 << (level)) - 1) + 3*(1 << (level)) - 2];
      const double xi_36 = xi_5*_data_edgeEdgeSrc[ctr_1 + face_id_0*(3*(1 << (level)) - 1) + (1 << (level)) - 1];
      const double xi_37 = xi_6*_data_edgeEdgeSrc[ctr_1];
      const double xi_38 = xi_7*_data_edgeEdgeSrc[ctr_1 + face_id_1*(3*(1 << (level)) - 1) + (1 << (level)) - 1];
      const double xi_39 = xi_8*_data_edgeEdgeSrc[ctr_1 - 1];
      const double xi_40 = xi_9*_data_edgeEdgeSrc[ctr_1 + face_id_0*(3*(1 << (level)) - 1) + 3*(1 << (level)) - 1];
      const double xi_23 = xi_10*_data_edgeEdgeSrc[cell_id*(7*(1 << (level)) - 7) + ctr_1 + num_neighbor_faces*(3*(1 << (level)) - 1) + 2*(1 << (level)) - 3];
      const double xi_24 = xi_11*_data_edgeEdgeSrc[ctr_1 + face_id_0*(3*(1 << (level)) - 1) + 3*(1 << (level)) - 2];
      const double xi_25 = xi_12*_data_edgeEdgeSrc[cell_id*(7*(1 << (level)) - 7) + ctr_1 + num_neighbor_faces*(3*(1 << (level)) - 1) + 5*(1 << (level)) - 6];
      const double xi_26 = xi_13*_data_edgeEdgeSrc[ctr_1 + face_id_1*(3*(1 << (level)) - 1) + 2*(1 << (level)) - 1];
      const double xi_27 = xi_14*_data_edgeEdgeSrc[ctr_1 + face_id_1*(3*(1 << (level)) - 1) + 2*(1 << (level)) - 2];
      const double xi_28 = xi_15*_data_edgeEdgeSrc[cell_id*(7*(1 << (level)) - 7) + ctr_1 + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 7];
      const double xi_29 = xi_16*_data_edgeEdgeSrc[cell_id*(7*(1 << (level)) - 7) + ctr_1 + num_neighbor_faces*(3*(1 << (level)) - 1) + 7*(1 << (level)) - 8];
      const double xi_30 = xi_17*_data_edgeEdgeSrc[ctr_1 + face_id_0*(3*(1 << (level)) - 1) + 2*(1 << (level)) - 1];
      const double xi_31 = xi_18*_data_edgeEdgeSrc[cell_id*(7*(1 << (level)) - 7) + ctr_1 + num_neighbor_faces*(3*(1 << (level)) - 1) + 4*(1 << (level)) - 5];
      const double xi_32 = xi_19*_data_edgeEdgeSrc[ctr_1 + face_id_0*(3*(1 << (level)) - 1) + 2*(1 << (level)) - 2];
      _data_vertexEdgeDst[ctr_1] = xi_22 + xi_23 + xi_24 + xi_25 + xi_26 + xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41;
   }
}


void apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_130(double const * RESTRICT const _data_edgeEdgeSrc, double * RESTRICT _data_vertexEdgeDst, int64_t cell_id, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2v_cell_stencil, int64_t face_id_0, int64_t face_id_1, int32_t level, int64_t num_neighbor_faces)
{
    switch( level )
    {

    default:
        apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add_impl_130_level_any(_data_edgeEdgeSrc, _data_vertexEdgeDst, cell_id, e2v_cell_stencil, face_id_0, face_id_1, level, num_neighbor_faces);
        break;
    }
}
    

} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hhg

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_103.hpp"

namespace hhg {
namespace EdgeDoFToVertexDoF {
namespace generated {

static void apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_103_level_any(double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * RESTRICT const _data_edgeFaceSrc_gl0_X, double const * RESTRICT const _data_edgeFaceSrc_gl0_XY, double const * RESTRICT const _data_edgeFaceSrc_gl0_XYZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_XZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_Y, double const * RESTRICT const _data_edgeFaceSrc_gl0_YZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_Z, double * RESTRICT _data_vertexFaceDst, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2v_cell_stencil, int32_t level)
{
   const double xi_1 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_2 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_3 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_4 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_5 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_6 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_7 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_8 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_9 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_10 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_11 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_12 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_13 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_14 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_15 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_16 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_17 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_18 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_19 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_20 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_21 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_22 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_23 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_24 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_25 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_26 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_27 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_28 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_29 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_30 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_31 = e2v_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_65 = _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_34 = xi_1*_data_edgeFaceSrc_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_45 = xi_2*_data_edgeFaceSrc_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_56 = xi_3*_data_edgeFaceSrc_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1];
         const double xi_59 = xi_4*_data_edgeFaceSrc_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_60 = xi_5*_data_edgeFaceSrc_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_61 = xi_6*_data_edgeFaceSrc_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_62 = xi_7*_data_edgeFaceSrc_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_63 = xi_8*_data_edgeFaceSrc_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_64 = xi_9*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = xi_10*_data_edgeFaceSrc_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_36 = xi_11*_data_edgeFaceSrc_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_37 = xi_12*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_38 = xi_13*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_39 = xi_14*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_40 = xi_15*_data_edgeFaceSrc_gl0_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_41 = xi_16*_data_edgeFaceSrc_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_42 = xi_17*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_43 = xi_18*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_44 = xi_19*_data_edgeFaceSrc_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_46 = xi_20*_data_edgeFaceSrc_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = xi_21*_data_edgeFaceSrc_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_48 = xi_22*_data_edgeFaceSrc_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_49 = xi_23*_data_edgeFaceSrc_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_50 = xi_24*_data_edgeFaceSrc_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_51 = xi_25*_data_edgeFaceSrc_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_52 = xi_26*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_53 = xi_27*_data_edgeFaceSrc_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_54 = xi_28*_data_edgeFaceSrc_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
         const double xi_55 = xi_29*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_57 = xi_30*_data_edgeFaceSrc_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1];
         const double xi_58 = xi_31*_data_edgeFaceSrc_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1];
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65;
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_103(double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * RESTRICT const _data_edgeFaceSrc_gl0_X, double const * RESTRICT const _data_edgeFaceSrc_gl0_XY, double const * RESTRICT const _data_edgeFaceSrc_gl0_XYZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_XZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_Y, double const * RESTRICT const _data_edgeFaceSrc_gl0_YZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_Z, double * RESTRICT _data_vertexFaceDst, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2v_cell_stencil, int32_t level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_vertexdof_add_103_level_any(_data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_X, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, _data_vertexFaceDst, e2v_cell_stencil, level);
        break;
    }
}
    

} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hhg
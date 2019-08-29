
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "sor_3D_macroface_P2_update_vertexdofs_impl_132_203.hpp"

namespace hhg {
namespace P2 {
namespace macroface {
namespace generated {

static void sor_3D_macroface_P2_update_vertexdofs_impl_132_203_level_any(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double const * RESTRICT const _data_edgeFaceDst_gl1_X, double const * RESTRICT const _data_edgeFaceDst_gl1_XY, double const * RESTRICT const _data_edgeFaceDst_gl1_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl1_XZ, double const * RESTRICT const _data_edgeFaceDst_gl1_Y, double const * RESTRICT const _data_edgeFaceDst_gl1_YZ, double const * RESTRICT const _data_edgeFaceDst_gl1_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double * RESTRICT _data_vertexFaceDst_gl1, double const * RESTRICT const _data_vertexFaceRhs, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2v_cell_stencil_fused_face_0, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2v_cell_stencil_fused_face_1, int32_t level, double relax, std::map< hhg::indexing::IndexIncrement, double > v2v_cell_stencil_fused_face_0, std::map< hhg::indexing::IndexIncrement, double > v2v_cell_stencil_fused_face_1)
{
   const double xi_1 = v2v_cell_stencil_fused_face_0[{ 0, 0, 0 }];
   const double xi_2 = v2v_cell_stencil_fused_face_1[{ 0, 0, 0 }];
   const double xi_86 = 1 / (xi_1 + xi_2);
   const double xi_3 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_4 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_5 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_6 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_7 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_8 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_9 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_10 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_11 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_12 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_13 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_14 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_15 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_16 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_17 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_18 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_19 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_20 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_21 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_22 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_23 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_24 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_25 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_26 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_27 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_28 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_29 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_30 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_31 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_32 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_33 = e2v_cell_stencil_fused_face_0[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_34 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_35 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_36 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_37 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_38 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_39 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_40 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_41 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_42 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_43 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_44 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_45 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_46 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_47 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_48 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_49 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_50 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_51 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_52 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_53 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_54 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_55 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_56 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_57 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_58 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_59 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_60 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_61 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_62 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_63 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_64 = e2v_cell_stencil_fused_face_1[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   const double xi_65 = v2v_cell_stencil_fused_face_0[{ -1, 0, 0 }];
   const double xi_66 = v2v_cell_stencil_fused_face_0[{ -1, 0, 1 }];
   const double xi_67 = v2v_cell_stencil_fused_face_0[{ -1, 1, -1 }];
   const double xi_68 = v2v_cell_stencil_fused_face_0[{ -1, 1, 0 }];
   const double xi_69 = v2v_cell_stencil_fused_face_0[{ 0, -1, 0 }];
   const double xi_70 = v2v_cell_stencil_fused_face_0[{ 0, -1, 1 }];
   const double xi_71 = v2v_cell_stencil_fused_face_0[{ 0, 0, -1 }];
   const double xi_72 = v2v_cell_stencil_fused_face_0[{ 0, 1, -1 }];
   const double xi_73 = v2v_cell_stencil_fused_face_0[{ 1, -1, 0 }];
   const double xi_74 = v2v_cell_stencil_fused_face_0[{ 1, 0, -1 }];
   const double xi_75 = v2v_cell_stencil_fused_face_1[{ 0, -1, 0 }];
   const double xi_76 = v2v_cell_stencil_fused_face_1[{ 0, -1, 1 }];
   const double xi_77 = v2v_cell_stencil_fused_face_1[{ 0, 0, -1 }];
   const double xi_78 = v2v_cell_stencil_fused_face_1[{ 0, 0, 1 }];
   const double xi_79 = v2v_cell_stencil_fused_face_1[{ 0, 1, -1 }];
   const double xi_80 = v2v_cell_stencil_fused_face_1[{ 0, 1, 0 }];
   const double xi_81 = v2v_cell_stencil_fused_face_1[{ 1, -1, 0 }];
   const double xi_82 = v2v_cell_stencil_fused_face_1[{ 1, -1, 1 }];
   const double xi_83 = v2v_cell_stencil_fused_face_1[{ 1, 0, -1 }];
   const double xi_84 = v2v_cell_stencil_fused_face_1[{ 1, 0, 0 }];
   for (int ctr_2 = 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_86*(-xi_10*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_11*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_12*_data_edgeFaceDst_gl0_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_13*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_14*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_15*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_16*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_17*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_18*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_19*_data_edgeFaceDst_gl0_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_20*_data_edgeFaceDst_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_21*_data_edgeFaceDst_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_22*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_edgeFaceDst_gl0_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_24*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_25*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_26*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_27*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_28*_data_edgeFaceDst_gl0_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_29*_data_edgeFaceDst_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_3*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_30*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_31*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_32*_data_edgeFaceDst_gl0_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_33*_data_edgeFaceDst_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_34*_data_edgeFaceDst_gl1_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_35*_data_edgeFaceDst_gl1_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_36*_data_edgeFaceDst_gl1_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_37*_data_edgeFaceDst_gl1_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_edgeFaceDst_gl1_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_39*_data_edgeFaceDst_gl1_Z[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_4*_data_edgeFaceDst_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_40*_data_edgeFaceDst_gl1_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_41*_data_edgeFaceDst_gl1_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_42*_data_edgeFaceDst_gl1_YZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_43*_data_edgeFaceDst_gl1_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_44*_data_edgeFaceDst_gl1_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_45*_data_edgeFaceDst_gl1_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_46*_data_edgeFaceDst_gl1_XZ[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_47*_data_edgeFaceDst_gl1_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_48*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_49*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_5*_data_edgeFaceDst_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_50*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_51*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_52*_data_edgeFaceDst_gl1_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_53*_data_edgeFaceDst_gl1_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_54*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_55*_data_edgeFaceDst_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_56*_data_edgeFaceDst_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_57*_data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_58*_data_edgeFaceDst_gl1_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_59*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_edgeFaceDst_gl0_Y[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_60*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_61*_data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_62*_data_edgeFaceDst_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_63*_data_edgeFaceDst_gl1_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_64*_data_edgeFaceDst_gl1_XY[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1] - xi_65*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_66*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_67*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_68*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_69*_data_vertexFaceDst_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_7*_data_edgeFaceDst_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_70*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_71*_data_vertexFaceDst_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_72*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_73*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_74*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_75*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_76*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_77*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_78*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_79*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_8*_data_edgeFaceDst_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_80*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_81*_data_vertexFaceDst_gl1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_82*_data_vertexFaceDst_gl1[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_83*_data_vertexFaceDst_gl1[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))] - xi_84*_data_vertexFaceDst_gl1[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_9*_data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] + _data_vertexFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
}


void sor_3D_macroface_P2_update_vertexdofs_impl_132_203(double const * RESTRICT const _data_edgeFaceDst_X, double const * RESTRICT const _data_edgeFaceDst_XY, double const * RESTRICT const _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_X, double const * RESTRICT const _data_edgeFaceDst_gl0_XY, double const * RESTRICT const _data_edgeFaceDst_gl0_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl0_XZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Y, double const * RESTRICT const _data_edgeFaceDst_gl0_YZ, double const * RESTRICT const _data_edgeFaceDst_gl0_Z, double const * RESTRICT const _data_edgeFaceDst_gl1_X, double const * RESTRICT const _data_edgeFaceDst_gl1_XY, double const * RESTRICT const _data_edgeFaceDst_gl1_XYZ, double const * RESTRICT const _data_edgeFaceDst_gl1_XZ, double const * RESTRICT const _data_edgeFaceDst_gl1_Y, double const * RESTRICT const _data_edgeFaceDst_gl1_YZ, double const * RESTRICT const _data_edgeFaceDst_gl1_Z, double * RESTRICT _data_vertexFaceDst, double * RESTRICT _data_vertexFaceDst_gl0, double * RESTRICT _data_vertexFaceDst_gl1, double const * RESTRICT const _data_vertexFaceRhs, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2v_cell_stencil_fused_face_0, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2v_cell_stencil_fused_face_1, int32_t level, double relax, std::map< hhg::indexing::IndexIncrement, double > v2v_cell_stencil_fused_face_0, std::map< hhg::indexing::IndexIncrement, double > v2v_cell_stencil_fused_face_1)
{
    switch( level )
    {

    default:
        sor_3D_macroface_P2_update_vertexdofs_impl_132_203_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceDst_gl0_X, _data_edgeFaceDst_gl0_XY, _data_edgeFaceDst_gl0_XYZ, _data_edgeFaceDst_gl0_XZ, _data_edgeFaceDst_gl0_Y, _data_edgeFaceDst_gl0_YZ, _data_edgeFaceDst_gl0_Z, _data_edgeFaceDst_gl1_X, _data_edgeFaceDst_gl1_XY, _data_edgeFaceDst_gl1_XYZ, _data_edgeFaceDst_gl1_XZ, _data_edgeFaceDst_gl1_Y, _data_edgeFaceDst_gl1_YZ, _data_edgeFaceDst_gl1_Z, _data_vertexFaceDst, _data_vertexFaceDst_gl0, _data_vertexFaceDst_gl1, _data_vertexFaceRhs, e2v_cell_stencil_fused_face_0, e2v_cell_stencil_fused_face_1, level, relax, v2v_cell_stencil_fused_face_0, v2v_cell_stencil_fused_face_1);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hhg
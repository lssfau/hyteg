
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsEdgeToVertexMacroCell3D.hpp"

namespace hhg {
namespace EdgeDoFToVertexDoF {
namespace generated {

static void apply_3D_macrocell_edgedof_to_vertexdof_replace_level_2(double const * RESTRICT const _data_edgeCellSrc_X, double const * RESTRICT const _data_edgeCellSrc_XY, double const * RESTRICT const _data_edgeCellSrc_XYZ, double const * RESTRICT const _data_edgeCellSrc_XZ, double const * RESTRICT const _data_edgeCellSrc_Y, double const * RESTRICT const _data_edgeCellSrc_YZ, double const * RESTRICT const _data_edgeCellSrc_Z, double * RESTRICT _data_vertexCellDst, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2vStencilMap)
{
   const double xi_1 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_2 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_3 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_5 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_6 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_7 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_8 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_9 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_10 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_11 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_12 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_13 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_14 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_15 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_16 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_17 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_18 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_19 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_20 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_21 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_22 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_23 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_24 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_25 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_26 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_27 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_28 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_29 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_30 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_31 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_32 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_33 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_34 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_35 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_36 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_37 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_38 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_39 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_40 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_41 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_42 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_43 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_44 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_45 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_46 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_47 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_48 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_49 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_50 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   for (int ctr_3 = 1; ctr_3 < 4; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 4; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 4; ctr_1 += 1)
         {
            const double xi_53 = xi_1*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6)) - 1];
            const double xi_64 = xi_2*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6)) - 1];
            const double xi_75 = xi_3*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6)) - 1];
            const double xi_86 = xi_4*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 5) + ((60) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            const double xi_97 = xi_5*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 4) + ((60) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))];
            const double xi_99 = xi_6*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((60) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            const double xi_100 = xi_7*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6)) - 1];
            const double xi_101 = xi_8*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 4) + ((120) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6)) - 1];
            const double xi_102 = xi_9*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 6) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6)) - 1];
            const double xi_54 = xi_10*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6)) - 1];
            const double xi_55 = xi_11*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            const double xi_56 = xi_12*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 4) + ((120) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))];
            const double xi_57 = xi_13*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 6) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6))];
            const double xi_58 = xi_14*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            const double xi_59 = xi_15*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 6) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6)) - 1];
            const double xi_60 = xi_16*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6)) - 1];
            const double xi_61 = xi_17*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 6) + ((120) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6)) - 1];
            const double xi_62 = xi_18*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            const double xi_63 = xi_19*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 6) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6))];
            const double xi_65 = xi_20*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            const double xi_66 = xi_21*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6)) - 1];
            const double xi_67 = xi_22*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 4) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6)) - 1];
            const double xi_68 = xi_23*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 6) + ((120) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6)) - 1];
            const double xi_69 = xi_24*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 5) + ((120) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6)) - 1];
            const double xi_70 = xi_25*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            const double xi_71 = xi_26*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 4) + ((120) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))];
            const double xi_72 = xi_27*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 6) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6))];
            const double xi_73 = xi_28*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            const double xi_74 = xi_29*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 6) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6)) - 1];
            const double xi_76 = xi_30*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6)) - 1];
            const double xi_77 = xi_31*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 6) + ((120) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6))];
            const double xi_78 = xi_32*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            const double xi_79 = xi_33*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 6) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6))];
            const double xi_80 = xi_34*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            const double xi_81 = xi_35*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 6) + ((120) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6)) + 1];
            const double xi_82 = xi_36*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6)) + 1];
            const double xi_83 = xi_37*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6)) - 1];
            const double xi_84 = xi_38*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            const double xi_85 = xi_39*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 4) + ((120) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 3)*(-ctr_3 + 4)*(-ctr_3 + 5)) / (6))];
            const double xi_87 = xi_40*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 6) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6))];
            const double xi_88 = xi_41*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            const double xi_89 = xi_42*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6)) + 1];
            const double xi_90 = xi_43*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6)) - 1];
            const double xi_91 = xi_44*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 6) + ((120) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6)) - 1];
            const double xi_92 = xi_45*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            const double xi_93 = xi_46*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 6) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6))];
            const double xi_94 = xi_47*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6))];
            const double xi_95 = xi_48*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 6) + ((120) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6))];
            const double xi_96 = xi_49*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 5) + ((120) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 4)*(-ctr_3 + 5)*(-ctr_3 + 6)) / (6)) + 1];
            const double xi_98 = xi_50*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 6) + ((120) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6)) + 1];
            _data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + 6) + ((210) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 5)*(-ctr_3 + 6)*(-ctr_3 + 7)) / (6))] = xi_100 + xi_101 + xi_102 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68 + xi_69 + xi_70 + xi_71 + xi_72 + xi_73 + xi_74 + xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80 + xi_81 + xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96 + xi_97 + xi_98 + xi_99;
         }
      }
   }
}

static void apply_3D_macrocell_edgedof_to_vertexdof_replace_level_3(double const * RESTRICT const _data_edgeCellSrc_X, double const * RESTRICT const _data_edgeCellSrc_XY, double const * RESTRICT const _data_edgeCellSrc_XYZ, double const * RESTRICT const _data_edgeCellSrc_XZ, double const * RESTRICT const _data_edgeCellSrc_Y, double const * RESTRICT const _data_edgeCellSrc_YZ, double const * RESTRICT const _data_edgeCellSrc_Z, double * RESTRICT _data_vertexCellDst, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2vStencilMap)
{
   const double xi_1 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_2 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_3 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_5 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_6 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_7 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_8 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_9 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_10 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_11 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_12 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_13 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_14 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_15 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_16 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_17 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_18 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_19 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_20 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_21 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_22 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_23 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_24 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_25 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_26 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_27 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_28 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_29 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_30 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_31 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_32 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_33 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_34 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_35 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_36 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_37 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_38 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_39 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_40 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_41 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_42 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_43 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_44 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_45 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_46 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_47 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_48 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_49 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_50 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   for (int ctr_3 = 1; ctr_3 < 8; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 8; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 8; ctr_1 += 1)
         {
            const double xi_53 = xi_1*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6)) - 1];
            const double xi_64 = xi_2*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6)) - 1];
            const double xi_75 = xi_3*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6)) - 1];
            const double xi_86 = xi_4*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 9) + ((504) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            const double xi_97 = xi_5*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 8) + ((504) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))];
            const double xi_99 = xi_6*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((504) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            const double xi_100 = xi_7*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6)) - 1];
            const double xi_101 = xi_8*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 8) + ((720) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6)) - 1];
            const double xi_102 = xi_9*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 10) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6)) - 1];
            const double xi_54 = xi_10*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6)) - 1];
            const double xi_55 = xi_11*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            const double xi_56 = xi_12*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 8) + ((720) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))];
            const double xi_57 = xi_13*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 10) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6))];
            const double xi_58 = xi_14*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            const double xi_59 = xi_15*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 10) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6)) - 1];
            const double xi_60 = xi_16*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6)) - 1];
            const double xi_61 = xi_17*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 10) + ((720) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6)) - 1];
            const double xi_62 = xi_18*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            const double xi_63 = xi_19*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 10) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6))];
            const double xi_65 = xi_20*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            const double xi_66 = xi_21*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6)) - 1];
            const double xi_67 = xi_22*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 8) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6)) - 1];
            const double xi_68 = xi_23*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 10) + ((720) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6)) - 1];
            const double xi_69 = xi_24*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 9) + ((720) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6)) - 1];
            const double xi_70 = xi_25*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            const double xi_71 = xi_26*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 8) + ((720) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))];
            const double xi_72 = xi_27*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 10) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6))];
            const double xi_73 = xi_28*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            const double xi_74 = xi_29*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 10) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6)) - 1];
            const double xi_76 = xi_30*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6)) - 1];
            const double xi_77 = xi_31*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 10) + ((720) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6))];
            const double xi_78 = xi_32*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            const double xi_79 = xi_33*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 10) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6))];
            const double xi_80 = xi_34*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            const double xi_81 = xi_35*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 10) + ((720) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6)) + 1];
            const double xi_82 = xi_36*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6)) + 1];
            const double xi_83 = xi_37*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6)) - 1];
            const double xi_84 = xi_38*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            const double xi_85 = xi_39*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 8) + ((720) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 7)*(-ctr_3 + 8)*(-ctr_3 + 9)) / (6))];
            const double xi_87 = xi_40*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 10) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6))];
            const double xi_88 = xi_41*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            const double xi_89 = xi_42*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6)) + 1];
            const double xi_90 = xi_43*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6)) - 1];
            const double xi_91 = xi_44*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 10) + ((720) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6)) - 1];
            const double xi_92 = xi_45*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            const double xi_93 = xi_46*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 10) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6))];
            const double xi_94 = xi_47*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6))];
            const double xi_95 = xi_48*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 10) + ((720) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6))];
            const double xi_96 = xi_49*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 9) + ((720) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 8)*(-ctr_3 + 9)*(-ctr_3 + 10)) / (6)) + 1];
            const double xi_98 = xi_50*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 10) + ((720) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6)) + 1];
            _data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + 10) + ((990) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 9)*(-ctr_3 + 10)*(-ctr_3 + 11)) / (6))] = xi_100 + xi_101 + xi_102 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68 + xi_69 + xi_70 + xi_71 + xi_72 + xi_73 + xi_74 + xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80 + xi_81 + xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96 + xi_97 + xi_98 + xi_99;
         }
      }
   }
}

static void apply_3D_macrocell_edgedof_to_vertexdof_replace_level_4(double const * RESTRICT const _data_edgeCellSrc_X, double const * RESTRICT const _data_edgeCellSrc_XY, double const * RESTRICT const _data_edgeCellSrc_XYZ, double const * RESTRICT const _data_edgeCellSrc_XZ, double const * RESTRICT const _data_edgeCellSrc_Y, double const * RESTRICT const _data_edgeCellSrc_YZ, double const * RESTRICT const _data_edgeCellSrc_Z, double * RESTRICT _data_vertexCellDst, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2vStencilMap)
{
   const double xi_1 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_2 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_3 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_5 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_6 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_7 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_8 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_9 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_10 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_11 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_12 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_13 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_14 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_15 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_16 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_17 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_18 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_19 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_20 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_21 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_22 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_23 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_24 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_25 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_26 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_27 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_28 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_29 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_30 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_31 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_32 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_33 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_34 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_35 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_36 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_37 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_38 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_39 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_40 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_41 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_42 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_43 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_44 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_45 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_46 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_47 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_48 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_49 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_50 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   for (int ctr_3 = 1; ctr_3 < 16; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 16; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 16; ctr_1 += 1)
         {
            const double xi_53 = xi_1*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6)) - 1];
            const double xi_64 = xi_2*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6)) - 1];
            const double xi_75 = xi_3*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6)) - 1];
            const double xi_86 = xi_4*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 17) + ((4080) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            const double xi_97 = xi_5*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 16) + ((4080) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))];
            const double xi_99 = xi_6*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            const double xi_100 = xi_7*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6)) - 1];
            const double xi_101 = xi_8*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 16) + ((4896) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6)) - 1];
            const double xi_102 = xi_9*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 18) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6)) - 1];
            const double xi_54 = xi_10*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6)) - 1];
            const double xi_55 = xi_11*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            const double xi_56 = xi_12*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 16) + ((4896) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))];
            const double xi_57 = xi_13*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 18) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6))];
            const double xi_58 = xi_14*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            const double xi_59 = xi_15*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 18) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6)) - 1];
            const double xi_60 = xi_16*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6)) - 1];
            const double xi_61 = xi_17*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 18) + ((4896) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6)) - 1];
            const double xi_62 = xi_18*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            const double xi_63 = xi_19*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 18) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6))];
            const double xi_65 = xi_20*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            const double xi_66 = xi_21*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6)) - 1];
            const double xi_67 = xi_22*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 16) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6)) - 1];
            const double xi_68 = xi_23*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 18) + ((4896) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6)) - 1];
            const double xi_69 = xi_24*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 17) + ((4896) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6)) - 1];
            const double xi_70 = xi_25*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            const double xi_71 = xi_26*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 16) + ((4896) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))];
            const double xi_72 = xi_27*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 18) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6))];
            const double xi_73 = xi_28*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            const double xi_74 = xi_29*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 18) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6)) - 1];
            const double xi_76 = xi_30*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6)) - 1];
            const double xi_77 = xi_31*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 18) + ((4896) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6))];
            const double xi_78 = xi_32*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            const double xi_79 = xi_33*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 18) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6))];
            const double xi_80 = xi_34*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            const double xi_81 = xi_35*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 18) + ((4896) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6)) + 1];
            const double xi_82 = xi_36*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6)) + 1];
            const double xi_83 = xi_37*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6)) - 1];
            const double xi_84 = xi_38*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            const double xi_85 = xi_39*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 16) + ((4896) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 15)*(-ctr_3 + 16)*(-ctr_3 + 17)) / (6))];
            const double xi_87 = xi_40*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 18) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6))];
            const double xi_88 = xi_41*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            const double xi_89 = xi_42*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6)) + 1];
            const double xi_90 = xi_43*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6)) - 1];
            const double xi_91 = xi_44*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 18) + ((4896) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6)) - 1];
            const double xi_92 = xi_45*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            const double xi_93 = xi_46*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 18) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6))];
            const double xi_94 = xi_47*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6))];
            const double xi_95 = xi_48*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 18) + ((4896) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6))];
            const double xi_96 = xi_49*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 17) + ((4896) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 16)*(-ctr_3 + 17)*(-ctr_3 + 18)) / (6)) + 1];
            const double xi_98 = xi_50*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 18) + ((4896) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6)) + 1];
            _data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + 18) + ((5814) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 17)*(-ctr_3 + 18)*(-ctr_3 + 19)) / (6))] = xi_100 + xi_101 + xi_102 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68 + xi_69 + xi_70 + xi_71 + xi_72 + xi_73 + xi_74 + xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80 + xi_81 + xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96 + xi_97 + xi_98 + xi_99;
         }
      }
   }
}

static void apply_3D_macrocell_edgedof_to_vertexdof_replace_level_5(double const * RESTRICT const _data_edgeCellSrc_X, double const * RESTRICT const _data_edgeCellSrc_XY, double const * RESTRICT const _data_edgeCellSrc_XYZ, double const * RESTRICT const _data_edgeCellSrc_XZ, double const * RESTRICT const _data_edgeCellSrc_Y, double const * RESTRICT const _data_edgeCellSrc_YZ, double const * RESTRICT const _data_edgeCellSrc_Z, double * RESTRICT _data_vertexCellDst, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2vStencilMap)
{
   const double xi_1 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_2 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_3 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_5 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_6 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_7 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_8 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_9 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_10 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_11 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_12 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_13 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_14 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_15 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_16 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_17 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_18 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_19 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_20 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_21 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_22 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_23 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_24 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_25 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_26 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_27 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_28 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_29 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_30 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_31 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_32 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_33 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_34 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_35 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_36 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_37 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_38 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_39 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_40 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_41 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_42 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_43 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_44 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_45 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_46 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_47 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_48 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_49 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_50 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   for (int ctr_3 = 1; ctr_3 < 32; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 32; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 32; ctr_1 += 1)
         {
            const double xi_53 = xi_1*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6)) - 1];
            const double xi_64 = xi_2*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6)) - 1];
            const double xi_75 = xi_3*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6)) - 1];
            const double xi_86 = xi_4*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 33) + ((32736) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            const double xi_97 = xi_5*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 32) + ((32736) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))];
            const double xi_99 = xi_6*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((32736) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            const double xi_100 = xi_7*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6)) - 1];
            const double xi_101 = xi_8*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 32) + ((35904) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6)) - 1];
            const double xi_102 = xi_9*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 34) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6)) - 1];
            const double xi_54 = xi_10*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6)) - 1];
            const double xi_55 = xi_11*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            const double xi_56 = xi_12*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 32) + ((35904) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))];
            const double xi_57 = xi_13*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 34) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6))];
            const double xi_58 = xi_14*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            const double xi_59 = xi_15*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 34) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6)) - 1];
            const double xi_60 = xi_16*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6)) - 1];
            const double xi_61 = xi_17*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 34) + ((35904) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6)) - 1];
            const double xi_62 = xi_18*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            const double xi_63 = xi_19*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 34) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6))];
            const double xi_65 = xi_20*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            const double xi_66 = xi_21*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6)) - 1];
            const double xi_67 = xi_22*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 32) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6)) - 1];
            const double xi_68 = xi_23*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 34) + ((35904) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6)) - 1];
            const double xi_69 = xi_24*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 33) + ((35904) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6)) - 1];
            const double xi_70 = xi_25*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            const double xi_71 = xi_26*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 32) + ((35904) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))];
            const double xi_72 = xi_27*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 34) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6))];
            const double xi_73 = xi_28*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            const double xi_74 = xi_29*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 34) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6)) - 1];
            const double xi_76 = xi_30*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6)) - 1];
            const double xi_77 = xi_31*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 34) + ((35904) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6))];
            const double xi_78 = xi_32*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            const double xi_79 = xi_33*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 34) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6))];
            const double xi_80 = xi_34*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            const double xi_81 = xi_35*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 34) + ((35904) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6)) + 1];
            const double xi_82 = xi_36*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6)) + 1];
            const double xi_83 = xi_37*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6)) - 1];
            const double xi_84 = xi_38*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            const double xi_85 = xi_39*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 32) + ((35904) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 31)*(-ctr_3 + 32)*(-ctr_3 + 33)) / (6))];
            const double xi_87 = xi_40*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 34) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6))];
            const double xi_88 = xi_41*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            const double xi_89 = xi_42*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6)) + 1];
            const double xi_90 = xi_43*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6)) - 1];
            const double xi_91 = xi_44*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 34) + ((35904) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6)) - 1];
            const double xi_92 = xi_45*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            const double xi_93 = xi_46*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 34) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6))];
            const double xi_94 = xi_47*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6))];
            const double xi_95 = xi_48*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 34) + ((35904) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6))];
            const double xi_96 = xi_49*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 33) + ((35904) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 32)*(-ctr_3 + 33)*(-ctr_3 + 34)) / (6)) + 1];
            const double xi_98 = xi_50*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 34) + ((35904) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6)) + 1];
            _data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + 34) + ((39270) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 33)*(-ctr_3 + 34)*(-ctr_3 + 35)) / (6))] = xi_100 + xi_101 + xi_102 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68 + xi_69 + xi_70 + xi_71 + xi_72 + xi_73 + xi_74 + xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80 + xi_81 + xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96 + xi_97 + xi_98 + xi_99;
         }
      }
   }
}

static void apply_3D_macrocell_edgedof_to_vertexdof_replace_level_6(double const * RESTRICT const _data_edgeCellSrc_X, double const * RESTRICT const _data_edgeCellSrc_XY, double const * RESTRICT const _data_edgeCellSrc_XYZ, double const * RESTRICT const _data_edgeCellSrc_XZ, double const * RESTRICT const _data_edgeCellSrc_Y, double const * RESTRICT const _data_edgeCellSrc_YZ, double const * RESTRICT const _data_edgeCellSrc_Z, double * RESTRICT _data_vertexCellDst, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2vStencilMap)
{
   const double xi_1 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_2 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_3 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_5 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_6 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_7 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_8 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_9 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_10 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_11 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_12 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_13 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_14 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_15 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_16 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_17 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_18 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_19 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_20 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_21 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_22 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_23 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_24 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_25 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_26 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_27 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_28 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_29 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_30 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_31 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_32 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_33 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_34 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_35 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_36 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_37 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_38 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_39 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_40 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_41 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_42 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_43 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_44 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_45 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_46 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_47 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_48 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_49 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_50 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   for (int ctr_3 = 1; ctr_3 < 64; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 64; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 64; ctr_1 += 1)
         {
            const double xi_53 = xi_1*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6)) - 1];
            const double xi_64 = xi_2*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6)) - 1];
            const double xi_75 = xi_3*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6)) - 1];
            const double xi_86 = xi_4*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 65) + ((262080) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            const double xi_97 = xi_5*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 64) + ((262080) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))];
            const double xi_99 = xi_6*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((262080) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            const double xi_100 = xi_7*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6)) - 1];
            const double xi_101 = xi_8*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 64) + ((274560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6)) - 1];
            const double xi_102 = xi_9*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 66) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6)) - 1];
            const double xi_54 = xi_10*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6)) - 1];
            const double xi_55 = xi_11*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            const double xi_56 = xi_12*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 64) + ((274560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))];
            const double xi_57 = xi_13*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 66) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6))];
            const double xi_58 = xi_14*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            const double xi_59 = xi_15*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 66) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6)) - 1];
            const double xi_60 = xi_16*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6)) - 1];
            const double xi_61 = xi_17*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 66) + ((274560) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6)) - 1];
            const double xi_62 = xi_18*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            const double xi_63 = xi_19*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 66) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6))];
            const double xi_65 = xi_20*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            const double xi_66 = xi_21*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6)) - 1];
            const double xi_67 = xi_22*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 64) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6)) - 1];
            const double xi_68 = xi_23*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 66) + ((274560) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6)) - 1];
            const double xi_69 = xi_24*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 65) + ((274560) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6)) - 1];
            const double xi_70 = xi_25*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            const double xi_71 = xi_26*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 64) + ((274560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))];
            const double xi_72 = xi_27*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 66) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6))];
            const double xi_73 = xi_28*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            const double xi_74 = xi_29*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 66) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6)) - 1];
            const double xi_76 = xi_30*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6)) - 1];
            const double xi_77 = xi_31*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 66) + ((274560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6))];
            const double xi_78 = xi_32*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            const double xi_79 = xi_33*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 66) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6))];
            const double xi_80 = xi_34*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            const double xi_81 = xi_35*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 66) + ((274560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6)) + 1];
            const double xi_82 = xi_36*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6)) + 1];
            const double xi_83 = xi_37*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6)) - 1];
            const double xi_84 = xi_38*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            const double xi_85 = xi_39*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 64) + ((274560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 63)*(-ctr_3 + 64)*(-ctr_3 + 65)) / (6))];
            const double xi_87 = xi_40*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 66) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6))];
            const double xi_88 = xi_41*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            const double xi_89 = xi_42*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6)) + 1];
            const double xi_90 = xi_43*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6)) - 1];
            const double xi_91 = xi_44*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 66) + ((274560) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6)) - 1];
            const double xi_92 = xi_45*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            const double xi_93 = xi_46*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 66) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6))];
            const double xi_94 = xi_47*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6))];
            const double xi_95 = xi_48*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 66) + ((274560) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6))];
            const double xi_96 = xi_49*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 65) + ((274560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 64)*(-ctr_3 + 65)*(-ctr_3 + 66)) / (6)) + 1];
            const double xi_98 = xi_50*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 66) + ((274560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6)) + 1];
            _data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + 66) + ((287430) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 65)*(-ctr_3 + 66)*(-ctr_3 + 67)) / (6))] = xi_100 + xi_101 + xi_102 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68 + xi_69 + xi_70 + xi_71 + xi_72 + xi_73 + xi_74 + xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80 + xi_81 + xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96 + xi_97 + xi_98 + xi_99;
         }
      }
   }
}

static void apply_3D_macrocell_edgedof_to_vertexdof_replace_level_7(double const * RESTRICT const _data_edgeCellSrc_X, double const * RESTRICT const _data_edgeCellSrc_XY, double const * RESTRICT const _data_edgeCellSrc_XYZ, double const * RESTRICT const _data_edgeCellSrc_XZ, double const * RESTRICT const _data_edgeCellSrc_Y, double const * RESTRICT const _data_edgeCellSrc_YZ, double const * RESTRICT const _data_edgeCellSrc_Z, double * RESTRICT _data_vertexCellDst, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2vStencilMap)
{
   const double xi_1 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_2 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_3 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_5 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_6 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_7 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_8 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_9 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_10 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_11 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_12 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_13 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_14 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_15 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_16 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_17 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_18 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_19 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_20 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_21 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_22 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_23 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_24 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_25 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_26 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_27 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_28 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_29 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_30 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_31 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_32 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_33 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_34 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_35 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_36 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_37 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_38 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_39 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_40 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_41 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_42 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_43 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_44 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_45 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_46 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_47 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_48 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_49 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_50 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   for (int ctr_3 = 1; ctr_3 < 128; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 128; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 128; ctr_1 += 1)
         {
            const double xi_53 = xi_1*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6)) - 1];
            const double xi_64 = xi_2*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6)) - 1];
            const double xi_75 = xi_3*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6)) - 1];
            const double xi_86 = xi_4*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 129) + ((2097024) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            const double xi_97 = xi_5*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 128) + ((2097024) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))];
            const double xi_99 = xi_6*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2097024) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            const double xi_100 = xi_7*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6)) - 1];
            const double xi_101 = xi_8*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 128) + ((2146560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6)) - 1];
            const double xi_102 = xi_9*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 130) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6)) - 1];
            const double xi_54 = xi_10*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6)) - 1];
            const double xi_55 = xi_11*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            const double xi_56 = xi_12*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 128) + ((2146560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))];
            const double xi_57 = xi_13*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 130) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6))];
            const double xi_58 = xi_14*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            const double xi_59 = xi_15*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 130) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6)) - 1];
            const double xi_60 = xi_16*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6)) - 1];
            const double xi_61 = xi_17*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 130) + ((2146560) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6)) - 1];
            const double xi_62 = xi_18*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            const double xi_63 = xi_19*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 130) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6))];
            const double xi_65 = xi_20*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            const double xi_66 = xi_21*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6)) - 1];
            const double xi_67 = xi_22*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 128) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6)) - 1];
            const double xi_68 = xi_23*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 130) + ((2146560) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6)) - 1];
            const double xi_69 = xi_24*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 129) + ((2146560) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6)) - 1];
            const double xi_70 = xi_25*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            const double xi_71 = xi_26*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 128) + ((2146560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))];
            const double xi_72 = xi_27*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 130) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6))];
            const double xi_73 = xi_28*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            const double xi_74 = xi_29*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 130) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6)) - 1];
            const double xi_76 = xi_30*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6)) - 1];
            const double xi_77 = xi_31*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 130) + ((2146560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6))];
            const double xi_78 = xi_32*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            const double xi_79 = xi_33*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 130) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6))];
            const double xi_80 = xi_34*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            const double xi_81 = xi_35*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 130) + ((2146560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6)) + 1];
            const double xi_82 = xi_36*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6)) + 1];
            const double xi_83 = xi_37*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6)) - 1];
            const double xi_84 = xi_38*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            const double xi_85 = xi_39*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 128) + ((2146560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 127)*(-ctr_3 + 128)*(-ctr_3 + 129)) / (6))];
            const double xi_87 = xi_40*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 130) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6))];
            const double xi_88 = xi_41*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            const double xi_89 = xi_42*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6)) + 1];
            const double xi_90 = xi_43*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6)) - 1];
            const double xi_91 = xi_44*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 130) + ((2146560) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6)) - 1];
            const double xi_92 = xi_45*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            const double xi_93 = xi_46*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 130) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6))];
            const double xi_94 = xi_47*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6))];
            const double xi_95 = xi_48*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 130) + ((2146560) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6))];
            const double xi_96 = xi_49*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 129) + ((2146560) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 128)*(-ctr_3 + 129)*(-ctr_3 + 130)) / (6)) + 1];
            const double xi_98 = xi_50*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 130) + ((2146560) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6)) + 1];
            _data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + 130) + ((2196870) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 129)*(-ctr_3 + 130)*(-ctr_3 + 131)) / (6))] = xi_100 + xi_101 + xi_102 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68 + xi_69 + xi_70 + xi_71 + xi_72 + xi_73 + xi_74 + xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80 + xi_81 + xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96 + xi_97 + xi_98 + xi_99;
         }
      }
   }
}

static void apply_3D_macrocell_edgedof_to_vertexdof_replace_level_8(double const * RESTRICT const _data_edgeCellSrc_X, double const * RESTRICT const _data_edgeCellSrc_XY, double const * RESTRICT const _data_edgeCellSrc_XYZ, double const * RESTRICT const _data_edgeCellSrc_XZ, double const * RESTRICT const _data_edgeCellSrc_Y, double const * RESTRICT const _data_edgeCellSrc_YZ, double const * RESTRICT const _data_edgeCellSrc_Z, double * RESTRICT _data_vertexCellDst, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2vStencilMap)
{
   const double xi_1 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_2 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_3 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_5 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_6 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_7 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_8 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_9 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_10 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_11 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_12 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_13 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_14 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_15 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_16 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_17 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_18 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_19 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_20 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_21 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_22 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_23 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_24 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_25 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_26 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_27 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_28 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_29 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_30 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_31 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_32 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_33 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_34 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_35 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_36 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_37 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_38 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_39 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_40 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_41 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_42 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_43 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_44 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_45 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_46 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_47 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_48 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_49 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_50 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   for (int ctr_3 = 1; ctr_3 < 256; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 256; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 256; ctr_1 += 1)
         {
            const double xi_53 = xi_1*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6)) - 1];
            const double xi_64 = xi_2*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6)) - 1];
            const double xi_75 = xi_3*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6)) - 1];
            const double xi_86 = xi_4*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 257) + ((16776960) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            const double xi_97 = xi_5*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 256) + ((16776960) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))];
            const double xi_99 = xi_6*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16776960) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            const double xi_100 = xi_7*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6)) - 1];
            const double xi_101 = xi_8*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 256) + ((16974336) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6)) - 1];
            const double xi_102 = xi_9*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 258) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6)) - 1];
            const double xi_54 = xi_10*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6)) - 1];
            const double xi_55 = xi_11*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            const double xi_56 = xi_12*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 256) + ((16974336) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))];
            const double xi_57 = xi_13*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 258) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6))];
            const double xi_58 = xi_14*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            const double xi_59 = xi_15*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 258) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6)) - 1];
            const double xi_60 = xi_16*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6)) - 1];
            const double xi_61 = xi_17*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 258) + ((16974336) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6)) - 1];
            const double xi_62 = xi_18*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            const double xi_63 = xi_19*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 258) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6))];
            const double xi_65 = xi_20*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            const double xi_66 = xi_21*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6)) - 1];
            const double xi_67 = xi_22*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 256) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6)) - 1];
            const double xi_68 = xi_23*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 258) + ((16974336) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6)) - 1];
            const double xi_69 = xi_24*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 257) + ((16974336) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6)) - 1];
            const double xi_70 = xi_25*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            const double xi_71 = xi_26*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 256) + ((16974336) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))];
            const double xi_72 = xi_27*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 258) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6))];
            const double xi_73 = xi_28*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            const double xi_74 = xi_29*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 258) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6)) - 1];
            const double xi_76 = xi_30*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6)) - 1];
            const double xi_77 = xi_31*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 258) + ((16974336) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6))];
            const double xi_78 = xi_32*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            const double xi_79 = xi_33*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 258) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6))];
            const double xi_80 = xi_34*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            const double xi_81 = xi_35*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 258) + ((16974336) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6)) + 1];
            const double xi_82 = xi_36*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6)) + 1];
            const double xi_83 = xi_37*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6)) - 1];
            const double xi_84 = xi_38*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            const double xi_85 = xi_39*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 256) + ((16974336) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 255)*(-ctr_3 + 256)*(-ctr_3 + 257)) / (6))];
            const double xi_87 = xi_40*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 258) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6))];
            const double xi_88 = xi_41*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            const double xi_89 = xi_42*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6)) + 1];
            const double xi_90 = xi_43*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6)) - 1];
            const double xi_91 = xi_44*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 258) + ((16974336) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6)) - 1];
            const double xi_92 = xi_45*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            const double xi_93 = xi_46*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 258) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6))];
            const double xi_94 = xi_47*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6))];
            const double xi_95 = xi_48*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 258) + ((16974336) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6))];
            const double xi_96 = xi_49*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 257) + ((16974336) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 256)*(-ctr_3 + 257)*(-ctr_3 + 258)) / (6)) + 1];
            const double xi_98 = xi_50*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 258) + ((16974336) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6)) + 1];
            _data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + 258) + ((17173254) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 257)*(-ctr_3 + 258)*(-ctr_3 + 259)) / (6))] = xi_100 + xi_101 + xi_102 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68 + xi_69 + xi_70 + xi_71 + xi_72 + xi_73 + xi_74 + xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80 + xi_81 + xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96 + xi_97 + xi_98 + xi_99;
         }
      }
   }
}

static void apply_3D_macrocell_edgedof_to_vertexdof_replace_level_9(double const * RESTRICT const _data_edgeCellSrc_X, double const * RESTRICT const _data_edgeCellSrc_XY, double const * RESTRICT const _data_edgeCellSrc_XYZ, double const * RESTRICT const _data_edgeCellSrc_XZ, double const * RESTRICT const _data_edgeCellSrc_Y, double const * RESTRICT const _data_edgeCellSrc_YZ, double const * RESTRICT const _data_edgeCellSrc_Z, double * RESTRICT _data_vertexCellDst, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2vStencilMap)
{
   const double xi_1 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_2 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_3 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_5 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_6 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_7 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_8 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_9 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_10 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_11 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_12 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_13 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_14 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_15 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_16 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_17 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_18 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_19 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_20 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_21 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_22 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_23 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_24 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_25 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_26 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_27 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_28 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_29 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_30 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_31 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_32 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_33 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_34 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_35 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_36 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_37 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_38 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_39 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_40 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_41 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_42 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_43 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_44 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_45 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_46 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_47 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_48 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_49 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_50 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   for (int ctr_3 = 1; ctr_3 < 512; ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + 512; ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + 512; ctr_1 += 1)
         {
            const double xi_53 = xi_1*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6)) - 1];
            const double xi_64 = xi_2*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6)) - 1];
            const double xi_75 = xi_3*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6)) - 1];
            const double xi_86 = xi_4*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 513) + ((134217216) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            const double xi_97 = xi_5*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 512) + ((134217216) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))];
            const double xi_99 = xi_6*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((134217216) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            const double xi_100 = xi_7*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6)) - 1];
            const double xi_101 = xi_8*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 512) + ((135005184) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6)) - 1];
            const double xi_102 = xi_9*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 514) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6)) - 1];
            const double xi_54 = xi_10*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6)) - 1];
            const double xi_55 = xi_11*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            const double xi_56 = xi_12*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 512) + ((135005184) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))];
            const double xi_57 = xi_13*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 514) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6))];
            const double xi_58 = xi_14*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            const double xi_59 = xi_15*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 514) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6)) - 1];
            const double xi_60 = xi_16*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6)) - 1];
            const double xi_61 = xi_17*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 514) + ((135005184) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6)) - 1];
            const double xi_62 = xi_18*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            const double xi_63 = xi_19*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 514) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6))];
            const double xi_65 = xi_20*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            const double xi_66 = xi_21*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6)) - 1];
            const double xi_67 = xi_22*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 512) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6)) - 1];
            const double xi_68 = xi_23*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 514) + ((135005184) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6)) - 1];
            const double xi_69 = xi_24*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 513) + ((135005184) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6)) - 1];
            const double xi_70 = xi_25*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            const double xi_71 = xi_26*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 512) + ((135005184) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))];
            const double xi_72 = xi_27*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 514) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6))];
            const double xi_73 = xi_28*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            const double xi_74 = xi_29*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 514) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6)) - 1];
            const double xi_76 = xi_30*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6)) - 1];
            const double xi_77 = xi_31*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 514) + ((135005184) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6))];
            const double xi_78 = xi_32*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            const double xi_79 = xi_33*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 514) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6))];
            const double xi_80 = xi_34*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            const double xi_81 = xi_35*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 514) + ((135005184) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6)) + 1];
            const double xi_82 = xi_36*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6)) + 1];
            const double xi_83 = xi_37*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6)) - 1];
            const double xi_84 = xi_38*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            const double xi_85 = xi_39*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 512) + ((135005184) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 511)*(-ctr_3 + 512)*(-ctr_3 + 513)) / (6))];
            const double xi_87 = xi_40*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 514) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6))];
            const double xi_88 = xi_41*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            const double xi_89 = xi_42*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6)) + 1];
            const double xi_90 = xi_43*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6)) - 1];
            const double xi_91 = xi_44*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 514) + ((135005184) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6)) - 1];
            const double xi_92 = xi_45*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            const double xi_93 = xi_46*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 514) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6))];
            const double xi_94 = xi_47*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6))];
            const double xi_95 = xi_48*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + 514) + ((135005184) / (6)) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6))];
            const double xi_96 = xi_49*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + 513) + ((135005184) / (6)) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + 512)*(-ctr_3 + 513)*(-ctr_3 + 514)) / (6)) + 1];
            const double xi_98 = xi_50*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + 514) + ((135005184) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6)) + 1];
            _data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + 514) + ((135796230) / (6)) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + 513)*(-ctr_3 + 514)*(-ctr_3 + 515)) / (6))] = xi_100 + xi_101 + xi_102 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68 + xi_69 + xi_70 + xi_71 + xi_72 + xi_73 + xi_74 + xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80 + xi_81 + xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96 + xi_97 + xi_98 + xi_99;
         }
      }
   }
}

static void apply_3D_macrocell_edgedof_to_vertexdof_replace_level_any(double const * RESTRICT const _data_edgeCellSrc_X, double const * RESTRICT const _data_edgeCellSrc_XY, double const * RESTRICT const _data_edgeCellSrc_XYZ, double const * RESTRICT const _data_edgeCellSrc_XZ, double const * RESTRICT const _data_edgeCellSrc_Y, double const * RESTRICT const _data_edgeCellSrc_YZ, double const * RESTRICT const _data_edgeCellSrc_Z, double * RESTRICT _data_vertexCellDst, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2vStencilMap, int64_t level)
{
   const double xi_1 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, -1, 0 }];
   const double xi_2 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, -1 }];
   const double xi_3 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_4 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, -1 }];
   const double xi_5 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_6 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_7 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 0 }];
   const double xi_8 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, -1, 1 }];
   const double xi_9 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, -1 }];
   const double xi_10 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_11 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 0 }];
   const double xi_12 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_13 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, -1 }];
   const double xi_14 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_15 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, -1 }];
   const double xi_16 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 0, 0 }];
   const double xi_17 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, -1 }];
   const double xi_18 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, -1, 0 }];
   const double xi_19 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, -1 }];
   const double xi_20 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_21 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 0 }];
   const double xi_22 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 0, 1 }];
   const double xi_23 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, -1 }];
   const double xi_24 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_25 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 0 }];
   const double xi_26 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, -1, 1 }];
   const double xi_27 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, -1 }];
   const double xi_28 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_29 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, -1 }];
   const double xi_30 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ -1, 0, 0 }];
   const double xi_31 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, -1 }];
   const double xi_32 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, -1, 0 }];
   const double xi_33 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_34 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_35 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, -1 }];
   const double xi_36 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_37 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ -1, 0, 0 }];
   const double xi_38 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 0 }];
   const double xi_39 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, -1, 1 }];
   const double xi_40 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, -1 }];
   const double xi_41 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_42 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Y][{ 1, -1, 0 }];
   const double xi_43 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 0 }];
   const double xi_44 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, -1 }];
   const double xi_45 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 0 }];
   const double xi_46 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, -1 }];
   const double xi_47 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_48 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_49 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 0 }];
   const double xi_50 = e2vStencilMap[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, -1 }];
   for (int ctr_3 = 1; ctr_3 < (1 << (level)); ctr_3 += 1)
   {
      for (int ctr_2 = 1; ctr_2 < -ctr_3 + (1 << (level)); ctr_2 += 1)
      {
         // cell (inner)
         for (int ctr_1 = 1; ctr_1 < -ctr_2 - ctr_3 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_53 = xi_1*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6)) - 1];
            const double xi_64 = xi_2*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6)) - 1];
            const double xi_75 = xi_3*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6)) - 1];
            const double xi_86 = xi_4*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
            const double xi_97 = xi_5*_data_edgeCellSrc_XYZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
            const double xi_99 = xi_6*_data_edgeCellSrc_XYZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) - 1)*((1 << (level)) + 1)*(1 << (level))) / (6))];
            const double xi_100 = xi_7*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_101 = xi_8*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_102 = xi_9*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_54 = xi_10*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_55 = xi_11*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_56 = xi_12*_data_edgeCellSrc_XY[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_57 = xi_13*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_58 = xi_14*_data_edgeCellSrc_XY[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_59 = xi_15*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_60 = xi_16*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_61 = xi_17*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_62 = xi_18*_data_edgeCellSrc_XZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_63 = xi_19*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_65 = xi_20*_data_edgeCellSrc_XZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_66 = xi_21*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_67 = xi_22*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_68 = xi_23*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_69 = xi_24*_data_edgeCellSrc_X[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_70 = xi_25*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_71 = xi_26*_data_edgeCellSrc_X[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_72 = xi_27*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_73 = xi_28*_data_edgeCellSrc_X[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_74 = xi_29*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_76 = xi_30*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_77 = xi_31*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_78 = xi_32*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_79 = xi_33*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_80 = xi_34*_data_edgeCellSrc_YZ[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_81 = xi_35*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) + 1];
            const double xi_82 = xi_36*_data_edgeCellSrc_YZ[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            const double xi_83 = xi_37*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_84 = xi_38*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_85 = xi_39*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level))) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) - 1)*(-ctr_3 + (1 << (level)) + 1)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_87 = xi_40*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_88 = xi_41*_data_edgeCellSrc_Y[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_89 = xi_42*_data_edgeCellSrc_Y[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            const double xi_90 = xi_43*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - 1];
            const double xi_91 = xi_44*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) - 1];
            const double xi_92 = xi_45*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_93 = xi_46*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_94 = xi_47*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6))];
            const double xi_95 = xi_48*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 + 1)*(-ctr_3 + (1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))];
            const double xi_96 = xi_49*_data_edgeCellSrc_Z[ctr_1 + (ctr_2 - 1)*(-ctr_3 + (1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) - (((-ctr_3 + (1 << (level)))*(-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)) / (6)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) + 1];
            const double xi_98 = xi_50*_data_edgeCellSrc_Z[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*(1 << (level))) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6)) + 1];
            _data_vertexCellDst[ctr_1 + ctr_2*(-ctr_3 + (1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*((1 << (level)) + 2)*((1 << (level)) + 3)) / (6)) - (((-ctr_3 + (1 << (level)) + 1)*(-ctr_3 + (1 << (level)) + 2)*(-ctr_3 + (1 << (level)) + 3)) / (6))] = xi_100 + xi_101 + xi_102 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68 + xi_69 + xi_70 + xi_71 + xi_72 + xi_73 + xi_74 + xi_75 + xi_76 + xi_77 + xi_78 + xi_79 + xi_80 + xi_81 + xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96 + xi_97 + xi_98 + xi_99;
         }
      }
   }
}


void apply_3D_macrocell_edgedof_to_vertexdof_replace(double const * RESTRICT const _data_edgeCellSrc_X, double const * RESTRICT const _data_edgeCellSrc_XY, double const * RESTRICT const _data_edgeCellSrc_XYZ, double const * RESTRICT const _data_edgeCellSrc_XZ, double const * RESTRICT const _data_edgeCellSrc_Y, double const * RESTRICT const _data_edgeCellSrc_YZ, double const * RESTRICT const _data_edgeCellSrc_Z, double * RESTRICT _data_vertexCellDst, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > e2vStencilMap, int64_t level)
{
    switch( level )
    {
    case 2:
        apply_3D_macrocell_edgedof_to_vertexdof_replace_level_2(_data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, _data_vertexCellDst, e2vStencilMap);
        break;
    case 3:
        apply_3D_macrocell_edgedof_to_vertexdof_replace_level_3(_data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, _data_vertexCellDst, e2vStencilMap);
        break;
    case 4:
        apply_3D_macrocell_edgedof_to_vertexdof_replace_level_4(_data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, _data_vertexCellDst, e2vStencilMap);
        break;
    case 5:
        apply_3D_macrocell_edgedof_to_vertexdof_replace_level_5(_data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, _data_vertexCellDst, e2vStencilMap);
        break;
    case 6:
        apply_3D_macrocell_edgedof_to_vertexdof_replace_level_6(_data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, _data_vertexCellDst, e2vStencilMap);
        break;
    case 7:
        apply_3D_macrocell_edgedof_to_vertexdof_replace_level_7(_data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, _data_vertexCellDst, e2vStencilMap);
        break;
    case 8:
        apply_3D_macrocell_edgedof_to_vertexdof_replace_level_8(_data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, _data_vertexCellDst, e2vStencilMap);
        break;
    case 9:
        apply_3D_macrocell_edgedof_to_vertexdof_replace_level_9(_data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, _data_vertexCellDst, e2vStencilMap);
        break;
    default:
        apply_3D_macrocell_edgedof_to_vertexdof_replace_level_any(_data_edgeCellSrc_X, _data_edgeCellSrc_XY, _data_edgeCellSrc_XYZ, _data_edgeCellSrc_XZ, _data_edgeCellSrc_Y, _data_edgeCellSrc_YZ, _data_edgeCellSrc_Z, _data_vertexCellDst, e2vStencilMap, level);
        break;
    }
}
    

} // namespace generated
} // namespace EdgeDoFToVertexDoF
} // namespace hhg

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "apply_3D_macroface_one_sided_vertexdof_to_edgedof_add_301.hpp"

namespace hhg {
namespace VertexDoFToEdgeDoF {
namespace generated {

static void apply_3D_macroface_one_sided_vertexdof_to_edgedof_add_301_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_vertexFaceSrc, double const * RESTRICT const _data_vertexFaceSrc_gl0, int32_t level, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > v2e_cell_stencil)
{
   const double xi_42 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_43 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 1 }];
   const double xi_44 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 1, -1 }];
   const double xi_45 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::X][{ 0, 1, 0 }];
   const double xi_46 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::X][{ 1, 0, -1 }];
   const double xi_47 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::X][{ 1, 0, 0 }];
   const double xi_101 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_102 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 1 }];
   const double xi_103 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 1, 0 }];
   const double xi_104 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][{ 1, 0, 0 }];
   const double xi_105 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][{ 1, 0, 1 }];
   const double xi_123 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 0, 1 }];
   const double xi_124 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, 0 }];
   const double xi_125 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_126 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 1 }];
   const double xi_127 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, 0 }];
   const double xi_128 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_56 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_50 = xi_42*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_51 = xi_43*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_52 = xi_44*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_53 = xi_45*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_54 = xi_46*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1];
            const double xi_55 = xi_47*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_114 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_120 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_109 = xi_42*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_110 = xi_43*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_111 = xi_44*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_112 = xi_45*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_113 = xi_46*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1];
            const double xi_108 = xi_47*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_115 = xi_101*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_116 = xi_102*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_117 = xi_103*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_118 = xi_104*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_119 = xi_105*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_115 + xi_116 + xi_117 + xi_118 + xi_119 + xi_120;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_71 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_66 = xi_101*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_67 = xi_102*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_68 = xi_103*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_69 = xi_104*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_70 = xi_105*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_66 + xi_67 + xi_68 + xi_69 + xi_70 + xi_71;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_143 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_150 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_137 = xi_123*_data_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_138 = xi_124*_data_vertexFaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_139 = xi_125*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_140 = xi_126*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_141 = xi_127*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_142 = xi_128*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_144 = xi_42*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_145 = xi_43*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_146 = xi_44*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_147 = xi_45*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_148 = xi_46*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1];
            const double xi_149 = xi_47*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_144 + xi_145 + xi_146 + xi_147 + xi_148 + xi_149 + xi_150;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_26 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_33 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_39 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = xi_123*_data_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_21 = xi_124*_data_vertexFaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_22 = xi_125*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_23 = xi_126*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_24 = xi_127*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_25 = xi_128*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_30 = xi_42*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_31 = xi_43*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_32 = xi_44*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_27 = xi_45*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_28 = xi_46*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1];
            const double xi_29 = xi_47*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_34 = xi_101*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_35 = xi_102*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_36 = xi_103*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_37 = xi_104*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_38 = xi_105*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32 + xi_33;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_172 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_178 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_166 = xi_123*_data_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_167 = xi_124*_data_vertexFaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_168 = xi_125*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_169 = xi_126*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_170 = xi_127*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_171 = xi_128*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_173 = xi_101*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_174 = xi_102*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_175 = xi_103*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_176 = xi_104*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_177 = xi_105*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_166 + xi_167 + xi_168 + xi_169 + xi_170 + xi_171 + xi_172;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_88 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_82 = xi_123*_data_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_83 = xi_124*_data_vertexFaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_84 = xi_125*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_85 = xi_126*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_86 = xi_127*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_87 = xi_128*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88;
         }
      }
   }
}


void apply_3D_macroface_one_sided_vertexdof_to_edgedof_add_301(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_vertexFaceSrc, double const * RESTRICT const _data_vertexFaceSrc_gl0, int32_t level, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > v2e_cell_stencil)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_vertexdof_to_edgedof_add_301_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_vertexFaceSrc, _data_vertexFaceSrc_gl0, level, v2e_cell_stencil);
        break;
    }
}
    

} // namespace generated
} // namespace VertexDoFToEdgeDoF
} // namespace hhg

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "apply_3D_macroface_one_sided_vertexdof_to_edgedof_add_320.hpp"

namespace hhg {
namespace VertexDoFToEdgeDoF {
namespace generated {

static void apply_3D_macroface_one_sided_vertexdof_to_edgedof_add_320_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_vertexFaceSrc, double const * RESTRICT const _data_vertexFaceSrc_gl0, int64_t level, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > v2e_cell_stencil)
{
   const double xi_42 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_43 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 1 }];
   const double xi_44 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 1, -1 }];
   const double xi_45 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 1, 0 }];
   const double xi_46 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Y][{ 1, 0, 0 }];
   const double xi_100 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, -1, 1 }];
   const double xi_101 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_102 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 1 }];
   const double xi_103 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, 0 }];
   const double xi_104 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, -1, 1 }];
   const double xi_105 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::Z][{ 1, 0, 0 }];
   const double xi_153 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_154 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 1 }];
   const double xi_155 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 1, 0 }];
   const double xi_156 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 1, 1 }];
   const double xi_157 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, 0 }];
   const double xi_158 = v2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, 1 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_54 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_49 = xi_42*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_50 = xi_43*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_51 = xi_44*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1];
            const double xi_52 = xi_45*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_53 = xi_46*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_113 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_120 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_108 = xi_42*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_109 = xi_43*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_110 = xi_44*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1];
            const double xi_111 = xi_45*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_112 = xi_46*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_114 = xi_100*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_115 = xi_101*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_116 = xi_102*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_117 = xi_103*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_118 = xi_104*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_119 = xi_105*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_114 + xi_115 + xi_116 + xi_117 + xi_118 + xi_119 + xi_120;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_71 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_65 = xi_100*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_66 = xi_101*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_67 = xi_102*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_68 = xi_103*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_69 = xi_104*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_70 = xi_105*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_65 + xi_66 + xi_67 + xi_68 + xi_69 + xi_70 + xi_71;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_172 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_178 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_166 = xi_153*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_167 = xi_154*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_168 = xi_155*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_169 = xi_156*_data_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_170 = xi_157*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_171 = xi_158*_data_vertexFaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
            const double xi_173 = xi_42*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_174 = xi_43*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_175 = xi_44*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1];
            const double xi_176 = xi_45*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_177 = xi_46*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_166 + xi_167 + xi_168 + xi_169 + xi_170 + xi_171 + xi_172;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_26 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_32 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_39 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = xi_153*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_21 = xi_154*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_22 = xi_155*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_23 = xi_156*_data_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_24 = xi_157*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_25 = xi_158*_data_vertexFaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
            const double xi_29 = xi_42*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_30 = xi_43*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_31 = xi_44*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1];
            const double xi_27 = xi_45*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_28 = xi_46*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_33 = xi_100*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_34 = xi_101*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_35 = xi_102*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_36 = xi_103*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_37 = xi_104*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_38 = xi_105*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_20 + xi_21 + xi_22 + xi_23 + xi_24 + xi_25 + xi_26;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_27 + xi_28 + xi_29 + xi_30 + xi_31 + xi_32;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_33 + xi_34 + xi_35 + xi_36 + xi_37 + xi_38 + xi_39;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_143 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_150 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_137 = xi_153*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_138 = xi_154*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_139 = xi_155*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_140 = xi_156*_data_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_141 = xi_157*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_142 = xi_158*_data_vertexFaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
            const double xi_144 = xi_100*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_145 = xi_101*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_146 = xi_102*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_147 = xi_103*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_148 = xi_104*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_149 = xi_105*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_144 + xi_145 + xi_146 + xi_147 + xi_148 + xi_149 + xi_150;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_88 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_82 = xi_153*_data_vertexFaceSrc[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_83 = xi_154*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_84 = xi_155*_data_vertexFaceSrc[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_85 = xi_156*_data_vertexFaceSrc[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_86 = xi_157*_data_vertexFaceSrc_gl0[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_87 = xi_158*_data_vertexFaceSrc_gl0[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88;
         }
      }
   }
}


void apply_3D_macroface_one_sided_vertexdof_to_edgedof_add_320(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_vertexFaceSrc, double const * RESTRICT const _data_vertexFaceSrc_gl0, int64_t level, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > v2e_cell_stencil)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_vertexdof_to_edgedof_add_320_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_vertexFaceSrc, _data_vertexFaceSrc_gl0, level, v2e_cell_stencil);
        break;
    }
}
    

} // namespace generated
} // namespace VertexDoFToEdgeDoF
} // namespace hhg
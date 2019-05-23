
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_321.hpp"

namespace hhg {
namespace edgedof {
namespace macroface {
namespace generated {

static void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_321_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > > e2e_cell_stencil, int32_t level)
{
   const double xi_71 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_72 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_73 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_74 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 1, -1 }];
   const double xi_75 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_76 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][hhg::edgedof::EdgeDoFOrientation::X][{ 0, 1, -1 }];
   const double xi_77 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_78 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_79 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, -1 }];
   const double xi_80 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_81 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_82 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XY][hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_165 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][hhg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_166 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][hhg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_167 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_168 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_169 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_170 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][hhg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_172 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::XZ][hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_245 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][hhg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_246 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_247 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][hhg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 1 }];
   const double xi_248 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][hhg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_249 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][hhg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, 0 }];
   const double xi_250 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][hhg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_251 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][hhg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_252 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][hhg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_253 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][hhg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_254 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][hhg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_255 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][hhg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, 0 }];
   const double xi_256 = e2e_cell_stencil[hhg::edgedof::EdgeDoFOrientation::YZ][hhg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_85 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_86 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_87 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_88 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_89 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_90 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_91 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_94 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_95 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_96 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_176 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_177 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_178 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_179 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_180 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_181 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_182 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_183 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_184 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_185 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_186 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_187 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_188 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_189 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_190 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_191 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_192 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_193 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = xi_171*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_195 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_196 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_176 + xi_177 + xi_178 + xi_179 + xi_180 + xi_181 + xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_118 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_110 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_111 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_112 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_113 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_114 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_115 = xi_171*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_116 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_117 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115 + xi_116 + xi_117 + xi_118;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_271 = xi_245*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_272 = xi_246*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_273 = xi_247*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_274 = xi_248*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_275 = xi_249*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_276 = xi_250*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_277 = xi_251*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_278 = xi_252*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_279 = xi_253*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_280 = xi_254*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_281 = xi_255*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_282 = xi_256*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_283 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_284 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_285 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_286 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_287 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_288 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_289 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_290 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_291 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_292 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_293 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_294 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_271 + xi_272 + xi_273 + xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_283 + xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293 + xi_294;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_36 = xi_245*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_40 = xi_246*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_41 = xi_247*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_42 = xi_248*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_43 = xi_249*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_44 = xi_250*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_45 = xi_251*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_46 = xi_252*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_47 = xi_253*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_37 = xi_254*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_38 = xi_255*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_39 = xi_256*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_48 = xi_71*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_49 = xi_72*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_50 = xi_73*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_51 = xi_74*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_52 = xi_75*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_53 = xi_76*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_54 = xi_77*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_55 = xi_78*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_56 = xi_79*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_57 = xi_80*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_58 = xi_81*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_59 = xi_82*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_60 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_61 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_62 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_63 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_64 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_65 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_66 = xi_171*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_67 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_68 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45 + xi_46 + xi_47;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58 + xi_59;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_222 = xi_245*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_223 = xi_246*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_224 = xi_247*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_225 = xi_248*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_226 = xi_249*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_227 = xi_250*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_228 = xi_251*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_229 = xi_252*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_230 = xi_253*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_231 = xi_254*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_232 = xi_255*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_233 = xi_256*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_234 = xi_165*_data_edgeFaceSrc_X[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_235 = xi_166*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
            const double xi_236 = xi_167*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_237 = xi_168*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_238 = xi_169*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_239 = xi_170*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_240 = xi_171*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_241 = xi_172*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_242 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_222 + xi_223 + xi_224 + xi_225 + xi_226 + xi_227 + xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_234 + xi_235 + xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_135 = xi_245*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 6*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_136 = xi_246*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_137 = xi_247*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_138 = xi_248*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_139 = xi_249*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_140 = xi_250*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_141 = xi_251*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_142 = xi_252*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 4*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_143 = xi_253*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_144 = xi_254*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 5*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            const double xi_145 = xi_255*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
            const double xi_146 = xi_256*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 3*((((1 << (level)) - 1)*(1 << (level))) / (2)) + 3*((((1 << (level)) + 1)*(1 << (level))) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140 + xi_141 + xi_142 + xi_143 + xi_144 + xi_145 + xi_146;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_321(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::edgedof::EdgeDoFOrientation, std::map< hhg::indexing::IndexIncrement, double > > > e2e_cell_stencil, int32_t level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl_321_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, e2e_cell_stencil, level);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace edgedof
} // namespace hhg
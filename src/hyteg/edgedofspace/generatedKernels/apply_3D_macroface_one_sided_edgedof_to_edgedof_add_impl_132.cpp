
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_132.hpp"

namespace hyteg {
namespace edgedof {
namespace macroface {
namespace generated {

static void apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_132_level_any(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * RESTRICT const _data_edgeFaceSrc_gl0_XY, double const * RESTRICT const _data_edgeFaceSrc_gl0_XYZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_XZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_Y, double const * RESTRICT const _data_edgeFaceSrc_gl0_YZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_Z, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::IndexIncrement, double > > > e2e_cell_stencil, int32_t level)
{
   const double xi_74 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ -1, 0, 0 }];
   const double xi_75 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 0 }];
   const double xi_76 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ -1, 0, 1 }];
   const double xi_77 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_78 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ -1, 1, 0 }];
   const double xi_79 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_80 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ -1, 1, 0 }];
   const double xi_81 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_82 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_83 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_84 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ -1, 1, 0 }];
   const double xi_85 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::YZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_171 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, 0, -1 }];
   const double xi_172 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_173 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_174 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 1, -1 }];
   const double xi_175 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_176 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 1, -1 }];
   const double xi_177 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, -1 }];
   const double xi_178 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_179 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, 0, -1 }];
   const double xi_180 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_181 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   const double xi_182 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XY][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 1, -1 }];
   const double xi_213 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XYZ][{ 0, -1, 0 }];
   const double xi_214 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, -1, 1 }];
   const double xi_215 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XY][{ 0, 0, 0 }];
   const double xi_216 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::XZ][{ 0, 0, 0 }];
   const double xi_217 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::X][{ 0, 0, 0 }];
   const double xi_218 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 0, 0, 0 }];
   const double xi_219 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::YZ][{ 1, -1, 0 }];
   const double xi_220 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Y][{ 0, 0, 0 }];
   const double xi_221 = e2e_cell_stencil[hyteg::edgedof::EdgeDoFOrientation::XZ][hyteg::edgedof::EdgeDoFOrientation::Z][{ 0, 0, 0 }];
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_100 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_88 = xi_74*_data_edgeFaceSrc_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_89 = xi_75*_data_edgeFaceSrc_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_90 = xi_76*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_91 = xi_77*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_92 = xi_78*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_93 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_94 = xi_80*_data_edgeFaceSrc_gl0_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_95 = xi_81*_data_edgeFaceSrc_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_96 = xi_82*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_97 = xi_83*_data_edgeFaceSrc_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_98 = xi_84*_data_edgeFaceSrc_gl0_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_99 = xi_85*_data_edgeFaceSrc_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_100 + xi_88 + xi_89 + xi_90 + xi_91 + xi_92 + xi_93 + xi_94 + xi_95 + xi_96 + xi_97 + xi_98 + xi_99;
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_197 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_210 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_185 = xi_74*_data_edgeFaceSrc_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_186 = xi_75*_data_edgeFaceSrc_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_187 = xi_76*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_188 = xi_77*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_189 = xi_78*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_190 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_191 = xi_80*_data_edgeFaceSrc_gl0_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_192 = xi_81*_data_edgeFaceSrc_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_193 = xi_82*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = xi_83*_data_edgeFaceSrc_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_195 = xi_84*_data_edgeFaceSrc_gl0_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_196 = xi_85*_data_edgeFaceSrc_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_198 = xi_171*_data_edgeFaceSrc_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_199 = xi_172*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_200 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_201 = xi_174*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_202 = xi_175*_data_edgeFaceSrc_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_203 = xi_176*_data_edgeFaceSrc_gl0_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_204 = xi_177*_data_edgeFaceSrc_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_205 = xi_178*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_206 = xi_179*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_207 = xi_180*_data_edgeFaceSrc_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_208 = xi_181*_data_edgeFaceSrc_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_209 = xi_182*_data_edgeFaceSrc_gl0_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190 + xi_191 + xi_192 + xi_193 + xi_194 + xi_195 + xi_196 + xi_197;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_198 + xi_199 + xi_200 + xi_201 + xi_202 + xi_203 + xi_204 + xi_205 + xi_206 + xi_207 + xi_208 + xi_209 + xi_210;
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
         {
            const double xi_129 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_117 = xi_171*_data_edgeFaceSrc_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_118 = xi_172*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_119 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_120 = xi_174*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_121 = xi_175*_data_edgeFaceSrc_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_122 = xi_176*_data_edgeFaceSrc_gl0_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_123 = xi_177*_data_edgeFaceSrc_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_124 = xi_178*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_125 = xi_179*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_126 = xi_180*_data_edgeFaceSrc_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_127 = xi_181*_data_edgeFaceSrc_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_128 = xi_182*_data_edgeFaceSrc_gl0_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_117 + xi_118 + xi_119 + xi_120 + xi_121 + xi_122 + xi_123 + xi_124 + xi_125 + xi_126 + xi_127 + xi_128 + xi_129;
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_245 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_258 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_236 = xi_213*_data_edgeFaceSrc_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))];
            const double xi_237 = xi_214*_data_edgeFaceSrc_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_238 = xi_215*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_239 = xi_216*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_240 = xi_217*_data_edgeFaceSrc_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_241 = xi_218*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_242 = xi_219*_data_edgeFaceSrc_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
            const double xi_243 = xi_220*_data_edgeFaceSrc_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_244 = xi_221*_data_edgeFaceSrc_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_246 = xi_74*_data_edgeFaceSrc_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_247 = xi_75*_data_edgeFaceSrc_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_248 = xi_76*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_249 = xi_77*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_250 = xi_78*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_251 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_252 = xi_80*_data_edgeFaceSrc_gl0_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_253 = xi_81*_data_edgeFaceSrc_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_254 = xi_82*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_255 = xi_83*_data_edgeFaceSrc_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_256 = xi_84*_data_edgeFaceSrc_gl0_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_257 = xi_85*_data_edgeFaceSrc_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_236 + xi_237 + xi_238 + xi_239 + xi_240 + xi_241 + xi_242 + xi_243 + xi_244 + xi_245;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_246 + xi_247 + xi_248 + xi_249 + xi_250 + xi_251 + xi_252 + xi_253 + xi_254 + xi_255 + xi_256 + xi_257 + xi_258;
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
         {
            const double xi_45 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_58 = _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_71 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_36 = xi_213*_data_edgeFaceSrc_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))];
            const double xi_37 = xi_214*_data_edgeFaceSrc_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_38 = xi_215*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_39 = xi_216*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_40 = xi_217*_data_edgeFaceSrc_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_41 = xi_218*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_42 = xi_219*_data_edgeFaceSrc_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
            const double xi_43 = xi_220*_data_edgeFaceSrc_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_44 = xi_221*_data_edgeFaceSrc_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_46 = xi_74*_data_edgeFaceSrc_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_47 = xi_75*_data_edgeFaceSrc_gl0_Y[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_48 = xi_76*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
            const double xi_49 = xi_77*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_50 = xi_78*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_51 = xi_79*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_52 = xi_80*_data_edgeFaceSrc_gl0_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_53 = xi_81*_data_edgeFaceSrc_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_54 = xi_82*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_55 = xi_83*_data_edgeFaceSrc_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_56 = xi_84*_data_edgeFaceSrc_gl0_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
            const double xi_57 = xi_85*_data_edgeFaceSrc_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_59 = xi_171*_data_edgeFaceSrc_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_60 = xi_172*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_61 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_174*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_63 = xi_175*_data_edgeFaceSrc_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_64 = xi_176*_data_edgeFaceSrc_gl0_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_65 = xi_177*_data_edgeFaceSrc_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_66 = xi_178*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_67 = xi_179*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_68 = xi_180*_data_edgeFaceSrc_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_69 = xi_181*_data_edgeFaceSrc_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_70 = xi_182*_data_edgeFaceSrc_gl0_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43 + xi_44 + xi_45;
            _data_edgeFaceDst_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54 + xi_55 + xi_56 + xi_57 + xi_58;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65 + xi_66 + xi_67 + xi_68 + xi_69 + xi_70 + xi_71;
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
         {
            const double xi_293 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_306 = _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_284 = xi_213*_data_edgeFaceSrc_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))];
            const double xi_285 = xi_214*_data_edgeFaceSrc_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_286 = xi_215*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_287 = xi_216*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_288 = xi_217*_data_edgeFaceSrc_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_289 = xi_218*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_290 = xi_219*_data_edgeFaceSrc_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
            const double xi_291 = xi_220*_data_edgeFaceSrc_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_292 = xi_221*_data_edgeFaceSrc_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_294 = xi_171*_data_edgeFaceSrc_gl0_XYZ[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_295 = xi_172*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_296 = xi_173*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_297 = xi_174*_data_edgeFaceSrc_X[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_298 = xi_175*_data_edgeFaceSrc_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_299 = xi_176*_data_edgeFaceSrc_gl0_Z[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            const double xi_300 = xi_177*_data_edgeFaceSrc_gl0_XY[ctr_1 + ctr_2*(1 << (level)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_301 = xi_178*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_302 = xi_179*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
            const double xi_303 = xi_180*_data_edgeFaceSrc_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_304 = xi_181*_data_edgeFaceSrc_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_305 = xi_182*_data_edgeFaceSrc_gl0_XZ[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_284 + xi_285 + xi_286 + xi_287 + xi_288 + xi_289 + xi_290 + xi_291 + xi_292 + xi_293;
            _data_edgeFaceDst_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_294 + xi_295 + xi_296 + xi_297 + xi_298 + xi_299 + xi_300 + xi_301 + xi_302 + xi_303 + xi_304 + xi_305 + xi_306;
         }
      }
      for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_152 = _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_143 = xi_213*_data_edgeFaceSrc_gl0_XYZ[ctr_1 + (ctr_2 - 1)*(1 << (level)) - ((ctr_2*(ctr_2 - 1)) / (2))];
            const double xi_144 = xi_214*_data_edgeFaceSrc_Y[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
            const double xi_145 = xi_215*_data_edgeFaceSrc_Y[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_146 = xi_216*_data_edgeFaceSrc_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_147 = xi_217*_data_edgeFaceSrc_gl0_Z[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_148 = xi_218*_data_edgeFaceSrc_XY[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_149 = xi_219*_data_edgeFaceSrc_XY[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2))];
            const double xi_150 = xi_220*_data_edgeFaceSrc_gl0_YZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_151 = xi_221*_data_edgeFaceSrc_gl0_XZ[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFaceDst_X[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = xi_143 + xi_144 + xi_145 + xi_146 + xi_147 + xi_148 + xi_149 + xi_150 + xi_151 + xi_152;
         }
      }
   }
}


void apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_132(double * RESTRICT _data_edgeFaceDst_X, double * RESTRICT _data_edgeFaceDst_XY, double * RESTRICT _data_edgeFaceDst_Y, double const * RESTRICT const _data_edgeFaceSrc_X, double const * RESTRICT const _data_edgeFaceSrc_XY, double const * RESTRICT const _data_edgeFaceSrc_Y, double const * RESTRICT const _data_edgeFaceSrc_gl0_XY, double const * RESTRICT const _data_edgeFaceSrc_gl0_XYZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_XZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_Y, double const * RESTRICT const _data_edgeFaceSrc_gl0_YZ, double const * RESTRICT const _data_edgeFaceSrc_gl0_Z, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::edgedof::EdgeDoFOrientation, std::map< hyteg::indexing::IndexIncrement, double > > > e2e_cell_stencil, int32_t level)
{
    switch( level )
    {

    default:
        apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl_132_level_any(_data_edgeFaceDst_X, _data_edgeFaceDst_XY, _data_edgeFaceDst_Y, _data_edgeFaceSrc_X, _data_edgeFaceSrc_XY, _data_edgeFaceSrc_Y, _data_edgeFaceSrc_gl0_XY, _data_edgeFaceSrc_gl0_XYZ, _data_edgeFaceSrc_gl0_XZ, _data_edgeFaceSrc_gl0_Y, _data_edgeFaceSrc_gl0_YZ, _data_edgeFaceSrc_gl0_Z, e2e_cell_stencil, level);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace edgedof
} // namespace hyteg
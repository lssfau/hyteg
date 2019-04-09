
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsP2MacroFace2D.hpp"

namespace hhg {
namespace P2 {
namespace macroface {
namespace generated {

static void sor_2D_macroface_P2_update_edgedofs_level_2(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceRhs, double const * const _data_edge_stencil_at_edge_x, double const * const _data_edge_stencil_at_edge_xy, double const * const _data_edge_stencil_at_edge_y, double * RESTRICT _data_vertexFaceDst, double const * const _data_vertex_stencil_at_edge_x, double const * const _data_vertex_stencil_at_edge_xy, double const * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_80 = 1.0;
   const double xi_81 = -relax;
   const double xi_67 = _data_edge_stencil_at_edge_xy[0];
   const double xi_78 = 1 / (xi_67);
   const double xi_68 = _data_edge_stencil_at_edge_xy[1];
   const double xi_69 = _data_edge_stencil_at_edge_xy[4];
   const double xi_70 = _data_edge_stencil_at_edge_xy[3];
   const double xi_71 = _data_edge_stencil_at_edge_xy[2];
   const double xi_72 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_73 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_74 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_155 = _data_edge_stencil_at_edge_y[0];
   const double xi_167 = 1 / (xi_155);
   const double xi_156 = _data_edge_stencil_at_edge_y[4];
   const double xi_157 = _data_edge_stencil_at_edge_y[3];
   const double xi_158 = _data_edge_stencil_at_edge_y[1];
   const double xi_159 = _data_edge_stencil_at_edge_y[2];
   const double xi_160 = _data_vertex_stencil_at_edge_y[3];
   const double xi_161 = _data_vertex_stencil_at_edge_y[0];
   const double xi_162 = _data_vertex_stencil_at_edge_y[2];
   const double xi_163 = _data_vertex_stencil_at_edge_y[1];
   const double xi_238 = _data_edge_stencil_at_edge_x[0];
   const double xi_258 = 1 / (xi_238);
   const double xi_239 = _data_edge_stencil_at_edge_x[1];
   const double xi_240 = _data_edge_stencil_at_edge_x[3];
   const double xi_241 = _data_edge_stencil_at_edge_x[4];
   const double xi_242 = _data_edge_stencil_at_edge_x[2];
   const double xi_243 = _data_vertex_stencil_at_edge_x[1];
   const double xi_244 = _data_vertex_stencil_at_edge_x[3];
   const double xi_245 = _data_vertex_stencil_at_edge_x[2];
   const double xi_246 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_90 = _data_edgeFaceRhs[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_82 = -xi_68*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = -xi_69*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_84 = -xi_70*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5];
         const double xi_85 = -xi_71*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_86 = -xi_72*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = -xi_73*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6];
         const double xi_88 = -xi_74*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_89 = -xi_75*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 7];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
      {
         const double xi_179 = _data_edgeFaceRhs[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_190 = _data_edgeFaceRhs[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_171 = -xi_68*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_172 = -xi_69*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = -xi_70*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5];
         const double xi_174 = -xi_71*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = -xi_72*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_176 = -xi_73*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6];
         const double xi_177 = -xi_74*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = -xi_75*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 7];
         const double xi_182 = -xi_156*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_183 = -xi_157*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4];
         const double xi_184 = -xi_158*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_185 = -xi_159*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_186 = -xi_160*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5];
         const double xi_187 = -xi_161*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_188 = -xi_162*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6];
         const double xi_189 = -xi_163*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178 + xi_179) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 3; ctr_1 < 4; ctr_1 += 1)
      {
         const double xi_115 = _data_edgeFaceRhs[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = -xi_156*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_109 = -xi_157*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4];
         const double xi_110 = -xi_158*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = -xi_159*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_112 = -xi_160*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5];
         const double xi_113 = -xi_161*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_114 = -xi_162*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6];
         const double xi_107 = -xi_163*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 3; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_271 = _data_edgeFaceRhs[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_282 = _data_edgeFaceRhs[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_263 = -xi_239*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_264 = -xi_240*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_265 = -xi_241*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_266 = -xi_242*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4];
         const double xi_267 = -xi_243*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_268 = -xi_244*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6];
         const double xi_269 = -xi_245*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_270 = -xi_246*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_274 = -xi_68*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_275 = -xi_69*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_276 = -xi_70*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5];
         const double xi_277 = -xi_71*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_278 = -xi_72*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_279 = -xi_73*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6];
         const double xi_280 = -xi_74*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_281 = -xi_75*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 7];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_263 + xi_264 + xi_265 + xi_266 + xi_267 + xi_268 + xi_269 + xi_270 + xi_271) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 3; ctr_1 += 1)
      {
         const double xi_43 = _data_edgeFaceRhs[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = _data_edgeFaceRhs[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = _data_edgeFaceRhs[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = -xi_239*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_36 = -xi_240*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_37 = -xi_241*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_242*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4];
         const double xi_39 = -xi_243*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_40 = -xi_244*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6];
         const double xi_41 = -xi_245*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_42 = -xi_246*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = -xi_68*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = -xi_69*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_48 = -xi_70*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5];
         const double xi_49 = -xi_71*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_50 = -xi_72*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_51 = -xi_73*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6];
         const double xi_52 = -xi_74*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_53 = -xi_75*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 7];
         const double xi_57 = -xi_156*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_58 = -xi_157*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4];
         const double xi_59 = -xi_158*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = -xi_159*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_61 = -xi_160*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5];
         const double xi_62 = -xi_161*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = -xi_162*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6];
         const double xi_64 = -xi_163*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 3; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
      {
         const double xi_225 = _data_edgeFaceRhs[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_236 = _data_edgeFaceRhs[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_217 = -xi_239*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_218 = -xi_240*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_219 = -xi_241*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_220 = -xi_242*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4];
         const double xi_221 = -xi_243*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_222 = -xi_244*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6];
         const double xi_223 = -xi_245*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_224 = -xi_246*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_228 = -xi_156*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_229 = -xi_157*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4];
         const double xi_230 = -xi_158*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_231 = -xi_159*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_232 = -xi_160*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5];
         const double xi_233 = -xi_161*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_234 = -xi_162*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6];
         const double xi_235 = -xi_163*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_217 + xi_218 + xi_219 + xi_220 + xi_221 + xi_222 + xi_223 + xi_224 + xi_225) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 3; ctr_2 < 4; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_140 = _data_edgeFaceRhs[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_132 = -xi_239*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_133 = -xi_240*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_134 = -xi_241*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = -xi_242*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4];
         const double xi_136 = -xi_243*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_137 = -xi_244*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6];
         const double xi_138 = -xi_245*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5];
         const double xi_139 = -xi_246*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_132 + xi_133 + xi_134 + xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_3(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceRhs, double const * const _data_edge_stencil_at_edge_x, double const * const _data_edge_stencil_at_edge_xy, double const * const _data_edge_stencil_at_edge_y, double * RESTRICT _data_vertexFaceDst, double const * const _data_vertex_stencil_at_edge_x, double const * const _data_vertex_stencil_at_edge_xy, double const * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_80 = 1.0;
   const double xi_81 = -relax;
   const double xi_67 = _data_edge_stencil_at_edge_xy[0];
   const double xi_78 = 1 / (xi_67);
   const double xi_68 = _data_edge_stencil_at_edge_xy[1];
   const double xi_69 = _data_edge_stencil_at_edge_xy[4];
   const double xi_70 = _data_edge_stencil_at_edge_xy[3];
   const double xi_71 = _data_edge_stencil_at_edge_xy[2];
   const double xi_72 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_73 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_74 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_155 = _data_edge_stencil_at_edge_y[0];
   const double xi_167 = 1 / (xi_155);
   const double xi_156 = _data_edge_stencil_at_edge_y[4];
   const double xi_157 = _data_edge_stencil_at_edge_y[3];
   const double xi_158 = _data_edge_stencil_at_edge_y[1];
   const double xi_159 = _data_edge_stencil_at_edge_y[2];
   const double xi_160 = _data_vertex_stencil_at_edge_y[3];
   const double xi_161 = _data_vertex_stencil_at_edge_y[0];
   const double xi_162 = _data_vertex_stencil_at_edge_y[2];
   const double xi_163 = _data_vertex_stencil_at_edge_y[1];
   const double xi_238 = _data_edge_stencil_at_edge_x[0];
   const double xi_258 = 1 / (xi_238);
   const double xi_239 = _data_edge_stencil_at_edge_x[1];
   const double xi_240 = _data_edge_stencil_at_edge_x[3];
   const double xi_241 = _data_edge_stencil_at_edge_x[4];
   const double xi_242 = _data_edge_stencil_at_edge_x[2];
   const double xi_243 = _data_vertex_stencil_at_edge_x[1];
   const double xi_244 = _data_vertex_stencil_at_edge_x[3];
   const double xi_245 = _data_vertex_stencil_at_edge_x[2];
   const double xi_246 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_90 = _data_edgeFaceRhs[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_82 = -xi_68*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = -xi_69*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_84 = -xi_70*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9];
         const double xi_85 = -xi_71*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_86 = -xi_72*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = -xi_73*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10];
         const double xi_88 = -xi_74*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_89 = -xi_75*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 11];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 7; ctr_1 += 1)
      {
         const double xi_179 = _data_edgeFaceRhs[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_190 = _data_edgeFaceRhs[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_171 = -xi_68*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_172 = -xi_69*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = -xi_70*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9];
         const double xi_174 = -xi_71*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = -xi_72*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_176 = -xi_73*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10];
         const double xi_177 = -xi_74*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = -xi_75*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 11];
         const double xi_182 = -xi_156*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_183 = -xi_157*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8];
         const double xi_184 = -xi_158*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_185 = -xi_159*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_186 = -xi_160*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9];
         const double xi_187 = -xi_161*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_188 = -xi_162*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10];
         const double xi_189 = -xi_163*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178 + xi_179) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 7; ctr_1 < 8; ctr_1 += 1)
      {
         const double xi_115 = _data_edgeFaceRhs[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = -xi_156*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_109 = -xi_157*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8];
         const double xi_110 = -xi_158*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = -xi_159*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_112 = -xi_160*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9];
         const double xi_113 = -xi_161*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_114 = -xi_162*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10];
         const double xi_107 = -xi_163*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 7; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_271 = _data_edgeFaceRhs[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_282 = _data_edgeFaceRhs[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_263 = -xi_239*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_264 = -xi_240*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_265 = -xi_241*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_266 = -xi_242*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8];
         const double xi_267 = -xi_243*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_268 = -xi_244*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10];
         const double xi_269 = -xi_245*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_270 = -xi_246*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_274 = -xi_68*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_275 = -xi_69*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_276 = -xi_70*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9];
         const double xi_277 = -xi_71*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_278 = -xi_72*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_279 = -xi_73*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10];
         const double xi_280 = -xi_74*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_281 = -xi_75*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 11];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_263 + xi_264 + xi_265 + xi_266 + xi_267 + xi_268 + xi_269 + xi_270 + xi_271) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 7; ctr_1 += 1)
      {
         const double xi_43 = _data_edgeFaceRhs[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = _data_edgeFaceRhs[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = _data_edgeFaceRhs[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = -xi_239*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_36 = -xi_240*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_37 = -xi_241*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_242*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8];
         const double xi_39 = -xi_243*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_40 = -xi_244*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10];
         const double xi_41 = -xi_245*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_42 = -xi_246*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = -xi_68*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = -xi_69*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_48 = -xi_70*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9];
         const double xi_49 = -xi_71*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_50 = -xi_72*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_51 = -xi_73*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10];
         const double xi_52 = -xi_74*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_53 = -xi_75*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 11];
         const double xi_57 = -xi_156*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_58 = -xi_157*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8];
         const double xi_59 = -xi_158*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = -xi_159*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_61 = -xi_160*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9];
         const double xi_62 = -xi_161*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = -xi_162*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10];
         const double xi_64 = -xi_163*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 7; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
      {
         const double xi_225 = _data_edgeFaceRhs[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_236 = _data_edgeFaceRhs[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_217 = -xi_239*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_218 = -xi_240*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_219 = -xi_241*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_220 = -xi_242*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8];
         const double xi_221 = -xi_243*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_222 = -xi_244*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10];
         const double xi_223 = -xi_245*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_224 = -xi_246*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_228 = -xi_156*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_229 = -xi_157*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8];
         const double xi_230 = -xi_158*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_231 = -xi_159*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_232 = -xi_160*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9];
         const double xi_233 = -xi_161*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_234 = -xi_162*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10];
         const double xi_235 = -xi_163*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_217 + xi_218 + xi_219 + xi_220 + xi_221 + xi_222 + xi_223 + xi_224 + xi_225) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 7; ctr_2 < 8; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_140 = _data_edgeFaceRhs[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_132 = -xi_239*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_133 = -xi_240*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_134 = -xi_241*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = -xi_242*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8];
         const double xi_136 = -xi_243*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_137 = -xi_244*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10];
         const double xi_138 = -xi_245*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9];
         const double xi_139 = -xi_246*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_132 + xi_133 + xi_134 + xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_4(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceRhs, double const * const _data_edge_stencil_at_edge_x, double const * const _data_edge_stencil_at_edge_xy, double const * const _data_edge_stencil_at_edge_y, double * RESTRICT _data_vertexFaceDst, double const * const _data_vertex_stencil_at_edge_x, double const * const _data_vertex_stencil_at_edge_xy, double const * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_80 = 1.0;
   const double xi_81 = -relax;
   const double xi_67 = _data_edge_stencil_at_edge_xy[0];
   const double xi_78 = 1 / (xi_67);
   const double xi_68 = _data_edge_stencil_at_edge_xy[1];
   const double xi_69 = _data_edge_stencil_at_edge_xy[4];
   const double xi_70 = _data_edge_stencil_at_edge_xy[3];
   const double xi_71 = _data_edge_stencil_at_edge_xy[2];
   const double xi_72 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_73 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_74 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_155 = _data_edge_stencil_at_edge_y[0];
   const double xi_167 = 1 / (xi_155);
   const double xi_156 = _data_edge_stencil_at_edge_y[4];
   const double xi_157 = _data_edge_stencil_at_edge_y[3];
   const double xi_158 = _data_edge_stencil_at_edge_y[1];
   const double xi_159 = _data_edge_stencil_at_edge_y[2];
   const double xi_160 = _data_vertex_stencil_at_edge_y[3];
   const double xi_161 = _data_vertex_stencil_at_edge_y[0];
   const double xi_162 = _data_vertex_stencil_at_edge_y[2];
   const double xi_163 = _data_vertex_stencil_at_edge_y[1];
   const double xi_238 = _data_edge_stencil_at_edge_x[0];
   const double xi_258 = 1 / (xi_238);
   const double xi_239 = _data_edge_stencil_at_edge_x[1];
   const double xi_240 = _data_edge_stencil_at_edge_x[3];
   const double xi_241 = _data_edge_stencil_at_edge_x[4];
   const double xi_242 = _data_edge_stencil_at_edge_x[2];
   const double xi_243 = _data_vertex_stencil_at_edge_x[1];
   const double xi_244 = _data_vertex_stencil_at_edge_x[3];
   const double xi_245 = _data_vertex_stencil_at_edge_x[2];
   const double xi_246 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_90 = _data_edgeFaceRhs[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_82 = -xi_68*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = -xi_69*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_84 = -xi_70*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17];
         const double xi_85 = -xi_71*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_86 = -xi_72*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = -xi_73*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18];
         const double xi_88 = -xi_74*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_89 = -xi_75*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 19];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
      {
         const double xi_179 = _data_edgeFaceRhs[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_190 = _data_edgeFaceRhs[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_171 = -xi_68*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_172 = -xi_69*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = -xi_70*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17];
         const double xi_174 = -xi_71*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = -xi_72*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_176 = -xi_73*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18];
         const double xi_177 = -xi_74*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = -xi_75*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 19];
         const double xi_182 = -xi_156*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_183 = -xi_157*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16];
         const double xi_184 = -xi_158*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_185 = -xi_159*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_186 = -xi_160*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17];
         const double xi_187 = -xi_161*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_188 = -xi_162*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18];
         const double xi_189 = -xi_163*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178 + xi_179) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 15; ctr_1 < 16; ctr_1 += 1)
      {
         const double xi_115 = _data_edgeFaceRhs[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = -xi_156*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_109 = -xi_157*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16];
         const double xi_110 = -xi_158*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = -xi_159*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_112 = -xi_160*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17];
         const double xi_113 = -xi_161*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_114 = -xi_162*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18];
         const double xi_107 = -xi_163*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 15; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_271 = _data_edgeFaceRhs[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_282 = _data_edgeFaceRhs[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_263 = -xi_239*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_264 = -xi_240*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_265 = -xi_241*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_266 = -xi_242*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16];
         const double xi_267 = -xi_243*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_268 = -xi_244*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18];
         const double xi_269 = -xi_245*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_270 = -xi_246*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_274 = -xi_68*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_275 = -xi_69*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_276 = -xi_70*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17];
         const double xi_277 = -xi_71*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_278 = -xi_72*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_279 = -xi_73*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18];
         const double xi_280 = -xi_74*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_281 = -xi_75*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 19];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_263 + xi_264 + xi_265 + xi_266 + xi_267 + xi_268 + xi_269 + xi_270 + xi_271) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 15; ctr_1 += 1)
      {
         const double xi_43 = _data_edgeFaceRhs[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = _data_edgeFaceRhs[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = _data_edgeFaceRhs[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = -xi_239*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_36 = -xi_240*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_37 = -xi_241*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_242*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16];
         const double xi_39 = -xi_243*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_40 = -xi_244*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18];
         const double xi_41 = -xi_245*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_42 = -xi_246*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = -xi_68*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = -xi_69*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_48 = -xi_70*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17];
         const double xi_49 = -xi_71*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_50 = -xi_72*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_51 = -xi_73*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18];
         const double xi_52 = -xi_74*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_53 = -xi_75*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 19];
         const double xi_57 = -xi_156*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_58 = -xi_157*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16];
         const double xi_59 = -xi_158*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = -xi_159*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_61 = -xi_160*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17];
         const double xi_62 = -xi_161*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = -xi_162*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18];
         const double xi_64 = -xi_163*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 15; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
      {
         const double xi_225 = _data_edgeFaceRhs[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_236 = _data_edgeFaceRhs[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_217 = -xi_239*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_218 = -xi_240*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_219 = -xi_241*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_220 = -xi_242*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16];
         const double xi_221 = -xi_243*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_222 = -xi_244*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18];
         const double xi_223 = -xi_245*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_224 = -xi_246*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_228 = -xi_156*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_229 = -xi_157*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16];
         const double xi_230 = -xi_158*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_231 = -xi_159*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_232 = -xi_160*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17];
         const double xi_233 = -xi_161*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_234 = -xi_162*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18];
         const double xi_235 = -xi_163*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_217 + xi_218 + xi_219 + xi_220 + xi_221 + xi_222 + xi_223 + xi_224 + xi_225) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 15; ctr_2 < 16; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_140 = _data_edgeFaceRhs[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_132 = -xi_239*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_133 = -xi_240*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_134 = -xi_241*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = -xi_242*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16];
         const double xi_136 = -xi_243*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_137 = -xi_244*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18];
         const double xi_138 = -xi_245*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17];
         const double xi_139 = -xi_246*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_132 + xi_133 + xi_134 + xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_5(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceRhs, double const * const _data_edge_stencil_at_edge_x, double const * const _data_edge_stencil_at_edge_xy, double const * const _data_edge_stencil_at_edge_y, double * RESTRICT _data_vertexFaceDst, double const * const _data_vertex_stencil_at_edge_x, double const * const _data_vertex_stencil_at_edge_xy, double const * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_80 = 1.0;
   const double xi_81 = -relax;
   const double xi_67 = _data_edge_stencil_at_edge_xy[0];
   const double xi_78 = 1 / (xi_67);
   const double xi_68 = _data_edge_stencil_at_edge_xy[1];
   const double xi_69 = _data_edge_stencil_at_edge_xy[4];
   const double xi_70 = _data_edge_stencil_at_edge_xy[3];
   const double xi_71 = _data_edge_stencil_at_edge_xy[2];
   const double xi_72 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_73 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_74 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_155 = _data_edge_stencil_at_edge_y[0];
   const double xi_167 = 1 / (xi_155);
   const double xi_156 = _data_edge_stencil_at_edge_y[4];
   const double xi_157 = _data_edge_stencil_at_edge_y[3];
   const double xi_158 = _data_edge_stencil_at_edge_y[1];
   const double xi_159 = _data_edge_stencil_at_edge_y[2];
   const double xi_160 = _data_vertex_stencil_at_edge_y[3];
   const double xi_161 = _data_vertex_stencil_at_edge_y[0];
   const double xi_162 = _data_vertex_stencil_at_edge_y[2];
   const double xi_163 = _data_vertex_stencil_at_edge_y[1];
   const double xi_238 = _data_edge_stencil_at_edge_x[0];
   const double xi_258 = 1 / (xi_238);
   const double xi_239 = _data_edge_stencil_at_edge_x[1];
   const double xi_240 = _data_edge_stencil_at_edge_x[3];
   const double xi_241 = _data_edge_stencil_at_edge_x[4];
   const double xi_242 = _data_edge_stencil_at_edge_x[2];
   const double xi_243 = _data_vertex_stencil_at_edge_x[1];
   const double xi_244 = _data_vertex_stencil_at_edge_x[3];
   const double xi_245 = _data_vertex_stencil_at_edge_x[2];
   const double xi_246 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_90 = _data_edgeFaceRhs[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_82 = -xi_68*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = -xi_69*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_84 = -xi_70*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33];
         const double xi_85 = -xi_71*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_86 = -xi_72*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = -xi_73*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34];
         const double xi_88 = -xi_74*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_89 = -xi_75*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 35];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
      {
         const double xi_179 = _data_edgeFaceRhs[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_190 = _data_edgeFaceRhs[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_171 = -xi_68*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_172 = -xi_69*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = -xi_70*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33];
         const double xi_174 = -xi_71*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = -xi_72*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_176 = -xi_73*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34];
         const double xi_177 = -xi_74*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = -xi_75*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 35];
         const double xi_182 = -xi_156*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_183 = -xi_157*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32];
         const double xi_184 = -xi_158*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_185 = -xi_159*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_186 = -xi_160*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33];
         const double xi_187 = -xi_161*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_188 = -xi_162*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34];
         const double xi_189 = -xi_163*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178 + xi_179) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 31; ctr_1 < 32; ctr_1 += 1)
      {
         const double xi_115 = _data_edgeFaceRhs[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = -xi_156*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_109 = -xi_157*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32];
         const double xi_110 = -xi_158*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = -xi_159*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_112 = -xi_160*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33];
         const double xi_113 = -xi_161*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_114 = -xi_162*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34];
         const double xi_107 = -xi_163*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 31; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_271 = _data_edgeFaceRhs[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_282 = _data_edgeFaceRhs[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_263 = -xi_239*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_264 = -xi_240*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_265 = -xi_241*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_266 = -xi_242*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32];
         const double xi_267 = -xi_243*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_268 = -xi_244*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34];
         const double xi_269 = -xi_245*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_270 = -xi_246*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_274 = -xi_68*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_275 = -xi_69*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_276 = -xi_70*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33];
         const double xi_277 = -xi_71*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_278 = -xi_72*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_279 = -xi_73*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34];
         const double xi_280 = -xi_74*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_281 = -xi_75*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 35];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_263 + xi_264 + xi_265 + xi_266 + xi_267 + xi_268 + xi_269 + xi_270 + xi_271) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 31; ctr_1 += 1)
      {
         const double xi_43 = _data_edgeFaceRhs[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = _data_edgeFaceRhs[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = _data_edgeFaceRhs[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = -xi_239*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_36 = -xi_240*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_37 = -xi_241*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_242*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32];
         const double xi_39 = -xi_243*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_40 = -xi_244*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34];
         const double xi_41 = -xi_245*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_42 = -xi_246*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = -xi_68*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = -xi_69*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_48 = -xi_70*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33];
         const double xi_49 = -xi_71*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_50 = -xi_72*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_51 = -xi_73*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34];
         const double xi_52 = -xi_74*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_53 = -xi_75*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 35];
         const double xi_57 = -xi_156*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_58 = -xi_157*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32];
         const double xi_59 = -xi_158*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = -xi_159*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_61 = -xi_160*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33];
         const double xi_62 = -xi_161*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = -xi_162*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34];
         const double xi_64 = -xi_163*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 31; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
      {
         const double xi_225 = _data_edgeFaceRhs[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_236 = _data_edgeFaceRhs[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_217 = -xi_239*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_218 = -xi_240*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_219 = -xi_241*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_220 = -xi_242*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32];
         const double xi_221 = -xi_243*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_222 = -xi_244*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34];
         const double xi_223 = -xi_245*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_224 = -xi_246*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_228 = -xi_156*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_229 = -xi_157*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32];
         const double xi_230 = -xi_158*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_231 = -xi_159*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_232 = -xi_160*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33];
         const double xi_233 = -xi_161*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_234 = -xi_162*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34];
         const double xi_235 = -xi_163*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_217 + xi_218 + xi_219 + xi_220 + xi_221 + xi_222 + xi_223 + xi_224 + xi_225) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 31; ctr_2 < 32; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_140 = _data_edgeFaceRhs[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_132 = -xi_239*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_133 = -xi_240*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_134 = -xi_241*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = -xi_242*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32];
         const double xi_136 = -xi_243*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_137 = -xi_244*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34];
         const double xi_138 = -xi_245*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33];
         const double xi_139 = -xi_246*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_132 + xi_133 + xi_134 + xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_6(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceRhs, double const * const _data_edge_stencil_at_edge_x, double const * const _data_edge_stencil_at_edge_xy, double const * const _data_edge_stencil_at_edge_y, double * RESTRICT _data_vertexFaceDst, double const * const _data_vertex_stencil_at_edge_x, double const * const _data_vertex_stencil_at_edge_xy, double const * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_80 = 1.0;
   const double xi_81 = -relax;
   const double xi_67 = _data_edge_stencil_at_edge_xy[0];
   const double xi_78 = 1 / (xi_67);
   const double xi_68 = _data_edge_stencil_at_edge_xy[1];
   const double xi_69 = _data_edge_stencil_at_edge_xy[4];
   const double xi_70 = _data_edge_stencil_at_edge_xy[3];
   const double xi_71 = _data_edge_stencil_at_edge_xy[2];
   const double xi_72 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_73 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_74 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_155 = _data_edge_stencil_at_edge_y[0];
   const double xi_167 = 1 / (xi_155);
   const double xi_156 = _data_edge_stencil_at_edge_y[4];
   const double xi_157 = _data_edge_stencil_at_edge_y[3];
   const double xi_158 = _data_edge_stencil_at_edge_y[1];
   const double xi_159 = _data_edge_stencil_at_edge_y[2];
   const double xi_160 = _data_vertex_stencil_at_edge_y[3];
   const double xi_161 = _data_vertex_stencil_at_edge_y[0];
   const double xi_162 = _data_vertex_stencil_at_edge_y[2];
   const double xi_163 = _data_vertex_stencil_at_edge_y[1];
   const double xi_238 = _data_edge_stencil_at_edge_x[0];
   const double xi_258 = 1 / (xi_238);
   const double xi_239 = _data_edge_stencil_at_edge_x[1];
   const double xi_240 = _data_edge_stencil_at_edge_x[3];
   const double xi_241 = _data_edge_stencil_at_edge_x[4];
   const double xi_242 = _data_edge_stencil_at_edge_x[2];
   const double xi_243 = _data_vertex_stencil_at_edge_x[1];
   const double xi_244 = _data_vertex_stencil_at_edge_x[3];
   const double xi_245 = _data_vertex_stencil_at_edge_x[2];
   const double xi_246 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_90 = _data_edgeFaceRhs[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_82 = -xi_68*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = -xi_69*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_84 = -xi_70*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65];
         const double xi_85 = -xi_71*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_86 = -xi_72*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = -xi_73*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66];
         const double xi_88 = -xi_74*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_89 = -xi_75*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 67];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
      {
         const double xi_179 = _data_edgeFaceRhs[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_190 = _data_edgeFaceRhs[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_171 = -xi_68*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_172 = -xi_69*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = -xi_70*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65];
         const double xi_174 = -xi_71*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = -xi_72*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_176 = -xi_73*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66];
         const double xi_177 = -xi_74*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = -xi_75*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 67];
         const double xi_182 = -xi_156*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_183 = -xi_157*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64];
         const double xi_184 = -xi_158*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_185 = -xi_159*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_186 = -xi_160*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65];
         const double xi_187 = -xi_161*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_188 = -xi_162*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66];
         const double xi_189 = -xi_163*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178 + xi_179) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 63; ctr_1 < 64; ctr_1 += 1)
      {
         const double xi_115 = _data_edgeFaceRhs[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = -xi_156*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_109 = -xi_157*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64];
         const double xi_110 = -xi_158*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = -xi_159*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_112 = -xi_160*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65];
         const double xi_113 = -xi_161*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_114 = -xi_162*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66];
         const double xi_107 = -xi_163*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 63; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_271 = _data_edgeFaceRhs[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_282 = _data_edgeFaceRhs[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_263 = -xi_239*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_264 = -xi_240*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_265 = -xi_241*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_266 = -xi_242*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64];
         const double xi_267 = -xi_243*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_268 = -xi_244*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66];
         const double xi_269 = -xi_245*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_270 = -xi_246*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_274 = -xi_68*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_275 = -xi_69*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_276 = -xi_70*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65];
         const double xi_277 = -xi_71*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_278 = -xi_72*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_279 = -xi_73*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66];
         const double xi_280 = -xi_74*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_281 = -xi_75*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 67];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_263 + xi_264 + xi_265 + xi_266 + xi_267 + xi_268 + xi_269 + xi_270 + xi_271) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 63; ctr_1 += 1)
      {
         const double xi_43 = _data_edgeFaceRhs[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = _data_edgeFaceRhs[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = _data_edgeFaceRhs[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = -xi_239*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_36 = -xi_240*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_37 = -xi_241*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_242*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64];
         const double xi_39 = -xi_243*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_40 = -xi_244*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66];
         const double xi_41 = -xi_245*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_42 = -xi_246*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = -xi_68*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = -xi_69*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_48 = -xi_70*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65];
         const double xi_49 = -xi_71*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_50 = -xi_72*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_51 = -xi_73*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66];
         const double xi_52 = -xi_74*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_53 = -xi_75*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 67];
         const double xi_57 = -xi_156*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_58 = -xi_157*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64];
         const double xi_59 = -xi_158*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = -xi_159*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_61 = -xi_160*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65];
         const double xi_62 = -xi_161*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = -xi_162*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66];
         const double xi_64 = -xi_163*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 63; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
      {
         const double xi_225 = _data_edgeFaceRhs[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_236 = _data_edgeFaceRhs[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_217 = -xi_239*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_218 = -xi_240*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_219 = -xi_241*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_220 = -xi_242*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64];
         const double xi_221 = -xi_243*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_222 = -xi_244*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66];
         const double xi_223 = -xi_245*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_224 = -xi_246*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_228 = -xi_156*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_229 = -xi_157*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64];
         const double xi_230 = -xi_158*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_231 = -xi_159*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_232 = -xi_160*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65];
         const double xi_233 = -xi_161*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_234 = -xi_162*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66];
         const double xi_235 = -xi_163*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_217 + xi_218 + xi_219 + xi_220 + xi_221 + xi_222 + xi_223 + xi_224 + xi_225) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 63; ctr_2 < 64; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_140 = _data_edgeFaceRhs[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_132 = -xi_239*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_133 = -xi_240*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_134 = -xi_241*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = -xi_242*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64];
         const double xi_136 = -xi_243*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_137 = -xi_244*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66];
         const double xi_138 = -xi_245*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65];
         const double xi_139 = -xi_246*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_132 + xi_133 + xi_134 + xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_7(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceRhs, double const * const _data_edge_stencil_at_edge_x, double const * const _data_edge_stencil_at_edge_xy, double const * const _data_edge_stencil_at_edge_y, double * RESTRICT _data_vertexFaceDst, double const * const _data_vertex_stencil_at_edge_x, double const * const _data_vertex_stencil_at_edge_xy, double const * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_80 = 1.0;
   const double xi_81 = -relax;
   const double xi_67 = _data_edge_stencil_at_edge_xy[0];
   const double xi_78 = 1 / (xi_67);
   const double xi_68 = _data_edge_stencil_at_edge_xy[1];
   const double xi_69 = _data_edge_stencil_at_edge_xy[4];
   const double xi_70 = _data_edge_stencil_at_edge_xy[3];
   const double xi_71 = _data_edge_stencil_at_edge_xy[2];
   const double xi_72 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_73 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_74 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_155 = _data_edge_stencil_at_edge_y[0];
   const double xi_167 = 1 / (xi_155);
   const double xi_156 = _data_edge_stencil_at_edge_y[4];
   const double xi_157 = _data_edge_stencil_at_edge_y[3];
   const double xi_158 = _data_edge_stencil_at_edge_y[1];
   const double xi_159 = _data_edge_stencil_at_edge_y[2];
   const double xi_160 = _data_vertex_stencil_at_edge_y[3];
   const double xi_161 = _data_vertex_stencil_at_edge_y[0];
   const double xi_162 = _data_vertex_stencil_at_edge_y[2];
   const double xi_163 = _data_vertex_stencil_at_edge_y[1];
   const double xi_238 = _data_edge_stencil_at_edge_x[0];
   const double xi_258 = 1 / (xi_238);
   const double xi_239 = _data_edge_stencil_at_edge_x[1];
   const double xi_240 = _data_edge_stencil_at_edge_x[3];
   const double xi_241 = _data_edge_stencil_at_edge_x[4];
   const double xi_242 = _data_edge_stencil_at_edge_x[2];
   const double xi_243 = _data_vertex_stencil_at_edge_x[1];
   const double xi_244 = _data_vertex_stencil_at_edge_x[3];
   const double xi_245 = _data_vertex_stencil_at_edge_x[2];
   const double xi_246 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_90 = _data_edgeFaceRhs[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_82 = -xi_68*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = -xi_69*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_84 = -xi_70*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129];
         const double xi_85 = -xi_71*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_86 = -xi_72*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = -xi_73*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130];
         const double xi_88 = -xi_74*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_89 = -xi_75*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 131];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 127; ctr_1 += 1)
      {
         const double xi_179 = _data_edgeFaceRhs[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_190 = _data_edgeFaceRhs[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_171 = -xi_68*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_172 = -xi_69*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = -xi_70*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129];
         const double xi_174 = -xi_71*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = -xi_72*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_176 = -xi_73*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130];
         const double xi_177 = -xi_74*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = -xi_75*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 131];
         const double xi_182 = -xi_156*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_183 = -xi_157*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128];
         const double xi_184 = -xi_158*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_185 = -xi_159*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_186 = -xi_160*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129];
         const double xi_187 = -xi_161*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_188 = -xi_162*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130];
         const double xi_189 = -xi_163*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178 + xi_179) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 127; ctr_1 < 128; ctr_1 += 1)
      {
         const double xi_115 = _data_edgeFaceRhs[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = -xi_156*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_109 = -xi_157*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128];
         const double xi_110 = -xi_158*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = -xi_159*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_112 = -xi_160*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129];
         const double xi_113 = -xi_161*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_114 = -xi_162*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130];
         const double xi_107 = -xi_163*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 127; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_271 = _data_edgeFaceRhs[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_282 = _data_edgeFaceRhs[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_263 = -xi_239*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_264 = -xi_240*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_265 = -xi_241*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_266 = -xi_242*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128];
         const double xi_267 = -xi_243*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_268 = -xi_244*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130];
         const double xi_269 = -xi_245*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_270 = -xi_246*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_274 = -xi_68*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_275 = -xi_69*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_276 = -xi_70*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129];
         const double xi_277 = -xi_71*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_278 = -xi_72*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_279 = -xi_73*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130];
         const double xi_280 = -xi_74*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_281 = -xi_75*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 131];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_263 + xi_264 + xi_265 + xi_266 + xi_267 + xi_268 + xi_269 + xi_270 + xi_271) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 127; ctr_1 += 1)
      {
         const double xi_43 = _data_edgeFaceRhs[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = _data_edgeFaceRhs[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = _data_edgeFaceRhs[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = -xi_239*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_36 = -xi_240*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_37 = -xi_241*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_242*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128];
         const double xi_39 = -xi_243*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_40 = -xi_244*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130];
         const double xi_41 = -xi_245*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_42 = -xi_246*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = -xi_68*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = -xi_69*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_48 = -xi_70*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129];
         const double xi_49 = -xi_71*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_50 = -xi_72*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_51 = -xi_73*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130];
         const double xi_52 = -xi_74*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_53 = -xi_75*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 131];
         const double xi_57 = -xi_156*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_58 = -xi_157*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128];
         const double xi_59 = -xi_158*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = -xi_159*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_61 = -xi_160*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129];
         const double xi_62 = -xi_161*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = -xi_162*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130];
         const double xi_64 = -xi_163*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 127; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
      {
         const double xi_225 = _data_edgeFaceRhs[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_236 = _data_edgeFaceRhs[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_217 = -xi_239*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_218 = -xi_240*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_219 = -xi_241*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_220 = -xi_242*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128];
         const double xi_221 = -xi_243*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_222 = -xi_244*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130];
         const double xi_223 = -xi_245*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_224 = -xi_246*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_228 = -xi_156*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_229 = -xi_157*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128];
         const double xi_230 = -xi_158*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_231 = -xi_159*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_232 = -xi_160*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129];
         const double xi_233 = -xi_161*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_234 = -xi_162*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130];
         const double xi_235 = -xi_163*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_217 + xi_218 + xi_219 + xi_220 + xi_221 + xi_222 + xi_223 + xi_224 + xi_225) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 127; ctr_2 < 128; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_140 = _data_edgeFaceRhs[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_132 = -xi_239*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_133 = -xi_240*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_134 = -xi_241*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = -xi_242*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128];
         const double xi_136 = -xi_243*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_137 = -xi_244*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130];
         const double xi_138 = -xi_245*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129];
         const double xi_139 = -xi_246*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_132 + xi_133 + xi_134 + xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_8(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceRhs, double const * const _data_edge_stencil_at_edge_x, double const * const _data_edge_stencil_at_edge_xy, double const * const _data_edge_stencil_at_edge_y, double * RESTRICT _data_vertexFaceDst, double const * const _data_vertex_stencil_at_edge_x, double const * const _data_vertex_stencil_at_edge_xy, double const * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_80 = 1.0;
   const double xi_81 = -relax;
   const double xi_67 = _data_edge_stencil_at_edge_xy[0];
   const double xi_78 = 1 / (xi_67);
   const double xi_68 = _data_edge_stencil_at_edge_xy[1];
   const double xi_69 = _data_edge_stencil_at_edge_xy[4];
   const double xi_70 = _data_edge_stencil_at_edge_xy[3];
   const double xi_71 = _data_edge_stencil_at_edge_xy[2];
   const double xi_72 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_73 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_74 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_155 = _data_edge_stencil_at_edge_y[0];
   const double xi_167 = 1 / (xi_155);
   const double xi_156 = _data_edge_stencil_at_edge_y[4];
   const double xi_157 = _data_edge_stencil_at_edge_y[3];
   const double xi_158 = _data_edge_stencil_at_edge_y[1];
   const double xi_159 = _data_edge_stencil_at_edge_y[2];
   const double xi_160 = _data_vertex_stencil_at_edge_y[3];
   const double xi_161 = _data_vertex_stencil_at_edge_y[0];
   const double xi_162 = _data_vertex_stencil_at_edge_y[2];
   const double xi_163 = _data_vertex_stencil_at_edge_y[1];
   const double xi_238 = _data_edge_stencil_at_edge_x[0];
   const double xi_258 = 1 / (xi_238);
   const double xi_239 = _data_edge_stencil_at_edge_x[1];
   const double xi_240 = _data_edge_stencil_at_edge_x[3];
   const double xi_241 = _data_edge_stencil_at_edge_x[4];
   const double xi_242 = _data_edge_stencil_at_edge_x[2];
   const double xi_243 = _data_vertex_stencil_at_edge_x[1];
   const double xi_244 = _data_vertex_stencil_at_edge_x[3];
   const double xi_245 = _data_vertex_stencil_at_edge_x[2];
   const double xi_246 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_90 = _data_edgeFaceRhs[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_82 = -xi_68*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = -xi_69*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_84 = -xi_70*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257];
         const double xi_85 = -xi_71*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_86 = -xi_72*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = -xi_73*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258];
         const double xi_88 = -xi_74*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_89 = -xi_75*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 259];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
      {
         const double xi_179 = _data_edgeFaceRhs[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_190 = _data_edgeFaceRhs[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_171 = -xi_68*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_172 = -xi_69*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = -xi_70*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257];
         const double xi_174 = -xi_71*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = -xi_72*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_176 = -xi_73*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258];
         const double xi_177 = -xi_74*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = -xi_75*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 259];
         const double xi_182 = -xi_156*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_183 = -xi_157*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256];
         const double xi_184 = -xi_158*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_185 = -xi_159*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_186 = -xi_160*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257];
         const double xi_187 = -xi_161*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_188 = -xi_162*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258];
         const double xi_189 = -xi_163*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178 + xi_179) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 255; ctr_1 < 256; ctr_1 += 1)
      {
         const double xi_115 = _data_edgeFaceRhs[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = -xi_156*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_109 = -xi_157*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256];
         const double xi_110 = -xi_158*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = -xi_159*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_112 = -xi_160*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257];
         const double xi_113 = -xi_161*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_114 = -xi_162*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258];
         const double xi_107 = -xi_163*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 255; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_271 = _data_edgeFaceRhs[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_282 = _data_edgeFaceRhs[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_263 = -xi_239*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_264 = -xi_240*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_265 = -xi_241*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_266 = -xi_242*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256];
         const double xi_267 = -xi_243*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_268 = -xi_244*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258];
         const double xi_269 = -xi_245*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_270 = -xi_246*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_274 = -xi_68*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_275 = -xi_69*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_276 = -xi_70*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257];
         const double xi_277 = -xi_71*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_278 = -xi_72*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_279 = -xi_73*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258];
         const double xi_280 = -xi_74*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_281 = -xi_75*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 259];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_263 + xi_264 + xi_265 + xi_266 + xi_267 + xi_268 + xi_269 + xi_270 + xi_271) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 255; ctr_1 += 1)
      {
         const double xi_43 = _data_edgeFaceRhs[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = _data_edgeFaceRhs[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = _data_edgeFaceRhs[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = -xi_239*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_36 = -xi_240*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_37 = -xi_241*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_242*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256];
         const double xi_39 = -xi_243*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_40 = -xi_244*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258];
         const double xi_41 = -xi_245*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_42 = -xi_246*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = -xi_68*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = -xi_69*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_48 = -xi_70*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257];
         const double xi_49 = -xi_71*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_50 = -xi_72*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_51 = -xi_73*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258];
         const double xi_52 = -xi_74*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_53 = -xi_75*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 259];
         const double xi_57 = -xi_156*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_58 = -xi_157*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256];
         const double xi_59 = -xi_158*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = -xi_159*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_61 = -xi_160*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257];
         const double xi_62 = -xi_161*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = -xi_162*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258];
         const double xi_64 = -xi_163*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 255; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
      {
         const double xi_225 = _data_edgeFaceRhs[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_236 = _data_edgeFaceRhs[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_217 = -xi_239*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_218 = -xi_240*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_219 = -xi_241*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_220 = -xi_242*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256];
         const double xi_221 = -xi_243*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_222 = -xi_244*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258];
         const double xi_223 = -xi_245*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_224 = -xi_246*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_228 = -xi_156*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_229 = -xi_157*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256];
         const double xi_230 = -xi_158*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_231 = -xi_159*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_232 = -xi_160*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257];
         const double xi_233 = -xi_161*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_234 = -xi_162*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258];
         const double xi_235 = -xi_163*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_217 + xi_218 + xi_219 + xi_220 + xi_221 + xi_222 + xi_223 + xi_224 + xi_225) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 255; ctr_2 < 256; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_140 = _data_edgeFaceRhs[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_132 = -xi_239*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_133 = -xi_240*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_134 = -xi_241*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = -xi_242*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256];
         const double xi_136 = -xi_243*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_137 = -xi_244*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258];
         const double xi_138 = -xi_245*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257];
         const double xi_139 = -xi_246*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_132 + xi_133 + xi_134 + xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_9(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceRhs, double const * const _data_edge_stencil_at_edge_x, double const * const _data_edge_stencil_at_edge_xy, double const * const _data_edge_stencil_at_edge_y, double * RESTRICT _data_vertexFaceDst, double const * const _data_vertex_stencil_at_edge_x, double const * const _data_vertex_stencil_at_edge_xy, double const * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_80 = 1.0;
   const double xi_81 = -relax;
   const double xi_67 = _data_edge_stencil_at_edge_xy[0];
   const double xi_78 = 1 / (xi_67);
   const double xi_68 = _data_edge_stencil_at_edge_xy[1];
   const double xi_69 = _data_edge_stencil_at_edge_xy[4];
   const double xi_70 = _data_edge_stencil_at_edge_xy[3];
   const double xi_71 = _data_edge_stencil_at_edge_xy[2];
   const double xi_72 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_73 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_74 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_155 = _data_edge_stencil_at_edge_y[0];
   const double xi_167 = 1 / (xi_155);
   const double xi_156 = _data_edge_stencil_at_edge_y[4];
   const double xi_157 = _data_edge_stencil_at_edge_y[3];
   const double xi_158 = _data_edge_stencil_at_edge_y[1];
   const double xi_159 = _data_edge_stencil_at_edge_y[2];
   const double xi_160 = _data_vertex_stencil_at_edge_y[3];
   const double xi_161 = _data_vertex_stencil_at_edge_y[0];
   const double xi_162 = _data_vertex_stencil_at_edge_y[2];
   const double xi_163 = _data_vertex_stencil_at_edge_y[1];
   const double xi_238 = _data_edge_stencil_at_edge_x[0];
   const double xi_258 = 1 / (xi_238);
   const double xi_239 = _data_edge_stencil_at_edge_x[1];
   const double xi_240 = _data_edge_stencil_at_edge_x[3];
   const double xi_241 = _data_edge_stencil_at_edge_x[4];
   const double xi_242 = _data_edge_stencil_at_edge_x[2];
   const double xi_243 = _data_vertex_stencil_at_edge_x[1];
   const double xi_244 = _data_vertex_stencil_at_edge_x[3];
   const double xi_245 = _data_vertex_stencil_at_edge_x[2];
   const double xi_246 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_90 = _data_edgeFaceRhs[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_82 = -xi_68*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = -xi_69*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_84 = -xi_70*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513];
         const double xi_85 = -xi_71*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_86 = -xi_72*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = -xi_73*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514];
         const double xi_88 = -xi_74*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_89 = -xi_75*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 515];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
      {
         const double xi_179 = _data_edgeFaceRhs[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_190 = _data_edgeFaceRhs[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_171 = -xi_68*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_172 = -xi_69*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = -xi_70*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513];
         const double xi_174 = -xi_71*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = -xi_72*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_176 = -xi_73*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514];
         const double xi_177 = -xi_74*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = -xi_75*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 515];
         const double xi_182 = -xi_156*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_183 = -xi_157*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512];
         const double xi_184 = -xi_158*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_185 = -xi_159*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_186 = -xi_160*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513];
         const double xi_187 = -xi_161*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_188 = -xi_162*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514];
         const double xi_189 = -xi_163*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178 + xi_179) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 511; ctr_1 < 512; ctr_1 += 1)
      {
         const double xi_115 = _data_edgeFaceRhs[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = -xi_156*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_109 = -xi_157*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512];
         const double xi_110 = -xi_158*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = -xi_159*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_112 = -xi_160*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513];
         const double xi_113 = -xi_161*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_114 = -xi_162*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514];
         const double xi_107 = -xi_163*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 511; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_271 = _data_edgeFaceRhs[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_282 = _data_edgeFaceRhs[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_263 = -xi_239*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_264 = -xi_240*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_265 = -xi_241*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_266 = -xi_242*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512];
         const double xi_267 = -xi_243*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_268 = -xi_244*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514];
         const double xi_269 = -xi_245*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_270 = -xi_246*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_274 = -xi_68*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_275 = -xi_69*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_276 = -xi_70*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513];
         const double xi_277 = -xi_71*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_278 = -xi_72*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_279 = -xi_73*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514];
         const double xi_280 = -xi_74*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_281 = -xi_75*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 515];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_263 + xi_264 + xi_265 + xi_266 + xi_267 + xi_268 + xi_269 + xi_270 + xi_271) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 511; ctr_1 += 1)
      {
         const double xi_43 = _data_edgeFaceRhs[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = _data_edgeFaceRhs[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = _data_edgeFaceRhs[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = -xi_239*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_36 = -xi_240*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_37 = -xi_241*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_242*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512];
         const double xi_39 = -xi_243*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_40 = -xi_244*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514];
         const double xi_41 = -xi_245*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_42 = -xi_246*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = -xi_68*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = -xi_69*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_48 = -xi_70*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513];
         const double xi_49 = -xi_71*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_50 = -xi_72*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_51 = -xi_73*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514];
         const double xi_52 = -xi_74*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_53 = -xi_75*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 515];
         const double xi_57 = -xi_156*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_58 = -xi_157*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512];
         const double xi_59 = -xi_158*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = -xi_159*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_61 = -xi_160*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513];
         const double xi_62 = -xi_161*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = -xi_162*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514];
         const double xi_64 = -xi_163*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 511; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
      {
         const double xi_225 = _data_edgeFaceRhs[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_236 = _data_edgeFaceRhs[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_217 = -xi_239*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_218 = -xi_240*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_219 = -xi_241*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_220 = -xi_242*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512];
         const double xi_221 = -xi_243*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_222 = -xi_244*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514];
         const double xi_223 = -xi_245*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_224 = -xi_246*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_228 = -xi_156*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_229 = -xi_157*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512];
         const double xi_230 = -xi_158*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_231 = -xi_159*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_232 = -xi_160*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513];
         const double xi_233 = -xi_161*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_234 = -xi_162*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514];
         const double xi_235 = -xi_163*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_217 + xi_218 + xi_219 + xi_220 + xi_221 + xi_222 + xi_223 + xi_224 + xi_225) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 511; ctr_2 < 512; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_140 = _data_edgeFaceRhs[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_132 = -xi_239*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_133 = -xi_240*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_134 = -xi_241*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = -xi_242*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512];
         const double xi_136 = -xi_243*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_137 = -xi_244*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514];
         const double xi_138 = -xi_245*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513];
         const double xi_139 = -xi_246*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_132 + xi_133 + xi_134 + xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_10(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceRhs, double const * const _data_edge_stencil_at_edge_x, double const * const _data_edge_stencil_at_edge_xy, double const * const _data_edge_stencil_at_edge_y, double * RESTRICT _data_vertexFaceDst, double const * const _data_vertex_stencil_at_edge_x, double const * const _data_vertex_stencil_at_edge_xy, double const * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_80 = 1.0;
   const double xi_81 = -relax;
   const double xi_67 = _data_edge_stencil_at_edge_xy[0];
   const double xi_78 = 1 / (xi_67);
   const double xi_68 = _data_edge_stencil_at_edge_xy[1];
   const double xi_69 = _data_edge_stencil_at_edge_xy[4];
   const double xi_70 = _data_edge_stencil_at_edge_xy[3];
   const double xi_71 = _data_edge_stencil_at_edge_xy[2];
   const double xi_72 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_73 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_74 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_155 = _data_edge_stencil_at_edge_y[0];
   const double xi_167 = 1 / (xi_155);
   const double xi_156 = _data_edge_stencil_at_edge_y[4];
   const double xi_157 = _data_edge_stencil_at_edge_y[3];
   const double xi_158 = _data_edge_stencil_at_edge_y[1];
   const double xi_159 = _data_edge_stencil_at_edge_y[2];
   const double xi_160 = _data_vertex_stencil_at_edge_y[3];
   const double xi_161 = _data_vertex_stencil_at_edge_y[0];
   const double xi_162 = _data_vertex_stencil_at_edge_y[2];
   const double xi_163 = _data_vertex_stencil_at_edge_y[1];
   const double xi_238 = _data_edge_stencil_at_edge_x[0];
   const double xi_258 = 1 / (xi_238);
   const double xi_239 = _data_edge_stencil_at_edge_x[1];
   const double xi_240 = _data_edge_stencil_at_edge_x[3];
   const double xi_241 = _data_edge_stencil_at_edge_x[4];
   const double xi_242 = _data_edge_stencil_at_edge_x[2];
   const double xi_243 = _data_vertex_stencil_at_edge_x[1];
   const double xi_244 = _data_vertex_stencil_at_edge_x[3];
   const double xi_245 = _data_vertex_stencil_at_edge_x[2];
   const double xi_246 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_90 = _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_82 = -xi_68*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = -xi_69*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_84 = -xi_70*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025];
         const double xi_85 = -xi_71*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_86 = -xi_72*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = -xi_73*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026];
         const double xi_88 = -xi_74*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_89 = -xi_75*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1027];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 1023; ctr_1 += 1)
      {
         const double xi_179 = _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_190 = _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_171 = -xi_68*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_172 = -xi_69*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = -xi_70*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025];
         const double xi_174 = -xi_71*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = -xi_72*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_176 = -xi_73*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026];
         const double xi_177 = -xi_74*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = -xi_75*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1027];
         const double xi_182 = -xi_156*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_183 = -xi_157*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024];
         const double xi_184 = -xi_158*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_185 = -xi_159*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_186 = -xi_160*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025];
         const double xi_187 = -xi_161*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_188 = -xi_162*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026];
         const double xi_189 = -xi_163*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178 + xi_179) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 1023; ctr_1 < 1024; ctr_1 += 1)
      {
         const double xi_115 = _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = -xi_156*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_109 = -xi_157*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024];
         const double xi_110 = -xi_158*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = -xi_159*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_112 = -xi_160*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025];
         const double xi_113 = -xi_161*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_114 = -xi_162*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026];
         const double xi_107 = -xi_163*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 1023; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_271 = _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_282 = _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_263 = -xi_239*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_264 = -xi_240*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_265 = -xi_241*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_266 = -xi_242*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024];
         const double xi_267 = -xi_243*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_268 = -xi_244*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026];
         const double xi_269 = -xi_245*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_270 = -xi_246*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_274 = -xi_68*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_275 = -xi_69*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_276 = -xi_70*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025];
         const double xi_277 = -xi_71*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_278 = -xi_72*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_279 = -xi_73*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026];
         const double xi_280 = -xi_74*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_281 = -xi_75*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1027];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_263 + xi_264 + xi_265 + xi_266 + xi_267 + xi_268 + xi_269 + xi_270 + xi_271) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1023; ctr_1 += 1)
      {
         const double xi_43 = _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = -xi_239*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_36 = -xi_240*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_37 = -xi_241*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_242*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024];
         const double xi_39 = -xi_243*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_40 = -xi_244*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026];
         const double xi_41 = -xi_245*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_42 = -xi_246*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = -xi_68*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = -xi_69*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_48 = -xi_70*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025];
         const double xi_49 = -xi_71*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_50 = -xi_72*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_51 = -xi_73*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026];
         const double xi_52 = -xi_74*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_53 = -xi_75*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1027];
         const double xi_57 = -xi_156*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_58 = -xi_157*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024];
         const double xi_59 = -xi_158*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = -xi_159*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_61 = -xi_160*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025];
         const double xi_62 = -xi_161*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = -xi_162*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026];
         const double xi_64 = -xi_163*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 1023; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
      {
         const double xi_225 = _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_236 = _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_217 = -xi_239*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_218 = -xi_240*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_219 = -xi_241*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_220 = -xi_242*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024];
         const double xi_221 = -xi_243*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_222 = -xi_244*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026];
         const double xi_223 = -xi_245*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_224 = -xi_246*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_228 = -xi_156*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_229 = -xi_157*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024];
         const double xi_230 = -xi_158*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_231 = -xi_159*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_232 = -xi_160*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025];
         const double xi_233 = -xi_161*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_234 = -xi_162*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026];
         const double xi_235 = -xi_163*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_217 + xi_218 + xi_219 + xi_220 + xi_221 + xi_222 + xi_223 + xi_224 + xi_225) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1023; ctr_2 < 1024; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_140 = _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_132 = -xi_239*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_133 = -xi_240*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_134 = -xi_241*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = -xi_242*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024];
         const double xi_136 = -xi_243*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_137 = -xi_244*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026];
         const double xi_138 = -xi_245*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025];
         const double xi_139 = -xi_246*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_132 + xi_133 + xi_134 + xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_11(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceRhs, double const * const _data_edge_stencil_at_edge_x, double const * const _data_edge_stencil_at_edge_xy, double const * const _data_edge_stencil_at_edge_y, double * RESTRICT _data_vertexFaceDst, double const * const _data_vertex_stencil_at_edge_x, double const * const _data_vertex_stencil_at_edge_xy, double const * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_80 = 1.0;
   const double xi_81 = -relax;
   const double xi_67 = _data_edge_stencil_at_edge_xy[0];
   const double xi_78 = 1 / (xi_67);
   const double xi_68 = _data_edge_stencil_at_edge_xy[1];
   const double xi_69 = _data_edge_stencil_at_edge_xy[4];
   const double xi_70 = _data_edge_stencil_at_edge_xy[3];
   const double xi_71 = _data_edge_stencil_at_edge_xy[2];
   const double xi_72 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_73 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_74 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_155 = _data_edge_stencil_at_edge_y[0];
   const double xi_167 = 1 / (xi_155);
   const double xi_156 = _data_edge_stencil_at_edge_y[4];
   const double xi_157 = _data_edge_stencil_at_edge_y[3];
   const double xi_158 = _data_edge_stencil_at_edge_y[1];
   const double xi_159 = _data_edge_stencil_at_edge_y[2];
   const double xi_160 = _data_vertex_stencil_at_edge_y[3];
   const double xi_161 = _data_vertex_stencil_at_edge_y[0];
   const double xi_162 = _data_vertex_stencil_at_edge_y[2];
   const double xi_163 = _data_vertex_stencil_at_edge_y[1];
   const double xi_238 = _data_edge_stencil_at_edge_x[0];
   const double xi_258 = 1 / (xi_238);
   const double xi_239 = _data_edge_stencil_at_edge_x[1];
   const double xi_240 = _data_edge_stencil_at_edge_x[3];
   const double xi_241 = _data_edge_stencil_at_edge_x[4];
   const double xi_242 = _data_edge_stencil_at_edge_x[2];
   const double xi_243 = _data_vertex_stencil_at_edge_x[1];
   const double xi_244 = _data_vertex_stencil_at_edge_x[3];
   const double xi_245 = _data_vertex_stencil_at_edge_x[2];
   const double xi_246 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_90 = _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_82 = -xi_68*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = -xi_69*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_84 = -xi_70*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049];
         const double xi_85 = -xi_71*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_86 = -xi_72*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = -xi_73*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050];
         const double xi_88 = -xi_74*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_89 = -xi_75*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2051];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
      {
         const double xi_179 = _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_190 = _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_171 = -xi_68*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_172 = -xi_69*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = -xi_70*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049];
         const double xi_174 = -xi_71*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = -xi_72*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_176 = -xi_73*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050];
         const double xi_177 = -xi_74*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = -xi_75*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2051];
         const double xi_182 = -xi_156*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_183 = -xi_157*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048];
         const double xi_184 = -xi_158*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_185 = -xi_159*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_186 = -xi_160*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049];
         const double xi_187 = -xi_161*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_188 = -xi_162*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050];
         const double xi_189 = -xi_163*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178 + xi_179) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 2047; ctr_1 < 2048; ctr_1 += 1)
      {
         const double xi_115 = _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = -xi_156*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_109 = -xi_157*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048];
         const double xi_110 = -xi_158*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = -xi_159*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_112 = -xi_160*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049];
         const double xi_113 = -xi_161*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_114 = -xi_162*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050];
         const double xi_107 = -xi_163*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 2047; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_271 = _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_282 = _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_263 = -xi_239*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_264 = -xi_240*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_265 = -xi_241*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_266 = -xi_242*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048];
         const double xi_267 = -xi_243*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_268 = -xi_244*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050];
         const double xi_269 = -xi_245*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_270 = -xi_246*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_274 = -xi_68*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_275 = -xi_69*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_276 = -xi_70*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049];
         const double xi_277 = -xi_71*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_278 = -xi_72*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_279 = -xi_73*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050];
         const double xi_280 = -xi_74*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_281 = -xi_75*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2051];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_263 + xi_264 + xi_265 + xi_266 + xi_267 + xi_268 + xi_269 + xi_270 + xi_271) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2047; ctr_1 += 1)
      {
         const double xi_43 = _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = -xi_239*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_36 = -xi_240*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_37 = -xi_241*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_242*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048];
         const double xi_39 = -xi_243*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_40 = -xi_244*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050];
         const double xi_41 = -xi_245*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_42 = -xi_246*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = -xi_68*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = -xi_69*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_48 = -xi_70*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049];
         const double xi_49 = -xi_71*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_50 = -xi_72*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_51 = -xi_73*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050];
         const double xi_52 = -xi_74*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_53 = -xi_75*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2051];
         const double xi_57 = -xi_156*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_58 = -xi_157*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048];
         const double xi_59 = -xi_158*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = -xi_159*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_61 = -xi_160*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049];
         const double xi_62 = -xi_161*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = -xi_162*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050];
         const double xi_64 = -xi_163*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 2047; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
      {
         const double xi_225 = _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_236 = _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_217 = -xi_239*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_218 = -xi_240*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_219 = -xi_241*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_220 = -xi_242*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048];
         const double xi_221 = -xi_243*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_222 = -xi_244*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050];
         const double xi_223 = -xi_245*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_224 = -xi_246*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_228 = -xi_156*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_229 = -xi_157*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048];
         const double xi_230 = -xi_158*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_231 = -xi_159*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_232 = -xi_160*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049];
         const double xi_233 = -xi_161*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_234 = -xi_162*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050];
         const double xi_235 = -xi_163*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_217 + xi_218 + xi_219 + xi_220 + xi_221 + xi_222 + xi_223 + xi_224 + xi_225) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 2047; ctr_2 < 2048; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_140 = _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_132 = -xi_239*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_133 = -xi_240*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_134 = -xi_241*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = -xi_242*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048];
         const double xi_136 = -xi_243*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_137 = -xi_244*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050];
         const double xi_138 = -xi_245*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049];
         const double xi_139 = -xi_246*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_132 + xi_133 + xi_134 + xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_12(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceRhs, double const * const _data_edge_stencil_at_edge_x, double const * const _data_edge_stencil_at_edge_xy, double const * const _data_edge_stencil_at_edge_y, double * RESTRICT _data_vertexFaceDst, double const * const _data_vertex_stencil_at_edge_x, double const * const _data_vertex_stencil_at_edge_xy, double const * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_80 = 1.0;
   const double xi_81 = -relax;
   const double xi_67 = _data_edge_stencil_at_edge_xy[0];
   const double xi_78 = 1 / (xi_67);
   const double xi_68 = _data_edge_stencil_at_edge_xy[1];
   const double xi_69 = _data_edge_stencil_at_edge_xy[4];
   const double xi_70 = _data_edge_stencil_at_edge_xy[3];
   const double xi_71 = _data_edge_stencil_at_edge_xy[2];
   const double xi_72 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_73 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_74 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_155 = _data_edge_stencil_at_edge_y[0];
   const double xi_167 = 1 / (xi_155);
   const double xi_156 = _data_edge_stencil_at_edge_y[4];
   const double xi_157 = _data_edge_stencil_at_edge_y[3];
   const double xi_158 = _data_edge_stencil_at_edge_y[1];
   const double xi_159 = _data_edge_stencil_at_edge_y[2];
   const double xi_160 = _data_vertex_stencil_at_edge_y[3];
   const double xi_161 = _data_vertex_stencil_at_edge_y[0];
   const double xi_162 = _data_vertex_stencil_at_edge_y[2];
   const double xi_163 = _data_vertex_stencil_at_edge_y[1];
   const double xi_238 = _data_edge_stencil_at_edge_x[0];
   const double xi_258 = 1 / (xi_238);
   const double xi_239 = _data_edge_stencil_at_edge_x[1];
   const double xi_240 = _data_edge_stencil_at_edge_x[3];
   const double xi_241 = _data_edge_stencil_at_edge_x[4];
   const double xi_242 = _data_edge_stencil_at_edge_x[2];
   const double xi_243 = _data_vertex_stencil_at_edge_x[1];
   const double xi_244 = _data_vertex_stencil_at_edge_x[3];
   const double xi_245 = _data_vertex_stencil_at_edge_x[2];
   const double xi_246 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_90 = _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_82 = -xi_68*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = -xi_69*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_84 = -xi_70*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097];
         const double xi_85 = -xi_71*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_86 = -xi_72*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = -xi_73*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098];
         const double xi_88 = -xi_74*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_89 = -xi_75*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4099];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
      {
         const double xi_179 = _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_190 = _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_171 = -xi_68*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_172 = -xi_69*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = -xi_70*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097];
         const double xi_174 = -xi_71*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = -xi_72*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_176 = -xi_73*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098];
         const double xi_177 = -xi_74*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = -xi_75*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4099];
         const double xi_182 = -xi_156*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_183 = -xi_157*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096];
         const double xi_184 = -xi_158*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_185 = -xi_159*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_186 = -xi_160*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097];
         const double xi_187 = -xi_161*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_188 = -xi_162*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098];
         const double xi_189 = -xi_163*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178 + xi_179) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 4095; ctr_1 < 4096; ctr_1 += 1)
      {
         const double xi_115 = _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = -xi_156*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_109 = -xi_157*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096];
         const double xi_110 = -xi_158*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = -xi_159*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_112 = -xi_160*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097];
         const double xi_113 = -xi_161*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_114 = -xi_162*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098];
         const double xi_107 = -xi_163*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 4095; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_271 = _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_282 = _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_263 = -xi_239*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_264 = -xi_240*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_265 = -xi_241*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_266 = -xi_242*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096];
         const double xi_267 = -xi_243*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_268 = -xi_244*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098];
         const double xi_269 = -xi_245*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_270 = -xi_246*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_274 = -xi_68*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_275 = -xi_69*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_276 = -xi_70*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097];
         const double xi_277 = -xi_71*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_278 = -xi_72*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_279 = -xi_73*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098];
         const double xi_280 = -xi_74*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_281 = -xi_75*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4099];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_263 + xi_264 + xi_265 + xi_266 + xi_267 + xi_268 + xi_269 + xi_270 + xi_271) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4095; ctr_1 += 1)
      {
         const double xi_43 = _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = -xi_239*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_36 = -xi_240*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_37 = -xi_241*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_242*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096];
         const double xi_39 = -xi_243*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_40 = -xi_244*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098];
         const double xi_41 = -xi_245*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_42 = -xi_246*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = -xi_68*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = -xi_69*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_48 = -xi_70*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097];
         const double xi_49 = -xi_71*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_50 = -xi_72*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_51 = -xi_73*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098];
         const double xi_52 = -xi_74*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_53 = -xi_75*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4099];
         const double xi_57 = -xi_156*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_58 = -xi_157*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096];
         const double xi_59 = -xi_158*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = -xi_159*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_61 = -xi_160*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097];
         const double xi_62 = -xi_161*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = -xi_162*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098];
         const double xi_64 = -xi_163*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 4095; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
      {
         const double xi_225 = _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_236 = _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_217 = -xi_239*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_218 = -xi_240*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_219 = -xi_241*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_220 = -xi_242*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096];
         const double xi_221 = -xi_243*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_222 = -xi_244*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098];
         const double xi_223 = -xi_245*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_224 = -xi_246*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_228 = -xi_156*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_229 = -xi_157*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096];
         const double xi_230 = -xi_158*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_231 = -xi_159*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_232 = -xi_160*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097];
         const double xi_233 = -xi_161*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_234 = -xi_162*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098];
         const double xi_235 = -xi_163*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_217 + xi_218 + xi_219 + xi_220 + xi_221 + xi_222 + xi_223 + xi_224 + xi_225) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 4095; ctr_2 < 4096; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_140 = _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_132 = -xi_239*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_133 = -xi_240*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_134 = -xi_241*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = -xi_242*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096];
         const double xi_136 = -xi_243*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_137 = -xi_244*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098];
         const double xi_138 = -xi_245*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097];
         const double xi_139 = -xi_246*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_132 + xi_133 + xi_134 + xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_13(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceRhs, double const * const _data_edge_stencil_at_edge_x, double const * const _data_edge_stencil_at_edge_xy, double const * const _data_edge_stencil_at_edge_y, double * RESTRICT _data_vertexFaceDst, double const * const _data_vertex_stencil_at_edge_x, double const * const _data_vertex_stencil_at_edge_xy, double const * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_80 = 1.0;
   const double xi_81 = -relax;
   const double xi_67 = _data_edge_stencil_at_edge_xy[0];
   const double xi_78 = 1 / (xi_67);
   const double xi_68 = _data_edge_stencil_at_edge_xy[1];
   const double xi_69 = _data_edge_stencil_at_edge_xy[4];
   const double xi_70 = _data_edge_stencil_at_edge_xy[3];
   const double xi_71 = _data_edge_stencil_at_edge_xy[2];
   const double xi_72 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_73 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_74 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_155 = _data_edge_stencil_at_edge_y[0];
   const double xi_167 = 1 / (xi_155);
   const double xi_156 = _data_edge_stencil_at_edge_y[4];
   const double xi_157 = _data_edge_stencil_at_edge_y[3];
   const double xi_158 = _data_edge_stencil_at_edge_y[1];
   const double xi_159 = _data_edge_stencil_at_edge_y[2];
   const double xi_160 = _data_vertex_stencil_at_edge_y[3];
   const double xi_161 = _data_vertex_stencil_at_edge_y[0];
   const double xi_162 = _data_vertex_stencil_at_edge_y[2];
   const double xi_163 = _data_vertex_stencil_at_edge_y[1];
   const double xi_238 = _data_edge_stencil_at_edge_x[0];
   const double xi_258 = 1 / (xi_238);
   const double xi_239 = _data_edge_stencil_at_edge_x[1];
   const double xi_240 = _data_edge_stencil_at_edge_x[3];
   const double xi_241 = _data_edge_stencil_at_edge_x[4];
   const double xi_242 = _data_edge_stencil_at_edge_x[2];
   const double xi_243 = _data_vertex_stencil_at_edge_x[1];
   const double xi_244 = _data_vertex_stencil_at_edge_x[3];
   const double xi_245 = _data_vertex_stencil_at_edge_x[2];
   const double xi_246 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_90 = _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_82 = -xi_68*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = -xi_69*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_84 = -xi_70*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193];
         const double xi_85 = -xi_71*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_86 = -xi_72*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = -xi_73*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194];
         const double xi_88 = -xi_74*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_89 = -xi_75*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8195];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 8191; ctr_1 += 1)
      {
         const double xi_179 = _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_190 = _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_171 = -xi_68*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_172 = -xi_69*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = -xi_70*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193];
         const double xi_174 = -xi_71*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = -xi_72*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_176 = -xi_73*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194];
         const double xi_177 = -xi_74*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = -xi_75*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8195];
         const double xi_182 = -xi_156*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_183 = -xi_157*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192];
         const double xi_184 = -xi_158*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_185 = -xi_159*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_186 = -xi_160*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193];
         const double xi_187 = -xi_161*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_188 = -xi_162*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194];
         const double xi_189 = -xi_163*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178 + xi_179) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 8191; ctr_1 < 8192; ctr_1 += 1)
      {
         const double xi_115 = _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = -xi_156*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_109 = -xi_157*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192];
         const double xi_110 = -xi_158*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = -xi_159*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_112 = -xi_160*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193];
         const double xi_113 = -xi_161*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_114 = -xi_162*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194];
         const double xi_107 = -xi_163*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 8191; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_271 = _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_282 = _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_263 = -xi_239*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_264 = -xi_240*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_265 = -xi_241*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_266 = -xi_242*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192];
         const double xi_267 = -xi_243*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_268 = -xi_244*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194];
         const double xi_269 = -xi_245*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_270 = -xi_246*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_274 = -xi_68*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_275 = -xi_69*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_276 = -xi_70*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193];
         const double xi_277 = -xi_71*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_278 = -xi_72*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_279 = -xi_73*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194];
         const double xi_280 = -xi_74*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_281 = -xi_75*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8195];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_263 + xi_264 + xi_265 + xi_266 + xi_267 + xi_268 + xi_269 + xi_270 + xi_271) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8191; ctr_1 += 1)
      {
         const double xi_43 = _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = -xi_239*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_36 = -xi_240*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_37 = -xi_241*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_242*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192];
         const double xi_39 = -xi_243*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_40 = -xi_244*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194];
         const double xi_41 = -xi_245*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_42 = -xi_246*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = -xi_68*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = -xi_69*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_48 = -xi_70*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193];
         const double xi_49 = -xi_71*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_50 = -xi_72*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_51 = -xi_73*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194];
         const double xi_52 = -xi_74*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_53 = -xi_75*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8195];
         const double xi_57 = -xi_156*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_58 = -xi_157*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192];
         const double xi_59 = -xi_158*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = -xi_159*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_61 = -xi_160*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193];
         const double xi_62 = -xi_161*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = -xi_162*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194];
         const double xi_64 = -xi_163*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 8191; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
      {
         const double xi_225 = _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_236 = _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_217 = -xi_239*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_218 = -xi_240*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_219 = -xi_241*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_220 = -xi_242*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192];
         const double xi_221 = -xi_243*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_222 = -xi_244*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194];
         const double xi_223 = -xi_245*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_224 = -xi_246*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_228 = -xi_156*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_229 = -xi_157*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192];
         const double xi_230 = -xi_158*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_231 = -xi_159*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_232 = -xi_160*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193];
         const double xi_233 = -xi_161*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_234 = -xi_162*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194];
         const double xi_235 = -xi_163*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_217 + xi_218 + xi_219 + xi_220 + xi_221 + xi_222 + xi_223 + xi_224 + xi_225) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 8191; ctr_2 < 8192; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_140 = _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_132 = -xi_239*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_133 = -xi_240*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_134 = -xi_241*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = -xi_242*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192];
         const double xi_136 = -xi_243*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_137 = -xi_244*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194];
         const double xi_138 = -xi_245*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193];
         const double xi_139 = -xi_246*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_132 + xi_133 + xi_134 + xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_14(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceRhs, double const * const _data_edge_stencil_at_edge_x, double const * const _data_edge_stencil_at_edge_xy, double const * const _data_edge_stencil_at_edge_y, double * RESTRICT _data_vertexFaceDst, double const * const _data_vertex_stencil_at_edge_x, double const * const _data_vertex_stencil_at_edge_xy, double const * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_80 = 1.0;
   const double xi_81 = -relax;
   const double xi_67 = _data_edge_stencil_at_edge_xy[0];
   const double xi_78 = 1 / (xi_67);
   const double xi_68 = _data_edge_stencil_at_edge_xy[1];
   const double xi_69 = _data_edge_stencil_at_edge_xy[4];
   const double xi_70 = _data_edge_stencil_at_edge_xy[3];
   const double xi_71 = _data_edge_stencil_at_edge_xy[2];
   const double xi_72 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_73 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_74 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_155 = _data_edge_stencil_at_edge_y[0];
   const double xi_167 = 1 / (xi_155);
   const double xi_156 = _data_edge_stencil_at_edge_y[4];
   const double xi_157 = _data_edge_stencil_at_edge_y[3];
   const double xi_158 = _data_edge_stencil_at_edge_y[1];
   const double xi_159 = _data_edge_stencil_at_edge_y[2];
   const double xi_160 = _data_vertex_stencil_at_edge_y[3];
   const double xi_161 = _data_vertex_stencil_at_edge_y[0];
   const double xi_162 = _data_vertex_stencil_at_edge_y[2];
   const double xi_163 = _data_vertex_stencil_at_edge_y[1];
   const double xi_238 = _data_edge_stencil_at_edge_x[0];
   const double xi_258 = 1 / (xi_238);
   const double xi_239 = _data_edge_stencil_at_edge_x[1];
   const double xi_240 = _data_edge_stencil_at_edge_x[3];
   const double xi_241 = _data_edge_stencil_at_edge_x[4];
   const double xi_242 = _data_edge_stencil_at_edge_x[2];
   const double xi_243 = _data_vertex_stencil_at_edge_x[1];
   const double xi_244 = _data_vertex_stencil_at_edge_x[3];
   const double xi_245 = _data_vertex_stencil_at_edge_x[2];
   const double xi_246 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_90 = _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_82 = -xi_68*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = -xi_69*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_84 = -xi_70*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385];
         const double xi_85 = -xi_71*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_86 = -xi_72*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = -xi_73*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386];
         const double xi_88 = -xi_74*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_89 = -xi_75*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16387];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
      {
         const double xi_179 = _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_190 = _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_171 = -xi_68*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_172 = -xi_69*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = -xi_70*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385];
         const double xi_174 = -xi_71*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_175 = -xi_72*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_176 = -xi_73*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386];
         const double xi_177 = -xi_74*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = -xi_75*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16387];
         const double xi_182 = -xi_156*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_183 = -xi_157*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384];
         const double xi_184 = -xi_158*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_185 = -xi_159*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_186 = -xi_160*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385];
         const double xi_187 = -xi_161*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_188 = -xi_162*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386];
         const double xi_189 = -xi_163*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178 + xi_179) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 16383; ctr_1 < 16384; ctr_1 += 1)
      {
         const double xi_115 = _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = -xi_156*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_109 = -xi_157*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384];
         const double xi_110 = -xi_158*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = -xi_159*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_112 = -xi_160*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385];
         const double xi_113 = -xi_161*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_114 = -xi_162*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386];
         const double xi_107 = -xi_163*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 16383; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_271 = _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_282 = _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_263 = -xi_239*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_264 = -xi_240*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_265 = -xi_241*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_266 = -xi_242*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384];
         const double xi_267 = -xi_243*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_268 = -xi_244*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386];
         const double xi_269 = -xi_245*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_270 = -xi_246*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_274 = -xi_68*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_275 = -xi_69*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_276 = -xi_70*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385];
         const double xi_277 = -xi_71*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_278 = -xi_72*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_279 = -xi_73*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386];
         const double xi_280 = -xi_74*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_281 = -xi_75*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16387];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_263 + xi_264 + xi_265 + xi_266 + xi_267 + xi_268 + xi_269 + xi_270 + xi_271) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16383; ctr_1 += 1)
      {
         const double xi_43 = _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_65 = _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_35 = -xi_239*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_36 = -xi_240*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_37 = -xi_241*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_38 = -xi_242*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384];
         const double xi_39 = -xi_243*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_40 = -xi_244*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386];
         const double xi_41 = -xi_245*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_42 = -xi_246*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = -xi_68*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = -xi_69*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_48 = -xi_70*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385];
         const double xi_49 = -xi_71*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_50 = -xi_72*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_51 = -xi_73*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386];
         const double xi_52 = -xi_74*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_53 = -xi_75*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16387];
         const double xi_57 = -xi_156*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_58 = -xi_157*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384];
         const double xi_59 = -xi_158*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = -xi_159*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_61 = -xi_160*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385];
         const double xi_62 = -xi_161*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = -xi_162*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386];
         const double xi_64 = -xi_163*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_78*(xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 16383; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
      {
         const double xi_225 = _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_236 = _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_217 = -xi_239*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_218 = -xi_240*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_219 = -xi_241*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_220 = -xi_242*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384];
         const double xi_221 = -xi_243*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_222 = -xi_244*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386];
         const double xi_223 = -xi_245*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_224 = -xi_246*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_228 = -xi_156*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1];
         const double xi_229 = -xi_157*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384];
         const double xi_230 = -xi_158*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_231 = -xi_159*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_232 = -xi_160*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385];
         const double xi_233 = -xi_161*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_234 = -xi_162*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386];
         const double xi_235 = -xi_163*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_217 + xi_218 + xi_219 + xi_220 + xi_221 + xi_222 + xi_223 + xi_224 + xi_225) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_167*(xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 16383; ctr_2 < 16384; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_140 = _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_132 = -xi_239*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_133 = -xi_240*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_134 = -xi_241*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_135 = -xi_242*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384];
         const double xi_136 = -xi_243*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_137 = -xi_244*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386];
         const double xi_138 = -xi_245*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385];
         const double xi_139 = -xi_246*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_132 + xi_133 + xi_134 + xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_any(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceRhs, double const * const _data_edge_stencil_at_edge_x, double const * const _data_edge_stencil_at_edge_xy, double const * const _data_edge_stencil_at_edge_y, double * RESTRICT _data_vertexFaceDst, double const * const _data_vertex_stencil_at_edge_x, double const * const _data_vertex_stencil_at_edge_xy, double const * const _data_vertex_stencil_at_edge_y, int64_t level, double relax)
{
   const double xi_80 = 1.0;
   const double xi_81 = -relax;
   const double xi_67 = _data_edge_stencil_at_edge_xy[0];
   const double xi_78 = 1 / (xi_67);
   const double xi_68 = _data_edge_stencil_at_edge_xy[1];
   const double xi_69 = _data_edge_stencil_at_edge_xy[4];
   const double xi_70 = _data_edge_stencil_at_edge_xy[3];
   const double xi_71 = _data_edge_stencil_at_edge_xy[2];
   const double xi_72 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_73 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_74 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_155 = _data_edge_stencil_at_edge_y[0];
   const double xi_167 = 1 / (xi_155);
   const double xi_156 = _data_edge_stencil_at_edge_y[4];
   const double xi_157 = _data_edge_stencil_at_edge_y[3];
   const double xi_158 = _data_edge_stencil_at_edge_y[1];
   const double xi_159 = _data_edge_stencil_at_edge_y[2];
   const double xi_160 = _data_vertex_stencil_at_edge_y[3];
   const double xi_161 = _data_vertex_stencil_at_edge_y[0];
   const double xi_162 = _data_vertex_stencil_at_edge_y[2];
   const double xi_163 = _data_vertex_stencil_at_edge_y[1];
   const double xi_238 = _data_edge_stencil_at_edge_x[0];
   const double xi_258 = 1 / (xi_238);
   const double xi_239 = _data_edge_stencil_at_edge_x[1];
   const double xi_240 = _data_edge_stencil_at_edge_x[3];
   const double xi_241 = _data_edge_stencil_at_edge_x[4];
   const double xi_242 = _data_edge_stencil_at_edge_x[2];
   const double xi_243 = _data_vertex_stencil_at_edge_x[1];
   const double xi_244 = _data_vertex_stencil_at_edge_x[3];
   const double xi_245 = _data_vertex_stencil_at_edge_x[2];
   const double xi_246 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_90 = _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_82 = -xi_68*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = -xi_69*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_84 = -xi_70*_data_edgeFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_85 = -xi_71*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
         const double xi_86 = -xi_72*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = -xi_73*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_88 = -xi_74*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_89 = -xi_75*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] = relax*xi_78*(xi_82 + xi_83 + xi_84 + xi_85 + xi_86 + xi_87 + xi_88 + xi_89 + xi_90) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
      {
         const double xi_179 = _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_190 = _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_171 = -xi_68*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_172 = -xi_69*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_173 = -xi_70*_data_edgeFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_174 = -xi_71*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
         const double xi_175 = -xi_72*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_176 = -xi_73*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_177 = -xi_74*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_178 = -xi_75*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1];
         const double xi_182 = -xi_156*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
         const double xi_183 = -xi_157*_data_edgeFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_184 = -xi_158*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_185 = -xi_159*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_186 = -xi_160*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_187 = -xi_161*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_188 = -xi_162*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_189 = -xi_163*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] = relax*xi_78*(xi_171 + xi_172 + xi_173 + xi_174 + xi_175 + xi_176 + xi_177 + xi_178 + xi_179) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] = relax*xi_167*(xi_182 + xi_183 + xi_184 + xi_185 + xi_186 + xi_187 + xi_188 + xi_189 + xi_190) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
      {
         const double xi_115 = _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_108 = -xi_156*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
         const double xi_109 = -xi_157*_data_edgeFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_110 = -xi_158*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = -xi_159*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_112 = -xi_160*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_113 = -xi_161*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_114 = -xi_162*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_107 = -xi_163*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] = relax*xi_167*(xi_107 + xi_108 + xi_109 + xi_110 + xi_111 + xi_112 + xi_113 + xi_114 + xi_115) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_271 = _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_282 = _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_263 = -xi_239*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_264 = -xi_240*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_265 = -xi_241*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_266 = -xi_242*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
         const double xi_267 = -xi_243*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_268 = -xi_244*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_269 = -xi_245*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_270 = -xi_246*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_274 = -xi_68*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_275 = -xi_69*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_276 = -xi_70*_data_edgeFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_277 = -xi_71*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
         const double xi_278 = -xi_72*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_279 = -xi_73*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_280 = -xi_74*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_281 = -xi_75*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_263 + xi_264 + xi_265 + xi_266 + xi_267 + xi_268 + xi_269 + xi_270 + xi_271) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] = relax*xi_78*(xi_274 + xi_275 + xi_276 + xi_277 + xi_278 + xi_279 + xi_280 + xi_281 + xi_282) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
      {
         const double xi_43 = _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_54 = _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_65 = _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_35 = -xi_239*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_36 = -xi_240*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_37 = -xi_241*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_38 = -xi_242*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
         const double xi_39 = -xi_243*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_40 = -xi_244*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_41 = -xi_245*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_42 = -xi_246*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_46 = -xi_68*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_47 = -xi_69*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_48 = -xi_70*_data_edgeFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_49 = -xi_71*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
         const double xi_50 = -xi_72*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_51 = -xi_73*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_52 = -xi_74*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_53 = -xi_75*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1];
         const double xi_57 = -xi_156*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
         const double xi_58 = -xi_157*_data_edgeFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_59 = -xi_158*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_60 = -xi_159*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_61 = -xi_160*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_62 = -xi_161*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_63 = -xi_162*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_64 = -xi_163*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_35 + xi_36 + xi_37 + xi_38 + xi_39 + xi_40 + xi_41 + xi_42 + xi_43) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] = relax*xi_78*(xi_46 + xi_47 + xi_48 + xi_49 + xi_50 + xi_51 + xi_52 + xi_53 + xi_54) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] = relax*xi_167*(xi_57 + xi_58 + xi_59 + xi_60 + xi_61 + xi_62 + xi_63 + xi_64 + xi_65) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         const double xi_225 = _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_236 = _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_217 = -xi_239*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_218 = -xi_240*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_219 = -xi_241*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_220 = -xi_242*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
         const double xi_221 = -xi_243*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_222 = -xi_244*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_223 = -xi_245*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_224 = -xi_246*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         const double xi_228 = -xi_156*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1];
         const double xi_229 = -xi_157*_data_edgeFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_230 = -xi_158*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_231 = -xi_159*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_232 = -xi_160*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1];
         const double xi_233 = -xi_161*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_234 = -xi_162*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_235 = -xi_163*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_217 + xi_218 + xi_219 + xi_220 + xi_221 + xi_222 + xi_223 + xi_224 + xi_225) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] = relax*xi_167*(xi_228 + xi_229 + xi_230 + xi_231 + xi_232 + xi_233 + xi_234 + xi_235 + xi_236) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
      }
   }
   for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_140 = _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_132 = -xi_239*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_133 = -xi_240*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_134 = -xi_241*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
         const double xi_135 = -xi_242*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1];
         const double xi_136 = -xi_243*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_137 = -xi_244*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))];
         const double xi_138 = -xi_245*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1];
         const double xi_139 = -xi_246*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_258*(xi_132 + xi_133 + xi_134 + xi_135 + xi_136 + xi_137 + xi_138 + xi_139 + xi_140) + (xi_80 + xi_81)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}


void sor_2D_macroface_P2_update_edgedofs(double * RESTRICT _data_edgeFaceDst, double * RESTRICT _data_edgeFaceRhs, double const * const _data_edge_stencil_at_edge_x, double const * const _data_edge_stencil_at_edge_xy, double const * const _data_edge_stencil_at_edge_y, double * RESTRICT _data_vertexFaceDst, double const * const _data_vertex_stencil_at_edge_x, double const * const _data_vertex_stencil_at_edge_xy, double const * const _data_vertex_stencil_at_edge_y, int64_t level, double relax)
{
    switch( level )
    {
    case 2:
        sor_2D_macroface_P2_update_edgedofs_level_2(_data_edgeFaceDst, _data_edgeFaceRhs, _data_edge_stencil_at_edge_x, _data_edge_stencil_at_edge_xy, _data_edge_stencil_at_edge_y, _data_vertexFaceDst, _data_vertex_stencil_at_edge_x, _data_vertex_stencil_at_edge_xy, _data_vertex_stencil_at_edge_y, relax);
        break;
    case 3:
        sor_2D_macroface_P2_update_edgedofs_level_3(_data_edgeFaceDst, _data_edgeFaceRhs, _data_edge_stencil_at_edge_x, _data_edge_stencil_at_edge_xy, _data_edge_stencil_at_edge_y, _data_vertexFaceDst, _data_vertex_stencil_at_edge_x, _data_vertex_stencil_at_edge_xy, _data_vertex_stencil_at_edge_y, relax);
        break;
    case 4:
        sor_2D_macroface_P2_update_edgedofs_level_4(_data_edgeFaceDst, _data_edgeFaceRhs, _data_edge_stencil_at_edge_x, _data_edge_stencil_at_edge_xy, _data_edge_stencil_at_edge_y, _data_vertexFaceDst, _data_vertex_stencil_at_edge_x, _data_vertex_stencil_at_edge_xy, _data_vertex_stencil_at_edge_y, relax);
        break;
    case 5:
        sor_2D_macroface_P2_update_edgedofs_level_5(_data_edgeFaceDst, _data_edgeFaceRhs, _data_edge_stencil_at_edge_x, _data_edge_stencil_at_edge_xy, _data_edge_stencil_at_edge_y, _data_vertexFaceDst, _data_vertex_stencil_at_edge_x, _data_vertex_stencil_at_edge_xy, _data_vertex_stencil_at_edge_y, relax);
        break;
    case 6:
        sor_2D_macroface_P2_update_edgedofs_level_6(_data_edgeFaceDst, _data_edgeFaceRhs, _data_edge_stencil_at_edge_x, _data_edge_stencil_at_edge_xy, _data_edge_stencil_at_edge_y, _data_vertexFaceDst, _data_vertex_stencil_at_edge_x, _data_vertex_stencil_at_edge_xy, _data_vertex_stencil_at_edge_y, relax);
        break;
    case 7:
        sor_2D_macroface_P2_update_edgedofs_level_7(_data_edgeFaceDst, _data_edgeFaceRhs, _data_edge_stencil_at_edge_x, _data_edge_stencil_at_edge_xy, _data_edge_stencil_at_edge_y, _data_vertexFaceDst, _data_vertex_stencil_at_edge_x, _data_vertex_stencil_at_edge_xy, _data_vertex_stencil_at_edge_y, relax);
        break;
    case 8:
        sor_2D_macroface_P2_update_edgedofs_level_8(_data_edgeFaceDst, _data_edgeFaceRhs, _data_edge_stencil_at_edge_x, _data_edge_stencil_at_edge_xy, _data_edge_stencil_at_edge_y, _data_vertexFaceDst, _data_vertex_stencil_at_edge_x, _data_vertex_stencil_at_edge_xy, _data_vertex_stencil_at_edge_y, relax);
        break;
    case 9:
        sor_2D_macroface_P2_update_edgedofs_level_9(_data_edgeFaceDst, _data_edgeFaceRhs, _data_edge_stencil_at_edge_x, _data_edge_stencil_at_edge_xy, _data_edge_stencil_at_edge_y, _data_vertexFaceDst, _data_vertex_stencil_at_edge_x, _data_vertex_stencil_at_edge_xy, _data_vertex_stencil_at_edge_y, relax);
        break;
    case 10:
        sor_2D_macroface_P2_update_edgedofs_level_10(_data_edgeFaceDst, _data_edgeFaceRhs, _data_edge_stencil_at_edge_x, _data_edge_stencil_at_edge_xy, _data_edge_stencil_at_edge_y, _data_vertexFaceDst, _data_vertex_stencil_at_edge_x, _data_vertex_stencil_at_edge_xy, _data_vertex_stencil_at_edge_y, relax);
        break;
    case 11:
        sor_2D_macroface_P2_update_edgedofs_level_11(_data_edgeFaceDst, _data_edgeFaceRhs, _data_edge_stencil_at_edge_x, _data_edge_stencil_at_edge_xy, _data_edge_stencil_at_edge_y, _data_vertexFaceDst, _data_vertex_stencil_at_edge_x, _data_vertex_stencil_at_edge_xy, _data_vertex_stencil_at_edge_y, relax);
        break;
    case 12:
        sor_2D_macroface_P2_update_edgedofs_level_12(_data_edgeFaceDst, _data_edgeFaceRhs, _data_edge_stencil_at_edge_x, _data_edge_stencil_at_edge_xy, _data_edge_stencil_at_edge_y, _data_vertexFaceDst, _data_vertex_stencil_at_edge_x, _data_vertex_stencil_at_edge_xy, _data_vertex_stencil_at_edge_y, relax);
        break;
    case 13:
        sor_2D_macroface_P2_update_edgedofs_level_13(_data_edgeFaceDst, _data_edgeFaceRhs, _data_edge_stencil_at_edge_x, _data_edge_stencil_at_edge_xy, _data_edge_stencil_at_edge_y, _data_vertexFaceDst, _data_vertex_stencil_at_edge_x, _data_vertex_stencil_at_edge_xy, _data_vertex_stencil_at_edge_y, relax);
        break;
    case 14:
        sor_2D_macroface_P2_update_edgedofs_level_14(_data_edgeFaceDst, _data_edgeFaceRhs, _data_edge_stencil_at_edge_x, _data_edge_stencil_at_edge_xy, _data_edge_stencil_at_edge_y, _data_vertexFaceDst, _data_vertex_stencil_at_edge_x, _data_vertex_stencil_at_edge_xy, _data_vertex_stencil_at_edge_y, relax);
        break;
    default:
        sor_2D_macroface_P2_update_edgedofs_level_any(_data_edgeFaceDst, _data_edgeFaceRhs, _data_edge_stencil_at_edge_x, _data_edge_stencil_at_edge_xy, _data_edge_stencil_at_edge_y, _data_vertexFaceDst, _data_vertex_stencil_at_edge_x, _data_vertex_stencil_at_edge_xy, _data_vertex_stencil_at_edge_y, level, relax);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hhg

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsP2MacroFace2D.hpp"

namespace hhg {
namespace P2 {
namespace macroface {
namespace generated {

static void sor_2D_macroface_P2_update_edgedofs_level_2(double * _data_edgeFaceDst, double * _data_edgeFaceRhs, double * const _data_edge_stencil_at_edge_x, double * const _data_edge_stencil_at_edge_xy, double * const _data_edge_stencil_at_edge_y, double * _data_vertexFaceDst, double * const _data_vertex_stencil_at_edge_x, double * const _data_vertex_stencil_at_edge_xy, double * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_32 = _data_edge_stencil_at_edge_xy[0];
   const double xi_42 = 1 / (xi_32);
   const double xi_33 = _data_edge_stencil_at_edge_xy[1];
   const double xi_34 = _data_edge_stencil_at_edge_xy[4];
   const double xi_35 = _data_edge_stencil_at_edge_xy[3];
   const double xi_36 = _data_edge_stencil_at_edge_xy[2];
   const double xi_37 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_38 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_39 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_40 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_70 = _data_edge_stencil_at_edge_xy[0];
   const double xi_89 = 1 / (xi_70);
   const double xi_71 = _data_edge_stencil_at_edge_xy[1];
   const double xi_72 = _data_edge_stencil_at_edge_xy[4];
   const double xi_73 = _data_edge_stencil_at_edge_xy[3];
   const double xi_74 = _data_edge_stencil_at_edge_xy[2];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_76 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_77 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_78 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_79 = _data_edge_stencil_at_edge_y[0];
   const double xi_90 = 1 / (xi_79);
   const double xi_80 = _data_edge_stencil_at_edge_y[4];
   const double xi_81 = _data_edge_stencil_at_edge_y[3];
   const double xi_82 = _data_edge_stencil_at_edge_y[1];
   const double xi_83 = _data_edge_stencil_at_edge_y[2];
   const double xi_84 = _data_vertex_stencil_at_edge_y[3];
   const double xi_85 = _data_vertex_stencil_at_edge_y[0];
   const double xi_86 = _data_vertex_stencil_at_edge_y[2];
   const double xi_87 = _data_vertex_stencil_at_edge_y[1];
   const double xi_44 = _data_edge_stencil_at_edge_y[0];
   const double xi_54 = 1 / (xi_44);
   const double xi_45 = _data_edge_stencil_at_edge_y[4];
   const double xi_46 = _data_edge_stencil_at_edge_y[3];
   const double xi_47 = _data_edge_stencil_at_edge_y[1];
   const double xi_48 = _data_edge_stencil_at_edge_y[2];
   const double xi_49 = _data_vertex_stencil_at_edge_y[3];
   const double xi_50 = _data_vertex_stencil_at_edge_y[0];
   const double xi_51 = _data_vertex_stencil_at_edge_y[2];
   const double xi_52 = _data_vertex_stencil_at_edge_y[1];
   const double xi_114 = _data_edge_stencil_at_edge_x[0];
   const double xi_133 = 1 / (xi_114);
   const double xi_115 = _data_edge_stencil_at_edge_x[1];
   const double xi_116 = _data_edge_stencil_at_edge_x[3];
   const double xi_117 = _data_edge_stencil_at_edge_x[4];
   const double xi_118 = _data_edge_stencil_at_edge_x[2];
   const double xi_119 = _data_vertex_stencil_at_edge_x[1];
   const double xi_120 = _data_vertex_stencil_at_edge_x[3];
   const double xi_121 = _data_vertex_stencil_at_edge_x[2];
   const double xi_122 = _data_vertex_stencil_at_edge_x[0];
   const double xi_123 = _data_edge_stencil_at_edge_xy[0];
   const double xi_134 = 1 / (xi_123);
   const double xi_124 = _data_edge_stencil_at_edge_xy[1];
   const double xi_125 = _data_edge_stencil_at_edge_xy[4];
   const double xi_126 = _data_edge_stencil_at_edge_xy[3];
   const double xi_127 = _data_edge_stencil_at_edge_xy[2];
   const double xi_128 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_129 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_130 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_131 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_0 = _data_edge_stencil_at_edge_x[0];
   const double xi_28 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_edge_x[1];
   const double xi_2 = _data_edge_stencil_at_edge_x[3];
   const double xi_3 = _data_edge_stencil_at_edge_x[4];
   const double xi_4 = _data_edge_stencil_at_edge_x[2];
   const double xi_5 = _data_vertex_stencil_at_edge_x[1];
   const double xi_6 = _data_vertex_stencil_at_edge_x[3];
   const double xi_7 = _data_vertex_stencil_at_edge_x[2];
   const double xi_8 = _data_vertex_stencil_at_edge_x[0];
   const double xi_9 = _data_edge_stencil_at_edge_xy[0];
   const double xi_29 = 1 / (xi_9);
   const double xi_10 = _data_edge_stencil_at_edge_xy[1];
   const double xi_11 = _data_edge_stencil_at_edge_xy[4];
   const double xi_12 = _data_edge_stencil_at_edge_xy[3];
   const double xi_13 = _data_edge_stencil_at_edge_xy[2];
   const double xi_14 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_15 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_16 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_17 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_18 = _data_edge_stencil_at_edge_y[0];
   const double xi_30 = 1 / (xi_18);
   const double xi_19 = _data_edge_stencil_at_edge_y[4];
   const double xi_20 = _data_edge_stencil_at_edge_y[3];
   const double xi_21 = _data_edge_stencil_at_edge_y[1];
   const double xi_22 = _data_edge_stencil_at_edge_y[2];
   const double xi_23 = _data_vertex_stencil_at_edge_y[3];
   const double xi_24 = _data_vertex_stencil_at_edge_y[0];
   const double xi_25 = _data_vertex_stencil_at_edge_y[2];
   const double xi_26 = _data_vertex_stencil_at_edge_y[1];
   const double xi_92 = _data_edge_stencil_at_edge_x[0];
   const double xi_111 = 1 / (xi_92);
   const double xi_93 = _data_edge_stencil_at_edge_x[1];
   const double xi_94 = _data_edge_stencil_at_edge_x[3];
   const double xi_95 = _data_edge_stencil_at_edge_x[4];
   const double xi_96 = _data_edge_stencil_at_edge_x[2];
   const double xi_97 = _data_vertex_stencil_at_edge_x[1];
   const double xi_98 = _data_vertex_stencil_at_edge_x[3];
   const double xi_99 = _data_vertex_stencil_at_edge_x[2];
   const double xi_100 = _data_vertex_stencil_at_edge_x[0];
   const double xi_101 = _data_edge_stencil_at_edge_y[0];
   const double xi_112 = 1 / (xi_101);
   const double xi_102 = _data_edge_stencil_at_edge_y[4];
   const double xi_103 = _data_edge_stencil_at_edge_y[3];
   const double xi_104 = _data_edge_stencil_at_edge_y[1];
   const double xi_105 = _data_edge_stencil_at_edge_y[2];
   const double xi_106 = _data_vertex_stencil_at_edge_y[3];
   const double xi_107 = _data_vertex_stencil_at_edge_y[0];
   const double xi_108 = _data_vertex_stencil_at_edge_y[2];
   const double xi_109 = _data_vertex_stencil_at_edge_y[1];
   const double xi_56 = _data_edge_stencil_at_edge_x[0];
   const double xi_66 = 1 / (xi_56);
   const double xi_57 = _data_edge_stencil_at_edge_x[1];
   const double xi_58 = _data_edge_stencil_at_edge_x[3];
   const double xi_59 = _data_edge_stencil_at_edge_x[4];
   const double xi_60 = _data_edge_stencil_at_edge_x[2];
   const double xi_61 = _data_vertex_stencil_at_edge_x[1];
   const double xi_62 = _data_vertex_stencil_at_edge_x[3];
   const double xi_63 = _data_vertex_stencil_at_edge_x[2];
   const double xi_64 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_42*(-xi_33*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_34*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_35*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] - xi_36*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6] - xi_39*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_40*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 7] + _data_edgeFaceRhs[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 3; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_89*(-xi_71*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_72*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_73*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] - xi_74*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_75*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_76*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6] - xi_77*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_78*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 7] + _data_edgeFaceRhs[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_90*(-xi_80*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_81*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4] - xi_82*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_83*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_84*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] - xi_85*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_86*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6] - xi_87*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 3; ctr_1 < 4; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_54*(-xi_45*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_46*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4] - xi_47*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_48*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_49*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] - xi_50*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_51*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6] - xi_52*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 3; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_133*(-xi_115*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] - xi_116*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_117*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_118*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4] - xi_119*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_120*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6] - xi_121*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] - xi_122*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_134*(-xi_124*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_125*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_126*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] - xi_127*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_128*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_129*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6] - xi_130*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_131*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 7] + _data_edgeFaceRhs[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 3; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_28*(-xi_1*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] - xi_2*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4] - xi_5*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6] - xi_7*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] - xi_8*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_29*(-xi_10*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] - xi_13*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_14*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6] - xi_16*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_17*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 7] + _data_edgeFaceRhs[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_30*(-xi_19*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4] - xi_21*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] - xi_24*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6] - xi_26*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 3; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_111*(-xi_100*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_93*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] - xi_94*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_95*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_96*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4] - xi_97*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_98*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6] - xi_99*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] + _data_edgeFaceRhs[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_112*(-xi_102*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_103*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4] - xi_104*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_105*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_106*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 5] - xi_107*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_108*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6] - xi_109*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 3; ctr_2 < 4; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_66*(-xi_57*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] - xi_58*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + ((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_59*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_60*_data_edgeFaceDst[ctr_1 + 5*ctr_2 + 2*((20) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4] - xi_61*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_62*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 6] - xi_63*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 5] - xi_64*_data_vertexFaceDst[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 5*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_3(double * _data_edgeFaceDst, double * _data_edgeFaceRhs, double * const _data_edge_stencil_at_edge_x, double * const _data_edge_stencil_at_edge_xy, double * const _data_edge_stencil_at_edge_y, double * _data_vertexFaceDst, double * const _data_vertex_stencil_at_edge_x, double * const _data_vertex_stencil_at_edge_xy, double * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_32 = _data_edge_stencil_at_edge_xy[0];
   const double xi_42 = 1 / (xi_32);
   const double xi_33 = _data_edge_stencil_at_edge_xy[1];
   const double xi_34 = _data_edge_stencil_at_edge_xy[4];
   const double xi_35 = _data_edge_stencil_at_edge_xy[3];
   const double xi_36 = _data_edge_stencil_at_edge_xy[2];
   const double xi_37 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_38 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_39 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_40 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_70 = _data_edge_stencil_at_edge_xy[0];
   const double xi_89 = 1 / (xi_70);
   const double xi_71 = _data_edge_stencil_at_edge_xy[1];
   const double xi_72 = _data_edge_stencil_at_edge_xy[4];
   const double xi_73 = _data_edge_stencil_at_edge_xy[3];
   const double xi_74 = _data_edge_stencil_at_edge_xy[2];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_76 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_77 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_78 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_79 = _data_edge_stencil_at_edge_y[0];
   const double xi_90 = 1 / (xi_79);
   const double xi_80 = _data_edge_stencil_at_edge_y[4];
   const double xi_81 = _data_edge_stencil_at_edge_y[3];
   const double xi_82 = _data_edge_stencil_at_edge_y[1];
   const double xi_83 = _data_edge_stencil_at_edge_y[2];
   const double xi_84 = _data_vertex_stencil_at_edge_y[3];
   const double xi_85 = _data_vertex_stencil_at_edge_y[0];
   const double xi_86 = _data_vertex_stencil_at_edge_y[2];
   const double xi_87 = _data_vertex_stencil_at_edge_y[1];
   const double xi_44 = _data_edge_stencil_at_edge_y[0];
   const double xi_54 = 1 / (xi_44);
   const double xi_45 = _data_edge_stencil_at_edge_y[4];
   const double xi_46 = _data_edge_stencil_at_edge_y[3];
   const double xi_47 = _data_edge_stencil_at_edge_y[1];
   const double xi_48 = _data_edge_stencil_at_edge_y[2];
   const double xi_49 = _data_vertex_stencil_at_edge_y[3];
   const double xi_50 = _data_vertex_stencil_at_edge_y[0];
   const double xi_51 = _data_vertex_stencil_at_edge_y[2];
   const double xi_52 = _data_vertex_stencil_at_edge_y[1];
   const double xi_114 = _data_edge_stencil_at_edge_x[0];
   const double xi_133 = 1 / (xi_114);
   const double xi_115 = _data_edge_stencil_at_edge_x[1];
   const double xi_116 = _data_edge_stencil_at_edge_x[3];
   const double xi_117 = _data_edge_stencil_at_edge_x[4];
   const double xi_118 = _data_edge_stencil_at_edge_x[2];
   const double xi_119 = _data_vertex_stencil_at_edge_x[1];
   const double xi_120 = _data_vertex_stencil_at_edge_x[3];
   const double xi_121 = _data_vertex_stencil_at_edge_x[2];
   const double xi_122 = _data_vertex_stencil_at_edge_x[0];
   const double xi_123 = _data_edge_stencil_at_edge_xy[0];
   const double xi_134 = 1 / (xi_123);
   const double xi_124 = _data_edge_stencil_at_edge_xy[1];
   const double xi_125 = _data_edge_stencil_at_edge_xy[4];
   const double xi_126 = _data_edge_stencil_at_edge_xy[3];
   const double xi_127 = _data_edge_stencil_at_edge_xy[2];
   const double xi_128 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_129 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_130 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_131 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_0 = _data_edge_stencil_at_edge_x[0];
   const double xi_28 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_edge_x[1];
   const double xi_2 = _data_edge_stencil_at_edge_x[3];
   const double xi_3 = _data_edge_stencil_at_edge_x[4];
   const double xi_4 = _data_edge_stencil_at_edge_x[2];
   const double xi_5 = _data_vertex_stencil_at_edge_x[1];
   const double xi_6 = _data_vertex_stencil_at_edge_x[3];
   const double xi_7 = _data_vertex_stencil_at_edge_x[2];
   const double xi_8 = _data_vertex_stencil_at_edge_x[0];
   const double xi_9 = _data_edge_stencil_at_edge_xy[0];
   const double xi_29 = 1 / (xi_9);
   const double xi_10 = _data_edge_stencil_at_edge_xy[1];
   const double xi_11 = _data_edge_stencil_at_edge_xy[4];
   const double xi_12 = _data_edge_stencil_at_edge_xy[3];
   const double xi_13 = _data_edge_stencil_at_edge_xy[2];
   const double xi_14 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_15 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_16 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_17 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_18 = _data_edge_stencil_at_edge_y[0];
   const double xi_30 = 1 / (xi_18);
   const double xi_19 = _data_edge_stencil_at_edge_y[4];
   const double xi_20 = _data_edge_stencil_at_edge_y[3];
   const double xi_21 = _data_edge_stencil_at_edge_y[1];
   const double xi_22 = _data_edge_stencil_at_edge_y[2];
   const double xi_23 = _data_vertex_stencil_at_edge_y[3];
   const double xi_24 = _data_vertex_stencil_at_edge_y[0];
   const double xi_25 = _data_vertex_stencil_at_edge_y[2];
   const double xi_26 = _data_vertex_stencil_at_edge_y[1];
   const double xi_92 = _data_edge_stencil_at_edge_x[0];
   const double xi_111 = 1 / (xi_92);
   const double xi_93 = _data_edge_stencil_at_edge_x[1];
   const double xi_94 = _data_edge_stencil_at_edge_x[3];
   const double xi_95 = _data_edge_stencil_at_edge_x[4];
   const double xi_96 = _data_edge_stencil_at_edge_x[2];
   const double xi_97 = _data_vertex_stencil_at_edge_x[1];
   const double xi_98 = _data_vertex_stencil_at_edge_x[3];
   const double xi_99 = _data_vertex_stencil_at_edge_x[2];
   const double xi_100 = _data_vertex_stencil_at_edge_x[0];
   const double xi_101 = _data_edge_stencil_at_edge_y[0];
   const double xi_112 = 1 / (xi_101);
   const double xi_102 = _data_edge_stencil_at_edge_y[4];
   const double xi_103 = _data_edge_stencil_at_edge_y[3];
   const double xi_104 = _data_edge_stencil_at_edge_y[1];
   const double xi_105 = _data_edge_stencil_at_edge_y[2];
   const double xi_106 = _data_vertex_stencil_at_edge_y[3];
   const double xi_107 = _data_vertex_stencil_at_edge_y[0];
   const double xi_108 = _data_vertex_stencil_at_edge_y[2];
   const double xi_109 = _data_vertex_stencil_at_edge_y[1];
   const double xi_56 = _data_edge_stencil_at_edge_x[0];
   const double xi_66 = 1 / (xi_56);
   const double xi_57 = _data_edge_stencil_at_edge_x[1];
   const double xi_58 = _data_edge_stencil_at_edge_x[3];
   const double xi_59 = _data_edge_stencil_at_edge_x[4];
   const double xi_60 = _data_edge_stencil_at_edge_x[2];
   const double xi_61 = _data_vertex_stencil_at_edge_x[1];
   const double xi_62 = _data_vertex_stencil_at_edge_x[3];
   const double xi_63 = _data_vertex_stencil_at_edge_x[2];
   const double xi_64 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_42*(-xi_33*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_34*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_35*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] - xi_36*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10] - xi_39*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_40*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 11] + _data_edgeFaceRhs[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 7; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_89*(-xi_71*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_72*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_73*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] - xi_74*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_75*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_76*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10] - xi_77*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_78*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 11] + _data_edgeFaceRhs[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_90*(-xi_80*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_81*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8] - xi_82*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_83*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_84*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] - xi_85*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_86*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10] - xi_87*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 7; ctr_1 < 8; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_54*(-xi_45*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_46*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8] - xi_47*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_48*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_49*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] - xi_50*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_51*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10] - xi_52*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 7; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_133*(-xi_115*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] - xi_116*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_117*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_118*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8] - xi_119*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_120*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10] - xi_121*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] - xi_122*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_134*(-xi_124*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_125*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_126*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] - xi_127*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_128*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_129*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10] - xi_130*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_131*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 11] + _data_edgeFaceRhs[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 7; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_28*(-xi_1*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] - xi_2*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8] - xi_5*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10] - xi_7*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] - xi_8*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_29*(-xi_10*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] - xi_13*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_14*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10] - xi_16*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_17*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 11] + _data_edgeFaceRhs[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_30*(-xi_19*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8] - xi_21*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] - xi_24*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10] - xi_26*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 7; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_111*(-xi_100*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_93*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] - xi_94*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_95*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_96*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8] - xi_97*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_98*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10] - xi_99*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] + _data_edgeFaceRhs[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_112*(-xi_102*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_103*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8] - xi_104*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_105*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_106*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 9] - xi_107*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_108*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10] - xi_109*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 7; ctr_2 < 8; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_66*(-xi_57*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] - xi_58*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + ((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_59*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_60*_data_edgeFaceDst[ctr_1 + 9*ctr_2 + 2*((72) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8] - xi_61*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_62*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 10] - xi_63*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 9] - xi_64*_data_vertexFaceDst[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 9*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_4(double * _data_edgeFaceDst, double * _data_edgeFaceRhs, double * const _data_edge_stencil_at_edge_x, double * const _data_edge_stencil_at_edge_xy, double * const _data_edge_stencil_at_edge_y, double * _data_vertexFaceDst, double * const _data_vertex_stencil_at_edge_x, double * const _data_vertex_stencil_at_edge_xy, double * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_32 = _data_edge_stencil_at_edge_xy[0];
   const double xi_42 = 1 / (xi_32);
   const double xi_33 = _data_edge_stencil_at_edge_xy[1];
   const double xi_34 = _data_edge_stencil_at_edge_xy[4];
   const double xi_35 = _data_edge_stencil_at_edge_xy[3];
   const double xi_36 = _data_edge_stencil_at_edge_xy[2];
   const double xi_37 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_38 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_39 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_40 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_70 = _data_edge_stencil_at_edge_xy[0];
   const double xi_89 = 1 / (xi_70);
   const double xi_71 = _data_edge_stencil_at_edge_xy[1];
   const double xi_72 = _data_edge_stencil_at_edge_xy[4];
   const double xi_73 = _data_edge_stencil_at_edge_xy[3];
   const double xi_74 = _data_edge_stencil_at_edge_xy[2];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_76 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_77 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_78 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_79 = _data_edge_stencil_at_edge_y[0];
   const double xi_90 = 1 / (xi_79);
   const double xi_80 = _data_edge_stencil_at_edge_y[4];
   const double xi_81 = _data_edge_stencil_at_edge_y[3];
   const double xi_82 = _data_edge_stencil_at_edge_y[1];
   const double xi_83 = _data_edge_stencil_at_edge_y[2];
   const double xi_84 = _data_vertex_stencil_at_edge_y[3];
   const double xi_85 = _data_vertex_stencil_at_edge_y[0];
   const double xi_86 = _data_vertex_stencil_at_edge_y[2];
   const double xi_87 = _data_vertex_stencil_at_edge_y[1];
   const double xi_44 = _data_edge_stencil_at_edge_y[0];
   const double xi_54 = 1 / (xi_44);
   const double xi_45 = _data_edge_stencil_at_edge_y[4];
   const double xi_46 = _data_edge_stencil_at_edge_y[3];
   const double xi_47 = _data_edge_stencil_at_edge_y[1];
   const double xi_48 = _data_edge_stencil_at_edge_y[2];
   const double xi_49 = _data_vertex_stencil_at_edge_y[3];
   const double xi_50 = _data_vertex_stencil_at_edge_y[0];
   const double xi_51 = _data_vertex_stencil_at_edge_y[2];
   const double xi_52 = _data_vertex_stencil_at_edge_y[1];
   const double xi_114 = _data_edge_stencil_at_edge_x[0];
   const double xi_133 = 1 / (xi_114);
   const double xi_115 = _data_edge_stencil_at_edge_x[1];
   const double xi_116 = _data_edge_stencil_at_edge_x[3];
   const double xi_117 = _data_edge_stencil_at_edge_x[4];
   const double xi_118 = _data_edge_stencil_at_edge_x[2];
   const double xi_119 = _data_vertex_stencil_at_edge_x[1];
   const double xi_120 = _data_vertex_stencil_at_edge_x[3];
   const double xi_121 = _data_vertex_stencil_at_edge_x[2];
   const double xi_122 = _data_vertex_stencil_at_edge_x[0];
   const double xi_123 = _data_edge_stencil_at_edge_xy[0];
   const double xi_134 = 1 / (xi_123);
   const double xi_124 = _data_edge_stencil_at_edge_xy[1];
   const double xi_125 = _data_edge_stencil_at_edge_xy[4];
   const double xi_126 = _data_edge_stencil_at_edge_xy[3];
   const double xi_127 = _data_edge_stencil_at_edge_xy[2];
   const double xi_128 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_129 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_130 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_131 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_0 = _data_edge_stencil_at_edge_x[0];
   const double xi_28 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_edge_x[1];
   const double xi_2 = _data_edge_stencil_at_edge_x[3];
   const double xi_3 = _data_edge_stencil_at_edge_x[4];
   const double xi_4 = _data_edge_stencil_at_edge_x[2];
   const double xi_5 = _data_vertex_stencil_at_edge_x[1];
   const double xi_6 = _data_vertex_stencil_at_edge_x[3];
   const double xi_7 = _data_vertex_stencil_at_edge_x[2];
   const double xi_8 = _data_vertex_stencil_at_edge_x[0];
   const double xi_9 = _data_edge_stencil_at_edge_xy[0];
   const double xi_29 = 1 / (xi_9);
   const double xi_10 = _data_edge_stencil_at_edge_xy[1];
   const double xi_11 = _data_edge_stencil_at_edge_xy[4];
   const double xi_12 = _data_edge_stencil_at_edge_xy[3];
   const double xi_13 = _data_edge_stencil_at_edge_xy[2];
   const double xi_14 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_15 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_16 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_17 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_18 = _data_edge_stencil_at_edge_y[0];
   const double xi_30 = 1 / (xi_18);
   const double xi_19 = _data_edge_stencil_at_edge_y[4];
   const double xi_20 = _data_edge_stencil_at_edge_y[3];
   const double xi_21 = _data_edge_stencil_at_edge_y[1];
   const double xi_22 = _data_edge_stencil_at_edge_y[2];
   const double xi_23 = _data_vertex_stencil_at_edge_y[3];
   const double xi_24 = _data_vertex_stencil_at_edge_y[0];
   const double xi_25 = _data_vertex_stencil_at_edge_y[2];
   const double xi_26 = _data_vertex_stencil_at_edge_y[1];
   const double xi_92 = _data_edge_stencil_at_edge_x[0];
   const double xi_111 = 1 / (xi_92);
   const double xi_93 = _data_edge_stencil_at_edge_x[1];
   const double xi_94 = _data_edge_stencil_at_edge_x[3];
   const double xi_95 = _data_edge_stencil_at_edge_x[4];
   const double xi_96 = _data_edge_stencil_at_edge_x[2];
   const double xi_97 = _data_vertex_stencil_at_edge_x[1];
   const double xi_98 = _data_vertex_stencil_at_edge_x[3];
   const double xi_99 = _data_vertex_stencil_at_edge_x[2];
   const double xi_100 = _data_vertex_stencil_at_edge_x[0];
   const double xi_101 = _data_edge_stencil_at_edge_y[0];
   const double xi_112 = 1 / (xi_101);
   const double xi_102 = _data_edge_stencil_at_edge_y[4];
   const double xi_103 = _data_edge_stencil_at_edge_y[3];
   const double xi_104 = _data_edge_stencil_at_edge_y[1];
   const double xi_105 = _data_edge_stencil_at_edge_y[2];
   const double xi_106 = _data_vertex_stencil_at_edge_y[3];
   const double xi_107 = _data_vertex_stencil_at_edge_y[0];
   const double xi_108 = _data_vertex_stencil_at_edge_y[2];
   const double xi_109 = _data_vertex_stencil_at_edge_y[1];
   const double xi_56 = _data_edge_stencil_at_edge_x[0];
   const double xi_66 = 1 / (xi_56);
   const double xi_57 = _data_edge_stencil_at_edge_x[1];
   const double xi_58 = _data_edge_stencil_at_edge_x[3];
   const double xi_59 = _data_edge_stencil_at_edge_x[4];
   const double xi_60 = _data_edge_stencil_at_edge_x[2];
   const double xi_61 = _data_vertex_stencil_at_edge_x[1];
   const double xi_62 = _data_vertex_stencil_at_edge_x[3];
   const double xi_63 = _data_vertex_stencil_at_edge_x[2];
   const double xi_64 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_42*(-xi_33*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_34*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_35*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] - xi_36*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18] - xi_39*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_40*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 19] + _data_edgeFaceRhs[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 15; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_89*(-xi_71*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_72*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_73*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] - xi_74*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_75*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_76*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18] - xi_77*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_78*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 19] + _data_edgeFaceRhs[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_90*(-xi_80*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_81*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16] - xi_82*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_83*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_84*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] - xi_85*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_86*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18] - xi_87*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 15; ctr_1 < 16; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_54*(-xi_45*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_46*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16] - xi_47*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_48*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_49*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] - xi_50*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_51*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18] - xi_52*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 15; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_133*(-xi_115*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] - xi_116*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_117*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_118*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16] - xi_119*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_120*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18] - xi_121*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] - xi_122*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_134*(-xi_124*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_125*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_126*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] - xi_127*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_128*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_129*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18] - xi_130*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_131*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 19] + _data_edgeFaceRhs[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 15; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_28*(-xi_1*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] - xi_2*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16] - xi_5*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18] - xi_7*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] - xi_8*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_29*(-xi_10*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] - xi_13*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_14*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18] - xi_16*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_17*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 19] + _data_edgeFaceRhs[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_30*(-xi_19*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16] - xi_21*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] - xi_24*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18] - xi_26*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 15; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_111*(-xi_100*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_93*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] - xi_94*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_95*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_96*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16] - xi_97*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_98*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18] - xi_99*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] + _data_edgeFaceRhs[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_112*(-xi_102*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_103*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16] - xi_104*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_105*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_106*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 17] - xi_107*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_108*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18] - xi_109*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 15; ctr_2 < 16; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_66*(-xi_57*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] - xi_58*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + ((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_59*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_60*_data_edgeFaceDst[ctr_1 + 17*ctr_2 + 2*((272) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16] - xi_61*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_62*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 18] - xi_63*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 17] - xi_64*_data_vertexFaceDst[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 17*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_5(double * _data_edgeFaceDst, double * _data_edgeFaceRhs, double * const _data_edge_stencil_at_edge_x, double * const _data_edge_stencil_at_edge_xy, double * const _data_edge_stencil_at_edge_y, double * _data_vertexFaceDst, double * const _data_vertex_stencil_at_edge_x, double * const _data_vertex_stencil_at_edge_xy, double * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_32 = _data_edge_stencil_at_edge_xy[0];
   const double xi_42 = 1 / (xi_32);
   const double xi_33 = _data_edge_stencil_at_edge_xy[1];
   const double xi_34 = _data_edge_stencil_at_edge_xy[4];
   const double xi_35 = _data_edge_stencil_at_edge_xy[3];
   const double xi_36 = _data_edge_stencil_at_edge_xy[2];
   const double xi_37 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_38 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_39 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_40 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_70 = _data_edge_stencil_at_edge_xy[0];
   const double xi_89 = 1 / (xi_70);
   const double xi_71 = _data_edge_stencil_at_edge_xy[1];
   const double xi_72 = _data_edge_stencil_at_edge_xy[4];
   const double xi_73 = _data_edge_stencil_at_edge_xy[3];
   const double xi_74 = _data_edge_stencil_at_edge_xy[2];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_76 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_77 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_78 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_79 = _data_edge_stencil_at_edge_y[0];
   const double xi_90 = 1 / (xi_79);
   const double xi_80 = _data_edge_stencil_at_edge_y[4];
   const double xi_81 = _data_edge_stencil_at_edge_y[3];
   const double xi_82 = _data_edge_stencil_at_edge_y[1];
   const double xi_83 = _data_edge_stencil_at_edge_y[2];
   const double xi_84 = _data_vertex_stencil_at_edge_y[3];
   const double xi_85 = _data_vertex_stencil_at_edge_y[0];
   const double xi_86 = _data_vertex_stencil_at_edge_y[2];
   const double xi_87 = _data_vertex_stencil_at_edge_y[1];
   const double xi_44 = _data_edge_stencil_at_edge_y[0];
   const double xi_54 = 1 / (xi_44);
   const double xi_45 = _data_edge_stencil_at_edge_y[4];
   const double xi_46 = _data_edge_stencil_at_edge_y[3];
   const double xi_47 = _data_edge_stencil_at_edge_y[1];
   const double xi_48 = _data_edge_stencil_at_edge_y[2];
   const double xi_49 = _data_vertex_stencil_at_edge_y[3];
   const double xi_50 = _data_vertex_stencil_at_edge_y[0];
   const double xi_51 = _data_vertex_stencil_at_edge_y[2];
   const double xi_52 = _data_vertex_stencil_at_edge_y[1];
   const double xi_114 = _data_edge_stencil_at_edge_x[0];
   const double xi_133 = 1 / (xi_114);
   const double xi_115 = _data_edge_stencil_at_edge_x[1];
   const double xi_116 = _data_edge_stencil_at_edge_x[3];
   const double xi_117 = _data_edge_stencil_at_edge_x[4];
   const double xi_118 = _data_edge_stencil_at_edge_x[2];
   const double xi_119 = _data_vertex_stencil_at_edge_x[1];
   const double xi_120 = _data_vertex_stencil_at_edge_x[3];
   const double xi_121 = _data_vertex_stencil_at_edge_x[2];
   const double xi_122 = _data_vertex_stencil_at_edge_x[0];
   const double xi_123 = _data_edge_stencil_at_edge_xy[0];
   const double xi_134 = 1 / (xi_123);
   const double xi_124 = _data_edge_stencil_at_edge_xy[1];
   const double xi_125 = _data_edge_stencil_at_edge_xy[4];
   const double xi_126 = _data_edge_stencil_at_edge_xy[3];
   const double xi_127 = _data_edge_stencil_at_edge_xy[2];
   const double xi_128 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_129 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_130 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_131 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_0 = _data_edge_stencil_at_edge_x[0];
   const double xi_28 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_edge_x[1];
   const double xi_2 = _data_edge_stencil_at_edge_x[3];
   const double xi_3 = _data_edge_stencil_at_edge_x[4];
   const double xi_4 = _data_edge_stencil_at_edge_x[2];
   const double xi_5 = _data_vertex_stencil_at_edge_x[1];
   const double xi_6 = _data_vertex_stencil_at_edge_x[3];
   const double xi_7 = _data_vertex_stencil_at_edge_x[2];
   const double xi_8 = _data_vertex_stencil_at_edge_x[0];
   const double xi_9 = _data_edge_stencil_at_edge_xy[0];
   const double xi_29 = 1 / (xi_9);
   const double xi_10 = _data_edge_stencil_at_edge_xy[1];
   const double xi_11 = _data_edge_stencil_at_edge_xy[4];
   const double xi_12 = _data_edge_stencil_at_edge_xy[3];
   const double xi_13 = _data_edge_stencil_at_edge_xy[2];
   const double xi_14 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_15 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_16 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_17 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_18 = _data_edge_stencil_at_edge_y[0];
   const double xi_30 = 1 / (xi_18);
   const double xi_19 = _data_edge_stencil_at_edge_y[4];
   const double xi_20 = _data_edge_stencil_at_edge_y[3];
   const double xi_21 = _data_edge_stencil_at_edge_y[1];
   const double xi_22 = _data_edge_stencil_at_edge_y[2];
   const double xi_23 = _data_vertex_stencil_at_edge_y[3];
   const double xi_24 = _data_vertex_stencil_at_edge_y[0];
   const double xi_25 = _data_vertex_stencil_at_edge_y[2];
   const double xi_26 = _data_vertex_stencil_at_edge_y[1];
   const double xi_92 = _data_edge_stencil_at_edge_x[0];
   const double xi_111 = 1 / (xi_92);
   const double xi_93 = _data_edge_stencil_at_edge_x[1];
   const double xi_94 = _data_edge_stencil_at_edge_x[3];
   const double xi_95 = _data_edge_stencil_at_edge_x[4];
   const double xi_96 = _data_edge_stencil_at_edge_x[2];
   const double xi_97 = _data_vertex_stencil_at_edge_x[1];
   const double xi_98 = _data_vertex_stencil_at_edge_x[3];
   const double xi_99 = _data_vertex_stencil_at_edge_x[2];
   const double xi_100 = _data_vertex_stencil_at_edge_x[0];
   const double xi_101 = _data_edge_stencil_at_edge_y[0];
   const double xi_112 = 1 / (xi_101);
   const double xi_102 = _data_edge_stencil_at_edge_y[4];
   const double xi_103 = _data_edge_stencil_at_edge_y[3];
   const double xi_104 = _data_edge_stencil_at_edge_y[1];
   const double xi_105 = _data_edge_stencil_at_edge_y[2];
   const double xi_106 = _data_vertex_stencil_at_edge_y[3];
   const double xi_107 = _data_vertex_stencil_at_edge_y[0];
   const double xi_108 = _data_vertex_stencil_at_edge_y[2];
   const double xi_109 = _data_vertex_stencil_at_edge_y[1];
   const double xi_56 = _data_edge_stencil_at_edge_x[0];
   const double xi_66 = 1 / (xi_56);
   const double xi_57 = _data_edge_stencil_at_edge_x[1];
   const double xi_58 = _data_edge_stencil_at_edge_x[3];
   const double xi_59 = _data_edge_stencil_at_edge_x[4];
   const double xi_60 = _data_edge_stencil_at_edge_x[2];
   const double xi_61 = _data_vertex_stencil_at_edge_x[1];
   const double xi_62 = _data_vertex_stencil_at_edge_x[3];
   const double xi_63 = _data_vertex_stencil_at_edge_x[2];
   const double xi_64 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_42*(-xi_33*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_34*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_35*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] - xi_36*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34] - xi_39*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_40*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 35] + _data_edgeFaceRhs[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 31; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_89*(-xi_71*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_72*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_73*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] - xi_74*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_75*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_76*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34] - xi_77*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_78*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 35] + _data_edgeFaceRhs[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_90*(-xi_80*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_81*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32] - xi_82*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_83*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_84*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] - xi_85*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_86*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34] - xi_87*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 31; ctr_1 < 32; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_54*(-xi_45*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_46*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32] - xi_47*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_48*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_49*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] - xi_50*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_51*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34] - xi_52*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 31; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_133*(-xi_115*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] - xi_116*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_117*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_118*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32] - xi_119*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_120*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34] - xi_121*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] - xi_122*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_134*(-xi_124*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_125*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_126*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] - xi_127*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_128*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_129*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34] - xi_130*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_131*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 35] + _data_edgeFaceRhs[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 31; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_28*(-xi_1*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] - xi_2*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32] - xi_5*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34] - xi_7*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] - xi_8*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_29*(-xi_10*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] - xi_13*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_14*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34] - xi_16*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_17*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 35] + _data_edgeFaceRhs[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_30*(-xi_19*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32] - xi_21*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] - xi_24*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34] - xi_26*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 31; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_111*(-xi_100*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_93*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] - xi_94*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_95*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_96*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32] - xi_97*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_98*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34] - xi_99*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] + _data_edgeFaceRhs[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_112*(-xi_102*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_103*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 32] - xi_104*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_105*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_106*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 33] - xi_107*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_108*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34] - xi_109*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 31; ctr_2 < 32; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_66*(-xi_57*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] - xi_58*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + ((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_59*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_60*_data_edgeFaceDst[ctr_1 + 33*ctr_2 + 2*((1056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 32] - xi_61*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_62*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 34] - xi_63*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 33] - xi_64*_data_vertexFaceDst[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 33*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_6(double * _data_edgeFaceDst, double * _data_edgeFaceRhs, double * const _data_edge_stencil_at_edge_x, double * const _data_edge_stencil_at_edge_xy, double * const _data_edge_stencil_at_edge_y, double * _data_vertexFaceDst, double * const _data_vertex_stencil_at_edge_x, double * const _data_vertex_stencil_at_edge_xy, double * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_32 = _data_edge_stencil_at_edge_xy[0];
   const double xi_42 = 1 / (xi_32);
   const double xi_33 = _data_edge_stencil_at_edge_xy[1];
   const double xi_34 = _data_edge_stencil_at_edge_xy[4];
   const double xi_35 = _data_edge_stencil_at_edge_xy[3];
   const double xi_36 = _data_edge_stencil_at_edge_xy[2];
   const double xi_37 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_38 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_39 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_40 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_70 = _data_edge_stencil_at_edge_xy[0];
   const double xi_89 = 1 / (xi_70);
   const double xi_71 = _data_edge_stencil_at_edge_xy[1];
   const double xi_72 = _data_edge_stencil_at_edge_xy[4];
   const double xi_73 = _data_edge_stencil_at_edge_xy[3];
   const double xi_74 = _data_edge_stencil_at_edge_xy[2];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_76 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_77 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_78 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_79 = _data_edge_stencil_at_edge_y[0];
   const double xi_90 = 1 / (xi_79);
   const double xi_80 = _data_edge_stencil_at_edge_y[4];
   const double xi_81 = _data_edge_stencil_at_edge_y[3];
   const double xi_82 = _data_edge_stencil_at_edge_y[1];
   const double xi_83 = _data_edge_stencil_at_edge_y[2];
   const double xi_84 = _data_vertex_stencil_at_edge_y[3];
   const double xi_85 = _data_vertex_stencil_at_edge_y[0];
   const double xi_86 = _data_vertex_stencil_at_edge_y[2];
   const double xi_87 = _data_vertex_stencil_at_edge_y[1];
   const double xi_44 = _data_edge_stencil_at_edge_y[0];
   const double xi_54 = 1 / (xi_44);
   const double xi_45 = _data_edge_stencil_at_edge_y[4];
   const double xi_46 = _data_edge_stencil_at_edge_y[3];
   const double xi_47 = _data_edge_stencil_at_edge_y[1];
   const double xi_48 = _data_edge_stencil_at_edge_y[2];
   const double xi_49 = _data_vertex_stencil_at_edge_y[3];
   const double xi_50 = _data_vertex_stencil_at_edge_y[0];
   const double xi_51 = _data_vertex_stencil_at_edge_y[2];
   const double xi_52 = _data_vertex_stencil_at_edge_y[1];
   const double xi_114 = _data_edge_stencil_at_edge_x[0];
   const double xi_133 = 1 / (xi_114);
   const double xi_115 = _data_edge_stencil_at_edge_x[1];
   const double xi_116 = _data_edge_stencil_at_edge_x[3];
   const double xi_117 = _data_edge_stencil_at_edge_x[4];
   const double xi_118 = _data_edge_stencil_at_edge_x[2];
   const double xi_119 = _data_vertex_stencil_at_edge_x[1];
   const double xi_120 = _data_vertex_stencil_at_edge_x[3];
   const double xi_121 = _data_vertex_stencil_at_edge_x[2];
   const double xi_122 = _data_vertex_stencil_at_edge_x[0];
   const double xi_123 = _data_edge_stencil_at_edge_xy[0];
   const double xi_134 = 1 / (xi_123);
   const double xi_124 = _data_edge_stencil_at_edge_xy[1];
   const double xi_125 = _data_edge_stencil_at_edge_xy[4];
   const double xi_126 = _data_edge_stencil_at_edge_xy[3];
   const double xi_127 = _data_edge_stencil_at_edge_xy[2];
   const double xi_128 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_129 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_130 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_131 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_0 = _data_edge_stencil_at_edge_x[0];
   const double xi_28 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_edge_x[1];
   const double xi_2 = _data_edge_stencil_at_edge_x[3];
   const double xi_3 = _data_edge_stencil_at_edge_x[4];
   const double xi_4 = _data_edge_stencil_at_edge_x[2];
   const double xi_5 = _data_vertex_stencil_at_edge_x[1];
   const double xi_6 = _data_vertex_stencil_at_edge_x[3];
   const double xi_7 = _data_vertex_stencil_at_edge_x[2];
   const double xi_8 = _data_vertex_stencil_at_edge_x[0];
   const double xi_9 = _data_edge_stencil_at_edge_xy[0];
   const double xi_29 = 1 / (xi_9);
   const double xi_10 = _data_edge_stencil_at_edge_xy[1];
   const double xi_11 = _data_edge_stencil_at_edge_xy[4];
   const double xi_12 = _data_edge_stencil_at_edge_xy[3];
   const double xi_13 = _data_edge_stencil_at_edge_xy[2];
   const double xi_14 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_15 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_16 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_17 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_18 = _data_edge_stencil_at_edge_y[0];
   const double xi_30 = 1 / (xi_18);
   const double xi_19 = _data_edge_stencil_at_edge_y[4];
   const double xi_20 = _data_edge_stencil_at_edge_y[3];
   const double xi_21 = _data_edge_stencil_at_edge_y[1];
   const double xi_22 = _data_edge_stencil_at_edge_y[2];
   const double xi_23 = _data_vertex_stencil_at_edge_y[3];
   const double xi_24 = _data_vertex_stencil_at_edge_y[0];
   const double xi_25 = _data_vertex_stencil_at_edge_y[2];
   const double xi_26 = _data_vertex_stencil_at_edge_y[1];
   const double xi_92 = _data_edge_stencil_at_edge_x[0];
   const double xi_111 = 1 / (xi_92);
   const double xi_93 = _data_edge_stencil_at_edge_x[1];
   const double xi_94 = _data_edge_stencil_at_edge_x[3];
   const double xi_95 = _data_edge_stencil_at_edge_x[4];
   const double xi_96 = _data_edge_stencil_at_edge_x[2];
   const double xi_97 = _data_vertex_stencil_at_edge_x[1];
   const double xi_98 = _data_vertex_stencil_at_edge_x[3];
   const double xi_99 = _data_vertex_stencil_at_edge_x[2];
   const double xi_100 = _data_vertex_stencil_at_edge_x[0];
   const double xi_101 = _data_edge_stencil_at_edge_y[0];
   const double xi_112 = 1 / (xi_101);
   const double xi_102 = _data_edge_stencil_at_edge_y[4];
   const double xi_103 = _data_edge_stencil_at_edge_y[3];
   const double xi_104 = _data_edge_stencil_at_edge_y[1];
   const double xi_105 = _data_edge_stencil_at_edge_y[2];
   const double xi_106 = _data_vertex_stencil_at_edge_y[3];
   const double xi_107 = _data_vertex_stencil_at_edge_y[0];
   const double xi_108 = _data_vertex_stencil_at_edge_y[2];
   const double xi_109 = _data_vertex_stencil_at_edge_y[1];
   const double xi_56 = _data_edge_stencil_at_edge_x[0];
   const double xi_66 = 1 / (xi_56);
   const double xi_57 = _data_edge_stencil_at_edge_x[1];
   const double xi_58 = _data_edge_stencil_at_edge_x[3];
   const double xi_59 = _data_edge_stencil_at_edge_x[4];
   const double xi_60 = _data_edge_stencil_at_edge_x[2];
   const double xi_61 = _data_vertex_stencil_at_edge_x[1];
   const double xi_62 = _data_vertex_stencil_at_edge_x[3];
   const double xi_63 = _data_vertex_stencil_at_edge_x[2];
   const double xi_64 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_42*(-xi_33*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_34*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_35*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] - xi_36*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66] - xi_39*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_40*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 67] + _data_edgeFaceRhs[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 63; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_89*(-xi_71*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_72*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_73*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] - xi_74*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_75*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_76*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66] - xi_77*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_78*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 67] + _data_edgeFaceRhs[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_90*(-xi_80*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_81*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64] - xi_82*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_83*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_84*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] - xi_85*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_86*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66] - xi_87*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 63; ctr_1 < 64; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_54*(-xi_45*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_46*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64] - xi_47*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_48*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_49*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] - xi_50*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_51*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66] - xi_52*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 63; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_133*(-xi_115*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] - xi_116*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_117*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_118*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64] - xi_119*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_120*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66] - xi_121*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] - xi_122*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_134*(-xi_124*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_125*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_126*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] - xi_127*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_128*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_129*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66] - xi_130*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_131*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 67] + _data_edgeFaceRhs[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 63; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_28*(-xi_1*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] - xi_2*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64] - xi_5*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66] - xi_7*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] - xi_8*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_29*(-xi_10*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] - xi_13*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_14*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66] - xi_16*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_17*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 67] + _data_edgeFaceRhs[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_30*(-xi_19*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64] - xi_21*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] - xi_24*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66] - xi_26*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 63; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_111*(-xi_100*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_93*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] - xi_94*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_95*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_96*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64] - xi_97*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_98*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66] - xi_99*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] + _data_edgeFaceRhs[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_112*(-xi_102*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_103*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 64] - xi_104*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_105*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_106*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 65] - xi_107*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_108*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66] - xi_109*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 63; ctr_2 < 64; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_66*(-xi_57*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] - xi_58*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + ((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_59*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_60*_data_edgeFaceDst[ctr_1 + 65*ctr_2 + 2*((4160) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 64] - xi_61*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_62*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 66] - xi_63*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 65] - xi_64*_data_vertexFaceDst[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 65*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_7(double * _data_edgeFaceDst, double * _data_edgeFaceRhs, double * const _data_edge_stencil_at_edge_x, double * const _data_edge_stencil_at_edge_xy, double * const _data_edge_stencil_at_edge_y, double * _data_vertexFaceDst, double * const _data_vertex_stencil_at_edge_x, double * const _data_vertex_stencil_at_edge_xy, double * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_32 = _data_edge_stencil_at_edge_xy[0];
   const double xi_42 = 1 / (xi_32);
   const double xi_33 = _data_edge_stencil_at_edge_xy[1];
   const double xi_34 = _data_edge_stencil_at_edge_xy[4];
   const double xi_35 = _data_edge_stencil_at_edge_xy[3];
   const double xi_36 = _data_edge_stencil_at_edge_xy[2];
   const double xi_37 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_38 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_39 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_40 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_70 = _data_edge_stencil_at_edge_xy[0];
   const double xi_89 = 1 / (xi_70);
   const double xi_71 = _data_edge_stencil_at_edge_xy[1];
   const double xi_72 = _data_edge_stencil_at_edge_xy[4];
   const double xi_73 = _data_edge_stencil_at_edge_xy[3];
   const double xi_74 = _data_edge_stencil_at_edge_xy[2];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_76 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_77 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_78 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_79 = _data_edge_stencil_at_edge_y[0];
   const double xi_90 = 1 / (xi_79);
   const double xi_80 = _data_edge_stencil_at_edge_y[4];
   const double xi_81 = _data_edge_stencil_at_edge_y[3];
   const double xi_82 = _data_edge_stencil_at_edge_y[1];
   const double xi_83 = _data_edge_stencil_at_edge_y[2];
   const double xi_84 = _data_vertex_stencil_at_edge_y[3];
   const double xi_85 = _data_vertex_stencil_at_edge_y[0];
   const double xi_86 = _data_vertex_stencil_at_edge_y[2];
   const double xi_87 = _data_vertex_stencil_at_edge_y[1];
   const double xi_44 = _data_edge_stencil_at_edge_y[0];
   const double xi_54 = 1 / (xi_44);
   const double xi_45 = _data_edge_stencil_at_edge_y[4];
   const double xi_46 = _data_edge_stencil_at_edge_y[3];
   const double xi_47 = _data_edge_stencil_at_edge_y[1];
   const double xi_48 = _data_edge_stencil_at_edge_y[2];
   const double xi_49 = _data_vertex_stencil_at_edge_y[3];
   const double xi_50 = _data_vertex_stencil_at_edge_y[0];
   const double xi_51 = _data_vertex_stencil_at_edge_y[2];
   const double xi_52 = _data_vertex_stencil_at_edge_y[1];
   const double xi_114 = _data_edge_stencil_at_edge_x[0];
   const double xi_133 = 1 / (xi_114);
   const double xi_115 = _data_edge_stencil_at_edge_x[1];
   const double xi_116 = _data_edge_stencil_at_edge_x[3];
   const double xi_117 = _data_edge_stencil_at_edge_x[4];
   const double xi_118 = _data_edge_stencil_at_edge_x[2];
   const double xi_119 = _data_vertex_stencil_at_edge_x[1];
   const double xi_120 = _data_vertex_stencil_at_edge_x[3];
   const double xi_121 = _data_vertex_stencil_at_edge_x[2];
   const double xi_122 = _data_vertex_stencil_at_edge_x[0];
   const double xi_123 = _data_edge_stencil_at_edge_xy[0];
   const double xi_134 = 1 / (xi_123);
   const double xi_124 = _data_edge_stencil_at_edge_xy[1];
   const double xi_125 = _data_edge_stencil_at_edge_xy[4];
   const double xi_126 = _data_edge_stencil_at_edge_xy[3];
   const double xi_127 = _data_edge_stencil_at_edge_xy[2];
   const double xi_128 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_129 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_130 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_131 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_0 = _data_edge_stencil_at_edge_x[0];
   const double xi_28 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_edge_x[1];
   const double xi_2 = _data_edge_stencil_at_edge_x[3];
   const double xi_3 = _data_edge_stencil_at_edge_x[4];
   const double xi_4 = _data_edge_stencil_at_edge_x[2];
   const double xi_5 = _data_vertex_stencil_at_edge_x[1];
   const double xi_6 = _data_vertex_stencil_at_edge_x[3];
   const double xi_7 = _data_vertex_stencil_at_edge_x[2];
   const double xi_8 = _data_vertex_stencil_at_edge_x[0];
   const double xi_9 = _data_edge_stencil_at_edge_xy[0];
   const double xi_29 = 1 / (xi_9);
   const double xi_10 = _data_edge_stencil_at_edge_xy[1];
   const double xi_11 = _data_edge_stencil_at_edge_xy[4];
   const double xi_12 = _data_edge_stencil_at_edge_xy[3];
   const double xi_13 = _data_edge_stencil_at_edge_xy[2];
   const double xi_14 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_15 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_16 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_17 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_18 = _data_edge_stencil_at_edge_y[0];
   const double xi_30 = 1 / (xi_18);
   const double xi_19 = _data_edge_stencil_at_edge_y[4];
   const double xi_20 = _data_edge_stencil_at_edge_y[3];
   const double xi_21 = _data_edge_stencil_at_edge_y[1];
   const double xi_22 = _data_edge_stencil_at_edge_y[2];
   const double xi_23 = _data_vertex_stencil_at_edge_y[3];
   const double xi_24 = _data_vertex_stencil_at_edge_y[0];
   const double xi_25 = _data_vertex_stencil_at_edge_y[2];
   const double xi_26 = _data_vertex_stencil_at_edge_y[1];
   const double xi_92 = _data_edge_stencil_at_edge_x[0];
   const double xi_111 = 1 / (xi_92);
   const double xi_93 = _data_edge_stencil_at_edge_x[1];
   const double xi_94 = _data_edge_stencil_at_edge_x[3];
   const double xi_95 = _data_edge_stencil_at_edge_x[4];
   const double xi_96 = _data_edge_stencil_at_edge_x[2];
   const double xi_97 = _data_vertex_stencil_at_edge_x[1];
   const double xi_98 = _data_vertex_stencil_at_edge_x[3];
   const double xi_99 = _data_vertex_stencil_at_edge_x[2];
   const double xi_100 = _data_vertex_stencil_at_edge_x[0];
   const double xi_101 = _data_edge_stencil_at_edge_y[0];
   const double xi_112 = 1 / (xi_101);
   const double xi_102 = _data_edge_stencil_at_edge_y[4];
   const double xi_103 = _data_edge_stencil_at_edge_y[3];
   const double xi_104 = _data_edge_stencil_at_edge_y[1];
   const double xi_105 = _data_edge_stencil_at_edge_y[2];
   const double xi_106 = _data_vertex_stencil_at_edge_y[3];
   const double xi_107 = _data_vertex_stencil_at_edge_y[0];
   const double xi_108 = _data_vertex_stencil_at_edge_y[2];
   const double xi_109 = _data_vertex_stencil_at_edge_y[1];
   const double xi_56 = _data_edge_stencil_at_edge_x[0];
   const double xi_66 = 1 / (xi_56);
   const double xi_57 = _data_edge_stencil_at_edge_x[1];
   const double xi_58 = _data_edge_stencil_at_edge_x[3];
   const double xi_59 = _data_edge_stencil_at_edge_x[4];
   const double xi_60 = _data_edge_stencil_at_edge_x[2];
   const double xi_61 = _data_vertex_stencil_at_edge_x[1];
   const double xi_62 = _data_vertex_stencil_at_edge_x[3];
   const double xi_63 = _data_vertex_stencil_at_edge_x[2];
   const double xi_64 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_42*(-xi_33*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_34*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_35*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] - xi_36*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130] - xi_39*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_40*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 131] + _data_edgeFaceRhs[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 127; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_89*(-xi_71*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_72*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_73*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] - xi_74*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_75*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_76*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130] - xi_77*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_78*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 131] + _data_edgeFaceRhs[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_90*(-xi_80*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_81*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128] - xi_82*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_83*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_84*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] - xi_85*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_86*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130] - xi_87*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 127; ctr_1 < 128; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_54*(-xi_45*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_46*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128] - xi_47*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_48*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_49*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] - xi_50*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_51*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130] - xi_52*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 127; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_133*(-xi_115*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] - xi_116*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_117*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_118*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128] - xi_119*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_120*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130] - xi_121*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] - xi_122*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_134*(-xi_124*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_125*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_126*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] - xi_127*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_128*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_129*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130] - xi_130*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_131*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 131] + _data_edgeFaceRhs[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 127; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_28*(-xi_1*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] - xi_2*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128] - xi_5*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130] - xi_7*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] - xi_8*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_29*(-xi_10*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] - xi_13*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_14*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130] - xi_16*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_17*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 131] + _data_edgeFaceRhs[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_30*(-xi_19*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128] - xi_21*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] - xi_24*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130] - xi_26*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 127; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_111*(-xi_100*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_93*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] - xi_94*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_95*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_96*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128] - xi_97*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_98*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130] - xi_99*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] + _data_edgeFaceRhs[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_112*(-xi_102*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_103*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 128] - xi_104*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_105*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_106*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 129] - xi_107*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_108*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130] - xi_109*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 127; ctr_2 < 128; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_66*(-xi_57*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] - xi_58*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + ((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_59*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_60*_data_edgeFaceDst[ctr_1 + 129*ctr_2 + 2*((16512) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 128] - xi_61*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_62*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 130] - xi_63*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 129] - xi_64*_data_vertexFaceDst[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 129*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_8(double * _data_edgeFaceDst, double * _data_edgeFaceRhs, double * const _data_edge_stencil_at_edge_x, double * const _data_edge_stencil_at_edge_xy, double * const _data_edge_stencil_at_edge_y, double * _data_vertexFaceDst, double * const _data_vertex_stencil_at_edge_x, double * const _data_vertex_stencil_at_edge_xy, double * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_32 = _data_edge_stencil_at_edge_xy[0];
   const double xi_42 = 1 / (xi_32);
   const double xi_33 = _data_edge_stencil_at_edge_xy[1];
   const double xi_34 = _data_edge_stencil_at_edge_xy[4];
   const double xi_35 = _data_edge_stencil_at_edge_xy[3];
   const double xi_36 = _data_edge_stencil_at_edge_xy[2];
   const double xi_37 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_38 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_39 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_40 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_70 = _data_edge_stencil_at_edge_xy[0];
   const double xi_89 = 1 / (xi_70);
   const double xi_71 = _data_edge_stencil_at_edge_xy[1];
   const double xi_72 = _data_edge_stencil_at_edge_xy[4];
   const double xi_73 = _data_edge_stencil_at_edge_xy[3];
   const double xi_74 = _data_edge_stencil_at_edge_xy[2];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_76 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_77 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_78 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_79 = _data_edge_stencil_at_edge_y[0];
   const double xi_90 = 1 / (xi_79);
   const double xi_80 = _data_edge_stencil_at_edge_y[4];
   const double xi_81 = _data_edge_stencil_at_edge_y[3];
   const double xi_82 = _data_edge_stencil_at_edge_y[1];
   const double xi_83 = _data_edge_stencil_at_edge_y[2];
   const double xi_84 = _data_vertex_stencil_at_edge_y[3];
   const double xi_85 = _data_vertex_stencil_at_edge_y[0];
   const double xi_86 = _data_vertex_stencil_at_edge_y[2];
   const double xi_87 = _data_vertex_stencil_at_edge_y[1];
   const double xi_44 = _data_edge_stencil_at_edge_y[0];
   const double xi_54 = 1 / (xi_44);
   const double xi_45 = _data_edge_stencil_at_edge_y[4];
   const double xi_46 = _data_edge_stencil_at_edge_y[3];
   const double xi_47 = _data_edge_stencil_at_edge_y[1];
   const double xi_48 = _data_edge_stencil_at_edge_y[2];
   const double xi_49 = _data_vertex_stencil_at_edge_y[3];
   const double xi_50 = _data_vertex_stencil_at_edge_y[0];
   const double xi_51 = _data_vertex_stencil_at_edge_y[2];
   const double xi_52 = _data_vertex_stencil_at_edge_y[1];
   const double xi_114 = _data_edge_stencil_at_edge_x[0];
   const double xi_133 = 1 / (xi_114);
   const double xi_115 = _data_edge_stencil_at_edge_x[1];
   const double xi_116 = _data_edge_stencil_at_edge_x[3];
   const double xi_117 = _data_edge_stencil_at_edge_x[4];
   const double xi_118 = _data_edge_stencil_at_edge_x[2];
   const double xi_119 = _data_vertex_stencil_at_edge_x[1];
   const double xi_120 = _data_vertex_stencil_at_edge_x[3];
   const double xi_121 = _data_vertex_stencil_at_edge_x[2];
   const double xi_122 = _data_vertex_stencil_at_edge_x[0];
   const double xi_123 = _data_edge_stencil_at_edge_xy[0];
   const double xi_134 = 1 / (xi_123);
   const double xi_124 = _data_edge_stencil_at_edge_xy[1];
   const double xi_125 = _data_edge_stencil_at_edge_xy[4];
   const double xi_126 = _data_edge_stencil_at_edge_xy[3];
   const double xi_127 = _data_edge_stencil_at_edge_xy[2];
   const double xi_128 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_129 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_130 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_131 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_0 = _data_edge_stencil_at_edge_x[0];
   const double xi_28 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_edge_x[1];
   const double xi_2 = _data_edge_stencil_at_edge_x[3];
   const double xi_3 = _data_edge_stencil_at_edge_x[4];
   const double xi_4 = _data_edge_stencil_at_edge_x[2];
   const double xi_5 = _data_vertex_stencil_at_edge_x[1];
   const double xi_6 = _data_vertex_stencil_at_edge_x[3];
   const double xi_7 = _data_vertex_stencil_at_edge_x[2];
   const double xi_8 = _data_vertex_stencil_at_edge_x[0];
   const double xi_9 = _data_edge_stencil_at_edge_xy[0];
   const double xi_29 = 1 / (xi_9);
   const double xi_10 = _data_edge_stencil_at_edge_xy[1];
   const double xi_11 = _data_edge_stencil_at_edge_xy[4];
   const double xi_12 = _data_edge_stencil_at_edge_xy[3];
   const double xi_13 = _data_edge_stencil_at_edge_xy[2];
   const double xi_14 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_15 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_16 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_17 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_18 = _data_edge_stencil_at_edge_y[0];
   const double xi_30 = 1 / (xi_18);
   const double xi_19 = _data_edge_stencil_at_edge_y[4];
   const double xi_20 = _data_edge_stencil_at_edge_y[3];
   const double xi_21 = _data_edge_stencil_at_edge_y[1];
   const double xi_22 = _data_edge_stencil_at_edge_y[2];
   const double xi_23 = _data_vertex_stencil_at_edge_y[3];
   const double xi_24 = _data_vertex_stencil_at_edge_y[0];
   const double xi_25 = _data_vertex_stencil_at_edge_y[2];
   const double xi_26 = _data_vertex_stencil_at_edge_y[1];
   const double xi_92 = _data_edge_stencil_at_edge_x[0];
   const double xi_111 = 1 / (xi_92);
   const double xi_93 = _data_edge_stencil_at_edge_x[1];
   const double xi_94 = _data_edge_stencil_at_edge_x[3];
   const double xi_95 = _data_edge_stencil_at_edge_x[4];
   const double xi_96 = _data_edge_stencil_at_edge_x[2];
   const double xi_97 = _data_vertex_stencil_at_edge_x[1];
   const double xi_98 = _data_vertex_stencil_at_edge_x[3];
   const double xi_99 = _data_vertex_stencil_at_edge_x[2];
   const double xi_100 = _data_vertex_stencil_at_edge_x[0];
   const double xi_101 = _data_edge_stencil_at_edge_y[0];
   const double xi_112 = 1 / (xi_101);
   const double xi_102 = _data_edge_stencil_at_edge_y[4];
   const double xi_103 = _data_edge_stencil_at_edge_y[3];
   const double xi_104 = _data_edge_stencil_at_edge_y[1];
   const double xi_105 = _data_edge_stencil_at_edge_y[2];
   const double xi_106 = _data_vertex_stencil_at_edge_y[3];
   const double xi_107 = _data_vertex_stencil_at_edge_y[0];
   const double xi_108 = _data_vertex_stencil_at_edge_y[2];
   const double xi_109 = _data_vertex_stencil_at_edge_y[1];
   const double xi_56 = _data_edge_stencil_at_edge_x[0];
   const double xi_66 = 1 / (xi_56);
   const double xi_57 = _data_edge_stencil_at_edge_x[1];
   const double xi_58 = _data_edge_stencil_at_edge_x[3];
   const double xi_59 = _data_edge_stencil_at_edge_x[4];
   const double xi_60 = _data_edge_stencil_at_edge_x[2];
   const double xi_61 = _data_vertex_stencil_at_edge_x[1];
   const double xi_62 = _data_vertex_stencil_at_edge_x[3];
   const double xi_63 = _data_vertex_stencil_at_edge_x[2];
   const double xi_64 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_42*(-xi_33*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_34*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_35*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] - xi_36*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258] - xi_39*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_40*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 259] + _data_edgeFaceRhs[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 255; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_89*(-xi_71*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_72*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_73*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] - xi_74*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_75*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_76*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258] - xi_77*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_78*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 259] + _data_edgeFaceRhs[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_90*(-xi_80*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_81*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256] - xi_82*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_83*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_84*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] - xi_85*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_86*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258] - xi_87*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 255; ctr_1 < 256; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_54*(-xi_45*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_46*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256] - xi_47*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_48*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_49*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] - xi_50*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_51*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258] - xi_52*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 255; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_133*(-xi_115*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] - xi_116*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_117*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_118*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256] - xi_119*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_120*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258] - xi_121*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] - xi_122*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_134*(-xi_124*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_125*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_126*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] - xi_127*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_128*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_129*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258] - xi_130*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_131*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 259] + _data_edgeFaceRhs[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 255; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_28*(-xi_1*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] - xi_2*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256] - xi_5*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258] - xi_7*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] - xi_8*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_29*(-xi_10*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] - xi_13*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_14*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258] - xi_16*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_17*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 259] + _data_edgeFaceRhs[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_30*(-xi_19*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256] - xi_21*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] - xi_24*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258] - xi_26*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 255; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_111*(-xi_100*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_93*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] - xi_94*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_95*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_96*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256] - xi_97*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_98*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258] - xi_99*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] + _data_edgeFaceRhs[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_112*(-xi_102*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_103*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 256] - xi_104*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_105*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_106*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 257] - xi_107*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_108*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258] - xi_109*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 255; ctr_2 < 256; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_66*(-xi_57*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] - xi_58*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + ((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_59*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_60*_data_edgeFaceDst[ctr_1 + 257*ctr_2 + 2*((65792) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 256] - xi_61*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_62*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 258] - xi_63*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 257] - xi_64*_data_vertexFaceDst[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 257*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_9(double * _data_edgeFaceDst, double * _data_edgeFaceRhs, double * const _data_edge_stencil_at_edge_x, double * const _data_edge_stencil_at_edge_xy, double * const _data_edge_stencil_at_edge_y, double * _data_vertexFaceDst, double * const _data_vertex_stencil_at_edge_x, double * const _data_vertex_stencil_at_edge_xy, double * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_32 = _data_edge_stencil_at_edge_xy[0];
   const double xi_42 = 1 / (xi_32);
   const double xi_33 = _data_edge_stencil_at_edge_xy[1];
   const double xi_34 = _data_edge_stencil_at_edge_xy[4];
   const double xi_35 = _data_edge_stencil_at_edge_xy[3];
   const double xi_36 = _data_edge_stencil_at_edge_xy[2];
   const double xi_37 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_38 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_39 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_40 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_70 = _data_edge_stencil_at_edge_xy[0];
   const double xi_89 = 1 / (xi_70);
   const double xi_71 = _data_edge_stencil_at_edge_xy[1];
   const double xi_72 = _data_edge_stencil_at_edge_xy[4];
   const double xi_73 = _data_edge_stencil_at_edge_xy[3];
   const double xi_74 = _data_edge_stencil_at_edge_xy[2];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_76 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_77 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_78 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_79 = _data_edge_stencil_at_edge_y[0];
   const double xi_90 = 1 / (xi_79);
   const double xi_80 = _data_edge_stencil_at_edge_y[4];
   const double xi_81 = _data_edge_stencil_at_edge_y[3];
   const double xi_82 = _data_edge_stencil_at_edge_y[1];
   const double xi_83 = _data_edge_stencil_at_edge_y[2];
   const double xi_84 = _data_vertex_stencil_at_edge_y[3];
   const double xi_85 = _data_vertex_stencil_at_edge_y[0];
   const double xi_86 = _data_vertex_stencil_at_edge_y[2];
   const double xi_87 = _data_vertex_stencil_at_edge_y[1];
   const double xi_44 = _data_edge_stencil_at_edge_y[0];
   const double xi_54 = 1 / (xi_44);
   const double xi_45 = _data_edge_stencil_at_edge_y[4];
   const double xi_46 = _data_edge_stencil_at_edge_y[3];
   const double xi_47 = _data_edge_stencil_at_edge_y[1];
   const double xi_48 = _data_edge_stencil_at_edge_y[2];
   const double xi_49 = _data_vertex_stencil_at_edge_y[3];
   const double xi_50 = _data_vertex_stencil_at_edge_y[0];
   const double xi_51 = _data_vertex_stencil_at_edge_y[2];
   const double xi_52 = _data_vertex_stencil_at_edge_y[1];
   const double xi_114 = _data_edge_stencil_at_edge_x[0];
   const double xi_133 = 1 / (xi_114);
   const double xi_115 = _data_edge_stencil_at_edge_x[1];
   const double xi_116 = _data_edge_stencil_at_edge_x[3];
   const double xi_117 = _data_edge_stencil_at_edge_x[4];
   const double xi_118 = _data_edge_stencil_at_edge_x[2];
   const double xi_119 = _data_vertex_stencil_at_edge_x[1];
   const double xi_120 = _data_vertex_stencil_at_edge_x[3];
   const double xi_121 = _data_vertex_stencil_at_edge_x[2];
   const double xi_122 = _data_vertex_stencil_at_edge_x[0];
   const double xi_123 = _data_edge_stencil_at_edge_xy[0];
   const double xi_134 = 1 / (xi_123);
   const double xi_124 = _data_edge_stencil_at_edge_xy[1];
   const double xi_125 = _data_edge_stencil_at_edge_xy[4];
   const double xi_126 = _data_edge_stencil_at_edge_xy[3];
   const double xi_127 = _data_edge_stencil_at_edge_xy[2];
   const double xi_128 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_129 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_130 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_131 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_0 = _data_edge_stencil_at_edge_x[0];
   const double xi_28 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_edge_x[1];
   const double xi_2 = _data_edge_stencil_at_edge_x[3];
   const double xi_3 = _data_edge_stencil_at_edge_x[4];
   const double xi_4 = _data_edge_stencil_at_edge_x[2];
   const double xi_5 = _data_vertex_stencil_at_edge_x[1];
   const double xi_6 = _data_vertex_stencil_at_edge_x[3];
   const double xi_7 = _data_vertex_stencil_at_edge_x[2];
   const double xi_8 = _data_vertex_stencil_at_edge_x[0];
   const double xi_9 = _data_edge_stencil_at_edge_xy[0];
   const double xi_29 = 1 / (xi_9);
   const double xi_10 = _data_edge_stencil_at_edge_xy[1];
   const double xi_11 = _data_edge_stencil_at_edge_xy[4];
   const double xi_12 = _data_edge_stencil_at_edge_xy[3];
   const double xi_13 = _data_edge_stencil_at_edge_xy[2];
   const double xi_14 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_15 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_16 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_17 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_18 = _data_edge_stencil_at_edge_y[0];
   const double xi_30 = 1 / (xi_18);
   const double xi_19 = _data_edge_stencil_at_edge_y[4];
   const double xi_20 = _data_edge_stencil_at_edge_y[3];
   const double xi_21 = _data_edge_stencil_at_edge_y[1];
   const double xi_22 = _data_edge_stencil_at_edge_y[2];
   const double xi_23 = _data_vertex_stencil_at_edge_y[3];
   const double xi_24 = _data_vertex_stencil_at_edge_y[0];
   const double xi_25 = _data_vertex_stencil_at_edge_y[2];
   const double xi_26 = _data_vertex_stencil_at_edge_y[1];
   const double xi_92 = _data_edge_stencil_at_edge_x[0];
   const double xi_111 = 1 / (xi_92);
   const double xi_93 = _data_edge_stencil_at_edge_x[1];
   const double xi_94 = _data_edge_stencil_at_edge_x[3];
   const double xi_95 = _data_edge_stencil_at_edge_x[4];
   const double xi_96 = _data_edge_stencil_at_edge_x[2];
   const double xi_97 = _data_vertex_stencil_at_edge_x[1];
   const double xi_98 = _data_vertex_stencil_at_edge_x[3];
   const double xi_99 = _data_vertex_stencil_at_edge_x[2];
   const double xi_100 = _data_vertex_stencil_at_edge_x[0];
   const double xi_101 = _data_edge_stencil_at_edge_y[0];
   const double xi_112 = 1 / (xi_101);
   const double xi_102 = _data_edge_stencil_at_edge_y[4];
   const double xi_103 = _data_edge_stencil_at_edge_y[3];
   const double xi_104 = _data_edge_stencil_at_edge_y[1];
   const double xi_105 = _data_edge_stencil_at_edge_y[2];
   const double xi_106 = _data_vertex_stencil_at_edge_y[3];
   const double xi_107 = _data_vertex_stencil_at_edge_y[0];
   const double xi_108 = _data_vertex_stencil_at_edge_y[2];
   const double xi_109 = _data_vertex_stencil_at_edge_y[1];
   const double xi_56 = _data_edge_stencil_at_edge_x[0];
   const double xi_66 = 1 / (xi_56);
   const double xi_57 = _data_edge_stencil_at_edge_x[1];
   const double xi_58 = _data_edge_stencil_at_edge_x[3];
   const double xi_59 = _data_edge_stencil_at_edge_x[4];
   const double xi_60 = _data_edge_stencil_at_edge_x[2];
   const double xi_61 = _data_vertex_stencil_at_edge_x[1];
   const double xi_62 = _data_vertex_stencil_at_edge_x[3];
   const double xi_63 = _data_vertex_stencil_at_edge_x[2];
   const double xi_64 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_42*(-xi_33*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_34*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_35*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] - xi_36*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514] - xi_39*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_40*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 515] + _data_edgeFaceRhs[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 511; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_89*(-xi_71*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_72*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_73*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] - xi_74*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_75*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_76*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514] - xi_77*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_78*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 515] + _data_edgeFaceRhs[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_90*(-xi_80*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_81*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512] - xi_82*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_83*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_84*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] - xi_85*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_86*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514] - xi_87*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 511; ctr_1 < 512; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_54*(-xi_45*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_46*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512] - xi_47*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_48*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_49*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] - xi_50*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_51*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514] - xi_52*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 511; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_133*(-xi_115*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] - xi_116*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_117*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_118*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512] - xi_119*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_120*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514] - xi_121*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] - xi_122*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_134*(-xi_124*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_125*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_126*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] - xi_127*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_128*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_129*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514] - xi_130*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_131*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 515] + _data_edgeFaceRhs[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 511; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_28*(-xi_1*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] - xi_2*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512] - xi_5*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514] - xi_7*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] - xi_8*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_29*(-xi_10*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] - xi_13*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_14*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514] - xi_16*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_17*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 515] + _data_edgeFaceRhs[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_30*(-xi_19*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512] - xi_21*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] - xi_24*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514] - xi_26*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 511; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_111*(-xi_100*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_93*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] - xi_94*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_95*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_96*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512] - xi_97*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_98*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514] - xi_99*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] + _data_edgeFaceRhs[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_112*(-xi_102*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_103*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 512] - xi_104*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_105*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_106*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 513] - xi_107*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_108*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514] - xi_109*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 511; ctr_2 < 512; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_66*(-xi_57*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] - xi_58*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + ((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_59*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_60*_data_edgeFaceDst[ctr_1 + 513*ctr_2 + 2*((262656) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 512] - xi_61*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_62*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 514] - xi_63*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 513] - xi_64*_data_vertexFaceDst[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 513*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_10(double * _data_edgeFaceDst, double * _data_edgeFaceRhs, double * const _data_edge_stencil_at_edge_x, double * const _data_edge_stencil_at_edge_xy, double * const _data_edge_stencil_at_edge_y, double * _data_vertexFaceDst, double * const _data_vertex_stencil_at_edge_x, double * const _data_vertex_stencil_at_edge_xy, double * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_32 = _data_edge_stencil_at_edge_xy[0];
   const double xi_42 = 1 / (xi_32);
   const double xi_33 = _data_edge_stencil_at_edge_xy[1];
   const double xi_34 = _data_edge_stencil_at_edge_xy[4];
   const double xi_35 = _data_edge_stencil_at_edge_xy[3];
   const double xi_36 = _data_edge_stencil_at_edge_xy[2];
   const double xi_37 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_38 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_39 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_40 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_70 = _data_edge_stencil_at_edge_xy[0];
   const double xi_89 = 1 / (xi_70);
   const double xi_71 = _data_edge_stencil_at_edge_xy[1];
   const double xi_72 = _data_edge_stencil_at_edge_xy[4];
   const double xi_73 = _data_edge_stencil_at_edge_xy[3];
   const double xi_74 = _data_edge_stencil_at_edge_xy[2];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_76 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_77 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_78 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_79 = _data_edge_stencil_at_edge_y[0];
   const double xi_90 = 1 / (xi_79);
   const double xi_80 = _data_edge_stencil_at_edge_y[4];
   const double xi_81 = _data_edge_stencil_at_edge_y[3];
   const double xi_82 = _data_edge_stencil_at_edge_y[1];
   const double xi_83 = _data_edge_stencil_at_edge_y[2];
   const double xi_84 = _data_vertex_stencil_at_edge_y[3];
   const double xi_85 = _data_vertex_stencil_at_edge_y[0];
   const double xi_86 = _data_vertex_stencil_at_edge_y[2];
   const double xi_87 = _data_vertex_stencil_at_edge_y[1];
   const double xi_44 = _data_edge_stencil_at_edge_y[0];
   const double xi_54 = 1 / (xi_44);
   const double xi_45 = _data_edge_stencil_at_edge_y[4];
   const double xi_46 = _data_edge_stencil_at_edge_y[3];
   const double xi_47 = _data_edge_stencil_at_edge_y[1];
   const double xi_48 = _data_edge_stencil_at_edge_y[2];
   const double xi_49 = _data_vertex_stencil_at_edge_y[3];
   const double xi_50 = _data_vertex_stencil_at_edge_y[0];
   const double xi_51 = _data_vertex_stencil_at_edge_y[2];
   const double xi_52 = _data_vertex_stencil_at_edge_y[1];
   const double xi_114 = _data_edge_stencil_at_edge_x[0];
   const double xi_133 = 1 / (xi_114);
   const double xi_115 = _data_edge_stencil_at_edge_x[1];
   const double xi_116 = _data_edge_stencil_at_edge_x[3];
   const double xi_117 = _data_edge_stencil_at_edge_x[4];
   const double xi_118 = _data_edge_stencil_at_edge_x[2];
   const double xi_119 = _data_vertex_stencil_at_edge_x[1];
   const double xi_120 = _data_vertex_stencil_at_edge_x[3];
   const double xi_121 = _data_vertex_stencil_at_edge_x[2];
   const double xi_122 = _data_vertex_stencil_at_edge_x[0];
   const double xi_123 = _data_edge_stencil_at_edge_xy[0];
   const double xi_134 = 1 / (xi_123);
   const double xi_124 = _data_edge_stencil_at_edge_xy[1];
   const double xi_125 = _data_edge_stencil_at_edge_xy[4];
   const double xi_126 = _data_edge_stencil_at_edge_xy[3];
   const double xi_127 = _data_edge_stencil_at_edge_xy[2];
   const double xi_128 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_129 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_130 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_131 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_0 = _data_edge_stencil_at_edge_x[0];
   const double xi_28 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_edge_x[1];
   const double xi_2 = _data_edge_stencil_at_edge_x[3];
   const double xi_3 = _data_edge_stencil_at_edge_x[4];
   const double xi_4 = _data_edge_stencil_at_edge_x[2];
   const double xi_5 = _data_vertex_stencil_at_edge_x[1];
   const double xi_6 = _data_vertex_stencil_at_edge_x[3];
   const double xi_7 = _data_vertex_stencil_at_edge_x[2];
   const double xi_8 = _data_vertex_stencil_at_edge_x[0];
   const double xi_9 = _data_edge_stencil_at_edge_xy[0];
   const double xi_29 = 1 / (xi_9);
   const double xi_10 = _data_edge_stencil_at_edge_xy[1];
   const double xi_11 = _data_edge_stencil_at_edge_xy[4];
   const double xi_12 = _data_edge_stencil_at_edge_xy[3];
   const double xi_13 = _data_edge_stencil_at_edge_xy[2];
   const double xi_14 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_15 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_16 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_17 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_18 = _data_edge_stencil_at_edge_y[0];
   const double xi_30 = 1 / (xi_18);
   const double xi_19 = _data_edge_stencil_at_edge_y[4];
   const double xi_20 = _data_edge_stencil_at_edge_y[3];
   const double xi_21 = _data_edge_stencil_at_edge_y[1];
   const double xi_22 = _data_edge_stencil_at_edge_y[2];
   const double xi_23 = _data_vertex_stencil_at_edge_y[3];
   const double xi_24 = _data_vertex_stencil_at_edge_y[0];
   const double xi_25 = _data_vertex_stencil_at_edge_y[2];
   const double xi_26 = _data_vertex_stencil_at_edge_y[1];
   const double xi_92 = _data_edge_stencil_at_edge_x[0];
   const double xi_111 = 1 / (xi_92);
   const double xi_93 = _data_edge_stencil_at_edge_x[1];
   const double xi_94 = _data_edge_stencil_at_edge_x[3];
   const double xi_95 = _data_edge_stencil_at_edge_x[4];
   const double xi_96 = _data_edge_stencil_at_edge_x[2];
   const double xi_97 = _data_vertex_stencil_at_edge_x[1];
   const double xi_98 = _data_vertex_stencil_at_edge_x[3];
   const double xi_99 = _data_vertex_stencil_at_edge_x[2];
   const double xi_100 = _data_vertex_stencil_at_edge_x[0];
   const double xi_101 = _data_edge_stencil_at_edge_y[0];
   const double xi_112 = 1 / (xi_101);
   const double xi_102 = _data_edge_stencil_at_edge_y[4];
   const double xi_103 = _data_edge_stencil_at_edge_y[3];
   const double xi_104 = _data_edge_stencil_at_edge_y[1];
   const double xi_105 = _data_edge_stencil_at_edge_y[2];
   const double xi_106 = _data_vertex_stencil_at_edge_y[3];
   const double xi_107 = _data_vertex_stencil_at_edge_y[0];
   const double xi_108 = _data_vertex_stencil_at_edge_y[2];
   const double xi_109 = _data_vertex_stencil_at_edge_y[1];
   const double xi_56 = _data_edge_stencil_at_edge_x[0];
   const double xi_66 = 1 / (xi_56);
   const double xi_57 = _data_edge_stencil_at_edge_x[1];
   const double xi_58 = _data_edge_stencil_at_edge_x[3];
   const double xi_59 = _data_edge_stencil_at_edge_x[4];
   const double xi_60 = _data_edge_stencil_at_edge_x[2];
   const double xi_61 = _data_vertex_stencil_at_edge_x[1];
   const double xi_62 = _data_vertex_stencil_at_edge_x[3];
   const double xi_63 = _data_vertex_stencil_at_edge_x[2];
   const double xi_64 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_42*(-xi_33*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_34*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_35*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] - xi_36*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026] - xi_39*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_40*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1027] + _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 1023; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_89*(-xi_71*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_72*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_73*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] - xi_74*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_75*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_76*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026] - xi_77*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_78*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1027] + _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_90*(-xi_80*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_81*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024] - xi_82*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_83*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_84*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] - xi_85*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_86*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026] - xi_87*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1023; ctr_1 < 1024; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_54*(-xi_45*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_46*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024] - xi_47*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_48*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_49*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] - xi_50*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_51*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026] - xi_52*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 1023; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_133*(-xi_115*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] - xi_116*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_117*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_118*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024] - xi_119*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_120*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026] - xi_121*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] - xi_122*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_134*(-xi_124*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_125*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_126*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] - xi_127*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_128*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_129*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026] - xi_130*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_131*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1027] + _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1023; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_28*(-xi_1*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] - xi_2*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024] - xi_5*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026] - xi_7*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] - xi_8*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_29*(-xi_10*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] - xi_13*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_14*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026] - xi_16*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_17*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1027] + _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_30*(-xi_19*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024] - xi_21*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] - xi_24*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026] - xi_26*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 1023; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_111*(-xi_100*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_93*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] - xi_94*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_95*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_96*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024] - xi_97*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_98*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026] - xi_99*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] + _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_112*(-xi_102*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_103*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1024] - xi_104*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_105*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_106*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1025] - xi_107*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_108*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026] - xi_109*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1023; ctr_2 < 1024; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_66*(-xi_57*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] - xi_58*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + ((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_59*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_60*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 + 2*((1049600) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 1024] - xi_61*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_62*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1026] - xi_63*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 1025] - xi_64*_data_vertexFaceDst[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 1025*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_11(double * _data_edgeFaceDst, double * _data_edgeFaceRhs, double * const _data_edge_stencil_at_edge_x, double * const _data_edge_stencil_at_edge_xy, double * const _data_edge_stencil_at_edge_y, double * _data_vertexFaceDst, double * const _data_vertex_stencil_at_edge_x, double * const _data_vertex_stencil_at_edge_xy, double * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_32 = _data_edge_stencil_at_edge_xy[0];
   const double xi_42 = 1 / (xi_32);
   const double xi_33 = _data_edge_stencil_at_edge_xy[1];
   const double xi_34 = _data_edge_stencil_at_edge_xy[4];
   const double xi_35 = _data_edge_stencil_at_edge_xy[3];
   const double xi_36 = _data_edge_stencil_at_edge_xy[2];
   const double xi_37 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_38 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_39 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_40 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_70 = _data_edge_stencil_at_edge_xy[0];
   const double xi_89 = 1 / (xi_70);
   const double xi_71 = _data_edge_stencil_at_edge_xy[1];
   const double xi_72 = _data_edge_stencil_at_edge_xy[4];
   const double xi_73 = _data_edge_stencil_at_edge_xy[3];
   const double xi_74 = _data_edge_stencil_at_edge_xy[2];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_76 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_77 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_78 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_79 = _data_edge_stencil_at_edge_y[0];
   const double xi_90 = 1 / (xi_79);
   const double xi_80 = _data_edge_stencil_at_edge_y[4];
   const double xi_81 = _data_edge_stencil_at_edge_y[3];
   const double xi_82 = _data_edge_stencil_at_edge_y[1];
   const double xi_83 = _data_edge_stencil_at_edge_y[2];
   const double xi_84 = _data_vertex_stencil_at_edge_y[3];
   const double xi_85 = _data_vertex_stencil_at_edge_y[0];
   const double xi_86 = _data_vertex_stencil_at_edge_y[2];
   const double xi_87 = _data_vertex_stencil_at_edge_y[1];
   const double xi_44 = _data_edge_stencil_at_edge_y[0];
   const double xi_54 = 1 / (xi_44);
   const double xi_45 = _data_edge_stencil_at_edge_y[4];
   const double xi_46 = _data_edge_stencil_at_edge_y[3];
   const double xi_47 = _data_edge_stencil_at_edge_y[1];
   const double xi_48 = _data_edge_stencil_at_edge_y[2];
   const double xi_49 = _data_vertex_stencil_at_edge_y[3];
   const double xi_50 = _data_vertex_stencil_at_edge_y[0];
   const double xi_51 = _data_vertex_stencil_at_edge_y[2];
   const double xi_52 = _data_vertex_stencil_at_edge_y[1];
   const double xi_114 = _data_edge_stencil_at_edge_x[0];
   const double xi_133 = 1 / (xi_114);
   const double xi_115 = _data_edge_stencil_at_edge_x[1];
   const double xi_116 = _data_edge_stencil_at_edge_x[3];
   const double xi_117 = _data_edge_stencil_at_edge_x[4];
   const double xi_118 = _data_edge_stencil_at_edge_x[2];
   const double xi_119 = _data_vertex_stencil_at_edge_x[1];
   const double xi_120 = _data_vertex_stencil_at_edge_x[3];
   const double xi_121 = _data_vertex_stencil_at_edge_x[2];
   const double xi_122 = _data_vertex_stencil_at_edge_x[0];
   const double xi_123 = _data_edge_stencil_at_edge_xy[0];
   const double xi_134 = 1 / (xi_123);
   const double xi_124 = _data_edge_stencil_at_edge_xy[1];
   const double xi_125 = _data_edge_stencil_at_edge_xy[4];
   const double xi_126 = _data_edge_stencil_at_edge_xy[3];
   const double xi_127 = _data_edge_stencil_at_edge_xy[2];
   const double xi_128 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_129 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_130 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_131 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_0 = _data_edge_stencil_at_edge_x[0];
   const double xi_28 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_edge_x[1];
   const double xi_2 = _data_edge_stencil_at_edge_x[3];
   const double xi_3 = _data_edge_stencil_at_edge_x[4];
   const double xi_4 = _data_edge_stencil_at_edge_x[2];
   const double xi_5 = _data_vertex_stencil_at_edge_x[1];
   const double xi_6 = _data_vertex_stencil_at_edge_x[3];
   const double xi_7 = _data_vertex_stencil_at_edge_x[2];
   const double xi_8 = _data_vertex_stencil_at_edge_x[0];
   const double xi_9 = _data_edge_stencil_at_edge_xy[0];
   const double xi_29 = 1 / (xi_9);
   const double xi_10 = _data_edge_stencil_at_edge_xy[1];
   const double xi_11 = _data_edge_stencil_at_edge_xy[4];
   const double xi_12 = _data_edge_stencil_at_edge_xy[3];
   const double xi_13 = _data_edge_stencil_at_edge_xy[2];
   const double xi_14 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_15 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_16 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_17 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_18 = _data_edge_stencil_at_edge_y[0];
   const double xi_30 = 1 / (xi_18);
   const double xi_19 = _data_edge_stencil_at_edge_y[4];
   const double xi_20 = _data_edge_stencil_at_edge_y[3];
   const double xi_21 = _data_edge_stencil_at_edge_y[1];
   const double xi_22 = _data_edge_stencil_at_edge_y[2];
   const double xi_23 = _data_vertex_stencil_at_edge_y[3];
   const double xi_24 = _data_vertex_stencil_at_edge_y[0];
   const double xi_25 = _data_vertex_stencil_at_edge_y[2];
   const double xi_26 = _data_vertex_stencil_at_edge_y[1];
   const double xi_92 = _data_edge_stencil_at_edge_x[0];
   const double xi_111 = 1 / (xi_92);
   const double xi_93 = _data_edge_stencil_at_edge_x[1];
   const double xi_94 = _data_edge_stencil_at_edge_x[3];
   const double xi_95 = _data_edge_stencil_at_edge_x[4];
   const double xi_96 = _data_edge_stencil_at_edge_x[2];
   const double xi_97 = _data_vertex_stencil_at_edge_x[1];
   const double xi_98 = _data_vertex_stencil_at_edge_x[3];
   const double xi_99 = _data_vertex_stencil_at_edge_x[2];
   const double xi_100 = _data_vertex_stencil_at_edge_x[0];
   const double xi_101 = _data_edge_stencil_at_edge_y[0];
   const double xi_112 = 1 / (xi_101);
   const double xi_102 = _data_edge_stencil_at_edge_y[4];
   const double xi_103 = _data_edge_stencil_at_edge_y[3];
   const double xi_104 = _data_edge_stencil_at_edge_y[1];
   const double xi_105 = _data_edge_stencil_at_edge_y[2];
   const double xi_106 = _data_vertex_stencil_at_edge_y[3];
   const double xi_107 = _data_vertex_stencil_at_edge_y[0];
   const double xi_108 = _data_vertex_stencil_at_edge_y[2];
   const double xi_109 = _data_vertex_stencil_at_edge_y[1];
   const double xi_56 = _data_edge_stencil_at_edge_x[0];
   const double xi_66 = 1 / (xi_56);
   const double xi_57 = _data_edge_stencil_at_edge_x[1];
   const double xi_58 = _data_edge_stencil_at_edge_x[3];
   const double xi_59 = _data_edge_stencil_at_edge_x[4];
   const double xi_60 = _data_edge_stencil_at_edge_x[2];
   const double xi_61 = _data_vertex_stencil_at_edge_x[1];
   const double xi_62 = _data_vertex_stencil_at_edge_x[3];
   const double xi_63 = _data_vertex_stencil_at_edge_x[2];
   const double xi_64 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_42*(-xi_33*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_34*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_35*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] - xi_36*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050] - xi_39*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_40*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2051] + _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 2047; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_89*(-xi_71*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_72*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_73*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] - xi_74*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_75*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_76*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050] - xi_77*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_78*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2051] + _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_90*(-xi_80*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_81*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048] - xi_82*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_83*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_84*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] - xi_85*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_86*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050] - xi_87*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 2047; ctr_1 < 2048; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_54*(-xi_45*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_46*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048] - xi_47*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_48*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_49*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] - xi_50*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_51*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050] - xi_52*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 2047; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_133*(-xi_115*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] - xi_116*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_117*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_118*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048] - xi_119*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_120*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050] - xi_121*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] - xi_122*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_134*(-xi_124*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_125*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_126*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] - xi_127*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_128*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_129*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050] - xi_130*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_131*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2051] + _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2047; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_28*(-xi_1*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] - xi_2*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048] - xi_5*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050] - xi_7*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] - xi_8*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_29*(-xi_10*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] - xi_13*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_14*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050] - xi_16*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_17*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2051] + _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_30*(-xi_19*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048] - xi_21*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] - xi_24*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050] - xi_26*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 2047; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_111*(-xi_100*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_93*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] - xi_94*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_95*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_96*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048] - xi_97*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_98*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050] - xi_99*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] + _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_112*(-xi_102*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_103*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2048] - xi_104*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_105*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_106*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2049] - xi_107*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_108*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050] - xi_109*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 2047; ctr_2 < 2048; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_66*(-xi_57*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] - xi_58*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + ((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_59*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_60*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 + 2*((4196352) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 2048] - xi_61*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_62*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 2050] - xi_63*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 2049] - xi_64*_data_vertexFaceDst[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 2049*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_12(double * _data_edgeFaceDst, double * _data_edgeFaceRhs, double * const _data_edge_stencil_at_edge_x, double * const _data_edge_stencil_at_edge_xy, double * const _data_edge_stencil_at_edge_y, double * _data_vertexFaceDst, double * const _data_vertex_stencil_at_edge_x, double * const _data_vertex_stencil_at_edge_xy, double * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_32 = _data_edge_stencil_at_edge_xy[0];
   const double xi_42 = 1 / (xi_32);
   const double xi_33 = _data_edge_stencil_at_edge_xy[1];
   const double xi_34 = _data_edge_stencil_at_edge_xy[4];
   const double xi_35 = _data_edge_stencil_at_edge_xy[3];
   const double xi_36 = _data_edge_stencil_at_edge_xy[2];
   const double xi_37 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_38 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_39 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_40 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_70 = _data_edge_stencil_at_edge_xy[0];
   const double xi_89 = 1 / (xi_70);
   const double xi_71 = _data_edge_stencil_at_edge_xy[1];
   const double xi_72 = _data_edge_stencil_at_edge_xy[4];
   const double xi_73 = _data_edge_stencil_at_edge_xy[3];
   const double xi_74 = _data_edge_stencil_at_edge_xy[2];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_76 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_77 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_78 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_79 = _data_edge_stencil_at_edge_y[0];
   const double xi_90 = 1 / (xi_79);
   const double xi_80 = _data_edge_stencil_at_edge_y[4];
   const double xi_81 = _data_edge_stencil_at_edge_y[3];
   const double xi_82 = _data_edge_stencil_at_edge_y[1];
   const double xi_83 = _data_edge_stencil_at_edge_y[2];
   const double xi_84 = _data_vertex_stencil_at_edge_y[3];
   const double xi_85 = _data_vertex_stencil_at_edge_y[0];
   const double xi_86 = _data_vertex_stencil_at_edge_y[2];
   const double xi_87 = _data_vertex_stencil_at_edge_y[1];
   const double xi_44 = _data_edge_stencil_at_edge_y[0];
   const double xi_54 = 1 / (xi_44);
   const double xi_45 = _data_edge_stencil_at_edge_y[4];
   const double xi_46 = _data_edge_stencil_at_edge_y[3];
   const double xi_47 = _data_edge_stencil_at_edge_y[1];
   const double xi_48 = _data_edge_stencil_at_edge_y[2];
   const double xi_49 = _data_vertex_stencil_at_edge_y[3];
   const double xi_50 = _data_vertex_stencil_at_edge_y[0];
   const double xi_51 = _data_vertex_stencil_at_edge_y[2];
   const double xi_52 = _data_vertex_stencil_at_edge_y[1];
   const double xi_114 = _data_edge_stencil_at_edge_x[0];
   const double xi_133 = 1 / (xi_114);
   const double xi_115 = _data_edge_stencil_at_edge_x[1];
   const double xi_116 = _data_edge_stencil_at_edge_x[3];
   const double xi_117 = _data_edge_stencil_at_edge_x[4];
   const double xi_118 = _data_edge_stencil_at_edge_x[2];
   const double xi_119 = _data_vertex_stencil_at_edge_x[1];
   const double xi_120 = _data_vertex_stencil_at_edge_x[3];
   const double xi_121 = _data_vertex_stencil_at_edge_x[2];
   const double xi_122 = _data_vertex_stencil_at_edge_x[0];
   const double xi_123 = _data_edge_stencil_at_edge_xy[0];
   const double xi_134 = 1 / (xi_123);
   const double xi_124 = _data_edge_stencil_at_edge_xy[1];
   const double xi_125 = _data_edge_stencil_at_edge_xy[4];
   const double xi_126 = _data_edge_stencil_at_edge_xy[3];
   const double xi_127 = _data_edge_stencil_at_edge_xy[2];
   const double xi_128 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_129 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_130 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_131 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_0 = _data_edge_stencil_at_edge_x[0];
   const double xi_28 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_edge_x[1];
   const double xi_2 = _data_edge_stencil_at_edge_x[3];
   const double xi_3 = _data_edge_stencil_at_edge_x[4];
   const double xi_4 = _data_edge_stencil_at_edge_x[2];
   const double xi_5 = _data_vertex_stencil_at_edge_x[1];
   const double xi_6 = _data_vertex_stencil_at_edge_x[3];
   const double xi_7 = _data_vertex_stencil_at_edge_x[2];
   const double xi_8 = _data_vertex_stencil_at_edge_x[0];
   const double xi_9 = _data_edge_stencil_at_edge_xy[0];
   const double xi_29 = 1 / (xi_9);
   const double xi_10 = _data_edge_stencil_at_edge_xy[1];
   const double xi_11 = _data_edge_stencil_at_edge_xy[4];
   const double xi_12 = _data_edge_stencil_at_edge_xy[3];
   const double xi_13 = _data_edge_stencil_at_edge_xy[2];
   const double xi_14 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_15 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_16 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_17 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_18 = _data_edge_stencil_at_edge_y[0];
   const double xi_30 = 1 / (xi_18);
   const double xi_19 = _data_edge_stencil_at_edge_y[4];
   const double xi_20 = _data_edge_stencil_at_edge_y[3];
   const double xi_21 = _data_edge_stencil_at_edge_y[1];
   const double xi_22 = _data_edge_stencil_at_edge_y[2];
   const double xi_23 = _data_vertex_stencil_at_edge_y[3];
   const double xi_24 = _data_vertex_stencil_at_edge_y[0];
   const double xi_25 = _data_vertex_stencil_at_edge_y[2];
   const double xi_26 = _data_vertex_stencil_at_edge_y[1];
   const double xi_92 = _data_edge_stencil_at_edge_x[0];
   const double xi_111 = 1 / (xi_92);
   const double xi_93 = _data_edge_stencil_at_edge_x[1];
   const double xi_94 = _data_edge_stencil_at_edge_x[3];
   const double xi_95 = _data_edge_stencil_at_edge_x[4];
   const double xi_96 = _data_edge_stencil_at_edge_x[2];
   const double xi_97 = _data_vertex_stencil_at_edge_x[1];
   const double xi_98 = _data_vertex_stencil_at_edge_x[3];
   const double xi_99 = _data_vertex_stencil_at_edge_x[2];
   const double xi_100 = _data_vertex_stencil_at_edge_x[0];
   const double xi_101 = _data_edge_stencil_at_edge_y[0];
   const double xi_112 = 1 / (xi_101);
   const double xi_102 = _data_edge_stencil_at_edge_y[4];
   const double xi_103 = _data_edge_stencil_at_edge_y[3];
   const double xi_104 = _data_edge_stencil_at_edge_y[1];
   const double xi_105 = _data_edge_stencil_at_edge_y[2];
   const double xi_106 = _data_vertex_stencil_at_edge_y[3];
   const double xi_107 = _data_vertex_stencil_at_edge_y[0];
   const double xi_108 = _data_vertex_stencil_at_edge_y[2];
   const double xi_109 = _data_vertex_stencil_at_edge_y[1];
   const double xi_56 = _data_edge_stencil_at_edge_x[0];
   const double xi_66 = 1 / (xi_56);
   const double xi_57 = _data_edge_stencil_at_edge_x[1];
   const double xi_58 = _data_edge_stencil_at_edge_x[3];
   const double xi_59 = _data_edge_stencil_at_edge_x[4];
   const double xi_60 = _data_edge_stencil_at_edge_x[2];
   const double xi_61 = _data_vertex_stencil_at_edge_x[1];
   const double xi_62 = _data_vertex_stencil_at_edge_x[3];
   const double xi_63 = _data_vertex_stencil_at_edge_x[2];
   const double xi_64 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_42*(-xi_33*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_34*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_35*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] - xi_36*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098] - xi_39*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_40*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4099] + _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 4095; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_89*(-xi_71*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_72*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_73*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] - xi_74*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_75*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_76*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098] - xi_77*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_78*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4099] + _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_90*(-xi_80*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_81*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096] - xi_82*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_83*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_84*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] - xi_85*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_86*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098] - xi_87*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 4095; ctr_1 < 4096; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_54*(-xi_45*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_46*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096] - xi_47*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_48*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_49*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] - xi_50*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_51*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098] - xi_52*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 4095; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_133*(-xi_115*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] - xi_116*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_117*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_118*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096] - xi_119*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_120*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098] - xi_121*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] - xi_122*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_134*(-xi_124*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_125*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_126*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] - xi_127*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_128*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_129*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098] - xi_130*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_131*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4099] + _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4095; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_28*(-xi_1*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] - xi_2*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096] - xi_5*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098] - xi_7*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] - xi_8*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_29*(-xi_10*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] - xi_13*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_14*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098] - xi_16*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_17*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4099] + _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_30*(-xi_19*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096] - xi_21*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] - xi_24*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098] - xi_26*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 4095; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_111*(-xi_100*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_93*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] - xi_94*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_95*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_96*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096] - xi_97*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_98*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098] - xi_99*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] + _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_112*(-xi_102*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_103*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4096] - xi_104*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_105*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_106*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4097] - xi_107*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_108*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098] - xi_109*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 4095; ctr_2 < 4096; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_66*(-xi_57*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] - xi_58*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + ((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_59*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_60*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 + 2*((16781312) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 4096] - xi_61*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_62*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 4098] - xi_63*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 4097] - xi_64*_data_vertexFaceDst[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 4097*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_13(double * _data_edgeFaceDst, double * _data_edgeFaceRhs, double * const _data_edge_stencil_at_edge_x, double * const _data_edge_stencil_at_edge_xy, double * const _data_edge_stencil_at_edge_y, double * _data_vertexFaceDst, double * const _data_vertex_stencil_at_edge_x, double * const _data_vertex_stencil_at_edge_xy, double * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_32 = _data_edge_stencil_at_edge_xy[0];
   const double xi_42 = 1 / (xi_32);
   const double xi_33 = _data_edge_stencil_at_edge_xy[1];
   const double xi_34 = _data_edge_stencil_at_edge_xy[4];
   const double xi_35 = _data_edge_stencil_at_edge_xy[3];
   const double xi_36 = _data_edge_stencil_at_edge_xy[2];
   const double xi_37 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_38 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_39 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_40 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_70 = _data_edge_stencil_at_edge_xy[0];
   const double xi_89 = 1 / (xi_70);
   const double xi_71 = _data_edge_stencil_at_edge_xy[1];
   const double xi_72 = _data_edge_stencil_at_edge_xy[4];
   const double xi_73 = _data_edge_stencil_at_edge_xy[3];
   const double xi_74 = _data_edge_stencil_at_edge_xy[2];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_76 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_77 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_78 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_79 = _data_edge_stencil_at_edge_y[0];
   const double xi_90 = 1 / (xi_79);
   const double xi_80 = _data_edge_stencil_at_edge_y[4];
   const double xi_81 = _data_edge_stencil_at_edge_y[3];
   const double xi_82 = _data_edge_stencil_at_edge_y[1];
   const double xi_83 = _data_edge_stencil_at_edge_y[2];
   const double xi_84 = _data_vertex_stencil_at_edge_y[3];
   const double xi_85 = _data_vertex_stencil_at_edge_y[0];
   const double xi_86 = _data_vertex_stencil_at_edge_y[2];
   const double xi_87 = _data_vertex_stencil_at_edge_y[1];
   const double xi_44 = _data_edge_stencil_at_edge_y[0];
   const double xi_54 = 1 / (xi_44);
   const double xi_45 = _data_edge_stencil_at_edge_y[4];
   const double xi_46 = _data_edge_stencil_at_edge_y[3];
   const double xi_47 = _data_edge_stencil_at_edge_y[1];
   const double xi_48 = _data_edge_stencil_at_edge_y[2];
   const double xi_49 = _data_vertex_stencil_at_edge_y[3];
   const double xi_50 = _data_vertex_stencil_at_edge_y[0];
   const double xi_51 = _data_vertex_stencil_at_edge_y[2];
   const double xi_52 = _data_vertex_stencil_at_edge_y[1];
   const double xi_114 = _data_edge_stencil_at_edge_x[0];
   const double xi_133 = 1 / (xi_114);
   const double xi_115 = _data_edge_stencil_at_edge_x[1];
   const double xi_116 = _data_edge_stencil_at_edge_x[3];
   const double xi_117 = _data_edge_stencil_at_edge_x[4];
   const double xi_118 = _data_edge_stencil_at_edge_x[2];
   const double xi_119 = _data_vertex_stencil_at_edge_x[1];
   const double xi_120 = _data_vertex_stencil_at_edge_x[3];
   const double xi_121 = _data_vertex_stencil_at_edge_x[2];
   const double xi_122 = _data_vertex_stencil_at_edge_x[0];
   const double xi_123 = _data_edge_stencil_at_edge_xy[0];
   const double xi_134 = 1 / (xi_123);
   const double xi_124 = _data_edge_stencil_at_edge_xy[1];
   const double xi_125 = _data_edge_stencil_at_edge_xy[4];
   const double xi_126 = _data_edge_stencil_at_edge_xy[3];
   const double xi_127 = _data_edge_stencil_at_edge_xy[2];
   const double xi_128 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_129 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_130 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_131 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_0 = _data_edge_stencil_at_edge_x[0];
   const double xi_28 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_edge_x[1];
   const double xi_2 = _data_edge_stencil_at_edge_x[3];
   const double xi_3 = _data_edge_stencil_at_edge_x[4];
   const double xi_4 = _data_edge_stencil_at_edge_x[2];
   const double xi_5 = _data_vertex_stencil_at_edge_x[1];
   const double xi_6 = _data_vertex_stencil_at_edge_x[3];
   const double xi_7 = _data_vertex_stencil_at_edge_x[2];
   const double xi_8 = _data_vertex_stencil_at_edge_x[0];
   const double xi_9 = _data_edge_stencil_at_edge_xy[0];
   const double xi_29 = 1 / (xi_9);
   const double xi_10 = _data_edge_stencil_at_edge_xy[1];
   const double xi_11 = _data_edge_stencil_at_edge_xy[4];
   const double xi_12 = _data_edge_stencil_at_edge_xy[3];
   const double xi_13 = _data_edge_stencil_at_edge_xy[2];
   const double xi_14 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_15 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_16 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_17 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_18 = _data_edge_stencil_at_edge_y[0];
   const double xi_30 = 1 / (xi_18);
   const double xi_19 = _data_edge_stencil_at_edge_y[4];
   const double xi_20 = _data_edge_stencil_at_edge_y[3];
   const double xi_21 = _data_edge_stencil_at_edge_y[1];
   const double xi_22 = _data_edge_stencil_at_edge_y[2];
   const double xi_23 = _data_vertex_stencil_at_edge_y[3];
   const double xi_24 = _data_vertex_stencil_at_edge_y[0];
   const double xi_25 = _data_vertex_stencil_at_edge_y[2];
   const double xi_26 = _data_vertex_stencil_at_edge_y[1];
   const double xi_92 = _data_edge_stencil_at_edge_x[0];
   const double xi_111 = 1 / (xi_92);
   const double xi_93 = _data_edge_stencil_at_edge_x[1];
   const double xi_94 = _data_edge_stencil_at_edge_x[3];
   const double xi_95 = _data_edge_stencil_at_edge_x[4];
   const double xi_96 = _data_edge_stencil_at_edge_x[2];
   const double xi_97 = _data_vertex_stencil_at_edge_x[1];
   const double xi_98 = _data_vertex_stencil_at_edge_x[3];
   const double xi_99 = _data_vertex_stencil_at_edge_x[2];
   const double xi_100 = _data_vertex_stencil_at_edge_x[0];
   const double xi_101 = _data_edge_stencil_at_edge_y[0];
   const double xi_112 = 1 / (xi_101);
   const double xi_102 = _data_edge_stencil_at_edge_y[4];
   const double xi_103 = _data_edge_stencil_at_edge_y[3];
   const double xi_104 = _data_edge_stencil_at_edge_y[1];
   const double xi_105 = _data_edge_stencil_at_edge_y[2];
   const double xi_106 = _data_vertex_stencil_at_edge_y[3];
   const double xi_107 = _data_vertex_stencil_at_edge_y[0];
   const double xi_108 = _data_vertex_stencil_at_edge_y[2];
   const double xi_109 = _data_vertex_stencil_at_edge_y[1];
   const double xi_56 = _data_edge_stencil_at_edge_x[0];
   const double xi_66 = 1 / (xi_56);
   const double xi_57 = _data_edge_stencil_at_edge_x[1];
   const double xi_58 = _data_edge_stencil_at_edge_x[3];
   const double xi_59 = _data_edge_stencil_at_edge_x[4];
   const double xi_60 = _data_edge_stencil_at_edge_x[2];
   const double xi_61 = _data_vertex_stencil_at_edge_x[1];
   const double xi_62 = _data_vertex_stencil_at_edge_x[3];
   const double xi_63 = _data_vertex_stencil_at_edge_x[2];
   const double xi_64 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_42*(-xi_33*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_34*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_35*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] - xi_36*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194] - xi_39*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_40*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8195] + _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 8191; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_89*(-xi_71*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_72*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_73*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] - xi_74*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_75*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_76*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194] - xi_77*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_78*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8195] + _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_90*(-xi_80*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_81*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192] - xi_82*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_83*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_84*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] - xi_85*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_86*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194] - xi_87*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 8191; ctr_1 < 8192; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_54*(-xi_45*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_46*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192] - xi_47*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_48*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_49*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] - xi_50*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_51*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194] - xi_52*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 8191; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_133*(-xi_115*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] - xi_116*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_117*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_118*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192] - xi_119*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_120*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194] - xi_121*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] - xi_122*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_134*(-xi_124*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_125*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_126*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] - xi_127*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_128*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_129*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194] - xi_130*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_131*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8195] + _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8191; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_28*(-xi_1*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] - xi_2*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192] - xi_5*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194] - xi_7*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] - xi_8*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_29*(-xi_10*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] - xi_13*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_14*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194] - xi_16*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_17*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8195] + _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_30*(-xi_19*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192] - xi_21*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] - xi_24*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194] - xi_26*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 8191; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_111*(-xi_100*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_93*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] - xi_94*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_95*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_96*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192] - xi_97*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_98*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194] - xi_99*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] + _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_112*(-xi_102*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_103*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8192] - xi_104*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_105*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_106*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8193] - xi_107*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_108*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194] - xi_109*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 8191; ctr_2 < 8192; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_66*(-xi_57*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] - xi_58*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + ((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_59*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_60*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 + 2*((67117056) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 8192] - xi_61*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_62*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 8194] - xi_63*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 8193] - xi_64*_data_vertexFaceDst[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 8193*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_14(double * _data_edgeFaceDst, double * _data_edgeFaceRhs, double * const _data_edge_stencil_at_edge_x, double * const _data_edge_stencil_at_edge_xy, double * const _data_edge_stencil_at_edge_y, double * _data_vertexFaceDst, double * const _data_vertex_stencil_at_edge_x, double * const _data_vertex_stencil_at_edge_xy, double * const _data_vertex_stencil_at_edge_y, double relax)
{
   const double xi_32 = _data_edge_stencil_at_edge_xy[0];
   const double xi_42 = 1 / (xi_32);
   const double xi_33 = _data_edge_stencil_at_edge_xy[1];
   const double xi_34 = _data_edge_stencil_at_edge_xy[4];
   const double xi_35 = _data_edge_stencil_at_edge_xy[3];
   const double xi_36 = _data_edge_stencil_at_edge_xy[2];
   const double xi_37 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_38 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_39 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_40 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_70 = _data_edge_stencil_at_edge_xy[0];
   const double xi_89 = 1 / (xi_70);
   const double xi_71 = _data_edge_stencil_at_edge_xy[1];
   const double xi_72 = _data_edge_stencil_at_edge_xy[4];
   const double xi_73 = _data_edge_stencil_at_edge_xy[3];
   const double xi_74 = _data_edge_stencil_at_edge_xy[2];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_76 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_77 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_78 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_79 = _data_edge_stencil_at_edge_y[0];
   const double xi_90 = 1 / (xi_79);
   const double xi_80 = _data_edge_stencil_at_edge_y[4];
   const double xi_81 = _data_edge_stencil_at_edge_y[3];
   const double xi_82 = _data_edge_stencil_at_edge_y[1];
   const double xi_83 = _data_edge_stencil_at_edge_y[2];
   const double xi_84 = _data_vertex_stencil_at_edge_y[3];
   const double xi_85 = _data_vertex_stencil_at_edge_y[0];
   const double xi_86 = _data_vertex_stencil_at_edge_y[2];
   const double xi_87 = _data_vertex_stencil_at_edge_y[1];
   const double xi_44 = _data_edge_stencil_at_edge_y[0];
   const double xi_54 = 1 / (xi_44);
   const double xi_45 = _data_edge_stencil_at_edge_y[4];
   const double xi_46 = _data_edge_stencil_at_edge_y[3];
   const double xi_47 = _data_edge_stencil_at_edge_y[1];
   const double xi_48 = _data_edge_stencil_at_edge_y[2];
   const double xi_49 = _data_vertex_stencil_at_edge_y[3];
   const double xi_50 = _data_vertex_stencil_at_edge_y[0];
   const double xi_51 = _data_vertex_stencil_at_edge_y[2];
   const double xi_52 = _data_vertex_stencil_at_edge_y[1];
   const double xi_114 = _data_edge_stencil_at_edge_x[0];
   const double xi_133 = 1 / (xi_114);
   const double xi_115 = _data_edge_stencil_at_edge_x[1];
   const double xi_116 = _data_edge_stencil_at_edge_x[3];
   const double xi_117 = _data_edge_stencil_at_edge_x[4];
   const double xi_118 = _data_edge_stencil_at_edge_x[2];
   const double xi_119 = _data_vertex_stencil_at_edge_x[1];
   const double xi_120 = _data_vertex_stencil_at_edge_x[3];
   const double xi_121 = _data_vertex_stencil_at_edge_x[2];
   const double xi_122 = _data_vertex_stencil_at_edge_x[0];
   const double xi_123 = _data_edge_stencil_at_edge_xy[0];
   const double xi_134 = 1 / (xi_123);
   const double xi_124 = _data_edge_stencil_at_edge_xy[1];
   const double xi_125 = _data_edge_stencil_at_edge_xy[4];
   const double xi_126 = _data_edge_stencil_at_edge_xy[3];
   const double xi_127 = _data_edge_stencil_at_edge_xy[2];
   const double xi_128 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_129 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_130 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_131 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_0 = _data_edge_stencil_at_edge_x[0];
   const double xi_28 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_edge_x[1];
   const double xi_2 = _data_edge_stencil_at_edge_x[3];
   const double xi_3 = _data_edge_stencil_at_edge_x[4];
   const double xi_4 = _data_edge_stencil_at_edge_x[2];
   const double xi_5 = _data_vertex_stencil_at_edge_x[1];
   const double xi_6 = _data_vertex_stencil_at_edge_x[3];
   const double xi_7 = _data_vertex_stencil_at_edge_x[2];
   const double xi_8 = _data_vertex_stencil_at_edge_x[0];
   const double xi_9 = _data_edge_stencil_at_edge_xy[0];
   const double xi_29 = 1 / (xi_9);
   const double xi_10 = _data_edge_stencil_at_edge_xy[1];
   const double xi_11 = _data_edge_stencil_at_edge_xy[4];
   const double xi_12 = _data_edge_stencil_at_edge_xy[3];
   const double xi_13 = _data_edge_stencil_at_edge_xy[2];
   const double xi_14 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_15 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_16 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_17 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_18 = _data_edge_stencil_at_edge_y[0];
   const double xi_30 = 1 / (xi_18);
   const double xi_19 = _data_edge_stencil_at_edge_y[4];
   const double xi_20 = _data_edge_stencil_at_edge_y[3];
   const double xi_21 = _data_edge_stencil_at_edge_y[1];
   const double xi_22 = _data_edge_stencil_at_edge_y[2];
   const double xi_23 = _data_vertex_stencil_at_edge_y[3];
   const double xi_24 = _data_vertex_stencil_at_edge_y[0];
   const double xi_25 = _data_vertex_stencil_at_edge_y[2];
   const double xi_26 = _data_vertex_stencil_at_edge_y[1];
   const double xi_92 = _data_edge_stencil_at_edge_x[0];
   const double xi_111 = 1 / (xi_92);
   const double xi_93 = _data_edge_stencil_at_edge_x[1];
   const double xi_94 = _data_edge_stencil_at_edge_x[3];
   const double xi_95 = _data_edge_stencil_at_edge_x[4];
   const double xi_96 = _data_edge_stencil_at_edge_x[2];
   const double xi_97 = _data_vertex_stencil_at_edge_x[1];
   const double xi_98 = _data_vertex_stencil_at_edge_x[3];
   const double xi_99 = _data_vertex_stencil_at_edge_x[2];
   const double xi_100 = _data_vertex_stencil_at_edge_x[0];
   const double xi_101 = _data_edge_stencil_at_edge_y[0];
   const double xi_112 = 1 / (xi_101);
   const double xi_102 = _data_edge_stencil_at_edge_y[4];
   const double xi_103 = _data_edge_stencil_at_edge_y[3];
   const double xi_104 = _data_edge_stencil_at_edge_y[1];
   const double xi_105 = _data_edge_stencil_at_edge_y[2];
   const double xi_106 = _data_vertex_stencil_at_edge_y[3];
   const double xi_107 = _data_vertex_stencil_at_edge_y[0];
   const double xi_108 = _data_vertex_stencil_at_edge_y[2];
   const double xi_109 = _data_vertex_stencil_at_edge_y[1];
   const double xi_56 = _data_edge_stencil_at_edge_x[0];
   const double xi_66 = 1 / (xi_56);
   const double xi_57 = _data_edge_stencil_at_edge_x[1];
   const double xi_58 = _data_edge_stencil_at_edge_x[3];
   const double xi_59 = _data_edge_stencil_at_edge_x[4];
   const double xi_60 = _data_edge_stencil_at_edge_x[2];
   const double xi_61 = _data_vertex_stencil_at_edge_x[1];
   const double xi_62 = _data_vertex_stencil_at_edge_x[3];
   const double xi_63 = _data_vertex_stencil_at_edge_x[2];
   const double xi_64 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_42*(-xi_33*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_34*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_35*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] - xi_36*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_37*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386] - xi_39*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_40*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16387] + _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < 16383; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_89*(-xi_71*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_72*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_73*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] - xi_74*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_75*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_76*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386] - xi_77*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_78*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16387] + _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_90*(-xi_80*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_81*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384] - xi_82*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_83*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_84*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] - xi_85*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_86*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386] - xi_87*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 16383; ctr_1 < 16384; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_54*(-xi_45*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_46*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384] - xi_47*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_48*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_49*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] - xi_50*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_51*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386] - xi_52*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 16383; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_133*(-xi_115*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] - xi_116*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_117*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_118*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384] - xi_119*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_120*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386] - xi_121*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] - xi_122*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_134*(-xi_124*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_125*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_126*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] - xi_127*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_128*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_129*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386] - xi_130*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_131*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16387] + _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16383; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_28*(-xi_1*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] - xi_2*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_3*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_4*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384] - xi_5*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386] - xi_7*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] - xi_8*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_29*(-xi_10*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] - xi_13*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_14*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386] - xi_16*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_17*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16387] + _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_30*(-xi_19*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_20*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384] - xi_21*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_23*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] - xi_24*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386] - xi_26*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      for (int ctr_1 = -ctr_2 + 16383; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_111*(-xi_100*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_93*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] - xi_94*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_95*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_96*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384] - xi_97*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_98*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386] - xi_99*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] + _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_112*(-xi_102*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2)) - 1] - xi_103*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16384] - xi_104*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_105*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_106*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16385] - xi_107*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_108*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386] - xi_109*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 16383; ctr_2 < 16384; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_66*(-xi_57*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] - xi_58*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + ((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_59*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_60*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 + 2*((268451840) / (2)) - ((ctr_2*(ctr_2 - 1)) / (2)) - 16384] - xi_61*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_62*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 16386] - xi_63*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 - 1)) / (2)) - 16385] - xi_64*_data_vertexFaceDst[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + 16385*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void sor_2D_macroface_P2_update_edgedofs_level_any(double * _data_edgeFaceDst, double * _data_edgeFaceRhs, double * const _data_edge_stencil_at_edge_x, double * const _data_edge_stencil_at_edge_xy, double * const _data_edge_stencil_at_edge_y, double * _data_vertexFaceDst, double * const _data_vertex_stencil_at_edge_x, double * const _data_vertex_stencil_at_edge_xy, double * const _data_vertex_stencil_at_edge_y, int64_t level, double relax)
{
   const double xi_32 = _data_edge_stencil_at_edge_xy[0];
   const double xi_42 = 1 / (xi_32);
   const double xi_33 = _data_edge_stencil_at_edge_xy[1];
   const double xi_34 = _data_edge_stencil_at_edge_xy[4];
   const double xi_35 = _data_edge_stencil_at_edge_xy[3];
   const double xi_36 = _data_edge_stencil_at_edge_xy[2];
   const double xi_37 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_38 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_39 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_40 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_70 = _data_edge_stencil_at_edge_xy[0];
   const double xi_89 = 1 / (xi_70);
   const double xi_71 = _data_edge_stencil_at_edge_xy[1];
   const double xi_72 = _data_edge_stencil_at_edge_xy[4];
   const double xi_73 = _data_edge_stencil_at_edge_xy[3];
   const double xi_74 = _data_edge_stencil_at_edge_xy[2];
   const double xi_75 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_76 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_77 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_78 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_79 = _data_edge_stencil_at_edge_y[0];
   const double xi_90 = 1 / (xi_79);
   const double xi_80 = _data_edge_stencil_at_edge_y[4];
   const double xi_81 = _data_edge_stencil_at_edge_y[3];
   const double xi_82 = _data_edge_stencil_at_edge_y[1];
   const double xi_83 = _data_edge_stencil_at_edge_y[2];
   const double xi_84 = _data_vertex_stencil_at_edge_y[3];
   const double xi_85 = _data_vertex_stencil_at_edge_y[0];
   const double xi_86 = _data_vertex_stencil_at_edge_y[2];
   const double xi_87 = _data_vertex_stencil_at_edge_y[1];
   const double xi_44 = _data_edge_stencil_at_edge_y[0];
   const double xi_54 = 1 / (xi_44);
   const double xi_45 = _data_edge_stencil_at_edge_y[4];
   const double xi_46 = _data_edge_stencil_at_edge_y[3];
   const double xi_47 = _data_edge_stencil_at_edge_y[1];
   const double xi_48 = _data_edge_stencil_at_edge_y[2];
   const double xi_49 = _data_vertex_stencil_at_edge_y[3];
   const double xi_50 = _data_vertex_stencil_at_edge_y[0];
   const double xi_51 = _data_vertex_stencil_at_edge_y[2];
   const double xi_52 = _data_vertex_stencil_at_edge_y[1];
   const double xi_114 = _data_edge_stencil_at_edge_x[0];
   const double xi_133 = 1 / (xi_114);
   const double xi_115 = _data_edge_stencil_at_edge_x[1];
   const double xi_116 = _data_edge_stencil_at_edge_x[3];
   const double xi_117 = _data_edge_stencil_at_edge_x[4];
   const double xi_118 = _data_edge_stencil_at_edge_x[2];
   const double xi_119 = _data_vertex_stencil_at_edge_x[1];
   const double xi_120 = _data_vertex_stencil_at_edge_x[3];
   const double xi_121 = _data_vertex_stencil_at_edge_x[2];
   const double xi_122 = _data_vertex_stencil_at_edge_x[0];
   const double xi_123 = _data_edge_stencil_at_edge_xy[0];
   const double xi_134 = 1 / (xi_123);
   const double xi_124 = _data_edge_stencil_at_edge_xy[1];
   const double xi_125 = _data_edge_stencil_at_edge_xy[4];
   const double xi_126 = _data_edge_stencil_at_edge_xy[3];
   const double xi_127 = _data_edge_stencil_at_edge_xy[2];
   const double xi_128 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_129 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_130 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_131 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_0 = _data_edge_stencil_at_edge_x[0];
   const double xi_28 = 1 / (xi_0);
   const double xi_1 = _data_edge_stencil_at_edge_x[1];
   const double xi_2 = _data_edge_stencil_at_edge_x[3];
   const double xi_3 = _data_edge_stencil_at_edge_x[4];
   const double xi_4 = _data_edge_stencil_at_edge_x[2];
   const double xi_5 = _data_vertex_stencil_at_edge_x[1];
   const double xi_6 = _data_vertex_stencil_at_edge_x[3];
   const double xi_7 = _data_vertex_stencil_at_edge_x[2];
   const double xi_8 = _data_vertex_stencil_at_edge_x[0];
   const double xi_9 = _data_edge_stencil_at_edge_xy[0];
   const double xi_29 = 1 / (xi_9);
   const double xi_10 = _data_edge_stencil_at_edge_xy[1];
   const double xi_11 = _data_edge_stencil_at_edge_xy[4];
   const double xi_12 = _data_edge_stencil_at_edge_xy[3];
   const double xi_13 = _data_edge_stencil_at_edge_xy[2];
   const double xi_14 = _data_vertex_stencil_at_edge_xy[3];
   const double xi_15 = _data_vertex_stencil_at_edge_xy[2];
   const double xi_16 = _data_vertex_stencil_at_edge_xy[0];
   const double xi_17 = _data_vertex_stencil_at_edge_xy[1];
   const double xi_18 = _data_edge_stencil_at_edge_y[0];
   const double xi_30 = 1 / (xi_18);
   const double xi_19 = _data_edge_stencil_at_edge_y[4];
   const double xi_20 = _data_edge_stencil_at_edge_y[3];
   const double xi_21 = _data_edge_stencil_at_edge_y[1];
   const double xi_22 = _data_edge_stencil_at_edge_y[2];
   const double xi_23 = _data_vertex_stencil_at_edge_y[3];
   const double xi_24 = _data_vertex_stencil_at_edge_y[0];
   const double xi_25 = _data_vertex_stencil_at_edge_y[2];
   const double xi_26 = _data_vertex_stencil_at_edge_y[1];
   const double xi_92 = _data_edge_stencil_at_edge_x[0];
   const double xi_111 = 1 / (xi_92);
   const double xi_93 = _data_edge_stencil_at_edge_x[1];
   const double xi_94 = _data_edge_stencil_at_edge_x[3];
   const double xi_95 = _data_edge_stencil_at_edge_x[4];
   const double xi_96 = _data_edge_stencil_at_edge_x[2];
   const double xi_97 = _data_vertex_stencil_at_edge_x[1];
   const double xi_98 = _data_vertex_stencil_at_edge_x[3];
   const double xi_99 = _data_vertex_stencil_at_edge_x[2];
   const double xi_100 = _data_vertex_stencil_at_edge_x[0];
   const double xi_101 = _data_edge_stencil_at_edge_y[0];
   const double xi_112 = 1 / (xi_101);
   const double xi_102 = _data_edge_stencil_at_edge_y[4];
   const double xi_103 = _data_edge_stencil_at_edge_y[3];
   const double xi_104 = _data_edge_stencil_at_edge_y[1];
   const double xi_105 = _data_edge_stencil_at_edge_y[2];
   const double xi_106 = _data_vertex_stencil_at_edge_y[3];
   const double xi_107 = _data_vertex_stencil_at_edge_y[0];
   const double xi_108 = _data_vertex_stencil_at_edge_y[2];
   const double xi_109 = _data_vertex_stencil_at_edge_y[1];
   const double xi_56 = _data_edge_stencil_at_edge_x[0];
   const double xi_66 = 1 / (xi_56);
   const double xi_57 = _data_edge_stencil_at_edge_x[1];
   const double xi_58 = _data_edge_stencil_at_edge_x[3];
   const double xi_59 = _data_edge_stencil_at_edge_x[4];
   const double xi_60 = _data_edge_stencil_at_edge_x[2];
   const double xi_61 = _data_vertex_stencil_at_edge_x[1];
   const double xi_62 = _data_vertex_stencil_at_edge_x[3];
   const double xi_63 = _data_vertex_stencil_at_edge_x[2];
   const double xi_64 = _data_vertex_stencil_at_edge_x[0];
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] = relax*xi_42*(-xi_33*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_34*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_35*_data_edgeFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_36*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1] - xi_37*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_38*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_39*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_40*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < (1 << (level)) - 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] = relax*xi_89*(-xi_71*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_72*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_73*_data_edgeFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_74*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1] - xi_75*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_76*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_77*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_78*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] = relax*xi_90*(-xi_80*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1] - xi_81*_data_edgeFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_82*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_83*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_84*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_85*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_86*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_87*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
      }
      for (int ctr_1 = (1 << (level)) - 1; ctr_1 < (1 << (level)); ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] = relax*xi_54*(-xi_45*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1] - xi_46*_data_edgeFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_47*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_48*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_49*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_50*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_51*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_52*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < (1 << (level)) - 1; ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_133*(-xi_115*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_116*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_117*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_118*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1] - xi_119*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_120*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_121*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_122*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] = relax*xi_134*(-xi_124*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_125*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_126*_data_edgeFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_127*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1] - xi_128*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_129*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_130*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_131*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
      }
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (level)) - 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_28*(-xi_1*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_2*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_3*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_4*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1] - xi_5*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_6*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_7*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_8*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] = relax*xi_29*(-xi_10*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_11*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_12*_data_edgeFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_13*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1] - xi_14*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_15*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_16*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_17*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] = relax*xi_30*(-xi_19*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1] - xi_20*_data_edgeFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_21*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_22*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_23*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_24*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_25*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_26*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
      }
      for (int ctr_1 = -ctr_2 + (1 << (level)) - 1; ctr_1 < -ctr_2 + (1 << (level)); ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_111*(-xi_100*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] - xi_93*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_94*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_95*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_96*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1] - xi_97*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_98*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_99*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] = relax*xi_112*(-xi_102*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2)) - 1] - xi_103*_data_edgeFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 1) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_104*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_105*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_106*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2)) - 1] - xi_107*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_108*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_109*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))];
      }
   }
   for (int ctr_2 = (1 << (level)) - 1; ctr_2 < (1 << (level)); ctr_2 += 1)
   {
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         _data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))] = relax*xi_66*(-xi_57*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_58*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + ((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_59*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2))] - xi_60*_data_edgeFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 1) - ((ctr_2*(ctr_2 - 1)) / (2)) + 2*((((1 << (level)) + 1)*(1 << (level))) / (2)) + 1] - xi_61*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))] - xi_62*_data_vertexFaceDst[ctr_1 + (ctr_2 + 1)*((1 << (level)) + 2) - (((ctr_2 + 1)*(ctr_2 + 2)) / (2))] - xi_63*_data_vertexFaceDst[ctr_1 + (ctr_2 - 1)*((1 << (level)) + 2) - ((ctr_2*(ctr_2 - 1)) / (2)) + 1] - xi_64*_data_vertexFaceDst[ctr_1 + ctr_2*((1 << (level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2)) + 1] + _data_edgeFaceRhs[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))]) + (-relax + 1.0)*_data_edgeFaceDst[ctr_1 + ctr_2*((1 << (level)) + 1) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}


void sor_2D_macroface_P2_update_edgedofs(double * _data_edgeFaceDst, double * _data_edgeFaceRhs, double * const _data_edge_stencil_at_edge_x, double * const _data_edge_stencil_at_edge_xy, double * const _data_edge_stencil_at_edge_y, double * _data_vertexFaceDst, double * const _data_vertex_stencil_at_edge_x, double * const _data_vertex_stencil_at_edge_xy, double * const _data_vertex_stencil_at_edge_y, int64_t level, double relax)
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
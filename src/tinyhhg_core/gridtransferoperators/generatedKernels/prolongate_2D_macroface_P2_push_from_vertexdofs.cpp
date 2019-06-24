
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "prolongate_2D_macroface_P2_push_from_vertexdofs.hpp"

namespace hhg {
namespace P2 {
namespace macroface {
namespace generated {

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_any(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, int32_t coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_54 = 1 / (num_neighbor_faces_edge0);
   const double xi_55 = 1 / (num_neighbor_faces_edge2);
   const double xi_56 = 1 / (num_neighbor_faces_vertex0);
   const double xi_74 = 1 / (num_neighbor_faces_edge1);
   const double xi_75 = 1 / (num_neighbor_faces_vertex1);
   const double xi_94 = 1 / (num_neighbor_faces_vertex2);
   {
      for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
      {
         // bottom left vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_59 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_61 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_63 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_65 = _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_68 = -0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_67 = _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_69 = _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_58 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_64 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_58 + xi_59;
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_60 + xi_61;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_62 + xi_63;
            _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_64 + xi_65;
            _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_67 + xi_68;
            _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_68 + xi_69;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_56*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < (1 << (coarse_level)); ctr_1 += 1)
         {
            const double xi_118 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_120 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_122 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_124 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_143 = -0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_126 = _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_128 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_130 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_132 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_134 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_136 = _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_138 = _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_141 = 0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_140 = _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_142 = _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_144 = _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_123 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_121 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_118 + xi_123;
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_120 + xi_121;
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_121 + xi_122;
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_123 + xi_124;
            _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2] = xi_126 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2] = xi_128 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1] = xi_130 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1] = xi_132 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_134 + xi_143;
            _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_136 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_138 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_140 + xi_141;
            _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_141 + xi_142;
            _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_143 + xi_144;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_54*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = (1 << (coarse_level)); ctr_1 < (1 << (coarse_level)) + 1; ctr_1 += 1)
         {
            const double xi_78 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_80 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_82 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_84 = _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_87 = -0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_86 = _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_88 = _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_77 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_79 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_81 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_83 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_77 + xi_78;
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_79 + xi_80;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2] = xi_81 + xi_82;
            _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
            _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2] = xi_86 + xi_87;
            _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_87 + xi_88;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_75*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < (1 << (coarse_level)); ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_151 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_153 = _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_155 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_157 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_176 = -0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_159 = _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_166 = 0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_161 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_163 = _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_165 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_167 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_169 = _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_171 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_173 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_175 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_177 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_156 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_154 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_151 + xi_156;
            _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_153 + xi_154;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_154 + xi_155;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))] = xi_156 + xi_157;
            _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_159 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_161 + xi_166;
            _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_163 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_165 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_166 + xi_167;
            _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_169 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_171 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))] = xi_173 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1] = xi_175 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1] = xi_176 + xi_177;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_74*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (coarse_level)); ctr_1 += 1)
         {
            const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_4 = _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_6 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_8 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_10 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1];
            const double xi_12 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_14 = _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_16 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_18 = _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_22 = _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_24 = _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_26 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_28 = _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_30 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_32 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1];
            const double xi_34 = _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1];
            const double xi_36 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_38 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_40 = _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_42 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_44 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_46 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_48 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_50 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2] = xi_4 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2] = xi_49 + xi_6;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1] = xi_49 + xi_8;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 1] = xi_10 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_12 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_14 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1] = xi_32 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1] = xi_34 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_36 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_37 + xi_38;
            _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_40 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_42 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))] = xi_44 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))] = xi_46 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1] = xi_48 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1] = xi_49 + xi_50;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = -ctr_2 + (1 << (coarse_level)); ctr_1 < -ctr_2 + (1 << (coarse_level)) + 1; ctr_1 += 1)
         {
            const double xi_184 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_186 = _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_188 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_190 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_209 = -0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_192 = _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_194 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_196 = _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_203 = 0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_198 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_200 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1];
            const double xi_202 = _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1];
            const double xi_204 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_206 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_208 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_210 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_189 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_187 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2] = xi_184 + xi_189;
            _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_186 + xi_187;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_187 + xi_188;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1] = xi_189 + xi_190;
            _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2] = xi_192 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_194 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_196 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_198 + xi_203;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1] = xi_200 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1] = xi_202 + xi_209;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_203 + xi_204;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))] = xi_206 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))] = xi_208 + xi_209;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1] = xi_209 + xi_210;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_55*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = (1 << (coarse_level)); ctr_2 < (1 << (coarse_level)) + 1; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_97 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_99 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_101 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_103 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_106 = -0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_105 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_107 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_96 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_98 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_100 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_102 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_96 + xi_97;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1] = xi_98 + xi_99;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_100 + xi_101;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))] = xi_102 + xi_103;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))] = xi_105 + xi_106;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1] = xi_106 + xi_107;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_94*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}


void prolongate_2D_macroface_P2_push_from_vertexdofs(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, int32_t coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
    switch( coarse_level )
    {

    default:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_any(_data_edgeFineDst_X, _data_edgeFineDst_XY, _data_edgeFineDst_Y, _data_vertexCoarseSrc, _data_vertexFineDst, coarse_level, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hhg

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsP2MacroFace2D.hpp"

namespace hhg {
namespace P2 {
namespace macroface {
namespace generated {

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_2(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
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
            const double xi_59 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_61 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_63 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
            const double xi_65 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_68 = -0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_67 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
            const double xi_69 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_58 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_64 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_58 + xi_59;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_60 + xi_61;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9] = xi_62 + xi_63;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_64 + xi_65;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9] = xi_67 + xi_68;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_68 + xi_69;
            _data_vertexFineDst[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_56*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < 4; ctr_1 += 1)
         {
            const double xi_118 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_120 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_122 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_124 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_143 = -0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_126 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7];
            const double xi_128 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7];
            const double xi_130 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
            const double xi_132 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
            const double xi_134 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
            const double xi_136 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
            const double xi_138 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_141 = 0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_140 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_142 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_144 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_123 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_121 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_118 + xi_123;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_120 + xi_121;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_121 + xi_122;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_123 + xi_124;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7] = xi_126 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7] = xi_128 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8] = xi_130 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8] = xi_132 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9] = xi_134 + xi_143;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9] = xi_136 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_138 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_140 + xi_141;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_141 + xi_142;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_143 + xi_144;
            _data_vertexFineDst[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_54*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = 4; ctr_1 < 5; ctr_1 += 1)
         {
            const double xi_78 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_80 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_82 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7];
            const double xi_84 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_87 = -0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_86 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7];
            const double xi_88 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_77 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_79 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_81 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_83 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_77 + xi_78;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_79 + xi_80;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7] = xi_81 + xi_82;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7] = xi_86 + xi_87;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_87 + xi_88;
            _data_vertexFineDst[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_75*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < 4; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_184 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
            const double xi_186 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_188 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9];
            const double xi_190 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18];
            const double xi_209 = -0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_192 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
            const double xi_199 = 0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_196 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_198 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_200 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9];
            const double xi_202 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
            const double xi_204 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
            const double xi_206 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18];
            const double xi_208 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17];
            const double xi_210 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17];
            const double xi_189 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_187 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9] = xi_184 + xi_189;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_186 + xi_187;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9] = xi_187 + xi_188;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18] = xi_189 + xi_190;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9] = xi_192 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_194 + xi_199;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_196 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_198 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9] = xi_199 + xi_200;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8] = xi_202 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8] = xi_204 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18] = xi_206 + xi_209;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17] = xi_208 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17] = xi_209 + xi_210;
            _data_vertexFineDst[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_55*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < 4 - ctr_2; ctr_1 += 1)
         {
            const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_4 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7];
            const double xi_6 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7];
            const double xi_8 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
            const double xi_10 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
            const double xi_12 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
            const double xi_14 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
            const double xi_16 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_18 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_22 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_24 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_26 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_28 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_30 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_32 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 10];
            const double xi_34 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 10];
            const double xi_36 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9];
            const double xi_38 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9];
            const double xi_40 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
            const double xi_42 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
            const double xi_44 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18];
            const double xi_46 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18];
            const double xi_48 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17];
            const double xi_50 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17];
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7] = xi_4 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7] = xi_49 + xi_6;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8] = xi_49 + xi_8;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8] = xi_10 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9] = xi_12 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9] = xi_14 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 10] = xi_32 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 10] = xi_34 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9] = xi_36 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9] = xi_37 + xi_38;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8] = xi_40 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8] = xi_42 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18] = xi_44 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18] = xi_46 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17] = xi_48 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17] = xi_49 + xi_50;
            _data_vertexFineDst[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = 4 - ctr_2; ctr_1 < 5 - ctr_2; ctr_1 += 1)
         {
            const double xi_151 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7];
            const double xi_153 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_155 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9];
            const double xi_157 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17];
            const double xi_176 = -0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_159 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7];
            const double xi_161 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_163 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_170 = 0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_165 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_167 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 10];
            const double xi_169 = _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 10];
            const double xi_171 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9];
            const double xi_173 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18];
            const double xi_175 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18];
            const double xi_177 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17];
            const double xi_156 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_154 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7] = xi_151 + xi_156;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_153 + xi_154;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9] = xi_154 + xi_155;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17] = xi_156 + xi_157;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7] = xi_159 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_161 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_163 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_165 + xi_170;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 10] = xi_167 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 10] = xi_169 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9] = xi_170 + xi_171;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18] = xi_173 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18] = xi_175 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17] = xi_176 + xi_177;
            _data_vertexFineDst[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_74*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 4; ctr_2 < 5; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_97 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9];
            const double xi_99 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17];
            const double xi_101 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9];
            const double xi_103 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18];
            const double xi_106 = -0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_105 = _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18];
            const double xi_107 = _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17];
            const double xi_96 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_98 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_100 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_102 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9] = xi_96 + xi_97;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17] = xi_98 + xi_99;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9] = xi_100 + xi_101;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18] = xi_102 + xi_103;
            _data_edgeFineDst_XY[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18] = xi_105 + xi_106;
            _data_edgeFineDst_Y[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17] = xi_106 + xi_107;
            _data_vertexFineDst[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_94*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_3(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
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
            const double xi_59 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_61 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_63 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
            const double xi_65 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_68 = -0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_67 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
            const double xi_69 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_58 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_64 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_58 + xi_59;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_60 + xi_61;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17] = xi_62 + xi_63;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_64 + xi_65;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17] = xi_67 + xi_68;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_68 + xi_69;
            _data_vertexFineDst[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_56*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < 8; ctr_1 += 1)
         {
            const double xi_118 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_120 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_122 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_124 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_143 = -0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_126 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15];
            const double xi_128 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15];
            const double xi_130 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
            const double xi_132 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
            const double xi_134 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
            const double xi_136 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
            const double xi_138 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_141 = 0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_140 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_142 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_144 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_123 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_121 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_118 + xi_123;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_120 + xi_121;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_121 + xi_122;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_123 + xi_124;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15] = xi_126 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15] = xi_128 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16] = xi_130 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16] = xi_132 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17] = xi_134 + xi_143;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17] = xi_136 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_138 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_140 + xi_141;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_141 + xi_142;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_143 + xi_144;
            _data_vertexFineDst[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_54*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = 8; ctr_1 < 9; ctr_1 += 1)
         {
            const double xi_78 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_80 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_82 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15];
            const double xi_84 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_87 = -0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_86 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15];
            const double xi_88 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_77 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_79 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_81 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_83 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_77 + xi_78;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_79 + xi_80;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15] = xi_81 + xi_82;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15] = xi_86 + xi_87;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_87 + xi_88;
            _data_vertexFineDst[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_75*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < 8; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_184 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
            const double xi_186 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_188 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17];
            const double xi_190 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34];
            const double xi_209 = -0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_192 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
            const double xi_199 = 0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_196 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_198 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_200 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17];
            const double xi_202 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
            const double xi_204 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
            const double xi_206 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34];
            const double xi_208 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33];
            const double xi_210 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33];
            const double xi_189 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_187 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17] = xi_184 + xi_189;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_186 + xi_187;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17] = xi_187 + xi_188;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34] = xi_189 + xi_190;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17] = xi_192 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_194 + xi_199;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_196 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_198 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17] = xi_199 + xi_200;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16] = xi_202 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16] = xi_204 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34] = xi_206 + xi_209;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33] = xi_208 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33] = xi_209 + xi_210;
            _data_vertexFineDst[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_55*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < 8 - ctr_2; ctr_1 += 1)
         {
            const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_4 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15];
            const double xi_6 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15];
            const double xi_8 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
            const double xi_10 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
            const double xi_12 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
            const double xi_14 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
            const double xi_16 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_18 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_22 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_24 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_26 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_28 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_30 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_32 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 18];
            const double xi_34 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 18];
            const double xi_36 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17];
            const double xi_38 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17];
            const double xi_40 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
            const double xi_42 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
            const double xi_44 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34];
            const double xi_46 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34];
            const double xi_48 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33];
            const double xi_50 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33];
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15] = xi_4 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15] = xi_49 + xi_6;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16] = xi_49 + xi_8;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16] = xi_10 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17] = xi_12 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17] = xi_14 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 18] = xi_32 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 18] = xi_34 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17] = xi_36 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17] = xi_37 + xi_38;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16] = xi_40 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16] = xi_42 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34] = xi_44 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34] = xi_46 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33] = xi_48 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33] = xi_49 + xi_50;
            _data_vertexFineDst[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = 8 - ctr_2; ctr_1 < 9 - ctr_2; ctr_1 += 1)
         {
            const double xi_151 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15];
            const double xi_153 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_155 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17];
            const double xi_157 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33];
            const double xi_176 = -0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_159 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15];
            const double xi_161 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_163 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_170 = 0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_165 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_167 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 18];
            const double xi_169 = _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 18];
            const double xi_171 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17];
            const double xi_173 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34];
            const double xi_175 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34];
            const double xi_177 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33];
            const double xi_156 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_154 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15] = xi_151 + xi_156;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_153 + xi_154;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17] = xi_154 + xi_155;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33] = xi_156 + xi_157;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15] = xi_159 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_161 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_163 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_165 + xi_170;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 18] = xi_167 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 18] = xi_169 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17] = xi_170 + xi_171;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34] = xi_173 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34] = xi_175 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33] = xi_176 + xi_177;
            _data_vertexFineDst[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_74*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 8; ctr_2 < 9; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_97 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17];
            const double xi_99 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33];
            const double xi_101 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17];
            const double xi_103 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34];
            const double xi_106 = -0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_105 = _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34];
            const double xi_107 = _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33];
            const double xi_96 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_98 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_100 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_102 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17] = xi_96 + xi_97;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33] = xi_98 + xi_99;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17] = xi_100 + xi_101;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34] = xi_102 + xi_103;
            _data_edgeFineDst_XY[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34] = xi_105 + xi_106;
            _data_edgeFineDst_Y[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33] = xi_106 + xi_107;
            _data_vertexFineDst[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_94*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_4(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
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
            const double xi_59 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_61 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_63 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
            const double xi_65 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_68 = -0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_67 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
            const double xi_69 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_58 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_64 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_58 + xi_59;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_60 + xi_61;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33] = xi_62 + xi_63;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_64 + xi_65;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33] = xi_67 + xi_68;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_68 + xi_69;
            _data_vertexFineDst[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_56*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < 16; ctr_1 += 1)
         {
            const double xi_118 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_120 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_122 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_124 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_143 = -0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_126 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31];
            const double xi_128 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31];
            const double xi_130 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
            const double xi_132 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
            const double xi_134 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
            const double xi_136 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
            const double xi_138 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_141 = 0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_140 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_142 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_144 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_123 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_121 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_118 + xi_123;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_120 + xi_121;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_121 + xi_122;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_123 + xi_124;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31] = xi_126 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31] = xi_128 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32] = xi_130 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32] = xi_132 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33] = xi_134 + xi_143;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33] = xi_136 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_138 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_140 + xi_141;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_141 + xi_142;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_143 + xi_144;
            _data_vertexFineDst[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_54*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = 16; ctr_1 < 17; ctr_1 += 1)
         {
            const double xi_78 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_80 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_82 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31];
            const double xi_84 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_87 = -0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_86 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31];
            const double xi_88 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_77 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_79 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_81 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_83 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_77 + xi_78;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_79 + xi_80;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31] = xi_81 + xi_82;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31] = xi_86 + xi_87;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_87 + xi_88;
            _data_vertexFineDst[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_75*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < 16; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_184 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
            const double xi_186 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_188 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33];
            const double xi_190 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66];
            const double xi_209 = -0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_192 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
            const double xi_199 = 0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_196 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_198 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_200 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33];
            const double xi_202 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
            const double xi_204 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
            const double xi_206 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66];
            const double xi_208 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65];
            const double xi_210 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65];
            const double xi_189 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_187 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33] = xi_184 + xi_189;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_186 + xi_187;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33] = xi_187 + xi_188;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66] = xi_189 + xi_190;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33] = xi_192 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_194 + xi_199;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_196 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_198 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33] = xi_199 + xi_200;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32] = xi_202 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32] = xi_204 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66] = xi_206 + xi_209;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65] = xi_208 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65] = xi_209 + xi_210;
            _data_vertexFineDst[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_55*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < 16 - ctr_2; ctr_1 += 1)
         {
            const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_4 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31];
            const double xi_6 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31];
            const double xi_8 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
            const double xi_10 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
            const double xi_12 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
            const double xi_14 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
            const double xi_16 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_18 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_22 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_24 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_26 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_28 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_30 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_32 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 34];
            const double xi_34 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 34];
            const double xi_36 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33];
            const double xi_38 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33];
            const double xi_40 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
            const double xi_42 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
            const double xi_44 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66];
            const double xi_46 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66];
            const double xi_48 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65];
            const double xi_50 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65];
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31] = xi_4 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31] = xi_49 + xi_6;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32] = xi_49 + xi_8;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32] = xi_10 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33] = xi_12 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33] = xi_14 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 34] = xi_32 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 34] = xi_34 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33] = xi_36 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33] = xi_37 + xi_38;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32] = xi_40 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32] = xi_42 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66] = xi_44 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66] = xi_46 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65] = xi_48 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65] = xi_49 + xi_50;
            _data_vertexFineDst[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = 16 - ctr_2; ctr_1 < 17 - ctr_2; ctr_1 += 1)
         {
            const double xi_151 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31];
            const double xi_153 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_155 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33];
            const double xi_157 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65];
            const double xi_176 = -0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_159 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31];
            const double xi_161 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_163 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_170 = 0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_165 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_167 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 34];
            const double xi_169 = _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 34];
            const double xi_171 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33];
            const double xi_173 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66];
            const double xi_175 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66];
            const double xi_177 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65];
            const double xi_156 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_154 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31] = xi_151 + xi_156;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_153 + xi_154;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33] = xi_154 + xi_155;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65] = xi_156 + xi_157;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31] = xi_159 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_161 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_163 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_165 + xi_170;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 34] = xi_167 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 34] = xi_169 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33] = xi_170 + xi_171;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66] = xi_173 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66] = xi_175 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65] = xi_176 + xi_177;
            _data_vertexFineDst[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_74*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 16; ctr_2 < 17; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_97 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33];
            const double xi_99 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65];
            const double xi_101 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33];
            const double xi_103 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66];
            const double xi_106 = -0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_105 = _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66];
            const double xi_107 = _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65];
            const double xi_96 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_98 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_100 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_102 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33] = xi_96 + xi_97;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65] = xi_98 + xi_99;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33] = xi_100 + xi_101;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66] = xi_102 + xi_103;
            _data_edgeFineDst_XY[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66] = xi_105 + xi_106;
            _data_edgeFineDst_Y[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65] = xi_106 + xi_107;
            _data_vertexFineDst[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_94*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_5(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
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
            const double xi_59 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_61 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_63 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
            const double xi_65 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_68 = -0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_67 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
            const double xi_69 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_58 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_64 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_58 + xi_59;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_60 + xi_61;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65] = xi_62 + xi_63;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_64 + xi_65;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65] = xi_67 + xi_68;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_68 + xi_69;
            _data_vertexFineDst[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_56*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < 32; ctr_1 += 1)
         {
            const double xi_118 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_120 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_122 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_124 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_143 = -0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_126 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63];
            const double xi_128 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63];
            const double xi_130 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
            const double xi_132 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
            const double xi_134 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
            const double xi_136 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
            const double xi_138 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_141 = 0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_140 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_142 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_144 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_123 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_121 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_118 + xi_123;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_120 + xi_121;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_121 + xi_122;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_123 + xi_124;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63] = xi_126 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63] = xi_128 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64] = xi_130 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64] = xi_132 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65] = xi_134 + xi_143;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65] = xi_136 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_138 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_140 + xi_141;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_141 + xi_142;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_143 + xi_144;
            _data_vertexFineDst[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_54*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = 32; ctr_1 < 33; ctr_1 += 1)
         {
            const double xi_78 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_80 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_82 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63];
            const double xi_84 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_87 = -0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_86 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63];
            const double xi_88 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_77 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_79 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_81 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_83 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_77 + xi_78;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_79 + xi_80;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63] = xi_81 + xi_82;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63] = xi_86 + xi_87;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_87 + xi_88;
            _data_vertexFineDst[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_75*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < 32; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_184 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
            const double xi_186 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_188 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65];
            const double xi_190 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130];
            const double xi_209 = -0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_192 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
            const double xi_199 = 0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_196 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_198 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_200 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65];
            const double xi_202 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
            const double xi_204 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
            const double xi_206 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130];
            const double xi_208 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129];
            const double xi_210 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129];
            const double xi_189 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_187 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65] = xi_184 + xi_189;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_186 + xi_187;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65] = xi_187 + xi_188;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130] = xi_189 + xi_190;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65] = xi_192 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_194 + xi_199;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_196 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_198 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65] = xi_199 + xi_200;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64] = xi_202 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64] = xi_204 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130] = xi_206 + xi_209;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129] = xi_208 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129] = xi_209 + xi_210;
            _data_vertexFineDst[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_55*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < 32 - ctr_2; ctr_1 += 1)
         {
            const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_4 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63];
            const double xi_6 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63];
            const double xi_8 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
            const double xi_10 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
            const double xi_12 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
            const double xi_14 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
            const double xi_16 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_18 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_22 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_24 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_26 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_28 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_30 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_32 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 66];
            const double xi_34 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 66];
            const double xi_36 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65];
            const double xi_38 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65];
            const double xi_40 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
            const double xi_42 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
            const double xi_44 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130];
            const double xi_46 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130];
            const double xi_48 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129];
            const double xi_50 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129];
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63] = xi_4 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63] = xi_49 + xi_6;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64] = xi_49 + xi_8;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64] = xi_10 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65] = xi_12 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65] = xi_14 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 66] = xi_32 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 66] = xi_34 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65] = xi_36 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65] = xi_37 + xi_38;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64] = xi_40 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64] = xi_42 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130] = xi_44 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130] = xi_46 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129] = xi_48 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129] = xi_49 + xi_50;
            _data_vertexFineDst[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = 32 - ctr_2; ctr_1 < 33 - ctr_2; ctr_1 += 1)
         {
            const double xi_151 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63];
            const double xi_153 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_155 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65];
            const double xi_157 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129];
            const double xi_176 = -0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_159 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63];
            const double xi_161 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_163 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_170 = 0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_165 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_167 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 66];
            const double xi_169 = _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 66];
            const double xi_171 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65];
            const double xi_173 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130];
            const double xi_175 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130];
            const double xi_177 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129];
            const double xi_156 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_154 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63] = xi_151 + xi_156;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_153 + xi_154;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65] = xi_154 + xi_155;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129] = xi_156 + xi_157;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63] = xi_159 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_161 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_163 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_165 + xi_170;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 66] = xi_167 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 66] = xi_169 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65] = xi_170 + xi_171;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130] = xi_173 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130] = xi_175 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129] = xi_176 + xi_177;
            _data_vertexFineDst[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_74*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 32; ctr_2 < 33; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_97 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65];
            const double xi_99 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129];
            const double xi_101 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65];
            const double xi_103 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130];
            const double xi_106 = -0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_105 = _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130];
            const double xi_107 = _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129];
            const double xi_96 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_98 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_100 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_102 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65] = xi_96 + xi_97;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129] = xi_98 + xi_99;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65] = xi_100 + xi_101;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130] = xi_102 + xi_103;
            _data_edgeFineDst_XY[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130] = xi_105 + xi_106;
            _data_edgeFineDst_Y[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129] = xi_106 + xi_107;
            _data_vertexFineDst[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_94*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_6(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
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
            const double xi_59 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_61 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_63 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
            const double xi_65 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_68 = -0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_67 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
            const double xi_69 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_58 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_64 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_58 + xi_59;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_60 + xi_61;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129] = xi_62 + xi_63;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_64 + xi_65;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129] = xi_67 + xi_68;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_68 + xi_69;
            _data_vertexFineDst[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_56*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < 64; ctr_1 += 1)
         {
            const double xi_118 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_120 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_122 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_124 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_143 = -0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_126 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127];
            const double xi_128 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127];
            const double xi_130 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
            const double xi_132 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
            const double xi_134 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
            const double xi_136 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
            const double xi_138 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_141 = 0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_140 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_142 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_144 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_123 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_121 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_118 + xi_123;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_120 + xi_121;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_121 + xi_122;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_123 + xi_124;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127] = xi_126 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127] = xi_128 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128] = xi_130 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128] = xi_132 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129] = xi_134 + xi_143;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129] = xi_136 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_138 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_140 + xi_141;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_141 + xi_142;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_143 + xi_144;
            _data_vertexFineDst[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_54*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = 64; ctr_1 < 65; ctr_1 += 1)
         {
            const double xi_78 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_80 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_82 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127];
            const double xi_84 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_87 = -0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_86 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127];
            const double xi_88 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_77 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_79 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_81 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_83 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_77 + xi_78;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_79 + xi_80;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127] = xi_81 + xi_82;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127] = xi_86 + xi_87;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_87 + xi_88;
            _data_vertexFineDst[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_75*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < 64; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_184 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
            const double xi_186 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_188 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129];
            const double xi_190 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258];
            const double xi_209 = -0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_192 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
            const double xi_199 = 0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_196 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_198 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_200 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129];
            const double xi_202 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
            const double xi_204 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
            const double xi_206 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258];
            const double xi_208 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257];
            const double xi_210 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257];
            const double xi_189 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_187 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129] = xi_184 + xi_189;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_186 + xi_187;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129] = xi_187 + xi_188;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258] = xi_189 + xi_190;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129] = xi_192 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_194 + xi_199;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_196 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_198 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129] = xi_199 + xi_200;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128] = xi_202 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128] = xi_204 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258] = xi_206 + xi_209;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257] = xi_208 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257] = xi_209 + xi_210;
            _data_vertexFineDst[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_55*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < 64 - ctr_2; ctr_1 += 1)
         {
            const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_4 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127];
            const double xi_6 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127];
            const double xi_8 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
            const double xi_10 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
            const double xi_12 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
            const double xi_14 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
            const double xi_16 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_18 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_22 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_24 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_26 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_28 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_30 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_32 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 130];
            const double xi_34 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 130];
            const double xi_36 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129];
            const double xi_38 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129];
            const double xi_40 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
            const double xi_42 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
            const double xi_44 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258];
            const double xi_46 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258];
            const double xi_48 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257];
            const double xi_50 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257];
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127] = xi_4 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127] = xi_49 + xi_6;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128] = xi_49 + xi_8;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128] = xi_10 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129] = xi_12 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129] = xi_14 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 130] = xi_32 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 130] = xi_34 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129] = xi_36 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129] = xi_37 + xi_38;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128] = xi_40 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128] = xi_42 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258] = xi_44 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258] = xi_46 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257] = xi_48 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257] = xi_49 + xi_50;
            _data_vertexFineDst[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = 64 - ctr_2; ctr_1 < 65 - ctr_2; ctr_1 += 1)
         {
            const double xi_151 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127];
            const double xi_153 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_155 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129];
            const double xi_157 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257];
            const double xi_176 = -0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_159 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127];
            const double xi_161 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_163 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_170 = 0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_165 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_167 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 130];
            const double xi_169 = _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 130];
            const double xi_171 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129];
            const double xi_173 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258];
            const double xi_175 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258];
            const double xi_177 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257];
            const double xi_156 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_154 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127] = xi_151 + xi_156;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_153 + xi_154;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129] = xi_154 + xi_155;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257] = xi_156 + xi_157;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127] = xi_159 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_161 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_163 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_165 + xi_170;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 130] = xi_167 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 130] = xi_169 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129] = xi_170 + xi_171;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258] = xi_173 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258] = xi_175 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257] = xi_176 + xi_177;
            _data_vertexFineDst[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_74*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 64; ctr_2 < 65; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_97 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129];
            const double xi_99 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257];
            const double xi_101 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129];
            const double xi_103 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258];
            const double xi_106 = -0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_105 = _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258];
            const double xi_107 = _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257];
            const double xi_96 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_98 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_100 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_102 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129] = xi_96 + xi_97;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257] = xi_98 + xi_99;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129] = xi_100 + xi_101;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258] = xi_102 + xi_103;
            _data_edgeFineDst_XY[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258] = xi_105 + xi_106;
            _data_edgeFineDst_Y[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257] = xi_106 + xi_107;
            _data_vertexFineDst[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_94*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_7(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
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
            const double xi_59 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_61 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_63 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
            const double xi_65 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_68 = -0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_67 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
            const double xi_69 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_58 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_64 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_58 + xi_59;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_60 + xi_61;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257] = xi_62 + xi_63;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_64 + xi_65;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257] = xi_67 + xi_68;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_68 + xi_69;
            _data_vertexFineDst[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_56*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < 128; ctr_1 += 1)
         {
            const double xi_118 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_120 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_122 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_124 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_143 = -0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_126 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255];
            const double xi_128 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255];
            const double xi_130 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
            const double xi_132 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
            const double xi_134 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
            const double xi_136 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
            const double xi_138 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_141 = 0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_140 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_142 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_144 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_123 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_121 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_118 + xi_123;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_120 + xi_121;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_121 + xi_122;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_123 + xi_124;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255] = xi_126 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255] = xi_128 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256] = xi_130 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256] = xi_132 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257] = xi_134 + xi_143;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257] = xi_136 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_138 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_140 + xi_141;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_141 + xi_142;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_143 + xi_144;
            _data_vertexFineDst[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_54*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = 128; ctr_1 < 129; ctr_1 += 1)
         {
            const double xi_78 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_80 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_82 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255];
            const double xi_84 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_87 = -0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_86 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255];
            const double xi_88 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_77 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_79 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_81 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_83 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_77 + xi_78;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_79 + xi_80;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255] = xi_81 + xi_82;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255] = xi_86 + xi_87;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_87 + xi_88;
            _data_vertexFineDst[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_75*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < 128; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_184 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
            const double xi_186 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_188 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257];
            const double xi_190 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514];
            const double xi_209 = -0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_192 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
            const double xi_199 = 0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_196 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_198 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_200 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257];
            const double xi_202 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
            const double xi_204 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
            const double xi_206 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514];
            const double xi_208 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513];
            const double xi_210 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513];
            const double xi_189 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_187 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257] = xi_184 + xi_189;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_186 + xi_187;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257] = xi_187 + xi_188;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514] = xi_189 + xi_190;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257] = xi_192 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_194 + xi_199;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_196 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_198 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257] = xi_199 + xi_200;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256] = xi_202 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256] = xi_204 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514] = xi_206 + xi_209;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513] = xi_208 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513] = xi_209 + xi_210;
            _data_vertexFineDst[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_55*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < 128 - ctr_2; ctr_1 += 1)
         {
            const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_4 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255];
            const double xi_6 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255];
            const double xi_8 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
            const double xi_10 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
            const double xi_12 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
            const double xi_14 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
            const double xi_16 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_18 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_22 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_24 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_26 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_28 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_30 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_32 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 258];
            const double xi_34 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 258];
            const double xi_36 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257];
            const double xi_38 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257];
            const double xi_40 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
            const double xi_42 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
            const double xi_44 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514];
            const double xi_46 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514];
            const double xi_48 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513];
            const double xi_50 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513];
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255] = xi_4 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255] = xi_49 + xi_6;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256] = xi_49 + xi_8;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256] = xi_10 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257] = xi_12 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257] = xi_14 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 258] = xi_32 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 258] = xi_34 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257] = xi_36 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257] = xi_37 + xi_38;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256] = xi_40 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256] = xi_42 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514] = xi_44 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514] = xi_46 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513] = xi_48 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513] = xi_49 + xi_50;
            _data_vertexFineDst[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = 128 - ctr_2; ctr_1 < 129 - ctr_2; ctr_1 += 1)
         {
            const double xi_151 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255];
            const double xi_153 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_155 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257];
            const double xi_157 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513];
            const double xi_176 = -0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_159 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255];
            const double xi_161 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_163 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_170 = 0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_165 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_167 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 258];
            const double xi_169 = _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 258];
            const double xi_171 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257];
            const double xi_173 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514];
            const double xi_175 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514];
            const double xi_177 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513];
            const double xi_156 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_154 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255] = xi_151 + xi_156;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_153 + xi_154;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257] = xi_154 + xi_155;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513] = xi_156 + xi_157;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255] = xi_159 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_161 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_163 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_165 + xi_170;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 258] = xi_167 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 258] = xi_169 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257] = xi_170 + xi_171;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514] = xi_173 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514] = xi_175 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513] = xi_176 + xi_177;
            _data_vertexFineDst[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_74*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 128; ctr_2 < 129; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_97 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257];
            const double xi_99 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513];
            const double xi_101 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257];
            const double xi_103 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514];
            const double xi_106 = -0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_105 = _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514];
            const double xi_107 = _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513];
            const double xi_96 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_98 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_100 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_102 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257] = xi_96 + xi_97;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513] = xi_98 + xi_99;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257] = xi_100 + xi_101;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514] = xi_102 + xi_103;
            _data_edgeFineDst_XY[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514] = xi_105 + xi_106;
            _data_edgeFineDst_Y[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513] = xi_106 + xi_107;
            _data_vertexFineDst[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_94*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_8(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
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
            const double xi_59 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_61 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_63 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
            const double xi_65 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_68 = -0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_67 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
            const double xi_69 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_58 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_64 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_58 + xi_59;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_60 + xi_61;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513] = xi_62 + xi_63;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_64 + xi_65;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513] = xi_67 + xi_68;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_68 + xi_69;
            _data_vertexFineDst[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_56*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < 256; ctr_1 += 1)
         {
            const double xi_118 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_120 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_122 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_124 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_143 = -0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_126 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511];
            const double xi_128 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511];
            const double xi_130 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
            const double xi_132 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
            const double xi_134 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
            const double xi_136 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
            const double xi_138 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_141 = 0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_140 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_142 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_144 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_123 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_121 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_118 + xi_123;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_120 + xi_121;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_121 + xi_122;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_123 + xi_124;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511] = xi_126 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511] = xi_128 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512] = xi_130 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512] = xi_132 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513] = xi_134 + xi_143;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513] = xi_136 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_138 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_140 + xi_141;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_141 + xi_142;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_143 + xi_144;
            _data_vertexFineDst[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_54*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = 256; ctr_1 < 257; ctr_1 += 1)
         {
            const double xi_78 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_80 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_82 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511];
            const double xi_84 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_87 = -0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_86 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511];
            const double xi_88 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_77 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_79 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_81 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_83 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_77 + xi_78;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_79 + xi_80;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511] = xi_81 + xi_82;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511] = xi_86 + xi_87;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_87 + xi_88;
            _data_vertexFineDst[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_75*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < 256; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_184 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
            const double xi_186 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_188 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513];
            const double xi_190 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026];
            const double xi_209 = -0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_192 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
            const double xi_199 = 0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_196 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_198 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_200 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513];
            const double xi_202 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
            const double xi_204 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
            const double xi_206 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026];
            const double xi_208 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025];
            const double xi_210 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025];
            const double xi_189 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_187 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513] = xi_184 + xi_189;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_186 + xi_187;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513] = xi_187 + xi_188;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026] = xi_189 + xi_190;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513] = xi_192 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_194 + xi_199;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_196 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_198 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513] = xi_199 + xi_200;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512] = xi_202 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512] = xi_204 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026] = xi_206 + xi_209;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025] = xi_208 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025] = xi_209 + xi_210;
            _data_vertexFineDst[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_55*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < 256 - ctr_2; ctr_1 += 1)
         {
            const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_4 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511];
            const double xi_6 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511];
            const double xi_8 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
            const double xi_10 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
            const double xi_12 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
            const double xi_14 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
            const double xi_16 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_18 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_22 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_24 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_26 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_28 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_30 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_32 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 514];
            const double xi_34 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 514];
            const double xi_36 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513];
            const double xi_38 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513];
            const double xi_40 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
            const double xi_42 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
            const double xi_44 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026];
            const double xi_46 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026];
            const double xi_48 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025];
            const double xi_50 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025];
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511] = xi_4 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511] = xi_49 + xi_6;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512] = xi_49 + xi_8;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512] = xi_10 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513] = xi_12 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513] = xi_14 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 514] = xi_32 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 514] = xi_34 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513] = xi_36 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513] = xi_37 + xi_38;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512] = xi_40 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512] = xi_42 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026] = xi_44 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026] = xi_46 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025] = xi_48 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025] = xi_49 + xi_50;
            _data_vertexFineDst[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = 256 - ctr_2; ctr_1 < 257 - ctr_2; ctr_1 += 1)
         {
            const double xi_151 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511];
            const double xi_153 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_155 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513];
            const double xi_157 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025];
            const double xi_176 = -0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_159 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511];
            const double xi_161 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_163 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_170 = 0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_165 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_167 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 514];
            const double xi_169 = _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 514];
            const double xi_171 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513];
            const double xi_173 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026];
            const double xi_175 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026];
            const double xi_177 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025];
            const double xi_156 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_154 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511] = xi_151 + xi_156;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_153 + xi_154;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513] = xi_154 + xi_155;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025] = xi_156 + xi_157;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511] = xi_159 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_161 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_163 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_165 + xi_170;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 514] = xi_167 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 514] = xi_169 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513] = xi_170 + xi_171;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026] = xi_173 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026] = xi_175 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025] = xi_176 + xi_177;
            _data_vertexFineDst[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_74*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 256; ctr_2 < 257; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_97 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513];
            const double xi_99 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025];
            const double xi_101 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513];
            const double xi_103 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026];
            const double xi_106 = -0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_105 = _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026];
            const double xi_107 = _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025];
            const double xi_96 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_98 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_100 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_102 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513] = xi_96 + xi_97;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025] = xi_98 + xi_99;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513] = xi_100 + xi_101;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026] = xi_102 + xi_103;
            _data_edgeFineDst_XY[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026] = xi_105 + xi_106;
            _data_edgeFineDst_Y[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025] = xi_106 + xi_107;
            _data_vertexFineDst[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_94*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_9(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
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
            const double xi_59 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_61 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_63 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
            const double xi_65 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_68 = -0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_67 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
            const double xi_69 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_58 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_64 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_58 + xi_59;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_60 + xi_61;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025] = xi_62 + xi_63;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_64 + xi_65;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025] = xi_67 + xi_68;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_68 + xi_69;
            _data_vertexFineDst[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_56*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < 512; ctr_1 += 1)
         {
            const double xi_118 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_120 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_122 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_124 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_143 = -0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_126 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023];
            const double xi_128 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023];
            const double xi_130 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
            const double xi_132 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
            const double xi_134 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
            const double xi_136 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
            const double xi_138 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_141 = 0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_140 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_142 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_144 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_123 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_121 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_118 + xi_123;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_120 + xi_121;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_121 + xi_122;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_123 + xi_124;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023] = xi_126 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023] = xi_128 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024] = xi_130 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024] = xi_132 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025] = xi_134 + xi_143;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025] = xi_136 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_138 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_140 + xi_141;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_141 + xi_142;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_143 + xi_144;
            _data_vertexFineDst[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_54*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = 512; ctr_1 < 513; ctr_1 += 1)
         {
            const double xi_78 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_80 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_82 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023];
            const double xi_84 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_87 = -0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_86 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023];
            const double xi_88 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_77 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_79 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_81 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_83 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_77 + xi_78;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_79 + xi_80;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023] = xi_81 + xi_82;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023] = xi_86 + xi_87;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_87 + xi_88;
            _data_vertexFineDst[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_75*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < 512; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_184 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
            const double xi_186 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_188 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025];
            const double xi_190 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050];
            const double xi_209 = -0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_192 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
            const double xi_199 = 0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_196 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_198 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_200 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025];
            const double xi_202 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
            const double xi_204 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
            const double xi_206 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050];
            const double xi_208 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049];
            const double xi_210 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049];
            const double xi_189 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_187 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025] = xi_184 + xi_189;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_186 + xi_187;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025] = xi_187 + xi_188;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050] = xi_189 + xi_190;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025] = xi_192 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_194 + xi_199;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_196 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_198 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025] = xi_199 + xi_200;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024] = xi_202 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024] = xi_204 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050] = xi_206 + xi_209;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049] = xi_208 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049] = xi_209 + xi_210;
            _data_vertexFineDst[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_55*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < 512 - ctr_2; ctr_1 += 1)
         {
            const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_4 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023];
            const double xi_6 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023];
            const double xi_8 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
            const double xi_10 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
            const double xi_12 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
            const double xi_14 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
            const double xi_16 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_18 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_22 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_24 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_26 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_28 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_30 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_32 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1026];
            const double xi_34 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1026];
            const double xi_36 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025];
            const double xi_38 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025];
            const double xi_40 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
            const double xi_42 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
            const double xi_44 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050];
            const double xi_46 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050];
            const double xi_48 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049];
            const double xi_50 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049];
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023] = xi_4 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023] = xi_49 + xi_6;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024] = xi_49 + xi_8;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024] = xi_10 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025] = xi_12 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025] = xi_14 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1026] = xi_32 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1026] = xi_34 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025] = xi_36 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025] = xi_37 + xi_38;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024] = xi_40 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024] = xi_42 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050] = xi_44 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050] = xi_46 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049] = xi_48 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049] = xi_49 + xi_50;
            _data_vertexFineDst[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = 512 - ctr_2; ctr_1 < 513 - ctr_2; ctr_1 += 1)
         {
            const double xi_151 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023];
            const double xi_153 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_155 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025];
            const double xi_157 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049];
            const double xi_176 = -0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_159 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023];
            const double xi_161 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_163 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_170 = 0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_165 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_167 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1026];
            const double xi_169 = _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1026];
            const double xi_171 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025];
            const double xi_173 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050];
            const double xi_175 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050];
            const double xi_177 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049];
            const double xi_156 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_154 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023] = xi_151 + xi_156;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_153 + xi_154;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025] = xi_154 + xi_155;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049] = xi_156 + xi_157;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023] = xi_159 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_161 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_163 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_165 + xi_170;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1026] = xi_167 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1026] = xi_169 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025] = xi_170 + xi_171;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050] = xi_173 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050] = xi_175 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049] = xi_176 + xi_177;
            _data_vertexFineDst[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_74*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 512; ctr_2 < 513; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_97 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025];
            const double xi_99 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049];
            const double xi_101 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025];
            const double xi_103 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050];
            const double xi_106 = -0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_105 = _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050];
            const double xi_107 = _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049];
            const double xi_96 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_98 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_100 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_102 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025] = xi_96 + xi_97;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049] = xi_98 + xi_99;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025] = xi_100 + xi_101;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050] = xi_102 + xi_103;
            _data_edgeFineDst_XY[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050] = xi_105 + xi_106;
            _data_edgeFineDst_Y[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049] = xi_106 + xi_107;
            _data_vertexFineDst[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_94*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_10(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
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
            const double xi_59 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_61 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_63 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
            const double xi_65 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_68 = -0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_67 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
            const double xi_69 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_58 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_64 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_58 + xi_59;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_60 + xi_61;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049] = xi_62 + xi_63;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_64 + xi_65;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049] = xi_67 + xi_68;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_68 + xi_69;
            _data_vertexFineDst[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_56*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < 1024; ctr_1 += 1)
         {
            const double xi_118 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_120 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_122 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_124 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_143 = -0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_126 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047];
            const double xi_128 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047];
            const double xi_130 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
            const double xi_132 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
            const double xi_134 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
            const double xi_136 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
            const double xi_138 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_141 = 0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_140 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_142 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_144 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_123 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_121 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_118 + xi_123;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_120 + xi_121;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_121 + xi_122;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_123 + xi_124;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047] = xi_126 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047] = xi_128 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048] = xi_130 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048] = xi_132 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049] = xi_134 + xi_143;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049] = xi_136 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_138 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_140 + xi_141;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_141 + xi_142;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_143 + xi_144;
            _data_vertexFineDst[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_54*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = 1024; ctr_1 < 1025; ctr_1 += 1)
         {
            const double xi_78 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_80 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_82 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047];
            const double xi_84 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_87 = -0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_86 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047];
            const double xi_88 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_77 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_79 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_81 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_83 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_77 + xi_78;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_79 + xi_80;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047] = xi_81 + xi_82;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047] = xi_86 + xi_87;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_87 + xi_88;
            _data_vertexFineDst[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_75*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < 1024; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_184 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
            const double xi_186 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_188 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049];
            const double xi_190 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098];
            const double xi_209 = -0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_192 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
            const double xi_199 = 0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_196 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_198 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_200 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049];
            const double xi_202 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
            const double xi_204 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
            const double xi_206 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098];
            const double xi_208 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097];
            const double xi_210 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097];
            const double xi_189 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_187 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049] = xi_184 + xi_189;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_186 + xi_187;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049] = xi_187 + xi_188;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098] = xi_189 + xi_190;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049] = xi_192 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_194 + xi_199;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_196 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_198 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049] = xi_199 + xi_200;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048] = xi_202 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048] = xi_204 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098] = xi_206 + xi_209;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097] = xi_208 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097] = xi_209 + xi_210;
            _data_vertexFineDst[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_55*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < 1024 - ctr_2; ctr_1 += 1)
         {
            const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_4 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047];
            const double xi_6 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047];
            const double xi_8 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
            const double xi_10 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
            const double xi_12 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
            const double xi_14 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
            const double xi_16 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_18 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_22 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_24 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_26 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_28 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_30 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_32 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2050];
            const double xi_34 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2050];
            const double xi_36 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049];
            const double xi_38 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049];
            const double xi_40 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
            const double xi_42 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
            const double xi_44 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098];
            const double xi_46 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098];
            const double xi_48 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097];
            const double xi_50 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097];
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047] = xi_4 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047] = xi_49 + xi_6;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048] = xi_49 + xi_8;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048] = xi_10 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049] = xi_12 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049] = xi_14 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2050] = xi_32 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2050] = xi_34 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049] = xi_36 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049] = xi_37 + xi_38;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048] = xi_40 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048] = xi_42 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098] = xi_44 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098] = xi_46 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097] = xi_48 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097] = xi_49 + xi_50;
            _data_vertexFineDst[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = 1024 - ctr_2; ctr_1 < 1025 - ctr_2; ctr_1 += 1)
         {
            const double xi_151 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047];
            const double xi_153 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_155 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049];
            const double xi_157 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097];
            const double xi_176 = -0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_159 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047];
            const double xi_161 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_163 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_170 = 0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_165 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_167 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2050];
            const double xi_169 = _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2050];
            const double xi_171 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049];
            const double xi_173 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098];
            const double xi_175 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098];
            const double xi_177 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097];
            const double xi_156 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_154 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047] = xi_151 + xi_156;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_153 + xi_154;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049] = xi_154 + xi_155;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097] = xi_156 + xi_157;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047] = xi_159 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_161 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_163 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_165 + xi_170;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2050] = xi_167 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2050] = xi_169 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049] = xi_170 + xi_171;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098] = xi_173 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098] = xi_175 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097] = xi_176 + xi_177;
            _data_vertexFineDst[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_74*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1024; ctr_2 < 1025; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_97 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049];
            const double xi_99 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097];
            const double xi_101 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049];
            const double xi_103 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098];
            const double xi_106 = -0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_105 = _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098];
            const double xi_107 = _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097];
            const double xi_96 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_98 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_100 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_102 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049] = xi_96 + xi_97;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097] = xi_98 + xi_99;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049] = xi_100 + xi_101;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098] = xi_102 + xi_103;
            _data_edgeFineDst_XY[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098] = xi_105 + xi_106;
            _data_edgeFineDst_Y[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097] = xi_106 + xi_107;
            _data_vertexFineDst[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_94*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_11(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
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
            const double xi_59 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_61 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_63 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
            const double xi_65 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_68 = -0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_67 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
            const double xi_69 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_58 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_64 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_58 + xi_59;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_60 + xi_61;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097] = xi_62 + xi_63;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_64 + xi_65;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097] = xi_67 + xi_68;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_68 + xi_69;
            _data_vertexFineDst[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_56*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < 2048; ctr_1 += 1)
         {
            const double xi_118 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_120 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_122 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_124 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_143 = -0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_126 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095];
            const double xi_128 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095];
            const double xi_130 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
            const double xi_132 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
            const double xi_134 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
            const double xi_136 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
            const double xi_138 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_141 = 0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_140 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_142 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_144 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_123 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_121 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_118 + xi_123;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_120 + xi_121;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_121 + xi_122;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_123 + xi_124;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095] = xi_126 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095] = xi_128 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096] = xi_130 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096] = xi_132 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097] = xi_134 + xi_143;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097] = xi_136 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_138 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_140 + xi_141;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_141 + xi_142;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_143 + xi_144;
            _data_vertexFineDst[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_54*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = 2048; ctr_1 < 2049; ctr_1 += 1)
         {
            const double xi_78 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_80 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_82 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095];
            const double xi_84 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_87 = -0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_86 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095];
            const double xi_88 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_77 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_79 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_81 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_83 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_77 + xi_78;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_79 + xi_80;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095] = xi_81 + xi_82;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095] = xi_86 + xi_87;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_87 + xi_88;
            _data_vertexFineDst[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_75*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < 2048; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_184 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
            const double xi_186 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_188 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097];
            const double xi_190 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194];
            const double xi_209 = -0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_192 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
            const double xi_199 = 0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_196 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_198 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_200 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097];
            const double xi_202 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
            const double xi_204 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
            const double xi_206 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194];
            const double xi_208 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193];
            const double xi_210 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193];
            const double xi_189 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_187 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097] = xi_184 + xi_189;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_186 + xi_187;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097] = xi_187 + xi_188;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194] = xi_189 + xi_190;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097] = xi_192 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_194 + xi_199;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_196 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_198 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097] = xi_199 + xi_200;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096] = xi_202 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096] = xi_204 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194] = xi_206 + xi_209;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193] = xi_208 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193] = xi_209 + xi_210;
            _data_vertexFineDst[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_55*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < 2048 - ctr_2; ctr_1 += 1)
         {
            const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_4 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095];
            const double xi_6 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095];
            const double xi_8 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
            const double xi_10 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
            const double xi_12 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
            const double xi_14 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
            const double xi_16 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_18 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_22 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_24 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_26 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_28 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_30 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_32 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4098];
            const double xi_34 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4098];
            const double xi_36 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097];
            const double xi_38 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097];
            const double xi_40 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
            const double xi_42 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
            const double xi_44 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194];
            const double xi_46 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194];
            const double xi_48 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193];
            const double xi_50 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193];
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095] = xi_4 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095] = xi_49 + xi_6;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096] = xi_49 + xi_8;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096] = xi_10 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097] = xi_12 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097] = xi_14 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4098] = xi_32 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4098] = xi_34 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097] = xi_36 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097] = xi_37 + xi_38;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096] = xi_40 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096] = xi_42 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194] = xi_44 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194] = xi_46 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193] = xi_48 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193] = xi_49 + xi_50;
            _data_vertexFineDst[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = 2048 - ctr_2; ctr_1 < 2049 - ctr_2; ctr_1 += 1)
         {
            const double xi_151 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095];
            const double xi_153 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_155 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097];
            const double xi_157 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193];
            const double xi_176 = -0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_159 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095];
            const double xi_161 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_163 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_170 = 0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_165 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_167 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4098];
            const double xi_169 = _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4098];
            const double xi_171 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097];
            const double xi_173 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194];
            const double xi_175 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194];
            const double xi_177 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193];
            const double xi_156 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_154 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095] = xi_151 + xi_156;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_153 + xi_154;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097] = xi_154 + xi_155;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193] = xi_156 + xi_157;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095] = xi_159 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_161 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_163 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_165 + xi_170;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4098] = xi_167 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4098] = xi_169 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097] = xi_170 + xi_171;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194] = xi_173 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194] = xi_175 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193] = xi_176 + xi_177;
            _data_vertexFineDst[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_74*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 2048; ctr_2 < 2049; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_97 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097];
            const double xi_99 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193];
            const double xi_101 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097];
            const double xi_103 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194];
            const double xi_106 = -0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_105 = _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194];
            const double xi_107 = _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193];
            const double xi_96 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_98 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_100 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_102 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097] = xi_96 + xi_97;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193] = xi_98 + xi_99;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097] = xi_100 + xi_101;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194] = xi_102 + xi_103;
            _data_edgeFineDst_XY[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194] = xi_105 + xi_106;
            _data_edgeFineDst_Y[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193] = xi_106 + xi_107;
            _data_vertexFineDst[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_94*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_12(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
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
            const double xi_59 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_61 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_63 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
            const double xi_65 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_68 = -0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_67 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
            const double xi_69 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_58 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_64 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_58 + xi_59;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_60 + xi_61;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193] = xi_62 + xi_63;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_64 + xi_65;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193] = xi_67 + xi_68;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_68 + xi_69;
            _data_vertexFineDst[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_56*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < 4096; ctr_1 += 1)
         {
            const double xi_118 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_120 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_122 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_124 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_143 = -0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_126 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191];
            const double xi_128 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191];
            const double xi_130 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
            const double xi_132 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
            const double xi_134 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
            const double xi_136 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
            const double xi_138 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_141 = 0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_140 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_142 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_144 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_123 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_121 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_118 + xi_123;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_120 + xi_121;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_121 + xi_122;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_123 + xi_124;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191] = xi_126 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191] = xi_128 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192] = xi_130 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192] = xi_132 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193] = xi_134 + xi_143;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193] = xi_136 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_138 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_140 + xi_141;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_141 + xi_142;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_143 + xi_144;
            _data_vertexFineDst[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_54*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = 4096; ctr_1 < 4097; ctr_1 += 1)
         {
            const double xi_78 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_80 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_82 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191];
            const double xi_84 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_87 = -0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_86 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191];
            const double xi_88 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_77 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_79 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_81 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_83 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_77 + xi_78;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_79 + xi_80;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191] = xi_81 + xi_82;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191] = xi_86 + xi_87;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_87 + xi_88;
            _data_vertexFineDst[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_75*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < 4096; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_184 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
            const double xi_186 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_188 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193];
            const double xi_190 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386];
            const double xi_209 = -0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_192 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
            const double xi_199 = 0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_196 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_198 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_200 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193];
            const double xi_202 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
            const double xi_204 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
            const double xi_206 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386];
            const double xi_208 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385];
            const double xi_210 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385];
            const double xi_189 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_187 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193] = xi_184 + xi_189;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_186 + xi_187;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193] = xi_187 + xi_188;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386] = xi_189 + xi_190;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193] = xi_192 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_194 + xi_199;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_196 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_198 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193] = xi_199 + xi_200;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192] = xi_202 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192] = xi_204 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386] = xi_206 + xi_209;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385] = xi_208 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385] = xi_209 + xi_210;
            _data_vertexFineDst[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_55*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < 4096 - ctr_2; ctr_1 += 1)
         {
            const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_4 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191];
            const double xi_6 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191];
            const double xi_8 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
            const double xi_10 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
            const double xi_12 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
            const double xi_14 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
            const double xi_16 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_18 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_22 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_24 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_26 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_28 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_30 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_32 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8194];
            const double xi_34 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8194];
            const double xi_36 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193];
            const double xi_38 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193];
            const double xi_40 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
            const double xi_42 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
            const double xi_44 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386];
            const double xi_46 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386];
            const double xi_48 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385];
            const double xi_50 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385];
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191] = xi_4 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191] = xi_49 + xi_6;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192] = xi_49 + xi_8;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192] = xi_10 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193] = xi_12 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193] = xi_14 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8194] = xi_32 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8194] = xi_34 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193] = xi_36 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193] = xi_37 + xi_38;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192] = xi_40 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192] = xi_42 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386] = xi_44 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386] = xi_46 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385] = xi_48 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385] = xi_49 + xi_50;
            _data_vertexFineDst[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = 4096 - ctr_2; ctr_1 < 4097 - ctr_2; ctr_1 += 1)
         {
            const double xi_151 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191];
            const double xi_153 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_155 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193];
            const double xi_157 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385];
            const double xi_176 = -0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_159 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191];
            const double xi_161 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_163 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_170 = 0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_165 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_167 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8194];
            const double xi_169 = _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8194];
            const double xi_171 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193];
            const double xi_173 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386];
            const double xi_175 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386];
            const double xi_177 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385];
            const double xi_156 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_154 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191] = xi_151 + xi_156;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_153 + xi_154;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193] = xi_154 + xi_155;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385] = xi_156 + xi_157;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191] = xi_159 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_161 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_163 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_165 + xi_170;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8194] = xi_167 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8194] = xi_169 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193] = xi_170 + xi_171;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386] = xi_173 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386] = xi_175 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385] = xi_176 + xi_177;
            _data_vertexFineDst[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_74*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 4096; ctr_2 < 4097; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_97 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193];
            const double xi_99 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385];
            const double xi_101 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193];
            const double xi_103 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386];
            const double xi_106 = -0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_105 = _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386];
            const double xi_107 = _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385];
            const double xi_96 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_98 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_100 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_102 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193] = xi_96 + xi_97;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385] = xi_98 + xi_99;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193] = xi_100 + xi_101;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386] = xi_102 + xi_103;
            _data_edgeFineDst_XY[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386] = xi_105 + xi_106;
            _data_edgeFineDst_Y[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385] = xi_106 + xi_107;
            _data_vertexFineDst[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_94*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_13(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
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
            const double xi_59 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_61 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_63 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
            const double xi_65 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_68 = -0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_67 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
            const double xi_69 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_58 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_64 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_58 + xi_59;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_60 + xi_61;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385] = xi_62 + xi_63;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_64 + xi_65;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385] = xi_67 + xi_68;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_68 + xi_69;
            _data_vertexFineDst[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_56*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < 8192; ctr_1 += 1)
         {
            const double xi_118 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_120 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_122 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_124 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_143 = -0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_126 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383];
            const double xi_128 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383];
            const double xi_130 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
            const double xi_132 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
            const double xi_134 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
            const double xi_136 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
            const double xi_138 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_141 = 0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_140 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_142 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_144 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_123 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_121 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_118 + xi_123;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_120 + xi_121;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_121 + xi_122;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_123 + xi_124;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383] = xi_126 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383] = xi_128 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384] = xi_130 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384] = xi_132 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385] = xi_134 + xi_143;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385] = xi_136 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_138 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_140 + xi_141;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_141 + xi_142;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_143 + xi_144;
            _data_vertexFineDst[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_54*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = 8192; ctr_1 < 8193; ctr_1 += 1)
         {
            const double xi_78 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_80 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_82 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383];
            const double xi_84 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_87 = -0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_86 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383];
            const double xi_88 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_77 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_79 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_81 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_83 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_77 + xi_78;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_79 + xi_80;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383] = xi_81 + xi_82;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383] = xi_86 + xi_87;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_87 + xi_88;
            _data_vertexFineDst[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_75*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < 8192; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_184 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
            const double xi_186 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_188 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385];
            const double xi_190 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770];
            const double xi_209 = -0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_192 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
            const double xi_199 = 0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_196 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_198 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_200 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385];
            const double xi_202 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
            const double xi_204 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
            const double xi_206 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770];
            const double xi_208 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769];
            const double xi_210 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769];
            const double xi_189 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_187 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385] = xi_184 + xi_189;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_186 + xi_187;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385] = xi_187 + xi_188;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770] = xi_189 + xi_190;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385] = xi_192 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_194 + xi_199;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_196 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_198 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385] = xi_199 + xi_200;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384] = xi_202 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384] = xi_204 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770] = xi_206 + xi_209;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769] = xi_208 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769] = xi_209 + xi_210;
            _data_vertexFineDst[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_55*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < 8192 - ctr_2; ctr_1 += 1)
         {
            const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_4 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383];
            const double xi_6 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383];
            const double xi_8 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
            const double xi_10 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
            const double xi_12 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
            const double xi_14 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
            const double xi_16 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_18 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_22 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_24 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_26 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_28 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_30 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_32 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16386];
            const double xi_34 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16386];
            const double xi_36 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385];
            const double xi_38 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385];
            const double xi_40 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
            const double xi_42 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
            const double xi_44 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770];
            const double xi_46 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770];
            const double xi_48 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769];
            const double xi_50 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769];
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383] = xi_4 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383] = xi_49 + xi_6;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384] = xi_49 + xi_8;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384] = xi_10 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385] = xi_12 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385] = xi_14 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16386] = xi_32 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16386] = xi_34 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385] = xi_36 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385] = xi_37 + xi_38;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384] = xi_40 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384] = xi_42 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770] = xi_44 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770] = xi_46 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769] = xi_48 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769] = xi_49 + xi_50;
            _data_vertexFineDst[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = 8192 - ctr_2; ctr_1 < 8193 - ctr_2; ctr_1 += 1)
         {
            const double xi_151 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383];
            const double xi_153 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_155 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385];
            const double xi_157 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769];
            const double xi_176 = -0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_159 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383];
            const double xi_161 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_163 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_170 = 0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_165 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_167 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16386];
            const double xi_169 = _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16386];
            const double xi_171 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385];
            const double xi_173 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770];
            const double xi_175 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770];
            const double xi_177 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769];
            const double xi_156 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_154 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383] = xi_151 + xi_156;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_153 + xi_154;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385] = xi_154 + xi_155;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769] = xi_156 + xi_157;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383] = xi_159 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_161 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_163 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_165 + xi_170;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16386] = xi_167 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16386] = xi_169 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385] = xi_170 + xi_171;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770] = xi_173 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770] = xi_175 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769] = xi_176 + xi_177;
            _data_vertexFineDst[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_74*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 8192; ctr_2 < 8193; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_97 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385];
            const double xi_99 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769];
            const double xi_101 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385];
            const double xi_103 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770];
            const double xi_106 = -0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_105 = _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770];
            const double xi_107 = _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769];
            const double xi_96 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_98 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_100 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_102 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385] = xi_96 + xi_97;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769] = xi_98 + xi_99;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385] = xi_100 + xi_101;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770] = xi_102 + xi_103;
            _data_edgeFineDst_XY[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770] = xi_105 + xi_106;
            _data_edgeFineDst_Y[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769] = xi_106 + xi_107;
            _data_vertexFineDst[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_94*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_14(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
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
            const double xi_59 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_61 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_63 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
            const double xi_65 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_68 = -0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_67 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
            const double xi_69 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_58 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_60 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_64 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_58 + xi_59;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_60 + xi_61;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769] = xi_62 + xi_63;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_64 + xi_65;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769] = xi_67 + xi_68;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_68 + xi_69;
            _data_vertexFineDst[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_56*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom edge
         for (int ctr_1 = 1; ctr_1 < 16384; ctr_1 += 1)
         {
            const double xi_118 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_120 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_122 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_124 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_143 = -0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_126 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767];
            const double xi_128 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767];
            const double xi_130 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
            const double xi_132 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
            const double xi_134 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
            const double xi_136 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
            const double xi_138 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_141 = 0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_140 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_142 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_144 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_123 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_121 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_118 + xi_123;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_120 + xi_121;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_121 + xi_122;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_123 + xi_124;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767] = xi_126 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767] = xi_128 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768] = xi_130 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768] = xi_132 + xi_143;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769] = xi_134 + xi_143;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769] = xi_136 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_138 + xi_143;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_140 + xi_141;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_141 + xi_142;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_143 + xi_144;
            _data_vertexFineDst[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_54*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // bottom right vertex
         for (int ctr_1 = 16384; ctr_1 < 16385; ctr_1 += 1)
         {
            const double xi_78 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_80 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_82 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767];
            const double xi_84 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_87 = -0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_86 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767];
            const double xi_88 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_77 = xi_54*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_79 = xi_54*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_81 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_83 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_77 + xi_78;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_79 + xi_80;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767] = xi_81 + xi_82;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767] = xi_86 + xi_87;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_87 + xi_88;
            _data_vertexFineDst[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_75*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 1; ctr_2 < 16384; ctr_2 += 1)
      {
         // left edge
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_184 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
            const double xi_186 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_188 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769];
            const double xi_190 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538];
            const double xi_209 = -0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_192 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
            const double xi_199 = 0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_196 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_198 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_200 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769];
            const double xi_202 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
            const double xi_204 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
            const double xi_206 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538];
            const double xi_208 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537];
            const double xi_210 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537];
            const double xi_189 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_187 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769] = xi_184 + xi_189;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_186 + xi_187;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769] = xi_187 + xi_188;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538] = xi_189 + xi_190;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769] = xi_192 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_194 + xi_199;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_196 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_198 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769] = xi_199 + xi_200;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768] = xi_202 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768] = xi_204 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538] = xi_206 + xi_209;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537] = xi_208 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537] = xi_209 + xi_210;
            _data_vertexFineDst[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_55*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // inner triangle
         for (int ctr_1 = 1; ctr_1 < 16384 - ctr_2; ctr_1 += 1)
         {
            const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_4 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767];
            const double xi_6 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767];
            const double xi_8 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
            const double xi_10 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
            const double xi_12 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
            const double xi_14 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
            const double xi_16 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_18 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_20 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_22 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_24 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_26 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_28 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_30 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_32 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32770];
            const double xi_34 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32770];
            const double xi_36 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769];
            const double xi_38 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769];
            const double xi_40 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
            const double xi_42 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
            const double xi_44 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538];
            const double xi_46 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538];
            const double xi_48 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537];
            const double xi_50 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537];
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767] = xi_4 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767] = xi_49 + xi_6;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768] = xi_49 + xi_8;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768] = xi_10 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769] = xi_12 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769] = xi_14 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32770] = xi_32 + xi_49;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32770] = xi_34 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769] = xi_36 + xi_37;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769] = xi_37 + xi_38;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768] = xi_40 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768] = xi_42 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538] = xi_44 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538] = xi_46 + xi_49;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537] = xi_48 + xi_49;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537] = xi_49 + xi_50;
            _data_vertexFineDst[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
         // diagonal edge
         for (int ctr_1 = 16384 - ctr_2; ctr_1 < 16385 - ctr_2; ctr_1 += 1)
         {
            const double xi_151 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767];
            const double xi_153 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_155 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769];
            const double xi_157 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537];
            const double xi_176 = -0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_159 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767];
            const double xi_161 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_163 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_170 = 0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_165 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_167 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32770];
            const double xi_169 = _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32770];
            const double xi_171 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769];
            const double xi_173 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538];
            const double xi_175 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538];
            const double xi_177 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537];
            const double xi_156 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_154 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767] = xi_151 + xi_156;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_153 + xi_154;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769] = xi_154 + xi_155;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537] = xi_156 + xi_157;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767] = xi_159 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_161 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_163 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_165 + xi_170;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32770] = xi_167 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32770] = xi_169 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769] = xi_170 + xi_171;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538] = xi_173 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538] = xi_175 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537] = xi_176 + xi_177;
            _data_vertexFineDst[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_74*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
      for (int ctr_2 = 16384; ctr_2 < 16385; ctr_2 += 1)
      {
         // top vertex
         for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
         {
            const double xi_97 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769];
            const double xi_99 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537];
            const double xi_101 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769];
            const double xi_103 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538];
            const double xi_106 = -0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_105 = _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538];
            const double xi_107 = _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537];
            const double xi_96 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_98 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_100 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_102 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769] = xi_96 + xi_97;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537] = xi_98 + xi_99;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769] = xi_100 + xi_101;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538] = xi_102 + xi_103;
            _data_edgeFineDst_XY[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538] = xi_105 + xi_106;
            _data_edgeFineDst_Y[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537] = xi_106 + xi_107;
            _data_vertexFineDst[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_94*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         }
      }
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_any(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, int64_t coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
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
            const double xi_184 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_186 = _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_188 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_190 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_209 = -0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_192 = _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
            const double xi_199 = 0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_194 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
            const double xi_196 = _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_198 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
            const double xi_200 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_202 = _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_204 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_206 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_208 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_210 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_189 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_187 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_184 + xi_189;
            _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_186 + xi_187;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_187 + xi_188;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))] = xi_189 + xi_190;
            _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_192 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_194 + xi_199;
            _data_edgeFineDst_Y[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_196 + xi_209;
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_198 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_199 + xi_200;
            _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_202 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_204 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))] = xi_206 + xi_209;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1] = xi_208 + xi_209;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1] = xi_209 + xi_210;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_55*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
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
            const double xi_151 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_153 = _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_155 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_157 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_176 = -0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_159 = _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
            const double xi_161 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_163 = _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
            const double xi_170 = 0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_165 = _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
            const double xi_167 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1];
            const double xi_169 = _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1];
            const double xi_171 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))];
            const double xi_173 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_175 = _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))];
            const double xi_177 = _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1];
            const double xi_156 = xi_74*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            const double xi_154 = xi_74*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2] = xi_151 + xi_156;
            _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_153 + xi_154;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_154 + xi_155;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1] = xi_156 + xi_157;
            _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2] = xi_159 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_161 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_163 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_165 + xi_170;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1] = xi_167 + xi_176;
            _data_edgeFineDst_X[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1] = xi_169 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2))] = xi_170 + xi_171;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))] = xi_173 + xi_176;
            _data_edgeFineDst_XY[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2))] = xi_175 + xi_176;
            _data_edgeFineDst_Y[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 1] = xi_176 + xi_177;
            _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_74*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
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


void prolongate_2D_macroface_P2_push_from_vertexdofs(double * RESTRICT _data_edgeFineDst_X, double * RESTRICT _data_edgeFineDst_XY, double * RESTRICT _data_edgeFineDst_Y, double * RESTRICT _data_vertexCoarseSrc, double * RESTRICT _data_vertexFineDst, int64_t coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
    switch( coarse_level )
    {
    case 2:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_2(_data_edgeFineDst_X, _data_edgeFineDst_XY, _data_edgeFineDst_Y, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 3:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_3(_data_edgeFineDst_X, _data_edgeFineDst_XY, _data_edgeFineDst_Y, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 4:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_4(_data_edgeFineDst_X, _data_edgeFineDst_XY, _data_edgeFineDst_Y, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 5:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_5(_data_edgeFineDst_X, _data_edgeFineDst_XY, _data_edgeFineDst_Y, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 6:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_6(_data_edgeFineDst_X, _data_edgeFineDst_XY, _data_edgeFineDst_Y, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 7:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_7(_data_edgeFineDst_X, _data_edgeFineDst_XY, _data_edgeFineDst_Y, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 8:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_8(_data_edgeFineDst_X, _data_edgeFineDst_XY, _data_edgeFineDst_Y, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 9:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_9(_data_edgeFineDst_X, _data_edgeFineDst_XY, _data_edgeFineDst_Y, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 10:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_10(_data_edgeFineDst_X, _data_edgeFineDst_XY, _data_edgeFineDst_Y, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 11:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_11(_data_edgeFineDst_X, _data_edgeFineDst_XY, _data_edgeFineDst_Y, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 12:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_12(_data_edgeFineDst_X, _data_edgeFineDst_XY, _data_edgeFineDst_Y, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 13:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_13(_data_edgeFineDst_X, _data_edgeFineDst_XY, _data_edgeFineDst_Y, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 14:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_14(_data_edgeFineDst_X, _data_edgeFineDst_XY, _data_edgeFineDst_Y, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    default:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_any(_data_edgeFineDst_X, _data_edgeFineDst_XY, _data_edgeFineDst_Y, _data_vertexCoarseSrc, _data_vertexFineDst, coarse_level, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hhg
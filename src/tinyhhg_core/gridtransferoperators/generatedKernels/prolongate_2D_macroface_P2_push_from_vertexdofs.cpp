
//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#include "GeneratedKernelsP2MacroFace2D.hpp"

namespace hhg {
namespace P2 {
namespace macroface {
namespace generated {

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_2(double * _data_edgeFineDst, double * _data_vertexCoarseSrc, double * _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_55 = 1 / (num_neighbor_faces_edge0);
   const double xi_57 = 1 / (num_neighbor_faces_edge2);
   const double xi_58 = 1 / (num_neighbor_faces_vertex0);
   const double xi_125 = 1 / (num_neighbor_faces_edge0);
   const double xi_76 = 1 / (num_neighbor_faces_edge0);
   const double xi_78 = 1 / (num_neighbor_faces_edge1);
   const double xi_79 = 1 / (num_neighbor_faces_vertex1);
   const double xi_199 = 1 / (num_neighbor_faces_edge2);
   const double xi_162 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_99 = 1 / (num_neighbor_faces_edge2);
   const double xi_100 = 1 / (num_neighbor_faces_vertex2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_63 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_65 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_67 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_70 = -0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_69 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_71 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_60 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_57*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = xi_57*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_60 + xi_61;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_62 + xi_63;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9] = xi_64 + xi_65;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_66 + xi_67;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9] = xi_69 + xi_70;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_70 + xi_71;
         _data_vertexFineDst[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_58*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 4; ctr_1 += 1)
      {
         const double xi_128 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_130 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_134 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_153 = -0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7];
         const double xi_138 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7];
         const double xi_140 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
         const double xi_142 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
         const double xi_144 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_146 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_148 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_151 = 0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_150 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_152 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_154 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_133 = xi_125*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_131 = xi_125*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_128 + xi_133;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_130 + xi_131;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_131 + xi_132;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_133 + xi_134;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7] = xi_136 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7] = xi_138 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8] = xi_140 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8] = xi_142 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9] = xi_144 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9] = xi_146 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_148 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_150 + xi_151;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_151 + xi_152;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_153 + xi_154;
         _data_vertexFineDst[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 4; ctr_1 < 5; ctr_1 += 1)
      {
         const double xi_82 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_84 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_86 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7];
         const double xi_88 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_91 = -0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_90 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7];
         const double xi_92 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_81 = xi_76*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = xi_76*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_85 = xi_78*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = xi_78*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_81 + xi_82;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7] = xi_85 + xi_86;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_87 + xi_88;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7] = xi_90 + xi_91;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_91 + xi_92;
         _data_vertexFineDst[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_79*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 4; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_202 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_204 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_206 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9];
         const double xi_208 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18];
         const double xi_227 = -0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_210 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_217 = 0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_212 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_216 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_218 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9];
         const double xi_220 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
         const double xi_222 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
         const double xi_224 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18];
         const double xi_226 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17];
         const double xi_228 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17];
         const double xi_207 = xi_199*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_205 = xi_199*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9] = xi_202 + xi_207;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_204 + xi_205;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9] = xi_205 + xi_206;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18] = xi_207 + xi_208;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9] = xi_210 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_212 + xi_217;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_214 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_216 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9] = xi_217 + xi_218;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8] = xi_220 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8] = xi_222 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18] = xi_224 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17] = xi_226 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17] = xi_227 + xi_228;
         _data_vertexFineDst[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_199*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4; ctr_1 += 1)
      {
         const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7];
         const double xi_6 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7];
         const double xi_8 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
         const double xi_10 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8];
         const double xi_12 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_14 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9];
         const double xi_16 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_18 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_22 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_26 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_32 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 10];
         const double xi_34 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 10];
         const double xi_36 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9];
         const double xi_38 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9];
         const double xi_40 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
         const double xi_42 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8];
         const double xi_44 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18];
         const double xi_46 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18];
         const double xi_48 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17];
         const double xi_50 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17];
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7] = xi_4 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7] = xi_49 + xi_6;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8] = xi_49 + xi_8;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8] = xi_10 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9] = xi_12 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 9] = xi_14 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 10] = xi_32 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 10] = xi_34 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9] = xi_36 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9] = xi_37 + xi_38;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8] = xi_40 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8] = xi_42 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18] = xi_44 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18] = xi_46 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17] = xi_48 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17] = xi_49 + xi_50;
         _data_vertexFineDst[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 4; ctr_1 < -ctr_2 + 5; ctr_1 += 1)
      {
         const double xi_165 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7];
         const double xi_167 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_169 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9];
         const double xi_171 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17];
         const double xi_190 = -0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7];
         const double xi_175 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_177 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_184 = 0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_179 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_181 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 10];
         const double xi_183 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 10];
         const double xi_185 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9];
         const double xi_187 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18];
         const double xi_189 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18];
         const double xi_191 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17];
         const double xi_170 = xi_162*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_168 = xi_162*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7] = xi_165 + xi_170;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_167 + xi_168;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9] = xi_168 + xi_169;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17] = xi_170 + xi_171;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 7] = xi_173 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_175 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_177 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_179 + xi_184;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 10] = xi_181 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 10] = xi_183 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9] = xi_184 + xi_185;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18] = xi_187 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18] = xi_189 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17] = xi_190 + xi_191;
         _data_vertexFineDst[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_162*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 4; ctr_2 < 5; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_103 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9];
         const double xi_105 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17];
         const double xi_107 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9];
         const double xi_109 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18];
         const double xi_112 = -0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18];
         const double xi_113 = _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17];
         const double xi_102 = xi_97*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_104 = xi_97*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_106 = xi_99*0.375*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = xi_99*-0.125*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9] = xi_102 + xi_103;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17] = xi_104 + xi_105;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 9] = xi_106 + xi_107;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18] = xi_108 + xi_109;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + ((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 18] = xi_111 + xi_112;
         _data_edgeFineDst[2*ctr_1 + 18*ctr_2 + 2*((72) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 17] = xi_112 + xi_113;
         _data_vertexFineDst[2*ctr_1 + 20*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_100*_data_vertexCoarseSrc[ctr_1 + 6*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_3(double * _data_edgeFineDst, double * _data_vertexCoarseSrc, double * _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_55 = 1 / (num_neighbor_faces_edge0);
   const double xi_57 = 1 / (num_neighbor_faces_edge2);
   const double xi_58 = 1 / (num_neighbor_faces_vertex0);
   const double xi_125 = 1 / (num_neighbor_faces_edge0);
   const double xi_76 = 1 / (num_neighbor_faces_edge0);
   const double xi_78 = 1 / (num_neighbor_faces_edge1);
   const double xi_79 = 1 / (num_neighbor_faces_vertex1);
   const double xi_199 = 1 / (num_neighbor_faces_edge2);
   const double xi_162 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_99 = 1 / (num_neighbor_faces_edge2);
   const double xi_100 = 1 / (num_neighbor_faces_vertex2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_63 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_65 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_67 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_70 = -0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_69 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_71 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_60 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_57*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = xi_57*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_60 + xi_61;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_62 + xi_63;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17] = xi_64 + xi_65;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_66 + xi_67;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17] = xi_69 + xi_70;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_70 + xi_71;
         _data_vertexFineDst[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_58*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 8; ctr_1 += 1)
      {
         const double xi_128 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_130 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_134 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_153 = -0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15];
         const double xi_138 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15];
         const double xi_140 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
         const double xi_142 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
         const double xi_144 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_146 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_148 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_151 = 0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_150 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_152 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_154 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_133 = xi_125*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_131 = xi_125*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_128 + xi_133;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_130 + xi_131;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_131 + xi_132;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_133 + xi_134;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15] = xi_136 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15] = xi_138 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16] = xi_140 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16] = xi_142 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17] = xi_144 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17] = xi_146 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_148 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_150 + xi_151;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_151 + xi_152;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_153 + xi_154;
         _data_vertexFineDst[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 8; ctr_1 < 9; ctr_1 += 1)
      {
         const double xi_82 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_84 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_86 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15];
         const double xi_88 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_91 = -0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_90 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15];
         const double xi_92 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_81 = xi_76*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = xi_76*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_85 = xi_78*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = xi_78*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_81 + xi_82;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15] = xi_85 + xi_86;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_87 + xi_88;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15] = xi_90 + xi_91;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_91 + xi_92;
         _data_vertexFineDst[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_79*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 8; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_202 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_204 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_206 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17];
         const double xi_208 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34];
         const double xi_227 = -0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_210 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_217 = 0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_212 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_216 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_218 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17];
         const double xi_220 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
         const double xi_222 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
         const double xi_224 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34];
         const double xi_226 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33];
         const double xi_228 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33];
         const double xi_207 = xi_199*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_205 = xi_199*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17] = xi_202 + xi_207;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_204 + xi_205;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17] = xi_205 + xi_206;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34] = xi_207 + xi_208;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17] = xi_210 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_212 + xi_217;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_214 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_216 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17] = xi_217 + xi_218;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16] = xi_220 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16] = xi_222 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34] = xi_224 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33] = xi_226 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33] = xi_227 + xi_228;
         _data_vertexFineDst[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_199*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8; ctr_1 += 1)
      {
         const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15];
         const double xi_6 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15];
         const double xi_8 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
         const double xi_10 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16];
         const double xi_12 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_14 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17];
         const double xi_16 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_18 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_22 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_26 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_32 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 18];
         const double xi_34 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 18];
         const double xi_36 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17];
         const double xi_38 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17];
         const double xi_40 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
         const double xi_42 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16];
         const double xi_44 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34];
         const double xi_46 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34];
         const double xi_48 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33];
         const double xi_50 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33];
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15] = xi_4 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15] = xi_49 + xi_6;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16] = xi_49 + xi_8;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16] = xi_10 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17] = xi_12 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 17] = xi_14 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 18] = xi_32 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 18] = xi_34 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17] = xi_36 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17] = xi_37 + xi_38;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16] = xi_40 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16] = xi_42 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34] = xi_44 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34] = xi_46 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33] = xi_48 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33] = xi_49 + xi_50;
         _data_vertexFineDst[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 8; ctr_1 < -ctr_2 + 9; ctr_1 += 1)
      {
         const double xi_165 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15];
         const double xi_167 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_169 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17];
         const double xi_171 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33];
         const double xi_190 = -0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15];
         const double xi_175 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_177 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_184 = 0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_179 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_181 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 18];
         const double xi_183 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 18];
         const double xi_185 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17];
         const double xi_187 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34];
         const double xi_189 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34];
         const double xi_191 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33];
         const double xi_170 = xi_162*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_168 = xi_162*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15] = xi_165 + xi_170;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_167 + xi_168;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17] = xi_168 + xi_169;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33] = xi_170 + xi_171;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 15] = xi_173 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_175 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_177 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_179 + xi_184;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 18] = xi_181 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 18] = xi_183 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17] = xi_184 + xi_185;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34] = xi_187 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34] = xi_189 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33] = xi_190 + xi_191;
         _data_vertexFineDst[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_162*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 8; ctr_2 < 9; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_103 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17];
         const double xi_105 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33];
         const double xi_107 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17];
         const double xi_109 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34];
         const double xi_112 = -0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34];
         const double xi_113 = _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33];
         const double xi_102 = xi_97*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_104 = xi_97*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_106 = xi_99*0.375*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = xi_99*-0.125*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17] = xi_102 + xi_103;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33] = xi_104 + xi_105;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 17] = xi_106 + xi_107;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34] = xi_108 + xi_109;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + ((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 34] = xi_111 + xi_112;
         _data_edgeFineDst[2*ctr_1 + 34*ctr_2 + 2*((272) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 33] = xi_112 + xi_113;
         _data_vertexFineDst[2*ctr_1 + 36*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_100*_data_vertexCoarseSrc[ctr_1 + 10*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_4(double * _data_edgeFineDst, double * _data_vertexCoarseSrc, double * _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_55 = 1 / (num_neighbor_faces_edge0);
   const double xi_57 = 1 / (num_neighbor_faces_edge2);
   const double xi_58 = 1 / (num_neighbor_faces_vertex0);
   const double xi_125 = 1 / (num_neighbor_faces_edge0);
   const double xi_76 = 1 / (num_neighbor_faces_edge0);
   const double xi_78 = 1 / (num_neighbor_faces_edge1);
   const double xi_79 = 1 / (num_neighbor_faces_vertex1);
   const double xi_199 = 1 / (num_neighbor_faces_edge2);
   const double xi_162 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_99 = 1 / (num_neighbor_faces_edge2);
   const double xi_100 = 1 / (num_neighbor_faces_vertex2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_63 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_65 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_67 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_70 = -0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_69 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_71 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_60 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_57*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = xi_57*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_60 + xi_61;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_62 + xi_63;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33] = xi_64 + xi_65;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_66 + xi_67;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33] = xi_69 + xi_70;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_70 + xi_71;
         _data_vertexFineDst[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_58*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 16; ctr_1 += 1)
      {
         const double xi_128 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_130 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_134 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_153 = -0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31];
         const double xi_138 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31];
         const double xi_140 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
         const double xi_142 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
         const double xi_144 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_146 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_148 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_151 = 0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_150 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_152 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_154 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_133 = xi_125*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_131 = xi_125*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_128 + xi_133;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_130 + xi_131;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_131 + xi_132;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_133 + xi_134;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31] = xi_136 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31] = xi_138 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32] = xi_140 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32] = xi_142 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33] = xi_144 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33] = xi_146 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_148 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_150 + xi_151;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_151 + xi_152;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_153 + xi_154;
         _data_vertexFineDst[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 16; ctr_1 < 17; ctr_1 += 1)
      {
         const double xi_82 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_84 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_86 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31];
         const double xi_88 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_91 = -0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_90 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31];
         const double xi_92 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_81 = xi_76*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = xi_76*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_85 = xi_78*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = xi_78*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_81 + xi_82;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31] = xi_85 + xi_86;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_87 + xi_88;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31] = xi_90 + xi_91;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_91 + xi_92;
         _data_vertexFineDst[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_79*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 16; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_202 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_204 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_206 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33];
         const double xi_208 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66];
         const double xi_227 = -0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_210 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_217 = 0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_212 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_216 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_218 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33];
         const double xi_220 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
         const double xi_222 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
         const double xi_224 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66];
         const double xi_226 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65];
         const double xi_228 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65];
         const double xi_207 = xi_199*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_205 = xi_199*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33] = xi_202 + xi_207;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_204 + xi_205;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33] = xi_205 + xi_206;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66] = xi_207 + xi_208;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33] = xi_210 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_212 + xi_217;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_214 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_216 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33] = xi_217 + xi_218;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32] = xi_220 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32] = xi_222 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66] = xi_224 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65] = xi_226 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65] = xi_227 + xi_228;
         _data_vertexFineDst[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_199*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16; ctr_1 += 1)
      {
         const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31];
         const double xi_6 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31];
         const double xi_8 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
         const double xi_10 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32];
         const double xi_12 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_14 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33];
         const double xi_16 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_18 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_22 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_26 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_32 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 34];
         const double xi_34 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 34];
         const double xi_36 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33];
         const double xi_38 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33];
         const double xi_40 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
         const double xi_42 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32];
         const double xi_44 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66];
         const double xi_46 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66];
         const double xi_48 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65];
         const double xi_50 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65];
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31] = xi_4 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31] = xi_49 + xi_6;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32] = xi_49 + xi_8;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32] = xi_10 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33] = xi_12 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 33] = xi_14 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 34] = xi_32 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 34] = xi_34 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33] = xi_36 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33] = xi_37 + xi_38;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32] = xi_40 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32] = xi_42 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66] = xi_44 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66] = xi_46 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65] = xi_48 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65] = xi_49 + xi_50;
         _data_vertexFineDst[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 16; ctr_1 < -ctr_2 + 17; ctr_1 += 1)
      {
         const double xi_165 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31];
         const double xi_167 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_169 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33];
         const double xi_171 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65];
         const double xi_190 = -0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31];
         const double xi_175 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_177 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_184 = 0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_179 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_181 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 34];
         const double xi_183 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 34];
         const double xi_185 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33];
         const double xi_187 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66];
         const double xi_189 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66];
         const double xi_191 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65];
         const double xi_170 = xi_162*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_168 = xi_162*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31] = xi_165 + xi_170;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_167 + xi_168;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33] = xi_168 + xi_169;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65] = xi_170 + xi_171;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 31] = xi_173 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_175 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_177 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_179 + xi_184;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 34] = xi_181 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 34] = xi_183 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33] = xi_184 + xi_185;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66] = xi_187 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66] = xi_189 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65] = xi_190 + xi_191;
         _data_vertexFineDst[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_162*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 16; ctr_2 < 17; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_103 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33];
         const double xi_105 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65];
         const double xi_107 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33];
         const double xi_109 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66];
         const double xi_112 = -0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66];
         const double xi_113 = _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65];
         const double xi_102 = xi_97*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_104 = xi_97*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_106 = xi_99*0.375*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = xi_99*-0.125*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33] = xi_102 + xi_103;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65] = xi_104 + xi_105;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 33] = xi_106 + xi_107;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66] = xi_108 + xi_109;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + ((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 66] = xi_111 + xi_112;
         _data_edgeFineDst[2*ctr_1 + 66*ctr_2 + 2*((1056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65] = xi_112 + xi_113;
         _data_vertexFineDst[2*ctr_1 + 68*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_100*_data_vertexCoarseSrc[ctr_1 + 18*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_5(double * _data_edgeFineDst, double * _data_vertexCoarseSrc, double * _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_55 = 1 / (num_neighbor_faces_edge0);
   const double xi_57 = 1 / (num_neighbor_faces_edge2);
   const double xi_58 = 1 / (num_neighbor_faces_vertex0);
   const double xi_125 = 1 / (num_neighbor_faces_edge0);
   const double xi_76 = 1 / (num_neighbor_faces_edge0);
   const double xi_78 = 1 / (num_neighbor_faces_edge1);
   const double xi_79 = 1 / (num_neighbor_faces_vertex1);
   const double xi_199 = 1 / (num_neighbor_faces_edge2);
   const double xi_162 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_99 = 1 / (num_neighbor_faces_edge2);
   const double xi_100 = 1 / (num_neighbor_faces_vertex2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_63 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_65 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_67 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_70 = -0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_69 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_71 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_60 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_57*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = xi_57*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_60 + xi_61;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_62 + xi_63;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65] = xi_64 + xi_65;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_66 + xi_67;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65] = xi_69 + xi_70;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_70 + xi_71;
         _data_vertexFineDst[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_58*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 32; ctr_1 += 1)
      {
         const double xi_128 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_130 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_134 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_153 = -0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63];
         const double xi_138 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63];
         const double xi_140 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
         const double xi_142 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
         const double xi_144 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_146 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_148 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_151 = 0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_150 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_152 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_154 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_133 = xi_125*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_131 = xi_125*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_128 + xi_133;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_130 + xi_131;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_131 + xi_132;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_133 + xi_134;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63] = xi_136 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63] = xi_138 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64] = xi_140 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64] = xi_142 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65] = xi_144 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65] = xi_146 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_148 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_150 + xi_151;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_151 + xi_152;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_153 + xi_154;
         _data_vertexFineDst[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 32; ctr_1 < 33; ctr_1 += 1)
      {
         const double xi_82 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_84 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_86 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63];
         const double xi_88 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_91 = -0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_90 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63];
         const double xi_92 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_81 = xi_76*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = xi_76*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_85 = xi_78*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = xi_78*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_81 + xi_82;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63] = xi_85 + xi_86;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_87 + xi_88;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63] = xi_90 + xi_91;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_91 + xi_92;
         _data_vertexFineDst[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_79*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 32; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_202 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_204 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_206 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65];
         const double xi_208 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130];
         const double xi_227 = -0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_210 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_217 = 0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_212 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_216 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_218 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65];
         const double xi_220 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
         const double xi_222 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
         const double xi_224 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130];
         const double xi_226 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129];
         const double xi_228 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129];
         const double xi_207 = xi_199*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_205 = xi_199*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65] = xi_202 + xi_207;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_204 + xi_205;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65] = xi_205 + xi_206;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130] = xi_207 + xi_208;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65] = xi_210 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_212 + xi_217;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_214 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_216 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65] = xi_217 + xi_218;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64] = xi_220 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64] = xi_222 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130] = xi_224 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129] = xi_226 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129] = xi_227 + xi_228;
         _data_vertexFineDst[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_199*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 32; ctr_1 += 1)
      {
         const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63];
         const double xi_6 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63];
         const double xi_8 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
         const double xi_10 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64];
         const double xi_12 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_14 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65];
         const double xi_16 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_18 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_22 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_26 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_32 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 66];
         const double xi_34 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 66];
         const double xi_36 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65];
         const double xi_38 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65];
         const double xi_40 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
         const double xi_42 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64];
         const double xi_44 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130];
         const double xi_46 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130];
         const double xi_48 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129];
         const double xi_50 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129];
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63] = xi_4 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63] = xi_49 + xi_6;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64] = xi_49 + xi_8;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 64] = xi_10 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65] = xi_12 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 65] = xi_14 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 66] = xi_32 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 66] = xi_34 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65] = xi_36 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65] = xi_37 + xi_38;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64] = xi_40 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 64] = xi_42 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130] = xi_44 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130] = xi_46 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129] = xi_48 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129] = xi_49 + xi_50;
         _data_vertexFineDst[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 32; ctr_1 < -ctr_2 + 33; ctr_1 += 1)
      {
         const double xi_165 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63];
         const double xi_167 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_169 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65];
         const double xi_171 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129];
         const double xi_190 = -0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63];
         const double xi_175 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_177 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_184 = 0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_179 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_181 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 66];
         const double xi_183 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 66];
         const double xi_185 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65];
         const double xi_187 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130];
         const double xi_189 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130];
         const double xi_191 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129];
         const double xi_170 = xi_162*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_168 = xi_162*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63] = xi_165 + xi_170;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_167 + xi_168;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65] = xi_168 + xi_169;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129] = xi_170 + xi_171;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 63] = xi_173 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_175 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_177 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_179 + xi_184;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 66] = xi_181 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 66] = xi_183 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65] = xi_184 + xi_185;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130] = xi_187 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130] = xi_189 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129] = xi_190 + xi_191;
         _data_vertexFineDst[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_162*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 32; ctr_2 < 33; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_103 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65];
         const double xi_105 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129];
         const double xi_107 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65];
         const double xi_109 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130];
         const double xi_112 = -0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130];
         const double xi_113 = _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129];
         const double xi_102 = xi_97*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_104 = xi_97*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_106 = xi_99*0.375*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = xi_99*-0.125*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65] = xi_102 + xi_103;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129] = xi_104 + xi_105;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 65] = xi_106 + xi_107;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130] = xi_108 + xi_109;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + ((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 130] = xi_111 + xi_112;
         _data_edgeFineDst[2*ctr_1 + 130*ctr_2 + 2*((4160) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 129] = xi_112 + xi_113;
         _data_vertexFineDst[2*ctr_1 + 132*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_100*_data_vertexCoarseSrc[ctr_1 + 34*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_6(double * _data_edgeFineDst, double * _data_vertexCoarseSrc, double * _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_55 = 1 / (num_neighbor_faces_edge0);
   const double xi_57 = 1 / (num_neighbor_faces_edge2);
   const double xi_58 = 1 / (num_neighbor_faces_vertex0);
   const double xi_125 = 1 / (num_neighbor_faces_edge0);
   const double xi_76 = 1 / (num_neighbor_faces_edge0);
   const double xi_78 = 1 / (num_neighbor_faces_edge1);
   const double xi_79 = 1 / (num_neighbor_faces_vertex1);
   const double xi_199 = 1 / (num_neighbor_faces_edge2);
   const double xi_162 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_99 = 1 / (num_neighbor_faces_edge2);
   const double xi_100 = 1 / (num_neighbor_faces_vertex2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_63 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_65 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_67 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_70 = -0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_69 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_71 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_60 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_57*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = xi_57*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_60 + xi_61;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_62 + xi_63;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129] = xi_64 + xi_65;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_66 + xi_67;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129] = xi_69 + xi_70;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_70 + xi_71;
         _data_vertexFineDst[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_58*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 64; ctr_1 += 1)
      {
         const double xi_128 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_130 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_134 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_153 = -0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127];
         const double xi_138 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127];
         const double xi_140 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
         const double xi_142 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
         const double xi_144 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_146 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_148 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_151 = 0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_150 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_152 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_154 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_133 = xi_125*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_131 = xi_125*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_128 + xi_133;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_130 + xi_131;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_131 + xi_132;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_133 + xi_134;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127] = xi_136 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127] = xi_138 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128] = xi_140 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128] = xi_142 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129] = xi_144 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129] = xi_146 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_148 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_150 + xi_151;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_151 + xi_152;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_153 + xi_154;
         _data_vertexFineDst[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 64; ctr_1 < 65; ctr_1 += 1)
      {
         const double xi_82 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_84 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_86 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127];
         const double xi_88 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_91 = -0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_90 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127];
         const double xi_92 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_81 = xi_76*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = xi_76*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_85 = xi_78*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = xi_78*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_81 + xi_82;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127] = xi_85 + xi_86;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_87 + xi_88;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127] = xi_90 + xi_91;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_91 + xi_92;
         _data_vertexFineDst[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_79*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 64; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_202 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_204 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_206 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129];
         const double xi_208 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258];
         const double xi_227 = -0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_210 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_217 = 0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_212 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_216 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_218 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129];
         const double xi_220 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
         const double xi_222 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
         const double xi_224 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258];
         const double xi_226 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257];
         const double xi_228 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257];
         const double xi_207 = xi_199*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_205 = xi_199*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129] = xi_202 + xi_207;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_204 + xi_205;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129] = xi_205 + xi_206;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258] = xi_207 + xi_208;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129] = xi_210 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_212 + xi_217;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_214 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_216 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129] = xi_217 + xi_218;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128] = xi_220 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128] = xi_222 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258] = xi_224 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257] = xi_226 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257] = xi_227 + xi_228;
         _data_vertexFineDst[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_199*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 64; ctr_1 += 1)
      {
         const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127];
         const double xi_6 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127];
         const double xi_8 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
         const double xi_10 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128];
         const double xi_12 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_14 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129];
         const double xi_16 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_18 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_22 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_26 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_32 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 130];
         const double xi_34 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 130];
         const double xi_36 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129];
         const double xi_38 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129];
         const double xi_40 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
         const double xi_42 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128];
         const double xi_44 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258];
         const double xi_46 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258];
         const double xi_48 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257];
         const double xi_50 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257];
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127] = xi_4 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127] = xi_49 + xi_6;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128] = xi_49 + xi_8;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 128] = xi_10 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129] = xi_12 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 129] = xi_14 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 130] = xi_32 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 130] = xi_34 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129] = xi_36 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129] = xi_37 + xi_38;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128] = xi_40 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 128] = xi_42 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258] = xi_44 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258] = xi_46 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257] = xi_48 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257] = xi_49 + xi_50;
         _data_vertexFineDst[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 64; ctr_1 < -ctr_2 + 65; ctr_1 += 1)
      {
         const double xi_165 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127];
         const double xi_167 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_169 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129];
         const double xi_171 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257];
         const double xi_190 = -0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127];
         const double xi_175 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_177 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_184 = 0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_179 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_181 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 130];
         const double xi_183 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 130];
         const double xi_185 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129];
         const double xi_187 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258];
         const double xi_189 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258];
         const double xi_191 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257];
         const double xi_170 = xi_162*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_168 = xi_162*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127] = xi_165 + xi_170;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_167 + xi_168;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129] = xi_168 + xi_169;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257] = xi_170 + xi_171;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 127] = xi_173 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_175 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_177 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_179 + xi_184;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 130] = xi_181 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 130] = xi_183 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129] = xi_184 + xi_185;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258] = xi_187 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258] = xi_189 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257] = xi_190 + xi_191;
         _data_vertexFineDst[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_162*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 64; ctr_2 < 65; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_103 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129];
         const double xi_105 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257];
         const double xi_107 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129];
         const double xi_109 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258];
         const double xi_112 = -0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258];
         const double xi_113 = _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257];
         const double xi_102 = xi_97*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_104 = xi_97*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_106 = xi_99*0.375*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = xi_99*-0.125*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129] = xi_102 + xi_103;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257] = xi_104 + xi_105;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 129] = xi_106 + xi_107;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258] = xi_108 + xi_109;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + ((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 258] = xi_111 + xi_112;
         _data_edgeFineDst[2*ctr_1 + 258*ctr_2 + 2*((16512) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 257] = xi_112 + xi_113;
         _data_vertexFineDst[2*ctr_1 + 260*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_100*_data_vertexCoarseSrc[ctr_1 + 66*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_7(double * _data_edgeFineDst, double * _data_vertexCoarseSrc, double * _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_55 = 1 / (num_neighbor_faces_edge0);
   const double xi_57 = 1 / (num_neighbor_faces_edge2);
   const double xi_58 = 1 / (num_neighbor_faces_vertex0);
   const double xi_125 = 1 / (num_neighbor_faces_edge0);
   const double xi_76 = 1 / (num_neighbor_faces_edge0);
   const double xi_78 = 1 / (num_neighbor_faces_edge1);
   const double xi_79 = 1 / (num_neighbor_faces_vertex1);
   const double xi_199 = 1 / (num_neighbor_faces_edge2);
   const double xi_162 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_99 = 1 / (num_neighbor_faces_edge2);
   const double xi_100 = 1 / (num_neighbor_faces_vertex2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_63 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_65 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_67 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_70 = -0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_69 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_71 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_60 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_57*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = xi_57*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_60 + xi_61;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_62 + xi_63;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257] = xi_64 + xi_65;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_66 + xi_67;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257] = xi_69 + xi_70;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_70 + xi_71;
         _data_vertexFineDst[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_58*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 128; ctr_1 += 1)
      {
         const double xi_128 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_130 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_134 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_153 = -0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255];
         const double xi_138 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255];
         const double xi_140 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
         const double xi_142 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
         const double xi_144 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_146 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_148 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_151 = 0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_150 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_152 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_154 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_133 = xi_125*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_131 = xi_125*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_128 + xi_133;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_130 + xi_131;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_131 + xi_132;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_133 + xi_134;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255] = xi_136 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255] = xi_138 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256] = xi_140 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256] = xi_142 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257] = xi_144 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257] = xi_146 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_148 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_150 + xi_151;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_151 + xi_152;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_153 + xi_154;
         _data_vertexFineDst[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 128; ctr_1 < 129; ctr_1 += 1)
      {
         const double xi_82 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_84 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_86 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255];
         const double xi_88 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_91 = -0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_90 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255];
         const double xi_92 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_81 = xi_76*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = xi_76*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_85 = xi_78*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = xi_78*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_81 + xi_82;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255] = xi_85 + xi_86;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_87 + xi_88;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255] = xi_90 + xi_91;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_91 + xi_92;
         _data_vertexFineDst[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_79*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 128; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_202 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_204 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_206 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257];
         const double xi_208 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514];
         const double xi_227 = -0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_210 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_217 = 0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_212 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_216 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_218 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257];
         const double xi_220 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
         const double xi_222 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
         const double xi_224 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514];
         const double xi_226 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513];
         const double xi_228 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513];
         const double xi_207 = xi_199*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_205 = xi_199*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257] = xi_202 + xi_207;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_204 + xi_205;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257] = xi_205 + xi_206;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514] = xi_207 + xi_208;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257] = xi_210 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_212 + xi_217;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_214 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_216 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257] = xi_217 + xi_218;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256] = xi_220 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256] = xi_222 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514] = xi_224 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513] = xi_226 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513] = xi_227 + xi_228;
         _data_vertexFineDst[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_199*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 128; ctr_1 += 1)
      {
         const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255];
         const double xi_6 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255];
         const double xi_8 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
         const double xi_10 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256];
         const double xi_12 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_14 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257];
         const double xi_16 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_18 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_22 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_26 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_32 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 258];
         const double xi_34 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 258];
         const double xi_36 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257];
         const double xi_38 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257];
         const double xi_40 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
         const double xi_42 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256];
         const double xi_44 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514];
         const double xi_46 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514];
         const double xi_48 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513];
         const double xi_50 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513];
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255] = xi_4 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255] = xi_49 + xi_6;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256] = xi_49 + xi_8;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 256] = xi_10 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257] = xi_12 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 257] = xi_14 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 258] = xi_32 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 258] = xi_34 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257] = xi_36 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257] = xi_37 + xi_38;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256] = xi_40 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 256] = xi_42 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514] = xi_44 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514] = xi_46 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513] = xi_48 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513] = xi_49 + xi_50;
         _data_vertexFineDst[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 128; ctr_1 < -ctr_2 + 129; ctr_1 += 1)
      {
         const double xi_165 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255];
         const double xi_167 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_169 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257];
         const double xi_171 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513];
         const double xi_190 = -0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255];
         const double xi_175 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_177 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_184 = 0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_179 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_181 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 258];
         const double xi_183 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 258];
         const double xi_185 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257];
         const double xi_187 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514];
         const double xi_189 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514];
         const double xi_191 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513];
         const double xi_170 = xi_162*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_168 = xi_162*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255] = xi_165 + xi_170;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_167 + xi_168;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257] = xi_168 + xi_169;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513] = xi_170 + xi_171;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 255] = xi_173 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_175 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_177 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_179 + xi_184;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 258] = xi_181 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 258] = xi_183 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257] = xi_184 + xi_185;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514] = xi_187 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514] = xi_189 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513] = xi_190 + xi_191;
         _data_vertexFineDst[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_162*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 128; ctr_2 < 129; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_103 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257];
         const double xi_105 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513];
         const double xi_107 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257];
         const double xi_109 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514];
         const double xi_112 = -0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514];
         const double xi_113 = _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513];
         const double xi_102 = xi_97*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_104 = xi_97*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_106 = xi_99*0.375*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = xi_99*-0.125*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257] = xi_102 + xi_103;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513] = xi_104 + xi_105;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 257] = xi_106 + xi_107;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514] = xi_108 + xi_109;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + ((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 514] = xi_111 + xi_112;
         _data_edgeFineDst[2*ctr_1 + 514*ctr_2 + 2*((65792) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 513] = xi_112 + xi_113;
         _data_vertexFineDst[2*ctr_1 + 516*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_100*_data_vertexCoarseSrc[ctr_1 + 130*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_8(double * _data_edgeFineDst, double * _data_vertexCoarseSrc, double * _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_55 = 1 / (num_neighbor_faces_edge0);
   const double xi_57 = 1 / (num_neighbor_faces_edge2);
   const double xi_58 = 1 / (num_neighbor_faces_vertex0);
   const double xi_125 = 1 / (num_neighbor_faces_edge0);
   const double xi_76 = 1 / (num_neighbor_faces_edge0);
   const double xi_78 = 1 / (num_neighbor_faces_edge1);
   const double xi_79 = 1 / (num_neighbor_faces_vertex1);
   const double xi_199 = 1 / (num_neighbor_faces_edge2);
   const double xi_162 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_99 = 1 / (num_neighbor_faces_edge2);
   const double xi_100 = 1 / (num_neighbor_faces_vertex2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_63 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_65 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_67 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_70 = -0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_69 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_71 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_60 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_57*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = xi_57*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_60 + xi_61;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_62 + xi_63;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513] = xi_64 + xi_65;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_66 + xi_67;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513] = xi_69 + xi_70;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_70 + xi_71;
         _data_vertexFineDst[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_58*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 256; ctr_1 += 1)
      {
         const double xi_128 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_130 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_134 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_153 = -0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511];
         const double xi_138 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511];
         const double xi_140 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
         const double xi_142 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
         const double xi_144 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_146 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_148 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_151 = 0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_150 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_152 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_154 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_133 = xi_125*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_131 = xi_125*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_128 + xi_133;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_130 + xi_131;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_131 + xi_132;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_133 + xi_134;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511] = xi_136 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511] = xi_138 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512] = xi_140 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512] = xi_142 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513] = xi_144 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513] = xi_146 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_148 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_150 + xi_151;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_151 + xi_152;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_153 + xi_154;
         _data_vertexFineDst[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 256; ctr_1 < 257; ctr_1 += 1)
      {
         const double xi_82 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_84 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_86 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511];
         const double xi_88 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_91 = -0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_90 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511];
         const double xi_92 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_81 = xi_76*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = xi_76*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_85 = xi_78*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = xi_78*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_81 + xi_82;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511] = xi_85 + xi_86;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_87 + xi_88;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511] = xi_90 + xi_91;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_91 + xi_92;
         _data_vertexFineDst[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_79*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 256; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_202 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_204 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_206 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513];
         const double xi_208 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026];
         const double xi_227 = -0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_210 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_217 = 0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_212 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_216 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_218 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513];
         const double xi_220 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
         const double xi_222 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
         const double xi_224 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026];
         const double xi_226 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025];
         const double xi_228 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025];
         const double xi_207 = xi_199*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_205 = xi_199*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513] = xi_202 + xi_207;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_204 + xi_205;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513] = xi_205 + xi_206;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026] = xi_207 + xi_208;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513] = xi_210 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_212 + xi_217;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_214 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_216 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513] = xi_217 + xi_218;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512] = xi_220 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512] = xi_222 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026] = xi_224 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025] = xi_226 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025] = xi_227 + xi_228;
         _data_vertexFineDst[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_199*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 256; ctr_1 += 1)
      {
         const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511];
         const double xi_6 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511];
         const double xi_8 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
         const double xi_10 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512];
         const double xi_12 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_14 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513];
         const double xi_16 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_18 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_22 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_26 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_32 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 514];
         const double xi_34 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 514];
         const double xi_36 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513];
         const double xi_38 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513];
         const double xi_40 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
         const double xi_42 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512];
         const double xi_44 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026];
         const double xi_46 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026];
         const double xi_48 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025];
         const double xi_50 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025];
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511] = xi_4 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511] = xi_49 + xi_6;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512] = xi_49 + xi_8;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 512] = xi_10 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513] = xi_12 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 513] = xi_14 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 514] = xi_32 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 514] = xi_34 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513] = xi_36 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513] = xi_37 + xi_38;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512] = xi_40 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 512] = xi_42 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026] = xi_44 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026] = xi_46 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025] = xi_48 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025] = xi_49 + xi_50;
         _data_vertexFineDst[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 256; ctr_1 < -ctr_2 + 257; ctr_1 += 1)
      {
         const double xi_165 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511];
         const double xi_167 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_169 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513];
         const double xi_171 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025];
         const double xi_190 = -0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511];
         const double xi_175 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_177 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_184 = 0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_179 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_181 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 514];
         const double xi_183 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 514];
         const double xi_185 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513];
         const double xi_187 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026];
         const double xi_189 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026];
         const double xi_191 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025];
         const double xi_170 = xi_162*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_168 = xi_162*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511] = xi_165 + xi_170;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_167 + xi_168;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513] = xi_168 + xi_169;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025] = xi_170 + xi_171;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 511] = xi_173 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_175 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_177 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_179 + xi_184;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 514] = xi_181 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 514] = xi_183 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513] = xi_184 + xi_185;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026] = xi_187 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026] = xi_189 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025] = xi_190 + xi_191;
         _data_vertexFineDst[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_162*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 256; ctr_2 < 257; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_103 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513];
         const double xi_105 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025];
         const double xi_107 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513];
         const double xi_109 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026];
         const double xi_112 = -0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026];
         const double xi_113 = _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025];
         const double xi_102 = xi_97*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_104 = xi_97*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_106 = xi_99*0.375*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = xi_99*-0.125*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513] = xi_102 + xi_103;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025] = xi_104 + xi_105;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 513] = xi_106 + xi_107;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026] = xi_108 + xi_109;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + ((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1026] = xi_111 + xi_112;
         _data_edgeFineDst[2*ctr_1 + 1026*ctr_2 + 2*((262656) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 1025] = xi_112 + xi_113;
         _data_vertexFineDst[2*ctr_1 + 1028*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_100*_data_vertexCoarseSrc[ctr_1 + 258*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_9(double * _data_edgeFineDst, double * _data_vertexCoarseSrc, double * _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_55 = 1 / (num_neighbor_faces_edge0);
   const double xi_57 = 1 / (num_neighbor_faces_edge2);
   const double xi_58 = 1 / (num_neighbor_faces_vertex0);
   const double xi_125 = 1 / (num_neighbor_faces_edge0);
   const double xi_76 = 1 / (num_neighbor_faces_edge0);
   const double xi_78 = 1 / (num_neighbor_faces_edge1);
   const double xi_79 = 1 / (num_neighbor_faces_vertex1);
   const double xi_199 = 1 / (num_neighbor_faces_edge2);
   const double xi_162 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_99 = 1 / (num_neighbor_faces_edge2);
   const double xi_100 = 1 / (num_neighbor_faces_vertex2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_63 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_65 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_67 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_70 = -0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_69 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_71 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_60 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_57*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = xi_57*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_60 + xi_61;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_62 + xi_63;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025] = xi_64 + xi_65;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_66 + xi_67;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025] = xi_69 + xi_70;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_70 + xi_71;
         _data_vertexFineDst[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_58*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 512; ctr_1 += 1)
      {
         const double xi_128 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_130 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_134 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_153 = -0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023];
         const double xi_138 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023];
         const double xi_140 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
         const double xi_142 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
         const double xi_144 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_146 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_148 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_151 = 0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_150 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_152 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_154 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_133 = xi_125*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_131 = xi_125*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_128 + xi_133;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_130 + xi_131;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_131 + xi_132;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_133 + xi_134;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023] = xi_136 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023] = xi_138 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024] = xi_140 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024] = xi_142 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025] = xi_144 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025] = xi_146 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_148 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_150 + xi_151;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_151 + xi_152;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_153 + xi_154;
         _data_vertexFineDst[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 512; ctr_1 < 513; ctr_1 += 1)
      {
         const double xi_82 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_84 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_86 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023];
         const double xi_88 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_91 = -0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_90 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023];
         const double xi_92 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_81 = xi_76*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = xi_76*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_85 = xi_78*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = xi_78*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_81 + xi_82;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023] = xi_85 + xi_86;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_87 + xi_88;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023] = xi_90 + xi_91;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_91 + xi_92;
         _data_vertexFineDst[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_79*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 512; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_202 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_204 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_206 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025];
         const double xi_208 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050];
         const double xi_227 = -0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_210 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_217 = 0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_212 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_216 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_218 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025];
         const double xi_220 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
         const double xi_222 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
         const double xi_224 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050];
         const double xi_226 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049];
         const double xi_228 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049];
         const double xi_207 = xi_199*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_205 = xi_199*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025] = xi_202 + xi_207;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_204 + xi_205;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025] = xi_205 + xi_206;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050] = xi_207 + xi_208;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025] = xi_210 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_212 + xi_217;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_214 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_216 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025] = xi_217 + xi_218;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024] = xi_220 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024] = xi_222 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050] = xi_224 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049] = xi_226 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049] = xi_227 + xi_228;
         _data_vertexFineDst[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_199*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 512; ctr_1 += 1)
      {
         const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023];
         const double xi_6 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023];
         const double xi_8 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
         const double xi_10 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024];
         const double xi_12 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_14 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025];
         const double xi_16 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_18 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_22 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_26 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_32 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1026];
         const double xi_34 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1026];
         const double xi_36 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025];
         const double xi_38 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025];
         const double xi_40 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
         const double xi_42 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024];
         const double xi_44 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050];
         const double xi_46 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050];
         const double xi_48 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049];
         const double xi_50 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049];
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023] = xi_4 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023] = xi_49 + xi_6;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024] = xi_49 + xi_8;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1024] = xi_10 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025] = xi_12 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1025] = xi_14 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1026] = xi_32 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1026] = xi_34 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025] = xi_36 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025] = xi_37 + xi_38;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024] = xi_40 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1024] = xi_42 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050] = xi_44 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050] = xi_46 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049] = xi_48 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049] = xi_49 + xi_50;
         _data_vertexFineDst[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 512; ctr_1 < -ctr_2 + 513; ctr_1 += 1)
      {
         const double xi_165 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023];
         const double xi_167 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_169 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025];
         const double xi_171 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049];
         const double xi_190 = -0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023];
         const double xi_175 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_177 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_184 = 0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_179 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_181 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1026];
         const double xi_183 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1026];
         const double xi_185 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025];
         const double xi_187 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050];
         const double xi_189 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050];
         const double xi_191 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049];
         const double xi_170 = xi_162*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_168 = xi_162*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023] = xi_165 + xi_170;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_167 + xi_168;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025] = xi_168 + xi_169;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049] = xi_170 + xi_171;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 1023] = xi_173 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_175 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_177 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_179 + xi_184;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1026] = xi_181 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1026] = xi_183 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025] = xi_184 + xi_185;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050] = xi_187 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050] = xi_189 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049] = xi_190 + xi_191;
         _data_vertexFineDst[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_162*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 512; ctr_2 < 513; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_103 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025];
         const double xi_105 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049];
         const double xi_107 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025];
         const double xi_109 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050];
         const double xi_112 = -0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050];
         const double xi_113 = _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049];
         const double xi_102 = xi_97*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_104 = xi_97*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_106 = xi_99*0.375*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = xi_99*-0.125*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025] = xi_102 + xi_103;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049] = xi_104 + xi_105;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1025] = xi_106 + xi_107;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050] = xi_108 + xi_109;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + ((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2050] = xi_111 + xi_112;
         _data_edgeFineDst[2*ctr_1 + 2050*ctr_2 + 2*((1049600) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 2049] = xi_112 + xi_113;
         _data_vertexFineDst[2*ctr_1 + 2052*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_100*_data_vertexCoarseSrc[ctr_1 + 514*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_10(double * _data_edgeFineDst, double * _data_vertexCoarseSrc, double * _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_55 = 1 / (num_neighbor_faces_edge0);
   const double xi_57 = 1 / (num_neighbor_faces_edge2);
   const double xi_58 = 1 / (num_neighbor_faces_vertex0);
   const double xi_125 = 1 / (num_neighbor_faces_edge0);
   const double xi_76 = 1 / (num_neighbor_faces_edge0);
   const double xi_78 = 1 / (num_neighbor_faces_edge1);
   const double xi_79 = 1 / (num_neighbor_faces_vertex1);
   const double xi_199 = 1 / (num_neighbor_faces_edge2);
   const double xi_162 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_99 = 1 / (num_neighbor_faces_edge2);
   const double xi_100 = 1 / (num_neighbor_faces_vertex2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_63 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_65 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_67 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_70 = -0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_69 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_71 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_60 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_57*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = xi_57*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_60 + xi_61;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_62 + xi_63;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049] = xi_64 + xi_65;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_66 + xi_67;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049] = xi_69 + xi_70;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_70 + xi_71;
         _data_vertexFineDst[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_58*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 1024; ctr_1 += 1)
      {
         const double xi_128 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_130 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_134 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_153 = -0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047];
         const double xi_138 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047];
         const double xi_140 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
         const double xi_142 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
         const double xi_144 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_146 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_148 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_151 = 0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_150 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_152 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_154 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_133 = xi_125*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_131 = xi_125*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_128 + xi_133;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_130 + xi_131;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_131 + xi_132;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_133 + xi_134;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047] = xi_136 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047] = xi_138 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048] = xi_140 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048] = xi_142 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049] = xi_144 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049] = xi_146 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_148 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_150 + xi_151;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_151 + xi_152;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_153 + xi_154;
         _data_vertexFineDst[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 1024; ctr_1 < 1025; ctr_1 += 1)
      {
         const double xi_82 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_84 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_86 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047];
         const double xi_88 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_91 = -0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_90 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047];
         const double xi_92 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_81 = xi_76*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = xi_76*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_85 = xi_78*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = xi_78*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_81 + xi_82;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047] = xi_85 + xi_86;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_87 + xi_88;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047] = xi_90 + xi_91;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_91 + xi_92;
         _data_vertexFineDst[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_79*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 1024; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_202 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_204 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_206 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049];
         const double xi_208 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098];
         const double xi_227 = -0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_210 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_217 = 0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_212 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_216 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_218 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049];
         const double xi_220 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
         const double xi_222 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
         const double xi_224 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098];
         const double xi_226 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097];
         const double xi_228 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097];
         const double xi_207 = xi_199*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_205 = xi_199*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049] = xi_202 + xi_207;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_204 + xi_205;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049] = xi_205 + xi_206;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098] = xi_207 + xi_208;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049] = xi_210 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_212 + xi_217;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_214 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_216 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049] = xi_217 + xi_218;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048] = xi_220 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048] = xi_222 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098] = xi_224 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097] = xi_226 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097] = xi_227 + xi_228;
         _data_vertexFineDst[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_199*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 1024; ctr_1 += 1)
      {
         const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047];
         const double xi_6 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047];
         const double xi_8 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
         const double xi_10 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048];
         const double xi_12 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_14 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049];
         const double xi_16 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_18 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_22 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_26 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_32 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2050];
         const double xi_34 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2050];
         const double xi_36 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049];
         const double xi_38 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049];
         const double xi_40 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
         const double xi_42 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048];
         const double xi_44 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098];
         const double xi_46 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098];
         const double xi_48 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097];
         const double xi_50 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097];
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047] = xi_4 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047] = xi_49 + xi_6;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048] = xi_49 + xi_8;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2048] = xi_10 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049] = xi_12 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2049] = xi_14 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2050] = xi_32 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2050] = xi_34 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049] = xi_36 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049] = xi_37 + xi_38;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048] = xi_40 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2048] = xi_42 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098] = xi_44 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098] = xi_46 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097] = xi_48 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097] = xi_49 + xi_50;
         _data_vertexFineDst[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 1024; ctr_1 < -ctr_2 + 1025; ctr_1 += 1)
      {
         const double xi_165 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047];
         const double xi_167 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_169 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049];
         const double xi_171 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097];
         const double xi_190 = -0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047];
         const double xi_175 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_177 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_184 = 0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_179 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_181 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2050];
         const double xi_183 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2050];
         const double xi_185 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049];
         const double xi_187 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098];
         const double xi_189 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098];
         const double xi_191 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097];
         const double xi_170 = xi_162*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_168 = xi_162*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047] = xi_165 + xi_170;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_167 + xi_168;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049] = xi_168 + xi_169;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097] = xi_170 + xi_171;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2047] = xi_173 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_175 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_177 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_179 + xi_184;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2050] = xi_181 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2050] = xi_183 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049] = xi_184 + xi_185;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098] = xi_187 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098] = xi_189 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097] = xi_190 + xi_191;
         _data_vertexFineDst[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_162*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1024; ctr_2 < 1025; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_103 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049];
         const double xi_105 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097];
         const double xi_107 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049];
         const double xi_109 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098];
         const double xi_112 = -0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098];
         const double xi_113 = _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097];
         const double xi_102 = xi_97*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_104 = xi_97*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_106 = xi_99*0.375*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = xi_99*-0.125*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049] = xi_102 + xi_103;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097] = xi_104 + xi_105;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 2049] = xi_106 + xi_107;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098] = xi_108 + xi_109;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + ((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4098] = xi_111 + xi_112;
         _data_edgeFineDst[2*ctr_1 + 4098*ctr_2 + 2*((4196352) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 4097] = xi_112 + xi_113;
         _data_vertexFineDst[2*ctr_1 + 4100*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_100*_data_vertexCoarseSrc[ctr_1 + 1026*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_11(double * _data_edgeFineDst, double * _data_vertexCoarseSrc, double * _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_55 = 1 / (num_neighbor_faces_edge0);
   const double xi_57 = 1 / (num_neighbor_faces_edge2);
   const double xi_58 = 1 / (num_neighbor_faces_vertex0);
   const double xi_125 = 1 / (num_neighbor_faces_edge0);
   const double xi_76 = 1 / (num_neighbor_faces_edge0);
   const double xi_78 = 1 / (num_neighbor_faces_edge1);
   const double xi_79 = 1 / (num_neighbor_faces_vertex1);
   const double xi_199 = 1 / (num_neighbor_faces_edge2);
   const double xi_162 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_99 = 1 / (num_neighbor_faces_edge2);
   const double xi_100 = 1 / (num_neighbor_faces_vertex2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_63 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_65 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_67 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_70 = -0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_69 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_71 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_60 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_57*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = xi_57*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_60 + xi_61;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_62 + xi_63;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097] = xi_64 + xi_65;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_66 + xi_67;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097] = xi_69 + xi_70;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_70 + xi_71;
         _data_vertexFineDst[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_58*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 2048; ctr_1 += 1)
      {
         const double xi_128 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_130 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_134 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_153 = -0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095];
         const double xi_138 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095];
         const double xi_140 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
         const double xi_142 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
         const double xi_144 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_146 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_148 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_151 = 0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_150 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_152 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_154 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_133 = xi_125*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_131 = xi_125*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_128 + xi_133;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_130 + xi_131;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_131 + xi_132;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_133 + xi_134;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095] = xi_136 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095] = xi_138 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096] = xi_140 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096] = xi_142 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097] = xi_144 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097] = xi_146 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_148 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_150 + xi_151;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_151 + xi_152;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_153 + xi_154;
         _data_vertexFineDst[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 2048; ctr_1 < 2049; ctr_1 += 1)
      {
         const double xi_82 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_84 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_86 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095];
         const double xi_88 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_91 = -0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_90 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095];
         const double xi_92 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_81 = xi_76*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = xi_76*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_85 = xi_78*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = xi_78*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_81 + xi_82;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095] = xi_85 + xi_86;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_87 + xi_88;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095] = xi_90 + xi_91;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_91 + xi_92;
         _data_vertexFineDst[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_79*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 2048; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_202 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_204 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_206 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097];
         const double xi_208 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194];
         const double xi_227 = -0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_210 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_217 = 0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_212 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_216 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_218 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097];
         const double xi_220 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
         const double xi_222 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
         const double xi_224 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194];
         const double xi_226 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193];
         const double xi_228 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193];
         const double xi_207 = xi_199*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_205 = xi_199*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097] = xi_202 + xi_207;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_204 + xi_205;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097] = xi_205 + xi_206;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194] = xi_207 + xi_208;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097] = xi_210 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_212 + xi_217;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_214 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_216 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097] = xi_217 + xi_218;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096] = xi_220 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096] = xi_222 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194] = xi_224 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193] = xi_226 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193] = xi_227 + xi_228;
         _data_vertexFineDst[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_199*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 2048; ctr_1 += 1)
      {
         const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095];
         const double xi_6 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095];
         const double xi_8 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
         const double xi_10 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096];
         const double xi_12 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_14 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097];
         const double xi_16 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_18 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_22 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_26 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_32 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4098];
         const double xi_34 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4098];
         const double xi_36 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097];
         const double xi_38 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097];
         const double xi_40 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
         const double xi_42 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096];
         const double xi_44 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194];
         const double xi_46 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194];
         const double xi_48 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193];
         const double xi_50 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193];
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095] = xi_4 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095] = xi_49 + xi_6;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096] = xi_49 + xi_8;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4096] = xi_10 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097] = xi_12 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4097] = xi_14 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4098] = xi_32 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4098] = xi_34 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097] = xi_36 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097] = xi_37 + xi_38;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096] = xi_40 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4096] = xi_42 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194] = xi_44 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194] = xi_46 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193] = xi_48 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193] = xi_49 + xi_50;
         _data_vertexFineDst[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 2048; ctr_1 < -ctr_2 + 2049; ctr_1 += 1)
      {
         const double xi_165 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095];
         const double xi_167 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_169 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097];
         const double xi_171 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193];
         const double xi_190 = -0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095];
         const double xi_175 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_177 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_184 = 0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_179 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_181 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4098];
         const double xi_183 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4098];
         const double xi_185 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097];
         const double xi_187 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194];
         const double xi_189 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194];
         const double xi_191 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193];
         const double xi_170 = xi_162*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_168 = xi_162*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095] = xi_165 + xi_170;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_167 + xi_168;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097] = xi_168 + xi_169;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193] = xi_170 + xi_171;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 4095] = xi_173 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_175 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_177 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_179 + xi_184;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4098] = xi_181 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4098] = xi_183 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097] = xi_184 + xi_185;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194] = xi_187 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194] = xi_189 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193] = xi_190 + xi_191;
         _data_vertexFineDst[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_162*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 2048; ctr_2 < 2049; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_103 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097];
         const double xi_105 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193];
         const double xi_107 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097];
         const double xi_109 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194];
         const double xi_112 = -0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194];
         const double xi_113 = _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193];
         const double xi_102 = xi_97*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_104 = xi_97*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_106 = xi_99*0.375*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = xi_99*-0.125*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097] = xi_102 + xi_103;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193] = xi_104 + xi_105;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 4097] = xi_106 + xi_107;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194] = xi_108 + xi_109;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + ((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8194] = xi_111 + xi_112;
         _data_edgeFineDst[2*ctr_1 + 8194*ctr_2 + 2*((16781312) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 8193] = xi_112 + xi_113;
         _data_vertexFineDst[2*ctr_1 + 8196*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_100*_data_vertexCoarseSrc[ctr_1 + 2050*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_12(double * _data_edgeFineDst, double * _data_vertexCoarseSrc, double * _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_55 = 1 / (num_neighbor_faces_edge0);
   const double xi_57 = 1 / (num_neighbor_faces_edge2);
   const double xi_58 = 1 / (num_neighbor_faces_vertex0);
   const double xi_125 = 1 / (num_neighbor_faces_edge0);
   const double xi_76 = 1 / (num_neighbor_faces_edge0);
   const double xi_78 = 1 / (num_neighbor_faces_edge1);
   const double xi_79 = 1 / (num_neighbor_faces_vertex1);
   const double xi_199 = 1 / (num_neighbor_faces_edge2);
   const double xi_162 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_99 = 1 / (num_neighbor_faces_edge2);
   const double xi_100 = 1 / (num_neighbor_faces_vertex2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_63 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_65 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_67 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_70 = -0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_69 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_71 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_60 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_57*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = xi_57*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_60 + xi_61;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_62 + xi_63;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193] = xi_64 + xi_65;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_66 + xi_67;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193] = xi_69 + xi_70;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_70 + xi_71;
         _data_vertexFineDst[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_58*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 4096; ctr_1 += 1)
      {
         const double xi_128 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_130 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_134 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_153 = -0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191];
         const double xi_138 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191];
         const double xi_140 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
         const double xi_142 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
         const double xi_144 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_146 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_148 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_151 = 0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_150 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_152 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_154 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_133 = xi_125*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_131 = xi_125*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_128 + xi_133;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_130 + xi_131;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_131 + xi_132;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_133 + xi_134;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191] = xi_136 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191] = xi_138 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192] = xi_140 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192] = xi_142 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193] = xi_144 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193] = xi_146 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_148 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_150 + xi_151;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_151 + xi_152;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_153 + xi_154;
         _data_vertexFineDst[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 4096; ctr_1 < 4097; ctr_1 += 1)
      {
         const double xi_82 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_84 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_86 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191];
         const double xi_88 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_91 = -0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_90 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191];
         const double xi_92 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_81 = xi_76*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = xi_76*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_85 = xi_78*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = xi_78*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_81 + xi_82;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191] = xi_85 + xi_86;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_87 + xi_88;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191] = xi_90 + xi_91;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_91 + xi_92;
         _data_vertexFineDst[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_79*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 4096; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_202 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_204 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_206 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193];
         const double xi_208 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386];
         const double xi_227 = -0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_210 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_217 = 0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_212 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_216 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_218 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193];
         const double xi_220 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
         const double xi_222 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
         const double xi_224 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386];
         const double xi_226 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385];
         const double xi_228 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385];
         const double xi_207 = xi_199*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_205 = xi_199*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193] = xi_202 + xi_207;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_204 + xi_205;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193] = xi_205 + xi_206;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386] = xi_207 + xi_208;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193] = xi_210 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_212 + xi_217;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_214 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_216 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193] = xi_217 + xi_218;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192] = xi_220 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192] = xi_222 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386] = xi_224 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385] = xi_226 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385] = xi_227 + xi_228;
         _data_vertexFineDst[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_199*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 4096; ctr_1 += 1)
      {
         const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191];
         const double xi_6 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191];
         const double xi_8 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
         const double xi_10 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192];
         const double xi_12 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_14 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193];
         const double xi_16 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_18 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_22 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_26 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_32 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8194];
         const double xi_34 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8194];
         const double xi_36 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193];
         const double xi_38 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193];
         const double xi_40 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
         const double xi_42 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192];
         const double xi_44 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386];
         const double xi_46 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386];
         const double xi_48 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385];
         const double xi_50 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385];
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191] = xi_4 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191] = xi_49 + xi_6;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192] = xi_49 + xi_8;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8192] = xi_10 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193] = xi_12 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8193] = xi_14 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8194] = xi_32 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8194] = xi_34 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193] = xi_36 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193] = xi_37 + xi_38;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192] = xi_40 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8192] = xi_42 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386] = xi_44 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386] = xi_46 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385] = xi_48 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385] = xi_49 + xi_50;
         _data_vertexFineDst[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 4096; ctr_1 < -ctr_2 + 4097; ctr_1 += 1)
      {
         const double xi_165 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191];
         const double xi_167 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_169 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193];
         const double xi_171 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385];
         const double xi_190 = -0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191];
         const double xi_175 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_177 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_184 = 0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_179 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_181 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8194];
         const double xi_183 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8194];
         const double xi_185 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193];
         const double xi_187 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386];
         const double xi_189 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386];
         const double xi_191 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385];
         const double xi_170 = xi_162*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_168 = xi_162*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191] = xi_165 + xi_170;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_167 + xi_168;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193] = xi_168 + xi_169;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385] = xi_170 + xi_171;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 8191] = xi_173 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_175 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_177 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_179 + xi_184;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8194] = xi_181 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8194] = xi_183 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193] = xi_184 + xi_185;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386] = xi_187 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386] = xi_189 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385] = xi_190 + xi_191;
         _data_vertexFineDst[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_162*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 4096; ctr_2 < 4097; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_103 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193];
         const double xi_105 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385];
         const double xi_107 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193];
         const double xi_109 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386];
         const double xi_112 = -0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386];
         const double xi_113 = _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385];
         const double xi_102 = xi_97*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_104 = xi_97*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_106 = xi_99*0.375*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = xi_99*-0.125*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193] = xi_102 + xi_103;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385] = xi_104 + xi_105;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 8193] = xi_106 + xi_107;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386] = xi_108 + xi_109;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + ((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16386] = xi_111 + xi_112;
         _data_edgeFineDst[2*ctr_1 + 16386*ctr_2 + 2*((67117056) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 16385] = xi_112 + xi_113;
         _data_vertexFineDst[2*ctr_1 + 16388*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_100*_data_vertexCoarseSrc[ctr_1 + 4098*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_13(double * _data_edgeFineDst, double * _data_vertexCoarseSrc, double * _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_55 = 1 / (num_neighbor_faces_edge0);
   const double xi_57 = 1 / (num_neighbor_faces_edge2);
   const double xi_58 = 1 / (num_neighbor_faces_vertex0);
   const double xi_125 = 1 / (num_neighbor_faces_edge0);
   const double xi_76 = 1 / (num_neighbor_faces_edge0);
   const double xi_78 = 1 / (num_neighbor_faces_edge1);
   const double xi_79 = 1 / (num_neighbor_faces_vertex1);
   const double xi_199 = 1 / (num_neighbor_faces_edge2);
   const double xi_162 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_99 = 1 / (num_neighbor_faces_edge2);
   const double xi_100 = 1 / (num_neighbor_faces_vertex2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_63 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_65 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_67 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_70 = -0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_69 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_71 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_60 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_57*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = xi_57*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_60 + xi_61;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_62 + xi_63;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385] = xi_64 + xi_65;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_66 + xi_67;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385] = xi_69 + xi_70;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_70 + xi_71;
         _data_vertexFineDst[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_58*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 8192; ctr_1 += 1)
      {
         const double xi_128 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_130 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_134 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_153 = -0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383];
         const double xi_138 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383];
         const double xi_140 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
         const double xi_142 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
         const double xi_144 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_146 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_148 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_151 = 0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_150 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_152 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_154 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_133 = xi_125*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_131 = xi_125*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_128 + xi_133;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_130 + xi_131;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_131 + xi_132;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_133 + xi_134;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383] = xi_136 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383] = xi_138 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384] = xi_140 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384] = xi_142 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385] = xi_144 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385] = xi_146 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_148 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_150 + xi_151;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_151 + xi_152;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_153 + xi_154;
         _data_vertexFineDst[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 8192; ctr_1 < 8193; ctr_1 += 1)
      {
         const double xi_82 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_84 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_86 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383];
         const double xi_88 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_91 = -0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_90 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383];
         const double xi_92 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_81 = xi_76*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = xi_76*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_85 = xi_78*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = xi_78*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_81 + xi_82;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383] = xi_85 + xi_86;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_87 + xi_88;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383] = xi_90 + xi_91;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_91 + xi_92;
         _data_vertexFineDst[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_79*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 8192; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_202 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_204 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_206 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385];
         const double xi_208 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770];
         const double xi_227 = -0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_210 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_217 = 0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_212 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_216 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_218 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385];
         const double xi_220 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
         const double xi_222 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
         const double xi_224 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770];
         const double xi_226 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769];
         const double xi_228 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769];
         const double xi_207 = xi_199*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_205 = xi_199*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385] = xi_202 + xi_207;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_204 + xi_205;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385] = xi_205 + xi_206;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770] = xi_207 + xi_208;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385] = xi_210 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_212 + xi_217;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_214 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_216 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385] = xi_217 + xi_218;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384] = xi_220 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384] = xi_222 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770] = xi_224 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769] = xi_226 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769] = xi_227 + xi_228;
         _data_vertexFineDst[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_199*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 8192; ctr_1 += 1)
      {
         const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383];
         const double xi_6 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383];
         const double xi_8 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
         const double xi_10 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384];
         const double xi_12 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_14 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385];
         const double xi_16 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_18 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_22 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_26 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_32 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16386];
         const double xi_34 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16386];
         const double xi_36 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385];
         const double xi_38 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385];
         const double xi_40 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
         const double xi_42 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384];
         const double xi_44 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770];
         const double xi_46 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770];
         const double xi_48 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769];
         const double xi_50 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769];
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383] = xi_4 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383] = xi_49 + xi_6;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384] = xi_49 + xi_8;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16384] = xi_10 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385] = xi_12 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16385] = xi_14 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16386] = xi_32 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16386] = xi_34 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385] = xi_36 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385] = xi_37 + xi_38;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384] = xi_40 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16384] = xi_42 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770] = xi_44 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770] = xi_46 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769] = xi_48 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769] = xi_49 + xi_50;
         _data_vertexFineDst[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 8192; ctr_1 < -ctr_2 + 8193; ctr_1 += 1)
      {
         const double xi_165 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383];
         const double xi_167 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_169 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385];
         const double xi_171 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769];
         const double xi_190 = -0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383];
         const double xi_175 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_177 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_184 = 0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_179 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_181 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16386];
         const double xi_183 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16386];
         const double xi_185 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385];
         const double xi_187 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770];
         const double xi_189 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770];
         const double xi_191 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769];
         const double xi_170 = xi_162*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_168 = xi_162*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383] = xi_165 + xi_170;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_167 + xi_168;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385] = xi_168 + xi_169;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769] = xi_170 + xi_171;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 16383] = xi_173 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_175 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_177 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_179 + xi_184;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16386] = xi_181 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16386] = xi_183 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385] = xi_184 + xi_185;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770] = xi_187 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770] = xi_189 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769] = xi_190 + xi_191;
         _data_vertexFineDst[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_162*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 8192; ctr_2 < 8193; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_103 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385];
         const double xi_105 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769];
         const double xi_107 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385];
         const double xi_109 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770];
         const double xi_112 = -0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770];
         const double xi_113 = _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769];
         const double xi_102 = xi_97*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_104 = xi_97*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_106 = xi_99*0.375*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = xi_99*-0.125*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385] = xi_102 + xi_103;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769] = xi_104 + xi_105;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 16385] = xi_106 + xi_107;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770] = xi_108 + xi_109;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + ((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32770] = xi_111 + xi_112;
         _data_edgeFineDst[2*ctr_1 + 32770*ctr_2 + 2*((268451840) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 32769] = xi_112 + xi_113;
         _data_vertexFineDst[2*ctr_1 + 32772*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_100*_data_vertexCoarseSrc[ctr_1 + 8194*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_14(double * _data_edgeFineDst, double * _data_vertexCoarseSrc, double * _data_vertexFineDst, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_55 = 1 / (num_neighbor_faces_edge0);
   const double xi_57 = 1 / (num_neighbor_faces_edge2);
   const double xi_58 = 1 / (num_neighbor_faces_vertex0);
   const double xi_125 = 1 / (num_neighbor_faces_edge0);
   const double xi_76 = 1 / (num_neighbor_faces_edge0);
   const double xi_78 = 1 / (num_neighbor_faces_edge1);
   const double xi_79 = 1 / (num_neighbor_faces_vertex1);
   const double xi_199 = 1 / (num_neighbor_faces_edge2);
   const double xi_162 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_99 = 1 / (num_neighbor_faces_edge2);
   const double xi_100 = 1 / (num_neighbor_faces_vertex2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_63 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_65 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_67 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_70 = -0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_69 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_71 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_60 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_57*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = xi_57*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_60 + xi_61;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_62 + xi_63;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769] = xi_64 + xi_65;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_66 + xi_67;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769] = xi_69 + xi_70;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_70 + xi_71;
         _data_vertexFineDst[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_58*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < 16384; ctr_1 += 1)
      {
         const double xi_128 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_130 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_134 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_153 = -0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767];
         const double xi_138 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767];
         const double xi_140 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
         const double xi_142 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
         const double xi_144 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_146 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_148 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_151 = 0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_150 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_152 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_154 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_133 = xi_125*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_131 = xi_125*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_128 + xi_133;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_130 + xi_131;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_131 + xi_132;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_133 + xi_134;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767] = xi_136 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767] = xi_138 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768] = xi_140 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768] = xi_142 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769] = xi_144 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769] = xi_146 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_148 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_150 + xi_151;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_151 + xi_152;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_153 + xi_154;
         _data_vertexFineDst[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = 16384; ctr_1 < 16385; ctr_1 += 1)
      {
         const double xi_82 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_84 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_86 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767];
         const double xi_88 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_91 = -0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_90 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767];
         const double xi_92 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_81 = xi_76*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = xi_76*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_85 = xi_78*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = xi_78*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_81 + xi_82;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767] = xi_85 + xi_86;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_87 + xi_88;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767] = xi_90 + xi_91;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_91 + xi_92;
         _data_vertexFineDst[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_79*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < 16384; ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_202 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_204 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_206 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769];
         const double xi_208 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538];
         const double xi_227 = -0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_210 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_217 = 0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_212 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_216 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_218 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769];
         const double xi_220 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
         const double xi_222 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
         const double xi_224 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538];
         const double xi_226 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537];
         const double xi_228 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537];
         const double xi_207 = xi_199*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_205 = xi_199*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769] = xi_202 + xi_207;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_204 + xi_205;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769] = xi_205 + xi_206;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538] = xi_207 + xi_208;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769] = xi_210 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_212 + xi_217;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_214 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_216 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769] = xi_217 + xi_218;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768] = xi_220 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768] = xi_222 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538] = xi_224 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537] = xi_226 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537] = xi_227 + xi_228;
         _data_vertexFineDst[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_199*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + 16384; ctr_1 += 1)
      {
         const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767];
         const double xi_6 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767];
         const double xi_8 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
         const double xi_10 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768];
         const double xi_12 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_14 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769];
         const double xi_16 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_18 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_22 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_24 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_26 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_30 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_32 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32770];
         const double xi_34 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32770];
         const double xi_36 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769];
         const double xi_38 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769];
         const double xi_40 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
         const double xi_42 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768];
         const double xi_44 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538];
         const double xi_46 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538];
         const double xi_48 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537];
         const double xi_50 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537];
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767] = xi_4 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767] = xi_49 + xi_6;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768] = xi_49 + xi_8;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32768] = xi_10 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769] = xi_12 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32769] = xi_14 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_18 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_22 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_24 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_28 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32770] = xi_32 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32770] = xi_34 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769] = xi_36 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769] = xi_37 + xi_38;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768] = xi_40 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32768] = xi_42 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538] = xi_44 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538] = xi_46 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537] = xi_48 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537] = xi_49 + xi_50;
         _data_vertexFineDst[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + 16384; ctr_1 < -ctr_2 + 16385; ctr_1 += 1)
      {
         const double xi_165 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767];
         const double xi_167 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_169 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769];
         const double xi_171 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537];
         const double xi_190 = -0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767];
         const double xi_175 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_177 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_184 = 0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_179 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_181 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32770];
         const double xi_183 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32770];
         const double xi_185 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769];
         const double xi_187 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538];
         const double xi_189 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538];
         const double xi_191 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537];
         const double xi_170 = xi_162*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_168 = xi_162*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767] = xi_165 + xi_170;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_167 + xi_168;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769] = xi_168 + xi_169;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537] = xi_170 + xi_171;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 32767] = xi_173 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_175 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_177 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_179 + xi_184;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32770] = xi_181 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32770] = xi_183 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769] = xi_184 + xi_185;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538] = xi_187 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538] = xi_189 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537] = xi_190 + xi_191;
         _data_vertexFineDst[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_162*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 16384; ctr_2 < 16385; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_103 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769];
         const double xi_105 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537];
         const double xi_107 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769];
         const double xi_109 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538];
         const double xi_112 = -0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538];
         const double xi_113 = _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537];
         const double xi_102 = xi_97*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_104 = xi_97*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_106 = xi_99*0.375*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = xi_99*-0.125*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769] = xi_102 + xi_103;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537] = xi_104 + xi_105;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 32769] = xi_106 + xi_107;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538] = xi_108 + xi_109;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + ((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65538] = xi_111 + xi_112;
         _data_edgeFineDst[2*ctr_1 + 65538*ctr_2 + 2*((1073774592) / (2)) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) - 65537] = xi_112 + xi_113;
         _data_vertexFineDst[2*ctr_1 + 65540*ctr_2 - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_100*_data_vertexCoarseSrc[ctr_1 + 16386*ctr_2 - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}

static void prolongate_2D_macroface_P2_push_from_vertexdofs_level_any(double * _data_edgeFineDst, double * _data_vertexCoarseSrc, double * _data_vertexFineDst, int64_t coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
   const double xi_55 = 1 / (num_neighbor_faces_edge0);
   const double xi_57 = 1 / (num_neighbor_faces_edge2);
   const double xi_58 = 1 / (num_neighbor_faces_vertex0);
   const double xi_125 = 1 / (num_neighbor_faces_edge0);
   const double xi_76 = 1 / (num_neighbor_faces_edge0);
   const double xi_78 = 1 / (num_neighbor_faces_edge1);
   const double xi_79 = 1 / (num_neighbor_faces_vertex1);
   const double xi_199 = 1 / (num_neighbor_faces_edge2);
   const double xi_162 = 1 / (num_neighbor_faces_edge1);
   const double xi_97 = 1 / (num_neighbor_faces_edge1);
   const double xi_99 = 1 / (num_neighbor_faces_edge2);
   const double xi_100 = 1 / (num_neighbor_faces_vertex2);
   for (int ctr_2 = 0; ctr_2 < 1; ctr_2 += 1)
   {
      // bottom left vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_61 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_63 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_65 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_67 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_70 = -0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_69 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_71 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_60 = xi_55*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_62 = xi_55*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_64 = xi_57*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_66 = xi_57*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_60 + xi_61;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_62 + xi_63;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_64 + xi_65;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_66 + xi_67;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_69 + xi_70;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1] = xi_70 + xi_71;
         _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_58*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom edge
      for (int ctr_1 = 1; ctr_1 < (1 << (coarse_level)); ctr_1 += 1)
      {
         const double xi_128 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_130 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_132 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_134 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_153 = -0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_136 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
         const double xi_138 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 2];
         const double xi_140 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_142 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_144 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_146 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_148 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 2];
         const double xi_151 = 0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_150 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_152 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_154 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_133 = xi_125*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_131 = xi_125*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_128 + xi_133;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_130 + xi_131;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_131 + xi_132;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_133 + xi_134;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2] = xi_136 + xi_153;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 2] = xi_138 + xi_153;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1] = xi_140 + xi_153;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1] = xi_142 + xi_153;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_144 + xi_153;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_146 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 2] = xi_148 + xi_153;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1] = xi_150 + xi_151;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_151 + xi_152;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1] = xi_153 + xi_154;
         _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // bottom right vertex
      for (int ctr_1 = (1 << (coarse_level)); ctr_1 < (1 << (coarse_level)) + 1; ctr_1 += 1)
      {
         const double xi_82 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_84 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_86 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 2];
         const double xi_88 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_91 = -0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_90 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
         const double xi_92 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 2];
         const double xi_81 = xi_76*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_83 = xi_76*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_85 = xi_78*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_87 = xi_78*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_81 + xi_82;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_83 + xi_84;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 2] = xi_85 + xi_86;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1] = xi_87 + xi_88;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2] = xi_90 + xi_91;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 2] = xi_91 + xi_92;
         _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_79*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = 1; ctr_2 < (1 << (coarse_level)); ctr_2 += 1)
   {
      // left edge
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_202 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_204 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_206 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_208 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_227 = -0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_210 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_217 = 0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_212 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_214 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_216 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_218 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_220 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
         const double xi_222 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_224 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_226 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_228 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_207 = xi_199*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_205 = xi_199*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_202 + xi_207;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_204 + xi_205;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_205 + xi_206;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_207 + xi_208;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_210 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_212 + xi_217;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1] = xi_214 + xi_227;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_216 + xi_227;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_217 + xi_218;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_220 + xi_227;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1] = xi_222 + xi_227;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_224 + xi_227;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1] = xi_226 + xi_227;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1] = xi_227 + xi_228;
         _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_199*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // inner triangle
      for (int ctr_1 = 1; ctr_1 < -ctr_2 + (1 << (coarse_level)); ctr_1 += 1)
      {
         const double xi_49 = -0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_4 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
         const double xi_6 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 2];
         const double xi_8 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_10 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_12 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_14 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))];
         const double xi_16 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_18 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 2];
         const double xi_37 = 0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_20 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_22 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_24 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_26 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))];
         const double xi_28 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_30 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1];
         const double xi_32 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_34 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1];
         const double xi_36 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_38 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_40 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1];
         const double xi_42 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_44 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_46 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_48 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_50 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2] = xi_4 + xi_49;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 2] = xi_49 + xi_6;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1] = xi_49 + xi_8;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1] = xi_10 + xi_49;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_12 + xi_49;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2))] = xi_14 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_16 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 2] = xi_18 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_20 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1] = xi_22 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_24 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = xi_26 + xi_37;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1] = xi_28 + xi_49;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + 1] = xi_30 + xi_49;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1] = xi_32 + xi_49;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1] = xi_34 + xi_49;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_36 + xi_37;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_37 + xi_38;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 1] = xi_40 + xi_49;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1] = xi_42 + xi_49;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_44 + xi_49;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_46 + xi_49;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1] = xi_48 + xi_49;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1] = xi_49 + xi_50;
         _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = _data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
      // diagonal edge
      for (int ctr_1 = -ctr_2 + (1 << (coarse_level)); ctr_1 < -ctr_2 + (1 << (coarse_level)) + 1; ctr_1 += 1)
      {
         const double xi_165 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 2];
         const double xi_167 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_169 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_171 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_190 = -0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_173 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2];
         const double xi_175 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2];
         const double xi_177 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 2];
         const double xi_184 = 0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_179 = _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1];
         const double xi_181 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1];
         const double xi_183 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1];
         const double xi_185 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_187 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_189 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_191 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_170 = xi_162*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_168 = xi_162*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 2] = xi_165 + xi_170;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1] = xi_167 + xi_168;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_168 + xi_169;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1] = xi_170 + xi_171;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 + 1)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 + 1)*(2*ctr_2 + 2)) / (2)) - 2] = xi_173 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 2] = xi_175 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 2] = xi_177 + xi_190;
         _data_edgeFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 + 1)) / (2)) - 1] = xi_179 + xi_184;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) - 1] = xi_181 + xi_190;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) - 1] = xi_183 + xi_190;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_184 + xi_185;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_187 + xi_190;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_189 + xi_190;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1] = xi_190 + xi_191;
         _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_162*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   for (int ctr_2 = (1 << (coarse_level)); ctr_2 < (1 << (coarse_level)) + 1; ctr_2 += 1)
   {
      // top vertex
      for (int ctr_1 = 0; ctr_1 < 1; ctr_1 += 1)
      {
         const double xi_103 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_105 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_107 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_109 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_112 = -0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_111 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))];
         const double xi_113 = _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1];
         const double xi_102 = xi_97*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_104 = xi_97*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_106 = xi_99*0.375*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         const double xi_108 = xi_99*-0.125*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_102 + xi_103;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1] = xi_104 + xi_105;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 1)*((1 << (coarse_level + 1)) + 1) - ((2*ctr_2*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_106 + xi_107;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_108 + xi_109;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + ((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2))] = xi_111 + xi_112;
         _data_edgeFineDst[2*ctr_1 + (2*ctr_2 - 2)*((1 << (coarse_level + 1)) + 1) - (((2*ctr_2 - 2)*(2*ctr_2 - 1)) / (2)) + 2*((((1 << (coarse_level + 1)) + 1)*(1 << (coarse_level + 1))) / (2)) + 1] = xi_112 + xi_113;
         _data_vertexFineDst[2*ctr_1 + 2*ctr_2*((1 << (coarse_level + 1)) + 2) - ((2*ctr_2*(2*ctr_2 + 1)) / (2))] = 1.0*xi_100*_data_vertexCoarseSrc[ctr_1 + ctr_2*((1 << (coarse_level)) + 2) - ((ctr_2*(ctr_2 + 1)) / (2))];
      }
   }
   {
      
   }
}


void prolongate_2D_macroface_P2_push_from_vertexdofs(double * _data_edgeFineDst, double * _data_vertexCoarseSrc, double * _data_vertexFineDst, int64_t coarse_level, double num_neighbor_faces_edge0, double num_neighbor_faces_edge1, double num_neighbor_faces_edge2, double num_neighbor_faces_vertex0, double num_neighbor_faces_vertex1, double num_neighbor_faces_vertex2)
{
    switch( coarse_level )
    {
    case 2:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_2(_data_edgeFineDst, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 3:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_3(_data_edgeFineDst, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 4:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_4(_data_edgeFineDst, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 5:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_5(_data_edgeFineDst, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 6:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_6(_data_edgeFineDst, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 7:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_7(_data_edgeFineDst, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 8:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_8(_data_edgeFineDst, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 9:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_9(_data_edgeFineDst, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 10:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_10(_data_edgeFineDst, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 11:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_11(_data_edgeFineDst, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 12:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_12(_data_edgeFineDst, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 13:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_13(_data_edgeFineDst, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    case 14:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_14(_data_edgeFineDst, _data_vertexCoarseSrc, _data_vertexFineDst, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    default:
        prolongate_2D_macroface_P2_push_from_vertexdofs_level_any(_data_edgeFineDst, _data_vertexCoarseSrc, _data_vertexFineDst, coarse_level, num_neighbor_faces_edge0, num_neighbor_faces_edge1, num_neighbor_faces_edge2, num_neighbor_faces_vertex0, num_neighbor_faces_vertex1, num_neighbor_faces_vertex2);
        break;
    }
}
    

} // namespace generated
} // namespace macroface
} // namespace P2
} // namespace hhg